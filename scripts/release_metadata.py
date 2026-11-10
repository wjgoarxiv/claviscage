#!/usr/bin/env python3
"""Create and validate content-addressed metadata for native release artifacts."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SOURCE_FILES = ("Cargo.toml", "Cargo.lock", "pyproject.toml", "build.rs", "LICENSE", "NOTICE")


class ReleaseMetadataError(ValueError):
    __slots__ = ("detail",)

    def __init__(self, detail: str) -> None:
        self.detail = detail

    def __str__(self) -> str:
        return self.detail


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_paths(root: Path) -> list[Path]:
    paths = [root / name for name in SOURCE_FILES if (root / name).is_file()]
    source_root = root / "src"
    if source_root.is_dir():
        paths.extend(path for path in source_root.rglob("*") if path.is_file())
    return sorted(paths, key=lambda path: path.relative_to(root).as_posix())


def source_tree_sha256(root: Path) -> str:
    """Hash only native build inputs, with paths and bytes framed unambiguously."""
    root = root.resolve()
    digest = hashlib.sha256()
    for path in _source_paths(root):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        content = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


def _spdx_id(name: str, version: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9.-]+", "-", f"{name}-{version}").strip("-")
    return f"SPDXRef-Package-{slug}"


def _creation_timestamp(engine_root: Path) -> str:
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch is None:
        epoch = subprocess.run(
            ["git", "show", "-s", "--format=%ct", "HEAD"],
            cwd=engine_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    return datetime.fromtimestamp(int(epoch), timezone.utc).isoformat().replace("+00:00", "Z")


def generate_spdx(
    engine_root: Path,
    wheels: list[Path],
    output: Path,
    source_commit: str,
) -> dict[str, Any]:
    metadata = json.loads(
        subprocess.run(
            ["cargo", "metadata", "--locked", "--format-version", "1", "--features", "pyo3"],
            cwd=engine_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )
    packages = []
    package_by_id: dict[str, dict[str, Any]] = {}
    for package in metadata["packages"]:
        spdx_id = _spdx_id(package["name"], package["version"])
        item = {
            "SPDXID": spdx_id,
            "name": package["name"],
            "versionInfo": package["version"],
            "downloadLocation": "NOASSERTION",
            "filesAnalyzed": False,
            "licenseConcluded": "NOASSERTION",
            "licenseDeclared": package.get("license") or "NOASSERTION",
            "supplier": "NOASSERTION",
        }
        packages.append(item)
        package_by_id[package["id"]] = item
    root_package = next(
        item for item in packages if item["name"] == "clavis-wasm" and item["versionInfo"] == "0.3.17"
    )
    relationships = [
        {
            "spdxElementId": "SPDXRef-DOCUMENT",
            "relationshipType": "DESCRIBES",
            "relatedSpdxElement": root_package["SPDXID"],
        }
    ]
    root_id = next(package["id"] for package in metadata["packages"] if package["name"] == "clavis-wasm")
    resolve = metadata["resolve"] or {"nodes": []}
    for node in resolve["nodes"]:
        source_package = package_by_id[node["id"]]
        for dependency in node["dependencies"]:
            relationships.append(
                {
                    "spdxElementId": source_package["SPDXID"],
                    "relationshipType": "DEPENDS_ON",
                    "relatedSpdxElement": package_by_id[dependency]["SPDXID"],
                }
            )
    tree_digest = source_tree_sha256(engine_root)
    files = [
        {
            "SPDXID": f"SPDXRef-File-{index}",
            "fileName": wheel.name,
            "checksums": [{"algorithm": "SHA256", "checksumValue": _sha256(wheel)}],
            "licenseConcluded": (
                "GPL-3.0-or-later" if wheel.name.startswith("clavis_engine-") else "NOASSERTION"
            ),
            "copyrightText": "NOASSERTION",
        }
        for index, wheel in enumerate(sorted(wheels), start=1)
    ]
    document = {
        "spdxVersion": "SPDX-2.3",
        "dataLicense": "CC0-1.0",
        "SPDXID": "SPDXRef-DOCUMENT",
        "name": "clavis-engine-0.3.17-release-artifacts",
        "documentNamespace": f"https://claviscage.app/spdx/0.3.17/{source_commit}/{tree_digest}",
        "creationInfo": {
            "created": _creation_timestamp(engine_root),
            "creators": ["Tool: ClavisCage release_metadata.py"],
        },
        "comment": f"Release build from {source_commit} plus source tree {tree_digest}.",
        "packages": packages,
        "files": files,
        "relationships": relationships,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return document


def validate_spdx(sbom: Path, wheels: list[Path], engine_version: str) -> None:
    document = json.loads(sbom.read_text(encoding="utf-8"))
    if document.get("spdxVersion") != "SPDX-2.3":
        raise ReleaseMetadataError("SBOM must use SPDX-2.3")
    packages = document.get("packages")
    if not isinstance(packages, list) or not any(
        item.get("name") == "clavis-wasm" and item.get("versionInfo") == engine_version
        for item in packages
        if isinstance(item, dict)
    ):
        raise ReleaseMetadataError("SBOM has no exact engine package")
    relationships = document.get("relationships")
    if not isinstance(relationships, list) or not relationships:
        raise ReleaseMetadataError("SBOM has no relationships")
    files = document.get("files")
    if not isinstance(files, list):
        raise ReleaseMetadataError("SBOM has no files")
    actual = {wheel.name: _sha256(wheel) for wheel in wheels}
    declared: dict[str, str] = {}
    for item in files:
        if not isinstance(item, dict) or not isinstance(item.get("fileName"), str):
            continue
        checksums = item.get("checksums")
        if not isinstance(checksums, list):
            continue
        sha = next(
            (
                checksum.get("checksumValue")
                for checksum in checksums
                if isinstance(checksum, dict) and checksum.get("algorithm") == "SHA256"
            ),
            None,
        )
        if isinstance(sha, str):
            declared[item["fileName"]] = sha
    if declared != actual:
        raise ReleaseMetadataError(
            f"wheel hash set does not match SBOM: declared={declared} actual={actual}"
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    identity = subparsers.add_parser("identity")
    identity.add_argument("--engine-root", type=Path, required=True)
    sbom = subparsers.add_parser("sbom")
    sbom.add_argument("--engine-root", type=Path, required=True)
    sbom.add_argument("--wheel", type=Path, action="append", required=True)
    sbom.add_argument("--output", type=Path, required=True)
    sbom.add_argument("--source-commit", required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("--sbom", type=Path, required=True)
    validate.add_argument("--wheel", type=Path, action="append", required=True)
    validate.add_argument("--engine-version", required=True)
    args = parser.parse_args()
    if args.command == "identity":
        print(json.dumps({"sourceTreeSha256": source_tree_sha256(args.engine_root)}, sort_keys=True))
        return 0
    if args.command == "sbom":
        generate_spdx(args.engine_root, args.wheel, args.output, args.source_commit)
        return 0
    validate_spdx(args.sbom, args.wheel, args.engine_version)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
