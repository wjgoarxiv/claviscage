from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "release_metadata.py"
ENGINE_ROOT = SCRIPT.parents[1]
SPEC = importlib.util.spec_from_file_location("release_metadata", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
release_metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release_metadata)


def test_source_tree_digest_is_content_addressed_and_ignores_build_output(tmp_path: Path):
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "lib.rs").write_text("pub fn value() -> u8 { 1 }\n")
    (tmp_path / "Cargo.toml").write_text("[package]\nname='x'\nversion='0.1.0'\n")
    first = release_metadata.source_tree_sha256(tmp_path)
    (tmp_path / "target").mkdir()
    (tmp_path / "target" / "artifact").write_bytes(b"ignored")
    assert release_metadata.source_tree_sha256(tmp_path) == first
    (tmp_path / "src" / "lib.rs").write_text("pub fn value() -> u8 { 2 }\n")
    assert release_metadata.source_tree_sha256(tmp_path) != first


def test_spdx_validation_rejects_a_wheel_hash_mismatch(tmp_path: Path):
    wheel = tmp_path / "clavis_engine-0.3.16-cp311-abi3-test.whl"
    wheel.write_bytes(b"wheel")
    sbom = tmp_path / "engine.spdx.json"
    sbom.write_text(
        json.dumps(
            {
                "spdxVersion": "SPDX-2.3",
                "packages": [{"name": "clavis-wasm", "versionInfo": "0.3.16"}],
                "files": [
                    {
                        "fileName": wheel.name,
                        "checksums": [{"algorithm": "SHA256", "checksumValue": "0" * 64}],
                    }
                ],
                "relationships": [{"spdxElementId": "SPDXRef-Package-clavis-wasm"}],
            }
        )
    )
    with pytest.raises(ValueError, match="wheel hash"):
        release_metadata.validate_spdx(sbom, [wheel], "0.3.16")


def test_spdx_validation_accepts_exact_wheel_hash(tmp_path: Path):
    wheel = tmp_path / "clavis_engine-0.3.16-cp311-abi3-test.whl"
    wheel.write_bytes(b"wheel")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    sbom = tmp_path / "engine.spdx.json"
    sbom.write_text(
        json.dumps(
            {
                "spdxVersion": "SPDX-2.3",
                "packages": [{"name": "clavis-wasm", "versionInfo": "0.3.16"}],
                "files": [
                    {
                        "fileName": wheel.name,
                        "checksums": [{"algorithm": "SHA256", "checksumValue": digest}],
                    }
                ],
                "relationships": [{"spdxElementId": "SPDXRef-Package-clavis-wasm"}],
            }
        )
    )
    release_metadata.validate_spdx(sbom, [wheel], "0.3.16")


def test_spdx_generation_includes_pyo3_graph_and_every_supplied_wheel(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    engine_wheel = tmp_path / "clavis_engine-0.3.16-cp311-abi3-test.whl"
    dependency_wheel = tmp_path / "numpy-2.5.1-cp313-cp313-test.whl"
    engine_wheel.write_bytes(b"engine")
    dependency_wheel.write_bytes(b"dependency")
    output = tmp_path / "engine.spdx.json"
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")

    document = release_metadata.generate_spdx(
        ENGINE_ROOT,
        [engine_wheel, dependency_wheel],
        output,
        "test-commit",
    )

    package_names = {item["name"] for item in document["packages"]}
    package_ids = {item["name"]: item["SPDXID"] for item in document["packages"]}
    declared_files = {item["fileName"] for item in document["files"]}
    declared_licenses = {
        item["fileName"]: item["licenseConcluded"]
        for item in document["files"]
    }
    dependency_edges = {
        (item["spdxElementId"], item["relatedSpdxElement"])
        for item in document["relationships"]
        if item["relationshipType"] == "DEPENDS_ON"
    }
    assert {"pyo3", "pyo3-ffi", "pyo3-macros"} <= package_names
    assert (package_ids["clavis-wasm"], package_ids["pyo3"]) in dependency_edges
    assert (package_ids["pyo3"], package_ids["pyo3-ffi"]) in dependency_edges
    assert declared_files == {engine_wheel.name, dependency_wheel.name}
    assert declared_licenses[engine_wheel.name] == "GPL-3.0-or-later"
    assert declared_licenses[dependency_wheel.name] == "NOASSERTION"
    assert document["creationInfo"]["created"] == "1970-01-01T00:00:00Z"
