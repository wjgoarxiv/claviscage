"""Correctness tests for the native `clavis_engine` Python bindings (Slice 0).

The bridge is a thin pass-through to the exact same pure-Rust functions that back
`clavis-cli` and the WASM build, so faithfulness is proven by:

  * exact parity with `clavis-cli` on the DETERMINISTIC analyses (chill, f3f4);
  * self-consistency + structural sanity on the analyses whose cage assembly has
    pre-existing run-to-run ordering sensitivity (trace, htr) — for these, even two
    `clavis-cli` runs disagree, so exact equality is not a valid oracle;
  * smoke coverage (runs, returns well-formed JSON) for hbond/occupancy/grade/
    guest_centric, which `clavis-cli` does not expose.

Run:  .venv-py311/bin/python -m pytest pytests/ -q
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

import clavis_engine as ce

ENGINE_DIR = Path(__file__).resolve().parents[1]
CLI = ENGINE_DIR / "target" / "release" / "clavis-cli"
GRO = ENGINE_DIR / "tests" / "fixtures" / "crystals" / "222_S2.gro"

DETERMINISTIC = {"chill": 0.35, "f3f4": 0.35}
CAGE_CUTOFF = 0.36


@pytest.fixture(scope="module")
def content() -> str:
    return GRO.read_text()


def _cli(kind: str, cutoff: float) -> dict:
    res = subprocess.run(
        [str(CLI), kind, str(GRO), str(cutoff)], capture_output=True, text=True
    )
    assert res.returncode == 0, f"clavis-cli {kind} failed: {res.stderr}"
    return json.loads(res.stdout)


def test_version():
    assert ce.version() == "0.3.18"


def test_source_identity_is_exposed_by_the_loaded_native_engine():
    identity = json.loads(ce.source_identity())
    # Release builds embed real digests; local builds report the placeholder from src/py.rs.
    placeholder = "development-unverified"
    commit, tree = identity["sourceCommit"], identity["sourceTreeSha256"]
    assert commit == placeholder or re.fullmatch(r"[0-9a-f]{40}", commit)
    assert tree == placeholder or re.fullmatch(r"[0-9a-f]{64}", tree)


def test_parse_gro_shape(content):
    parsed = json.loads(ce.parse_gro(content))
    assert set(parsed) == {"waters", "guests", "box_dims", "box_vectors"}
    assert len(parsed["box_dims"]) == 3
    assert len(parsed["box_vectors"]) in (3, 9)
    assert len(parsed["waters"]) > 0
    # Every water carries O + 2 H, coordinates in nm.
    w0 = parsed["waters"][0]
    assert set(w0) == {"oxygen", "hydrogen1", "hydrogen2"}
    assert len(w0["oxygen"]) == 3


@pytest.mark.skipif(not CLI.exists(), reason="clavis-cli oracle not built")
@pytest.mark.parametrize("kind", sorted(DETERMINISTIC))
def test_deterministic_exact_parity_with_cli(content, kind):
    """chill / f3f4 are deterministic → the bridge must match clavis-cli byte-for-byte."""
    cutoff = DETERMINISTIC[kind]
    assert json.loads(ce.analyze_gro(content, kind, cutoff)) == _cli(kind, cutoff)


def test_hbond_symmetric_and_self_consistent(content):
    """H-bond adjacency is deterministic (the handoff's correctness oracle)."""
    adj1 = json.loads(ce.analyze_gro(content, "hbond", CAGE_CUTOFF))
    adj2 = json.loads(ce.analyze_gro(content, "hbond", CAGE_CUTOFF))
    assert adj1 == adj2, "hbond network must be deterministic run-to-run"
    # Undirected: j in adj[i]  <=>  i in adj[j].
    for i, nbrs in enumerate(adj1):
        for j in nbrs:
            assert i in adj1[j]


@pytest.mark.parametrize("kind", ["trace", "htr"])
def test_cage_analyses_wellformed(content, kind):
    """trace/htr have pre-existing run-to-run cage-assembly variance, so we assert
    structural sanity (parallel arrays, index bounds) rather than exact equality."""
    n_waters = len(json.loads(ce.parse_gro(content))["waters"])
    out = json.loads(ce.analyze_gro(content, kind, CAGE_CUTOFF))
    if kind == "trace":
        assert len(out["cages_all"]) == len(out["cage_centers"]) == len(out["cage_types"])
        for cage in out["cages_all"]:
            assert all(0 <= v < n_waters for v in cage)
    else:  # htr
        assert {"rings", "cages", "n_water"} <= set(out)
        for ring in out["rings"]:
            assert ring["size"] in (4, 5, 6)
            assert len(ring["vertices"]) == ring["size"]


@pytest.mark.parametrize("kind", ["occupancy", "grade", "guest_centric"])
def test_extra_analyses_smoke(content, kind):
    """hbond/occupancy/grade/guest_centric have no CLI oracle — assert they run and
    return well-formed JSON with the expected top-level keys."""
    out = json.loads(ce.analyze_gro(content, kind, CAGE_CUTOFF))
    if kind == "grade":
        assert {
            "cages_all",
            "cage_centers",
            "cage_type_labels",
            "filled_512",
            "filled_51262",
            "filled_51264",
        } <= set(out)
        assert len(out["cages_all"]) == len(out["cage_centers"]) == len(
            out["cage_type_labels"]
        )
    elif kind == "occupancy":
        assert {"assignments", "overall_occupancy", "total_cages"} <= set(out)
    else:  # guest_centric
        assert {"results", "n_liquid", "n_small", "n_large"} <= set(out)


def test_analyze_waters_matches_analyze_gro_for_deterministic(content):
    """The per-frame coords-in path must agree with the parse+analyze path (chill)."""
    parsed = json.loads(ce.parse_gro(content))
    waters_json = json.dumps(parsed["waters"])
    guests_json = json.dumps(parsed["guests"])
    via_gro = json.loads(ce.analyze_gro(content, "chill", 0.35))
    via_waters = json.loads(
        ce.analyze_waters(waters_json, guests_json, parsed["box_dims"], "chill", 0.35)
    )
    assert via_gro == via_waters


def test_unknown_analysis_raises(content):
    with pytest.raises(ValueError):
        ce.analyze_gro(content, "nonsense", 0.35)


def test_ortho_3v_equals_6v_zero_tilt(content):
    """A 6-value box with zero tilt must reproduce the 3-value orthorhombic result
    exactly — proves the cell path is a faithful superset of the ortho path."""
    parsed = json.loads(ce.parse_gro(content))
    w, g = json.dumps(parsed["waters"]), json.dumps(parsed["guests"])
    lx, ly, lz = parsed["box_dims"]
    via_3 = json.loads(ce.analyze_waters(w, g, [lx, ly, lz], "chill", 0.35))
    via_6 = json.loads(ce.analyze_waters(w, g, [lx, ly, lz, 0.0, 0.0, 0.0], "chill", 0.35))
    assert via_3 == via_6


def test_tilt_changes_result():
    """Decisive minimal case: two oxygens separated across the y-boundary are neighbors
    (O-O 0.2 nm) under an orthorhombic box, but the xy tilt shifts the y-image out to
    ~0.54 nm, breaking the CHILL O-O edge. Proves the tilt is genuinely applied."""
    def water(o):
        return {"oxygen": o, "hydrogen1": [o[0] + 0.01, o[1], o[2]], "hydrogen2": [o[0], o[1] + 0.01, o[2]]}
    w = json.dumps([water([0.5, 0.1, 0.5]), water([0.5, 0.9, 0.5])])
    ortho = json.loads(ce.analyze_waters(w, "[]", [1.0, 1.0, 1.0], "chill", 0.35))
    tilted = json.loads(ce.analyze_waters(w, "[]", [1.0, 1.0, 1.0, 0.5, 0.0, 0.0], "chill", 0.35))
    assert len(ortho["edges"]) == 1, "expected an O-O edge across the y-boundary (ortho)"
    assert len(tilted["edges"]) == 0, "xy tilt should have broken the y-boundary edge"


def test_invalid_box_length_raises():
    with pytest.raises(ValueError):
        ce.analyze_waters("[]", "[]", [1.0, 1.0, 1.0, 0.1, 0.0], "chill", 0.35)  # 5 values
