"""Steinhardt q-bar and MCG through the Python bindings, and the method registry.

The expected numbers are the ones the Rust unit tests pin (src/steinhardt.rs,
src/mcg.rs), so a Python/Rust contract drift fails here.
"""
from __future__ import annotations

import json
from pathlib import Path

import clavis_engine as ce

FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures"


def _analyze(rel: str, kind: str) -> dict:
    return json.loads(ce.analyze_gro((FIXTURES / rel).read_text(), kind))


def test_steinhardt_matches_freud_reference():
    out = _analyze("crystals/222_S2.gro", "steinhardt")
    ref = json.loads((FIXTURES / "reference" / "steinhardt_222_S2.json").read_text())
    assert out["ls"] == [4, 6]
    assert out["n_waters"] == ref["n_waters"] == 1088
    for k, l in enumerate((4, 6)):
        worst = max(abs(a - b) for a, b in zip(out["qbar"][k], ref[f"qbar{l}"]))
        assert worst < 1e-4


def test_mcg_full_crystals_form_one_cluster():
    for rel, n in (("guests/sI_CH4_222.gro", 64), ("guests/sII_CH4_222.gro", 192)):
        out = _analyze(rel, "mcg")
        assert out["n_guests"] == n
        assert out["largest_cluster"] == n
        assert out["cluster_sizes"] == [n]
        assert out["parameters"] == {
            "r_guest_nm": 0.9, "r_water_nm": 0.6, "phi_deg": 45.0, "n_water": 5, "n_count": 3,
        }


def test_list_methods_tiers():
    methods = json.loads(ce.list_methods())
    tiers = {m["id"]: m["tier"] for m in methods}
    assert tiers == {
        "trace": "author_code", "htr": "author_code", "grade": "author_code",
        "chill": "paper", "f3f4": "paper",
        "steinhardt": "reference_implementation", "mcg": "paper", "tetra": "reference_implementation",
    }
    assert all(m["doi"].startswith("10.") for m in methods)


def test_descriptors_match_the_shared_snapshot():
    # The same snapshot is checked against the WASM build by its own test suite.
    snapshot = json.loads((FIXTURES / "methods-descriptor.json").read_text())
    assert json.loads(ce.list_methods()) == snapshot
    for m in snapshot:
        assert json.loads(ce.describe_method(m["id"])) == m


def _close(a, b, tol=1e-12) -> bool:
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_close(a[k], b[k], tol) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(_close(x, y, tol) for x, y in zip(a, b))
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_run_method_matches_the_wasm_parity_file():
    # Written from the WASM build by its run_method parity test. Native and wasm32 use
    # different libm builds, so values may differ in the last bit (measured: 1 ulp on F4);
    # keys, parameters and numbers to 1e-12 relative must agree.
    parity = json.loads((FIXTURES / "run-method-parity.json").read_text())
    cases = {
        "example_H2O_frame0": "crystals/example_H2O_frame0.gro",
        "sI_CH4_222": "guests/sI_CH4_222.gro",
        "sH_222": "benchmark/sH_222.gro",  # tilted (9-value) cell
    }
    assert set(parity) == set(cases)
    for case, rel in cases.items():
        text = (FIXTURES / rel).read_text()
        parsed = json.loads(ce.parse_gro(text))
        # Every number on the box line: parse_gro keeps only the first three.
        cell = [float(v) for v in text.rstrip().splitlines()[-1].split()]
        for method, expected in parity[case].items():
            out = json.loads(ce.run_method(
                method, json.dumps(parsed["waters"]), json.dumps(parsed["guests"]), cell, None, True,
            ))
            assert out["id"] == method
            got = {"series": out["series"], "molecule": out["molecule"], "params": out["params"]}
            assert got["params"] == expected["params"], (case, method)
            assert _close(got, expected), (case, method)


def test_run_method_errors_carry_codes():
    import pytest

    parsed = json.loads(ce.parse_gro((FIXTURES / "crystals/example_H2O_frame0.gro").read_text()))
    args = (json.dumps(parsed["waters"]), json.dumps(parsed["guests"]), parsed["box_dims"])
    with pytest.raises(ValueError, match="^unknown_method"):
        ce.run_method("nope", *args)
    with pytest.raises(ValueError, match="^missing_input: guests"):
        ce.run_method("mcg", *args)
    with pytest.raises(ValueError, match="^param_out_of_range"):
        ce.run_method("f3f4", *args, json.dumps({"cutoff_nm": 2.0}))


def _gmx_value(name: str) -> float:
    lines = (FIXTURES / "reference" / "tetra" / name).read_text().splitlines()
    return float([l for l in lines if l and l[0] not in "#@"][-1].split()[1])


def test_tetra_matches_gmx_order():
    # Same reference as src/tetra.rs: GROMACS 2026.2 `gmx order -Sg -Sk`, single precision, six decimals.
    for rel, stem in (("crystals/example_H2O_frame0.gro", "example_H2O_frame0"), ("crystals/222_S2.gro", "222_S2")):
        parsed = json.loads(ce.parse_gro((FIXTURES / rel).read_text()))
        out = json.loads(ce.run_method(
            "tetra", json.dumps(parsed["waters"]), json.dumps(parsed["guests"]), parsed["box_dims"], None, True,
        ))
        sg, sk = _gmx_value(f"{stem}_sg.xvg"), _gmx_value(f"{stem}_sk.xvg")
        assert abs(out["series"]["sg"] - sg) < 2e-6, rel
        assert abs(out["series"]["sk"] - sk) < 2e-6, rel
        assert abs(out["series"]["q"] - (1 - 4 * out["series"]["sg"])) < 1e-12
        assert len(out["molecule"]["q"]) == len(parsed["waters"])


def test_describe_unknown_method_raises():
    import pytest

    with pytest.raises(ValueError, match="unknown method"):
        ce.describe_method("nope")
