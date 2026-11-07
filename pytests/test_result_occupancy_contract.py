"""Failing-first contracts for result-owned PyO3 cage occupancy."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import clavis_engine as ce

ENGINE_DIR = Path(__file__).resolve().parents[1]
GRO = ENGINE_DIR / "tests" / "fixtures" / "crystals" / "222_S2.gro"
CAGE_CUTOFF = 0.36


@pytest.fixture(scope="module")
def content() -> str:
    return GRO.read_text()


@pytest.mark.parametrize("surface", ["analyze_gro", "analyze_waters"])
def test_legacy_occupancy_is_trace_derived(content, surface):
    """Legacy occupancy keeps its historical TRACE-owned cage selection."""
    # Given a non-empty sII fixture whose TRACE and HTR cage labels are distinct.
    parsed = json.loads(ce.parse_gro(content))
    prefixes = {
        "analyze_gro": (content,),
        "analyze_waters": (
            json.dumps(parsed["waters"]),
            json.dumps(parsed["guests"]),
            parsed["box_dims"],
        ),
    }
    analyze = getattr(ce, surface)

    # When TRACE, HTR, and legacy occupancy run through the same bridge surface.
    trace = json.loads(analyze(*prefixes[surface], "trace", CAGE_CUTOFF))
    htr = json.loads(analyze(*prefixes[surface], "htr", CAGE_CUTOFF))
    occupancy = json.loads(analyze(*prefixes[surface], "occupancy", CAGE_CUTOFF))

    # Then occupancy exposes TRACE labels, not HTR's generalized signatures.
    occupancy_types = {row["cage_type"] for row in occupancy["by_cage_type"]}
    trace_types = set(trace["cage_types"])
    htr_types = {row["cage_type"] for row in htr["cage_types"]}
    assert occupancy["total_cages"] > 0
    assert occupancy_types == trace_types
    assert occupancy_types.isdisjoint(htr_types)


@pytest.mark.parametrize(
    ("algorithm", "analysis_result"),
    [
        ("trace", {"cage_centers": [[0.5, 0.5, 0.5]], "cage_types": ["512"]}),
        ("htr", {"cages": [{"center": [0.5, 0.5, 0.5], "cage_type": "0-12-0"}]}),
        ("grade", {"cage_centers": [[0.5, 0.5, 0.5]], "cage_type_labels": ["512"]}),
    ],
)
def test_result_owned_occupancy_uses_selected_algorithm_shape(algorithm, analysis_result):
    # Given one complete cage in the selected algorithm's native result shape.
    guests_json = json.dumps([{"name": "CH4", "position": [0.5, 0.5, 0.5]}])

    # When occupancy consumes that already-computed result.
    out = json.loads(
        ce.calculate_occupancy_from_result(
            algorithm, json.dumps(analysis_result), guests_json, [1.0, 1.0, 1.0]
        )
    )

    # Then the guest occupies that cage and HTR's canonical signature is normalized.
    assert (out["total_cages"], out["total_occupied"]) == (1, 1)
    assert out["by_cage_type"][0]["cage_type"] == "512"


def test_result_owned_occupancy_rejects_parallel_array_length_mismatch():
    # Given a TRACE result whose centers and types are not parallel arrays.
    result_json = json.dumps({"cage_centers": [[0.5, 0.5, 0.5]], "cage_types": []})

    # When the result boundary parses the malformed scientific contract.
    # Then it raises a Python value error naming the parallel-array mismatch.
    with pytest.raises(ValueError, match="center.*type|type.*center|length|parallel"):
        ce.calculate_occupancy_from_result("trace", result_json, "[]", [1.0, 1.0, 1.0])


def test_result_owned_occupancy_rejects_malformed_json():
    # Given malformed analysis-result JSON.
    # When the result boundary parses it.
    # Then it reports a Python value error.
    with pytest.raises(ValueError):
        ce.calculate_occupancy_from_result("trace", "{", "[]", [1.0, 1.0, 1.0])


def test_result_owned_occupancy_rejects_unsupported_algorithm():
    # Given a syntactically valid result with an unsupported algorithm discriminator.
    result_json = json.dumps({"cage_centers": [], "cage_types": []})

    # When dispatch selects that algorithm.
    # Then it fails instead of silently falling back to TRACE.
    with pytest.raises(ValueError, match="algorithm|unsupported|nonsense"):
        ce.calculate_occupancy_from_result("nonsense", result_json, "[]", [1.0, 1.0, 1.0])


def test_result_owned_occupancy_rejects_shape_from_another_algorithm():
    # Given an HTR result mislabeled as TRACE.
    htr_json = json.dumps({"cages": [{"center": [0.5, 0.5, 0.5], "cage_type": "0-12-0"}]})

    # When TRACE parsing consumes the HTR shape.
    # Then selected-algorithm parsing rejects it.
    with pytest.raises(ValueError):
        ce.calculate_occupancy_from_result("trace", htr_json, "[]", [1.0, 1.0, 1.0])


def test_result_owned_occupancy_applies_triclinic_minimum_image():
    # Given a tilted cell where the nearest image is 0.2 nm from a 5^12 center.
    result_json = json.dumps({"cage_centers": [[0.1, 0.1, 0.5]], "cage_types": ["512"]})
    guests_json = json.dumps([{"name": "CH4", "position": [0.6, 0.9, 0.5]}])

    # When occupancy uses the restricted-triclinic six-value box.
    out = json.loads(
        ce.calculate_occupancy_from_result(
            "trace", result_json, guests_json, [1.0, 1.0, 1.0, 0.5, 0.0, 0.0]
        )
    )

    # Then MIC assigns the guest to the cage across the tilted boundary.
    assert out["total_occupied"] == 1
    assert out["assignments"][0]["distance_to_center"] == pytest.approx(0.2)
