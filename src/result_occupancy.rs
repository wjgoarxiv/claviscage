// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Result-owned occupancy extraction for the native Python bridge.

use serde::Deserialize;

use crate::geometry::Cell;
use crate::occupancy::{self, OccupancyResult};
use crate::Guest;

#[derive(Deserialize)]
struct TraceResultInput {
    cage_centers: Vec<[f64; 3]>,
    cage_types: Vec<String>,
}

#[derive(Deserialize)]
struct HtrCageInput {
    center: [f64; 3],
    cage_type: String,
}

#[derive(Deserialize)]
struct HtrResultInput {
    cages: Vec<HtrCageInput>,
}

#[derive(Deserialize)]
struct GradeResultInput {
    cage_centers: Vec<[f64; 3]>,
    cage_type_labels: Vec<String>,
}

fn parse_result(
    algorithm: &str,
    result_json: &str,
) -> Result<(Vec<[f64; 3]>, Vec<String>), String> {
    match algorithm {
        "trace" => {
            let result: TraceResultInput = serde_json::from_str(result_json)
                .map_err(|error| format!("invalid TRACE result JSON: {error}"))?;
            Ok((result.cage_centers, result.cage_types))
        }
        "htr" => {
            let result: HtrResultInput = serde_json::from_str(result_json)
                .map_err(|error| format!("invalid HTR result JSON: {error}"))?;
            let mut centers = Vec::with_capacity(result.cages.len());
            let mut cage_types = Vec::with_capacity(result.cages.len());
            for cage in result.cages {
                centers.push(cage.center);
                cage_types.push(normalize_htr_type(&cage.cage_type));
            }
            Ok((centers, cage_types))
        }
        "grade" => {
            let result: GradeResultInput = serde_json::from_str(result_json)
                .map_err(|error| format!("invalid GRADE result JSON: {error}"))?;
            Ok((result.cage_centers, result.cage_type_labels))
        }
        unsupported => Err(format!(
            "unsupported occupancy algorithm '{unsupported}' (expected trace, htr, or grade)"
        )),
    }
}

fn normalize_htr_type(cage_type: &str) -> String {
    match cage_type {
        "0-12-0" => "512".to_owned(),
        "0-12-2" => "51262".to_owned(),
        "0-12-4" => "51264".to_owned(),
        other => other.to_owned(),
    }
}

/// Calculate occupancy from cages already owned by the selected algorithm result.
///
/// # Errors
/// Returns an error when JSON parsing, selected-result shape, geometry, or cell validation fails.
pub fn calculate_from_result(
    algorithm: &str,
    result_json: &str,
    guests_json: &str,
    box_dims: &[f64],
) -> Result<OccupancyResult, String> {
    let (cage_centers, cage_types) = parse_result(algorithm, result_json)?;
    if cage_centers.len() != cage_types.len() {
        return Err(format!(
            "cage center/type arrays must have equal length; got {} centers and {} types",
            cage_centers.len(),
            cage_types.len()
        ));
    }
    if !cage_centers.iter().flatten().all(|value| value.is_finite()) {
        return Err("cage centers must contain only finite coordinates".to_owned());
    }

    let guests: Vec<Guest> = serde_json::from_str(guests_json)
        .map_err(|error| format!("invalid guest JSON: {error}"))?;
    if !guests
        .iter()
        .flat_map(|guest| guest.position)
        .all(f64::is_finite)
    {
        return Err("guest positions must contain only finite coordinates".to_owned());
    }

    let cell = Cell::from_box_values(box_dims)?;
    Ok(occupancy::calculate_occupancy_cell(
        &cage_centers,
        &cage_types,
        &guests,
        &cell,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ORTHO_BOX: &[f64] = &[1.0, 1.0, 1.0];
    const CENTERED_GUEST: &str = r#"[{"name":"CH4","position":[0.5,0.5,0.5]}]"#;

    #[test]
    fn result_occupancy_uses_trace_shape() {
        let result = calculate_from_result(
            "trace",
            r#"{"cage_centers":[[0.5,0.5,0.5]],"cage_types":["512"]}"#,
            CENTERED_GUEST,
            ORTHO_BOX,
        );
        assert_eq!(result.map(|value| value.total_occupied), Ok(1));
    }

    #[test]
    fn result_occupancy_normalizes_htr_types_and_preserves_unknown_types() {
        let result = calculate_from_result(
            "htr",
            r#"{"cages":[{"center":[0.5,0.5,0.5],"cage_type":"0-12-0"},{"center":[0.1,0.1,0.1],"cage_type":"1-2-3"}]}"#,
            CENTERED_GUEST,
            ORTHO_BOX,
        );
        let cage_types = result.map(|value| {
            value
                .by_cage_type
                .into_iter()
                .map(|row| row.cage_type)
                .collect::<Vec<_>>()
        });
        assert_eq!(cage_types, Ok(vec!["1-2-3".to_owned(), "512".to_owned()]));
    }

    #[test]
    fn result_occupancy_uses_grade_shape() {
        let result = calculate_from_result(
            "grade",
            r#"{"cage_centers":[[0.5,0.5,0.5]],"cage_type_labels":["51264"]}"#,
            CENTERED_GUEST,
            ORTHO_BOX,
        );
        assert_eq!(result.map(|value| value.total_occupied), Ok(1));
    }

    #[test]
    fn result_occupancy_rejects_invalid_boundaries() {
        let mismatched = calculate_from_result(
            "trace",
            r#"{"cage_centers":[[0.5,0.5,0.5]],"cage_types":[]}"#,
            "[]",
            ORTHO_BOX,
        );
        let wrong_shape = calculate_from_result("trace", r#"{"cages":[]}"#, "[]", ORTHO_BOX);
        let invalid_guests = calculate_from_result(
            "trace",
            r#"{"cage_centers":[],"cage_types":[]}"#,
            "{",
            ORTHO_BOX,
        );
        let invalid_box = calculate_from_result(
            "trace",
            r#"{"cage_centers":[],"cage_types":[]}"#,
            "[]",
            &[1.0, 1.0],
        );
        assert!(mismatched.is_err());
        assert!(wrong_shape.is_err());
        assert!(invalid_guests.is_err());
        assert!(invalid_box.is_err());
    }

    #[test]
    fn result_occupancy_rejects_unsupported_algorithm_and_non_finite_center() {
        let unsupported = calculate_from_result("other", "{}", "[]", ORTHO_BOX);
        let non_finite = calculate_from_result(
            "trace",
            r#"{"cage_centers":[[1e400,0.0,0.0]],"cage_types":["512"]}"#,
            "[]",
            ORTHO_BOX,
        );
        assert!(unsupported.is_err());
        assert!(non_finite.is_err());
    }

    #[test]
    fn result_occupancy_applies_triclinic_minimum_image() {
        let result = calculate_from_result(
            "trace",
            r#"{"cage_centers":[[0.1,0.1,0.5]],"cage_types":["512"]}"#,
            r#"[{"name":"CH4","position":[0.6,0.9,0.5]}]"#,
            &[1.0, 1.0, 1.0, 0.5, 0.0, 0.0],
        );
        let assignment = result.map(|value| {
            (
                value.total_occupied,
                value.assignments[0].distance_to_center,
            )
        });
        assert_eq!(assignment.as_ref().map(|value| value.0), Ok(1));
        assert!((assignment.as_ref().map(|value| value.1).unwrap_or_default() - 0.2).abs() < 1e-12);
    }
}
