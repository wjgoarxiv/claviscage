// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Native CPython bindings exposing the ClavisCage analysis engine to Python.
//!
//! This is the integration surface for the Blender extension. It reuses the exact
//! same pure-Rust analysis functions that back `clavis-cli` and the WASM build, so
//! results are identical across all three surfaces. Enabled only under `--features pyo3`
//! on host targets; it is never compiled into the wasm (cdylib) build.
//!
//! Every function returns a JSON string (the same serde serialization `clavis-cli`
//! prints), which the Python caller decodes with `json.loads`. Water/guest inputs are
//! likewise passed as JSON so this module needs no `pyo3` derives on the engine structs.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::geometry::Cell;
use crate::{cage, chill, grade, hbond, htr, occupancy, order_params, parser, result_occupancy};
use crate::{GradeResult, Guest, Water};

/// Analyses supported by [`run`]. Cutoff defaults mirror `clavis-cli` / `lib.rs`.
const ANALYSES: &[&str] = &[
    "chill",
    "trace",
    "htr",
    "hbond",
    "occupancy",
    "guest_centric",
    "f3f4",
    "grade",
];

/// Dispatch one analysis over already-parsed waters/guests and serialize it to JSON.
///
/// `box_dims` may be 3 values (orthorhombic `[lx, ly, lz]`), 6 values (restricted
/// triclinic `[lx, ly, lz, xy, xz, yz]`), or 9 values (GROMACS box order). A 3-value box
/// takes the original orthorhombic code path unchanged; 6/9-value boxes build a general
/// periodic `Cell` and use the engine's `*_cell` variants so tilted cells are handled
/// with correct minimum-image geometry.
fn run(
    waters: &[Water],
    guests: &[Guest],
    box_dims: &[f64],
    kind: &str,
    cutoff: Option<f64>,
    theta: Option<f64>,
) -> Result<String, String> {
    match box_dims.len() {
        3 => run_ortho(
            waters,
            guests,
            [box_dims[0], box_dims[1], box_dims[2]],
            kind,
            cutoff,
            theta,
        ),
        6 | 9 => {
            let cell = Cell::from_box_values(box_dims)?;
            run_cell(waters, guests, &cell, kind, cutoff, theta)
        }
        n => Err(format!("box must have 3, 6, or 9 values; got {n}")),
    }
}

/// Orthorhombic path — unchanged behavior, byte-identical to `clavis-cli`.
fn run_ortho(
    waters: &[Water],
    guests: &[Guest],
    box3: [f64; 3],
    kind: &str,
    cutoff: Option<f64>,
    theta: Option<f64>,
) -> Result<String, String> {
    let hbond_theta = theta.unwrap_or(35.0);

    let json = match kind {
        "chill" => serde_json::to_string(&chill::calculate_chill(
            waters,
            &box3,
            cutoff.unwrap_or(0.35),
        )),
        "f3f4" => serde_json::to_string(&order_params::calculate_f3_f4(
            waters,
            &box3,
            cutoff.unwrap_or(0.35),
        )),
        "trace" => serde_json::to_string(&cage::detect_cages(
            waters,
            &box3,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "htr" => serde_json::to_string(&htr::detect_htr_rings(
            waters,
            &box3,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "hbond" => serde_json::to_string(&hbond::build_hbond_network(
            waters,
            &box3,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "occupancy" => {
            let cages = cage::detect_cages(waters, &box3, cutoff.unwrap_or(0.36), hbond_theta);
            serde_json::to_string(&occupancy::calculate_occupancy(
                &cages.cage_centers,
                &cages.cage_types,
                guests,
                &box3,
            ))
        }
        "guest_centric" => serde_json::to_string(&occupancy::calculate_guest_centric(
            waters,
            guests,
            &box3,
            cutoff.unwrap_or(0.36),
            hbond_theta,
            0.55,
        )),
        "grade" => serde_json::to_string(&run_grade(
            waters,
            guests,
            &box3,
            cutoff.unwrap_or(0.35),
            theta.unwrap_or(45.0),
        )?),
        other => {
            return Err(format!(
                "unknown analysis '{}' (expected one of {})",
                other,
                ANALYSES.join(", ")
            ))
        }
    };
    json.map_err(|e| e.to_string())
}

/// General periodic-cell path — used for 6/9-value (triclinic) boxes. Mirrors
/// [`run_ortho`] but calls the engine's `*_cell` variants for correct tilted-cell MIC.
fn run_cell(
    waters: &[Water],
    guests: &[Guest],
    cell: &Cell,
    kind: &str,
    cutoff: Option<f64>,
    theta: Option<f64>,
) -> Result<String, String> {
    let hbond_theta = theta.unwrap_or(35.0);

    let json = match kind {
        "chill" => serde_json::to_string(&chill::calculate_chill_cell(
            waters,
            cell,
            cutoff.unwrap_or(0.35),
        )),
        "f3f4" => serde_json::to_string(&order_params::calculate_f3_f4_cell(
            waters,
            cell,
            cutoff.unwrap_or(0.35),
        )),
        "trace" => serde_json::to_string(&cage::detect_cages_cell(
            waters,
            cell,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "htr" => serde_json::to_string(&htr::detect_htr_rings_cell(
            waters,
            cell,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "hbond" => serde_json::to_string(&hbond::build_hbond_network_cell(
            waters,
            cell,
            cutoff.unwrap_or(0.36),
            hbond_theta,
        )),
        "occupancy" => {
            let cages = cage::detect_cages_cell(waters, cell, cutoff.unwrap_or(0.36), hbond_theta);
            serde_json::to_string(&occupancy::calculate_occupancy_cell(
                &cages.cage_centers,
                &cages.cage_types,
                guests,
                cell,
            ))
        }
        "guest_centric" => serde_json::to_string(&occupancy::calculate_guest_centric_cell(
            waters,
            guests,
            cell,
            cutoff.unwrap_or(0.36),
            hbond_theta,
            0.55,
        )),
        "grade" => {
            let parameters = grade::GradeParameters {
                rcut: cutoff.unwrap_or(0.35),
                theta_deg: theta.unwrap_or(45.0),
                ..grade::GradeParameters::default()
            };
            let guest_positions: Vec<_> = guests.iter().map(|guest| guest.position).collect();
            let native = grade::detect_grade_cell(waters, &guest_positions, cell, parameters)?;
            serde_json::to_string(&GradeResult::from_native(native, parameters))
        }
        other => {
            return Err(format!(
                "unknown analysis '{}' (expected one of {})",
                other,
                ANALYSES.join(", ")
            ))
        }
    };
    json.map_err(|e| e.to_string())
}

/// Native GRADE cage detection with the original 0.35 nm / 45 degree defaults.
fn run_grade(
    waters: &[Water],
    guests: &[Guest],
    box3: &[f64; 3],
    rcut: f64,
    theta: f64,
) -> Result<GradeResult, String> {
    let parameters = grade::GradeParameters {
        rcut,
        theta_deg: theta,
        ..grade::GradeParameters::default()
    };
    let guest_positions: Vec<_> = guests.iter().map(|guest| guest.position).collect();
    let native = grade::detect_grade(waters, &guest_positions, box3, parameters)?;
    Ok(GradeResult::from_native(native, parameters))
}

/// Engine version (matches the crate version and the WASM `version()`).
#[pyfunction]
fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// Exact source identity injected by the isolated wheel build.
#[pyfunction]
fn source_identity() -> String {
    serde_json::json!({
        "sourceCommit": option_env!("CLAVIS_SOURCE_COMMIT").unwrap_or("development-unverified"),
        "sourceTreeSha256": option_env!("CLAVIS_SOURCE_TREE_SHA256").unwrap_or("development-unverified"),
    })
    .to_string()
}

/// Parse GRO text into `{waters, guests, box_dims}` (JSON string). Coordinates are nm.
#[pyfunction]
fn parse_gro(content: &str) -> PyResult<String> {
    let parsed = parser::parse_gro_content(content).map_err(PyValueError::new_err)?;
    serde_json::to_string(&parsed).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse GRO text and run `kind`, returning the analysis as a JSON string.
/// This mirrors `clavis-cli <kind> <file.gro>` exactly.
#[pyfunction]
#[pyo3(signature = (content, kind, cutoff=None, theta=None))]
fn analyze_gro(
    content: &str,
    kind: &str,
    cutoff: Option<f64>,
    theta: Option<f64>,
) -> PyResult<String> {
    let parsed = parser::parse_gro_content(content).map_err(PyValueError::new_err)?;
    run(
        &parsed.waters,
        &parsed.guests,
        &parsed.box_vectors,
        kind,
        cutoff,
        theta,
    )
    .map_err(PyValueError::new_err)
}

/// Run `kind` over pre-parsed waters/guests (JSON strings) — the per-frame trajectory
/// path that skips re-parsing GRO text every frame. `waters_json` / `guests_json` are
/// JSON arrays matching the engine's `Water` / `Guest` shapes; `box_dims` is `[lx,ly,lz]`.
#[pyfunction]
#[pyo3(signature = (waters_json, guests_json, box_dims, kind, cutoff=None, theta=None))]
fn analyze_waters(
    waters_json: &str,
    guests_json: &str,
    box_dims: Vec<f64>,
    kind: &str,
    cutoff: Option<f64>,
    theta: Option<f64>,
) -> PyResult<String> {
    let waters: Vec<Water> =
        serde_json::from_str(waters_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let guests: Vec<Guest> =
        serde_json::from_str(guests_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    run(&waters, &guests, &box_dims, kind, cutoff, theta).map_err(PyValueError::new_err)
}

/// Calculate occupancy from an already-computed TRACE, HTR, or GRADE result.
#[pyfunction]
fn calculate_occupancy_from_result(
    algorithm: &str,
    result_json: &str,
    guests_json: &str,
    box_dims: Vec<f64>,
) -> PyResult<String> {
    let result =
        result_occupancy::calculate_from_result(algorithm, result_json, guests_json, &box_dims)
            .map_err(PyValueError::new_err)?;
    serde_json::to_string(&result).map_err(|error| PyValueError::new_err(error.to_string()))
}

#[pymodule]
fn clavis_engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(source_identity, m)?)?;
    m.add_function(wrap_pyfunction!(parse_gro, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_gro, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_waters, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_occupancy_from_result, m)?)?;
    Ok(())
}
