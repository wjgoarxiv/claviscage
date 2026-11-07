// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Native GRADE cage detection.
//!
//! This module reproduces the GRADE v1.00 reference program step by
//! step: an oxygen-distance graph, planar five- and six-membered rings, `5^6`
//! and `6^1 5^6` cups, and cup pairing into `5^12`, `5^12 6^2` and `5^12 6^4`
//! cages. It does not reuse TRACE cage labels or infer an sI/sII/sH verdict.
//!
//! Parity with the reference is the goal, including its defects: planarity is
//! judged on raw coordinates, repeated hexagons can survive, and hexagon and
//! pentagon indices are compared as if they shared one numbering. The
//! per-step notes in the submodules say where each behaviour comes from.

mod cages;
mod cups;
mod faces;
mod metric;
mod planarity;
mod rings;

use crate::geometry::Cell;
use crate::Water;
use cages::{cage_vertices, filled_count, group_cups_64512, pair_cups, reference_center};
use cups::{find_cups_512, find_cups_62512, Topology};
use faces::{link_rings, remove_duplicates};
use metric::Metric;
use rings::find_rings;

const MIN_OO_DISTANCE: f64 = 0.18;
const FILLED_DISTANCE: f64 = 0.20;

/// Numerical parameters exposed by the GRADE v1.00 command line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradeParameters {
    pub rcut: f64,
    pub theta_deg: f64,
    pub pentagon_diagonal_delta: f64,
    pub hexagon_diagonal_delta: f64,
}

impl Default for GradeParameters {
    fn default() -> Self {
        Self {
            rcut: 0.35,
            theta_deg: 45.0,
            pentagon_diagonal_delta: 0.18,
            hexagon_diagonal_delta: 0.26,
        }
    }
}

/// The only cage families defined by GRADE.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GradeCageType {
    Cage512,
    Cage51262,
    Cage51264,
}

impl GradeCageType {
    pub fn from_face_counts(n4: usize, n5: usize, n6: usize) -> Option<Self> {
        match (n4, n5, n6) {
            (0, 12, 0) => Some(Self::Cage512),
            (0, 12, 2) => Some(Self::Cage51262),
            (0, 12, 4) => Some(Self::Cage51264),
            _ => None,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Cage512 => "512",
            Self::Cage51262 => "51262",
            Self::Cage51264 => "51264",
        }
    }
}

/// Canonical cage lists and the corresponding GRADE `.xvg` filled counts.
#[derive(Debug, Clone, PartialEq)]
pub struct NativeGradeResult {
    pub cages_512: Vec<Vec<usize>>,
    pub cages_51262: Vec<Vec<usize>>,
    pub cages_51264: Vec<Vec<usize>>,
    pub centers_512: Vec<[f64; 3]>,
    pub centers_51262: Vec<[f64; 3]>,
    pub centers_51264: Vec<[f64; 3]>,
    pub filled_512: usize,
    pub filled_51262: usize,
    pub filled_51264: usize,
}

impl NativeGradeResult {
    fn empty() -> Self {
        Self {
            cages_512: Vec::new(),
            cages_51262: Vec::new(),
            cages_51264: Vec::new(),
            centers_512: Vec::new(),
            centers_51262: Vec::new(),
            centers_51264: Vec::new(),
            filled_512: 0,
            filled_51262: 0,
            filled_51264: 0,
        }
    }
}

/// Detect native GRADE cages in an orthorhombic box.
///
/// `guest_positions` is intentionally explicit: the filled columns in the
/// reference `.xvg` count unique solute sites within 0.20 nm of a cage center.
pub fn detect_grade(
    waters: &[Water],
    guest_positions: &[[f64; 3]],
    box_dims: &[f64; 3],
    parameters: GradeParameters,
) -> Result<NativeGradeResult, String> {
    let cell = Cell::orthorhombic(*box_dims)?;
    detect_grade_cell(waters, guest_positions, &cell, parameters)
}

/// Detect native GRADE cages using a general periodic cell.
pub fn detect_grade_cell(
    waters: &[Water],
    guest_positions: &[[f64; 3]],
    cell: &Cell,
    parameters: GradeParameters,
) -> Result<NativeGradeResult, String> {
    validate_parameters(parameters)?;
    if waters.len() < 20 {
        return Ok(NativeGradeResult::empty());
    }
    let positions: Vec<[f64; 3]> = waters.iter().map(|water| water.oxygen).collect();
    let metric = Metric::new(cell);
    let found = find_rings(&positions, &metric, parameters);
    Ok(assemble_grade_result(
        &found,
        &positions,
        guest_positions,
        cell,
        &metric,
    ))
}

fn validate_parameters(parameters: GradeParameters) -> Result<(), String> {
    if !parameters.rcut.is_finite() || parameters.rcut <= MIN_OO_DISTANCE {
        return Err("GRADE rcut must be finite and greater than 0.18 nm".to_string());
    }
    if !parameters.theta_deg.is_finite()
        || parameters.theta_deg <= 0.0
        || parameters.theta_deg >= 180.0
    {
        return Err("GRADE planarity angle must be between 0 and 180 degrees".to_string());
    }
    for (name, value) in [
        (
            "pentagon diagonal delta",
            parameters.pentagon_diagonal_delta,
        ),
        ("hexagon diagonal delta", parameters.hexagon_diagonal_delta),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(format!("GRADE {name} must be finite and non-negative"));
        }
    }
    Ok(())
}

/// Cups, cages, display centres and filled counts for one set of rings.
fn assemble_grade_result(
    found: &rings::Rings,
    positions: &[[f64; 3]],
    guests: &[[f64; 3]],
    cell: &Cell,
    metric: &Metric,
) -> NativeGradeResult {
    let links = link_rings(&found.pentagons, &found.hexagons);
    let topology = Topology {
        pentagons: &found.pentagons,
        hexagons: &found.hexagons,
        links: &links,
    };
    let cups_5 = remove_duplicates(find_cups_512(&topology), 0, 0);
    let cups_6 = remove_duplicates(find_cups_62512(&topology), 0, 0);
    let rows_512 = pair_cups(&cups_5, &topology);
    let rows_51262 = pair_cups(&cups_6, &topology);
    let rows_51264 = group_cups_64512(&cups_6);

    let describe = |rows: &[cages::CageRings], geometric_center: bool| {
        let vertices: Vec<Vec<usize>> = rows.iter().map(|r| cage_vertices(r, &topology)).collect();
        let reference: Vec<[f64; 3]> = rows
            .iter()
            .enumerate()
            .map(|(i, r)| reference_center(r, i == 0, &topology, positions))
            .collect();
        let displayed = if geometric_center {
            vertices
                .iter()
                .map(|cage| {
                    let points: Vec<_> = cage.iter().map(|&v| positions[v]).collect();
                    cell.average(&points)
                })
                .collect()
        } else {
            reference.clone()
        };
        (
            vertices,
            displayed,
            filled_count(&reference, guests, metric),
        )
    };
    let (cages_512, centers_512, filled_512) = describe(&rows_512, true);
    let (cages_51262, centers_51262, filled_51262) = describe(&rows_51262, true);
    let (cages_51264, centers_51264, filled_51264) = describe(&rows_51264, false);

    NativeGradeResult {
        cages_512,
        cages_51262,
        cages_51264,
        centers_512,
        centers_51262,
        centers_51264,
        filled_512,
        filled_51262,
        filled_51264,
    }
}

#[cfg(test)]
mod tests;
