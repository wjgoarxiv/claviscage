// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! ClavisCage WASM Engine
//!
//! High-performance hydrate cage analysis running directly in the browser.
//! All computation happens on the user's machine using WebAssembly.

// Make modules public for integration tests
pub mod cage;
pub mod chill;
pub mod grade;
pub mod geometry;
pub mod hbond;
pub mod htr;
pub mod occupancy;
pub mod order_params;
pub mod parser;

// Native CPython bindings for the Blender extension. Only compiled with `--features pyo3`
// on host targets; never part of the wasm build.
#[cfg(feature = "pyo3")]
pub mod py;

#[cfg(feature = "pyo3")]
pub mod result_occupancy;

use geometry::Cell;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// Initialize panic hook for better error messages in console
#[wasm_bindgen(start)]
pub fn init() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

/// Water molecule with oxygen and hydrogen positions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Water {
    pub oxygen: [f64; 3],
    pub hydrogen1: [f64; 3],
    pub hydrogen2: [f64; 3],
}

/// Guest molecule (CO2, CH4, etc.)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Guest {
    pub name: String,
    pub position: [f64; 3],
}

/// Result of F3/F4 analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct F3F4Result {
    pub f3_overall: f64,
    pub f4_overall: f64,
    pub f3_per_atom: Vec<f64>,
    pub f4_per_atom: Vec<f64>,
    pub n_waters: usize,
}

/// Incomplete or non-SEC cage info
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IncompleteCageInfo {
    pub vertices: Vec<usize>,
    pub vertex_positions: Vec<[f64; 3]>,
    pub edges: Vec<[usize; 2]>,
    pub center: [f64; 3],
    pub n_squares: usize, // count of 4-membered rings
    pub n_pentagons: usize,
    pub n_hexagons: usize,
    pub n_faces_total: usize,
    pub quality: u8,            // 0=IC, 1=non-SEC, 2=SEC
    pub estimated_type: String, // cage type label or "unknown"
    pub completeness: f64,      // ratio of actual faces to expected faces (0.0-1.0)
}

/// Result of cage detection
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CageResult {
    pub cages_512: Vec<Vec<usize>>,   // 5^12 cages (small) - legacy
    pub cages_51262: Vec<Vec<usize>>, // 5^12 6^2 cages (medium) - legacy
    pub cages_51264: Vec<Vec<usize>>, // 5^12 6^4 cages (large) - legacy
    pub cages_all: Vec<Vec<usize>>, // all cage vertex lists (parallel with cage_types/cage_centers)
    pub cage_centers: Vec<[f64; 3]>,
    pub cage_types: Vec<String>,
    pub incomplete_cages: Vec<IncompleteCageInfo>,
}

/// HTR ring info
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HtrRingInfo {
    pub size: usize,
    pub vertices: Vec<usize>,
    pub center: Vec<f64>,
    pub vertex_positions: Vec<Vec<f64>>,
    pub edges: Vec<Vec<usize>>,
}

/// HTR cage info
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HtrCageInfo {
    pub cage_type: String, // "n4-n5-n6" format
    pub ring_indices: Vec<usize>,
    pub vertices: Vec<usize>,
    pub center: [f64; 3],
    pub vertex_positions: Vec<Vec<f64>>,
    pub edges: Vec<Vec<usize>>,
    pub n4: usize,
    pub n5: usize,
    pub n6: usize,
}

/// Cage type count
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CageTypeCount {
    pub cage_type: String,
    pub count: usize,
}

/// HTR analysis result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HtrResult {
    pub n_water: usize,
    pub n_4ring: usize,
    pub n_5ring: usize,
    pub n_6ring: usize,
    pub rings: Vec<HtrRingInfo>,
    pub cages: Vec<HtrCageInfo>,
    pub cage_types: Vec<CageTypeCount>,
    pub box_dims: Vec<f64>,
    pub incomplete_cages: Vec<IncompleteCageInfo>,
}

/// Native GRADE cage-detection result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GradeResult {
    pub cages_512: Vec<Vec<usize>>,
    pub cages_51262: Vec<Vec<usize>>,
    pub cages_51264: Vec<Vec<usize>>,
    pub filled_512: usize,
    pub filled_51262: usize,
    pub filled_51264: usize,
    pub cages_all: Vec<Vec<usize>>,
    pub cage_centers: Vec<[f64; 3]>,
    pub cage_type_labels: Vec<String>,
    pub rcut: f64,
    pub planarity_theta: f64,
    pub pentagon_diagonal_delta: f64,
    pub hexagon_diagonal_delta: f64,
}

impl GradeResult {
    pub fn from_native(result: grade::NativeGradeResult, parameters: grade::GradeParameters) -> Self {
        let mut cages_all = Vec::with_capacity(
            result.cages_512.len() + result.cages_51262.len() + result.cages_51264.len(),
        );
        let mut cage_centers = Vec::with_capacity(cages_all.capacity());
        let mut cage_type_labels = Vec::with_capacity(cages_all.capacity());

        cages_all.extend(result.cages_512.iter().cloned());
        cage_centers.extend(result.centers_512.iter().copied());
        cage_type_labels.extend((0..result.cages_512.len()).map(|_| "512".to_string()));

        cages_all.extend(result.cages_51262.iter().cloned());
        cage_centers.extend(result.centers_51262.iter().copied());
        cage_type_labels.extend((0..result.cages_51262.len()).map(|_| "51262".to_string()));

        cages_all.extend(result.cages_51264.iter().cloned());
        cage_centers.extend(result.centers_51264.iter().copied());
        cage_type_labels.extend((0..result.cages_51264.len()).map(|_| "51264".to_string()));

        Self {
            cages_512: result.cages_512,
            cages_51262: result.cages_51262,
            cages_51264: result.cages_51264,
            filled_512: result.filled_512,
            filled_51262: result.filled_51262,
            filled_51264: result.filled_51264,
            cages_all,
            cage_centers,
            cage_type_labels,
            rcut: parameters.rcut,
            planarity_theta: parameters.theta_deg,
            pentagon_diagonal_delta: parameters.pentagon_diagonal_delta,
            hexagon_diagonal_delta: parameters.hexagon_diagonal_delta,
        }
    }
}

/// Parse a GRO file and return waters and guests
#[wasm_bindgen]
pub fn parse_gro(content: &str) -> Result<JsValue, JsValue> {
    let result = parser::parse_gro_content(content).map_err(|e| JsValue::from_str(&e))?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Count the number of frames in a multi-frame GRO file
#[wasm_bindgen]
pub fn count_frames(content: &str) -> usize {
    parser::count_frames(content)
}

/// Parse a specific frame from a multi-frame GRO file
#[wasm_bindgen]
pub fn parse_gro_frame(content: &str, frame_idx: usize) -> Result<JsValue, JsValue> {
    let result = parser::parse_gro_frame(content, frame_idx).map_err(|e| JsValue::from_str(&e))?;

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Calculate F3 and F4 order parameters
///
/// # Arguments
/// * `waters_js` - Array of water molecules (from parse_gro)
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `cutoff` - O-O cutoff distance in nm (default 0.35)
#[wasm_bindgen]
pub fn calculate_f3_f4(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];
    let cutoff = cutoff.unwrap_or(0.35);

    let result = order_params::calculate_f3_f4(&waters, &box_arr, cutoff);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Calculate F3 and F4 order parameters with a 3/6/9-value periodic cell.
///
/// Cell formats: `[lx,ly,lz]`, restricted triclinic `[lx,ly,lz,xy,xz,yz]`,
/// or 9-value GROMACS `.gro` box order.
#[wasm_bindgen]
pub fn calculate_f3_f4_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let cutoff = cutoff.unwrap_or(0.35);

    let result = order_params::calculate_f3_f4_cell(&waters, &cell, cutoff);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// CHILL+ structure identification: classify each water into 6 classes
/// (0 other/liquid, 1 hexagonal ice, 2 cubic ice, 3 interfacial ice,
/// 4 gas hydrate, 5 interfacial hydrate) from the q3 bond-order correlation.
///
/// # Arguments
/// * `waters_js` - Array of water molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `cutoff` - O-O cutoff distance in nm (default 0.35)
#[wasm_bindgen]
pub fn calculate_chill(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cutoff = cutoff.unwrap_or(chill::DEFAULT_CUTOFF);

    let result = chill::calculate_chill(&waters, &box_dims, cutoff);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// CHILL+ structure identification with a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn calculate_chill_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let cutoff = cutoff.unwrap_or(chill::DEFAULT_CUTOFF);

    let result = chill::calculate_chill_cell(&waters, &cell, cutoff);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect hydrate cages using TRACE algorithm
///
/// # Arguments
/// * `waters_js` - Array of water molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `rcut` - O-O cutoff for H-bond (default 0.36 nm)
/// * `theta` - H-bond angle threshold in degrees (default 35.0)
#[wasm_bindgen]
pub fn detect_cages(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let result = cage::detect_cages(&waters, &box_arr, rcut, theta);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect hydrate cages using TRACE algorithm with a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn detect_cages_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let result = cage::detect_cages_cell(&waters, &cell, rcut, theta);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect primitive rings using HTR algorithm
///
/// Returns ring counts (4/5/6-membered), ring geometries, and assembled cages.
///
/// # Arguments
/// * `waters_js` - Array of water molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `rcut` - O-O cutoff for H-bond (default 0.36 nm)
/// * `theta` - H-bond angle threshold in degrees (default 35.0)
#[wasm_bindgen]
pub fn detect_htr(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let result = htr::detect_htr_rings(&waters, &box_arr, rcut, theta);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect primitive rings using HTR algorithm with a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn detect_htr_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let result = htr::detect_htr_rings_cell(&waters, &cell, rcut, theta);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect the three cage families supported by the native GRADE algorithm.
///
/// # Arguments
/// * `waters_js` - Array of water molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `rcut` - O-O cutoff (default 0.35 nm)
/// * `theta` - ring planarity threshold in degrees (default 45.0)
/// * `guests_js` - optional guest positions used for the native filled-cage counts
#[wasm_bindgen]
pub fn classify_grade(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
    guests_js: Option<JsValue>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let guests: Vec<Guest> = match guests_js {
        Some(value) if !value.is_null() && !value.is_undefined() => serde_wasm_bindgen::from_value(value)
            .map_err(|e| JsValue::from_str(&e.to_string()))?,
        _ => Vec::new(),
    };

    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];
    let parameters = grade::GradeParameters {
        rcut: rcut.unwrap_or(0.35),
        theta_deg: theta.unwrap_or(45.0),
        ..grade::GradeParameters::default()
    };
    let guest_positions: Vec<_> = guests.iter().map(|guest| guest.position).collect();
    let native = grade::detect_grade(&waters, &guest_positions, &box_arr, parameters)
        .map_err(|e| JsValue::from_str(&e))?;
    let result = GradeResult::from_native(native, parameters);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Detect native GRADE cages using a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn classify_grade_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
    guests_js: Option<JsValue>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;
    let guests: Vec<Guest> = match guests_js {
        Some(value) if !value.is_null() && !value.is_undefined() => serde_wasm_bindgen::from_value(value)
            .map_err(|e| JsValue::from_str(&e.to_string()))?,
        _ => Vec::new(),
    };

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let parameters = grade::GradeParameters {
        rcut: rcut.unwrap_or(0.35),
        theta_deg: theta.unwrap_or(45.0),
        ..grade::GradeParameters::default()
    };
    let guest_positions: Vec<_> = guests.iter().map(|guest| guest.position).collect();
    let native = grade::detect_grade_cell(&waters, &guest_positions, &cell, parameters)
        .map_err(|e| JsValue::from_str(&e))?;
    let result = GradeResult::from_native(native, parameters);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Build H-bond network and return adjacency list
#[wasm_bindgen]
pub fn build_hbond_network(
    waters_js: JsValue,
    box_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let adj_list = hbond::build_hbond_network(&waters, &box_arr, rcut, theta);

    serde_wasm_bindgen::to_value(&adj_list).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Build H-bond network and return adjacency list using a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn build_hbond_network_cell(
    waters_js: JsValue,
    cell_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> =
        serde_wasm_bindgen::from_value(waters_js).map_err(|e| JsValue::from_str(&e.to_string()))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;
    let rcut = rcut.unwrap_or(0.36);
    let theta = theta.unwrap_or(35.0);

    let adj_list = hbond::build_hbond_network_cell(&waters, &cell, rcut, theta);

    serde_wasm_bindgen::to_value(&adj_list).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Calculate cage occupancy (Mode A: cage-first)
///
/// Given pre-detected cage centers/types and guest positions, assigns each guest
/// to its nearest cage within a type-specific radius cutoff.
///
/// # Arguments
/// * `cage_result_js` - CageResult from detect_cages or classify_grade
/// * `guests_js` - Array of Guest molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
#[wasm_bindgen]
pub fn calculate_occupancy(
    cage_result_js: JsValue,
    guests_js: JsValue,
    box_dims: Vec<f64>,
) -> Result<JsValue, JsValue> {
    let cage_result: CageResult = serde_wasm_bindgen::from_value(cage_result_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse cage_result: {}", e)))?;

    let guests: Vec<Guest> = serde_wasm_bindgen::from_value(guests_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse guests: {}", e)))?;

    if box_dims.len() < 3 {
        return Err(JsValue::from_str("box_dims must have at least 3 elements"));
    }
    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];

    let result = occupancy::calculate_occupancy(
        &cage_result.cage_centers,
        &cage_result.cage_types,
        &guests,
        &box_arr,
    );

    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
}

/// Calculate cage occupancy (Mode A: cage-first) with a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn calculate_occupancy_cell(
    cage_result_js: JsValue,
    guests_js: JsValue,
    cell_dims: Vec<f64>,
) -> Result<JsValue, JsValue> {
    let cage_result: CageResult = serde_wasm_bindgen::from_value(cage_result_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse cage_result: {}", e)))?;

    let guests: Vec<Guest> = serde_wasm_bindgen::from_value(guests_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse guests: {}", e)))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;

    let result = occupancy::calculate_occupancy_cell(
        &cage_result.cage_centers,
        &cage_result.cage_types,
        &guests,
        &cell,
    );

    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
}

/// Calculate guest-centric cage classification (Mode B)
///
/// For each guest molecule, builds a hydration shell, computes local order (F3),
/// counts primitive rings, and classifies the local environment.
///
/// # Arguments
/// * `waters_js` - Array of water molecules
/// * `guests_js` - Array of Guest molecules
/// * `box_dims` - Box dimensions [x, y, z] in nm
/// * `rcut` - H-bond O-O cutoff in nm (default 0.35)
/// * `theta` - H-bond angle threshold in degrees (default 35.0)
/// * `shell_cutoff` - Shell radius in nm (default 0.55)
#[wasm_bindgen]
pub fn calculate_guest_centric(
    waters_js: JsValue,
    guests_js: JsValue,
    box_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
    shell_cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> = serde_wasm_bindgen::from_value(waters_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse waters: {}", e)))?;

    let guests: Vec<Guest> = serde_wasm_bindgen::from_value(guests_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse guests: {}", e)))?;

    if box_dims.len() < 3 {
        return Err(JsValue::from_str("box_dims must have at least 3 elements"));
    }
    let box_arr = [box_dims[0], box_dims[1], box_dims[2]];

    let hbond_rcut = rcut.unwrap_or(0.35);
    let hbond_theta = theta.unwrap_or(35.0);
    let s_cutoff = shell_cutoff.unwrap_or(0.55);

    let result = occupancy::calculate_guest_centric(
        &waters,
        &guests,
        &box_arr,
        hbond_rcut,
        hbond_theta,
        s_cutoff,
    );

    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
}

/// Calculate guest-centric cage classification (Mode B) with a 3/6/9-value periodic cell.
#[wasm_bindgen]
pub fn calculate_guest_centric_cell(
    waters_js: JsValue,
    guests_js: JsValue,
    cell_dims: Vec<f64>,
    rcut: Option<f64>,
    theta: Option<f64>,
    shell_cutoff: Option<f64>,
) -> Result<JsValue, JsValue> {
    let waters: Vec<Water> = serde_wasm_bindgen::from_value(waters_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse waters: {}", e)))?;

    let guests: Vec<Guest> = serde_wasm_bindgen::from_value(guests_js)
        .map_err(|e| JsValue::from_str(&format!("Failed to parse guests: {}", e)))?;

    let cell = Cell::from_box_values(&cell_dims).map_err(|e| JsValue::from_str(&e))?;

    let hbond_rcut = rcut.unwrap_or(0.35);
    let hbond_theta = theta.unwrap_or(35.0);
    let s_cutoff = shell_cutoff.unwrap_or(0.55);

    let result = occupancy::calculate_guest_centric_cell(
        &waters,
        &guests,
        &cell,
        hbond_rcut,
        hbond_theta,
        s_cutoff,
    );

    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsValue::from_str(&format!("Serialization error: {}", e)))
}

/// Get version info
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}
