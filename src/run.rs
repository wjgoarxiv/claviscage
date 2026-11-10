// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! One entry point for every method that returns per-frame numbers: the app runs any such method
//! by id with parameters named in its descriptor (`methods.rs`), without method-specific code.
//! Each call returns the descriptor's series, optionally its per-molecule values, and the method's
//! own result unchanged. Cage detectors keep their own entry points.

use std::collections::BTreeMap;
use std::fmt;

use serde::Serialize;

use crate::geometry::Cell;
use crate::methods::{self, MethodInfo, View};
use crate::{chill, mcg, order_params, steinhardt, tetra, F3F4Result, Guest, Water};

#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "code", content = "detail", rename_all = "snake_case")]
pub enum MethodError {
    UnknownMethod(String),
    /// The method has no generic runner (cage detectors).
    NotSeriesMethod(String),
    UnknownParam(String),
    ParamOutOfRange { name: String, value: f64, min: f64, max: f64 },
    NotInteger { name: String, value: f64 },
    /// The structure lacks something the method needs (e.g. guests for MCG).
    MissingInput(String),
    Invalid(String),
}

impl fmt::Display for MethodError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownMethod(id) => write!(f, "unknown_method: {id}"),
            Self::NotSeriesMethod(id) => write!(f, "not_series_method: {id} has no generic runner"),
            Self::UnknownParam(n) => write!(f, "unknown_param: {n}"),
            Self::ParamOutOfRange { name, value, min, max } => {
                write!(f, "param_out_of_range: {name} = {value} (allowed {min}..{max})")
            }
            Self::NotInteger { name, value } => write!(f, "not_integer: {name} = {value}"),
            Self::MissingInput(what) => write!(f, "missing_input: {what}"),
            Self::Invalid(msg) => write!(f, "invalid: {msg}"),
        }
    }
}

/// The method's own result type, serialised exactly as its dedicated entry point does.
#[derive(Debug, Clone, Serialize)]
#[serde(untagged)]
pub enum MethodResult {
    F3F4(F3F4Result),
    Steinhardt(steinhardt::SteinhardtResult),
    Mcg(mcg::McgResult),
    Chill(chill::ChillResult),
    Tetra(tetra::TetraResult),
}

#[derive(Debug, Clone, Serialize)]
pub struct MethodOutput {
    pub id: &'static str,
    pub descriptor_version: u32,
    /// The parameters used, defaults filled in.
    pub params: BTreeMap<String, f64>,
    /// One number per descriptor `series` key.
    pub series: BTreeMap<String, f64>,
    /// One array per descriptor `molecule_values` key; empty unless requested.
    pub molecule: BTreeMap<String, Vec<f64>>,
    pub result: MethodResult,
}

/// Fill defaults and check every supplied parameter against the descriptor.
fn resolve_params(m: &MethodInfo, supplied: &BTreeMap<String, f64>) -> Result<BTreeMap<String, f64>, MethodError> {
    if let Some(name) = supplied.keys().find(|k| !m.params.iter().any(|p| p.name == k.as_str())) {
        return Err(MethodError::UnknownParam(name.clone()));
    }
    let mut out = BTreeMap::new();
    for p in m.params {
        let value = supplied.get(p.name).copied().unwrap_or(p.default);
        if !value.is_finite() || value < p.min || value > p.max {
            return Err(MethodError::ParamOutOfRange { name: p.name.into(), value, min: p.min, max: p.max });
        }
        if p.integer && value.fract() != 0.0 {
            return Err(MethodError::NotInteger { name: p.name.into(), value });
        }
        out.insert(p.name.to_string(), value);
    }
    Ok(out)
}

fn keyed(keys: &[methods::SeriesSpec], values: impl IntoIterator<Item = f64>) -> BTreeMap<String, f64> {
    keys.iter().map(|k| k.key.to_string()).zip(values).collect()
}

/// Run one series method on one frame. `cell_values` takes the 3/6/9-value box formats.
pub fn run_method(
    id: &str,
    waters: &[Water],
    guests: &[Guest],
    cell_values: &[f64],
    supplied: &BTreeMap<String, f64>,
    with_molecules: bool,
) -> Result<MethodOutput, MethodError> {
    let m = methods::describe_method(id).ok_or_else(|| MethodError::UnknownMethod(id.into()))?;
    if m.views.contains(&View::CageView) {
        return Err(MethodError::NotSeriesMethod(id.into()));
    }
    let params = resolve_params(&m, supplied)?;
    // The parser keeps only waters with two hydrogens, so an oxygen-only file arrives empty.
    if waters.is_empty() {
        return Err(MethodError::MissingInput("waters".into()));
    }
    let finite = |p: &[f64; 3]| p.iter().all(|x| x.is_finite());
    if !waters.iter().all(|w| finite(&w.oxygen) && finite(&w.hydrogen1) && finite(&w.hydrogen2))
        || !guests.iter().all(|g| finite(&g.position))
    {
        return Err(MethodError::Invalid("a coordinate is not a finite number".into()));
    }
    if m.inputs.guests && guests.is_empty() {
        return Err(MethodError::MissingInput("guests".into()));
    }
    let cell = Cell::from_box_values(cell_values).map_err(MethodError::Invalid)?;
    let p = |name: &str| params[name];

    let (series, molecule, result) = match m.id {
        "f3f4" => {
            let r = order_params::calculate_f3_f4_cell(waters, &cell, p("cutoff_nm"));
            let molecule = [r.f3_per_atom.clone(), r.f4_per_atom.clone()];
            (keyed(m.series, [r.f3_overall, r.f4_overall]), molecule.to_vec(), MethodResult::F3F4(r))
        }
        "steinhardt" => {
            let ls = [4, 6];
            steinhardt::check_inputs(&cell, p("cutoff_nm"), &ls).map_err(MethodError::Invalid)?;
            let r = steinhardt::calculate_steinhardt_cell(waters, &cell, p("cutoff_nm"), &ls);
            let molecule = r.qbar.clone();
            (keyed(m.series, r.qbar_mean.iter().copied()), molecule, MethodResult::Steinhardt(r))
        }
        "mcg" => {
            let mp = mcg::McgParameters {
                r_guest_nm: p("r_guest_nm"),
                r_water_nm: p("r_water_nm"),
                phi_deg: p("phi_deg"),
                n_water: p("n_water") as usize,
                n_count: p("n_count") as usize,
            };
            mp.check(&cell).map_err(MethodError::Invalid)?;
            let r = mcg::calculate_mcg_cell(waters, guests, &cell, mp);
            (keyed(m.series, [r.largest_cluster as f64]), Vec::new(), MethodResult::Mcg(r))
        }
        "chill" => {
            let r = chill::calculate_chill_cell(waters, &cell, p("cutoff_nm"));
            let molecule = vec![r.per_atom_class.iter().map(|&c| f64::from(c)).collect()];
            (keyed(m.series, r.fractions), molecule, MethodResult::Chill(r))
        }
        "tetra" => {
            let r = tetra::calculate_tetra_cell(waters, &cell).map_err(MethodError::Invalid)?;
            let molecule = vec![r.q_per_water.clone()];
            (keyed(m.series, [r.q_mean, r.sg_mean, r.sk_mean]), molecule, MethodResult::Tetra(r))
        }
        other => return Err(MethodError::NotSeriesMethod(other.into())),
    };
    let molecule = if with_molecules {
        m.molecule_values.iter().map(|k| k.key.to_string()).zip(molecule).collect()
    } else {
        BTreeMap::new()
    };
    Ok(MethodOutput { id: m.id, descriptor_version: m.descriptor_version, params, series, molecule, result })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_gro_content;

    fn fixture(rel: &str) -> crate::parser::ParseResult {
        let text = std::fs::read_to_string(format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), rel)).unwrap();
        parse_gro_content(&text).unwrap()
    }

    fn none() -> BTreeMap<String, f64> {
        BTreeMap::new()
    }

    #[test]
    fn generic_runs_equal_the_direct_functions() {
        for rel in ["crystals/example_H2O_frame0.gro", "crystals/222_S2.gro", "guests/sI_CH4_222.gro"] {
            let d = fixture(rel);
            let cell = Cell::from_box_values(&d.box_dims).unwrap();

            let out = run_method("f3f4", &d.waters, &d.guests, &d.box_dims, &none(), true).unwrap();
            let direct = order_params::calculate_f3_f4_cell(&d.waters, &cell, 0.35);
            let MethodResult::F3F4(r) = &out.result else { panic!() };
            assert_eq!(r.f3_per_atom, direct.f3_per_atom, "{rel}");
            assert_eq!(r.f4_per_atom, direct.f4_per_atom, "{rel}");
            assert_eq!(out.series["f3"], direct.f3_overall);
            assert_eq!(out.series["f4"], direct.f4_overall);
            assert_eq!(out.molecule["f4"], direct.f4_per_atom);

            let out = run_method("steinhardt", &d.waters, &d.guests, &d.box_dims, &none(), true).unwrap();
            let direct = steinhardt::calculate_steinhardt_cell(&d.waters, &cell, 0.35, &[4, 6]);
            let MethodResult::Steinhardt(r) = &out.result else { panic!() };
            assert_eq!(r.qbar, direct.qbar, "{rel}");
            assert_eq!(out.series["qbar6"], direct.qbar_mean[1]);
            assert_eq!(out.molecule["qbar4"], direct.qbar[0]);

            let out = run_method("chill", &d.waters, &d.guests, &d.box_dims, &none(), false).unwrap();
            let direct = chill::calculate_chill_cell(&d.waters, &cell, 0.35);
            let MethodResult::Chill(r) = &out.result else { panic!() };
            assert_eq!(r.per_atom_class, direct.per_atom_class, "{rel}");
            assert_eq!(out.series["hydrate"], direct.fractions[4]);
            assert!(out.molecule.is_empty());
        }
        let d = fixture("crystals/222_S2.gro");
        let out = run_method("tetra", &d.waters, &d.guests, &d.box_dims, &none(), true).unwrap();
        let direct = tetra::calculate_tetra_cell(&d.waters, &Cell::from_box_values(&d.box_dims).unwrap()).unwrap();
        assert_eq!(out.series["q"], direct.q_mean);
        assert_eq!(out.series["sk"], direct.sk_mean);
        assert_eq!(out.molecule["q"], direct.q_per_water);
        let d = fixture("guests/sI_CH4_222.gro");
        let out = run_method("mcg", &d.waters, &d.guests, &d.box_dims, &none(), true).unwrap();
        assert_eq!(out.series["largest_cluster"], 64.0);
        assert!(out.molecule.is_empty());
    }

    #[test]
    fn params_are_resolved_and_checked() {
        let d = fixture("guests/sI_CH4_222.gro");
        let run = |id: &str, pairs: &[(&str, f64)]| {
            let p = pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect();
            run_method(id, &d.waters, &d.guests, &d.box_dims, &p, false)
        };
        let out = run("f3f4", &[("cutoff_nm", 0.32)]).unwrap();
        assert_eq!(out.params["cutoff_nm"], 0.32);
        assert!(matches!(run("f3f4", &[("cutoff", 0.3)]), Err(MethodError::UnknownParam(_))));
        assert!(matches!(run("f3f4", &[("cutoff_nm", 0.9)]), Err(MethodError::ParamOutOfRange { .. })));
        assert!(matches!(run("f3f4", &[("cutoff_nm", f64::NAN)]), Err(MethodError::ParamOutOfRange { .. })));
        assert!(matches!(run("mcg", &[("n_count", 2.5)]), Err(MethodError::NotInteger { .. })));
        let out = run("mcg", &[("n_count", 4.0)]).unwrap();
        let MethodResult::Mcg(r) = &out.result else { panic!() };
        assert_eq!(r.parameters.n_count, 4);
    }

    #[test]
    fn typed_errors_instead_of_panics() {
        let d = fixture("crystals/example_H2O_frame0.gro");
        let run = |id: &str, cell: &[f64]| run_method(id, &d.waters, &d.guests, cell, &none(), false);
        assert_eq!(run("nope", &d.box_dims).unwrap_err(), MethodError::UnknownMethod("nope".into()));
        assert!(matches!(run("trace", &d.box_dims), Err(MethodError::NotSeriesMethod(_))));
        assert_eq!(run("mcg", &d.box_dims).unwrap_err(), MethodError::MissingInput("guests".into()));
        assert!(matches!(run("f3f4", &[1.0, 2.0]), Err(MethodError::Invalid(_))));
        assert!(matches!(run("f3f4", &[0.0, 0.0, 0.0]), Err(MethodError::Invalid(_))));

        for id in ["f3f4", "steinhardt", "chill", "tetra"] {
            let empty = run_method(id, &[], &d.guests, &d.box_dims, &none(), false);
            assert_eq!(empty.unwrap_err(), MethodError::MissingInput("waters".into()), "{id}");
        }
        let mut bad = d.waters.clone();
        bad[3].oxygen[1] = f64::INFINITY;
        for id in ["f3f4", "steinhardt", "chill", "tetra"] {
            let out = run_method(id, &bad, &d.guests, &d.box_dims, &none(), false);
            assert!(matches!(out, Err(MethodError::Invalid(_))), "{id}");
        }
    }
}
