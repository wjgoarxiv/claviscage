// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Steinhardt bond-orientational order q_l and its locally averaged form q̄_l.
//!
//! References:
//! - P. J. Steinhardt, D. R. Nelson, M. Ronchetti, Phys. Rev. B 28 (1983) 784.
//! - W. Lechner, C. Dellago, J. Chem. Phys. 129 (2008) 114707,
//!   doi:10.1063/1.2977970 (averaged form, Eqs. 1-3).
//!
//! For water oxygen i with neighbour set N(i) (other oxygens closer than `cutoff`
//! under the minimum-image convention):
//!
//! ```text
//! q_lm(i)  = 1/|N(i)|        Σ_{j∈N(i)}      Y_lm(r̂_ij)
//! q̄_lm(i)  = 1/(|N(i)| + 1)  Σ_{k∈{i}∪N(i)}  q_lm(k)
//! q_l(i)   = sqrt(4π/(2l+1) Σ_m |q_lm(i)|²),   q̄_l(i) likewise from q̄_lm(i)
//! ```
//!
//! A water with no neighbours has q_lm = 0. Results are checked against freud's
//! `Steinhardt(average=True)` on the fixtures in `tests/fixtures/reference/steinhardt_*.json`.

use crate::geometry::Cell;
use crate::Water;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Default O–O neighbour cutoff (nm), the same first-shell cutoff used by F3/F4.
pub const DEFAULT_CUTOFF_NM: f64 = 0.35;

/// Per-water and frame-mean Steinhardt parameters for each requested l.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SteinhardtResult {
    /// The degrees l, in request order; every other field is indexed the same way.
    pub ls: Vec<u32>,
    pub cutoff_nm: f64,
    pub n_waters: usize,
    /// `q[k][i]` is q_l for l = `ls[k]` and water i.
    pub q: Vec<Vec<f64>>,
    /// `qbar[k][i]` is q̄_l for l = `ls[k]` and water i.
    pub qbar: Vec<Vec<f64>>,
    pub q_mean: Vec<f64>,
    pub qbar_mean: Vec<f64>,
}

/// Largest degree accepted; higher l gains nothing for water and overflows the recurrence.
pub const MAX_L: u32 = 20;

/// Rejects inputs the calculation cannot answer correctly: a non-positive cutoff, a degree
/// outside 1..=MAX_L, or a cell narrower than twice the cutoff (neighbours would be counted
/// through only one of several periodic images).
pub fn check_inputs(cell: &Cell, cutoff: f64, ls: &[u32]) -> Result<(), String> {
    if !(cutoff.is_finite() && cutoff > 0.0) {
        return Err(format!("cutoff must be a positive length in nm, got {cutoff}"));
    }
    if let Some(l) = ls.iter().find(|&&l| l == 0 || l > MAX_L) {
        return Err(format!("degree l must be within 1..={MAX_L}, got {l}"));
    }
    let width = cell.min_width();
    if width < 2.0 * cutoff {
        return Err(format!(
            "cell is {width:.3} nm across; Steinhardt needs at least twice the cutoff ({:.3} nm)",
            2.0 * cutoff
        ));
    }
    Ok(())
}

/// Orthorhombic-box entry point.
pub fn calculate_steinhardt(
    waters: &[Water],
    box_dims: &[f64; 3],
    cutoff: f64,
    ls: &[u32],
) -> SteinhardtResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    calculate_steinhardt_cell(waters, &cell, cutoff, ls)
}

/// General periodic-cell entry point.
pub fn calculate_steinhardt_cell(
    waters: &[Water],
    cell: &Cell,
    cutoff: f64,
    ls: &[u32],
) -> SteinhardtResult {
    let n = waters.len();
    let neighbors = neighbor_list(waters, cell, cutoff);
    let mut q = Vec::with_capacity(ls.len());
    let mut qbar = Vec::with_capacity(ls.len());

    for &l in ls {
        let l = l as usize;
        // qlm[i][m] for m = 0..=l; negative m follow from Y_l,-m = (-1)^m conj(Y_lm),
        // so |q_l,-m| = |q_lm| and the sum over m doubles the m > 0 terms.
        let mut qlm = vec![vec![(0.0f64, 0.0f64); l + 1]; n];
        for i in 0..n {
            if neighbors[i].is_empty() {
                continue;
            }
            for &j in &neighbors[i] {
                let d = cell.diff(&waters[j].oxygen, &waters[i].oxygen);
                let ylm = spherical_harmonics(l, d);
                for m in 0..=l {
                    qlm[i][m].0 += ylm[m].0;
                    qlm[i][m].1 += ylm[m].1;
                }
            }
            let inv = 1.0 / neighbors[i].len() as f64;
            for m in 0..=l {
                qlm[i][m].0 *= inv;
                qlm[i][m].1 *= inv;
            }
        }

        let mut ql = vec![0.0; n];
        let mut qbl = vec![0.0; n];
        for i in 0..n {
            ql[i] = magnitude(l, &qlm[i]);
            let mut avg = qlm[i].clone();
            for &k in &neighbors[i] {
                for m in 0..=l {
                    avg[m].0 += qlm[k][m].0;
                    avg[m].1 += qlm[k][m].1;
                }
            }
            let inv = 1.0 / (neighbors[i].len() + 1) as f64;
            for v in avg.iter_mut() {
                v.0 *= inv;
                v.1 *= inv;
            }
            qbl[i] = magnitude(l, &avg);
        }
        q.push(ql);
        qbar.push(qbl);
    }

    let mean = |v: &Vec<f64>| if v.is_empty() { 0.0 } else { v.iter().sum::<f64>() / v.len() as f64 };
    SteinhardtResult {
        ls: ls.to_vec(),
        cutoff_nm: cutoff,
        n_waters: n,
        q_mean: q.iter().map(mean).collect(),
        qbar_mean: qbar.iter().map(mean).collect(),
        q,
        qbar,
    }
}

fn magnitude(l: usize, qlm: &[(f64, f64)]) -> f64 {
    let mut s = qlm[0].0 * qlm[0].0 + qlm[0].1 * qlm[0].1;
    for v in &qlm[1..] {
        s += 2.0 * (v.0 * v.0 + v.1 * v.1);
    }
    (4.0 * PI / (2 * l + 1) as f64 * s).sqrt()
}

fn neighbor_list(waters: &[Water], cell: &Cell, cutoff: f64) -> Vec<Vec<usize>> {
    let n = waters.len();
    let mut out = vec![Vec::new(); n];
    let c2 = cutoff * cutoff;
    for i in 0..n {
        for j in (i + 1)..n {
            let d2 = cell.distance_sq(&waters[i].oxygen, &waters[j].oxygen);
            if d2 < c2 && d2 > 0.0 {
                out[i].push(j);
                out[j].push(i);
            }
        }
    }
    out
}

/// Orthonormal spherical harmonics Y_lm(r̂) for m = 0..=l as (re, im).
fn spherical_harmonics(l: usize, d: [f64; 3]) -> Vec<(f64, f64)> {
    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let x = (d[2] / r).clamp(-1.0, 1.0);
    let phi = d[1].atan2(d[0]);
    let p = associated_legendre(l, x);
    (0..=l)
        .map(|m| {
            let k = ((2 * l + 1) as f64 / (4.0 * PI) * factorial_ratio(l - m, l + m)).sqrt();
            let a = k * p[m];
            let mp = m as f64 * phi;
            (a * mp.cos(), a * mp.sin())
        })
        .collect()
}

/// P_l^m(x) for m = 0..=l (Condon–Shortley phase), by the standard recurrences.
fn associated_legendre(l: usize, x: f64) -> Vec<f64> {
    let s = (1.0 - x * x).max(0.0).sqrt();
    let mut out = vec![0.0; l + 1];
    for m in 0..=l {
        // P_m^m
        let mut pmm = 1.0;
        let mut fact = 1.0;
        for _ in 0..m {
            pmm *= -fact * s;
            fact += 2.0;
        }
        if l == m {
            out[m] = pmm;
            continue;
        }
        let mut pm1 = x * (2 * m + 1) as f64 * pmm;
        if l == m + 1 {
            out[m] = pm1;
            continue;
        }
        let mut pll = 0.0;
        for ll in (m + 2)..=l {
            pll = (x * (2 * ll - 1) as f64 * pm1 - (ll + m - 1) as f64 * pmm) / (ll - m) as f64;
            pmm = pm1;
            pm1 = pll;
        }
        out[m] = pll;
    }
    out
}

/// (a!)/(b!) for a <= b.
fn factorial_ratio(a: usize, b: usize) -> f64 {
    let mut r = 1.0;
    for k in (a + 1)..=b {
        r /= k as f64;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_gro_content;

    fn fixture(rel: &str) -> String {
        std::fs::read_to_string(format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), rel))
            .expect("fixture must exist")
    }

    /// Pull the numeric array stored under `key` out of a reference JSON file.
    fn json_array(json: &str, key: &str) -> Vec<f64> {
        let start = json.find(&format!("\"{key}\": [")).expect("key present") + key.len() + 5;
        let end = start + json[start..].find(']').expect("array end");
        json[start..end].split(',').map(|v| v.trim().parse::<f64>().unwrap()).collect()
    }

    fn check_against_freud(gro: &str, reference: &str) -> f64 {
        let data = parse_gro_content(&fixture(gro)).expect("parse");
        let r = calculate_steinhardt(&data.waters, &data.box_dims, DEFAULT_CUTOFF_NM, &[4, 6]);
        let json = fixture(reference);
        // freud computes in float32 (positions f32, harmonics complex64), so agreement is
        // limited to single precision: per-water |Δ| below 1e-4 and mean |Δ| below 1e-5.
        // A wrong neighbour set or normalisation shows up as |Δ| of order 1e-2 or more.
        let mut worst = 0.0f64;
        for (k, l) in [4u32, 6].iter().enumerate() {
            for (field, ours) in [("q", &r.q[k]), ("qbar", &r.qbar[k])] {
                let theirs = json_array(&json, &format!("{field}{l}"));
                assert_eq!(theirs.len(), ours.len(), "{reference} {field}{l} length");
                let mut sum = 0.0;
                for (a, b) in ours.iter().zip(&theirs) {
                    worst = worst.max((a - b).abs());
                    sum += (a - b).abs();
                }
                let mean = sum / ours.len() as f64;
                assert!(mean < 1e-5, "{reference} {field}{l}: mean |Δ| = {mean:e}");
            }
        }
        assert!(worst < 1e-4, "{reference}: max |Δ| = {worst:e}");
        worst
    }

    #[test]
    fn steinhardt_matches_freud_si_crystal() {
        check_against_freud("crystals/444_S1.gro", "reference/steinhardt_444_S1.json");
    }

    #[test]
    fn steinhardt_matches_freud_sii_crystal() {
        check_against_freud("crystals/222_S2.gro", "reference/steinhardt_222_S2.json");
    }

    #[test]
    fn steinhardt_matches_freud_si_with_guests() {
        check_against_freud("guests/sI_CH4_222.gro", "reference/steinhardt_sI_CH4_222.json");
    }

    #[test]
    fn steinhardt_matches_freud_disordered_frame() {
        check_against_freud("reference/grade_prod.gro", "reference/steinhardt_grade_prod.json");
    }

    #[test]
    fn crystal_is_more_ordered_than_disordered_frame() {
        let crystal = parse_gro_content(&fixture("crystals/444_S1.gro")).unwrap();
        let disordered = parse_gro_content(&fixture("reference/grade_prod.gro")).unwrap();
        let a = calculate_steinhardt(&crystal.waters, &crystal.box_dims, DEFAULT_CUTOFF_NM, &[6]);
        let b = calculate_steinhardt(&disordered.waters, &disordered.box_dims, DEFAULT_CUTOFF_NM, &[6]);
        assert!(a.qbar_mean[0] > b.qbar_mean[0] + 0.03, "{} vs {}", a.qbar_mean[0], b.qbar_mean[0]);
    }

    #[test]
    fn rejects_bad_degree_cutoff_and_small_cell() {
        let cell = Cell::orthorhombic([3.0; 3]).unwrap();
        assert!(check_inputs(&cell, 0.35, &[4, 6]).is_ok());
        assert!(check_inputs(&cell, 0.35, &[0]).is_err());
        assert!(check_inputs(&cell, 0.35, &[MAX_L + 1]).is_err());
        assert!(check_inputs(&cell, 0.0, &[6]).is_err());
        assert!(check_inputs(&Cell::orthorhombic([0.6; 3]).unwrap(), 0.35, &[6]).is_err());
    }

    #[test]
    fn isolated_water_is_zero() {
        let w = Water { oxygen: [0.5, 0.5, 0.5], hydrogen1: [0.6, 0.5, 0.5], hydrogen2: [0.5, 0.6, 0.5] };
        let r = calculate_steinhardt(&[w], &[3.0, 3.0, 3.0], DEFAULT_CUTOFF_NM, &[4, 6]);
        assert_eq!(r.q, vec![vec![0.0], vec![0.0]]);
        assert_eq!(r.qbar, vec![vec![0.0], vec![0.0]]);
    }
}
