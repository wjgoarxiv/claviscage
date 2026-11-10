// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Tetrahedral order of each water from its four nearest oxygens (minimum image, no cutoff):
//! the angular and distance parameters S_g and S_k of Chau & Hardwick (Mol. Phys. 93 (1998) 511,
//! S_g eq. 3) and the orientational order q of Errington & Debenedetti (Nature 409 (2001) 318),
//! which is q = 1 - 4 S_g. The normalisation follows GROMACS `gmx order -Sg -Sk`, the reference
//! these values are checked against (`tests/fixtures/reference/tetra/`).

use serde::Serialize;

use crate::geometry::{dot3, Cell};
use crate::Water;

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct TetraResult {
    pub q_per_water: Vec<f64>,
    pub sg_per_water: Vec<f64>,
    pub sk_per_water: Vec<f64>,
    pub q_mean: f64,
    pub sg_mean: f64,
    pub sk_mean: f64,
}

/// The four nearest oxygens of water `i`, nearest first. Ties keep the lower index.
fn nearest_four(i: usize, oxygens: &[[f64; 3]], cell: &Cell) -> [(f64, usize); 4] {
    let mut nn = [(f64::INFINITY, usize::MAX); 4];
    for (j, o) in oxygens.iter().enumerate() {
        if j == i {
            continue;
        }
        let r2 = cell.distance_sq(&oxygens[i], o);
        if r2 < nn[3].0 {
            let at = nn.iter().position(|&(d, _)| r2 < d).unwrap();
            nn.copy_within(at..3, at + 1);
            nn[at] = (r2, j);
        }
    }
    nn
}

/// The neighbour search compares every pair of waters, so its cost grows with the square of the
/// count; larger frames are refused until it uses a cell list.
pub const MAX_WATERS: usize = 20_000;

/// Needs at least five waters, so that every water has four neighbours, and at most `MAX_WATERS`.
pub fn calculate_tetra_cell(waters: &[Water], cell: &Cell) -> Result<TetraResult, String> {
    if waters.len() < 5 {
        return Err(format!("tetrahedral order needs at least 5 waters, got {}", waters.len()));
    }
    if waters.len() > MAX_WATERS {
        return Err(format!(
            "tetrahedral order is limited to {MAX_WATERS} waters per frame for now, got {}",
            waters.len()
        ));
    }
    let oxygens: Vec<[f64; 3]> = waters.iter().map(|w| w.oxygen).collect();
    let n = oxygens.len();
    let (mut q, mut sg, mut sk) = (Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        let nn = nearest_four(i, &oxygens, cell);
        if !(nn[0].0 > 0.0 && nn[3].0.is_finite()) {
            return Err(format!("water {i} shares its oxygen position with another water or has no four neighbours"));
        }
        let units: Vec<[f64; 3]> = nn
            .iter()
            .map(|&(r2, j)| {
                let d = cell.diff(&oxygens[j], &oxygens[i]);
                let r = r2.sqrt();
                [d[0] / r, d[1] / r, d[2] / r]
            })
            .collect();
        let mut sum = 0.0;
        for j in 0..3 {
            for k in j + 1..4 {
                let c = dot3(&units[j], &units[k]) + 1.0 / 3.0;
                sum += c * c;
            }
        }
        let sg_i = 3.0 * sum / 32.0;
        let r: Vec<f64> = nn.iter().map(|&(r2, _)| r2.sqrt()).collect();
        let mean = r.iter().sum::<f64>() / 4.0;
        let sk_i = r.iter().map(|x| (mean - x) * (mean - x)).sum::<f64>() / (12.0 * mean * mean);
        q.push(1.0 - 4.0 * sg_i);
        sg.push(sg_i);
        sk.push(sk_i);
    }
    let avg = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    Ok(TetraResult { q_mean: avg(&q), sg_mean: avg(&sg), sk_mean: avg(&sk), q_per_water: q, sg_per_water: sg, sk_per_water: sk })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_gro_content;

    fn fixture(rel: &str) -> (Vec<Water>, Cell) {
        let text = std::fs::read_to_string(format!("{}/tests/fixtures/{rel}", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let last = text.lines().rev().find(|l| !l.trim().is_empty()).unwrap();
        let values: Vec<f64> = last.split_whitespace().map(|v| v.parse().unwrap()).collect();
        (parse_gro_content(&text).unwrap().waters, Cell::from_box_values(&values).unwrap())
    }

    /// The frame value from a GROMACS `.xvg` (last data line, second column).
    fn gmx(rel: &str) -> f64 {
        let text = std::fs::read_to_string(format!("{}/tests/fixtures/reference/tetra/{rel}", env!("CARGO_MANIFEST_DIR"))).unwrap();
        let line = text.lines().filter(|l| !l.starts_with(['#', '@'])).last().unwrap();
        line.split_whitespace().nth(1).unwrap().parse().unwrap()
    }

    /// GROMACS runs in single precision and prints six decimals.
    const TOL: f64 = 2e-6;

    #[test]
    fn matches_gmx_order() {
        for (rel, stem) in [
            ("crystals/example_H2O_frame0.gro", "example_H2O_frame0"),
            ("crystals/222_S2.gro", "222_S2"),
            ("benchmark/Ih_332.gro", "Ih_332"),
            ("benchmark/sH_222.gro", "sH_222"),
        ] {
            let (w, cell) = fixture(rel);
            let r = calculate_tetra_cell(&w, &cell).unwrap();
            let (sg, sk) = (gmx(&format!("{stem}_sg.xvg")), gmx(&format!("{stem}_sk.xvg")));
            assert!((r.sg_mean - sg).abs() < TOL, "{rel}: S_g {} vs gmx {sg}", r.sg_mean);
            assert!((r.sk_mean - sk).abs() < TOL, "{rel}: S_k {} vs gmx {sk}", r.sk_mean);
            assert!((r.q_mean - (1.0 - 4.0 * sg)).abs() < 4.0 * TOL, "{rel}: q");
        }
    }

    #[test]
    fn regular_tetrahedron_is_perfect() {
        // A water at the centre of a regular tetrahedron of four waters, in a large box.
        let c = [5.0, 5.0, 5.0];
        let s = 0.276 / 3f64.sqrt();
        let pos = [c, [c[0] + s, c[1] + s, c[2] + s], [c[0] + s, c[1] - s, c[2] - s], [c[0] - s, c[1] + s, c[2] - s], [c[0] - s, c[1] - s, c[2] + s]];
        let waters: Vec<Water> = pos.iter().map(|&o| Water { oxygen: o, hydrogen1: o, hydrogen2: o }).collect();
        let r = calculate_tetra_cell(&waters, &Cell::orthorhombic([10.0; 3]).unwrap()).unwrap();
        assert!(r.sg_per_water[0].abs() < 1e-12 && (r.q_per_water[0] - 1.0).abs() < 1e-12);
        assert!(r.sk_per_water[0].abs() < 1e-12);
    }

    #[test]
    fn too_few_waters_is_an_error() {
        let o = [1.0, 1.0, 1.0];
        let w = vec![Water { oxygen: o, hydrogen1: o, hydrogen2: o }; 4];
        assert!(calculate_tetra_cell(&w, &Cell::orthorhombic([3.0; 3]).unwrap()).is_err());
    }

    fn at(o: [f64; 3]) -> Water {
        Water { oxygen: o, hydrogen1: o, hydrogen2: o }
    }

    #[test]
    fn bad_coordinates_are_errors_not_panics() {
        let cell = Cell::orthorhombic([3.0; 3]).unwrap();
        let mut w: Vec<Water> = (0..6).map(|i| at([0.3 * i as f64, 1.0, 1.0])).collect();
        w[5] = at(w[0].oxygen);
        assert!(calculate_tetra_cell(&w, &cell).is_err(), "coincident oxygens");
        w[5] = at([f64::NAN, 1.0, 1.0]);
        assert!(calculate_tetra_cell(&w, &cell).is_err(), "NaN coordinate");
    }

    #[test]
    fn frames_above_the_limit_are_refused() {
        let w = vec![at([1.0, 1.0, 1.0]); MAX_WATERS + 1];
        let err = calculate_tetra_cell(&w, &Cell::orthorhombic([3.0; 3]).unwrap()).unwrap_err();
        assert!(err.contains("limited to"), "{err}");
    }
}
