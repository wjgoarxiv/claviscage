// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

// CHILL+ structure identification (Nguyen & Molinero, J. Phys. Chem. B 2015, 119:9369).
//
// Classifies each water oxygen into 6 classes from the ten-Wolde q3 bond-orientational
// order parameter and the staggered/eclipsed character of its O-O cutoff-neighbor bonds.
// Implemented from the paper + OVITO's documented thresholds.
//
// OVITO integer codes: 0 Other/liquid, 1 Hexagonal ice, 2 Cubic ice, 3 Interfacial ice,
// 4 Gas hydrate (clathrate), 5 Interfacial hydrate.
//
// Reuses the PBC toolbox (geometry::Cell). Units are GROMACS nm; the default O-O cutoff
// is 0.35 nm (3.5 A).

use crate::geometry::{norm3, Cell};
use crate::Water;
use serde::{Deserialize, Serialize};

/// Default O-O candidate-shell cutoff in nm (3.5 A).
pub const DEFAULT_CUTOFF: f64 = 0.35;
/// Bond is "staggered" (ice-like) when the q3 correlation is below this.
const STAGGERED_MAX: f64 = -0.8;
/// Bond is "eclipsed" (clathrate-like) when the q3 correlation is within [ECLIPSED_MIN, ECLIPSED_MAX].
const ECLIPSED_MIN: f64 = -0.35;
const ECLIPSED_MAX: f64 = 0.25;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChillResult {
    /// OVITO class code (0..5) per analyzed water, in input order.
    pub per_atom_class: Vec<u8>,
    /// Count of waters in each class, indexed by class code 0..5.
    pub class_counts: [usize; 6],
    /// Fraction of analyzed waters in each class (sums to 1 when n_water > 0).
    pub fractions: [f64; 6],
    /// Number of waters analyzed (equal to the number of input waters).
    pub n_water: usize,
    /// Unique O-O neighbor edges (`[i, j]` with i < j) from the CHILL+ O-O cutoff search, into
    /// the input water order. Used to draw the classified-water network overlay (the ice/hydrate
    /// lattice) in the 3D scene, colored by class.
    pub edges: Vec<[u32; 2]>,
}

pub const CLASS_LABELS: [&str; 6] = [
    "Other/liquid",
    "Hexagonal ice",
    "Cubic ice",
    "Interfacial ice",
    "Gas hydrate",
    "Interfacial hydrate",
];

// --- minimal complex arithmetic (no num-complex dependency) ---
#[derive(Clone, Copy)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    #[inline]
    fn zero() -> C {
        C { re: 0.0, im: 0.0 }
    }
    #[inline]
    fn add(self, o: C) -> C {
        C {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
    #[inline]
    fn scale(self, s: f64) -> C {
        C {
            re: self.re * s,
            im: self.im * s,
        }
    }
    #[inline]
    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    /// Re(self * conj(other)).
    #[inline]
    fn re_dot_conj(self, o: C) -> f64 {
        self.re * o.re + self.im * o.im
    }
}

/// Complex spherical harmonics Y_3^m evaluated on a unit vector (x, y, z), for m = 0..3.
/// Returned as [Y30, Y31, Y32, Y33]; the m<0 terms follow from Y_3^-m = (-1)^m conj(Y_3^m).
/// Condon-Shortley phase. Prefactors are kept (they set the relative weights between m and
/// cancel in the normalized correlation).
fn y3(x: f64, y: f64, z: f64) -> [C; 4] {
    use std::f64::consts::PI;
    let c0 = 0.25 * (7.0 / PI).sqrt();
    let c1 = -0.125 * (21.0 / PI).sqrt();
    let c2 = 0.25 * (105.0 / (2.0 * PI)).sqrt();
    let c3 = -0.125 * (35.0 / PI).sqrt();

    let z2 = z * z;
    // (x + iy)^2 and (x + iy)^3
    let p2 = C {
        re: x * x - y * y,
        im: 2.0 * x * y,
    };
    let p3 = C {
        re: x * x * x - 3.0 * x * y * y,
        im: 3.0 * x * x * y - y * y * y,
    };

    let y30 = C {
        re: c0 * (5.0 * z2 * z - 3.0 * z),
        im: 0.0,
    };
    let f1 = c1 * (5.0 * z2 - 1.0);
    let y31 = C {
        re: f1 * x,
        im: f1 * y,
    };
    let y32 = p2.scale(c2 * z);
    let y33 = p3.scale(c3);
    [y30, y31, y32, y33]
}

/// Local ten-Wolde q3m for one water, averaged over its neighbor bond directions.
/// Returns [q30, q31, q32, q33]; q30 is real.
fn local_q3(center: [f64; 3], neighbors: &[[f64; 3]], cell: &Cell) -> [C; 4] {
    let mut acc = [C::zero(); 4];
    if neighbors.is_empty() {
        return acc;
    }
    for &nb in neighbors {
        let d = cell.diff(&nb, &center); // r_ij = neighbor - center (MIC)
        let n = norm3(&d);
        if n <= 0.0 {
            continue;
        }
        let (ux, uy, uz) = (d[0] / n, d[1] / n, d[2] / n);
        let ylm = y3(ux, uy, uz);
        for m in 0..4 {
            acc[m] = acc[m].add(ylm[m]);
        }
    }
    let inv = 1.0 / neighbors.len() as f64;
    for m in 0..4 {
        acc[m] = acc[m].scale(inv);
    }
    acc
}

/// |q3|^2 = sum_{m=-3}^{3} |q3m|^2 = q30^2 + 2*(|q31|^2 + |q32|^2 + |q33|^2).
#[inline]
fn q3_norm_sq(q: &[C; 4]) -> f64 {
    q[0].norm_sq() + 2.0 * (q[1].norm_sq() + q[2].norm_sq() + q[3].norm_sq())
}

/// Real q3 correlation c(i,j) in [-1, 1] between two waters' q3m vectors.
fn correlation(qi: &[C; 4], nsi: f64, qj: &[C; 4], nsj: f64) -> f64 {
    if nsi <= 0.0 || nsj <= 0.0 {
        return 0.0;
    }
    // sum_{m=-3}^{3} q3m(i) conj(q3m(j)) is real by the (-1)^m conjugate symmetry:
    // = q30(i)q30(j) + 2*sum_{m=1}^{3} Re(q3m(i) conj(q3m(j)))
    let num = qi[0].re * qj[0].re
        + 2.0 * (qi[1].re_dot_conj(qj[1]) + qi[2].re_dot_conj(qj[2]) + qi[3].re_dot_conj(qj[3]));
    num / (nsi.sqrt() * nsj.sqrt())
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum BondKind {
    Staggered,
    Eclipsed,
    Neither,
}

#[inline]
fn classify_bond(c: f64) -> BondKind {
    if c < STAGGERED_MAX {
        BondKind::Staggered
    } else if c >= ECLIPSED_MIN && c <= ECLIPSED_MAX {
        BondKind::Eclipsed
    } else {
        BondKind::Neither
    }
}

/// CHILL+ molecule class code (0..5) from the staggered/eclipsed bond counts of the 4 bonds.
fn classify_molecule(n_stag: usize, n_ecl: usize) -> u8 {
    if n_stag == 4 {
        2 // cubic ice
    } else if n_stag == 3 && n_ecl == 1 {
        1 // hexagonal ice
    } else if n_ecl == 4 {
        4 // gas hydrate
    } else if n_ecl == 3 {
        5 // interfacial hydrate
    } else if n_stag == 2 || n_stag == 3 {
        3 // interfacial ice
    } else {
        0 // other / liquid
    }
}

/// Sort neighbor candidates by distance, breaking ties by atom index. Shared by both the
/// brute-force and cell-list paths so their `nearest` lists are byte-identical.
fn sort_candidates(cand: &mut [(f64, usize)]) {
    cand.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });
}

/// Orthorhombic box lengths if the cell is (near-)orthorhombic, else None (triclinic).
fn orthorhombic_lengths(cell: &Cell) -> Option<[f64; 3]> {
    let v = cell.vectors();
    let off_diag =
        v[0][1].abs() + v[0][2].abs() + v[1][0].abs() + v[1][2].abs() + v[2][0].abs() + v[2][1].abs();
    if off_diag > 1e-9 {
        return None;
    }
    let lengths = [v[0][0].abs(), v[1][1].abs(), v[2][2].abs()];
    if lengths.iter().any(|&l| l <= 0.0) {
        return None;
    }
    Some(lengths)
}

/// All O neighbors within `cutoff` for each water, sorted by distance then index.
/// Uses an orthorhombic cell list (~O(n)) when the box is orthorhombic, else an O(n^2) scan.
/// Both paths use the same MIC `distance_sq` and the same tie-break, so their output is identical.
fn build_nearest(waters: &[Water], cell: &Cell, cutoff: f64, cutoff_sq: f64) -> Vec<Vec<usize>> {
    match orthorhombic_lengths(cell) {
        Some(lengths) => build_nearest_cell_list(waters, cell, cutoff, cutoff_sq, lengths),
        None => build_nearest_bruteforce(waters, cell, cutoff_sq),
    }
}

/// Reference O(n^2) all-pairs neighbor search. Kept as the correctness oracle for the cell list.
fn build_nearest_bruteforce(waters: &[Water], cell: &Cell, cutoff_sq: f64) -> Vec<Vec<usize>> {
    let n = waters.len();
    let mut nearest: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        let mut cand: Vec<(f64, usize)> = Vec::new();
        for j in 0..n {
            if i == j {
                continue;
            }
            let d2 = cell.distance_sq(&waters[i].oxygen, &waters[j].oxygen);
            if d2 <= cutoff_sq {
                cand.push((d2, j));
            }
        }
        sort_candidates(&mut cand);
        nearest[i] = cand.into_iter().map(|(_, j)| j).collect();
    }
    nearest
}

/// Orthorhombic uniform cell list. Bins oxygens into cells of edge >= cutoff, then scans only
/// each atom's own cell + its face/edge/corner neighbors. Because cell edges are >= cutoff, every
/// neighbor within the cutoff is guaranteed to fall in the searched 27-cell (or fewer, wrapped)
/// block. Distances still go through `cell.distance_sq` (MIC), so results match the brute force.
fn build_nearest_cell_list(
    waters: &[Water],
    cell: &Cell,
    cutoff: f64,
    cutoff_sq: f64,
    lengths: [f64; 3],
) -> Vec<Vec<usize>> {
    let n = waters.len();
    let dims: [usize; 3] = [
        ((lengths[0] / cutoff).floor() as usize).max(1),
        ((lengths[1] / cutoff).floor() as usize).max(1),
        ((lengths[2] / cutoff).floor() as usize).max(1),
    ];
    let flat = |c: [usize; 3]| (c[0] * dims[1] + c[1]) * dims[2] + c[2];
    let cell_of = |pos: &[f64; 3]| -> [usize; 3] {
        let mut idx = [0usize; 3];
        for k in 0..3 {
            let l = lengths[k];
            let mut w = pos[k] - l * (pos[k] / l).floor(); // wrap into [0, l)
            if !(w >= 0.0) {
                w = 0.0;
            }
            if w >= l {
                w = 0.0;
            }
            let c = ((w / l) * dims[k] as f64).floor() as usize;
            idx[k] = c.min(dims[k] - 1);
        }
        idx
    };

    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
    let mut atom_cell: Vec<[usize; 3]> = Vec::with_capacity(n);
    for i in 0..n {
        let c = cell_of(&waters[i].oxygen);
        buckets[flat(c)].push(i);
        atom_cell.push(c);
    }

    // Unique neighbor cell indices along one axis. When a dimension has <=2 cells, wrapping
    // would revisit cells, so scan every cell in that dimension exactly once instead.
    let neigh = |c: usize, ndim: usize| -> Vec<usize> {
        if ndim <= 2 {
            (0..ndim).collect()
        } else {
            vec![(c + ndim - 1) % ndim, c, (c + 1) % ndim]
        }
    };

    let mut nearest: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        let ci = atom_cell[i];
        let nx = neigh(ci[0], dims[0]);
        let ny = neigh(ci[1], dims[1]);
        let nz = neigh(ci[2], dims[2]);
        let mut cand: Vec<(f64, usize)> = Vec::new();
        for &x in &nx {
            for &y in &ny {
                for &z in &nz {
                    for &j in &buckets[flat([x, y, z])] {
                        if i == j {
                            continue;
                        }
                        let d2 = cell.distance_sq(&waters[i].oxygen, &waters[j].oxygen);
                        if d2 <= cutoff_sq {
                            cand.push((d2, j));
                        }
                    }
                }
            }
        }
        sort_candidates(&mut cand);
        nearest[i] = cand.into_iter().map(|(_, j)| j).collect();
    }
    nearest
}

/// CHILL+ classification of every water. `box_dims` accepts the 3/6/9-value GROMACS box.
pub fn calculate_chill(waters: &[Water], box_dims: &[f64], cutoff: f64) -> ChillResult {
    let cell = Cell::from_box_values(box_dims)
        .unwrap_or_else(|_| Cell::orthorhombic([1.0, 1.0, 1.0]).unwrap());
    calculate_chill_cell(waters, &cell, cutoff)
}

pub fn calculate_chill_cell(waters: &[Water], cell: &Cell, cutoff: f64) -> ChillResult {
    let n = waters.len();
    let cutoff_sq = cutoff * cutoff;

    // All O-O cutoff neighbors for each water. Orthorhombic cell list
    // (~O(n)) with an O(n^2) fallback for triclinic boxes; both paths return byte-identical
    // neighbor sets, so classification is unchanged (see the cell_list identity tests).
    let nearest = build_nearest(waters, cell, cutoff, cutoff_sq);

    // Local q3 vectors and their |q3|^2 for each water.
    let mut q: Vec<[C; 4]> = Vec::with_capacity(n);
    let mut qns: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let nb: Vec<[f64; 3]> = nearest[i].iter().map(|&j| waters[j].oxygen).collect();
        let qi = local_q3(waters[i].oxygen, &nb, cell);
        qns.push(q3_norm_sq(&qi));
        q.push(qi);
    }

    let mut class_counts = [0usize; 6];
    let mut per_atom_class = vec![0u8; n];
    for i in 0..n {
        // CHILL+/OVITO classify bulk-like waters with exactly 4 O-O cutoff neighbors.
        // Under- and over-coordinated waters remain Other instead of being forced through a
        // nearest-4 truncation, which overclassifies dense liquid environments as hydrate-like.
        let class = if nearest[i].len() != 4 {
            0
        } else {
            let mut n_stag = 0usize;
            let mut n_ecl = 0usize;
            for &j in &nearest[i] {
                let c = correlation(&q[i], qns[i], &q[j], qns[j]);
                match classify_bond(c) {
                    BondKind::Staggered => n_stag += 1,
                    BondKind::Eclipsed => n_ecl += 1,
                    BondKind::Neither => {}
                }
            }
            classify_molecule(n_stag, n_ecl)
        };
        per_atom_class[i] = class;
        class_counts[class as usize] += 1;
    }

    let mut fractions = [0.0f64; 6];
    if n > 0 {
        for k in 0..6 {
            fractions[k] = class_counts[k] as f64 / n as f64;
        }
    }

    // Unique O-O neighbor edges (i < j) for the network overlay. Deterministic (sort + dedup).
    let mut edges: Vec<[u32; 2]> = Vec::new();
    for i in 0..n {
        for &j in &nearest[i] {
            let a = i.min(j) as u32;
            let b = i.max(j) as u32;
            edges.push([a, b]);
        }
    }
    edges.sort_unstable();
    edges.dedup();

    ChillResult {
        per_atom_class,
        class_counts,
        fractions,
        n_water: n,
        edges,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic PRNG (LCG) so the identity fixtures are reproducible without a rand crate.
    fn lcg_next(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn random_ortho_waters(n: usize, box_l: [f64; 3], seed: u64) -> Vec<Water> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                // Some coordinates fall outside [0, L) on purpose to exercise wrapping.
                let o = [
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[0],
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[1],
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[2],
                ];
                Water { oxygen: o, hydrogen1: o, hydrogen2: o }
            })
            .collect()
    }

    fn water_at(oxygen: [f64; 3]) -> Water {
        Water { oxygen, hydrogen1: oxygen, hydrogen2: oxygen }
    }

    #[test]
    fn cutoff_neighbor_search_keeps_overcoordinated_candidates() {
        let cutoff = DEFAULT_CUTOFF;
        let cutoff_sq = cutoff * cutoff;
        let cell = Cell::orthorhombic([3.0, 3.0, 3.0]).unwrap();
        let waters = vec![
            water_at([1.50, 1.50, 1.50]),
            water_at([1.72, 1.50, 1.50]),
            water_at([1.28, 1.50, 1.50]),
            water_at([1.50, 1.72, 1.50]),
            water_at([1.50, 1.28, 1.50]),
            water_at([1.50, 1.50, 1.72]),
        ];

        let bf = build_nearest_bruteforce(&waters, &cell, cutoff_sq);
        let cl = build_nearest(&waters, &cell, cutoff, cutoff_sq);

        assert_eq!(bf, cl, "cell list must keep the same cutoff-neighbor set as brute force");
        assert_eq!(bf[0].len(), 5, "over-coordinated waters must retain all cutoff neighbors");
    }

    #[test]
    fn triclinic_chill_cutoff_edges_cross_tilted_boundary() {
        let cutoff = DEFAULT_CUTOFF;
        let cutoff_sq = cutoff * cutoff;
        let cell = Cell::restricted_triclinic(2.0, 2.0, 2.0, 0.5, 0.0, 0.0).unwrap();
        let waters = vec![
            water_at([0.10, 0.10, 0.10]),
            // Across the tilted b-boundary: raw delta [0.45, 1.80, 0.0], MIC [-0.05, -0.20, 0.0].
            water_at([0.55, 1.90, 0.10]),
            water_at([0.30, 0.10, 0.10]),
            water_at([0.10, 0.30, 0.10]),
            water_at([0.10, 0.10, 0.30]),
            water_at([1.90, 0.10, 0.10]),
        ];

        let neighbors = build_nearest_bruteforce(&waters, &cell, cutoff_sq);
        assert!(
            neighbors[0].contains(&1),
            "triclinic MIC should connect the tilted-boundary neighbor within the cutoff"
        );
        assert_eq!(
            neighbors[0].len(),
            5,
            "central water should retain all five triclinic cutoff neighbors"
        );

        let result = calculate_chill_cell(&waters, &cell, cutoff);
        assert!(
            result.edges.contains(&[0, 1]),
            "CHILL+ network edges should include the tilted-boundary cutoff pair"
        );
        assert_eq!(
            result.per_atom_class[0],
            0,
            "over-coordinated triclinic waters should remain Other instead of nearest-4 classified"
        );
    }

    #[test]
    fn cell_list_matches_bruteforce_neighbors() {
        // The cell list must reproduce the O(n^2) neighbor sets byte-for-byte across grid sizes,
        // including a box small enough that a dimension has a single cell (dims == 1).
        let cutoff = DEFAULT_CUTOFF;
        let cutoff_sq = cutoff * cutoff;
        for (idx, &box_l) in [[2.0, 2.0, 2.0], [1.5, 2.5, 3.0], [0.6, 0.6, 0.6]].iter().enumerate() {
            let cell = Cell::orthorhombic(box_l).unwrap();
            let waters = random_ortho_waters(500, box_l, 0x00C0FFEE + idx as u64);
            let bf = build_nearest_bruteforce(&waters, &cell, cutoff_sq);
            let cl = build_nearest(&waters, &cell, cutoff, cutoff_sq);
            assert_eq!(bf, cl, "cell list must match brute-force for box {:?}", box_l);
        }
    }

    #[test]
    #[ignore] // manual bench: cargo test --lib chill_cell_list_speedup -- --ignored --nocapture
    fn chill_cell_list_speedup() {
        use std::time::Instant;
        let cutoff = DEFAULT_CUTOFF;
        let cutoff_sq = cutoff * cutoff;
        // ~16k waters in a realistic ~7.9 nm box (~0.032 waters/A^3 liquid density).
        let n = 16_000;
        let box_l = [7.9, 7.9, 7.9];
        let cell = Cell::orthorhombic(box_l).unwrap();
        let waters = random_ortho_waters(n, box_l, 0xBEEF);
        let t0 = Instant::now();
        let bf = build_nearest_bruteforce(&waters, &cell, cutoff_sq);
        let bf_ms = t0.elapsed().as_secs_f64() * 1e3;
        let t1 = Instant::now();
        let cl = build_nearest(&waters, &cell, cutoff, cutoff_sq);
        let cl_ms = t1.elapsed().as_secs_f64() * 1e3;
        assert_eq!(bf, cl);
        eprintln!(
            "CHILL+ neighbor search n={}: bruteforce {:.1} ms, cell-list {:.1} ms, speedup {:.0}x",
            n, bf_ms, cl_ms, bf_ms / cl_ms
        );
    }

    #[test]
    fn cell_list_matches_bruteforce_on_ideal_sii() {
        use crate::parser::parse_gro_content;
        let content = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("222_S2 fixture present");
        let parsed = parse_gro_content(&content).expect("parse 222_S2");
        let cell = Cell::from_box_values(&parsed.box_dims).unwrap();
        let cutoff = DEFAULT_CUTOFF;
        let cutoff_sq = cutoff * cutoff;
        let bf = build_nearest_bruteforce(&parsed.waters, &cell, cutoff_sq);
        let cl = build_nearest(&parsed.waters, &cell, cutoff, cutoff_sq);
        assert_eq!(bf, cl, "cell list must match brute-force on the real sII lattice");
    }

    #[test]
    fn y3_canonical_along_z() {
        // Along +z only Y30 is nonzero: (1/2)sqrt(7/pi) ~ 0.74635.
        let y = y3(0.0, 0.0, 1.0);
        assert!((y[0].re - 0.5 * (7.0 / std::f64::consts::PI).sqrt()).abs() < 1e-9);
        assert!(y[0].im.abs() < 1e-12);
        for m in 1..4 {
            assert!(y[m].norm_sq() < 1e-18, "Y3{} should vanish along z", m);
        }
    }

    #[test]
    fn correlation_self_is_one() {
        // c(i,i) = 1 for any nonzero q3 vector.
        let cell = Cell::orthorhombic([10.0, 10.0, 10.0]).unwrap();
        let center = [5.0, 5.0, 5.0];
        // an arbitrary asymmetric neighbor set so q3 is nonzero
        let nb = vec![
            [5.30, 5.05, 5.00],
            [4.75, 5.28, 5.02],
            [5.02, 4.72, 5.26],
            [4.98, 5.00, 4.71],
        ];
        let qi = local_q3(center, &nb, &cell);
        let ns = q3_norm_sq(&qi);
        assert!(ns > 0.0);
        let c = correlation(&qi, ns, &qi, ns);
        assert!((c - 1.0).abs() < 1e-9, "self-correlation should be 1, got {}", c);
    }

    #[test]
    fn chill_on_ideal_sii_hydrate() {
        use crate::parser::parse_gro_content;
        let content = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("222_S2 fixture present");
        let parsed = parse_gro_content(&content).expect("parse 222_S2");
        let result = calculate_chill(&parsed.waters, &parsed.box_dims, DEFAULT_CUTOFF);
        eprintln!(
            "CHILL+ 222_S2 (n={}) counts {:?} = [other,hex,cubic,ifaceice,hydrate,ifacehyd]",
            result.n_water, result.class_counts
        );
        // Every water of an ideal sII lattice is class 4 (gas hydrate).
        assert_eq!(result.n_water, 1088);
        assert_eq!(result.class_counts, [0, 0, 0, 0, 1088, 0]);
    }

    #[test]
    #[ignore] // manual probe: CHILL_PROBE_GRO=<file.gro> cargo test chill_probe -- --ignored --nocapture
    fn chill_probe() {
        let path = match std::env::var("CHILL_PROBE_GRO") {
            Ok(p) => p,
            Err(_) => return,
        };
        let content = std::fs::read_to_string(&path).expect("probe gro");
        let parsed = crate::parser::parse_gro_content(&content).expect("parse");
        let r = calculate_chill(&parsed.waters, &parsed.box_dims, DEFAULT_CUTOFF);
        eprintln!(
            "CHILL+ probe: n={} counts={:?} other={:.3} ice={:.3} hydrate={:.3}",
            r.n_water,
            r.class_counts,
            r.fractions[0],
            r.fractions[1] + r.fractions[2] + r.fractions[3],
            r.fractions[4] + r.fractions[5],
        );
    }

    #[test]
    fn bond_thresholds() {
        assert_eq!(classify_bond(-0.9), BondKind::Staggered);
        assert_eq!(classify_bond(-0.2), BondKind::Eclipsed);
        assert_eq!(classify_bond(0.6), BondKind::Neither);
        assert_eq!(classify_molecule(4, 0), 2); // cubic
        assert_eq!(classify_molecule(3, 1), 1); // hex
        assert_eq!(classify_molecule(0, 4), 4); // hydrate
        assert_eq!(classify_molecule(0, 3), 5); // interfacial hydrate
        assert_eq!(classify_molecule(2, 0), 3); // interfacial ice
        assert_eq!(classify_molecule(0, 0), 0); // other
    }
}
