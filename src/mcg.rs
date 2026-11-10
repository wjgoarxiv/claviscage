// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Mutually Coordinated Guest (MCG) order parameter.
//!
//! Reference: B. C. Barnes, G. T. Beckham, D. T. Wu, A. K. Sum, J. Chem. Phys. 140 (2014)
//! 164506, doi:10.1063/1.4871898, Sec. II.
//!
//! 1. For a candidate guest g, every other guest h closer than `r_guest` is considered.
//! 2. The pair (g, h) qualifies when at least `n_water` water oxygens lie within `r_water`
//!    of both guests and inside both cones of half-angle `phi_deg` projected along the
//!    guest–guest vector (from g towards h and from h towards g).
//! 3. g is an MCG monomer when it has at least `n_count` qualifying neighbours.
//! 4. Monomers joined by qualifying pairs form clusters; the order parameter is the number
//!    of guests in the largest cluster.
//!
//! The paper's MCG-3 parameters are the defaults; `n_count = 1` gives MCG-1. Positions are
//! in nm (the paper's 9 Å and 6 Å are 0.9 nm and 0.6 nm). A water is located at its oxygen.

use crate::geometry::Cell;
use crate::{Guest, Water};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct McgParameters {
    pub r_guest_nm: f64,
    pub r_water_nm: f64,
    pub phi_deg: f64,
    pub n_water: usize,
    pub n_count: usize,
}

impl Default for McgParameters {
    /// MCG-3 as published: Rgcut 9 Å, Rwcut 6 Å, φ 45°, Nw 5, Nc 3.
    fn default() -> Self {
        Self { r_guest_nm: 0.9, r_water_nm: 0.6, phi_deg: 45.0, n_water: 5, n_count: 3 }
    }
}

impl McgParameters {
    /// Rejects parameters outside the method's meaning and cells too small for the pair
    /// search to see a single periodic image of each guest (width below 2 × r_guest).
    pub fn check(&self, cell: &Cell) -> Result<(), String> {
        let lengths = [self.r_guest_nm, self.r_water_nm];
        if lengths.iter().any(|r| !(r.is_finite() && *r > 0.0)) {
            return Err("r_guest_nm and r_water_nm must be positive lengths in nm".into());
        }
        if !(self.phi_deg > 0.0 && self.phi_deg <= 180.0) {
            return Err(format!("phi_deg must be within (0, 180], got {}", self.phi_deg));
        }
        if self.n_water == 0 || self.n_count == 0 {
            return Err("n_water and n_count must be at least 1".into());
        }
        let need = 2.0 * self.r_guest_nm.max(self.r_water_nm);
        let width = cell.min_width();
        if width < need {
            return Err(format!("cell is {width:.3} nm across; MCG needs at least {need:.3} nm"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct McgResult {
    pub parameters: McgParameters,
    pub n_guests: usize,
    /// Number of guests that passed step 3.
    pub n_monomers: usize,
    /// The order parameter: guests in the largest MCG cluster.
    pub largest_cluster: usize,
    /// Cluster sizes, largest first.
    pub cluster_sizes: Vec<usize>,
    /// Cluster index per guest (`None` for guests that are not monomers).
    pub guest_cluster: Vec<Option<usize>>,
}

pub fn calculate_mcg(
    waters: &[Water],
    guests: &[Guest],
    box_dims: &[f64; 3],
    params: McgParameters,
) -> McgResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    calculate_mcg_cell(waters, guests, &cell, params)
}

pub fn calculate_mcg_cell(
    waters: &[Water],
    guests: &[Guest],
    cell: &Cell,
    params: McgParameters,
) -> McgResult {
    let ng = guests.len();
    let rg2 = params.r_guest_nm * params.r_guest_nm;
    let rw2 = params.r_water_nm * params.r_water_nm;
    let cos_phi = params.phi_deg.to_radians().cos();

    // Waters within r_water of each guest, as MIC displacement vectors from the guest.
    let shells: Vec<Vec<[f64; 3]>> = guests
        .iter()
        .map(|g| {
            waters
                .iter()
                .map(|w| cell.diff(&w.oxygen, &g.position))
                .filter(|d| dot(d, d) < rw2)
                .collect()
        })
        .collect();

    let mut qualifying: Vec<Vec<usize>> = vec![Vec::new(); ng];
    for a in 0..ng {
        for b in (a + 1)..ng {
            let gab = cell.diff(&guests[b].position, &guests[a].position);
            let dab2 = dot(&gab, &gab);
            if dab2 >= rg2 || dab2 == 0.0 {
                continue;
            }
            let dab = dab2.sqrt();
            let mut count = 0usize;
            for &wa in &shells[a] {
                // Same water seen from b, using the a-frame vector so both share one image.
                let wb = [wa[0] - gab[0], wa[1] - gab[1], wa[2] - gab[2]];
                if dot(&wb, &wb) >= rw2 {
                    continue;
                }
                let na = dot(&wa, &wa).sqrt();
                let nb = dot(&wb, &wb).sqrt();
                if na == 0.0 || nb == 0.0 {
                    continue;
                }
                let cos_a = dot(&wa, &gab) / (na * dab);
                let cos_b = -dot(&wb, &gab) / (nb * dab);
                if cos_a >= cos_phi && cos_b >= cos_phi {
                    count += 1;
                    if count >= params.n_water {
                        break;
                    }
                }
            }
            if count >= params.n_water {
                qualifying[a].push(b);
                qualifying[b].push(a);
            }
        }
    }

    let monomer: Vec<bool> = qualifying.iter().map(|q| q.len() >= params.n_count).collect();
    let mut guest_cluster = vec![None; ng];
    let mut sizes = Vec::new();
    for start in 0..ng {
        if !monomer[start] || guest_cluster[start].is_some() {
            continue;
        }
        let id = sizes.len();
        let mut stack = vec![start];
        guest_cluster[start] = Some(id);
        let mut size = 0;
        while let Some(g) = stack.pop() {
            size += 1;
            for &h in &qualifying[g] {
                if monomer[h] && guest_cluster[h].is_none() {
                    guest_cluster[h] = Some(id);
                    stack.push(h);
                }
            }
        }
        sizes.push(size);
    }

    // Renumber clusters so that index 0 is the largest.
    let mut order: Vec<usize> = (0..sizes.len()).collect();
    order.sort_by(|&x, &y| sizes[y].cmp(&sizes[x]).then(x.cmp(&y)));
    let mut rank = vec![0; sizes.len()];
    for (r, &old) in order.iter().enumerate() {
        rank[old] = r;
    }
    for c in guest_cluster.iter_mut().flatten() {
        *c = rank[*c];
    }
    let cluster_sizes: Vec<usize> = order.iter().map(|&i| sizes[i]).collect();

    McgResult {
        parameters: params,
        n_guests: ng,
        n_monomers: monomer.iter().filter(|&&m| m).count(),
        largest_cluster: cluster_sizes.first().copied().unwrap_or(0),
        cluster_sizes,
        guest_cluster,
    }
}

#[inline]
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_gro_content;

    fn fixture(rel: &str) -> String {
        std::fs::read_to_string(format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), rel))
            .expect("fixture must exist")
    }

    fn water(o: [f64; 3]) -> Water {
        Water { oxygen: o, hydrogen1: o, hydrogen2: o }
    }

    fn guest(p: [f64; 3]) -> Guest {
        Guest { name: "CH4".into(), position: p }
    }

    /// A pentagon of waters in the mid-plane between two guests 0.6 nm apart.
    fn shared_face(center: [f64; 3], radius: f64) -> Vec<Water> {
        (0..5)
            .map(|k| {
                let t = k as f64 * std::f64::consts::TAU / 5.0;
                water([center[0], center[1] + radius * t.cos(), center[2] + radius * t.sin()])
            })
            .collect()
    }

    #[test]
    fn pair_with_shared_face_qualifies_and_wide_ring_does_not() {
        let mcg1 = McgParameters { n_count: 1, ..Default::default() };
        let g = vec![guest([1.2, 1.5, 1.5]), guest([1.8, 1.5, 1.5])];
        // Face radius 0.235 nm at 0.3 nm from each guest: 38° from the axis, inside 45°.
        let near = calculate_mcg(&shared_face([1.5, 1.5, 1.5], 0.235), &g, &[3.0; 3], mcg1);
        assert_eq!((near.n_monomers, near.largest_cluster), (2, 2));
        // Radius 0.4 nm: 53° from the axis, outside both cones.
        let wide = calculate_mcg(&shared_face([1.5, 1.5, 1.5], 0.4), &g, &[3.0; 3], mcg1);
        assert_eq!((wide.n_monomers, wide.largest_cluster), (0, 0));
        // Four waters only: below Nw = 5.
        let four = &shared_face([1.5, 1.5, 1.5], 0.235)[..4];
        let r = calculate_mcg(four, &g, &[3.0; 3], mcg1);
        assert_eq!(r.largest_cluster, 0);
    }

    #[test]
    fn pair_across_periodic_boundary_qualifies() {
        let mcg1 = McgParameters { n_count: 1, ..Default::default() };
        let g = vec![guest([2.85, 1.5, 1.5]), guest([0.45, 1.5, 1.5])];
        let r = calculate_mcg(&shared_face([0.15, 1.5, 1.5], 0.235), &g, &[3.0; 3], mcg1);
        assert_eq!(r.largest_cluster, 2);
    }

    #[test]
    fn mcg3_needs_three_qualifying_neighbours() {
        let g = vec![guest([1.2, 1.5, 1.5]), guest([1.8, 1.5, 1.5])];
        let r = calculate_mcg(&shared_face([1.5, 1.5, 1.5], 0.235), &g, &[3.0; 3], McgParameters::default());
        assert_eq!(r.largest_cluster, 0);
    }

    #[test]
    fn mcg_full_crystal_single_cluster() {
        for (rel, n) in [("guests/sI_CH4_222.gro", 64), ("guests/sII_CH4_222.gro", 192)] {
            let d = parse_gro_content(&fixture(rel)).unwrap();
            let r = calculate_mcg(&d.waters, &d.guests, &d.box_dims, McgParameters::default());
            assert_eq!(r.n_guests, n, "{rel}");
            assert_eq!(r.largest_cluster, n, "{rel}: {:?}", r.cluster_sizes);
            assert_eq!(r.cluster_sizes.len(), 1, "{rel}");
        }
    }

    #[test]
    fn rejects_degenerate_parameters_and_small_cells() {
        let cell = Cell::orthorhombic([3.0; 3]).unwrap();
        assert!(McgParameters::default().check(&cell).is_ok());
        for bad in [
            McgParameters { n_count: 0, ..Default::default() },
            McgParameters { n_water: 0, ..Default::default() },
            McgParameters { r_guest_nm: -1.0, ..Default::default() },
            McgParameters { phi_deg: 0.0, ..Default::default() },
        ] {
            assert!(bad.check(&cell).is_err(), "{bad:?}");
        }
        let small = Cell::orthorhombic([1.2; 3]).unwrap();
        assert!(McgParameters::default().check(&small).is_err());
    }

    #[test]
    fn no_waters_no_clusters() {
        let d = parse_gro_content(&fixture("guests/sI_CH4_222.gro")).unwrap();
        let r = calculate_mcg(&[], &d.guests, &d.box_dims, McgParameters::default());
        assert_eq!((r.n_monomers, r.largest_cluster), (0, 0));
        assert!(r.cluster_sizes.is_empty());
    }

    #[test]
    fn disordered_frame_runs_and_is_below_guest_count() {
        let d = parse_gro_content(&fixture("reference/grade_prod.gro")).unwrap();
        let r = calculate_mcg(&d.waters, &d.guests, &d.box_dims, McgParameters::default());
        assert_eq!(r.n_guests, 790);
        assert!(r.largest_cluster < r.n_guests);
        assert_eq!(r.cluster_sizes.iter().sum::<usize>(), r.n_monomers);
    }
}
