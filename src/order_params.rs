// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! F3 and F4 order parameter calculations for hydrate structure analysis.
//!
//! ## F4 Order Parameter
//!
//! F4 = ⟨cos(3φ)⟩ where φ is the H-O...O-H dihedral angle between hydrogen-bonded
//! water pairs. This parameter distinguishes hydrate crystals from liquid water and ice.
//!
//! ### Expected Values
//!
//! | Structure | F4 Range | Notes |
//! |-----------|----------|-------|
//! | Perfect hydrate crystal (static, 0K) | 0.90 - 0.97 | Ideal dihedral angles, no thermal motion |
//! | Hydrate crystal (MD, 245-275K) | 0.65 - 0.80 | Thermal fluctuations reduce order |
//! | Liquid water | -0.1 to +0.1 | Random dihedral distribution |
//! | Ice Ih | -0.4 to -0.5 | Staggered H-bond geometry |
//!
//! ### Why Static Crystals Have Higher F4
//!
//! In a perfect hydrate crystal (e.g., .gro file at t=0 from crystallography):
//! - H-O...O-H dihedral angles are near their ideal eclipsed values (φ ≈ 0°, ±120°)
//! - cos(3φ) ≈ 1.0 for these ideal angles
//! - Result: F4 ≈ 0.9-0.97
//!
//! In MD simulations at finite temperature (245-275K):
//! - Thermal motion causes dihedral fluctuations of ±15-30° around ideal values
//! - These fluctuations reduce ⟨cos(3φ)⟩
//! - Result: F4 ≈ 0.65-0.80 (as commonly cited in literature)
//!
//! ### Literature References
//!
//! - Moon, Hawtin, Rodger (2007) Faraday Discussions 136:367 - Original F4 definition
//! - Walsh et al. (2011) J. Phys. Chem. C 115:21241 - F4 in nucleation studies
//! - Rodger et al. (1996) Fluid Phase Equilibria 116:326 - "bulk hydrate ≈ 0.7" at finite T
//! - GRADE: Mahmoudinobar & Dias (2019) Comput. Phys. Commun. 244:385
//!
//! ### Implementation
//!
//! This implementation follows the CalF4 reference executable:
//! - For each water pair within O-O cutoff distance
//! - Select the pair with the greatest H-H distance under PBC
//! - Calculate H-O...O-H dihedral angle φ
//! - Compute cos(3φ) and average over all pairs
//!
//! Compact hydrate fixtures are also checked against published F4SPC values.

use crate::geometry::{apply_mic, dot3, norm3, Cell};
use crate::F3F4Result;
use crate::Water;

/// Calculate F3 and F4 order parameters for hydrate structure analysis.
///
/// # F3 (Three-body)
/// Measures tetrahedral ordering of the oxygen network around each water molecule.
/// Lower values indicate more tetrahedral (ice/hydrate-like) local structure.
///
/// # F4 (Four-body/Torsional)
/// Measures H-bond torsional alignment. See module-level documentation for expected
/// values and interpretation.
///
/// # Arguments
/// * `waters` - Slice of Water structs with oxygen and hydrogen positions (nm)
/// * `box_dims` - Periodic box dimensions [x, y, z] in nm
/// * `cutoff` - O-O distance cutoff for neighbor detection (typically 0.35 nm)
///
/// # Returns
/// F3F4Result containing overall and per-atom F3/F4 values.
pub fn calculate_f3_f4(waters: &[Water], box_dims: &[f64; 3], cutoff: f64) -> F3F4Result {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    calculate_f3_f4_cell(waters, &cell, cutoff)
}

/// Calculate F3 and F4 order parameters using a general periodic cell.
pub fn calculate_f3_f4_cell(waters: &[Water], cell: &Cell, cutoff: f64) -> F3F4Result {
    let n_waters = waters.len();

    if n_waters < 2 {
        return F3F4Result {
            f3_overall: 0.0,
            f4_overall: 0.0,
            f3_per_atom: vec![],
            f4_per_atom: vec![],
            n_waters,
        };
    }

    // Build neighbor list (using H-bond network but with distance-only cutoff)
    let neighbors = build_neighbor_list_cell(waters, cell, cutoff);

    // Ideal tetrahedral angle: cos(109.47°) = -1/3
    let cos_tet = -1.0 / 3.0;
    let cos_tet_sq = cos_tet * cos_tet;

    // Calculate F3
    let mut f3_per_atom = vec![0.0; n_waters];
    let mut f3_count = vec![0usize; n_waters];
    let mut f3_total_sum = 0.0;
    let mut f3_total_triplets = 0usize;

    for i in 0..n_waters {
        let neigh_i = &neighbors[i];
        let m = neigh_i.len();

        if m < 2 {
            continue;
        }

        for idx_j in 0..m {
            let j = neigh_i[idx_j];

            // Compute rij with MIC
            let rij = cell.diff(&waters[j].oxygen, &waters[i].oxygen);
            let rij_norm = norm3(&rij);

            if rij_norm < 1e-10 {
                continue;
            }

            for idx_k in (idx_j + 1)..m {
                let k = neigh_i[idx_k];

                // Compute rik with MIC
                let rik = cell.diff(&waters[k].oxygen, &waters[i].oxygen);
                let rik_norm = norm3(&rik);

                if rik_norm < 1e-10 {
                    continue;
                }

                // Calculate cos(theta) between rij and rik
                let cos_theta = dot3(&rij, &rik) / (rij_norm * rik_norm);
                let cos_theta = cos_theta.max(-1.0).min(1.0);

                // F3 term: (cos(theta) * |cos(theta)| + cos²(109.47°))²
                let term = if cos_theta >= 0.0 {
                    cos_theta * cos_theta + cos_tet_sq
                } else {
                    -cos_theta * cos_theta + cos_tet_sq
                };

                let contribution = term * term;
                f3_per_atom[i] += contribution;
                f3_count[i] += 1;
                f3_total_sum += contribution;
                f3_total_triplets += 1;
            }
        }
    }

    // Normalize F3
    for i in 0..n_waters {
        if f3_count[i] > 0 {
            f3_per_atom[i] /= f3_count[i] as f64;
        }
    }

    // Calculate F4 using CalF4's maximum H-H distance selection rule.
    let mut f4_per_atom = vec![0.0; n_waters];
    let mut f4_count = vec![0usize; n_waters];
    let mut f4_total_sum = 0.0;
    let mut f4_total_pairs = 0usize;

    for i in 0..n_waters {
        for &j in &neighbors[i] {
            if j <= i {
                continue;
            }

            let hydrogen_pairs = [
                (
                    &waters[i].hydrogen1,
                    &waters[j].hydrogen1,
                    cell.distance_sq(&waters[i].hydrogen1, &waters[j].hydrogen1),
                ),
                (
                    &waters[i].hydrogen1,
                    &waters[j].hydrogen2,
                    cell.distance_sq(&waters[i].hydrogen1, &waters[j].hydrogen2),
                ),
                (
                    &waters[i].hydrogen2,
                    &waters[j].hydrogen1,
                    cell.distance_sq(&waters[i].hydrogen2, &waters[j].hydrogen1),
                ),
                (
                    &waters[i].hydrogen2,
                    &waters[j].hydrogen2,
                    cell.distance_sq(&waters[i].hydrogen2, &waters[j].hydrogen2),
                ),
            ];
            let (h_i_far, h_j_far, _) = hydrogen_pairs
                .into_iter()
                .max_by(|a, b| a.2.total_cmp(&b.2))
                .expect("a water pair always has four hydrogen pairs");

            // Dihedral: H_i_far - O_i - O_j - H_j_far
            let phi = cell.dihedral(h_i_far, &waters[i].oxygen, &waters[j].oxygen, h_j_far);

            // F4 = cos(3*phi)
            let cos_3phi = (3.0 * phi).cos();

            f4_per_atom[i] += cos_3phi;
            f4_count[i] += 1;
            f4_per_atom[j] += cos_3phi;
            f4_count[j] += 1;

            f4_total_sum += cos_3phi;
            f4_total_pairs += 1;
        }
    }

    // Normalize per-atom F4 for visualization
    for i in 0..n_waters {
        if f4_count[i] > 0 {
            f4_per_atom[i] /= f4_count[i] as f64;
        }
    }

    // F3 = average over all neighbor triplets. Per-atom means are intentionally
    // not averaged here because waters can contribute different triplet counts.
    let f3_overall = if f3_total_triplets > 0 {
        f3_total_sum / f3_total_triplets as f64
    } else {
        0.0
    };

    // F4 = average cos(3φ) over all pairs (bounded in [-1, 1])
    let f4_overall = if f4_total_pairs > 0 {
        f4_total_sum / f4_total_pairs as f64
    } else {
        0.0
    };

    F3F4Result {
        f3_overall,
        f4_overall,
        f3_per_atom,
        f4_per_atom,
        n_waters,
    }
}

/// Build neighbor list based on O-O distance cutoff
#[allow(dead_code)]
fn build_neighbor_list(waters: &[Water], box_dims: &[f64; 3], cutoff: f64) -> Vec<Vec<usize>> {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    build_neighbor_list_cell(waters, &cell, cutoff)
}

fn build_neighbor_list_cell(waters: &[Water], cell: &Cell, cutoff: f64) -> Vec<Vec<usize>> {
    let n_water = waters.len();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n_water];
    let cutoff_sq = cutoff * cutoff;

    for i in 0..n_water {
        for j in (i + 1)..n_water {
            let dist_sq = cell.distance_sq(&waters[i].oxygen, &waters[j].oxygen);

            if dist_sq < cutoff_sq {
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }

    neighbors
}

/// Calculate distance squared with PBC
#[allow(dead_code)]
fn pbc_distance_sq(a: &[f64; 3], b: &[f64; 3], box_dims: &[f64; 3]) -> f64 {
    let delta = apply_mic([b[0] - a[0], b[1] - a[1], b[2] - a[2]], box_dims);
    delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]
}

/// Compute F3 order parameter for a subset of waters by index.
///
/// This is used by guest-centric analysis to evaluate the local ordering
/// of the hydration shell around a guest molecule.
///
/// # Arguments
/// * `waters` - Full water array
/// * `indices` - Subset of water indices to consider
/// * `box_dims` - Periodic box dimensions [x, y, z] in nm
/// * `cutoff_nm` - O-O neighbor cutoff (e.g., 0.35 nm)
pub fn f3_for_subset(
    waters: &[Water],
    indices: &[usize],
    box_dims: &[f64; 3],
    cutoff_nm: f64,
) -> f64 {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    f3_for_subset_cell(waters, indices, &cell, cutoff_nm)
}

/// Compute F3 order parameter for a subset of waters by index using a general cell.
pub fn f3_for_subset_cell(waters: &[Water], indices: &[usize], cell: &Cell, cutoff_nm: f64) -> f64 {
    let cutoff_sq = cutoff_nm * cutoff_nm;
    // Reference angle 104.5° (H-O-H bond angle), matching the reference implementation's convention.
    // NOTE: This differs from the standard tetrahedral angle (109.47°) used in calculate_f3_f4.
    // The hydrate-like threshold of F3 ≤ 0.04 is calibrated to this angle; do not change
    // the angle without recalibrating the threshold.
    let cos_ref_sq = {
        let angle_rad = 104.5f64.to_radians();
        angle_rad.cos().powi(2)
    };

    let mut sum_f3 = 0.0f64;
    let mut n_triplets = 0usize;

    for &i in indices {
        // Find neighbors of i within the subset
        let neighbors: Vec<usize> = indices
            .iter()
            .filter(|&&j| {
                if j == i {
                    return false;
                }
                let dsq = cell.distance_sq(&waters[i].oxygen, &waters[j].oxygen);
                dsq <= cutoff_sq
            })
            .copied()
            .collect();

        if neighbors.len() < 2 {
            continue;
        }

        for a in 0..neighbors.len() {
            for b in (a + 1)..neighbors.len() {
                let j = neighbors[a];
                let k = neighbors[b];
                // Compute angle j-i-k
                let vij = cell.diff(&waters[j].oxygen, &waters[i].oxygen);
                let vik = cell.diff(&waters[k].oxygen, &waters[i].oxygen);
                let dot = vij[0] * vik[0] + vij[1] * vik[1] + vij[2] * vik[2];
                let norm_ij = (vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2]).sqrt();
                let norm_ik = (vik[0] * vik[0] + vik[1] * vik[1] + vik[2] * vik[2]).sqrt();
                if norm_ij < 1e-10 || norm_ik < 1e-10 {
                    continue;
                }
                let cos_phi = (dot / (norm_ij * norm_ik)).clamp(-1.0, 1.0);
                let term = cos_phi * cos_phi.abs() + cos_ref_sq;
                sum_f3 += term * term;
                n_triplets += 1;
            }
        }
    }

    if n_triplets == 0 {
        return f64::NAN;
    }
    sum_f3 / n_triplets as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_f3_f4_basic() {
        let waters = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.1, 0.0, 0.0],
                hydrogen2: [0.0, 0.1, 0.0],
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.38, 0.0, 0.0],
                hydrogen2: [0.28, 0.1, 0.0],
            },
        ];

        let box_dims = [3.0, 3.0, 3.0];
        let result = calculate_f3_f4(&waters, &box_dims, 0.35);

        assert_eq!(result.n_waters, 2);
    }

    #[test]
    fn test_f3_overall_weights_each_triplet() {
        fn water_at(oxygen: [f64; 3]) -> Water {
            Water {
                oxygen,
                hydrogen1: [oxygen[0], oxygen[1] + 0.0957, oxygen[2]],
                hydrogen2: [oxygen[0], oxygen[1], oxygen[2] + 0.0957],
            }
        }

        // The first cluster contributes one 180-degree triplet. The second
        // contributes three 120-degree triplets. Overall F3 is defined over
        // all four triplets, not as an unweighted mean of the two centers.
        let waters = vec![
            water_at([1.0, 1.0, 1.0]),
            water_at([0.7, 1.0, 1.0]),
            water_at([1.3, 1.0, 1.0]),
            water_at([3.0, 3.0, 3.0]),
            water_at([3.3, 3.0, 3.0]),
            water_at([2.85, 3.259_807_621, 3.0]),
            water_at([2.85, 2.740_192_379, 3.0]),
        ];

        let result = calculate_f3_f4(&waters, &[10.0, 10.0, 10.0], 0.35);
        let expected = 1099.0 / 5184.0;

        assert!(
            (result.f3_overall - expected).abs() < 1e-12,
            "F3 overall must weight every neighbor triplet; got {}, expected {}",
            result.f3_overall,
            expected
        );
    }

    #[test]
    fn test_f4_selects_maximum_hydrogen_pair_with_pbc() {
        // Compact pair from prod.gro. The waters straddle the x boundary, and
        // independently choosing each H farthest from the opposite O selects
        // a different torsion than CalF4's maximum H-H pair rule.
        let waters = vec![
            Water {
                oxygen: [0.081, 0.486, 8.232],
                hydrogen1: [0.006, 0.480, 8.292],
                hydrogen2: [0.041, 0.500, 8.146],
            },
            Water {
                oxygen: [3.578, 0.612, 7.976],
                hydrogen1: [3.497, 0.646, 7.938],
                hydrogen2: [3.643, 0.621, 7.907],
            },
        ];

        let result = calculate_f3_f4(&waters, &[3.58839, 3.58839, 12.66492], 0.35);
        let expected = -0.995_154_676_890_555;

        assert!(
            (result.f4_overall - expected).abs() < 1e-12,
            "F4 must use the MIC-farthest H-H pair; got {}, expected {}",
            result.f4_overall,
            expected
        );
    }

    #[test]
    fn test_f4_eclipsed_geometry() {
        // Eclipsed geometry (hydrate-like): dihedral ≈ 0° → F4 ≈ +1
        // Two waters with H-O...O-H dihedral near 0°
        let waters = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.0, 0.0957, 0.0], // H1 pointing +Y (not in bond)
                hydrogen2: [0.0, -0.0478, 0.0829], // H2 pointing toward O2 (in bond)
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.28, 0.0957, 0.0], // H1 pointing +Y (not in bond) - ECLIPSED with H1 of mol1
                hydrogen2: [0.28, -0.0478, -0.0829], // H2 away from O1
            },
        ];

        let box_dims = [3.0, 3.0, 3.0];
        let result = calculate_f3_f4(&waters, &box_dims, 0.35);

        println!("Eclipsed test: F4 = {}", result.f4_overall);
        // Eclipsed: expect cos(3*0) = 1, so F4 should be high (> 0.5)
        assert!(
            result.f4_overall > 0.5,
            "Eclipsed geometry should give F4 > 0.5, got {}",
            result.f4_overall
        );
    }

    #[test]
    fn test_f4_staggered_geometry() {
        let sin_60 = 3.0f64.sqrt() / 2.0;
        let waters = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.0, 0.0957, 0.0],
                hydrogen2: [0.0957, 0.0, 0.0],
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.28, 0.0957 * 0.5, 0.0957 * sin_60],
                hydrogen2: [0.28 - 0.0957, 0.0, 0.0],
            },
        ];

        let result = calculate_f3_f4(&waters, &[3.0, 3.0, 3.0], 0.35);

        assert!(
            (result.f4_overall + 1.0).abs() < 1e-12,
            "A 60-degree staggered torsion must give F4 = -1, got {}",
            result.f4_overall
        );
    }

    #[test]
    fn test_f4_signed_dihedral() {
        // Verify that dihedral calculation gives SIGNED angles
        use crate::geometry::calculate_dihedral;

        let box_dims = [10.0, 10.0, 10.0];

        // Points defining a dihedral with positive angle
        let p1 = [0.0, 1.0, 0.0];
        let p2 = [0.0, 0.0, 0.0];
        let p3 = [1.0, 0.0, 0.0];
        let p4 = [1.0, 1.0, 0.0]; // Same plane as p1-p2-p3, dihedral = 0

        let phi = calculate_dihedral(&p1, &p2, &p3, &p4, &box_dims);
        println!("Planar dihedral: {} rad = {}°", phi, phi * 180.0 / PI);
        assert!(
            (phi.abs() - 0.0).abs() < 0.1,
            "Planar points should give dihedral ≈ 0, got {}",
            phi
        );

        // Rotate p4 out of plane to get non-zero dihedral
        let p4_pos = [1.0, 0.707, 0.707]; // Rotated 45° around X-axis
        let phi_pos = calculate_dihedral(&p1, &p2, &p3, &p4_pos, &box_dims);
        println!(
            "Positive dihedral: {} rad = {}°",
            phi_pos,
            phi_pos * 180.0 / PI
        );

        let p4_neg = [1.0, 0.707, -0.707]; // Rotated -45° around X-axis
        let phi_neg = calculate_dihedral(&p1, &p2, &p3, &p4_neg, &box_dims);
        println!(
            "Negative dihedral: {} rad = {}°",
            phi_neg,
            phi_neg * 180.0 / PI
        );

        // With atan2, phi_pos and phi_neg should have opposite signs
        assert!(
            phi_pos * phi_neg < 0.0,
            "Positive and negative rotations should give opposite sign dihedrals. Got {} and {}",
            phi_pos,
            phi_neg
        );
    }

    /// Validate F4 calculation against F4SPC Python reference implementation
    /// Reference: https://github.com/adibifard/F4SPC
    ///
    /// Test data: First 10 waters from sII hydrate (222_S2.gro)
    /// F4SPC gives F4 = 0.926698 for these 10 waters with cutoff 0.35nm
    #[test]
    fn test_f4_matches_f4spc_reference() {
        // First 10 waters from 222_S2.gro (sII hydrate structure)
        let waters = vec![
            Water {
                oxygen: [0.592, 0.592, 0.592],
                hydrogen1: [0.576, 0.684, 0.569],
                hydrogen2: [0.684, 0.578, 0.570],
            },
            Water {
                oxygen: [0.532, 0.532, 0.857],
                hydrogen1: [0.554, 0.553, 0.766],
                hydrogen2: [0.613, 0.495, 0.893],
            },
            Water {
                oxygen: [0.766, 0.424, 0.965],
                hydrogen1: [0.803, 0.460, 1.045],
                hydrogen2: [0.746, 0.333, 0.987],
            },
            Water {
                oxygen: [0.965, 0.424, 0.766],
                hydrogen1: [0.895, 0.425, 0.832],
                hydrogen2: [0.925, 0.461, 0.688],
            },
            Water {
                oxygen: [0.857, 0.532, 0.532],
                hydrogen1: [0.858, 0.467, 0.462],
                hydrogen2: [0.894, 0.611, 0.493],
            },
            Water {
                oxygen: [0.424, 0.766, 0.965],
                hydrogen1: [0.460, 0.686, 0.928],
                hydrogen2: [0.333, 0.744, 0.985],
            },
            Water {
                oxygen: [0.420, 0.964, 0.765],
                hydrogen1: [0.423, 0.898, 0.835],
                hydrogen2: [0.459, 1.042, 0.804],
            },
            Water {
                oxygen: [0.526, 0.855, 0.530],
                hydrogen1: [0.487, 0.893, 0.609],
                hydrogen2: [0.457, 0.857, 0.465],
            },
            Water {
                oxygen: [0.759, 0.969, 0.425],
                hydrogen1: [0.679, 0.930, 0.460],
                hydrogen2: [0.739, 0.988, 0.334],
            },
            Water {
                oxygen: [0.966, 0.767, 0.426],
                hydrogen1: [0.900, 0.836, 0.428],
                hydrogen2: [1.045, 0.808, 0.462],
            },
        ];

        let box_dims = [3.462, 3.462, 3.462];
        let cutoff = 0.35;

        let result = calculate_f3_f4(&waters, &box_dims, cutoff);

        // F4SPC reference value: 0.926698 (12 pairs found)
        let f4spc_reference = 0.926698;
        let tolerance = 0.01; // 1% tolerance

        println!("F4 from WASM: {:.6}", result.f4_overall);
        println!("F4 from F4SPC: {:.6}", f4spc_reference);
        println!(
            "Difference: {:.6}",
            (result.f4_overall - f4spc_reference).abs()
        );

        assert!(
            (result.f4_overall - f4spc_reference).abs() < tolerance,
            "F4 should match F4SPC reference within {}. Got {}, expected {}",
            tolerance,
            result.f4_overall,
            f4spc_reference
        );

        // Also verify reasonable range for crystalline hydrate
        assert!(
            result.f4_overall > 0.8,
            "sII hydrate should have F4 > 0.8 (crystalline order), got {}",
            result.f4_overall
        );
    }

    /// Validate F4 for sH hydrate structure
    /// Test data: First 15 waters from 222_SH.gro
    /// F4SPC gives F4 = 0.917329 for these 15 waters with cutoff 0.35nm (14 pairs)
    #[test]
    fn test_f4_sh_hydrate() {
        let waters = vec![
            Water {
                oxygen: [0.137, 1.063, 0.65],
                hydrogen1: [0.17, 0.981, 0.686],
                hydrogen2: [0.142, 1.05, 0.555],
            },
            Water {
                oxygen: [0.229, 0.836, 0.782],
                hydrogen1: [0.151, 0.792, 0.815],
                hydrogen2: [0.276, 0.864, 0.861],
            },
            Water {
                oxygen: [0.005, 0.702, 0.873],
                hydrogen1: [-0.074, 0.748, 0.846],
                hydrogen2: [0.002, 0.704, 0.968],
            },
            Water {
                oxygen: [0.226, 1.285, 0.787],
                hydrogen1: [0.198, 1.206, 0.741],
                hydrogen2: [0.274, 1.252, 0.863],
            },
            Water {
                oxygen: [0.001, 1.415, 0.88],
                hydrogen1: [-0.075, 1.367, 0.848],
                hydrogen2: [0.076, 1.367, 0.848],
            },
            Water {
                oxygen: [0.006, 1.677, 0.79],
                hydrogen1: [0.004, 1.586, 0.82],
                hydrogen2: [-0.076, 1.687, 0.741],
            },
            Water {
                oxygen: [0.143, 1.062, 0.371],
                hydrogen1: [0.171, 0.985, 0.321],
                hydrogen2: [0.048, 1.059, 0.369],
            },
            Water {
                oxygen: [0.005, 0.708, 0.134],
                hydrogen1: [-0.075, 0.749, 0.167],
                hydrogen2: [0.001, 0.618, 0.167],
            },
            Water {
                oxygen: [0.005, 0.447, 0.227],
                hydrogen1: [0.001, 0.389, 0.151],
                hydrogen2: [-0.078, 0.432, 0.273],
            },
            Water {
                oxygen: [0.364, 0.915, 0.0],
                hydrogen1: [0.452, 0.878, 0.0],
                hydrogen2: [0.376, 1.009, 0.0],
            },
            Water {
                oxygen: [0.229, 0.836, 0.232],
                hydrogen1: [0.151, 0.792, 0.199],
                hydrogen2: [0.276, 0.864, 0.153],
            },
            Water {
                oxygen: [0.367, 1.193, -0.005],
                hydrogen1: [0.32, 1.222, 0.073],
                hydrogen2: [0.45, 1.242, -0.003],
            },
            Water {
                oxygen: [0.226, 1.279, 0.219],
                hydrogen1: [0.2, 1.207, 0.276],
                hydrogen2: [0.276, 1.338, 0.276],
            },
            Water {
                oxygen: [0.468, 2.11, 0.643],
                hydrogen1: [0.441, 2.188, 0.693],
                hydrogen2: [0.564, 2.113, 0.645],
            },
            Water {
                oxygen: [0.386, 1.893, 0.795],
                hydrogen1: [0.412, 1.966, 0.738],
                hydrogen2: [0.336, 1.834, 0.738],
            },
        ];

        let box_dims = [2.442, 4.23, 2.028];
        let cutoff = 0.35;

        let result = calculate_f3_f4(&waters, &box_dims, cutoff);

        let f4spc_reference = 0.917329;
        let tolerance = 0.01;

        println!("sH F4 from WASM: {:.6}", result.f4_overall);
        println!("sH F4 from F4SPC: {:.6}", f4spc_reference);

        assert!(
            (result.f4_overall - f4spc_reference).abs() < tolerance,
            "sH F4 should match F4SPC. Got {}, expected {}",
            result.f4_overall,
            f4spc_reference
        );
    }
}
