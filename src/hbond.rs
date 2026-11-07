// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! H-bond network construction

use crate::geometry::Cell;
use crate::Water;

/// Build hydrogen bond network between water molecules
///
/// H-bond criteria:
/// 1. O-O distance <= rcut (default 0.36 nm)
/// 2. O1-O2-H angle < theta (default 35 degrees)
pub fn build_hbond_network(
    waters: &[Water],
    box_dims: &[f64; 3],
    rcut: f64,
    theta: f64,
) -> Vec<Vec<usize>> {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    build_hbond_network_cell(waters, &cell, rcut, theta)
}

/// Build hydrogen bond network between water molecules using a general periodic cell.
pub fn build_hbond_network_cell(
    waters: &[Water],
    cell: &Cell,
    rcut: f64,
    theta: f64,
) -> Vec<Vec<usize>> {
    match orthorhombic_lengths(cell) {
        Some(lengths) if rcut.is_finite() && rcut > 0.0 => {
            build_hbond_network_cell_list(waters, cell, rcut, theta, lengths)
        }
        _ => build_hbond_network_bruteforce(waters, cell, rcut, theta),
    }
}

/// Reference O(n²) all-pairs H-bond builder. Kept as the correctness oracle for optimized paths.
fn build_hbond_network_bruteforce(
    waters: &[Water],
    cell: &Cell,
    rcut: f64,
    theta: f64,
) -> Vec<Vec<usize>> {
    let n_water = waters.len();
    let mut adj_list: Vec<Vec<usize>> = vec![Vec::new(); n_water];

    if n_water < 2 {
        return adj_list;
    }

    // Simple O(n²) approach - can optimize with spatial grid later
    for i in 0..n_water {
        for j in (i + 1)..n_water {
            let dist = cell.distance(&waters[i].oxygen, &waters[j].oxygen);

            if dist > rcut {
                continue;
            }

            // Check H-bond angles
            let is_hbond = check_hbond_angles_cell(&waters[i], &waters[j], cell, theta);

            if is_hbond {
                adj_list[i].push(j);
                adj_list[j].push(i);
            }
        }
    }

    // Sort adjacency lists
    for list in adj_list.iter_mut() {
        list.sort_unstable();
    }

    adj_list
}

/// Orthorhombic box lengths if the cell has no tilt, else None.
fn orthorhombic_lengths(cell: &Cell) -> Option<[f64; 3]> {
    let v = cell.vectors();
    let off_diag =
        v[0][1].abs() + v[0][2].abs() + v[1][0].abs() + v[1][2].abs() + v[2][0].abs() + v[2][1].abs();
    if off_diag > 1e-9 {
        return None;
    }
    let lengths = [v[0][0].abs(), v[1][1].abs(), v[2][2].abs()];
    if lengths.iter().any(|&l| l <= 0.0 || !l.is_finite()) {
        return None;
    }
    Some(lengths)
}

/// Orthorhombic uniform cell list. It prunes candidate O-O pairs only; exact H-bond
/// acceptance still uses the same MIC distance and angle checks as the brute-force oracle.
fn build_hbond_network_cell_list(
    waters: &[Water],
    cell: &Cell,
    rcut: f64,
    theta: f64,
    lengths: [f64; 3],
) -> Vec<Vec<usize>> {
    let n_water = waters.len();
    let mut adj_list: Vec<Vec<usize>> = vec![Vec::new(); n_water];
    if n_water < 2 {
        return adj_list;
    }

    let mut dims: [usize; 3] = [
        ((lengths[0] / rcut).floor() as usize).max(1),
        ((lengths[1] / rcut).floor() as usize).max(1),
        ((lengths[2] / rcut).floor() as usize).max(1),
    ];
    let max_buckets = (n_water.saturating_mul(64)).max(1);
    while dims[0].saturating_mul(dims[1]).saturating_mul(dims[2]) > max_buckets {
        let axis = if dims[0] >= dims[1] && dims[0] >= dims[2] {
            0
        } else if dims[1] >= dims[2] {
            1
        } else {
            2
        };
        if dims[axis] <= 1 {
            break;
        }
        dims[axis] = dims[axis].div_ceil(2);
    }

    let flat = |c: [usize; 3]| (c[0] * dims[1] + c[1]) * dims[2] + c[2];
    let cell_of = |pos: &[f64; 3]| -> [usize; 3] {
        let mut idx = [0usize; 3];
        for k in 0..3 {
            let l = lengths[k];
            let mut w = pos[k] - l * (pos[k] / l).floor();
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
    let mut atom_cell: Vec<[usize; 3]> = Vec::with_capacity(n_water);
    for i in 0..n_water {
        let c = cell_of(&waters[i].oxygen);
        buckets[flat(c)].push(i);
        atom_cell.push(c);
    }

    let neigh = |c: usize, ndim: usize| -> Vec<usize> {
        if ndim <= 2 {
            (0..ndim).collect()
        } else {
            vec![(c + ndim - 1) % ndim, c, (c + 1) % ndim]
        }
    };

    for i in 0..n_water {
        let ci = atom_cell[i];
        let nx = neigh(ci[0], dims[0]);
        let ny = neigh(ci[1], dims[1]);
        let nz = neigh(ci[2], dims[2]);
        for &x in &nx {
            for &y in &ny {
                for &z in &nz {
                    for &j in &buckets[flat([x, y, z])] {
                        if j <= i {
                            continue;
                        }
                        let dist = cell.distance(&waters[i].oxygen, &waters[j].oxygen);
                        if dist > rcut {
                            continue;
                        }
                        if check_hbond_angles_cell(&waters[i], &waters[j], cell, theta) {
                            adj_list[i].push(j);
                            adj_list[j].push(i);
                        }
                    }
                }
            }
        }
    }

    for list in adj_list.iter_mut() {
        list.sort_unstable();
    }

    adj_list
}

/// Check if any H-bond angle criterion is satisfied
///
/// Uses the standard geometric criterion: the angle O_donor - H - O_acceptor
/// should be greater than (180 - theta) degrees. Equivalently, the angle
/// H - O_donor - O_acceptor should be less than theta degrees.
///
/// We check 4 possible donations (2 hydrogens per water, 2 directions).
#[allow(dead_code)]
fn check_hbond_angles(water1: &Water, water2: &Water, box_dims: &[f64; 3], theta: f64) -> bool {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    check_hbond_angles_cell(water1, water2, &cell, theta)
}

fn check_hbond_angles_cell(water1: &Water, water2: &Water, cell: &Cell, theta: f64) -> bool {
    // Water1 donating H1 to O2: angle H1-O1-O2 < theta
    let angle1 = cell.angle(&water1.oxygen, &water1.hydrogen1, &water2.oxygen);
    if angle1 < theta {
        return true;
    }

    // Water1 donating H2 to O2: angle H2-O1-O2 < theta
    let angle2 = cell.angle(&water1.oxygen, &water1.hydrogen2, &water2.oxygen);
    if angle2 < theta {
        return true;
    }

    // Water2 donating H1 to O1: angle H1-O2-O1 < theta
    let angle3 = cell.angle(&water2.oxygen, &water2.hydrogen1, &water1.oxygen);
    if angle3 < theta {
        return true;
    }

    // Water2 donating H2 to O1: angle H2-O2-O1 < theta
    let angle4 = cell.angle(&water2.oxygen, &water2.hydrogen2, &water1.oxygen);
    if angle4 < theta {
        return true;
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deterministic PRNG (LCG) so identity fixtures are reproducible without a rand crate.
    fn lcg_next(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn random_waters(n: usize, box_l: [f64; 3], seed: u64) -> Vec<Water> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                // Some coordinates fall outside [0, L) on purpose to exercise wrapping.
                let o = [
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[0],
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[1],
                    (lcg_next(&mut s) * 1.4 - 0.2) * box_l[2],
                ];
                Water {
                    oxygen: o,
                    hydrogen1: [o[0] + 0.096, o[1], o[2]],
                    hydrogen2: [o[0], o[1] + 0.096, o[2]],
                }
            })
            .collect()
    }


    #[test]
    fn test_hbond_network() {
        // Two waters close together with H1 of water1 pointing toward O of water2
        let waters = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.096, 0.0, 0.0], // H pointing toward water2
                hydrogen2: [0.0, 0.096, 0.0],
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.184, 0.0, 0.0], // H pointing toward water1
                hydrogen2: [0.28, 0.096, 0.0],
            },
        ];

        let box_dims = [3.0, 3.0, 3.0];
        let adj_list = build_hbond_network(&waters, &box_dims, 0.36, 35.0);

        // They should be bonded (H1-O1-O2 angle is ~0° for water1 donating H1)
        assert!(!adj_list[0].is_empty(), "Expected H-bond between waters");
        assert!(adj_list[0].contains(&1));
        assert!(adj_list[1].contains(&0));
    }

    #[test]
    fn test_hbond_no_bond_bad_angle() {
        // Two waters close together but H pointing away
        let waters = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [-0.096, 0.0, 0.0], // H pointing away from water2
                hydrogen2: [0.0, -0.096, 0.0], // H pointing away
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.376, 0.0, 0.0], // H pointing away from water1
                hydrogen2: [0.28, 0.096, 0.0], // H pointing perpendicular
            },
        ];

        let box_dims = [3.0, 3.0, 3.0];
        let adj_list = build_hbond_network(&waters, &box_dims, 0.36, 35.0);

        // No H-bond because all H atoms point away
        assert!(
            adj_list[0].is_empty(),
            "Expected no H-bond when H atoms point away"
        );
    }

    #[test]
    fn hbond_network_matches_bruteforce_random_orthorhombic() {
        let cell = Cell::orthorhombic([4.0, 4.0, 4.0]).unwrap();
        let waters = random_waters(128, [4.0, 4.0, 4.0], 0x0BAD_F00D);

        let expected = build_hbond_network_bruteforce(&waters, &cell, 0.36, 35.0);
        let actual = build_hbond_network_cell(&waters, &cell, 0.36, 35.0);

        assert_eq!(actual, expected);
    }

    #[test]
    fn hbond_network_matches_bruteforce_on_ideal_sii_fixture() {
        let gro_content = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("222_S2.gro fixture not found");
        let parsed = crate::parser::parse_gro_content(&gro_content).expect("parse 222_S2");
        let cell = Cell::orthorhombic(parsed.box_dims).unwrap();

        let expected = build_hbond_network_bruteforce(&parsed.waters, &cell, 0.36, 35.0);
        let actual = build_hbond_network_cell(&parsed.waters, &cell, 0.36, 35.0);

        assert_eq!(actual, expected);
    }

}
