// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Cage occupancy analysis
//!
//! Two complementary modes:
//! - **Mode A (cage-first):** Given pre-detected cage centers and types, assign guest
//!   molecules to their nearest cage using type-specific radius cutoffs.
//! - **Mode B (guest-centric):** For each guest, build a hydration shell from nearby
//!   waters, compute local order (F3), count primitive rings, and classify the local
//!   environment as liquid / partial / small_cage / large_cage.

use crate::geometry::Cell;
use crate::hbond::build_hbond_network_cell;
use crate::order_params::f3_for_subset_cell;
use crate::{Guest, Water};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashSet, VecDeque};

// =========================================================================
// Mode A — cage-first occupancy
// =========================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestAssignment {
    pub guest_idx: usize,
    pub guest_name: String,
    pub guest_position: [f64; 3],
    pub cage_idx: Option<usize>,
    pub cage_type: Option<String>,
    pub distance_to_center: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CageTypeOccupancy {
    pub cage_type: String,
    pub total: usize,
    pub occupied: usize,
    pub occupancy_ratio: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OccupancyResult {
    pub assignments: Vec<GuestAssignment>,
    pub by_cage_type: Vec<CageTypeOccupancy>,
    pub total_cages: usize,
    pub total_occupied: usize,
    pub overall_occupancy: f64,
    pub theta_s: f64,
    pub theta_l: f64,
}

/// Squared radius cutoff (nm^2) for each cage type.
fn cage_cutoff_sq(cage_type: &str) -> f64 {
    match cage_type {
        "512" => 0.395 * 0.395,
        "51262" => 0.433 * 0.433,
        "51264" => 0.473 * 0.473,
        _ => 0.473 * 0.473,
    }
}

/// Mode A: assign guests to pre-detected cage centers.
///
/// For each guest the nearest cage center (PBC distance) is found.  The guest is
/// assigned only if the distance is below the type-specific cutoff **and** the cage
/// is not already occupied by a closer guest.
pub fn calculate_occupancy(
    cage_centers: &[[f64; 3]],
    cage_types: &[String],
    guests: &[Guest],
    box_dims: &[f64; 3],
) -> OccupancyResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    calculate_occupancy_cell(cage_centers, cage_types, guests, &cell)
}

/// Mode A: assign guests to pre-detected cage centers using a general cell.
pub fn calculate_occupancy_cell(
    cage_centers: &[[f64; 3]],
    cage_types: &[String],
    guests: &[Guest],
    cell: &Cell,
) -> OccupancyResult {
    let n_cages = cage_centers.len();

    // Track which guest occupies each cage and at what distance
    let mut cage_occupant: Vec<Option<usize>> = vec![None; n_cages];
    let mut cage_occupant_dist: Vec<f64> = vec![f64::MAX; n_cages];

    // Each guest proposes to eligible cages from nearest to farthest. A cage keeps
    // the closer proposal; a displaced guest continues with its next candidate.
    // This preserves the documented closer-guest conflict rule without silently
    // dropping a guest that can still occupy another cage.
    let mut candidates: Vec<Vec<(usize, f64)>> = Vec::with_capacity(guests.len());
    for guest in guests {
        let mut eligible = Vec::new();
        for ci in 0..n_cages {
            let dsq = cell.distance_sq(&guest.position, &cage_centers[ci]);
            let cutoff_sq = cage_cutoff_sq(&cage_types[ci]);
            if dsq < cutoff_sq {
                eligible.push((ci, dsq));
            }
        }
        eligible.sort_by(|(left_ci, left_dist), (right_ci, right_dist)| {
            left_dist
                .total_cmp(right_dist)
                .then_with(|| left_ci.cmp(right_ci))
        });
        candidates.push(eligible);
    }

    let mut next_candidate = vec![0usize; guests.len()];
    let mut pending: VecDeque<usize> = (0..guests.len()).collect();
    while let Some(gi) = pending.pop_front() {
        while let Some(&(ci, dsq)) = candidates[gi].get(next_candidate[gi]) {
            next_candidate[gi] += 1;
            let incumbent = cage_occupant[ci];
            let incumbent_dist = cage_occupant_dist[ci];
            let proposal_wins = incumbent.is_none()
                || dsq < incumbent_dist
                || (dsq == incumbent_dist && gi < incumbent.unwrap_or(usize::MAX));
            if proposal_wins {
                cage_occupant[ci] = Some(gi);
                cage_occupant_dist[ci] = dsq;
                if let Some(displaced) = incumbent {
                    pending.push_back(displaced);
                }
                break;
            }
        }
    }

    // Rebuild guest assignments from cage_occupant (source of truth)
    let mut final_cage_for_guest: Vec<Option<usize>> = vec![None; guests.len()];
    for (ci, occupant) in cage_occupant.iter().enumerate() {
        if let Some(gi) = occupant {
            final_cage_for_guest[*gi] = Some(ci);
        }
    }

    // Build assignment records
    let assignments: Vec<GuestAssignment> = guests
        .iter()
        .enumerate()
        .map(|(gi, guest)| {
            let (cage_idx, cage_type, dist) = match final_cage_for_guest[gi] {
                Some(ci) => (
                    Some(ci),
                    Some(cage_types[ci].clone()),
                    cage_occupant_dist[ci].sqrt(),
                ),
                None => (None, None, f64::MAX),
            };
            GuestAssignment {
                guest_idx: gi,
                guest_name: guest.name.clone(),
                guest_position: guest.position,
                cage_idx,
                cage_type,
                distance_to_center: dist,
            }
        })
        .collect();

    // Per-type statistics
    let mut type_set: Vec<String> = cage_types
        .iter()
        .cloned()
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    type_set.sort();

    let by_cage_type: Vec<CageTypeOccupancy> = type_set
        .iter()
        .map(|ct| {
            let total = cage_types.iter().filter(|t| *t == ct).count();
            let occupied = cage_occupant
                .iter()
                .enumerate()
                .filter(|(ci, occ)| occ.is_some() && cage_types[*ci] == *ct)
                .count();
            CageTypeOccupancy {
                cage_type: ct.clone(),
                total,
                occupied,
                occupancy_ratio: if total > 0 {
                    occupied as f64 / total as f64
                } else {
                    0.0
                },
            }
        })
        .collect();

    let total_occupied = cage_occupant.iter().filter(|o| o.is_some()).count();

    // theta_S: occupancy of 512 cages
    let n_512_total = cage_types.iter().filter(|t| t.as_str() == "512").count();
    let n_512_occ = cage_occupant
        .iter()
        .enumerate()
        .filter(|(ci, occ)| occ.is_some() && cage_types[*ci] == "512")
        .count();
    let theta_s = if n_512_total > 0 {
        n_512_occ as f64 / n_512_total as f64
    } else {
        0.0
    };

    // theta_L: occupancy of 51262 + 51264 combined
    let n_large_total = cage_types
        .iter()
        .filter(|t| t.as_str() == "51262" || t.as_str() == "51264")
        .count();
    let n_large_occ = cage_occupant
        .iter()
        .enumerate()
        .filter(|(ci, occ)| {
            occ.is_some() && (cage_types[*ci] == "51262" || cage_types[*ci] == "51264")
        })
        .count();
    let theta_l = if n_large_total > 0 {
        n_large_occ as f64 / n_large_total as f64
    } else {
        0.0
    };

    OccupancyResult {
        assignments,
        by_cage_type,
        total_cages: n_cages,
        total_occupied,
        overall_occupancy: if n_cages > 0 {
            total_occupied as f64 / n_cages as f64
        } else {
            0.0
        },
        theta_s,
        theta_l,
    }
}

// =========================================================================
// Mode B — guest-centric analysis
// =========================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestCentricResult {
    pub guest_idx: usize,
    pub guest_name: String,
    pub n_shell: usize,
    pub shell_f3: f64,
    pub is_hydrate_like: bool,
    pub n_rule_class: String,
    pub ring5_in_shell: usize,
    pub ring6_in_shell: usize,
    pub final_class: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuestCentricAnalysis {
    pub results: Vec<GuestCentricResult>,
    pub n_liquid: usize,
    pub n_partial: usize,
    pub n_small: usize,
    pub n_large: usize,
}

/// Classify environment by coordination number (n-rule).
///
/// Boundaries match the reference classification exactly (strict inequalities):
/// - n < 10 or F3 > 0.04           → LIQUID
/// - 10 < n < 18  (11-17)          → PARTIAL
/// - 18 < n < 21  (19-20)          → SMALL  (5¹² cage: ~20 coordination waters)
/// - n > 27       (28+)            → LARGE  (5¹²6⁴ cage: ~28 coordination waters)
/// - else (n=10, 18, 21-27)        → PARTIAL (conservative fallthrough by design)
pub fn classify_cage_nrule(n_shell: usize, is_hydrate_like: bool) -> &'static str {
    if !is_hydrate_like || n_shell < 10 {
        return "LIQUID";
    }
    if n_shell > 10 && n_shell < 18 {
        return "PARTIAL";
    }
    if n_shell > 18 && n_shell < 21 {
        return "SMALL";
    }
    if n_shell > 27 {
        return "LARGE";
    }
    "PARTIAL"
}

/// Refine the n-rule class using ring counts.
pub fn finalize_class(n_rule: &str, ring5: usize, ring6: usize) -> &'static str {
    if n_rule == "LARGE" && ring5 == 12 && ring6 == 4 {
        return "large_cage";
    }
    if n_rule == "SMALL" && ring5 == 12 {
        return "small_cage";
    }
    if n_rule == "LIQUID" {
        return "liquid_water";
    }
    "partial_cage"
}

/// Find all primitive k-cycles in the H-bond adjacency graph via DFS.
///
/// For each start vertex, DFS with maximum depth k is performed. A cycle is
/// emitted when the path length equals k and the last vertex is adjacent to
/// the start. Cycles are deduplicated using sorted vertex sets.
fn find_k_cycles_prim(adj: &[Vec<usize>], k: usize) -> Vec<BTreeSet<usize>> {
    let n = adj.len();
    let mut seen: HashSet<BTreeSet<usize>> = HashSet::new();
    let mut result: Vec<BTreeSet<usize>> = Vec::new();

    for start in 0..n {
        let mut path = Vec::with_capacity(k);
        path.push(start);
        let mut visited = vec![false; n];
        visited[start] = true;
        cycle_dfs(
            start,
            &mut path,
            &mut visited,
            adj,
            k,
            &mut seen,
            &mut result,
        );
    }

    result
}

fn cycle_dfs(
    start: usize,
    path: &mut Vec<usize>,
    visited: &mut Vec<bool>,
    adj: &[Vec<usize>],
    k: usize,
    seen: &mut HashSet<BTreeSet<usize>>,
    result: &mut Vec<BTreeSet<usize>>,
) {
    let current = *path.last().unwrap();
    let depth = path.len();

    if depth == k {
        // Check if current is adjacent to start (closing the cycle)
        if adj[current].contains(&start) {
            let cycle_set: BTreeSet<usize> = path.iter().copied().collect();
            if cycle_set.len() == k && !seen.contains(&cycle_set) {
                seen.insert(cycle_set.clone());
                result.push(cycle_set);
            }
        }
        return;
    }

    for &next in &adj[current] {
        // Allow returning to start only at depth == k (handled above)
        if next == start && depth < k - 1 {
            continue;
        }
        if visited[next] {
            continue;
        }
        // Pruning: only explore vertices > start to avoid duplicate enumeration
        // when depth == 1 (first step from start)
        if depth == 1 && next < start {
            continue;
        }
        visited[next] = true;
        path.push(next);
        cycle_dfs(start, path, visited, adj, k, seen, result);
        path.pop();
        visited[next] = false;
    }
}

/// Mode B: guest-centric cage classification.
///
/// For each guest molecule:
/// 1. Find the hydration shell (waters within `shell_cutoff` of guest)
/// 2. Compute F3 order parameter for the shell
/// 3. Count 5- and 6-membered rings with all vertices in the shell
/// 4. Apply n-rule then ring-rule to classify the environment
pub fn calculate_guest_centric(
    waters: &[Water],
    guests: &[Guest],
    box_dims: &[f64; 3],
    hbond_rcut: f64,
    hbond_theta_deg: f64,
    shell_cutoff: f64,
) -> GuestCentricAnalysis {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    calculate_guest_centric_cell(
        waters,
        guests,
        &cell,
        hbond_rcut,
        hbond_theta_deg,
        shell_cutoff,
    )
}

/// Mode B: guest-centric cage classification using a general periodic cell.
pub fn calculate_guest_centric_cell(
    waters: &[Water],
    guests: &[Guest],
    cell: &Cell,
    hbond_rcut: f64,
    hbond_theta_deg: f64,
    shell_cutoff: f64,
) -> GuestCentricAnalysis {
    // Build H-bond adjacency list once
    let adj = build_hbond_network_cell(waters, cell, hbond_rcut, hbond_theta_deg);

    // Find primitive 5- and 6-cycles
    let rings5 = find_k_cycles_prim(&adj, 5);
    let rings6 = find_k_cycles_prim(&adj, 6);

    let shell_cutoff_sq = shell_cutoff * shell_cutoff;

    let mut results = Vec::with_capacity(guests.len());
    let mut n_liquid = 0usize;
    let mut n_partial = 0usize;
    let mut n_small = 0usize;
    let mut n_large = 0usize;

    for (gi, guest) in guests.iter().enumerate() {
        // Find shell: water indices within cutoff of guest position
        let shell_indices: Vec<usize> = (0..waters.len())
            .filter(|&i| cell.distance_sq(&waters[i].oxygen, &guest.position) <= shell_cutoff_sq)
            .collect();

        let n_shell = shell_indices.len();
        let shell_set: HashSet<usize> = shell_indices.iter().copied().collect();

        // Compute F3 for the shell
        let shell_f3 = if n_shell >= 2 {
            f3_for_subset_cell(waters, &shell_indices, cell, 0.35)
        } else {
            f64::NAN
        };

        let is_hydrate_like = shell_f3.is_finite() && shell_f3 <= 0.04;

        // Count rings fully contained in shell
        let ring5_in_shell = rings5
            .iter()
            .filter(|ring| ring.iter().all(|v| shell_set.contains(v)))
            .count();
        let ring6_in_shell = rings6
            .iter()
            .filter(|ring| ring.iter().all(|v| shell_set.contains(v)))
            .count();

        let n_rule_class = classify_cage_nrule(n_shell, is_hydrate_like).to_string();
        let final_class = finalize_class(&n_rule_class, ring5_in_shell, ring6_in_shell).to_string();

        match final_class.as_str() {
            "liquid_water" => n_liquid += 1,
            "partial_cage" => n_partial += 1,
            "small_cage" => n_small += 1,
            "large_cage" => n_large += 1,
            _ => {}
        }

        results.push(GuestCentricResult {
            guest_idx: gi,
            guest_name: guest.name.clone(),
            n_shell,
            shell_f3,
            is_hydrate_like,
            n_rule_class,
            ring5_in_shell,
            ring6_in_shell,
            final_class,
        });
    }

    GuestCentricAnalysis {
        results,
        n_liquid,
        n_partial,
        n_small,
        n_large,
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Guest, Water};

    fn dummy_water(pos: [f64; 3]) -> Water {
        Water {
            oxygen: pos,
            hydrogen1: [pos[0] + 0.1, pos[1], pos[2]],
            hydrogen2: [pos[0] - 0.1, pos[1], pos[2]],
        }
    }

    #[test]
    fn test_n_rule_boundaries() {
        // n_shell < 10 -> LIQUID
        assert_eq!(classify_cage_nrule(9, true), "LIQUID");
        // F3 > 0.04 -> LIQUID
        assert_eq!(classify_cage_nrule(20, false), "LIQUID");
        // 11-17 -> PARTIAL
        assert_eq!(classify_cage_nrule(15, true), "PARTIAL");
        // 18-21 -> SMALL
        assert_eq!(classify_cage_nrule(20, true), "SMALL");
        // >27 -> LARGE
        assert_eq!(classify_cage_nrule(28, true), "LARGE");
    }

    #[test]
    fn test_ring_rule() {
        assert_eq!(finalize_class("SMALL", 12, 0), "small_cage");
        assert_eq!(finalize_class("LARGE", 12, 4), "large_cage");
        assert_eq!(finalize_class("LIQUID", 0, 0), "liquid_water");
        assert_eq!(finalize_class("SMALL", 10, 0), "partial_cage"); // wrong ring count
        assert_eq!(finalize_class("PARTIAL", 8, 2), "partial_cage");
    }

    #[test]
    fn test_occupancy_mode_a_simple() {
        // One 512 cage, one guest at its center
        let cage_centers = vec![[0.5, 0.5, 0.5]];
        let cage_types = vec!["512".to_string()];
        let guests = vec![Guest {
            name: "CH4".to_string(),
            position: [0.5, 0.5, 0.5],
        }];
        let box_dims = [3.0, 3.0, 3.0];

        let result = calculate_occupancy(&cage_centers, &cage_types, &guests, &box_dims);

        assert_eq!(result.total_occupied, 1);
        assert_eq!(result.total_cages, 1);
        assert!((result.theta_s - 1.0).abs() < 1e-9);
        assert_eq!(result.assignments[0].cage_idx, Some(0));
        assert_eq!(result.assignments[0].cage_type, Some("512".to_string()));
    }

    #[test]
    fn test_occupancy_mode_a_unassigned() {
        // Guest far from cage (> 0.395 nm)
        let cage_centers = vec![[0.5, 0.5, 0.5]];
        let cage_types = vec!["512".to_string()];
        let guests = vec![Guest {
            name: "CH4".to_string(),
            position: [2.0, 2.0, 2.0],
        }];
        let box_dims = [5.0, 5.0, 5.0];

        let result = calculate_occupancy(&cage_centers, &cage_types, &guests, &box_dims);

        assert_eq!(result.total_occupied, 0);
        assert_eq!(result.assignments[0].cage_idx, None);
    }

    #[test]
    fn test_occupancy_pbc_wrapping() {
        // Guest near box edge, cage near opposite edge — should be close via PBC
        let cage_centers = vec![[0.01, 0.5, 0.5]];
        let cage_types = vec!["51264".to_string()];
        let guests = vec![Guest {
            name: "CO2".to_string(),
            position: [2.99, 0.5, 0.5],
        }];
        let box_dims = [3.0, 3.0, 3.0];

        let result = calculate_occupancy(&cage_centers, &cage_types, &guests, &box_dims);

        // PBC distance = 0.02 nm, well within 0.473 nm cutoff
        assert_eq!(result.total_occupied, 1);
        assert_eq!(result.assignments[0].cage_idx, Some(0));
    }

    #[test]
    fn test_occupancy_conflict_resolution() {
        // Two guests competing for the same cage — closer one wins
        let cage_centers = vec![[0.5, 0.5, 0.5]];
        let cage_types = vec!["51264".to_string()];
        let guests = vec![
            Guest {
                name: "CH4".to_string(),
                position: [0.5, 0.5, 0.5], // distance = 0
            },
            Guest {
                name: "CO2".to_string(),
                position: [0.6, 0.5, 0.5], // distance = 0.1
            },
        ];
        let box_dims = [3.0, 3.0, 3.0];

        let result = calculate_occupancy(&cage_centers, &cage_types, &guests, &box_dims);

        assert_eq!(result.total_occupied, 1);
        // Guest 0 is closer (distance = 0)
        assert_eq!(result.assignments[0].cage_idx, Some(0));
        assert_eq!(result.assignments[1].cage_idx, None);
    }

    #[test]
    fn test_occupancy_displaced_guest_retries_next_nearest_cage() {
        // Both guests prefer cage 0. The closer guest keeps it, while the displaced
        // guest must retry cage 1 instead of being silently dropped.
        let cage_centers = vec![[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let cage_types = vec!["512".to_string(), "512".to_string()];
        let guests = vec![
            Guest {
                name: "CH4".to_string(),
                position: [0.10, 0.0, 0.0],
            },
            Guest {
                name: "CH4".to_string(),
                position: [0.05, 0.0, 0.0],
            },
        ];

        let result = calculate_occupancy(&cage_centers, &cage_types, &guests, &[3.0; 3]);

        assert_eq!(result.total_occupied, 2);
        assert_eq!(result.assignments[0].cage_idx, Some(1));
        assert_eq!(result.assignments[1].cage_idx, Some(0));
    }

    #[test]
    fn test_find_k_cycles_triangle() {
        // A simple triangle (3-cycle) should not appear in k=5 search
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let cycles5 = find_k_cycles_prim(&adj, 5);
        assert!(cycles5.is_empty());
    }

    #[test]
    fn test_find_k_cycles_pentagon() {
        // Pentagon: 0-1-2-3-4-0
        let adj = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        let cycles5 = find_k_cycles_prim(&adj, 5);
        assert_eq!(cycles5.len(), 1);
        let expected: BTreeSet<usize> = (0..5).collect();
        assert_eq!(cycles5[0], expected);
    }
}
