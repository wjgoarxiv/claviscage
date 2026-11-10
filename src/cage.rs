// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from TRACE, Copyright (c) 2025 Jun Wei Hsu and Shiang-Tai Lin, GPL-3.0-or-later; modified.

//! Cage detection using TRACE algorithm
//!
//! Detects all 17 known hydrate cage types from H-bond network.
//!
//! Algorithm:
//! 1. Build H-bond network from water positions
//! 2. Find 4, 5, and 6-membered rings via DFS with geometric validation:
//!    - PBC closure check (prevents ghost rings spanning periodic boundaries)
//!    - Diagonal distance check (prevents collapsed/spurious rings)
//! 3. Build edge-to-ring index for fast adjacency lookups
//! 4. Assemble rings into closed polyhedra using constraint propagation + bounded DFS
//! 5. Validate cages with Euler formula: F - E + V = 2,
//!    every vertex in exactly 3 faces, every edge in exactly 2 faces
//! 6. Classify by face composition (n4, n5, n6) into 17 known types + "other"
//!
//! References:
//! - TRACE: Hsu & Lin, J. Chem. Theory Comput. 2025, 21, 11436-54, doi:10.1021/acs.jctc.5c01459
//! - GRADE diagonal validation: Mahmoudinobar & Dias, Comput. Phys. Commun. (2019)

use crate::geometry::{cross3, dot3, norm3, Cell};
use crate::hbond::build_hbond_network_cell;
use crate::Water;
use crate::{CageResult, IncompleteCageInfo};
use std::collections::{HashMap, HashSet};

/// Maximum faces in any known hydrate cage (5^12·6^8 = 20)
const MAX_CAGE_FACES: usize = 20;

/// Detect hydrate cages from water molecule positions
pub fn detect_cages(waters: &[Water], box_dims: &[f64; 3], rcut: f64, theta: f64) -> CageResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    detect_cages_cell(waters, &cell, rcut, theta)
}

/// Detect hydrate cages from water molecule positions using a general periodic cell.
pub fn detect_cages_cell(waters: &[Water], cell: &Cell, rcut: f64, theta: f64) -> CageResult {
    let n_water = waters.len();
    if n_water < 20 {
        return empty_result();
    }

    // Step 1: Build H-bond network
    let adj_list = build_hbond_network_cell(waters, cell, rcut, theta);

    // Step 2: Reproduce TRACE's native 4-10-member ring walk and pruning.
    let mut rings = find_trace_rings_cell(&adj_list, waters, cell, 10, 90.0, 20.0);
    if rings.is_empty() {
        return empty_result();
    }
    rings.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));

    // Step 3: Build edge-to-ring index
    let edge_to_rings = build_edge_ring_index(&rings);

    // Step 4: Assemble cages (first pass - strict SEC)
    let cage_face_sets = assemble_cages(&rings, &edge_to_rings);

    // Step 5-6: Validate and classify SEC cages
    let mut result = classify_and_build_result_cell(&cage_face_sets, &rings, waters, cell);

    // Step 7: Reproduce TRACE's cup-based T and T+2 passes for nSEC and IC cages.
    // The first pass marks cups consumed by SEC/nSEC cages; the second pass reuses
    // that state so those cups cannot be reported again as incomplete cages.
    let mut cups = find_trace_cups(&rings, &edge_to_rings);
    let complete = find_trace_cage_faces(&mut cups, &rings, 0);
    let incomplete = find_trace_cage_faces(&mut cups, &rings, 2);
    result.incomplete_cages = complete
        .into_iter()
        .filter(|cage| cage.quality == 1)
        .chain(incomplete.into_iter().filter(|cage| cage.quality == 0))
        .map(|cage| build_incomplete_cage_info(&cage.faces, cage.quality, &rings, waters, cell))
        .collect();

    result
}

fn empty_result() -> CageResult {
    CageResult {
        cages_512: vec![],
        cages_51262: vec![],
        cages_51264: vec![],
        cages_all: vec![],
        cage_centers: vec![],
        cage_types: vec![],
        incomplete_cages: vec![],
    }
}

// =========================================================================
// Ring Finding with Geometric Validation
// =========================================================================

/// Find all geometrically valid rings of size 4, 5 and 6
#[allow(dead_code)]
fn find_validated_rings(
    adj_list: &[Vec<usize>],
    waters: &[Water],
    box_dims: &[f64; 3],
    rcut: f64,
) -> Vec<Vec<usize>> {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    find_validated_rings_cell(adj_list, waters, &cell, rcut)
}

fn find_validated_rings_cell(
    adj_list: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
    rcut: f64,
) -> Vec<Vec<usize>> {
    let mut ring_set: HashSet<Vec<usize>> = HashSet::new();
    let mut valid_rings = Vec::new();

    for target_size in [4, 5, 6] {
        let raw_rings = find_rings_topo(adj_list, target_size, &mut ring_set);
        for ring in raw_rings {
            if is_valid_ring_cell(&ring, waters, cell, rcut) {
                valid_rings.push(ring);
            }
        }
    }

    valid_rings
}

fn find_trace_rings_cell(
    adj_list: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
    max_ring_size: usize,
    max_planarity_deg: f64,
    interior_angle_tolerance_deg: f64,
) -> Vec<Vec<usize>> {
    let mut rings = Vec::new();
    let mut visited = vec![false; adj_list.len()];

    for start in 0..adj_list.len() {
        let mut ring = Vec::with_capacity(max_ring_size);
        let mut interior_angles = Vec::with_capacity(max_ring_size.saturating_sub(2));
        let mut seen_vertex_sets = HashSet::new();
        let mut found_on_branch = false;
        trace_ring_dfs(
            start,
            start,
            adj_list,
            waters,
            cell,
            max_ring_size,
            max_planarity_deg,
            interior_angle_tolerance_deg,
            &mut visited,
            &mut ring,
            &mut interior_angles,
            None,
            &mut seen_vertex_sets,
            &mut found_on_branch,
            &mut rings,
        );
    }

    rings
}

#[allow(clippy::too_many_arguments)]
fn trace_ring_dfs(
    current: usize,
    start: usize,
    adj_list: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
    max_ring_size: usize,
    max_planarity_deg: f64,
    interior_angle_tolerance_deg: f64,
    visited: &mut [bool],
    ring: &mut Vec<usize>,
    interior_angles: &mut Vec<f64>,
    reference_normal: Option<[f64; 3]>,
    seen_vertex_sets: &mut HashSet<Vec<usize>>,
    found_on_branch: &mut bool,
    rings: &mut Vec<Vec<usize>>,
) {
    visited[current] = true;
    ring.push(current);
    let level = ring.len() - 1;
    let mut normal = reference_normal;

    if level >= 2 {
        let a = waters[ring[level - 2]].oxygen;
        let b = waters[ring[level - 1]].oxygen;
        let c = waters[ring[level]].oxygen;
        let angle = cell.angle(&b, &a, &c);
        let first = cross3(&cell.diff(&b, &a), &cell.diff(&c, &b));

        if level == 2 {
            normal = Some(first);
            interior_angles.push(angle);
        } else {
            let current_normal = first;
            let base_normal = normal.expect("TRACE ring normal starts at level 2");
            let denominator = norm3(&base_normal) * norm3(&current_normal);
            let mut plane_angle = if denominator <= f64::EPSILON {
                0.0
            } else {
                (dot3(&base_normal, &current_normal) / denominator)
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees()
            };
            if plane_angle > 90.0 {
                plane_angle = 180.0 - plane_angle;
            }
            if plane_angle > max_planarity_deg {
                ring.pop();
                visited[current] = false;
                return;
            }

            if adj_list[current].binary_search(&start).is_ok() {
                if !trace_ring_is_invalid(ring, adj_list, waters, cell) {
                    let mut key = ring.clone();
                    key.sort_unstable();
                    if seen_vertex_sets.insert(key) {
                        rings.push(ring.clone());
                    }
                    *found_on_branch = true;
                }
                ring.pop();
                visited[current] = false;
                return;
            }

            if level > 4 {
                let sides = (level + 2) as f64;
                let ideal_angle = (sides - 2.0) * 180.0 / sides;
                if interior_angles
                    .iter()
                    .any(|value| (value - ideal_angle).abs() > interior_angle_tolerance_deg)
                {
                    ring.pop();
                    visited[current] = false;
                    return;
                }
            }
            interior_angles.push(angle);
        }
    }

    if ring.len() == max_ring_size {
        ring.pop();
        if level >= 2 {
            interior_angles.pop();
        }
        visited[current] = false;
        return;
    }

    for &next in &adj_list[current] {
        if !visited[next] && next > start {
            if !*found_on_branch {
                trace_ring_dfs(
                    next,
                    start,
                    adj_list,
                    waters,
                    cell,
                    max_ring_size,
                    max_planarity_deg,
                    interior_angle_tolerance_deg,
                    visited,
                    ring,
                    interior_angles,
                    normal,
                    seen_vertex_sets,
                    found_on_branch,
                    rings,
                );
            }
            if *found_on_branch && level > 1 {
                if level == 2 {
                    *found_on_branch = false;
                }
                ring.pop();
                interior_angles.pop();
                visited[current] = false;
                return;
            }
        }
    }

    ring.pop();
    if level >= 2 {
        interior_angles.pop();
    }
    visited[current] = false;
}

fn trace_ring_is_invalid(
    ring: &[usize],
    adj_list: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
) -> bool {
    if ring.len() < 7 {
        return false;
    }

    for (i, &vertex) in ring.iter().enumerate() {
        let bonded_in_ring = ring
            .iter()
            .enumerate()
            .filter(|&(j, candidate)| i != j && adj_list[vertex].binary_search(candidate).is_ok())
            .count();
        if bonded_in_ring > 2 {
            return true;
        }
    }

    let mut closure = [0.0; 3];
    for i in 0..ring.len() {
        let from = waters[ring[i]].oxygen;
        let to = waters[ring[(i + 1) % ring.len()]].oxygen;
        let delta = cell.diff(&to, &from);
        closure[0] += delta[0];
        closure[1] += delta[1];
        closure[2] += delta[2];
    }
    closure.iter().any(|value| value.abs() > 0.1)
}

/// Topological ring finding via DFS
fn find_rings_topo(
    adj_list: &[Vec<usize>],
    ring_size: usize,
    ring_set: &mut HashSet<Vec<usize>>,
) -> Vec<Vec<usize>> {
    let n = adj_list.len();
    let mut rings = Vec::new();

    for start in 0..n {
        let mut path = vec![start];
        let mut visited = HashSet::new();
        visited.insert(start);

        ring_dfs(
            start,
            start,
            &mut path,
            &mut visited,
            adj_list,
            ring_size,
            &mut rings,
            ring_set,
        );
    }

    rings
}

fn ring_dfs(
    current: usize,
    start: usize,
    path: &mut Vec<usize>,
    visited: &mut HashSet<usize>,
    adj_list: &[Vec<usize>],
    target_size: usize,
    rings: &mut Vec<Vec<usize>>,
    ring_set: &mut HashSet<Vec<usize>>,
) {
    if path.len() == target_size {
        if adj_list[current].contains(&start) {
            let normalized = normalize_ring(path);
            if !ring_set.contains(&normalized) {
                ring_set.insert(normalized.clone());
                rings.push(normalized);
            }
        }
        return;
    }

    for &next in &adj_list[current] {
        if next == start && path.len() < target_size {
            continue;
        }
        if visited.contains(&next) {
            continue;
        }
        // Pruning: first edge goes to higher index only (avoids duplicates)
        if path.len() == 1 && next < start {
            continue;
        }

        visited.insert(next);
        path.push(next);
        ring_dfs(
            next,
            start,
            path,
            visited,
            adj_list,
            target_size,
            rings,
            ring_set,
        );
        path.pop();
        visited.remove(&next);
    }
}

/// Normalize ring to canonical form (smallest vertex first, smaller neighbor second)
fn normalize_ring(ring: &[usize]) -> Vec<usize> {
    let n = ring.len();
    if n == 0 {
        return vec![];
    }

    let min_pos = ring
        .iter()
        .enumerate()
        .min_by_key(|&(_, v)| *v)
        .map(|(i, _)| i)
        .unwrap();

    let prev = ring[(min_pos + n - 1) % n];
    let next = ring[(min_pos + 1) % n];

    let mut result = Vec::with_capacity(n);
    if prev < next {
        for i in 0..n {
            result.push(ring[(min_pos + n - i) % n]);
        }
    } else {
        for i in 0..n {
            result.push(ring[(min_pos + i) % n]);
        }
    }

    result
}

/// Geometric ring validation (pub(crate) so HTR module can reuse)
///
/// Two checks:
/// 1. PBC closure: sum of MIC displacement vectors around ring ≈ 0
///    (rejects ghost rings that wrap around periodic boundaries)
/// 2. Diagonal distance: non-adjacent vertices must be far enough apart
///    (rejects collapsed/spurious rings where distant atoms appear close via PBC)
#[allow(dead_code)]
pub(crate) fn is_valid_ring(
    ring: &[usize],
    waters: &[Water],
    box_dims: &[f64; 3],
    rcut: f64,
) -> bool {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    is_valid_ring_cell(ring, waters, &cell, rcut)
}

pub(crate) fn is_valid_ring_cell(ring: &[usize], waters: &[Water], cell: &Cell, rcut: f64) -> bool {
    let n = ring.len();
    let positions: Vec<[f64; 3]> = ring.iter().map(|&v| waters[v].oxygen).collect();

    // Check 1: PBC closure
    // For a valid ring, sum of MIC displacements around the ring = 0
    // Non-zero sum means the ring wraps around the periodic box
    let mut pbc_sum = [0.0f64; 3];
    for i in 0..n {
        let j = (i + 1) % n;
        let delta = cell.diff(&positions[j], &positions[i]);
        pbc_sum[0] += delta[0];
        pbc_sum[1] += delta[1];
        pbc_sum[2] += delta[2];
    }
    // Tolerance: 0.1 nm per axis (from TRACE reference)
    if pbc_sum[0].abs() > 0.1 || pbc_sum[1].abs() > 0.1 || pbc_sum[2].abs() > 0.1 {
        return false;
    }

    // Check 2: Diagonal distances (from GRADE v1.14)
    // Pentagon: non-adjacent vertices must be > Lrange5
    // Hexagon: opposite vertices must be > Lrange6
    match n {
        4 => {
            // Square rings: diagonals must be above minimum
            let lrange4 = 1.2 * rcut - 0.08;
            for i in 0..4 {
                let j = (i + 2) % 4; // diagonal vertex
                let d = cell.distance(&positions[i], &positions[j]);
                if d < lrange4 {
                    return false;
                }
            }
        }
        5 => {
            // Lrange5 = 1.6 * rcut - 0.18 (≈0.38 nm for rcut=0.35)
            let lrange5 = 1.6 * rcut - 0.18;
            // Check all non-adjacent pairs in the pentagon
            for i in 0..n {
                for j in (i + 2)..n {
                    // Skip pair (0, n-1) which ARE adjacent in the ring
                    if i == 0 && j == n - 1 {
                        continue;
                    }
                    let d = cell.distance(&positions[i], &positions[j]);
                    if d < lrange5 {
                        return false;
                    }
                }
            }
        }
        6 => {
            // Lrange6 = 2.0 * rcut - 0.26 (≈0.44 nm for rcut=0.35)
            let lrange6 = 2.0 * rcut - 0.26;
            // Check opposite vertex pairs: (0,3), (1,4), (2,5)
            for i in 0..3 {
                let d = cell.distance(&positions[i], &positions[i + 3]);
                if d < lrange6 {
                    return false;
                }
            }
        }
        _ => {}
    }

    true
}

// =========================================================================
// Edge-Ring Index
// =========================================================================

/// Build index: edge → list of ring indices containing that edge
fn build_edge_ring_index(rings: &[Vec<usize>]) -> HashMap<(usize, usize), Vec<usize>> {
    let mut index: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    for (ri, ring) in rings.iter().enumerate() {
        let n = ring.len();
        for i in 0..n {
            let edge = make_edge(ring[i], ring[(i + 1) % n]);
            index.entry(edge).or_default().push(ri);
        }
    }

    index
}

#[derive(Clone)]
struct TraceCup {
    faces: Vec<usize>,
    vertices: Vec<usize>,
    found: bool,
}

fn find_trace_cups(
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
) -> Vec<TraceCup> {
    let mut cups = Vec::new();

    for (base_index, base) in rings.iter().enumerate() {
        let mut lateral_groups = Vec::with_capacity(base.len());
        let mut complete = true;
        for i in 0..base.len() {
            let edge = make_edge(base[i], base[(i + 1) % base.len()]);
            let candidates: Vec<usize> = edge_to_rings
                .get(&edge)
                .into_iter()
                .flatten()
                .copied()
                .filter(|&candidate| candidate != base_index)
                .filter(|&candidate| shared_vertex_count(base, &rings[candidate]) == 2)
                .collect();
            if candidates.is_empty() {
                complete = false;
                break;
            }
            lateral_groups.push(candidates);
        }
        if !complete {
            continue;
        }

        for &first_lateral in &lateral_groups[0] {
            let mut current = Vec::with_capacity(lateral_groups.len());
            current.push(first_lateral);
            let Some(laterals) =
                first_closed_lateral_combination(&lateral_groups, rings, 1, &mut current)
            else {
                continue;
            };
            let mut vertices = HashSet::new();
            for &ring_index in &laterals {
                vertices.extend(rings[ring_index].iter().copied());
            }
            let mut vertices: Vec<usize> = vertices.into_iter().collect();
            vertices.sort_unstable();
            let mut faces = Vec::with_capacity(laterals.len() + 1);
            faces.push(base_index);
            faces.extend(laterals);
            cups.push(TraceCup {
                faces,
                vertices,
                found: false,
            });
        }
    }

    cups.sort_by(|left, right| {
        left.vertices
            .len()
            .cmp(&right.vertices.len())
            .then_with(|| left.vertices.cmp(&right.vertices))
    });
    cups
}

fn first_closed_lateral_combination(
    groups: &[Vec<usize>],
    rings: &[Vec<usize>],
    index: usize,
    current: &mut Vec<usize>,
) -> Option<Vec<usize>> {
    if index == groups.len() {
        return (!current.is_empty()
            && rings_share_edge(&rings[current[0]], &rings[*current.last().unwrap()]))
        .then(|| current.clone());
    }

    for &candidate in &groups[index] {
        if let Some(&previous) = current.last() {
            if !rings_share_edge(&rings[previous], &rings[candidate]) {
                continue;
            }
        }
        current.push(candidate);
        if let Some(combination) =
            first_closed_lateral_combination(groups, rings, index + 1, current)
        {
            current.pop();
            return Some(combination);
        }
        current.pop();
    }
    None
}

fn shared_vertex_count(left: &[usize], right: &[usize]) -> usize {
    left.iter().filter(|vertex| right.contains(vertex)).count()
}

fn rings_share_edge(left: &[usize], right: &[usize]) -> bool {
    for i in 0..left.len() {
        let edge = make_edge(left[i], left[(i + 1) % left.len()]);
        for j in 0..right.len() {
            if edge == make_edge(right[j], right[(j + 1) % right.len()]) {
                return true;
            }
        }
    }
    false
}

#[derive(Clone)]
struct TraceCageFaces {
    faces: Vec<usize>,
    quality: i8,
}

fn find_trace_cage_faces(
    cups: &mut [TraceCup],
    rings: &[Vec<usize>],
    tolerance: usize,
) -> Vec<TraceCageFaces> {
    let mut cups_by_base: HashMap<usize, Vec<usize>> = HashMap::new();
    for (index, cup) in cups.iter().enumerate() {
        cups_by_base.entry(cup.faces[0]).or_default().push(index);
    }

    let mut cages = Vec::new();
    for base_index in 0..cups.len() {
        if cups[base_index].found {
            continue;
        }
        let base = cups[base_index].clone();
        let mut faces: HashSet<usize> = base.faces.iter().copied().collect();
        let mut vertices: HashSet<usize> = base.vertices.iter().copied().collect();
        let remaining = base.faces[1..].to_vec();
        let previous = vec![base.faces[0]; remaining.len()];
        trace_cage_search(
            2,
            &mut faces,
            &mut vertices,
            &remaining,
            &previous,
            &cups_by_base,
            cups,
            rings,
            tolerance,
            &mut cages,
        );
    }
    cages
}

#[allow(clippy::too_many_arguments)]
fn trace_cage_search(
    layer: usize,
    faces: &mut HashSet<usize>,
    vertices: &mut HashSet<usize>,
    remaining: &[usize],
    previous: &[usize],
    cups_by_base: &HashMap<usize, Vec<usize>>,
    cups: &mut [TraceCup],
    rings: &[Vec<usize>],
    tolerance: usize,
    cages: &mut Vec<TraceCageFaces>,
) -> i8 {
    if layer >= 3 {
        let total_face_vertices: usize = faces.iter().map(|&face| rings[face].len()).sum();
        if vertices.len() * 3 <= total_face_vertices + tolerance {
            let face_list: Vec<usize> = faces.iter().copied().collect();
            let quality = classify_cage_quality(&face_list, rings);
            if quality > 0 {
                cages.push(TraceCageFaces {
                    faces: face_list,
                    quality,
                });
                return 1;
            }
            if quality == 0 {
                if tolerance > 0 {
                    cages.push(TraceCageFaces {
                        faces: face_list,
                        quality,
                    });
                }
                return 2;
            }
            return 2;
        }
    }
    if layer >= 4 || remaining.is_empty() || remaining.len() > 12 {
        return 0;
    }

    let mut lateral_groups: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut seen_bases = HashSet::new();
    for (&face, &previous_face) in remaining.iter().zip(previous) {
        if !seen_bases.insert(face) {
            continue;
        }
        let candidates: Vec<usize> = cups_by_base
            .get(&face)
            .into_iter()
            .flatten()
            .copied()
            .filter(|&cup_index| cups[cup_index].faces.contains(&previous_face))
            .collect();
        if !candidates.is_empty() {
            lateral_groups.push((face, candidates));
        }
    }
    if layer == 2 && lateral_groups.len() < 4 {
        return 0;
    }

    let mut total_product = 1usize;
    for (_, candidates) in &lateral_groups {
        total_product = total_product.saturating_mul(candidates.len());
        if total_product > 10 {
            return 0;
        }
    }

    let mut combinations = Vec::new();
    trace_cup_combinations(&lateral_groups, cups, 0, &mut Vec::new(), &mut combinations);
    for combination in combinations {
        let saved_faces = faces.clone();
        let saved_vertices = vertices.clone();
        let mut next_remaining = Vec::new();
        let mut next_previous = Vec::new();

        for (group_index, &cup_index) in combination.iter().enumerate() {
            let key = lateral_groups[group_index].0;
            let cup = &cups[cup_index];
            for &face in &cup.faces {
                if faces.insert(face) {
                    vertices.extend(cup.vertices.iter().copied());
                    next_remaining.push(face);
                    next_previous.push(key);
                }
            }
        }

        let status = trace_cage_search(
            layer + 1,
            faces,
            vertices,
            &next_remaining,
            &next_previous,
            cups_by_base,
            cups,
            rings,
            tolerance,
            cages,
        );
        if status != 0 {
            if status == 1 || tolerance > 0 {
                for cup_index in combination {
                    cups[cup_index].found = true;
                }
            }
            return status;
        }
        *faces = saved_faces;
        *vertices = saved_vertices;
    }
    0
}

fn trace_cup_combinations(
    groups: &[(usize, Vec<usize>)],
    cups: &[TraceCup],
    index: usize,
    current: &mut Vec<usize>,
    combinations: &mut Vec<Vec<usize>>,
) {
    if index == groups.len() {
        combinations.push(current.clone());
        return;
    }
    for &cup_index in &groups[index].1 {
        if cups[cup_index].found {
            continue;
        }
        current.push(cup_index);
        trace_cup_combinations(groups, cups, index + 1, current, combinations);
        current.pop();
    }
}

#[inline]
fn make_edge(u: usize, v: usize) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

// =========================================================================
// Cage Assembly via Deterministic Vertex-Seeded Growth
// =========================================================================

/// Assemble rings into cages using constraint propagation + bounded DFS.
///
/// Algorithm:
/// 1. For each ring, try to grow a cage containing it
/// 2. Growth uses constraint propagation (forced moves) + DFS (branching)
/// 3. Validate closed polyhedra with Euler formula
/// 4. Deduplicate by sorted face-index key
/// 5. Filter to keep only physically valid cage types (512, 51262, 51264)
///
/// Note: In a hydrate crystal, rings (faces) CAN be shared between adjacent
/// cages. So we do NOT enforce ring partitioning.
fn assemble_cages(
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
) -> Vec<Vec<usize>> {
    let mut found_cages: Vec<Vec<usize>> = Vec::new();
    let mut found_cage_keys: HashSet<Vec<usize>> = HashSet::new();

    for seed_ri in 0..rings.len() {
        // Every closure this ring participates in, not just the first one found. A face is
        // shared by two cages, so a seed that stops at its first closure reports one of them
        // and hides the other -- which is how 4 of the 64 large cages in an sII crystal went
        // missing: every hexagon of those cages closed into the neighbour first.
        for cage_faces in grow_cage_deterministic(&[seed_ri], rings, edge_to_rings) {
            let mut key = cage_faces.clone();
            key.sort();
            if found_cage_keys.insert(key) {
                found_cages.push(cage_faces);
            }
        }
    }

    found_cages
}

/// Grow a cage from seed rings using constraint propagation + bounded DFS.
///
/// Phase 1 (Propagation): Iteratively close open edges that have exactly 1
/// valid candidate ring. This is deterministic and fast.
///
/// Phase 2 (Branching): When all remaining open edges have 2+ candidates,
/// pick the most-constrained edge (fewest candidates) and try each candidate
/// via recursive DFS with backtracking. Bounded to prevent explosion.
///
/// Returns EVERY valid cage the seed can close into.
fn grow_cage_deterministic(
    seed_rings: &[usize],
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
) -> Vec<Vec<usize>> {
    let mut faces: Vec<usize> = seed_rings.to_vec();
    let mut face_set: HashSet<usize> = faces.iter().cloned().collect();
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();

    // Initialize counts from seed rings
    for &ri in &faces {
        add_ring_to_counts(ri, rings, &mut edge_count, &mut vertex_count);
    }

    // Check initial validity
    if edge_count.values().any(|&c| c > 2) || vertex_count.values().any(|&c| c > 3) {
        return vec![];
    }

    let mut out = Vec::new();
    grow_cage_recursive(
        &mut faces,
        &mut face_set,
        &mut edge_count,
        &mut vertex_count,
        rings,
        edge_to_rings,
        &mut out,
    );
    out
}

/// Add a ring's edges and vertices to the count maps
#[inline]
fn add_ring_to_counts(
    ri: usize,
    rings: &[Vec<usize>],
    edge_count: &mut HashMap<(usize, usize), usize>,
    vertex_count: &mut HashMap<usize, usize>,
) {
    let ring = &rings[ri];
    let n = ring.len();
    for i in 0..n {
        let edge = make_edge(ring[i], ring[(i + 1) % n]);
        *edge_count.entry(edge).or_insert(0) += 1;
    }
    for &v in ring {
        *vertex_count.entry(v).or_insert(0) += 1;
    }
}

/// Remove a ring's edges and vertices from the count maps
#[inline]
fn remove_ring_from_counts(
    ri: usize,
    rings: &[Vec<usize>],
    edge_count: &mut HashMap<(usize, usize), usize>,
    vertex_count: &mut HashMap<usize, usize>,
) {
    let ring = &rings[ri];
    let n = ring.len();
    for i in 0..n {
        let edge = make_edge(ring[i], ring[(i + 1) % n]);
        let c = edge_count.get_mut(&edge).unwrap();
        *c -= 1;
        if *c == 0 {
            edge_count.remove(&edge);
        }
    }
    for &v in ring {
        let c = vertex_count.get_mut(&v).unwrap();
        *c -= 1;
        if *c == 0 {
            vertex_count.remove(&v);
        }
    }
}

/// Check if adding ring ri would violate SEC constraints
#[inline]
fn ring_addition_valid(
    ri: usize,
    rings: &[Vec<usize>],
    edge_count: &HashMap<(usize, usize), usize>,
    vertex_count: &HashMap<usize, usize>,
) -> bool {
    let ring = &rings[ri];
    let n = ring.len();
    // Edge constraint: no edge would exceed 2
    let edge_ok = (0..n).all(|i| {
        let e = make_edge(ring[i], ring[(i + 1) % n]);
        edge_count.get(&e).copied().unwrap_or(0) < 2
    });
    if !edge_ok {
        return false;
    }
    // Vertex constraint: no vertex would exceed 3
    ring.iter()
        .all(|v| vertex_count.get(v).copied().unwrap_or(0) < 3)
}

/// Get valid candidates for an open edge
fn get_candidates(
    edge: (usize, usize),
    face_set: &HashSet<usize>,
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    edge_count: &HashMap<(usize, usize), usize>,
    vertex_count: &HashMap<usize, usize>,
) -> Vec<usize> {
    edge_to_rings
        .get(&edge)
        .map(|r| {
            r.iter()
                .filter(|ri| !face_set.contains(ri))
                .filter(|&&ri| ring_addition_valid(ri, rings, edge_count, vertex_count))
                .cloned()
                .collect()
        })
        .unwrap_or_default()
}

/// Recursive cage growth with constraint propagation and bounded DFS.
fn grow_cage_recursive(
    faces: &mut Vec<usize>,
    face_set: &mut HashSet<usize>,
    edge_count: &mut HashMap<(usize, usize), usize>,
    vertex_count: &mut HashMap<usize, usize>,
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    out: &mut Vec<Vec<usize>>,
) {
    // Every exit unwinds `propagated`, so the caller's state is exactly as it was. That matters
    // now that the search does not stop at its first cage: a branch that closes must be undone
    // so its siblings are explored from the same state.
    let mut propagated: Vec<usize> = Vec::new();

    macro_rules! unwind {
        () => {{
            for &ri in propagated.iter().rev() {
                faces.pop();
                face_set.remove(&ri);
                remove_ring_from_counts(ri, rings, edge_count, vertex_count);
            }
            return;
        }};
    }

    loop {
        // Sorted, not hash order: the walk below closes the FIRST open edge that admits a
        // single candidate, so iteration order decides which cage gets grown -- and Rust
        // reseeds every HashMap, which made the cage count differ between two runs on one
        // unchanged structure.
        let mut open_edges: Vec<(usize, usize)> = edge_count
            .iter()
            .filter(|(_, &c)| c == 1)
            .map(|(&e, _)| e)
            .collect();
        open_edges.sort_unstable();

        if open_edges.is_empty() {
            if validate_euler(faces, rings) {
                out.push(faces.clone());
            }
            unwind!();
        }

        if faces.len() >= MAX_CAGE_FACES {
            unwind!();
        }

        // Phase 1: constraint propagation — close any edge with exactly one candidate
        let mut forced_any = false;
        for &edge in &open_edges {
            if edge_count.get(&edge).copied().unwrap_or(0) != 1 {
                continue;
            }
            let candidates = get_candidates(
                edge,
                face_set,
                rings,
                edge_to_rings,
                edge_count,
                vertex_count,
            );
            if candidates.is_empty() {
                unwind!(); // dead end: an open edge no ring can close
            }
            if candidates.len() == 1 {
                let ri = candidates[0];
                faces.push(ri);
                face_set.insert(ri);
                add_ring_to_counts(ri, rings, edge_count, vertex_count);
                propagated.push(ri);
                forced_any = true;
                break; // counts changed — rescan
            }
        }
        if forced_any {
            continue;
        }

        // Phase 2: nothing forced. Branch on the most constrained edge, and take EVERY
        // candidate — the first one that closes is not the only cage here.
        let mut best_candidates: Vec<usize> = Vec::new();
        let mut best_count = usize::MAX;

        for &edge in &open_edges {
            if edge_count.get(&edge).copied().unwrap_or(0) != 1 {
                continue;
            }
            let candidates = get_candidates(
                edge,
                face_set,
                rings,
                edge_to_rings,
                edge_count,
                vertex_count,
            );
            if candidates.is_empty() {
                unwind!();
            }
            if candidates.len() < best_count {
                best_count = candidates.len();
                best_candidates = candidates;
            }
        }

        if best_candidates.is_empty() {
            unwind!();
        }

        for candidate in best_candidates {
            faces.push(candidate);
            face_set.insert(candidate);
            add_ring_to_counts(candidate, rings, edge_count, vertex_count);

            grow_cage_recursive(
                faces,
                face_set,
                edge_count,
                vertex_count,
                rings,
                edge_to_rings,
                out,
            );

            faces.pop();
            face_set.remove(&candidate);
            remove_ring_from_counts(candidate, rings, edge_count, vertex_count);
        }

        unwind!();
    }
}

// =========================================================================
// Euler Validation
// =========================================================================

/// Validate a cage using Euler's formula for convex polyhedra.
///
/// Matching TRACE's SEC_cage_identify():
/// 1. F - E + V = 2
/// 2. Every edge shared by exactly 2 faces
/// 3. Every vertex in exactly 3 faces
///
/// Returns true only for "perfect SEC cages" (type 2 in TRACE).
pub(crate) fn validate_euler(face_indices: &[usize], rings: &[Vec<usize>]) -> bool {
    let f = face_indices.len();
    if f < 4 {
        return false;
    }

    let mut edge_face_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_face_count: HashMap<usize, usize> = HashMap::new();

    for &fi in face_indices {
        let ring = &rings[fi];
        let n = ring.len();
        for i in 0..n {
            *vertex_face_count.entry(ring[i]).or_insert(0) += 1;
            let edge = make_edge(ring[i], ring[(i + 1) % n]);
            *edge_face_count.entry(edge).or_insert(0) += 1;
        }
    }

    let e = edge_face_count.len();
    let v = vertex_face_count.len();

    // Euler's formula
    if (f as i64) - (e as i64) + (v as i64) != 2 {
        return false;
    }

    // Every edge in exactly 2 faces
    if edge_face_count.values().any(|&c| c != 2) {
        return false;
    }

    // Every vertex in exactly 3 faces (SEC criterion)
    if vertex_face_count.values().any(|&c| c != 3) {
        return false;
    }

    true
}

// =========================================================================
// Classification
// =========================================================================

/// Classify cages by face composition and build the final result
#[allow(dead_code)]
fn classify_and_build_result(
    cage_face_sets: &[Vec<usize>],
    rings: &[Vec<usize>],
    waters: &[Water],
    box_dims: &[f64; 3],
) -> CageResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    classify_and_build_result_cell(cage_face_sets, rings, waters, &cell)
}

fn classify_and_build_result_cell(
    cage_face_sets: &[Vec<usize>],
    rings: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
) -> CageResult {
    let mut cages_512 = Vec::new();
    let mut cages_51262 = Vec::new();
    let mut cages_51264 = Vec::new();
    let mut cages_all = Vec::new();
    let mut cage_centers = Vec::new();
    let mut cage_types = Vec::new();

    for face_set in cage_face_sets {
        // Count faces by size
        let mut n4 = 0usize;
        let mut n5 = 0usize;
        let mut n6 = 0usize;
        for &fi in face_set {
            match rings[fi].len() {
                4 => n4 += 1,
                5 => n5 += 1,
                6 => n6 += 1,
                _ => {}
            }
        }

        // Face-count classification: all 17 known hydrate cage types
        let cage_type = match (n4, n5, n6) {
            (0, 12, 0) => "512",
            (0, 12, 2) => "51262",
            (0, 12, 3) => "51263",
            (0, 12, 4) => "51264",
            (0, 12, 5) => "51265",
            (0, 12, 6) => "51266",
            (0, 12, 8) => "51268",
            (1, 10, 2) => "4151062",
            (1, 10, 3) => "4151063",
            (1, 10, 4) => "4151064",
            (1, 10, 5) => "4151065",
            (2, 8, 1) => "425861",
            (2, 8, 2) => "425862",
            (2, 8, 3) => "425863",
            (2, 8, 4) => "425864",
            (3, 6, 3) => "435663",
            (3, 6, 4) => "435664",
            _ => "other",
        };

        // Collect cage vertices
        let mut vertex_set: HashSet<usize> = HashSet::new();
        for &fi in face_set {
            vertex_set.extend(rings[fi].iter());
        }
        let mut vertices: Vec<usize> = vertex_set.into_iter().collect();
        vertices.sort_unstable();

        // Calculate cage center
        let positions: Vec<[f64; 3]> = vertices.iter().map(|&i| waters[i].oxygen).collect();
        let center = cell.average(&positions);

        // Legacy fields for backwards compatibility
        match cage_type {
            "512" => cages_512.push(vertices.clone()),
            "51262" => cages_51262.push(vertices.clone()),
            "51264" => cages_51264.push(vertices.clone()),
            _ => {}
        }
        cages_all.push(vertices);
        cage_centers.push(center);
        cage_types.push(cage_type.to_string());
    }

    CageResult {
        cages_512,
        cages_51262,
        cages_51264,
        cages_all,
        cage_centers,
        cage_types,
        incomplete_cages: vec![], // populated by second pass in detect_cages
    }
}

// =========================================================================
// Incomplete Cage Detection (Second Pass)
// =========================================================================

/// Classify a polyhedron following TRACE's SEC_cage_identify logic.
///
/// Returns:
///  2 = SEC (every edge in exactly 2 faces, every vertex in exactly 3 faces, F-E+V=2)
///  1 = non-SEC (F-E+V=2, no vertex with count=1, no edge<2, but has vertex>3)
///  0 = IC (F-E+V=2, no vertex with count=1, no edge<2, but has vertex<3)
/// -1 = not a cage
pub(crate) fn classify_cage_quality(face_indices: &[usize], rings: &[Vec<usize>]) -> i8 {
    let f = face_indices.len();
    if f < 4 {
        return -1;
    }

    let mut edge_face_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_face_count: HashMap<usize, usize> = HashMap::new();

    for &fi in face_indices {
        let ring = &rings[fi];
        let n = ring.len();
        for i in 0..n {
            *vertex_face_count.entry(ring[i]).or_insert(0) += 1;
            let edge = make_edge(ring[i], ring[(i + 1) % n]);
            *edge_face_count.entry(edge).or_insert(0) += 1;
        }
    }

    let e = edge_face_count.len();
    let v = vertex_face_count.len();

    // Euler's formula
    if (f as i64) - (e as i64) + (v as i64) != 2 {
        return -1;
    }

    let has_one_vertex = vertex_face_count.values().any(|&c| c == 1);
    let has_low_edge = edge_face_count.values().any(|&c| c < 2);

    // Vertex with count=1 or edge<2 → not a cage
    if has_one_vertex || has_low_edge {
        return -1;
    }

    let has_low_vertex = vertex_face_count.values().any(|&c| c < 3);
    let has_high_vertex = vertex_face_count.values().any(|&c| c > 3);

    if has_low_vertex {
        return 0; // IC
    }
    if has_high_vertex {
        return 1; // non-SEC
    }
    2 // SEC
}

/// Relaxed Euler validation for incomplete cage assembly.
///
/// Matching TRACE's tolerance: `uni_vertice.size() * 3 <= total_vertex + tolerance`
/// where tolerance=2 for the IC pass.
///
/// This allows some vertices to have valency 2 instead of 3 (boundary vertices).
/// Still requires F-E+V=2, no vertex with count=1, no edge with count<2.
fn validate_relaxed(face_indices: &[usize], rings: &[Vec<usize>]) -> bool {
    let f = face_indices.len();
    if f < 4 {
        return false;
    }

    let mut edge_face_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_face_count: HashMap<usize, usize> = HashMap::new();

    for &fi in face_indices {
        let ring = &rings[fi];
        let n = ring.len();
        for i in 0..n {
            *vertex_face_count.entry(ring[i]).or_insert(0) += 1;
            let edge = make_edge(ring[i], ring[(i + 1) % n]);
            *edge_face_count.entry(edge).or_insert(0) += 1;
        }
    }

    let e = edge_face_count.len();
    let v = vertex_face_count.len();

    // Euler's formula must still hold
    if (f as i64) - (e as i64) + (v as i64) != 2 {
        return false;
    }

    // No vertex appearing in only 1 face
    if vertex_face_count.values().any(|&c| c == 1) {
        return false;
    }

    // No edge in fewer than 2 faces
    if edge_face_count.values().any(|&c| c < 2) {
        return false;
    }

    // TRACE tolerance: vertex_count * 3 <= total_face_vertex_sum + 2
    // This allows up to 2 vertices to have valency 2 instead of 3
    let total_vertex_sum: usize = vertex_face_count.values().sum();
    v * 3 <= total_vertex_sum + 2
}

/// Estimate cage type based on square/pentagon/hexagon face counts.
/// Returns (type_name, completeness) where completeness is 0.0-1.0.
fn estimate_cage_type(n4: usize, n5: usize, n6: usize) -> (String, f64) {
    // Known cage types: (name, expected_n4, expected_n5, expected_n6)
    let candidates: &[(&str, usize, usize, usize)] = &[
        ("512", 0, 12, 0),
        ("51262", 0, 12, 2),
        ("51263", 0, 12, 3),
        ("51264", 0, 12, 4),
        ("51265", 0, 12, 5),
        ("51266", 0, 12, 6),
        ("51268", 0, 12, 8),
        ("4151062", 1, 10, 2),
        ("4151063", 1, 10, 3),
        ("4151064", 1, 10, 4),
        ("4151065", 1, 10, 5),
        ("425861", 2, 8, 1),
        ("425862", 2, 8, 2),
        ("425863", 2, 8, 3),
        ("425864", 2, 8, 4),
        ("435663", 3, 6, 3),
        ("435664", 3, 6, 4),
    ];

    let mut best_type = "unknown";
    let mut best_completeness = 0.0f64;

    for &(name, exp4, exp5, exp6) in candidates {
        let expected_total = exp4 + exp5 + exp6;
        if expected_total == 0 {
            continue;
        }

        let m4 = if exp4 > 0 {
            n4.min(exp4) as f64 / exp4 as f64
        } else {
            if n4 == 0 {
                1.0
            } else {
                0.0
            }
        };
        let m5 = if exp5 > 0 {
            n5.min(exp5) as f64 / exp5 as f64
        } else {
            if n5 == 0 {
                1.0
            } else {
                0.0
            }
        };
        let m6 = if exp6 > 0 {
            n6.min(exp6) as f64 / exp6 as f64
        } else {
            if n6 == 0 {
                1.0
            } else {
                0.0
            }
        };

        let completeness = (m4 + m5 + m6) / 3.0;

        // Penalize if actual faces exceed expected
        let excess =
            (n4.saturating_sub(exp4) + n5.saturating_sub(exp5) + n6.saturating_sub(exp6)) as f64;
        let penalty = excess / expected_total as f64;
        let adjusted = (completeness - penalty * 0.5).max(0.0);

        if adjusted > best_completeness {
            best_completeness = adjusted;
            best_type = name;
        }
    }

    (best_type.to_string(), best_completeness)
}

fn build_incomplete_cage_info(
    cage_faces: &[usize],
    quality: i8,
    rings: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
) -> IncompleteCageInfo {
    let mut n4 = 0;
    let mut n5 = 0;
    let mut n6 = 0;
    let mut vertex_set = HashSet::new();
    for &face in cage_faces {
        vertex_set.extend(rings[face].iter().copied());
        match rings[face].len() {
            4 => n4 += 1,
            5 => n5 += 1,
            6 => n6 += 1,
            _ => {}
        }
    }

    let mut vertices: Vec<usize> = vertex_set.into_iter().collect();
    vertices.sort_unstable();
    let vertex_positions: Vec<[f64; 3]> = vertices
        .iter()
        .map(|&vertex| waters[vertex].oxygen)
        .collect();
    let vertex_index: HashMap<usize, usize> = vertices
        .iter()
        .enumerate()
        .map(|(local, &global)| (global, local))
        .collect();
    let mut edge_set = HashSet::new();
    for &face in cage_faces {
        let ring = &rings[face];
        for i in 0..ring.len() {
            edge_set.insert(make_edge(
                vertex_index[&ring[i]],
                vertex_index[&ring[(i + 1) % ring.len()]],
            ));
        }
    }
    let mut edges: Vec<[usize; 2]> = edge_set
        .into_iter()
        .map(|(left, right)| [left, right])
        .collect();
    edges.sort_unstable();
    let center = cell.average(&vertex_positions);
    let (estimated_type, completeness) = estimate_cage_type(n4, n5, n6);

    IncompleteCageInfo {
        vertices,
        vertex_positions,
        edges,
        center,
        n_squares: n4,
        n_pentagons: n5,
        n_hexagons: n6,
        n_faces_total: cage_faces.len(),
        quality: quality as u8,
        estimated_type,
        completeness,
    }
}

/// Second pass: find incomplete and non-SEC cages from rings not used by SEC cages.
///
/// Uses the same constraint propagation + bounded DFS but with relaxed validation
/// (TRACE tolerance T+2). Only considers rings NOT consumed by the first SEC pass.
#[allow(dead_code)]
fn find_incomplete_cages(
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    sec_used_rings: &HashSet<usize>,
    waters: &[Water],
    cell: &Cell,
) -> Vec<IncompleteCageInfo> {
    let mut incomplete_cages = Vec::new();
    let mut found_cage_keys: HashSet<Vec<usize>> = HashSet::new();

    // Try each unused ring as a seed
    for seed_ri in 0..rings.len() {
        if sec_used_rings.contains(&seed_ri) {
            continue;
        }

        let seed = vec![seed_ri];
        if let Some(cage_faces) = grow_cage_relaxed(&seed, rings, edge_to_rings, sec_used_rings) {
            let mut key = cage_faces.clone();
            key.sort();
            if found_cage_keys.contains(&key) {
                continue;
            }

            // Classify quality using TRACE's SEC_cage_identify logic
            let quality = classify_cage_quality(&cage_faces, rings);
            if quality < 0 {
                continue; // Not even an IC
            }
            // Skip SEC cages found in second pass (shouldn't happen but be safe)
            if quality == 2 {
                continue;
            }

            // Count face types
            let mut n4 = 0usize;
            let mut n5 = 0usize;
            let mut n6 = 0usize;
            for &fi in &cage_faces {
                match rings[fi].len() {
                    4 => n4 += 1,
                    5 => n5 += 1,
                    6 => n6 += 1,
                    _ => {}
                }
            }

            // Require minimum faces for a recognizable partial cage
            let total_faces = cage_faces.len();
            if total_faces < 8 {
                continue;
            }

            let (estimated_type, completeness) = estimate_cage_type(n4, n5, n6);

            // Require at least 60% completeness
            if completeness < 0.60 {
                continue;
            }

            // Collect vertices
            let mut vertex_set: HashSet<usize> = HashSet::new();
            for &fi in &cage_faces {
                vertex_set.extend(rings[fi].iter());
            }
            let mut vertices: Vec<usize> = vertex_set.into_iter().collect();
            vertices.sort_unstable();

            // Vertex positions
            let vertex_positions: Vec<[f64; 3]> =
                vertices.iter().map(|&i| waters[i].oxygen).collect();

            // Build local vertex index map for edges
            let vertex_index: HashMap<usize, usize> =
                vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();

            // Collect edges (mapped to local indices)
            let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
            for &fi in &cage_faces {
                let ring = &rings[fi];
                let n = ring.len();
                for i in 0..n {
                    let u = vertex_index[&ring[i]];
                    let v = vertex_index[&ring[(i + 1) % n]];
                    let edge = if u < v { (u, v) } else { (v, u) };
                    edge_set.insert(edge);
                }
            }
            let mut edges: Vec<[usize; 2]> = edge_set.into_iter().map(|(u, v)| [u, v]).collect();
            edges.sort_unstable(); // hash order otherwise

            // Calculate center
            let positions: Vec<[f64; 3]> = vertices.iter().map(|&i| waters[i].oxygen).collect();
            let center = cell.average(&positions);

            found_cage_keys.insert(key);
            incomplete_cages.push(IncompleteCageInfo {
                vertices,
                vertex_positions,
                edges,
                center,
                n_squares: n4,
                n_pentagons: n5,
                n_hexagons: n6,
                n_faces_total: total_faces,
                quality: quality as u8,
                estimated_type,
                completeness,
            });
        }
    }

    incomplete_cages
}

/// Grow a cage from seed rings using relaxed validation.
/// Excludes rings already consumed by SEC cages.
/// Uses same constraint propagation + bounded DFS as strict pass.
fn grow_cage_relaxed(
    seed_rings: &[usize],
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    excluded_rings: &HashSet<usize>,
) -> Option<Vec<usize>> {
    let mut faces: Vec<usize> = seed_rings.to_vec();
    let mut face_set: HashSet<usize> = faces.iter().cloned().collect();
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();

    // Initialize counts from seed rings
    for &ri in &faces {
        add_ring_to_counts(ri, rings, &mut edge_count, &mut vertex_count);
    }

    // Check initial validity (still enforce edge<=2, vertex<=3)
    if edge_count.values().any(|&c| c > 2) || vertex_count.values().any(|&c| c > 3) {
        return None;
    }

    grow_cage_relaxed_recursive(
        &mut faces,
        &mut face_set,
        &mut edge_count,
        &mut vertex_count,
        rings,
        edge_to_rings,
        excluded_rings,
    )
}

/// Recursive cage growth with relaxed validation for incomplete cage detection.
fn grow_cage_relaxed_recursive(
    faces: &mut Vec<usize>,
    face_set: &mut HashSet<usize>,
    edge_count: &mut HashMap<(usize, usize), usize>,
    vertex_count: &mut HashMap<usize, usize>,
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    excluded_rings: &HashSet<usize>,
) -> Option<Vec<usize>> {
    let mut propagated = Vec::new();
    loop {
        // Sorted, not hash order: the walk below closes the FIRST open edge that admits a
        // single candidate, so iteration order decides which cage gets grown -- and Rust
        // reseeds every HashMap, which made the cage count differ between two runs on one
        // unchanged structure.
        let mut open_edges: Vec<(usize, usize)> = edge_count
            .iter()
            .filter(|(_, &c)| c == 1)
            .map(|(&e, _)| e)
            .collect();
        open_edges.sort_unstable();

        if open_edges.is_empty() {
            // All edges closed — validate with relaxed criteria
            if validate_relaxed(faces, rings) {
                return Some(faces.clone());
            }
            for &ri in propagated.iter().rev() {
                faces.pop();
                face_set.remove(&ri);
                remove_ring_from_counts(ri, rings, edge_count, vertex_count);
            }
            return None;
        }

        if faces.len() >= MAX_CAGE_FACES {
            for &ri in propagated.iter().rev() {
                faces.pop();
                face_set.remove(&ri);
                remove_ring_from_counts(ri, rings, edge_count, vertex_count);
            }
            return None;
        }

        let mut forced_any = false;
        let mut dead_end = false;

        for &edge in &open_edges {
            if edge_count.get(&edge).copied().unwrap_or(0) != 1 {
                continue;
            }
            let candidates = get_candidates_excluding(
                edge,
                face_set,
                rings,
                edge_to_rings,
                edge_count,
                vertex_count,
                excluded_rings,
            );
            if candidates.is_empty() {
                dead_end = true;
                break;
            }
            if candidates.len() == 1 {
                let ri = candidates[0];
                faces.push(ri);
                face_set.insert(ri);
                add_ring_to_counts(ri, rings, edge_count, vertex_count);
                propagated.push(ri);
                forced_any = true;
                break;
            }
        }

        if dead_end {
            for &ri in propagated.iter().rev() {
                faces.pop();
                face_set.remove(&ri);
                remove_ring_from_counts(ri, rings, edge_count, vertex_count);
            }
            return None;
        }

        if !forced_any {
            let mut best_edge = None;
            let mut best_candidates: Vec<usize> = Vec::new();
            let mut best_count = usize::MAX;

            for &edge in &open_edges {
                if edge_count.get(&edge).copied().unwrap_or(0) != 1 {
                    continue;
                }
                let candidates = get_candidates_excluding(
                    edge,
                    face_set,
                    rings,
                    edge_to_rings,
                    edge_count,
                    vertex_count,
                    excluded_rings,
                );
                if candidates.is_empty() {
                    for &ri in propagated.iter().rev() {
                        faces.pop();
                        face_set.remove(&ri);
                        remove_ring_from_counts(ri, rings, edge_count, vertex_count);
                    }
                    return None;
                }
                if candidates.len() < best_count {
                    best_count = candidates.len();
                    best_candidates = candidates;
                    best_edge = Some(edge);
                }
            }

            if best_edge.is_none() {
                for &ri in propagated.iter().rev() {
                    faces.pop();
                    face_set.remove(&ri);
                    remove_ring_from_counts(ri, rings, edge_count, vertex_count);
                }
                return None;
            }

            for candidate in best_candidates {
                faces.push(candidate);
                face_set.insert(candidate);
                add_ring_to_counts(candidate, rings, edge_count, vertex_count);

                if let Some(result) = grow_cage_relaxed_recursive(
                    faces,
                    face_set,
                    edge_count,
                    vertex_count,
                    rings,
                    edge_to_rings,
                    excluded_rings,
                ) {
                    faces.pop();
                    face_set.remove(&candidate);
                    remove_ring_from_counts(candidate, rings, edge_count, vertex_count);
                    for &ri in propagated.iter().rev() {
                        faces.pop();
                        face_set.remove(&ri);
                        remove_ring_from_counts(ri, rings, edge_count, vertex_count);
                    }
                    return Some(result);
                }

                faces.pop();
                face_set.remove(&candidate);
                remove_ring_from_counts(candidate, rings, edge_count, vertex_count);
            }

            for &ri in propagated.iter().rev() {
                faces.pop();
                face_set.remove(&ri);
                remove_ring_from_counts(ri, rings, edge_count, vertex_count);
            }
            return None;
        }
    }
}

/// Get valid candidates for an open edge, excluding rings used by SEC cages
fn get_candidates_excluding(
    edge: (usize, usize),
    face_set: &HashSet<usize>,
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    edge_count: &HashMap<(usize, usize), usize>,
    vertex_count: &HashMap<usize, usize>,
    excluded_rings: &HashSet<usize>,
) -> Vec<usize> {
    edge_to_rings
        .get(&edge)
        .map(|r| {
            r.iter()
                .filter(|ri| !face_set.contains(ri))
                .filter(|ri| !excluded_rings.contains(ri))
                .filter(|&&ri| ring_addition_valid(ri, rings, edge_count, vertex_count))
                .cloned()
                .collect()
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_ring() {
        let ring = vec![3, 1, 2, 5, 4];
        let normalized = normalize_ring(&ring);
        assert_eq!(normalized[0], 1);
    }

    #[test]
    fn test_find_rings_simple_5cycle() {
        let adj_list = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        let mut ring_set = HashSet::new();
        let rings = find_rings_topo(&adj_list, 5, &mut ring_set);
        assert_eq!(rings.len(), 1);
    }

    #[test]
    fn test_find_rings_6cycle() {
        let adj_list = vec![
            vec![1, 5],
            vec![0, 2],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
            vec![4, 0],
        ];
        let mut ring_set = HashSet::new();
        let rings = find_rings_topo(&adj_list, 6, &mut ring_set);
        assert_eq!(rings.len(), 1);
    }

    #[test]
    fn test_euler_dodecahedron() {
        // A dodecahedron: 12 pentagonal faces, 30 edges, 20 vertices
        // F - E + V = 12 - 30 + 20 = 2
        let rings: Vec<Vec<usize>> = vec![
            vec![0, 1, 2, 3, 4],
            vec![0, 4, 9, 14, 5],
            vec![0, 5, 10, 6, 1],
            vec![1, 6, 11, 7, 2],
            vec![2, 7, 12, 8, 3],
            vec![3, 8, 13, 9, 4],
            vec![5, 14, 19, 15, 10],
            vec![6, 10, 15, 16, 11],
            vec![7, 11, 16, 17, 12],
            vec![8, 12, 17, 18, 13],
            vec![9, 13, 18, 19, 14],
            vec![15, 19, 18, 17, 16],
        ];
        let face_indices: Vec<usize> = (0..12).collect();
        assert!(validate_euler(&face_indices, &rings));
    }

    #[test]
    fn test_euler_invalid_open() {
        // Two pentagons sharing one edge — NOT a valid polyhedron
        let rings: Vec<Vec<usize>> = vec![vec![0, 1, 2, 3, 4], vec![0, 4, 5, 6, 7]];
        let face_indices = vec![0, 1];
        assert!(!validate_euler(&face_indices, &rings));
    }

    #[test]
    fn test_euler_invalid_vertex_valency() {
        // Construct faces where a vertex appears in 4 faces
        // Tetrahedron has vertices in 3 faces each but uses triangles (not pentagons)
        // Let's just verify the vertex check catches violations
        let rings: Vec<Vec<usize>> = vec![
            vec![0, 1, 2, 3, 4],
            vec![0, 4, 5, 6, 7],
            vec![0, 7, 8, 9, 1],
            vec![0, 1, 9, 10, 4], // vertex 0 now in 4 faces
        ];
        let face_indices = vec![0, 1, 2, 3];
        // Even if Euler holds, vertex valency check should fail
        assert!(!validate_euler(&face_indices, &rings));
    }

    #[test]
    fn test_pbc_closure_valid_ring() {
        // A ring that fits within the box — PBC closure sum should be ~0
        // Regular pentagon with ~0.28 nm edges, centered at (0.5, 0.5, 0.5)
        let r = 0.238; // circumradius for edge length ~0.28
        let cx = 0.5;
        let cy = 0.5;
        let waters: Vec<Water> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / 5.0;
                Water {
                    oxygen: [cx + r * angle.cos(), cy + r * angle.sin(), 0.5],
                    hydrogen1: [0.0; 3],
                    hydrogen2: [0.0; 3],
                }
            })
            .collect();
        let box_dims = [3.0, 3.0, 3.0];
        let ring = vec![0, 1, 2, 3, 4];
        assert!(is_valid_ring(&ring, &waters, &box_dims, 0.36));
    }

    #[test]
    fn test_pbc_closure_ghost_ring() {
        // Atoms spread across the box — PBC sum should be nonzero
        let waters: Vec<Water> = (0..5)
            .map(|i| Water {
                oxygen: [i as f64 * 0.28, 0.0, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            })
            .collect();
        let box_dims = [1.0, 1.0, 1.0];
        let ring = vec![0, 1, 2, 3, 4];
        // PBC sum in x = 4*0.28 + mic(0.0 - 1.12) = 1.12 + (-0.12) = 1.0 > 0.1
        assert!(!is_valid_ring(&ring, &waters, &box_dims, 0.36));
    }

    #[test]
    fn test_diagonal_check_collapsed_pentagon() {
        // A "collapsed" pentagon where non-adjacent vertices are very close
        let waters: Vec<Water> = vec![
            Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [0.28, 0.0, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [0.14, 0.05, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [0.28, 0.10, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [0.0, 0.10, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
        ];
        let box_dims = [3.0, 3.0, 3.0];
        let ring = vec![0, 1, 2, 3, 4];
        // Non-adjacent distances like d(0,2)=~0.15 < lrange5=~0.396
        assert!(!is_valid_ring(&ring, &waters, &box_dims, 0.36));
    }

    #[test]
    fn test_build_edge_ring_index() {
        let rings = vec![vec![0, 1, 2, 3, 4], vec![0, 4, 5, 6, 7]];
        let index = build_edge_ring_index(&rings);

        // Edge (0,4) should appear in both rings
        assert_eq!(index.get(&(0, 4)).unwrap().len(), 2);
        // Edge (1,2) should appear in only ring 0
        assert_eq!(index.get(&(1, 2)).unwrap().len(), 1);
    }

    #[test]
    fn test_sii_hydrate_cage_detection() {
        // Integration test: detect cages in a real sII hydrate structure
        // 222_S2.gro: 1088 waters, box 3.462 nm, 8 unit cells
        // Expected per unit cell: 16x512 + 8x51264
        // Expected for 8 unit cells: 128x512 + 64x51264 = 192 total
        let gro_content = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("Failed to read test GRO file");

        let parsed =
            crate::parser::parse_gro_content(&gro_content).expect("Failed to parse GRO file");

        assert_eq!(
            parsed.waters.len(),
            1088,
            "Expected 1088 waters in sII hydrate"
        );

        let result = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);

        let total = result.cages_512.len() + result.cages_51262.len() + result.cages_51264.len();

        eprintln!(
            "sII hydrate cage detection: 512={}, 51262={}, 51264={}, total={}",
            result.cages_512.len(),
            result.cages_51262.len(),
            result.cages_51264.len(),
            total
        );

        // sII structure: only 512 and 51264 cages, no 51262
        assert_eq!(
            result.cages_51262.len(),
            0,
            "sII should have no 51262 cages"
        );

        // 8 sII unit cells: exactly 128 small cages. This is pinned, not bounded -- a `>= 100`
        // lower bound is what let the count drift between runs without any test noticing.
        assert_eq!(
            result.cages_512.len(),
            128,
            "8 sII unit cells hold exactly 128 512 cages"
        );
        // 8 sII unit cells hold exactly 64 large cages. Pinned, like the 512 count: a lower
        // bound is what let 4 of them go missing unnoticed while the total drifted between runs.
        assert_eq!(
            result.cages_51264.len(),
            64,
            "8 sII unit cells hold exactly 64 51264 cages"
        );
        assert_eq!(
            total, 192,
            "128 512 + 64 51264 = the full sII 2x2x2 crystal"
        );
    }

    /// sI, as the counterpart to the sII test above: 4x4x4 = 64 unit cells, each holding
    /// 2 x 512 + 6 x 51262, so 128 and 384. Guards the enumeration change against the other
    /// hydrate structure -- sI was already complete, and must stay that way.
    #[test]
    fn test_si_hydrate_cage_detection() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/444_S1.gro"))
            .expect("Failed to read test GRO file");
        let parsed = crate::parser::parse_gro_content(&gro).expect("Failed to parse GRO file");

        let result = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        let total = result.cages_512.len() + result.cages_51262.len() + result.cages_51264.len();
        eprintln!(
            "sI hydrate cage detection: 512={}, 51262={}, 51264={}, total={}",
            result.cages_512.len(),
            result.cages_51262.len(),
            result.cages_51264.len(),
            total
        );

        assert_eq!(
            result.cages_512.len(),
            128,
            "64 sI cells hold 128 512 cages"
        );
        assert_eq!(
            result.cages_51262.len(),
            384,
            "64 sI cells hold 384 51262 cages"
        );
        assert_eq!(result.cages_51264.len(), 0, "sI has no 51264 cages");
        assert_eq!(total, 512, "128 + 384 = the full sI 4x4x4 crystal");
    }

    /// sH is hexagonal, so its GRO box line carries a tilt (v2x = -a/2). GenIce2's sH cell
    /// holds two unit cells; 2x2x2 of it is 16 unit cells, each 3 x 512 + 2 x 435663 +
    /// 1 x 51268. With only the three box lengths the minimum image is wrong across the
    /// tilted face and cages there are lost (41 instead of 48 512 cages).
    #[test]
    fn test_sh_triclinic_cage_detection() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_SH.gro"))
            .expect("Failed to read test GRO file");
        let parsed = crate::parser::parse_gro_content(&gro).expect("Failed to parse GRO file");
        assert_eq!(parsed.box_vectors.len(), 9, "sH box line keeps its tilt");
        let cell = Cell::from_box_values(&parsed.box_vectors).unwrap();

        let result = detect_cages_cell(&parsed.waters, &cell, 0.36, 35.0);
        let count = |t: &str| result.cage_types.iter().filter(|c| c.as_str() == t).count();
        assert_eq!(count("512"), 48, "16 sH cells hold 48 512 cages");
        assert_eq!(count("435663"), 32, "16 sH cells hold 32 435663 cages");
        assert_eq!(count("51268"), 16, "16 sH cells hold 16 51268 cages");
        assert_eq!(result.cage_types.len(), 96);
    }

    /// Cage assembly walks the open edges of a partially-grown polyhedron and closes the first
    /// one that admits a single candidate. Reading those edges out of a HashMap made that "first"
    /// depend on hash order, which Rust reseeds per HashMap instance -- so the same structure
    /// yielded 189, 190 or 191 cages depending on the run. The loose lower bounds in
    /// test_sii_hydrate_cage_detection (>= 100, >= 150) were wide enough to hide it.
    #[test]
    fn cage_detection_is_deterministic() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("Failed to read test GRO file");
        let parsed = crate::parser::parse_gro_content(&gro).expect("Failed to parse GRO file");

        let first = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        let counts = |r: &CageResult| (r.cages_512.len(), r.cages_51262.len(), r.cages_51264.len());
        for run in 1..8 {
            let again = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
            assert_eq!(
                counts(&first),
                counts(&again),
                "run {run} disagreed with run 0 on an unchanged structure"
            );
        }
    }

    #[test]
    fn test_classify_512() {
        // 12 pentagons = 512 cage
        let rings: Vec<Vec<usize>> = vec![
            vec![0, 1, 2, 3, 4],
            vec![0, 4, 9, 14, 5],
            vec![0, 5, 10, 6, 1],
            vec![1, 6, 11, 7, 2],
            vec![2, 7, 12, 8, 3],
            vec![3, 8, 13, 9, 4],
            vec![5, 14, 19, 15, 10],
            vec![6, 10, 15, 16, 11],
            vec![7, 11, 16, 17, 12],
            vec![8, 12, 17, 18, 13],
            vec![9, 13, 18, 19, 14],
            vec![15, 19, 18, 17, 16],
        ];
        let face_set = vec![(0..12).collect::<Vec<usize>>()];

        // Create dummy waters at origin (positions don't matter for classification)
        let waters: Vec<Water> = (0..20)
            .map(|_| Water {
                oxygen: [0.0, 0.0, 0.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            })
            .collect();
        let box_dims = [3.0, 3.0, 3.0];

        let result = classify_and_build_result(&face_set, &rings, &waters, &box_dims);
        assert_eq!(result.cages_512.len(), 1);
        assert_eq!(result.cages_51262.len(), 0);
        assert_eq!(result.cages_51264.len(), 0);
        assert_eq!(result.cage_types[0], "512");
    }

    #[test]
    fn trace_reference_222_s1_exact_counts() {
        // TRACE reference program (commit 927903b), example/222_S1/inputH2O.gro, TRACE defaults:
        // frame 0 contains 16 x 5^12 and 48 x 5^12 6^2 standard cages.
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S1_inputH2O.gro"))
            .expect("Failed to read the TRACE 222_S1 reference fixture");
        let parsed = crate::parser::parse_gro_content(&gro)
            .expect("Failed to parse the TRACE 222_S1 reference fixture");

        let result = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        let mut type_counts: HashMap<&str, usize> = HashMap::new();
        for cage_type in &result.cage_types {
            *type_counts.entry(cage_type.as_str()).or_default() += 1;
        }

        assert_eq!(type_counts.get("512"), Some(&16));
        assert_eq!(type_counts.get("51262"), Some(&48));
        assert_eq!(result.cage_types.len(), 64);
        assert!(result.incomplete_cages.is_empty());
    }

    #[test]
    fn trace_reference_example_h2o_ring_counts() {
        // TRACE reference program (commit 927903b), default -mr 10. The native detail_cage.txt
        // reports these ring counts before cup/cage assembly for frame 0.
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/example_H2O_frame0.gro"))
            .expect("Failed to read the TRACE example_H2O reference fixture");
        let parsed = crate::parser::parse_gro_content(&gro)
            .expect("Failed to parse the TRACE example_H2O reference fixture");
        let cell = Cell::orthorhombic(parsed.box_dims).unwrap();
        let adjacency = build_hbond_network_cell(&parsed.waters, &cell, 0.36, 35.0);
        let rings = find_trace_rings_cell(&adjacency, &parsed.waters, &cell, 10, 90.0, 20.0);

        let mut counts = [0usize; 11];
        for ring in rings {
            counts[ring.len()] += 1;
        }
        assert_eq!(&counts[4..=10], &[284, 1400, 1348, 153, 35, 1, 0]);
    }

    #[test]
    fn trace_reference_example_h2o_cup_count() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/example_H2O_frame0.gro"))
            .expect("Failed to read the TRACE example_H2O reference fixture");
        let parsed = crate::parser::parse_gro_content(&gro)
            .expect("Failed to parse the TRACE example_H2O reference fixture");
        let cell = Cell::orthorhombic(parsed.box_dims).unwrap();
        let adjacency = build_hbond_network_cell(&parsed.waters, &cell, 0.36, 35.0);
        let mut rings = find_trace_rings_cell(&adjacency, &parsed.waters, &cell, 10, 90.0, 20.0);
        rings.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));
        let edge_to_rings = build_edge_ring_index(&rings);

        let mut cups = find_trace_cups(&rings, &edge_to_rings);
        assert_eq!(cups.len(), 1366);

        let complete = find_trace_cage_faces(&mut cups, &rings, 0);
        assert_eq!(complete.iter().filter(|cage| cage.quality == 2).count(), 10);
        assert_eq!(complete.iter().filter(|cage| cage.quality == 1).count(), 0);

        let incomplete = find_trace_cage_faces(&mut cups, &rings, 2);
        assert_eq!(
            incomplete.iter().filter(|cage| cage.quality == 0).count(),
            12
        );
    }

    #[test]
    fn trace_reference_example_h2o_cage_counts() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/example_H2O_frame0.gro"))
            .expect("Failed to read the TRACE example_H2O reference fixture");
        let parsed = crate::parser::parse_gro_content(&gro)
            .expect("Failed to parse the TRACE example_H2O reference fixture");
        let result = detect_cages(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        let mut type_counts: HashMap<&str, usize> = HashMap::new();
        for cage_type in &result.cage_types {
            *type_counts.entry(cage_type.as_str()).or_default() += 1;
        }

        assert_eq!(type_counts.get("512"), Some(&1));
        assert_eq!(type_counts.get("51262"), Some(&2));
        assert_eq!(type_counts.get("4151062"), Some(&6));
        assert_eq!(type_counts.get("other"), Some(&1));
        assert_eq!(result.cage_types.len(), 10);
        assert_eq!(result.incomplete_cages.len(), 12);
        assert!(result.incomplete_cages.iter().all(|cage| cage.quality == 0));
    }

    #[test]
    fn trace_quality_contract_classifies_standard_edge_saturated_cage() {
        // Tetrahedron: F=4, E=6, V=4; every edge belongs to two faces and every
        // vertex belongs to three faces. This is TRACE's SEC quality 2.
        let rings = vec![vec![0, 1, 2], vec![0, 3, 1], vec![1, 3, 2], vec![2, 3, 0]];
        assert_eq!(classify_cage_quality(&[0, 1, 2, 3], &rings), 2);
    }

    #[test]
    fn trace_quality_contract_classifies_nonstandard_edge_saturated_cage() {
        // Square pyramid: F=5, E=8, V=5 and all edges are saturated, but the
        // apex belongs to four faces. TRACE classifies that high-valency cage as nSEC.
        let rings = vec![
            vec![0, 1, 2, 3],
            vec![0, 4, 1],
            vec![1, 4, 2],
            vec![2, 4, 3],
            vec![3, 4, 0],
        ];
        assert_eq!(classify_cage_quality(&[0, 1, 2, 3, 4], &rings), 1);
    }

    #[test]
    fn trace_quality_contract_classifies_incomplete_cage() {
        // Subdivide one tetrahedron edge in both incident faces. Euler and edge
        // saturation still hold, while the inserted vertex belongs to only two faces.
        // TRACE classifies that low-valency cage as IC.
        let rings = vec![
            vec![0, 4, 1, 2],
            vec![0, 3, 1, 4],
            vec![1, 3, 2],
            vec![2, 3, 0],
        ];
        assert_eq!(classify_cage_quality(&[0, 1, 2, 3], &rings), 0);
    }

    #[test]
    fn trace_preserves_arbitrary_435663_face_signature() {
        // Native TRACE records 4(3)5(6)6(3) as a named SEC type rather than
        // restricting results to the three canonical GRADE cage families.
        let mut rings = Vec::new();
        rings.extend((0..3).map(|_| vec![0, 1, 2, 3]));
        rings.extend((0..6).map(|_| vec![0, 1, 2, 3, 4]));
        rings.extend((0..3).map(|_| vec![0, 1, 2, 3, 4, 5]));
        let face_set = vec![(0..rings.len()).collect::<Vec<_>>()];
        let waters: Vec<Water> = (0..6)
            .map(|_| Water {
                oxygen: [0.0; 3],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            })
            .collect();

        let result = classify_and_build_result(&face_set, &rings, &waters, &[3.0; 3]);
        assert_eq!(result.cage_types, vec!["435663"]);
        assert_eq!(result.cages_all.len(), 1);
        assert!(result.cages_512.is_empty());
        assert!(result.cages_51262.is_empty());
        assert!(result.cages_51264.is_empty());
    }
}
