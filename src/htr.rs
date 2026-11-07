// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
//! HTR (Hydrogen-bonded Topological Ring) detection
//!
//! Finds primitive 4/5/6-membered rings in the H-bond network
//! and assembles them into cage-like polyhedra.
//!
//! References: Liu et al., Nanotechnol. Rev. 11 (2022) 699 (HTR);
//! Shi et al., J. Phys.: Condens. Matter 36 (2024) 365901 (HTR+).
//! Implements HTR as described in these papers; results are checked against the
//! authors' published program (https://gitee.com/TxxLewis/HTR).

use crate::cage::{classify_cage_quality, is_valid_ring_cell, validate_euler};
use crate::geometry::Cell;
use crate::{CageTypeCount, HtrCageInfo, HtrResult, HtrRingInfo, IncompleteCageInfo, Water};
use std::collections::{HashMap, HashSet};

/// Fixed O--O cutoff used by the reference HTR implementation (3.5 A).
const HTR_OXYGEN_CUTOFF_NM: f64 = 0.35;
/// The reference stores at most six neighbors for each oxygen.
const HTR_MAX_NEIGHBORS: usize = 6;

/// Maximum faces in any known hydrate cage (5^12·6^8 = 20)
const MAX_CAGE_FACES: usize = 20;

/// Detect primitive rings (4/5/6-membered) and assemble into cages
pub fn detect_htr_rings(waters: &[Water], box_dims: &[f64; 3], rcut: f64, theta: f64) -> HtrResult {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    let mut result = detect_htr_rings_cell(waters, &cell, rcut, theta);
    result.box_dims = box_dims.to_vec();
    result
}

/// Detect primitive rings and assemble cages using a general periodic cell.
pub fn detect_htr_rings_cell(waters: &[Water], cell: &Cell, _rcut: f64, _theta: f64) -> HtrResult {
    let n_water = waters.len();
    let box_dims = cell.vectors();

    if n_water < 4 {
        return HtrResult {
            n_water,
            n_4ring: 0,
            n_5ring: 0,
            n_6ring: 0,
            rings: vec![],
            cages: vec![],
            cage_types: vec![],
            box_dims: box_dims.iter().flatten().copied().collect(),
            incomplete_cages: vec![],
        };
    }

    // The reference HTR topology is an oxygen proximity graph.  Hydrogen
    // orientations and the TRACE H-bond angle criterion are not part of HTR.
    let adj_list = build_htr_topology_cell(waters, cell);

    // Find primitive rings of size 4, 5, 6 with geometric validation
    let mut ring_set: HashSet<Vec<usize>> = HashSet::new();
    let mut all_rings: Vec<Vec<usize>> = Vec::new();

    for size in [4, 5, 6] {
        let rings = find_primitive_rings_cell(
            &adj_list,
            size,
            &mut ring_set,
            waters,
            cell,
            HTR_OXYGEN_CUTOFF_NM,
        );
        all_rings.extend(rings);
    }

    // Count by size
    let n_4ring = all_rings.iter().filter(|r| r.len() == 4).count();
    let n_5ring = all_rings.iter().filter(|r| r.len() == 5).count();
    let n_6ring = all_rings.iter().filter(|r| r.len() == 6).count();

    // Build ring info with geometry
    let rings: Vec<HtrRingInfo> = all_rings
        .iter()
        .map(|ring| {
            let vertex_positions: Vec<Vec<f64>> =
                ring.iter().map(|&v| waters[v].oxygen.to_vec()).collect();

            let positions: Vec<[f64; 3]> = ring.iter().map(|&v| waters[v].oxygen).collect();
            let center = cell.average(&positions);

            // Edges are consecutive pairs in the ring
            let n = ring.len();
            let edges: Vec<Vec<usize>> = (0..n).map(|i| vec![i, (i + 1) % n]).collect();

            HtrRingInfo {
                size: ring.len(),
                vertices: ring.clone(),
                center: center.to_vec(),
                vertex_positions,
                edges,
            }
        })
        .collect();

    // Assemble rings into cages (first pass - strict SEC)
    let (cages, sec_used_rings) = assemble_cages_with_tracking(&all_rings, waters, cell);

    // Count cage types
    let mut type_counts: HashMap<String, usize> = HashMap::new();
    for cage in &cages {
        *type_counts.entry(cage.cage_type.clone()).or_insert(0) += 1;
    }
    let mut cage_types: Vec<CageTypeCount> = type_counts
        .into_iter()
        .map(|(cage_type, count)| CageTypeCount { cage_type, count })
        .collect();
    cage_types.sort_by(|a, b| a.cage_type.cmp(&b.cage_type)); // hash order otherwise: reported list reshuffles between runs

    // Second pass: find incomplete cages from unused rings
    let incomplete_cages = find_incomplete_cages_htr(&all_rings, &sec_used_rings, waters, cell);

    HtrResult {
        n_water,
        n_4ring,
        n_5ring,
        n_6ring,
        rings,
        cages,
        cage_types,
        box_dims: box_dims.iter().flatten().copied().collect(),
        incomplete_cages,
    }
}

fn build_htr_topology_cell(waters: &[Water], cell: &Cell) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); waters.len()];
    for i in 0..waters.len() {
        for j in (i + 1)..waters.len() {
            if cell.distance(&waters[i].oxygen, &waters[j].oxygen) <= HTR_OXYGEN_CUTOFF_NM {
                if adjacency[i].len() < HTR_MAX_NEIGHBORS {
                    adjacency[i].push(j);
                }
                if adjacency[j].len() < HTR_MAX_NEIGHBORS {
                    adjacency[j].push(i);
                }
            }
        }
    }
    adjacency
}

/// Find all primitive rings of given size.
/// A primitive ring has no chord (no shortcut between non-adjacent vertices).
/// Also applies geometric validation (PBC closure + diagonal distance checks).
#[allow(dead_code)]
fn find_primitive_rings(
    adj_list: &[Vec<usize>],
    ring_size: usize,
    ring_set: &mut HashSet<Vec<usize>>,
    waters: &[Water],
    box_dims: &[f64; 3],
    rcut: f64,
) -> Vec<Vec<usize>> {
    let cell = Cell::orthorhombic(*box_dims).expect("orthorhombic box must be non-degenerate");
    find_primitive_rings_cell(adj_list, ring_size, ring_set, waters, &cell, rcut)
}

fn find_primitive_rings_cell(
    adj_list: &[Vec<usize>],
    ring_size: usize,
    ring_set: &mut HashSet<Vec<usize>>,
    waters: &[Water],
    cell: &Cell,
    rcut: f64,
) -> Vec<Vec<usize>> {
    let n = adj_list.len();
    let mut raw_rings = Vec::new();

    for start in 0..n {
        let mut path = vec![start];
        let mut visited = HashSet::new();
        visited.insert(start);

        find_rings_dfs(
            start,
            start,
            &mut path,
            &mut visited,
            adj_list,
            ring_size,
            &mut raw_rings,
            ring_set,
        );
    }

    // Filter for primitive rings (no chords) AND geometric validity.
    // NOTE (intentional divergence from TRACE): HTR restricts to PRIMITIVE rings — the
    // topologically-correct ring set for ring-network analysis — whereas TRACE
    // (cage.rs::find_rings_topo) does not apply the primitive filter. As a result, on
    // thermally disordered / near-dissociation systems HTR is the more conservative
    // estimator and may report fewer large (hexagon-bearing) cages than TRACE. This is by
    // design, not a bug; keep the two policies distinct.
    raw_rings
        .into_iter()
        .filter(|ring| is_primitive(ring, adj_list) && is_valid_ring_cell(ring, waters, cell, rcut))
        .collect()
}

fn find_rings_dfs(
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
        // Check if we can close the ring back to start
        if adj_list[current].contains(&start) {
            let normalized = normalize_ring(path);
            if !ring_set.contains(&normalized) {
                ring_set.insert(normalized.clone());
                rings.push(normalized);
            }
        }
        return;
    }

    if path.len() > target_size {
        return;
    }

    for &next in &adj_list[current] {
        if next == start && path.len() < target_size {
            continue;
        }
        if visited.contains(&next) {
            continue;
        }
        // Pruning: only proceed if next > start to avoid duplicates
        if path.len() == 1 && next < start {
            continue;
        }

        visited.insert(next);
        path.push(next);
        find_rings_dfs(
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

/// Check if a ring is primitive (no chords).
/// A chord is an edge between two non-adjacent vertices in the ring.
fn is_primitive(ring: &[usize], adj_list: &[Vec<usize>]) -> bool {
    let n = ring.len();
    for i in 0..n {
        for j in i + 2..n {
            // Skip the edge between first and last vertex (they ARE adjacent)
            if i == 0 && j == n - 1 {
                continue;
            }
            // Check if non-adjacent vertices are connected (chord)
            if adj_list[ring[i]].contains(&ring[j]) {
                return false;
            }
        }
    }
    true
}

// =========================================================================
// Edge helpers
// =========================================================================

/// Create a canonical edge (sorted pair)
#[inline]
fn make_edge(u: usize, v: usize) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

/// Build index: edge -> list of ring indices containing that edge
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

// =========================================================================
// Cage Assembly via Constraint Propagation + Bounded DFS
// =========================================================================

/// Assemble rings into cage-like polyhedra using constraint propagation + bounded DFS.
/// Also returns the set of ring indices used in SEC cages for the second pass.
fn assemble_cages_with_tracking(
    rings: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
) -> (Vec<HtrCageInfo>, HashSet<usize>) {
    if rings.is_empty() {
        return (vec![], HashSet::new());
    }

    let edge_to_rings = build_edge_ring_index(rings);

    let mut found_cages: Vec<Vec<usize>> = Vec::new();
    let mut found_cage_keys: HashSet<Vec<usize>> = HashSet::new();

    for seed_ri in 0..rings.len() {
        // Every closure this ring participates in, not just the first. A face is shared by two
        // cages, so stopping at the first closure reports one of them and hides the other.
        for cage_faces in grow_cage_deterministic(&[seed_ri], rings, &edge_to_rings) {
            let mut key = cage_faces.clone();
            key.sort();
            if found_cage_keys.insert(key) {
                found_cages.push(cage_faces);
            }
        }
    }

    // Track rings used in valid SEC cages
    let mut sec_used: HashSet<usize> = HashSet::new();
    let mut cages = Vec::new();
    for face_set in &found_cages {
        if validate_euler(face_set, rings) {
            sec_used.extend(face_set.iter());
        }
        if let Some(cage) = build_cage_info(face_set, rings, waters, cell) {
            cages.push(cage);
        }
    }
    (cages, sec_used)
}

/// Second pass: find incomplete cages from unused rings (HTR variant).
///
/// Uses the same relaxed validation as TRACE's cage.rs second pass.
fn find_incomplete_cages_htr(
    rings: &[Vec<usize>],
    sec_used_rings: &HashSet<usize>,
    waters: &[Water],
    cell: &Cell,
) -> Vec<IncompleteCageInfo> {
    if rings.is_empty() {
        return vec![];
    }

    let edge_to_rings = build_edge_ring_index(rings);
    let mut incomplete_cages = Vec::new();
    let mut found_cage_keys: HashSet<Vec<usize>> = HashSet::new();

    for seed_ri in 0..rings.len() {
        if sec_used_rings.contains(&seed_ri) {
            continue;
        }

        let seed = vec![seed_ri];
        if let Some(cage_faces) =
            grow_cage_relaxed_htr(&seed, rings, &edge_to_rings, sec_used_rings)
        {
            let mut key = cage_faces.clone();
            key.sort();
            if found_cage_keys.contains(&key) {
                continue;
            }

            let quality = classify_cage_quality(&cage_faces, rings);
            if quality < 0 || quality == 2 {
                continue;
            }

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

            let total_faces = cage_faces.len();
            if total_faces < 8 {
                continue;
            }

            let (estimated_type, completeness) = estimate_cage_type_htr(n4, n5, n6);
            if completeness < 0.60 {
                continue;
            }

            let mut vertex_set: HashSet<usize> = HashSet::new();
            for &fi in &cage_faces {
                vertex_set.extend(rings[fi].iter());
            }
            let mut vertices: Vec<usize> = vertex_set.into_iter().collect();
            vertices.sort_unstable();

            let vertex_positions: Vec<[f64; 3]> =
                vertices.iter().map(|&i| waters[i].oxygen).collect();

            let vertex_index: HashMap<usize, usize> =
                vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();

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
            edges.sort_unstable(); // hash order otherwise: the cage mesh's edge list reshuffles

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

/// Estimate cage type based on square/pentagon/hexagon face counts (same logic as cage.rs)
fn estimate_cage_type_htr(n4: usize, n5: usize, n6: usize) -> (String, f64) {
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

/// Relaxed Euler validation for HTR incomplete cages (same as cage.rs)
fn validate_relaxed_htr(face_indices: &[usize], rings: &[Vec<usize>]) -> bool {
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

    if (f as i64) - (e as i64) + (v as i64) != 2 {
        return false;
    }
    if vertex_face_count.values().any(|&c| c == 1) {
        return false;
    }
    if edge_face_count.values().any(|&c| c < 2) {
        return false;
    }

    let total_vertex_sum: usize = vertex_face_count.values().sum();
    v * 3 <= total_vertex_sum + 2
}

/// Grow a cage from seed rings using relaxed validation, excluding SEC-used rings.
fn grow_cage_relaxed_htr(
    seed_rings: &[usize],
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
    excluded_rings: &HashSet<usize>,
) -> Option<Vec<usize>> {
    let mut faces: Vec<usize> = seed_rings.to_vec();
    let mut face_set: HashSet<usize> = faces.iter().cloned().collect();
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();

    for &ri in &faces {
        add_ring_to_counts(ri, rings, &mut edge_count, &mut vertex_count);
    }

    if edge_count.values().any(|&c| c > 2) || vertex_count.values().any(|&c| c > 3) {
        return None;
    }

    grow_cage_relaxed_recursive_htr(
        &mut faces,
        &mut face_set,
        &mut edge_count,
        &mut vertex_count,
        rings,
        edge_to_rings,
        excluded_rings,
    )
}

/// Recursive relaxed cage growth for HTR.
fn grow_cage_relaxed_recursive_htr(
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
            if validate_relaxed_htr(faces, rings) {
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
            let candidates = get_candidates_excluding_htr(
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
                let candidates = get_candidates_excluding_htr(
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

                if let Some(result) = grow_cage_relaxed_recursive_htr(
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

/// Get valid candidates excluding SEC-used rings (HTR variant)
fn get_candidates_excluding_htr(
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

/// Grow a cage from seed rings using constraint propagation + bounded DFS.
///
/// Phase 1 (Propagation): Iteratively close open edges that have exactly 1
/// valid candidate ring. This is deterministic and fast.
///
/// Phase 2 (Branching): When all remaining open edges have 2+ candidates,
/// pick the most-constrained edge (fewest candidates) and try each candidate
/// via recursive DFS with backtracking.
///
/// Returns the first valid cage found, or None.
fn grow_cage_deterministic(
    seed_rings: &[usize],
    rings: &[Vec<usize>],
    edge_to_rings: &HashMap<(usize, usize), Vec<usize>>,
) -> Vec<Vec<usize>> {
    let mut faces: Vec<usize> = seed_rings.to_vec();
    let mut face_set: HashSet<usize> = faces.iter().cloned().collect();
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();

    for &ri in &faces {
        add_ring_to_counts(ri, rings, &mut edge_count, &mut vertex_count);
    }

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
        // Sorted, not hash order -- see the note in cage.rs: hash order made the cage count
        // differ between two runs on one unchanged structure.
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

        // Phase 1: constraint propagation -- close any edge with exactly one candidate
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
                break; // counts changed -- rescan
            }
        }
        if forced_any {
            continue;
        }

        // Phase 2: nothing forced. Branch on the most constrained edge, and take EVERY
        // candidate -- the first one that closes is not the only cage here.
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
// Cage Info Construction
// =========================================================================

/// Build HtrCageInfo from a set of face (ring) indices
fn build_cage_info(
    face_set: &[usize],
    rings: &[Vec<usize>],
    waters: &[Water],
    cell: &Cell,
) -> Option<HtrCageInfo> {
    let ring_indices = face_set.to_vec();

    // Count ring sizes
    let n4 = ring_indices
        .iter()
        .filter(|&&ri| rings[ri].len() == 4)
        .count();
    let n5 = ring_indices
        .iter()
        .filter(|&&ri| rings[ri].len() == 5)
        .count();
    let n6 = ring_indices
        .iter()
        .filter(|&&ri| rings[ri].len() == 6)
        .count();

    let cage_type = format!("{}-{}-{}", n4, n5, n6);

    // Collect unique vertices
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &ri in &ring_indices {
        vertex_set.extend(rings[ri].iter());
    }
    let mut vertices: Vec<usize> = vertex_set.into_iter().collect();
    vertices.sort_unstable();

    // Vertex positions and a PBC-aware center owned by the same cell used for HTR topology.
    let positions: Vec<[f64; 3]> = vertices.iter().map(|&v| waters[v].oxygen).collect();
    let center = cell.average(&positions);
    let vertex_positions: Vec<Vec<f64>> =
        positions.iter().map(|position| position.to_vec()).collect();

    // Build local vertex index map
    let vertex_index: HashMap<usize, usize> =
        vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();

    // Collect edges (mapped to local indices)
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    for &ri in &ring_indices {
        let ring = &rings[ri];
        let n = ring.len();
        for i in 0..n {
            let u = vertex_index[&ring[i]];
            let v = vertex_index[&ring[(i + 1) % n]];
            let edge = if u < v { (u, v) } else { (v, u) };
            edge_set.insert(edge);
        }
    }
    let mut edges: Vec<Vec<usize>> = edge_set.into_iter().map(|(u, v)| vec![u, v]).collect();
    edges.sort_unstable(); // hash order otherwise

    Some(HtrCageInfo {
        cage_type,
        ring_indices,
        vertices,
        center,
        vertex_positions,
        edges,
        n4,
        n5,
        n6,
    })
}

/// Classify hydrate structure type based on cage composition (GRADE algorithm).
///
/// Reference: Mahmoudinobar & Dias, Comput. Phys. Commun. (2019)
///
/// Structure types:
/// - sI:  2×5^12 + 6×5^12·6^2 per unit cell (ratio ~1:3)
/// - sII: 16×5^12 + 8×5^12·6^4 per unit cell (ratio ~2:1)
/// - sH:  3×5^12 + 2×4^3·5^6·6^3 + 1×5^12·6^8 per unit cell
pub fn classify_structure(
    n_512: usize,
    n_51262: usize,
    n_51264: usize,
    n_435663: usize,
    n_51268: usize,
) -> (String, f64) {
    let total = n_512 + n_51262 + n_51264 + n_435663 + n_51268;
    if total == 0 {
        return ("amorphous".to_string(), 0.0);
    }

    let has_51262 = n_51262 > 0;
    let has_51264 = n_51264 > 0;
    let has_435663 = n_435663 > 0;
    let has_51268 = n_51268 > 0;

    // sH detection: requires 4^3·5^6·6^3 (medium) and/or 5^12·6^8 (large)
    // Ideal sH ratio: 3×5^12 : 2×4^3·5^6·6^3 : 1×5^12·6^8
    if has_435663 || has_51268 {
        let sh_cages = n_435663 + n_51268;
        if sh_cages > 0 && n_512 > 0 {
            // Calculate confidence based on ratio deviation from ideal 3:2:1
            let total_sh = n_512 + n_435663 + n_51268;
            let r_512 = n_512 as f64 / total_sh as f64; // ideal: 0.5 (3/6)
            let r_435663 = n_435663 as f64 / total_sh as f64; // ideal: 0.333 (2/6)
            let r_51268 = n_51268 as f64 / total_sh as f64; // ideal: 0.167 (1/6)

            let dev =
                ((r_512 - 0.5).abs() + (r_435663 - 0.333).abs() + (r_51268 - 0.167).abs()) / 3.0;
            let confidence = (1.0 - dev * 2.0).max(0.2).min(1.0);

            // If also has sI/sII cages, it might be mixed
            if has_51262 || has_51264 {
                let non_sh = n_51262 + n_51264;
                if non_sh as f64 > total as f64 * 0.3 {
                    return ("mixed".to_string(), 0.5);
                }
            }
            return ("sH".to_string(), confidence);
        }
        // Has sH-specific cages but no 5^12 → still likely sH
        if sh_cages > 0 {
            return ("sH".to_string(), 0.4);
        }
    }

    if has_51262 && !has_51264 {
        // Only 5^12 and 5^12·6^2 → sI
        // Ideal sI ratio: 5^12·6^2 / 5^12 = 3.0
        let confidence = if n_512 > 0 {
            let ratio = n_51262 as f64 / n_512 as f64;
            (1.0 - ((ratio - 3.0).abs() / 3.0)).max(0.1).min(1.0)
        } else {
            0.5 // Only large cages, still likely sI
        };
        ("sI".to_string(), confidence)
    } else if has_51264 && !has_51262 {
        // Only 5^12 and 5^12·6^4 → sII
        // Ideal sII ratio: 5^12 / 5^12·6^4 = 2.0
        let confidence = if n_51264 > 0 {
            let ratio = n_512 as f64 / n_51264 as f64;
            (1.0 - ((ratio - 2.0).abs() / 2.0)).max(0.1).min(1.0)
        } else {
            0.5
        };
        ("sII".to_string(), confidence)
    } else if has_51262 && has_51264 {
        // Both present → likely sH or mixed
        ("mixed".to_string(), 0.5)
    } else {
        // Only 5^12 cages
        if n_512 > 0 {
            ("unknown".to_string(), 0.3)
        } else {
            ("amorphous".to_string(), 0.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_htr_uses_fixed_oxygen_cutoff_without_hydrogen_angles() {
        // The reference HTR program builds its topology from O--O distances only,
        // using a fixed 3.5 A cutoff.  This square crosses both periodic boundaries;
        // coincident hydrogen coordinates deliberately cannot satisfy an angle test.
        let waters = vec![
            Water {
                oxygen: [2.85, 2.85, 1.5],
                hydrogen1: [2.85, 2.85, 1.5],
                hydrogen2: [2.85, 2.85, 1.5],
            },
            Water {
                oxygen: [0.15, 2.85, 1.5],
                hydrogen1: [0.15, 2.85, 1.5],
                hydrogen2: [0.15, 2.85, 1.5],
            },
            Water {
                oxygen: [0.15, 0.15, 1.5],
                hydrogen1: [0.15, 0.15, 1.5],
                hydrogen2: [0.15, 0.15, 1.5],
            },
            Water {
                oxygen: [2.85, 0.15, 1.5],
                hydrogen1: [2.85, 0.15, 1.5],
                hydrogen2: [2.85, 0.15, 1.5],
            },
        ];

        let result = detect_htr_rings(&waters, &[3.0, 3.0, 3.0], 0.10, 0.0);

        assert_eq!(result.n_4ring, 1);
        assert_eq!(result.rings[0].vertices, vec![0, 1, 2, 3]);
    }

    #[test]
    fn reference_htr_caps_each_oxygen_at_six_neighbors() {
        let waters: Vec<Water> = (0..8)
            .map(|i| Water {
                oxygen: [1.0 + i as f64 * 0.01, 1.0, 1.0],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            })
            .collect();
        let cell = Cell::orthorhombic([3.0, 3.0, 3.0]).unwrap();

        let adjacency = build_htr_topology_cell(&waters, &cell);

        assert!(adjacency.iter().all(|neighbors| neighbors.len() <= 6));
    }

    #[test]
    fn matches_bundled_reference_htr_1242_complete_oracle() {
        let manifest = include_str!("../tests/fixtures/reference/htr_1242.json");
        assert!(
            manifest.contains("7a7872f3a76644acb361e0dca6d62ade173ae7275ad6e92594c07fd4733ab122")
        );
        let gro = include_str!("../tests/fixtures/reference/htr_1242.gro");
        let parsed = crate::parser::parse_gro_content(&gro).unwrap();

        let result = detect_htr_rings(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        let actual: Vec<(String, usize)> = result
            .cage_types
            .iter()
            .map(|count| (count.cage_type.clone(), count.count))
            .collect();

        assert_eq!(result.n_water, 1242);
        assert_eq!(result.n_4ring, 0);
        assert_eq!(result.n_5ring, 1296);
        assert_eq!(result.n_6ring, 162);
        assert_eq!(result.rings.len(), 1458);
        assert_eq!(result.cages.len(), 216);
        assert_eq!(
            actual,
            vec![("0-12-0".to_string(), 54), ("0-12-2".to_string(), 162)]
        );

        let hbonds =
            crate::hbond::build_hbond_network(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        assert_eq!(hbonds.iter().map(Vec::len).sum::<usize>() / 2, 2484);
    }

    /// Same defect as cage.rs: the open-edge walk read its edges out of a HashMap, so which edge
    /// got closed first depended on hash order and the cage count wobbled run to run.
    #[test]
    fn htr_detection_is_deterministic() {
        let gro = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("Failed to read test GRO file");
        let parsed = crate::parser::parse_gro_content(&gro).expect("Failed to parse GRO file");

        let types = |r: &HtrResult| -> Vec<(String, usize)> {
            r.cage_types
                .iter()
                .map(|t| (t.cage_type.clone(), t.count))
                .collect()
        };

        let first = detect_htr_rings(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        for run in 1..8 {
            let again = detect_htr_rings(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
            assert_eq!(
                first.cages.len(),
                again.cages.len(),
                "run {run} found a different number of cages on an unchanged structure"
            );
            assert_eq!(
                types(&first),
                types(&again),
                "run {run} disagreed on the cage types"
            );
            assert_eq!(
                first.incomplete_cages.len(),
                again.incomplete_cages.len(),
                "run {run} disagreed on the incomplete cages"
            );
        }
    }

    #[test]
    fn test_is_primitive() {
        // 5-cycle: 0-1-2-3-4 with no chords
        let adj_list = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        let ring = vec![0, 1, 2, 3, 4];
        assert!(is_primitive(&ring, &adj_list));
    }

    #[test]
    fn test_is_not_primitive() {
        // 5-cycle with a chord: 0-2 edge
        let adj_list = vec![
            vec![1, 2, 4], // 0 connected to 1, 2, 4
            vec![0, 2],
            vec![0, 1, 3],
            vec![2, 4],
            vec![3, 0],
        ];
        let ring = vec![0, 1, 2, 3, 4];
        assert!(!is_primitive(&ring, &adj_list));
    }

    #[test]
    fn test_find_primitive_rings_simple() {
        // Simple 5-cycle graph with no chords
        // Place waters in a regular pentagon with ~0.28nm edges
        let r = 0.238;
        let waters: Vec<Water> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * (i as f64) / 5.0;
                Water {
                    oxygen: [1.5 + r * angle.cos(), 1.5 + r * angle.sin(), 1.5],
                    hydrogen1: [0.0; 3],
                    hydrogen2: [0.0; 3],
                }
            })
            .collect();
        let box_dims = [3.0, 3.0, 3.0];
        let adj_list = vec![vec![1, 4], vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 0]];
        let mut ring_set = HashSet::new();
        let rings = find_primitive_rings(&adj_list, 5, &mut ring_set, &waters, &box_dims, 0.36);
        assert_eq!(rings.len(), 1);
    }

    #[test]
    fn test_find_primitive_4_rings() {
        // Square graph: 0-1-2-3-0
        // Place waters in a square with ~0.28nm edges
        let waters: Vec<Water> = vec![
            Water {
                oxygen: [1.5, 1.5, 1.5],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [1.78, 1.5, 1.5],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [1.78, 1.78, 1.5],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
            Water {
                oxygen: [1.5, 1.78, 1.5],
                hydrogen1: [0.0; 3],
                hydrogen2: [0.0; 3],
            },
        ];
        let box_dims = [3.0, 3.0, 3.0];
        let adj_list = vec![vec![1, 3], vec![0, 2], vec![1, 3], vec![2, 0]];
        let mut ring_set = HashSet::new();
        let rings = find_primitive_rings(&adj_list, 4, &mut ring_set, &waters, &box_dims, 0.36);
        assert_eq!(rings.len(), 1);
    }

    #[test]
    fn test_classify_si() {
        let (stype, conf) = classify_structure(2, 6, 0, 0, 0);
        assert_eq!(stype, "sI");
        assert!(conf > 0.5);
    }

    #[test]
    fn test_classify_sii() {
        let (stype, conf) = classify_structure(16, 0, 8, 0, 0);
        assert_eq!(stype, "sII");
        assert!(conf > 0.5);
    }

    #[test]
    fn test_classify_amorphous() {
        let (stype, _) = classify_structure(0, 0, 0, 0, 0);
        assert_eq!(stype, "amorphous");
    }

    #[test]
    fn test_classify_mixed() {
        let (stype, _) = classify_structure(5, 3, 2, 0, 0);
        assert_eq!(stype, "mixed");
    }

    #[test]
    fn test_classify_sh() {
        // sH: 3×5^12 + 2×4^3·5^6·6^3 + 1×5^12·6^8
        let (stype, conf) = classify_structure(3, 0, 0, 2, 1);
        assert_eq!(stype, "sH");
        assert!(conf > 0.5, "sH confidence should be > 0.5, got {}", conf);

        // Partial sH: has medium cage only
        let (stype2, _) = classify_structure(2, 0, 0, 1, 0);
        assert_eq!(stype2, "sH");

        // Partial sH: has large cage only
        let (stype3, _) = classify_structure(3, 0, 0, 0, 1);
        assert_eq!(stype3, "sH");
    }

    #[test]
    fn test_normalize_ring() {
        let ring = vec![3, 1, 2, 5, 4];
        let normalized = normalize_ring(&ring);
        assert_eq!(normalized[0], 1);
    }

    #[test]
    fn test_sii_hydrate_htr() {
        let gro_content = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/crystals/222_S2.gro"))
            .expect("Failed to read test GRO file");
        let parsed =
            crate::parser::parse_gro_content(&gro_content).expect("Failed to parse GRO file");
        let result = detect_htr_rings(&parsed.waters, &parsed.box_dims, 0.36, 35.0);
        eprintln!(
            "HTR: 4-rings={}, 5-rings={}, 6-rings={}, cages={}",
            result.n_4ring,
            result.n_5ring,
            result.n_6ring,
            result.cages.len()
        );
        for ct in &result.cage_types {
            eprintln!("  {} x {}", ct.count, ct.cage_type);
        }
        // Ideal sII 2x2x2: 128 dodecahedra (0-12-0) and 64 hexakaidecahedra (0-12-4).
        assert_eq!((result.n_4ring, result.n_5ring, result.n_6ring), (0, 1152, 128));
        assert_eq!(result.cages.len(), 192);
        let count = |label: &str| {
            result.cage_types.iter().find(|ct| ct.cage_type == label).map_or(0, |ct| ct.count)
        };
        assert_eq!(count("0-12-0"), 128);
        assert_eq!(count("0-12-4"), 64);
    }
}
