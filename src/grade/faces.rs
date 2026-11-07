// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Ring rows and the edge predicates GRADE uses to link them.

use std::collections::{BTreeSet, HashMap, HashSet};

/// A ring as the reference stores it: path order with the first vertex
/// repeated at the end, so `windows(2)` yields every edge of the cycle.
pub(super) type Row = Vec<usize>;

fn edge(pair: &[usize]) -> (usize, usize) {
    if pair[0] < pair[1] {
        (pair[0], pair[1])
    } else {
        (pair[1], pair[0])
    }
}

fn has_edge(row: &[usize], wanted: (usize, usize)) -> bool {
    row.windows(2).any(|pair| edge(pair) == wanted)
}

/// `compare`: the two rings share at least one edge.
pub(super) fn share_edge(a: &[usize], b: &[usize]) -> bool {
    a.windows(2).any(|pair| has_edge(b, edge(pair)))
}

/// `compare_adjacant`: `a` and `b` share an edge and no shared edge lies on
/// the base ring.
pub(super) fn adjacent_off_base(a: &[usize], b: &[usize], base: &[usize]) -> bool {
    let mut shared = false;
    for pair in a.windows(2) {
        let wanted = edge(pair);
        if has_edge(b, wanted) {
            if has_edge(base, wanted) {
                return false;
            }
            shared = true;
        }
    }
    shared
}

/// `std::set_intersection` size over ranges that are not sorted, which is how
/// the cup finders call it. The merge walk can miss matches, so the result
/// is a lower bound of the true overlap.
pub(super) fn merge_overlap(a: &[usize], b: &[usize]) -> usize {
    let (mut i, mut j, mut common) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            i += 1;
        } else if b[j] < a[i] {
            j += 1;
        } else {
            common += 1;
            i += 1;
            j += 1;
        }
    }
    common
}

/// Ring-to-ring neighbour lists (`My_neigh_ring5`, `My_neigh_ring6_ring5`).
pub(super) struct Links {
    pub(super) pentagon: Vec<Vec<usize>>,
    pub(super) hexagon_pentagon: Vec<Vec<usize>>,
}

/// Link rings that share an edge, each list sorted and free of repeats.
///
/// The reference skips the pair when the two loop indices are equal. That is
/// right for pentagon-pentagon links, but for a hexagon `i` and a pentagon `k`
/// it compares indices of two different lists, so pentagon `i` is never a
/// neighbour of hexagon `i`. The same skip is applied here.
pub(super) fn link_rings(pentagons: &[Row], hexagons: &[Row]) -> Links {
    let mut by_edge: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (index, ring) in pentagons.iter().enumerate() {
        for pair in ring.windows(2) {
            by_edge.entry(edge(pair)).or_default().push(index);
        }
    }
    let neighbours = |ring: &Row, own_index: usize| -> Vec<usize> {
        let found: BTreeSet<usize> = ring
            .windows(2)
            .filter_map(|pair| by_edge.get(&edge(pair)))
            .flatten()
            .copied()
            .filter(|&other| other != own_index)
            .collect();
        found.into_iter().collect()
    };
    Links {
        pentagon: pentagons
            .iter()
            .enumerate()
            .map(|(i, ring)| neighbours(ring, i))
            .collect(),
        hexagon_pentagon: hexagons
            .iter()
            .enumerate()
            .map(|(i, ring)| neighbours(ring, i))
            .collect(),
    }
}

/// `remove_duplicates_map` / `remove_duplicates_map_rings`: keep the first row
/// for each sorted-value key, preserving order.
///
/// A row is dropped when its sorted key holds more than `tolerated` equal
/// adjacent pairs, counted from position `first_pair`. Keys mix whatever
/// indices the row holds, so a hexagon index equal to a pentagon index in the
/// same cup or cage counts as a repeat.
pub(super) fn remove_duplicates(rows: Vec<Row>, first_pair: usize, tolerated: usize) -> Vec<Row> {
    let mut seen = HashSet::new();
    rows.into_iter()
        .filter(|row| {
            let mut key = row.clone();
            key.sort_unstable();
            let repeats = (first_pair..key.len().saturating_sub(1))
                .filter(|&j| key[j] == key[j + 1])
                .count();
            repeats <= tolerated && seen.insert(key)
        })
        .collect()
}
