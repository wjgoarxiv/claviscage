// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Cage assembly (`cage_Finder`, `cage_Finder_64512`) and guest counting.

use std::collections::BTreeSet;

use super::cups::Topology;
use super::faces::{remove_duplicates, Row};
use super::metric::Metric;
use super::FILLED_DISTANCE;

/// A cage as the rows of its rings: hexagon and pentagon indices in the
/// reference's layout (12 pentagons, 14 rings with hexagons at 0 and 7, or
/// 28 rings with four leading hexagons).
pub(super) type CageRings = Vec<usize>;

/// `5^12` and `5^12 6^2` cages: two cups whose lateral pentagons touch each
/// other exactly twice per ring.
pub(super) fn pair_cups(cups: &[Row], topology: &Topology) -> Vec<CageRings> {
    let mut rows = Vec::new();
    for (i, cup) in cups.iter().enumerate() {
        for other in &cups[i + 1..] {
            let matched = cup[1..].iter().all(|&ring| {
                let neighbours = &topology.links.pentagon[ring];
                other[1..].iter().filter(|s| neighbours.contains(s)).count() == 2
            });
            if matched {
                rows.push([cup.as_slice(), other.as_slice()].concat());
            }
        }
    }
    remove_duplicates(rows, 0, 0)
}

/// `5^12 6^4` cages: a cup together with the later cups it shares at least
/// two lateral pentagons with. The reference appends such a cup once for every
/// shared pentagon beyond the first, and accepts exactly four entries.
pub(super) fn group_cups_64512(cups: &[Row]) -> Vec<CageRings> {
    let mut rows = Vec::new();
    for (i, cup) in cups.iter().enumerate() {
        let mut group = vec![i];
        for (j, other) in cups.iter().enumerate().skip(i + 1) {
            let shared = cup[1..].iter().filter(|s| other[1..].contains(s)).count();
            group.extend(std::iter::repeat_n(j, shared.saturating_sub(1)));
        }
        if group.len() == 4 {
            let hexagons = group.iter().map(|&g| cups[g][0]);
            let pentagons = group.iter().flat_map(|&g| cups[g][1..].iter().copied());
            rows.push(hexagons.chain(pentagons).collect());
        }
    }
    rows
}

fn ring_vertices<'a>(row: &CageRings, position: usize, topology: &'a Topology) -> &'a [usize] {
    let hexagon = match row.len() {
        14 => position == 0 || position == 7,
        28 => position < 4,
        _ => false,
    };
    if hexagon {
        &topology.hexagons[row[position]][..6]
    } else {
        &topology.pentagons[row[position]][..5]
    }
}

/// Unique oxygens of a cage, ascending.
pub(super) fn cage_vertices(row: &CageRings, topology: &Topology) -> Vec<usize> {
    let set: BTreeSet<usize> = (0..row.len())
        .flat_map(|p| ring_vertices(row, p, topology).iter().copied())
        .collect();
    set.into_iter().collect()
}

/// The reference's cage "centre of mass": raw (unwrapped) oxygen coordinates
/// summed once per ring membership, divided by that count. The counter starts
/// at one for the first cage of each family, as `oxygen_counter` does.
pub(super) fn reference_center(
    row: &CageRings,
    first_cage: bool,
    topology: &Topology,
    positions: &[[f64; 3]],
) -> [f64; 3] {
    let mut sum = [0.0; 3];
    let mut count = usize::from(first_cage);
    for position in 0..row.len() {
        for &vertex in ring_vertices(row, position, topology) {
            for (axis, total) in sum.iter_mut().enumerate() {
                *total += positions[vertex][axis];
            }
            count += 1;
        }
    }
    sum.map(|total| total / count as f64)
}

/// Unique solutes within 0.20 nm of any reference centre.
pub(super) fn filled_count(centers: &[[f64; 3]], guests: &[[f64; 3]], metric: &Metric) -> usize {
    guests
        .iter()
        .filter(|guest| {
            centers
                .iter()
                .any(|center| metric.distance(guest, center) <= FILLED_DISTANCE)
        })
        .count()
}
