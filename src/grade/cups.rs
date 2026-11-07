// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Cup search of `cup_512_Finder` and `cup_62512_Finder`.
//!
//! A cup is a base ring followed by its lateral pentagons, as ring indices.
//! The reference mixes pentagon and hexagon indices in these comparisons; the
//! mixing is kept because it decides which cups survive.

use super::faces::{adjacent_off_base, merge_overlap, share_edge, Links, Row};

pub(super) struct Topology<'a> {
    pub(super) pentagons: &'a [Row],
    pub(super) hexagons: &'a [Row],
    pub(super) links: &'a Links,
}

impl Topology<'_> {
    fn pentagon_degree(&self, index: usize) -> usize {
        self.links.pentagon[index].len()
    }

    /// Check the overlap of consecutive rings (and, for hexagon cups, the
    /// last ring against the base) after a complete cup is found.
    fn overlaps_allowed(&self, base: &[usize], sides: &[usize], wrap: bool) -> bool {
        let mut rings: Vec<&[usize]> = vec![base];
        rings.extend(sides.iter().map(|&s| &self.pentagons[s][..5]));
        let mut pairs: Vec<(usize, usize)> = (0..rings.len() - 1).map(|i| (i, i + 1)).collect();
        if wrap {
            pairs.push((rings.len() - 1, 0));
        }
        pairs
            .into_iter()
            .all(|(a, b)| merge_overlap(rings[a], rings[b]) <= 2)
    }
}

/// `5^6` cups: a base pentagon plus five lateral pentagons.
pub(super) fn find_cups_512(topology: &Topology) -> Vec<Vec<usize>> {
    let mut cups = Vec::new();
    for first in 0..topology.pentagons.len() {
        if topology.pentagon_degree(first) > 4 {
            extend_512(topology, &mut vec![first], &mut cups);
        }
    }
    cups
}

fn extend_512(topology: &Topology, chosen: &mut Vec<usize>, cups: &mut Vec<Vec<usize>>) {
    let depth = chosen.len();
    let first = chosen[0];
    let base = &topology.pentagons[first];
    for &candidate in &topology.links.pentagon[chosen[depth - 1]] {
        let ring = &topology.pentagons[candidate];
        let eligible = match depth {
            1 => topology.pentagon_degree(candidate) > 2,
            2 => {
                !chosen.contains(&candidate)
                    && share_edge(base, ring)
                    && topology.pentagon_degree(candidate) > 2
            }
            3 | 4 => {
                !chosen.contains(&candidate)
                    && share_edge(base, ring)
                    && adjacent_off_base(&topology.pentagons[chosen[depth - 1]], ring, base)
                    && topology.pentagon_degree(candidate) > 2
            }
            _ => {
                !chosen.contains(&candidate)
                    && share_edge(&topology.pentagons[chosen[1]], ring)
                    && adjacent_off_base(&topology.pentagons[chosen[4]], ring, base)
                    && share_edge(base, ring)
            }
        };
        if !eligible {
            continue;
        }
        chosen.push(candidate);
        if depth == 5 {
            if topology.overlaps_allowed(&base[..5], &chosen[1..], false) {
                cups.push(chosen.clone());
            }
        } else {
            extend_512(topology, chosen, cups);
        }
        chosen.pop();
    }
}

/// `6^1 5^6` cups: a base hexagon plus six lateral pentagons.
pub(super) fn find_cups_62512(topology: &Topology) -> Vec<Vec<usize>> {
    let mut cups = Vec::new();
    for first in 0..topology.hexagons.len() {
        if topology.links.hexagon_pentagon[first].len() > 5 {
            extend_62512(topology, &mut vec![first], &mut cups);
        }
    }
    cups
}

fn extend_62512(topology: &Topology, chosen: &mut Vec<usize>, cups: &mut Vec<Vec<usize>>) {
    let depth = chosen.len();
    let first = chosen[0];
    let base = &topology.hexagons[first];
    let pool = if depth == 1 {
        &topology.links.hexagon_pentagon[first]
    } else {
        &topology.links.pentagon[chosen[depth - 1]]
    };
    for &candidate in pool {
        let ring = &topology.pentagons[candidate];
        let eligible = match depth {
            1 => candidate != first && topology.pentagon_degree(candidate) > 3,
            2..=5 => {
                !chosen.contains(&candidate)
                    && share_edge(base, ring)
                    && adjacent_off_base(&topology.pentagons[chosen[depth - 1]], ring, base)
                    && topology.pentagon_degree(candidate) > 3
            }
            _ => {
                !chosen.contains(&candidate)
                    && share_edge(&topology.pentagons[chosen[1]], ring)
                    && adjacent_off_base(&topology.pentagons[chosen[5]], ring, base)
                    && topology.links.hexagon_pentagon[first].contains(&candidate)
                    && share_edge(base, ring)
            }
        };
        if !eligible {
            continue;
        }
        chosen.push(candidate);
        if depth == 6 {
            if topology.overlaps_allowed(&base[..6], &chosen[1..], true) {
                cups.push(chosen.clone());
            }
        } else {
            extend_62512(topology, chosen, cups);
        }
        chosen.pop();
    }
}
