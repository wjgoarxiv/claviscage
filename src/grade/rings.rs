// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Pentagon and hexagon search of `ring_Finder` in the GRADE reference.

use std::collections::HashSet;

use super::faces::{remove_duplicates, Row};
use super::metric::Metric;
use super::planarity::is_planar;
use super::{GradeParameters, MIN_OO_DISTANCE};

/// Planar, de-duplicated rings in the order the reference reports them.
pub(super) struct Rings {
    pub(super) pentagons: Vec<Row>,
    pub(super) hexagons: Vec<Row>,
}

pub(super) fn find_rings(
    positions: &[[f64; 3]],
    metric: &Metric,
    parameters: GradeParameters,
) -> Rings {
    let neighbours = oxygen_neighbours(positions, metric, parameters.rcut);
    let (pentagons, hexagons) = candidate_rings(&neighbours, positions, metric, parameters);
    let planar = |rows: Vec<Row>| -> Vec<Row> {
        let kept = rows
            .into_iter()
            .filter(|row| is_planar(&row[..row.len() - 1], positions, parameters.theta_deg))
            .collect();
        remove_duplicates(kept, 1, 1)
    };
    Rings {
        pentagons: planar(pentagons),
        hexagons: planar(hexagons),
    }
}

/// Ascending oxygen neighbours with `0.18 < d <= rcut`.
fn oxygen_neighbours(positions: &[[f64; 3]], metric: &Metric, rcut: f64) -> Vec<Vec<usize>> {
    let mut neighbours = vec![Vec::new(); positions.len()];
    for i in 0..positions.len() {
        for j in 0..positions.len() {
            if i == j {
                continue;
            }
            let distance = metric.distance(&positions[i], &positions[j]);
            if distance <= rcut && distance > MIN_OO_DISTANCE {
                neighbours[i].push(j);
            }
        }
    }
    neighbours
}

/// Ring candidates with the reference's two de-duplication behaviours.
///
/// Pentagons are kept once per vertex set. For hexagons the reference tracks
/// its history with a counter `k6` that advances only when a new ring is
/// stored, while the history list also grows by one on every sighting made
/// while `k6 == 0` and by a second entry for the second distinct ring. The
/// membership scan covers `m6[0..=k6]`, so the newest few history entries are
/// never compared and a later sighting of those rings is stored again.
#[derive(Default)]
pub(super) struct Candidates {
    pentagons: Vec<Row>,
    pentagon_sets: HashSet<[usize; 5]>,
    pub(super) hexagons: Vec<Row>,
    hexagon_history: Vec<[usize; 6]>,
    hexagon_counter: usize,
}

impl Candidates {
    fn add_pentagon(&mut self, ring: [usize; 5]) {
        let mut key = ring;
        key.sort_unstable();
        if self.pentagon_sets.insert(key) {
            self.pentagons.push(closed(&ring));
        }
    }

    pub(super) fn add_hexagon(&mut self, ring: [usize; 6]) {
        let mut key = ring;
        key.sort_unstable();
        if self.hexagon_counter == 0 {
            self.hexagon_history.push(key);
            self.hexagons.push(closed(&ring));
        }
        let scanned = &self.hexagon_history[..=self.hexagon_counter];
        if !scanned.contains(&key) {
            self.hexagon_history.push(key);
            self.hexagon_counter += 1;
            self.hexagons.push(closed(&ring));
        }
    }
}

fn closed(ring: &[usize]) -> Row {
    let mut row = ring.to_vec();
    row.push(ring[0]);
    row
}

fn candidate_rings(
    neighbours: &[Vec<usize>],
    positions: &[[f64; 3]],
    metric: &Metric,
    parameters: GradeParameters,
) -> (Vec<Row>, Vec<Row>) {
    let pentagon_min = 1.6 * parameters.rcut - parameters.pentagon_diagonal_delta;
    let hexagon_min = 2.0 * parameters.rcut - parameters.hexagon_diagonal_delta;
    let apart =
        |i: usize, j: usize, minimum: f64| metric.distance(&positions[i], &positions[j]) > minimum;
    let mut found = Candidates::default();
    for a1 in 0..neighbours.len() {
        // The reference closes a ring by scanning `My_neigh[a1][0..N-1]`: its
        // dummy leading zero shifts the list and drops the largest neighbour.
        let closing = &neighbours[a1][..neighbours[a1].len().saturating_sub(1)];
        for &a2 in &neighbours[a1] {
            for &a3 in neighbours[a2].iter().filter(|&&v| v != a1) {
                for &a4 in neighbours[a3].iter().filter(|&&v| v != a2 && v != a1) {
                    for &a5 in neighbours[a4]
                        .iter()
                        .filter(|&&v| v != a3 && v != a2 && v != a1)
                    {
                        if closing.binary_search(&a5).is_ok()
                            && apart(a1, a4, pentagon_min)
                            && apart(a1, a3, pentagon_min)
                            && apart(a2, a5, pentagon_min)
                        {
                            found.add_pentagon([a1, a2, a3, a4, a5]);
                        }
                        for &a6 in &neighbours[a5] {
                            if [a1, a2, a3, a4].contains(&a6) || closing.binary_search(&a6).is_err()
                            {
                                continue;
                            }
                            if apart(a1, a4, hexagon_min)
                                && apart(a2, a5, hexagon_min)
                                && apart(a3, a6, hexagon_min)
                            {
                                found.add_hexagon([a1, a2, a3, a4, a5, a6]);
                            }
                        }
                    }
                }
            }
        }
    }
    (found.pentagons, found.hexagons)
}
