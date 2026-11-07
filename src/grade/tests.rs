// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

use super::cages::{filled_count, group_cups_64512, reference_center};
use super::cups::Topology;
use super::faces::{link_rings, merge_overlap, remove_duplicates, Row};
use super::metric::Metric;
use super::planarity::is_planar;
use super::rings::Rings;
use super::*;

fn closed(face: &[usize]) -> Row {
    let mut row = face.to_vec();
    row.push(face[0]);
    row
}

fn line_positions(n: usize) -> Vec<[f64; 3]> {
    (0..n).map(|i| [i as f64 * 0.01, 0.0, 0.0]).collect()
}

fn dodecahedron_faces() -> Vec<Row> {
    [
        [0, 19, 3, 2, 1],
        [0, 10, 11, 18, 19],
        [0, 1, 8, 9, 10],
        [1, 2, 6, 7, 8],
        [2, 3, 4, 5, 6],
        [3, 19, 18, 17, 4],
        [4, 17, 16, 15, 5],
        [5, 15, 14, 7, 6],
        [7, 14, 13, 9, 8],
        [9, 13, 12, 11, 10],
        [11, 12, 16, 17, 18],
        [12, 13, 14, 15, 16],
    ]
    .iter()
    .map(|face| closed(face))
    .collect()
}

#[test]
fn grade_defaults_match_reference_cli() {
    assert_eq!(
        GradeParameters::default(),
        GradeParameters {
            rcut: 0.35,
            theta_deg: 45.0,
            pentagon_diagonal_delta: 0.18,
            hexagon_diagonal_delta: 0.26,
        }
    );
}

#[test]
fn grade_rejects_irregular_trace_face_signature() {
    assert_eq!(GradeCageType::from_face_counts(3, 6, 3), None);
    assert_eq!(GradeCageType::from_face_counts(0, 12, 3), None);
    assert_eq!(GradeCageType::from_face_counts(0, 12, 8), None);
}

#[test]
fn invalid_parameters_are_rejected_at_api_boundary() {
    let parameters = GradeParameters {
        rcut: 0.18,
        ..GradeParameters::default()
    };
    let error = detect_grade(&[], &[], &[1.0, 1.0, 1.0], parameters).unwrap_err();
    assert!(error.contains("greater than 0.18"));
}

#[test]
fn cup_pairing_finds_one_512_and_no_noncanonical_type() {
    let found = Rings {
        pentagons: dodecahedron_faces(),
        hexagons: Vec::new(),
    };
    let cell = Cell::orthorhombic([10.0, 10.0, 10.0]).unwrap();
    let result =
        assemble_grade_result(&found, &line_positions(20), &[], &cell, &Metric::new(&cell));
    assert_eq!(result.cages_512.len(), 1);
    assert_eq!(result.cages_512[0], (0..20).collect::<Vec<_>>());
    assert!(result.cages_51262.is_empty());
    assert!(result.cages_51264.is_empty());
}

#[test]
fn filled_count_matches_reference_unique_guest_semantics() {
    let cell = Cell::orthorhombic([2.0, 2.0, 2.0]).unwrap();
    let centers = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
    let guests = [[0.05, 0.0, 0.0], [1.0, 1.0, 1.0]];
    // One guest is close to both cages, but the reference `.xvg` column
    // counts that solute only once for a given cage family.
    assert_eq!(filled_count(&centers, &guests, &Metric::new(&cell)), 1);
}

#[test]
fn reference_center_counts_every_ring_membership_and_starts_first_cage_at_one() {
    let found = Rings {
        pentagons: vec![closed(&[0, 1, 2, 3, 4])],
        hexagons: vec![closed(&[0, 1, 2, 3, 4, 5])],
    };
    let links = link_rings(&found.pentagons, &found.hexagons);
    let topology = Topology {
        pentagons: &found.pentagons,
        hexagons: &found.hexagons,
        links: &links,
    };
    // A 14-ring row: hexagons at positions 0 and 7, pentagons elsewhere.
    let row = vec![0; 14];
    let center = reference_center(&row, true, &topology, &line_positions(6));
    // Two hexagons (sum 0.15 each) and twelve pentagons (sum 0.10 each) over
    // 2 * 6 + 12 * 5 + 1 vertex slots.
    assert!((center[0] - (2.0 * 0.15 + 12.0 * 0.10) / 73.0).abs() < 1e-12);
}

#[test]
fn planarity_rejects_a_ring_in_a_symmetry_plane_like_the_reference() {
    // sII 2x2x2 pentagon lying exactly in the plane x = z. Its second normal
    // angle has cosine 1.0000000000000002, so the reference's acos gives NaN
    // and the ring fails `theta < THETA`.
    let positions = [
        [0.965, 1.307, 0.965],
        [1.025, 1.571, 1.025],
        [0.865, 1.731, 0.865],
        [0.706, 1.571, 0.706],
        [0.766, 1.307, 0.766],
    ];
    assert!(!is_planar(&[0, 1, 2, 3, 4], &positions, 45.0));
}

#[test]
fn planarity_uses_raw_coordinates_so_a_wrapped_vertex_fails() {
    let flat = [
        [0.0, 0.0, 0.0],
        [0.27, 0.0, 0.01],
        [0.35, 0.26, 0.0],
        [0.13, 0.42, 0.01],
        [-0.08, 0.25, 0.0],
    ];
    assert!(is_planar(&[0, 1, 2, 3, 4], &flat, 45.0));
    // The same ring with one oxygen wrapped by a 2 nm box length along z.
    let mut wrapped = flat;
    wrapped[2][2] += 2.0;
    assert!(!is_planar(&[0, 1, 2, 3, 4], &wrapped, 45.0));
}

#[test]
fn hexagon_history_lag_stores_recent_repeats_again() {
    let (a, b, c) = ([0, 1, 2, 3, 4, 5], [1, 2, 3, 4, 5, 6], [2, 3, 4, 5, 6, 7]);
    let mut seen = rings::Candidates::default();
    for ring in [a, a, b, b, c, c] {
        seen.add_hexagon(ring);
    }
    // The reference stores A and B twice (via the k6 == 0 branch) and
    // re-stores B and C because the newest history entries are never scanned.
    let order: Vec<usize> = seen.hexagons.iter().map(|row| row[0]).collect();
    assert_eq!(order, [0, 0, 1, 1, 1, 2, 2]);
}

#[test]
fn duplicate_removal_compares_hexagon_and_pentagon_indices_as_one_numbering() {
    let kept = remove_duplicates(
        vec![
            vec![3, 1, 2, 3, 4, 5, 6],
            vec![7, 1, 2, 3, 4, 5, 6],
            vec![7, 6, 5, 4, 3, 2, 1],
        ],
        0,
        0,
    );
    // Row 0 repeats the index 3 and is dropped; row 2 has the key of row 1.
    assert_eq!(kept, vec![vec![7, 1, 2, 3, 4, 5, 6]]);
}

#[test]
fn hexagon_never_links_to_the_pentagon_with_its_own_index() {
    let pentagons = vec![closed(&[0, 1, 2, 3, 4]), closed(&[0, 1, 5, 6, 7])];
    let hexagons = vec![closed(&[0, 1, 8, 9, 10, 11])];
    let links = link_rings(&pentagons, &hexagons);
    assert_eq!(links.pentagon[0], vec![1]);
    // Both pentagons share the edge 0-1 with the hexagon, but index 0 is skipped.
    assert_eq!(links.hexagon_pentagon[0], vec![1]);
}

#[test]
fn four_cup_group_counts_each_shared_pentagon_beyond_the_first() {
    let cup = |hexagon: usize, sides: [usize; 6]| [vec![hexagon], sides.to_vec()].concat();
    let accepted = [
        cup(0, [1, 2, 3, 4, 5, 6]),
        cup(1, [1, 2, 3, 20, 21, 22]),
        cup(2, [4, 5, 30, 31, 32, 33]),
    ];
    // Cup 1 shares three pentagons (appended twice), cup 2 shares two (once).
    assert_eq!(group_cups_64512(&accepted).len(), 1);
    let rejected = [
        cup(0, [1, 2, 3, 4, 5, 6]),
        cup(1, [1, 2, 3, 4, 21, 22]),
        cup(2, [5, 6, 30, 31, 32, 33]),
    ];
    // Four shared pentagons append cup 1 three times, so the group has five.
    assert!(group_cups_64512(&rejected).is_empty());
}

#[test]
fn merge_overlap_misses_matches_in_unsorted_rings() {
    assert_eq!(merge_overlap(&[1, 2, 3], &[3, 2, 1]), 1);
    assert_eq!(merge_overlap(&[1, 2, 3], &[1, 2, 3]), 3);
}

fn crystal_counts(gro: &str) -> (usize, usize, usize) {
    let parsed = crate::parser::parse_gro_content(gro).unwrap();
    let result = detect_grade(
        &parsed.waters,
        &[],
        &parsed.box_dims,
        GradeParameters::default(),
    )
    .unwrap();
    (
        result.cages_512.len(),
        result.cages_51262.len(),
        result.cages_51264.len(),
    )
}

// Native GRADE v1.00 counts for the GenIce2 crystals below, produced with
// `tests/fixtures/reference/grade_crystals.json`. They are the reference's own
// output, not the ideal cage counts (16 + 48, 128 + 64, 54 + 162).
#[test]
fn matches_native_grade_on_222_s1_crystal() {
    let gro = include_str!("../../tests/fixtures/reference/grade_222_s1.gro");
    assert_eq!(crystal_counts(gro), (6, 44, 0));
}

#[test]
fn matches_native_grade_on_222_s2_crystal() {
    let gro = include_str!("../../tests/fixtures/reference/grade_222_s2.gro");
    assert_eq!(crystal_counts(gro), (46, 0, 24));
}

#[test]
fn matches_native_grade_on_htr_1242_crystal() {
    let gro = include_str!("../../tests/fixtures/reference/htr_1242.gro");
    assert_eq!(crystal_counts(gro), (35, 159, 0));
}

#[test]
fn matches_bundled_grade_prod_complete_oracle() {
    let manifest = include_str!("../../tests/fixtures/reference/grade_prod.json");
    assert!(manifest.contains("7b06be697cb4fc01e8ac280a12adc85ac7e5f0b15a0023c0b60e636c4f76df22"));
    let gro = include_str!("../../tests/fixtures/reference/grade_prod.gro");
    let parsed = crate::parser::parse_gro_content(gro).unwrap();
    let guest_positions: Vec<_> = parsed.guests.iter().map(|guest| guest.position).collect();

    let result = detect_grade(
        &parsed.waters,
        &guest_positions,
        &parsed.box_dims,
        GradeParameters::default(),
    )
    .unwrap();

    assert_eq!(parsed.waters.len(), 4500);
    assert_eq!(guest_positions.len(), 790);
    assert_eq!(result.cages_512.len(), 64);
    assert_eq!(result.filled_512, 63);
    assert_eq!(result.cages_51262.len(), 23);
    assert_eq!(result.filled_51262, 23);
    assert_eq!(result.cages_51264.len(), 4);
    assert_eq!(result.filled_51264, 4);
}
