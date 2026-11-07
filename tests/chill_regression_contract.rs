// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Frozen v0.3.15 CHILL+ regression coverage.
//!
//! This is a same-engine regression fixture, not an independent reference implementation.

use clavis_wasm::{chill, parser};

#[test]
fn sii_222_frozen_chill_contract() {
    let gro = include_str!("fixtures/crystals/222_S2.gro");
    let parsed = parser::parse_gro_content(gro).expect("parse frozen sII fixture");
    let result = chill::calculate_chill(&parsed.waters, &parsed.box_dims, 0.35);

    assert_eq!(result.n_water, 1088);
    assert_eq!(result.class_counts, [0, 0, 0, 0, 1088, 0]);
    assert_eq!(result.fractions, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    assert_eq!(result.edges.len(), 2176);
    assert!(result.per_atom_class.iter().all(|class| *class == 4));
}

#[test]
fn si_222_frozen_chill_contract() {
    let gro = include_str!("fixtures/crystals/222_S1_inputH2O.gro");
    let parsed = parser::parse_gro_content(gro).expect("parse frozen sI 2x2x2 fixture");
    let result = chill::calculate_chill(&parsed.waters, &parsed.box_dims, 0.35);

    assert_eq!(result.n_water, 368);
    assert_eq!(result.class_counts, [0, 0, 0, 0, 368, 0]);
    assert_eq!(result.fractions, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    assert_eq!(result.edges.len(), 736);
    assert!(result.per_atom_class.iter().all(|class| *class == 4));
}

#[test]
fn si_444_frozen_chill_contract() {
    let gro = include_str!("fixtures/crystals/444_S1.gro");
    let parsed = parser::parse_gro_content(gro).expect("parse frozen sI 4x4x4 fixture");
    let result = chill::calculate_chill(&parsed.waters, &parsed.box_dims, 0.35);

    assert_eq!(result.n_water, 2944);
    assert_eq!(result.class_counts, [0, 0, 0, 0, 2944, 0]);
    assert_eq!(result.fractions, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    assert_eq!(result.edges.len(), 5888);
    assert!(result.per_atom_class.iter().all(|class| *class == 4));
}
