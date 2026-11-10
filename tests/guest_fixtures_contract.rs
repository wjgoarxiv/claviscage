// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Guest-containing crystal fixtures (GenIce2, fully occupied CH4): water and guest counts.

use clavis_wasm::parser::parse_gro_content;

#[test]
fn guests_fixture_counts() {
    let si = parse_gro_content(include_str!("fixtures/guests/sI_CH4_222.gro")).unwrap();
    assert_eq!((si.waters.len(), si.guests.len()), (368, 64));
    let sii = parse_gro_content(include_str!("fixtures/guests/sII_CH4_222.gro")).unwrap();
    assert_eq!((sii.waters.len(), sii.guests.len()), (1088, 192));
    assert!(si.guests.iter().chain(&sii.guests).all(|g| g.name == "CH4"));
}
