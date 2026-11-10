// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Registry of the published methods the engine implements: what each needs, what it returns,
//! how the app may show it, and how it was verified. The app builds its analysis blocks from this.
//!
//! Verification tiers, strongest first (`tier` keeps the original names; `level` is the T0–T3 scale):
//! - T3 `author_code`: results match the authors' own program on the reference inputs.
//! - T2 `reference_implementation`: results match an established open implementation.
//! - T1 `paper` with fixtures: implemented from the paper and checked on cases with known answers
//!   (crystals, hand-built systems).
//! - T0 `paper` without fixtures: implemented from the paper only.

use serde::{Deserialize, Serialize};

/// Bumped when a method's parameters, outputs or numbers change, so saved recipes can notice.
pub const DESCRIPTOR_VERSION: u32 = 1;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum VerificationTier {
    AuthorCode,
    ReferenceImplementation,
    Paper,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Level {
    T0,
    T1,
    T2,
    T3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MethodKind {
    CageDetector,
    WaterClassifier,
    OrderParameter,
}

/// How the app may show a method's result.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum View {
    /// Per-frame metrics as a chart and CSV. Every method has it.
    MetricChart,
    /// Each water coloured in 3D by a per-molecule value or class.
    MoleculeColour,
    /// The 3D cage and occupancy view.
    CageView,
}

/// Which parts of the structure a method reads.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
pub struct InputNeeds {
    pub hydrogens: bool,
    pub guests: bool,
}

#[derive(Debug, Clone, Copy, Serialize, PartialEq)]
pub struct ParamSpec {
    pub name: &'static str,
    pub label: &'static str,
    pub unit: &'static str,
    pub default: f64,
    pub min: f64,
    pub max: f64,
    pub step: f64,
    /// Whole numbers only (counts).
    pub integer: bool,
}

/// One per-frame number (or one per-molecule array) a method returns.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq)]
pub struct SeriesSpec {
    pub key: &'static str,
    pub label: &'static str,
    pub unit: &'static str,
}

#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct MethodInfo {
    pub id: &'static str,
    pub name: &'static str,
    pub kind: MethodKind,
    pub citation: &'static str,
    pub doi: &'static str,
    pub tier: VerificationTier,
    /// What the tier rests on.
    pub verified_against: &'static str,
    pub descriptor_version: u32,
    pub level: Level,
    /// Licence of this implementation.
    pub licence: &'static str,
    /// The program this implementation was ported from and its licence, when it is a port.
    pub ported_from: Option<&'static str>,
    pub inputs: InputNeeds,
    pub params: &'static [ParamSpec],
    /// Per-frame numbers, shown by `MetricChart`.
    pub series: &'static [SeriesSpec],
    /// Per-molecule arrays, shown by `MoleculeColour`.
    pub molecule_values: &'static [SeriesSpec],
    pub views: &'static [View],
    /// Test inputs (under `tests/fixtures/`) the tier claim rests on.
    pub fixtures: &'static [&'static str],
    /// Details the paper leaves open and the choice made for each.
    pub spec_gaps: &'static [&'static str],
}

const fn level_of(tier: VerificationTier, fixtures: &[&str]) -> Level {
    match tier {
        VerificationTier::AuthorCode => Level::T3,
        VerificationTier::ReferenceImplementation => Level::T2,
        VerificationTier::Paper if fixtures.is_empty() => Level::T0,
        VerificationTier::Paper => Level::T1,
    }
}

const O_ONLY: InputNeeds = InputNeeds { hydrogens: false, guests: false };
const O_H: InputNeeds = InputNeeds { hydrogens: true, guests: false };

const fn param(name: &'static str, label: &'static str, unit: &'static str, default: f64, min: f64, max: f64, step: f64) -> ParamSpec {
    ParamSpec { name, label, unit, default, min, max, step, integer: false }
}

const fn count_param(name: &'static str, label: &'static str, default: f64, min: f64, max: f64) -> ParamSpec {
    ParamSpec { name, label, unit: "", default, min, max, step: 1.0, integer: true }
}

const fn series(key: &'static str, label: &'static str, unit: &'static str) -> SeriesSpec {
    SeriesSpec { key, label, unit }
}

const HBOND_PARAMS: &[ParamSpec] = &[
    param("cutoff_nm", "H-bond O–O cutoff", "nm", 0.36, 0.25, 0.45, 0.01),
    param("theta_deg", "H-bond angle", "deg", 35.0, 20.0, 45.0, 1.0),
];
const ORDER_CUTOFF: ParamSpec = param("cutoff_nm", "Neighbour cutoff", "nm", 0.35, 0.3, 0.4, 0.005);
const CAGE_SERIES: &[SeriesSpec] = &[series("n_cages", "Cages", "")];
const CAGE_VIEWS: &[View] = &[View::MetricChart, View::CageView];

const METHODS: &[MethodInfo] = &[
    MethodInfo {
        id: "trace",
        name: "TRACE",
        kind: MethodKind::CageDetector,
        citation: "Hsu & Lin, J. Chem. Theory Comput. 21 (2025) 11436",
        doi: "10.1021/acs.jctc.5c01459",
        tier: VerificationTier::AuthorCode,
        verified_against: "Port of the authors' TRACE program; exact cage counts on sI/sII crystals",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(VerificationTier::AuthorCode, &["crystals/222_S2.gro", "crystals/444_S1.gro"]),
        licence: "GPL-3.0-or-later",
        ported_from: Some("TRACE (GPL-3.0-or-later)"),
        inputs: O_H,
        params: HBOND_PARAMS,
        series: CAGE_SERIES,
        molecule_values: &[],
        views: CAGE_VIEWS,
        fixtures: &["crystals/222_S2.gro", "crystals/444_S1.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "htr",
        name: "HTR",
        kind: MethodKind::CageDetector,
        citation: "Liu et al., Nanotechnol. Rev. 11 (2022) 699",
        doi: "10.1515/ntrev-2022-0044",
        tier: VerificationTier::AuthorCode,
        verified_against: "Checked against the authors' published HTR program",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(VerificationTier::AuthorCode, &["reference/htr_1242.gro", "crystals/222_S2.gro"]),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: O_H,
        params: HBOND_PARAMS,
        series: CAGE_SERIES,
        molecule_values: &[],
        views: CAGE_VIEWS,
        fixtures: &["reference/htr_1242.gro", "crystals/222_S2.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "grade",
        name: "GRADE",
        kind: MethodKind::CageDetector,
        citation: "Mahmoudinobar & Dias, Comput. Phys. Commun. 244 (2019) 385",
        doi: "10.1016/j.cpc.2019.06.004",
        tier: VerificationTier::AuthorCode,
        verified_against: "Step-by-step port of GRADE v1.00; equal counts on six reference inputs",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(
            VerificationTier::AuthorCode,
            &["reference/grade_222_s1.gro", "reference/grade_222_s2.gro", "reference/grade_prod.gro"],
        ),
        licence: "GPL-3.0-or-later",
        ported_from: Some("GRADE v1.00 (GPL-3.0-or-later)"),
        inputs: O_ONLY,
        params: &[
            param("cutoff_nm", "H-bond O–O cutoff", "nm", 0.35, 0.25, 0.45, 0.01),
            param("theta_deg", "H-bond angle", "deg", 45.0, 20.0, 45.0, 1.0),
        ],
        series: CAGE_SERIES,
        molecule_values: &[],
        views: CAGE_VIEWS,
        fixtures: &["reference/grade_222_s1.gro", "reference/grade_222_s2.gro", "reference/grade_prod.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "chill",
        name: "CHILL+",
        kind: MethodKind::WaterClassifier,
        citation: "Nguyen & Molinero, J. Phys. Chem. B 119 (2015) 9369",
        doi: "10.1021/jp510289t",
        tier: VerificationTier::Paper,
        verified_against: "Paper definitions and documented thresholds; every water of the reference crystals is classed as hydrate",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(VerificationTier::Paper, &["crystals/222_S2.gro"]),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: O_ONLY,
        params: &[ORDER_CUTOFF],
        series: &[
            series("liquid", "Liquid / other", "fraction"),
            series("ice_ih", "Hexagonal ice", "fraction"),
            series("ice_ic", "Cubic ice", "fraction"),
            series("interfacial_ice", "Interfacial ice", "fraction"),
            series("hydrate", "Hydrate", "fraction"),
            series("interfacial_hydrate", "Interfacial hydrate", "fraction"),
        ],
        molecule_values: &[series("class", "CHILL+ class", "class")],
        views: &[View::MetricChart, View::MoleculeColour],
        fixtures: &["crystals/222_S2.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "f3f4",
        name: "F3/F4",
        kind: MethodKind::OrderParameter,
        citation: "Rodger et al., Fluid Phase Equilib. 116 (1996) 326; Moon et al., Faraday Discuss. 136 (2007) 367",
        doi: "10.1039/b618194p",
        tier: VerificationTier::Paper,
        verified_against: "Published definitions; F4 pair rule follows the CalF4 program distributed with GRADE",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(VerificationTier::Paper, &["inline: sI and sH water clusters vs F4SPC values"]),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: O_H,
        params: &[ORDER_CUTOFF],
        series: &[series("f3", "F3", ""), series("f4", "F4", "")],
        molecule_values: &[series("f3", "F3", ""), series("f4", "F4", "")],
        views: &[View::MetricChart, View::MoleculeColour],
        fixtures: &["inline: sI and sH water clusters vs F4SPC values"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "steinhardt",
        name: "Averaged Steinhardt q4/q6",
        kind: MethodKind::OrderParameter,
        citation: "Lechner & Dellago, J. Chem. Phys. 129 (2008) 114707",
        doi: "10.1063/1.2977970",
        tier: VerificationTier::ReferenceImplementation,
        verified_against: "freud Steinhardt(average=True), per-water agreement within its float32 precision",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(
            VerificationTier::ReferenceImplementation,
            &["crystals/444_S1.gro", "crystals/222_S2.gro", "guests/sI_CH4_222.gro", "reference/grade_prod.gro"],
        ),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: O_ONLY,
        params: &[ORDER_CUTOFF],
        series: &[series("qbar4", "q̄4", ""), series("qbar6", "q̄6", "")],
        molecule_values: &[series("qbar4", "q̄4", ""), series("qbar6", "q̄6", "")],
        views: &[View::MetricChart, View::MoleculeColour],
        fixtures: &["crystals/444_S1.gro", "crystals/222_S2.gro", "guests/sI_CH4_222.gro", "reference/grade_prod.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "mcg",
        name: "MCG",
        kind: MethodKind::OrderParameter,
        citation: "Barnes, Beckham, Wu & Sum, J. Chem. Phys. 140 (2014) 164506",
        doi: "10.1063/1.4871898",
        tier: VerificationTier::Paper,
        verified_against: "Paper definitions (MCG-3); hand-built pairs and fully occupied sI/sII crystals",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(VerificationTier::Paper, &["guests/sI_CH4_222.gro", "guests/sII_CH4_222.gro"]),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: InputNeeds { hydrogens: false, guests: true },
        params: &[
            param("r_guest_nm", "Guest–guest distance", "nm", 0.9, 0.5, 1.5, 0.05),
            param("r_water_nm", "Guest–water distance", "nm", 0.6, 0.3, 1.0, 0.05),
            param("phi_deg", "Cone angle", "deg", 45.0, 15.0, 90.0, 1.0),
            count_param("n_water", "Shared waters", 5.0, 1.0, 10.0),
            count_param("n_count", "MCG-n", 3.0, 1.0, 6.0),
        ],
        series: &[series("largest_cluster", "Largest MCG cluster", "guests")],
        molecule_values: &[],
        views: &[View::MetricChart],
        fixtures: &["guests/sI_CH4_222.gro", "guests/sII_CH4_222.gro"],
        spec_gaps: &[],
    },
    MethodInfo {
        id: "tetra",
        name: "Tetrahedral order",
        kind: MethodKind::OrderParameter,
        citation: "Errington & Debenedetti, Nature 409 (2001) 318; Chau & Hardwick, Mol. Phys. 93 (1998) 511",
        doi: "10.1038/35053024",
        tier: VerificationTier::ReferenceImplementation,
        verified_against: "GROMACS 2026.2 gmx order -Sg -Sk: frame means of S_g and S_k within 2e-6 on a liquid frame and three crystals",
        descriptor_version: DESCRIPTOR_VERSION,
        level: level_of(
            VerificationTier::ReferenceImplementation,
            &["crystals/example_H2O_frame0.gro", "crystals/222_S2.gro", "benchmark/Ih_332.gro", "benchmark/sH_222.gro"],
        ),
        licence: "GPL-3.0-or-later",
        ported_from: None,
        inputs: O_ONLY,
        params: &[],
        series: &[
            series("q", "Tetrahedral order q", ""),
            series("sg", "Angular order S_g", ""),
            series("sk", "Distance order S_k", ""),
        ],
        molecule_values: &[series("q", "Tetrahedral order q", "")],
        views: &[View::MetricChart, View::MoleculeColour],
        fixtures: &["crystals/example_H2O_frame0.gro", "crystals/222_S2.gro", "benchmark/Ih_332.gro", "benchmark/sH_222.gro"],
        spec_gaps: &[
            "Neighbours: the papers say the four nearest neighbours; taken as the four nearest oxygens under the minimum image, with no cutoff, ties keeping the lower index (as gmx order)",
            "S_k normalisation: the mean distance is that of the same four neighbours, giving the sum of squared deviations over 12 r_mean^2 (as gmx order)",
        ],
    },
];

pub fn list_methods() -> Vec<MethodInfo> {
    METHODS.to_vec()
}

pub fn describe_method(id: &str) -> Option<MethodInfo> {
    METHODS.iter().find(|m| m.id == id).cloned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_ids_are_unique_and_complete() {
        let m = list_methods();
        let mut ids: Vec<_> = m.iter().map(|x| x.id).collect();
        ids.sort();
        ids.dedup();
        assert_eq!(ids, ["chill", "f3f4", "grade", "htr", "mcg", "steinhardt", "tetra", "trace"]);
        assert!(m.iter().all(|x| x.doi.starts_with("10.") && !x.citation.is_empty()));
    }

    #[test]
    fn descriptors_are_complete_and_consistent() {
        for m in list_methods() {
            assert!(m.descriptor_version >= 1, "{}", m.id);
            assert!(!m.licence.is_empty(), "{}", m.id);
            // Every method shows its metrics; 3D views only where the output supports them.
            assert!(m.views.contains(&View::MetricChart), "{}", m.id);
            assert!(!m.series.is_empty(), "{}", m.id);
            if m.views.contains(&View::CageView) {
                assert_eq!(m.kind, MethodKind::CageDetector, "{}", m.id);
            }
            assert_eq!(m.views.contains(&View::MoleculeColour), !m.molecule_values.is_empty(), "{}", m.id);
            for p in m.params {
                assert!(p.min <= p.default && p.default <= p.max && p.step > 0.0, "{}.{}", m.id, p.name);
            }
            let mut names: Vec<_> = m.params.iter().map(|p| p.name).collect();
            names.sort();
            names.dedup();
            assert_eq!(names.len(), m.params.len(), "{} has duplicate params", m.id);
            // A tier above T0 must name the fixtures it rests on.
            if m.level != Level::T0 {
                assert!(!m.fixtures.is_empty(), "{}", m.id);
            }
            // Spec-gap rule: a shipped paper-only method has at most two documented open details.
            assert!(m.spec_gaps.len() <= 2, "{}", m.id);
        }
    }

    #[test]
    fn levels_follow_the_old_tiers() {
        for m in list_methods() {
            let expected = match m.tier {
                VerificationTier::AuthorCode => Level::T3,
                VerificationTier::ReferenceImplementation => Level::T2,
                VerificationTier::Paper if m.fixtures.is_empty() => Level::T0,
                VerificationTier::Paper => Level::T1,
            };
            assert_eq!(m.level, expected, "{}", m.id);
        }
        assert!(describe_method("steinhardt").is_some());
        assert!(describe_method("nope").is_none());
    }

    /// The descriptor snapshot both bindings are checked against (pytest for PyO3, the WASM build's own test suite for WASM).
    #[test]
    fn descriptor_snapshot_is_current() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/methods-descriptor.json");
        let text = format!("{}\n", serde_json::to_string_pretty(&list_methods()).unwrap());
        if std::env::var("UPDATE_DESCRIPTOR_SNAPSHOT").as_deref() == Ok("1") || !std::path::Path::new(path).exists() {
            std::fs::write(path, &text).unwrap();
        }
        assert_eq!(std::fs::read_to_string(path).unwrap(), text, "re-run with UPDATE_DESCRIPTOR_SNAPSHOT=1");
    }

    #[test]
    fn fixture_files_exist() {
        for m in list_methods() {
            for f in m.fixtures.iter().filter(|f| !f.starts_with("inline:")) {
                let path = format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), f);
                assert!(std::path::Path::new(&path).exists(), "{}: {}", m.id, f);
            }
        }
    }
}
