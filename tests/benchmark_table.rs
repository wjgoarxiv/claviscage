// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! The benchmark table: every registered method on its own fixtures and on the GenIce2 crystals
//! in `tests/fixtures/benchmark/` (see PROVENANCE.md there), written to `benchmarks/results.json`.
//! Rows whose answer is fixed by the crystal's construction must match exactly; the rest are
//! recorded so a change in any number shows up in the diff.

use std::collections::BTreeMap;

use clavis_wasm::geometry::Cell;
use clavis_wasm::methods;
use clavis_wasm::parser::{parse_gro_content, ParseResult};
use clavis_wasm::{cage, grade, htr, run};
use serde::Serialize;

const BENCHMARK: &[&str] = &[
    "benchmark/sI_222.gro",
    "benchmark/sII_222.gro",
    "benchmark/sH_222.gro",
    "benchmark/Ih_332.gro",
    "benchmark/Ic_222.gro",
];

#[derive(Serialize)]
struct Row {
    method: &'static str,
    fixture: String,
    waters: usize,
    /// Cage counts by face signature `n4-n5-n6`, or the method's series values.
    observed: BTreeMap<String, f64>,
    expected: Option<BTreeMap<String, f64>>,
    /// `exact`, `mismatch`, `recorded` (no exact answer) or `error: ...`.
    verdict: String,
}

/// Rows that miss their exact answer for a known reason. The table must show exactly these.
const KNOWN_MISMATCH: &[(&str, &str, &str)] = &[
    ("grade", "benchmark/sI_222.gro", GRADE_BOUNDARY),
    ("grade", "benchmark/sII_222.gro", GRADE_BOUNDARY),
    ("grade", "benchmark/sH_222.gro", GRADE_BOUNDARY),
];
const GRADE_BOUNDARY: &str = "GRADE v1.00 tests ring planarity on raw coordinates without the minimum image, so rings \
    that straddle the box edge are dropped and the count depends on where the edge falls; the port keeps this";

struct Fixture {
    data: ParseResult,
    /// The full box line (3 or 9 values); `ParseResult::box_dims` keeps only the diagonal.
    cell: Vec<f64>,
}

fn read(rel: &str) -> Fixture {
    let path = format!("{}/tests/fixtures/{rel}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let last = text.lines().rev().find(|l| !l.trim().is_empty()).unwrap();
    let cell = last.split_whitespace().map(|v| v.parse().unwrap()).collect();
    Fixture { data: parse_gro_content(&text).unwrap(), cell }
}

fn counts(pairs: &[(&str, usize)]) -> BTreeMap<String, f64> {
    pairs.iter().map(|(k, v)| (k.to_string(), *v as f64)).collect()
}

/// Answers fixed by construction. Detectors that type only 512/51262/51264 are held to those
/// types; the sH medium and large cages are outside them.
fn expected(method: &str, fixture: &str) -> Option<BTreeMap<String, f64>> {
    let three = |a: usize, b: usize, c: usize| counts(&[("0-12-0", a), ("0-12-2", b), ("0-12-4", c)]);
    let cages = match fixture {
        "crystals/444_S1.gro" => Some(three(128, 384, 0)),
        "crystals/222_S2.gro" | "benchmark/sII_222.gro" => Some(three(128, 0, 64)),
        "benchmark/sI_222.gro" => Some(three(16, 48, 0)),
        "benchmark/sH_222.gro" if method == "htr" => Some(counts(&[("0-12-0", 48), ("3-6-3", 32), ("0-12-8", 16)])),
        "benchmark/sH_222.gro" => Some(three(48, 0, 0)),
        _ => None,
    };
    let one = |key: &str| Some(BTreeMap::from([(key.to_string(), 1.0)]));
    match method {
        "trace" | "htr" | "grade" => cages,
        "chill" => match fixture {
            "crystals/222_S2.gro" | "benchmark/sI_222.gro" | "benchmark/sII_222.gro" => one("hydrate"),
            "benchmark/Ih_332.gro" => one("ice_ih"),
            "benchmark/Ic_222.gro" => one("ice_ic"),
            _ => None,
        },
        "mcg" => match fixture {
            "guests/sI_CH4_222.gro" => Some(counts(&[("largest_cluster", 64)])),
            "guests/sII_CH4_222.gro" => Some(counts(&[("largest_cluster", 192)])),
            _ => None,
        },
        _ => None,
    }
}

fn label(n4: usize, n5: usize, n6: usize) -> String {
    format!("{n4}-{n5}-{n6}")
}

fn param(id: &str, name: &str) -> f64 {
    methods::describe_method(id).unwrap().params.iter().find(|p| p.name == name).unwrap().default
}

fn observe(id: &str, f: &Fixture) -> Result<BTreeMap<String, f64>, String> {
    let d = &f.data;
    let cell = Cell::from_box_values(&f.cell)?;
    let hbond = || (param(id, "cutoff_nm"), param(id, "theta_deg"));
    let tally = |labels: Vec<String>| {
        let mut out = BTreeMap::new();
        for l in labels {
            *out.entry(l).or_insert(0.0) += 1.0;
        }
        out
    };
    match id {
        "trace" => {
            let (rcut, theta) = hbond();
            let r = cage::detect_cages_cell(&d.waters, &cell, rcut, theta);
            let mut out = counts(&[("0-12-0", r.cages_512.len()), ("0-12-2", r.cages_51262.len()), ("0-12-4", r.cages_51264.len())]);
            out.retain(|_, v| *v > 0.0);
            Ok(out)
        }
        "htr" => {
            let (rcut, theta) = hbond();
            let r = htr::detect_htr_rings_cell(&d.waters, &cell, rcut, theta);
            Ok(r.cage_types.iter().map(|ct| (ct.cage_type.clone(), ct.count as f64)).collect())
        }
        "grade" => {
            let (rcut, theta) = hbond();
            let p =grade::GradeParameters { rcut, theta_deg: theta, ..grade::GradeParameters::default() };
            let guests: Vec<_> = d.guests.iter().map(|g| g.position).collect();
            let r = grade::detect_grade_cell(&d.waters, &guests, &cell, p)?;
            Ok(tally(
                std::iter::repeat_n(label(0, 12, 0), r.cages_512.len())
                    .chain(std::iter::repeat_n(label(0, 12, 2), r.cages_51262.len()))
                    .chain(std::iter::repeat_n(label(0, 12, 4), r.cages_51264.len()))
                    .collect(),
            ))
        }
        _ => run::run_method(id, &d.waters, &d.guests, &f.cell, &BTreeMap::new(), false)
            .map(|out| out.series)
            .map_err(|e| e.to_string()),
    }
}

/// Exact when every key either side names agrees; a key one side omits counts as 0.
fn verdict(observed: &BTreeMap<String, f64>, expected: &Option<BTreeMap<String, f64>>) -> String {
    match expected {
        None => "recorded".into(),
        Some(want) => {
            let value = |m: &BTreeMap<String, f64>, k: &String| m.get(k).copied().unwrap_or(0.0);
            let same = want.keys().chain(observed.keys()).all(|k| value(observed, k) == value(want, k));
            if same { "exact".into() } else { "mismatch".into() }
        }
    }
}

#[test]
fn benchmark_table() {
    let mut rows = Vec::new();
    for m in methods::list_methods() {
        let own: Vec<&str> = m.fixtures.iter().copied().filter(|f| f.ends_with(".gro")).collect();
        // A method's own fixture may already be a benchmark crystal; list it once.
        let extra = BENCHMARK.iter().copied().filter(|b| !own.contains(b));
        for fixture in own.iter().copied().chain(extra) {
            let f = read(fixture);
            let want = expected(m.id, fixture);
            let known = KNOWN_MISMATCH.iter().find(|(id, fx, _)| *id == m.id && *fx == fixture);
            let (observed, verdict) = match observe(m.id, &f) {
                Ok(obs) => {
                    let v = match (verdict(&obs, &want), known) {
                        (v, Some((_, _, why))) if v == "mismatch" => format!("known mismatch: {why}"),
                        (v, _) => v,
                    };
                    (obs, v)
                }
                Err(e) => (BTreeMap::new(), format!("error: {e}")),
            };
            rows.push(Row { method: m.id, fixture: fixture.into(), waters: f.data.waters.len(), observed, expected: want, verdict });
        }
    }

    let table = serde_json::json!({
        "engine": env!("CARGO_PKG_VERSION"),
        "descriptor_versions": methods::list_methods().iter().map(|m| (m.id, m.descriptor_version)).collect::<BTreeMap<_, _>>(),
        "rows": rows,
    });
    let dir = format!("{}/benchmarks", env!("CARGO_MANIFEST_DIR"));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(format!("{dir}/results.json"), serde_json::to_string_pretty(&table).unwrap() + "\n").unwrap();

    let bad: Vec<String> = rows
        .iter()
        .filter(|r| r.verdict == "mismatch" || (r.expected.is_some() && r.verdict.starts_with("error")))
        .map(|r| format!("{} on {}: observed {:?}, expected {:?} ({})", r.method, r.fixture, r.observed, r.expected, r.verdict))
        .collect();
    assert!(bad.is_empty(), "exact-by-construction rows failed:\n{}", bad.join("\n"));
    // A known mismatch that starts matching is news too: the list must be updated with it.
    let known: Vec<_> = rows.iter().filter(|r| r.verdict.starts_with("known")).map(|r| (r.method, r.fixture.as_str())).collect();
    assert_eq!(known, KNOWN_MISMATCH.iter().map(|(m, f, _)| (*m, *f)).collect::<Vec<_>>());
    // The sI 512 and sII 192 counts of the reference crystals are rows of the table.
    let total = |method: &str, fixture: &str| -> f64 {
        rows.iter().find(|r| r.method == method && r.fixture == fixture).map_or(0.0, |r| r.observed.values().sum())
    };
    assert_eq!(total("trace", "crystals/444_S1.gro"), 512.0);
    assert_eq!(total("trace", "crystals/222_S2.gro"), 192.0);
}
