// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// clavis-cli — headless native runner over the ClavisCage analysis engine.
// Reads a .gro, runs an analysis, prints JSON.
// Build/run: cargo run --release --features cli --bin clavis-cli -- <analysis> <file.gro> [cutoff_nm]
// analysis: chill | trace | htr | grade | f3f4

use std::env;
use std::fs;
use std::process;

use clavis_wasm::geometry::Cell;
use clavis_wasm::{cage, chill, grade, htr, order_params, parser, GradeResult};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: clavis-cli <chill|trace|htr|grade|f3f4> <file.gro> [cutoff_nm]");
        process::exit(2);
    }
    let analysis = args[1].as_str();
    let path = &args[2];
    let cutoff: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.35);

    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error: cannot read {}: {}", path, e);
            process::exit(1);
        }
    };
    let parsed = match parser::parse_gro_content(&content) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("error: parse failed: {}", e);
            process::exit(1);
        }
    };

    let box3 = parsed.box_dims;
    let theta = 35.0;
    // CHILL+/F3-F4 use an O-O cutoff (default 0.35 nm); TRACE/HTR use an H-bond rcut (0.36 nm).
    let rcut = args
        .get(3)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.36);

    let grade_parameters = grade::GradeParameters {
        rcut: args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.35),
        ..grade::GradeParameters::default()
    };
    let guest_positions: Vec<_> = parsed.guests.iter().map(|guest| guest.position).collect();
    let known = ["chill", "trace", "htr", "grade", "f3f4"];
    if !known.contains(&analysis) {
        eprintln!("error: unknown analysis '{}' (expected chill|trace|htr|grade|f3f4)", analysis);
        process::exit(2);
    }

    // A tilted box (9 values on the GRO box line) takes the engine's general-cell path.
    let json: Result<String, String> = if parsed.box_vectors.len() == 9 {
        let cell = match Cell::from_box_values(&parsed.box_vectors) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("error: invalid box: {}", e);
                process::exit(1);
            }
        };
        match analysis {
            "chill" => to_json(&chill::calculate_chill_cell(&parsed.waters, &cell, cutoff)),
            "trace" => to_json(&cage::detect_cages_cell(&parsed.waters, &cell, rcut, theta)),
            "htr" => to_json(&htr::detect_htr_rings_cell(&parsed.waters, &cell, rcut, theta)),
            "grade" => grade::detect_grade_cell(&parsed.waters, &guest_positions, &cell, grade_parameters)
                .and_then(|native| to_json(&GradeResult::from_native(native, grade_parameters))),
            _ => to_json(&order_params::calculate_f3_f4_cell(&parsed.waters, &cell, cutoff)),
        }
    } else {
        match analysis {
            "chill" => to_json(&chill::calculate_chill(&parsed.waters, &box3, cutoff)),
            "trace" => to_json(&cage::detect_cages(&parsed.waters, &box3, rcut, theta)),
            "htr" => to_json(&htr::detect_htr_rings(&parsed.waters, &box3, rcut, theta)),
            "grade" => grade::detect_grade(&parsed.waters, &guest_positions, &box3, grade_parameters)
                .and_then(|native| to_json(&GradeResult::from_native(native, grade_parameters))),
            _ => to_json(&order_params::calculate_f3_f4(&parsed.waters, &box3, cutoff)),
        }
    };

    match json {
        Ok(s) => println!("{}", s),
        Err(e) => {
            eprintln!("error: serialize failed: {}", e);
            process::exit(1);
        }
    }
}

fn to_json<T: serde::Serialize>(value: &T) -> Result<String, String> {
    serde_json::to_string_pretty(value).map_err(|e| e.to_string())
}
