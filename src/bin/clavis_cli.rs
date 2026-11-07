// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// clavis-cli — headless native runner over the ClavisCage analysis engine.
// Reads a .gro, runs an analysis, prints JSON.
// Build/run: cargo run --release --features cli --bin clavis-cli -- <analysis> <file.gro> [cutoff_nm]
// analysis: chill | trace | htr | grade | f3f4

use std::env;
use std::fs;
use std::process;

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

    let json: Result<String, String> = match analysis {
        "chill" => serde_json::to_string_pretty(&chill::calculate_chill(&parsed.waters, &box3, cutoff))
            .map_err(|e| e.to_string()),
        "trace" => serde_json::to_string_pretty(&cage::detect_cages(&parsed.waters, &box3, rcut, theta))
            .map_err(|e| e.to_string()),
        "htr" => serde_json::to_string_pretty(&htr::detect_htr_rings(&parsed.waters, &box3, rcut, theta))
            .map_err(|e| e.to_string()),
        "grade" => {
            let parameters = grade::GradeParameters {
                rcut: args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.35),
                ..grade::GradeParameters::default()
            };
            let guest_positions: Vec<_> = parsed.guests.iter().map(|guest| guest.position).collect();
            grade::detect_grade(&parsed.waters, &guest_positions, &box3, parameters)
                .map(|native| GradeResult::from_native(native, parameters))
                .and_then(|result| serde_json::to_string_pretty(&result).map_err(|e| e.to_string()))
        }
        "f3f4" => serde_json::to_string_pretty(&order_params::calculate_f3_f4(&parsed.waters, &box3, cutoff))
            .map_err(|e| e.to_string()),
        other => {
            eprintln!("error: unknown analysis '{}' (expected chill|trace|htr|grade|f3f4)", other);
            process::exit(2);
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
