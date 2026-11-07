// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! GRO file parser

use crate::{Guest, Water};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParseResult {
    pub waters: Vec<Water>,
    pub guests: Vec<Guest>,
    pub box_dims: [f64; 3],
}

/// Parse GRO file content into waters, guests, and box dimensions
pub fn parse_gro_content(content: &str) -> Result<ParseResult, String> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.len() < 3 {
        return Err("GRO file too short".to_string());
    }

    // Line 1: title (ignored)
    // Line 2: number of atoms
    let n_atoms: usize = lines[1].trim().parse().map_err(|_| "Invalid atom count")?;

    if lines.len() < n_atoms + 3 {
        return Err("GRO file truncated".to_string());
    }

    let mut waters: Vec<Water> = Vec::new();
    let mut guests: Vec<Guest> = Vec::new();

    // Temporary storage for building water molecules
    let mut current_water_id: Option<i32> = None;
    let mut current_oxygen: Option<[f64; 3]> = None;
    let mut current_hydrogens: Vec<[f64; 3]> = Vec::new();

    // Temporary storage for building a guest molecule: (res_num, res_name, atom positions).
    // The guest position is the centroid of all its atoms, finalized (MIC-aware) once the
    // box is known.
    let mut current_guest: Option<(i32, String, Vec<[f64; 3]>)> = None;
    let mut completed_guests: Vec<(String, Vec<[f64; 3]>)> = Vec::new();

    // Parse atoms (lines 3 to n_atoms+2)
    for i in 2..(2 + n_atoms) {
        let line = lines[i];
        if line.len() < 44 {
            continue;
        }

        // GRO format (fixed width):
        // Columns 1-5: residue number
        // Columns 6-10: residue name
        // Columns 11-15: atom name
        // Columns 16-20: atom number
        // Columns 21-28: x position
        // Columns 29-36: y position
        // Columns 37-44: z position

        let res_num: i32 = line[0..5].trim().parse().unwrap_or(0);
        let res_name = line[5..10].trim();
        let atom_name = line[10..15].trim();

        let x: f64 = line[20..28].trim().parse().unwrap_or(0.0);
        let y: f64 = line[28..36].trim().parse().unwrap_or(0.0);
        let z: f64 = line[36..44].trim().parse().unwrap_or(0.0);
        let pos = [x, y, z];

        // A pending guest residue ends as soon as an atom of a different residue appears.
        if let Some((gnum, gname, _)) = &current_guest {
            if *gnum != res_num || gname.as_str() != res_name {
                if let Some((_, name, positions)) = current_guest.take() {
                    if !positions.is_empty() {
                        completed_guests.push((name, positions));
                    }
                }
            }
        }

        // Detect water molecules
        if is_water_residue(res_name) {
            if is_oxygen_atom(atom_name) {
                // Save previous water if complete
                if let Some(_) = current_water_id {
                    if let Some(oxy) = current_oxygen {
                        if current_hydrogens.len() >= 2 {
                            waters.push(Water {
                                oxygen: oxy,
                                hydrogen1: current_hydrogens[0],
                                hydrogen2: current_hydrogens[1],
                            });
                        }
                    }
                }

                // Start new water
                current_water_id = Some(res_num);
                current_oxygen = Some(pos);
                current_hydrogens.clear();
            } else if is_hydrogen_atom(atom_name) {
                if current_water_id == Some(res_num) && current_hydrogens.len() < 2 {
                    current_hydrogens.push(pos);
                }
            }
        } else if is_guest_residue(res_name) {
            // Accumulate every atom of the guest residue; the guest position is the centroid.
            match &mut current_guest {
                Some((_, _, positions)) => positions.push(pos),
                None => current_guest = Some((res_num, res_name.to_string(), vec![pos])),
            }
        }
    }

    // Save last water
    if let Some(oxy) = current_oxygen {
        if current_hydrogens.len() >= 2 {
            waters.push(Water {
                oxygen: oxy,
                hydrogen1: current_hydrogens[0],
                hydrogen2: current_hydrogens[1],
            });
        }
    }

    // Save last guest
    if let Some((_, name, positions)) = current_guest.take() {
        if !positions.is_empty() {
            completed_guests.push((name, positions));
        }
    }

    // Parse box dimensions (last line)
    let box_line = lines[2 + n_atoms];
    let box_parts: Vec<f64> = box_line
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();

    let box_dims = if box_parts.len() >= 3 {
        [box_parts[0], box_parts[1], box_parts[2]]
    } else {
        return Err("Invalid box dimensions".to_string());
    };

    // Finalize guest positions with MIC-aware centroids now that the box is known.
    for (name, positions) in completed_guests {
        guests.push(Guest {
            name,
            position: mic_centroid(&positions, &box_dims),
        });
    }

    Ok(ParseResult {
        waters,
        guests,
        box_dims,
    })
}

/// Count the number of frames in a multi-frame GRO file
pub fn count_frames(content: &str) -> usize {
    if content.is_empty() {
        return 0;
    }

    let lines: Vec<&str> = content.lines().collect();
    let mut frame_count = 0;
    let mut idx = 0;

    while idx < lines.len() {
        // Need at least: title, atom_count, box line = 3 lines minimum
        if idx + 2 >= lines.len() {
            break;
        }

        // Skip title line
        idx += 1;

        // Read atom count
        if idx >= lines.len() {
            break;
        }
        let n_atoms: usize = match lines[idx].trim().parse() {
            Ok(n) => n,
            Err(_) => break,
        };
        idx += 1;

        // Skip atom lines + box line
        let lines_to_skip = n_atoms + 1;
        if idx + lines_to_skip > lines.len() {
            break;
        }
        idx += lines_to_skip;

        frame_count += 1;
    }

    frame_count
}

/// Parse a specific frame from a multi-frame GRO file
pub fn parse_gro_frame(content: &str, frame_idx: usize) -> Result<ParseResult, String> {
    if content.is_empty() {
        return Err("Empty content".to_string());
    }

    let lines: Vec<&str> = content.lines().collect();
    let mut current_frame = 0;
    let mut idx = 0;

    // Scan to the target frame
    while current_frame < frame_idx {
        if idx + 2 >= lines.len() {
            return Err(format!("Frame {} not found (only {} frames)", frame_idx, current_frame));
        }

        // Skip title
        idx += 1;

        // Read atom count
        if idx >= lines.len() {
            return Err(format!("Frame {} not found (truncated)", frame_idx));
        }
        let n_atoms: usize = lines[idx].trim().parse()
            .map_err(|_| format!("Invalid atom count at frame {}", current_frame))?;
        idx += 1;

        // Skip atom lines + box line
        let lines_to_skip = n_atoms + 1;
        if idx + lines_to_skip > lines.len() {
            return Err(format!("Frame {} truncated", current_frame));
        }
        idx += lines_to_skip;

        current_frame += 1;
    }

    // Now extract the target frame
    if idx >= lines.len() {
        return Err(format!("Frame {} not found (only {} frames)", frame_idx, current_frame));
    }

    let frame_start = idx;

    // Title line
    idx += 1;

    // Atom count line
    if idx >= lines.len() {
        return Err(format!("Frame {} truncated (no atom count)", frame_idx));
    }
    let n_atoms: usize = lines[idx].trim().parse()
        .map_err(|_| format!("Invalid atom count at frame {}", frame_idx))?;
    idx += 1;

    // Atom lines + box line
    let frame_end = idx + n_atoms + 1;
    if frame_end > lines.len() {
        return Err(format!("Frame {} truncated (expected {} atoms)", frame_idx, n_atoms));
    }

    // Extract frame lines and reconstruct as single string
    let frame_content = lines[frame_start..frame_end].join("\n") + "\n";

    // Delegate to existing parser
    parse_gro_content(&frame_content)
}

fn is_water_residue(name: &str) -> bool {
    matches!(
        name.to_uppercase().as_str(),
        "SOL"
            | "WAT"
            | "HOH"
            | "H2O"
            | "TIP"
            | "TIP3"
            | "TIP3P"
            | "TIP4"
            | "TIP4P"
            | "TIP5"
            | "TIP5P"
            | "SPC"
            | "SPCE"
            | "SPC/E"
            | "HSL"
            | "ICE"
            | "HYD"
            | "HYW"
            | "CAGE"
            | "WCL"
            | "WATC"
    )
}

fn is_oxygen_atom(name: &str) -> bool {
    let upper = name.to_uppercase();
    upper.starts_with('O') || upper.contains("OW")
}

fn is_hydrogen_atom(name: &str) -> bool {
    let upper = name.to_uppercase();
    upper.starts_with('H') && !upper.contains("HE")
}

/// Minimum-image-aware centroid: unwrap every atom relative to the first against the
/// (orthorhombic) box, average, then wrap back in. Matches the frontend
/// `calculateGuestCentroid`, so a guest split across the periodic boundary is located
/// correctly instead of at the box center.
fn mic_centroid(positions: &[[f64; 3]], box_dims: &[f64; 3]) -> [f64; 3] {
    let anchor = positions[0];
    let mut sum = [0.0f64; 3];
    for p in positions {
        for k in 0..3 {
            let mut d = p[k] - anchor[k];
            if box_dims[k] > 0.0 {
                d -= box_dims[k] * (d / box_dims[k]).round();
            }
            sum[k] += anchor[k] + d;
        }
    }
    let n = positions.len() as f64;
    let mut c = [sum[0] / n, sum[1] / n, sum[2] / n];
    for k in 0..3 {
        if box_dims[k] > 0.0 {
            c[k] = c[k].rem_euclid(box_dims[k]);
        }
    }
    c
}

/// Ions, counter-ions, and metal/slab residues that are never clathrate guests.
fn is_non_guest_residue(name: &str) -> bool {
    matches!(
        name.to_uppercase().as_str(),
        // monatomic ions / counter-ions
        "NA" | "NA+" | "CL" | "CL-" | "K" | "K+" | "LI" | "LI+" | "F" | "F-"
            | "BR" | "BR-" | "I" | "I-"
            | "MG" | "MG2+" | "CA" | "CA2+" | "ZN" | "ZN2+" | "SOD" | "CLA" | "POT" | "IB" | "IB+"
        // metals / slabs
            | "FE" | "CR" | "NI" | "CU" | "AL" | "SLAB" | "STEEL" | "METAL"
    )
}

/// A guest is any non-water residue that is not an ion/metal. The caller only reaches
/// this for non-water residues. Deliberately a denylist (not an allowlist) so custom
/// guests are not silently dropped and this agrees with the frontend guest filter.
fn is_guest_residue(name: &str) -> bool {
    !is_non_guest_residue(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_gro() {
        let content = r#"Test system
6
    1SOL     OW    1   0.100   0.200   0.300
    1SOL    HW1    2   0.150   0.200   0.300
    1SOL    HW2    3   0.050   0.200   0.300
    2SOL     OW    4   0.400   0.500   0.600
    2SOL    HW1    5   0.450   0.500   0.600
    2SOL    HW2    6   0.350   0.500   0.600
   1.00000   1.00000   1.00000
"#;

        let result = parse_gro_content(content).unwrap();
        assert_eq!(result.waters.len(), 2);
        assert_eq!(result.box_dims, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_parse_hydrate_water_alias_gro() {
        let content = r#"Hydrate alias system
7
    1HSL     OW    1   0.100   0.200   0.300
    1HSL    HW1    2   0.150   0.200   0.300
    1HSL    HW2    3   0.050   0.200   0.300
    2CH4      C    4   0.400   0.500   0.600
    2CH4     H1    5   0.450   0.500   0.600
    2CH4     H2    6   0.350   0.500   0.600
    3SLAB    Fe    7   0.700   0.800   0.900
   1.00000   1.00000   1.00000
"#;

        let result = parse_gro_content(content).unwrap();
        assert_eq!(result.waters.len(), 1);
        assert_eq!(result.guests.len(), 1);
        assert_eq!(result.guests[0].name, "CH4");
    }

    #[test]
    fn test_parse_cyclopentane_guest_and_excludes_ions() {
        // Cyclopentane (CPN) is a real sII hydrate guest; NA is a counter-ion that must
        // NOT be treated as a guest. The guest position is the centroid of the residue.
        let content = r#"Guest parser test
9
    1SOL     OW    1   0.100   0.200   0.300
    1SOL    HW1    2   0.150   0.200   0.300
    1SOL    HW2    3   0.050   0.200   0.300
    2CPN     C1    4   0.100   0.500   0.500
    2CPN     C2    5   0.200   0.500   0.500
    2CPN     C3    6   0.300   0.500   0.500
    2CPN     C4    7   0.400   0.500   0.500
    2CPN     C5    8   0.500   0.500   0.500
    3NA      NA    9   0.900   0.900   0.900
   1.00000   1.00000   1.00000
"#;

        let result = parse_gro_content(content).unwrap();
        assert_eq!(result.waters.len(), 1);
        assert_eq!(result.guests.len(), 1, "CPN must be a guest; NA must be excluded");
        assert_eq!(result.guests[0].name, "CPN");
        let p = result.guests[0].position;
        assert!((p[0] - 0.3).abs() < 1e-4 && (p[1] - 0.5).abs() < 1e-4 && (p[2] - 0.5).abs() < 1e-4,
            "guest position should be the residue centroid, got {:?}", p);
    }

    #[test]
    fn test_parse_propane_guest_not_dropped_as_biomolecule() {
        // Residue name "PRO" is propane here, not proline; it must be kept as a guest.
        let content = r#"Propane guest test
3
    1PRO     C1    1   0.100   0.500   0.500
    1PRO     C2    2   0.200   0.500   0.500
    1PRO     C3    3   0.300   0.500   0.500
  10.00000  10.00000  10.00000
"#;
        let result = parse_gro_content(content).unwrap();
        assert_eq!(result.guests.len(), 1);
        assert_eq!(result.guests[0].name, "PRO");
    }

    #[test]
    fn test_guest_centroid_is_minimum_image_aware() {
        // A guest split across the periodic boundary must be located near the boundary,
        // not at the box center (which a plain arithmetic mean would give).
        let content = r#"PBC guest test
2
    1GUE     C1    1   9.900   5.000   5.000
    1GUE     C2    2   0.100   5.000   5.000
  10.00000  10.00000  10.00000
"#;
        let result = parse_gro_content(content).unwrap();
        assert_eq!(result.guests.len(), 1);
        let p = result.guests[0].position;
        assert!(p[0] < 0.11 || p[0] > 9.89,
            "MIC centroid should sit near the boundary, not the box center; got {:?}", p);
    }

    #[test]
    fn test_count_frames_single() {
        let content = "Test system\n6\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.150   0.200   0.300\n    1SOL    HW2    3   0.050   0.200   0.300\n    2SOL     OW    4   0.400   0.500   0.600\n    2SOL    HW1    5   0.450   0.500   0.600\n    2SOL    HW2    6   0.350   0.500   0.600\n   1.00000   1.00000   1.00000\n";
        assert_eq!(count_frames(content), 1);
    }

    #[test]
    fn test_count_frames_multi() {
        let content = "Frame 1\n3\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.150   0.200   0.300\n    1SOL    HW2    3   0.050   0.200   0.300\n   1.00000   1.00000   1.00000\nFrame 2\n3\n    1SOL     OW    1   0.200   0.300   0.400\n    1SOL    HW1    2   0.250   0.300   0.400\n    1SOL    HW2    3   0.150   0.300   0.400\n   2.00000   2.00000   2.00000\n";
        assert_eq!(count_frames(content), 2);
    }

    #[test]
    fn test_count_frames_empty() {
        assert_eq!(count_frames(""), 0);
    }

    #[test]
    fn test_parse_gro_frame_first() {
        let content = "Frame 1\n3\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.150   0.200   0.300\n    1SOL    HW2    3   0.050   0.200   0.300\n   1.00000   1.00000   1.00000\nFrame 2\n3\n    1SOL     OW    1   0.200   0.300   0.400\n    1SOL    HW1    2   0.250   0.300   0.400\n    1SOL    HW2    3   0.150   0.300   0.400\n   2.00000   2.00000   2.00000\n";
        let result = parse_gro_frame(content, 0).unwrap();
        assert_eq!(result.box_dims, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_parse_gro_frame_second() {
        let content = "Frame 1\n3\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.150   0.200   0.300\n    1SOL    HW2    3   0.050   0.200   0.300\n   1.00000   1.00000   1.00000\nFrame 2\n3\n    1SOL     OW    1   0.200   0.300   0.400\n    1SOL    HW1    2   0.250   0.300   0.400\n    1SOL    HW2    3   0.150   0.300   0.400\n   2.00000   2.00000   2.00000\n";
        let result = parse_gro_frame(content, 1).unwrap();
        assert_eq!(result.box_dims, [2.0, 2.0, 2.0]);
    }

    #[test]
    fn test_parse_gro_frame_out_of_bounds() {
        let content = "Test\n3\n    1SOL     OW    1   0.100   0.200   0.300\n    1SOL    HW1    2   0.150   0.200   0.300\n    1SOL    HW2    3   0.050   0.200   0.300\n   1.00000   1.00000   1.00000\n";
        assert!(parse_gro_frame(content, 5).is_err());
    }
}
