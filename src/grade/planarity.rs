// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Planarity test of `coplanar_Points` in the GRADE reference.

/// The reference converts radians with this truncated constant; the exact
/// value shifts an angle by about 5e-8 degrees at the 45 degree cutoff.
#[allow(clippy::approx_constant)]
const PI_REFERENCE: f64 = 3.14159265;

type Vec3 = [f64; 3];

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn magnitude(a: Vec3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Angle between two normals in degrees, with no clamping: a cosine that
/// rounds just above one gives NaN, which fails every `< theta` test.
fn angle(a: Vec3, b: Vec3) -> f64 {
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    (dot / (magnitude(a) * magnitude(b))).acos() * 180.0 / PI_REFERENCE
}

/// True when the open vertex list (five or six oxygens) is planar.
///
/// The reference measures the angle between the normals of the triangle fan
/// rooted at the first vertex, on raw coordinates: no minimum-image shift is
/// applied, so a ring that straddles a periodic boundary is judged on its
/// wrapped positions and normally fails.
pub(super) fn is_planar(ring: &[usize], positions: &[[f64; 3]], theta_deg: f64) -> bool {
    let a = positions[ring[0]];
    let spokes: Vec<Vec3> = ring[1..].iter().map(|&v| sub(positions[v], a)).collect();
    let normals: Vec<Vec3> = spokes.windows(2).map(|w| cross(w[0], w[1])).collect();
    normals
        .windows(2)
        .all(|pair| angle(pair[0], pair[1]) < theta_deg)
}
