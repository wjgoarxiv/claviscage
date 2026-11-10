// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo

//! Geometry utilities with periodic boundary conditions

/// Periodic simulation cell represented by three Cartesian basis vectors.
///
/// The vectors are stored as `a`, `b`, `c` where Cartesian coordinates are
/// reconstructed from fractional coordinates as `f0*a + f1*b + f2*c`.
/// This supports orthorhombic cells and restricted triclinic cells such as
/// LAMMPS/GROMACS boxes with `a=(lx,0,0)`, `b=(xy,ly,0)`, `c=(xz,yz,lz)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    vectors: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
}

impl Cell {
    /// Build an orthorhombic cell from legacy `[x, y, z]` box lengths.
    pub fn orthorhombic(box_dims: [f64; 3]) -> Result<Self, String> {
        Self::from_vectors([
            [box_dims[0], 0.0, 0.0],
            [0.0, box_dims[1], 0.0],
            [0.0, 0.0, box_dims[2]],
        ])
    }

    /// Build a restricted triclinic cell from lengths and tilt factors.
    pub fn restricted_triclinic(
        lx: f64,
        ly: f64,
        lz: f64,
        xy: f64,
        xz: f64,
        yz: f64,
    ) -> Result<Self, String> {
        Self::from_vectors([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]])
    }

    /// Build a cell from Cartesian basis vectors `a`, `b`, `c`.
    pub fn from_vectors(vectors: [[f64; 3]; 3]) -> Result<Self, String> {
        if !vectors.iter().flatten().all(|v| v.is_finite()) {
            return Err("cell vectors must be finite".to_string());
        }

        let m = [
            [vectors[0][0], vectors[1][0], vectors[2][0]],
            [vectors[0][1], vectors[1][1], vectors[2][1]],
            [vectors[0][2], vectors[1][2], vectors[2][2]],
        ];
        let det = determinant3(&m);
        if det.abs() < 1e-12 {
            return Err("cell vectors must form a non-degenerate cell".to_string());
        }

        Ok(Self {
            vectors,
            inverse: invert3(&m, det),
        })
    }

    /// Parse a browser/WASM cell array.
    ///
    /// Supported forms:
    /// - 3 values: orthorhombic `[lx, ly, lz]`
    /// - 6 values: restricted triclinic `[lx, ly, lz, xy, xz, yz]`
    /// - 9 values: GROMACS `.gro` box order
    ///   `[v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y]`
    pub fn from_box_values(values: &[f64]) -> Result<Self, String> {
        match values.len() {
            3 => Self::orthorhombic([values[0], values[1], values[2]]),
            6 => Self::restricted_triclinic(
                values[0], values[1], values[2], values[3], values[4], values[5],
            ),
            9 => Self::from_vectors([
                [values[0], values[3], values[4]],
                [values[5], values[1], values[6]],
                [values[7], values[8], values[2]],
            ]),
            _ => Err("cell must contain 3, 6, or 9 values".to_string()),
        }
    }

    #[inline]
    pub fn vectors(&self) -> [[f64; 3]; 3] {
        self.vectors
    }

    /// Smallest distance between opposite cell faces. A spherical cutoff below half of it
    /// sees at most one periodic image of any neighbour.
    pub fn min_width(&self) -> f64 {
        let [a, b, c] = self.vectors;
        let volume = dot3(&a, &cross3(&b, &c)).abs();
        [cross3(&b, &c), cross3(&c, &a), cross3(&a, &b)]
            .iter()
            .map(|n| volume / norm3(n))
            .fold(f64::INFINITY, f64::min)
    }

    /// Apply minimum image convention to a Cartesian displacement vector.
    #[inline]
    pub fn mic_delta(&self, delta: [f64; 3]) -> [f64; 3] {
        let mut frac = self.cart_to_frac(delta);
        for v in &mut frac {
            *v -= v.round();
        }
        self.frac_to_cart(frac)
    }

    #[inline]
    pub fn diff(&self, a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
        self.mic_delta([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
    }

    #[inline]
    pub fn distance_sq(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let delta = self.mic_delta([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        dot3(&delta, &delta)
    }

    #[inline]
    pub fn distance(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        self.distance_sq(a, b).sqrt()
    }

    /// Calculate angle `a-center-b` in degrees with cell-aware MIC vectors.
    pub fn angle(&self, center: &[f64; 3], a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let v1 = self.diff(a, center);
        let v2 = self.diff(b, center);

        let mag1 = norm3(&v1);
        let mag2 = norm3(&v2);

        if mag1 < 1e-10 || mag2 < 1e-10 {
            return 180.0;
        }

        let cos_angle = dot3(&v1, &v2) / (mag1 * mag2);
        cos_angle.clamp(-1.0, 1.0).acos() * (180.0 / std::f64::consts::PI)
    }

    /// Calculate signed dihedral angle in radians with cell-aware MIC vectors.
    pub fn dihedral(&self, p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> f64 {
        let b1 = self.diff(p2, p1);
        let b2 = self.diff(p3, p2);
        let b3 = self.diff(p4, p3);

        signed_dihedral_from_bonds(&b1, &b2, &b3)
    }

    /// Calculate a PBC-aware average position and wrap it into the primary cell.
    pub fn average(&self, points: &[[f64; 3]]) -> [f64; 3] {
        if points.is_empty() {
            return [0.0, 0.0, 0.0];
        }
        if points.len() == 1 {
            return self.wrap_point(points[0]);
        }

        let ref_point = points[0];
        let mut sum = [0.0, 0.0, 0.0];
        for p in points {
            let delta = self.diff(p, &ref_point);
            sum[0] += delta[0];
            sum[1] += delta[1];
            sum[2] += delta[2];
        }

        let n = points.len() as f64;
        self.wrap_point([
            ref_point[0] + sum[0] / n,
            ref_point[1] + sum[1] / n,
            ref_point[2] + sum[2] / n,
        ])
    }

    #[inline]
    fn cart_to_frac(&self, cart: [f64; 3]) -> [f64; 3] {
        mat_vec_mul(&self.inverse, cart)
    }

    #[inline]
    fn frac_to_cart(&self, frac: [f64; 3]) -> [f64; 3] {
        [
            frac[0] * self.vectors[0][0]
                + frac[1] * self.vectors[1][0]
                + frac[2] * self.vectors[2][0],
            frac[0] * self.vectors[0][1]
                + frac[1] * self.vectors[1][1]
                + frac[2] * self.vectors[2][1],
            frac[0] * self.vectors[0][2]
                + frac[1] * self.vectors[1][2]
                + frac[2] * self.vectors[2][2],
        ]
    }

    #[inline]
    fn wrap_point(&self, point: [f64; 3]) -> [f64; 3] {
        let mut frac = self.cart_to_frac(point);
        for v in &mut frac {
            *v -= v.floor();
            if *v >= 1.0 {
                *v -= 1.0;
            }
        }
        self.frac_to_cart(frac)
    }
}

#[inline]
fn determinant3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[inline]
fn invert3(m: &[[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

#[inline]
fn mat_vec_mul(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Apply minimum image convention to a displacement vector
#[inline]
pub fn apply_mic(delta: [f64; 3], box_dims: &[f64; 3]) -> [f64; 3] {
    Cell::orthorhombic(*box_dims)
        .expect("orthorhombic box must be non-degenerate")
        .mic_delta(delta)
}

/// Calculate distance squared between two points with PBC
#[inline]
pub fn pbc_distance_sq(a: &[f64; 3], b: &[f64; 3], box_dims: &[f64; 3]) -> f64 {
    Cell::orthorhombic(*box_dims)
        .expect("orthorhombic box must be non-degenerate")
        .distance_sq(a, b)
}

/// Calculate distance between two points with PBC
#[inline]
pub fn pbc_distance(a: &[f64; 3], b: &[f64; 3], box_dims: &[f64; 3]) -> f64 {
    pbc_distance_sq(a, b, box_dims).sqrt()
}

/// Calculate vector norm
#[inline]
pub fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Dot product of two 3D vectors
#[inline]
pub fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3D vectors
#[inline]
pub fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Calculate angle between vectors a-center-b in degrees with PBC
pub fn pbc_angle(center: &[f64; 3], a: &[f64; 3], b: &[f64; 3], box_dims: &[f64; 3]) -> f64 {
    Cell::orthorhombic(*box_dims)
        .expect("orthorhombic box must be non-degenerate")
        .angle(center, a, b)
}

/// Calculate dihedral angle between four points with PBC
/// Returns angle in radians
pub fn calculate_dihedral(
    p1: &[f64; 3],
    p2: &[f64; 3],
    p3: &[f64; 3],
    p4: &[f64; 3],
    box_dims: &[f64; 3],
) -> f64 {
    Cell::orthorhombic(*box_dims)
        .expect("orthorhombic box must be non-degenerate")
        .dihedral(p1, p2, p3, p4)
}

fn signed_dihedral_from_bonds(b1: &[f64; 3], b2: &[f64; 3], b3: &[f64; 3]) -> f64 {
    let n1 = cross3(&b1, &b2);
    let n2 = cross3(&b2, &b3);

    let n1_norm = norm3(&n1);
    let n2_norm = norm3(&n2);
    let b2_norm = norm3(&b2);

    if n1_norm < 1e-10 || n2_norm < 1e-10 || b2_norm < 1e-10 {
        // Degenerate case: collinear atoms. Return π/2 so cos(3φ) ≈ 0 (neutral contribution)
        return std::f64::consts::FRAC_PI_2;
    }

    // Use atan2 for SIGNED dihedral angle (-π to +π)
    // Reference: F4SPC (https://github.com/adibifard/F4SPC)
    // Formula: x = dot(b2, cross(n1, n2)), y = norm(b2) * dot(n1, n2)
    // theta = atan2(x, y)
    //
    // This is critical for F4 order parameter:
    // - Hydrates have eclipsed geometry (φ ≈ 0°) → cos(3φ) ≈ +1
    // - Ice has staggered geometry (φ ≈ 60°) → cos(3φ) ≈ -1
    // - Liquid has random angles → cos(3φ) ≈ 0
    let n1_cross_n2 = cross3(&n1, &n2);
    let x = dot3(&b2, &n1_cross_n2);
    let y = b2_norm * dot3(&n1, &n2);

    x.atan2(y)
}

/// Compute MIC displacement vector from point b to point a: a - b with PBC wrapping
#[inline]
pub fn apply_mic_diff(a: &[f64; 3], b: &[f64; 3], box_dims: &[f64; 3]) -> [f64; 3] {
    apply_mic([a[0] - b[0], a[1] - b[1], a[2] - b[2]], box_dims)
}

/// Calculate average position with PBC
pub fn pbc_average(points: &[[f64; 3]], box_dims: &[f64; 3]) -> [f64; 3] {
    if points.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    if points.len() == 1 {
        return points[0];
    }

    Cell::orthorhombic(*box_dims)
        .expect("orthorhombic box must be non-degenerate")
        .average(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-9,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn restricted_triclinic_mic_uses_fractional_wrapping() {
        // Restricted triclinic cell: a=(2,0,0), b=(0.5,2,0), c=(0.25,0.4,2).
        let cell = Cell::restricted_triclinic(2.0, 2.0, 2.0, 0.5, 0.25, 0.4).unwrap();
        let delta = [1.30, 0.0, 0.0];

        // Orthorhombic MIC would be -0.70 along x. In the triclinic basis the nearest
        // image wraps by one a-vector, giving the same x displacement here and proving
        // the new path operates through Cell rather than box length division.
        let mic = cell.mic_delta(delta);
        assert_close(mic[0], -0.70);
        assert_close(mic[1], 0.0);
        assert_close(mic[2], 0.0);
    }

    #[test]
    fn restricted_triclinic_distance_can_cross_tilted_boundary() {
        let cell = Cell::restricted_triclinic(2.0, 2.0, 2.0, 0.5, 0.0, 0.0).unwrap();
        let a = [0.10, 0.10, 0.10];
        let b = [0.55, 1.90, 0.10];

        // Raw delta is [0.45, 1.8, 0]. Wrapping by one b-vector subtracts [0.5, 2, 0],
        // producing [-0.05, -0.2, 0] with distance sqrt(0.0425).
        assert_close(
            cell.distance(&a, &b),
            (0.05_f64 * 0.05 + 0.2_f64 * 0.2).sqrt(),
        );
    }

    #[test]
    fn cell_average_wraps_cluster_across_triclinic_boundary() {
        let cell = Cell::restricted_triclinic(2.0, 2.0, 2.0, 0.5, 0.0, 0.0).unwrap();
        let points = [[0.10, 0.10, 0.0], [0.55, 1.90, 0.0]];

        let avg = cell.average(&points);
        assert_close(avg[0], 0.075);
        assert_close(avg[1], 0.0);
        assert_close(avg[2], 0.0);
    }

    #[test]
    fn cell_angle_and_dihedral_match_orthorhombic_legacy_helpers() {
        let cell = Cell::orthorhombic([10.0, 10.0, 10.0]).unwrap();
        let center = [0.0, 0.0, 0.0];
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        assert_close(cell.angle(&center, &x, &y), 90.0);

        let p1 = [0.0, 1.0, 0.0];
        let p2 = [0.0, 0.0, 0.0];
        let p3 = [1.0, 0.0, 0.0];
        let p4 = [1.0, 1.0, 0.0];
        assert_close(cell.dihedral(&p1, &p2, &p3, &p4), 0.0);
    }
}
