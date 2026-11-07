// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
// Ported from GRADE v1.00, Copyright (c) 2018 Farbod Mahmoudinobar and Cristiano L. Dias, GPL-3.0-or-later; modified.

//! Distance rule used by the GRADE reference.

use crate::geometry::Cell;

/// Periodic distance with the reference's own arithmetic.
///
/// `calc_Distance` in the reference `MyFunctions.cpp` folds each component with
/// `if (|d| >= L / 2) d = L - |d|` and reads only the three box lengths.
/// For an orthorhombic cell the same expression is evaluated here so that
/// cutoff comparisons agree bit for bit. Tilted cells have no reference
/// behaviour and fall back to the general minimum-image distance.
#[derive(Debug, Clone, Copy)]
pub(super) struct Metric {
    cell: Cell,
    lengths: Option<[f64; 3]>,
}

impl Metric {
    pub(super) fn new(cell: &Cell) -> Self {
        let v = cell.vectors();
        let diagonal = (0..3).all(|i| (0..3).all(|j| i == j || v[i][j] == 0.0));
        Self {
            cell: *cell,
            lengths: diagonal.then(|| [v[0][0], v[1][1], v[2][2]]),
        }
    }

    pub(super) fn distance(&self, a: &[f64; 3], b: &[f64; 3]) -> f64 {
        let Some(lengths) = self.lengths else {
            return self.cell.distance(a, b);
        };
        let fold = |k: usize| {
            let d = b[k] - a[k];
            if d.abs() >= lengths[k] * 0.5 {
                lengths[k] - d.abs()
            } else {
                d
            }
        };
        let (dx, dy, dz) = (fold(0), fold(1), fold(2));
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
