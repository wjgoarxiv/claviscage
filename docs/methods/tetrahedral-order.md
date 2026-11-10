# Tetrahedral order (q, S_g, S_k)

Tetrahedral order of each water from its four nearest oxygens: the orientational order q of
Errington and Debenedetti, and the angular and distance parameters S_g and S_k of Chau and
Hardwick. It is not a cage detector.

Method id: `tetra`. Verification level: **T2** (matches an established open implementation).

## What it computes

For each water, the four nearest oxygens under the minimum image (no cutoff; ties keep the lower
index) give S_g, S_k and q = 1 - 4 S_g. Results are per-water values (`q`, `sg`, `sk`) and their
frame means. The method has no parameters and needs only oxygen positions.

## Choices where the papers leave a detail open

- Neighbours: the papers say the four nearest neighbours; the engine takes the four nearest oxygens
  under the minimum image, with no cutoff, as `gmx order` does.
- S_k normalisation: the mean distance is that of the same four neighbours, giving the sum of
  squared deviations over 12 r_mean^2, as `gmx order` does.

## Limit

The neighbour search compares every pair of waters, so a frame is limited to 20,000 waters for now;
a larger frame returns an error that says so.

## Verification

Frame means of S_g and S_k match GROMACS 2026.2 `gmx order -Sg -Sk` within 2e-6 on a liquid frame,
sII, hexagonal ice and the triclinic sH crystal (`tests/fixtures/reference/tetra/`). GROMACS
(LGPL-2.1-or-later) was used to generate the test values only; no GROMACS code is included. See
`NOTICE`.

## References

- Errington and Debenedetti, *Nature* 409 (2001) 318. https://doi.org/10.1038/35053024
- Chau and Hardwick, *Mol. Phys.* 93 (1998) 511. https://doi.org/10.1080/002689798169195
