# Averaged Steinhardt q4 / q6

Locally averaged bond-orientational order of each water oxygen, after Lechner and Dellago. It
separates crystalline from liquid-like waters; it is not a cage detector.

Method id: `steinhardt`. Verification level: **T2** (matches an open reference implementation).

## What it computes

For every oxygen, the neighbours inside the cutoff define the bond-orientational vector for
l = 4 and l = 6; the vectors are then averaged over the oxygen and its neighbours, giving the
per-water values `qbar4` and `qbar6`. Distances use the minimum image, in orthorhombic and
triclinic cells.

| Parameter | Default | Range |
|---|---|---|
| Neighbour cutoff (`cutoff_nm`) | 0.35 nm | 0.30 to 0.40 nm |

Only oxygen positions are needed. The Python bindings (analysis kind `steinhardt`) take the published
defaults (l = 4 and 6); the WebAssembly surface also accepts overrides.

## Verification

Per-water values match freud's `Steinhardt(average=True)` within its float32 precision on the sI and
sII crystals, a CH4-filled sI crystal and a disordered 4,500-water frame
(`tests/fixtures/reference/steinhardt_*.json`). freud (BSD-3-Clause) was used to generate the test
values only; no freud code is included. See `NOTICE`.

## Reference

Lechner and Dellago, *J. Chem. Phys.* 129 (2008) 114707.
https://doi.org/10.1063/1.2977970
