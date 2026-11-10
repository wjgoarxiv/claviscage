# Validation summary

Numbers the ClavisCage engine reproduces, as recorded in the project changelog and docs. The inputs
and expected outputs behind the engine tests are in `tests/fixtures/`; other inputs are not
redistributed here, and blank cells are things we have not published.

| Check | Result | Not claimed |
|---|---|---|
| TRACE, bundled `example_H2O` reference | 1,366 cups; SEC 10, nSEC 0, IC 12 (exact match) | Agreement on other inputs without their own reference comparison |
| HTR, 1,242-water oracle | 216 cages = 54 `0-12-0` + 162 `0-12-2` (exact match) | An exact full-trajectory reproduction: the original 1,001-frame input is unavailable |
| HTR, production UI check on the same 1,242-water file | 216 cages, 1,458 rings | Per-size ring split |
| TRACE, same UI check | Same 54 / 162 canonical split | SEC / nSEC / IC for that file |
| GRADE, same UI check | 35 `5^12`, 96 `5^12·6^2` | `5^12·6^4` total and the three filled counts |
| TRACE and HTR on perfect GenIce2 crystals | sII 2×2×2: 128 `5^12` + 64 `5^12·6^4` = 192; sI: 128 `5^12` + 384 `5^12·6^2` = 512; identical over repeated runs | |
| F3 / F4, supplied `prod.gro` | F3 0.03089055687939763, F4 0.5403493788565187; match an independent calculation and CalF4 | Water model, cutoff and frame are not restated in this record |
| HTR occupancy, 4,500-water `prod.gro` | 164 complete cages with finite periodic centres, 161 occupied, θ_S 98.18 %, θ_L 100 % | |
| CHILL+, GenIce2 ice Ih, ice Ic, sI, sII | 100 % expected class | Per-class counts |
| Averaged q4 / q6 (engine 0.3.18) | Per-water values match freud `Steinhardt(average=True)` within its float32 precision on sI, sII, CH4-filled sI and a disordered 4,500-water frame (level T2) | Agreement with other implementations |
| Tetrahedral order (engine 0.3.18) | Frame means of S_g and S_k match GROMACS 2026.2 `gmx order -Sg -Sk` within 2e-6 on a liquid frame, sII, hexagonal ice and the triclinic sH crystal (level T2) | Frames of more than 20,000 waters |
| MCG-3 (engine 0.3.18) | One cluster of 64 guests in CH4-filled sI and 192 in CH4-filled sII (level T1) | Agreement with a program of the authors; none is available |
| TRACE and HTR on a GenIce2 sH 2×2×2 crystal | 48 `5^12`, 32 `4^3 5^6 6^3`, 16 `5^12 6^8` (engine test with the full triclinic cell) | |
| `.claviscage.zip` bundle | Schema v1 and v2 round-trip in the web app and Blender extension regression tests | Those tests and their fixtures are not in this repository |

Reference tools named above (TRACE, HTR, GRADE, CalF4) belong to their authors and are cited in
`CITATION.cff`; they are not included in this repository.

## Benchmark on ideal crystals (engine 0.3.18)

`benchmarks/results.json` is written by `cargo test --release --test benchmark_table`. It runs every
registered method on its own fixtures and on five GenIce2 2.2.13.1 crystals (seed 1000, TIP3P water;
commands in `tests/fixtures/benchmark/PROVENANCE.md`): sI 2×2×2 (368 waters; 16 `5^12` + 48
`5^12·6^2`), sII 2×2×2 (1,088; 128 + 64), sH 2×2×2 (544; 48 `5^12`, 32 `4^3 5^6 6^3`, 16 `5^12 6^8`),
hexagonal ice Ih (288) and cubic ice Ic (64). Rows whose answer is known by construction must match
exactly; every other row is recorded as it comes out, so a later change shows up in the diff.

| Method | Exact on | Notes |
|---|---|---|
| TRACE | sI, sII, and the `5^12` count of sH | The table records TRACE's `5^12`, `5^12·6^2` and `5^12·6^4` counts; no cages in either ice |
| HTR | sI, sII, sH including `4^3 5^6 6^3` and `5^12 6^8` | No cages in either ice |
| GRADE | none | Known mismatch, below |
| CHILL+ | every water classed hydrate in sI and sII, Ih in Ih, Ic in Ic | sH is recorded, not asserted |
| F3/F4, averaged q4/q6, tetrahedral order | recorded values | No count to be exact against |
| MCG | 64 guests in CH4-filled sI, 192 in CH4-filled sII | Reports a missing input on the water-only crystals |

**Known mismatch, GRADE on periodic crystals.** GRADE v1.00 tests ring planarity on raw coordinates
without the minimum image, so a ring that crosses the box edge is dropped and the count depends on
where the edge falls: shifting and re-wrapping the same sI crystal gives between 2 and 9 of its 16
small cages. The engine's GRADE matches the authors' program, so it keeps this behaviour and the rows
are listed as known mismatches. Use TRACE or HTR to count cages in a periodic crystal.
