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
| `.claviscage.zip` bundle | Schema v1 and v2 round-trip in the web app and Blender extension regression tests | Those tests and their fixtures are not in this repository |

Reference tools named above (TRACE, HTR, GRADE, CalF4) belong to their authors and are cited in
`CITATION.cff`; they are not included in this repository.
