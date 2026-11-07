# GRADE

A cage detector with a deliberately narrow output: three cage families only.
Full page: https://docs.claviscage.app/analysis/grade

| Machine label | Face signature |
|---|---|
| `512` | `5^12` |
| `51262` | `5^12·6^2` |
| `51264` | `5^12·6^4` |

GRADE builds an oxygen-distance graph, finds planar 5- and 6-membered rings, assembles cups and
pairs compatible cups into those families. It does not run TRACE, does not return arbitrary
signatures, and gives no sI/sII/sH verdict or confidence score.

| Parameter | Default |
|---|---|
| O–O cutoff | 0.35 nm |
| Ring-planarity threshold | 45° |
| Pentagon diagonal tolerance (`d1`) | 0.18 nm |
| Hexagon diagonal tolerance (`d2`) | 0.26 nm |

With guest sites, GRADE also reports `filled_512`, `filled_51262` and `filled_51264`: filled-cage
counts, not a structure metric. The H-bond cutoff control does not override GRADE.

The engine's GRADE is a step-by-step port of GRADE v1.00 and returns the same cage counts as the
reference program on the test crystals and production frame in `tests/fixtures/reference/`. GRADE
v1.00 does not return the ideal counts of perfect periodic crystals (for example 46 + 24 instead of
128 + 64 on sII 2×2×2), so use TRACE or HTR when exact crystal counts matter.

## Reference

F. Mahmoudinobar and C. L. Dias, "GRADE: A code to determine clathrate hydrate structures",
*Computer Physics Communications* 244 (2019) 385–391. https://doi.org/10.1016/j.cpc.2019.06.004
(The GRADE program is distributed by its authors; it is not included here.)
