# Benchmark crystals

Ideal crystals written by [GenIce2](https://github.com/vitroid/GenIce) (MIT licence), version **2.2.13.1**, installed with `uv pip install genice2==2.2.13.1` in a temporary virtual environment that was deleted after use. Every file uses seed 1000 and TIP3P water geometry. GenIce2 places the oxygens on the ideal lattice and picks one hydrogen arrangement that obeys the ice rules. The seed fixes that arrangement, so rerunning these commands writes the same files.

| File | Command | Waters | Box |
|------|---------|-------:|-----|
| `sI_222.gro` | `genice2 CS1 --rep 2 2 2 --seed 1000 --format gromacs --water tip3p` | 368 | cubic |
| `sII_222.gro` | `genice2 CS2 --rep 2 2 2 --seed 1000 --format gromacs --water tip3p` | 1088 | cubic |
| `sH_222.gro` | `genice2 sH --rep 2 2 2 --seed 1000 --format gromacs --water tip3p` | 544 | hexagonal (triclinic, 9-value box line) |
| `Ih_332.gro` | `genice2 1h --rep 3 3 2 --seed 1000 --format gromacs --water tip3p` | 288 | orthorhombic |
| `Ic_222.gro` | `genice2 1c --rep 2 2 2 --seed 1000 --format gromacs --water tip3p` | 64 | cubic |

## Answers fixed by construction

These counts follow from the crystal structure itself, not from any program:

- **sI**: each unit cell holds 2 small (5¹²) and 6 large (5¹²6²) cages. 8 cells give 16 + 48.
- **sII**: each unit cell holds 16 small (5¹²) and 8 large (5¹²6⁴) cages. 8 cells give 128 + 64.
- **sH**: each hexagonal cell of 34 waters holds 3 small (5¹²), 2 medium (4³5⁶6³) and 1 large (5¹²6⁸) cage. The 544 waters make 16 such cells, giving 48 + 32 + 16. TRACE and GRADE only type 5¹², 5¹²6² and 5¹²6⁴ cages, so they are held to 48 small cages and none of the other two types.
- **Ih / Ic**: every water of hexagonal (cubic) ice is hexagonal (cubic) ice under CHILL+.

## Results table

`cargo test --release --test benchmark_table` runs every registered method on its own descriptor fixtures and on these crystals. It writes `benchmarks/results.json`, and two runs write the same bytes. Rows with an answer above must match exactly. Any other row is recorded as it comes out, so a later change to that number shows up in the diff. MCG rows on these crystals record `missing_input: guests`, because they contain no guests.

Two points from building the table:

- **GRADE** misses its exact answer on the three GenIce hydrates. The GRADE v1.00 program tests ring planarity on raw coordinates, with no minimum image. A ring that crosses the box edge therefore fails, and the count depends on where the edge falls. Shifting the same sI crystal and wrapping it back into the box gives between 2 and 9 small cages instead of 16. The port reproduces the program on purpose (its T3 level rests on equal counts with GRADE v1.00), so these rows are listed as known mismatches, not loosened.
- The table reads the full box line of each file, which matters for the triclinic sH cell. Since 0.3.17 the `.gro` parser also keeps it (`box_vectors`); before that fix it kept only the diagonal, and TRACE found 41 small cages in sH instead of 48.
