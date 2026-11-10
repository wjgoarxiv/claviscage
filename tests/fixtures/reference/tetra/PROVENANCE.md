# Tetrahedral order reference values

Output of GROMACS **2026.2** (mixed precision, Homebrew build) `gmx order -Sg -Sk` on four fixtures. This folder holds only the outputs and the command log (`commands.txt`); no GROMACS code is in this repository.

| Fixture | Waters | S_g (`*_sg.xvg`) | S_k (`*_sk.xvg`) |
|---------|-------:|-----:|-----:|
| `crystals/example_H2O_frame0.gro` (liquid/interface slab) | 3000 | 0.057853 | 0.001457 |
| `crystals/222_S2.gro` (sII) | 1088 | 0.002704 | 0.000028 |
| `benchmark/Ih_332.gro` (ice Ih) | 288 | 0.000490 | 0.000007 |
| `benchmark/sH_222.gro` (sH, triclinic) | 544 | 0.007654 | 0.000007 |

How each value was made: `gmx order` needs a run input (`-s`), so `gmx grompp` built one from the same coordinates. The topology was TIP3P water from `amber99sb-ildn.ff`, and the `.mdp` had zero steps and 0.6 nm cut-offs. The fixture itself was read as the trajectory (`-f`). An index group listed the oxygens, 15 numbers per line. A single long line is cut at 4096 characters by GROMACS and corrupts the group. That cut is why a first attempt gave `nan` for the liquid frame. The run input was built in a temporary folder that was deleted afterwards. The `.mdp` and topology were not kept; they do not affect these numbers, which `gmx order -Sg -Sk` computes from the group's coordinates and the box alone. In each `.xvg` header the `Working dir:` line, which named that local temporary folder, was replaced by a placeholder; nothing else was edited.

What GROMACS computes (`src/gromacs/gmxana/gmx_order.cpp`, read for the definitions only):

- the four nearest atoms of the group under the minimum image, with no cutoff;
- S_g = 3/32 · Σ_{j<k} (cos ψ_jk + 1/3)², Chau & Hardwick eq. 3;
- S_k = Σ_k (r̄ − r_k)² / (12 r̄²);
- the frame value is the mean over the group, printed with six decimals.

The engine's `tetra` reports these and q = 1 − 4 S_g (Errington & Debenedetti). `src/tetra.rs` and `pytests/test_new_methods_parity.py` check the frame means within 2·10⁻⁶.
