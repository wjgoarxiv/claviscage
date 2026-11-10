# Guest-containing crystal fixtures

Fully occupied CH4 hydrates used by the MCG and Steinhardt tests.

| File | Structure | Waters | CH4 | sha256 |
| --- | --- | ---: | ---: | --- |
| `sI_CH4_222.gro` | sI (CS1), 2×2×2 | 368 | 64 | `30c462a13b5a3c94a527324994b08b0923649967bc4cfa2a33a4957351ccebcd` |
| `sII_CH4_222.gro` | sII (CS2), 2×2×2 | 1088 | 192 | `380049b3be548c01608a988aef2f773e4dee0d01159a0788baf2f31ad0aad7dc` |

Generated with GenIce2 2.2.13.3 (MIT, https://github.com/genice-dev/GenIce2):

    genice2 CS1 -r 2 2 2 -g 12=me -g 14=me -s 1 -f gromacs
    genice2 CS2 -r 2 2 2 -g 12=me -g 16=me -s 1 -f gromacs

The united-atom guest residue `Me`/atom `Me` was renamed to residue `CH4`, atom `C`;
coordinates are unchanged. The Steinhardt reference JSON next to the other references
is produced by `scripts/make_detector_fixtures.py` (freud, test data only).
