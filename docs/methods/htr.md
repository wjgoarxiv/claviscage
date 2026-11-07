# HTR

Primitive-ring and generalized cage analysis on the oxygen topology. HTR is a separate algorithm
from TRACE; equal counts on an ideal crystal do not make the two equivalent.
Full page: https://docs.claviscage.app/analysis/htr

## Topology contract

1. Oxygen-only proximity graph with minimum-image O–O distances.
2. Connect oxygens at or below a fixed **0.35 nm** cutoff.
3. Keep at most six neighbours per oxygen.
4. Enumerate primitive 4-, 5- and 6-membered rings.
5. Assemble cages and report face counts as `n4-n5-n6`.

HTR does not use hydrogen donor orientation, so the TRACE angle setting never changes HTR topology.

## Output

`0-12-0` is a `5^12` cage and `0-12-2` is a `5^12·6^2` cage; HTR is not limited to canonical
labels. Ring and cage overlays come from the same result as the counts. Incomplete cages are shown
but are not occupiable; occupancy uses complete HTR cages only.

## Reference

Y. Liu, K. Xu, Y. Xu, J. Liu, J. Wu and Z. Zhang, "HTR: An ultra-high speed algorithm for cage
recognition of clathrate hydrates", *Nanotechnology Reviews* 11 (1) (2022) 699–711.
https://doi.org/10.1515/ntrev-2022-0044

HTR+ extends the method to hydrate types and polycrystals: Q. Shi, Z. Lin, Y. Qu, J. Wu and Z. Zhang,
"HTR+: a novel algorithm for identifying type and polycrystal of gas hydrates", *Journal of Physics:
Condensed Matter* 36 (2024) 365901. https://doi.org/10.1088/1361-648X/ad52df
