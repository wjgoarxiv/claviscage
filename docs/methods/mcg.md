# MCG-3 (mutually coordinated guest)

Size of the largest cluster of guests that are mutually coordinated through shared water, after
Barnes, Beckham, Wu and Sum. It follows hydrate nucleation and growth from the guests' side; it is
not a cage detector.

Method id: `mcg`. Verification level: **T1** (implemented from the paper and checked on cases with
known answers).

## What it computes

1. For a candidate guest g, every other guest h closer than the guest-guest distance is considered.
2. The pair (g, h) qualifies when at least `n_water` water oxygens lie within the guest-water
   distance of both guests and inside both cones of half-angle `phi_deg` projected along the
   guest-guest vector.
3. g is an MCG monomer when it has at least `n_count` qualifying neighbours.
4. Monomers joined by qualifying pairs form clusters; the result is the number of guests in the
   largest cluster (`largest_cluster`).

| Parameter | Default (MCG-3) | Range |
|---|---|---|
| Guest-guest distance (`r_guest_nm`) | 0.9 nm | 0.5 to 1.5 nm |
| Guest-water distance (`r_water_nm`) | 0.6 nm | 0.3 to 1.0 nm |
| Cone angle (`phi_deg`) | 45 degrees | 15 to 90 degrees |
| Shared waters (`n_water`) | 5 | 1 to 10 |
| MCG-n (`n_count`) | 3 | 1 to 6 |

`n_count = 1` gives MCG-1. A water is located at its oxygen. The method needs guests; on a frame
without guests it reports a missing input instead of zero.

## Verification

Hand-built pairs, and fully occupied CH4 sI and sII crystals (`tests/fixtures/guests/`), where one
cluster holds all 64 and 192 guests. There is no program of the authors to compare with, so this is
level T1.

## Reference

Barnes, Beckham, Wu and Sum, *J. Chem. Phys.* 140 (2014) 164506.
https://doi.org/10.1063/1.4871898
