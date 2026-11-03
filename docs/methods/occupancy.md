# Guest occupancy

How many detected cages hold a guest. Full page: https://docs.claviscage.app/analysis/occupancy

- **θ_S**: occupancy of small cages (`5^12`).
- **θ_L**: occupancy of large cages such as `5^12·6^2` and `5^12·6^4`.
- **Overall**: occupied cages divided by all detected cages in the current result.

## How it works

After TRACE, GRADE or complete-cage HTR analysis, the app measures the periodic minimum-image
distance from each candidate guest to each cage centre. A cage-type-specific radius admits nearby
candidates, and a one-guest-per-cage rule keeps the closest. This is a centre-and-radius model, not
a point-in-polyhedron test. HTR incomplete cages are not occupiable. Cage centres are computed in
the periodic cell, so cages that cross the box boundary get the right centre.

Occupied cages stay attached to the algorithm that found them; switching algorithms restores that
algorithm's own result.

## Reference

No paper record for the occupancy definitions in our sources yet.
