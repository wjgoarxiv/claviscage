# TRACE

General hydrogen-bond topology search for hydrate cages, including unusual face signatures.
Full page: https://docs.claviscage.app/analysis/trace

## How it works

1. Build a periodic (minimum-image) hydrogen-bond graph from an O–O distance and a donor-angle criterion.
2. Enumerate and geometrically validate 4-, 5- and 6-membered rings.
3. Assemble closed polyhedra from shared ring edges.
4. Classify complete and incomplete candidates by topology and face signature.

Default thresholds: **O–O 0.36 nm**, **donor angle 35°**. These are TRACE settings; HTR and GRADE
use their own contracts.

## Output

Canonical cages such as `5^12`, `5^12·6^2` and `5^12·6^4`, and arbitrary signatures such as
`4^3·5^6·6^3`.

| Quality | Meaning |
|---|---|
| SEC | Complete standard edge-saturated cage (closed polyhedron, `F − E + V = 2`) |
| nSEC | Closed cage with non-standard topology |
| IC | Incomplete cage candidate |

An arbitrary signature is valid TRACE output; it is never presented as GRADE output.

## Reference

J.-W. Hsu and S.-T. Lin, "TRACE: A Topological Algorithm for Detecting Additive-Coordinated
Hydrate Cages", *Journal of Chemical Theory and Computation* 21 (22) (2025) 11436–11454.
https://doi.org/10.1021/acs.jctc.5c01459
