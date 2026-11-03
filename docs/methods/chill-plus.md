# CHILL+

Per-water classification of the local orientational environment. It complements cage finders:
TRACE reports closed cage topology, CHILL+ shows where ordered ice-like and hydrate-like water is.
Full page: https://docs.claviscage.app/analysis/chill-plus

## Classes

Other / liquid, hexagonal ice, cubic ice, interfacial ice, gas hydrate, interfacial hydrate.
The result reports the count and fraction of waters in each class.

Only recognized water residues are classified; guests, ions and other residues stay ordinary scene
geometry. The CHILL+ network drawn in the scene is the classifier's O–O neighbour lattice; the
H-bond overlay is separate and uses its own cutoff (default 0.35 nm) and angle.

## Interpreting

Treat class fractions as a structural descriptor, compare them over time, and use a cage-topology
method (TRACE, HTR, GRADE) when you need explicit cage counts.

## Reference

A. H. Nguyen and V. Molinero, "Identification of Clathrate Hydrates, Hexagonal Ice, Cubic Ice, and
Liquid Water in Simulations: the CHILL+ Algorithm", *J. Phys. Chem. B* 119 (29) (2015) 9369–9376.
https://doi.org/10.1021/jp510289t
