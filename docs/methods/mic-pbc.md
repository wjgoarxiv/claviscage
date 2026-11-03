# Periodic boundaries (minimum image)

Every distance the app measures across a periodic cell uses the minimum-image convention:

- O–O distances for H-bonds, rings and cages (TRACE, HTR, GRADE)
- guest-to-cage-centre distances for occupancy
- H–H separations and O···O neighbours for F3/F4
- bonds and overlay segments that cross the box edge, which are drawn unwrapped instead of
  stretching across the box

Orthorhombic and triclinic cells follow the same rule. We do not cite a specific paper for the
convention itself.
