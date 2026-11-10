# Changelog

User-visible changes to the ClavisCage web app and to this documentation. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The analysis engine is in this
repository; the web app itself is not (see `NOTICE.md`).

## Unreleased — analysis engine

### Added
- The ClavisCage analysis engine (Rust crate, v0.3.17) under GPL-3.0-or-later: `src/`, reference
  tests and fixtures in `tests/`, Python binding tests in `pytests/`, `LICENSE.txt` and `NOTICE`.
- README section on building, testing and running the command-line tool.

### Fixed
- Triclinic `.gro` boxes keep their tilt: the parser returns all nine box values (`box_vectors`),
  and the Python binding and `clavis-cli` use the full cell. TRACE and HTR now find 48/32/16
  cages on a GenIce2 sH 2×2×2 crystal (fixture `tests/fixtures/crystals/222_SH.gro`).

### Changed
- The documentation licence file is now `LICENSE-docs.md` (CC BY 4.0, unchanged terms).
- GRADE method note: the engine's GRADE is a step-by-step port of GRADE v1.00 and matches the
  reference program's counts.

## 2026-11-03 — public documentation repository

### Added
- This repository: method notes, a validation summary, the `.claviscage.zip` bundle schema,
  citation metadata, issue templates, a cover loop and app screenshots.

### Current app (https://claviscage.app)
- Free, no account, runs in the browser. TRACE, HTR, GRADE, CHILL+, F3/F4 and guest occupancy on
  GRO, PDB/mmCIF, GROMACS (GRO + XTC/TRR/TPR) and LAMMPS (data + dump) inputs.
- Light and dark themes in a computational-brutalist style; the 3D scene follows the theme.
- Path-traced PNG export renders the number of samples you choose.
- MP4 export and the 3D scene work offline: the ffmpeg.wasm core and the scene's environment map
  are served by the site itself, with no CDN.
- Large `.gro` trajectories (100 MB and up) stream frame by frame.
