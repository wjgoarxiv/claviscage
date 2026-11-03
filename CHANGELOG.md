# Changelog

User-visible changes to the ClavisCage web app and to this documentation. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The application itself is not in this
repository; see `NOTICE.md`.

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
