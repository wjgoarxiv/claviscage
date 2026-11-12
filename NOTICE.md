# Notice

This repository holds the ClavisCage analysis engine and its documentation. ClavisCage is a tool for
analysing gas hydrate cages; the browser application runs at https://claviscage.app.

- **Engine source** (`src/`, `tests/`, `pytests/`, `scripts/`, `Cargo.toml`, `pyproject.toml`):
  GPL-3.0-or-later, see `LICENSE.txt`. The engine includes ports of TRACE and GRADE (both
  GPL-3.0-or-later) and implements HTR from its paper, checked against the authors' program, whose
  code is not included; `NOTICE` lists the original authors and
  the origin of the test fixtures.
- **Documentation, method notes, diagrams and the bundle schema**: CC BY 4.0, see
  `LICENSE-docs.md`.
- **The ClavisCage web application and the Blender extension** are separate programs built on the
  engine. Their source is not in this repository, and nothing here grants a licence to them.

Papers and third-party tools mentioned here are cited, not redistributed; their authors' own
licences apply to them.
