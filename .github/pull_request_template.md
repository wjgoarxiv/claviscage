## Engine or docs change

This repository holds the ClavisCage analysis engine (`src/`, `tests/`) and its documentation; pull
requests for either are welcome. Pull requests that add web app or Blender extension source,
WebAssembly or other binaries, Python wheels, reference-tool code or simulation trajectories will
be closed.

- [ ] Every new citation has a DOI
- [ ] No trajectory or private structure is attached, and any new test fixture states where it came from
- [ ] If this touches `src/`: `cargo test --lib` passes

What changed and why:

Reference-count impact (for `src/` changes: which cage, ring or order-parameter counts move, or "none"):
