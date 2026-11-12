# Contributing a method

ClavisCage aims to offer every published method for finding hydrate cages and classifying water
structure, each one verified. Every method lives in this engine, so the web app, the Python package
and the Blender extension run the same implementation. A method describes itself through one descriptor in
[`src/methods.rs`](src/methods.rs): its citation, parameters, outputs, verification level and how it
can be shown. A new per-frame (series) method therefore appears, once released, in the web app's
block shelf, Results, CSV export and recipes and in Python's `run_method` without host-specific code.
Tetrahedral order was added this way. A new cage detector, and use from the command-line tool or the
Blender extension, still needs code in that host.

Bug reports, documentation fixes and method mismatches have their own issue forms. This page covers
new methods.

## Requesting a method

Open a [method request](https://github.com/wjgoarxiv/claviscage/issues/new?template=method_request.yml)
with:

| Field | What to give |
|-------|--------------|
| Paper | Full reference and DOI. If the method was refined later, give the version you want. |
| Inputs | What the method reads: oxygens only, hydrogens, guests, a trajectory rather than one frame. |
| Outputs | Per-frame numbers, per-water values or classes, or cages. |
| A fixture with its expected output | A small structure and the answer the method must give on it, ideally one that does not depend on any program (an ideal crystal, a hand-built cluster). |
| Code, if any | Repository URL and its licence. Say if there is no public code. |

A request is easier to accept when the expected output comes from the authors' own program or from
an established open implementation, because that sets the verification level.

## Verification levels

Each method carries a level, shown in the app and in its descriptor:

| Level | Meaning |
|-------|---------|
| **T3** | Matches the authors' own program on reference inputs (TRACE, HTR, GRADE). |
| **T2** | Matches an established open implementation within a stated tolerance (averaged Steinhardt vs freud; tetrahedral order vs GROMACS `gmx order`). |
| **T1** | Implemented from the paper and checked on cases whose answer is known by construction (CHILL+, F3/F4, MCG). |
| **T0** | Implemented from the paper only. |

A method's level rests on the fixtures listed in its descriptor and on the engine tests that check
them. The [benchmark](docs/validation/benchmarks.md) runs every method on ideal GenIce2 crystals.

## The spec-gap rule

Papers often leave details open: which neighbours count, how ties are broken, how a sum is
normalised. Before a paper-only method is written, its spec lists each open detail and the choice
made for it. These are the descriptor's `spec_gaps`. **A method with more than two open details is
not shipped** until the authors' code, an open reference or a correspondence with the authors closes
them.

## The licence rule

The engine is GPL-3.0-or-later, and contributions are accepted under the same licence.

- Code under a licence compatible with GPL-3.0-or-later (for example GPL-3.0, GPL-2.0-or-later,
  LGPL-2.1-or-later, MIT, 3-clause BSD or Apache-2.0) may be ported into the engine. Its notice is
  kept in the source file and listed in [`NOTICE`](NOTICE).
- Code with no licence, or with a licence that forbids it, is **never ported or redistributed**. It
  may be cited, and its published outputs may serve as reference values with their source named.
  HTR's public code, for example, declares no licence; the engine's HTR is checked against the
  authors' published program and does not include its code.
- Model weights follow the same rule as code.

## What a method needs to land

1. A descriptor in `src/methods.rs` with the citation, DOI, level, inputs, parameters with ranges,
   series, views, fixtures and spec gaps.
2. The implementation, reachable through `run_method`, with Rust tests on its fixtures
   (`cargo test --lib`).
3. Reference data under `tests/fixtures/reference/` with a provenance note (program, version,
   command).
4. A row in the benchmark table.
5. Python parity in `pytests/`, and the descriptor snapshot regenerated
   (`UPDATE_DESCRIPTOR_SNAPSHOT=1 cargo test --lib descriptor_snapshot`).

Pull requests follow the [pull request template](.github/pull_request_template.md). The method notes
in [`docs/methods/`](docs/methods/) show the expected level of detail.
