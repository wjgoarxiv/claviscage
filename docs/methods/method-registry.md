# Method registry and verification levels

The engine describes every published method it implements. The registry is the same in
WebAssembly and Python:

- `list_methods()` returns all descriptors.
- `describe_method(id)` returns one descriptor, or an error for an unknown id.
- `run_method(id, ...)` runs a series method by id with parameters named in its descriptor.

A descriptor holds the method's id and name, citation and DOI, inputs it needs (hydrogens, guests),
parameters with default, range and step, the per-frame series and per-water values it returns, how
the result can be shown (metric chart, per-water colouring, cage view), the fixtures its
verification rests on, details the paper leaves open and the choice made for each, the licence of
this implementation, and a `descriptor_version` that is raised when a method's parameters, outputs
or numbers change.

`run_method` covers the series methods: `chill`, `f3f4`, `steinhardt`, `mcg` and `tetra`. Cage
detectors (`trace`, `htr`, `grade`) keep their own entry points. A call returns the descriptor's
series, optionally the per-water values, and the method's own result. Errors start with a code:
`unknown_method`, `not_series_method`, `unknown_param`, `param_out_of_range`, `not_integer`,
`missing_input` or `invalid`.

```python
import json, clavis_engine
methods = json.loads(clavis_engine.list_methods())
info = json.loads(clavis_engine.describe_method("tetra"))
# result = clavis_engine.run_method("tetra", waters_json, guests_json, box, params_json=None)
```

## Verification levels

| Level | Meaning |
|---|---|
| **T3** | Results match the authors' own program on the reference inputs. |
| **T2** | Results match an established open implementation. |
| **T1** | Implemented from the paper and checked on cases with known answers (crystals, hand-built systems). |
| **T0** | Implemented from the paper only. |

A level above T0 names the fixtures it rests on. The current levels:

| Method | Id | Level | Verified against |
|---|---|---|---|
| TRACE | `trace` | T3 | Port of the authors' TRACE program; exact cage counts on sI and sII crystals |
| HTR | `htr` | T3 | Implemented from the paper and checked against the authors' published program (its code is not included here) |
| GRADE | `grade` | T3 | Step-by-step port of GRADE v1.00; equal counts on six reference inputs |
| CHILL+ | `chill` | T1 | Paper definitions and documented thresholds; every water of the reference crystals is classed as hydrate |
| F3/F4 | `f3f4` | T1 | Published definitions; F4 pair rule follows the CalF4 program distributed with GRADE |
| Averaged Steinhardt q4/q6 | `steinhardt` | T2 | freud `Steinhardt(average=True)`, within its float32 precision |
| MCG-3 | `mcg` | T1 | Paper definitions; hand-built pairs and fully occupied sI and sII crystals |
| Tetrahedral order | `tetra` | T2 | GROMACS 2026.2 `gmx order -Sg -Sk`, frame means within 2e-6 |

The benchmark that runs every method on GenIce2 sI, sII, sH, Ih and Ic crystals is described in
[`../validation/benchmarks.md`](../validation/benchmarks.md).
