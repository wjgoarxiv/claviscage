# `.claviscage.zip` bundle

The web app's **Blender bundle** export is a single-frame handoff to the ClavisCage Blender
extension. This page describes its shape so you can audit an export; it is not an editing API.
Machine-readable draft: `schema/claviscage-bundle.v2.schema.json`.

## Archive

- File name ends in `.claviscage.zip`.
- Exactly two members: `manifest.json` and `scene.json`. Importers reject missing members, extra
  members and unsupported schema versions.
- One displayed frame. Not a trajectory, not a `.blend`, not a rendered image.

## `manifest.json`

| Field | Value |
|---|---|
| `format` | `"claviscage.bundle"` |
| `schemaVersion` | `2` for files the web app writes now |
| `appVersion` | string |
| `createdAt` | string (timestamp) |
| `warnings` | array of strings; analyses that belonged to another frame are omitted and listed here |

## `scene.json`

- `frame.title`
- `frame.waters`, `frame.guests`, `frame.molecules` (`name`, `elements`, `positions`)
- `frame.boxDims` (3 numbers), `frame.cell` (6 numbers), `frame.isTriclinic`
- `analyses.chill` or `null`
- `analyses.adjacency`: neighbour index lists
- `analyses.cageAnalysis` or `null`

`cageAnalysis`, when present:

- `algorithm`: `"trace"`, `"htr"` or `"grade"` (CHILL+ is never this field)
- `frameIndex`: zero-based source frame
- `totalCages`, `visibleCages`
- `parameters`: `oxygenCutoffNm`, and where relevant `donorAngleDeg`, `planarityDeg`, `d1Nm`, `d2Nm`
- `cages`: geometry owned by that algorithm and that frame

## Versions

Write **v2**. The Blender extension also reads legacy **v1**, whose cage field is HTR-specific and
cannot record TRACE or GRADE provenance.
