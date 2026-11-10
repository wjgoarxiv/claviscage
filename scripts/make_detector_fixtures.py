# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Woojin Go, Junghoon Mok and Yongwon Seo
"""Regenerate reference data for the averaged Steinhardt (q-bar) tests.

The reference values come from freud (BSD-3-Clause), an established open
implementation of the Lechner-Dellago averaged Steinhardt order parameter.
freud is a test-data tool only; the engine does not depend on it.

Usage (from the engine directory, in a throwaway environment):

    python -m venv .venv-ref && .venv-ref/bin/pip install freud-analysis numpy
    .venv-ref/bin/python scripts/make_detector_fixtures.py

Neighbour rule (identical to the engine): every other water oxygen closer
than 0.35 nm under the minimum-image convention; q-bar averages q_lm over the
water itself and those neighbours (Lechner & Dellago 2008, Eq. 3).
"""

import hashlib
import json
from pathlib import Path

import freud
import numpy as np

ENGINE = Path(__file__).resolve().parents[1]
CUTOFF_NM = 0.35
LS = (4, 6)
WATER = {"SOL", "WAT", "HOH", "H2O", "TIP3", "TIP3P", "TIP4", "TIP4P", "SPC", "SPCE"}
INPUTS = [
    ("crystals/444_S1.gro", "steinhardt_444_S1.json"),
    ("crystals/222_S2.gro", "steinhardt_222_S2.json"),
    ("guests/sI_CH4_222.gro", "steinhardt_sI_CH4_222.json"),
    ("reference/grade_prod.gro", "steinhardt_grade_prod.json"),
]


def read_oxygens(path):
    lines = path.read_text().splitlines()
    n = int(lines[1])
    oxygens = []
    for line in lines[2 : 2 + n]:
        res = line[5:10].strip().upper()
        atom = line[10:15].strip().upper()
        if res in WATER and (atom.startswith("O") or "OW" in atom):
            oxygens.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    box = [float(v) for v in lines[2 + n].split()[:3]]
    return np.asarray(oxygens, dtype=np.float64), box


def main():
    for rel, out in INPUTS:
        src = ENGINE / "tests" / "fixtures" / rel
        pts, box = read_oxygens(src)
        fbox = freud.box.Box(*box)
        wrapped = fbox.wrap(pts - np.asarray(box) / 2.0)
        result = {
            "reference": f"freud {freud.__version__} freud.order.Steinhardt(average=True)",
            "input": f"tests/fixtures/{rel}",
            "input_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
            "cutoff_nm": CUTOFF_NM,
            "n_waters": int(len(pts)),
        }
        for l in LS:
            q = freud.order.Steinhardt(l, average=False)
            q.compute((fbox, wrapped), neighbors={"r_max": CUTOFF_NM, "exclude_ii": True})
            qb = freud.order.Steinhardt(l, average=True)
            qb.compute((fbox, wrapped), neighbors={"r_max": CUTOFF_NM, "exclude_ii": True})
            result[f"q{l}"] = [round(float(v), 10) for v in np.nan_to_num(q.particle_order)]
            result[f"qbar{l}"] = [round(float(v), 10) for v in np.nan_to_num(qb.particle_order)]
        dst = ENGINE / "tests" / "fixtures" / "reference" / out
        dst.write_text(json.dumps(result) + "\n")
        print(f"{out}: {len(pts)} waters, mean qbar6 {np.mean(result['qbar6']):.4f}")


if __name__ == "__main__":
    main()
