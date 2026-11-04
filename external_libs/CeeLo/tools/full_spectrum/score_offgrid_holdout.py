#!/usr/bin/env python3
# CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
# and X-ray detector efficiency - developed as part of InterSpec.
#
# Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

"""Score an assembled source-escape payload on the off-grid tau holdout.

The off-grid bank's optical depths all sit BETWEEN the secondary fit's
training nodes and it is never trained on, so it is the only measure of
between-node interpolation.  Two metrics, each a relative error against the
bank's CeeLo truth with a z against the truth cell's own MC sigma (the model
side is deterministic):

1. collided-only leakage -- every escaping group except the analytic,
   exact-by-construction uncollided term, which would otherwise dominate thin
   cells and flatter the score;
2. the bremsstrahlung block (the three electron-brems channels plus in-flight
   annihilation) at tau >= 3, where the escape kernel does its work;
3. the fluorescence atoms (shield K-alpha, K-beta and L x-rays, summed) in
   cells whose truth holds at least 25 of them; the Compton-scattered x-ray
   continuum the model does not produce is reported as a share of them.

Usage:
    score_offgrid_holdout.py --offgrid DIR --payload FILE [--evaluator PATH]
        [--output JSON] [--workers N]
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from source_escape_groups import BREMS, FLUORESCENCE, UNMODELLED, read_rows, totals
import bank_io

REPO = Path(__file__).resolve().parents[2]
DEEP_TAU = 3.0


def score_cell(path: Path, payload: Path, evaluator: Path) -> dict:
    rows = read_rows(path)
    first = rows[0]
    leak, leak_sigma, _ = totals(rows)
    brems, brems_sigma, _ = totals(rows, BREMS)
    fluorescence, fluorescence_sigma, _ = totals(rows, FLUORESCENCE)
    completed = subprocess.run(
        [str(evaluator), "--payload", str(payload),
         "--energy-kev", first["energy_keV"], "--z", first["Z"],
         "--areal-density", first["areal_density_g_cm2"],
         "--hydrogen-mass-fraction",
         first.get("hydrogen_mass_fraction") or "0",
         "--bins", "256"],
        stdout=subprocess.PIPE, text=True, check=True)
    model_rows = read_rows(completed.stdout)
    return {
        "identity": bank_io.cell_id(path),
        "energy_keV": float(first["energy_keV"]),
        "Z": int(float(first["Z"])),
        "tau": float(first["tau"]),
        "leak": leak, "leak_sigma": leak_sigma,
        "model_leak": totals(model_rows)[0],
        "brems": brems, "brems_sigma": brems_sigma,
        "model_brems": totals(model_rows, BREMS)[0],
        "fluorescence": fluorescence, "fluorescence_sigma": fluorescence_sigma,
        "model_fluorescence": totals(model_rows, FLUORESCENCE)[0],
        "fluorescence_events": fluorescence * float(first["events"]),
        "fluorescence_continuum": totals(rows, UNMODELLED)[0],
    }


def summary(records: list[dict], key: str) -> dict:
    errors, pulls = [], []
    for record in records:
        truth, sigma = record[key], record[f"{key}_sigma"]
        difference = record[f"model_{key}"] - truth
        errors.append(abs(difference) / truth)
        if sigma > 0.0:
            pulls.append(abs(difference) / sigma)
    if not errors:
        return {"cells": 0}
    return {
        "cells": len(errors),
        "median": float(np.median(errors)),
        "p95": float(np.percentile(errors, 95)),
        "max": float(np.max(errors)),
        "abs_z_median": float(np.median(pulls)) if pulls else math.nan,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--offgrid", type=Path, required=True,
                        help="off-grid bank directory (cells/*_ceelo.csv.gz)")
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument(
        "--evaluator", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_source_model_evaluate")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    paths = bank_io.cell_files(args.offgrid / "cells")
    if not paths:
        print(f"no off-grid cells under {args.offgrid}", file=sys.stderr)
        return 2
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(
            lambda path: score_cell(path, args.payload, args.evaluator),
            paths))
    leaking = [r for r in records if r["leak"] > 0.0]
    deep = [r for r in records if r["tau"] >= DEEP_TAU and r["brems"] > 0.0]
    fluorescent = [r for r in records if r["fluorescence_events"] >= 25.0]
    report = {
        "cells": len(records),
        "leakage": summary(leaking, "leak"),
        "deep_tau_brems": summary(deep, "brems"),
        "fluorescence": summary(fluorescent, "fluorescence"),
        "fluorescence_continuum_share": (
            {"cells": len(fluorescent),
             "median": float(np.median([r["fluorescence_continuum"] / r["fluorescence"]
                                        for r in fluorescent])),
             "max": float(np.max([r["fluorescence_continuum"] / r["fluorescence"]
                                  for r in fluorescent]))}
            if fluorescent else {"cells": 0}),
        "deep_tau": DEEP_TAU,
        "records": records,
    }
    for label, key in (("collided-only leakage", "leakage"),
                       (f"brems block at tau >= {DEEP_TAU:g}",
                        "deep_tau_brems"),
                       ("fluorescence atoms", "fluorescence")):
        values = report[key]
        if not values.get("cells"):
            print(f"{label:34s} no cells")
            continue
        print(f"{label:34s} n={values['cells']:3d}  "
              f"median {100 * values['median']:6.2f}%  "
              f"p95 {100 * values['p95']:7.2f}%  "
              f"max {100 * values['max']:7.2f}%  "
              f"|z| median {values['abs_z_median']:5.1f}")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
