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

"""Measure the source-secondary model's escaping yield BEYOND its trained taus.

An opaque shield emits nothing, so past the deepest trained optical depth
(tau = 8) the escaping secondary yield must fall to ~0.  A kernel field left
alone extrapolates flat off the deepest populated node instead; the fit
counters that with Poisson-bound anchor nodes and the tail-decay envelope
exp(-kappa (tau - tau_last)).  This check evaluates the model at optical
depths far past the grid and reports the worst spurious secondary yield per
primary -- an ABSOLUTE yield, since the truth there is ~0 and a relative error
would divide by nothing.  The off-grid holdout cannot see this defect: all its
cells sit inside the trained range.

Usage:
    check_deep_tau_tail.py MODEL [MODEL ...] [--evaluator PATH] [--output JSON]

MODEL is an assembled CEELO_SOURCE_ESCAPE payload or a loose
CEELO_SOURCE_SECONDARY_SURROGATE file.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from source_escape_groups import SECONDARY, read_rows, totals

REPO = Path(__file__).resolve().parents[2]
# 8 is the deepest trained node; everything past it is extrapolation.
TAUS = (8.0, 12.0, 20.0, 50.0, 100.0)
CELLS = ((250.0, 82), (662.0, 26), (3000.0, 13), (3000.0, 82))


def model_arguments(path: Path) -> list[str]:
    with path.open() as stream:
        magic = stream.readline().split()
    if magic and magic[0] == "CEELO_SOURCE_ESCAPE":
        return ["--payload", str(path)]
    if magic and magic[0] == "CEELO_SOURCE_SECONDARY_SURROGATE":
        return ["--secondary-surrogate", str(path)]
    raise ValueError(f"{path} is neither a payload nor a secondary surrogate")


def secondary_yield(evaluator: Path, model: list[str], energy: float,
                    z_value: int, tau: float) -> float:
    completed = subprocess.run(
        [str(evaluator), *model, "--energy-kev", repr(energy),
         "--z", str(z_value), "--tau", repr(tau), "--bins", "128"],
        stdout=subprocess.PIPE, text=True, check=True)
    return totals(read_rows(completed.stdout), SECONDARY)[0]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model", type=Path, nargs="+")
    parser.add_argument(
        "--evaluator", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_source_model_evaluate")
    parser.add_argument("--output", type=Path,
                        help="write the table and the worst yields as JSON")
    args = parser.parse_args()
    for path in (args.evaluator, *args.model):
        if not path.exists():
            print(f"missing: {path}", file=sys.stderr)
            return 2

    names = [str(path) for path in args.model]
    models = [model_arguments(path) for path in args.model]
    print("escaping secondary yield per primary (truth ~ 0 for tau >> 8)\n")
    print(f"{'cell':>14s} {'tau':>6s} " + "".join(
        f"{f'model {index}':>14s}" for index in range(len(models))))
    table = []
    worst = [0.0] * len(models)
    for energy, z_value in CELLS:
        for tau in TAUS:
            values = [secondary_yield(args.evaluator, model, energy,
                                      z_value, tau) for model in models]
            for index, value in enumerate(values):
                if tau > 8.0:
                    worst[index] = max(worst[index], value)
            table.append({"energy_keV": energy, "Z": z_value, "tau": tau,
                          "secondary_yield": values})
            print(f"{f'E{energy:g} Z{z_value}':>14s} {tau:6.0f} "
                  + "".join(f"{value:14.3e}" for value in values))
    print("\nworst spurious secondary yield past the trained grid (tau > 8):")
    for index, name in enumerate(names):
        print(f"  model {index}: {worst[index]:.3e} photons/primary  ({name})")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({
            "models": names,
            "worst_spurious_yield": worst,
            "cells": table,
        }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
