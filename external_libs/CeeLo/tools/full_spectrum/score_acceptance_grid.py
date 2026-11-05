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

"""Score a source-escape payload on the acceptance grid, against CeeLo and Geant4.

The grid (written by `run_source_anchor.py --acceptance-grid`) holds a CeeLo
truth CSV and a Geant4 CSV per cell.  Only the MODEL leg is computed here, so
scoring a payload costs minutes rather than a Geant4 campaign; model spectra
are cached per payload (keyed by its SHA-256), so re-scoring the same payload
is free.

The metric is the collided-only leakage per primary: every escaping group of
`source_escape_groups.GROUPS` except the analytic, exact-by-construction
uncollided term (--with-uncollided restores it).  Cells whose truth leakage is
below MINIMUM_SCORED_LEAKAGE are skipped: that removes the unshielded cells and
the opaque ones, where a relative error is undefined or measures only the
truth's own noise.  The floor is a fixed leakage, so the scored cells do not
change with the grid's statistics.

Usage:
    score_acceptance_grid.py --grid DIR --payload FILE [--evaluator PATH]
        [--output JSON] [--workers N] [--with-uncollided]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from source_escape_groups import KNOWN, read_rows, totals
import bank_io

REPO = Path(__file__).resolve().parents[2]
IDENTITY = re.compile(
    r"^E(?P<energy>[\dp]+)_Z(?P<z>\d+)_rhoA(?P<rho>[\dp]+)$")
# Collided leakage per primary below which a cell is not scored: 25 counts at
# 500k histories.
MINIMUM_SCORED_LEAKAGE = 5.0e-5


def number(text: str) -> float:
    return float(text.replace("p", "."))


def score_cell(identity: str, grid: Path, payload: Path, evaluator: Path,
               cache: Path, labels) -> dict | None:
    match = IDENTITY.match(identity)
    ceelo_path = bank_io.cell_path(grid / "cells", identity)
    if match is None or ceelo_path is None:
        return None
    energy = number(match.group("energy"))
    z_value = int(match.group("z"))
    rho = number(match.group("rho"))
    ceelo_rows = read_rows(ceelo_path)
    # Like-for-like binning: the model uses the grid CSVs' own bin count.
    bins = max(64, len({
        (row["bin_lo_keV"], row["bin_hi_keV"]) for row in ceelo_rows
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])}))
    model_path = cache / f"{identity}_model.csv"
    if not model_path.exists():
        completed = subprocess.run(
            [str(evaluator), "--payload", str(payload),
             "--energy-kev", repr(energy), "--z", str(z_value),
             "--areal-density", repr(rho), "--bins", str(bins)],
            stdout=subprocess.PIPE, text=True, check=True)
        temporary = model_path.with_suffix(".csv.tmp")
        temporary.write_text(completed.stdout)
        temporary.replace(model_path)
    model, _, ungrouped = totals(read_rows(model_path), labels)
    ceelo, ceelo_sigma, _ = totals(ceelo_rows, labels)
    record = {
        "identity": identity, "energy_keV": energy, "Z": z_value,
        "areal_density_g_cm2": rho, "model": model,
        "ceelo": ceelo, "ceelo_sigma": ceelo_sigma,
        "ungrouped_model_labels": sorted(ungrouped),
    }
    geant4_path = bank_io.cell_path(grid / "cells", identity, "geant4")
    if geant4_path is not None:
        geant4, geant4_sigma, _ = totals(read_rows(geant4_path), labels)
        record["geant4"] = geant4
        record["geant4_sigma"] = geant4_sigma
    return record


def summarize(records: list[dict], key: str) -> dict:
    usable = [r for r in records
              if r.get(key, 0.0) >= MINIMUM_SCORED_LEAKAGE]
    if not usable:
        return {"cells": 0}
    errors = np.asarray([(r["model"] - r[key]) / r[key] for r in usable])
    pulls = [abs(r["model"] - r[key]) / r[f"{key}_sigma"]
             for r in usable if r.get(f"{key}_sigma", 0.0) > 0.0]
    return {
        "cells": len(usable),
        "median": float(np.median(errors)),
        "abs_median": float(np.median(np.abs(errors))),
        "p95": float(np.percentile(np.abs(errors), 95)),
        "max": float(np.max(np.abs(errors))),
        "abs_z_median": float(np.median(pulls)) if pulls else None,
    }


def line(title: str, values: dict) -> str:
    if not values.get("cells"):
        return f"  {title}: no usable cells"
    z = values.get("abs_z_median")
    return (f"  {title}: n={values['cells']:3d}  "
            f"median {100 * values['median']:+7.2f}%  "
            f"|median| {100 * values['abs_median']:6.2f}%  "
            f"p95 {100 * values['p95']:7.2f}%  "
            f"max {100 * values['max']:8.2f}%"
            + (f"  |z| median {z:5.1f}" if z is not None else ""))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument(
        "--evaluator", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_source_model_evaluate")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--model-cache", type=Path,
        help="where model spectra are kept (default: "
             "GRID/models/<payload sha256 prefix>)")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--with-uncollided", action="store_true",
        help="include the analytic uncollided term (dilutes the metric)")
    args = parser.parse_args()

    digest = hashlib.sha256(args.payload.read_bytes()).hexdigest()
    cache = args.model_cache or args.grid / "models" / digest[:16]
    cache.mkdir(parents=True, exist_ok=True)
    labels = KNOWN if args.with_uncollided else None
    identities = sorted({
        bank_io.cell_id(path)
        for path in bank_io.cell_files(args.grid / "cells")})
    if not identities:
        print(f"no cells under {args.grid}", file=sys.stderr)
        return 2
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = [record for record in pool.map(
            lambda identity: score_cell(
                identity, args.grid, args.payload, args.evaluator, cache,
                labels),
            identities) if record is not None]

    ungrouped = sorted({name for r in records
                        for name in r["ungrouped_model_labels"]})
    if ungrouped:
        print(f"WARNING: model components outside GROUPS (not counted): "
              f"{ungrouped}")
    metric = ("total incl. uncollided" if args.with_uncollided
              else "collided-only leakage")
    report = {
        "payload": str(args.payload),
        "payload_sha256": digest,
        "metric": metric,
        "cells": len(records),
        "vs_ceelo": summarize(records, "ceelo"),
        "vs_geant4": summarize(records, "geant4"),
        "by_areal_density_vs_ceelo": {
            f"{rho:g}": summarize(
                [r for r in records if r["areal_density_g_cm2"] == rho],
                "ceelo")
            for rho in sorted({r["areal_density_g_cm2"] for r in records})},
        "records": records,
    }
    print(f"acceptance grid: {len(records)} cells, payload sha256 "
          f"{digest[:16]}; metric: {metric}")
    print(line("model vs CeeLo ", report["vs_ceelo"]))
    print(line("model vs Geant4", report["vs_geant4"]))
    print("  by areal density (vs CeeLo):")
    for rho, values in report["by_areal_density_vs_ceelo"].items():
        print("  " + line(f"{rho:>4s} g/cm2", values))
    ranked = sorted((r for r in records
                     if r["ceelo"] >= MINIMUM_SCORED_LEAKAGE),
                    key=lambda r: abs(r["model"] / r["ceelo"] - 1.0),
                    reverse=True)[:10]
    print("  worst cells vs CeeLo:")
    for r in ranked:
        sigma = r["ceelo_sigma"]
        z = (r["model"] - r["ceelo"]) / sigma if sigma > 0 else float("nan")
        print(f"    {r['identity']:20s} model {r['model']:.4e} truth "
              f"{r['ceelo']:.4e}  {100 * (r['model'] / r['ceelo'] - 1):+9.2f}%"
              f"  z={z:+7.1f}")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
