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

"""Adaptively increase MC statistics in populated multi-scatter cells."""

from __future__ import annotations

import argparse
import csv
import bank_io
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from validate_multiple_scatter import (
    TRUTH_COMPONENTS,
    run_cell,
)


def inspect_cell(path: Path) -> dict:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty campaign cell: {path}")
    continuous = [
        row for row in rows
        if row["component"] in TRUTH_COMPONENTS
        and float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])]
    histories = sum(float(row["sum_w"]) for row in continuous)
    energy = float(rows[0]["energy_keV"])
    bins = sum(
        row["component"] == TRUTH_COMPONENTS[0]
        and float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])
        for row in rows)
    return {
        "energy": energy,
        "Z": int(float(rows[0]["Z"])),
        "tau": float(rows[0]["tau"]),
        "events": int(float(rows[0]["events"])),
        "histories": histories,
        "bins": bins,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--target-histories", type=float, default=1000.0)
    parser.add_argument(
        "--minimum-pilot-histories", type=float, default=20.0)
    parser.add_argument("--max-events", type=int, default=2000000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--z-values", type=int, nargs="*")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--source-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_escape_campaign"))
    parser.add_argument(
        "--model-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_model_evaluate"))
    args = parser.parse_args()

    selected = []
    for path in bank_io.cell_files(args.campaign / "cells"):
        cell = inspect_cell(path)
        if args.z_values and cell["Z"] not in args.z_values:
            continue
        if not (
            args.minimum_pilot_histories <= cell["histories"]
            < args.target_histories
        ):
            continue
        estimated_yield = cell["histories"] / cell["events"]
        requested = min(
            args.max_events,
            max(
                cell["events"],
                int(math.ceil(
                    1.2 * args.target_histories / estimated_yield))))
        if requested <= cell["events"]:
            # Already at the history cap: re-running would reproduce the
            # same cell, so a resumed refinement skips it.
            continue
        selected.append((cell, requested))

    total_events = sum(events for _, events in selected)
    print(
        f"refining {len(selected)} cells; "
        f"{total_events} capped deterministic histories")
    if args.dry_run:
        for cell, events in selected:
            print(
                f"E={cell['energy']:g} Z={cell['Z']} "
                f"tau={cell['tau']:g}: "
                f"{cell['histories']:.0f} -> {events} events")
        return 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(
                run_cell, args.campaign,
                cell["energy"], cell["Z"], cell["tau"],
                cell["bins"], events,
                args.source_driver, args.model_driver,
                True, False, None, False)
            for cell, events in selected
        ]
        for future in as_completed(futures):
            identity, state = future.result()
            print(state, identity)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
