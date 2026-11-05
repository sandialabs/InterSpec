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

"""Fit the Rayleigh-only ("elastic line") escape constant.

The model in SourceGeomFluxModel::evaluate is

    elastic = s * exp(-tau_abs) * (1 - exp(-tau_rs)) * exp(-d * tau_cs)

with `s` the payload's `elastic_escape_scale` and `d` its Compton damping
(`elastic_compton_damping`).  The line has a closed form and truth in every
bank cell, so it is fitted directly rather than enveloped.

Taking logs makes it exactly linear in (log s, d):

    log y + tau_abs - log(1 - exp(-tau_rs)) = log s - d * tau_cs

so this is a weighted linear least squares, not an optimizer. Weights are the
cells' own Poisson counts, so cells with a handful of escaping photons do not
outvote well-measured ones.

The free fit lands on d ~ 1, the physically natural value: a Compton
interaction removes the photon from the elastic line outright.  The payload
therefore pins d = 1 and takes `s` from the pinned solve; the free fit is
reported alongside as the check that pinning costs nothing.

Usage:
    fit_elastic_escape.py <campaign-dir> [<campaign-dir> ...]
        [--output elastic_escape.json] [--cross-section-helper PATH]
"""

from __future__ import annotations

import argparse
import csv
import bank_io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
# Matches the evaluator: the elastic line is only meaningful where Rayleigh
# actually escapes, and cells with almost no counts carry no information.
MINIMUM_COUNTS = 25.0


def load(campaigns: list[Path]):
    cells = []
    for campaign in campaigns:
        paths = bank_io.cell_files(campaign / "cells", recursive=True)
        if not paths:
            paths = bank_io.cell_files(campaign / "cells")
        for path in paths:
            with bank_io.open_cell(path) as stream:
                rows = list(csv.DictReader(stream))
            if not rows:
                continue
            first = rows[0]
            events = float(first["events"])
            total = 0.0
            for row in rows:
                if (row["component"] == "elastic_line"
                        and float(row["bin_hi_keV"])
                        <= float(row["bin_lo_keV"])):
                    total += float(row["sum_w"]) / events
            cells.append({
                "energy_keV": float(first["energy_keV"]),
                "Z": int(first["Z"]),
                "areal": float(first["areal_density_g_cm2"]),
                "f_H": float(
                    first.get("hydrogen_mass_fraction", 0.0) or 0.0),
                "events": events,
                "yield": total,
            })
    return cells


def attach_mu(cells, helper: Path) -> None:
    text = "".join(
        f"{c['energy_keV']} {c['Z']} {c['areal']} {c['f_H']}\n"
        for c in cells)
    completed = subprocess.run(
        [str(helper), "--with-hydrogen-column", "--smooth-z-blend"],
        input=text, text=True, stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    assert len(rows) == len(cells)
    for cell, row in zip(cells, rows):
        t = cell["areal"]
        cell["tau_abs"] = t * (float(row["mu_pe"]) + float(row["mu_pp"]))
        cell["tau_rs"] = t * float(row["mu_rs"])
        cell["tau_cs"] = t * float(row["mu_cs"])


def summary(relative: np.ndarray) -> dict[str, float]:
    return {
        "median": float(np.median(relative)),
        "abs_median": float(np.median(np.abs(relative))),
        "p95": float(np.percentile(np.abs(relative), 95)),
        "max": float(np.max(np.abs(relative))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign", type=Path, nargs="+",
                        help="bank directories holding cells/*_ceelo.csv.gz")
    parser.add_argument(
        "--output", type=Path,
        help="write the pinned scale, the free fit and the metrics as JSON")
    parser.add_argument(
        "--cross-section-helper", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_source_cross_section_features")
    args = parser.parse_args()
    cells = load(args.campaign)
    if not cells:
        print("no cells found", file=sys.stderr)
        return 2
    attach_mu(cells, args.cross_section_helper)

    usable = [
        c for c in cells
        if c["yield"] * c["events"] >= MINIMUM_COUNTS
        and c["tau_rs"] > 1e-9]
    print(f"{len(usable)} usable cells of {len(cells)}")
    if len(usable) < 2:
        print("too few usable cells to fit", file=sys.stderr)
        return 2

    counts = np.asarray([c["yield"] * c["events"] for c in usable])
    tau_abs = np.asarray([c["tau_abs"] for c in usable])
    tau_rs = np.asarray([c["tau_rs"] for c in usable])
    tau_cs = np.asarray([c["tau_cs"] for c in usable])
    observed = np.asarray([c["yield"] for c in usable])

    target = (
        np.log(observed) + tau_abs - np.log(-np.expm1(-tau_rs)))
    design = np.column_stack([np.ones_like(tau_cs), -tau_cs])
    # A log-space residual has variance ~ 1/counts, so weight by sqrt(counts).
    weight = np.sqrt(counts)
    solution, *_ = np.linalg.lstsq(
        design * weight[:, None], target * weight, rcond=None)
    scale = float(np.exp(solution[0]))
    damping = float(solution[1])
    # Pinned d = 1: one fewer free constant for the same accuracy.
    pinned_target = target + tau_cs
    pinned_scale = float(np.exp(
        np.sum(weight * weight * pinned_target)
        / np.sum(weight * weight)))

    def relative_error(s: float, d: float, mask=None) -> np.ndarray:
        select = np.ones_like(observed, dtype=bool) if mask is None else mask
        predicted = (
            s * np.exp(-tau_abs[select]) * (-np.expm1(-tau_rs[select]))
            * np.exp(-d * tau_cs[select]))
        return (predicted - observed[select]) / observed[select]

    deep = (tau_abs + tau_rs + tau_cs) > 6.0
    fits = {"free": (scale, damping), "pinned": (pinned_scale, 1.0)}
    metrics = {}
    print("\nelastic-line accuracy vs CeeLo truth:")
    for label, (s, d) in fits.items():
        metrics[label] = summary(relative_error(s, d))
        m = metrics[label]
        print(
            f"  {label:7s} s={s:.5f} d={d:.5f}  "
            f"median {m['median']*100:+7.2f}%  "
            f"|median| {m['abs_median']*100:6.2f}%  "
            f"p95 {m['p95']*100:7.2f}%  max {m['max']*100:8.1f}%")
        if deep.any():
            metrics[label]["deep_tau"] = summary(relative_error(s, d, deep))
            print(
                f"          tau>6 ({int(deep.sum())} cells) median "
                f"{metrics[label]['deep_tau']['median']*100:+8.2f}%")

    print(f"\nfree fit      : elastic_escape_scale {scale:.17g} "
          f"damping {damping:.17g}")
    print(f"pinned d=1    : elastic_escape_scale {pinned_scale:.17g}")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({
            "elastic_escape_scale": pinned_scale,
            "elastic_compton_damping": 1.0,
            "usable_cells": len(usable),
            "cells": len(cells),
            "deep_tau_cells": int(deep.sum()),
            "free_fit": {"elastic_escape_scale": scale,
                         "elastic_compton_damping": damping},
            "metrics": metrics,
        }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
