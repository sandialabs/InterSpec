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

"""Pre-fit gate on the source-secondary mechanistic bases.

The secondary fit stores log(truth / base) per node and multiplies the base
back in at evaluation, so the base's job is to carry the optical-depth
falloff.  What is left -- the conditional yield truth / base -- is what the
fitted kernel field must interpolate between tau nodes several bandwidths
apart.  This measures it on the truth banks BEFORE the fit is spent:

  span  = max/min of the conditional yield along one (E, Z, f_H) ray over tau
  slope = d log(conditional yield) / d tau, least squares along the ray

reported as medians over rays, per component.  The bases come from
`ceelo_source_cross_section_features --emit-component-bases`, the same
detail::source_component_bases the fitter divides out and the evaluator
multiplies back in, so this measures exactly the base the fit will use.

`pair_511_atom` calibrates the scale: its base is a full depth x exit-cosine
quadrature and its fitted field is flat, so its span/slope is what the banks'
own statistics allow.  A component passes at span <= 5 and |slope| <= 0.25.
`pair_annihilation_continuum` is the one component known to sit above that
bar (its base carries only one Compton scatter at 511 keV); it is held to its
own bound instead.  These bars apply to hydrogen-free rays.  Hydrogen-mixture
rays (f_H > 0) are gated separately, at the same bar except for the two
components whose bases are weakest there (`compton_electron_brems`, whose
parent energy is the maximum Compton recoil, and
`pair_annihilation_continuum`); a component with fewer than
HYDROGEN_MINIMUM_RAYS hydrogen rays is reported, not gated.

Usage:
    check_kernel_span.py --bank DIR [--bank DIR ...] [--helper PATH]
        [--output JSON]
Exit status 1 when a component is outside its bound.
"""

from __future__ import annotations

import argparse
import csv
import bank_io
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
COMPONENTS = ("pair_511_atom", "pair_annihilation_continuum",
              "pair_annihilation_in_flight", "photoelectron_brems",
              "compton_electron_brems", "pair_electron_brems",
              "fluorescence_k_alpha_atom", "fluorescence_k_beta_atom",
              "fluorescence_l_atom")
# The banks label the atoms without the `_atom` suffix the surrogate uses.
BANK_LABEL = {"pair_511_atom": "pair_511",
              "fluorescence_k_alpha_atom": "fluorescence_k_alpha",
              "fluorescence_k_beta_atom": "fluorescence_k_beta",
              "fluorescence_l_atom": "fluorescence_l"}
# (maximum median span, maximum |median slope|, origin of the bound)
DEFAULT_BOUND = (5.0, 0.25, "kernel acceptance bar: span <= 5, |slope| <= 0.25")
BOUNDS = {
    "pair_annihilation_continuum": (
        11.5, 0.36,
        "last accepted payload: span 9.2, slope +0.291 (x1.25)"),
    # Accepted at +0.249 against the 0.25 bar, so a new realization of the
    # same truth can land either side of it.
    "pair_annihilation_in_flight": (
        5.0, 0.31,
        "last accepted payload: span 4.10, slope +0.249 (x1.25)"),
    # The shield K x-ray bases follow a scattered primary through one Compton
    # collision; deep in a high-Z shield the multiply scattered photons that
    # are finally absorbed add K vacancies the base does not count.
    "fluorescence_k_alpha_atom": (
        5.0, 0.36,
        "span bank (Z 32-92, 59 keV-12 MeV): span 2.00, slope +0.285 (x1.25)"),
    "fluorescence_k_beta_atom": (
        5.0, 0.39,
        "span bank (Z 32-92, 59 keV-12 MeV): span 1.80, slope +0.314 (x1.25)"),
}
# Hydrogen-mixture rays (f_H > 0) are held separately: the bar above was set
# on hydrogen-free rays, and the bases are weakest where hydrogen makes the
# Compton recoil dominate.  A component with fewer rays than this is reported,
# not gated.
HYDROGEN_MINIMUM_RAYS = 10
HYDROGEN_BOUNDS = {
    "compton_electron_brems": (
        11.84, 0.35,
        "first gated with the hydrogen bank: span 9.47, slope +0.277 (x1.25)"),
    "pair_annihilation_continuum": (
        15.51, 0.39,
        "first gated with the hydrogen bank: "
        "span 12.41, slope +0.308 (x1.25)"),
}


def load_cell(path: Path):
    wanted = {BANK_LABEL.get(name, name): name for name in COMPONENTS}
    totals = dict.fromkeys(COMPONENTS, 0.0)
    variances = dict.fromkeys(COMPONENTS, 0.0)
    first = None
    with bank_io.open_cell(path) as stream:
        for row in csv.DictReader(stream):
            first = first or row
            name = wanted.get(row["component"])
            if name is not None:
                totals[name] += float(row["sum_w"])
                variances[name] += float(row["sum_w2"])
    if first is None:
        return None
    events = float(first["events"])
    return {
        "energy_keV": float(first["energy_keV"]),
        "Z": int(float(first["Z"])),
        "tau": float(first["tau"]),
        "areal_density": float(first["areal_density_g_cm2"]),
        "f_H": float(first.get("hydrogen_mass_fraction") or 0.0),
        "events": events,
        "yield": {k: v / events for k, v in totals.items()},
        "sigma": {k: math.sqrt(v) / events for k, v in variances.items()},
    }


def attach_bases(cells: list[dict], helper: Path) -> None:
    text = "".join(
        f"{c['energy_keV']!r} {c['Z']} {c['areal_density']!r} {c['f_H']!r}\n"
        for c in cells)
    completed = subprocess.run(
        [str(helper), "--emit-component-bases"], input=text, text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(COMPONENTS) * len(cells):
        raise RuntimeError("component-base helper returned wrong row count")
    for index, cell in enumerate(cells):
        block = rows[index * len(COMPONENTS):(index + 1) * len(COMPONENTS)]
        cell["base"] = {row["component"]: float(row["base"]) for row in block}


def ray_statistics(rays: dict) -> tuple[float, float, int]:
    spans, slopes = [], []
    for points in rays.values():
        if len(points) < 4:
            continue
        taus = np.asarray([tau for tau, _ in points])
        values = np.asarray([value for _, value in points])
        spans.append(float(values.max() / values.min()))
        slopes.append(float(np.polyfit(taus, np.log(values), 1)[0]))
    if not spans:
        return math.nan, math.nan, 0
    return float(np.median(spans)), float(np.median(slopes)), len(spans)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bank", type=Path, action="append", required=True,
                        help="bank directory holding cells/*_ceelo.csv.gz")
    parser.add_argument(
        "--helper", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_source_cross_section_features")
    parser.add_argument(
        "--max-relative-sigma", type=float, default=0.15,
        help="drop rows noisier than this; the deep-tau rows are a few "
             "counts per million and would otherwise dominate every span")
    parser.add_argument("--minimum-counts", type=float, default=25.0,
                        help="matches the fitter's --minimum-counts")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    cells = []
    for bank in args.bank:
        paths = bank_io.cell_files(bank / "cells")
        if not paths:
            print(f"no cells under {bank}", file=sys.stderr)
            return 2
        cells += [cell for cell in map(load_cell, paths) if cell is not None]
    attach_bases(cells, args.helper)

    # rays[hydrogen][name]: hydrogen-free (False) and hydrogen-mixture (True)
    rays: dict[bool, dict[str, dict]] = {
        hydrogen: {name: {} for name in COMPONENTS}
        for hydrogen in (False, True)}
    for cell in cells:
        for name in COMPONENTS:
            truth, base = cell["yield"][name], cell["base"][name]
            if truth <= 0.0 or base <= 0.0:
                continue
            if truth * cell["events"] < args.minimum_counts:
                continue
            if cell["sigma"][name] > args.max_relative_sigma * truth:
                continue
            key = (cell["energy_keV"], cell["Z"], cell["f_H"])
            rays[cell["f_H"] > 0.0][name].setdefault(key, []).append(
                (cell["tau"], truth / base))

    print(f"{len(cells)} cells; rows kept at sigma/yield <= "
          f"{args.max_relative_sigma:.0%} and >= {args.minimum_counts:g} "
          "counts; rays need >= 4 optical depths\n")
    report = {"cells": len(cells), "components": {},
              "hydrogen_components": {}}
    failures = 0
    for hydrogen, section in ((False, "components"),
                              (True, "hydrogen_components")):
        title = "hydrogen rays (f_H > 0)" if hydrogen else "hydrogen-free rays"
        print(f"{title:30s} {'rays':>5s} {'span':>7s} {'bound':>6s} "
              f"{'slope':>7s} {'bound':>6s}  verdict")
        for name in COMPONENTS:
            span, slope, count = ray_statistics(rays[hydrogen][name])
            bounds = HYDROGEN_BOUNDS if hydrogen else BOUNDS
            max_span, max_slope, origin = bounds.get(name, DEFAULT_BOUND)
            if count == 0:
                verdict = "REPORTED" if hydrogen else "NO DATA"
            elif hydrogen and count < HYDROGEN_MINIMUM_RAYS:
                verdict = "REPORTED"
            elif span <= max_span and abs(slope) <= max_slope:
                verdict = "PASS"
            else:
                verdict = "FAIL"
            failures += verdict not in ("PASS", "REPORTED")
            print(f"{name:30s} {count:5d} {span:7.2f} {max_span:6.2f} "
                  f"{slope:+7.3f} {max_slope:6.2f}  {verdict}")
            report[section][name] = {
                "rays": count, "median_span": span, "median_slope": slope,
                "max_span": max_span, "max_abs_slope": max_slope,
                "bound_origin": origin, "verdict": verdict,
                "pass": verdict in ("PASS", "REPORTED")}
        print()
    report["pass"] = failures == 0
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print("\nPASS" if failures == 0 else f"\nFAIL ({failures} component(s))")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
