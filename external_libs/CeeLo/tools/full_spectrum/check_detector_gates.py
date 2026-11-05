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

"""Acceptance gates for a baked HPGe detector-response payload.

Reads the reports `fit_detector_response.py` wrote for the candidate (scored
on the held-out validation bank) and for the audit banks, plus the direct-MC
spot checks, and prints one table.  Exit status 0 means every REQUIRED gate
passed; 1 means at least one failed.

Required gates (validation bank):
  detected yield      relative error median <= 0.15%, p95 <= 0.5%
  pulse-shape TV      median <= 1.3%, p95 <= 4%
  full-energy atom    relative error median <= 0.5%, p95 <= 1.5%
  spot checks         |model FEP / direct-MC FEP - 1| <= 0.5% at every
                      spot energy (normal incidence; the truth's own sigma
                      is printed and must be well below the tolerance)

Required gates (arrival bank, `--arrival`, score_arrival_bank.py's report),
in every arrival class (iso, front, back, side) and energy band (40-100,
100-400, 400-1000 keV):
  detected yield      |relative error| median <= 0.15%; max <= 0.5% beyond
                      2 sigma of each cell's MC noise
  continuum TV        (no zero-deposit or full-energy atom) median <= 3x its
                      MC noise floor
  valley density      |relative error| median <= 3%, max <= 10%
  envelope area       |z| max <= 5 (MC boundary against the payload envelope)
The arrival full-energy atom is reported, not gated, until the fitter's
full-energy misfit is fixed (TODO.md).

The audit banks are reported, not gated: they sit off every training grid,
including near-grazing incidence, and are for judging a candidate against the
payload it would replace.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

# (metric path in the report's "validation" section, statistic, bound)
REQUIRED = (
    (("detected_yield_relative_error",), "median", 0.0015),
    (("detected_yield_relative_error",), "p95", 0.005),
    (("pulse_shape_total_variation",), "median", 0.013),
    (("pulse_shape_total_variation",), "p95", 0.04),
    (("component_yield_relative_error", "full_energy"), "median", 0.005),
    (("component_yield_relative_error", "full_energy"), "p95", 0.015),
)
SPOT_TOLERANCE = 0.005
# (metric in each arrival group, statistic, bound)
ARRIVAL_REQUIRED = (
    ("detected_yield_relative_error", "median", 0.0015),
    ("detected_yield_excess", "max", 0.005),
    ("continuum_tv_over_noise", "median", 3.0),
    ("valley_density_relative_error", "median", 0.03),
    ("valley_density_relative_error", "max", 0.10),
    ("envelope_area_z", "max", 5.0),
)
ARRIVAL_REPORTED = ("full_energy_relative_error", "pulse_shape_total_variation")
# Informational rows printed for the candidate and every audit bank.
REPORTED = (
    ("detected_yield_relative_error",),
    ("pulse_shape_total_variation",),
    ("continuum_shape_total_variation",),
    ("component_yield_relative_error", "full_energy"),
)


def lookup(section: dict, path: tuple[str, ...]) -> dict | None:
    node = section
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def validation_section(report_path: Path) -> dict:
    report = json.loads(report_path.read_text())
    section = report.get("validation")
    if not isinstance(section, dict):
        raise SystemExit(f"{report_path}: no 'validation' section")
    return section


def full_energy(path: Path) -> tuple[float, float]:
    """Full-energy yield per incident photon and its MC sigma."""
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["component"] == "full_energy":
                events = float(row["events"])
                value = float(row["sum_w"]) / events
                sigma = math.sqrt(max(float(row["sum_w2"]), 0.0)) / events
                return value, sigma
    raise SystemExit(f"{path}: no full_energy row")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--candidate", type=Path, required=True,
                        help="report.json of the candidate on the validation bank")
    parser.add_argument("--audit", action="append", default=[],
                        metavar="NAME=REPORT",
                        help="report.json of the candidate on an audit bank")
    parser.add_argument("--spot", type=Path, required=True,
                        help="directory of truth_<E>.csv / model_<E>.csv pairs")
    parser.add_argument("--arrival", type=Path,
                        help="score_arrival_bank.py report.json of the candidate")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    checks = []
    validation = validation_section(args.candidate)
    for path, statistic, bound in REQUIRED:
        stats = lookup(validation, path) or {}
        value = stats.get(statistic)
        passed = value is not None and value <= bound
        checks.append({
            "gate": f"validation {'.'.join(path)} {statistic}",
            "value": value, "bound": bound, "passed": passed})

    truths = sorted(args.spot.glob("truth_*.csv"))
    if not truths:
        raise SystemExit(f"no spot checks under {args.spot}")
    for truth_path in truths:
        energy = truth_path.stem.removeprefix("truth_")
        model_path = args.spot / f"model_{energy}.csv"
        truth, sigma = full_energy(truth_path)
        model, _ = full_energy(model_path)
        error = model / truth - 1.0
        checks.append({
            "gate": f"spot check FEP {energy} keV (truth sigma "
                    f"{sigma / truth:.3%})",
            "value": abs(error), "signed": error, "bound": SPOT_TOLERANCE,
            "truth": truth, "model": model,
            "passed": abs(error) <= SPOT_TOLERANCE})

    arrival_reported = []
    if args.arrival is not None:
        groups = json.loads(args.arrival.read_text()).get("groups") or []
        if not groups:
            raise SystemExit(f"{args.arrival}: no arrival groups")
        for group in groups:
            lo, hi = group["energy_band_keV"]
            label = (f"arrival {group['incidence']} {lo:g}-{hi:g} keV"
                     if hi is not None else
                     f"arrival {group['incidence']} >= {lo:g} keV")
            for metric, statistic, bound in ARRIVAL_REQUIRED:
                value = (group.get(metric) or {}).get(statistic)
                checks.append({
                    "gate": f"{label} {metric} {statistic}",
                    "value": value, "bound": bound,
                    "ratio": metric.endswith("_over_noise")
                    or metric.endswith("_z"),
                    "passed": value is not None and value <= bound})
            arrival_reported.append({
                "group": label, "cells": group["cells"],
                **{metric: group.get(metric) for metric in ARRIVAL_REPORTED}})

    reported = {"validation": {}}
    banks = [("validation", validation)]
    for item in args.audit:
        name, _, report = item.partition("=")
        section = validation_section(Path(report))
        reported[name] = {}
        banks.append((name, section))
    for name, section in banks:
        for path in REPORTED:
            stats = lookup(section, path) or {}
            reported[name][".".join(path)] = {
                key: stats.get(key) for key in ("count", "median", "p95")}

    print(f"{'gate':76s} {'value':>9s} {'bound':>8s}  result")
    for check in checks:
        value = check["value"]
        if check.get("ratio"):
            text = "missing" if value is None else f"{value:.3f}"
            bound = f"{check['bound']:.2f}"
        else:
            text = "missing" if value is None else f"{value:.4%}"
            bound = f"{check['bound']:.2%}"
        print(f"{check['gate']:76s} {text:>9s} {bound:>8s}  "
              f"{'PASS' if check['passed'] else 'FAIL'}")
    print("\nreported (median / p95):")
    for name, metrics in reported.items():
        for metric, stats in metrics.items():
            median, p95 = stats.get("median"), stats.get("p95")
            if median is None:
                continue
            print(f"  {name:11s} {metric:45s} {median:8.3%} / {p95:8.3%}"
                  f"  (n={stats.get('count')})")

    if arrival_reported:
        print("\narrival, reported (median / max):")
        for row in arrival_reported:
            for metric in ARRIVAL_REPORTED:
                stats = row.get(metric) or {}
                if "median" in stats:
                    print(f"  {row['group']:34s} {metric:30s} "
                          f"{stats['median']:8.3%} / {stats['max']:8.3%}"
                          f"  (n={stats['count']})")

    passed = all(check["passed"] for check in checks)
    args.output.write_text(json.dumps(
        {"passed": passed, "checks": checks, "reported": reported,
         "arrival_reported": arrival_reported},
        indent=2) + "\n")
    print(f"\n{'ALL REQUIRED GATES PASSED' if passed else 'GATES FAILED'}"
          f" -- {args.output}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
