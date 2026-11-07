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

"""Acceptance gates for a source-escape bake: one table, one JSON, one verdict.

Reads the reports the bake's scoring steps wrote, applies every bound, runs
the C++/Python lockstep guards against the build that produced the payload,
prints a table and writes `gates.json`.  Exit status 1 when any gate fails.

Each bound carries its origin.  Where a documented bar exists it is used;
otherwise the bound is the last accepted payload's value with headroom,
max(1.25 x last, last + 0.5 percentage points), so a statistically equivalent
re-bake passes and a real regression does not.

Photonuclear and the GEANT4 leg of the acceptance grid are optional parts of
the bake: without --photonuclear its gates are skipped, and --without-geant4
skips the acceptance grid's GEANT4 comparison.  Both are explicit, so a run
that meant to include them and lost their data fails instead.

Usage:
    check_source_escape_gates.py --offgrid J --deep-tau-tail J
        --multiple-scatter J [--photonuclear J] --kernel-span J
        --acceptance J [--without-geant4] --build DIR --output gates.json
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def headroom(last: float) -> float:
    """max(1.25 x last, last + 0.5 percentage points) for a fraction."""
    return max(1.25 * last, last + 0.005)


# (report, key path, bound, origin); every value is a fraction unless noted.
BOUNDS = [
    ("offgrid", ("deep_tau_brems", "median"), 0.10,
     "bar: single digits; last accepted payload: 6.8%"),
    ("offgrid", ("leakage", "median"), headroom(0.0048),
     "last accepted payload: 0.48%"),
    ("offgrid", ("fluorescence", "median"), headroom(0.0305),
     "last accepted payload: 3.05%"),
    ("deep_tau_tail", ("worst_spurious_yield", 0), 1.25 * 9.61e-8,
     "photons/primary; last accepted payload: 9.61e-8 (x1.25)"),
    ("multiple_scatter", ("test", "yield_error_median_abs"),
     headroom(0.0126), "last accepted payload: 1.26%"),
    ("multiple_scatter", ("test", "yield_error_p95_abs"),
     headroom(0.084), "last accepted payload: 8.4%"),
    ("multiple_scatter", ("test", "shape_tv_median"),
     headroom(0.0497), "last accepted payload: 4.97%"),
    ("multiple_scatter", ("test", "shape_tv_p95"),
     headroom(0.1146), "last accepted payload: 11.46%"),
    ("photonuclear", ("system_impact", "heldout",
                      "max_added_percentage_points"), 2.0,
     "percentage points of the escaping spectrum; bar: <= 2 pp added "
     "to the CeeLo/Geant4 discrepancy"),
    ("photonuclear", ("system_impact", "offgrid",
                      "max_added_percentage_points"), 2.0,
     "percentage points of the escaping spectrum; bar: <= 2 pp added"),
    ("acceptance", ("vs_ceelo", "abs_median"), headroom(0.0083),
     "last accepted payload: 0.83%"),
    ("acceptance", ("vs_ceelo", "p95"), headroom(0.1049),
     "last accepted payload: 10.49%"),
    ("acceptance", ("vs_geant4", "abs_median"), headroom(0.0098),
     "last accepted payload: 0.98%"),
    ("acceptance", ("vs_geant4", "p95"), headroom(0.0843),
     "last accepted payload: 8.43%"),
]
# Reported next to the gates but not gated: known to sit above any bar the
# current model can meet.
REPORTED = [
    ("offgrid", ("deep_tau_brems", "p95"),
     "known unmet; last accepted payload: 51.5%"),
    ("offgrid", ("leakage", "p95"),
     "a multiple-scatter residual; last accepted payload: 5.85%"),
    ("offgrid", ("fluorescence", "p95"),
     "fluorescence atoms between training nodes; last accepted payload: 35.5%"),
    ("offgrid", ("fluorescence_continuum_share", "median"),
     "Compton-scattered shield x-rays (not modelled) / the atoms"),
    ("photonuclear", ("metrics", "heldout", "yield_median"),
     "held-out element; last accepted payload: 27.2%"),
    ("photonuclear", ("metrics", "offgrid", "yield_median"),
     "off-grid; last accepted payload: 13.4%"),
]
LOCKSTEP = ("check_base_lockstep.py", "check_feature_lockstep.py",
            "check_basis_lockstep.py")


def lookup(report: dict | None, path: tuple):
    value = report
    for key in path:
        if value is None:
            return None
        try:
            value = value[key]
        except (KeyError, IndexError, TypeError):
            return None
    return value


def show(value, bound_is_pp: bool) -> str:
    if value is None:
        return "missing"
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        return str(value)
    if bound_is_pp:
        return f"{value:.3f}"
    if abs(value) < 1e-3:
        return f"{value:.3e}"
    return f"{100 * value:.2f}%"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("offgrid", "deep-tau-tail", "multiple-scatter",
                 "kernel-span", "acceptance"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument(
        "--photonuclear", type=Path,
        help="photonuclear fit report; omit when photonuclear is not baked")
    parser.add_argument(
        "--without-geant4", action="store_true",
        help="the acceptance grid has no GEANT4 leg: skip its GEANT4 bounds")
    parser.add_argument("--build", type=Path, required=True,
                        help="build tree for the lockstep guards")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    reports = {}
    for name in ("offgrid", "deep_tau_tail", "multiple_scatter",
                 "photonuclear", "kernel_span", "acceptance"):
        path = getattr(args, name)
        reports[name] = (json.loads(path.read_text())
                         if path is not None and path.exists() else None)

    def applies(report: str, path: tuple) -> bool:
        if report == "photonuclear":
            return args.photonuclear is not None
        if report == "acceptance" and path[0] == "vs_geant4":
            return not args.without_geant4
        return True

    gates = []
    for report, path, bound, origin in BOUNDS:
        if not applies(report, path):
            continue
        value = lookup(reports[report], path)
        passed = (isinstance(value, (int, float)) and math.isfinite(value)
                  and abs(value) <= bound)
        gates.append({"gate": f"{report}." + ".".join(map(str, path)),
                      "value": value, "bound": bound, "origin": origin,
                      "pass": passed})
    kernel = reports["kernel_span"] or {}
    gates.append({
        "gate": "kernel_span.pass", "value": kernel.get("pass"),
        "bound": True, "pass": kernel.get("pass") is True,
        "origin": "every component within its span/slope bound "
                  "(check_kernel_span.py)"})
    for script in LOCKSTEP:
        command = [sys.executable, str(HERE / script)]
        if script != "check_basis_lockstep.py":
            command += ["--build", str(args.build)]
        completed = subprocess.run(command, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        tail = completed.stdout.strip().splitlines()[-1:] or [""]
        gates.append({
            "gate": f"lockstep.{script.removesuffix('.py')}",
            "value": completed.returncode, "bound": 0,
            "pass": completed.returncode == 0,
            "origin": f"C++/Python parity guard; {tail[0][:80]}"})

    width = max(len(g["gate"]) for g in gates)
    print(f"{'gate':{width}s} {'value':>11s} {'bound':>11s}  verdict  origin")
    for gate in gates:
        pp = "percentage points" in gate["origin"]
        bound = gate["bound"]
        bound_text = (show(bound, pp) if isinstance(bound, float)
                      else str(bound))
        print(f"{gate['gate']:{width}s} {show(gate['value'], pp):>11s} "
              f"{bound_text:>11s}  {'PASS' if gate['pass'] else 'FAIL':7s}"
              f"  {gate['origin']}")
    reported = []
    print("\nreported, not gated:")
    for report, path, note in REPORTED:
        if not applies(report, path):
            continue
        value = lookup(reports[report], path)
        reported.append({"metric": f"{report}." + ".".join(map(str, path)),
                         "value": value, "note": note})
        print(f"  {reported[-1]['metric']:{width - 2}s} "
              f"{show(value, False):>11s}  {note}")
    failed = [gate["gate"] for gate in gates if not gate["pass"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({
        "pass": not failed, "failed": failed, "gates": gates,
        "reported": reported,
        "inputs": {name: str(getattr(args, name)) for name in (
            "offgrid", "deep_tau_tail", "multiple_scatter", "photonuclear",
            "kernel_span", "acceptance")},
        "photonuclear_baked": args.photonuclear is not None,
        "geant4_acceptance_leg": not args.without_geant4,
    }, indent=2) + "\n")
    print(f"\n{'ALL GATES PASS' if not failed else 'FAILED: ' + ', '.join(failed)}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
