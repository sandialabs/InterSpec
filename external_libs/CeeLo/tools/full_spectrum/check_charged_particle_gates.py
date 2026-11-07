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

"""Acceptance gates for a baked charged-particle source payload.

Reads the joint-shape payload, its score report and its C++/Python parity
record (fit_charged_particle_joint_shape.py score / parity), prints one
table and writes the gates JSON.  Exit status 0 means every gate passed; 1
means at least one failed.

Gates:
  format      header `CEELO_CHARGED_PARTICLE_SOURCE 1`; every component a
              spline field; shape rank <= 16 (0 for the yield-only
              components); every grid dimension in [4, 32]; archetypes
              non-negative with rows summing to 1 within 1e-12; all values
              finite
  consistency the report and the parity record were computed for this
              payload (SHA-256)
  parity      C++ evaluator vs the Python reference evaluator on whole-Z
              holdout cells: max relative difference <= 5e-3 (both use the
              smooth effective-Z cross-section blend; the residual is the two
              implementations of the absorption envelope)
  accuracy    |yield error| and shape total variation, median and p95, on
              the whole-Z holdout planes (the matched holdout), the
              total-spectrum shape TV on the two off-grid banks, and the
              positron components on the off-grid final test

Each accuracy bound is the last accepted payload's value with headroom
max(1.25 x last, last + 0.5 percentage points), rounded up to 0.01 pp; the
comment on each bound gives the last accepted value.  The off-grid electron
yield errors are reported, not gated.  A component is scored in a cell when
its truth yield is at least fit_charged_particle_source.MINIMUM_SCORED_YIELD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_charged_particle_source as fit  # noqa: E402
import fit_charged_particle_source_structured as structured  # noqa: E402

PARITY_TOLERANCE = 5e-3  # last accepted payload: 2.8e-3
MAXIMUM_RANK = 16
GRID_RANGE = (4, 32)
ROW_SUM_TOLERANCE = 1e-12

YIELD = "yield_relative_error"
TV = "shape_total_variation"
# (bank, component, metric, median bound, p95 bound)
ACCURACY = (
    # off-grid electron banks, total spectrum
    ("offgrid_validation", "electron_total", TV, 0.0218, 0.0843),  # last 1.675% / 6.737%
    ("offgrid_final_test", "electron_total", TV, 0.0171, 0.0700),  # last 1.200% / 5.593%
    # whole-Z holdout planes Z = 35, 64
    ("test", "electron_brems", YIELD, 0.0166, 0.0990),  # last 1.151% / 7.913%
    ("test", "electron_brems", TV, 0.0323, 0.1883),  # last 2.576% / 15.064%
    ("test", "positron_brems", YIELD, 0.0207, 0.1774),  # last 1.562% / 14.187%
    ("test", "positron_brems", TV, 0.0297, 0.1398),  # last 2.375% / 11.182%
    ("test", "positron_containment", YIELD, 0.0052, 0.0864),  # last 0.010% / 6.911%
    ("test", "positron_511_atom", YIELD, 0.0116, 0.1085),  # last 0.659% / 8.679%
    ("test", "positron_annihilation_continuum", YIELD, 0.0224, 0.1933),  # last 1.733% / 15.459%
    ("test", "positron_annihilation_continuum", TV, 0.0353, 0.2156),  # last 2.817% / 17.243%
    ("test", "positron_annihilation_in_flight_direct", YIELD, 0.0290, 0.2053),  # last 2.320% / 16.416%
    ("test", "positron_annihilation_in_flight_direct", TV, 0.0357, 0.2018),  # last 2.852% / 16.141%
    ("test", "positron_annihilation_in_flight_continuum", YIELD, 0.0538, 0.3184),  # last 4.303% / 25.470%
    ("test", "positron_annihilation_in_flight_continuum", TV, 0.1198, 0.3980),  # last 9.580% / 31.832%
    ("test", "positron_total", YIELD, 0.0103, 0.0679),  # last 0.528% / 5.427%
    ("test", "positron_total", TV, 0.0070, 0.0235),  # last 0.197% / 1.841%
    # off-grid final test, positrons (T, Z 16/47/79, r off the training grid;
    # last accepted payload scored on this bank, its first with positrons)
    ("offgrid_final_test", "positron_brems", YIELD, 0.0218, 0.1705),  # last 1.674% / 13.632%
    ("offgrid_final_test", "positron_brems", TV, 0.0191, 0.0820),  # last 1.409% / 6.556%
    ("offgrid_final_test", "positron_containment", YIELD, 0.0051, 0.2436),  # last 0.000% / 19.484%
    ("offgrid_final_test", "positron_511_atom", YIELD, 0.0094, 0.2501),  # last 0.436% / 20.003%
    ("offgrid_final_test", "positron_annihilation_continuum", YIELD, 0.0187, 0.3144),  # last 1.361% / 25.150%
    ("offgrid_final_test", "positron_annihilation_continuum", TV, 0.0137, 0.0505),  # last 0.862% / 4.036%
    ("offgrid_final_test", "positron_annihilation_in_flight_direct", YIELD, 0.0188, 0.0912),  # last 1.377% / 7.294%
    ("offgrid_final_test", "positron_annihilation_in_flight_direct", TV, 0.0155, 0.0404),  # last 1.043% / 3.230%
    ("offgrid_final_test", "positron_annihilation_in_flight_continuum", YIELD, 0.0452, 0.1493),  # last 3.615% / 11.937%
    ("offgrid_final_test", "positron_annihilation_in_flight_continuum", TV, 0.0401, 0.0921),  # last 3.206% / 7.365%
    ("offgrid_final_test", "positron_total", YIELD, 0.0105, 0.1445),  # last 0.541% / 11.553%
    ("offgrid_final_test", "positron_total", TV, 0.0077, 0.0350),  # last 0.267% / 2.795%
    # shield K x-rays, set by the first bake with them (holdout planes and the
    # off-grid final test; the p95 tails are low-count cells)
    ("test", "electron_fluorescence_k_alpha", YIELD, 0.0930, 0.6901),  # first bake: 7.432% / 55.203% (n=661)
    ("test", "electron_fluorescence_k_beta", YIELD, 0.0688, 0.3355),  # first bake: 5.498% / 26.835% (n=553)
    ("test", "positron_fluorescence_k_alpha", YIELD, 0.0711, 0.6809),  # first bake: 5.684% / 54.466% (n=316)
    ("test", "positron_fluorescence_k_beta", YIELD, 0.0585, 0.1932),  # first bake: 4.678% / 15.453% (n=265)
    ("offgrid_final_test", "electron_fluorescence_k_alpha", YIELD, 0.0468, 0.8318),  # first bake: 3.741% / 66.540% (n=117)
    ("offgrid_final_test", "electron_fluorescence_k_beta", YIELD, 0.0347, 0.1569),  # first bake: 2.773% / 12.546% (n=84)
    ("offgrid_final_test", "positron_fluorescence_k_alpha", YIELD, 0.0658, 0.5632),  # first bake: 5.259% / 45.050% (n=152)
    ("offgrid_final_test", "positron_fluorescence_k_beta", YIELD, 0.0608, 0.6067),  # first bake: 4.857% / 48.529% (n=109)
)
# (bank, component, metric), printed for information
REPORTED = (
    ("offgrid_validation", "electron_total", YIELD),
    ("offgrid_final_test", "electron_total", YIELD),
    ("test", "electron_total", "yield_relative_uncertainty"),
    ("test", "positron_total", "yield_relative_uncertainty"),
    # Shield L x-rays: 5-6 off-grid cells per particle and none on the holdout
    # planes, too few to bound.
    *((bank, f"{particle}_fluorescence_l", YIELD)
      for bank in ("test", "offgrid_final_test")
      for particle in ("electron", "positron")),
)


def format_problems(path: Path) -> list[str]:
    """Every violated format invariant, as text; empty when the payload is
    well formed."""
    lines = path.read_text().splitlines()
    if not lines or lines[0].split() != ["CEELO_CHARGED_PARTICLE_SOURCE", "1"]:
        return ["header is not 'CEELO_CHARGED_PARTICLE_SOURCE 1'"]
    for line in lines:
        for token in line.split()[1:]:
            try:
                value = float(token)
            except ValueError:
                continue
            if not math.isfinite(value):
                return [f"non-finite value in record {line.split()[0]}"]
    try:
        payload = structured.load_payload(path)
    except (AssertionError, IndexError, ValueError) as error:
        return [f"does not parse: {error!r}"]
    problems = []
    for name in fit.MODEL_COMPONENTS:
        component = payload.components.get(name)
        if component is None:
            problems.append(f"{name}: no spline field")
            continue
        grids = [component.yield_grid]
        if component.rank >= 2:
            grids.append(component.weight_grid)
        if any(not GRID_RANGE[0] <= size <= GRID_RANGE[1]
               for grid in grids for size in grid):
            problems.append(f"{name}: grid dimension outside {GRID_RANGE}")
        expected_zero = name in structured.YIELD_ONLY_COMPONENTS
        if component.rank > MAXIMUM_RANK or (
                expected_zero != (component.rank == 0)):
            problems.append(f"{name}: shape rank {component.rank}")
        if component.rank >= 1:
            archetypes = component.archetypes
            if (archetypes < 0.0).any():
                problems.append(f"{name}: negative archetype value")
            worst = float(abs(archetypes.sum(axis=1) - 1.0).max())
            if worst > ROW_SUM_TOLERANCE:
                problems.append(
                    f"{name}: archetype row sums off by {worst:.2e}")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True,
                        help="report JSON of the score subcommand")
    parser.add_argument("--parity", type=Path, required=True,
                        help="JSON of the parity subcommand")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    digest = hashlib.sha256(args.payload.read_bytes()).hexdigest()
    report = json.loads(args.report.read_text())
    parity = json.loads(args.parity.read_text())
    checks = []

    problems = format_problems(args.payload)
    checks.append({"gate": "payload format invariants",
                   "value": None, "bound": None,
                   "detail": "; ".join(problems) or "ok",
                   "passed": not problems})
    for label, record in (("score report", report), ("parity record", parity)):
        same = record.get("payload_sha256") == digest
        checks.append({
            "gate": f"{label} is for this payload", "value": None,
            "bound": None, "passed": same,
            "detail": "ok" if same else
            "stale: re-run the score and parity steps for this payload"})
    worst = parity.get("worst_relative_difference")
    checks.append({
        "gate": "C++ vs Python parity (max relative difference)",
        "value": worst, "bound": PARITY_TOLERANCE,
        "detail": json.dumps(parity.get("worst_at", {})),
        "passed": worst is not None and worst <= PARITY_TOLERANCE})
    banks = report.get("banks", {})
    for bank, component, metric, median_bound, p95_bound in ACCURACY:
        stats = (banks.get(bank, {}).get("components", {})
                 .get(component, {}).get(metric, {}))
        for statistic, bound in (("median", median_bound),
                                 ("p95", p95_bound)):
            value = stats.get(statistic)
            checks.append({
                "gate": f"{bank} {component} {metric} {statistic}",
                "value": value, "bound": bound, "count": stats.get("count"),
                "passed": value is not None and value <= bound})
    reported = {}
    for bank, component, metric in REPORTED:
        stats = (banks.get(bank, {}).get("components", {})
                 .get(component, {}).get(metric, {}))
        reported[f"{bank} {component} {metric}"] = {
            key: stats.get(key) for key in ("count", "median", "p95")}

    width = max(len(check["gate"]) for check in checks)
    print(f"{'gate':{width}s} {'value':>9s} {'bound':>8s}  result")
    for check in checks:
        value, bound = check["value"], check["bound"]
        if bound is None:
            text, bound_text = check["detail"][:9], ""
        elif bound < 0.001:
            text = "missing" if value is None else f"{value:.2e}"
            bound_text = f"{bound:.0e}"
        else:
            text = "missing" if value is None else f"{value:.3%}"
            bound_text = f"{bound:.2%}"
        print(f"{check['gate']:{width}s} {text:>9s} {bound_text:>8s}  "
              f"{'PASS' if check['passed'] else 'FAIL'}")
        if not check["passed"] and check.get("detail") not in (None, "ok"):
            print(f"    {check['detail']}")
    print("\nreported (median / p95):")
    for key, stats in reported.items():
        if stats.get("median") is not None:
            print(f"  {key:78s} {stats['median']:8.3%} / {stats['p95']:8.3%}"
                  f"  (n={stats['count']})")

    passed = all(check["passed"] for check in checks)
    args.output.write_text(json.dumps(
        {"passed": passed, "payload_sha256": digest, "checks": checks,
         "reported": reported}, indent=2) + "\n")
    print(f"\n{'ALL GATES PASSED' if passed else 'GATES FAILED'} -- "
          f"{args.output}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
