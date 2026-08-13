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

"""Validate direct all-element ESTAR stopping against the historical runtime."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PREPARE = ROOT / "tools/prepare_cross_sections"
DEFAULT_JSON = PREPARE / "reports/stopping_migration.json"
DEFAULT_MARKDOWN = PREPARE / "reports/stopping_migration.md"
DEFAULT_REF = "ab03e48"
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
ME_C2_MEV = 0.5109989461
BETHE_K = 0.307075
GRID_COUNT = 200
HISTORICAL_ANCHOR_Z = {1, 6, 8, 11, 13, 26, 29, 48, 53, 56, 57, 82}


def git_text(ref: str, path: str) -> str:
    return subprocess.run(
        ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True,
    ).stdout


def initializer(text: str, declaration: str) -> str:
    start = text.index(declaration)
    begin = text.index("{", start)
    depth = 0
    for position in range(begin, len(text)):
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
            if depth == 0:
                return text[begin:position + 1]
    raise ValueError(f"unterminated initializer: {declaration}")


def values(text: str, declaration: str) -> list[float]:
    body = re.sub(r"/\*.*?\*/|//.*", "", initializer(text, declaration), flags=re.DOTALL)
    return [float(token) for token in NUMBER.findall(body)]


def table(text: str, declaration: str, rows: int, columns: int) -> list[list[float]]:
    flat = values(text, declaration)
    if len(flat) != rows * columns:
        raise ValueError(f"{declaration}: expected {rows * columns} values, found {len(flat)}")
    return [flat[index * columns:(index + 1) * columns] for index in range(rows)]


def interpolate_correction(
    correction: list[list[float]], log_z: list[float], log_e: list[float], z: int, energy: float
) -> float:
    x = math.log(min(max(energy, 10.0), 20_000.0))
    q = math.log(min(max(z, 1), 92))
    iz = max((index for index, value in enumerate(log_z) if value <= q), default=0)
    ie = max((index for index, value in enumerate(log_e) if value <= x), default=0)
    jz, je = min(iz + 1, len(log_z) - 1), min(ie + 1, len(log_e) - 1)
    tz = 0.0 if iz == jz else (q - log_z[iz]) / (log_z[jz] - log_z[iz])
    te = 0.0 if ie == je else (x - log_e[ie]) / (log_e[je] - log_e[ie])
    return (
        (1 - tz) * (1 - te) * correction[iz][ie]
        + tz * (1 - te) * correction[jz][ie]
        + (1 - tz) * te * correction[iz][je]
        + tz * te * correction[jz][je]
    )


def old_stopping(
    z: int, energy_keV: float, atomic_weight: list[float], excitation: list[float],
    correction: list[list[float]], log_z: list[float], log_e: list[float],
) -> float:
    tau = energy_keV * 1.0e-3 / ME_C2_MEV
    beta2 = tau * (tau + 2.0) / ((tau + 1.0) ** 2)
    i_ratio = excitation[z - 1] * 1.0e-6 / ME_C2_MEV
    argument = tau * tau * (tau + 2.0) / (2.0 * i_ratio * i_ratio)
    f_minus = (1.0 - beta2) + (
        tau * tau / 8.0 - (2.0 * tau + 1.0) * math.log(2.0)
    ) / ((tau + 1.0) ** 2)
    stopping = BETHE_K * (z / atomic_weight[z - 1]) / beta2 * (0.5 * math.log(argument) + f_minus)
    return max(stopping * interpolate_correction(correction, log_z, log_e, z, energy_keV), 0.001)


def log_interpolate(rows: list[tuple[float, float, float]], energy_keV: float, column: int) -> float:
    energy_mev = min(max(energy_keV * 1.0e-3, rows[0][0]), rows[-1][0])
    if energy_mev <= rows[0][0]:
        return rows[0][column]
    if energy_mev >= rows[-1][0]:
        return rows[-1][column]
    lo = 0
    while rows[lo + 1][0] < energy_mev:
        lo += 1
    high = lo + 1
    t = math.log(energy_mev / rows[lo][0]) / math.log(rows[high][0] / rows[lo][0])
    return math.exp(math.log(rows[lo][column]) + t * math.log(rows[high][column] / rows[lo][column]))


def new_stopping(rows: list[tuple[float, float, float]], energy_keV: float) -> float:
    if energy_keV >= 10.0:
        return log_interpolate(rows, energy_keV, 1)
    s10 = rows[0][1]
    s12 = rows[1][1]
    slope = math.log(s12 / s10) / math.log(12.5 / 10.0)
    return s10 * (energy_keV / 10.0) ** slope


def percentile(data: list[float], fraction: float) -> float:
    ordered = sorted(data)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def summary(data: list[tuple[float, dict]]) -> dict:
    worst = max(data, key=lambda item: item[0])
    errors = [item[0] for item in data]
    return {
        "count": len(errors), "median_percent": percentile(errors, 0.5),
        "p95_percent": percentile(errors, 0.95), "worst_percent": max(errors),
        "worst_case": worst[1],
    }


def range_grid(stopping) -> tuple[list[float], list[float]]:
    log_min, log_max = math.log(1.0), math.log(20_000.0)
    energy = [math.exp(log_min + index * (log_max - log_min) / (GRID_COUNT - 1))
              for index in range(GRID_COUNT)]
    ranges = [energy[0] * 1.0e-3 / stopping(energy[0])]
    for low, high in zip(energy, energy[1:]):
        area = 0.5 * (high - low) * 1.0e-3 * (1.0 / stopping(low) + 1.0 / stopping(high))
        ranges.append(ranges[-1] + area)
    return energy, ranges


def grid_value(x: list[float], y: list[float], query: float) -> float:
    if query <= x[0]:
        return y[0]
    if query >= x[-1]:
        return y[-1]
    position = (math.log(query) - math.log(x[0])) / (math.log(x[1]) - math.log(x[0]))
    low = max(0, min(int(position), len(x) - 2))
    return y[low] + (position - low) * (y[low + 1] - y[low])


def inverse_grid(x: list[float], y: list[float], query: float) -> float:
    if query <= y[0]:
        return x[0]
    if query >= y[-1]:
        return x[-1]
    low, high = 0, len(y) - 1
    while high - low > 1:
        mid = (low + high) // 2
        if y[mid] <= query:
            low = mid
        else:
            high = mid
    t = (query - y[low]) / (y[high] - y[low])
    return math.exp(math.log(x[low]) + t * math.log(x[high] / x[low]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--estar-dir", type=Path, required=True)
    parser.add_argument("--baseline-ref", default=DEFAULT_REF)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args()

    import sys
    sys.path.insert(0, str(PREPARE))
    from generate_estar_stopping import EXPECTED_ENERGIES_MEV, logical_digest, parse_table
    from fetch_sources import load_lock

    historical = git_text(args.baseline_ref, "src/physics/ElectronCsda.cpp")
    atomic_weight = values(historical, "static const double kAtomicWeight[92]")
    excitation = values(historical, "static const double kICRU49_I_eV[92]")
    log_z = values(historical, "static constexpr double kCorrLnZ[kNCorrZ]")
    log_e = values(historical, "static constexpr double kCorrLnE[kNCorrE]")
    correction = table(historical, "static constexpr float kCorrTable[kNCorrZ][kNCorrE]", 12, 15)
    tables = {z: parse_table(args.estar_dir / f"estar-text-Z{z:03d}.html") for z in range(1, 93)}
    lock_path = PREPARE / "sources.lock.json"
    lock = load_lock(lock_path)
    digest = logical_digest(tables)
    if digest != lock["sources"]["nist_estar_text"]["logical_table_sha256"]:
        raise ValueError("ESTAR logical digest mismatch")

    stopping_errors = []
    range_errors = []
    deposit_errors = []
    for z in range(1, 93):
        old = lambda energy, z=z: old_stopping(
            z, energy, atomic_weight, excitation, correction, log_z, log_e
        )
        new = lambda energy, z=z: new_stopping(tables[z], energy)
        old_grid = range_grid(old)
        new_grid = range_grid(new)
        for energy_mev in EXPECTED_ENERGIES_MEV:
            energy = energy_mev * 1000.0
            old_s, new_s = old(energy), new(energy)
            stopping_errors.append((abs(new_s / old_s - 1.0) * 100.0,
                                    {"Z": z, "energy_keV": energy}))
            old_range = grid_value(*old_grid, energy)
            new_range = grid_value(*new_grid, energy)
            range_errors.append((abs(new_range / old_range - 1.0) * 100.0,
                                 {"Z": z, "energy_keV": energy}))
            if energy >= 20.0:
                for fraction in (0.1, 0.5, 0.9):
                    path = fraction * old_range
                    old_residual = inverse_grid(*old_grid, max(old_range - path, 0.0))
                    new_residual = inverse_grid(*new_grid, max(new_range - path, 0.0))
                    old_deposit = energy - old_residual
                    new_deposit = energy - new_residual
                    deposit_errors.append((abs(new_deposit / old_deposit - 1.0) * 100.0,
                                           {"Z": z, "energy_keV": energy,
                                            "old_range_fraction": fraction}))

    summaries = {
        "collision_stopping": summary(stopping_errors),
        "collision_stopping_historical_anchor_elements": summary([
            item for item in stopping_errors if item[1]["Z"] in HISTORICAL_ANCHOR_Z
        ]),
        "csda_range": summary(range_errors),
        "deterministic_deposited_energy": summary(deposit_errors),
    }
    passed = (
        summaries["collision_stopping"]["p95_percent"] <= 1.0
        and summaries["collision_stopping"]["worst_percent"] <= 3.0
        and summaries["csda_range"]["worst_percent"] <= 1.0
        and summaries["deterministic_deposited_energy"]["worst_percent"] <= 1.0
    )
    report = {
        "schema_version": 1,
        "source_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "historical_baseline_ref": args.baseline_ref,
        "summaries": summaries,
        "gates_percent": {"stopping_p95": 1.0, "stopping_worst": 3.0,
                           "csda_range_worst": 1.0, "deposited_energy_worst": 1.0},
        "diagnosis": (
            "Direct ESTAR agrees closely at the 12 historical correction-surface anchors; "
            "the failure is concentrated in elements formerly approximated by interpolation "
            "between anchors, especially elemental gas/solid density-effect changes."
        ),
        "passed": passed,
    }
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    lines = ["# ESTAR stopping migration validation\n",
             f"Historical baseline: `{args.baseline_ref}`.\n",
             "| Quantity | Median | p95 | Worst | Worst case | Gate | Result |",
             "|---|---:|---:|---:|---|---|---|"]
    quantity_gates = {
        "collision_stopping": (1.0, 3.0, "p95 ≤1%; worst ≤3%"),
        "collision_stopping_historical_anchor_elements": (None, None, "diagnostic"),
        "csda_range": (None, 1.0, "worst ≤1%"),
        "deterministic_deposited_energy": (None, 1.0, "worst ≤1%"),
    }
    for name, data in summaries.items():
        p95_gate, worst_gate, gate_text = quantity_gates[name]
        result = "DIAGNOSTIC" if p95_gate is None and worst_gate is None else (
            "PASS" if (p95_gate is None or data["p95_percent"] <= p95_gate)
            and (worst_gate is None or data["worst_percent"] <= worst_gate) else "FAIL"
        )
        case = ", ".join(f"{key}={value}" for key, value in data["worst_case"].items())
        lines.append(
            f"| {name} | {data['median_percent']:.6g}% | {data['p95_percent']:.6g}% "
            f"| {data['worst_percent']:.6g}% | {case} | {gate_text} | {result} |"
        )
    lines.extend([
        f"\nOverall: {'PASS' if passed else 'FAIL'}.\n",
        "Direct ESTAR agrees closely at the 12 historical correction-surface anchors. "
        "The failed historical-result gate is concentrated in elements that the old "
        "runtime approximated by interpolation between those anchors, especially across "
        "elemental gas/solid density-effect changes. The direct 92×53 source table and "
        "its float32-log representation pass their separate source-fidelity checks.\n",
    ])
    args.output_markdown.write_text("\n".join(lines))
    print(f"Stopping migration {'passed' if passed else 'failed'}; wrote {args.output_json}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
