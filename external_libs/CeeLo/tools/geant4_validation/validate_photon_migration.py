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

"""Compare direct EPICS2023 runtime curves with the historical main results.

The historical cache is read transiently from a pinned git object. The report
contains aggregate errors and case locations only; it never reproduces raw G4
cross-section columns.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import io
import json
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CPP = ROOT / "src/cross_sections/photon_epics_data.cpp"
DEFAULT_JSON = ROOT / "tools/prepare_cross_sections/reports/photon_migration.json"
DEFAULT_MARKDOWN = ROOT / "tools/prepare_cross_sections/reports/photon_migration.md"
DEFAULT_BASELINE = "ab03e48"
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
PROCESS = {
    "photoelectric": "pe_barn",
    "compton": "compt_barn",
    "rayleigh": "rayl_barn",
    "pair_production": "conv_barn",
}
GATES = {
    "compton": (0.25, 0.30),
    "pair_production": (0.05, 0.10),
    "photoelectric": (1.25, 4.00),
    "rayleigh": (0.75, 6.50),
}
EDGE_PAIR_RELATIVE_SPACING = 5.0e-4
MATERIALS = {
    "Ge": [(32, 1.0)],
    "NaI": [(11, 0.1533750517), (53, 0.8466249483)],
    "CsI": [(55, 0.5115488686), (53, 0.4884511314)],
    "CZT": [(48, 0.4299534318), (30, 0.0277845537), (52, 0.5422620145)],
    "Al": [(13, 1.0)],
    "SS304": [(26, 0.70), (24, 0.19), (28, 0.10), (25, 0.01)],
    "Pb": [(82, 1.0)],
}
ATOMIC_WEIGHT = {
    11: 22.989769, 13: 26.981538, 24: 51.9961, 25: 54.938044,
    26: 55.845, 28: 58.6934, 30: 65.38, 32: 72.630, 48: 112.414,
    52: 127.60, 53: 126.90447, 55: 132.90545, 82: 207.2,
}


def git_text(ref: str, path: str) -> str:
    return subprocess.run(
        ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True,
    ).stdout


def parse_array(text: str, name: str) -> list[float]:
    match = re.search(
        rf"(?:static|extern) const (?:float|uint16_t) {re.escape(name)}\[[^]]+\]\s*=\s*\{{(.*?)\}};",
        text, re.DOTALL,
    )
    if not match:
        raise ValueError(f"missing generated array {name}")
    return [float(token) for token in NUMBER.findall(re.sub(r"//.*", "", match.group(1)))]


def direct_curves(text: str) -> dict[tuple[int, str], tuple[list[float], list[float]]]:
    rayleigh_energy = parse_array(text, "g_rayleigh_log_energy")
    rayleigh_grid_offset = [
        int(value) for value in parse_array(text, "g_rayleigh_group_grid_offset")
    ]
    rayleigh_value_offset = [
        int(value) for value in parse_array(text, "g_rayleigh_group_value_offset")
    ]
    rayleigh_packed = [
        int(value) for value in parse_array(text, "g_rayleigh_log_value_q")
    ]
    energy_pool = parse_array(text, "g_photon_energy_pool")
    grid_index = [int(value) for value in parse_array(text, "g_photon_process_grid_index")]
    packed_value = [int(value) for value in parse_array(text, "g_photon_process_log_value_q")]
    descriptor_match = re.search(
        r"g_photon_epics_data\[92\]\s*=\s*\{(.*?)\n\};",
        text, re.DOTALL,
    )
    if not descriptor_match:
        raise ValueError("missing generated photon descriptor table")
    descriptor_lines = [
        line for line in descriptor_match.group(1).splitlines() if "// Z=" in line
    ]
    if len(descriptor_lines) != 92:
        raise ValueError("generated photon descriptor table is not 92 elements")
    process_slot = {
        "compton": 0,
        "pair_production": 1,
        "photoelectric": 2,
    }
    curves = {}
    for z in range(1, 93):
        fields = [
            float(token) for token in NUMBER.findall(
                re.sub(r"//.*", "", descriptor_lines[z - 1])
            )
        ]
        if len(fields) != 22:
            raise ValueError(f"Z={z}: malformed compact descriptor")
        group = (z - 1) // 4
        lane = (z - 1) % 4
        grid_begin, grid_end = (
            rayleigh_grid_offset[group], rayleigh_grid_offset[group + 1]
        )
        size = grid_end - grid_begin
        value_begin = rayleigh_value_offset[group] + lane * size
        offset, scale = fields[0], fields[1]
        x = rayleigh_energy[grid_begin:grid_end]
        packed = rayleigh_packed[value_begin:value_begin + size]
        y = [-30.0 if code == 0 else offset + (code - 1) * scale for code in packed]
        curves[z, "rayleigh"] = (x, y)
        for process in PROCESS:
            if process == "rayleigh":
                continue
            base = 2 + process_slot[process] * 4
            size = int(fields[base])
            data_offset = int(fields[base + 1])
            offset, scale = fields[base + 2], fields[base + 3]
            x = [
                energy_pool[grid_index[data_offset + index]]
                for index in range(size)
            ]
            packed = packed_value[data_offset:data_offset + size]
            y = [-30.0 if code == 0 else offset + (code - 1) * scale for code in packed]
            if len(x) != len(y) or len(x) < 2:
                raise ValueError(f"Z={z} {process}: malformed direct curve")
            curves[z, process] = (x, y)
    return curves


def interpolate(curve: tuple[list[float], list[float]], energy_keV: float) -> float:
    x, y = curve
    query = math.log10(energy_keV * 1.0e-3)
    if query <= x[0]:
        return 0.0 if y[0] <= -30.0 else 10.0**y[0]
    if query >= x[-1]:
        return 0.0 if y[-1] <= -30.0 else 10.0**y[-1]
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x[mid] <= query:
            lo = mid
        else:
            hi = mid
    fraction = (query - x[lo]) / (x[hi] - x[lo])
    if y[lo] <= -30.0 or y[hi] <= -30.0:
        low = 0.0 if y[lo] <= -30.0 else 10.0**y[lo]
        high = 0.0 if y[hi] <= -30.0 else 10.0**y[hi]
        return low + fraction * (high - low)
    return 10.0 ** (y[lo] + fraction * (y[hi] - y[lo]))


def historical_interpolate(curve: tuple[list[float], list[float]], energy_keV: float) -> float:
    energy, value = curve
    if energy_keV <= energy[0]:
        return value[0]
    if energy_keV >= energy[-1]:
        return value[-1]
    lo, hi = 0, len(energy) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if energy[mid] <= energy_keV:
            lo = mid
        else:
            hi = mid
    if value[lo] <= 0.0 or value[hi] <= 0.0:
        fraction = (energy_keV - energy[lo]) / (energy[hi] - energy[lo])
        return value[lo] + fraction * (value[hi] - value[lo])
    fraction = math.log(energy_keV / energy[lo]) / math.log(energy[hi] / energy[lo])
    return math.exp(math.log(value[lo]) + fraction * math.log(value[hi] / value[lo]))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def summarize(cases: list[dict]) -> dict:
    values = [case["difference_percent"] for case in cases]
    if not values:
        return {"count": 0, "median_percent": 0.0, "p95_percent": 0.0,
                "worst_percent": 0.0, "worst_case": None}
    worst = max(cases, key=lambda item: item["difference_percent"])
    return {
        "count": len(values),
        "median_percent": percentile(values, 0.5),
        "p95_percent": percentile(values, 0.95),
        "worst_percent": max(values),
        "worst_case": {key: worst[key] for key in ("Z", "energy_keV", "difference_percent")},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", default=DEFAULT_BASELINE)
    parser.add_argument("--direct-cpp", type=Path, default=DEFAULT_CPP)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args()

    old_csv = git_text(args.baseline_ref, "tools/extract_g4_xs/g4_xs_cache.csv")
    rows = list(csv.DictReader(io.StringIO(old_csv)))
    curves = direct_curves(args.direct_cpp.read_text(encoding="utf-8"))
    by_z: dict[int, list[float]] = defaultdict(list)
    by_z_pe: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for row in rows:
        energy = float(row["E_keV"])
        if 10.0 <= energy <= 10_000.0:
            by_z[int(row["Z"])].append(energy)
            by_z_pe[int(row["Z"])].append((energy, float(row["pe_barn"])))
    edge_nodes = set()
    edge_intervals: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for z, energies in by_z.items():
        energies.sort()
        for low, high in zip(energies, energies[1:]):
            if high / low - 1.0 <= EDGE_PAIR_RELATIVE_SPACING:
                edge_nodes.add((z, low))
                edge_nodes.add((z, high))
                edge_intervals[z].append((low, high))
        # The historical shared grid sometimes represented an absorption jump
        # as a broad interpolation interval instead of a tight flank pair.
        # Classify the whole rising interval as an edge neighborhood; assessing
        # smooth-curve migration error inside a smeared discontinuity is not
        # meaningful and is explicitly outside the away-from-edge gate.
        samples = sorted(by_z_pe[z])
        for (low_e, low_v), (high_e, high_v) in zip(samples, samples[1:]):
            if low_v > 0.0 and high_v / low_v > 1.10:
                edge_nodes.add((z, low_e))
                edge_nodes.add((z, high_e))
                edge_intervals[z].append((low_e, high_e))

    cases = {name: {"away_from_edges": [], "edge_neighborhoods": []} for name in PROCESS}
    historical: dict[tuple[int, str], tuple[list[float], list[float]]] = {}
    for z in range(1, 93):
        z_rows = [row for row in rows if int(row["Z"]) == z and 10.0 <= float(row["E_keV"]) <= 10_000.0]
        for process, column in PROCESS.items():
            historical[z, process] = (
                [float(row["E_keV"]) for row in z_rows],
                [float(row[column]) for row in z_rows],
            )
        direct_x, direct_y = curves[z, "photoelectric"]
        old_energy = historical[z, "photoelectric"][0]
        for low_x, high_x, low_y, high_y in zip(
            direct_x, direct_x[1:], direct_y, direct_y[1:]
        ):
            if high_y - low_y <= math.log10(1.10):
                continue
            edge_energy = 1000.0 * 10.0 ** (0.5 * (low_x + high_x))
            upper = bisect.bisect_right(old_energy, edge_energy)
            if 0 < upper < len(old_energy):
                low, high = old_energy[upper - 1], old_energy[upper]
                edge_nodes.add((z, low))
                edge_nodes.add((z, high))
                edge_intervals[z].append((low, high))
    for row in rows:
        z = int(row["Z"])
        energy = float(row["E_keV"])
        if not 10.0 <= energy <= 10_000.0:
            continue
        edge = (z, energy) in edge_nodes
        bucket = "edge_neighborhoods" if edge else "away_from_edges"
        for process, column in PROCESS.items():
            old = float(row[column])
            new = 0.0 if process == "pair_production" and energy < 1022.0 else interpolate(
                curves[z, process], energy
            )
            if old > 0.0 and not (process == "pair_production" and energy <= 1022.0):
                cases[process][bucket].append({
                    "Z": z,
                    "energy_keV": energy,
                    "difference_percent": abs(new / old - 1.0) * 100.0,
                })

    # Contribution-weighted attenuation is assessed for the actual migration
    # material set. Microscopic element cross sections are weighted by mass
    # fraction / atomic weight before summation, so low-contribution processes
    # and trace constituents do not dominate a percentile merely by case count.
    attenuation = {"away_from_edges": [], "edge_neighborhoods": []}
    # These are the committed downstream material-comparison energies within
    # the restored historical photon-table upper bound of 10 MeV.
    material_energies = [20.0, 50.0, 100.0, 500.0, 1000.0, 3000.0, 10_000.0]
    for material, composition in MATERIALS.items():
        for energy in material_energies:
            edge = any(
                any(low <= energy <= high for low, high in edge_intervals[z])
                for z, _ in composition
            )
            bucket = "edge_neighborhoods" if edge else "away_from_edges"
            old_total = 0.0
            new_total = 0.0
            for z, mass_fraction in composition:
                weight = mass_fraction / ATOMIC_WEIGHT[z]
                for process in PROCESS:
                    old = (0.0 if process == "pair_production" and energy < 1022.0
                           else historical_interpolate(historical[z, process], energy))
                    new = (0.0 if process == "pair_production" and energy < 1022.0
                           else interpolate(curves[z, process], energy))
                    old_total += weight * old
                    new_total += weight * new
            attenuation[bucket].append({
                "material": material,
                "Z": composition[0][0],
                "energy_keV": energy,
                "difference_percent": abs(new_total / old_total - 1.0) * 100.0,
            })

    process_report = {}
    errors = []
    for process, split in cases.items():
        away = summarize(split["away_from_edges"])
        edge = summarize(split["edge_neighborhoods"])
        p95_gate, worst_gate = GATES[process]
        passed = away["p95_percent"] <= p95_gate and away["worst_percent"] <= worst_gate
        if not passed:
            errors.append(
                f"{process}: p95/worst {away['p95_percent']:.6g}%/"
                f"{away['worst_percent']:.6g}% exceeds {p95_gate}%/{worst_gate}%"
            )
        process_report[process] = {
            "away_from_edges": away,
            "edge_neighborhoods_reported_separately": edge,
            "gate_percent": {"p95": p95_gate, "worst": worst_gate},
            "passed": passed,
        }

    attenuation_report = {
        "away_from_edges": summarize(attenuation["away_from_edges"]),
        "edge_neighborhoods_reported_separately": summarize(attenuation["edge_neighborhoods"]),
        "away_from_edges_by_material": {
            material: summarize([
                case for case in attenuation["away_from_edges"]
                if case["material"] == material
            ])
            for material in MATERIALS
        },
        "gate_percent": {"p95": 0.5, "worst": 2.0},
    }
    attenuation_report["passed"] = (
        attenuation_report["away_from_edges"]["p95_percent"] <= 0.5
        and attenuation_report["away_from_edges"]["worst_percent"] <= 2.0
    )
    if not attenuation_report["passed"]:
        away = attenuation_report["away_from_edges"]
        errors.append(
            f"attenuation: p95/worst {away['p95_percent']:.6g}%/"
            f"{away['worst_percent']:.6g}% exceeds 0.5%/2%"
        )

    lock = ROOT / "tools/prepare_cross_sections/sources.lock.json"
    report = {
        "schema_version": 1,
        "source_lock_sha256": hashlib.sha256(lock.read_bytes()).hexdigest(),
        "historical_baseline_ref": args.baseline_ref,
        "historical_input": "transient git object; raw columns are not emitted",
        "historical_comparison_range_keV": [10.0, 10_000.0],
        "runtime_photon_range_keV": [10.0, 10_000.0],
        "attenuation_materials": list(MATERIALS),
        "attenuation_energies_keV": material_energies,
        "edge_definition": (
            "Both nodes of historical adjacent flank pairs separated by at most 0.05% in energy, "
            "plus every historical interpolation interval that contains a greater-than-10% "
            "photoelectric rise in either the historical or direct EPICS curve; the full smeared "
            "discontinuity interval is reported separately."
        ),
        "processes": process_report,
        "contribution_weighted_total_attenuation": attenuation_report,
        "passed": not errors,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [
        "# Direct EPICS2023 photon migration validation\n",
        f"Historical baseline: `{args.baseline_ref}`. Raw historical columns are read transiently and are not reproduced.\n",
        "Absorption-edge neighborhoods include the full historical interpolation bracket around direct or historical discontinuities and are reported separately.\n",
        "| Process | Away p95 | Away worst | Gate p95/worst | Edge p95/worst | Result |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for process, values in process_report.items():
        away = values["away_from_edges"]
        edge = values["edge_neighborhoods_reported_separately"]
        gate = values["gate_percent"]
        lines.append(
            f"| {process} | {away['p95_percent']:.6g}% | {away['worst_percent']:.6g}% "
            f"| {gate['p95']:.6g}% / {gate['worst']:.6g}% "
            f"| {edge['p95_percent']:.6g}% / {edge['worst_percent']:.6g}% "
            f"| {'PASS' if values['passed'] else 'FAIL'} |"
        )
    away = attenuation_report["away_from_edges"]
    edge = attenuation_report["edge_neighborhoods_reported_separately"]
    lines.extend([
        "\n## Contribution-weighted total attenuation\n",
        f"Away from edges: p95 {away['p95_percent']:.6g}%, worst {away['worst_percent']:.6g}% "
        f"(gate 0.5% / 2%): {'PASS' if attenuation_report['passed'] else 'FAIL'}.\n",
        f"Edge neighborhoods (reported separately): p95 {edge['p95_percent']:.6g}%, "
        f"worst {edge['worst_percent']:.6g}%.\n",
    ])
    args.output_markdown.write_text("\n".join(lines), encoding="utf-8")
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"Photon migration gates passed; wrote {args.output_json} and {args.output_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
