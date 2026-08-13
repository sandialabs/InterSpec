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

"""Generate before/after coefficient, accuracy, and compiled resident-data accounting."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPORTS = Path(__file__).resolve().parent / "reports"
BASELINES = {
    "main": "ab03e4810c548126daced9d7c91d62b0467a3d6b",
    "feature/FullSpecApprox": "7e8713b517c13999a148bea2578cc3f26d99e67f",
    "feature/env-scatter": "4977fccff48fd26bdcc9bb873577c72ff0973183",
}
OBJECTS = {
    "photon_epics_data.cpp.o": "src/CMakeFiles/CeeLo.dir/cross_sections/photon_epics_data.cpp.o",
    "element_data.cpp.o": "src/CMakeFiles/CeeLo.dir/cross_sections/element_data.cpp.o",
    "relaxation_epics_data.cpp.o": "src/CMakeFiles/CeeLo.dir/cross_sections/relaxation_epics_data.cpp.o",
    "estar_stopping_data.cpp.o": "src/CMakeFiles/CeeLo.dir/physics/estar_stopping_data.cpp.o",
}
RUNTIME_OBJECT = "src/CMakeFiles/CeeLo.dir/cross_sections/CrossSectionData.cpp.o"
PHOTON_NAMES = (
    "rayleigh", "compton", "pair_production", "photoelectric",
    "k_photoelectric", "form_factor", "scattering_function",
)
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def git_text(ref: str, path: str) -> str:
    return subprocess.run(
        ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True,
    ).stdout


def array_counts(text: str, stem: str) -> list[int]:
    counts = []
    for z in range(1, 93):
        match = re.search(
            rf"static const (?:float|uint16_t) Z{z}_{stem}(?:\[[^]]*\])?\s*=\s*\{{(.*?)\}};",
            text, re.DOTALL,
        )
        if not match:
            raise ValueError(f"missing historical Z{z}_{stem}")
        body = re.sub(r"//.*", "", match.group(1))
        counts.append(len(NUMBER.findall(body)))
    return counts


def percentile(values: list[int], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return float(ordered[low])
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def node_stats(values: list[int]) -> dict[str, float]:
    return {"min": min(values), "median": percentile(values, 0.5),
            "p95": percentile(values, 0.95), "max": max(values)}


def compiled_constant_bytes(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(f"build the CeeLo target before accounting: {path}")
    # Darwin's size(1) identifies read-only constants and any initialized data
    # explicitly. Both remain resident.
    darwin = subprocess.run(
        ["size", "-m", str(path)], check=False, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True,
    )
    matches = re.findall(
        r"Section \([^,]+, __(?:const|data)\):\s+(\d+)", darwin.stdout
    )
    if matches:
        return sum(int(value) for value in matches)
    # GNU size: include .rodata, initialized data, and relocated read-only
    # descriptor sections.
    portable = subprocess.run(
        ["size", "-A", str(path)], check=True, stdout=subprocess.PIPE, text=True,
    ).stdout
    total = 0
    for line in portable.splitlines():
        fields = line.split()
        if len(fields) >= 2 and (
            fields[0].startswith(".rodata")
            or fields[0].startswith(".data.rel.ro")
            or fields[0] == ".data"
        ):
            total += int(fields[1], 0)
    if total <= 0:
        raise ValueError(f"could not identify compiled constant sections in {path}")
    return total


def measured_section_bytes(path: Path, darwin_name: str, gnu_name: str) -> int:
    darwin = subprocess.run(
        ["size", "-m", str(path)], check=False, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True,
    )
    match = re.search(
        rf"Section \([^,]+, {re.escape(darwin_name)}\):\s+(\d+)", darwin.stdout
    )
    if match:
        return int(match.group(1))
    # A successful Darwin report with no named section means the object has no
    # such allocation (not that GNU-style parsing is required).
    if darwin.returncode == 0 and "Segment" in darwin.stdout:
        return 0
    portable = subprocess.run(
        ["size", "-A", str(path)], check=True, stdout=subprocess.PIPE, text=True,
    ).stdout
    return sum(
        int(fields[1], 0) for line in portable.splitlines()
        if len((fields := line.split())) >= 2 and fields[0] == gnu_name
    )


def row(
    baseline: str,
    process: str,
    before_count: int,
    before_encoding: str,
    before_coefficient_bytes: int,
    before_grid_bytes: int,
    before_metadata_bytes: int,
    before_nodes: dict[str, float],
    after_count: int,
    native_after: int,
    after_encoding: str,
    after_coefficient_bytes: int,
    after_grid_bytes: int,
    after_metadata_bytes: int,
    after_nodes: dict[str, float],
    p95_error: float,
    max_error: float,
) -> dict:
    before_total = before_coefficient_bytes + before_grid_bytes + before_metadata_bytes
    after_total = after_coefficient_bytes + after_grid_bytes + after_metadata_bytes
    change = after_total - before_total
    return {
        "baseline": baseline,
        "process": process,
        "native_coefficients_before": before_count,
        "retained_coefficients_before": before_count,
        "scalar_encoding_before": before_encoding,
        "coefficient_bytes_before": before_coefficient_bytes,
        "grid_bytes_before": before_grid_bytes,
        "metadata_bytes_before": before_metadata_bytes,
        "compiled_rodata_bytes_before": before_total,
        "nodes_min_before": before_nodes["min"],
        "nodes_median_before": before_nodes["median"],
        "nodes_p95_before": before_nodes["p95"],
        "nodes_max_before": before_nodes["max"],
        "native_coefficients_after": native_after,
        "retained_coefficients_after": after_count,
        "scalar_encoding_after": after_encoding,
        "coefficient_bytes_after": after_coefficient_bytes,
        "grid_bytes_after": after_grid_bytes,
        "metadata_bytes_after": after_metadata_bytes,
        "compiled_rodata_bytes_after": after_total,
        "nodes_min_after": after_nodes["min"],
        "nodes_median_after": after_nodes["median"],
        "nodes_p95_after": after_nodes["p95"],
        "nodes_max_after": after_nodes["max"],
        "absolute_byte_change": change,
        "percentage_byte_change": 100.0 * change / before_total if before_total else None,
        "p95_interpolation_error_percent": p95_error,
        "maximum_interpolation_error_percent": max_error,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, default=REPORTS)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build-rel")
    args = parser.parse_args()
    photon = json.loads((args.report_dir / "photon_epics.json").read_text())
    brems = json.loads((args.report_dir / "bremsstrahlung_epq.json").read_text())
    stopping = json.loads((args.report_dir / "stopping_estar.json").read_text())
    relaxation = json.loads((args.report_dir / "relaxation_eadl.json").read_text())
    direct = {entry["process"]: entry for entry in photon["processes"]}
    rows: list[dict] = []

    for baseline, ref in BASELINES.items():
        historical = git_text(ref, "src/cross_sections/element_data.cpp")
        energy_nodes = array_counts(historical, "log_energy")
        angular_nodes = [18] * 92
        for process in PHOTON_NAMES:
            after = direct[process]
            before_nodes = angular_nodes if process in ("form_factor", "scattering_function") else energy_nodes
            before_count = sum(before_nodes)
            # Historical photon energy grids were shared across five processes;
            # angular grids were each one shared 18-node array.
            before_grid = 0
            after_nodes = {key.removeprefix("nodes_"): after[key]
                           for key in ("nodes_min", "nodes_median", "nodes_p95", "nodes_max")}
            rows.append(row(
                baseline, process, before_count, "float32",
                before_count * 4, before_grid, 0, node_stats(before_nodes),
                after["retained_coefficients"], after["native_coefficients"],
                after["encoding"], after["coefficient_bytes"], after["grid_bytes"],
                after["metadata_bytes"] + after.get("legacy_compatibility_bytes", 0),
                after_nodes, after["p95_interpolation_error_percent"],
                after["maximum_interpolation_error_percent"],
            ))

        quantized = "static const uint16_t Z1_sb_chi" in historical
        sb_nodes = array_counts(historical, "sb_chi")
        before_coeff = sum(sb_nodes)
        before_scalar = "uint16" if quantized else "float32"
        before_bytes = before_coeff * (2 if quantized else 4)
        before_grid = (32 + len(sb_nodes[0:1]) * 0 + (26 if quantized else 27)) * 4
        before_metadata = 92 * 4 + 4 if quantized else 4
        storage = brems["storage"]
        packing = brems["validation"]["uint16_packing"]
        packing_p95 = max(
            packing[name]["p95_percent"] for name in (
                "row_rms_percent", "integrated_number_rate_difference_percent",
                "mean_kappa_difference_percent",
            )
        )
        packing_max = max(
            packing[name]["maximum_percent"] for name in (
                "row_rms_percent", "integrated_number_rate_difference_percent",
                "mean_kappa_difference_percent",
            )
        )
        rows.append(row(
            baseline, "bremsstrahlung", before_coeff, before_scalar, before_bytes,
            before_grid, before_metadata, node_stats([count // 32 for count in sb_nodes]),
            storage["retained_coefficients"], storage["native_coefficients_in_runtime_range"],
            storage["scalar_encoding"], storage["coefficient_bytes"],
            storage["energy_grid_bytes"] + storage["kappa_grid_bytes"], storage["metadata_bytes"],
            node_stats([27] * 92), packing_p95, packing_max,
        ))

        # The former stopping representation was one shared 12x15 correction
        # surface per process; the new representation retains all 92x53 values.
        for process in ("collision_stopping", "radiative_stopping"):
            after = stopping["processes"][process]
            rows.append(row(
                baseline, process, 180, "float32 correction ratio", 720, 0, 0,
                node_stats([15] * 12), after["retained_coefficients"],
                after["native_coefficients"], stopping["encoding"],
                after["coefficient_bytes"], 0, after["metadata_bytes"], node_stats([53] * 92),
                after["p95_native_representation_error_percent"],
                after["maximum_native_representation_error_percent"],
            ))

        # Relaxation stores only populated line arrays. Pointer/count/alignment
        # overhead is reported separately from the retained float coefficients.
        relaxation_elements = relaxation["maximum_atomic_number"]
        relaxation_slots = relaxation_elements + 1  # direct Z indexing includes Z=0
        k_relaxation = relaxation["processes"]["K_relaxation"]
        l_relaxation = relaxation["processes"]["L_relaxation"]
        relaxation_p95_error = max(
            relaxation["p95_float32_yield_error_percent"],
            relaxation["p95_float32_mean_emitted_energy_error_percent"],
        )
        relaxation_max_error = max(
            relaxation["maximum_float32_yield_error_percent"],
            relaxation["maximum_float32_mean_emitted_energy_error_percent"],
        )
        rows.append(row(
            baseline, "K_relaxation", 92 * 14, "float32", 92 * 14 * 4,
            0, 92 * 4, node_stats([6] * 92),
            k_relaxation["retained_coefficients"], k_relaxation["native_coefficients"],
            "float32 variable-length", k_relaxation["coefficient_bytes"], 0,
            k_relaxation["metadata_bytes"],
            node_stats([item["line_count"] for item in relaxation["subshells"] if item["shell"] == "K"]),
            relaxation_p95_error, relaxation_max_error,
        ))
        rows.append(row(
            baseline, "L_relaxation", 100 * 43, "float32", 100 * 43 * 4,
            0, 100 * 12, node_stats([6] * 100),
            l_relaxation["retained_coefficients"], l_relaxation["native_coefficients"],
            "float32 variable-length", l_relaxation["coefficient_bytes"], 0,
            l_relaxation["metadata_bytes"],
            node_stats([item["line_count"] for item in relaxation["subshells"] if item["shell"] != "K"]),
            relaxation_p95_error, relaxation_max_error,
        ))

    shared = []
    for baseline, ref in BASELINES.items():
        historical = git_text(ref, "src/cross_sections/element_data.cpp")
        energy_nodes = array_counts(historical, "log_energy")
        shared.extend([
            {"baseline": baseline, "grid": "historical_photon_energy_grids",
             "bytes_before": sum(energy_nodes) * 4, "bytes_after": 0,
             "notes": "Previously shared by five photon processes; never double-counted above."},
            {"baseline": baseline, "grid": "historical_scattering_function_grid",
             "bytes_before": 18 * 4, "bytes_after": 0,
             "notes": "Direct EPDL now uses per-element process-specific grids."},
            {"baseline": baseline, "grid": "historical_form_factor_grid",
             "bytes_before": 18 * 4, "bytes_after": 0,
             "notes": "Direct EPDL now uses per-element process-specific grids."},
            {"baseline": baseline, "grid": "ESTAR_energy_grid",
             "bytes_before": 0, "bytes_after": stopping["shared_grid"]["bytes"],
             "notes": "Shared by collision and radiative stopping; counted once."},
            {"baseline": baseline, "grid": "EPDL_photon_energy_coordinate_pool",
             "bytes_before": 0,
             "bytes_after": photon["shared_energy_coordinate_pool"]["bytes"],
             "notes": "Bit-exact float32 pool shared by five process-specific uint16 index streams."},
            {"baseline": baseline, "grid": "Rayleigh_offline_sampling_CDF",
             "bytes_before": 0,
             "bytes_after": (
                 photon["rayleigh_sampling_cdf"]["coefficient_bytes"]
                 + photon["rayleigh_sampling_cdf"]["shared_grid_bytes"]
             ),
             "notes": "Normalized uint16 CDF for all elements plus one shared float32 x grid; generated offline."},
        ])

    compiled = {
        name: compiled_constant_bytes(args.build_dir / relative)
        for name, relative in OBJECTS.items()
    }
    runtime_object = args.build_dir / RUNTIME_OBJECT
    singleton_static = compiled_constant_bytes(runtime_object)
    singleton_bss = measured_section_bytes(runtime_object, "__bss", ".bss")
    per_thread_bss = measured_section_bytes(runtime_object, "__thread_bss", ".tbss")
    runtime = {
        "singleton_static_data_measured": singleton_static,
        "singleton_and_guard_bss_measured": singleton_bss,
        "thread_local_bss_per_thread_measured": per_thread_bss,
    }
    runtime_bytes = singleton_static + singleton_bss + per_thread_bss
    combined_bytes = sum(compiled.values()) + runtime_bytes
    report = {
        "schema_version": 1,
        "source_lock_sha256": hashlib.sha256(
            (Path(__file__).resolve().parent / "sources.lock.json").read_bytes()
        ).hexdigest(),
        "baselines": BASELINES,
        "rows": rows,
        "shared_grids": shared,
        "compiled_object_constant_bytes_after": compiled,
        "compiled_constant_bytes_total_after": sum(compiled.values()),
        "initialized_runtime_bytes_after": runtime,
        "initialized_runtime_bytes_total_after_one_thread": runtime_bytes,
        "compiled_plus_initialized_bytes_total_after": combined_bytes,
        "memory_gate_bytes": 1024 * 1024,
        "passed": combined_bytes < 1024 * 1024,
        "notes": [
            "Compiled resident data includes the immutable compact tables; measured Mach-O section totals are listed separately.",
            "Shared historical grids and the shared ESTAR grid are reported separately and not charged to individual process rows.",
            "Every omitted EPDL node remains a validation probe; retained nodes are selected by accuracy rather than a fixed memory target.",
            "The under-1-MiB gate includes measured singleton BSS and one thread's measured TLS cache; each additional active thread adds the separately reported TLS bytes.",
        ],
    }
    json_path = args.report_dir / "coefficient_sizes.json"
    csv_path = args.report_dir / "coefficient_sizes.csv"
    md_path = args.report_dir / "coefficient_sizes.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    lines = [
        "# Cross-section coefficient and resident-data accounting\n",
        "Shared grids are listed separately. Size increases are accepted only where the committed accuracy report passes its gate.\n",
        "| Baseline | Process | Before coeffs / encoding | After retained / native | Before bytes | After bytes | Change | p95 / max error |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for item in rows:
        before_total = item["compiled_rodata_bytes_before"]
        after_total = item["compiled_rodata_bytes_after"]
        lines.append(
            f"| {item['baseline']} | {item['process']} | {item['retained_coefficients_before']:,} / {item['scalar_encoding_before']} "
            f"| {item['retained_coefficients_after']:,} / {item['native_coefficients_after']:,} "
            f"| {before_total:,} | {after_total:,} | {item['absolute_byte_change']:+,} "
            f"({item['percentage_byte_change']:+.1f}%) | {item['p95_interpolation_error_percent']:.6g}% / "
            f"{item['maximum_interpolation_error_percent']:.6g}% |"
        )
    lines.extend(["\n## Shared grids\n", "| Baseline | Grid | Before bytes | After bytes | Note |",
                  "|---|---|---:|---:|---|"])
    for item in shared:
        lines.append(f"| {item['baseline']} | {item['grid']} | {item['bytes_before']:,} | "
                     f"{item['bytes_after']:,} | {item['notes']} |")
    lines.extend(["\n## Measured compiled resident-data sections\n",
                  "| Object | Bytes |", "|---|---:|"])
    for name, size in report["compiled_object_constant_bytes_after"].items():
        lines.append(f"| `{name}` | {size:,} |")
    lines.extend([
        f"\nCompiled resident data total: {report['compiled_constant_bytes_total_after']:,} bytes.",
        f"Initialized singleton + one thread TLS: {report['initialized_runtime_bytes_total_after_one_thread']:,} bytes.",
        f"Combined for one active thread: {report['compiled_plus_initialized_bytes_total_after']:,} bytes "
        f"({report['compiled_plus_initialized_bytes_total_after'] / 1048576.0:.3f} MiB).",
    ])
    md_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {csv_path}, {json_path}, and {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
