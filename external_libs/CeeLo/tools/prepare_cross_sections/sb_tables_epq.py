#!/usr/bin/env python3
# CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
# and X-ray detector efficiency - developed as part of InterSpec.
#
# Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License, version 2.1 or later.

"""Generate the G4-grid-aligned bremsstrahlung table from public-domain NIST EPQ.

This production generator reads no Geant4 or G4EMLOW files.  Two missing/revised
low-kappa samples are obtained by one element- and energy-independent polynomial
interpolation rule.  All other samples are copied from equal-kappa EPQ nodes.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import struct
import sys
from pathlib import Path
from typing import Any, Iterable

from fetch_sources import epq_table_set_sha256, load_lock
from generation_utils import (
    ELEMENT_SYMBOLS,
    OUT_KAPPA,
    OUT_T_keV,
    Z_MAX,
    as_float32,
    replace_section,
    write_section,
)


HERE = Path(__file__).resolve().parent
DEFAULT_LOCK = HERE / "sources.lock.json"
DEFAULT_REPORT_JSON = HERE / "reports" / "bremsstrahlung_epq.json"
DEFAULT_REPORT_MARKDOWN = HERE / "reports" / "bremsstrahlung_epq.md"

EPQ_KAPPA = [
    0.0, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925,
    0.95, 0.97, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1.0,
]
EPQ_NATIVE_ENERGIES = 57
NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?")

# Constants needed only to turn the tabulated first moment back into the
# dimensionless phi_rad definition (Seltzer-Berger 1986, Eq. 5).
ALPHA = 1.0 / 137.035999084
CLASSICAL_ELECTRON_RADIUS_M = 2.8179403262e-15
ELECTRON_REST_ENERGY_EV = 510_998.950
MILLIBARN_M2 = 1.0e-31


class EpqTable:
    def __init__(self, z: int, energies_eV: list[float], rows: list[list[float]], phi_rad: list[float]):
        self.z = z
        self.energies_eV = energies_eV
        self.rows = rows
        self.phi_rad = phi_rad


def parse_epq_table(path: Path, expected_z: int) -> EpqTable:
    values = [float(token) for token in NUMBER.findall(path.read_text(encoding="ascii"))]
    expected = 1 + EPQ_NATIVE_ENERGIES * (1 + len(EPQ_KAPPA) + 1)
    if len(values) != expected:
        raise ValueError(f"{path}: expected {expected} numeric fields, found {len(values)}")
    z = int(values[0])
    if z != expected_z or values[0] != z:
        raise ValueError(f"{path}: element number is {values[0]}, expected {expected_z}")

    energies: list[float] = []
    rows: list[list[float]] = []
    phi_rad: list[float] = []
    width = 1 + len(EPQ_KAPPA) + 1
    for index in range(EPQ_NATIVE_ENERGIES):
        offset = 1 + index * width
        energy = values[offset]
        row = values[offset + 1:offset + 1 + len(EPQ_KAPPA)]
        phi = values[offset + 1 + len(EPQ_KAPPA)]
        if energy <= 0.0 or (energies and energy <= energies[-1]):
            raise ValueError(f"{path}: non-increasing energy at native row {index}")
        if any(not math.isfinite(value) or value <= 0.0 for value in row):
            raise ValueError(f"{path}: non-positive/non-finite cross section at row {index}")
        if not math.isfinite(phi) or phi <= 0.0:
            raise ValueError(f"{path}: invalid phi_rad at row {index}")
        energies.append(energy)
        rows.append(row)
        phi_rad.append(phi)
    return EpqTable(z, energies, rows, phi_rad)


def lagrange(xs: Iterable[float], ys: Iterable[float], x: float) -> float:
    x_values, y_values = list(xs), list(ys)
    if len(x_values) != len(y_values) or not x_values:
        raise ValueError("Lagrange interpolation requires equal non-empty inputs")
    result = 0.0
    for i, yi in enumerate(y_values):
        basis = 1.0
        for j, xj in enumerate(x_values):
            if i != j:
                basis *= (x - xj) / (x_values[i] - xj)
        result += yi * basis
    return result


def transform_row(row: list[float]) -> list[float]:
    """Apply the fixed EPQ-to-CeeLo kappa-grid transformation."""
    if len(row) != len(EPQ_KAPPA):
        raise ValueError(f"expected {len(EPQ_KAPPA)} kappa values, got {len(row)}")
    by_kappa = {kappa: value for kappa, value in zip(EPQ_KAPPA, row)}
    output: list[float] = []
    for kappa in OUT_KAPPA:
        if kappa == 0.025:
            value = lagrange(
                (0.0, 0.05, 0.1),
                (by_kappa[0.0], by_kappa[0.05], by_kappa[0.1]),
                kappa,
            )
        elif kappa == 0.075:
            value = lagrange(
                (0.0, 0.05, 0.1, 0.15),
                (by_kappa[0.0], by_kappa[0.05], by_kappa[0.1], by_kappa[0.15]),
                kappa,
            )
        else:
            value = by_kappa[kappa]
        value = as_float32(value)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"low-kappa transformation produced invalid value {value}")
        output.append(value)
    return output


def trapezoid(xs: list[float], ys: list[float]) -> float:
    return sum(
        0.5 * (ys[index] + ys[index + 1]) * (xs[index + 1] - xs[index])
        for index in range(len(xs) - 1)
    )


def implied_phi_rad(energy_eV: float, kappa: list[float], row: list[float]) -> float:
    moment = trapezoid(kappa, row)
    ratio = (energy_eV + ELECTRON_REST_ENERGY_EV) / (
        energy_eV + 2.0 * ELECTRON_REST_ENERGY_EV
    )
    return (
        moment
        * MILLIBARN_M2
        / (ALPHA * CLASSICAL_ELECTRON_RADIUS_M**2)
        * ratio
    )


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("cannot take percentile of an empty sequence")
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def summarize_percent(values: list[float]) -> dict[str, float]:
    return {
        "median_percent": percentile(values, 0.50),
        "p95_percent": percentile(values, 0.95),
        "p99_percent": percentile(values, 0.99),
        "maximum_percent": max(values),
    }


def table_payload_sha256(tables: dict[int, list[list[float]]]) -> str:
    digest = hashlib.sha256()
    for z in range(1, Z_MAX + 1):
        for row in tables[z]:
            for value in row:
                digest.update(struct.pack("<f", as_float32(value)))
    return digest.hexdigest()


def quantize_tables(
    tables: dict[int, list[list[float]]]
) -> tuple[dict[int, list[list[float]]], str]:
    decoded: dict[int, list[list[float]]] = {}
    digest = hashlib.sha256()
    for z in range(1, Z_MAX + 1):
        maximum = max(value for row in tables[z] for value in row)
        scale = as_float32(maximum / 65535.0) if maximum > 0.0 else 0.0
        digest.update(struct.pack("<f", scale))
        rows = []
        for row in tables[z]:
            packed = [
                min(65535, max(0, 0 if scale == 0.0 else round(value / scale)))
                for value in row
            ]
            for code in packed:
                digest.update(struct.pack("<H", code))
            rows.append([code * scale for code in packed])
        decoded[z] = rows
    return decoded, digest.hexdigest()


def row_observables(row: list[float], energy_keV: float) -> tuple[float, float] | None:
    kmin = 10.0 / energy_keV
    if kmin >= 1.0:
        return None

    def value_at(kappa: float) -> float:
        import bisect
        high = min(max(bisect.bisect_right(OUT_KAPPA, kappa), 1), len(OUT_KAPPA) - 1)
        low = high - 1
        fraction = (kappa - OUT_KAPPA[low]) / (OUT_KAPPA[high] - OUT_KAPPA[low])
        return row[low] + fraction * (row[high] - row[low])

    x = [kmin] + [kappa for kappa in OUT_KAPPA if kappa > kmin]
    y = [value_at(kappa) for kappa in x]
    number = sum(
        (x[index + 1] - x[index]) * 0.5
        * (y[index] / x[index] + y[index + 1] / x[index + 1])
        for index in range(len(x) - 1)
    )
    energy = sum(
        (x[index + 1] - x[index]) * 0.5 * (y[index] + y[index + 1])
        for index in range(len(x) - 1)
    )
    return number, energy / number


def quantization_validation(tables: dict[int, list[list[float]]]) -> tuple[dict[str, Any], str]:
    decoded, payload_hash = quantize_tables(tables)
    row_rms, rate_error, mean_error = [], [], []
    for z in range(1, Z_MAX + 1):
        for energy, reference, represented in zip(OUT_T_keV, tables[z], decoded[z]):
            row_rms.append(math.sqrt(sum(
                (candidate / source - 1.0) ** 2
                for source, candidate in zip(reference, represented) if source > 0.0
            ) / sum(source > 0.0 for source in reference)) * 100.0)
            source_observables = row_observables(reference, energy)
            represented_observables = row_observables(represented, energy)
            if source_observables and represented_observables:
                rate_error.append(abs(represented_observables[0] / source_observables[0] - 1.0) * 100.0)
                mean_error.append(abs(represented_observables[1] / source_observables[1] - 1.0) * 100.0)
    validation = {
        "row_rms_percent": summarize_percent(row_rms),
        "integrated_number_rate_difference_percent": summarize_percent(rate_error),
        "mean_kappa_difference_percent": summarize_percent(mean_error),
        "gates_percent": {
            "row_rms_p95": 0.02, "row_rms_maximum": 0.25,
            "integrated_rate_p95": 0.005, "integrated_rate_maximum": 0.01,
            "mean_kappa_p95": 0.005, "mean_kappa_maximum": 0.01,
        },
    }
    gates = validation["gates_percent"]
    validation["passed"] = (
        validation["row_rms_percent"]["p95_percent"] <= gates["row_rms_p95"]
        and validation["row_rms_percent"]["maximum_percent"] <= gates["row_rms_maximum"]
        and validation["integrated_number_rate_difference_percent"]["p95_percent"] <= gates["integrated_rate_p95"]
        and validation["integrated_number_rate_difference_percent"]["maximum_percent"] <= gates["integrated_rate_maximum"]
        and validation["mean_kappa_difference_percent"]["p95_percent"] <= gates["mean_kappa_p95"]
        and validation["mean_kappa_difference_percent"]["maximum_percent"] <= gates["mean_kappa_maximum"]
    )
    return validation, payload_hash


def epq_provenance(source: dict[str, Any]) -> str:
    return "".join(
        [
            "// chi(Z, T, kappa) = (beta^2/Z^2) * kappa * dSigma/dkappa\n",
            f"// {len(OUT_KAPPA)} kappa values x {len(OUT_T_keV)} native electron energies ",
            f"({OUT_T_keV[0]:g} keV - {OUT_T_keV[-1] / 1000.0:g} MeV)\n",
            "//\n",
            "// SOURCE AND TERMS\n",
            "//   NIST EPQ BergerSeltzerBrem tables, public domain in the United States.\n",
            f"//   Pinned EPQ commit: {source['version'].removeprefix('git-')}\n",
            f"//   Verified 92-table SHA-256: {source['table_set_sha256']}\n",
            "//   Retain the upstream NIST notice; neither NIST nor the authors endorse\n",
            "//   derived products. See tools/prepare_cross_sections/sources.lock.json.\n",
            "//   Underlying reference: S.M. Seltzer and M.J. Berger, Atomic Data and\n",
            "//   Nuclear Data Tables 35 (1986) 345-418. EPQ notes PENELOPE provenance.\n",
            "//\n",
            "// FIXED, UNIFORM TRANSFORMATION\n",
            "//   * kappa=0.025: quadratic interpolation through EPQ {0, 0.05, 0.1}.\n",
            "//   * kappa=0.075: cubic interpolation through EPQ {0, 0.05, 0.1, 0.15}.\n",
            "//   * Copy every other equal-kappa node; omit EPQ kappa=0.125.\n",
            "//   * Select EPQ's 27 native incident energies from 10 keV through 20 MeV.\n",
            "//   * Preserve the differential magnitude. phi_rad is checked but is not\n",
            "//     used to normalize rows. Values use per-element uint16 scales.\n",
            "//\n",
            "// This is a general interpolation of public-domain EPQ data, validated\n",
            "// against Geant4 behavior. It is not a transcription, fit, or reproduction\n",
            "// of G4EMLOW coefficients. This generator reads no Geant4/G4EMLOW data.\n",
            "//\n",
            "// GENERATED by tools/prepare_cross_sections/sb_tables_epq.py -- do not\n",
            "// hand-edit. Regenerate through the verified preparation cache; normal\n",
            "// builds remain offline.\n",
            "// =====================================================================\n\n",
        ]
    )


def build_tables(
    directory: Path,
) -> tuple[dict[int, list[list[float]]], list[float], list[float]]:
    tables: dict[int, list[list[float]]] = {}
    source_phi_errors: list[float] = []
    output_phi_errors: list[float] = []
    target_eV = [energy * 1000.0 for energy in OUT_T_keV]

    for z in range(1, Z_MAX + 1):
        parsed = parse_epq_table(directory / f"pdebr{z:02d}.tab", z)
        energy_index = {energy: index for index, energy in enumerate(parsed.energies_eV)}
        missing = [energy for energy in target_eV if energy not in energy_index]
        if missing:
            raise ValueError(f"Z={z}: requested energies are not native EPQ nodes: {missing}")

        output_rows: list[list[float]] = []
        for energy in target_eV:
            index = energy_index[energy]
            source_row = parsed.rows[index]
            output_row = transform_row(source_row)
            phi = parsed.phi_rad[index]
            source_phi_errors.append(
                abs(implied_phi_rad(energy, EPQ_KAPPA, source_row) / phi - 1.0) * 100.0
            )
            output_phi_errors.append(
                abs(implied_phi_rad(energy, OUT_KAPPA, output_row) / phi - 1.0) * 100.0
            )
            output_rows.append(output_row)
        tables[z] = output_rows
        print(f"\r  transformed Z={z:3d}/{Z_MAX}", end="", file=sys.stderr)
    print(file=sys.stderr)
    return tables, source_phi_errors, output_phi_errors


def make_report(
    lock_path: Path,
    source: dict[str, Any],
    tables: dict[int, list[list[float]]],
    source_phi_errors: list[float],
    output_phi_errors: list[float],
) -> dict[str, Any]:
    output_count = Z_MAX * len(OUT_T_keV) * len(OUT_KAPPA)
    packing_validation, packed_hash = quantization_validation(tables)
    return {
        "schema_version": 1,
        "process": "bremsstrahlung_differential_scaled_chi",
        "generator": "tools/prepare_cross_sections/sb_tables_epq.py",
        "source_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "source": {
            "name": "nist_epq",
            "version": source["version"],
            "archive_sha256": source["sha256"],
            "table_set_sha256": source["table_set_sha256"],
        },
        "grid": {
            "elements": Z_MAX,
            "source_incident_energy_nodes_all": EPQ_NATIVE_ENERGIES,
            "native_incident_energy_nodes_in_runtime_range": len(OUT_T_keV),
            "retained_incident_energy_nodes": len(OUT_T_keV),
            "source_kappa_nodes": len(EPQ_KAPPA),
            "retained_kappa_nodes": len(OUT_KAPPA),
            "energy_min_keV": OUT_T_keV[0],
            "energy_max_keV": OUT_T_keV[-1],
        },
        "storage": {
            "scalar_encoding": "uint16 with one float32 scale per element",
            "source_coefficients_all_energies": Z_MAX * EPQ_NATIVE_ENERGIES * len(EPQ_KAPPA),
            "native_coefficients_in_runtime_range": output_count,
            "retained_coefficients": output_count,
            "coefficient_bytes": output_count * 2,
            "kappa_grid_bytes": len(OUT_KAPPA) * 4,
            "energy_grid_bytes": len(OUT_T_keV) * 4,
            "metadata_bytes": 4 + Z_MAX * 4,
            "logical_total_bytes": output_count * 2 + len(OUT_KAPPA) * 4
                + len(OUT_T_keV) * 4 + 4 + Z_MAX * 4,
        },
        "source_derived_payload_sha256_float32_little_endian": table_payload_sha256(tables),
        "packed_payload_sha256_scale_float32_then_uint16_little_endian": packed_hash,
        "validation": {
            "epq_phi_rad_source_grid_absolute_relative_error": summarize_percent(source_phi_errors),
            "epq_phi_rad_output_grid_absolute_relative_error": summarize_percent(output_phi_errors),
            "uint16_packing": packing_validation,
            "note": "phi_rad is an independent first-moment checksum and is not a row-normalization target.",
        },
        "passed": packing_validation["passed"],
    }


def markdown_report(report: dict[str, Any]) -> str:
    storage = report["storage"]
    grid = report["grid"]
    validation = report["validation"]
    source_stats = validation["epq_phi_rad_source_grid_absolute_relative_error"]
    output_stats = validation["epq_phi_rad_output_grid_absolute_relative_error"]
    packing = validation["uint16_packing"]
    return f"""# G4-aligned NIST EPQ bremsstrahlung generation report

This report is generated by `tools/prepare_cross_sections/sb_tables_epq.py` from
the pinned source lock. The production generator reads no Geant4/G4EMLOW data.

| Field | Value |
|---|---:|
| Source | `{report['source']['version']}` |
| Verified EPQ table-set SHA-256 | `{report['source']['table_set_sha256']}` |
| Elements | {grid['elements']} |
| Native energies retained (10 keV–20 MeV) | {grid['retained_incident_energy_nodes']} / {grid['source_incident_energy_nodes_all']} |
| Source/output κ nodes | {grid['source_kappa_nodes']} / {grid['retained_kappa_nodes']} |
| Retained coefficients | {storage['retained_coefficients']:,} |
| Encoding | {storage['scalar_encoding']} |
| Coefficient bytes | {storage['coefficient_bytes']:,} ({storage['coefficient_bytes'] / 1024.0:.1f} KiB) |
| Separate grid + metadata bytes | {storage['kappa_grid_bytes'] + storage['energy_grid_bytes'] + storage['metadata_bytes']:,} |
| Pre-packing float32 payload SHA-256 | `{report['source_derived_payload_sha256_float32_little_endian']}` |
| Packed payload SHA-256 | `{report['packed_payload_sha256_scale_float32_then_uint16_little_endian']}` |

## EPQ `phi_rad` first-moment checksum

`phi_rad` is reported as an independent checksum. It is deliberately not used to
renormalize differential rows.

| Grid | Median | p95 | p99 | Maximum |
|---|---:|---:|---:|---:|
| Original EPQ | {source_stats['median_percent']:.6f}% | {source_stats['p95_percent']:.6f}% | {source_stats['p99_percent']:.6f}% | {source_stats['maximum_percent']:.6f}% |
| G4-aligned interpolation | {output_stats['median_percent']:.6f}% | {output_stats['p95_percent']:.6f}% | {output_stats['p99_percent']:.6f}% | {output_stats['maximum_percent']:.6f}% |

## uint16 storage error

| Metric | Median | p95 | Maximum |
|---|---:|---:|---:|
| Row RMS | {packing['row_rms_percent']['median_percent']:.6f}% | {packing['row_rms_percent']['p95_percent']:.6f}% | {packing['row_rms_percent']['maximum_percent']:.6f}% |
| Integrated number rate | {packing['integrated_number_rate_difference_percent']['median_percent']:.6f}% | {packing['integrated_number_rate_difference_percent']['p95_percent']:.6f}% | {packing['integrated_number_rate_difference_percent']['maximum_percent']:.6f}% |
| Mean kappa | {packing['mean_kappa_difference_percent']['median_percent']:.6f}% | {packing['mean_kappa_difference_percent']['p95_percent']:.6f}% | {packing['mean_kappa_difference_percent']['maximum_percent']:.6f}% |

The only synthesized samples are κ=0.025 (quadratic through 0, 0.05, 0.1) and
κ=0.075 (cubic through 0, 0.05, 0.1, 0.15). Every element and energy uses the
same rule; all other output samples are direct equal-node EPQ values.
"""


def write_report(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--epq", type=Path, required=True, help="directory containing pdebr01.tab ... pdebr92.tab")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--output", type=Path, help="element_data.cpp to patch (default: section on stdout)")
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT_JSON)
    parser.add_argument("--report-markdown", type=Path, default=DEFAULT_REPORT_MARKDOWN)
    parser.add_argument("--no-report", action="store_true")
    args = parser.parse_args()

    lock = load_lock(args.lock)
    source = lock["sources"]["nist_epq"]
    actual_digest = epq_table_set_sha256(args.epq, source["table_glob"])
    if actual_digest != source["table_set_sha256"]:
        raise ValueError(
            f"EPQ table-set SHA-256 mismatch: expected {source['table_set_sha256']}, "
            f"got {actual_digest}"
        )

    tables, source_phi_errors, output_phi_errors = build_tables(args.epq)
    provenance = epq_provenance(source)
    if args.output:
        replace_section(args.output, tables, "uint16", provenance)
        print(f"Updated {args.output} with {Z_MAX * len(OUT_T_keV) * len(OUT_KAPPA):,} uint16 coefficients", file=sys.stderr)
    else:
        write_section(sys.stdout, tables, "uint16", provenance)

    if not args.no_report:
        report = make_report(args.lock, source, tables, source_phi_errors, output_phi_errors)
        write_report(args.report_json, json.dumps(report, indent=2, sort_keys=True) + "\n")
        write_report(args.report_markdown, markdown_report(report))
        print(f"Wrote {args.report_json} and {args.report_markdown}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
