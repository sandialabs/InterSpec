#!/usr/bin/env python3
"""Generate K/L relaxation fixtures directly from locked EPICS2023 EADL."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

from epics_endf import EndfFile, control, endf_float, fields
from fetch_sources import load_lock, sha256_file
from generation_utils import as_float32, format_float32

CEELO_LICENCE_HEADER = """\
/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC.
 This library is distributed under the GNU Lesser General Public License,
 version 2.1 or (at your option) any later version.
 */
"""



HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DEFAULT_LOCK = HERE / "sources.lock.json"
DEFAULT_OUTPUT = ROOT / "src" / "cross_sections" / "relaxation_epics_data.cpp"
DEFAULT_REPORT = HERE / "reports" / "relaxation_eadl.json"
MAX_LINES = 24
MAX_RELAXATION_Z = 99
L_TRANSPORT_CUT_EV = 10_000.0


@dataclass(frozen=True)
class Transition:
    filling_shell: int
    emitted_shell: int
    energy_eV: float
    probability: float


@dataclass(frozen=True)
class Subshell:
    designator: int
    binding_eV: float
    electrons: float
    transitions: tuple[Transition, ...]


def parse_relaxation(endf: EndfFile, z: int) -> dict[int, Subshell]:
    lines = endf.section(z, 28, 533)
    number_subshells = control(lines[0])[4]
    position = 1
    result: dict[int, Subshell] = {}
    for _ in range(number_subshells):
        c1, _, _, _, npl, number_transitions = control(lines[position])
        designator = int(c1)
        position += 1
        data_lines = (npl + 5) // 6
        values: list[float] = []
        for line in lines[position:position + data_lines]:
            values.extend(endf_float(field) for field in fields(line))
        position += data_lines
        if npl != 6 * (number_transitions + 1):
            raise ValueError(
                f"Z={z} shell={designator}: NPL={npl} does not encode "
                f"{number_transitions} six-field transitions"
            )
        transitions = tuple(
            Transition(
                int(values[6 + 6 * index]),
                int(values[7 + 6 * index]),
                values[8 + 6 * index],
                values[9 + 6 * index],
            )
            for index in range(number_transitions)
        )
        if any(t.energy_eV < 0.0 or t.probability < 0.0 for t in transitions):
            raise ValueError(f"Z={z} shell={designator}: invalid transition")
        probability_sum = sum(t.probability for t in transitions)
        if number_transitions > 0 and abs(probability_sum - 1.0) > 2.0e-4:
            raise ValueError(
                f"Z={z} shell={designator}: transition probabilities sum to {probability_sum}"
            )
        result[designator] = Subshell(designator, values[0], values[1], transitions)
    if len(result) != number_subshells:
        raise ValueError(f"Z={z}: duplicate/missing EADL subshell blocks")
    return result


def radiative_lines(shell: Subshell, minimum_energy_eV: float = 0.0) -> list[Transition]:
    lines = [
        transition for transition in shell.transitions
        if transition.emitted_shell == 0
        and transition.energy_eV >= minimum_energy_eV
        and transition.probability > 0.0
    ]
    if len(lines) > MAX_LINES:
        raise ValueError(
            f"shell {shell.designator}: {len(lines)} radiative lines exceeds storage {MAX_LINES}"
        )
    return sorted(lines, key=lambda transition: transition.probability, reverse=True)


def line_statistics(lines: list[Transition]) -> tuple[float, float]:
    emitted_yield = sum(line.probability for line in lines)
    mean_eV = (
        sum(line.energy_eV * line.probability for line in lines) / emitted_yield
        if emitted_yield > 0.0 else 0.0
    )
    return emitted_yield, mean_eV


def line_payload(
    shell: Subshell | None,
    cut_eV: float,
    prefix: str,
    arrays: list[tuple[str, list[float]]],
) -> tuple[str, str, int, dict]:
    if shell is None:
        return (
            "nullptr", "nullptr", 0,
            {"yield": 0.0, "mean_energy_keV": 0.0, "line_count": 0},
        )
    lines = radiative_lines(shell, cut_eV)
    emitted_yield, mean_eV = line_statistics(lines)
    probabilities = [line.probability / emitted_yield for line in lines] if emitted_yield else []
    if lines:
        energy_name = prefix + "_energy_keV"
        probability_name = prefix + "_probability"
        arrays.append((energy_name, [line.energy_eV / 1000.0 for line in lines]))
        arrays.append((probability_name, probabilities))
    else:
        energy_name = probability_name = "nullptr"
    return energy_name, probability_name, len(lines), {
        "yield": emitted_yield,
        "mean_energy_keV": mean_eV / 1000.0,
        "line_count": len(lines),
    }


def fluorescence_initializer(
    shell: Subshell | None,
    cut_eV: float,
    prefix: str,
    arrays: list[tuple[str, list[float]]],
) -> tuple[str, dict]:
    energy_name, probability_name, count, stats = line_payload(
        shell, cut_eV, prefix, arrays
    )
    binding = shell.binding_eV / 1000.0 if shell else 0.0
    initializer = (
        "{" + format_float32(binding) + ", " + format_float32(stats["yield"])
        + f", {count}, {energy_name}, {probability_name}" + "}"
    )
    return initializer, stats


def coster_kronig(shell: Subshell | None, destination: int) -> float:
    if shell is None:
        return 0.0
    # A Coster-Kronig transition fills the initial L vacancy from a higher L
    # subshell and ejects an outer electron, leaving the destination L vacancy.
    return sum(
        transition.probability for transition in shell.transitions
        if transition.emitted_shell != 0 and transition.filling_shell == destination
    )


def l_initializer(
    subshells: dict[int, Subshell],
    z: int,
    arrays: list[tuple[str, list[float]]],
) -> tuple[str, list[dict]]:
    l1, l2, l3 = (subshells.get(index) for index in (2, 3, 4))
    ck = (
        coster_kronig(l1, 3),
        coster_kronig(l1, 4),
        coster_kronig(l2, 4),
    )
    sub_initializers = []
    stats = []
    for label, shell in zip(("L1", "L2", "L3"), (l1, l2, l3)):
        energy_name, probability_name, count, stat = line_payload(
            shell, L_TRANSPORT_CUT_EV, f"Z{z}_{label}", arrays
        )
        sub_initializers.append(
            "{" + format_float32(stat["yield"]) + f", {count}, "
            + energy_name + ", " + probability_name + "}"
        )
        stats.append(stat)
    l3_edge = l3.binding_eV / 1000.0 if l3 else 0.0
    initializer = (
        "{" + format_float32(l3_edge) + ", {"
        + ", ".join(format_float32(value) for value in ck) + "}, {"
        + ", ".join(sub_initializers) + "}}"
    )
    return initializer, stats


def represented_statistics(initializer_stats: dict) -> tuple[float, float]:
    return initializer_stats["yield"], initializer_stats["mean_energy_keV"]


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eadl", type=Path, required=True)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    lock = load_lock(args.lock)
    source = lock["sources"]["eadl2023"]
    actual = sha256_file(args.eadl)
    if actual != source["sha256"]:
        raise ValueError(f"EADL SHA-256 mismatch: expected {source['sha256']}, got {actual}")
    endf = EndfFile(args.eadl)

    k_initializers = ["{0.0f, 0.0f, 0, nullptr, nullptr}" for _ in range(100)]
    l_empty = "{0.0f, {0.0f, 0.0f, 0.0f}, {{}, {}, {}}}"
    l_initializers = [l_empty for _ in range(100)]
    line_arrays: list[tuple[str, list[float]]] = []
    all_stats = []
    yield_errors = []
    mean_errors = []
    maximum_yield_error = 0.0
    maximum_mean_error = 0.0

    # Photon/electron transport is limited to Z<=92, but radioactive-decay
    # daughters queried by the cascade model extend through Z=99 (for example,
    # Am-241 -> Np-237, Z=93). EPICS2023 EADL has MF=28/MT=533 sections through
    # Z=99, so retain that complete relaxation-only domain.
    for z in range(1, MAX_RELAXATION_Z + 1):
        subshells = parse_relaxation(endf, z)
        k_initializer, k_stats = fluorescence_initializer(
            subshells.get(1), 0.0, f"Z{z}_K", line_arrays
        )
        k_initializers[z] = k_initializer
        l_initializers[z], l_stats = l_initializer(subshells, z, line_arrays)
        for label, stats in [("K", k_stats)] + [
            (name, stats) for name, stats in zip(("L1", "L2", "L3"), l_stats)
        ]:
            # Validate the exact float32 payload, including normalized line probabilities.
            shell = subshells.get({"K": 1, "L1": 2, "L2": 3, "L3": 4}[label])
            cut = 0.0 if label == "K" else L_TRANSPORT_CUT_EV
            lines = radiative_lines(shell, cut) if shell else []
            source_yield, source_mean_eV = line_statistics(lines)
            stored_yield = as_float32(source_yield)
            stored_prob = [as_float32(line.probability / source_yield) for line in lines] if source_yield else []
            stored_energy = [as_float32(line.energy_eV / 1000.0) for line in lines]
            norm = sum(stored_prob)
            stored_mean_keV = (
                sum(e * p for e, p in zip(stored_energy, stored_prob)) / norm if norm else 0.0
            )
            yield_error = abs(stored_yield / source_yield - 1.0) if source_yield else 0.0
            mean_error = abs(stored_mean_keV / (source_mean_eV / 1000.0) - 1.0) if source_mean_eV else 0.0
            maximum_yield_error = max(maximum_yield_error, yield_error)
            maximum_mean_error = max(maximum_mean_error, mean_error)
            yield_errors.append(yield_error)
            mean_errors.append(mean_error)
            all_stats.append({"Z": z, "shell": label, **stats})

    if maximum_yield_error > 1.0e-3 or maximum_mean_error > 1.0e-3:
        raise ValueError(
            f"relaxation float32 gates failed: yield={maximum_yield_error:.6%}, "
            f"mean={maximum_mean_error:.6%}"
        )

    with args.output.open("w", encoding="utf-8") as out:
        out.write(CEELO_LICENCE_HEADER)
        out.write("""/* AUTO-GENERATED by tools/prepare_cross_sections/generate_eadl_relaxation.py.
 Direct source: EPICS2023 EADL MF=28 MT=533, verified by sources.lock.json.
 No Geant4/G4EMLOW data is read or emitted. */
#include "cross_sections/relaxation_epics_data.h"
namespace ceelo {
""")
        for name, values in line_arrays:
            out.write(f"static const float {name}[{len(values)}] = {{")
            out.write(", ".join(format_float32(value) for value in values))
            out.write("};\n")
        out.write("\nextern const FluorescenceData g_epics_k_relaxation[kRelaxationMaxZ + 1] = {\n")
        for z, initializer in enumerate(k_initializers):
            out.write(f"    {initializer}{',' if z < 99 else ''} // Z={z}\n")
        out.write("};\nextern const LFluorescenceData g_epics_l_relaxation[kRelaxationMaxZ + 1] = {\n")
        for z, initializer in enumerate(l_initializers):
            out.write(f"    {initializer}{',' if z < 99 else ''} // Z={z}\n")
        out.write("};\n} // namespace ceelo\n")

    k_stats = [stats for stats in all_stats if stats["shell"] == "K"]
    l_stats = [stats for stats in all_stats if stats["shell"] != "K"]
    relaxation_slots = MAX_RELAXATION_Z + 1
    k_native_coefficients = sum(2 + 2 * stats["line_count"] for stats in k_stats)
    l_native_coefficients = (
        MAX_RELAXATION_Z * 7
        + 2 * sum(stats["line_count"] for stats in l_stats)
    )
    report = {
        "schema_version": 1,
        "source": "EPICS2023 EADL MF=28 MT=533",
        "source_sha256": source["sha256"],
        "source_lock_sha256": hashlib.sha256(args.lock.read_bytes()).hexdigest(),
        "maximum_float32_yield_error_percent": maximum_yield_error * 100.0,
        "maximum_float32_mean_emitted_energy_error_percent": maximum_mean_error * 100.0,
        "p95_float32_yield_error_percent": percentile(yield_errors, 0.95) * 100.0,
        "p95_float32_mean_emitted_energy_error_percent": percentile(mean_errors, 0.95) * 100.0,
        "acceptance_gate_percent": 0.1,
        "passed": maximum_yield_error <= 1.0e-3 and maximum_mean_error <= 1.0e-3,
        "maximum_lines_per_subshell": max(stats["line_count"] for stats in all_stats),
        "maximum_atomic_number": MAX_RELAXATION_Z,
        "l_transport_cut_keV": L_TRANSPORT_CUT_EV / 1000.0,
        "processes": {
            "K_relaxation": {
                "native_coefficients": k_native_coefficients,
                "retained_coefficients": 2 * relaxation_slots
                    + 2 * sum(stats["line_count"] for stats in k_stats),
                "coefficient_bytes": (2 * relaxation_slots
                    + 2 * sum(stats["line_count"] for stats in k_stats)) * 4,
                "metadata_bytes": relaxation_slots * 24,
            },
            "L_relaxation": {
                "native_coefficients": l_native_coefficients,
                "retained_coefficients": 7 * relaxation_slots
                    + 2 * sum(stats["line_count"] for stats in l_stats),
                "coefficient_bytes": (7 * relaxation_slots
                    + 2 * sum(stats["line_count"] for stats in l_stats)) * 4,
                "metadata_bytes": relaxation_slots * 72,
            },
        },
        "validation_definition": (
            "For every K/L subshell, compare the float32 runtime transported-x-ray "
            "yield and probability-weighted mean energy against the same direct EADL "
            "transitions after applying the documented 10 keV L-photon transport cut."
        ),
        "subshells": all_stats,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {args.output} and {args.report}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
