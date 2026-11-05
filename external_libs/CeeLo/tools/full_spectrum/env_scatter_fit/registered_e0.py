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

"""Generate and structurally recompute the registered-energy E0 gate."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Mapping


SCHEMA = "CEEFIT-REGISTERED-E0-GATE-v1"
GENERATOR_SCHEMA = "CEEFIT-REGISTERED-E0-GENERATOR-v1"
GRID_ENERGIES_KEV = (59.5, 122.0, 239.0, 478.0, 955.0, 1836.0, 2614.0)
MIDPOINT_ENERGIES_KEV = tuple(
    math.sqrt(left * right)
    for left, right in zip(GRID_ENERGIES_KEV, GRID_ENERGIES_KEV[1:])
)
_SCHEMES = (
    ("lin(s) vs lin(x)", False, False, False),
    ("lin(s) vs log(x)", True, False, False),
    ("log(s) vs lin(x)", False, True, False),
    ("log(s) vs log(x)", True, True, False),
    ("PCHIP s  vs lin(x)", False, False, True),
    ("PCHIP s  vs log(x)", True, False, True),
    ("PCHIP ln s vs lin(x)", False, True, True),
    ("PCHIP ln s vs log(x)", True, True, True),
)


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"registered-E0 {name} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"registered-E0 {name} must be finite")
    return result


def _edge_slope(h0: float, h1: float, m0: float, m1: float) -> float:
    slope = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    sign = lambda value: (value > 0.0) - (value < 0.0)
    if sign(slope) != sign(m0):
        return 0.0
    if sign(m0) != sign(m1) and abs(slope) > 3.0 * abs(m0):
        return 3.0 * m0
    return slope


def _pchip(xs: list[float], ys: list[float], query: float) -> float:
    """Dependency-free mirror of ceelo::Pchip for evidence recomputation."""

    if len(xs) != len(ys) or len(xs) < 2 or any(
        right <= left for left, right in zip(xs, xs[1:])
    ):
        raise ValueError("registered-E0 raw interpolation nodes are invalid")
    if query <= xs[0]:
        return ys[0]
    if query >= xs[-1]:
        return ys[-1]
    interval = 0
    while interval + 1 < len(xs) and xs[interval + 1] <= query:
        interval += 1
    widths = [right - left for left, right in zip(xs, xs[1:])]
    slopes = [
        (right - left) / width
        for left, right, width in zip(ys, ys[1:], widths)
    ]
    derivatives = [0.0] * len(xs)
    if len(xs) == 2:
        derivatives[:] = [slopes[0], slopes[0]]
    else:
        for index in range(1, len(xs) - 1):
            before, after = slopes[index - 1], slopes[index]
            if before == 0.0 or after == 0.0 or (before > 0.0) != (after > 0.0):
                derivatives[index] = 0.0
            else:
                w1 = 2.0 * widths[index] + widths[index - 1]
                w2 = widths[index] + 2.0 * widths[index - 1]
                derivatives[index] = (w1 + w2) / (w1 / before + w2 / after)
        derivatives[0] = _edge_slope(
            widths[0], widths[1], slopes[0], slopes[1]
        )
        derivatives[-1] = _edge_slope(
            widths[-1], widths[-2], slopes[-1], slopes[-2]
        )
    width = widths[interval]
    t = (query - xs[interval]) / width
    t2, t3 = t * t, t * t * t
    return (
        (2.0 * t3 - 3.0 * t2 + 1.0) * ys[interval]
        + (t3 - 2.0 * t2 + t) * width * derivatives[interval]
        + (-2.0 * t3 + 3.0 * t2) * ys[interval + 1]
        + (t3 - t2) * width * derivatives[interval + 1]
    )


def _interpolate(
    scheme: tuple[str, bool, bool, bool],
    xs: list[float],
    ys: list[float],
    query: float,
    *,
    force_linear_y: bool = False,
) -> float:
    _, log_x, log_y, pchip = scheme
    log_y = log_y and not force_linear_y
    transformed_x = [math.log(value) if log_x else value for value in xs]
    transformed_y = [
        math.log(max(value, 1.0e-300)) if log_y else value for value in ys
    ]
    transformed_query = math.log(query) if log_x else query
    if pchip:
        value = _pchip(transformed_x, transformed_y, transformed_query)
    else:
        index = 0
        while (
            index + 2 < len(transformed_x)
            and transformed_x[index + 1] < transformed_query
        ):
            index += 1
        fraction = (
            (transformed_query - transformed_x[index])
            / (transformed_x[index + 1] - transformed_x[index])
        )
        value = transformed_y[index] + fraction * (
            transformed_y[index + 1] - transformed_y[index]
        )
    return math.exp(value) if log_y else value


def _require_close(actual: Any, expected: float, name: str) -> float:
    value = _finite(actual, name)
    if not math.isclose(value, expected, rel_tol=2.0e-10, abs_tol=2.0e-13):
        raise ValueError(
            f"registered-E0 {name} differs from independently recomputed raw nodes"
        )
    return value


def validate_registered_e0_report(record: Mapping[str, Any]) -> bool:
    """Independently reconstruct every acceptance metric from raw node maps."""

    if record.get("schema") != SCHEMA or type(record.get("passed")) is not bool:
        raise ValueError("registered-E0 report schema/pass flag is invalid")
    gate = record.get("gate_configuration")
    if not isinstance(gate, Mapping):
        raise ValueError("registered-E0 report lacks its gate configuration")
    geometry_tier = gate.get("geometry_tier")
    material = gate.get("material")
    if (
        geometry_tier not in {"air", "field"}
        or (geometry_tier == "air" and material != "Air")
        or (geometry_tier == "field" and material not in {"Soil", "Concrete"})
        or _finite(gate.get("distance_cm"), "gate distance") != 500.0
        or _finite(gate.get("source_height_cm"), "gate source height") != 100.0
        or _finite(gate.get("detector_height_cm"), "gate detector height") != 100.0
        or _finite(gate.get("scoring_radius_cm"), "gate scoring radius") != 4.0
        or _finite(gate.get("energy_cut_keV"), "gate energy cutoff") != 10.0
    ):
        raise ValueError("registered-E0 gate configuration is noncanonical")
    grid = record.get("grid")
    if not isinstance(grid, Mapping):
        raise ValueError("registered-E0 report lacks its grid definition")
    expected_registered = [index / 12 for index in range(13)]
    expected_mu = [-1.0 + index / 4 for index in range(9)]
    if (
        grid.get("native_registered_energy_bins") != 384
        or grid.get("native_mu_bins") != 64
        or grid.get("registered_energy_bands") != 12
        or grid.get("mu_bands") != 8
        or grid.get("registered_energy_edges") != expected_registered
        or grid.get("mu_edges") != expected_mu
    ):
        raise ValueError("registered-E0 report does not use the canonical 384x64 -> 12x8 grid")

    targets = record.get("stopping_targets")
    if not isinstance(targets, Mapping):
        raise ValueError("registered-E0 report lacks stopping targets")
    total_target = _finite(targets.get("total_relative_se"), "total_relative_se")
    absolute_target = _finite(
        targets.get("absolute_se_fraction_of_total"),
        "absolute_se_fraction_of_total",
    )
    significant_fraction = _finite(
        targets.get("significant_fraction_of_total"),
        "significant_fraction_of_total",
    )
    if not (
        0.0 < total_target <= 0.0025
        and 0.0 < absolute_target <= 0.001
        and 0.0 <= significant_fraction <= 0.01
        and isinstance(targets.get("batch_histories"), int)
        and targets["batch_histories"] > 0
        and isinstance(targets.get("max_histories_per_node"), int)
        and targets["max_histories_per_node"] >= targets["batch_histories"]
    ):
        raise ValueError("registered-E0 stopping targets are noncanonical")

    thresholds = record.get("acceptance_thresholds")
    expected_thresholds = {
        "max_integrated_relative_error": 0.01,
        "max_joint_tv": 0.01,
        "max_energy_marginal_tv": 0.005,
        "max_angular_marginal_tv": 0.005,
    }
    if not isinstance(thresholds, Mapping) or any(
        _finite(thresholds.get(name), f"acceptance {name}") != value
        for name, value in expected_thresholds.items()
    ):
        raise ValueError("registered-E0 acceptance thresholds are noncanonical")
    if record.get("observable") != "total_continuum; exact lines remain atoms":
        raise ValueError("registered-E0 observable definition is noncanonical")

    nodes = record.get("nodes")
    if not isinstance(nodes, list) or len(nodes) != 13:
        raise ValueError("registered-E0 report must contain seven grid and six held-out nodes")
    roles = [node.get("role") for node in nodes if isinstance(node, Mapping)]
    if roles.count("grid") != 7 or roles.count("heldout") != 6:
        raise ValueError("registered-E0 node role inventory is incomplete")
    node_targets_met = True
    raw_nodes: dict[str, list[dict[str, Any]]] = {"grid": [], "heldout": []}
    for index, node in enumerate(nodes):
        if not isinstance(node, Mapping):
            raise ValueError("registered-E0 node record is malformed")
        role = str(node["role"])
        histories = node.get("histories")
        energy = _finite(node.get("energy_keV"), f"nodes[{index}].energy_keV")
        total_mass = _finite(
            node.get("total_continuum_s_hat"),
            f"nodes[{index}].total_continuum_s_hat",
        )
        joint = node.get("joint_cell_masses")
        if (
            total_mass <= 0.0
            or not isinstance(joint, list)
            or len(joint) != 96
        ):
            raise ValueError("registered-E0 raw node maps are incomplete")
        joint_values = [
            _finite(value, f"nodes[{index}].joint_cell_masses[{cell}]")
            for cell, value in enumerate(joint)
        ]
        if any(value < 0.0 for value in joint_values):
            raise ValueError("registered-E0 raw node maps contain negative fluence")
        total_se = _finite(node.get("total_relative_se"), f"nodes[{index}].total_relative_se")
        absolute_se = _finite(
            node.get("max_significant_absolute_se_fraction"),
            f"nodes[{index}].max_significant_absolute_se_fraction",
        )
        flags = (
            node.get("total_target_met"),
            node.get("absolute_target_met"),
            node.get("relative_band_target_met"),
        )
        if any(type(flag) is not bool for flag in flags):
            raise ValueError("registered-E0 node stopping flags must be boolean")
        enough_rounds = bool(
            isinstance(histories, int)
            and histories >= 12 * targets["batch_histories"]
            and histories <= targets["max_histories_per_node"]
            and histories % targets["batch_histories"] == 0
        )
        expected_flags = (
            enough_rounds and 0.0 <= total_se <= total_target,
            enough_rounds and 0.0 <= absolute_se <= absolute_target,
            enough_rounds,
        )
        if tuple(flags) != expected_flags:
            raise ValueError("registered-E0 node stopping flags differ from raw targets")
        derived = bool(enough_rounds and all(expected_flags))
        node_targets_met = node_targets_met and derived
        raw_nodes[role].append(
            {
                "energy": energy,
                "total": total_mass,
                "joint": joint_values,
            }
        )
    if record.get("all_node_stopping_targets_met") is not node_targets_met:
        raise ValueError("registered-E0 stopping summary differs from detailed nodes")

    raw_nodes["grid"].sort(key=lambda item: item["energy"])
    raw_nodes["heldout"].sort(key=lambda item: item["energy"])
    for supplied, expected, role in (
        ([item["energy"] for item in raw_nodes["grid"]], GRID_ENERGIES_KEV, "grid"),
        ([item["energy"] for item in raw_nodes["heldout"]], MIDPOINT_ENERGIES_KEV, "heldout"),
    ):
        if len(supplied) != len(expected) or any(
            not math.isclose(left, right, rel_tol=1.0e-13, abs_tol=1.0e-12)
            for left, right in zip(supplied, expected)
        ):
            raise ValueError(f"registered-E0 {role} energy ladder is noncanonical")

    midpoints = record.get("midpoints")
    if not isinstance(midpoints, list) or len(midpoints) != 6:
        raise ValueError("registered-E0 report must contain all six held-out midpoint metrics")
    metric_fields = (
        "integrated_relative_error",
        "joint_tv",
        "energy_marginal_tv",
        "angular_marginal_tv",
    )
    supplied_midpoints: dict[float, Mapping[str, Any]] = {}
    for index, midpoint in enumerate(midpoints):
        if not isinstance(midpoint, Mapping):
            raise ValueError("registered-E0 midpoint record is malformed")
        energy = _finite(midpoint.get("energy_keV"), f"midpoints[{index}].energy_keV")
        if energy in supplied_midpoints:
            raise ValueError("registered-E0 midpoint energies are repeated")
        supplied_midpoints[energy] = midpoint
        for name in metric_fields:
            value = _finite(midpoint.get(name), f"midpoints[{index}].{name}")
            if value < 0.0:
                raise ValueError("registered-E0 errors cannot be negative")

    xs = [item["energy"] for item in raw_nodes["grid"]]
    totals = [item["total"] for item in raw_nodes["grid"]]
    best_scheme = min(
        _SCHEMES,
        key=lambda scheme: max(
            abs(
                _interpolate(scheme, xs, totals, truth["energy"])
                / truth["total"]
                - 1.0
            )
            for truth in raw_nodes["heldout"]
        ),
    )
    if record.get("best_scheme") != best_scheme[0]:
        raise ValueError("registered-E0 best scheme differs from raw node totals")
    maxima = {name: 0.0 for name in metric_fields}
    matched_midpoint_ids: set[int] = set()
    for truth in raw_nodes["heldout"]:
        candidates = [
            value
            for energy, value in supplied_midpoints.items()
            if math.isclose(energy, truth["energy"], rel_tol=1.0e-13, abs_tol=1.0e-12)
        ]
        if len(candidates) != 1:
            raise ValueError("registered-E0 midpoint details do not match raw heldout nodes")
        supplied = candidates[0]
        matched_midpoint_ids.add(id(supplied))
        predicted = [
            max(
                0.0,
                _interpolate(
                    best_scheme,
                    xs,
                    [item["joint"][cell] for item in raw_nodes["grid"]],
                    truth["energy"],
                    force_linear_y=True,
                ),
            )
            for cell in range(96)
        ]
        predicted_total = sum(predicted)
        truth_total = sum(truth["joint"])
        predicted_integral = _interpolate(
            best_scheme, xs, totals, truth["energy"]
        )
        derived = {
            "integrated_relative_error": abs(
                predicted_integral - truth["total"]
            ) / max(truth["total"], 1.0e-300),
            "joint_tv": 1.0,
            "energy_marginal_tv": 1.0,
            "angular_marginal_tv": 1.0,
        }
        if predicted_total > 0.0 and truth_total > 0.0:
            predicted_energy = [sum(predicted[8 * e : 8 * (e + 1)]) for e in range(12)]
            truth_energy = [sum(truth["joint"][8 * e : 8 * (e + 1)]) for e in range(12)]
            predicted_mu = [sum(predicted[8 * e + m] for e in range(12)) for m in range(8)]
            truth_mu = [sum(truth["joint"][8 * e + m] for e in range(12)) for m in range(8)]
            derived.update(
                {
                    "joint_tv": 0.5 * sum(
                        abs(left / predicted_total - right / truth_total)
                        for left, right in zip(predicted, truth["joint"])
                    ),
                    "energy_marginal_tv": 0.5 * sum(
                        abs(left / predicted_total - right / truth_total)
                        for left, right in zip(predicted_energy, truth_energy)
                    ),
                    "angular_marginal_tv": 0.5 * sum(
                        abs(left / predicted_total - right / truth_total)
                        for left, right in zip(predicted_mu, truth_mu)
                    ),
                }
            )
        for name, expected in derived.items():
            _require_close(supplied.get(name), expected, f"midpoint {truth['energy']} {name}")
            maxima[name] = max(maxima[name], expected)
        _require_close(
            supplied.get("continuum_fraction_of_scattered"),
            truth_total / max(truth["total"], 1.0e-300),
            f"midpoint {truth['energy']} continuum coverage",
        )
    if len(matched_midpoint_ids) != 6:
        raise ValueError("registered-E0 midpoint inventory is incomplete")
    summary_names = {
        "integrated_relative_error": "max_integrated_relative_error",
        "joint_tv": "max_joint_tv",
        "energy_marginal_tv": "max_energy_marginal_tv",
        "angular_marginal_tv": "max_angular_marginal_tv",
    }
    for detail, summary in summary_names.items():
        _require_close(record.get(summary), maxima[detail], summary)

    computed = bool(
        record.get("predefined_index_bands") is True
        and record.get("statistically_empty_bands_included") is True
        and _finite(record.get("shape_coverage_fraction"), "shape_coverage_fraction") == 1.0
        and node_targets_met
        and maxima["integrated_relative_error"] <= 0.01
        and maxima["joint_tv"] <= 0.01
        and maxima["energy_marginal_tv"] <= 0.005
        and maxima["angular_marginal_tv"] <= 0.005
    )
    if record["passed"] != computed:
        raise ValueError("registered-E0 passed flag differs from recomputed detailed gates")
    return computed


def generate_registered_e0_report(
    campaign: Mapping[str, Any],
    executable_path: str | Path,
    output_path: str | Path,
    *,
    batch_histories: int = 1_000_000,
    max_histories: int = 128_000_000,
    threads: int = 1,
) -> dict[str, Any]:
    """Run the registered-axis producer and record its exact invocation."""

    executable = Path(executable_path).resolve()
    output = Path(output_path).resolve()
    if not executable.is_file():
        raise ValueError("registered-E0 executable does not exist")
    if batch_histories <= 0 or max_histories < batch_histories or threads <= 0:
        raise ValueError("registered-E0 history/thread parameters are invalid")
    if output == executable:
        raise ValueError("registered-E0 output must not overwrite the executable")
    output.parent.mkdir(parents=True, exist_ok=True)
    config = campaign["config"]
    physics = config["physics"]
    material_argument = str(config["material"]).lower()
    if config["mode"] == "air" and material_argument not in {"soil", "concrete"}:
        # The producer retains a field-material parser even though its
        # air branch does not instantiate a slab. Keep that inert argument
        # canonical and record the real campaign material in input_parameters.
        material_argument = "soil"
    with tempfile.NamedTemporaryFile(
        prefix=output.name + ".raw.", suffix=".json", dir=output.parent, delete=False
    ) as temporary:
        raw_path = Path(temporary.name)
    argv = [
        str(executable), "--axis", "e0", "--escale", "landmark",
        "--mode", str(config["mode"]), "--material", material_argument,
        "--precision", "0.0025", "--absolute-precision", "0.001",
        "--absolute-min-frac", "0.01", "--native-energy-bins", "384",
        "--native-mu-bins", "64", "--e-cut", str(float(physics["energy_cut_keV"])),
        "--scoring-radius", str(float(physics["scoring_radius_cm"])),
        "--batch", str(batch_histories), "--max-histories", str(max_histories),
        "--threads", str(threads), "--config-sha256", str(config["sha256"]),
        "--json-report", str(raw_path),
    ]
    try:
        subprocess.run(
            argv,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        report = json.loads(raw_path.read_text(encoding="utf-8"))
    finally:
        raw_path.unlink(missing_ok=True)
    report["generator"] = {
        "schema": GENERATOR_SCHEMA,
        "argv": argv,
        "input_parameters": {
            "config_sha256": config["sha256"],
            "canonical_configuration_json": config["canonical_json"],
            "canonical_configuration": json.loads(config["canonical_json"]),
            "gate_configuration": report["gate_configuration"],
            "batch_histories": batch_histories,
            "max_histories_per_node": max_histories,
            "threads": threads,
        },
        "output_parameters": {
            "observable": "total continuum on 12x8 registered-energy/mu bands",
            "acceptance_metrics": [
                "integrated_relative_error", "joint_tv",
                "energy_marginal_tv", "angular_marginal_tv",
            ],
        },
    }
    validate_registered_e0_report(report)
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    temporary_output = output.with_name(output.name + f".tmp.{os.getpid()}")
    with temporary_output.open("w", encoding="utf-8") as destination:
        destination.write(encoded)
        destination.flush()
        os.fsync(destination.fileno())
    os.replace(temporary_output, output)
    return report
