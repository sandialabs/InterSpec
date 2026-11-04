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

"""Validate the source multiple-scatter model over whole transport cells.

Truth is the union of exactly two Compton scatters and three-or-more Compton
scatters. Rayleigh remains active in transport but does not increment the
interaction order. Source-electron bremsstrahlung and pair-annihilation
photons are excluded.
"""

from __future__ import annotations

import argparse
import csv
import bank_io
import hashlib
import json
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from validate_single_scatter import calculate_metrics, cell_name


DEFAULT_ENERGIES = [100.0, 300.0, 661.657, 1332.5]
DEFAULT_Z = [6, 26, 50, 82]
DEFAULT_TAU = [0.2, 1.0, 4.0]

# The multiple-scatter surrogate's banks (`--grid`).  Training: twelve anchor
# energies plus the 30 interior points of a 32-point log grid over
# 10 keV-12 MeV, 20 elements, 7 optical depths (5880 cells).  Test: energies at
# the midpoints of a 31-interval log grid, elements and optical depths off the
# training grid (792 cells), so every test cell is a genuine off-grid query.
TRAINING_ANCHOR_ENERGIES = [
    10.0, 20.0, 40.0, 80.0, 150.0, 300.0, 511.0, 1000.0, 2000.0, 5000.0,
    8000.0, 12000.0]
TRAINING_Z = [1, 2, 4, 6, 8, 10, 13, 16, 20, 26, 29, 35, 42, 50, 56, 64,
              69, 74, 82, 92]
TRAINING_TAU = [0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
TEST_ENERGY_INDICES = [0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 28, 30]
TEST_Z = [2, 10, 16, 23, 32, 38, 46, 60, 69, 78, 87]
TEST_TAU = [0.3, 0.75, 1.5, 3.0, 6.0, 10.0]


def preset_grid(name: str) -> tuple[list[float], list[int], list[float]]:
    if name == "training":
        interior = [float(value) for value in np.geomspace(10.0, 12000.0, 32)[1:-1]]
        return (TRAINING_ANCHOR_ENERGIES + interior, TRAINING_Z, TRAINING_TAU)
    if name == "test":
        energies = [10.0 * 1200.0 ** ((k + 0.5) / 31.0)
                    for k in TEST_ENERGY_INDICES]
        return (energies, TEST_Z, TEST_TAU)
    raise ValueError(f"unknown grid preset: {name}")
TRUTH_COMPONENTS = (
    "double_compton",
    "double_compton_rayleigh",
    "higher_compton",
    "higher_compton_rayleigh",
)
FINE_TO_AGGREGATE = {
    "triple_compton": "higher_compton",
    "quadruple_compton": "higher_compton",
    "five_plus_compton": "higher_compton",
    "triple_compton_rayleigh": "higher_compton_rayleigh",
    "quadruple_compton_rayleigh": "higher_compton_rayleigh",
    "five_plus_compton_rayleigh": "higher_compton_rayleigh",
}


def stable_seed(energy: float, z_value: int, tau: float) -> int:
    identity = (
        f"multiple-scatter-compton-order:"
        f"{energy:.12g}:{z_value}:{tau:.12g}")
    digest = hashlib.sha256(identity.encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def run_cell(
    output: Path,
    energy: float,
    z_value: int,
    tau: float,
    bins: int,
    events: int,
    source_driver: Path,
    model_driver: Path,
    force: bool,
    force_model: bool,
    multiple_surrogate: Path | None,
    fine_interaction_orders: bool,
) -> tuple[str, str]:
    identity = cell_name(energy, z_value, tau)
    truth = bank_io.written_path(output / "cells", identity, "ceelo")
    model = bank_io.written_path(output / "cells", identity, "model")
    truth.parent.mkdir(parents=True, exist_ok=True)
    common = [
        "--energy-kev", str(energy), "--z", str(z_value),
        "--tau", str(tau), "--bins", str(bins),
    ]
    def is_rayleigh_enabled(path: Path) -> bool:
        path = bank_io.existing(path)
        if path is None:
            return False
        with bank_io.open_cell(path) as stream:
            first = next(csv.DictReader(stream), None)
        if first is None:
            return False
        value = first.get("rayleigh_enabled")
        return (
            value is None
            or value.strip().lower() in {"1", "true", "yes"})

    if force or not is_rayleigh_enabled(truth):
        truth_command = [
            str(source_driver), *common, "--events", str(events),
            "--seed", str(stable_seed(energy, z_value, tau)),
            "--interaction-orders"]
        if fine_interaction_orders:
            truth_command.append("--fine-interaction-orders")
        result = subprocess.run(
            truth_command,
            check=True, text=True, stdout=subprocess.PIPE)
        bank_io.write_cell(truth, result.stdout)
    if force or force_model or not is_rayleigh_enabled(model):
        model_command = [str(model_driver), *common]
        if multiple_surrogate is not None:
            model_command += [
                "--multiple-surrogate", str(multiple_surrogate)]
        result = subprocess.run(
            model_command,
            check=True, text=True, stdout=subprocess.PIPE)
        bank_io.write_cell(model, result.stdout)
    return identity, "ready"


def load_components(
    path: Path,
    bins: int,
    truth_file: bool,
) -> tuple[np.ndarray, np.ndarray, dict[str, float], dict[str, np.ndarray]]:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty cell: {path}")
    rayleigh_enabled = rows[0].get("rayleigh_enabled")
    if (rayleigh_enabled is not None
            and rayleigh_enabled.strip().lower()
            not in {"1", "true", "yes"}):
        raise ValueError(
            f"{path} is not a Rayleigh-enabled campaign cell")
    events = float(rows[0]["events"])
    energy = float(rows[0]["energy_keV"])
    wanted = TRUTH_COMPONENTS
    components = {name: np.zeros(bins) for name in wanted}
    component_variance = {name: np.zeros(bins) for name in wanted}
    fine_order_mode = any(
        row["component"] in FINE_TO_AGGREGATE
        and float(row["sum_w"]) != 0.0
        for row in rows)
    for row in rows:
        name = row["component"]
        if fine_order_mode:
            name = FINE_TO_AGGREGATE.get(name, name)
        if name not in components:
            continue
        lo = float(row["bin_lo_keV"])
        hi = float(row["bin_hi_keV"])
        if not hi > lo:
            continue
        index = min(bins - 1, int(round(lo / energy * bins)))
        components[name][index] += float(row["sum_w"]) / events
        component_variance[name][index] += (
            float(row["sum_w2"]) / (events * events))
    total = sum(components.values(), np.zeros(bins))
    variance = sum(component_variance.values(), np.zeros(bins))
    metadata = {
        "energy_keV": energy,
        "Z": float(rows[0]["Z"]),
        "tau": float(rows[0]["tau"]),
        "areal_density_g_cm2": float(rows[0]["areal_density_g_cm2"]),
        "events": events,
        "rayleigh_enabled": 1.0,
    }
    return total, variance, metadata, components


def write_plots(
    output: Path,
    records: list[dict],
    loaded: dict[str, tuple[
        np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]],
    worst: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    energies = sorted({record["energy_keV"] for record in records})
    z_values = sorted({record["Z"] for record in records})
    figure, axes = plt.subplots(
        2, len(energies), figsize=(4.2 * len(energies), 7.0),
        sharex="col", constrained_layout=True)
    if len(energies) == 1:
        axes = axes[:, None]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(z_values)))
    for column, energy in enumerate(energies):
        for color, z_value in zip(colors, z_values):
            cells = sorted(
                (record for record in records
                 if record["energy_keV"] == energy
                 and record["Z"] == z_value
                 and math.isfinite(record["yield_relative_error"])
                 and math.isfinite(record["shape_total_variation"])),
                key=lambda record: record["tau"])
            axes[0, column].plot(
                [cell["tau"] for cell in cells],
                [100.0 * cell["yield_relative_error"] for cell in cells],
                marker="o", color=color, label=f"Z={z_value}")
            axes[1, column].plot(
                [cell["tau"] for cell in cells],
                [100.0 * cell["shape_total_variation"] for cell in cells],
                marker="o", color=color)
        axes[0, column].axhspan(-2.0, 2.0, color="#59a14f", alpha=0.12)
        axes[0, column].axhline(0.0, color="#303030", linewidth=0.8)
        axes[0, column].set_title(f"{energy:g} keV")
        axes[0, column].set_ylabel("Multiple-scatter yield error (%)")
        axes[1, column].set_ylabel("Normalized-shape TV distance (%)")
        axes[1, column].set_xlabel("Process optical depth")
        for axis in axes[:, column]:
            axis.set_xscale("log")
            axis.grid(alpha=0.2)
    axes[0, -1].legend(fontsize=8)
    figure.suptitle(
        "Aggregate multiple-scatter validation across E, Z, and shield")
    figure.savefig(output / "multiple_scatter_grid_summary.png", dpi=180)
    plt.close(figure)

    ranked = sorted(
        (record for record in records
         if math.isfinite(record["yield_relative_error"])
         and math.isfinite(record["shape_total_variation"])),
        key=lambda record: max(
            abs(record["yield_relative_error"]),
            record["shape_total_variation"]),
        reverse=True)
    with PdfPages(output / "multiple_scatter_worst_cells.pdf") as pdf:
        for record in ranked[:worst]:
            identity = record["cell_id"]
            truth, variance, model, components = loaded[identity]
            energy = record["energy_keV"]
            edges = np.linspace(0.0, energy, truth.size + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            sigma = np.sqrt(np.maximum(variance, 0.0))
            figure, (shape_axis, ratio_axis) = plt.subplots(
                2, 1, figsize=(11.0, 7.0), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]},
                constrained_layout=True)
            shape_axis.stairs(
                truth, edges, color="#202020", linewidth=1.4,
                label="CeeLo aggregate multiple scatter")
            colors = (
                "#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756")
            for name, color in zip(TRUTH_COMPONENTS, colors):
                shape_axis.stairs(
                    components[name], edges, color=color, linewidth=0.9,
                    label=name.replace("_", " "))
            shape_axis.stairs(
                model, edges, color="#d1495b", linewidth=1.2,
                label="functional model")
            shape_axis.fill_between(
                centers, np.maximum(0.0, truth - sigma), truth + sigma,
                step="mid", color="#4c78a8", alpha=0.15,
                label="MC uncertainty")
            shape_axis.set_yscale("log")
            shape_axis.set_ylim(
                max(1e-12, truth.max() * 1e-6),
                max(truth.max(), model.max()) * 2.0)
            shape_axis.set_ylabel("Yield / emitted photon / bin")
            shape_axis.legend(fontsize=8, ncol=2)
            shape_axis.grid(alpha=0.2)

            valid = truth > 0.0
            ratio = np.full_like(truth, np.nan)
            ratio[valid] = model[valid] / truth[valid]
            ratio_axis.plot(centers, ratio, color="#d1495b", linewidth=0.9)
            ratio_axis.axhline(1.0, color="#202020", linewidth=0.8)
            ratio_axis.set_ylim(0.0, 3.0)
            ratio_axis.set_ylabel("model / truth")
            ratio_axis.set_xlabel("Escaped photon energy (keV)")
            ratio_axis.grid(alpha=0.2)
            figure.suptitle(
                f"{energy:g} keV, Z={record['Z']}, "
                f"tau={record['tau']:g}; yield error "
                f"{100.0 * record['yield_relative_error']:+.1f}%, "
                f"shape TV {100.0 * record['shape_total_variation']:.1f}%")
            pdf.savefig(figure)
            plt.close(figure)


def parse_cells(args: argparse.Namespace) -> list[tuple[float, int, float]]:
    if args.grid is not None and not args.cell:
        energies, z_values, tau_values = preset_grid(args.grid)
        return [
            (energy, z_value, tau)
            for energy in energies
            for z_value in z_values
            for tau in tau_values
        ]
    if args.cell:
        cells = []
        for text in args.cell:
            fields = text.split(",")
            if len(fields) != 3:
                raise ValueError(
                    f"--cell must be ENERGY,Z,TAU, received: {text}")
            cells.append((float(fields[0]), int(fields[1]), float(fields[2])))
        return cells
    return [
        (energy, z_value, tau)
        for energy in args.energies
        for z_value in args.z_values
        for tau in args.tau_values
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action", choices=["run", "analyze", "all"], default="all", nargs="?")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--energies", type=float, nargs="+",
                        default=DEFAULT_ENERGIES)
    parser.add_argument("--z-values", type=int, nargs="+", default=DEFAULT_Z)
    parser.add_argument("--tau-values", type=float, nargs="+",
                        default=DEFAULT_TAU)
    parser.add_argument("--cell", action="append", default=[],
                        metavar="ENERGY,Z,TAU")
    parser.add_argument(
        "--grid", choices=("training", "test"),
        help="the surrogate's training or test bank grid instead of "
             "--energies/--z-values/--tau-values")
    parser.add_argument("--bins", type=int, default=256)
    parser.add_argument("--metric-bins", type=int, default=64)
    parser.add_argument("--events", type=int, default=1000000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--worst-plots", type=int, default=12)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--force-model", action="store_true",
        help="regenerate only functional-model cells, preserving MC truth")
    parser.add_argument(
        "--fine-interaction-orders", action="store_true",
        help="split MC truth into exactly 3, exactly 4, and 5+ Compton groups")
    parser.add_argument(
        "--multiple-surrogate", type=Path,
        help="functional high-order residual payload")
    parser.add_argument(
        "--source-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_escape_campaign"))
    parser.add_argument(
        "--model-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_model_evaluate"))
    args = parser.parse_args()
    cells = parse_cells(args)

    if args.action in {"run", "all"}:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [
                pool.submit(
                    run_cell, args.output, *cell, args.bins, args.events,
                    args.source_driver, args.model_driver,
                    args.force, args.force_model,
                    args.multiple_surrogate,
                    args.fine_interaction_orders)
                for cell in cells
            ]
            for future in as_completed(futures):
                identity, state = future.result()
                print(state, identity)
        if args.action == "run":
            return 0

    records = []
    loaded = {}
    for energy, z_value, tau in cells:
        identity = cell_name(energy, z_value, tau)
        truth, variance, metadata, components = load_components(
            bank_io.existing(bank_io.written_path(args.output / "cells", identity, "ceelo")),
            args.bins, True)
        model, _, _, _ = load_components(
            bank_io.existing(bank_io.written_path(args.output / "cells", identity, "model")),
            args.bins, False)
        record = calculate_metrics(
            truth, variance, model, metadata, args.metric_bins)
        record["cell_id"] = identity
        for name, values in components.items():
            record[f"{name}_yield"] = float(values.sum())
        records.append(record)
        loaded[identity] = (truth, variance, model, components)

    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "multiple_scatter_metrics.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    finite = [
        record for record in records
        if math.isfinite(record["weighted_rms_relative_error"])]
    report = {
        "schema": 1,
        "cells": len(records),
        "statistically_evaluable_cells": len(finite),
        "energies_keV": sorted({cell[0] for cell in cells}),
        "Z_values": sorted({cell[1] for cell in cells}),
        "tau_values": sorted({cell[2] for cell in cells}),
        "events_per_cell": args.events,
        "bins": args.bins,
        "metric_bins": args.metric_bins,
        "strict_truth_label": (
            "two-or-more Compton, with Rayleigh active but excluded from "
            "the interaction-order count; "
            "pair and brems photons excluded"),
        "source_rayleigh": True,
        "optical_depth_basis": "mu_pe + mu_cs + mu_rs + mu_pp",
        "yield_error_median_abs": float(np.median([
            abs(record["yield_relative_error"]) for record in finite])),
        "yield_error_p95_abs": float(np.percentile([
            abs(record["yield_relative_error"]) for record in finite], 95)),
        "shape_tv_median": float(np.median([
            record["shape_total_variation"] for record in finite])),
        "shape_tv_p95": float(np.percentile([
            record["shape_total_variation"] for record in finite], 95)),
        "weighted_rms_median": float(np.median([
            record["weighted_rms_relative_error"] for record in finite])),
        "weighted_rms_p95": float(np.percentile([
            record["weighted_rms_relative_error"] for record in finite], 95)),
        "bin_p95_error_median": float(np.median([
            record["p95_relative_error"] for record in finite])),
        "bin_p95_error_p95": float(np.percentile([
            record["p95_relative_error"] for record in finite], 95)),
        "maximum_absolute_yield_pull": max(
            abs(record["yield_pull"]) for record in finite),
    }
    (args.output / "multiple_scatter_report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    write_plots(args.output, records, loaded, args.worst_plots)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
