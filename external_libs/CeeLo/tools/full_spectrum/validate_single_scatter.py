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

"""Validate the physics single-Compton source model over whole transport cells."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


DEFAULT_ENERGIES = [100.0, 300.0, 661.657, 1332.5]
DEFAULT_Z = [6, 26, 50, 82]
DEFAULT_TAU = [0.2, 1.0, 4.0]


def cell_name(energy: float, z_value: int, tau: float) -> str:
    return f"E{energy:g}_Z{z_value}_tau{tau:g}".replace(".", "p")


def stable_seed(energy: float, z_value: int, tau: float) -> int:
    text = (
        f"single-scatter-compton-order:"
        f"{energy:.12g}:{z_value}:{tau:.12g}")
    digest = hashlib.sha256(text.encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def is_matching_physics_cell(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open(newline="") as stream:
        first = next(csv.DictReader(stream), None)
    if first is None:
        return False
    rayleigh = first.get("rayleigh_enabled")
    doppler = first.get("compton_doppler_enabled")
    return (
        rayleigh is not None
        and rayleigh.strip().lower() in {"1", "true", "yes"}
        and doppler is not None
        and doppler.strip().lower() in {"1", "true", "yes"})


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
) -> tuple[str, str]:
    identity = cell_name(energy, z_value, tau)
    truth = output / "cells" / f"{identity}_ceelo.csv"
    model = output / "cells" / f"{identity}_model.csv"
    truth.parent.mkdir(parents=True, exist_ok=True)
    common = [
        "--energy-kev", str(energy), "--z", str(z_value),
        "--tau", str(tau), "--bins", str(bins),
    ]
    if force or not is_matching_physics_cell(truth):
        result = subprocess.run(
            [str(source_driver), *common, "--events", str(events),
             "--seed", str(stable_seed(energy, z_value, tau))],
            check=True, text=True, stdout=subprocess.PIPE)
        truth.write_text(result.stdout)
    if force or force_model or not is_matching_physics_cell(model):
        result = subprocess.run(
            [str(model_driver), *common],
            check=True, text=True, stdout=subprocess.PIPE)
        model.write_text(result.stdout)
    return identity, "ready"


def load_single(path: Path, bins: int) -> tuple[
        np.ndarray, np.ndarray, dict[str, float]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty cell: {path}")
    rayleigh_enabled = rows[0].get("rayleigh_enabled")
    if (rayleigh_enabled is not None
            and rayleigh_enabled.strip().lower()
            not in {"1", "true", "yes"}):
        raise ValueError(
            f"{path} is not a Rayleigh-enabled source cell")
    doppler_enabled = rows[0].get("compton_doppler_enabled")
    if (doppler_enabled is None
            or doppler_enabled.strip().lower()
            not in {"1", "true", "yes"}):
        raise ValueError(
            f"{path} is not a Doppler-enabled source cell")
    events = float(rows[0]["events"])
    values = np.zeros(bins)
    variance = np.zeros(bins)
    energy = float(rows[0]["energy_keV"])
    for row in rows:
        if row["component"] not in {
                "single_compton", "single_compton_rayleigh"}:
            continue
        lo = float(row["bin_lo_keV"])
        hi = float(row["bin_hi_keV"])
        if not hi > lo:
            continue
        index = min(bins - 1, int(round(lo / energy * bins)))
        values[index] += float(row["sum_w"]) / events
        variance[index] += float(row["sum_w2"]) / (events * events)
    metadata = {
        "energy_keV": energy,
        "Z": float(rows[0]["Z"]),
        "tau": float(rows[0]["tau"]),
        "areal_density_g_cm2": float(rows[0]["areal_density_g_cm2"]),
        "events": events,
        "rayleigh_enabled": 1.0,
    }
    return values, variance, metadata


def calculate_metrics(
    truth: np.ndarray,
    variance: np.ndarray,
    model: np.ndarray,
    metadata: dict[str, float],
    metric_bins: int,
) -> dict[str, float | int]:
    if truth.size % metric_bins != 0:
        raise ValueError("metric bin count must divide simulation bin count")
    group = truth.size // metric_bins
    truth = truth.reshape(metric_bins, group).sum(axis=1)
    variance = variance.reshape(metric_bins, group).sum(axis=1)
    model = model.reshape(metric_bins, group).sum(axis=1)
    truth_yield = float(truth.sum())
    model_yield = float(model.sum())
    truth_yield_sigma = math.sqrt(float(variance.sum()))
    yield_relative_error = (
        (model_yield - truth_yield) / truth_yield
        if truth_yield > 0.0 else math.inf)
    if truth_yield > 0.0 and model_yield > 0.0:
        total_variation = 0.5 * float(np.abs(
            truth / truth_yield - model / model_yield).sum())
    else:
        total_variation = math.inf
    relative_uncertainty = np.sqrt(np.maximum(variance, 0.0)) \
        / np.maximum(truth, 1e-30)
    valid = (truth > 1e-12) & (relative_uncertainty <= 0.05)
    relative_error = np.abs(model - truth) / np.maximum(truth, 1e-30)
    if np.any(valid):
        weighted_rms = math.sqrt(float(
            np.sum(truth[valid] * relative_error[valid] ** 2)
            / np.sum(truth[valid])))
        p95 = float(np.percentile(relative_error[valid], 95))
    else:
        weighted_rms = math.inf
        p95 = math.inf

    energy = metadata["energy_keV"]
    minimum = energy / (1.0 + 2.0 * energy / 510.99895)
    edges = np.linspace(0.0, energy, metric_bins + 1)
    outside = edges[1:] <= minimum
    below_free_edge_yield = float(model[outside].sum())
    return {
        "energy_keV": energy,
        "Z": int(metadata["Z"]),
        "tau": metadata["tau"],
        "areal_density_g_cm2": metadata["areal_density_g_cm2"],
        "events": int(metadata["events"]),
        "truth_yield": truth_yield,
        "truth_yield_sigma": truth_yield_sigma,
        "model_yield": model_yield,
        "yield_relative_error": yield_relative_error,
        "yield_pull": (
            (model_yield - truth_yield) / truth_yield_sigma
            if truth_yield_sigma > 0.0 else math.inf),
        "shape_total_variation": total_variation,
        "weighted_rms_relative_error": weighted_rms,
        "p95_relative_error": p95,
        "below_free_edge_yield": below_free_edge_yield,
        "valid_bins": int(valid.sum()),
    }


def write_plots(
    output: Path,
    records: list[dict],
    loaded: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
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
                 and record["Z"] == z_value),
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
        axes[0, column].set_ylabel("Single-scatter yield error (%)")
        axes[1, column].set_ylabel("Normalized-shape TV distance (%)")
        axes[1, column].set_xlabel("Process optical depth")
        for axis in axes[:, column]:
            axis.set_xscale("log")
            axis.grid(alpha=0.2)
    axes[0, -1].legend(fontsize=8)
    figure.suptitle("Physics single-Compton validation across E, Z, and shield")
    figure.savefig(output / "single_scatter_grid_summary.png", dpi=180)
    plt.close(figure)

    ranked = sorted(
        records,
        key=lambda record: max(
            abs(record["yield_relative_error"]),
            record["shape_total_variation"]),
        reverse=True)
    with PdfPages(output / "single_scatter_worst_cells.pdf") as pdf:
        for record in ranked[:worst]:
            identity = record["cell_id"]
            truth, variance, model = loaded[identity]
            energy = record["energy_keV"]
            edges = np.linspace(0.0, energy, truth.size + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            sigma = np.sqrt(np.maximum(variance, 0.0))
            figure, (shape_axis, ratio_axis) = plt.subplots(
                2, 1, figsize=(11.0, 7.0), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]},
                constrained_layout=True)
            shape_axis.stairs(
                truth, edges, color="#202020", linewidth=1.3,
                label="CeeLo exactly-one-Compton")
            shape_axis.stairs(
                model, edges, color="#d1495b", linewidth=1.1,
                label="physics model")
            shape_axis.fill_between(
                centers, np.maximum(0.0, truth - sigma), truth + sigma,
                step="mid", color="#4c78a8", alpha=0.20,
                label="MC uncertainty")
            shape_axis.set_ylabel("Yield / emitted photon / bin")
            shape_axis.legend(fontsize=8)
            shape_axis.grid(alpha=0.2)

            valid = truth > 0.0
            ratio = np.full_like(truth, np.nan)
            ratio[valid] = model[valid] / truth[valid]
            ratio_axis.plot(centers, ratio, color="#d1495b", linewidth=0.9)
            ratio_axis.axhline(1.0, color="#202020", linewidth=0.8)
            ratio_axis.set_ylim(0.0, 2.0)
            ratio_axis.set_ylabel("model / truth")
            ratio_axis.set_xlabel("Escaped photon energy (keV)")
            ratio_axis.grid(alpha=0.2)
            minimum = energy / (1.0 + 2.0 * energy / 510.99895)
            for axis in (shape_axis, ratio_axis):
                axis.axvline(
                    minimum, color="#777777", linestyle="--",
                    linewidth=0.8)
            figure.suptitle(
                f"{energy:g} keV, Z={record['Z']}, "
                f"tau={record['tau']:g}; yield error "
                f"{100.0 * record['yield_relative_error']:+.2f}%, "
                f"shape TV {100.0 * record['shape_total_variation']:.2f}%")
            pdf.savefig(figure)
            plt.close(figure)


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
    parser.add_argument(
        "--cell", action="append", default=[],
        metavar="ENERGY,Z,TAU",
        help="explicit validation cell; repeat to replace the Cartesian grid")
    parser.add_argument("--bins", type=int, default=256)
    parser.add_argument(
        "--metric-bins", type=int, default=64,
        help="coarser binning used for stable shape-error metrics")
    parser.add_argument("--events", type=int, default=1000000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--worst-plots", type=int, default=12)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--force-model", action="store_true",
        help="regenerate only functional-model cells, retaining MC truth")
    parser.add_argument(
        "--source-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_escape_campaign"))
    parser.add_argument(
        "--model-driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_model_evaluate"))
    args = parser.parse_args()
    if args.cell:
        cells = []
        for text in args.cell:
            fields = text.split(",")
            if len(fields) != 3:
                raise ValueError(
                    f"--cell must be ENERGY,Z,TAU, received: {text}")
            cells.append(
                (float(fields[0]), int(fields[1]), float(fields[2])))
    else:
        cells = [
            (energy, z_value, tau)
            for energy in args.energies
            for z_value in args.z_values
            for tau in args.tau_values
        ]

    if args.action in {"run", "all"}:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [
                pool.submit(
                    run_cell, args.output, *cell, args.bins, args.events,
                    args.source_driver, args.model_driver, args.force,
                    args.force_model)
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
        truth, variance, metadata = load_single(
            args.output / "cells" / f"{identity}_ceelo.csv", args.bins)
        model, _, _ = load_single(
            args.output / "cells" / f"{identity}_model.csv", args.bins)
        record = calculate_metrics(
            truth, variance, model, metadata, args.metric_bins)
        record["cell_id"] = identity
        records.append(record)
        loaded[identity] = (truth, variance, model)

    args.output.mkdir(parents=True, exist_ok=True)
    fieldnames = list(records[0])
    with (args.output / "single_scatter_metrics.csv").open(
            "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    finite = [
        record for record in records
        if math.isfinite(record["weighted_rms_relative_error"])]
    def finite_summary(key: str, percentile: float) -> float | None:
        if not finite:
            return None
        values = [record[key] for record in finite]
        return float(
            np.median(values)
            if percentile == 50.0
            else np.percentile(values, percentile))

    report = {
        "schema": 1,
        "cells": len(records),
        "energies_keV": sorted({cell[0] for cell in cells}),
        "Z_values": sorted({cell[1] for cell in cells}),
        "tau_values": sorted({cell[2] for cell in cells}),
        "events_per_cell": args.events,
        "bins": args.bins,
        "metric_bins": args.metric_bins,
        "source_rayleigh": True,
        "source_compton_doppler": True,
        "optical_depth_basis": "mu_pe + mu_cs + mu_rs + mu_pp",
        "strict_truth_label": (
            "exactly one Compton, with any number of Rayleigh scatters"),
        "yield_error_median_abs": float(np.median([
            abs(record["yield_relative_error"]) for record in records])),
        "yield_error_p95_abs": float(np.percentile([
            abs(record["yield_relative_error"]) for record in records], 95)),
        "shape_tv_median": float(np.median([
            record["shape_total_variation"] for record in records])),
        "shape_tv_p95": float(np.percentile([
            record["shape_total_variation"] for record in records], 95)),
        "weighted_rms_median": finite_summary(
            "weighted_rms_relative_error", 50.0),
        "weighted_rms_p95": finite_summary(
            "weighted_rms_relative_error", 95.0),
        "bin_p95_error_median": finite_summary(
            "p95_relative_error", 50.0),
        "bin_p95_error_p95": finite_summary(
            "p95_relative_error", 95.0),
        "maximum_absolute_yield_pull": max(
            abs(record["yield_pull"]) for record in records),
        "maximum_below_free_edge_yield": max(
            record["below_free_edge_yield"] for record in records),
    }
    (args.output / "single_scatter_report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    write_plots(args.output, records, loaded, args.worst_plots)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
