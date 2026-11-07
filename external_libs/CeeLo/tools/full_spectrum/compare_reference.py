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

"""Rebin, resolve and compare CeeLo, Geant4 or InterSpec/GADRAS spectra.

Accepted inputs are campaign-cell CSVs, the weighted Geant4 histogram CSV
emitted by ceelo_g4val, or a simple two-column energy/count CSV export.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np


def noncomment_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines()
            if line and not line.startswith("#")
            and not line.startswith("RealTime:")]


def campaign_spectrum(rows: list[dict[str, str]],
                      selected_component: str | None = None
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    events = int(rows[0]["events"])
    kind = rows[0]["cell_kind"]
    excluded = {"absorbed"} if kind == "source" else {"zero_deposit"}
    boundaries = sorted({
        float(value)
        for row in rows
        for value in (row["bin_lo_keV"], row["bin_hi_keV"])
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])
    })
    edges = np.asarray(boundaries)
    values = np.zeros(edges.size - 1)
    variance = np.zeros(edges.size - 1)
    for row in rows:
        if row["component"] in excluded:
            continue
        if selected_component is not None \
                and row["component"] != selected_component:
            continue
        lo, hi = float(row["bin_lo_keV"]), float(row["bin_hi_keV"])
        if hi > lo:
            index = int(np.searchsorted(edges, lo))
        else:
            index = min(
                values.size - 1,
                int(np.searchsorted(edges, lo, side="right") - 1))
        if 0 <= index < values.size:
            values[index] += float(row["sum_w"]) / events
            variance[index] += float(row["sum_w2"]) / (events * events)
    return edges, values, variance


def weighted_g4_spectrum(lines: list[str], step_hint: float | None = None,
                         selected_component: str | None = None
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = list(csv.DictReader(lines))
    component = selected_component or "total"
    rows = [row for row in rows if row["component"] == component]
    lower = np.asarray([float(row["energy_lower_keV"]) for row in rows])
    if lower.size < 2 and step_hint is None:
        raise ValueError("Geant4 weighted histogram has too few populated bins")
    step = step_hint or float(np.min(np.diff(np.unique(lower))))
    max_edge = float(lower.max() + step)
    edges = np.arange(0.0, max_edge + 0.5 * step, step)
    values = np.zeros(edges.size - 1)
    variance = np.zeros(edges.size - 1)
    # Geant4 output is weighted counts. Normalization may be performed by the
    # caller with --scale when comparing efficiencies.
    for row in rows:
        index = int(round(float(row["energy_lower_keV"]) / step))
        values[index] = float(row["sum_w"])
        variance[index] = float(row["sum_w2"])
    return edges, values, variance


def simple_spectrum(lines: list[str]
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = list(csv.DictReader(lines))
    names = rows[0].keys()
    energy_name = next(
        name for name in names if "energy" in name.lower()
        or "lower" in name.lower())
    count_name = next(
        name for name in names if name != energy_name
        and ("count" in name.lower() or "yield" in name.lower()
             or "value" in name.lower()))
    energy = np.asarray([float(row[energy_name]) for row in rows])
    values = np.asarray([float(row[count_name]) for row in rows])
    order = np.argsort(energy)
    energy, values = energy[order], values[order]
    if energy.size < 2:
        raise ValueError("Simple spectrum needs at least two energy rows")
    widths = np.diff(energy)
    final_width = float(np.median(widths))
    edges = np.concatenate((energy, [energy[-1] + final_width]))
    return edges, values, np.zeros_like(values)


def load(path: Path, selected_component: str | None = None
         ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = noncomment_lines(path)
    if not lines:
        raise ValueError(f"empty spectrum: {path}")
    header = lines[0]
    if "cell_kind" in header:
        return campaign_spectrum(
            list(csv.DictReader(lines)), selected_component)
    if "component" in header and "sum_w2" in header:
        step_hint = None
        for raw in path.read_text().splitlines():
            match = re.match(
                r"#\s*Bin width:\s*([0-9.eE+-]+)\s*keV", raw)
            if match:
                step_hint = float(match.group(1))
                break
        return weighted_g4_spectrum(
            lines, step_hint, selected_component)
    if selected_component is not None:
        raise ValueError(
            "--component is only available for component-labelled inputs")
    return simple_spectrum(lines)


def load_components(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Load named components on the file's native energy grid."""
    lines = noncomment_lines(path)
    if not lines:
        raise ValueError(f"empty spectrum: {path}")
    header = lines[0]
    edges, total, _ = load(path)
    components: dict[str, np.ndarray] = {}
    if "cell_kind" in header:
        rows = list(csv.DictReader(lines))
        events = int(rows[0]["events"])
        kind = rows[0]["cell_kind"]
        excluded = {"absorbed"} if kind == "source" else {"zero_deposit"}
        for row in rows:
            name = row["component"]
            if name in excluded:
                continue
            values = components.setdefault(name, np.zeros(edges.size - 1))
            lo, hi = float(row["bin_lo_keV"]), float(row["bin_hi_keV"])
            if hi > lo:
                index = int(np.searchsorted(edges, lo))
            else:
                index = int(np.searchsorted(edges, lo, side="right") - 1)
            index = min(values.size - 1, max(0, index))
            values[index] += float(row["sum_w"]) / events
    elif "component" in header and "sum_w2" in header:
        rows = list(csv.DictReader(lines))
        total_rows = [row for row in rows if row["component"] == "total"]
        if total_rows and edges.size >= 2:
            step = float(edges[1] - edges[0])
            for row in rows:
                name = row["component"]
                if name == "total":
                    continue
                values = components.setdefault(
                    name, np.zeros(edges.size - 1))
                index = int(round(
                    float(row["energy_lower_keV"]) / step))
                if 0 <= index < values.size:
                    values[index] += float(row["sum_w"])
    if not components:
        components["total"] = total
    return edges, components


def spectrum_metadata(path: Path) -> tuple[str | None, float | None]:
    """Return transport kind and incident energy when the format carries it."""
    raw_lines = path.read_text().splitlines()
    lines = noncomment_lines(path)
    if lines and "cell_kind" in lines[0]:
        first = next(csv.DictReader(lines))
        return first["cell_kind"], float(first["energy_keV"])
    for line in raw_lines:
        match = re.match(
            r"#\s*Primary energy:\s*([0-9.eE+-]+)\s*keV", line)
        if match:
            return "detector", float(match.group(1))
    return None, None


def rebin(old_edges: np.ndarray, old_values: np.ndarray,
          new_edges: np.ndarray) -> np.ndarray:
    cumulative = np.concatenate(([0.0], np.cumsum(old_values)))
    # Piecewise-linear cumulative integral of a piecewise-constant density.
    sampled = np.interp(new_edges, old_edges, cumulative)
    return np.diff(sampled)


def rebin_variance(
    old_edges: np.ndarray, old_variance: np.ndarray, new_edges: np.ndarray
) -> np.ndarray:
    """Propagate independent bin variances through overlap-weighted rebinning."""
    result = np.zeros(new_edges.size - 1)
    old_index = 0
    new_index = 0
    while old_index < old_variance.size and new_index < result.size:
        lo = max(old_edges[old_index], new_edges[new_index])
        hi = min(old_edges[old_index + 1], new_edges[new_index + 1])
        if hi > lo:
            fraction = (hi - lo) / (
                old_edges[old_index + 1] - old_edges[old_index])
            result[new_index] += old_variance[old_index] * fraction * fraction
        if old_edges[old_index + 1] <= new_edges[new_index + 1]:
            old_index += 1
        else:
            new_index += 1
    return result


def resolve(edges: np.ndarray, values: np.ndarray,
            a0: float, a1: float, a2: float) -> np.ndarray:
    if a0 == 0.0 and a1 == 0.0 and a2 == 0.0:
        return values.copy()
    result = np.zeros_like(values)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = float(np.median(np.diff(edges)))
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    for index, (energy, weight) in enumerate(zip(centers, values)):
        if weight == 0.0:
            continue
        fwhm2 = a0 + a1 * energy + a2 * energy * energy
        if fwhm2 <= 0.0:
            raise ValueError("resolution FWHM squared is nonpositive")
        sigma = math.sqrt(fwhm2 / (8.0 * math.log(2.0)))
        radius = max(1, int(math.ceil(6.0 * sigma / width)))
        first, last = max(0, index - radius), min(values.size, index + radius + 1)
        window_zlo = (edges[first] - energy) * inv_sqrt2 / sigma
        window_zhi = (edges[last] - energy) * inv_sqrt2 / sigma
        normalization = max(
            0.5 * (math.erf(window_zhi) - math.erf(window_zlo)), 1e-300)
        for target in range(first, last):
            zlo = (edges[target] - energy) * inv_sqrt2 / sigma
            zhi = (edges[target + 1] - energy) * inv_sqrt2 / sigma
            result[target] += weight * 0.5 * (
                math.erf(zhi) - math.erf(zlo)) / normalization
    return result


def resolve_variance(
    edges: np.ndarray, variance: np.ndarray,
    a0: float, a1: float, a2: float
) -> np.ndarray:
    """Propagate independent input-bin variances through Gaussian smearing."""
    if a0 == 0.0 and a1 == 0.0 and a2 == 0.0:
        return variance.copy()
    result = np.zeros_like(variance)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = float(np.median(np.diff(edges)))
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    for index, (energy, input_variance) in enumerate(zip(centers, variance)):
        if input_variance == 0.0:
            continue
        fwhm2 = a0 + a1 * energy + a2 * energy * energy
        if fwhm2 <= 0.0:
            raise ValueError("resolution FWHM squared is nonpositive")
        sigma = math.sqrt(fwhm2 / (8.0 * math.log(2.0)))
        radius = max(1, int(math.ceil(6.0 * sigma / width)))
        first = max(0, index - radius)
        last = min(variance.size, index + radius + 1)
        window_zlo = (edges[first] - energy) * inv_sqrt2 / sigma
        window_zhi = (edges[last] - energy) * inv_sqrt2 / sigma
        normalization = max(
            0.5 * (math.erf(window_zhi) - math.erf(window_zlo)), 1e-300)
        for target in range(first, last):
            zlo = (edges[target] - energy) * inv_sqrt2 / sigma
            zhi = (edges[target + 1] - energy) * inv_sqrt2 / sigma
            probability = 0.5 * (
                math.erf(zhi) - math.erf(zlo)) / normalization
            result[target] += input_variance * probability * probability
    return result


def prepare_components(
    path: Path,
    target_edges: np.ndarray,
    scale: float,
    resolution: tuple[float, float, float],
) -> dict[str, np.ndarray]:
    edges, components = load_components(path)
    return {
        name: resolve(
            target_edges,
            rebin(edges, values, target_edges) * scale,
            *resolution,
        )
        for name, values in components.items()
    }


def physics_landmarks(
    incident_keV: float | None, spectrum_kind: str | None
) -> list[tuple[str, float]]:
    if incident_keV is None or incident_keV <= 0.0:
        return []
    electron_mass_keV = 510.99895
    backscatter = incident_keV / (
        1.0 + 2.0 * incident_keV / electron_mass_keV)
    if spectrum_kind == "source":
        landmarks = [
            ("single-scatter minimum", backscatter),
            ("unscattered / elastic", incident_keV),
        ]
    else:
        landmarks = [
            ("backscatter", backscatter),
            ("Compton edge", incident_keV - backscatter),
            ("FEP / E0", incident_keV),
        ]
    if incident_keV > 1021.9979:
        landmarks.append(("511", electron_mass_keV))
        if spectrum_kind != "source":
            landmarks.append(
                ("single escape", incident_keV - electron_mass_keV))
    if incident_keV > 1532.99685 and spectrum_kind != "source":
        landmarks.append(
            ("double escape", incident_keV - 2.0 * electron_mass_keV))
    return sorted(landmarks, key=lambda item: item[1])


def write_comparison_plot(
    output: Path,
    edges: np.ndarray,
    truth: np.ndarray,
    candidate: np.ndarray,
    variance: np.ndarray,
    truth_components: dict[str, np.ndarray],
    candidate_components: dict[str, np.ndarray],
    truth_label: str,
    candidate_label: str,
    title: str,
    incident_keV: float | None,
    spectrum_kind: str | None,
    linear_y: bool,
    show_components: bool,
    dpi: int,
) -> None:
    # Plotting remains an optional research dependency.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centers = 0.5 * (edges[:-1] + edges[1:])
    sigma = np.sqrt(np.maximum(variance, 0.0))
    positive = np.concatenate(
        (truth[truth > 0.0], candidate[candidate > 0.0]))
    floor = max(
        float(np.max(positive)) * 1e-9 if positive.size else 1e-30,
        1e-30,
    )
    valid = truth > floor
    ratio = np.full_like(truth, np.nan)
    relative = np.full_like(truth, np.nan)
    ratio_uncertainty = np.full_like(truth, np.nan)
    ratio[valid] = candidate[valid] / truth[valid]
    relative[valid] = 100.0 * (
        candidate[valid] - truth[valid]) / truth[valid]
    ratio_uncertainty[valid] = sigma[valid] / truth[valid]

    rows = 4 if show_components else 3
    height_ratios = [2.4, 1.8, 1.0, 1.0] if show_components \
        else [2.5, 1.0, 1.0]
    figure, axes = plt.subplots(
        rows, 1, figsize=(13.0, 10.5 if show_components else 8.5),
        sharex=True, gridspec_kw={"height_ratios": height_ratios},
        constrained_layout=True,
    )
    total_ax = axes[0]
    component_ax = axes[1] if show_components else None
    ratio_ax = axes[-2]
    residual_ax = axes[-1]

    total_ax.stairs(
        np.maximum(truth, floor), edges, label=truth_label,
        color="#202020", linewidth=1.5)
    total_ax.stairs(
        np.maximum(candidate, floor), edges, label=candidate_label,
        color="#d1495b", linewidth=1.25)
    total_ax.fill_between(
        centers,
        np.maximum(truth - sigma, floor),
        np.maximum(truth + sigma, floor),
        color="#4c78a8", alpha=0.22, step="mid",
        label=f"{truth_label} MC uncertainty",
    )
    if not linear_y:
        total_ax.set_yscale("log")
    total_ax.set_ylabel("Expected counts / bin")
    total_ax.legend(loc="best", ncols=3, fontsize=8)
    total_ax.grid(alpha=0.2)

    if component_ax is not None:
        color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        component_names = sorted(
            set(truth_components) | set(candidate_components))
        for index, name in enumerate(component_names):
            color = color_cycle[index % len(color_cycle)]
            if name in truth_components and np.any(truth_components[name] > 0):
                component_ax.stairs(
                    np.maximum(truth_components[name], floor), edges,
                    color=color, linewidth=1.0, label=f"{name} truth")
            if (name in candidate_components
                    and np.any(candidate_components[name] > 0)):
                component_ax.stairs(
                    np.maximum(candidate_components[name], floor), edges,
                    color=color, linewidth=1.0, linestyle=":",
                    label=f"{name} candidate")
        if not linear_y:
            component_ax.set_yscale("log")
        component_ax.set_ylabel("Components")
        component_ax.legend(
            loc="upper right", ncols=3, fontsize=7,
            framealpha=0.85)
        component_ax.grid(alpha=0.2)

    ratio_ax.plot(centers, ratio, color="#d1495b", linewidth=0.8)
    ratio_ax.fill_between(
        centers, 1.0 - ratio_uncertainty, 1.0 + ratio_uncertainty,
        color="#4c78a8", alpha=0.22, step="mid")
    ratio_ax.axhline(1.0, color="#202020", linewidth=0.8)
    finite_ratio = ratio[np.isfinite(ratio)]
    if finite_ratio.size:
        limit = min(4.0, max(
            0.10, float(np.percentile(np.abs(finite_ratio - 1.0), 99.0))
            * 1.15))
        ratio_ax.set_ylim(max(0.0, 1.0 - limit), 1.0 + limit)
    ratio_ax.set_ylabel("candidate / truth")
    ratio_ax.grid(alpha=0.2)

    residual_ax.plot(
        centers, relative, color="#d1495b", linewidth=0.8)
    residual_ax.fill_between(
        centers, -100.0 * ratio_uncertainty,
        100.0 * ratio_uncertainty,
        color="#4c78a8", alpha=0.22, step="mid")
    residual_ax.axhline(0.0, color="#202020", linewidth=0.8)
    residual_ax.axhspan(-3.0, 3.0, color="#59a14f", alpha=0.10)
    finite_relative = np.abs(relative[np.isfinite(relative)])
    if finite_relative.size:
        limit = min(500.0, max(10.0, float(
            np.percentile(finite_relative, 98.0)) * 1.15))
        residual_ax.set_ylim(-limit, limit)
    residual_ax.set_ylabel("relative residual (%)")
    residual_ax.set_xlabel("Energy (keV)")
    residual_ax.grid(alpha=0.2)
    if incident_keV is not None and incident_keV > 0.0:
        residual_ax.set_xlim(
            max(0.0, edges[0]), min(edges[-1], 1.05 * incident_keV))

    landmarks = physics_landmarks(incident_keV, spectrum_kind)
    for axis in axes:
        for _, energy in landmarks:
            if edges[0] <= energy <= edges[-1]:
                axis.axvline(
                    energy, color="#7f7f7f", linewidth=0.7,
                    alpha=0.45, linestyle="--")
    for name, energy in landmarks:
        if edges[0] <= energy <= edges[-1]:
            total_ax.annotate(
                name, xy=(energy, 1.0), xycoords=("data", "axes fraction"),
                xytext=(2, -2), textcoords="offset points",
                rotation=90, va="top", ha="left", fontsize=7,
                color="#555555")

    figure.suptitle(title)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("truth", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--channels", type=int, default=4096)
    parser.add_argument("--max-kev", type=float, default=12000.0)
    parser.add_argument("--resolution", type=float, nargs=3,
                        metavar=("A0", "A1", "A2"),
                        default=(0.0, 0.0, 0.0))
    parser.add_argument("--truth-scale", type=float, default=1.0)
    parser.add_argument("--candidate-scale", type=float, default=1.0)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--plot", type=Path,
        help="write an overlay/component/ratio/residual comparison (PNG or PDF)")
    parser.add_argument("--title")
    parser.add_argument("--truth-label")
    parser.add_argument("--candidate-label")
    parser.add_argument(
        "--component",
        help="compare only this labelled component (for example, "
             "single_compton)")
    parser.add_argument(
        "--incident-kev", type=float,
        help="incident energy for Compton/backscatter/escape landmarks")
    parser.add_argument("--linear-y", action="store_true")
    parser.add_argument("--no-components", action="store_true")
    parser.add_argument("--plot-dpi", type=int, default=180)
    args = parser.parse_args()

    target_edges = np.linspace(0.0, args.max_kev, args.channels + 1)
    truth_edges, truth_values, truth_variance = load(
        args.truth, args.component)
    candidate_edges, candidate_values, _ = load(
        args.candidate, args.component)
    truth = rebin(truth_edges, truth_values, target_edges) * args.truth_scale
    candidate = rebin(
        candidate_edges, candidate_values, target_edges) * args.candidate_scale
    variance = rebin_variance(
        truth_edges, truth_variance, target_edges) * args.truth_scale ** 2
    truth = resolve(target_edges, truth, *args.resolution)
    candidate = resolve(target_edges, candidate, *args.resolution)
    variance = resolve_variance(
        target_edges, variance, *args.resolution)

    rel_unc = np.sqrt(np.maximum(variance, 0.0)) / np.maximum(truth, 1e-30)
    valid = (truth > max(1e-12, truth.max() * 1e-8)) \
        & ((variance == 0.0) | (rel_unc <= 0.05))
    rel = np.abs(candidate - truth) / np.maximum(truth, 1e-30)
    weights = np.where(valid, truth, 0.0)
    report = {
        "truth_integral": float(truth.sum()),
        "candidate_integral": float(candidate.sum()),
        "integral_relative_error": float(
            abs(candidate.sum() - truth.sum()) / max(truth.sum(), 1e-30)),
        "weighted_rms_relative_error": float(math.sqrt(
            np.sum(weights * rel * rel) / max(np.sum(weights), 1e-30))),
        "p95_relative_error": float(np.percentile(rel[valid], 95))
        if np.any(valid) else None,
        "valid_bins": int(valid.sum()),
        "channels": args.channels,
        "resolution": list(args.resolution),
    }
    if args.plot:
        truth_kind, truth_incident = spectrum_metadata(args.truth)
        candidate_kind, candidate_incident = spectrum_metadata(args.candidate)
        incident = args.incident_kev
        if incident is None:
            incident = truth_incident or candidate_incident
        truth_components = prepare_components(
            args.truth, target_edges, args.truth_scale,
            tuple(args.resolution))
        candidate_components = prepare_components(
            args.candidate, target_edges, args.candidate_scale,
            tuple(args.resolution))
        truth_label = args.truth_label or args.truth.stem
        candidate_label = args.candidate_label or args.candidate.stem
        inferred_kind = truth_kind or candidate_kind or "spectrum"
        title = args.title or (
            f"{inferred_kind.capitalize()} response comparison"
            + (f" at {incident:g} keV" if incident is not None else ""))
        write_comparison_plot(
            args.plot, target_edges, truth, candidate, variance,
            truth_components, candidate_components,
            truth_label, candidate_label, title, incident,
            inferred_kind,
            args.linear_y, not args.no_components, args.plot_dpi)
        report["plot"] = str(args.plot)
        report["incident_keV"] = incident
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
