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

"""Fit positive source-secondary photon measures.

The fitted components deliberately exclude fluorescence:

* the unscattered 511-keV annihilation atom;
* the Compton-scattered at-rest annihilation continuum;
* the positron in-flight annihilation continuum;
* bremsstrahlung rooted in photoelectric, Compton, and pair electrons.

Each continuum uses a positive cubic B-spline in the same inverse-energy warp
as the high-order Compton model. Smooth Gaussian RBF fields map
(log E0, Z, log total optical depth) to component yields and normalized shape
coefficients. The serialized payload contains coefficients only, never source
histograms.
"""

from __future__ import annotations

import argparse
import csv
import bank_io
import hashlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import nnls


ELECTRON_MASS_KEV = 510.99895
BASIS_SIZE = 32
WARP_SCALE = 2.0

# The 7-D process coordinate; must match
# ceelo::full_spectrum::kSourceRbfMaxFeatures.
PROCESS_FEATURE_SIZE = 7
# 90% CL Poisson upper bound on the mean given zero observed counts.
POISSON_ZERO_UPPER_BOUND = 2.3
CONTINUOUS_COMPONENTS = (
    "pair_annihilation_continuum",
    "pair_annihilation_in_flight",
    "photoelectron_brems",
    "compton_electron_brems",
    "pair_electron_brems",
)
# Line components: one yield per cell, a flat shape. The fluorescence atoms are
# the shield's K-alpha, K-beta and L x-ray families (the C++ evaluator spreads
# each family's yield over its lines).
FLUORESCENCE_COMPONENTS = (
    "fluorescence_k_alpha_atom",
    "fluorescence_k_beta_atom",
    "fluorescence_l_atom",
)
ATOM_COMPONENTS = ("pair_511_atom",) + FLUORESCENCE_COMPONENTS
# Bank row label of each atom component.
ATOM_BANK_LABEL = {
    "pair_511_atom": "pair_511",
    "fluorescence_k_alpha_atom": "fluorescence_k_alpha",
    "fluorescence_k_beta_atom": "fluorescence_k_beta",
    "fluorescence_l_atom": "fluorescence_l",
}
# Order of SourceSecondarySurrogatePayload::components() in the C++ model.
ALL_COMPONENTS = ("pair_511_atom",) + CONTINUOUS_COMPONENTS + FLUORESCENCE_COMPONENTS
# Populated by fit_nodes; reported so the anchor count is never silent.
ANCHOR_COUNTS: dict[str, int] = {}


@dataclass
class Cell:
    identity: str
    path: Path
    energy_keV: float
    z_value: int
    tau: float
    areal_density: float
    events: float
    edges: np.ndarray
    values: dict[str, np.ndarray]
    variance: dict[str, np.ndarray]
    atom_yield: dict[str, float]
    atom_variance: dict[str, float]
    feature: np.ndarray
    split: str
    hydrogen_mass_fraction: float = 0.0
    # Both built by the C++ helper, never recomputed here -- see
    # attach_process_features and attach_component_bases.
    process_features: dict[str, np.ndarray] | None = None
    component_bases: dict[str, float] | None = None


@dataclass
class Node:
    feature: np.ndarray
    log_conditional_yield: float
    shape_coefficients: np.ndarray


def feature(energy_keV: float, z_value: int, tau: float) -> np.ndarray:
    return np.asarray([
        math.log(energy_keV / 10.0) / math.log(1200.0),
        z_value / 92.0,
        math.log1p(tau) / math.log(25.0),
    ])


def split_for(identity: str, z_value: int, holdout_z: set[int]) -> str:
    if z_value in holdout_z:
        return "test"
    return "train"


def load_cell(path: Path, holdout_z: set[int]) -> Cell:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty source cell: {path}")
    first = rows[0]
    for switch in ("fluorescence_enabled", "impact_ionization_enabled"):
        if first.get(switch, "0").strip().lower() not in {"1", "true", "yes"}:
            raise ValueError(f"{switch} must be on (the bank must hold the shield's "
                             f"x-rays for the fluorescence atoms): {path}")
    energy = float(first["energy_keV"])
    events = float(first["events"])
    continuous_rows = [
        row for row in rows
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])
    ]
    boundaries = sorted({
        float(value)
        for row in continuous_rows
        for value in (row["bin_lo_keV"], row["bin_hi_keV"])
    })
    if len(boundaries) < 2:
        raise ValueError(f"no continuous grid in {path}")
    edges = np.asarray(boundaries)
    values = {
        name: np.zeros(edges.size - 1)
        for name in CONTINUOUS_COMPONENTS}
    variance = {
        name: np.zeros(edges.size - 1)
        for name in CONTINUOUS_COMPONENTS}
    atom_of_label = {label: name for name, label in ATOM_BANK_LABEL.items()}
    atom_yield = {name: 0.0 for name in ATOM_COMPONENTS}
    atom_variance = {name: 0.0 for name in ATOM_COMPONENTS}
    for row in rows:
        name = row["component"]
        lo = float(row["bin_lo_keV"])
        hi = float(row["bin_hi_keV"])
        if name in atom_of_label and not hi > lo:
            atom = atom_of_label[name]
            atom_yield[atom] += float(row["sum_w"]) / events
            atom_variance[atom] += float(row["sum_w2"]) / (events * events)
        elif name in values and hi > lo:
            index = int(np.searchsorted(edges, lo, side="right") - 1)
            if 0 <= index < values[name].size:
                values[name][index] += float(row["sum_w"]) / events
                variance[name][index] += (
                    float(row["sum_w2"]) / (events * events))
    identity = bank_io.cell_id(path)
    z_value = int(first["Z"])
    tau = float(first["tau"])
    hydrogen = float(first["hydrogen_mass_fraction"])
    return Cell(
        identity=identity,
        path=path,
        energy_keV=energy,
        z_value=z_value,
        tau=tau,
        areal_density=float(first["areal_density_g_cm2"]),
        events=events,
        edges=edges,
        values=values,
        variance=variance,
        atom_yield=atom_yield,
        atom_variance=atom_variance,
        feature=feature(energy, z_value, tau),
        split=split_for(identity, z_value, holdout_z),
        hydrogen_mass_fraction=hydrogen,
    )


def attach_component_bases(cells: list[Cell], helper: Path) -> None:
    """Read the finished mechanistic baselines from the C++ helper.

    Deliberately NOT recomputed here, for the same reason as the process
    features one layer up: the fit stores log(truth / base) and the evaluator
    multiplies the base back in, so any C++/Python disagreement is a purely
    multiplicative bias that no metric in report.json can see.
    """
    lines = [
        f"{cell.energy_keV} {cell.z_value} {cell.areal_density} "
        f"{cell.hydrogen_mass_fraction}\n"
        for cell in cells]
    completed = subprocess.run(
        [str(helper), "--emit-component-bases"],
        input="".join(lines), text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(ALL_COMPONENTS) * len(cells):
        raise RuntimeError(
            "component-base helper returned wrong row count: "
            f"{len(rows)} for {len(cells)} cells")
    for index, cell in enumerate(cells):
        block = rows[index * len(ALL_COMPONENTS):
                     (index + 1) * len(ALL_COMPONENTS)]
        bases = {row["component"]: float(row["base"]) for row in block}
        if set(bases) != set(ALL_COMPONENTS):
            raise RuntimeError(
                f"component-base helper missed a component for {cell.identity}")
        cell.component_bases = bases


def attach_process_features(cells: list[Cell], helper: Path) -> None:
    """Read the finished 7-D process feature vectors from the C++ helper.

    Deliberately NOT recomputed here. The helper calls
    detail::source_process_features, the same function the evaluator uses, so
    fitter/evaluator lockstep on this coordinate is structural rather than a
    convention two files have to keep agreeing on by hand.
    """
    lines = [
        f"{cell.energy_keV} {cell.z_value} {cell.areal_density} "
        f"{cell.hydrogen_mass_fraction}\n"
        for cell in cells]
    completed = subprocess.run(
        [str(helper), "--emit-process-features"],
        input="".join(lines), text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(ALL_COMPONENTS) * len(cells):
        raise RuntimeError(
            "process-feature helper returned wrong row count: "
            f"{len(rows)} for {len(cells)} cells")
    axes = [f"f{index}" for index in range(PROCESS_FEATURE_SIZE)]
    for index, cell in enumerate(cells):
        block = rows[index * len(ALL_COMPONENTS):
                     (index + 1) * len(ALL_COMPONENTS)]
        features = {}
        for row in block:
            features[row["component"]] = np.asarray(
                [float(row[axis]) for axis in axes])
        if set(features) != set(ALL_COMPONENTS):
            raise RuntimeError(
                f"process-feature helper missed a component for {cell.identity}")
        cell.process_features = features


def component_feature(cell: Cell, name: str) -> np.ndarray:
    assert cell.process_features is not None
    return cell.process_features[name]


def component_base(cell: Cell, name: str) -> float:
    """The mechanistic baseline the fitted conditional yield rides on.

    A lookup, not a calculation.  detail::source_component_bases in
    FullSpectrumModel.h is the single implementation; the helper hands the
    finished numbers over in attach_component_bases.  Mirror of
    component_feature one layer down.
    """
    assert cell.component_bases is not None
    return cell.component_bases[name]


_BSPLINE_CACHE: dict[tuple, np.ndarray] = {}


def warped_bspline_matrix_cached(
    cell_identity: str,
    output_edges: np.ndarray,
    incident_energy_keV: float,
) -> np.ndarray:
    key = (cell_identity, round(incident_energy_keV, 9))
    cached = _BSPLINE_CACHE.get(key)
    if cached is None:
        cached = warped_bspline_matrix(output_edges, incident_energy_keV)
        _BSPLINE_CACHE[key] = cached
    return cached


def warped_bspline_matrix(
    output_edges: np.ndarray,
    incident_energy_keV: float,
) -> np.ndarray:
    degree = 3
    internal_count = BASIS_SIZE - degree - 1
    knots = np.concatenate((
        np.zeros(degree + 1),
        np.linspace(0.0, 1.0, internal_count + 2)[1:-1],
        np.ones(degree + 1),
    ))

    clipped = np.clip(output_edges, 0.0, incident_energy_keV)

    def warped(energy: np.ndarray) -> np.ndarray:
        safe = np.maximum(energy, 1e-100)
        loss = 0.5 * ELECTRON_MASS_KEV * (
            1.0 / safe - 1.0 / incident_energy_keV)
        return np.clip(loss / (WARP_SCALE + loss), 0.0, 1.0)

    upper = warped(clipped[:-1])
    lower = warped(clipped[1:])
    upper[output_edges[:-1] <= 0.0] = 1.0
    valid = output_edges[:-1] < incident_energy_keV
    answer = []
    for index in range(BASIS_SIZE):
        coefficients = np.zeros(BASIS_SIZE)
        coefficients[index] = 1.0
        antiderivative = BSpline(
            knots, coefficients, degree,
            extrapolate=False).antiderivative()
        column = np.zeros(output_edges.size - 1)
        column[valid] = (
            antiderivative(upper[valid])
            - antiderivative(lower[valid]))
        normalization = float(column.sum())
        if normalization > 0.0:
            column /= normalization
        answer.append(column)
    return np.column_stack(answer)


def truth_for(cell: Cell, name: str) -> tuple[np.ndarray, np.ndarray]:
    if name in ATOM_COMPONENTS:
        return (
            np.asarray([cell.atom_yield[name]]),
            np.asarray([cell.atom_variance[name]]))
    return cell.values[name], cell.variance[name]


def fit_shape_coefficients(
    cell: Cell,
    name: str,
    truth: np.ndarray,
    variance: np.ndarray,
    total: float,
) -> np.ndarray | None:
    """Poisson-weighted nonnegative least squares for the normalized shape.

    The weights are 1/sigma per bin with the zero-count bins floored at ONE
    count (1/events^2). Without that floor an empty bin gets variance ~0 and
    therefore the LARGEST weight in the fit -- the same defect found and fixed
    in fit_photonuclear.py.
    """
    support = (
        ELECTRON_MASS_KEV
        if name == "pair_annihilation_continuum"
        else cell.energy_keV)
    basis = warped_bspline_matrix_cached(cell.identity, cell.edges, support)
    shape = truth / total
    one_count = 1.0 / (cell.events * cell.events)
    floored = np.maximum(variance, one_count) / (total * total)
    weight = 1.0 / np.sqrt(np.maximum(floored, 1e-300))
    # scipy's nnls IS the exact Lawson-Hanson active-set algorithm.
    coefficients, _ = nnls(basis * weight[:, None], shape * weight)
    if coefficients.sum() <= 0.0:
        return None
    return coefficients / coefficients.sum()


def fit_nodes(
    cells: list[Cell],
    name: str,
    minimum_counts: float,
    anchor_starved: bool = True,
    anchor_first_only: bool = False,
) -> list[Node]:
    """One node per training cell.

    Cells whose truth is too sparse to fit a shape become ANCHOR nodes
    carrying the 90% CL Poisson upper bound (2.3 counts / N) on the yield, so
    the field is pulled down by data. Dropping them instead would leave the
    deep-tau region without nodes, and the kernel field would extrapolate flat
    off the last populated node, keeping spurious escaping yield in fully
    opaque shields. The shape of an
    anchor is taken from the cell's own mechanistic support (flat over the
    basis), which is irrelevant at these magnitudes but keeps the node valid.
    """
    nodes = []
    anchor_indices: list[int] = []
    anchors = 0
    # Anchoring EVERY starved cell puts a long run of near-zero nodes at the
    # deep end of each (E, Z, f_H) ray. A local-linear slope fitted through
    # that run extrapolates a steeper falloff than the truth, which is what
    # over-suppressed the brems channels between tau = 4 and 8. Keeping only
    # the shallowest starved cell per ray still pins where the field must come
    # down, without stacking the slope.
    allowed_anchor: set[str] | None = None
    if anchor_first_only:
        by_ray: dict[tuple, list[tuple[float, str]]] = {}
        for cell in cells:
            if cell.split != "train":
                continue
            truth, _ = truth_for(cell, name)
            if float(truth.sum()) * cell.events >= minimum_counts:
                continue
            key = (cell.energy_keV, cell.z_value, cell.hydrogen_mass_fraction)
            by_ray.setdefault(key, []).append((cell.tau, cell.identity))
        allowed_anchor = {
            sorted(entries)[0][1] for entries in by_ray.values()}
    for cell in cells:
        if cell.split != "train":
            continue
        truth, variance = truth_for(cell, name)
        total = float(truth.sum())
        base = component_base(cell, name)
        if base <= 0.0:
            continue
        starved = total * cell.events < minimum_counts
        if starved and not anchor_starved:
            continue
        if (starved and allowed_anchor is not None
                and cell.identity not in allowed_anchor):
            continue
        if starved:
            # The larger of what was actually seen and the zero-count bound:
            # never below the observation, never above the 90% CL ceiling a
            # genuinely empty cell allows.
            bound = max(total, POISSON_ZERO_UPPER_BOUND / cell.events)
            # A shield x-ray family whose base is below the bound (just under
            # an edge, or lines barely above the 10 keV floor) learns nothing
            # from an empty cell: the anchor would store bound / base >> 1 and
            # raise the field around it instead of pulling it down.
            if name in FLUORESCENCE_COMPONENTS and bound >= base:
                continue
            coefficients = np.full(BASIS_SIZE, 1.0 / BASIS_SIZE)
            anchor_indices.append(len(nodes))
            nodes.append(Node(
                feature=component_feature(cell, name),
                log_conditional_yield=math.log(max(bound / base, 1e-300)),
                shape_coefficients=coefficients,
            ))
            anchors += 1
            continue
        if name in ATOM_COMPONENTS:
            coefficients = np.full(BASIS_SIZE, 1.0 / BASIS_SIZE)
        else:
            fitted = fit_shape_coefficients(
                cell, name, truth, variance, total)
            if fitted is None:
                continue
            coefficients = fitted
        nodes.append(Node(
            feature=component_feature(cell, name),
            log_conditional_yield=math.log(max(total / base, 1e-300)),
            shape_coefficients=coefficients,
        ))
    ANCHOR_COUNTS[name] = anchors
    # An anchor constrains the YIELD field only. Leaving it with a flat shape
    # would drag the shape kernel toward a spectrum no cell actually has, so
    # each anchor inherits the shape of the nearest populated node instead.
    # That keeps the payload format unchanged and makes anchors shape-neutral.
    populated = [
        (index, node) for index, node in enumerate(nodes)
        if index not in anchor_indices]
    if populated:
        populated_features = np.asarray(
            [node.feature for _, node in populated])
        for index in anchor_indices:
            distance = np.sum(
                (populated_features - nodes[index].feature) ** 2, axis=1)
            nearest = populated[int(np.argmin(distance))][1]
            nodes[index].shape_coefficients = (
                nearest.shape_coefficients.copy())
    return nodes


def _local_linear(
    weight: np.ndarray,
    offset: np.ndarray,
    value: np.ndarray,
    ridge: float,
) -> float:
    """Intercept of a weighted linear fit of `value` against `offset`."""
    weight_sum = float(weight.sum())
    offset_sum = float(np.dot(weight, offset))
    offset2_sum = float(np.dot(weight, offset * offset))
    value_sum = float(np.dot(weight, value))
    value_offset_sum = float(np.dot(weight, offset * value))
    stabilized = offset2_sum + ridge * weight_sum
    determinant = (
        weight_sum * stabilized - offset_sum * offset_sum
        + np.finfo(float).tiny)
    return float(
        (value_sum * stabilized - value_offset_sum * offset_sum)
        / determinant)


def rbf_predict(
    cell: Cell,
    name: str,
    nodes: list[Node],
    yield_bandwidth: np.ndarray,
    shape_bandwidth: np.ndarray,
    local_linear_ridge: float = 0.0,
    local_linear_feature_index: int = 2,
    tail_decay_kappa: float = 0.0,
    tail_decay_tau: float = 0.0,
) -> np.ndarray:
    if not nodes:
        size = 1 if name in ATOM_COMPONENTS else cell.edges.size - 1
        return np.zeros(size)
    features = np.asarray([node.feature for node in nodes])
    query_feature = component_feature(cell, name)
    dy = (features - query_feature) / yield_bandwidth
    ds = (features - query_feature) / shape_bandwidth
    wy = np.exp(-0.5 * np.sum(dy * dy, axis=1))
    ws = np.exp(-0.5 * np.sum(ds * ds, axis=1))
    if wy.sum() <= 1e-300 or ws.sum() <= 1e-300:
        size = 1 if name in ATOM_COMPONENTS else cell.edges.size - 1
        return np.zeros(size)
    log_yields = np.asarray(
        [node.log_conditional_yield for node in nodes])
    shapes = np.asarray([node.shape_coefficients for node in nodes])
    if local_linear_ridge > 0.0:
        # Ridge-stabilized local LINEAR regression along the steep axis,
        # centred on the query so the intercept IS the prediction there.
        # Identical algebra to the C++ evaluator.
        axis = min(local_linear_feature_index, features.shape[1] - 1)
        uy = dy[:, axis]
        us = ds[:, axis]
        log_yield = _local_linear(wy, uy, log_yields, local_linear_ridge)
        coefficients = np.asarray([
            _local_linear(ws, us, shapes[:, k], local_linear_ridge)
            for k in range(shapes.shape[1])])
        # Local-linear coefficients can go negative where a kernel mean
        # cannot; the same smooth softplus the evaluator applies.
        width = 1.0e-12
        coefficients = 0.5 * (
            coefficients + np.sqrt(coefficients * coefficients + width * width))
    else:
        log_yield = float(np.dot(wy, log_yields) / wy.sum())
        coefficients = np.average(shapes, axis=0, weights=ws)
    if tail_decay_kappa > 0.0 and tail_decay_tau > 0.0:
        log_yield -= tail_decay_kappa * max(0.0, cell.tau - tail_decay_tau)
    total = component_base(cell, name) * math.exp(log_yield)
    if name in ATOM_COMPONENTS:
        return np.asarray([total])
    coefficients = np.maximum(coefficients, 0.0)
    coefficients /= max(float(coefficients.sum()), 1e-300)
    support = (
        ELECTRON_MASS_KEV
        if name == "pair_annihilation_continuum"
        else cell.energy_keV)
    return total * (
        warped_bspline_matrix_cached(cell.identity, cell.edges, support)
        @ coefficients)


def grouped_metrics(
    truth: np.ndarray,
    variance: np.ndarray,
    prediction: np.ndarray,
    groups: int = 32,
) -> dict[str, float]:
    usable_groups = min(groups, truth.size)
    boundaries = np.linspace(
        0, truth.size, usable_groups + 1, dtype=int)
    combine = lambda value: np.asarray([
        value[boundaries[i]:boundaries[i + 1]].sum()
        for i in range(usable_groups)])
    truth_g = combine(truth)
    variance_g = combine(variance)
    prediction_g = combine(prediction)
    truth_yield = float(truth_g.sum())
    predicted_yield = float(prediction_g.sum())
    if truth_yield <= 0.0:
        return {
            "truth_yield": truth_yield,
            "predicted_yield": predicted_yield,
            "yield_relative_error": math.inf,
            "shape_total_variation": math.inf,
            "weighted_rms_relative_error": math.inf,
            "p95_relative_error": math.inf,
        }
    relative_uncertainty = (
        np.sqrt(np.maximum(variance_g, 0.0))
        / np.maximum(truth_g, 1e-300))
    valid = (truth_g > 0.0) & (relative_uncertainty <= 0.05)
    relative = (
        np.abs(prediction_g - truth_g)
        / np.maximum(truth_g, 1e-300))
    if np.any(valid):
        wrms = math.sqrt(float(
            np.sum(truth_g[valid] * relative[valid] ** 2)
            / np.sum(truth_g[valid])))
        p95 = float(np.percentile(relative[valid], 95))
    else:
        wrms = math.inf
        p95 = math.inf
    predicted_norm = prediction_g / max(predicted_yield, 1e-300)
    return {
        "truth_yield": truth_yield,
        "predicted_yield": predicted_yield,
        "yield_relative_error": (
            predicted_yield - truth_yield) / truth_yield,
        "shape_total_variation": 0.5 * float(np.abs(
            truth_g / truth_yield - predicted_norm).sum()),
        "weighted_rms_relative_error": wrms,
        "p95_relative_error": p95,
    }


def evaluate(
    cells: list[Cell],
    node_map: dict[str, list[Node]],
    yield_bandwidth: np.ndarray,
    shape_bandwidth: np.ndarray,
    knobs: dict[str, float] | None = None,
    include_train: bool = False,
) -> tuple[list[dict[str, float | str]], dict[tuple[str, str], np.ndarray]]:
    knobs = knobs or {}
    tail_tau = {
        name: max((float(node.feature[2]) for node in nodes), default=0.0)
        for name, nodes in node_map.items()}
    scale = math.log(25.0)
    records = []
    predictions = {}
    for cell in cells:
        if cell.split == "train" and not include_train:
            continue
        for name in ALL_COMPONENTS:
            truth, variance = truth_for(cell, name)
            prediction = rbf_predict(
                cell, name, node_map[name],
                yield_bandwidth, shape_bandwidth,
                local_linear_ridge=knobs.get("local_linear_ridge", 0.0),
                local_linear_feature_index=int(
                    knobs.get("local_linear_feature_index", 2)),
                tail_decay_kappa=knobs.get("tail_decay_kappa", 0.0),
                tail_decay_tau=(
                    math.expm1(tail_tau[name] * scale)
                    if tail_tau[name] > 0.0 else 0.0))
            predictions[(cell.identity, name)] = prediction
            if float(truth.sum()) * cell.events < 25.0:
                continue
            metrics = grouped_metrics(truth, variance, prediction)
            records.append({
                "identity": cell.identity,
                "split": cell.split,
                "component": name,
                "energy_keV": cell.energy_keV,
                "Z": cell.z_value,
                "tau": cell.tau,
                **metrics,
            })
    return records, predictions


# Interior tau nodes of the training grid.  Holding one out and predicting it
# from the rest is a BETWEEN-node criterion built entirely from training data,
# which is what the bandwidth scan needs and what it did not have.
CROSS_VALIDATION_TAU = (0.5, 2.0)


def tau_cross_validation_folds(
    cells: list[Cell],
    minimum_counts: float,
    anchor_starved: bool,
    anchor_first_only: bool,
) -> list[tuple[dict[str, list[Node]], list[Cell]]]:
    """(node_map without one interior tau, the cells at that tau) per fold."""
    folds = []
    for tau in CROSS_VALIDATION_TAU:
        kept = [cell for cell in cells if cell.tau != tau]
        held = [cell for cell in cells if cell.tau == tau]
        if not held or not kept:
            continue
        node_map = {
            name: fit_nodes(kept, name, minimum_counts,
                            anchor_starved=anchor_starved,
                            anchor_first_only=anchor_first_only)
            for name in ALL_COMPONENTS}
        folds.append((node_map, held))
    return folds


def scan_bandwidths(
    cells: list[Cell],
    node_map: dict[str, list[Node]],
    base_yield: np.ndarray,
    base_shape: np.ndarray,
    knobs: dict[str, float],
    folds: list[tuple[dict[str, list[Node]], list[Cell]]],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    """Pick the RBF bandwidths instead of trusting hand-set constants.

    Scored on BOTH between-node generalization and on-grid training-cell
    reproduction, combined with a max so neither can be traded away: scoring
    interpolation consistency alone once let a degenerate over-wide kernel win,
    because two equally-wrong predictions agree with each other perfectly.

    The generalization term is a leave-one-interior-tau-node-out cross
    validation built from the TRAINING grid. Z-holdout cells sit ON the
    training tau grid, so they cannot measure between-node interpolation; the
    off-grid-tau bank can, which is why it stays report-only rather than being
    spent on model selection.

    The yield and shape kernels are independent in the estimator (separate
    weight sets), so they are scanned separately -- 2*N evaluations rather
    than N^2.
    """
    global_scales = (0.35, 0.5, 0.75, 1.0, 1.5, 2.0)
    process_scales = (0.5, 1.0, 2.0)
    trace: list[dict[str, float]] = []

    def candidates(base: np.ndarray) -> list[np.ndarray]:
        answer = []
        for scale in global_scales:
            for process in process_scales:
                trial = base * scale
                if trial.size > 3:
                    trial = trial.copy()
                    trial[3:] *= process
                answer.append(trial)
                if trial.size <= 3:
                    break
        return answer

    def score(yield_bandwidth: np.ndarray, shape_bandwidth: np.ndarray,
              metric: str) -> float:
        train = [
            abs(float(r[metric]))
            for r in evaluate(cells, node_map, yield_bandwidth,
                              shape_bandwidth, knobs, include_train=True)[0]
            if r["split"] == "train" and math.isfinite(float(r[metric]))]
        between = []
        for fold_nodes, held in folds:
            between += [
                abs(float(r[metric]))
                for r in evaluate(held, fold_nodes, yield_bandwidth,
                                  shape_bandwidth, knobs,
                                  include_train=True)[0]
                if math.isfinite(float(r[metric]))]
        if not between or not train:
            return math.inf
        return max(float(np.median(between)), float(np.median(train)))

    best_yield, best_yield_score = base_yield, math.inf
    for trial in candidates(base_yield):
        value = score(trial, base_shape, "yield_relative_error")
        trace.append({
            "kind": "yield", "bandwidth": [float(v) for v in trial],
            "score": value})
        if value < best_yield_score:
            best_yield, best_yield_score = trial, value

    best_shape, best_shape_score = base_shape, math.inf
    for trial in candidates(base_shape):
        value = score(best_yield, trial, "shape_total_variation")
        trace.append({
            "kind": "shape", "bandwidth": [float(v) for v in trial],
            "score": value})
        if value < best_shape_score:
            best_shape, best_shape_score = trial, value

    return best_yield, best_shape, trace


def finite_summary(values: list[float]) -> dict[str, float | int]:
    finite = np.asarray([value for value in values if math.isfinite(value)])
    if finite.size == 0:
        return {"count": 0}
    return {
        "count": int(finite.size),
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95)),
        "maximum": float(np.max(finite)),
    }


def write_payload(
    path: Path,
    node_map: dict[str, list[Node]],
    yield_bandwidth: np.ndarray,
    shape_bandwidth: np.ndarray,
    local_linear_ridge: float = 0.0,
    local_linear_feature_index: int = 2,
    tail_decay_kappa: float = 0.0,
) -> None:
    with path.open("w") as stream:
        stream.write("CEELO_SOURCE_SECONDARY_SURROGATE 1\n")
        stream.write(f"inverse_energy_warp_scale {WARP_SCALE:.17g}\n")
        stream.write("absorption_envelope 0 3.5 10\n")
        stream.write("basis cubic_bspline\n")
        stream.write(
            f"local_linear {float(local_linear_ridge):.17g} "
            f"{int(local_linear_feature_index)}\n")
        stream.write(f"tail_decay {float(tail_decay_kappa):.17g}\n")
        for name in ALL_COMPONENTS:
            nodes = node_map[name]
            stream.write(
                f"{name} {BASIS_SIZE - 1} "
                + " ".join(f"{value:.17g}" for value in yield_bandwidth)
                + " "
                + " ".join(f"{value:.17g}" for value in shape_bandwidth)
                + f" {len(nodes)}\n")
            for node in nodes:
                stream.write(
                    f"node {name} "
                    + " ".join(f"{value:.17g}" for value in node.feature)
                    + f" {node.log_conditional_yield:.17g} "
                    + " ".join(
                        f"{value:.17g}"
                        for value in node.shape_coefficients)
                    + "\n")
        stream.write("END\n")


def write_metrics(path: Path, records: list[dict[str, float | str]]) -> None:
    if not records:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_plots(
    output: Path,
    cells: list[Cell],
    records: list[dict[str, float | str]],
    predictions: dict[tuple[str, str], np.ndarray],
    worst: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    by_identity = {cell.identity: cell for cell in cells}
    ranked = sorted(
        records,
        key=lambda record: abs(float(record["yield_relative_error"]))
        + float(record["shape_total_variation"]),
        reverse=True)
    with PdfPages(output / "source_secondary_worst_comparisons.pdf") as pdf:
        for record in ranked[:worst]:
            cell = by_identity[str(record["identity"])]
            name = str(record["component"])
            truth, variance = truth_for(cell, name)
            prediction = predictions[(cell.identity, name)]
            fig, axes = plt.subplots(
                2, 1, figsize=(9.0, 6.8), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]})
            if name in ATOM_COMPONENTS:
                # A fluorescence family has many lines; one marker stands
                # for their sum.
                atom_energy = (ELECTRON_MASS_KEV if name == "pair_511_atom"
                               else cell.energy_keV / 2.0)
                axes[0].bar(
                    [atom_energy], truth,
                    width=max(cell.energy_keV / 128.0, 2.0),
                    alpha=0.45, label="CeeLo truth")
                axes[0].scatter(
                    [atom_energy], prediction,
                    color="tab:orange", label="functional model")
                axes[1].axhline(0.0, color="black", linewidth=0.8)
                axes[1].scatter(
                    [atom_energy],
                    [(prediction[0] - truth[0])
                     / max(truth[0], 1e-300)])
            else:
                centers = 0.5 * (cell.edges[:-1] + cell.edges[1:])
                widths = np.diff(cell.edges)
                axes[0].step(
                    centers, truth / widths, where="mid",
                    label="CeeLo truth")
                axes[0].fill_between(
                    centers,
                    np.maximum(
                        0.0, truth - np.sqrt(np.maximum(variance, 0.0)))
                    / widths,
                    (truth + np.sqrt(np.maximum(variance, 0.0))) / widths,
                    step="mid", alpha=0.2)
                axes[0].step(
                    centers, prediction / widths, where="mid",
                    label="functional model")
                residual = (
                    prediction - truth) / np.maximum(truth, 1e-12)
                axes[1].step(
                    centers, residual, where="mid",
                    color="tab:red")
                axes[1].set_ylim(-1.5, 1.5)
            axes[0].set_yscale("log")
            axes[0].set_ylabel("escaping photons / keV / primary")
            axes[0].legend()
            axes[1].set_ylabel("(model-truth)/truth")
            axes[1].set_xlabel("escaping photon energy (keV)")
            fig.suptitle(
                f"{name}: E0={cell.energy_keV:g} keV, "
                f"Z={cell.z_value}, tau={cell.tau:g}\n"
                f"yield error={100*float(record['yield_relative_error']):.2f}%, "
                f"TV={100*float(record['shape_total_variation']):.2f}%")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    # A second PDF shows the physically relevant sum of all non-fluorescence
    # secondary photons, including the 511-keV atom placed in its source bin.
    candidates = [
        cell for cell in cells
        if cell.split != "train"
        and any((cell.identity, name) in predictions for name in ALL_COMPONENTS)
    ]
    total_ranked = []
    for cell in candidates:
        truth = sum(cell.values.values(), np.zeros(cell.edges.size - 1))
        prediction = np.zeros_like(truth)
        for name in CONTINUOUS_COMPONENTS:
            prediction += predictions.get(
                (cell.identity, name), np.zeros_like(truth))
        pair_index = int(np.searchsorted(
            cell.edges, ELECTRON_MASS_KEV, side="right") - 1)
        if 0 <= pair_index < truth.size:
            truth[pair_index] += cell.atom_yield["pair_511_atom"]
            atom_prediction = predictions.get(
                (cell.identity, "pair_511_atom"), np.zeros(1))
            prediction[pair_index] += atom_prediction[0]
        score = float(np.abs(prediction - truth).sum()) / max(
            float(truth.sum()), 1e-300)
        total_ranked.append((score, cell, truth, prediction))
    total_ranked.sort(key=lambda value: value[0], reverse=True)
    with PdfPages(output / "source_secondary_sum_comparisons.pdf") as pdf:
        for score, cell, truth, prediction in total_ranked[:worst]:
            centers = 0.5 * (cell.edges[:-1] + cell.edges[1:])
            widths = np.diff(cell.edges)
            fig, axes = plt.subplots(
                2, 1, figsize=(9.0, 6.8), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]})
            axes[0].step(
                centers, truth / widths, where="mid", label="CeeLo truth")
            axes[0].step(
                centers, prediction / widths, where="mid",
                label="functional sum")
            axes[0].set_yscale("log")
            axes[0].set_ylabel("escaping photons / keV / primary")
            axes[0].legend()
            residual = (
                prediction - truth) / np.maximum(truth, 1e-12)
            axes[1].step(
                centers, residual, where="mid", color="tab:red")
            axes[1].set_ylim(-1.5, 1.5)
            axes[1].set_ylabel("relative residual")
            axes[1].set_xlabel("escaping photon energy (keV)")
            fig.suptitle(
                f"all non-fluorescence secondaries: "
                f"E0={cell.energy_keV:g} keV, Z={cell.z_value}, "
                f"tau={cell.tau:g}; L1/yield={100*score:.2f}%")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--offgrid-campaign", type=Path, action="append", default=[],
        help="cells held out entirely and scored as their own split; every "
             "tau in them sits BETWEEN training nodes, which whole-cell Z "
             "holdouts (sharing the training tau grid) cannot probe")
    parser.add_argument(
        "--anchor-first-only", action="store_true",
        help="anchor only the shallowest starved cell of each (E, Z, f_H) "
             "ray instead of every one")
    parser.add_argument(
        "--no-anchor-starved", action="store_true",
        help="drop sparse cells instead of anchoring them at the 90%% CL "
             "Poisson upper bound")
    parser.add_argument("--local-linear-ridge", type=float, default=0.0)
    parser.add_argument(
        "--local-linear-feature-index", type=int, default=2)
    parser.add_argument("--tail-decay-kappa", type=float, default=0.0)
    parser.add_argument(
        "--cross-section-helper", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_cross_section_features"))
    parser.add_argument("--holdout-z", default="13,35,64")
    parser.add_argument("--minimum-counts", type=float, default=25.0)
    parser.add_argument(
        "--yield-bandwidth",
        default="0.08,0.12,0.08,0.08,0.10,0.08,0.06")
    parser.add_argument(
        "--shape-bandwidth",
        default="0.12,0.14,0.12,0.12,0.14,0.12,0.10")
    parser.add_argument("--plot-worst", type=int, default=20)
    parser.add_argument(
        "--scan-bandwidths", action="store_true",
        help="select the RBF bandwidths by scan instead of using the "
             "hand-set defaults")
    args = parser.parse_args()

    holdout_z = {
        int(value) for value in args.holdout_z.split(",") if value}
    yield_bandwidth = np.asarray([
        float(value) for value in args.yield_bandwidth.split(",")])
    shape_bandwidth = np.asarray([
        float(value) for value in args.shape_bandwidth.split(",")])
    if (yield_bandwidth.size != PROCESS_FEATURE_SIZE
            or shape_bandwidth.size != PROCESS_FEATURE_SIZE):
        raise ValueError(
            f"RBF bandwidths need {PROCESS_FEATURE_SIZE} comma-separated "
            "values")

    def collect(root: Path) -> list[Path]:
        # An empty return means the bank is not there.  Callers MUST raise on
        # it rather than fitting a smaller problem: the off-grid bank was
        # deleted once and this silently dropped the only between-node holdout
        # instead of failing.
        found = bank_io.cell_files(root / "cells", recursive=True)
        return found or bank_io.cell_files(root / "cells")

    cells = []
    for root in args.campaign:
        found = collect(root)
        if not found:
            raise ValueError(f"no source cells under {root}")
        cells.extend(load_cell(path, holdout_z) for path in found)
    offgrid_cells = []
    for root in args.offgrid_campaign:
        found = collect(root)
        if not found:
            raise ValueError(f"no off-grid cells under {root}")
        for path in found:
            cell = load_cell(path, holdout_z)
            cell.split = "offgrid"
            offgrid_cells.append(cell)
    every_cell = cells + offgrid_cells
    attach_component_bases(every_cell, args.cross_section_helper)
    attach_process_features(every_cell, args.cross_section_helper)
    node_map = {
        name: fit_nodes(
            cells, name, args.minimum_counts,
            anchor_starved=not args.no_anchor_starved,
            anchor_first_only=args.anchor_first_only)
        for name in ALL_COMPONENTS}
    knobs = {
        "local_linear_ridge": args.local_linear_ridge,
        "local_linear_feature_index": args.local_linear_feature_index,
        "tail_decay_kappa": args.tail_decay_kappa,
    }
    bandwidth_trace: list[dict[str, float]] = []
    if args.scan_bandwidths:
        folds = tau_cross_validation_folds(
            cells, args.minimum_counts,
            anchor_starved=not args.no_anchor_starved,
            anchor_first_only=args.anchor_first_only)
        print(
            "bandwidth scan: leave-one-interior-tau-out folds at "
            + ", ".join(f"tau={tau:g} ({len(held)} cells)"
                        for (_, held), tau in zip(folds, CROSS_VALIDATION_TAU))
            + "; the off-grid bank is report-only")
        yield_bandwidth, shape_bandwidth, bandwidth_trace = scan_bandwidths(
            cells, node_map, yield_bandwidth, shape_bandwidth, knobs, folds)
        print(
            "selected yield bandwidth "
            + ",".join(f"{v:.4g}" for v in yield_bandwidth)
            + "  shape bandwidth "
            + ",".join(f"{v:.4g}" for v in shape_bandwidth))
    records, predictions = evaluate(
        every_cell, node_map, yield_bandwidth, shape_bandwidth, knobs)
    ongrid_records, _ = evaluate(
        cells, node_map, yield_bandwidth, shape_bandwidth, knobs,
        include_train=True)
    ongrid_records = [
        record for record in ongrid_records if record["split"] == "train"]

    args.output.mkdir(parents=True, exist_ok=True)
    write_payload(
        args.output / "secondary_surrogate.txt",
        node_map, yield_bandwidth, shape_bandwidth,
        local_linear_ridge=args.local_linear_ridge,
        local_linear_feature_index=args.local_linear_feature_index,
        tail_decay_kappa=args.tail_decay_kappa)
    write_metrics(args.output / "metrics.csv", records + ongrid_records)
    if args.plot_worst > 0:
        write_plots(
            args.output, every_cell, records, predictions, args.plot_worst)

    report = {
        "model": "positive_inverse_energy_bspline_rbf",
        "fluorescence_included": True,
        "basis_size": BASIS_SIZE,
        "yield_bandwidth": yield_bandwidth.tolist(),
        "shape_bandwidth": shape_bandwidth.tolist(),
        "cells": len(every_cell),
        "local_linear_ridge": args.local_linear_ridge,
        "local_linear_feature_index": args.local_linear_feature_index,
        "tail_decay_kappa": args.tail_decay_kappa,
        "anchor_starved_cells": not args.no_anchor_starved,
        "anchor_first_only": args.anchor_first_only,
        "anchor_nodes": dict(ANCHOR_COUNTS),
        "bandwidth_scan": bandwidth_trace,
        "splits": {
            split: sum(cell.split == split for cell in every_cell)
            for split in ("train", "test", "offgrid")},
        "nodes": {
            name: len(nodes) for name, nodes in node_map.items()},
        "metrics": {},
        "accuracy": {},
    }
    for name in ALL_COMPONENTS:
        selected = [
            record for record in records if record["component"] == name]
        report["metrics"][name] = {
            metric: finite_summary([
                abs(float(record[metric])) for record in selected])
            for metric in (
                "yield_relative_error",
                "shape_total_variation",
                "weighted_rms_relative_error",
                "p95_relative_error")
        }
    # The three-way accuracy split the payload README must quote: on-grid
    # (training cells reproduced), off-grid tau, and held-out element.
    for label, subset in (
            ("ongrid", ongrid_records),
            ("offgrid_tau", [
                r for r in records if r["split"] == "offgrid"]),
            ("heldout_element", [
                r for r in records if r["split"] == "test"])):
        report["accuracy"][label] = {
            metric: finite_summary([
                abs(float(record[metric])) for record in subset])
            for metric in (
                "yield_relative_error", "shape_total_variation")
        }
    (args.output / "report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
