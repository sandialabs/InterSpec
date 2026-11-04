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

"""Fit positive functional residuals for high-order source scattering.

Rayleigh remains active in source transport, but interaction order is defined
only by the number of Compton collisions. The exact low-order C++ transport
handles two-Compton histories and this script fits only the remaining
three-or-more-Compton residual. Each
spectrum is represented by either a
positive cubic B-spline or positive Bernstein density in an additive
inverse-energy Compton-loss coordinate. A broadened rational warp resolves
the absorption-shaped low-energy turn-on while preserving exact support down
to zero. Smooth normalized RBF functions map
(log E, effective Z, log optical depth) to a log yield and positive basis
coefficients; no histogram is interpolated at evaluation time.
"""

from __future__ import annotations

import argparse
import csv
import bank_io
import hashlib
import io
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import nnls
from scipy.interpolate import BSpline
from scipy.special import betainc

from compare_reference import (
    rebin,
    rebin_variance,
    resolve,
    resolve_variance,
)
from validate_single_scatter import calculate_metrics


ELECTRON_MASS_KEV = 510.99895
INVERSE_ENERGY_WARP_SCALE = 2.0
ABSORPTION_ENVELOPE_STRENGTH = 0.0
ABSORPTION_ENVELOPE_POWER = 3.5
ABSORPTION_ENVELOPE_MAX_RATIO = 10.0
RESIDUAL_COMPONENTS = ("higher_aggregate",)
RAW_RESIDUAL_COMPONENTS = (
    "higher_compton",
    "higher_compton_rayleigh",
)
FALLBACK_SCALES = {
    "higher_compton": 0.61,
    "higher_compton_rayleigh": 0.20,
}
SOLVED_COMPONENTS = (
    "double_compton",
    "double_compton_rayleigh",
)


@dataclass
class Cell:
    identity: str
    energy_keV: float
    z_value: int
    tau: float
    areal_density: float
    events: int
    feature: np.ndarray
    basis: np.ndarray
    truth: dict[str, np.ndarray]
    variance: dict[str, np.ndarray]
    model: dict[str, np.ndarray]
    split: str


def feature(energy_keV: float, z_value: int, tau: float) -> np.ndarray:
    return np.array([
        math.log(energy_keV / 10.0) / math.log(1200.0),
        z_value / 92.0,
        math.log1p(tau) / math.log(25.0),
    ])


def process_feature(
    energy_keV: float,
    z_value: int,
    tau: float,
    areal_density: float,
    cross_sections: dict[str, float],
) -> np.ndarray:
    """Physics coordinates for the changing multiple-scatter mixture."""
    scale = math.log(25.0)
    return np.asarray([
        math.log(energy_keV / 10.0) / math.log(1200.0),
        z_value / 92.0,
        math.log1p(tau) / scale,
        math.log1p(areal_density * cross_sections["mu_cs"]) / scale,
        math.log1p(areal_density * cross_sections["mu_rs"]) / scale,
        math.log1p(cross_sections["tau_abs3"]) / scale,
        math.log1p(cross_sections["tau_abs8"]) / scale,
        math.log1p(areal_density * cross_sections["mu_pp"]) / scale,
    ])


def yield_feature(value: np.ndarray, tau: float) -> np.ndarray:
    """Remove fitted thickness from process coordinates for yield curves."""
    if value.size == 3:
        return value[:2]
    if value.size != 8 or tau <= 0.0:
        raise ValueError("invalid process feature for a tau-yield node")
    scale = math.log(25.0)
    process_tau = np.expm1(value[3:] * scale)
    process_per_total_tau = process_tau / tau
    return np.concatenate((
        value[:2],
        np.log1p(process_per_total_tau) / scale))


def warped_bernstein_matrix(
    energy_keV: float,
    bins: int,
    degree: int,
) -> np.ndarray:
    edges = np.linspace(0.0, energy_keV, bins + 1)

    def warped(energy: np.ndarray) -> np.ndarray:
        safe = np.maximum(energy, 1e-100)
        loss = 0.5 * ELECTRON_MASS_KEV * (
            1.0 / safe - 1.0 / energy_keV)
        return loss / (INVERSE_ENERGY_WARP_SCALE + loss)

    upper_warp = warped(edges[:-1])
    lower_warp = warped(edges[1:])
    upper_warp[0] = 1.0
    return np.column_stack([
        betainc(index + 1, degree - index + 1, upper_warp)
        - betainc(index + 1, degree - index + 1, lower_warp)
        for index in range(degree + 1)
    ])


def warped_cubic_bspline_matrix(
    energy_keV: float,
    bins: int,
    basis_size: int,
) -> np.ndarray:
    spline_degree = 3
    if basis_size <= spline_degree:
        raise ValueError("cubic B-spline basis is too small")
    edges = np.linspace(0.0, energy_keV, bins + 1)

    def warped(energy: np.ndarray) -> np.ndarray:
        safe = np.maximum(energy, 1e-100)
        loss = 0.5 * ELECTRON_MASS_KEV * (
            1.0 / safe - 1.0 / energy_keV)
        return loss / (INVERSE_ENERGY_WARP_SCALE + loss)

    upper_warp = warped(edges[:-1])
    lower_warp = warped(edges[1:])
    upper_warp[0] = 1.0
    internal_count = basis_size - spline_degree - 1
    knots = np.concatenate((
        np.zeros(spline_degree + 1),
        np.linspace(0.0, 1.0, internal_count + 2)[1:-1],
        np.ones(spline_degree + 1),
    ))
    answer = []
    for index in range(basis_size):
        coefficients = np.zeros(basis_size)
        coefficients[index] = 1.0
        antiderivative = BSpline(
            knots, coefficients, spline_degree,
            extrapolate=False).antiderivative()
        column = (
            antiderivative(upper_warp)
            - antiderivative(lower_warp))
        column /= max(float(column.sum()), 1e-300)
        answer.append(column)
    return np.column_stack(answer)


def warped_basis_matrix(
    energy_keV: float,
    bins: int,
    degree: int,
    basis_type: str,
) -> np.ndarray:
    if basis_type == "bernstein":
        return warped_bernstein_matrix(
            energy_keV, bins, degree)
    if basis_type == "cubic_bspline":
        return warped_cubic_bspline_matrix(
            energy_keV, bins, degree + 1)
    raise ValueError(f"unsupported basis type: {basis_type}")


def load_histograms(
    path: Path,
    component_names: tuple[str, ...],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, float]]:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty cell: {path}")
    rayleigh_enabled = rows[0].get("rayleigh_enabled")
    if (rayleigh_enabled is not None
            and rayleigh_enabled.strip().lower()
            not in {"1", "true", "yes"}):
        raise ValueError(
            f"{path} is not a Rayleigh-enabled source-study cell")
    energy = float(rows[0]["energy_keV"])
    events = float(rows[0]["events"])
    continuous = [
        row for row in rows
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])]
    bins = max(
        int(round(float(row["bin_hi_keV"]) / energy))
        for row in continuous)
    # The preceding expression gives a fraction, not the actual number of
    # bins. Count the rows for any complete continuous component instead.
    bin_counts = {
        name: sum(row["component"] == name for row in continuous)
        for name in component_names}
    bins = max(bin_counts.values(), default=0)
    if bins == 0:
        raise ValueError(f"no requested components in {path}")
    values = {name: np.zeros(bins) for name in component_names}
    variance = {name: np.zeros(bins) for name in component_names}
    fine_to_aggregate = {
        "triple_compton": "higher_compton",
        "quadruple_compton": "higher_compton",
        "five_plus_compton": "higher_compton",
        "triple_compton_rayleigh": "higher_compton_rayleigh",
        "quadruple_compton_rayleigh": "higher_compton_rayleigh",
        "five_plus_compton_rayleigh": "higher_compton_rayleigh",
    }
    # Fine-order diagnostics replace, rather than duplicate, the aggregate 3+
    # labels in the campaign driver.  Consolidate them on input so
    # research campaigns with diagnostics enabled remain numerically
    # compatible with the aggregate model components.
    fine_order_mode = any(
        row["component"] in fine_to_aggregate
        and float(row["sum_w"]) != 0.0
        for row in continuous)
    for row in continuous:
        name = row["component"]
        if fine_order_mode:
            name = fine_to_aggregate.get(name, name)
        if name not in values:
            continue
        index = min(
            bins - 1,
            int(round(float(row["bin_lo_keV"]) / energy * bins)))
        values[name][index] += float(row["sum_w"]) / events
        variance[name][index] += float(row["sum_w2"]) / (events * events)
    metadata = {
        "energy_keV": energy,
        "Z": float(rows[0]["Z"]),
        "tau": float(rows[0]["tau"]),
        "areal_density_g_cm2": float(rows[0]["areal_density_g_cm2"]),
        "events": events,
        "bins": bins,
        "rayleigh_enabled": 1.0,
    }
    return values, variance, metadata


def split_for(identity: str, z_value: int, holdout_z: set[int]) -> str:
    if z_value in holdout_z:
        return "test"
    bucket = int.from_bytes(
        hashlib.sha256(identity.encode()).digest()[:4], "big") % 5
    return "validation" if bucket == 0 else "train"


def load_cells(
    campaign: Path,
    degree: int,
    holdout_z: set[int],
    cross_section_helper: Path | None = None,
    absorption_envelope_strength: float = 0.0,
    absorption_envelope_power: float = ABSORPTION_ENVELOPE_POWER,
    absorption_envelope_max_ratio: float =
        ABSORPTION_ENVELOPE_MAX_RATIO,
    basis_type: str = "bernstein",
    feature_mode: str = "raw",
) -> list[Cell]:
    names = SOLVED_COMPONENTS + RAW_RESIDUAL_COMPONENTS
    cells = []
    for truth_path in bank_io.cell_files(campaign / "cells"):
        model_path = bank_io.companion(truth_path, "ceelo", "model")
        if model_path is None:
            continue
        truth, variance, metadata = load_histograms(truth_path, names)
        model, _, model_metadata = load_histograms(model_path, names)
        if int(metadata["bins"]) != int(model_metadata["bins"]):
            raise ValueError(f"bin mismatch for {truth_path}")
        identity = bank_io.cell_id(truth_path)
        z_value = int(metadata["Z"])
        cells.append(Cell(
            identity=identity,
            energy_keV=metadata["energy_keV"],
            z_value=z_value,
            tau=metadata["tau"],
            areal_density=metadata["areal_density_g_cm2"],
            events=int(metadata["events"]),
            feature=feature(
                metadata["energy_keV"], z_value, metadata["tau"]),
            basis=warped_basis_matrix(
                metadata["energy_keV"], int(metadata["bins"]),
                degree, basis_type),
            truth=truth,
            variance=variance,
            model=model,
            split=split_for(identity, z_value, holdout_z),
        ))
    if not cells:
        raise ValueError(f"no paired cells under {campaign}")
    if absorption_envelope_strength > 0.0 or feature_mode == "process":
        if cross_section_helper is None:
            raise ValueError(
                "process features or an absorption envelope require a "
                "cross-section helper")
        keys = sorted({
            (cell.energy_keV, cell.z_value) for cell in cells
        })
        result = subprocess.run(
            [str(cross_section_helper)],
            input="".join(
                f"{energy:.17g} {z_value} 1\n"
                for energy, z_value in keys),
            check=True, text=True, stdout=subprocess.PIPE)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        if len(rows) != len(keys):
            raise ValueError(
                "cross-section helper returned the wrong row count")
        cross_section_rows = {}
        for key, row in zip(keys, rows):
            if int(row["Z"]) != key[1] or not math.isclose(
                    float(row["energy_keV"]), key[0],
                    rel_tol=1e-12, abs_tol=1e-10):
                raise ValueError(
                    "cross-section helper changed row ordering")
            cross_section_rows[key] = {
                name: float(row[name])
                for name in (
                    "mu_pe", "mu_cs", "mu_rs", "mu_pp",
                    "tau_abs3", "tau_scatter3",
                    "tau_abs8", "tau_scatter8")}
        for cell in cells:
            row = cross_section_rows[
                (cell.energy_keV, cell.z_value)]
            if feature_mode == "process":
                cell.feature = process_feature(
                    cell.energy_keV, cell.z_value, cell.tau,
                    cell.areal_density, row)
            if absorption_envelope_strength <= 0.0:
                continue
            centers = (
                np.arange(cell.basis.shape[0], dtype=float) + 0.5
            ) * cell.energy_keV / cell.basis.shape[0]
            loss_factor = (
                np.minimum(
                    absorption_envelope_max_ratio,
                    cell.energy_keV / centers)
                ** absorption_envelope_power - 1.0)
            envelope = np.exp(-np.minimum(
                745.0,
                absorption_envelope_strength
                * cell.areal_density
                * row["mu_pe"]
                * loss_factor))
            cell.basis *= envelope[:, None]
            cell.basis /= np.maximum(
                cell.basis.sum(axis=0, keepdims=True), 1e-300)
    return cells


def per_cell_coefficients(
    cell: Cell,
    component: str,
    minimum_histories: int,
) -> tuple[float, np.ndarray] | None:
    if component != "higher_aggregate":
        raise ValueError(f"unsupported fitted component: {component}")
    values = sum(
        (cell.truth[name] for name in RAW_RESIDUAL_COMPONENTS),
        np.zeros_like(cell.truth[RAW_RESIDUAL_COMPONENTS[0]]))
    histories = values.sum() * cell.events
    base_yield = sum(
        cell.model[name].sum() / FALLBACK_SCALES[name]
        for name in RAW_RESIDUAL_COMPONENTS)
    if histories < minimum_histories or not base_yield > 0.0:
        return None
    coefficients = nnls(cell.basis, values)[0]
    fitted_yield = coefficients.sum()
    if not fitted_yield > 0.0:
        return None
    shape = coefficients / fitted_yield
    floor = max(1e-10 / coefficients.size, 1e-30)
    return (
        math.log(fitted_yield / base_yield),
        np.maximum(shape, floor))


@dataclass
class ComponentFit:
    component: str
    degree: int
    yield_bandwidth: np.ndarray
    shape_bandwidth: np.ndarray
    centers: np.ndarray
    log_conditional_yields: np.ndarray
    shape_coefficients: np.ndarray

    def predict_conditional(self, value: np.ndarray) -> np.ndarray:
        yield_distance = 0.5 * np.sum(
            ((self.centers - value) / self.yield_bandwidth) ** 2,
            axis=1)
        shape_distance = 0.5 * np.sum(
            ((self.centers - value) / self.shape_bandwidth) ** 2,
            axis=1)
        # Subtracting the nearest-node distance leaves normalized Gaussian
        # weights unchanged and prevents underflow in strongly absorbing
        # process-feature regimes.
        yield_weights = np.exp(-(yield_distance - yield_distance.min()))
        shape_weights = np.exp(-(shape_distance - shape_distance.min()))
        log_yield = (
            yield_weights @ self.log_conditional_yields
            / yield_weights.sum())
        # A convex average preserves positivity and interpolates through sparse
        # Bernstein coefficients without the artificial narrow spikes caused
        # by geometric (log-coefficient) interpolation near genuine zeros.
        shape = (
            shape_weights @ self.shape_coefficients
            / shape_weights.sum())
        shape /= shape.sum()
        return math.exp(log_yield) * shape

    def predict_shape(self, value: np.ndarray) -> np.ndarray:
        distance = 0.5 * np.sum(
            ((self.centers - value) / self.shape_bandwidth) ** 2,
            axis=1)
        weights = np.exp(-(distance - distance.min()))
        shape = weights @ self.shape_coefficients / weights.sum()
        return shape / shape.sum()


@dataclass
class TauYieldModel:
    degree: int
    bandwidth: np.ndarray
    centers: np.ndarray
    coefficients: np.ndarray
    minimum_tau: float = 0.2
    maximum_tau: float = 12.0
    use_process_features: bool = False

    def tau_coordinate(self, tau: float) -> float:
        lo = math.log1p(self.minimum_tau)
        hi = math.log1p(self.maximum_tau)
        return 2.0 * (math.log1p(tau) - lo) / (hi - lo) - 1.0

    def predict(self, value: np.ndarray, tau: float) -> float:
        value = (
            yield_feature(value, tau)
            if self.use_process_features else value[:2])
        distance = 0.5 * np.sum(
            ((self.centers - value) / self.bandwidth) ** 2,
            axis=1)
        weights = np.exp(-(distance - distance.min()))
        coefficients = weights @ self.coefficients / weights.sum()
        log_reduced_yield = np.polynomial.chebyshev.chebval(
            self.tau_coordinate(tau), coefficients)
        return tau ** 3 * math.exp(log_reduced_yield)


def fit_tau_yield_node(
    cells: list[Cell],
    degree: int,
    minimum_histories: int,
) -> np.ndarray | None:
    usable = []
    for cell in cells:
        values = sum(
            (cell.truth[name] for name in RAW_RESIDUAL_COMPONENTS),
            np.zeros(cell.basis.shape[0]))
        histories = values.sum() * cell.events
        # A zero-count cell is an upper-limit measurement, not missing data.
        # Jeffreys-like half-count regularization keeps its logarithm finite;
        # the square-root weight remains tiny beside populated cells. This
        # prevents a high-absorption (E,Z) curve from disappearing entirely
        # and borrowing a grossly larger yield from a neighboring element.
        estimated_yield = (
            histories + 0.5) / (cell.events + 1.0)
        usable.append((cell.tau, estimated_yield, histories))
    if len(usable) < degree + 1:
        return None
    tau_lo = math.log1p(0.2)
    tau_hi = math.log1p(12.0)
    coordinate = np.asarray([
        2.0 * (math.log1p(tau) - tau_lo) / (tau_hi - tau_lo)
        - 1.0
        for tau, _, _ in usable])
    design = np.polynomial.chebyshev.chebvander(coordinate, degree)
    target = np.asarray([
        math.log(yield_value / tau ** 3)
        for tau, yield_value, _ in usable])
    weights = np.sqrt(np.minimum(
        np.maximum(
            np.asarray([
                histories for _, _, histories in usable]), 0.5),
        max(1000.0, float(minimum_histories))))
    return np.linalg.lstsq(
        design * weights[:, None],
        target * weights,
        rcond=None)[0]


def fit_tau_yield_model(
    cells: list[Cell],
    minimum_histories: int,
    feature_mode: str,
) -> tuple[TauYieldModel, dict]:
    campaign_cells = [cell for cell in cells if cell.split != "test"]
    groups: dict[tuple[float, int], list[Cell]] = {}
    for cell in campaign_cells:
        groups.setdefault(
            (cell.energy_keV, cell.z_value), []).append(cell)
    if feature_mode == "raw":
        bandwidths = [
            np.asarray([energy, z_value])
            for energy in (0.03, 0.05, 0.08, 0.12, 0.18)
            for z_value in (0.01, 0.015, 0.025, 0.04, 0.06)
        ]
    else:
        bandwidths = [
            np.asarray(
                [energy, z_value, process, process, process,
                 process, process])
            for energy in (0.03, 0.05, 0.08)
            for z_value in (0.015, 0.03, 0.06)
            for process in (0.03, 0.06, 0.10)]
    candidates = []
    for degree in (3, 4, 5):
        nodes = []
        for (energy, z_value), group_cells in groups.items():
            coefficients = fit_tau_yield_node(
                group_cells, degree, minimum_histories)
            if coefficients is None:
                continue
            identity = f"yield_E{energy:.12g}_Z{z_value}"
            node_split = split_for(identity, z_value, set())
            nodes.append((
                (
                    group_cells[0].feature[:2]
                    if feature_mode == "raw"
                    else yield_feature(
                        group_cells[0].feature, group_cells[0].tau)
                ),
                coefficients, node_split, group_cells))
        training = [node for node in nodes if node[2] == "train"]
        validation_cells = [
            cell
            for _, _, split, group_cells in nodes
            if split == "validation"
            for cell in group_cells]
        if not training or not validation_cells:
            continue
        centers = np.stack([node[0] for node in training])
        coefficients = np.stack([node[1] for node in training])
        for bandwidth in bandwidths:
            model = TauYieldModel(
                degree, bandwidth, centers, coefficients,
                use_process_features=feature_mode == "process")
            errors = []
            for cell in validation_cells:
                truth = sum(
                    cell.truth[name].sum()
                    for name in RAW_RESIDUAL_COMPONENTS)
                histories = truth * cell.events
                if histories < max(100.0, minimum_histories):
                    continue
                prediction = model.predict(cell.feature, cell.tau)
                errors.append(abs(prediction / truth - 1.0))
            if not errors:
                continue
            record = {
                "degree": degree,
                "bandwidth": bandwidth.tolist(),
                "nodes": len(training),
                "validation_cells": len(errors),
                "yield_median": float(np.median(errors)),
                "yield_p95": float(np.percentile(errors, 95)),
            }
            record["score"] = (
                record["yield_median"] + 0.25 * record["yield_p95"])
            candidates.append(record)
    if not candidates:
        raise ValueError("insufficient cells for tau-yield model")
    selected = min(candidates, key=lambda record: record["score"])

    # Refit the selected per-(E,Z) curves and 2-D coefficient RBF on all
    # non-test campaign groups.
    centers = []
    coefficients = []
    for (energy, z_value), group_cells in groups.items():
        value = fit_tau_yield_node(
            group_cells, selected["degree"], minimum_histories)
        if value is None:
            continue
        centers.append(
            group_cells[0].feature[:2]
            if feature_mode == "raw"
            else yield_feature(
                group_cells[0].feature, group_cells[0].tau))
        coefficients.append(value)
    model = TauYieldModel(
        selected["degree"], np.asarray(selected["bandwidth"]),
        np.stack(centers), np.stack(coefficients),
        use_process_features=feature_mode == "process")
    return model, {
        "refit_nodes": len(centers),
        "selected": selected,
        "candidates": candidates,
    }


def fit_component(
    cells: list[Cell],
    component: str,
    degree: int,
    minimum_histories: int,
    bandwidths: list[np.ndarray],
    refit_validation: bool,
) -> tuple[ComponentFit, dict]:
    usable = [
        (cell, per_cell_coefficients(
            cell, component, minimum_histories))
        for cell in cells]
    train = [
        (cell, coefficients)
        for cell, coefficients in usable
        if cell.split == "train" and coefficients is not None]
    validation = [
        cell for cell, coefficients in usable
        if cell.split == "validation" and coefficients is not None]
    if not train or not validation:
        raise ValueError(
            f"insufficient train/validation cells for {component}")
    centers = np.stack([cell.feature for cell, _ in train])
    log_yields = np.asarray([value[0] for _, value in train])
    shape_coefficients = np.stack([value[1] for _, value in train])
    candidates = []
    for bandwidth in bandwidths:
        fit = ComponentFit(
            component, degree, bandwidth, bandwidth, centers,
            log_yields, shape_coefficients)
        yield_errors = []
        shape_tv = []
        for cell in validation:
            prediction = cell.basis @ fit.predict_conditional(cell.feature)
            prediction *= sum(
                cell.model[name].sum() / FALLBACK_SCALES[name]
                for name in RAW_RESIDUAL_COMPONENTS)
            truth = sum(
                (cell.truth[name] for name in RAW_RESIDUAL_COMPONENTS),
                np.zeros(cell.basis.shape[0]))
            yield_errors.append(abs(prediction.sum() / truth.sum() - 1.0))
            # Select shape smoothing for the observable multiple-scatter
            # continuum, not for the normalized sparse 3+ residual alone.
            # The exact low-order truth is common to every candidate and
            # naturally downweights residual-shape fluctuations when that
            # residual is a small part of the emitted continuum.
            prediction *= truth.sum() / prediction.sum()
            low_order = sum(
                (cell.truth[name] for name in SOLVED_COMPONENTS),
                np.zeros(cell.basis.shape[0]))
            aggregate_prediction = low_order + prediction
            aggregate_truth = low_order + truth
            shape_tv.append(0.5 * np.abs(
                aggregate_prediction / aggregate_prediction.sum()
                - aggregate_truth / aggregate_truth.sum()).sum())
        candidates.append({
            "bandwidth": bandwidth.tolist(),
            "yield_median": float(np.median(yield_errors)),
            "yield_p95": float(np.percentile(yield_errors, 95)),
            "shape_tv_median": float(np.median(shape_tv)),
            "shape_tv_p95": float(np.percentile(shape_tv, 95)),
            "yield_score": float(
                np.median(yield_errors)
                + 0.25 * np.percentile(yield_errors, 95)),
            "shape_score": float(
                np.median(shape_tv)
                + 0.25 * np.percentile(shape_tv, 95)),
        })
    selected_yield = min(
        candidates, key=lambda record: record["yield_score"])
    selected_shape = min(
        candidates, key=lambda record: record["shape_score"])

    # Hyperparameters are selected without the validation cells. Refit the
    # final payload on every non-test campaign cell so the intentionally
    # sparse whole-cell validation split does not leave holes in the RBF bank.
    refit = [
        (cell, coefficients)
        for cell, coefficients in usable
        if (
            cell.split == "train"
            or (refit_validation and cell.split == "validation")
        ) and coefficients is not None]
    centers = np.stack([cell.feature for cell, _ in refit])
    log_yields = np.asarray([value[0] for _, value in refit])
    shape_coefficients = np.stack([value[1] for _, value in refit])
    fit = ComponentFit(
        component, degree,
        np.asarray(selected_yield["bandwidth"]),
        np.asarray(selected_shape["bandwidth"]),
        centers, log_yields, shape_coefficients)
    return fit, {
        "training_cells": len(train),
        "validation_cells": len(validation),
        "refit_cells": len(refit),
        "selected_yield": selected_yield,
        "selected_shape": selected_shape,
        "bandwidth_candidates": candidates,
    }


def evaluate_cell(
    cell: Cell,
    fits: dict[str, ComponentFit],
    yield_model: TauYieldModel,
    metric_bins: int,
) -> tuple[dict, np.ndarray]:
    truth = sum(cell.truth.values(), np.zeros(cell.basis.shape[0]))
    variance = sum(cell.variance.values(), np.zeros(cell.basis.shape[0]))
    prediction = sum(
        (cell.model[name] for name in SOLVED_COMPONENTS),
        np.zeros(cell.basis.shape[0]))
    component_metrics = {}
    for component, fit in fits.items():
        component_prediction = (
            cell.basis @ fit.predict_shape(cell.feature)
            * yield_model.predict(cell.feature, cell.tau))
        prediction += component_prediction
        component_truth = sum(
            (cell.truth[name] for name in RAW_RESIDUAL_COMPONENTS),
            np.zeros(cell.basis.shape[0]))
        if component_truth.sum() > 0.0:
            component_metrics[f"{component}_yield_error"] = float(
                component_prediction.sum() / component_truth.sum() - 1.0)
            component_metrics[f"{component}_shape_tv"] = float(
                0.5 * np.abs(
                    component_prediction / component_prediction.sum()
                    - component_truth / component_truth.sum()).sum())
    metadata = {
        "energy_keV": cell.energy_keV,
        "Z": float(cell.z_value),
        "tau": cell.tau,
        "areal_density_g_cm2": cell.areal_density,
        "events": float(cell.events),
    }
    metrics = calculate_metrics(
        truth, variance, prediction, metadata, metric_bins)
    metrics.update(component_metrics)
    metrics["cell_id"] = cell.identity
    metrics["split"] = cell.split
    return metrics, prediction


def resolved_metrics(
    cell: Cell,
    prediction: np.ndarray,
    channels: int,
    metric_bins: int,
    resolution: tuple[float, float, float],
) -> tuple[dict, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    if channels % metric_bins != 0:
        raise ValueError(
            "resolution metric bins must divide resolution channels")
    truth = sum(cell.truth.values(), np.zeros(cell.basis.shape[0]))
    variance = sum(cell.variance.values(), np.zeros(cell.basis.shape[0]))
    input_edges = np.linspace(
        0.0, cell.energy_keV, truth.size + 1)
    output_edges = np.linspace(
        0.0, cell.energy_keV, channels + 1)
    # Resolve on the native MC grid before the 4096-channel rebin. The
    # campaign stores only native-bin integrals, so subdividing first would
    # manufacture within-bin information and make the convolution needlessly
    # expensive.
    truth = resolve(input_edges, truth, *resolution)
    variance = resolve_variance(input_edges, variance, *resolution)
    prediction = resolve(input_edges, prediction, *resolution)
    truth = rebin(input_edges, truth, output_edges)
    variance = rebin_variance(input_edges, variance, output_edges)
    prediction = rebin(input_edges, prediction, output_edges)
    metadata = {
        "energy_keV": cell.energy_keV,
        "Z": float(cell.z_value),
        "tau": cell.tau,
        "areal_density_g_cm2": cell.areal_density,
        "events": float(cell.events),
    }
    metrics = calculate_metrics(
        truth, variance, prediction, metadata, metric_bins)
    return metrics, (output_edges, truth, variance, prediction)


def summarize(
    records: list[dict],
    split: str,
    prefix: str = "",
    eligibility_prefix: str | None = None,
) -> dict:
    def name(field: str) -> str:
        return prefix + field

    eligibility = (
        prefix if eligibility_prefix is None else eligibility_prefix)
    selected = [
        record for record in records
        if record["split"] == split
        and math.isfinite(record[name("weighted_rms_relative_error")])
        and math.isfinite(
            record[eligibility + "weighted_rms_relative_error"])
        and record[eligibility + "truth_yield"] > 0.0
        and record[eligibility + "truth_yield_sigma"]
            / record[eligibility + "truth_yield"] <= 0.10]
    if not selected:
        return {"cells": 0}
    return {
        "cells": len(selected),
        "yield_error_median_abs": float(np.median([
            abs(record[name("yield_relative_error")])
            for record in selected])),
        "yield_error_p95_abs": float(np.percentile([
            abs(record[name("yield_relative_error")])
            for record in selected], 95)),
        "shape_tv_median": float(np.median([
            record[name("shape_total_variation")]
            for record in selected])),
        "shape_tv_p95": float(np.percentile([
            record[name("shape_total_variation")]
            for record in selected], 95)),
        "weighted_rms_median": float(np.median([
            record[name("weighted_rms_relative_error")]
            for record in selected])),
        "weighted_rms_p95": float(np.percentile([
            record[name("weighted_rms_relative_error")]
            for record in selected], 95)),
        "bin_p95_error_median": float(np.median([
            record[name("p95_relative_error")]
            for record in selected])),
        "bin_p95_error_p95": float(np.percentile([
            record[name("p95_relative_error")]
            for record in selected], 95)),
    }


def payload_dict(fit: ComponentFit) -> dict:
    return {
        "degree": fit.degree,
        "yield_bandwidth": fit.yield_bandwidth.tolist(),
        "shape_bandwidth": fit.shape_bandwidth.tolist(),
        "centers": fit.centers.tolist(),
        "log_conditional_yields": fit.log_conditional_yields.tolist(),
        "shape_coefficients": fit.shape_coefficients.tolist(),
    }

def write_ceelo_payload(
    path: Path,
    fits: dict[str, ComponentFit],
    yield_model: TauYieldModel,
    absorption_envelope_strength: float,
    absorption_envelope_power: float,
    absorption_envelope_max_ratio: float,
    basis_type: str,
) -> None:
    names = {
        "higher_aggregate": "higher_aggregate",
    }
    lines = [
        "CEELO_MULTIPLE_SURROGATE 1",
        f"inverse_energy_warp_scale "
        f"{INVERSE_ENERGY_WARP_SCALE:.17g}",
        "absorption_envelope "
        f"{absorption_envelope_strength:.17g} "
        f"{absorption_envelope_power:.17g} "
        f"{absorption_envelope_max_ratio:.17g}",
        f"basis {basis_type}",
        "higher_mode merged",
        f"yield_model {yield_model.degree} "
        + " ".join(
            f"{value:.17g}" for value in yield_model.bandwidth)
        + f" {yield_model.minimum_tau:.17g} "
        f"{yield_model.maximum_tau:.17g} "
        f"{yield_model.centers.shape[0]}",
    ]
    for center, coefficients in zip(
            yield_model.centers, yield_model.coefficients):
        lines.append(
            "yield_node "
            + " ".join(
                f"{value:.17g}"
                for value in np.concatenate((center, coefficients))))
    for component in RESIDUAL_COMPONENTS:
        fit = fits[component]
        label = names[component]
        lines.append(
            f"{label} {fit.degree} "
            + " ".join(
                f"{value:.17g}" for value in fit.yield_bandwidth)
            + " "
            + " ".join(
                f"{value:.17g}" for value in fit.shape_bandwidth)
            + " "
            f"{fit.centers.shape[0]}")
        for center, log_yield, coefficients in zip(
                fit.centers, fit.log_conditional_yields,
                fit.shape_coefficients):
            values = np.concatenate((center, [log_yield], coefficients))
            lines.append(
                "node " + label + " "
                + " ".join(f"{value:.17g}" for value in values))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_plots(
    output: Path,
    cells: list[Cell],
    records: list[dict],
    predictions: dict[str, np.ndarray],
    resolved: dict[
        str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    worst: int,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    by_id = {cell.identity: cell for cell in cells}
    ranked = sorted(
        (record for record in records
         if record["split"] != "train"
         and math.isfinite(record["yield_relative_error"])
         and math.isfinite(record["shape_total_variation"])
         and record["truth_yield"] * record["events"] >= 100.0),
        key=lambda record: max(
            abs(record["yield_relative_error"]),
            record["shape_total_variation"]),
        reverse=True)
    with PdfPages(output / "multiple_scatter_surrogate_worst.pdf") as pdf:
        for record in ranked[:worst]:
            cell = by_id[record["cell_id"]]
            truth = sum(
                cell.truth.values(), np.zeros(cell.basis.shape[0]))
            truth_low_order = sum(
                (cell.truth[name] for name in SOLVED_COMPONENTS),
                np.zeros(cell.basis.shape[0]))
            truth_residual = sum(
                (cell.truth[name] for name in RAW_RESIDUAL_COMPONENTS),
                np.zeros(cell.basis.shape[0]))
            model_low_order = sum(
                (cell.model[name] for name in SOLVED_COMPONENTS),
                np.zeros(cell.basis.shape[0]))
            prediction = predictions[cell.identity]
            model_residual = np.maximum(
                0.0, prediction - model_low_order)
            variance = sum(
                cell.variance.values(),
                np.zeros(cell.basis.shape[0]))
            sigma = np.sqrt(np.maximum(variance, 0.0))
            edges = np.linspace(
                0.0, cell.energy_keV, truth.size + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            figure, (axis, ratio_axis) = plt.subplots(
                2, 1, figsize=(11, 7), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]},
                constrained_layout=True)
            axis.stairs(
                truth, edges, color="#202020", linewidth=1.3,
                label="CeeLo aggregate multiple scatter")
            axis.stairs(
                prediction, edges, color="#d1495b", linewidth=1.1,
                label="mechanistic low orders + functional residual")
            axis.stairs(
                truth_low_order, edges, color="#4c78a8",
                linewidth=0.8, alpha=0.85,
                label="CeeLo exact low-order subset")
            axis.stairs(
                truth_residual, edges, color="#59a14f",
                linewidth=0.8, alpha=0.85,
                label="CeeLo 3+ residual")
            axis.stairs(
                model_residual, edges, color="#f28e2b",
                linewidth=0.8, alpha=0.9,
                label="functional 3+ residual")
            axis.fill_between(
                centers, np.maximum(0.0, truth - sigma),
                truth + sigma, step="mid",
                color="#4c78a8", alpha=0.12,
                label="CeeLo 1σ")
            axis.set_yscale("log")
            axis.set_ylim(
                max(1e-14, truth.max() * 1e-6),
                max(truth.max(), prediction.max()) * 2)
            axis.set_ylabel("Yield / emitted photon / bin")
            axis.legend()
            axis.grid(alpha=0.2)
            ratio = np.full_like(truth, np.nan)
            valid = (
                truth > 0.0) & (
                sigma / np.maximum(truth, 1e-300) <= 0.20)
            ratio[valid] = prediction[valid] / truth[valid]
            ratio_axis.plot(centers, ratio, color="#d1495b")
            ratio_axis.axhline(1.0, color="#202020", linewidth=0.8)
            ratio_axis.set_ylim(0, 2)
            ratio_axis.set_xlabel("Escaped photon energy (keV)")
            ratio_axis.set_ylabel("model / truth")
            ratio_axis.grid(alpha=0.2)
            figure.suptitle(
                f"{record['split']}: {cell.energy_keV:g} keV, "
                f"Z={cell.z_value}, tau={cell.tau:g}; "
                f"yield {100*record['yield_relative_error']:+.1f}%, "
                f"TV {100*record['shape_total_variation']:.1f}%")
            pdf.savefig(figure)
            plt.close(figure)

    resolved_ranked = sorted(
        (record for record in records
         if record["split"] != "train"
         and math.isfinite(
             record["resolved_yield_relative_error"])
         and math.isfinite(
             record["resolved_shape_total_variation"])
         and record["resolved_truth_yield"]
             * record["events"] >= 100.0),
        key=lambda record: max(
            abs(record["resolved_yield_relative_error"]),
            record["resolved_shape_total_variation"]),
        reverse=True)
    with PdfPages(
            output
            / "multiple_scatter_surrogate_worst_resolved.pdf") as pdf:
        for record in resolved_ranked[:worst]:
            edges, truth, variance, prediction = resolved[
                record["cell_id"]]
            centers = 0.5 * (edges[:-1] + edges[1:])
            sigma = np.sqrt(np.maximum(variance, 0.0))
            figure, (axis, ratio_axis) = plt.subplots(
                2, 1, figsize=(11, 7), sharex=True,
                gridspec_kw={"height_ratios": [2.2, 1.0]},
                constrained_layout=True)
            axis.stairs(
                truth, edges, color="#202020", linewidth=1.3,
                label="CeeLo after resolution")
            axis.stairs(
                prediction, edges, color="#d1495b", linewidth=1.1,
                label="functional model after resolution")
            axis.fill_between(
                centers, np.maximum(0.0, truth - sigma),
                truth + sigma, step="mid",
                color="#4c78a8", alpha=0.12,
                label="propagated CeeLo 1σ")
            axis.set_yscale("log")
            axis.set_ylim(
                max(1e-14, truth.max() * 1e-6),
                max(truth.max(), prediction.max()) * 2)
            axis.set_ylabel("Resolved yield / channel")
            axis.legend()
            axis.grid(alpha=0.2)
            ratio = np.full_like(truth, np.nan)
            valid = (
                truth > max(1e-14, truth.max() * 1e-8)
            ) & (
                sigma / np.maximum(truth, 1e-300) <= 0.20)
            ratio[valid] = prediction[valid] / truth[valid]
            ratio_axis.plot(centers, ratio, color="#d1495b")
            ratio_axis.axhline(1.0, color="#202020", linewidth=0.8)
            ratio_axis.set_ylim(0, 2)
            ratio_axis.set_xlabel("Escaped photon energy (keV)")
            ratio_axis.set_ylabel("model / truth")
            ratio_axis.grid(alpha=0.2)
            figure.suptitle(
                f"{record['split']}: {record['energy_keV']:g} keV, "
                f"Z={record['Z']}, tau={record['tau']:g}; "
                "resolved yield "
                f"{100*record['resolved_yield_relative_error']:+.1f}%, "
                "TV "
                f"{100*record['resolved_shape_total_variation']:.1f}%")
            pdf.savefig(figure)
            plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path)
    parser.add_argument(
        "--test-campaign", type=Path,
        help="disjoint whole-cell campaign used only for final testing")
    parser.add_argument(
        "--validation-campaign", type=Path,
        help="structured whole-cell campaign used only for model selection")
    parser.add_argument(
        "--additional-training-campaign", type=Path, action="append",
        default=[],
        help="adaptive whole-cell campaign promoted to training after a "
             "completed diagnostic round; may be repeated")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--degree", type=int, default=31)
    parser.add_argument("--metric-bins", type=int, default=32)
    parser.add_argument("--minimum-histories", type=int, default=200)
    parser.add_argument(
        "--minimum-yield-histories", type=int, default=40,
        help="minimum aggregate histories at a tau node for yield curves")
    parser.add_argument(
        "--resolution", type=float, nargs=3,
        metavar=("A0", "A1", "A2"),
        default=(4.0, 1.5, 2.0e-4),
        help="observable-spectrum FWHM^2=a0+a1*E+a2*E^2")
    parser.add_argument(
        "--resolution-channels", type=int, default=4096)
    parser.add_argument(
        "--resolution-metric-bins", type=int, default=64)
    parser.add_argument(
        "--cross-section-helper", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_cross_section_features"))
    parser.add_argument(
        "--absorption-envelope-strength", type=float,
        default=ABSORPTION_ENVELOPE_STRENGTH)
    parser.add_argument(
        "--absorption-envelope-power", type=float,
        default=ABSORPTION_ENVELOPE_POWER)
    parser.add_argument(
        "--absorption-envelope-max-ratio", type=float,
        default=ABSORPTION_ENVELOPE_MAX_RATIO)
    parser.add_argument(
        "--basis", choices=("cubic_bspline", "bernstein"),
        default="cubic_bspline")
    parser.add_argument(
        "--feature-mode", choices=("raw", "process"), default="raw",
        help="RBF coordinates for continuum-shape interpolation")
    parser.add_argument(
        "--yield-feature-mode", choices=("raw", "process"), default="raw",
        help="coordinates used to interpolate thickness-yield curves; "
             "process is retained as an ablation")
    parser.add_argument("--holdout-z", type=int, nargs="*", default=[])
    parser.add_argument("--worst-plots", type=int, default=20)
    args = parser.parse_args()
    output = args.output or args.campaign / "multiple_surrogate"
    output.mkdir(parents=True, exist_ok=True)

    primary_cells = load_cells(
        args.campaign, args.degree, set(args.holdout_z),
        args.cross_section_helper,
        args.absorption_envelope_strength,
        args.absorption_envelope_power,
        args.absorption_envelope_max_ratio,
        args.basis,
        args.feature_mode)
    cells = list(primary_cells)
    for training_campaign in args.additional_training_campaign:
        adaptive_cells = load_cells(
            training_campaign, args.degree, set(),
            args.cross_section_helper,
            args.absorption_envelope_strength,
            args.absorption_envelope_power,
            args.absorption_envelope_max_ratio,
            args.basis,
            args.feature_mode)
        for cell in adaptive_cells:
            cell.split = "train"
        primary_cells.extend(adaptive_cells)
        cells.extend(adaptive_cells)
    if args.validation_campaign is not None:
        for cell in primary_cells:
            cell.split = "train"
        validation_cells = load_cells(
            args.validation_campaign, args.degree, set(),
            args.cross_section_helper,
            args.absorption_envelope_strength,
            args.absorption_envelope_power,
            args.absorption_envelope_max_ratio,
            args.basis,
            args.feature_mode)
        for cell in validation_cells:
            cell.split = "validation"
        cells.extend(validation_cells)
    if args.test_campaign is not None:
        test_cells = load_cells(
            args.test_campaign, args.degree, set(),
            args.cross_section_helper,
            args.absorption_envelope_strength,
            args.absorption_envelope_power,
            args.absorption_envelope_max_ratio,
            args.basis,
            args.feature_mode)
        for cell in test_cells:
            cell.split = "test"
        training_identities = {
            cell.identity for cell in primary_cells}
        cells.extend(
            cell for cell in test_cells
            if cell.identity not in training_identities)
    if args.feature_mode == "raw":
        bandwidths = [
            np.array([energy, z_value, tau])
            for energy in (0.03, 0.05, 0.08, 0.12)
            for z_value in (0.015, 0.025, 0.04, 0.06)
            for tau in (0.015, 0.025, 0.04, 0.06)
        ]
    else:
        bandwidths = [
            np.asarray(
                [energy, z_value, tau, scatter, scatter,
                 absorption, absorption, 0.10])
            for energy in (0.03, 0.05, 0.08)
            for z_value in (0.015, 0.03)
            for tau in (0.025, 0.05, 0.08)
            for scatter in (0.06, 0.12)
            for absorption in (0.04, 0.08, 0.15)
        ]
    fits = {}
    fit_reports = {}
    yield_model, yield_report = fit_tau_yield_model(
        primary_cells, args.minimum_yield_histories,
        args.yield_feature_mode)
    for component in RESIDUAL_COMPONENTS:
        fits[component], fit_reports[component] = fit_component(
            cells, component, args.degree,
            args.minimum_histories, bandwidths,
            args.validation_campaign is None)

    records = []
    predictions = {}
    resolved = {}
    for cell in cells:
        record, prediction = evaluate_cell(
            cell, fits, yield_model, args.metric_bins)
        resolved_record, resolved_arrays = resolved_metrics(
            cell, prediction,
            args.resolution_channels,
            args.resolution_metric_bins,
            tuple(args.resolution))
        for name, value in resolved_record.items():
            if name not in {
                "energy_keV", "Z", "tau",
                "areal_density_g_cm2", "events"
            }:
                record[f"resolved_{name}"] = value
        records.append(record)
        predictions[cell.identity] = prediction
        resolved[cell.identity] = resolved_arrays
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    payload = {
        "schema": 1,
        "coordinate": (
            "w=q/(s+q), q=(m_e/2)*(1/E-1/E0)"),
        "basis": args.basis,
        "inverse_energy_warp_scale": INVERSE_ENERGY_WARP_SCALE,
        "absorption_envelope": {
            "strength": args.absorption_envelope_strength,
            "power": args.absorption_envelope_power,
            "maximum_energy_ratio":
                args.absorption_envelope_max_ratio,
            "form": (
                "exp(-strength*areal_density*mu_pe(E0)"
                "*(min(maximum_energy_ratio,E0/E)^power-1))"),
        },
        "higher_mode": "merged",
        "feature_mode": args.feature_mode,
        "features": (
            [
                "log(E0/10)/log(1200)", "Z/92",
                "log(1+tau)/log(25)",
            ]
            if args.feature_mode == "raw"
            else [
                "log(E0/10)/log(1200)", "Z/92",
                "log(1+tau_total)/log(25)",
                "log(1+tau_cs_incident)/log(25)",
                "log(1+tau_rs_incident)/log(25)",
                "log(1+tau_abs_at_reference3)/log(25)",
                "log(1+tau_abs_at_reference8)/log(25)",
                "log(1+tau_pp_incident)/log(25)",
            ]),
        "yield_model": {
            "form": (
                "tau^3*exp(Chebyshev(log1p(tau))), "
                "2-D RBF coefficients over log(E),Z"),
            "degree": yield_model.degree,
            "bandwidth": yield_model.bandwidth.tolist(),
            "minimum_tau": yield_model.minimum_tau,
            "maximum_tau": yield_model.maximum_tau,
            "centers": yield_model.centers.tolist(),
            "coefficients": yield_model.coefficients.tolist(),
        },
        "components": {
            name: payload_dict(fit) for name, fit in fits.items()},
    }
    (output / "payload.json").write_text(
        json.dumps(payload, separators=(",", ":")) + "\n")
    write_ceelo_payload(
        output / "surrogate.txt", fits, yield_model,
        args.absorption_envelope_strength,
        args.absorption_envelope_power,
        args.absorption_envelope_max_ratio,
        args.basis)
    report = {
        "schema": 1,
        "campaign": str(args.campaign),
        "additional_training_campaigns": [
            str(path) for path in args.additional_training_campaign],
        "validation_campaign": (
            str(args.validation_campaign)
            if args.validation_campaign is not None else None),
        "test_campaign": (
            str(args.test_campaign)
            if args.test_campaign is not None else None),
        "cells": len(cells),
        "degree": args.degree,
        "minimum_component_histories": args.minimum_histories,
        "minimum_yield_histories": args.minimum_yield_histories,
        "absorption_envelope_strength":
            args.absorption_envelope_strength,
        "absorption_envelope_power":
            args.absorption_envelope_power,
        "absorption_envelope_max_ratio":
            args.absorption_envelope_max_ratio,
        "basis": args.basis,
        "feature_mode": args.feature_mode,
        "yield_feature_mode": args.yield_feature_mode,
        "yield_fit": yield_report,
        "fit": fit_reports,
        "training": summarize(records, "train"),
        "validation": summarize(records, "validation"),
        "test": summarize(records, "test"),
        "resolved": {
            "channels": args.resolution_channels,
            "metric_bins": args.resolution_metric_bins,
            "fwhm2_parameters": list(args.resolution),
            "training": summarize(
                records, "train", "resolved_", ""),
            "validation": summarize(
                records, "validation", "resolved_", ""),
            "test": summarize(
                records, "test", "resolved_", ""),
        },
        "payload_bytes": (
            output / "payload.json").stat().st_size,
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    write_plots(
        output, cells, records, predictions, resolved,
        args.worst_plots)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
