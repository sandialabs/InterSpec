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

"""Fit the charged-particle source payload.

Per component:
  - a tensor-product clamped-uniform cubic B-spline yield field in log space
    on the normalized features (freeze-clamped domain), and
  - a rank-R softmax mixture of nonnegative unit-sum 32-vector shape
    archetypes whose R-1 log-odds are tensor B-spline fields.

The representation, the payload text format (`CEELO_CHARGED_PARTICLE_SOURCE
1`) and the evaluation semantics are documented on
ChargedParticleSplineField / ChargedParticleSourcePayload in
src/full_spectrum/FullSpectrumModel.h and must stay in exact lockstep with the
C++ evaluator. This script is the Python parity reference: `load_payload` +
`payload_predict_fn` evaluate a serialized payload exactly as C++ does.

Cell loading, features, mechanistic bases, per-cell NNLS shape targets and
metrics come from fit_charged_particle_source.py.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")

import argparse
import bank_io
import hashlib
import json
import math
import pickle
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_charged_particle_source as fit  # noqa: E402

BASIS_SIZE = fit.BASIS_SIZE
MODEL_COMPONENTS = fit.MODEL_COMPONENTS
YIELD_ONLY_COMPONENTS = {"positron_containment", "positron_511_atom",
                         *fit.FLUORESCENCE_COMPONENTS}
DEFAULT_PAPF = 2.0  # validation-selected photon absorption path fraction
LOG_ODDS_EPS = 0.01
WEIGHT_FIELD_COUNT_CAP = 1.0e4
RIDGE_RELATIVE = 1.0e-6


# ---------------------------------------------------------------------------
# Generalized clamped-uniform cubic B-spline basis (spec section "Basis")
# ---------------------------------------------------------------------------

def bspline_knots(n: int) -> np.ndarray:
    degree = 3
    knot_count = n + degree + 1
    knots = np.empty(knot_count)
    for index in range(knot_count):
        if index <= degree:
            knots[index] = 0.0
        elif index >= n:
            knots[index] = 1.0
        else:
            knots[index] = (index - degree) / (n - degree)
    return knots


def bspline_basis_matrix(u: np.ndarray, n: int) -> np.ndarray:
    """Rows of cubic B-spline basis values; identical knot convention to
    fit.cubic_bspline_values (verified for n == 32 in --selftest)."""
    u = np.atleast_1d(np.asarray(u, dtype=float))
    degree = 3
    knots = bspline_knots(n)
    work = n + degree
    previous = np.zeros((u.size, work))
    for index in range(knots.size - 1):
        previous[:, index] = (knots[index] <= u) & (u < knots[index + 1])
    for order in range(1, degree + 1):
        current = np.zeros((u.size, work))
        for index in range(work - order):
            left = knots[index + order] - knots[index]
            if left > 0.0:
                current[:, index] += (
                    previous[:, index] * (u - knots[index]) / left)
            right = knots[index + order + 1] - knots[index + 1]
            if right > 0.0:
                current[:, index] += (
                    previous[:, index + 1]
                    * (knots[index + order + 1] - u) / right)
        previous = current
    answer = previous[:, :n]
    at_one = u >= 1.0
    answer[at_one] = 0.0
    answer[at_one, n - 1] = 1.0
    return answer


def clamp01(u: np.ndarray) -> np.ndarray:
    return np.clip(u, 0.0, 1.0)


def normalized_axes(
    features: np.ndarray, domain: np.ndarray,
) -> list[np.ndarray]:
    axes = []
    for axis in range(domain.shape[0]):
        lo, hi = domain[axis]
        axes.append(clamp01((features[:, axis] - lo) / (hi - lo)))
    return axes


def tensor_design(
    features: np.ndarray,
    domain: np.ndarray,
    grid: tuple[int, ...],
) -> scipy.sparse.csr_matrix:
    """Sparse design matrix; each row is the outer product of per-axis basis
    rows (<= 4 nonzeros per axis for cubic splines)."""
    axes_u = normalized_axes(features, domain)
    count = features.shape[0]
    windows = []
    values = []
    for u, n in zip(axes_u, grid):
        dense = bspline_basis_matrix(u, n)
        start = np.clip(
            np.floor(u * (n - 3)).astype(int), 0, n - 4)
        cols = start[:, None] + np.arange(4)[None, :]
        vals = np.take_along_axis(dense, cols, axis=1)
        windows.append(cols)
        values.append(vals)
    n0, n1, n2 = grid
    data = (
        values[0][:, :, None, None]
        * values[1][:, None, :, None]
        * values[2][:, None, None, :]).reshape(count, 64)
    cols = (
        (windows[0][:, :, None, None] * n1
         + windows[1][:, None, :, None]) * n2
        + windows[2][:, None, None, :]).reshape(count, 64)
    rows = np.repeat(np.arange(count), 64)
    return scipy.sparse.csr_matrix(
        (data.ravel(), (rows, cols.ravel())),
        shape=(count, n0 * n1 * n2))


def second_difference_penalty(
    grid: tuple[int, ...], lambdas: tuple[float, ...],
) -> scipy.sparse.csr_matrix:
    """sum_axis lambda_axis * (n_axis - 3)^3 * ||D2_axis c||^2 as a sparse
    quadratic form (the (n-3)^3 factor approximates an integrated squared
    second derivative in u, keeping lambda comparable across grid sizes)."""
    total = None
    eyes = [scipy.sparse.identity(n, format="csr") for n in grid]
    for axis, (n, lam) in enumerate(zip(grid, lambdas)):
        if n < 3 or lam <= 0.0:
            continue
        diff = scipy.sparse.diags(
            [np.ones(n - 2), -2.0 * np.ones(n - 2), np.ones(n - 2)],
            offsets=[0, 1, 2], shape=(n - 2, n), format="csr")
        quad = (diff.T @ diff) * (lam * (n - 3) ** 3)
        factors = [
            quad if index == axis else eyes[index]
            for index in range(len(grid))]
        term = scipy.sparse.kron(
            scipy.sparse.kron(factors[0], factors[1]), factors[2],
            format="csr")
        total = term if total is None else total + term
    if total is None:
        size = int(np.prod(grid))
        total = scipy.sparse.csr_matrix((size, size))
    return total


def solve_penalized(
    design: scipy.sparse.csr_matrix,
    target: np.ndarray,
    weight: np.ndarray,
    penalty: scipy.sparse.csr_matrix,
) -> np.ndarray:
    weighted = design.T.multiply(weight)
    gram = (weighted @ design).tocsr()
    rhs = weighted @ target
    scale = gram.diagonal().sum() / max(gram.shape[0], 1)
    ridge = RIDGE_RELATIVE * (scale if scale > 0.0 else 1.0)
    system = (
        gram + penalty
        + ridge * scipy.sparse.identity(gram.shape[0], format="csr"))
    answer = scipy.sparse.linalg.spsolve(system.tocsc(), rhs)
    if not np.all(np.isfinite(answer)):
        raise RuntimeError("penalized solve produced non-finite output")
    return np.asarray(answer)


def tensor_evaluate(
    features: np.ndarray,
    domain: np.ndarray,
    grid: tuple[int, ...],
    coefficients: np.ndarray,
) -> np.ndarray:
    design = tensor_design(features, domain, grid)
    return design @ coefficients


# ---------------------------------------------------------------------------
# Weighted NMF shape archetypes
# ---------------------------------------------------------------------------

def weighted_nmf(
    shapes: np.ndarray,
    sample_weight: np.ndarray,
    rank: int,
    iterations: int = 400,
    seed: int = 20260731,
) -> tuple[np.ndarray, np.ndarray, float]:
    rng = np.random.default_rng(seed)
    count = shapes.shape[0]
    mixing = rng.uniform(0.5, 1.5, size=(count, rank)) / rank
    archetypes = rng.uniform(0.5, 1.5, size=(rank, BASIS_SIZE))
    archetypes *= shapes.mean(axis=0)[None, :] / max(
        archetypes.mean(), 1e-300)
    s = sample_weight[:, None]
    sx = s * shapes
    floor = 1e-12
    for _ in range(iterations):
        model = mixing @ archetypes
        archetypes *= (mixing.T @ sx) / np.maximum(
            mixing.T @ (s * model), floor)
        archetypes = np.maximum(archetypes, floor)
        model = mixing @ archetypes
        mixing *= (sx @ archetypes.T) / np.maximum(
            (s * model) @ archetypes.T, floor)
        mixing = np.maximum(mixing, floor)
    if not (np.all(np.isfinite(archetypes))
            and np.all(np.isfinite(mixing))):
        raise RuntimeError("weighted NMF diverged")
    row_sums = archetypes.sum(axis=1)
    row_sums[row_sums <= 0.0] = 1.0
    archetypes /= row_sums[:, None]
    mixing *= row_sums[None, :]
    residual = shapes - mixing @ archetypes
    error = float(
        np.sum(sample_weight * np.sum(residual * residual, axis=1))
        / max(np.sum(sample_weight * np.sum(shapes * shapes, axis=1)),
              1e-300))
    return mixing, archetypes, error


def robust_nnls(matrix: np.ndarray, target: np.ndarray) -> np.ndarray:
    """scipy NNLS, falling back to fit_charged_particle_source's
    Lawson-Hanson when the LAPACK path returns non-finite output
    (ill-conditioned attenuated bases)."""
    try:
        solution, _ = scipy.optimize.nnls(matrix, target)
        if np.all(np.isfinite(solution)):
            return solution
    except Exception:
        pass
    solution, _ = fit.nonnegative_least_squares(matrix, target)
    if not np.all(np.isfinite(solution)):
        raise RuntimeError("NNLS produced non-finite output")
    return solution


def mixture_weights(
    shapes: np.ndarray, archetypes: np.ndarray,
) -> np.ndarray:
    rank = archetypes.shape[0]
    basis = archetypes.T
    answer = np.empty((shapes.shape[0], rank))
    for index in range(shapes.shape[0]):
        solution = robust_nnls(basis, shapes[index])
        total = solution.sum()
        answer[index] = (
            solution / total if total > 0.0
            else np.full(rank, 1.0 / rank))
    return answer


# ---------------------------------------------------------------------------
# Model containers
# ---------------------------------------------------------------------------

@dataclass
class StructuredComponent:
    name: str
    domain: np.ndarray | None = None            # (3, 2)
    yield_grid: tuple[int, int, int] | None = None
    yield_coefficients: np.ndarray | None = None
    rank: int = 0
    archetypes: np.ndarray | None = None        # (R, 32)
    weight_grid: tuple[int, int, int] | None = None
    log_odds: np.ndarray | None = None          # (R-1, m0*m1*m2)
    # bookkeeping
    lambdas: tuple[float, float, float] | None = None
    weight_lambdas: tuple[float, float, float] | None = None

    def predict_yield(self, features: np.ndarray) -> np.ndarray:
        return np.exp(tensor_evaluate(
            features, self.domain, self.yield_grid,
            self.yield_coefficients))

    def predict_shape(self, features: np.ndarray) -> np.ndarray:
        count = features.shape[0]
        if self.rank == 0:
            return np.full(
                (count, BASIS_SIZE), 1.0 / BASIS_SIZE)
        if self.rank == 1:
            shape = self.archetypes[0]
            shape = shape / max(shape.sum(), 1e-300)
            return np.tile(shape, (count, 1))
        logits = np.zeros((count, self.rank))
        for index in range(self.rank - 1):
            logits[:, index] = tensor_evaluate(
                features, self.domain, self.weight_grid,
                self.log_odds[index])
        # max-shift softmax (spec)
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        weights /= weights.sum(axis=1, keepdims=True)
        shapes = weights @ self.archetypes
        shapes /= np.maximum(
            shapes.sum(axis=1, keepdims=True), 1e-300)
        return shapes


@dataclass
class StructuredPayload:
    photon_absorption_path_fraction: float
    components: dict[str, StructuredComponent] = field(
        default_factory=dict)


def payload_predict_fn(payload: StructuredPayload):
    """The predictor fit.evaluate calls, evaluating `payload`."""

    def predict(
        cell, name, containment_probability=None, basis_cache=None,
        photon_absorption_path_fraction=0.0,
    ):
        component = payload.components.get(name)
        output_size = (
            1 if name in YIELD_ONLY_COMPONENTS
            else cell.edges.size - 1)
        if component is None or not fit.component_applies(cell, name):
            return np.zeros(output_size)
        feature = fit.cell_feature(cell, name)[None, :]
        conditional = float(component.predict_yield(feature)[0])
        total = fit.component_base(
            cell, name, containment_probability) * conditional
        if name == "positron_containment":
            return np.asarray([min(1.0, total)])
        if name in fit.ATOM_COMPONENTS:
            return np.asarray([total])
        shape = component.predict_shape(feature)[0]
        basis_key = (
            name, cell.energy_keV, cell.z_value,
            cell.areal_density, cell.edges.size,
            float(cell.edges[-1]),
            photon_absorption_path_fraction)
        basis = (
            basis_cache.get(basis_key)
            if basis_cache is not None else None)
        if basis is None:
            basis = fit.component_basis(
                cell, name, photon_absorption_path_fraction)
            if basis_cache is not None:
                basis_cache[basis_key] = basis
        return total * (basis @ shape)

    return predict


# ---------------------------------------------------------------------------
# Payload serialization (CEELO_CHARGED_PARTICLE_SOURCE 1)
# ---------------------------------------------------------------------------

def format_values(values) -> str:
    return " ".join(f"{float(value):.17g}" for value in values)


def write_payload(path: Path, payload: StructuredPayload) -> None:
    with path.open("w") as stream:
        stream.write("CEELO_CHARGED_PARTICLE_SOURCE 1\n")
        stream.write("validity 10 12000 10 1 92 0 60\n")
        stream.write(
            f"features 0.41999999999999998 "
            f"{fit.RANGE_FEATURE_SCALE:.17g} "
            f"{payload.photon_absorption_path_fraction:.17g}\n")
        for name in MODEL_COMPONENTS:
            component = payload.components.get(name)
            if component is None:
                stream.write(f"{name} none\n")
                continue
            stream.write(f"{name} spline\n")
            stream.write(
                f"{name}_domain 3 "
                + format_values(component.domain.ravel()) + "\n")
            n0, n1, n2 = component.yield_grid
            stream.write(f"{name}_yield_grid 3 {n0} {n1} {n2}\n")
            stream.write(
                f"{name}_yield_coefficients "
                + format_values(component.yield_coefficients) + "\n")
            stream.write(f"{name}_shape_rank {component.rank}\n")
            for row in range(component.rank):
                stream.write(
                    f"{name}_shape_basis "
                    + format_values(component.archetypes[row]) + "\n")
            if component.rank >= 2:
                m0, m1, m2 = component.weight_grid
                stream.write(
                    f"{name}_weight_grid 3 {m0} {m1} {m2}\n")
                for row in range(component.rank - 1):
                    stream.write(
                        f"{name}_weight_log_odds "
                        + format_values(component.log_odds[row])
                        + "\n")
        stream.write("END\n")


def load_payload(path: Path) -> StructuredPayload:
    tokens_by_line = [
        line.split() for line in path.read_text().splitlines()
        if line.strip()]
    cursor = 0

    def next_line():
        nonlocal cursor
        line = tokens_by_line[cursor]
        cursor += 1
        return line

    header = next_line()
    if header != ["CEELO_CHARGED_PARTICLE_SOURCE", "1"]:
        raise ValueError(
            "not a version-1 charged-particle source payload")
    next_line()  # validity
    features_line = next_line()
    if len(features_line) != 4 or features_line[0] != "features":
        raise ValueError(f"bad features record: {features_line}")
    payload = StructuredPayload(
        photon_absorption_path_fraction=float(features_line[3]))
    for name in MODEL_COMPONENTS:
        tag = next_line()
        if tag[0] != name or tag[1] not in {"spline", "none"}:
            raise ValueError(f"bad component tag: {tag}")
        if tag[1] == "none":
            continue
        domain_line = next_line()
        assert domain_line[0] == f"{name}_domain"
        naxes = int(domain_line[1])
        domain = np.asarray(
            [float(v) for v in domain_line[2:2 + 2 * naxes]]
        ).reshape(naxes, 2)
        grid_line = next_line()
        assert grid_line[0] == f"{name}_yield_grid"
        grid = tuple(int(v) for v in grid_line[2:2 + naxes])
        coeff_line = next_line()
        assert coeff_line[0] == f"{name}_yield_coefficients"
        coefficients = np.asarray([float(v) for v in coeff_line[1:]])
        if coefficients.size != int(np.prod(grid)):
            raise ValueError(f"{name}: coefficient count mismatch")
        rank_line = next_line()
        assert rank_line[0] == f"{name}_shape_rank"
        rank = int(rank_line[1])
        archetypes = None
        weight_grid = None
        log_odds = None
        if rank > 0:
            rows = []
            for _ in range(rank):
                row = next_line()
                assert row[0] == f"{name}_shape_basis"
                rows.append([float(v) for v in row[1:]])
            archetypes = np.asarray(rows)
            if archetypes.shape != (rank, BASIS_SIZE):
                raise ValueError(f"{name}: archetype shape mismatch")
            if np.any(archetypes < 0.0):
                raise ValueError(f"{name}: negative archetype")
        if rank >= 2:
            wg_line = next_line()
            assert wg_line[0] == f"{name}_weight_grid"
            weight_grid = tuple(int(v) for v in wg_line[2:2 + naxes])
            rows = []
            for _ in range(rank - 1):
                row = next_line()
                assert row[0] == f"{name}_weight_log_odds"
                rows.append([float(v) for v in row[1:]])
            log_odds = np.asarray(rows)
            if log_odds.shape[1] != int(np.prod(weight_grid)):
                raise ValueError(f"{name}: log-odds count mismatch")
        payload.components[name] = StructuredComponent(
            name=name, domain=domain,
            yield_grid=grid, yield_coefficients=coefficients,
            rank=rank, archetypes=archetypes,
            weight_grid=weight_grid, log_odds=log_odds)
    end_line = next_line()
    if end_line[0] != "END":
        raise ValueError("missing END")
    return payload


# ---------------------------------------------------------------------------
# Training-data extraction (mirrors fit.fit_nodes qualifying rules)
# ---------------------------------------------------------------------------

@dataclass
class ComponentData:
    identities: list[str]
    features: np.ndarray       # (n, 3)
    log_yield: np.ndarray      # (n,)
    sigma_rel: np.ndarray      # (n,)
    counts: np.ndarray         # (n,)
    splits: np.ndarray         # (n,) object array of split strings
    shapes: np.ndarray | None  # (n, 32) or None for yield-only
    base: np.ndarray | None = None       # mechanistic base per cell
    raw_total: np.ndarray | None = None  # unsmoothed truth total


def component_training_data(
    cells: list[fit.Cell],
    name: str,
    basis_cache: dict,
    papf: float,
    minimum_counts: float = 25.0,
    shape_cache: dict | None = None,
) -> ComponentData:
    identities = []
    features = []
    log_yields = []
    sigma_rels = []
    counts_list = []
    splits = []
    bases = []
    raw_totals = []
    shapes = [] if name not in YIELD_ONLY_COMPONENTS else None
    for cell in cells:
        if not fit.component_applies(cell, name):
            continue
        truth, variance = fit.truth_for(cell, name)
        total = float(truth.sum())
        raw_total = total
        if name == "positron_containment":
            total = (total * cell.events + 0.5) / (cell.events + 1.0)
        base = fit.component_base(cell, name)
        if (
            name != "positron_containment"
            and raw_total * cell.events < minimum_counts
        ) or base <= 0.0:
            continue
        if shapes is not None:
            key = (cell.identity, name, papf)
            coefficients = (
                shape_cache.get(key) if shape_cache is not None
                else None)
            if coefficients is None:
                basis_key = (
                    name, cell.energy_keV, cell.z_value,
                    cell.areal_density, cell.edges.size,
                    float(cell.edges[-1]), papf)
                basis = basis_cache.get(basis_key)
                if basis is None:
                    basis = fit.component_basis(cell, name, papf)
                    basis_cache[basis_key] = basis
                coefficients = robust_nnls(basis, truth / total)
                if shape_cache is not None:
                    shape_cache[key] = coefficients
            coefficient_sum = float(np.sum(coefficients))
            if coefficient_sum <= 0.0:
                continue
            shapes.append(
                np.asarray(coefficients) / coefficient_sum)
        variance_sum = float(np.sum(variance))
        sigma_rel = (
            math.sqrt(max(variance_sum, 0.0)) / total
            if total > 0.0 else math.inf)
        if name == "positron_containment":
            # The stored escaped-primary variance is the uncentered sum_w2
            # estimator; for the complement 1 - escaped (escaped ~ 1) it
            # overstates the noise by orders of magnitude (missing -p^2/N).
            # Use the binomial sigma of the smoothed contained fraction.
            sigma_rel = math.sqrt(
                max(total * (1.0 - total), 0.0) / cell.events
                + 1.0 / (cell.events * cell.events)) / total
        identities.append(cell.identity)
        features.append(fit.cell_feature(cell, name))
        log_yields.append(math.log(max(total / base, 1e-300)))
        sigma_rels.append(sigma_rel)
        counts_list.append(raw_total * cell.events)
        splits.append(cell.split)
        bases.append(base)
        raw_totals.append(raw_total)
    return ComponentData(
        identities=identities,
        features=np.asarray(features).reshape(-1, 3),
        log_yield=np.asarray(log_yields),
        sigma_rel=np.asarray(sigma_rels),
        counts=np.asarray(counts_list),
        splits=np.asarray(splits, dtype=object),
        shapes=(
            np.asarray(shapes).reshape(-1, BASIS_SIZE)
            if shapes is not None else None),
        base=np.asarray(bases),
        raw_total=np.asarray(raw_totals))


def training_domain(features: np.ndarray) -> np.ndarray:
    lo = features.min(axis=0)
    hi = features.max(axis=0)
    span = np.maximum(hi - lo, 1e-6)
    pad = 0.01 * span
    return np.stack([lo - pad, hi + pad], axis=1)


def component_domain(
    name: str, data: ComponentData, mask: np.ndarray,
) -> np.ndarray:
    """Fit domain per axis. For positron_containment the f2 (gate-depth)
    axis is truncated at the saturation point: beyond it every training
    cell is fully contained (truth ~ 1, log-yield ~ 0), so the freeze
    clamp makes the field constant there and all knots concentrate in the
    sharp turn-on region instead of the long saturated tail."""
    domain = training_domain(data.features[mask])
    if name == "positron_containment":
        unsaturated = mask & (data.raw_total < 0.9995)
        if np.any(unsaturated) and np.any(mask & ~unsaturated):
            span = domain[2, 1] - domain[2, 0]
            hi = float(data.features[unsaturated, 2].max())
            domain[2, 1] = min(
                domain[2, 1], hi + 0.02 * span)
    return domain


def yield_weights(data: ComponentData, mask: np.ndarray) -> np.ndarray:
    sigma = np.maximum(data.sigma_rel[mask], 0.01)
    weight = 1.0 / (sigma * sigma)
    return weight / weight.mean()


def fit_yield_field(
    data: ComponentData,
    mask: np.ndarray,
    domain: np.ndarray,
    grid: tuple[int, int, int],
    lam: float,
) -> np.ndarray:
    design = tensor_design(data.features[mask], domain, grid)
    penalty = second_difference_penalty(grid, (lam, lam, lam))
    return solve_penalized(
        design, data.log_yield[mask], yield_weights(data, mask),
        penalty)


def fit_shape_model(
    data: ComponentData,
    mask: np.ndarray,
    domain: np.ndarray,
    rank: int,
    weight_grid: tuple[int, int, int],
    weight_lam: float,
    seed: int = 20260731,
) -> tuple[np.ndarray, np.ndarray | None, float]:
    shapes = data.shapes[mask]
    counts = np.clip(data.counts[mask], 1.0, WEIGHT_FIELD_COUNT_CAP)
    if rank == 1:
        mean_shape = (
            (counts[:, None] * shapes).sum(axis=0) / counts.sum())
        mean_shape /= mean_shape.sum()
        return mean_shape[None, :], None, 0.0
    mixing, archetypes, nmf_error = weighted_nmf(
        shapes, counts, rank, seed=seed)
    design = tensor_design(
        data.features[mask], domain, weight_grid)
    penalty = second_difference_penalty(
        weight_grid, (weight_lam, weight_lam, weight_lam))
    weight_norm = counts / counts.mean()

    def fit_fields(weights: np.ndarray) -> np.ndarray:
        log_odds_targets = np.log(
            (weights[:, :-1] + LOG_ODDS_EPS)
            / (weights[:, -1:] + LOG_ODDS_EPS))
        fields = np.empty((rank - 1, int(np.prod(weight_grid))))
        for index in range(rank - 1):
            fields[index] = solve_penalized(
                design, log_odds_targets[:, index], weight_norm,
                penalty)
        return fields

    def softmax_weights(fields: np.ndarray) -> np.ndarray:
        logits = np.zeros((design.shape[0], rank))
        for index in range(rank - 1):
            logits[:, index] = design @ fields[index]
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        return weights / weights.sum(axis=1, keepdims=True)

    weights = mixture_weights(shapes, archetypes)
    # Reference archetype (gauge l_{R-1} = 0) = most-used archetype.
    usage = (counts[:, None] * weights).sum(axis=0)
    order = np.argsort(usage)  # ascending; last = most used
    archetypes = archetypes[order]
    weights = weights[:, order]
    fields = fit_fields(weights)
    # One field-consistent alternation: refit the archetypes against the
    # SMOOTHED softmax weights (what the payload actually mixes with),
    # renormalize rows to unit sum (spec), then refit the fields against
    # NNLS weights on the refined archetypes. Directly reduces the
    # deployed-shape TV rather than the NNLS-weight surrogate.
    smoothed = softmax_weights(fields)
    scaled = smoothed * np.sqrt(counts)[:, None]
    refined = np.empty_like(archetypes)
    for bin_index in range(BASIS_SIZE):
        refined[:, bin_index] = robust_nnls(
            scaled, shapes[:, bin_index] * np.sqrt(counts))
    row_sums = refined.sum(axis=1)
    degenerate = row_sums <= 1e-12
    refined[degenerate] = archetypes[degenerate]
    row_sums = refined.sum(axis=1)
    archetypes = refined / row_sums[:, None]
    weights = mixture_weights(shapes, archetypes)
    usage = (counts[:, None] * weights).sum(axis=0)
    order = np.argsort(usage)
    archetypes = archetypes[order]
    weights = weights[:, order]
    fields = fit_fields(weights)
    return archetypes, fields, nmf_error


# ---------------------------------------------------------------------------
# Selection metrics (fast paths on cached arrays + spectra)
# ---------------------------------------------------------------------------

def eval_yield_metrics(
    component: StructuredComponent,
    data: ComponentData,
    mask: np.ndarray,
) -> dict:
    """|relative yield error| summaries on masked cells, from log space:
    err = exp(model - truth_log) - 1 (base cancels)."""
    if not np.any(mask):
        return {"count": 0}
    conditional = component.predict_yield(data.features[mask])
    if component.name == "positron_containment":
        # deployment clamps the probability at 1 and is judged against the
        # raw (unsmoothed) containment truth, matching fit.evaluate
        predicted = np.minimum(1.0, data.base[mask] * conditional)
        truth = np.maximum(data.raw_total[mask], 1e-300)
        errors = np.abs(predicted - data.raw_total[mask]) / truth
    else:
        predicted = data.base[mask] * conditional
        truth = np.maximum(data.raw_total[mask], 1e-300)
        errors = np.abs(predicted - data.raw_total[mask]) / truth
    return fit.finite_summary([float(v) for v in errors])


def eval_shape_tv(
    component: StructuredComponent,
    data: ComponentData,
    mask: np.ndarray,
    groups: int = 32,
) -> dict:
    """Grouped total variation between per-cell NNLS truth shapes and the
    model shapes, both pushed through the 32-term basis coefficients
    directly (basis columns are unit-sum, so coefficient-space TV grouped
    over coefficients matches spectrum TV closely and is fast)."""
    if not np.any(mask):
        return {"count": 0}
    predicted = component.predict_shape(data.features[mask])
    truth = data.shapes[mask]
    tv = 0.5 * np.abs(predicted - truth).sum(axis=1)
    return fit.finite_summary([float(v) for v in tv])


def selection_objective(summary: dict) -> float:
    if summary.get("count", 0) == 0:
        return math.inf
    return summary["median"] + 0.3 * summary["p95"]


# ---------------------------------------------------------------------------
# Bank loading with caches
# ---------------------------------------------------------------------------

def bank_signature(paths: list[Path]) -> list[tuple[str, int, int]]:
    """Name, size and modification time of every cell file: a cache built
    from a different set of cells (a regenerated or extended bank) never
    matches."""
    answer = []
    for path in paths:
        status = path.stat()
        answer.append((path.name, status.st_size, status.st_mtime_ns))
    return answer


def read_cell_cache(path: Path) -> list[fit.Cell]:
    """Cells parsed by `load_bank`, read back without the bank directory
    (the joint shape fit reads the structured fit's caches this way)."""
    if not path.exists():
        raise SystemExit(
            f"{path} is missing: it is written by "
            "fit_charged_particle_source_structured.py, run that first")
    with path.open("rb") as stream:
        blob = pickle.load(stream)
    if not isinstance(blob, dict) or "cells" not in blob:
        raise SystemExit(
            f"{path} is not a cell cache from this fitter; delete it and "
            "re-run fit_charged_particle_source_structured.py")
    return blob["cells"]


def load_bank(
    campaign: Path,
    helper: Path,
    cache_dir: Path,
    cache_name: str,
    holdout_z: set[int],
) -> tuple[list[fit.Cell], bool]:
    """Parse a bank's cells (with cross sections attached), reusing
    cache/<cache_name>.pkl when it was built from exactly these cell files
    with the same cross-section helper binary and query options (the cached
    cells carry what the helper returned).  Returns the cells and whether
    they were parsed afresh."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{cache_name}.pkl"
    paths = bank_io.cell_files(campaign / "cells")
    if not paths:
        raise SystemExit(f"no campaign cells in {campaign}/cells")
    # The fluorescence bases come from the charged evaluator next to the
    # cross-section helper.
    evaluator = Path(helper).with_name("ceelo_charged_particle_model_evaluate")
    signature = {
        "cells": bank_signature(paths),
        "helper_sha256":
            hashlib.sha256(Path(helper).read_bytes()).hexdigest(),
        "evaluator_sha256":
            hashlib.sha256(evaluator.read_bytes()).hexdigest(),
        "helper_options": list(fit.HELPER_OPTIONS)}
    cells = None
    if cache_path.exists():
        with cache_path.open("rb") as stream:
            blob = pickle.load(stream)
        if isinstance(blob, dict) and blob.get("signature") == signature:
            cells = blob["cells"]
    parsed = cells is None
    if parsed:
        cells = [fit.load_cell(path, set()) for path in paths]
        fit.attach_properties(cells, helper)
        fit.attach_photon_absorption(cells, helper)
        fit.attach_fluorescence_bases(cells, evaluator)
        with cache_path.open("wb") as stream:
            pickle.dump(
                {"signature": signature, "cells": cells}, stream,
                protocol=4)
    for cell in cells:
        cell.split = fit.split_for(cell.identity, cell.z_value, holdout_z)
    return cells, parsed


def load_shape_cache(path: Path) -> dict:
    if path.exists():
        with path.open("rb") as stream:
            return pickle.load(stream)
    return {}


def save_shape_cache(path: Path, cache: dict) -> None:
    with path.open("wb") as stream:
        pickle.dump(cache, stream, protocol=4)


# ---------------------------------------------------------------------------
# Self tests
# ---------------------------------------------------------------------------

def selftest() -> None:
    rng = np.random.default_rng(1)
    for n in (4, 6, 10, 14, 20, 24, 32):
        u = np.concatenate([
            rng.uniform(0.0, 1.0, 257), [0.0, 0.5, 1.0]])
        basis = bspline_basis_matrix(u, n)
        assert np.allclose(basis.sum(axis=1), 1.0, atol=1e-12), n
        assert np.all(basis >= -1e-15), n
        # window support: <= 4 consecutive nonzeros
        start = np.clip(
            np.floor(u * (n - 3)).astype(int), 0, n - 4)
        for row in range(u.size):
            outside = np.ones(n, dtype=bool)
            outside[start[row]:start[row] + 4] = False
            assert np.all(np.abs(basis[row][outside]) < 1e-14), n
    # exact agreement with fit_charged_particle_source's basis at n == 32
    for u in np.concatenate([rng.uniform(0, 1, 64), [0.0, 1.0, 0.999]]):
        ours = bspline_basis_matrix(np.asarray([u]), 32)[0]
        theirs = fit.cubic_bspline_values(float(u))
        assert np.allclose(ours, theirs, atol=1e-14), u
    # tensor design partition of unity
    features = rng.uniform(-0.2, 1.2, size=(300, 3))
    domain = np.asarray([[0.0, 1.0]] * 3)
    design = tensor_design(features, domain, (8, 6, 10))
    assert np.allclose(
        np.asarray(design.sum(axis=1)).ravel(), 1.0, atol=1e-12)
    # penalized solve recovers a smooth function
    truth = (
        np.sin(2.0 * clamp01(features[:, 0]))
        + clamp01(features[:, 1]) ** 2
        + 0.5 * clamp01(features[:, 2]))
    coefficients = solve_penalized(
        design, truth, np.ones(truth.size),
        second_difference_penalty((8, 6, 10), (1e-5,) * 3))
    fitted = design @ coefficients
    assert float(np.max(np.abs(fitted - truth))) < 0.01
    print("selftest passed", flush=True)


# ---------------------------------------------------------------------------
# Selection + fitting driver
# ---------------------------------------------------------------------------

YIELD_N0 = (8, 12, 16)
YIELD_N1 = (6, 8, 10)
YIELD_N2 = (14, 20, 26)
YIELD_N2_CONTAINMENT = (20, 26, 32)
YIELD_LAMBDA = (1e-5, 1e-4, 1e-3)
YIELD_LAMBDA_CONTAINMENT = (1e-7, 1e-6, 1e-5)
SHAPE_RANKS = (2, 3, 4, 6, 8)
WEIGHT_GRIDS = ((8, 6, 10), (10, 6, 12), (12, 8, 14))
WEIGHT_LAMBDA = (1e-3, 1e-2)
PARSIMONY_SLACK = 0.02
# selection-only quality gates: score structure candidates on cells whose
# truth is statistically resolved, so the objective responds to structure
# instead of the MC-noise floor (reports use the standard gates)
SELECT_MAX_SIGMA_REL = 0.05
SELECT_MIN_SHAPE_COUNTS = 200.0


def select_yield_structure(
    name: str,
    data: ComponentData,
    train_mask: np.ndarray,
    val_mask: np.ndarray,
    offgrid_data: ComponentData | None,
    log,
) -> tuple[tuple[int, int, int], float, dict]:
    domain = component_domain(name, data, train_mask)
    select_val = val_mask & (data.sigma_rel <= SELECT_MAX_SIGMA_REL)
    if not np.any(select_val):
        select_val = val_mask
    offgrid_select = None
    if offgrid_data is not None and offgrid_data.features.size:
        offgrid_select = (
            offgrid_data.sigma_rel <= SELECT_MAX_SIGMA_REL)
        if not np.any(offgrid_select):
            offgrid_select = np.ones(
                offgrid_data.log_yield.size, dtype=bool)
    n2_candidates = (
        YIELD_N2_CONTAINMENT if name == "positron_containment"
        else YIELD_N2)
    lambda_candidates = (
        YIELD_LAMBDA_CONTAINMENT if name == "positron_containment"
        else YIELD_LAMBDA)
    current = [12, 6, n2_candidates[1], lambda_candidates[1]]

    def score(n0, n1, n2, lam):
        grid = (int(n0), int(n1), int(n2))
        coefficients = fit_yield_field(
            data, train_mask, domain, grid, lam)
        component = StructuredComponent(
            name=name, domain=domain,
            yield_grid=grid, yield_coefficients=coefficients)
        errors = []
        summary = eval_yield_metrics(component, data, select_val)
        if summary.get("count", 0):
            errors.append((summary, summary["count"]))
        if offgrid_select is not None:
            off = eval_yield_metrics(
                component, offgrid_data, offgrid_select)
            if off.get("count", 0):
                errors.append((off, off["count"]))
        if not errors:
            return math.inf, grid
        objective = sum(
            selection_objective(summary) * count
            for summary, count in errors
        ) / sum(count for _, count in errors)
        return objective, grid

    best = math.inf
    history = {}
    for sweep in range(2):
        improved = False
        for axis, candidates in (
            (2, n2_candidates), (3, lambda_candidates),
            (0, YIELD_N0), (1, YIELD_N1),
        ):
            results = []
            for value in candidates:
                trial = list(current)
                trial[axis] = value
                key = tuple(trial)
                if key not in history:
                    history[key], _ = score(*trial)
                results.append((history[key], value))
            # parsimony: smallest grid within slack of the best
            best_obj = min(result[0] for result in results)
            chosen = None
            for objective, value in sorted(
                results,
                key=lambda r: (
                    r[1] if axis < 3 else -math.log10(r[1]))):
                if objective <= best_obj * (1.0 + PARSIMONY_SLACK):
                    chosen = value
                    break
            if chosen is None:
                chosen = min(results)[1]
            if current[axis] != chosen:
                improved = True
            current[axis] = chosen
            best = history[tuple(current)]
        log(f"  {name} yield sweep {sweep}: "
            f"grid=({current[0]},{current[1]},{current[2]}) "
            f"lambda={current[3]:g} objective={best:.4f}")
        if not improved:
            break
    grid = (int(current[0]), int(current[1]), int(current[2]))
    return grid, float(current[3]), {"objective": best}


def select_shape_structure(
    name: str,
    data: ComponentData,
    train_mask: np.ndarray,
    val_mask: np.ndarray,
    offgrid_data: ComponentData | None,
    domain: np.ndarray,
    log,
) -> tuple[int, tuple[int, int, int], float, dict]:
    select_val = val_mask & (data.counts >= SELECT_MIN_SHAPE_COUNTS)
    if not np.any(select_val):
        select_val = val_mask
    offgrid_select = None
    if offgrid_data is not None and offgrid_data.features.size:
        offgrid_select = (
            offgrid_data.counts >= SELECT_MIN_SHAPE_COUNTS)
        if not np.any(offgrid_select):
            offgrid_select = np.ones(
                offgrid_data.log_yield.size, dtype=bool)

    def score(rank, weight_grid, weight_lam):
        archetypes, fields, _ = fit_shape_model(
            data, train_mask, domain, rank, weight_grid, weight_lam)
        component = StructuredComponent(
            name=name, domain=domain,
            rank=rank, archetypes=archetypes,
            weight_grid=weight_grid, log_odds=fields)
        summaries = []
        local = eval_shape_tv(component, data, select_val)
        if local.get("count", 0):
            summaries.append((local, local["count"]))
        if offgrid_select is not None:
            off = eval_shape_tv(
                component, offgrid_data, offgrid_select)
            if off.get("count", 0):
                summaries.append((off, off["count"]))
        if not summaries:
            return math.inf
        return sum(
            selection_objective(summary) * count
            for summary, count in summaries
        ) / sum(count for _, count in summaries)

    current = [4, (10, 6, 12), 1e-2]
    history = {}

    def cached(rank, weight_grid, weight_lam):
        key = (rank, weight_grid, weight_lam)
        if key not in history:
            history[key] = score(rank, weight_grid, weight_lam)
        return history[key]

    for sweep in range(2):
        improved = False
        # rank
        results = [
            (cached(rank, current[1], current[2]), rank)
            for rank in SHAPE_RANKS]
        best_obj = min(result[0] for result in results)
        for objective, rank in sorted(results, key=lambda r: r[1]):
            if objective <= best_obj + PARSIMONY_SLACK * max(
                best_obj, 0.01):
                if current[0] != rank:
                    improved = True
                current[0] = rank
                break
        # weight grid (only matters for rank >= 2)
        if current[0] >= 2:
            results = [
                (cached(current[0], grid, current[2]), grid)
                for grid in WEIGHT_GRIDS]
            best_obj = min(result[0] for result in results)
            for objective, grid in sorted(
                results, key=lambda r: int(np.prod(r[1]))):
                if objective <= best_obj + PARSIMONY_SLACK * max(
                    best_obj, 0.01):
                    if current[1] != grid:
                        improved = True
                    current[1] = grid
                    break
            results = [
                (cached(current[0], current[1], lam), lam)
                for lam in WEIGHT_LAMBDA]
            chosen = min(results)[1]
            if current[2] != chosen:
                improved = True
            current[2] = chosen
        objective = cached(*current)
        log(f"  {name} shape sweep {sweep}: R={current[0]} "
            f"grid={current[1]} lambda={current[2]:g} "
            f"objective={objective:.4f}")
        if not improved:
            break
    return current[0], current[1], float(current[2]), {
        "objective": cached(*current)}


def build_component(
    name: str,
    data: ComponentData,
    fit_mask: np.ndarray,
    structure: dict,
) -> StructuredComponent:
    domain = component_domain(name, data, fit_mask)
    grid = tuple(structure["yield_grid"])
    lam = structure["yield_lambda"]
    coefficients = fit_yield_field(
        data, fit_mask, domain, grid, lam)
    component = StructuredComponent(
        name=name, domain=domain,
        yield_grid=grid, yield_coefficients=coefficients,
        lambdas=(lam, lam, lam))
    if name in YIELD_ONLY_COMPONENTS:
        component.rank = 0
        return component
    rank = structure["shape_rank"]
    weight_grid = tuple(structure["weight_grid"])
    weight_lam = structure["weight_lambda"]
    archetypes, fields, nmf_error = fit_shape_model(
        data, fit_mask, domain, rank, weight_grid, weight_lam)
    component.rank = rank
    component.archetypes = archetypes
    component.weight_grid = weight_grid if rank >= 2 else None
    component.log_odds = fields
    component.weight_lambdas = (weight_lam,) * 3
    return component


def summarize_records(records: list[dict], components) -> dict:
    answer = {}
    for name in components:
        selected = [
            record for record in records
            if record["component"] == name]
        answer[name] = {
            split: {
                metric: fit.finite_summary([
                    abs(float(record[metric]))
                    for record in selected
                    if record["split"] == split])
                for metric in (
                    "yield_relative_uncertainty",
                    "yield_relative_error",
                    "shape_total_variation",
                    "weighted_rms_relative_error",
                    "p95_relative_error")
            }
            for split in sorted({
                record["split"] for record in selected})
        }
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cross-section-helper", type=Path,
        default=Path(
            "build-rel/tools/full_spectrum/"
            "ceelo_source_cross_section_features"))
    parser.add_argument("--holdout-z", default="35,64")
    parser.add_argument(
        "--offgrid-validation", type=Path,
        help="off-grid validation bank, used with the bank's validation "
             "split for structure selection and scored separately")
    parser.add_argument(
        "--offgrid-final-test", type=Path,
        help="off-grid FINAL TEST bank, scored once after everything is "
             "frozen")
    parser.add_argument(
        "--photon-absorption-path-fraction", type=float,
        default=DEFAULT_PAPF)
    parser.add_argument(
        "--structure-json", type=Path,
        help="skip selection; use this per-component structure")
    parser.add_argument(
        "--no-refit-validation", action="store_true",
        help="keep hash-validation cells out of the final fit")
    parser.add_argument(
        "--skip-final-test", action="store_true",
        help="do not score the off-grid FINAL TEST bank")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    if args.selftest:
        selftest()
        return 0

    log = lambda msg: print(msg, flush=True)  # noqa: E731
    papf = args.photon_absorption_path_fraction
    holdout_z = {
        int(value) for value in args.holdout_z.split(",") if value}
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    cache_dir = output / "cache"

    log("loading training bank")
    cells, reparsed = load_bank(
        args.campaign, args.cross_section_helper, cache_dir,
        "bank_cells", holdout_z)
    log(f"loaded {len(cells)} cells "
        f"(train {sum(c.split == 'train' for c in cells)}, "
        f"validation {sum(c.split == 'validation' for c in cells)}, "
        f"test {sum(c.split == 'test' for c in cells)})")
    offgrid_val_cells = []
    if args.offgrid_validation is not None:
        offgrid_val_cells, parsed = load_bank(
            args.offgrid_validation, args.cross_section_helper,
            cache_dir, "offgrid_validation_cells", set())
        reparsed = reparsed or parsed
        for cell in offgrid_val_cells:
            cell.split = "offgrid_validation"
    log(f"loaded {len(offgrid_val_cells)} off-grid validation cells")

    basis_cache: dict = {}
    shape_cache_path = cache_dir / f"shape_nnls_papf{papf:g}.pkl"
    # Per-cell NNLS shapes are keyed by cell identity; a re-parsed bank may
    # hold different truth under the same identities, so start afresh.
    shape_cache = {} if reparsed else load_shape_cache(shape_cache_path)
    shape_cache_size = -1 if reparsed else len(shape_cache)

    log("extracting per-component training data (NNLS shapes)")
    data_by_component: dict[str, ComponentData] = {}
    offgrid_by_component: dict[str, ComponentData] = {}
    for name in MODEL_COMPONENTS:
        data_by_component[name] = component_training_data(
            cells, name, basis_cache, papf,
            shape_cache=shape_cache)
        offgrid_by_component[name] = component_training_data(
            offgrid_val_cells, name, basis_cache, papf,
            shape_cache=shape_cache)
        log(f"  {name}: {data_by_component[name].log_yield.size} "
            f"bank cells, "
            f"{offgrid_by_component[name].log_yield.size} off-grid")
    if len(shape_cache) != shape_cache_size:
        save_shape_cache(shape_cache_path, shape_cache)

    # ------------------------------------------------------------------
    # Structure selection (train mask = 'train'; select on 'validation'
    # union the off-grid validation bank; never the final-test bank)
    # ------------------------------------------------------------------
    structures: dict[str, dict] = {}
    selection_notes: dict[str, dict] = {}
    if args.structure_json:
        structures = json.loads(args.structure_json.read_text())
        log("using structure from " + str(args.structure_json))
    else:
        for name in MODEL_COMPONENTS:
            data = data_by_component[name]
            if data.log_yield.size == 0:
                log(f"no training cells for {name}: left out of the payload")
                continue
            train_mask = data.splits == "train"
            val_mask = data.splits == "validation"
            offgrid = offgrid_by_component[name]
            offgrid = offgrid if offgrid.log_yield.size else None
            log(f"selecting structure for {name}")
            grid, lam, yield_note = select_yield_structure(
                name, data, train_mask, val_mask, offgrid, log)
            structure = {
                "yield_grid": list(grid),
                "yield_lambda": lam,
            }
            note = {"yield": yield_note}
            if name not in YIELD_ONLY_COMPONENTS:
                domain = component_domain(name, data, train_mask)
                rank, weight_grid, weight_lam, shape_note = (
                    select_shape_structure(
                        name, data, train_mask, val_mask, offgrid,
                        domain, log))
                structure.update({
                    "shape_rank": rank,
                    "weight_grid": list(weight_grid),
                    "weight_lambda": weight_lam,
                })
                note["shape"] = shape_note
            structures[name] = structure
            selection_notes[name] = note
    (output / "selected_structure.json").write_text(
        json.dumps(structures, indent=2) + "\n")

    # ------------------------------------------------------------------
    # Final fit (validation folded into train unless disabled)
    # ------------------------------------------------------------------
    log("final fit")
    payload = StructuredPayload(photon_absorption_path_fraction=papf)
    fit_splits = (
        {"train"} if args.no_refit_validation
        else {"train", "validation"})
    for name in MODEL_COMPONENTS:
        data = data_by_component[name]
        fit_mask = np.isin(data.splits, sorted(fit_splits))
        if not np.any(fit_mask) or name not in structures:
            # A component with no populated cell (e.g. shield x-rays in a bank
            # without a heavy element) is written as "none".
            log(f"  {name}: no populated cells, left out")
            continue
        payload.components[name] = build_component(
            name, data, fit_mask, structures[name])
        log(f"  fitted {name} on {int(fit_mask.sum())} cells")
    payload_path = output / "charged_particle_source_payload.txt"
    write_payload(payload_path, payload)
    payload_bytes = payload_path.stat().st_size
    log(f"payload written: {payload_bytes} bytes")

    # round-trip parity: parse the serialized payload back and require
    # predictions to match the in-memory model to float precision
    parsed = load_payload(payload_path)
    rng = np.random.default_rng(7)
    for name in MODEL_COMPONENTS:
        if name not in payload.components:
            continue
        left = payload.components[name]
        right = parsed.components[name]
        probe = rng.uniform(-0.2, 1.2, size=(64, 3)) * (
            left.domain[:, 1] - left.domain[:, 0]) + left.domain[:, 0]
        dy = np.abs(
            left.predict_yield(probe) - right.predict_yield(probe))
        scale = np.abs(left.predict_yield(probe)) + 1e-300
        assert float((dy / scale).max()) < 1e-12, name
        ds = np.abs(
            left.predict_shape(probe) - right.predict_shape(probe))
        assert float(ds.max()) < 1e-12, name
    log("round-trip parity check passed")

    # ------------------------------------------------------------------
    # Scoring through the shared metrics machinery
    # ------------------------------------------------------------------
    predict_fn = payload_predict_fn(parsed)

    log("scoring bank test cells (final model)")
    test_cells = [cell for cell in cells if cell.split == "test"]
    records, _ = fit.evaluate(
        test_cells, predict_fn, basis_cache, papf)
    if not args.no_refit_validation:
        # validation cells were folded into the final fit; score them with
        # a train-only model so the reported validation numbers stay honest
        log("scoring bank validation cells (train-only model)")
        trainonly = StructuredPayload(photon_absorption_path_fraction=papf)
        for name in MODEL_COMPONENTS:
            data = data_by_component[name]
            if name not in structures or not np.any(data.splits == "train"):
                continue
            trainonly.components[name] = build_component(
                name, data, data.splits == "train", structures[name])
        validation_cells = [
            cell for cell in cells if cell.split == "validation"]
        validation_records, _ = fit.evaluate(
            validation_cells, payload_predict_fn(trainonly), basis_cache, papf)
        records = validation_records + records
    else:
        log("scoring bank validation cells (same model)")
        validation_cells = [
            cell for cell in cells if cell.split == "validation"]
        validation_records, _ = fit.evaluate(
            validation_cells, predict_fn, basis_cache, papf)
        records = validation_records + records
    fit.write_metrics(output / "metrics.csv", records)

    log("scoring off-grid validation bank")
    offgrid_val_records, _ = fit.evaluate(
        offgrid_val_cells, predict_fn, basis_cache, papf)
    fit.write_metrics(
        output / "offgrid_validation_metrics.csv",
        offgrid_val_records)

    offgrid_final_summary = None
    if not args.skip_final_test and args.offgrid_final_test is not None:
        log("scoring off-grid FINAL TEST bank (single, final scoring)")
        offgrid_final_cells, _ = load_bank(
            args.offgrid_final_test, args.cross_section_helper,
            cache_dir, "offgrid_final_test_cells", set())
        for cell in offgrid_final_cells:
            cell.split = "offgrid_final_test"
        offgrid_final_records, _ = fit.evaluate(
            offgrid_final_cells, predict_fn, basis_cache, papf)
        fit.write_metrics(
            output / "offgrid_final_test_metrics.csv",
            offgrid_final_records)
        offgrid_final_summary = summarize_records(
            offgrid_final_records, fit.REPORT_COMPONENTS)

    structure_block = {}
    for name in MODEL_COMPONENTS:
        if name not in payload.components:
            continue
        component = payload.components[name]
        entry = {
            "domain": component.domain.tolist(),
            "yield_grid": list(component.yield_grid),
            "yield_lambda": structures[name]["yield_lambda"],
            "shape_rank": component.rank,
        }
        if component.rank >= 2:
            entry["weight_grid"] = list(component.weight_grid)
            entry["weight_lambda"] = structures[name]["weight_lambda"]
        entry["parameter_count"] = (
            int(np.prod(component.yield_grid))
            + component.rank * BASIS_SIZE
            + (component.rank - 1) * int(np.prod(component.weight_grid))
            if component.rank >= 2 else
            int(np.prod(component.yield_grid))
            + component.rank * BASIS_SIZE)
        structure_block[name] = entry

    report = {
        "model": "structured_tensor_bspline_yield_softmax_nmf_shapes",
        "format_version": 1,
        "basis_size": BASIS_SIZE,
        "photon_absorption_path_fraction": papf,
        "payload_bytes": payload_bytes,
        "cells": len(cells),
        "holdout_z": sorted(holdout_z),
        "splits": {
            split: sum(cell.split == split for cell in cells)
            for split in ("train", "validation", "test")},
        "refit_validation": not args.no_refit_validation,
        "structure": structure_block,
        "selection_notes": selection_notes,
        "metrics": summarize_records(records, fit.REPORT_COMPONENTS),
        "offgrid": {
            "validation": summarize_records(
                offgrid_val_records, fit.REPORT_COMPONENTS),
            "final_test": offgrid_final_summary,
        },
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps(
        {k: report[k] for k in (
            "payload_bytes", "structure")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
