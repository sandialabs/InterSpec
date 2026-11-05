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

"""Positive structured continuum fits and portable runtime payloads.

This module is an intentionally compact reference implementation of the fit
contract.  It requires NumPy but not SciPy or scikit-learn.  The transport and
campaign-management modules do not import it, keeping worker installations
dependency-free.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import io
import json
import math
from pathlib import Path
import struct
from typing import Any, Iterable, Mapping, Sequence

try:  # Deliberately optional; see require_numpy().
    import numpy as np
except ImportError:  # pragma: no cover - exercised in the bare-Python CLI test
    np = None  # type: ignore[assignment]

from .axes import registered_edges_from_physical
from .design import PhysicalNode, geometry_transform


PAYLOAD_MAGIC = b"CEEFIT-PAYLOAD1\0"
PAYLOAD_VERSION = 1
PAYLOAD_ENDIAN_TAG = 0x01020304
PAYLOAD_NORMALIZATION_S_HAT = 1
PAYLOAD_HEADER = struct.Struct("<16sIIQ32s")
_AXIS_KIND_TO_ID = {
    "primary": 1,
    "annihilation": 2,
    "brems": 3,
    "line_e0": 4,
    "line_511": 5,
}
_AXIS_ID_TO_KIND = {value: key for key, value in _AXIS_KIND_TO_ID.items()}
PAIR_THRESHOLD_KEV = 2.0 * 510.998950
_MODEL_NAME_TO_AXIS = {
    "primary_order1": "primary",
    "primary_order2plus": "primary",
    "annihilation_order1": "annihilation",
    "annihilation_order2plus": "annihilation",
    "brems_order0": "brems",
    "brems_order1": "brems",
    "brems_order2plus": "brems",
    "elastic_E0": "line_e0",
    "annihilation_511": "line_511",
}
_LINE_AXIS_KINDS = frozenset(("line_e0", "line_511"))


def _validate_model_identity(name: str, axis_kind: str, *, line: bool | None = None) -> None:
    """Reject aliases and unknown taxonomy before they reach a runtime payload."""

    expected = _MODEL_NAME_TO_AXIS.get(name)
    if expected is None:
        raise ValueError(f"unknown runtime model name {name!r}")
    if axis_kind not in _AXIS_KIND_TO_ID:
        raise ValueError(f"unknown runtime axis kind {axis_kind!r}")
    if axis_kind != expected:
        raise ValueError(
            f"runtime model {name!r} requires axis {expected!r}, not {axis_kind!r}"
        )
    if line is not None and (axis_kind in _LINE_AXIS_KINDS) != line:
        collection = "line" if line else "continuum"
        raise ValueError(f"runtime model {name!r} does not belong in the {collection} inventory")


def require_numpy() -> Any:
    if np is None:
        raise RuntimeError(
            "structured fitting requires NumPy; install the fit extras with "
            "`python -m pip install -r tools/full_spectrum/env_scatter_fit/requirements-fit.txt`"
        )
    return np


def _as_2d(values: Any, name: str) -> Any:
    numpy = require_numpy()
    array = numpy.asarray(values, dtype=float)
    if array.ndim != 2 or not numpy.all(numpy.isfinite(array)):
        raise ValueError(f"{name} must be a finite 2-D array")
    return array


def _open_uniform_knots(n_basis: int, degree: int = 3) -> Any:
    numpy = require_numpy()
    if n_basis < degree + 1:
        raise ValueError("not enough spline basis functions for requested degree")
    n_internal = n_basis - degree - 1
    return numpy.concatenate(
        (
            numpy.zeros(degree + 1),
            numpy.linspace(0.0, 1.0, n_internal + 2)[1:-1],
            numpy.ones(degree + 1),
        )
    )


def _bspline_basis_from_knots(values: Any, knots: Any, degree: int = 3) -> Any:
    numpy = require_numpy()
    knot_array = numpy.asarray(knots, dtype=float)
    n_basis = knot_array.size - degree - 1
    if n_basis < degree + 1 or numpy.any(numpy.diff(knot_array) < 0.0):
        raise ValueError("invalid clamped spline knot vector")
    x = numpy.clip(numpy.asarray(values, dtype=float), 0.0, 1.0)
    knots = knot_array
    basis = numpy.zeros((x.size, len(knots) - 1), dtype=float)
    for i in range(len(knots) - 1):
        basis[:, i] = ((x >= knots[i]) & (x < knots[i + 1])).astype(float)
    basis[x == 1.0, n_basis - 1] = 1.0
    for order in range(1, degree + 1):
        next_basis = numpy.zeros((x.size, len(knots) - order - 1), dtype=float)
        for i in range(next_basis.shape[1]):
            left_denom = knots[i + order] - knots[i]
            right_denom = knots[i + order + 1] - knots[i + 1]
            if left_denom > 0.0:
                next_basis[:, i] += (x - knots[i]) / left_denom * basis[:, i]
            if right_denom > 0.0:
                next_basis[:, i] += (knots[i + order + 1] - x) / right_denom * basis[:, i + 1]
        basis = next_basis
    return basis[:, :n_basis]


def _bspline_basis_1d(values: Any, n_basis: int, degree: int = 3) -> Any:
    return _bspline_basis_from_knots(values, _open_uniform_knots(n_basis, degree), degree)


def _unit_pdf_basis(values: Any, knots: Any, degree: int = 3) -> Any:
    """Evaluate B-spline PDFs whose individual integrals on [0,1] are one."""

    numpy = require_numpy()
    knot_array = numpy.asarray(knots, dtype=float)
    basis = _bspline_basis_from_knots(values, knot_array, degree)
    n_basis = knot_array.size - degree - 1
    integrals = numpy.asarray(
        [(knot_array[i + degree + 1] - knot_array[i]) / (degree + 1) for i in range(n_basis)]
    )
    if numpy.any(integrals <= 0.0):
        raise ValueError("spline basis contains a zero-integral function")
    return basis / integrals[None, :]


def _interval_basis_masses(edges: Any, knots: Any, degree: int = 3) -> Any:
    """Integrate unit-normalized cubic basis PDFs over arbitrary intervals."""

    numpy = require_numpy()
    values = numpy.asarray(edges, dtype=float)
    if values.ndim != 1 or values.size < 2 or numpy.any(numpy.diff(values) < 0.0):
        raise ValueError("integration edges must be an ordered vector")
    # Eight-point Gauss-Legendre integrates each polynomial knot segment exactly;
    # split at every knot so no quadrature interval crosses a derivative break.
    nodes, weights = numpy.polynomial.legendre.leggauss(8)
    n_basis = len(knots) - degree - 1
    result = numpy.zeros((values.size - 1, n_basis), dtype=float)
    unique_knots = sorted(set(float(value) for value in knots if 0.0 < value < 1.0))
    for cell, (raw_low, raw_high) in enumerate(zip(values, values[1:])):
        low, high = max(0.0, float(raw_low)), min(1.0, float(raw_high))
        if high <= low:
            continue
        boundaries = [low, *(k for k in unique_knots if low < k < high), high]
        for left, right in zip(boundaries, boundaries[1:]):
            points = 0.5 * (right - left) * nodes + 0.5 * (right + left)
            result[cell] += 0.5 * (right - left) * (
                weights[:, None] * _unit_pdf_basis(points, knots, degree)
            ).sum(axis=0)
    result[abs(result) < 1.0e-15] = 0.0
    return result


def _project_archetypes_to_splines(
    cell_archetypes: Any,
    energy_edges: Any,
    mu_edges: Any,
    energy_knots: Any,
    mu_knots: Any,
    *,
    cell_weights: Any | None = None,
    iterations: int = 300,
) -> Any:
    """Nonnegative multiplicative NNLS projection onto a tensor PDF basis."""

    numpy = require_numpy()
    archetypes = _as_2d(cell_archetypes, "cell_archetypes")
    e_mass = _interval_basis_masses(energy_edges, energy_knots)
    mu_t_edges = 0.5 * (numpy.asarray(mu_edges, dtype=float) + 1.0)
    m_mass = _interval_basis_masses(mu_t_edges, mu_knots)
    expected = e_mass.shape[0] * m_mass.shape[0]
    if archetypes.shape[1] != expected:
        raise ValueError("cell archetypes do not match projection axes")
    if cell_weights is None:
        weight_grid = numpy.ones((e_mass.shape[0], m_mass.shape[0]))
    else:
        weight_grid = numpy.asarray(cell_weights, dtype=float).reshape(
            e_mass.shape[0], m_mass.shape[0]
        )
        positive = weight_grid[weight_grid > 0.0]
        median = float(numpy.median(positive)) if positive.size else 1.0
        weight_grid = numpy.clip(weight_grid, median * 1.0e-3, median * 1.0e3)
    result = numpy.empty((archetypes.shape[0], e_mass.shape[1], m_mass.shape[1]))
    eps = 1.0e-18
    for rank_index, flat in enumerate(archetypes):
        target = flat.reshape(e_mass.shape[0], m_mass.shape[0])
        coefficients = numpy.full(
            (e_mass.shape[1], m_mass.shape[1]),
            1.0 / (e_mass.shape[1] * m_mass.shape[1]),
        )
        for _ in range(iterations):
            predicted = e_mass @ coefficients @ m_mass.T
            numerator = e_mass.T @ (weight_grid * target) @ m_mass
            denominator = e_mass.T @ (weight_grid * predicted) @ m_mass
            updated = coefficients * numerator / numpy.maximum(denominator, eps)
            norm = float(numpy.sum(updated))
            if norm <= 0.0:
                raise RuntimeError("spline projection produced an empty archetype")
            updated /= norm
            if numpy.max(numpy.abs(updated - coefficients)) < 1.0e-9:
                coefficients = updated
                break
            coefficients = updated
        result[rank_index] = coefficients
    return result


def _tensor_spline_features(normalized_q: Any, level: int) -> Any:
    numpy = require_numpy()
    q = _as_2d(normalized_q, "normalized_q")
    if q.shape[1] != 4:
        raise ValueError("geometry coordinates must have four columns")
    n_basis = _geometry_basis_count(level)
    knots = _geometry_knots(level)
    features = numpy.ones((q.shape[0], 1), dtype=float)
    for dim in range(4):
        one_dim = _bspline_basis_from_knots(q[:, dim], knots)
        features = (features[:, :, None] * one_dim[:, None, :]).reshape(q.shape[0], -1)
    return features


def _geometry_basis_count(level: int) -> int:
    """Return the 1-D basis size for a nested dyadic geometry level.

    Level one is the cubic Bernstein space.  Each subsequent level bisects
    every knot interval, so the spaces (and not merely their dimensions) are
    nested: 1, 2, 4, and 8 intervals give 4, 5, 7, and 11 basis functions.
    """

    if level not in {1, 2, 3, 4}:
        raise ValueError("geometry spline level must be one of 1,2,3,4")
    return (1 << (level - 1)) + 3


def _geometry_feature_count(level: int) -> int:
    return _geometry_basis_count(level) ** 4


def _geometry_knots(level: int) -> Any:
    return _open_uniform_knots(_geometry_basis_count(level))


def _cubic_second_derivative_basis(values: Any, knots: Any) -> Any:
    """Evaluate exact second derivatives of a clamped cubic B-spline basis."""

    numpy = require_numpy()
    knot_array = numpy.asarray(knots, dtype=float)
    n_basis = knot_array.size - 4
    linear = _bspline_basis_from_knots(values, knot_array, degree=1)
    result = numpy.zeros((linear.shape[0], n_basis), dtype=float)
    for i in range(n_basis):
        cubic_left = knot_array[i + 3] - knot_array[i]
        cubic_right = knot_array[i + 4] - knot_array[i + 1]
        if cubic_left > 0.0:
            quadratic_left = knot_array[i + 2] - knot_array[i]
            quadratic_right = knot_array[i + 3] - knot_array[i + 1]
            if quadratic_left > 0.0:
                result[:, i] += (
                    6.0 / (cubic_left * quadratic_left) * linear[:, i]
                )
            if quadratic_right > 0.0:
                result[:, i] -= (
                    6.0 / (cubic_left * quadratic_right) * linear[:, i + 1]
                )
        if cubic_right > 0.0:
            quadratic_left = knot_array[i + 3] - knot_array[i + 1]
            quadratic_right = knot_array[i + 4] - knot_array[i + 2]
            if quadratic_left > 0.0:
                result[:, i] -= (
                    6.0 / (cubic_right * quadratic_left) * linear[:, i + 1]
                )
            if quadratic_right > 0.0:
                result[:, i] += (
                    6.0 / (cubic_right * quadratic_right) * linear[:, i + 2]
                )
    return result


def _geometry_roughness_matrix(level: int) -> Any:
    """Exact 1-D integral of B_i'' B_j'' for the requested dyadic level."""

    numpy = require_numpy()
    knots = _geometry_knots(level)
    nodes, weights = numpy.polynomial.legendre.leggauss(4)
    result = numpy.zeros((_geometry_basis_count(level),) * 2, dtype=float)
    unique = sorted(set(map(float, knots)))
    for low, high in zip(unique, unique[1:]):
        if high <= low:
            continue
        half = 0.5 * (high - low)
        points = half * nodes + 0.5 * (high + low)
        derivative = _cubic_second_derivative_basis(points, knots)
        result += half * (derivative.T * weights) @ derivative
    result = 0.5 * (result + result.T)
    # Enforce the analytically exact constant/affine nullspace after floating-
    # point quadrature.  This prevents roundoff-scale curvature from growing
    # under four tensor products at the finest level.
    greville = numpy.asarray(
        [numpy.mean(knots[i + 1 : i + 4]) for i in range(result.shape[0])]
    )
    nullspace = numpy.column_stack((numpy.ones(result.shape[0]), greville))
    projector = numpy.eye(result.shape[0]) - nullspace @ numpy.linalg.solve(
        nullspace.T @ nullspace, nullspace.T
    )
    result = projector @ result @ projector
    eigenvalues, eigenvectors = numpy.linalg.eigh(0.5 * (result + result.T))
    scale = max(float(numpy.max(eigenvalues)), 1.0)
    eigenvalues[eigenvalues < scale * 1.0e-12] = 0.0
    result = (eigenvectors * eigenvalues) @ eigenvectors.T
    result = projector @ result @ projector
    return 0.5 * (result + result.T)


def _geometry_gram_matrix(level: int) -> Any:
    """Exact 1-D L2 Gram matrix of the geometry B-spline basis."""

    numpy = require_numpy()
    knots = _geometry_knots(level)
    nodes, weights = numpy.polynomial.legendre.leggauss(4)
    result = numpy.zeros((_geometry_basis_count(level),) * 2, dtype=float)
    unique = sorted(set(map(float, knots)))
    for low, high in zip(unique, unique[1:]):
        if high <= low:
            continue
        half = 0.5 * (high - low)
        points = half * nodes + 0.5 * (high + low)
        basis = _bspline_basis_from_knots(points, knots)
        result += half * (basis.T * weights) @ basis
    return 0.5 * (result + result.T)


def _mode_matrix_multiply(tensor: Any, matrix: Any, axis: int) -> Any:
    numpy = require_numpy()
    moved = numpy.moveaxis(tensor, axis, 0)
    transformed = (matrix @ moved.reshape(matrix.shape[1], -1)).reshape(
        (matrix.shape[0],) + moved.shape[1:]
    )
    return numpy.moveaxis(transformed, 0, axis)


def _apply_geometry_roughness(coefficients: Any, level: int) -> Any:
    """Apply the four-axis tensor-spline curvature penalty without forming it."""

    numpy = require_numpy()
    values = numpy.asarray(coefficients, dtype=float)
    was_vector = values.ndim == 1
    if was_vector:
        values = values[:, None]
    feature_count = _geometry_feature_count(level)
    if values.ndim != 2 or values.shape[0] != feature_count:
        raise ValueError("geometry coefficient rows do not match spline level")
    n_basis = _geometry_basis_count(level)
    tensor = values.reshape((n_basis,) * 4 + (values.shape[1],))
    penalty = _geometry_roughness_matrix(level)
    gram = _geometry_gram_matrix(level)
    result = numpy.zeros_like(tensor)
    for axis in range(4):
        term = tensor
        for other_axis in range(4):
            term = _mode_matrix_multiply(
                term, penalty if other_axis == axis else gram, other_axis
            )
        result += term
    flat = result.reshape(values.shape)
    return flat[:, 0] if was_vector else flat


def _normalized_uncertainty_weights(
    variances: Any,
    yields: Any,
    yield_variances: Any | None = None,
) -> tuple[Any, Any]:
    """Return dimensionless shape and yield weights.

    Archive variances are in squared fluence units.  Both fitting targets are
    dimensionless (normalized cell shape and relative yield), so their floors
    must be computed *after* the same normalization.  Mixing a raw-mass floor
    with normalized variances makes otherwise identical relative-precision
    nodes receive different weights solely because their yields differ.
    """

    numpy = require_numpy()
    variance = numpy.asarray(variances, dtype=float)
    total = numpy.asarray(yields, dtype=float).reshape(-1)
    if variance.ndim != 2 or variance.shape[0] != total.size:
        raise ValueError("variance rows must match yields")
    if numpy.any(variance < 0.0) or numpy.any(total <= 0.0):
        raise ValueError("uncertainty weights require nonnegative variance and positive yields")

    shape_variance = variance / total[:, None] ** 2
    positive_shape = shape_variance[shape_variance > 0.0]
    if positive_shape.size:
        shape_floor = max(float(numpy.median(positive_shape)) * 1.0e-9, 1.0e-30)
        shape_weights = 1.0 / numpy.maximum(shape_variance, shape_floor)
        shape_weights /= float(numpy.median(shape_weights))
    else:
        shape_weights = numpy.ones_like(shape_variance)

    if yield_variances is None:
        # Fallback for callers without replica statistics (the single-component
        # NPZ `fit` command).  Archive fitting supplies replica-derived yield
        # variances, because summing cell variances silently discards the
        # positive covariance created when one history scores many cells.
        absolute_yield_variance = numpy.sum(variance, axis=1)
    else:
        absolute_yield_variance = numpy.asarray(yield_variances, dtype=float).reshape(-1)
        if absolute_yield_variance.shape != total.shape or numpy.any(
            absolute_yield_variance < 0.0
        ):
            raise ValueError(
                "yield_variances must be nonnegative and match yields"
            )
    relative_yield_variance = absolute_yield_variance / total**2
    positive_yield = relative_yield_variance[relative_yield_variance > 0.0]
    if positive_yield.size:
        yield_floor = max(float(numpy.median(positive_yield)) * 1.0e-9, 1.0e-30)
        sample_weights = 1.0 / numpy.maximum(relative_yield_variance, yield_floor)
        sample_weights /= float(numpy.median(sample_weights))
    else:
        sample_weights = numpy.ones(total.size, dtype=float)
    return shape_weights, sample_weights


class _GeometryFieldSolver:
    """Reusable exact smoother for all yield/logit targets at one level."""

    def __init__(
        self,
        normalized_q: Any,
        level: int,
        weights: Any,
        roughness_penalty: float,
    ):
        numpy = require_numpy()
        q = _as_2d(normalized_q, "normalized_q")
        if q.shape[1] != 4 or numpy.any(q < 0.0) or numpy.any(q > 1.0):
            raise ValueError("normalized geometry coordinates must lie in [0,1]^4")
        sample_weights = numpy.asarray(weights, dtype=float).reshape(-1)
        if (
            sample_weights.shape[0] != q.shape[0]
            or numpy.any(~numpy.isfinite(sample_weights))
            or numpy.any(sample_weights <= 0.0)
        ):
            raise ValueError("geometry weights must be finite, positive, and match rows")
        if not math.isfinite(roughness_penalty) or roughness_penalty <= 0.0:
            raise ValueError("geometry roughness penalty must be finite and positive")

        self.level = level
        self.row_count = q.shape[0]
        self.n_basis = _geometry_basis_count(level)
        self.sqrt_weights = numpy.sqrt(
            sample_weights / float(numpy.sum(sample_weights))
        )

        # The generalized K v = d G v eigenbasis diagonalizes the exact
        # integral of squared second derivatives over the full 4-D domain.
        gram = _geometry_gram_matrix(level)
        curvature = _geometry_roughness_matrix(level)
        cholesky = numpy.linalg.cholesky(gram)
        left = numpy.linalg.solve(cholesky, curvature)
        whitened = numpy.linalg.solve(cholesky, left.T).T
        eigenvalues, eigenvectors = numpy.linalg.eigh(
            0.5 * (whitened + whitened.T)
        )
        scale = max(float(numpy.max(eigenvalues)), 1.0)
        eigenvalues[eigenvalues < scale * 1.0e-11] = 0.0
        transform = numpy.linalg.solve(cholesky.T, eigenvectors)
        self.transform = transform

        transformed_design = numpy.ones((q.shape[0], 1), dtype=float)
        knots = _geometry_knots(level)
        for dim in range(4):
            one_dim = _bspline_basis_from_knots(q[:, dim], knots) @ transform
            transformed_design = (
                transformed_design[:, :, None] * one_dim[:, None, :]
            ).reshape(q.shape[0], -1)
        roughness_eigenvalues = (
            eigenvalues[:, None, None, None]
            + eigenvalues[None, :, None, None]
            + eigenvalues[None, None, :, None]
            + eigenvalues[None, None, None, :]
        ).reshape(-1)
        # This L2-field ridge is a numerical tie-breaker only in the 16-D
        # multi-affine nullspace.  It is three orders below the curvature
        # coefficient and at most O(1e-8) in the default averaged data loss.
        nullspace_ridge = max(1.0e-12, roughness_penalty * 1.0e-3)
        self.sqrt_inverse_penalty = 1.0 / numpy.sqrt(
            roughness_penalty * roughness_eigenvalues + nullspace_ridge
        )
        self.whitened_design = (
            self.sqrt_weights[:, None]
            * transformed_design
            * self.sqrt_inverse_penalty[None, :]
        )
        normal = numpy.eye(q.shape[0]) + self.whitened_design @ self.whitened_design.T
        self.data_cholesky = numpy.linalg.cholesky(0.5 * (normal + normal.T))

    def solve(self, targets: Any) -> Any:
        numpy = require_numpy()
        y = numpy.asarray(targets, dtype=float)
        if y.ndim == 1:
            y = y[:, None]
        if y.ndim != 2 or y.shape[0] != self.row_count or numpy.any(~numpy.isfinite(y)):
            raise ValueError("geometry targets must be finite and match design rows")
        weighted_targets = self.sqrt_weights[:, None] * y
        intermediate = numpy.linalg.solve(self.data_cholesky, weighted_targets)
        dual = numpy.linalg.solve(self.data_cholesky.T, intermediate)
        spectral = self.sqrt_inverse_penalty[:, None] * (
            self.whitened_design.T @ dual
        )
        tensor = spectral.reshape((self.n_basis,) * 4 + (y.shape[1],))
        for axis in range(4):
            tensor = _mode_matrix_multiply(tensor, self.transform, axis)
        result = tensor.reshape(_geometry_feature_count(self.level), y.shape[1])
        if numpy.any(~numpy.isfinite(result)):
            raise RuntimeError("tensor-spline roughness solve produced non-finite coefficients")
        return result


def _fit_geometry_fields(
    normalized_q: Any,
    level: int,
    targets: Any,
    weights: Any,
    roughness_penalty: float,
) -> Any:
    """Fit fields with exact tensor curvature and a reusable NumPy dual solve."""

    return _GeometryFieldSolver(
        normalized_q, level, weights, roughness_penalty
    ).solve(targets)


@dataclass(frozen=True)
class NMFResult:
    coefficients: Any
    archetypes: Any
    loss_history: tuple[float, ...]


def weighted_nmf(
    values: Any,
    rank: int,
    *,
    weights: Any | None = None,
    iterations: int = 400,
    tolerance: float = 1.0e-7,
    seed: int = 0xCEE10,
) -> NMFResult:
    """Fit a nonnegative basis with weighted multiplicative updates."""

    numpy = require_numpy()
    x = _as_2d(values, "values")
    if numpy.any(x < 0.0):
        raise ValueError("NMF values must be nonnegative")
    if not (1 <= rank <= min(x.shape)):
        raise ValueError("rank must be between one and min(values.shape)")
    if iterations <= 0 or tolerance <= 0.0:
        raise ValueError("iterations and tolerance must be positive")
    omega = numpy.ones_like(x) if weights is None else numpy.asarray(weights, dtype=float)
    if omega.shape != x.shape or numpy.any(~numpy.isfinite(omega)) or numpy.any(omega <= 0.0):
        raise ValueError("NMF weights must be finite, positive, and match values")
    # Extreme inverse-variance weights make multiplicative NMF brittle and let a
    # handful of bins define the entire basis.  Clipping is deterministic and
    # retains six orders of meaningful relative weighting.
    positive = omega[omega > 0.0]
    median = float(numpy.median(positive))
    omega = numpy.clip(omega, median * 1.0e-3, median * 1.0e3)
    rng = numpy.random.default_rng(seed)
    scale = max(float(numpy.mean(x)), 1.0e-12)
    coeff = rng.uniform(0.5, 1.5, size=(x.shape[0], rank)) * math.sqrt(scale / rank)
    basis = rng.uniform(0.5, 1.5, size=(rank, x.shape[1])) * math.sqrt(scale / rank)
    eps = 1.0e-15
    losses: list[float] = []
    previous = math.inf
    for iteration in range(iterations):
        reconstruction = coeff @ basis
        basis *= (coeff.T @ (omega * x)) / (coeff.T @ (omega * reconstruction) + eps)
        reconstruction = coeff @ basis
        coeff *= ((omega * x) @ basis.T) / ((omega * reconstruction) @ basis.T + eps)
        if iteration % 10 == 0 or iteration == iterations - 1:
            residual = omega * (x - coeff @ basis) ** 2
            loss = float(numpy.sum(residual) / numpy.sum(omega))
            losses.append(loss)
            if previous < math.inf and abs(previous - loss) <= tolerance * max(previous, eps):
                break
            previous = loss
    norms = numpy.sum(basis, axis=1)
    if numpy.any(norms <= 0.0):
        raise RuntimeError("NMF produced an empty archetype")
    basis = basis / norms[:, None]
    coeff = coeff * norms[None, :]
    return NMFResult(coeff, basis, tuple(losses))


@dataclass(frozen=True)
class _ShapeDecomposition:
    """Rank-specific shape state shared by every geometry spline level."""

    rank: int
    row_count: int
    energy_knots: Any
    mu_knots: Any
    logits: Any
    archetype_coefficients: Any


def _fit_shape_decomposition(
    shapes: Any,
    axis_kind: str,
    registered_energy_edges: Sequence[float],
    mu_edges: Sequence[float],
    nmf_weights: Any,
    *,
    rank: int,
    nmf_iterations: int,
    seed: int,
    energy_basis_count: int,
    mu_basis_count: int,
    projection_iterations: int,
) -> _ShapeDecomposition:
    """Fit NMF/archetypes once; geometry level changes only their fields."""

    numpy = require_numpy()
    normalized_shapes = _as_2d(shapes, "shapes")
    weights = numpy.asarray(nmf_weights, dtype=float)
    if weights.shape != normalized_shapes.shape:
        raise ValueError("NMF weights must match normalized shapes")
    nmf = weighted_nmf(
        normalized_shapes,
        rank,
        weights=weights,
        iterations=nmf_iterations,
        seed=seed,
    )
    mu_knots = _open_uniform_knots(mu_basis_count)
    is_line = axis_kind in {"line_e0", "line_511"}
    if is_line:
        energy_knots = numpy.asarray([], dtype=float)
        mu_mass = _interval_basis_masses(
            0.5 * (numpy.asarray(mu_edges, dtype=float) + 1.0), mu_knots
        )
        projected = numpy.empty((rank, mu_basis_count))
        projection_weight = numpy.mean(weights, axis=0)
        projection_weight /= max(float(numpy.median(projection_weight)), 1.0e-30)
        for index, target in enumerate(nmf.archetypes):
            coefficients = numpy.full(mu_basis_count, 1.0 / mu_basis_count)
            for _ in range(projection_iterations):
                predicted = mu_mass @ coefficients
                update = coefficients * (
                    mu_mass.T @ (projection_weight * target)
                ) / numpy.maximum(
                    mu_mass.T @ (projection_weight * predicted), 1.0e-18
                )
                update /= max(float(numpy.sum(update)), 1.0e-18)
                if numpy.max(numpy.abs(update - coefficients)) < 1.0e-9:
                    coefficients = update
                    break
                coefficients = update
            projected[index] = coefficients
    else:
        energy_knots = _open_uniform_knots(energy_basis_count)
        projected = _project_archetypes_to_splines(
            nmf.archetypes,
            registered_energy_edges,
            mu_edges,
            energy_knots,
            mu_knots,
            cell_weights=numpy.mean(weights, axis=0),
            iterations=projection_iterations,
        )
    mixture = nmf.coefficients / numpy.maximum(
        numpy.sum(nmf.coefficients, axis=1, keepdims=True), 1.0e-15
    )
    logits = numpy.log(numpy.maximum(mixture, 1.0e-12))
    logits -= numpy.mean(logits, axis=1, keepdims=True)
    return _ShapeDecomposition(
        rank,
        normalized_shapes.shape[0],
        energy_knots,
        mu_knots,
        logits,
        projected,
    )


@dataclass
class StructuredPositiveModel:
    """Positive yield times softmax mixtures of tensor cubic-spline PDFs."""

    component: str
    axis_kind: str
    energy_knots: list[float]
    mu_knots: list[float]
    q_min: list[float]
    q_span: list[float]
    rank: int
    level: int
    roughness_penalty: float
    log_yield_coefficients: Any
    logit_coefficients: Any
    archetype_coefficients: Any
    e_cut_keV: float = 10.0
    pair_threshold_keV: float | None = None

    @property
    def parameter_count(self) -> int:
        """Number of serialized scalar fit coefficients.

        This is the model-selection complexity measure.  In particular, the
        nested four-dimensional geometry field grows with dyadic refinement;
        comparing ``rank * level`` substantially understates higher levels.
        """

        feature_count = _geometry_feature_count(self.level)
        archetype_count = int(require_numpy().asarray(self.archetype_coefficients).size)
        return feature_count * (self.rank + 1) + archetype_count

    @classmethod
    def fit(
        cls,
        component: str,
        axis_kind: str,
        q: Any,
        cell_masses: Any,
        *,
        registered_energy_edges: Sequence[float],
        mu_edges: Sequence[float],
        variances: Any | None = None,
        yield_variances: Any | None = None,
        rank: int = 8,
        level: int = 2,
        roughness_penalty: float = 1.0e-5,
        nmf_iterations: int = 400,
        seed: int = 0xCEE10,
        pair_threshold_keV: float | None = None,
        q_domain: tuple[Sequence[float], Sequence[float]] | None = None,
        energy_basis_count: int = 32,
        mu_basis_count: int = 12,
        projection_iterations: int = 300,
        _field_solver: Any | None = None,
        _shape_decomposition: Any | None = None,
    ) -> "StructuredPositiveModel":
        numpy = require_numpy()
        _validate_model_identity(component, axis_kind)
        coords = _as_2d(q, "q")
        masses = _as_2d(cell_masses, "cell_masses")
        if coords.shape[1] != 4 or coords.shape[0] != masses.shape[0]:
            raise ValueError("q must be N x 4 and match the mass rows")
        if numpy.any(masses < 0.0):
            raise ValueError("cell masses must be nonnegative")
        is_line = axis_kind in {"line_e0", "line_511"}
        expected_cells = (
            len(mu_edges) - 1
            if is_line
            else (len(registered_energy_edges) - 1) * (len(mu_edges) - 1)
        )
        if masses.shape[1] != expected_cells:
            raise ValueError("mass columns do not match registered energy/mu axes")
        yields = numpy.sum(masses, axis=1)
        if numpy.any(yields <= 0.0):
            raise ValueError("every training row must have positive total yield")
        shapes = masses / yields[:, None]
        if q_domain is None:
            q_min = numpy.min(coords, axis=0)
            q_max = numpy.max(coords, axis=0)
        else:
            q_min = numpy.asarray(q_domain[0], dtype=float)
            q_max = numpy.asarray(q_domain[1], dtype=float)
            if q_min.shape != (4,) or q_max.shape != (4,) or numpy.any(q_max < q_min):
                raise ValueError("q_domain must contain two ordered four-vectors")
            if numpy.any(coords < q_min - 1.0e-12) or numpy.any(coords > q_max + 1.0e-12):
                raise ValueError("training q lies outside q_domain")
        q_span = numpy.maximum(q_max - q_min, 1.0e-12)
        normalized_q = (coords - q_min) / q_span

        if variances is None:
            nmf_weights = numpy.ones_like(masses)
            sample_weights = numpy.ones(coords.shape[0])
        else:
            variance = numpy.asarray(variances, dtype=float)
            if variance.shape != masses.shape or numpy.any(variance < 0.0):
                raise ValueError("variances must be nonnegative and match masses")
            nmf_weights, sample_weights = _normalized_uncertainty_weights(
                variance, yields, yield_variances
            )

        if _shape_decomposition is None:
            shape_decomposition = _fit_shape_decomposition(
                shapes,
                axis_kind,
                registered_energy_edges,
                mu_edges,
                nmf_weights,
                rank=rank,
                nmf_iterations=nmf_iterations,
                seed=seed,
                energy_basis_count=energy_basis_count,
                mu_basis_count=mu_basis_count,
                projection_iterations=projection_iterations,
            )
        else:
            if (
                not isinstance(_shape_decomposition, _ShapeDecomposition)
                or _shape_decomposition.rank != rank
                or _shape_decomposition.row_count != coords.shape[0]
            ):
                raise ValueError("reused shape decomposition does not match training data")
            shape_decomposition = _shape_decomposition
        field_targets = numpy.column_stack(
            (numpy.log(yields), shape_decomposition.logits)
        )
        if _field_solver is None:
            fields = _fit_geometry_fields(
                normalized_q,
                level,
                field_targets,
                sample_weights,
                roughness_penalty,
            )
        else:
            if (
                not isinstance(_field_solver, _GeometryFieldSolver)
                or _field_solver.level != level
                or _field_solver.row_count != coords.shape[0]
            ):
                raise ValueError("reused geometry solver does not match training design")
            fields = _field_solver.solve(field_targets)
        yield_coef = fields[:, :1]
        logit_coef = fields[:, 1:]
        return cls(
            component=component,
            axis_kind=axis_kind,
            energy_knots=list(map(float, shape_decomposition.energy_knots)),
            mu_knots=list(map(float, shape_decomposition.mu_knots)),
            q_min=q_min.tolist(),
            q_span=q_span.tolist(),
            rank=rank,
            level=level,
            roughness_penalty=roughness_penalty,
            log_yield_coefficients=yield_coef,
            logit_coefficients=logit_coef,
            archetype_coefficients=shape_decomposition.archetype_coefficients,
            pair_threshold_keV=pair_threshold_keV,
        )

    def _features(self, q: Any) -> Any:
        numpy = require_numpy()
        coords = _as_2d(q, "q")
        normalized = (coords - numpy.asarray(self.q_min)) / numpy.asarray(self.q_span)
        tolerance = 1.0e-12
        if numpy.any(normalized < -tolerance) or numpy.any(normalized > 1.0 + tolerance):
            raise ValueError(f"query lies outside fitted q domain for {self.component}")
        return _tensor_spline_features(numpy.clip(normalized, 0.0, 1.0), self.level)

    def predict_coefficients(self, q: Any) -> tuple[Any, Any, Any]:
        """Return yields, mixture weights, and mixed spline-PDF coefficients."""

        numpy = require_numpy()
        features = self._features(q)
        yields = numpy.exp(features @ numpy.asarray(self.log_yield_coefficients)).reshape(-1)
        logits = features @ numpy.asarray(self.logit_coefficients)
        logits -= numpy.max(logits, axis=1, keepdims=True)
        mixture = numpy.exp(logits)
        mixture /= numpy.sum(mixture, axis=1, keepdims=True)
        archetypes = numpy.asarray(self.archetype_coefficients)
        mixed = mixture @ archetypes.reshape(self.rank, -1)
        mixed /= numpy.sum(mixed, axis=1, keepdims=True)
        return yields, mixture, mixed

    def predict_cell_masses(
        self,
        q: Any,
        registered_energy_edges: Sequence[float],
        mu_edges: Sequence[float],
    ) -> tuple[Any, Any, Any]:
        """Integrate the spline PDFs over caller-provided registered cells."""

        numpy = require_numpy()
        yields, _, mixed = self.predict_coefficients(q)
        mu_mass = _interval_basis_masses(
            0.5 * (numpy.asarray(mu_edges, dtype=float) + 1.0), self.mu_knots
        )
        if self.axis_kind in {"line_e0", "line_511"}:
            shapes = mixed @ mu_mass.T
        else:
            energy_mass = _interval_basis_masses(
                registered_energy_edges, self.energy_knots
            )
            coefficients = mixed.reshape(
                mixed.shape[0], len(self.energy_knots) - 4, len(self.mu_knots) - 4
            )
            shapes = numpy.stack(
                [(energy_mass @ value @ mu_mass.T).reshape(-1) for value in coefficients]
            )
        return yields, shapes, yields[:, None] * shapes

    def to_dict(self) -> dict[str, Any]:
        numpy = require_numpy()
        return {
            "component": self.component,
            "axis_kind": self.axis_kind,
            "energy_knots": self.energy_knots,
            "mu_knots": self.mu_knots,
            "q_min": self.q_min,
            "q_span": self.q_span,
            "rank": self.rank,
            "level": self.level,
            "roughness_penalty": self.roughness_penalty,
            "e_cut_keV": self.e_cut_keV,
            "pair_threshold_keV": self.pair_threshold_keV,
            "log_yield_coefficients": numpy.asarray(self.log_yield_coefficients).tolist(),
            "logit_coefficients": numpy.asarray(self.logit_coefficients).tolist(),
            "archetype_coefficients": numpy.asarray(self.archetype_coefficients).tolist(),
        }

    @classmethod
    def from_dict(cls, record: Mapping[str, Any]) -> "StructuredPositiveModel":
        numpy = require_numpy()
        model = cls(
            component=str(record["component"]),
            axis_kind=str(record["axis_kind"]),
            energy_knots=list(map(float, record["energy_knots"])),
            mu_knots=list(map(float, record["mu_knots"])),
            q_min=list(map(float, record["q_min"])),
            q_span=list(map(float, record["q_span"])),
            rank=int(record["rank"]),
            level=int(record["level"]),
            roughness_penalty=float(
                record.get("roughness_penalty", record.get("ridge_penalty", 0.0))
            ),
            e_cut_keV=float(record.get("e_cut_keV", 10.0)),
            pair_threshold_keV=(
                None
                if record.get("pair_threshold_keV") is None
                else float(record["pair_threshold_keV"])
            ),
            log_yield_coefficients=numpy.asarray(record["log_yield_coefficients"], dtype=float),
            logit_coefficients=numpy.asarray(record["logit_coefficients"], dtype=float),
            archetype_coefficients=numpy.asarray(record["archetype_coefficients"], dtype=float),
        )
        _validate_model_identity(model.component, model.axis_kind)
        return model


@dataclass(frozen=True)
class ModelThresholds:
    max_yield_relative_error: float = 0.005
    max_joint_tv: float = 0.01


@dataclass
class CandidateFit:
    model: StructuredPositiveModel
    max_yield_relative_error: float
    max_joint_tv: float
    mean_joint_tv: float
    passed: bool


def score_model(
    model: StructuredPositiveModel,
    q: Any,
    truth_masses: Any,
    registered_energy_edges: Sequence[float],
    mu_edges: Sequence[float],
    *,
    thresholds: ModelThresholds = ModelThresholds(),
    validation_reference_yields: Any | None = None,
    validation_material_mask: Any | None = None,
    additive_baseline_masses: Any | None = None,
    weak_absolute_yield_fraction: float = 0.001,
    weak_absolute_shape_fraction: float = 0.001,
) -> tuple[float, float, float]:
    """Score validation rows in the currency relevant to the final result.

    The standalone API remains relative-yield plus normalized-shape TV.  An
    archive fit may additionally supply a total-reference yield and a
    conservative materiality mask.  Material rows retain the relative gates;
    weak rows are scored by their absolute yield and angular/joint L1 effect on
    the corresponding total.  ``additive_baseline_masses`` is used for exact
    atoms such as environmental elastic E0, whose delivered observable also
    contains the deterministic direct atom.
    """

    numpy = require_numpy()
    truth = _as_2d(truth_masses, "truth_masses")
    truth_yield = numpy.sum(truth, axis=1)
    if numpy.any(truth_yield <= 0.0):
        raise ValueError("validation truth rows must have positive yield")
    predicted_yield, predicted_shape, _ = model.predict_cell_masses(
        q, registered_energy_edges, mu_edges
    )
    predicted = predicted_yield[:, None] * predicted_shape
    evaluated_truth = truth
    evaluated_prediction = predicted
    if additive_baseline_masses is not None:
        baseline = _as_2d(additive_baseline_masses, "additive_baseline_masses")
        if baseline.shape != truth.shape or numpy.any(~numpy.isfinite(baseline)) or numpy.any(baseline < 0.0):
            raise ValueError("additive validation baseline must be finite, nonnegative, and match truth")
        evaluated_truth = truth + baseline
        evaluated_prediction = predicted + baseline
    evaluated_truth_yield = numpy.sum(evaluated_truth, axis=1)
    evaluated_prediction_yield = numpy.sum(evaluated_prediction, axis=1)
    truth_shape = evaluated_truth / evaluated_truth_yield[:, None]
    predicted_evaluated_shape = (
        evaluated_prediction / evaluated_prediction_yield[:, None]
    )
    yield_error = numpy.abs(
        evaluated_prediction_yield / evaluated_truth_yield - 1.0
    )
    tv = 0.5 * numpy.sum(
        numpy.abs(predicted_evaluated_shape - truth_shape), axis=1
    )

    if validation_reference_yields is not None:
        reference = numpy.asarray(validation_reference_yields, dtype=float).reshape(-1)
        material = numpy.asarray(validation_material_mask, dtype=bool).reshape(-1)
        if (
            reference.shape != truth_yield.shape
            or material.shape != truth_yield.shape
            or numpy.any(~numpy.isfinite(reference))
            or numpy.any(reference <= 0.0)
            or not math.isfinite(weak_absolute_yield_fraction)
            or weak_absolute_yield_fraction <= 0.0
            or not math.isfinite(weak_absolute_shape_fraction)
            or weak_absolute_shape_fraction <= 0.0
        ):
            raise ValueError("contribution-aware validation inputs are invalid")
        weak_yield = (
            numpy.abs(predicted_yield - truth_yield)
            / reference
            * thresholds.max_yield_relative_error
            / weak_absolute_yield_fraction
        )
        weak_shape = (
            0.5 * numpy.sum(numpy.abs(predicted - truth), axis=1)
            / reference
            * thresholds.max_joint_tv
            / weak_absolute_shape_fraction
        )
        yield_error = numpy.where(material, yield_error, weak_yield)
        tv = numpy.where(material, tv, weak_shape)
    elif validation_material_mask is not None:
        raise ValueError("a validation materiality mask requires reference yields")
    return float(numpy.max(yield_error)), float(numpy.max(tv)), float(numpy.mean(tv))


def select_rank_and_level(
    component: str,
    axis_kind: str,
    train_q: Any,
    train_masses: Any,
    validation_q: Any,
    validation_masses: Any,
    *,
    registered_energy_edges: Sequence[float],
    mu_edges: Sequence[float],
    train_variances: Any | None = None,
    train_yield_variances: Any | None = None,
    ranks: Iterable[int] = (4, 8, 12, 16, 24),
    levels: Iterable[int] = (1, 2, 3, 4),
    thresholds: ModelThresholds = ModelThresholds(),
    roughness_penalty: float = 1.0e-5,
    nmf_iterations: int = 400,
    seed: int = 0xCEE10,
    source_energy_threshold_keV: float | None = None,
    q_domain: tuple[Sequence[float], Sequence[float]] | None = None,
    validation_reference_yields: Any | None = None,
    validation_material_mask: Any | None = None,
    additive_baseline_masses: Any | None = None,
    weak_absolute_yield_fraction: float = 0.001,
    weak_absolute_shape_fraction: float = 0.001,
) -> tuple[CandidateFit, tuple[CandidateFit, ...]]:
    """Choose the least-complex rank/level candidate that passes validation."""

    results: list[CandidateFit] = []
    _validate_model_identity(component, axis_kind)
    train = _as_2d(train_masses, "train_masses")
    train_coords = _as_2d(train_q, "train_q")
    numpy = require_numpy()
    all_q = numpy.vstack((train_coords, _as_2d(validation_q, "validation_q")))
    if source_energy_threshold_keV is None:
        if axis_kind in {"annihilation", "line_511"}:
            source_energy_threshold_keV = PAIR_THRESHOLD_KEV
        elif axis_kind == "brems":
            raise ValueError(
                "brems fitting requires the campaign brems activation threshold"
            )
    if q_domain is None:
        q_domain = (numpy.min(all_q, axis=0), numpy.max(all_q, axis=0))
    feasible_ranks = [rank for rank in sorted(set(ranks)) if rank <= min(train.shape)]
    if not feasible_ranks:
        raise ValueError("no requested rank is feasible for the training matrix")
    q_low = numpy.asarray(q_domain[0], dtype=float)
    q_high = numpy.asarray(q_domain[1], dtype=float)
    if q_low.shape != (4,) or q_high.shape != (4,) or numpy.any(q_high < q_low):
        raise ValueError("q_domain must contain two ordered four-vectors")
    q_span = numpy.maximum(q_high - q_low, 1.0e-12)
    normalized_train_q = (train_coords - q_low) / q_span
    training_yields = numpy.sum(train, axis=1)
    if numpy.any(training_yields <= 0.0):
        raise ValueError("every training row must have positive total yield")
    if train_variances is None:
        nmf_weights = numpy.ones_like(train)
        field_weights = numpy.ones(train.shape[0])
    else:
        variance = numpy.asarray(train_variances, dtype=float)
        if variance.shape != train.shape or numpy.any(variance < 0.0):
            raise ValueError("train_variances must be nonnegative and match train_masses")
        nmf_weights, field_weights = _normalized_uncertainty_weights(
            variance, training_yields, train_yield_variances
        )
    normalized_shapes = train / training_yields[:, None]
    shape_decompositions = {
        rank: _fit_shape_decomposition(
            normalized_shapes,
            axis_kind,
            registered_energy_edges,
            mu_edges,
            nmf_weights,
            rank=rank,
            nmf_iterations=nmf_iterations,
            seed=seed + rank * 101,
            energy_basis_count=32,
            mu_basis_count=12,
            projection_iterations=300,
        )
        for rank in feasible_ranks
    }
    for level in sorted(set(levels)):
        field_solver = _GeometryFieldSolver(
            normalized_train_q, level, field_weights, roughness_penalty
        )
        for rank in feasible_ranks:
            model = StructuredPositiveModel.fit(
                component,
                axis_kind,
                train_coords,
                train,
                registered_energy_edges=registered_energy_edges,
                mu_edges=mu_edges,
                variances=train_variances,
                yield_variances=train_yield_variances,
                rank=rank,
                level=level,
                roughness_penalty=roughness_penalty,
                nmf_iterations=nmf_iterations,
                seed=seed + rank * 101,
                q_domain=q_domain,
                pair_threshold_keV=source_energy_threshold_keV,
                _field_solver=field_solver,
                _shape_decomposition=shape_decompositions[rank],
            )
            yield_error, joint_tv, mean_tv = score_model(
                model,
                validation_q,
                validation_masses,
                registered_energy_edges,
                mu_edges,
                thresholds=thresholds,
                validation_reference_yields=validation_reference_yields,
                validation_material_mask=validation_material_mask,
                additive_baseline_masses=additive_baseline_masses,
                weak_absolute_yield_fraction=weak_absolute_yield_fraction,
                weak_absolute_shape_fraction=weak_absolute_shape_fraction,
            )
            results.append(
                CandidateFit(
                    model,
                    yield_error,
                    joint_tv,
                    mean_tv,
                    yield_error <= thresholds.max_yield_relative_error
                    and joint_tv <= thresholds.max_joint_tv,
                )
            )
    passing = [item for item in results if item.passed]
    if passing:
        selected = min(
            passing,
            key=lambda item: (
                item.model.parameter_count,
                item.model.rank,
                item.model.level,
            ),
        )
    else:
        selected = min(
            results,
            key=lambda item: (
                item.max_yield_relative_error / thresholds.max_yield_relative_error
                + item.max_joint_tv / thresholds.max_joint_tv,
                item.model.parameter_count,
                item.model.rank,
                item.model.level,
            ),
        )
    return selected, tuple(results)


@dataclass
class RuntimePayload:
    """Compact fitted continuum payload plus exact-energy angular atoms."""

    components: dict[str, StructuredPositiveModel]
    lines: dict[str, StructuredPositiveModel] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    schema: str = "CEEFIT-PAYLOAD-v1"
    accepted_for_publication: bool = False
    model_content_sha256: str = ""

    @staticmethod
    def _put_string(writer: io.BytesIO, value: str) -> None:
        raw = value.encode("utf-8")
        writer.write(struct.pack("<Q", len(raw)))
        writer.write(raw)

    @staticmethod
    def _put_vector(writer: io.BytesIO, values: Any) -> None:
        numpy = require_numpy()
        flat = numpy.asarray(values, dtype="<f8").reshape(-1)
        if numpy.any(~numpy.isfinite(flat)):
            raise ValueError("binary payload arrays must be finite")
        writer.write(struct.pack("<Q", flat.size))
        writer.write(flat.tobytes(order="C"))

    @classmethod
    def _put_model(cls, writer: io.BytesIO, model: StructuredPositiveModel) -> None:
        numpy = require_numpy()
        _validate_model_identity(model.component, model.axis_kind)
        axis_id = _AXIS_KIND_TO_ID[model.axis_kind]
        expected_features = _geometry_feature_count(model.level)
        yield_coefficients = numpy.asarray(model.log_yield_coefficients, dtype=float).reshape(-1)
        logits = numpy.asarray(model.logit_coefficients, dtype=float)
        archetypes = numpy.asarray(model.archetype_coefficients, dtype=float)
        n_energy_basis = len(model.energy_knots) - 4
        n_mu_basis = len(model.mu_knots) - 4
        if model.rank <= 0 or model.level not in {1, 2, 3, 4}:
            raise ValueError(f"{model.component} rank/spline level is invalid")
        if yield_coefficients.size != expected_features:
            raise ValueError(f"{model.component} yield coefficient count disagrees with spline level")
        if logits.shape != (expected_features, model.rank):
            raise ValueError(f"{model.component} logit coefficient shape is invalid")
        expected_shape = (
            (model.rank, n_mu_basis)
            if model.axis_kind in _LINE_AXIS_KINDS
            else (model.rank, n_energy_basis, n_mu_basis)
        )
        if archetypes.shape != expected_shape:
            raise ValueError(f"{model.component} archetype shape is invalid")
        if numpy.any(~numpy.isfinite(archetypes)) or numpy.any(archetypes < 0.0) or not numpy.allclose(
            numpy.sum(archetypes.reshape(model.rank, -1), axis=1),
            1.0,
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ValueError(f"{model.component} archetypes are not nonnegative unit masses")
        cls._put_string(writer, model.component)
        writer.write(
            struct.pack(
                "<IdB",
                axis_id,
                model.e_cut_keV,
                model.pair_threshold_keV is not None,
            )
        )
        writer.write(struct.pack("<dII", model.pair_threshold_keV or 0.0, model.rank, model.level))
        cls._put_vector(writer, model.energy_knots)
        cls._put_vector(writer, model.mu_knots)
        writer.write(struct.pack("<4d", *model.q_min))
        writer.write(struct.pack("<4d", *model.q_span))
        cls._put_vector(writer, yield_coefficients)
        cls._put_vector(writer, logits)
        cls._put_vector(writer, archetypes)

    def _canonical_inventory(self) -> tuple[tuple[int, str, str], ...]:
        record = self.metadata.get("model_inventory")
        if not isinstance(record, Mapping) or set(record) != {"continuum", "lines"}:
            raise ValueError(
                "runtime payload metadata needs explicit continuum/line model_inventory"
            )
        entries: list[tuple[int, str, str]] = []
        for collection_id, key in enumerate(("continuum", "lines")):
            values = record[key]
            if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
                raise ValueError(f"model_inventory {key} must be a sequence")
            for value in values:
                if not isinstance(value, Mapping) or set(value) != {"name", "axis_kind"}:
                    raise ValueError(f"model_inventory {key} entry is invalid")
                name = str(value["name"])
                axis_kind = str(value["axis_kind"])
                _validate_model_identity(name, axis_kind, line=bool(collection_id))
                entries.append((collection_id, name, axis_kind))
        if len(entries) != len({(collection, name) for collection, name, _ in entries}):
            raise ValueError("model_inventory contains duplicate names")
        return tuple(sorted(entries, key=lambda item: (item[0], item[1])))

    def _validate_runtime_contract(self) -> tuple[tuple[int, str, str], ...]:
        if type(self.accepted_for_publication) is not bool:
            raise ValueError("accepted_for_publication must be an explicit Boolean")
        inventory = self._canonical_inventory()
        actual: list[tuple[int, str, str]] = []
        for collection_id, models in enumerate((self.components, self.lines)):
            if not isinstance(models, Mapping):
                raise ValueError("runtime payload model collections must be mappings")
            for key, model in models.items():
                if key != model.component:
                    raise ValueError(
                        f"runtime model dictionary key {key!r} disagrees with {model.component!r}"
                    )
                _validate_model_identity(model.component, model.axis_kind, line=bool(collection_id))
                actual.append((collection_id, model.component, model.axis_kind))
        if tuple(sorted(actual, key=lambda item: (item[0], item[1]))) != inventory:
            raise ValueError("runtime model inventory does not exactly match fitted models")

        try:
            brems_threshold = float(self.metadata["brems_activation_threshold_keV"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                "runtime payload metadata needs brems_activation_threshold_keV"
            ) from error
        if not math.isfinite(brems_threshold) or brems_threshold < 0.0:
            raise ValueError("brems activation threshold must be finite and nonnegative")
        for _, model in [*self.components.items(), *self.lines.items()]:
            threshold = model.pair_threshold_keV
            if model.axis_kind == "brems":
                expected = brems_threshold if brems_threshold > 0.0 else None
            elif model.axis_kind in {"annihilation", "line_511"}:
                expected = PAIR_THRESHOLD_KEV
            else:
                expected = None
            if expected is None:
                if threshold is not None:
                    raise ValueError(f"{model.component} must not have an activation threshold")
            elif threshold is None or not math.isclose(
                threshold, expected, rel_tol=0.0, abs_tol=1.0e-9
            ):
                raise ValueError(
                    f"{model.component} activation threshold disagrees with payload physics metadata"
                )
        return inventory

    def _binary_identity_fields(
        self,
    ) -> tuple[str, str, str, tuple[float, ...], tuple[float, ...], float, tuple[tuple[int, str, str], ...]]:
        inventory = self._validate_runtime_contract()
        config_hash = str(self.metadata.get("config_sha256", ""))
        if len(config_hash) != 64 or any(c not in "0123456789abcdef" for c in config_hash):
            raise ValueError("binary payload metadata needs a lowercase 64-character config_sha256")
        material_name = str(self.metadata.get("material_name", ""))
        geometry_tier = str(self.metadata.get("geometry_tier", ""))
        if not material_name or not geometry_tier:
            raise ValueError("binary payload metadata needs material_name and geometry_tier")
        if self.metadata.get("normalization") != "s_hat_to_direct_fluence":
            raise ValueError("binary payload normalization must be s_hat_to_direct_fluence")
        physical_domain = self.metadata.get("physical_domain")
        if (
            not isinstance(physical_domain, Mapping)
            or len(physical_domain.get("min", ())) != 4
            or len(physical_domain.get("max", ())) != 4
        ):
            raise ValueError("binary payload metadata needs four-dimensional physical_domain")
        physical_min = tuple(map(float, physical_domain["min"]))
        physical_max = tuple(map(float, physical_domain["max"]))
        if any(not math.isfinite(value) for value in (*physical_min, *physical_max)) or any(
            high < low for low, high in zip(physical_min, physical_max)
        ):
            raise ValueError("binary payload physical_domain is invalid")
        if physical_min[0] <= 0.0 or physical_min[1] <= 0.0 or physical_min[2] < 0.0 or physical_min[3] < 0.0:
            raise ValueError("binary payload physical_domain has invalid physical minima")
        return (
            config_hash,
            material_name,
            geometry_tier,
            physical_min,
            physical_max,
            float(self.metadata["brems_activation_threshold_keV"]),
            inventory,
        )

    def _model_content_bytes(self) -> bytes:
        (
            config_hash,
            material_name,
            geometry_tier,
            physical_min,
            physical_max,
            brems_threshold,
            inventory,
        ) = self._binary_identity_fields()
        writer = io.BytesIO()
        self._put_string(writer, config_hash)
        self._put_string(writer, material_name)
        self._put_string(writer, geometry_tier)
        writer.write(struct.pack("<I", PAYLOAD_NORMALIZATION_S_HAT))
        writer.write(struct.pack("<d", brems_threshold))
        writer.write(struct.pack("<4d", *physical_min))
        writer.write(struct.pack("<4d", *physical_max))
        writer.write(struct.pack("<Q", len(inventory)))
        for collection_id, name, axis_kind in inventory:
            writer.write(struct.pack("<B", collection_id))
            self._put_string(writer, name)
            writer.write(struct.pack("<I", _AXIS_KIND_TO_ID[axis_kind]))
        writer.write(struct.pack("<Q", len(self.components)))
        for _, model in sorted(self.components.items()):
            self._put_model(writer, model)
        writer.write(struct.pack("<Q", len(self.lines)))
        for _, model in sorted(self.lines.items()):
            self._put_model(writer, model)
        return writer.getvalue()

    def compute_model_content_sha256(self) -> str:
        """Hash all runtime model state while excluding validation/publication state."""

        return hashlib.sha256(self._model_content_bytes()).hexdigest()

    def evaluate(
        self,
        energy_keV: float,
        distance_m: float,
        source_height_m: float,
        detector_height_m: float,
        target_energy_edges_keV: Sequence[float],
        target_mu_edges: Sequence[float],
        *,
        direct_fluence_cm2: float,
    ) -> dict[str, Any]:
        numpy = require_numpy()
        if not math.isfinite(direct_fluence_cm2) or direct_fluence_cm2 <= 0.0:
            raise ValueError("runtime evaluation requires positive analytic direct_fluence_cm2")
        node = PhysicalNode(energy_keV, distance_m, source_height_m, detector_height_m)
        node.validate()
        physical_domain = self.metadata.get("physical_domain")
        if physical_domain is not None:
            physical = (energy_keV, distance_m, source_height_m, detector_height_m)
            low = physical_domain["min"]
            high = physical_domain["max"]
            if any(value < lower or value > upper for value, lower, upper in zip(physical, low, high)):
                raise ValueError("query lies outside runtime payload physical domain")
        q = numpy.asarray([geometry_transform(node)], dtype=float)
        total = [
            [0.0 for _ in range(len(target_mu_edges) - 1)]
            for _ in range(len(target_energy_edges_keV) - 1)
        ]
        components_s_hat: dict[str, tuple[tuple[float, ...], ...]] = {}
        for name, model in self.components.items():
            if model.pair_threshold_keV is not None and energy_keV < model.pair_threshold_keV:
                rebinned = tuple(tuple(0.0 for _ in total[0]) for _ in total)
            else:
                registered_edges = registered_edges_from_physical(
                    model.axis_kind,
                    target_energy_edges_keV,
                    energy_keV,
                    e_cut_keV=model.e_cut_keV,
                )
                _, _, flat = model.predict_cell_masses(
                    q,
                    registered_edges,
                    target_mu_edges,
                )
                rebinned = tuple(
                    tuple(float(value) for value in row)
                    for row in flat[0].reshape(len(target_energy_edges_keV) - 1, -1)
                )
            components_s_hat[name] = rebinned
            for i, row in enumerate(rebinned):
                for j, value in enumerate(row):
                    total[i][j] += value

        exact_lines: list[dict[str, Any]] = []
        for name, model in self.lines.items():
            if model.pair_threshold_keV is not None and energy_keV < model.pair_threshold_keV:
                continue
            _, _, line_mass = model.predict_cell_masses(
                q,
                (0.0, 1.0),
                target_mu_edges,
            )
            angular_rebinned = tuple(float(value) for value in line_mass[0])
            exact_lines.append(
                {
                    "name": name,
                    "energy_keV": energy_keV if model.axis_kind == "line_e0" else 510.998950,
                    "s_hat_mu_bin_masses": angular_rebinned,
                    "mu_bin_masses": tuple(
                        value * direct_fluence_cm2 for value in angular_rebinned
                    ),
                }
            )
        continuum_s_hat = tuple(tuple(row) for row in total)
        return {
            "continuum_s_hat": continuum_s_hat,
            "continuum": tuple(
                tuple(value * direct_fluence_cm2 for value in row)
                for row in continuum_s_hat
            ),
            "components_s_hat": components_s_hat,
            "components": {
                name: tuple(
                    tuple(value * direct_fluence_cm2 for value in row) for row in masses
                )
                for name, masses in components_s_hat.items()
            },
            "lines": exact_lines,
        }

    def _content_dict(self) -> dict[str, Any]:
        computed_model_hash = self.compute_model_content_sha256()
        if self.model_content_sha256 and self.model_content_sha256 != computed_model_hash:
            raise ValueError("stored model-content SHA-256 does not match runtime model state")
        return {
            "schema": self.schema,
            "accepted_for_publication": self.accepted_for_publication,
            "model_content_sha256": computed_model_hash,
            "metadata": self.metadata,
            "components": {key: value.to_dict() for key, value in sorted(self.components.items())},
            "lines": {key: value.to_dict() for key, value in sorted(self.lines.items())},
        }

    def to_dict(self) -> dict[str, Any]:
        content = self._content_dict()
        canonical = json.dumps(content, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        return {**content, "sha256": hashlib.sha256(canonical).hexdigest()}

    def save(self, path: str | Path) -> None:
        destination = Path(path)
        encoded = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False)
        temporary = destination.with_name(destination.name + ".tmp")
        temporary.write_text(encoded + "\n", encoding="utf-8")
        temporary.replace(destination)

    def to_binary_bytes(self) -> bytes:
        """Return the exact dependency-free CEEFIT-PAYLOAD-v1 file bytes.

        Publication preflight must measure the artifact that would actually be
        written, including its header, checksum, publication bit, and validation
        metadata.  Keeping that operation beside the encoder also prevents a
        parameter-count or JSON-size proxy from drifting away from the runtime
        representation.
        """

        (
            config_hash,
            material_name,
            geometry_tier,
            physical_min,
            physical_max,
            brems_threshold,
            inventory,
        ) = self._binary_identity_fields()
        computed_model_hash = self.compute_model_content_sha256()
        if self.model_content_sha256 and self.model_content_sha256 != computed_model_hash:
            raise ValueError("stored model-content SHA-256 does not match runtime model state")
        validation = self.metadata.get("validation_metadata", {})
        q_domain = self.metadata.get("q_domain")
        if (
            not isinstance(q_domain, Mapping)
            or len(q_domain.get("min", ())) != 4
            or len(q_domain.get("max", ())) != 4
        ):
            raise ValueError("binary payload metadata needs four-dimensional q_domain")
        expected_q_min = tuple(map(float, q_domain["min"]))
        expected_q_max = tuple(map(float, q_domain["max"]))
        if any(
            not math.isfinite(value) for value in (*expected_q_min, *expected_q_max)
        ) or any(high <= low for low, high in zip(expected_q_min, expected_q_max)):
            raise ValueError("binary payload q_domain is non-finite or empty")
        physical_q_energy_min = math.log(physical_min[0])
        physical_q_energy_max = math.log(physical_max[0])
        if (
            expected_q_min[0] > physical_q_energy_min + 1.0e-12
            or expected_q_max[0] < physical_q_energy_max - 1.0e-12
        ):
            raise ValueError(
                "binary payload q_domain does not cover the physical energy domain"
            )
        for model in (*self.components.values(), *self.lines.values()):
            model_q_max = tuple(low + span for low, span in zip(model.q_min, model.q_span))
            if any(
                not math.isclose(a, b, rel_tol=1.0e-14, abs_tol=1.0e-14)
                for a, b in zip(model.q_min, expected_q_min)
            ) or any(
                not math.isclose(a, b, rel_tol=1.0e-14, abs_tol=1.0e-14)
                for a, b in zip(model_q_max, expected_q_max)
            ):
                raise ValueError(f"model {model.component} does not cover the campaign q_domain")
        validation_json = json.dumps(
            validation, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        writer = io.BytesIO()
        self._put_string(writer, config_hash)
        self._put_string(writer, material_name)
        self._put_string(writer, geometry_tier)
        writer.write(struct.pack("<I", PAYLOAD_NORMALIZATION_S_HAT))
        writer.write(struct.pack("<B", self.accepted_for_publication))
        writer.write(struct.pack("<d", brems_threshold))
        self._put_string(writer, computed_model_hash)
        writer.write(struct.pack("<4d", *physical_min))
        writer.write(struct.pack("<4d", *physical_max))
        writer.write(struct.pack("<Q", len(inventory)))
        for collection_id, name, axis_kind in inventory:
            writer.write(struct.pack("<B", collection_id))
            self._put_string(writer, name)
            writer.write(struct.pack("<I", _AXIS_KIND_TO_ID[axis_kind]))
        self._put_string(writer, validation_json)
        writer.write(struct.pack("<Q", len(self.components)))
        for _, model in sorted(self.components.items()):
            self._put_model(writer, model)
        writer.write(struct.pack("<Q", len(self.lines)))
        for _, model in sorted(self.lines.items()):
            self._put_model(writer, model)
        payload = writer.getvalue()
        digest = hashlib.sha256(payload).digest()
        return PAYLOAD_HEADER.pack(
            PAYLOAD_MAGIC,
            PAYLOAD_VERSION,
            PAYLOAD_ENDIAN_TAG,
            len(payload),
            digest,
        ) + payload

    def save_binary(self, path: str | Path) -> None:
        """Atomically write the dependency-free CEEFIT-PAYLOAD-v1 format."""

        encoded = self.to_binary_bytes()
        destination = Path(path)
        temporary = destination.with_name(destination.name + ".tmp")
        with open(temporary, "wb") as output:
            output.write(encoded)
        temporary.replace(destination)

    @classmethod
    def load(
        cls, path: str | Path, *, allow_candidate: bool = False
    ) -> "RuntimePayload":
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        checksum = record.pop("sha256", None)
        canonical = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        if checksum != hashlib.sha256(canonical).hexdigest():
            raise ValueError("runtime payload checksum mismatch")
        if record.get("schema") != "CEEFIT-PAYLOAD-v1":
            raise ValueError(f"unsupported runtime payload schema {record.get('schema')!r}")
        components = {
            key: StructuredPositiveModel.from_dict(value)
            for key, value in record.get("components", {}).items()
        }
        lines = {
            key: StructuredPositiveModel.from_dict(value)
            for key, value in record.get("lines", {}).items()
        }
        accepted = record.get("accepted_for_publication")
        if type(accepted) is not bool:
            raise ValueError("runtime payload lacks an explicit publication flag")
        if not accepted and not allow_candidate:
            raise ValueError(
                "runtime payload is a candidate; pass allow_candidate=True for validation use"
            )
        result = cls(
            components=components,
            lines=lines,
            metadata=dict(record.get("metadata", {})),
            schema=record["schema"],
            accepted_for_publication=accepted,
            model_content_sha256=str(record.get("model_content_sha256", "")),
        )
        result._validate_runtime_contract()
        if result.compute_model_content_sha256() != result.model_content_sha256:
            raise ValueError("runtime payload model-content SHA-256 mismatch")
        return result

    @classmethod
    def load_binary(
        cls, path: str | Path, *, allow_candidate: bool = False
    ) -> "RuntimePayload":
        """Read and strictly validate CEEFIT-PAYLOAD-v1."""

        numpy = require_numpy()
        raw = Path(path).read_bytes()
        if len(raw) < PAYLOAD_HEADER.size:
            raise ValueError("runtime payload is shorter than its header")
        magic, version, endian, payload_size, expected_sha = PAYLOAD_HEADER.unpack_from(raw)
        if magic != PAYLOAD_MAGIC or version != PAYLOAD_VERSION or endian != PAYLOAD_ENDIAN_TAG:
            raise ValueError("unsupported runtime payload header")
        payload = raw[PAYLOAD_HEADER.size :]
        if len(payload) != payload_size:
            raise ValueError("runtime payload length mismatch")
        if hashlib.sha256(payload).digest() != expected_sha:
            raise ValueError("runtime payload checksum mismatch")
        reader = io.BytesIO(payload)

        def take(size: int) -> bytes:
            value = reader.read(size)
            if len(value) != size:
                raise ValueError("truncated runtime payload")
            return value

        def get_string() -> str:
            count = struct.unpack("<Q", take(8))[0]
            if count > len(payload) - reader.tell():
                raise ValueError("runtime payload string length is invalid")
            try:
                return take(count).decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError("runtime payload contains invalid UTF-8") from error

        def get_vector() -> Any:
            count = struct.unpack("<Q", take(8))[0]
            byte_count = count * 8
            if byte_count > len(payload) - reader.tell():
                raise ValueError("runtime payload vector length is invalid")
            return numpy.frombuffer(take(byte_count), dtype="<f8").astype(float, copy=True)

        def get_model() -> StructuredPositiveModel:
            name = get_string()
            axis_id, e_cut, has_threshold = struct.unpack("<IdB", take(13))
            threshold, rank, level = struct.unpack("<dII", take(16))
            if has_threshold not in (0, 1) or axis_id not in _AXIS_ID_TO_KIND:
                raise ValueError("runtime payload model flags are invalid")
            energy_knots = get_vector()
            mu_knots = get_vector()
            q_min = list(struct.unpack("<4d", take(32)))
            q_span = list(struct.unpack("<4d", take(32)))
            yields = get_vector()
            logits = get_vector()
            archetypes = get_vector()
            features = _geometry_feature_count(level)
            n_energy_basis = energy_knots.size - 4
            n_mu_basis = mu_knots.size - 4
            is_line = axis_id in (4, 5)
            coefficient_count = n_mu_basis if is_line else n_energy_basis * n_mu_basis
            if (
                rank == 0
                or level not in {1, 2, 3, 4}
                or n_mu_basis < 4
                or (not is_line and n_energy_basis < 4)
                or (is_line and energy_knots.size != 0)
            ):
                raise ValueError("runtime payload model dimensions are invalid")
            if yields.size != features or logits.size != features * rank or archetypes.size != rank * coefficient_count:
                raise ValueError("runtime payload model vector dimensions are invalid")
            if any(not math.isfinite(value) for value in (*q_min, *q_span)) or any(
                value <= 0.0 for value in q_span
            ):
                raise ValueError("runtime payload q domain is invalid")
            if numpy.any(~numpy.isfinite(archetypes)) or numpy.any(archetypes < 0.0):
                raise ValueError("runtime payload archetypes are invalid")
            shaped_archetypes = archetypes.reshape(
                (rank, n_mu_basis) if is_line else (rank, n_energy_basis, n_mu_basis)
            )
            if not numpy.allclose(
                numpy.sum(shaped_archetypes.reshape(rank, -1), axis=1),
                1.0,
                rtol=1.0e-10,
                atol=1.0e-12,
            ):
                raise ValueError("runtime payload archetypes are not unit masses")
            return StructuredPositiveModel(
                component=name,
                axis_kind=_AXIS_ID_TO_KIND[axis_id],
                energy_knots=energy_knots.tolist(),
                mu_knots=mu_knots.tolist(),
                q_min=q_min,
                q_span=q_span,
                rank=rank,
                level=level,
                roughness_penalty=0.0,
                e_cut_keV=e_cut,
                pair_threshold_keV=threshold if has_threshold else None,
                log_yield_coefficients=yields.reshape(features, 1),
                logit_coefficients=logits.reshape(features, rank),
                archetype_coefficients=shaped_archetypes,
            )

        config_hash = get_string()
        material_name = get_string()
        geometry_tier = get_string()
        normalization_id = struct.unpack("<I", take(4))[0]
        if normalization_id != PAYLOAD_NORMALIZATION_S_HAT:
            raise ValueError("unsupported runtime payload normalization")
        publication_flag = struct.unpack("<B", take(1))[0]
        if publication_flag not in (0, 1):
            raise ValueError("runtime payload publication flag is invalid")
        brems_threshold = struct.unpack("<d", take(8))[0]
        model_content_sha256 = get_string()
        physical_min = list(struct.unpack("<4d", take(32)))
        physical_max = list(struct.unpack("<4d", take(32)))
        inventory_count = struct.unpack("<Q", take(8))[0]
        if inventory_count > len(_MODEL_NAME_TO_AXIS):
            raise ValueError("runtime payload model inventory is too large")
        inventory: dict[str, list[dict[str, str]]] = {
            "continuum": [],
            "lines": [],
        }
        for _ in range(inventory_count):
            collection_id = struct.unpack("<B", take(1))[0]
            name = get_string()
            axis_id = struct.unpack("<I", take(4))[0]
            if collection_id not in (0, 1) or axis_id not in _AXIS_ID_TO_KIND:
                raise ValueError("runtime payload model inventory entry is invalid")
            inventory[("continuum", "lines")[collection_id]].append(
                {"name": name, "axis_kind": _AXIS_ID_TO_KIND[axis_id]}
            )
        validation_json = get_string()
        try:
            validation = json.loads(validation_json)
        except json.JSONDecodeError as error:
            raise ValueError("runtime validation metadata is invalid JSON") from error
        continuum_count = struct.unpack("<Q", take(8))[0]
        components: dict[str, StructuredPositiveModel] = {}
        for _ in range(continuum_count):
            model = get_model()
            if model.component in components:
                raise ValueError("duplicate continuum model name")
            components[model.component] = model
        line_count = struct.unpack("<Q", take(8))[0]
        lines: dict[str, StructuredPositiveModel] = {}
        for _ in range(line_count):
            model = get_model()
            if model.component in lines:
                raise ValueError("duplicate line model name")
            lines[model.component] = model
        if reader.read(1):
            raise ValueError("runtime payload has trailing bytes")
        all_models = [*components.values(), *lines.values()]
        q_metadata = (
            {
                "min": list(all_models[0].q_min),
                "max": [
                    low + span
                    for low, span in zip(all_models[0].q_min, all_models[0].q_span)
                ],
            }
            if all_models
            else None
        )
        result = cls(
            components,
            lines,
            {
                "config_sha256": config_hash,
                "material_name": material_name,
                "geometry_tier": geometry_tier,
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": {"min": physical_min, "max": physical_max},
                "q_domain": q_metadata,
                "brems_activation_threshold_keV": brems_threshold,
                "model_inventory": inventory,
                "validation_metadata": validation,
            },
            accepted_for_publication=bool(publication_flag),
            model_content_sha256=model_content_sha256,
        )
        result._validate_runtime_contract()
        if result.compute_model_content_sha256() != model_content_sha256:
            raise ValueError("runtime payload model-content SHA-256 mismatch")
        if not result.accepted_for_publication and not allow_candidate:
            raise ValueError(
                "runtime payload is a candidate; pass allow_candidate=True for validation use"
            )
        return result
