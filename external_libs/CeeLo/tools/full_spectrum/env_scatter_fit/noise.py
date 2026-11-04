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

"""Contribution-aware replica noise metrics for archived joint shapes.

The archive deliberately retains the native 384 x Nmu maps.  Requiring the
empirical L1/TV noise of all of those cells to fall below 0.25% is a different
thing from avoiding relative-error gates on tiny cells: the absolute sampling
noise from thousands of resolved cells still adds in L1.  History allocation
therefore uses a fixed, mass-conserving 12 x 8 sentinel grid.  Native-grid TV
remains available as a diagnostic and native masses remain the fit/acceptance
truth.
"""

from __future__ import annotations

import math
from typing import Any


JOINT_SHAPE_SENTINEL_ENERGY_BINS = 12
JOINT_SHAPE_SENTINEL_MU_BINS = 8

# A weak exact atom is controlled by its absolute effect on the complete
# result.  A material atom additionally receives the stricter standalone
# relative-yield and angular-shape gates.
EXACT_LINE_MATERIAL_FRACTION = 0.01
EXACT_LINE_SIGNIFICANCE_SIGMAS = 3.0
EXACT_LINE_ABSOLUTE_NOISE_FRACTION = 0.001
EXACT_LINE_RELATIVE_NOISE_FRACTION = 0.0025
EXACT_LINE_ABSOLUTE_ERROR_FRACTION = 0.001
EXACT_LINE_RELATIVE_ERROR_FRACTION = 0.01

MODEL_SELECTION_POLICY_SCHEMA = (
    "CEEFIT-CONTRIBUTION-AWARE-MODEL-SELECTION-v1"
)


def contribution_aware_model_selection_policy() -> dict[str, Any]:
    """Return the versioned acceptance semantics bound into fit reports.

    This is deliberately a function rather than a mutable module-level mapping:
    callers receive an independent JSON-compatible record that can be serialized
    verbatim and compared exactly when stale reports or payloads are loaded.
    """

    return {
        "schema": MODEL_SELECTION_POLICY_SCHEMA,
        "component_material_fraction_3sigma": 0.05,
        "exact_line_material_fraction_3sigma": EXACT_LINE_MATERIAL_FRACTION,
        "significance_sigmas": EXACT_LINE_SIGNIFICANCE_SIGMAS,
        "weak_absolute_yield_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
        "weak_absolute_shape_l1_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
        "elastic_e0": "combined-direct-plus-environmental-atom",
    }


def exact_line_materiality(
    line_yield: float,
    line_yield_se: float,
    total_yield: float,
) -> tuple[bool, float]:
    """Classify a line using a conservative three-sigma contribution bound."""

    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (line_yield, line_yield_se, total_yield)
    ) or total_yield <= 0.0:
        return True, math.inf
    upper_fraction = (
        line_yield + EXACT_LINE_SIGNIFICANCE_SIGMAS * line_yield_se
    ) / total_yield
    return upper_fraction >= EXACT_LINE_MATERIAL_FRACTION, upper_fraction


def direct_line_map(numpy: Any, direct_yield: float, mu_bins: int) -> Any:
    """Represent the analytic head-on direct atom on a binned mu axis."""

    if not math.isfinite(direct_yield) or direct_yield < 0.0 or mu_bins <= 0:
        raise ValueError("direct exact-line inputs are invalid")
    result = numpy.zeros(mu_bins, dtype=float)
    result[-1] = direct_yield
    return result


def sentinel_shape(native_energy_bins: int, native_mu_bins: int) -> tuple[int, int]:
    """Return the fixed allocation grid supported by a native joint map."""

    if native_energy_bins <= 0 or native_mu_bins <= 0:
        raise ValueError("native joint-shape dimensions must be positive")
    energy_bins = min(JOINT_SHAPE_SENTINEL_ENERGY_BINS, native_energy_bins)
    mu_bins = min(JOINT_SHAPE_SENTINEL_MU_BINS, native_mu_bins)
    if native_energy_bins % energy_bins or native_mu_bins % mu_bins:
        raise ValueError(
            "native joint-shape dimensions must divide the 12x8 sentinel grid"
        )
    return energy_bins, mu_bins


def rebin_joint_masses(
    numpy: Any,
    values: Any,
    *,
    native_shape: tuple[int, int],
    target_shape: tuple[int, int] | None = None,
) -> Any:
    """Mass-conservingly sum the last, flattened joint-map dimension."""

    native_energy, native_mu = map(int, native_shape)
    target_energy, target_mu = (
        sentinel_shape(native_energy, native_mu)
        if target_shape is None
        else tuple(map(int, target_shape))
    )
    if (
        target_energy <= 0
        or target_mu <= 0
        or native_energy % target_energy
        or native_mu % target_mu
    ):
        raise ValueError("target joint-shape dimensions must divide the native grid")
    array = numpy.asarray(values, dtype=float)
    expected = native_energy * native_mu
    if array.ndim < 1 or array.shape[-1] != expected:
        raise ValueError(
            f"joint masses must have {expected} cells on their last dimension"
        )
    energy_factor = native_energy // target_energy
    mu_factor = native_mu // target_mu
    reshaped = array.reshape(
        *array.shape[:-1],
        target_energy,
        energy_factor,
        target_mu,
        mu_factor,
    )
    return reshaped.sum(axis=(-3, -1)).reshape(
        *array.shape[:-1], target_energy * target_mu
    )


def replica_shape_tv_noise(
    numpy: Any,
    replica_s1: Any,
    replica_histories: Any,
    truth: Any,
    *,
    native_shape: tuple[int, int],
    target_shape: tuple[int, int] | None = None,
) -> float:
    """Estimate aggregate normalized-shape TV noise from logical replicas.

    Each logical replica is normalized independently, so yield fluctuations do
    not masquerade as shape noise.  The RMS replica-to-aggregate TV is scaled
    by ``sqrt(R)`` to estimate the noise of their aggregate.  Rebinning occurs
    by exact summation before normalization/TV; no cell is dropped or assigned
    a relative-precision requirement.
    """

    raw = numpy.asarray(replica_s1, dtype=float)
    if raw.ndim != 2:
        raise ValueError("replica maps must be a two-dimensional array")
    counts = numpy.asarray(replica_histories, dtype=float).reshape(-1)
    if raw.shape[0] != counts.size:
        raise ValueError("replica maps and history counts disagree")
    valid = counts > 0.0
    if int(numpy.count_nonzero(valid)) < 2:
        return math.inf

    truth_array = numpy.asarray(truth, dtype=float).reshape(-1)
    native_energy, native_mu = map(int, native_shape)
    if truth_array.size != native_energy * native_mu:
        raise ValueError("truth map does not match its declared native shape")
    estimates = raw[valid] / counts[valid, None]
    if target_shape is not None or sentinel_shape(native_energy, native_mu) != native_shape:
        estimates = rebin_joint_masses(
            numpy,
            estimates,
            native_shape=native_shape,
            target_shape=target_shape,
        )
        truth_array = rebin_joint_masses(
            numpy,
            truth_array,
            native_shape=native_shape,
            target_shape=target_shape,
        )

    truth_yield = float(numpy.sum(truth_array))
    if truth_yield <= 0.0:
        return 0.0 if not numpy.any(raw) else math.inf
    yields = numpy.sum(estimates, axis=1)
    positive = yields > 0.0
    if int(numpy.count_nonzero(positive)) < 2:
        return math.inf
    truth_shape = truth_array / truth_yield
    shapes = estimates[positive] / yields[positive, None]
    replica_tv = 0.5 * numpy.sum(
        numpy.abs(shapes - truth_shape[None, :]), axis=1
    )
    return float(
        math.sqrt(float(numpy.mean(replica_tv * replica_tv)) / shapes.shape[0])
    )
