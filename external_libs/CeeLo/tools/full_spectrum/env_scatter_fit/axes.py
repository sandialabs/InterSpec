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

"""Registered component axes and conservative histogram rebinning."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


ELECTRON_REST_KEV = 510.998950
ANNIHILATION_BACKSCATTER_KEV = ELECTRON_REST_KEV / 3.0


@dataclass(frozen=True)
class AxisDefinition:
    component: str
    physical_energy_edges_keV: tuple[float, ...]
    registered_edges: tuple[float, ...]
    mu_edges: tuple[float, ...]

    @property
    def shape(self) -> tuple[int, int]:
        return (len(self.physical_energy_edges_keV) - 1, len(self.mu_edges) - 1)


def _log_edges(low: float, high: float, bins: int) -> list[float]:
    if not (low > 0.0 and high > low and bins > 0):
        raise ValueError("log edges require 0 < low < high and bins > 0")
    ratio = high / low
    result = [low * ratio ** (i / bins) for i in range(bins + 1)]
    result[0], result[-1] = low, high
    return result


def uniform_mu_edges(n_mu: int = 64) -> tuple[float, ...]:
    if n_mu <= 0:
        raise ValueError("n_mu must be positive")
    edges = [-1.0 + 2.0 * i / n_mu for i in range(n_mu + 1)]
    edges[0], edges[-1] = -1.0, 1.0
    return tuple(edges)


def _two_block_axis(
    component: str,
    low: float,
    landmark: float,
    high: float,
    n_per_block: int,
    n_mu: int,
) -> AxisDefinition:
    if not (low < landmark < high):
        raise ValueError("registered landmark must lie strictly inside the energy range")
    first = _log_edges(low, landmark, n_per_block)
    second = _log_edges(landmark, high, n_per_block)
    physical = tuple(first + second[1:])
    registered = tuple(
        [0.5 * i / n_per_block for i in range(n_per_block + 1)]
        + [0.5 + 0.5 * i / n_per_block for i in range(1, n_per_block + 1)]
    )
    return AxisDefinition(component, physical, registered, uniform_mu_edges(n_mu))


def primary_registered_axis(
    energy_keV: float,
    *,
    e_cut_keV: float = 10.0,
    n_per_block: int = 192,
    n_mu: int = 64,
) -> AxisDefinition:
    backscatter = energy_keV / (1.0 + 2.0 * energy_keV / ELECTRON_REST_KEV)
    return _two_block_axis(
        "primary", e_cut_keV, backscatter, energy_keV, n_per_block, n_mu
    )


def annihilation_registered_axis(
    *, e_cut_keV: float = 10.0, n_per_block: int = 192, n_mu: int = 64
) -> AxisDefinition:
    return _two_block_axis(
        "annihilation",
        e_cut_keV,
        ANNIHILATION_BACKSCATTER_KEV,
        ELECTRON_REST_KEV,
        n_per_block,
        n_mu,
    )


def brems_registered_axis(
    energy_keV: float,
    *,
    e_cut_keV: float = 10.0,
    n_energy: int = 384,
    n_mu: int = 64,
) -> AxisDefinition:
    physical = tuple(_log_edges(e_cut_keV, energy_keV, n_energy))
    registered = tuple(i / n_energy for i in range(n_energy + 1))
    return AxisDefinition("brems", physical, registered, uniform_mu_edges(n_mu))


def physical_edges_from_registered(
    component: str,
    registered_edges: Sequence[float],
    energy_keV: float,
    *,
    e_cut_keV: float = 10.0,
) -> tuple[float, ...]:
    """Map canonical registered edges back into physical energy."""

    if component == "primary":
        landmark, high = (
            energy_keV / (1.0 + 2.0 * energy_keV / ELECTRON_REST_KEV),
            energy_keV,
        )
    elif component == "annihilation":
        landmark, high = ANNIHILATION_BACKSCATTER_KEV, ELECTRON_REST_KEV
    elif component == "brems":
        landmark, high = None, energy_keV
    else:
        raise ValueError(f"unknown component axis {component!r}")
    result: list[float] = []
    for u in registered_edges:
        if not 0.0 <= u <= 1.0:
            raise ValueError("registered edges must lie in [0, 1]")
        if component == "brems":
            value = e_cut_keV * (high / e_cut_keV) ** u
        elif u <= 0.5:
            value = e_cut_keV * (landmark / e_cut_keV) ** (2.0 * u)  # type: ignore[operator]
        else:
            value = landmark * (high / landmark) ** (2.0 * u - 1.0)  # type: ignore[operator]
        result.append(value)
    if result:
        result[0] = e_cut_keV if registered_edges[0] == 0.0 else result[0]
        result[-1] = high if registered_edges[-1] == 1.0 else result[-1]
    return tuple(result)


def registered_edges_from_physical(
    component: str,
    physical_edges_keV: Sequence[float],
    energy_keV: float,
    *,
    e_cut_keV: float = 10.0,
) -> tuple[float, ...]:
    """Map physical energy edges into a component coordinate (unclipped)."""

    if component == "primary":
        landmark, high = (
            energy_keV / (1.0 + 2.0 * energy_keV / ELECTRON_REST_KEV),
            energy_keV,
        )
    elif component == "annihilation":
        landmark, high = ANNIHILATION_BACKSCATTER_KEV, ELECTRON_REST_KEV
    elif component == "brems":
        landmark, high = None, energy_keV
    else:
        raise ValueError(f"unknown component axis {component!r}")
    result = []
    for physical in physical_edges_keV:
        if physical <= 0.0:
            result.append(-math.inf)
        elif component == "brems":
            result.append(math.log(physical / e_cut_keV) / math.log(high / e_cut_keV))
        elif physical <= landmark:  # type: ignore[operator]
            result.append(0.5 * math.log(physical / e_cut_keV) / math.log(landmark / e_cut_keV))  # type: ignore[operator]
        else:
            result.append(0.5 + 0.5 * math.log(physical / landmark) / math.log(high / landmark))  # type: ignore[operator]
    return tuple(result)


def _validate_edges(edges: Sequence[float], name: str) -> None:
    if len(edges) < 2:
        raise ValueError(f"{name} needs at least two edges")
    if not all(math.isfinite(x) for x in edges):
        raise ValueError(f"{name} contains a non-finite edge")
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError(f"{name} must be strictly increasing")


def overlap_matrix(
    source_edges: Sequence[float],
    target_edges: Sequence[float],
    *,
    require_full_coverage: bool = True,
) -> tuple[tuple[float, ...], ...]:
    """Fractions of each source bin assigned to each target bin."""

    _validate_edges(source_edges, "source_edges")
    _validate_edges(target_edges, "target_edges")
    if require_full_coverage and (
        target_edges[0] > source_edges[0] or target_edges[-1] < source_edges[-1]
    ):
        raise ValueError("target edges must cover the full source range")
    n_src, n_dst = len(source_edges) - 1, len(target_edges) - 1
    matrix = [[0.0] * n_src for _ in range(n_dst)]
    i = j = 0
    while i < n_src and j < n_dst:
        low = max(source_edges[i], target_edges[j])
        high = min(source_edges[i + 1], target_edges[j + 1])
        if high > low:
            matrix[j][i] = (high - low) / (source_edges[i + 1] - source_edges[i])
        if source_edges[i + 1] <= target_edges[j + 1]:
            i += 1
        else:
            j += 1
    return tuple(tuple(row) for row in matrix)


def conservative_rebin_2d(
    masses: Sequence[Sequence[float]],
    source_energy_edges: Sequence[float],
    source_mu_edges: Sequence[float],
    target_energy_edges: Sequence[float],
    target_mu_edges: Sequence[float],
    *,
    require_full_coverage: bool = True,
) -> tuple[tuple[float, ...], ...]:
    """Conservatively rebin cell-integrated masses in energy and angle.

    Uniform density is assumed only within each source cell.  The operation is
    exactly mass-conserving (up to floating-point summation) whenever the target
    axes cover the source axes.
    """

    n_e, n_mu = len(source_energy_edges) - 1, len(source_mu_edges) - 1
    if len(masses) != n_e or any(len(row) != n_mu for row in masses):
        raise ValueError("mass array shape does not match source axes")
    e_overlap = overlap_matrix(
        source_energy_edges,
        target_energy_edges,
        require_full_coverage=require_full_coverage,
    )
    mu_overlap = overlap_matrix(
        source_mu_edges,
        target_mu_edges,
        require_full_coverage=require_full_coverage,
    )
    result = [[0.0] * len(mu_overlap) for _ in range(len(e_overlap))]
    for de, e_row in enumerate(e_overlap):
        e_terms = [(se, fraction) for se, fraction in enumerate(e_row) if fraction]
        for dm, mu_row in enumerate(mu_overlap):
            mu_terms = [(sm, fraction) for sm, fraction in enumerate(mu_row) if fraction]
            result[de][dm] = math.fsum(
                masses[se][sm] * ef * mf
                for se, ef in e_terms
                for sm, mf in mu_terms
            )
    return tuple(tuple(row) for row in result)
