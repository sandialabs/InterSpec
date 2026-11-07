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

"""Shared training-data and basis library for the charged-particle source fit.

Loads e-/e+ campaign cells, builds the mechanistic yield baselines (CSDA range,
radiation length, and 511-keV photon optical depths), the positive cubic
B-spline shape bases in log photon energy, and the scoring metrics. The fitter
itself is fit_charged_particle_source_structured.py; the detector-response
fitter reuses the B-spline helpers.
"""

from __future__ import annotations

import csv
import bank_io
import hashlib
import math
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


ELECTRON_MASS_KEV = 510.99895
BASIS_SIZE = 32
MINIMUM_PHOTON_KEV = 10.0
RANGE_FEATURE_SCALE = 50.0
# Options for every cross-section helper query: the evaluator's escape envelope
# uses the smooth effective-Z blend, not the exact element, and next to a K
# edge the two differ by a few percent.
HELPER_OPTIONS = ("--smooth-z-blend",)
# Scored rows: a component whose truth gives at least 25 counts in the bank's
# smallest cells (60k events). A fixed yield, so the scored population does not
# grow with a cell's statistics.
MINIMUM_SCORED_YIELD = 25.0 / 60000.0
SPECTRAL_COMPONENTS = (
    "electron_brems",
    "positron_brems",
    "positron_511_atom",
    "positron_annihilation_continuum",
    "positron_annihilation_in_flight_direct",
    "positron_annihilation_in_flight_continuum",
)
# Shield K-alpha, K-beta and L x-rays per particle: yield-only atoms whose
# bases (a + contained x b) come from the C++ charged_fluorescence_bases
# through the evaluator's --fluorescence-bases mode.
FLUORESCENCE_FAMILIES = ("k_alpha", "k_beta", "l")
FLUORESCENCE_COMPONENTS = tuple(
    f"{particle}_fluorescence_{family}"
    for particle in ("electron", "positron") for family in FLUORESCENCE_FAMILIES)
MODEL_COMPONENTS = (
    "electron_brems",
    "positron_brems",
    "positron_containment",
    "positron_511_atom",
    "positron_annihilation_continuum",
    "positron_annihilation_in_flight_direct",
    "positron_annihilation_in_flight_continuum",
) + FLUORESCENCE_COMPONENTS
ATOM_COMPONENTS = ("positron_511_atom",) + FLUORESCENCE_COMPONENTS
REPORT_COMPONENTS = MODEL_COMPONENTS + (
    "electron_total",
    "positron_total",
)
CONTINUOUS_COMPONENTS = tuple(
    name for name in SPECTRAL_COMPONENTS
    if name != "positron_511_atom")


@dataclass
class Cell:
    identity: str
    path: Path
    energy_keV: float
    z_value: int
    areal_density: float
    is_positron: bool
    events: float
    edges: np.ndarray
    values: dict[str, np.ndarray]
    variance: dict[str, np.ndarray]
    atom_yield: float
    atom_variance: float
    escaped_primary_yield: float
    escaped_primary_variance: float
    split: str
    csda_range: float = 0.0
    electron_csda_range: float = 0.0
    radiation_length: float = 0.0
    mu511_pe: float = 0.0
    mu511_cs: float = 0.0
    mu511_pp: float = 0.0
    photon_absorption_mu: np.ndarray | None = None
    # Fluorescence truth per family, and the bases' a and b per family.
    fluorescence_yield: dict = field(default_factory=dict)
    fluorescence_variance: dict = field(default_factory=dict)
    fluorescence_a: dict = field(default_factory=dict)
    fluorescence_b: dict = field(default_factory=dict)


def split_for(
    identity: str,
    z_value: int,
    holdout_z: set[int],
    validation_energies_keV: set[float] | None = None,
    energy_keV: float | None = None,
    validation_z: set[int] | None = None,
) -> str:
    if z_value in holdout_z:
        return "test"
    if validation_z and z_value in validation_z:
        return "validation"
    if (
        energy_keV is not None
        and validation_energies_keV
        and energy_keV in validation_energies_keV
    ):
        return "validation"
    digest = hashlib.sha256(identity.encode()).digest()
    return "validation" if digest[0] % 5 == 0 else "train"


def load_cell(
    path: Path,
    holdout_z: set[int],
    validation_energies_keV: set[float] | None = None,
    validation_z: set[int] | None = None,
) -> Cell:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty campaign cell: {path}")
    first = rows[0]
    source = first["source"].strip()
    if source not in {"e-", "e+"}:
        raise ValueError(f"unsupported charged source {source}: {path}")
    is_positron = source == "e+"
    boundaries = sorted({
        float(value)
        for row in rows
        for value in (row["bin_lo_keV"], row["bin_hi_keV"])
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])
    })
    edges = np.asarray(boundaries)
    values = {
        name: np.zeros(edges.size - 1)
        for name in CONTINUOUS_COMPONENTS}
    variance = {
        name: np.zeros(edges.size - 1)
        for name in CONTINUOUS_COMPONENTS}
    atom_yield = 0.0
    atom_variance = 0.0
    escaped_primary_yield = 0.0
    escaped_primary_variance = 0.0
    fluorescence_yield = dict.fromkeys(FLUORESCENCE_FAMILIES, 0.0)
    fluorescence_variance = dict.fromkeys(FLUORESCENCE_FAMILIES, 0.0)
    events = float(first["normalization_events"])
    for row in rows:
        value = float(row["sum_w"]) / events
        value_variance = float(row["sum_w2"]) / (events * events)
        if row["component"] == "escaped_primary":
            escaped_primary_yield += value
            escaped_primary_variance += value_variance
            continue
        if row["particle"] != "photon":
            continue
        name = row["component"]
        if name == "secondary_charged_brems":
            # Delta-electron descendants are part of the escaping
            # charged-particle brems measure. Their separate MC truth label is
            # diagnostic; retaining it as a separate surrogate did not improve
            # the observable spectrum enough to justify another component.
            name = (
                "positron_brems" if is_positron else "electron_brems")
        if name == "positron_511_atom":
            atom_yield += value
            atom_variance += value_variance
            continue
        if name.startswith("fluorescence_"):
            family = name.removeprefix("fluorescence_")
            if family in fluorescence_yield:
                fluorescence_yield[family] += value
                fluorescence_variance[family] += value_variance
            continue
        if name not in values:
            continue
        lo = float(row["bin_lo_keV"])
        index = int(np.searchsorted(edges, lo, side="right") - 1)
        if 0 <= index < values[name].size:
            values[name][index] += value
            variance[name][index] += value_variance
    identity = bank_io.cell_id(path)
    z_value = int(first["Z"])
    energy_keV = float(first["primary_energy_keV"])
    return Cell(
        identity=identity,
        path=path,
        energy_keV=energy_keV,
        z_value=z_value,
        areal_density=float(first["areal_density_g_cm2"]),
        is_positron=is_positron,
        events=events,
        edges=edges,
        values=values,
        variance=variance,
        atom_yield=atom_yield,
        atom_variance=atom_variance,
        escaped_primary_yield=escaped_primary_yield,
        escaped_primary_variance=escaped_primary_variance,
        split=split_for(
            identity, z_value, holdout_z,
            validation_energies_keV, energy_keV, validation_z),
        fluorescence_yield=fluorescence_yield,
        fluorescence_variance=fluorescence_variance,
    )


def attach_properties(cells: list[Cell], helper: Path) -> None:
    queries = []
    for cell in cells:
        queries.append(
            f"{cell.energy_keV:.17g} {cell.z_value} "
            f"{cell.areal_density:.17g}\n")
        queries.append(
            f"{ELECTRON_MASS_KEV:.17g} {cell.z_value} "
            f"{cell.areal_density:.17g}\n")
    completed = subprocess.run(
        [str(helper), *HELPER_OPTIONS], input="".join(queries), text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != 2 * len(cells):
        raise RuntimeError("cross-section helper returned the wrong row count")
    for index, cell in enumerate(cells):
        charged = rows[2 * index]
        photon511 = rows[2 * index + 1]
        cell.csda_range = float(charged[
            "positron_range_g_cm2"
            if cell.is_positron else "electron_range_g_cm2"])
        cell.electron_csda_range = float(
            charged["electron_range_g_cm2"])
        cell.radiation_length = float(
            charged["radiation_length_g_cm2"])
        cell.mu511_pe = float(photon511["mu_pe"])
        cell.mu511_cs = float(photon511["mu_cs"])
        cell.mu511_pp = float(photon511["mu_pp"])


def attach_fluorescence_bases(cells: list[Cell], evaluator: Path) -> None:
    """a and b of each fluorescence family, from the C++ evaluator's own
    charged_fluorescence_bases (never recomputed here)."""
    queries = "".join(
        f"{cell.energy_keV:.17g} {cell.z_value} {cell.areal_density:.17g} 0 "
        f"{int(cell.is_positron)}\n" for cell in cells)
    completed = subprocess.run(
        [str(evaluator), "--fluorescence-bases"], input=queries, text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(FLUORESCENCE_FAMILIES) * len(cells):
        raise RuntimeError("fluorescence-bases helper returned the wrong row count")
    for index, cell in enumerate(cells):
        block = rows[index * len(FLUORESCENCE_FAMILIES):
                     (index + 1) * len(FLUORESCENCE_FAMILIES)]
        cell.fluorescence_a = {row["family"]: float(row["a"]) for row in block}
        cell.fluorescence_b = {row["family"]: float(row["b"]) for row in block}


def fluorescence_family(name: str) -> str:
    return name.split("_fluorescence_", 1)[1]


def attach_photon_absorption(cells: list[Cell], helper: Path) -> None:
    signatures: dict[
        tuple[int, float, int, float], tuple[np.ndarray, list[int]]
    ] = {}
    for cell in cells:
        key = (
            cell.z_value, cell.energy_keV,
            cell.edges.size, float(cell.edges[-1]))
        if key in signatures:
            continue
        valid = np.flatnonzero(
            (cell.edges[:-1] < cell.energy_keV)
            & (cell.edges[1:] > MINIMUM_PHOTON_KEV))
        lo = np.maximum(
            cell.edges[valid], MINIMUM_PHOTON_KEV)
        hi = np.minimum(
            cell.edges[valid + 1], cell.energy_keV)
        signatures[key] = (np.sqrt(lo * hi), valid.tolist())
    queries = []
    ordered = []
    for key, (energies, indices) in signatures.items():
        z_value = key[0]
        for energy, index in zip(energies, indices):
            queries.append(f"{energy:.17g} {z_value} 0\n")
            ordered.append((key, index))
    completed = subprocess.run(
        [str(helper), *HELPER_OPTIONS], input="".join(queries), text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(ordered):
        raise RuntimeError(
            "cross-section helper returned the wrong absorption row count")
    absorption: dict[tuple[int, float, int, float], np.ndarray] = {
        key: np.zeros(key[2] - 1) for key in signatures}
    for (key, index), row in zip(ordered, rows):
        absorption[key][index] = (
            float(row["mu_pe"]) + float(row["mu_pp"]))
    for cell in cells:
        key = (
            cell.z_value, cell.energy_keV,
            cell.edges.size, float(cell.edges[-1]))
        cell.photon_absorption_mu = absorption[key]


def positron_gate_range_ratio(cell: Cell) -> float:
    """Radial depth in the empirical source-escape gate coordinate.

    CeeLo's validated gate uses an electron extrapolated range even for a
    positron exit.  Its lower detour clamp produces real, sharp containment
    landmarks in the MC bank as energy and Z vary.  Warping the containment
    interpolation by this coordinate aligns those landmarks instead of asking
    a smooth raw-Z RBF to straddle them.

    The smooth clamps match the Jet-compatible evaluator.  Their 0.002 width is
    narrow compared with the empirical model uncertainty while avoiding a
    derivative discontinuity at a fitted effective Z.
    """
    energy_mev = cell.energy_keV * 1.0e-3
    energy_modulation = 1.0 / (
        1.0 + (energy_mev / 0.35) ** 1.6)
    raw_detour = 1.0 - 0.0645 * cell.z_value * energy_modulation
    smoothing = 0.002
    detour = 0.5 * (
        raw_detour + 0.10
        + math.sqrt((raw_detour - 0.10) ** 2 + smoothing ** 2))
    detour = 0.5 * (
        detour + 0.98
        - math.sqrt((detour - 0.98) ** 2 + smoothing ** 2))
    extrapolated_range = max(
        cell.electron_csda_range * detour, 1.0e-300)
    return cell.areal_density / extrapolated_range


def cell_feature(cell: Cell, component: str | None = None) -> np.ndarray:
    range_ratio = cell.areal_density / cell.csda_range
    if component == "positron_containment":
        range_ratio = positron_gate_range_ratio(cell)
    return np.asarray([
        math.log(cell.energy_keV / 10.0) / math.log(1200.0),
        cell.z_value / 92.0,
        math.log1p(range_ratio) / math.log1p(RANGE_FEATURE_SCALE),
    ])


ELECTRON_COMPONENTS = {"electron_brems"} | {
    name for name in FLUORESCENCE_COMPONENTS if name.startswith("electron_")}


def component_applies(cell: Cell, name: str) -> bool:
    return (name in ELECTRON_COMPONENTS) != cell.is_positron


def positron_containment_truth(cell: Cell) -> float:
    return min(1.0, max(0.0, 1.0 - cell.escaped_primary_yield))


def positron_containment_turn_on(cell: Cell) -> float:
    range_ratio = cell.areal_density / cell.csda_range
    return -math.expm1(-100.0 * range_ratio)


def annihilation_escape_averages(cell: Cell) -> tuple[float, float]:
    range_ratio = cell.areal_density / cell.csda_range
    radius_fraction = math.exp(
        -0.125 * math.log1p(range_ratio ** 8))
    nodes, weights = np.polynomial.legendre.leggauss(8)
    tau_absorb = cell.areal_density * (
        cell.mu511_pe + cell.mu511_pp)
    tau_compton = cell.areal_density * cell.mu511_cs
    atom_escape = 0.0
    compton_escape = 0.0
    for cosine, weight in zip(nodes, weights):
        transverse2 = (
            radius_fraction * radius_fraction
            * (1.0 - cosine * cosine))
        path_fraction = (
            -radius_fraction * cosine
            + math.sqrt(max(0.0, 1.0 - transverse2)))
        survives_absorption = math.exp(
            -tau_absorb * path_fraction)
        survives_without_compton = math.exp(
            -tau_compton * path_fraction)
        atom_escape += (
            0.5 * weight * survives_absorption
            * survives_without_compton)
        compton_escape += (
            0.5 * weight * survives_absorption
            * (1.0 - survives_without_compton))
    return atom_escape, compton_escape


def component_base(
    cell: Cell,
    name: str,
    containment_probability: float | None = None,
) -> float:
    range_ratio = cell.areal_density / cell.csda_range
    if name in {"electron_brems", "positron_brems"}:
        effective_path = (
            cell.csda_range * (-math.expm1(-range_ratio)))
        return (
            effective_path / cell.radiation_length
            * math.log(cell.energy_keV / MINIMUM_PHOTON_KEV))
    if name == "positron_containment":
        # Keep the exact zero-thickness constraint in the mechanistic base.
        # The fitted conditional residual must divide this factor out; applying
        # it only after fitting biases every finite-thickness training node.
        return positron_containment_turn_on(cell)
    contained = (
        positron_containment_truth(cell)
        if containment_probability is None
        else containment_probability)
    if name in FLUORESCENCE_COMPONENTS:
        family = fluorescence_family(name)
        b = cell.fluorescence_b.get(family, 0.0) if cell.is_positron else 0.0
        return cell.fluorescence_a.get(family, 0.0) + contained * b
    atom_escape, compton_escape = annihilation_escape_averages(cell)
    if name == "positron_511_atom":
        return 2.0 * contained * atom_escape
    if name == "positron_annihilation_continuum":
        return 2.0 * contained * compton_escape
    if name in {
        "positron_annihilation_in_flight_direct",
        "positron_annihilation_in_flight_continuum",
    }:
        energy_onset = 1.0 - math.exp(
            -cell.energy_keV / ELECTRON_MASS_KEV)
        return (
            2.0 * energy_onset * (-math.expm1(-range_ratio)))
    raise KeyError(name)


def component_support(cell: Cell, name: str) -> float:
    if name == "positron_annihilation_continuum":
        return ELECTRON_MASS_KEV
    if name in {
        "positron_annihilation_in_flight_direct",
        "positron_annihilation_in_flight_continuum",
    }:
        return cell.energy_keV + ELECTRON_MASS_KEV
    return cell.energy_keV


def truth_for(
    cell: Cell, name: str,
) -> tuple[np.ndarray, np.ndarray]:
    if name == "positron_containment":
        return (
            np.asarray([positron_containment_truth(cell)]),
            np.asarray([cell.escaped_primary_variance]))
    if name == "positron_511_atom":
        return (
            np.asarray([cell.atom_yield]),
            np.asarray([cell.atom_variance]))
    if name in FLUORESCENCE_COMPONENTS:
        family = fluorescence_family(name)
        return (
            np.asarray([cell.fluorescence_yield.get(family, 0.0)]),
            np.asarray([cell.fluorescence_variance.get(family, 0.0)]))
    return cell.values[name], cell.variance[name]


def cubic_bspline_values(coordinate: float) -> np.ndarray:
    degree = 3
    knot_count = BASIS_SIZE + degree + 1
    knots = np.zeros(knot_count)
    for index in range(knot_count):
        if index <= degree:
            knots[index] = 0.0
        elif index >= BASIS_SIZE:
            knots[index] = 1.0
        else:
            knots[index] = (
                (index - degree) / (BASIS_SIZE - degree))
    answer = np.zeros(BASIS_SIZE)
    if coordinate >= 1.0:
        answer[-1] = 1.0
        return answer
    work_size = BASIS_SIZE + degree
    previous = np.zeros(work_size)
    for index in range(knot_count - 1):
        if knots[index] <= coordinate < knots[index + 1]:
            previous[index] = 1.0
    for order in range(1, degree + 1):
        current = np.zeros(work_size)
        count = work_size - order
        for index in range(count):
            left_denominator = (
                knots[index + order] - knots[index])
            if left_denominator > 0.0:
                current[index] += (
                    previous[index]
                    * (coordinate - knots[index])
                    / left_denominator)
            right_denominator = (
                knots[index + order + 1] - knots[index + 1])
            if right_denominator > 0.0:
                current[index] += (
                    previous[index + 1]
                    * (knots[index + order + 1] - coordinate)
                    / right_denominator)
        previous = current
    return previous[:BASIS_SIZE]


def cubic_bspline_integrals(
    coordinate_lo: float, coordinate_hi: float,
) -> np.ndarray:
    answer = np.zeros(BASIS_SIZE)
    if not coordinate_hi > coordinate_lo:
        return answer
    degree = 3
    intervals = BASIS_SIZE - degree
    inverse_sqrt_three = 0.5773502691896258
    segment_lo = coordinate_lo
    for boundary in range(1, intervals + 1):
        knot = boundary / intervals
        segment_hi = min(coordinate_hi, knot)
        if segment_hi > segment_lo:
            midpoint = 0.5 * (segment_lo + segment_hi)
            half_width = 0.5 * (segment_hi - segment_lo)
            answer += half_width * (
                cubic_bspline_values(
                    midpoint - half_width * inverse_sqrt_three)
                + cubic_bspline_values(
                    midpoint + half_width * inverse_sqrt_three))
        if knot >= coordinate_hi:
            break
        segment_lo = max(segment_lo, knot)
    for index in range(BASIS_SIZE):
        knot_lo = index
        knot_hi = index + degree + 1

        def knot_value(knot_index: int) -> float:
            if knot_index <= degree:
                return 0.0
            if knot_index >= BASIS_SIZE:
                return 1.0
            return (knot_index - degree) / intervals

        normalization = (
            (knot_value(knot_hi) - knot_value(knot_lo))
            / (degree + 1))
        answer[index] /= normalization
    return answer


def bspline_matrix_from_coordinates(
    coordinate_lo: np.ndarray,
    coordinate_hi: np.ndarray,
    valid: np.ndarray,
) -> np.ndarray:
    answer = np.zeros((coordinate_lo.size, BASIS_SIZE))
    for index in np.flatnonzero(valid):
        answer[index] = cubic_bspline_integrals(
            coordinate_lo[index], coordinate_hi[index])
    column_normalization = answer.sum(axis=0)
    populated = column_normalization > 0.0
    answer[:, populated] /= column_normalization[populated]
    return answer


def log_bspline_matrix(
    output_edges: np.ndarray,
    maximum_energy_keV: float,
) -> np.ndarray:

    valid = (
        (output_edges[:-1] < maximum_energy_keV)
        & (output_edges[1:] > MINIMUM_PHOTON_KEV))
    lo_energy = np.maximum(
        output_edges[:-1], MINIMUM_PHOTON_KEV)
    hi_energy = np.minimum(
        output_edges[1:], maximum_energy_keV)
    inverse_span = 1.0 / math.log(
        maximum_energy_keV / MINIMUM_PHOTON_KEV)
    lo = np.clip(
        np.log(lo_energy / MINIMUM_PHOTON_KEV) * inverse_span,
        0.0, 1.0)
    hi = np.clip(
        np.log(hi_energy / MINIMUM_PHOTON_KEV) * inverse_span,
        0.0, 1.0)
    return bspline_matrix_from_coordinates(lo, hi, valid)


def in_flight_direct_bspline_matrix(
    output_edges: np.ndarray,
    kinetic_energy_keV: float,
) -> np.ndarray:
    minimum = ELECTRON_MASS_KEV
    maximum = kinetic_energy_keV + ELECTRON_MASS_KEV
    valid = (
        (output_edges[:-1] < maximum)
        & (output_edges[1:] > minimum))
    lo = np.clip(
        (np.maximum(output_edges[:-1], minimum) - minimum)
        / kinetic_energy_keV, 0.0, 1.0)
    hi = np.clip(
        (np.minimum(output_edges[1:], maximum) - minimum)
        / kinetic_energy_keV, 0.0, 1.0)
    return bspline_matrix_from_coordinates(lo, hi, valid)


_UNATTENUATED_BASIS_CACHE: dict[
    tuple[str, float, int, float], np.ndarray] = {}


def component_basis(
    cell: Cell, name: str,
    photon_absorption_path_fraction: float = 0.0,
) -> np.ndarray:
    key = (
        name, cell.energy_keV, cell.edges.size,
        float(cell.edges[-1]))
    cached = _UNATTENUATED_BASIS_CACHE.get(key)
    if cached is not None:
        basis = cached
    elif name == "positron_annihilation_in_flight_direct":
        basis = in_flight_direct_bspline_matrix(
            cell.edges, cell.energy_keV)
        _UNATTENUATED_BASIS_CACHE[key] = basis
    else:
        basis = log_bspline_matrix(
            cell.edges, component_support(cell, name))
        _UNATTENUATED_BASIS_CACHE[key] = basis
    if (
        name not in {"electron_brems", "positron_brems"}
        or not photon_absorption_path_fraction > 0.0
    ):
        return basis
    if cell.photon_absorption_mu is None:
        raise RuntimeError("photon absorption cross sections are missing")
    envelope = np.exp(
        -photon_absorption_path_fraction
        * cell.areal_density * cell.photon_absorption_mu)
    attenuated = basis * envelope[:, None]
    normalization = attenuated.sum(axis=0)
    populated = normalization > 1e-300
    attenuated[:, populated] /= normalization[populated]
    return attenuated


def least_squares_solve(matrix: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Least-squares solve that tolerates a non-converging LAPACK SVD.

    Apple's Accelerate LAPACK occasionally reports "SVD did not converge" from
    the divide-and-conquer driver even on well-conditioned inputs. Fall back to
    the normal equations with a relative ridge far below the fitted signal,
    which yields the same solution to well past fit precision.
    """
    try:
        answer, *_ = np.linalg.lstsq(matrix, target, rcond=None)
        if np.all(np.isfinite(answer)):
            return answer
    except np.linalg.LinAlgError:
        pass
    gram = matrix.T @ matrix
    scale = float(np.trace(gram)) / max(gram.shape[0], 1)
    ridge = 1.0e-12 * (scale if scale > 0.0 else 1.0)
    return np.linalg.solve(
        gram + ridge * np.eye(gram.shape[0]), matrix.T @ target)


def nonnegative_least_squares(
    matrix: np.ndarray, target: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Small Lawson-Hanson NNLS implementation for the 32-term basis."""
    columns = matrix.shape[1]
    solution = np.zeros(columns)
    passive = np.zeros(columns, dtype=bool)
    tolerance = (
        10.0 * np.finfo(float).eps
        * np.linalg.norm(matrix, ord=1)
        * max(np.linalg.norm(target), 1.0))
    gradient = matrix.T @ (target - matrix @ solution)
    iterations = 0
    while (
        np.any((~passive) & (gradient > tolerance))
        and iterations < 5 * columns
    ):
        masked_gradient = np.where(~passive, gradient, -np.inf)
        passive[int(np.argmax(masked_gradient))] = True
        while True:
            candidate = np.zeros(columns)
            candidate[passive] = least_squares_solve(
                matrix[:, passive], target)
            if np.all(candidate[passive] > tolerance):
                solution = candidate
                break
            offending = passive & (candidate <= tolerance)
            denominators = solution[offending] - candidate[offending]
            valid = denominators > 0.0
            alpha = (
                np.min(solution[offending][valid] / denominators[valid])
                if np.any(valid) else 0.0)
            solution += alpha * (candidate - solution)
            passive[passive & (solution <= tolerance)] = False
            solution[solution <= tolerance] = 0.0
        gradient = matrix.T @ (target - matrix @ solution)
        iterations += 1
    residual = float(np.linalg.norm(matrix @ solution - target))
    return solution, residual


def grouped_metrics(
    truth: np.ndarray,
    variance: np.ndarray,
    prediction: np.ndarray,
    groups: int = 32,
) -> dict[str, float]:
    usable_groups = min(groups, truth.size)
    boundaries = np.linspace(
        0, truth.size, usable_groups + 1, dtype=int)

    def combine(value: np.ndarray) -> np.ndarray:
        return np.asarray([
            value[boundaries[index]:boundaries[index + 1]].sum()
            for index in range(usable_groups)])

    truth_grouped = combine(truth)
    variance_grouped = combine(variance)
    prediction_grouped = combine(prediction)
    truth_yield = float(truth_grouped.sum())
    predicted_yield = float(prediction_grouped.sum())
    if truth_yield <= 0.0:
        return {
            "truth_yield": truth_yield,
            "predicted_yield": predicted_yield,
            "yield_relative_uncertainty": math.inf,
            "yield_relative_error": math.inf,
            "shape_total_variation": math.inf,
            "weighted_rms_relative_error": math.inf,
            "p95_relative_error": math.inf,
        }
    uncertainty = (
        np.sqrt(np.maximum(variance_grouped, 0.0))
        / np.maximum(truth_grouped, 1e-300))
    valid = (truth_grouped > 0.0) & (uncertainty <= 0.05)
    relative = (
        np.abs(prediction_grouped - truth_grouped)
        / np.maximum(truth_grouped, 1e-300))
    if np.any(valid):
        weighted_rms = math.sqrt(float(
            np.sum(
                truth_grouped[valid] * relative[valid] ** 2)
            / np.sum(truth_grouped[valid])))
        p95 = float(np.percentile(relative[valid], 95))
    else:
        weighted_rms = math.inf
        p95 = math.inf
    truth_shape = truth_grouped / truth_yield
    prediction_shape = (
        prediction_grouped / max(predicted_yield, 1e-300))
    return {
        "truth_yield": truth_yield,
        "predicted_yield": predicted_yield,
        "yield_relative_uncertainty": (
            math.sqrt(float(variance_grouped.sum())) / truth_yield),
        "yield_relative_error": (
            predicted_yield - truth_yield) / truth_yield,
        "shape_total_variation": 0.5 * float(
            np.abs(truth_shape - prediction_shape).sum()),
        "weighted_rms_relative_error": weighted_rms,
        "p95_relative_error": p95,
    }


def evaluate(
    cells: list[Cell],
    predict,
    basis_cache: dict[tuple, np.ndarray] | None = None,
    photon_absorption_path_fraction: float = 0.0,
) -> tuple[list[dict], dict[tuple[str, str], np.ndarray]]:
    """Score every non-training cell against its truth.

    ``predict(cell, name, containment_probability, basis_cache,
    photon_absorption_path_fraction)`` returns one component's prediction:
    a one-element array for the yield-only components, otherwise the binned
    spectrum on the cell's edges.
    """
    records = []
    predictions = {}
    for cell in cells:
        if cell.split == "train":
            continue
        containment_probability = None
        if cell.is_positron:
            containment_prediction = predict(
                cell, "positron_containment", None, basis_cache,
                photon_absorption_path_fraction)
            containment_probability = float(
                containment_prediction[0])
            predictions[
                (cell.identity, "positron_containment")
            ] = containment_prediction
        for name in MODEL_COMPONENTS:
            if not component_applies(cell, name):
                continue
            truth, variance = truth_for(cell, name)
            if name == "positron_containment":
                prediction = predictions[(cell.identity, name)]
            else:
                prediction = predict(
                    cell, name, containment_probability, basis_cache,
                    photon_absorption_path_fraction)
                predictions[(cell.identity, name)] = prediction
            if float(truth.sum()) < MINIMUM_SCORED_YIELD:
                continue
            records.append({
                "identity": cell.identity,
                "split": cell.split,
                "component": name,
                "energy_keV": cell.energy_keV,
                "Z": cell.z_value,
                "areal_density_g_cm2": cell.areal_density,
                "range_ratio": cell.areal_density / cell.csda_range,
                **grouped_metrics(truth, variance, prediction),
            })
        total_truth = np.zeros(cell.edges.size - 1)
        total_variance = np.zeros(cell.edges.size - 1)
        total_prediction = np.zeros(cell.edges.size - 1)
        total_name = (
            "positron_total" if cell.is_positron else "electron_total")
        continuous_names = (
            tuple(
                name for name in CONTINUOUS_COMPONENTS
                if component_applies(cell, name)))
        for name in continuous_names:
            total_truth += cell.values[name]
            total_variance += cell.variance[name]
            total_prediction += predictions[(cell.identity, name)]
        if cell.is_positron:
            atom_bin = int(
                np.searchsorted(
                    cell.edges, ELECTRON_MASS_KEV, side="right") - 1)
            if 0 <= atom_bin < total_truth.size:
                total_truth[atom_bin] += cell.atom_yield
                total_variance[atom_bin] += cell.atom_variance
                total_prediction[atom_bin] += float(
                    predictions[
                        (cell.identity, "positron_511_atom")][0])
        predictions[(cell.identity, total_name)] = total_prediction
        if float(total_truth.sum()) >= MINIMUM_SCORED_YIELD:
            records.append({
                "identity": cell.identity,
                "split": cell.split,
                "component": total_name,
                "energy_keV": cell.energy_keV,
                "Z": cell.z_value,
                "areal_density_g_cm2": cell.areal_density,
                "range_ratio": cell.areal_density / cell.csda_range,
                **grouped_metrics(
                    total_truth, total_variance, total_prediction),
            })
    return records, predictions


def finite_summary(values: list[float]) -> dict[str, float | int]:
    finite = np.asarray([
        value for value in values if math.isfinite(value)])
    if finite.size == 0:
        return {"count": 0}
    return {
        "count": int(finite.size),
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95)),
        "maximum": float(np.max(finite)),
    }


def write_metrics(path: Path, records: list[dict]) -> None:
    if not records:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


