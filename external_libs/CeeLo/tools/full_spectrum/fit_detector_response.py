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

"""Fit a positive functional HPGe pulse-response model.

Each process yield is a smooth positive RBF field in incident energy and
angle.  Continuous process shapes are positive cubic B-spline expansions in a
moving, process-specific deposited-energy coordinate.  The serialized payload
contains smooth coefficients, not detector-response histograms.
"""

from __future__ import annotations

import argparse
import atexit
import csv
import bank_io
import json
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import fit_charged_particle_source as basis_helpers


BASIS_SIZE = 96
RATE_ENERGY_BASIS_SIZE = 32
RATE_ANGLE_BASIS_SIZE = 24
# The passive-proxy stretch of the energy coordinate buys sub-60 keV
# resolution when there is no surface-entry anchor. With the anchor carrying
# the turn-on, alpha = 0 measures better everywhere; un-anchored fits should
# pass alpha = 3.
RATE_ENERGY_WARP_ALPHA = 0.0
# Angle coordinate reshaping for the detected-rate surface. A plain polar
# angle (exponent 1) leaves the rate surface unable to follow how fast the
# detected rate turns over with angle: its error then oscillates between
# training angles, reaching +5% near end-on incidence. Passing the centred
# polar angle through x -> x*(x^2+eps^2)^((q-1)/2) concentrates resolution
# where it is needed; eps > 0 keeps the derivative at cos(theta) = 0 finite,
# which a plain |x|^q would not.
RATE_ANGLE_WARP_EXPONENT = 0.5
RATE_ANGLE_WARP_EPSILON = 0.03
RATE_SMOOTHING_RIDGE = 1.0e-5
# Extra second-difference smoothing applied ONLY to the top of the energy
# axis. Above ~2 MeV the training grid is sparse (9000/10000/12000 keV) while
# 32 Bernstein terms span the range, so the fit has more effective resolution
# than the data supports and bows between the top two nodes: undamped, the
# anchored surface runs +0.40% (z = +8.6 on the 8-cell mean) at 11 MeV.
# Damping just those coefficients removes ~45% of that bias while
# improving every other held-out group and the leave-one-energy-out score,
# and leaves the sub-130 keV region untouched by construction. This is a
# fitting regularizer only — it never enters the payload. Zero reproduces the
# unpenalized fit exactly.
RATE_HIGH_ENERGY_RIDGE_KAPPA = 100.0
RATE_HIGH_ENERGY_RIDGE_CUT = 0.75
SHAPE_ROUGHNESS = 0.0
SHAPE_ROUGHNESS_REFERENCE_COUNTS = 10_000.0
PAIR_YIELD_BANDWIDTH_SCALE = 3.0
# Residual iterations turning each yield field's node values into
# coefficients its kernel average reproduces (interpolate_yield_nodes).
YIELD_NODE_ITERATIONS = 16
# Pair rates need broader interpolation than the ordinary interaction classes,
# but broadening the independently aligned shapes washed out escape shoulders.
# The held-out bandwidth ablation therefore keeps the native shape bandwidth.
PAIR_SHAPE_BANDWIDTH_SCALE = 1.0
# The charged-source helper implements the same normalized cubic B-spline
# algebra with its module-level basis size. Override it only in this detector
# fitter process; the compiled source/charged payload remains 32 terms.
basis_helpers.BASIS_SIZE = BASIS_SIZE
PAIR_THRESHOLD_KEV = 2.0 * basis_helpers.ELECTRON_MASS_KEV
PAIR_THRESHOLD_SCALE_KEV = 50.0
PAIR_RATE_ONSET_SCALE_KEV = basis_helpers.ELECTRON_MASS_KEV
COMPONENTS = (
    "zero_deposit",
    "full_energy",
    "pair_fraction",
    "pair_full_energy",
    "single_escape",
    "double_escape",
    "annihilation_511",
    "k_alpha_escape",
    "k_beta_escape",
    "k_alpha_absorbed",
    "k_beta_absorbed",
    "photoelectric_escape",
    "single_compton",
    "multiple_compton",
    "pair_no_compton",
    "pair_single_compton",
    "pair_multiple_compton",
)
K_ESCAPE_COMPONENTS = (
    "k_alpha_escape", "k_beta_escape",
    "k_alpha_absorbed", "k_beta_absorbed")
# The detector the banks were made with: its label, its geometry (the
# CEELO_GEOMETRY text) and what the anchor helper derives from it -- the
# passive-layer landmark of the features' third coordinate, and the crystal K
# edge gating the K-escape atoms (the escape-line energies themselves come
# from the campaign atom rows, so they always match what transport sampled).
# Set once per run by set_detector(), from the first bank's manifest or a
# loaded payload, before any cell feature is computed.
DETECTOR_NAME: str | None = None
GEOMETRY_TEXT: str | None = None
PASSIVE_LANDMARK_KEV = 55.0
K_EDGE_KEV = 11.103
K_LINE_ELEMENT = 32
ANCHOR_EXECUTABLE: Path | None = None
# A loaded payload's finite_distance record (text), written back unchanged; a
# fresh fit has none ("finite_distance 0") until the bake measures it.
FINITE_DISTANCE_TEXT = "finite_distance 0\n"
_GEOMETRY_FILE: Path | None = None
OUTPUT_COMPONENTS = tuple(
    name for name in COMPONENTS if name != "pair_fraction")
CONTINUOUS = {
    "photoelectric_escape", "single_compton", "multiple_compton",
    "pair_no_compton", "pair_single_compton",
    "pair_multiple_compton"}
PAIR_OUTCOMES = {
    "pair_full_energy",
    "single_escape", "double_escape", "annihilation_511",
    "pair_no_compton", "pair_single_compton",
    "pair_multiple_compton"}
NONPAIR_OUTCOMES = {
    "full_energy", "photoelectric_escape",
    "k_alpha_escape", "k_beta_escape",
    "k_alpha_absorbed", "k_beta_absorbed",
    "single_compton", "multiple_compton"}
PAIR_COMPONENTS = PAIR_OUTCOMES | {"pair_fraction"}
_BASIS_CACHE: dict[tuple[str, str], np.ndarray] = {}
_SHAPE_CACHE: dict[tuple[str, str], np.ndarray] = {}
COMPONENT_BASIS_SIZES: dict[str, int] = {
    name: (BASIS_SIZE if name in CONTINUOUS else 1)
    for name in COMPONENTS
}
# Block partition of each continuum's deposited-energy support.  Block
# boundaries are the component's moving physical landmarks; the coefficient
# count per block is fixed so a coefficient keeps its meaning as incident
# energy changes and the incident-energy interpolation stays well posed.
#
# A single C2 spline over one warped coordinate cannot represent the step at
# the free-electron Compton edge, and with nonnegative coefficients the wide
# below-edge bases whose support crosses the edge get suppressed rather than
# cancelled - which produced a below-edge deficit and an edge spike.  Giving
# each side its own basis removes that coupling.  Allocations were chosen by
# single-cell NNLS chi-square over the training banks.
COMPTON_BLOCK_COMPONENTS = frozenset({"single_compton", "multiple_compton"})
PAIR_BLOCK_COMPONENTS = frozenset({
    "pair_no_compton", "pair_single_compton", "pair_multiple_compton"})
# Fraction of a Compton component's coefficients placed below the edge. The
# above-edge region is narrow at high incident energy and wide at low incident
# energy, so a two-thirds split is the compromise that scored best over the
# training banks at both ends.
COMPTON_BELOW_EDGE_FRACTION = 2.0 / 3.0
MINIMUM_BLOCK_SIZE = 4
def component_block_sizes(name: str, basis_size: int) -> tuple[int, ...]:
    """Per-block coefficient counts, or () when the component is unpartitioned."""
    if name in COMPTON_BLOCK_COMPONENTS:
        below = int(round(basis_size * COMPTON_BELOW_EDGE_FRACTION))
        below = max(MINIMUM_BLOCK_SIZE,
                    min(basis_size - MINIMUM_BLOCK_SIZE, below))
        if basis_size < 2 * MINIMUM_BLOCK_SIZE:
            return ()
        return (below, basis_size - below)
    if name in PAIR_BLOCK_COMPONENTS:
        if basis_size < 4 * MINIMUM_BLOCK_SIZE:
            return ()
        per = basis_size // 4
        return (per, per, per, basis_size - 3 * per)
    return ()


@dataclass
class Cell:
    identity: str
    path: Path
    energy_keV: float
    cos_theta: float
    events: float
    edges: np.ndarray
    values: dict[str, np.ndarray]
    variances: dict[str, np.ndarray]
    yields: dict[str, float]
    yield_variances: dict[str, float]
    # Escape-line offsets Ein - E_atom taken from the campaign's K-escape
    # atom rows (None when the cell is below the crystal K edge).
    k_alpha_offset_keV: float | None = None
    k_beta_offset_keV: float | None = None

    @property
    def feature(self) -> np.ndarray:
        energy = (
            math.log(self.energy_keV / 10.0)
            / math.log(12000.0 / 10.0))
        passive_ratio = self.energy_keV / PASSIVE_LANDMARK_KEV
        return np.asarray([
            energy,
            math.acos(0.99999999 * self.cos_theta) / math.pi,
            1.0 / (1.0 + passive_ratio ** 3),
        ])

    @property
    def pair_feature(self) -> np.ndarray:
        feature = self.feature
        if self.energy_keV > PAIR_THRESHOLD_KEV:
            feature[0] = (
                math.log1p(
                    (self.energy_keV - PAIR_THRESHOLD_KEV)
                    / PAIR_THRESHOLD_SCALE_KEV)
                / math.log1p(
                    (12000.0 - PAIR_THRESHOLD_KEV)
                    / PAIR_THRESHOLD_SCALE_KEV))
        return feature


@dataclass
class Node:
    feature: np.ndarray
    log_yield: float
    shape: np.ndarray


@dataclass
class ComponentModel:
    yield_bandwidth: np.ndarray
    shape_bandwidth: np.ndarray
    basis_size: int
    nodes: list[Node]


def anchor_helper_path(executable: Path | None = None) -> Path:
    """Locate `ceelo_rate_surface_anchor`.

    The path is deliberately NOT serialized into the payload (it is
    machine-specific), so a loaded anchor has to find the helper itself: an
    explicit path, ANCHOR_EXECUTABLE or CEELO_RATE_SURFACE_ANCHOR first, then
    the usual build directories relative to the working directory and to this
    file's repository root, so evaluating a payload does not require being in
    the repo root.
    """
    candidates: list[Path] = []
    for explicit in (executable, ANCHOR_EXECUTABLE):
        if explicit is not None:
            candidates.append(Path(explicit))
    override = os.environ.get("CEELO_RATE_SURFACE_ANCHOR")
    if override:
        candidates.append(Path(override))
    relative = Path("tools/full_spectrum/ceelo_rate_surface_anchor")
    roots = [Path.cwd(), Path(__file__).resolve().parents[2]]
    for root in roots:
        for build in ("build-rel", "build", "build-rel-ceres"):
            candidates.append(root / build / relative)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "cannot find ceelo_rate_surface_anchor; build it or set "
        "CEELO_RATE_SURFACE_ANCHOR. Tried: "
        + ", ".join(str(candidate) for candidate in candidates))


def set_detector(name: str, geometry_text: str,
                 landmark_keV: float | None = None) -> None:
    """Adopt the detector every bank and payload of this run describes.

    The helper derives the passive landmark (unless a loaded payload already
    fixes it), the crystal K edge and its K-line element from the geometry.
    """
    global DETECTOR_NAME, GEOMETRY_TEXT, PASSIVE_LANDMARK_KEV, K_EDGE_KEV
    global K_LINE_ELEMENT, _GEOMETRY_FILE
    handle = tempfile.NamedTemporaryFile(
        "w", suffix=".geometry", delete=False)
    with handle:
        handle.write(geometry_text)
    path = Path(handle.name)
    atexit.register(path.unlink, missing_ok=True)
    completed = subprocess.run(
        [str(anchor_helper_path()), "--geometry", str(path), "--info"],
        capture_output=True, text=True, check=True)
    info = json.loads(completed.stdout)
    DETECTOR_NAME = name
    GEOMETRY_TEXT = geometry_text
    PASSIVE_LANDMARK_KEV = (
        info["passive_landmark_keV"] if landmark_keV is None
        else landmark_keV)
    K_EDGE_KEV = info["k_edge_keV"]
    K_LINE_ELEMENT = int(info["k_line_element"])
    _GEOMETRY_FILE = path
    _ANCHOR_CACHE.clear()


def adopt_detector(name: str, geometry_text: str,
                   landmark_keV: float | None = None,
                   source: str = "") -> None:
    """set_detector() on first use; afterwards, insist on the same detector."""
    if GEOMETRY_TEXT is None:
        set_detector(name, geometry_text, landmark_keV)
    elif geometry_text != GEOMETRY_TEXT or name != DETECTOR_NAME:
        raise ValueError(
            f"{source} describes detector {name!r}, but this run's banks "
            f"and payloads are of {DETECTOR_NAME!r} with another geometry")


def adopt_campaign_detector(root: Path) -> None:
    manifest_path = root / "campaign_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    geometry_text = manifest.get("geometry_text")
    if geometry_text is None:
        raise ValueError(
            f"{manifest_path} records no detector geometry; re-run that "
            "bank's step (finished cells are reused) to record it")
    adopt_detector(manifest["detector"], geometry_text,
                   source=str(manifest_path))


@dataclass
class RateSurfaceAnchor:
    """Deterministic surface-entry transmission anchor on the log odds.

    The anchor has no fitted parameters -- only the detector geometry, which
    is the run's (set_detector) -- and is evaluated by the C++ helper
    `ceelo_rate_surface_anchor`, which calls the same
    `rate_surface_anchor_log_odds` the C++ evaluator uses. Sharing one
    implementation is deliberate: the fitter subtracts exactly what the
    evaluator adds back, so a re-implementation here would put any mismatch
    straight onto the detected rate.
    """
    executable: Path | None = None
    blend_keV: float = 150.0
    blend_power: float = 2.0
    azimuth_quadrature_points: int = 48

    def options(self) -> list[str]:
        if _GEOMETRY_FILE is None:
            raise ValueError("the rate-surface anchor needs set_detector() first")
        return [
            "--geometry", str(_GEOMETRY_FILE),
            "--blend-kev", repr(self.blend_keV),
            "--blend-power", repr(self.blend_power),
            "--azimuth-points", str(self.azimuth_quadrature_points),
        ]

    def payload_record(self) -> str:
        return (
            f"rate_surface_anchor {self.blend_keV!r} {self.blend_power!r} "
            f"{self.azimuth_quadrature_points}")

    def resolved_executable(self) -> Path:
        return anchor_helper_path(self.executable)

    def cache_key(self) -> tuple:
        """Identity for caching: the geometry, never the object address."""
        return (
            GEOMETRY_TEXT, self.blend_keV, self.blend_power,
            self.azimuth_quadrature_points)

    def log_odds(self, pairs: list[tuple[float, float]]) -> np.ndarray:
        """L0 for a list of (energy_keV, cos_theta) pairs."""
        if not pairs:
            return np.zeros(0)
        # Keyed on the geometry, NOT id(self): CPython reuses addresses, so
        # an id-keyed cache silently hands one detector's anchor to another.
        key = (self.cache_key(), tuple(pairs))
        cached = _ANCHOR_CACHE.get(key)
        if cached is not None:
            return cached
        query = "\n".join(
            f"{float(energy)!r} {float(cosine)!r}"
            for energy, cosine in pairs)
        result = subprocess.run(
            [str(self.resolved_executable()), *self.options()],
            input=query, capture_output=True, text=True, check=True)
        if result.stderr:
            print(
                f"rate-surface anchor helper: {result.stderr.strip()}",
                file=sys.stderr)
        values = np.asarray(
            [float(token) for token in result.stdout.split()])
        if values.size != len(pairs):
            raise ValueError(
                "rate-surface anchor helper returned "
                f"{values.size} values for {len(pairs)} queries")
        _ANCHOR_CACHE[key] = values
        return values


_ANCHOR_CACHE: dict[tuple, np.ndarray] = {}


@dataclass
class RateSurface:
    energy_warp_alpha: float
    coefficients: np.ndarray
    anchor: RateSurfaceAnchor | None = None


def rate_energy_coordinate(
    energy_keV: float,
    alpha: float,
) -> float:
    energy_u = math.log(energy_keV / 10.0) / math.log(1200.0)
    passive_at_ten = 1.0 / (1.0 + (10.0 / PASSIVE_LANDMARK_KEV) ** 3)
    passive = 1.0 / (1.0 + (energy_keV / PASSIVE_LANDMARK_KEV) ** 3)
    return (
        energy_u + alpha * (passive_at_ten - passive)
    ) / (1.0 + alpha * passive_at_ten)


def bernstein_values(coordinate: float, size: int) -> np.ndarray:
    degree = size - 1
    return np.asarray([
        math.comb(degree, index)
        * coordinate ** index
        * (1.0 - coordinate) ** (degree - index)
        for index in range(size)
    ])


def rate_angle_coordinate(cos_theta: float) -> float:
    """Reshaped polar angle in [0, 1].

    Mirrors the angle warp in DetectorPulseResponseModel::evaluate; the two
    must stay in lockstep.
    """
    polar = math.acos(0.99999999 * cos_theta) / math.pi
    centred = 2.0 * polar - 1.0
    power = 0.5 * (RATE_ANGLE_WARP_EXPONENT - 1.0)
    epsilon2 = RATE_ANGLE_WARP_EPSILON * RATE_ANGLE_WARP_EPSILON
    shaped = centred * (centred * centred + epsilon2) ** power
    shaped /= (1.0 + epsilon2) ** power
    return 0.5 * (1.0 + shaped)


def rate_surface_row(
    energy_keV: float,
    cos_theta: float,
    alpha: float,
) -> np.ndarray:
    energy = bernstein_values(
        rate_energy_coordinate(energy_keV, alpha),
        RATE_ENERGY_BASIS_SIZE)
    angle = bernstein_values(
        rate_angle_coordinate(cos_theta), RATE_ANGLE_BASIS_SIZE)
    return np.kron(angle, energy)


def fit_rate_surface(
    cells: list[Cell],
    anchor: RateSurfaceAnchor | None = None,
) -> RateSurface:
    populated = []
    for cell in cells:
        detected = max(1.0 - cell.yields["zero_deposit"], 1e-300)
        if detected * cell.events >= 25.0:
            populated.append((cell, detected))
    matrix = np.vstack([
        rate_surface_row(
            cell.energy_keV, cell.cos_theta,
            RATE_ENERGY_WARP_ALPHA)
        for cell, _ in populated
    ])
    target = np.asarray([
        math.log(detected / max(1.0 - detected, 1e-300))
        for _, detected in populated
    ])
    if anchor is not None:
        # Fit only the residual above the deterministic anchor; the evaluator
        # adds the identical quantity back.
        target = target - anchor.log_odds(
            [(cell.energy_keV, cell.cos_theta) for cell, _ in populated])
    surface_size = RATE_ENERGY_BASIS_SIZE * RATE_ANGLE_BASIS_SIZE
    penalty_rows = []
    for angle in range(RATE_ANGLE_BASIS_SIZE):
        for energy in range(RATE_ENERGY_BASIS_SIZE - 2):
            row = np.zeros(surface_size)
            offset = angle * RATE_ENERGY_BASIS_SIZE + energy
            fraction = (energy + 1) / (RATE_ENERGY_BASIS_SIZE - 1)
            scale = math.sqrt(
                1.0 + RATE_HIGH_ENERGY_RIDGE_KAPPA
                if fraction > RATE_HIGH_ENERGY_RIDGE_CUT else 1.0)
            row[offset:offset + 3] = (scale, -2.0 * scale, scale)
            penalty_rows.append(row)
    for angle in range(RATE_ANGLE_BASIS_SIZE - 2):
        for energy in range(RATE_ENERGY_BASIS_SIZE):
            row = np.zeros(surface_size)
            row[angle * RATE_ENERGY_BASIS_SIZE + energy] = 1.0
            row[(angle + 1) * RATE_ENERGY_BASIS_SIZE + energy] = -2.0
            row[(angle + 2) * RATE_ENERGY_BASIS_SIZE + energy] = 1.0
            penalty_rows.append(row)
    penalty = np.vstack(penalty_rows)
    normal = matrix.T @ matrix
    normal += RATE_SMOOTHING_RIDGE * (penalty.T @ penalty)
    normal += 1.0e-12 * np.eye(surface_size)
    coefficients = np.linalg.solve(normal, matrix.T @ target)
    return RateSurface(
        energy_warp_alpha=RATE_ENERGY_WARP_ALPHA,
        coefficients=coefficients,
        anchor=anchor)


def evaluate_rate_surface(
    cell: Cell,
    surface: RateSurface,
) -> float:
    log_odds = float(
        rate_surface_row(
            cell.energy_keV, cell.cos_theta,
            surface.energy_warp_alpha)
        @ surface.coefficients)
    if surface.anchor is not None:
        log_odds += float(surface.anchor.log_odds(
            [(cell.energy_keV, cell.cos_theta)])[0])
    return 1.0 / (1.0 + math.exp(-log_odds))


def load_cell(path: Path) -> Cell:
    with bank_io.open_cell(path) as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty cell: {path}")
    first = rows[0]
    if first["cell_kind"] != "detector":
        raise ValueError(f"not a detector response cell: {path}")
    boundaries = sorted({
        float(value)
        for row in rows
        for value in (row["bin_lo_keV"], row["bin_hi_keV"])
        if float(row["bin_hi_keV"]) > float(row["bin_lo_keV"])
    })
    edges = np.asarray(boundaries)
    values = {
        name: np.zeros(edges.size - 1) for name in CONTINUOUS}
    variances = {
        name: np.zeros(edges.size - 1) for name in CONTINUOUS}
    yields = {name: 0.0 for name in COMPONENTS}
    yield_variances = {name: 0.0 for name in COMPONENTS}
    events = float(first["events"])
    energy_keV = float(first["energy_keV"])
    k_alpha_offset = None
    k_beta_offset = None
    for row in rows:
        name = row["component"]
        if name not in yields:
            continue
        value = float(row["sum_w"]) / events
        variance = float(row["sum_w2"]) / (events * events)
        yields[name] += value
        yield_variances[name] += variance
        if name == "k_alpha_escape":
            k_alpha_offset = energy_keV - float(row["bin_lo_keV"])
        elif name == "k_beta_escape":
            k_beta_offset = energy_keV - float(row["bin_lo_keV"])
        if name not in CONTINUOUS:
            continue
        lo = float(row["bin_lo_keV"])
        index = int(np.searchsorted(edges, lo, side="right") - 1)
        if 0 <= index < edges.size - 1:
            values[name][index] += value
            variances[name][index] += variance
    yields["pair_fraction"] = sum(
        yields[name] for name in PAIR_OUTCOMES)
    yield_variances["pair_fraction"] = sum(
        yield_variances[name] for name in PAIR_OUTCOMES)
    return Cell(
        identity=bank_io.cell_id(path),
        path=path,
        energy_keV=energy_keV,
        cos_theta=float(first["cos_theta"]),
        events=events,
        edges=edges,
        values=values,
        variances=variances,
        yields=yields,
        yield_variances=yield_variances,
        k_alpha_offset_keV=k_alpha_offset,
        k_beta_offset_keV=k_beta_offset,
    )


def load_campaign(root: Path) -> list[Cell]:
    adopt_campaign_detector(root)
    return [
        load_cell(path)
        for path in bank_io.cell_files(root / "cells")
    ]


def prefer_highest_statistics(cells: list[Cell]) -> list[Cell]:
    """Deduplicate whole transport cells, retaining the largest event bank."""
    selected: dict[tuple[float, float], Cell] = {}
    for cell in cells:
        key = (cell.energy_keV, cell.cos_theta)
        current = selected.get(key)
        if current is None or cell.events > current.events:
            selected[key] = cell
    return sorted(
        selected.values(),
        key=lambda cell: (cell.energy_keV, cell.cos_theta))


def component_support(cell: Cell, name: str) -> float:
    # Interaction-order labels span the complete composite detector assembly.
    # A photon can Compton scatter and subsequently transfer more energy after
    # crossing active, dead, or passive regions, so the free-electron Compton
    # edge is a moving feature, not a hard support limit.
    return cell.energy_keV


def block_record(name: str, basis_size: int) -> str:
    """Serialized block partition: count then per-block coefficient counts."""
    blocks = component_block_sizes(name, basis_size)
    if not blocks:
        return " 0"
    return f" {len(blocks)}" + "".join(f" {size}" for size in blocks)


def block_offsets(name: str, basis_size: int) -> list[tuple[int, int]]:
    """(offset, size) of each block, or one whole-basis block if unpartitioned."""
    blocks = component_block_sizes(name, basis_size)
    if not blocks:
        return [(0, basis_size)]
    offsets = []
    offset = 0
    for size in blocks:
        offsets.append((offset, size))
        offset += size
    return offsets


def block_bounds(cell: Cell, name: str) -> list[float]:
    """Moving landmark boundaries of a block-partitioned component.

    Mirrors detector_compton_block_layout / detector_pair_block_layout in
    src/full_spectrum/FullSpectrumModel.h; the two must stay in lockstep.
    """
    energy = cell.energy_keV
    if name in COMPTON_BLOCK_COMPONENTS:
        edge = energy * (
            1.0 - 1.0
            / (1.0 + 2.0 * energy / basis_helpers.ELECTRON_MASS_KEV))
        return [0.0, edge, energy]
    points = sorted([
        0.0,
        energy - PAIR_THRESHOLD_KEV,
        basis_helpers.ELECTRON_MASS_KEV,
        energy - basis_helpers.ELECTRON_MASS_KEV,
        energy,
    ])
    return [min(max(point, 0.0), energy) for point in points]


def block_partitioned_basis(
    cell: Cell,
    name: str,
    blocks: tuple[int, ...],
) -> np.ndarray:
    """Concatenated per-block cubic B-spline integral matrices.

    A block that the training histogram does not resolve - fewer than four bin
    edges strictly inside it, which happens where landmarks nearly coincide
    just above the pair threshold or at the 511 / Ein-1022 crossing - is
    emitted as zero columns. Its basis functions would otherwise be mutually
    indistinguishable over the bins they touch, making the least-squares
    problem exactly collinear rather than merely underdetermined. The
    coefficient count is unchanged, so those coefficients simply stay zero and
    the evaluator reproduces the same (vanishing) contribution.
    """
    bounds = block_bounds(cell, name)
    columns = []
    for index, size in enumerate(blocks):
        low, high = bounds[index], bounds[index + 1]
        interior = int(np.count_nonzero(
            (cell.edges > low + 1.0e-12) & (cell.edges < high - 1.0e-12)))
        if high - low < 1.0e-9 or interior < 4:
            columns.append(np.zeros((cell.edges.size - 1, size)))
            continue
        width = high - low
        lo = np.clip((cell.edges[:-1] - low) / width, 0.0, 1.0)
        hi = np.clip((cell.edges[1:] - low) / width, 0.0, 1.0)
        valid = (cell.edges[1:] > low) & (cell.edges[:-1] < high)
        previous_basis_size = basis_helpers.BASIS_SIZE
        basis_helpers.BASIS_SIZE = size
        try:
            columns.append(
                basis_helpers.bspline_matrix_from_coordinates(lo, hi, valid))
        finally:
            basis_helpers.BASIS_SIZE = previous_basis_size
    return np.hstack(columns)


def component_basis(
    cell: Cell,
    name: str,
    basis_size: int | None = None,
) -> np.ndarray:
    if basis_size is None:
        basis_size = COMPONENT_BASIS_SIZES[name]
    key = (str(cell.path), f"{name}:{basis_size}")
    cached = _BASIS_CACHE.get(key)
    if cached is not None:
        return cached
    support = component_support(cell, name)
    blocks = component_block_sizes(name, basis_size)
    if blocks:
        answer = block_partitioned_basis(cell, name, blocks)
        _BASIS_CACHE[key] = answer
        return answer
    # Anything without a block partition (photoelectric escape) uses a plain
    # linear-energy coordinate.
    valid = (cell.edges[:-1] < support) & (cell.edges[1:] > 0.0)
    lo = np.clip(cell.edges[:-1] / support, 0.0, 1.0)
    hi = np.clip(cell.edges[1:] / support, 0.0, 1.0)
    previous_basis_size = basis_helpers.BASIS_SIZE
    basis_helpers.BASIS_SIZE = basis_size
    try:
        answer = basis_helpers.bspline_matrix_from_coordinates(
            lo, hi, valid)
    finally:
        basis_helpers.BASIS_SIZE = previous_basis_size
    _BASIS_CACHE[key] = answer
    return answer


def shape_coefficients(cell: Cell, name: str) -> np.ndarray:
    key = (str(cell.path), name)
    cached = _SHAPE_CACHE.get(key)
    if cached is not None:
        return cached
    if name not in CONTINUOUS:
        answer = np.zeros(COMPONENT_BASIS_SIZES[name])
        answer[0] = 1.0
        _SHAPE_CACHE[key] = answer
        return answer
    total = cell.yields[name]
    if not total > 0.0:
        size = COMPONENT_BASIS_SIZES[name]
        answer = np.full(size, 1.0 / size)
        _SHAPE_CACHE[key] = answer
        return answer
    basis = component_basis(cell, name)
    target = cell.values[name] / total
    basis_size = COMPONENT_BASIS_SIZES[name]
    # Columns for a block the histogram cannot resolve are identically zero
    # (see block_partitioned_basis). Solve without them and scatter the answer
    # back, so the coefficient count stays fixed and the least-squares problem
    # stays well posed. Determined before any penalty rows are stacked on.
    supported = basis.any(axis=0)
    if SHAPE_ROUGHNESS > 0.0 and basis_size > 2:
        # Curvature is penalized inside each block only. A stencil spanning a
        # block boundary would smooth across a physical step (the Compton
        # edge, the escape edges) - exactly the structure the partition exists
        # to represent.
        rows = []
        for offset, size in block_offsets(name, basis_size):
            for index in range(size - 2):
                row = np.zeros(basis_size)
                row[offset + index:offset + index + 3] = (1.0, -2.0, 1.0)
                rows.append(row)
        penalty = np.vstack(rows) if rows else np.zeros((0, basis_size))
        component_counts = max(total * cell.events, 1.0)
        effective_roughness = SHAPE_ROUGHNESS * max(
            1.0,
            SHAPE_ROUGHNESS_REFERENCE_COUNTS / component_counts)
        basis = np.vstack((
            basis,
            math.sqrt(effective_roughness) * penalty))
        target = np.concatenate((target, np.zeros(penalty.shape[0])))
    if supported.all():
        coefficients, _ = basis_helpers.nonnegative_least_squares(
            basis, target)
    else:
        coefficients = np.zeros(basis_size)
        solved, _ = basis_helpers.nonnegative_least_squares(
            basis[:, supported], target)
        coefficients[supported] = solved
    coefficient_sum = coefficients.sum()
    if not coefficient_sum > 0.0:
        answer = np.full(basis_size, 1.0 / basis_size)
        _SHAPE_CACHE[key] = answer
        return answer
    answer = coefficients / coefficient_sum
    _SHAPE_CACHE[key] = answer
    return answer


def interpolate_yield_nodes(model: ComponentModel, iterations: int) -> None:
    """Replace the yield nodes' training values with coefficients the kernel
    average reproduces at the nodes.

    The yield field is a normalized Gaussian-kernel average of node log yields
    (here and in the C++ evaluator).  With the training values y as its
    coefficients it is a smoother: at node i it returns (S y)_i, biased
    wherever the field curves within a bandwidth or the nodes crowd on one
    side -- the full-energy atom read up to +3% at its own training cells.
    Damped residual iteration v <- v + (y - S v) (Landweber; one step is
    Tukey's twicing) drives S v toward y.  S's eigenvalues lie in [0, 1], so
    after k steps the field at the nodes is the data filtered by
    1 - (1 - lambda)^(k + 1) <= 1: the bias goes mode by mode, and the noise
    is never amplified beyond the training cells' own.
    """
    if iterations <= 0 or len(model.nodes) < 2:
        return
    centers = np.vstack([node.feature for node in model.nodes])
    scaled = centers / model.yield_bandwidth
    distance2 = np.sum(
        (scaled[:, None, :] - scaled[None, :, :]) ** 2, axis=2)
    smoother = np.exp(
        -0.5 * (distance2 - np.min(distance2, axis=1, keepdims=True)))
    smoother /= smoother.sum(axis=1, keepdims=True)
    target = np.asarray([node.log_yield for node in model.nodes])
    coefficients = target.copy()
    for _ in range(iterations):
        coefficients += target - smoother @ coefficients
    for node, value in zip(model.nodes, coefficients):
        node.log_yield = float(value)


def fit_components(
    cells: list[Cell],
    yield_bandwidth: np.ndarray,
    shape_bandwidth: np.ndarray,
    minimum_component_counts: float,
    component_yield_scales: dict[str, float] | None = None,
    yield_node_iterations: int = 0,
) -> dict[str, ComponentModel]:
    models = {}
    for name in COMPONENTS:
        nodes = []
        for cell in cells:
            detected_yield = max(
                1.0 - cell.yields["zero_deposit"], 1e-300)
            if name == "zero_deposit":
                # Store the log odds of a nonzero deposit. Interpolating a
                # near-unity zero-deposit probability and subtracting it from
                # one catastrophically amplifies low-energy rate errors.
                modeled_yield = (
                    detected_yield
                    / max(1.0 - detected_yield, 1e-300))
                expected_counts = detected_yield * cell.events
            elif name == "pair_fraction":
                pair_yield = cell.yields["pair_fraction"]
                nonpair_yield = max(
                    detected_yield - pair_yield, 1e-300)
                pair_odds = pair_yield / nonpair_yield
                onset_odds = max(
                    (cell.energy_keV - PAIR_THRESHOLD_KEV)
                    / PAIR_RATE_ONSET_SCALE_KEV,
                    1e-300)
                # Fit only the positive residual to a mechanistic threshold
                # onset. Evaluation multiplies this residual by onset_odds,
                # enforcing exactly zero pair rate at threshold.
                modeled_yield = pair_odds / onset_odds
                expected_counts = pair_yield * cell.events
            elif name in PAIR_OUTCOMES:
                pair_yield = max(
                    cell.yields["pair_fraction"], 1e-300)
                modeled_yield = cell.yields[name] / pair_yield
                expected_counts = cell.yields[name] * cell.events
            else:
                nonpair_yield = max(
                    detected_yield
                    - cell.yields["pair_fraction"], 1e-300)
                modeled_yield = cell.yields[name] / nonpair_yield
                expected_counts = cell.yields[name] * cell.events
            if (
                name in PAIR_COMPONENTS
                and cell.energy_keV <= PAIR_THRESHOLD_KEV
            ):
                continue
            if expected_counts < minimum_component_counts:
                continue
            nodes.append(Node(
                feature=(
                    cell.pair_feature
                    if name in PAIR_COMPONENTS else cell.feature),
                log_yield=math.log(max(modeled_yield, 1e-300)),
                shape=shape_coefficients(cell, name),
            ))
        if not nodes:
            if name in K_ESCAPE_COMPONENTS:
                # A bank generated before the K-escape labels (or below the
                # crystal K edge throughout) simply has no atoms; write an
                # empty component so the evaluator keeps them disabled.
                models[name] = ComponentModel(
                    yield_bandwidth=yield_bandwidth,
                    shape_bandwidth=shape_bandwidth,
                    basis_size=COMPONENT_BASIS_SIZES[name],
                    nodes=[],
                )
                continue
            raise ValueError(f"no populated nodes for {name}")
        models[name] = ComponentModel(
            yield_bandwidth=(
                yield_bandwidth
                * (
                    PAIR_YIELD_BANDWIDTH_SCALE
                    if name in PAIR_COMPONENTS else 1.0)
                * (component_yield_scales or {}).get(name, 1.0)),
            shape_bandwidth=(
                shape_bandwidth
                * (
                    PAIR_SHAPE_BANDWIDTH_SCALE
                    if name in PAIR_COMPONENTS else 1.0)),
            basis_size=COMPONENT_BASIS_SIZES[name],
            nodes=nodes,
        )
        # At the component's final bandwidth: the coefficients reproduce the
        # nodes only for the kernel they were solved with.
        interpolate_yield_nodes(models[name], yield_node_iterations)
    return models


def interpolate(
    model: ComponentModel,
    feature: np.ndarray,
    shape: bool,
) -> np.ndarray | float:
    bandwidth = (
        model.shape_bandwidth if shape else model.yield_bandwidth)
    centers = np.vstack([node.feature for node in model.nodes])
    distance2 = np.sum(
        ((centers - feature) / bandwidth) ** 2, axis=1)
    weights = np.exp(-0.5 * (distance2 - np.min(distance2)))
    weights /= weights.sum()
    if shape:
        coefficients = weights @ np.vstack([
            node.shape for node in model.nodes])
        return coefficients / max(coefficients.sum(), 1e-300)
    return math.exp(float(weights @ np.asarray([
        node.log_yield for node in model.nodes])))


def k_escape_applies(
    cell: Cell,
    name: str,
    models: dict[str, ComponentModel],
    k_escape_lines: tuple[float, float, float, float] | None,
) -> bool:
    """Whether a K-escape atom exists for this cell.

    Must agree with `detector_spectrum`, which places the atom only when the
    cell carries the campaign row, and with the C++ evaluator, which gates on
    the payload's edge and support limit. A mismatch would allocate yield to
    an atom that never lands in the modeled spectrum.
    """
    if not models[name].nodes:
        return False
    if cell.k_alpha_offset_keV is None or cell.k_beta_offset_keV is None:
        return False
    edge = (
        k_escape_lines[2]
        if k_escape_lines and k_escape_lines[2] > 0.0 else K_EDGE_KEV)
    if cell.energy_keV <= edge:
        return False
    # Each component has its own support; they differ by up to 2x.
    if not k_escape_lines:
        return True
    support = k_escape_lines[3 + K_ESCAPE_COMPONENTS.index(name)]
    return support <= 0.0 or cell.energy_keV <= support


def predict(
    cell: Cell,
    models: dict[str, ComponentModel],
    rate_surface: RateSurface | None = None,
    k_escape_lines: tuple[float, float, float, float] | None = None,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    raw_yields = {
        name: float(interpolate(
            models[name],
            cell.pair_feature
            if name in PAIR_COMPONENTS else cell.feature,
            False))
        for name in COMPONENTS
        if name != "zero_deposit"
        and not (
            name in PAIR_COMPONENTS
            and cell.energy_keV <= PAIR_THRESHOLD_KEV)
        and not (
            name in K_ESCAPE_COMPONENTS
            and not k_escape_applies(cell, name, models, k_escape_lines))
    }
    if rate_surface is not None:
        detected_yield = evaluate_rate_surface(cell, rate_surface)
    else:
        detected_odds = float(interpolate(
            models["zero_deposit"], cell.feature, False))
        detected_yield = detected_odds / (1.0 + detected_odds)
    if cell.energy_keV > PAIR_THRESHOLD_KEV:
        pair_odds = (
            raw_yields["pair_fraction"]
            * (cell.energy_keV - PAIR_THRESHOLD_KEV)
            / PAIR_RATE_ONSET_SCALE_KEV)
        pair_probability = pair_odds / (1.0 + pair_odds)
    else:
        pair_probability = 0.0
    nonpair_probability = 1.0 - pair_probability
    nonpair_normalization = max(
        sum(raw_yields.get(name, 0.0) for name in NONPAIR_OUTCOMES),
        1e-300)
    pair_normalization = max(
        sum(raw_yields.get(name, 0.0) for name in PAIR_OUTCOMES),
        1e-300)
    yields = {
        "zero_deposit": 1.0 - detected_yield,
        "pair_fraction": detected_yield * pair_probability,
    }
    yields.update({
        name: (
            detected_yield * nonpair_probability
            * raw_yields.get(name, 0.0) / nonpair_normalization)
        for name in NONPAIR_OUTCOMES
    })
    yields.update({
        name: (
            detected_yield * pair_probability
            * raw_yields.get(name, 0.0) / pair_normalization)
        for name in PAIR_OUTCOMES
    })
    spectra = {}
    for name in CONTINUOUS:
        coefficients = np.asarray(
            interpolate(
                models[name],
                cell.pair_feature
                if name in PAIR_COMPONENTS else cell.feature,
                True))
        spectra[name] = (
            component_basis(
                cell, name, models[name].basis_size)
            @ coefficients * yields[name])
    return yields, spectra


def detector_spectrum(
    cell: Cell,
    yields: dict[str, float],
    spectra: dict[str, np.ndarray],
) -> np.ndarray:
    result = sum(
        (spectra[name] for name in CONTINUOUS),
        np.zeros(cell.edges.size - 1))
    atom_energies = {
        "full_energy": cell.energy_keV,
        "pair_full_energy": cell.energy_keV,
        "single_escape": cell.energy_keV
            - basis_helpers.ELECTRON_MASS_KEV,
        "double_escape": cell.energy_keV
            - 2.0 * basis_helpers.ELECTRON_MASS_KEV,
        "annihilation_511": basis_helpers.ELECTRON_MASS_KEV,
    }
    if cell.k_alpha_offset_keV is not None:
        atom_energies["k_alpha_escape"] = (
            cell.energy_keV - cell.k_alpha_offset_keV)
        # Absorbed fluorescence deposits the line energy itself.
        atom_energies["k_alpha_absorbed"] = cell.k_alpha_offset_keV
    if cell.k_beta_offset_keV is not None:
        atom_energies["k_beta_escape"] = (
            cell.energy_keV - cell.k_beta_offset_keV)
        atom_energies["k_beta_absorbed"] = cell.k_beta_offset_keV
    for name, energy in atom_energies.items():
        if energy < 0.0:
            continue
        index = int(np.searchsorted(
            cell.edges, energy, side="right") - 1)
        if result.size:
            index = min(max(index, 0), result.size - 1)
            result[index] += yields[name]
    return result


def truth_spectrum(cell: Cell) -> np.ndarray:
    return detector_spectrum(cell, cell.yields, cell.values)


def finite_summary(values: list[float]) -> dict:
    array = np.asarray([value for value in values if math.isfinite(value)])
    if not array.size:
        return {"count": 0}
    return {
        "count": int(array.size),
        "median": float(np.median(array)),
        "p90": float(np.percentile(array, 90)),
        "p95": float(np.percentile(array, 95)),
        "maximum": float(np.max(array)),
    }


def metrics(
    cells: list[Cell],
    models: dict[str, ComponentModel],
    rate_surface: RateSurface | None = None,
    k_escape_lines: tuple[float, float, float, float] | None = None,
) -> tuple[dict, list[dict]]:
    rows = []
    for cell in cells:
        yields, spectra = predict(
            cell, models, rate_surface, k_escape_lines)
        truth = truth_spectrum(cell)
        estimate = detector_spectrum(cell, yields, spectra)
        truth_detected = max(1.0 - cell.yields["zero_deposit"], 1e-15)
        estimate_detected = max(1.0 - yields["zero_deposit"], 1e-15)
        truth_shape = truth / max(truth.sum(), 1e-300)
        estimate_shape = estimate / max(estimate.sum(), 1e-300)
        continuum_truth = sum(
            (cell.values[name] for name in CONTINUOUS),
            np.zeros(truth.size))
        continuum_estimate = sum(
            (spectra[name] for name in CONTINUOUS),
            np.zeros(truth.size))
        continuum_counts = continuum_truth.sum() * cell.events
        if continuum_truth.sum() > 0.0 and continuum_estimate.sum() > 0.0:
            continuum_tv = 0.5 * float(np.sum(np.abs(
                continuum_truth / continuum_truth.sum()
                - continuum_estimate / continuum_estimate.sum())))
        else:
            continuum_tv = 0.0
        detected_counts = truth_detected * cell.events
        row = {
            "identity": cell.identity,
            "energy_keV": cell.energy_keV,
            "cos_theta": cell.cos_theta,
            "detected_yield_relative_error": (
                estimate_detected / truth_detected - 1.0
                if detected_counts >= 400.0 else math.nan),
            "pulse_shape_total_variation": (
                0.5 * float(np.sum(
                    np.abs(estimate_shape - truth_shape)))
                if detected_counts >= 400.0 else math.nan),
            "continuum_shape_total_variation": (
                continuum_tv
                if continuum_counts >= 400.0 else math.nan),
            "detected_mc_counts": detected_counts,
            "continuum_mc_counts": continuum_counts,
        }
        for name in COMPONENTS:
            row[f"{name}_yield_error"] = (
                yields[name] / cell.yields[name] - 1.0
                if cell.yields[name] * cell.events >= 400.0
                else math.nan)
        for name in CONTINUOUS:
            truth_component = cell.values[name]
            estimate_component = spectra[name]
            truth_area = float(truth_component.sum())
            estimate_area = float(estimate_component.sum())
            row[f"{name}_shape_total_variation"] = (
                0.5 * float(np.sum(np.abs(
                    truth_component / truth_area
                    - estimate_component / estimate_area)))
                if truth_area * cell.events >= 400.0
                and estimate_area > 0.0
                else math.nan)
        rows.append(row)
    report = {}
    for metric in (
        "detected_yield_relative_error",
        "pulse_shape_total_variation",
        "continuum_shape_total_variation",
    ):
        report[metric] = finite_summary([
            abs(float(row[metric])) for row in rows])
    report["component_yield_relative_error"] = {
        name: finite_summary([
            abs(float(row[f"{name}_yield_error"]))
            for row in rows])
        for name in COMPONENTS
    }
    report["component_shape_total_variation"] = {
        name: finite_summary([
            abs(float(row[f"{name}_shape_total_variation"]))
            for row in rows])
        for name in CONTINUOUS
    }
    return report, rows


def k_escape_line_energies(
    training_cells: list[Cell],
    minimum_counts: float = 25.0,
) -> tuple[float, float, float, float]:
    """(E_K_alpha, E_K_beta, K_edge) + a support limit PER COMPONENT.

    Each limit is the highest training energy that produced enough counts to
    place a node (the same `minimum_counts` rule `fit_components` uses), so
    the gate matches each component's own node cloud. They differ materially
    — 320 / 175 / 280 / 150 keV for the fitted HPGe payload — so one shared
    limit would let the shorter-lived components extrapolate flat across most
    of a decade, which is exactly the defect the gate exists to prevent. The fitted yield fields are normalized kernel averages,
    so above their node cloud they level off at the nearest node instead of
    decaying; without the limit the model emits a spurious escape satellite
    beside the full-energy peak all the way to 12 MeV. Zeros disable the
    components entirely (bank without K-escape rows).
    """
    lines = None
    support = {name: 0.0 for name in K_ESCAPE_COMPONENTS}
    for cell in training_cells:
        if (cell.k_alpha_offset_keV is None
                or cell.k_beta_offset_keV is None):
            continue
        if lines is None:
            lines = (cell.k_alpha_offset_keV, cell.k_beta_offset_keV)
        for name in K_ESCAPE_COMPONENTS:
            if cell.yields[name] * cell.events >= minimum_counts:
                support[name] = max(support[name], cell.energy_keV)
    if lines is None:
        return (0.0, 0.0, 0.0) + (0.0,) * len(K_ESCAPE_COMPONENTS)
    return (lines[0], lines[1], K_EDGE_KEV) + tuple(
        support[name] for name in K_ESCAPE_COMPONENTS)


def write_payload(
    path: Path,
    models: dict[str, ComponentModel],
    training_cells: list[Cell],
    validity_keV: tuple[float, float] | None = None,
    rate_surface: RateSurface | None = None,
    precision: int = 12,
    k_escape_lines: tuple[float, float, float] | None = None,
) -> None:
    if validity_keV is None:
        validity_keV = (
            min(cell.energy_keV for cell in training_cells),
            max(cell.energy_keV for cell in training_cells),
        )
    if k_escape_lines is None:
        k_escape_lines = k_escape_line_energies(training_cells)
    if len(k_escape_lines) != 3 + len(K_ESCAPE_COMPONENTS):
        raise ValueError(
            "k_escape_lines needs three energies plus one support limit "
            "per K component")
    number = lambda value: f"{value:.{precision}g}"
    if GEOMETRY_TEXT is None:
        raise ValueError("write_payload needs the run's detector (set_detector)")
    geometry_lines = GEOMETRY_TEXT.splitlines()
    with path.open("w") as output:
        output.write("CEELO_DETECTOR_RESPONSE 1\n")
        # The family picks the evaluator's pilot preset and resolution window:
        # germanium crystals take the HPGe one.
        family = "hpge_coax" if K_LINE_ELEMENT == 32 else "nai3x3"
        output.write(f"family {family}\n")
        output.write(f"detector {DETECTOR_NAME}\n")
        output.write(f"geometry {len(geometry_lines)}\n")
        output.write("".join(line + "\n" for line in geometry_lines))
        output.write(
            f"validity {number(validity_keV[0])} "
            f"{number(validity_keV[1])} -1 1\n")
        output.write(
            "feature_map log_energy_angle_passive_pair_threshold\n")
        # At full precision: the features were computed with this value.
        output.write(f"passive_landmark_keV {PASSIVE_LANDMARK_KEV!r}\n")
        output.write(f"basis_size {max(COMPONENT_BASIS_SIZES.values())}\n")
        output.write(
            f"rate_surface_enabled {int(rate_surface is not None)}\n")
        if rate_surface is not None:
            output.write(
                f"rate_surface {RATE_ENERGY_BASIS_SIZE} "
                f"{RATE_ANGLE_BASIS_SIZE} "
                f"{number(rate_surface.energy_warp_alpha)} "
                f"{number(RATE_ANGLE_WARP_EXPONENT)} "
                f"{number(RATE_ANGLE_WARP_EPSILON)}\n")
            output.write(
                "rate_surface_coefficients"
                + "".join(
                    f" {number(value)}"
                    for value in rate_surface.coefficients)
                + "\n")
            anchor = rate_surface.anchor
            output.write(
                f"rate_surface_anchor_enabled {int(anchor is not None)}\n")
            if anchor is not None:
                output.write(anchor.payload_record() + "\n")
        for name in COMPONENTS:
            if name == "k_alpha_escape":
                output.write(
                    "k_escape_lines "
                    + " ".join(
                        number(value) for value in k_escape_lines[:3])
                    + "\n")
                output.write(
                    "k_escape_support "
                    + " ".join(
                        number(value) for value in k_escape_lines[3:])
                    + "\n")
            model = models[name]
            nodes = (
                [] if name == "zero_deposit"
                and rate_surface is not None else model.nodes)
            output.write(
                name
                + "".join(
                    f" {number(value)}" for value in model.yield_bandwidth)
                + "".join(
                    f" {number(value)}" for value in model.shape_bandwidth)
                + f" {model.basis_size} {len(nodes)}"
                + block_record(name, model.basis_size)
                + "\n")
            for node in nodes:
                output.write(
                    f"{name}_node"
                    + "".join(f" {number(value)}" for value in node.feature)
                    + f" {number(node.log_yield)}"
                    + "".join(f" {number(value)}" for value in node.shape)
                    + "\n")
        output.write(FINITE_DISTANCE_TEXT)
        output.write("END\n")


def read_payload(
    path: Path,
) -> tuple[dict[str, ComponentModel], RateSurface | None,
           tuple[float, float, float]]:
    """Read a compact detector payload without refitting its basis nodes.

    Returns (models, rate surface, K-escape line energies).
    """
    with path.open() as stream:
        lines = iter(stream)
        header = next(lines).split()
        if header != ["CEELO_DETECTOR_RESPONSE", "1"]:
            raise ValueError(
                f"unsupported detector payload (expected version 1): {path}")
        family = next(lines).split()
        if len(family) != 2 or family[0] != "family" or family[1] not in {
                "hpge_coax", "nai3x3"}:
            raise ValueError(f"unsupported detector family: {family}")
        record = next(lines).split()
        if len(record) != 2 or record[0] != "detector":
            raise ValueError(f"{path}: missing the detector label and geometry")
        geometry = next(lines).split()
        if len(geometry) != 2 or geometry[0] != "geometry":
            raise ValueError(f"{path}: invalid detector geometry block")
        geometry_text = "".join(next(lines) for _ in range(int(geometry[1])))
        validity = next(lines).split()
        if len(validity) != 5 or validity[0] != "validity":
            raise ValueError("invalid detector validity record")
        feature_map = next(lines).split()
        if feature_map != [
            "feature_map", "log_energy_angle_passive_pair_threshold",
        ]:
            raise ValueError(
                f"unsupported detector feature map: {feature_map}")
        landmark = next(lines).split()
        if len(landmark) != 2 or landmark[0] != "passive_landmark_keV":
            raise ValueError(f"{path}: missing the passive landmark")
        adopt_detector(record[1], geometry_text, float(landmark[1]),
                       source=str(path))
        basis = next(lines).split()
        if (
            len(basis) != 2
            or basis[0] != "basis_size"
            or int(basis[1]) < 1
        ):
            raise ValueError(f"unsupported detector basis: {basis}")
        rate_surface = None
        enabled = next(lines).split()
        if (
            len(enabled) != 2
            or enabled[0] != "rate_surface_enabled"
            or enabled[1] not in {"0", "1"}
        ):
            raise ValueError(
                f"invalid rate-surface enable flag: {enabled}")
        if enabled[1] == "1":
            surface = next(lines).split()
            if (
                len(surface) != 6
                or surface[0] != "rate_surface"
                or int(surface[1]) != RATE_ENERGY_BASIS_SIZE
                or int(surface[2]) != RATE_ANGLE_BASIS_SIZE
            ):
                raise ValueError(
                    f"unsupported detector rate surface: {surface}")
            coefficients = next(lines).split()
            expected_size = (
                RATE_ENERGY_BASIS_SIZE * RATE_ANGLE_BASIS_SIZE)
            if (
                len(coefficients) != expected_size + 1
                or coefficients[0] != "rate_surface_coefficients"
            ):
                raise ValueError(
                    "invalid detector rate-surface coefficients")
            anchor = None
            flag = next(lines).split()
            if (
                len(flag) != 2
                or flag[0] != "rate_surface_anchor_enabled"
                or flag[1] not in {"0", "1"}
            ):
                raise ValueError(
                    f"invalid rate-surface anchor flag: {flag}")
            if flag[1] == "1":
                anchor = _parse_anchor_record(
                    ["rate_surface_anchor_enabled 1", next(lines)], 0)
            rate_surface = RateSurface(
                energy_warp_alpha=float(surface[3]),
                coefficients=np.asarray([
                    float(value) for value in coefficients[1:]]),
                anchor=anchor,
            )
        models = {}
        k_escape_lines = (0.0, 0.0, 0.0) + (0.0,) * len(
            K_ESCAPE_COMPONENTS)
        for name in COMPONENTS:
            if name == "k_alpha_escape":
                marker = next(lines).split()
                if len(marker) != 4 or marker[0] != "k_escape_lines":
                    raise ValueError(
                        f"invalid detector K-escape line record: {marker}")
                support = next(lines).split()
                if (len(support) != 1 + len(K_ESCAPE_COMPONENTS)
                        or support[0] != "k_escape_support"):
                    raise ValueError(
                        f"invalid detector K-escape support record: {support}")
                k_escape_lines = tuple(
                    [float(value) for value in marker[1:]]
                    + [float(value) for value in support[1:]])
            record = next(lines).split()
            if len(record) < 10 or record[0] != name:
                raise ValueError(f"invalid detector component {name}")
            block_count = int(record[9])
            block_sizes = tuple(
                int(value) for value in record[10:10 + block_count])
            if len(record) != 10 + block_count:
                raise ValueError(f"invalid detector component {name}")
            model_basis_size = int(record[7])
            node_count = int(record[8])
            if block_sizes and sum(block_sizes) != model_basis_size:
                raise ValueError(
                    f"detector block sizes must sum to basis size: {name}")
            nodes = []
            for _ in range(node_count):
                node = next(lines).split()
                if (
                    len(node) != 5 + model_basis_size
                    or node[0] != f"{name}_node"
                ):
                    raise ValueError(f"invalid detector node for {name}")
                nodes.append(Node(
                    feature=np.asarray([
                        float(value) for value in node[1:4]]),
                    log_yield=float(node[4]),
                    shape=np.asarray([
                        float(value) for value in node[5:]]),
                ))
            models[name] = ComponentModel(
                yield_bandwidth=np.asarray([
                    float(value) for value in record[1:4]]),
                shape_bandwidth=np.asarray([
                    float(value) for value in record[4:7]]),
                basis_size=model_basis_size,
                nodes=nodes,
            )
        # The finite-distance depths (ceelo_detector_finite_distance),
        # kept verbatim so a rewritten payload carries them on.
        record = next(lines)
        fields = record.split()
        if len(fields) != 2 or fields[0] != "finite_distance":
            raise ValueError(f"{path}: missing the finite_distance record")
        rows = [next(lines) for _ in range(int(fields[1]))]
        globals()["FINITE_DISTANCE_TEXT"] = record + "".join(rows)
        if next(lines).strip() != "END":
            raise ValueError("missing detector payload END marker")
    return models, rate_surface, k_escape_lines


def _parse_anchor_record(lines: list[str], start: int
                         ) -> RateSurfaceAnchor | None:
    """Parse the optional anchor records following the coefficients line."""
    if start >= len(lines):
        return None
    flag = lines[start].split()
    if not flag or flag[0] != "rate_surface_anchor_enabled":
        return None
    if len(flag) != 2 or flag[1] not in {"0", "1"}:
        raise ValueError(
            f"invalid rate-surface anchor flag: {flag}")
    if flag[1] == "0":
        return None
    if start + 1 >= len(lines):
        raise ValueError("truncated rate-surface anchor record")
    record = lines[start + 1].split()
    if len(record) != 4 or record[0] != "rate_surface_anchor":
        raise ValueError("invalid rate-surface anchor record")
    # Its geometry is the payload's, adopted when the payload was read.
    return RateSurfaceAnchor(
        blend_keV=float(record[1]), blend_power=float(record[2]),
        azimuth_quadrature_points=int(record[3]))


def read_rate_surface_only(path: Path) -> RateSurface:
    """Load just the global detected-rate field from any detector payload."""
    lines = path.read_text().splitlines(keepends=True)
    # The anchor is evaluated on the payload's geometry, so the surface is
    # only reusable for the same detector.
    for index, line in enumerate(lines):
        values = line.split()
        if values and values[0] == "detector":
            count = int(lines[index + 1].split()[1])
            adopt_detector(
                values[1], "".join(lines[index + 2:index + 2 + count]),
                source=str(path))
            break
    else:
        raise ValueError(f"{path} records no detector geometry")
    for index, line in enumerate(lines):
        values = line.split()
        if not values or values[0] != "rate_surface":
            continue
        if (
            len(values) != 6
            or int(values[1]) != RATE_ENERGY_BASIS_SIZE
            or int(values[2]) != RATE_ANGLE_BASIS_SIZE
            or index + 1 >= len(lines)
        ):
            raise ValueError(f"invalid rate surface in {path}")
        coefficients = lines[index + 1].split()
        expected_size = (
            RATE_ENERGY_BASIS_SIZE * RATE_ANGLE_BASIS_SIZE)
        if (
            len(coefficients) != expected_size + 1
            or coefficients[0] != "rate_surface_coefficients"
        ):
            raise ValueError(
                f"invalid rate-surface coefficients in {path}")
        return RateSurface(
            energy_warp_alpha=float(values[3]),
            coefficients=np.asarray([
                float(value) for value in coefficients[1:]]),
            anchor=_parse_anchor_record(lines, index + 2),
        )
    raise ValueError(f"no detected-rate surface in {path}")


def write_metrics(
    path: Path,
    rows: list[dict],
) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def append_detected_yield_nodes(
    models: dict[str, ComponentModel],
    cells: list[Cell],
    minimum_counts: float,
) -> int:
    """Promote cells into the scalar detection-rate field only.

    This is useful when an adaptive campaign finds a rate landmark but the
    process-conditional pulse shapes already meet their gate. The nonzero
    shape attached to this atom-only component is a serialization placeholder
    and is not evaluated by the detector model. The appended nodes carry
    training values beside the field's solved coefficients
    (interpolate_yield_nodes), so the field then reproduces neither exactly;
    refit when the kernel-average rate matters (a payload with a detected-rate
    surface does not evaluate this field).
    """
    model = models["zero_deposit"]
    added = 0
    for cell in cells:
        detected_yield = max(
            1.0 - cell.yields["zero_deposit"], 1e-300)
        if detected_yield * cell.events < minimum_counts:
            continue
        if any(np.allclose(
            node.feature, cell.feature, rtol=0.0, atol=1e-14)
               for node in model.nodes):
            continue
        shape = np.zeros(model.basis_size)
        shape[0] = 1.0
        model.nodes.append(Node(
            feature=cell.feature,
            log_yield=math.log(
                detected_yield
                / max(1.0 - detected_yield, 1e-300)),
            shape=shape,
        ))
        added += 1
    model.nodes.sort(key=lambda node: tuple(node.feature))
    return added


def bandwidth_values(text: str) -> list[np.ndarray]:
    answer = []
    for candidate in text.split(";"):
        values = np.asarray([
            float(value) for value in candidate.split(",")])
        if values.size != 3 or np.any(values <= 0.0):
            raise argparse.ArgumentTypeError(
                "bandwidths must be semicolon-separated triples")
        answer.append(values)
    return answer


def component_scale(text: str) -> tuple[str, float]:
    name, separator, scale_text = text.partition("=")
    if separator != "=" or name not in COMPONENTS:
        raise argparse.ArgumentTypeError(
            "component scale must be COMPONENT=FACTOR")
    try:
        scale = float(scale_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "component scale factor must be numeric") from error
    if not math.isfinite(scale) or scale <= 0.0:
        raise argparse.ArgumentTypeError(
            "component scale factor must be positive and finite")
    return name, scale


def selection_score(report: dict) -> float:
    return (
        report["detected_yield_relative_error"]["p95"]
        + report["pulse_shape_total_variation"]["p95"]
        + report["continuum_shape_total_variation"]["p95"])


def main() -> int:
    global BASIS_SIZE, SHAPE_ROUGHNESS
    global SHAPE_ROUGHNESS_REFERENCE_COUNTS
    global COMPONENT_BASIS_SIZES
    global PAIR_YIELD_BANDWIDTH_SCALE, PAIR_SHAPE_BANDWIDTH_SCALE
    parser = argparse.ArgumentParser()
    parser.add_argument("training", type=Path, nargs="?")
    parser.add_argument(
        "--additional-training", action="append", type=Path, default=[])
    parser.add_argument("--validation", type=Path)
    parser.add_argument(
        "--evaluate-payload", type=Path,
        help="audit an existing payload without refitting")
    parser.add_argument(
        "--augment-detected-yield", action="store_true",
        help=(
            "with --evaluate-payload, promote --additional-training cells "
            "into only the detected-yield field and write an updated payload"))
    parser.add_argument(
        "--fit-detected-rate-surface", action="store_true",
        help=(
            f"fit the {RATE_ENERGY_BASIS_SIZE}x{RATE_ANGLE_BASIS_SIZE} "
            "Bernstein log-odds surface; with "
            "--evaluate-payload this replaces only the total-rate field"))
    parser.add_argument(
        "--rate-energy-warp-alpha", type=float,
        default=RATE_ENERGY_WARP_ALPHA,
        help=(
            "passive-proxy stretch of the rate-surface energy coordinate. "
            "The default suits anchored fits; un-anchored fits should "
            "use 3"))
    parser.add_argument(
        "--rate-smoothing-ridge", type=float,
        default=RATE_SMOOTHING_RIDGE,
        help="second-difference penalty weight for the rate surface")
    parser.add_argument(
        "--rate-high-energy-ridge-kappa", type=float,
        default=RATE_HIGH_ENERGY_RIDGE_KAPPA,
        help=(
            "extra smoothing on rate-surface energy coefficients above the "
            "cut; 0 reproduces an unpenalized fit"))
    parser.add_argument(
        "--rate-high-energy-ridge-cut", type=float,
        default=RATE_HIGH_ENERGY_RIDGE_CUT,
        help="energy-coordinate fraction above which the extra smoothing acts")
    parser.add_argument(
        "--rate-surface-anchor", action="store_true",
        help=(
            "anchor the detected-rate log odds on the deterministic "
            "surface-entry transmission of the passive layers (the "
            "sub-60 keV turn-on the Bernstein surface cannot interpolate "
            "from sparse training energies), evaluated on the geometry the "
            "banks record."))
    parser.add_argument(
        "--rate-surface-anchor-executable", type=Path,
        default=Path(
            "build-rel/tools/full_spectrum/ceelo_rate_surface_anchor"),
        help=(
            "helper that evaluates the anchor and derives the landmark "
            "(shares the C++ evaluator)"))
    parser.add_argument(
        "--rate-surface-anchor-blend-kev", type=float, default=150.0)
    parser.add_argument(
        "--rate-surface-anchor-blend-power", type=float, default=2.0)
    parser.add_argument(
        "--detected-rate-training", action="append", type=Path, default=[],
        help="campaign used only to fit the detected-rate surface")
    parser.add_argument(
        "--detected-rate-payload", type=Path,
        help=(
            "reuse the independently selected detected-rate surface from "
            "another detector payload"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--yield-bandwidths", type=bandwidth_values,
        default=bandwidth_values(
            "0.035,0.10,0.10;0.05,0.15,0.15;"
            "0.07,0.20,0.20;0.10,0.30,0.30"))
    parser.add_argument(
        "--shape-bandwidths", type=bandwidth_values,
        default=bandwidth_values(
            "0.035,0.10,0.10;0.05,0.15,0.15;"
            "0.07,0.20,0.20;0.10,0.30,0.30"))
    parser.add_argument("--minimum-component-counts", type=float, default=25.0)
    parser.add_argument(
        "--basis-size", type=int, default=BASIS_SIZE,
        help=(
            "cubic B-spline terms per continuum component; values below the "
            "compiled production size are intended for model-selection "
            "ablations"))
    parser.add_argument(
        "--shape-roughness", type=float, default=SHAPE_ROUGHNESS,
        help=(
            "second-difference ridge applied to each continuum node; it is "
            "automatically strengthened for statistically sparse components"))
    parser.add_argument(
        "--shape-roughness-reference-counts", type=float,
        default=SHAPE_ROUGHNESS_REFERENCE_COUNTS)
    parser.add_argument(
        "--pair-yield-bandwidth-scale", type=float,
        default=PAIR_YIELD_BANDWIDTH_SCALE,
        help="scale all RBF bandwidth coordinates for pair-history rates")
    parser.add_argument(
        "--pair-shape-bandwidth-scale", type=float,
        default=PAIR_SHAPE_BANDWIDTH_SCALE,
        help="scale all RBF bandwidth coordinates for pair-history shapes")
    parser.add_argument(
        "--rescale-loaded-pair-yield-bandwidth", type=float, default=1.0,
        help=(
            "with --evaluate-payload, multiply existing pair-component "
            "yield bandwidths without refitting node coefficients"))
    parser.add_argument(
        "--rescale-loaded-pair-shape-bandwidth", type=float, default=1.0,
        help=(
            "with --evaluate-payload, multiply existing pair-component "
            "shape bandwidths without refitting node coefficients"))
    parser.add_argument(
        "--component-yield-bandwidth-scale",
        action="append", type=component_scale, default=[],
        metavar="COMPONENT=FACTOR",
        help=(
            "when fitting, multiply one component's yield bandwidth (after "
            "the pair scale); repeat per component"))
    parser.add_argument(
        "--yield-node-iterations", type=int, default=YIELD_NODE_ITERATIONS,
        help=(
            "residual iterations that make each yield field reproduce its "
            "training nodes (0: the plain kernel average of the training "
            "values)"))
    parser.add_argument(
        "--rescale-loaded-component-yield-bandwidth",
        action="append", type=component_scale, default=[],
        metavar="COMPONENT=FACTOR",
        help=(
            "with --evaluate-payload, additionally multiply one component's "
            "yield bandwidth; repeat for independently selected pair "
            "fractions"))
    parser.add_argument(
        "--photoelectric-basis-size", type=int,
        help="override the photoelectric-escape continuum basis size")
    parser.add_argument(
        "--single-compton-basis-size", type=int,
        help="override the single-Compton continuum basis size")
    parser.add_argument(
        "--multiple-compton-basis-size", type=int,
        help="override the multiple-Compton continuum basis size")
    parser.add_argument(
        "--pair-continuum-basis-size", type=int,
        help=(
            "override all three hierarchical pair-continuum basis sizes"))
    parser.add_argument(
        "--pair-no-compton-basis-size", type=int)
    parser.add_argument(
        "--pair-single-compton-basis-size", type=int)
    parser.add_argument(
        "--pair-multiple-compton-basis-size", type=int)
    parser.add_argument(
        "--payload-precision", type=int, default=12,
        help=(
            "significant decimal digits in the portable text payload; "
            "12 is well below the Monte Carlo/model uncertainty while "
            "keeping the selected payload compact"))
    args = parser.parse_args()
    globals()["RATE_ENERGY_WARP_ALPHA"] = args.rate_energy_warp_alpha
    globals()["RATE_SMOOTHING_RIDGE"] = args.rate_smoothing_ridge
    globals()["RATE_HIGH_ENERGY_RIDGE_KAPPA"] = (
        args.rate_high_energy_ridge_kappa)
    globals()["RATE_HIGH_ENERGY_RIDGE_CUT"] = args.rate_high_energy_ridge_cut
    if args.rate_surface_anchor_executable.is_file():
        globals()["ANCHOR_EXECUTABLE"] = args.rate_surface_anchor_executable
    rate_surface_anchor = None
    if args.rate_surface_anchor:
        rate_surface_anchor = RateSurfaceAnchor(
            executable=args.rate_surface_anchor_executable,
            blend_keV=args.rate_surface_anchor_blend_kev,
            blend_power=args.rate_surface_anchor_blend_power)
    if args.basis_size < 4:
        parser.error("--basis-size must be at least 4")
    if args.shape_roughness < 0.0:
        parser.error("--shape-roughness cannot be negative")
    if args.shape_roughness_reference_counts <= 0.0:
        parser.error(
            "--shape-roughness-reference-counts must be positive")
    if args.pair_yield_bandwidth_scale <= 0.0:
        parser.error("--pair-yield-bandwidth-scale must be positive")
    if args.yield_node_iterations < 0:
        parser.error("--yield-node-iterations cannot be negative")
    if args.pair_shape_bandwidth_scale <= 0.0:
        parser.error("--pair-shape-bandwidth-scale must be positive")
    if not 9 <= args.payload_precision <= 17:
        parser.error("--payload-precision must be between 9 and 17")
    if args.rescale_loaded_pair_yield_bandwidth <= 0.0:
        parser.error(
            "--rescale-loaded-pair-yield-bandwidth must be positive")
    if args.rescale_loaded_pair_shape_bandwidth <= 0.0:
        parser.error(
            "--rescale-loaded-pair-shape-bandwidth must be positive")
    BASIS_SIZE = args.basis_size
    SHAPE_ROUGHNESS = args.shape_roughness
    SHAPE_ROUGHNESS_REFERENCE_COUNTS = (
        args.shape_roughness_reference_counts)
    PAIR_YIELD_BANDWIDTH_SCALE = args.pair_yield_bandwidth_scale
    PAIR_SHAPE_BANDWIDTH_SCALE = args.pair_shape_bandwidth_scale
    basis_helpers.BASIS_SIZE = BASIS_SIZE
    COMPONENT_BASIS_SIZES = {
        name: (BASIS_SIZE if name in CONTINUOUS else 1)
        for name in COMPONENTS
    }
    for name, value in (
        ("photoelectric_escape", args.photoelectric_basis_size),
        ("single_compton", args.single_compton_basis_size),
        ("multiple_compton", args.multiple_compton_basis_size),
        ("pair_no_compton",
         args.pair_no_compton_basis_size
         or args.pair_continuum_basis_size),
        ("pair_single_compton",
         args.pair_single_compton_basis_size
         or args.pair_continuum_basis_size),
        ("pair_multiple_compton",
         args.pair_multiple_compton_basis_size
         or args.pair_continuum_basis_size),
    ):
        if value is not None:
            if value < 4:
                parser.error(
                    f"--{name.replace('_', '-')}-basis-size must be "
                    "at least 4")
            COMPONENT_BASIS_SIZES[name] = value
    _BASIS_CACHE.clear()
    _SHAPE_CACHE.clear()

    args.output.mkdir(parents=True, exist_ok=True)
    if args.evaluate_payload:
        if args.validation is None:
            parser.error("--evaluate-payload requires --validation")
        models, rate_surface, loaded_k_escape_lines = read_payload(
            args.evaluate_payload)
        for name in PAIR_COMPONENTS:
            models[name].yield_bandwidth *= (
                args.rescale_loaded_pair_yield_bandwidth)
            models[name].shape_bandwidth *= (
                args.rescale_loaded_pair_shape_bandwidth)
        for name, scale in args.rescale_loaded_component_yield_bandwidth:
            models[name].yield_bandwidth *= scale
        if args.detected_rate_payload:
            rate_surface = read_rate_surface_only(
                args.detected_rate_payload)
        added_nodes = 0
        payload_training_cells = []
        if args.augment_detected_yield:
            promoted_cells = []
            for campaign in args.additional_training:
                promoted_cells.extend(load_campaign(campaign))
            if not promoted_cells:
                parser.error(
                    "--augment-detected-yield requires "
                    "--additional-training")
            added_nodes = append_detected_yield_nodes(
                models, promoted_cells, args.minimum_component_counts)
            payload_training_cells.extend(promoted_cells)
        rate_training_cell_count = 0
        rate_training_populated_cell_count = 0
        if (
            args.fit_detected_rate_surface
            and args.detected_rate_payload is None
        ):
            rate_training_cells = []
            for campaign in args.detected_rate_training:
                rate_training_cells.extend(load_campaign(campaign))
            if not rate_training_cells:
                parser.error(
                    "--fit-detected-rate-surface requires "
                    "--detected-rate-training")
            rate_surface = fit_rate_surface(
                rate_training_cells, rate_surface_anchor)
            rate_training_cell_count = len(rate_training_cells)
            rate_training_populated_cell_count = sum(
                (1.0 - cell.yields["zero_deposit"]) * cell.events
                >= 25.0
                for cell in rate_training_cells)
            payload_training_cells.extend(rate_training_cells)
        if (
            args.augment_detected_yield
            or args.fit_detected_rate_surface
            or args.detected_rate_payload
            or args.rescale_loaded_pair_yield_bandwidth != 1.0
            or args.rescale_loaded_pair_shape_bandwidth != 1.0
            or args.rescale_loaded_component_yield_bandwidth
        ):
            write_payload(
                args.output / "detector_response_payload.txt",
                models, payload_training_cells,
                validity_keV=(10.0, 12000.0),
                rate_surface=rate_surface,
                precision=args.payload_precision,
                k_escape_lines=loaded_k_escape_lines)
        validation = load_campaign(args.validation)
        validation_report, validation_rows = metrics(
            validation, models, rate_surface, loaded_k_escape_lines)
        report = {
            "model": (
                "bernstein-detected-rate-plus-"
                "positive-physics-coordinate-bspline-rbf"),
            "component_basis_sizes": {
                name: model.basis_size
                for name, model in models.items()
            },
            "payload_precision_significant_digits":
                args.payload_precision,
            "payload": str(args.evaluate_payload),
            "validation_cell_count": len(validation),
            "added_detected_yield_nodes": added_nodes,
            "component_yield_bandwidth_rescaling": dict(
                args.rescale_loaded_component_yield_bandwidth),
            "detected_rate_training_cell_count":
                rate_training_cell_count,
            "detected_rate_populated_cell_count":
                rate_training_populated_cell_count,
            "validation": validation_report,
        }
        (args.output / "report.json").write_text(
            json.dumps(report, indent=2) + "\n")
        write_metrics(
            args.output / "validation_metrics.csv", validation_rows)
        print(json.dumps(validation_report, indent=2))
        return 0
    if args.training is None:
        parser.error("training campaign is required when fitting")

    training = load_campaign(args.training)
    for campaign in args.additional_training:
        training.extend(load_campaign(campaign))
    training = prefer_highest_statistics(training)
    validation = (
        load_campaign(args.validation) if args.validation else training)
    rate_training = list(training)
    for campaign in args.detected_rate_training:
        rate_training.extend(load_campaign(campaign))
    rate_training = prefer_highest_statistics(rate_training)
    rate_surface = (
        read_rate_surface_only(args.detected_rate_payload)
        if args.detected_rate_payload
        else fit_rate_surface(rate_training, rate_surface_anchor)
        if args.fit_detected_rate_surface
        else None)
    # One K-escape line/support record for the whole run: the fit, the
    # metrics and the written payload must all gate identically.
    fitted_k_escape_lines = k_escape_line_energies(
        training, args.minimum_component_counts)
    candidates = []
    best = None
    for yield_bandwidth in args.yield_bandwidths:
        for shape_bandwidth in args.shape_bandwidths:
            models = fit_components(
                training, yield_bandwidth, shape_bandwidth,
                args.minimum_component_counts,
                dict(args.component_yield_bandwidth_scale),
                args.yield_node_iterations)
            report, _ = metrics(
                validation, models, rate_surface, fitted_k_escape_lines)
            candidate = {
                "yield_bandwidth": yield_bandwidth.tolist(),
                "shape_bandwidth": shape_bandwidth.tolist(),
                "component_yield_bandwidth_scale": dict(
                    args.component_yield_bandwidth_scale),
                "yield_node_iterations": args.yield_node_iterations,
                "validation": report,
                "score": selection_score(report),
            }
            candidates.append(candidate)
            if best is None or candidate["score"] < best[0]:
                best = (candidate["score"], models, candidate)
    assert best is not None
    selected_models = best[1]
    train_report, train_rows = metrics(
        training, selected_models, rate_surface, fitted_k_escape_lines)
    validation_report, validation_rows = metrics(
        validation, selected_models, rate_surface, fitted_k_escape_lines)
    report = {
        "model": "positive-physics-coordinate-bspline-rbf",
        "component_basis_sizes": {
            name: model.basis_size
            for name, model in selected_models.items()
        },
        "payload_precision_significant_digits": args.payload_precision,
        "training_cell_count": len(training),
        "detected_rate_training_cell_count": len(rate_training),
        "validation_cell_count": len(validation),
        "selected": best[2],
        "training": train_report,
        "validation": validation_report,
        "candidates": candidates,
    }
    write_payload(
        args.output / "detector_response_payload.txt",
        selected_models, training, rate_surface=rate_surface,
        precision=args.payload_precision,
        k_escape_lines=fitted_k_escape_lines)
    (args.output / "report.json").write_text(
        json.dumps(report, indent=2) + "\n")
    for name, rows in (
        ("training_metrics.csv", train_rows),
        ("validation_metrics.csv", validation_rows),
    ):
        write_metrics(args.output / name, rows)
    print(json.dumps(report["selected"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
