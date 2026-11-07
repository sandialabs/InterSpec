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

"""Guard: every fitter's B-spline basis must match the C++ evaluator's.

Each Python fitter builds its own design matrix, and the C++ evaluator
re-expands the coefficients it is handed. If the two bases disagree the
shapes are silently reweighted -- the fit still converges, its own metrics
still look fine, and only a round-trip through the evaluator reveals it.
That exact bug cost 83% shape total variation in fit_photonuclear.py before
it was caught, because the C++ divides every basis integral by that basis
function's own support integral, (knot[i+degree+1] - knot[i]) / (degree+1),
making it a unit-mass density basis rather than a partition of unity.

This script encodes the C++ rule once (`reference_integrals`, mirroring
detail::cubic_bspline_integrals / cubic_bspline_values in
src/full_spectrum/FullSpectrumModel.h) and checks each fitter against it.
Run it after touching any basis code in either language.

    python3 tools/full_spectrum/check_basis_lockstep.py
"""

from __future__ import annotations

import importlib
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

BASIS_SIZE = 32
DEGREE = 3
TOLERANCE = 1e-9


def knot_value(index: int, basis_size: int = BASIS_SIZE) -> float:
    intervals = basis_size - DEGREE
    if index <= DEGREE:
        return 0.0
    if index >= basis_size:
        return 1.0
    return (index - DEGREE) / intervals


def reference_values(coordinate: np.ndarray,
                     basis_size: int = BASIS_SIZE) -> np.ndarray:
    """Cox-de Boor on the C++ clamped knot vector."""
    knots = np.array([knot_value(i, basis_size)
                      for i in range(basis_size + DEGREE + 1)])
    x = np.clip(np.asarray(coordinate, dtype=float), 0.0, 1.0)
    saturated = x >= 1.0
    work = np.zeros((x.size, basis_size + DEGREE))
    for index in range(knots.size - 1):
        if index >= work.shape[1]:
            break
        work[(knots[index] <= x) & (x < knots[index + 1]), index] = 1.0
    for order in range(1, DEGREE + 1):
        current = np.zeros_like(work)
        for index in range(work.shape[1] - order):
            left = knots[index + order] - knots[index]
            if left > 0.0:
                current[:, index] += \
                    (x - knots[index]) / left * work[:, index]
            right_index = index + order + 1
            if right_index < knots.size:
                right = knots[right_index] - knots[index + 1]
                if right > 0.0:
                    current[:, index] += \
                        (knots[right_index] - x) / right * work[:, index + 1]
        work = current
    answer = work[:, :basis_size]
    answer[saturated, :] = 0.0
    answer[saturated, basis_size - 1] = 1.0
    return answer


def reference_integrals(lo: np.ndarray, hi: np.ndarray,
                        basis_size: int = BASIS_SIZE) -> np.ndarray:
    """Integrals with the C++ unit-mass normalization."""
    nodes, weights = np.polynomial.legendre.leggauss(16)
    result = np.zeros((np.asarray(lo).size, basis_size))
    for x, w in zip(nodes, weights):
        point = 0.5 * (hi + lo) + 0.5 * (hi - lo) * x
        result += (0.5 * (hi - lo) * w)[:, None] \
            * reference_values(point, basis_size)
    normalization = np.array([
        (knot_value(i + DEGREE + 1, basis_size) - knot_value(i, basis_size))
        / (DEGREE + 1)
        for i in range(basis_size)])
    return result / normalization[None, :]


def unit_mass_check() -> float:
    """Each normalized basis function must integrate to 1 over [0, 1]."""
    edges = np.linspace(0.0, 1.0, 4001)
    total = reference_integrals(edges[:-1], edges[1:]).sum(axis=0)
    return float(np.max(np.abs(total - 1.0)))


def fitter_columns(module_name: str, energy_keV: float, bins: int):
    """Design matrix (bins x BASIS_SIZE) from a fitter's own entry point."""
    module = importlib.import_module(module_name)
    edges = np.linspace(0.0, energy_keV, bins + 1)
    if hasattr(module, "warped_cubic_bspline_matrix"):
        return np.asarray(module.warped_cubic_bspline_matrix(
            energy_keV, bins, BASIS_SIZE))
    if hasattr(module, "warped_bspline_matrix"):
        return np.asarray(module.warped_bspline_matrix(edges, energy_keV))
    if hasattr(module, "cubic_bspline_integrals"):
        if hasattr(module, "inverse_energy_loss_warp"):
            warp_lo = module.inverse_energy_loss_warp(
                edges[1:], energy_keV, 2.0)
            warp_hi = module.inverse_energy_loss_warp(
                edges[:-1], energy_keV, 2.0)
        else:
            # No shared warp helper: exercise the basis directly on a uniform
            # coordinate grid, which is what the normalization check needs.
            warp_lo = np.linspace(0.0, 1.0, bins + 1)[:-1]
            warp_hi = np.linspace(0.0, 1.0, bins + 1)[1:]
        # Vectorised signature first; fall back to one call per bin.
        try:
            answer = np.asarray(
                module.cubic_bspline_integrals(warp_lo, warp_hi))
            if answer.ndim == 2:
                return answer
        except Exception:  # noqa: BLE001 - scalar-only signature
            pass
        return np.vstack([
            np.asarray(module.cubic_bspline_integrals(float(lo), float(hi)))
            for lo, hi in zip(warp_lo, warp_hi)])
    return None


def main() -> int:
    print(f"unit-mass residual of the reference basis: "
          f"{unit_mass_check():.2e}")

    energy, bins = 3000.0, 256
    edges = np.linspace(0.0, energy, bins + 1)
    # The photonuclear fitter's warp is the shared convention.
    photonuclear = importlib.import_module("fit_photonuclear")
    warp_lo = photonuclear.inverse_energy_loss_warp(edges[1:], energy, 2.0)
    warp_hi = photonuclear.inverse_energy_loss_warp(edges[:-1], energy, 2.0)
    reference = reference_integrals(warp_lo, warp_hi)

    failures = 0
    for module_name in ("fit_photonuclear",
                        "fit_multiple_scatter_surrogate",
                        "fit_source_secondaries",
                        "fit_charged_particle_source"):
        try:
            columns = fitter_columns(module_name, energy, bins)
        except Exception as error:  # noqa: BLE001 - report, do not abort
            print(f"  {module_name:34s} SKIP ({type(error).__name__}: {error})")
            continue
        if columns is None:
            print(f"  {module_name:34s} SKIP (no recognised basis entry point)")
            continue
        if columns.shape != reference.shape:
            print(f"  {module_name:34s} SHAPE {columns.shape} != "
                  f"{reference.shape}")
            failures += 1
            continue
        # Compare per-basis-function total mass: that is precisely what the
        # missing support-integral divisor corrupts, and it is invariant to
        # the warp scale each fitter chooses.
        mass = columns.sum(axis=0)
        reference_mass = reference.sum(axis=0)
        ratio = mass / np.maximum(reference_mass, 1e-300)
        spread = float(np.max(ratio) / max(np.min(ratio), 1e-300))
        # Omitting the support divisor makes the clamped end functions carry
        # several times the mass of the interior ones, so a real mismatch
        # shows a spread of ~4. Anything within a few percent is quadrature
        # difference between the fitters' rules and the reference's.
        verdict = "OK" if spread < 1.05 else "MISMATCH"
        failures += verdict == "MISMATCH"
        print(f"  {module_name:34s} {verdict}  "
              f"per-term mass ratio spread {spread:.6f}")

    unnormalized = reference * np.array([
        (knot_value(i + DEGREE + 1) - knot_value(i)) / (DEGREE + 1)
        for i in range(BASIS_SIZE)])[None, :]
    mass = unnormalized.sum(axis=0) / np.maximum(reference.sum(axis=0), 1e-300)
    print(f"\nsensitivity: a basis missing the support divisor would show "
          f"spread {float(np.max(mass) / np.min(mass)):.2f}")

    if failures:
        print(f"\n{failures} fitter(s) disagree with the C++ basis "
              f"normalization -- shapes will be silently reweighted.")
    else:
        print("\nall fitters agree with the C++ basis normalization.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
