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

"""Fit the escaping-photonuclear (E, Z, tau) yield and shape surrogate.

CeeLo transport has no photonuclear physics, so this component's truth is
Geant4-only (the bank written by run_photonuclear_campaign.py). The payload keeps
the mechanistic removal baseline `1 - exp(-tau_pn)` from the Geant4
attenuation table and this fit supplies the positive conditional correction
plus the emitted shape:

    yield(E, Z, tau) = (1 - exp(-tau_pn)) * exp(log_conditional_yield)
    shape            = positive cubic B-splines over the inverse-energy-loss
                       warp, normalized to unit area

Both fields are normalized-Gaussian RBF interpolants over the same three
features the other source surrogates use (log-energy warp, Z/92, log1p(tau)).
That choice is deliberate: outside the trained hull a normalized-Gaussian
weight renormalizes to the nearest nodes, so the field FREEZES rather than
extrapolating polynomially -- the failure mode that produced the deep-tau
exp() overflow in the multiple-scatter tau-yield model.

Validation is whole-cell: Z holdouts never seen in training, and (when the
off-grid bank is present) cells at optical depths and Z values BETWEEN the
training nodes. On-grid holdouts alone hide interpolation error -- that is
how the mid-tau defect in the secondary surrogate went unnoticed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import statistics
import sys

import numpy as np

BASIS_SIZE = 32
ELECTRON_MASS_KEV = 510.99895


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bank", type=pathlib.Path, required=True,
        help="training cells written by run_photonuclear_campaign.py")
    parser.add_argument(
        "--offgrid", type=pathlib.Path,
        help="optional off-grid (tau, Z) validation cells")
    parser.add_argument(
        "--photonuclear-xs", type=pathlib.Path, required=True,
        help="natural-element photonuclear attenuation table "
             "(g4_photonuclear_xs output)")
    parser.add_argument(
        "--output", type=pathlib.Path, required=True)
    parser.add_argument("--holdout-z", type=int, nargs="+", default=[47, 73])
    parser.add_argument(
        "--warp-scale", type=float, nargs="+",
        default=[0.5, 1.0, 2.0, 4.0, 8.0],
        help="inverse-energy-loss warp scales to scan")
    parser.add_argument(
        "--yield-bandwidth", type=float, nargs=3, default=None,
        help="override the scanned yield bandwidths")
    parser.add_argument(
        "--shape-bandwidth", type=float, nargs=3, default=None)
    parser.add_argument(
        "--escape-kappa", type=float, default=0.8,
        help="secondary/primary mass-attenuation ratio in the escape factor")
    return parser.parse_args()


def escape_factor(tau: float, kappa: float) -> float:
    """Mirror of detail::photonuclear_escape_factor."""
    x = (1.0 - kappa) * tau
    if abs(x) < 1e-6:
        ratio = 1.0 - x / 2.0 + x * x / 6.0
    else:
        ratio = (1.0 - math.exp(-x)) / x
    return math.exp(-kappa * tau) * ratio


# ---- features (must match SourceGeomFluxModel::evaluate) --------------------

def feature(energy_keV: float, z_value: float, tau: float) -> np.ndarray:
    return np.array([
        math.log(energy_keV / 10.0) / math.log(1200.0),
        z_value / 92.0,
        math.log1p(tau) / math.log(25.0),
    ])


def inverse_energy_loss_warp(
        energy_keV: np.ndarray, incident_keV: float,
        scale: float) -> np.ndarray:
    """Mirror of detail::inverse_energy_loss_warp."""
    safe = np.maximum(energy_keV, 1e-9)
    loss = 0.5 * ELECTRON_MASS_KEV * (1.0 / safe - 1.0 / incident_keV)
    return np.clip(loss / (scale + loss), 0.0, 1.0)


def cubic_bspline_integrals(
        warp_lo: np.ndarray, warp_hi: np.ndarray) -> np.ndarray:
    """Integrals of the clamped cubic B-spline basis, matching the C++ basis.

    Evaluated by dense quadrature: this is fitting-side only, and the
    round-trip against the evaluator is checked by the caller.
    """
    nodes = 8
    gl_x, gl_w = np.polynomial.legendre.leggauss(nodes)
    result = np.zeros((warp_lo.size, BASIS_SIZE))
    for x, w in zip(gl_x, gl_w):
        point = 0.5 * (warp_hi + warp_lo) + 0.5 * (warp_hi - warp_lo) * x
        result += (0.5 * (warp_hi - warp_lo) * w)[:, None] \
            * cubic_bspline_values(point)
    # The C++ basis divides each integral by the basis function's own support
    # integral, (knot[i+degree+1] - knot[i]) / (degree+1), so every term
    # carries unit mass over [0, 1] -- a density basis, not a partition of
    # unity. The fitter MUST mirror this or the shapes it emits are silently
    # reweighted by the evaluator.
    return result / _basis_normalization()[None, :]


def _knot_value(index: int) -> float:
    degree = 3
    intervals = BASIS_SIZE - degree
    if index <= degree:
        return 0.0
    if index >= BASIS_SIZE:
        return 1.0
    return (index - degree) / intervals


def _basis_normalization() -> np.ndarray:
    degree = 3
    return np.array([
        (_knot_value(i + degree + 1) - _knot_value(i)) / (degree + 1)
        for i in range(BASIS_SIZE)])


def cubic_bspline_values(coordinate: np.ndarray) -> np.ndarray:
    """Clamped uniform cubic B-spline basis on [0, 1] with BASIS_SIZE terms."""
    degree = 3
    interior = BASIS_SIZE - degree - 1
    knots = np.concatenate([
        np.zeros(degree + 1),
        np.arange(1, interior + 1) / (interior + 1),
        np.ones(degree + 1)])
    x = np.clip(coordinate, 0.0, 1.0)
    saturated = x >= 1.0
    work = np.zeros((x.size, BASIS_SIZE + degree))
    for index in range(knots.size - 1):
        if index >= work.shape[1]:
            break
        inside = (knots[index] <= x) & (x < knots[index + 1])
        work[inside, index] = 1.0
    for order in range(1, degree + 1):
        current = np.zeros_like(work)
        for index in range(work.shape[1] - order):
            left_den = knots[index + order] - knots[index] \
                if index + order < knots.size else 0.0
            if left_den > 0.0:
                current[:, index] += (x - knots[index]) / left_den \
                    * work[:, index]
            right_index = index + order + 1
            if right_index < knots.size:
                right_den = knots[right_index] - knots[index + 1]
                if right_den > 0.0:
                    current[:, index] += \
                        (knots[right_index] - x) / right_den \
                        * work[:, index + 1]
        work = current
    answer = work[:, :BASIS_SIZE]
    # Match the C++ early return: at the clamped upper end the last basis
    # function carries everything.
    answer[saturated, :] = 0.0
    answer[saturated, BASIS_SIZE - 1] = 1.0
    return answer


# ---- bank loading ----------------------------------------------------------

class Cell:
    __slots__ = ("energy", "z", "tau", "areal", "yield_", "sigma",
                 "edges", "counts", "variance", "identity", "events",
                 "escaping")


def load_cell(path: pathlib.Path) -> Cell | None:
    lo_list: list[float] = []
    hi_list: list[float] = []
    values: list[float] = []
    variances: list[float] = []
    total = 0.0
    total_var = 0.0
    energy = z = tau = areal = None
    escaping = 0.0
    with open(path) as stream:
        for row in csv.DictReader(stream):
            # Every escaping photon, uncollided included: the denominator of
            # the photonuclear share of the escaping spectrum.
            if row["component"] != "absorbed":
                escaping += float(row["sum_w"]) / float(row["events"])
            if row["component"] != "photonuclear_secondary":
                continue
            if energy is None:
                energy = float(row["energy_keV"])
                z = int(row["Z"])
                tau = float(row["tau"])
                areal = float(row["areal_density_g_cm2"])
            events = float(row["events"])
            weight = float(row["sum_w"]) / events
            variance = float(row["sum_w2"]) / (events * events)
            total += weight
            total_var += variance
            lower = float(row["bin_lo_keV"])
            upper = float(row["bin_hi_keV"])
            if upper > lower:
                lo_list.append(lower)
                hi_list.append(upper)
                values.append(weight)
                variances.append(variance)
    if energy is None or total <= 0.0:
        return None
    cell = Cell()
    cell.energy, cell.z, cell.tau, cell.areal = energy, z, tau, areal
    cell.events = events
    cell.yield_, cell.sigma = total, math.sqrt(total_var)
    cell.edges = (np.asarray(lo_list), np.asarray(hi_list))
    cell.counts = np.asarray(values)
    cell.variance = np.asarray(variances)
    cell.identity = path.stem[: -len("_geant4")]
    cell.escaping = escaping
    return cell


def load_photonuclear_mu(path: pathlib.Path) -> dict:
    table: dict[float, dict[int, float]] = {}
    with open(path) as stream:
        for row in csv.DictReader(stream):
            table.setdefault(float(row["energy_keV"]), {})[
                int(row["Z"])] = float(row["mu_over_rho_cm2_g"])
    return table


def mu_photonuclear(table: dict, energy_keV: float, z_value: int) -> float:
    energies = sorted(table)
    if energy_keV <= energies[0]:
        return table[energies[0]].get(z_value, 0.0)
    if energy_keV >= energies[-1]:
        return table[energies[-1]].get(z_value, 0.0)
    for lower, upper in zip(energies, energies[1:]):
        if lower <= energy_keV <= upper:
            f = (energy_keV - lower) / (upper - lower)
            return (1.0 - f) * table[lower].get(z_value, 0.0) \
                + f * table[upper].get(z_value, 0.0)
    return 0.0


# ---- per-cell shape fit ----------------------------------------------------

def fit_cell_shape(cell: Cell, warp_scale: float) -> np.ndarray:
    """Nonnegative least squares of the cell spectrum onto the warped basis."""
    lo, hi = cell.edges
    # The warp is decreasing in energy, so the basis integral runs hi -> lo.
    warp_lo = inverse_energy_loss_warp(hi, cell.energy, warp_scale)
    warp_hi = inverse_energy_loss_warp(lo, cell.energy, warp_scale)
    design = cubic_bspline_integrals(warp_lo, warp_hi)
    # Poisson chi-square weights. An EMPTY bin has zero sampled variance, so
    # floor the variance at one count: without the floor a zero-count bin
    # takes the largest weight in the fit -- the opposite of its information
    # content -- and drags the shape toward the sparse tail.
    one_count = 1.0 / (cell.events * cell.events)
    weights = 1.0 / np.sqrt(np.maximum(cell.variance, one_count))
    a = design * weights[:, None]
    b = cell.counts * weights
    coefficients = nnls(a, b)
    total = coefficients.sum()
    if total <= 0.0:
        coefficients = np.ones(BASIS_SIZE)
        total = coefficients.sum()
    return coefficients / total


def reconstruction_tv(cell: Cell, coefficients: np.ndarray,
                      warp_scale: float) -> float:
    """Total variation between the fitted shape and the truth histogram.

    This is the check that the BASIS can represent the cell at all. Comparing
    one fitted coefficient vector against another (the interpolation metric)
    is blind to a warp scale that collapses the whole spectrum into a corner
    of the basis domain: both fits are then equally wrong and agree perfectly.
    """
    lo, hi = cell.edges
    warp_lo = inverse_energy_loss_warp(hi, cell.energy, warp_scale)
    warp_hi = inverse_energy_loss_warp(lo, cell.energy, warp_scale)
    predicted = cubic_bspline_integrals(warp_lo, warp_hi) @ coefficients
    truth_total = cell.counts.sum()
    predicted_total = predicted.sum()
    if truth_total <= 0.0 or predicted_total <= 0.0:
        return 1.0
    return 0.5 * float(np.abs(
        predicted / predicted_total - cell.counts / truth_total).sum())


def nnls(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Lawson-Hanson active-set NNLS (scipy is not a dependency here).

    Exact for this problem size (32 unknowns), which matters: an iterative
    projected-gradient solve was still 10 percentage points of shape total
    variation away from the optimum after 50k iterations, and that slack
    reads as basis inadequacy when it is really solver error.
    """
    n = a.shape[1]
    ata = a.T @ a
    atb = a.T @ b
    x = np.zeros(n)
    passive = np.zeros(n, dtype=bool)
    for _ in range(3 * n):
        gradient = atb - ata @ x
        candidates = (~passive) & (gradient > 1e-12)
        if not candidates.any():
            break
        passive[int(np.argmax(np.where(candidates, gradient, -np.inf)))] = True
        for _ in range(3 * n):
            index = np.flatnonzero(passive)
            try:
                solution = np.linalg.solve(
                    ata[np.ix_(index, index)], atb[index])
            except np.linalg.LinAlgError:
                solution, _, _, _ = np.linalg.lstsq(
                    ata[np.ix_(index, index)], atb[index], rcond=None)
            if (solution > 0.0).all():
                x = np.zeros(n)
                x[index] = solution
                break
            # Back off along the segment to the first zero crossing.
            ratios = np.where(solution <= 0.0,
                              x[index] / np.maximum(x[index] - solution, 1e-300),
                              np.inf)
            alpha = float(np.min(ratios))
            x[index] = x[index] + alpha * (solution - x[index])
            passive[index[x[index] <= 1e-14]] = False
            x[~passive] = 0.0
        else:
            break
    return x


# ---- RBF evaluation (mirrors interpolate_residual) -------------------------

def evaluate_nodes(nodes, yield_bw, shape_bw, query):
    yield_weight_sum = 0.0
    shape_weight_sum = 0.0
    log_yield_sum = 0.0
    shape_sum = np.zeros(BASIS_SIZE)
    tiny = sys.float_info.min
    for node_feature, log_yield, shape in nodes:
        dy = (query - node_feature) / yield_bw
        ds = (query - node_feature) / shape_bw
        wy = math.exp(-0.5 * float(dy @ dy))
        ws = math.exp(-0.5 * float(ds @ ds))
        yield_weight_sum += wy
        shape_weight_sum += ws
        log_yield_sum += wy * log_yield
        shape_sum += ws * shape
    conditional = math.exp(log_yield_sum / (yield_weight_sum + tiny))
    shape = shape_sum / (shape_weight_sum + tiny)
    total = shape.sum()
    return conditional, shape / (total + tiny)


def shape_total_variation(a: np.ndarray, b: np.ndarray) -> float:
    return 0.5 * float(np.abs(a - b).sum())


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    mu_table = load_photonuclear_mu(args.photonuclear_xs)
    cells = [c for c in (load_cell(p) for p in sorted(args.bank.glob("*_geant4.csv")))
             if c is not None]
    offgrid = []
    if args.offgrid is not None:
        if not args.offgrid.is_dir():
            raise SystemExit(f"off-grid bank not found: {args.offgrid}")
        offgrid = [c for c in (load_cell(p)
                               for p in sorted(args.offgrid.glob("*_geant4.csv")))
                   if c is not None]
    print(f"{len(cells)} bank cells, {len(offgrid)} off-grid cells")

    def baseline(cell: Cell) -> float:
        # Production (Geant4 attenuation table) times the mechanistic escape
        # of the secondaries through the remaining shield. Without the escape
        # term the fitted field would have to carry a ~1700x tau falloff, and
        # a kernel-mean interpolant over-predicts by 60-75% between the tau
        # training nodes.
        tau_pn = cell.areal * mu_photonuclear(mu_table, cell.energy, cell.z)
        return max((1.0 - math.exp(-tau_pn))
                   * escape_factor(cell.tau, args.escape_kappa), 1e-30)

    train = [c for c in cells if c.z not in args.holdout_z]
    heldout = [c for c in cells if c.z in args.holdout_z]
    print(f"train {len(train)}, held-out Z {args.holdout_z} {len(heldout)}, "
          f"off-grid {len(offgrid)}")

    # Node spacing in feature units sets the scale: energy ~0.017 per MeV
    # step, Z/92 as little as 0.011 (73 vs 74), and log1p(tau)/log(25)
    # ~0.16-0.28 between the tau nodes. The conditional yield falls by more
    # than an order of magnitude from tau=4 to tau=8, so the tau kernel must
    # stay well inside that spacing or it averages across the falloff.
    candidates_yield = [args.yield_bandwidth] if args.yield_bandwidth else [
        [e, z, t]
        for e in (0.006, 0.010, 0.016)
        for z in (0.020, 0.035, 0.055)
        for t in (0.035, 0.055, 0.080, 0.120)]
    candidates_shape = [args.shape_bandwidth] if args.shape_bandwidth else [
        [e, z, t]
        for e in (0.008, 0.014, 0.022)
        for z in (0.030, 0.050, 0.080)
        for t in (0.055, 0.090, 0.140)]

    best = None
    for warp_scale in args.warp_scale:
        shapes = {c.identity: fit_cell_shape(c, warp_scale)
                  for c in cells + offgrid}
        recon = [reconstruction_tv(c, shapes[c.identity], warp_scale)
                 for c in cells]
        recon_median = statistics.median(recon)
        recon_p95 = sorted(recon)[int(0.95 * len(recon)) - 1]
        print(f"  warp_scale={warp_scale}: basis reconstruction TV "
              f"med {100*recon_median:.2f}%  p95 {100*recon_p95:.2f}%")
        train_nodes = [
            (feature(c.energy, c.z, c.tau),
             math.log(max(c.yield_ / baseline(c), 1e-300)),
             shapes[c.identity])
            for c in train]
        for yield_bw in candidates_yield:
            for shape_bw in candidates_shape:
                ybw = np.asarray(yield_bw)
                sbw = np.asarray(shape_bw)
                metrics = {}
                for label, group in (("ongrid", train),
                                     ("heldout", heldout),
                                     ("offgrid", offgrid)):
                    yerr, tvs = [], []
                    for cell in group:
                        query = feature(cell.energy, cell.z, cell.tau)
                        conditional, shape = evaluate_nodes(
                            train_nodes, ybw, sbw, query)
                        predicted = conditional * baseline(cell)
                        yerr.append(abs(predicted / cell.yield_ - 1.0))
                        tvs.append(shape_total_variation(
                            shape, shapes[cell.identity]))
                    if yerr:
                        metrics[label] = (
                            statistics.median(yerr),
                            sorted(yerr)[int(0.95 * len(yerr)) - 1],
                            statistics.median(tvs),
                            sorted(tvs)[int(0.95 * len(tvs)) - 1])
                if "heldout" not in metrics:
                    continue
                # Select on the worst of the two validation families so a
                # configuration cannot win by overfitting the on-grid holdout.
                # Score on the worst validation family AND on training-cell
                # reproduction. A kernel-mean field is a weighted average, not
                # an interpolator: too wide a kernel biases even the cells it
                # was built from, and a holdout-only score is blind to that.
                score = (recon_median
                         + metrics["ongrid"][0] + metrics["ongrid"][2]
                         + max(metrics["heldout"][0] + metrics["heldout"][2],
                               (metrics.get("offgrid", metrics["heldout"])[0]
                                + metrics.get("offgrid",
                                              metrics["heldout"])[2])))
                if best is None or score < best["score"]:
                    best = {"score": score, "warp_scale": warp_scale,
                            "yield_bandwidth": list(yield_bw),
                            "shape_bandwidth": list(shape_bw),
                            "metrics": metrics, "shapes": shapes,
                            "reconstruction": (recon_median, recon_p95)}

    if best is None:
        raise SystemExit(
            "no configuration could be scored: the bank has no held-out-Z "
            f"cells (--holdout-z {args.holdout_z})")
    m = best["metrics"]
    print(f"\nselected warp_scale={best['warp_scale']} "
          f"yield_bw={best['yield_bandwidth']} "
          f"shape_bw={best['shape_bandwidth']}")
    print(f"  basis reconstruction TV med/p95 "
          f"{100*best['reconstruction'][0]:.2f}%/"
          f"{100*best['reconstruction'][1]:.2f}%")
    for label, values in m.items():
        print(f"  {label:8s} yield med/p95 {100*values[0]:.2f}%/"
              f"{100*values[1]:.2f}%   shape TV med/p95 "
              f"{100*values[2]:.2f}%/{100*values[3]:.2f}%")

    # Final payload: every cell becomes a node (the bank is small, and the
    # holdout scan above already established the interpolation quality).
    shapes = best["shapes"]
    # repr() of a numpy scalar is "np.float64(x)" under numpy 2, which the
    # C++ reader cannot parse: force plain Python floats everywhere.
    def num(value) -> str:
        return repr(float(value))

    lines = ["CEELO_PHOTONUCLEAR_SURROGATE 1",
             f"inverse_energy_warp_scale {num(best['warp_scale'])}",
             "basis cubic_bspline",
             f"escape_kappa {num(args.escape_kappa)}"]
    header = "photonuclear 31"
    header += "".join(f" {num(v)}" for v in best["yield_bandwidth"])
    header += "".join(f" {num(v)}" for v in best["shape_bandwidth"])
    header += f" {len(cells)}"
    lines.append(header)
    for cell in cells:
        f = feature(cell.energy, cell.z, cell.tau)
        log_yield = math.log(max(cell.yield_ / baseline(cell), 1e-300))
        row = "node photonuclear"
        row += "".join(f" {num(v)}" for v in f)
        row += f" {num(log_yield)}"
        row += "".join(f" {num(v)}" for v in shapes[cell.identity])
        lines.append(row)
    lines.append("END")
    surrogate = args.output / "photonuclear_surrogate.txt"
    surrogate.write_text("\n".join(lines) + "\n")

    # System impact of the photonuclear error on the escaping spectrum:
    # |model - truth| of the photonuclear yield in percentage points of each
    # held-out / off-grid cell's total escaping yield, evaluated with the
    # training nodes only (the cells stay held out).
    selected_nodes = [
        (feature(c.energy, c.z, c.tau),
         math.log(max(c.yield_ / baseline(c), 1e-300)),
         shapes[c.identity])
        for c in train]
    impact = {}
    for label, group in (("heldout", heldout), ("offgrid", offgrid)):
        rows = []
        for cell in group:
            conditional, _ = evaluate_nodes(
                selected_nodes, np.asarray(best["yield_bandwidth"]),
                np.asarray(best["shape_bandwidth"]),
                feature(cell.energy, cell.z, cell.tau))
            predicted = conditional * baseline(cell)
            rows.append({
                "identity": cell.identity,
                "photonuclear_truth": cell.yield_,
                "photonuclear_sigma": cell.sigma,
                "photonuclear_model": predicted,
                "escaping_total": cell.escaping,
                "added_percentage_points": 100.0 * abs(
                    predicted - cell.yield_) / max(cell.escaping, 1e-300),
            })
        if rows:
            impact[label] = {
                "max_added_percentage_points": max(
                    r["added_percentage_points"] for r in rows),
                "cells": rows,
            }
    for label, values in impact.items():
        print(f"  {label:8s} worst photonuclear error "
              f"{values['max_added_percentage_points']:.3f} percentage "
              "points of the escaping spectrum")

    report = {
        "cells": len(cells), "train": len(train),
        "holdout_z": args.holdout_z,
        "offgrid_cells": len(offgrid),
        "warp_scale": best["warp_scale"],
        "escape_kappa": args.escape_kappa,
        "yield_bandwidth": best["yield_bandwidth"],
        "shape_bandwidth": best["shape_bandwidth"],
        "reconstruction_tv_median": best["reconstruction"][0],
        "reconstruction_tv_p95": best["reconstruction"][1],
        "metrics": {k: {"yield_median": v[0], "yield_p95": v[1],
                        "shape_tv_median": v[2], "shape_tv_p95": v[3]}
                    for k, v in m.items()},
        "system_impact": impact,
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=1))
    print(f"\nsurrogate -> {surrogate}\nreport -> {args.output/'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
