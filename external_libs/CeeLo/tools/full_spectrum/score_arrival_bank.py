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

"""Score a detector payload on an arrival bank (photons from every side).

Each cell of the bank is Monte Carlo truth for an isotropic field over a band
of direction cosines (`ceelo_detector_response_campaign --incidence`); the
payload is averaged over the same band by `ceelo_detector_model_evaluate
--incidence` (Gauss-Legendre in the polar angle per hemisphere, projected-area
weighted).  This is the response a scattered field is folded with, which the
fixed-angle banks never test as a whole.  Per cell, relative errors with z
against the truth's own MC sigma (the model side is deterministic):

  detected yield       1 - zero_deposit, per envelope photon; also its
                       excess over 2 sigma of the truth's noise (cells far
                       below the passive landmark detect a few percent of
                       their photons, so 4M events leave ~0.2% noise there)
  full-energy atom     per envelope photon
  continuum TV         total variation of the deposited spectrum without the
                       zero-deposit and full-energy atoms (atoms in their
                       bins), normalized; with its MC noise floor (the mean TV
                       of multinomial redraws of the truth against itself)
  valley density       continuum in [E_C + (E - E_C)/4, E - max(3 keV, 2% E)]
                       below the full-energy line, per envelope photon
  envelope area        the boundary's mean projected area over the band, MC
                       (bounding disk x acceptance) against the payload's
                       envelope; a mismatch would bias every fold

The report groups cells by arrival class and energy band (40-100, 100-400,
400-1000 keV, and >= 1000 keV if present); check_detector_gates.py --arrival
applies the bounds.

    score_arrival_bank.py BANK --payload P [--evaluator PATH] --output DIR
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import bank_io

REPO = Path(__file__).resolve().parents[2]
ELECTRON_MASS_KEV = 510.99895
ENERGY_BANDS_KEV = ((40.0, 100.0), (100.0, 400.0), (400.0, 1000.0),
                    (1000.0, math.inf))
NOISE_DRAWS = 40
# The campaign's opt-in refined pair labels re-split rows already counted.
REFINED_PAIR_LABELS = frozenset((
    "pair_escape0", "pair_escape1_unscattered", "pair_escape1_scattered",
    "pair_escape2_unscattered", "pair_escape2_scattered", "pair_in_flight"))


def read_cell(stream) -> dict:
    """One arrival cell: energy, class, events, atoms {name: (sum_w, sum_w2,
    energy)} and the continuum summed over components {bin lower edge:
    (sum_w, sum_w2)} on the campaign's uniform grid over [0, E)."""
    cell: dict = {"atoms": {}, "bins": {}}
    for row in csv.DictReader(stream):
        if row["cell_kind"] != "detector_arrival":
            raise ValueError(f"cell_kind {row['cell_kind']!r} is not an "
                             "arrival cell")
        if row["component"] in REFINED_PAIR_LABELS:
            continue
        energy = float(row["energy_keV"])
        cell.setdefault("energy_keV", energy)
        cell.setdefault("incidence", row["incidence"])
        cell.setdefault("events", float(row["events"]))
        cell.setdefault("area", float(row["mean_projected_area_cm2"]))
        cell.setdefault("area_sigma",
                        float(row["mean_projected_area_sigma_cm2"]))
        lo, hi = float(row["bin_lo_keV"]), float(row["bin_hi_keV"])
        value = (float(row["sum_w"]), float(row["sum_w2"]))
        if lo == hi:
            old = cell["atoms"].get(row["component"], (0.0, 0.0, lo))
            cell["atoms"][row["component"]] = (old[0] + value[0],
                                               old[1] + value[1], lo)
        else:
            old = cell["bins"].get(lo, (0.0, 0.0))
            cell["bins"][lo] = (old[0] + value[0], old[1] + value[1])
    return cell


def bin_index(energy_keV: float, bins: int, lower_keV: float) -> int:
    # Truth and model print the edges E b / bins at different precisions, so
    # bins are matched by index, not by their printed edges.
    return min(bins - 1, max(0, int(round(lower_keV / energy_keV * bins))))


def atom(cell: dict, name: str) -> float:
    return cell["atoms"].get(name, (0.0, 0.0))[0]


def spectrum(cell: dict, energy_keV: float, bins: int,
             with_fep: bool) -> np.ndarray:
    """The deposited spectrum on `bins` bins over [0, E): continuum plus atoms
    in their bins, the zero-deposit atom (and, without FEP, the full-energy
    atom) left out."""
    values = np.zeros(bins)
    for lo, (w, _) in cell["bins"].items():
        values[bin_index(energy_keV, bins, lo)] += w
    for name, (w, _, energy) in cell["atoms"].items():
        if name == "zero_deposit" or (name == "full_energy" and not with_fep):
            continue
        values[min(bins - 1, int(energy / energy_keV * bins))] += w
    return values


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    if a.sum() <= 0.0 or b.sum() <= 0.0:
        return float("nan")
    return 0.5 * float(np.abs(a / a.sum() - b / b.sum()).sum())


def noise_floor(counts: np.ndarray, seed: int) -> float:
    """Mean TV of multinomial redraws of an unweighted truth against itself."""
    total = int(round(counts.sum()))
    if total <= 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    p = counts / counts.sum()
    return float(np.mean([
        total_variation(rng.multinomial(total, p).astype(float), counts)
        for _ in range(NOISE_DRAWS)]))


def valley_window(energy_keV: float) -> tuple[float, float]:
    ratio = 2.0 * energy_keV / ELECTRON_MASS_KEV
    compton_edge = energy_keV * ratio / (1.0 + ratio)
    return (compton_edge + 0.25 * (energy_keV - compton_edge),
            energy_keV - max(3.0, 0.02 * energy_keV))


def valley(cell: dict, energy_keV: float, bins: int) -> tuple[float, float]:
    """Continuum in the valley window: (sum_w, sum_w2), whole bins inside."""
    lo, hi = valley_window(energy_keV)
    first = math.ceil(lo / energy_keV * bins - 1e-9)
    last = math.floor(hi / energy_keV * bins + 1e-9) - 1
    w = w2 = 0.0
    for b_lo, (sw, sw2) in cell["bins"].items():
        if first <= bin_index(energy_keV, bins, b_lo) <= last:
            w += sw
            w2 += sw2
    return w, w2


def score_cell(path: Path, payload: Path, evaluator: Path,
               nodes: int) -> dict:
    with bank_io.open_cell(path) as stream:
        truth = read_cell(stream)
    energy, incidence = truth["energy_keV"], truth["incidence"]
    bins = len(truth["bins"])
    with tempfile.TemporaryDirectory() as scratch:
        output = Path(scratch) / "model.csv"
        subprocess.run(
            [str(evaluator), "--payload", str(payload),
             "--energy-kev", f"{energy:.17g}", "--incidence", incidence,
             "--nodes-per-hemisphere", str(nodes),
             "--bins", str(bins), "--output", str(output)],
            check=True, stdout=subprocess.DEVNULL)
        with output.open(newline="") as stream:
            model = read_cell(stream)

    n = truth["events"]
    seed = int.from_bytes(
        hashlib.sha256(bank_io.cell_id(path).encode()).digest()[:8], "big")

    def relative(model_value: float, truth_sum: float, truth_sum2: float):
        if truth_sum <= 0.0:
            return float("nan"), float("nan")
        truth_value = truth_sum / n
        sigma = math.sqrt(max(truth_sum2, 0.0)) / n
        return (model_value / truth_value - 1.0,
                (model_value - truth_value) / sigma if sigma > 0 else float("nan"))

    zero = atom(truth, "zero_deposit")
    detected_truth = 1.0 - zero / n
    detected_sigma = math.sqrt(max(detected_truth * (1 - detected_truth), 0) / n)
    detected_model = 1.0 - atom(model, "zero_deposit")
    detected_error = detected_model / detected_truth - 1.0
    detected_sigma_rel = detected_sigma / detected_truth
    fep_w, fep_w2, _ = truth["atoms"].get("full_energy", (0.0, 0.0, 0.0))
    fep_error, fep_z = relative(atom(model, "full_energy"), fep_w, fep_w2)
    valley_w, valley_w2 = valley(truth, energy, bins)
    valley_error, valley_z = relative(valley(model, energy, bins)[0],
                                      valley_w, valley_w2)
    truth_full = spectrum(truth, energy, bins, True)
    truth_continuum = spectrum(truth, energy, bins, False)
    tv_continuum = total_variation(spectrum(model, energy, bins, False),
                                   truth_continuum)
    tv_continuum_noise = noise_floor(truth_continuum, seed)
    return {
        "cell": bank_io.cell_id(path),
        "energy_keV": energy,
        "incidence": incidence,
        "events": n,
        "detected_truth": detected_truth,
        "detected_yield_relative_error": detected_error,
        "detected_yield_sigma": detected_sigma_rel,
        "detected_yield_excess": max(0.0, abs(detected_error)
                                     - 2.0 * detected_sigma_rel),
        "detected_yield_z": (detected_model - detected_truth) / detected_sigma
        if detected_sigma > 0 else float("nan"),
        "full_energy_per_detected_truth": fep_w / n / detected_truth,
        "full_energy_relative_error": fep_error,
        "full_energy_z": fep_z,
        "pulse_shape_total_variation": total_variation(
            spectrum(model, energy, bins, True), truth_full),
        "pulse_shape_total_variation_noise": noise_floor(truth_full, seed + 1),
        "continuum_total_variation": tv_continuum,
        "continuum_total_variation_noise": tv_continuum_noise,
        "continuum_tv_over_noise": tv_continuum / tv_continuum_noise
        if tv_continuum_noise > 0 else float("nan"),
        "valley_truth_per_photon": valley_w / n,
        "valley_density_relative_error": valley_error,
        "valley_density_z": valley_z,
        "envelope_area_relative_error": model["area"] / truth["area"] - 1.0,
        "envelope_area_z": (model["area"] - truth["area"]) / truth["area_sigma"]
        if truth["area_sigma"] > 0 else float("nan"),
    }


GROUP_METRICS = (
    ("detected_yield_relative_error", True),
    ("detected_yield_excess", False),
    ("full_energy_relative_error", True),
    ("continuum_tv_over_noise", False),
    ("pulse_shape_total_variation", False),
    ("valley_density_relative_error", True),
    ("envelope_area_z", True),
)


def statistics(values: list[float]) -> dict:
    finite = np.array([v for v in values if math.isfinite(v)])
    if finite.size == 0:
        return {"count": 0}
    return {"count": int(finite.size), "median": float(np.median(finite)),
            "p95": float(np.percentile(finite, 95)),
            "max": float(finite.max())}


def groups(records: list[dict]) -> list[dict]:
    out = []
    classes = sorted({r["incidence"] for r in records})
    for incidence in classes:
        for lo, hi in ENERGY_BANDS_KEV:
            members = [r for r in records if r["incidence"] == incidence
                       and lo <= r["energy_keV"] < hi]
            if not members:
                continue
            group = {"incidence": incidence,
                     "energy_band_keV": [lo, hi if math.isfinite(hi) else None],
                     "cells": len(members)}
            for metric, absolute in GROUP_METRICS:
                values = [abs(r[metric]) if absolute else r[metric]
                          for r in members]
                group[metric] = statistics(values)
            out.append(group)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bank", type=Path,
                        help="arrival bank directory (cells/*_ceelo.csv.gz)")
    parser.add_argument("--payload", type=Path, required=True)
    parser.add_argument(
        "--evaluator", type=Path,
        default=REPO / "build-rel/tools/full_spectrum/"
                       "ceelo_detector_model_evaluate")
    parser.add_argument("--nodes-per-hemisphere", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    paths = bank_io.cell_files(args.bank / "cells")
    if not paths:
        print(f"no arrival cells under {args.bank}", file=sys.stderr)
        return 2
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(
            lambda path: score_cell(path, args.payload, args.evaluator,
                                    args.nodes_per_hemisphere), paths))
    records.sort(key=lambda r: (r["incidence"], r["energy_keV"]))
    report = {
        "payload": str(args.payload),
        "nodes_per_hemisphere": args.nodes_per_hemisphere,
        "groups": groups(records),
        "records": records,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "report.json").write_text(json.dumps(report, indent=1) + "\n")

    print(f"{'class':6s} {'E keV':>7s} {'events':>8s} {'det err':>8s} {'z':>6s}"
          f" {'FEP err':>8s} {'z':>6s} {'TVc/noise':>9s} {'valley':>7s} {'z':>6s}"
          f" {'area z':>6s}")
    for r in records:
        print(f"{r['incidence']:6s} {r['energy_keV']:7.2f} {r['events']:8.0f}"
              f" {r['detected_yield_relative_error']:8.3%} {r['detected_yield_z']:6.1f}"
              f" {r['full_energy_relative_error']:8.3%} {r['full_energy_z']:6.1f}"
              f" {r['continuum_tv_over_noise']:9.2f}"
              f" {r['valley_density_relative_error']:7.2%} {r['valley_density_z']:6.1f}"
              f" {r['envelope_area_z']:6.1f}")
    print(f"\nreport: {args.output / 'report.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
