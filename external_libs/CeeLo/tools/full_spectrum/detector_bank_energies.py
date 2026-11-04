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

"""Place a detector bake's bank energies for one detector.

Below 100 keV the response turns on through the passive layers, so energies
there are written for the reference detector (hpge_coax, whose passive-layer
landmark is 55 keV) and scaled by landmark/55 for the detector being baked.
Held-out banks must stay off the training grid: a held-out energy that lands
within 2% of a training energy is moved, by the smallest multiple of 0.5%,
until it clears every one.  Every energy is rounded to 0.05 keV, so a detector
whose landmark is within ~0.05% of 55 keV keeps the reference energies exactly
(and its finished cells).

    detector_bank_energies.py --landmark 63.44 \\
        --training training=10,20,40,... --training-scaled lowE=35,37.5,... \\
        --held-out-scaled audit_a=42.5,52.5,...

prints one `NAME=E1,E2,...` line per bank for a shell to read; moved energies
are reported on stderr.
"""

from __future__ import annotations

import argparse
import sys

REFERENCE_LANDMARK_KEV = 55.0
SCALED_BELOW_KEV = 100.0
ROUNDING_KEV = 0.05
MIN_SEPARATION = 0.02
SHIFT_STEP = 0.005
MAX_SHIFT = 0.2


def rounded(energy_keV: float) -> float:
    # An integer over 20 is the double nearest the 0.05 keV value, so equal
    # placements give equal identities.
    return round(energy_keV / ROUNDING_KEV) / round(1.0 / ROUNDING_KEV)


def landmark_scaled(energies_keV: list[float], landmark_keV: float) -> list[float]:
    """Energies below 100 keV scaled by landmark/55; the rest unchanged."""
    if not landmark_keV > 0.0:
        raise ValueError("the landmark must be positive")
    factor = landmark_keV / REFERENCE_LANDMARK_KEV
    return [rounded(e * factor) if e < SCALED_BELOW_KEV else e
            for e in energies_keV]


def separation(a: float, b: float) -> float:
    return abs(a - b) / min(a, b)


def clear_of(energy_keV: float, training_keV: list[float]) -> bool:
    return all(separation(energy_keV, t) >= MIN_SEPARATION for t in training_keV)


def held_out(energies_keV: list[float], training_keV: list[float]
             ) -> tuple[list[float], list[tuple[float, float]]]:
    """Each energy moved off the training grid if it is on it, and the moves.

    Candidates are E (1 + k 0.5%) for k = 1, -1, 2, -2, ...: the nearest
    position at least 2% from every training energy wins (up to 20%).
    """
    placed: list[float] = []
    moves: list[tuple[float, float]] = []
    for energy in energies_keV:
        if clear_of(energy, training_keV):
            placed.append(energy)
            continue
        steps = int(round(MAX_SHIFT / SHIFT_STEP))
        for k in (s * sign for s in range(1, steps + 1) for sign in (1, -1)):
            candidate = rounded(energy * (1.0 + k * SHIFT_STEP))
            if clear_of(candidate, training_keV):
                placed.append(candidate)
                moves.append((energy, candidate))
                break
        else:
            raise ValueError(
                f"no held-out energy within {MAX_SHIFT:.0%} of {energy} keV "
                "is 2% clear of the training grid")
    if len(set(placed)) != len(placed):
        raise ValueError(f"held-out energies coincide after placement: {placed}")
    return placed, moves


def named_list(text: str) -> tuple[str, list[float]]:
    name, sep, values = text.partition("=")
    if not sep or not name:
        raise argparse.ArgumentTypeError(f"{text!r} is not NAME=E1,E2,...")
    return name, [float(v) for v in values.split(",") if v.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--landmark", type=float, required=True,
                        help="the detector's passive-layer landmark, keV")
    for flag, help_text in (
            ("--training", "training bank, energies used as given"),
            ("--training-scaled", "training bank, landmark-relative below 100 keV"),
            ("--held-out-scaled", "held-out bank, landmark-relative below "
                                  "100 keV and kept off the training grid")):
        parser.add_argument(flag, type=named_list, action="append", default=[],
                            metavar="NAME=E1,E2,...", help=help_text)
    args = parser.parse_args()

    banks: list[tuple[str, list[float]]] = []
    training: list[float] = []
    for name, energies in args.training:
        banks.append((name, energies))
        training += energies
    for name, energies in args.training_scaled:
        scaled = landmark_scaled(energies, args.landmark)
        banks.append((name, scaled))
        training += scaled
    for name, energies in args.held_out_scaled:
        placed, moves = held_out(landmark_scaled(energies, args.landmark), training)
        for before, after in moves:
            print(f"{name}: {before:g} keV is on the training grid; held out at "
                  f"{after:g} keV", file=sys.stderr)
        banks.append((name, placed))
    for name, energies in banks:
        print(f"{name}=" + ",".join(f"{e:.10g}" for e in energies))
    return 0


if __name__ == "__main__":
    sys.exit(main())
