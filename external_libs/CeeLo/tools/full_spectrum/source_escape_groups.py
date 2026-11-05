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

"""Physical component groups of a source-escape CSV, shared by the scorers.

CeeLo truth, the functional model and the Geant4 harness label escaping
photons with overlapping vocabularies (the model merges the Compton orders the
truth splits, the truth splits what the model merges), so every total is
summed over these physical GROUPS.  A naive sum over the raw rows would
double-count.  `absorbed` is not an escaping photon.
"""

from __future__ import annotations

import csv
import bank_io
import math
from collections import OrderedDict
from pathlib import Path

GROUPS = OrderedDict([
    ("direct", {"uncollided"}),
    ("elastic", {"elastic_line"}),
    ("single Compton", {
        "single_compton", "compton_rayleigh",
        "single_compton_rayleigh",
    }),
    ("multiple Compton", {
        "multiple_compton", "double_compton", "higher_compton",
        "double_compton_rayleigh", "higher_compton_rayleigh",
        "multiple_compton_rayleigh", "triple_compton",
        "quadruple_compton", "five_plus_compton",
        "triple_compton_rayleigh", "quadruple_compton_rayleigh",
        "five_plus_compton_rayleigh",
    }),
    ("511-keV atom", {"pair_511"}),
    ("at-rest annihilation scatter", {"pair_annihilation_continuum"}),
    ("in-flight annihilation", {"pair_annihilation_in_flight"}),
    ("photoelectron brems", {"photoelectron_brems"}),
    ("Compton-electron brems", {"compton_electron_brems"}),
    ("pair-electron brems", {"pair_electron_brems"}),
    ("fluorescence K-alpha", {"fluorescence_k_alpha"}),
    ("fluorescence K-beta", {"fluorescence_k_beta"}),
    ("fluorescence L", {"fluorescence_l"}),
    # Compton-scattered shield x-rays: tallied by the campaign, not modelled.
    ("fluorescence continuum", {"fluorescence_continuum"}),
    ("photonuclear", {"photonuclear_secondary"}),
])
EXCLUDED = {"absorbed"}
# Tallied by the campaign but not produced by the model: kept out of the
# default (leakage) total and reported on its own.
UNMODELLED = {"fluorescence_continuum"}
KNOWN = {label for labels in GROUPS.values() for label in labels}
# The uncollided term is analytic exp(-sum tau) and exact by construction, so
# scores exclude it: at thin shields it IS most of the spectrum and would
# flatter every number ("collided-only leakage").
DIRECT = set(GROUPS["direct"])
# The four channels the bremsstrahlung escape kernel carries.
BREMS = {"photoelectron_brems", "compton_electron_brems",
         "pair_electron_brems", "pair_annihilation_in_flight"}
# Every component of the source-secondary surrogate.
FLUORESCENCE = {"fluorescence_k_alpha", "fluorescence_k_beta", "fluorescence_l"}
SECONDARY = BREMS | FLUORESCENCE | {"pair_511", "pair_annihilation_continuum"}


def read_rows(source: Path | str) -> list[dict[str, str]]:
    """Rows of a campaign CSV given as a path or as the CSV text itself."""
    if isinstance(source, Path):
        with bank_io.open_cell(source) as stream:
            return list(csv.DictReader(stream))
    return list(csv.DictReader(source.splitlines()))


def totals(rows: list[dict[str, str]], labels: set[str] | None = None
           ) -> tuple[float, float, set[str]]:
    """(yield per primary, its MC sigma, ungrouped labels seen).

    `labels` selects the component labels summed; by default the
    collided-only leakage (every grouped label except the uncollided one and
    the UNMODELLED ones).
    The model's rows carry sum_w2 = 0, so its sigma is 0.
    """
    if not rows:
        return 0.0, 0.0, set()
    events = float(rows[0]["events"])
    total = variance = 0.0
    ungrouped: set[str] = set()
    for row in rows:
        name = row["component"]
        if name in EXCLUDED:
            continue
        if name not in KNOWN:
            ungrouped.add(name)
            continue
        if labels is None:
            if name in DIRECT or name in UNMODELLED:
                continue
        elif name not in labels:
            continue
        total += float(row["sum_w"]) / events
        variance += float(row.get("sum_w2") or 0.0) / (events * events)
    return total, math.sqrt(max(variance, 0.0)), ungrouped
