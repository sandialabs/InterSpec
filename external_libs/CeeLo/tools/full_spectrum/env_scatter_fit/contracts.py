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

"""Canonical, dependency-free contracts shared by campaign and archive code."""

from __future__ import annotations

import json
from typing import Any


STOPPING_VECTOR_SCHEMA = "CEEFIT-STOPPING-v1"
STOPPING_VECTOR_DIMENSION = 16
STOPPING_ENERGY_BANDS = 12
FORMAL_BLOCK_HISTORIES = 262_144
# Blocks 1..RETAINED_RAW_BLOCKS of every node stay in the archive as raw
# records (the uncertainty calibration pairs them); every later block is folded
# into the node's running sums as it is ingested and keeps only a ledger entry.
RETAINED_RAW_BLOCKS = 4


def stopping_vector_definition(
    energy_axis_index: int,
    energy_axis_name: str,
    energy_bins: int,
) -> dict[str, Any]:
    """Return the exact semantic definition of the 16 stopping observables."""

    if energy_axis_index < 0 or not energy_axis_name or energy_bins <= 0:
        raise ValueError("invalid total-energy axis for stopping-vector definition")
    entries: list[dict[str, Any]] = []
    # FluenceTally assigns cell i to floor(i * 12 / n_energy).  These ceil
    # boundaries describe that exact index partition, including empty bands in
    # small non-production test axes.  Production's 384 bins give 32/bin band.
    for band in range(STOPPING_ENERGY_BANDS):
        begin = (band * energy_bins + STOPPING_ENERGY_BANDS - 1) // STOPPING_ENERGY_BANDS
        end = (
            (band + 1) * energy_bins + STOPPING_ENERGY_BANDS - 1
        ) // STOPPING_ENERGY_BANDS
        entries.append(
            {
                "energy_axis_index": energy_axis_index,
                "energy_axis_name": energy_axis_name,
                "energy_bin_begin": begin,
                "energy_bin_end": end,
                "kind": "registered_energy_index_band",
                "name": f"registered_band_{band:02d}",
            }
        )
    entries.extend(
        {
            "kind": "legendre_angular_moment",
            "legendre_order": order,
            "name": f"legendre_l{order}",
        }
        for order in range(1, 5)
    )
    return {
        "angular_basis": "standard Legendre P_l(mu), mu=cos(direction_of_travel,+detector_axis)",
        "dimension": STOPPING_VECTOR_DIMENSION,
        "entries": entries,
        "normalization": "raw per-history scattered continuum-plus-lines fluence moments",
        "version": STOPPING_VECTOR_SCHEMA,
    }


def stopping_vector_definition_json(
    energy_axis_index: int,
    energy_axis_name: str,
    energy_bins: int,
) -> str:
    return json.dumps(
        stopping_vector_definition(energy_axis_index, energy_axis_name, energy_bins),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
