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

"""Run the CeeLo truth banks for the source-secondary fit.

Three presets, each a separate bank directory of `cells/<id>_ceelo.csv.gz` plus
`campaign_manifest.json`:

* the training grid (default): 11 energies x 17 Z x 6 tau = 1122 cells, the
  Z = 13, 35, 64 slices marked "test" (whole-element holdouts);
* `--offgrid`: 75 cells at optical depths BETWEEN the training nodes, never
  trained on (the only between-node holdout);
* `--hydrogen`: 672 two-component shields (hydrogen mass fraction f_H plus
  element Z), the training data for the f_H coordinate.

Every cell runs `ceelo_source_escape_campaign` with Rayleigh and bound-electron
Compton Doppler on (the campaign defaults, passed explicitly) and
fluorescence off.  Seeds derive from each cell's coordinates, so a bank is
reproducible cell by cell; existing cells are skipped, so an interrupted run
resumes where it stopped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import bank_io
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# The 59/100/250 keV rows extend the bank down to the acceptance grid's 59 keV
# floor; without them the surrogate would extrapolate ~15 yield-bandwidths in
# energy below its support.
TRAINING_ENERGIES_KEV = [
    59.0, 100.0, 250.0,
    300.0, 661.657, 1000.0, 1045.0, 1500.0,
    3000.0, 6000.0, 12000.0,
]
TRAIN_Z = [6, 10, 16, 20, 26, 32, 38, 42, 50, 60, 70, 74, 82, 92]
HOLDOUT_Z = [13, 35, 64]
# tau 0.5 and 2.0 close the mid-tau interpolation defect between the
# {0.2, 1, 4, 8} nodes.
TRAINING_TAU = [0.2, 0.5, 1.0, 2.0, 4.0, 8.0]

# Held out of training entirely: every value sits between training tau nodes, so
# on-grid metrics cannot see between-node behaviour.  Scored as its own column.
OFFGRID_TAU = [0.35, 0.7, 1.4, 3.0, 6.0]
OFFGRID_Z = [13, 26, 35, 64, 82]
OFFGRID_ENERGIES_KEV = [100.0, 661.657, 3000.0]

# Hydrogenous shields.  0.1119 is water's hydrogen mass fraction and 0.1437
# polyethylene's; (Z 8, water) and (Z 6, polyethylene) are the exact-material
# anchors, and Z 35 at the water fraction is a whole-element holdout.
HYDROGEN_ENERGIES_KEV = [100.0, 200.0, 300.0, 662.0, 1045.0, 1500.0,
                         3000.0, 6000.0]
HYDROGEN_Z = [6, 13, 26, 50, 82]
HYDROGEN_TAU = [0.2, 1.0, 4.0, 8.0]
WATER_FH = 0.1119
POLYETHYLENE_FH = 0.1437
HYDROGEN_FH = [0.05, WATER_FH, POLYETHYLENE_FH, 0.25]


def hydrogen_cells() -> list[tuple[float, int, float, float]]:
    """The 672 unique hydrogen-bank cells (the polyethylene anchors repeat
    grid cells and are kept once)."""
    cells: list[tuple[float, int, float, float]] = []
    for energy in HYDROGEN_ENERGIES_KEV:
        for z_value in HYDROGEN_Z:
            for tau in HYDROGEN_TAU:
                for f_h in HYDROGEN_FH:
                    cells.append((energy, z_value, tau, f_h))
        for z_value, f_h in ((8, WATER_FH), (6, POLYETHYLENE_FH),
                             (35, WATER_FH)):
            for tau in (1.0, 4.0):
                cells.append((energy, z_value, tau, f_h))
    return list(dict.fromkeys(cells))


def identity(energy: float, z_value: int, tau: float, f_h: float = 0.0) -> str:
    e_text = f"{energy:g}".replace(".", "p")
    t_text = f"{tau:g}".replace(".", "p")
    text = f"E{e_text}_Z{z_value}_tau{t_text}"
    if f_h > 0.0:
        text += "_fH" + f"{f_h:g}".replace(".", "p")
    return text


def seed_for(energy: float, z_value: int, tau: float, f_h: float = 0.0) -> int:
    text = f"source-secondary-v1|{energy:.12g}|{z_value}|{tau:.12g}"
    if f_h > 0.0:
        text = f"source-hydrogen|{energy:.12g}|{z_value}|{tau:.12g}|{f_h:.12g}"
    return int.from_bytes(
        hashlib.sha256(text.encode()).digest()[:8], "big") & ((1 << 63) - 1)


def run_cell(
    driver: Path,
    output: Path,
    cell: tuple[float, int, float, float],
    events: int,
    bins: int,
    force: bool,
) -> tuple[str, str]:
    energy, z_value, tau, f_h = cell
    cell_id = identity(*cell)
    path = bank_io.written_path(output / "cells", cell_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    if bank_io.existing(path) is not None and not force:
        return cell_id, "cached"
    command = [
        str(driver),
        "--energy-kev", str(energy),
        "--z", str(z_value),
        "--tau", str(tau),
        "--events", str(events),
        "--bins", str(bins),
        "--seed", str(seed_for(*cell)),
        "--source-rayleigh",
        "--source-compton-doppler",
    ]
    if f_h > 0.0:
        command += ["--hydrogen-mass-fraction", str(f_h)]
    completed = subprocess.run(
        command, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bank_io.write_cell(path, completed.stdout)
    return cell_id, "ran"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--driver", type=Path,
        default=Path("build-rel/tools/full_spectrum/"
                     "ceelo_source_escape_campaign"))
    parser.add_argument("--events", type=int, default=100000)
    parser.add_argument("--bins", type=int, default=256)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--cell", action="append", default=[],
        metavar="ENERGY,Z,TAU[,F_H]",
        help="run explicit cells instead of the preset grid; repeatable")
    preset = parser.add_mutually_exclusive_group()
    preset.add_argument(
        "--offgrid", action="store_true",
        help="the off-grid tau validation set instead of the training grid")
    preset.add_argument(
        "--hydrogen", action="store_true",
        help="the hydrogen-fraction bank instead of the training grid")
    args = parser.parse_args()

    if args.cell:
        cells = []
        for text in args.cell:
            fields = text.split(",")
            if len(fields) not in (3, 4):
                raise ValueError(
                    f"--cell must be ENERGY,Z,TAU[,F_H], received: {text}")
            cells.append((
                float(fields[0]), int(fields[1]), float(fields[2]),
                float(fields[3]) if len(fields) == 4 else 0.0))
    elif args.offgrid:
        cells = [
            (energy, z_value, tau, 0.0)
            for energy in OFFGRID_ENERGIES_KEV
            for z_value in OFFGRID_Z
            for tau in OFFGRID_TAU
        ]
    elif args.hydrogen:
        cells = hydrogen_cells()
    else:
        cells = [
            (energy, z_value, tau, 0.0)
            for energy in TRAINING_ENERGIES_KEV
            for z_value in TRAIN_Z + HOLDOUT_Z
            for tau in TRAINING_TAU
        ]
    cells = list(dict.fromkeys(cells))
    args.output.mkdir(parents=True, exist_ok=True)
    campaign = ("source-secondary-offgrid-tau" if args.offgrid
                else "source-secondary-hydrogen" if args.hydrogen
                else "source-secondary-training")
    manifest = {
        "schema": 1,
        "campaign": campaign,
        "energies_keV": sorted({cell[0] for cell in cells}),
        "training_Z": sorted({
            cell[1] for cell in cells if cell[1] not in HOLDOUT_Z}),
        "holdout_Z": sorted({
            cell[1] for cell in cells if cell[1] in HOLDOUT_Z}),
        "tau": sorted({cell[2] for cell in cells}),
        "hydrogen_mass_fractions": sorted({cell[3] for cell in cells}),
        "events_per_cell": args.events,
        "bins": args.bins,
        "physics": {
            "rayleigh": True,
            "compton_incoherent_scattering_function": True,
            "compton_doppler": True,
            "pair_annihilation": True,
            "positron_transport": True,
            "positron_in_flight_annihilation": True,
            "annihilation_components": [
                "511_atom", "scattered_at_rest", "in_flight"],
            "source_electron_bremsstrahlung": True,
            "fluorescence": True,
            "impact_ionization": True,
        },
        "cells": [
            {
                "id": identity(*cell),
                "energy_keV": cell[0],
                "Z": cell[1],
                "tau": cell[2],
                "hydrogen_mass_fraction": cell[3],
                "seed": seed_for(*cell),
                "split": (
                    "offgrid" if args.offgrid
                    else "test" if cell[1] in HOLDOUT_Z
                    else "candidate"),
            }
            for cell in cells
        ],
    }
    (args.output / "campaign_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(f"{campaign}: {len(cells)} cells at {args.events} events")

    counts = {"ran": 0, "cached": 0}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                run_cell, args.driver, args.output, cell,
                args.events, args.bins, args.force)
            for cell in cells
        ]
        try:
            for index, future in enumerate(as_completed(futures), 1):
                cell_id, status = future.result()
                counts[status] += 1
                if index % 20 == 0 or index == len(futures):
                    print(
                        f"{index}/{len(futures)} cells; "
                        f"ran={counts['ran']} cached={counts['cached']}; "
                        f"latest={cell_id}")
        except BaseException:
            # Stop queued cells at once; finished cells stay on disk and a
            # re-run resumes after them.
            executor.shutdown(wait=True, cancel_futures=True)
            raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
