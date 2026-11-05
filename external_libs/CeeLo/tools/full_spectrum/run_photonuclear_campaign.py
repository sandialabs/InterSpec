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

"""Photonuclear (E, Z, tau) Geant4 bank runner.

Runs `g4_source_escape_campaign` over an energy x Z x tau grid at 4M
histories per cell with deterministic hashed seeds; each cell is written as
`<output>/<id>_geant4.csv`, and cells already on disk are skipped, so an
interrupted run resumes where it stopped.  Requires the Geant4 environment
(source geant4.sh) before running.

Presets:

* default: the training bank, 152 cells.  Train Z {13, 26, 29, 42, 50, 74,
  82, 92}, whole-cell holdouts {47, 73}, tau {0.2, 1, 4, 8}, and per Z the
  energies 8-12 MeV above an importance threshold (cells where photonuclear
  exceeds ~0.2% of escaping photons);
* `--offgrid`: 11 validation cells at optical depths BETWEEN the training
  nodes and an element (Z = 60) outside the trained set.

`--grid-json` or the explicit grid options override a preset.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]

DEFAULT_ENERGIES_KEV = [8000, 9000, 10000, 11000, 12000]
DEFAULT_ZS = [13, 26, 29, 42, 47, 50, 73, 74, 82, 92]
DEFAULT_TAUS = [0.2, 1.0, 4.0, 8.0]
# Lowest energy per Z at which photonuclear secondaries exceed ~0.2% of the
# escaping photons (the importance threshold that sets the bank's extent).
MIN_ENERGY_BY_Z_KEV = {
    13: 11000.0, 26: 11000.0,
    29: 9000.0, 42: 9000.0, 47: 9000.0, 50: 9000.0, 73: 9000.0, 74: 9000.0,
    82: 8000.0, 92: 8000.0,
}
# (energy_keV, Z, tau): taus between the training nodes, and Z = 60.
OFFGRID_CELLS = [
    (9000.0, 50, 0.5), (9000.0, 50, 2.0), (9000.0, 82, 0.5),
    (9000.0, 82, 2.0), (12000.0, 50, 0.5), (12000.0, 50, 2.0),
    (12000.0, 82, 0.5), (12000.0, 82, 2.0),
    (11000.0, 60, 1.0), (11000.0, 60, 4.0), (10000.0, 60, 2.0),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--harness", type=pathlib.Path,
        default=REPO / "build-rel/tools/full_spectrum/g4_source_escape_campaign")
    parser.add_argument(
        "--output", type=pathlib.Path, required=True)
    parser.add_argument("--events", type=int, default=4000000)
    parser.add_argument("--bins", type=int, default=512)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument(
        "--offgrid", action="store_true",
        help="run the 11 off-grid validation cells instead of the bank")
    parser.add_argument(
        "--energies-kev", type=float, nargs="+", default=None)
    parser.add_argument("--zs", type=int, nargs="+", default=None)
    parser.add_argument("--taus", type=float, nargs="+", default=None)
    parser.add_argument(
        "--grid-json", type=pathlib.Path, default=None,
        help="JSON [{energy_keV, Z, tau}, ...] overriding the product grid")
    parser.add_argument(
        "--min-energy-by-z-json", type=pathlib.Path, default=None,
        help="JSON {Z: min_energy_keV} replacing the built-in importance "
             "threshold")
    return parser.parse_args()


def tag(value: float) -> str:
    return (f"{value:g}").replace(".", "p")


def seed_for(name: str, salt: str) -> int:
    digest = hashlib.sha256(f"{salt}:{name}".encode()).digest()
    return int.from_bytes(digest[:6], "big")


def main() -> int:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    explicit = (args.energies_kev is not None or args.zs is not None
                or args.taus is not None)
    if args.grid_json:
        grid = [
            (float(cell["energy_keV"]), int(cell["Z"]), float(cell["tau"]))
            for cell in json.loads(args.grid_json.read_text())]
    elif args.offgrid and not explicit:
        grid = list(OFFGRID_CELLS)
    else:
        minimum_by_z = dict(MIN_ENERGY_BY_Z_KEV) if not explicit else {}
        if args.min_energy_by_z_json:
            minimum_by_z = {
                int(z): float(e)
                for z, e in json.loads(
                    args.min_energy_by_z_json.read_text()).items()}
        grid = [
            (float(energy), z_value, tau)
            for energy in (args.energies_kev or DEFAULT_ENERGIES_KEV)
            for z_value in (args.zs or DEFAULT_ZS)
            for tau in (args.taus or DEFAULT_TAUS)
            if energy >= minimum_by_z.get(z_value, 0.0)]
    salt = "photonuclear-offgrid" if args.offgrid else "photonuclear-bank"
    print(f"{len(grid)} cells at {args.events} histories")

    def run_cell(cell) -> str:
        energy, z_value, tau = cell
        name = f"E{tag(energy)}_Z{z_value}_tau{tag(tau)}"
        output = args.output / f"{name}_geant4.csv"
        if output.exists():
            return f"skip {name}"
        log = args.output / f"{name}.log"
        # The harness writes its CSV at the end of the run; write it under a
        # temporary name so an interrupted cell is never mistaken for a
        # finished one.
        partial = args.output / f"{name}_geant4.csv.tmp"
        command = [
            str(args.harness),
            "--energy-kev", str(energy), "--z", str(z_value),
            "--tau", str(tau), "--events", str(args.events),
            "--bins", str(args.bins), "--seed", str(seed_for(name, salt)),
            "--output", str(partial)]
        with open(log, "w") as handle:
            result = subprocess.run(
                command, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode != 0 or not partial.exists():
            partial.unlink(missing_ok=True)
            return f"FAILED {name} (rc={result.returncode}; see {log})"
        partial.replace(output)
        return f"done {name}"

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for message in pool.map(run_cell, grid):
            print(message, flush=True)
            failures += message.startswith("FAILED")
    print(f"bank complete: {len(grid)} cells, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
