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

"""Run CeeLo and Geant4 source-escape truth cells (the acceptance anchors).

Each anchor is a centered isotropic point source in a homogeneous
equivalent-element sphere.  Per cell this writes, under OUTPUT/cells/:

* `<id>_ceelo.csv.gz` -- `ceelo_source_escape_campaign` (Rayleigh and Compton
  Doppler on, fine Compton-order labels);
* `<id>_geant4.csv` -- `g4_source_escape_campaign` (Livermore EM physics,
  photonuclear on), the independent reference; skipped with --no-geant4;
* `<id>_model.csv.gz` -- only with --payload: the functional model's spectrum;

plus `<id>_manifest.json`.  Seeds derive from the cell, and finished legs are
skipped, so an interrupted grid resumes where it stopped.

`--acceptance-grid` runs the standing 150-cell grid, {59, 100, 250, 661,
1000, 3000} keV x Z {6, 13, 26, 48, 82} x {0, 1, 5, 15, 50} g/cm^2, which
no fit selects on: it is the independent end-to-end arbiter of a payload.

Usage:
    run_source_anchor.py --output DIR --acceptance-grid [--workers N]
    run_source_anchor.py --output DIR --energy-kev E --z Z
        (--tau T | --areal-density G) [--payload FILE]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import bank_io
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "build-rel/tools/full_spectrum"
ACCEPTANCE_ENERGIES_KEV = (59.0, 100.0, 250.0, 661.0, 1000.0, 3000.0)
ACCEPTANCE_Z = (6, 13, 26, 48, 82)
ACCEPTANCE_AREAL_DENSITY = (0.0, 1.0, 5.0, 15.0, 50.0)


def atomic_stdout(command: list[str], output: Path) -> None:
    completed = subprocess.run(
        command, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bank_io.write_cell(output, completed.stdout)


def stable_seed(label: str, energy: float, z_value: int, amount: float) -> int:
    encoded = (f"source-anchor|{label}|{energy:.12g}|"
               f"{z_value}|{amount:.12g}").encode()
    return (int.from_bytes(hashlib.sha256(encoded).digest()[:8], "big")
            & ((1 << 31) - 1))


def identity(energy: float, z_value: int, amount_name: str,
             amount: float) -> str:
    kind = "tau" if amount_name == "--tau" else "rhoA"
    return f"E{energy:g}_Z{z_value}_{kind}{amount:g}".replace(".", "p")


def run_anchor(args: argparse.Namespace, energy: float, z_value: int,
               amount_name: str, amount: float) -> str:
    cells = args.output / "cells"
    cells.mkdir(parents=True, exist_ok=True)
    name = identity(energy, z_value, amount_name, amount)
    ceelo_path = bank_io.written_path(cells, name, "ceelo")
    model_path = bank_io.written_path(cells, name, "model")
    geant4_path = cells / f"{name}_geant4.csv"
    common = ["--energy-kev", repr(energy), "--z", str(z_value),
              amount_name, repr(amount), "--bins", str(args.bins)]
    if args.force or bank_io.existing(ceelo_path) is None:
        atomic_stdout([
            str(args.ceelo_driver), *common,
            "--events", str(args.ceelo_events),
            "--seed", str(stable_seed("ceelo", energy, z_value, amount)),
            "--fine-interaction-orders", "--source-compton-doppler",
        ], ceelo_path)
    if args.payload is not None and (args.force or bank_io.existing(model_path) is None):
        atomic_stdout([
            str(args.model_driver), *common, "--payload", str(args.payload),
        ], model_path)
    if args.no_geant4:
        pass
    elif args.geant4_input is not None:
        if args.force or not geant4_path.exists():
            shutil.copyfile(args.geant4_input, geant4_path)
    elif args.force or not geant4_path.exists():
        partial = cells / f"{name}_geant4.csv.tmp"
        with (cells / f"{name}_geant4.log").open("w") as log:
            subprocess.run([
                str(args.geant4_driver), *common,
                "--events", str(args.geant4_events),
                "--seed", str(stable_seed("geant4", energy, z_value, amount)),
                # The CeeLo banks carry shield K/L x-rays with impact
                # ionization at every energy; GEANT4 needs PIXE for the latter.
                "--fluorescence", "--pixe",
                "--output", str(partial),
            ], check=True, stdout=log, stderr=subprocess.STDOUT)
        partial.replace(geant4_path)
    manifest = {
        "schema": 1,
        "identity": name,
        "energy_keV": energy,
        "Z": z_value,
        amount_name.lstrip("-").replace("-", "_"): amount,
        "bins": args.bins,
        "ceelo_events": args.ceelo_events,
        "geant4_events": (None if args.geant4_input is not None
                          or args.no_geant4 else args.geant4_events),
        "physics": {
            "rayleigh_transport": True,
            "scatter_order_counts_compton_only": True,
            "compton_doppler": True,
            "positron_transport": True,
            "fluorescence": True,
            "impact_ionization": True,
            "geant4_em": "Livermore",
            "geant4_photonuclear": True,
        },
        "payload": str(args.payload) if args.payload is not None else None,
    }
    (args.output / f"{name}_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    return name


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--energy-kev", type=float)
    parser.add_argument("--z", type=int)
    amount = parser.add_mutually_exclusive_group()
    amount.add_argument("--tau", type=float)
    amount.add_argument("--areal-density", type=float)
    parser.add_argument(
        "--acceptance-grid", action="store_true",
        help="run the 150-cell acceptance grid")
    parser.add_argument(
        "--cell", action="append", default=[], metavar="ENERGY,Z,AREAL",
        help="run explicit (energy, Z, areal density) cells; repeatable")
    parser.add_argument("--bins", type=int, default=1024)
    parser.add_argument("--ceelo-events", type=int, default=500000)
    parser.add_argument("--geant4-events", type=int, default=100000)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--no-geant4", action="store_true",
        help="run only the CeeLo (and model) legs; no GEANT4 reference")
    parser.add_argument(
        "--geant4-input", type=Path,
        help="use this Geant4 CSV instead of running the harness "
             "(single-cell mode)")
    parser.add_argument(
        "--ceelo-driver", type=Path,
        default=BIN / "ceelo_source_escape_campaign")
    parser.add_argument(
        "--model-driver", type=Path,
        default=BIN / "ceelo_source_model_evaluate")
    parser.add_argument(
        "--geant4-driver", type=Path,
        default=BIN / "g4_source_escape_campaign")
    parser.add_argument(
        "--payload", type=Path,
        help="also evaluate this source-escape payload in every cell")
    args = parser.parse_args()

    jobs: list[tuple[float, int, str, float]] = []
    if args.acceptance_grid:
        jobs = [(energy, z_value, "--areal-density", rho)
                for energy in ACCEPTANCE_ENERGIES_KEV
                for z_value in ACCEPTANCE_Z
                for rho in ACCEPTANCE_AREAL_DENSITY]
    for text in args.cell:
        fields = text.split(",")
        if len(fields) != 3:
            parser.error(f"--cell must be ENERGY,Z,AREAL, received: {text}")
        jobs.append((float(fields[0]), int(fields[1]), "--areal-density",
                     float(fields[2])))
    if args.energy_kev is not None:
        if args.z is None or (args.tau is None and args.areal_density is None):
            parser.error("a single cell needs --z and --tau or "
                         "--areal-density")
        jobs.append((args.energy_kev, args.z,
                     "--tau" if args.tau is not None else "--areal-density",
                     args.tau if args.tau is not None
                     else args.areal_density))
    if not jobs:
        parser.error("give --acceptance-grid, --cell, or a single cell")
    if args.geant4_input is not None and len(jobs) != 1:
        parser.error("--geant4-input needs exactly one cell")

    args.output.mkdir(parents=True, exist_ok=True)
    geant4 = ("no GEANT4 leg" if args.no_geant4
              else f"GEANT4 {args.geant4_events}")
    print(f"{len(jobs)} anchor cells: CeeLo {args.ceelo_events} histories, "
          f"{geant4}, {args.bins} bins")
    failures = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_anchor, args, *job) for job in jobs]
        for index, future in enumerate(futures, 1):
            try:
                name = future.result()
            except subprocess.CalledProcessError as error:
                failures += 1
                print(f"FAILED {' '.join(map(str, error.cmd[:8]))} ... "
                      f"(exit {error.returncode})", flush=True)
                continue
            if index % 10 == 0 or index == len(futures):
                print(f"{index}/{len(futures)} cells; latest {name}",
                      flush=True)
    if failures:
        print(f"{failures} cell(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
