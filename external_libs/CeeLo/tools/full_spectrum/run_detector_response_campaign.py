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

"""Run a reproducible HPGe parallel-beam detector-response campaign.

The detector is a built-in (`--detector hpge_coax`) or any geometry file
(`--geometry FILE --detector LABEL`); the manifest records its label and
geometry text, which the fitter embeds in the payload, so a bank cannot be fit
as another detector. `--incidences iso,front,back,side` writes an arrival bank
instead: each cell is an isotropic field over a band of direction cosines
(cell_kind `detector_arrival`, which the fitter refuses, so it can only score). Cells are written under the output directory (keep it
outside the source tree, or in the git-ignored bake work area). This driver
covers only the detector kernel, so source-escape banks do not have to be
regenerated while the detector representation is refined.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import bank_io
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_TRAIN_ENERGIES_KEV = (
    10.0, 20.0, 40.0, 59.54, 80.0, 122.0, 200.0, 356.013, 511.0,
    661.657, 1000.0, 1023.0, 1332.492, 2000.0, 2614.511, 4000.0,
    6000.0, 9000.0, 12000.0,
)
DEFAULT_TRAIN_ANGLES = (
    -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0,
)


def number_list(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


ARRIVAL_CLASSES = ("iso", "front", "back", "side")


def identity(energy_keV: float, cosine: float) -> str:
    cosine_text = f"{cosine:+.5f}".replace("+", "p").replace("-", "m")
    return f"hpge_E{energy_keV:.6f}_c{cosine_text}"


def arrival_identity(energy_keV: float, incidence: str) -> str:
    return f"hpge_E{energy_keV:.6f}_{incidence}"


def incidence_list(text: str) -> list[str]:
    names = [name.strip() for name in text.split(",") if name.strip()]
    for name in names:
        if name not in ARRIVAL_CLASSES:
            raise argparse.ArgumentTypeError(
                f"{name!r} is not one of {', '.join(ARRIVAL_CLASSES)}")
    return names


def seed_for(identity_text: str, salt: str) -> int:
    digest = hashlib.sha256(f"{salt}:{identity_text}".encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def detector_options(detector: str, geometry: Path | None) -> list[str]:
    options = ["--detector", detector]
    if geometry is not None:
        options += ["--geometry", str(geometry)]
    return options


def run_cell(
    executable: Path,
    detector: list[str],
    output: Path,
    energy_keV: float,
    cosine: float | None,
    incidence: str | None,
    bins: int,
    min_events: int,
    max_events: int,
    seed: int,
    force: bool,
    extra_flags: list[str] | None = None,
) -> str:
    if bank_io.existing(output) is not None and not force:
        return "cached"
    completed = subprocess.run(
        [
            str(executable),
            "--energy-kev", f"{energy_keV:.17g}",
            *(["--incidence", incidence] if incidence is not None
              else ["--cos-theta", f"{cosine:.17g}"]),
            *detector,
            "--bins", str(bins),
            "--min-events", str(min_events),
            "--max-events", str(max_events),
            "--seed", str(seed),
            *(extra_flags or []),
        ],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    bank_io.write_cell(output, completed.stdout)
    return "ran"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--executable", type=Path,
        default=Path(
            "build-rel/tools/full_spectrum/"
            "ceelo_detector_response_campaign"))
    parser.add_argument(
        "--energies", type=number_list,
        default=list(DEFAULT_TRAIN_ENERGIES_KEV))
    angles = parser.add_mutually_exclusive_group()
    angles.add_argument(
        "--angles", type=number_list,
        default=list(DEFAULT_TRAIN_ANGLES))
    angles.add_argument(
        "--incidences", type=incidence_list,
        help="arrival classes (iso,front,back,side) instead of beam angles")
    parser.add_argument(
        "--detector", default="hpge_coax",
        help="built-in detector, or the label of the --geometry detector")
    parser.add_argument(
        "--geometry", type=Path,
        help="CEELO_GEOMETRY file describing the detector")
    parser.add_argument("--bins", type=int, default=256)
    parser.add_argument("--min-events", type=int, default=200_000)
    parser.add_argument("--max-events", type=int, default=500_000)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--seed-salt", default="hpge-training-v1")
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--detector-brems", action="store_true",
        help="enable detector-side bremsstrahlung (matches Geant4)")
    parser.add_argument(
        "--detector-fluorescence", action="store_true",
        help="enable detector-side fluorescence (matches Geant4)")
    parser.add_argument(
        "--detector-spectrum-fidelity", action="store_true",
        help="Compton vacancy x-rays, the Detailed electron walk and its "
             "electron-impact ionization x-rays")
    args = parser.parse_args()
    detector = detector_options(args.detector, args.geometry)
    geometry_text = subprocess.run(
        [str(args.executable), *detector, "--print-geometry"],
        check=True, text=True, stdout=subprocess.PIPE).stdout

    if args.incidences:
        cells = [
            {
                "identity": arrival_identity(energy, incidence),
                "energy_keV": energy,
                "incidence": incidence,
                "seed": seed_for(arrival_identity(energy, incidence),
                                 args.seed_salt),
            }
            for energy in args.energies
            for incidence in args.incidences
        ]
    else:
        cells = [
            {
                "identity": identity(energy, cosine),
                "energy_keV": energy,
                "cos_theta": cosine,
                "seed": seed_for(identity(energy, cosine), args.seed_salt),
            }
            for energy in args.energies
            for cosine in args.angles
        ]
    args.output.mkdir(parents=True, exist_ok=True)
    cell_directory = args.output / "cells"
    cell_directory.mkdir(exist_ok=True)
    manifest_path = args.output / "campaign_manifest.json"
    if manifest_path.exists() and not args.force:
        # Finished cells are reused, so they must be of this detector.
        previous = json.loads(manifest_path.read_text())
        if (previous.get("geometry_text", geometry_text) != geometry_text
                or previous.get("detector") != args.detector):
            parser.error(
                f"{args.output} holds cells of another detector "
                f"({previous.get('detector')}); use a new directory or --force")
    # The fitter and the scorers read every cell in the directory, so a cell
    # left from another energy grid would join this bank unseen -- and could
    # sit on an energy another bank holds out.
    planned = {cell["identity"] for cell in cells}
    stray = sorted(
        bank_io.cell_id(path) for path in bank_io.cell_files(cell_directory)
        if bank_io.cell_id(path) not in planned)
    if stray:
        parser.error(
            f"{cell_directory} holds {len(stray)} cell(s) this campaign does "
            f"not plan (e.g. {stray[0]}); move them out or use a new "
            "directory")
    manifest = {
        "schema": 1,
        "campaign": "FullSpecApprox HPGe detector response",
        "detector": args.detector,
        "geometry_text": geometry_text,
        "normalization": "per photon incident on detector envelope",
        "energy_range_keV": [
            min(args.energies), max(args.energies)],
        "energies_keV": args.energies,
        "cos_theta": [] if args.incidences else args.angles,
        "incidences": args.incidences or [],
        "bins": args.bins,
        "min_events": args.min_events,
        "max_events": args.max_events,
        "seed_salt": args.seed_salt,
        "physics": {
            "CeeLo": True,
            "fluorescence": False,
            "pair_production": True,
            "rayleigh": True,
            "bremsstrahlung": False,
        },
        "cells": cells,
        "detector_brems": args.detector_brems,
        "detector_fluorescence": args.detector_fluorescence,
        "detector_spectrum_fidelity": args.detector_spectrum_fidelity,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    extra_flags = []
    if args.detector_brems:
        extra_flags.append("--detector-brems")
    if args.detector_fluorescence:
        extra_flags.append("--detector-fluorescence")
    if args.detector_spectrum_fidelity:
        extra_flags.append("--detector-spectrum-fidelity")

    def submit(cell: dict) -> tuple[str, str]:
        destination = bank_io.written_path(cell_directory, cell['identity'])
        status = run_cell(
            args.executable, detector, destination,
            cell["energy_keV"], cell.get("cos_theta"),
            cell.get("incidence"), args.bins,
            args.min_events, args.max_events, cell["seed"], args.force,
            extra_flags)
        return cell["identity"], status

    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(submit, cell) for cell in cells]
        for future in as_completed(futures):
            name, status = future.result()
            completed_count += 1
            print(f"[{completed_count}/{len(cells)}] {name}: {status}",
                  flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
