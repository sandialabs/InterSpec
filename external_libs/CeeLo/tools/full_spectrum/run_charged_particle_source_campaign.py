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

"""Generate a reusable CeeLo mono-electron/positron source-escape bank.

One call runs one tensor product energy x element x range ratio x particle;
a bank made of several sub-grids is several calls into the same --output.
Thickness is sampled in units of the charged particle's CSDA range. This
places nodes around the physically important containment transition while the
CSV retains the corresponding areal density in g/cm2.

Each cell is `cells/T<T>_Z<Z>_r<r>_<particle>_ceelo.csv.gz`, written through a
temporary file so an interrupted run never leaves a truncated cell behind;
--skip-existing then resumes a bank cell by cell.  The CPU time of every cell
run is appended to `timing.csv` in the output directory.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import bank_io
import hashlib
import os
import subprocess
import sys
import threading
from pathlib import Path


def parse_numbers(text: str, conversion: type) -> list:
    return [conversion(value) for value in text.split(",") if value]


def stable_seed(identity: str) -> int:
    digest = hashlib.sha256(identity.encode()).digest()
    return int.from_bytes(digest[:8], "little") & 0x7FFFFFFFFFFFFFFF


def load_ranges(
    helper: Path,
    energies: list[float],
    elements: list[int],
) -> dict[tuple[float, int], tuple[float, float]]:
    queries = "".join(
        f"{energy:.17g} {z_value} 0\n"
        for energy in energies for z_value in elements)
    completed = subprocess.run(
        [str(helper)], input=queries, text=True,
        stdout=subprocess.PIPE, check=True)
    rows = list(csv.DictReader(completed.stdout.splitlines()))
    if len(rows) != len(energies) * len(elements):
        raise RuntimeError("cross-section helper returned the wrong row count")
    answer = {}
    for row in rows:
        answer[(float(row["energy_keV"]), int(row["Z"]))] = (
            float(row["electron_range_g_cm2"]),
            float(row["positron_range_g_cm2"]),
        )
    return answer


def energy_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--executable", type=Path,
        default=Path(
            "build-rel/tools/full_spectrum/"
            "ceelo_source_electron_campaign"))
    parser.add_argument(
        "--cross-section-helper", type=Path,
        default=Path(
            "build-rel/tools/full_spectrum/"
            "ceelo_source_cross_section_features"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--energies-kev",
        default="100,200,300,500,1000,2290,5000,12000")
    parser.add_argument(
        "--elements", default="6,13,26,35,48,64,82,92")
    parser.add_argument(
        "--range-ratios", default="0.05,0.2,0.5,1,2,5")
    parser.add_argument("--events", type=int, default=30000)
    parser.add_argument(
        "--seed-salt", default="",
        help="optional deterministic salt for an independent campaign seed")
    parser.add_argument("--bins", type=int, default=512)
    parser.add_argument("--minimum-photon-kev", type=float, default=10.0)
    parser.add_argument(
        "--jobs", type=int, default=1,
        help="number of independent campaign cells to run concurrently")
    parser.add_argument("--electrons-only", action="store_true")
    parser.add_argument("--positrons-only", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    if args.electrons_only and args.positrons_only:
        raise ValueError("particle filters are mutually exclusive")
    if args.jobs < 1:
        raise ValueError("jobs must be positive")
    energies = parse_numbers(args.energies_kev, float)
    elements = parse_numbers(args.elements, int)
    ratios = parse_numbers(args.range_ratios, float)
    ranges = load_ranges(args.cross_section_helper, energies, elements)
    particles = (
        ["electron"] if args.electrons_only
        else ["positron"] if args.positrons_only
        else ["electron", "positron"])
    cells = args.output / "cells"
    cells.mkdir(parents=True, exist_ok=True)

    total = len(energies) * len(elements) * len(ratios) * len(particles)
    commands: list[tuple[str, list[str], Path]] = []
    for energy in energies:
        for z_value in elements:
            electron_range, positron_range = ranges[(energy, z_value)]
            for particle in particles:
                particle_range = (
                    positron_range
                    if particle == "positron" else electron_range)
                for ratio in ratios:
                    identity = (
                        f"T{energy_label(energy)}_Z{z_value}_"
                        f"r{energy_label(ratio)}_{particle}")
                    output = bank_io.written_path(cells, identity)
                    if (
                        args.skip_existing
                        and (finished := bank_io.existing(output)) is not None
                        and finished.stat().st_size > 100
                    ):
                        commands.append((identity, [], output))
                        continue
                    areal_density = ratio * particle_range
                    maximum_energy = energy + (
                        510.99895 if particle == "positron" else 0.0)
                    command = [
                        str(args.executable),
                        "--energy-kev", f"{energy:.17g}",
                        "--z", str(z_value),
                        "--areal-density", f"{areal_density:.17g}",
                        "--events", str(args.events),
                        "--seed", str(stable_seed(
                            identity
                            if not args.seed_salt
                            else f"{args.seed_salt}|{identity}")),
                        "--bins", str(args.bins),
                        "--max-energy-kev", f"{maximum_energy:.17g}",
                        "--min-photon-kev",
                        f"{args.minimum_photon_kev:.17g}",
                        "--output", str(output) + ".part",
                    ]
                    if particle == "positron":
                        command.append("--positron")
                    commands.append((identity, command, output))

    timing_path = args.output / "timing.csv"
    timing_lock = threading.Lock()
    totals = {"run": 0, "reused": 0, "cpu_s": 0.0}

    def run_cell(task: tuple[int, tuple[str, list[str], Path]]) -> str:
        index, (identity, command, output) = task
        if not command:
            with timing_lock:
                totals["reused"] += 1
            return f"[{index}/{total}] {identity} (exists, reused)"
        process = subprocess.Popen(command)
        # wait4 reaps this child alone, so its CPU time is exact even with
        # other cells running concurrently.
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
        if process.returncode != 0:
            Path(str(output) + ".part").unlink(missing_ok=True)
            raise RuntimeError(
                f"cell {identity} failed (exit {process.returncode}); "
                "command:\n  " + " ".join(command))
        part = Path(str(output) + ".part")
        bank_io.write_cell(output, part.read_text())
        part.unlink()
        cpu_s = usage.ru_utime + usage.ru_stime
        with timing_lock:
            totals["run"] += 1
            totals["cpu_s"] += cpu_s
            new_file = not timing_path.exists()
            with timing_path.open("a") as stream:
                if new_file:
                    stream.write("identity,events,cpu_s\n")
                stream.write(f"{identity},{args.events},{cpu_s:.4f}\n")
        return f"[{index}/{total}] {identity} ({cpu_s:.2f} s CPU)"

    indexed = list(enumerate(commands, start=1))
    try:
        if args.jobs == 1:
            for task in indexed:
                print(run_cell(task), flush=True)
        else:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.jobs
            ) as executor:
                futures = [
                    executor.submit(run_cell, task) for task in indexed]
                try:
                    for future in concurrent.futures.as_completed(futures):
                        print(future.result(), flush=True)
                except RuntimeError:
                    # Fail fast: cells already running finish (and are
                    # kept), queued ones are dropped.
                    for future in futures:
                        future.cancel()
                    raise
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    rate = (
        f", {totals['run'] * args.events / totals['cpu_s']:.0f} "
        "events per CPU-second" if totals["cpu_s"] > 0.0 else "")
    print(
        f"{args.output}: {total} cells in this grid, {totals['run']} run "
        f"({totals['cpu_s']:.1f} s CPU{rate}), {totals['reused']} reused",
        flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
