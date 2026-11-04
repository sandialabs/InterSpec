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

"""Generate K-shell electron-impact ionization cross sections from locked EPICS2023 EEDL.

For every element Z = 1..92, writes EEDL's K-subshell electroionization cross section
(MF=23 MT=534) at its native nodes from threshold (the K binding energy) through the first
node at or above 20 MeV. EEDL declares linear-linear interpolation for these tables; the
runtime interpolates the stored float32 nodes the same way, and the generator gates that
representation against EEDL's own interpolation on 16 probes per interval. The detailed
electron walk samples K vacancies from it (electron-impact ionization x-rays).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

from epics_endf import EndfFile
from fetch_sources import load_lock, sha256_file
from generation_utils import ELECTRON_Z_MAX, ELEMENT_SYMBOLS, as_float32, format_float32

CEELO_LICENCE_HEADER = """\
/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC.
 This library is distributed under the GNU Lesser General Public License,
 version 2.1 or (at your option) any later version.
 */
"""

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DEFAULT_LOCK = HERE / "sources.lock.json"
DEFAULT_OUTPUT = ROOT / "src" / "cross_sections" / "electron_ionization_data.cpp"
DEFAULT_REPORT = HERE / "reports" / "electron_ionization_eedl.json"

MAX_ENERGY_EV = 20.0e6
PROBES_PER_INTERVAL = 16
MAX_RELATIVE_ERROR = 1.0e-4


def lin_lin(x: list[float], y: list[float], e: float) -> float:
    """Runtime interpolation of the stored nodes (mirrors electron_ionization_data.h)."""
    if e <= x[0]:
        return 0.0
    for i in range(len(x) - 1):
        if e <= x[i + 1]:
            t = (e - x[i]) / (x[i + 1] - x[i])
            return y[i] + t * (y[i + 1] - y[i])
    return y[-1]


def element_nodes(endf: EndfFile, z: int) -> tuple[list[float], list[float]]:
    table = endf.tab1(z, 23, 534)
    if any(law != 2 for _, law in table.regions):
        raise ValueError(f"Z={z}: MT=534 interpolation laws {table.regions} are not lin-lin")
    xs, ys = [], []
    for x, y in zip(table.x, table.y):
        xs.append(x)
        ys.append(y)
        if x >= MAX_ENERGY_EV:
            break
    if xs[-1] < MAX_ENERGY_EV:
        raise ValueError(f"Z={z}: MT=534 stops at {xs[-1]} eV")
    if ys[0] != 0.0:
        raise ValueError(f"Z={z}: MT=534 does not start from zero at threshold")
    return xs, ys


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--eedl", type=Path, required=True, help="EEDL2023.ALL from the source cache")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    lock = load_lock(args.lock)
    source = lock["sources"]["eedl2023"]
    digest = sha256_file(args.eedl)
    if digest != source["sha256"]:
        print(f"EEDL hash {digest} does not match the lock", file=sys.stderr)
        return 1

    endf = EndfFile(args.eedl)
    elements, report_rows = [], {}
    maximum_error = 0.0
    for z in range(1, ELECTRON_Z_MAX + 1):
        xs_eV, ys = element_nodes(endf, z)
        table = endf.tab1(z, 23, 534)
        energy_keV = [as_float32(x / 1000.0) for x in xs_eV]
        sigma_b = [as_float32(y) for y in ys]
        element_error = 0.0
        for i in range(len(xs_eV) - 1):
            for k in range(1, PROBES_PER_INTERVAL + 1):
                e_eV = xs_eV[i] + (xs_eV[i + 1] - xs_eV[i]) * k / (PROBES_PER_INTERVAL + 1)
                if e_eV > MAX_ENERGY_EV:
                    break
                truth = table.evaluate(e_eV)
                stored = lin_lin(energy_keV, sigma_b, e_eV / 1000.0)
                if truth > 0.0:
                    element_error = max(element_error, abs(stored / truth - 1.0))
        maximum_error = max(maximum_error, element_error)
        elements.append((z, energy_keV, sigma_b))
        report_rows[ELEMENT_SYMBOLS[z - 1]] = {
            "threshold_keV": energy_keV[0],
            "nodes": len(energy_keV),
            "peak_sigma_barn": max(sigma_b),
            "maximum_relative_error": element_error,
        }
    if maximum_error > MAX_RELATIVE_ERROR:
        raise ValueError(f"float32 representation error {maximum_error} exceeds {MAX_RELATIVE_ERROR}")

    with args.output.open("w", encoding="utf-8") as out:
        out.write(CEELO_LICENCE_HEADER)
        out.write(
            "/* AUTO-GENERATED by tools/prepare_cross_sections/generate_eedl_ionization.py.\n"
            " Direct source: EPICS2023 EEDL MF=23 MT=534 (K-subshell electroionization),\n"
            " verified by sources.lock.json. No Geant4/G4EMLOW data is read or emitted. */\n")
        out.write('#include "cross_sections/electron_ionization_data.h"\n')
        out.write("namespace ceelo {\nnamespace eedl_ionization {\n")
        out.write(f"static_assert(kElementCount == {ELECTRON_Z_MAX}, "
                  "\"regenerate electron_ionization_data.cpp\");\n")
        for z, energy, sigma in elements:
            symbol = ELEMENT_SYMBOLS[z - 1]
            out.write(f"static const float {symbol}_energy_keV[{len(energy)}] = {{"
                      + ", ".join(format_float32(v) for v in energy) + "};\n")
            out.write(f"static const float {symbol}_sigma_barn[{len(sigma)}] = {{"
                      + ", ".join(format_float32(v) for v in sigma) + "};\n")
        out.write("extern const KShellIonization kKShell[kElementCount] = {\n")
        for z, energy, _ in elements:
            symbol = ELEMENT_SYMBOLS[z - 1]
            out.write(f"    {{{len(energy)}, {symbol}_energy_keV, {symbol}_sigma_barn}}"
                      f"{',' if z < ELECTRON_Z_MAX else ''} // Z={z}\n")
        out.write("};\n} // namespace eedl_ionization\n} // namespace ceelo\n")

    node_count = sum(len(energy) for _, energy, _ in elements)
    report = {
        "source": "eedl2023",
        "source_sha256": digest,
        "source_lock_sha256": hashlib.sha256(args.lock.read_bytes()).hexdigest(),
        "reactions": ["MF=23 MT=534"],
        "interpolation": "linear-linear, as EEDL declares; zero below threshold",
        "maximum_energy_keV": MAX_ENERGY_EV / 1000.0,
        "maximum_relative_error": maximum_error,
        "acceptance_relative_error": MAX_RELATIVE_ERROR,
        "passed": maximum_error <= MAX_RELATIVE_ERROR,
        "retained_coefficients": 2 * node_count,
        "coefficient_bytes": 8 * node_count,
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "elements": report_rows,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output} and {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
