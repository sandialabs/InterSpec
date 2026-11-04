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

"""Generate electron elastic-scattering transport moments from locked EPICS2023 EEDL.

For every element Z = 1..92 and a log grid of electron energies, writes

    sigma_1 = integral (1 - mu) dsigma            (first transport cross section)
    R       = sigma_2 / sigma_1,  sigma_2 = integral (3/2)(1 - mu^2) dsigma

from the EEDL large-angle elastic cross section (MF=23 MT=525) and its angular
distributions (MF=26 MT=525). The detailed electron walk fits a Wentzel
(screened-Rutherford-shaped) distribution to (sigma_1, R) per material, which
reproduces the first two transport moments of the evaluated cross section. (A
relativistic large-angle suppression fitted to EEDL's backward cross section was
tried as a third parameter: with the moments matched it moved electron-beam
spectra by nothing measurable, so it was dropped.)

EEDL's small-angle part (mu > 0.999999) is excluded: each such collision moves
(1 - mu) by under 1e-6, so its share of sigma_1 is negligible.

The angular distributions are tabulated at only ~2 energies per decade, with a
gap from 256 keV to 10 MeV. Between them each transport cross section is
interpolated as its ratio to the screened-Rutherford one (Moliere screening),
ln(sigma_k / sigma_k^SR) linear in ln E: that ratio is in effect a Mott
correction factor and varies slowly (Ge: 1.00 at 256 keV, 0.97 at 10 MeV), where
the cross sections themselves fall by three decades. (Interpolating <1 - mu>
(pc)^2 instead put sigma_1 30-46% high at 1-3 MeV against GEANT4's
Goudsmit-Saunderson model; this reproduces it within 1% for Ge.)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

from epics_endf import EndfFile, control, endf_float, fields
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
DEFAULT_OUTPUT = ROOT / "src" / "cross_sections" / "electron_elastic_data.cpp"
DEFAULT_REPORT = HERE / "reports" / "electron_elastic_eedl.json"

MEC2_EV = 510998.95
GRID_MIN_KEV = 1.0
GRID_MAX_KEV = 20000.0
GRID_POINTS = 36


def grid_keV() -> list[float]:
    ratio = math.log(GRID_MAX_KEV / GRID_MIN_KEV) / (GRID_POINTS - 1)
    return [GRID_MIN_KEV * math.exp(ratio * k) for k in range(GRID_POINTS)]


def angular_tables(endf: EndfFile, z: int) -> list[tuple[float, list[float], list[float]]]:
    """(energy eV, mu, p(mu)) for each MF=26 MT=525 incident energy."""
    lines = endf.section(z, 26, 525)
    # HEAD, then the product-yield TAB1 (control, interpolation, data).
    _, _, _, _, nr, n_points = control(lines[1])
    position = 2 + (2 * nr + 5) // 6 + (2 * n_points + 5) // 6
    # TAB2 over incident energy (LAW=2 discrete two-body angular tables).
    _, _, _, _, nr2, n_energies = control(lines[position])
    position += 1 + (2 * nr2 + 5) // 6
    tables = []
    for _ in range(n_energies):
        _, energy_eV, _, _, n_words, n_pairs = control(lines[position])
        position += 1
        n_lines = (n_words + 5) // 6
        values: list[float] = []
        for line in lines[position:position + n_lines]:
            values.extend(endf_float(field) for field in fields(line))
        position += n_lines
        if n_words != 2 * n_pairs:
            raise ValueError(f"Z={z} E={energy_eV}: NW={n_words} != 2*NL={2 * n_pairs}")
        tables.append((energy_eV, values[0:2 * n_pairs:2], values[1:2 * n_pairs:2]))
    return tables


def moments(mu: list[float], p: list[float]) -> tuple[float, float, float]:
    """Norm, <1 - mu> and <1 - mu^2> of a lin-lin tabulated p(mu) (Simpson: exact here)."""
    norm = m1 = m2 = 0.0
    for j in range(len(mu) - 1):
        a, b = mu[j], mu[j + 1]
        pa, pb = p[j], p[j + 1]
        c, pc = 0.5 * (a + b), 0.5 * (pa + pb)
        h = (b - a) / 6.0
        norm += h * (pa + 4.0 * pc + pb)
        m1 += h * (pa * (1 - a) + 4.0 * pc * (1 - c) + pb * (1 - b))
        m2 += h * (pa * (1 - a * a) + 4.0 * pc * (1 - c * c) + pb * (1 - b * b))
    return norm, m1 / norm, m2 / norm


ALPHA = 1.0 / 137.035999
RE2_BARN = 7.9407877e-26 * 1e24


def screened_rutherford(z: int, energy_eV: float) -> tuple[float, float]:
    """sigma_1, sigma_2 (barn) of screened-Rutherford scattering, Moliere screening.

    dsigma/dOmega = Z(Z+1) r_e^2 (m c^2 E)^2 / (pc)^4 / (x + 2 eta)^2, x = 1 - mu;
    only the smooth energy dependence matters (it is divided out and back in).
    """
    t_keV, me_keV = energy_eV / 1e3, MEC2_EV / 1e3
    pc = math.sqrt(t_keV * (t_keV + 2.0 * me_keV))
    e_tot = t_keV + me_keV
    b2 = (pc / e_tot) ** 2
    eta = (2.1067 * z ** (1.0 / 3.0) / pc) ** 2 * (1.13 + 3.76 * (ALPHA * z) ** 2 / b2)
    c = 2.0 * math.pi * z * (z + 1) * RE2_BARN * (me_keV * e_tot) ** 2 / pc ** 4
    a = 2.0 * eta
    s1 = c * (math.log((2.0 + a) / a) + a / (2.0 + a) - 1.0)
    s2 = c * 1.5 * (-2.0 + (2.0 + 2.0 * a) * math.log((2.0 + a) / a)
                    - a * (2.0 + a) * (1.0 / a - 1.0 / (2.0 + a)))
    return s1, s2


def element_rows(endf: EndfFile, z: int, energies_keV: list[float]):
    sigma_la = endf.tab1(z, 23, 525)
    nodes = []
    for energy_eV, mu, p in angular_tables(endf, z):
        if energy_eV < 100.0:
            continue  # below every grid energy; keeps the interpolation local
        norm, m1, m2 = moments(mu, p)
        if abs(norm - 1.0) > 0.02:
            raise ValueError(f"Z={z} E={energy_eV}: angular table norm {norm}")
        R = 1.5 * m2 / m1
        s1 = sigma_la.evaluate(energy_eV) * m1
        sr1, sr2 = screened_rutherford(z, energy_eV)
        nodes.append((math.log(energy_eV), math.log(s1 / sr1), math.log(s1 * R / sr2)))
    sigma1, ratio = [], []
    for e_keV in energies_keV:
        e_eV = e_keV * 1e3
        x = math.log(e_eV)
        k = max(0, min(len(nodes) - 2, next((i for i in range(len(nodes) - 1)
                                             if nodes[i + 1][0] >= x), len(nodes) - 2)))
        (x0, f0, g0), (x1, f1, g1) = nodes[k], nodes[k + 1]
        t = (x - x0) / (x1 - x0)
        sr1, sr2 = screened_rutherford(z, e_eV)
        s1 = math.exp(f0 + t * (f1 - f0)) * sr1
        s2 = math.exp(g0 + t * (g1 - g0)) * sr2
        sigma1.append(s1)  # barns
        ratio.append(s2 / s1)
    return sigma1, ratio, nodes


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
    energies = grid_keV()
    sigma1_rows, ratio_rows, report_rows = [], [], {}
    for z in range(1, ELECTRON_Z_MAX + 1):
        s1, r, nodes = element_rows(endf, z, energies)
        sigma1_rows.append([as_float32(v) for v in s1])
        ratio_rows.append([as_float32(v) for v in r])
        report_rows[ELEMENT_SYMBOLS[z - 1]] = {
            "nodes_keV": [round(math.exp(x) / 1e3, 6) for x, _, _ in nodes],
            "ratio_range": [min(r), max(r)],
            "sigma1_barn_at_100keV": s1[min(range(len(energies)), key=lambda i: abs(energies[i] - 100))],
        }

    with args.output.open("w", encoding="utf-8") as out:
        out.write(CEELO_LICENCE_HEADER)
        out.write(
            "/* AUTO-GENERATED by tools/prepare_cross_sections/generate_eedl_elastic.py.\n"
            " Direct source: EPICS2023 EEDL MF=23/MF=26 MT=525 (large-angle elastic),\n"
            " verified by sources.lock.json. No Geant4/G4EMLOW data is read or emitted. */\n")
        out.write('#include "cross_sections/electron_elastic_data.h"\n')
        out.write("namespace ceelo {\nnamespace eedl_elastic {\n")
        out.write(f"static_assert(kElementCount == {ELECTRON_Z_MAX} && kEnergyCount == {GRID_POINTS},\n"
                  "              \"regenerate electron_elastic_data.cpp\");\n")
        out.write("extern const float kEnergy_keV[kEnergyCount] = {")
        out.write(", ".join(format_float32(as_float32(e)) for e in energies))
        out.write("};\n")
        for name, rows in (("kSigma1_barn", sigma1_rows), ("kRatio21", ratio_rows)):
            out.write(f"extern const float {name}[kElementCount][kEnergyCount] = {{\n")
            for z, row in enumerate(rows, start=1):
                out.write("    {" + ", ".join(format_float32(v) for v in row)
                          + f"}}{',' if z < ELECTRON_Z_MAX else ''} // {ELEMENT_SYMBOLS[z - 1]}\n")
            out.write("};\n")
        out.write("} // namespace eedl_elastic\n} // namespace ceelo\n")

    report = {
        "source": "eedl2023",
        "source_sha256": digest,
        "source_lock_sha256": hashlib.sha256(args.lock.read_bytes()).hexdigest(),
        "reactions": ["MF=23 MT=525", "MF=26 MT=525"],
        "model": "Wentzel fitted per material to sigma_1 and sigma_2",
        "grid_keV": [GRID_MIN_KEV, GRID_MAX_KEV, GRID_POINTS],
        "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "elements": report_rows,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output} and {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
