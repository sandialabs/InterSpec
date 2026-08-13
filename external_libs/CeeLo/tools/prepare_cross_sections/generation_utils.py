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

"""Small dependency-free helpers shared by locked cross-section generators."""

from __future__ import annotations

import io
import math
import os
import struct
import tempfile
from pathlib import Path
from typing import TextIO


OUT_T_keV = [
    10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0,
    300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 3000.0,
    4000.0, 5000.0, 6000.0, 8000.0, 10000.0, 15000.0, 20000.0,
]
OUT_KAPPA = [
    0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.97,
    0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1.0,
]
Z_MAX = 92
SECTION_START = (
    "// =====================================================================\n"
    "// Seltzer-Berger bremsstrahlung spectral data\n"
)
SECTION_END = (
    "// =====================================================================\n"
    "// Master element data array: g_element_data[92]\n"
)
ELEMENT_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
]


def as_float32(value: float) -> float:
    return struct.unpack("f", struct.pack("f", value))[0]


def format_float32(value: float) -> str:
    """Return the shortest stable C++ literal that round-trips a float32."""
    text = format(as_float32(value), ".9g")
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text + "f"


def format_rows(values, per_line, indent, formatter) -> str:
    lines = []
    for index in range(0, len(values), per_line):
        chunk = values[index:index + per_line]
        suffix = "," if index + per_line < len(values) else ""
        lines.append(indent + ", ".join(formatter(value) for value in chunk) + suffix)
    return "\n".join(lines)


def write_section(
    out: TextIO,
    tables: dict[int, list[list[float]]],
    encoding: str,
    provenance: str,
) -> None:
    if encoding not in ("float32", "uint16"):
        raise ValueError(f"unsupported bremsstrahlung encoding: {encoding}")
    n_kappa, n_energy = len(OUT_KAPPA), len(OUT_T_keV)
    out.write(SECTION_START)
    out.write(provenance)
    out.write(f"extern const uint16_t kSB_n_kappa = {n_kappa};\n")
    out.write(f"extern const uint16_t kSB_n_energy = {n_energy};\n\n")
    out.write(f"extern const float kSB_kappa[{n_kappa}] = {{\n"
              f"{format_rows(OUT_KAPPA, 6, '    ', format_float32)}\n}};\n\n")
    log_energy = [math.log10(value) for value in OUT_T_keV]
    out.write(f"extern const float kSB_log_E_keV[{n_energy}] = {{\n"
              f"{format_rows(log_energy, 6, '    ', format_float32)}\n}};\n\n")
    count = n_energy * n_kappa
    scales: dict[int, float] = {}
    if encoding == "uint16":
        for z in range(1, Z_MAX + 1):
            maximum = max(value for row in tables[z] for value in row)
            scales[z] = as_float32(maximum / 65535.0) if maximum > 0.0 else 0.0
        out.write(f"extern const float kSB_chi_scale[{Z_MAX}] = {{\n")
        out.write(format_rows([scales[z] for z in range(1, Z_MAX + 1)], 6,
                              "    ", format_float32))
        out.write("\n};\n\n")
    for z in range(1, Z_MAX + 1):
        flat = [value for row in tables[z] for value in row]
        if encoding == "uint16":
            scale = scales[z]
            packed = [
                min(65535, max(0, 0 if scale == 0.0 else round(value / scale)))
                for value in flat
            ]
            out.write(f"static const uint16_t Z{z}_sb_chi[{count}] = {{  // {ELEMENT_SYMBOLS[z - 1]}\n")
            out.write(format_rows(packed, 12, "    ", str))
        else:
            out.write(f"static const float Z{z}_sb_chi[{count}] = {{  // {ELEMENT_SYMBOLS[z - 1]}\n")
            out.write(format_rows(flat, 6, "    ", format_float32))
        out.write("\n};\n")
    out.write("\n")


def replace_section(
    path: Path,
    tables: dict[int, list[list[float]]],
    encoding: str,
    provenance: str,
) -> None:
    original = path.read_text(encoding="utf-8")
    start, end = original.find(SECTION_START), original.find(SECTION_END)
    if start < 0 or end <= start:
        raise ValueError(f"Seltzer-Berger markers not found in {path}")
    buffer = io.StringIO()
    write_section(buffer, tables, encoding, provenance)
    updated = original[:start] + buffer.getvalue() + original[end:]
    mode = path.stat().st_mode
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(updated)
        os.chmod(temporary, mode & 0o7777)
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
