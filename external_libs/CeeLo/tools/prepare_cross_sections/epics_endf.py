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

"""Small, dependency-free ENDF-6 readers used by the EPICS generators.

Only the TAB1 records needed for EPDL MF=23/27 and the LIST records used by
EADL MF=28 are implemented. Keeping the parser here makes source preparation
reproducible without xraylib and keeps normal builds entirely offline.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from pathlib import Path


def endf_float(field: str) -> float:
    text = field.strip().replace("D", "E").replace("d", "e")
    if not text:
        return 0.0
    try:
        return float(text)
    except ValueError:
        for index in range(len(text) - 1, 0, -1):
            if text[index] in "+-" and (text[index - 1].isdigit() or text[index - 1] == "."):
                return float(text[:index] + "e" + text[index:])
        raise ValueError(f"invalid ENDF floating field {field!r}")


def endf_int(field: str) -> int:
    text = field.strip()
    return int(text) if text else 0


def fields(line: str) -> list[str]:
    return [line[index * 11:(index + 1) * 11] for index in range(6)]


def control(line: str) -> tuple[float, float, int, int, int, int]:
    values = fields(line)
    return (
        endf_float(values[0]),
        endf_float(values[1]),
        endf_int(values[2]),
        endf_int(values[3]),
        endf_int(values[4]),
        endf_int(values[5]),
    )


@dataclass(frozen=True)
class Tab1:
    mat: int
    mf: int
    mt: int
    x: tuple[float, ...]
    y: tuple[float, ...]
    regions: tuple[tuple[int, int], ...]  # (last one-based point, interpolation law)

    def law_for_interval(self, lower_index: int) -> int:
        upper_one_based = lower_index + 2
        for last, law in self.regions:
            if upper_one_based <= last:
                return law
        return self.regions[-1][1]

    def evaluate(self, argument: float) -> float:
        if argument <= self.x[0]:
            return self.y[0]
        if argument >= self.x[-1]:
            return self.y[-1]
        # bisect_right selects the upper side of an exactly duplicated edge.
        lower = bisect.bisect_right(self.x, argument) - 1
        if lower >= len(self.x) - 1:
            return self.y[-1]
        x0, x1 = self.x[lower], self.x[lower + 1]
        y0, y1 = self.y[lower], self.y[lower + 1]
        if x1 <= x0:
            return y1
        law = self.law_for_interval(lower)
        if law == 1:
            return y0
        if law == 2:
            t = (argument - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
        if law == 3 and x0 > 0.0 and argument > 0.0:
            t = math.log(argument / x0) / math.log(x1 / x0)
            return y0 + t * (y1 - y0)
        if law == 4 and y0 > 0.0 and y1 > 0.0:
            t = (argument - x0) / (x1 - x0)
            return math.exp(math.log(y0) + t * math.log(y1 / y0))
        if law == 5 and x0 > 0.0 and argument > 0.0 and y0 > 0.0 and y1 > 0.0:
            t = math.log(argument / x0) / math.log(x1 / x0)
            return math.exp(math.log(y0) + t * math.log(y1 / y0))
        t = (argument - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)


class EndfFile:
    def __init__(self, path: Path):
        self.path = path
        self._sections: dict[tuple[int, int, int], list[str]] = {}
        with path.open(encoding="latin-1") as stream:
            for raw in stream:
                line = raw.rstrip("\r\n")
                if len(line) < 75:
                    continue
                try:
                    mat = endf_int(line[66:70])
                    mf = endf_int(line[70:72])
                    mt = endf_int(line[72:75])
                except ValueError:
                    continue
                if mat > 0 and mf > 0 and mt > 0:
                    self._sections.setdefault((mat, mf, mt), []).append(line)

    def section(self, z: int, mf: int, mt: int) -> list[str]:
        key = (z * 100, mf, mt)
        try:
            return self._sections[key]
        except KeyError as error:
            raise KeyError(f"{self.path}: missing MAT/MF/MT={key}") from error

    def tab1(self, z: int, mf: int, mt: int) -> Tab1:
        lines = self.section(z, mf, mt)
        if len(lines) < 3:
            raise ValueError(f"Z={z} MF={mf} MT={mt}: truncated TAB1")
        _, _, _, _, nr, np = control(lines[1])
        if nr <= 0 or np <= 0:
            raise ValueError(f"Z={z} MF={mf} MT={mt}: invalid NR={nr}, NP={np}")
        interpolation_fields: list[int] = []
        interpolation_lines = (2 * nr + 5) // 6
        for line in lines[2:2 + interpolation_lines]:
            interpolation_fields.extend(endf_int(field) for field in fields(line))
        regions = tuple(
            (interpolation_fields[2 * index], interpolation_fields[2 * index + 1])
            for index in range(nr)
        )
        if regions[-1][0] != np or any(law not in (1, 2, 3, 4, 5) for _, law in regions):
            raise ValueError(f"Z={z} MF={mf} MT={mt}: malformed interpolation regions {regions}")

        data: list[float] = []
        for line in lines[2 + interpolation_lines:]:
            data.extend(endf_float(field) for field in fields(line))
        if len(data) < 2 * np:
            raise ValueError(f"Z={z} MF={mf} MT={mt}: expected {2*np} TAB1 values, got {len(data)}")
        x = tuple(data[index] for index in range(0, 2 * np, 2))
        y = tuple(data[index] for index in range(1, 2 * np, 2))
        if any(x[index] > x[index + 1] for index in range(np - 1)):
            raise ValueError(f"Z={z} MF={mf} MT={mt}: decreasing TAB1 argument")
        return Tab1(z * 100, mf, mt, x, y, regions)
