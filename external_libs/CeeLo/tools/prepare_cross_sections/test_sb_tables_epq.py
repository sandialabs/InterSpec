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

"""Focused tests for the source-safe EPQ bremsstrahlung transformation."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import sb_tables_epq as epq


class EpqTransformationTest(unittest.TestCase):
    def test_lagrange_reproduces_polynomials(self) -> None:
        quadratic = lambda x: 2.0 - 3.0 * x + 7.0 * x * x
        xs2 = [0.0, 0.05, 0.1]
        self.assertAlmostEqual(
            epq.lagrange(xs2, [quadratic(x) for x in xs2], 0.025),
            quadratic(0.025),
            places=14,
        )

        cubic = lambda x: 3.0 + 2.0 * x - 5.0 * x * x + 11.0 * x**3
        xs3 = [0.0, 0.05, 0.1, 0.15]
        self.assertAlmostEqual(
            epq.lagrange(xs3, [cubic(x) for x in xs3], 0.075),
            cubic(0.075),
            places=14,
        )

    def test_transform_only_synthesizes_two_nodes(self) -> None:
        row = [10.0 + 3.0 * kappa + 2.0 * kappa**2 for kappa in epq.EPQ_KAPPA]
        transformed = epq.transform_row(row)
        self.assertEqual(len(transformed), 32)
        output = dict(zip(epq.OUT_KAPPA, transformed))
        source = dict(zip(epq.EPQ_KAPPA, row))
        for kappa in epq.OUT_KAPPA:
            if kappa not in (0.025, 0.075):
                self.assertEqual(output[kappa], epq.as_float32(source[kappa]))
        self.assertNotIn(0.125, output)
        self.assertAlmostEqual(output[0.025], 10.0 + 3.0 * 0.025 + 2.0 * 0.025**2, places=6)
        self.assertAlmostEqual(output[0.075], 10.0 + 3.0 * 0.075 + 2.0 * 0.075**2, places=6)

    def test_parser_rejects_malformed_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pdebr01.tab"
            path.write_text("1\n1000 2 3\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "numeric fields"):
                epq.parse_epq_table(path, 1)

    def test_source_lock_pins_epq_table_set(self) -> None:
        lock = json.loads(epq.DEFAULT_LOCK.read_text(encoding="utf-8"))
        source = lock["sources"]["nist_epq"]
        self.assertEqual(source["version"], "git-ce1e5895892b17854d857e596250e4ae70d6b274")
        self.assertRegex(source["sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(source["table_set_sha256"], r"^[0-9a-f]{64}$")
        self.assertIn("No row normalization", source["transformation"])


if __name__ == "__main__":
    unittest.main()
