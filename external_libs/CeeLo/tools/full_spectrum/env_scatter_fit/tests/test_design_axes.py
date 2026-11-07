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

import math
import unittest

from tools.full_spectrum.env_scatter_fit.axes import (
    conservative_rebin_2d,
    overlap_matrix,
    primary_registered_axis,
    uniform_mu_edges,
)
from tools.full_spectrum.env_scatter_fit.design import (
    DesignBounds,
    PhysicalNode,
    generate_candidate_pool,
    geometry_transform,
    scrambled_sobol,
    select_disjoint_maximin,
    stable_node_id_from_cm,
    transformed_domain,
)


class DesignTests(unittest.TestCase):
    def test_node_id_cross_language_fixture(self):
        config_hash = "6b9e67799f16dd979ae2913272145feffb9ac83bf5fd3bf5eebb71fe86e11530"
        node = PhysicalNode(662.0, 10.0, 1.0, 1.0)
        self.assertEqual(node.stable_id_hex(config_hash), "e00a5bb3240138d0")

    def test_node_id_uses_stored_centimetres_without_round_trip(self):
        config_hash = "d465fe820aac54179527eadc3caded1e922276bed97b90ab98a9709c8679ff2e"
        self.assertEqual(
            stable_node_id_from_cm(
                config_hash,
                2337.8461284396781,
                360.20485631067271,
                51.505621833298378,
                401.7725328733622,
            ),
            0xCA3EF7B30817B0CF,
        )

    def test_sobol_is_deterministic_and_bounded(self):
        first = scrambled_sobol(32, seed=1234)
        self.assertEqual(first, scrambled_sobol(32, seed=1234))
        self.assertNotEqual(first, scrambled_sobol(32, seed=1235))
        self.assertEqual(len(set(first)), 32)
        self.assertTrue(all(0.0 <= value < 1.0 for row in first for value in row))

    def test_candidate_validity_and_nonreciprocal_transform(self):
        bounds = DesignBounds(rho_min_m=0.02)
        pool = generate_candidate_pool(128, bounds=bounds, seed=7)
        self.assertEqual(len(pool), 128)
        q_min, q_max = transformed_domain(bounds)
        for node in pool:
            self.assertGreater(node.distance_m, abs(node.source_height_m - node.detector_height_m))
            q = geometry_transform(node)
            self.assertGreaterEqual(math.exp(q[1]), bounds.rho_min_m)
            self.assertTrue(all(low <= value <= high for value, low, high in zip(q, q_min, q_max)))
        a = geometry_transform(PhysicalNode(662.0, 5.0, 0.5, 2.0))
        b = geometry_transform(PhysicalNode(662.0, 5.0, 2.0, 0.5))
        self.assertEqual(a[:2], b[:2])
        self.assertNotEqual(a[2:], b[2:])

    def test_transformed_domain_is_explicit_not_sample_extrema(self):
        bounds = DesignBounds(rho_min_m=0.015)
        expected = transformed_domain(bounds)
        for seed in (1, 999):
            pool = generate_candidate_pool(32, bounds=bounds, seed=seed)
            sample_q = [geometry_transform(node) for node in pool]
            self.assertGreater(min(row[1] for row in sample_q), expected[0][1])
            self.assertLess(max(row[1] for row in sample_q), expected[1][1])

        air_min, air_max = transformed_domain(bounds, mode="air")
        air = generate_candidate_pool(32, bounds=bounds, seed=5, mode="air")
        self.assertTrue(all(node.source_height_m == node.detector_height_m == 0.0 for node in air))
        self.assertEqual(air_min[1], math.log(bounds.distance_m[0]))
        self.assertEqual(air_max[1], math.log(bounds.distance_m[1]))

    def test_disjoint_maximin_keeps_forced_anchor(self):
        candidates = generate_candidate_pool(80, seed=11)
        anchor = PhysicalNode(662.0, 5.0, 1.0, 1.0)
        selected = select_disjoint_maximin(
            candidates,
            sizes={"train": 12, "validation": 6, "final": 5},
            forced={"train": (anchor,)},
        )
        self.assertIn(anchor, selected.train)
        all_nodes = selected.train + selected.validation + selected.final
        self.assertEqual(len(all_nodes), len(set(all_nodes)))


class AxisTests(unittest.TestCase):
    def test_registered_landmark_and_nested_rebin_conserve_mass(self):
        fine = primary_registered_axis(662.0, n_per_block=4, n_mu=4)
        coarse = primary_registered_axis(662.0, n_per_block=2, n_mu=2)
        self.assertEqual(fine.physical_energy_edges_keV[4], 662.0 / (1.0 + 2.0 * 662.0 / 510.998950))
        masses = tuple(
            tuple(float(1 + ie * 4 + imu) for imu in range(4)) for ie in range(8)
        )
        rebinned = conservative_rebin_2d(
            masses,
            fine.physical_energy_edges_keV,
            fine.mu_edges,
            coarse.physical_energy_edges_keV,
            coarse.mu_edges,
        )
        self.assertTrue(math.isclose(sum(map(sum, masses)), sum(map(sum, rebinned)), rel_tol=1e-14))

    def test_production_384x64_to_192x32_conserves_mass_and_lines_stay_atoms(self):
        fine = primary_registered_axis(2614.0, n_per_block=192, n_mu=64)
        coarse = primary_registered_axis(2614.0, n_per_block=96, n_mu=32)
        masses = tuple(
            tuple(
                (1.0 + energy_bin) * (1.0 + 0.125 * mu_bin)
                for mu_bin in range(64)
            )
            for energy_bin in range(384)
        )
        rebinned = conservative_rebin_2d(
            masses,
            fine.physical_energy_edges_keV,
            fine.mu_edges,
            coarse.physical_energy_edges_keV,
            coarse.mu_edges,
        )
        self.assertEqual((len(rebinned), len(rebinned[0])), (192, 32))
        self.assertTrue(
            math.isclose(
                math.fsum(map(math.fsum, masses)),
                math.fsum(map(math.fsum, rebinned)),
                rel_tol=2.0e-15,
            )
        )

        # Lines carry only an exact energy label plus angular bin masses.  A
        # mu rebin must not turn either atom into an energy histogram.
        line_energies_keV = (2614.0, 510.998950)
        line_mu_masses = tuple(
            tuple(float(1 + line + mu_bin) for mu_bin in range(64))
            for line in range(2)
        )
        angular_overlap = overlap_matrix(fine.mu_edges, coarse.mu_edges)
        rebinned_lines = tuple(
            tuple(
                math.fsum(
                    source_mass * fraction
                    for source_mass, fraction in zip(line, target_row)
                )
                for target_row in angular_overlap
            )
            for line in line_mu_masses
        )
        self.assertEqual(line_energies_keV, (2614.0, 510.998950))
        for source, target in zip(line_mu_masses, rebinned_lines):
            self.assertTrue(
                math.isclose(math.fsum(source), math.fsum(target), rel_tol=1.0e-15)
            )

    def test_rebin_requires_target_coverage(self):
        with self.assertRaisesRegex(ValueError, "cover"):
            conservative_rebin_2d(
                ((1.0,),),
                (0.0, 1.0),
                (-1.0, 1.0),
                (0.1, 1.0),
                (-1.0, 1.0),
            )


if __name__ == "__main__":
    unittest.main()
