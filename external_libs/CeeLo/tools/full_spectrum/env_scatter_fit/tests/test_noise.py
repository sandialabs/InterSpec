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

import importlib.util
import math
import unittest

from tools.full_spectrum.env_scatter_fit.noise import (
    EXACT_LINE_MATERIAL_FRACTION,
    direct_line_map,
    exact_line_materiality,
    rebin_joint_masses,
    replica_shape_tv_noise,
    sentinel_shape,
)


NUMPY_AVAILABLE = importlib.util.find_spec("numpy") is not None


@unittest.skipUnless(NUMPY_AVAILABLE, "numpy is not installed")
class JointShapeNoiseTests(unittest.TestCase):
    def test_exact_line_materiality_uses_three_sigma_upper_bound(self):
        material, upper = exact_line_materiality(0.004, 0.0019, 1.0)
        self.assertFalse(material)
        self.assertLess(upper, EXACT_LINE_MATERIAL_FRACTION)

        material, upper = exact_line_materiality(0.004, 0.002, 1.0)
        self.assertTrue(material)
        self.assertAlmostEqual(upper, EXACT_LINE_MATERIAL_FRACTION)

    def test_invalid_line_uncertainty_fails_conservatively(self):
        material, upper = exact_line_materiality(0.0, math.inf, 1.0)
        self.assertTrue(material)
        self.assertTrue(math.isinf(upper))

    def test_direct_atom_occupies_head_on_mu_bin(self):
        import numpy as np

        np.testing.assert_array_equal(direct_line_map(np, 3.0, 4), [0, 0, 0, 3])

    def test_rebin_conserves_every_replica_mass(self):
        import numpy as np

        values = np.arange(2 * 8 * 8, dtype=float).reshape(2, -1)
        rebinned = rebin_joint_masses(
            np, values, native_shape=(8, 8), target_shape=(2, 4)
        )
        self.assertEqual(rebinned.shape, (2, 8))
        np.testing.assert_allclose(rebinned.sum(axis=1), values.sum(axis=1))

    def test_subcell_jitter_does_not_force_more_histories(self):
        import numpy as np

        # Every replica has the exact same 2x2 parent masses, but alternates
        # between the four native cells inside each parent.  That is genuine
        # native-bin sampling noise, not uncertainty in the stopping-scale
        # joint shape.
        replicas = np.zeros((32, 16), dtype=float)
        for replica in range(32):
            for energy_parent in range(2):
                for mu_parent in range(2):
                    energy = 2 * energy_parent + replica % 2
                    mu = 2 * mu_parent + (replica // 2) % 2
                    replicas[replica, energy * 4 + mu] = 0.25
        truth = replicas.mean(axis=0)
        native = replica_shape_tv_noise(
            np,
            replicas,
            np.ones(32),
            truth,
            native_shape=(4, 4),
            target_shape=(4, 4),
        )
        stopping = replica_shape_tv_noise(
            np,
            replicas,
            np.ones(32),
            truth,
            native_shape=(4, 4),
            target_shape=(2, 2),
        )
        self.assertGreater(native, 0.05)
        self.assertAlmostEqual(stopping, 0.0)

    def test_tiny_cell_uses_absolute_not_relative_noise(self):
        import numpy as np

        replicas = np.tile(np.asarray([1.0, 0.0]), (32, 1))
        replicas[::2] = [1.0 - 1.0e-4, 1.0e-4]
        truth = replicas.mean(axis=0)
        noise = replica_shape_tv_noise(
            np,
            replicas,
            np.ones(32),
            truth,
            native_shape=(1, 2),
            target_shape=(1, 2),
        )
        # The second cell has 100% relative fluctuations, but only 5e-5 of
        # the total mass; its aggregate TV noise remains correspondingly tiny.
        self.assertLess(noise, 1.0e-5)
        self.assertTrue(math.isfinite(noise))

    def test_production_native_grid_selects_12_by_8(self):
        self.assertEqual(sentinel_shape(384, 64), (12, 8))
        self.assertEqual(sentinel_shape(384, 32), (12, 8))
        self.assertEqual(sentinel_shape(1, 64), (1, 8))


if __name__ == "__main__":
    unittest.main()
