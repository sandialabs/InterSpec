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

"""The detector fitter's yield-node correction."""

import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import fit_detector_response as fit  # noqa: E402


def curved_model(spacing: float) -> fit.ComponentModel:
    """A yield field log y = sin(3 x0) + x1^2 on an uneven grid, one shape
    coefficient per node."""
    rng = np.random.default_rng(3)
    nodes = []
    for x0 in np.arange(0.0, 1.0 + 1e-9, spacing):
        for x1 in (0.0, 0.3, 0.35, 0.4, 1.0):  # crowded in the middle
            feature = np.array([x0 + 0.01 * rng.random(), x1, 0.5])
            nodes.append(fit.Node(feature=feature,
                                  log_yield=math.sin(3.0 * feature[0])
                                  + feature[1] ** 2,
                                  shape=np.array([1.0])))
    return fit.ComponentModel(
        yield_bandwidth=np.array([0.08, 0.15, 0.1]),
        shape_bandwidth=np.array([0.08, 0.15, 0.1]),
        basis_size=1, nodes=nodes)


def node_error(model: fit.ComponentModel, truth: list[float]) -> float:
    return max(abs(math.log(fit.interpolate(model, node.feature, False)) - y)
               for node, y in zip(model.nodes, truth))


class YieldNodeCorrection(unittest.TestCase):
    def test_zero_iterations_keep_the_training_values(self):
        model = curved_model(0.1)
        before = [node.log_yield for node in model.nodes]
        fit.interpolate_yield_nodes(model, 0)
        self.assertEqual([node.log_yield for node in model.nodes], before)

    def test_iterations_make_the_kernel_average_reproduce_the_nodes(self):
        truth = [node.log_yield for node in curved_model(0.1).nodes]
        errors = []
        for iterations in (0, 2, 8, 32):
            model = curved_model(0.1)
            fit.interpolate_yield_nodes(model, iterations)
            errors.append(node_error(model, truth))
        # The plain kernel average misses its own nodes (smoothing bias);
        # every step shrinks the miss.
        self.assertGreater(errors[0], 0.02)
        self.assertTrue(all(b < a for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], 0.1 * errors[0])

    def test_noise_is_not_amplified_at_the_nodes(self):
        # Pure noise on a flat field: the corrected field at the nodes never
        # strays further from zero than the noisy values themselves.
        model = curved_model(0.05)
        rng = np.random.default_rng(11)
        for node in model.nodes:
            node.log_yield = 0.01 * rng.standard_normal()
        noisy = np.array([node.log_yield for node in model.nodes])
        fit.interpolate_yield_nodes(model, 16)
        fitted = np.array([math.log(fit.interpolate(model, node.feature, False))
                           for node in model.nodes])
        self.assertLessEqual(np.sqrt(np.mean(fitted ** 2)),
                             np.sqrt(np.mean(noisy ** 2)) * (1 + 1e-9))


if __name__ == "__main__":
    unittest.main()
