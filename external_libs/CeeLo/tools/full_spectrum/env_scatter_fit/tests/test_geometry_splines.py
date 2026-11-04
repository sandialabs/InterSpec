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

import io
import unittest
from unittest import mock

try:
    import numpy as np
except ImportError:  # pragma: no cover - dependency-free worker environment
    np = None


@unittest.skipIf(np is None, "NumPy is optional")
class HierarchicalGeometrySplineTests(unittest.TestCase):
    def test_levels_are_nested_dyadic_spaces(self):
        from tools.full_spectrum.env_scatter_fit.fit import (
            _geometry_basis_count,
            _geometry_feature_count,
            _geometry_knots,
            _tensor_spline_features,
        )

        previous = set()
        for level, expected_basis in enumerate((4, 5, 7, 11), start=1):
            knots = _geometry_knots(level)
            interior = {
                round(float(value), 15) for value in knots if 0.0 < value < 1.0
            }
            self.assertTrue(previous.issubset(interior))
            self.assertEqual(_geometry_basis_count(level), expected_basis)
            self.assertEqual(_geometry_feature_count(level), expected_basis**4)
            features = _tensor_spline_features(
                np.asarray(((0.13, 0.37, 0.61, 0.83),)), level
            )
            self.assertEqual(features.shape, (1, expected_basis**4))
            self.assertAlmostEqual(float(features.sum()), 1.0, places=13)
            previous = interior

    def test_curvature_has_affine_nullspace_and_penalizes_rough_fields(self):
        from tools.full_spectrum.env_scatter_fit.fit import (
            _apply_geometry_roughness,
            _geometry_basis_count,
            _geometry_knots,
        )

        for level in range(1, 5):
            knots = _geometry_knots(level)
            n_basis = _geometry_basis_count(level)
            greville = np.asarray(
                [np.mean(knots[index + 1 : index + 4]) for index in range(n_basis)]
            )
            grid = np.meshgrid(greville, greville, greville, greville, indexing="ij")
            affine = 0.7 + sum((axis + 1.0) * values for axis, values in enumerate(grid))
            affine_penalty = _apply_geometry_roughness(affine.reshape(-1), level)
            checkerboard = np.indices((n_basis,) * 4).sum(axis=0) % 2
            rough = 2.0 * checkerboard.astype(float) - 1.0
            rough_penalty = _apply_geometry_roughness(rough.reshape(-1), level)
            self.assertLess(
                float(np.linalg.norm(affine_penalty)),
                2.0e-10 * max(1.0, float(np.linalg.norm(affine))),
            )
            self.assertGreater(float(rough.reshape(-1) @ rough_penalty), 1.0)

    def test_fine_level_roughness_fit_is_finite_and_stable(self):
        from tools.full_spectrum.env_scatter_fit.fit import (
            _fit_geometry_fields,
            _tensor_spline_features,
        )

        rng = np.random.default_rng(0xCEE10)
        q = rng.random((96, 4))
        targets = np.column_stack(
            (
                0.3 + q @ np.asarray((0.2, -0.1, 0.05, 0.15)),
                np.sin(2.0 * np.pi * q[:, 0]),
            )
        )
        coefficients = _fit_geometry_fields(
            q, 4, targets, np.ones(q.shape[0]), 1.0e-5
        )
        prediction = _tensor_spline_features(q, 4) @ coefficients
        self.assertTrue(np.all(np.isfinite(coefficients)))
        self.assertLess(float(np.sqrt(np.mean((prediction[:, 0] - targets[:, 0]) ** 2))), 1e-6)
        self.assertLess(float(np.sqrt(np.mean((prediction[:, 1] - targets[:, 1]) ** 2))), 0.04)

    def test_affine_feature_evaluation_and_serialized_dimensions(self):
        from tools.full_spectrum.env_scatter_fit.fit import (
            RuntimePayload,
            StructuredPositiveModel,
            _geometry_basis_count,
            _geometry_feature_count,
            _geometry_knots,
            _tensor_spline_features,
        )

        point = np.asarray(((0.13, 0.37, 0.61, 0.83),))
        slopes = np.asarray((0.2, -0.3, 0.4, 0.1))
        expected = 0.7 + float(point[0] @ slopes)
        archetype = np.asarray([1.0] + [0.0] * 15).reshape(1, 4, 4)
        knots = [0.0] * 4 + [1.0] * 4
        for level in range(1, 5):
            geometry_knots = _geometry_knots(level)
            n_basis = _geometry_basis_count(level)
            greville = np.asarray(
                [
                    np.mean(geometry_knots[index + 1 : index + 4])
                    for index in range(n_basis)
                ]
            )
            grid = np.meshgrid(greville, greville, greville, greville, indexing="ij")
            coefficients = 0.7 + sum(
                slopes[axis] * values for axis, values in enumerate(grid)
            )
            evaluated = float(
                (
                    _tensor_spline_features(point, level)
                    @ coefficients.reshape(-1)
                ).item()
            )
            self.assertAlmostEqual(evaluated, expected, places=12)
            model = StructuredPositiveModel(
                "primary_order1",
                "primary",
                knots,
                knots,
                [0.0] * 4,
                [1.0] * 4,
                1,
                level,
                1.0e-5,
                coefficients.reshape(-1, 1),
                np.zeros((_geometry_feature_count(level), 1)),
                archetype,
            )
            self.assertEqual(
                model.parameter_count, 2 * _geometry_feature_count(level) + 16
            )
            RuntimePayload._put_model(io.BytesIO(), model)

        # A nonnested level-three coefficient count (6^4) must not be accepted
        # as a level-three payload; such a candidate fails loudly rather than
        # being evaluated with a different spline space.
        invalid = model
        invalid.level = 3
        invalid.log_yield_coefficients = np.zeros((6**4, 1))
        invalid.logit_coefficients = np.zeros((6**4, 1))
        with self.assertRaisesRegex(ValueError, "coefficient count|shape"):
            RuntimePayload._put_model(io.BytesIO(), invalid)

    def test_rank_shape_decomposition_is_reused_across_levels(self):
        import tools.full_spectrum.env_scatter_fit.fit as fit_module

        rng = np.random.default_rng(42)
        train_q = rng.random((8, 4))
        validation_q = rng.random((2, 4))
        masses = rng.uniform(0.1, 1.0, size=(10, 4))
        original = fit_module._fit_shape_decomposition
        with mock.patch.object(
            fit_module, "_fit_shape_decomposition", wraps=original
        ) as decomposition:
            _, candidates = fit_module.select_rank_and_level(
                "primary_order1",
                "primary",
                train_q,
                masses[:8],
                validation_q,
                masses[8:],
                registered_energy_edges=(0.0, 0.5, 1.0),
                mu_edges=(-1.0, 0.0, 1.0),
                ranks=(1,),
                levels=(1, 2),
                nmf_iterations=10,
                q_domain=(
                    np.minimum(train_q.min(axis=0), validation_q.min(axis=0)),
                    np.maximum(train_q.max(axis=0), validation_q.max(axis=0)),
                ),
            )
        self.assertEqual(decomposition.call_count, 1)
        self.assertEqual(len(candidates), 2)
        np.testing.assert_array_equal(
            candidates[0].model.archetype_coefficients,
            candidates[1].model.archetype_coefficients,
        )


if __name__ == "__main__":
    unittest.main()
