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
from pathlib import Path
import tempfile
import unittest


NUMPY_AVAILABLE = importlib.util.find_spec("numpy") is not None


class OptionalDependencyTests(unittest.TestCase):
    def test_fit_dependency_is_actionable_or_available(self):
        from tools.full_spectrum.env_scatter_fit.fit import require_numpy

        if NUMPY_AVAILABLE:
            self.assertIsNotNone(require_numpy())
        else:
            with self.assertRaisesRegex(RuntimeError, "requirements-fit"):
                require_numpy()

    def test_archive_dependency_is_actionable_or_available(self):
        from tools.full_spectrum.env_scatter_fit.archive import require_hdf5

        if importlib.util.find_spec("h5py") is not None:
            self.assertEqual(len(require_hdf5()), 2)
        else:
            with self.assertRaisesRegex(RuntimeError, "requirements-archive"):
                require_hdf5()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy fit extra is not installed")
class FitTests(unittest.TestCase):
    def test_weak_component_validation_uses_absolute_total_effect(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import score_model

        class FixedModel:
            def predict_cell_masses(self, q, energy_edges, mu_edges):
                masses = np.asarray([[0.0042, 0.0063]])
                yields = masses.sum(axis=1)
                return yields, masses / yields[:, None], masses

        truth = np.asarray([[0.004, 0.006]])
        standalone_yield, _, _ = score_model(
            FixedModel(), np.zeros((1, 4)), truth, (0.0, 1.0), (-1.0, 1.0)
        )
        effective_yield, effective_shape, _ = score_model(
            FixedModel(),
            np.zeros((1, 4)),
            truth,
            (0.0, 1.0),
            (-1.0, 1.0),
            validation_reference_yields=np.asarray([1.0]),
            validation_material_mask=np.asarray([False]),
        )

        self.assertAlmostEqual(standalone_yield, 0.05)
        self.assertAlmostEqual(effective_yield, 0.0025)
        self.assertAlmostEqual(effective_shape, 0.0025)

    def test_elastic_line_validation_includes_direct_atom(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import score_model

        class FixedModel:
            def predict_cell_masses(self, q, energy_edges, mu_edges):
                masses = np.asarray([[0.0, 0.002]])
                yields = masses.sum(axis=1)
                return yields, masses / yields[:, None], masses

        truth = np.asarray([[0.0, 0.001]])
        baseline = np.asarray([[0.0, 1.0]])
        standalone_yield, _, _ = score_model(
            FixedModel(), np.zeros((1, 4)), truth, (0.0, 1.0), (-1.0, 1.0)
        )
        combined_yield, combined_shape, _ = score_model(
            FixedModel(),
            np.zeros((1, 4)),
            truth,
            (0.0, 1.0),
            (-1.0, 1.0),
            additive_baseline_masses=baseline,
        )

        self.assertAlmostEqual(standalone_yield, 1.0)
        self.assertLess(combined_yield, 0.001)
        self.assertAlmostEqual(combined_shape, 0.0)

    def test_uncertainty_weights_are_invariant_to_yield_units(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import _normalized_uncertainty_weights

        masses = np.asarray([[0.2, 0.8], [20.0, 80.0], [2000.0, 8000.0]])
        yields = np.sum(masses, axis=1)
        # Identical cell-relative uncertainty at yields spanning four orders
        # of magnitude must produce identical dimensionless fit weights.
        variances = (0.02 * masses) ** 2
        shape, sample = _normalized_uncertainty_weights(variances, yields)
        np.testing.assert_allclose(shape[0], shape[1], rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(shape[1], shape[2], rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(sample, np.ones(3), rtol=1e-14, atol=1e-14)

    def test_yield_weights_use_covariance_preserving_variance_when_supplied(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import _normalized_uncertainty_weights

        variances = np.ones((2, 2))
        yields = np.ones(2)
        _, fallback = _normalized_uncertainty_weights(variances, yields)
        _, covariance_preserving = _normalized_uncertainty_weights(
            variances, yields, np.asarray([2.0, 8.0])
        )
        np.testing.assert_allclose(fallback, np.ones(2))
        np.testing.assert_allclose(covariance_preserving, np.asarray([1.6, 0.4]))

    def test_replica_yield_variance_sums_cells_before_variance(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import _replica_yield_variance_of_mean

        # The two cells move together within a replica.  Summing marginal
        # variances would miss the positive covariance and be half as large.
        offsets = np.arange(32, dtype=float) - 15.5
        scores = np.stack((offsets, offsets), axis=1)[None, :, :]
        histories = np.ones((1, 32), dtype=float)
        observed = _replica_yield_variance_of_mean(np, scores, histories)[0]
        expected = np.var(2.0 * offsets, ddof=1) / 32.0
        self.assertAlmostEqual(observed, expected)

    def test_archive_replica_yield_variance_reads_bounded_row_batches(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import (
            _archive_component_has_nonzero,
            _replica_yield_variance_from_archive,
            _replica_yield_variance_of_mean,
        )

        class BoundedDataset:
            def __init__(self, values):
                self.values = values
                self.shape = values.shape
                self.ndim = values.ndim
                self.maximum_rows = 0

            def __getitem__(self, key):
                rows = key[0]
                if not isinstance(rows, slice) or rows.start is None or rows.stop is None:
                    raise AssertionError("replica archive read must use a bounded row slice")
                self.maximum_rows = max(self.maximum_rows, rows.stop - rows.start)
                return self.values[key]

        rng = np.random.default_rng(1234)
        scores = rng.uniform(0.0, 1.0, size=(19, 32, 13)).astype("f4")
        histories = np.full((19, 32), 8192.0)
        dataset = BoundedDataset(scores)
        expected = _replica_yield_variance_of_mean(
            np, scores[:, :, 3:11], histories
        )
        observed = _replica_yield_variance_from_archive(
            np, dataset, histories, 3, 11, row_batch=4
        )

        np.testing.assert_array_equal(observed, expected)
        self.assertLessEqual(dataset.maximum_rows, 4)

        zeros = BoundedDataset(np.zeros_like(scores))
        self.assertFalse(
            _archive_component_has_nonzero(
                np, zeros, 3, 11, row_batch=4
            )
        )
        self.assertLessEqual(zeros.maximum_rows, 4)
        scores[17, 4, 7] = 1.0
        nonzero = BoundedDataset(scores)
        self.assertTrue(
            _archive_component_has_nonzero(
                np, nonzero, 3, 11, row_batch=4
            )
        )
        self.assertLessEqual(nonzero.maximum_rows, 4)

    def test_archive_replica_reduction_is_bitwise_exact_at_production_width(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import (
            _replica_yield_variance_from_archive,
            _replica_yield_variance_of_mean,
        )

        class Dataset:
            def __init__(self, values):
                self.values = values
                self.shape = values.shape
                self.ndim = values.ndim

            def __getitem__(self, key):
                return self.values[key]

        rows = 3
        component_cells = 384 * 64
        rng = np.random.default_rng(4321)
        scores = rng.uniform(
            0.0,
            1.0,
            size=(rows, 32, 3 * component_cells),
        ).astype("f4")
        histories = np.full((rows, 32), 8192.0)
        dataset = Dataset(scores)
        for component in range(3):
            start = component * component_cells
            stop = start + component_cells
            expected = _replica_yield_variance_of_mean(
                np, scores[:, :, start:stop], histories
            )
            for row_batch in (1, 2, 8):
                with self.subTest(component=component, row_batch=row_batch):
                    observed = _replica_yield_variance_from_archive(
                        np,
                        dataset,
                        histories,
                        start,
                        stop,
                        row_batch=row_batch,
                    )
                    np.testing.assert_array_equal(observed, expected)

    def test_stopping_vector_uses_distinct_active_band_and_legendre_targets(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import _stopping_noise_metrics

        count = 100
        first = np.zeros(16)
        first[0] = 2.0  # mean 0.02: activates the >=1%-of-total band gate
        first[1] = 0.5  # mean 0.005: deliberately inactive
        outer = np.zeros((16, 16))

        def set_standard_error(index, standard_error):
            outer[index, index] = (
                first[index] ** 2 / count
                + standard_error**2 * count * (count - 1)
            )

        set_standard_error(0, 0.0009)
        set_standard_error(1, 0.5)
        set_standard_error(12, 0.0024)
        band, legendre = _stopping_noise_metrics(
            np, count, first, outer, total_mean=1.0
        )
        self.assertAlmostEqual(band, 0.0009)
        self.assertAlmostEqual(legendre, 0.0024)

    def test_component_thresholds_are_applied_with_explicit_q_domain(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.design import PhysicalNode, geometry_transform
        from tools.full_spectrum.env_scatter_fit.fit import RuntimePayload, select_rank_and_level

        nodes = [
            PhysicalNode(1100.0 + 40.0 * index, 5.0 + 0.1 * index, 1.0, 1.1)
            for index in range(6)
        ]
        q = np.asarray([geometry_transform(node) for node in nodes])
        q_domain = (q.min(axis=0), q.max(axis=0))
        masses = np.tile(np.asarray((0.2, 0.3, 0.1, 0.4)), (4, 1))
        for component, axis_kind, expected_threshold in (
            ("annihilation_order1", "annihilation", 2.0 * 510.998950),
            ("brems_order0", "brems", 150.0),
        ):
            selected, _ = select_rank_and_level(
                component,
                axis_kind,
                q[:4],
                masses,
                q[4:],
                masses[:2],
                registered_energy_edges=(0.0, 0.5, 1.0),
                mu_edges=(-1.0, 0.0, 1.0),
                ranks=(1,),
                levels=(1,),
                nmf_iterations=10,
                source_energy_threshold_keV=(
                    expected_threshold if axis_kind == "brems" else None
                ),
                q_domain=q_domain,
            )
            self.assertAlmostEqual(selected.model.pair_threshold_keV, expected_threshold)
            below = 100.0 if axis_kind == "brems" else 1000.0
            # Use the model directly at an in-domain q coordinate while checking
            # that the runtime support gate keys off physical source energy.
            payload = RuntimePayload(
                {component: selected.model},
                metadata={
                    "physical_domain": {
                        "min": [below, nodes[0].distance_m, 1.0, 1.1],
                        "max": [nodes[-1].energy_keV, nodes[-1].distance_m, 1.0, 1.1],
                    }
                },
            )
            # Query geometry must remain inside the fitted q box. Alter only the
            # physical E threshold check by supplying the low source energy and
            # a q box broad enough to cover it for this support test.
            selected.model.q_min[0] = min(selected.model.q_min[0], np.log(below))
            selected.model.q_span[0] = max(q_domain[1][0], np.log(below)) - selected.model.q_min[0]
            evaluated = payload.evaluate(
                below,
                nodes[0].distance_m,
                1.0,
                1.1,
                (0.0, 10.0, below),
                (-1.0, 0.0, 1.0),
                direct_fluence_cm2=1.0,
            )
            self.assertEqual(float(np.sum(evaluated["continuum_s_hat"])), 0.0)

    def test_cpp_binary_payload_fixture(self):
        import hashlib
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import RuntimePayload, StructuredPositiveModel

        knots = [0.0] * 4 + [1.0] * 4
        model = StructuredPositiveModel(
            "primary_order1",
            "primary",
            knots,
            knots,
            [0.0] * 4,
            [1.0] * 4,
            1,
            1,
            0.0,
            np.zeros((256, 1)),
            np.zeros((256, 1)),
            np.asarray([1.0] + [0.0] * 15).reshape(1, 4, 4),
        )
        self.assertEqual(model.parameter_count, 256 * 2 + 16)
        payload = RuntimePayload(
            {"primary_order1": model},
            metadata={
                "config_sha256": "0" * 64,
                "material_name": "Air",
                "geometry_tier": "air",
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": {"min": [1, 1, 0, 0], "max": [2, 2, 1, 1]},
                "q_domain": {"min": [0] * 4, "max": [1] * 4},
                "brems_activation_threshold_keV": 150.0,
                "model_inventory": {
                    "continuum": [
                        {"name": "primary_order1", "axis_kind": "primary"}
                    ],
                    "lines": [],
                },
                "validation_metadata": {},
            },
            accepted_for_publication=True,
        )
        encoded = payload.to_binary_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "golden.ceefitp"
            payload.save_binary(path)
            raw = path.read_bytes()
        self.assertEqual(raw, encoded)
        self.assertEqual(len(raw), 4875)
        self.assertEqual(
            hashlib.sha256(raw).hexdigest(),
            "7a947f72ddc89c54337afe8274f390b7443f1cc83a23b504bb32f4824b399808",
        )
        from tools.full_spectrum.env_scatter_fit.workflow import _bind_payload_size_report

        size_report = _bind_payload_size_report(payload)
        self.assertTrue(size_report["passed"])
        self.assertEqual(
            size_report["proposed_publication_bytes"],
            len(payload.to_binary_bytes()),
        )

    def test_positive_model_payload_round_trip(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.design import PhysicalNode, geometry_transform
        from tools.full_spectrum.env_scatter_fit.fit import RuntimePayload, StructuredPositiveModel

        q = []
        masses = []
        for index in range(12):
            node = PhysicalNode(100.0 + 40.0 * index, 5.0 + 0.2 * index, 1.0, 1.2)
            q.append(geometry_transform(node))
            shape = np.array((0.2 + 0.005 * index, 0.3, 0.1, 0.4 - 0.005 * index))
            shape /= shape.sum()
            masses.append(np.exp(0.1 * q[-1][0]) * shape)
        model = StructuredPositiveModel.fit(
            "primary_order1",
            "primary",
            np.asarray(q),
            np.asarray(masses),
            registered_energy_edges=(0.0, 0.5, 1.0),
            mu_edges=(-1.0, 0.0, 1.0),
            rank=2,
            level=1,
            nmf_iterations=40,
        )
        yields, shapes, prediction = model.predict_cell_masses(
            np.asarray(q), (0.0, 0.5, 1.0), (-1.0, 0.0, 1.0)
        )
        self.assertTrue(np.all(yields > 0.0))
        self.assertTrue(np.all(prediction >= 0.0))
        np.testing.assert_allclose(shapes.sum(axis=1), 1.0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "payload.json"
            metadata = {
                "config_sha256": "0" * 64,
                "material_name": "Soil",
                "geometry_tier": "field",
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": {
                    "min": [100.0, 5.0, 1.0, 1.2],
                    "max": [540.0, 7.2, 1.0, 1.2],
                },
                "q_domain": {
                    "min": model.q_min,
                    "max": [
                        low + span for low, span in zip(model.q_min, model.q_span)
                    ],
                },
                "brems_activation_threshold_keV": 150.0,
                "model_inventory": {
                    "continuum": [
                        {"name": "primary_order1", "axis_kind": "primary"}
                    ],
                    "lines": [],
                },
                "validation_metadata": {
                    "passed": True,
                    "fit_parameters": {
                        "input_parameters": {"fit": {"seed": 7}},
                        "output_parameters": {"observable": "fixture"},
                    },
                },
            }
            RuntimePayload(
                {"primary_order1": model},
                metadata=metadata,
                accepted_for_publication=True,
            ).save(path)
            restored = RuntimePayload.load(path)
            restored_prediction = restored.components["primary_order1"].predict_cell_masses(
                np.asarray(q), (0.0, 0.5, 1.0), (-1.0, 0.0, 1.0)
            )[2]
            np.testing.assert_allclose(restored_prediction, prediction, rtol=1e-12, atol=1e-14)
            binary = Path(temporary) / "payload.ceefitp"
            RuntimePayload(
                {"primary_order1": model},
                metadata=metadata,
                accepted_for_publication=True,
            ).save_binary(binary)
            loaded_binary = RuntimePayload.load_binary(binary)
            binary_model = loaded_binary.components["primary_order1"]
            self.assertEqual(
                RuntimePayload.load_binary(binary).metadata[
                    "brems_activation_threshold_keV"
                ],
                150.0,
            )
            self.assertEqual(
                loaded_binary.metadata["validation_metadata"]["fit_parameters"],
                metadata["validation_metadata"]["fit_parameters"],
            )
            binary_prediction = binary_model.predict_cell_masses(
                np.asarray(q), (0.0, 0.5, 1.0), (-1.0, 0.0, 1.0)
            )[2]
            np.testing.assert_allclose(binary_prediction, prediction, rtol=1e-12, atol=1e-14)

            candidate = RuntimePayload(
                {"primary_order1": model}, metadata=dict(metadata)
            )
            model_hash = candidate.compute_model_content_sha256()
            candidate.metadata["validation_metadata"] = {"passed": False, "round": 2}
            self.assertEqual(candidate.compute_model_content_sha256(), model_hash)
            candidate_json = Path(temporary) / "candidate.json"
            candidate.save(candidate_json)
            with self.assertRaisesRegex(ValueError, "candidate"):
                RuntimePayload.load(candidate_json)
            self.assertFalse(
                RuntimePayload.load(
                    candidate_json, allow_candidate=True
                ).accepted_for_publication
            )
            candidate_path = Path(temporary) / "candidate.ceefitp"
            candidate.save_binary(candidate_path)
            with self.assertRaisesRegex(ValueError, "candidate"):
                RuntimePayload.load_binary(candidate_path)
            restored_candidate = RuntimePayload.load_binary(
                candidate_path, allow_candidate=True
            )
            self.assertFalse(restored_candidate.accepted_for_publication)
            self.assertEqual(restored_candidate.model_content_sha256, model_hash)

            bad_inventory = RuntimePayload(
                {"primary_order1": model},
                metadata={
                    **metadata,
                    "model_inventory": {"continuum": [], "lines": []},
                },
            )
            with self.assertRaisesRegex(ValueError, "inventory"):
                bad_inventory.compute_model_content_sha256()

    def test_nontrivial_python_cpp_evaluator_golden_vector(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.design import PhysicalNode, geometry_transform
        from tools.full_spectrum.env_scatter_fit.fit import RuntimePayload, StructuredPositiveModel

        node = PhysicalNode(662.0, 8.0, 1.5, 0.8)
        q0 = np.asarray(geometry_transform(node))
        offsets = np.asarray((0.23, 0.37, 0.11, 0.19))
        q_span = np.asarray((1.1, 0.9, 0.6, 0.8))
        q_min = q0 - offsets
        indices = np.indices((4, 4, 4, 4)).reshape(4, -1).T
        i0, i1, i2, i3 = indices.T
        yield_coefficients = (
            0.4
            + 0.01 * (i0 - 1.5)
            + 0.02 * (i1 - 1.5)
            - 0.015 * (i2 - 1.5)
            + 0.005 * (i3 - 1.5)
        )[:, None]
        logits = np.column_stack(
            (
                -0.3 + 0.05 * i0 - 0.04 * i1 + 0.03 * i2 - 0.02 * i3,
                0.2 - 0.02 * i0 + 0.06 * i1 - 0.01 * i2 + 0.04 * i3,
            )
        )
        energy_knots = [0.0] * 4 + [0.4] + [1.0] * 4
        mu_knots = [0.0] * 4 + [0.3, 0.7] + [1.0] * 4
        energy_index, mu_index = np.indices((5, 6))
        archetypes = np.stack(
            (
                (energy_index + 1.0) * (mu_index + 2.0),
                (6.0 - energy_index) * (7.0 - mu_index)
                + 0.3 * energy_index * mu_index,
            )
        )
        archetypes /= archetypes.sum(axis=(1, 2), keepdims=True)
        continuum = StructuredPositiveModel(
            "brems_order0",
            "brems",
            energy_knots,
            mu_knots,
            q_min.tolist(),
            q_span.tolist(),
            2,
            1,
            0.0,
            yield_coefficients,
            logits,
            archetypes,
        )
        line_archetypes = np.stack(
            (np.arange(1.0, 7.0), np.arange(6.0, 0.0, -1.0))
        )
        line_archetypes /= line_archetypes.sum(axis=1, keepdims=True)
        line = StructuredPositiveModel(
            "elastic_E0",
            "line_e0",
            [],
            mu_knots,
            q_min.tolist(),
            q_span.tolist(),
            2,
            1,
            0.0,
            yield_coefficients - 1.1,
            logits,
            line_archetypes,
        )
        payload = RuntimePayload(
            {continuum.component: continuum},
            {line.component: line},
            {
                "physical_domain": {
                    "min": [662.0, 8.0, 1.5, 0.8],
                    "max": [662.0, 8.0, 1.5, 0.8],
                }
            },
        )
        result = payload.evaluate(
            662.0,
            8.0,
            1.5,
            0.8,
            (0.0, 10.0, 17.0, 83.0, 251.0, 500.0, 662.0, 900.0),
            (-1.0, -0.73, -0.1, 0.37, 1.0),
            direct_fluence_cm2=0.2,
        )
        expected_continuum = np.asarray(
            (
                0, 0, 0, 0,
                0.08332336753366282, 0.08478522302174381, 0.046526312153207856, 0.07838387164298673,
                0.12131648531579858, 0.13280593285509656, 0.0796843340401528, 0.15671859710585198,
                0.06857324135921093, 0.08270895357362763, 0.05479029268231385, 0.12335331946729351,
                0.04509909119202966, 0.05932341261672018, 0.0423212377411697, 0.10354784535377856,
                0.020488893946752323, 0.028309544066437688, 0.020960043780111244, 0.05322364672284197,
                0, 0, 0, 0,
            )
        )
        expected_line = np.asarray(
            (0.12409460886151082, 0.1348171647338879, 0.08019496002549564, 0.15562079951936106)
        )
        np.testing.assert_allclose(
            np.asarray(result["continuum_s_hat"]).reshape(-1),
            expected_continuum,
            rtol=0.0,
            atol=1.0e-10,
        )
        np.testing.assert_allclose(
            result["lines"][0]["s_hat_mu_bin_masses"],
            expected_line,
            rtol=0.0,
            atol=1.0e-10,
        )
        np.testing.assert_allclose(
            np.asarray(result["continuum"]).reshape(-1),
            0.2 * expected_continuum,
            rtol=0.0,
            atol=2.0e-11,
        )
        np.testing.assert_allclose(
            result["lines"][0]["mu_bin_masses"],
            0.2 * expected_line,
            rtol=0.0,
            atol=2.0e-11,
        )


if __name__ == "__main__":
    unittest.main()
