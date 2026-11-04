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

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest


NUMPY_HDF5_AVAILABLE = (
    importlib.util.find_spec("numpy") is not None
    and importlib.util.find_spec("h5py") is not None
)


@unittest.skipUnless(NUMPY_HDF5_AVAILABLE, "NumPy/HDF5 archive extras are not installed")
class DetectorFoldReportTests(unittest.TestCase):
    def _model(self, name, axis_kind, *, yield_value):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.fit import StructuredPositiveModel

        knots = [0.0] * 4 + [1.0] * 4
        features = 4**4
        yield_coefficients = np.full(features, math.log(yield_value))
        if axis_kind.startswith("line_"):
            archetypes = np.zeros((1, 4))
            archetypes[0, 0] = 1.0
            energy_knots = []
        else:
            archetypes = np.zeros((1, 4, 4))
            archetypes[0, 0, 0] = 1.0
            energy_knots = knots
        return StructuredPositiveModel(
            name,
            axis_kind,
            energy_knots,
            knots,
            [math.log(100.0), 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            1,
            1,
            0.0,
            yield_coefficients,
            np.zeros((features, 1)),
            archetypes,
        )

    def _write_response(
        self,
        path,
        detector_name,
        *,
        energy_edges=(10.0, 50.0, 100.0),
        mu_edges=(-1.0, 0.0, 1.0),
        energy_variation=False,
        angular_variation=False,
    ):
        import numpy as np

        n_energy = len(energy_edges) - 1
        n_mu = len(mu_edges) - 1
        continuum = np.zeros((n_energy, n_mu, 2), dtype=np.float64)
        for energy_index in range(n_energy):
            continuum[
                energy_index, :, min(energy_index, 1) if energy_variation else 0
            ] = 0.8
        line = np.zeros((2, n_mu, 2), dtype=np.float64)
        line[0, :, 0] = 0.8
        line[1, :, 1] = 0.8
        if angular_variation:
            continuum[:, 1::2] = continuum[:, 1::2, ::-1]
            line[:, 1::2] = line[:, 1::2, ::-1]
        geometry = {
            "schema": "CEEFIT-DETECTOR-GEOMETRY-v1",
            "detector_name": detector_name,
            "geometry_id": f"representative-{detector_name.lower()}-fixture-v1",
            "description": f"Representative synthetic {detector_name} validation detector",
            "active_material": "HPGe" if detector_name == "HPGe" else "NaI(Tl)",
            "dimensions_cm": {"length": 7.0, "radius": 3.0},
            "mu_definition": "cosine-to-detector-axis;mu=1=head-on-v1",
            "source_surface": "far-field-planar-angular-source",
        }
        provenance = {
            "schema": "CEEFIT-DETECTOR-RESPONSE-PROVENANCE-v1",
            "transport_code": "synthetic-test",
            "transport_version": "1",
            "physics_configuration": {"electromagnetic": "fixture-v1"},
            "histories_per_grid_point": 1000000,
            "rng": "counter-fixture-v1",
            "created_utc": "2026-08-08T00:00:00Z",
            "input_sha256": {"geometry": hashlib.sha256(detector_name.encode()).hexdigest()},
        }
        np.savez(
            path,
            schema=np.asarray("CEEFIT-DETECTOR-RESPONSE-v1"),
            detector_name=np.asarray(detector_name),
            response_semantics=np.asarray("probability-per-incident-cell-mass-v1"),
            line_interpolation=np.asarray("linear-exact-energy-v1"),
            geometry_json=np.asarray(
                json.dumps(geometry, sort_keys=True, separators=(",", ":"))
            ),
            provenance_json=np.asarray(
                json.dumps(provenance, sort_keys=True, separators=(",", ":"))
            ),
            generator_name=np.asarray("synthetic-response-generator"),
            generator_version=np.asarray("1.0"),
            generator_sha256=np.asarray(hashlib.sha256(b"synthetic-generator").hexdigest()),
            continuum_energy_edges_keV=np.asarray(energy_edges, dtype=np.float64),
            mu_edges=np.asarray(mu_edges, dtype=np.float64),
            spectrum_edges_keV=np.asarray((0.0, 50.0, 100.0), dtype=np.float64),
            continuum_response=continuum,
            line_energy_grid_keV=np.asarray((50.0, 600.0), dtype=np.float64),
            line_response=line,
        )

    def _fixture(self, root, *, truth_scale=1.0):
        from dataclasses import asdict
        import h5py

        from tools.full_spectrum.env_scatter_fit.archive import ingest_blocks
        from tools.full_spectrum.env_scatter_fit.block import read_fit_block
        from tools.full_spectrum.env_scatter_fit.campaign import PhysicsConfig
        from tools.full_spectrum.env_scatter_fit.design import PhysicalNode
        from tools.full_spectrum.env_scatter_fit.fit import RuntimePayload
        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256
        from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block

        root.mkdir(parents=True, exist_ok=True)
        block_paths = []
        nodes = [PhysicalNode(100.0, 1.0 + index * 0.001, 0.0, 0.0) for index in range(64)]
        for index, node in enumerate(nodes):
            direct = 1.0 / (4.0 * math.pi * (100.0 * node.distance_m) ** 2)
            path = root / f"final_{index:02d}.ceefit"
            path.write_bytes(
                fixture_block(
                    mu_bins=64,
                    role="final",
                    physical_node=node,
                    scale=direct * truth_scale,
                )
            )
            block_paths.append(path)
        first_block = read_fit_block(block_paths[0])
        canonical = first_block.manifest["canonical_config_json"]
        config_hash = first_block.manifest["config_sha256"]
        final_ids = [node.stable_id_hex(config_hash) for node in nodes]
        archive_path = root / "archive.h5"
        ingest_blocks(block_paths, archive_path)

        campaign = {
            "schema": "CEEFIT-CAMPAIGN-v1",
            "config": {
                "canonical_json": canonical,
                "sha256": config_hash,
                "mode": "air",
                "material": "soil",
                "physics": asdict(PhysicsConfig(energy_bins=2, mu_bins=64)),
            },
            "allocation": {"maximum_final_generations": 3},
            "campaign_state": {"final_generation": 0},
            "design": {
                "nodes": {
                    "train": [],
                    "validation": [],
                    "final": [
                        {"node_id": int(node_id, 16), "node_id_hex": node_id, "role": "final"}
                        for node_id in final_ids
                    ],
                }
            },
        }
        campaign_path = root / "campaign.json"
        campaign_path.write_text(json.dumps(campaign), encoding="utf-8")

        payload = RuntimePayload(
            {"primary_order1": self._model("primary_order1", "primary", yield_value=1.0)},
            {"elastic_E0": self._model("elastic_E0", "line_e0", yield_value=0.2)},
            {
                "config_sha256": config_hash,
                "material_name": "Air",
                "geometry_tier": "air",
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": {
                    "min": [100.0, 1.0, 0.0, 0.0],
                    "max": [100.0, nodes[-1].distance_m, 0.0, 0.0],
                },
                "q_domain": {
                    "min": [math.log(100.0), 0.0, 0.0, 0.0],
                    "max": [math.log(100.0) + 1.0, 1.0, 1.0, 1.0],
                },
                "brems_activation_threshold_keV": 150.0,
                "model_inventory": {
                    "continuum": [{"name": "primary_order1", "axis_kind": "primary"}],
                    "lines": [{"name": "elastic_E0", "axis_kind": "line_e0"}],
                },
            },
        )
        payload.model_content_sha256 = payload.compute_model_content_sha256()
        with h5py.File(archive_path, "r") as archive:
            snapshot_hash = archive_snapshot_sha256(archive)
            total_axis = int(archive["axes"].attrs["total_energy_axis"])
            response_energy_edges = tuple(
                float(value)
                for value in archive[
                    f"axes/energy/{total_axis}/physical_edges_keV"
                ][0]
            )
        payload.metadata["validation_metadata"] = {
            "model_content_sha256": payload.model_content_sha256,
            "archive_snapshot_sha256": snapshot_hash,
            "final": {
                "physics_passed": True,
                "nodes": [{"node_id": node_id} for node_id in final_ids],
            },
        }
        candidate_path = root / "payload.candidate"
        payload.save_binary(candidate_path)
        hpge = root / "hpge.npz"
        nai = root / "nai.npz"
        self._write_response(hpge, "HPGe", energy_edges=response_energy_edges)
        self._write_response(nai, "NaI", energy_edges=response_energy_edges)
        return candidate_path, archive_path, campaign_path, hpge, nai, final_ids

    def _mu_fixture(self, root, *, angular_variation=False, response_mu_edges=None):
        import h5py

        from tools.full_spectrum.env_scatter_fit.archive import ingest_blocks
        from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block

        root.mkdir(parents=True, exist_ok=True)
        block = root / "pilot.ceefit"
        block.write_bytes(fixture_block(mu_bins=64, role="pilot"))
        archive = root / "pilot.h5"
        ingest_blocks((block,), archive)
        with h5py.File(archive, "r") as stored:
            total_axis = int(stored["axes"].attrs["total_energy_axis"])
            energy_edges = tuple(
                float(value)
                for value in stored[
                    f"axes/energy/{total_axis}/physical_edges_keV"
                ][0]
            )
            archive_mu_edges = tuple(float(value) for value in stored["axes/mu_edges"][:])
        mu_edges = response_mu_edges or archive_mu_edges
        hpge = root / "hpge.npz"
        nai = root / "nai.npz"
        self._write_response(
            hpge,
            "HPGe",
            energy_edges=energy_edges,
            mu_edges=mu_edges,
            angular_variation=angular_variation,
        )
        self._write_response(
            nai,
            "NaI",
            energy_edges=energy_edges,
            mu_edges=mu_edges,
            angular_variation=angular_variation,
        )
        return archive, hpge, nai

    def test_real_candidate_and_synthetic_64_node_archive_pass_consumer(self):
        from tools.full_spectrum.env_scatter_fit.detector_fold import generate_detector_fold_report
        from tools.full_spectrum.env_scatter_fit.workflow import _external_detector_validation

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, archive, campaign, hpge, nai, final_ids = self._fixture(root)
            output = root / "folds.json"
            report = generate_detector_fold_report(
                candidate, archive, campaign, hpge, nai, output
            )
            self.assertTrue(report["passed"])
            self.assertEqual(report["final_node_ids"], final_ids)
            self.assertEqual(len(report["detectors"]["HPGe"]["nodes"]), 64)
            self.assertEqual(
                report["detectors"]["HPGe"]["response_sha256"],
                hashlib.sha256(hpge.read_bytes()).hexdigest(),
            )
            self.assertFalse(output.with_name(output.name + ".tmp").exists())
            consumed = _external_detector_validation(
                output,
                config_hash=report["config_sha256"],
                model_content_sha256=report["model_content_sha256"],
                archive_snapshot_hash=report["archive_snapshot_sha256"],
                final_node_ids=final_ids,
            )
            self.assertTrue(consumed["passed"])

    def test_metric_failure_writes_a_bound_failed_report(self):
        from tools.full_spectrum.env_scatter_fit.detector_fold import generate_detector_fold_report

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, archive, campaign, hpge, nai, _ = self._fixture(
                root, truth_scale=1.02
            )
            output = root / "failed.json"
            report = generate_detector_fold_report(
                candidate, archive, campaign, hpge, nai, output
            )
            self.assertFalse(report["passed"])
            self.assertGreater(
                report["detectors"]["HPGe"]["maxima"][
                    "continuum_total_relative_error"
                ],
                0.01,
            )
            self.assertEqual(json.loads(output.read_text())["passed"], False)

    def test_candidate_tamper_and_axis_coverage_fail_without_replacing_output(self):
        from tools.full_spectrum.env_scatter_fit.detector_fold import generate_detector_fold_report

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate, archive, campaign, hpge, nai, _ = self._fixture(root)
            output = root / "folds.json"
            output.write_text("sentinel\n", encoding="utf-8")
            raw = bytearray(candidate.read_bytes())
            raw[-1] ^= 1
            candidate.write_bytes(raw)
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                generate_detector_fold_report(
                    candidate, archive, campaign, hpge, nai, output
                )
            self.assertEqual(output.read_text(encoding="utf-8"), "sentinel\n")

            candidate, archive, campaign, hpge, nai, _ = self._fixture(root / "coverage")
            self._write_response(hpge, "HPGe", energy_edges=(20.0, 50.0, 100.0))
            with self.assertRaisesRegex(ValueError, "does not cover"):
                generate_detector_fold_report(
                    candidate, archive, campaign, hpge, nai, output
                )
            self.assertEqual(output.read_text(encoding="utf-8"), "sentinel\n")

    def test_exact_line_uses_explicit_energy_interpolation(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.detector_fold import _line_kernel, load_detector_response

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hpge.npz"
            self._write_response(path, "HPGe")
            response = load_detector_response(path, expected_detector_name="HPGe")
            kernel = _line_kernel(np, response, 325.0)
            np.testing.assert_allclose(
                kernel,
                0.5 * response.line_response[0] + 0.5 * response.line_response[1],
            )

    def test_response_dimensions_and_nonfinite_values_are_rejected(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.detector_fold import load_detector_response

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hpge.npz"
            self._write_response(path, "HPGe")
            with np.load(path, allow_pickle=False) as stored:
                record = {name: stored[name].copy() for name in stored.files}
            record["continuum_response"] = record["continuum_response"][:, :1, :]
            np.savez(path, **record)
            with self.assertRaisesRegex(ValueError, "shape"):
                load_detector_response(path, expected_detector_name="HPGe")

            self._write_response(path, "HPGe")
            with np.load(path, allow_pickle=False) as stored:
                record = {name: stored[name].copy() for name in stored.files}
            record["line_response"][0, 0, 0] = np.nan
            np.savez(path, **record)
            with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                load_detector_response(path, expected_detector_name="HPGe")

    def test_mu_resolution_generator_passes_consumer_and_binds_responses(self):
        from tools.full_spectrum.env_scatter_fit.detector_fold import (
            generate_mu_resolution_validation_report,
        )
        from tools.full_spectrum.env_scatter_fit.workflow import (
            _mu_detector_validation,
            mu_resolution_report,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, hpge, nai = self._mu_fixture(root)
            output = root / "mu_folds.json"
            report = generate_mu_resolution_validation_report(
                archive, hpge, nai, output
            )
            self.assertTrue(report["passed"])
            self.assertEqual(report["source_bins"], 64)
            self.assertEqual(report["candidate_bins"], 32)
            self.assertEqual(len(report["pilot_node_ids"]), 1)
            self.assertEqual(
                report["detectors"]["NaI"]["response_sha256"],
                hashlib.sha256(nai.read_bytes()).hexdigest(),
            )
            consumed = _mu_detector_validation(
                output,
                config_sha256=report["config_sha256"],
                archive_snapshot_hash=report["archive_snapshot_sha256"],
                pilot_node_ids=report["pilot_node_ids"],
                tolerance=0.002,
            )
            self.assertTrue(consumed["passed"])
            decision = mu_resolution_report(archive, detector_validation_path=output)
            self.assertTrue(decision["detector_validation"]["passed"])

    def test_mu_resolution_metrics_and_discretization_bound_fail_closed(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.detector_fold import (
            generate_mu_resolution_validation_report,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, hpge, nai = self._mu_fixture(root / "angular", angular_variation=True)
            angular_report = generate_mu_resolution_validation_report(
                archive, hpge, nai, root / "angular.json"
            )
            self.assertFalse(angular_report["passed"])
            self.assertGreater(
                angular_report["detectors"]["HPGe"]["maxima"][
                    "max_detector_folded_spectral_tv"
                ],
                0.002,
            )
            self.assertGreater(
                angular_report["detectors"]["HPGe"]["fold_discretization"][
                    "max_continuum_spectral_tv_bound"
                ],
                0.002,
            )

            fine_response_mu = tuple(np.linspace(-1.0, 1.0, 129))
            archive, hpge, nai = self._mu_fixture(
                root / "split",
                angular_variation=True,
                response_mu_edges=fine_response_mu,
            )
            split_report = generate_mu_resolution_validation_report(
                archive, hpge, nai, root / "split.json"
            )
            self.assertFalse(split_report["passed"])
            self.assertGreater(
                split_report["detectors"]["HPGe"]["fold_discretization"][
                    "max_continuum_spectral_tv_bound"
                ],
                0.002,
            )

    def test_response_provenance_is_mandatory_and_geometry_bound(self):
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.detector_fold import load_detector_response

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hpge.npz"
            self._write_response(path, "HPGe")
            with np.load(path, allow_pickle=False) as stored:
                record = {name: stored[name].copy() for name in stored.files}
            del record["generator_sha256"]
            np.savez(path, **record)
            with self.assertRaisesRegex(ValueError, "key inventory"):
                load_detector_response(path, expected_detector_name="HPGe")

            self._write_response(path, "HPGe")
            with np.load(path, allow_pickle=False) as stored:
                record = {name: stored[name].copy() for name in stored.files}
            geometry = json.loads(str(record["geometry_json"].item()))
            geometry["active_material"] = "NaI(Tl)"
            record["geometry_json"] = np.asarray(
                json.dumps(geometry, sort_keys=True, separators=(",", ":"))
            )
            np.savez(path, **record)
            with self.assertRaisesRegex(ValueError, "active material"):
                load_detector_response(path, expected_detector_name="HPGe")

    def test_cli_exposes_mu_fold_report(self):
        from tools.full_spectrum.env_scatter_fit.cli import build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "mu-fold-report",
                "--archive",
                "pilot.h5",
                "--hpge-response",
                "hpge.npz",
                "--nai-response",
                "nai.npz",
                "--output",
                "mu.json",
            ]
        )
        self.assertEqual(args.tolerance, 0.002)


if __name__ == "__main__":
    unittest.main()
