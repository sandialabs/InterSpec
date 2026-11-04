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
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tools.full_spectrum.env_scatter_fit.campaign import (
    _detector_fold_noise_criteria,
    generate_campaign_manifest,
)


def _response_evidence(detector_name, character):
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
        "input_sha256": {"geometry": character * 64},
    }

    def canonical_digest(value):
        return hashlib.sha256(
            json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

    return {
        "response_sha256": character * 64,
        "response_schema": "CEEFIT-DETECTOR-RESPONSE-v1",
        "response_semantics": "probability-per-incident-cell-mass-v1",
        "line_interpolation": "linear-exact-energy-v1",
        "geometry": geometry,
        "geometry_id": geometry["geometry_id"],
        "geometry_sha256": canonical_digest(geometry),
        "provenance": provenance,
        "provenance_sha256": canonical_digest(provenance),
        "generator": {
            "name": "synthetic-response-generator",
            "version": "1.0",
            "sha256": hashlib.sha256(detector_name.encode("ascii")).hexdigest(),
        },
        "dimensions": {
            "continuum_energy_bins": 4,
            "mu_bins": 4,
            "spectrum_bins": 4,
            "line_energy_grid_points": 2,
        },
        "fold_discretization": {
            "method": "archive-cell-response-kernel-diameter-v1",
            "gate": 0.01,
            "max_continuum_total_relative_bound": 0.0,
            "max_continuum_spectral_tv_bound": 0.0,
            "max_exact_line_relative_bound": 0.0,
        },
    }


class DetectorAllocationTests(unittest.TestCase):
    def _manifest_and_report(self, directory: Path):
        manifest = generate_campaign_manifest(
            candidate_count=256,
            train_count=32,
            validation_count=16,
            final_count=64,
        )
        final_ids = [
            item["node_id_hex"] for item in manifest["design"]["nodes"]["final"]
        ]
        nodes = [
            {
                "node_id": node_id,
                "continuum_total_relative_error": 0.0,
                "spectral_tv": 0.0,
                "exact_line_angular_fold_relative_error": 0.0,
                "replica_relative_se": 0.003,
                "replica_yield_relative_se": 0.003,
                "replica_spectral_tv_se": 0.002,
                "continuum_fold_discretization_total_relative_bound": 0.0,
                "continuum_fold_discretization_spectral_tv_bound": 0.0,
                "exact_line_fold_discretization_relative_bound": 0.0,
            }
            for node_id in final_ids
        ]
        report = {
            "schema": "CEEFIT-DETECTOR-FOLD-VALIDATION-v1",
            "config_sha256": manifest["config"]["sha256"],
            "archive_snapshot_sha256": "a" * 64,
            "model_content_sha256": "b" * 64,
            "final_node_ids": final_ids,
            "detectors": {
                name: {
                    **_response_evidence(name, character),
                    "nodes": nodes,
                }
                for name, character in (("HPGe", "c"), ("NaI", "d"))
            },
        }
        path = directory / "detector.json"
        path.write_text(json.dumps(report), encoding="utf-8")
        return manifest, report, path

    def test_archive_bound_detector_replica_noise_becomes_allocation_criteria(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest, _report, path = self._manifest_and_report(Path(temporary))
            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.publication.archive_snapshot_sha256",
                return_value="a" * 64,
            ):
                criteria, status = _detector_fold_noise_criteria(
                    path, root=object(), manifest=manifest
                )
        self.assertTrue(status["assessed"])
        self.assertFalse(status["all_noise_targets_met"])
        self.assertEqual(len(criteria), 64)
        for values in criteria.values():
            self.assertEqual([item.name for item in values], [
                "detector_fold:HPGe",
                "detector_fold:NaI",
            ])
            self.assertTrue(all(item.active and not item.passed for item in values))

    def test_detector_replica_summary_must_close_to_yield_and_shape_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            manifest, report, path = self._manifest_and_report(Path(temporary))
            report["detectors"]["HPGe"]["nodes"][0]["replica_relative_se"] = 0.001
            path.write_text(json.dumps(report), encoding="utf-8")
            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.publication.archive_snapshot_sha256",
                return_value="a" * 64,
            ):
                with self.assertRaisesRegex(ValueError, "does not close"):
                    _detector_fold_noise_criteria(path, root=object(), manifest=manifest)


if __name__ == "__main__":
    unittest.main()
