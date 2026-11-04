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
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from tools.full_spectrum.env_scatter_fit.campaign import (
    MODEL_SELECTION_HISTORY_SCHEMA,
    MODEL_SELECTION_REPORT_SCHEMA,
    canonical_fit_config_json,
)
from tools.full_spectrum.env_scatter_fit.publication import (
    assemble_publication_evidence,
    validate_publication_evidence,
)
from tools.full_spectrum.env_scatter_fit.noise import contribution_aware_model_selection_policy
from tools.full_spectrum.env_scatter_fit.workflow import (
    PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES,
    _external_detector_validation,
    _mu_detector_validation,
    _write_fitted_payload,
    promote_candidate,
)


REGISTERED_CANONICAL_CONFIG_JSON = canonical_fit_config_json("field", "soil")


def digest(label):
    payload = REGISTERED_CANONICAL_CONFIG_JSON if label == "config" else label
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def detector_response_evidence(detector_name, gate):
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
        "input_sha256": {"geometry": digest(f"{detector_name}-geometry")},
    }

    def canonical_digest(value):
        return hashlib.sha256(
            json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

    return {
        "response_path": f"raw/{detector_name}.npz",
        "response_sha256": digest(detector_name),
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
            "sha256": digest(f"{detector_name}-generator"),
        },
        "dimensions": {
            "continuum_energy_bins": 4,
            "mu_bins": 4,
            "spectrum_bins": 4,
            "line_energy_grid_points": 2,
        },
        "fold_discretization": {
            "method": "archive-cell-response-kernel-diameter-v1",
            "gate": gate,
            "max_continuum_total_relative_bound": 0.0,
            "max_continuum_spectral_tv_bound": 0.0,
            "max_exact_line_relative_bound": 0.0,
        },
    }


def geant4_report(config, archive):
    radius = 4.0
    volume = (4.0 / 3.0) * math.pi * radius**3
    arms = [
        {
            "path": "raw/g4-arm-1.csv",
            "sha256": digest("g4-arm-1"),
            "events": 1000,
            "e0_keV": 100.0,
            "sphere_radius_cm": radius,
            "gdml_scoring_radius_cm": radius,
            "energy_bins": 48,
            "mu_bins": 32,
            "cell_count": 1536,
            "energy_min_keV": 10.0,
            "energy_max_keV": 3000.0,
            "mu_min": -1.0,
            "mu_max": 1.0,
            "observable": "scattered-track-length-flux-per-primary-v1",
            "physics_list": "G4EmLivermorePhysics",
            "world_half_extent_m": 60.0,
            "primary_sampling": "gps-isotropic-4pi",
            "rng_engine": "MixMaxRng",
            "rng_seed": 101,
            "production_cut_mm": 0.7,
            "gdml_sha256": digest("g4-gdml"),
            "macro_sha256": digest("g4-macro-1"),
            "inputs": {
                "gdml": {"path": "raw/env.gdml", "sha256": digest("g4-gdml")},
                "macro": {"path": "raw/arm-1.mac", "sha256": digest("g4-macro-1")},
            },
            "geometry_tier": "field",
            "material": "soil",
            "distance_cm": 1000.0,
            "source_height_cm": 100.0,
            "detector_height_cm": 100.0,
            "scattered_total_moment_schema": "per-event-track-length-v1",
            "scattered_mean": 0.99,
            "scattered_standard_error": 0.01,
        },
        {
            "path": "raw/g4-arm-2.csv",
            "sha256": digest("g4-arm-2"),
            "events": 3000,
            "e0_keV": 100.0,
            "sphere_radius_cm": radius,
            "gdml_scoring_radius_cm": radius,
            "energy_bins": 48,
            "mu_bins": 32,
            "cell_count": 1536,
            "energy_min_keV": 10.0,
            "energy_max_keV": 3000.0,
            "mu_min": -1.0,
            "mu_max": 1.0,
            "observable": "scattered-track-length-flux-per-primary-v1",
            "physics_list": "G4EmLivermorePhysics",
            "world_half_extent_m": 60.0,
            "primary_sampling": "gps-isotropic-4pi",
            "rng_engine": "MixMaxRng",
            "rng_seed": 202,
            "production_cut_mm": 0.7,
            "gdml_sha256": digest("g4-gdml"),
            "macro_sha256": digest("g4-macro-2"),
            "inputs": {
                "gdml": {"path": "raw/env.gdml", "sha256": digest("g4-gdml")},
                "macro": {"path": "raw/arm-2.mac", "sha256": digest("g4-macro-2")},
            },
            "geometry_tier": "field",
            "material": "soil",
            "distance_cm": 1000.0,
            "source_height_cm": 100.0,
            "detector_height_cm": 100.0,
            "scattered_total_moment_schema": "per-event-track-length-v1",
            "scattered_mean": 1.0033333333333334,
            "scattered_standard_error": 0.01,
        },
    ]
    for item in arms:
        n = item["events"]
        first = item["scattered_mean"] * n * volume
        second = first * first / n + (
            item["scattered_standard_error"] * volume
        ) ** 2 * n * (n - 1)
        item["scattered_total_N"] = n
        item["scattered_total_S1_cm"] = first
        item["scattered_total_S2_cm2"] = second
        item["scattered_mean"] = first / (n * volume)
        item["scattered_standard_error"] = math.sqrt(
            max(0.0, second - first * first / n) / (n * (n - 1))
        ) / volume
    events = sum(item["events"] for item in arms)
    total_s1 = math.fsum(item["scattered_total_S1_cm"] for item in arms)
    total_s2 = math.fsum(item["scattered_total_S2_cm2"] for item in arms)
    reference_mean = total_s1 / (events * volume)
    reference_se = math.sqrt(
        max(0.0, total_s2 - total_s1 * total_s1 / events)
        / (events * (events - 1))
    ) / volume
    weights = [float(item["events"]) for item in arms]
    weight_sum = math.fsum(weights)
    squared_weight_sum = math.fsum(weight * weight for weight in weights)
    between_variance = math.fsum(
        weight * (item["scattered_mean"] - reference_mean) ** 2
        for weight, item in zip(weights, arms)
    ) / (weight_sum - squared_weight_sum / weight_sum)
    between_arm_se = math.sqrt(
        between_variance / (weight_sum * weight_sum / squared_weight_sum)
    )
    combined_reference_se = max(reference_se, between_arm_se)
    archive_mean = 1.002
    archive_se = 0.005
    difference = archive_mean - reference_mean
    combined_se = math.hypot(archive_se, combined_reference_se)
    standardized = difference / combined_se
    relative = abs(difference) / reference_mean
    anchor = {
        "anchor_id": "field-soil-100keV",
        "node_id": "00000000000000a1",
        "coordinates": [100.0, 1000.0, 100.0, 100.0, 6000.0],
        "archive_role": "diagnostic",
        "archive_status": "complete",
        "world_half_extent_m": 60.0,
        "reference": {
            "events": events,
            "arm_count": len(arms),
            "e0_keV": 100.0,
            "sphere_radius_cm": radius,
            "gdml_scoring_radius_cm": radius,
            "energy_bins": 48,
            "mu_bins": 32,
            "energy_min_keV": 10.0,
            "energy_max_keV": 3000.0,
            "mu_min": -1.0,
            "mu_max": 1.0,
            "observable": "scattered-track-length-flux-per-primary-v1",
            "physics_list": "G4EmLivermorePhysics",
            "world_half_extent_m": 60.0,
            "primary_sampling": "gps-isotropic-4pi",
            "rng_engine": "MixMaxRng",
            "production_cut_mm": 0.7,
            "geometry_tier": "field",
            "material": "soil",
            "distance_cm": 1000.0,
            "source_height_cm": 100.0,
            "detector_height_cm": 100.0,
            "mean": reference_mean,
            "scattered_total_moment_schema": "per-event-track-length-v1",
            "scattered_total_N": events,
            "scattered_total_S1_cm": total_s1,
            "scattered_total_S2_cm2": total_s2,
            "per_history_standard_error": reference_se,
            "between_arm_standard_error": between_arm_se,
            "standard_error": combined_reference_se,
            "files": arms,
        },
        "archive": {
            "histories": 1048576,
            "mean": archive_mean,
            "standard_error": archive_se,
        },
        "difference": difference,
        "absolute_relative_difference": relative,
        "combined_standard_error": combined_se,
        "standardized_difference": standardized,
        "passed": True,
    }
    return {
        "schema": "CEEFIT-GEANT4-ANCHOR-AGREEMENT-v1",
        "config_sha256": config,
        "archive_snapshot_sha256": archive,
        "passed": True,
        "anchor_count": 1,
        "within_combined_uncertainty": True,
        "physics_list": "G4EmLivermorePhysics",
        "comparison_sigma": 2.0,
        "minimum_field_world_half_extent_m": 60.0,
        "manifest_path": "g4_manifest.json",
        "manifest_sha256": digest("g4-manifest"),
        "reference_file_count": len(arms),
        "reference_event_count": events,
        "max_abs_standardized_difference": abs(standardized),
        "max_absolute_relative_difference": relative,
        "anchors": [anchor],
    }


def detailed_reports(config, archive, model):
    grid_energies = (59.5, 122.0, 239.0, 478.0, 955.0, 1836.0, 2614.0)
    heldout_energies = tuple(
        math.sqrt(left * right)
        for left, right in zip(grid_energies, grid_energies[1:])
    )
    midpoint_metrics = [
        {
            "energy_keV": energy,
            "integrated_relative_error": 0.0,
            "joint_tv": 0.0,
            "energy_marginal_tv": 0.0,
            "angular_marginal_tv": 0.0,
            "continuum_fraction_of_scattered": 1.0,
        }
        for energy in heldout_energies
    ]
    stopping_nodes = [
        {
            "role": "grid" if index < 7 else "heldout",
            "energy_keV": (*grid_energies, *heldout_energies)[index],
            "histories": 12_000_000,
            "total_continuum_s_hat": 1.0,
            "joint_cell_masses": [1.0] + [0.0] * 95,
            "total_relative_se": 0.002,
            "max_significant_absolute_se_fraction": 0.0008,
            "total_target_met": True,
            "absolute_target_met": True,
            "relative_band_target_met": True,
        }
        for index in range(13)
    ]
    return {
        "registered_e0_gate": {
            "schema": "CEEFIT-REGISTERED-E0-GATE-v1",
            "config_sha256": config,
            "passed": True,
            "predefined_index_bands": True,
            "statistically_empty_bands_included": True,
            "shape_coverage_fraction": 1.0,
            "max_integrated_relative_error": 0.0,
            "max_joint_tv": 0.0,
            "max_energy_marginal_tv": 0.0,
            "max_angular_marginal_tv": 0.0,
            "all_node_stopping_targets_met": True,
            "best_scheme": "lin(s) vs lin(x)",
            "observable": "total_continuum; exact lines remain atoms",
            "gate_configuration": {
                "geometry_tier": "field",
                "material": "Soil",
                "distance_cm": 500.0,
                "source_height_cm": 100.0,
                "detector_height_cm": 100.0,
                "scoring_radius_cm": 4.0,
                "energy_cut_keV": 10.0,
            },
            "grid": {
                "native_registered_energy_bins": 384,
                "native_mu_bins": 64,
                "registered_energy_bands": 12,
                "registered_energy_edges": [index / 12 for index in range(13)],
                "mu_bands": 8,
                "mu_edges": [-1.0 + index / 4 for index in range(9)],
            },
            "stopping_targets": {
                "total_relative_se": 0.0025,
                "absolute_se_fraction_of_total": 0.001,
                "significant_fraction_of_total": 0.01,
                "batch_histories": 1_000_000,
                "max_histories_per_node": 128_000_000,
            },
            "acceptance_thresholds": {
                "max_integrated_relative_error": 0.01,
                "max_joint_tv": 0.01,
                "max_energy_marginal_tv": 0.005,
                "max_angular_marginal_tv": 0.005,
            },
            "midpoints": midpoint_metrics,
            "nodes": stopping_nodes,
            "generator": {
                "schema": "CEEFIT-REGISTERED-E0-GENERATOR-v1",
                "argv": ["env_registered_e0_gate", "--axis", "e0"],
                "input_parameters": {
                    "config_sha256": config,
                    "canonical_configuration_json": REGISTERED_CANONICAL_CONFIG_JSON,
                    "canonical_configuration": json.loads(
                        REGISTERED_CANONICAL_CONFIG_JSON
                    ),
                    "gate_configuration": {
                        "geometry_tier": "field",
                        "material": "Soil",
                        "distance_cm": 500.0,
                        "source_height_cm": 100.0,
                        "detector_height_cm": 100.0,
                        "scoring_radius_cm": 4.0,
                        "energy_cut_keV": 10.0,
                    },
                    "batch_histories": 1_000_000,
                },
                "output_parameters": {"observable": "total continuum"},
            },
        },
        "uncertainty_calibration": {
            "schema": "CEEFIT-UNCERTAINTY-CALIBRATION-v1",
            "config_sha256": config,
            "archive_snapshot_sha256": archive,
            "passed": True,
            "repeated_null_runs": 8,
            "s2_passed": True,
            "coarse_covariance_passed": True,
            "replica_covariance_passed": True,
            "family_coverage_passed": True,
            "criteria": {
                "minimum_independent_comparisons": 8,
                "required_moment_and_replica_tensors": [
                    "continuum_components",
                    "continuum_total",
                    "line_components",
                    "line_total",
                ],
            },
            "replication": {
                "sufficient": True,
                "independent": True,
                "comparison_count": 8,
                "blocks_used": 16,
            },
            "raw_total_s2": {"passed": True},
            "cell_s2": {
                "passed": True,
                "family_inventory_complete": True,
                "required_tensors": [
                    "continuum_components",
                    "continuum_total",
                    "line_components",
                    "line_total",
                ],
                "tensors": {
                    name: {"passed": True, "segments": [{"passed": True}]}
                    for name in (
                        "continuum_components",
                        "continuum_total",
                        "line_components",
                        "line_total",
                    )
                },
            },
            "stopping_outer_moments": {
                "passed": True,
                "dimensions": [
                    {"index": index, "passed": True} for index in range(16)
                ],
            },
            "replica_covariance": {
                "passed": True,
                "family_inventory_complete": True,
                "families": {
                    name: {"passed": True}
                    for name in (
                        "continuum_components",
                        "continuum_total",
                        "line_components",
                        "line_total",
                    )
                },
            },
        },
        "geant4_anchor_agreement": geant4_report(config, archive),
        "python_cpp_parity": {
            "schema": "CEEFIT-PYTHON-CPP-PARITY-v1",
            "config_sha256": config,
            "model_content_sha256": model,
            "passed": True,
            "candidate_payload_sha256": digest("candidate"),
            "artifact_paths": {"candidate_payload": "raw/candidate.bin"},
            "protocol": "CEEFIT-PARITY-REQUEST-v1",
            "request_sha256": digest("request"),
            "cpp_result_sha256": digest("result"),
            "tolerance": 1.0e-9,
            "query_count": 15,
            "coverage": {
                "interior": 1,
                "domain_face": 8,
                "brems_threshold": 3,
                "pair_threshold": 3,
            },
            "compared_value_count": 1234,
            "max_absolute_difference": 8.0e-10,
        },
    }


def write_detailed_reports(root, config, archive, model):
    Path(root).mkdir(parents=True, exist_ok=True)
    raw = Path(root) / "raw"
    raw.mkdir(exist_ok=True)
    (Path(root) / "g4_manifest.json").write_bytes(b"g4-manifest")
    for name in ("g4-arm-1", "g4-arm-2", "g4-gdml", "g4-macro-1", "g4-macro-2"):
        suffix = ".csv" if name.startswith("g4-arm") else ".gdml" if name == "g4-gdml" else ".mac"
        filename = (
            "env.gdml"
            if name == "g4-gdml"
            else f"arm-{name[-1]}.mac"
            if name.startswith("g4-macro")
            else f"{name}.csv"
        )
        (raw / filename).write_bytes(name.encode("ascii"))
    (raw / "candidate.bin").write_bytes(b"candidate")
    reports = detailed_reports(config, archive, model)
    paths = {}
    for name, record in reports.items():
        path = Path(root) / f"{name}.json"
        path.write_text(json.dumps(record, sort_keys=True) + "\n", encoding="utf-8")
        paths[name] = path
    return paths


class PublicationEvidenceTests(unittest.TestCase):
    def test_missing_evidence_is_a_candidate_not_an_exception(self):
        result = validate_publication_evidence(
            None,
            config_sha256=digest("config"),
            archive_snapshot_sha256=digest("archive"),
            model_content_sha256=digest("model"),
            final_node_ids=("0000000000000001",),
        )
        self.assertFalse(result["supplied"])
        self.assertFalse(result["passed"])

    def test_registered_e0_summary_is_recomputed_and_raw_generator_is_required(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = write_detailed_reports(root / "reports", config, archive, model)
            registered = json.loads(
                sources["registered_e0_gate"].read_text(encoding="utf-8")
            )
            registered["max_joint_tv"] = 0.001
            sources["registered_e0_gate"].write_text(
                json.dumps(registered, sort_keys=True) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "independently recomputed raw nodes"):
                assemble_publication_evidence(
                    root / "evidence.json",
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=["0000000000000001"],
                    report_paths=sources,
                )

            for mutate, message in (
                (
                    lambda generator: generator.pop("argv"),
                    "exact generator arguments",
                ),
                (
                    lambda generator: generator.__setitem__("schema", "unbound"),
                    "integrated generator record",
                ),
            ):
                sources = write_detailed_reports(
                    root / "reports", config, archive, model
                )
                registered = json.loads(
                    sources["registered_e0_gate"].read_text(encoding="utf-8")
                )
                mutate(registered["generator"])
                sources["registered_e0_gate"].write_text(
                    json.dumps(registered, sort_keys=True) + "\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(ValueError, message):
                    assemble_publication_evidence(
                        root / "evidence.json",
                        config_sha256=config,
                        archive_snapshot_sha256=archive,
                        model_content_sha256=model,
                        final_node_ids=["0000000000000001"],
                        report_paths=sources,
                    )

    def test_all_evidence_is_identity_bound_and_metric_checked(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        final_ids = ["0000000000000002", "0000000000000001"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "evidence.json"
            sources = write_detailed_reports(root / "reports", config, archive, model)
            assembled = assemble_publication_evidence(
                path,
                config_sha256=config,
                archive_snapshot_sha256=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
                report_paths=sources,
            )
            self.assertTrue(assembled["passed"])
            result = validate_publication_evidence(
                path,
                config_sha256=config,
                archive_snapshot_sha256=archive,
                model_content_sha256=model,
                final_node_ids=tuple(reversed(final_ids)),
            )
            self.assertTrue(result["passed"])

            aggregate = json.loads(path.read_text(encoding="utf-8"))
            aggregate["checks"]["python_cpp_parity"]["max_absolute_difference"] = 0.0
            path.write_text(json.dumps(aggregate, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "differs from its source"):
                validate_publication_evidence(
                    path,
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                )
            assemble_publication_evidence(
                path,
                config_sha256=config,
                archive_snapshot_sha256=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
                report_paths=sources,
            )

            parity = json.loads(sources["python_cpp_parity"].read_text(encoding="utf-8"))
            parity["max_absolute_difference"] = 1.1e-9
            sources["python_cpp_parity"].write_text(
                json.dumps(parity, sort_keys=True) + "\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "SHA-256 differs"):
                validate_publication_evidence(
                    path,
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                )
            result = assemble_publication_evidence(
                path,
                config_sha256=config,
                archive_snapshot_sha256=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
                report_paths=sources,
            )
            self.assertFalse(result["passed"])

            sources["python_cpp_parity"].unlink()
            with self.assertRaisesRegex(ValueError, "missing"):
                validate_publication_evidence(
                    path,
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                )

            with self.assertRaisesRegex(ValueError, "does not match"):
                validate_publication_evidence(
                    path,
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=digest("different-model"),
                    final_node_ids=final_ids,
                )

    def test_geant4_source_requires_recomputable_detailed_anchor_evidence(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        final_ids = ["0000000000000001"]

        def remove_anchors(record):
            record.pop("anchors")

        mutations = (
            ("summary only", remove_anchors, "lacks 'anchors'"),
            (
                "anchor count",
                lambda record: record.__setitem__("anchor_count", 2),
                "anchor_count differs",
            ),
            (
                "physics list",
                lambda record: record.__setitem__(
                    "physics_list", "G4EmStandardPhysics_option4"
                ),
                "must use G4EmLivermorePhysics",
            ),
            (
                "comparison sigma",
                lambda record: record.__setitem__("comparison_sigma", 3.0),
                "comparison_sigma must equal 2",
            ),
            (
                "world extent",
                lambda record: record["anchors"][0].__setitem__(
                    "world_half_extent_m", 59.0
                ),
                "inadequate field-world extent",
            ),
            (
                "source hash",
                lambda record: record["anchors"][0]["reference"]["files"][0].__setitem__(
                    "sha256", "not-a-sha256"
                ),
                "lacks a source SHA-256",
            ),
            (
                "duplicate source content",
                lambda record: record["anchors"][0]["reference"]["files"][1].__setitem__(
                    "sha256",
                    record["anchors"][0]["reference"]["files"][0]["sha256"],
                ),
                "identical source content more than once",
            ),
            (
                "merged reference",
                lambda record: record["anchors"][0]["reference"].__setitem__(
                    "mean", 1.1
                ),
                "reference.mean is inconsistent",
            ),
            (
                "missing exact total moment",
                lambda record: record["anchors"][0]["reference"]["files"][0].pop(
                    "scattered_total_S2_cm2"
                ),
                "scattered_total_S2_cm2 must be .*finite",
            ),
            (
                "biased sampling",
                lambda record: record["anchors"][0]["reference"]["files"][0].__setitem__(
                    "primary_sampling", "cone-biased"
                ),
                "uses biased sampling",
            ),
            (
                "standardized residual",
                lambda record: record["anchors"][0].__setitem__(
                    "standardized_difference", 0.0
                ),
                "standardized_difference is inconsistent",
            ),
            (
                "anchor pass",
                lambda record: record["anchors"][0].__setitem__("passed", False),
                r"anchors\[0\].passed is inconsistent",
            ),
            (
                "top maximum",
                lambda record: record.__setitem__(
                    "max_abs_standardized_difference", 0.0
                ),
                "max_abs_standardized_difference is inconsistent",
            ),
            (
                "top pass",
                lambda record: record.__setitem__("passed", False),
                "passed flag is inconsistent",
            ),
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "evidence.json"
            sources = write_detailed_reports(root / "reports", config, archive, model)
            for label, mutate, message in mutations:
                with self.subTest(label=label):
                    record = geant4_report(config, archive)
                    mutate(record)
                    sources["geant4_anchor_agreement"].write_text(
                        json.dumps(record, sort_keys=True) + "\n", encoding="utf-8"
                    )
                    with self.assertRaisesRegex(ValueError, message):
                        assemble_publication_evidence(
                            output,
                            config_sha256=config,
                            archive_snapshot_sha256=archive,
                            model_content_sha256=model,
                            final_node_ids=final_ids,
                            report_paths=sources,
                        )

    def test_detector_folds_are_bound_to_candidate_archive_response_and_nodes(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        final_ids = ["0000000000000001", "0000000000000002"]
        node_metrics = [
            {
                "node_id": node_id,
                "continuum_total_relative_error": 0.004,
                "spectral_tv": 0.006,
                "exact_line_angular_fold_relative_error": 0.007,
                "replica_relative_se": 0.002,
                "continuum_fold_discretization_total_relative_bound": 0.0,
                "continuum_fold_discretization_spectral_tv_bound": 0.0,
                "exact_line_fold_discretization_relative_bound": 0.0,
            }
            for node_id in final_ids
        ]
        record = {
            "schema": "CEEFIT-DETECTOR-FOLD-VALIDATION-v1",
            "config_sha256": config,
            "archive_snapshot_sha256": archive,
            "model_content_sha256": model,
            "final_node_ids": final_ids,
            "detectors": {
                name: {
                    **detector_response_evidence(name, 0.01),
                    "nodes": node_metrics,
                }
                for name in ("HPGe", "NaI")
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "raw").mkdir()
            for name in ("HPGe", "NaI"):
                (root / "raw" / f"{name}.npz").write_bytes(name.encode("ascii"))
            path = root / "detectors.json"
            path.write_text(json.dumps(record), encoding="utf-8")
            result = _external_detector_validation(
                path,
                config_hash=config,
                archive_snapshot_hash=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
            )
            self.assertTrue(result["passed"])
            self.assertTrue(result["physics_passed"])
            self.assertTrue(result["noise_passed"])
            self.assertEqual(result["noise_failing_node_ids"], [])
            self.assertEqual(result["detectors"]["HPGe"]["node_count"], 2)

            node_metrics[1]["replica_relative_se"] = 0.003
            path.write_text(json.dumps(record), encoding="utf-8")
            result = _external_detector_validation(
                path,
                config_hash=config,
                archive_snapshot_hash=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
            )
            self.assertFalse(result["passed"])
            self.assertTrue(result["physics_passed"])
            self.assertFalse(result["noise_passed"])
            self.assertEqual(result["noise_failing_node_ids"], [final_ids[1]])

            node_metrics[1]["replica_relative_se"] = 0.002
            node_metrics[1]["spectral_tv"] = 0.02
            path.write_text(json.dumps(record), encoding="utf-8")
            result = _external_detector_validation(
                path,
                config_hash=config,
                archive_snapshot_hash=archive,
                model_content_sha256=model,
                final_node_ids=final_ids,
            )
            self.assertFalse(result["passed"])
            self.assertFalse(result["physics_passed"])
            self.assertTrue(result["noise_passed"])

            node_metrics[1]["spectral_tv"] = 0.006
            record["archive_snapshot_sha256"] = digest("stale archive")
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "different fit archive"):
                _external_detector_validation(
                    path,
                    config_hash=config,
                    archive_snapshot_hash=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                )

            record["archive_snapshot_sha256"] = archive
            generator = record["detectors"]["HPGe"].pop("generator")
            path.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "generator"):
                _external_detector_validation(
                    path,
                    config_hash=config,
                    archive_snapshot_hash=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                )
            record["detectors"]["HPGe"]["generator"] = generator

    def test_uncertainty_source_requires_eight_comparisons_and_family_closure(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        final_ids = ["0000000000000001"]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = write_detailed_reports(root / "reports", config, archive, model)
            calibration_path = sources["uncertainty_calibration"]
            calibration = json.loads(calibration_path.read_text(encoding="utf-8"))
            calibration["repeated_null_runs"] = 7
            calibration["replication"].update(
                sufficient=False, comparison_count=7, blocks_used=14
            )
            calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "passed flag is inconsistent"):
                assemble_publication_evidence(
                    root / "evidence.json",
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                    report_paths=sources,
                )

            calibration = detailed_reports(config, archive, model)[
                "uncertainty_calibration"
            ]
            calibration["cell_s2"]["tensors"]["line_total"]["passed"] = False
            calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cell-S2 passed flag is inconsistent"):
                assemble_publication_evidence(
                    root / "evidence.json",
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                    report_paths=sources,
                )

            calibration = detailed_reports(config, archive, model)[
                "uncertainty_calibration"
            ]
            calibration["failures"] = [
                {"gate": "family_coverage", "reason": "injected contradiction"}
            ]
            calibration_path.write_text(json.dumps(calibration), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cannot contain failure records"):
                assemble_publication_evidence(
                    root / "evidence.json",
                    config_sha256=config,
                    archive_snapshot_sha256=archive,
                    model_content_sha256=model,
                    final_node_ids=final_ids,
                    report_paths=sources,
                )

    def test_mu_fold_decision_requires_every_bound_detector_sentinel(self):
        config, archive = digest("config"), digest("archive")
        pilot_ids = ["0000000000000010", "0000000000000020"]
        node_metrics = [
            {
                "node_id": node_id,
                "detector_folded_total_relative_change": 0.001,
                "detector_folded_spectral_tv": 0.0015,
                "exact_line_angular_fold_relative_change": 0.0019,
                "continuum_fold_discretization_total_relative_bound": 0.0,
                "continuum_fold_discretization_spectral_tv_bound": 0.0,
                "exact_line_fold_discretization_relative_bound": 0.0,
            }
            for node_id in pilot_ids
        ]
        record = {
            "schema": "CEEFIT-MU-RESOLUTION-VALIDATION-v1",
            "config_sha256": config,
            "archive_snapshot_sha256": archive,
            "source_bins": 64,
            "candidate_bins": 32,
            "tolerance": 0.002,
            "pilot_node_ids": pilot_ids,
            "detectors": {
                name: {
                    **detector_response_evidence(name, 0.002),
                    "nodes": node_metrics,
                }
                for name in ("HPGe", "NaI")
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "raw").mkdir()
            for name in ("HPGe", "NaI"):
                (root / "raw" / f"{name}.npz").write_bytes(name.encode("ascii"))
            path = root / "mu.json"
            path.write_text(json.dumps(record), encoding="utf-8")
            result = _mu_detector_validation(
                path,
                config_sha256=config,
                archive_snapshot_hash=archive,
                pilot_node_ids=pilot_ids,
                tolerance=0.002,
            )
            self.assertTrue(result["passed"])

            record["detectors"]["NaI"]["nodes"][0][
                "detector_folded_spectral_tv"
            ] = 0.0021
            path.write_text(json.dumps(record), encoding="utf-8")
            result = _mu_detector_validation(
                path,
                config_sha256=config,
                archive_snapshot_hash=archive,
                pilot_node_ids=pilot_ids,
                tolerance=0.002,
            )
            self.assertFalse(result["passed"])

    def test_fit_acceptance_reports_oversize_without_rejecting_valid_model(self):
        class SizedEncoding:
            def __init__(self, size):
                self.size = size

            def __len__(self):
                return self.size

        class Payload:
            accepted_for_publication = False

            def __init__(self):
                self.metadata = {}
                self.saved = None

            def to_binary_bytes(self):
                metadata_bytes = len(
                    json.dumps(
                        self.metadata, sort_keys=True, separators=(",", ":")
                    ).encode("utf-8")
                )
                return SizedEncoding(PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES + metadata_bytes)

            def save_binary(self, path):
                self.saved = Path(path)

        payload = Payload()
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "model.ceefitp"
            report = _write_fitted_payload(
                payload,
                output,
                {"policy": {"action": "publish_payload"}},
                acceptance=True,
            )
        self.assertTrue(report["acceptance_passed"])
        self.assertTrue(payload.accepted_for_publication)
        self.assertEqual(payload.saved.name, "model.ceefitp")
        self.assertTrue(report["runtime_payload_size"]["passed"])
        self.assertFalse(report["runtime_payload_size"]["within_advisory_limit"])
        self.assertEqual(
            report["runtime_payload_size"]["exact_bound_candidate_bytes"],
            len(payload.to_binary_bytes()),
        )

    def test_candidate_promotion_reuses_locked_model_without_refitting(self):
        config, archive, model = digest("config"), digest("archive"), digest("model")
        final_ids = [f"{value:016x}" for value in range(1, 65)]

        class FakePayload:
            accepted_for_publication = False
            model_content_sha256 = model

            def __init__(self):
                self.metadata = {
                    "config_sha256": config,
                    "validation_metadata": {
                        "model_selection_policy": (
                            contribution_aware_model_selection_policy()
                        ),
                        "model_content_sha256": model,
                        "archive_snapshot_sha256": archive,
                        "selection": {"all": {"validation_passed": True}},
                        "final": {
                            "physics_passed": True,
                            "passed": False,
                            "component_acceptance_semantics": (
                                "continuum-origin-order-plus-exact-atom-sentinel3-v1"
                            ),
                            "transport_status": {"passed": True},
                            "gates": {"detector_fold_validation": False},
                            "nodes": [{"node_id": value} for value in final_ids],
                        },
                    },
                }
                self.saved = None

            def compute_model_content_sha256(self):
                return model

            def to_binary_bytes(self):
                return b"candidate-payload-fixture"

            def save_binary(self, path):
                self.saved = Path(path)

        payload = FakePayload()

        class Dataset:
            def __init__(self, values):
                self.values = values

            def __getitem__(self, key):
                return self.values[key]

        class TextDataset(Dataset):
            def asstr(self):
                return self

        class Root:
            attrs = {"schema": "CEEFIT-HDF5-v1"}

            def __getitem__(self, key):
                if key == "manifest":
                    return SimpleNamespace(attrs={"config_sha256": config})
                if key == "nodes/id":
                    return Dataset(list(range(1, 65)))
                if key == "nodes/status":
                    return TextDataset(["complete"] * 64)
                if key == "nodes/transport_failures":
                    return Dataset([0] * 64)
                if key == "nodes/capped":
                    return Dataset([0] * 64)
                if key == "nodes/histories":
                    return Dataset([4 * 262144] * 64)
                raise KeyError(key)

        class OpenArchive:
            def __enter__(self):
                return Root()

            def __exit__(self, *unused):
                return False

        campaign = {
            "config": {"sha256": config},
            "design": {
                "nodes": {
                    "validation": [{"node_id_hex": "f" * 16}],
                    "final": [{"node_id_hex": value} for value in final_ids],
                }
            },
            "campaign_state": {
                "round": 2,
                "final_generation": 0,
                "stopped": True,
                "consecutive_validation_passes": 2,
                "validation_history": [
                    {
                        "schema": MODEL_SELECTION_HISTORY_SCHEMA,
                        "kind": "bound_model_selection_report",
                        "passed": True,
                        "report_schema": MODEL_SELECTION_REPORT_SCHEMA,
                        "report_canonical_sha256": digest(f"report-{round_number}"),
                        "config_sha256": config,
                        "archive_snapshot_sha256": archive,
                        "model_content_sha256": model,
                        "eligible_validation_node_ids": ["f" * 16],
                        "campaign_round": round_number,
                        "final_generation": 0,
                        "campaign_round_identity_sha256": digest(
                            f"round-identity-{round_number}"
                        ),
                    }
                    for round_number in (1, 2)
                ],
            },
        }
        with tempfile.TemporaryDirectory() as temporary, patch(
            "tools.full_spectrum.env_scatter_fit.workflow.RuntimePayload.load_binary",
            return_value=payload,
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.load_manifest", return_value=campaign
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_hdf5",
            return_value=(SimpleNamespace(File=lambda *unused: OpenArchive()), object()),
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.archive_snapshot_sha256", return_value=archive
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_production_configuration",
            return_value={"passed": True, "schema": "CEEFIT-PRODUCTION-CONFIGURATION-v1"},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow._external_detector_validation",
            return_value={"supplied": True, "passed": True},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.validate_publication_evidence",
            return_value={"supplied": True, "passed": True},
        ):
            root = Path(temporary)
            (root / "model.candidate").write_bytes(payload.to_binary_bytes())
            report = promote_candidate(
                root / "model.candidate",
                root / "archive.h5",
                root / "campaign.json",
                root / "detector.json",
                root / "evidence.json",
                root / "model.ceefitp",
            )
        self.assertTrue(report["acceptance_passed"])
        self.assertTrue(payload.accepted_for_publication)
        self.assertEqual(payload.saved.name, "model.ceefitp")
        self.assertEqual(report["policy"]["action"], "publish_payload_without_refit")
        self.assertTrue(report["runtime_payload_size"]["passed"])
        self.assertEqual(
            report["runtime_payload_size"]["exact_input_candidate_bytes"],
            len(payload.to_binary_bytes()),
        )

        class FailedRoot(Root):
            def __getitem__(self, key):
                if key == "nodes/status":
                    return TextDataset(["capped", *(["complete"] * 63)])
                if key == "nodes/capped":
                    return Dataset([1, *([0] * 63)])
                return super().__getitem__(key)

        class OpenFailedArchive:
            def __enter__(self):
                return FailedRoot()

            def __exit__(self, *unused):
                return False

        failed_payload = FakePayload()
        with tempfile.TemporaryDirectory() as temporary, patch(
            "tools.full_spectrum.env_scatter_fit.workflow.RuntimePayload.load_binary",
            return_value=failed_payload,
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.load_manifest", return_value=campaign
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_hdf5",
            return_value=(
                SimpleNamespace(File=lambda *unused: OpenFailedArchive()),
                object(),
            ),
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.archive_snapshot_sha256",
            return_value=archive,
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_production_configuration",
            return_value={"passed": True, "schema": "CEEFIT-PRODUCTION-CONFIGURATION-v1"},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow._external_detector_validation",
            side_effect=AssertionError("invalid transport must not reach detector validation"),
        ):
            root = Path(temporary)
            failed = promote_candidate(
                root / "model.candidate",
                root / "archive.h5",
                root / "campaign.json",
                root / "detector.json",
                root / "evidence.json",
                root / "model.ceefitp",
            )
        self.assertFalse(failed["acceptance_passed"])
        self.assertEqual(
            failed["policy"]["action"],
            "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
        )
        self.assertEqual(failed["transport_status"]["failed_node_ids"], [final_ids[0]])

        class SizedEncoding:
            def __len__(self):
                return PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES + 1

        class OversizedPayload(FakePayload):
            def __init__(self):
                super().__init__()
                self.saved = None

            def to_binary_bytes(self):
                return SizedEncoding()

            def save_binary(self, path):
                self.saved = Path(path)

        oversized_payload = OversizedPayload()
        with tempfile.TemporaryDirectory() as temporary, patch(
            "tools.full_spectrum.env_scatter_fit.workflow.RuntimePayload.load_binary",
            return_value=oversized_payload,
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.load_manifest", return_value=campaign
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_hdf5",
            return_value=(SimpleNamespace(File=lambda *unused: OpenArchive()), object()),
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.archive_snapshot_sha256", return_value=archive
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_production_configuration",
            return_value={"passed": True, "schema": "CEEFIT-PRODUCTION-CONFIGURATION-v1"},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow._external_detector_validation",
            return_value={"supplied": True, "passed": True},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.validate_publication_evidence",
            return_value={"supplied": True, "passed": True},
        ):
            root = Path(temporary)
            candidate_path = root / "oversized.candidate"
            candidate_path.write_bytes(b"small-input-candidate")
            oversized = promote_candidate(
                candidate_path,
                root / "archive.h5",
                root / "campaign.json",
                root / "detector.json",
                root / "evidence.json",
                root / "oversized.ceefitp",
            )
        self.assertTrue(oversized["acceptance_passed"])
        self.assertEqual(oversized_payload.saved.name, "oversized.ceefitp")
        self.assertTrue(oversized["runtime_payload_size"]["passed"])
        self.assertFalse(
            oversized["runtime_payload_size"]["within_advisory_limit"]
        )
        self.assertEqual(
            oversized["runtime_payload_size"]["proposed_publication_bytes"],
            PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES + 1,
        )


if __name__ == "__main__":
    unittest.main()
