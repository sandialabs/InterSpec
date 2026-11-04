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

import json
import hashlib
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tools.full_spectrum.env_scatter_fit.allocation import (
    BLOCK_HISTORIES,
    LOGICAL_REPLICAS,
    MIN_INITIAL_BLOCKS,
    AcquisitionCandidate,
    ContributionMoment,
    NodeStatistics,
    RawMoment,
    SentinelMoment,
    ValidationRoundState,
    allocate_actions,
    assess_node,
    compare_adaptive_batches,
    continuation_block_count,
    continuation_loss_reduction_components,
    effective_history_cap,
    history_cap_blocks,
    predict_candidate_utilities,
)
from tools.full_spectrum.env_scatter_fit.campaign import (
    MODEL_SELECTION_REPORT_SCHEMA,
    PhysicsConfig,
    _load_adaptive_utility_report,
    advance_campaign_round,
    campaign_round_identity,
    canonical_fit_config_json,
    config_sha256,
    default_forced_anchors,
    generate_campaign_manifest,
    load_manifest,
    replace_final_set_from_reserve,
    replica_seed,
    save_manifest,
)
from tools.full_spectrum.env_scatter_fit.design import DesignBounds
from tools.full_spectrum.env_scatter_fit.design import PhysicalNode
from tools.full_spectrum.env_scatter_fit.noise import contribution_aware_model_selection_policy


def moment(mean, standard_deviation, histories):
    return RawMoment(
        histories * mean,
        (histories - 1) * standard_deviation**2 + histories * mean**2,
        histories,
    )


def replica_moment_with_se(mean, standard_error):
    """Raw moments of 32 replica means with a requested SE of their mean."""

    return moment(mean, standard_error * math.sqrt(LOGICAL_REPLICAS), LOGICAL_REPLICAS)


def model_selection_report(
    manifest, *, passed=True, archive_hash="a" * 64, model_hash="b" * 64
):
    identity = campaign_round_identity(manifest)
    return {
        "schema": MODEL_SELECTION_REPORT_SCHEMA,
        "model_selection_policy": contribution_aware_model_selection_policy(),
        "config_sha256": manifest["config"]["sha256"],
        "archive_snapshot_sha256": archive_hash,
        "model_content_sha256": model_hash,
        "campaign_round": identity["campaign_round"],
        "final_generation": identity["final_generation"],
        "campaign_round_identity_sha256": identity["identity_sha256"],
        "eligible_validation_node_ids": identity["validation_node_ids"],
        "passed": passed,
        "selection": {
            "transport_status": {
                "kind": "transport_status",
                "validation_passed": passed,
            },
            "mc_noise": {"kind": "mc_noise", "validation_passed": passed},
        },
    }


def final_failure_report(manifest):
    report = model_selection_report(manifest)
    final_ids = [
        item["node_id_hex"] for item in manifest["design"]["nodes"]["final"]
    ]
    report["final"] = {
        "evaluated": True,
        "passed": False,
        "physics_passed": False,
        "transport_status": {"passed": True},
        "gates": {"mc_noise_eligible": True},
        "detector_validation": {"supplied": False},
        "nodes": [{"node_id": node_id} for node_id in final_ids],
    }
    report["policy"] = {
        "action": (
            "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve"
        ),
        "retired_node_ids": final_ids,
    }
    return report


class AllocationTests(unittest.TestCase):
    def test_builtin_acquisition_predictor_targets_cost_weighted_validation_loss(self):
        result = predict_candidate_utilities(
            training_q=((0.0, 0.0, 0.0, 0.0),),
            training_seconds_per_block=(2.0,),
            validation_q=(
                (0.8, 0.0, 0.0, 0.0),
                (0.0, 0.8, 0.0, 0.0),
            ),
            validation_losses=(10.0, 0.1),
            candidates={
                "near_failure": (0.75, 0.0, 0.0, 0.0),
                "near_good": (0.0, 0.75, 0.0, 0.0),
                "already_trained": (0.0, 0.0, 0.0, 0.0),
            },
            q_min=(0.0, 0.0, 0.0, 0.0),
            q_max=(1.0, 1.0, 1.0, 1.0),
        )
        self.assertGreater(
            result["near_failure"]["utility"], result["near_good"]["utility"]
        )
        self.assertEqual(result["already_trained"]["utility"], 0.0)
        self.assertEqual(result["near_failure"]["predicted_cpu_seconds"], 8.0)

    def test_raw_moment_formula_and_merge(self):
        value = moment(2.0, 4.0, 100)
        self.assertAlmostEqual(value.mean, 2.0)
        self.assertAlmostEqual(value.standard_error, 0.4)
        merged = value.merge(value)
        self.assertEqual(merged.histories, 200)
        self.assertAlmostEqual(merged.mean, 2.0)

    def test_absolute_contribution_gate_ignores_tiny_channel(self):
        histories = 4 * 262_144
        stats = NodeStatistics(
            "n",
            "train",
            histories,
            moment(1.0, 2.56, histories),  # SE = 0.25%, exactly on target.
            (
                ContributionMoment("material", moment(0.02, 0.5, histories)),
                ContributionMoment("tiny", moment(0.001, 100.0, histories)),
            ),
            completed_blocks=4,
        )
        result = assess_node(stats)
        checks = {item.name: item for item in result.criteria}
        self.assertAlmostEqual(checks["total_yield"].ratio, 1.0)
        self.assertTrue(checks["total_yield"].passed)
        self.assertTrue(checks["contribution:tiny"].passed)
        self.assertFalse(checks["contribution:tiny"].active)
        self.assertTrue(checks["contribution:material"].active)

    def test_exact_total_contribution_and_replica_sentinel_gates(self):
        histories = 4 * BLOCK_HISTORIES
        total_mean = 1.0
        total_target = 0.0025 * total_mean
        contribution_target = 0.001 * total_mean
        stats = NodeStatistics(
            "exact-gates",
            "train",
            histories,
            moment(total_mean, 0.9 * total_target * math.sqrt(histories), histories),
            (
                # Exactly 1% is active (the plan says "at least 1%").
                ContributionMoment(
                    "component-at-threshold",
                    moment(
                        0.01,
                        0.9 * contribution_target * math.sqrt(histories),
                        histories,
                    ),
                ),
                # A negligible cell/band has no relative-precision gate even
                # when its observed uncertainty is enormous.
                ContributionMoment(
                    "negligible-cell",
                    moment(0.009, 100.0, histories),
                ),
            ),
            (
                SentinelMoment(
                    "hpge-fold",
                    replica_moment_with_se(2.0, 0.0049),
                    2.0,  # target = 0.25% * 2 = 0.005
                ),
                SentinelMoment(
                    "joint-shape",
                    replica_moment_with_se(1.0, 0.0026),
                    1.0,  # deliberately above the 0.0025 target
                ),
            ),
            completed_blocks=4,
        )
        result = assess_node(stats)
        checks = {item.name: item for item in result.criteria}
        self.assertTrue(checks["total_yield"].passed)
        self.assertTrue(checks["contribution:component-at-threshold"].active)
        self.assertTrue(checks["contribution:component-at-threshold"].passed)
        self.assertFalse(checks["contribution:negligible-cell"].active)
        self.assertTrue(checks["contribution:negligible-cell"].passed)
        self.assertTrue(checks["sentinel:hpge-fold"].passed)
        self.assertFalse(checks["sentinel:joint-shape"].passed)
        self.assertFalse(result.passed)

    def test_sentinel_uncertainty_must_be_from_32_logical_replicas(self):
        histories = 4 * BLOCK_HISTORIES
        per_history_sentinel = NodeStatistics(
            "wrong-sentinel-n",
            "validation",
            histories,
            moment(1.0, 1.0, histories),
            sentinels=(
                SentinelMoment(
                    "nai-fold",
                    moment(1.0, 1.0, histories),
                    1.0,
                ),
            ),
            completed_blocks=4,
        )
        with self.assertRaisesRegex(ValueError, "32 logical replica means"):
            assess_node(per_history_sentinel)

    def test_initial_four_blocks_are_required_even_if_variance_is_zero(self):
        histories = BLOCK_HISTORIES
        stats = NodeStatistics(
            "interrupted-initial",
            "train",
            histories,
            moment(1.0, 0.0, histories),
            completed_blocks=1,
        )
        assessment = assess_node(stats)
        self.assertFalse(assessment.minimum_blocks_satisfied)
        self.assertFalse(assessment.passed)
        self.assertEqual(assessment.recommended_additional_blocks, 3)

    def test_role_caps_are_floored_to_complete_blocks(self):
        self.assertEqual(history_cap_blocks("train"), 244)
        self.assertEqual(effective_history_cap("train"), 63_963_136)
        self.assertEqual(history_cap_blocks("pilot"), 244)
        self.assertEqual(history_cap_blocks("diagnostic"), 244)
        self.assertEqual(history_cap_blocks("validation"), 488)
        self.assertEqual(history_cap_blocks("final"), 488)
        self.assertEqual(effective_history_cap("final"), 127_926_272)

        validation_blocks = history_cap_blocks("validation") - 1
        validation_histories = validation_blocks * BLOCK_HISTORIES
        one_block_left = assess_node(
            NodeStatistics(
                "validation-one-block-left",
                "validation",
                validation_histories,
                moment(1.0, 100.0, validation_histories),
                completed_blocks=validation_blocks,
            )
        )
        self.assertFalse(one_block_left.history_budget_exhausted)
        self.assertEqual(one_block_left.recommended_additional_blocks, 1)

        with self.assertRaisesRegex(ValueError, "integral number of blocks"):
            assess_node(
                NodeStatistics(
                    "partial-cap",
                    "train",
                    64_000_000,
                    moment(1.0, 1.0, 64_000_000),
                )
            )

    def test_history_budget_cap_remains_usable_but_transport_failure_does_not(self):
        train_histories = effective_history_cap("train")
        noisy_at_budget = NodeStatistics(
            "budget-cap",
            "train",
            train_histories,
            moment(1.0, 100.0, train_histories),
            completed_blocks=history_cap_blocks("train"),
        )
        budget = assess_node(noisy_at_budget)
        self.assertTrue(budget.history_budget_exhausted)
        self.assertTrue(budget.history_budget_capped)
        self.assertTrue(budget.usable)
        self.assertFalse(budget.transport_failed)
        self.assertFalse(budget.passed)
        self.assertEqual(budget.recommended_additional_blocks, 0)
        self.assertEqual(budget.stop_reason, "history_budget_cap")

        transport_failure = NodeStatistics(
            "transport-failure",
            "train",
            train_histories,
            moment(1.0, 100.0, train_histories),
            completed_blocks=history_cap_blocks("train"),
            transport_failed=True,
        )
        failed = assess_node(transport_failure)
        self.assertTrue(failed.transport_failed)
        self.assertTrue(failed.history_budget_exhausted)
        self.assertFalse(failed.history_budget_capped)
        self.assertFalse(failed.usable)
        self.assertFalse(failed.passed)
        self.assertEqual(failed.recommended_additional_blocks, 0)
        self.assertEqual(failed.stop_reason, "transport_failure")

    def test_scheduler_preserves_loss_reduction_per_cpu_priority(self):
        histories = 4 * BLOCK_HISTORIES
        # total SE = 0.5%, so the active total gate has ratio 2.
        continuation = NodeStatistics(
            "continue",
            "train",
            histories,
            moment(1.0, 0.005 * math.sqrt(histories), histories),
            wall_seconds=8.0,
            # Missing metadata is recoverable from the integral history count;
            # CPU utility must not silently fall back to one second per block.
            completed_blocks=0,
        )
        failed_transport = NodeStatistics(
            "failed",
            "train",
            histories,
            moment(1.0, 100.0, histories),
            wall_seconds=4.0,
            completed_blocks=4,
            transport_failed=True,
        )
        actions = allocate_actions(
            (continuation, failed_transport),
            (AcquisitionCandidate("new", 3.0, 2.0),),
            continuation_local_losses={"continue": 10.0},
        )
        self.assertEqual([item.node_id for item in actions], ["new", "continue"])
        self.assertAlmostEqual(actions[0].utility, 1.5)
        # Ratio 2 projects 16 blocks; the pass doubles the four it has, and
        # both the reduction and the CPU are priced for those four blocks.
        self.assertEqual(actions[1].blocks, 4)
        expected_reduction = 10.0 * 0.5 * (1.0 - math.sqrt(4.0 / 8.0))
        self.assertAlmostEqual(actions[1].predicted_loss_reduction, expected_reduction)
        self.assertAlmostEqual(actions[1].predicted_cpu_seconds, 8.0)
        self.assertAlmostEqual(actions[1].utility, expected_reduction / 8.0)
        self.assertNotIn("failed", {item.node_id for item in actions})

        zero_residual = allocate_actions(
            (continuation,),
            (),
            continuation_local_losses={"continue": 0.0},
        )
        self.assertGreater(zero_residual[0].predicted_loss_reduction, 0.0)
        components = continuation_loss_reduction_components(
            0.0, assess_node(continuation), 4
        )
        self.assertEqual(components.local_model_loss, 0.0)
        self.assertEqual(components.precision_guard_loss, 1.0)
        self.assertEqual(components.predicted_model_loss_reduction, 0.0)
        self.assertGreater(components.predicted_precision_guard_reduction, 0.0)

    def test_continuation_doubles_without_overshooting_the_projection(self):
        # Baseline completion first, then at most a doubling, never more than
        # the projected (cap-limited) need.
        self.assertEqual(continuation_block_count(0, 100), MIN_INITIAL_BLOCKS)
        self.assertEqual(continuation_block_count(2, 100), MIN_INITIAL_BLOCKS - 2)
        self.assertEqual(continuation_block_count(4, 100), 4)
        self.assertEqual(continuation_block_count(64, 100), 64)
        self.assertEqual(continuation_block_count(64, 3), 3)
        self.assertEqual(continuation_block_count(128, 0), 0)
        # A node needing 60x its histories reaches them in six passes.
        blocks, passes = 4, 0
        need = 60 * 4
        while blocks < need:
            blocks += continuation_block_count(blocks, need - blocks)
            passes += 1
        self.assertEqual((blocks, passes), (need, 6))
        components = continuation_loss_reduction_components(
            1.0,
            assess_node(
                NodeStatistics(
                    "n",
                    "train",
                    4 * BLOCK_HISTORIES,
                    moment(1.0, 0.005 * math.sqrt(4 * BLOCK_HISTORIES), 4 * BLOCK_HISTORIES),
                    completed_blocks=4,
                )
            ),
            4,
            4,
        )
        self.assertEqual(components.blocks, 4)
        self.assertAlmostEqual(components.noise_reduction, 1.0 - math.sqrt(0.5))

    def test_histories_and_new_node_batch_use_one_loss_per_cpu_currency(self):
        decision = compare_adaptive_batches(
            [(2.0, 4.0), (1.0, 2.0)],
            [(8.0, 8.0), (4.0, 4.0)],
        )
        self.assertEqual(decision.selected_branch, "new_node_batch")
        self.assertAlmostEqual(decision.continuation.utility, 0.5)
        self.assertAlmostEqual(decision.new_node_batch.utility, 1.0)

        forced = compare_adaptive_batches(
            [(2.0, 4.0)],
            [(100.0, 1.0)],
            forced_continuation=True,
        )
        self.assertEqual(forced.selected_branch, "continuation")

    def test_two_consecutive_validation_rounds(self):
        state = ValidationRoundState()
        self.assertFalse(state.update(True))
        self.assertFalse(state.update(False))
        self.assertFalse(state.update(True))
        self.assertTrue(state.update(True))


class CampaignTests(unittest.TestCase):
    def test_32_mu_production_requires_bound_passing_pilot_decision(self):
        production = PhysicsConfig(mu_bins=32)
        with self.assertRaisesRegex(ValueError, "64-to-32 decision"):
            generate_campaign_manifest(
                physics=production,
                candidate_count=96,
                train_count=24,
                validation_count=4,
                final_count=4,
            )

        pilot_hash = config_sha256(
            canonical_fit_config_json("field", "soil", PhysicsConfig(mu_bins=64))
        )
        decision = {
            "schema": "CEEFIT-MU-DECISION-v1",
            "config_sha256": pilot_hash,
            "source_bins": 64,
            "candidate_bins": 32,
            "tolerance": 0.002,
            "max_joint_tv": 0.001,
            "archive_snapshot_sha256": "a" * 64,
            "source_archive": {
                "path": "/local/pilot-64.h5",
                "archive_snapshot_sha256": "a" * 64,
            },
            "source_detector_validation": {
                "path": "/local/pilot-mu-fold.json",
                "report_sha256": "f" * 64,
            },
            "pilot_node_ids": ["1" * 16],
            "detector_validation": {
                "passed": True,
                "detectors": {
                    name: {
                        "response_sha256": character * 64,
                        "response_schema": "CEEFIT-DETECTOR-RESPONSE-v1",
                        "response_semantics": "probability-per-incident-cell-mass-v1",
                        "line_interpolation": "linear-exact-energy-v1",
                        "geometry_id": f"representative-{name.lower()}-fixture-v1",
                        "geometry_sha256": "d" * 64,
                        "provenance_sha256": "e" * 64,
                        "generator": {
                            "name": "synthetic-response-generator",
                            "version": "1.0",
                            "sha256": "f" * 64,
                        },
                        "dimensions": {
                            "continuum_energy_bins": 4,
                            "mu_bins": 4,
                            "spectrum_bins": 4,
                            "line_energy_grid_points": 2,
                        },
                        "fold_discretization": {
                            "method": "archive-cell-response-kernel-diameter-v1",
                            "gate": 0.002,
                            "max_continuum_total_relative_bound": 0.0,
                            "max_continuum_spectral_tv_bound": 0.0,
                            "max_exact_line_relative_bound": 0.0,
                        },
                        "node_count": 1,
                        "max_detector_folded_total_relative_change": 0.001,
                        "max_detector_folded_spectral_tv": 0.001,
                        "max_exact_line_angular_fold_relative_change": 0.001,
                        "passed": True,
                    }
                    for name, character in (("HPGe", "b"), ("NaI", "c"))
                },
            },
            "decision": "use_32",
            "passed": True,
        }
        with patch(
            "tools.full_spectrum.env_scatter_fit.campaign._verify_mu_resolution_sources"
        ) as verify_sources:
            campaign = generate_campaign_manifest(
                physics=production,
                mu_resolution_decision=decision,
                candidate_count=96,
                train_count=24,
                validation_count=4,
                final_count=4,
            )
        verify_sources.assert_called_once_with(decision)
        self.assertEqual(campaign["config"]["physics"]["mu_bins"], 32)
        self.assertEqual(
            campaign["mu_resolution_decision"]["report"]["config_sha256"],
            pilot_hash,
        )
        tampered = json.loads(json.dumps(campaign))
        tampered["mu_resolution_decision"]["report"]["max_joint_tv"] = 0.1
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "tampered.json"
            path.write_text(json.dumps(tampered), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum"):
                load_manifest(path)

    def test_default_canonical_config_hash_is_stable(self):
        canonical = canonical_fit_config_json("field", "soil")
        self.assertEqual(
            config_sha256(canonical),
            "d465fe820aac54179527eadc3caded1e922276bed97b90ab98a9709c8679ff2e",
        )

    def test_cpp_node_and_seed_vector(self):
        config_hash = "0" * 64
        node_id = PhysicalNode(122.0, 1.0, 0.0, 0.0).stable_id(config_hash)
        self.assertEqual(f"{node_id:016x}", "c5dbb8b7b2b1ad9e")
        self.assertEqual(
            f"{replica_seed(config_hash, node_id, 1, 0):016x}",
            "8fa2f48da2f95952",
        )

    def test_small_manifest_is_disjoint_and_continuable(self):
        anchors = (PhysicalNode(662.0, 5.0, 1.0, 1.0),)
        manifest = generate_campaign_manifest(
            candidate_count=160,
            train_count=32,
            validation_count=4,
            final_count=4,
            anchors=anchors,
        )
        self.assertEqual(len(manifest["work"]), 40 * 4)
        nodes = [
            *manifest["design"]["nodes"]["train"],
            *manifest["design"]["nodes"]["validation"],
            *manifest["design"]["nodes"]["final"],
        ]
        self.assertEqual(len({item["node_id"] for item in nodes}), 40)
        self.assertEqual(len(manifest["design"]["candidate_pool"]), 160)
        self.assertEqual(len(manifest["worker_groups"]), 40)
        self.assertEqual(manifest["worker_groups"][0]["block_count"], 4)
        self.assertEqual(
            manifest["worker_groups"][0]["argv"][-2:], ["--block-count", "4"]
        )
        first_work = manifest["work"][0]
        self.assertEqual(first_work["block_id"], 1)
        self.assertEqual(first_work["histories"], 262_144)
        self.assertEqual(first_work["replicas"][0]["history_count"], 8192)
        self.assertIn("--fit-block", first_work["argv"])
        self.assertIn("--role", first_work["argv"])
        self.assertEqual(
            first_work["replicas"][0]["rng_seed"],
            replica_seed(
                manifest["config"]["sha256"],
                first_work["node_id"],
                first_work["block_id"],
                0,
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.json"
            save_manifest(manifest, path)
            self.assertEqual(load_manifest(path)["config"], manifest["config"])

    def test_pair_threshold_forced_anchors_bracket_exact_threshold(self):
        energies = {node.energy_keV for node in default_forced_anchors(DesignBounds())}
        threshold = 2.0 * 510.998950
        self.assertIn(threshold - 0.1, energies)
        self.assertIn(threshold + 0.1, energies)

    def test_adaptive_round_state_batch_cap_and_two_pass_stop(self):
        manifest = generate_campaign_manifest(
            candidate_count=192,
            train_count=32,
            validation_count=8,
            final_count=8,
        )
        failed = advance_campaign_round(
            manifest,
            diagnostic_validation_failed=True,
        )
        self.assertEqual(
            failed["campaign_state"]["last_action"],
            "await_unified_adaptive_allocation",
        )
        self.assertEqual(len(failed["campaign_state"]["training_node_ids"]), 32)
        self.assertEqual(len(failed["campaign_state"]["acquired_node_ids"]), 0)
        self.assertEqual(len(failed["worker_groups"]), 48)

        first_pass = advance_campaign_round(
            manifest,
            validation_report=model_selection_report(manifest),
            expected_archive_snapshot_sha256="a" * 64,
            expected_model_content_sha256="b" * 64,
        )
        self.assertFalse(first_pass["campaign_state"]["stopped"])
        self.assertEqual(
            first_pass["campaign_state"]["last_action"],
            "await_confirmation_validation_round",
        )
        stopped = advance_campaign_round(
            first_pass,
            validation_report=model_selection_report(first_pass),
            expected_archive_snapshot_sha256="a" * 64,
            expected_model_content_sha256="b" * 64,
        )
        self.assertTrue(stopped["campaign_state"]["stopped"])
        self.assertEqual(stopped["campaign_state"]["consecutive_validation_passes"], 2)
        self.assertEqual(
            stopped["campaign_state"]["validation_history"][-1]["kind"],
            "bound_model_selection_report",
        )
        self.assertEqual(
            len(stopped["campaign_state"]["validation_history"][-1]["report_canonical_sha256"]),
            64,
        )

    def test_passing_round_rejects_unbound_or_replayed_validation_evidence(self):
        manifest = generate_campaign_manifest(
            candidate_count=96,
            train_count=24,
            validation_count=4,
            final_count=4,
        )
        report = model_selection_report(manifest)
        with self.assertRaisesRegex(ValueError, "independently opened archive"):
            advance_campaign_round(manifest, validation_report=report)
        with self.assertRaisesRegex(ValueError, "another archive"):
            advance_campaign_round(
                manifest,
                validation_report=report,
                expected_archive_snapshot_sha256="c" * 64,
                expected_model_content_sha256="b" * 64,
            )

    def test_adaptive_utility_bytes_are_bound_to_the_failed_validation_round(self):
        manifest = generate_campaign_manifest(
            candidate_count=96,
            train_count=24,
            validation_count=4,
            final_count=4,
        )
        identity = campaign_round_identity(manifest)
        utility = {
            "schema": "CEEFIT-ADAPTIVE-UTILITY-v1",
            "config_sha256": manifest["config"]["sha256"],
            "archive_snapshot_sha256": "a" * 64,
            "model_content_sha256": "b" * 64,
            "campaign_round": identity["campaign_round"],
            "campaign_round_identity_sha256": identity["identity_sha256"],
            "new_nodes": {},
            "continuation_local_validation_losses": {},
            "input_parameters": {},
            "output_parameters": {},
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "utilities.json"
            path.write_text(
                json.dumps(utility, sort_keys=True) + "\n", encoding="utf-8"
            )
            utility_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            report = model_selection_report(manifest)
            report["selection"]["mc_noise"]["validation_passed"] = False
            report["passed"] = False
            report["acquisition"] = {
                "utilities_artifact_sha256": utility_hash
            }
            failed = advance_campaign_round(
                manifest,
                validation_report=report,
                expected_archive_snapshot_sha256="a" * 64,
                expected_model_content_sha256="b" * 64,
            )
            loaded = _load_adaptive_utility_report(
                path, failed, archive_snapshot_sha256="a" * 64
            )
            self.assertEqual(loaded["report_sha256"], utility_hash)

            path.write_text(
                json.dumps({**utility, "objective": "stale"}, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "bytes differ"):
                _load_adaptive_utility_report(
                    path, failed, archive_snapshot_sha256="a" * 64
                )
        with self.assertRaisesRegex(ValueError, "another model content"):
            advance_campaign_round(
                manifest,
                validation_report=report,
                expected_archive_snapshot_sha256="a" * 64,
                expected_model_content_sha256="c" * 64,
            )

        for field, replacement, message in (
            ("config_sha256", "f" * 64, "config_sha256 binding differs"),
            ("final_generation", 12, "final_generation binding differs"),
            (
                "eligible_validation_node_ids",
                report["eligible_validation_node_ids"][:-1],
                "exact eligible validation nodes",
            ),
        ):
            tampered = json.loads(json.dumps(report))
            tampered[field] = replacement
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                advance_campaign_round(
                    manifest,
                    validation_report=tampered,
                    expected_archive_snapshot_sha256="a" * 64,
                    expected_model_content_sha256="b" * 64,
                )

        first = advance_campaign_round(
            manifest,
            validation_report=report,
            expected_archive_snapshot_sha256="a" * 64,
            expected_model_content_sha256="b" * 64,
        )
        with self.assertRaisesRegex(ValueError, "campaign_round binding differs"):
            advance_campaign_round(
                first,
                validation_report=report,
                expected_archive_snapshot_sha256="a" * 64,
                expected_model_content_sha256="b" * 64,
            )

        forged = model_selection_report(manifest)
        forged["passed"] = True
        forged["selection"]["mc_noise"]["validation_passed"] = False
        with self.assertRaisesRegex(ValueError, "pass flag differs"):
            advance_campaign_round(
                manifest,
                validation_report=forged,
                expected_archive_snapshot_sha256="a" * 64,
                expected_model_content_sha256="b" * 64,
            )

        stale_policy = model_selection_report(manifest)
        stale_policy.pop("model_selection_policy")
        with self.assertRaisesRegex(ValueError, "contribution-aware policy"):
            advance_campaign_round(
                manifest,
                validation_report=stale_policy,
                expected_archive_snapshot_sha256="a" * 64,
                expected_model_content_sha256="b" * 64,
            )

    def test_failed_final_is_replaced_from_unused_sobol_reserve(self):
        manifest = generate_campaign_manifest(
            candidate_count=256,
            train_count=32,
            validation_count=16,
            final_count=64,
        )
        first = advance_campaign_round(
            manifest,
            validation_report=model_selection_report(manifest),
            expected_archive_snapshot_sha256="a" * 64,
            expected_model_content_sha256="b" * 64,
        )
        manifest = advance_campaign_round(
            first,
            validation_report=model_selection_report(first),
            expected_archive_snapshot_sha256="a" * 64,
            expected_model_content_sha256="b" * 64,
        )
        old_ids = {item["node_id_hex"] for item in manifest["design"]["nodes"]["final"]}
        failure = final_failure_report(manifest)
        tampered = json.loads(json.dumps(failure))
        tampered["final_generation"] += 1
        with self.assertRaisesRegex(ValueError, "final_generation binding differs"):
            replace_final_set_from_reserve(manifest, failure_report=tampered)

        replaced = replace_final_set_from_reserve(
            manifest, failure_report=failure
        )
        new_ids = {item["node_id_hex"] for item in replaced["design"]["nodes"]["final"]}
        self.assertEqual(len(new_ids), 64)
        self.assertTrue(old_ids.isdisjoint(new_ids))
        self.assertEqual(
            set(replaced["campaign_state"]["retired_final_node_ids"]), old_ids
        )
        self.assertTrue(
            old_ids.issubset(set(replaced["campaign_state"]["training_node_ids"]))
        )
        self.assertTrue(
            all(item["role"] == "diagnostic" for item in replaced["design"]["nodes"]["train"][-64:])
        )
        self.assertEqual(replaced["campaign_state"]["final_generation"], 1)
        self.assertFalse(replaced["campaign_state"]["stopped"])
        self.assertEqual(replaced["campaign_state"]["consecutive_validation_passes"], 0)
        self.assertEqual(replaced["campaign_state"]["validation_history"], [])
        self.assertEqual(
            replaced["campaign_state"]["validation_epochs"][0]["validation_history"],
            manifest["campaign_state"]["validation_history"],
        )
        self.assertEqual(
            len(
                replaced["campaign_state"]["validation_epochs"][0][
                    "failure_report_sha256"
                ]
            ),
            64,
        )

        capped = json.loads(json.dumps(manifest))
        capped["campaign_state"]["final_generation"] = 2
        with self.assertRaisesRegex(ValueError, "generation cap reached"):
            replace_final_set_from_reserve(
                capped, failure_report=final_failure_report(capped)
            )


if __name__ == "__main__":
    unittest.main()
