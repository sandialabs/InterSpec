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

import contextlib
import importlib.util
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from tools.full_spectrum.env_scatter_fit.allocation import (
    BLOCK_HISTORIES,
    CriterionResult,
    NodeStatistics,
    RawMoment,
    assess_node,
    history_cap_blocks,
)
from tools.full_spectrum.env_scatter_fit.campaign import (
    continuation_shell_commands,
    generate_campaign_manifest,
    replica_seed,
    save_manifest,
    schedule_archive_continuation,
)
from tools.full_spectrum.env_scatter_fit.cli import main
from tools.full_spectrum.env_scatter_fit.workflow import (
    _score_final_set,
    _transport_status_report,
    _validation_mc_noise_report,
)


H5PY_AVAILABLE = importlib.util.find_spec("h5py") is not None


class AdditionalCriterionTests(unittest.TestCase):
    def test_archive_derived_criterion_participates_in_block_recommendation(self):
        histories = 4 * BLOCK_HISTORIES
        stats = NodeStatistics(
            "archive-derived",
            "train",
            histories,
            RawMoment(float(histories), float(histories), histories),
            completed_blocks=4,
            additional_criteria=(
                CriterionResult("joint_continuum_shape", 0.005, 0.0025, True),
            ),
        )
        assessment = assess_node(stats)
        self.assertFalse(assessment.passed)
        self.assertEqual(assessment.worst_ratio, 2.0)
        self.assertEqual(assessment.recommended_additional_blocks, 12)


def _one_node_campaign():
    campaign = generate_campaign_manifest(
        mode="air",
        candidate_count=96,
        train_count=24,
        validation_count=4,
        final_count=4,
    )
    node = campaign["design"]["nodes"]["train"][0]
    campaign["design"]["nodes"] = {"train": [node], "validation": [], "final": []}
    campaign["campaign_state"]["training_node_ids"] = [node["node_id_hex"]]
    campaign["campaign_state"]["acquired_node_ids"] = []
    campaign["work"] = [
        item for item in campaign["work"] if item["node_id_hex"] == node["node_id_hex"]
    ]
    campaign["worker_groups"] = [
        item
        for item in campaign["worker_groups"]
        if item["node_id_hex"] == node["node_id_hex"]
    ]
    return campaign, node


def _raw_second(first, count, standard_error):
    return first * first / count + standard_error**2 * count * (count - 1)


def _write_archive(
    path,
    campaign,
    node,
    block_ids,
    *,
    total_relative_se=0.005,
    status="complete",
    transport_failures=0,
    capped=False,
    archive_role=None,
):
    import h5py
    import numpy as np

    block_ids = tuple(block_ids)
    histories = len(block_ids) * BLOCK_HISTORIES
    replicas = 32
    replica_histories = histories // replicas if histories else 0
    text_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as root:
        root.attrs["schema"] = "CEEFIT-HDF5-v1"
        manifest = root.create_group("manifest")
        manifest.attrs["config_sha256"] = campaign["config"]["sha256"]

        nodes = root.create_group("nodes")
        nodes.create_dataset("id", data=[int(node["node_id"])], dtype="u8")
        nodes.create_dataset(
            "coordinates",
            data=[
                [
                    node["energy_keV"],
                    node["distance_m"] * 100.0,
                    node["source_height_m"] * 100.0,
                    node["detector_height_m"] * 100.0,
                    0.0,
                ]
            ],
            dtype="f8",
        )
        nodes.create_dataset(
            "role", data=[archive_role or node["role"]], dtype=text_dtype
        )
        nodes.create_dataset("status", data=[status], dtype=text_dtype)
        nodes.create_dataset("histories", data=[histories], dtype="u8")
        # The exact elastic environmental atom is assessed together with the
        # deterministic direct E0 atom, as it is in the delivered spectrum.
        nodes.create_dataset("direct_fluence", data=[1.0], dtype="f8")
        nodes.create_dataset("wall_seconds", data=[8.0 * len(block_ids)], dtype="f8")
        nodes.create_dataset("cpu_seconds", data=[4.0 * len(block_ids)], dtype="f8")
        nodes.create_dataset(
            "transport_failures", data=[transport_failures], dtype="u8"
        )
        nodes.create_dataset("capped", data=[int(capped)], dtype="u1")

        axes = root.create_group("axes")
        axes.create_dataset("mu_edges", data=[-1.0, 0.0, 1.0], dtype="f8")
        axes.attrs["line_definitions_json"] = json.dumps(
            [
                {"name": "elastic_E0", "origin": 0, "post_creation_order": 3},
                {
                    "name": "annihilation_511",
                    "origin": 1,
                    "post_creation_order": 3,
                },
            ]
        )
        energy = axes.create_group("energy").create_group("0")
        energy.create_dataset(
            "physical_edges_keV",
            data=[[10.0, 50.0, float(node["energy_keV"])]],
            dtype="f8",
        )

        components = root.create_group("components")
        components.attrs["descriptors_json"] = json.dumps(
            [{"name": "primary_order1", "energy_axis": 0}]
        )

        continuum = root.create_group("continuum")
        total = continuum.create_group("total")
        cell_n = np.full((1, 4), histories, dtype=np.uint64)
        cell_s1 = np.asarray([[0.25 * histories] * 4], dtype=float)
        replica_map = np.full((1, replicas, 4), 0.25 * replica_histories, dtype=float)
        total.create_dataset("N", data=cell_n)
        total.create_dataset("S1", data=cell_s1)
        total.create_dataset("replica_S1", data=replica_map)
        component = continuum.create_group("components")
        component.create_dataset("N", data=cell_n)
        component.create_dataset("S1", data=cell_s1)
        component.create_dataset("replica_S1", data=replica_map)

        lines = root.create_group("lines")
        lines.create_dataset("support_mask", data=[[0, 0]], dtype="u1")
        line_total = lines.create_group("total")
        line_total.create_dataset("N", data=np.full((1, 4), histories, dtype=np.uint64))
        line_total.create_dataset("S1", data=np.zeros((1, 4), dtype=float))
        line_total.create_dataset(
            "replica_S1", data=np.zeros((1, replicas, 4), dtype=float)
        )
        line_components = lines.create_group("components")
        line_components.create_dataset(
            "N", data=np.full((1, 4), histories, dtype=np.uint64)
        )
        line_components.create_dataset("S1", data=np.zeros((1, 4), dtype=float))
        line_components.create_dataset(
            "replica_S1", data=np.zeros((1, replicas, 4), dtype=float)
        )

        replica_group = root.create_group("replicas")
        replica_group.create_dataset(
            "histories",
            data=np.full((1, replicas), replica_histories, dtype=np.uint64),
        )

        stopping = root.create_group("stopping")
        stopping.create_dataset("N", data=[histories], dtype="u8")
        stopping_s1 = np.zeros((1, 16), dtype=float)
        stopping_s1[0, 0] = 0.02 * histories
        stopping.create_dataset("S1", data=stopping_s1)
        outer = np.zeros((1, 16, 16), dtype=float)
        if histories:
            outer[0, 0, 0] = stopping_s1[0, 0] ** 2 / histories
        stopping.create_dataset("outer_S2", data=outer)

        diagnostics = root.create_group("diagnostics")
        first = float(histories)
        diagnostics.create_dataset("scattered_total_N", data=[histories], dtype="u8")
        diagnostics.create_dataset("scattered_total_S1", data=[first], dtype="f8")
        diagnostics.create_dataset(
            "scattered_total_S2",
            data=[
                _raw_second(first, histories, total_relative_se)
                if histories
                else 0.0
            ],
            dtype="f8",
        )

        blocks = root.create_group("blocks").create_group(node["node_id_hex"])
        histories_per_replica = BLOCK_HISTORIES // replicas
        for block_id in block_ids:
            block = blocks.create_group(str(block_id))
            block.attrs["histories"] = BLOCK_HISTORIES
            block.attrs["payload_sha256"] = f"{block_id:064x}"
            block.create_dataset("id", data=list(range(replicas)), dtype="u8")
            block.create_dataset(
                "history_begin",
                data=[
                    (block_id - 1) * BLOCK_HISTORIES
                    + replica * histories_per_replica
                    for replica in range(replicas)
                ],
                dtype="u8",
            )
            block.create_dataset(
                "history_count",
                data=[histories_per_replica] * replicas,
                dtype="u8",
            )


@unittest.skipUnless(H5PY_AVAILABLE, "h5py archive extra is not installed")
class ArchiveContinuationTests(unittest.TestCase):
    def test_max_vertex_capped_final_is_rejected_before_residual_evaluation(self):
        with tempfile.TemporaryDirectory() as temporary:
            import h5py
            import numpy as np

            campaign, node = _one_node_campaign()
            archive = Path(temporary) / "capped-final.h5"
            _write_archive(
                archive,
                campaign,
                node,
                (1, 2, 3, 4),
                status="capped",
                capped=True,
                archive_role="final",
            )

            class PayloadMustNotRun:
                def evaluate(self, *args, **kwargs):
                    raise AssertionError("invalid final truth must not be scored")

            with h5py.File(archive, "r") as root:
                status = _transport_status_report(root, np.asarray([0]))
                report = _score_final_set(
                    root,
                    PayloadMustNotRun(),
                    np.asarray([0]),
                    [{"name": "primary_order1", "energy_axis": 0}],
                    [0, 4],
                    detector_validation_path=None,
                    config_hash=campaign["config"]["sha256"],
                    model_content_sha256="a" * 64,
                    archive_snapshot_hash="b" * 64,
                )

            self.assertFalse(status["passed"])
            self.assertEqual(status["failed_node_ids"], [node["node_id_hex"]])
            self.assertFalse(report["evaluated"])
            self.assertFalse(report["gates"]["transport_complete"])
            self.assertFalse(report["detector_validation"]["supplied"])

    def test_component_noise_significance_includes_its_exact_line_atom(self):
        with tempfile.TemporaryDirectory() as temporary:
            import h5py
            import numpy as np

            campaign, node = _one_node_campaign()
            archive = Path(temporary) / "line-component.h5"
            _write_archive(archive, campaign, node, (1, 2, 3, 4))
            histories = 4 * BLOCK_HISTORIES
            per_replica = histories // 32
            descriptor = {
                "name": "primary_order0",
                "origin": 0,
                "post_creation_order": 0,
                "energy_axis": 0,
            }
            with h5py.File(archive, "r+") as root:
                root["components"].attrs["descriptors_json"] = json.dumps([descriptor])
                root["continuum/components/S1"][0] = 0.0
                root["continuum/components/replica_S1"][0] = 0.0
                root["lines/support_mask"][0] = [1, 0]
                line_first = np.asarray([0.5, 0.5, 0.0, 0.0]) * histories
                root["lines/components/S1"][0] = line_first
                root["lines/total/S1"][0] = line_first
                replica_first = np.tile(
                    np.asarray([0.5, 0.5, 0.0, 0.0]) * per_replica,
                    (32, 1),
                )
                root["lines/components/replica_S1"][0] = replica_first
                root["lines/total/replica_S1"][0] = replica_first
                report = _validation_mc_noise_report(
                    root,
                    np.asarray([0]),
                    [descriptor],
                    [0, 4],
                )
            scheduled = schedule_archive_continuation(campaign, archive, max_nodes=1)

            continuum = report["nodes"][0]["physical_components"][0]
            component = report["nodes"][0]["exact_line_atoms"][0]
            criteria = {
                item["name"]: item
                for item in scheduled["campaign_state"]["continuation_schedule"][
                    "assessments"
                ][0]["criteria"]
            }
            self.assertEqual(
                report["component_acceptance_semantics"],
                "continuum-origin-order-plus-exact-atom-sentinel3-v1",
            )
            self.assertEqual(continuum["truth_yield"], 0.0)
            self.assertEqual(continuum["post_creation_order"], 0)
            self.assertEqual(component["post_creation_order"], 3)
            self.assertAlmostEqual(component["truth_yield"], 1.0)
            self.assertTrue(component["yield_noise_gate_active"])
            self.assertTrue(component["shape_noise_gate_active"])
            self.assertEqual(component["absolute_replica_se_fraction"], 0.0)
            self.assertFalse(criteria["component_yield:primary_order0"]["active"])
            self.assertTrue(
                criteria["component_yield:exact_line:elastic_E0"]["active"]
            )
            self.assertTrue(
                criteria["line_combined_yield:elastic_E0"]["active"]
            )
            self.assertTrue(
                criteria["line_combined_shape_2:elastic_E0"]["active"]
            )

    def test_noisy_node_at_most_doubles_and_resume_is_idempotent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            archive = root / "fit.h5"
            _write_archive(archive, campaign, node, (1, 2, 3, 4))

            updated = schedule_archive_continuation(campaign, archive, max_nodes=1)
            schedule = updated["campaign_state"]["continuation_schedule"]
            # Four archived blocks: the pass doubles the node, no further.
            self.assertEqual(
                schedule["pending_block_ids_by_node"], {node["node_id_hex"]: [5, 6, 7, 8]}
            )
            self.assertEqual(len(schedule["pending_work"]), 4)
            work = schedule["pending_work"][0]
            self.assertEqual(work["block_id"], 5)
            self.assertEqual(
                work["replicas"][0]["history_begin"], 4 * BLOCK_HISTORIES
            )
            self.assertEqual(
                work["replicas"][0]["rng_seed"],
                replica_seed(
                    campaign["config"]["sha256"], int(node["node_id"]), 5, 0
                ),
            )
            self.assertEqual(len(updated["work"]), 8)
            self.assertEqual(len(continuation_shell_commands(updated)), 1)

            resumed = schedule_archive_continuation(updated, archive, max_nodes=1)
            resumed_schedule = resumed["campaign_state"]["continuation_schedule"]
            self.assertEqual(resumed_schedule["generation"], schedule["generation"])
            self.assertEqual(resumed_schedule["pending_block_keys"], schedule["pending_block_keys"])
            self.assertEqual(len(resumed["work"]), 8)
            self.assertEqual(
                resumed["campaign_state"]["last_action"],
                "archive_continuation_unchanged",
            )

    def test_interrupted_initial_allocation_reuses_unarchived_disjoint_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            archive = root / "partial.h5"
            _write_archive(archive, campaign, node, (1, 3))

            updated = schedule_archive_continuation(campaign, archive, max_nodes=1)
            schedule = updated["campaign_state"]["continuation_schedule"]
            self.assertEqual(
                schedule["archived_block_ids_by_node"][node["node_id_hex"]],
                [1, 3],
            )
            self.assertEqual(
                schedule["pending_block_ids_by_node"][node["node_id_hex"]],
                [2, 4],
            )
            self.assertEqual(
                [item["block_id"] for item in schedule["pending_work"]], [2, 4]
            )
            self.assertEqual(len(schedule["pending_worker_groups"]), 2)

    def test_training_history_cap_is_persisted_without_more_work(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            cap = history_cap_blocks("train")
            archive = root / "capped.h5"
            _write_archive(archive, campaign, node, range(1, cap + 1))

            updated = schedule_archive_continuation(campaign, archive, max_nodes=1)
            schedule = updated["campaign_state"]["continuation_schedule"]
            self.assertEqual(schedule["pending_work"], [])
            assessment = schedule["assessments"][0]
            self.assertTrue(assessment["history_budget_exhausted"])
            self.assertTrue(assessment["history_budget_capped"])
            self.assertEqual(assessment["stop_reason"], "history_budget_cap")
            self.assertEqual(
                updated["campaign_state"]["last_action"],
                "archive_precision_complete_or_capped",
            )

    def test_transport_capped_node_is_never_rescheduled(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            archive = root / "transport-capped.h5"
            _write_archive(
                archive,
                campaign,
                node,
                (1, 2, 3, 4),
                status="capped",
                capped=True,
            )

            updated = schedule_archive_continuation(campaign, archive, max_nodes=1)
            schedule = updated["campaign_state"]["continuation_schedule"]
            self.assertEqual(schedule["pending_work"], [])
            assessment = schedule["assessments"][0]
            self.assertTrue(assessment["transport_failed"])
            self.assertEqual(assessment["stop_reason"], "transport_failure")

    def test_retired_final_above_training_cap_remains_usable_without_more_work(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            node["role"] = "diagnostic"
            cap = history_cap_blocks("diagnostic")
            archive = root / "retired-final.h5"
            _write_archive(
                archive,
                campaign,
                node,
                range(1, cap + 2),
                archive_role="final",
            )

            updated = schedule_archive_continuation(campaign, archive, max_nodes=1)
            schedule = updated["campaign_state"]["continuation_schedule"]
            self.assertEqual(schedule["pending_work"], [])
            assessment = schedule["assessments"][0]
            self.assertEqual(assessment["role"], "diagnostic")
            self.assertEqual(assessment["archive_role"], "final")
            self.assertTrue(assessment["usable"])
            self.assertTrue(assessment["history_budget_exhausted"])
            self.assertTrue(assessment["history_budget_capped"])

    def test_cli_writes_manifest_and_only_pending_commands(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            manifest = root / "campaign.json"
            archive = root / "fit.h5"
            output = root / "continued.json"
            commands = root / "commands.sh"
            save_manifest(campaign, manifest)
            _write_archive(archive, campaign, node, (1, 2, 3, 4))

            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                result = main(
                    [
                        "continue-archive",
                        "--manifest",
                        str(manifest),
                        "--archive",
                        str(archive),
                        "--output",
                        str(output),
                        "--commands-output",
                        str(commands),
                        "--max-nodes",
                        "1",
                    ]
                )
            self.assertEqual(result, 0)
            self.assertTrue(output.is_file())
            command_text = commands.read_text(encoding="utf-8")
            self.assertIn("--block-id 5", command_text)
            self.assertIn("--block-count 4", command_text)
            self.assertNotIn("--block-id 1 ", command_text)
            self.assertIn("blocks=4", stderr.getvalue())

    def test_projection_sums_every_node_and_writes_nothing_else(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign, node = _one_node_campaign()
            manifest = root / "campaign.json"
            archive = root / "fit.h5"
            output = root / "projection.json"
            save_manifest(campaign, manifest)
            _write_archive(archive, campaign, node, (1, 2, 3, 4))
            before = sorted(path.name for path in root.iterdir())

            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                result = main(
                    [
                        "project-continuation",
                        "--manifest",
                        str(manifest),
                        "--archive",
                        str(archive),
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(result, 0)
            report = json.loads(output.read_text(encoding="utf-8"))
            scheduled = schedule_archive_continuation(campaign, archive)
            expected = sum(
                int(row["recommended_additional_blocks"])
                for row in scheduled["campaign_state"]["continuation_schedule"]["assessments"]
            )
            # The whole remaining recommendation, not one pass's single block.
            self.assertEqual(report["additional_blocks"], expected)
            self.assertGreater(expected, 1)
            # Blocks 5.. fold into the node's state: no raw records, a state
            # bounded by FOLD_STATE_TO_RAW_RECORD_BOUND of the node's records.
            growth = report["archive_growth"]
            self.assertEqual(growth["raw_record_bytes"], 0.0)
            self.assertGreater(growth["fold_state_bytes"], 0.0)
            self.assertEqual(report["additional_archive_bytes"], growth["total_bytes"])
            self.assertLessEqual(growth["expected_bytes"], growth["total_bytes"])
            self.assertIn("at most +", stdout.getvalue())
            self.assertIn(f"total +{expected} blocks", stdout.getvalue())
            self.assertEqual(
                sorted(path.name for path in root.iterdir()),
                sorted(before + [output.name]),
            )


if __name__ == "__main__":
    unittest.main()
