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
from unittest import mock

from tools.full_spectrum.env_scatter_fit.archive import _aligned_fit_cell_chunk, ingest_blocks
from tools.full_spectrum.env_scatter_fit.design import PhysicalNode
from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block


H5PY_AVAILABLE = importlib.util.find_spec("h5py") is not None


@unittest.skipUnless(H5PY_AVAILABLE, "h5py archive extra is not installed")
class ArchiveHdf5Tests(unittest.TestCase):
    def test_fit_cell_chunks_are_bounded_and_logically_aligned(self):
        self.assertEqual(_aligned_fit_cell_chunk(384 * 64), 8192)
        self.assertEqual(_aligned_fit_cell_chunk(320 * 64), 5120)
        for width in (64, 511, 20_480, 24_576, 24_577):
            chunk = _aligned_fit_cell_chunk(width)
            self.assertLessEqual(chunk, 8192)
            self.assertEqual(width % chunk, 0)

    def test_new_archive_preflight_is_atomic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first.ceefit"
            mismatch = root / "mismatch.ceefit"
            first.write_bytes(fixture_block(mu_bins=32))
            mismatch.write_bytes(fixture_block(mu_bins=64, block_id=2))
            output = root / "archive.h5"
            with self.assertRaisesRegex(ValueError, "layout/configuration"):
                ingest_blocks((first, mismatch), output)
            self.assertFalse(output.exists())

            with self.assertRaisesRegex(ValueError, "duplicate block"):
                ingest_blocks((first, first), output)
            self.assertFalse(output.exists())

    def test_verified_ingest_can_remove_redundant_worker_block(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "block.ceefit"
            source.write_bytes(fixture_block())
            archive = root / "archive.h5"

            result = ingest_blocks(
                (source,), archive, delete_source_blocks=True
            )

            self.assertEqual(result["sources_deleted"], 1)
            self.assertFalse(source.exists())
            with h5py.File(archive, "r") as stored:
                node_hex = f"{int(stored['nodes/id'][0]):016x}"
                self.assertIn("1", stored[f"blocks/{node_hex}"])

    def test_failed_post_write_verification_preserves_worker_block(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "block.ceefit"
            source.write_bytes(fixture_block())
            archive = root / "archive.h5"

            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.archive._verify_materialized_node",
                side_effect=ValueError("injected post-write verification failure"),
            ):
                with self.assertRaisesRegex(ValueError, "injected post-write"):
                    ingest_blocks(
                        (source,), archive, delete_source_blocks=True
                    )

            self.assertTrue(source.exists())
            self.assertTrue(archive.exists())

    def _snapshot(self, path):
        import h5py
        import numpy as np

        result = {}
        with h5py.File(path, "r") as root:
            def collect(name, value):
                if not isinstance(value, h5py.Dataset) or name.startswith("blocks/"):
                    return
                if h5py.check_string_dtype(value.dtype) is not None:
                    result[name] = tuple(value.asstr()[...].reshape(-1).tolist())
                else:
                    result[name] = np.asarray(value[...]).copy()

            root.visititems(collect)
        return result

    def test_incremental_permutations_rematerialize_identical_science(self):
        import h5py
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for block_id, scale in ((1, 1.0e8), (2, 1.0), (3, 1.0e-4)):
                path = root / f"block_{block_id}.ceefit"
                path.write_bytes(fixture_block(block_id=block_id, scale=scale))
                paths.append(path)
            archive_a = root / "a.h5"
            archive_b = root / "b.h5"
            ingest_blocks((paths[2], paths[0]), archive_a)
            ingest_blocks((paths[1],), archive_a)
            ingest_blocks((paths[0],), archive_b)
            ingest_blocks((paths[1], paths[2]), archive_b)

            first = self._snapshot(archive_a)
            second = self._snapshot(archive_b)
            self.assertEqual(first.keys(), second.keys())
            for name in first:
                if isinstance(first[name], tuple):
                    self.assertEqual(first[name], second[name], name)
                else:
                    np.testing.assert_array_equal(first[name], second[name], err_msg=name)
            with h5py.File(archive_a, "r") as first_archive, h5py.File(
                archive_b, "r"
            ) as second_archive:
                self.assertEqual(
                    archive_snapshot_sha256(first_archive),
                    archive_snapshot_sha256(second_archive),
                )
            with h5py.File(archive_a, "r") as archive:
                self.assertEqual(int(archive["nodes/histories"][0]), 3 * 262_144)
                self.assertEqual(tuple(archive["lines/support_mask"][0]), (1, 0))
                block_group = archive[f"blocks/{int(archive['nodes/id'][0]):016x}/1"]
                self.assertIn("moments", block_group)
                self.assertIn("coordinates", block_group)
                self.assertIn("q", block_group)
                self.assertIn("physical_energy_edges_keV", block_group)
                self.assertIn("component_support_mask", block_group)
                self.assertIn("line_energy_keV", block_group)
                self.assertEqual(
                    block_group.attrs["record_schema"],
                    "CEEFIT-HDF5-BLOCK-RECORD-v1",
                )
                self.assertEqual(len(block_group.attrs["scientific_record_sha256"]), 64)
                self.assertEqual(
                    block_group["moments/continuum_components/S1"].compression,
                    "gzip",
                )
                self.assertEqual(
                    block_group[
                        "moments/continuum_components/S1"
                    ].compression_opts,
                    1,
                )
                self.assertFalse(
                    block_group["moments/continuum_components/S1"].shuffle
                )
                self.assertEqual(
                    block_group["replica_s1/continuum_components"].chunks[0], 1
                )
                self.assertLessEqual(
                    archive["continuum/components/replica_S1"].chunks[-1],
                    8192,
                )
                self.assertEqual(
                    archive["continuum/components/replica_S1"].compression_opts,
                    1,
                )
                component_cells = (
                    archive["continuum/components/replica_S1"].shape[-1] // 9
                )
                self.assertEqual(
                    component_cells
                    % archive["continuum/components/replica_S1"].chunks[-1],
                    0,
                )
                self.assertEqual(archive["nodes/status"].asstr()[0], "complete")
                self.assertEqual(
                    archive["stopping"].attrs["definition_json"],
                    block_group.attrs["stopping_definition_json"],
                )
                self.assertIn(
                    '"name":"registered_band_00"',
                    archive["stopping"].attrs["definition_json"],
                )
                self.assertIn(
                    '"name":"legendre_l4"',
                    archive["stopping"].attrs["definition_json"],
                )
                stopping_n = int(archive["stopping/N"][0])
                stopping_s1 = np.asarray(archive["stopping/S1"][0])
                stopping_outer = np.asarray(archive["stopping/outer_S2"][0])
                expected_se = np.sqrt(
                    np.maximum(
                        0.0,
                        np.diag(stopping_outer) - stopping_s1**2 / stopping_n,
                    )
                    / (stopping_n * (stopping_n - 1))
                )
                np.testing.assert_array_equal(
                    archive["stopping/standard_error"][0], expected_se
                )
                total_mean = (
                    float(archive["diagnostics/scattered_total_S1"][0])
                    / int(archive["diagnostics/scattered_total_N"][0])
                )
                np.testing.assert_array_equal(
                    archive["stopping/target_standard_error"][0],
                    np.asarray([0.001 * total_mean] * 12 + [0.0025 * total_mean] * 4),
                )
                active = np.asarray(archive["stopping/active"][0])
                expected_active = np.zeros(16, dtype="u1")
                expected_active[:12] = (
                    stopping_s1[:12] / stopping_n >= 0.01 * total_mean
                )
                expected_active[12:] = 1
                np.testing.assert_array_equal(active, expected_active)
                np.testing.assert_array_equal(
                    archive["stopping/passed"][0],
                    np.logical_or(active == 0, expected_se <= archive["stopping/target_standard_error"][0]),
                )
            with self.assertRaisesRegex(ValueError, "already archived"):
                ingest_blocks((paths[0],), archive_a)
            mismatched = root / "mismatched_mu.ceefit"
            mismatched.write_bytes(fixture_block(mu_bins=64))
            with self.assertRaisesRegex(ValueError, "configuration/layout"):
                ingest_blocks((mismatched,), archive_a)
            with h5py.File(archive_a, "r") as archive:
                self.assertEqual(len(archive["nodes/id"]), 1)
                self.assertEqual(len(archive["blocks"][next(iter(archive["blocks"]))]), 3)

    def test_snapshot_rejects_materialized_or_raw_mutation(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            block = root / "block.ceefit"
            block.write_bytes(fixture_block())
            archive_path = root / "archive.h5"
            ingest_blocks((block,), archive_path)

            with h5py.File(archive_path, "r+") as archive:
                baseline = archive_snapshot_sha256(archive)

                aggregate = archive["continuum/total/S1"]
                aggregate_value = float(aggregate[0, 0])
                aggregate[0, 0] = aggregate_value + 0.125
                with self.assertRaisesRegex(ValueError, "does not close to raw blocks"):
                    archive_snapshot_sha256(archive)
                aggregate[0, 0] = aggregate_value
                self.assertEqual(baseline, archive_snapshot_sha256(archive))

                node_hex = f"{int(archive['nodes/id'][0]):016x}"
                raw = archive[f"blocks/{node_hex}/1/moments/continuum_total/S1"]
                raw_value = float(raw[0])
                raw[0] = raw_value + 0.125
                with self.assertRaisesRegex(ValueError, "scientific checksum"):
                    archive_snapshot_sha256(archive)

            orphan_path = root / "orphan.h5"
            ingest_blocks((block,), orphan_path)
            with h5py.File(orphan_path, "r+") as orphan:
                node_hex = f"{int(orphan['nodes/id'][0]):016x}"
                orphan["blocks"].move(node_hex, "0000000000000001")
                with self.assertRaisesRegex(
                    ValueError, "raw-block and materialized node inventories differ"
                ):
                    archive_snapshot_sha256(orphan)

    def test_snapshot_closure_covers_metadata_moments_replicas_and_diagnostics(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            block = root / "block.ceefit"
            block.write_bytes(fixture_block())
            archive_path = root / "archive.h5"
            ingest_blocks((block,), archive_path)
            with h5py.File(archive_path, "r+") as archive:
                baseline = archive_snapshot_sha256(archive)
                numeric_mutations = (
                    ("nodes/coordinates", (0, 0), 0.25),
                    ("nodes/q", (0, 1), 0.25),
                    ("axes/energy/0/physical_edges_keV", (0, 1), 0.25),
                    ("components/support_mask", (0, 0), 1),
                    ("lines/energy_keV", (0, 0), 0.25),
                    ("continuum/components/N", (0, 0), 1),
                    ("continuum/total/S2", (0, 0), 0.25),
                    ("continuum/total/replica_S1", (0, 0, 0), 0.25),
                    ("stopping/standard_error", (0, 0), 0.25),
                    ("diagnostics/interaction_counts", (0, 0), 1),
                    (
                        "diagnostics/path_channels/air_only/continuum/S1",
                        (0, 0),
                        0.25,
                    ),
                )
                for path, location, increment in numeric_mutations:
                    with self.subTest(path=path):
                        dataset = archive[path]
                        original = dataset[location]
                        dataset[location] = original + increment
                        with self.assertRaisesRegex(
                            ValueError, "does not close to raw blocks"
                        ):
                            archive_snapshot_sha256(archive)
                        dataset[location] = original
                        self.assertEqual(baseline, archive_snapshot_sha256(archive))

                definition = archive["stopping"].attrs["definition_json"]
                archive["stopping"].attrs["definition_json"] = "{}"
                with self.assertRaisesRegex(ValueError, "stopping/definition_json"):
                    archive_snapshot_sha256(archive)
                archive["stopping"].attrs["definition_json"] = definition
                self.assertEqual(baseline, archive_snapshot_sha256(archive))

                role = archive["nodes/role"]
                original_role = role.asstr()[0]
                role[0] = "final"
                with self.assertRaisesRegex(ValueError, "does not close to raw blocks"):
                    archive_snapshot_sha256(archive)
                role[0] = original_role
                self.assertEqual(baseline, archive_snapshot_sha256(archive))

                mu_edges = archive["axes/mu_edges"]
                original_edge = float(mu_edges[1])
                mu_edges[1] = original_edge + 0.01
                with self.assertRaisesRegex(ValueError, "does not close to raw blocks"):
                    archive_snapshot_sha256(archive)

    def test_incomplete_raw_record_requires_reingest(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            block = root / "block.ceefit"
            block.write_bytes(fixture_block())
            archive_path = root / "archive.h5"
            ingest_blocks((block,), archive_path)
            with h5py.File(archive_path, "r+") as archive:
                node_hex = f"{int(archive['nodes/id'][0]):016x}"
                del archive[f"blocks/{node_hex}/1"].attrs["record_schema"]
                with self.assertRaisesRegex(ValueError, "re-ingest"):
                    archive_snapshot_sha256(archive)

    def test_snapshot_is_invariant_to_two_node_ingest_order(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_node = root / "near.ceefit"
            second_node = root / "far.ceefit"
            first_node.write_bytes(
                fixture_block(physical_node=PhysicalNode(100.0, 1.0, 0.0, 0.0))
            )
            second_node.write_bytes(
                fixture_block(physical_node=PhysicalNode(100.0, 2.0, 0.0, 0.0))
            )
            forward = root / "forward.h5"
            reverse = root / "reverse.h5"
            ingest_blocks((first_node, second_node), forward)
            ingest_blocks((second_node, first_node), reverse)

            with h5py.File(forward, "r") as forward_archive, h5py.File(
                reverse, "r"
            ) as reverse_archive:
                self.assertNotEqual(
                    list(forward_archive["nodes/id"][:]),
                    list(reverse_archive["nodes/id"][:]),
                )
                self.assertEqual(
                    archive_snapshot_sha256(forward_archive),
                    archive_snapshot_sha256(reverse_archive),
                )

    def test_mu_resolution_workflow_runs_on_real_hdf5(self):
        import h5py
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import _transport_complete_mask
        from tools.full_spectrum.env_scatter_fit.workflow import mu_resolution_report

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            block = root / "pilot.ceefit"
            block.write_bytes(fixture_block(mu_bins=64, role="pilot"))
            archive = root / "pilot.h5"
            ingest_blocks((block,), archive)
            with h5py.File(archive, "r+") as stored:
                self.assertTrue(bool(_transport_complete_mask(stored, np)[0]))
                stored["nodes/status"][0] = "capped"
                self.assertFalse(bool(_transport_complete_mask(stored, np)[0]))
                stored["nodes/status"][0] = "complete"
            report = mu_resolution_report(archive)
        self.assertEqual(report["source_bins"], 64)
        self.assertEqual(report["candidate_bins"], 32)
        self.assertIn(report["decision"], {"use_32", "retain_64"})

    def test_capped_continuation_status_is_preserved_and_excluded_from_truth(self):
        import h5py
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.workflow import _transport_complete_mask

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            capped = root / "capped.ceefit"
            continued = root / "continued.ceefit"
            capped.write_bytes(
                fixture_block(
                    block_id=1,
                    status="capped",
                    capped=True,
                    failure_reason="maximum vertex cap reached",
                )
            )
            continued.write_bytes(fixture_block(block_id=2))
            archive = root / "status.h5"
            ingest_blocks((continued, capped), archive)
            with h5py.File(archive, "r") as stored:
                self.assertEqual(stored["nodes/status"].asstr()[0], "capped")
                self.assertIn("maximum vertex", stored["nodes/failure_reason"].asstr()[0])
                self.assertFalse(bool(_transport_complete_mask(stored, np)[0]))


if __name__ == "__main__":
    unittest.main()
