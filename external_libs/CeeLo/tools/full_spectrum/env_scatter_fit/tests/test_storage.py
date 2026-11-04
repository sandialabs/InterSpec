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
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from tools.full_spectrum.env_scatter_fit.archive import ingest_blocks
from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block


H5PY_AVAILABLE = importlib.util.find_spec("h5py") is not None


class StorageCliAndOrchestrationTests(unittest.TestCase):
    def test_parser_exposes_bounded_storage_workflow(self):
        import argparse

        from tools.full_spectrum.env_scatter_fit.cli import build_parser

        parser = build_parser()
        subparsers = next(
            action
            for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        for command in (
            "run-streaming",
            "seal-archive",
            "verify-archive-seal",
            "fit-cache",
            "verify-fit-cache",
        ):
            with self.subTest(command=command):
                self.assertIn(command, subparsers.choices)

    def test_streaming_runs_worker_then_verified_background_pack(self):
        from tools.full_spectrum.env_scatter_fit.storage import run_streaming_campaign

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            spool = root / "spool"
            archive = root / "archive.h5"
            manifest = {
                "work": [
                    {
                        "node_id": 7,
                        "block_id": 1,
                        "argv": ["worker", "--fit-block", "unused.ceefit"],
                    }
                ],
                "worker_groups": [
                    {
                        "node_id": 7,
                        "role": "train",
                        "first_block_id": 1,
                        "block_count": 1,
                        "argv": ["worker", "--fit-block", "unused.ceefit"],
                    }
                ],
            }

            def fake_worker(argv, check):
                self.assertTrue(check)
                output = Path(argv[argv.index("--fit-block") + 1])
                output = Path(str(output).format(block_id=1))
                output.write_bytes(b"worker result")

            def fake_ingest(paths, destination, *, delete_source_blocks):
                self.assertEqual(Path(destination), archive.resolve())
                self.assertEqual(len(paths), 1)
                self.assertTrue(paths[0].is_file())
                if delete_source_blocks:
                    paths[0].unlink()
                return {
                    "blocks_ingested": 1,
                    "nodes_total": 1,
                    "sources_deleted": int(delete_source_blocks),
                }

            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.storage.subprocess.run",
                side_effect=fake_worker,
            ), mock.patch(
                "tools.full_spectrum.env_scatter_fit.archive.ingest_blocks",
                side_effect=fake_ingest,
            ):
                result = run_streaming_campaign(
                    manifest,
                    archive,
                    spool,
                    delete_packed_blocks=True,
                    seal_output=False,
                )

            self.assertEqual(result["blocks_generated"], 1)
            self.assertEqual(result["blocks_packed"], 1)
            self.assertEqual(result["source_blocks_deleted"], 1)
            self.assertFalse(any(spool.rglob("*.ceefit")))


@unittest.skipUnless(H5PY_AVAILABLE, "h5py archive extra is not installed")
class StorageWorkflowTests(unittest.TestCase):
    def _archive(self, root: Path) -> Path:
        block = root / "block.ceefit"
        block.write_bytes(fixture_block())
        archive = root / "canonical.h5"
        ingest_blocks((block,), archive)
        return archive

    def test_seal_reuses_full_scientific_snapshot_for_read_only_consumers(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256
        from tools.full_spectrum.env_scatter_fit.storage import create_archive_seal

        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(Path(temporary))
            seal = create_archive_seal(archive)
            self.assertIsNone(seal["archive_file_sha256"])

            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.archive.verify_archive_materialization",
                side_effect=AssertionError("full verification should be cached"),
            ):
                with h5py.File(archive, "r") as root:
                    self.assertEqual(
                        archive_snapshot_sha256(root),
                        seal["archive_snapshot_sha256"],
                    )

    def _two_node_archive(self, root: Path) -> tuple[Path, Path]:
        from tools.full_spectrum.env_scatter_fit.design import PhysicalNode

        near = root / "near_1.ceefit"
        far = root / "far_1.ceefit"
        near.write_bytes(fixture_block(physical_node=PhysicalNode(100.0, 1.0, 0.0, 0.0)))
        far.write_bytes(fixture_block(physical_node=PhysicalNode(100.0, 2.0, 0.0, 0.0)))
        archive = root / "canonical.h5"
        ingest_blocks((near, far), archive)
        continuation = root / "near_2.ceefit"
        continuation.write_bytes(
            fixture_block(block_id=2, physical_node=PhysicalNode(100.0, 1.0, 0.0, 0.0))
        )
        return archive, continuation

    def test_incremental_seal_reverifies_only_changed_nodes(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit import archive as archive_module
        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256
        from tools.full_spectrum.env_scatter_fit.storage import (
            FULL_VERIFICATION,
            INCREMENTAL_VERIFICATION,
            cached_archive_snapshot,
            create_archive_seal,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, continuation = self._two_node_archive(root)
            first = create_archive_seal(archive)
            self.assertEqual(first["verification"], FULL_VERIFICATION)
            self.assertEqual(first["nodes_verified"], 2)

            ingest_blocks((continuation,), archive)
            verified = []
            original = archive_module._verify_materialized_node

            def counting(root_group, index, node_hex, numpy):
                verified.append(node_hex)
                return original(root_group, index, node_hex, numpy)

            with mock.patch.object(
                archive_module, "_verify_materialized_node", side_effect=counting
            ):
                incremental = create_archive_seal(archive, incremental=True)
            self.assertEqual(incremental["verification"], INCREMENTAL_VERIFICATION)
            self.assertEqual(incremental["nodes_verified"], 1)
            self.assertEqual(len(verified), 1)
            with h5py.File(archive, "r") as handle:
                # The same identity a full closure assigns.
                self.assertEqual(
                    archive_snapshot_sha256(handle, use_cached=False),
                    incremental["archive_snapshot_sha256"],
                )
                self.assertEqual(
                    cached_archive_snapshot(handle), incremental["archive_snapshot_sha256"]
                )
                # Publication does not rest on it.
                self.assertIsNone(
                    cached_archive_snapshot(handle, require_full_verification=True)
                )

    def test_incremental_seal_still_catches_a_changed_trusted_row(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.storage import create_archive_seal

        with tempfile.TemporaryDirectory() as temporary:
            archive, _ = self._two_node_archive(Path(temporary))
            create_archive_seal(archive)
            with h5py.File(archive, "r+") as handle:
                handle["nodes/wall_seconds"][1] = handle["nodes/wall_seconds"][1] + 1.0
            with self.assertRaisesRegex(ValueError, "does not close to raw blocks"):
                create_archive_seal(archive, incremental=True)

    def test_archive_mutation_invalidates_default_seal(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.storage import (
            create_archive_seal,
            verify_archive_seal,
        )

        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(Path(temporary))
            create_archive_seal(archive)
            with h5py.File(archive, "r+") as root:
                root.attrs["test_mutation"] = 1
            with self.assertRaisesRegex(ValueError, "stat identity"):
                verify_archive_seal(archive)

    def test_fit_cache_is_exact_block_free_and_source_bound(self):
        import h5py
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256
        from tools.full_spectrum.env_scatter_fit.storage import (
            FIT_CACHE_SCHEMA,
            create_archive_seal,
            create_fit_cache,
            verify_fit_cache,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = self._archive(root)
            canonical_seal = create_archive_seal(archive)
            cache = root / "fit-cache.h5"
            result = create_fit_cache(archive, cache)

            self.assertLess(cache.stat().st_size, archive.stat().st_size)
            self.assertIsNone(result["cache_file_sha256"])
            self.assertEqual(verify_fit_cache(cache)["schema"], result["schema"])
            with h5py.File(archive, "r") as source, h5py.File(cache, "r") as fitted:
                self.assertEqual(fitted.attrs["schema"], FIT_CACHE_SCHEMA)
                self.assertNotIn("blocks", fitted)
                np.testing.assert_array_equal(
                    source["continuum/total/S1"][...],
                    fitted["continuum/total/S1"][...],
                )
                self.assertEqual(
                    archive_snapshot_sha256(fitted),
                    canonical_seal["archive_snapshot_sha256"],
                )

            with h5py.File(cache, "r+") as fitted:
                fitted.attrs["test_mutation"] = 1
            with self.assertRaisesRegex(ValueError, "stat identity"):
                verify_fit_cache(cache)

    def test_ingest_invalidates_a_prior_default_seal(self):
        from tools.full_spectrum.env_scatter_fit.storage import (
            create_archive_seal,
            default_archive_seal_path,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = self._archive(root)
            create_archive_seal(archive)
            seal = default_archive_seal_path(archive)
            self.assertTrue(seal.exists())
            continuation = root / "continuation.ceefit"
            continuation.write_bytes(fixture_block(block_id=2))
            ingest_blocks((continuation,), archive)
            self.assertFalse(seal.exists())

    def test_resume_inventory_is_incremental_and_avoids_raw_group_walk(self):
        from tools.full_spectrum.env_scatter_fit.storage import (
            _archive_inventory,
            cached_archive_inventory,
            default_archive_inventory_path,
        )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = self._archive(root)
            initial = cached_archive_inventory(archive)
            self.assertIsNotNone(initial)
            self.assertEqual(len(initial), 1)
            continuation = root / "continuation.ceefit"
            continuation.write_bytes(fixture_block(block_id=2))
            ingest_blocks((continuation,), archive)
            self.assertTrue(default_archive_inventory_path(archive).is_file())
            self.assertEqual(len(cached_archive_inventory(archive)), 2)

            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.archive.require_hdf5",
                side_effect=AssertionError("valid resume index should avoid HDF5 walk"),
            ):
                self.assertEqual(len(_archive_inventory(archive)), 2)

    def test_seal_records_canonical_input_configuration(self):
        from tools.full_spectrum.env_scatter_fit.storage import create_archive_seal

        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(Path(temporary))
            record = create_archive_seal(archive)
            configuration = json.loads(record["canonical_config_json"])
            self.assertEqual(configuration["block_histories"], 262_144)
