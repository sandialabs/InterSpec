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

"""Folding blocks beyond RETAINED_RAW_BLOCKS into a per-node state."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

from tools.full_spectrum.env_scatter_fit.archive import (
    RETAINED_RAW_BLOCKS,
    _block_record,
    _raw_record_sha256,
    archived_block_inventory,
    ingest_blocks,
    verify_archive_materialization,
)
from tools.full_spectrum.env_scatter_fit.block import read_fit_block
from tools.full_spectrum.env_scatter_fit.design import PhysicalNode
from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block


H5PY_AVAILABLE = importlib.util.find_spec("h5py") is not None
NEAR = PhysicalNode(100.0, 1.0, 0.0, 0.0)
FAR = PhysicalNode(100.0, 2.0, 0.0, 0.0)
# Widely separated magnitudes make the compensated sums order-sensitive, so a
# fold resumed in the wrong state or order would not reproduce the rows.
SCALES = (1.0e8, 1.0, 1.0e-4, 3.0, 7.0e-3, 2.0e5)


def _write_blocks(directory: Path, node: PhysicalNode, label: str, block_ids) -> dict[int, Path]:
    paths = {}
    for block_id in block_ids:
        path = directory / f"{label}_{block_id}.ceefit"
        path.write_bytes(
            fixture_block(
                block_id=block_id,
                physical_node=node,
                scale=SCALES[block_id % len(SCALES)],
            )
        )
        paths[block_id] = path
    return paths


def _materialized(path: Path) -> dict:
    import h5py
    import numpy as np

    values = {}
    with h5py.File(path, "r") as root:
        def collect(name, value):
            if not isinstance(value, h5py.Dataset) or name.split("/", 1)[0] in {"blocks", "folds"}:
                return
            if h5py.check_string_dtype(value.dtype) is not None:
                values[name] = tuple(value.asstr()[...].reshape(-1).tolist())
            else:
                values[name] = np.asarray(value[...]).copy()

        root.visititems(collect)
    return values


@unittest.skipUnless(H5PY_AVAILABLE, "h5py archive extra is not installed")
class ArchiveFoldTests(unittest.TestCase):
    def assertSameScience(self, first: Path, second: Path) -> None:
        import h5py
        import numpy as np

        from tools.full_spectrum.env_scatter_fit.publication import archive_snapshot_sha256

        left, right = _materialized(first), _materialized(second)
        self.assertEqual(left.keys(), right.keys())
        for name in left:
            if isinstance(left[name], tuple):
                self.assertEqual(left[name], right[name], name)
            else:
                self.assertEqual(left[name].dtype, right[name].dtype, name)
                self.assertEqual(left[name].tobytes(), right[name].tobytes(), name)
        with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
            self.assertEqual(archive_snapshot_sha256(a), archive_snapshot_sha256(b))
            np.testing.assert_equal(
                {node: {block: (entry["payload_sha256"], entry["scientific_record_sha256"])
                        for block, entry in blocks.items()}
                 for node, blocks in archived_block_inventory(a).items()},
                {node: {block: (entry["payload_sha256"], entry["scientific_record_sha256"])
                        for block, entry in blocks.items()}
                 for node, blocks in archived_block_inventory(b).items()},
            )

    def test_in_memory_record_seal_equals_the_stored_record_seal(self):
        import h5py
        import numpy as np

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            variants = (
                {},
                {"block_id": 2, "scale": 1.0e-4},
                {"block_id": 3, "capped": True, "status": "capped", "failure_reason": "maximum vertex cap reached"},
                {"block_id": 4, "physical_node": PhysicalNode(1500.0, 1.0, 0.0, 0.0)},
            )
            for index, kwargs in enumerate(variants):
                path = root / f"block_{index}.ceefit"
                path.write_bytes(fixture_block(**kwargs))
                archive = root / f"archive_{index}.h5"
                ingest_blocks((path,), archive, fold_blocks=False)
                record = _block_record(read_fit_block(path), np)
                with h5py.File(archive, "r") as stored:
                    node_hex = f"{int(stored['nodes/id'][0]):016x}"
                    group = stored[f"blocks/{node_hex}/{kwargs.get('block_id', 1)}"]
                    self.assertEqual(
                        record.attrs[""]["scientific_record_sha256"],
                        _raw_record_sha256(group, np),
                    )
                    self.assertEqual(
                        group.attrs["scientific_record_sha256"],
                        record.attrs[""]["scientific_record_sha256"],
                    )

    def test_folded_archive_reproduces_raw_archive_bit_for_bit(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 11))
            far = _write_blocks(root, FAR, "far", range(1, 6))
            raw = root / "raw.h5"
            folded = root / "folded.h5"
            ingest_blocks([*near.values(), *far.values()], raw, fold_blocks=False)
            # Rounds: the baseline, then continuations that each fold on arrival.
            ingest_blocks([near[i] for i in range(1, 5)] + [far[i] for i in range(1, 5)], folded)
            ingest_blocks([near[i] for i in range(5, 8)], folded)
            ingest_blocks([far[5]] + [near[i] for i in range(8, 11)], folded)

            self.assertSameScience(raw, folded)
            with h5py.File(folded, "r") as archive:
                near_hex = f"{NEAR.stable_id(archive['manifest'].attrs['config_sha256']):016x}"
                self.assertEqual(
                    sorted(int(name) for name in archive[f"blocks/{near_hex}"]),
                    list(range(1, RETAINED_RAW_BLOCKS + 1)),
                )
                self.assertEqual(archive[f"folds/{near_hex}"].attrs["current"], "10")
                self.assertEqual(list(archive[f"folds/{near_hex}"]), ["10"])
                self.assertEqual(
                    sorted(archived_block_inventory(archive)[near_hex]), list(range(1, 11))
                )
                strategy = archive.id.get_create_plist().get_file_space_strategy()
                self.assertTrue(strategy[1], "persistent free-space tracking is off")
                verify_archive_materialization(archive)

    def test_a_new_node_folds_within_one_batch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 9))
            raw = root / "raw.h5"
            folded = root / "folded.h5"
            ingest_blocks(list(near.values()), raw, fold_blocks=False)
            result = ingest_blocks(list(near.values()), folded)
            self.assertEqual(result["blocks_folded"], 8 - RETAINED_RAW_BLOCKS)
            self.assertSameScience(raw, folded)

    def test_blocks_beyond_a_gap_stay_raw_until_the_gap_is_filled(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 9))
            raw = root / "raw.h5"
            folded = root / "folded.h5"
            ingest_blocks(list(near.values()), raw, fold_blocks=False)
            ingest_blocks([near[i] for i in (1, 2, 3, 4, 5)], folded)
            result = ingest_blocks([near[7], near[8]], folded)
            self.assertEqual(result["blocks_folded"], 0)
            with h5py.File(folded, "r") as archive:
                node_hex = next(iter(archive["folds"]))
                self.assertEqual(sorted(int(n) for n in archive[f"blocks/{node_hex}"]), [1, 2, 3, 4, 7, 8])
                verify_archive_materialization(archive)
            result = ingest_blocks([near[6]], folded)
            self.assertEqual(result["blocks_folded"], 1)
            with h5py.File(folded, "r") as archive:
                self.assertEqual(sorted(int(n) for n in archive[f"blocks/{node_hex}"]), [1, 2, 3, 4])
                self.assertEqual(archive[f"folds/{node_hex}"].attrs["current"], "8")
            self.assertSameScience(raw, folded)

    def test_folded_blocks_are_rejected_as_replays(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 7))
            archive = root / "folded.h5"
            ingest_blocks(list(near.values()), archive)
            with self.assertRaisesRegex(ValueError, "already archived"):
                ingest_blocks([near[6]], archive)

    def test_source_blocks_and_the_old_generation_go_only_after_verification(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 9))
            archive = root / "folded.h5"
            ingest_blocks([near[i] for i in range(1, 7)], archive, delete_source_blocks=True)
            self.assertFalse(any(near[i].exists() for i in range(1, 7)))
            ingest_blocks([near[7], near[8]], archive, delete_source_blocks=True)
            with h5py.File(archive, "r") as stored:
                node_hex = next(iter(stored["folds"]))
                self.assertEqual(list(stored[f"folds/{node_hex}"]), ["8"])

            near = _write_blocks(root, NEAR, "again", range(9, 11))
            from unittest import mock

            with mock.patch(
                "tools.full_spectrum.env_scatter_fit.archive._verify_materialized_node",
                side_effect=ValueError("injected post-write verification failure"),
            ):
                with self.assertRaisesRegex(ValueError, "injected"):
                    ingest_blocks([near[9], near[10]], archive, delete_source_blocks=True)
            self.assertTrue(near[9].exists() and near[10].exists())
            with h5py.File(archive, "r") as stored:
                # The new generation is current; the one it replaced is kept.
                self.assertEqual(stored[f"folds/{node_hex}"].attrs["current"], "10")
                self.assertEqual(sorted(stored[f"folds/{node_hex}"]), ["10", "8"])
                verify_archive_materialization(stored)

    def test_fold_state_and_retained_raw_corruption_are_detected(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 7))
            archive = root / "folded.h5"
            ingest_blocks(list(near.values()), archive)
            with h5py.File(archive, "a") as stored:
                node_hex = next(iter(stored["folds"]))
                total = stored[f"folds/{node_hex}/6/float/moments/continuum_total/S1/total"]
                values = total[...]
                values[0] += 1.0
                total[...] = values
            with h5py.File(archive, "r") as stored:
                with self.assertRaisesRegex(ValueError, "fold state of node .* failed its checksum"):
                    verify_archive_materialization(stored)

            archive = root / "missing_raw.h5"
            ingest_blocks(list(near.values()), archive)
            with h5py.File(archive, "a") as stored:
                del stored[f"blocks/{node_hex}/2"]
            with h5py.File(archive, "r") as stored:
                with self.assertRaisesRegex(ValueError, "lacks its retained raw blocks"):
                    verify_archive_materialization(stored)

    def test_rebuild_refuses_to_drop_histories(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 4))
            archive = root / "raw.h5"
            ingest_blocks([near[1], near[2]], archive)
            with h5py.File(archive, "a") as stored:
                node_hex = next(iter(stored["blocks"]))
                del stored[f"blocks/{node_hex}/1"]
                del stored[f"blocks/{node_hex}/2"]
            with self.assertRaisesRegex(ValueError, "would drop histories"):
                ingest_blocks([near[3]], archive)

    def test_calibration_pairs_only_the_raw_records(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.calibration import _pair_inventory

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 11))
            archive = root / "folded.h5"
            ingest_blocks(list(near.values()), archive)
            with h5py.File(archive, "r") as stored:
                pairs, unused = _pair_inventory(stored)
            self.assertEqual([pair[1:] for pair in pairs], [("1", "2"), ("3", "4")])
            self.assertEqual(unused, [])

    def test_a_damaged_fold_state_is_never_extended(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 8))
            archive = root / "folded.h5"
            ingest_blocks([near[i] for i in range(1, 7)], archive)
            with h5py.File(archive, "a") as stored:
                node_hex = next(iter(stored["folds"]))
                total = stored[f"folds/{node_hex}/6/float/moments/continuum_total/S1/total"]
                values = total[...]
                values[0] += 1.0
                total[...] = values
            with self.assertRaisesRegex(ValueError, "failed its checksum"):
                ingest_blocks([near[7]], archive, delete_source_blocks=True)
            self.assertTrue(near[7].exists())

    def test_failed_verifications_never_drop_the_last_verified_generation(self):
        import h5py
        from unittest import mock

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 12))
            archive = root / "folded.h5"
            ingest_blocks([near[i] for i in range(1, 9)], archive)
            for batch in ((9,), (10,)):
                with mock.patch(
                    "tools.full_spectrum.env_scatter_fit.archive._verify_materialized_node",
                    side_effect=ValueError("injected post-write verification failure"),
                ):
                    with self.assertRaisesRegex(ValueError, "injected"):
                        ingest_blocks([near[i] for i in batch], archive)
            with h5py.File(archive, "r") as stored:
                node_hex = next(iter(stored["folds"]))
                self.assertEqual(sorted(stored[f"folds/{node_hex}"], key=int), ["8", "9", "10"])
            ingest_blocks([near[11]], archive)
            with h5py.File(archive, "r") as stored:
                self.assertEqual(list(stored[f"folds/{node_hex}"]), ["11"])
                verify_archive_materialization(stored)

    def test_an_interrupted_first_fold_write_is_recoverable(self):
        import h5py
        from unittest import mock

        from tools.full_spectrum.env_scatter_fit import archive as archive_module

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 7))
            raw = root / "raw.h5"
            archive = root / "folded.h5"
            ingest_blocks(list(near.values()), raw, fold_blocks=False)
            ingest_blocks([near[i] for i in range(1, 5)], archive)
            write = archive_module._Record.write

            def interrupted(self, parent, name):
                if "folds" in parent.name:
                    parent.create_group(name)
                    raise KeyboardInterrupt("interrupted generation write")
                return write(self, parent, name)

            with mock.patch.object(archive_module._Record, "write", interrupted):
                with self.assertRaises(KeyboardInterrupt):
                    ingest_blocks([near[5], near[6]], archive)
            with h5py.File(archive, "r") as stored:
                self.assertEqual(sorted(archived_block_inventory(stored).popitem()[1]), [1, 2, 3, 4])
            ingest_blocks([near[5], near[6]], archive)
            self.assertSameScience(raw, archive)

    def test_blocks_ahead_of_the_retained_ones_fold_once_those_arrive(self):
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 9))
            raw = root / "raw.h5"
            archive = root / "folded.h5"
            ingest_blocks(list(near.values()), raw, fold_blocks=False)
            ingest_blocks([near[i] for i in (8, 5, 7, 6)], archive)
            ingest_blocks([near[i] for i in (3, 1, 4, 2)], archive)
            with h5py.File(archive, "r") as stored:
                node_hex = next(iter(stored["folds"]))
                self.assertEqual(sorted(int(n) for n in stored[f"blocks/{node_hex}"]), [1, 2, 3, 4])
                self.assertEqual(stored[f"folds/{node_hex}"].attrs["current"], "8")
            self.assertSameScience(raw, archive)

    def test_a_folded_node_calibrates_on_its_retained_blocks_only(self):
        import h5py

        from tools.full_spectrum.env_scatter_fit.calibration import _pair_inventory

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            near = _write_blocks(root, NEAR, "near", range(1, 9))
            archive = root / "folded.h5"
            ingest_blocks([near[i] for i in range(1, 6)], archive)
            ingest_blocks([near[7], near[8]], archive)
            with h5py.File(archive, "r") as stored:
                node_hex = next(iter(stored["folds"]))
                self.assertEqual(sorted(int(n) for n in stored[f"blocks/{node_hex}"]), [1, 2, 3, 4, 7, 8])
                pairs, unused = _pair_inventory(stored)
            self.assertEqual([pair[1:] for pair in pairs], [("1", "2"), ("3", "4")])
            self.assertEqual(unused, [])


if __name__ == "__main__":
    unittest.main()
