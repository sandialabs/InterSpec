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

"""Lossless storage orchestration for fit-first Monte Carlo campaigns.

The canonical HDF5 archive remains the scientific record.  Worker blocks and
fit caches are deliberately disposable: a block is removed only after its
normalized HDF5 record and affected node aggregates have been verified, while
a fit cache is accepted only when it is bound to a sealed canonical snapshot.
"""

from __future__ import annotations

from concurrent.futures import Future, ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Iterable, Mapping, Sequence


ARCHIVE_SEAL_SCHEMA = "CEEFIT-ARCHIVE-SEAL-v1"
ARCHIVE_INVENTORY_SCHEMA = "CEEFIT-ARCHIVE-INVENTORY-v1"
ARCHIVE_VERIFIED_NODES_SCHEMA = "CEEFIT-ARCHIVE-VERIFIED-NODES-v1"
FULL_VERIFICATION = "full-raw-record-and-materialized-science-closure-v1"
INCREMENTAL_VERIFICATION = "incremental-raw-record-and-materialized-science-closure-v1"
FIT_CACHE_SCHEMA = "CEEFIT-FIT-CACHE-v1"
FIT_CACHE_SEAL_SCHEMA = "CEEFIT-FIT-CACHE-SEAL-v1"
_FIT_CACHE_GROUPS = (
    "manifest",
    "axes",
    "components",
    "nodes",
    "continuum",
    "lines",
    "diagnostics",
    "replicas",
    "stopping",
)


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _fsync_file(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        # Some removable filesystems do not expose directory fsync.  The data
        # file itself is still explicitly flushed before transactional cleanup.
        pass
    finally:
        os.close(descriptor)


def _atomic_json(path: Path, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
        _fsync_directory(path.parent)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _stat_record(path: Path) -> dict[str, int]:
    status = path.stat()
    return {
        "size_bytes": int(status.st_size),
        "mtime_ns": int(status.st_mtime_ns),
        "device": int(status.st_dev),
        "inode": int(status.st_ino),
    }


def _stat_matches(path: Path, record: Mapping[str, Any]) -> bool:
    try:
        current = _stat_record(path)
    except OSError:
        return False
    return all(record.get(name) == value for name, value in current.items())


def default_archive_seal_path(archive_path: str | Path) -> Path:
    archive = Path(archive_path)
    return archive.with_name(archive.name + ".seal.json")


def default_archive_inventory_path(archive_path: str | Path) -> Path:
    archive = Path(archive_path)
    return archive.with_name(archive.name + ".inventory.json")


def default_archive_verified_nodes_path(archive_path: str | Path) -> Path:
    archive = Path(archive_path)
    return archive.with_name(archive.name + ".verified-nodes.json")


def _trusted_archive_nodes(archive: Path, config_sha256: str) -> dict[str, str] | None:
    """The node closure identities the last seal verified, or ``None``.

    The ledger survives archive mutation on purpose (the seal does not): it
    names nodes whose raw-block closure was verified, keyed by an identity
    that changes whenever the node's blocks or materialized row change.
    """

    try:
        record = json.loads(
            default_archive_verified_nodes_path(archive).read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    nodes = record.get("nodes")
    if (
        record.get("schema") != ARCHIVE_VERIFIED_NODES_SCHEMA
        or record.get("config_sha256") != config_sha256
        or not isinstance(nodes, dict)
        or not all(_is_sha256(value) for value in nodes.values())
    ):
        return None
    return {str(key): str(value) for key, value in nodes.items()}


def default_fit_cache_seal_path(cache_path: str | Path) -> Path:
    cache = Path(cache_path)
    return cache.with_name(cache.name + ".seal.json")


def cached_archive_inventory(
    archive_path: str | Path,
) -> set[tuple[int, int]] | None:
    """Return a stat-bound resume inventory, or ``None`` when it is stale."""

    archive = Path(archive_path)
    if not archive.is_file():
        return set()
    sidecar = default_archive_inventory_path(archive)
    try:
        record = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if record.get("schema") != ARCHIVE_INVENTORY_SCHEMA or not _stat_matches(
        archive, record.get("archive_stat", {})
    ):
        return None
    raw_pairs = record.get("blocks")
    if not isinstance(raw_pairs, list):
        return None
    result: set[tuple[int, int]] = set()
    try:
        for pair in raw_pairs:
            if (
                not isinstance(pair, list)
                or len(pair) != 2
                or not isinstance(pair[0], str)
                or len(pair[0]) != 16
                or type(pair[1]) is not int
                or pair[1] <= 0
            ):
                return None
            node_id = int(pair[0], 16)
            if f"{node_id:016x}" != pair[0]:
                return None
            result.add((node_id, pair[1]))
    except ValueError:
        return None
    if len(result) != len(raw_pairs):
        return None
    return result


def write_archive_inventory(
    archive_path: str | Path,
    blocks: Iterable[tuple[int, int]],
    config_sha256: str,
) -> Path:
    """Atomically write a disposable fast-resume index for one archive."""

    archive = Path(archive_path)
    if not archive.is_file() or not _is_sha256(config_sha256):
        raise ValueError("cannot index a missing archive or invalid configuration")
    normalized = sorted(
        {(int(node_id), int(block_id)) for node_id, block_id in blocks}
    )
    if any(node_id < 0 or block_id <= 0 for node_id, block_id in normalized):
        raise ValueError("archive inventory contains an invalid block identity")
    destination = default_archive_inventory_path(archive)
    _atomic_json(
        destination,
        {
            "schema": ARCHIVE_INVENTORY_SCHEMA,
            "archive_name": archive.name,
            "archive_stat": _stat_record(archive),
            "config_sha256": config_sha256,
            "blocks": [
                [f"{node_id:016x}", block_id]
                for node_id, block_id in normalized
            ],
        },
    )
    return destination


def invalidate_archive_seal(archive_path: str | Path) -> None:
    """Remove a derived default seal after the canonical archive changes."""

    seal = default_archive_seal_path(archive_path)
    if seal.exists():
        seal.unlink()
        _fsync_directory(seal.parent)


def cached_archive_snapshot(
    root: Any, *, require_full_verification: bool = False
) -> str | None:
    """Return a stat-bound cached canonical identity, or ``None``.

    A stale/missing seal never blesses data: callers fall back to the complete
    scientific closure calculation.  Only read-only HDF5 handles may use this
    shortcut, so mutation tests and repair tools always exercise full closure.
    With ``require_full_verification`` an incremental seal does not count.
    """

    if str(root.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
        return None
    if str(getattr(root, "mode", "")) != "r":
        return None
    filename = getattr(root, "filename", None)
    if not filename:
        return None
    archive = Path(filename)
    seal_path = default_archive_seal_path(archive)
    if not seal_path.is_file():
        return None
    try:
        record = json.loads(seal_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if record.get("schema") != ARCHIVE_SEAL_SCHEMA:
        return None
    if require_full_verification and record.get("verification") != FULL_VERIFICATION:
        return None
    if not _stat_matches(archive, record.get("archive_stat", {})):
        return None
    snapshot = record.get("archive_snapshot_sha256")
    config = record.get("config_sha256")
    if not _is_sha256(snapshot) or config != str(
        root["manifest"].attrs.get("config_sha256", "")
    ):
        return None
    return str(snapshot)


def create_archive_seal(
    archive_path: str | Path,
    output_path: str | Path | None = None,
    *,
    full_file_hash: bool = False,
    incremental: bool = False,
) -> dict[str, Any]:
    """Verify an archive and persist its reusable local identity.

    A full seal verifies every node's raw-block closure. An ``incremental`` one
    re-verifies only the nodes whose closure identity changed since the last
    seal recorded them (in the ``.verified-nodes.json`` ledger every seal
    writes), and runs every archive-wide check; it yields the same snapshot
    identity, but publication requires a full seal. Without a usable ledger
    an incremental seal is a full one.
    """

    from .archive import archived_block_inventory, require_hdf5
    from .publication import archive_snapshot_sha256

    h5py, _ = require_hdf5()
    archive = Path(archive_path).resolve()
    if not archive.is_file():
        raise ValueError(f"archive does not exist: {archive}")
    identities: dict[str, str] = {}
    with h5py.File(archive, "r") as root:
        if str(root.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
            raise ValueError("only canonical CEEFIT HDF5 archives may be sealed")
        config = str(root["manifest"].attrs["config_sha256"])
        trusted = _trusted_archive_nodes(archive, config) if incremental else None
        snapshot = archive_snapshot_sha256(
            root,
            use_cached=False,
            trusted_nodes=trusted,
            identities_out=identities,
        )
        canonical_config_json = str(
            root["manifest"].attrs.get("canonical_config_json", "")
        )
        node_count = len(root["nodes/id"])
        block_count = sum(len(blocks) for blocks in archived_block_inventory(root).values())
    _fsync_file(archive)
    reverified = sum(
        1
        for node, identity in identities.items()
        if trusted is None or trusted.get(node) != identity
    )
    record = {
        "schema": ARCHIVE_SEAL_SCHEMA,
        "archive_schema": "CEEFIT-HDF5-v1",
        "archive_name": archive.name,
        "archive_stat": _stat_record(archive),
        "archive_file_sha256": _file_sha256(archive) if full_file_hash else None,
        "archive_snapshot_sha256": snapshot,
        "config_sha256": config,
        "canonical_config_json": canonical_config_json,
        "node_count": node_count,
        "block_count": block_count,
        "verification": FULL_VERIFICATION if trusted is None else INCREMENTAL_VERIFICATION,
        "nodes_verified": reverified,
    }
    destination = (
        default_archive_seal_path(archive)
        if output_path is None
        else Path(output_path).resolve()
    )
    _atomic_json(destination, record)
    if output_path is None:
        _atomic_json(
            default_archive_verified_nodes_path(archive),
            {
                "schema": ARCHIVE_VERIFIED_NODES_SCHEMA,
                "archive_name": archive.name,
                "config_sha256": config,
                "nodes": identities,
            },
        )
    return {**record, "output": str(destination)}


def verify_archive_seal(
    archive_path: str | Path,
    seal_path: str | Path | None = None,
    *,
    full_file_hash: bool = False,
) -> dict[str, Any]:
    """Validate a cached seal without decompressing scientific datasets."""

    archive = Path(archive_path).resolve()
    source = (
        default_archive_seal_path(archive)
        if seal_path is None
        else Path(seal_path)
    )
    record = json.loads(source.read_text(encoding="utf-8"))
    if record.get("schema") != ARCHIVE_SEAL_SCHEMA:
        raise ValueError("unsupported archive seal schema")
    if not _stat_matches(archive, record.get("archive_stat", {})):
        raise ValueError("archive stat identity differs from its seal")
    from .archive import require_hdf5

    h5py, _ = require_hdf5()
    with h5py.File(archive, "r") as root:
        if str(root.attrs.get("schema", "")) != record.get("archive_schema"):
            raise ValueError("archive schema differs from its seal")
        if str(root["manifest"].attrs.get("config_sha256", "")) != record.get(
            "config_sha256"
        ):
            raise ValueError("archive configuration differs from its seal")
    if full_file_hash:
        expected = record.get("archive_file_sha256")
        if not _is_sha256(expected):
            raise ValueError(
                "archive seal does not contain an optional whole-file hash"
            )
        if _file_sha256(archive) != expected:
            raise ValueError("archive file checksum differs from its seal")
    return record


def create_fit_cache(
    archive_path: str | Path,
    output_path: str | Path,
    *,
    full_file_hash: bool = False,
) -> dict[str, Any]:
    """Create an aggregate-only, disposable cache for model fitting.

    The cache preserves every aggregate dataset used by fitting and validation,
    including full-resolution component/total replica maps.  It omits only the
    embedded immutable block records used for calibration and publication
    closure, reducing the measured soil cache from about 13 GiB to 2.9 GiB.
    """

    from .archive import require_hdf5
    from .publication import archive_snapshot_sha256

    h5py, _ = require_hdf5()
    archive = Path(archive_path).resolve()
    output = Path(output_path).resolve()
    if archive == output:
        raise ValueError("fit cache output must differ from the canonical archive")
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.cache-", suffix=".tmp", dir=output.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with h5py.File(archive, "r") as source:
            if str(source.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
                raise ValueError("fit cache source is not a canonical CEEFIT archive")
            snapshot = archive_snapshot_sha256(source)
            config = str(source["manifest"].attrs["config_sha256"])
            with h5py.File(temporary, "w") as destination:
                for name, value in source.attrs.items():
                    destination.attrs[name] = value
                destination.attrs["schema"] = FIT_CACHE_SCHEMA
                destination.attrs["source_archive_snapshot_sha256"] = snapshot
                destination.attrs["source_archive_config_sha256"] = config
                destination.attrs["source_archive_path"] = str(archive)
                for name in _FIT_CACHE_GROUPS:
                    if name not in source:
                        raise ValueError(f"canonical archive lacks required group {name!r}")
                    source.copy(name, destination, name=name)
                destination.flush()
        _fsync_file(temporary)
        os.replace(temporary, output)
        _fsync_directory(output.parent)
        temporary = None
        record = {
            "schema": FIT_CACHE_SEAL_SCHEMA,
            "fit_cache_schema": FIT_CACHE_SCHEMA,
            "cache_name": output.name,
            "cache_stat": _stat_record(output),
            "cache_file_sha256": _file_sha256(output) if full_file_hash else None,
            "source_archive_path": str(archive),
            "source_archive_snapshot_sha256": snapshot,
            "config_sha256": config,
            "groups": list(_FIT_CACHE_GROUPS),
        }
        seal = default_fit_cache_seal_path(output)
        _atomic_json(seal, record)
        return {**record, "output": str(output), "seal": str(seal)}
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def fit_cache_source_snapshot(root: Any) -> str:
    """Validate a fit cache's stat-bound seal and return its source identity."""

    if str(root.attrs.get("schema", "")) != FIT_CACHE_SCHEMA:
        raise ValueError("input is not a CEEFIT fit cache")
    if str(getattr(root, "mode", "")) != "r":
        raise ValueError("fit caches must be opened read-only")
    filename = getattr(root, "filename", None)
    if not filename:
        raise ValueError("fit cache has no filesystem identity")
    cache = Path(filename)
    seal_path = default_fit_cache_seal_path(cache)
    try:
        record = json.loads(seal_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("fit cache lacks a valid seal sidecar") from error
    if record.get("schema") != FIT_CACHE_SEAL_SCHEMA:
        raise ValueError("unsupported fit-cache seal schema")
    if not _stat_matches(cache, record.get("cache_stat", {})):
        raise ValueError("fit cache stat identity differs from its seal")
    snapshot = record.get("source_archive_snapshot_sha256")
    config = str(root["manifest"].attrs.get("config_sha256", ""))
    if (
        not _is_sha256(snapshot)
        or record.get("config_sha256") != config
        or str(root.attrs.get("source_archive_snapshot_sha256", "")) != snapshot
    ):
        raise ValueError("fit cache source identity is malformed or inconsistent")
    return str(snapshot)


def verify_fit_cache(
    cache_path: str | Path, *, full_file_hash: bool = False
) -> dict[str, Any]:
    cache = Path(cache_path).resolve()
    seal = default_fit_cache_seal_path(cache)
    record = json.loads(seal.read_text(encoding="utf-8"))
    if record.get("schema") != FIT_CACHE_SEAL_SCHEMA:
        raise ValueError("unsupported fit-cache seal schema")
    if not _stat_matches(cache, record.get("cache_stat", {})):
        raise ValueError("fit cache stat identity differs from its seal")
    from .archive import require_hdf5

    h5py, _ = require_hdf5()
    with h5py.File(cache, "r") as root:
        snapshot = fit_cache_source_snapshot(root)
        if snapshot != record.get("source_archive_snapshot_sha256"):
            raise ValueError("fit cache source snapshot differs from its seal")
    if full_file_hash:
        expected = record.get("cache_file_sha256")
        if not _is_sha256(expected):
            raise ValueError(
                "fit-cache seal does not contain an optional whole-file hash"
            )
        if _file_sha256(cache) != expected:
            raise ValueError("fit cache file checksum differs from its seal")
    return record


def _replace_fit_block_output(argv: Sequence[str], output: Path) -> list[str]:
    result = list(argv)
    try:
        index = result.index("--fit-block") + 1
    except ValueError as error:
        raise ValueError("worker command lacks --fit-block") from error
    result[index] = str(output)
    return result


def _archive_inventory(archive_path: Path) -> set[tuple[int, int]]:
    if not archive_path.exists():
        return set()
    cached = cached_archive_inventory(archive_path)
    if cached is not None:
        return cached
    from .archive import archived_block_inventory, require_hdf5

    h5py, _ = require_hdf5()
    result: set[tuple[int, int]] = set()
    with h5py.File(archive_path, "r") as root:
        config = str(root["manifest"].attrs.get("config_sha256", ""))
        for node_hex, blocks in archived_block_inventory(root).items():
            node_id = int(node_hex, 16)
            result.update((node_id, block_id) for block_id in blocks)
    write_archive_inventory(archive_path, result, config)
    return result


# A node's fold state measured at most this many of its own raw block records
# on the field/soil campaign (<= 2.0 at 4 blocks, <= 4.3 at the 488-block cap,
# folding pseudo-blocks resampled from each node's baseline): the expected
# size.  Its ceiling is the state's uncompressed size, which real independent
# blocks (filling cells the baseline never hit) can approach.
FOLD_STATE_TO_RAW_RECORD_ESTIMATE = 5.0
# Per folded block: ledger replica arrays (4 x 32 u8) and its JSON row.
_LEDGER_BYTES_PER_BLOCK = 2048
# HDF5 object headers and attributes a raw record carries beyond its datasets.
_RAW_RECORD_METADATA = 1.03


def _storage_bytes(group: Any) -> int:
    total = 0

    def visit(_: str, value: Any) -> None:
        nonlocal total
        if hasattr(value, "id") and hasattr(value.id, "get_storage_size"):
            total += int(value.id.get_storage_size())

    group.visititems(visit)
    return total


def projected_archive_growth(
    archive_path: str | Path,
    work: Iterable[tuple[int, int]],
    *,
    fit_cache: bool = False,
) -> dict[str, Any]:
    """What archiving ``work`` ((node id, block id) pairs) adds to the archive.

    Blocks 1..RETAINED_RAW_BLOCKS become raw records, priced at the node's own
    mean record size (the archive's mean for a new node).  Later blocks fold
    into the node's state.  ``total_bytes`` is a ceiling: every continued
    node's state grows to its uncompressed size, and while a node is packed
    the generation its new state replaces still exists.  ``expected_bytes``
    uses the measured FOLD_STATE_TO_RAW_RECORD_ESTIMATE instead.  New nodes
    add an aggregate row (and a fit-cache row).  An archive created without
    persistent free-space tracking never reuses the space of a replaced
    generation, so there a continued node's whole new state is new space.
    """

    from .archive import RETAINED_RAW_BLOCKS, _fold_generation, _node_major_paths, require_hdf5

    h5py, _ = require_hdf5()
    by_node: dict[int, list[int]] = {}
    for node_id, block_id in work:
        by_node.setdefault(int(node_id), []).append(int(block_id))
    with h5py.File(Path(archive_path), "r") as root:
        reuses_space = bool(root.id.get_create_plist().get_file_space_strategy()[1])
        raw_bytes: dict[int, float] = {}
        state_bytes: dict[int, int] = {}
        for node_hex, node_group in root["blocks"].items():
            records = len(node_group)
            if records:
                raw_bytes[int(node_hex, 16)] = _storage_bytes(node_group) / records
            generation = _fold_generation(root, node_hex)
            if generation is not None:
                state_bytes[int(node_hex, 16)] = _storage_bytes(generation)
        # A fold state holds at most one float64 (total, compensation) pair per
        # element of the node's materialized row.
        ceiling = 0
        for path in _node_major_paths(root):
            elements = 1
            for size in root[path].shape[1:]:
                elements *= int(size)
            ceiling += 16 * elements
        node_count = len(root["nodes/id"])
        aggregate = sum(
            _storage_bytes(root[name])
            for name in root
            if name not in {"blocks", "folds"} and hasattr(root[name], "visititems")
        )
    mean_raw = sum(raw_bytes.values()) / len(raw_bytes) if raw_bytes else 0.0
    per_node_aggregate = aggregate / node_count if node_count else 0.0
    raw_growth = state_ceiling = state_expected = new_rows = 0.0
    largest_ceiling = largest_expected = 0.0
    for node_id, block_ids in by_node.items():
        record = raw_bytes.get(node_id, mean_raw)
        raw_growth += _RAW_RECORD_METADATA * record * sum(
            1 for value in block_ids if value <= RETAINED_RAW_BLOCKS
        )
        if node_id not in raw_bytes:
            new_rows += 1
        if not any(value > RETAINED_RAW_BLOCKS for value in block_ids):
            continue
        present = 0 if not reuses_space else state_bytes.get(node_id, 0)
        node_ceiling = ceiling + _LEDGER_BYTES_PER_BLOCK * max(block_ids)
        node_expected = FOLD_STATE_TO_RAW_RECORD_ESTIMATE * record
        state_ceiling += max(0.0, node_ceiling - present)
        state_expected += max(0.0, node_expected - present)
        largest_ceiling = max(largest_ceiling, node_ceiling)
        largest_expected = max(largest_expected, node_expected)
    aggregate_growth = new_rows * per_node_aggregate * (2 if fit_cache else 1)
    return {
        "blocks": sum(len(values) for values in by_node.values()),
        "reuses_freed_space": reuses_space,
        "raw_record_bytes": raw_growth,
        "fold_state_bytes": state_ceiling,
        "expected_fold_state_bytes": state_expected,
        "aggregate_bytes": aggregate_growth,
        "superseded_generation_bytes": largest_ceiling,
        "total_bytes": raw_growth + state_ceiling + aggregate_growth + largest_ceiling,
        "expected_bytes": raw_growth + state_expected + aggregate_growth + largest_expected,
        "mean_raw_record_bytes": mean_raw,
    }


def run_streaming_campaign(
    manifest: Mapping[str, Any],
    archive_path: str | Path,
    spool_path: str | Path,
    *,
    roles: Iterable[str] | None = None,
    spool_max_bytes: int = 4 * 1024**3,
    delete_packed_blocks: bool = False,
    seal_output: bool = True,
    incremental_seal: bool = False,
    max_groups: int | None = None,
) -> dict[str, Any]:
    """Run worker groups while one background writer drains a bounded spool."""

    from .archive import ingest_blocks

    if spool_max_bytes <= 0:
        raise ValueError("spool_max_bytes must be positive")
    archive = Path(archive_path).resolve()
    spool = Path(spool_path).resolve()
    if archive == spool or spool in archive.parents:
        raise ValueError("spool must not contain the canonical archive")
    spool.mkdir(parents=True, exist_ok=True)
    selected_roles = None if roles is None else {str(value) for value in roles}
    logical_work = {
        (int(item["node_id"]), int(item["block_id"])): item
        for item in manifest.get("work", ())
    }
    groups = [
        item
        for item in manifest.get("worker_groups", ())
        if selected_roles is None or str(item["role"]) in selected_roles
    ]
    if max_groups is not None:
        if max_groups <= 0:
            raise ValueError("max_groups must be positive")
        groups = groups[:max_groups]
    archived = _archive_inventory(archive)
    pending: list[tuple[Future[dict[str, int]], int]] = []
    pending_bytes = 0
    generated_blocks = 0
    skipped_blocks = 0
    packed_blocks = 0
    deleted_blocks = 0

    def finish_oldest() -> None:
        nonlocal pending_bytes, packed_blocks, deleted_blocks
        future, byte_count = pending.pop(0)
        result = future.result()
        pending_bytes -= byte_count
        packed_blocks += int(result["blocks_ingested"])
        deleted_blocks += int(result.get("sources_deleted", 0))

    with ThreadPoolExecutor(max_workers=1, thread_name_prefix="ceefit-packer") as pool:
        for group in groups:
            node_id = int(group["node_id"])
            role = str(group["role"])
            first = int(group["first_block_id"])
            count = int(group["block_count"])
            block_ids = list(range(first, first + count))
            needed = [value for value in block_ids if (node_id, value) not in archived]
            skipped_blocks += count - len(needed)
            if not needed:
                continue
            node_dir = spool / role / f"{node_id:016x}"
            node_dir.mkdir(parents=True, exist_ok=True)
            paths = {value: node_dir / f"block_{value}.ceefit" for value in needed}
            missing = [value for value, path in paths.items() if not path.is_file()]
            if missing == block_ids:
                template = node_dir / "block_{block_id}.ceefit"
                subprocess.run(
                    _replace_fit_block_output(group["argv"], template), check=True
                )
                generated_blocks += len(missing)
            else:
                for block_id in missing:
                    item = logical_work.get((node_id, block_id))
                    if item is None:
                        raise ValueError(
                            f"manifest lacks logical work for {node_id:016x}/{block_id}"
                        )
                    subprocess.run(
                        _replace_fit_block_output(item["argv"], paths[block_id]),
                        check=True,
                    )
                    generated_blocks += 1
            ordered_paths = [paths[value] for value in needed]
            byte_count = sum(path.stat().st_size for path in ordered_paths)
            future = pool.submit(
                ingest_blocks,
                ordered_paths,
                archive,
                delete_source_blocks=delete_packed_blocks,
            )
            pending.append((future, byte_count))
            pending_bytes += byte_count
            while pending and pending_bytes >= spool_max_bytes:
                finish_oldest()
        while pending:
            finish_oldest()

    seal = (
        create_archive_seal(archive, incremental=incremental_seal)
        if seal_output and archive.exists()
        else None
    )
    return {
        "groups_considered": len(groups),
        "blocks_generated": generated_blocks,
        "blocks_skipped_as_archived": skipped_blocks,
        "blocks_packed": packed_blocks,
        "source_blocks_deleted": deleted_blocks,
        "archive": str(archive),
        "spool": str(spool),
        "archive_seal": None if seal is None else seal["output"],
    }
