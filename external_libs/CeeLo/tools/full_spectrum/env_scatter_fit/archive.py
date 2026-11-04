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

"""Convert validated CEEFIT worker blocks into a chunked HDF5 fit archive."""

from __future__ import annotations

import json
import hashlib
import math
import os
from pathlib import Path
import struct
import tempfile
from typing import Any, Iterable, Mapping, Sequence

from .block import FitBlock, read_fit_block
from .contracts import FORMAL_BLOCK_HISTORIES, RETAINED_RAW_BLOCKS


ARCHIVE_SCHEMA = "CEEFIT-HDF5-v1"
RAW_RECORD_SCHEMA = "CEEFIT-HDF5-BLOCK-RECORD-v1"
FOLD_STATE_SCHEMA = "CEEFIT-HDF5-FOLD-STATE-v1"
# Free space left by deleted records and replaced fold generations is tracked
# in the file and reused; without it HDF5 never reuses space freed in an
# earlier session and the archive grows as if nothing were deleted.
_FILE_SPACE = {"fs_strategy": "fsm", "fs_persist": True, "fs_threshold": 1}
_CANONICAL_LINES = (
    {"name": "elastic_E0", "origin": 0, "post_creation_order": 3},
    {"name": "annihilation_511", "origin": 1, "post_creation_order": 3},
)


def require_hdf5() -> tuple[Any, Any]:
    try:
        import h5py
        import numpy
    except ImportError as error:
        raise RuntimeError(
            "HDF5 conversion requires NumPy and h5py; install the archive extras "
            "with `python -m pip install -r "
            "tools/full_spectrum/env_scatter_fit/requirements-archive.txt`"
        ) from error
    return h5py, numpy


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _normal_hdf_value(value: Any) -> Any:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, bytes):
        return {"bytes_hex": value.hex()}
    if isinstance(value, (tuple, list)):
        return [_normal_hdf_value(item) for item in value]
    if isinstance(value, (str, bool, int, float)) or value is None:
        if isinstance(value, float) and not math.isfinite(value):
            return {"float_hex": value.hex()}
        return value
    if hasattr(value, "item"):
        return _normal_hdf_value(value.item())
    raise ValueError(f"unsupported HDF5 attribute value {type(value).__name__}")


def _hash_frame(digest: Any, label: str, payload: bytes) -> None:
    encoded = label.encode("utf-8")
    digest.update(struct.pack("<QQ", len(encoded), len(payload)))
    digest.update(encoded)
    digest.update(payload)


_RECORD_SEAL_EXCLUDED = frozenset({"source_path", "scientific_record_sha256"})
_FOLD_SEAL_EXCLUDED = frozenset({"state_sha256"})


def _hash_attributes(digest: Any, relative: str, attributes: Mapping[str, Any], excluded: frozenset[str]) -> None:
    normalized = {
        str(name): _normal_hdf_value(value)
        for name, value in attributes.items()
        if str(name) not in excluded
    }
    _hash_frame(
        digest,
        f"group:{relative}:attributes",
        _canonical(normalized).encode("utf-8"),
    )


def _hash_dataset(digest: Any, path: str, values: Any, numpy: Any) -> None:
    _hash_frame(
        digest,
        f"dataset:{path}:shape",
        _canonical(list(values.shape)).encode("ascii"),
    )
    kind = values.dtype.kind
    if kind in {"O", "S", "U"}:
        payload = bytearray()
        for item in values.reshape(-1):
            encoded = item if isinstance(item, bytes) else str(item).encode("utf-8")
            payload.extend(struct.pack("<Q", len(encoded)))
            payload.extend(encoded)
        dtype = b"utf8"
        raw = bytes(payload)
    elif kind in {"b", "i", "u", "f", "c"}:
        canonical_dtype = values.dtype.newbyteorder("<")
        normalized = numpy.ascontiguousarray(
            values.astype(canonical_dtype, copy=False)
        )
        dtype = f"{kind}{values.dtype.itemsize}-little-endian".encode("ascii")
        raw = normalized.tobytes(order="C")
    else:
        raise ValueError(f"unsupported raw-record dataset dtype {values.dtype}")
    _hash_frame(digest, f"dataset:{path}:dtype", dtype)
    _hash_frame(digest, f"dataset:{path}:value", raw)


def _raw_record_sha256(
    group: Any, numpy: Any, excluded: frozenset[str] = _RECORD_SEAL_EXCLUDED
) -> str:
    """Seal one normalized raw record independently of its source pathname."""

    digest = hashlib.sha256()

    def visit(current: Any, relative: str) -> None:
        _hash_attributes(digest, relative, current.attrs, excluded)
        for name in sorted(current.keys()):
            value = current[name]
            path = f"{relative}/{name}" if relative else str(name)
            if hasattr(value, "keys"):
                visit(value, path)
                continue
            _hash_dataset(digest, path, numpy.asarray(value[...]), numpy)

    visit(group, "")
    return digest.hexdigest()


class _Record:
    """An HDF5 group tree held in memory: group attributes, datasets, filters.

    Raw block records and fold generations are built this way, sealed with
    the same framing `_raw_record_sha256` applies to the stored group, then
    written.  A block that is folded on arrival is never written, but its
    seal is exactly the one its raw record would carry.
    """

    def __init__(self) -> None:
        self.attrs: dict[str, dict[str, Any]] = {"": {}}
        self.datasets: dict[str, Any] = {}
        self.compressed: set[str] = set()

    def group(self, path: str) -> dict[str, Any]:
        parts = path.split("/")
        for depth in range(1, len(parts) + 1):
            self.attrs.setdefault("/".join(parts[:depth]), {})
        return self.attrs[path]

    def dataset(self, path: str, values: Any, dtype: Any, numpy: Any, *, compressed: bool) -> None:
        if "/" in path:
            self.group(path.rsplit("/", 1)[0])
        self.datasets[path] = numpy.asarray(values, dtype=dtype)
        if compressed and self.datasets[path].size:
            self.compressed.add(path)

    def sha256(self, numpy: Any, excluded: frozenset[str]) -> str:
        digest = hashlib.sha256()
        children: dict[str, list[str]] = {}
        for path in [*self.attrs, *self.datasets]:
            if path:
                parent = path.rsplit("/", 1)[0] if "/" in path else ""
                children.setdefault(parent, []).append(path)

        def visit(relative: str) -> None:
            _hash_attributes(digest, relative, self.attrs[relative], excluded)
            for path in sorted(children.get(relative, ()), key=lambda value: value.rsplit("/", 1)[-1]):
                if path in self.attrs:
                    visit(path)
                else:
                    _hash_dataset(digest, path, self.datasets[path], numpy)

        visit("")
        return digest.hexdigest()

    def write(self, parent: Any, name: str) -> Any:
        target = parent.create_group(name)
        for path in sorted(self.attrs, key=lambda value: value.count("/")):
            group = target if not path else target.require_group(path)
            for key, value in self.attrs[path].items():
                group.attrs[key] = value
        for path, values in self.datasets.items():
            kwargs: dict[str, Any] = {}
            if path in self.compressed:
                if values.ndim == 1:
                    chunks: Any = (min(values.shape[0], 131_072),)
                elif values.ndim == 2:
                    chunks = (1, min(values.shape[1], 131_072))
                else:
                    chunks = True
                kwargs.update(
                    chunks=chunks,
                    compression="gzip",
                    compression_opts=_HDF5_GZIP_LEVEL,
                )
            dtype: Any = values.dtype
            if dtype.kind == "O":
                h5py, _ = require_hdf5()
                dtype = h5py.string_dtype(encoding="utf-8")
            target.create_dataset(path, data=values, dtype=dtype, **kwargs)
        return target


def _layout_record(block: FitBlock) -> dict[str, Any]:
    return {
        "manifest": block.manifest,
        "energy_axes": [
            {
                "name": axis["name"],
                "transform_version": axis["transform_version"],
                "registered_edges": axis["registered_edges"],
            }
            for axis in block.axes["energy"]
        ],
        "total_energy_axis": block.axes["total_energy_axis"],
        "mu_edges": block.axes["mu_edges"],
        "lines": _CANONICAL_LINES,
        "components": [
            {key: component[key] for key in ("name", "origin", "post_creation_order", "energy_axis")}
            for component in block.components
        ],
        "stopping_definition_json": block.stopping["definition_json"],
        "diagnostic_channels": [channel["name"] for channel in block.diagnostics["channels"]],
        "interaction_names": block.diagnostics["interaction_names"],
    }


def _same_node(a: FitBlock, b: FitBlock) -> bool:
    keys = (
        "energy_keV",
        "distance_cm",
        "source_height_cm",
        "detector_height_cm",
        "extent_cm",
        "q",
        "role",
        "direct_transmission",
        "direct_fluence",
        "target_total_rel_se",
    )
    if any(a.node[key] != b.node[key] for key in keys):
        return False
    return a.axes["lines"] == b.axes["lines"] and all(
        aa["physical_edges_keV"] == bb["physical_edges_keV"]
        for aa, bb in zip(a.axes["energy"], b.axes["energy"])
    )


def _validate_batch(blocks: Sequence[FitBlock]) -> None:
    if not blocks:
        raise ValueError("no CEEFIT blocks were supplied")
    layout = _canonical(_layout_record(blocks[0]))
    identities: set[tuple[int, int]] = set()
    by_node: dict[int, list[FitBlock]] = {}
    for block in blocks:
        if _canonical(_layout_record(block)) != layout:
            raise ValueError("blocks do not share one archive layout/configuration")
        identity = (block.node_id, block.block_id)
        if identity in identities:
            raise ValueError(f"duplicate block node={block.node_id:016x} id={block.block_id}")
        identities.add(identity)
        by_node.setdefault(block.node_id, []).append(block)
    for node_id, node_blocks in by_node.items():
        first = node_blocks[0]
        if any(not _same_node(first, block) for block in node_blocks[1:]):
            raise ValueError(f"blocks for node {node_id:016x} disagree on node metadata/axes")
        ranges: list[tuple[int, int, int]] = []
        seeds: set[int] = set()
        for block in node_blocks:
            for replica in block.replicas:
                if replica["rng_seed"] in seeds:
                    raise ValueError(f"RNG seed replay within node {node_id:016x}")
                seeds.add(replica["rng_seed"])
                ranges.append(
                    (
                        replica["history_begin"],
                        replica["history_begin"] + replica["history_count"],
                        block.block_id,
                    )
                )
        ranges.sort()
        if any(current[0] < previous[1] for previous, current in zip(ranges, ranges[1:])):
            raise ValueError(f"history ranges overlap within node {node_id:016x}")


def _node_identity_json(block: FitBlock) -> str:
    """Compact exact node identity retained during streaming preflight."""

    keys = (
        "energy_keV",
        "distance_cm",
        "source_height_cm",
        "detector_height_cm",
        "extent_cm",
        "q",
        "role",
        "direct_transmission",
        "direct_fluence",
        "target_total_rel_se",
    )
    return _canonical(
        {
            "node": {key: block.node[key] for key in keys},
            "lines": block.axes["lines"],
            "physical_edges_keV": [
                axis["physical_edges_keV"] for axis in block.axes["energy"]
            ],
        }
    )


def _preflight_block_paths(
    paths: Sequence[Path], output: Path, h5py: Any
) -> dict[int, int]:
    """Validate a batch without retaining dense scientific arrays in memory.

    This deliberately performs a complete first pass, so every input is
    validated before anything is written, while peak memory stays at one
    decoded worker block plus compact identities and history ranges.
    """

    if not paths:
        raise ValueError("no CEEFIT blocks were supplied")
    layout: str | None = None
    identities: set[tuple[int, int]] = set()
    node_identities: dict[int, str] = {}
    ranges_by_node: dict[int, list[tuple[int, int]]] = {}
    seeds_by_node: dict[int, set[int]] = {}
    block_counts: dict[int, int] = {}

    existing = h5py.File(output, "r") if output.exists() else None
    try:
        for path in paths:
            block = read_fit_block(path)
            block_layout = _canonical(_layout_record(block))
            if layout is None:
                layout = block_layout
            elif block_layout != layout:
                raise ValueError("blocks do not share one archive layout/configuration")
            if existing is not None:
                _check_against_archive(existing, block)

            identity = (block.node_id, block.block_id)
            if identity in identities:
                raise ValueError(
                    f"duplicate block node={block.node_id:016x} id={block.block_id}"
                )
            identities.add(identity)
            block_counts[block.node_id] = block_counts.get(block.node_id, 0) + 1

            node_identity = _node_identity_json(block)
            previous_identity = node_identities.setdefault(block.node_id, node_identity)
            if previous_identity != node_identity:
                raise ValueError(
                    f"blocks for node {block.node_id:016x} disagree on node metadata/axes"
                )

            ranges = ranges_by_node.setdefault(block.node_id, [])
            seeds = seeds_by_node.setdefault(block.node_id, set())
            for replica in block.replicas:
                seed = replica["rng_seed"]
                if seed in seeds:
                    raise ValueError(f"RNG seed replay within node {block.node_id:016x}")
                seeds.add(seed)
                start = replica["history_begin"]
                ranges.append((start, start + replica["history_count"]))
            del block
    finally:
        if existing is not None:
            existing.close()

    for node_id, ranges in ranges_by_node.items():
        ranges.sort()
        if any(current[0] < previous[1] for previous, current in zip(ranges, ranges[1:])):
            raise ValueError(f"history ranges overlap within node {node_id:016x}")
    return block_counts


_FIT_CELL_CHUNK = 8192
_HDF5_GZIP_LEVEL = 1


def _aligned_fit_cell_chunk(segment_width: int) -> int:
    """Choose the largest bounded chunk that exactly tiles a logical segment."""

    if type(segment_width) is not int or segment_width <= 0:
        raise ValueError("fit archive segment width must be a positive integer")
    best = 1
    for divisor in range(1, math.isqrt(segment_width) + 1):
        if segment_width % divisor:
            continue
        complement = segment_width // divisor
        if divisor <= _FIT_CELL_CHUNK:
            best = max(best, divisor)
        if complement <= _FIT_CELL_CHUNK:
            best = max(best, complement)
    return best


def _create_extendible(
    group: Any,
    name: str,
    tail: tuple[int, ...],
    dtype: Any,
    *,
    chunks: tuple[int, ...] | None = None,
) -> Any:
    if chunks is None:
        chunks = (1, *tail) if tail else (256,)
    kwargs: dict[str, Any] = {"chunks": chunks}
    if tail:
        kwargs.update(
            compression="gzip",
            compression_opts=_HDF5_GZIP_LEVEL,
            shuffle=True,
        )
    return group.create_dataset(
        name,
        shape=(0, *tail),
        maxshape=(None, *tail),
        dtype=dtype,
        **kwargs,
    )


def _initialize_archive(root: Any, block: FitBlock, h5py: Any, numpy: Any) -> None:
    root.attrs["schema"] = ARCHIVE_SCHEMA
    root.attrs["layout_json"] = _canonical(_layout_record(block))
    manifest = root.create_group("manifest")
    for key, value in block.manifest.items():
        if isinstance(value, (str, bool, int, float)):
            manifest.attrs[key] = value

    axes = root.create_group("axes")
    axes.attrs["total_energy_axis"] = block.axes["total_energy_axis"]
    axes.attrs["line_definitions_json"] = _canonical(_CANONICAL_LINES)
    axes.create_dataset("mu_edges", data=block.axes["mu_edges"], dtype="f8")
    energy_group = axes.create_group("energy")
    for index, axis in enumerate(block.axes["energy"]):
        group = energy_group.create_group(str(index))
        group.attrs["name"] = axis["name"]
        group.attrs["transform_version"] = axis["transform_version"]
        group.create_dataset("registered_edges", data=axis["registered_edges"], dtype="f8")
        _create_extendible(group, "physical_edges_keV", (len(axis["physical_edges_keV"]),), "f8")

    components = root.create_group("components")
    components.attrs["descriptors_json"] = _canonical(
        [
            {key: item[key] for key in ("name", "origin", "post_creation_order", "energy_axis")}
            for item in block.components
        ]
    )
    support_size = sum(len(item["support_mask"]) for item in block.components)
    offsets = [0]
    for item in block.components:
        offsets.append(offsets[-1] + len(item["support_mask"]))
    components.create_dataset("support_offsets", data=offsets, dtype="u8")
    _create_extendible(components, "support_mask", (support_size,), "u1")

    nodes = root.create_group("nodes")
    _create_extendible(nodes, "id", (), "u8")
    _create_extendible(nodes, "coordinates", (5,), "f8")
    _create_extendible(nodes, "q", (4,), "f8")
    text_dtype = h5py.string_dtype(encoding="utf-8")
    _create_extendible(nodes, "role", (), text_dtype)
    _create_extendible(nodes, "status", (), text_dtype)
    _create_extendible(nodes, "failure_reason", (), text_dtype)
    for name, dtype in (
        ("direct_transmission", "f8"),
        ("direct_fluence", "f8"),
        ("histories", "u8"),
        ("wall_seconds", "f8"),
        ("cpu_seconds", "f8"),
        ("target_total_rel_se", "f8"),
        ("achieved_total_rel_se", "f8"),
        ("transport_failures", "u8"),
        ("max_history_contribution", "f8"),
        ("capped", "u1"),
    ):
        _create_extendible(nodes, name, (), dtype)

    continuum_root = root.create_group("continuum")
    line_root = root.create_group("lines")
    _create_extendible(line_root, "energy_keV", (len(_CANONICAL_LINES),), "f8")
    _create_extendible(line_root, "support_mask", (len(_CANONICAL_LINES),), "u1")
    for name, moment in block.moments.items():
        family, label = name.split("_", 1)
        group = (continuum_root if family == "continuum" else line_root).create_group(label)
        if name == "line_components":
            size = len(block.components) * len(_CANONICAL_LINES) * (len(block.axes["mu_edges"]) - 1)
        elif name == "line_total":
            size = len(_CANONICAL_LINES) * (len(block.axes["mu_edges"]) - 1)
        else:
            size = len(moment["N"])
        if name == "continuum_components":
            if size % len(block.components):
                raise ValueError("continuum component moments are not component-major")
            segment_width = size // len(block.components)
        elif name in {"line_components", "line_total"}:
            segment_width = len(block.axes["mu_edges"]) - 1
        else:
            segment_width = size
        cell_chunk = _aligned_fit_cell_chunk(segment_width)
        for field_name, dtype in (("N", "u8"), ("S1", "f8"), ("S2", "f8")):
            _create_extendible(
                group,
                field_name,
                (size,),
                dtype,
                chunks=(1, cell_chunk),
            )
        _create_extendible(
            group,
            "replica_S1",
            (32, size),
            "f4",
            chunks=(1, 32, cell_chunk),
        )
    replica_root = root.create_group("replicas")
    _create_extendible(replica_root, "histories", (32,), "u8")

    stopping = root.create_group("stopping")
    stopping.attrs["definition_json"] = block.stopping["definition_json"]
    _create_extendible(stopping, "N", (), "u8")
    _create_extendible(stopping, "S1", (16,), "f8")
    _create_extendible(stopping, "outer_S2", (16, 16), "f8")
    _create_extendible(stopping, "standard_error", (16,), "f8")
    _create_extendible(stopping, "target_standard_error", (16,), "f8")
    _create_extendible(stopping, "active", (16,), "u1")
    _create_extendible(stopping, "passed", (16,), "u1")

    diagnostics = root.create_group("diagnostics")
    diagnostics.attrs["interaction_names_json"] = _canonical(
        block.diagnostics["interaction_names"]
    )
    if block.diagnostics["interaction_counts"]:
        _create_extendible(
            diagnostics,
            "interaction_counts",
            (len(block.diagnostics["interaction_counts"]),),
            "u8",
        )
    else:
        diagnostics.create_dataset("interaction_counts", shape=(0,), dtype="u8")
        diagnostics["interaction_counts"].attrs["empty_per_node"] = True
    _create_extendible(diagnostics, "dropped_N", (), "u8")
    _create_extendible(diagnostics, "dropped_S1", (), "f8")
    _create_extendible(diagnostics, "dropped_S2", (), "f8")
    _create_extendible(diagnostics, "scattered_total_N", (), "u8")
    _create_extendible(diagnostics, "scattered_total_S1", (), "f8")
    _create_extendible(diagnostics, "scattered_total_S2", (), "f8")
    channels = diagnostics.create_group("path_channels")
    for channel in block.diagnostics["channels"]:
        group = channels.create_group(channel["name"])
        for tensor_name in ("continuum", "lines"):
            tensor = group.create_group(tensor_name)
            size = (
                len(_CANONICAL_LINES) * (len(block.axes["mu_edges"]) - 1)
                if tensor_name == "lines"
                else len(channel[tensor_name]["N"])
            )
            for field_name, dtype in (("N", "u8"), ("S1", "f8"), ("S2", "f8")):
                _create_extendible(tensor, field_name, (size,), dtype)
    root.create_group("blocks")


def _resize_node_axis(root: Any, new_size: int) -> None:
    def walk(group: Any) -> None:
        for value in group.values():
            if hasattr(value, "shape"):
                if value.maxshape and value.maxshape[0] is None:
                    value.resize(new_size, axis=0)
            else:
                walk(value)

    # Block records and fold states are fixed-size and live below /blocks and
    # /folds; do not resize.
    for name, group in root.items():
        if name not in {"blocks", "folds", "manifest"}:
            if hasattr(group, "values"):
                walk(group)


def _node_index(root: Any, node_id: int) -> int | None:
    values = root["nodes/id"][:]
    matches = [i for i, value in enumerate(values) if int(value) == node_id]
    if len(matches) > 1:
        raise RuntimeError(f"archive contains duplicate node id {node_id:016x}")
    return matches[0] if matches else None


_LEDGER_ARRAYS = ("id", "history_begin", "history_count", "rng_seed")


def _fold_generation(root: Any, node_hex: str) -> Any | None:
    """The current fold generation of a node, or None if nothing is folded."""

    folds = root.get("folds")
    node_group = None if folds is None else folds.get(node_hex)
    if node_group is None or "current" not in node_group.attrs:
        # No generation has become current yet (a first write was interrupted).
        return None
    if str(node_group.attrs.get("fold_schema", "")) != FOLD_STATE_SCHEMA:
        raise ValueError(f"fold state of node {node_hex} lacks the {FOLD_STATE_SCHEMA} schema")
    current = str(node_group.attrs.get("current", ""))
    if current not in node_group:
        raise ValueError(f"fold state of node {node_hex} names a missing generation {current!r}")
    return node_group[current]


def _ledger_rows(generation: Any) -> list[dict[str, Any]]:
    """Per-block rows of a fold generation: scalars, checksums, replica arrays."""

    rows = json.loads(generation["ledger_json"][()])
    arrays = {name: generation[f"ledger/{name}"][...] for name in _LEDGER_ARRAYS}
    for index, row in enumerate(rows):
        for name in _LEDGER_ARRAYS:
            row[name] = tuple(int(value) for value in arrays[name][index])
    return rows


def node_block_inventory(root: Any, node_hex: str) -> dict[int, dict[str, Any]]:
    """Every archived block of one node: raw records and folded ledger rows.

    Each entry carries the block's checksums, history count and replica
    id/begin/count/seed arrays, and whether a raw record and/or a ledger row
    holds it (blocks 1..RETAINED_RAW_BLOCKS of a folded node have both, and
    both must carry the same checksums).
    """

    inventory: dict[int, dict[str, Any]] = {}
    generation = _fold_generation(root, node_hex)
    if generation is not None:
        for row in _ledger_rows(generation):
            inventory[int(row["block_id"])] = {
                "block_id": int(row["block_id"]),
                "histories": int(row["histories"]),
                "payload_sha256": str(row["payload_sha256"]),
                "scientific_record_sha256": str(row["scientific_record_sha256"]),
                **{name: row[name] for name in _LEDGER_ARRAYS},
                "raw": False,
                "folded": True,
            }
    node_blocks = root["blocks"].get(node_hex)
    for name, group in ({} if node_blocks is None else node_blocks).items():
        block_id = int(name)
        entry = {
            "block_id": block_id,
            "histories": int(group.attrs["histories"]),
            "payload_sha256": str(group.attrs.get("payload_sha256", "")),
            "scientific_record_sha256": str(group.attrs.get("scientific_record_sha256", "")),
            **{
                key: tuple(int(value) for value in group[key][:]) if key in group else ()
                for key in _LEDGER_ARRAYS
            },
            "raw": True,
            "folded": False,
        }
        folded = inventory.get(block_id)
        if folded is not None:
            if any(folded[key] != entry[key] for key in ("payload_sha256", "scientific_record_sha256", "histories", *_LEDGER_ARRAYS)):
                raise ValueError(
                    f"raw block {node_hex}/{block_id} differs from its fold-ledger row"
                )
            entry["folded"] = True
        inventory[block_id] = entry
    return dict(sorted(inventory.items()))


def archived_block_inventory(root: Any) -> dict[str, dict[int, dict[str, Any]]]:
    """node_block_inventory() for every node that has archived blocks."""

    names = set(root["blocks"])
    if "folds" in root:
        names |= set(root["folds"])
    return {
        node_hex: node_block_inventory(root, node_hex)
        for node_hex in sorted(names, key=lambda value: int(value, 16))
    }


def _check_against_archive(root: Any, block: FitBlock) -> None:
    if root.attrs.get("schema") != ARCHIVE_SCHEMA:
        raise ValueError("output exists but is not a supported CEEFIT HDF5 archive")
    if root.attrs.get("layout_json") != _canonical(_layout_record(block)):
        raise ValueError("block configuration/layout does not match the HDF5 archive")
    node_hex = f"{block.node_id:016x}"
    inventory = node_block_inventory(root, node_hex)
    if block.block_id in inventory:
        raise ValueError(f"block {block.block_id} for node {node_hex} is already archived")
    old_ranges = [
        (begin, begin + count)
        for entry in inventory.values()
        for begin, count in zip(entry["history_begin"], entry["history_count"])
    ]
    old_seeds = {seed for entry in inventory.values() for seed in entry["rng_seed"]}
    for replica in block.replicas:
        start = replica["history_begin"]
        end = start + replica["history_count"]
        if any(start < old_end and old_start < end for old_start, old_end in old_ranges):
            raise ValueError(f"block history range overlaps archived node {node_hex}")
        if replica["rng_seed"] in old_seeds:
            raise ValueError(f"block RNG seed replays archived node {node_hex}")


def _block_record(block: FitBlock, numpy: Any) -> _Record:
    """The normalized raw record of one worker block, sealed, not yet written."""

    record = _Record()
    attrs = record.attrs[""]
    attrs["record_schema"] = RAW_RECORD_SCHEMA
    attrs["node_id"] = block.node_id
    attrs["block_id"] = block.block_id
    attrs["layout_json"] = _canonical(_layout_record(block))
    attrs["payload_sha256"] = block.payload_sha256
    attrs["stopping_definition_json"] = block.stopping["definition_json"]
    attrs["source_path"] = block.source_path or ""
    attrs["role"] = block.node["role"]
    attrs["status"] = block.node["status"]
    attrs["failure_reason"] = block.diagnostics["failure_reason"]
    attrs["capped"] = block.diagnostics["capped"]
    attrs["transport_failures"] = block.diagnostics["transport_failures"]
    attrs["wall_seconds"] = block.node["wall_seconds"]
    attrs["cpu_seconds"] = block.node["cpu_seconds"]
    attrs["histories"] = block.node["histories"]
    attrs["direct_transmission"] = block.node["direct_transmission"]
    attrs["direct_fluence"] = block.node["direct_fluence"]
    attrs["target_total_rel_se"] = block.node["target_total_rel_se"]
    attrs["max_history_contribution"] = block.diagnostics["max_history_contribution"]
    attrs["achieved_total_rel_se"] = block.node["achieved_total_rel_se"]
    attrs["scattered_total_N"] = block.diagnostics["scattered_total"]["N"]
    attrs["scattered_total_S1"] = block.diagnostics["scattered_total"]["S1"]
    attrs["scattered_total_S2"] = block.diagnostics["scattered_total"]["S2"]
    attrs["dropped_N"] = block.diagnostics["dropped"]["N"]
    attrs["dropped_S1"] = block.diagnostics["dropped"]["S1"]
    attrs["dropped_S2"] = block.diagnostics["dropped"]["S2"]

    def scientific(path: str, values: Any, dtype: str) -> None:
        record.dataset(path, values, dtype, numpy, compressed=True)

    def compact(path: str, values: Any, dtype: str) -> None:
        record.dataset(path, values, dtype, numpy, compressed=False)

    # These compact deterministic records make the immutable worker blocks a
    # complete source of truth for every node-varying HDF5 coordinate/axis.
    compact(
        "coordinates",
        (
            block.node["energy_keV"],
            block.node["distance_cm"],
            block.node["source_height_cm"],
            block.node["detector_height_cm"],
            block.node["extent_cm"],
        ),
        "f8",
    )
    compact("q", block.node["q"], "f8")
    record.group("physical_energy_edges_keV")
    for axis_index, axis in enumerate(block.axes["energy"]):
        scientific(f"physical_energy_edges_keV/{axis_index}", axis["physical_edges_keV"], "f8")
    scientific(
        "component_support_mask",
        tuple(
            value
            for component in block.components
            for value in component["support_mask"]
        ),
        "u1",
    )
    support_offsets = [0]
    for component in block.components:
        support_offsets.append(support_offsets[-1] + len(component["support_mask"]))
    compact("component_support_offsets", support_offsets, "u8")
    line_energy = {line["name"]: line["energy_keV"] for line in block.axes["lines"]}
    compact(
        "line_energy_keV",
        (
            line_energy.get("elastic_E0", block.node["energy_keV"]),
            line_energy.get("annihilation_511", 510.998950),
        ),
        "f8",
    )
    compact(
        "line_support_mask",
        tuple(int(descriptor["name"] in line_energy) for descriptor in _CANONICAL_LINES),
        "u1",
    )
    compact("interaction_counts", block.diagnostics["interaction_counts"], "u8")
    for key in ("id", "history_begin", "history_count", "rng_seed"):
        compact(key, [item[key] for item in block.replicas], "u8")
    record.group("moments")
    record.group("replica_s1")
    for tensor_name in block.moments:
        for field_name, dtype in (("N", "u8"), ("S1", "f8"), ("S2", "f8")):
            scientific(
                f"moments/{tensor_name}/{field_name}",
                _expanded_moment_values(block, tensor_name, field_name),
                dtype,
            )
        replica_values = (
            _expanded_replica_lines(block, tensor_name, numpy)
            if tensor_name.startswith("line_")
            else numpy.asarray(block.replica_s1[tensor_name], dtype="f4").reshape(32, -1)
        )
        scientific(f"replica_s1/{tensor_name}", replica_values, "f4")
    record.group("stopping")["N"] = block.stopping["N"]
    compact("stopping/S1", block.stopping["S1"], "f8")
    compact(
        "stopping/outer_S2",
        numpy.asarray(block.stopping["outer_S2"], dtype="f8").reshape(16, 16),
        "f8",
    )
    record.group("path_channels")
    for channel in block.diagnostics["channels"]:
        for tensor_name in ("continuum", "lines"):
            record.group(f"path_channels/{channel['name']}/{tensor_name}")
            for field_name, dtype in (("N", "u8"), ("S1", "f8"), ("S2", "f8")):
                values = channel[tensor_name][field_name]
                if tensor_name == "lines":
                    values = _expand_line_values(
                        block,
                        values,
                        components=False,
                        missing_value=(block.node["histories"] if field_name == "N" else 0.0),
                    )
                scientific(f"path_channels/{channel['name']}/{tensor_name}/{field_name}", values, dtype)

    # Seal the normalized HDF5 values, rather than trusting only the copied
    # worker payload checksum attribute.  The source path is deliberately not
    # part of this seal so archives remain invariant under staging layout.
    attrs["scientific_record_sha256"] = record.sha256(numpy, _RECORD_SEAL_EXCLUDED)
    return record


def _write_block_record(root: Any, block: FitBlock, numpy: Any, record: _Record | None = None) -> None:
    node_group = root["blocks"].require_group(f"{block.node_id:016x}")
    (record or _block_record(block, numpy)).write(node_group, str(block.block_id))


def _merged_node_status(
    previous_status: str,
    incoming_status: str,
    *,
    transport_failures: int,
    capped: bool,
    failure_reason: str,
) -> str:
    """Merge status without misclassifying a documented cap as a failure."""

    if transport_failures or "failed" in {previous_status, incoming_status}:
        return "failed"
    if capped or "capped" in {previous_status, incoming_status}:
        return "capped"
    if failure_reason:
        return "failed"
    if "converged" in {previous_status, incoming_status}:
        return "converged"
    return "complete"


def _expand_line_values(
    block: FitBlock,
    values: Sequence[Any],
    *,
    components: bool,
    missing_value: Any,
) -> tuple[Any, ...]:
    """Expand threshold-dependent block line tensors onto the fixed two-line axis."""

    n_mu = len(block.axes["mu_edges"]) - 1
    source_names = [line["name"] for line in block.axes["lines"]]
    n_component = len(block.components) if components else 1
    source_per_component = len(source_names) * n_mu
    output: list[Any] = []
    for component in range(n_component):
        base = component * source_per_component
        for descriptor in _CANONICAL_LINES:
            if descriptor["name"] in source_names:
                source_line = source_names.index(descriptor["name"])
                start = base + source_line * n_mu
                output.extend(values[start : start + n_mu])
            else:
                output.extend([missing_value] * n_mu)
    return tuple(output)


def _expanded_moment_values(
    block: FitBlock, tensor_name: str, field_name: str
) -> tuple[Any, ...]:
    values = block.moments[tensor_name][field_name]
    if not tensor_name.startswith("line_"):
        return values
    missing = block.node["histories"] if field_name == "N" else 0.0
    return _expand_line_values(
        block,
        values,
        components=tensor_name == "line_components",
        missing_value=missing,
    )


def _expanded_replica_lines(block: FitBlock, tensor_name: str, numpy: Any) -> Any:
    values = block.replica_s1[tensor_name]
    source_base = len(values) // 32
    rows = []
    for replica in range(32):
        rows.append(
            _expand_line_values(
                block,
                values[replica * source_base : (replica + 1) * source_base],
                components=tensor_name == "line_components",
                missing_value=0.0,
            )
        )
    return numpy.asarray(rows, dtype="f4")


# ---------------------------------------------------------------------------
# Node folds.  A node's materialized row is a fold over its blocks in
# ascending block-ID order: compensated (Kahan) f8 sums of every float tensor,
# exact u8 sums of the counts, and fsum / max / join merges of the per-block
# scalars.  The fold state after blocks 1..k -- every (total, compensation)
# pair, every integer sum and one ledger row per block -- is kept as the
# generation /folds/<node>/<k>.  Resuming the fold from that state reproduces
# the fold over all blocks bit for bit, so a block beyond RETAINED_RAW_BLOCKS
# is folded as it is ingested and never stored as a raw record.
# ---------------------------------------------------------------------------

_ROW_SCALARS = (
    "block_id",
    "payload_sha256",
    "scientific_record_sha256",
    "histories",
    "wall_seconds",
    "cpu_seconds",
    "transport_failures",
    "capped",
    "failure_reason",
    "max_history_contribution",
    "status",
    "scattered_total_N",
    "scattered_total_S1",
    "scattered_total_S2",
    "dropped_N",
    "dropped_S1",
    "dropped_S2",
)
_DETERMINISTIC_ATTRIBUTES = (
    "role",
    "direct_transmission",
    "direct_fluence",
    "target_total_rel_se",
)
_DETERMINISTIC_DATASETS = (
    "coordinates",
    "q",
    "component_support_mask",
    "component_support_offsets",
    "line_energy_keV",
    "line_support_mask",
)
# Materialized datasets written once, when a node first appears.
_INGEST_TIME_PATHS = frozenset(
    {
        "nodes/id",
        "nodes/coordinates",
        "nodes/q",
        "nodes/role",
        "components/support_mask",
        "lines/energy_keV",
        "lines/support_mask",
    }
)


class _BlockView:
    """Read access to one block for folding: a stored raw record or a _Record."""

    def __init__(self, source: Any) -> None:
        self.source = source
        self.in_memory = isinstance(source, _Record)

    def attrs(self) -> Mapping[str, Any]:
        return self.source.attrs[""] if self.in_memory else self.source.attrs

    def array(self, path: str) -> Any:
        return self.source.datasets[path] if self.in_memory else self.source[path][...]

    def names(self, group: str) -> list[str]:
        if not self.in_memory:
            return sorted(self.source[group])
        prefix = f"{group}/"
        return sorted(
            {
                path[len(prefix):].split("/", 1)[0]
                for path in (*self.source.attrs, *self.source.datasets)
                if path.startswith(prefix)
            }
        )

    def stopping_n(self) -> int:
        if self.in_memory:
            return int(self.source.attrs["stopping"]["N"])
        return int(self.source["stopping"].attrs["N"])


def _fold_paths(view: _BlockView, per_node_interactions: bool) -> tuple[list[str], list[str]]:
    floats: list[str] = []
    integers: list[str] = []
    for tensor in view.names("moments"):
        integers.append(f"moments/{tensor}/N")
        floats += [f"moments/{tensor}/S1", f"moments/{tensor}/S2", f"replica_s1/{tensor}"]
    integers.append("history_count")
    floats += ["stopping/S1", "stopping/outer_S2"]
    if per_node_interactions:
        integers.append("interaction_counts")
    for channel in view.names("path_channels"):
        for tensor in ("continuum", "lines"):
            base = f"path_channels/{channel}/{tensor}"
            integers.append(f"{base}/N")
            floats += [f"{base}/S1", f"{base}/S2"]
    return floats, integers


def _view_row(view: _BlockView) -> dict[str, Any]:
    attrs = view.attrs()
    row = {name: _normal_hdf_value(attrs[name]) for name in _ROW_SCALARS}
    row["stopping_N"] = view.stopping_n()
    for name in _LEDGER_ARRAYS:
        row[name] = tuple(int(value) for value in view.array(name))
    return row


def _new_fold_state(view: _BlockView, per_node_interactions: bool, numpy: Any) -> dict[str, Any]:
    floats, integers = _fold_paths(view, per_node_interactions)
    state: dict[str, Any] = {"through": 0, "float": {}, "int": {}, "rows": []}
    for path in floats:
        total = numpy.zeros_like(numpy.asarray(view.array(path), dtype="f8"), dtype="f8")
        state["float"][path] = [total, numpy.zeros_like(total, dtype="f8")]
    for path in integers:
        state["int"][path] = numpy.zeros_like(
            numpy.asarray(view.array(path), dtype="u8"), dtype="u8"
        )
    return state


def _copy_fold_state(state: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "through": state["through"],
        "float": {path: [pair[0].copy(), pair[1].copy()] for path, pair in state["float"].items()},
        "int": {path: total.copy() for path, total in state["int"].items()},
        "rows": [dict(row) for row in state["rows"]],
    }


def _fold_view(state: dict[str, Any], view: _BlockView, numpy: Any) -> None:
    """Add one block, exactly as the ascending-order compensated fold does."""

    for path, pair in state["float"].items():
        value = numpy.asarray(view.array(path), dtype="f8")
        corrected = value - pair[1]
        updated = pair[0] + corrected
        pair[1] = (updated - pair[0]) - corrected
        pair[0] = updated
    for path, total in state["int"].items():
        numpy.add(total, numpy.asarray(view.array(path), dtype="u8"), out=total)
    row = _view_row(view)
    state["rows"].append(row)
    state["through"] = int(row["block_id"])


def _ledger_chain_sha256(rows: Sequence[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256(b"CEEFIT-FOLD-LEDGER-v1\0").digest()
    for row in rows:
        link = _canonical(
            [int(row["block_id"]), str(row["payload_sha256"]), str(row["scientific_record_sha256"])]
        ).encode("utf-8")
        digest = hashlib.sha256(digest + link).digest()
    return digest.hex()


def _fold_generation_record(state: Mapping[str, Any], numpy: Any) -> _Record:
    record = _Record()
    attrs = record.attrs[""]
    attrs["folded_through"] = int(state["through"])
    attrs["paths_json"] = _canonical({"float": sorted(state["float"]), "int": sorted(state["int"])})
    attrs["ledger_chain_sha256"] = _ledger_chain_sha256(state["rows"])
    for path, (total, compensation) in state["float"].items():
        record.dataset(f"float/{path}/total", total, "f8", numpy, compressed=True)
        record.dataset(f"float/{path}/compensation", compensation, "f8", numpy, compressed=True)
    for path, total in state["int"].items():
        record.dataset(f"int/{path}", total, "u8", numpy, compressed=True)
    scalars = [
        {key: value for key, value in row.items() if key not in _LEDGER_ARRAYS}
        for row in state["rows"]
    ]
    record.dataset("ledger_json", _canonical(scalars), object, numpy, compressed=False)
    for name in _LEDGER_ARRAYS:
        record.dataset(
            f"ledger/{name}", [row[name] for row in state["rows"]], "u8", numpy, compressed=False
        )
    attrs["state_sha256"] = record.sha256(numpy, _FOLD_SEAL_EXCLUDED)
    return record


def _load_fold_state(generation: Any, numpy: Any) -> dict[str, Any]:
    """A stored fold state, checked against its seal before anything uses it:
    re-sealing the loaded state must reproduce `state_sha256` (which also
    covers the ledger and its hash chain), so a damaged state is never
    extended and re-sealed as if it were valid."""

    paths = json.loads(str(generation.attrs["paths_json"]))
    state = {
        "through": int(generation.attrs["folded_through"]),
        "float": {
            path: [
                numpy.asarray(generation[f"float/{path}/total"][...], dtype="f8"),
                numpy.asarray(generation[f"float/{path}/compensation"][...], dtype="f8"),
            ]
            for path in paths["float"]
        },
        "int": {
            path: numpy.asarray(generation[f"int/{path}"][...], dtype="u8")
            for path in paths["int"]
        },
        "rows": _ledger_rows(generation),
    }
    sealed = str(generation.attrs.get("state_sha256", ""))
    if _fold_generation_record(state, numpy).attrs[""]["state_sha256"] != sealed:
        raise ValueError(f"fold state {generation.name} failed its checksum")
    return state


def _store_fold_generation(root: Any, node_hex: str, state: Mapping[str, Any], numpy: Any) -> str | None:
    """Write a new current generation; return the one it replaces (kept until verified)."""

    node_group = root.require_group("folds").require_group(node_hex)
    previous = str(node_group.attrs.get("current", "")) or None
    name = str(int(state["through"]))
    if name == previous:
        raise RuntimeError(f"fold state of node {node_hex} did not advance")
    if name in node_group:
        # Left by an interrupted earlier ingest that never became current.
        del node_group[name]
    _fold_generation_record(state, numpy).write(node_group, name)
    # The pointer moves only once the generation is complete.
    node_group.attrs["fold_schema"] = FOLD_STATE_SCHEMA
    node_group.attrs["current"] = name
    return previous


def _drop_stale_fold_generations(root: Any, node_hex: str) -> int:
    folds = root.get("folds")
    node_group = None if folds is None else folds.get(node_hex)
    if node_group is None:
        return 0
    current = str(node_group.attrs.get("current", ""))
    stale = [name for name in node_group if name != current]
    for name in stale:
        del node_group[name]
    return len(stale)


def _check_fold_consistency(
    node_hex: str, state: Mapping[str, Any], raw: Mapping[int, Any]
) -> None:
    """A fold covers exactly blocks 1..k, keeps blocks 1..RETAINED_RAW_BLOCKS
    raw, and agrees with every raw record it covers."""

    through = int(state["through"])
    ids = [int(row["block_id"]) for row in state["rows"]]
    if ids != list(range(1, through + 1)):
        raise ValueError(f"fold ledger of node {node_hex} is not blocks 1..{through}")
    missing = [
        block_id
        for block_id in range(1, min(through, RETAINED_RAW_BLOCKS) + 1)
        if block_id not in raw
    ]
    if missing:
        raise ValueError(f"node {node_hex} lacks its retained raw blocks {missing}")
    for row in state["rows"]:
        group = raw.get(int(row["block_id"]))
        if group is not None and (
            str(group.attrs.get("payload_sha256", "")) != row["payload_sha256"]
            or str(group.attrs.get("scientific_record_sha256", "")) != row["scientific_record_sha256"]
        ):
            raise ValueError(
                f"raw block {node_hex}/{row['block_id']} differs from its fold-ledger row"
            )


def _fold_node(
    root: Any,
    node_hex: str,
    numpy: Any,
    *,
    state: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], _BlockView]:
    """The node's full fold: its stored (or given) fold state extended by the
    raw records beyond it, and a view of its lowest raw block."""

    node_blocks = root["blocks"].get(node_hex)
    raw = {} if node_blocks is None else {int(name): group for name, group in node_blocks.items()}
    if not raw:
        raise RuntimeError(f"node {node_hex} has no archived blocks")
    first = _BlockView(raw[min(raw)])
    if state is None:
        generation = _fold_generation(root, node_hex)
        state = None if generation is None else _load_fold_state(generation, numpy)
    if state is not None:
        _check_fold_consistency(node_hex, state, raw)
    through = 0 if state is None else int(state["through"])
    tail = [raw[block_id] for block_id in sorted(raw) if block_id > through]
    if state is None:
        state = _new_fold_state(
            first, root["diagnostics/interaction_counts"].ndim > 1, numpy
        )
    elif tail:
        state = _copy_fold_state(state)
    for group in tail:
        _fold_view(state, _BlockView(group), numpy)
    return state, first  # type: ignore[return-value]


def _node_values(
    state: Mapping[str, Any], first: _BlockView, node_hex: str, numpy: Any
) -> dict[str, Any]:
    """Every materialized node-major value a fold state determines."""

    rows = state["rows"]

    def total(path: str) -> Any:
        return state["float"][path][0]

    attrs = first.attrs()
    values: dict[str, Any] = {
        "nodes/id": int(node_hex, 16),
        "nodes/coordinates": first.array("coordinates"),
        "nodes/q": first.array("q"),
        "nodes/role": str(attrs["role"]),
        "nodes/direct_transmission": float(attrs["direct_transmission"]),
        "nodes/direct_fluence": float(attrs["direct_fluence"]),
        "nodes/target_total_rel_se": float(attrs["target_total_rel_se"]),
        "nodes/histories": sum(int(row["histories"]) for row in rows),
        "nodes/wall_seconds": math.fsum(float(row["wall_seconds"]) for row in rows),
        "nodes/cpu_seconds": math.fsum(float(row["cpu_seconds"]) for row in rows),
        "components/support_mask": first.array("component_support_mask"),
        "lines/energy_keV": first.array("line_energy_keV"),
        "lines/support_mask": first.array("line_support_mask"),
    }
    for axis_name in first.names("physical_energy_edges_keV"):
        values[f"axes/energy/{axis_name}/physical_edges_keV"] = first.array(
            f"physical_energy_edges_keV/{axis_name}"
        )

    failures = sum(int(row["transport_failures"]) for row in rows)
    capped = any(bool(row["capped"]) for row in rows)
    reasons = [str(row["failure_reason"]) for row in rows if str(row["failure_reason"])]
    failure_reason = "; ".join(reasons)
    values["nodes/transport_failures"] = failures
    values["nodes/capped"] = int(capped)
    values["nodes/failure_reason"] = failure_reason
    values["nodes/max_history_contribution"] = max(
        float(row["max_history_contribution"]) for row in rows
    )
    statuses = {str(row["status"]) for row in rows}
    representative = (
        "failed"
        if "failed" in statuses
        else "capped"
        if "capped" in statuses
        else "converged"
        if "converged" in statuses
        else "complete"
    )
    values["nodes/status"] = _merged_node_status(
        "",
        representative,
        transport_failures=failures,
        capped=capped,
        failure_reason=failure_reason,
    )

    for tensor_name in first.names("moments"):
        family, label = tensor_name.split("_", 1)
        destination = f"{'lines' if family == 'line' else family}/{label}"
        values[f"{destination}/N"] = state["int"][f"moments/{tensor_name}/N"]
        values[f"{destination}/S1"] = total(f"moments/{tensor_name}/S1")
        values[f"{destination}/S2"] = total(f"moments/{tensor_name}/S2")
        values[f"{destination}/replica_S1"] = total(f"replica_s1/{tensor_name}")
    values["replicas/histories"] = state["int"]["history_count"]
    stopping_n = sum(int(row["stopping_N"]) for row in rows)
    stopping_s1 = total("stopping/S1")
    stopping_outer = total("stopping/outer_S2")
    values["stopping/N"] = stopping_n
    values["stopping/S1"] = stopping_s1
    values["stopping/outer_S2"] = stopping_outer
    if "interaction_counts" in state["int"]:
        values["diagnostics/interaction_counts"] = state["int"]["interaction_counts"]
    for prefix in ("scattered_total", "dropped"):
        values[f"diagnostics/{prefix}_N"] = sum(int(row[f"{prefix}_N"]) for row in rows)
        values[f"diagnostics/{prefix}_S1"] = math.fsum(float(row[f"{prefix}_S1"]) for row in rows)
        values[f"diagnostics/{prefix}_S2"] = math.fsum(float(row[f"{prefix}_S2"]) for row in rows)
    for channel_name in first.names("path_channels"):
        for tensor_name in ("continuum", "lines"):
            base = f"path_channels/{channel_name}/{tensor_name}"
            destination = f"diagnostics/{base}"
            values[f"{destination}/N"] = state["int"][f"{base}/N"]
            values[f"{destination}/S1"] = total(f"{base}/S1")
            values[f"{destination}/S2"] = total(f"{base}/S2")

    total_n = int(values["diagnostics/scattered_total_N"])
    total_s1 = float(values["diagnostics/scattered_total_S1"])
    total_s2 = float(values["diagnostics/scattered_total_S2"])
    if total_n >= 2 and total_s1 != 0.0:
        variance = max(0.0, total_s2 - total_s1 * total_s1 / total_n) / (
            total_n * (total_n - 1)
        )
        achieved = math.sqrt(variance) / abs(total_s1 / total_n)
    else:
        achieved = 0.0 if total_s1 == 0.0 and total_s2 == 0.0 else math.inf
    values["nodes/achieved_total_rel_se"] = achieved

    if stopping_n >= 2:
        diagonal = numpy.diag(stopping_outer)
        numerator = numpy.maximum(0.0, diagonal - stopping_s1 * stopping_s1 / stopping_n)
        stopping_se = numpy.sqrt(numerator / (stopping_n * (stopping_n - 1)))
    else:
        stopping_se = numpy.full(16, math.inf, dtype="f8")
    total_mean = total_s1 / total_n if total_n else 0.0
    target = numpy.concatenate(
        (
            numpy.full(12, 0.001 * total_mean, dtype="f8"),
            numpy.full(4, 0.0025 * total_mean, dtype="f8"),
        )
    )
    active = numpy.zeros(16, dtype="u1")
    if total_mean > 0.0 and stopping_n > 0:
        active[:12] = (stopping_s1[:12] / stopping_n >= 0.01 * total_mean).astype("u1")
        active[12:] = 1
    values["stopping/standard_error"] = stopping_se
    values["stopping/target_standard_error"] = target
    values["stopping/active"] = active
    values["stopping/passed"] = numpy.logical_or(active == 0, stopping_se <= target).astype("u1")
    return values


def _rematerialize_node(
    root: Any,
    index: int,
    node_hex: str,
    numpy: Any,
    *,
    state: Mapping[str, Any] | None = None,
) -> None:
    """Rebuild one node's row from its fold state and the raw records beyond it.

    Per-block sufficient statistics are the source of truth, folded in
    ascending block-ID order, so the HDF aggregate is bitwise independent of
    ingest order and of how many blocks were folded before.
    """

    full, first = _fold_node(root, node_hex, numpy, state=state)
    values = _node_values(full, first, node_hex, numpy)
    existing = int(root["nodes/histories"][index])
    if int(values["nodes/histories"]) < existing:
        raise ValueError(
            f"rebuilding node {node_hex} would drop histories ({existing} archived, "
            f"{values['nodes/histories']} in its blocks): its block records are incomplete"
        )
    for path, value in values.items():
        if path not in _INGEST_TIME_PATHS and not path.startswith("axes/"):
            root[path][index] = value


def _verify_materialized_node(
    root: Any, index: int, node_hex: str, numpy: Any
) -> set[str]:
    """Verify one materialized row exactly against its raw blocks and fold state."""

    node_blocks = root[f"blocks/{node_hex}"]
    blocks = [
        node_blocks[name] for name in sorted(node_blocks, key=lambda value: int(value))
    ]
    if not blocks:
        raise ValueError(f"node {node_hex} has no archived blocks")
    for block in blocks:
        block_name = block.name.rsplit("/", 1)[-1]
        schema = str(block.attrs.get("record_schema", ""))
        seal = str(block.attrs.get("scientific_record_sha256", ""))
        if schema != RAW_RECORD_SCHEMA or len(seal) != 64:
            raise ValueError(
                f"raw block record {node_hex}/{block_name} lacks the "
                f"{RAW_RECORD_SCHEMA} schema or its scientific-record seal; "
                "re-ingest the original CEEFIT blocks into a new archive"
            )
        if seal != _raw_record_sha256(block, numpy):
            raise ValueError(
                f"raw block record {node_hex}/{block_name} failed its scientific "
                "checksum; re-ingest the original CEEFIT block"
            )
        if int(block.attrs.get("node_id", -1)) != int(node_hex, 16):
            raise ValueError(f"raw block record {node_hex}/{block_name} has wrong node id")
        if int(block.attrs.get("block_id", -1)) != int(block_name):
            raise ValueError(f"raw block record {node_hex}/{block_name} has wrong block id")
        if int(block.attrs.get("histories", -1)) != FORMAL_BLOCK_HISTORIES:
            raise ValueError(
                f"formal raw block {node_hex}/{block_name} does not contain "
                f"{FORMAL_BLOCK_HISTORIES} histories"
            )
        if str(block.attrs.get("layout_json", "")) != str(
            root.attrs.get("layout_json", "")
        ):
            raise ValueError(
                f"raw block record {node_hex}/{block_name} does not match archive layout"
            )

    first = blocks[0]
    for block in blocks[1:]:
        _require_same_node(_BlockView(first), _BlockView(block), node_hex)

    generation = _fold_generation(root, node_hex)
    if generation is not None:
        if str(generation.attrs.get("state_sha256", "")) != _raw_record_sha256(
            generation, numpy, _FOLD_SEAL_EXCLUDED
        ):
            raise ValueError(f"fold state of node {node_hex} failed its checksum")
        rows = _ledger_rows(generation)
        if str(generation.attrs.get("ledger_chain_sha256", "")) != _ledger_chain_sha256(rows):
            raise ValueError(f"fold ledger of node {node_hex} failed its hash chain")
        if any(int(row["histories"]) != FORMAL_BLOCK_HISTORIES for row in rows):
            raise ValueError(
                f"fold ledger of node {node_hex} holds a block without "
                f"{FORMAL_BLOCK_HISTORIES} histories"
            )

    full, view = _fold_node(root, node_hex, numpy)
    expected = _node_values(full, view, node_hex, numpy)
    for path, wanted in expected.items():
        dataset = root[path]
        actual = dataset.asstr()[index] if hasattr(dataset, "asstr") and dataset.dtype.kind in {"O", "S", "U"} else dataset[index]
        if isinstance(wanted, str):
            same = str(actual) == wanted
        else:
            normalized = numpy.asarray(wanted, dtype=dataset.dtype)
            same = numpy.array_equal(numpy.asarray(actual), normalized)
        if not same:
            raise ValueError(
                f"materialized archive value {path}[{index}] does not close to raw blocks"
            )
    return set(expected)


def _require_same_node(first: _BlockView, block: _BlockView, node_hex: str) -> None:
    """Every block of a node must carry the node's exact deterministic values."""

    for name in _DETERMINISTIC_ATTRIBUTES:
        if _normal_hdf_value(block.attrs()[name]) != _normal_hdf_value(first.attrs()[name]):
            raise ValueError(f"raw blocks for node {node_hex} disagree on {name}")
    for name in _DETERMINISTIC_DATASETS:
        if not _same_array(block.array(name), first.array(name)):
            raise ValueError(f"raw blocks for node {node_hex} disagree on {name}")
    axes = first.names("physical_energy_edges_keV")
    if block.names("physical_energy_edges_keV") != axes:
        raise ValueError(f"raw blocks for node {node_hex} disagree on physical axes")
    for axis_name in axes:
        path = f"physical_energy_edges_keV/{axis_name}"
        if not _same_array(block.array(path), first.array(path)):
            raise ValueError(f"raw blocks for node {node_hex} disagree on physical axes")


def _same_array(left: Any, right: Any) -> bool:
    import numpy

    return numpy.array_equal(numpy.asarray(left), numpy.asarray(right))


def _prepare_node_row(root: Any, block: FitBlock) -> int:
    """The node's row index, appending and initializing the row for a new node."""

    index = _node_index(root, block.node_id)
    is_new = index is None
    if is_new:
        index = len(root["nodes/id"])
        _resize_node_axis(root, index + 1)
        root["nodes/id"][index] = block.node_id
        root["nodes/coordinates"][index] = (
            block.node["energy_keV"],
            block.node["distance_cm"],
            block.node["source_height_cm"],
            block.node["detector_height_cm"],
            block.node["extent_cm"],
        )
        root["nodes/q"][index] = block.node["q"]
        root["nodes/role"][index] = block.node["role"]
        root["nodes/direct_transmission"][index] = block.node["direct_transmission"]
        root["nodes/direct_fluence"][index] = block.node["direct_fluence"]
        root["nodes/target_total_rel_se"][index] = block.node["target_total_rel_se"]
        for axis_index, axis in enumerate(block.axes["energy"]):
            root[f"axes/energy/{axis_index}/physical_edges_keV"][index] = axis[
                "physical_edges_keV"
            ]
        root["components/support_mask"][index] = tuple(
            value for component in block.components for value in component["support_mask"]
        )
        line_energy = {line["name"]: line["energy_keV"] for line in block.axes["lines"]}
        root["lines/energy_keV"][index] = (
            line_energy.get("elastic_E0", block.node["energy_keV"]),
            line_energy.get("annihilation_511", 510.998950),
        )
        root["lines/support_mask"][index] = tuple(
            int(descriptor["name"] in line_energy) for descriptor in _CANONICAL_LINES
        )
    else:
        coordinates = tuple(float(x) for x in root["nodes/coordinates"][index])
        expected = (
            block.node["energy_keV"],
            block.node["distance_cm"],
            block.node["source_height_cm"],
            block.node["detector_height_cm"],
            block.node["extent_cm"],
        )
        if coordinates != expected or root["nodes/role"].asstr()[index] != block.node["role"]:
            raise ValueError(f"node metadata differs for existing node {block.node_id:016x}")
        for axis_index, axis in enumerate(block.axes["energy"]):
            if tuple(root[f"axes/energy/{axis_index}/physical_edges_keV"][index]) != axis[
                "physical_edges_keV"
            ]:
                raise ValueError(f"physical axes differ for existing node {block.node_id:016x}")
        stored_line_energy = tuple(float(value) for value in root["lines/energy_keV"][index])
        expected_line_energy = (
            block.node["energy_keV"],
            next(
                (
                    line["energy_keV"]
                    for line in block.axes["lines"]
                    if line["name"] == "annihilation_511"
                ),
                510.998950,
            ),
        )
        if stored_line_energy != expected_line_energy:
            raise ValueError(f"line energies differ for existing node {block.node_id:016x}")
        for name in ("direct_transmission", "direct_fluence", "target_total_rel_se"):
            if not math.isclose(
                float(root[f"nodes/{name}"][index]),
                float(block.node[name]),
                rel_tol=1.0e-13,
                abs_tol=0.0,
            ):
                raise ValueError(f"deterministic node value {name} changed across blocks")

    assert index is not None
    return index


def _open_node_fold(root: Any, node_hex: str, numpy: Any) -> dict[str, Any] | None:
    """The node's fold, ready to extend, or None while any of blocks
    1..RETAINED_RAW_BLOCKS is not yet archived (later blocks then stay raw).

    Generations other than the current one are left alone: they are dropped
    only after a verified ingest, so a generation that failed its post-write
    verification never costs the last verified one."""

    node_blocks = root["blocks"].get(node_hex)
    raw = {} if node_blocks is None else {int(name): group for name, group in node_blocks.items()}
    generation = _fold_generation(root, node_hex)
    if generation is not None:
        state = _load_fold_state(generation, numpy)
        _check_fold_consistency(node_hex, state, raw)
    elif all(block_id in raw for block_id in range(1, RETAINED_RAW_BLOCKS + 1)):
        state = _new_fold_state(
            _BlockView(raw[1]), root["diagnostics/interaction_counts"].ndim > 1, numpy
        )
        for block_id in range(1, RETAINED_RAW_BLOCKS + 1):
            _fold_view(state, _BlockView(raw[block_id]), numpy)
    else:
        return None
    fold = {
        "state": state,
        "first": _BlockView(raw[min(raw)]),
        "changed": False,
        # Raw records the fold covers beyond the retained blocks; removed only
        # after the generation that covers them has been verified.
        "covered": [
            block_id
            for block_id in sorted(raw)
            if RETAINED_RAW_BLOCKS < block_id <= int(state["through"])
        ],
    }
    _absorb_raw_tail(root, node_hex, fold, numpy)
    return fold


def _absorb_raw_tail(root: Any, node_hex: str, fold: dict[str, Any], numpy: Any) -> None:
    """Fold the raw records that directly follow the fold (blocks once held
    back by a gap in the block sequence)."""

    node_blocks = root["blocks"][node_hex]
    while str(int(fold["state"]["through"]) + 1) in node_blocks:
        block_id = int(fold["state"]["through"]) + 1
        view = _BlockView(node_blocks[str(block_id)])
        _require_same_node(fold["first"], view, node_hex)
        _fold_view(fold["state"], view, numpy)
        fold["changed"] = True
        if block_id > RETAINED_RAW_BLOCKS:
            fold["covered"].append(block_id)


def _ingest_one(
    root: Any, block: FitBlock, numpy: Any, *, rematerialize: bool = True
) -> None:
    """Store one block as a raw record (no folding)."""

    index = _prepare_node_row(root, block)
    _write_block_record(root, block, numpy)
    if rematerialize:
        _rematerialize_node(root, index, f"{block.node_id:016x}", numpy)


def ingest_blocks(
    block_paths: Iterable[str | Path],
    output_path: str | Path,
    *,
    delete_source_blocks: bool = False,
    fold_blocks: bool = True,
) -> dict[str, int]:
    """Validate then merge blocks into an appendable, chunked HDF5 archive.

    Blocks 1..RETAINED_RAW_BLOCKS of a node are stored as raw records.  With
    ``fold_blocks`` every later block that directly extends the node's fold
    is folded into it and keeps only a ledger row; a block beyond a gap in the
    node's block sequence is stored raw until the gap is filled.

    When ``delete_source_blocks`` is requested, deletion is transactional at
    the scientific-record level: every newly embedded record is re-hashed,
    every newly folded block's ledger row is checked, and every affected node
    is checked against its materialized aggregates after the archive has been
    closed and fsynced.  A failed check leaves every source worker block in
    place.  A replaced fold generation, and raw records a new fold covers,
    are removed only after that check.

    Memory: one in-memory fold state (up to ~140 MB in production) per node
    whose blocks interleave with another node's in the batch; run-streaming
    passes one node per batch.
    """

    h5py, numpy = require_hdf5()
    paths = tuple(Path(path) for path in block_paths)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    blocks_remaining = _preflight_block_paths(paths, output, h5py)
    source_records: list[tuple[Path, int, int, str, str, bool]] = []
    affected_nodes: set[int] = set()
    superseded: dict[str, dict[str, Any]] = {}
    new_archive = not output.exists()
    # Preserve a valid fast-resume inventory across the mutation. Missing or
    # stale sidecars are harmless and will be rebuilt by the streaming runner.
    # Remove a cached scientific identity before opening the file for mutation,
    # so no concurrent read can mistake an in-progress append for sealed state.
    from .storage import (
        cached_archive_inventory,
        invalidate_archive_seal,
        write_archive_inventory,
    )

    prior_inventory = set() if new_archive else cached_archive_inventory(output)
    invalidate_archive_seal(output)
    temporary: Path | None = None
    destination = output
    if new_archive:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{output.name}.ingest-",
            suffix=".tmp",
            dir=output.parent,
        )
        os.close(descriptor)
        temporary = Path(temporary_name)
        destination = temporary
    try:
        with h5py.File(
            destination, "w" if new_archive else "a", **(_FILE_SPACE if new_archive else {})
        ) as root:
            reuses_space = bool(root.id.get_create_plist().get_file_space_strategy()[1])
            folds: dict[int, dict[str, Any] | None] = {}
            for index, path in enumerate(paths):
                block = read_fit_block(path)
                if new_archive and index == 0:
                    _initialize_archive(root, block, h5py, numpy)
                _check_against_archive(root, block)
                _prepare_node_row(root, block)
                node_hex = f"{block.node_id:016x}"
                record = _block_record(block, numpy)
                folded = False
                if fold_blocks and block.block_id > RETAINED_RAW_BLOCKS:
                    fold = folds.get(block.node_id)
                    if fold is None:
                        # Opened again while it is None: blocks 1..RETAINED_RAW_BLOCKS
                        # of a new node arrive earlier in the same batch.
                        fold = folds[block.node_id] = _open_node_fold(root, node_hex, numpy)
                    if fold is not None and block.block_id == int(fold["state"]["through"]) + 1:
                        view = _BlockView(record)
                        _require_same_node(fold["first"], view, node_hex)
                        _fold_view(fold["state"], view, numpy)
                        fold["changed"] = True
                        _absorb_raw_tail(root, node_hex, fold, numpy)
                        folded = True
                if not folded:
                    _write_block_record(root, block, numpy, record)
                source_records.append(
                    (
                        path,
                        block.node_id,
                        block.block_id,
                        block.payload_sha256,
                        str(record.attrs[""]["scientific_record_sha256"]),
                        folded,
                    )
                )
                affected_nodes.add(block.node_id)
                blocks_remaining[block.node_id] -= 1
                if blocks_remaining[block.node_id] == 0:
                    node_index = _node_index(root, block.node_id)
                    assert node_index is not None
                    fold = folds.pop(block.node_id, None)
                    if fold is None and fold_blocks and any(
                        int(name) > RETAINED_RAW_BLOCKS for name in root[f"blocks/{node_hex}"]
                    ):
                        # Blocks beyond the fourth that arrived before blocks
                        # 1..RETAINED_RAW_BLOCKS were complete, or before a gap
                        # below them was filled.
                        fold = _open_node_fold(root, node_hex, numpy)
                    state = None
                    if fold is not None and fold["changed"]:
                        previous = _store_fold_generation(root, node_hex, fold["state"], numpy)
                        superseded[node_hex] = {"generation": previous, "raw": list(fold["covered"])}
                        state = fold["state"]
                    _rematerialize_node(root, node_index, node_hex, numpy, state=state)
                del block, record
            root.flush()
            node_count = len(root["nodes/id"])
            config_sha256 = str(root["manifest"].attrs["config_sha256"])
        if superseded and not reuses_space:
            import sys

            print(
                f"warning: {output} predates persistent free-space tracking, so the "
                "space of replaced fold states is not reused; a new archive avoids it",
                file=sys.stderr,
            )
        if temporary is not None:
            os.replace(temporary, output)
            temporary = None
        _fsync_path(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

    sources_deleted = 0
    if delete_source_blocks or superseded:
        # Reopen after the durable close.  This is intentionally narrower than
        # a publication audit: only the newly embedded records, the newly
        # folded ledger rows and the affected nodes need to close before
        # redundant worker files or superseded fold generations are removed.
        with h5py.File(output, "r") as root:
            inventories: dict[int, dict[int, dict[str, Any]]] = {}
            for _, node_id, block_id, payload_sha256, seal, folded in source_records:
                if folded:
                    if node_id not in inventories:
                        inventories[node_id] = node_block_inventory(root, f"{node_id:016x}")
                    entry = inventories[node_id].get(block_id)
                    if (
                        entry is None
                        or not entry["folded"]
                        or entry["payload_sha256"] != payload_sha256
                        or entry["scientific_record_sha256"] != seal
                    ):
                        raise ValueError(
                            f"folded block {node_id:016x}/{block_id} has no matching "
                            "fold-ledger row"
                        )
                    continue
                group = root[f"blocks/{node_id:016x}/{block_id}"]
                if str(group.attrs.get("payload_sha256", "")) != payload_sha256:
                    raise ValueError(
                        f"embedded block {node_id:016x}/{block_id} payload identity "
                        "differs from its worker block"
                    )
                stored = str(group.attrs.get("scientific_record_sha256", ""))
                if stored != seal or stored != _raw_record_sha256(group, numpy):
                    raise ValueError(
                        f"embedded block {node_id:016x}/{block_id} failed its "
                        "scientific checksum"
                    )
            for node_id in sorted(affected_nodes):
                node_index = _node_index(root, node_id)
                assert node_index is not None
                _verify_materialized_node(
                    root, node_index, f"{node_id:016x}", numpy
                )
        if superseded:
            with h5py.File(output, "a") as root:
                for node_hex, item in superseded.items():
                    _drop_stale_fold_generations(root, node_hex)
                    for block_id in item["raw"]:
                        del root[f"blocks/{node_hex}/{block_id}"]
            _fsync_path(output)

    if prior_inventory is not None:
        write_archive_inventory(
            output,
            prior_inventory
            | {(node_id, block_id) for _, node_id, block_id, _, _, _ in source_records},
            config_sha256,
        )

    if delete_source_blocks:
        source_directories: set[Path] = set()
        for path, _, _, _, _, _ in source_records:
            path.unlink()
            source_directories.add(path.parent)
            sources_deleted += 1
        for directory in source_directories:
            try:
                descriptor = os.open(directory, os.O_RDONLY)
                try:
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
            except OSError:
                pass
    return {
        "blocks_ingested": len(paths),
        "blocks_folded": sum(1 for record in source_records if record[5]),
        "nodes_total": node_count,
        "sources_deleted": sources_deleted,
    }


def _fsync_path(output: Path) -> None:
    descriptor = os.open(output, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    try:
        directory_descriptor = os.open(output.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except OSError:
        # Some removable filesystems do not implement directory fsync.
        pass


def _node_major_paths(root: Any) -> set[str]:
    """Every materialized dataset with one row per node (the fitter's view)."""

    paths: set[str] = set()

    def collect(group: Any, prefix: str) -> None:
        for name, value in group.items():
            path = f"{prefix}/{name}" if prefix else str(name)
            if path.split("/", 1)[0] in {"blocks", "folds"}:
                continue
            if hasattr(value, "items"):
                collect(value, path)
            elif (
                value.maxshape is not None
                and len(value.shape) >= 1
                and value.maxshape[0] is None
            ):
                paths.add(path)

    collect(root, "")
    return paths


def node_closure_identity(
    root: Any, index: int, node_hex: str, node_major_paths: set[str], numpy: Any
) -> str:
    """Cheap identity of what one node's closure covers.

    Hashes the node's block inventory as its record seals and fold ledger
    describe it (block ids, payload and scientific-record checksums), its
    current fold generation and that generation's seal, and every value of
    its materialized row, without decompressing raw records or fold states.
    Equal identities mean neither the node's blocks nor its fitter-facing
    values have changed since a closure that produced the identity, so an
    incremental seal may skip re-verifying that node.
    """

    digest = hashlib.sha256()
    _hash_frame(digest, "node", node_hex.encode("ascii"))
    _hash_frame(digest, "paths", _canonical(sorted(node_major_paths)).encode("utf-8"))
    for block_id, entry in node_block_inventory(root, node_hex).items():
        _hash_frame(
            digest,
            f"block:{block_id}",
            _canonical(
                [entry["payload_sha256"], entry["scientific_record_sha256"]]
            ).encode("utf-8"),
        )
    generation = _fold_generation(root, node_hex)
    if generation is not None:
        _hash_frame(
            digest,
            "fold",
            _canonical(
                [
                    generation.name.rsplit("/", 1)[-1],
                    str(generation.attrs.get("state_sha256", "")),
                ]
            ).encode("utf-8"),
        )
    for path in sorted(node_major_paths):
        dataset = root[path]
        if dataset.dtype.kind in {"O", "S", "U"}:
            payload = str(dataset.asstr()[index]).encode("utf-8")
        else:
            values = numpy.ascontiguousarray(numpy.asarray(dataset[index]))
            payload = (
                _canonical([list(values.shape), values.dtype.str]).encode("ascii")
                + values.tobytes()
            )
        _hash_frame(digest, f"row:{path}", payload)
    return digest.hexdigest()


def verify_archive_materialization(
    root: Any,
    *,
    trusted_nodes: Mapping[str, str] | None = None,
    identities_out: dict[str, str] | None = None,
) -> None:
    """Fail unless every materialized scientific value closes to raw blocks.

    This function is deliberately read-only.  It is the publication boundary:
    corrupt or stale aggregates are never silently assigned a fresh snapshot
    identity.  Raw records without the complete normalized metadata fail with
    a re-ingestion instruction instead of receiving weaker trust.

    ``trusted_nodes`` maps node ids to the node_closure_identity() a previous
    closure verified; a node whose identity is unchanged is not re-verified
    (an incremental seal).  Every other node, and every archive-wide check,
    runs in full.  ``identities_out`` receives every node's identity.
    """

    _, numpy = require_hdf5()
    if str(root.attrs.get("schema", "")) != ARCHIVE_SCHEMA:
        raise ValueError("output is not a supported CEEFIT HDF5 archive")
    layout_text = str(root.attrs.get("layout_json", ""))
    try:
        layout = json.loads(layout_text)
    except json.JSONDecodeError as error:
        raise ValueError("fit archive has invalid canonical layout metadata") from error
    if _canonical(layout) != layout_text:
        raise ValueError("fit archive layout metadata is not canonical")

    block_root = root.get("blocks")
    if block_root is None or not len(block_root):
        raise ValueError("fit archive contains no raw block records")
    node_ids = [int(value) for value in root["nodes/id"][...]]
    if len(set(node_ids)) != len(node_ids):
        raise ValueError("fit archive contains duplicate materialized node ids")
    raw_ids: list[int] = []
    for node_name in block_root:
        try:
            node_id = int(node_name, 16)
        except ValueError as error:
            raise ValueError(f"invalid raw-block node id {node_name!r}") from error
        if f"{node_id:016x}" != node_name:
            raise ValueError(f"raw-block node id {node_name!r} is not canonical")
        raw_ids.append(node_id)
    if len(set(raw_ids)) != len(raw_ids) or set(raw_ids) != set(node_ids):
        raise ValueError("fit archive raw-block and materialized node inventories differ")
    for node_name in root.get("folds", {}):
        if node_name not in block_root:
            raise ValueError(f"fold state {node_name!r} belongs to no archived node")

    def require_attribute(group: Any, name: str, wanted: Any, label: str) -> None:
        if name not in group.attrs or _normal_hdf_value(group.attrs[name]) != _normal_hdf_value(wanted):
            raise ValueError(f"materialized archive metadata {label} does not close to raw blocks")

    def require_dataset(path: str, wanted: Any) -> None:
        if path not in root:
            raise ValueError(f"materialized archive lacks required dataset {path}")
        dataset = root[path]
        actual = numpy.asarray(dataset[...])
        expected = numpy.asarray(wanted, dtype=dataset.dtype)
        if not numpy.array_equal(actual, expected):
            raise ValueError(
                f"materialized archive dataset {path} does not close to raw blocks"
            )

    # Static archive metadata also comes from the per-block canonical layout.
    manifest_record = layout.get("manifest")
    if not isinstance(manifest_record, dict):
        raise ValueError("fit archive layout lacks its worker manifest")
    try:
        canonical_configuration = json.loads(manifest_record["canonical_config_json"])
        configured_block_histories = int(canonical_configuration["block_histories"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError("fit archive lacks a canonical block-size contract") from error
    if configured_block_histories != FORMAL_BLOCK_HISTORIES:
        raise ValueError(
            f"formal fit archives require {FORMAL_BLOCK_HISTORIES} histories per block"
        )
    for name, wanted in manifest_record.items():
        if not isinstance(wanted, (str, bool, int, float)):
            continue
        require_attribute(root["manifest"], name, wanted, f"manifest/{name}")

    require_attribute(
        root["axes"],
        "total_energy_axis",
        layout["total_energy_axis"],
        "axes/total_energy_axis",
    )
    require_attribute(
        root["axes"],
        "line_definitions_json",
        _canonical(layout["lines"]),
        "axes/line_definitions_json",
    )
    require_dataset("axes/mu_edges", layout["mu_edges"])
    energy_axes = layout["energy_axes"]
    if set(root["axes/energy"]) != {str(index) for index in range(len(energy_axes))}:
        raise ValueError("materialized archive energy-axis inventory is inconsistent")
    for index, axis in enumerate(energy_axes):
        group = root[f"axes/energy/{index}"]
        require_attribute(group, "name", axis["name"], f"axes/energy/{index}/name")
        require_attribute(
            group,
            "transform_version",
            axis["transform_version"],
            f"axes/energy/{index}/transform_version",
        )
        require_dataset(
            f"axes/energy/{index}/registered_edges", axis["registered_edges"]
        )
    require_attribute(
        root["components"],
        "descriptors_json",
        _canonical(layout["components"]),
        "components/descriptors_json",
    )
    require_attribute(
        root["diagnostics"],
        "interaction_names_json",
        _canonical(layout["interaction_names"]),
        "diagnostics/interaction_names_json",
    )
    stopping_definition_json = layout.get("stopping_definition_json")
    if not isinstance(stopping_definition_json, str) or not stopping_definition_json:
        raise ValueError("fit archive layout lacks its stopping-vector definition")
    require_attribute(
        root["stopping"],
        "definition_json",
        stopping_definition_json,
        "stopping/definition_json",
    )

    first_node = f"{min(raw_ids):016x}"
    first_block_name = min(root[f"blocks/{first_node}"], key=int)
    first_block = root[f"blocks/{first_node}/{first_block_name}"]
    require_attribute(
        first_block,
        "stopping_definition_json",
        stopping_definition_json,
        "raw-block/stopping_definition_json",
    )
    if "component_support_offsets" not in first_block:
        raise ValueError(
            "fit archive raw records lack component support offsets; "
            "re-ingest the original CEEFIT blocks into a new archive"
        )
    require_dataset(
        "components/support_offsets", first_block["component_support_offsets"][...]
    )

    node_major_paths = _node_major_paths(root)
    want_identity = trusted_nodes is not None or identities_out is not None
    for node_id in sorted(node_ids):
        index = node_ids.index(node_id)
        node_hex = f"{node_id:016x}"
        identity = (
            node_closure_identity(root, index, node_hex, node_major_paths, numpy)
            if want_identity
            else None
        )
        if identities_out is not None:
            identities_out[node_hex] = identity
        if trusted_nodes is not None and trusted_nodes.get(node_hex) == identity:
            continue
        paths = _verify_materialized_node(root, index, node_hex, numpy)
        if paths != node_major_paths:
            missing = sorted(node_major_paths - paths)
            extra = sorted(paths - node_major_paths)
            raise ValueError(
                "raw-block closure does not cover every materialized node dataset "
                f"(unverified={missing}, absent={extra})"
            )
