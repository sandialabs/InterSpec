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

"""Strict decoder and validator for the dependency-free CEEFIT block format."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import struct
from typing import Any, BinaryIO

from .axes import (
    annihilation_registered_axis,
    brems_registered_axis,
    primary_registered_axis,
    uniform_mu_edges,
)
from .contracts import FORMAL_BLOCK_HISTORIES, stopping_vector_definition_json
from .design import PhysicalNode, geometry_transform, stable_node_id_from_cm


MAGIC = b"CEEFIT-BLOCK-v1\0"
VERSION = 1
HEADER = struct.Struct("<16sIIQ32s")
MAX_PAYLOAD_BYTES = 8 * 1024 * 1024 * 1024
MAX_VECTOR_ITEMS = 1 << 31
_CANONICAL_AXIS_NAMES = (
    "primary_registered",
    "annihilation_511_registered",
    "brems_log_relative",
)
_CANONICAL_COMPONENTS = tuple(
    (f"{origin_name}_{order_name}", origin, order, origin)
    for origin, origin_name in enumerate(("primary", "annihilation", "brems"))
    for order, order_name in enumerate(("order0", "order1", "order2plus"))
)
_PAIR_THRESHOLD_KEV = 2.0 * 510.998950
_CANONICAL_DIAGNOSTIC_CHANNELS = (
    "air_only",
    "one_slab_bounce",
    "multi_same_slab",
    "multi_distinct_slab",
)
_CANONICAL_INTERACTIONS = (
    "photoelectric",
    "compton",
    "rayleigh",
    "pair_production",
    "brems_emission",
)
_RNG_ALGORITHM = "mt19937_64/sha256-tuple-stream-v1"
_UNITS = "energy=keV,length=cm,fluence=cm^-2"
_OBSERVABLE = "bin-integrated scalar incident fluence per source photon"
_NORMALIZATION = "raw S1/S2; s_hat divisor=T(E0)/(4*pi*d^2)"


class BlockFormatError(ValueError):
    """The block is corrupt, unsupported, or structurally inconsistent."""


@dataclass(frozen=True)
class FitBlock:
    node_id: int
    block_id: int
    manifest: dict[str, Any]
    axes: dict[str, Any]
    components: tuple[dict[str, Any], ...]
    node: dict[str, Any]
    replicas: tuple[dict[str, int], ...]
    moments: dict[str, dict[str, tuple[Any, ...]]]
    replica_s1: dict[str, tuple[float, ...]]
    stopping: dict[str, Any]
    diagnostics: dict[str, Any]
    payload_sha256: str
    source_path: str | None = None


class _Reader:
    def __init__(self, payload: bytes):
        self.stream = io.BytesIO(payload)
        self.size = len(payload)

    @property
    def offset(self) -> int:
        return self.stream.tell()

    def _unpack(self, fmt: str) -> tuple[Any, ...]:
        spec = struct.Struct("<" + fmt)
        raw = self.stream.read(spec.size)
        if len(raw) != spec.size:
            raise BlockFormatError(f"truncated payload at byte {self.offset}")
        return spec.unpack(raw)

    def u8(self) -> int:
        return int(self._unpack("B")[0])

    def boolean(self, name: str) -> bool:
        value = self.u8()
        if value not in (0, 1):
            raise BlockFormatError(f"{name} is not encoded as canonical 0/1")
        return bool(value)

    def u32(self) -> int:
        return int(self._unpack("I")[0])

    def u64(self) -> int:
        return int(self._unpack("Q")[0])

    def f64(self) -> float:
        return float(self._unpack("d")[0])

    def string(self) -> str:
        count = self.u64()
        if count > self.size - self.offset:
            raise BlockFormatError("string length exceeds remaining payload")
        raw = self.stream.read(count)
        try:
            return raw.decode("utf-8")
        except UnicodeDecodeError as error:
            raise BlockFormatError("block contains invalid UTF-8") from error

    def vector(self, code: str) -> tuple[Any, ...]:
        count = self.u64()
        if count > MAX_VECTOR_ITEMS:
            raise BlockFormatError(f"vector length {count} exceeds safety limit")
        spec = struct.Struct("<" + code)
        byte_count = count * spec.size
        if byte_count > self.size - self.offset:
            raise BlockFormatError("vector length exceeds remaining payload")
        raw = self.stream.read(byte_count)
        if not raw:
            return ()
        # iter_unpack avoids constructing a gigantic format string.
        return tuple(value[0] for value in struct.iter_unpack("<" + code, raw))

    def done(self) -> None:
        if self.offset != self.size:
            raise BlockFormatError(f"payload has {self.size - self.offset} trailing bytes")


def _read_moments(reader: _Reader) -> dict[str, tuple[Any, ...]]:
    return {
        "N": reader.vector("Q"),
        "S1": reader.vector("d"),
        "S2": reader.vector("d"),
    }


def _validate_moments(name: str, moments: dict[str, tuple[Any, ...]], expected: int) -> None:
    for field in ("N", "S1", "S2"):
        if len(moments[field]) != expected:
            raise BlockFormatError(
                f"{name}.{field} has {len(moments[field])} values; expected {expected}"
            )


def _validate_numeric_moments(
    name: str,
    moments: dict[str, tuple[Any, ...]],
    histories: int,
) -> None:
    for index, (count, first, second) in enumerate(
        zip(moments["N"], moments["S1"], moments["S2"])
    ):
        if count != histories:
            raise BlockFormatError(f"{name}.N[{index}] does not equal node histories")
        if not math.isfinite(first) or not math.isfinite(second):
            raise BlockFormatError(f"{name} contains a non-finite sufficient statistic")
        if first < 0.0 or second < 0.0:
            raise BlockFormatError(f"{name} contains a negative fluence moment")
        minimum = first * first / histories
        tolerance = 2.0e-12 * max(1.0, abs(minimum), abs(second))
        if second + tolerance < minimum:
            raise BlockFormatError(f"{name}.S2[{index}] violates S2 >= S1^2/N")


def _close(actual: float, expected: float, *, float32: bool = False) -> bool:
    relative = 2.0e-5 if float32 else 2.0e-11
    return math.isclose(actual, expected, rel_tol=relative, abs_tol=1.0e-12)


def _valid_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _validate_replica_closure(
    name: str,
    values: tuple[float, ...],
    first_moments: tuple[Any, ...],
    replicas: int,
) -> None:
    cells = len(first_moments)
    for value in values:
        if not math.isfinite(value) or value < 0.0:
            raise BlockFormatError(f"replica {name} contains non-finite/negative fluence")
    for cell, expected in enumerate(first_moments):
        actual = math.fsum(values[replica * cells + cell] for replica in range(replicas))
        if not _close(actual, float(expected), float32=True):
            raise BlockFormatError(f"replica {name} does not close to block S1 at cell {cell}")


def _decode_payload(payload: bytes, payload_sha256: str, source_path: str | None) -> FitBlock:
    reader = _Reader(payload)
    node_id, block_id = reader.u64(), reader.u64()
    manifest = {
        "canonical_config_json": reader.string(),
        "config_sha256": reader.string(),
    }
    for name in (
        "material_name",
        "material_definition_json",
        "geometry_tier",
        "rng_algorithm",
        "estimator_flags",
        "units",
        "observable",
        "normalization",
    ):
        manifest[name] = reader.string()

    energy_axes = []
    for _ in range(reader.u64()):
        energy_axes.append(
            {
                "name": reader.string(),
                "transform_version": reader.u32(),
                "registered_edges": reader.vector("d"),
                "physical_edges_keV": reader.vector("d"),
            }
        )
    total_energy_axis = reader.u32()
    mu_edges = reader.vector("d")
    lines = []
    for _ in range(reader.u64()):
        lines.append(
            {
                "name": reader.string(),
                "energy_keV": reader.f64(),
                "origin": reader.u32(),
                "post_creation_order": reader.u32(),
            }
        )
    axes = {
        "energy": tuple(energy_axes),
        "total_energy_axis": total_energy_axis,
        "mu_edges": mu_edges,
        "lines": tuple(lines),
    }

    components = []
    for _ in range(reader.u64()):
        components.append(
            {
                "name": reader.string(),
                "origin": reader.u32(),
                "post_creation_order": reader.u32(),
                "energy_axis": reader.u32(),
                "support_mask": reader.vector("B"),
            }
        )

    coordinates = [reader.f64() for _ in range(5)]
    transformed_q = tuple(reader.f64() for _ in range(4))
    node = {
        "energy_keV": coordinates[0],
        "distance_cm": coordinates[1],
        "source_height_cm": coordinates[2],
        "detector_height_cm": coordinates[3],
        "extent_cm": coordinates[4],
        "q": transformed_q,
        "role": reader.string(),
        "status": reader.string(),
        "direct_transmission": reader.f64(),
        "direct_fluence": reader.f64(),
        "histories": reader.u64(),
        "wall_seconds": reader.f64(),
        "cpu_seconds": reader.f64(),
        "target_total_rel_se": reader.f64(),
        "achieved_total_rel_se": reader.f64(),
    }

    replica_count = reader.u32()
    replicas = tuple(
        {
            "id": reader.u32(),
            "history_begin": reader.u64(),
            "history_count": reader.u64(),
            "rng_seed": reader.u64(),
        }
        for _ in range(replica_count)
    )
    moments = {
        "continuum_components": _read_moments(reader),
        "continuum_total": _read_moments(reader),
        "line_components": _read_moments(reader),
        "line_total": _read_moments(reader),
    }
    replica_s1 = {
        "continuum_components": reader.vector("f"),
        "continuum_total": reader.vector("f"),
        "line_components": reader.vector("f"),
        "line_total": reader.vector("f"),
    }
    stopping_definition_json = reader.string()
    stopping_n = reader.u64()
    stopping = {
        "definition_json": stopping_definition_json,
        "N": stopping_n,
        "S1": tuple(reader.f64() for _ in range(16)),
        "outer_S2": tuple(reader.f64() for _ in range(16 * 16)),
    }

    diagnostic_channels = []
    for _ in range(reader.u64()):
        diagnostic_channels.append(
            {
                "name": reader.string(),
                "continuum": _read_moments(reader),
                "lines": _read_moments(reader),
            }
        )
    interaction_names = tuple(reader.string() for _ in range(reader.u64()))
    diagnostics = {
        "channels": tuple(diagnostic_channels),
        "interaction_names": interaction_names,
        "interaction_counts": reader.vector("Q"),
        "scattered_total": {"N": reader.u64(), "S1": reader.f64(), "S2": reader.f64()},
        "dropped": {"N": reader.u64(), "S1": reader.f64(), "S2": reader.f64()},
        "transport_failures": reader.u64(),
        "max_history_contribution": reader.f64(),
        "capped": reader.boolean("diagnostics.capped"),
        "failure_reason": reader.string(),
    }
    reader.done()

    if not _valid_sha256(manifest["config_sha256"]):
        raise BlockFormatError("manifest config_sha256 is not lowercase SHA-256")
    required_manifest_fields = (
        "material_name",
        "material_definition_json",
        "geometry_tier",
        "rng_algorithm",
        "estimator_flags",
        "units",
        "observable",
        "normalization",
    )
    if any(not manifest[name] for name in required_manifest_fields):
        raise BlockFormatError(
            "manifest is missing a required material, estimator, or observable field"
        )

    if replica_count != 32:
        raise BlockFormatError(f"replica_count is {replica_count}; CEEFIT-v1 requires 32")
    replica_ids = [item["id"] for item in replicas]
    if replica_ids != list(range(32)):
        raise BlockFormatError("replica IDs must be exactly 0..31 in order")
    replica_counts = {item["history_count"] for item in replicas}
    if len(replica_counts) != 1 or 0 in replica_counts:
        raise BlockFormatError("replica history counts must be equal and nonzero")
    if len({item["rng_seed"] for item in replicas}) != replica_count:
        raise BlockFormatError("replica RNG seeds must be unique")
    ranges = sorted(
        (item["history_begin"], item["history_begin"] + item["history_count"])
        for item in replicas
    )
    if any(start < previous_end for (_, previous_end), (start, _) in zip(ranges, ranges[1:])):
        raise BlockFormatError("replica history ranges overlap")
    if sum(item["history_count"] for item in replicas) != node["histories"]:
        raise BlockFormatError("replica history counts do not sum to node histories")
    for replica in replicas:
        seed_record = (
            b"CEEFIT-RNG-v1\0"
            + manifest["config_sha256"].encode("ascii")
            + struct.pack("<QQI", node_id, block_id, replica["id"])
        )
        expected_seed = int.from_bytes(hashlib.sha256(seed_record).digest()[:8], "big")
        if replica["rng_seed"] != expected_seed:
            raise BlockFormatError(f"replica {replica['id']} RNG seed violates the CEEFIT identity")

    canonical = manifest["canonical_config_json"].encode("utf-8")
    if hashlib.sha256(canonical).hexdigest() != manifest["config_sha256"]:
        raise BlockFormatError("canonical configuration SHA-256 mismatch")
    try:
        configuration = json.loads(manifest["canonical_config_json"])
        physics_configuration = configuration["physics"]
        if any(
            not isinstance(physics_configuration[name], bool)
            for name in ("include_annihilation", "include_brems", "include_rayleigh")
        ):
            raise TypeError("physics switches must be JSON booleans")
        configuration_schema = str(configuration["schema"])
        configuration_tier = str(configuration["geometry_tier"])
        configured_energy_bins = int(configuration["axis"]["energy_bins"])
        configured_e_cut = float(configuration["axis"]["energy_cut_keV"])
        configured_mu_bins = int(configuration["axis"]["mu_bins"])
        configured_transform = str(configuration["axis"]["transform"])
        configured_block_histories = int(configuration["block_histories"])
        brems_floor = float(physics_configuration["brems_electron_floor_keV"])
        include_annihilation = physics_configuration["include_annihilation"]
        include_brems = physics_configuration["include_brems"]
        include_rayleigh = physics_configuration["include_rayleigh"]
        material_definition = (
            configuration["air"]
            if configuration_tier == "air"
            else configuration["material"]
        )
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise BlockFormatError("canonical configuration lacks the CEEFIT axis/physics record") from error
    if (
        configuration_schema != "ceelo-env-fit-config-v1"
        or configuration_tier not in {"air", "field"}
        or configured_transform != "registered-component-v1"
        or configured_energy_bins <= 0
        or configured_energy_bins % 2
        or configured_e_cut <= 0.0
        or configured_mu_bins not in (32, 64)
        or configured_block_histories != FORMAL_BLOCK_HISTORIES
        or brems_floor < 0.0
        or not include_annihilation
        or not include_brems
        or not include_rayleigh
    ):
        raise BlockFormatError(
            "canonical configuration has invalid fit axes/threshold/block size"
        )
    try:
        manifest_material = json.loads(manifest["material_definition_json"])
    except json.JSONDecodeError as error:
        raise BlockFormatError("manifest material definition is not valid JSON") from error
    expected_estimator_flags = (
        "next-event;rayleigh="
        + ("on" if include_rayleigh else "off")
        + ";brems="
        + ("on" if include_brems else "off")
        + ";annihilation="
        + ("on" if include_annihilation else "off")
    )
    if (
        manifest["geometry_tier"] != configuration_tier
        or manifest_material != material_definition
        or not isinstance(material_definition, dict)
        or manifest["material_name"] != material_definition.get("name")
        or manifest["rng_algorithm"] != _RNG_ALGORITHM
        or manifest["estimator_flags"] != expected_estimator_flags
        or manifest["units"] != _UNITS
        or manifest["observable"] != _OBSERVABLE
        or manifest["normalization"] != _NORMALIZATION
    ):
        raise BlockFormatError("manifest fields disagree with the canonical physics configuration")
    if node["histories"] != configured_block_histories:
        raise BlockFormatError("node histories disagree with canonical block_histories")
    histories_per_replica = configured_block_histories // 32
    for replica in replicas:
        expected_begin = (
            (block_id - 1) * configured_block_histories
            + replica["id"] * histories_per_replica
        )
        if replica["history_count"] != histories_per_replica or replica["history_begin"] != expected_begin:
            raise BlockFormatError("replica history count/range violates canonical block slots")
    if len(energy_axes) != 3 or tuple(axis["name"] for axis in energy_axes) != _CANONICAL_AXIS_NAMES:
        raise BlockFormatError("energy-axis inventory is not the canonical three-component taxonomy")
    if total_energy_axis != 0:
        raise BlockFormatError("the common total continuum must use the primary registered axis")
    expected_axis_records = (
        primary_registered_axis(
            node["energy_keV"],
            e_cut_keV=configured_e_cut,
            n_per_block=configured_energy_bins // 2,
            n_mu=configured_mu_bins,
        ),
        annihilation_registered_axis(
            e_cut_keV=configured_e_cut,
            n_per_block=configured_energy_bins // 2,
            n_mu=configured_mu_bins,
        ),
        brems_registered_axis(
            node["energy_keV"],
            e_cut_keV=configured_e_cut,
            n_energy=configured_energy_bins,
            n_mu=configured_mu_bins,
        ),
    )
    for axis, expected_axis in zip(energy_axes, expected_axis_records):
        registered = axis["registered_edges"]
        physical = axis["physical_edges_keV"]
        if axis["transform_version"] != 1:
            raise BlockFormatError(f"axis {axis['name']!r} has an unsupported transform")
        if len(registered) != configured_energy_bins + 1 or len(physical) != len(registered):
            raise BlockFormatError(f"axis {axis['name']!r} disagrees with configured energy bins")
        if (
            registered[0] != 0.0
            or registered[-1] != 1.0
            or any(not math.isfinite(value) for value in (*registered, *physical))
            or any(right <= left for left, right in zip(registered, registered[1:]))
            or any(right <= left for left, right in zip(physical, physical[1:]))
            or physical[0] <= 0.0
        ):
            raise BlockFormatError(f"axis {axis['name']!r} has invalid ordered edges")
        for actual, expected in zip(registered, expected_axis.registered_edges):
            if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=2.0e-15):
                raise BlockFormatError(f"axis {axis['name']!r} registered edges are noncanonical")
        for actual, expected in zip(physical, expected_axis.physical_energy_edges_keV):
            if not math.isclose(actual, expected, rel_tol=2.0e-13, abs_tol=1.0e-12):
                raise BlockFormatError(f"axis {axis['name']!r} physical transform is noncanonical")
    if total_energy_axis >= len(energy_axes):
        raise BlockFormatError("total energy-axis index is out of range")
    expected_stopping_definition = stopping_vector_definition_json(
        total_energy_axis,
        energy_axes[total_energy_axis]["name"],
        configured_energy_bins,
    )
    if stopping["definition_json"] != expected_stopping_definition:
        raise BlockFormatError(
            "stopping-vector definition is absent, malformed, or noncanonical"
        )
    if (
        len(mu_edges) != configured_mu_bins + 1
        or mu_edges[0] != -1.0
        or mu_edges[-1] != 1.0
        or any(not math.isfinite(value) for value in mu_edges)
        or any(right <= left for left, right in zip(mu_edges, mu_edges[1:]))
    ):
        raise BlockFormatError("mu axis disagrees with the canonical ordered configuration")
    if any(
        not math.isclose(actual, expected, rel_tol=0.0, abs_tol=2.0e-15)
        for actual, expected in zip(mu_edges, uniform_mu_edges(configured_mu_bins))
    ):
        raise BlockFormatError("mu edges are not the canonical uniform partition")
    n_mu = len(mu_edges) - 1
    component_inventory = tuple(
        (
            component["name"],
            component["origin"],
            component["post_creation_order"],
            component["energy_axis"],
        )
        for component in components
    )
    if component_inventory != _CANONICAL_COMPONENTS:
        raise BlockFormatError("component inventory is not the canonical 3x3 taxonomy")
    total_cells = (len(energy_axes[total_energy_axis]["physical_edges_keV"]) - 1) * n_mu
    component_cells = 0
    for component_index, component in enumerate(components):
        axis_index = component["energy_axis"]
        if axis_index >= len(energy_axes):
            raise BlockFormatError(f"component {component['name']!r} axis is out of range")
        cells = (len(energy_axes[axis_index]["physical_edges_keV"]) - 1) * n_mu
        energy_bins = len(energy_axes[axis_index]["physical_edges_keV"]) - 1
        if len(component["support_mask"]) != energy_bins:
            raise BlockFormatError(f"component {component['name']!r} support mask has wrong size")
        if any(value not in (0, 1) for value in component["support_mask"]):
            raise BlockFormatError(f"component {component['name']!r} support mask is not binary")
        origin, order = divmod(component_index, 3)
        expected_support = (
            order != 0
            if origin == 0
            else node["energy_keV"] >= _PAIR_THRESHOLD_KEV and order != 0
            if origin == 1
            else node["energy_keV"] >= brems_floor
        )
        if any(bool(value) != expected_support for value in component["support_mask"]):
            raise BlockFormatError(f"component {component['name']!r} support violates physics threshold")
        component_cells += cells
    expected_lines = [
        ("elastic_E0", node["energy_keV"], 0, 3),
    ]
    if node["energy_keV"] >= _PAIR_THRESHOLD_KEV:
        expected_lines.append(("annihilation_511", 510.998950, 1, 3))
    actual_lines = [
        (line["name"], line["energy_keV"], line["origin"], line["post_creation_order"])
        for line in lines
    ]
    if len(actual_lines) != len(expected_lines):
        raise BlockFormatError("exact-line inventory disagrees with the source-energy threshold")
    for actual, expected in zip(actual_lines, expected_lines):
        if actual[0] != expected[0] or actual[2:] != expected[2:] or not math.isclose(
            actual[1], expected[1], rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise BlockFormatError("exact-line taxonomy/energy is noncanonical")
    line_cells = len(lines) * n_mu
    line_component_cells = len(components) * line_cells
    _validate_moments("continuum_components", moments["continuum_components"], component_cells)
    _validate_moments("continuum_total", moments["continuum_total"], total_cells)
    _validate_moments("line_components", moments["line_components"], line_component_cells)
    _validate_moments("line_total", moments["line_total"], line_cells)
    for name, value in moments.items():
        _validate_numeric_moments(name, value, node["histories"])
    for name, base in (
        ("continuum_components", component_cells),
        ("continuum_total", total_cells),
        ("line_components", line_component_cells),
        ("line_total", line_cells),
    ):
        if len(replica_s1[name]) != replica_count * base:
            raise BlockFormatError(f"replica {name} tensor has wrong size")
        _validate_replica_closure(name, replica_s1[name], moments[name]["S1"], replica_count)

    continuum_component_total = math.fsum(moments["continuum_components"]["S1"])
    continuum_total = math.fsum(moments["continuum_total"]["S1"])
    if not _close(continuum_component_total, continuum_total):
        raise BlockFormatError("continuum component S1 values do not close to total continuum")
    for line_cell in range(line_cells):
        component_sum = math.fsum(
            moments["line_components"]["S1"][component * line_cells + line_cell]
            for component in range(len(components))
        )
        if not _close(component_sum, moments["line_total"]["S1"][line_cell]):
            raise BlockFormatError("line component S1 values do not close cellwise")
    for replica in range(replica_count):
        component_start = replica * component_cells
        total_start = replica * total_cells
        if not _close(
            math.fsum(replica_s1["continuum_components"][component_start : component_start + component_cells]),
            math.fsum(replica_s1["continuum_total"][total_start : total_start + total_cells]),
            float32=True,
        ):
            raise BlockFormatError("replica continuum components do not close to total")
        line_component_start = replica * line_component_cells
        line_total_start = replica * line_cells
        for line_cell in range(line_cells):
            component_sum = math.fsum(
                replica_s1["line_components"][
                    line_component_start + component * line_cells + line_cell
                ]
                for component in range(len(components))
            )
            if not _close(
                component_sum,
                replica_s1["line_total"][line_total_start + line_cell],
                float32=True,
            ):
                raise BlockFormatError("replica line components do not close cellwise")
    if len(diagnostics["interaction_counts"]) != len(interaction_names):
        raise BlockFormatError("interaction names/counts have different lengths")
    if tuple(interaction_names) != _CANONICAL_INTERACTIONS:
        raise BlockFormatError("interaction diagnostic taxonomy is noncanonical")
    channel_names = [channel["name"] for channel in diagnostic_channels]
    if tuple(channel_names) != _CANONICAL_DIAGNOSTIC_CHANNELS:
        raise BlockFormatError("path diagnostic channel taxonomy is invalid")
    for diagnostic_name in ("scattered_total", "dropped"):
        _validate_numeric_moments(
            f"diagnostics.{diagnostic_name}",
            {key: (diagnostics[diagnostic_name][key],) for key in ("N", "S1", "S2")},
            node["histories"],
        )
    for channel in diagnostic_channels:
        _validate_moments(f"diagnostic:{channel['name']}:continuum", channel["continuum"], total_cells)
        _validate_moments(f"diagnostic:{channel['name']}:lines", channel["lines"], line_cells)
        _validate_numeric_moments(
            f"diagnostic:{channel['name']}:continuum", channel["continuum"], node["histories"]
        )
        _validate_numeric_moments(
            f"diagnostic:{channel['name']}:lines", channel["lines"], node["histories"]
        )
    for cell in range(total_cells):
        channel_sum = math.fsum(
            channel["continuum"]["S1"][cell] for channel in diagnostic_channels
        )
        if not _close(channel_sum, moments["continuum_total"]["S1"][cell]):
            raise BlockFormatError("path diagnostic continuum channels do not close cellwise")
    for cell in range(line_cells):
        channel_sum = math.fsum(
            channel["lines"]["S1"][cell] for channel in diagnostic_channels
        )
        if not _close(channel_sum, moments["line_total"]["S1"][cell]):
            raise BlockFormatError("path diagnostic line channels do not close cellwise")
    if not _close(
        diagnostics["scattered_total"]["S1"],
        continuum_total + math.fsum(moments["line_total"]["S1"]),
    ):
        raise BlockFormatError("scattered scalar does not close to continuum plus exact lines")
    if stopping["N"] != node["histories"]:
        raise BlockFormatError("stopping-vector N does not match node histories")
    if len(stopping["S1"]) != 16 or len(stopping["outer_S2"]) != 256:
        raise BlockFormatError("stopping vector/matrix dimensions are invalid")
    if any(not math.isfinite(value) for value in (*stopping["S1"], *stopping["outer_S2"])):
        raise BlockFormatError("stopping moments contain non-finite values")
    for row in range(16):
        diagonal = stopping["outer_S2"][row * 16 + row]
        if diagonal < -1.0e-14 or diagonal + 2.0e-12 * max(1.0, abs(diagonal)) < stopping["S1"][row] ** 2 / node["histories"]:
            raise BlockFormatError("stopping outer moment has an invalid diagonal")
        for column in range(row):
            left = stopping["outer_S2"][row * 16 + column]
            right = stopping["outer_S2"][column * 16 + row]
            if not _close(left, right):
                raise BlockFormatError("stopping outer moment is not symmetric")
    band_total_s1 = math.fsum(stopping["S1"][:12])
    band_total_s2 = math.fsum(
        stopping["outer_S2"][row * 16 + column]
        for row in range(12)
        for column in range(12)
    )
    if not _close(band_total_s1, diagnostics["scattered_total"]["S1"]):
        raise BlockFormatError("registered stopping-band S1 does not close to scattered total")
    if not _close(band_total_s2, diagnostics["scattered_total"]["S2"]):
        raise BlockFormatError("registered stopping-band S2 does not close to scattered total")

    if node["role"] not in {"pilot", "train", "validation", "final", "diagnostic"}:
        raise BlockFormatError("node role is not recognized")
    if node["status"] not in {"complete", "converged", "capped", "failed"}:
        raise BlockFormatError("node status is not recognized")
    for name in (
        "direct_transmission",
        "direct_fluence",
        "wall_seconds",
        "cpu_seconds",
        "target_total_rel_se",
        "achieved_total_rel_se",
    ):
        if not math.isfinite(node[name]) or node[name] < 0.0:
            raise BlockFormatError(f"node scalar {name} is non-finite or negative")
    if (
        not math.isfinite(node["energy_keV"])
        or node["energy_keV"] <= configured_e_cut
        or not math.isfinite(node["distance_cm"])
        or node["distance_cm"] <= 0.0
        or not math.isfinite(node["source_height_cm"])
        or node["source_height_cm"] < 0.0
        or not math.isfinite(node["detector_height_cm"])
        or node["detector_height_cm"] < 0.0
        or not math.isfinite(node["extent_cm"])
        or node["extent_cm"] < 0.0
        or node["distance_cm"]
        <= abs(node["source_height_cm"] - node["detector_height_cm"])
    ):
        raise BlockFormatError("node coordinates are outside the physical fit domain")
    if configuration_tier == "air" and any(
        node[name] != 0.0
        for name in ("source_height_cm", "detector_height_cm", "extent_cm")
    ):
        raise BlockFormatError("air-tier node has nonzero field-geometry coordinates")
    expected_direct_fluence = node["direct_transmission"] / (
        4.0 * math.pi * node["distance_cm"] * node["distance_cm"]
    )
    if (
        node["direct_transmission"] <= 0.0
        or node["direct_transmission"] > 1.0 + 2.0e-12
        or node["direct_fluence"] <= 0.0
        or not math.isclose(
            node["direct_fluence"],
            expected_direct_fluence,
            rel_tol=2.0e-12,
            abs_tol=1.0e-300,
        )
        or node["target_total_rel_se"] <= 0.0
    ):
        raise BlockFormatError("direct-beam normalization or stopping target is noncanonical")
    expected_status = (
        "failed"
        if diagnostics["transport_failures"]
        else "capped"
        if diagnostics["capped"]
        else node["status"]
    )
    if (
        node["status"] != expected_status
        or (node["status"] in {"complete", "converged"} and diagnostics["failure_reason"])
        or (node["status"] in {"capped", "failed"} and not diagnostics["failure_reason"])
    ):
        raise BlockFormatError("node status is inconsistent with transport diagnostics")
    if (
        not math.isfinite(diagnostics["max_history_contribution"])
        or diagnostics["max_history_contribution"] < 0.0
    ):
        raise BlockFormatError("maximum history contribution is invalid")
    expected_node_id = stable_node_id_from_cm(
        manifest["config_sha256"],
        node["energy_keV"],
        node["distance_cm"],
        node["source_height_cm"],
        node["detector_height_cm"],
        extent_cm=node["extent_cm"],
    )
    if node_id != expected_node_id:
        raise BlockFormatError("node ID does not match configuration and coordinates")
    expected_q = geometry_transform(
        PhysicalNode(
            node["energy_keV"],
            node["distance_cm"] * 0.01,
            node["source_height_cm"] * 0.01,
            node["detector_height_cm"] * 0.01,
        )
    )
    # Near rho_min, rho is obtained by subtracting nearly equal squared
    # distances.  Correct C++ and Python binary64 evaluations can consequently
    # differ at roughly 1e-12 in log(rho).  Keep the comparison strict enough
    # to catch altered metadata without rejecting that conditioned roundoff.
    if any(
        not math.isclose(actual, expected, rel_tol=2.0e-12, abs_tol=2.0e-12)
        for actual, expected in zip(node["q"], expected_q)
    ):
        raise BlockFormatError("stored q transform does not match physical coordinates")
    if block_id == 0:
        raise BlockFormatError("block ID zero is reserved for in-memory aggregates")

    return FitBlock(
        node_id,
        block_id,
        manifest,
        axes,
        tuple(components),
        node,
        replicas,
        moments,
        replica_s1,
        stopping,
        diagnostics,
        payload_sha256,
        source_path,
    )


def read_fit_block(source: str | Path | BinaryIO) -> FitBlock:
    """Read, checksum, fully consume, and structurally validate one block."""

    source_path: str | None = None
    close = False
    if isinstance(source, (str, Path)):
        source_path = str(source)
        stream: BinaryIO = open(source, "rb")
        close = True
    else:
        stream = source
    try:
        header = stream.read(HEADER.size)
        if len(header) != HEADER.size:
            raise BlockFormatError("file is shorter than the 64-byte CEEFIT header")
        magic, version, flags, payload_bytes, expected_sha = HEADER.unpack(header)
        if magic != MAGIC:
            raise BlockFormatError("CEEFIT magic mismatch")
        if version != VERSION:
            raise BlockFormatError(f"unsupported CEEFIT version {version}")
        if flags != 0:
            raise BlockFormatError(f"unsupported CEEFIT header flags 0x{flags:x}")
        if payload_bytes > MAX_PAYLOAD_BYTES:
            raise BlockFormatError("declared payload exceeds safety limit")
        payload = stream.read(payload_bytes)
        if len(payload) != payload_bytes:
            raise BlockFormatError("truncated CEEFIT payload")
        if stream.read(1):
            raise BlockFormatError("bytes follow the declared CEEFIT payload")
        actual_sha = hashlib.sha256(payload).digest()
        if actual_sha != expected_sha:
            raise BlockFormatError("CEEFIT payload SHA-256 mismatch")
        return _decode_payload(payload, actual_sha.hex(), source_path)
    finally:
        if close:
            stream.close()
