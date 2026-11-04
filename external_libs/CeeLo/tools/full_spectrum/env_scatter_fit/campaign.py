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

"""Tracked campaign manifests and one-node/one-block worker commands."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import struct
from typing import Any, Iterable, Mapping, Sequence

from .allocation import BLOCK_HISTORIES, MIN_INITIAL_BLOCKS
from .design import (
    DesignBounds,
    DesignSelection,
    PhysicalNode,
    generate_candidate_pool,
    geometry_transform,
    select_disjoint_maximin,
    transformed_domain,
)
from .noise import (
    EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
    EXACT_LINE_RELATIVE_NOISE_FRACTION,
    contribution_aware_model_selection_policy,
    direct_line_map,
    exact_line_materiality,
    replica_shape_tv_noise,
    sentinel_shape,
)


CAMPAIGN_SCHEMA = "CEEFIT-CAMPAIGN-v1"
MODEL_SELECTION_REPORT_SCHEMA = "CEEFIT-MODEL-SELECTION-VALIDATION-v1"
MODEL_SELECTION_ROUND_IDENTITY_SCHEMA = "CEEFIT-CAMPAIGN-ROUND-INPUT-v1"
MODEL_SELECTION_HISTORY_SCHEMA = "CEEFIT-VALIDATION-HISTORY-v1"
DIAGNOSTIC_FAILURE_SCHEMA = "CEEFIT-DIAGNOSTIC-VALIDATION-FAILURE-v1"
ACQUISITION_BATCH_SIZE = 64
MAX_TRAINING_NODES = 2_048
MAXIMUM_FINAL_GENERATIONS = 3
# The bake build directory convention: every worker binary is built into
# <CeeLo root>/build-rel/tools/full_spectrum/, and campaign commands run from
# the CeeLo root.
DEFAULT_WORKER_EXECUTABLE = "build-rel/tools/full_spectrum/env_kernel_bake"


_MATERIALS: dict[str, tuple[str, float, tuple[tuple[int, float], ...]]] = {
    "air": ("Air", 0.001205, ((7, 0.7553), (8, 0.2318), (18, 0.0128), (6, 0.0001))),
    "soil": (
        "Soil",
        1.5,
        (
            (8, 0.460),
            (14, 0.276),
            (13, 0.083),
            (26, 0.050),
            (20, 0.037),
            (12, 0.021),
            (11, 0.021),
            (19, 0.026),
            (22, 0.006),
            (1, 0.020),
        ),
    ),
    "concrete": (
        "Concrete",
        2.3,
        (
            (1, 0.0100),
            (6, 0.0010),
            (8, 0.5299),
            (11, 0.0160),
            (12, 0.0020),
            (13, 0.0339),
            (14, 0.3374),
            (19, 0.0130),
            (20, 0.0440),
            (26, 0.0140),
        ),
    ),
}


@dataclass(frozen=True)
class PhysicsConfig:
    energy_bins: int = 384
    energy_cut_keV: float = 10.0
    mu_bins: int = 64
    brems_electron_floor_keV: float = 150.0
    doppler: bool = False
    free_electron_compton: bool = False
    include_annihilation: bool = True
    include_brems: bool = True
    include_rayleigh: bool = True
    max_vertices: int = 50
    positron_range: bool = True
    roulette_tau: float = 8.0
    walk_annihilation: bool = True
    walk_brems: bool = True
    weight_floor: float = 1.0e-5
    scoring_radius_cm: float = 4.0

    def validate(self) -> None:
        if self.energy_bins != 384:
            raise ValueError("CEEFIT-v1 requires 384 registered energy bins")
        if self.energy_cut_keV < 10.0 or self.mu_bins not in {32, 64}:
            raise ValueError("invalid fit axis: Ecut >= 10 and mu bins 32 or 64 are required")
        if self.max_vertices <= 0 or self.brems_electron_floor_keV < 0.0:
            raise ValueError("invalid transport termination configuration")
        if self.scoring_radius_cm <= 0.0:
            raise ValueError("scoring radius must be positive")


def _number(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("canonical configuration cannot contain non-finite numbers")
    return format(value, ".17g")


def _json_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


def material_json(key: str) -> str:
    try:
        name, density, composition = _MATERIALS[key]
    except KeyError as error:
        raise ValueError(f"unsupported material {key!r}") from error
    parts = ",".join(
        '{"Z":' + str(z) + ',"mass_fraction":' + _number(fraction) + "}"
        for z, fraction in composition
    )
    return (
        '{"composition":['
        + parts
        + '],"density_g_cm3":'
        + _number(density)
        + ',"name":'
        + _json_string(name)
        + "}"
    )


def canonical_fit_config_json(
    mode: str,
    material: str,
    physics: PhysicsConfig = PhysicsConfig(),
) -> str:
    """Reproduce env_kernel_bake's byte-exact canonical configuration."""

    physics.validate()
    if mode not in {"air", "field"}:
        raise ValueError("initial fit campaign supports only air or field")
    if material not in {"soil", "concrete"}:
        raise ValueError("field material must be soil or concrete")
    boolean = lambda value: "true" if value else "false"
    return (
        '{"air":'
        + material_json("air")
        + ',"axis":{"energy_bins":'
        + str(physics.energy_bins)
        + ',"energy_cut_keV":'
        + _number(physics.energy_cut_keV)
        + ',"mu_bins":'
        + str(physics.mu_bins)
        + ',"transform":"registered-component-v1"}'
        + ',"block_histories":262144'
        + ',"geometry_tier":'
        + _json_string(mode)
        + ',"material":'
        + ("null" if mode == "air" else material_json(material))
        + ',"physics":{"brems_electron_floor_keV":'
        + _number(physics.brems_electron_floor_keV)
        + ',"doppler":'
        + boolean(physics.doppler)
        + ',"free_electron_compton":'
        + boolean(physics.free_electron_compton)
        + ',"include_annihilation":'
        + boolean(physics.include_annihilation)
        + ',"include_brems":'
        + boolean(physics.include_brems)
        + ',"include_rayleigh":'
        + boolean(physics.include_rayleigh)
        + ',"max_vertices":'
        + str(physics.max_vertices)
        + ',"positron_range":'
        + boolean(physics.positron_range)
        + ',"roulette_tau":'
        + _number(physics.roulette_tau)
        + ',"walk_annihilation":'
        + boolean(physics.walk_annihilation)
        + ',"walk_brems":'
        + boolean(physics.walk_brems)
        + ',"weight_floor":'
        + _number(physics.weight_floor)
        + '},"scoring_radius_cm":'
        + _number(physics.scoring_radius_cm)
        + ',"schema":"ceelo-env-fit-config-v1"}'
    )


def config_sha256(canonical_config_json: str) -> str:
    return hashlib.sha256(canonical_config_json.encode("utf-8")).hexdigest()


def replica_seed(config_hash: str, node_id: int, block_id: int, replica_id: int) -> int:
    if block_id <= 0 or not 0 <= replica_id < 32:
        raise ValueError("invalid block or replica id")
    raw = (
        b"CEEFIT-RNG-v1\0"
        + config_hash.encode("ascii")
        + struct.pack("<QQI", node_id, block_id, replica_id)
    )
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def _inside(node: PhysicalNode, bounds: DesignBounds) -> bool:
    dz = node.source_height_m - node.detector_height_m
    rho2 = node.distance_m * node.distance_m - dz * dz
    return (
        bounds.energy_keV[0] <= node.energy_keV <= bounds.energy_keV[1]
        and bounds.distance_m[0] <= node.distance_m <= bounds.distance_m[1]
        and bounds.source_height_m[0] <= node.source_height_m <= bounds.source_height_m[1]
        and bounds.detector_height_m[0] <= node.detector_height_m <= bounds.detector_height_m[1]
        and abs(node.source_height_m - node.detector_height_m) < node.distance_m
        and rho2 >= bounds.rho_min_m * bounds.rho_min_m
    )


def default_forced_anchors(
    bounds: DesignBounds,
    *,
    mode: str = "field",
) -> tuple[PhysicalNode, ...]:
    """G4 anchors, thresholds, and valid domain-face sentinels."""

    e_mid = math.sqrt(bounds.energy_keV[0] * bounds.energy_keV[1])
    d_mid = math.sqrt(bounds.distance_m[0] * bounds.distance_m[1])
    hs_mid = sum(bounds.source_height_m) / 2.0
    hd_mid = sum(bounds.detector_height_m) / 2.0
    if mode == "air":
        raw = [
            PhysicalNode(energy, distance, 0.0, 0.0)
            for energy in (
                bounds.energy_keV[0],
                122.0,
                149.0,
                151.0,
                662.0,
                2.0 * 510.998950 - 0.1,
                2.0 * 510.998950 + 0.1,
                bounds.energy_keV[1],
            )
            for distance in (bounds.distance_m[0], d_mid, bounds.distance_m[1])
        ]
        # Air deliberately collapses the two physically irrelevant height axes.
        return tuple(dict.fromkeys(raw))
    raw = [
        PhysicalNode(662.0, 5.0, 1.0, 1.0),
        PhysicalNode(662.0, 5.0, 0.3, 1.0),
        PhysicalNode(122.0, 5.0, 1.0, 1.0),
        PhysicalNode(2614.0, 5.0, 1.0, 1.0),
    ]
    raw.extend(
        PhysicalNode(energy, d_mid, hs_mid, hd_mid)
        for energy in (
            bounds.energy_keV[0],
            149.0,
            151.0,
            2.0 * 510.998950 - 0.1,
            2.0 * 510.998950 + 0.1,
            bounds.energy_keV[1],
        )
    )
    raw.extend(
        (
            PhysicalNode(e_mid, bounds.distance_m[0], bounds.source_height_m[0], bounds.detector_height_m[0]),
            PhysicalNode(e_mid, bounds.distance_m[0], bounds.source_height_m[1], bounds.detector_height_m[1]),
            PhysicalNode(e_mid, bounds.distance_m[1], bounds.source_height_m[0], bounds.detector_height_m[1]),
            PhysicalNode(e_mid, bounds.distance_m[1], bounds.source_height_m[1], bounds.detector_height_m[0]),
            PhysicalNode(e_mid, d_mid, bounds.source_height_m[0], hd_mid),
            PhysicalNode(e_mid, d_mid, bounds.source_height_m[1], hd_mid),
            PhysicalNode(e_mid, d_mid, hs_mid, bounds.detector_height_m[0]),
            PhysicalNode(e_mid, d_mid, hs_mid, bounds.detector_height_m[1]),
        )
    )
    return tuple(dict.fromkeys(node for node in raw if _inside(node, bounds)))


def _worker_command(
    executable: str,
    output_root: str,
    node: Mapping[str, Any],
    block_id: int,
    *,
    mode: str,
    material: str,
    physics: PhysicsConfig,
) -> list[str]:
    output = (
        Path(output_root)
        / str(node["role"])
        / str(node["node_id_hex"])
        / f"block_{block_id}.ceefit"
    )
    command = [
        executable,
        "--mode",
        mode,
        "--energies",
        _number(float(node["energy_keV"])),
        "--distances",
        _number(float(node["distance_m"]) * 100.0),
        "--histories",
        str(BLOCK_HISTORIES),
        "--registered",
        str(physics.energy_bins),
        "--e-floor",
        _number(physics.energy_cut_keV),
        "--mu-bins",
        str(physics.mu_bins),
        "--brems-floor",
        _number(physics.brems_electron_floor_keV),
        "--max-vertices",
        str(physics.max_vertices),
        "--scoring-radius",
        _number(physics.scoring_radius_cm),
        "--fit-block",
        str(output),
        "--block-id",
        str(block_id),
        "--role",
        str(node["role"]),
    ]
    if mode == "field":
        command += [
            "--material",
            material,
            "--h-src",
            _number(float(node["source_height_m"]) * 100.0),
            "--h-det",
            _number(float(node["detector_height_m"]) * 100.0),
        ]
    else:
        command += ["--extent", "0"]
    # Explicitly represent any non-default physics setting supported by the C++
    # executable. The initial campaign intentionally keeps all production
    # physics enabled.
    if not physics.include_rayleigh:
        command.append("--no-rayleigh")
    if not physics.positron_range:
        command.append("--no-positron-range")
    if not physics.walk_annihilation:
        command.append("--no-walk-annihilation")
    if not physics.include_brems:
        command.append("--no-brems")
    if not physics.walk_brems:
        command.append("--no-walk-brems")
    if physics.free_electron_compton:
        command.append("--free-electron")
    return command


def _build_work_for_nodes(
    nodes: Sequence[Mapping[str, Any]],
    *,
    first_block_id: int,
    block_count: int,
    executable: str,
    block_output_root: str,
    mode: str,
    material: str,
    physics: PhysicsConfig,
    config_hash: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Materialize exact logical work and grouped persistent-worker commands."""

    if first_block_id <= 0 or block_count <= 0:
        raise ValueError("work generation requires positive block IDs/counts")
    histories_per_replica = BLOCK_HISTORIES // 32
    work: list[dict[str, Any]] = []
    worker_groups: list[dict[str, Any]] = []
    for node in nodes:
        node_work: list[dict[str, Any]] = []
        for block_id in range(first_block_id, first_block_id + block_count):
            replicas = [
                {
                    "replica_id": replica,
                    "history_begin": (block_id - 1) * BLOCK_HISTORIES
                    + replica * histories_per_replica,
                    "history_count": histories_per_replica,
                    "rng_seed": replica_seed(
                        config_hash, int(node["node_id"]), block_id, replica
                    ),
                }
                for replica in range(32)
            ]
            command = _worker_command(
                executable,
                block_output_root,
                node,
                block_id,
                mode=mode,
                material=material,
                physics=physics,
            )
            output = command[command.index("--fit-block") + 1]
            item = {
                "node_id": node["node_id"],
                "node_id_hex": node["node_id_hex"],
                "role": node["role"],
                "block_id": block_id,
                "histories": BLOCK_HISTORIES,
                "replicas": replicas,
                "output": output,
                "argv": command,
                "shell_command": shlex.join(
                    ["mkdir", "-p", str(Path(output).parent)]
                )
                + " && "
                + shlex.join(command),
            }
            node_work.append(item)
            work.append(item)
        first = node_work[0]
        grouped = list(first["argv"])
        output_index = grouped.index("--fit-block") + 1
        output_template = str(Path(first["output"]).with_name("block_{block_id}.ceefit"))
        grouped[output_index] = output_template
        grouped += ["--block-count", str(block_count)]
        worker_groups.append(
            {
                "node_id": node["node_id"],
                "node_id_hex": node["node_id_hex"],
                "role": node["role"],
                "first_block_id": first_block_id,
                "block_count": block_count,
                "argv": grouped,
                "shell_command": shlex.join(
                    ["mkdir", "-p", str(Path(output_template).parent)]
                )
                + " && "
                + shlex.join(grouped),
            }
        )
    return work, worker_groups


def _coordinate_key(node: PhysicalNode | Mapping[str, Any]) -> tuple[str, str, str, str]:
    if isinstance(node, PhysicalNode):
        values = (
            node.energy_keV,
            node.distance_m,
            node.source_height_m,
            node.detector_height_m,
        )
    else:
        values = tuple(
            float(node[name])
            for name in (
                "energy_keV",
                "distance_m",
                "source_height_m",
                "detector_height_m",
            )
        )
    return tuple(value.hex() for value in values)  # type: ignore[return-value]


def _encode_node(node: PhysicalNode, role: str, config_hash: str) -> dict[str, Any]:
    return {
        **node.as_dict(),
        "q": list(geometry_transform(node)),
        "role": role,
        "node_id": node.stable_id(config_hash),
        "node_id_hex": node.stable_id_hex(config_hash),
    }


def _physical_node(record: Mapping[str, Any]) -> PhysicalNode:
    return PhysicalNode(
        float(record["energy_keV"]),
        float(record["distance_m"]),
        float(record["source_height_m"]),
        float(record["detector_height_m"]),
    )


def _validate_manifest_record(manifest: Mapping[str, Any]) -> None:
    if manifest.get("schema") != CAMPAIGN_SCHEMA:
        raise ValueError(f"unsupported campaign schema {manifest.get('schema')!r}")
    canonical = manifest["config"]["canonical_json"]
    if config_sha256(canonical) != manifest["config"]["sha256"]:
        raise ValueError("campaign canonical configuration checksum mismatch")
    allocation = manifest.get("allocation")
    state = manifest.get("campaign_state")
    if (
        not isinstance(allocation, Mapping)
        or allocation.get("maximum_final_generations")
        != MAXIMUM_FINAL_GENERATIONS
        or not isinstance(state, Mapping)
        or type(state.get("final_generation")) is not int
        or not 0 <= state["final_generation"] < MAXIMUM_FINAL_GENERATIONS
    ):
        raise ValueError("campaign final-generation policy is invalid")
    physics = PhysicsConfig(**manifest["config"]["physics"])
    decision = manifest.get("mu_resolution_decision")
    if physics.mu_bins == 32:
        _validated_mu_resolution_decision(
            decision,
            mode=str(manifest["config"]["mode"]),
            material=str(manifest["config"]["material"]),
            production_physics=physics,
        )
    elif decision is not None:
        raise ValueError("a 64-bin campaign must not carry a 32-bin production decision")


def _validated_mu_resolution_decision(
    decision: Any,
    *,
    mode: str,
    material: str,
    production_physics: PhysicsConfig,
    verify_sources: bool = False,
) -> dict[str, Any]:
    """Validate the identity-bound 64→32 pilot decision embedded in a campaign."""

    if not isinstance(decision, Mapping):
        raise ValueError("a 32-bin campaign requires a passing 64-to-32 decision")
    report = decision.get("report")
    digest = decision.get("canonical_sha256")
    if not isinstance(report, Mapping) or not isinstance(digest, str):
        raise ValueError("mu-resolution decision record is malformed")
    try:
        canonical_report = json.dumps(
            report, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ValueError("mu-resolution decision is not finite canonical JSON") from error
    if hashlib.sha256(canonical_report).hexdigest() != digest:
        raise ValueError("mu-resolution decision canonical checksum differs")
    if (
        report.get("schema") != "CEEFIT-MU-DECISION-v1"
        or report.get("source_bins") != 64
        or report.get("candidate_bins") != 32
        or report.get("passed") is not True
        or report.get("decision") != "use_32"
    ):
        raise ValueError("mu-resolution decision does not authorize 32 production bins")
    tolerance = float(report.get("tolerance", math.inf))
    max_tv = float(report.get("max_joint_tv", math.inf))
    if not (
        math.isfinite(tolerance)
        and 0.0 < tolerance <= 0.002
        and math.isfinite(max_tv)
        and 0.0 <= max_tv <= tolerance
    ):
        raise ValueError("mu-resolution decision exceeds the 0.2% shape tolerance")
    archive_hash = report.get("archive_snapshot_sha256")
    if not (
        isinstance(archive_hash, str)
        and len(archive_hash) == 64
        and all(character in "0123456789abcdef" for character in archive_hash)
    ):
        raise ValueError("mu-resolution decision lacks an archive snapshot SHA-256")
    source_archive = report.get("source_archive")
    if not (
        isinstance(source_archive, Mapping)
        and isinstance(source_archive.get("path"), str)
        and source_archive["path"]
        and source_archive.get("archive_snapshot_sha256") == archive_hash
    ):
        raise ValueError("mu-resolution decision lacks its reopenable source archive")
    source_detector = report.get("source_detector_validation")
    if not (
        isinstance(source_detector, Mapping)
        and isinstance(source_detector.get("path"), str)
        and source_detector["path"]
        and _is_sha256(source_detector.get("report_sha256"))
    ):
        raise ValueError(
            "mu-resolution decision lacks its reopenable detector-fold report"
        )
    pilot_ids = [str(value).lower().removeprefix("0x") for value in report.get("pilot_node_ids", ())]
    if not pilot_ids or len(pilot_ids) != len(set(pilot_ids)) or any(
        len(value) != 16 or any(character not in "0123456789abcdef" for character in value)
        for value in pilot_ids
    ):
        raise ValueError("mu-resolution decision has invalid pilot node identities")
    detector = report.get("detector_validation")
    if not isinstance(detector, Mapping) or detector.get("passed") is not True:
        raise ValueError("mu-resolution decision lacks passing bound HPGe/NaI folds")
    detector_records = detector.get("detectors")
    if not isinstance(detector_records, Mapping):
        raise ValueError("mu-resolution decision lacks detector-fold details")
    metric_names = (
        "max_detector_folded_total_relative_change",
        "max_detector_folded_spectral_tv",
        "max_exact_line_angular_fold_relative_change",
    )
    for detector_name in ("HPGe", "NaI"):
        detector_record = detector_records.get(detector_name)
        if not isinstance(detector_record, Mapping):
            raise ValueError(f"mu-resolution decision lacks {detector_name} folds")
        if not _is_sha256(detector_record.get("response_sha256")):
            raise ValueError(
                f"mu-resolution decision has an invalid {detector_name} response hash"
            )
        if (
            detector_record.get("passed") is not True
            or detector_record.get("node_count") != len(pilot_ids)
        ):
            raise ValueError(
                f"mu-resolution decision does not cover every {detector_name} pilot fold"
            )
        for metric_name in metric_names:
            value = float(detector_record.get(metric_name, math.inf))
            if not math.isfinite(value) or not 0.0 <= value <= tolerance:
                raise ValueError(
                    f"mu-resolution decision {detector_name} exceeds the 0.2% fold tolerance"
                )

    pilot_physics_record = asdict(production_physics)
    pilot_physics_record["mu_bins"] = 64
    expected_config = config_sha256(
        canonical_fit_config_json(mode, material, PhysicsConfig(**pilot_physics_record))
    )
    if report.get("config_sha256") != expected_config:
        raise ValueError("mu-resolution decision is bound to another 64-bin configuration")
    if verify_sources:
        _verify_mu_resolution_sources(report)
    return {"report": deepcopy(dict(report)), "canonical_sha256": digest}


def _verify_mu_resolution_sources(report: Mapping[str, Any]) -> None:
    """Reopen the pilot archive and detector report named by a 64→32 decision."""

    from .archive import require_hdf5
    from .publication import archive_snapshot_sha256, file_sha256

    archive_record = report["source_archive"]
    detector_record = report["source_detector_validation"]
    archive_path = Path(str(archive_record["path"]))
    detector_path = Path(str(detector_record["path"]))
    if not archive_path.is_file():
        raise ValueError("mu-resolution source archive is missing")
    if not detector_path.is_file():
        raise ValueError("mu-resolution detector-fold source report is missing")
    h5py, _ = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        if archive_snapshot_sha256(root) != archive_record["archive_snapshot_sha256"]:
            raise ValueError("mu-resolution source archive snapshot is stale")
    if file_sha256(detector_path) != detector_record["report_sha256"]:
        raise ValueError("mu-resolution detector-fold source report is stale")


def _execution_arguments(manifest: Mapping[str, Any]) -> dict[str, Any]:
    execution = manifest["execution"]
    return {
        "executable": str(execution["executable"]),
        "block_output_root": str(execution["block_output_root"]),
        "mode": str(manifest["config"]["mode"]),
        "material": str(manifest["config"]["material"]),
        "physics": PhysicsConfig(**manifest["config"]["physics"]),
        "config_hash": str(manifest["config"]["sha256"]),
    }


def generate_campaign_manifest(
    *,
    mode: str = "field",
    material: str = "soil",
    bounds: DesignBounds = DesignBounds(),
    physics: PhysicsConfig = PhysicsConfig(),
    candidate_count: int = 10_000,
    train_count: int = 256,
    validation_count: int = 64,
    final_count: int = 64,
    seed: int = 0xCEE10,
    anchors: Iterable[PhysicalNode] | None = None,
    initial_blocks: int = MIN_INITIAL_BLOCKS,
    executable: str = DEFAULT_WORKER_EXECUTABLE,
    block_output_root: str = "ceefit_blocks",
    mu_resolution_decision: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Create a decision-complete deterministic pilot campaign manifest."""

    if initial_blocks < MIN_INITIAL_BLOCKS:
        raise ValueError(f"initial campaign needs at least {MIN_INITIAL_BLOCKS} blocks per node")
    if train_count > MAX_TRAINING_NODES:
        raise ValueError(f"initial training design exceeds the {MAX_TRAINING_NODES}-node cap")
    canonical = canonical_fit_config_json(mode, material, physics)
    config_hash = config_sha256(canonical)
    embedded_mu_decision = None
    if physics.mu_bins == 32:
        if mu_resolution_decision is None:
            raise ValueError(
                "32-bin production requires a passing identity-bound 64-to-32 decision"
            )
        report = deepcopy(dict(mu_resolution_decision))
        encoded_report = json.dumps(
            report, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        embedded_mu_decision = _validated_mu_resolution_decision(
            {
                "report": report,
                "canonical_sha256": hashlib.sha256(encoded_report).hexdigest(),
            },
            mode=mode,
            material=material,
            production_physics=physics,
            verify_sources=True,
        )
    elif mu_resolution_decision is not None:
        raise ValueError("mu-resolution decision applies only to 32-bin production")
    candidates = generate_candidate_pool(
        candidate_count, bounds=bounds, seed=seed, mode=mode
    )
    default_anchors = default_forced_anchors(bounds, mode=mode)
    additional_anchors = tuple(anchors or ())
    if mode == "air":
        additional_anchors = tuple(
            PhysicalNode(node.energy_keV, node.distance_m, 0.0, 0.0)
            for node in additional_anchors
        )
    if any(
        not (
            bounds.energy_keV[0] <= node.energy_keV <= bounds.energy_keV[1]
            and bounds.distance_m[0] <= node.distance_m <= bounds.distance_m[1]
        )
        if mode == "air"
        else not _inside(node, bounds)
        for node in additional_anchors
    ):
        raise ValueError("a forced campaign anchor lies outside the valid design bounds")
    forced_train = tuple(
        dict.fromkeys((*default_anchors, *additional_anchors))
    )
    q_min_tuple, q_max_tuple = transformed_domain(bounds, mode=mode)
    q_min = list(q_min_tuple)
    q_max = list(q_max_tuple)
    selection: DesignSelection = select_disjoint_maximin(
        candidates,
        sizes={"train": train_count, "validation": validation_count, "final": final_count},
        forced={"train": forced_train},
        q_domain=(q_min, q_max),
    )
    design = selection.as_dict(config_hash)
    selection_provenance = {
        "method": "greedy-maximin-q-v1",
        "forced_train": [
            _encode_node(node, "train", config_hash) for node in forced_train
        ],
        "initial_role_node_ids": {
            role: [str(item["node_id_hex"]) for item in design[role]]
            for role in ("train", "validation", "final")
        },
    }
    selection_provenance["canonical_sha256"] = _canonical_json_sha256(
        selection_provenance
    )
    nodes = [*design["train"], *design["validation"], *design["final"]]
    work, worker_groups = _build_work_for_nodes(
        nodes,
        first_block_id=1,
        block_count=initial_blocks,
        executable=executable,
        block_output_root=block_output_root,
        mode=mode,
        material=material,
        physics=physics,
        config_hash=config_hash,
    )
    selected_roles = {
        _coordinate_key(item): role
        for role in ("train", "validation", "final")
        for item in design[role]
    }
    candidate_pool = []
    unused_candidate_indices = []
    for index, candidate in enumerate(candidates):
        selected_role = selected_roles.get(_coordinate_key(candidate))
        candidate_pool.append(
            {
                "candidate_index": index,
                **_encode_node(candidate, "candidate", config_hash),
                "selected_role": selected_role,
            }
        )
        if selected_role is None:
            unused_candidate_indices.append(index)
    return {
        "schema": CAMPAIGN_SCHEMA,
        "mu_resolution_decision": embedded_mu_decision,
        "config": {
            "canonical_json": canonical,
            "sha256": config_hash,
            "mode": mode,
            "material": material,
            "physics": asdict(physics),
        },
        "execution": {
            "executable": executable,
            "block_output_root": block_output_root,
        },
        "design": {
            "bounds": asdict(bounds),
            "sobol_seed": seed,
            "candidate_count": candidate_count,
            "forced_train_count": len(forced_train),
            "selection": selection_provenance,
            "q_domain": {"min": q_min, "max": q_max},
            "nodes": design,
            "candidate_pool": candidate_pool,
        },
        "allocation": {
            "block_histories": BLOCK_HISTORIES,
            "logical_replicas": 32,
            "initial_blocks_per_node": initial_blocks,
            "total_relative_se_target": 0.0025,
            "material_contribution_fraction": 0.01,
            "contribution_absolute_se_fraction": 0.001,
            "sentinel_absolute_se_fraction": 0.0025,
            "training_history_cap": 64_000_000,
            "validation_history_cap": 128_000_000,
            "acquisition_batch_size": ACQUISITION_BATCH_SIZE,
            "maximum_training_nodes": MAX_TRAINING_NODES,
            "maximum_final_generations": MAXIMUM_FINAL_GENERATIONS,
        },
        "campaign_state": {
            "round": 0,
            "validation_history": [],
            "consecutive_validation_passes": 0,
            "stopped": False,
            "training_node_ids": [item["node_id_hex"] for item in design["train"]],
            "acquired_node_ids": [],
            "unused_candidate_indices": unused_candidate_indices,
            "retired_final_node_ids": [],
            "final_generation": 0,
            "acquisition_history": [],
            "last_action": "initial_design",
        },
        "work": work,
        "worker_groups": worker_groups,
    }


def _canonical_json_sha256(record: Mapping[str, Any]) -> str:
    try:
        encoded = json.dumps(
            record, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ValueError("validation evidence must be finite JSON") from error
    return hashlib.sha256(encoded).hexdigest()


def _is_sha256(value: Any) -> bool:
    return bool(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def campaign_round_identity(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Return the content-addressed input identity for the next fit round.

    Besides the requested round/final generation, bind the complete role
    overlay.  This prevents a report generated before a training acquisition,
    validation-set edit, or final-set refresh from being replayed later.
    """

    state = manifest["campaign_state"]
    validation_ids = [
        str(item["node_id_hex"]).lower()
        for item in manifest["design"]["nodes"]["validation"]
    ]
    final_ids = [
        str(item["node_id_hex"]).lower()
        for item in manifest["design"]["nodes"]["final"]
    ]
    training_ids = [str(value).lower() for value in state["training_node_ids"]]
    for name, values in (
        ("training", training_ids),
        ("validation", validation_ids),
        ("final", final_ids),
    ):
        if len(values) != len(set(values)) or any(
            len(value) != 16
            or any(character not in "0123456789abcdef" for character in value)
            for value in values
        ):
            raise ValueError(f"campaign has invalid or duplicate {name} node IDs")
    if set(training_ids) & set(validation_ids) or set(training_ids) & set(final_ids):
        raise ValueError("campaign role overlay contains overlapping node IDs")
    if set(validation_ids) & set(final_ids):
        raise ValueError("campaign validation and final node IDs overlap")
    identity = {
        "schema": MODEL_SELECTION_ROUND_IDENTITY_SCHEMA,
        "config_sha256": str(manifest["config"]["sha256"]),
        "campaign_round": int(state["round"]) + 1,
        "final_generation": int(state["final_generation"]),
        "training_node_ids": training_ids,
        "validation_node_ids": validation_ids,
        "final_node_ids": final_ids,
    }
    return {**identity, "identity_sha256": _canonical_json_sha256(identity)}


def validate_model_selection_report(
    manifest: Mapping[str, Any],
    report: Mapping[str, Any],
    *,
    expected_archive_snapshot_sha256: str | None = None,
    expected_model_content_sha256: str | None = None,
) -> dict[str, Any]:
    """Validate and reduce a model-selection report to an audit record.

    A passing report is accepted only when its archive and serialized model
    identities have also been independently resolved by the caller.  Failed
    reports may be recorded without those artifacts because they can only
    reset the consecutive-pass counter.
    """

    if not isinstance(report, Mapping):
        raise ValueError("model-selection validation report must be a JSON object")
    canonical_report_sha256 = _canonical_json_sha256(report)
    if report.get("schema") != MODEL_SELECTION_REPORT_SCHEMA:
        raise ValueError(
            f"unsupported model-selection report schema {report.get('schema')!r}"
        )
    if report.get("model_selection_policy") != (
        contribution_aware_model_selection_policy()
    ):
        raise ValueError(
            "model-selection report lacks the current contribution-aware policy"
        )
    identity = campaign_round_identity(manifest)
    scalar_bindings = {
        "config_sha256": identity["config_sha256"],
        "campaign_round": identity["campaign_round"],
        "final_generation": identity["final_generation"],
        "campaign_round_identity_sha256": identity["identity_sha256"],
    }
    for field, expected in scalar_bindings.items():
        if report.get(field) != expected:
            raise ValueError(f"model-selection report {field} binding differs")
    eligible_ids = report.get("eligible_validation_node_ids")
    if eligible_ids != identity["validation_node_ids"]:
        raise ValueError(
            "model-selection report does not bind the exact eligible validation nodes"
        )
    archive_hash = report.get("archive_snapshot_sha256")
    model_hash = report.get("model_content_sha256")
    if not _is_sha256(archive_hash):
        raise ValueError("model-selection report lacks a valid archive snapshot SHA-256")
    if not _is_sha256(model_hash):
        raise ValueError("model-selection report lacks a valid model-content SHA-256")

    selection = report.get("selection")
    if not isinstance(selection, Mapping) or not selection:
        raise ValueError("model-selection report lacks component validation results")
    selection_flags: list[bool] = []
    for name, item in selection.items():
        if not isinstance(name, str) or not isinstance(item, Mapping):
            raise ValueError("model-selection component records are malformed")
        flag = item.get("validation_passed")
        if type(flag) is not bool:
            raise ValueError(
                "each model-selection component must carry a boolean validation_passed"
            )
        selection_flags.append(flag)
    derived_passed = all(selection_flags)
    if type(report.get("passed")) is not bool or report["passed"] != derived_passed:
        raise ValueError("model-selection report pass flag differs from its gate records")
    adaptive_hash = (
        report.get("acquisition", {}).get("utilities_artifact_sha256")
        if isinstance(report.get("acquisition"), Mapping)
        else None
    )
    if adaptive_hash is not None and not _is_sha256(adaptive_hash):
        raise ValueError("model-selection report adaptive utility hash is invalid")

    if derived_passed:
        if expected_archive_snapshot_sha256 is None:
            raise ValueError(
                "a passing validation report requires an independently opened archive"
            )
        if expected_model_content_sha256 is None:
            raise ValueError(
                "a passing validation report requires an independently loaded candidate"
            )
    if (
        expected_archive_snapshot_sha256 is not None
        and archive_hash != expected_archive_snapshot_sha256
    ):
        raise ValueError("model-selection report is bound to another archive snapshot")
    if (
        expected_model_content_sha256 is not None
        and model_hash != expected_model_content_sha256
    ):
        raise ValueError("model-selection report is bound to another model content")

    return {
        "schema": MODEL_SELECTION_HISTORY_SCHEMA,
        "kind": "bound_model_selection_report",
        "passed": derived_passed,
        "report_schema": MODEL_SELECTION_REPORT_SCHEMA,
        "report_canonical_sha256": canonical_report_sha256,
        "config_sha256": identity["config_sha256"],
        "archive_snapshot_sha256": archive_hash,
        "model_content_sha256": model_hash,
        "eligible_validation_node_ids": list(eligible_ids),
        "campaign_round": identity["campaign_round"],
        "final_generation": identity["final_generation"],
        "campaign_round_identity_sha256": identity["identity_sha256"],
        "adaptive_utility_report_sha256": adaptive_hash,
    }


def advance_campaign_round(
    manifest: Mapping[str, Any],
    *,
    validation_report: Mapping[str, Any] | None = None,
    diagnostic_validation_failed: bool = False,
    expected_archive_snapshot_sha256: str | None = None,
    expected_model_content_sha256: str | None = None,
) -> dict[str, Any]:
    """Persist one adaptive validation/acquisition round.

    Failed rounds are recorded but do not choose work here. The archive-bound
    unified scheduler subsequently compares a 64-node acquisition batch with
    continuations (each at most doubling its node) in the same predicted
    validation-loss/CPU units.
    A first passing round schedules no acquisition, and the second consecutive
    pass locks the campaign.
    """

    _validate_manifest_record(manifest)
    updated = deepcopy(dict(manifest))
    state = updated["campaign_state"]
    if state["stopped"]:
        raise ValueError("campaign has already stopped after two passing rounds")
    if (validation_report is None) == (not diagnostic_validation_failed):
        raise ValueError(
            "provide exactly one of a bound validation report or the explicit "
            "diagnostic-validation-failed path"
        )
    identity = campaign_round_identity(updated)
    if validation_report is not None:
        history_record = validate_model_selection_report(
            updated,
            validation_report,
            expected_archive_snapshot_sha256=expected_archive_snapshot_sha256,
            expected_model_content_sha256=expected_model_content_sha256,
        )
    else:
        history_record = {
            "schema": MODEL_SELECTION_HISTORY_SCHEMA,
            "kind": "diagnostic_manual_failure",
            "diagnostic_schema": DIAGNOSTIC_FAILURE_SCHEMA,
            "passed": False,
            "config_sha256": identity["config_sha256"],
            "campaign_round": identity["campaign_round"],
            "final_generation": identity["final_generation"],
            "campaign_round_identity_sha256": identity["identity_sha256"],
        }
    validation_gates_passed = history_record["passed"] is True
    state["round"] = identity["campaign_round"]
    state["validation_history"].append(history_record)
    consecutive = (
        int(state["consecutive_validation_passes"]) + 1
        if validation_gates_passed
        else 0
    )
    should_stop = consecutive >= 2
    state["consecutive_validation_passes"] = consecutive
    state["stopped"] = should_stop
    if should_stop:
        state["last_action"] = "stop_after_two_consecutive_validation_passes"
        return updated
    if validation_gates_passed:
        state["last_action"] = "await_confirmation_validation_round"
        return updated
    state["last_action"] = "await_unified_adaptive_allocation"
    return updated


def _validate_final_refresh_report(
    manifest: Mapping[str, Any], report: Mapping[str, Any]
) -> str:
    """Require one bound, adequately resolved failed-final adjudication."""

    if not isinstance(report, Mapping) or report.get("schema") != MODEL_SELECTION_REPORT_SCHEMA:
        raise ValueError("final refresh requires a CEEFIT model-selection report")
    identity = campaign_round_identity(manifest)
    for name, expected in (
        ("config_sha256", identity["config_sha256"]),
        ("campaign_round", identity["campaign_round"]),
        ("final_generation", identity["final_generation"]),
        ("campaign_round_identity_sha256", identity["identity_sha256"]),
    ):
        if report.get(name) != expected:
            raise ValueError(f"final-refresh report {name} binding differs")
    if not _is_sha256(report.get("archive_snapshot_sha256")) or not _is_sha256(
        report.get("model_content_sha256")
    ):
        raise ValueError("final-refresh report lacks archive/model identities")
    if manifest["campaign_state"].get("stopped") is not True:
        raise ValueError("final set cannot be refreshed before the two-pass stop")

    current_ids = [
        str(item["node_id_hex"]).lower()
        for item in manifest["design"]["nodes"]["final"]
    ]
    final = report.get("final")
    policy = report.get("policy")
    if not isinstance(final, Mapping) or not isinstance(policy, Mapping):
        raise ValueError("final-refresh report lacks final adjudication and policy")
    reported_ids = [
        str(item.get("node_id", "")).lower().removeprefix("0x")
        for item in final.get("nodes", ())
        if isinstance(item, Mapping)
    ]
    if len(reported_ids) != len(current_ids) or set(reported_ids) != set(current_ids):
        raise ValueError("final-refresh report does not cover the current final set")
    retired_ids = [
        str(value).lower().removeprefix("0x")
        for value in policy.get("retired_node_ids", ())
    ]
    if (
        policy.get("action")
        != "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve"
        or len(retired_ids) != len(current_ids)
        or set(retired_ids) != set(current_ids)
    ):
        raise ValueError("final-refresh report does not authorize this exact retirement")
    if final.get("evaluated") is not True or final.get("passed") is not False:
        raise ValueError("final refresh requires an evaluated failed final set")

    transport_failed = final.get("transport_status", {}).get("passed") is False
    noise_eligible = final.get("gates", {}).get("mc_noise_eligible") is True
    physics_failed = noise_eligible and final.get("physics_passed") is False
    capped_noise = bool(policy.get("precision_capped_node_ids"))
    detector = final.get("detector_validation", {})
    detector_failed = bool(
        isinstance(detector, Mapping)
        and detector.get("supplied") is True
        and detector.get("noise_passed") is True
        and detector.get("physics_passed") is False
    )
    if not (transport_failed or physics_failed or capped_noise or detector_failed):
        raise ValueError(
            "final-refresh report has no adjudicated transport, physics, or capped-noise failure"
        )
    return _canonical_json_sha256(report)


def replace_final_set_from_reserve(
    manifest: Mapping[str, Any],
    *,
    failure_report: Mapping[str, Any],
    count: int = 64,
) -> dict[str, Any]:
    """Retire a failed final set and select a disjoint fresh set from Sobol reserve."""

    _validate_manifest_record(manifest)
    failure_report_sha256 = _validate_final_refresh_report(manifest, failure_report)
    if count != 64:
        raise ValueError("CEEFIT-v1 fresh final sets contain exactly 64 nodes")
    updated = deepcopy(dict(manifest))
    state = updated["campaign_state"]
    maximum_final_generations = int(
        updated["allocation"].get("maximum_final_generations", 0)
    )
    if (
        maximum_final_generations != MAXIMUM_FINAL_GENERATIONS
        or int(state["final_generation"]) + 1 >= maximum_final_generations
    ):
        raise ValueError(
            "fresh-final generation cap reached; start a new campaign revision"
        )
    design_nodes = updated["design"]["nodes"]
    old_final = list(design_nodes["final"])
    if len(old_final) != count:
        raise ValueError("the current final set does not contain 64 nodes")
    pool = updated["design"]["candidate_pool"]
    available_records = [pool[int(index)] for index in state["unused_candidate_indices"]]
    if len(available_records) < count:
        raise ValueError("Sobol reserve cannot supply a fresh final set")
    train_nodes = tuple(_physical_node(item) for item in design_nodes["train"])
    validation_nodes = tuple(_physical_node(item) for item in design_nodes["validation"])
    available_nodes = tuple(_physical_node(item) for item in available_records)
    selection = select_disjoint_maximin(
        available_nodes,
        sizes={
            "train": len(train_nodes),
            "validation": len(validation_nodes),
            "final": count,
        },
        forced={"train": train_nodes, "validation": validation_nodes},
        q_domain=(
            updated["design"]["q_domain"]["min"],
            updated["design"]["q_domain"]["max"],
        ),
    )
    config_hash = str(updated["config"]["sha256"])
    new_final = [_encode_node(node, "final", config_hash) for node in selection.final]
    record_by_key = {_coordinate_key(record): record for record in available_records}
    selected_indices = {
        int(record_by_key[_coordinate_key(node)]["candidate_index"])
        for node in selection.final
    }
    for index in selected_indices:
        pool[index]["selected_role"] = "final"
    for record in pool:
        if record.get("selected_role") == "final" and int(record["candidate_index"]) not in selected_indices:
            record["selected_role"] = "diagnostic"
    state["unused_candidate_indices"] = [
        index for index in state["unused_candidate_indices"] if int(index) not in selected_indices
    ]
    state["retired_final_node_ids"].extend(item["node_id_hex"] for item in old_final)
    state["training_node_ids"].extend(item["node_id_hex"] for item in old_final)
    # The failed final truth is now training/diagnostic information, so the
    # fitted model changes. Preserve the completed epoch for audit, then lock
    # the fresh final set behind two new consecutive validation passes.
    state.setdefault("validation_epochs", []).append(
        {
            "final_generation": int(state["final_generation"]),
            "validation_history": list(state.get("validation_history", [])),
            "consecutive_validation_passes": int(
                state.get("consecutive_validation_passes", 0)
            ),
            "stopped": bool(state.get("stopped", False)),
            "retired_final_node_ids": [item["node_id_hex"] for item in old_final],
            "closed_reason": "failed_final_retired_into_training",
            "failure_report_sha256": failure_report_sha256,
            "archive_snapshot_sha256": str(failure_report["archive_snapshot_sha256"]),
            "model_content_sha256": str(failure_report["model_content_sha256"]),
        }
    )
    state["validation_history"] = []
    state["consecutive_validation_passes"] = 0
    state["stopped"] = False
    state.pop("continuation_schedule", None)
    state["final_generation"] = int(state["final_generation"]) + 1
    state["last_action"] = "replace_failed_final_from_unused_sobol_reserve"
    design_nodes["final"] = new_final
    design_nodes["train"].extend(
        [{**item, "role": "diagnostic"} for item in old_final]
    )
    work, groups = _build_work_for_nodes(
        new_final,
        first_block_id=1,
        block_count=int(updated["allocation"]["initial_blocks_per_node"]),
        **_execution_arguments(updated),
    )
    updated["work"].extend(work)
    updated["worker_groups"].extend(groups)
    return updated


def save_manifest(manifest: Mapping[str, Any], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        stream.write(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, destination)


def load_manifest(path: str | Path) -> dict[str, Any]:
    manifest = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate_manifest_record(manifest)
    return manifest


CONTINUATION_SCHEMA = "CEEFIT-CONTINUATION-v1"
ADAPTIVE_UTILITY_SCHEMA = "CEEFIT-ADAPTIVE-UTILITY-v1"


def _load_adaptive_utility_report(
    path: str | Path,
    manifest: Mapping[str, Any],
    *,
    archive_snapshot_sha256: str,
) -> dict[str, Any]:
    """Validate a fit-produced common-currency adaptive decision input."""

    report = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(report, dict) or report.get("schema") != ADAPTIVE_UTILITY_SCHEMA:
        raise ValueError("adaptive utility input is not CEEFIT-ADAPTIVE-UTILITY-v1")
    state = manifest["campaign_state"]
    history = state.get("validation_history", [])
    if not history or history[-1].get("passed") is not False:
        raise ValueError("adaptive allocation requires the latest validation round to fail")
    latest = history[-1]
    expected = {
        "config_sha256": manifest["config"]["sha256"],
        "archive_snapshot_sha256": archive_snapshot_sha256,
        "campaign_round": state["round"],
        "campaign_round_identity_sha256": latest["campaign_round_identity_sha256"],
    }
    for name, value in expected.items():
        if report.get(name) != value:
            raise ValueError(f"adaptive utility {name} is stale or bound to another run")
    expected_report_hash = latest.get("adaptive_utility_report_sha256")
    actual_report_hash = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if not _is_sha256(expected_report_hash) or actual_report_hash != expected_report_hash:
        raise ValueError("adaptive utility bytes differ from the report recorded at validation")
    if latest.get("kind") == "bound_model_selection_report" and report.get(
        "model_content_sha256"
    ) != latest.get("model_content_sha256"):
        raise ValueError("adaptive utility model identity differs from validation")
    if not isinstance(report.get("input_parameters"), Mapping) or not isinstance(
        report.get("output_parameters"), Mapping
    ):
        raise ValueError("adaptive utility report lacks input/output parameters")

    new_nodes = report.get("new_nodes")
    local_losses = report.get("continuation_local_validation_losses")
    if not isinstance(new_nodes, Mapping) or not isinstance(local_losses, Mapping):
        raise ValueError("adaptive utility report lacks both allocation branches")
    normalized_new: dict[str, dict[str, float]] = {}
    for raw_id, raw in new_nodes.items():
        node_id = str(raw_id).lower().removeprefix("0x").zfill(16)
        if (
            len(node_id) != 16
            or any(character not in "0123456789abcdef" for character in node_id)
            or node_id in normalized_new
        ):
            raise ValueError("adaptive new-node identity is invalid or duplicated")
        if not isinstance(raw, Mapping):
            raise ValueError("adaptive new-node utility record is malformed")
        reduction = float(raw.get("predicted_loss_reduction", math.nan))
        cpu = float(raw.get("predicted_cpu_seconds", math.nan))
        utility = float(raw.get("utility", math.nan))
        if not (
            math.isfinite(reduction)
            and reduction >= 0.0
            and math.isfinite(cpu)
            and cpu > 0.0
            and math.isfinite(utility)
            and math.isclose(utility, reduction / cpu, rel_tol=1.0e-12, abs_tol=1.0e-15)
        ):
            raise ValueError("adaptive new-node loss/CPU utility is inconsistent")
        normalized_new[node_id] = {
            "predicted_loss_reduction": reduction,
            "predicted_cpu_seconds": cpu,
            "utility": utility,
        }
    normalized_losses: dict[str, float] = {}
    for raw_id, raw_loss in local_losses.items():
        node_id = str(raw_id).lower().removeprefix("0x").zfill(16)
        if (
            len(node_id) != 16
            or any(character not in "0123456789abcdef" for character in node_id)
            or node_id in normalized_losses
        ):
            raise ValueError("adaptive continuation identity is invalid or duplicated")
        loss = float(raw_loss)
        if not math.isfinite(loss) or loss < 0.0:
            raise ValueError("adaptive continuation local loss is invalid")
        normalized_losses[node_id] = loss
    return {
        **report,
        "new_nodes": normalized_new,
        "continuation_local_validation_losses": normalized_losses,
        "report_sha256": actual_report_hash,
    }


def _active_campaign_nodes(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return the current role overlay without exposing the Sobol reserve."""

    design_nodes = manifest["design"]["nodes"]
    records: dict[str, dict[str, Any]] = {}
    ordered_ids: list[str] = []
    for role in ("train", "validation", "final"):
        for raw in design_nodes.get(role, []):
            record = dict(raw)
            node_hex = str(record["node_id_hex"]).lower().removeprefix("0x").zfill(16)
            if node_hex in records:
                raise ValueError(f"campaign role overlay repeats node {node_hex}")
            if int(record["node_id"]) != int(node_hex, 16):
                raise ValueError(f"campaign node identity disagrees for {node_hex}")
            record["node_id_hex"] = node_hex
            record["role"] = str(record.get("role", role))
            if record["role"] not in {
                "pilot",
                "train",
                "validation",
                "final",
                "diagnostic",
            }:
                raise ValueError(f"campaign node {node_hex} has an invalid active role")
            records[node_hex] = record
            ordered_ids.append(node_hex)

    state = manifest.get("campaign_state", {})
    training_ids = [
        str(value).lower().removeprefix("0x").zfill(16)
        for value in state.get("training_node_ids", ())
    ]
    missing = sorted(set(training_ids) - set(records))
    if missing:
        raise ValueError(
            "campaign training-node state references nodes absent from the design overlay: "
            + ", ".join(missing[:4])
        )
    return [records[node_id] for node_id in ordered_ids]


def _archive_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _archive_json(value: Any) -> Any:
    return json.loads(_archive_text(value))


def _standard_error(count: int, first: float, second: float) -> float:
    if count < 2:
        return math.inf
    numerator = max(0.0, second - first * first / count)
    return math.sqrt(numerator / (count * (count - 1)))


def _replica_uncertainties(
    numpy: Any,
    replica_s1: Any,
    replica_histories: Any,
    truth: Any,
    *,
    native_shape: tuple[int, int],
) -> tuple[float, float]:
    """Return absolute yield SE and sentinel-grid shape-TV uncertainty."""

    raw = numpy.asarray(replica_s1, dtype=float).reshape(32, -1)
    counts = numpy.asarray(replica_histories, dtype=float).reshape(32)
    if numpy.any(counts <= 0.0):
        return math.inf, math.inf
    estimates = raw / counts[:, None]
    truth_array = numpy.asarray(truth, dtype=float).reshape(-1)
    truth_yield = float(numpy.sum(truth_array))
    yields = numpy.sum(estimates, axis=1)
    yield_se = float(numpy.std(yields, ddof=1) / math.sqrt(len(yields)))
    if truth_yield <= 0.0:
        if numpy.any(raw):
            return math.inf, math.inf
        return yield_se, 0.0
    shape_se = replica_shape_tv_noise(
        numpy,
        raw,
        counts,
        truth_array,
        native_shape=native_shape,
    )
    return yield_se, shape_se


def _archive_component_offsets(root: Any, descriptors: Sequence[Mapping[str, Any]]) -> list[int]:
    n_mu = len(root["axes/mu_edges"]) - 1
    offsets = [0]
    for descriptor in descriptors:
        axis = int(descriptor["energy_axis"])
        n_energy = root[f"axes/energy/{axis}/physical_edges_keV"].shape[1] - 1
        offsets.append(offsets[-1] + n_energy * n_mu)
    return offsets


def _cell_means(numpy: Any, count: Any, first: Any) -> Any:
    count_array = numpy.asarray(count, dtype=float)
    first_array = numpy.asarray(first, dtype=float)
    return numpy.divide(
        first_array,
        count_array,
        out=numpy.zeros_like(first_array),
        where=count_array > 0.0,
    )


def _archive_additional_criteria(
    root: Any,
    index: int,
    *,
    total_mean: float,
    histories: int,
    numpy: Any,
) -> tuple[Any, ...]:
    """Build exact band/moment and replica-covariance allocation checks."""

    from .allocation import CriterionResult

    result: list[CriterionResult] = []
    stopping_n = int(root["stopping/N"][index])
    stopping_s1 = numpy.asarray(root["stopping/S1"][index], dtype=float).reshape(16)
    stopping_s2 = numpy.asarray(
        root["stopping/outer_S2"][index], dtype=float
    ).reshape(16, 16)
    if stopping_n != histories:
        raise ValueError("archive stopping-vector N disagrees with node histories")
    for band in range(12):
        mean = stopping_s1[band] / stopping_n if stopping_n else 0.0
        result.append(
            CriterionResult(
                f"registered_band_{band:02d}",
                _standard_error(
                    stopping_n,
                    float(stopping_s1[band]),
                    float(stopping_s2[band, band]),
                ),
                0.001 * total_mean,
                total_mean > 0.0 and mean >= 0.01 * total_mean,
            )
        )
    for moment in range(4):
        slot = 12 + moment
        result.append(
            CriterionResult(
                f"legendre_l{moment + 1}",
                _standard_error(
                    stopping_n,
                    float(stopping_s1[slot]),
                    float(stopping_s2[slot, slot]),
                ),
                0.0025 * total_mean,
                total_mean > 0.0,
            )
        )

    replica_histories = numpy.asarray(root["replicas/histories"][index], dtype=float)
    n_mu = len(root["axes/mu_edges"]) - 1
    total_truth = _cell_means(
        numpy,
        root["continuum/total/N"][index],
        root["continuum/total/S1"][index],
    )
    _, total_shape_se = _replica_uncertainties(
        numpy,
        root["continuum/total/replica_S1"][index],
        replica_histories,
        total_truth,
        native_shape=(total_truth.size // n_mu, n_mu),
    )
    total_shape_resolution = sentinel_shape(total_truth.size // n_mu, n_mu)
    result.append(
        CriterionResult(
            f"joint_continuum_shape:{total_shape_resolution[0]}x{total_shape_resolution[1]}",
            total_shape_se,
            0.0025,
            total_mean > 0.0,
        )
    )

    descriptors = _archive_json(root["components"].attrs["descriptors_json"])
    offsets = _archive_component_offsets(root, descriptors)
    continuum_n = root["continuum/components/N"][index]
    continuum_s1 = root["continuum/components/S1"][index]
    continuum_replica = root["continuum/components/replica_S1"][index]
    n_lines = root["lines/support_mask"].shape[1]
    for component, descriptor in enumerate(descriptors):
        start, stop = offsets[component], offsets[component + 1]
        truth = _cell_means(
            numpy, continuum_n[start:stop], continuum_s1[start:stop]
        )
        component_mean = float(numpy.sum(truth))
        yield_se, shape_se = _replica_uncertainties(
            numpy,
            continuum_replica[:, start:stop],
            replica_histories,
            truth,
            native_shape=((stop - start) // n_mu, n_mu),
        )
        shape_resolution = sentinel_shape((stop - start) // n_mu, n_mu)
        name = str(descriptor["name"])
        result.append(
            CriterionResult(
                f"component_yield:{name}",
                yield_se,
                0.001 * total_mean,
                total_mean > 0.0 and component_mean >= 0.01 * total_mean,
            )
        )
        result.append(
            CriterionResult(
                f"component_shape_{shape_resolution[0]}x{shape_resolution[1]}:{name}",
                shape_se,
                0.0025,
                total_mean > 0.0 and component_mean >= 0.05 * total_mean,
            )
        )

    line_support = numpy.asarray(root["lines/support_mask"][index], dtype=bool)
    total_line_n = numpy.asarray(root["lines/total/N"][index]).reshape(n_lines, n_mu)
    total_line_s1 = numpy.asarray(root["lines/total/S1"][index]).reshape(n_lines, n_mu)
    total_line_replica = numpy.asarray(root["lines/total/replica_S1"][index]).reshape(
        32, n_lines, n_mu
    )
    line_definitions = _archive_json(root["axes"].attrs["line_definitions_json"])
    direct_fluence = float(root["nodes/direct_fluence"][index])
    for line in range(n_lines):
        if not bool(line_support[line]):
            continue
        truth = _cell_means(numpy, total_line_n[line], total_line_s1[line])
        line_mean = float(numpy.sum(truth))
        yield_se, shape_se = _replica_uncertainties(
            numpy,
            total_line_replica[:, line, :],
            replica_histories,
            truth,
            native_shape=(1, n_mu),
        )
        line_shape_resolution = sentinel_shape(1, n_mu)
        descriptor = line_definitions[line]
        name = str(descriptor["name"])
        if int(descriptor["post_creation_order"]) != 3:
            raise ValueError(f"exact line {name} lacks sentinel order 3")
        result.extend(
            (
                CriterionResult(
                    f"component_yield:exact_line:{name}",
                    yield_se,
                    EXACT_LINE_ABSOLUTE_NOISE_FRACTION * total_mean,
                    total_mean > 0.0 and line_mean >= 0.01 * total_mean,
                ),
                CriterionResult(
                    f"component_shape_{line_shape_resolution[0]}x{line_shape_resolution[1]}:exact_line:{name}",
                    shape_se,
                    EXACT_LINE_RELATIVE_NOISE_FRACTION,
                    total_mean > 0.0 and line_mean >= 0.05 * total_mean,
                ),
            )
        )
        if name == "elastic_E0":
            direct_map = direct_line_map(numpy, direct_fluence, n_mu)
            combined_truth = truth + direct_map
            combined_replicas = numpy.asarray(
                total_line_replica[:, line, :], dtype=float
            ) + replica_histories[:, None] * direct_map[None, :]
            combined_yield_se, combined_shape_se = _replica_uncertainties(
                numpy,
                combined_replicas,
                replica_histories,
                combined_truth,
                native_shape=(1, n_mu),
            )
            combined_mean = direct_fluence + line_mean
            result.extend(
                (
                    CriterionResult(
                        "line_combined_yield:elastic_E0",
                        combined_yield_se,
                        EXACT_LINE_RELATIVE_NOISE_FRACTION * combined_mean,
                        combined_mean > 0.0,
                    ),
                    CriterionResult(
                        f"line_combined_shape_{line_shape_resolution[1]}:elastic_E0",
                        combined_shape_se,
                        EXACT_LINE_RELATIVE_NOISE_FRACTION,
                        combined_mean > 0.0,
                    ),
                )
            )
        else:
            material, _ = exact_line_materiality(line_mean, yield_se, total_mean)
            result.extend(
                (
                    CriterionResult(
                        f"line_absolute_yield:{name}",
                        yield_se,
                        EXACT_LINE_ABSOLUTE_NOISE_FRACTION * total_mean,
                        total_mean > 0.0,
                    ),
                    CriterionResult(
                        f"line_relative_yield:{name}",
                        yield_se,
                        EXACT_LINE_RELATIVE_NOISE_FRACTION * line_mean,
                        material and line_mean > 0.0,
                    ),
                    CriterionResult(
                        f"line_shape_{line_shape_resolution[1]}:{name}",
                        shape_se,
                        EXACT_LINE_RELATIVE_NOISE_FRACTION,
                        material and line_mean > 0.0,
                    ),
                )
            )
    return tuple(result)


def _detector_fold_noise_criteria(
    path: str | Path | None,
    *,
    root: Any,
    manifest: Mapping[str, Any],
) -> tuple[dict[str, tuple[Any, ...]], dict[str, Any]]:
    """Load archive-bound detector replica-noise criteria for final nodes.

    Detector response matrices intentionally remain outside the physics
    archive.  Once the real fold report exists, however, its per-node replica
    uncertainty is an allocation observable: an imprecise final node must get
    more histories before its observed fit residual can retire the final set.
    """

    from .allocation import CriterionResult
    from .publication import archive_snapshot_sha256, file_sha256

    if path is None:
        return {}, {
            "supplied": False,
            "assessed": False,
            "reason": "no archive-bound HPGe/NaI detector-fold report supplied",
        }
    report_path = Path(path)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("cannot parse detector-fold allocation report") from error
    if not isinstance(report, Mapping):
        raise ValueError("detector-fold allocation report must be a JSON object")
    if report.get("schema") != "CEEFIT-DETECTOR-FOLD-VALIDATION-v1":
        raise ValueError("unsupported detector-fold allocation report schema")
    config_hash = str(manifest["config"]["sha256"])
    if report.get("config_sha256") != config_hash:
        raise ValueError("detector-fold allocation report configuration differs")
    snapshot_hash = archive_snapshot_sha256(root)
    if report.get("archive_snapshot_sha256") != snapshot_hash:
        raise ValueError("detector-fold allocation report is stale for this archive")
    model_hash = report.get("model_content_sha256")
    if not _is_sha256(model_hash):
        raise ValueError("detector-fold allocation report lacks a model-content hash")

    expected_ids = [
        str(item["node_id_hex"]).lower()
        for item in manifest["design"]["nodes"]["final"]
    ]
    supplied_ids = [
        str(value).lower().removeprefix("0x")
        for value in report.get("final_node_ids", ())
    ]
    if (
        len(expected_ids) != 64
        or len(set(expected_ids)) != 64
        or supplied_ids != expected_ids
    ):
        raise ValueError(
            "detector-fold allocation report does not bind the exact ordered final set"
        )

    by_node: dict[str, list[CriterionResult]] = {
        node_id: [] for node_id in expected_ids
    }
    response_hashes: dict[str, str] = {}
    for detector_name in ("HPGe", "NaI"):
        detector = report.get("detectors", {}).get(detector_name)
        if not isinstance(detector, Mapping):
            raise ValueError(f"detector-fold allocation report lacks {detector_name}")
        response_hash = detector.get("response_sha256")
        if not _is_sha256(response_hash):
            raise ValueError(f"{detector_name} detector response hash is invalid")
        response_hashes[detector_name] = str(response_hash)
        nodes = detector.get("nodes")
        if not isinstance(nodes, list):
            raise ValueError(f"{detector_name} detector fold lacks node metrics")
        records: dict[str, Mapping[str, Any]] = {}
        for item in nodes:
            if not isinstance(item, Mapping):
                raise ValueError(f"{detector_name} detector node metric is malformed")
            node_id = str(item.get("node_id", "")).lower().removeprefix("0x")
            if node_id in records:
                raise ValueError(f"{detector_name} detector fold repeats node {node_id}")
            records[node_id] = item
        if set(records) != set(expected_ids):
            raise ValueError(f"{detector_name} detector fold misses final nodes")
        for node_id in expected_ids:
            item = records[node_id]
            achieved = float(item.get("replica_relative_se", math.inf))
            yield_se = float(item.get("replica_yield_relative_se", math.inf))
            shape_se = float(item.get("replica_spectral_tv_se", math.inf))
            if not all(
                math.isfinite(value) and value >= 0.0
                for value in (achieved, yield_se, shape_se)
            ):
                raise ValueError(
                    f"{detector_name} detector replica uncertainty is invalid at {node_id}"
                )
            expected = max(yield_se, shape_se)
            if not math.isclose(achieved, expected, rel_tol=1.0e-12, abs_tol=1.0e-15):
                raise ValueError(
                    f"{detector_name} detector replica scalar does not close at {node_id}"
                )
            by_node[node_id].append(
                CriterionResult(
                    f"detector_fold:{detector_name}",
                    achieved,
                    0.0025,
                    True,
                )
            )
    criteria = {node_id: tuple(values) for node_id, values in by_node.items()}
    return criteria, {
        "supplied": True,
        "assessed": True,
        "report_sha256": file_sha256(report_path),
        "archive_snapshot_sha256": snapshot_hash,
        "model_content_sha256": model_hash,
        "response_sha256": response_hashes,
        "all_noise_targets_met": all(
            criterion.passed
            for values in criteria.values()
            for criterion in values
        ),
    }


def _archive_block_ids(
    root: Any,
    node_hex: str,
    *,
    block_histories: int,
    maximum_block_id: int,
) -> tuple[int, ...]:
    """Validate and return the exact logical block inventory for one node."""

    from .archive import node_block_inventory

    for name in root["blocks"].get(node_hex, ()):
        if str(name) != str(int(name)):
            raise ValueError(f"archive node {node_hex} has invalid block ID {name!r}")
    inventory = node_block_inventory(root, node_hex)
    if not inventory:
        return ()
    block_ids: list[int] = []
    histories_per_replica = block_histories // 32
    for block_id, entry in inventory.items():
        if not 1 <= block_id <= maximum_block_id:
            raise ValueError(f"archive node {node_hex} has invalid block ID {block_id!r}")
        if int(entry["histories"]) != block_histories:
            raise ValueError(f"archive node {node_hex} block {block_id} has wrong history count")
        replica_ids = tuple(entry["id"])
        begins = tuple(entry["history_begin"])
        counts = tuple(entry["history_count"])
        if replica_ids != tuple(range(32)) or counts != (histories_per_replica,) * 32:
            raise ValueError(f"archive node {node_hex} block {block_id} has invalid replicas")
        expected_begins = tuple(
            (block_id - 1) * block_histories + replica * histories_per_replica
            for replica in range(32)
        )
        if begins != expected_begins:
            raise ValueError(
                f"archive node {node_hex} block {block_id} has noncanonical history ranges"
            )
        block_ids.append(block_id)
    if len(block_ids) != len(set(block_ids)):
        raise ValueError(f"archive node {node_hex} repeats a block ID")
    return tuple(sorted(block_ids))


def _scheduled_block_ids(manifest: Mapping[str, Any]) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    for item in manifest.get("work", ()):
        node_hex = str(item["node_id_hex"]).lower().removeprefix("0x").zfill(16)
        block_id = int(item["block_id"])
        if block_id <= 0:
            raise ValueError("campaign work contains a non-positive block ID")
        result.setdefault(node_hex, set()).add(block_id)
    continuation = manifest.get("campaign_state", {}).get("continuation_schedule", {})
    for item in continuation.get("pending_work", ()):
        node_hex = str(item["node_id_hex"]).lower().removeprefix("0x").zfill(16)
        result.setdefault(node_hex, set()).add(int(item["block_id"]))
    return result


def _finite_number(value: float) -> float | None:
    return float(value) if math.isfinite(value) else None


def _assessment_record(
    stats: Any,
    assessment: Any,
    *,
    block_ids: Sequence[int],
    archive_role: str | None,
    status: str,
) -> dict[str, Any]:
    return {
        "node_id": stats.node_id,
        "role": stats.role,
        "archive_role": archive_role,
        "archive_status": status,
        "histories": stats.histories,
        "completed_block_ids": list(block_ids),
        "passed": bool(assessment.passed),
        "usable": bool(assessment.usable),
        "stop_reason": assessment.stop_reason,
        "history_budget_exhausted": bool(assessment.history_budget_exhausted),
        "history_budget_capped": bool(assessment.history_budget_capped),
        "transport_failed": bool(assessment.transport_failed),
        "recommended_additional_blocks": assessment.recommended_additional_blocks,
        "worst_ratio": _finite_number(assessment.worst_ratio),
        "criteria": [
            {
                "name": value.name,
                "achieved_se": _finite_number(value.achieved_se),
                "target_se": _finite_number(value.target_se),
                "active": bool(value.active),
                "ratio": _finite_number(value.ratio),
                "passed": bool(value.passed),
            }
            for value in assessment.criteria
        ],
    }


def _contiguous_runs(values: Sequence[int]) -> list[tuple[int, int]]:
    if not values:
        return []
    ordered = sorted(set(int(value) for value in values))
    result: list[tuple[int, int]] = []
    first = previous = ordered[0]
    for value in ordered[1:]:
        if value != previous + 1:
            result.append((first, previous - first + 1))
            first = value
        previous = value
    result.append((first, previous - first + 1))
    return result


def _archive_inventory_sha256(root: Any, node_ids: Sequence[str]) -> str:
    from .archive import node_block_inventory

    records: list[dict[str, Any]] = []
    for node_hex in sorted(node_ids):
        for block_id, entry in node_block_inventory(root, node_hex).items():
            records.append(
                {
                    "node_id": node_hex,
                    "block_id": block_id,
                    "payload_sha256": entry["payload_sha256"],
                }
            )
    encoded = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def schedule_archive_continuation(
    manifest: Mapping[str, Any],
    archive_path: str | Path,
    *,
    max_nodes: int = 64,
    detector_validation_path: str | Path | None = None,
    adaptive_utility_path: str | Path | None = None,
) -> dict[str, Any]:
    """Assess archived uncertainty and issue only the next disjoint node work.

    The archive is the completion authority. Previously issued commands that
    have not been ingested are re-emitted only when the current uncertainty
    assessment still calls for them. Existing precision continuations at most
    double their histories per pass (continuation_block_count), and their
    utility is priced for that many blocks; interrupted nodes receive enough
    blocks to complete the four-block initial allocation.
    """

    from .allocation import (
        NodeStatistics,
        RawMoment,
        assess_node,
        compare_adaptive_batches,
        continuation_block_count,
        continuation_loss_reduction_components,
        history_cap_blocks,
    )
    from .archive import ARCHIVE_SCHEMA, require_hdf5

    _validate_manifest_record(manifest)
    if max_nodes <= 0:
        raise ValueError("archive continuation max_nodes must be positive")
    block_histories = int(manifest["allocation"]["block_histories"])
    if block_histories != BLOCK_HISTORIES or block_histories % 32:
        raise ValueError("campaign continuation requires canonical 262144-history blocks")
    nodes = _active_campaign_nodes(manifest)
    node_by_id = {str(item["node_id_hex"]): item for item in nodes}
    scheduled = _scheduled_block_ids(manifest)
    h5py, numpy = require_hdf5()

    assessments: list[dict[str, Any]] = []
    scheduling_rows: list[
        tuple[float, float, float, str, dict[str, Any], Any, Any, tuple[int, ...], str]
    ] = []
    continuation_utility_components: dict[str, dict[str, float]] = {}
    archived_ids_by_node: dict[str, list[int]] = {}
    with h5py.File(archive_path, "r") as root:
        if _archive_text(root.attrs.get("schema", "")) != ARCHIVE_SCHEMA:
            raise ValueError("continuation input is not a supported CEEFIT HDF5 archive")
        archive_hash = _archive_text(root["manifest"].attrs["config_sha256"])
        if archive_hash != manifest["config"]["sha256"]:
            raise ValueError("campaign and HDF5 archive configuration hashes differ")
        adaptive_report = None
        if adaptive_utility_path is not None:
            from .publication import archive_snapshot_sha256

            adaptive_report = _load_adaptive_utility_report(
                adaptive_utility_path,
                manifest,
                archive_snapshot_sha256=archive_snapshot_sha256(root),
            )
        adaptive_required = (
            manifest.get("campaign_state", {}).get("last_action")
            == "await_unified_adaptive_allocation"
        )
        if adaptive_required and adaptive_report is None:
            raise ValueError(
                "failed validation requires --adaptive-utilities so histories and "
                "new nodes are compared in one decision"
            )
        local_losses = (
            adaptive_report["continuation_local_validation_losses"]
            if adaptive_report is not None
            else {}
        )
        archive_index = {
            f"{int(value):016x}": index for index, value in enumerate(root["nodes/id"][:])
        }
        if len(archive_index) != len(root["nodes/id"]):
            raise ValueError("archive contains duplicate node IDs")
        detector_criteria, detector_status = _detector_fold_noise_criteria(
            detector_validation_path,
            root=root,
            manifest=manifest,
        )

        for node in nodes:
            node_hex = str(node["node_id_hex"])
            role = str(node["role"])
            cap_blocks = history_cap_blocks(role, block_histories)
            index = archive_index.get(node_hex)
            archive_role: str | None = None
            status = "not_archived"
            command_role = role
            if index is None:
                block_ids: tuple[int, ...] = ()
                stats = NodeStatistics(
                    node_hex,
                    role,
                    0,
                    RawMoment(0.0, 0.0, 0),
                    completed_blocks=0,
                )
            else:
                coordinates = tuple(float(value) for value in root["nodes/coordinates"][index])
                expected = (
                    float(node["energy_keV"]),
                    float(node["distance_m"]) * 100.0,
                    float(node["source_height_m"]) * 100.0,
                    float(node["detector_height_m"]) * 100.0,
                )
                if any(
                    not math.isclose(actual, wanted, rel_tol=1.0e-13, abs_tol=1.0e-12)
                    for actual, wanted in zip(coordinates[:4], expected)
                ):
                    raise ValueError(f"archive coordinates disagree for campaign node {node_hex}")
                archive_role = _archive_text(root["nodes/role"].asstr()[index])
                if archive_role != role and not (
                    role == "diagnostic" and archive_role in {"train", "final"}
                ):
                    raise ValueError(
                        f"archive role {archive_role!r} disagrees with campaign role "
                        f"{role!r} for node {node_hex}"
                    )
                command_role = archive_role
                inventory_cap_blocks = max(
                    cap_blocks,
                    history_cap_blocks(archive_role, block_histories),
                )
                block_ids = _archive_block_ids(
                    root,
                    node_hex,
                    block_histories=block_histories,
                    maximum_block_id=inventory_cap_blocks,
                )
                histories = int(root["nodes/histories"][index])
                if histories != len(block_ids) * block_histories:
                    raise ValueError(
                        f"archive history total disagrees with block inventory for node {node_hex}"
                    )
                status = _archive_text(root["nodes/status"].asstr()[index])
                transport_failed = bool(
                    status not in {"complete", "converged"}
                    or int(root["nodes/transport_failures"][index]) != 0
                    or bool(root["nodes/capped"][index])
                )
                total = RawMoment(
                    float(root["diagnostics/scattered_total_S1"][index]),
                    float(root["diagnostics/scattered_total_S2"][index]),
                    int(root["diagnostics/scattered_total_N"][index]),
                )
                if total.histories != histories:
                    raise ValueError(
                        f"archive total moment N disagrees with histories for node {node_hex}"
                    )
                total_mean = abs(total.mean) if histories else 0.0
                additional = _archive_additional_criteria(
                    root,
                    index,
                    total_mean=total_mean,
                    histories=histories,
                    numpy=numpy,
                )
                additional = (
                    *additional,
                    *detector_criteria.get(node_hex, ()),
                )
                timing = float(root["nodes/cpu_seconds"][index])
                if not math.isfinite(timing) or timing <= 0.0:
                    timing = float(root["nodes/wall_seconds"][index])
                stats = NodeStatistics(
                    node_hex,
                    role,
                    histories,
                    total,
                    wall_seconds=max(0.0, timing),
                    completed_blocks=len(block_ids),
                    transport_failed=transport_failed,
                    additional_criteria=additional,
                )
            archived_ids_by_node[node_hex] = list(block_ids)
            assessment = assess_node(
                stats,
                block_histories=block_histories,
                # A failed final set can legitimately contain more histories
                # than the training/diagnostic cap before it is reclassified.
                # Preserve that valid truth, mark its new budget exhausted,
                # and never allocate another block to it.
                allow_existing_over_cap=(
                    archive_role is not None
                    and role == "diagnostic"
                    and len(block_ids) > cap_blocks
                ),
            )
            assessments.append(
                _assessment_record(
                    stats,
                    assessment,
                    block_ids=block_ids,
                    archive_role=archive_role,
                    status=status,
                )
            )
            if assessment.recommended_additional_blocks <= 0:
                continue
            observed_blocks = len(block_ids)
            issue_blocks = continuation_block_count(
                observed_blocks, assessment.recommended_additional_blocks
            )
            if observed_blocks < MIN_INITIAL_BLOCKS:
                priority = math.inf
                predicted_reduction = math.inf
                predicted_cpu = max(1.0, stats.wall_seconds / max(1, observed_blocks))
            else:
                seconds_per_block = (
                    stats.wall_seconds / observed_blocks
                    if stats.wall_seconds > 0.0 and observed_blocks
                    else 1.0
                )
                if adaptive_required and role not in {"train", "validation"}:
                    continue
                if adaptive_report is not None:
                    if node_hex not in local_losses:
                        if adaptive_required and role in {"train", "validation"}:
                            raise ValueError(
                                f"adaptive utility report lacks active continuation node {node_hex}"
                            )
                        continue
                    reduction = continuation_loss_reduction_components(
                        float(local_losses[node_hex]),
                        assessment,
                        observed_blocks,
                        issue_blocks,
                    )
                    predicted_reduction = reduction.total
                    continuation_utility_components[node_hex] = reduction.as_dict()
                    predicted_cpu = seconds_per_block * issue_blocks
                    priority = predicted_reduction / predicted_cpu
                else:
                    # Precision-only scheduling remains available for final-set
                    # noise resolution, where no unused final set may be inspected
                    # or compared until the adaptive campaign has stopped.
                    predicted_reduction = max(0.0, assessment.worst_ratio - 1.0)
                    predicted_cpu = seconds_per_block * issue_blocks
                    priority = predicted_reduction / seconds_per_block
            scheduling_rows.append(
                (
                    priority,
                    predicted_reduction,
                    predicted_cpu,
                    node_hex,
                    node,
                    stats,
                    assessment,
                    block_ids,
                    command_role,
                )
            )

        inventory_hash = _archive_inventory_sha256(root, tuple(node_by_id))

    scheduling_rows.sort(
        key=lambda item: (
            -item[0] if math.isfinite(item[0]) else -math.inf,
            item[3],
        )
    )
    selected = scheduling_rows[:max_nodes]
    updated = deepcopy(dict(manifest))
    unified_decision: dict[str, Any] = {
        "required": adaptive_required,
        "selected_branch": "precision_only_continuation",
        "adaptive_utility_report_sha256": (
            adaptive_report.get("report_sha256") if adaptive_report is not None else None
        ),
    }
    selected_new_records: list[dict[str, Any]] = []
    if adaptive_required:
        batch_size = int(manifest["allocation"]["acquisition_batch_size"])
        if max_nodes != batch_size:
            raise ValueError(
                f"unified adaptive allocation requires exactly {batch_size} node jobs per decision"
            )
        forced_initial = any(not math.isfinite(row[0]) for row in selected)
        forced_validation_truth_precision = any(
            str(row[4]["role"]) == "validation" for row in scheduling_rows
        )
        forced_continuation = forced_initial or forced_validation_truth_precision
        continuation_rows = [row for row in scheduling_rows if math.isfinite(row[0])][
            :batch_size
        ]

        state_record = manifest["campaign_state"]
        maximum = int(manifest["allocation"]["maximum_training_nodes"])
        capacity = maximum - len(state_record["training_node_ids"])
        if 0 < capacity < batch_size:
            raise ValueError("training-node cap would require a partial acquisition batch")
        available = {
            int(index): manifest["design"]["candidate_pool"][int(index)]
            for index in state_record["unused_candidate_indices"]
        }
        new_predictions = adaptive_report["new_nodes"]
        eligible_new = [
            record
            for record in available.values()
            if str(record["node_id_hex"]) in new_predictions
        ]
        eligible_new.sort(
            key=lambda record: (
                -new_predictions[str(record["node_id_hex"])]["utility"],
                int(record["candidate_index"]),
            )
        )
        new_batch = eligible_new[:batch_size] if capacity >= batch_size else []
        if capacity >= batch_size and len(new_batch) != batch_size:
            raise ValueError("adaptive report does not cover a complete 64-node candidate batch")
        branch_decision = compare_adaptive_batches(
            [(row[1], row[2]) for row in continuation_rows],
            [
                (
                    new_predictions[str(record["node_id_hex"])]["predicted_loss_reduction"],
                    new_predictions[str(record["node_id_hex"])]["predicted_cpu_seconds"],
                )
                for record in new_batch
            ],
            forced_continuation=forced_continuation,
        )
        choose_new = branch_decision.selected_branch == "new_node_batch"
        unified_decision.update(
            {
                "criterion": "predicted_validation_loss_reduction_per_cpu_second",
                "continuation": {
                    "action_count": branch_decision.continuation.action_count,
                    "predicted_loss_reduction": branch_decision.continuation.predicted_loss_reduction,
                    "predicted_cpu_seconds": branch_decision.continuation.predicted_cpu_seconds,
                    "utility": branch_decision.continuation.utility,
                    "node_components": {
                        row[3]: continuation_utility_components[row[3]]
                        for row in continuation_rows
                        if row[3] in continuation_utility_components
                    },
                },
                "new_node_batch": {
                    "action_count": branch_decision.new_node_batch.action_count,
                    "predicted_loss_reduction": branch_decision.new_node_batch.predicted_loss_reduction,
                    "predicted_cpu_seconds": branch_decision.new_node_batch.predicted_cpu_seconds,
                    "utility": branch_decision.new_node_batch.utility,
                },
                "forced_initial_block_completion": forced_initial,
                "forced_validation_truth_precision": forced_validation_truth_precision,
                "selected_branch": "new_node_batch" if choose_new else "continuation",
            }
        )
        if choose_new:
            selected = []
            selected_new_records = new_batch
        elif forced_continuation:
            selected = sorted(
                scheduling_rows,
                key=lambda row: (
                    0
                    if not math.isfinite(row[0])
                    else 1
                    if str(row[4]["role"]) == "validation"
                    else 2,
                    -row[0] if math.isfinite(row[0]) else 0.0,
                    row[3],
                ),
            )[:batch_size]
        else:
            selected = continuation_rows
    pending_ids_by_node: dict[str, list[int]] = {}
    pending_work: list[dict[str, Any]] = []
    pending_groups: list[dict[str, Any]] = []
    cumulative_keys = {
        (str(item["node_id_hex"]), int(item["block_id"]))
        for item in manifest.get("work", ())
    }
    new_cumulative_work: list[dict[str, Any]] = []
    new_cumulative_groups: list[dict[str, Any]] = []
    execution = _execution_arguments(manifest)
    for _, _, _, node_hex, node, stats, assessment, archived, command_role in selected:
        current_blocks = len(archived)
        issue_count = continuation_block_count(
            current_blocks, int(assessment.recommended_additional_blocks)
        )
        cap_blocks = history_cap_blocks(str(node["role"]), block_histories)
        already_scheduled = set(scheduled.get(node_hex, ()))
        outstanding = sorted(already_scheduled - set(archived))
        chosen = outstanding[:issue_count]
        unavailable = set(archived) | already_scheduled
        next_id = 1
        while len(chosen) < issue_count:
            while next_id in unavailable and next_id <= cap_blocks:
                next_id += 1
            if next_id > cap_blocks:
                raise RuntimeError(
                    f"no disjoint block ID remains below the history cap for node {node_hex}"
                )
            chosen.append(next_id)
            unavailable.add(next_id)
        pending_ids_by_node[node_hex] = chosen
        command_node = {**node, "role": command_role}
        for first, count in _contiguous_runs(chosen):
            run_work, run_groups = _build_work_for_nodes(
                [command_node],
                first_block_id=first,
                block_count=count,
                **execution,
            )
            pending_work.extend(run_work)
            pending_groups.extend(run_groups)
        for block_id in chosen:
            key = (node_hex, block_id)
            if key in cumulative_keys:
                continue
            one_work, one_group = _build_work_for_nodes(
                [command_node],
                first_block_id=block_id,
                block_count=1,
                **execution,
            )
            new_cumulative_work.extend(one_work)
            new_cumulative_groups.extend(one_group)
            cumulative_keys.add(key)

    if selected_new_records:
        new_nodes = [
            {
                key: record[key]
                for key in (
                    "energy_keV",
                    "distance_m",
                    "source_height_m",
                    "detector_height_m",
                    "q",
                    "node_id",
                    "node_id_hex",
                )
            }
            | {"role": "train"}
            for record in selected_new_records
        ]
        chosen_indices = {
            int(record["candidate_index"]) for record in selected_new_records
        }
        pool = updated["design"]["candidate_pool"]
        for record in selected_new_records:
            pool[int(record["candidate_index"])]["selected_role"] = "train"
        state_for_acquisition = updated["campaign_state"]
        state_for_acquisition["unused_candidate_indices"] = [
            index
            for index in state_for_acquisition["unused_candidate_indices"]
            if int(index) not in chosen_indices
        ]
        chosen_ids = [str(record["node_id_hex"]) for record in selected_new_records]
        state_for_acquisition["training_node_ids"].extend(chosen_ids)
        state_for_acquisition["acquired_node_ids"].extend(chosen_ids)
        state_for_acquisition["acquisition_history"].append(
            {
                "round": state_for_acquisition["round"],
                "criterion": "unified_predicted_validation_loss_reduction_per_cpu_second",
                "selected_branch": "new_node_batch",
                "node_ids": chosen_ids,
                "branch_comparison": unified_decision,
            }
        )
        updated["design"]["nodes"]["train"].extend(new_nodes)
        initial_blocks = int(updated["allocation"]["initial_blocks_per_node"])
        run_work, run_groups = _build_work_for_nodes(
            new_nodes,
            first_block_id=1,
            block_count=initial_blocks,
            **execution,
        )
        pending_work.extend(run_work)
        pending_groups.extend(run_groups)
        new_cumulative_work.extend(run_work)
        new_cumulative_groups.extend(run_groups)
        for node in new_nodes:
            pending_ids_by_node[str(node["node_id_hex"])] = list(
                range(1, initial_blocks + 1)
            )

    pending_signature = [
        [str(item["node_id_hex"]), int(item["block_id"])] for item in pending_work
    ]
    updated.setdefault("work", []).extend(new_cumulative_work)
    updated.setdefault("worker_groups", []).extend(new_cumulative_groups)
    state = updated.setdefault("campaign_state", {})
    previous = state.get("continuation_schedule", {})
    assessment_identity = hashlib.sha256(
        json.dumps(
            {
                "archive_inventory_sha256": inventory_hash,
                "detector_report_sha256": detector_status.get("report_sha256"),
                "adaptive_utility_report_sha256": unified_decision.get(
                    "adaptive_utility_report_sha256"
                ),
                "unified_decision": unified_decision,
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    unchanged = bool(
        previous.get("assessment_identity_sha256") == assessment_identity
        and previous.get("pending_block_keys") == pending_signature
    )
    generation = int(previous.get("generation", 0)) + (0 if unchanged else 1)
    schedule = {
        "schema": CONTINUATION_SCHEMA,
        "generation": generation,
        "archive_inventory_sha256": inventory_hash,
        "assessment_identity_sha256": assessment_identity,
        "archived_block_ids_by_node": archived_ids_by_node,
        "pending_block_ids_by_node": pending_ids_by_node,
        "pending_block_keys": pending_signature,
        "pending_work": pending_work,
        "pending_worker_groups": pending_groups,
        "assessment_count": len(assessments),
        "assessments": assessments,
        "detector_fold_status": detector_status,
        "unified_allocation_decision": unified_decision,
    }
    state["continuation_schedule"] = schedule
    if not unchanged:
        state.setdefault("continuation_history", []).append(
            {
                "generation": generation,
                "archive_inventory_sha256": inventory_hash,
                "assessment_identity_sha256": assessment_identity,
                "detector_report_sha256": detector_status.get("report_sha256"),
                "adaptive_utility_report_sha256": unified_decision.get(
                    "adaptive_utility_report_sha256"
                ),
                "selected_branch": unified_decision.get("selected_branch"),
                "pending_block_keys": pending_signature,
                "assessed_nodes": len(assessments),
                "scheduled_nodes": len(pending_ids_by_node),
                "scheduled_blocks": len(pending_work),
            }
        )
    state["last_action"] = (
        "archive_continuation_unchanged"
        if unchanged
        else "schedule_adaptive_new_node_batch"
        if pending_work and unified_decision.get("selected_branch") == "new_node_batch"
        else "schedule_archive_continuation"
        if pending_work
        else "archive_precision_complete_or_capped"
    )
    return updated


def continuation_shell_commands(
    manifest: Mapping[str, Any],
    *,
    roles: Iterable[str] | None = None,
    grouped: bool = True,
) -> tuple[str, ...]:
    """Return only the latest archive-reconciled continuation commands."""

    schedule = manifest.get("campaign_state", {}).get("continuation_schedule")
    if not isinstance(schedule, Mapping) or schedule.get("schema") != CONTINUATION_SCHEMA:
        raise ValueError("campaign manifest has no archive continuation schedule")
    allowed = set(roles) if roles is not None else None
    source = schedule["pending_worker_groups" if grouped else "pending_work"]
    return tuple(
        str(item["shell_command"])
        for item in source
        if allowed is None or str(item["role"]) in allowed
    )
