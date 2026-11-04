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

"""Publication evidence and immutable archive-snapshot identities.

Fitting a model and publishing it are deliberately separate operations.  This
module binds external validation reports to the exact fitted model, archive
block inventory, configuration, and untouched final-node set.  A missing or
stale report therefore leaves a payload in candidate state.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import struct
from pathlib import Path
from typing import Any, Mapping, Sequence


PUBLICATION_EVIDENCE_SCHEMA = "CEEFIT-PUBLICATION-EVIDENCE-v1"
PRODUCTION_CONFIGURATION_SCHEMA = "CEEFIT-PRODUCTION-CONFIGURATION-v1"
_GEANT4_REQUIRED_PHYSICS_LIST = "G4EmLivermorePhysics"
_GEANT4_COMPARISON_SIGMA = 2.0
_GEANT4_MINIMUM_FIELD_WORLD_HALF_EXTENT_M = 60.0
_GEANT4_REFERENCE_ENERGY_BINS = 48
_GEANT4_REFERENCE_MU_BINS = 32
_GEANT4_REQUIRED_OBSERVABLE = "scattered-track-length-flux-per-primary-v1"
_GEANT4_REQUIRED_SAMPLING = "gps-isotropic-4pi"
_GEANT4_TOTAL_MOMENT_SCHEMA = "per-event-track-length-v1"
_CALIBRATION_MINIMUM_COMPARISONS = 8
_CALIBRATION_MOMENT_FAMILIES = (
    "continuum_components",
    "continuum_total",
    "line_components",
    "line_total",
)
_MANDATORY_CHECKS = (
    "registered_e0_gate",
    "uncertainty_calibration",
    "geant4_anchor_agreement",
    "python_cpp_parity",
)
_PRODUCTION_AXIS_NAMES = (
    "primary_registered",
    "annihilation_511_registered",
    "brems_log_relative",
)
_PRODUCTION_COMPONENTS = tuple(
    (f"{origin_name}_{order_name}", origin, order, origin)
    for origin, origin_name in enumerate(("primary", "annihilation", "brems"))
    for order, order_name in enumerate(("order0", "order1", "order2plus"))
)
_PRODUCTION_REQUIRED_TRUE_SWITCHES = (
    "include_annihilation",
    "include_brems",
    "include_rayleigh",
    "positron_range",
    "walk_annihilation",
    "walk_brems",
)
_PRODUCTION_ESTIMATOR_FLAGS = (
    "next-event;rayleigh=on;brems=on;annihilation=on"
)
_PRODUCTION_CANDIDATE_COUNT = 10_000
_PRODUCTION_INITIAL_DESIGN_COUNTS = {
    "train": 256,
    "validation": 64,
    "final": 64,
}
_PRODUCTION_INITIAL_BLOCKS_PER_NODE = 4
_PRODUCTION_ALLOCATION = {
    "block_histories": 262_144,
    "logical_replicas": 32,
    "initial_blocks_per_node": _PRODUCTION_INITIAL_BLOCKS_PER_NODE,
    "total_relative_se_target": 0.0025,
    "material_contribution_fraction": 0.01,
    "contribution_absolute_se_fraction": 0.001,
    "sentinel_absolute_se_fraction": 0.0025,
    "training_history_cap": 64_000_000,
    "validation_history_cap": 128_000_000,
    "acquisition_batch_size": 64,
    "maximum_training_nodes": 2_048,
    "maximum_final_generations": 3,
}
_PRODUCTION_BASELINE_PHYSICS = {
    "brems_electron_floor_keV": 150.0,
    "doppler": False,
    "max_vertices": 50,
    "roulette_tau": 8.0,
    "weight_floor": 1.0e-5,
}
_PRODUCTION_SCORING_RADIUS_CM = 4.0

_SOURCE_REPORT_CONTRACTS = {
    "registered_e0_gate": {
        "schema": "CEEFIT-REGISTERED-E0-GATE-v1",
        "bindings": ("config_sha256",),
        "fields": (
            "predefined_index_bands",
            "statistically_empty_bands_included",
            "shape_coverage_fraction",
            "max_integrated_relative_error",
            "max_joint_tv",
            "max_energy_marginal_tv",
            "max_angular_marginal_tv",
            "all_node_stopping_targets_met",
            "generator",
        ),
    },
    "uncertainty_calibration": {
        "schema": "CEEFIT-UNCERTAINTY-CALIBRATION-v1",
        "bindings": ("config_sha256", "archive_snapshot_sha256"),
        "fields": (
            "repeated_null_runs",
            "s2_passed",
            "coarse_covariance_passed",
            "replica_covariance_passed",
            "family_coverage_passed",
            "criteria",
        ),
    },
    "geant4_anchor_agreement": {
        "schema": "CEEFIT-GEANT4-ANCHOR-AGREEMENT-v1",
        "bindings": ("config_sha256", "archive_snapshot_sha256"),
        "fields": (
            "anchor_count",
            "within_combined_uncertainty",
            "physics_list",
            "comparison_sigma",
            "minimum_field_world_half_extent_m",
            "manifest_path",
            "manifest_sha256",
            "reference_file_count",
            "reference_event_count",
            "max_abs_standardized_difference",
            "max_absolute_relative_difference",
            "anchors",
        ),
    },
    "python_cpp_parity": {
        "schema": "CEEFIT-PYTHON-CPP-PARITY-v1",
        "bindings": ("config_sha256", "model_content_sha256"),
        "fields": (
            "candidate_payload_sha256",
            "artifact_paths",
            "protocol",
            "request_sha256",
            "cpp_result_sha256",
            "tolerance",
            "query_count",
            "coverage",
            "compared_value_count",
            "max_absolute_difference",
        ),
    },
}


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def file_sha256(path: str | Path) -> str:
    """Hash the exact bytes of a detailed validation report or data artifact."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _hdf_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def production_configuration_gate(root: Any, campaign: Mapping[str, Any]) -> dict[str, Any]:
    """Assess whether an archive/campaign pair is eligible for production.

    The block and archive readers intentionally accept reduced axes for small,
    fast diagnostic fixtures.  That flexibility must stop at publication: an
    accepted payload needs the fit campaign's 384-cell component-native energy
    axes, the 10-keV transport cutoff, the baseline transport configuration,
    the formal 10,000-point adaptive design, four initial blocks per node, and
    a valid 64-bin pilot or identity-bound 64-to-32 decision.

    Adaptive acquisition and failed-final retirement may grow the current
    training overlay.  Its initial size is recovered from the identity-tracked
    acquired and retired node sets rather than requiring the current training
    set to remain at 256 nodes.
    """

    issues: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            issues.append(message)

    configured_mu_bins: int | None = None
    scoring_radius_cm: float | None = None
    candidate_count: int | None = None
    initial_design_counts: dict[str, int | None] = {
        "train": None,
        "validation": None,
        "final": None,
    }
    initial_blocks_per_node: int | None = None
    try:
        input_schema = _hdf_text(root.attrs.get("schema", ""))
        require(
            input_schema in {"CEEFIT-HDF5-v1", "CEEFIT-FIT-CACHE-v1"},
            "input schema is neither a canonical archive nor a fit cache",
        )
        if input_schema == "CEEFIT-FIT-CACHE-v1":
            # A cache is eligible only when its stat-bound sidecar carries a
            # previously verified canonical snapshot through the fit.
            from .storage import fit_cache_source_snapshot

            fit_cache_source_snapshot(root)
        manifest_attrs = root["manifest"].attrs
        archive_hash = _hdf_text(manifest_attrs["config_sha256"])
        canonical_json = _hdf_text(manifest_attrs["canonical_config_json"])
        require(_is_sha256(archive_hash), "archive configuration hash is invalid")
        require(
            hashlib.sha256(canonical_json.encode("utf-8")).hexdigest() == archive_hash,
            "archive canonical configuration checksum differs",
        )
        configuration = json.loads(canonical_json)
        campaign_config = campaign["config"]
        require(
            campaign_config.get("sha256") == archive_hash,
            "campaign and archive configuration hashes differ",
        )
        require(
            campaign_config.get("canonical_json") == canonical_json,
            "campaign and archive canonical configurations differ",
        )

        from .campaign import (
            PhysicsConfig,
            _validated_mu_resolution_decision,
            canonical_fit_config_json,
            default_forced_anchors,
        )
        from .design import (
            DesignBounds,
            PhysicalNode,
            generate_candidate_pool,
            geometry_transform,
            select_disjoint_maximin,
            transformed_domain,
        )

        design = campaign["design"]
        design_nodes = design["nodes"]
        campaign_state = campaign["campaign_state"]
        allocation = campaign["allocation"]
        raw_bounds = design["bounds"]
        production_bounds = DesignBounds(
            energy_keV=tuple(float(value) for value in raw_bounds["energy_keV"]),
            distance_m=tuple(float(value) for value in raw_bounds["distance_m"]),
            source_height_m=tuple(
                float(value) for value in raw_bounds["source_height_m"]
            ),
            detector_height_m=tuple(
                float(value) for value in raw_bounds["detector_height_m"]
            ),
            rho_min_m=float(raw_bounds["rho_min_m"]),
        )
        production_bounds.validate()
        require(
            production_bounds == DesignBounds(),
            "production design bounds differ from the validated campaign domain",
        )
        expected_q_domain = transformed_domain(
            production_bounds, mode=str(campaign_config["mode"])
        )
        raw_q_domain = design["q_domain"]
        q_domain_matches = all(
            len(raw_q_domain[name]) == len(expected)
            and all(
                not isinstance(actual, bool) and float(actual) == target
                for actual, target in zip(raw_q_domain[name], expected)
            )
            for name, expected in zip(("min", "max"), expected_q_domain)
        )
        require(
            q_domain_matches,
            "production transformed q domain is not canonical for its bounds and tier",
        )
        sobol_seed = design["sobol_seed"]
        require(
            type(sobol_seed) is int,
            "production Sobol seed must be an integer",
        )
        candidate_count = design.get("candidate_count")
        require(
            type(candidate_count) is int
            and candidate_count == _PRODUCTION_CANDIDATE_COUNT,
            "production requires exactly 10000 Sobol candidates",
        )
        candidate_pool = design.get("candidate_pool")
        candidate_pool_is_sized = (
            isinstance(candidate_pool, Sequence)
            and not isinstance(candidate_pool, (str, bytes))
            and len(candidate_pool) == _PRODUCTION_CANDIDATE_COUNT
        )
        require(
            candidate_pool_is_sized,
            "production candidate pool must contain exactly 10000 records",
        )

        selection_record = design.get("selection")
        selection_record_valid = isinstance(selection_record, Mapping)
        if selection_record_valid:
            selection_without_hash = {
                key: value
                for key, value in selection_record.items()
                if key != "canonical_sha256"
            }
            try:
                selection_record_valid = bool(
                    selection_record.get("method") == "greedy-maximin-q-v1"
                    and selection_record.get("canonical_sha256")
                    == hashlib.sha256(_canonical(selection_without_hash)).hexdigest()
                )
            except (TypeError, ValueError):
                selection_record_valid = False
        require(
            selection_record_valid,
            "production design lacks checksummed maximin-selection provenance",
        )

        train_records = design_nodes["train"]
        validation_records = design_nodes["validation"]
        final_records = design_nodes["final"]
        train_ids = [str(item["node_id_hex"]).lower() for item in train_records]
        acquired_ids = [
            str(value).lower() for value in campaign_state["acquired_node_ids"]
        ]
        retired_final_ids = [
            str(value).lower()
            for value in campaign_state["retired_final_node_ids"]
        ]
        training_state_ids = [
            str(value).lower() for value in campaign_state["training_node_ids"]
        ]
        require(
            len(train_ids) == len(set(train_ids)),
            "production training overlay contains duplicate node IDs",
        )
        require(
            len(acquired_ids) == len(set(acquired_ids)),
            "production acquired-node history contains duplicate node IDs",
        )
        require(
            len(retired_final_ids) == len(set(retired_final_ids)),
            "production retired-final history contains duplicate node IDs",
        )
        require(
            not (set(acquired_ids) & set(retired_final_ids)),
            "production acquired and retired-final node histories overlap",
        )
        require(
            set(training_state_ids) == set(train_ids)
            and len(training_state_ids) == len(train_ids),
            "campaign training-node state differs from the design overlay",
        )
        require(
            set(acquired_ids).issubset(train_ids)
            and set(retired_final_ids).issubset(train_ids),
            "adaptive or retired training nodes are absent from the design overlay",
        )
        acquisition_history = campaign_state["acquisition_history"]
        flattened_acquisitions: list[str] = []
        acquisition_batches_are_canonical = (
            isinstance(acquisition_history, Sequence)
            and not isinstance(acquisition_history, (str, bytes))
        )
        if acquisition_batches_are_canonical:
            for record in acquisition_history:
                if not isinstance(record, Mapping):
                    acquisition_batches_are_canonical = False
                    break
                batch_ids = [str(value).lower() for value in record.get("node_ids", ())]
                if len(batch_ids) != 64 or len(set(batch_ids)) != 64:
                    acquisition_batches_are_canonical = False
                    break
                flattened_acquisitions.extend(batch_ids)
        require(
            acquisition_batches_are_canonical
            and flattened_acquisitions == acquired_ids
            and len(acquired_ids) % 64 == 0,
            "production acquisitions must be identity-tracked exact 64-node batches",
        )
        require(
            len(retired_final_ids) % 64 == 0
            and len(retired_final_ids) == 64 * int(campaign_state["final_generation"]),
            "production retired-final history must contain one 64-node batch per generation",
        )
        adaptive_ids = set(acquired_ids) | set(retired_final_ids)
        initial_design_counts = {
            "train": sum(node_id not in adaptive_ids for node_id in train_ids),
            "validation": len(validation_records),
            "final": len(final_records),
        }
        for role, expected_count in _PRODUCTION_INITIAL_DESIGN_COUNTS.items():
            require(
                initial_design_counts[role] == expected_count,
                f"production initial {role} design requires exactly {expected_count} nodes",
            )

        initial_blocks_per_node = allocation.get("initial_blocks_per_node")
        for name, expected_value in _PRODUCTION_ALLOCATION.items():
            actual_value = allocation.get(name)
            if isinstance(expected_value, int):
                allocation_matches = (
                    type(actual_value) is int and actual_value == expected_value
                )
            else:
                allocation_matches = (
                    not isinstance(actual_value, bool)
                    and float(actual_value) == expected_value
                )
            require(
                allocation_matches,
                f"production allocation {name} must equal {expected_value!r}",
            )

        role_by_node_id: dict[str, str] = {}
        overlay_is_disjoint = True
        overlay_records_are_authentic = True
        mode = str(campaign_config["mode"])
        for container_role, records in (
            ("train", train_records),
            ("validation", validation_records),
            ("final", final_records),
        ):
            for record in records:
                node_id = str(record["node_id_hex"]).lower()
                record_role = str(record.get("role", container_role))
                allowed_roles = (
                    {"train", "diagnostic"}
                    if container_role == "train"
                    else {container_role}
                )
                if record_role not in allowed_roles or node_id in role_by_node_id:
                    overlay_is_disjoint = False
                role_by_node_id[node_id] = record_role
                try:
                    physical_node = PhysicalNode(
                        float(record["energy_keV"]),
                        float(record["distance_m"]),
                        float(record["source_height_m"]),
                        float(record["detector_height_m"]),
                    )
                    physical_node.validate()
                    expected_node_id = physical_node.stable_id_hex(archive_hash)
                    expected_q = geometry_transform(physical_node)
                    if mode == "air":
                        coordinate_is_in_domain = (
                            physical_node.source_height_m == 0.0
                            and physical_node.detector_height_m == 0.0
                            and production_bounds.energy_keV[0]
                            <= physical_node.energy_keV
                            <= production_bounds.energy_keV[1]
                            and production_bounds.distance_m[0]
                            <= physical_node.distance_m
                            <= production_bounds.distance_m[1]
                        )
                    else:
                        rho = math.sqrt(
                            physical_node.distance_m * physical_node.distance_m
                            - (
                                physical_node.source_height_m
                                - physical_node.detector_height_m
                            )
                            ** 2
                        )
                        coordinate_is_in_domain = (
                            production_bounds.energy_keV[0]
                            <= physical_node.energy_keV
                            <= production_bounds.energy_keV[1]
                            and production_bounds.distance_m[0]
                            <= physical_node.distance_m
                            <= production_bounds.distance_m[1]
                            and production_bounds.source_height_m[0]
                            <= physical_node.source_height_m
                            <= production_bounds.source_height_m[1]
                            and production_bounds.detector_height_m[0]
                            <= physical_node.detector_height_m
                            <= production_bounds.detector_height_m[1]
                            and rho >= production_bounds.rho_min_m
                        )
                    record_is_authentic = (
                        coordinate_is_in_domain
                        and not any(
                            isinstance(record[name], bool)
                            for name in (
                                "energy_keV",
                                "distance_m",
                                "source_height_m",
                                "detector_height_m",
                            )
                        )
                        and int(record["node_id"]) == int(expected_node_id, 16)
                        and node_id == expected_node_id
                        and len(record["q"]) == 4
                        and all(
                            not isinstance(actual, bool) and float(actual) == target
                            for actual, target in zip(record["q"], expected_q)
                        )
                    )
                except (KeyError, TypeError, ValueError, OverflowError):
                    record_is_authentic = False
                if not record_is_authentic:
                    overlay_records_are_authentic = False
        require(
            overlay_is_disjoint,
            "production design role overlay is invalid or contains duplicate node IDs",
        )
        require(
            overlay_records_are_authentic,
            "production design overlay contains an invalid coordinate, transform, or node identity",
        )

        if candidate_pool_is_sized and type(sobol_seed) is int:
            regenerated_pool = generate_candidate_pool(
                _PRODUCTION_CANDIDATE_COUNT,
                bounds=production_bounds,
                seed=sobol_seed,
                mode=str(campaign_config["mode"]),
            )
            pool_node_ids: set[str] = set()
            candidate_inventory_matches = True
            candidate_inventory_issue = ""
            for index, (record, expected_node) in enumerate(
                zip(candidate_pool, regenerated_pool)
            ):
                expected_node_id = expected_node.stable_id_hex(archive_hash)
                expected_q = geometry_transform(expected_node)
                expected_role = role_by_node_id.get(expected_node_id)
                try:
                    record_matches = (
                        isinstance(record, Mapping)
                        and type(record.get("candidate_index")) is int
                        and record["candidate_index"] == index
                        and record.get("role") == "candidate"
                        and int(record["node_id"]) == int(expected_node_id, 16)
                        and str(record["node_id_hex"]).lower() == expected_node_id
                        and all(
                            not isinstance(record[name], bool)
                            and float(record[name]) == expected_value
                            for name, expected_value in expected_node.as_dict().items()
                        )
                        and len(record["q"]) == 4
                        and all(
                            not isinstance(actual, bool) and float(actual) == target
                            for actual, target in zip(record["q"], expected_q)
                        )
                        and record.get("selected_role") == expected_role
                    )
                except (KeyError, TypeError, ValueError, OverflowError):
                    record_matches = False
                if not record_matches:
                    candidate_inventory_matches = False
                    candidate_inventory_issue = (
                        "production candidate pool differs from regenerated Sobol "
                        f"record {index}"
                    )
                    break
                pool_node_ids.add(expected_node_id)
            require(candidate_inventory_matches, candidate_inventory_issue)
            require(
                set(acquired_ids).issubset(pool_node_ids)
                and set(retired_final_ids).issubset(pool_node_ids)
                and {
                    str(item["node_id_hex"]).lower() for item in validation_records
                }.issubset(pool_node_ids)
                and {
                    str(item["node_id_hex"]).lower() for item in final_records
                }.issubset(pool_node_ids),
                "validation/final/adaptive nodes must originate in the Sobol pool",
            )

            if selection_record_valid:
                try:
                    forced_records = selection_record["forced_train"]
                    forced_nodes = tuple(
                        PhysicalNode(
                            float(item["energy_keV"]),
                            float(item["distance_m"]),
                            float(item["source_height_m"]),
                            float(item["detector_height_m"]),
                        )
                        for item in forced_records
                    )
                    required_defaults = default_forced_anchors(
                        production_bounds, mode=str(campaign_config["mode"])
                    )
                    forced_keys = {
                        tuple(value.hex() for value in node.as_dict().values())
                        for node in forced_nodes
                    }
                    defaults_present = all(
                        tuple(value.hex() for value in node.as_dict().values())
                        in forced_keys
                        for node in required_defaults
                    )
                    regenerated_selection = select_disjoint_maximin(
                        regenerated_pool,
                        sizes=_PRODUCTION_INITIAL_DESIGN_COUNTS,
                        forced={"train": forced_nodes},
                        q_domain=expected_q_domain,
                    ).as_dict(archive_hash)
                    expected_initial_ids = {
                        role: [item["node_id_hex"] for item in regenerated_selection[role]]
                        for role in ("train", "validation", "final")
                    }
                    recorded_initial_ids = selection_record["initial_role_node_ids"]
                    maximin_selection_matches = all(
                        list(recorded_initial_ids.get(role, ()))
                        == expected_initial_ids[role]
                        for role in ("train", "validation", "final")
                    )
                    initial_history_matches = bool(
                        set(expected_initial_ids["train"]).issubset(train_ids)
                        and expected_initial_ids["validation"]
                        == [str(item["node_id_hex"]).lower() for item in validation_records]
                        and (
                            expected_initial_ids["final"]
                            == [str(item["node_id_hex"]).lower() for item in final_records]
                            if int(campaign_state["final_generation"]) == 0
                            else set(expected_initial_ids["final"]).issubset(
                                retired_final_ids
                            )
                        )
                    )
                except (KeyError, TypeError, ValueError, OverflowError):
                    defaults_present = False
                    maximin_selection_matches = False
                    initial_history_matches = False
                require(
                    defaults_present,
                    "production forced-anchor inventory omits a required anchor",
                )
                require(
                    maximin_selection_matches,
                    "production initial roles do not reproduce the recorded maximin design",
                )
                require(
                    initial_history_matches,
                    "production initial maximin roles do not survive in campaign history",
                )

        initial_group_records = [
            record
            for record in campaign["worker_groups"]
            if int(record.get("first_block_id", 0)) == 1
        ]
        initial_group_ids = [
            str(record["node_id_hex"]).lower() for record in initial_group_records
        ]
        require(
            len(initial_group_ids) == len(set(initial_group_ids))
            and set(initial_group_ids) == set(role_by_node_id)
            and all(
                type(record.get("block_count")) is int
                and record["block_count"] == _PRODUCTION_INITIAL_BLOCKS_PER_NODE
                for record in initial_group_records
            ),
            "every production node must have exactly one four-block initial worker group",
        )

        # Re-emitting the byte-exact configuration catches extra/missing fields
        # and prevents a self-consistent but noncanonical JSON record from
        # acquiring publication status.
        production_physics = PhysicsConfig(**campaign_config["physics"])
        expected_json = canonical_fit_config_json(
            str(campaign_config["mode"]),
            str(campaign_config["material"]),
            production_physics,
        )
        require(
            canonical_json == expected_json,
            "configuration is not the byte-exact campaign configuration",
        )

        axis = configuration["axis"]
        physics = configuration["physics"]
        configured_mu_bins = int(axis["mu_bins"])
        require(type(axis.get("energy_bins")) is int and axis["energy_bins"] == 384,
                "production requires exactly 384 energy bins")
        require(
            not isinstance(axis.get("energy_cut_keV"), bool)
            and float(axis["energy_cut_keV"]) == 10.0,
            "production requires an exact 10-keV cutoff",
        )
        require(
            configured_mu_bins in (32, 64),
            "production requires 64 mu bins or an authorized 32-bin reduction",
        )
        require(
            axis.get("transform") == "registered-component-v1",
            "production registered-axis transform is wrong",
        )
        require(
            type(configuration.get("block_histories")) is int
            and configuration["block_histories"] == 262_144,
            "production requires 262144-history mergeable blocks",
        )
        for name in _PRODUCTION_REQUIRED_TRUE_SWITCHES:
            require(
                physics.get(name) is True,
                f"production physics switch {name} must be enabled",
            )
        require(
            physics.get("free_electron_compton") is False,
            "free-electron Compton is diagnostic-only and must be disabled",
        )
        for name, expected_value in _PRODUCTION_BASELINE_PHYSICS.items():
            actual_value = physics.get(name)
            if isinstance(expected_value, bool):
                matches_baseline = actual_value is expected_value
            elif isinstance(expected_value, int):
                matches_baseline = (
                    type(actual_value) is int and actual_value == expected_value
                )
            else:
                matches_baseline = (
                    not isinstance(actual_value, bool)
                    and float(actual_value) == expected_value
                )
            require(
                matches_baseline,
                f"production physics {name} must equal the validated baseline "
                f"{expected_value!r}",
            )
        scoring_radius_cm = float(configuration["scoring_radius_cm"])
        require(
            math.isfinite(scoring_radius_cm)
            and scoring_radius_cm == _PRODUCTION_SCORING_RADIUS_CM,
            "production scoring radius must equal the validated 4-cm baseline",
        )
        require(
            _hdf_text(manifest_attrs.get("estimator_flags", ""))
            == _PRODUCTION_ESTIMATOR_FLAGS,
            "archive estimator flags do not preserve next-event/Rayleigh/brems/annihilation",
        )

        decision = campaign.get("mu_resolution_decision")
        if configured_mu_bins == 32:
            try:
                _validated_mu_resolution_decision(
                    decision,
                    mode=str(campaign_config["mode"]),
                    material=str(campaign_config["material"]),
                    production_physics=production_physics,
                    verify_sources=True,
                )
            except (KeyError, TypeError, ValueError, OverflowError) as error:
                issues.append(f"32-bin production decision is invalid: {error}")
        else:
            require(
                decision is None,
                "a 64-bin production campaign must not carry a reduction decision",
            )

        axes = root["axes"]
        require(
            int(axes.attrs["total_energy_axis"]) == 0,
            "common total continuum does not use the primary registered axis",
        )
        energy_group = axes["energy"]
        require(
            set(energy_group.keys()) == {"0", "1", "2"},
            "archive does not contain exactly three component-native energy axes",
        )
        node_count = len(root["nodes/id"])
        for index, expected_name in enumerate(_PRODUCTION_AXIS_NAMES):
            group = energy_group[str(index)]
            require(
                _hdf_text(group.attrs.get("name", "")) == expected_name,
                f"energy axis {index} is not {expected_name}",
            )
            registered = group["registered_edges"]
            physical = group["physical_edges_keV"]
            require(
                tuple(registered.shape) == (385,),
                f"energy axis {expected_name} is not 384 bins",
            )
            require(
                tuple(physical.shape) == (node_count, 385),
                f"physical axis {expected_name} is not node-major 384-bin data",
            )
            if tuple(registered.shape) == (385,):
                registered_values = [float(value) for value in registered[...]]
                require(
                    all(
                        math.isclose(value, edge / 384.0, rel_tol=0.0, abs_tol=2.0e-15)
                        for edge, value in enumerate(registered_values)
                    ),
                    f"energy axis {expected_name} has noncanonical registered edges",
                )
                if index < 2:
                    require(
                        registered_values[192] == 0.5,
                        f"energy axis {expected_name} lacks the 192+192 landmark split",
                    )
            if len(physical.shape) == 2 and physical.shape[0] == node_count:
                require(
                    all(float(value) == 10.0 for value in physical[:, 0]),
                    f"physical axis {expected_name} does not start at 10 keV",
                )

        mu_edges = axes["mu_edges"]
        expected_mu_edges = configured_mu_bins + 1
        require(
            tuple(mu_edges.shape) == (expected_mu_edges,),
            "archive mu axis disagrees with the configured resolution",
        )
        if tuple(mu_edges.shape) == (expected_mu_edges,):
            require(
                all(
                    math.isclose(
                        float(value),
                        -1.0 + 2.0 * edge / configured_mu_bins,
                        rel_tol=0.0,
                        abs_tol=2.0e-15,
                    )
                    for edge, value in enumerate(mu_edges[...])
                ),
                "archive mu edges are not the canonical uniform partition",
            )

        descriptors = json.loads(_hdf_text(root["components"].attrs["descriptors_json"]))
        inventory = tuple(
            (
                item["name"],
                int(item["origin"]),
                int(item["post_creation_order"]),
                int(item["energy_axis"]),
            )
            for item in descriptors
        )
        require(
            inventory == _PRODUCTION_COMPONENTS,
            "archive component taxonomy is not the canonical origin/order inventory",
        )
        component_cells = 9 * 384 * configured_mu_bins
        total_cells = 384 * configured_mu_bins
        require(
            tuple(root["continuum/components/N"].shape) == (node_count, component_cells),
            "component continuum tensor is not component-native 384-by-mu data",
        )
        require(
            tuple(root["continuum/total/N"].shape) == (node_count, total_cells),
            "common total continuum tensor is not 384-by-mu data",
        )

        expected_lines = [
            {"name": "elastic_E0", "origin": 0, "post_creation_order": 3},
            {"name": "annihilation_511", "origin": 1, "post_creation_order": 3},
        ]
        line_definitions = json.loads(
            _hdf_text(axes.attrs["line_definitions_json"])
        )
        require(line_definitions == expected_lines, "archive exact-line taxonomy is incomplete")
        require(
            tuple(root["lines/energy_keV"].shape) == (node_count, 2)
            and tuple(root["lines/support_mask"].shape) == (node_count, 2),
            "archive does not preserve both exact line atoms",
        )
    except (KeyError, TypeError, ValueError, OverflowError, AttributeError) as error:
        issues.append(f"production configuration record is malformed: {error}")

    return {
        "schema": PRODUCTION_CONFIGURATION_SCHEMA,
        "passed": not issues,
        "energy_bins": 384,
        "energy_cut_keV": 10.0,
        "mu_bins": configured_mu_bins,
        "mu_policy": (
            "64-bin-pilot-or-production"
            if configured_mu_bins == 64
            else "identity-bound-64-to-32-decision"
            if configured_mu_bins == 32
            else "invalid"
        ),
        "required_true_physics_switches": list(_PRODUCTION_REQUIRED_TRUE_SWITCHES),
        "baseline_physics": dict(_PRODUCTION_BASELINE_PHYSICS),
        "candidate_count": candidate_count,
        "initial_design_counts": initial_design_counts,
        "initial_blocks_per_node": initial_blocks_per_node,
        "scoring_radius_cm": scoring_radius_cm,
        "scoring_radius_policy": "validated-4-cm-baseline-v1",
        "issues": issues,
    }


def require_production_configuration(
    root: Any, campaign: Mapping[str, Any]
) -> dict[str, Any]:
    """Return the production gate or reject a diagnostic-only configuration."""

    report = production_configuration_gate(root, campaign)
    if not report["passed"]:
        raise ValueError(
            "candidate is diagnostic-only; production configuration gate failed: "
            + "; ".join(report["issues"])
        )
    return report


def _strict_json(path: str | Path) -> Mapping[str, Any]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"non-finite JSON constant {value!r}")

    def unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON member {key!r}")
            result[key] = value
        return result

    try:
        value = json.loads(
            Path(path).read_text(encoding="utf-8"),
            parse_constant=reject_constant,
            object_pairs_hook=unique_object,
        )
    except (json.JSONDecodeError, UnicodeError) as error:
        raise ValueError(f"cannot parse validation report {path}") from error
    if not isinstance(value, Mapping):
        raise ValueError(f"validation report {path} is not a JSON object")
    return value


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _normal_node_id(value: Any) -> str:
    result = str(value).lower().removeprefix("0x")
    if len(result) != 16 or any(character not in "0123456789abcdef" for character in result):
        raise ValueError(f"invalid final node id {value!r}")
    return result


def _hash_frame(digest: Any, label: str, payload: bytes) -> None:
    """Add one unambiguous labelled byte frame to ``digest``."""

    encoded_label = label.encode("utf-8")
    digest.update(struct.pack("<QQ", len(encoded_label), len(payload)))
    digest.update(encoded_label)
    digest.update(payload)


def _normal_hdf_value(value: Any) -> Any:
    """Convert small HDF5 attribute values to a canonical JSON value."""

    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, bytes):
        return {"bytes_hex": value.hex()}
    if isinstance(value, tuple):
        return [_normal_hdf_value(item) for item in value]
    if isinstance(value, list):
        return [_normal_hdf_value(item) for item in value]
    if isinstance(value, Mapping):
        return {
            str(key): _normal_hdf_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (str, bool, int, float)) or value is None:
        if isinstance(value, float) and not math.isfinite(value):
            return {"float_hex": value.hex()}
        return value
    if hasattr(value, "item"):
        return _normal_hdf_value(value.item())
    raise ValueError(f"unsupported HDF5 attribute value {type(value).__name__}")


def _dataset_payload(array: Any, numpy: Any) -> tuple[bytes, bytes]:
    """Return canonical dtype metadata and bytes for one dataset slice."""

    values = numpy.asarray(array)
    kind = values.dtype.kind
    if kind in {"O", "S", "U"}:
        payload = bytearray()
        for value in values.reshape(-1):
            if isinstance(value, bytes):
                encoded = value
            else:
                encoded = str(value).encode("utf-8")
            payload.extend(struct.pack("<Q", len(encoded)))
            payload.extend(encoded)
        return b"utf8", bytes(payload)
    if kind not in {"b", "i", "u", "f", "c"}:
        raise ValueError(f"unsupported HDF5 scientific dataset dtype {values.dtype}")
    canonical_dtype = values.dtype.newbyteorder("<")
    canonical = numpy.ascontiguousarray(values.astype(canonical_dtype, copy=False))
    metadata = f"{kind}{values.dtype.itemsize}-little-endian".encode("ascii")
    return metadata, canonical.tobytes(order="C")


def _scientific_tree_sha256(
    group: Any,
    *,
    node_order: Sequence[int] | None = None,
    skip_children: frozenset[str] = frozenset(),
    skip_attributes: frozenset[str] = frozenset(),
) -> str:
    """Hash actual HDF5 attributes and datasets, independent of container bytes.

    Extendible node-major datasets are visited in stable node-ID order.  This
    preserves merge-order invariance while binding the identity to the exact
    arrays consumed by fitting and validation.
    """

    try:
        import numpy
    except ImportError as error:  # pragma: no cover - an HDF5 archive implies NumPy
        raise RuntimeError("archive identity requires NumPy") from error

    digest = hashlib.sha256()

    def visit(current: Any, relative: str) -> None:
        attributes = {
            str(name): _normal_hdf_value(value)
            for name, value in current.attrs.items()
            if str(name) not in skip_attributes
        }
        _hash_frame(digest, f"group:{relative}:attributes", _canonical(attributes))
        for name in sorted(current.keys()):
            if not relative and name in skip_children:
                continue
            value = current[name]
            path = f"{relative}/{name}" if relative else str(name)
            if hasattr(value, "keys"):
                visit(value, path)
                continue
            shape = tuple(int(size) for size in value.shape)
            is_node_major = (
                node_order is not None
                and value.maxshape is not None
                and len(shape) >= 1
                and value.maxshape[0] is None
                and shape[0] == len(node_order)
            )
            metadata = {
                "shape": list(shape),
                "node_ordered": is_node_major,
            }
            _hash_frame(digest, f"dataset:{path}:metadata", _canonical(metadata))
            if is_node_major:
                for logical_index, physical_index in enumerate(node_order or ()):
                    dtype_metadata, payload = _dataset_payload(
                        value[int(physical_index)], numpy
                    )
                    _hash_frame(
                        digest,
                        f"dataset:{path}:row:{logical_index}:dtype",
                        dtype_metadata,
                    )
                    _hash_frame(
                        digest,
                        f"dataset:{path}:row:{logical_index}:value",
                        payload,
                    )
            else:
                dtype_metadata, payload = _dataset_payload(value[...], numpy)
                _hash_frame(digest, f"dataset:{path}:dtype", dtype_metadata)
                _hash_frame(digest, f"dataset:{path}:value", payload)

    visit(group, "")
    return digest.hexdigest()


def archive_snapshot_sha256(
    root: Any,
    *,
    use_cached: bool = True,
    require_full_verification: bool = False,
    trusted_nodes: Mapping[str, str] | None = None,
    identities_out: dict[str, str] | None = None,
) -> str:
    """Hash the checked raw records and materialized scientific archive state.

    HDF5 container bytes are not a stable scientific identity: insertion order,
    free-space bookkeeping, and library versions may change them.  This identity
    therefore hashes canonical dataset values instead.  It deliberately does
    not trust the stored worker-checksum attributes alone: actual per-block
    records and the node-major arrays consumed by fitters are both covered.

    The identity does not depend on how it was verified.  A cached seal from an
    incremental verification (``trusted_nodes``, see
    verify_archive_materialization) serves every consumer except those that
    pass ``require_full_verification``, such as publication, which then
    re-verifies everything when the seal is not a full one.
    """

    schema = str(root.attrs.get("schema", ""))
    if use_cached:
        # Fit caches contain the exact fitter-facing aggregates but deliberately
        # omit embedded raw records.  Their stat-bound sidecar carries the
        # already-verified canonical archive identity through fitting reports.
        if schema == "CEEFIT-FIT-CACHE-v1":
            from .storage import fit_cache_source_snapshot

            return fit_cache_source_snapshot(root)
        if schema == "CEEFIT-HDF5-v1":
            from .storage import cached_archive_snapshot

            cached = cached_archive_snapshot(
                root, require_full_verification=require_full_verification
            )
            if cached is not None:
                return cached

    # A snapshot is a scientific blessing, not merely a content hash.  Verify
    # the fitter-facing aggregates exactly against sealed raw worker records
    # before assigning any identity.
    from .archive import verify_archive_materialization

    verify_archive_materialization(
        root, trusted_nodes=trusted_nodes, identities_out=identities_out
    )

    if schema != "CEEFIT-HDF5-v1":
        raise ValueError(f"unsupported fit archive schema {schema!r}")
    try:
        config_hash = str(root["manifest"].attrs["config_sha256"])
        layout_json = str(root.attrs["layout_json"])
        block_root = root["blocks"]
    except (KeyError, TypeError) as error:
        raise ValueError("fit archive lacks snapshot identity metadata") from error
    if not _is_sha256(config_hash):
        raise ValueError("fit archive has an invalid configuration hash")

    try:
        node_ids = [int(value) for value in root["nodes/id"][...]]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("fit archive lacks canonical node identities") from error
    if len(set(node_ids)) != len(node_ids):
        raise ValueError("fit archive contains duplicate node identities")
    node_order = sorted(range(len(node_ids)), key=node_ids.__getitem__)

    from .archive import node_block_inventory

    inventory: list[dict[str, Any]] = []
    raw_node_ids: list[int] = []
    for node_name in sorted(block_root, key=lambda value: int(value, 16)):
        node_id = _normal_node_id(node_name)
        raw_node_ids.append(int(node_id, 16))
        # Raw records and fold-ledger rows alike: a folded block keeps the
        # checksums its raw record would carry, so folding does not change the
        # snapshot.
        for block_id, entry in node_block_inventory(root, node_name).items():
            if block_id <= 0:
                raise ValueError("fit archive contains a non-positive block id")
            checksum = entry["payload_sha256"]
            if not _is_sha256(checksum):
                raise ValueError("fit archive block lacks a valid payload checksum")
            scientific_record = entry["scientific_record_sha256"]
            if not _is_sha256(scientific_record):
                raise ValueError("fit archive block lacks a valid scientific-record seal")
            inventory.append(
                {
                    "node_id": node_id,
                    "block_id": block_id,
                    "payload_sha256": checksum,
                    # verify_archive_materialization recomputed a raw
                    # record's seal from its datasets; a folded block's seal
                    # is protected by its fold generation's state_sha256.
                    "scientific_record_sha256": scientific_record,
                }
            )
    if len(set(raw_node_ids)) != len(raw_node_ids):
        raise ValueError("fit archive contains duplicate raw-block node identities")
    if set(raw_node_ids) != set(node_ids):
        raise ValueError("fit archive raw-block and materialized node inventories differ")
    if not inventory:
        raise ValueError("fit archive contains no raw blocks")
    aggregate_sha256 = _scientific_tree_sha256(
        root,
        node_order=node_order,
        skip_children=frozenset({"blocks", "folds"}),
    )
    identity = {
        "schema": "CEEFIT-ARCHIVE-SNAPSHOT-v1",
        "archive_schema": schema,
        "config_sha256": config_hash,
        "layout_json": layout_json,
        "materialized_science_sha256": aggregate_sha256,
        "blocks": inventory,
    }
    return hashlib.sha256(_canonical(identity)).hexdigest()


def _strict_nonnegative_integer(value: Any, name: str, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"GEANT4 anchor report {name} must be an integer >= {minimum}")
    return value


def _strict_finite(
    value: Any,
    name: str,
    *,
    minimum: float | None = None,
) -> float:
    if isinstance(value, bool):
        raise ValueError(f"GEANT4 anchor report {name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"GEANT4 anchor report {name} must be a finite number") from error
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        qualifier = f" >= {minimum:g}" if minimum is not None else ""
        raise ValueError(f"GEANT4 anchor report {name} must be finite{qualifier}")
    return result


def _require_computed_number(actual: Any, expected: float, name: str) -> float:
    """Require one serialized derived value to agree with its raw inputs.

    JSON's binary64 round trip is exact.  A few ulps are nevertheless allowed
    so a report generated on a platform with a different libm implementation
    is not rejected merely because ``hypot`` or ``sqrt`` rounded differently.
    """

    result = _strict_finite(actual, name)
    tolerance = 16.0 * math.ulp(max(abs(result), abs(expected)))
    if abs(result - expected) > tolerance:
        raise ValueError(
            f"GEANT4 anchor report {name} is inconsistent with its detailed records"
        )
    return result


def _strict_report_path(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"GEANT4 anchor report {name} must be a relative path")
    path = Path(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"GEANT4 anchor report {name} must be a contained relative path")
    return value


def _validate_geant4_anchor_report(record: Mapping[str, Any]) -> bool:
    """Recompute a GEANT4 report from its detailed, source-hashed records.

    The exact source report is retained and byte-hashed by the surrounding
    publication contract.  This function additionally prevents a handcrafted
    ``anchor_count``/``passed`` summary from becoming publication evidence: all
    arm merges, residuals, uncertainty combinations, per-anchor decisions, and
    top-level counts/maxima must close to the detailed records.
    """

    if record.get("physics_list") != _GEANT4_REQUIRED_PHYSICS_LIST:
        raise ValueError(
            "GEANT4 anchor report must use " + _GEANT4_REQUIRED_PHYSICS_LIST
        )
    comparison_sigma = _strict_finite(
        record.get("comparison_sigma"), "comparison_sigma", minimum=0.0
    )
    if comparison_sigma != _GEANT4_COMPARISON_SIGMA:
        raise ValueError(
            f"GEANT4 anchor report comparison_sigma must equal "
            f"{_GEANT4_COMPARISON_SIGMA:g}"
        )
    minimum_world_extent = _strict_finite(
        record.get("minimum_field_world_half_extent_m"),
        "minimum_field_world_half_extent_m",
        minimum=0.0,
    )
    if minimum_world_extent != _GEANT4_MINIMUM_FIELD_WORLD_HALF_EXTENT_M:
        raise ValueError(
            "GEANT4 anchor report has the wrong minimum field-world extent"
        )
    if not _is_sha256(record.get("manifest_sha256")):
        raise ValueError("GEANT4 anchor report lacks a valid manifest SHA-256")

    anchors = record.get("anchors")
    if not isinstance(anchors, list) or not anchors:
        raise ValueError("GEANT4 anchor report needs nonempty detailed anchors")
    anchor_count = _strict_nonnegative_integer(
        record.get("anchor_count"), "anchor_count", minimum=1
    )
    if anchor_count != len(anchors):
        raise ValueError("GEANT4 anchor report anchor_count differs from its anchors")

    seen_anchor_ids: set[str] = set()
    seen_node_ids: set[str] = set()
    seen_reference_paths: set[str] = set()
    seen_reference_hashes: set[str] = set()
    expected_file_count = 0
    expected_event_count = 0
    standardized_magnitudes: list[float] = []
    relative_differences: list[float] = []
    computed_anchor_passes: list[bool] = []

    for anchor_index, anchor in enumerate(anchors):
        prefix = f"anchors[{anchor_index}]"
        if not isinstance(anchor, Mapping):
            raise ValueError(f"GEANT4 anchor report {prefix} must be an object")
        anchor_id = anchor.get("anchor_id")
        if not isinstance(anchor_id, str) or not anchor_id or anchor_id in seen_anchor_ids:
            raise ValueError("GEANT4 anchor report anchor_id values must be nonempty and unique")
        seen_anchor_ids.add(anchor_id)
        try:
            node_id = _normal_node_id(anchor.get("node_id"))
        except ValueError as error:
            raise ValueError(f"GEANT4 anchor report {prefix} has an invalid node_id") from error
        if node_id in seen_node_ids:
            raise ValueError("GEANT4 anchor report repeats an archive node")
        seen_node_ids.add(node_id)

        coordinates = anchor.get("coordinates")
        if not isinstance(coordinates, list) or len(coordinates) != 5:
            raise ValueError(f"GEANT4 anchor report {prefix}.coordinates must have length five")
        coordinate_values = [
            _strict_finite(value, f"{prefix}.coordinates[{index}]")
            for index, value in enumerate(coordinates)
        ]
        if (
            coordinate_values[0] <= 0.0
            or coordinate_values[1] <= 0.0
            or coordinate_values[1]
            <= abs(coordinate_values[2] - coordinate_values[3])
        ):
            raise ValueError(f"GEANT4 anchor report {prefix} has invalid physical coordinates")
        if any(value < 0.0 for value in coordinate_values[2:]):
            raise ValueError(f"GEANT4 anchor report {prefix} has invalid physical coordinates")

        role = anchor.get("archive_role")
        if role not in {"pilot", "train", "validation", "final", "diagnostic"}:
            raise ValueError(f"GEANT4 anchor report {prefix} has an invalid archive role")
        status = anchor.get("archive_status")
        if status not in {"complete", "converged", "capped", "failed"}:
            raise ValueError(f"GEANT4 anchor report {prefix} has an invalid archive status")
        world_extent = _strict_finite(
            anchor.get("world_half_extent_m"),
            f"{prefix}.world_half_extent_m",
            minimum=0.0,
        )
        if world_extent < minimum_world_extent:
            raise ValueError(
                f"GEANT4 anchor report {prefix} has an inadequate field-world extent"
            )

        reference = anchor.get("reference")
        if not isinstance(reference, Mapping):
            raise ValueError(f"GEANT4 anchor report {prefix}.reference must be an object")
        files = reference.get("files")
        if not isinstance(files, list) or len(files) < 2:
            raise ValueError(
                f"GEANT4 anchor report {prefix} needs at least two source-hashed "
                "reference arms"
            )
        arm_records: list[dict[str, Any]] = []
        for file_index, source in enumerate(files):
            source_prefix = f"{prefix}.reference.files[{file_index}]"
            if not isinstance(source, Mapping):
                raise ValueError(f"GEANT4 anchor report {source_prefix} must be an object")
            source_path = _strict_report_path(source.get("path"), f"{source_prefix}.path")
            if source_path in seen_reference_paths:
                raise ValueError("GEANT4 anchor report assigns one reference arm more than once")
            seen_reference_paths.add(source_path)
            source_sha256 = source.get("sha256")
            if not _is_sha256(source_sha256):
                raise ValueError(f"GEANT4 anchor report {source_prefix} lacks a source SHA-256")
            if source_sha256 in seen_reference_hashes:
                raise ValueError(
                    "GEANT4 anchor report counts identical source content more than once"
                )
            seen_reference_hashes.add(source_sha256)
            events = _strict_nonnegative_integer(
                source.get("events"), f"{source_prefix}.events", minimum=2
            )
            cell_count = _strict_nonnegative_integer(
                source.get("cell_count"), f"{source_prefix}.cell_count", minimum=1
            )
            energy_bins = _strict_nonnegative_integer(
                source.get("energy_bins"), f"{source_prefix}.energy_bins", minimum=1
            )
            mu_bins = _strict_nonnegative_integer(
                source.get("mu_bins"), f"{source_prefix}.mu_bins", minimum=1
            )
            if (
                energy_bins != _GEANT4_REFERENCE_ENERGY_BINS
                or mu_bins != _GEANT4_REFERENCE_MU_BINS
                or cell_count != energy_bins * mu_bins
            ):
                raise ValueError(
                    f"GEANT4 anchor report {source_prefix} has the wrong scorer layout"
                )
            source_e0 = _strict_finite(
                source.get("e0_keV"), f"{source_prefix}.e0_keV", minimum=0.0
            )
            source_radius = _strict_finite(
                source.get("sphere_radius_cm"),
                f"{source_prefix}.sphere_radius_cm",
                minimum=0.0,
            )
            if source_e0 <= 0.0 or source_radius <= 0.0:
                raise ValueError(
                    f"GEANT4 anchor report {source_prefix} has invalid scorer metadata"
                )
            gdml_radius = _strict_finite(
                source.get("gdml_scoring_radius_cm"),
                f"{source_prefix}.gdml_scoring_radius_cm",
                minimum=0.0,
            )
            if not math.isclose(
                gdml_radius, source_radius, rel_tol=1.0e-12, abs_tol=1.0e-12
            ):
                raise ValueError(
                    f"GEANT4 anchor report {source_prefix} scorer/GDML radii differ"
                )
            energy_min = _strict_finite(
                source.get("energy_min_keV"),
                f"{source_prefix}.energy_min_keV",
                minimum=0.0,
            )
            energy_max = _strict_finite(
                source.get("energy_max_keV"),
                f"{source_prefix}.energy_max_keV",
                minimum=0.0,
            )
            mu_min = _strict_finite(source.get("mu_min"), f"{source_prefix}.mu_min")
            mu_max = _strict_finite(source.get("mu_max"), f"{source_prefix}.mu_max")
            if energy_min <= 0.0 or energy_max < source_e0 or energy_max <= energy_min:
                raise ValueError(f"GEANT4 anchor report {source_prefix} has invalid energy support")
            if mu_min != -1.0 or mu_max != 1.0:
                raise ValueError(f"GEANT4 anchor report {source_prefix} has invalid mu support")
            if source.get("observable") != _GEANT4_REQUIRED_OBSERVABLE:
                raise ValueError(f"GEANT4 anchor report {source_prefix} has wrong observable")
            if source.get("physics_list") != _GEANT4_REQUIRED_PHYSICS_LIST:
                raise ValueError(f"GEANT4 anchor report {source_prefix} has wrong physics list")
            source_world = _strict_finite(
                source.get("world_half_extent_m"),
                f"{source_prefix}.world_half_extent_m",
                minimum=minimum_world_extent,
            )
            if source.get("primary_sampling") != _GEANT4_REQUIRED_SAMPLING:
                raise ValueError(f"GEANT4 anchor report {source_prefix} uses biased sampling")
            rng_engine = source.get("rng_engine")
            if not isinstance(rng_engine, str) or not rng_engine:
                raise ValueError(f"GEANT4 anchor report {source_prefix} lacks RNG engine")
            rng_seed = _strict_nonnegative_integer(
                source.get("rng_seed"), f"{source_prefix}.rng_seed", minimum=1
            )
            production_cut = _strict_finite(
                source.get("production_cut_mm"),
                f"{source_prefix}.production_cut_mm",
                minimum=0.0,
            )
            if not 0.0 < production_cut <= 0.7:
                raise ValueError(f"GEANT4 anchor report {source_prefix} has invalid production cut")
            geometry_tier = source.get("geometry_tier")
            material = source.get("material")
            if geometry_tier not in {"air", "field"} or (
                geometry_tier == "air" and material != "none"
            ) or (
                geometry_tier == "field" and material not in {"soil", "concrete"}
            ):
                raise ValueError(f"GEANT4 anchor report {source_prefix} has invalid geometry identity")
            distance_cm = _strict_finite(
                source.get("distance_cm"), f"{source_prefix}.distance_cm", minimum=0.0
            )
            source_height_cm = _strict_finite(
                source.get("source_height_cm"),
                f"{source_prefix}.source_height_cm",
                minimum=0.0,
            )
            detector_height_cm = _strict_finite(
                source.get("detector_height_cm"),
                f"{source_prefix}.detector_height_cm",
                minimum=0.0,
            )
            if distance_cm <= abs(source_height_cm - detector_height_cm):
                raise ValueError(f"GEANT4 anchor report {source_prefix} has invalid coordinates")
            gdml_sha256 = source.get("gdml_sha256")
            macro_sha256 = source.get("macro_sha256")
            if not _is_sha256(gdml_sha256) or not _is_sha256(macro_sha256):
                raise ValueError(f"GEANT4 anchor report {source_prefix} lacks input hashes")
            inputs = source.get("inputs")
            if not isinstance(inputs, Mapping) or set(inputs) != {"gdml", "macro"}:
                raise ValueError(f"GEANT4 anchor report {source_prefix} lacks bound inputs")
            for input_name, wanted_hash in (
                ("gdml", gdml_sha256),
                ("macro", macro_sha256),
            ):
                input_record = inputs[input_name]
                if not isinstance(input_record, Mapping):
                    raise ValueError(f"GEANT4 anchor report {source_prefix} input is invalid")
                _strict_report_path(
                    input_record.get("path"), f"{source_prefix}.inputs.{input_name}.path"
                )
                if input_record.get("sha256") != wanted_hash:
                    raise ValueError(f"GEANT4 anchor report {source_prefix} input hash does not close")

            if source.get("scattered_total_moment_schema") != _GEANT4_TOTAL_MOMENT_SCHEMA:
                raise ValueError(f"GEANT4 anchor report {source_prefix} lacks exact total moments")
            moment_n = _strict_nonnegative_integer(
                source.get("scattered_total_N"),
                f"{source_prefix}.scattered_total_N",
                minimum=2,
            )
            if moment_n != events:
                raise ValueError(f"GEANT4 anchor report {source_prefix} moment count does not close")
            moment_s1 = _strict_finite(
                source.get("scattered_total_S1_cm"),
                f"{source_prefix}.scattered_total_S1_cm",
                minimum=0.0,
            )
            moment_s2 = _strict_finite(
                source.get("scattered_total_S2_cm2"),
                f"{source_prefix}.scattered_total_S2_cm2",
                minimum=0.0,
            )
            centered = moment_s2 - moment_s1 * moment_s1 / events
            centered_scale = max(abs(moment_s2), abs(moment_s1 * moment_s1 / events), 1.0)
            if centered < -64.0 * math.ulp(centered_scale):
                raise ValueError(f"GEANT4 anchor report {source_prefix} moments have negative variance")
            volume = (4.0 / 3.0) * math.pi * source_radius**3
            arm_mean = moment_s1 / (events * volume)
            arm_se = math.sqrt(max(0.0, centered) / (events * (events - 1))) / volume
            _require_computed_number(
                source.get("scattered_mean"), arm_mean, f"{source_prefix}.scattered_mean"
            )
            _require_computed_number(
                source.get("scattered_standard_error"),
                arm_se,
                f"{source_prefix}.scattered_standard_error",
            )
            arm_records.append(
                {
                    "events": events,
                    "mean": arm_mean,
                    "standard_error": arm_se,
                    "S1": moment_s1,
                    "S2": moment_s2,
                    "e0_keV": source_e0,
                    "radius_cm": source_radius,
                    "energy_min_keV": energy_min,
                    "energy_max_keV": energy_max,
                    "world_half_extent_m": source_world,
                    "rng_engine": rng_engine,
                    "rng_seed": rng_seed,
                    "production_cut_mm": production_cut,
                    "geometry_tier": geometry_tier,
                    "material": material,
                    "distance_cm": distance_cm,
                    "source_height_cm": source_height_cm,
                    "detector_height_cm": detector_height_cm,
                }
            )

        cell_counts = {
            int(source["cell_count"])
            for source in files
        }
        if len(cell_counts) != 1:
            raise ValueError(
                f"GEANT4 anchor report {prefix} reference arms have different cell layouts"
            )
        invariant_fields = (
            "e0_keV", "radius_cm", "energy_min_keV", "energy_max_keV",
            "world_half_extent_m", "rng_engine", "production_cut_mm",
            "geometry_tier", "material", "distance_cm", "source_height_cm",
            "detector_height_cm",
        )
        if any(len({arm[field] for arm in arm_records}) != 1 for field in invariant_fields):
            raise ValueError(
                f"GEANT4 anchor report {prefix} reference arms have different scorer metadata"
            )
        if len({arm["rng_seed"] for arm in arm_records}) != len(arm_records):
            raise ValueError(f"GEANT4 anchor report {prefix} reuses an RNG seed")
        first_arm = arm_records[0]
        if not math.isclose(
            first_arm["e0_keV"], coordinate_values[0], rel_tol=1.0e-8, abs_tol=1.0e-12
        ):
            raise ValueError(
                f"GEANT4 anchor report {prefix} E0 does not match its archive coordinates"
            )
        for coordinate_index, field in enumerate(
            ("e0_keV", "distance_cm", "source_height_cm", "detector_height_cm")
        ):
            if not math.isclose(
                float(first_arm[field]), coordinate_values[coordinate_index],
                rel_tol=1.0e-9, abs_tol=1.0e-9,
            ):
                raise ValueError(f"GEANT4 anchor report {prefix} coordinates do not close")
        if not math.isclose(world_extent, float(first_arm["world_half_extent_m"])):
            raise ValueError(f"GEANT4 anchor report {prefix} world extent does not close")
        total_events = sum(int(value["events"]) for value in arm_records)
        total_s1 = math.fsum(float(value["S1"]) for value in arm_records)
        total_s2 = math.fsum(float(value["S2"]) for value in arm_records)
        volume = (4.0 / 3.0) * math.pi * float(first_arm["radius_cm"]) ** 3
        reference_mean = total_s1 / (total_events * volume)
        pooled_centered = total_s2 - total_s1 * total_s1 / total_events
        pooled_scale = max(abs(total_s2), abs(total_s1 * total_s1 / total_events), 1.0)
        if pooled_centered < -64.0 * math.ulp(pooled_scale):
            raise ValueError(f"GEANT4 anchor report {prefix} merged moments have negative variance")
        per_history_se = math.sqrt(
            max(0.0, pooled_centered) / (total_events * (total_events - 1))
        ) / volume
        weights = [float(value["events"]) for value in arm_records]
        means = [float(value["mean"]) for value in arm_records]
        weight_sum = math.fsum(weights)
        squared_weight_sum = math.fsum(weight * weight for weight in weights)
        between_denominator = weight_sum - squared_weight_sum / weight_sum
        effective_arms = weight_sum * weight_sum / squared_weight_sum
        if between_denominator <= 0.0 or effective_arms <= 1.0:
            raise ValueError(
                f"GEANT4 anchor report {prefix} has insufficient independent arm weight"
            )
        between_variance = math.fsum(
            weight * (mean - reference_mean) ** 2
            for weight, mean in zip(weights, means)
        ) / between_denominator
        between_arm_se = math.sqrt(max(0.0, between_variance) / effective_arms)
        reference_se = max(per_history_se, between_arm_se)
        if _strict_nonnegative_integer(
            reference.get("events"), f"{prefix}.reference.events", minimum=1
        ) != total_events:
            raise ValueError(f"GEANT4 anchor report {prefix}.reference event count does not close")
        if _strict_nonnegative_integer(
            reference.get("arm_count"), f"{prefix}.reference.arm_count", minimum=2
        ) != len(arm_records):
            raise ValueError(f"GEANT4 anchor report {prefix}.reference arm count does not close")
        if _strict_nonnegative_integer(
            reference.get("energy_bins"), f"{prefix}.reference.energy_bins", minimum=1
        ) != _GEANT4_REFERENCE_ENERGY_BINS or _strict_nonnegative_integer(
            reference.get("mu_bins"), f"{prefix}.reference.mu_bins", minimum=1
        ) != _GEANT4_REFERENCE_MU_BINS:
            raise ValueError(f"GEANT4 anchor report {prefix}.reference layout does not close")
        _require_computed_number(
            reference.get("e0_keV"),
            float(first_arm["e0_keV"]),
            f"{prefix}.reference.e0_keV",
        )
        _require_computed_number(
            reference.get("sphere_radius_cm"),
            float(first_arm["radius_cm"]),
            f"{prefix}.reference.sphere_radius_cm",
        )
        _require_computed_number(
            reference.get("gdml_scoring_radius_cm"),
            float(first_arm["radius_cm"]),
            f"{prefix}.reference.gdml_scoring_radius_cm",
        )
        for field in (
            "energy_min_keV",
            "energy_max_keV",
            "world_half_extent_m",
            "production_cut_mm",
            "distance_cm",
            "source_height_cm",
            "detector_height_cm",
        ):
            _require_computed_number(
                reference.get(field),
                float(first_arm[field]),
                f"{prefix}.reference.{field}",
            )
        _require_computed_number(
            reference.get("mu_min"), -1.0, f"{prefix}.reference.mu_min"
        )
        _require_computed_number(
            reference.get("mu_max"), 1.0, f"{prefix}.reference.mu_max"
        )
        for field, expected in (
            ("observable", _GEANT4_REQUIRED_OBSERVABLE),
            ("physics_list", _GEANT4_REQUIRED_PHYSICS_LIST),
            ("primary_sampling", _GEANT4_REQUIRED_SAMPLING),
            ("rng_engine", first_arm["rng_engine"]),
            ("geometry_tier", first_arm["geometry_tier"]),
            ("material", first_arm["material"]),
            ("scattered_total_moment_schema", _GEANT4_TOTAL_MOMENT_SCHEMA),
        ):
            if reference.get(field) != expected:
                raise ValueError(
                    f"GEANT4 anchor report {prefix}.reference.{field} does not close"
                )
        if _strict_nonnegative_integer(
            reference.get("scattered_total_N"),
            f"{prefix}.reference.scattered_total_N",
            minimum=2,
        ) != total_events:
            raise ValueError(
                f"GEANT4 anchor report {prefix} total moment count does not close"
            )
        _require_computed_number(
            reference.get("scattered_total_S1_cm"),
            total_s1,
            f"{prefix}.reference.scattered_total_S1_cm",
        )
        _require_computed_number(
            reference.get("scattered_total_S2_cm2"),
            total_s2,
            f"{prefix}.reference.scattered_total_S2_cm2",
        )
        _require_computed_number(
            reference.get("mean"), reference_mean, f"{prefix}.reference.mean"
        )
        _require_computed_number(
            reference.get("per_history_standard_error"),
            per_history_se,
            f"{prefix}.reference.per_history_standard_error",
        )
        _require_computed_number(
            reference.get("between_arm_standard_error"),
            between_arm_se,
            f"{prefix}.reference.between_arm_standard_error",
        )
        _require_computed_number(
            reference.get("standard_error"),
            reference_se,
            f"{prefix}.reference.standard_error",
        )
        expected_file_count += len(arm_records)
        expected_event_count += total_events

        archive = anchor.get("archive")
        if not isinstance(archive, Mapping):
            raise ValueError(f"GEANT4 anchor report {prefix}.archive must be an object")
        _strict_nonnegative_integer(
            archive.get("histories"), f"{prefix}.archive.histories", minimum=2
        )
        archive_mean = _strict_finite(
            archive.get("mean"), f"{prefix}.archive.mean", minimum=0.0
        )
        archive_se = _strict_finite(
            archive.get("standard_error"),
            f"{prefix}.archive.standard_error",
            minimum=0.0,
        )
        difference = archive_mean - reference_mean
        combined_se = math.hypot(archive_se, reference_se)
        if combined_se == 0.0:
            standardized = (
                0.0 if difference == 0.0 else math.copysign(1.0e100, difference)
            )
        else:
            standardized = difference / combined_se
        relative_difference = (
            abs(difference) / abs(reference_mean)
            if reference_mean != 0.0
            else (0.0 if difference == 0.0 else 1.0e100)
        )
        _require_computed_number(anchor.get("difference"), difference, f"{prefix}.difference")
        _require_computed_number(
            anchor.get("combined_standard_error"),
            combined_se,
            f"{prefix}.combined_standard_error",
        )
        _require_computed_number(
            anchor.get("standardized_difference"),
            standardized,
            f"{prefix}.standardized_difference",
        )
        _require_computed_number(
            anchor.get("absolute_relative_difference"),
            relative_difference,
            f"{prefix}.absolute_relative_difference",
        )
        expected_pass = (
            status in {"complete", "converged"}
            and abs(standardized) <= comparison_sigma
        )
        if type(anchor.get("passed")) is not bool or anchor["passed"] != expected_pass:
            raise ValueError(f"GEANT4 anchor report {prefix}.passed is inconsistent")
        computed_anchor_passes.append(expected_pass)
        standardized_magnitudes.append(abs(standardized))
        relative_differences.append(relative_difference)

    if _strict_nonnegative_integer(
        record.get("reference_file_count"), "reference_file_count", minimum=1
    ) != expected_file_count:
        raise ValueError("GEANT4 anchor report reference_file_count does not close")
    if _strict_nonnegative_integer(
        record.get("reference_event_count"), "reference_event_count", minimum=1
    ) != expected_event_count:
        raise ValueError("GEANT4 anchor report reference_event_count does not close")
    _require_computed_number(
        record.get("max_abs_standardized_difference"),
        max(standardized_magnitudes),
        "max_abs_standardized_difference",
    )
    _require_computed_number(
        record.get("max_absolute_relative_difference"),
        max(relative_differences),
        "max_absolute_relative_difference",
    )
    computed_passed = all(computed_anchor_passes)
    if (
        type(record.get("within_combined_uncertainty")) is not bool
        or record["within_combined_uncertainty"] != computed_passed
    ):
        raise ValueError(
            "GEANT4 anchor report within_combined_uncertainty is inconsistent"
        )
    if type(record.get("passed")) is not bool or record["passed"] != computed_passed:
        raise ValueError("GEANT4 anchor report passed flag is inconsistent")
    return computed_passed


def _validate_uncertainty_calibration_report(record: Mapping[str, Any]) -> bool:
    """Recompute the calibration summary from its detailed family records."""

    def explicit_bool(value: Any, name: str) -> bool:
        if type(value) is not bool:
            raise ValueError(f"uncertainty calibration {name} must be boolean")
        return bool(value)

    repeated = _strict_nonnegative_integer(
        record.get("repeated_null_runs"), "repeated_null_runs", minimum=0
    )
    failures = record.get("failures", [])
    if not isinstance(failures, list) or any(
        not isinstance(failure, Mapping) for failure in failures
    ):
        raise ValueError("uncertainty calibration failures must be a list of records")
    criteria = record.get("criteria")
    if not isinstance(criteria, Mapping):
        raise ValueError("uncertainty calibration lacks its criteria record")
    if _strict_nonnegative_integer(
        criteria.get("minimum_independent_comparisons"),
        "criteria.minimum_independent_comparisons",
        minimum=1,
    ) != _CALIBRATION_MINIMUM_COMPARISONS:
        raise ValueError("uncertainty calibration uses the wrong comparison minimum")
    if criteria.get("required_moment_and_replica_tensors") != list(
        _CALIBRATION_MOMENT_FAMILIES
    ):
        raise ValueError("uncertainty calibration uses the wrong moment-family inventory")

    replication = record.get("replication")
    if not isinstance(replication, Mapping):
        raise ValueError("uncertainty calibration lacks replication diagnostics")
    if _strict_nonnegative_integer(
        replication.get("comparison_count"),
        "replication.comparison_count",
        minimum=0,
    ) != repeated:
        raise ValueError("uncertainty calibration comparison count does not close")
    if _strict_nonnegative_integer(
        replication.get("blocks_used"), "replication.blocks_used", minimum=0
    ) != 2 * repeated:
        raise ValueError("uncertainty calibration block count does not close")
    sufficient = repeated >= _CALIBRATION_MINIMUM_COMPARISONS
    if explicit_bool(replication.get("sufficient"), "replication.sufficient") != sufficient:
        raise ValueError("uncertainty calibration replication sufficiency is inconsistent")
    independent = explicit_bool(
        replication.get("independent"), "replication.independent"
    )

    raw_total = record.get("raw_total_s2")
    cell_s2 = record.get("cell_s2")
    stopping = record.get("stopping_outer_moments")
    replica = record.get("replica_covariance")
    if not all(isinstance(value, Mapping) for value in (raw_total, cell_s2, stopping, replica)):
        raise ValueError("uncertainty calibration lacks detailed estimator families")

    raw_total_passed = explicit_bool(raw_total.get("passed"), "raw_total_s2.passed")
    cell_passed = explicit_bool(cell_s2.get("passed"), "cell_s2.passed")
    cell_inventory = explicit_bool(
        cell_s2.get("family_inventory_complete"),
        "cell_s2.family_inventory_complete",
    )
    if cell_s2.get("required_tensors") != list(_CALIBRATION_MOMENT_FAMILIES):
        raise ValueError("uncertainty calibration cell-S2 inventory differs")
    cell_tensors = cell_s2.get("tensors")
    if not isinstance(cell_tensors, Mapping) or set(cell_tensors) != set(
        _CALIBRATION_MOMENT_FAMILIES
    ):
        raise ValueError("uncertainty calibration cell-S2 tensors are incomplete")
    tensor_passes: list[bool] = []
    segment_passes: list[bool] = []
    for family in _CALIBRATION_MOMENT_FAMILIES:
        tensor = cell_tensors[family]
        if not isinstance(tensor, Mapping):
            raise ValueError(f"uncertainty calibration cell-S2 family {family!r} is malformed")
        tensor_passes.append(
            explicit_bool(tensor.get("passed"), f"cell_s2.tensors.{family}.passed")
        )
        segments = tensor.get("segments")
        if not isinstance(segments, list) or not segments:
            raise ValueError(f"uncertainty calibration cell-S2 family {family!r} lacks segments")
        for index, segment in enumerate(segments):
            if not isinstance(segment, Mapping):
                raise ValueError("uncertainty calibration cell-S2 segment is malformed")
            segment_passes.append(
                explicit_bool(
                    segment.get("passed"),
                    f"cell_s2.tensors.{family}.segments[{index}].passed",
                )
            )
    if cell_passed != bool(cell_inventory and all(tensor_passes)):
        raise ValueError("uncertainty calibration cell-S2 passed flag is inconsistent")

    stopping_passed = explicit_bool(
        stopping.get("passed"), "stopping_outer_moments.passed"
    )
    dimensions = stopping.get("dimensions")
    if not isinstance(dimensions, list) or len(dimensions) != 16:
        raise ValueError("uncertainty calibration needs all 16 stopping dimensions")
    dimension_indices: set[int] = set()
    dimension_passes: list[bool] = []
    for position, dimension in enumerate(dimensions):
        if not isinstance(dimension, Mapping):
            raise ValueError("uncertainty calibration stopping dimension is malformed")
        index = _strict_nonnegative_integer(
            dimension.get("index"),
            f"stopping_outer_moments.dimensions[{position}].index",
            minimum=0,
        )
        dimension_indices.add(index)
        dimension_passes.append(
            explicit_bool(
                dimension.get("passed"),
                f"stopping_outer_moments.dimensions[{position}].passed",
            )
        )
    if dimension_indices != set(range(16)):
        raise ValueError("uncertainty calibration stopping indices are incomplete")

    replica_passed = explicit_bool(replica.get("passed"), "replica_covariance.passed")
    replica_inventory = explicit_bool(
        replica.get("family_inventory_complete"),
        "replica_covariance.family_inventory_complete",
    )
    families = replica.get("families")
    if not isinstance(families, Mapping) or set(families) != set(
        _CALIBRATION_MOMENT_FAMILIES
    ):
        raise ValueError("uncertainty calibration replica families are incomplete")
    replica_family_passes: list[bool] = []
    for family in _CALIBRATION_MOMENT_FAMILIES:
        value = families[family]
        if not isinstance(value, Mapping):
            raise ValueError(f"uncertainty calibration replica family {family!r} is malformed")
        replica_family_passes.append(
            explicit_bool(
                value.get("passed"),
                f"replica_covariance.families.{family}.passed",
            )
        )
    if replica_passed != bool(replica_inventory and all(replica_family_passes)):
        raise ValueError(
            "uncertainty calibration replica-covariance passed flag is inconsistent"
        )

    expected_s2 = raw_total_passed and cell_passed
    expected_family_coverage = bool(
        cell_inventory
        and all(segment_passes)
        and replica_inventory
        and all(dimension_passes)
    )
    summary_checks = {
        "s2_passed": expected_s2,
        "coarse_covariance_passed": stopping_passed,
        "replica_covariance_passed": replica_passed,
        "family_coverage_passed": expected_family_coverage,
    }
    for name, expected in summary_checks.items():
        if explicit_bool(record.get(name), name) != expected:
            raise ValueError(f"uncertainty calibration {name} is inconsistent")
    computed_passed = bool(
        sufficient
        and independent
        and expected_s2
        and stopping_passed
        and replica_passed
        and expected_family_coverage
        and not failures
    )
    reported_passed = explicit_bool(record.get("passed"), "passed")
    if reported_passed and failures:
        raise ValueError(
            "passing uncertainty calibration cannot contain failure records"
        )
    if reported_passed != computed_passed:
        raise ValueError("uncertainty calibration passed flag is inconsistent")
    return computed_passed


def _source_report(
    path: Path,
    name: str,
    expected_bindings: Mapping[str, str],
) -> Mapping[str, Any]:
    contract = _SOURCE_REPORT_CONTRACTS[name]
    record = _strict_json(path)
    if record.get("schema") != contract["schema"]:
        raise ValueError(f"publication source report {name!r} has the wrong schema")
    if type(record.get("passed")) is not bool:
        raise ValueError(f"publication source report {name!r} needs an explicit passed flag")
    for binding in contract["bindings"]:
        if record.get(binding) != expected_bindings[binding]:
            raise ValueError(
                f"publication source report {name!r} {binding} does not match the candidate"
            )
    for field in contract["fields"]:
        if field not in record:
            raise ValueError(f"publication source report {name!r} lacks {field!r}")
    if name == "registered_e0_gate":
        from .registered_e0 import GENERATOR_SCHEMA, validate_registered_e0_report

        validate_registered_e0_report(record)
        generator = record.get("generator")
        if not isinstance(generator, Mapping) or generator.get("schema") != GENERATOR_SCHEMA:
            raise ValueError("registered-E0 report lacks its integrated generator record")
        if not isinstance(generator.get("argv"), list) or not generator["argv"]:
            raise ValueError("registered-E0 report lacks the exact generator arguments")
        if not isinstance(generator.get("input_parameters"), Mapping) or not isinstance(
            generator.get("output_parameters"), Mapping
        ):
            raise ValueError("registered-E0 report lacks input/output parameters")
        inputs = generator["input_parameters"]
        canonical_json = inputs.get("canonical_configuration_json")
        if not isinstance(canonical_json, str):
            raise ValueError("registered-E0 report lacks canonical configuration JSON")
        try:
            canonical_record = json.loads(canonical_json)
        except json.JSONDecodeError as error:
            raise ValueError("registered-E0 canonical configuration JSON is invalid") from error
        if (
            hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()
            != record.get("config_sha256")
            or inputs.get("config_sha256") != record.get("config_sha256")
            or inputs.get("canonical_configuration") != canonical_record
            or inputs.get("gate_configuration") != record.get("gate_configuration")
        ):
            raise ValueError("registered-E0 input configuration binding does not close")
        gate = record["gate_configuration"]
        canonical_material = canonical_record.get("material")
        expected_material = (
            "Air"
            if canonical_record.get("geometry_tier") == "air"
            else canonical_material.get("name")
            if isinstance(canonical_material, Mapping)
            else None
        )
        if (
            canonical_record.get("geometry_tier") != gate.get("geometry_tier")
            or expected_material != gate.get("material")
            or canonical_record.get("scoring_radius_cm")
            != gate.get("scoring_radius_cm")
            or canonical_record.get("axis", {}).get("energy_cut_keV")
            != gate.get("energy_cut_keV")
        ):
            raise ValueError("registered-E0 gate slice differs from canonical configuration")
    elif name == "geant4_anchor_agreement":
        _validate_geant4_anchor_report(record)
        manifest_relative = record.get("manifest_path")
        if not isinstance(manifest_relative, str) or not manifest_relative:
            raise ValueError("GEANT4 anchor report lacks its source manifest path")
        manifest_path = (path.parent / manifest_relative).resolve()
        if not manifest_path.is_file():
            raise ValueError("GEANT4 anchor source manifest is missing")
        if file_sha256(manifest_path) != record.get("manifest_sha256"):
            raise ValueError("GEANT4 anchor source manifest SHA-256 differs")
        for anchor in record["anchors"]:
            for source in anchor["reference"]["files"]:
                source_path = (manifest_path.parent / source["path"]).resolve()
                if not source_path.is_file() or file_sha256(source_path) != source["sha256"]:
                    raise ValueError("GEANT4 raw reference SHA-256 differs from its file")
                for input_record in source["inputs"].values():
                    input_path = (manifest_path.parent / input_record["path"]).resolve()
                    if (
                        not input_path.is_file()
                        or file_sha256(input_path) != input_record["sha256"]
                    ):
                        raise ValueError("GEANT4 GDML/macro SHA-256 differs from its file")
    elif name == "uncertainty_calibration":
        _validate_uncertainty_calibration_report(record)
    elif name == "python_cpp_parity":
        artifact_paths = record.get("artifact_paths")
        if not isinstance(artifact_paths, Mapping) or set(artifact_paths) != {
            "candidate_payload"
        }:
            raise ValueError("Python/C++ parity report lacks its candidate payload path")
        relative = artifact_paths["candidate_payload"]
        if not isinstance(relative, str) or not relative:
            raise ValueError("Python/C++ parity candidate payload path is invalid")
        candidate_path = (path.parent / relative).resolve()
        if (
            not candidate_path.is_file()
            or file_sha256(candidate_path) != record["candidate_payload_sha256"]
        ):
            raise ValueError(
                "Python/C++ parity candidate_payload SHA-256 differs from its file"
            )
    return record


def _check_record(
    checks: Mapping[str, Any],
    name: str,
    *,
    evidence_path: Path,
    expected_bindings: Mapping[str, str],
) -> Mapping[str, Any]:
    value = checks.get(name)
    if not isinstance(value, Mapping):
        raise ValueError(f"publication evidence is missing {name!r}")
    if not _is_sha256(value.get("report_sha256")):
        raise ValueError(f"publication check {name!r} lacks a report SHA-256")
    if type(value.get("passed")) is not bool:
        raise ValueError(f"publication check {name!r} needs an explicit passed flag")
    relative = value.get("report_path")
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError(f"publication check {name!r} needs a relative source-report path")
    source_path = evidence_path.parent / relative
    if not source_path.is_file():
        raise ValueError(f"publication source report {name!r} is missing")
    if file_sha256(source_path) != value["report_sha256"]:
        raise ValueError(f"publication source report {name!r} SHA-256 differs")
    source = _source_report(source_path, name, expected_bindings)
    contract = _SOURCE_REPORT_CONTRACTS[name]
    if source["passed"] != value["passed"]:
        raise ValueError(f"publication check {name!r} passed flag differs from its source")
    for field in contract["fields"]:
        if source[field] != value.get(field):
            raise ValueError(f"publication check {name!r} field {field!r} differs from its source")
    return {
        "report_path": relative,
        "report_sha256": value["report_sha256"],
        "source_schema": source["schema"],
        "passed": source["passed"],
        **{field: source[field] for field in contract["fields"]},
    }


def _metric_passes(name: str, details: Mapping[str, Any]) -> bool:
    if name == "registered_e0_gate":
        coverage = float(details.get("shape_coverage_fraction", -1.0))
        integrated = float(details.get("max_integrated_relative_error", math.inf))
        joint = float(details.get("max_joint_tv", math.inf))
        energy = float(details.get("max_energy_marginal_tv", math.inf))
        angular = float(details.get("max_angular_marginal_tv", math.inf))
        return bool(
            details.get("predefined_index_bands") is True
            and details.get("statistically_empty_bands_included") is True
            and details.get("all_node_stopping_targets_met") is True
            and all(math.isfinite(value) for value in (coverage, integrated, joint, energy, angular))
            and 1.0 - 1.0e-12 <= coverage <= 1.0 + 1.0e-12
            and 0.0 <= integrated <= 0.01
            and 0.0 <= joint <= 0.01
            and 0.0 <= energy <= 0.005
            and 0.0 <= angular <= 0.005
        )
    if name == "uncertainty_calibration":
        criteria = details.get("criteria")
        return bool(
            int(details.get("repeated_null_runs", 0))
            >= _CALIBRATION_MINIMUM_COMPARISONS
            and details.get("s2_passed") is True
            and details.get("coarse_covariance_passed") is True
            and details.get("replica_covariance_passed") is True
            and details.get("family_coverage_passed") is True
            and isinstance(criteria, Mapping)
            and criteria.get("minimum_independent_comparisons")
            == _CALIBRATION_MINIMUM_COMPARISONS
            and criteria.get("required_moment_and_replica_tensors")
            == list(_CALIBRATION_MOMENT_FAMILIES)
        )
    if name == "geant4_anchor_agreement":
        return _validate_geant4_anchor_report(details)
    if name == "python_cpp_parity":
        parity_error = float(details.get("max_absolute_difference", math.inf))
        tolerance = float(details.get("tolerance", math.inf))
        coverage = details.get("coverage")
        return bool(
            math.isfinite(parity_error)
            and 0.0 <= parity_error <= 1.0e-9
            and math.isfinite(tolerance)
            and 0.0 < tolerance <= 1.0e-9
            and details.get("protocol") == "CEEFIT-PARITY-REQUEST-v1"
            and _is_sha256(details.get("candidate_payload_sha256"))
            and _is_sha256(details.get("request_sha256"))
            and _is_sha256(details.get("cpp_result_sha256"))
            and int(details.get("query_count", 0)) >= 11
            and int(details.get("compared_value_count", 0)) > 0
            and isinstance(coverage, Mapping)
            and int(coverage.get("interior", 0)) >= 1
            and int(coverage.get("domain_face", 0)) >= 4
            and int(coverage.get("brems_threshold", 0)) >= 3
            and int(coverage.get("pair_threshold", 0)) >= 3
        )
    raise ValueError(f"unknown publication check {name!r}")


def validate_publication_evidence(
    path: str | Path | None,
    *,
    config_sha256: str,
    archive_snapshot_sha256: str,
    model_content_sha256: str,
    final_node_ids: Sequence[str],
) -> dict[str, Any]:
    """Validate external production evidence against one candidate identity.

    The result is report-friendly and never treats absence as success.  A
    structurally malformed or identity-mismatched report raises; a complete
    report containing a failed scientific check returns ``passed=False``.
    """

    expected_ids = tuple(sorted(_normal_node_id(value) for value in final_node_ids))
    if len(expected_ids) != len(set(expected_ids)):
        raise ValueError("untouched final node IDs are not unique")
    if path is None:
        return {
            "supplied": False,
            "passed": False,
            "reason": "publication evidence was not supplied",
            "required_checks": list(_MANDATORY_CHECKS),
        }

    evidence_path = Path(path)
    record = _strict_json(evidence_path)
    if not isinstance(record, Mapping) or record.get("schema") != PUBLICATION_EVIDENCE_SCHEMA:
        raise ValueError("unsupported publication-evidence schema")
    expected_bindings = {
        "config_sha256": config_sha256,
        "archive_snapshot_sha256": archive_snapshot_sha256,
        "model_content_sha256": model_content_sha256,
    }
    for name, expected in expected_bindings.items():
        if not _is_sha256(expected) or record.get(name) != expected:
            raise ValueError(f"publication evidence {name} does not match the candidate")
    supplied_ids = tuple(
        sorted(_normal_node_id(value) for value in record.get("final_node_ids", ()))
    )
    if supplied_ids != expected_ids:
        raise ValueError("publication evidence does not cover exactly the untouched final set")

    checks = record.get("checks")
    if not isinstance(checks, Mapping):
        raise ValueError("publication evidence lacks its checks object")
    details = {
        name: dict(
            _check_record(
                checks,
                name,
                evidence_path=evidence_path,
                expected_bindings=expected_bindings,
            )
        )
        for name in _MANDATORY_CHECKS
    }
    for name in _MANDATORY_CHECKS:
        details[name]["metric_passed"] = _metric_passes(name, details[name])

    passed = all(
        bool(details[name]["passed"] and details[name]["metric_passed"])
        for name in _MANDATORY_CHECKS
    )
    return {
        "supplied": True,
        "passed": passed,
        "schema": PUBLICATION_EVIDENCE_SCHEMA,
        "bindings": expected_bindings,
        "final_node_ids": list(expected_ids),
        "checks": details,
    }


def assemble_publication_evidence(
    output_path: str | Path,
    *,
    config_sha256: str,
    archive_snapshot_sha256: str,
    model_content_sha256: str,
    final_node_ids: Sequence[str],
    report_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    """Create a portable aggregate from exact, independently retained reports.

    Source paths are stored relative to the aggregate.  Validation reopens each
    source and recomputes its exact byte SHA, so deleting, editing, or replacing
    a detailed report makes promotion fail closed.
    """

    output = Path(output_path)
    bindings = {
        "config_sha256": config_sha256,
        "archive_snapshot_sha256": archive_snapshot_sha256,
        "model_content_sha256": model_content_sha256,
    }
    if any(not _is_sha256(value) for value in bindings.values()):
        raise ValueError("publication-evidence bindings need lowercase SHA-256 values")
    final_ids = tuple(sorted(_normal_node_id(value) for value in final_node_ids))
    if len(final_ids) != len(set(final_ids)) or not final_ids:
        raise ValueError("publication evidence needs a nonempty unique final-node set")
    if set(report_paths) != set(_MANDATORY_CHECKS):
        raise ValueError("publication evidence needs exactly the four mandatory source reports")

    output.parent.mkdir(parents=True, exist_ok=True)
    checks: dict[str, Any] = {}
    for name in _MANDATORY_CHECKS:
        source_path = Path(report_paths[name])
        if source_path.resolve() == output.resolve():
            raise ValueError("aggregate publication evidence cannot overwrite a source report")
        if not source_path.is_file():
            raise ValueError(f"publication source report {name!r} is missing")
        source = _source_report(source_path, name, bindings)
        relative = os.path.relpath(source_path.resolve(), output.parent.resolve())
        contract = _SOURCE_REPORT_CONTRACTS[name]
        checks[name] = {
            "report_path": relative,
            "report_sha256": file_sha256(source_path),
            "passed": source["passed"],
            **{field: source[field] for field in contract["fields"]},
        }

    record = {
        "schema": PUBLICATION_EVIDENCE_SCHEMA,
        **bindings,
        "final_node_ids": list(final_ids),
        "checks": checks,
    }
    encoded = json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + "\n"
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(encoded, encoding="utf-8")
    temporary.replace(output)
    return validate_publication_evidence(
        output,
        config_sha256=config_sha256,
        archive_snapshot_sha256=archive_snapshot_sha256,
        model_content_sha256=model_content_sha256,
        final_node_ids=final_ids,
    )


def assemble_candidate_publication_evidence(
    candidate_payload: str | Path,
    archive_path: str | Path,
    campaign_manifest_path: str | Path,
    output_path: str | Path,
    *,
    report_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    """Derive immutable bindings and assemble evidence for one candidate."""

    from .archive import require_hdf5
    from .campaign import load_manifest
    from .fit import RuntimePayload

    output = Path(output_path)
    protected = {
        Path(candidate_payload).resolve(),
        Path(archive_path).resolve(),
        Path(campaign_manifest_path).resolve(),
    }
    if output.resolve() in protected:
        raise ValueError("publication evidence must not overwrite candidate, archive, or campaign")
    payload = RuntimePayload.load_binary(candidate_payload, allow_candidate=True)
    model_hash = payload.compute_model_content_sha256()
    if payload.model_content_sha256 != model_hash:
        raise ValueError("candidate model-content hash is inconsistent")
    campaign = load_manifest(campaign_manifest_path)
    final_ids = tuple(
        str(item["node_id_hex"]).lower().removeprefix("0x")
        for item in campaign["design"]["nodes"]["final"]
    )
    if len(final_ids) != 64 or len(set(final_ids)) != 64:
        raise ValueError("campaign must identify exactly 64 fresh final nodes")

    h5py, _ = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        config_hash = str(root["manifest"].attrs["config_sha256"])
        if campaign["config"]["sha256"] != config_hash:
            raise ValueError("campaign and archive configuration hashes differ")
        if payload.metadata.get("config_sha256") != config_hash:
            raise ValueError("candidate and archive configuration hashes differ")
        require_production_configuration(root, campaign)
        archive_ids = {f"{int(value):016x}" for value in root["nodes/id"][:]}
        if not set(final_ids).issubset(archive_ids):
            raise ValueError("archive does not contain the campaign's fresh final set")
        # Publication never rests on an incremental seal.
        snapshot_hash = archive_snapshot_sha256(root, require_full_verification=True)

    return assemble_publication_evidence(
        output_path,
        config_sha256=config_hash,
        archive_snapshot_sha256=snapshot_hash,
        model_content_sha256=model_hash,
        final_node_ids=final_ids,
        report_paths=report_paths,
    )
