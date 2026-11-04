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

"""End-to-end archive fitting, pilot resolution, and acceptance reporting."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .allocation import BLOCK_HISTORIES, effective_history_cap, predict_candidate_utilities
from .archive import require_hdf5
from .axes import registered_edges_from_physical
from .campaign import (
    MAXIMUM_FINAL_GENERATIONS,
    MODEL_SELECTION_HISTORY_SCHEMA,
    MODEL_SELECTION_REPORT_SCHEMA,
    campaign_round_identity,
    load_manifest,
)
from .fit import (
    ModelThresholds,
    RuntimePayload,
    score_model,
    select_rank_and_level,
)
from .noise import (
    EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
    EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
    EXACT_LINE_MATERIAL_FRACTION,
    EXACT_LINE_RELATIVE_ERROR_FRACTION,
    EXACT_LINE_RELATIVE_NOISE_FRACTION,
    EXACT_LINE_SIGNIFICANCE_SIGMAS,
    contribution_aware_model_selection_policy,
    direct_line_map,
    exact_line_materiality,
    replica_shape_tv_noise,
    sentinel_shape,
)
from .publication import (
    archive_snapshot_sha256,
    file_sha256,
    production_configuration_gate,
    require_production_configuration,
    validate_publication_evidence,
)


_STRUCTURAL_CONTINUUM_EXCLUSIONS = {
    "primary_order0": "analytic direct beam, not a fitted continuum",
    "annihilation_order0": "exact 511-keV atom, not a fitted continuum",
}
_ENERGY_AXIS_TO_RUNTIME_KIND = {
    "primary_registered": "primary",
    "annihilation_511_registered": "annihilation",
    "brems_log_relative": "brems",
}
_LINE_TO_RUNTIME_KIND = {
    "elastic_E0": "line_e0",
    "annihilation_511": "line_511",
}


def _apply_final_generation_cap(
    policy: Mapping[str, Any], campaign: Mapping[str, Any]
) -> dict[str, Any]:
    result = dict(policy)
    if result.get("action") != (
        "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve"
    ):
        return result
    generation = int(campaign["campaign_state"]["final_generation"])
    allocation = campaign.get("allocation", {})
    maximum = int(
        allocation.get("maximum_final_generations", MAXIMUM_FINAL_GENERATIONS)
    )
    if generation + 1 < maximum:
        return result
    return {
        "action": "start_new_campaign_revision",
        "reason": (
            "the preregistered three-generation final-test allowance is exhausted"
        ),
        "final_generation": generation,
        "maximum_final_generations": maximum,
        "unresolved_final_policy": result,
    }


_DETECTOR_RESPONSE_SCHEMA = "CEEFIT-DETECTOR-RESPONSE-v1"
_DETECTOR_RESPONSE_SEMANTICS = "probability-per-incident-cell-mass-v1"
_DETECTOR_LINE_INTERPOLATION = "linear-exact-energy-v1"
_DETECTOR_GEOMETRY_SCHEMA = "CEEFIT-DETECTOR-GEOMETRY-v1"
_DETECTOR_PROVENANCE_SCHEMA = "CEEFIT-DETECTOR-RESPONSE-PROVENANCE-v1"
_DETECTOR_MU_DEFINITION = "cosine-to-detector-axis;mu=1=head-on-v1"
_DETECTOR_DISCRETIZATION_METHOD = "archive-cell-response-kernel-diameter-v1"
_PHYSICAL_COMPONENT_ACCEPTANCE_SEMANTICS = (
    "continuum-origin-order-plus-exact-atom-sentinel3-v1"
)
PUBLICATION_PAYLOAD_SIZE_SCHEMA = "CEEFIT-RUNTIME-PAYLOAD-SIZE-v1"
# Advisory deployment-size target, reported only as a diagnostic.  Scientific
# acceptance is based on the fitted model and untouched-test gates, never on
# the size of this local/intermediate serialization.
PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES = 30 * 1024 * 1024


def _json_attr(value: Any) -> Any:
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return json.loads(value)


def _is_lower_sha256(value: Any) -> bool:
    return bool(
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _exact_binary_payload_size(payload: RuntimePayload, *, accepted: bool) -> int:
    """Measure the exact encoded artifact without writing or renaming a file."""

    previous = payload.accepted_for_publication
    payload.accepted_for_publication = accepted
    try:
        return len(payload.to_binary_bytes())
    finally:
        payload.accepted_for_publication = previous


def _bind_payload_size_report(
    payload: RuntimePayload,
    *,
    exact_input_candidate_bytes: int | None = None,
) -> dict[str, Any]:
    """Bind a fixed-point size report to the exact candidate/publication bytes.

    The report itself is serialized into the payload, so one preliminary size
    measurement would under-count it.  Iterate until both recorded byte counts
    describe the exact representation containing those same counts.  Diagnostic
    candidates are always serializable; callers decide whether a failed size
    report prevents setting the publication bit.
    """

    validation = payload.metadata.get("validation_metadata")
    if not isinstance(validation, dict):
        raise ValueError("runtime payload lacks mutable validation metadata")
    if exact_input_candidate_bytes is not None and (
        type(exact_input_candidate_bytes) is not int or exact_input_candidate_bytes < 0
    ):
        raise ValueError("serialized input-candidate byte count is invalid")
    record: dict[str, Any] = {
        "schema": PUBLICATION_PAYLOAD_SIZE_SCHEMA,
        "measurement": "exact-ceefit-payload-file-bytes-v1",
        "advisory_bytes": PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES,
        "advisory_mib": 30,
        "acceptance_policy": "diagnostic-only;not-a-scientific-gate",
    }
    validation["runtime_payload_size"] = record
    for _ in range(32):
        encoded_candidate_bytes = _exact_binary_payload_size(payload, accepted=False)
        proposed_publication_bytes = _exact_binary_payload_size(payload, accepted=True)
        input_candidate_bytes = (
            encoded_candidate_bytes
            if exact_input_candidate_bytes is None
            else exact_input_candidate_bytes
        )
        updated = {
            "schema": PUBLICATION_PAYLOAD_SIZE_SCHEMA,
            "measurement": "exact-ceefit-payload-file-bytes-v1",
            "advisory_bytes": PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES,
            "advisory_mib": 30,
            "acceptance_policy": "diagnostic-only;not-a-scientific-gate",
            "exact_input_candidate_bytes": input_candidate_bytes,
            "exact_bound_candidate_bytes": encoded_candidate_bytes,
            "proposed_publication_bytes": proposed_publication_bytes,
            "within_advisory_limit": bool(
                input_candidate_bytes <= PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES
                and encoded_candidate_bytes <= PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES
                and proposed_publication_bytes <= PUBLICATION_PAYLOAD_SIZE_ADVISORY_BYTES
            ),
            "passed": True,
        }
        if updated == record:
            return dict(record)
        record.clear()
        record.update(updated)
    raise RuntimeError("runtime payload size report did not reach a serialization fixed point")


def _write_fitted_payload(
    payload: RuntimePayload,
    output_payload: str | Path,
    validation_metadata: dict[str, Any],
    *,
    acceptance: bool,
) -> dict[str, Any]:
    """Record exact size, then write according to scientific acceptance only."""

    artifact_path = (
        Path(output_payload)
        if acceptance
        else Path(str(output_payload) + ".candidate")
    )
    validation_metadata["acceptance_passed"] = acceptance
    validation_metadata["payload_artifact"] = str(artifact_path)
    payload.metadata["validation_metadata"] = validation_metadata
    payload.accepted_for_publication = acceptance
    payload_size = _bind_payload_size_report(payload)
    validation_metadata["runtime_payload_size"] = payload_size
    payload.accepted_for_publication = acceptance
    payload.save_binary(artifact_path)
    return validation_metadata


def _mean_and_variance_of_mean(numpy: Any, n: Any, s1: Any, s2: Any) -> tuple[Any, Any]:
    count = numpy.asarray(n, dtype=float)
    first = numpy.asarray(s1, dtype=float)
    second = numpy.asarray(s2, dtype=float)
    mean = numpy.divide(first, count, out=numpy.zeros_like(first), where=count > 0.0)
    numerator = numpy.maximum(0.0, second - numpy.divide(first * first, count, out=numpy.zeros_like(first), where=count > 0.0))
    variance = numpy.divide(
        numerator,
        count * numpy.maximum(count - 1.0, 1.0),
        out=numpy.full_like(first, numpy.inf),
        where=count > 1.0,
    )
    return mean, variance


def _replica_yield_variance_of_mean(
    numpy: Any, replica_s1: Any, replica_histories: Any
) -> Any:
    """Estimate yield variance without discarding between-cell covariance.

    Each of the 32 logical replicas is an independent estimate of the complete
    rebinned yield.  Summing its cells before taking the across-replica sample
    variance preserves arbitrary same-history covariance in the stored map.
    """

    scores = numpy.asarray(replica_s1, dtype=float)
    histories = numpy.asarray(replica_histories, dtype=float)
    if scores.ndim < 3 or scores.shape[:2] != histories.shape:
        raise ValueError("replica yield statistics have incompatible shapes")
    if scores.shape[1] != 32 or numpy.any(histories <= 0.0):
        raise ValueError("yield variance requires 32 nonempty logical replicas")
    flattened = scores.reshape(scores.shape[0], scores.shape[1], -1)
    estimates = numpy.sum(flattened, axis=2) / histories
    return numpy.var(estimates, axis=1, ddof=1) / estimates.shape[1]


def _replica_yield_variance_from_archive(
    numpy: Any,
    replica_s1: Any,
    replica_histories: Any,
    start: int,
    stop: int,
    *,
    row_batch: int = 8,
) -> Any:
    """Compute replica yield variance while bounding archive read memory.

    Production replica maps are float32 and much larger than every other fit
    input.  Only their per-replica cell sums are needed for the component-yield
    covariance.  Read complete component rows in small batches so the summation
    and its floating-point order are identical to the in-memory calculation,
    without ever materializing the full node/component/replica tensor.
    """

    histories = numpy.asarray(replica_histories, dtype=float)
    if histories.ndim != 2 or histories.shape[1] != 32:
        raise ValueError("replica histories must be an N x 32 array")
    if replica_s1.ndim != 3 or replica_s1.shape[:2] != histories.shape:
        raise ValueError("archive replica yield statistics have incompatible shapes")
    if not (0 <= start < stop <= replica_s1.shape[2]):
        raise ValueError("archive replica component slice is invalid")
    if row_batch <= 0:
        raise ValueError("archive replica row batch must be positive")

    totals = numpy.empty(histories.shape, dtype=float)
    for begin in range(0, histories.shape[0], row_batch):
        end = min(begin + row_batch, histories.shape[0])
        values = numpy.asarray(
            replica_s1[begin:end, :, start:stop], dtype=float
        )
        totals[begin:end] = numpy.sum(values, axis=2)
    return _replica_yield_variance_of_mean(
        numpy, totals[:, :, None], histories
    )


def _archive_component_has_nonzero(
    numpy: Any,
    dataset: Any,
    start: int,
    stop: int,
    *,
    row_batch: int = 8,
) -> bool:
    """Return whether a component slice is nonzero without a full-tensor read."""

    if dataset.ndim not in {2, 3}:
        raise ValueError("archive component dataset must have two or three dimensions")
    if not (0 <= start < stop <= dataset.shape[-1]):
        raise ValueError("archive component slice is invalid")
    if row_batch <= 0:
        raise ValueError("archive component row batch must be positive")
    for begin in range(0, dataset.shape[0], row_batch):
        end = min(begin + row_batch, dataset.shape[0])
        selection = (
            dataset[begin:end, :, start:stop]
            if dataset.ndim == 3
            else dataset[begin:end, start:stop]
        )
        if numpy.any(selection != 0.0):
            return True
    return False


def _roles(root: Any) -> list[str]:
    return list(root["nodes/role"].asstr()[:])


def _campaign_stop_status(campaign: Mapping[str, Any]) -> dict[str, Any]:
    """Return the auditable two-pass model-selection stop state.

    The untouched final set must not be inspected until two consecutive
    model-selection rounds have passed.  Treat all missing, malformed, or
    internally inconsistent state as not stopped so a hand-edited manifest can
    never bypass that boundary.
    """

    state = campaign.get("campaign_state")
    if not isinstance(state, Mapping):
        return {
            "passed": False,
            "stopped": False,
            "consecutive_validation_passes": 0,
            "validation_history": [],
            "reason": "campaign manifest lacks adaptive validation state",
        }
    history = state.get("validation_history", [])
    if not isinstance(history, list) or any(not isinstance(value, Mapping) for value in history):
        return {
            "passed": False,
            "stopped": False,
            "consecutive_validation_passes": 0,
            "validation_history": [],
            "reason": "campaign validation history is malformed",
        }
    try:
        config_hash = str(campaign["config"]["sha256"])
        validation_ids = [
            str(item["node_id_hex"]).lower()
            for item in campaign["design"]["nodes"]["validation"]
        ]
        final_generation = int(state["final_generation"])
        round_number = int(state["round"])
    except (KeyError, TypeError, ValueError, OverflowError):
        return {
            "passed": False,
            "stopped": False,
            "consecutive_validation_passes": 0,
            "validation_history": [],
            "reason": "campaign validation identity is malformed",
        }

    def valid_sha256(value: Any) -> bool:
        return bool(
            isinstance(value, str)
            and len(value) == 64
            and all(character in "0123456789abcdef" for character in value)
        )

    flags: list[bool] = []
    first_round = round_number - len(history) + 1
    for offset, record in enumerate(history):
        flag = record.get("passed")
        kind = record.get("kind")
        common_valid = bool(
            record.get("schema") == MODEL_SELECTION_HISTORY_SCHEMA
            and type(flag) is bool
            and record.get("config_sha256") == config_hash
            and record.get("campaign_round") == first_round + offset
            and record.get("final_generation") == final_generation
            and valid_sha256(record.get("campaign_round_identity_sha256"))
        )
        if kind == "bound_model_selection_report":
            common_valid = bool(
                common_valid
                and record.get("report_schema") == MODEL_SELECTION_REPORT_SCHEMA
                and valid_sha256(record.get("report_canonical_sha256"))
                and valid_sha256(record.get("archive_snapshot_sha256"))
                and valid_sha256(record.get("model_content_sha256"))
                and record.get("eligible_validation_node_ids") == validation_ids
            )
        elif kind == "diagnostic_manual_failure":
            common_valid = common_valid and flag is False
        else:
            common_valid = False
        if not common_valid:
            return {
                "passed": False,
                "stopped": False,
                "consecutive_validation_passes": 0,
                "validation_history": list(history),
                "reason": "campaign validation history is malformed or unbound",
            }
        flags.append(flag)
    try:
        consecutive = int(state.get("consecutive_validation_passes", 0))
    except (TypeError, ValueError, OverflowError):
        consecutive = 0
    stopped = state.get("stopped") is True
    actual_consecutive = 0
    for flag in reversed(flags):
        if not flag:
            break
        actual_consecutive += 1
    passed = bool(
        stopped
        and consecutive >= 2
        and consecutive == actual_consecutive
        and len(flags) >= 2
        and flags[-2:] == [True, True]
    )
    return {
        "passed": passed,
        "stopped": stopped,
        "consecutive_validation_passes": consecutive,
        "validation_history": list(history),
        "reason": (
            "two consecutive model-selection rounds passed"
            if passed
            else "two consecutive passing model-selection rounds are required"
        ),
    }


def _transport_complete_mask(root: Any, numpy: Any) -> Any:
    statuses = numpy.asarray(root["nodes/status"].asstr()[:], dtype=object)
    failures = numpy.asarray(root["nodes/transport_failures"][:], dtype=int)
    capped = numpy.asarray(root["nodes/capped"][:], dtype=bool)
    return (
        numpy.isin(statuses, ("complete", "converged"))
        & (failures == 0)
        & ~capped
    )


def _transport_status_report(root: Any, indices: Any) -> dict[str, Any]:
    """Return a fail-closed transport assessment without reading fit truth.

    A maximum-vertex termination is an unusable transport result, not a
    statistical history-budget cap.  Keeping that distinction here prevents
    final-set policy from asking the scheduler to continue a node which can
    never become valid by appending histories.
    """

    statuses = list(root["nodes/status"].asstr()[:])
    failures = root["nodes/transport_failures"][:]
    capped = root["nodes/capped"][:]
    records: list[dict[str, Any]] = []
    for raw_index in indices:
        index = int(raw_index)
        status = str(statuses[index])
        transport_failures = int(failures[index])
        maximum_vertex_capped = bool(capped[index])
        eligible = bool(
            status in {"complete", "converged"}
            and transport_failures == 0
            and not maximum_vertex_capped
        )
        records.append(
            {
                "node_id": f"{int(root['nodes/id'][index]):016x}",
                "status": status,
                "transport_failures": transport_failures,
                "maximum_vertex_capped": maximum_vertex_capped,
                "eligible": eligible,
            }
        )
    failed_ids = [item["node_id"] for item in records if not item["eligible"]]
    return {
        "kind": "transport_status",
        "passed": bool(records) and not failed_ids,
        "node_count": len(records),
        "failed_node_ids": failed_ids,
        "nodes": records,
        "reason": (
            "all final nodes have transport-complete truth"
            if records and not failed_ids
            else "failed/max-vertex-capped final truth must be retired and replaced"
        ),
    }


def _component_offsets(root: Any, descriptors: Sequence[dict[str, Any]]) -> list[int]:
    n_mu = len(root["axes/mu_edges"]) - 1
    offsets = [0]
    for descriptor in descriptors:
        n_energy = root[
            f"axes/energy/{descriptor['energy_axis']}/physical_edges_keV"
        ].shape[1] - 1
        offsets.append(offsets[-1] + n_energy * n_mu)
    return offsets


def _accumulate_validation_model_loss(
    numpy: Any,
    losses: Any,
    indices: Any,
    model: Any,
    q: Any,
    truth_masses: Any,
    registered_energy_edges: Sequence[float],
    mu_edges: Sequence[float],
    *,
    reference_yields: Any | None = None,
    material_mask: Any | None = None,
    additive_baseline_masses: Any | None = None,
    weak_absolute_yield_fraction: float = 0.001,
    weak_absolute_shape_fraction: float = 0.001,
) -> None:
    """Accumulate dimensionless per-node loss for adaptive acquisition."""

    if not len(indices):
        return
    truth = numpy.asarray(truth_masses, dtype=float).reshape(len(indices), -1)
    truth_yield = numpy.sum(truth, axis=1)
    positive = truth_yield > 0.0
    if not bool(numpy.any(positive)):
        return
    predicted_yield, predicted_shape, _ = model.predict_cell_masses(
        q[positive], registered_energy_edges, mu_edges
    )
    positive_truth = truth[positive]
    predicted = predicted_yield[:, None] * predicted_shape
    evaluated_truth = positive_truth
    evaluated_prediction = predicted
    if additive_baseline_masses is not None:
        baseline = numpy.asarray(additive_baseline_masses, dtype=float).reshape(
            len(indices), -1
        )[positive]
        evaluated_truth = evaluated_truth + baseline
        evaluated_prediction = evaluated_prediction + baseline
    evaluated_truth_yield = numpy.sum(evaluated_truth, axis=1)
    evaluated_prediction_yield = numpy.sum(evaluated_prediction, axis=1)
    truth_shape = evaluated_truth / evaluated_truth_yield[:, None]
    predicted_evaluated_shape = (
        evaluated_prediction / evaluated_prediction_yield[:, None]
    )
    yield_ratio = numpy.abs(
        evaluated_prediction_yield / evaluated_truth_yield - 1.0
    ) / 0.005
    shape_ratio = 0.5 * numpy.sum(
        numpy.abs(predicted_evaluated_shape - truth_shape), axis=1
    ) / 0.01
    if reference_yields is not None:
        reference = numpy.asarray(reference_yields, dtype=float).reshape(-1)[positive]
        material = numpy.asarray(material_mask, dtype=bool).reshape(-1)[positive]
        weak_yield_ratio = (
            numpy.abs(predicted_yield - truth_yield[positive])
            / reference
            / weak_absolute_yield_fraction
        )
        weak_shape_ratio = (
            0.5 * numpy.sum(numpy.abs(predicted - positive_truth), axis=1)
            / reference
            / weak_absolute_shape_fraction
        )
        yield_ratio = numpy.where(material, yield_ratio, weak_yield_ratio)
        shape_ratio = numpy.where(material, shape_ratio, weak_shape_ratio)
    node_loss = numpy.maximum(yield_ratio, shape_ratio) ** 2
    selected_indices = numpy.asarray(indices)[positive]
    losses[selected_indices] = numpy.maximum(losses[selected_indices], node_loss)


def fit_archive(
    archive_path: str | Path,
    campaign_manifest_path: str | Path,
    output_payload: str | Path,
    *,
    ranks: Iterable[int] = (4, 8, 12, 16, 24),
    levels: Iterable[int] = (1, 2, 3, 4),
    nmf_iterations: int = 400,
    seed: int = 0xCEE10,
    detector_validation_path: str | Path | None = None,
    publication_evidence_path: str | Path | None = None,
    acquisition_utilities_path: str | Path | None = None,
) -> dict[str, Any]:
    """Fit every supported component/line by role and emit a binary payload."""

    h5py, numpy = require_hdf5()
    campaign = load_manifest(campaign_manifest_path)
    campaign_stop = _campaign_stop_status(campaign)
    brems_threshold = float(
        campaign["config"]["physics"]["brems_electron_floor_keV"]
    )
    q_record = campaign["design"]["q_domain"]
    q_domain = (
        numpy.asarray(q_record["min"], dtype=float),
        numpy.asarray(q_record["max"], dtype=float),
    )
    bounds = campaign["design"]["bounds"]
    physical_domain = {
        "min": [
            bounds["energy_keV"][0],
            bounds["distance_m"][0],
            0.0 if campaign["config"]["mode"] == "air" else bounds["source_height_m"][0],
            0.0 if campaign["config"]["mode"] == "air" else bounds["detector_height_m"][0],
        ],
        "max": [
            bounds["energy_keV"][1],
            bounds["distance_m"][1],
            (0.0 if campaign["config"]["mode"] == "air" else bounds["source_height_m"][1]),
            (0.0 if campaign["config"]["mode"] == "air" else bounds["detector_height_m"][1]),
        ],
    }
    with h5py.File(archive_path, "r") as root:
        archive_hash = str(root["manifest"].attrs["config_sha256"])
        if archive_hash != campaign["config"]["sha256"]:
            raise ValueError("campaign and HDF5 archive configuration hashes differ")
        production_configuration = production_configuration_gate(root, campaign)
        q_all = numpy.asarray(root["nodes/q"][:], dtype=float)
        if numpy.any(q_all < q_domain[0] - 1.0e-12) or numpy.any(q_all > q_domain[1] + 1.0e-12):
            raise ValueError("archive contains a node outside campaign q_domain")
        roles = _roles(root)
        archive_ids = {
            f"{int(value):016x}": index for index, value in enumerate(root["nodes/id"][:])
        }
        state = campaign.get("campaign_state", {})
        train_ids = list(
            state.get(
                "training_node_ids",
                [item["node_id_hex"] for item in campaign["design"]["nodes"]["train"]],
            )
        )
        validation_ids = [
            item["node_id_hex"] for item in campaign["design"]["nodes"]["validation"]
        ]
        final_ids = [item["node_id_hex"] for item in campaign["design"]["nodes"]["final"]]
        if len(final_ids) != len(set(final_ids)) or not final_ids:
            raise ValueError("current campaign manifest has empty or duplicate final nodes")
        if production_configuration["passed"] and len(final_ids) != 64:
            raise ValueError(
                "production campaign manifest must identify exactly 64 fresh final nodes"
            )
        if set(train_ids) & set(validation_ids) or set(train_ids) & set(final_ids) or set(validation_ids) & set(final_ids):
            raise ValueError("campaign role overlay contains overlapping node IDs")
        required_ids = set(train_ids) | set(validation_ids) | set(final_ids)
        missing = sorted(required_ids - set(archive_ids))
        if missing:
            raise ValueError(
                f"archive is missing {len(missing)} nodes required by the current campaign role overlay"
            )
        requested_train_indices = numpy.asarray([archive_ids[value] for value in train_ids], dtype=int)
        requested_validation_indices = numpy.asarray(
            [archive_ids[value] for value in validation_ids], dtype=int
        )
        final_indices = numpy.asarray([archive_ids[value] for value in final_ids], dtype=int)
        transport_complete = _transport_complete_mask(root, numpy)
        train_indices = requested_train_indices[transport_complete[requested_train_indices]]
        validation_indices = requested_validation_indices[
            transport_complete[requested_validation_indices]
        ]
        if not train_indices.size or not validation_indices.size or not final_indices.size:
            raise ValueError("archive must contain train, validation, and untouched final roles")
        validation_model_loss = numpy.zeros(len(roles), dtype=float)

        descriptors = _json_attr(root["components"].attrs["descriptors_json"])
        continuum_inventory = []
        component_axis_kinds: dict[str, str] = {}
        for descriptor in descriptors:
            if descriptor["name"] in _STRUCTURAL_CONTINUUM_EXCLUSIONS:
                continue
            axis_index = descriptor["energy_axis"]
            axis_name = str(root[f"axes/energy/{axis_index}"].attrs["name"])
            try:
                axis_kind = _ENERGY_AXIS_TO_RUNTIME_KIND[axis_name]
            except KeyError as error:
                raise ValueError(f"unknown archive energy axis {axis_name!r}") from error
            component_axis_kinds[descriptor["name"]] = axis_kind
            continuum_inventory.append(
                {"name": descriptor["name"], "axis_kind": axis_kind}
            )
        offsets = _component_offsets(root, descriptors)
        support_offsets = root["components/support_offsets"][:]
        support = root["components/support_mask"][:]
        n_mu = len(root["axes/mu_edges"]) - 1
        mu_edges = numpy.asarray(root["axes/mu_edges"][:], dtype=float)
        continuum_n = root["continuum/components/N"]
        continuum_s1 = root["continuum/components/S1"]
        continuum_s2 = root["continuum/components/S2"]
        continuum_replica = root["continuum/components/replica_S1"]
        replica_histories = root["replicas/histories"][:]
        direct_fluence = numpy.asarray(root["nodes/direct_fluence"][:], dtype=float)
        if numpy.any(direct_fluence <= 0.0):
            raise ValueError("archive contains a non-positive deterministic direct fluence")
        scattered_n = numpy.asarray(root["diagnostics/scattered_total_N"][:], dtype=float)
        scattered_s1 = numpy.asarray(root["diagnostics/scattered_total_S1"][:], dtype=float)
        scattered_reference_yields = numpy.divide(
            scattered_s1,
            scattered_n,
            out=numpy.zeros_like(scattered_s1),
            where=scattered_n > 0.0,
        ) / direct_fluence
        if numpy.any(~numpy.isfinite(scattered_reference_yields)) or numpy.any(
            scattered_reference_yields <= 0.0
        ):
            raise ValueError("archive contains a non-positive scattered reference yield")
        component_models = {}
        excluded_train_ids = [
            f"{int(root['nodes/id'][index]):016x}"
            for index in requested_train_indices
            if not transport_complete[index]
        ]
        excluded_validation_ids = [
            f"{int(root['nodes/id'][index]):016x}"
            for index in requested_validation_indices
            if not transport_complete[index]
        ]
        selection_report: dict[str, Any] = {
            "transport_status": {
                "kind": "transport_status",
                "validation_passed": not excluded_validation_ids,
                "excluded_train_node_ids": excluded_train_ids,
                "excluded_validation_node_ids": excluded_validation_ids,
                "reason": "failed/max-vertex-capped transport is excluded from fit truth",
            }
        }
        for component_index, descriptor in enumerate(descriptors):
            if descriptor["name"] in _STRUCTURAL_CONTINUUM_EXCLUSIONS:
                start, stop = offsets[component_index], offsets[component_index + 1]
                nonzero = []
                for name, dataset in (
                    ("S1", continuum_s1),
                    ("S2", continuum_s2),
                    ("replica_S1", continuum_replica),
                ):
                    if _archive_component_has_nonzero(
                        numpy, dataset, start, stop
                    ):
                        nonzero.append(name)
                if nonzero:
                    raise ValueError(
                        f"structurally excluded continuum {descriptor['name']!r} "
                        f"contains nonzero archive data in {', '.join(nonzero)}"
                    )
                selection_report[descriptor["name"]] = {
                    "kind": "structural_exclusion",
                    "validation_passed": True,
                    "archive_zero_verified": True,
                    "reason": _STRUCTURAL_CONTINUUM_EXCLUSIONS[descriptor["name"]],
                }
                continue
            support_slice = support[
                :, support_offsets[component_index] : support_offsets[component_index + 1]
            ]
            supported = numpy.any(support_slice != 0, axis=1)
            start, stop = offsets[component_index], offsets[component_index + 1]
            component_n = continuum_n[:, start:stop]
            component_s1 = continuum_s1[:, start:stop]
            component_s2 = continuum_s2[:, start:stop]
            means, variances = _mean_and_variance_of_mean(
                numpy,
                component_n,
                component_s1,
                component_s2,
            )
            del component_n, component_s1, component_s2
            means = means / direct_fluence[:, None]
            variances = variances / direct_fluence[:, None] ** 2
            yield_variances = _replica_yield_variance_from_archive(
                numpy,
                continuum_replica,
                replica_histories,
                start,
                stop,
            ) / direct_fluence**2
            positive = numpy.sum(means, axis=1) > 0.0
            train = train_indices[supported[train_indices] & positive[train_indices]]
            validation = validation_indices[
                supported[validation_indices] & positive[validation_indices]
            ]
            if train.size < 2 or validation.size < 1:
                validation_model_loss[
                    validation_indices[supported[validation_indices]]
                ] = numpy.maximum(
                    validation_model_loss[
                        validation_indices[supported[validation_indices]]
                    ],
                    4.0,
                )
                selection_report[descriptor["name"]] = {
                    "kind": "continuum",
                    "validation_passed": False,
                    "reason": "insufficient supported positive train/validation nodes",
                    "train_nodes": int(train.size),
                    "validation_nodes": int(validation.size),
                }
                continue
            axis_index = descriptor["energy_axis"]
            registered_edges = root[f"axes/energy/{axis_index}/registered_edges"][:]
            axis_kind = component_axis_kinds[descriptor["name"]]
            validation_component_yields = numpy.sum(means[validation], axis=1)
            validation_component_se = numpy.sqrt(
                numpy.maximum(yield_variances[validation], 0.0)
            )
            validation_reference = scattered_reference_yields[validation]
            validation_material = (
                validation_component_yields
                + EXACT_LINE_SIGNIFICANCE_SIGMAS * validation_component_se
                >= 0.05 * validation_reference
            )
            selected, candidates = select_rank_and_level(
                descriptor["name"],
                axis_kind,
                q_all[train],
                means[train],
                q_all[validation],
                means[validation],
                registered_energy_edges=registered_edges,
                mu_edges=mu_edges,
                train_variances=variances[train],
                train_yield_variances=yield_variances[train],
                ranks=ranks,
                levels=levels,
                nmf_iterations=nmf_iterations,
                seed=seed + component_index * 1009,
                source_energy_threshold_keV=(
                    brems_threshold if axis_kind == "brems" else None
                ),
                q_domain=q_domain,
                validation_reference_yields=validation_reference,
                validation_material_mask=validation_material,
                weak_absolute_yield_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                weak_absolute_shape_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
            )
            component_models[descriptor["name"]] = selected.model
            _accumulate_validation_model_loss(
                numpy,
                validation_model_loss,
                validation,
                selected.model,
                q_all[validation],
                means[validation],
                registered_edges,
                mu_edges,
                reference_yields=validation_reference,
                material_mask=validation_material,
                weak_absolute_yield_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                weak_absolute_shape_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
            )
            selection_report[descriptor["name"]] = {
                "kind": "continuum",
                "rank": selected.model.rank,
                "level": selected.model.level,
                "parameter_count": selected.model.parameter_count,
                "validation_passed": selected.passed,
                "max_yield_relative_error": selected.max_yield_relative_error,
                "max_joint_tv": selected.max_joint_tv,
                "candidate_count": len(candidates),
                "validation_policy": {
                    "kind": "component-contribution-aware-v1",
                    "material_fraction_3sigma": 0.05,
                    "material_node_count": int(numpy.count_nonzero(validation_material)),
                    "weak_absolute_yield_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                    "weak_absolute_shape_l1_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                },
            }

        line_definitions = _json_attr(root["axes"].attrs["line_definitions_json"])
        line_inventory = []
        for descriptor in line_definitions:
            try:
                axis_kind = _LINE_TO_RUNTIME_KIND[descriptor["name"]]
            except KeyError as error:
                raise ValueError(
                    f"unknown exact-line model {descriptor['name']!r}"
                ) from error
            line_inventory.append(
                {"name": descriptor["name"], "axis_kind": axis_kind}
            )
        line_support = root["lines/support_mask"][:]
        line_n = root["lines/total/N"]
        line_s1 = root["lines/total/S1"]
        line_s2 = root["lines/total/S2"]
        line_replica = root["lines/total/replica_S1"]
        line_models = {}
        for line_index, descriptor in enumerate(line_definitions):
            line_start = line_index * n_mu
            line_stop = line_start + n_mu
            line_n_values = line_n[:, line_start:line_stop]
            line_s1_values = line_s1[:, line_start:line_stop]
            line_s2_values = line_s2[:, line_start:line_stop]
            means, variances = _mean_and_variance_of_mean(
                numpy,
                line_n_values,
                line_s1_values,
                line_s2_values,
            )
            del line_n_values, line_s1_values, line_s2_values
            means = means / direct_fluence[:, None]
            variances = variances / direct_fluence[:, None] ** 2
            yield_variances = _replica_yield_variance_from_archive(
                numpy,
                line_replica,
                replica_histories,
                line_start,
                line_stop,
            ) / direct_fluence**2
            positive = numpy.sum(means, axis=1) > 0.0
            train = train_indices[(line_support[train_indices, line_index] != 0) & positive[train_indices]]
            validation = validation_indices[
                (line_support[validation_indices, line_index] != 0) & positive[validation_indices]
            ]
            if train.size < 2 or validation.size < 1:
                supported_validation = validation_indices[
                    line_support[validation_indices, line_index] != 0
                ]
                validation_model_loss[supported_validation] = numpy.maximum(
                    validation_model_loss[supported_validation], 4.0
                )
                selection_report[descriptor["name"]] = {
                    "kind": "line",
                    "validation_passed": False,
                    "reason": "insufficient supported positive train/validation nodes",
                    "train_nodes": int(train.size),
                    "validation_nodes": int(validation.size),
                }
                continue
            axis_kind = _LINE_TO_RUNTIME_KIND[descriptor["name"]]
            threshold = None if axis_kind == "line_e0" else 2.0 * 510.998950
            feasible_ranks = tuple(rank for rank in ranks if rank <= min(train.size, n_mu))
            line_name = str(descriptor["name"])
            validation_line_yields = numpy.sum(means[validation], axis=1)
            validation_line_se = numpy.sqrt(
                numpy.maximum(yield_variances[validation], 0.0)
            )
            validation_reference = None
            validation_material = None
            validation_baseline = None
            if line_name == "elastic_E0":
                validation_baseline = numpy.zeros_like(means[validation])
                validation_baseline[:, -1] = 1.0
                line_validation_policy = {
                    "kind": "combined-direct-plus-elastic-relative-v1",
                    "material_node_count": int(validation.size),
                }
            else:
                validation_reference = scattered_reference_yields[validation]
                validation_material = (
                    validation_line_yields
                    + EXACT_LINE_SIGNIFICANCE_SIGMAS * validation_line_se
                    >= EXACT_LINE_MATERIAL_FRACTION * validation_reference
                )
                line_validation_policy = {
                    "kind": "exact-line-contribution-aware-v1",
                    "material_fraction_3sigma": EXACT_LINE_MATERIAL_FRACTION,
                    "material_node_count": int(
                        numpy.count_nonzero(validation_material)
                    ),
                    "weak_absolute_yield_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                    "weak_absolute_angular_l1_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                }
            selected, candidates = select_rank_and_level(
                descriptor["name"],
                axis_kind,
                q_all[train],
                means[train],
                q_all[validation],
                means[validation],
                registered_energy_edges=(0.0, 1.0),
                mu_edges=mu_edges,
                train_variances=variances[train],
                train_yield_variances=yield_variances[train],
                ranks=feasible_ranks,
                levels=levels,
                nmf_iterations=nmf_iterations,
                seed=seed + 50_000 + line_index,
                source_energy_threshold_keV=threshold,
                q_domain=q_domain,
                validation_reference_yields=validation_reference,
                validation_material_mask=validation_material,
                additive_baseline_masses=validation_baseline,
                weak_absolute_yield_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                weak_absolute_shape_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
            )
            line_models[descriptor["name"]] = selected.model
            _accumulate_validation_model_loss(
                numpy,
                validation_model_loss,
                validation,
                selected.model,
                q_all[validation],
                means[validation],
                (0.0, 1.0),
                mu_edges,
                reference_yields=validation_reference,
                material_mask=validation_material,
                additive_baseline_masses=validation_baseline,
                weak_absolute_yield_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
                weak_absolute_shape_fraction=EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
            )
            selection_report[descriptor["name"]] = {
                "kind": "line",
                "rank": selected.model.rank,
                "level": selected.model.level,
                "parameter_count": selected.model.parameter_count,
                "validation_passed": selected.passed,
                "max_yield_relative_error": selected.max_yield_relative_error,
                "max_joint_tv": selected.max_joint_tv,
                "candidate_count": len(candidates),
                "validation_policy": line_validation_policy,
            }

        selection_report["mc_noise"] = _validation_mc_noise_report(
            root,
            validation_indices,
            descriptors,
            offsets,
        )

        pool = campaign["design"]["candidate_pool"]
        unused_indices = state.get(
            "unused_candidate_indices",
            [
                int(item["candidate_index"])
                for item in pool
                if item.get("selected_role") == "candidate"
            ],
        )
        unused_candidates = {
            str(pool[int(index)]["node_id_hex"]): pool[int(index)]["q"]
            for index in unused_indices
        }
        histories_by_node = numpy.asarray(root["nodes/histories"][:], dtype=float)
        cpu_by_node = numpy.asarray(root["nodes/cpu_seconds"][:], dtype=float)
        wall_by_node = numpy.asarray(root["nodes/wall_seconds"][:], dtype=float)
        seconds_per_block = []
        for index in train_indices:
            blocks = histories_by_node[index] / BLOCK_HISTORIES
            timing = cpu_by_node[index] / blocks if blocks > 0.0 else 0.0
            if not math.isfinite(timing) or timing <= 0.0:
                timing = wall_by_node[index] / blocks if blocks > 0.0 else 0.0
            seconds_per_block.append(timing)
        positive_timings = [
            value for value in seconds_per_block if math.isfinite(value) and value > 0.0
        ]
        fallback_timing = (
            float(numpy.median(positive_timings)) if positive_timings else 1.0
        )
        seconds_per_block = [
            value if math.isfinite(value) and value > 0.0 else fallback_timing
            for value in seconds_per_block
        ]
        acquisition_predictions = predict_candidate_utilities(
            q_all[train_indices].tolist(),
            seconds_per_block,
            q_all[validation_indices].tolist(),
            validation_model_loss[validation_indices].tolist(),
            unused_candidates,
            q_min=q_domain[0].tolist(),
            q_max=q_domain[1].tolist(),
        )
        train_local_predictions = predict_candidate_utilities(
            q_all[train_indices].tolist(),
            seconds_per_block,
            q_all[validation_indices].tolist(),
            validation_model_loss[validation_indices].tolist(),
            {
                f"{int(root['nodes/id'][index]):016x}": q_all[index].tolist()
                for index in train_indices
            },
            q_min=q_domain[0].tolist(),
            q_max=q_domain[1].tolist(),
        )
        continuation_local_losses = {
            node_id: values["local_validation_loss"]
            for node_id, values in train_local_predictions.items()
        }
        continuation_local_losses.update(
            {
                f"{int(root['nodes/id'][index]):016x}": float(
                    validation_model_loss[index]
                )
                for index in validation_indices
            }
        )
        acquisition_utilities = {
            node_id: values["utility"]
            for node_id, values in acquisition_predictions.items()
        }
        ranked_acquisition = sorted(
            acquisition_predictions.items(),
            key=lambda item: (-item[1]["utility"], item[0]),
        )

        # The fit settings and output semantics needed to reproduce and consume
        # this result; the configuration itself is bound by config_sha256.
        fit_parameters = {
            "input_parameters": {
                "canonical_configuration": _json_attr(
                    root["manifest"].attrs["canonical_config_json"]
                ),
                "campaign_design_bounds": campaign["design"]["bounds"],
                "campaign_sobol_seed": campaign["design"]["sobol_seed"],
                "fit": {
                    "ranks": list(ranks),
                    "levels": list(levels),
                    "nmf_iterations": nmf_iterations,
                    "seed": seed,
                },
            },
            "output_parameters": {
                "observable": (
                    "incident-continuum-cell-masses-plus-exact-line-atoms"
                ),
                "normalization": "s_hat_to_direct_fluence",
                "energy_cut_keV": float(
                    _json_attr(root["manifest"].attrs["canonical_config_json"])[
                        "axis"
                    ]["energy_cut_keV"]
                ),
                "mu_definition": "cosine-to-detector-axis;mu=1=head-on-v1",
                "requested_edges": (
                    "runtime evaluator integrates onto caller-supplied energy "
                    "and mu edges"
                ),
                "models": {
                    name: {
                        "axis_kind": model.axis_kind,
                        "rank": model.rank,
                        "geometry_spline_level": model.level,
                        "parameter_count": model.parameter_count,
                    }
                    for name, model in sorted(
                        {**component_models, **line_models}.items()
                    )
                },
            },
        }
        payload = RuntimePayload(
            component_models,
            line_models,
            {
                "config_sha256": archive_hash,
                "material_name": str(root["manifest"].attrs["material_name"]),
                "geometry_tier": str(root["manifest"].attrs["geometry_tier"]),
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": physical_domain,
                "q_domain": {"min": q_domain[0].tolist(), "max": q_domain[1].tolist()},
                "brems_activation_threshold_keV": brems_threshold,
                "model_inventory": {
                    "continuum": continuum_inventory,
                    "lines": line_inventory,
                },
            },
        )
        candidate_model_hash = payload.compute_model_content_sha256()
        payload.model_content_sha256 = candidate_model_hash
        archive_snapshot_hash = archive_snapshot_sha256(root)
        round_identity = campaign_round_identity(campaign)
        adaptive_utility_report_sha256 = None
        if acquisition_utilities_path is not None:
            destination = Path(acquisition_utilities_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            utility_report = {
                "schema": "CEEFIT-ADAPTIVE-UTILITY-v1",
                "config_sha256": archive_hash,
                "archive_snapshot_sha256": archive_snapshot_hash,
                "model_content_sha256": candidate_model_hash,
                "campaign_round": round_identity["campaign_round"],
                "campaign_round_identity_sha256": round_identity["identity_sha256"],
                "objective": (
                    "threshold-normalized validation model loss; predicted loss "
                    "reduction per measured CPU second"
                ),
                "new_nodes": acquisition_predictions,
                "continuation_local_validation_losses": continuation_local_losses,
                "input_parameters": {
                    "q_domain": {
                        "min": q_domain[0].tolist(),
                        "max": q_domain[1].tolist(),
                    },
                    "initial_blocks_per_new_node": 4,
                    "training_node_ids": round_identity["training_node_ids"],
                    "validation_node_ids": round_identity["validation_node_ids"],
                },
                "output_parameters": {
                    "new_node_fields": sorted(next(iter(acquisition_predictions.values())).keys())
                    if acquisition_predictions
                    else [],
                    "continuation_value": "local_validation_loss",
                },
            }
            temporary = destination.with_name(destination.name + ".tmp")
            with temporary.open("w", encoding="utf-8") as stream:
                json.dump(utility_report, stream, indent=2, sort_keys=True, allow_nan=False)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            temporary.replace(destination)
            adaptive_utility_report_sha256 = hashlib.sha256(destination.read_bytes()).hexdigest()
        if campaign_stop["passed"]:
            final_report = _score_final_set(
                root,
                payload,
                final_indices,
                descriptors,
                offsets,
                detector_validation_path=detector_validation_path,
                config_hash=archive_hash,
                model_content_sha256=candidate_model_hash,
                archive_snapshot_hash=archive_snapshot_hash,
            )
            if (
                final_report["transport_status"]["passed"]
                and production_configuration["passed"]
            ):
                publication_evidence = validate_publication_evidence(
                    publication_evidence_path,
                    config_sha256=archive_hash,
                    archive_snapshot_sha256=archive_snapshot_hash,
                    model_content_sha256=candidate_model_hash,
                    final_node_ids=final_ids,
                )
            elif not production_configuration["passed"]:
                publication_evidence = {
                    "supplied": False,
                    "passed": False,
                    "reason": (
                        "publication evidence is not consumed for a "
                        "diagnostic-only configuration"
                    ),
                }
            else:
                publication_evidence = {
                    "supplied": False,
                    "passed": False,
                    "reason": (
                        "publication evidence is not consumed for a final set "
                        "with invalid transport truth"
                    ),
                }
        else:
            # Looking at final-node residuals before the adaptive stop would
            # turn the nominally untouched set into another validation set.
            final_report = {
                "evaluated": False,
                "passed": False,
                "physics_passed": False,
                "nodes": [],
                "detector_validation": {
                    "supplied": False,
                    "passed": False,
                    "reason": "final set is locked until the two-pass adaptive stop",
                },
                "reason": "untouched final set was not evaluated",
            }
            publication_evidence = {
                "supplied": False,
                "passed": False,
                "reason": "publication evidence is checked only after the adaptive stop",
            }
        selection_passed = bool(selection_report) and all(
            item["validation_passed"] for item in selection_report.values()
        )
        acceptance = bool(
            production_configuration["passed"]
            and selection_passed
            and campaign_stop["passed"]
            and final_report["passed"]
            and publication_evidence["passed"]
        )
        if acceptance:
            policy = {"action": "publish_payload"}
        elif not production_configuration["passed"]:
            policy = {
                "action": "retain_diagnostic_candidate_only",
                "reason": (
                    "compact or altered configurations may be fitted for diagnostics "
                    "but cannot become runtime production payloads"
                ),
                "issues": production_configuration["issues"],
            }
        elif (
            campaign_stop["passed"]
            and not final_report.get("transport_status", {}).get("passed", False)
        ):
            policy = {
                "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                "reason": (
                    "failed/max-vertex-capped transport is unusable and cannot "
                    "be repaired by same-node continuation"
                ),
                "retired_node_ids": [
                    f"{int(root['nodes/id'][index]):016x}" for index in final_indices
                ],
                "transport_failed_node_ids": final_report.get(
                    "transport_status", {}
                ).get("failed_node_ids", []),
                "fresh_final_count": 64,
                "campaign_command": "refresh-final",
            }
        elif not selection_passed:
            policy = {
                "action": "continue_adaptive_training_or_histories",
                "reason": "model-selection or validation-MC gates have not passed",
                "acquisition_utilities": (
                    None
                    if acquisition_utilities_path is None
                    else str(acquisition_utilities_path)
                ),
            }
        elif not campaign_stop["passed"]:
            policy = {
                "action": "record_passing_validation_round",
                "reason": campaign_stop["reason"],
                "campaign_command": (
                    "campaign-round --validation-report <fit-report.json> "
                    "--archive <archive.h5> --candidate <payload.candidate>"
                ),
                "final_set_evaluated": False,
            }
        elif not final_report.get("gates", {}).get("mc_noise_eligible", False):
            noisy_nodes = [
                item
                for item in final_report.get("nodes", ())
                if item.get("mc_noise_eligible") is not True
            ]
            precision_capped = [
                str(item["node_id"])
                for item in noisy_nodes
                if item.get("history_budget_exhausted") is True
            ]
            if precision_capped:
                policy = {
                    "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                    "reason": (
                        "transport-complete truth remains usable, but its statistical "
                        "history cap was reached before the MC-noise gate"
                    ),
                    "retired_node_ids": [
                        f"{int(root['nodes/id'][index]):016x}" for index in final_indices
                    ],
                    "precision_capped_node_ids": precision_capped,
                    "fresh_final_count": 64,
                    "campaign_command": "refresh-final",
                }
            else:
                policy = {
                    "action": "continue_final_histories_for_mc_noise",
                    "reason": (
                        "untouched-final residuals cannot be adjudicated until all "
                        "archive uncertainty targets pass"
                    ),
                    "node_ids": [str(item["node_id"]) for item in noisy_nodes],
                    "campaign_command": "continue-archive --roles final",
                }
        elif not final_report["physics_passed"]:
            policy = {
                "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                "retired_node_ids": [f"{int(root['nodes/id'][index]):016x}" for index in final_indices],
                "fresh_final_count": 64,
                "campaign_command": "refresh-final",
            }
        elif not final_report["detector_validation"]["supplied"]:
            policy = {
                "action": "validate_candidate_detector_folds",
                "candidate_model_content_sha256": candidate_model_hash,
                "archive_snapshot_sha256": archive_snapshot_hash,
            }
        elif not final_report["detector_validation"].get("noise_passed", False):
            detector_noise_ids = final_report["detector_validation"].get(
                "noise_failing_node_ids", []
            )
            capped_ids = {
                str(item["node_id"])
                for item in final_report.get("nodes", ())
                if item.get("history_budget_exhausted") is True
            }
            precision_capped = sorted(set(detector_noise_ids) & capped_ids)
            if precision_capped:
                policy = {
                    "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                    "reason": (
                        "detector-fold replica noise remains above tolerance at "
                        "the statistical history cap"
                    ),
                    "retired_node_ids": [
                        f"{int(root['nodes/id'][index]):016x}" for index in final_indices
                    ],
                    "precision_capped_node_ids": precision_capped,
                    "fresh_final_count": 64,
                    "campaign_command": "refresh-final",
                }
            else:
                policy = {
                    "action": "continue_final_histories_for_detector_fold_noise",
                    "reason": (
                        "detector-fold residuals cannot retire the final set while "
                        "their replica uncertainty exceeds 0.25%"
                    ),
                    "node_ids": detector_noise_ids,
                    "campaign_command": (
                        "continue-archive --roles final --detector-validation "
                        "<detector-fold-report.json>"
                    ),
                }
        elif not final_report["detector_validation"]["passed"]:
            policy = {
                "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                "retired_node_ids": [f"{int(root['nodes/id'][index]):016x}" for index in final_indices],
                "fresh_final_count": 64,
                "campaign_command": "refresh-final",
            }
        else:
            policy = {
                "action": "collect_or_correct_publication_evidence",
                "candidate_model_content_sha256": candidate_model_hash,
                "archive_snapshot_sha256": archive_snapshot_hash,
                "required_checks": publication_evidence.get("required_checks", []),
            }
        policy = _apply_final_generation_cap(policy, campaign)
        validation_metadata = {
            "schema": MODEL_SELECTION_REPORT_SCHEMA,
            "model_selection_policy": contribution_aware_model_selection_policy(),
            "fit_parameters": fit_parameters,
            "config_sha256": archive_hash,
            "archive_snapshot_sha256": archive_snapshot_hash,
            "model_content_sha256": candidate_model_hash,
            "campaign_round": round_identity["campaign_round"],
            "final_generation": round_identity["final_generation"],
            "campaign_round_identity_sha256": round_identity["identity_sha256"],
            "eligible_validation_node_ids": round_identity["validation_node_ids"],
            "passed": selection_passed,
            "selection": selection_report,
            "production_configuration": production_configuration,
            "campaign_stop": campaign_stop,
            "final": final_report,
            "acquisition": {
                "schema": "CEEFIT-ACQUISITION-UTILITY-v1",
                "criterion": "matern32_validation_loss_reduction_per_cpu_second",
                "candidate_count": len(acquisition_predictions),
                "utilities_artifact": (
                    None
                    if acquisition_utilities_path is None
                    else str(acquisition_utilities_path)
                ),
                "utilities_artifact_sha256": adaptive_utility_report_sha256,
                "top_candidates": [
                    {"node_id": node_id, **values}
                    for node_id, values in ranked_acquisition[:64]
                ],
            },
            "acceptance_passed": acceptance,
            "publication_evidence": publication_evidence,
            "payload_artifact": (
                str(output_payload)
                if acceptance
                else str(output_payload) + ".candidate"
            ),
            "policy": policy,
        }
        return _write_fitted_payload(
            payload,
            output_payload,
            validation_metadata,
            acceptance=acceptance,
        )


def promote_candidate(
    candidate_payload: str | Path,
    archive_path: str | Path,
    campaign_manifest_path: str | Path,
    detector_validation_path: str | Path,
    publication_evidence_path: str | Path,
    output_payload: str | Path,
) -> dict[str, Any]:
    """Promote a validated candidate without repeating the expensive fit.

    Only validation/publication metadata and the publication bit change.  The
    serialized model-content hash excludes both and is rechecked against every
    external report before the accepted artifact is written atomically.
    """

    candidate = Path(candidate_payload)
    output = Path(output_payload)
    if candidate.resolve() == output.resolve():
        raise ValueError("candidate and accepted output paths must be distinct")
    payload = RuntimePayload.load_binary(candidate, allow_candidate=True)
    if payload.accepted_for_publication:
        raise ValueError("input payload is already accepted for publication")
    model_hash = payload.compute_model_content_sha256()
    if payload.model_content_sha256 != model_hash:
        raise ValueError("candidate model-content hash is inconsistent")

    campaign = load_manifest(campaign_manifest_path)
    campaign_stop = _campaign_stop_status(campaign)
    if not campaign_stop["passed"]:
        raise ValueError(
            "candidate cannot be promoted before two consecutive passing "
            "model-selection rounds"
        )
    final_ids = [
        str(item["node_id_hex"]).lower().removeprefix("0x")
        for item in campaign["design"]["nodes"]["final"]
    ]
    if len(final_ids) != 64 or len(set(final_ids)) != 64:
        raise ValueError("campaign must identify exactly 64 fresh final nodes")

    h5py, _ = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        config_hash = str(root["manifest"].attrs["config_sha256"])
        if config_hash != campaign["config"]["sha256"]:
            raise ValueError("campaign and archive configuration hashes differ")
        if payload.metadata.get("config_sha256") != config_hash:
            raise ValueError("candidate and archive configuration hashes differ")
        production_configuration = require_production_configuration(root, campaign)
        archive_ids = {
            f"{int(value):016x}": index
            for index, value in enumerate(root["nodes/id"][:])
        }
        if not set(final_ids).issubset(archive_ids):
            raise ValueError("archive no longer contains the candidate's fresh final set")
        final_indices = [archive_ids[node_id] for node_id in final_ids]
        transport_status = _transport_status_report(root, final_indices)
        final_histories = {
            node_id: int(root["nodes/histories"][archive_ids[node_id]])
            for node_id in final_ids
        }
        snapshot_hash = archive_snapshot_sha256(root)

    if not transport_status["passed"]:
        policy = _apply_final_generation_cap(
            {
                "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                "reason": (
                    "failed/max-vertex-capped transport is unusable and cannot "
                    "be repaired by same-node continuation"
                ),
                "retired_node_ids": final_ids,
                "transport_failed_node_ids": transport_status["failed_node_ids"],
                "campaign_command": "refresh-final",
            },
            campaign,
        )
        return {
            "acceptance_passed": False,
            "model_content_sha256": model_hash,
            "archive_snapshot_sha256": snapshot_hash,
            "transport_status": transport_status,
            "policy": policy,
        }

    validation = deepcopy(payload.metadata.get("validation_metadata"))
    if not isinstance(validation, dict):
        raise ValueError("candidate lacks locked validation metadata")
    if validation.get("model_selection_policy") != (
        contribution_aware_model_selection_policy()
    ):
        raise ValueError(
            "candidate lacks the current contribution-aware model-selection policy"
        )
    if validation.get("model_content_sha256") != model_hash:
        raise ValueError("candidate validation metadata names another model")
    if validation.get("archive_snapshot_sha256") != snapshot_hash:
        raise ValueError("archive changed after the candidate was fitted")
    selection = validation.get("selection")
    final = validation.get("final")
    if not isinstance(selection, dict) or not selection or not all(
        isinstance(item, dict) and item.get("validation_passed") is True
        for item in selection.values()
    ):
        raise ValueError("candidate did not pass every model-selection/validation gate")
    if not isinstance(final, dict) or final.get("physics_passed") is not True:
        raise ValueError("candidate did not pass the untouched-final physics gates")
    if (
        final.get("component_acceptance_semantics")
        != _PHYSICAL_COMPONENT_ACCEPTANCE_SEMANTICS
    ):
        raise ValueError(
            "candidate final report lacks physical continuum-component/exact-atom gates"
        )
    locked_transport = final.get("transport_status")
    if not isinstance(locked_transport, Mapping) or locked_transport.get("passed") is not True:
        raise ValueError("candidate final report lacks a passing transport-status gate")
    reported_final_ids = {
        str(item.get("node_id", "")).lower().removeprefix("0x")
        for item in final.get("nodes", ())
        if isinstance(item, dict)
    }
    if reported_final_ids != set(final_ids) or len(final.get("nodes", ())) != 64:
        raise ValueError("candidate validation is not bound to the campaign's fresh final set")
    detector = _external_detector_validation(
        detector_validation_path,
        config_hash=config_hash,
        model_content_sha256=model_hash,
        archive_snapshot_hash=snapshot_hash,
        final_node_ids=final_ids,
    )
    evidence = validate_publication_evidence(
        publication_evidence_path,
        config_sha256=config_hash,
        archive_snapshot_sha256=snapshot_hash,
        model_content_sha256=model_hash,
        final_node_ids=final_ids,
    )
    if not detector["passed"] or not evidence["passed"]:
        if detector.get("supplied") and not detector.get("noise_passed", False):
            noise_ids = detector.get("noise_failing_node_ids", [])
            precision_capped = sorted(
                node_id
                for node_id in noise_ids
                if final_histories.get(node_id, 0)
                >= effective_history_cap("final", BLOCK_HISTORIES)
            )
            if precision_capped:
                policy = {
                    "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                    "reason": (
                        "detector-fold replica noise remains above tolerance at "
                        "the statistical history cap"
                    ),
                    "retired_node_ids": final_ids,
                    "precision_capped_node_ids": precision_capped,
                    "campaign_command": "refresh-final",
                }
            else:
                policy = {
                    "action": "continue_final_histories_for_detector_fold_noise",
                    "node_ids": noise_ids,
                    "campaign_command": (
                        "continue-archive --roles final --detector-validation "
                        "<detector-fold-report.json>"
                    ),
                }
        elif detector.get("supplied") and not detector.get("physics_passed", False):
            policy = {
                "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve",
                "retired_node_ids": final_ids,
                "campaign_command": "refresh-final",
            }
        else:
            policy = {"action": "correct_external_validation_evidence"}
        policy = _apply_final_generation_cap(policy, campaign)
        return {
            "acceptance_passed": False,
            "model_content_sha256": model_hash,
            "archive_snapshot_sha256": snapshot_hash,
            "detector_validation": detector,
            "publication_evidence": evidence,
            "policy": policy,
        }

    final["detector_validation"] = detector
    final.setdefault("gates", {})["detector_fold_validation"] = True
    final["passed"] = True
    validation.update(
        {
            "final": final,
            "production_configuration": production_configuration,
            "publication_evidence": evidence,
            "acceptance_passed": True,
            "model_content_sha256": model_hash,
            "archive_snapshot_sha256": snapshot_hash,
            "payload_artifact": str(output),
            "policy": {"action": "publish_payload_without_refit"},
        }
    )
    payload.metadata["validation_metadata"] = validation
    payload.accepted_for_publication = True
    payload_size = _bind_payload_size_report(
        payload,
        exact_input_candidate_bytes=candidate.stat().st_size,
    )
    validation["runtime_payload_size"] = payload_size
    payload.save_binary(output)
    return validation


_FINITE_FAILURE = 1.0e300


def _replica_noise_metrics(
    numpy: Any,
    replica_s1: Any,
    histories: Any,
    truth: Any,
    *,
    native_shape: tuple[int, int],
    native_grid: bool = False,
) -> tuple[float, float]:
    """Return yield-relative SE and shape-TV noise on the requested grid.

    Allocation and MC eligibility use the fixed, mass-conserving 12x8
    sentinel (or 1x8 for a line).  ``native_grid`` is diagnostic only; the
    native maps are still used for fitting and physical acceptance metrics.
    """

    raw = numpy.asarray(replica_s1, dtype=float).reshape(32, -1)
    counts = numpy.asarray(histories, dtype=float).reshape(32)
    valid = counts > 0.0
    if int(numpy.count_nonzero(valid)) < 2:
        return _FINITE_FAILURE, _FINITE_FAILURE
    estimates = raw[valid] / counts[valid, None]
    yields = numpy.sum(estimates, axis=1)
    truth_array = numpy.asarray(truth, dtype=float).reshape(-1)
    truth_yield = float(numpy.sum(truth_array))
    if truth_yield <= 0.0:
        return (0.0, 0.0) if not numpy.any(raw) else (_FINITE_FAILURE, _FINITE_FAILURE)
    replicas = estimates.shape[0]
    yield_relative_se = float(numpy.std(yields, ddof=1) / math.sqrt(replicas) / truth_yield)
    target_shape = native_shape if native_grid else None
    shape_tv_se = replica_shape_tv_noise(
        numpy,
        raw,
        counts,
        truth_array,
        native_shape=native_shape,
        target_shape=target_shape,
    )
    return yield_relative_se, shape_tv_se


def _stopping_noise_metrics(
    numpy: Any,
    count: int,
    first: Any,
    outer_second: Any,
    total_mean: float,
) -> tuple[float, float]:
    """Return the two distinct stopping-vector uncertainty criteria.

    Entries 0..11 are nonnegative registered-energy bands.  Only bands that
    carry at least one percent of the total activate the tighter absolute
    0.1%-of-total criterion.  Entries 12..15 are signed Legendre moments and
    use the 0.25%-of-total joint-shape sentinel criterion.
    """

    if count < 2 or not math.isfinite(total_mean) or total_mean <= 0.0:
        return _FINITE_FAILURE, _FINITE_FAILURE
    first_array = numpy.asarray(first, dtype=float).reshape(16)
    outer_array = numpy.asarray(outer_second, dtype=float).reshape(16, 16)
    diagonal = numpy.diag(outer_array)
    numerator = numpy.maximum(0.0, diagonal - first_array * first_array / count)
    standard_error = numpy.sqrt(numerator / (count * (count - 1)))
    band_means = first_array[:12] / count
    active_bands = band_means >= 0.01 * total_mean
    band_fraction = (
        float(numpy.max(standard_error[:12][active_bands]) / total_mean)
        if bool(numpy.any(active_bands))
        else 0.0
    )
    legendre_fraction = float(numpy.max(standard_error[12:]) / total_mean)
    if not math.isfinite(band_fraction) or not math.isfinite(legendre_fraction):
        return _FINITE_FAILURE, _FINITE_FAILURE
    return band_fraction, legendre_fraction


def _validation_mc_noise_report(
    root: Any,
    indices: Any,
    descriptors: Sequence[dict[str, Any]],
    offsets: Sequence[int],
) -> dict[str, Any]:
    """Apply the campaign's absolute MC-noise gates to validation truth."""

    _, numpy = require_hdf5()
    n_mu = len(root["axes/mu_edges"]) - 1
    total_n = root["continuum/total/N"]
    n_total_energy = total_n.shape[1] // n_mu
    total_noise_resolution = sentinel_shape(n_total_energy, n_mu)
    total_s1 = root["continuum/total/S1"]
    total_replica = root["continuum/total/replica_S1"]
    replica_histories = root["replicas/histories"][:]
    scattered_n = root["diagnostics/scattered_total_N"][:]
    scattered_s1 = root["diagnostics/scattered_total_S1"][:]
    scattered_s2 = root["diagnostics/scattered_total_S2"][:]
    component_n = root["continuum/components/N"]
    component_s1 = root["continuum/components/S1"]
    component_replica = root["continuum/components/replica_S1"]
    stopping_n = root["stopping/N"][:]
    stopping_s1 = root["stopping/S1"][:]
    stopping_outer = root["stopping/outer_S2"][:]
    line_definitions = _json_attr(root["axes"].attrs["line_definitions_json"])
    n_lines = len(line_definitions)
    line_support = root["lines/support_mask"][:]
    line_n = root["lines/total/N"]
    line_s1 = root["lines/total/S1"]
    line_replica = root["lines/total/replica_S1"]
    direct_fluence = root["nodes/direct_fluence"][:]

    node_reports: list[dict[str, Any]] = []
    for raw_index in indices:
        index = int(raw_index)
        total_n_row = numpy.asarray(total_n[index])
        total_s1_row = numpy.asarray(total_s1[index])
        total_replica_row = numpy.asarray(total_replica[index])
        component_n_row = numpy.asarray(component_n[index])
        component_s1_row = numpy.asarray(component_s1[index])
        component_replica_row = numpy.asarray(component_replica[index])
        line_n_row = numpy.asarray(line_n[index]).reshape(n_lines, n_mu)
        line_s1_row = numpy.asarray(line_s1[index]).reshape(n_lines, n_mu)
        line_replica_row = numpy.asarray(line_replica[index]).reshape(
            32, n_lines, n_mu
        )
        count = int(scattered_n[index])
        total_mean = float(scattered_s1[index] / count) if count else 0.0
        if count > 1 and total_mean > 0.0:
            numerator = max(
                0.0,
                float(scattered_s2[index])
                - float(scattered_s1[index]) ** 2 / count,
            )
            total_relative_se = math.sqrt(numerator / (count * (count - 1))) / total_mean
        else:
            total_relative_se = _FINITE_FAILURE
        band_fraction, legendre_fraction = _stopping_noise_metrics(
            numpy,
            int(stopping_n[index]),
            stopping_s1[index],
            stopping_outer[index],
            total_mean,
        )
        total_truth = numpy.divide(
            total_s1_row,
            total_n_row,
            out=numpy.zeros_like(total_s1_row),
            where=total_n_row > 0,
        )
        _, joint_shape_se = _replica_noise_metrics(
            numpy,
            total_replica_row,
            replica_histories[index],
            total_truth,
            native_shape=(n_total_energy, n_mu),
        )
        _, native_joint_shape_noise = _replica_noise_metrics(
            numpy,
            total_replica_row,
            replica_histories[index],
            total_truth,
            native_shape=(n_total_energy, n_mu),
            native_grid=True,
        )

        component_absolute_yield_noise: list[float] = []
        significant_component_shape_noise: list[float] = []
        component_reports: list[dict[str, Any]] = []
        for component_index, descriptor in enumerate(descriptors):
            start, stop = offsets[component_index], offsets[component_index + 1]
            truth_map = numpy.divide(
                component_s1_row[start:stop],
                component_n_row[start:stop],
                out=numpy.zeros(stop - start),
                where=component_n_row[start:stop] > 0,
            )
            component_mean = float(numpy.sum(truth_map))
            if total_mean <= 0.0 or component_mean <= 0.0:
                component_reports.append(
                    {
                        "kind": "continuum_component",
                        "name": str(descriptor["name"]),
                        "origin": int(descriptor["origin"]),
                        "post_creation_order": int(descriptor["post_creation_order"]),
                        "truth_yield": component_mean,
                        "truth_fraction": 0.0,
                        "yield_noise_gate_active": False,
                        "shape_noise_gate_active": False,
                        "structural_exclusion": _STRUCTURAL_CONTINUUM_EXCLUSIONS.get(
                            str(descriptor["name"])
                        ),
                    }
                )
                continue
            yield_relative_se, shape_se = _replica_noise_metrics(
                numpy,
                component_replica_row[:, start:stop],
                replica_histories[index],
                truth_map,
                native_shape=((stop - start) // n_mu, n_mu),
            )
            _, native_shape_noise = _replica_noise_metrics(
                numpy,
                component_replica_row[:, start:stop],
                replica_histories[index],
                truth_map,
                native_shape=((stop - start) // n_mu, n_mu),
                native_grid=True,
            )
            yield_active = component_mean >= 0.01 * total_mean
            shape_active = component_mean >= 0.05 * total_mean
            absolute_yield_fraction = yield_relative_se * component_mean / total_mean
            if yield_active:
                component_absolute_yield_noise.append(
                    absolute_yield_fraction
                )
            if shape_active:
                significant_component_shape_noise.append(shape_se)
            component_reports.append(
                {
                    "kind": "continuum_component",
                    "name": str(descriptor["name"]),
                    "origin": int(descriptor["origin"]),
                    "post_creation_order": int(descriptor["post_creation_order"]),
                    "truth_yield": component_mean,
                    "truth_fraction": component_mean / total_mean,
                    "yield_noise_gate_active": yield_active,
                    "shape_noise_gate_active": shape_active,
                    "absolute_replica_se_fraction": absolute_yield_fraction,
                    "shape_replica_tv_se": shape_se,
                    "native_shape_replica_tv_noise": native_shape_noise,
                    "shape_noise_resolution": list(
                        sentinel_shape((stop - start) // n_mu, n_mu)
                    ),
                }
            )

        line_noise_ratios: list[float] = []
        line_replica_diagnostics: list[float] = []
        exact_line_reports: list[dict[str, Any]] = []
        for line_index in range(n_lines):
            line_descriptor = line_definitions[line_index]
            line_order = int(line_descriptor["post_creation_order"])
            if line_order != 3:
                raise ValueError(
                    f"exact line {line_descriptor['name']} lacks sentinel order 3"
                )
            if not bool(line_support[index, line_index]):
                exact_line_reports.append(
                    {
                        "kind": "exact_line_atom",
                        "name": str(line_descriptor["name"]),
                        "origin": int(line_descriptor["origin"]),
                        "post_creation_order": line_order,
                        "supported": False,
                        "truth_yield": 0.0,
                        "truth_fraction": 0.0,
                        "yield_noise_gate_active": False,
                        "shape_noise_gate_active": False,
                        "exact_line_noise_gate_active": False,
                        "exact_line_noise_passed": True,
                        "taxonomy_note": (
                            "all-Rayleigh atom; sentinel order 3 is not continuum order0"
                        ),
                    }
                )
                continue
            truth_map = numpy.divide(
                line_s1_row[line_index],
                line_n_row[line_index],
                out=numpy.zeros(n_mu),
                where=line_n_row[line_index] > 0,
            )
            line_yield = float(numpy.sum(truth_map))
            yield_relative_se, shape_se = _replica_noise_metrics(
                numpy,
                line_replica_row[:, line_index, :],
                replica_histories[index],
                truth_map,
                native_shape=(1, n_mu),
            )
            _, native_shape_noise = _replica_noise_metrics(
                numpy,
                line_replica_row[:, line_index, :],
                replica_histories[index],
                truth_map,
                native_shape=(1, n_mu),
                native_grid=True,
            )
            yield_active = total_mean > 0.0 and line_yield >= 0.01 * total_mean
            shape_active = total_mean > 0.0 and line_yield >= 0.05 * total_mean
            absolute_yield_fraction = (
                yield_relative_se * line_yield / total_mean
                if total_mean > 0.0
                else _FINITE_FAILURE
            )
            if yield_active:
                component_absolute_yield_noise.append(absolute_yield_fraction)
            if shape_active:
                significant_component_shape_noise.append(shape_se)
            line_yield_se = yield_relative_se * line_yield
            line_name = str(line_descriptor["name"])
            material, upper_fraction = exact_line_materiality(
                line_yield, line_yield_se, total_mean
            )
            combined_yield = None
            combined_yield_relative_se = None
            combined_shape_se = None
            if line_name == "elastic_E0":
                direct_map = direct_line_map(
                    numpy, float(direct_fluence[index]), n_mu
                )
                combined_truth = truth_map + direct_map
                combined_replicas = (
                    line_replica_row[:, line_index, :]
                    + replica_histories[index, :, None] * direct_map[None, :]
                )
                combined_yield_relative_se, combined_shape_se = _replica_noise_metrics(
                    numpy,
                    combined_replicas,
                    replica_histories[index],
                    combined_truth,
                    native_shape=(1, n_mu),
                )
                combined_yield = float(numpy.sum(combined_truth))
                noise_ratio = max(
                    combined_yield_relative_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                    combined_shape_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                )
                material = True
            else:
                noise_ratio = absolute_yield_fraction / EXACT_LINE_ABSOLUTE_NOISE_FRACTION
                if material:
                    noise_ratio = max(
                        noise_ratio,
                        yield_relative_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                        shape_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                    )
            line_noise_ratios.append(noise_ratio)
            line_replica_diagnostics.append(max(yield_relative_se, shape_se))
            exact_line_reports.append(
                {
                    "kind": "exact_line_atom",
                    "name": line_name,
                    "origin": int(line_descriptor["origin"]),
                    "post_creation_order": line_order,
                    "supported": True,
                    "truth_yield": line_yield,
                    "truth_fraction": (
                        line_yield / total_mean if total_mean > 0.0 else _FINITE_FAILURE
                    ),
                    "yield_noise_gate_active": yield_active,
                    "shape_noise_gate_active": shape_active,
                    "absolute_replica_se_fraction": absolute_yield_fraction,
                    "shape_replica_tv_se": shape_se,
                    "native_shape_replica_tv_noise": native_shape_noise,
                    "shape_noise_resolution": list(sentinel_shape(1, n_mu)),
                    "exact_line_noise_gate_active": True,
                    "line_replica_se": max(yield_relative_se, shape_se),
                    "line_yield_standard_error": line_yield_se,
                    "material_line_gate_active": material,
                    "materiality_upper_fraction_3sigma": upper_fraction,
                    "combined_direct_line_yield": combined_yield,
                    "combined_direct_line_yield_relative_se": combined_yield_relative_se,
                    "combined_direct_line_shape_replica_tv_se": combined_shape_se,
                    "exact_line_noise_ratio": noise_ratio,
                    "exact_line_noise_passed": noise_ratio <= 1.0,
                    "taxonomy_note": (
                        "exact all-Rayleigh atom; sentinel order 3 is not continuum order0"
                    ),
                }
            )

        component_yield_fraction = max(component_absolute_yield_noise, default=0.0)
        component_shape_se = max(significant_component_shape_noise, default=0.0)
        line_noise_ratio = max(line_noise_ratios, default=0.0)
        line_se = max(line_replica_diagnostics, default=0.0)
        passed = bool(
            total_relative_se <= 0.0025
            and band_fraction <= 0.001
            and legendre_fraction <= 0.0025
            and joint_shape_se <= 0.0025
            and component_yield_fraction <= 0.001
            and component_shape_se <= 0.0025
            and line_noise_ratio <= 1.0
        )
        node_reports.append(
            {
                "node_id": f"{int(root['nodes/id'][index]):016x}",
                "passed": passed,
                "total_relative_se": total_relative_se,
                "active_band_absolute_se_fraction": band_fraction,
                "legendre_absolute_se_fraction": legendre_fraction,
                "joint_shape_replica_tv_se": joint_shape_se,
                "native_joint_shape_replica_tv_noise": native_joint_shape_noise,
                "joint_shape_noise_resolution": list(total_noise_resolution),
                "material_component_absolute_replica_se_fraction": component_yield_fraction,
                "significant_component_shape_replica_tv_se": component_shape_se,
                "line_replica_se": line_se,
                "exact_line_noise_ratio": line_noise_ratio,
                "physical_components": [*component_reports, *exact_line_reports],
                "exact_line_atoms": exact_line_reports,
            }
        )
    return {
        "kind": "mc_noise",
        "joint_shape_noise_resolution": list(total_noise_resolution),
        "native_joint_shape_noise_is_diagnostic": True,
        "component_acceptance_semantics": _PHYSICAL_COMPONENT_ACCEPTANCE_SEMANTICS,
        "validation_passed": bool(node_reports) and all(item["passed"] for item in node_reports),
        "nodes": node_reports,
        "detector_fold_noise": (
            "requires external detector-response matrices and is enforced by the "
            "detector-fold validation report on the untouched final set"
        ),
    }


def _canonical_detector_metadata_sha256(value: Any, label: str) -> str:
    if not isinstance(value, Mapping):
        raise ValueError(f"detector response {label} must be a JSON object")
    try:
        canonical = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ValueError(f"detector response {label} is not finite JSON") from error
    return hashlib.sha256(canonical).hexdigest()


def _detector_response_evidence(
    values: Mapping[str, Any],
    detector_name: str,
    *,
    discretization_gate: float,
) -> dict[str, Any]:
    """Validate response provenance and its conservative fold-error bound."""

    response_hash = values.get("response_sha256")
    if not _is_lower_sha256(response_hash):
        raise ValueError(f"{detector_name} detector fold lacks a response-matrix SHA-256")
    if values.get("response_schema") != _DETECTOR_RESPONSE_SCHEMA:
        raise ValueError(f"{detector_name} detector response schema is missing or unsupported")
    if values.get("response_semantics") != _DETECTOR_RESPONSE_SEMANTICS:
        raise ValueError(f"{detector_name} detector response semantics are unsupported")
    if values.get("line_interpolation") != _DETECTOR_LINE_INTERPOLATION:
        raise ValueError(f"{detector_name} exact-line interpolation is unsupported")

    geometry = values.get("geometry")
    geometry_keys = {
        "schema",
        "detector_name",
        "geometry_id",
        "description",
        "active_material",
        "dimensions_cm",
        "mu_definition",
        "source_surface",
    }
    if not isinstance(geometry, Mapping) or set(geometry) != geometry_keys:
        raise ValueError(f"{detector_name} detector geometry evidence is malformed")
    if (
        geometry.get("schema") != _DETECTOR_GEOMETRY_SCHEMA
        or geometry.get("detector_name") != detector_name
        or geometry.get("active_material")
        != {"HPGe": "HPGe", "NaI": "NaI(Tl)"}[detector_name]
        or geometry.get("mu_definition") != _DETECTOR_MU_DEFINITION
    ):
        raise ValueError(f"{detector_name} detector geometry is not representative")
    for name in ("geometry_id", "description", "source_surface"):
        if not isinstance(geometry.get(name), str) or not str(geometry[name]).strip():
            raise ValueError(f"{detector_name} detector geometry {name} is empty")
    if values.get("geometry_id") != geometry["geometry_id"]:
        raise ValueError(f"{detector_name} detector geometry identity does not close")
    dimensions_cm = geometry.get("dimensions_cm")
    if not isinstance(dimensions_cm, Mapping) or not dimensions_cm:
        raise ValueError(f"{detector_name} detector geometry dimensions are missing")
    for name, raw_value in dimensions_cm.items():
        if not isinstance(name, str) or not name or isinstance(raw_value, bool):
            raise ValueError(f"{detector_name} detector geometry dimension is malformed")
        value = float(raw_value)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{detector_name} detector geometry dimensions are invalid")
    geometry_hash = _canonical_detector_metadata_sha256(geometry, "geometry")
    if values.get("geometry_sha256") != geometry_hash:
        raise ValueError(f"{detector_name} detector geometry hash does not close")

    provenance = values.get("provenance")
    provenance_keys = {
        "schema",
        "transport_code",
        "transport_version",
        "physics_configuration",
        "histories_per_grid_point",
        "rng",
        "created_utc",
        "input_sha256",
    }
    if not isinstance(provenance, Mapping) or set(provenance) != provenance_keys:
        raise ValueError(f"{detector_name} detector response provenance is malformed")
    if provenance.get("schema") != _DETECTOR_PROVENANCE_SCHEMA:
        raise ValueError(f"{detector_name} detector response provenance schema is unsupported")
    for name in ("transport_code", "transport_version", "rng", "created_utc"):
        if not isinstance(provenance.get(name), str) or not str(provenance[name]).strip():
            raise ValueError(f"{detector_name} detector response provenance {name} is empty")
    histories = provenance.get("histories_per_grid_point")
    if isinstance(histories, bool) or not isinstance(histories, int) or histories <= 0:
        raise ValueError(f"{detector_name} detector response history count is invalid")
    physics = provenance.get("physics_configuration")
    if not isinstance(physics, Mapping) or not physics:
        raise ValueError(f"{detector_name} detector response physics provenance is empty")
    input_hashes = provenance.get("input_sha256")
    if not isinstance(input_hashes, Mapping) or not input_hashes or any(
        not isinstance(name, str) or not name or not _is_lower_sha256(value)
        for name, value in input_hashes.items()
    ):
        raise ValueError(f"{detector_name} detector response input hashes are invalid")
    provenance_hash = _canonical_detector_metadata_sha256(provenance, "provenance")
    if values.get("provenance_sha256") != provenance_hash:
        raise ValueError(f"{detector_name} detector response provenance hash does not close")

    generator = values.get("generator")
    if not isinstance(generator, Mapping) or set(generator) != {"name", "version", "sha256"}:
        raise ValueError(f"{detector_name} detector response generator is malformed")
    if (
        not isinstance(generator.get("name"), str)
        or not str(generator["name"]).strip()
        or not isinstance(generator.get("version"), str)
        or not str(generator["version"]).strip()
        or not _is_lower_sha256(generator.get("sha256"))
    ):
        raise ValueError(f"{detector_name} detector response generator is invalid")

    dimensions = values.get("dimensions")
    dimension_names = {
        "continuum_energy_bins",
        "mu_bins",
        "spectrum_bins",
        "line_energy_grid_points",
    }
    if not isinstance(dimensions, Mapping) or set(dimensions) != dimension_names or any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in dimensions.values()
    ):
        raise ValueError(f"{detector_name} detector response dimensions are invalid")

    fold = values.get("fold_discretization")
    bound_names = (
        "max_continuum_total_relative_bound",
        "max_continuum_spectral_tv_bound",
        "max_exact_line_relative_bound",
    )
    if not isinstance(fold, Mapping) or set(fold) != {
        "method",
        "gate",
        *bound_names,
    }:
        raise ValueError(f"{detector_name} fold discretization evidence is malformed")
    gate = float(fold.get("gate", math.inf))
    if (
        fold.get("method") != _DETECTOR_DISCRETIZATION_METHOD
        or not math.isfinite(gate)
        or not math.isclose(
            gate, discretization_gate, rel_tol=0.0, abs_tol=1.0e-15
        )
    ):
        raise ValueError(f"{detector_name} fold discretization gate is invalid")
    bounds = {name: float(fold.get(name, math.inf)) for name in bound_names}
    if not all(
        math.isfinite(value) and 0.0 <= value <= discretization_gate
        for value in bounds.values()
    ):
        raise ValueError(f"{detector_name} fold discretization bound exceeds its gate")
    return {
        "response_sha256": response_hash,
        "response_schema": values["response_schema"],
        "response_semantics": values["response_semantics"],
        "line_interpolation": values["line_interpolation"],
        "geometry_id": geometry["geometry_id"],
        "geometry_sha256": geometry_hash,
        "provenance_sha256": provenance_hash,
        "generator": dict(generator),
        "dimensions": dict(dimensions),
        "fold_discretization": {"method": fold["method"], "gate": gate, **bounds},
    }


def _close_detector_discretization_evidence(
    records_by_id: Mapping[str, Mapping[str, Any]],
    ordered_node_ids: Sequence[str],
    evidence: Mapping[str, Any],
    detector_name: str,
) -> None:
    node_to_summary = {
        "continuum_fold_discretization_total_relative_bound": (
            "max_continuum_total_relative_bound"
        ),
        "continuum_fold_discretization_spectral_tv_bound": (
            "max_continuum_spectral_tv_bound"
        ),
        "exact_line_fold_discretization_relative_bound": "max_exact_line_relative_bound",
    }
    fold = evidence["fold_discretization"]
    for node_name, summary_name in node_to_summary.items():
        node_values = [
            float(records_by_id[node_id].get(node_name, math.inf))
            for node_id in ordered_node_ids
        ]
        if not all(math.isfinite(value) and value >= 0.0 for value in node_values):
            raise ValueError(
                f"{detector_name} detector fold lacks finite per-node discretization bounds"
            )
        maximum = max(node_values, default=math.inf)
        if not math.isclose(
            maximum,
            float(fold[summary_name]),
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        ):
            raise ValueError(
                f"{detector_name} detector fold discretization summary does not close"
            )


def _external_detector_validation(
    path: str | Path | None,
    *,
    config_hash: str,
    model_content_sha256: str,
    archive_snapshot_hash: str,
    final_node_ids: Sequence[str],
) -> dict[str, Any]:
    if path is None:
        return {
            "supplied": False,
            "passed": False,
            "physics_passed": False,
            "noise_passed": False,
            "noise_failing_node_ids": [],
            "reason": "HPGe/NaI detector-fold validation report was not supplied",
        }
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    if record.get("schema") != "CEEFIT-DETECTOR-FOLD-VALIDATION-v1":
        raise ValueError("unsupported detector-fold validation schema")
    if record.get("config_sha256") != config_hash:
        raise ValueError("detector-fold validation configuration hash differs")
    if record.get("model_content_sha256") != model_content_sha256:
        raise ValueError("detector-fold validation is bound to a different fitted candidate")
    if record.get("archive_snapshot_sha256") != archive_snapshot_hash:
        raise ValueError("detector-fold validation is bound to a different fit archive snapshot")
    supplied_id_list = [
        str(value).lower().removeprefix("0x") for value in record.get("final_node_ids", ())
    ]
    supplied_ids = set(supplied_id_list)
    if len(supplied_id_list) != len(supplied_ids) or supplied_ids != set(final_node_ids):
        raise ValueError("detector-fold report does not cover exactly the untouched final set")
    detectors = record.get("detectors", {})
    details: dict[str, Any] = {}
    passed = True
    physics_passed = True
    noise_passed = True
    noise_failing_node_ids: set[str] = set()
    for name in ("HPGe", "NaI"):
        values = detectors.get(name)
        if not isinstance(values, dict):
            details[name] = {
                "passed": False,
                "physics_passed": False,
                "noise_passed": False,
                "reason": "missing detector record",
            }
            passed = False
            physics_passed = False
            noise_passed = False
            continue
        response_path = values.get("response_path")
        if not isinstance(response_path, str) or not response_path:
            raise ValueError(f"{name} detector fold lacks its response file path")
        resolved_response = (Path(path).parent / response_path).resolve()
        if not resolved_response.is_file():
            raise ValueError(f"{name} detector response file is missing")
        from .publication import file_sha256

        if file_sha256(resolved_response) != values.get("response_sha256"):
            raise ValueError(f"{name} detector response SHA-256 differs from its file")
        response_evidence = _detector_response_evidence(
            values, name, discretization_gate=0.01
        )
        node_values = values.get("nodes")
        if not isinstance(node_values, list):
            raise ValueError(f"{name} detector fold lacks per-final-node metrics")
        records_by_id: dict[str, dict[str, Any]] = {}
        metric_names = (
            "continuum_total_relative_error",
            "spectral_tv",
            "exact_line_angular_fold_relative_error",
            "replica_relative_se",
        )
        for item in node_values:
            if not isinstance(item, dict):
                raise ValueError(f"{name} detector fold has a malformed node record")
            node_id = str(item.get("node_id", "")).lower().removeprefix("0x")
            if node_id in records_by_id:
                raise ValueError(f"{name} detector fold repeats final node {node_id}")
            records_by_id[node_id] = item
        if set(records_by_id) != supplied_ids:
            raise ValueError(f"{name} detector fold does not score every final node exactly once")
        _close_detector_discretization_evidence(
            records_by_id, supplied_id_list, response_evidence, name
        )
        metrics = {}
        for metric_name in metric_names:
            metric_values = [
                float(records_by_id[node_id].get(metric_name, _FINITE_FAILURE))
                for node_id in supplied_id_list
            ]
            metrics[f"max_{metric_name}"] = max(metric_values, default=_FINITE_FAILURE)
        physics_metric_names = (
            "max_continuum_total_relative_error",
            "max_spectral_tv",
            "max_exact_line_angular_fold_relative_error",
        )
        physics_finite = all(
            math.isfinite(metrics[name]) and metrics[name] >= 0.0
            for name in physics_metric_names
        )
        noise_finite = bool(
            math.isfinite(metrics["max_replica_relative_se"])
            and metrics["max_replica_relative_se"] >= 0.0
        )
        detector_physics_passed = bool(
            physics_finite
            and metrics["max_continuum_total_relative_error"] <= 0.01
            and metrics["max_spectral_tv"] <= 0.01
            and metrics["max_exact_line_angular_fold_relative_error"] <= 0.01
        )
        detector_noise_passed = bool(
            noise_finite and metrics["max_replica_relative_se"] <= 0.0025
        )
        detector_passed = detector_physics_passed and detector_noise_passed
        detector_noise_failures = sorted(
            node_id
            for node_id in supplied_id_list
            if not (
                math.isfinite(
                    float(records_by_id[node_id].get("replica_relative_se", _FINITE_FAILURE))
                )
                and 0.0
                <= float(
                    records_by_id[node_id].get("replica_relative_se", _FINITE_FAILURE)
                )
                <= 0.0025
            )
        )
        noise_failing_node_ids.update(detector_noise_failures)
        details[name] = {
            **response_evidence,
            "node_count": len(records_by_id),
            **metrics,
            "physics_passed": detector_physics_passed,
            "noise_passed": detector_noise_passed,
            "noise_failing_node_ids": detector_noise_failures,
            "passed": detector_passed,
        }
        passed = passed and detector_passed
        physics_passed = physics_passed and detector_physics_passed
        noise_passed = noise_passed and detector_noise_passed
    return {
        "supplied": True,
        "passed": bool(passed),
        "physics_passed": bool(physics_passed),
        "noise_passed": bool(noise_passed),
        "noise_failing_node_ids": sorted(noise_failing_node_ids),
        "detectors": details,
    }


def _score_final_set(
    root: Any,
    payload: RuntimePayload,
    final_indices: Any,
    descriptors: Sequence[dict[str, Any]],
    offsets: Sequence[int],
    *,
    detector_validation_path: str | Path | None,
    config_hash: str,
    model_content_sha256: str,
    archive_snapshot_hash: str,
) -> dict[str, Any]:
    _, numpy = require_hdf5()
    transport_status = _transport_status_report(root, final_indices)
    if not transport_status["passed"]:
        return {
            "evaluated": False,
            "passed": False,
            "physics_passed": False,
            "component_acceptance_semantics": _PHYSICAL_COMPONENT_ACCEPTANCE_SEMANTICS,
            "transport_status": transport_status,
            "gates": {"transport_complete": False},
            "maxima": {},
            # Status metadata is safe to inspect, but residuals from the
            # compromised final set are deliberately never evaluated.
            "nodes": transport_status["nodes"],
            "detector_validation": {
                "supplied": False,
                "passed": False,
                "physics_passed": False,
                "noise_passed": False,
                "reason": "detector folds are not evaluated for invalid transport truth",
            },
            "reason": transport_status["reason"],
        }
    q = numpy.asarray(root["nodes/q"][:], dtype=float)
    coordinates = root["nodes/coordinates"][:]
    mu_edges = numpy.asarray(root["axes/mu_edges"][:], dtype=float)
    n_mu = len(mu_edges) - 1
    total_n = root["continuum/total/N"]
    total_s1 = root["continuum/total/S1"]
    total_replica = root["continuum/total/replica_S1"]
    n_total_energy = total_n.shape[1] // n_mu
    total_noise_resolution = sentinel_shape(n_total_energy, n_mu)
    replica_histories = root["replicas/histories"][:]
    line_definitions = _json_attr(root["axes"].attrs["line_definitions_json"])
    line_names = [item["name"] for item in line_definitions]
    line_support = root["lines/support_mask"][:]
    line_n = root["lines/total/N"]
    line_s1 = root["lines/total/S1"]
    line_replica = root["lines/total/replica_S1"]
    scattered_n = root["diagnostics/scattered_total_N"][:]
    scattered_s1 = root["diagnostics/scattered_total_S1"][:]
    scattered_s2 = root["diagnostics/scattered_total_S2"][:]
    component_n = root["continuum/components/N"]
    component_s1 = root["continuum/components/S1"]
    component_replica = root["continuum/components/replica_S1"]
    stopping_n = root["stopping/N"][:]
    stopping_s1 = root["stopping/S1"][:]
    stopping_outer = root["stopping/outer_S2"][:]
    direct = root["nodes/direct_fluence"][:]
    histories = root["nodes/histories"][:]
    statuses = list(root["nodes/status"].asstr()[:])
    failures = root["nodes/transport_failures"][:]
    final_history_cap = effective_history_cap("final", BLOCK_HISTORIES)

    metrics: dict[str, list[float]] = {
        key: []
        for key in (
            "integrated_error",
            "fit_only_proxy_error",
            "joint_tv",
            "energy_marginal_tv",
            "mu_marginal_tv",
            "significant_component_yield_error",
            "significant_component_joint_tv",
            "line_yield_error",
            "line_angular_tv",
            "material_line_yield_relative_error",
            "material_line_angular_tv",
            "weak_line_absolute_yield_error_fraction",
            "weak_line_absolute_angular_l1_fraction",
            "exact_line_noise_ratio",
            "total_relative_mc_se",
            "active_band_absolute_mc_se_fraction",
            "legendre_absolute_mc_se_fraction",
            "joint_shape_replica_tv_se",
            "native_joint_shape_replica_tv_noise",
            "material_component_absolute_replica_se_fraction",
            "significant_component_shape_replica_tv_se",
            "significant_component_replica_se",
            "line_replica_se",
            "support_leakage_fraction",
            "most_negative_mass",
        )
    }
    node_reports = []
    final_node_ids = [f"{int(root['nodes/id'][index]):016x}" for index in final_indices]
    for raw_index in final_indices:
        index = int(raw_index)
        total_n_row = numpy.asarray(total_n[index])
        total_s1_row = numpy.asarray(total_s1[index])
        total_replica_row = numpy.asarray(total_replica[index])
        component_n_row = numpy.asarray(component_n[index])
        component_s1_row = numpy.asarray(component_s1[index])
        component_replica_row = numpy.asarray(component_replica[index])
        line_n_row = numpy.asarray(line_n[index]).reshape(len(line_names), n_mu)
        line_s1_row = numpy.asarray(line_s1[index]).reshape(len(line_names), n_mu)
        line_replica_row = numpy.asarray(line_replica[index]).reshape(
            32, len(line_names), n_mu
        )
        energy, distance_cm, hsrc_cm, hdet_cm, _ = map(float, coordinates[index])
        physical_edges = numpy.asarray(
            root[
                f"axes/energy/{int(root['axes'].attrs['total_energy_axis'])}/physical_edges_keV"
            ][index],
            dtype=float,
        )
        evaluated = payload.evaluate(
            energy,
            distance_cm * 0.01,
            hsrc_cm * 0.01,
            hdet_cm * 0.01,
            physical_edges,
            mu_edges,
            direct_fluence_cm2=float(direct[index]),
        )
        predicted_continuum = numpy.asarray(evaluated["continuum"], dtype=float)
        truth_continuum = numpy.divide(
            total_s1_row,
            total_n_row,
            out=numpy.zeros_like(total_s1_row),
            where=total_n_row > 0,
        ).reshape(predicted_continuum.shape)
        pred_cont_yield = float(numpy.sum(predicted_continuum))
        truth_cont_yield = float(numpy.sum(truth_continuum))
        evaluated_lines = {item["name"]: item for item in evaluated["lines"]}
        truth_line_maps = numpy.divide(
            line_s1_row,
            line_n_row,
            out=numpy.zeros_like(line_s1_row),
            where=line_n_row > 0,
        )
        predicted_line_yield = float(
            sum(sum(item["mu_bin_masses"]) for item in evaluated_lines.values())
        )
        truth_line_yield = float(numpy.sum(truth_line_maps))
        truth_total = float(scattered_s1[index] / scattered_n[index]) if scattered_n[index] else 0.0
        if truth_total <= 0.0:
            total_error = _FINITE_FAILURE
        else:
            total_error = abs((pred_cont_yield + predicted_line_yield) / truth_total - 1.0)
        metrics["integrated_error"].append(total_error)
        metrics["fit_only_proxy_error"].append(total_error)

        if pred_cont_yield > 0.0 and truth_cont_yield > 0.0:
            pred_shape = predicted_continuum / pred_cont_yield
            truth_shape = truth_continuum / truth_cont_yield
            jtv = 0.5 * float(numpy.sum(numpy.abs(pred_shape - truth_shape)))
            etv = 0.5 * float(
                numpy.sum(
                    numpy.abs(
                        numpy.sum(pred_shape, axis=1) - numpy.sum(truth_shape, axis=1)
                    )
                )
            )
            mtv = 0.5 * float(
                numpy.sum(
                    numpy.abs(
                        numpy.sum(pred_shape, axis=0) - numpy.sum(truth_shape, axis=0)
                    )
                )
            )
        else:
            jtv = etv = mtv = _FINITE_FAILURE
        metrics["joint_tv"].append(jtv)
        metrics["energy_marginal_tv"].append(etv)
        metrics["mu_marginal_tv"].append(mtv)

        total_count = int(scattered_n[index])
        if total_count > 1 and truth_total > 0.0:
            total_variance = max(
                0.0,
                float(scattered_s2[index])
                - float(scattered_s1[index]) ** 2 / total_count,
            ) / (total_count * (total_count - 1))
            total_relative_se = math.sqrt(total_variance) / truth_total
        else:
            total_relative_se = _FINITE_FAILURE
        metrics["total_relative_mc_se"].append(total_relative_se)
        band_noise_fraction, legendre_noise_fraction = _stopping_noise_metrics(
            numpy,
            int(stopping_n[index]),
            stopping_s1[index],
            stopping_outer[index],
            truth_total,
        )
        metrics["active_band_absolute_mc_se_fraction"].append(band_noise_fraction)
        metrics["legendre_absolute_mc_se_fraction"].append(legendre_noise_fraction)
        _, total_shape_replica_se = _replica_noise_metrics(
            numpy,
            total_replica_row,
            replica_histories[index],
            truth_continuum,
            native_shape=(n_total_energy, n_mu),
        )
        _, native_total_shape_replica_noise = _replica_noise_metrics(
            numpy,
            total_replica_row,
            replica_histories[index],
            truth_continuum,
            native_shape=(n_total_energy, n_mu),
            native_grid=True,
        )
        metrics["joint_shape_replica_tv_se"].append(total_shape_replica_se)
        metrics["native_joint_shape_replica_tv_noise"].append(
            native_total_shape_replica_noise
        )

        node_component_yield_errors: list[float] = []
        node_component_shape_errors: list[float] = []
        node_component_absolute_yield_noise: list[float] = []
        node_component_shape_noise: list[float] = []
        node_component_noise: list[float] = []
        node_component_reports: list[dict[str, Any]] = []
        for component_index, descriptor in enumerate(descriptors):
            start, stop = offsets[component_index], offsets[component_index + 1]
            truth_map = numpy.divide(
                component_s1_row[start:stop],
                component_n_row[start:stop],
                out=numpy.zeros(stop - start),
                where=component_n_row[start:stop] > 0,
            )
            truth_yield = float(numpy.sum(truth_map))
            if truth_total <= 0.0 or truth_yield <= 0.0:
                node_component_reports.append(
                    {
                        "kind": "continuum_component",
                        "name": str(descriptor["name"]),
                        "origin": int(descriptor["origin"]),
                        "post_creation_order": int(descriptor["post_creation_order"]),
                        "truth_yield": truth_yield,
                        "truth_fraction": 0.0,
                        "yield_gate_active": False,
                        "shape_gate_active": False,
                        "structural_exclusion": _STRUCTURAL_CONTINUUM_EXCLUSIONS.get(
                            str(descriptor["name"])
                        ),
                    }
                )
                continue
            replica_yield_se, replica_shape_se = _replica_noise_metrics(
                numpy,
                component_replica_row[:, start:stop],
                replica_histories[index],
                truth_map,
                native_shape=((stop - start) // n_mu, n_mu),
            )
            _, native_replica_shape_noise = _replica_noise_metrics(
                numpy,
                component_replica_row[:, start:stop],
                replica_histories[index],
                truth_map,
                native_shape=((stop - start) // n_mu, n_mu),
                native_grid=True,
            )
            truth_fraction = (
                truth_yield / truth_total if truth_total > 0.0 else _FINITE_FAILURE
            )
            yield_noise_active = (
                truth_total > 0.0 and truth_yield >= 0.01 * truth_total
            )
            acceptance_active = (
                truth_total > 0.0 and truth_yield >= 0.05 * truth_total
            )
            absolute_yield_noise = replica_yield_se * truth_fraction
            if yield_noise_active:
                node_component_absolute_yield_noise.append(
                    absolute_yield_noise
                )

            prediction = numpy.zeros(stop - start, dtype=float)
            model = payload.components.get(descriptor["name"])
            if model is not None:
                axis = str(descriptor["energy_axis"])
                registered_edges = root[f"axes/energy/{axis}/registered_edges"][:]
                prediction = (
                    model.predict_cell_masses(
                        q[index : index + 1], registered_edges, mu_edges
                    )[2][0]
                    * float(direct[index])
                )
            predicted_yield = float(numpy.sum(prediction))
            yield_error = (
                abs(predicted_yield / truth_yield - 1.0)
                if predicted_yield >= 0.0 and truth_yield > 0.0
                else _FINITE_FAILURE
            )
            shape_error = (
                0.5
                * float(
                    numpy.sum(
                        numpy.abs(prediction / predicted_yield - truth_map / truth_yield)
                    )
                )
                if predicted_yield > 0.0
                else _FINITE_FAILURE
            )
            if acceptance_active:
                node_component_yield_errors.append(yield_error)
                node_component_shape_errors.append(shape_error)
                node_component_shape_noise.append(replica_shape_se)
                node_component_noise.append(max(replica_yield_se, replica_shape_se))
            node_component_reports.append(
                {
                    "kind": "continuum_component",
                    "name": str(descriptor["name"]),
                    "origin": int(descriptor["origin"]),
                    "post_creation_order": int(descriptor["post_creation_order"]),
                    "truth_yield": truth_yield,
                    "truth_fraction": truth_fraction,
                    "predicted_yield": predicted_yield,
                    "yield_noise_gate_active": yield_noise_active,
                    "yield_gate_active": acceptance_active,
                    "shape_gate_active": acceptance_active,
                    "yield_error": yield_error,
                    "joint_tv": shape_error,
                    "absolute_replica_se_fraction": absolute_yield_noise,
                    "shape_replica_tv_se": replica_shape_se,
                    "native_shape_replica_tv_noise": native_replica_shape_noise,
                    "shape_noise_resolution": list(
                        sentinel_shape((stop - start) // n_mu, n_mu)
                    ),
                }
            )
        node_line_yield_errors: list[float] = []
        node_line_shape_errors: list[float] = []
        node_material_line_yield_errors: list[float] = []
        node_material_line_shape_errors: list[float] = []
        node_weak_line_absolute_yield_errors: list[float] = []
        node_weak_line_absolute_angular_errors: list[float] = []
        node_line_noise_ratios: list[float] = []
        node_line_noise_diagnostics: list[float] = []
        node_exact_line_reports: list[dict[str, Any]] = []
        leakage = 0.0
        exact_atom_energies = True
        for line_index, line_name in enumerate(line_names):
            line_descriptor = line_definitions[line_index]
            line_order = int(line_descriptor["post_creation_order"])
            if line_order != 3:
                raise ValueError(f"exact line {line_name} lacks sentinel order 3")
            supported = bool(line_support[index, line_index])
            predicted_record = evaluated_lines.get(line_name)
            truth_map = truth_line_maps[line_index]
            truth_yield = float(numpy.sum(truth_map))
            if not supported:
                if predicted_record is not None:
                    leakage += abs(float(sum(predicted_record["mu_bin_masses"])))
                node_exact_line_reports.append(
                    {
                        "kind": "exact_line_atom",
                        "name": line_name,
                        "origin": int(line_descriptor["origin"]),
                        "post_creation_order": line_order,
                        "supported": False,
                        "truth_yield": truth_yield,
                        "truth_fraction": 0.0,
                        "yield_gate_active": False,
                        "shape_gate_active": False,
                        "exact_line_yield_gate_active": False,
                        "exact_line_angular_gate_active": False,
                        "exact_line_noise_gate_active": False,
                        "taxonomy_note": (
                            "all-Rayleigh atom; sentinel order 3 is not continuum order0"
                        ),
                    }
                )
                continue
            if predicted_record is None:
                node_material_line_yield_errors.append(_FINITE_FAILURE)
                node_material_line_shape_errors.append(_FINITE_FAILURE)
                node_line_noise_ratios.append(_FINITE_FAILURE)
                node_exact_line_reports.append(
                    {
                        "kind": "exact_line_atom",
                        "name": line_name,
                        "origin": int(line_descriptor["origin"]),
                        "post_creation_order": line_order,
                        "supported": True,
                        "truth_yield": truth_yield,
                        "truth_fraction": (
                            truth_yield / truth_total if truth_total > 0.0 else 0.0
                        ),
                        "yield_gate_active": (
                            truth_total > 0.0 and truth_yield >= 0.05 * truth_total
                        ),
                        "shape_gate_active": (
                            truth_total > 0.0 and truth_yield >= 0.05 * truth_total
                        ),
                        "standalone_yield_error": _FINITE_FAILURE,
                        "standalone_joint_tv": _FINITE_FAILURE,
                        "exact_line_yield_gate_active": True,
                        "exact_line_angular_gate_active": True,
                        "exact_line_noise_gate_active": True,
                        "exact_line_noise_ratio": _FINITE_FAILURE,
                        "exact_line_gates_passed": False,
                        "taxonomy_note": (
                            "all-Rayleigh atom; sentinel order 3 is not continuum order0"
                        ),
                    }
                )
                continue
            prediction = numpy.asarray(predicted_record["mu_bin_masses"], dtype=float)
            predicted_yield = float(numpy.sum(prediction))
            yield_error = (
                abs(predicted_yield / truth_yield - 1.0)
                if truth_yield > 0.0
                else (0.0 if predicted_yield == 0.0 else _FINITE_FAILURE)
            )
            shape_error = (
                0.5
                * float(
                    numpy.sum(
                        numpy.abs(prediction / predicted_yield - truth_map / truth_yield)
                    )
                )
                if predicted_yield > 0.0 and truth_yield > 0.0
                else (0.0 if predicted_yield == truth_yield == 0.0 else _FINITE_FAILURE)
            )
            node_line_yield_errors.append(yield_error)
            node_line_shape_errors.append(shape_error)
            replica_yield_se, replica_shape_se = _replica_noise_metrics(
                numpy,
                line_replica_row[:, line_index, :],
                replica_histories[index],
                truth_map,
                native_shape=(1, n_mu),
            )
            _, native_replica_shape_noise = _replica_noise_metrics(
                numpy,
                line_replica_row[:, line_index, :],
                replica_histories[index],
                truth_map,
                native_shape=(1, n_mu),
                native_grid=True,
            )
            line_replica_value = max(replica_yield_se, replica_shape_se)
            node_line_noise_diagnostics.append(line_replica_value)
            truth_fraction = (
                truth_yield / truth_total if truth_total > 0.0 else _FINITE_FAILURE
            )
            yield_noise_active = (
                truth_total > 0.0 and truth_yield >= 0.01 * truth_total
            )
            acceptance_active = (
                truth_total > 0.0 and truth_yield >= 0.05 * truth_total
            )
            absolute_yield_noise = replica_yield_se * truth_fraction
            if yield_noise_active:
                node_component_absolute_yield_noise.append(absolute_yield_noise)
            if acceptance_active:
                node_component_yield_errors.append(yield_error)
                node_component_shape_errors.append(shape_error)
                node_component_shape_noise.append(replica_shape_se)
                node_component_noise.append(max(replica_yield_se, replica_shape_se))
            line_yield_se = replica_yield_se * truth_yield
            material, upper_fraction = exact_line_materiality(
                truth_yield, line_yield_se, truth_total
            )
            combined_truth_yield = None
            combined_predicted_yield = None
            combined_yield_error = None
            combined_shape_error = None
            combined_replica_yield_se = None
            combined_replica_shape_se = None
            weak_absolute_yield_error = None
            weak_absolute_angular_error = None
            if line_name == "elastic_E0":
                direct_map = direct_line_map(numpy, float(direct[index]), n_mu)
                combined_truth = truth_map + direct_map
                combined_prediction = prediction + direct_map
                combined_truth_yield = float(numpy.sum(combined_truth))
                combined_predicted_yield = float(numpy.sum(combined_prediction))
                combined_yield_error = abs(
                    combined_predicted_yield / combined_truth_yield - 1.0
                )
                combined_shape_error = 0.5 * float(
                    numpy.sum(
                        numpy.abs(
                            combined_prediction / combined_predicted_yield
                            - combined_truth / combined_truth_yield
                        )
                    )
                )
                combined_replicas = (
                    line_replica_row[:, line_index, :]
                    + replica_histories[index, :, None] * direct_map[None, :]
                )
                (
                    combined_replica_yield_se,
                    combined_replica_shape_se,
                ) = _replica_noise_metrics(
                    numpy,
                    combined_replicas,
                    replica_histories[index],
                    combined_truth,
                    native_shape=(1, n_mu),
                )
                policy_yield_error = combined_yield_error
                policy_shape_error = combined_shape_error
                node_material_line_yield_errors.append(policy_yield_error)
                node_material_line_shape_errors.append(policy_shape_error)
                noise_ratio = max(
                    combined_replica_yield_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                    combined_replica_shape_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                )
                material = True
            elif material:
                policy_yield_error = yield_error
                policy_shape_error = shape_error
                node_material_line_yield_errors.append(policy_yield_error)
                node_material_line_shape_errors.append(policy_shape_error)
                noise_ratio = max(
                    absolute_yield_noise / EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
                    replica_yield_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                    replica_shape_se / EXACT_LINE_RELATIVE_NOISE_FRACTION,
                )
            else:
                weak_absolute_yield_error = (
                    abs(predicted_yield - truth_yield) / truth_total
                )
                weak_absolute_angular_error = (
                    0.5 * float(numpy.sum(numpy.abs(prediction - truth_map)))
                    / truth_total
                )
                node_weak_line_absolute_yield_errors.append(
                    weak_absolute_yield_error
                )
                node_weak_line_absolute_angular_errors.append(
                    weak_absolute_angular_error
                )
                policy_yield_error = weak_absolute_yield_error
                policy_shape_error = weak_absolute_angular_error
                noise_ratio = max(
                    absolute_yield_noise / EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
                    replica_shape_se * truth_fraction
                    / EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
                )
            node_line_noise_ratios.append(noise_ratio)
            physics_passed = bool(
                (
                    policy_yield_error <= EXACT_LINE_RELATIVE_ERROR_FRACTION
                    and policy_shape_error <= EXACT_LINE_RELATIVE_ERROR_FRACTION
                )
                if material
                else (
                    policy_yield_error <= EXACT_LINE_ABSOLUTE_ERROR_FRACTION
                    and policy_shape_error <= EXACT_LINE_ABSOLUTE_ERROR_FRACTION
                )
            )
            node_exact_line_reports.append(
                {
                    "kind": "exact_line_atom",
                    "name": line_name,
                    "origin": int(line_descriptor["origin"]),
                    "post_creation_order": line_order,
                    "supported": True,
                    "truth_yield": truth_yield,
                    "truth_fraction": truth_fraction,
                    "predicted_yield": predicted_yield,
                    "yield_noise_gate_active": yield_noise_active,
                    "yield_gate_active": acceptance_active,
                    "shape_gate_active": acceptance_active,
                    "standalone_yield_error": yield_error,
                    "standalone_joint_tv": shape_error,
                    "material_line_gate_active": material,
                    "materiality_upper_fraction_3sigma": upper_fraction,
                    "policy_yield_error": policy_yield_error,
                    "policy_angular_error": policy_shape_error,
                    "weak_absolute_yield_error_fraction": weak_absolute_yield_error,
                    "weak_absolute_angular_l1_fraction": weak_absolute_angular_error,
                    "combined_direct_line_truth_yield": combined_truth_yield,
                    "combined_direct_line_predicted_yield": combined_predicted_yield,
                    "combined_direct_line_yield_error": combined_yield_error,
                    "combined_direct_line_angular_tv": combined_shape_error,
                    "exact_line_yield_gate_active": True,
                    "exact_line_angular_gate_active": True,
                    "exact_line_noise_gate_active": True,
                    "line_replica_se": line_replica_value,
                    "combined_direct_line_yield_relative_se": combined_replica_yield_se,
                    "combined_direct_line_shape_replica_tv_se": combined_replica_shape_se,
                    "exact_line_noise_ratio": noise_ratio,
                    "exact_line_gates_passed": physics_passed and noise_ratio <= 1.0,
                    "absolute_replica_se_fraction": absolute_yield_noise,
                    "shape_replica_tv_se": replica_shape_se,
                    "native_shape_replica_tv_noise": native_replica_shape_noise,
                    "shape_noise_resolution": list(sentinel_shape(1, n_mu)),
                    "taxonomy_note": (
                        "all-Rayleigh atom; sentinel order 3 is not continuum order0"
                    ),
                }
            )
            expected_energy = energy if line_name == "elastic_E0" else 510.998950
            exact_atom_energies = exact_atom_energies and math.isclose(
                float(predicted_record["energy_keV"]), expected_energy, rel_tol=0.0, abs_tol=1.0e-12
            )
        metrics["line_yield_error"].append(max(node_line_yield_errors, default=0.0))
        metrics["line_angular_tv"].append(max(node_line_shape_errors, default=0.0))
        metrics["material_line_yield_relative_error"].append(
            max(node_material_line_yield_errors, default=0.0)
        )
        metrics["material_line_angular_tv"].append(
            max(node_material_line_shape_errors, default=0.0)
        )
        metrics["weak_line_absolute_yield_error_fraction"].append(
            max(node_weak_line_absolute_yield_errors, default=0.0)
        )
        metrics["weak_line_absolute_angular_l1_fraction"].append(
            max(node_weak_line_absolute_angular_errors, default=0.0)
        )
        metrics["exact_line_noise_ratio"].append(
            max(node_line_noise_ratios, default=0.0)
        )
        metrics["line_replica_se"].append(
            max(node_line_noise_diagnostics, default=0.0)
        )
        metrics["significant_component_yield_error"].append(
            max(node_component_yield_errors, default=0.0)
        )
        metrics["significant_component_joint_tv"].append(
            max(node_component_shape_errors, default=0.0)
        )
        metrics["significant_component_replica_se"].append(
            max(node_component_noise, default=0.0)
        )
        metrics["material_component_absolute_replica_se_fraction"].append(
            max(node_component_absolute_yield_noise, default=0.0)
        )
        metrics["significant_component_shape_replica_tv_se"].append(
            max(node_component_shape_noise, default=0.0)
        )

        all_predicted = [predicted_continuum.reshape(-1)]
        all_predicted.extend(
            numpy.asarray(value, dtype=float).reshape(-1)
            for value in evaluated["components"].values()
        )
        all_predicted.extend(
            numpy.asarray(value["mu_bin_masses"], dtype=float).reshape(-1)
            for value in evaluated_lines.values()
        )
        finite_nonnegative = all(
            numpy.all(numpy.isfinite(array)) and float(numpy.min(array, initial=0.0)) >= -1.0e-14
            for array in all_predicted
        )
        most_negative = min(
            (float(numpy.min(array, initial=0.0)) for array in all_predicted), default=0.0
        )
        metrics["most_negative_mass"].append(most_negative if finite_nonnegative else -_FINITE_FAILURE)
        probe_high = 1.05 * max(energy, 510.998950) + 1.0
        for name, model in payload.components.items():
            endpoint = 510.998950 if model.axis_kind == "annihilation" else energy
            probe_edges = tuple(sorted({0.0, model.e_cut_keV, endpoint, probe_high}))
            probe = payload.evaluate(
                energy,
                distance_cm * 0.01,
                hsrc_cm * 0.01,
                hdet_cm * 0.01,
                probe_edges,
                mu_edges,
                direct_fluence_cm2=float(direct[index]),
            )
            cells = numpy.asarray(probe["components"][name], dtype=float)
            leakage += float(numpy.sum(numpy.abs(cells[0])))
            leakage += float(numpy.sum(numpy.abs(cells[-1])))
            if model.pair_threshold_keV is not None and energy < model.pair_threshold_keV:
                leakage += float(numpy.sum(numpy.abs(cells)))
        if not exact_atom_energies:
            leakage = _FINITE_FAILURE
        leakage_fraction = leakage / max(truth_total, 1.0e-300)
        metrics["support_leakage_fraction"].append(min(leakage_fraction, _FINITE_FAILURE))

        status_eligible = statuses[index] in {"complete", "converged"} and int(failures[index]) == 0
        noise_eligible = bool(
            status_eligible
            and total_relative_se <= 0.0025
            and band_noise_fraction <= 0.001
            and legendre_noise_fraction <= 0.0025
            and total_shape_replica_se <= 0.0025
            and max(node_component_absolute_yield_noise, default=0.0) <= 0.001
            and max(node_component_shape_noise, default=0.0) <= 0.0025
            and max(node_line_noise_ratios, default=0.0) <= 1.0
        )
        node_reports.append(
            {
                "node_id": f"{int(root['nodes/id'][index]):016x}",
                "status": statuses[index],
                "transport_eligible": status_eligible,
                "histories": int(histories[index]),
                "history_budget_exhausted": int(histories[index]) >= final_history_cap,
                "integrated_error": total_error,
                "joint_tv": jtv,
                "energy_marginal_tv": etv,
                "mu_marginal_tv": mtv,
                "significant_component_yield_error": max(node_component_yield_errors, default=0.0),
                "significant_component_joint_tv": max(node_component_shape_errors, default=0.0),
                "line_yield_error": max(node_line_yield_errors, default=0.0),
                "line_angular_tv": max(node_line_shape_errors, default=0.0),
                "material_line_yield_relative_error": max(
                    node_material_line_yield_errors, default=0.0
                ),
                "material_line_angular_tv": max(
                    node_material_line_shape_errors, default=0.0
                ),
                "weak_line_absolute_yield_error_fraction": max(
                    node_weak_line_absolute_yield_errors, default=0.0
                ),
                "weak_line_absolute_angular_l1_fraction": max(
                    node_weak_line_absolute_angular_errors, default=0.0
                ),
                "exact_line_noise_ratio": max(node_line_noise_ratios, default=0.0),
                "total_relative_mc_se": total_relative_se,
                "active_band_absolute_mc_se_fraction": band_noise_fraction,
                "legendre_absolute_mc_se_fraction": legendre_noise_fraction,
                "joint_shape_replica_tv_se": total_shape_replica_se,
                "native_joint_shape_replica_tv_noise": native_total_shape_replica_noise,
                "joint_shape_noise_resolution": list(total_noise_resolution),
                "material_component_absolute_replica_se_fraction": max(
                    node_component_absolute_yield_noise, default=0.0
                ),
                "significant_component_shape_replica_tv_se": max(
                    node_component_shape_noise, default=0.0
                ),
                "physical_components": [
                    *node_component_reports,
                    *node_exact_line_reports,
                ],
                "exact_line_atoms": node_exact_line_reports,
                "mc_noise_eligible": noise_eligible,
                "support_leakage_fraction": min(leakage_fraction, _FINITE_FAILURE),
                "nonnegative_finite": finite_nonnegative,
            }
        )

    maxima = {
        key: max(values, default=_FINITE_FAILURE)
        for key, values in metrics.items()
        if key != "most_negative_mass"
    }
    maxima["most_negative_mass"] = min(metrics["most_negative_mass"], default=-_FINITE_FAILURE)
    detector = _external_detector_validation(
        detector_validation_path,
        config_hash=config_hash,
        model_content_sha256=model_content_sha256,
        archive_snapshot_hash=archive_snapshot_hash,
        final_node_ids=final_node_ids,
    )
    physics_gates = {
        "transport_complete": bool(transport_status["passed"]),
        "integrated_le_1pct": maxima["integrated_error"] <= 0.01,
        # A single noisy MC truth realization cannot identify the fit-only
        # contribution.  This conservative residual contains both fit and MC
        # error; the separate null-calibration publication check quantifies the
        # latter before an accepted payload can be written.
        "fit_plus_mc_residual_target_le_0p5pct": maxima["fit_only_proxy_error"] <= 0.005,
        "joint_tv_le_1pct": maxima["joint_tv"] <= 0.01,
        "energy_marginal_tv_le_0p5pct": maxima["energy_marginal_tv"] <= 0.005,
        "mu_marginal_tv_le_0p5pct": maxima["mu_marginal_tv"] <= 0.005,
        "significant_component_yield_le_1pct": maxima["significant_component_yield_error"] <= 0.01,
        "significant_component_shape_le_1pct": maxima["significant_component_joint_tv"] <= 0.01,
        "material_exact_line_yield_le_1pct": (
            maxima["material_line_yield_relative_error"]
            <= EXACT_LINE_RELATIVE_ERROR_FRACTION
        ),
        "material_exact_line_angular_shape_le_1pct": (
            maxima["material_line_angular_tv"]
            <= EXACT_LINE_RELATIVE_ERROR_FRACTION
        ),
        "weak_exact_line_absolute_yield_le_0p1pct_total": (
            maxima["weak_line_absolute_yield_error_fraction"]
            <= EXACT_LINE_ABSOLUTE_ERROR_FRACTION
        ),
        "weak_exact_line_absolute_angular_l1_le_0p1pct_total": (
            maxima["weak_line_absolute_angular_l1_fraction"]
            <= EXACT_LINE_ABSOLUTE_ERROR_FRACTION
        ),
        "mc_noise_eligible": all(item["mc_noise_eligible"] for item in node_reports),
        "no_support_leakage": maxima["support_leakage_fraction"] <= 1.0e-12,
        "nonnegative_finite": maxima["most_negative_mass"] >= -1.0e-14,
    }
    gates = {
        **physics_gates,
        "detector_fold_validation": bool(detector["passed"]),
    }
    return {
        "evaluated": True,
        "passed": bool(node_reports) and all(gates.values()),
        "physics_passed": bool(node_reports) and all(physics_gates.values()),
        "joint_shape_noise_resolution": list(total_noise_resolution),
        "native_joint_shape_noise_is_diagnostic": True,
        "component_acceptance_semantics": _PHYSICAL_COMPONENT_ACCEPTANCE_SEMANTICS,
        "transport_status": transport_status,
        "gates": gates,
        "maxima": maxima,
        "nodes": node_reports,
        "detector_validation": detector,
    }


def _mu_detector_validation(
    path: str | Path | None,
    *,
    config_sha256: str,
    archive_snapshot_hash: str,
    pilot_node_ids: Sequence[str],
    tolerance: float,
) -> dict[str, Any]:
    if path is None:
        return {
            "supplied": False,
            "passed": False,
            "reason": "HPGe/NaI 64-to-32 fold comparison was not supplied",
        }
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    if record.get("schema") != "CEEFIT-MU-RESOLUTION-VALIDATION-v1":
        raise ValueError("unsupported mu-resolution detector report schema")
    if record.get("config_sha256") != config_sha256:
        raise ValueError("mu-resolution report configuration hash differs")
    if record.get("archive_snapshot_sha256") != archive_snapshot_hash:
        raise ValueError("mu-resolution report is bound to another archive snapshot")
    if record.get("source_bins") != 64 or record.get("candidate_bins") != 32:
        raise ValueError("mu-resolution report does not compare 64 to 32 bins")
    report_tolerance = float(record.get("tolerance", math.inf))
    if (
        not math.isfinite(report_tolerance)
        or not 0.0 < report_tolerance <= 0.002
        or not math.isclose(report_tolerance, tolerance, rel_tol=0.0, abs_tol=1.0e-15)
    ):
        raise ValueError("mu-resolution detector report tolerance differs or exceeds 0.2%")
    supplied = [str(value).lower().removeprefix("0x") for value in record.get("pilot_node_ids", ())]
    if len(supplied) != len(set(supplied)) or set(supplied) != set(pilot_node_ids):
        raise ValueError("mu-resolution report does not cover every eligible pilot node")
    details: dict[str, Any] = {}
    passed = True
    for detector_name in ("HPGe", "NaI"):
        detector = record.get("detectors", {}).get(detector_name)
        if not isinstance(detector, dict):
            raise ValueError(f"mu-resolution report lacks {detector_name} folds")
        response_path = detector.get("response_path")
        if not isinstance(response_path, str) or not response_path:
            raise ValueError(
                f"mu-resolution report lacks {detector_name} response file path"
            )
        resolved_response = (Path(path).parent / response_path).resolve()
        if not resolved_response.is_file():
            raise ValueError(f"{detector_name} response file is missing")
        from .publication import file_sha256

        if file_sha256(resolved_response) != detector.get("response_sha256"):
            raise ValueError(
                f"{detector_name} response SHA-256 differs from its file"
            )
        response_evidence = _detector_response_evidence(
            detector, detector_name, discretization_gate=tolerance
        )
        nodes = detector.get("nodes")
        if not isinstance(nodes, list):
            raise ValueError(f"{detector_name} mu-resolution fold lacks per-node metrics")
        by_id: dict[str, dict[str, Any]] = {}
        for item in nodes:
            if not isinstance(item, dict):
                raise ValueError(f"{detector_name} mu-resolution node record is malformed")
            node_id = str(item.get("node_id", "")).lower().removeprefix("0x")
            if node_id in by_id:
                raise ValueError(f"{detector_name} mu-resolution report repeats {node_id}")
            by_id[node_id] = item
        if set(by_id) != set(pilot_node_ids):
            raise ValueError(f"{detector_name} mu-resolution report misses pilot nodes")
        _close_detector_discretization_evidence(
            by_id, supplied, response_evidence, detector_name
        )
        metric_names = (
            "detector_folded_total_relative_change",
            "detector_folded_spectral_tv",
            "exact_line_angular_fold_relative_change",
        )
        maxima = {
            f"max_{name}": max(
                (float(by_id[node_id].get(name, _FINITE_FAILURE)) for node_id in pilot_node_ids),
                default=_FINITE_FAILURE,
            )
            for name in metric_names
        }
        detector_passed = all(
            math.isfinite(value) and 0.0 <= value <= tolerance for value in maxima.values()
        )
        details[detector_name] = {
            **response_evidence,
            "node_count": len(by_id),
            **maxima,
            "passed": detector_passed,
        }
        passed = passed and detector_passed
    return {"supplied": True, "passed": bool(passed), "detectors": details}


def mu_resolution_report(
    archive_path: str | Path,
    *,
    tolerance: float = 0.002,
    detector_validation_path: str | Path | None = None,
) -> dict[str, Any]:
    """Execute the 64→32 pilot sentinel decision on total continuum maps."""

    h5py, numpy = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        mu_edges = root["axes/mu_edges"][:]
        if len(mu_edges) - 1 != 64:
            raise ValueError("mu pilot decision requires a 64-bin archive")
        if not math.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("mu-resolution tolerance must be finite and positive")
        roles = _roles(root)
        complete = _transport_complete_mask(root, numpy)
        indices = [i for i, role in enumerate(roles) if role == "pilot" and complete[i]]
        if not indices:
            indices = [i for i, role in enumerate(roles) if role == "train" and complete[i]]
        if not indices:
            raise ValueError("mu pilot decision has no transport-complete pilot/train nodes")
        n = root["continuum/total/N"][indices]
        s1 = root["continuum/total/S1"][indices]
        means = numpy.divide(s1, n, out=numpy.zeros_like(s1), where=n > 0)
        n_energy = means.shape[1] // 64
        means = means.reshape(len(indices), n_energy, 64)
        coarse = means.reshape(len(indices), n_energy, 32, 2).sum(axis=3)
        expanded = numpy.repeat(coarse * 0.5, 2, axis=2)
        totals = numpy.maximum(numpy.sum(means, axis=(1, 2)), 1.0e-300)
        tv = 0.5 * numpy.sum(numpy.abs(means - expanded), axis=(1, 2)) / totals
        worst_tv = float(numpy.max(tv))
        pilot_node_ids = [f"{int(root['nodes/id'][index]):016x}" for index in indices]
        snapshot_hash = archive_snapshot_sha256(root)
        detector = _mu_detector_validation(
            detector_validation_path,
            config_sha256=str(root["manifest"].attrs["config_sha256"]),
            archive_snapshot_hash=snapshot_hash,
            pilot_node_ids=pilot_node_ids,
            tolerance=tolerance,
        )
        passed = worst_tv <= tolerance and detector["passed"]
        return {
            "schema": "CEEFIT-MU-DECISION-v1",
            "config_sha256": str(root["manifest"].attrs["config_sha256"]),
            "source_bins": 64,
            "candidate_bins": 32,
            "tolerance": tolerance,
            "max_joint_tv": worst_tv,
            "archive_snapshot_sha256": snapshot_hash,
            "source_archive": {
                "path": str(Path(archive_path).resolve()),
                "archive_snapshot_sha256": snapshot_hash,
            },
            "source_detector_validation": (
                None
                if detector_validation_path is None
                else {
                    "path": str(Path(detector_validation_path).resolve()),
                    "report_sha256": file_sha256(detector_validation_path),
                }
            ),
            "pilot_node_ids": pilot_node_ids,
            "detector_validation": detector,
            "decision": "use_32" if passed else "retain_64",
            "passed": passed,
            "note": "32 bins are selected only when joint-shape and bound HPGe/NaI folds all pass.",
        }
