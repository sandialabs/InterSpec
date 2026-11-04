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

"""Command-line entry points for the fit-first campaign."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from .allocation import (
    AcquisitionCandidate,
    ContributionMoment,
    NodeStatistics,
    RawMoment,
    SentinelMoment,
    allocate_actions,
    assess_node,
)
from .block import read_fit_block
from .contracts import RETAINED_RAW_BLOCKS
from .campaign import (
    DEFAULT_WORKER_EXECUTABLE,
    PhysicsConfig,
    advance_campaign_round,
    generate_campaign_manifest,
    load_manifest,
    replace_final_set_from_reserve,
    save_manifest,
)
from .design import DesignBounds, PhysicalNode


def _csv_ints(text: str) -> tuple[int, ...]:
    try:
        return tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error


def _load_anchors(path: str | None) -> tuple[PhysicalNode, ...] | None:
    if path is None:
        return None
    records = json.loads(Path(path).read_text(encoding="utf-8"))
    return tuple(
        PhysicalNode(
            float(item["energy_keV"]),
            float(item["distance_m"]),
            float(item["source_height_m"]),
            float(item["detector_height_m"]),
        )
        for item in records
    )


def _design_command(args: argparse.Namespace) -> int:
    bounds = DesignBounds(
        (args.energy_min, args.energy_max),
        (args.distance_min, args.distance_max),
        (args.source_height_min, args.source_height_max),
        (args.detector_height_min, args.detector_height_max),
        args.rho_min,
    )
    physics = PhysicsConfig(mu_bins=args.mu_bins)
    mu_decision = (
        None
        if args.mu_decision is None
        else json.loads(Path(args.mu_decision).read_text(encoding="utf-8"))
    )
    manifest = generate_campaign_manifest(
        mode=args.mode,
        material=args.material,
        bounds=bounds,
        physics=physics,
        candidate_count=args.candidates,
        train_count=args.train,
        validation_count=args.validation,
        final_count=args.final,
        seed=args.seed,
        anchors=_load_anchors(args.anchors),
        executable=args.executable,
        block_output_root=args.block_output_root,
        mu_resolution_decision=mu_decision,
    )
    save_manifest(manifest, args.output)
    print(
        f"wrote {args.output}: {args.train} train, {args.validation} validation, "
        f"{args.final} final nodes; {len(manifest['work'])} initial blocks"
    )
    return 0


def _commands_command(args: argparse.Namespace) -> int:
    manifest = load_manifest(args.manifest)
    roles = set(args.roles.split(",")) if args.roles else None
    source = manifest["work"] if args.ungrouped else manifest.get("worker_groups", manifest["work"])
    commands = [
        item["shell_command"]
        for item in source
        if roles is None or item["role"] in roles
    ]
    output = "\n".join(commands) + ("\n" if commands else "")
    if args.output:
        Path(args.output).write_text(output, encoding="utf-8")
    else:
        sys.stdout.write(output)
    return 0


def _campaign_round_command(args: argparse.Namespace) -> int:
    manifest = load_manifest(args.manifest)
    validation_report = None
    archive_snapshot_hash = None
    model_content_hash = None
    if args.validation_report:
        if not args.archive or not args.candidate:
            raise ValueError(
                "--validation-report requires --archive and --candidate for "
                "independent identity checks"
            )
        validation_report = json.loads(
            Path(args.validation_report).read_text(encoding="utf-8")
        )
        if not isinstance(validation_report, dict):
            raise ValueError("model-selection validation report must be a JSON object")
        from .archive import require_hdf5
        from .fit import RuntimePayload
        from .publication import archive_snapshot_sha256

        h5py, _ = require_hdf5()
        with h5py.File(args.archive, "r") as root:
            archive_snapshot_hash = archive_snapshot_sha256(root)
        candidate = RuntimePayload.load_binary(args.candidate, allow_candidate=True)
        model_content_hash = candidate.compute_model_content_sha256()
        if candidate.model_content_sha256 != model_content_hash:
            raise ValueError("candidate model-content hash is inconsistent")
        if candidate.metadata.get("config_sha256") != manifest["config"]["sha256"]:
            raise ValueError("candidate is bound to another campaign configuration")
    updated = advance_campaign_round(
        manifest,
        validation_report=validation_report,
        diagnostic_validation_failed=args.diagnostic_validation_failed,
        expected_archive_snapshot_sha256=archive_snapshot_hash,
        expected_model_content_sha256=model_content_hash,
    )
    save_manifest(updated, args.output)
    state = updated["campaign_state"]
    print(
        f"round={state['round']} action={state['last_action']} "
        f"train={len(state['training_node_ids'])} stopped={state['stopped']}"
    )
    return 0


def _refresh_final_command(args: argparse.Namespace) -> int:
    report = json.loads(Path(args.failure_report).read_text(encoding="utf-8"))
    if not isinstance(report, dict):
        raise ValueError("final failure report must be a JSON object")
    updated = replace_final_set_from_reserve(
        load_manifest(args.manifest), failure_report=report
    )
    save_manifest(updated, args.output)
    state = updated["campaign_state"]
    print(
        f"final_generation={state['final_generation']} "
        f"retired={len(state['retired_final_node_ids'])}"
    )
    return 0


def _registered_e0_command(args: argparse.Namespace) -> int:
    from .registered_e0 import generate_registered_e0_report

    report = generate_registered_e0_report(
        load_manifest(args.campaign),
        args.executable,
        args.output,
        batch_histories=args.batch_histories,
        max_histories=args.max_histories,
        threads=args.threads,
    )
    print(
        f"registered-E0 {'PASS' if report['passed'] else 'FAIL'}: "
        f"max joint TV={report['max_joint_tv']:.6g}; wrote {args.output}"
    )
    return 0 if report["passed"] else 1


def _moment(record: dict[str, Any]) -> RawMoment:
    return RawMoment(float(record["S1"]), float(record["S2"]), int(record["N"]))


def _allocation_command(args: argparse.Namespace) -> int:
    record = json.loads(Path(args.input).read_text(encoding="utf-8"))
    nodes = []
    assessments = []
    for item in record.get("nodes", []):
        node = NodeStatistics(
            str(item["node_id"]),
            str(item["role"]),
            int(item["histories"]),
            _moment(item["total"]),
            tuple(
                ContributionMoment(str(value["name"]), _moment(value))
                for value in item.get("contributions", [])
            ),
            tuple(
                SentinelMoment(
                    str(value["name"]),
                    _moment(value),
                    float(value["normalizer_mean"]),
                )
                for value in item.get("sentinels", [])
            ),
            float(item.get("wall_seconds", 0.0)),
            int(item.get("completed_blocks", 0)),
        )
        nodes.append(node)
        assessment = assess_node(node)
        assessments.append(
            {
                "node_id": node.node_id,
                "passed": assessment.passed,
                "history_budget_capped": assessment.history_budget_capped,
                "recommended_additional_blocks": assessment.recommended_additional_blocks,
                "criteria": [
                    {
                        "name": value.name,
                        "achieved_se": value.achieved_se,
                        "target_se": value.target_se,
                        "active": value.active,
                        "ratio": value.ratio,
                        "passed": value.passed,
                    }
                    for value in assessment.criteria
                ],
            }
        )
    candidates = tuple(
        AcquisitionCandidate(
            str(item["node_id"]),
            float(item["predicted_loss_reduction"]),
            float(item["predicted_cpu_seconds"]),
            blocks=int(item.get("blocks", 4)),
        )
        for item in record.get("new_candidates", [])
    )
    local_losses = {
        str(key): float(value)
        for key, value in record.get("continuation_local_validation_losses", {}).items()
    }
    actions = allocate_actions(
        nodes,
        candidates,
        max_actions=args.max_actions,
        continuation_local_losses=local_losses,
    )
    output = {
        "schema": "CEEFIT-ALLOCATION-v1",
        "assessments": assessments,
        "actions": [
            {
                "node_id": action.node_id,
                "kind": action.kind,
                "blocks": action.blocks,
                "utility": action.utility,
                "reason": action.reason,
                "predicted_loss_reduction": action.predicted_loss_reduction,
                "predicted_cpu_seconds": action.predicted_cpu_seconds,
            }
            for action in actions
        ],
    }
    encoded = json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        Path(args.output).write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    return 0


def _continue_archive_command(args: argparse.Namespace) -> int:
    from .campaign import continuation_shell_commands, schedule_archive_continuation

    updated = schedule_archive_continuation(
        load_manifest(args.manifest),
        args.archive,
        max_nodes=args.max_nodes,
        detector_validation_path=args.detector_validation,
        adaptive_utility_path=args.adaptive_utilities,
    )
    save_manifest(updated, args.output)
    roles = tuple(value for value in args.roles.split(",") if value) if args.roles else None
    commands = continuation_shell_commands(
        updated,
        roles=roles,
        grouped=not args.ungrouped,
    )
    encoded = "\n".join(commands) + ("\n" if commands else "")
    if args.commands_output:
        Path(args.commands_output).write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    schedule = updated["campaign_state"]["continuation_schedule"]
    print(
        f"continuation generation={schedule['generation']} "
        f"assessed={schedule['assessment_count']} "
        f"nodes={len(schedule['pending_block_ids_by_node'])} "
        f"blocks={len(schedule['pending_work'])}; wrote {args.output}",
        file=sys.stderr,
    )
    return 0


def _project_continuation_command(args: argparse.Namespace) -> int:
    """What bringing every node to its noise target would add, read-only.

    The same assessment `continue-archive` makes, summed over every node
    instead of one round's worth, so the checkpoint can size the adaptive
    rounds' archive before any is run. Writes nothing.
    """
    from .allocation import BLOCK_HISTORIES, history_cap_blocks
    from .campaign import schedule_archive_continuation
    from .storage import projected_archive_growth

    updated = schedule_archive_continuation(
        load_manifest(args.manifest),
        args.archive,
        adaptive_utility_path=args.adaptive_utilities,
    )
    rows = updated["campaign_state"]["continuation_schedule"]["assessments"]
    work: list[tuple[int, int]] = []
    for row in rows:
        completed = {int(value) for value in row["completed_block_ids"]}
        candidate = 1
        for _ in range(int(row["recommended_additional_blocks"])):
            while candidate in completed:
                candidate += 1
            work.append((int(str(row["node_id"]), 16), candidate))
            candidate += 1
    growth = projected_archive_growth(args.archive, work)
    roles: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        roles.setdefault(str(row["archive_role"] or row["role"]), []).append(row)
    report: dict[str, Any] = {"archive_growth": growth, "roles": {}}
    total_blocks = 0
    for role, members in sorted(roles.items()):
        cap = history_cap_blocks(role)
        extra = sum(int(row["recommended_additional_blocks"]) for row in members)
        at_cap = sum(
            1 for row in members
            if int(row["recommended_additional_blocks"])
            >= cap - len(row["completed_block_ids"])
        )
        ratios = sorted(
            float(row["worst_ratio"]) for row in members if row["worst_ratio"] is not None
        )
        total_blocks += extra
        report["roles"][role] = {
            "nodes": len(members),
            "noise_passed": sum(1 for row in members if row["passed"]),
            "history_cap_blocks": cap,
            "nodes_needing_the_cap": at_cap,
            "additional_blocks": extra,
            "worst_ratio_median": ratios[len(ratios) // 2] if ratios else None,
            "worst_ratio_max": ratios[-1] if ratios else None,
        }
        print(
            f"{role:10s} {len(members):4d} nodes: {report['roles'][role]['noise_passed']} at "
            f"their noise target, {at_cap} need the full {cap}-block cap; "
            f"+{extra} blocks"
        )
    report["additional_blocks"] = total_blocks
    report["additional_archive_bytes"] = growth["total_bytes"]
    print(
        f"total +{total_blocks} blocks ({total_blocks * BLOCK_HISTORIES / 1e9:.1f} G histories): "
        f"archive about +{growth['expected_bytes'] / 1e9:.1f} GB, at most "
        f"+{growth['total_bytes'] / 1e9:.1f} GB (blocks beyond the first "
        f"{RETAINED_RAW_BLOCKS} of a node fold into its running sums; raw records "
        f"+{growth['raw_record_bytes'] / 1e9:.1f} GB)"
    )
    if not growth["reuses_freed_space"]:
        print(
            "this archive predates persistent free-space tracking: space freed by "
            "replaced fold states is never reused (start a new archive to avoid it)"
        )
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0


def _inspect_block_command(args: argparse.Namespace) -> int:
    block = read_fit_block(args.block)
    summary = {
        "node_id": f"{block.node_id:016x}",
        "block_id": block.block_id,
        "config_sha256": block.manifest["config_sha256"],
        "role": block.node["role"],
        "histories": block.node["histories"],
        "components": [item["name"] for item in block.components],
        "energy_axes": [item["name"] for item in block.axes["energy"]],
        "lines": [item["name"] for item in block.axes["lines"]],
        "payload_sha256": block.payload_sha256,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def _ingest_command(args: argparse.Namespace) -> int:
    from .archive import ingest_blocks
    from .storage import create_archive_seal

    if args.full_file_hash and not args.seal:
        raise ValueError("--full-file-hash requires --seal")
    result = ingest_blocks(
        args.blocks,
        args.output,
        delete_source_blocks=args.delete_source_blocks,
        fold_blocks=not args.keep_raw_blocks,
    )
    if args.seal:
        result["archive_seal"] = create_archive_seal(
            args.output, full_file_hash=args.full_file_hash
        )["output"]
    print(json.dumps(result, sort_keys=True))
    return 0


def _seal_archive_command(args: argparse.Namespace) -> int:
    from .storage import create_archive_seal

    result = create_archive_seal(
        args.archive,
        args.output,
        full_file_hash=args.full_file_hash,
        incremental=args.incremental,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


def _verify_archive_seal_command(args: argparse.Namespace) -> int:
    from .storage import verify_archive_seal

    result = verify_archive_seal(
        args.archive,
        args.seal,
        full_file_hash=args.full_file_hash,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


def _fit_cache_command(args: argparse.Namespace) -> int:
    from .storage import create_fit_cache

    result = create_fit_cache(
        args.archive,
        args.output,
        full_file_hash=args.full_file_hash,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


def _verify_fit_cache_command(args: argparse.Namespace) -> int:
    from .storage import verify_fit_cache

    result = verify_fit_cache(args.cache, full_file_hash=args.full_file_hash)
    print(json.dumps(result, sort_keys=True))
    return 0


def _run_streaming_command(args: argparse.Namespace) -> int:
    from .storage import run_streaming_campaign

    manifest = load_manifest(args.manifest)
    roles = None if not args.roles else tuple(args.roles.split(","))
    result = run_streaming_campaign(
        manifest,
        args.archive,
        args.spool,
        roles=roles,
        spool_max_bytes=int(args.spool_max_gib * 1024**3),
        delete_packed_blocks=args.delete_packed_blocks,
        seal_output=not args.no_seal,
        incremental_seal=args.incremental_seal,
        max_groups=args.max_groups,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


def _fit_command(args: argparse.Namespace) -> int:
    from .fit import RuntimePayload, select_rank_and_level, require_numpy

    numpy = require_numpy()
    data = numpy.load(args.input)
    campaign = load_manifest(args.campaign)
    q_record = campaign["design"]["q_domain"]
    brems_threshold = float(
        campaign["config"]["physics"]["brems_electron_floor_keV"]
    )
    variances = data["train_variances"] if "train_variances" in data else None
    selected, candidates = select_rank_and_level(
        args.component,
        args.axis_kind,
        data["train_q"],
        data["train_masses"],
        data["validation_q"],
        data["validation_masses"],
        registered_energy_edges=data["registered_energy_edges"],
        mu_edges=data["mu_edges"],
        train_variances=variances,
        ranks=args.ranks,
        levels=args.levels,
        nmf_iterations=args.nmf_iterations,
        seed=args.seed,
        source_energy_threshold_keV=(
            brems_threshold if args.axis_kind == "brems" else None
        ),
        q_domain=(q_record["min"], q_record["max"]),
    )
    bounds = campaign["design"]["bounds"]
    payload = RuntimePayload(
        {args.component: selected.model},
        metadata={
            "config_sha256": campaign["config"]["sha256"],
            "material_name": campaign["config"]["material"].title(),
            "geometry_tier": campaign["config"]["mode"],
            "normalization": "s_hat_to_direct_fluence",
            "physical_domain": {
                "min": [bounds["energy_keV"][0], bounds["distance_m"][0], bounds["source_height_m"][0], bounds["detector_height_m"][0]],
                "max": [bounds["energy_keV"][1], bounds["distance_m"][1], bounds["source_height_m"][1], bounds["detector_height_m"][1]],
            },
            "q_domain": q_record,
            "brems_activation_threshold_keV": brems_threshold,
            "model_inventory": {
                "continuum": [
                    {"name": args.component, "axis_kind": args.axis_kind}
                ],
                "lines": [],
            },
            "validation_metadata": {
                "selection_passed": selected.passed,
                "validation_max_yield_relative_error": selected.max_yield_relative_error,
                "validation_max_joint_tv": selected.max_joint_tv,
                "candidates": [
                {
                    "rank": item.model.rank,
                    "level": item.model.level,
                    "parameter_count": item.model.parameter_count,
                    "passed": item.passed,
                    "max_yield_relative_error": item.max_yield_relative_error,
                    "max_joint_tv": item.max_joint_tv,
                }
                for item in candidates
                ],
            },
        },
    )
    payload.save_binary(args.output)
    print(
        f"selected rank={selected.model.rank} level={selected.model.level}; "
        f"validation pass={selected.passed}; wrote candidate {args.output}"
    )
    return 0


def _fit_archive_command(args: argparse.Namespace) -> int:
    from .workflow import fit_archive

    report = fit_archive(
        args.archive,
        args.campaign,
        args.output,
        ranks=args.ranks,
        levels=args.levels,
        nmf_iterations=args.nmf_iterations,
        seed=args.seed,
        detector_validation_path=args.detector_validation,
        publication_evidence_path=args.publication_evidence,
        acquisition_utilities_path=args.acquisition_utilities,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.report:
        Path(args.report).write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    return 0 if report.get("acceptance_passed") else 2


def _promote_candidate_command(args: argparse.Namespace) -> int:
    from .workflow import promote_candidate

    report = promote_candidate(
        args.candidate,
        args.archive,
        args.campaign,
        args.detector_validation,
        args.publication_evidence,
        args.output,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.report:
        Path(args.report).write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    return 0 if report.get("acceptance_passed") else 2


def _mu_decision_command(args: argparse.Namespace) -> int:
    from .workflow import mu_resolution_report

    report = mu_resolution_report(
        args.archive,
        tolerance=args.tolerance,
        detector_validation_path=args.detector_validation,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        Path(args.output).write_text(encoded, encoding="utf-8")
    else:
        sys.stdout.write(encoded)
    return 0


def _parity_report_command(args: argparse.Namespace) -> int:
    from .parity import generate_parity_report

    report = generate_parity_report(
        args.candidate,
        args.cpp_evaluator,
        args.output,
        tolerance=args.tolerance,
        timeout_seconds=args.timeout,
    )
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "query_count": report["query_count"],
                "compared_value_count": report["compared_value_count"],
                "max_absolute_difference": report["max_absolute_difference"],
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def _detector_fold_report_command(args: argparse.Namespace) -> int:
    from .detector_fold import generate_detector_fold_report

    report = generate_detector_fold_report(
        args.candidate,
        args.archive,
        args.campaign,
        args.hpge_response,
        args.nai_response,
        args.output,
    )
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "final_node_count": len(report["final_node_ids"]),
                "response_sha256": {
                    name: value["response_sha256"]
                    for name, value in report["detectors"].items()
                },
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def _mu_fold_report_command(args: argparse.Namespace) -> int:
    from .detector_fold import generate_mu_resolution_validation_report

    report = generate_mu_resolution_validation_report(
        args.archive,
        args.hpge_response,
        args.nai_response,
        args.output,
        tolerance=args.tolerance,
    )
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "pilot_node_count": len(report["pilot_node_ids"]),
                "response_sha256": {
                    name: value["response_sha256"]
                    for name, value in report["detectors"].items()
                },
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def _calibrate_uncertainty_command(args: argparse.Namespace) -> int:
    from .calibration import calibrate_uncertainty_archive

    report = calibrate_uncertainty_archive(args.archive, args.output)
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "repeated_null_runs": report["repeated_null_runs"],
                "s2_passed": report["s2_passed"],
                "coarse_covariance_passed": report["coarse_covariance_passed"],
                "replica_covariance_passed": report["replica_covariance_passed"],
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def _geant4_anchor_report_command(args: argparse.Namespace) -> int:
    from .geant4_anchor import generate_geant4_anchor_report

    report = generate_geant4_anchor_report(
        args.archive,
        args.manifest,
        args.output,
    )
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "anchor_count": report["anchor_count"],
                "reference_file_count": report["reference_file_count"],
                "reference_event_count": report["reference_event_count"],
                "max_abs_standardized_difference": report[
                    "max_abs_standardized_difference"
                ],
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def _assemble_evidence_command(args: argparse.Namespace) -> int:
    from .publication import assemble_candidate_publication_evidence

    report = assemble_candidate_publication_evidence(
        args.candidate,
        args.archive,
        args.campaign,
        args.output,
        report_paths={
            "registered_e0_gate": args.registered_e0_report,
            "uncertainty_calibration": args.uncertainty_calibration_report,
            "geant4_anchor_agreement": args.geant4_anchor_report,
            "python_cpp_parity": args.python_cpp_parity_report,
        },
    )
    print(
        json.dumps(
            {
                "output": args.output,
                "passed": report["passed"],
                "source_reports": {
                    name: value["report_sha256"]
                    for name, value in report["checks"].items()
                },
            },
            sort_keys=True,
        )
    )
    return 0 if report["passed"] else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    design = subparsers.add_parser("design", help="write a deterministic campaign manifest")
    design.add_argument("--output", required=True)
    design.add_argument("--mode", choices=("air", "field"), default="field")
    design.add_argument("--material", choices=("soil", "concrete"), default="soil")
    design.add_argument("--candidates", type=int, default=10_000)
    design.add_argument("--train", type=int, default=256)
    design.add_argument("--validation", type=int, default=64)
    design.add_argument("--final", type=int, default=64)
    design.add_argument("--seed", type=int, default=0xCEE10)
    design.add_argument("--anchors", help="JSON list of additional/replacement forced train nodes")
    design.add_argument("--mu-bins", type=int, choices=(32, 64), default=64)
    design.add_argument(
        "--mu-decision",
        help="passing CEEFIT-MU-DECISION-v1 report required for 32 production bins",
    )
    design.add_argument("--energy-min", type=float, default=59.5)
    design.add_argument("--energy-max", type=float, default=2614.0)
    design.add_argument("--distance-min", type=float, default=0.5)
    design.add_argument("--distance-max", type=float, default=49.87)
    design.add_argument("--source-height-min", type=float, default=0.3)
    design.add_argument("--source-height-max", type=float, default=4.16)
    design.add_argument("--detector-height-min", type=float, default=0.3)
    design.add_argument("--detector-height-max", type=float, default=4.16)
    design.add_argument(
        "--rho-min",
        type=float,
        default=0.01,
        help="minimum supported transverse separation in metres for log(rho)",
    )
    design.add_argument(
        "--executable",
        default=DEFAULT_WORKER_EXECUTABLE,
        help=f"env_kernel_bake worker binary (default: {DEFAULT_WORKER_EXECUTABLE})",
    )
    design.add_argument("--block-output-root", default="ceefit_blocks")
    design.set_defaults(handler=_design_command)

    commands = subparsers.add_parser("commands", help="emit shell-safe worker commands")
    commands.add_argument("manifest")
    commands.add_argument("--roles", help="comma-separated role filter")
    commands.add_argument("--output")
    commands.add_argument(
        "--ungrouped",
        action="store_true",
        help="emit one process per block instead of persistent four-block node jobs",
    )
    commands.set_defaults(handler=_commands_command)

    campaign_round = subparsers.add_parser(
        "campaign-round",
        help="persist a bound validation round before unified adaptive allocation",
    )
    campaign_round.add_argument("--manifest", required=True)
    campaign_round.add_argument("--output", required=True)
    round_outcome = campaign_round.add_mutually_exclusive_group(required=True)
    round_outcome.add_argument(
        "--validation-report",
        help=(
            "CEEFIT-MODEL-SELECTION-VALIDATION-v1 report; independently "
            "checked against --archive and --candidate"
        ),
    )
    round_outcome.add_argument(
        "--diagnostic-validation-failed",
        action="store_true",
        help="record an explicit failure when no complete fit report exists",
    )
    campaign_round.add_argument(
        "--archive", help="HDF5 archive used to independently verify a validation report"
    )
    campaign_round.add_argument(
        "--candidate", help="candidate payload used to independently verify model content"
    )
    campaign_round.set_defaults(handler=_campaign_round_command)

    refresh_final = subparsers.add_parser(
        "refresh-final", help="retire failed final nodes and draw a fresh Sobol-reserve set"
    )
    refresh_final.add_argument("--manifest", required=True)
    refresh_final.add_argument(
        "--failure-report",
        required=True,
        help="bound fit report whose final policy authorizes this retirement",
    )
    refresh_final.add_argument("--output", required=True)
    refresh_final.set_defaults(handler=_refresh_final_command)

    allocate = subparsers.add_parser("allocate", help="rank node continuation/acquisition work")
    allocate.add_argument("--input", required=True)
    allocate.add_argument("--output")
    allocate.add_argument("--max-actions", type=int, default=64)
    allocate.set_defaults(handler=_allocation_command)

    continue_archive = subparsers.add_parser(
        "continue-archive",
        help="choose and emit archive-bound continuation or 64-node acquisition work",
    )
    continue_archive.add_argument("--manifest", required=True)
    continue_archive.add_argument("--archive", required=True)
    continue_archive.add_argument("--output", required=True)
    continue_archive.add_argument(
        "--detector-validation",
        help=(
            "archive-bound CEEFIT detector-fold report; its per-final-node "
            "replica SE participates in continuation allocation"
        ),
    )
    continue_archive.add_argument(
        "--adaptive-utilities",
        help=(
            "archive/model/round-bound CEEFIT-ADAPTIVE-UTILITY-v1 report; "
            "required after failed model-selection validation"
        ),
    )
    continue_archive.add_argument(
        "--commands-output",
        help="write pending commands here instead of standard output",
    )
    continue_archive.add_argument(
        "--roles", help="comma-separated archive roles to emit after assessing every active node"
    )
    continue_archive.add_argument(
        "--max-nodes",
        type=int,
        default=64,
        help="maximum node jobs to issue in this scheduling pass (default: 64)",
    )
    continue_archive.add_argument(
        "--ungrouped",
        action="store_true",
        help="emit one process per block instead of one process per contiguous node run",
    )
    continue_archive.set_defaults(handler=_continue_archive_command)

    project = subparsers.add_parser(
        "project-continuation",
        help="read-only: blocks and archive growth to bring every node to its noise target",
    )
    project.add_argument("--manifest", required=True)
    project.add_argument("--archive", required=True)
    project.add_argument("--adaptive-utilities")
    project.add_argument("--output", help="also write the projection as JSON")
    project.set_defaults(handler=_project_continuation_command)

    inspect_block = subparsers.add_parser("inspect-block", help="validate and summarize a block")
    inspect_block.add_argument("block")
    inspect_block.set_defaults(handler=_inspect_block_command)

    ingest = subparsers.add_parser("ingest", help="merge blocks into chunked HDF5")
    ingest.add_argument("--output", required=True)
    ingest.add_argument(
        "--delete-source-blocks",
        action="store_true",
        help=(
            "delete worker blocks only after their embedded records and affected "
            "node aggregates pass a durable post-write verification"
        ),
    )
    ingest.add_argument(
        "--keep-raw-blocks",
        action="store_true",
        help=(
            "store every block as a raw record instead of folding blocks beyond "
            "each node's first four into its running sums (tests, diagnostics)"
        ),
    )
    ingest.add_argument(
        "--seal",
        action="store_true",
        help="perform one full scientific closure audit and write the default seal",
    )
    ingest.add_argument(
        "--full-file-hash",
        action="store_true",
        help="also compute a redundant bytewise archive hash when sealing",
    )
    ingest.add_argument("blocks", nargs="+")
    ingest.set_defaults(handler=_ingest_command)

    streaming = subparsers.add_parser(
        "run-streaming",
        help="run deterministic workers through a bounded local spool into HDF5",
    )
    streaming.add_argument("--manifest", required=True)
    streaming.add_argument("--archive", required=True)
    streaming.add_argument("--spool", required=True)
    streaming.add_argument(
        "--roles", help="comma-separated roles to run (default: every manifest role)"
    )
    streaming.add_argument(
        "--spool-max-gib",
        type=float,
        default=4.0,
        help="maximum generated-but-not-packed spool volume (default: 4 GiB)",
    )
    streaming.add_argument(
        "--delete-packed-blocks",
        action="store_true",
        help="transactionally remove redundant worker blocks after verified packing",
    )
    streaming.add_argument(
        "--no-seal",
        action="store_true",
        help="skip the final full archive audit/seal (useful between campaign rounds)",
    )
    streaming.add_argument(
        "--incremental-seal",
        action="store_true",
        help=(
            "seal by re-verifying only the nodes whose blocks or rows changed since "
            "the last seal (adaptive rounds); publication still requires a full seal"
        ),
    )
    streaming.add_argument(
        "--max-groups",
        type=int,
        help="run only this many worker groups (diagnostic/resume testing)",
    )
    streaming.set_defaults(handler=_run_streaming_command)

    seal_archive = subparsers.add_parser(
        "seal-archive",
        help="fully audit a canonical archive once and cache its scientific identity",
    )
    seal_archive.add_argument("--archive", required=True)
    seal_archive.add_argument("--output")
    seal_archive.add_argument("--full-file-hash", action="store_true")
    seal_archive.add_argument(
        "--incremental",
        action="store_true",
        help="re-verify only nodes whose closure identity changed since the last seal",
    )
    seal_archive.set_defaults(handler=_seal_archive_command)

    verify_seal = subparsers.add_parser(
        "verify-archive-seal",
        help="quickly validate a stat-bound archive seal",
    )
    verify_seal.add_argument("--archive", required=True)
    verify_seal.add_argument("--seal")
    verify_seal.add_argument("--full-file-hash", action="store_true")
    verify_seal.set_defaults(handler=_verify_archive_seal_command)

    fit_cache = subparsers.add_parser(
        "fit-cache",
        help="copy exact fitter aggregates into a disposable block-free HDF5 cache",
    )
    fit_cache.add_argument("--archive", required=True)
    fit_cache.add_argument("--output", required=True)
    fit_cache.add_argument("--full-file-hash", action="store_true")
    fit_cache.set_defaults(handler=_fit_cache_command)

    verify_cache = subparsers.add_parser(
        "verify-fit-cache", help="validate a fit cache and its source-snapshot binding"
    )
    verify_cache.add_argument("--cache", required=True)
    verify_cache.add_argument("--full-file-hash", action="store_true")
    verify_cache.set_defaults(handler=_verify_fit_cache_command)

    fit = subparsers.add_parser("fit", help="fit one positive continuum component from NPZ")
    fit.add_argument("--input", required=True)
    fit.add_argument("--campaign", required=True)
    fit.add_argument("--output", required=True)
    fit.add_argument("--component", required=True)
    fit.add_argument("--axis-kind", choices=("primary", "annihilation", "brems"), required=True)
    fit.add_argument("--ranks", type=_csv_ints, default=(4, 8, 12, 16, 24))
    fit.add_argument("--levels", type=_csv_ints, default=(1, 2, 3, 4))
    fit.add_argument("--nmf-iterations", type=int, default=400)
    fit.add_argument("--seed", type=int, default=0xCEE10)
    fit.set_defaults(handler=_fit_command)

    fit_archive_parser = subparsers.add_parser(
        "fit-archive", help="fit all supported components/lines from an HDF5 archive"
    )
    fit_archive_parser.add_argument("--archive", required=True)
    fit_archive_parser.add_argument("--campaign", required=True)
    fit_archive_parser.add_argument("--output", required=True)
    fit_archive_parser.add_argument("--report")
    fit_archive_parser.add_argument(
        "--detector-validation",
        help="external CEEFIT-DETECTOR-FOLD-VALIDATION-v1 HPGe/NaI report",
    )
    fit_archive_parser.add_argument(
        "--publication-evidence",
        help="identity-bound CEEFIT-PUBLICATION-EVIDENCE-v1 report",
    )
    fit_archive_parser.add_argument(
        "--acquisition-utilities",
        help="write built-in validation-loss-reduction/CPU utilities for campaign-round",
    )
    fit_archive_parser.add_argument("--ranks", type=_csv_ints, default=(4, 8, 12, 16, 24))
    fit_archive_parser.add_argument("--levels", type=_csv_ints, default=(1, 2, 3, 4))
    fit_archive_parser.add_argument("--nmf-iterations", type=int, default=400)
    fit_archive_parser.add_argument("--seed", type=int, default=0xCEE10)
    fit_archive_parser.set_defaults(handler=_fit_archive_command)

    promote = subparsers.add_parser(
        "promote-candidate",
        help="validate external evidence and publish an existing candidate without refitting",
    )
    promote.add_argument("--candidate", required=True)
    promote.add_argument("--archive", required=True)
    promote.add_argument("--campaign", required=True)
    promote.add_argument("--detector-validation", required=True)
    promote.add_argument("--publication-evidence", required=True)
    promote.add_argument("--output", required=True)
    promote.add_argument("--report")
    promote.set_defaults(handler=_promote_candidate_command)

    mu_decision = subparsers.add_parser(
        "mu-decision", help="run the 64-to-32 pilot shape/fold decision"
    )
    mu_decision.add_argument("--archive", required=True)
    mu_decision.add_argument("--tolerance", type=float, default=0.002)
    mu_decision.add_argument(
        "--detector-validation",
        help="CEEFIT-MU-RESOLUTION-VALIDATION-v1 HPGe/NaI comparison",
    )
    mu_decision.add_argument("--output")
    mu_decision.set_defaults(handler=_mu_decision_command)

    parity = subparsers.add_parser(
        "parity-report",
        help="compare one serialized candidate through Python and C++ evaluators",
    )
    parity.add_argument("--candidate", required=True)
    parity.add_argument(
        "--cpp-evaluator",
        required=True,
        help="env_fit_payload_eval binary (build-rel/tools/full_spectrum/)",
    )
    parity.add_argument("--output", required=True)
    parity.add_argument("--tolerance", type=float, default=1.0e-9)
    parity.add_argument("--timeout", type=float, default=120.0)
    parity.set_defaults(handler=_parity_report_command)

    detector_fold = subparsers.add_parser(
        "detector-fold-report",
        help="fold the exact candidate/final set through bound HPGe and NaI responses",
    )
    detector_fold.add_argument("--candidate", required=True)
    detector_fold.add_argument("--archive", required=True)
    detector_fold.add_argument("--campaign", required=True)
    detector_fold.add_argument("--hpge-response", required=True)
    detector_fold.add_argument("--nai-response", required=True)
    detector_fold.add_argument("--output", required=True)
    detector_fold.set_defaults(handler=_detector_fold_report_command)

    mu_fold = subparsers.add_parser(
        "mu-fold-report",
        help="compare 64-bin pilot folds with conservative 32-bin HPGe/NaI folds",
    )
    mu_fold.add_argument("--archive", required=True)
    mu_fold.add_argument("--hpge-response", required=True)
    mu_fold.add_argument("--nai-response", required=True)
    mu_fold.add_argument("--output", required=True)
    mu_fold.add_argument("--tolerance", type=float, default=0.002)
    mu_fold.set_defaults(handler=_mu_fold_report_command)

    calibration = subparsers.add_parser(
        "calibrate-uncertainty",
        help="calibrate S2 and covariance from disjoint same-node block pairs",
    )
    calibration.add_argument("--archive", required=True)
    calibration.add_argument("--output", required=True)
    calibration.set_defaults(handler=_calibrate_uncertainty_command)

    geant4_anchor = subparsers.add_parser(
        "geant4-anchor-report",
        help="compute archive-bound agreement from raw Livermore reference arms",
    )
    geant4_anchor.add_argument("--archive", required=True)
    geant4_anchor.add_argument(
        "--manifest",
        required=True,
        help="CEEFIT-GEANT4-ANCHOR-MANIFEST-v1 input with source CSV hashes",
    )
    geant4_anchor.add_argument("--output", required=True)
    geant4_anchor.set_defaults(handler=_geant4_anchor_report_command)

    registered_e0 = subparsers.add_parser(
        "registered-e0-report",
        help="run the full-coverage registered-E0 gate and bind its invocation",
    )
    registered_e0.add_argument("--campaign", required=True)
    registered_e0.add_argument(
        "--executable",
        required=True,
        help="env_registered_e0_gate binary (build-rel/tools/full_spectrum/)",
    )
    registered_e0.add_argument("--output", required=True)
    registered_e0.add_argument("--batch-histories", type=int, default=1_000_000)
    registered_e0.add_argument("--max-histories", type=int, default=128_000_000)
    registered_e0.add_argument("--threads", type=int, default=1)
    registered_e0.set_defaults(handler=_registered_e0_command)

    evidence = subparsers.add_parser(
        "assemble-evidence",
        help="bind and hash the four detailed publication reports",
    )
    evidence.add_argument("--candidate", required=True)
    evidence.add_argument("--archive", required=True)
    evidence.add_argument("--campaign", required=True)
    evidence.add_argument("--registered-e0-report", required=True)
    evidence.add_argument("--uncertainty-calibration-report", required=True)
    evidence.add_argument("--geant4-anchor-report", required=True)
    evidence.add_argument("--python-cpp-parity-report", required=True)
    evidence.add_argument("--output", required=True)
    evidence.set_defaults(handler=_assemble_evidence_command)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (ValueError, RuntimeError, OSError) as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
