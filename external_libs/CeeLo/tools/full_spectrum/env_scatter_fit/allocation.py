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

"""Noise-aware history and adaptive-node allocation decisions."""

from __future__ import annotations

from dataclasses import dataclass, field
import heapq
import math
from typing import Iterable, Mapping, Sequence


BLOCK_HISTORIES = 262_144
LOGICAL_REPLICAS = 32
MIN_INITIAL_BLOCKS = 4
TRAIN_HISTORY_CAP = 64_000_000
VALIDATION_HISTORY_CAP = 128_000_000


@dataclass(frozen=True)
class RawMoment:
    """Mergeable scalar sufficient statistics.

    For transport observables, ``histories`` is the number of primary
    histories and ``S2`` contains squared per-history contributions.  A
    :class:`SentinelMoment` instead uses the same algebra on the 32 logical
    *replica means*: in that context ``histories`` is exactly 32.  Keeping the
    raw ``N/S1/S2`` form makes both uncertainties auditable without storing
    rounded standard errors.
    """

    s1: float
    s2: float
    histories: int

    def validate(self) -> None:
        if self.histories < 0:
            raise ValueError("history count must be non-negative")
        if not math.isfinite(self.s1) or not math.isfinite(self.s2):
            raise ValueError("moments must be finite")
        if self.s2 < 0.0:
            raise ValueError("S2 must be non-negative")
        if self.histories == 0:
            if self.s1 != 0.0 or self.s2 != 0.0:
                raise ValueError("zero-count moments must have zero S1 and S2")
            return
        minimum_s2 = self.s1 * self.s1 / self.histories
        tolerance = 2.0e-12 * max(1.0, abs(self.s2), abs(minimum_s2))
        if self.s2 + tolerance < minimum_s2:
            raise ValueError("S2 is inconsistent with N and S1")

    @property
    def mean(self) -> float:
        self.validate()
        return self.s1 / self.histories if self.histories else math.nan

    @property
    def standard_error(self) -> float:
        self.validate()
        if self.histories < 2:
            return math.inf
        numerator = max(0.0, self.s2 - self.s1 * self.s1 / self.histories)
        return math.sqrt(numerator / (self.histories * (self.histories - 1)))

    def merge(self, other: "RawMoment") -> "RawMoment":
        self.validate()
        other.validate()
        return RawMoment(
            self.s1 + other.s1,
            self.s2 + other.s2,
            self.histories + other.histories,
        )


@dataclass(frozen=True)
class ContributionMoment:
    name: str
    moment: RawMoment


@dataclass(frozen=True)
class SentinelMoment:
    """Detector-fold or joint-shape uncertainty from 32 logical replicas.

    ``moment`` must be formed from the 32 per-replica observable means, not
    from individual primary histories.  Its standard error therefore retains
    arbitrary rebinning/fold covariance exactly at the replica resolution.
    ``normalizer_mean`` is the corresponding folded or joint-shape total.
    """

    name: str
    moment: RawMoment
    normalizer_mean: float


@dataclass(frozen=True)
class AllocationCriteria:
    total_relative_se: float = 0.0025
    material_contribution_fraction: float = 0.01
    contribution_absolute_se_fraction: float = 0.001
    sentinel_absolute_se_fraction: float = 0.0025

    def validate(self) -> None:
        for value in (
            self.total_relative_se,
            self.material_contribution_fraction,
            self.contribution_absolute_se_fraction,
            self.sentinel_absolute_se_fraction,
        ):
            if not (math.isfinite(value) and value > 0.0):
                raise ValueError("allocation tolerances must be finite and positive")


@dataclass(frozen=True)
class CriterionResult:
    name: str
    achieved_se: float
    target_se: float
    active: bool

    @property
    def ratio(self) -> float:
        if not self.active:
            return 0.0
        if self.target_se <= 0.0:
            return math.inf
        return self.achieved_se / self.target_se

    @property
    def passed(self) -> bool:
        return not self.active or self.ratio <= 1.0


@dataclass(frozen=True)
class NodeStatistics:
    node_id: str
    role: str
    histories: int
    total: RawMoment
    contributions: tuple[ContributionMoment, ...] = ()
    sentinels: tuple[SentinelMoment, ...] = ()
    wall_seconds: float = 0.0
    completed_blocks: int = 0
    # A transport failure/max-vertex termination is not a history-budget cap:
    # it makes the node unusable and must never be cured by scheduling blocks.
    transport_failed: bool = False
    # Archive-derived replica and stopping-vector diagnostics sometimes have
    # exact uncertainty estimates without a scalar per-history RawMoment (for
    # example a normalized joint-shape TV).  Keep those checks explicit rather
    # than manufacturing pseudo histories merely to pass them through the
    # scalar contribution/sentinel interfaces.
    additional_criteria: tuple[CriterionResult, ...] = ()

    def validate(self) -> None:
        if self.role not in {"pilot", "train", "validation", "final", "diagnostic"}:
            raise ValueError(f"unknown node role {self.role!r}")
        if self.histories != self.total.histories:
            raise ValueError("node and total history counts disagree")
        if self.completed_blocks < 0 or self.wall_seconds < 0.0:
            raise ValueError("blocks and wall time must be non-negative")
        if not isinstance(self.transport_failed, bool):
            raise ValueError("transport_failed must be boolean")
        self.total.validate()
        for item in self.contributions:
            item.moment.validate()
            if item.moment.histories != self.histories:
                raise ValueError(f"contribution {item.name!r} history count disagrees")
        for item in self.sentinels:
            item.moment.validate()
            if item.moment.histories != LOGICAL_REPLICAS:
                raise ValueError(
                    f"sentinel {item.name!r} must contain exactly "
                    f"{LOGICAL_REPLICAS} logical replica means"
                )
            if not math.isfinite(item.normalizer_mean) or item.normalizer_mean < 0.0:
                raise ValueError(f"sentinel {item.name!r} has invalid normalizer")
        names: set[str] = set()
        for item in self.additional_criteria:
            if not item.name or item.name in names:
                raise ValueError("additional allocation criterion names must be unique")
            names.add(item.name)
            if item.achieved_se < 0.0 or item.target_se < 0.0:
                raise ValueError(
                    f"additional criterion {item.name!r} has a negative uncertainty"
                )
            if math.isnan(item.achieved_se) or math.isnan(item.target_se):
                raise ValueError(
                    f"additional criterion {item.name!r} has a NaN uncertainty"
                )


@dataclass(frozen=True)
class NodeAssessment:
    node_id: str
    criteria: tuple[CriterionResult, ...]
    recommended_additional_blocks: int
    minimum_blocks_satisfied: bool
    history_budget_exhausted: bool
    history_budget_capped: bool
    transport_failed: bool

    @property
    def passed(self) -> bool:
        return (
            not self.transport_failed
            and self.minimum_blocks_satisfied
            and all(item.passed for item in self.criteria)
        )

    @property
    def usable(self) -> bool:
        """A precision-capped node stays usable with its stored uncertainty."""

        return not self.transport_failed

    @property
    def stop_reason(self) -> str:
        if self.transport_failed:
            return "transport_failure"
        if self.history_budget_capped:
            return "history_budget_cap"
        if self.passed:
            return "precision_met"
        return "continue"

    @property
    def worst_ratio(self) -> float:
        return max((item.ratio for item in self.criteria), default=0.0)


def history_cap(role: str) -> int:
    if role in {"pilot", "train", "diagnostic"}:
        return TRAIN_HISTORY_CAP
    if role in {"validation", "final"}:
        return VALIDATION_HISTORY_CAP
    raise ValueError(f"unknown node role {role!r}")


def history_cap_blocks(role: str, block_histories: int = BLOCK_HISTORIES) -> int:
    """Maximum complete blocks whose histories do not exceed the role cap."""

    if block_histories <= 0:
        raise ValueError("block_histories must be positive")
    return history_cap(role) // block_histories


def effective_history_cap(role: str, block_histories: int = BLOCK_HISTORIES) -> int:
    """Reachable history cap after enforcing indivisible worker blocks."""

    return history_cap_blocks(role, block_histories) * block_histories


def assess_node(
    stats: NodeStatistics,
    criteria: AllocationCriteria = AllocationCriteria(),
    *,
    block_histories: int = BLOCK_HISTORIES,
    allow_existing_over_cap: bool = False,
) -> NodeAssessment:
    """Evaluate absolute uncertainty gates and estimate the next block count."""

    stats.validate()
    criteria.validate()
    if block_histories <= 0:
        raise ValueError("block_histories must be positive")
    if stats.histories % block_histories != 0:
        raise ValueError("node histories must be an integral number of blocks")
    current_blocks = stats.histories // block_histories
    if stats.completed_blocks not in {0, current_blocks}:
        raise ValueError("completed block count does not match node histories")

    total_mean = abs(stats.total.mean) if stats.histories else 0.0
    total_se = stats.total.standard_error
    checks: list[CriterionResult] = [
        CriterionResult(
            "total_yield",
            total_se,
            criteria.total_relative_se * total_mean,
            True,
        )
    ]
    for contribution in stats.contributions:
        contribution_mean = abs(contribution.moment.mean)
        active = total_mean > 0.0 and (
            contribution_mean >= criteria.material_contribution_fraction * total_mean
        )
        checks.append(
            CriterionResult(
                f"contribution:{contribution.name}",
                contribution.moment.standard_error,
                criteria.contribution_absolute_se_fraction * total_mean,
                active,
            )
        )
    for sentinel in stats.sentinels:
        scale = abs(sentinel.normalizer_mean)
        checks.append(
            CriterionResult(
                f"sentinel:{sentinel.name}",
                sentinel.moment.standard_error,
                criteria.sentinel_absolute_se_fraction * scale,
                scale > 0.0,
            )
        )
    checks.extend(stats.additional_criteria)

    minimum_blocks = max(0, MIN_INITIAL_BLOCKS - current_blocks)
    worst = max((check.ratio for check in checks if check.active), default=0.0)
    if not math.isfinite(worst):
        desired_total_blocks = max(MIN_INITIAL_BLOCKS, current_blocks + 1)
    else:
        # Standard errors scale as N^-1/2.  This is the whole projected need;
        # a scheduling pass issues at most a doubling of it
        # (continuation_block_count), and the estimate is recomputed after
        # every pass rather than committing a large run from a noisy pilot.
        desired_total_blocks = math.ceil(max(MIN_INITIAL_BLOCKS, current_blocks * worst * worst))
    cap_blocks = history_cap_blocks(stats.role, block_histories)
    if current_blocks > cap_blocks and not allow_existing_over_cap:
        raise ValueError("node histories exceed the integral role history cap")
    remaining_blocks = max(0, cap_blocks - current_blocks)
    requested_blocks = max(minimum_blocks, desired_total_blocks - current_blocks)
    additional = 0 if stats.transport_failed else min(requested_blocks, remaining_blocks)
    uncertainty_failed = any(not check.passed for check in checks)
    budget_exhausted = current_blocks >= cap_blocks
    budget_capped = (
        not stats.transport_failed and budget_exhausted and uncertainty_failed
    )
    return NodeAssessment(
        stats.node_id,
        tuple(checks),
        additional,
        current_blocks >= MIN_INITIAL_BLOCKS,
        budget_exhausted,
        budget_capped,
        stats.transport_failed,
    )


@dataclass(frozen=True)
class AcquisitionCandidate:
    node_id: str
    predicted_loss_reduction: float
    predicted_cpu_seconds: float
    kind: str = "new_node"
    blocks: int = MIN_INITIAL_BLOCKS

    @property
    def utility(self) -> float:
        if self.predicted_cpu_seconds <= 0.0:
            return math.inf if self.predicted_loss_reduction > 0.0 else 0.0
        return max(0.0, self.predicted_loss_reduction) / self.predicted_cpu_seconds


def predict_candidate_utilities(
    training_q: Sequence[Sequence[float]],
    training_seconds_per_block: Sequence[float],
    validation_q: Sequence[Sequence[float]],
    validation_losses: Sequence[float],
    candidates: Mapping[str, Sequence[float]],
    *,
    q_min: Sequence[float],
    q_max: Sequence[float],
    initial_blocks: int = MIN_INITIAL_BLOCKS,
) -> dict[str, dict[str, float]]:
    """Predict validation-loss reduction per CPU second for unused nodes.

    This is a deterministic, inexpensive active-learning approximation to the
    conditional-variance reduction used in computer-experiment design. A
    Matérn-3/2 kernel transfers measured validation loss to nearby candidates;
    the candidate's novelty relative to existing training nodes supplies the
    conditional-variance factor. CPU cost is an inverse-distance estimate from
    the eight nearest completed training nodes. No transport observable is
    predicted or altered by this scheduler.
    """

    if not training_q or len(training_q) != len(training_seconds_per_block):
        raise ValueError("utility prediction needs CPU timing for every training node")
    if not validation_q or len(validation_q) != len(validation_losses):
        raise ValueError("utility prediction needs a loss for every validation node")
    if initial_blocks <= 0 or len(q_min) != 4 or len(q_max) != 4:
        raise ValueError("utility prediction has invalid domain or block count")
    span = [float(high) - float(low) for low, high in zip(q_min, q_max)]
    if any(not math.isfinite(value) or value <= 0.0 for value in span):
        raise ValueError("utility prediction q domain must have four positive spans")

    def normalized(values: Sequence[float]) -> tuple[float, float, float, float]:
        if len(values) != 4:
            raise ValueError("utility prediction q vectors must be four-dimensional")
        result = tuple(
            (float(value) - float(low)) / width
            for value, low, width in zip(values, q_min, span)
        )
        if any(not math.isfinite(value) for value in result):
            raise ValueError("utility prediction q vectors must be finite")
        return result  # type: ignore[return-value]

    train = [normalized(value) for value in training_q]
    validation = [normalized(value) for value in validation_q]
    losses = [float(value) for value in validation_losses]
    timings = [float(value) for value in training_seconds_per_block]
    if any(not math.isfinite(value) or value < 0.0 for value in losses):
        raise ValueError("validation losses must be finite and nonnegative")
    if any(not math.isfinite(value) or value <= 0.0 for value in timings):
        raise ValueError("training CPU seconds per block must be finite and positive")

    def distance(a: Sequence[float], b: Sequence[float]) -> float:
        return math.sqrt(sum((left - right) ** 2 for left, right in zip(a, b)))

    nearest_validation_distances = [
        min(distance(point, existing) for existing in train) for point in validation
    ]
    ordered = sorted(nearest_validation_distances)
    middle = len(ordered) // 2
    median = (
        ordered[middle]
        if len(ordered) % 2
        else 0.5 * (ordered[middle - 1] + ordered[middle])
    )
    length_scale = min(0.5, max(0.05, 1.5 * median))
    root_three = math.sqrt(3.0)

    def kernel(separation: float) -> float:
        scaled = root_three * separation / length_scale
        return (1.0 + scaled) * math.exp(-scaled)

    result: dict[str, dict[str, float]] = {}
    for node_id, raw_q in candidates.items():
        point = normalized(raw_q)
        neighbors = heapq.nsmallest(
            min(8, len(train)),
            (
                (distance(point, existing), timing)
                for existing, timing in zip(train, timings)
            ),
        )
        nearest = neighbors[0][0]
        novelty = max(0.0, 1.0 - kernel(nearest) ** 2)
        local_loss = sum(
            loss * kernel(distance(point, validation_point)) ** 2
            for loss, validation_point in zip(losses, validation)
        ) / len(validation)
        if neighbors[0][0] <= 1.0e-14:
            seconds_per_block = neighbors[0][1]
        else:
            weights = [1.0 / max(separation, 1.0e-8) for separation, _ in neighbors]
            seconds_per_block = sum(
                weight * timing for weight, (_, timing) in zip(weights, neighbors)
            ) / sum(weights)
        predicted_cpu = initial_blocks * seconds_per_block
        predicted_reduction = novelty * local_loss
        result[str(node_id)] = {
            "predicted_loss_reduction": predicted_reduction,
            "predicted_cpu_seconds": predicted_cpu,
            "utility": predicted_reduction / predicted_cpu,
            "local_validation_loss": local_loss,
            "novelty": novelty,
            "length_scale": length_scale,
        }
    return result


@dataclass(frozen=True)
class AllocationAction:
    node_id: str
    kind: str
    blocks: int
    utility: float
    reason: str
    predicted_loss_reduction: float = 0.0
    predicted_cpu_seconds: float = 0.0


@dataclass(frozen=True)
class BranchUtility:
    action_count: int
    predicted_loss_reduction: float
    predicted_cpu_seconds: float

    @property
    def utility(self) -> float:
        return (
            self.predicted_loss_reduction / self.predicted_cpu_seconds
            if self.predicted_cpu_seconds > 0.0
            else 0.0
        )


@dataclass(frozen=True)
class AdaptiveBatchDecision:
    selected_branch: str
    continuation: BranchUtility
    new_node_batch: BranchUtility
    forced_continuation: bool


def continuation_block_count(current_blocks: int, recommended_blocks: int) -> int:
    """Blocks one scheduling pass gives an existing node.

    A node first completes its MIN_INITIAL_BLOCKS baseline.  After that a pass
    at most doubles its histories, and never issues more than its projected
    need (``recommended_blocks``, already capped by the role's history cap).
    The noise estimate is refreshed between these geometric steps, so a noisy
    pilot never commits a long run, while a node that needs N times its
    histories reaches them in about log2(N) passes instead of N.
    """

    if current_blocks < 0 or recommended_blocks < 0:
        raise ValueError("block counts must be non-negative")
    if current_blocks < MIN_INITIAL_BLOCKS:
        return min(recommended_blocks, MIN_INITIAL_BLOCKS - current_blocks)
    return min(recommended_blocks, current_blocks)


@dataclass(frozen=True)
class ContinuationLossReduction:
    """Explicit model-residual and mandatory-precision utility components."""

    local_model_loss: float
    precision_guard_loss: float
    noise_fraction: float
    blocks: int
    noise_reduction: float

    @property
    def predicted_model_loss_reduction(self) -> float:
        return (
            self.local_model_loss
            * self.noise_fraction
            * self.noise_reduction
        )

    @property
    def predicted_precision_guard_reduction(self) -> float:
        return (
            self.precision_guard_loss
            * self.noise_fraction
            * self.noise_reduction
        )

    @property
    def total(self) -> float:
        return (
            self.predicted_model_loss_reduction
            + self.predicted_precision_guard_reduction
        )

    def as_dict(self) -> dict[str, float]:
        return {
            "local_model_loss": self.local_model_loss,
            "precision_guard_loss": self.precision_guard_loss,
            "noise_fraction": self.noise_fraction,
            "blocks": self.blocks,
            "noise_reduction": self.noise_reduction,
            "predicted_model_loss_reduction": self.predicted_model_loss_reduction,
            "predicted_precision_guard_reduction": (
                self.predicted_precision_guard_reduction
            ),
            "predicted_total_loss_reduction": self.total,
        }


def compare_adaptive_batches(
    continuations: Sequence[tuple[float, float]],
    new_nodes: Sequence[tuple[float, float]],
    *,
    forced_continuation: bool = False,
) -> AdaptiveBatchDecision:
    """Choose histories or nodes using one summed loss-reduction/CPU currency."""

    def aggregate(values: Sequence[tuple[float, float]]) -> BranchUtility:
        reduction = 0.0
        cpu = 0.0
        for raw_reduction, raw_cpu in values:
            item_reduction = float(raw_reduction)
            item_cpu = float(raw_cpu)
            if not (
                math.isfinite(item_reduction)
                and item_reduction >= 0.0
                and math.isfinite(item_cpu)
                and item_cpu > 0.0
            ):
                raise ValueError("adaptive branch actions need finite loss reduction and CPU")
            reduction += item_reduction
            cpu += item_cpu
        return BranchUtility(len(values), reduction, cpu)

    continuation = aggregate(continuations)
    new_batch = aggregate(new_nodes)
    choose_new = bool(
        not forced_continuation
        and new_batch.action_count
        and (
            not continuation.action_count
            or new_batch.utility > continuation.utility
        )
    )
    return AdaptiveBatchDecision(
        "new_node_batch" if choose_new else "continuation",
        continuation,
        new_batch,
        bool(forced_continuation),
    )


def continuation_loss_reduction_components(
    local_validation_loss: float,
    assessment: NodeAssessment,
    current_blocks: int,
    blocks: int = 1,
) -> ContinuationLossReduction:
    """Predict ``blocks`` more blocks' reduction in the loss used for new nodes.

    The validation objectives are amplitude-like errors (relative error and
    total variation), so their MC-noise part scales as ``N^-1/2``.  The worst
    stopping-ratio estimates how much of the local validation loss can still be
    attributed to noise above its allowed target.  Multiplying the two gives a
    conservative loss reduction in the same units as a new-node model-loss
    reduction, rather than a bare uncertainty overshoot.
    """

    if (
        not math.isfinite(local_validation_loss)
        or local_validation_loss < 0.0
        or current_blocks <= 0
        or blocks <= 0
    ):
        raise ValueError("continuation utility inputs must be finite and positive")
    worst = assessment.worst_ratio
    if not math.isfinite(worst):
        noise_fraction = 1.0
    elif worst <= 1.0:
        noise_fraction = 0.0
    else:
        noise_fraction = min(1.0, max(0.0, 1.0 - 1.0 / worst))
    noise_reduction = 1.0 - math.sqrt(current_blocks / (current_blocks + blocks))
    # Keep the measured model residual and the mandatory precision guard
    # explicit.  Below the model acceptance threshold the guard supplies only
    # the missing fraction up to one threshold-normalized loss unit; it does not
    # masquerade as observed model error in the audit record.
    return ContinuationLossReduction(
        local_model_loss=local_validation_loss,
        precision_guard_loss=max(0.0, 1.0 - local_validation_loss),
        noise_fraction=noise_fraction,
        blocks=blocks,
        noise_reduction=noise_reduction,
    )


def continuation_loss_reduction(
    local_validation_loss: float,
    assessment: NodeAssessment,
    current_blocks: int,
    blocks: int = 1,
) -> float:
    return continuation_loss_reduction_components(
        local_validation_loss, assessment, current_blocks, blocks
    ).total


def allocate_actions(
    existing: Iterable[NodeStatistics],
    new_candidates: Iterable[AcquisitionCandidate],
    *,
    max_actions: int = 64,
    criteria: AllocationCriteria = AllocationCriteria(),
    default_seconds_per_block: float = 1.0,
    continuation_local_losses: Mapping[str, float] | None = None,
) -> tuple[AllocationAction, ...]:
    """Rank continuation and new-node work by expected loss reduction/CPU.

    Existing nodes receive at most a doubling per scheduling pass
    (continuation_block_count), so observations update the noise estimate
    before more transport is committed.
    """

    if max_actions < 0 or default_seconds_per_block <= 0.0:
        raise ValueError("invalid scheduler limits")
    ranked: list[AllocationAction] = []
    local_losses = continuation_local_losses or {}
    for stats in existing:
        assessment = assess_node(stats, criteria)
        if (
            assessment.passed
            or assessment.history_budget_capped
            or assessment.recommended_additional_blocks <= 0
        ):
            continue
        observed_blocks = stats.completed_blocks or (
            stats.histories // BLOCK_HISTORIES
        )
        seconds_per_block = (
            stats.wall_seconds / observed_blocks
            if observed_blocks > 0 and stats.wall_seconds > 0.0
            else default_seconds_per_block
        )
        if stats.node_id not in local_losses:
            raise ValueError(
                f"continuation node {stats.node_id!r} lacks a predicted local validation loss"
            )
        blocks = continuation_block_count(
            observed_blocks, assessment.recommended_additional_blocks
        )
        reduction = continuation_loss_reduction(
            float(local_losses[stats.node_id]), assessment, observed_blocks, blocks
        )
        cpu_seconds = blocks * seconds_per_block
        ranked.append(
            AllocationAction(
                stats.node_id,
                "continue",
                blocks,
                reduction / cpu_seconds,
                (
                    f"{blocks}-block MC-noise reduction in threshold-normalized "
                    f"validation loss; worst uncertainty ratio {assessment.worst_ratio:.6g}"
                ),
                reduction,
                cpu_seconds,
            )
        )
    for candidate in new_candidates:
        if candidate.kind != "new_node" or candidate.blocks <= 0:
            raise ValueError("new candidate must request positive new_node blocks")
        ranked.append(
            AllocationAction(
                candidate.node_id,
                candidate.kind,
                candidate.blocks,
                candidate.utility,
                "predicted validation-loss reduction per CPU second",
                max(0.0, candidate.predicted_loss_reduction),
                candidate.predicted_cpu_seconds,
            )
        )
    ranked.sort(key=lambda item: (-item.utility, item.kind, item.node_id))
    return tuple(ranked[:max_actions])


@dataclass
class ValidationRoundState:
    """Stop only after every gate passes in two consecutive rounds."""

    consecutive_passes: int = 0
    history: list[bool] = field(default_factory=list)

    def update(self, all_validation_gates_pass: bool) -> bool:
        self.history.append(bool(all_validation_gates_pass))
        self.consecutive_passes = self.consecutive_passes + 1 if all_validation_gates_pass else 0
        return self.consecutive_passes >= 2
