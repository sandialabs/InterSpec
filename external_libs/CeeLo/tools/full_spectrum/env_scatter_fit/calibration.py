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

"""Repeated-null calibration of uncertainty estimators in CEEFIT archives.

The aggregate node arrays in a fit archive are deliberately *not* used here.
Every comparison is made between two immutable worker-block records from the
same node.  Blocks are sorted by numeric ID and consumed once, in adjacent
pairs, which makes the result independent of archive insertion/merge order.

The complete uncertainty record is checked:

* the raw per-history ``S2`` of total scattered fluence and every stored
  continuum-component, total-continuum, line-component, and total-line cell;
* the exact outer moments of the 16-element stopping vector; and
* covariance estimated from all four 32-replica map families after
  deterministic component/line yields, coarse registered-energy rebinnings,
  angular moments, and smooth response-like folds.

The scalar total and stopping-vector omnibus tests use central-limit chi-square
diagnostics; cell marginals retain explicit standardized-residual coverage
without pretending that unrecorded cell covariances are independent. Replica
covariance is estimated rather than known, so it is checked with the
finite-sample two-sample Hotelling statistic.  The report records coverage as
well as the
explicit gates; a report cannot pass with fewer than eight disjoint
comparisons.  Physically supported component/line segments must be represented
by stochastic observations rather than being silently declared calibrated
from empty cells.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

from .archive import RETAINED_RAW_BLOCKS, _fold_generation, require_hdf5
from .publication import archive_snapshot_sha256


CALIBRATION_SCHEMA = "CEEFIT-UNCERTAINTY-CALIBRATION-v1"
PROJECTION_SCHEMA = "CEEFIT-REPLICA-CALIBRATION-PROJECTIONS-v1"
MINIMUM_COMPARISONS = 8
OMNIBUS_TWO_SIDED_ALPHA = 1.0e-3
INDIVIDUAL_REPLICA_TWO_SIDED_ALPHA = 1.0e-4
MAX_ABS_STANDARDIZED_RESIDUAL = 5.0
_PSD_RELATIVE_TOLERANCE = 1.0e-10
# Replica maps intentionally use float32.  The line angular projections can
# combine 64 signed mu cells, so their worst-case closure roundoff is bounded by
# 64 float32 ulps.  This is about 8 ppm, far below the 0.25% MC-noise target.
_FLOAT32_REPLICA_NULL_TOLERANCE = 64.0 * 2.0**-23
_REQUIRED_MOMENT_TENSORS = (
    "continuum_components",
    "continuum_total",
    "line_components",
    "line_total",
)
_MINIMUM_WITHIN_1_96_FRACTION = 0.70
_MINIMUM_WITHIN_3_FRACTION = 0.90


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _atomic_json(path: str | Path, value: Mapping[str, Any]) -> None:
    """Write strict JSON in the destination directory, then atomically replace."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode(
        "ascii"
    )
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=destination.name + ".",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            temporary.write(encoded)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                Path(temporary_name).unlink()
            except FileNotFoundError:
                pass


def _regularized_gamma_p(shape: float, value: float) -> float:
    """Regularized lower incomplete gamma, using stable standard expansions."""

    if not math.isfinite(shape) or shape <= 0.0 or value < 0.0:
        raise ValueError("invalid regularized-gamma arguments")
    if value == 0.0:
        return 0.0
    if math.isinf(value):
        return 1.0
    log_scale = -value + shape * math.log(value) - math.lgamma(shape)
    epsilon = 2.0e-15
    tiny = 1.0e-300
    if value < shape + 1.0:
        term = 1.0 / shape
        total = term
        denominator = shape
        for _ in range(10_000):
            denominator += 1.0
            term *= value / denominator
            total += term
            if abs(term) <= abs(total) * epsilon:
                return min(1.0, max(0.0, total * math.exp(log_scale)))
        raise RuntimeError("regularized-gamma series did not converge")

    b = value + 1.0 - shape
    c = 1.0 / tiny
    d = 1.0 / max(abs(b), tiny)
    if b < 0.0:
        d = -d
    fraction = d
    for iteration in range(1, 10_001):
        coefficient = -iteration * (iteration - shape)
        b += 2.0
        d = coefficient * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + coefficient / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) <= epsilon:
            upper = math.exp(log_scale) * fraction
            return min(1.0, max(0.0, 1.0 - upper))
    raise RuntimeError("regularized-gamma continued fraction did not converge")


def _beta_continued_fraction(a: float, b: float, value: float) -> float:
    epsilon = 2.0e-15
    tiny = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * value / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    result = d
    for iteration in range(1, 10_001):
        twice = 2 * iteration
        coefficient = iteration * (b - iteration) * value / (
            (qam + twice) * (a + twice)
        )
        d = 1.0 + coefficient * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + coefficient / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        result *= d * c

        coefficient = -(a + iteration) * (qab + iteration) * value / (
            (a + twice) * (qap + twice)
        )
        d = 1.0 + coefficient * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + coefficient / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) <= epsilon:
            return result
    raise RuntimeError("regularized-beta continued fraction did not converge")


def _regularized_beta(value: float, a: float, b: float) -> float:
    if not (a > 0.0 and b > 0.0 and 0.0 <= value <= 1.0):
        raise ValueError("invalid regularized-beta arguments")
    if value == 0.0:
        return 0.0
    if value == 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(value)
        + b * math.log1p(-value)
    )
    if value < (a + 1.0) / (a + b + 2.0):
        result = front * _beta_continued_fraction(a, b, value) / a
    else:
        result = 1.0 - front * _beta_continued_fraction(b, a, 1.0 - value) / b
    return min(1.0, max(0.0, result))


def _chi_square_cdf(value: float, degrees_of_freedom: int) -> float:
    if degrees_of_freedom <= 0 or not math.isfinite(value) or value < 0.0:
        raise ValueError("invalid chi-square statistic")
    return _regularized_gamma_p(0.5 * degrees_of_freedom, 0.5 * value)


def _f_cdf(value: float, numerator_df: int, denominator_df: int) -> float:
    if (
        numerator_df <= 0
        or denominator_df <= 0
        or not math.isfinite(value)
        or value < 0.0
    ):
        raise ValueError("invalid F statistic")
    transformed = numerator_df * value / (
        numerator_df * value + denominator_df
    )
    return _regularized_beta(
        transformed, 0.5 * numerator_df, 0.5 * denominator_df
    )


def _two_sided_probability(cdf: float) -> float:
    return min(1.0, max(0.0, 2.0 * min(cdf, 1.0 - cdf)))


def _wilson_interval(successes: int, trials: int) -> list[float] | None:
    if trials <= 0:
        return None
    z = 1.959963984540054
    fraction = successes / trials
    scale = 1.0 + z * z / trials
    center = (fraction + z * z / (2.0 * trials)) / scale
    half = z * math.sqrt(
        fraction * (1.0 - fraction) / trials + z * z / (4.0 * trials * trials)
    ) / scale
    return [max(0.0, center - half), min(1.0, center + half)]


def _new_coverage_state() -> dict[str, Any]:
    return {
        "count": 0,
        "within_1sigma": 0,
        "within_1_96sigma": 0,
        "within_3sigma": 0,
        "max_abs_standardized_residual": 0.0,
    }


def _update_coverage_state(
    state: dict[str, Any], values: Any, numpy: Any | None = None
) -> None:
    if numpy is not None:
        array = numpy.asarray(values, dtype="f8").reshape(-1)
        if not numpy.all(numpy.isfinite(array)):
            raise ValueError("coverage contains a non-finite standardized residual")
        absolute = numpy.abs(array)
        state["count"] += int(array.size)
        state["within_1sigma"] += int(numpy.count_nonzero(absolute <= 1.0))
        state["within_1_96sigma"] += int(
            numpy.count_nonzero(absolute <= 1.959963984540054)
        )
        state["within_3sigma"] += int(numpy.count_nonzero(absolute <= 3.0))
        state["max_abs_standardized_residual"] = max(
            float(state["max_abs_standardized_residual"]),
            float(numpy.max(absolute)) if absolute.size else 0.0,
        )
        return

    for value in values:
        normalized = float(value)
        if not math.isfinite(normalized):
            raise ValueError("coverage contains a non-finite standardized residual")
        absolute = abs(normalized)
        state["count"] += 1
        state["within_1sigma"] += int(absolute <= 1.0)
        state["within_1_96sigma"] += int(absolute <= 1.959963984540054)
        state["within_3sigma"] += int(absolute <= 3.0)
        state["max_abs_standardized_residual"] = max(
            float(state["max_abs_standardized_residual"]), absolute
        )


def _coverage_from_state(state: Mapping[str, Any]) -> dict[str, Any]:
    trials = int(state["count"])
    result: dict[str, Any] = {
        "count": trials,
        "expected_fractions": {
            "within_1sigma": 0.6826894921370859,
            "within_1_96sigma": 0.95,
            "within_3sigma": 0.9973002039367398,
        },
    }
    for label in ("within_1sigma", "within_1_96sigma", "within_3sigma"):
        successes = int(state[label])
        result[label] = {
            "count": successes,
            "fraction": successes / trials if trials else 0.0,
            "wilson_95_interval": _wilson_interval(successes, trials),
        }
    result["max_abs_standardized_residual"] = float(
        state["max_abs_standardized_residual"]
    )
    return result


def _coverage(values: Sequence[float]) -> dict[str, Any]:
    state = _new_coverage_state()
    _update_coverage_state(state, values)
    return _coverage_from_state(state)


def _pair_inventory(root: Any) -> tuple[list[tuple[str, str, str]], list[dict[str, Any]]]:
    pairs: list[tuple[str, str, str]] = []
    unused: list[dict[str, Any]] = []
    block_root = root.get("blocks")
    if block_root is None:
        return pairs, unused
    for node_name in sorted(block_root, key=lambda value: int(value, 16)):
        names = sorted(block_root[node_name], key=int)
        if _fold_generation(root, node_name) is not None:
            # A folded node's calibration records are its retained blocks;
            # any other raw record is transient (held back by a gap, or
            # covered by the fold and not yet removed).
            names = [name for name in names if int(name) <= RETAINED_RAW_BLOCKS]
        for index in range(0, len(names) - 1, 2):
            pairs.append((node_name.lower(), names[index], names[index + 1]))
        if len(names) % 2:
            unused.append({"node_id": node_name.lower(), "block_id": int(names[-1])})
    return pairs, unused


def _comparison_label(pair: tuple[str, str, str]) -> dict[str, Any]:
    return {
        "node_id": pair[0],
        "block_ids": [int(pair[1]), int(pair[2])],
    }


def _independence_diagnostics(
    root: Any, pairs: Sequence[tuple[str, str, str]], numpy: Any
) -> dict[str, Any]:
    """Recheck disjoint history ranges and streams instead of trusting group IDs."""

    failures: list[dict[str, Any]] = []
    ranges_by_node: dict[str, list[tuple[int, int, int]]] = {}
    used_seeds: set[int] = set()
    for pair_index, pair in enumerate(pairs):
        label = _comparison_label(pair)
        try:
            for block_name in pair[1:]:
                group = root[f"blocks/{pair[0]}/{block_name}"]
                status = str(group.attrs.get("status", ""))
                failures_count = int(group.attrs.get("transport_failures", -1))
                histories = int(group.attrs.get("histories", -1))
                if status not in {"complete", "converged", "capped"}:
                    raise ValueError(f"block {block_name} has unusable status {status!r}")
                if failures_count != 0:
                    raise ValueError(f"block {block_name} has transport failures")
                ids = numpy.asarray(group["id"][...], dtype="u8")
                begins = numpy.asarray(group["history_begin"][...], dtype="u8")
                counts = numpy.asarray(group["history_count"][...], dtype="u8")
                seeds = numpy.asarray(group["rng_seed"][...], dtype="u8")
                if any(values.shape != (32,) for values in (ids, begins, counts, seeds)):
                    raise ValueError(f"block {block_name} lacks exactly 32 replica records")
                if not numpy.array_equal(ids, numpy.arange(32, dtype="u8")):
                    raise ValueError(f"block {block_name} has noncanonical replica IDs")
                if numpy.any(counts == 0) or int(numpy.sum(counts, dtype="u8")) != histories:
                    raise ValueError(f"block {block_name} has invalid replica history accounting")

                node_ranges = ranges_by_node.setdefault(pair[0], [])
                for begin_value, count_value, seed_value in zip(begins, counts, seeds):
                    begin = int(begin_value)
                    count = int(count_value)
                    end = begin + count
                    if end <= begin:
                        raise ValueError(f"block {block_name} has an overflowing history range")
                    if any(begin < old_end and old_begin < end for old_begin, old_end, _ in node_ranges):
                        raise ValueError(
                            f"block {block_name} overlaps a history range used by another comparison"
                        )
                    seed = int(seed_value)
                    if seed in used_seeds:
                        raise ValueError(
                            f"block {block_name} reuses an RNG seed from another comparison"
                        )
                    node_ranges.append((begin, end, pair_index))
                    used_seeds.add(seed)
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            failures.append({**label, "reason": str(error)})
    return {
        "passed": not failures,
        "checked_comparisons": len(pairs),
        "checked_blocks": 2 * len(pairs),
        "failures": failures,
    }


def _total_moments(group: Any) -> tuple[int, float, float]:
    try:
        count = int(group.attrs["scattered_total_N"])
        s1 = float(group.attrs["scattered_total_S1"])
        s2 = float(group.attrs["scattered_total_S2"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("block lacks total scattered-fluence moments") from error
    if count < 2 or not (math.isfinite(s1) and math.isfinite(s2)):
        raise ValueError("total moments require finite values and N >= 2")
    centered = s2 - s1 * s1 / count
    tolerance = 64.0 * math.ulp(1.0) * max(
        float.fromhex("0x1.0p-1022"), abs(s2), abs(s1 * s1 / count)
    )
    if centered < -tolerance:
        raise ValueError("total S2 is smaller than S1^2/N")
    centered = max(0.0, centered)
    return count, s1 / count, centered / (count * (count - 1))


def _s2_diagnostics(root: Any, pairs: Sequence[tuple[str, str, str]]) -> dict[str, Any]:
    residuals: list[float] = []
    comparisons: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for pair in pairs:
        label = _comparison_label(pair)
        try:
            left = _total_moments(root[f"blocks/{pair[0]}/{pair[1]}"])
            right = _total_moments(root[f"blocks/{pair[0]}/{pair[2]}"])
            predicted_variance = left[2] + right[2]
            difference = left[1] - right[1]
            if not math.isfinite(predicted_variance) or predicted_variance <= 0.0:
                raise ValueError("total repeated-null comparison has no positive variance")
            z = difference / math.sqrt(predicted_variance)
            if not math.isfinite(z):
                raise ValueError("total standardized residual is non-finite")
            residuals.append(z)
            comparisons.append(
                {
                    **label,
                    "difference": difference,
                    "predicted_standard_error": math.sqrt(predicted_variance),
                    "standardized_residual": z,
                }
            )
        except (KeyError, ValueError, FloatingPointError) as error:
            failures.append({**label, "reason": str(error)})

    statistic = math.fsum(value * value for value in residuals)
    degrees = len(residuals)
    omnibus = (
        _two_sided_probability(_chi_square_cdf(statistic, degrees)) if degrees else 0.0
    )
    coverage = _coverage(residuals)
    passed = bool(
        len(pairs) >= MINIMUM_COMPARISONS
        and degrees == len(pairs)
        and not failures
        and omnibus >= OMNIBUS_TWO_SIDED_ALPHA
        and coverage["max_abs_standardized_residual"]
        <= MAX_ABS_STANDARDIZED_RESIDUAL
    )
    return {
        "passed": passed,
        "comparison_count": len(comparisons),
        "chi_square": statistic,
        "degrees_of_freedom": degrees,
        "omnibus_two_sided_probability": omnibus,
        "coverage": coverage,
        "comparisons": comparisons,
        "failures": failures,
    }


def _json_attribute(group: Any, name: str, label: str) -> Any:
    try:
        value = group.attrs[name]
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        return json.loads(str(value))
    except (KeyError, TypeError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"archive lacks valid {label}") from error


def _moment_segments(
    root: Any, tensor_name: str, cell_count: int
) -> list[dict[str, Any]]:
    """Return the fixed logical segments for one flattened raw moment tensor."""

    components = _json_attribute(
        root["components"], "descriptors_json", "component descriptors"
    )
    lines = _json_attribute(
        root["axes"], "line_definitions_json", "line definitions"
    )
    n_mu = len(root["axes/mu_edges"]) - 1
    if not isinstance(components, list) or not components or n_mu <= 0:
        raise ValueError("archive component/mu inventory is empty")
    if not isinstance(lines, list) or not lines:
        raise ValueError("archive line inventory is empty")

    result: list[dict[str, Any]] = []
    if tensor_name == "continuum_components":
        if cell_count % len(components):
            raise ValueError("continuum-component moment shape is not component-major")
        stride = cell_count // len(components)
        if stride % n_mu:
            raise ValueError("continuum-component cells do not match the mu axis")
        for index, component in enumerate(components):
            result.append(
                {
                    "name": str(component["name"]),
                    "begin": index * stride,
                    "end": (index + 1) * stride,
                    "component_index": index,
                    "kind": "continuum_component",
                    "globally_required": bool(
                        int(component["origin"]) == 2
                        or int(component["post_creation_order"]) != 0
                    ),
                }
            )
    elif tensor_name == "continuum_total":
        if cell_count % n_mu:
            raise ValueError("total-continuum cells do not match the mu axis")
        result.append(
            {
                "name": "total_continuum",
                "begin": 0,
                "end": cell_count,
                "kind": "continuum_total",
                "globally_required": True,
            }
        )
    elif tensor_name == "line_components":
        divisor = len(components) * len(lines)
        if cell_count != divisor * n_mu:
            raise ValueError("line-component moment shape is not component/line/mu major")
        for component_index, component in enumerate(components):
            for line_index, line in enumerate(lines):
                offset = (component_index * len(lines) + line_index) * n_mu
                origin_matches = int(component["origin"]) == int(line["origin"])
                component_order = int(component["post_creation_order"])
                # The next-event score is made at the interaction that sends a
                # photon toward the detector.  Consequently an all-Rayleigh
                # primary line first appears in order 1, never order 0.  A
                # 511-keV photon also has an order-0 birth score and may retain
                # its exact energy through later Rayleigh interactions.
                physically_assigned = bool(
                    origin_matches
                    and (
                        (int(line["origin"]) == 0 and component_order >= 1)
                        or int(line["origin"]) == 1
                    )
                )
                result.append(
                    {
                        "name": f"{component['name']}:{line['name']}",
                        "begin": offset,
                        "end": offset + n_mu,
                        "component_index": component_index,
                        "line_index": line_index,
                        "kind": "line_component",
                        "physically_assigned": physically_assigned,
                        # Component-line requirements are finalized from
                        # observed nonzero support below; a static order-0
                        # predicate would make this gate impossible for
                        # genuine elastic-E0 transport records.
                        "globally_required": False,
                    }
                )
    elif tensor_name == "line_total":
        if cell_count != len(lines) * n_mu:
            raise ValueError("total-line moment shape is not line/mu major")
        for line_index, line in enumerate(lines):
            result.append(
                {
                    "name": str(line["name"]),
                    "begin": line_index * n_mu,
                    "end": (line_index + 1) * n_mu,
                    "line_index": line_index,
                    "kind": "line_total",
                    "globally_required": True,
                }
            )
    else:  # pragma: no cover - all callers use the fixed v1 inventory
        raise ValueError(f"unknown moment tensor {tensor_name!r}")
    return result


def _segment_supported(group: Any, segment: Mapping[str, Any], numpy: Any) -> bool:
    kind = str(segment["kind"])
    if kind == "continuum_total":
        return True
    if kind == "continuum_component":
        offsets = numpy.asarray(group["component_support_offsets"][...], dtype="u8")
        support = numpy.asarray(group["component_support_mask"][...], dtype="u1")
        index = int(segment["component_index"])
        if offsets.ndim != 1 or index + 1 >= len(offsets):
            raise ValueError("raw block has invalid component-support offsets")
        begin, end = int(offsets[index]), int(offsets[index + 1])
        return bool(begin < end and numpy.any(support[begin:end] != 0))
    support = numpy.asarray(group["line_support_mask"][...], dtype="u1")
    line_index = int(segment["line_index"])
    if support.ndim != 1 or line_index >= len(support):
        raise ValueError("raw block has invalid exact-line support mask")
    line_supported = bool(support[line_index])
    if kind == "line_component":
        return line_supported and bool(segment["physically_assigned"])
    return line_supported


def _cell_moments(group: Any, tensor_name: str, numpy: Any) -> tuple[Any, Any]:
    """Return per-history means and their standard errors.

    Keep the uncertainty as an SE rather than an SE squared.  For the most
    strongly attenuated cells, the centered S2 is a representable subnormal
    while division by N(N-1) would underflow.  Taking its square root first
    retains the statistically meaningful scale.
    """
    try:
        moments = group[f"moments/{tensor_name}"]
        count = numpy.asarray(moments["N"][...], dtype="f8")
        s1 = numpy.asarray(moments["S1"][...], dtype="f8")
        s2 = numpy.asarray(moments["S2"][...], dtype="f8")
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"block lacks {tensor_name} cell moments") from error
    if count.ndim != 1 or s1.shape != count.shape or s2.shape != count.shape:
        raise ValueError(f"{tensor_name} cell moments have incompatible shapes")
    if count.size == 0 or numpy.any(count < 2.0):
        raise ValueError(f"{tensor_name} cell moments require N >= 2")
    if not (
        numpy.all(numpy.isfinite(count))
        and numpy.all(numpy.isfinite(s1))
        and numpy.all(numpy.isfinite(s2))
        and numpy.all(s1 >= 0.0)
        and numpy.all(s2 >= 0.0)
    ):
        raise ValueError(f"{tensor_name} cell moments contain invalid values")
    mean_square_sum = s1 * s1 / count
    centered = s2 - mean_square_sum
    # These moments are incident fluences and routinely live at 1e-8--1e-16
    # scales.  A unit absolute floor would classify every production-scale
    # variance as roundoff.  Bound cancellation by the two operands of the
    # subtraction itself; exact zero operands consequently have zero tolerance.
    tolerance = 128.0 * numpy.finfo("f8").eps * numpy.maximum(
        numpy.abs(s2), numpy.abs(mean_square_sum)
    )
    if numpy.any(centered < -tolerance):
        raise ValueError(f"{tensor_name} S2 is smaller than S1^2/N")
    centered = numpy.where(centered <= tolerance, 0.0, centered)
    standard_error = numpy.sqrt(centered) / (
        numpy.sqrt(count) * numpy.sqrt(count - 1.0)
    )
    return s1 / count, standard_error


def _cell_s2_diagnostics(
    root: Any, pairs: Sequence[tuple[str, str, str]], numpy: Any
) -> dict[str, Any]:
    """Calibrate every stored cell marginal without inventing cell covariance."""

    tensor_reports: dict[str, Any] = {}
    overall_coverage_state = _new_coverage_state()
    all_failures: list[dict[str, Any]] = []
    for tensor_name in _REQUIRED_MOMENT_TENSORS:
        coverage_state = _new_coverage_state()
        failures: list[dict[str, Any]] = []
        comparisons: list[dict[str, Any]] = []
        segment_accumulators: dict[str, dict[str, Any]] | None = None
        for pair in pairs:
            label = _comparison_label(pair)
            try:
                left_group = root[f"blocks/{pair[0]}/{pair[1]}"]
                right_group = root[f"blocks/{pair[0]}/{pair[2]}"]
                left_mean, left_standard_error = _cell_moments(
                    left_group, tensor_name, numpy
                )
                right_mean, right_standard_error = _cell_moments(
                    right_group, tensor_name, numpy
                )
                if left_mean.shape != right_mean.shape:
                    raise ValueError(f"paired {tensor_name} cell shapes differ")
                segments = _moment_segments(root, tensor_name, left_mean.size)
                if segment_accumulators is None:
                    segment_accumulators = {
                        str(segment["name"]): {
                            "name": str(segment["name"]),
                            "stored_cells": int(segment["end"] - segment["begin"]),
                            "globally_required": bool(segment["globally_required"]),
                            "support_comparisons": 0,
                            "stochastic_comparisons": 0,
                            "nonzero_mean_comparisons": 0,
                            "stochastic_residual_count": 0,
                            "deterministic_cell_comparisons": 0,
                            "max_abs_standardized_residual": 0.0,
                        }
                        for segment in segments
                    }
                elif tuple(segment_accumulators) != tuple(
                    str(segment["name"]) for segment in segments
                ):
                    raise ValueError(f"{tensor_name} segment inventory changed between blocks")

                difference = left_mean - right_mean
                standard_error = numpy.hypot(
                    left_standard_error, right_standard_error
                )
                scale = numpy.maximum(numpy.abs(left_mean), numpy.abs(right_mean))
                # If W^2 is below the smallest representable subnormal, S2
                # cannot carry its variance even though S1 can still carry W.
                # This is an IEEE representation limit (~1e-162 in W), not a
                # physics precision cutoff.  Above it, every unsupported
                # difference remains a hard calibration failure.
                square_underflow = math.sqrt(
                    float(
                        numpy.nextafter(
                            numpy.float64(0.0), numpy.float64(1.0)
                        )
                    )
                )
                deterministic_tolerance = numpy.maximum(
                    256.0 * numpy.finfo("f8").eps * scale,
                    2.0 * square_underflow,
                )
                # Each cell carries its own raw marginal S2.  A relative cutoff
                # against the largest cell would incorrectly discard weak but
                # scientifically real components, so every strictly positive
                # finite variance is calibrated.
                stochastic = standard_error > 0.0
                if numpy.any(numpy.logical_and(~stochastic, numpy.abs(difference) > deterministic_tolerance)):
                    raise ValueError(
                        f"{tensor_name} has a repeated-null difference with no S2 support"
                    )
                z = difference[stochastic] / standard_error[stochastic]
                if not numpy.all(numpy.isfinite(z)):
                    raise ValueError(f"{tensor_name} has non-finite standardized residuals")
                _update_coverage_state(coverage_state, z, numpy)
                _update_coverage_state(overall_coverage_state, z, numpy)

                for segment in segments:
                    name = str(segment["name"])
                    accumulator = segment_accumulators[name]
                    begin, end = int(segment["begin"]), int(segment["end"])
                    supported = _segment_supported(left_group, segment, numpy) or _segment_supported(
                        right_group, segment, numpy
                    )
                    accumulator["support_comparisons"] += int(supported)
                    segment_scale = numpy.maximum.reduce(
                        (
                            numpy.ones(end - begin, dtype="f8"),
                            numpy.abs(left_mean[begin:end]),
                            numpy.abs(right_mean[begin:end]),
                        )
                    )
                    accumulator["nonzero_mean_comparisons"] += int(
                        numpy.any(
                            numpy.maximum(
                                numpy.abs(left_mean[begin:end]),
                                numpy.abs(right_mean[begin:end]),
                            )
                            > 256.0 * numpy.finfo("f8").eps * segment_scale
                        )
                    )
                    segment_stochastic = stochastic[begin:end]
                    accumulator["stochastic_comparisons"] += int(
                        numpy.any(segment_stochastic)
                    )
                    segment_z = (
                        difference[begin:end][segment_stochastic]
                        / standard_error[begin:end][segment_stochastic]
                    )
                    accumulator["stochastic_residual_count"] += int(
                        numpy.count_nonzero(segment_stochastic)
                    )
                    accumulator["deterministic_cell_comparisons"] += int(
                        (end - begin) - numpy.count_nonzero(segment_stochastic)
                    )
                    accumulator["max_abs_standardized_residual"] = max(
                        float(accumulator["max_abs_standardized_residual"]),
                        float(numpy.max(numpy.abs(segment_z)))
                        if segment_z.size
                        else 0.0,
                    )
                comparisons.append(
                    {
                        **label,
                        "stored_cell_count": int(left_mean.size),
                        "stochastic_residual_count": int(numpy.count_nonzero(stochastic)),
                        "deterministic_cell_count": int(left_mean.size - numpy.count_nonzero(stochastic)),
                        "max_abs_standardized_residual": (
                            float(numpy.max(numpy.abs(z))) if z.size else 0.0
                        ),
                    }
                )
            except (KeyError, TypeError, ValueError, FloatingPointError) as error:
                failures.append({**label, "reason": str(error)})

        segment_records = list((segment_accumulators or {}).values())
        for record in segment_records:
            required = bool(
                record["globally_required"]
                or int(record["nonzero_mean_comparisons"]) > 0
            )
            record["required"] = required
            record["passed"] = bool(
                not required
                or int(record["stochastic_comparisons"]) >= MINIMUM_COMPARISONS
            )
        coverage = _coverage_from_state(coverage_state)
        familywise_maximum_gate = max(
            MAX_ABS_STANDARDIZED_RESIDUAL,
            math.sqrt(
                2.0
                * math.log(
                    2.0
                    * max(1, int(coverage["count"]))
                    / OMNIBUS_TWO_SIDED_ALPHA
                )
            ),
        )
        fractions_pass = bool(
            coverage["count"] >= MINIMUM_COMPARISONS
            and coverage["within_1_96sigma"]["fraction"]
            >= _MINIMUM_WITHIN_1_96_FRACTION
            and coverage["within_3sigma"]["fraction"] >= _MINIMUM_WITHIN_3_FRACTION
        )
        passed = bool(
            len(pairs) >= MINIMUM_COMPARISONS
            and len(comparisons) == len(pairs)
            and not failures
            and segment_records
            and all(bool(record["passed"]) for record in segment_records)
            and fractions_pass
            and coverage["max_abs_standardized_residual"]
            <= familywise_maximum_gate
        )
        tensor_report = {
            "passed": passed,
            "comparison_count": len(comparisons),
            "coverage": coverage,
            "maximum_abs_standardized_residual_gate": familywise_maximum_gate,
            "segments": segment_records,
            "comparisons": comparisons,
            "failures": failures,
        }
        tensor_reports[tensor_name] = tensor_report
        all_failures.extend(
            {"tensor": tensor_name, **failure} for failure in failures
        )

    family_inventory_complete = set(tensor_reports) == set(_REQUIRED_MOMENT_TENSORS)
    overall_coverage = _coverage_from_state(overall_coverage_state)
    return {
        "passed": bool(
            family_inventory_complete
            and all(tensor_reports[name]["passed"] for name in _REQUIRED_MOMENT_TENSORS)
        ),
        "required_tensors": list(_REQUIRED_MOMENT_TENSORS),
        "family_inventory_complete": family_inventory_complete,
        "coverage": overall_coverage,
        "tensors": tensor_reports,
        "failures": all_failures,
    }


def _stopping_moments(group: Any, numpy: Any) -> tuple[Any, Any]:
    try:
        stopping = group["stopping"]
        count = int(stopping.attrs["N"])
        s1 = numpy.asarray(stopping["S1"][...], dtype="f8")
        outer = numpy.asarray(stopping["outer_S2"][...], dtype="f8")
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("block lacks exact stopping-vector moments") from error
    if count < 2 or s1.shape != (16,) or outer.shape != (16, 16):
        raise ValueError("stopping moments require N >= 2 and exact 16-vector shapes")
    if not (numpy.all(numpy.isfinite(s1)) and numpy.all(numpy.isfinite(outer))):
        raise ValueError("stopping moments contain non-finite values")
    asymmetry = float(numpy.max(numpy.abs(outer - outer.T)))
    scale = max(float(numpy.finfo("f8").tiny), float(numpy.max(numpy.abs(outer))))
    if asymmetry > 1.0e-11 * scale:
        raise ValueError("stopping outer moment is not symmetric")
    centered = 0.5 * (outer + outer.T) - numpy.outer(s1, s1) / count
    covariance = centered / (count * (count - 1))
    covariance = _positive_semidefinite(covariance, numpy, "stopping covariance")
    return s1 / count, covariance


def _positive_semidefinite(matrix: Any, numpy: Any, name: str) -> Any:
    symmetric = 0.5 * (matrix + matrix.T)
    eigenvalues, eigenvectors = numpy.linalg.eigh(symmetric)
    scale = max(1.0e-300, float(numpy.max(numpy.abs(eigenvalues))))
    if float(numpy.min(eigenvalues)) < -_PSD_RELATIVE_TOLERANCE * scale:
        raise ValueError(f"{name} is not positive semidefinite")
    clipped = numpy.maximum(eigenvalues, 0.0)
    return (eigenvectors * clipped) @ eigenvectors.T


def _mahalanobis(
    difference: Any,
    covariance: Any,
    numpy: Any,
    *,
    null_relative_tolerance: float | None = None,
) -> tuple[float, int, list[float]]:
    covariance = _positive_semidefinite(covariance, numpy, "difference covariance")
    eigenvalues, eigenvectors = numpy.linalg.eigh(covariance)
    maximum = max(0.0, float(numpy.max(eigenvalues)))
    # Do not reuse the deliberately loose PSD-negativity tolerance as a rank
    # cutoff.  The stopping vector can have well-resolved, low-variance modes
    # 1e-12--1e-14 below its leading eigenvalue.  Discarding those modes turns
    # ordinary sub-sigma differences into false "outside covariance support"
    # failures.  The conventional matrix-rank scale retains every eigenmode
    # distinguishable in binary64 while still excluding numerical zero modes.
    rank_relative_tolerance = (
        float(numpy.finfo("f8").eps) * max(covariance.shape)
    )
    threshold = max(
        float(numpy.finfo("f8").tiny), rank_relative_tolerance * maximum
    )
    active = eigenvalues > threshold
    rank = int(numpy.count_nonzero(active))
    difference_norm = float(numpy.linalg.norm(difference))
    if rank:
        coordinates = eigenvectors[:, active].T @ difference
        statistic = float(numpy.sum(coordinates * coordinates / eigenvalues[active]))
        projected = eigenvectors[:, active] @ coordinates
        null_residual = float(numpy.linalg.norm(difference - projected))
    else:
        statistic = 0.0
        null_residual = float(numpy.linalg.norm(difference))
    # Exact closure relations can leave roundoff-sized components in algebraic
    # null spaces after projection.  Replica callers supply the float32 storage
    # tolerance explicitly; exact stopping moments retain the tighter float64
    # / covariance-rank tolerance.
    relative_tolerance = max(
        128.0 * float(numpy.finfo("f8").eps),
        _PSD_RELATIVE_TOLERANCE,
        0.0 if null_relative_tolerance is None else null_relative_tolerance,
    )
    numerical_tolerance = max(
        relative_tolerance * max(
            difference_norm,
            math.sqrt(maximum),
            float(numpy.finfo("f8").tiny),
        ),
        # Directions discarded at the binary64 rank threshold are not known
        # to have exactly zero variance.  Permit the same familywise 5-sigma
        # excursion used by the calibration gate at the largest unresolved
        # standard deviation.  Truly deterministic mismatches larger than the
        # numerical covariance resolution still fail below.
        MAX_ABS_STANDARDIZED_RESIDUAL * math.sqrt(threshold),
        float(numpy.finfo("f8").tiny),
    )
    if null_residual > numerical_tolerance:
        raise ValueError("observed difference lies outside covariance support")

    diagonal = numpy.diag(covariance)
    standardized: list[float] = []
    for index, variance in enumerate(diagonal):
        if variance > threshold:
            standardized.append(float(difference[index] / math.sqrt(float(variance))))
        elif abs(float(difference[index])) > numerical_tolerance:
            raise ValueError("deterministic observable differs between null blocks")
    return statistic, rank, standardized


def _coarse_diagnostics(
    root: Any, pairs: Sequence[tuple[str, str, str]], numpy: Any
) -> dict[str, Any]:
    residuals: list[float] = []
    comparisons: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    statistic = 0.0
    degrees = 0
    dimension_records = [
        {
            "index": index,
            "name": (
                f"registered_energy_band_{index:02d}"
                if index < 12
                else f"legendre_l{index - 11}"
            ),
            "stochastic_comparisons": 0,
            "deterministic_comparisons": 0,
            "nonzero_mean_comparisons": 0,
        }
        for index in range(16)
    ]
    for pair in pairs:
        label = _comparison_label(pair)
        try:
            left_mean, left_covariance = _stopping_moments(
                root[f"blocks/{pair[0]}/{pair[1]}"], numpy
            )
            right_mean, right_covariance = _stopping_moments(
                root[f"blocks/{pair[0]}/{pair[2]}"], numpy
            )
            difference = left_mean - right_mean
            difference_covariance = left_covariance + right_covariance
            value, rank, z_values = _mahalanobis(
                difference, difference_covariance, numpy
            )
            if rank <= 0:
                raise ValueError("stopping comparison has no stochastic covariance rank")
            probability = _two_sided_probability(_chi_square_cdf(value, rank))
            statistic += value
            degrees += rank
            residuals.extend(z_values)
            diagonal = numpy.diag(difference_covariance)
            diagonal_threshold = max(
                float(numpy.finfo("f8").tiny),
                _PSD_RELATIVE_TOLERANCE
                * max(0.0, float(numpy.max(numpy.abs(diagonal)))),
            )
            mean_scale = numpy.maximum.reduce(
                (
                    numpy.ones(16, dtype="f8"),
                    numpy.abs(left_mean),
                    numpy.abs(right_mean),
                )
            )
            for index, record in enumerate(dimension_records):
                stochastic = float(diagonal[index]) > diagonal_threshold
                record["stochastic_comparisons"] += int(stochastic)
                record["deterministic_comparisons"] += int(not stochastic)
                record["nonzero_mean_comparisons"] += int(
                    max(abs(float(left_mean[index])), abs(float(right_mean[index])))
                    > 256.0 * float(numpy.finfo("f8").eps) * float(mean_scale[index])
                )
            comparisons.append(
                {
                    **label,
                    "mahalanobis_squared": value,
                    "covariance_rank": rank,
                    "two_sided_probability": probability,
                    "active_marginal_count": len(z_values),
                }
            )
        except (KeyError, ValueError, FloatingPointError) as error:
            failures.append({**label, "reason": str(error)})

    omnibus = (
        _two_sided_probability(_chi_square_cdf(statistic, degrees)) if degrees else 0.0
    )
    coverage = _coverage(residuals)
    for record in dimension_records:
        record["required"] = bool(record["nonzero_mean_comparisons"])
        record["passed"] = bool(
            not record["required"]
            or int(record["stochastic_comparisons"]) >= MINIMUM_COMPARISONS
        )
    passed = bool(
        len(pairs) >= MINIMUM_COMPARISONS
        and len(comparisons) == len(pairs)
        and not failures
        and degrees > 0
        and omnibus >= OMNIBUS_TWO_SIDED_ALPHA
        and coverage["max_abs_standardized_residual"]
        <= MAX_ABS_STANDARDIZED_RESIDUAL
        and all(bool(record["passed"]) for record in dimension_records)
    )
    return {
        "passed": passed,
        "comparison_count": len(comparisons),
        "chi_square": statistic,
        "degrees_of_freedom": degrees,
        "omnibus_two_sided_probability": omnibus,
        "coverage": coverage,
        "dimensions": dimension_records,
        "comparisons": comparisons,
        "failures": failures,
    }


def _total_continuum_projection_matrix(
    root: Any, cell_count: int, numpy: Any
) -> tuple[Any, dict[str, Any]]:
    try:
        mu_edges = numpy.asarray(root["axes/mu_edges"][...], dtype="f8")
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("archive lacks mu edges for replica projections") from error
    if (
        mu_edges.ndim != 1
        or len(mu_edges) < 2
        or not numpy.all(numpy.isfinite(mu_edges))
        or numpy.any(numpy.diff(mu_edges) <= 0.0)
        or not math.isclose(float(mu_edges[0]), -1.0, abs_tol=1.0e-12)
        or not math.isclose(float(mu_edges[-1]), 1.0, abs_tol=1.0e-12)
    ):
        raise ValueError("archive mu edges are invalid")
    n_mu = len(mu_edges) - 1
    if cell_count <= 0 or cell_count % n_mu:
        raise ValueError("replica total-continuum cells do not match the mu axis")
    n_energy = cell_count // n_mu
    if n_energy < 12:
        raise ValueError("replica calibration requires at least 12 energy bins")

    projection = numpy.zeros((16, n_energy, n_mu), dtype="f8")
    bands: list[dict[str, Any]] = []
    for band in range(12):
        begin = band * n_energy // 12
        end = (band + 1) * n_energy // 12
        projection[band, begin:end, :] = 1.0
        bands.append(
            {
                "name": f"registered_energy_band_{band:02d}",
                "energy_index_begin": begin,
                "energy_index_end": end,
                "mu_index_begin": 0,
                "mu_index_end": n_mu,
            }
        )

    energy = (numpy.arange(n_energy, dtype="f8") + 0.5) / n_energy
    mu = 0.5 * (mu_edges[:-1] + mu_edges[1:])
    projection[12] = (0.2 + 0.8 * energy[:, None]) * (0.75 + 0.25 * mu[None, :])
    projection[13] = numpy.exp(-1.7 * energy[:, None]) * (
        0.75 - 0.25 * mu[None, :]
    )
    projection[14] = numpy.sqrt(energy[:, None]) * (
        0.5 + 0.5 * mu[None, :] * mu[None, :]
    )
    projection[15] = numpy.exp(-0.5 * ((energy[:, None] - 0.55) / 0.18) ** 2) * (
        0.8 + 0.2 * mu[None, :]
    )
    definitions: dict[str, Any] = {
        "schema": PROJECTION_SCHEMA,
        "tensor": "continuum_total",
        "flattening": "C-order [registered_energy_index,mu_index]",
        "energy_bins": n_energy,
        "mu_edges": [float(value) for value in mu_edges],
        "coarse_rebinnings": bands,
        "smooth_nonnegative_folds": [
            "(0.2+0.8*x)*(0.75+0.25*mu)",
            "exp(-1.7*x)*(0.75-0.25*mu)",
            "sqrt(x)*(0.5+0.5*mu^2)",
            "exp(-0.5*((x-0.55)/0.18)^2)*(0.8+0.2*mu)",
        ],
    }
    definitions["projection_sha256"] = hashlib.sha256(_canonical(definitions)).hexdigest()
    return projection.reshape(16, -1), definitions


def _replica_projection_matrix(
    root: Any, tensor_name: str, cell_count: int, numpy: Any
) -> tuple[Any, dict[str, Any]]:
    """Build a deterministic, low-dimensional projection for each map family."""

    if tensor_name == "continuum_total":
        return _total_continuum_projection_matrix(root, cell_count, numpy)

    components = _json_attribute(
        root["components"], "descriptors_json", "component descriptors"
    )
    lines = _json_attribute(
        root["axes"], "line_definitions_json", "line definitions"
    )
    mu_edges = numpy.asarray(root["axes/mu_edges"][...], dtype="f8")
    n_component = len(components)
    n_line = len(lines)
    n_mu = len(mu_edges) - 1
    if n_component <= 0 or n_line <= 0 or n_mu <= 0:
        raise ValueError("archive taxonomy is empty for replica projections")
    mu = 0.5 * (mu_edges[:-1] + mu_edges[1:])
    rows: list[Any] = []
    observables: list[dict[str, Any]] = []

    def append(name: str, values: Any, kind: str) -> None:
        flattened = numpy.asarray(values, dtype="f8").reshape(-1)
        if flattened.shape != (cell_count,):
            raise ValueError(f"replica projection {name!r} has the wrong shape")
        rows.append(flattened)
        observables.append({"name": name, "kind": kind})

    if tensor_name == "continuum_components":
        if cell_count % (n_component * n_mu):
            raise ValueError("component replica maps do not match component/mu axes")
        n_energy = cell_count // (n_component * n_mu)
        shape = (n_component, n_energy, n_mu)
        for component_index, component in enumerate(components):
            values = numpy.zeros(shape, dtype="f8")
            values[component_index, :, :] = 1.0
            append(str(component["name"]), values, "component_yield")
        for origin in sorted({int(component["origin"]) for component in components}):
            values = numpy.zeros(shape, dtype="f8")
            for component_index, component in enumerate(components):
                if int(component["origin"]) == origin:
                    values[component_index, : max(1, n_energy // 2), :] = 1.0
            append(f"origin_{origin}_low_registered_energy", values, "coarse_rebin")
        energy = (numpy.arange(n_energy, dtype="f8") + 0.5) / n_energy
        fold_definitions = (
            ("forward_hemisphere", numpy.broadcast_to((mu >= 0.0), (n_energy, n_mu))),
            ("upper_registered_half", numpy.broadcast_to((energy[:, None] >= 0.5), (n_energy, n_mu))),
            (
                "smooth_fold_0",
                (0.2 + 0.8 * energy[:, None]) * (0.75 + 0.25 * mu[None, :]),
            ),
            (
                "smooth_fold_1",
                numpy.exp(-1.7 * energy[:, None]) * (0.75 - 0.25 * mu[None, :]),
            ),
        )
        for name, fold in fold_definitions:
            append(
                name,
                numpy.broadcast_to(fold, shape),
                "coarse_rebin" if "half" in name or "hemisphere" in name else "smooth_fold",
            )
    elif tensor_name in {"line_components", "line_total"}:
        if tensor_name == "line_components":
            if cell_count != n_component * n_line * n_mu:
                raise ValueError("line-component replica maps have the wrong shape")
            shape = (n_component, n_line, n_mu)
            for component_index, component in enumerate(components):
                for line_index, line in enumerate(lines):
                    values = numpy.zeros(shape, dtype="f8")
                    values[component_index, line_index, :] = 1.0
                    append(
                        f"{component['name']}:{line['name']}",
                        values,
                        "component_line_yield",
                    )
            for line_index, line in enumerate(lines):
                for order in range(1, 5):
                    values = numpy.zeros(shape, dtype="f8")
                    values[:, line_index, :] = numpy.polynomial.legendre.legval(
                        mu, [0.0] * order + [1.0]
                    )
                    append(
                        f"{line['name']}:legendre_l{order}",
                        values,
                        "line_angular_moment",
                    )
        else:
            if cell_count != n_line * n_mu:
                raise ValueError("total-line replica maps have the wrong shape")
            shape = (n_line, n_mu)
            for line_index, line in enumerate(lines):
                values = numpy.zeros(shape, dtype="f8")
                values[line_index, :] = 1.0
                append(str(line["name"]), values, "line_yield")
            for line_index, line in enumerate(lines):
                for order in range(1, 5):
                    values = numpy.zeros(shape, dtype="f8")
                    values[line_index, :] = numpy.polynomial.legendre.legval(
                        mu, [0.0] * order + [1.0]
                    )
                    append(
                        f"{line['name']}:legendre_l{order}",
                        values,
                        "line_angular_moment",
                    )
    else:  # pragma: no cover - fixed caller inventory
        raise ValueError(f"unknown replica tensor {tensor_name!r}")

    if not rows:
        raise ValueError(f"{tensor_name} has no replica calibration projections")
    definitions: dict[str, Any] = {
        "schema": PROJECTION_SCHEMA,
        "tensor": tensor_name,
        "flattening": (
            "C-order [component,registered_energy,mu]"
            if tensor_name == "continuum_components"
            else "C-order [component,line,mu]"
            if tensor_name == "line_components"
            else "C-order [line,mu]"
        ),
        "input_cells": cell_count,
        "projection_count": len(rows),
        "observables": observables,
    }
    definitions["projection_sha256"] = hashlib.sha256(_canonical(definitions)).hexdigest()
    return numpy.stack(rows, axis=0), definitions


def _replica_observations(
    group: Any,
    root: Any,
    tensor_name: str,
    numpy: Any,
    projection_cache: dict[tuple[str, int], tuple[Any, dict[str, Any]]] | None = None,
) -> tuple[Any, dict[str, Any]]:
    try:
        counts = numpy.asarray(group["history_count"][...], dtype="f8")
        s1 = numpy.asarray(group[f"replica_s1/{tensor_name}"][...], dtype="f8")
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"block lacks {tensor_name} replica maps/counts") from error
    if counts.shape != (32,) or s1.ndim != 2 or s1.shape[0] != 32:
        raise ValueError("replica calibration requires exactly 32 maps and counts")
    if (
        not numpy.all(numpy.isfinite(s1))
        or numpy.any(s1 < 0.0)
        or not numpy.all(numpy.isfinite(counts))
        or numpy.any(counts <= 0.0)
        or not numpy.all(counts == counts[0])
    ):
        raise ValueError("replica maps/counts are invalid or unequal within a block")
    cache_key = (tensor_name, int(s1.shape[1]))
    cached = projection_cache.get(cache_key) if projection_cache is not None else None
    if cached is None:
        cached = _replica_projection_matrix(
            root, tensor_name, s1.shape[1], numpy
        )
        if projection_cache is not None:
            projection_cache[cache_key] = cached
    projection, definitions = cached
    per_history = s1 / counts[:, None]
    return per_history @ projection.T, definitions


def _hotelling(
    left: Any, right: Any, numpy: Any
) -> tuple[float, int, int, float, list[float]]:
    n_left, dimension = left.shape
    n_right, right_dimension = right.shape
    if dimension != right_dimension or n_left < 2 or n_right < 2:
        raise ValueError("replica observations have incompatible dimensions")
    left_covariance = numpy.cov(left, rowvar=False, ddof=1)
    right_covariance = numpy.cov(right, rowvar=False, ddof=1)
    pooled = (
        (n_left - 1) * left_covariance + (n_right - 1) * right_covariance
    ) / (n_left + n_right - 2)
    difference = numpy.mean(left, axis=0) - numpy.mean(right, axis=0)
    statistic, rank, _ = _mahalanobis(
        difference,
        pooled,
        numpy,
        null_relative_tolerance=_FLOAT32_REPLICA_NULL_TOLERANCE,
    )
    if rank <= 0:
        raise ValueError("replica comparison has no stochastic covariance rank")
    denominator_df = n_left + n_right - rank - 1
    if denominator_df <= 0:
        raise ValueError("replica projection rank is too large for Hotelling calibration")
    # _mahalanobis used pooled covariance, so scale to the conventional T^2.
    t_squared = n_left * n_right / (n_left + n_right) * statistic
    f_statistic = denominator_df * t_squared / (
        rank * (n_left + n_right - 2)
    )
    probability = _two_sided_probability(_f_cdf(f_statistic, rank, denominator_df))

    mean_covariance = left_covariance / n_left + right_covariance / n_right
    diagonal = numpy.diag(mean_covariance)
    threshold = max(
        float(numpy.finfo("f8").tiny),
        _PSD_RELATIVE_TOLERANCE * max(0.0, float(numpy.max(diagonal))),
    )
    standardized = [
        float(difference[index] / math.sqrt(float(variance)))
        for index, variance in enumerate(diagonal)
        if variance > threshold
    ]
    return f_statistic, rank, denominator_df, probability, standardized


def _replica_diagnostics(
    root: Any, pairs: Sequence[tuple[str, str, str]], numpy: Any
) -> dict[str, Any]:
    family_reports: dict[str, Any] = {}
    projection_families: dict[str, Any] = {}
    all_residuals: list[float] = []
    all_failures: list[dict[str, Any]] = []
    projection_cache: dict[
        tuple[str, int], tuple[Any, dict[str, Any]]
    ] = {}
    for tensor_name in _REQUIRED_MOMENT_TENSORS:
        residuals: list[float] = []
        comparisons: list[dict[str, Any]] = []
        failures: list[dict[str, Any]] = []
        probabilities: list[float] = []
        projection_record: dict[str, Any] | None = None
        for pair in pairs:
            label = _comparison_label(pair)
            try:
                left_group = root[f"blocks/{pair[0]}/{pair[1]}"]
                right_group = root[f"blocks/{pair[0]}/{pair[2]}"]
                left, left_projection = _replica_observations(
                    left_group, root, tensor_name, numpy, projection_cache
                )
                right, right_projection = _replica_observations(
                    right_group, root, tensor_name, numpy, projection_cache
                )
                if left_projection != right_projection:
                    raise ValueError(
                        f"{tensor_name} replica projection changed between blocks"
                    )
                if projection_record is None:
                    projection_record = left_projection
                elif projection_record != left_projection:
                    raise ValueError(
                        f"{tensor_name} replica projection changed between nodes"
                    )
                left_counts = numpy.asarray(
                    left_group["history_count"][...], dtype="u8"
                )
                right_counts = numpy.asarray(
                    right_group["history_count"][...], dtype="u8"
                )
                if not numpy.array_equal(left_counts, right_counts):
                    raise ValueError(
                        "paired blocks need equal replica histories for the Hotelling null test"
                    )
                f_value, rank, denominator, probability, z_values = _hotelling(
                    left, right, numpy
                )
                probabilities.append(probability)
                residuals.extend(z_values)
                comparisons.append(
                    {
                        **label,
                        "hotelling_f": f_value,
                        "numerator_degrees_of_freedom": rank,
                        "denominator_degrees_of_freedom": denominator,
                        "two_sided_probability": probability,
                        "active_marginal_count": len(z_values),
                    }
                )
            except (KeyError, TypeError, ValueError, FloatingPointError) as error:
                failures.append({**label, "reason": str(error)})

        fisher = -2.0 * math.fsum(
            math.log(max(value, 1.0e-300)) for value in probabilities
        )
        fisher_degrees = 2 * len(probabilities)
        fisher_upper_probability = (
            1.0 - _chi_square_cdf(fisher, fisher_degrees)
            if fisher_degrees
            else 0.0
        )
        coverage = _coverage(residuals)
        passed = bool(
            len(pairs) >= MINIMUM_COMPARISONS
            and len(comparisons) == len(pairs)
            and not failures
            and projection_record is not None
            and min(probabilities, default=0.0)
            >= INDIVIDUAL_REPLICA_TWO_SIDED_ALPHA
            and fisher_upper_probability >= OMNIBUS_TWO_SIDED_ALPHA
            and coverage["count"] >= MINIMUM_COMPARISONS
            and coverage["max_abs_standardized_residual"]
            <= MAX_ABS_STANDARDIZED_RESIDUAL
        )
        family_reports[tensor_name] = {
            "passed": passed,
            "comparison_count": len(comparisons),
            "fisher_combined_statistic": fisher,
            "fisher_degrees_of_freedom": fisher_degrees,
            "fisher_upper_tail_probability": fisher_upper_probability,
            "minimum_individual_two_sided_probability": min(
                probabilities, default=0.0
            ),
            "coverage": coverage,
            "projection_definition": projection_record,
            "comparisons": comparisons,
            "failures": failures,
        }
        if projection_record is not None:
            projection_families[tensor_name] = projection_record
        all_residuals.extend(residuals)
        all_failures.extend(
            {"tensor": tensor_name, **failure} for failure in failures
        )

    projection_definition: dict[str, Any] = {
        "schema": PROJECTION_SCHEMA,
        "families": projection_families,
    }
    projection_definition["projection_sha256"] = hashlib.sha256(
        _canonical(projection_definition)
    ).hexdigest()
    inventory_complete = set(projection_families) == set(_REQUIRED_MOMENT_TENSORS)
    return {
        "passed": bool(
            inventory_complete
            and all(family_reports[name]["passed"] for name in _REQUIRED_MOMENT_TENSORS)
        ),
        "comparison_count": len(pairs),
        "family_test_count": sum(
            int(report["comparison_count"]) for report in family_reports.values()
        ),
        "family_inventory_complete": inventory_complete,
        "coverage": _coverage(all_residuals),
        "projection_definition": projection_definition,
        "families": family_reports,
        "failures": all_failures,
    }


def analyze_uncertainty_calibration(archive_path: str | Path) -> dict[str, Any]:
    """Analyze one immutable-block archive and return a publication-ready report."""

    h5py, numpy = require_hdf5()
    archive = Path(archive_path)
    with h5py.File(archive, "r") as root:
        if str(root.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
            raise ValueError("unsupported CEEFIT fit-archive schema")
        try:
            config_hash = str(root["manifest"].attrs["config_sha256"])
        except (KeyError, TypeError) as error:
            raise ValueError("fit archive lacks config_sha256") from error
        if not _is_sha256(config_hash):
            raise ValueError("fit archive config_sha256 is invalid")
        snapshot_hash = archive_snapshot_sha256(root)
        pairs, unused = _pair_inventory(root)
        independence = _independence_diagnostics(root, pairs, numpy)
        s2 = _s2_diagnostics(root, pairs)
        cell_s2 = _cell_s2_diagnostics(root, pairs, numpy)
        coarse = _coarse_diagnostics(root, pairs, numpy)
        replica = _replica_diagnostics(root, pairs, numpy)

    comparison_inventory = [_comparison_label(pair) for pair in pairs]
    sufficient_replication = len(pairs) >= MINIMUM_COMPARISONS
    independent_replication = bool(independence["passed"])
    if not independent_replication:
        s2["passed"] = False
        cell_s2["passed"] = False
        coarse["passed"] = False
        replica["passed"] = False
    s2_passed = bool(s2["passed"] and cell_s2["passed"])
    family_coverage_passed = bool(
        cell_s2["family_inventory_complete"]
        and all(
            bool(segment["passed"])
            for tensor in cell_s2["tensors"].values()
            for segment in tensor["segments"]
        )
        and replica["family_inventory_complete"]
        and all(
            bool(dimension["passed"])
            for dimension in coarse["dimensions"]
        )
    )
    failures: list[dict[str, Any]] = []
    if not sufficient_replication:
        failures.append(
            {
                "gate": "independent_replication",
                "reason": (
                    f"calibration has {len(pairs)} independent comparisons; "
                    f"at least {MINIMUM_COMPARISONS} are required"
                ),
            }
        )
    if not independent_replication:
        failures.append(
            {
                "gate": "independent_replication",
                "reason": "block/replica identity checks did not establish independent null runs",
            }
        )
    if not cell_s2["family_inventory_complete"]:
        failures.append(
            {
                "gate": "family_coverage",
                "reason": "cell-S2 calibration lacks one or more required tensor families",
            }
        )
    for tensor_name, tensor in cell_s2["tensors"].items():
        for segment in tensor["segments"]:
            if not bool(segment["required"]) or bool(segment["passed"]):
                continue
            name = str(segment["name"])
            reason = (
                "required uncertainty family has no stochastic support"
                if int(segment["stochastic_comparisons"]) == 0
                else "required uncertainty family has too few stochastic comparisons"
            )
            if "annihilation" in name:
                reason += (
                    f"; provide {MINIMUM_COMPARISONS} independent block-pair "
                    "comparisons at or above the 1021.9979-keV pair threshold "
                    f"({2 * MINIMUM_COMPARISONS} blocks total)"
                )
            failures.append(
                {
                    "gate": "family_coverage",
                    "tensor": tensor_name,
                    "family": name,
                    "stochastic_comparisons": int(
                        segment["stochastic_comparisons"]
                    ),
                    "required_comparisons": MINIMUM_COMPARISONS,
                    "reason": reason,
                }
            )
    if not replica["family_inventory_complete"]:
        failures.append(
            {
                "gate": "family_coverage",
                "reason": "replica calibration lacks one or more required tensor families",
            }
        )
    for section_name, section in (
        ("raw_total_s2", s2),
        ("cell_s2", cell_s2),
        ("stopping_outer_moments", coarse),
        ("replica_covariance", replica),
    ):
        failures.extend(
            {"gate": section_name, **failure}
            for failure in section.get("failures", ())
        )
    for gate, passed, reason in (
        ("s2", s2_passed, "raw or per-cell S2 calibration did not pass"),
        (
            "stopping_outer_moments",
            bool(coarse["passed"]),
            "stopping-vector outer-moment calibration did not pass",
        ),
        (
            "replica_covariance",
            bool(replica["passed"]),
            "32-replica covariance calibration did not pass",
        ),
        (
            "family_coverage",
            family_coverage_passed,
            "one or more required uncertainty families lack calibrated support",
        ),
    ):
        if not passed:
            failures.append({"gate": gate, "reason": reason})
    passed = bool(
        sufficient_replication
        and independent_replication
        and family_coverage_passed
        and s2_passed
        and coarse["passed"]
        and replica["passed"]
        and not failures
    )
    report: dict[str, Any] = {
        "schema": CALIBRATION_SCHEMA,
        "config_sha256": config_hash,
        "archive_snapshot_sha256": snapshot_hash,
        "passed": passed,
        "repeated_null_runs": len(pairs),
        "s2_passed": s2_passed,
        "coarse_covariance_passed": bool(coarse["passed"]),
        "replica_covariance_passed": bool(replica["passed"]),
        "family_coverage_passed": family_coverage_passed,
        "failures": failures,
        "criteria": {
            "minimum_independent_comparisons": MINIMUM_COMPARISONS,
            "pairing": (
                "numeric block IDs paired adjacently within each node; every block "
                "is used in at most one comparison"
            ),
            "required_moment_and_replica_tensors": list(
                _REQUIRED_MOMENT_TENSORS
            ),
            "omnibus_two_sided_alpha": OMNIBUS_TWO_SIDED_ALPHA,
            "minimum_individual_replica_two_sided_probability": (
                INDIVIDUAL_REPLICA_TWO_SIDED_ALPHA
            ),
            "maximum_absolute_standardized_residual": (
                MAX_ABS_STANDARDIZED_RESIDUAL
            ),
            "minimum_cell_fraction_within_1_96sigma": (
                _MINIMUM_WITHIN_1_96_FRACTION
            ),
            "minimum_cell_fraction_within_3sigma": _MINIMUM_WITHIN_3_FRACTION,
            "cell_maximum_residual_gate": (
                "max(5,sqrt(2*log(2*n/omnibus_alpha))) per tensor"
            ),
        },
        "replication": {
            "sufficient": sufficient_replication,
            "independent": independent_replication,
            "comparison_count": len(pairs),
            "blocks_used": 2 * len(pairs),
            "comparison_inventory": comparison_inventory,
            "unpaired_blocks": unused,
            "independence_validation": independence,
        },
        "raw_total_s2": s2,
        "cell_s2": cell_s2,
        "stopping_outer_moments": coarse,
        "replica_covariance": replica,
    }
    # Exercise strict serialization before returning even when no output is requested.
    _canonical(report)
    return report


def calibrate_uncertainty_archive(
    archive_path: str | Path, output_path: str | Path
) -> dict[str, Any]:
    """Analyze ``archive_path`` and atomically write strict calibration JSON."""

    archive = Path(archive_path)
    output = Path(output_path)
    if archive.resolve() == output.resolve():
        raise ValueError(
            "uncertainty-calibration output must be distinct from the HDF5 archive"
        )
    report = analyze_uncertainty_calibration(archive_path)
    _atomic_json(output_path, report)
    return report
