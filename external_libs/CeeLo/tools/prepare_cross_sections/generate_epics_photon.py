#!/usr/bin/env python3
"""Generate compact process-specific photon curves from locked EPICS2023 EPDL.

Endpoints, thresholds, and both sides of duplicated absorption edges are
mandatory.  All other native source nodes are validation samples rather than
automatically retained.  The adaptive reducer checks every source node plus 16
logarithmic probes in each retained interval and inserts the worst failing
sample until the accuracy gates pass.  Values are then stored with a per-curve
uint16 affine encoding; the committed errors include that decode step.
Rayleigh total cross sections jointly refine one grid for each four adjacent
elements, while preserving the same per-element error tests.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from epics_endf import EndfFile, Tab1
from fetch_sources import load_lock, sha256_file
from generation_utils import ELEMENT_SYMBOLS, as_float32, format_float32, format_rows


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DEFAULT_LOCK = HERE / "sources.lock.json"
DEFAULT_OUTPUT = ROOT / "src" / "cross_sections" / "photon_epics_data.cpp"
DEFAULT_REPORT_DIR = HERE / "reports"

ENERGY_MIN_EV = 10_000.0
ENERGY_MAX_EV = 10_000_000.0
ANGULAR_MIN = 1.0e-2
ANGULAR_REFERENCE_MIN = 1.0e-7
ANGULAR_MAX = 1.0e9
PAIR_THRESHOLD_EV = 1_022_000.0
EDGE_RELATIVE_WIDTH = 2.0e-6
PROBES_PER_INTERVAL = 16
MAX_RELATIVE_ERROR = 1.0e-3
MAX_INTEGRATED_ERROR = 5.0e-4
MAX_ANGULAR_CDF_ERROR = 2.5e-3
MAX_ANGULAR_MEAN_COS_ERROR = 2.0e-3
MAX_RAYLEIGH_CDF_QUANTIZATION_ERROR = 8.0e-6
MAX_RAYLEIGH_CDF_DISCRETIZATION_ERROR = 2.5e-3
RAYLEIGH_CDF_N = 256
RAYLEIGH_XS_GROUP_SIZE = 4
LOG_FLOOR = -30.0
HC_KEV_ANGSTROM = 12.398
MEC2_KEV = 510.998950

PROCESS_SPECS = (
    ("rayleigh", 23, 502, "log_log"),
    ("compton", 23, 504, "log_log"),
    ("pair_production", 23, 516, "log_log"),
    ("photoelectric", 23, 522, "log_log"),
    ("k_photoelectric", 23, 534, "log_log"),
)
ANGULAR_SPECS = (
    ("form_factor", 27, 502, "log_log"),
    ("scattering_function", 27, 504, "log_linear"),
)
RUNTIME_ANGULAR_SPECS = (
    ("scattering_function", 27, 504, "log_linear"),
)


@dataclass
class EncodedCurve:
    argument: list[float]
    packed_value: list[int]
    value_offset: float
    value_scale: float
    native_count: int
    maximum_error: float
    p95_error: float
    integrated_error: float


@dataclass
class CompactStorage:
    rayleigh_log_energy: list[float]
    rayleigh_value: list[int]
    rayleigh_grid_offset: list[int]
    rayleigh_value_offset: list[int]
    energy_pool: list[float]
    process_grid_index: list[int]
    process_value: list[int]
    process_offset: dict[tuple[int, str], int]
    angular_log_x: list[float]
    angular_value: list[int]
    angular_offset: dict[tuple[int, str], int]
    rayleigh_sampling_x: list[float]
    rayleigh_cdf_q: list[list[int]]


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = int(math.floor(position)), int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def safe_log10(value: float) -> float:
    return max(math.log10(value), LOG_FLOOR) if value > 0.0 and math.isfinite(value) else LOG_FLOOR


def log_probe(a: float, b: float, index: int, count: int = PROBES_PER_INTERVAL) -> float:
    fraction = index / (count + 1)
    return math.exp(math.log(a) + fraction * math.log(b / a))


def source_points(
    table: Tab1, lower: float, upper: float
) -> tuple[list[tuple[float, float]], list[float], int]:
    """Return mandatory points and the complete native validation abscissae."""
    mandatory: dict[float, float] = {lower: table.evaluate(lower)}
    validation_arguments = [lower]
    native_count = 0
    index = 0
    while index < len(table.x):
        end = index + 1
        while end < len(table.x) and table.x[end] == table.x[index]:
            end += 1
        argument = table.x[index]
        if lower <= argument <= upper:
            native_count += end - index
            validation_arguments.append(argument)
            values = table.y[index:end]
            if len(values) > 1:
                flank = max(lower, argument * (1.0 - EDGE_RELATIVE_WIDTH))
                mandatory[flank] = values[0]
                mandatory[argument] = values[-1]
                validation_arguments.append(flank)
        index = end
    mandatory[upper] = table.evaluate(upper)
    validation_arguments.append(upper)
    return (
        sorted(mandatory.items()),
        sorted(set(validation_arguments)),
        native_count,
    )


def encoded_value_cached(
    points: list[tuple[float, float]],
    arguments: list[float],
    argument: float,
    mode: str,
) -> float:
    if argument <= points[0][0]:
        return points[0][1]
    if argument >= points[-1][0]:
        return points[-1][1]
    import bisect
    upper = bisect.bisect_right(arguments, argument)
    lower = upper - 1
    x0, y0 = points[lower]
    x1, y1 = points[upper]
    t = math.log(argument / x0) / math.log(x1 / x0)
    if mode == "log_log" and y0 > 0.0 and y1 > 0.0:
        return math.exp(math.log(y0) + t * math.log(y1 / y0))
    return y0 + t * (y1 - y0)


def encoded_value(points: list[tuple[float, float]], argument: float, mode: str) -> float:
    return encoded_value_cached(points, [point[0] for point in points], argument, mode)


def pack_values(stored: list[float], mode: str) -> tuple[list[int], float, float, list[float]]:
    """Affine uint16 encoding, reserving zero for log-domain zero/floor."""
    if mode == "log_log":
        finite = [value for value in stored if value > LOG_FLOOR]
        if not finite:
            return [0] * len(stored), 0.0, 0.0, [LOG_FLOOR] * len(stored)
        offset = as_float32(min(finite))
        span = max(finite) - offset
        scale = as_float32(span / 65534.0) if span > 0.0 else 0.0
        packed = []
        recovered = []
        for value in stored:
            if value <= LOG_FLOOR:
                packed.append(0)
                recovered.append(LOG_FLOOR)
                continue
            code = 1 if scale == 0.0 else 1 + round((value - offset) / scale)
            code = min(65535, max(1, code))
            packed.append(code)
            recovered.append(offset + (code - 1) * scale)
        return packed, offset, scale, recovered

    offset = as_float32(min(stored))
    span = max(stored) - offset
    scale = as_float32(span / 65535.0) if span > 0.0 else 0.0
    packed = [
        min(65535, max(0, 0 if scale == 0.0 else round((value - offset) / scale)))
        for value in stored
    ]
    recovered = [offset + code * scale for code in packed]
    return packed, offset, scale, recovered


def quantize_points(
    points: list[tuple[float, float]], mode: str
) -> tuple[list[tuple[float, float]], list[int], float, float]:
    coordinates: list[float] = []
    stored: list[float] = []
    for argument, value in points:
        log_argument = as_float32(math.log10(argument))
        if mode in ("log_log", "log_linear"):
            encoded_value = as_float32(safe_log10(value))
        else:
            encoded_value = as_float32(value)
        recovered_argument = 10.0**log_argument
        if coordinates and recovered_argument <= coordinates[-1]:
            # Some EPICS source nodes are closer than a float32 log argument can
            # distinguish. Keep the later node deterministically. Deliberate edge
            # flanks are wider than this and therefore remain distinct.
            stored[-1] = encoded_value
            continue
        coordinates.append(recovered_argument)
        stored.append(encoded_value)
    packed, offset, scale, recovered_stored = pack_values(
        stored, "log_log" if mode in ("log_log", "log_linear") else mode
    )
    recovered = [
        10.0**value if mode in ("log_log", "log_linear") and value > LOG_FLOOR else
        (0.0 if mode in ("log_log", "log_linear") else value)
        for value in recovered_stored
    ]
    return list(zip(coordinates, recovered)), packed, offset, scale


def relative_error(reference: float, represented: float, material_floor: float) -> float:
    if reference <= material_floor:
        return 0.0
    return abs(represented / reference - 1.0)


def curve_errors(
    table: Tab1,
    points: list[tuple[float, float]],
    source_arguments: list[float],
    mode: str,
    material_floor: float,
    error_gate: float = MAX_RELATIVE_ERROR,
) -> tuple[list[float], tuple[float, float, int] | None, list[tuple[float, float]]]:
    import bisect

    encoded, _, _, _ = quantize_points(points, mode)
    encoded_arguments = [point[0] for point in encoded]
    errors: list[float] = []
    worst: tuple[float, float, int] | None = None
    refinements: list[tuple[float, float]] = []
    for interval in range(len(encoded) - 1):
        a, b = encoded[interval][0], encoded[interval + 1][0]
        # A duplicated source edge is represented by this deliberately narrow
        # transition. It is tested as an edge neighborhood, not as a smooth span.
        if b / a <= 1.0 + 2.0 * EDGE_RELATIVE_WIDTH:
            continue
        interval_worst: tuple[float, float] | None = None
        source_begin = bisect.bisect_left(source_arguments, a)
        source_end = bisect.bisect_right(source_arguments, b)
        candidates = source_arguments[source_begin:source_end]
        candidates.extend(
            log_probe(a, b, probe_index)
            for probe_index in range(1, PROBES_PER_INTERVAL + 1)
        )
        for argument in candidates:
            reference = table.evaluate(argument)
            represented = encoded_value_cached(encoded, encoded_arguments, argument, mode)
            error = relative_error(reference, represented, material_floor)
            errors.append(error)
            if worst is None or error > worst[0]:
                worst = (error, argument, interval)
            if interval_worst is None or error > interval_worst[0]:
                interval_worst = (error, argument)
        if interval_worst is not None and interval_worst[0] > error_gate:
            refinements.append(interval_worst)
    return errors, worst, refinements


def integrated_error(
    table: Tab1,
    points: list[tuple[float, float]],
    mode: str,
) -> tuple[float, int]:
    encoded, _, _, _ = quantize_points(points, mode)
    encoded_arguments = [point[0] for point in encoded]
    true_area = 0.0
    represented_area = 0.0
    interval_differences: list[float] = []
    for interval in range(len(encoded) - 1):
        a, b = encoded[interval][0], encoded[interval + 1][0]
        previous_x = a
        previous_true = table.evaluate(a)
        previous_represented = encoded_value_cached(encoded, encoded_arguments, a, mode)
        interval_true = 0.0
        interval_represented = 0.0
        for probe_index in range(1, PROBES_PER_INTERVAL + 2):
            x = b if probe_index == PROBES_PER_INTERVAL + 1 else log_probe(a, b, probe_index)
            true = table.evaluate(x)
            represented = encoded_value_cached(encoded, encoded_arguments, x, mode)
            width = x - previous_x
            true_piece = 0.5 * width * (previous_true + true)
            represented_piece = 0.5 * width * (previous_represented + represented)
            true_area += true_piece
            represented_area += represented_piece
            interval_true += true_piece
            interval_represented += represented_piece
            previous_x, previous_true, previous_represented = x, true, represented
        interval_differences.append(abs(interval_represented - interval_true))
    error = abs(represented_area / true_area - 1.0) if true_area > 0.0 else 0.0
    worst = max(range(len(interval_differences)), key=interval_differences.__getitem__)
    return error, worst


def build_curve(table: Tab1, lower: float, upper: float, mode: str) -> EncodedCurve:
    points, source_arguments, native_count = source_points(table, lower, upper)
    peak = max(table.evaluate(argument) for argument in source_arguments)
    # Pair production rises from an exact zero at threshold. Relative error is
    # undefined in that vanishing onset; test it once it reaches one ppm of the
    # curve peak, far below any material contribution. Other processes use a
    # substantially stricter ten-orders-below-peak material floor.
    material_fraction = 1.0e-6 if table.mt == 516 else 1.0e-10
    material_floor = max(peak * material_fraction, 1.0e-30)

    # The migration comparison has a tighter 0.05% p95 gate for pair
    # production. Its direct EPICS source agrees with the historical values to
    # printed precision, so use half the general representation tolerance and
    # leave margin for float32/runtime interpolation.
    relative_gate = 5.0e-4 if table.mt == 516 else MAX_RELATIVE_ERROR
    for _ in range(64):
        errors, worst, refinements = curve_errors(
            table, points, source_arguments, mode, material_floor, relative_gate
        )
        if worst is None or worst[0] <= relative_gate:
            break
        encoded_logs = {as_float32(math.log10(argument)) for argument, _ in points}
        added = []
        for _, argument in refinements:
            encoded_log = as_float32(math.log10(argument))
            if encoded_log not in encoded_logs:
                added.append(argument)
                encoded_logs.add(encoded_log)
        if not added:
            raise ValueError(
                f"MF={table.mf} MT={table.mt}: uint16/float32 encoding floor "
                f"prevents the {relative_gate:.4%} pointwise gate; worst={worst[0]:.4%}"
            )
        points.extend((argument, table.evaluate(argument)) for argument in added)
        points.sort()
    else:
        raise ValueError(f"MF={table.mf} MT={table.mt}: interval refinement did not converge")

    area_error, worst_area_interval = integrated_error(table, points, mode)
    for _ in range(64):
        if area_error <= MAX_INTEGRATED_ERROR:
            break
        # Split the interval contributing the largest absolute quadrature
        # discrepancy; this also sharpens an edge flank if that dominates.
        encoded, _, _, _ = quantize_points(points, mode)
        argument = math.sqrt(
            encoded[worst_area_interval][0] * encoded[worst_area_interval + 1][0]
        )
        points.append((argument, table.evaluate(argument)))
        points.sort()
        area_error, worst_area_interval = integrated_error(table, points, mode)
    if area_error > MAX_INTEGRATED_ERROR:
        raise ValueError(
            f"MF={table.mf} MT={table.mt}: integrated error {area_error:.4%} exceeds 0.05%"
        )

    # Adding an integration midpoint slightly changes a curve's affine uint16
    # scale. Re-close both gates together so packing can never invalidate the
    # pointwise result established above.
    for _ in range(64):
        errors, worst, refinements = curve_errors(
            table, points, source_arguments, mode, material_floor, relative_gate
        )
        area_error, worst_area_interval = integrated_error(table, points, mode)
        if ((worst is None or worst[0] <= relative_gate)
                and area_error <= MAX_INTEGRATED_ERROR):
            break
        if worst is not None and worst[0] > relative_gate:
            encoded_logs = {as_float32(math.log10(argument)) for argument, _ in points}
            added = []
            for _, argument in refinements:
                encoded_log = as_float32(math.log10(argument))
                if encoded_log not in encoded_logs:
                    added.append(argument)
                    encoded_logs.add(encoded_log)
            if not added:
                raise ValueError(
                    f"MF={table.mf} MT={table.mt}: packed pointwise gate cannot close"
                )
            points.extend((argument, table.evaluate(argument)) for argument in added)
        else:
            encoded, _, _, _ = quantize_points(points, mode)
            argument = math.sqrt(
                encoded[worst_area_interval][0] * encoded[worst_area_interval + 1][0]
            )
            points.append((argument, table.evaluate(argument)))
        points.sort()
    else:
        raise ValueError(f"MF={table.mf} MT={table.mt}: combined packed gates did not converge")

    encoded, packed, value_offset, value_scale = quantize_points(points, mode)
    arguments = [as_float32(math.log10(point[0])) for point in encoded]
    return EncodedCurve(
        arguments,
        packed,
        value_offset,
        value_scale,
        native_count,
        max(errors, default=0.0),
        percentile(errors, 0.95) if errors else 0.0,
        area_error,
    )


def build_shared_rayleigh_group(tables: list[Tab1]) -> list[EncodedCurve]:
    """Build one adaptive energy grid shared by adjacent-Z Rayleigh curves."""
    source_arguments: list[list[float]] = []
    native_counts: list[int] = []
    material_floors: list[float] = []
    for table in tables:
        _, arguments, native_count = source_points(
            table, ENERGY_MIN_EV, ENERGY_MAX_EV
        )
        source_arguments.append(arguments)
        native_counts.append(native_count)
        peak = max(table.evaluate(argument) for argument in arguments)
        material_floors.append(max(peak * 1.0e-10, 1.0e-30))

    grid = [ENERGY_MIN_EV, ENERGY_MAX_EV]
    final_errors: list[list[float]] = []
    final_areas: list[float] = []
    for _ in range(64):
        refinements: dict[int, tuple[float, float]] = {}
        final_errors = []
        final_areas = []
        all_pass = True
        for table, arguments, material_floor in zip(
            tables, source_arguments, material_floors
        ):
            points = [(argument, table.evaluate(argument)) for argument in grid]
            errors, worst, requested = curve_errors(
                table, points, arguments, "log_log", material_floor
            )
            area_error, worst_area_interval = integrated_error(
                table, points, "log_log"
            )
            final_errors.append(errors)
            final_areas.append(area_error)
            if worst is not None and worst[0] > MAX_RELATIVE_ERROR:
                all_pass = False
                for error, argument in requested:
                    interval = bisect.bisect_right(grid, argument) - 1
                    previous = refinements.get(interval)
                    if previous is None or error > previous[0]:
                        refinements[interval] = (error, argument)
            if area_error > MAX_INTEGRATED_ERROR:
                all_pass = False
                encoded, _, _, _ = quantize_points(points, "log_log")
                argument = math.sqrt(
                    encoded[worst_area_interval][0]
                    * encoded[worst_area_interval + 1][0]
                )
                interval = bisect.bisect_right(grid, argument) - 1
                refinements.setdefault(interval, (area_error, argument))
        if all_pass:
            break
        encoded_logs = {as_float32(math.log10(argument)) for argument in grid}
        added: list[float] = []
        for _, argument in refinements.values():
            encoded_log = as_float32(math.log10(argument))
            if encoded_log not in encoded_logs:
                added.append(argument)
                encoded_logs.add(encoded_log)
        if not added:
            raise ValueError(
                "shared Rayleigh grid cannot close the packed accuracy gates"
            )
        grid.extend(added)
        grid.sort()
    else:
        raise ValueError("shared Rayleigh grid refinement did not converge")

    curves: list[EncodedCurve] = []
    shared_arguments: list[float] | None = None
    for table, native_count, errors, area_error in zip(
        tables, native_counts, final_errors, final_areas
    ):
        encoded, packed, value_offset, value_scale = quantize_points(
            [(argument, table.evaluate(argument)) for argument in grid], "log_log"
        )
        arguments = [as_float32(math.log10(point[0])) for point in encoded]
        if shared_arguments is None:
            shared_arguments = arguments
        elif arguments != shared_arguments:
            raise ValueError("Rayleigh group did not retain an identical grid")
        curves.append(EncodedCurve(
            arguments, packed, value_offset, value_scale, native_count,
            max(errors, default=0.0),
            percentile(errors, 0.95) if errors else 0.0,
            area_error,
        ))
    return curves


def curve_runtime(curve: EncodedCurve, argument: float, mode: str) -> float:
    points = []
    for log_argument, code in zip(curve.argument, curve.packed_value):
        stored = (
            LOG_FLOOR if mode in ("log_log", "log_linear") and code == 0 else
            curve.value_offset + (code - 1) * curve.value_scale
            if mode in ("log_log", "log_linear") else
            curve.value_offset + code * curve.value_scale
        )
        value = 10.0**stored if mode in ("log_log", "log_linear") and stored > LOG_FLOOR else stored
        if mode in ("log_log", "log_linear") and stored <= LOG_FLOOR:
            value = 0.0
        points.append((10.0**log_argument, value))
    return encoded_value(points, argument, mode)


def angular_distribution_errors(
    z: int,
    form_table: Tab1,
    form_curve: EncodedCurve,
    scatter_table: Tab1,
    scatter_curve: EncodedCurve,
) -> tuple[float, float]:
    maximum_cdf_error = 0.0
    maximum_mean_cos_error = 0.0
    energies = (20.0, 50.0, 100.0, 500.0, 1000.0, 3000.0, 10000.0, 20000.0)
    panels = 1024
    encoded_curves: dict[str, tuple[list[tuple[float, float]], list[float]]] = {}
    for kind, curve, mode in (
        ("rayleigh", form_curve, "log_log"),
        ("compton", scatter_curve, "log_linear"),
    ):
        points = []
        for log_argument, code in zip(curve.argument, curve.packed_value):
            stored = (
                LOG_FLOOR if mode in ("log_log", "log_linear") and code == 0 else
                curve.value_offset + (code - 1) * curve.value_scale
                if mode in ("log_log", "log_linear") else
                curve.value_offset + code * curve.value_scale
            )
            value = 10.0**stored if mode in ("log_log", "log_linear") and stored > LOG_FLOOR else stored
            if mode in ("log_log", "log_linear") and stored <= LOG_FLOOR:
                value = 0.0
            points.append((10.0**log_argument, value))
        encoded_curves[kind] = (points, [point[0] for point in points])
    for energy in energies:
        for kind, table, curve in (
            ("rayleigh", form_table, form_curve),
            ("compton", scatter_table, scatter_curve),
        ):
            true_weights: list[float] = []
            represented_weights: list[float] = []
            mus: list[float] = []
            alpha = energy / MEC2_KEV
            for index in range(panels + 1):
                s = index / panels
                mu = 1.0 - 2.0 * s * s
                # Compare the stored 0.01 inv-A clamp against the underlying
                # EPDL behavior down to the former 1e-7 lower bound. This makes
                # the deliberate low-x truncation part of the angular gate.
                x = max(ANGULAR_REFERENCE_MIN, (energy / HC_KEV_ANGSTROM) * s)
                true_factor = table.evaluate(x)
                mode = "log_log" if kind == "rayleigh" else "log_linear"
                encoded_points, encoded_arguments = encoded_curves[kind]
                represented_factor = encoded_value_cached(
                    encoded_points, encoded_arguments, x, mode
                )
                if kind == "rayleigh":
                    polarization = 0.5 * (1.0 + mu * mu)
                    true = true_factor * true_factor * polarization * 4.0 * s
                    represented = represented_factor * represented_factor * polarization * 4.0 * s
                else:
                    epsilon = 1.0 / (1.0 + alpha * (1.0 - mu))
                    kn = epsilon * epsilon * (epsilon + 1.0 / epsilon - (1.0 - mu * mu))
                    true = true_factor * kn * 4.0 * s
                    represented = represented_factor * kn * 4.0 * s
                true_weights.append(true)
                represented_weights.append(represented)
                mus.append(mu)

            def distribution(weights: list[float]) -> tuple[list[float], float]:
                cumulative = [0.0]
                moment = 0.0
                for index in range(panels):
                    area = 0.5 * (weights[index] + weights[index + 1]) / panels
                    cumulative.append(cumulative[-1] + area)
                    moment += 0.5 * (
                        weights[index] * mus[index] + weights[index + 1] * mus[index + 1]
                    ) / panels
                total = cumulative[-1]
                return [value / total for value in cumulative], moment / total

            true_cdf, true_mean = distribution(true_weights)
            represented_cdf, represented_mean = distribution(represented_weights)
            maximum_cdf_error = max(
                maximum_cdf_error,
                max(abs(a - b) for a, b in zip(true_cdf, represented_cdf)),
            )
            maximum_mean_cos_error = max(maximum_mean_cos_error, abs(true_mean - represented_mean))
    return maximum_cdf_error, maximum_mean_cos_error


def rayleigh_runtime_cdf(
    curve: EncodedCurve,
) -> tuple[list[float], list[int], float, float]:
    """Build the packed runtime inverse CDF and measure its representation error."""
    count = RAYLEIGH_CDF_N
    arguments = [
        as_float32(10.0 ** (-3.0 + index * 7.0 / (count - 1)))
        for index in range(count)
    ]
    x0 = arguments[0]
    f0 = curve_runtime(curve, x0, "log_log")
    cumulative = [f0 * f0 * x0 * x0]
    previous_x = x0
    previous_integrand = f0 * f0 * 2.0 * x0
    for argument in arguments[1:]:
        factor = curve_runtime(curve, argument, "log_log")
        integrand = factor * factor * 2.0 * argument
        cumulative.append(
            cumulative[-1]
            + 0.5 * (previous_integrand + integrand) * (argument - previous_x)
        )
        previous_x, previous_integrand = argument, integrand
    source_total = cumulative[-1]
    packed = [
        min(65535, max(0, round(value / source_total * 65535.0)))
        for value in cumulative
    ]
    packed[-1] = 65535
    storage_error = max(
        abs(source / source_total - stored / 65535.0)
        for source, stored in zip(cumulative, packed)
    )

    # Compare both CDF interpolations used by the sampler with a much denser
    # trapezoidal integration of the same decoded form-factor curve. Linear-x
    # interpolation is implied by inverse sampling within a grid cell, while
    # log-x interpolation determines G(x_max) for a truncated distribution.
    dense_count = 8193
    dense_arguments = [
        10.0 ** (-3.0 + index * 7.0 / (dense_count - 1))
        for index in range(dense_count)
    ]
    dense_x0 = dense_arguments[0]
    dense_f0 = curve_runtime(curve, dense_x0, "log_log")
    dense_cumulative = [dense_f0 * dense_f0 * dense_x0 * dense_x0]
    previous_x = dense_x0
    previous_integrand = dense_f0 * dense_f0 * 2.0 * dense_x0
    for argument in dense_arguments[1:]:
        factor = curve_runtime(curve, argument, "log_log")
        integrand = factor * factor * 2.0 * argument
        dense_cumulative.append(
            dense_cumulative[-1]
            + 0.5 * (previous_integrand + integrand) * (argument - previous_x)
        )
        previous_x, previous_integrand = argument, integrand
    dense_total = dense_cumulative[-1]
    discretization_error = 0.0
    for argument, reference in zip(dense_arguments, dense_cumulative):
        position = (math.log10(argument) + 3.0) * (count - 1) / 7.0
        low = min(count - 2, max(0, int(math.floor(position))))
        log_fraction = position - low
        linear_fraction = (
            (argument - arguments[low]) / (arguments[low + 1] - arguments[low])
        )
        for fraction in (log_fraction, linear_fraction):
            represented = (
                packed[low] + fraction * (packed[low + 1] - packed[low])
            ) / 65535.0
            discretization_error = max(
                discretization_error,
                abs(represented - reference / dense_total),
            )
    return arguments, packed, storage_error, discretization_error


def emit_external_array(out, ctype: str, name: str, values: list, formatter,
                        width: int) -> None:
    out.write(f"extern const {ctype} {name}[{len(values)}] = {{\n")
    out.write(format_rows(values, width, "    ", formatter))
    out.write("\n};\n\n")


def emit_external_uint16_rows(
    out, name: str, rows: list[list[int]], width: int,
) -> None:
    out.write(
        f"extern const uint16_t {name}[{len(rows)}][{len(rows[0])}] = {{\n"
    )
    for index, values in enumerate(rows):
        out.write("    {\n")
        out.write(format_rows(values, width, "        ", str))
        out.write("\n    }" + ("," if index + 1 < len(rows) else "") + "\n")
    out.write("};\n\n")


def compact_storage(
    processes: dict[int, dict[str, EncodedCurve]],
    angular: dict[int, dict[str, EncodedCurve]],
    rayleigh_sampling_x: list[float],
    rayleigh_cdf_q: list[list[int]],
) -> CompactStorage:
    rayleigh_log_energy: list[float] = []
    rayleigh_value: list[int] = []
    rayleigh_grid_offset = [0]
    rayleigh_value_offset = [0]
    for first_z in range(1, 93, RAYLEIGH_XS_GROUP_SIZE):
        last_z = min(first_z + RAYLEIGH_XS_GROUP_SIZE, 93)
        group = [processes[z]["rayleigh"] for z in range(first_z, last_z)]
        arguments = group[0].argument
        if any(curve.argument != arguments for curve in group[1:]):
            raise ValueError("Rayleigh curves in a storage group do not share a grid")
        rayleigh_log_energy.extend(
            as_float32(argument - 6.0) for argument in arguments
        )
        for curve in group:
            rayleigh_value.extend(curve.packed_value)
        rayleigh_grid_offset.append(len(rayleigh_log_energy))
        rayleigh_value_offset.append(len(rayleigh_value))

    # Photon coordinates are stored only once and addressed by exact
    # float32-bit identity. No contiguous raw-grid compatibility slices remain.
    energy_pool: list[float] = []
    pool_index: dict[bytes, int] = {}
    for name, _, _, _ in PROCESS_SPECS:
        if name == "rayleigh":
            continue
        for z in range(1, 93):
            for argument in processes[z][name].argument:
                value = as_float32(argument - 6.0)
                key = struct.pack("<f", value)
                if key not in pool_index:
                    pool_index[key] = len(energy_pool)
                    energy_pool.append(value)

    process_grid_index: list[int] = []
    process_value: list[int] = []
    process_offset: dict[tuple[int, str], int] = {}
    for name, _, _, _ in PROCESS_SPECS:
        if name == "rayleigh":
            continue
        for z in range(1, 93):
            curve = processes[z][name]
            process_offset[z, name] = len(process_value)
            for argument in curve.argument:
                value = as_float32(argument - 6.0)
                process_grid_index.append(pool_index[struct.pack("<f", value)])
            process_value.extend(curve.packed_value)

    angular_log_x: list[float] = []
    angular_value: list[int] = []
    angular_offset: dict[tuple[int, str], int] = {}
    for name, _, _, _ in RUNTIME_ANGULAR_SPECS:
        for z in range(1, 93):
            curve = angular[z][name]
            angular_offset[z, name] = len(angular_value)
            angular_log_x.extend(curve.argument)
            angular_value.extend(curve.packed_value)
    limits = {
        "Rayleigh grid offset": len(rayleigh_log_energy),
        "Rayleigh value offset": len(rayleigh_value),
        "energy-pool index": len(energy_pool),
        "process data offset": len(process_value),
        "angular data offset": len(angular_value),
    }
    for label, count in limits.items():
        if count > 65535:
            raise ValueError(f"{label} count {count} exceeds uint16 capacity")
    if len(process_grid_index) != len(process_value):
        raise ValueError("process grid/value streams are not parallel")
    if len(angular_log_x) != len(angular_value):
        raise ValueError("angular grid/value streams are not parallel")
    if len(rayleigh_sampling_x) != RAYLEIGH_CDF_N:
        raise ValueError("Rayleigh sampling grid has the wrong size")
    if len(rayleigh_cdf_q) != 92 or any(
        len(row) != RAYLEIGH_CDF_N for row in rayleigh_cdf_q
    ):
        raise ValueError("Rayleigh sampling CDF has the wrong shape")
    if any(
        row[-1] != 65535
        or any(low > high for low, high in zip(row, row[1:]))
        for row in rayleigh_cdf_q
    ):
        raise ValueError("Rayleigh sampling CDF is not normalized and monotonic")
    return CompactStorage(
        rayleigh_log_energy, rayleigh_value,
        rayleigh_grid_offset, rayleigh_value_offset,
        energy_pool, process_grid_index, process_value, process_offset,
        angular_log_x, angular_value, angular_offset,
        rayleigh_sampling_x, rayleigh_cdf_q,
    )


def emit_cpp(
    output: Path,
    source: dict,
    processes: dict[int, dict[str, EncodedCurve]],
    angular: dict[int, dict[str, EncodedCurve]],
    storage: CompactStorage,
) -> None:
    with output.open("w", encoding="utf-8") as out:
        out.write("""/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC.
 This library is distributed under the GNU Lesser General Public License,
 version 2.1 or (at your option) any later version.
 */

/// AUTO-GENERATED by tools/prepare_cross_sections/generate_epics_photon.py.
/// Direct source: EPICS2023 EPDL, verified by sources.lock.json.
/// Native nodes from 10 keV--10 MeV are accuracy probes; adaptive reduction
/// retains endpoints, thresholds, edge flanks, and only the additional points
/// needed by the committed gates. Values use per-curve uint16 packing; the
/// Rayleigh total cross sections share one adaptive energy grid per four
/// adjacent elements; other photon processes use bit-exact pooled coordinates.
/// The normalized Rayleigh sampling CDF is generated offline as uint16, and its
/// source form-factor curve is not emitted at runtime. No Geant4/G4EMLOW data
/// is read or emitted.

#include "cross_sections/photon_epics_data.h"

namespace ceelo {

""")
        emit_external_array(
            out, "float", "g_rayleigh_log_energy", storage.rayleigh_log_energy,
            format_float32, 8,
        )
        emit_external_array(
            out, "uint16_t", "g_rayleigh_log_value_q", storage.rayleigh_value,
            str, 16,
        )
        emit_external_array(
            out, "uint16_t", "g_rayleigh_group_grid_offset",
            storage.rayleigh_grid_offset, str, 16,
        )
        emit_external_array(
            out, "uint16_t", "g_rayleigh_group_value_offset",
            storage.rayleigh_value_offset, str, 16,
        )
        emit_external_array(
            out, "float", "g_photon_energy_pool", storage.energy_pool,
            format_float32, 8,
        )
        emit_external_array(
            out, "uint16_t", "g_photon_process_grid_index",
            storage.process_grid_index, str, 16,
        )
        emit_external_array(
            out, "uint16_t", "g_photon_process_log_value_q",
            storage.process_value, str, 16,
        )
        emit_external_array(
            out, "float", "g_photon_angular_log_x", storage.angular_log_x,
            format_float32, 8,
        )
        emit_external_array(
            out, "uint16_t", "g_photon_angular_log_value_q",
            storage.angular_value, str, 16,
        )
        emit_external_array(
            out, "float", "g_rayleigh_sampling_x",
            storage.rayleigh_sampling_x, format_float32, 8,
        )
        emit_external_uint16_rows(
            out, "g_rayleigh_sampling_cdf_q", storage.rayleigh_cdf_q, 16,
        )

        out.write("extern const PhotonEpicsElementData g_photon_epics_data[92] = {\n")
        for z in range(1, 93):
            rayleigh = processes[z]["rayleigh"]
            refs = [
                f"{{{format_float32(rayleigh.value_offset)}, "
                f"{format_float32(rayleigh.value_scale)}}}"
            ]
            for name, _, _, _ in PROCESS_SPECS:
                if name == "rayleigh":
                    continue
                size = len(processes[z][name].argument)
                curve = processes[z][name]
                refs.append(
                    f"{{{size}, {storage.process_offset[z, name]}, "
                    f"{format_float32(curve.value_offset)}, "
                    f"{format_float32(curve.value_scale)}}}"
                )
            for name, _, _, _ in RUNTIME_ANGULAR_SPECS:
                size = len(angular[z][name].argument)
                curve = angular[z][name]
                refs.append(
                    f"{{{size}, {storage.angular_offset[z, name]}, "
                    f"{format_float32(curve.value_offset)}, "
                    f"{format_float32(curve.value_scale)}}}"
                )
            comma = "," if z < 92 else ""
            out.write(f"    {{{', '.join(refs)}}}{comma} // Z={z}\n")
        out.write("};\n\n} // namespace ceelo\n")


def report_rows(
    processes: dict[int, dict[str, EncodedCurve]],
    angular: dict[int, dict[str, EncodedCurve]],
    storage: CompactStorage,
) -> list[dict]:
    rows: list[dict] = []
    for name, _, _, _ in PROCESS_SPECS:
        curves = [processes[z][name] for z in range(1, 93)]
        counts = [len(curve.argument) for curve in curves]
        native = [curve.native_count for curve in curves]
        if name == "rayleigh":
            encoding = (
                f"one adaptive float32 log10(argument) grid per "
                f"{RAYLEIGH_XS_GROUP_SIZE} adjacent elements, per-curve affine "
                "uint16 log10(value)"
            )
            grid_bytes = len(storage.rayleigh_log_energy) * 4
            metadata_bytes = (
                92 * 8
                + len(storage.rayleigh_grid_offset) * 2
                + len(storage.rayleigh_value_offset) * 2
            )
        else:
            encoding = "shared float32 log10(argument) pool with uint16 indices, per-curve affine uint16 log10(value)"
            grid_bytes = sum(counts) * 2
            metadata_bytes = 92 * 12
        rows.append({
            "process": name,
            "runtime_emitted": True,
            "encoding": encoding,
            "native_coefficients": sum(native),
            "retained_coefficients": sum(counts),
            "coefficient_bytes": sum(counts) * 2,
            "grid_bytes": grid_bytes,
            "metadata_bytes": metadata_bytes,
            "legacy_compatibility_bytes": 0,
            "nodes_min": min(counts),
            "nodes_median": percentile(counts, 0.5),
            "nodes_p95": percentile(counts, 0.95),
            "nodes_max": max(counts),
            "p95_interpolation_error_percent": percentile([c.p95_error for c in curves], 0.95) * 100.0,
            "maximum_interpolation_error_percent": max(c.maximum_error for c in curves) * 100.0,
            "maximum_integrated_error_percent": max(c.integrated_error for c in curves) * 100.0,
        })
    for name, _, _, mode in ANGULAR_SPECS:
        curves = [angular[z][name] for z in range(1, 93)]
        counts = [len(curve.argument) for curve in curves]
        runtime_emitted = name != "form_factor"
        runtime_counts = counts if runtime_emitted else [0] * len(counts)
        rows.append({
            "process": name,
            "runtime_emitted": runtime_emitted,
            "encoding": (
                "generator-only adaptive EPDL curve; not emitted at runtime"
                if not runtime_emitted else
                "float32 log10(argument), per-curve affine uint16 log10(value)"
                if mode == "log_log"
                else "float32 log10(argument), per-curve affine uint16 log10(value), linear-value interpolation"
            ),
            "native_coefficients": sum(curve.native_count for curve in curves),
            "retained_coefficients": sum(runtime_counts),
            "coefficient_bytes": sum(runtime_counts) * 2,
            "grid_bytes": sum(runtime_counts) * 4,
            "metadata_bytes": 92 * 12 if runtime_emitted else 0,
            "legacy_compatibility_bytes": 0,
            "nodes_min": min(runtime_counts),
            "nodes_median": percentile(runtime_counts, 0.5),
            "nodes_p95": percentile(runtime_counts, 0.95),
            "nodes_max": max(runtime_counts),
            "p95_interpolation_error_percent": percentile([c.p95_error for c in curves], 0.95) * 100.0,
            "maximum_interpolation_error_percent": max(c.maximum_error for c in curves) * 100.0,
            "maximum_integrated_error_percent": max(c.integrated_error for c in curves) * 100.0,
        })
    return rows


def write_reports(
    directory: Path,
    lock_path: Path,
    source: dict,
    rows: list[dict],
    storage: CompactStorage,
    angular_cdf_error: float,
    angular_mean_cos_error: float,
    rayleigh_cdf_quantization_error: float,
    rayleigh_cdf_discretization_error: float,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "source": "EPICS2023 EPDL",
        "source_sha256": source["sha256"],
        "source_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "probe_positions_per_interval": PROBES_PER_INTERVAL,
        "maximum_allowed_interpolation_error_percent": MAX_RELATIVE_ERROR * 100.0,
        "maximum_allowed_integrated_error_percent": MAX_INTEGRATED_ERROR * 100.0,
        "maximum_allowed_angular_cdf_error_percent": MAX_ANGULAR_CDF_ERROR * 100.0,
        "maximum_allowed_mean_cosine_error": MAX_ANGULAR_MEAN_COS_ERROR,
        "maximum_angular_cdf_error_percent": angular_cdf_error * 100.0,
        "maximum_mean_cosine_error": angular_mean_cos_error,
        "maximum_rayleigh_runtime_cdf_uint16_error": rayleigh_cdf_quantization_error,
        "maximum_allowed_rayleigh_runtime_cdf_uint16_error": MAX_RAYLEIGH_CDF_QUANTIZATION_ERROR,
        "rayleigh_runtime_cdf_nodes": RAYLEIGH_CDF_N,
        "maximum_rayleigh_runtime_cdf_discretization_error_percent": rayleigh_cdf_discretization_error * 100.0,
        "maximum_allowed_rayleigh_runtime_cdf_discretization_error_percent": MAX_RAYLEIGH_CDF_DISCRETIZATION_ERROR * 100.0,
        "shared_energy_coordinate_pool": {
            "stored_float32_coordinates": len(storage.energy_pool),
            "bytes": len(storage.energy_pool) * 4,
            "note": "Every coordinate is bit-exact deduplicated; no raw legacy grid views are retained.",
        },
        "rayleigh_cross_section_grids": {
            "elements_per_group": RAYLEIGH_XS_GROUP_SIZE,
            "groups": len(storage.rayleigh_grid_offset) - 1,
            "stored_float32_coordinates": len(storage.rayleigh_log_energy),
            "bytes": len(storage.rayleigh_log_energy) * 4,
            "offset_bytes": (
                len(storage.rayleigh_grid_offset)
                + len(storage.rayleigh_value_offset)
            ) * 2,
            "note": "Each adjacent-Z group uses one jointly refined grid; every element and probe must pass the unchanged gates.",
        },
        "rayleigh_sampling_cdf": {
            "nodes_per_element": RAYLEIGH_CDF_N,
            "elements": 92,
            "encoding": "normalized uint16",
            "coefficient_bytes": 92 * RAYLEIGH_CDF_N * 2,
            "shared_grid_bytes": RAYLEIGH_CDF_N * 4,
            "generated_offline": True,
        },
        "processes": rows,
    }
    (directory / "photon_epics.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    with (directory / "photon_epics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# EPICS2023 process-specific photon table report\n",
        "Every native source node in 10 keV–10 MeV is validated, but only endpoints, "
        "thresholds, duplicated-edge flank/edge pairs, and adaptively required nodes are retained. "
        "The reported errors include uint16 value decoding and 16 logarithmic probes per interval.\n",
        "The form-factor curve is generator-only: its error columns validate the curve used to build "
        "the offline Rayleigh CDF, while its runtime retained-node and byte counts are zero.\n",
        f"The shared energy-coordinate pool stores {len(storage.energy_pool):,} float32 values "
        f"({len(storage.energy_pool) * 4:,} bytes); process grid columns below are uint16 pool indices.\n",
        f"Rayleigh total cross sections use {len(storage.rayleigh_grid_offset) - 1} shared "
        f"adjacent-Z adaptive grids containing {len(storage.rayleigh_log_energy):,} total "
        f"float32 coordinates ({len(storage.rayleigh_log_energy) * 4:,} bytes).\n",
        "| Process | Native | Retained | Coeff bytes | Grid bytes | Nodes min/med/p95/max | p95/max error | Max integrated |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['process']} | {row['native_coefficients']:,} | {row['retained_coefficients']:,} | "
            f"{row['coefficient_bytes']:,} | {row['grid_bytes']:,} | "
            f"{row['nodes_min']}/{row['nodes_median']:.0f}/{row['nodes_p95']:.0f}/{row['nodes_max']} | "
            f"{row['p95_interpolation_error_percent']:.5f}%/{row['maximum_interpolation_error_percent']:.5f}% | "
            f"{row['maximum_integrated_error_percent']:.5f}% |"
        )
    lines.extend([
        "\nAngular distribution gates:",
        f"maximum CDF error = {angular_cdf_error * 100.0:.6f}% "
        f"(limit {MAX_ANGULAR_CDF_ERROR * 100.0:.3g}%); maximum mean-cosine error = "
        f"{angular_mean_cos_error:.7f} (limit {MAX_ANGULAR_MEAN_COS_ERROR:.4g}).",
        f"Offline Rayleigh inverse-CDF uint16 quantization error = "
        f"{rayleigh_cdf_quantization_error:.3g} "
        f"(limit {MAX_RAYLEIGH_CDF_QUANTIZATION_ERROR:.3g}).",
        f"Runtime Rayleigh {RAYLEIGH_CDF_N}-node inverse-CDF discretization error = "
        f"{rayleigh_cdf_discretization_error * 100.0:.6f}% "
        f"(limit {MAX_RAYLEIGH_CDF_DISCRETIZATION_ERROR * 100.0:.3g}%).",
    ])
    (directory / "photon_epics.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epdl", type=Path, required=True)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    args = parser.parse_args()

    lock = load_lock(args.lock)
    source = lock["sources"]["epdl2023"]
    actual = sha256_file(args.epdl)
    if actual != source["sha256"]:
        raise ValueError(f"EPDL SHA-256 mismatch: expected {source['sha256']}, got {actual}")
    epdl = EndfFile(args.epdl)
    processes: dict[int, dict[str, EncodedCurve]] = {}
    angular: dict[int, dict[str, EncodedCurve]] = {}
    angular_cdf_error = 0.0
    angular_mean_cos_error = 0.0
    rayleigh_cdf_quantization_error = 0.0
    rayleigh_cdf_discretization_error = 0.0
    rayleigh_sampling_x: list[float] = []
    rayleigh_cdf_q: list[list[int]] = []

    for z in range(1, 93):
        processes[z] = {}
        for name, mf, mt, mode in PROCESS_SPECS:
            table = epdl.tab1(z, mf, mt)
            if name == "rayleigh":
                continue
            if name == "pair_production":
                lower = max(PAIR_THRESHOLD_EV, table.x[0])
            elif name == "k_photoelectric":
                lower = max(ENERGY_MIN_EV, table.x[0])
            else:
                lower = ENERGY_MIN_EV
            processes[z][name] = build_curve(table, lower, ENERGY_MAX_EV, mode)
        angular[z] = {}
        angular_tables = {}
        for name, mf, mt, mode in ANGULAR_SPECS:
            table = epdl.tab1(z, mf, mt)
            angular_tables[name] = table
            angular[z][name] = build_curve(table, ANGULAR_MIN, ANGULAR_MAX, mode)
        cdf_error, mean_error = angular_distribution_errors(
            z,
            angular_tables["form_factor"], angular[z]["form_factor"],
            angular_tables["scattering_function"], angular[z]["scattering_function"],
        )
        angular_cdf_error = max(angular_cdf_error, cdf_error)
        angular_mean_cos_error = max(angular_mean_cos_error, mean_error)
        sampling_x, cdf_q, storage_error, discretization_error = rayleigh_runtime_cdf(
            angular[z]["form_factor"]
        )
        if not rayleigh_sampling_x:
            rayleigh_sampling_x = sampling_x
        elif sampling_x != rayleigh_sampling_x:
            raise ValueError("Rayleigh sampling grids differ between elements")
        rayleigh_cdf_q.append(cdf_q)
        rayleigh_cdf_quantization_error = max(
            rayleigh_cdf_quantization_error, storage_error
        )
        rayleigh_cdf_discretization_error = max(
            rayleigh_cdf_discretization_error, discretization_error
        )
        print(f"\r  generated direct EPDL Z={z:3d}/92", end="", file=sys.stderr)
    print(file=sys.stderr)
    for first_z in range(1, 93, RAYLEIGH_XS_GROUP_SIZE):
        last_z = min(first_z + RAYLEIGH_XS_GROUP_SIZE, 93)
        curves = build_shared_rayleigh_group([
            epdl.tab1(z, 23, 502) for z in range(first_z, last_z)
        ])
        for z, curve in zip(range(first_z, last_z), curves):
            processes[z]["rayleigh"] = curve
        print(
            f"\r  generated shared Rayleigh groups through Z={last_z - 1:3d}/92",
            end="", file=sys.stderr,
        )
    print(file=sys.stderr)
    if (angular_cdf_error > MAX_ANGULAR_CDF_ERROR
            or angular_mean_cos_error > MAX_ANGULAR_MEAN_COS_ERROR):
        raise ValueError(
            f"angular gates failed: CDF={angular_cdf_error:.6%}, mean-cos={angular_mean_cos_error:.6g}"
        )
    if rayleigh_cdf_quantization_error > MAX_RAYLEIGH_CDF_QUANTIZATION_ERROR:
        raise ValueError(
            "Rayleigh uint16 CDF quantization gate failed: "
            f"{rayleigh_cdf_quantization_error:.6g}"
        )
    if rayleigh_cdf_discretization_error > MAX_RAYLEIGH_CDF_DISCRETIZATION_ERROR:
        raise ValueError(
            "Rayleigh inverse-CDF discretization gate failed: "
            f"{rayleigh_cdf_discretization_error:.6g}"
        )

    storage = compact_storage(
        processes, angular, rayleigh_sampling_x, rayleigh_cdf_q
    )
    emit_cpp(args.output, source, processes, angular, storage)
    rows = report_rows(processes, angular, storage)
    write_reports(
        args.report_dir, args.lock, source, rows, storage, angular_cdf_error,
        angular_mean_cos_error, rayleigh_cdf_quantization_error,
        rayleigh_cdf_discretization_error,
    )
    print(f"Wrote {args.output} and photon_epics CSV/JSON/Markdown reports", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
