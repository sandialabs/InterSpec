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

"""End-to-end Python/C++ parity for one serialized runtime candidate."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence

from .fit import PAIR_THRESHOLD_KEV, RuntimePayload
from .publication import file_sha256


PARITY_REPORT_SCHEMA = "CEEFIT-PYTHON-CPP-PARITY-v1"
PARITY_REQUEST_SCHEMA = "CEEFIT-PARITY-REQUEST-v1"
PARITY_RESULT_SCHEMA = "CEEFIT-PARITY-RESULT-v1"
PARITY_TOLERANCE = 1.0e-9


def _atomic_json(path: str | Path, value: Mapping[str, Any]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.write_text(encoded, encoding="utf-8")
    temporary.replace(destination)


@dataclass(frozen=True)
class ParityQuery:
    label: str
    category: str
    energy_keV: float
    distance_m: float
    source_height_m: float
    detector_height_m: float
    direct_fluence_cm2: float
    energy_edges_keV: tuple[float, ...]
    mu_edges: tuple[float, ...]


def _strict_edges(values: Iterable[float]) -> tuple[float, ...]:
    result: list[float] = []
    for raw in sorted(float(value) for value in values if math.isfinite(float(value))):
        if raw < 0.0:
            continue
        if not result or raw > result[-1] + 1.0e-13 * max(1.0, abs(raw)):
            result.append(raw)
    if len(result) < 2:
        raise ValueError("parity query needs at least two distinct energy edges")
    return tuple(result)


def _energy_edges(payload: RuntimePayload, energy_keV: float) -> tuple[float, ...]:
    e_cuts = [float(model.e_cut_keV) for model in payload.components.values()]
    e_cut = min(e_cuts, default=10.0)
    return _strict_edges(
        (
            0.0,
            0.5 * e_cut,
            e_cut,
            0.073 * energy_keV,
            0.211 * energy_keV,
            170.3,
            0.499 * energy_keV,
            510.998950,
            0.811 * energy_keV,
            math.nextafter(energy_keV, 0.0),
            energy_keV,
            math.nextafter(energy_keV, math.inf),
            max(energy_keV + 1.0, 1.07 * energy_keV),
        )
    )


_PARITY_MU_EDGES = (-1.0, -0.91, -0.42, -0.03, 0.23, 0.77, 1.0)


def _direct_fluence(distance_m: float) -> float:
    distance_cm = 100.0 * distance_m
    return 1.0 / (4.0 * math.pi * distance_cm * distance_cm)


def deterministic_parity_queries(payload: RuntimePayload) -> tuple[ParityQuery, ...]:
    """Cover interior, physical faces, and both activation thresholds.

    Every proposed point is evaluated once in Python before inclusion.  This
    rejects geometrically invalid corners of the rectangular physical bounds,
    while a missing mandatory coverage category remains a hard error.
    """

    domain = payload.metadata.get("physical_domain")
    if not isinstance(domain, Mapping):
        raise ValueError("runtime payload lacks its physical domain")
    try:
        low = tuple(float(value) for value in domain["min"])
        high = tuple(float(value) for value in domain["max"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("runtime payload physical domain is malformed") from error
    if len(low) != 4 or len(high) != 4 or any(
        not math.isfinite(value) for value in (*low, *high)
    ) or any(upper < lower for lower, upper in zip(low, high)):
        raise ValueError("runtime payload physical domain is invalid")
    if low[0] <= 0.0 or low[1] <= 0.0:
        raise ValueError("runtime payload energy and distance domains must be positive")

    energy_center = math.sqrt(low[0] * high[0])
    distance_center = math.sqrt(low[1] * high[1])
    source_center = 0.5 * (low[2] + high[2])
    detector_center = 0.5 * (low[3] + high[3])
    if abs(source_center - detector_center) >= distance_center:
        distance_center = high[1]
    overlap_low = max(low[2], low[3])
    overlap_high = min(high[2], high[3])
    if overlap_low <= overlap_high:
        equal_source = equal_detector = 0.5 * (overlap_low + overlap_high)
    elif high[2] < low[3]:
        equal_source, equal_detector = high[2], low[3]
    else:
        equal_source, equal_detector = low[2], high[3]

    proposed: list[tuple[str, str, tuple[float, float, float, float]]] = [
        ("interior", "interior", (energy_center, distance_center, source_center, detector_center)),
        ("energy_min_face", "domain_face", (low[0], distance_center, source_center, detector_center)),
        ("energy_max_face", "domain_face", (high[0], distance_center, source_center, detector_center)),
        ("distance_min_face", "domain_face", (energy_center, low[1], equal_source, equal_detector)),
        ("distance_max_face", "domain_face", (energy_center, high[1], source_center, detector_center)),
        ("source_height_min_face", "domain_face", (energy_center, high[1], low[2], detector_center)),
        ("source_height_max_face", "domain_face", (energy_center, high[1], high[2], detector_center)),
        ("detector_height_min_face", "domain_face", (energy_center, high[1], source_center, low[3])),
        ("detector_height_max_face", "domain_face", (energy_center, high[1], source_center, high[3])),
    ]
    thresholds = (
        ("brems_threshold", float(payload.metadata["brems_activation_threshold_keV"])),
        ("pair_threshold", PAIR_THRESHOLD_KEV),
    )
    for category, threshold in thresholds:
        if low[0] <= threshold <= high[0]:
            delta = max(1.0e-8, 1.0e-10 * threshold)
            values = (
                ("below", max(low[0], threshold - delta)),
                ("exact", threshold),
                ("above", min(high[0], threshold + delta)),
            )
            for suffix, energy in values:
                proposed.append(
                    (
                        f"{category}_{suffix}",
                        category,
                        (energy, distance_center, source_center, detector_center),
                    )
                )

    result: list[ParityQuery] = []
    seen: set[tuple[float, float, float, float, str]] = set()
    for label, category, coordinates in proposed:
        key = (*coordinates, category)
        if key in seen:
            continue
        seen.add(key)
        energy, distance, source, detector = coordinates
        direct = _direct_fluence(distance)
        energy_edges = _energy_edges(payload, energy)
        try:
            payload.evaluate(
                energy,
                distance,
                source,
                detector,
                energy_edges,
                _PARITY_MU_EDGES,
                direct_fluence_cm2=direct,
            )
        except ValueError:
            continue
        result.append(
            ParityQuery(
                label,
                category,
                energy,
                distance,
                source,
                detector,
                direct,
                energy_edges,
                _PARITY_MU_EDGES,
            )
        )

    coverage = {category: 0 for category in (
        "interior", "domain_face", "brems_threshold", "pair_threshold"
    )}
    for query in result:
        coverage[query.category] += 1
    required = {
        "interior": 1,
        "domain_face": 4,
        "brems_threshold": 3 if low[0] <= thresholds[0][1] <= high[0] else 0,
        "pair_threshold": 3 if low[0] <= thresholds[1][1] <= high[0] else 0,
    }
    missing = [name for name, count in required.items() if coverage[name] < count]
    if missing:
        raise ValueError(
            "fitted q domain cannot cover mandatory parity queries: " + ", ".join(missing)
        )
    return tuple(result)


def _request_text(queries: Sequence[ParityQuery]) -> str:
    lines = [PARITY_REQUEST_SCHEMA, f"N {len(queries)}"]
    for index, query in enumerate(queries):
        lines.append(
            " ".join(
                (
                    "Q",
                    str(index),
                    *(format(value, ".17g") for value in (
                        query.energy_keV,
                        query.distance_m,
                        query.source_height_m,
                        query.detector_height_m,
                        query.direct_fluence_cm2,
                    )),
                    str(len(query.energy_edges_keV)),
                    str(len(query.mu_edges)),
                )
            )
        )
        lines.append("E " + " ".join(format(value, ".17g") for value in query.energy_edges_keV))
        lines.append("M " + " ".join(format(value, ".17g") for value in query.mu_edges))
    lines.append("END")
    return "\n".join(lines) + "\n"


class _Tokens:
    def __init__(self, text: str):
        self._values = text.split()
        self._index = 0

    def take(self, expected: str | None = None) -> str:
        if self._index >= len(self._values):
            raise ValueError("truncated C++ parity result")
        value = self._values[self._index]
        self._index += 1
        if expected is not None and value != expected:
            raise ValueError(f"expected C++ parity token {expected!r}, got {value!r}")
        return value

    def count(self, name: str, maximum: int = 1_000_000) -> int:
        try:
            value = int(self.take(), 10)
        except ValueError as error:
            raise ValueError(f"invalid C++ parity {name}") from error
        if value < 0 or value > maximum:
            raise ValueError(f"invalid C++ parity {name}")
        return value

    def number(self, name: str) -> float:
        try:
            value = float(self.take())
        except ValueError as error:
            raise ValueError(f"invalid C++ parity {name}") from error
        if not math.isfinite(value):
            raise ValueError(f"non-finite C++ parity {name}")
        return value

    def vector(self, name: str) -> list[float]:
        count = self.count(f"{name} length")
        return [self.number(name) for _ in range(count)]

    def finish(self) -> None:
        if self._index != len(self._values):
            raise ValueError("C++ parity result has trailing data")


def _parse_cpp_result(text: str, expected_queries: int) -> tuple[str, list[dict[str, Any]]]:
    tokens = _Tokens(text)
    tokens.take(PARITY_RESULT_SCHEMA)
    tokens.take("MODEL")
    model_hash = tokens.take()
    if len(model_hash) != 64 or any(character not in "0123456789abcdef" for character in model_hash):
        raise ValueError("C++ parity result has an invalid model hash")
    tokens.take("N")
    if tokens.count("query count", 4096) != expected_queries:
        raise ValueError("C++ parity result query count differs")
    results: list[dict[str, Any]] = []
    for index in range(expected_queries):
        tokens.take("Q")
        if tokens.count("query index", 4095) != index:
            raise ValueError("C++ parity result query order differs")
        n_energy = tokens.count("energy bin count", 4096)
        n_mu = tokens.count("mu bin count", 4096)
        n_components = tokens.count("component count", 4096)
        n_lines = tokens.count("line count", 4096)
        direct = tokens.number("direct fluence")
        tokens.take("TOTAL")
        total = tokens.vector("total")
        tokens.take("INCIDENT_TOTAL")
        incident_total = tokens.vector("incident total")
        if len(total) != n_energy * n_mu or len(incident_total) != len(total):
            raise ValueError("C++ parity total dimensions differ")
        components: dict[str, tuple[list[float], list[float]]] = {}
        for _ in range(n_components):
            tokens.take("COMP")
            name = tokens.take()
            if name in components:
                raise ValueError(f"C++ parity result repeats component {name}")
            components[name] = (tokens.vector(name), tokens.vector(name))
        lines: dict[str, tuple[float, list[float], list[float]]] = {}
        for _ in range(n_lines):
            tokens.take("LINE")
            name = tokens.take()
            if name in lines:
                raise ValueError(f"C++ parity result repeats line {name}")
            energy = tokens.number("line energy")
            lines[name] = (energy, tokens.vector(name), tokens.vector(name))
        tokens.take("ENDQ")
        results.append(
            {
                "n_energy": n_energy,
                "n_mu": n_mu,
                "direct": direct,
                "total": total,
                "incident_total": incident_total,
                "components": components,
                "lines": lines,
            }
        )
    tokens.take("END")
    tokens.finish()
    return model_hash, results


def _flatten(values: Any) -> list[float]:
    result: list[float] = []
    for value in values:
        if isinstance(value, (tuple, list)):
            result.extend(_flatten(value))
        else:
            number = float(value)
            if not math.isfinite(number):
                raise ValueError("Python evaluator produced a non-finite mass")
            result.append(number)
    return result


def _difference(left: Sequence[float], right: Sequence[float], name: str) -> tuple[float, int]:
    if len(left) != len(right):
        raise ValueError(f"Python/C++ {name} dimensions differ")
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0), len(left)


def generate_parity_report(
    candidate_path: str | Path,
    cpp_evaluator_path: str | Path,
    output_path: str | Path,
    *,
    tolerance: float = PARITY_TOLERANCE,
    timeout_seconds: float = 120.0,
) -> dict[str, Any]:
    """Evaluate one serialized candidate through both production implementations."""

    if not math.isfinite(tolerance) or tolerance <= 0.0 or tolerance > PARITY_TOLERANCE:
        raise ValueError("parity tolerance must be positive and no larger than 1e-9")
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0.0:
        raise ValueError("parity timeout must be finite and positive")
    candidate = Path(candidate_path)
    executable = Path(cpp_evaluator_path)
    output = Path(output_path)
    if output.resolve() in {candidate.resolve(), executable.resolve()}:
        raise ValueError("parity report must not overwrite its candidate or C++ evaluator")
    payload = RuntimePayload.load_binary(candidate, allow_candidate=True)
    model_hash = payload.compute_model_content_sha256()
    if payload.model_content_sha256 != model_hash:
        raise ValueError("candidate model hash is inconsistent")
    queries = deterministic_parity_queries(payload)
    request_text = _request_text(queries)
    try:
        process = subprocess.run(
            [str(executable), str(candidate)],
            input=request_text,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout_seconds,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError(f"cannot execute C++ fitted-payload evaluator: {error}") from error
    if process.returncode != 0:
        detail = process.stderr.strip() or f"exit status {process.returncode}"
        raise RuntimeError(f"C++ fitted-payload evaluator failed: {detail}")
    cpp_hash, cpp_results = _parse_cpp_result(process.stdout, len(queries))
    if cpp_hash != model_hash:
        raise ValueError("C++ evaluator loaded a different fitted model")

    maximum = 0.0
    compared = 0
    query_reports: list[dict[str, Any]] = []
    coverage: dict[str, int] = {}
    for query, cpp in zip(queries, cpp_results):
        python = payload.evaluate(
            query.energy_keV,
            query.distance_m,
            query.source_height_m,
            query.detector_height_m,
            query.energy_edges_keV,
            query.mu_edges,
            direct_fluence_cm2=query.direct_fluence_cm2,
        )
        query_maximum = abs(query.direct_fluence_cm2 - float(cpp["direct"]))
        query_compared = 1
        for name, left, right in (
            ("continuum total", _flatten(python["continuum_s_hat"]), cpp["total"]),
            ("incident continuum total", _flatten(python["continuum"]), cpp["incident_total"]),
        ):
            difference, count = _difference(left, right, name)
            query_maximum = max(query_maximum, difference)
            query_compared += count

        python_components = set(python["components_s_hat"])
        if python_components != set(cpp["components"]):
            raise ValueError("Python/C++ component inventories differ")
        for name in sorted(python_components):
            for label, left, right in (
                (f"{name} s_hat", _flatten(python["components_s_hat"][name]), cpp["components"][name][0]),
                (f"{name} incident", _flatten(python["components"][name]), cpp["components"][name][1]),
            ):
                difference, count = _difference(left, right, label)
                query_maximum = max(query_maximum, difference)
                query_compared += count

        python_lines = {item["name"]: item for item in python["lines"]}
        if set(python_lines) != set(cpp["lines"]):
            raise ValueError("Python/C++ exact-line inventories differ")
        for name in sorted(python_lines):
            item = python_lines[name]
            cpp_energy, cpp_s_hat, cpp_incident = cpp["lines"][name]
            query_maximum = max(query_maximum, abs(float(item["energy_keV"]) - cpp_energy))
            query_compared += 1
            for label, left, right in (
                (f"{name} s_hat", _flatten(item["s_hat_mu_bin_masses"]), cpp_s_hat),
                (f"{name} incident", _flatten(item["mu_bin_masses"]), cpp_incident),
            ):
                difference, count = _difference(left, right, label)
                query_maximum = max(query_maximum, difference)
                query_compared += count

        maximum = max(maximum, query_maximum)
        compared += query_compared
        coverage[query.category] = coverage.get(query.category, 0) + 1
        query_reports.append(
            {
                "label": query.label,
                "category": query.category,
                "coordinates": [
                    query.energy_keV,
                    query.distance_m,
                    query.source_height_m,
                    query.detector_height_m,
                ],
                "energy_edge_count": len(query.energy_edges_keV),
                "mu_edge_count": len(query.mu_edges),
                "compared_value_count": query_compared,
                "max_absolute_difference": query_maximum,
            }
        )

    report = {
        "schema": PARITY_REPORT_SCHEMA,
        "config_sha256": str(payload.metadata["config_sha256"]),
        "model_content_sha256": model_hash,
        "candidate_payload_sha256": file_sha256(candidate),
        "artifact_paths": {
            "candidate_payload": Path(
                os.path.relpath(candidate.resolve(), output.parent.resolve())
            ).as_posix(),
        },
        "protocol": PARITY_REQUEST_SCHEMA,
        "request_sha256": hashlib.sha256(request_text.encode("ascii")).hexdigest(),
        "cpp_result_sha256": hashlib.sha256(process.stdout.encode("utf-8")).hexdigest(),
        "tolerance": tolerance,
        "query_count": len(queries),
        "coverage": coverage,
        "compared_value_count": compared,
        "max_absolute_difference": maximum,
        "queries": query_reports,
        "passed": bool(compared > 0 and maximum <= tolerance),
    }
    _atomic_json(output, report)
    return report
