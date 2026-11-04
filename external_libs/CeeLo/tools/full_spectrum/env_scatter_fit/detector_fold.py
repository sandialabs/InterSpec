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

"""Generate identity-bound HPGe/NaI detector-fold validation evidence.

Detector responses use the dependency-safe NumPy ``.npz`` contract documented
in ``README.md``.  Continuum entries are cell kernels on incident energy and
mu.  Exact atoms use a separate line-energy grid and are linearly interpolated
at their exact energy; they are never assigned to continuum-bin centres.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from .archive import require_hdf5
from .campaign import load_manifest
from .fit import RuntimePayload
from .noise import (
    EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
    EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
    EXACT_LINE_RELATIVE_ERROR_FRACTION,
    EXACT_LINE_RELATIVE_NOISE_FRACTION,
    exact_line_materiality,
)
from .publication import archive_snapshot_sha256


RESPONSE_SCHEMA = "CEEFIT-DETECTOR-RESPONSE-v1"
REPORT_SCHEMA = "CEEFIT-DETECTOR-FOLD-VALIDATION-v1"
RESPONSE_SEMANTICS = "probability-per-incident-cell-mass-v1"
LINE_INTERPOLATION = "linear-exact-energy-v1"
GEOMETRY_SCHEMA = "CEEFIT-DETECTOR-GEOMETRY-v1"
PROVENANCE_SCHEMA = "CEEFIT-DETECTOR-RESPONSE-PROVENANCE-v1"
MU_DEFINITION = "cosine-to-detector-axis;mu=1=head-on-v1"
DETECTOR_NAMES = ("HPGe", "NaI")
_RESPONSE_KEYS = frozenset(
    {
        "schema",
        "detector_name",
        "response_semantics",
        "line_interpolation",
        "geometry_json",
        "provenance_json",
        "generator_name",
        "generator_version",
        "generator_sha256",
        "continuum_energy_edges_keV",
        "mu_edges",
        "spectrum_edges_keV",
        "continuum_response",
        "line_energy_grid_keV",
        "line_response",
    }
)
_LINE_NAMES = ("elastic_E0", "annihilation_511")
_FINITE_FAILURE = 1.0e100


@dataclass(frozen=True)
class DetectorResponse:
    """Validated detector response loaded from one exact NPZ byte stream."""

    detector_name: str
    continuum_energy_edges_keV: Any
    mu_edges: Any
    spectrum_edges_keV: Any
    continuum_response: Any
    line_energy_grid_keV: Any
    line_response: Any
    geometry: Mapping[str, Any]
    provenance: Mapping[str, Any]
    generator_name: str
    generator_version: str
    generator_sha256: str
    sha256: str


def _scalar_text(array: Any, name: str) -> str:
    if array.shape != () or array.dtype.kind not in {"U", "S"}:
        raise ValueError(f"detector response {name} must be a scalar string")
    value = array.item()
    if isinstance(value, bytes):
        try:
            value = value.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError(f"detector response {name} is not UTF-8") from error
    return str(value)


def _sha256_text(value: str, name: str) -> str:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"detector response {name} must be a lowercase SHA-256")
    return value


def _canonical_json_mapping(text: str, name: str) -> Mapping[str, Any]:
    try:
        record = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError(f"detector response {name} is not valid JSON") from error
    if not isinstance(record, Mapping):
        raise ValueError(f"detector response {name} must encode a JSON object")
    canonical = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if text != canonical:
        raise ValueError(f"detector response {name} must use canonical JSON encoding")
    return record


def _validate_response_metadata(
    geometry: Mapping[str, Any],
    provenance: Mapping[str, Any],
    *,
    detector_name: str,
    generator_name: str,
    generator_version: str,
    generator_sha256: str,
) -> None:
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
    if set(geometry) != geometry_keys or geometry.get("schema") != GEOMETRY_SCHEMA:
        raise ValueError("detector response geometry_json has an invalid schema/inventory")
    if geometry.get("detector_name") != detector_name:
        raise ValueError("detector response geometry names another detector")
    expected_material = {"HPGe": "HPGe", "NaI": "NaI(Tl)"}[detector_name]
    if geometry.get("active_material") != expected_material:
        raise ValueError(
            f"{detector_name} response active material must be {expected_material!r}"
        )
    if geometry.get("mu_definition") != MU_DEFINITION:
        raise ValueError("detector response geometry uses an incompatible mu definition")
    for name in ("geometry_id", "description", "source_surface"):
        if not isinstance(geometry.get(name), str) or not str(geometry[name]).strip():
            raise ValueError(f"detector response geometry {name} must be nonempty")
    dimensions = geometry.get("dimensions_cm")
    if not isinstance(dimensions, Mapping) or not dimensions:
        raise ValueError("detector response geometry dimensions_cm must be a nonempty object")
    for name, raw_value in dimensions.items():
        if not isinstance(name, str) or not name or isinstance(raw_value, bool):
            raise ValueError("detector response geometry dimension is malformed")
        value = float(raw_value)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("detector response geometry dimensions must be finite and positive")

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
    if set(provenance) != provenance_keys or provenance.get("schema") != PROVENANCE_SCHEMA:
        raise ValueError("detector response provenance_json has an invalid schema/inventory")
    for name in ("transport_code", "transport_version", "rng", "created_utc"):
        if not isinstance(provenance.get(name), str) or not str(provenance[name]).strip():
            raise ValueError(f"detector response provenance {name} must be nonempty")
    histories = provenance.get("histories_per_grid_point")
    if isinstance(histories, bool) or not isinstance(histories, int) or histories <= 0:
        raise ValueError("detector response provenance histories_per_grid_point is invalid")
    physics = provenance.get("physics_configuration")
    if not isinstance(physics, Mapping) or not physics:
        raise ValueError("detector response provenance physics_configuration is empty")
    hashes = provenance.get("input_sha256")
    if not isinstance(hashes, Mapping) or not hashes:
        raise ValueError("detector response provenance input_sha256 is empty")
    for name, value in hashes.items():
        if not isinstance(name, str) or not name or not isinstance(value, str):
            raise ValueError("detector response provenance input hash is malformed")
        _sha256_text(value, f"provenance input_sha256/{name}")
    if not generator_name.strip() or not generator_version.strip():
        raise ValueError("detector response generator name/version must be nonempty")
    _sha256_text(generator_sha256, "generator_sha256")


def _strict_edges(numpy: Any, values: Any, name: str, *, positive: bool = False) -> Any:
    result = numpy.asarray(values)
    if result.dtype != numpy.dtype("float64"):
        raise ValueError(f"detector response {name} must have dtype float64")
    if result.ndim != 1 or result.size < 2:
        raise ValueError(f"detector response {name} must be a one-dimensional edge vector")
    if numpy.any(~numpy.isfinite(result)) or numpy.any(numpy.diff(result) <= 0.0):
        raise ValueError(f"detector response {name} must be finite and strictly increasing")
    if positive and float(result[0]) < 0.0:
        raise ValueError(f"detector response {name} cannot begin below zero")
    return result


def _response_tensor(numpy: Any, values: Any, shape: tuple[int, ...], name: str) -> Any:
    result = numpy.asarray(values)
    if result.dtype != numpy.dtype("float64"):
        raise ValueError(f"detector response {name} must have dtype float64")
    if result.shape != shape:
        raise ValueError(
            f"detector response {name} shape {result.shape} does not match {shape}"
        )
    if numpy.any(~numpy.isfinite(result)) or numpy.any(result < 0.0):
        raise ValueError(f"detector response {name} must be finite and nonnegative")
    efficiency = numpy.sum(result, axis=-1)
    if numpy.any(efficiency > 1.0 + 1.0e-12):
        raise ValueError(f"detector response {name} has probability greater than one")
    if not bool(numpy.any(efficiency > 0.0)):
        raise ValueError(f"detector response {name} contains no detector sensitivity")
    return result


def load_detector_response(
    path: str | Path,
    *,
    expected_detector_name: str,
) -> DetectorResponse:
    """Read and strictly validate one response from the bytes that are hashed."""

    if expected_detector_name not in DETECTOR_NAMES:
        raise ValueError(f"unsupported detector name {expected_detector_name!r}")
    numpy = require_hdf5()[1]
    raw = Path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    try:
        with numpy.load(io.BytesIO(raw), allow_pickle=False) as archive:
            keys = frozenset(archive.files)
            if keys != _RESPONSE_KEYS:
                missing = sorted(_RESPONSE_KEYS - keys)
                extra = sorted(keys - _RESPONSE_KEYS)
                raise ValueError(
                    f"detector response key inventory differs; missing={missing}, extra={extra}"
                )
            schema = _scalar_text(archive["schema"], "schema")
            detector_name = _scalar_text(archive["detector_name"], "detector_name")
            semantics = _scalar_text(archive["response_semantics"], "response_semantics")
            interpolation = _scalar_text(archive["line_interpolation"], "line_interpolation")
            geometry = _canonical_json_mapping(
                _scalar_text(archive["geometry_json"], "geometry_json"),
                "geometry_json",
            )
            provenance = _canonical_json_mapping(
                _scalar_text(archive["provenance_json"], "provenance_json"),
                "provenance_json",
            )
            generator_name = _scalar_text(archive["generator_name"], "generator_name")
            generator_version = _scalar_text(
                archive["generator_version"], "generator_version"
            )
            generator_sha256 = _scalar_text(
                archive["generator_sha256"], "generator_sha256"
            )
            energy_edges = _strict_edges(
                numpy,
                archive["continuum_energy_edges_keV"],
                "continuum_energy_edges_keV",
                positive=True,
            ).copy()
            mu_edges = _strict_edges(numpy, archive["mu_edges"], "mu_edges").copy()
            spectrum_edges = _strict_edges(
                numpy,
                archive["spectrum_edges_keV"],
                "spectrum_edges_keV",
                positive=True,
            ).copy()
            line_energy = _strict_edges(
                numpy,
                archive["line_energy_grid_keV"],
                "line_energy_grid_keV",
                positive=True,
            ).copy()
            n_energy = energy_edges.size - 1
            n_mu = mu_edges.size - 1
            n_spectrum = spectrum_edges.size - 1
            continuum = _response_tensor(
                numpy,
                archive["continuum_response"],
                (n_energy, n_mu, n_spectrum),
                "continuum_response",
            ).copy()
            lines = _response_tensor(
                numpy,
                archive["line_response"],
                (line_energy.size, n_mu, n_spectrum),
                "line_response",
            ).copy()
    except (OSError, EOFError) as error:
        raise ValueError(f"cannot decode detector response NPZ {path}") from error

    if schema != RESPONSE_SCHEMA:
        raise ValueError(f"unsupported detector response schema {schema!r}")
    if detector_name != expected_detector_name:
        raise ValueError(
            f"detector response names {detector_name!r}, expected {expected_detector_name!r}"
        )
    if semantics != RESPONSE_SEMANTICS:
        raise ValueError(f"unsupported detector response semantics {semantics!r}")
    if interpolation != LINE_INTERPOLATION:
        raise ValueError(f"unsupported exact-line interpolation {interpolation!r}")
    _validate_response_metadata(
        geometry,
        provenance,
        detector_name=detector_name,
        generator_name=generator_name,
        generator_version=generator_version,
        generator_sha256=generator_sha256,
    )
    if not (
        math.isclose(float(mu_edges[0]), -1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(float(mu_edges[-1]), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    ):
        raise ValueError("detector response mu edges must span exactly [-1, 1]")
    return DetectorResponse(
        detector_name,
        energy_edges,
        mu_edges,
        spectrum_edges,
        continuum,
        line_energy,
        lines,
        geometry,
        provenance,
        generator_name,
        generator_version,
        generator_sha256,
        digest,
    )


def _overlap_matrix(numpy: Any, source_edges: Any, target_edges: Any, name: str) -> Any:
    source = numpy.asarray(source_edges, dtype=float)
    target = numpy.asarray(target_edges, dtype=float)
    if source.ndim != 1 or source.size < 2 or numpy.any(~numpy.isfinite(source)) or numpy.any(
        numpy.diff(source) <= 0.0
    ):
        raise ValueError(f"archive {name} source edges are invalid")
    scale = max(1.0, abs(float(source[0])), abs(float(source[-1])))
    tolerance = 1.0e-12 * scale
    if float(target[0]) > float(source[0]) + tolerance or float(target[-1]) < float(
        source[-1]
    ) - tolerance:
        raise ValueError(f"detector response {name} does not cover the archive source axis")
    low = numpy.maximum(target[:-1, None], source[None, :-1])
    high = numpy.minimum(target[1:, None], source[None, 1:])
    overlap = numpy.maximum(0.0, high - low) / numpy.diff(source)[None, :]
    if not numpy.allclose(
        numpy.sum(overlap, axis=0), 1.0, rtol=0.0, atol=5.0e-12
    ):
        raise ValueError(f"detector response {name} loses archive cell mass")
    return overlap


def _rebin_2d(numpy: Any, masses: Any, energy_overlap: Any, mu_overlap: Any) -> Any:
    source = numpy.asarray(masses, dtype=float)
    if source.shape != (energy_overlap.shape[1], mu_overlap.shape[1]):
        raise ValueError("archive continuum tensor does not match its declared axes")
    result = energy_overlap @ source @ mu_overlap.T
    if not math.isclose(
        float(numpy.sum(result)),
        float(numpy.sum(source)),
        rel_tol=5.0e-12,
        abs_tol=1.0e-300,
    ):
        raise ValueError("continuum conservative rebinning failed mass closure")
    return result


def _rebin_replicas(numpy: Any, masses: Any, energy_overlap: Any, mu_overlap: Any) -> Any:
    source = numpy.asarray(masses, dtype=float)
    if source.ndim != 3 or source.shape[1:] != (
        energy_overlap.shape[1],
        mu_overlap.shape[1],
    ):
        raise ValueError("archive continuum replicas do not match their declared axes")
    result = numpy.einsum(
        "ei,rij,mj->rem", energy_overlap, source, mu_overlap, optimize=True
    )
    source_total = numpy.sum(source, axis=(1, 2))
    target_total = numpy.sum(result, axis=(1, 2))
    if not numpy.allclose(target_total, source_total, rtol=5.0e-12, atol=1.0e-300):
        raise ValueError("replica conservative rebinning failed mass closure")
    return result


def _fold_discretization_bounds(
    numpy: Any,
    source_masses: Any,
    energy_overlap: Any,
    mu_overlap: Any,
    response_tensor: Any,
    approximate_spectrum: Any,
) -> dict[str, float]:
    """Bound response variation hidden inside archive cells.

    The archive stores only each source-cell mass. If a response boundary
    splits that cell, conservative rebinning uses overlap fractions. The true
    sub-cell allocation can be any convex combination of the overlapping
    response kernels. The kernel diameter around the overlap-weighted kernel
    therefore gives a distribution-free L1/yield bound.
    """

    masses = numpy.asarray(source_masses, dtype=float)
    response = numpy.asarray(response_tensor, dtype=float)
    spectrum = numpy.asarray(approximate_spectrum, dtype=float)
    if masses.shape != (energy_overlap.shape[1], mu_overlap.shape[1]):
        raise ValueError("fold-discretization source masses do not match overlap axes")
    if response.shape[:2] != (energy_overlap.shape[0], mu_overlap.shape[0]):
        raise ValueError("fold-discretization response tensor does not match target axes")
    l1_absolute = 0.0
    yield_absolute = 0.0
    for source_energy in range(masses.shape[0]):
        target_energy = numpy.flatnonzero(energy_overlap[:, source_energy] > 0.0)
        for source_mu in range(masses.shape[1]):
            mass = float(masses[source_energy, source_mu])
            if mass <= 0.0:
                continue
            target_mu = numpy.flatnonzero(mu_overlap[:, source_mu] > 0.0)
            if target_energy.size == 1 and target_mu.size == 1:
                continue
            weights = (
                energy_overlap[target_energy, source_energy, None]
                * mu_overlap[target_mu, source_mu][None, :]
            )
            kernels = response[numpy.ix_(target_energy, target_mu)]
            average = numpy.sum(kernels * weights[:, :, None], axis=(0, 1))
            l1_distance = numpy.sum(numpy.abs(kernels - average), axis=2)
            efficiencies = numpy.sum(kernels, axis=2)
            average_efficiency = float(numpy.sum(average))
            l1_absolute += mass * float(numpy.max(l1_distance))
            yield_absolute += mass * float(
                numpy.max(numpy.abs(efficiencies - average_efficiency))
            )
    approximate_total = float(numpy.sum(spectrum))
    lower_total = approximate_total - yield_absolute
    if approximate_total <= 0.0 or lower_total <= 0.0:
        return {
            "l1_absolute": l1_absolute,
            "yield_absolute": yield_absolute,
            "total_relative_bound": _FINITE_FAILURE,
            "spectral_tv_bound": _FINITE_FAILURE,
        }
    return {
        "l1_absolute": l1_absolute,
        "yield_absolute": yield_absolute,
        "total_relative_bound": yield_absolute / lower_total,
        "spectral_tv_bound": min(
            1.0, (l1_absolute + yield_absolute) / (2.0 * approximate_total)
        ),
    }


def _line_discretization_bounds(
    numpy: Any,
    source_masses: Any,
    mu_overlap: Any,
    line_kernel: Any,
    approximate_spectrum: Any,
) -> dict[str, float]:
    masses = numpy.asarray(source_masses, dtype=float).reshape(-1)
    kernel = numpy.asarray(line_kernel, dtype=float)
    spectrum = numpy.asarray(approximate_spectrum, dtype=float)
    if masses.shape != (mu_overlap.shape[1],) or kernel.shape[0] != mu_overlap.shape[0]:
        raise ValueError("line fold-discretization axes are inconsistent")
    l1_absolute = 0.0
    yield_absolute = 0.0
    for source_mu, mass_value in enumerate(masses):
        mass = float(mass_value)
        if mass <= 0.0:
            continue
        targets = numpy.flatnonzero(mu_overlap[:, source_mu] > 0.0)
        if targets.size == 1:
            continue
        weights = mu_overlap[targets, source_mu]
        kernels = kernel[targets]
        average = numpy.sum(kernels * weights[:, None], axis=0)
        l1_absolute += mass * float(
            numpy.max(numpy.sum(numpy.abs(kernels - average), axis=1))
        )
        yield_absolute += mass * float(
            numpy.max(numpy.abs(numpy.sum(kernels, axis=1) - numpy.sum(average)))
        )
    approximate_total = float(numpy.sum(spectrum))
    lower_total = approximate_total - yield_absolute
    if approximate_total <= 0.0 or lower_total <= 0.0:
        relative_l1 = _FINITE_FAILURE
    else:
        relative_l1 = l1_absolute / lower_total
    return {
        "l1_absolute": l1_absolute,
        "yield_absolute": yield_absolute,
        "relative_l1_bound": relative_l1,
    }


def _robust_total_relative_change(
    reference_total: float,
    comparison_total: float,
    reference_yield_absolute_bound: float,
    comparison_yield_absolute_bound: float = 0.0,
) -> float:
    reference_low = reference_total - reference_yield_absolute_bound
    reference_high = reference_total + reference_yield_absolute_bound
    comparison_low = comparison_total - comparison_yield_absolute_bound
    comparison_high = comparison_total + comparison_yield_absolute_bound
    if reference_low <= 0.0 or comparison_low < 0.0 or not all(
        math.isfinite(value)
        for value in (reference_low, reference_high, comparison_low, comparison_high)
    ):
        return _FINITE_FAILURE
    return max(
        abs(comparison_low / reference_high - 1.0),
        abs(comparison_high / reference_low - 1.0),
    )


def _total_ratio_discretization_bound(
    reference_total: float,
    comparison_total: float,
    reference_yield_absolute_bound: float,
    comparison_yield_absolute_bound: float,
) -> float:
    reference_low = reference_total - reference_yield_absolute_bound
    reference_high = reference_total + reference_yield_absolute_bound
    comparison_low = comparison_total - comparison_yield_absolute_bound
    comparison_high = comparison_total + comparison_yield_absolute_bound
    if (
        reference_total <= 0.0
        or reference_low <= 0.0
        or comparison_low < 0.0
        or not all(
            math.isfinite(value)
            for value in (
                reference_high,
                comparison_total,
                comparison_low,
                comparison_high,
            )
        )
    ):
        return _FINITE_FAILURE
    nominal_ratio = comparison_total / reference_total
    return max(
        abs(comparison_low / reference_high - nominal_ratio),
        abs(comparison_high / reference_low - nominal_ratio),
    )


def _line_kernel(numpy: Any, response: DetectorResponse, energy_keV: float) -> Any:
    grid = response.line_energy_grid_keV
    scale = max(1.0, abs(float(energy_keV)))
    if energy_keV < float(grid[0]) - 1.0e-12 * scale or energy_keV > float(
        grid[-1]
    ) + 1.0e-12 * scale:
        raise ValueError(
            f"{response.detector_name} exact-line response grid does not cover {energy_keV} keV"
        )
    if energy_keV <= float(grid[0]):
        return response.line_response[0]
    if energy_keV >= float(grid[-1]):
        return response.line_response[-1]
    upper = int(numpy.searchsorted(grid, energy_keV, side="right"))
    lower = upper - 1
    fraction = (energy_keV - float(grid[lower])) / float(grid[upper] - grid[lower])
    return (1.0 - fraction) * response.line_response[lower] + fraction * response.line_response[
        upper
    ]


def _fold_continuum(numpy: Any, masses: Any, response: DetectorResponse) -> Any:
    values = numpy.asarray(masses, dtype=float)
    expected = response.continuum_response.shape[:2]
    if values.shape != expected:
        raise ValueError(f"continuum mass shape {values.shape} does not match response {expected}")
    return numpy.einsum("em,ems->s", values, response.continuum_response, optimize=True)


def _replica_uncertainty(
    numpy: Any, spectra: Any, truth_spectrum: Any
) -> tuple[float, float, float]:
    values = numpy.asarray(spectra, dtype=float)
    truth = numpy.asarray(truth_spectrum, dtype=float)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1:] != truth.shape:
        raise ValueError("detector-fold replicas have an invalid shape")
    truth_total = float(numpy.sum(truth))
    if not math.isfinite(truth_total) or truth_total <= 0.0:
        raise ValueError("detector-fold truth has no positive continuum response")
    yields = numpy.sum(values, axis=1)
    if numpy.any(~numpy.isfinite(values)) or numpy.any(values < 0.0):
        raise ValueError("detector-fold replicas must be finite and nonnegative")
    replicas = values.shape[0]
    yield_se = float(numpy.std(yields, ddof=1) / math.sqrt(replicas) / truth_total)
    truth_shape = truth / truth_total
    positive = yields > 0.0
    if int(numpy.count_nonzero(positive)) < 2:
        return max(yield_se, _FINITE_FAILURE), yield_se, _FINITE_FAILURE
    shapes = values[positive] / yields[positive, None]
    distances = 0.5 * numpy.sum(numpy.abs(shapes - truth_shape[None, :]), axis=1)
    shape_se = float(math.sqrt(float(numpy.mean(distances * distances)) / shapes.shape[0]))
    return max(yield_se, shape_se), yield_se, shape_se


def _normal_node_id(value: Any) -> str:
    result = str(value).lower().removeprefix("0x")
    if len(result) != 16 or any(character not in "0123456789abcdef" for character in result):
        raise ValueError(f"invalid final node id {value!r}")
    return result


def _candidate_bindings(
    payload: RuntimePayload,
    *,
    config_hash: str,
    archive_snapshot_hash: str,
    final_node_ids: Sequence[str],
) -> str:
    if payload.accepted_for_publication:
        raise ValueError("detector-fold validation requires an unpublished candidate payload")
    model_hash = payload.compute_model_content_sha256()
    if payload.model_content_sha256 != model_hash:
        raise ValueError("candidate model-content SHA-256 is inconsistent")
    if payload.metadata.get("config_sha256") != config_hash:
        raise ValueError("candidate and archive configuration hashes differ")
    validation = payload.metadata.get("validation_metadata")
    if not isinstance(validation, Mapping):
        raise ValueError("candidate lacks locked validation metadata")
    if validation.get("model_content_sha256") != model_hash:
        raise ValueError("candidate validation metadata names another model")
    if validation.get("archive_snapshot_sha256") != archive_snapshot_hash:
        raise ValueError("archive changed after the candidate was serialized")
    final = validation.get("final")
    if not isinstance(final, Mapping) or final.get("physics_passed") is not True:
        raise ValueError("candidate has not passed untouched-final physics validation")
    reported = [
        _normal_node_id(item.get("node_id", ""))
        for item in final.get("nodes", ())
        if isinstance(item, Mapping)
    ]
    if len(reported) != 64 or set(reported) != set(final_node_ids):
        raise ValueError("candidate is not bound to the campaign's exact fresh final set")
    return model_hash


def _archive_final_indices(root: Any, final_node_ids: Sequence[str]) -> list[int]:
    ids = [_normal_node_id(f"{int(value):016x}") for value in root["nodes/id"][:]]
    if len(ids) != len(set(ids)):
        raise ValueError("fit archive contains duplicate node IDs")
    by_id = {value: index for index, value in enumerate(ids)}
    missing = sorted(set(final_node_ids) - set(by_id))
    if missing:
        raise ValueError(f"fit archive is missing {len(missing)} current final nodes")
    indices = [by_id[value] for value in final_node_ids]
    roles = [str(root["nodes/role"].asstr()[index]) for index in indices]
    statuses = [str(root["nodes/status"].asstr()[index]) for index in indices]
    failures = [int(root["nodes/transport_failures"][index]) for index in indices]
    if any(role != "final" for role in roles):
        raise ValueError("current campaign final node is not tagged final in the archive")
    if any(status not in {"complete", "converged"} for status in statuses) or any(
        int(value) != 0 for value in failures
    ):
        raise ValueError("current final set contains incomplete or failed transport")
    return indices


def _require_archive_paths(root: Any) -> None:
    required = (
        "manifest",
        "axes/mu_edges",
        "axes/energy",
        "nodes/id",
        "nodes/role",
        "nodes/status",
        "nodes/transport_failures",
        "nodes/histories",
        "nodes/coordinates",
        "nodes/direct_fluence",
        "continuum/total/N",
        "continuum/total/S1",
        "continuum/total/replica_S1",
        "diagnostics/scattered_total_N",
        "diagnostics/scattered_total_S1",
        "replicas/histories",
        "lines/support_mask",
        "lines/energy_keV",
        "lines/total/N",
        "lines/total/S1",
        "lines/total/replica_S1",
        "blocks",
    )
    missing = [path for path in required if path not in root]
    if missing:
        raise ValueError(f"fit archive lacks detector-fold datasets: {missing}")
    if "total_energy_axis" not in root["axes"].attrs:
        raise ValueError("fit archive lacks its total continuum energy-axis identity")
    if "line_definitions_json" not in root["axes"].attrs:
        raise ValueError("fit archive lacks exact-line definitions")
    if "config_sha256" not in root["manifest"].attrs:
        raise ValueError("fit archive lacks its configuration hash")
    try:
        total_axis = int(root["axes"].attrs["total_energy_axis"])
    except (TypeError, ValueError) as error:
        raise ValueError("fit archive total continuum energy-axis identity is invalid") from error
    physical_axis = f"axes/energy/{total_axis}/physical_edges_keV"
    if physical_axis not in root:
        raise ValueError(f"fit archive lacks total continuum physical axes at {physical_axis}")


def _score_detector(
    root: Any,
    payload: RuntimePayload,
    response: DetectorResponse,
    final_indices: Sequence[int],
    final_node_ids: Sequence[str],
) -> dict[str, Any]:
    _, numpy = require_hdf5()
    archive_mu = numpy.asarray(root["axes/mu_edges"][:], dtype=float)
    mu_overlap = _overlap_matrix(numpy, archive_mu, response.mu_edges, "mu")
    n_source_mu = archive_mu.size - 1
    total_axis = int(root["axes"].attrs["total_energy_axis"])
    source_energy_all = root[f"axes/energy/{total_axis}/physical_edges_keV"]
    continuum_n = root["continuum/total/N"]
    continuum_s1 = root["continuum/total/S1"]
    continuum_replica_s1 = root["continuum/total/replica_S1"]
    replica_histories = root["replicas/histories"]
    raw_line_definitions = root["axes"].attrs["line_definitions_json"]
    if isinstance(raw_line_definitions, bytes):
        raw_line_definitions = raw_line_definitions.decode("utf-8")
    line_definitions = json.loads(str(raw_line_definitions))
    line_names = tuple(str(item["name"]) for item in line_definitions)
    if line_names != _LINE_NAMES:
        raise ValueError("fit archive exact-line inventory is not canonical")
    n_lines = len(line_names)
    line_support = root["lines/support_mask"]
    line_energy = root["lines/energy_keV"]
    line_n = root["lines/total/N"]
    line_s1 = root["lines/total/S1"]
    line_replica_s1 = root["lines/total/replica_S1"]
    scattered_n = root["diagnostics/scattered_total_N"]
    scattered_s1 = root["diagnostics/scattered_total_S1"]
    direct_fluence = root["nodes/direct_fluence"]
    coordinates = root["nodes/coordinates"]
    node_histories = root["nodes/histories"]

    nodes: list[dict[str, Any]] = []
    for index, node_id in zip(final_indices, final_node_ids):
        source_energy = numpy.asarray(source_energy_all[index], dtype=float)
        energy_overlap = _overlap_matrix(
            numpy, source_energy, response.continuum_energy_edges_keV, "continuum energy"
        )
        n_source_energy = source_energy.size - 1
        source_shape = (n_source_energy, n_source_mu)
        counts = numpy.asarray(continuum_n[index], dtype=float).reshape(source_shape)
        first = numpy.asarray(continuum_s1[index], dtype=float).reshape(source_shape)
        if (
            numpy.any(~numpy.isfinite(first))
            or numpy.any(first < 0.0)
            or numpy.any(counts <= 0.0)
            or not numpy.all(counts == float(node_histories[index]))
        ):
            raise ValueError(f"final node {node_id} has invalid continuum sufficient statistics")
        truth_source = first / counts
        truth = _rebin_2d(numpy, truth_source, energy_overlap, mu_overlap)
        truth_spectrum = _fold_continuum(numpy, truth, response)
        truth_total = float(numpy.sum(truth_spectrum))
        if not math.isfinite(truth_total) or truth_total <= 0.0:
            raise ValueError(f"{response.detector_name} has no continuum sensitivity at {node_id}")
        continuum_discretization = _fold_discretization_bounds(
            numpy,
            truth_source,
            energy_overlap,
            mu_overlap,
            response.continuum_response,
            truth_spectrum,
        )

        histories = numpy.asarray(replica_histories[index], dtype=float)
        if histories.shape != (32,) or numpy.any(histories <= 0.0) or not math.isclose(
            float(numpy.sum(histories)), float(node_histories[index]), rel_tol=0.0, abs_tol=0.0
        ):
            raise ValueError(f"final node {node_id} has invalid replica history accounting")
        replica_source = numpy.asarray(continuum_replica_s1[index], dtype=float).reshape(
            32, n_source_energy, n_source_mu
        )
        if numpy.any(~numpy.isfinite(replica_source)) or numpy.any(replica_source < 0.0):
            raise ValueError(f"final node {node_id} has invalid continuum replicas")
        replica_source = replica_source / histories[:, None, None]
        replica_rebinned = _rebin_replicas(
            numpy, replica_source, energy_overlap, mu_overlap
        )
        replica_spectra = numpy.einsum(
            "rem,ems->rs", replica_rebinned, response.continuum_response, optimize=True
        )
        (
            continuum_replica_se,
            continuum_replica_yield_se,
            continuum_replica_shape_se,
        ) = _replica_uncertainty(numpy, replica_spectra, truth_spectrum)

        energy_keV, distance_cm, source_height_cm, detector_height_cm, _ = map(
            float, coordinates[index]
        )
        direct = float(direct_fluence[index])
        if not all(
            math.isfinite(value)
            for value in (
                energy_keV,
                distance_cm,
                source_height_cm,
                detector_height_cm,
                direct,
            )
        ) or direct <= 0.0:
            raise ValueError(f"final node {node_id} has invalid physical coordinates")
        direct_detector_spectrum = direct * _line_kernel(
            numpy, response, energy_keV
        )[-1]
        detector_total_reference = truth_total + float(
            numpy.sum(direct_detector_spectrum)
        )
        incident_count = int(scattered_n[index])
        incident_total = (
            float(scattered_s1[index]) / incident_count if incident_count else 0.0
        )
        if detector_total_reference <= 0.0 or incident_total <= 0.0:
            raise ValueError(f"final node {node_id} has invalid total line normalizers")
        evaluated = payload.evaluate(
            energy_keV,
            distance_cm * 0.01,
            source_height_cm * 0.01,
            detector_height_cm * 0.01,
            response.continuum_energy_edges_keV,
            response.mu_edges,
            direct_fluence_cm2=direct,
        )
        predicted = numpy.asarray(evaluated["continuum"], dtype=float)
        if predicted.shape != response.continuum_response.shape[:2]:
            raise ValueError(f"candidate continuum dimensions differ at final node {node_id}")
        if numpy.any(~numpy.isfinite(predicted)) or numpy.any(predicted < 0.0):
            raise ValueError(
                f"candidate continuum is nonfinite or negative at final node {node_id}"
            )
        predicted_spectrum = _fold_continuum(numpy, predicted, response)
        predicted_total = float(numpy.sum(predicted_spectrum))
        continuum_error = _robust_total_relative_change(
            truth_total,
            predicted_total,
            continuum_discretization["yield_absolute"],
        )
        observed_spectral_tv = (
            0.5
            * float(
                numpy.sum(
                    numpy.abs(
                        predicted_spectrum / predicted_total - truth_spectrum / truth_total
                    )
                )
            )
            if predicted_total > 0.0
            else 1.0
        )
        spectral_tv = min(
            1.0,
            observed_spectral_tv + continuum_discretization["spectral_tv_bound"],
        )

        predicted_lines = {str(item["name"]): item for item in evaluated["lines"]}
        if len(predicted_lines) != len(evaluated["lines"]):
            raise ValueError(f"candidate repeats an exact line at final node {node_id}")
        support = numpy.asarray(line_support[index], dtype=bool)
        energies = numpy.asarray(line_energy[index], dtype=float)
        line_counts = numpy.asarray(line_n[index], dtype=float).reshape(n_lines, n_source_mu)
        line_first = numpy.asarray(line_s1[index], dtype=float).reshape(n_lines, n_source_mu)
        line_replicas = numpy.asarray(line_replica_s1[index], dtype=float).reshape(
            32, n_lines, n_source_mu
        )
        if (
            numpy.any(~numpy.isfinite(energies))
            or numpy.any(~numpy.isfinite(line_first))
            or numpy.any(line_first < 0.0)
            or numpy.any(~numpy.isfinite(line_replicas))
            or numpy.any(line_replicas < 0.0)
            or numpy.any(line_counts <= 0.0)
            or not numpy.all(line_counts == float(node_histories[index]))
        ):
            raise ValueError(f"final node {node_id} has invalid exact-line statistics")
        line_errors: list[float] = []
        line_discretization_bounds: list[float] = []
        line_replica_metrics: list[tuple[float, float, float]] = []
        line_policy_reports: list[dict[str, Any]] = []
        for line_index, line_name in enumerate(line_names):
            record = predicted_lines.get(line_name)
            if not bool(support[line_index]):
                if record is not None and abs(float(sum(record["mu_bin_masses"]))) > 1.0e-15:
                    raise ValueError(f"candidate leaks unsupported line {line_name} at {node_id}")
                continue
            if record is None:
                raise ValueError(f"candidate omits supported line {line_name} at {node_id}")
            exact_energy = float(energies[line_index])
            if not math.isclose(
                float(record["energy_keV"]), exact_energy, rel_tol=0.0, abs_tol=1.0e-9
            ):
                raise ValueError(f"candidate changes exact line energy {line_name} at {node_id}")
            kernel = _line_kernel(numpy, response, exact_energy)
            truth_mu = mu_overlap @ (line_first[line_index] / line_counts[line_index])
            prediction_mu = numpy.asarray(record["mu_bin_masses"], dtype=float)
            if prediction_mu.shape != (response.mu_edges.size - 1,):
                raise ValueError(f"candidate line angular dimensions differ at {node_id}")
            if numpy.any(~numpy.isfinite(prediction_mu)) or numpy.any(prediction_mu < 0.0):
                raise ValueError(f"candidate exact line is nonfinite or negative at {node_id}")
            truth_line_spectrum = truth_mu @ kernel
            predicted_line_spectrum = prediction_mu @ kernel
            standalone_truth_line_spectrum = truth_line_spectrum
            standalone_predicted_line_spectrum = predicted_line_spectrum
            folded_truth_line = float(numpy.sum(standalone_truth_line_spectrum))
            if folded_truth_line <= 0.0:
                raise ValueError(
                    f"{response.detector_name} has no sensitivity to supported "
                    f"{line_name} at {node_id}"
                )
            line_discretization = _line_discretization_bounds(
                numpy,
                line_first[line_index] / line_counts[line_index],
                mu_overlap,
                kernel,
                truth_line_spectrum,
            )
            line_discretization_bounds.append(
                line_discretization["relative_l1_bound"]
            )
            observed_difference = float(
                numpy.sum(
                    numpy.abs(
                        standalone_predicted_line_spectrum
                        - standalone_truth_line_spectrum
                    )
                )
            )
            replica_mu = numpy.einsum(
                "mj,rj->rm",
                mu_overlap,
                line_replicas[:, line_index, :] / histories[:, None],
                optimize=True,
            )
            replica_line_spectra = replica_mu @ kernel
            incident_truth_mu = line_first[line_index] / line_counts[line_index]
            incident_line_yield = float(numpy.sum(incident_truth_mu))
            replica_incident_yields = numpy.sum(
                line_replicas[:, line_index, :] / histories[:, None], axis=1
            )
            incident_line_yield_se = float(
                numpy.std(replica_incident_yields, ddof=1) / math.sqrt(32)
            )
            material, upper_fraction = exact_line_materiality(
                incident_line_yield, incident_line_yield_se, incident_total
            )
            if line_name == "elastic_E0":
                truth_line_spectrum = (
                    standalone_truth_line_spectrum + direct_detector_spectrum
                )
                predicted_line_spectrum = (
                    standalone_predicted_line_spectrum + direct_detector_spectrum
                )
                replica_line_spectra = (
                    replica_line_spectra + direct_detector_spectrum[None, :]
                )
                denominator = float(numpy.sum(truth_line_spectrum)) - line_discretization[
                    "yield_absolute"
                ]
                effective_error = (
                    (observed_difference + line_discretization["l1_absolute"])
                    / denominator
                    if denominator > 0.0
                    else _FINITE_FAILURE
                )
                effective_discretization = (
                    line_discretization["l1_absolute"] / denominator
                    if denominator > 0.0
                    else _FINITE_FAILURE
                )
                effective_replica = _replica_uncertainty(
                    numpy, replica_line_spectra, truth_line_spectrum
                )
                material = True
                policy_kind = "combined-direct-plus-elastic-relative"
            elif material:
                denominator = folded_truth_line - line_discretization["yield_absolute"]
                effective_error = (
                    (observed_difference + line_discretization["l1_absolute"])
                    / denominator
                    if denominator > 0.0
                    else _FINITE_FAILURE
                )
                effective_discretization = (
                    line_discretization["relative_l1_bound"]
                )
                effective_replica = _replica_uncertainty(
                    numpy, replica_line_spectra, truth_line_spectrum
                )
                policy_kind = "material-line-relative"
            else:
                absolute_error_fraction = (
                    observed_difference + line_discretization["l1_absolute"]
                ) / detector_total_reference
                absolute_discretization_fraction = (
                    line_discretization["l1_absolute"] / detector_total_reference
                )
                raw_replica = _replica_uncertainty(
                    numpy, replica_line_spectra, truth_line_spectrum
                )
                folded_fraction = folded_truth_line / detector_total_reference
                absolute_yield_noise = raw_replica[1] * folded_fraction
                absolute_shape_noise = raw_replica[2] * folded_fraction
                # Express absolute-policy quantities on the relative scalar
                # scales, so the shared <=1% physics and <=0.25% noise gates
                # apply to every line policy unchanged.
                error_scale = (
                    EXACT_LINE_RELATIVE_ERROR_FRACTION
                    / EXACT_LINE_ABSOLUTE_ERROR_FRACTION
                )
                noise_scale = (
                    EXACT_LINE_RELATIVE_NOISE_FRACTION
                    / EXACT_LINE_ABSOLUTE_NOISE_FRACTION
                )
                effective_error = absolute_error_fraction * error_scale
                effective_discretization = (
                    absolute_discretization_fraction * error_scale
                )
                effective_replica = (
                    max(absolute_yield_noise, absolute_shape_noise) * noise_scale,
                    absolute_yield_noise * noise_scale,
                    absolute_shape_noise * noise_scale,
                )
                policy_kind = "weak-line-absolute"
            line_errors.append(effective_error)
            line_discretization_bounds[-1] = effective_discretization
            line_replica_metrics.append(effective_replica)
            line_policy_reports.append(
                {
                    "name": line_name,
                    "policy": policy_kind,
                    "materiality_upper_fraction_3sigma": upper_fraction,
                    "incident_truth_yield": incident_line_yield,
                    "incident_yield_standard_error": incident_line_yield_se,
                    "standalone_folded_truth_yield": folded_truth_line,
                    "effective_physics_error": effective_error,
                    "effective_replica_relative_se": effective_replica[0],
                    "effective_discretization_bound": effective_discretization,
                }
            )
        unexpected = set(predicted_lines) - {
            name for name, active in zip(line_names, support) if bool(active)
        }
        if unexpected:
            raise ValueError(f"candidate returns unsupported exact lines {sorted(unexpected)}")
        if not line_errors:
            raise ValueError(f"final node {node_id} has no supported exact line")

        replica_se = max(
            [continuum_replica_se, *(value[0] for value in line_replica_metrics)]
        )
        replica_yield_se = max(
            [continuum_replica_yield_se, *(value[1] for value in line_replica_metrics)]
        )
        replica_shape_se = max(
            [continuum_replica_shape_se, *(value[2] for value in line_replica_metrics)]
        )

        metric_values = (
            continuum_error,
            spectral_tv,
            max(line_errors),
            replica_se,
            replica_yield_se,
            replica_shape_se,
            continuum_discretization["total_relative_bound"],
            continuum_discretization["spectral_tv_bound"],
            max(line_discretization_bounds, default=_FINITE_FAILURE),
        )
        if any(not math.isfinite(value) or value < 0.0 for value in metric_values):
            raise ValueError(f"detector-fold metrics are invalid at final node {node_id}")
        nodes.append(
            {
                "node_id": node_id,
                "continuum_total_relative_error": continuum_error,
                "spectral_tv": spectral_tv,
                "exact_line_angular_fold_relative_error": max(line_errors),
                # This conservative scalar enforces both detector-folded yield
                # and normalized spectral-shape replica targets.
                "replica_relative_se": replica_se,
                "replica_yield_relative_se": replica_yield_se,
                "replica_spectral_tv_se": replica_shape_se,
                "continuum_fold_discretization_total_relative_bound": continuum_discretization[
                    "total_relative_bound"
                ],
                "continuum_fold_discretization_spectral_tv_bound": continuum_discretization[
                    "spectral_tv_bound"
                ],
                "exact_line_fold_discretization_relative_bound": max(
                    line_discretization_bounds, default=_FINITE_FAILURE
                ),
                "exact_line_policy": line_policy_reports,
            }
        )

    maxima = {
        name: max(float(item[name]) for item in nodes)
        for name in (
            "continuum_total_relative_error",
            "spectral_tv",
            "exact_line_angular_fold_relative_error",
            "replica_relative_se",
        )
    }
    passed = bool(
        maxima["continuum_total_relative_error"] <= 0.01
        and maxima["spectral_tv"] <= 0.01
        and maxima["exact_line_angular_fold_relative_error"] <= 0.01
        and maxima["replica_relative_se"] <= 0.0025
        and max(
            float(item["continuum_fold_discretization_total_relative_bound"])
            for item in nodes
        )
        <= 0.01
        and max(
            float(item["continuum_fold_discretization_spectral_tv_bound"])
            for item in nodes
        )
        <= 0.01
        and max(
            float(item["exact_line_fold_discretization_relative_bound"])
            for item in nodes
        )
        <= 0.01
    )
    geometry_json = json.dumps(
        response.geometry, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    provenance_json = json.dumps(
        response.provenance, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return {
        "response_sha256": response.sha256,
        "response_schema": RESPONSE_SCHEMA,
        "response_semantics": RESPONSE_SEMANTICS,
        "line_interpolation": LINE_INTERPOLATION,
        "exact_line_acceptance_policy": "contribution-aware-v1",
        "geometry": response.geometry,
        "geometry_id": response.geometry["geometry_id"],
        "geometry_sha256": hashlib.sha256(geometry_json.encode("utf-8")).hexdigest(),
        "provenance": response.provenance,
        "provenance_sha256": hashlib.sha256(
            provenance_json.encode("utf-8")
        ).hexdigest(),
        "generator": {
            "name": response.generator_name,
            "version": response.generator_version,
            "sha256": response.generator_sha256,
        },
        "fold_discretization": {
            "method": "archive-cell-response-kernel-diameter-v1",
            "gate": 0.01,
            "max_continuum_total_relative_bound": max(
                float(item["continuum_fold_discretization_total_relative_bound"])
                for item in nodes
            ),
            "max_continuum_spectral_tv_bound": max(
                float(item["continuum_fold_discretization_spectral_tv_bound"])
                for item in nodes
            ),
            "max_exact_line_relative_bound": max(
                float(item["exact_line_fold_discretization_relative_bound"])
                for item in nodes
            ),
        },
        "dimensions": {
            "continuum_energy_bins": int(response.continuum_energy_edges_keV.size - 1),
            "mu_bins": int(response.mu_edges.size - 1),
            "spectrum_bins": int(response.spectrum_edges_keV.size - 1),
            "line_energy_grid_points": int(response.line_energy_grid_keV.size),
        },
        "nodes": nodes,
        "maxima": maxima,
        "passed": passed,
    }


def _eligible_mu_pilot_indices(root: Any) -> tuple[list[int], list[str]]:
    roles = list(root["nodes/role"].asstr()[:])
    statuses = list(root["nodes/status"].asstr()[:])
    failures = [int(value) for value in root["nodes/transport_failures"][:]]
    complete = [
        status in {"complete", "converged"} and failure == 0
        for status, failure in zip(statuses, failures)
    ]
    indices = [
        index
        for index, (role, eligible) in enumerate(zip(roles, complete))
        if role == "pilot" and eligible
    ]
    if not indices:
        indices = [
            index
            for index, (role, eligible) in enumerate(zip(roles, complete))
            if role == "train" and eligible
        ]
    if not indices:
        raise ValueError("mu-resolution validation has no transport-complete pilot/train nodes")
    node_ids = [f"{int(root['nodes/id'][index]):016x}" for index in indices]
    if len(node_ids) != len(set(node_ids)):
        raise ValueError("mu-resolution pilot inventory repeats a node ID")
    return indices, node_ids


def _score_mu_resolution_detector(
    root: Any,
    response: DetectorResponse,
    pilot_indices: Sequence[int],
    pilot_node_ids: Sequence[str],
    *,
    tolerance: float,
) -> dict[str, Any]:
    _, numpy = require_hdf5()
    fine_mu_edges = numpy.asarray(root["axes/mu_edges"][:], dtype=float)
    if fine_mu_edges.shape != (65,):
        raise ValueError("mu-resolution detector validation requires exactly 64 source bins")
    coarse_mu_edges = fine_mu_edges[::2]
    if coarse_mu_edges.shape != (33,) or coarse_mu_edges[-1] != fine_mu_edges[-1]:
        raise ValueError("64-bin archive mu edges cannot be paired conservatively into 32 bins")
    fine_mu_overlap = _overlap_matrix(
        numpy, fine_mu_edges, response.mu_edges, "fine mu"
    )
    coarse_mu_overlap = _overlap_matrix(
        numpy, coarse_mu_edges, response.mu_edges, "coarse mu"
    )
    total_axis = int(root["axes"].attrs["total_energy_axis"])
    source_energy_all = root[f"axes/energy/{total_axis}/physical_edges_keV"]
    continuum_n = root["continuum/total/N"]
    continuum_s1 = root["continuum/total/S1"]
    node_histories = root["nodes/histories"]
    raw_line_definitions = root["axes"].attrs["line_definitions_json"]
    if isinstance(raw_line_definitions, bytes):
        raw_line_definitions = raw_line_definitions.decode("utf-8")
    line_definitions = json.loads(str(raw_line_definitions))
    line_names = tuple(str(item["name"]) for item in line_definitions)
    if line_names != _LINE_NAMES:
        raise ValueError("fit archive exact-line inventory is not canonical")
    line_support = root["lines/support_mask"]
    line_energy = root["lines/energy_keV"]
    line_n = root["lines/total/N"]
    line_s1 = root["lines/total/S1"]

    nodes: list[dict[str, Any]] = []
    for index, node_id in zip(pilot_indices, pilot_node_ids):
        source_energy = numpy.asarray(source_energy_all[index], dtype=float)
        energy_overlap = _overlap_matrix(
            numpy, source_energy, response.continuum_energy_edges_keV, "continuum energy"
        )
        n_energy = source_energy.size - 1
        counts = numpy.asarray(continuum_n[index], dtype=float).reshape(n_energy, 64)
        first = numpy.asarray(continuum_s1[index], dtype=float).reshape(n_energy, 64)
        if (
            numpy.any(~numpy.isfinite(first))
            or numpy.any(first < 0.0)
            or numpy.any(counts <= 0.0)
            or not numpy.all(counts == float(node_histories[index]))
        ):
            raise ValueError(f"pilot node {node_id} has invalid continuum sufficient statistics")
        fine_source = first / counts
        coarse_source = fine_source.reshape(n_energy, 32, 2).sum(axis=2)
        fine_rebinned = _rebin_2d(
            numpy, fine_source, energy_overlap, fine_mu_overlap
        )
        coarse_rebinned = _rebin_2d(
            numpy, coarse_source, energy_overlap, coarse_mu_overlap
        )
        fine_spectrum = _fold_continuum(numpy, fine_rebinned, response)
        coarse_spectrum = _fold_continuum(numpy, coarse_rebinned, response)
        fine_total = float(numpy.sum(fine_spectrum))
        coarse_total = float(numpy.sum(coarse_spectrum))
        if fine_total <= 0.0 or not all(
            math.isfinite(value) for value in (fine_total, coarse_total)
        ):
            raise ValueError(f"{response.detector_name} has no continuum pilot sensitivity")
        fine_discretization = _fold_discretization_bounds(
            numpy,
            fine_source,
            energy_overlap,
            fine_mu_overlap,
            response.continuum_response,
            fine_spectrum,
        )
        coarse_discretization = _fold_discretization_bounds(
            numpy,
            coarse_source,
            energy_overlap,
            coarse_mu_overlap,
            response.continuum_response,
            coarse_spectrum,
        )
        total_change = _robust_total_relative_change(
            fine_total,
            coarse_total,
            fine_discretization["yield_absolute"],
            coarse_discretization["yield_absolute"],
        )
        total_discretization_bound = _total_ratio_discretization_bound(
            fine_total,
            coarse_total,
            fine_discretization["yield_absolute"],
            coarse_discretization["yield_absolute"],
        )
        observed_tv = (
            0.5
            * float(
                numpy.sum(
                    numpy.abs(
                        coarse_spectrum / coarse_total - fine_spectrum / fine_total
                    )
                )
            )
            if coarse_total > 0.0
            else 1.0
        )
        spectral_discretization_bound = min(
            1.0,
            fine_discretization["spectral_tv_bound"]
            + coarse_discretization["spectral_tv_bound"],
        )
        spectral_tv = min(1.0, observed_tv + spectral_discretization_bound)

        support = numpy.asarray(line_support[index], dtype=bool)
        energies = numpy.asarray(line_energy[index], dtype=float)
        line_counts = numpy.asarray(line_n[index], dtype=float).reshape(2, 64)
        line_first = numpy.asarray(line_s1[index], dtype=float).reshape(2, 64)
        if (
            numpy.any(~numpy.isfinite(energies))
            or numpy.any(~numpy.isfinite(line_first))
            or numpy.any(line_first < 0.0)
            or numpy.any(line_counts <= 0.0)
            or not numpy.all(line_counts == float(node_histories[index]))
        ):
            raise ValueError(f"pilot node {node_id} has invalid exact-line statistics")
        line_changes: list[float] = []
        line_bounds: list[float] = []
        for line_index, line_name in enumerate(line_names):
            if not bool(support[line_index]):
                continue
            exact_energy = float(energies[line_index])
            kernel = _line_kernel(numpy, response, exact_energy)
            fine_line_source = line_first[line_index] / line_counts[line_index]
            coarse_line_source = fine_line_source.reshape(32, 2).sum(axis=1)
            fine_line_mu = fine_mu_overlap @ fine_line_source
            coarse_line_mu = coarse_mu_overlap @ coarse_line_source
            fine_line_spectrum = fine_line_mu @ kernel
            coarse_line_spectrum = coarse_line_mu @ kernel
            fine_line_total = float(numpy.sum(fine_line_spectrum))
            if fine_line_total <= 0.0:
                raise ValueError(
                    f"{response.detector_name} has no sensitivity to supported "
                    f"{line_name} at pilot {node_id}"
                )
            fine_line_discretization = _line_discretization_bounds(
                numpy,
                fine_line_source,
                fine_mu_overlap,
                kernel,
                fine_line_spectrum,
            )
            coarse_line_discretization = _line_discretization_bounds(
                numpy,
                coarse_line_source,
                coarse_mu_overlap,
                kernel,
                coarse_line_spectrum,
            )
            denominator = (
                fine_line_total - fine_line_discretization["yield_absolute"]
            )
            observed_l1 = float(
                numpy.sum(numpy.abs(coarse_line_spectrum - fine_line_spectrum))
            )
            line_discretization_l1 = (
                fine_line_discretization["l1_absolute"]
                + coarse_line_discretization["l1_absolute"]
            )
            line_changes.append(
                (observed_l1 + line_discretization_l1) / denominator
                if denominator > 0.0
                else _FINITE_FAILURE
            )
            line_bounds.append(
                line_discretization_l1 / denominator
                if denominator > 0.0
                else _FINITE_FAILURE
            )
        if not line_changes:
            raise ValueError(f"pilot node {node_id} has no supported exact line")

        metric_values = (
            total_change,
            spectral_tv,
            max(line_changes),
            total_discretization_bound,
            spectral_discretization_bound,
            max(line_bounds),
        )
        if any(not math.isfinite(value) or value < 0.0 for value in metric_values):
            raise ValueError(f"mu-resolution detector metrics are invalid at pilot {node_id}")
        nodes.append(
            {
                "node_id": node_id,
                "detector_folded_total_relative_change": total_change,
                "detector_folded_spectral_tv": spectral_tv,
                "exact_line_angular_fold_relative_change": max(line_changes),
                "continuum_fold_discretization_total_relative_bound": (
                    total_discretization_bound
                ),
                "continuum_fold_discretization_spectral_tv_bound": (
                    spectral_discretization_bound
                ),
                "exact_line_fold_discretization_relative_bound": max(line_bounds),
            }
        )

    metric_names = (
        "detector_folded_total_relative_change",
        "detector_folded_spectral_tv",
        "exact_line_angular_fold_relative_change",
    )
    maxima = {
        f"max_{name}": max(float(item[name]) for item in nodes)
        for name in metric_names
    }
    discretization_maxima = {
        "max_continuum_total_relative_bound": max(
            float(item["continuum_fold_discretization_total_relative_bound"])
            for item in nodes
        ),
        "max_continuum_spectral_tv_bound": max(
            float(item["continuum_fold_discretization_spectral_tv_bound"])
            for item in nodes
        ),
        "max_exact_line_relative_bound": max(
            float(item["exact_line_fold_discretization_relative_bound"])
            for item in nodes
        ),
    }
    passed = all(
        value <= tolerance
        for value in (*maxima.values(), *discretization_maxima.values())
    )
    geometry_json = json.dumps(
        response.geometry, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    provenance_json = json.dumps(
        response.provenance, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return {
        "response_sha256": response.sha256,
        "response_schema": RESPONSE_SCHEMA,
        "response_semantics": RESPONSE_SEMANTICS,
        "line_interpolation": LINE_INTERPOLATION,
        "geometry": response.geometry,
        "geometry_id": response.geometry["geometry_id"],
        "geometry_sha256": hashlib.sha256(geometry_json.encode("utf-8")).hexdigest(),
        "provenance": response.provenance,
        "provenance_sha256": hashlib.sha256(
            provenance_json.encode("utf-8")
        ).hexdigest(),
        "generator": {
            "name": response.generator_name,
            "version": response.generator_version,
            "sha256": response.generator_sha256,
        },
        "dimensions": {
            "continuum_energy_bins": int(response.continuum_energy_edges_keV.size - 1),
            "mu_bins": int(response.mu_edges.size - 1),
            "spectrum_bins": int(response.spectrum_edges_keV.size - 1),
            "line_energy_grid_points": int(response.line_energy_grid_keV.size),
        },
        "fold_discretization": {
            "method": "archive-cell-response-kernel-diameter-v1",
            "gate": tolerance,
            **discretization_maxima,
        },
        "nodes": nodes,
        "maxima": maxima,
        "passed": bool(passed),
    }


def generate_mu_resolution_validation_report(
    archive_path: str | Path,
    hpge_response_path: str | Path,
    nai_response_path: str | Path,
    output_path: str | Path,
    *,
    tolerance: float = 0.002,
) -> dict[str, Any]:
    """Generate the bound HPGe/NaI 64-to-32 detector-fold comparison."""

    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("mu-resolution tolerance must be finite and positive")
    destination = Path(output_path)
    inputs = tuple(
        Path(value) for value in (archive_path, hpge_response_path, nai_response_path)
    )
    if any(destination.resolve() == value.resolve() for value in inputs):
        raise ValueError("mu-resolution report output must be distinct from every input")
    responses = {
        "HPGe": load_detector_response(
            hpge_response_path, expected_detector_name="HPGe"
        ),
        "NaI": load_detector_response(nai_response_path, expected_detector_name="NaI"),
    }
    h5py, _ = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        if str(root.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
            raise ValueError("unsupported fit archive schema")
        _require_archive_paths(root)
        if len(root["axes/mu_edges"]) - 1 != 64:
            raise ValueError("mu-resolution detector validation requires a 64-bin archive")
        snapshot_hash = archive_snapshot_sha256(root)
        config_hash = str(root["manifest"].attrs["config_sha256"])
        pilot_indices, pilot_node_ids = _eligible_mu_pilot_indices(root)
        detectors = {
            name: _score_mu_resolution_detector(
                root,
                response,
                pilot_indices,
                pilot_node_ids,
                tolerance=tolerance,
            )
            for name, response in responses.items()
        }
    response_paths = {"HPGe": Path(hpge_response_path), "NaI": Path(nai_response_path)}
    for name, response_path in response_paths.items():
        detectors[name]["response_path"] = Path(
            os.path.relpath(response_path.resolve(), destination.parent.resolve())
        ).as_posix()
    report = {
        "schema": "CEEFIT-MU-RESOLUTION-VALIDATION-v1",
        "config_sha256": config_hash,
        "archive_snapshot_sha256": snapshot_hash,
        "source_bins": 64,
        "candidate_bins": 32,
        "tolerance": tolerance,
        "pilot_node_ids": pilot_node_ids,
        "detectors": detectors,
        "passed": all(value["passed"] for value in detectors.values()),
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)
    return report


def generate_detector_fold_report(
    candidate_payload: str | Path,
    archive_path: str | Path,
    campaign_manifest_path: str | Path,
    hpge_response_path: str | Path,
    nai_response_path: str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    """Fold the exact current final set and atomically write publication evidence."""

    destination = Path(output_path)
    inputs = tuple(
        Path(value)
        for value in (
            candidate_payload,
            archive_path,
            campaign_manifest_path,
            hpge_response_path,
            nai_response_path,
        )
    )
    if any(destination.resolve() == value.resolve() for value in inputs):
        raise ValueError("detector-fold report output must be distinct from every input artifact")
    campaign = load_manifest(campaign_manifest_path)
    try:
        raw_final_nodes = campaign["design"]["nodes"]["final"]
        final_node_ids = [
            _normal_node_id(item["node_id_hex"])
            for item in raw_final_nodes
        ]
    except (KeyError, TypeError) as error:
        raise ValueError("campaign manifest lacks its current final-node inventory") from error
    if len(final_node_ids) != 64 or len(set(final_node_ids)) != 64:
        raise ValueError("detector-fold validation requires exactly 64 fresh final nodes")
    candidate_path = Path(candidate_payload)
    candidate_bytes = candidate_path.read_bytes()
    candidate_sha256 = hashlib.sha256(candidate_bytes).hexdigest()
    payload = RuntimePayload.load_binary(candidate_path, allow_candidate=True)
    if hashlib.sha256(candidate_path.read_bytes()).hexdigest() != candidate_sha256:
        raise ValueError("candidate payload changed while detector folds were being prepared")
    responses = {
        "HPGe": load_detector_response(
            hpge_response_path, expected_detector_name="HPGe"
        ),
        "NaI": load_detector_response(nai_response_path, expected_detector_name="NaI"),
    }
    h5py, _ = require_hdf5()
    with h5py.File(archive_path, "r") as root:
        if str(root.attrs.get("schema", "")) != "CEEFIT-HDF5-v1":
            raise ValueError("unsupported fit archive schema")
        _require_archive_paths(root)
        config_hash = str(root["manifest"].attrs["config_sha256"])
        if config_hash != campaign["config"]["sha256"]:
            raise ValueError("campaign and fit archive configuration hashes differ")
        snapshot_hash = archive_snapshot_sha256(root)
        model_hash = _candidate_bindings(
            payload,
            config_hash=config_hash,
            archive_snapshot_hash=snapshot_hash,
            final_node_ids=final_node_ids,
        )
        final_indices = _archive_final_indices(root, final_node_ids)
        detectors = {
            name: _score_detector(root, payload, response, final_indices, final_node_ids)
            for name, response in responses.items()
        }
    response_paths = {"HPGe": Path(hpge_response_path), "NaI": Path(nai_response_path)}
    for name, response_path in response_paths.items():
        detectors[name]["response_path"] = Path(
            os.path.relpath(response_path.resolve(), destination.parent.resolve())
        ).as_posix()

    report = {
        "schema": REPORT_SCHEMA,
        "config_sha256": config_hash,
        "archive_snapshot_sha256": snapshot_hash,
        "model_content_sha256": model_hash,
        "candidate_payload_sha256": candidate_sha256,
        "final_node_ids": final_node_ids,
        "thresholds": {
            "continuum_total_relative_error": 0.01,
            "spectral_tv": 0.01,
            "exact_line_angular_fold_relative_error": 0.01,
            "replica_relative_se": 0.0025,
            "exact_line_policy": "contribution-aware-v1",
            "weak_line_absolute_error_fraction": EXACT_LINE_ABSOLUTE_ERROR_FRACTION,
            "weak_line_absolute_noise_fraction": EXACT_LINE_ABSOLUTE_NOISE_FRACTION,
        },
        "detectors": detectors,
        "passed": all(value["passed"] for value in detectors.values()),
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)
    return report
