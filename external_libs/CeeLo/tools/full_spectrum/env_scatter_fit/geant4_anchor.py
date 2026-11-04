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

"""Build GEANT4-anchor evidence from raw reference arms and archive moments.

The publication gate must not trust a hand-written ``passed`` flag.  This
module hashes and merges the raw GEANT4 sphere-flux CSV arms, obtains the CeeLo
mean and standard error directly from the identity-bound fit archive, and
computes every comparison in the emitted report.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence
import xml.etree.ElementTree as ElementTree

from .archive import require_hdf5
from .publication import archive_snapshot_sha256, file_sha256


MANIFEST_SCHEMA = "CEEFIT-GEANT4-ANCHOR-MANIFEST-v1"
REPORT_SCHEMA = "CEEFIT-GEANT4-ANCHOR-AGREEMENT-v1"
REQUIRED_PHYSICS_LIST = "G4EmLivermorePhysics"
REQUIRED_OBSERVABLE = "scattered-track-length-flux-per-primary-v1"
REQUIRED_PRIMARY_SAMPLING = "gps-isotropic-4pi"
DEFAULT_COMPARISON_SIGMA = 2.0
MINIMUM_FIELD_WORLD_HALF_EXTENT_M = 60.0
REFERENCE_ENERGY_BINS = 48
REFERENCE_MU_BINS = 32


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
        raise ValueError(f"invalid GEANT4 anchor manifest: {error}") from error
    if not isinstance(value, Mapping):
        raise ValueError("GEANT4 anchor manifest must be a JSON object")
    return value


def _atomic_json(path: str | Path, value: Mapping[str, Any]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = (
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode("ascii")
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


def _sha256(value: Any, name: str) -> str:
    result = str(value)
    if not (
        len(result) == 64
        and all(character in "0123456789abcdef" for character in result)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256")
    return result


def _node_hex(value: Any) -> str:
    result = str(value).lower().removeprefix("0x")
    if not (
        len(result) == 16
        and all(character in "0123456789abcdef" for character in result)
    ):
        raise ValueError("GEANT4 anchor node_id must contain exactly 16 hex digits")
    return result


def _finite_nonnegative(value: Any, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return result


def _length_cm(value: str, unit: str, name: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{name} has an invalid length") from error
    scales = {"mm": 0.1, "cm": 1.0, "m": 100.0}
    if unit not in scales or not math.isfinite(number):
        raise ValueError(f"{name} has an unsupported length unit")
    return number * scales[unit]


def _parse_exported_inputs(gdml_path: Path, macro_path: Path) -> dict[str, Any]:
    """Independently reconstruct the identity of one exported G4 run."""

    gdml_text = gdml_path.read_text(encoding="utf-8")
    marker = "CEEFIT-GEOMETRY-v1"
    begin = gdml_text.find(marker)
    end = gdml_text.find("-->", begin + len(marker))
    if begin < 0 or end < 0:
        raise ValueError(f"{gdml_path}: missing CEEFIT-GEOMETRY-v1 metadata")
    metadata: dict[str, str] = {}
    for token in gdml_text[begin + len(marker) : end].split():
        key, separator, value = token.partition("=")
        if not separator:
            continue
        if key in metadata:
            raise ValueError(f"{gdml_path}: repeated CEEFIT geometry field {key!r}")
        metadata[key] = value
    required = {
        "geometry_tier",
        "material",
        "energy_keV",
        "distance_cm",
        "source_height_cm",
        "detector_height_cm",
        "scoring_radius_cm",
        "world_half_extent_cm",
    }
    if set(metadata) != required:
        raise ValueError(f"{gdml_path}: incomplete CEEFIT geometry metadata")

    try:
        xml_root = ElementTree.fromstring(gdml_text)
    except ElementTree.ParseError as error:
        raise ValueError(f"{gdml_path}: invalid GDML XML") from error
    solids: dict[str, tuple[str, Mapping[str, str]]] = {}
    for element in xml_root.iter():
        tag = element.tag.rsplit("}", 1)[-1]
        if tag in {"box", "sphere"} and "name" in element.attrib:
            solids[element.attrib["name"]] = (tag, element.attrib)
    try:
        world_tag, world = solids["world_s"]
        score_tag, score = solids["score_s"]
    except KeyError as error:
        raise ValueError(f"{gdml_path}: missing world_s or score_s solid") from error
    if world_tag != "box" or score_tag != "sphere":
        raise ValueError(f"{gdml_path}: world_s/score_s have the wrong solid type")
    world_unit = str(world.get("lunit", "cm"))
    world_half_cm = 0.5 * min(
        _length_cm(str(world[axis]), world_unit, f"{gdml_path}: world_s")
        for axis in ("x", "y", "z")
    )
    score_radius_cm = _length_cm(
        str(score.get("rmax", "")),
        str(score.get("lunit", "cm")),
        f"{gdml_path}: score_s",
    )

    volumes = {
        element.attrib["name"]: element
        for element in xml_root.iter()
        if element.tag.rsplit("}", 1)[-1] == "volume" and "name" in element.attrib
    }
    try:
        active_volume = volumes["active_crystal"]
        world_volume = volumes["world"]
    except KeyError as error:
        raise ValueError(f"{gdml_path}: missing active_crystal or world volume") from error

    def child_reference(element: Any, tag_name: str) -> str:
        matches = [
            child.attrib.get("ref", "")
            for child in element
            if child.tag.rsplit("}", 1)[-1] == tag_name
        ]
        if len(matches) != 1 or not matches[0]:
            raise ValueError(f"{gdml_path}: invalid {tag_name} reference")
        return str(matches[0])

    if (
        child_reference(active_volume, "materialref") != "G4_AIR"
        or child_reference(active_volume, "solidref") != "score_s"
        or child_reference(world_volume, "solidref") != "world_s"
    ):
        raise ValueError(f"{gdml_path}: scoring/world volume identity is invalid")

    placements: dict[str, tuple[float, float, float]] = {}
    for physical in world_volume:
        if physical.tag.rsplit("}", 1)[-1] != "physvol":
            continue
        volume_refs = [
            child.attrib.get("ref", "")
            for child in physical
            if child.tag.rsplit("}", 1)[-1] == "volumeref"
        ]
        positions = [
            child
            for child in physical
            if child.tag.rsplit("}", 1)[-1] == "position"
        ]
        if len(volume_refs) != 1 or len(positions) != 1:
            raise ValueError(f"{gdml_path}: invalid physical-volume placement")
        position = positions[0]
        unit = str(position.attrib.get("unit", "cm"))
        placements[str(volume_refs[0])] = tuple(
            _length_cm(str(position.attrib.get(axis, "0")), unit, f"{gdml_path}: position")
            for axis in ("x", "y", "z")
        )
    if "active_crystal" not in placements:
        raise ValueError(f"{gdml_path}: scoring sphere is not placed in the world")
    score_position_cm = placements["active_crystal"]

    geometry_tier = metadata["geometry_tier"]
    material = metadata["material"]
    numeric: dict[str, float] = {}
    for key in required - {"geometry_tier", "material"}:
        try:
            numeric[key] = float(metadata[key])
        except ValueError as error:
            raise ValueError(f"{gdml_path}: invalid CEEFIT geometry field {key}") from error
        if not math.isfinite(numeric[key]):
            raise ValueError(f"{gdml_path}: non-finite CEEFIT geometry field {key}")
    if not math.isclose(
        numeric["scoring_radius_cm"], score_radius_cm, rel_tol=1.0e-12, abs_tol=1.0e-12
    ) or not math.isclose(
        numeric["world_half_extent_cm"], world_half_cm, rel_tol=1.0e-12, abs_tol=1.0e-12
    ):
        raise ValueError(f"{gdml_path}: CEEFIT geometry does not close to GDML solids")

    macro_values: dict[str, list[str]] = {}
    for raw_line in macro_path.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split()
        if not fields or fields[0].startswith("#"):
            continue
        command = fields[0]
        if command.startswith("/random/"):
            raise ValueError(f"{macro_path}: formal anchor macro overrides the RNG")
        if command in {
            "/gps/energy",
            "/gps/ene/mono",
            "/gps/pos/centre",
            "/gps/ang/type",
            "/run/beamOn",
        }:
            if command in macro_values:
                raise ValueError(f"{macro_path}: repeated {command} command")
            macro_values[command] = fields[1:]
    energy_fields = macro_values.get("/gps/energy") or macro_values.get("/gps/ene/mono")
    position_fields = macro_values.get("/gps/pos/centre")
    angle_fields = macro_values.get("/gps/ang/type")
    event_fields = macro_values.get("/run/beamOn")
    if (
        energy_fields is None
        or len(energy_fields) != 2
        or position_fields is None
        or len(position_fields) != 4
        or angle_fields != ["iso"]
        or event_fields is None
        or len(event_fields) != 1
    ):
        raise ValueError(f"{macro_path}: incomplete isotropic monoenergetic GPS source")
    energy_scale = {"eV": 0.001, "keV": 1.0, "MeV": 1000.0}.get(energy_fields[1])
    if energy_scale is None:
        raise ValueError(f"{macro_path}: unsupported source-energy unit")
    try:
        energy_keV = float(energy_fields[0]) * energy_scale
        events = int(event_fields[0])
    except ValueError as error:
        raise ValueError(f"{macro_path}: invalid source energy or event count") from error
    source_cm = tuple(
        _length_cm(value, position_fields[3], f"{macro_path}: source position")
        for value in position_fields[:3]
    )
    source_offset_cm = tuple(
        source - score for source, score in zip(source_cm, score_position_cm)
    )
    distance_cm = math.sqrt(math.fsum(value * value for value in source_offset_cm))
    if geometry_tier == "air":
        derived_material = "none"
        derived_detector_height_cm = 0.0
        derived_source_height_cm = 0.0
        if any(abs(value) > 1.0e-12 for value in score_position_cm) or abs(
            source_offset_cm[1]
        ) > 1.0e-12:
            raise ValueError(f"{gdml_path}: canonical air geometry is not at 0/0 height")
    elif geometry_tier == "field":
        try:
            ground_volume = volumes["ground"]
            ground_position_cm = placements["ground"]
            ground_solid = child_reference(ground_volume, "solidref")
            ground_material = child_reference(ground_volume, "materialref")
            ground_tag, ground = solids[ground_solid]
        except KeyError as error:
            raise ValueError(f"{gdml_path}: field geometry lacks ground slab") from error
        if ground_tag != "box":
            raise ValueError(f"{gdml_path}: ground slab is not a box")
        ground_half_y_cm = 0.5 * _length_cm(
            str(ground.get("y", "")),
            str(ground.get("lunit", "cm")),
            f"{gdml_path}: ground slab",
        )
        ground_surface_y_cm = ground_position_cm[1] + ground_half_y_cm
        derived_detector_height_cm = score_position_cm[1] - ground_surface_y_cm
        derived_source_height_cm = source_cm[1] - ground_surface_y_cm
        derived_material = {
            "MatSoil": "soil",
            "MatConcrete": "concrete",
        }.get(ground_material, "")
    else:
        raise ValueError(f"{gdml_path}: formal anchor geometry must be air or field")
    if derived_material != material:
        raise ValueError(f"{gdml_path}: material metadata does not match GDML")
    if events < 2 or not math.isclose(
        energy_keV, numeric["energy_keV"], rel_tol=1.0e-12, abs_tol=1.0e-12
    ) or not math.isclose(
        distance_cm, numeric["distance_cm"], rel_tol=1.0e-9, abs_tol=1.0e-9
    ) or not math.isclose(
        derived_source_height_cm,
        numeric["source_height_cm"],
        rel_tol=1.0e-9,
        abs_tol=1.0e-9,
    ) or not math.isclose(
        derived_detector_height_cm,
        numeric["detector_height_cm"],
        rel_tol=1.0e-9,
        abs_tol=1.0e-9,
    ):
        raise ValueError(f"{macro_path}: source does not close to CEEFIT geometry")
    return {
        "geometry_tier": geometry_tier,
        "material": material,
        **numeric,
        "gdml_scoring_radius_cm": score_radius_cm,
        "world_half_extent_m": world_half_cm * 0.01,
        "events": events,
        "e0_keV": energy_keV,
        "primary_sampling": REQUIRED_PRIMARY_SAMPLING,
    }


def _parse_reference_arm(path: Path) -> dict[str, Any]:
    """Read one sphere-flux arm and verify its per-event total moments.

    Cell standard errors cannot be added in quadrature for the integrated
    sphere flux: one primary can score in many cells.  Formal anchor evidence
    therefore requires the scorer's raw moments of the *whole per-event*
    scattered track length and derives the integrated uncertainty from them.
    """

    events: int | None = None
    e0_keV: float | None = None
    sphere_radius_cm: float | None = None
    gdml_scoring_radius_cm: float | None = None
    energy_bins: int | None = None
    mu_bins: int | None = None
    energy_min_keV: float | None = None
    energy_max_keV: float | None = None
    mu_min: float | None = None
    mu_max: float | None = None
    observable: str | None = None
    physics_list: str | None = None
    world_half_extent_m: float | None = None
    primary_sampling: str | None = None
    rng_engine: str | None = None
    rng_seed: int | None = None
    production_cut_mm: float | None = None
    gdml_sha256: str | None = None
    macro_sha256: str | None = None
    geometry_tier: str | None = None
    material: str | None = None
    distance_cm: float | None = None
    source_height_cm: float | None = None
    detector_height_cm: float | None = None
    moment_schema: str | None = None
    scattered_total_n: int | None = None
    scattered_total_s1_cm: float | None = None
    scattered_total_s2_cm2: float | None = None
    stated_total_mean: float | None = None
    stated_total_standard_error: float | None = None
    rows: set[tuple[int, int]] = set()
    energy_lows: dict[int, float] = {}
    mu_lows: dict[int, float] = {}
    scattered_values: list[float] = []
    header_seen = False
    with path.open("r", encoding="utf-8", newline="") as source:
        for line_number, raw_line in enumerate(source, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                if body.startswith("events:"):
                    if events is not None:
                        raise ValueError(f"{path}: repeated events header")
                    try:
                        parsed = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid events header") from error
                    if not math.isfinite(parsed) or parsed <= 0.0 or not parsed.is_integer():
                        raise ValueError(f"{path}: events must be a positive integer")
                    events = int(parsed)
                elif body.startswith("E0_keV:"):
                    if e0_keV is not None:
                        raise ValueError(f"{path}: repeated E0_keV header")
                    try:
                        e0_keV = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid E0_keV header") from error
                    if not math.isfinite(e0_keV) or e0_keV <= 0.0:
                        raise ValueError(f"{path}: E0_keV must be finite and positive")
                elif body.startswith("sphere_radius_cm:"):
                    if sphere_radius_cm is not None:
                        raise ValueError(f"{path}: repeated sphere_radius_cm header")
                    try:
                        sphere_radius_cm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid sphere_radius_cm header") from error
                    if not math.isfinite(sphere_radius_cm) or sphere_radius_cm <= 0.0:
                        raise ValueError(
                            f"{path}: sphere_radius_cm must be finite and positive"
                        )
                elif body.startswith("gdml_scoring_radius_cm:"):
                    if gdml_scoring_radius_cm is not None:
                        raise ValueError(
                            f"{path}: repeated gdml_scoring_radius_cm header"
                        )
                    try:
                        gdml_scoring_radius_cm = float(
                            body.partition(":")[2].strip()
                        )
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid gdml_scoring_radius_cm header"
                        ) from error
                elif body.startswith("n_energy:"):
                    if energy_bins is not None:
                        raise ValueError(f"{path}: repeated n_energy header")
                    try:
                        energy_bins = int(body.partition(":")[2].strip().split()[0])
                    except (ValueError, IndexError) as error:
                        raise ValueError(f"{path}: invalid n_energy header") from error
                elif body.startswith("n_mu:"):
                    if mu_bins is not None:
                        raise ValueError(f"{path}: repeated n_mu header")
                    try:
                        mu_bins = int(body.partition(":")[2].strip().split()[0])
                    except (ValueError, IndexError) as error:
                        raise ValueError(f"{path}: invalid n_mu header") from error
                elif body.startswith("energy_min_keV:"):
                    if energy_min_keV is not None:
                        raise ValueError(f"{path}: repeated energy_min_keV header")
                    try:
                        energy_min_keV = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid energy_min_keV header") from error
                elif body.startswith("energy_max_keV:"):
                    if energy_max_keV is not None:
                        raise ValueError(f"{path}: repeated energy_max_keV header")
                    try:
                        energy_max_keV = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid energy_max_keV header") from error
                elif body.startswith("mu_min:"):
                    if mu_min is not None:
                        raise ValueError(f"{path}: repeated mu_min header")
                    try:
                        mu_min = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid mu_min header") from error
                elif body.startswith("mu_max:"):
                    if mu_max is not None:
                        raise ValueError(f"{path}: repeated mu_max header")
                    try:
                        mu_max = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid mu_max header") from error
                elif body.startswith("observable:"):
                    if observable is not None:
                        raise ValueError(f"{path}: repeated observable header")
                    observable = body.partition(":")[2].strip()
                elif body.startswith("physics_list:"):
                    if physics_list is not None:
                        raise ValueError(f"{path}: repeated physics_list header")
                    physics_list = body.partition(":")[2].strip()
                elif body.startswith("world_half_extent_m:"):
                    if world_half_extent_m is not None:
                        raise ValueError(f"{path}: repeated world_half_extent_m header")
                    try:
                        world_half_extent_m = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid world_half_extent_m header"
                        ) from error
                elif body.startswith("primary_sampling:"):
                    if primary_sampling is not None:
                        raise ValueError(f"{path}: repeated primary_sampling header")
                    primary_sampling = body.partition(":")[2].strip()
                elif body.startswith("rng_engine:"):
                    if rng_engine is not None:
                        raise ValueError(f"{path}: repeated rng_engine header")
                    rng_engine = body.partition(":")[2].strip()
                elif body.startswith("rng_seed:"):
                    if rng_seed is not None:
                        raise ValueError(f"{path}: repeated rng_seed header")
                    try:
                        parsed = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid rng_seed header") from error
                    if not math.isfinite(parsed) or parsed <= 0.0 or not parsed.is_integer():
                        raise ValueError(f"{path}: rng_seed must be a positive integer")
                    rng_seed = int(parsed)
                elif body.startswith("production_cut_mm:"):
                    if production_cut_mm is not None:
                        raise ValueError(f"{path}: repeated production_cut_mm header")
                    try:
                        production_cut_mm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid production_cut_mm header"
                        ) from error
                elif body.startswith("gdml_sha256:"):
                    if gdml_sha256 is not None:
                        raise ValueError(f"{path}: repeated gdml_sha256 header")
                    gdml_sha256 = body.partition(":")[2].strip()
                elif body.startswith("macro_sha256:"):
                    if macro_sha256 is not None:
                        raise ValueError(f"{path}: repeated macro_sha256 header")
                    macro_sha256 = body.partition(":")[2].strip()
                elif body.startswith("geometry_tier:"):
                    if geometry_tier is not None:
                        raise ValueError(f"{path}: repeated geometry_tier header")
                    geometry_tier = body.partition(":")[2].strip()
                elif body.startswith("material:"):
                    if material is not None:
                        raise ValueError(f"{path}: repeated material header")
                    material = body.partition(":")[2].strip()
                elif body.startswith("distance_cm:"):
                    if distance_cm is not None:
                        raise ValueError(f"{path}: repeated distance_cm header")
                    try:
                        distance_cm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(f"{path}: invalid distance_cm header") from error
                elif body.startswith("source_height_cm:"):
                    if source_height_cm is not None:
                        raise ValueError(f"{path}: repeated source_height_cm header")
                    try:
                        source_height_cm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid source_height_cm header"
                        ) from error
                elif body.startswith("detector_height_cm:"):
                    if detector_height_cm is not None:
                        raise ValueError(f"{path}: repeated detector_height_cm header")
                    try:
                        detector_height_cm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid detector_height_cm header"
                        ) from error
                elif body.startswith("scattered_total_moment_schema:"):
                    if moment_schema is not None:
                        raise ValueError(
                            f"{path}: repeated scattered_total_moment_schema header"
                        )
                    moment_schema = body.partition(":")[2].strip()
                elif body.startswith("scattered_total_N:"):
                    if scattered_total_n is not None:
                        raise ValueError(f"{path}: repeated scattered_total_N header")
                    try:
                        parsed = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid scattered_total_N header"
                        ) from error
                    if not math.isfinite(parsed) or parsed < 2.0 or not parsed.is_integer():
                        raise ValueError(
                            f"{path}: scattered_total_N must be an integer of at least two"
                        )
                    scattered_total_n = int(parsed)
                elif body.startswith("scattered_total_S1_cm:"):
                    if scattered_total_s1_cm is not None:
                        raise ValueError(f"{path}: repeated scattered_total_S1_cm header")
                    try:
                        scattered_total_s1_cm = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid scattered_total_S1_cm header"
                        ) from error
                elif body.startswith("scattered_total_S2_cm2:"):
                    if scattered_total_s2_cm2 is not None:
                        raise ValueError(f"{path}: repeated scattered_total_S2_cm2 header")
                    try:
                        scattered_total_s2_cm2 = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid scattered_total_S2_cm2 header"
                        ) from error
                elif body.startswith("scattered_total_mean_flux_cm-2:"):
                    if stated_total_mean is not None:
                        raise ValueError(
                            f"{path}: repeated scattered_total_mean_flux_cm-2 header"
                        )
                    try:
                        stated_total_mean = float(body.partition(":")[2].strip())
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid scattered total mean header"
                        ) from error
                elif body.startswith(
                    "scattered_total_standard_error_flux_cm-2:"
                ):
                    if stated_total_standard_error is not None:
                        raise ValueError(
                            f"{path}: repeated scattered total standard-error header"
                        )
                    try:
                        stated_total_standard_error = float(
                            body.partition(":")[2].strip()
                        )
                    except ValueError as error:
                        raise ValueError(
                            f"{path}: invalid scattered total standard-error header"
                        ) from error
                continue
            fields = next(csv.reader([raw_line]))
            if fields == [
                "energy_bin",
                "mu_bin",
                "e_low_keV",
                "mu_low",
                "phi_direct",
                "phi_scattered",
                "phi_scattered_sigma",
            ]:
                if header_seen:
                    raise ValueError(f"{path}: repeated column header")
                header_seen = True
                continue
            if not header_seen:
                raise ValueError(f"{path}:{line_number}: data precedes the column header")
            if len(fields) != 7:
                # Origin-decomposition rows are a distinct trailing record and
                # are irrelevant to the integrated anchor comparison.
                if len(fields) == 4:
                    continue
                raise ValueError(f"{path}:{line_number}: expected seven CSV fields")
            try:
                energy_bin = int(fields[0])
                mu_bin = int(fields[1])
                numeric = [float(field) for field in fields[2:]]
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid numeric field") from error
            if (
                energy_bin < 0
                or energy_bin >= REFERENCE_ENERGY_BINS
                or mu_bin < 0
                or mu_bin >= REFERENCE_MU_BINS
                or any(not math.isfinite(value) for value in numeric)
            ):
                raise ValueError(f"{path}:{line_number}: invalid cell coordinates or value")
            if numeric[0] <= 0.0 or numeric[1] < -1.0 or numeric[1] >= 1.0:
                raise ValueError(f"{path}:{line_number}: invalid energy or mu lower edge")
            if any(value < 0.0 for value in numeric[2:]):
                raise ValueError(f"{path}:{line_number}: fluence and uncertainty must be nonnegative")
            cell = (energy_bin, mu_bin)
            if cell in rows:
                raise ValueError(f"{path}:{line_number}: duplicate energy/mu cell")
            rows.add(cell)
            if energy_bin in energy_lows and energy_lows[energy_bin] != numeric[0]:
                raise ValueError(f"{path}:{line_number}: inconsistent energy lower edge")
            if mu_bin in mu_lows and mu_lows[mu_bin] != numeric[1]:
                raise ValueError(f"{path}:{line_number}: inconsistent mu lower edge")
            energy_lows[energy_bin] = numeric[0]
            mu_lows[mu_bin] = numeric[1]
            scattered_values.append(numeric[3])
    if (
        events is None
        or e0_keV is None
        or sphere_radius_cm is None
        or gdml_scoring_radius_cm is None
        or energy_bins is None
        or mu_bins is None
        or energy_min_keV is None
        or energy_max_keV is None
        or mu_min is None
        or mu_max is None
        or observable is None
        or physics_list is None
        or world_half_extent_m is None
        or primary_sampling is None
        or rng_engine is None
        or rng_seed is None
        or production_cut_mm is None
        or gdml_sha256 is None
        or macro_sha256 is None
        or geometry_tier is None
        or material is None
        or distance_cm is None
        or source_height_cm is None
        or detector_height_cm is None
        or moment_schema is None
        or scattered_total_n is None
        or scattered_total_s1_cm is None
        or scattered_total_s2_cm2 is None
        or stated_total_mean is None
        or stated_total_standard_error is None
        or not header_seen
        or not rows
    ):
        raise ValueError(f"{path}: incomplete GEANT4 sphere-flux record")
    if energy_bins != REFERENCE_ENERGY_BINS or mu_bins != REFERENCE_MU_BINS:
        raise ValueError(
            f"{path}: GEANT4 sphere-flux layout must be "
            f"{REFERENCE_ENERGY_BINS}x{REFERENCE_MU_BINS}"
        )
    if observable != REQUIRED_OBSERVABLE:
        raise ValueError(f"{path}: wrong sphere-flux observable")
    if physics_list != REQUIRED_PHYSICS_LIST:
        raise ValueError(f"{path}: physics_list must be {REQUIRED_PHYSICS_LIST}")
    if world_half_extent_m < MINIMUM_FIELD_WORLD_HALF_EXTENT_M:
        raise ValueError(
            f"{path}: world half extent must be at least "
            f"{MINIMUM_FIELD_WORLD_HALF_EXTENT_M:g} m"
        )
    if primary_sampling != REQUIRED_PRIMARY_SAMPLING:
        raise ValueError(f"{path}: sphere-flux sampling must be unbiased isotropic 4pi")
    if not rng_engine:
        raise ValueError(f"{path}: rng_engine must be nonempty")
    if not math.isfinite(production_cut_mm) or not 0.0 < production_cut_mm <= 0.7:
        raise ValueError(f"{path}: production_cut_mm must be in (0, 0.7]")
    _sha256(gdml_sha256, "sphere-flux GDML SHA-256")
    _sha256(macro_sha256, "sphere-flux macro SHA-256")
    if geometry_tier not in {"air", "field"}:
        raise ValueError(f"{path}: formal anchors support only air or field geometry")
    if (geometry_tier == "air" and material != "none") or (
        geometry_tier == "field" and material not in {"soil", "concrete"}
    ):
        raise ValueError(f"{path}: geometry tier/material metadata is inconsistent")
    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (distance_cm, source_height_cm, detector_height_cm)
    ) or distance_cm <= abs(source_height_cm - detector_height_cm):
        raise ValueError(f"{path}: invalid sphere-flux geometry coordinates")
    if not math.isclose(
        sphere_radius_cm,
        gdml_scoring_radius_cm,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ValueError(f"{path}: scorer radius does not match the hashed GDML radius")
    if (
        not math.isfinite(energy_min_keV)
        or not math.isfinite(energy_max_keV)
        or energy_min_keV <= 0.0
        or energy_max_keV < e0_keV
        or energy_max_keV <= energy_min_keV
        or mu_min != -1.0
        or mu_max != 1.0
    ):
        raise ValueError(f"{path}: invalid sphere-flux energy/mu domain")
    expected_rows = {
        (energy_bin, mu_bin)
        for energy_bin in range(REFERENCE_ENERGY_BINS)
        for mu_bin in range(REFERENCE_MU_BINS)
    }
    if rows != expected_rows:
        raise ValueError(f"{path}: GEANT4 sphere-flux grid is incomplete")
    for index in range(REFERENCE_ENERGY_BINS):
        expected = energy_min_keV * (
            energy_max_keV / energy_min_keV
        ) ** (index / REFERENCE_ENERGY_BINS)
        if not math.isclose(
            energy_lows[index], expected, rel_tol=5.0e-12, abs_tol=0.0
        ):
            raise ValueError(f"{path}: GEANT4 sphere-flux energy grid is inconsistent")
    for index in range(REFERENCE_MU_BINS):
        expected = mu_min + (mu_max - mu_min) * index / REFERENCE_MU_BINS
        if not math.isclose(mu_lows[index], expected, rel_tol=0.0, abs_tol=5.0e-15):
            raise ValueError(f"{path}: GEANT4 sphere-flux mu grid is inconsistent")
    if any(
        right <= left
        for left, right in zip(
            (energy_lows[index] for index in range(REFERENCE_ENERGY_BINS - 1)),
            (energy_lows[index] for index in range(1, REFERENCE_ENERGY_BINS)),
        )
    ) or any(
        right <= left
        for left, right in zip(
            (mu_lows[index] for index in range(REFERENCE_MU_BINS - 1)),
            (mu_lows[index] for index in range(1, REFERENCE_MU_BINS)),
        )
    ):
        raise ValueError(f"{path}: GEANT4 sphere-flux edges are not strictly increasing")
    if moment_schema != "per-event-track-length-v1":
        raise ValueError(f"{path}: unsupported scattered total moment schema")
    if scattered_total_n != events:
        raise ValueError(f"{path}: scattered_total_N does not equal events")
    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (
            scattered_total_s1_cm,
            scattered_total_s2_cm2,
            stated_total_mean,
            stated_total_standard_error,
        )
    ):
        raise ValueError(f"{path}: scattered total moments must be finite and nonnegative")
    centered = (
        scattered_total_s2_cm2
        - scattered_total_s1_cm * scattered_total_s1_cm / scattered_total_n
    )
    centered_scale = max(
        abs(scattered_total_s2_cm2),
        abs(scattered_total_s1_cm * scattered_total_s1_cm / scattered_total_n),
        1.0,
    )
    if centered < -64.0 * math.ulp(centered_scale):
        raise ValueError(f"{path}: scattered total moments have negative variance")
    volume = (4.0 / 3.0) * math.pi * sphere_radius_cm**3
    scattered_mean = scattered_total_s1_cm / (scattered_total_n * volume)
    scattered_standard_error = math.sqrt(
        max(0.0, centered) / (scattered_total_n * (scattered_total_n - 1))
    ) / volume

    def require_close(stated: float, computed: float, name: str) -> None:
        absolute_tolerance = 64.0 * math.ulp(max(abs(stated), abs(computed)))
        if not math.isclose(
            stated,
            computed,
            rel_tol=5.0e-12,
            abs_tol=absolute_tolerance,
        ):
            raise ValueError(f"{path}: {name} does not close to raw N/S1/S2")

    require_close(stated_total_mean, scattered_mean, "scattered total mean")
    require_close(
        stated_total_standard_error,
        scattered_standard_error,
        "scattered total standard error",
    )
    cell_sum = math.fsum(scattered_values)
    require_close(cell_sum, scattered_mean, "scattered cell sum")
    return {
        "events": events,
        "e0_keV": e0_keV,
        "sphere_radius_cm": sphere_radius_cm,
        "gdml_scoring_radius_cm": gdml_scoring_radius_cm,
        "energy_bins": energy_bins,
        "mu_bins": mu_bins,
        "energy_min_keV": energy_min_keV,
        "energy_max_keV": energy_max_keV,
        "mu_min": mu_min,
        "mu_max": mu_max,
        "observable": observable,
        "physics_list": physics_list,
        "world_half_extent_m": world_half_extent_m,
        "primary_sampling": primary_sampling,
        "rng_engine": rng_engine,
        "rng_seed": rng_seed,
        "production_cut_mm": production_cut_mm,
        "gdml_sha256": gdml_sha256,
        "macro_sha256": macro_sha256,
        "geometry_tier": geometry_tier,
        "material": material,
        "distance_cm": distance_cm,
        "source_height_cm": source_height_cm,
        "detector_height_cm": detector_height_cm,
        "cell_count": len(rows),
        "scattered_total_moment_schema": moment_schema,
        "scattered_total_N": scattered_total_n,
        "scattered_total_S1_cm": scattered_total_s1_cm,
        "scattered_total_S2_cm2": scattered_total_s2_cm2,
        "scattered_mean": scattered_mean,
        "scattered_standard_error": scattered_standard_error,
    }


def _merge_reference_arms(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if len(records) < 2:
        raise ValueError(
            "GEANT4 anchor calibration requires at least two independent reference arms"
        )
    if len({int(record["cell_count"]) for record in records}) != 1:
        raise ValueError("GEANT4 reference arms do not share one cell layout")
    for field in (
        "e0_keV",
        "sphere_radius_cm",
        "gdml_scoring_radius_cm",
        "energy_bins",
        "mu_bins",
        "energy_min_keV",
        "energy_max_keV",
        "mu_min",
        "mu_max",
        "observable",
        "physics_list",
        "world_half_extent_m",
        "primary_sampling",
        "rng_engine",
        "production_cut_mm",
        "geometry_tier",
        "material",
        "distance_cm",
        "source_height_cm",
        "detector_height_cm",
    ):
        if len({record[field] for record in records}) != 1:
            raise ValueError(f"GEANT4 reference arms disagree on {field}")
    if len({int(record["rng_seed"]) for record in records}) != len(records):
        raise ValueError("GEANT4 reference arms must use distinct RNG seeds")
    total_events = sum(int(record["events"]) for record in records)
    if total_events <= 0:
        raise ValueError("GEANT4 anchor has no reference events")
    total_s1_cm = math.fsum(float(record["scattered_total_S1_cm"]) for record in records)
    total_s2_cm2 = math.fsum(float(record["scattered_total_S2_cm2"]) for record in records)
    radius_cm = float(records[0]["sphere_radius_cm"])
    volume = (4.0 / 3.0) * math.pi * radius_cm**3
    pooled_centered = total_s2_cm2 - total_s1_cm * total_s1_cm / total_events
    centered_scale = max(
        abs(total_s2_cm2), abs(total_s1_cm * total_s1_cm / total_events), 1.0
    )
    if pooled_centered < -64.0 * math.ulp(centered_scale):
        raise ValueError("merged GEANT4 reference moments have negative variance")
    reference_mean = total_s1_cm / (total_events * volume)
    per_history_se = math.sqrt(
        max(0.0, pooled_centered) / (total_events * (total_events - 1))
    ) / volume
    weights = [float(record["events"]) for record in records]
    means = [float(record["scattered_mean"]) for record in records]
    weight_sum = math.fsum(weights)
    squared_weight_sum = math.fsum(weight * weight for weight in weights)
    denominator = weight_sum - squared_weight_sum / weight_sum
    effective_arms = weight_sum * weight_sum / squared_weight_sum
    if denominator <= 0.0 or effective_arms <= 1.0:
        raise ValueError("GEANT4 reference arms have insufficient independent weight")
    between_variance = math.fsum(
        weight * (mean - reference_mean) ** 2
        for weight, mean in zip(weights, means)
    ) / denominator
    between_arm_se = math.sqrt(max(0.0, between_variance) / effective_arms)
    return {
        "events": total_events,
        "arm_count": len(records),
        "e0_keV": float(records[0]["e0_keV"]),
        "sphere_radius_cm": float(records[0]["sphere_radius_cm"]),
        "gdml_scoring_radius_cm": float(records[0]["gdml_scoring_radius_cm"]),
        "energy_bins": int(records[0]["energy_bins"]),
        "mu_bins": int(records[0]["mu_bins"]),
        "energy_min_keV": float(records[0]["energy_min_keV"]),
        "energy_max_keV": float(records[0]["energy_max_keV"]),
        "mu_min": float(records[0]["mu_min"]),
        "mu_max": float(records[0]["mu_max"]),
        "observable": str(records[0]["observable"]),
        "physics_list": str(records[0]["physics_list"]),
        "world_half_extent_m": float(records[0]["world_half_extent_m"]),
        "primary_sampling": str(records[0]["primary_sampling"]),
        "rng_engine": str(records[0]["rng_engine"]),
        "production_cut_mm": float(records[0]["production_cut_mm"]),
        "geometry_tier": str(records[0]["geometry_tier"]),
        "material": str(records[0]["material"]),
        "distance_cm": float(records[0]["distance_cm"]),
        "source_height_cm": float(records[0]["source_height_cm"]),
        "detector_height_cm": float(records[0]["detector_height_cm"]),
        "mean": reference_mean,
        "scattered_total_moment_schema": "per-event-track-length-v1",
        "scattered_total_N": total_events,
        "scattered_total_S1_cm": total_s1_cm,
        "scattered_total_S2_cm2": total_s2_cm2,
        "per_history_standard_error": per_history_se,
        "between_arm_standard_error": between_arm_se,
        "standard_error": max(per_history_se, between_arm_se),
    }


def _archive_scalar(root: Any, index: int) -> dict[str, float | int]:
    count = int(root["diagnostics/scattered_total_N"][index])
    first = float(root["diagnostics/scattered_total_S1"][index])
    second = float(root["diagnostics/scattered_total_S2"][index])
    if count < 2 or not all(math.isfinite(value) for value in (first, second)):
        raise ValueError("GEANT4 anchor archive node lacks usable scalar moments")
    numerator = second - first * first / count
    tolerance = 64.0 * math.ulp(max(abs(second), abs(first * first / count), 1.0))
    if numerator < -tolerance:
        raise ValueError("GEANT4 anchor archive node has negative centered variance")
    standard_error = math.sqrt(max(0.0, numerator) / (count * (count - 1)))
    return {"histories": count, "mean": first / count, "standard_error": standard_error}


def generate_geant4_anchor_report(
    archive_path: str | Path,
    manifest_path: str | Path,
    output_path: str | Path,
) -> dict[str, Any]:
    """Generate an identity-bound, independently recomputable anchor report."""

    archive_file = Path(archive_path).resolve()
    manifest_file = Path(manifest_path).resolve()
    output_file = Path(output_path).resolve()
    if output_file in {archive_file, manifest_file}:
        raise ValueError("GEANT4 anchor output must not overwrite an input artifact")
    manifest = _strict_json(manifest_file)
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unsupported GEANT4 anchor manifest schema")
    if (
        "physics_list" in manifest
        and manifest.get("physics_list") != REQUIRED_PHYSICS_LIST
    ):
        raise ValueError(
            f"GEANT4 anchor physics_list must be {REQUIRED_PHYSICS_LIST}"
        )
    comparison_sigma = _finite_nonnegative(
        manifest.get("comparison_sigma", DEFAULT_COMPARISON_SIGMA),
        "comparison_sigma",
    )
    if comparison_sigma != DEFAULT_COMPARISON_SIGMA:
        raise ValueError(f"comparison_sigma must equal {DEFAULT_COMPARISON_SIGMA}")
    anchors = manifest.get("anchors")
    if not isinstance(anchors, list) or not anchors:
        raise ValueError("GEANT4 anchor manifest must contain at least one anchor")

    h5py, _ = require_hdf5()
    report_anchors: list[dict[str, Any]] = []
    seen_anchor_ids: set[str] = set()
    seen_node_ids: set[str] = set()
    seen_reference_files: set[Path] = set()
    seen_reference_hashes: set[str] = set()
    with h5py.File(archive_file, "r") as root:
        config_hash = _sha256(
            root["manifest"].attrs.get("config_sha256", ""), "archive config_sha256"
        )
        snapshot_hash = archive_snapshot_sha256(root)
        canonical_text = root["manifest"].attrs.get("canonical_config_json")
        if not isinstance(canonical_text, str):
            raise ValueError("fit archive lacks canonical configuration JSON")
        try:
            archive_configuration = json.loads(canonical_text)
            archive_axis = archive_configuration["axis"]
            configured_energy_cut_keV = float(archive_axis["energy_cut_keV"])
            configured_scoring_radius_cm = float(
                archive_configuration["scoring_radius_cm"]
            )
            configured_geometry_tier = str(
                archive_configuration["geometry_tier"]
            )
            configured_material_record = archive_configuration["material"]
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise ValueError(
                "fit archive configuration lacks formal GEANT4 observable metadata"
            ) from error
        if configured_geometry_tier not in {"air", "field"}:
            raise ValueError("formal GEANT4 anchors support only air or field archives")
        if configured_geometry_tier == "air":
            configured_material = "none"
            if configured_material_record is not None:
                raise ValueError("air archive has unexpected field material")
        else:
            if not isinstance(configured_material_record, Mapping):
                raise ValueError("field archive lacks a material definition")
            material_name = str(configured_material_record.get("name", "")).lower()
            if "concrete" in material_name:
                configured_material = "concrete"
            elif "soil" in material_name:
                configured_material = "soil"
            else:
                raise ValueError("field archive has unsupported GEANT4 material")
        node_index = {
            f"{int(value):016x}": index
            for index, value in enumerate(root["nodes/id"][:])
        }
        if len(node_index) != len(root["nodes/id"]):
            raise ValueError("fit archive repeats a node id")
        for raw_anchor in anchors:
            if not isinstance(raw_anchor, Mapping):
                raise ValueError("GEANT4 anchor entries must be JSON objects")
            anchor_id = str(raw_anchor.get("anchor_id", ""))
            if not anchor_id or anchor_id in seen_anchor_ids:
                raise ValueError("GEANT4 anchor_id values must be nonempty and unique")
            seen_anchor_ids.add(anchor_id)
            node_id = _node_hex(raw_anchor.get("node_id", ""))
            if node_id in seen_node_ids:
                raise ValueError("GEANT4 anchor manifest repeats an archive node")
            seen_node_ids.add(node_id)
            if node_id not in node_index:
                raise ValueError(f"GEANT4 anchor node {node_id} is absent from the archive")
            references = raw_anchor.get("reference_files")
            if not isinstance(references, list) or not references:
                raise ValueError(f"GEANT4 anchor {anchor_id!r} has no reference files")
            arm_records: list[dict[str, Any]] = []
            source_records: list[dict[str, Any]] = []
            for raw_reference in references:
                if not isinstance(raw_reference, Mapping):
                    raise ValueError("GEANT4 reference file entries must be JSON objects")
                relative = Path(str(raw_reference.get("path", "")))
                if not str(relative) or relative.is_absolute():
                    raise ValueError("GEANT4 reference paths must be nonempty and relative")
                reference_file = (manifest_file.parent / relative).resolve()
                try:
                    relative_to_manifest = reference_file.relative_to(manifest_file.parent)
                except ValueError as error:
                    raise ValueError("GEANT4 reference path escapes the manifest directory") from error
                if reference_file == output_file:
                    raise ValueError(
                        "GEANT4 anchor output must not overwrite a reference artifact"
                    )
                if reference_file in seen_reference_files:
                    raise ValueError("a GEANT4 reference arm is assigned more than once")
                seen_reference_files.add(reference_file)
                actual_hash = file_sha256(reference_file)
                expected_hash = _sha256(
                    raw_reference.get("sha256", ""), "GEANT4 reference SHA-256"
                )
                if actual_hash != expected_hash:
                    raise ValueError(f"GEANT4 reference hash differs for {relative_to_manifest}")
                if actual_hash in seen_reference_hashes:
                    raise ValueError(
                        "GEANT4 reference arms must contain independent source bytes"
                    )
                seen_reference_hashes.add(actual_hash)
                arm = _parse_reference_arm(reference_file)
                bound_inputs: dict[str, dict[str, str]] = {}
                input_paths: dict[str, Path] = {}
                for input_name in ("gdml", "macro"):
                    path_key = f"{input_name}_path"
                    hash_key = f"{input_name}_sha256"
                    relative_input = Path(str(raw_reference.get(path_key, "")))
                    if not str(relative_input) or relative_input.is_absolute():
                        raise ValueError(
                            f"GEANT4 {input_name} paths must be nonempty and relative"
                        )
                    input_file = (manifest_file.parent / relative_input).resolve()
                    try:
                        relative_input = input_file.relative_to(manifest_file.parent)
                    except ValueError as error:
                        raise ValueError(
                            f"GEANT4 {input_name} path escapes the manifest directory"
                        ) from error
                    actual_input_hash = file_sha256(input_file)
                    expected_input_hash = _sha256(
                        raw_reference.get(hash_key, ""),
                        f"GEANT4 {input_name} SHA-256",
                    )
                    if actual_input_hash != expected_input_hash:
                        raise ValueError(
                            f"GEANT4 {input_name} hash differs for {relative_input}"
                        )
                    if actual_input_hash != arm[hash_key]:
                        raise ValueError(
                            f"GEANT4 sphere-flux header does not bind the exact {input_name}"
                        )
                    input_paths[input_name] = input_file
                    bound_inputs[input_name] = {
                        "path": relative_input.as_posix(),
                        "sha256": actual_input_hash,
                    }
                exported = _parse_exported_inputs(
                    input_paths["gdml"], input_paths["macro"]
                )
                for field in (
                    "geometry_tier",
                    "material",
                    "events",
                    "e0_keV",
                    "distance_cm",
                    "source_height_cm",
                    "detector_height_cm",
                    "gdml_scoring_radius_cm",
                    "world_half_extent_m",
                    "primary_sampling",
                ):
                    left = arm[field]
                    right = exported[field]
                    if isinstance(right, float):
                        matches = math.isclose(
                            float(left), right, rel_tol=1.0e-10, abs_tol=1.0e-10
                        )
                    else:
                        matches = left == right
                    if not matches:
                        raise ValueError(
                            f"GEANT4 sphere-flux {field} does not close to hashed inputs"
                        )
                arm_records.append(arm)
                source_records.append(
                    {
                        "path": relative_to_manifest.as_posix(),
                        "sha256": actual_hash,
                        "inputs": bound_inputs,
                        **arm,
                    }
                )
            reference = _merge_reference_arms(arm_records)
            index = node_index[node_id]
            archive_e0_keV = float(root["nodes/coordinates"][index][0])
            if not math.isclose(
                float(reference["e0_keV"]),
                archive_e0_keV,
                rel_tol=1.0e-8,
                abs_tol=1.0e-12,
            ):
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} E0 does not match its archive node"
                )
            if reference["geometry_tier"] != configured_geometry_tier:
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} geometry tier does not match archive"
                )
            if reference["material"] != configured_material:
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} material does not match archive"
                )
            if not math.isclose(
                float(reference["energy_min_keV"]),
                configured_energy_cut_keV,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ):
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} energy cutoff does not match archive"
                )
            if not math.isclose(
                float(reference["sphere_radius_cm"]),
                configured_scoring_radius_cm,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ):
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} scoring radius does not match archive"
                )
            archive_coordinates = [
                float(value) for value in root["nodes/coordinates"][index]
            ]
            for coordinate_index, field in enumerate(
                ("e0_keV", "distance_cm", "source_height_cm", "detector_height_cm")
            ):
                if not math.isclose(
                    float(reference[field]),
                    archive_coordinates[coordinate_index],
                    rel_tol=1.0e-9,
                    abs_tol=1.0e-9,
                ):
                    raise ValueError(
                        f"GEANT4 anchor {anchor_id!r} {field} does not match archive node"
                    )
            archive_mu_edges = root["axes/mu_edges"]
            if not math.isclose(float(reference["mu_min"]), float(archive_mu_edges[0])) or not math.isclose(
                float(reference["mu_max"]), float(archive_mu_edges[-1])
            ):
                raise ValueError(
                    f"GEANT4 anchor {anchor_id!r} angular domain does not match archive"
                )
            archive_result = _archive_scalar(root, index)
            archive_mean = float(archive_result["mean"])
            reference_mean = float(reference["mean"])
            difference = archive_mean - reference_mean
            combined_se = math.hypot(
                float(archive_result["standard_error"]),
                float(reference["standard_error"]),
            )
            if combined_se == 0.0:
                standardized = 0.0 if difference == 0.0 else math.copysign(1.0e100, difference)
            else:
                standardized = difference / combined_se
            relative_difference = (
                abs(difference) / abs(reference_mean)
                if reference_mean != 0.0
                else (0.0 if difference == 0.0 else 1.0e100)
            )
            passed = abs(standardized) <= comparison_sigma
            role_value = root["nodes/role"].asstr()[index]
            status_value = root["nodes/status"].asstr()[index]
            # These are the actual materialized archive states.  A capped node
            # remains usable with its measured uncertainty; only transport
            # failure makes an otherwise statistically consistent anchor
            # ineligible.
            if (
                status_value not in {"complete", "converged"}
                or int(root["nodes/capped"][index]) != 0
                or int(root["nodes/transport_failures"][index]) != 0
                or str(root["nodes/failure_reason"].asstr()[index])
            ):
                passed = False
            report_anchors.append(
                {
                    "anchor_id": anchor_id,
                    "node_id": node_id,
                    "coordinates": archive_coordinates,
                    "archive_role": str(role_value),
                    "archive_status": str(status_value),
                    "world_half_extent_m": float(reference["world_half_extent_m"]),
                    "reference": {**reference, "files": source_records},
                    "archive": archive_result,
                    "difference": difference,
                    "absolute_relative_difference": relative_difference,
                    "combined_standard_error": combined_se,
                    "standardized_difference": standardized,
                    "passed": bool(passed),
                }
            )

    within = all(anchor["passed"] for anchor in report_anchors)
    report = {
        "schema": REPORT_SCHEMA,
        "config_sha256": config_hash,
        "archive_snapshot_sha256": snapshot_hash,
        "passed": bool(report_anchors) and within,
        "anchor_count": len(report_anchors),
        "within_combined_uncertainty": bool(report_anchors) and within,
        "physics_list": REQUIRED_PHYSICS_LIST,
        "comparison_sigma": comparison_sigma,
        "minimum_field_world_half_extent_m": MINIMUM_FIELD_WORLD_HALF_EXTENT_M,
        "manifest_path": Path(
            os.path.relpath(manifest_file, output_file.parent)
        ).as_posix(),
        "manifest_sha256": file_sha256(manifest_file),
        "reference_file_count": sum(
            len(anchor["reference"]["files"]) for anchor in report_anchors
        ),
        "reference_event_count": sum(
            int(anchor["reference"]["events"]) for anchor in report_anchors
        ),
        "max_abs_standardized_difference": max(
            (abs(float(anchor["standardized_difference"])) for anchor in report_anchors),
            default=1.0e100,
        ),
        "max_absolute_relative_difference": max(
            (float(anchor["absolute_relative_difference"]) for anchor in report_anchors),
            default=1.0e100,
        ),
        "anchors": report_anchors,
    }
    _atomic_json(output_file, report)
    return report
