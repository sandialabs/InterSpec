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

"""Deterministic scattered designs for environmental-fluence campaigns."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import struct
from typing import Iterable, Mapping, Sequence


def stable_node_id_from_cm(
    config_hash: str,
    energy_keV: float,
    distance_cm: float,
    source_height_cm: float,
    detector_height_cm: float,
    *,
    extent_cm: float = 0.0,
) -> int:
    """Return the CEEFIT-NODE-v1 id from canonical stored coordinates.

    Worker blocks store lengths in centimetres and hash those exact binary64
    values.  Archive readers must therefore hash the stored centimetre values
    directly: converting cm -> m -> cm can change the final bit and spuriously
    reject an otherwise canonical block.
    """

    try:
        encoded_hash = config_hash.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("config_hash must be ASCII") from error
    coordinates = (
        energy_keV,
        distance_cm,
        source_height_cm,
        detector_height_cm,
        extent_cm,
    )
    if not all(math.isfinite(value) for value in coordinates):
        raise ValueError("node coordinates must be finite")
    if extent_cm < 0.0:
        raise ValueError("extent_cm must be non-negative")
    packed = struct.pack("<5d", *coordinates)
    digest = hashlib.sha256(b"CEEFIT-NODE-v1\0" + encoded_hash + packed).digest()
    return int.from_bytes(digest[:8], "big")


@dataclass(frozen=True)
class PhysicalNode:
    """One physical simulation coordinate."""

    energy_keV: float
    distance_m: float
    source_height_m: float
    detector_height_m: float

    def validate(self) -> None:
        values = (
            self.energy_keV,
            self.distance_m,
            self.source_height_m,
            self.detector_height_m,
        )
        if not all(math.isfinite(v) for v in values):
            raise ValueError("node coordinates must be finite")
        if self.energy_keV <= 0.0 or self.distance_m <= 0.0:
            raise ValueError("energy and source-detector distance must be positive")
        if self.source_height_m < 0.0 or self.detector_height_m < 0.0:
            raise ValueError("heights must be non-negative")
        if abs(self.source_height_m - self.detector_height_m) >= self.distance_m:
            raise ValueError("invalid geometry: distance must exceed |hsrc-hdet|")

    def as_dict(self) -> dict[str, float]:
        return {
            "energy_keV": self.energy_keV,
            "distance_m": self.distance_m,
            "source_height_m": self.source_height_m,
            "detector_height_m": self.detector_height_m,
        }

    def stable_id(self, config_hash: str, *, extent_m: float = 0.0) -> int:
        """Return the CEEFIT-NODE-v1 unsigned-64 identifier.

        The byte contract matches the C++ block writer: a domain separator,
        ASCII configuration SHA-256, then five little-endian IEEE-754 doubles.
        Geometry coordinates are converted from this API's metres to the C++
        transport interface's centimetres.
        """

        if not math.isfinite(extent_m) or extent_m < 0.0:
            raise ValueError("extent_m must be finite and non-negative")
        return stable_node_id_from_cm(
            config_hash,
            self.energy_keV,
            self.distance_m * 100.0,
            self.source_height_m * 100.0,
            self.detector_height_m * 100.0,
            extent_cm=extent_m * 100.0,
        )

    def stable_id_hex(self, config_hash: str, *, extent_m: float = 0.0) -> str:
        return f"{self.stable_id(config_hash, extent_m=extent_m):016x}"


@dataclass(frozen=True)
class DesignBounds:
    """Campaign bounds in physical space plus a finite transform margin.

    The valid physical set is open at ``rho = 0`` while the fitted coordinate
    contains ``log(rho)``.  A finite campaign therefore has to declare the
    smallest transverse separation it supports.  Keeping that margin explicit
    prevents a sampled point cloud from accidentally defining the runtime
    domain.
    """

    energy_keV: tuple[float, float] = (59.5, 2614.0)
    distance_m: tuple[float, float] = (0.5, 49.87)
    source_height_m: tuple[float, float] = (0.3, 4.16)
    detector_height_m: tuple[float, float] = (0.3, 4.16)
    rho_min_m: float = 0.01

    def validate(self) -> None:
        for name, pair in (
            ("energy_keV", self.energy_keV),
            ("distance_m", self.distance_m),
            ("source_height_m", self.source_height_m),
            ("detector_height_m", self.detector_height_m),
        ):
            if len(pair) != 2 or not all(math.isfinite(v) for v in pair):
                raise ValueError(f"{name} must contain two finite bounds")
            if pair[0] >= pair[1]:
                raise ValueError(f"{name} lower bound must be below upper bound")
        if self.energy_keV[0] <= 0.0 or self.distance_m[0] <= 0.0:
            raise ValueError("energy and distance bounds must be positive")
        if self.source_height_m[0] < 0.0 or self.detector_height_m[0] < 0.0:
            raise ValueError("height bounds must be non-negative")
        if (
            not math.isfinite(self.rho_min_m)
            or self.rho_min_m <= 0.0
            or self.rho_min_m >= self.distance_m[1]
        ):
            raise ValueError("rho_min_m must lie strictly between zero and maximum distance")


@dataclass(frozen=True)
class DesignSelection:
    train: tuple[PhysicalNode, ...]
    validation: tuple[PhysicalNode, ...]
    final: tuple[PhysicalNode, ...]

    def as_dict(self, config_hash: str) -> dict[str, list[dict[str, object]]]:
        def encode(role: str, nodes: Sequence[PhysicalNode]) -> list[dict[str, object]]:
            return [
                {
                    **node.as_dict(),
                    "q": list(geometry_transform(node)),
                    "role": role,
                    "node_id": node.stable_id(config_hash),
                    "node_id_hex": node.stable_id_hex(config_hash),
                }
                for node in nodes
            ]

        return {
            "train": encode("train", self.train),
            "validation": encode("validation", self.validation),
            "final": encode("final", self.final),
        }


def geometry_transform(node: PhysicalNode) -> tuple[float, float, float, float]:
    """Map physical inputs to the non-reciprocal fit coordinate ``q``."""

    node.validate()
    dz = node.source_height_m - node.detector_height_m
    rho2 = node.distance_m * node.distance_m - dz * dz
    # Validation makes rho2 positive; max protects against a one-ulp roundoff.
    rho = math.sqrt(max(rho2, math.ulp(node.distance_m * node.distance_m)))
    return (
        math.log(node.energy_keV),
        math.log(rho),
        math.atan2(node.source_height_m, rho),
        math.atan2(node.detector_height_m, rho),
    )


def transformed_domain(
    bounds: DesignBounds, *, mode: str = "field"
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Return a conservative, explicit bounding box for the fitted ``q`` domain."""

    bounds.validate()
    if mode not in {"air", "field"}:
        raise ValueError("transformed campaign domain supports air or field geometry")
    if mode == "air":
        return (
            (math.log(bounds.energy_keV[0]), math.log(bounds.distance_m[0]), 0.0, 0.0),
            (
                math.log(bounds.energy_keV[1]),
                math.log(bounds.distance_m[1]),
                1.0e-12,
                1.0e-12,
            ),
        )
    rho_low = bounds.rho_min_m
    rho_high = bounds.distance_m[1]
    return (
        (
            math.log(bounds.energy_keV[0]),
            math.log(rho_low),
            math.atan2(bounds.source_height_m[0], rho_high),
            math.atan2(bounds.detector_height_m[0], rho_high),
        ),
        (
            math.log(bounds.energy_keV[1]),
            math.log(rho_high),
            math.atan2(bounds.source_height_m[1], rho_low),
            math.atan2(bounds.detector_height_m[1], rho_low),
        ),
    )


# Primitive-polynomial parameters for Sobol dimensions 2..4.  Dimension one
# has the van-der-Corput direction numbers.  These are the Bratley-Fox/Joe-Kuo
# initial values and are sufficient for this fixed four-dimensional campaign.
_SOBOL_PARAMS = (
    None,
    (1, 0, (1,)),
    (2, 1, (1, 3)),
    (3, 1, (1, 3, 1)),
)
_SOBOL_BITS = 52
_SOBOL_SCALE = float(1 << _SOBOL_BITS)


def _splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return value ^ (value >> 31)


def _direction_numbers(dim: int) -> list[list[int]]:
    if dim < 1 or dim > 4:
        raise ValueError("the dependency-free Sobol generator supports 1..4 dimensions")
    directions = [[0] * (_SOBOL_BITS + 1) for _ in range(dim)]
    for j in range(1, _SOBOL_BITS + 1):
        directions[0][j] = 1 << (_SOBOL_BITS - j)
    for d in range(1, dim):
        s, a, initial = _SOBOL_PARAMS[d]  # type: ignore[misc]
        for j in range(1, s + 1):
            directions[d][j] = initial[j - 1] << (_SOBOL_BITS - j)
        for j in range(s + 1, _SOBOL_BITS + 1):
            value = directions[d][j - s] ^ (directions[d][j - s] >> s)
            for k in range(1, s):
                if (a >> (s - 1 - k)) & 1:
                    value ^= directions[d][j - k]
            directions[d][j] = value
    return directions


def scrambled_sobol(count: int, dim: int = 4, seed: int = 0) -> list[tuple[float, ...]]:
    """Generate a reproducible digitally shifted Sobol sequence.

    A digital XOR shift derived through SplitMix64 supplies the scramble.  The
    implementation is local so a design is bit-for-bit stable whether or not
    SciPy happens to be installed.
    """

    if count < 0:
        raise ValueError("count must be non-negative")
    directions = _direction_numbers(dim)
    mask = (1 << _SOBOL_BITS) - 1
    shifts = [(_splitmix64(seed + d) & mask) for d in range(dim)]
    state = [0] * dim
    result: list[tuple[float, ...]] = []
    for index in range(count):
        if index:
            # Position of the rightmost zero bit of index-1, one-based.
            c = ((index & -index).bit_length() - 1) + 1
            for d in range(dim):
                state[d] ^= directions[d][c]
        result.append(tuple(((state[d] ^ shifts[d]) & mask) / _SOBOL_SCALE for d in range(dim)))
    return result


def _lerp_log(pair: tuple[float, float], u: float) -> float:
    return math.exp(math.log(pair[0]) + u * math.log(pair[1] / pair[0]))


def _lerp(pair: tuple[float, float], u: float) -> float:
    return pair[0] + u * (pair[1] - pair[0])


def generate_candidate_pool(
    count: int = 10_000,
    *,
    bounds: DesignBounds = DesignBounds(),
    seed: int = 0xCEE10,
    mode: str = "field",
) -> tuple[PhysicalNode, ...]:
    """Return exactly ``count`` nodes sampled in the valid transformed domain.

    Sobol coordinates map first to ``(log E, log rho, alpha_src,
    alpha_det)``.  The angular ranges are conditioned on rho so both heights
    are in bounds by construction; only the coupled distance constraint needs
    rejection.  This avoids the table-first behavior of sampling a rectangular
    physical box and merely transforming it afterward.
    """

    if count <= 0:
        raise ValueError("count must be positive")
    bounds.validate()
    if mode not in {"air", "field"}:
        raise ValueError("candidate generation supports air or field geometry")
    accepted: list[PhysicalNode] = []
    # Generate by index so rejection never makes the result batch-size dependent.
    stream_count = max(count, 1024)
    while len(accepted) < count:
        points = scrambled_sobol(stream_count, seed=seed)
        accepted.clear()
        for u0, u1, u2, u3 in points:
            energy = _lerp_log(bounds.energy_keV, u0)
            if mode == "air":
                node = PhysicalNode(
                    energy, _lerp_log(bounds.distance_m, u1), 0.0, 0.0
                )
                accepted.append(node)
                if len(accepted) == count:
                    break
                continue

            rho = _lerp_log((bounds.rho_min_m, bounds.distance_m[1]), u1)
            src_angle = _lerp(
                (
                    math.atan2(bounds.source_height_m[0], rho),
                    math.atan2(bounds.source_height_m[1], rho),
                ),
                u2,
            )
            det_angle = _lerp(
                (
                    math.atan2(bounds.detector_height_m[0], rho),
                    math.atan2(bounds.detector_height_m[1], rho),
                ),
                u3,
            )
            h_src = rho * math.tan(src_angle)
            h_det = rho * math.tan(det_angle)
            distance = math.hypot(rho, h_src - h_det)
            if bounds.distance_m[0] <= distance <= bounds.distance_m[1]:
                node = PhysicalNode(energy, distance, h_src, h_det)
                accepted.append(node)
                if len(accepted) == count:
                    break
        stream_count *= 2
        if stream_count > (1 << 24) and len(accepted) < count:
            raise RuntimeError("bounds have too little valid geometry volume")
    return tuple(accepted)


def _node_key(node: PhysicalNode) -> tuple[str, str, str, str]:
    return tuple(
        v.hex()
        for v in (
            node.energy_keV,
            node.distance_m,
            node.source_height_m,
            node.detector_height_m,
        )
    )  # type: ignore[return-value]


def _normalizers(nodes: Sequence[PhysicalNode]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    q = [geometry_transform(node) for node in nodes]
    lo = tuple(min(row[j] for row in q) for j in range(4))
    hi = tuple(max(row[j] for row in q) for j in range(4))
    span = tuple(max(hi[j] - lo[j], 1.0e-15) for j in range(4))
    return lo, span


def _normalized_q(
    node: PhysicalNode, lo: Sequence[float], span: Sequence[float]
) -> tuple[float, float, float, float]:
    q = geometry_transform(node)
    return tuple((q[j] - lo[j]) / span[j] for j in range(4))  # type: ignore[return-value]


def _distance2(a: Sequence[float], b: Sequence[float]) -> float:
    return sum((x - y) * (x - y) for x, y in zip(a, b))


def select_disjoint_maximin(
    candidates: Sequence[PhysicalNode],
    *,
    sizes: Mapping[str, int] | None = None,
    forced: Mapping[str, Iterable[PhysicalNode]] | None = None,
    q_domain: tuple[Sequence[float], Sequence[float]] | None = None,
) -> DesignSelection:
    """Select disjoint train/validation/final sets by greedy maximin distance.

    Forced nodes count toward their role's size and need not be members of the
    candidate pool.  Distances for later roles include all earlier roles, which
    prevents a validation or final point from hugging training data.
    """

    sizes = dict(sizes or {"train": 256, "validation": 64, "final": 64})
    forced = {key: tuple(value) for key, value in (forced or {}).items()}
    roles = ("train", "validation", "final")
    if not candidates:
        raise ValueError("candidate pool is empty")
    if any(sizes.get(role, 0) < 0 for role in roles):
        raise ValueError("role sizes must be non-negative")
    all_forced = [node for role in roles for node in forced.get(role, ())]
    for node in [*candidates, *all_forced]:
        node.validate()
    keys = [_node_key(node) for node in candidates]
    if len(set(keys)) != len(keys):
        raise ValueError("candidate pool contains duplicate coordinates")
    forced_keys = [_node_key(node) for node in all_forced]
    if len(set(forced_keys)) != len(forced_keys):
        raise ValueError("forced nodes overlap across roles")
    required_candidates = sum(sizes.get(role, 0) for role in roles) - len(all_forced)
    if required_candidates > len(candidates):
        raise ValueError("not enough candidates for requested design")
    for role in roles:
        if len(forced.get(role, ())) > sizes.get(role, 0):
            raise ValueError(f"forced {role} nodes exceed requested size")

    if q_domain is None:
        lo, span = _normalizers([*candidates, *all_forced])
    else:
        lo = tuple(map(float, q_domain[0]))
        high = tuple(map(float, q_domain[1]))
        if len(lo) != 4 or len(high) != 4 or any(b <= a for a, b in zip(lo, high)):
            raise ValueError("q_domain must contain four nonempty bounds")
        span = tuple(b - a for a, b in zip(lo, high))
    candidate_q = [_normalized_q(node, lo, span) for node in candidates]
    selected_q = [_normalized_q(node, lo, span) for node in all_forced]
    unavailable = set(forced_keys)
    available = [i for i, key in enumerate(keys) if key not in unavailable]
    if required_candidates > len(available):
        raise ValueError("not enough non-forced candidates for requested design")
    result: dict[str, list[PhysicalNode]] = {
        role: list(forced.get(role, ())) for role in roles
    }

    if selected_q:
        min_d2 = [min(_distance2(q, chosen) for chosen in selected_q) for q in candidate_q]
    else:
        center = (0.5, 0.5, 0.5, 0.5)
        min_d2 = [_distance2(q, center) for q in candidate_q]

    for role in roles:
        needed = sizes.get(role, 0) - len(result[role])
        for _ in range(needed):
            # max() preserves the first index on a tie through the negative-key
            # secondary ordering, making selection stable across Python versions.
            best = max(available, key=lambda i: (min_d2[i], -i))
            available.remove(best)
            node = candidates[best]
            result[role].append(node)
            chosen_q = candidate_q[best]
            selected_q.append(chosen_q)
            for i in available:
                min_d2[i] = min(min_d2[i], _distance2(candidate_q[i], chosen_q))

    return DesignSelection(
        tuple(result["train"]),
        tuple(result["validation"]),
        tuple(result["final"]),
    )
