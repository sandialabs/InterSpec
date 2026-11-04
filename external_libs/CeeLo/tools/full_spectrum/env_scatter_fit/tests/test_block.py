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

import hashlib
import io
import json
import math
import random
import struct
import unittest

from tools.full_spectrum.env_scatter_fit.axes import (
    annihilation_registered_axis,
    brems_registered_axis,
    primary_registered_axis,
    uniform_mu_edges,
)
from tools.full_spectrum.env_scatter_fit.block import BlockFormatError, HEADER, MAGIC, read_fit_block
from tools.full_spectrum.env_scatter_fit.archive import _merged_node_status
from tools.full_spectrum.env_scatter_fit.campaign import replica_seed
from tools.full_spectrum.env_scatter_fit.contracts import (
    STOPPING_VECTOR_SCHEMA,
    stopping_vector_definition_json,
)
from tools.full_spectrum.env_scatter_fit.design import PhysicalNode, geometry_transform


def string(value):
    raw = value.encode()
    return struct.pack("<Q", len(raw)) + raw


def vector(code, values):
    values = tuple(values)
    return struct.pack("<Q", len(values)) + b"".join(struct.pack("<" + code, x) for x in values)


def moment_values(means, histories=262_144, *, variances=None, corrupt_s2=False):
    means = list(means)
    variances = [0.0] * len(means) if variances is None else list(variances)
    if len(variances) != len(means):
        raise ValueError("moment means/variances have different lengths")
    second = [
        (histories - 1) * variance + histories * value * value
        for value, variance in zip(means, variances)
    ]
    if corrupt_s2:
        nonzero = next(index for index, value in enumerate(means) if value)
        second[nonzero] = 0.0
    return (
        vector("Q", [histories] * len(means))
        + vector("d", [histories * value for value in means])
        + vector("d", second)
    )


def fixture_block(
    *,
    block_id=1,
    mu_bins=32,
    role="train",
    scale=1.0,
    status="complete",
    capped=False,
    failure_reason="",
    transport_failures=0,
    corrupt_taxonomy=False,
    corrupt_replica=False,
    corrupt_s2=False,
    corrupt_manifest=False,
    corrupt_path=False,
    corrupt_stopping=False,
    corrupt_stopping_definition=False,
    physical_node=None,
    histories=262_144,
    energy_bins=2,
    scoring_radius_cm=None,
    calibration_variant=None,
):
    if histories <= 0 or histories % 32:
        raise ValueError("fixture histories must be positive and divisible by 32")
    if energy_bins <= 0 or energy_bins % 2:
        raise ValueError("fixture energy bins must be positive and even")
    e_cut = 10.0
    air_definition = {
        "name": "Air",
        "density_g_cm3": 0.001205,
        "composition": [],
    }
    configuration = {
        "air": air_definition,
        "axis": {
            "energy_bins": energy_bins,
            "energy_cut_keV": e_cut,
            "mu_bins": mu_bins,
            "transform": "registered-component-v1",
        },
        "block_histories": histories,
        "geometry_tier": "air",
        "material": None,
        "physics": {
            "brems_electron_floor_keV": 150.0,
            "include_annihilation": True,
            "include_brems": True,
            "include_rayleigh": True,
        },
        "schema": "ceelo-env-fit-config-v1",
    }
    if scoring_radius_cm is not None:
        configuration["scoring_radius_cm"] = float(scoring_radius_cm)
    canonical = json.dumps(
        configuration,
        sort_keys=True,
        separators=(",", ":"),
    )
    config_hash = hashlib.sha256(canonical.encode()).hexdigest()
    node = physical_node or PhysicalNode(100.0, 1.0, 0.0, 0.0)
    node.validate()
    node_id = node.stable_id(config_hash)
    mu_edges = uniform_mu_edges(mu_bins)
    axes = (
        ("primary_registered", primary_registered_axis(node.energy_keV, e_cut_keV=e_cut, n_per_block=energy_bins // 2, n_mu=mu_bins)),
        ("annihilation_511_registered", annihilation_registered_axis(e_cut_keV=e_cut, n_per_block=energy_bins // 2, n_mu=mu_bins)),
        ("brems_log_relative", brems_registered_axis(node.energy_keV, e_cut_keV=e_cut, n_energy=energy_bins, n_mu=mu_bins)),
    )
    cells_per_component = energy_bins * mu_bins
    include_annihilation_line = node.energy_keV >= 1021.9979
    line_count = 2 if include_annihilation_line else 1
    line_cells = line_count * mu_bins
    component_supported = []
    for origin in range(3):
        for order in range(3):
            component_supported.append(
                (origin == 0 and order != 0)
                or (origin == 1 and include_annihilation_line and order != 0)
                or (origin == 2 and node.energy_keV >= 150.0)
            )

    def shifted_mean(base, variance, cell, family_seed):
        pair_index = (block_id - 1) // 2
        generator = random.Random(
            (node_id ^ (pair_index << 20) ^ (family_seed << 40) ^ cell)
            & ((1 << 64) - 1)
        )
        standardized = max(-3.5, min(3.5, generator.gauss(0.0, 1.0)))
        side = 1.0 if block_id % 2 else -1.0
        return base + side * 0.5 * standardized * math.sqrt(2.0 * variance / histories)

    def replica_maps(means, variances, family_seed):
        histories_per_replica = histories // 32
        replica_scale = 0.025 if calibration_variant == "bad_replica_covariance" else 1.0
        output = [0.0] * (32 * len(means))
        for cell, (mean, variance) in enumerate(zip(means, variances)):
            if variance <= 0.0:
                for replica in range(32):
                    output[replica * len(means) + cell] = mean * histories_per_replica
                continue
            generator = random.Random(
                (node_id ^ (block_id << 16) ^ (family_seed << 40) ^ cell)
                & ((1 << 64) - 1)
            )
            samples = [generator.gauss(0.0, 1.0) for _ in range(32)]
            center = math.fsum(samples) / 32.0
            samples = [value - center for value in samples]
            sample_variance = math.fsum(value * value for value in samples) / 31.0
            amplitude = replica_scale * math.sqrt(
                variance / (histories_per_replica * sample_variance)
            )
            for replica, value in enumerate(samples):
                replica_mean = mean + amplitude * value
                if replica_mean < 0.0:
                    raise ValueError("calibration fixture generated a negative replica mean")
                output[replica * len(means) + cell] = (
                    replica_mean * histories_per_replica
                )
        return output

    if calibration_variant is None:
        component_means = [0.0] * (9 * cells_per_component)
        component_means[1 * cells_per_component] = 0.6 * scale
        component_means[2 * cells_per_component + 1] = 0.4 * scale
        component_variances = [0.0] * len(component_means)
        total_means = [0.0] * cells_per_component
        total_means[0], total_means[1] = 0.6 * scale, 0.4 * scale
        total_variances = [0.0] * len(total_means)
        line_component_means = [0.0] * (9 * line_cells)
        line_component_means[0] = 0.2 * scale
        line_component_variances = [0.0] * len(line_component_means)
        line_total_means = [0.0] * line_cells
        line_total_means[0] = 0.2 * scale
        line_total_variances = [0.0] * len(line_total_means)
        histories_per_replica = histories // 32
        replica_component = [
            value * histories_per_replica for _ in range(32) for value in component_means
        ]
        replica_total = [
            value * histories_per_replica for _ in range(32) for value in total_means
        ]
        replica_line_component = [
            value * histories_per_replica
            for _ in range(32)
            for value in line_component_means
        ]
        replica_line_total = [
            value * histories_per_replica for _ in range(32) for value in line_total_means
        ]
    else:
        true_component_variances = [0.0] * (9 * cells_per_component)
        component_means = [0.0] * len(true_component_variances)
        component_bases = [0.0] * len(true_component_variances)
        for component in range(9):
            if not component_supported[component] or (
                calibration_variant == "empty_supported_component" and component == 8
            ):
                continue
            for cell in range(cells_per_component):
                index = component * cells_per_component + cell
                shape = 0.8 + 0.4 * ((cell % mu_bins) + 0.5) / mu_bins
                base = scale * (0.8 + 0.08 * component) * shape
                variance = scale * scale * (0.6 + 0.03 * component)
                component_bases[index] = base
                true_component_variances[index] = variance
                component_means[index] = shifted_mean(base, variance, index, 1)
        total_means = [0.0] * cells_per_component
        true_total_variances = [0.0] * cells_per_component
        for component in range(9):
            for cell in range(cells_per_component):
                index = component * cells_per_component + cell
                total_means[cell] += component_means[index]
                true_total_variances[cell] += true_component_variances[index]
        component_variances = list(true_component_variances)
        total_variances = list(true_total_variances)
        if calibration_variant == "bad_component_s2":
            component_variances = [value * 0.001 for value in component_variances]
        if calibration_variant == "bad_total_cell_s2":
            total_variances = [value * 0.001 for value in total_variances]

        line_component_means = [0.0] * (9 * line_cells)
        line_component_bases = [0.0] * len(line_component_means)
        true_line_component_variances = [0.0] * len(line_component_means)
        # The primary elastic atom is first scored at the Rayleigh interaction
        # that redirects it toward the detector, hence primary order 1.
        active_line_components = [(1, 0)]
        if include_annihilation_line:
            active_line_components.append((3, 1))
        for component, line in active_line_components:
            for mu_index in range(mu_bins):
                cell = component * line_cells + line * mu_bins + mu_index
                base = scale * (0.45 + 0.15 * line) * (
                    0.85 + 0.3 * (mu_index + 0.5) / mu_bins
                )
                variance = scale * scale * (0.35 + 0.1 * line)
                line_component_bases[cell] = base
                true_line_component_variances[cell] = variance
                line_component_means[cell] = shifted_mean(base, variance, cell, 2)
        line_total_means = [0.0] * line_cells
        true_line_total_variances = [0.0] * line_cells
        for component in range(9):
            for cell in range(line_cells):
                index = component * line_cells + cell
                line_total_means[cell] += line_component_means[index]
                true_line_total_variances[cell] += true_line_component_variances[index]
        line_component_variances = list(true_line_component_variances)
        line_total_variances = list(true_line_total_variances)
        if calibration_variant == "bad_line_s2":
            line_component_variances = [
                value * 0.001 for value in line_component_variances
            ]
            line_total_variances = [value * 0.001 for value in line_total_variances]

        replica_component = replica_maps(
            component_means, true_component_variances, 11
        )
        replica_line_component = replica_maps(
            line_component_means, true_line_component_variances, 12
        )
        histories_per_replica = histories // 32
        replica_total = [0.0] * (32 * cells_per_component)
        replica_line_total = [0.0] * (32 * line_cells)
        for replica in range(32):
            for component in range(9):
                for cell in range(cells_per_component):
                    replica_total[replica * cells_per_component + cell] += (
                        replica_component[
                            replica * len(component_means)
                            + component * cells_per_component
                            + cell
                        ]
                    )
                for cell in range(line_cells):
                    replica_line_total[replica * line_cells + cell] += (
                        replica_line_component[
                            replica * len(line_component_means)
                            + component * line_cells
                            + cell
                        ]
                    )
    if corrupt_replica:
        replica_component[0] += 0.1
    payload = bytearray()
    payload += struct.pack("<QQ", node_id, block_id)
    payload += string(canonical) + string(config_hash)
    for value in (
        "Air",
        json.dumps(air_definition, sort_keys=True, separators=(",", ":")),
        "air",
        "mt19937_64/sha256-tuple-stream-v1",
        (
            "incorrect-estimator-contract"
            if corrupt_manifest
            else "next-event;rayleigh=on;brems=on;annihilation=on"
        ),
        "energy=keV,length=cm,fluence=cm^-2",
        "bin-integrated scalar incident fluence per source photon",
        "raw S1/S2; s_hat divisor=T(E0)/(4*pi*d^2)",
    ):
        payload += string(value)
    payload += struct.pack("<Q", 3)
    for name, axis in axes:
        payload += string(name) + struct.pack("<I", 1)
        payload += vector("d", axis.registered_edges)
        payload += vector("d", axis.physical_energy_edges_keV)
    payload += struct.pack("<I", 0) + vector("d", mu_edges)
    payload += struct.pack("<Q", line_count)
    payload += string("elastic_E0") + struct.pack("<dII", node.energy_keV, 0, 3)
    if include_annihilation_line:
        payload += string("annihilation_511") + struct.pack(
            "<dII", 510.998950, 1, 3
        )
    payload += struct.pack("<Q", 9)
    for origin, origin_name in enumerate(("primary", "annihilation", "brems")):
        for order, order_name in enumerate(("order0", "order1", "order2plus")):
            name = f"{origin_name}_{order_name}"
            if corrupt_taxonomy and origin == 2 and order == 2:
                name = "unknown_component"
            supported = component_supported[origin * 3 + order]
            payload += string(name) + struct.pack("<III", origin, order, origin)
            payload += vector("B", [int(supported)] * energy_bins)
    payload += struct.pack(
        "<5d",
        node.energy_keV,
        100.0 * node.distance_m,
        100.0 * node.source_height_m,
        100.0 * node.detector_height_m,
        0.0,
    )
    payload += struct.pack("<4d", *geometry_transform(node))
    payload += string(role) + string(status)
    direct_fluence = 1.0 / (4.0 * math.pi * (100.0 * node.distance_m) ** 2)
    payload += struct.pack(
        "<ddQdddd", 1.0, direct_fluence, histories, 1.0, 2.0, 0.0025, 0.01
    )
    payload += struct.pack("<I", 32)
    histories_per_replica = histories // 32
    for replica in range(32):
        payload += struct.pack(
            "<IQQQ",
            replica,
            (block_id - 1) * histories + replica * histories_per_replica,
            histories_per_replica,
            replica_seed(config_hash, node_id, block_id, replica),
        )
    payload += moment_values(
        component_means,
        histories,
        variances=component_variances,
        corrupt_s2=corrupt_s2,
    )
    payload += moment_values(total_means, histories, variances=total_variances)
    payload += moment_values(
        line_component_means, histories, variances=line_component_variances
    )
    payload += moment_values(
        line_total_means, histories, variances=line_total_variances
    )
    payload += vector("f", replica_component)
    payload += vector("f", replica_total)
    payload += vector("f", replica_line_component)
    payload += vector("f", replica_line_total)

    scalar_mean = math.fsum(total_means) + math.fsum(line_total_means)
    if calibration_variant is None:
        scalar_variance = 0.0
        stopping_means = [0.0] * 16
        stopping_means[0] = 0.0 if corrupt_stopping else scalar_mean
        stopping_covariance = [[0.0] * 16 for _ in range(16)]
    else:
        scalar_variance = math.fsum(true_total_variances) + math.fsum(
            true_line_total_variances
        )
        # Recover the pair center without duplicating the large cell arrays.
        pair_scalar_mean = math.fsum(component_bases) + math.fsum(
            line_component_bases
        )
        side = 1.0 if block_id % 2 else -1.0
        scalar_pair_difference = 2.0 * side * (scalar_mean - pair_scalar_mean)
        pair_index = (block_id - 1) // 2
        band_generator = random.Random(
            (node_id ^ (pair_index << 24) ^ 0xBADD5EED) & ((1 << 64) - 1)
        )
        band_modes = [band_generator.gauss(0.0, 1.0) for _ in range(12)]
        band_center = math.fsum(band_modes) / 12.0
        band_scale = math.sqrt(2.0 * scalar_variance / (12.0 * histories))
        band_differences = [
            scalar_pair_difference / 12.0
            + band_scale * (value - band_center)
            for value in band_modes
        ]
        stopping_means = [
            pair_scalar_mean / 12.0 + side * 0.5 * difference
            for difference in band_differences
        ]
        for order in range(4):
            variance = scalar_variance * (0.04 + 0.01 * order)
            generator = random.Random(
                (node_id ^ (pair_index << 24) ^ (0x1E600000 + order))
                & ((1 << 64) - 1)
            )
            difference = max(-3.5, min(3.5, generator.gauss(0.0, 1.0))) * math.sqrt(
                2.0 * variance / histories
            )
            stopping_means.append(
                0.02 * (order + 1) * pair_scalar_mean + side * 0.5 * difference
            )
        stopping_covariance = [[0.0] * 16 for _ in range(16)]
        for index in range(12):
            stopping_covariance[index][index] = scalar_variance / 12.0
        for order in range(4):
            stopping_covariance[12 + order][12 + order] = scalar_variance * (
                0.04 + 0.01 * order
            )
        if calibration_variant in {"bad_s2", "bad_stopping_covariance"}:
            scalar_variance *= 0.001
            for row in range(16):
                for column in range(16):
                    stopping_covariance[row][column] *= 0.001

    stopping_s1 = [histories * value for value in stopping_means]
    stopping_outer = [0.0] * 256
    for row in range(16):
        for column in range(16):
            stopping_outer[row * 16 + column] = (
                (histories - 1) * stopping_covariance[row][column]
                + histories * stopping_means[row] * stopping_means[column]
            )
    stopping_definition = stopping_vector_definition_json(
        0, "primary_registered", energy_bins
    )
    if corrupt_stopping_definition:
        stopping_definition = "{}"
    payload += string(stopping_definition)
    payload += struct.pack("<Q", histories) + struct.pack("<16d", *stopping_s1)
    payload += struct.pack("<256d", *stopping_outer)
    channels = (
        "air_only",
        "one_slab_bounce",
        "multi_same_slab",
        "multi_distinct_slab",
    )
    payload += struct.pack("<Q", len(channels))
    for channel_index, name in enumerate(channels):
        payload += string(name)
        channel_continuum = (
            list(total_means)
            if channel_index == 0
            else [0.0] * cells_per_component
        )
        if corrupt_path and channel_index == 0:
            channel_continuum[0] += 0.1
        payload += moment_values(channel_continuum, histories)
        payload += moment_values(
            line_total_means if channel_index == 0 else [0.0] * line_cells,
            histories,
        )
    interactions = ("photoelectric", "compton", "rayleigh", "pair_production", "brems_emission")
    payload += struct.pack("<Q", len(interactions))
    for name in interactions:
        payload += string(name)
    payload += vector("Q", [0] * len(interactions))
    scalar_second = (
        (histories - 1) * scalar_variance + histories * scalar_mean * scalar_mean
    )
    payload += struct.pack("<Qdd", histories, histories * scalar_mean, scalar_second)
    payload += struct.pack("<Qdd", histories, 0.0, 0.0)
    payload += struct.pack(
        "<QdB", transport_failures, 1.2 * scale, int(capped)
    ) + string(failure_reason)
    digest = hashlib.sha256(payload).digest()
    return HEADER.pack(MAGIC, 1, 0, len(payload), digest) + payload


class BlockTests(unittest.TestCase):
    def test_capped_reason_remains_usable_but_transport_failure_wins(self):
        self.assertEqual(
            _merged_node_status(
                "complete",
                "capped",
                transport_failures=0,
                capped=True,
                failure_reason="maximum vertex cap reached",
            ),
            "capped",
        )
        self.assertEqual(
            _merged_node_status(
                "capped",
                "complete",
                transport_failures=1,
                capped=True,
                failure_reason="transport failure",
            ),
            "failed",
        )

    def test_decode_and_validate_reference_fixture(self):
        block = read_fit_block(io.BytesIO(fixture_block()))
        self.assertEqual(block.block_id, 1)
        self.assertEqual(block.node["histories"], 262_144)
        self.assertEqual(len(block.replicas), 32)
        self.assertEqual(len(block.moments["continuum_components"]["S1"]), 9 * 2 * 32)
        self.assertEqual(block.diagnostics["scattered_total"]["N"], 262_144)
        definition = json.loads(block.stopping["definition_json"])
        self.assertEqual(definition["version"], STOPPING_VECTOR_SCHEMA)
        self.assertEqual(len(definition["entries"]), 16)
        self.assertEqual(definition["entries"][12]["legendre_order"], 1)

    def test_semantic_corruption_is_rejected_after_valid_checksum(self):
        with self.assertRaisesRegex(BlockFormatError, "taxonomy"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_taxonomy=True)))
        with self.assertRaisesRegex(BlockFormatError, "replica continuum_components"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_replica=True)))
        with self.assertRaisesRegex(BlockFormatError, "S2"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_s2=True)))
        with self.assertRaisesRegex(BlockFormatError, "canonical physics configuration"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_manifest=True)))
        with self.assertRaisesRegex(BlockFormatError, "path diagnostic continuum"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_path=True)))
        with self.assertRaisesRegex(BlockFormatError, "stopping-band"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_stopping=True)))
        with self.assertRaisesRegex(BlockFormatError, "stopping-vector definition"):
            read_fit_block(io.BytesIO(fixture_block(corrupt_stopping_definition=True)))

    def test_checksum_and_trailing_bytes_are_rejected(self):
        damaged = bytearray(fixture_block())
        damaged[-1] ^= 1
        with self.assertRaisesRegex(BlockFormatError, "SHA-256"):
            read_fit_block(io.BytesIO(damaged))
        with self.assertRaisesRegex(BlockFormatError, "follow"):
            read_fit_block(io.BytesIO(fixture_block() + b"x"))

    def test_formal_archive_rejects_noncanonical_block_size(self):
        with self.assertRaisesRegex(BlockFormatError, "block size"):
            read_fit_block(io.BytesIO(fixture_block(histories=32)))


if __name__ == "__main__":
    unittest.main()
