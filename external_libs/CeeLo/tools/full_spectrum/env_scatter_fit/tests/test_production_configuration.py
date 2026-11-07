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

from copy import deepcopy
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

from tools.full_spectrum.env_scatter_fit.campaign import (
    PhysicsConfig,
    canonical_fit_config_json,
    config_sha256,
    default_forced_anchors,
)
from tools.full_spectrum.env_scatter_fit.design import (
    DesignBounds,
    generate_candidate_pool,
    geometry_transform,
    select_disjoint_maximin,
    transformed_domain,
)
from tools.full_spectrum.env_scatter_fit.publication import (
    production_configuration_gate,
    require_production_configuration,
)
from tools.full_spectrum.env_scatter_fit.workflow import promote_candidate


class Dataset:
    def __init__(self, values=None, *, shape=None):
        self.values = values
        if shape is not None:
            self.shape = tuple(shape)
        elif values and isinstance(values[0], list):
            self.shape = (len(values), len(values[0]))
        else:
            self.shape = (len(values or ()),)

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, key):
        if key is Ellipsis:
            return self.values
        if isinstance(key, tuple) and key == (slice(None), 0):
            return [row[0] for row in self.values]
        return self.values[key]


class Group(dict):
    def __init__(self, *args, attrs=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.attrs = dict(attrs or {})


class Root(Group):
    def __getitem__(self, key):
        if "/" not in key:
            return super().__getitem__(key)
        value = self
        for part in key.split("/"):
            value = dict.__getitem__(value, part)
        return value


def component_descriptors():
    return [
        {
            "name": f"{origin_name}_{order_name}",
            "origin": origin,
            "post_creation_order": order,
            "energy_axis": origin,
        }
        for origin, origin_name in enumerate(("primary", "annihilation", "brems"))
        for order, order_name in enumerate(("order0", "order1", "order2plus"))
    ]


def campaign_and_root(physics=PhysicsConfig()):
    canonical = canonical_fit_config_json("field", "soil", physics)
    config_hash = config_sha256(canonical)
    bounds = DesignBounds()
    sobol_seed = 0xCEE10
    candidates = generate_candidate_pool(10_000, bounds=bounds, seed=sobol_seed)

    def node_record(node, role):
        return {
            **node.as_dict(),
            "q": list(geometry_transform(node)),
            "role": role,
            "node_id": node.stable_id(config_hash),
            "node_id_hex": node.stable_id_hex(config_hash),
        }

    forced = default_forced_anchors(bounds, mode="field")
    q_domain = transformed_domain(bounds, mode="field")
    selection = select_disjoint_maximin(
        candidates,
        sizes={"train": 256, "validation": 64, "final": 64},
        forced={"train": forced},
        q_domain=q_domain,
    )
    train = [node_record(node, "train") for node in selection.train]
    validation = [node_record(node, "validation") for node in selection.validation]
    final = [node_record(node, "final") for node in selection.final]
    selected_roles = {
        item["node_id_hex"]: role
        for role, records in (
            ("train", train),
            ("validation", validation),
            ("final", final),
        )
        for item in records
    }
    candidate_pool = [
        {
            "candidate_index": index,
            **node_record(node, "candidate"),
            "selected_role": selected_roles.get(node.stable_id_hex(config_hash)),
        }
        for index, node in enumerate(candidates)
    ]
    allocation = {
        "block_histories": 262_144,
        "logical_replicas": 32,
        "initial_blocks_per_node": 4,
        "total_relative_se_target": 0.0025,
        "material_contribution_fraction": 0.01,
        "contribution_absolute_se_fraction": 0.001,
        "sentinel_absolute_se_fraction": 0.0025,
        "training_history_cap": 64_000_000,
        "validation_history_cap": 128_000_000,
        "acquisition_batch_size": 64,
        "maximum_training_nodes": 2_048,
        "maximum_final_generations": 3,
    }
    all_nodes = [*train, *validation, *final]
    selection_provenance = {
        "method": "greedy-maximin-q-v1",
        "forced_train": [node_record(node, "train") for node in forced],
        "initial_role_node_ids": {
            "train": [item["node_id_hex"] for item in train],
            "validation": [item["node_id_hex"] for item in validation],
            "final": [item["node_id_hex"] for item in final],
        },
    }
    selection_provenance["canonical_sha256"] = hashlib.sha256(
        json.dumps(
            selection_provenance,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    campaign = {
        "mu_resolution_decision": None,
        "config": {
            "canonical_json": canonical,
            "sha256": config_hash,
            "mode": "field",
            "material": "soil",
            "physics": asdict(physics),
        },
        "design": {
            "bounds": asdict(bounds),
            "sobol_seed": sobol_seed,
            "q_domain": dict(
                zip(("min", "max"), transformed_domain(bounds, mode="field"))
            ),
            "candidate_count": 10_000,
            "candidate_pool": candidate_pool,
            "forced_train_count": len(forced),
            "selection": selection_provenance,
            "nodes": {
                "train": train,
                "validation": validation,
                "final": final,
            },
        },
        "allocation": allocation,
        "campaign_state": {
            "training_node_ids": [item["node_id_hex"] for item in train],
            "acquired_node_ids": [],
            "retired_final_node_ids": [],
            "acquisition_history": [],
            "final_generation": 0,
        },
        "worker_groups": [
            {
                "node_id_hex": item["node_id_hex"],
                "first_block_id": 1,
                "block_count": 4,
            }
            for item in all_nodes
        ],
    }
    registered = [index / 384.0 for index in range(385)]
    physical = [[10.0 + index for index in range(385)]]
    axis_names = (
        "primary_registered",
        "annihilation_511_registered",
        "brems_log_relative",
    )
    energy = Group(
        {
            str(index): Group(
                {
                    "registered_edges": Dataset(registered),
                    "physical_edges_keV": Dataset(physical),
                },
                attrs={"name": name, "transform_version": 1},
            )
            for index, name in enumerate(axis_names)
        }
    )
    mu_edges = [
        -1.0 + 2.0 * index / physics.mu_bins
        for index in range(physics.mu_bins + 1)
    ]
    line_definitions = [
        {"name": "elastic_E0", "origin": 0, "post_creation_order": 3},
        {"name": "annihilation_511", "origin": 1, "post_creation_order": 3},
    ]
    estimator_flags = (
        "next-event;rayleigh="
        + ("on" if physics.include_rayleigh else "off")
        + ";brems="
        + ("on" if physics.include_brems else "off")
        + ";annihilation="
        + ("on" if physics.include_annihilation else "off")
    )
    component_cells = 9 * 384 * physics.mu_bins
    total_cells = 384 * physics.mu_bins
    root = Root(
        {
            "manifest": Group(
                attrs={
                    "config_sha256": config_hash,
                    "canonical_config_json": canonical,
                    "estimator_flags": estimator_flags,
                }
            ),
            "axes": Group(
                {"energy": energy, "mu_edges": Dataset(mu_edges)},
                attrs={
                    "total_energy_axis": 0,
                    "line_definitions_json": json.dumps(line_definitions),
                },
            ),
            "nodes": Group({"id": Dataset([1])}),
            "components": Group(
                attrs={"descriptors_json": json.dumps(component_descriptors())}
            ),
            "continuum": Group(
                {
                    "components": Group(
                        {"N": Dataset(shape=(1, component_cells))}
                    ),
                    "total": Group({"N": Dataset(shape=(1, total_cells))}),
                }
            ),
            "lines": Group(
                {
                    "energy_keV": Dataset([[100.0, 510.998950]]),
                    "support_mask": Dataset([[1, 0]]),
                }
            ),
        },
        attrs={"schema": "CEEFIT-HDF5-v1"},
    )
    return campaign, root


class ProductionConfigurationTests(unittest.TestCase):
    def test_canonical_64_bin_archive_passes_validated_baseline(self):
        campaign, root = campaign_and_root(PhysicsConfig(mu_bins=64))
        report = production_configuration_gate(root, campaign)
        self.assertTrue(report["passed"], report["issues"])
        self.assertEqual(report["mu_policy"], "64-bin-pilot-or-production")
        self.assertEqual(report["candidate_count"], 10_000)
        self.assertEqual(
            report["initial_design_counts"],
            {"train": 256, "validation": 64, "final": 64},
        )
        self.assertEqual(report["initial_blocks_per_node"], 4)
        self.assertEqual(campaign["allocation"]["maximum_final_generations"], 3)
        self.assertEqual(report["scoring_radius_cm"], 4.0)
        self.assertEqual(report["scoring_radius_policy"], "validated-4-cm-baseline-v1")

    def test_adaptive_training_growth_preserves_initial_design_eligibility(self):
        campaign, root = campaign_and_root()
        pool_records = [
            record
            for record in campaign["design"]["candidate_pool"]
            if record["selected_role"] is None
        ][:64]
        acquired = []
        for pool_record in pool_records:
            pool_record["selected_role"] = "train"
            acquired.append(
                {
                    key: value
                    for key, value in pool_record.items()
                    if key not in {"candidate_index", "selected_role"}
                }
                | {"role": "train"}
            )
        campaign["design"]["nodes"]["train"].extend(acquired)
        acquired_ids = [
            item["node_id_hex"] for item in acquired
        ]
        campaign["campaign_state"]["acquired_node_ids"] = acquired_ids
        campaign["campaign_state"]["training_node_ids"].extend(acquired_ids)
        campaign["campaign_state"]["acquisition_history"].append(
            {"node_ids": acquired_ids}
        )
        campaign["worker_groups"].extend(
            {
                "node_id_hex": item["node_id_hex"],
                "first_block_id": 1,
                "block_count": 4,
            }
            for item in acquired
        )

        report = production_configuration_gate(root, campaign)
        self.assertTrue(report["passed"], report["issues"])
        self.assertEqual(report["initial_design_counts"]["train"], 256)

    def test_retired_final_growth_preserves_initial_design_eligibility(self):
        campaign, root = campaign_and_root()
        old_final = campaign["design"]["nodes"]["final"]
        retired_ids = [item["node_id_hex"] for item in old_final]
        campaign["design"]["nodes"]["train"].extend(
            {**item, "role": "diagnostic"} for item in old_final
        )
        pool_by_id = {
            record["node_id_hex"]: record
            for record in campaign["design"]["candidate_pool"]
        }
        for node_id in retired_ids:
            pool_by_id[node_id]["selected_role"] = "diagnostic"

        new_final = []
        replacement_records = [
            record
            for record in campaign["design"]["candidate_pool"]
            if record["selected_role"] is None
        ][:64]
        for pool_record in replacement_records:
            pool_record["selected_role"] = "final"
            new_final.append(
                {
                    key: value
                    for key, value in pool_record.items()
                    if key not in {"candidate_index", "selected_role"}
                }
                | {"role": "final"}
            )
        campaign["design"]["nodes"]["final"] = new_final
        campaign["campaign_state"]["training_node_ids"].extend(retired_ids)
        campaign["campaign_state"]["retired_final_node_ids"] = retired_ids
        campaign["campaign_state"]["final_generation"] = 1
        campaign["worker_groups"].extend(
            {
                "node_id_hex": item["node_id_hex"],
                "first_block_id": 1,
                "block_count": 4,
            }
            for item in new_final
        )

        report = production_configuration_gate(root, campaign)
        self.assertTrue(report["passed"], report["issues"])
        self.assertEqual(report["initial_design_counts"]["train"], 256)

    def test_altered_campaign_dimensions_remain_diagnostic_only(self):
        campaign, root = campaign_and_root()
        cases = {
            "candidate count": lambda item: item["design"].update(
                {"candidate_count": 9_999}
            ),
            "candidate pool": lambda item: item["design"]["candidate_pool"].pop(),
            "candidate index": lambda item: item["design"]["candidate_pool"][1].update(
                {"candidate_index": 7}
            ),
            "candidate coordinate": lambda item: item["design"]["candidate_pool"][1].update(
                {"energy_keV": 100.0}
            ),
            "bounds": lambda item: item["design"]["bounds"].update(
                {"distance_m": (0.5, 40.0)}
            ),
            "initial train": lambda item: item["design"]["nodes"]["train"].pop(),
            "validation": lambda item: item["design"]["nodes"]["validation"].pop(),
            "final": lambda item: item["design"]["nodes"]["final"].pop(),
            "initial blocks": lambda item: item["allocation"].update(
                {"initial_blocks_per_node": 5}
            ),
            "allocation target": lambda item: item["allocation"].update(
                {"total_relative_se_target": 0.005}
            ),
        }
        for name, mutate in cases.items():
            with self.subTest(name=name):
                altered = deepcopy(campaign)
                mutate(altered)
                if name == "initial train":
                    altered["campaign_state"]["training_node_ids"].pop()
                report = production_configuration_gate(root, altered)
                self.assertFalse(report["passed"], report)

    def test_altered_baseline_physics_remains_diagnostic_only(self):
        baseline = PhysicsConfig()
        cases = {
            "brems_electron_floor_keV": 149.0,
            "doppler": True,
            "max_vertices": 49,
            "roulette_tau": 7.0,
            "weight_floor": 2.0e-5,
            "scoring_radius_cm": 7.25,
        }
        for name, value in cases.items():
            with self.subTest(name=name):
                campaign, root = campaign_and_root(replace(baseline, **{name: value}))
                report = production_configuration_gate(root, campaign)
                self.assertFalse(report["passed"], report)
                needle = "scoring radius" if name == "scoring_radius_cm" else name
                self.assertTrue(
                    any(needle in issue for issue in report["issues"]),
                    report["issues"],
                )

    def test_compact_archive_remains_diagnostic_only(self):
        campaign, root = campaign_and_root()
        for group in root["axes/energy"].values():
            group["registered_edges"] = Dataset([index / 8.0 for index in range(9)])
            group["physical_edges_keV"] = Dataset(
                [[10.0 + index for index in range(9)]]
            )
        root["continuum/components"]["N"] = Dataset(shape=(1, 9 * 8 * 64))
        root["continuum/total"]["N"] = Dataset(shape=(1, 8 * 64))

        report = production_configuration_gate(root, campaign)
        self.assertFalse(report["passed"])
        self.assertTrue(any("not 384 bins" in issue for issue in report["issues"]))
        with self.assertRaisesRegex(ValueError, "diagnostic-only.*384"):
            require_production_configuration(root, campaign)

    def test_disabled_validated_physics_cannot_pass(self):
        campaign, root = campaign_and_root(PhysicsConfig(include_rayleigh=False))
        report = production_configuration_gate(root, campaign)
        self.assertFalse(report["passed"])
        self.assertTrue(
            any("include_rayleigh must be enabled" in issue for issue in report["issues"])
        )

    def test_32_bins_require_identity_bound_pilot_decision(self):
        campaign, root = campaign_and_root(PhysicsConfig(mu_bins=32))
        report = production_configuration_gate(root, campaign)
        self.assertFalse(report["passed"])
        self.assertTrue(
            any("32-bin production decision is invalid" in issue for issue in report["issues"])
        )

    def test_promote_candidate_rejects_compact_archive_before_external_evidence(self):
        campaign, root = campaign_and_root()
        for group in root["axes/energy"].values():
            group["registered_edges"] = Dataset([index / 8.0 for index in range(9)])
            group["physical_edges_keV"] = Dataset(
                [[10.0 + index for index in range(9)]]
            )
        class Payload:
            accepted_for_publication = False
            model_content_sha256 = "1" * 64
            metadata = {"config_sha256": campaign["config"]["sha256"]}

            def compute_model_content_sha256(self):
                return self.model_content_sha256

        class OpenArchive:
            def __enter__(self):
                return root

            def __exit__(self, *unused):
                return False

        with tempfile.TemporaryDirectory() as temporary, patch(
            "tools.full_spectrum.env_scatter_fit.workflow.RuntimePayload.load_binary",
            return_value=Payload(),
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.load_manifest", return_value=campaign
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow._campaign_stop_status",
            return_value={"passed": True},
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow.require_hdf5",
            return_value=(SimpleNamespace(File=lambda *unused: OpenArchive()), object()),
        ), patch(
            "tools.full_spectrum.env_scatter_fit.workflow._external_detector_validation",
            side_effect=AssertionError("compact archive must fail before external evidence"),
        ):
            directory = Path(temporary)
            with self.assertRaisesRegex(ValueError, "diagnostic-only.*384"):
                promote_candidate(
                    directory / "model.candidate",
                    directory / "archive.h5",
                    directory / "campaign.json",
                    directory / "detector.json",
                    directory / "publication.json",
                    directory / "model.ceefitp",
                )


if __name__ == "__main__":
    unittest.main()
