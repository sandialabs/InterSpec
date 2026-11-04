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

import json
from dataclasses import asdict
from pathlib import Path
import tempfile
import unittest

from tools.full_spectrum.env_scatter_fit.campaign import (
    PhysicsConfig,
    canonical_fit_config_json,
    config_sha256,
)
from tools.full_spectrum.env_scatter_fit.registered_e0 import (
    GENERATOR_SCHEMA,
    generate_registered_e0_report,
    validate_registered_e0_report,
)
from tools.full_spectrum.env_scatter_fit.tests.test_publication import (
    REGISTERED_CANONICAL_CONFIG_JSON,
    detailed_reports,
    digest,
)


ROOT = Path(__file__).resolve().parents[4]
EXECUTABLE = ROOT / "build-rel" / "tools" / "full_spectrum" / "env_registered_e0_gate"


def _campaign(canonical_json, mode, material, physics):
    return {
        "config": {
            "canonical_json": canonical_json,
            "sha256": config_sha256(canonical_json),
            "mode": mode,
            "material": material,
            "physics": asdict(physics),
        },
    }


class RegisteredE0EvidenceTests(unittest.TestCase):
    def test_raw_joint_maps_are_the_metric_authority(self):
        report = detailed_reports(
            digest("config"), digest("archive"), digest("model")
        )["registered_e0_gate"]
        self.assertTrue(validate_registered_e0_report(report))

        report["nodes"][7]["joint_cell_masses"] = [0.0, 1.0] + [0.0] * 94
        with self.assertRaisesRegex(ValueError, "independently recomputed raw nodes"):
            validate_registered_e0_report(report)

    def test_integrated_generator_records_real_subprocess_invocation(self):
        campaign = _campaign(
            REGISTERED_CANONICAL_CONFIG_JSON,
            "field",
            "soil",
            PhysicsConfig(),
        )
        config_hash = campaign["config"]["sha256"]
        raw = detailed_reports(config_hash, digest("archive"), digest("model"))[
            "registered_e0_gate"
        ]
        raw.pop("generator")

        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "registered-e0-fixture"
            output = Path(temporary) / "registered.json"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "import json, pathlib, sys\n"
                f"report = json.loads({json.dumps(raw)!r})\n"
                "path = pathlib.Path(sys.argv[sys.argv.index('--json-report') + 1])\n"
                "path.write_text(json.dumps(report), encoding='utf-8')\n"
                "print('synthetic registered-E0 subprocess completed')\n",
                encoding="utf-8",
            )
            executable.chmod(0o755)

            report = generate_registered_e0_report(
                campaign,
                executable,
                output,
                batch_histories=1_000_000,
                max_histories=128_000_000,
            )

            self.assertTrue(report["passed"])
            generator = report["generator"]
            self.assertEqual(generator["schema"], GENERATOR_SCHEMA)
            self.assertEqual(generator["argv"][0], str(executable.resolve()))
            argv = generator["argv"]
            self.assertEqual(argv[argv.index("--config-sha256") + 1], config_hash)
            self.assertEqual(generator["input_parameters"]["config_sha256"], config_hash)
            self.assertEqual(
                generator["input_parameters"]["gate_configuration"],
                raw["gate_configuration"],
            )
            self.assertEqual(json.loads(output.read_text(encoding="utf-8")), report)

    def test_built_cpp_generator_runs_transport_and_python_recomputes_report(self):
        if not EXECUTABLE.is_file():
            self.skipTest(f"env_registered_e0_gate has not been built at {EXECUTABLE}")
        physics = PhysicsConfig()
        canonical = canonical_fit_config_json("air", "soil", physics)
        campaign = _campaign(canonical, "air", "soil", physics)
        with tempfile.TemporaryDirectory() as temporary:
            report = generate_registered_e0_report(
                campaign,
                EXECUTABLE,
                Path(temporary) / "real-low-history.json",
                batch_histories=32,
                max_histories=32,
                threads=4,
            )
        self.assertFalse(report["passed"])
        self.assertEqual(len(report["nodes"]), 13)
        self.assertFalse(report["all_node_stopping_targets_met"])
        self.assertFalse(validate_registered_e0_report(report))


if __name__ == "__main__":
    unittest.main()
