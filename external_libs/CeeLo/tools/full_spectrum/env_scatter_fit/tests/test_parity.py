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
import json
import math
from pathlib import Path
import tempfile
import unittest

try:
    import numpy as np
except ImportError:  # pragma: no cover - standard-library-only worker test
    np = None

from tools.full_spectrum.env_scatter_fit.parity import (
    _parse_cpp_result,
    deterministic_parity_queries,
    generate_parity_report,
)


def digest(label):
    return hashlib.sha256(label.encode("ascii")).hexdigest()


@unittest.skipIf(np is None, "NumPy is optional")
class CandidateParityTests(unittest.TestCase):
    @staticmethod
    def payload():
        from tools.full_spectrum.env_scatter_fit.fit import (
            PAIR_THRESHOLD_KEV,
            RuntimePayload,
            StructuredPositiveModel,
        )

        knots = [0.0] * 4 + [1.0] * 4
        q_min = [math.log(59.5), math.log(0.1), 0.0, 0.0]
        q_max = [math.log(2614.0), math.log(60.0), 1.55, 1.55]
        q_span = [high - low for low, high in zip(q_min, q_max)]
        def continuum(name, axis, threshold=None, scale=1.0, level=1):
            features = ((1 << (level - 1)) + 3) ** 4
            archetype = np.zeros((1, 4, 4))
            archetype[0, 1, 2] = 1.0
            return StructuredPositiveModel(
                name,
                axis,
                knots,
                knots,
                q_min,
                q_span,
                1,
                level,
                0.0,
                np.full((features, 1), math.log(scale)),
                np.zeros((features, 1)),
                archetype,
                pair_threshold_keV=threshold,
            )

        def line(name, axis, threshold=None, scale=1.0, level=1):
            features = ((1 << (level - 1)) + 3) ** 4
            archetype = np.zeros((1, 4))
            archetype[0, 2] = 1.0
            return StructuredPositiveModel(
                name,
                axis,
                [],
                knots,
                q_min,
                q_span,
                1,
                level,
                0.0,
                np.full((features, 1), math.log(scale)),
                np.zeros((features, 1)),
                archetype,
                pair_threshold_keV=threshold,
            )

        models = {
            "primary_order1": continuum(
                "primary_order1", "primary", scale=1.2, level=1
            ),
            "brems_order0": continuum(
                "brems_order0", "brems", 150.0, 0.4, level=3
            ),
            "annihilation_order1": continuum(
                "annihilation_order1", "annihilation", PAIR_THRESHOLD_KEV, 0.2,
                level=4,
            ),
        }
        lines = {
            "elastic_E0": line("elastic_E0", "line_e0", scale=0.3, level=2),
            "annihilation_511": line(
                "annihilation_511", "line_511", PAIR_THRESHOLD_KEV, 0.1,
                level=4,
            ),
        }
        inventory = {
            "continuum": [
                {"name": model.component, "axis_kind": model.axis_kind}
                for model in models.values()
            ],
            "lines": [
                {"name": model.component, "axis_kind": model.axis_kind}
                for model in lines.values()
            ],
        }
        return RuntimePayload(
            models,
            lines,
            {
                "config_sha256": digest("parity-config"),
                "material_name": "Soil",
                "geometry_tier": "field",
                "normalization": "s_hat_to_direct_fluence",
                "physical_domain": {
                    "min": [59.5, 0.5, 0.3, 0.3],
                    "max": [2614.0, 49.87, 4.16, 4.16],
                },
                "q_domain": {"min": q_min, "max": q_max},
                "brems_activation_threshold_keV": 150.0,
                "model_inventory": inventory,
                "validation_metadata": {},
            },
        )

    def test_query_bank_covers_faces_and_both_physics_thresholds(self):
        queries = deterministic_parity_queries(self.payload())
        counts = {}
        for query in queries:
            counts[query.category] = counts.get(query.category, 0) + 1
        self.assertGreaterEqual(counts["interior"], 1)
        self.assertGreaterEqual(counts["domain_face"], 4)
        self.assertEqual(counts["brems_threshold"], 3)
        self.assertEqual(counts["pair_threshold"], 3)
        self.assertEqual(len({query.label for query in queries}), len(queries))

    def test_cpp_result_parser_rejects_trailing_or_dimension_mismatch(self):
        good = "\n".join(
            (
                "CEEFIT-PARITY-RESULT-v1",
                f"MODEL {digest('model')}",
                "N 1",
                "Q 0 1 1 0 0 0.5",
                "TOTAL 1 2.0",
                "INCIDENT_TOTAL 1 1.0",
                "ENDQ",
                "END",
            )
        )
        model, results = _parse_cpp_result(good, 1)
        self.assertEqual(model, digest("model"))
        self.assertEqual(results[0]["total"], [2.0])
        with self.assertRaisesRegex(ValueError, "trailing"):
            _parse_cpp_result(good + " EXTRA", 1)
        with self.assertRaisesRegex(ValueError, "dimensions"):
            _parse_cpp_result(good.replace("TOTAL 1 2.0", "TOTAL 0"), 1)

    def test_exact_serialized_candidate_matches_built_cpp_evaluator(self):
        repository = Path(__file__).resolve().parents[4]
        executable = (
            repository / "build-rel" / "tools" / "full_spectrum" / "env_fit_payload_eval"
        )
        if not executable.is_file():
            self.skipTest(f"C++ parity evaluator has not been built at {executable}")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate = root / "candidate.ceefitp"
            output = root / "parity.json"
            self.payload().save_binary(candidate)
            report = generate_parity_report(candidate, executable, output)
            self.assertTrue(report["passed"])
            self.assertGreater(report["compared_value_count"], 100)
            self.assertLessEqual(report["max_absolute_difference"], 1.0e-9)
            restored = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(restored["schema"], "CEEFIT-PYTHON-CPP-PARITY-v1")
            self.assertEqual(restored["model_content_sha256"], report["model_content_sha256"])


if __name__ == "__main__":
    unittest.main()
