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

import unittest

from tools.full_spectrum.env_scatter_fit.campaign import (
    MODEL_SELECTION_HISTORY_SCHEMA,
    MODEL_SELECTION_REPORT_SCHEMA,
)
from tools.full_spectrum.env_scatter_fit.workflow import (
    _apply_final_generation_cap,
    _campaign_stop_status,
)


def history_record(round_number, passed):
    return {
        "schema": MODEL_SELECTION_HISTORY_SCHEMA,
        "kind": "bound_model_selection_report",
        "passed": passed,
        "report_schema": MODEL_SELECTION_REPORT_SCHEMA,
        "report_canonical_sha256": "a" * 64,
        "config_sha256": "b" * 64,
        "archive_snapshot_sha256": "c" * 64,
        "model_content_sha256": "d" * 64,
        "eligible_validation_node_ids": ["1" * 16],
        "campaign_round": round_number,
        "final_generation": 0,
        "campaign_round_identity_sha256": "e" * 64,
    }


def campaign(history, *, consecutive, stopped):
    return {
        "config": {"sha256": "b" * 64},
        "design": {"nodes": {"validation": [{"node_id_hex": "1" * 16}]}},
        "campaign_state": {
            "round": len(history),
            "final_generation": 0,
            "stopped": stopped,
            "consecutive_validation_passes": consecutive,
            "validation_history": history,
        },
    }


class CampaignStopBoundaryTests(unittest.TestCase):
    def test_final_replacement_policy_stops_after_three_total_generations(self):
        policy = {
            "action": "retire_final_to_diagnostic_and_replace_from_unused_sobol_reserve"
        }
        current = {
            "allocation": {"maximum_final_generations": 3},
            "campaign_state": {"final_generation": 2},
        }
        bounded = _apply_final_generation_cap(policy, current)
        self.assertEqual(bounded["action"], "start_new_campaign_revision")
        self.assertEqual(bounded["unresolved_final_policy"], policy)

    def test_only_two_consecutive_recorded_passes_unlock_final_set(self):
        missing = _campaign_stop_status({})
        self.assertFalse(missing["passed"])

        one_pass = _campaign_stop_status(
            campaign(
                [history_record(1, False), history_record(2, True)],
                consecutive=1,
                stopped=False,
            )
        )
        self.assertFalse(one_pass["passed"])

        inconsistent = _campaign_stop_status(
            campaign(
                [history_record(1, True), history_record(2, False)],
                consecutive=2,
                stopped=True,
            )
        )
        self.assertFalse(inconsistent["passed"])

        stopped = _campaign_stop_status(
            campaign(
                [
                    history_record(1, False),
                    history_record(2, True),
                    history_record(3, True),
                ],
                consecutive=2,
                stopped=True,
            )
        )
        self.assertTrue(stopped["passed"])

    def test_malformed_history_fails_closed(self):
        result = _campaign_stop_status(
            {
                "config": {"sha256": "b" * 64},
                "design": {"nodes": {"validation": [{"node_id_hex": "1" * 16}]}},
                "campaign_state": {
                    "round": 2,
                    "final_generation": 0,
                    "stopped": True,
                    "consecutive_validation_passes": 2,
                    "validation_history": [1, 1],
                }
            }
        )
        self.assertFalse(result["passed"])
        self.assertIn("malformed", result["reason"])

    def test_unbound_boolean_history_cannot_unlock_final_set(self):
        result = _campaign_stop_status(
            {
                "config": {"sha256": "b" * 64},
                "design": {"nodes": {"validation": [{"node_id_hex": "1" * 16}]}},
                "campaign_state": {
                    "round": 2,
                    "final_generation": 0,
                    "stopped": True,
                    "consecutive_validation_passes": 2,
                    "validation_history": [True, True],
                },
            }
        )
        self.assertFalse(result["passed"])


if __name__ == "__main__":
    unittest.main()
