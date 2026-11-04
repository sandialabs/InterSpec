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
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

try:
    import h5py
except ImportError:  # pragma: no cover - dependency-free test environment
    h5py = None
try:
    import numpy
except ImportError:  # pragma: no cover - minimal dependency-free environment
    numpy = None

from tools.full_spectrum.env_scatter_fit.archive import ingest_blocks
from tools.full_spectrum.env_scatter_fit.calibration import (
    CALIBRATION_SCHEMA,
    _FLOAT32_REPLICA_NULL_TOLERANCE,
    MINIMUM_COMPARISONS,
    calibrate_uncertainty_archive,
    _cell_moments,
    _mahalanobis,
)
from tools.full_spectrum.env_scatter_fit.cli import build_parser
from tools.full_spectrum.env_scatter_fit.design import PhysicalNode
from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block


CALIBRATION_NODE = PhysicalNode(1500.0, 1.0, 0.0, 0.0)


def _build_archive(
    directory,
    *,
    variant="valid",
    block_count=2 * MINIMUM_COMPARISONS,
    reverse=False,
    name="archive.h5",
    scale=1.0,
    physical_node=CALIBRATION_NODE,
):
    # Every block stays a raw record: these tests exercise the estimators on
    # many repeated-null pairs at one node, where production folds blocks
    # beyond RETAINED_RAW_BLOCKS.
    paths = []
    for block_id in range(1, block_count + 1):
        path = Path(directory) / f"{name}.{block_id}.ceefit"
        path.write_bytes(
            fixture_block(
                block_id=block_id,
                physical_node=physical_node,
                histories=262_144,
                energy_bins=12,
                calibration_variant=variant,
                scale=scale,
            )
        )
        paths.append(path)
    archive = Path(directory) / name
    ingest_blocks(
        tuple(reversed(paths)) if reverse else tuple(paths), archive, fold_blocks=False
    )
    return archive


class CalibrationCliSurfaceTest(unittest.TestCase):
    def test_parser_exposes_archive_bound_calibration_command(self):
        args = build_parser().parse_args(
            [
                "calibrate-uncertainty",
                "--archive",
                "fit.h5",
                "--output",
                "calibration.json",
            ]
        )
        self.assertEqual(args.command, "calibrate-uncertainty")
        self.assertEqual(args.archive, "fit.h5")
        self.assertEqual(args.output, "calibration.json")

    def test_cli_module_invocation_executes_main(self):
        completed = subprocess.run(
            [sys.executable, "-m", "tools.full_spectrum.env_scatter_fit.cli", "--help"],
            cwd=Path(__file__).resolve().parents[4],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("calibrate-uncertainty", completed.stdout)

    @unittest.skipIf(numpy is None, "NumPy is required for moment arithmetic")
    def test_real_transport_scale_s2_remains_stochastic(self):
        group = {
            "moments/continuum_total": {
                "N": numpy.asarray([262_144], dtype="u8"),
                "S1": numpy.asarray([2.825141e-08], dtype="f8"),
                "S2": numpy.asarray([4.635103e-16], dtype="f8"),
            }
        }
        mean, standard_error = _cell_moments(group, "continuum_total", numpy)
        self.assertGreater(float(mean[0]), 0.0)
        self.assertGreater(float(standard_error[0]), 0.0)

    @unittest.skipIf(numpy is None, "NumPy is required for moment arithmetic")
    def test_subnormal_s2_retains_a_representable_standard_error(self):
        group = {
            "moments/continuum_total": {
                "N": numpy.asarray([262_144], dtype="u8"),
                "S1": numpy.asarray([5.199175329330277e-160], dtype="f8"),
                "S2": numpy.asarray([2.70313e-319], dtype="f8"),
            }
        }
        mean, standard_error = _cell_moments(group, "continuum_total", numpy)
        self.assertGreater(float(mean[0]), 0.0)
        self.assertGreater(float(standard_error[0]), 0.0)

    @unittest.skipIf(numpy is None, "NumPy is required for covariance arithmetic")
    def test_float32_replica_closure_uses_its_storage_precision(self):
        difference = numpy.asarray([1.0, 5.0e-7], dtype="f8")
        covariance = numpy.asarray([[1.0, 0.0], [0.0, 0.0]], dtype="f8")
        with self.assertRaisesRegex(ValueError, "outside covariance support"):
            _mahalanobis(difference, covariance, numpy)
        statistic, rank, _ = _mahalanobis(
            difference,
            covariance,
            numpy,
            null_relative_tolerance=_FLOAT32_REPLICA_NULL_TOLERANCE,
        )
        self.assertEqual(rank, 1)
        self.assertEqual(statistic, 1.0)

    @unittest.skipIf(numpy is None, "NumPy is required for covariance arithmetic")
    def test_resolved_low_variance_mode_is_not_treated_as_null_space(self):
        covariance = numpy.asarray([[1.0, 0.0], [0.0, 1.0e-13]], dtype="f8")
        difference = numpy.asarray([0.0, 0.8 * numpy.sqrt(1.0e-13)], dtype="f8")
        statistic, rank, standardized = _mahalanobis(
            difference, covariance, numpy
        )
        self.assertEqual(rank, 2)
        self.assertAlmostEqual(statistic, 0.64)
        self.assertEqual(len(standardized), 2)

    @unittest.skipIf(numpy is None, "NumPy is required for covariance arithmetic")
    def test_unresolved_covariance_mode_uses_rank_resolution_tolerance(self):
        covariance = numpy.asarray([[1.0, 0.0], [0.0, 1.0e-18]], dtype="f8")
        difference = numpy.asarray([1.0, 0.8e-9], dtype="f8")
        statistic, rank, _ = _mahalanobis(difference, covariance, numpy)
        self.assertEqual(rank, 1)
        self.assertAlmostEqual(statistic, 1.0)


@unittest.skipIf(h5py is None, "h5py and NumPy are optional archive dependencies")
class UncertaintyCalibrationTest(unittest.TestCase):
    def test_production_scale_fluence_does_not_collapse_s2_to_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            # Per-history means near 1e-13 produce block S1 near 1e-8, the
            # scale measured in real transport by REVIEW_e5f76aa.md.  The
            # standardized null experiment is scale invariant and must pass.
            archive = _build_archive(directory, scale=1.0e-13)
            output = Path(directory) / "production-scale.json"

            report = calibrate_uncertainty_archive(archive, output)

            self.assertTrue(report["passed"])
            self.assertTrue(report["s2_passed"])
            self.assertTrue(report["replica_covariance_passed"])

    def test_valid_repeated_null_archive_covers_every_family(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = _build_archive(directory)
            output = Path(directory) / "calibration.json"

            report = calibrate_uncertainty_archive(archive, output)

            self.assertEqual(report["schema"], CALIBRATION_SCHEMA)
            self.assertTrue(report["passed"])
            self.assertEqual(report["repeated_null_runs"], MINIMUM_COMPARISONS)
            self.assertTrue(report["s2_passed"])
            self.assertTrue(report["coarse_covariance_passed"])
            self.assertTrue(report["replica_covariance_passed"])
            self.assertTrue(report["family_coverage_passed"])
            self.assertEqual(report["failures"], [])
            self.assertEqual(
                set(report["cell_s2"]["tensors"]),
                {
                    "continuum_components",
                    "continuum_total",
                    "line_components",
                    "line_total",
                },
            )
            self.assertTrue(
                all(
                    tensor["passed"]
                    for tensor in report["cell_s2"]["tensors"].values()
                )
            )
            self.assertEqual(
                set(report["replica_covariance"]["families"]),
                set(report["cell_s2"]["tensors"]),
            )
            self.assertTrue(
                all(
                    family["passed"]
                    for family in report["replica_covariance"]["families"].values()
                )
            )
            self.assertEqual(report["replication"]["blocks_used"], 16)
            self.assertEqual(
                len(
                    report["replica_covariance"]["projection_definition"][
                        "projection_sha256"
                    ]
                ),
                64,
            )
            self.assertEqual(len(report["stopping_outer_moments"]["dimensions"]), 16)
            self.assertTrue(
                all(
                    dimension["passed"]
                    for dimension in report["stopping_outer_moments"]["dimensions"]
                )
            )
            encoded = output.read_text(encoding="ascii")
            self.assertNotIn("NaN", encoded)
            self.assertNotIn("Infinity", encoded)
            self.assertEqual(json.loads(encoded), report)
            self.assertEqual(list(Path(directory).glob("*.tmp")), [])

    def test_each_uncertainty_family_has_an_effective_corruption_gate(self):
        cases = (
            ("bad_s2", ("raw_total_s2",)),
            ("bad_component_s2", ("cell_s2", "continuum_components")),
            ("bad_total_cell_s2", ("cell_s2", "continuum_total")),
            ("bad_line_s2", ("cell_s2", "line_total")),
            ("bad_stopping_covariance", ("stopping_outer_moments",)),
            ("bad_replica_covariance", ("replica_covariance",)),
        )
        with tempfile.TemporaryDirectory() as directory:
            for variant, path in cases:
                with self.subTest(variant=variant):
                    archive = _build_archive(
                        directory, variant=variant, name=f"{variant}.h5"
                    )
                    output = Path(directory) / f"{variant}.json"
                    report = calibrate_uncertainty_archive(archive, output)
                    selected = report
                    for key in path:
                        if key in {"continuum_components", "continuum_total", "line_total"}:
                            selected = selected["tensors"][key]
                        else:
                            selected = selected[key]
                    self.assertFalse(report["passed"])
                    self.assertFalse(selected["passed"])
                    self.assertEqual(json.loads(output.read_text()), report)

    def test_supported_but_statistically_empty_component_fails_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = _build_archive(directory, variant="empty_supported_component")
            output = Path(directory) / "empty_component.json"

            report = calibrate_uncertainty_archive(archive, output)

            segments = {
                segment["name"]: segment
                for segment in report["cell_s2"]["tensors"][
                    "continuum_components"
                ]["segments"]
            }
            self.assertTrue(segments["brems_order2plus"]["required"])
            self.assertEqual(
                segments["brems_order2plus"]["stochastic_residual_count"], 0
            )
            self.assertFalse(segments["brems_order2plus"]["passed"])
            self.assertFalse(report["s2_passed"])
            self.assertFalse(report["family_coverage_passed"])
            self.assertFalse(report["passed"])

    def test_subthreshold_archive_explains_missing_annihilation_support(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = _build_archive(
                directory,
                physical_node=PhysicalNode(662.0, 1.0, 0.0, 0.0),
            )
            output = Path(directory) / "subthreshold.json"

            report = calibrate_uncertainty_archive(archive, output)

            self.assertFalse(report["passed"])
            annihilation = [
                failure
                for failure in report["failures"]
                if "annihilation" in str(failure.get("family", ""))
            ]
            self.assertTrue(annihilation)
            self.assertTrue(
                all(
                    "8 independent block-pair comparisons" in item["reason"]
                    and "1021.9979-keV pair threshold" in item["reason"]
                    for item in annihilation
                )
            )

    def test_seven_independent_comparisons_are_not_formal_calibration(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = _build_archive(
                directory, block_count=2 * (MINIMUM_COMPARISONS - 1)
            )
            output = Path(directory) / "insufficient.json"

            report = calibrate_uncertainty_archive(archive, output)

            self.assertFalse(report["passed"])
            self.assertEqual(
                report["repeated_null_runs"], MINIMUM_COMPARISONS - 1
            )
            self.assertFalse(report["replication"]["sufficient"])
            self.assertFalse(report["s2_passed"])
            self.assertFalse(report["coarse_covariance_passed"])
            self.assertFalse(report["replica_covariance_passed"])

    def test_output_cannot_overwrite_the_source_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = _build_archive(directory)

            with self.assertRaisesRegex(ValueError, "must be distinct"):
                calibrate_uncertainty_archive(archive, archive)

            with h5py.File(archive, "r") as root:
                self.assertEqual(root.attrs["schema"], "CEEFIT-HDF5-v1")

    def test_archive_merge_and_group_order_do_not_change_report(self):
        with tempfile.TemporaryDirectory() as directory:
            first_archive = _build_archive(directory, name="forward.h5")
            second_archive = _build_archive(
                directory, reverse=True, name="reverse.h5"
            )
            first_output = Path(directory) / "forward.json"
            second_output = Path(directory) / "reverse.json"

            first = calibrate_uncertainty_archive(first_archive, first_output)
            second = calibrate_uncertainty_archive(second_archive, second_output)

            self.assertEqual(
                first["archive_snapshot_sha256"],
                second["archive_snapshot_sha256"],
            )
            self.assertEqual(first, second)
            self.assertEqual(first_output.read_bytes(), second_output.read_bytes())


if __name__ == "__main__":
    unittest.main()
