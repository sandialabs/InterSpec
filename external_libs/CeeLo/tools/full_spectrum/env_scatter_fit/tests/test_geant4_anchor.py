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
    import h5py
except ImportError:  # pragma: no cover - optional archive dependency
    h5py = None

from tools.full_spectrum.env_scatter_fit.archive import ingest_blocks
from tools.full_spectrum.env_scatter_fit.geant4_anchor import (
    MANIFEST_SCHEMA,
    REPORT_SCHEMA,
    generate_geant4_anchor_report,
    _parse_reference_arm,
)
from tools.full_spectrum.env_scatter_fit.cli import build_parser
from tools.full_spectrum.env_scatter_fit.tests.test_block import fixture_block


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _refresh_reference_hashes(record, root):
    for source in record["anchors"][0]["reference_files"]:
        for prefix in ("", "gdml_", "macro_"):
            path_key = "path" if not prefix else f"{prefix}path"
            hash_key = "sha256" if not prefix else f"{prefix}sha256"
            source[hash_key] = _hash(root / source[path_key])


def _write_arm(
    path,
    *,
    events,
    scattered,
    standard_error,
    rng_seed,
    energy_keV=100.0,
    world_half_extent_m=60.0,
    sphere_radius_cm=4.0,
    energy_min_keV=10.0,
    cell_standard_error=0.0,
):
    path = Path(path)
    gdml = path.with_suffix(".gdml")
    macro = path.with_suffix(".mac")
    world_half_cm = 100.0 * world_half_extent_m
    gdml.write_text(
        "\n".join(
            [
                '<?xml version="1.0" encoding="UTF-8"?>',
                "<gdml>",
                "  <!-- CEEFIT-GEOMETRY-v1 geometry_tier=air material=none "
                f"energy_keV={energy_keV:.17g} distance_cm=100 "
                "source_height_cm=0 detector_height_cm=0 "
                f"scoring_radius_cm={sphere_radius_cm:.17g} "
                f"world_half_extent_cm={world_half_cm:.17g} -->",
                "  <solids>",
                f'    <box name="world_s" x="{2 * world_half_cm:.17g}" '
                f'y="{2 * world_half_cm:.17g}" z="{2 * world_half_cm:.17g}" lunit="cm"/>',
                f'    <sphere name="score_s" rmin="0" rmax="{sphere_radius_cm:.17g}" lunit="cm"/>',
                "  </solids>",
                "  <structure>",
                '    <volume name="active_crystal"><materialref ref="G4_AIR"/>'
                '<solidref ref="score_s"/></volume>',
                '    <volume name="world"><materialref ref="G4_AIR"/>'
                '<solidref ref="world_s"/><physvol>'
                '<volumeref ref="active_crystal"/>'
                '<position name="score_p" x="0" y="0" z="0" unit="cm"/>'
                "</physvol></volume>",
                "  </structure>",
                "</gdml>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    macro.write_text(
        "\n".join(
            [
                "/gps/particle gamma",
                f"/gps/energy {energy_keV:.17g} keV",
                "/gps/pos/type Point",
                "/gps/pos/centre 0 0 -100 cm",
                "/gps/ang/type iso",
                f"/run/beamOn {events}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    volume = (4.0 / 3.0) * math.pi * sphere_radius_cm**3
    first = scattered * events * volume
    second = first * first / events + (
        standard_error * volume
    ) ** 2 * events * (events - 1)
    lines = [
        "# GEANT4 environment-scatter sphere flux",
        f"# E0_keV: {energy_keV:.17g}",
        f"# sphere_radius_cm: {sphere_radius_cm:.17g}",
        f"# gdml_scoring_radius_cm: {sphere_radius_cm:.17g}",
        "# observable: scattered-track-length-flux-per-primary-v1",
        "# physics_list: G4EmLivermorePhysics",
        f"# world_half_extent_m: {world_half_extent_m:.17g}",
        "# primary_sampling: gps-isotropic-4pi",
        "# rng_engine: MixMaxRng",
        f"# rng_seed: {rng_seed}",
        "# production_cut_mm: 0.7",
        f"# gdml_sha256: {_hash(gdml)}",
        f"# macro_sha256: {_hash(macro)}",
        "# geometry_tier: air",
        "# material: none",
        "# distance_cm: 100",
        "# source_height_cm: 0",
        "# detector_height_cm: 0",
        f"# events: {events}",
        "# scattered_total_moment_schema: per-event-track-length-v1",
        f"# scattered_total_N: {events}",
        f"# scattered_total_S1_cm: {first:.17g}",
        f"# scattered_total_S2_cm2: {second:.17g}",
        f"# scattered_total_mean_flux_cm-2: {scattered:.17g}",
        f"# scattered_total_standard_error_flux_cm-2: {standard_error:.17g}",
        "# n_energy: 48",
        f"# energy_min_keV: {energy_min_keV:.17g}",
        "# energy_max_keV: 3000",
        "# n_mu: 32",
        "# mu_min: -1",
        "# mu_max: 1",
        "energy_bin,mu_bin,e_low_keV,mu_low,phi_direct,phi_scattered,"
        "phi_scattered_sigma",
    ]
    for energy_bin in range(48):
        e_low = energy_min_keV * (3000.0 / energy_min_keV) ** (energy_bin / 48.0)
        for mu_bin in range(32):
            mu_low = -1.0 + 2.0 * mu_bin / 32.0
            first_cell = energy_bin == 0 and mu_bin == 0
            lines.append(
                f"{energy_bin},{mu_bin},{e_low:.17g},{mu_low:.17g},0,"
                f"{(scattered if first_cell else 0.0):.17g},"
                f"{(cell_standard_error if first_cell else 0.0):.17g}"
            )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gdml, macro


class Geant4AnchorCliSurfaceTests(unittest.TestCase):
    def test_cli_exposes_raw_reference_report_generator(self):
        args = build_parser().parse_args(
            [
                "geant4-anchor-report",
                "--archive",
                "fit.h5",
                "--manifest",
                "anchors.json",
                "--output",
                "agreement.json",
            ]
        )
        self.assertEqual(args.command, "geant4-anchor-report")
        self.assertEqual(args.archive, "fit.h5")
        self.assertEqual(args.manifest, "anchors.json")
        self.assertEqual(args.output, "agreement.json")

    def test_raw_arm_records_exact_gdml_and_macro_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            arm = root / "arm.csv"
            gdml, macro = _write_arm(
                arm,
                events=1000,
                scattered=1.0e-6,
                standard_error=1.0e-8,
                rng_seed=101,
            )
            parsed = _parse_reference_arm(arm)
            self.assertEqual(parsed["gdml_sha256"], _hash(gdml))
            self.assertEqual(parsed["macro_sha256"], _hash(macro))


@unittest.skipIf(h5py is None, "h5py is an optional archive dependency")
class Geant4AnchorReportTests(unittest.TestCase):
    def _fixture(self, root, **block_options):
        block = root / "block.ceefit"
        block.write_bytes(
            fixture_block(scoring_radius_cm=4.0, **block_options)
        )
        archive = root / "archive.h5"
        ingest_blocks((block,), archive)
        with h5py.File(archive, "r") as handle:
            node_id = f"{int(handle['nodes/id'][0]):016x}"
            count = int(handle["diagnostics/scattered_total_N"][0])
            mean = float(handle["diagnostics/scattered_total_S1"][0]) / count
        first = root / "g4-arm-1.csv"
        second = root / "g4-arm-2.csv"
        first_gdml, first_macro = _write_arm(
            first,
            events=1000,
            scattered=mean * 0.999,
            standard_error=mean * 0.01,
            rng_seed=101,
            cell_standard_error=mean * 100.0,
        )
        second_gdml, second_macro = _write_arm(
            second,
            events=3000,
            scattered=mean * (1.0 + 1.0 / 3000.0),
            standard_error=mean * 0.01,
            rng_seed=202,
            cell_standard_error=mean * 100.0,
        )
        manifest = root / "anchors.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema": MANIFEST_SCHEMA,
                    "physics_list": "G4EmLivermorePhysics",
                    "comparison_sigma": 2.0,
                    "anchors": [
                        {
                            "anchor_id": "field-soil-100keV",
                            "node_id": node_id,
                            "reference_files": [
                                {
                                    "path": first.name,
                                    "sha256": _hash(first),
                                    "gdml_path": first_gdml.name,
                                    "gdml_sha256": _hash(first_gdml),
                                    "macro_path": first_macro.name,
                                    "macro_sha256": _hash(first_macro),
                                },
                                {
                                    "path": second.name,
                                    "sha256": _hash(second),
                                    "gdml_path": second_gdml.name,
                                    "gdml_sha256": _hash(second_gdml),
                                    "macro_path": second_macro.name,
                                    "macro_sha256": _hash(second_macro),
                                },
                            ],
                        }
                    ],
                },
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        return archive, manifest, first, second

    def test_report_is_computed_from_raw_arms_and_archive_moments(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, _, _ = self._fixture(root)
            output = root / "report.json"

            report = generate_geant4_anchor_report(archive, manifest, output)

            self.assertEqual(report["schema"], REPORT_SCHEMA)
            self.assertTrue(report["passed"])
            self.assertTrue(report["within_combined_uncertainty"])
            self.assertEqual(report["anchor_count"], 1)
            self.assertEqual(report["reference_file_count"], 2)
            self.assertEqual(report["reference_event_count"], 4000)
            anchor = report["anchors"][0]
            self.assertAlmostEqual(anchor["reference"]["mean"], anchor["archive"]["mean"])
            self.assertEqual(anchor["reference"]["events"], 4000)
            self.assertLess(
                anchor["reference"]["per_history_standard_error"],
                anchor["reference"]["mean"] * 0.02,
            )
            self.assertEqual(json.loads(output.read_text(encoding="ascii")), report)

    def test_changed_reference_bytes_and_inadequate_world_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, _ = self._fixture(root)
            output = root / "report.json"
            first.write_text(first.read_text() + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash differs"):
                generate_geant4_anchor_report(archive, manifest, output)

            record = json.loads(manifest.read_text())
            record["anchors"][0]["reference_files"][0]["sha256"] = _hash(first)
            _write_arm(
                first,
                events=1000,
                scattered=1.0,
                standard_error=0.01,
                rng_seed=101,
                world_half_extent_m=15.0,
            )
            second = root / record["anchors"][0]["reference_files"][1]["path"]
            _write_arm(
                second,
                events=3000,
                scattered=1.0,
                standard_error=0.01,
                rng_seed=202,
                world_half_extent_m=15.0,
            )
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "world half extent"):
                generate_geant4_anchor_report(archive, manifest, output)

    def test_disagreement_is_reported_not_overridden_by_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, second = self._fixture(root)
            _write_arm(
                first, events=1000, scattered=0.0,
                standard_error=1.0e-30, rng_seed=101
            )
            _write_arm(
                second, events=3000, scattered=0.0,
                standard_error=1.0e-30, rng_seed=202
            )
            record = json.loads(manifest.read_text())
            _refresh_reference_hashes(record, root)
            record["passed"] = True
            manifest.write_text(json.dumps(record), encoding="utf-8")

            report = generate_geant4_anchor_report(
                archive, manifest, root / "report.json"
            )

            self.assertFalse(report["passed"])
            self.assertFalse(report["within_combined_uncertainty"])
            self.assertFalse(report["anchors"][0]["passed"])

    def test_reference_arms_must_be_independent_source_content(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, second = self._fixture(root)
            second.write_bytes(first.read_bytes())
            record = json.loads(manifest.read_text())
            for source in record["anchors"][0]["reference_files"]:
                source["sha256"] = _hash(root / source["path"])
            manifest.write_text(json.dumps(record), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "independent source bytes"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

    def test_exact_total_moments_and_distinct_rng_seeds_are_mandatory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, second = self._fixture(root)
            first.write_text(
                "\n".join(
                    line
                    for line in first.read_text(encoding="utf-8").splitlines()
                    if not line.startswith("# scattered_total_S2_cm2:")
                )
                + "\n",
                encoding="utf-8",
            )
            record = json.loads(manifest.read_text())
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "incomplete GEANT4 sphere-flux"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

            _write_arm(
                first, events=1000, scattered=1.0,
                standard_error=0.01, rng_seed=101
            )
            _write_arm(
                second, events=3000, scattered=1.0,
                standard_error=0.01, rng_seed=101
            )
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "distinct RNG seeds"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

    def test_cutoff_radius_and_hashed_geometry_must_match_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, second = self._fixture(root)
            record = json.loads(manifest.read_text())
            _write_arm(
                first, events=1000, scattered=1.0,
                standard_error=0.01, rng_seed=101, energy_min_keV=30.0
            )
            _write_arm(
                second, events=3000, scattered=1.0,
                standard_error=0.01, rng_seed=202, energy_min_keV=30.0
            )
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "energy cutoff does not match"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

            _write_arm(
                first, events=1000, scattered=1.0,
                standard_error=0.01, rng_seed=101, sphere_radius_cm=15.0
            )
            _write_arm(
                second, events=3000, scattered=1.0,
                standard_error=0.01, rng_seed=202, sphere_radius_cm=15.0
            )
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "scoring radius does not match"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

            # Alter the exact macro supplied to evidence while leaving the
            # source CSV's run-time hash unchanged.
            macro = first.with_suffix(".mac")
            macro.write_text(
                macro.read_text(encoding="utf-8").replace("-100 cm", "-90 cm"),
                encoding="utf-8",
            )
            record["anchors"][0]["reference_files"][0]["macro_sha256"] = _hash(macro)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not bind the exact macro"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

    def test_capped_archive_node_cannot_pass_anchor_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, _, _ = self._fixture(
                root,
                status="capped",
                capped=True,
                failure_reason="history cap",
            )
            report = generate_geant4_anchor_report(
                archive, manifest, root / "report.json"
            )
            self.assertFalse(report["passed"])
            self.assertFalse(report["anchors"][0]["passed"])

    def test_reference_grid_and_energy_must_match_the_archive_node(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, first, second = self._fixture(root)
            lines = first.read_text(encoding="utf-8").splitlines()
            first.write_text(
                "\n".join(line for line in lines if not line.startswith("47,31,"))
                + "\n",
                encoding="utf-8",
            )
            record = json.loads(manifest.read_text())
            record["anchors"][0]["reference_files"][0]["sha256"] = _hash(first)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "grid is incomplete"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

            _write_arm(
                first, events=1000, scattered=1.0,
                standard_error=0.01, rng_seed=101, energy_keV=101.0
            )
            _write_arm(
                second, events=3000, scattered=1.0,
                standard_error=0.01, rng_seed=202, energy_keV=101.0
            )
            record = json.loads(manifest.read_text())
            _refresh_reference_hashes(record, root)
            manifest.write_text(json.dumps(record), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "E0 does not match"):
                generate_geant4_anchor_report(archive, manifest, root / "report.json")

    def test_output_cannot_overwrite_manifest_or_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, manifest, _, _ = self._fixture(root)
            with self.assertRaisesRegex(ValueError, "must not overwrite"):
                generate_geant4_anchor_report(archive, manifest, manifest)
            with self.assertRaisesRegex(ValueError, "must not overwrite"):
                generate_geant4_anchor_report(archive, manifest, archive)


if __name__ == "__main__":
    unittest.main()
