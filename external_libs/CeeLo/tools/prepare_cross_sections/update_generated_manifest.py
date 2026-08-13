#!/usr/bin/env python3
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

"""Record hashes for generated cross-section sources and coefficient reports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "tools" / "prepare_cross_sections" / "generated_manifest.json"
PATHS = (
    "tools/prepare_cross_sections/sources.lock.json",
    "tools/prepare_cross_sections/source_notices/NIST_EPQ_LicenseFile.txt",
    "tools/prepare_cross_sections/source_notices/NIST_STAR_Disclaimer.txt",
    "tools/prepare_cross_sections/source_notices/xraylib_license_all.txt",
    "tools/prepare_cross_sections/fetch_sources.py",
    "tools/prepare_cross_sections/generation_utils.py",
    "tools/prepare_cross_sections/epics_endf.py",
    "tools/prepare_cross_sections/sb_tables_epq.py",
    "tools/prepare_cross_sections/generate_element_support.py",
    "tools/prepare_cross_sections/generate_epics_photon.py",
    "tools/prepare_cross_sections/generate_eadl_relaxation.py",
    "tools/prepare_cross_sections/generate_estar_stopping.py",
    "tools/prepare_cross_sections/coefficient_accounting.py",
    "tools/prepare_cross_sections/check_provenance.py",
    "tools/prepare_cross_sections/update_generated_manifest.py",
    "tools/geant4_validation/generate_brems_review_bundle.py",
    "tools/geant4_validation/validate_photon_migration.py",
    "tools/geant4_validation/validate_stopping_migration.py",
    "src/cross_sections/element_data.cpp",
    "src/cross_sections/photon_epics_data.cpp",
    "src/cross_sections/photon_epics_data.h",
    "src/cross_sections/relaxation_epics_data.cpp",
    "src/cross_sections/relaxation_epics_data.h",
    "src/physics/estar_stopping_data.cpp",
    "src/physics/estar_stopping_data.h",
    "tools/prepare_cross_sections/reports/bremsstrahlung_epq.json",
    "tools/prepare_cross_sections/reports/bremsstrahlung_epq.md",
    "tools/prepare_cross_sections/reports/element_support.json",
    "tools/prepare_cross_sections/reports/photon_epics.csv",
    "tools/prepare_cross_sections/reports/photon_epics.json",
    "tools/prepare_cross_sections/reports/photon_epics.md",
    "tools/prepare_cross_sections/reports/photon_migration.json",
    "tools/prepare_cross_sections/reports/photon_migration.md",
    "tools/prepare_cross_sections/reports/relaxation_eadl.json",
    "tools/prepare_cross_sections/reports/stopping_estar.csv",
    "tools/prepare_cross_sections/reports/stopping_estar.json",
    "tools/prepare_cross_sections/reports/stopping_estar.md",
    "tools/prepare_cross_sections/reports/stopping_migration.json",
    "tools/prepare_cross_sections/reports/stopping_migration.md",
    "tools/prepare_cross_sections/reports/migration_acceptance.json",
    "tools/prepare_cross_sections/reports/coefficient_sizes.csv",
    "tools/prepare_cross_sections/reports/coefficient_sizes.json",
    "tools/prepare_cross_sections/reports/coefficient_sizes.md",
    "validation/bremsstrahlung_migration_review/README.md",
    "validation/bremsstrahlung_migration_review/metadata.json",
    "validation/bremsstrahlung_migration_review/review.html",
    "validation/bremsstrahlung_migration_review/spectra.csv",
    "validation/bremsstrahlung_migration_review/summary.csv",
    "validation/bremsstrahlung_migration_review/worst_20.csv",
    "validation/bremsstrahlung_migration_review/manual_review.json",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    missing = [relative for relative in PATHS if not (ROOT / relative).is_file()]
    if missing:
        raise FileNotFoundError("generated manifest inputs are missing: " + ", ".join(missing))
    manifest = {
        "schema_version": 1,
        "purpose": "CI stale-generation gate; update only after verified regeneration",
        "files": {relative: sha256(ROOT / relative) for relative in PATHS},
    }
    OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Updated {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
