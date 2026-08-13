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

"""CI gate for stale generation and prohibited historical G4 payloads."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import struct
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "generated_manifest.json"
ACCEPTANCE = HERE / "reports" / "migration_acceptance.json"
PROHIBITED_PATH_PARTS = (
    "tools/extract_g4_xs/",
    "sb_tables_g4emlow.py",
    "g4_xs_cache.csv",
    "g4_Kpe_cache.csv",
    "g4_fluor_cache.csv",
    "g4_fluor_l_cache.csv",
)
PROHIBITED_BANNERS = (
    b"*** WARNING: DERIVED " + b"FROM GEANT4 G4EMLOW",
    b"float32 restoration " + b"archive",
)
PROHIBITED_FILE_HASHES = {
    "0d06331df264c52aae789b3824393a9b9201bced103e11ce659075f2d37a8614",
    "27c72b79c3a87a5e5b3badbcf95d7fde440f0ca1e19ccfaffec9b70abe79c9a5",
    "f196ba9342ab69e53112954eb88e8543170ac6c8f0439f813312f63e1acb9449",
    "cbcc960f24c2aeb4d271db3efbb72d3b42af38b35d069a2a0e4a323bf992d5f6",
}
PROHIBITED_DECODED_SB_HASHES = {
    "68f552be7454fb57a536439fe80717b97830f1155bbf4563d53b4c57faf3540a",
    "f521b1338829e2bf352f3c2e6d6b403bfe81c6540458d4d3a231746b27ae1993",
}
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Directories that are never part of the distributed source and can be large;
# skipped by the non-git fallback so a vendored checkout does not walk build trees.
_SKIP_DIRS = {
    "build", "build_g4xs", "studies", "scratch", "docs", "venv",
    "__pycache__", "node_modules", "Testing",
}


def _skip_dir(name: str) -> bool:
    # Dot-directories are tool/VCS state, never distributed source.
    return name.startswith(".") or name in _SKIP_DIRS or name.startswith("build")


def working_files() -> list[str]:
    """Repo-relative paths to check.

    Prefers `git ls-files` so a development checkout honours .gitignore. Falls
    back to a filesystem walk when git is unavailable or ROOT is not a git
    repository -- the case for a vendored copy, an unpacked release tarball, or
    a `git archive` extraction, where this check still has to run.
    """
    try:
        output = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, text=True,
        ).stdout
        return sorted(set(output.splitlines()))
    except (OSError, subprocess.CalledProcessError):
        pass

    paths: list[str] = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not _skip_dir(d)]
        for name in filenames:
            paths.append(
                str(Path(dirpath, name).relative_to(ROOT).as_posix()))
    return sorted(set(paths))


def sb_payload_hash(path: Path) -> str | None:
    text = path.read_text(encoding="utf-8")
    is_float = "static const float Z1_sb_chi" in text
    is_quantized = "static const uint16_t Z1_sb_chi" in text
    if not is_float and not is_quantized:
        return None
    scales: list[float] = []
    if is_quantized:
        match = re.search(
            r"extern const float kSB_chi_scale\[[^]]*\]\s*=\s*\{(.*?)\};",
            text, re.DOTALL,
        )
        if not match:
            raise ValueError("cannot parse quantized SB scales")
        scales = [float(token) for token in NUMBER.findall(match.group(1))]
        if len(scales) != 92:
            raise ValueError("quantized SB scale count is not 92")
    digest = hashlib.sha256()
    for z in range(1, 93):
        match = re.search(
            rf"static const (?:float|uint16_t) Z{z}_sb_chi\[[^]]*\]\s*=\s*\{{(.*?)\}};",
            text, re.DOTALL,
        )
        if not match:
            raise ValueError(f"cannot parse Z={z} SB payload")
        body = re.sub(r"//.*", "", match.group(1))
        for token in NUMBER.findall(body):
            value = float(token) * scales[z - 1] if is_quantized else float(token)
            digest.update(struct.pack("<f", value))
    return digest.hexdigest()


def accepted_gate_failures(errors: list[str]) -> dict[str, dict[str, object]]:
    if not ACCEPTANCE.is_file():
        return {}
    try:
        payload = json.loads(ACCEPTANCE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"cannot read migration acceptance: {exc}")
        return {}
    if payload.get("status") != "accepted":
        errors.append("migration acceptance status is not accepted")
        return {}
    records = payload.get("accepted_reports")
    if not isinstance(records, dict):
        errors.append("migration acceptance has no accepted_reports mapping")
        return {}
    valid: dict[str, dict[str, object]] = {}
    for relative, record in records.items():
        if not isinstance(relative, str):
            errors.append("migration acceptance report path is not a string")
            continue
        if not isinstance(record, dict):
            errors.append(f"invalid migration acceptance record: {relative}")
            continue
        if not record.get("rationale"):
            errors.append(f"migration acceptance has no rationale: {relative}")
        if not record.get("sha256"):
            errors.append(f"migration acceptance has no report hash: {relative}")
        valid[relative] = record
    return valid


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-manual-approval", action="store_true")
    args = parser.parse_args()
    errors: list[str] = []
    accepted_failures = accepted_gate_failures(errors)
    paths = working_files()
    for relative in paths:
        path = ROOT / relative
        if not path.is_file():
            continue
        if any(part in relative for part in PROHIBITED_PATH_PARTS):
            errors.append(f"prohibited historical path remains: {relative}")
            continue
        data = path.read_bytes()
        if any(banner in data for banner in PROHIBITED_BANNERS):
            errors.append(f"prohibited historical banner remains: {relative}")
        if hashlib.sha256(data).hexdigest() in PROHIBITED_FILE_HASHES:
            errors.append(f"known prohibited cache payload remains: {relative}")

    element_path = ROOT / "src/cross_sections/element_data.cpp"
    payload_hash = sb_payload_hash(element_path)
    if payload_hash in PROHIBITED_DECODED_SB_HASHES:
        errors.append("element_data.cpp contains a known historical G4 SB payload")

    lock_hash = sha256(HERE / "sources.lock.json")
    reports: dict[str, tuple[Path, dict[str, object]]] = {}
    for report_path in sorted((HERE / "reports").glob("*.json")):
        report = json.loads(report_path.read_text(encoding="utf-8"))
        relative = str(report_path.relative_to(ROOT))
        reports[relative] = (report_path, report)
        recorded = report.get("source_lock_sha256")
        if recorded is not None and recorded != lock_hash:
            errors.append(f"stale source-lock hash in {relative}")
        if report.get("passed") is False:
            acceptance = accepted_failures.get(relative)
            if acceptance is None:
                errors.append(f"unaccepted numerical gate failure in {relative}")
            elif acceptance.get("sha256") != sha256(report_path):
                errors.append(f"stale migration acceptance for {relative}")

    for relative, acceptance in accepted_failures.items():
        entry = reports.get(relative)
        if entry is None:
            errors.append(f"accepted migration report is missing: {relative}")
        elif entry[1].get("passed") is not False:
            errors.append(f"obsolete migration acceptance for passing report: {relative}")
        elif acceptance.get("sha256") != sha256(entry[0]):
            # This can also be reported above, but keep this check for accepted
            # entries that are malformed or otherwise skipped by report parsing.
            if f"stale migration acceptance for {relative}" not in errors:
                errors.append(f"stale migration acceptance for {relative}")

    if not MANIFEST.is_file():
        errors.append("generated_manifest.json is missing")
    else:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for relative, expected in manifest.get("files", {}).items():
            path = ROOT / relative
            if not path.is_file():
                errors.append(f"manifest output missing: {relative}")
            elif sha256(path) != expected:
                errors.append(f"stale generated output/report: {relative}")

    review = ROOT / "validation/bremsstrahlung_migration_review/manual_review.json"
    if args.require_manual_approval:
        status = json.loads(review.read_text(encoding="utf-8")).get("status") if review.is_file() else None
        if status != "approved":
            errors.append("bremsstrahlung manual review is not approved")

    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("Cross-section provenance, prohibited-payload, and stale-generation checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
