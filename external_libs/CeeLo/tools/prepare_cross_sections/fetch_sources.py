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

"""Fetch and verify the preparation-only inputs in ``sources.lock.json``.

Normal CeeLo builds never invoke this script and never require network access.
Downloaded archives and extracted trees live under the gitignored preparation
cache.  Every byte used by a generator is checked against the committed lock.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import re
import struct
import sys
import tarfile
import tempfile
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
DEFAULT_LOCK = HERE / "sources.lock.json"
DEFAULT_CACHE = HERE / "data" / "source-cache"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_lock(path: Path) -> dict[str, Any]:
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("schema_version") != 1 or not isinstance(doc.get("sources"), dict):
        raise ValueError(f"unsupported or malformed source lock: {path}")
    return doc


def verify_file(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label} SHA-256 mismatch: expected {expected}, got {actual}")


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(fd, "wb") as out, urllib.request.urlopen(url) as response:
            shutil.copyfileobj(response, out, length=1024 * 1024)
        os.replace(temporary, destination)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def safe_extract_tar(archive: Path, destination: Path) -> None:
    """Extract a tar archive after rejecting path and link escapes."""
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive, "r:*") as source:
        for member in source.getmembers():
            target = (destination / member.name).resolve()
            if target != root and root not in target.parents:
                raise ValueError(f"archive member escapes cache: {member.name!r}")
            if member.issym() or member.islnk():
                link = (target.parent / member.linkname).resolve()
                if link != root and root not in link.parents:
                    raise ValueError(f"archive link escapes cache: {member.name!r}")
        source.extractall(destination)


def epq_table_set_sha256(directory: Path, pattern: str) -> str:
    """Hash the logical EPQ table set independent of archive metadata."""
    paths = sorted(directory.glob(pattern), key=lambda item: item.name)
    if len(paths) != 92:
        raise ValueError(f"expected 92 EPQ tables in {directory}, found {len(paths)}")
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode("ascii"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def verify_extracted(source_name: str, spec: dict[str, Any], root: Path) -> None:
    for relative, expected in spec.get("file_sha256", {}).items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"{source_name}: missing extracted file {relative}")
        verify_file(path, expected, f"{source_name}:{relative}")

    expected_tables = spec.get("table_set_sha256")
    if expected_tables:
        directory = root / spec["table_subdir"]
        actual = epq_table_set_sha256(directory, spec["table_glob"])
        if actual != expected_tables:
            raise ValueError(
                f"{source_name} table-set SHA-256 mismatch: expected "
                f"{expected_tables}, got {actual}"
            )


def parse_estar_text(path: Path) -> list[tuple[float, float, float]]:
    row_re = re.compile(
        r"\s+([0-9.]+E[+-][0-9]+)\s+" + r"\s+".join(
            [r"([0-9.]+E[+-][0-9]+)"] * 6
        ) + r"\s*<br>"
    )
    rows = [
        (float(match.group(1)), float(match.group(2)), float(match.group(3)))
        for match in row_re.finditer(path.read_text(encoding="ascii"))
        if 0.01 <= float(match.group(1)) <= 20.0
    ]
    if len(rows) != 53:
        raise ValueError(f"{path}: expected 53 ESTAR rows, found {len(rows)}")
    return rows


def prepare_estar_text_source(
    source_name: str, spec: dict[str, Any], cache: Path, offline: bool
) -> Path:
    destination = cache / "generated" / source_name
    destination.mkdir(parents=True, exist_ok=True)
    for z in range(1, 93):
        path = destination / spec["filename_pattern"].replace("NNN", f"{z:03d}")
        if not path.is_file():
            if offline:
                raise FileNotFoundError(f"offline mode: {path} is not cached")
            request = urllib.request.Request(
                spec["url"],
                data=urllib.parse.urlencode({"matno": f"{z:03d}", "ShowDefault": "on"}).encode("ascii"),
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=destination)
            try:
                with os.fdopen(fd, "wb") as out, urllib.request.urlopen(request) as response:
                    shutil.copyfileobj(response, out)
                os.replace(temporary, path)
            except Exception:
                try:
                    os.unlink(temporary)
                except FileNotFoundError:
                    pass
                raise
    digest = hashlib.sha256()
    for z in range(1, 93):
        path = destination / spec["filename_pattern"].replace("NNN", f"{z:03d}")
        digest.update(struct.pack("<H", z))
        for row in parse_estar_text(path):
            digest.update(struct.pack("<ddd", *row))
    actual = digest.hexdigest()
    if actual != spec["logical_table_sha256"]:
        raise ValueError(
            f"{source_name} logical SHA-256 mismatch: expected "
            f"{spec['logical_table_sha256']}, got {actual}"
        )
    print(f"Verified {source_name}: {destination}", file=sys.stderr)
    return destination


def prepare_source(
    source_name: str,
    spec: dict[str, Any],
    cache: Path,
    offline: bool,
    force_extract: bool,
) -> Path:
    if spec.get("filename_pattern") and spec.get("logical_table_sha256"):
        return prepare_estar_text_source(source_name, spec, cache, offline)
    archive = cache / "downloads" / spec["filename"]
    if not archive.exists():
        if offline:
            raise FileNotFoundError(f"offline mode: {archive} is not cached")
        print(f"Downloading {source_name}: {spec['url']}", file=sys.stderr)
        download(spec["url"], archive)
    verify_file(archive, spec["sha256"], source_name)

    if not spec.get("archive"):
        print(f"Verified {source_name}: {archive}", file=sys.stderr)
        return archive

    extracted = cache / "extracted" / source_name
    marker = extracted / ".source.sha256"
    if force_extract or not marker.is_file() or marker.read_text().strip() != spec["sha256"]:
        if extracted.exists():
            shutil.rmtree(extracted)
        safe_extract_tar(archive, extracted)
        marker.write_text(spec["sha256"] + "\n", encoding="ascii")
    verify_extracted(source_name, spec, extracted)
    print(f"Verified {source_name}: {extracted}", file=sys.stderr)
    return extracted


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--source",
        action="append",
        dest="sources",
        help="source-lock key to prepare; repeat as needed (default: all)",
    )
    parser.add_argument("--offline", action="store_true", help="never access the network")
    parser.add_argument("--force-extract", action="store_true")
    args = parser.parse_args()

    lock = load_lock(args.lock)
    selected = args.sources or list(lock["sources"])
    unknown = sorted(set(selected) - set(lock["sources"]))
    if unknown:
        parser.error("unknown source(s): " + ", ".join(unknown))

    for name in selected:
        prepare_source(name, lock["sources"][name], args.cache, args.offline, args.force_extract)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
