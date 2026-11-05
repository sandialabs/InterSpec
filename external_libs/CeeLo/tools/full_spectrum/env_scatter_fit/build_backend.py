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

"""Tiny dependency-free PEP 517 backend for the single-package tool bundle."""

from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import io
from pathlib import Path
import tarfile
import zipfile


NAME = "ceelo-env-scatter-fit"
VERSION = "0.1.0"
NORMALIZED = "ceelo_env_scatter_fit"
DIST_INFO = f"{NORMALIZED}-{VERSION}.dist-info"
WHEEL_NAME = f"{NORMALIZED}-{VERSION}-py3-none-any.whl"


def _metadata() -> bytes:
    return (
        "Metadata-Version: 2.1\n"
        f"Name: {NAME}\n"
        f"Version: {VERSION}\n"
        "Summary: Fit-first environmental-scatter campaign and surrogate tools\n"
        "Requires-Python: >=3.10\n"
        "Provides-Extra: fit\n"
        "Requires-Dist: numpy>=1.26; extra == 'fit'\n"
        "Provides-Extra: archive\n"
        "Requires-Dist: numpy>=1.26; extra == 'archive'\n"
        "Requires-Dist: h5py>=3.10; extra == 'archive'\n"
        "\n"
    ).encode("utf-8")


def _wheel_metadata() -> bytes:
    return (
        "Wheel-Version: 1.0\n"
        "Generator: ceefit-stdlib-backend 1\n"
        "Root-Is-Purelib: true\n"
        "Tag: py3-none-any\n\n"
    ).encode("utf-8")


def _entry_points() -> bytes:
    return b"[console_scripts]\nceefit = env_scatter_fit.cli:main\n"


def _package_files() -> list[tuple[str, bytes]]:
    root = Path(__file__).resolve().parent
    result = []
    for path in sorted(root.glob("*.py")):
        if path.name == "build_backend.py":
            continue
        result.append((f"env_scatter_fit/{path.name}", path.read_bytes()))
    return result


def _record_line(path: str, data: bytes) -> tuple[str, str, str]:
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode("ascii")
    return path, f"sha256={digest}", str(len(data))


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    del config_settings, metadata_directory
    destination = Path(wheel_directory) / WHEEL_NAME
    files = _package_files()
    files.extend(
        (
            (f"{DIST_INFO}/METADATA", _metadata()),
            (f"{DIST_INFO}/WHEEL", _wheel_metadata()),
            (f"{DIST_INFO}/entry_points.txt", _entry_points()),
        )
    )
    records = [_record_line(path, data) for path, data in files]
    record_stream = io.StringIO(newline="")
    writer = csv.writer(record_stream, lineterminator="\n")
    writer.writerows(records)
    writer.writerow((f"{DIST_INFO}/RECORD", "", ""))
    files.append((f"{DIST_INFO}/RECORD", record_stream.getvalue().encode("utf-8")))
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as wheel:
        for path, data in files:
            info = zipfile.ZipInfo(path, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            wheel.writestr(info, data)
    return WHEEL_NAME


def prepare_metadata_for_build_wheel(metadata_directory, config_settings=None):
    del config_settings
    directory = Path(metadata_directory) / DIST_INFO
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "METADATA").write_bytes(_metadata())
    (directory / "WHEEL").write_bytes(_wheel_metadata())
    (directory / "entry_points.txt").write_bytes(_entry_points())
    return DIST_INFO


def get_requires_for_build_wheel(config_settings=None):
    del config_settings
    return []


def get_requires_for_build_sdist(config_settings=None):
    del config_settings
    return []


def build_sdist(sdist_directory, config_settings=None):
    del config_settings
    root = Path(__file__).resolve().parent
    filename = f"{NORMALIZED}-{VERSION}.tar.gz"
    destination = Path(sdist_directory) / filename
    selected = [
        *sorted(root.glob("*.py")),
        root / "pyproject.toml",
        root / "README.md",
        root / "requirements-fit.txt",
        root / "requirements-archive.txt",
    ]
    with open(destination, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in selected:
                    data = path.read_bytes()
                    info = tarfile.TarInfo(f"{NORMALIZED}-{VERSION}/{path.name}")
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    archive.addfile(info, io.BytesIO(data))
    return filename
