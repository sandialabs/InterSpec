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

"""Monte Carlo bank cells: one CSV per cell, `cells/<identity>_<kind>.csv.gz`.

The campaign runners write every cell gzip-compressed (15x smaller; parsing
costs 2-7% more).  Every reader goes through `cell_files`, `cell_path`,
`open_cell` and `cell_id`, which also accept uncompressed `.csv` cells, so
banks generated before compression and in-progress work areas stay readable.
A cell's identity -- which seeds the fitters' train/validation hash splits --
and the order cells are read in are those of the uncompressed names, whatever
the compression.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

GZIP_LEVEL = 6


def _plain(path: Path) -> str:
    name = str(path)
    return name[:-3] if name.endswith(".gz") else name


def cell_id(path: Path, kind: str = "ceelo") -> str:
    """`<identity>` of `<identity>_<kind>.csv[.gz]`."""
    name = Path(_plain(path)).name
    suffix = f"_{kind}.csv"
    if not name.endswith(suffix):
        raise ValueError(f"{path} is not a {kind} cell file")
    return name[: -len(suffix)]


def cell_files(directory: Path, kind: str = "ceelo", recursive: bool = False) -> list[Path]:
    """Every `<identity>_<kind>.csv[.gz]` under `directory`, in the order the
    uncompressed paths sort.  A cell left in both forms (a rewrite interrupted
    before it removed the uncompressed copy) is read from its compressed file,
    which is the newer."""
    pattern = f"**/*_{kind}.csv" if recursive else f"*_{kind}.csv"
    chosen: dict[Path, Path] = {path: path for path in directory.glob(pattern)}
    for path in directory.glob(pattern + ".gz"):
        chosen[Path(_plain(path))] = path
    return [chosen[plain] for plain in sorted(chosen)]


def cell_path(directory: Path, identity: str, kind: str = "ceelo") -> Path | None:
    """The existing file of one cell, or None."""
    for suffix in (".csv.gz", ".csv"):
        path = directory / f"{identity}_{kind}{suffix}"
        if path.exists():
            return path
    return None


def companion(path: Path, kind_from: str, kind_to: str) -> Path | None:
    """The `<kind_to>` file of the cell `path` (a `<kind_from>` file) holds."""
    return cell_path(path.parent, cell_id(path, kind_from), kind_to)


def open_cell(path: Path):
    """A text stream for csv.reader / csv.DictReader."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def written_path(directory: Path, identity: str, kind: str = "ceelo") -> Path:
    """Where a writer puts a new cell: always compressed."""
    return directory / f"{identity}_{kind}.csv.gz"


def existing(path: Path) -> Path | None:
    """The finished cell a writer's `path` names, compressed or not, or None:
    a runner's resume check, so a re-run never redoes an uncompressed cell."""
    if path.exists():
        return path
    plain = Path(_plain(path))
    return plain if plain != path and plain.exists() else None


def write_cell(path: Path, text: str) -> None:
    """Compress `text` into `path` in one step (temporary file + rename), so an
    interrupted run never leaves a truncated cell that a resume would accept.
    The gzip header carries no name or time, so equal text gives equal bytes."""
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw,
                           compresslevel=GZIP_LEVEL, mtime=0) as stream:
            stream.write(text.encode())
    os.replace(temporary, path)
    plain = Path(_plain(path))
    if plain != path and plain.exists():
        # A forced re-run replaces an uncompressed cell; keep one form only.
        plain.unlink()
