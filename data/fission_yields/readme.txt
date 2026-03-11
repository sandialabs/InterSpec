This data retrieved from the IAEA "LiveChart of Nuclides" website at:
https://www-nds.iaea.org/livechart/

You can either search for nuclide wanted, then click the "Fission Yields" tab, 
and then download the CSV for "Independent Fission Yield", or use the API, via 
a URL like:
    https://www-nds.iaea.org/relnsd/v1/data?fields=independent_fy&parents=235u


Spontaneous-fission yields were added from ENDF/B-VIII.0 SFY files (MF=8, MT=454),
using below Python code.

Source archive used by that script:
  https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_sfy.zip

The script writes/updates the column:
  independent_spontaneous_fy

along with:
  unc_isf

Current SFY coverage from ENDF/B-VIII.0 includes:
  U238, Cm244, Cm246, Cm248, Cf250, Cf252, Es253, Fm254, Fm256


Python code used to extract SF yields:
```
#!/usr/bin/env python3
"""
Generate InterSpec fission-yield CSV files from ENDF spontaneous-yield files.

This script uses ENDF/B-VIII.0 SFY files and writes/updates files in:
  data/fission_yields/

It adds the columns:
  - independent_spontaneous_fy
  - unc_isf

to enable spontaneous-fission support in ReferenceLineInfo.
"""

from __future__ import annotations

import csv
import datetime as _dt
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import urllib.request
import zipfile

from endf_parserpy import EndfParserPy


ROOT = Path(__file__).resolve().parents[2]
SFY_DIR = ROOT / "support" / "fission_data" / "ENDF-B-VIII.0_sfy"
SFY_ZIP = ROOT / "support" / "fission_data" / "ENDF-B-VIII.0_sfy.zip"
SFY_URL = "https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_sfy.zip"
OUT_DIR = ROOT / "data" / "fission_yields"


# Requested by user for initial generation.
REQUESTED_PARENTS = ("cf252", "cm248", "cm244", "pu240", "u238")


PERIODIC = [None, *"H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()]
CSV_COLUMNS = "z_daughter,a_daughter,element_daughter,z_parent,a_parent,element_parent,daughter_level_idx,independent_thermal_fy,unc_it,independent_fast_fy,unc_if,independent_14mev_fy,unc_i14,independent_spontaneous_fy,unc_isf,Extraction_date".split(",")


def symbol_for_z(z: int) -> str:
    if z <= 0 or z >= len(PERIODIC) or PERIODIC[z] is None:
        raise ValueError(f"Unsupported atomic number Z={z}")
    return PERIODIC[z]


def parse_parent_from_za(za: float) -> Tuple[int, int]:
    za_i = int(round(za))
    z = za_i // 1000
    a = za_i - 1000 * z
    return z, a


def format_float(value: float) -> str:
    if value == 0.0:
        return "0"
    # Keep enough precision for tiny yields while avoiding gratuitous verbosity.
    return f"{value:.10g}"


def spontaneous_rows_from_endf(path: Path) -> Tuple[str, List[Dict[str, str]]]:
    parser = EndfParserPy()
    data = parser.parsefile(str(path))

    mf8 = data.get(8, {})
    mt454 = mf8.get(454)
    if not mt454:
        raise RuntimeError(f"Missing MF=8 MT=454 independent yields in {path.name}")

    z_parent, a_parent = parse_parent_from_za(mt454["ZA"])
    parent_label = f"{symbol_for_z(z_parent).lower()}{a_parent}"

    # SFY files in ENDF/B-VIII.0 have a single energy point at E=0.
    # Pick the lowest-index energy bin robustly.
    e_idx = sorted(mt454["E"].keys())[0]
    nfp = int(mt454["NFP"][e_idx])
    zafp = mt454["ZAFP"][e_idx]
    fps = mt454["FPS"][e_idx]
    yi = mt454["YI"][e_idx]
    dyi = mt454["DYI"][e_idx]

    extraction_date = _dt.date.today().isoformat()
    rows: List[Dict[str, str]] = []

    for idx in range(1, nfp + 1):
        daughter_za = int(round(zafp[idx]))
        z_daughter = daughter_za // 1000
        a_daughter = daughter_za - 1000 * z_daughter
        level_idx = int(round(fps[idx]))
        row = {
            "z_daughter": str(z_daughter),
            "a_daughter": str(a_daughter),
            "element_daughter": symbol_for_z(z_daughter),
            "z_parent": str(z_parent),
            "a_parent": str(a_parent),
            "element_parent": symbol_for_z(z_parent),
            "daughter_level_idx": str(level_idx),
            "independent_thermal_fy": "",
            "unc_it": "",
            "independent_fast_fy": "",
            "unc_if": "",
            "independent_14mev_fy": "",
            "unc_i14": "",
            "independent_spontaneous_fy": format_float(float(yi[idx])),
            "unc_isf": format_float(float(dyi[idx])),
            "Extraction_date": extraction_date,
        }
        rows.append(row)

    # Deterministic order, matching existing style as much as possible.
    rows.sort(
        key=lambda r: (
            int(r["z_daughter"]),
            int(r["a_daughter"]),
            int(r["daughter_level_idx"]),
        )
    )
    return parent_label, rows


def read_existing_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return []

        rows: List[Dict[str, str]] = []
        for raw in reader:
            # Support both legacy InterSpec files (14 columns) and regenerated
            # files with spontaneous columns (16 columns).
            if len(raw) < 14:
                continue

            row_map = {
                "z_daughter": raw[0].strip() if len(raw) > 0 else "",
                "a_daughter": raw[1].strip() if len(raw) > 1 else "",
                "element_daughter": raw[2].strip() if len(raw) > 2 else "",
                "z_parent": raw[3].strip() if len(raw) > 3 else "",
                "a_parent": raw[4].strip() if len(raw) > 4 else "",
                "element_parent": raw[5].strip() if len(raw) > 5 else "",
                "daughter_level_idx": raw[6].strip() if len(raw) > 6 else "",
                "independent_thermal_fy": raw[7].strip() if len(raw) > 7 else "",
                "unc_it": raw[8].strip() if len(raw) > 8 else "",
                "independent_fast_fy": raw[9].strip() if len(raw) > 9 else "",
                "unc_if": raw[10].strip() if len(raw) > 10 else "",
                "independent_14mev_fy": raw[11].strip() if len(raw) > 11 else "",
                "unc_i14": raw[12].strip() if len(raw) > 12 else "",
                "independent_spontaneous_fy": raw[13].strip() if len(raw) > 13 and len(raw) >= 16 else "",
                "unc_isf": raw[14].strip() if len(raw) > 14 and len(raw) >= 16 else "",
                "Extraction_date": raw[15].strip() if len(raw) > 15 and len(raw) >= 16 else (raw[13].strip() if len(raw) > 13 else ""),
            }
            # Some historic files can include an empty trailing row.
            if not row_map.get("z_daughter", "") or not row_map.get("a_daughter", ""):
                continue

            normalized = {col: row_map.get(col, "") for col in CSV_COLUMNS}
            rows.append(normalized)
    return rows


def upsert_spontaneous_column(existing: List[Dict[str, str]], spontaneous: List[Dict[str, str]]) -> List[Dict[str, str]]:
    index: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    for row in existing:
        key = (row["z_daughter"], row["a_daughter"], row["daughter_level_idx"])
        index[key] = row
        # Ensure columns exist in case file had legacy schema.
        if "independent_spontaneous_fy" not in row:
            row["independent_spontaneous_fy"] = ""
        if "unc_isf" not in row:
            row["unc_isf"] = ""

    for row in spontaneous:
        key = (row["z_daughter"], row["a_daughter"], row["daughter_level_idx"])
        if key in index:
            index[key]["independent_spontaneous_fy"] = row["independent_spontaneous_fy"]
            index[key]["unc_isf"] = row["unc_isf"]
            index[key]["Extraction_date"] = row["Extraction_date"]
        else:
            index[key] = row

    merged = list(index.values())
    merged.sort(
        key=lambda r: (
            int(r["z_daughter"]),
            int(r["a_daughter"]),
            int(r["daughter_level_idx"]),
        )
    )
    return merged


def write_rows(path: Path, rows: Iterable[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in CSV_COLUMNS})


def ensure_sfy_library() -> None:
    if SFY_DIR.is_dir() and any(SFY_DIR.glob("sfy-*_*.endf")):
        return

    SFY_ZIP.parent.mkdir(parents=True, exist_ok=True)
    if not SFY_ZIP.is_file():
        req = urllib.request.Request(
            SFY_URL,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req) as response:
            SFY_ZIP.write_bytes(response.read())

    with zipfile.ZipFile(SFY_ZIP, "r") as zf:
        zf.extractall(SFY_ZIP.parent)


def main() -> int:
    ensure_sfy_library()

    sfy_files = sorted(SFY_DIR.glob("sfy-*_*.endf"))
    if not sfy_files:
        raise RuntimeError(f"No SFY ENDF files found under {SFY_DIR}")

    available: Dict[str, Path] = {}
    extracted: Dict[str, List[Dict[str, str]]] = {}

    for sfy_file in sfy_files:
        parent_label, rows = spontaneous_rows_from_endf(sfy_file)
        available[parent_label] = sfy_file
        extracted[parent_label] = rows

    print("Available spontaneous-yield parents from ENDF-B-VIII.0 SFY:")
    for k in sorted(available.keys()):
        print(f"  - {k} ({available[k].name})")

    for parent in REQUESTED_PARENTS:
        out_path = OUT_DIR / f"{parent}_independent_fy.csv"
        if parent not in extracted:
            print(f"WARNING: no spontaneous-yield source found for {parent}; not writing {out_path.name}")
            continue

        sf_rows = extracted[parent]
        if out_path.exists():
            merged_rows = upsert_spontaneous_column(read_existing_rows(out_path), sf_rows)
            write_rows(out_path, merged_rows)
            print(f"Updated {out_path.name}: merged spontaneous yields into existing file ({len(merged_rows)} rows)")
        else:
            write_rows(out_path, sf_rows)
            print(f"Created {out_path.name}: spontaneous yields ({len(sf_rows)} rows)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```