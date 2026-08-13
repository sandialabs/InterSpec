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

"""
compare_results.py — Compare CeeLo simulation against GEANT4 reference data.

Usage:
    python compare_results.py <ceelo_csv> <geant4_csv> [--plot]

Both CSV files must have the format produced by ceelo_g4val:
    energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events

The script prints a table of energy-by-energy comparisons and optionally
produces a plot.

Exit code 0 = all within 2.5%; non-zero = at least one deviation > 2.5%.
"""

import sys
import csv
import math
import argparse


def read_csv(path):
    """Read a results CSV, skipping comment lines starting with #."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(
            (line for line in f if not line.startswith("#")),
            fieldnames=["energy_keV", "fep_efficiency", "fep_uncertainty",
                        "total_efficiency", "total_uncertainty", "num_events"],
        )
        for row in reader:
            if row["energy_keV"] == "energy_keV":
                continue  # skip header row if present
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def pct_diff(a, b):
    """Relative difference (a - b) / b in percent."""
    if b == 0.0:
        return float("inf") if a != 0.0 else 0.0
    return (a - b) / b * 100.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ceelo_csv", help="CeeLo simulation results CSV")
    parser.add_argument("geant4_csv",   help="GEANT4 reference results CSV")
    parser.add_argument("--plot", action="store_true", help="Generate comparison plot")
    parser.add_argument("--tolerance", type=float, default=2.5,
                        help="Acceptance tolerance in %% (default 2.5)")
    args = parser.parse_args()

    our  = read_csv(args.ceelo_csv)
    ref  = read_csv(args.geant4_csv)

    # Match rows by energy.
    ref_dict = {round(r["energy_keV"], 2): r for r in ref}

    header = (
        f"{'Energy':>10}  {'ours FEP':>12}  {'G4 FEP':>12}  {'Δ FEP %':>9}  "
        f"{'ours Tot':>12}  {'G4 Tot':>12}  {'Δ Tot %':>9}  {'Status':>8}"
    )
    print(header)
    print("-" * len(header))

    all_pass = True
    for row in our:
        e = round(row["energy_keV"], 2)
        g4 = ref_dict.get(e)
        if g4 is None:
            print(f"  {e:>8} keV  [no G4 reference at this energy]")
            continue

        dfep = pct_diff(row["fep_efficiency"], g4["fep_efficiency"])
        dtot = pct_diff(row["total_efficiency"], g4["total_efficiency"])
        ok   = abs(dfep) <= args.tolerance and abs(dtot) <= args.tolerance
        if not ok:
            all_pass = False

        print(
            f"{e:>10.1f}  "
            f"{row['fep_efficiency']:>12.5f}  {g4['fep_efficiency']:>12.5f}  "
            f"{dfep:>+9.2f}%  "
            f"{row['total_efficiency']:>12.5f}  {g4['total_efficiency']:>12.5f}  "
            f"{dtot:>+9.2f}%  "
            f"{'PASS' if ok else 'FAIL':>8}"
        )

    print()
    if all_pass:
        print(f"ALL PASS (tolerance = {args.tolerance}%)")
    else:
        print(f"SOME FAILURES (tolerance = {args.tolerance}%)")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            energies   = [r["energy_keV"] for r in our]
            fep_ratio  = []
            tot_ratio  = []
            for row in our:
                e  = round(row["energy_keV"], 2)
                g4 = ref_dict.get(e)
                if g4:
                    fep_ratio.append(row["fep_efficiency"] / g4["fep_efficiency"] if g4["fep_efficiency"] > 0 else None)
                    tot_ratio.append(row["total_efficiency"] / g4["total_efficiency"] if g4["total_efficiency"] > 0 else None)

            fig, ax = plt.subplots(figsize=(9, 5))
            ax.plot(energies, fep_ratio, "o-", label="ε_FEP ratio (ours/G4)")
            ax.plot(energies, tot_ratio, "s--", label="ε_total ratio (ours/G4)")
            ax.axhline(1.0, color="k", linewidth=0.8)
            ax.axhline(1.025, color="r", linewidth=0.6, linestyle=":")
            ax.axhline(0.975, color="r", linewidth=0.6, linestyle=":")
            ax.set_xlabel("Energy (keV)")
            ax.set_ylabel("Ratio (ours / GEANT4)")
            ax.set_title("CeeLo vs GEANT4 efficiency comparison")
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig("comparison.png", dpi=150)
            print("Plot saved to comparison.png")
        except ImportError:
            print("matplotlib not available; skipping plot.")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
