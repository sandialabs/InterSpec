#!/usr/bin/env bash
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

# run_all_benchmarks.sh — Run all §9.3 GEANT4 benchmark configurations and
# compare against our MC results.
#
# Usage (from the directory containing detector_*.gdml, our_*.csv, run_*.mac):
#   bash /path/to/run_all_benchmarks.sh
#
# Requirements:
#   - ceelo_g4val in PATH or set G4VAL variable
#   - Geant4 datasets in place (GEANT4_INSTALL_DATA was used)
#   - generate_g4_benchmarks was run first to produce .gdml, .mac, and our_*.csv files
#   - compare_results.py from tools/geant4_validation/

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
G4VAL="${G4VAL:-ceelo_g4val}"
COMPARE="${COMPARE:-python3 $(dirname "$0")/compare_results.py}"
TOLERANCE="${TOLERANCE:-2.5}"

# Datasets for each config (energy in keV, must match filenames from generate_g4_benchmarks)
declare -A CONFIGS
CONFIGS[1]="30 50 80 100 200 300 500 662 800 1000 1173 1332 2000 3000"
CONFIGS[2]="30 50 80 100 200 300 500 662 800 1000 1173 1332 2000 3000"
CONFIGS[3]="30 50 80 100 200 300 500 662 800 1000 1173 1332 2000 3000"
CONFIGS[5]="30 50 80 100 200 300 500 662 800 1000 1500"
CONFIGS[6]="100 200 300 500 662 1000 1173 1332 2000"
CONFIGS[7]="100 200 300 500 662 1000 1173 1332 2000 3000"

# ---------------------------------------------------------------------------
# Run G4 for each config × energy, collect into per-config CSV
# ---------------------------------------------------------------------------
echo "=== Running GEANT4 benchmarks ==="
echo ""

for cfg in $(echo "${!CONFIGS[@]}" | tr ' ' '\n' | sort -n); do
    gdml="detector_${cfg}.gdml"
    if [ ! -f "$gdml" ]; then
        echo "  [SKIP] Config $cfg: $gdml not found"
        continue
    fi

    # Write header for combined G4 CSV
    g4_csv="g4_${cfg}_multi.csv"
    echo "# ceelo_g4val reference data (config $cfg)" > "$g4_csv"
    echo "energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events" >> "$g4_csv"

    echo "Config $cfg ($gdml):"
    for E in ${CONFIGS[$cfg]}; do
        mac="run_${cfg}_${E}keV.mac"
        tmp_csv="g4_${cfg}_${E}keV_tmp.csv"

        if [ ! -f "$mac" ]; then
            echo "  [SKIP] $mac not found"
            continue
        fi

        printf "  E=%6d keV ... " "$E"
        "$G4VAL" "$gdml" "$mac" "$tmp_csv" 1 2>/dev/null
        # Extract the data row (skip comment and header lines)
        grep -v "^#" "$tmp_csv" | tail -n 1 >> "$g4_csv"
        echo "done"
        rm -f "$tmp_csv"
    done
    echo ""
done

# ---------------------------------------------------------------------------
# Compare our results to G4 for each config
# ---------------------------------------------------------------------------
echo "=== Comparing CeeLo vs GEANT4 ==="
echo ""

all_pass=true
for cfg in $(echo "${!CONFIGS[@]}" | tr ' ' '\n' | sort -n); do
    our_csv="our_${cfg}_multi.csv"
    g4_csv="g4_${cfg}_multi.csv"

    if [ ! -f "$our_csv" ] || [ ! -f "$g4_csv" ]; then
        echo "  [SKIP] Config $cfg: missing CSV files"
        continue
    fi

    echo "--- Config $cfg ---"
    if ! $COMPARE "$our_csv" "$g4_csv" --tolerance "$TOLERANCE"; then
        all_pass=false
    fi
    echo ""
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "========================================"
if $all_pass; then
    echo "OVERALL: ALL PASS (tolerance = ${TOLERANCE}%)"
    exit 0
else
    echo "OVERALL: SOME FAILURES (tolerance = ${TOLERANCE}%)"
    exit 1
fi
