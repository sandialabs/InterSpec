#!/bin/bash
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

# Run all diagnostic comparisons for FEP discrepancy isolation.
#
# Prerequisites:
#   1. Build the CeeLo project:  cmake --build build -j8
#   2. Build the G4 harness:
#        source <your-geant4-install>/bin/geant4.sh
#        cmake -B <build-dir> -S tools/geant4_validation \
#              -DGeant4_DIR=<your-geant4-install>/lib/cmake/Geant4
#        cmake --build <build-dir> -j4
#        # then either put ceelo_g4val on PATH or set G4VAL=<build-dir>/ceelo_g4val
#
# Usage:
#   cd /path/to/mc_det_eff_plan
#   source <your-geant4-install>/bin/geant4.sh
#   bash tools/geant4_validation/run_diagnostics.sh

set -e

PROJECT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
HARNESS="${G4VAL:-ceelo_g4val}"   # override: G4VAL=/path/to/ceelo_g4val
# Override with DIAG_SPECTRUM=/path/to/diag_spectrum (a Release build is strongly
# recommended for the 10M-event >=2 MeV rows).
DIAG_SPECTRUM="${DIAG_SPECTRUM:-$PROJECT_DIR/build/examples/diag_spectrum}"
GDML="$PROJECT_DIR/build/examples/detector_1.gdml"
MACROS="$PROJECT_DIR/tools/geant4_validation/macros"
OUTDIR="$PROJECT_DIR/build/diagnostics"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== Diagnostic runs for FEP discrepancy isolation ==="
echo "Output directory: $OUTDIR"
echo ""

# =============================================
# Part 1: Our MC — energy deposit spectra
# =============================================
echo "--- Our MC: energy deposit spectra ---"

# FEP efficiency falls to ~3e-3 by 3 MeV, so the >=2 MeV rows need far more
# events than the sub-MeV ones to resolve the brems/MSC shares.
for E in 662 1173 2614 3000; do
    if [ "$E" -ge 2000 ]; then N_EVENTS=10000000; else N_EVENTS=500000; fi
    for MODE in full no_moliere no_brems no_both; do
        echo "  Our MC: ${E} keV, mode=${MODE}..."
        "$DIAG_SPECTRUM" "$E" "$N_EVENTS" "$MODE"
        echo ""
    done
done

# =============================================
# Part 2: GEANT4 — baseline + process-disabled
# =============================================
if [ ! -x "$HARNESS" ]; then
    echo "WARNING: G4 harness not found at $HARNESS"
    echo "Build it first (see instructions above), then re-run."
    echo "Skipping G4 runs."
    exit 0
fi

if [ ! -f "$GDML" ]; then
    echo "WARNING: GDML file not found at $GDML"
    echo "Run generate_g4_benchmarks first."
    exit 1
fi

echo ""
echo "--- GEANT4: baseline + process isolation runs ---"

# noRayl exists only for 662/1173 (Rayleigh is negligible >=2 MeV); the loop
# below skips absent macros.
VARIANTS="baseline noEBrem noMSC noBoth noRayl"

for E in 662 1173 2614 3000; do
    for VAR in $VARIANTS; do
        MAC="$MACROS/diag_${E}keV_${VAR}.mac"
        if [ ! -f "$MAC" ]; then
            echo "  SKIP: $MAC not found"
            continue
        fi
        OUT="g4_diag_${E}_${VAR}.csv"
        echo "  G4: ${E} keV, variant=${VAR}..."
        "$HARNESS" "$GDML" "$MAC" "$OUT" --histogram --cone-bias
        echo ""
    done
done

echo ""
echo "=== All diagnostic runs complete ==="
echo "Results in: $OUTDIR"
echo ""
echo "Compare histograms:"
echo "  Our MC:  diag_spectrum_*keV_*.csv"
echo "  GEANT4:  g4_diag_*_*_histogram.csv"
echo ""
echo "Compare efficiencies:"
echo "  Our MC:  diag_spectrum_*keV_*_summary.csv"
echo "  GEANT4:  g4_diag_*_*.csv (FEP/total efficiency)"
