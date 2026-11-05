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

# Source-escape bake: CeeLo and GEANT4 truth banks -> the four fitted parts of
# the source-escape payload (secondary surrogate, elastic-line scale,
# photonuclear surrogate, multiple-scatter surrogate) -> assembly -> scoring on
# held-out banks and the 150-cell acceptance grid -> gates -> publish.
#
#   tools/full_spectrum/bake/source_escape.sh            # production
#   SMOKE=1 tools/full_spectrum/bake/source_escape.sh    # smoke run
#
# Optional parts, both GEANT4 runs (source <geant4>/bin/geant4.sh first):
#   PHOTONUCLEAR=1   also bake the photonuclear term: a GEANT4 truth bank at
#                    8-12 MeV (~450 CPU-hours) and the natural-element
#                    photonuclear cross sections.  Without it the payload has
#                    no photonuclear term, which evaluation leaves off by
#                    default anyway (it matters only above ~8 MeV).
#   GEANT4_CHECK=1   also run the GEANT4 reference leg of the 150-cell
#                    acceptance grid (a few CPU-hours) and gate the payload
#                    against it; without it the grid is scored against CeeLo.
#
# Settings: see common.sh.  Resumable: re-run the same
# command after an interruption; bank tools skip finished cells.  A step whose
# inputs changed since it ran (a re-run bank or fit) runs again automatically,
# and so does everything downstream of it, so REDO=<step> never leaves a stale
# payload or stale scores behind.
. "$(dirname "$0")/common.sh"
family source_escape

PHOTONUCLEAR=${PHOTONUCLEAR:-0}
GEANT4_CHECK=${GEANT4_CHECK:-0}
if [ "$PHOTONUCLEAR" = 1 ] || [ "$GEANT4_CHECK" = 1 ]; then
    need_geant4
fi
stamp "photonuclear: $([ "$PHOTONUCLEAR" = 1 ] && echo baked || echo 'not baked (PHOTONUCLEAR=1 bakes it)'); GEANT4 acceptance leg: $([ "$GEANT4_CHECK" = 1 ] && echo on || echo 'off (GEANT4_CHECK=1 runs it)')"
need_bin ceelo_source_escape_campaign ceelo_source_model_evaluate \
    ceelo_source_cross_section_features

BANKS=$FW/banks
HELPER=$BIN/ceelo_source_cross_section_features
EVALUATOR=$BIN/ceelo_source_model_evaluate
CAMPAIGN=$BIN/ceelo_source_escape_campaign
PAYLOAD=$FW/payload/source_escape_payload.txt

if [ "$SMOKE" = 1 ]; then
    # Every step on a handful of cells and 1e3-1e4 histories each; the
    # numbers mean nothing, the code paths are the production ones.
    SECONDARY_BANK="--events 20000"
    for e in 661.657 3000; do
        for z in 26 82; do
            for t in 0.2 0.5 1 2 4 8; do
                SECONDARY_BANK="$SECONDARY_BANK --cell $e,$z,$t"
            done
        done
        SECONDARY_BANK="$SECONDARY_BANK --cell $e,13,1 --cell $e,13,4"
    done
    OFFGRID_BANK="--events 20000 --cell 661.657,26,0.7 --cell 661.657,26,3
        --cell 3000,82,0.7 --cell 3000,82,3 --cell 3000,82,6"
    HYDROGEN_BANK="--events 20000 --cell 661.657,26,1,0.1119
        --cell 661.657,26,4,0.1119 --cell 3000,26,1,0.25
        --cell 3000,26,4,0.25 --cell 661.657,35,1,0.1119"
    # The fitter keeps only cells with >= 200 histories of three or more
    # Compton scatters, so these cells are thick and low-Z enough for that.
    MS_TRAINING="--energies 300 1000 --z-values 6 13 26 50
        --tau-values 1 2 4 8 --events 100000 --bins 64"
    MS_REFINE="--max-events 200000"
    MS_TEST="--energies 500 --z-values 10 20 --tau-values 3 6
        --events 100000 --bins 64"
    PN_XS="--energy-nodes 16"
    PN_BANK="--events 3000 --energies-kev 12000 --zs 82 47
        --taus 0.2 1 4"
    PN_OFFGRID="--events 3000 --energies-kev 12000 --zs 82 --taus 2"
    ACCEPTANCE="--cell 661,26,5 --cell 3000,82,1 --ceelo-events 20000
        --geant4-events 2000 --bins 128"
else
    SECONDARY_BANK="--events 1000000"
    OFFGRID_BANK="--offgrid --events 4000000"
    HYDROGEN_BANK="--hydrogen --events 1000000"
    MS_TRAINING="--grid training --events 100000"
    MS_REFINE=""
    MS_TEST="--grid test --events 200000"
    PN_XS="--energy-nodes 128"
    PN_BANK=""
    PN_OFFGRID="--offgrid"
    # 20M histories per cell: at 500k the grid's opaque cells carry 10-45%
    # truth noise, and the p95 they set measures that noise, not the payload.
    ACCEPTANCE="--acceptance-grid --ceelo-events 20000000"
fi

# stale STEP OUTPUT INPUT...: forget that STEP completed when OUTPUT is missing
# or any INPUT (a file, or an upstream step's done marker) is newer than it.
stale() {
    local name=$1 output=$2 input
    shift 2
    [ -f "$DONE/$name" ] || return 0
    if [ ! -e "$output" ]; then
        rm -f "$DONE/$name"
        return 0
    fi
    for input in "$@"; do
        if [ -e "$input" ] && [ "$input" -nt "$output" ]; then
            stamp "$name is out of date ($input changed); it runs again"
            rm -f "$DONE/$name"
            return 0
        fi
    done
}

# mode_file NAME VALUE: a file holding VALUE whose timestamp changes only when
# VALUE does, so a step can list a setting among its inputs.
mode_file() {
    local file=$FW/$1.mode
    if [ ! -f "$file" ] || [ "$(cat "$file")" != "$2" ]; then
        echo "$2" > "$file"
    fi
    echo "$file"
}

json_value() {
    $PY -c 'import json, sys; print(repr(json.load(open(sys.argv[1]))[sys.argv[2]]))' "$1" "$2"
}

# ---- C++/Python parity first: a disagreement wastes everything downstream.
lockstep() {
    $PY $TOOLS/check_base_lockstep.py --build "$BUILD"
    $PY $TOOLS/check_feature_lockstep.py --build "$BUILD"
    $PY $TOOLS/check_basis_lockstep.py
}

# ---- CeeLo truth banks for the secondary fit (Rayleigh, Compton Doppler,
# shield fluorescence and impact ionization on).
secondary_bank() {
    # shellcheck disable=SC2086
    $PY $TOOLS/run_source_secondary_campaign.py --output "$BANKS/$1" \
        --driver "$CAMPAIGN" --workers "$JOBS" $2
}

kernel_span() {
    $PY $TOOLS/check_kernel_span.py --bank "$BANKS/secondary" \
        --bank "$BANKS/hydrogen" --helper "$HELPER" \
        --output "$FW/kernel_span.json"
}

fit_elastic() {
    $PY $TOOLS/fit_elastic_escape.py "$BANKS/secondary" "$BANKS/hydrogen" \
        --cross-section-helper "$HELPER" --output "$FW/elastic/elastic_escape.json"
}

fit_secondary() {
    $PY $TOOLS/fit_source_secondaries.py "$BANKS/secondary" "$BANKS/hydrogen" \
        --offgrid-campaign "$BANKS/offgrid" --output "$FW/secondary" \
        --plot-worst 0 --scan-bandwidths --tail-decay-kappa 0.8 \
        --local-linear-ridge 0.05 --cross-section-helper "$HELPER"
}

# ---- GEANT4: natural-element photonuclear attenuation and the truth banks.
photonuclear_xs() {
    mkdir -p "$FW/photonuclear"
    # shellcheck disable=SC2086
    "$BIN/g4_photonuclear_xs" $PN_XS --output "$FW/photonuclear/xs.csv.tmp"
    mv "$FW/photonuclear/xs.csv.tmp" "$FW/photonuclear/photonuclear_xs.csv"
}

photonuclear_bank() {
    # shellcheck disable=SC2086
    $PY $TOOLS/run_photonuclear_campaign.py --output "$BANKS/$1" \
        --harness "$BIN/g4_source_escape_campaign" --workers "$JOBS" $2
}

fit_photonuclear() {
    $PY $TOOLS/fit_photonuclear.py --bank "$BANKS/photonuclear" \
        --offgrid "$BANKS/photonuclear_offgrid" \
        --photonuclear-xs "$FW/photonuclear/photonuclear_xs.csv" \
        --output "$FW/photonuclear"
}

# ---- Multiple scatter: each cell is a CeeLo truth plus the mechanistic model
# evaluated at the same cell (the fitter needs both).
ms_bank() {
    # shellcheck disable=SC2086
    $PY $TOOLS/validate_multiple_scatter.py run --output "$BANKS/$1" \
        --workers "$JOBS" --source-driver "$CAMPAIGN" \
        --model-driver "$EVALUATOR" $2
}

refine_ms() {
    # shellcheck disable=SC2086
    $PY $TOOLS/refine_multiple_scatter_cells.py "$BANKS/ms_training" \
        --workers "$JOBS" --source-driver "$CAMPAIGN" \
        --model-driver "$EVALUATOR" $MS_REFINE
}

fit_multiple_scatter() {
    $PY $TOOLS/fit_multiple_scatter_surrogate.py "$BANKS/ms_training" \
        --test-campaign "$BANKS/ms_test" --output "$FW/multiple_scatter" \
        --cross-section-helper "$HELPER"
}

# ---- Assembly: one payload from the four fits, then a round trip through
# the serialized form must reproduce the assembly-time spectrum exactly.
assemble() {
    local scale check="--energy-kev 3000 --z 82 --areal-density 5 --bins 64"
    local photonuclear=""
    scale=$(json_value "$FW/elastic/elastic_escape.json" elastic_escape_scale)
    if [ "$PHOTONUCLEAR" = 1 ]; then
        photonuclear="--photonuclear-xs $FW/photonuclear/photonuclear_xs.csv
            --photonuclear-surrogate $FW/photonuclear/photonuclear_surrogate.txt"
    fi
    mkdir -p "$FW/payload"
    # shellcheck disable=SC2086
    "$EVALUATOR" \
        --multiple-surrogate "$FW/multiple_scatter/surrogate.txt" \
        --secondary-surrogate "$FW/secondary/secondary_surrogate.txt" \
        $photonuclear \
        --elastic-escape-scale "$scale" --elastic-compton-damping 1.0 \
        --write-payload "$PAYLOAD.tmp" $check > "$FW/payload/assembled.csv"
    # shellcheck disable=SC2086
    "$EVALUATOR" --payload "$PAYLOAD.tmp" $check > "$FW/payload/reloaded.csv"
    cmp "$FW/payload/assembled.csv" "$FW/payload/reloaded.csv" \
        || die "the serialized payload does not reproduce the assembled model"
    mv "$PAYLOAD.tmp" "$PAYLOAD"
    echo "payload: $PAYLOAD (elastic_escape_scale $scale)"
    shasum -a 256 "$PAYLOAD"
}

# ---- Scoring on data no fit selected on.
# acceptance_bank [--no-geant4]: the grid's CeeLo legs (and GEANT4 legs
# unless --no-geant4); finished legs are skipped.
acceptance_bank() {
    # shellcheck disable=SC2086
    $PY $TOOLS/run_source_anchor.py --output "$BANKS/acceptance" \
        --workers "$JOBS" --ceelo-driver "$CAMPAIGN" \
        --geant4-driver "$BIN/g4_source_escape_campaign" $ACCEPTANCE "$@"
}

score_offgrid() {
    $PY $TOOLS/score_offgrid_holdout.py --offgrid "$BANKS/offgrid" \
        --payload "$PAYLOAD" --evaluator "$EVALUATOR" --workers "$JOBS" \
        --output "$FW/scores/offgrid.json"
}

score_deep_tau_tail() {
    $PY $TOOLS/check_deep_tau_tail.py "$PAYLOAD" --evaluator "$EVALUATOR" \
        --output "$FW/scores/deep_tau_tail.json"
}

score_acceptance() {
    $PY $TOOLS/score_acceptance_grid.py --grid "$BANKS/acceptance" \
        --payload "$PAYLOAD" --evaluator "$EVALUATOR" --workers "$JOBS" \
        --output "$FW/scores/acceptance.json"
}

gates() {
    local optional=""
    [ "$PHOTONUCLEAR" = 1 ] \
        && optional="--photonuclear $FW/photonuclear/report.json"
    [ "$GEANT4_CHECK" = 1 ] || optional="$optional --without-geant4"
    # shellcheck disable=SC2086
    $PY $TOOLS/check_source_escape_gates.py \
        --offgrid "$FW/scores/offgrid.json" \
        --deep-tau-tail "$FW/scores/deep_tau_tail.json" \
        --multiple-scatter "$FW/multiple_scatter/report.json" \
        --kernel-span "$FW/kernel_span.json" \
        --acceptance "$FW/scores/acceptance.json" \
        $optional --build "$BUILD" --output "$FW/gates.json"
}

gate lockstep lockstep

step bank_secondary secondary_bank secondary "$SECONDARY_BANK"
step bank_offgrid secondary_bank offgrid "$OFFGRID_BANK"
step bank_hydrogen secondary_bank hydrogen "$HYDROGEN_BANK"
# Accept the mechanistic bases before spending the fit on them.
gate kernel_span kernel_span

stale fit_elastic "$FW/elastic/elastic_escape.json" \
    "$DONE/bank_secondary" "$DONE/bank_hydrogen"
step fit_elastic fit_elastic
stale fit_secondary "$FW/secondary/secondary_surrogate.txt" \
    "$DONE/bank_secondary" "$DONE/bank_offgrid" "$DONE/bank_hydrogen"
step fit_secondary fit_secondary

PHOTONUCLEAR_INPUTS=$(mode_file photonuclear "$PHOTONUCLEAR")
if [ "$PHOTONUCLEAR" = 1 ]; then
    step photonuclear_xs photonuclear_xs
    step bank_photonuclear photonuclear_bank photonuclear "$PN_BANK"
    step bank_photonuclear_offgrid photonuclear_bank photonuclear_offgrid \
        "$PN_OFFGRID"
    stale fit_photonuclear "$FW/photonuclear/photonuclear_surrogate.txt" \
        "$DONE/photonuclear_xs" "$DONE/bank_photonuclear" \
        "$DONE/bank_photonuclear_offgrid"
    step fit_photonuclear fit_photonuclear
    PHOTONUCLEAR_INPUTS="$PHOTONUCLEAR_INPUTS
        $FW/photonuclear/photonuclear_surrogate.txt
        $FW/photonuclear/photonuclear_xs.csv"
fi

step bank_ms_training ms_bank ms_training "$MS_TRAINING"
stale refine_ms_training "$DONE/refine_ms_training" "$DONE/bank_ms_training"
step refine_ms_training refine_ms
step bank_ms_test ms_bank ms_test "$MS_TEST"
stale fit_multiple_scatter "$FW/multiple_scatter/surrogate.txt" \
    "$DONE/refine_ms_training" "$DONE/bank_ms_test"
step fit_multiple_scatter fit_multiple_scatter

# shellcheck disable=SC2086
stale assemble "$PAYLOAD" "$FW/elastic/elastic_escape.json" \
    "$FW/secondary/secondary_surrogate.txt" \
    "$FW/multiple_scatter/surrogate.txt" $PHOTONUCLEAR_INPUTS
step assemble assemble

step bank_acceptance acceptance_bank --no-geant4
if [ "$GEANT4_CHECK" = 1 ]; then
    step bank_acceptance_geant4 acceptance_bank
fi
stale score_offgrid "$FW/scores/offgrid.json" "$PAYLOAD" "$DONE/bank_offgrid"
step score_offgrid score_offgrid
stale score_deep_tau_tail "$FW/scores/deep_tau_tail.json" "$PAYLOAD"
step score_deep_tau_tail score_deep_tau_tail
stale score_acceptance "$FW/scores/acceptance.json" "$PAYLOAD" \
    "$DONE/bank_acceptance" "$DONE/bank_acceptance_geant4" \
    "$(mode_file geant4_check "$GEANT4_CHECK")"
step score_acceptance score_acceptance
gate gates gates

# A production run only gets here with every gate passed.
rm -rf "$FW/publish"
mkdir -p "$FW/publish"
cp "$PAYLOAD" "$FW/gates.json" "$FW/publish/"
cp "$FW/secondary/report.json" "$FW/publish/secondary_report.json"
cp "$FW/multiple_scatter/report.json" "$FW/publish/multiple_scatter_report.json"
if [ "$PHOTONUCLEAR" = 1 ]; then
    cp "$FW/photonuclear/report.json" "$FW/publish/photonuclear_report.json"
fi
cp "$FW/elastic/elastic_escape.json" "$FW/kernel_span.json" "$FW/publish/"
cp "$FW/scores/offgrid.json" "$FW/publish/offgrid_holdout.json"
cp "$FW/scores/deep_tau_tail.json" "$FW/publish/deep_tau_tail.json"
cp "$FW/scores/acceptance.json" "$FW/publish/acceptance_grid.json"
publish "$FW"/publish/* source_escape
stamp "SOURCE_ESCAPE_BAKE_DONE"
