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

# HPGe detector-response bake: truth banks -> fit -> bandwidth-rescaled
# candidate -> audit scoring -> arrival scoring -> direct-MC spot checks ->
# gates -> publish.
#
#   tools/full_spectrum/bake/detector.sh                        # hpge_coax
#   DETECTOR=detective_x tools/full_spectrum/bake/detector.sh   # another one
#   SMOKE=1 tools/full_spectrum/bake/detector.sh    # smoke run (minutes)
#   PRECISION=fast tools/full_spectrum/bake/detector.sh   # 1/8 the events
#
# DETECTOR names the detector: bake/config/$DETECTOR.geometry if that file
# exists (a CEELO_GEOMETRY file; tools/full_spectrum/drf_geometry.py writes one
# from an InterSpec DRF), else a campaign built-in (hpge_coax).  Each detector
# has its own work area ($WORK/detector_$DETECTOR), seeds and published
# payload ($RESULTS/${DETECTOR}_detector_response/).
#
# Energies below 100 keV are written for hpge_coax (passive-layer landmark
# 55 keV) and scaled to the detector's own landmark, and the held-out banks
# (audits, arrival) are kept 2% off every training energy
# (detector_bank_energies.py).  The placed recipe is recorded in
# $FW/recipe.txt; a work area baked with another recipe is refused unless REDO
# names the steps to run again.
#
# Other settings: see common.sh.  Resumable: re-run the same command after an
# interruption.  The banks are the expensive part (about 800 cells, 0.5-1.1 G
# events); every later step is minutes, so a failed gate means "fix the fitter
# and REDO=\"fit candidate finite_distance audit_a audit_b arrival\"", never
# "re-bake the banks".
. "$(dirname "$0")/common.sh"
DETECTOR=${DETECTOR:-hpge_coax}
precision_tiers
family "detector_$DETECTOR"
DETECTOR_OPTS="--detector $DETECTOR"
if [ -f "$TOOLS/bake/config/$DETECTOR.geometry" ]; then
    DETECTOR_OPTS="$DETECTOR_OPTS --geometry $TOOLS/bake/config/$DETECTOR.geometry"
fi
stamp "detector $DETECTOR: $DETECTOR_OPTS"

need_bin ceelo_detector_response_campaign ceelo_detector_model_evaluate \
    ceelo_rate_surface_anchor ceelo_detector_finite_distance
export CEELO_RATE_SURFACE_ANCHOR=$BIN/ceelo_rate_surface_anchor

# The passive-layer landmark places the low-energy banks.
# shellcheck disable=SC2086
"$BIN/ceelo_detector_response_campaign" $DETECTOR_OPTS --print-geometry \
    > "$FW/geometry.txt"
LANDMARK=$("$BIN/ceelo_rate_surface_anchor" --geometry "$FW/geometry.txt" \
    --info | $PY -c 'import json, sys; print(json.load(sys.stdin)["passive_landmark_keV"])')
stamp "passive-layer landmark $LANDMARK keV"

BANKS=$FW/banks
RUN="$PY $TOOLS/run_detector_response_campaign.py --executable $BIN/ceelo_detector_response_campaign $DETECTOR_OPTS"
CELL_FLAGS="--bins 256 --jobs $JOBS --detector-brems --detector-fluorescence --detector-spectrum-fidelity"
# Basis sizes per continuum component; fitting and re-evaluation must agree.
BASIS="--photoelectric-basis-size 16 --single-compton-basis-size 48
  --multiple-compton-basis-size 96 --pair-no-compton-basis-size 96
  --pair-single-compton-basis-size 64 --pair-multiple-compton-basis-size 96"

# bank NAME ENERGIES ANGLES MIN_EVENTS MAX_EVENTS SALT
# ANGLES is a list of beam cosines, or incidence=CLASS,... for an arrival
# bank (isotropic over each class's band of cosines).  The seed of every cell
# derives from $DETECTOR-SALT and the cell's (energy, angle or class), so a
# bank is reproducible cell by cell and banks with different salts (or
# detectors) are independent.  The events are the final tier's, scaled here.
bank() {
    local directions
    case $3 in
        incidence=*) directions="--incidences ${3#incidence=}" ;;
        *) directions="--angles $3" ;;
    esac
    # shellcheck disable=SC2086
    $RUN "$BANKS/$1" $CELL_FLAGS --energies "$2" $directions \
        --min-events "$(events "$4")" --max-events "$(events "$5")" \
        --seed-salt "$DETECTOR-$6"
}

if [ "$SMOKE" = 1 ]; then
    # Same banks and steps, a handful of cells each, 20-50k events.
    TRAIN_E=20,59.54,122,356.013,661.657,1332.492,2614.511,6000,12000
    DENSE_E=30,450,3200
    DENSE_A_E=100,1450
    HIGHSTAT_E=2000,9000
    LOWE_E=37.5,57.5
    VALIDATION_E=70,800,7500
    AUDIT_A_E=52.5,1000
    AUDIT_B_E=45,1500
    ARRIVAL_E=56.3,340
else
    TRAIN_E=10,20,40,59.54,80,122,200,356.013,511,661.657,1000,1023,1332.492,2000,2614.511,4000,6000,9000,12000
    DENSE_E=30,70,150,280,450,800,1150,1700,3200,7500
    DENSE_A_E=50,100,175,320,600,900,1450,2400,5000,10000,12000
    HIGHSTAT_E=1332.492,2000,2614.511,4000,5000,6000,9000,12000
    # ~5 keV spacing through the passive-layer turn-on, 0.64-1.32 landmarks.
    LOWE_E=35,37.5,47.5,50,57.5,60,65,72.5
    VALIDATION_E=30,70,150,275,450,800,1173.228,1600,3200,7500,11000
    # Held out: energies and angles off every training grid, including
    # near-grazing incidence (|cos| = 0.93-0.97).
    AUDIT_A_E=42.5,52.5,62.5,75,95,500,1000,2000,5000,10000
    AUDIT_B_E=45,55,70,85,120,300,700,1500,3000,7000
    # Held out: the response a scattered field is folded with, from 0.75
    # landmarks to 850 keV.
    ARRIVAL_E=41,48.5,56.3,67.5,77.5,90,110,135,162.5,190,237.5,340,400,640,850
fi
PLACED=$($PY $TOOLS/detector_bank_energies.py --landmark "$LANDMARK" \
    --training E_TRAIN="$TRAIN_E" --training E_DENSE="$DENSE_E" \
    --training E_DENSE_A="$DENSE_A_E" --training E_HIGHSTAT="$HIGHSTAT_E" \
    --training-scaled E_LOWE="$LOWE_E" \
    --held-out-scaled E_AUDIT_A="$AUDIT_A_E" \
    --held-out-scaled E_AUDIT_B="$AUDIT_B_E" \
    --held-out-scaled E_ARRIVAL="$ARRIVAL_E") \
    || die "bank energy placement failed"
eval "$PLACED"

# One bank per line: NAME ENERGIES ANGLES MIN_EVENTS MAX_EVENTS SALT.
if [ "$SMOKE" = 1 ]; then
    BANK_SPECS="training $E_TRAIN -1,0,1 20000 20000 training
dense_energy $E_DENSE -0.5,0.5 20000 20000 dense-energy
dense_angles $E_DENSE_A -0.75,-0.25,0.25,0.75 20000 20000 dense-angles
highstat $E_HIGHSTAT -0.5,0.5 50000 50000 highstat
lowE $E_LOWE -1,0,1 20000 20000 lowE
validation $VALIDATION_E -0.5,0.5 20000 20000 validation
audit_a $E_AUDIT_A -0.85,0.85 20000 20000 audit-a
audit_b $E_AUDIT_B -0.72,0.72 20000 20000 audit-b
arrival $E_ARRIVAL incidence=iso,side 20000 20000 arrival"
    SPOT_EVENTS=100000
    FD_OPTS="--energies 60,662,2614 --events 20000"
else
    # highstat is the slow bank: 40 cells at a fixed 3M events.
    BANK_SPECS="training $E_TRAIN -1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1 200000 500000 training
dense_energy $E_DENSE -0.875,-0.625,-0.375,-0.125,0.125,0.375,0.625,0.875 300000 600000 dense-energy
dense_angles $E_DENSE_A -0.9375,-0.8125,-0.6875,-0.5625,-0.4375,-0.3125,-0.1875,-0.0625,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.8125,0.9375 250000 600000 dense-angles
highstat $E_HIGHSTAT -1,-0.5,0,0.5,1 3000000 3000000 highstat
lowE $E_LOWE -1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1 600000 3000000 lowE
validation $VALIDATION_E -0.875,-0.625,-0.375,-0.125,0.125,0.375,0.625,0.875 250000 600000 validation
audit_a $E_AUDIT_A -0.97,-0.85,-0.45,-0.15,0.15,0.45,0.85,0.97 600000 3000000 audit-a
audit_b $E_AUDIT_B -0.93,-0.72,-0.55,-0.28,0.28,0.55,0.72,0.93 600000 3000000 audit-b
arrival $E_ARRIVAL incidence=iso,front,back,side 4000000 4000000 arrival"
    SPOT_EVENTS=4000000
    FD_OPTS="--events $(events 1000000)"
fi
SPOT_EVENTS=$(events "$SPOT_EVENTS")

# A work area is baked with one recipe.  A done step is never re-run, so
# resuming an area started with another recipe would mix banks; REDO must say
# what to redo -- the changed banks and every step after them.
RECIPE="landmark_keV $LANDMARK
precision $PRECISION
$BANK_SPECS"
if [ -f "$FW/recipe.txt" ] && [ "$(cat "$FW/recipe.txt")" != "$RECIPE" ]; then
    printf '%s\n' "$RECIPE" > "$FW/recipe.txt.new"
    diff "$FW/recipe.txt" "$FW/recipe.txt.new" || true
    [ -n "$REDO" ] || die "$FW was baked with another bank recipe ($FW/recipe.txt; this one is recipe.txt.new, diff above).  Re-run with REDO naming the changed banks and every later step, e.g. REDO=\"bank_lowE fit candidate finite_distance audit_a audit_b arrival spot\", or use a new WORK"
    stamp "bank recipe changed; running with REDO=\"$REDO\""
fi
printf '%s\n' "$RECIPE" > "$FW/recipe.txt"
rm -f "$FW/recipe.txt.new"

while read -r name energies angles min_events max_events salt <&3; do
    step "bank_$name" bank "$name" "$energies" "$angles" \
        "$min_events" "$max_events" "$salt"
done 3<<BANKS_END
$BANK_SPECS
BANKS_END

# Yield bandwidths per component: the pair-family and 511-keV fields are
# smoothed x4/3 wider (pair_no_compton and the annihilation atom x2/3
# narrower) than the shared triple.  Each yield field's node coefficients are
# then solved so its kernel average reproduces the training values
# (--yield-node-iterations; fit_detector_response.interpolate_yield_nodes).
YIELD_SCALES="pair_fraction=1.3333333333 pair_full_energy=1.3333333333
  single_escape=1.3333333333 double_escape=1.3333333333
  annihilation_511=0.6666666667 pair_no_compton=0.6666666667
  pair_single_compton=1.3333333333 pair_multiple_compton=1.3333333333"

fit() {
    local scales="" scale
    for scale in $YIELD_SCALES; do
        scales="$scales --component-yield-bandwidth-scale $scale"
    done
    # shellcheck disable=SC2086
    $PY $TOOLS/fit_detector_response.py "$BANKS/training" \
        --additional-training "$BANKS/dense_energy" \
        --additional-training "$BANKS/dense_angles" \
        --additional-training "$BANKS/highstat" \
        --additional-training "$BANKS/lowE" \
        --validation "$BANKS/validation" \
        --fit-detected-rate-surface --rate-surface-anchor \
        --rate-high-energy-ridge-kappa 100 --rate-high-energy-ridge-cut 0.75 \
        $BASIS \
        --yield-bandwidths 0.02,0.035,0.05 --shape-bandwidths 0.012,0.035,0.06 \
        $scales --yield-node-iterations 16 \
        --payload-precision 12 \
        --output "$FW/fit"
}

# The candidate is the fit's payload, scored on the validation bank on its
# own (the report the gates read).
candidate() {
    mkdir -p "$FW/candidate"
    cp "$FW/fit/detector_response_payload.txt" "$FW/candidate/"
    # shellcheck disable=SC2086
    $PY $TOOLS/fit_detector_response.py \
        --evaluate-payload "$FW/candidate/detector_response_payload.txt" \
        --validation "$BANKS/validation" \
        $BASIS --output "$FW/candidate"
}

# The point source's finite distance: virtual depths from point-source Monte
# Carlo of the detector's own geometry at 5 cm - 10 m from its reference point,
# written into the candidate payload (FiniteDistanceDepth; predict's near-field
# acceptance).
finite_distance() {
    # shellcheck disable=SC2086
    "$BIN/ceelo_detector_finite_distance" \
        --payload "$FW/candidate/detector_response_payload.txt" \
        --output "$FW/candidate/detector_response_payload.txt" \
        $FD_OPTS --threads "$JOBS"
}

audit() {
    # shellcheck disable=SC2086
    $PY $TOOLS/fit_detector_response.py \
        --evaluate-payload "$FW/candidate/detector_response_payload.txt" \
        --validation "$BANKS/$1" \
        $BASIS --output "$FW/$1"
}

# Direct-MC spot checks at normal incidence: a fresh high-statistics truth cell
# and the candidate evaluated at the same (energy, angle).
spot() {
    local energy
    mkdir -p "$FW/spot"
    for energy in 122 661.657 1332.492; do
        # shellcheck disable=SC2086
        "$BIN/ceelo_detector_response_campaign" $DETECTOR_OPTS \
            --energy-kev "$energy" --cos-theta 1 --bins 256 \
            --min-events "$SPOT_EVENTS" --max-events "$SPOT_EVENTS" \
            --seed 777 --detector-brems --detector-fluorescence --detector-spectrum-fidelity \
            > "$FW/spot/truth_$energy.csv.tmp"
        mv "$FW/spot/truth_$energy.csv.tmp" "$FW/spot/truth_$energy.csv"
        "$BIN/ceelo_detector_model_evaluate" \
            --payload "$FW/candidate/detector_response_payload.txt" \
            --energy-kev "$energy" --cos-theta 1 --bins 256 \
            --output "$FW/spot/model_$energy.csv"
        echo "spot check $energy keV written"
    done
}

# The candidate on the arrival bank: the payload averaged over each class's
# band of arrival directions against sampled-direction truth.
score_arrival() {
    $PY $TOOLS/score_arrival_bank.py "$BANKS/arrival" \
        --payload "$FW/candidate/detector_response_payload.txt" \
        --evaluator "$BIN/ceelo_detector_model_evaluate" \
        --workers "$JOBS" --output "$FW/arrival"
}

gates() {
    $PY $TOOLS/check_detector_gates.py \
        --candidate "$FW/candidate/report.json" \
        --audit audit_a="$FW/audit_a/report.json" \
        --audit audit_b="$FW/audit_b/report.json" \
        --arrival "$FW/arrival/report.json" \
        --spot "$FW/spot" \
        --output "$FW/gates.json"
}

step fit fit
step candidate candidate
step finite_distance finite_distance
step audit_a audit audit_a
step audit_b audit audit_b
step arrival score_arrival
step spot spot
gate gates gates

# A production run only gets here with every gate passed (a smoke or
# sub-final run publishes into its own $WORK/results).
publish "$FW/candidate/detector_response_payload.txt" \
    "$FW/candidate/report.json" "$FW/gates.json" "${DETECTOR}_detector_response"
stamp "DETECTOR_BAKE_DONE"
