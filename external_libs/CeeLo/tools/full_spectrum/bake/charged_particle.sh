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

# Charged-particle source bake: training bank (three sub-grids) and two
# off-grid electron banks -> structured fit with the frozen structure -> joint
# shape fit -> scoring -> C++/Python parity -> gates -> publish.
#
#   tools/full_spectrum/bake/charged_particle.sh            # production
#   SMOKE=1 tools/full_spectrum/bake/charged_particle.sh    # smoke run (minutes)
#
# Settings: see common.sh.  Resumable: re-run the same command after an
# interruption; a bank step skips the cells it already wrote.  The banks are
# 11,802 cells, 2.99 G events; every later step is minutes to tens of minutes.
# A failed gate therefore means "fix the fitter or its settings, then
# REDO=\"fit_structured joint_extract joint_build joint_score parity\"" (or
# from joint_build on, for a shape-only change), never "re-bake the banks".
# REDO of a bank step only fills in missing cells; to regenerate a bank,
# delete its directory under $WORK/charged_particle/banks.
#
# Frozen hyperparameters (selected once, not re-selected by the bake):
#   bake/config/charged_particle_structure.json       per-component yield
#       grids and smoothing, and the four-stage shape rank / weight grid
#   bake/config/charged_particle_joint_settings.json  joint shape fit: rank,
#       weight grid, field smoothness, residual norm, cell weighting
#       (fit_charged_particle_joint_shape.py reselect re-derives it)
. "$(dirname "$0")/common.sh"
family charged_particle

need_bin ceelo_source_electron_campaign ceelo_source_cross_section_features \
    ceelo_charged_particle_model_evaluate

HELPER=$BIN/ceelo_source_cross_section_features
BANKS=$FW/banks
CONFIG=$TOOLS/bake/config
PAYLOAD=charged_particle_source_payload.txt
# Whole-Z holdout planes: the fits never see Z 35 or 64 ("test" split).
HOLDOUT_Z=35,64

# grid BANK SALT EVENTS ENERGIES ELEMENTS RATIOS [--electrons-only | --positrons-only]
# One tensor-product sub-grid (T keV x Z x thickness in CSDA ranges x
# particle) of a bank; 512 bins, 10 keV photon floor.  Every cell's seed
# derives from SALT and the cell's identity, so a bank is reproducible cell by
# cell and banks with different salts are independent.
grid() {
    "$PY" "$TOOLS/run_charged_particle_source_campaign.py" \
        --executable "$BIN/ceelo_source_electron_campaign" \
        --cross-section-helper "$HELPER" \
        --output "$BANKS/$1" --seed-salt "$2" --events "$3" \
        --energies-kev "$4" --elements "$5" --range-ratios "$6" \
        --bins 512 --minimum-photon-kev 10 --jobs "$JOBS" --skip-existing \
        "${@:7}"
}

# expect_cells BANK COUNT: the bank must hold exactly COUNT cells.
expect_cells() {
    local count
    count=$(find "$BANKS/$1/cells" \( -name '*_ceelo.csv' -o -name '*_ceelo.csv.gz' \) \
        | sed 's/\.gz$//' | sort -u | wc -l | tr -d ' ')
    if [ "$count" -ne "$2" ]; then
        die "bank $1 holds $count cells, expected $2; more than expected means $BANKS/$1 holds cells of another grid (delete it and re-run)"
    fi
    echo "bank $1: $count cells, as expected"
}

if [ "$SMOKE" = 1 ]; then
    # The same sub-grids and steps on a handful of cells at 20k events: the
    # holdout planes, one sparse Z, both particles in every stage.
    FULL_Z=6,35,64,82
    A_T=100,1000,5000;  A_R=0.2,2;  A_Z=$FULL_Z;  A_EVENTS=20000
    SPARSE_Z=50;  SPARSE_EVENTS=20000
    B_T=300,2000;  B_R=0.5,2;  B_EVENTS=20000
    C_T=600;  C_R=0.4,1.5;  C_EVENTS=20000
    TRAINING_CELLS=92
    HOLDOUT_FACTOR=2
    V_T=225,2150;  V_Z=32,70;  V_R=0.47,1.9;  V_EVENTS=20000;  V_CELLS=8
    F_T=325,3800;  F_Z=16,79;  F_R=0.75,2.3;  F_EVENTS=20000
    F_POSITRON_EVENTS=20000;  F_CELLS=16
    ROUNDS=2
    PROBES=2
else
    FULL_Z=1,4,6,8,13,20,26,29,35,42,56,64,74,82,92
    # Sub-grid A, both particles: 13 T x 11 r over the full Z list plus Z 50
    # at 60k events, and Z 32, 70 at 120k.
    A_T=50,100,200,300,500,750,1000,1500,2290,3500,5000,8000,12000
    A_R=0.02,0.05,0.1,0.2,0.4,0.7,1,1.5,2,5,10
    A_Z=$FULL_Z,50;  A_EVENTS=60000
    SPARSE_Z=32,70;  SPARSE_EVENTS=120000
    # Sub-grid B, electrons only, 200k events: dense in T and through the
    # containment transition (r 0.2-6).
    B_T=120,150,175,250,400,450,550,700,800,1200,1300,1800,2000,2500,2800,4200,6500,10000
    B_R=0.2,0.3,0.4,0.55,0.7,0.85,1,1.25,1.5,1.75,2,2.5,3,3.5,5,6
    B_EVENTS=200000
    # Sub-grid C, both particles, 300k events: 240-900 keV, where the
    # annihilation components change fastest.
    C_T=240,270,330,375,600,675,825,900
    C_R=0.1,0.2,0.4,0.7,1,1.5,2
    C_EVENTS=300000
    # 13*11*2*16 + 13*11*2*2 + 18*16*15 + 8*7*2*15 = 4576 + 572 + 4320 + 1680
    TRAINING_CELLS=11148
    # The holdout planes' positron cells (sub-grids A and C at Z 35, 64) carry
    # 16x their sub-grid's events: at 1x the in-flight annihilation yields
    # sit at their own MC noise floor, so the holdout medians measure noise.
    HOLDOUT_FACTOR=16
    # Off-grid banks: T, Z and r all off the training grid.  The validation
    # bank is electrons; the final test adds positrons at 4x the electron
    # events, the only off-grid positron truth.
    V_T=130,225,275,475,625,875,1100,1400,2150,3200,5500,9000
    V_Z=32,70;  V_R=0.27,0.47,0.62,0.93,1.38,1.9,2.8,4.2
    V_EVENTS=300000;  V_CELLS=192
    F_T=160,325,525,775,975,1250,1650,2400,3800,7200,11000
    F_Z=16,47,79;  F_R=0.35,0.75,1.1,1.6,2.3,3.3,5.2
    F_EVENTS=400000;  F_POSITRON_EVENTS=1600000;  F_CELLS=462
    ROUNDS=16
    PROBES=6
fi

training_bank() {
    # The holdout planes' positron cells first, so the grids below skip them.
    grid training charged-training $(( HOLDOUT_FACTOR * A_EVENTS )) \
        "$A_T" "$HOLDOUT_Z" "$A_R" --positrons-only
    grid training charged-training $(( HOLDOUT_FACTOR * C_EVENTS )) \
        "$C_T" "$HOLDOUT_Z" "$C_R" --positrons-only
    grid training charged-training "$A_EVENTS" "$A_T" "$A_Z" "$A_R"
    grid training charged-training "$SPARSE_EVENTS" "$A_T" "$SPARSE_Z" "$A_R"
    grid training charged-training "$B_EVENTS" "$B_T" "$FULL_Z" "$B_R" \
        --electrons-only
    grid training charged-training "$C_EVENTS" "$C_T" "$FULL_Z" "$C_R"
    expect_cells training "$TRAINING_CELLS"
}

offgrid_validation_bank() {
    grid offgrid_validation charged-offgrid-validation "$V_EVENTS" \
        "$V_T" "$V_Z" "$V_R" --electrons-only
    expect_cells offgrid_validation "$V_CELLS"
}

# Scored once, after everything else is frozen; never used for selection.
offgrid_final_test_bank() {
    grid offgrid_final_test charged-offgrid-final-test "$F_EVENTS" \
        "$F_T" "$F_Z" "$F_R" --electrons-only
    grid offgrid_final_test charged-offgrid-final-test "$F_POSITRON_EVENTS" \
        "$F_T" "$F_Z" "$F_R" --positrons-only
    expect_cells offgrid_final_test "$F_CELLS"
}

# Yield fields (kept byte for byte in the final payload) and the four-stage
# shape model the joint fit starts from.  The bank's hash-validation split is
# folded into the final fit and scored with a train-only model.
fit_structured() {
    "$PY" "$TOOLS/fit_charged_particle_source_structured.py" \
        "$BANKS/training" --output "$FW/structured" \
        --offgrid-validation "$BANKS/offgrid_validation" \
        --offgrid-final-test "$BANKS/offgrid_final_test" \
        --structure-json "$CONFIG/charged_particle_structure.json" \
        --cross-section-helper "$HELPER" \
        --holdout-z "$HOLDOUT_Z" --photon-absorption-path-fraction 2
}

joint() {
    "$PY" "$TOOLS/fit_charged_particle_joint_shape.py" "$@"
}

joint_extract() {
    joint extract --fit-dir "$FW/structured" --output "$FW/joint/data"
}

# The final payload: joint shape fit on the training bank's train and
# validation splits, starting from the structured payload.
joint_build() {
    joint build --fit-dir "$FW/structured" --data "$FW/joint/data" \
        --settings "$CONFIG/charged_particle_joint_settings.json" \
        --rounds "$ROUNDS" --output "$FW/joint"
}

joint_score() {
    joint score --fit-dir "$FW/structured" --payload "$FW/joint/$PAYLOAD" \
        --banks test,offgrid_validation,offgrid_final_test \
        --output "$FW/joint/validation_report.json"
}

parity() {
    joint parity --fit-dir "$FW/structured" --payload "$FW/joint/$PAYLOAD" \
        --evaluator "$BIN/ceelo_charged_particle_model_evaluate" \
        --probes-per-particle "$PROBES" --output "$FW/joint/parity.json"
}

gates() {
    "$PY" "$TOOLS/check_charged_particle_gates.py" \
        --payload "$FW/joint/$PAYLOAD" \
        --report "$FW/joint/validation_report.json" \
        --parity "$FW/joint/parity.json" --output "$FW/gates.json"
}

step bank_training training_bank
step bank_offgrid_validation offgrid_validation_bank
step bank_offgrid_final_test offgrid_final_test_bank
step fit_structured fit_structured
step joint_extract joint_extract
step joint_build joint_build
step joint_score joint_score
step parity parity
gate gates gates

# A production run only gets here with every gate passed.
publish "$FW/joint/$PAYLOAD" "$FW/joint/validation_report.json" \
    "$FW/joint/parity.json" "$FW/joint/build_notes.json" "$FW/gates.json" \
    charged_particle_source
stamp "CHARGED_PARTICLE_BAKE_DONE"
