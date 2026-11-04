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

# Shared settings and helpers for the full-spectrum bake drivers.  Sourced by
# every driver in this directory; never run directly.  Written for bash 3.2
# (macOS /bin/bash) and newer.
#
# Settings are environment variables; every one has a default:
#   SMOKE=0         1 = smoke run: the same steps and code paths on tiny banks,
#                   minutes instead of hours.  Its numbers mean nothing, its
#                   gates are reported but not enforced, and it writes only
#                   under $WORK (never into data/full_spectrum/).
#   PRECISION=final fast | medium | final: the event budget, for a driver
#                   that offers the tiers (detector.sh): fast = 1/8, medium =
#                   1/2 of every final bank and spot check.  A fast or medium
#                   run is a real bake at lower precision -- to see where a
#                   detector stands in an hour -- so, like SMOKE, it writes
#                   only under its own $WORK (bake_work/precision_<tier>), and
#                   reports its gates without enforcing them.  A driver without
#                   the tiers refuses to run below final.
#   BUILD=build-rel Release build tree, relative to the CeeLo root.
#   PY=python3      Python >= 3.10 with numpy, scipy, h5py (and matplotlib).
#   JOBS=<cores-1>  concurrent single-threaded MC cells / worker processes.
#   WORK=bake_work  banks, fits, reports and logs (bake_work/smoke if SMOKE=1).
#   RESULTS=data/full_spectrum
#                   where accepted payloads are copied ($WORK/results if SMOKE=1).
#   REDO=""         space-separated step names to run again although they
#                   completed ("all" = every step of the driver).
#
# Every step is recorded in $WORK/<family>/done/<step> when it succeeds, so an
# interrupted driver is resumed by running the same command again.

set -euo pipefail

CEELO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
cd "$CEELO_ROOT"

SMOKE=${SMOKE:-0}
PRECISION=${PRECISION:-final}
BUILD=${BUILD:-build-rel}
PY=${PY:-python3}
if [ -z "${JOBS:-}" ]; then
    _cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)
    JOBS=$(( _cores > 1 ? _cores - 1 : 1 ))
fi
case $PRECISION in
    final|medium|fast) ;;
    *) printf 'ERROR: PRECISION must be fast, medium or final (not %s)\n' \
           "$PRECISION" >&2; exit 1 ;;
esac
if [ "$SMOKE" = 1 ] && [ "$PRECISION" != final ]; then
    printf 'ERROR: SMOKE=1 and PRECISION=%s are exclusive\n' "$PRECISION" >&2
    exit 1
fi
if [ "$SMOKE" = 1 ]; then
    WORK=${WORK:-bake_work/smoke}
    RESULTS=${RESULTS:-$WORK/results}
elif [ "$PRECISION" != final ]; then
    WORK=${WORK:-bake_work/precision_$PRECISION}
    RESULTS=${RESULTS:-$WORK/results}
else
    WORK=${WORK:-bake_work}
    RESULTS=${RESULTS:-data/full_spectrum}
fi
REDO=${REDO:-}
case $BUILD in
    /*) ;;
    *) BUILD=$CEELO_ROOT/$BUILD ;;
esac
TOOLS=tools/full_spectrum
BIN=$BUILD/tools/full_spectrum
# Progress lines from the Python runners should appear as they happen, not in
# 4 KiB bursts when their output is piped into a log.
export PYTHONUNBUFFERED=1

stamp() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

# precision_tiers: a driver that scales its event budgets with events()
# declares so before calling family; any other driver refuses PRECISION below
# final.
PRECISION_TIERS=0
precision_tiers() { PRECISION_TIERS=1; }

# events N: N events scaled to the PRECISION tier.
events() {
    case $PRECISION in
        fast) echo $(( $1 / 8 )) ;;
        medium) echo $(( $1 / 2 )) ;;
        *) echo "$1" ;;
    esac
}

# family NAME: select the per-family work area; call once, near the top.
family() {
    if [ "$PRECISION" != final ] && [ "$PRECISION_TIERS" != 1 ]; then
        die "$1 has no precision tiers; run it with PRECISION=final"
    fi
    FAMILY=$1
    FW=$WORK/$1
    LOGS=$FW/logs
    DONE=$FW/done
    mkdir -p "$FW" "$LOGS" "$DONE"
    stamp "$FAMILY bake: WORK=$WORK SMOKE=$SMOKE PRECISION=$PRECISION JOBS=$JOBS BUILD=$BUILD PY=$PY"
}

need_bin() {
    local b
    for b in "$@"; do
        [ -x "$BIN/$b" ] || die "missing $BIN/$b -- run tools/full_spectrum/bake/build.sh first"
    done
}

# GEANT4 11.4's geant4.sh exports GEANT4_DATA_DIR (not the per-dataset
# G4LEDATA etc. of older releases); accept either, and check the datasets are
# really there, since a missing mount only fails deep inside a GEANT4 run.
geant4_data_dir() {
    if [ -n "${GEANT4_DATA_DIR:-}" ]; then
        echo "$GEANT4_DATA_DIR"
    elif [ -n "${G4LEDATA:-}" ]; then
        dirname "$G4LEDATA"
    fi
}

need_geant4() {
    local data
    data=$(geant4_data_dir)
    [ -n "$data" ] || die "the GEANT4 environment is not set: 'source <geant4 install>/bin/geant4.sh' and re-run"
    ls -d "$data"/G4EMLOW* > /dev/null 2>&1 \
        || die "GEANT4 datasets not found under $data (install or mount the GEANT4 11.4 datasets there)"
    need_bin g4_source_escape_campaign g4_photonuclear_xs
}

_redo() {
    local item
    for item in $REDO; do
        if [ "$item" = all ] || [ "$item" = "$1" ]; then return 0; fi
    done
    return 1
}

# step NAME COMMAND [ARGS...]: run COMMAND (usually a shell function of the
# driver) once.  Output goes to the terminal and to $LOGS/NAME.log.  On success
# $DONE/NAME is written and a later run skips the step; on failure the driver
# stops with the log's path.
step() {
    local name=$1 status
    shift
    if [ -f "$DONE/$name" ] && ! _redo "$name"; then
        stamp "skip $name (already done; REDO=$name re-runs it)"
        return 0
    fi
    rm -f "$DONE/$name"
    stamp "run  $name"
    # errexit is re-enabled inside the subshell: a step body is a function, and
    # bash ignores `set -e` for functions run where their status is tested.
    set +e
    ( set -e; "$@" ) 2>&1 | tee "$LOGS/$name.log"
    status=${PIPESTATUS[0]}
    set -e
    if [ "$status" -ne 0 ]; then
        die "step $name failed (exit $status); its log is $LOGS/$name.log"
    fi
    date '+%Y-%m-%d %H:%M:%S' > "$DONE/$name"
    stamp "done $name"
}

# gate NAME COMMAND [ARGS...]: like step, for an acceptance check.  A failing
# gate stops a production bake; in a smoke run it is reported and the driver
# carries on, because smoke banks are far too small to pass anything, and so
# it is below PRECISION=final, whose noise floors are not the gates'.
gate() {
    local name=$1 status
    shift
    stamp "gate $name"
    set +e
    ( set -e; "$@" ) 2>&1 | tee "$LOGS/$name.log"
    status=${PIPESTATUS[0]}
    set -e
    if [ "$status" -eq 0 ]; then
        stamp "PASS $name"
        return 0
    fi
    if [ "$SMOKE" = 1 ]; then
        stamp "gate $name did not pass (expected in a smoke run: exit $status)"
        return 0
    fi
    if [ "$PRECISION" != final ]; then
        stamp "gate $name did not pass at PRECISION=$PRECISION (exit $status; enforced at final only)"
        return 0
    fi
    die "gate $name FAILED (exit $status); its log is $LOGS/$name.log"
}

# publish SRC... DEST_DIR: replace the contents of $RESULTS/DEST_DIR with the
# accepted outputs, so no file of an earlier bake is left beside them.
publish() {
    local dest=${!#}
    local count=$(( $# - 1 ))
    mkdir -p "$RESULTS/$dest"
    rm -f "$RESULTS/$dest"/*
    cp "${@:1:$count}" "$RESULTS/$dest/"
    stamp "published to $RESULTS/$dest/"
    ls -l "$RESULTS/$dest/"
}
