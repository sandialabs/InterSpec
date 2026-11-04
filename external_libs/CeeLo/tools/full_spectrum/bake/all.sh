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

# Every family, unattended, in one command, for a machine left to bake:
# preflight -> build -> timing probe -> each detector -> source escape ->
# charged particle -> environment scatter (baseline, then ENV_ROUNDS adaptive
# rounds) -> the status report in $WORK/bake_report.txt.
#
#   nohup caffeinate -is nice -n 20 tools/full_spectrum/bake/all.sh \
#       > bake_all.out 2>&1 &
#   tools/full_spectrum/bake/status.sh       # any time, from another shell
#
# A family that fails does not stop the others (they are independent); the
# report says which stopped where, and each driver's output is in
# $WORK/<stage>.out.  A failed preflight or build stops everything.  Re-run
# the same command to resume: every driver skips the steps it has done.
#
# Settings (besides common.sh's, which every driver inherits):
#   DETECTORS="hpge_coax detective_x"   detector payloads (detector.sh DETECTOR)
#   FAMILIES="detector source_escape charged_particle env_scatter"
#                     ("detector" alone below PRECISION=final: only it has the
#                     precision tiers)
#   ENV_ROUNDS=0      adaptive environment rounds after the baseline; the
#                     campaign stops by itself after two passing rounds, so a
#                     full unattended campaign is a large number, e.g. 40
#   ARCHIVE_DIR       the environment archive's directory (env_scatter.sh)
. "$(dirname "$0")/common.sh"
export SMOKE PRECISION BUILD PY JOBS WORK RESULTS
DETECTORS=${DETECTORS:-"hpge_coax detective_x"}
if [ "$PRECISION" = final ]; then
    FAMILIES=${FAMILIES:-"detector source_escape charged_particle env_scatter"}
else
    FAMILIES=${FAMILIES:-detector}
fi
ENV_ROUNDS=${ENV_ROUNDS:-0}
export DETECTORS FAMILIES
BAKE=tools/full_spectrum/bake
STAGES=$WORK/all/stages.txt
mkdir -p "$WORK/all"

# stage NAME COMMAND...: run one driver with its output in $WORK/NAME.out and
# its outcome appended to $STAGES.  Returns the driver's exit status.
stage() {
    local name=$1 start status
    shift
    start=$(date '+%Y-%m-%d %H:%M:%S')
    stamp "stage $name: $*"
    set +e
    "$@" > "$WORK/$name.out" 2>&1
    status=$?
    set -e
    printf '%-24s %s -> %s  exit %d  %s\n' "$name" "$start" \
        "$(date '+%Y-%m-%d %H:%M:%S')" "$status" \
        "$(grep -E '_DONE|ERROR' "$WORK/$name.out" | tail -1 | cut -c1-90)" >> "$STAGES"
    stamp "stage $name: exit $status"
    return "$status"
}

stamp "all.sh: WORK=$WORK SMOKE=$SMOKE PRECISION=$PRECISION DETECTORS=\"$DETECTORS\" FAMILIES=\"$FAMILIES\" ENV_ROUNDS=$ENV_ROUNDS"
echo "--- $(date '+%Y-%m-%d %H:%M:%S') all.sh started" >> "$STAGES"

stage preflight "$BAKE/preflight.sh" || die "preflight failed: $WORK/preflight.out"
stage build "$BAKE/build.sh" || die "build failed: $WORK/build.out"
if [ ! -f "$WORK/timing_probe/done/probe" ]; then
    stage timing_probe "$BAKE/timing_probe.sh" || true
fi

failed=""
for family in $FAMILIES; do
    case $family in
        detector)
            for detector in $DETECTORS; do
                stage "detector_$detector" env DETECTOR="$detector" "$BAKE/detector.sh" \
                    || failed="$failed detector_$detector"
            done ;;
        source_escape|charged_particle)
            stage "$family" "$BAKE/$family.sh" || failed="$failed $family" ;;
        env_scatter)
            stage env_scatter env ROUNDS="$ENV_ROUNDS" "$BAKE/env_scatter.sh" \
                || failed="$failed env_scatter" ;;
        *) die "unknown family $family" ;;
    esac
done

"$BAKE/status.sh" > "$WORK/bake_report.txt" 2>&1 || true
stamp "report: $WORK/bake_report.txt"
if [ -n "$failed" ]; then
    stamp "ALL_DONE with failures:$failed (each stage's output is $WORK/<stage>.out)"
    exit 1
fi
stamp "ALL_DONE: every stage finished"
