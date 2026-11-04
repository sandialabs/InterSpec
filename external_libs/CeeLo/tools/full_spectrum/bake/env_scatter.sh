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

# Environment-scatter bake for one geometry tier (field/soil by default).
#
# Baseline (round 1): design the campaign -> transport every node's initial
# blocks into the canonical HDF5 archive -> verify its seal -> calibrate the
# uncertainty estimators -> fit every component -> record the round.  Then
# STOP and review (README.md, "Environment scatter: checkpoint").
#
# Adaptive rounds: ROUNDS=N runs up to N more rounds after the baseline, each
# continue-archive -> transport the allocated blocks -> fit -> record, and stops
# early once the campaign is stopped (two consecutive passing rounds).
#
#   tools/full_spectrum/bake/env_scatter.sh               # baseline
#   ROUNDS=3 tools/full_spectrum/bake/env_scatter.sh      # + up to 3 rounds
#   SMOKE=1 tools/full_spectrum/bake/env_scatter.sh       # smoke (24 nodes)
#
# Settings (besides common.sh's):
#   MODE=field | air          MATERIAL=soil | concrete (field only)
#   ARCHIVE_DIR=<dir>         canonical archive (~15 GB per baseline; the
#                             default is the work area -- put it on a big disk)
#   SPOOL_GIB=4               local spool high-water mark for run-streaming
#   ROUNDS=0                  adaptive rounds to run after the baseline
#   DISK_MARGIN_GIB=20        free space every round must leave on the archive
#                             and work disks, beyond its projected growth
#   FIT_CACHE=0               1 = fit from an aggregate-only copy of the archive
#                             in the work area (for an archive on a slow disk);
#                             0 fits the archive directly (same model)
#
# The transport workers use every core of the machine; JOBS does not limit
# them.  Results do not depend on the thread count (seeds are per replica).
. "$(dirname "$0")/common.sh"

MODE=${MODE:-field}
MATERIAL=${MATERIAL:-soil}
SPOOL_GIB=${SPOOL_GIB:-4}
ROUNDS=${ROUNDS:-0}
DISK_MARGIN_GIB=${DISK_MARGIN_GIB:-20}
FIT_CACHE=${FIT_CACHE:-0}
if [ "$MODE" = air ]; then TAG=air; else TAG=${MODE}_${MATERIAL}; fi
family env_$TAG
ARCHIVE_DIR=${ARCHIVE_DIR:-$FW}
mkdir -p "$ARCHIVE_DIR"

need_bin env_kernel_bake env_fit_block_tool env_fit_payload_eval

CEEFIT="$PY -m tools.full_spectrum.env_scatter_fit"
MANIFEST=$FW/campaign.json
ARCHIVE=$ARCHIVE_DIR/${TAG}_archive.h5
CACHE=$FW/fit_cache.h5
SPOOL=$FW/spool

# Physical domain: E0 20-3000 keV; distance 0.25-25 m; source height 0-4 m;
# detector height 0.3-4 m.  Below 25 cm the environment adds <= 1.2% of the
# direct fluence; detector heights below ~0.3 m cost 4-116x the histories
# (next-event variance near the ground) for a <= 4% change.
DOMAIN="--energy-min 20 --energy-max 3000 --distance-min 0.25 --distance-max 25"
if [ "$MODE" != air ]; then
    DOMAIN="$DOMAIN --source-height-min 0 --source-height-max 4
        --detector-height-min 0.3 --detector-height-max 4"
fi
if [ "$SMOKE" = 1 ]; then
    # The smallest valid design: the 18 forced training anchors (GEANT4
    # anchors, thresholds, domain faces), 4 validation nodes -- two above the
    # pair threshold, so every component has validation support -- and 2 final
    # nodes.  24 nodes x 4 blocks x 262,144 histories.
    SIZES="--train 18 --validation 4 --final 2 --candidates 256"
    # Only three training nodes lie above the pair threshold, and an NMF rank
    # cannot exceed a component's supporting node count.
    FIT_OPTS="--ranks 2 --levels 1,2"
else
    SIZES="--train 256 --validation 64 --final 64 --candidates 10000"
    FIT_OPTS=""
fi

design() {
    # A campaign manifest binds every block ever transported for it: never
    # regenerate it mid-campaign, even with REDO.
    if [ -f "$MANIFEST" ]; then
        echo "reusing the existing campaign manifest $MANIFEST"
        return 0
    fi
    local material=""
    [ "$MODE" = air ] || material="--material $MATERIAL"
    # shellcheck disable=SC2086
    $CEEFIT design --mode "$MODE" $material $SIZES $DOMAIN \
        --executable "$BIN/env_kernel_bake" \
        --block-output-root "$SPOOL/blocks" \
        --output "$MANIFEST"
}

# streaming MANIFEST [--incremental-seal]: transport every block of MANIFEST
# that is not yet in the archive, pack it, verify it and seal the archive.  The
# adaptive rounds seal incrementally (only the nodes a round changed are
# re-verified); the baseline and the final seal before publishing are full.
# Resumes after interruption.
streaming() {
    $CEEFIT run-streaming --manifest "$1" --archive "$ARCHIVE" \
        --spool "$SPOOL" --spool-max-gib "$SPOOL_GIB" --delete-packed-blocks ${2:-}
    $CEEFIT verify-archive-seal --archive "$ARCHIVE"
}

# full_seal: make sure the archive carries a full seal (publication requires
# one), re-verifying every node's raw-block closure only when the current seal
# is incremental or stale.
full_seal() {
    if $CEEFIT verify-archive-seal --archive "$ARCHIVE" 2> /dev/null \
            | grep -q '"verification": "full-'; then
        echo "the archive already carries a full seal"
    else
        $CEEFIT seal-archive --archive "$ARCHIVE"
    fi
}

calibrate() {
    $CEEFIT calibrate-uncertainty --archive "$ARCHIVE" \
        --output "$FW/calibration.json"
}

# fit ROUND CAMPAIGN: every rank/level candidate, fitted from the sealed
# archive (or, with FIT_CACHE=1, from a fresh aggregate-only cache bound to the
# archive snapshot).  fit-archive exits nonzero when validation fails --
# normal before the adaptive rounds converge -- and then writes a .candidate
# payload instead; either way the round is recorded next.
fit() {
    local dir=$FW/round$1 source=$ARCHIVE
    mkdir -p "$dir"
    rm -f "$CACHE" "$CACHE.seal.json"
    if [ "$FIT_CACHE" = 1 ]; then
        $CEEFIT fit-cache --archive "$ARCHIVE" --output "$CACHE"
        source=$CACHE
    fi
    set +e
    # shellcheck disable=SC2086
    $CEEFIT fit-archive --archive "$source" --campaign "$2" $FIT_OPTS \
        --output "$dir/$TAG.ceefitp" --report "$dir/acceptance.json" \
        --acquisition-utilities "$dir/utilities.json"
    local status=$?
    set -e
    echo "fit-archive exit status $status"
    [ -f "$dir/acceptance.json" ] \
        || die "fit-archive wrote no acceptance report; read the error above"
    [ -f "$dir/$TAG.ceefitp" ] || [ -f "$dir/$TAG.ceefitp.candidate" ] \
        || die "fit-archive wrote neither a payload nor a candidate"
}

# record ROUND CAMPAIGN: persist the round (the validation report is checked
# independently against the archive and the candidate).
record() {
    local dir=$FW/round$1 candidate
    candidate=$dir/$TAG.ceefitp
    [ -f "$candidate" ] || candidate=$candidate.candidate
    $CEEFIT campaign-round --manifest "$2" \
        --validation-report "$dir/acceptance.json" \
        --archive "$ARCHIVE" --candidate "$candidate" \
        --output "$dir/manifest.json"
}

# allocate ROUND: choose the next round's work from the previous round.
allocate() {
    local previous=$FW/round$(( $1 - 1 )) dir=$FW/round$1
    mkdir -p "$dir"
    $CEEFIT continue-archive --manifest "$previous/manifest.json" \
        --archive "$ARCHIVE" --adaptive-utilities "$previous/utilities.json" \
        --output "$dir/allocated.json" --commands-output "$dir/commands.sh"
}

# disk_guard ROUND: before a round transports anything, project what it adds
# -- a ceiling from its scheduled blocks: raw records for blocks 1..4 of new
# nodes, continued nodes' fold states at their uncompressed size, new
# aggregate (and fit-cache) rows -- and stop, before writing, unless the
# archive's disk and the work area's keep DISK_MARGIN_GIB (plus the spool)
# free afterwards.  Free space and re-run to resume.
disk_guard() {
    "$PY" - "$FW/round$1/allocated.json" "$ARCHIVE" "$FIT_CACHE" "$ARCHIVE_DIR" \
        "$FW" "$DISK_MARGIN_GIB" "$SPOOL_GIB" <<'EOF'
import json, os, shutil, sys
from tools.full_spectrum.env_scatter_fit.storage import projected_archive_growth
allocated, archive, fit_cache, archive_dir, work, margin, spool = sys.argv[1:]
gib = 1024.0 ** 3
pending = json.load(open(allocated))["campaign_state"]["continuation_schedule"]["pending_work"]
projection = projected_archive_growth(
    archive,
    [(int(item["node_id_hex"], 16), int(item["block_id"])) for item in pending],
    fit_cache=fit_cache == "1",
)
growth = projection["total_bytes"]
# Summed per disk: the archive and the work area may share one.
need = {}
for path, bytes_needed in ((archive_dir, growth), (work, float(spool) * gib)):
    device = os.stat(path).st_dev
    need[device] = (path, need.get(device, (path, 0.0))[1] + bytes_needed)
print(f"round adds {len(pending)} blocks: archive about "
      f"+{projection['expected_bytes'] / gib:.1f} GiB, at most +{growth / gib:.1f} GiB "
      f"(raw records +{projection['raw_record_bytes'] / gib:.1f}, new rows "
      f"+{projection['aggregate_bytes'] / gib:.1f})")
if not projection["reuses_freed_space"]:
    print("  the archive predates persistent free-space tracking: replaced fold "
          "states are not reused (a new archive avoids it)")
short = False
for path, bytes_needed in need.values():
    free = shutil.disk_usage(path).free
    left = (free - bytes_needed) / gib
    print(f"  {path}: {free / gib:.1f} GiB free, {left:.1f} GiB left after the round")
    if left < float(margin):
        short = True
if short:
    print(f"disk guard: the round would leave less than {margin} GiB free; free space "
          "(or move ARCHIVE_DIR to a bigger disk) and re-run the same command",
          file=sys.stderr)
    sys.exit(1)
EOF
}

campaign_state() {
    "$PY" - "$1" <<'EOF'
import json, sys
state = json.load(open(sys.argv[1]))["campaign_state"]
print(f"round {state['round']}: last action {state['last_action']}, "
      f"{len(state['training_node_ids'])} training nodes, "
      f"stopped={state['stopped']}")
sys.exit(0 if state["stopped"] else 3)
EOF
}

# Baseline.
step design design
step streaming streaming "$MANIFEST"
step calibrate calibrate
step fit1 fit 1 "$MANIFEST"
step record1 record 1 "$MANIFEST"

# Adaptive rounds.
round=1
stopped=0
campaign_state "$FW/round1/manifest.json" || true
while [ "$round" -le "$ROUNDS" ]; do
    if campaign_state "$FW/round$round/manifest.json" > /dev/null; then
        stopped=1
        break
    fi
    next=$(( round + 1 ))
    step "allocate$next" allocate "$next"
    if [ ! -f "$DONE/streaming$next" ]; then
        disk_guard "$next" || die "disk guard: round $next not started"
    fi
    step "streaming$next" streaming "$FW/round$next/allocated.json" --incremental-seal
    step "fit$next" fit "$next" "$FW/round$next/allocated.json"
    step "record$next" record "$next" "$FW/round$next/allocated.json"
    round=$next
    campaign_state "$FW/round$round/manifest.json" || true
done
if campaign_state "$FW/round$round/manifest.json" > /dev/null; then
    stopped=1
fi

last=$FW/round$round
payload=$last/$TAG.ceefitp
[ -f "$payload" ] || payload=$payload.candidate
if [ "$SMOKE" = 1 ] || [ "$stopped" = 1 ]; then
    full_seal
    publish "$payload" "$last/acceptance.json" "$last/manifest.json" \
        "$FW/calibration.json" "env_scatter_$TAG"
fi
if [ "$stopped" = 1 ]; then
    stamp "ENV_SCATTER_BAKE_DONE ($TAG): the campaign stopped after round $round."
else
    stamp "ENV_SCATTER_ROUND_DONE ($TAG): round $round recorded; the campaign continues."
    echo "  Review before more rounds: $FW/calibration.json,"
    echo "  $last/acceptance.json and $last/utilities.json (README.md, checkpoint)."
    echo "  Then: ROUNDS=<n> $0"
fi
