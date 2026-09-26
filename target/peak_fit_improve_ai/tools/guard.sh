#!/bin/zsh
# guard.sh BIN TAG REF_TAG [harness args...] : the HPGe guard alone - Detective-X inject and the HPGe
# hand-fit corpus - and a per-peak bit-identity check against runs ${REF_TAG}_detx / ${REF_TAG}_manual.
# Matching summary lines are not enough: compare every fitted peak (bitident.py).  For the opposite
# direction (HPGe work, scintillators must not move) use guard_lowres.sh.
set -u
source ${0:A:h}/env.sh
BIN=$1; TAG=$2; REF=$3; shift 3
fpr_inject_run $BIN ${TAG}_detx Detective-X High 600 "$@"
fpr_hand_run $BIN ${TAG}_manual "$@"
for s in detx manual; do
  $FPR_PY $FPR_TOOLS/bitident.py $FPR_WORK/runs/${REF}_$s $FPR_WORK/runs/${TAG}_$s
done
