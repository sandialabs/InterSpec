#!/bin/zsh
# guard_lowres.sh BIN TAG REF_TAG [harness args...] : the scintillator guard for HPGe work - the five
# scintillator inject sets (R500, NGH and SAM-Eagle NaI, Radseeker LaBr3, Kromek GR1 CZT) in two
# parallel streams, then a per-peak bit-identity check of each against ${REF_TAG}_<set>.  Every
# scintillator fit shares the non-HPGe config, so a change that is meant for HPGe must leave all five
# untouched; that is the only way to know a shared code path did not move.  Takes about 1.5 h.
set -u
source ${0:A:h}/env.sh
BIN=$1; TAG=$2; REF=$3; shift 3
FPR_THREADS=$(( (FPR_THREADS + 1) / 2 ))
(
  fpr_inject_run $BIN ${TAG}_r500 IdentiFINDER-R500-NaI Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_sam SAM-Eagle-NaI-3x3 Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_czt Kromek-GR1-CZT CZT 1200 "$@"
) &
(
  fpr_inject_run $BIN ${TAG}_ngh IdentiFINDER-NGH Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_labr Radseeker-LaBr3 LaBr 1200 "$@"
) &
wait
for s in r500 ngh sam labr czt; do
  $FPR_PY $FPR_TOOLS/bitident.py $FPR_WORK/runs/${REF}_$s $FPR_WORK/runs/${TAG}_$s
done
