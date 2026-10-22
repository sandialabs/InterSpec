#!/bin/zsh
# full.sh BIN TAG [harness args...] : the full evaluation of one candidate, in two parallel streams of
# FPR_THREADS/2 threads each:
#   ${TAG}_r500, ${TAG}_ngh, ${TAG}_sam   NaI inject corpora (det-type Low, 1200 s timeout)
#   ${TAG}_labr, ${TAG}_czt               LaBr3 and CZT inject corpora (they get the non-HPGe config too)
#   ${TAG}_detx                           HPGe inject corpus (Detective-X, det-type High)
#   ${TAG}_manual                         the HPGe hand-fit corpus (the harness's default corpus)
# Extra harness args (e.g. --set field=value to run a flag variant) go to every run.
# With FPR_GUARD_REF=<tag> the runs named in FPR_GUARD_SETS - default "detx manual", the HPGe guard;
# "r500 ngh sam labr czt" when it is the scintillators that must not move - are then checked peak by
# peak for bit-identity against <tag>'s.  Takes 1.5-2 h on 10 cores.
set -u
source ${0:A:h}/env.sh
BIN=$1; TAG=$2; shift 2
FPR_THREADS=$(( (FPR_THREADS + 1) / 2 ))
(
  fpr_inject_run $BIN ${TAG}_r500 IdentiFINDER-R500-NaI Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_sam SAM-Eagle-NaI-3x3 Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_czt Kromek-GR1-CZT CZT 1200 "$@"
  fpr_hand_run $BIN ${TAG}_manual "$@"
) &
(
  fpr_inject_run $BIN ${TAG}_ngh IdentiFINDER-NGH Low 1200 "$@"
  fpr_inject_run $BIN ${TAG}_labr Radseeker-LaBr3 LaBr 1200 "$@"
  fpr_inject_run $BIN ${TAG}_detx Detective-X High 600 "$@"
) &
wait
if [[ -n "${FPR_GUARD_REF:-}" ]]; then
  for s in ${=${FPR_GUARD_SETS:-detx manual}}; do
    $FPR_PY $FPR_TOOLS/bitident.py $FPR_WORK/runs/${FPR_GUARD_REF}_$s $FPR_WORK/runs/${TAG}_$s
  done
fi
echo FULL_DONE
