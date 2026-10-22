#!/bin/zsh
# full_hpge.sh BIN TAG [harness args...] : the HPGe evaluation of one candidate, in two parallel streams
# of FPR_THREADS/2 threads each (det-type High, 600 s timeout, Livermore background):
#   ${TAG}_detx, ${TAG}_detx30, ${TAG}_detx1800   Detective-X inject at 300 / 30 / 1800 s
#   ${TAG}_detex                                  Detective-EX inject, 300 s
#   ${TAG}_planar                                 HPGe_Planar_50% inject, 300 s
#   ${TAG}_falcon                                 Falcon 5000 inject, 300 s
#   ${TAG}_fulcrum                                Fulcrum40h inject, 300 s (most truth files are empty)
#   ${TAG}_manual                                 the 168 Detective-X hand fits
# Extra harness args go to every run.  The scintillator guard is guard_lowres.sh.
set -u
source ${0:A:h}/env.sh
BIN=$1; TAG=$2; shift 2
FPR_THREADS=$(( (FPR_THREADS + 1) / 2 ))
(
  fpr_inject_run $BIN ${TAG}_detx Detective-X High 600 "$@"
  FPR_DWELL=30 fpr_inject_run $BIN ${TAG}_detx30 Detective-X High 600 "$@"
  FPR_DWELL=1800 fpr_inject_run $BIN ${TAG}_detx1800 Detective-X High 600 "$@"
  fpr_inject_run $BIN ${TAG}_fulcrum Fulcrum40h High 600 "$@"
) &
(
  fpr_hand_run $BIN ${TAG}_manual "$@"
  fpr_inject_run $BIN ${TAG}_detex Detective-EX High 600 "$@"
  fpr_inject_run $BIN ${TAG}_planar HPGe_Planar_50% High 600 "$@"
  fpr_inject_run $BIN ${TAG}_falcon "Falcon 5000" High 600 "$@"
) &
wait
echo HPGE_DONE
