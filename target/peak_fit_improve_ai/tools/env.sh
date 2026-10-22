#!/bin/zsh
# env.sh : shared settings for the peak-fit improvement scripts; sourced by each of them.  Override any
# of these in your environment before running a script.
#
#   FPR_WORK    work area: bins/ (snapshotted binaries + sources), runs/, ov/ (images), triage/, galleries/
#   FPR_INJ     GADRAS-inject corpus root (<detector>/<site>/<dwell>_seconds/*.pcf + *_truth.csv)
#   FPR_SITE    inject background site: Livermore (default), Baltimore or Denver
#   FPR_DWELL   inject dwell time in seconds: 300 (default), 30 or 1800
#   FPR_MANUAL  hand-fit corpus root holding the 17 R500 NaI hand fits.  The 168 Detective-X HPGe hand
#               fits are the harness's own default corpus (see harness/fit_peaks_corpus_eval.cpp), so
#               a run without --corpus uses them, with the inject truth attached.
#   FPR_DATA    InterSpec data directory
#   FPR_PY      a python3 with numpy + matplotlib (for the image scripts)
#   FPR_THREADS harness threads per run

FPR_TOOLS=${0:A:h}
FPR_REPO=${FPR_TOOLS:h:h:h}
: ${FPR_WORK:=$HOME/fit_peaks_work}
# The corpora default to a checkout of InterSpec_peak_fit_improve beside this repository.
: ${FPR_INJ:=${FPR_REPO:h}/InterSpec_peak_fit_improve/peak_fit_accuracy_inject_compact}
: ${FPR_SITE:=Livermore}
: ${FPR_DWELL:=300}
: ${FPR_MANUAL:=${FPR_REPO:h}/InterSpec_peak_fit_improve/manual_fits}
: ${FPR_DATA:=$FPR_REPO/data}
: ${FPR_PY:=python3}
: ${FPR_THREADS:=8}
mkdir -p $FPR_WORK/bins $FPR_WORK/runs $FPR_WORK/ov $FPR_WORK/triage $FPR_WORK/galleries

# fpr_inject_run BIN OUT DETECTOR DET_TYPE TIMEOUT [harness args...] : one inject-corpus run (site
# $FPR_SITE, dwell $FPR_DWELL) into $FPR_WORK/runs/OUT, log alongside, echoing the harness summary line.
fpr_inject_run(){
  local bin=$1 out=$2 det=$3 dt=$4 to=$5; shift 5
  rm -rf $FPR_WORK/runs/$out
  $bin --datadir=$FPR_DATA --out=$FPR_WORK/runs/$out --corpus-format=inject \
    --corpus="$FPR_INJ/$det/$FPR_SITE/${FPR_DWELL}_seconds" --det-type=$dt --no-structure --background=file \
    --threads=$FPR_THREADS --no-html --weight=min_scored_energy=20 --fit-timeout=$to "$@" \
    > $FPR_WORK/runs/$out.log 2>&1
  echo "$out: $(grep '^problems=' $FPR_WORK/runs/$out.log | tail -1 | cut -c1-160)"
}

# fpr_hand_run BIN OUT [harness args...] : the 168 Detective-X HPGe hand fits (the harness's default
# corpus, scored against the hand fits AND the inject truth) into $FPR_WORK/runs/OUT.
fpr_hand_run(){
  local bin=$1 out=$2; shift 2
  rm -rf $FPR_WORK/runs/$out
  $bin --datadir=$FPR_DATA --out=$FPR_WORK/runs/$out --background=file --threads=$FPR_THREADS \
    --fit-timeout=600 --no-html "$@" > $FPR_WORK/runs/$out.log 2>&1
  echo "$out: $(grep '^problems=' $FPR_WORK/runs/$out.log | tail -1 | cut -c1-160)"
}
