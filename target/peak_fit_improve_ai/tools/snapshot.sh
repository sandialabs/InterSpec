#!/bin/zsh
# snapshot.sh TAG [BUILD_DIR] : freeze the current fit_peaks_corpus_eval build and the sources it came
# from as $FPR_WORK/bins/eval_TAG and $FPR_WORK/bins/src_TAG.  The build tree is shared (other sessions
# rebuild it), so measure only snapshots and keep each candidate's sources to diff or restore from.
set -u
source ${0:A:h}/env.sh
TAG=$1
BUILD=${2:-$FPR_REPO/target/testing/build_ninja}
cp $BUILD/fit_peaks_corpus_eval $FPR_WORK/bins/eval_$TAG
mkdir -p $FPR_WORK/bins/src_$TAG
for f in src/FitPeaksForNuclides.cpp InterSpec/FitPeaksForNuclides.h src/RelActCalcAuto.cpp InterSpec/RelActCalcAuto.h \
         InterSpec/RelActCalcAuto_imp.hpp src/PeakFitLM.cpp src/PeakFitUtils.cpp InterSpec/PeakFitUtils.h; do
  cp $FPR_REPO/$f $FPR_WORK/bins/src_$TAG/
done
( cd $FPR_REPO && git diff > $FPR_WORK/bins/src_$TAG/working_tree.diff )
echo "$FPR_WORK/bins/eval_$TAG"
