#!/bin/zsh
# run_corpus.sh BIN OUT DETECTOR DET_TYPE [TIMEOUT] [harness args...] : one full inject-corpus run.
#   DET_TYPE is Low (NaI/CsI), LaBr, CZT or High - the harness defaults to High, so scintillators MUST
#   say so.  TIMEOUT defaults to 1200 s: U/Pu problems take 2-5 min and a 600 s cap times out under load.
#   e.g. run_corpus.sh $FPR_WORK/bins/eval_c52 c52_r500 IdentiFINDER-R500-NaI Low
set -u
source ${0:A:h}/env.sh
BIN=$1; OUT=$2; DET=$3; DT=$4; TO=${5:-1200}
(( $# >= 5 )) && shift 5 || shift $#
fpr_inject_run $BIN $OUT $DET $DT $TO "$@"
