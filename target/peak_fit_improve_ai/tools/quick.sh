#!/bin/zsh
# quick.sh BIN OUT DETECTOR DET_TYPE PROBLEMS [harness args...] : fit a comma-separated list of problems
# into $FPR_WORK/runs/OUT (plot data kept), then render whole-spectrum review images to $FPR_WORK/ov/OUT.
# With FPR_CMP=<run dir> the images draw that run's fit in BLUE under the new one in RED, so each image
# is a before/after comparison (use a full run of the previous candidate, or the baseline gallery run).
#   e.g. FPR_CMP=$FPR_WORK/runs/c50_r500 quick.sh $FPR_WORK/bins/eval_c52 q52 IdentiFINDER-R500-NaI Low Pu238_Unsh,Lu177m_Unsh
set -u
source ${0:A:h}/env.sh
BIN=$1; OUT=$2; DET=$3; DT=$4; PROBS=$5; shift 5
rm -rf $FPR_WORK/runs/$OUT
$BIN --datadir=$FPR_DATA --out=$FPR_WORK/runs/$OUT --corpus-format=inject \
  --corpus="$FPR_INJ/$DET/$FPR_SITE/${FPR_DWELL}_seconds" --det-type=$DT --no-structure --background=file \
  --threads=$FPR_THREADS --no-html --weight=min_scored_energy=20 --fit-timeout=1200 --problems "$PROBS" "$@" \
  > $FPR_WORK/runs/$OUT.log 2>&1
grep '^problems=' $FPR_WORK/runs/$OUT.log | tail -1 | cut -c1-160
rm -rf $FPR_WORK/ov/$OUT
if [[ -n "${FPR_CMP:-}" ]]; then
  $FPR_PY $FPR_TOOLS/render_all.py $FPR_WORK/runs/$OUT $FPR_WORK/ov/$OUT --cmp $FPR_CMP > /dev/null 2>&1
else
  $FPR_PY $FPR_TOOLS/render_all.py $FPR_WORK/runs/$OUT $FPR_WORK/ov/$OUT > /dev/null 2>&1
fi
echo "images: $FPR_WORK/ov/$OUT"
