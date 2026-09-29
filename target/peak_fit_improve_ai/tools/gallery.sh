#!/bin/zsh
# gallery.sh BIN NAME [DETECTOR:DET_TYPE ...] : HTML galleries (the harness's gallery.html, with the
# background drawn) for one binary, into $FPR_WORK/galleries/NAME/<detector>/.  Defaults to the three NaI
# inject detectors plus LaBr3 and CZT; e.g. `gallery.sh BIN c5 Detective-X:High manual` for HPGe, where
# `manual` is the 168 Detective-X hand fits (reference fits shown too).  Inject sets use $FPR_SITE and
# $FPR_DWELL.  Open with: open $FPR_WORK/galleries/NAME/*/gallery.html
set -u
source ${0:A:h}/env.sh
BIN=$1; NAME=$2; shift 2
DETS=( "$@" )
(( ${#DETS} )) || DETS=( IdentiFINDER-R500-NaI:Low IdentiFINDER-NGH:Low SAM-Eagle-NaI-3x3:Low Radseeker-LaBr3:LaBr Kromek-GR1-CZT:CZT )
OUT=$FPR_WORK/galleries/$NAME
mkdir -p $OUT
for pair in $DETS; do
  det=${pair%%:*}; dt=${pair##*:}
  rm -rf $OUT/$det
  if [[ $det == manual ]]; then
    $BIN --datadir=$FPR_DATA --out=$OUT/manual --background=file --threads=$FPR_THREADS --gallery-background \
      --gallery-max 400 --fit-timeout=600 > $OUT/manual.log 2>&1
    echo "manual: $(grep '^problems=' $OUT/manual.log | tail -1 | cut -c1-140)"
    continue
  fi
  $BIN --datadir=$FPR_DATA --out=$OUT/$det --corpus-format=inject --corpus="$FPR_INJ/$det/$FPR_SITE/${FPR_DWELL}_seconds" \
    --det-type=$dt --no-structure --background=file --threads=$FPR_THREADS --gallery-background --gallery-max 400 \
    --weight=min_scored_energy=20 --fit-timeout=1200 > $OUT/$det.log 2>&1
  echo "$det: $(grep '^problems=' $OUT/$det.log | tail -1 | cut -c1-140)"
done
echo "open $OUT/*/gallery.html"
