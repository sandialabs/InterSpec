#!/bin/zsh
# run_hand.sh BIN OUTNAME [harness args...] : fit the 17 hand-fit R500 NaI spectra (reference = the hand
# fit, truth from the inject corpus), merged into $FPR_WORK/runs/OUTNAME
set -u
source ${0:A:h}/env.sh
BIN=$1; OUT=$FPR_WORK/runs/$2; shift 2
rm -rf $OUT; mkdir -p $OUT/plot_data $OUT/parts
one(){ local p=$1 s=$2; shift 2
  $BIN --datadir=$FPR_DATA --out=$OUT/parts/$p \
    --corpus=$FPR_MANUAL/IdentiFINDER-R500-NaI_300_seconds \
    --truth-inject=$FPR_INJ/IdentiFINDER-R500-NaI/Livermore/300_seconds --problems $p --sources $s \
    --det-type=Low --no-structure --background=file --threads=1 --weight=min_scored_energy=20 --fit-timeout=1200 --no-html "$@" \
    > $OUT/parts/$p.log 2>&1
  cp $OUT/parts/$p/plot_data/*.json $OUT/plot_data/ 2>/dev/null
}
LIST=( Ac225_Unsh:Ac225 Am241_Sh:Am241,Pb Bi213_Phantom:Bi213 Br76_Phantom:Br76 Eu152_Unsh:Eu152 Ir192_Sh:Ir192,Pb
  Ir192_Unsh:Ir192 Lu177m_Sh:Lu177m,Pb Lu177m_Unsh:Lu177m Pu239_Sh:Pu,Pu238,Pu239,Pu240,Pu241,Pb Ra223_Phantom:Ra223
  Ra226_Unsh:Ra226 Th228_Unsh:Th228 Th232_Unsh:Th232 U233_Sh:U,U232,U233,U234,U235,U238,Pb
  Uore_Unsh:U,U232,U234,U235,U238,Ra226 Yb169_Phantom:Yb169 )
for e in $LIST; do
  one ${e%%:*} ${e#*:} "$@" &
  while [ $(jobs -r | wc -l) -ge 6 ]; do sleep 1; done
done
wait
echo HAND_DONE $OUT
