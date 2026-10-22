#!/bin/zsh
# search_all.sh BIN TAG [harness args...] : the automated peak search alone (fit_peaks_corpus_eval
# --search-only), scored against the inject truth, on every guard detector at 30, 300 and 1800 s:
#   HPGe  Detective-X, Detective-EX, HPGe_Planar_50%, Falcon 5000, Fulcrum40h
#   NaI   IdentiFINDER-R500-NaI, IdentiFINDER-NGH, SAM-Eagle-NaI-3x3
#   LaBr3 Radseeker-LaBr3;  CZT Kromek-GR1-CZT
# into $FPR_WORK/runs/${TAG}_<set>_<dwell>.  Extra args (e.g. --search-set detection_z_chi2_weights=1, or a
# --search-sweep) go to every run.  A few minutes in all; compare two tags with search_cmp.py.
set -u
source ${0:A:h}/env.sh
BIN=$1; TAG=$2; shift 2
typeset -A SETS
SETS=( detx "Detective-X|High" detex "Detective-EX|High" planar "HPGe_Planar_50%|High" falcon "Falcon 5000|High"
       fulcrum "Fulcrum40h|High" r500 "IdentiFINDER-R500-NaI|Low" ngh "IdentiFINDER-NGH|Low"
       sam "SAM-Eagle-NaI-3x3|Low" labr "Radseeker-LaBr3|LaBr" czt "Kromek-GR1-CZT|CZT" )
for dwell in 30 300 1800; do
  for s in detx detex planar falcon fulcrum r500 ngh sam labr czt; do
    det=${SETS[$s]%%|*}; dt=${SETS[$s]##*|}
    dir="$FPR_INJ/$det/$FPR_SITE/${dwell}_seconds"
    [[ -d $dir ]] || continue
    out=${TAG}_${s}_${dwell}
    rm -rf $FPR_WORK/runs/$out
    $BIN --datadir=$FPR_DATA --out=$FPR_WORK/runs/$out --corpus-format=inject --corpus="$dir" --det-type=$dt \
      --no-structure --background=file --threads=$FPR_THREADS --weight=min_scored_energy=20 --search-only "$@" \
      > $FPR_WORK/runs/$out.log 2>&1
    echo "$out: $(grep -E '^problems=|: problems=' $FPR_WORK/runs/$out.log | tail -1 | cut -c1-230)"
  done
done
echo SEARCH_ALL_DONE
