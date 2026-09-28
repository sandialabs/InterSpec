#!/usr/bin/env bash
# Checks InterSpec itself reads the peaks this page writes the same as from files InterSpec wrote.
#  For each InterSpec-written N42 below, light_cli loads it and exports N42-2012 and SPE; then
#  interspec_check (built against InterSpec's library) prints the peaks InterSpec reads from:
#    - the original N42, and the page's N42 export: must be the same;
#    - an SPE InterSpec writes of the displayed samples, and the page's SPE export: must be the same.
#  `run.sh --make-fixtures` instead (re)writes tests/data/interspec_sources.{n42,spe}: files InterSpec
#  writes, with reaction, element x-ray, escape, and annihilation peaks added to a NaI spectrum.
#  Optional: needs InterSpec built (INTERSPEC_BUILD: its build directory, default <repo>/build_vscode),
#  and ./build.sh run with NATIVE_PREFIX set.  InterSpec's dependency prefix is read from its
#  CMakeCache.txt (or set INTERSPEC_PREFIX).  The light build does not depend on this.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/../.." && pwd)"
REPO="$(cd "${ROOT}/../.." && pwd)"
IS_BUILD="${INTERSPEC_BUILD:-${REPO}/build_vscode}"
PREFIX="${INTERSPEC_PREFIX:-$(sed -n 's/^CMAKE_PREFIX_PATH:[A-Z]*=//p' "${IS_BUILD}/CMakeCache.txt" 2> /dev/null | head -1)}"
IS_LIB="$(ls "${IS_BUILD}"/InterSpec.dylib "${IS_BUILD}"/libInterSpec.dylib "${IS_BUILD}"/InterSpec.so \
             "${IS_BUILD}"/libInterSpec.so 2> /dev/null | head -1 || true)"
if [ -z "${IS_LIB}" ] || [ -z "${PREFIX}" ]; then
  echo "ERROR: need InterSpec built in ${IS_BUILD} (set INTERSPEC_BUILD), and its dependency prefix" >&2
  exit 1
fi
export INTERSPEC_DATA_DIR="${REPO}/data"
T="${REPO}/target/testing/test_data"
OUT="$(mktemp -d)"
trap 'rm -rf "${OUT}"' EXIT

"${CXX:-c++}" -std=c++20 -O1 -I"${REPO}" -I"${IS_BUILD}" -I"${IS_BUILD}/external_libs/SpecUtils" \
  -I"${REPO}/external_libs" -I"${REPO}/external_libs/SpecUtils" -I"${REPO}/external_libs/SpecUtils/3rdparty" \
  -I"${REPO}/external_libs/SandiaDecay" -isystem "${PREFIX}/include" "${HERE}/interspec_check.cpp" \
  -Wl,-rpath,"${IS_BUILD}" "${IS_LIB}" -L"${PREFIX}/lib" -Wl,-rpath,"${PREFIX}/lib" -lwt \
  -o "${OUT}/interspec_check"
cd "${OUT}"

if [ "${1:-}" = "--make-fixtures" ]; then
  ./interspec_check add-peaks "${T}/PeakFitLM/Ba133_Unshielded.n42" "${ROOT}/tests/data/interspec_sources.n42"
  ./interspec_check write-spe "${ROOT}/tests/data/interspec_sources.n42" "${ROOT}/tests/data/interspec_sources.spe"
  echo "Wrote ${ROOT}/tests/data/interspec_sources.n42 and .spe"
  exit 0
fi

# The peaks InterSpec sees, without its (timestamped) debug log lines - e.g., its SPE parser logs
#  lines it does not know, which includes its own peak sections
peaks_seen(){ ./interspec_check read "$1" 2>&1 | grep -v -E "^[0-9]{8}T[0-9.]+: |^Unrecognized line|^displayed|^$" || true; }

files=( "${T}/AnalystTests/Ba133_Cs137_RandomSummingOnly.n42" "${ROOT}/tests/data/interspec_sources.n42"
        "${T}/findCharacteristics/202204_example_problem_1.n42" "${T}/manual_rel_eff/spec184_235U_12.9543.n42"
        "${T}/AnalystTests/escape_peak_check_Co56_Shielded_Fulcrum40h.n42"
        "${T}/ExportSpecFile/gapped_samples_peaks_on_9.n42" "${T}/ExportSpecFile/fore_plus_back_peaks_on_1_2.n42" )
status=0
for f in "${files[@]}"; do
  name="$(basename "${f}" .n42)"
  printf '{"method":"loadFile","params":{"path":"%s","name":"%s.n42","type":"FOREGROUND"}}\n' "${f}" "${name}" > req.jsonl
  printf '{"method":"exportFile","params":{"format":"N42-2012","path":"%s_light.n42"}}\n' "${name}" >> req.jsonl
  printf '{"method":"exportFile","params":{"format":"SPE","path":"%s_light.spe"}}\n' "${name}" >> req.jsonl
  "${ROOT}/build_native/light_cli" req.jsonl "${ROOT}/build_native/gen/ref_lines.json" > resp.jsonl
  samples="$(sed -n 2p resp.jsonl | python3 -c 'import json,sys; print(",".join(map(str, json.load(sys.stdin)["files"]["FOREGROUND"]["samples"])))')"

  if diff <(peaks_seen "${f}") <(peaks_seen "${name}_light.n42") > n42.diff; then
    n42="same"
  else
    n42="DIFFERENT"; status=1
  fi

  spe="same"
  if [ -f "${name}_light.spe" ] && grep -q PEAK_INFO_CSV "${name}_light.spe"; then
    ./interspec_check write-spe "${f}" "${name}_is.spe" "${samples}" > /dev/null 2>&1
    if ! diff <(peaks_seen "${name}_is.spe") <(peaks_seen "${name}_light.spe") > spe.diff; then
      spe="DIFFERENT"; status=1
    fi
  else
    spe="(no peaks on displayed samples ${samples})"
  fi

  echo "${name}: N42 ${n42}; SPE ${spe}"
  if [ "${n42}" = DIFFERENT ]; then head -20 n42.diff; fi
  if [ "${spe}" = DIFFERENT ]; then head -20 spe.diff; fi
done
exit ${status}
