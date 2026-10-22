#!/usr/bin/env bash
# Runs each tests/*.jsonl request script through the native light_cli and the WASM module (in
#  Node), and checks both succeed with the same results, and meet any expectations the script has
#  (see check_expectations.py).  Scripts starting "# requires: peakFileIo" are skipped for builds
#  without that option.  "@REPO@" in scripts is replaced by the path of the InterSpec repository.
#  Run after ./build.sh with NATIVE_PREFIX set (which builds light_cli), and emsdk activated.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/.." && pwd)"
REPO="$(cd "${ROOT}/../.." && pwd)"
NODE="${EMSDK_NODE:-node}"
CLI="${ROOT}/build_native/light_cli"
REF="${ROOT}/build_native/gen/ref_lines.json"
if [ ! -x "${CLI}" ]; then
  echo "ERROR: ${CLI} not built; run ./build.sh with NATIVE_PREFIX set (see README.md)" >&2
  exit 1
fi
OUT="$(mktemp -d)"
trap 'rm -rf "${OUT}"' EXIT
cd "${OUT}"  #scripts write exports to relative paths

PEAK_FILE_IO="$(echo '{"method":"getState"}' | "${CLI}" \
                 | python3 -c 'import json,sys; print(json.load(sys.stdin)["features"]["peakFileIo"])')"

status=0
for source_script in "${HERE}"/*.jsonl; do
  name="$(basename "${source_script}" .jsonl)"
  if grep -q '^# requires: peakFileIo' "${source_script}" && [ "${PEAK_FILE_IO}" != True ]; then
    echo "SKIP ${name} (built without peak file support)"; continue
  fi
  script="${OUT}/${name}.jsonl"
  sed "s#@REPO@#${REPO}#g" "${source_script}" > "${script}"
  if ! "${CLI}" "${script}" "${REF}" > "${OUT}/${name}.native"; then
    echo "FAIL ${name}: native run had errors"; status=1; continue
  fi
  if ! "${NODE}" "${HERE}/node_run.mjs" "${ROOT}/build_wasm/light_wasm.js" "${script}" "${REF}" > "${OUT}/${name}.wasm"; then
    echo "FAIL ${name}: wasm run had errors"; status=1; continue
  fi
  if ! python3 "${HERE}/compare_outputs.py" "${OUT}/${name}.native" "${OUT}/${name}.wasm"; then
    echo "FAIL ${name}: native and wasm results differ"; status=1; continue
  fi
  if ! python3 "${HERE}/check_expectations.py" "${script}" "${OUT}/${name}.native" \
     || ! python3 "${HERE}/check_expectations.py" "${script}" "${OUT}/${name}.wasm"; then
    echo "FAIL ${name}: expectations not met"; status=1; continue
  fi
  echo "PASS ${name}"
done
exit ${status}
