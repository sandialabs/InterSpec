#!/usr/bin/env bash
# Builds dist/InterSpecLight.html (and an uncompressed variant):
#   1. The WASM build's dependencies (Eigen, Abseil, Ceres)         -> deps/prefix   (first run only)
#   2. Native tools: refgen, and light_cli if NATIVE_PREFIX is set   -> build_native/
#   3. Reference-line library from InterSpec's decay data            -> build_native/gen/ref_lines.json
#   4. The WASM module (optimized for size)                          -> build_wasm/light_wasm.{js,wasm}
#   5. The self-contained pages   -> dist/InterSpecLight.html, dist/InterSpecLight_uncompressed.html
# Builds are limited to 4 cores.
#
# Environment:
#   EMSDK          required - an activated emsdk: `source <emsdk>/emsdk_env.sh`
#   NATIVE_PREFIX  optional - the dependency prefix InterSpec is built with (Eigen, Ceres): also builds
#                  light_cli, the test driver (needed by tests/run_tests.sh)
# Options:
#   --without-peak-files   leave out reading/writing peaks in N42 and SPE files (InterSpec-compatible);
#                          saves ~27 kB in the compressed page (~70 kB uncompressed).
set -euo pipefail

PEAK_FILE_IO=ON
for arg in "$@"; do
  case "${arg}" in
    --without-peak-files) PEAK_FILE_IO=OFF ;;
    *) echo "Unknown option '${arg}'" >&2; exit 1 ;;
  esac
done
REFGEN_ARGS=()
if [ "${PEAK_FILE_IO}" = ON ]; then
  REFGEN_ARGS+=( --transition-children )
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"
NCORES=4
cd "${HERE}"

if [ -z "${EMSDK:-}" ]; then
  echo "ERROR: activate emsdk first: source <emsdk>/emsdk_env.sh" >&2
  exit 1
fi
if ! command -v emcmake > /dev/null; then
  # shellcheck disable=SC1091
  source "${EMSDK}/emsdk_env.sh" > /dev/null 2>&1
fi

"${HERE}/deps/build_deps.sh"

NATIVE_ARGS=( -DCMAKE_BUILD_TYPE=Release -DLIGHT_PEAK_FILE_IO=${PEAK_FILE_IO} )
if [ -n "${NATIVE_PREFIX:-}" ]; then
  NATIVE_ARGS+=( -DCMAKE_PREFIX_PATH="${NATIVE_PREFIX}" )
fi
cmake -S . -B build_native -G Ninja "${NATIVE_ARGS[@]}" > /dev/null
cmake --build build_native -- -j${NCORES}

if [ ! -f refgen/sources.txt ]; then
  python3 refgen/collect_sources.py "${REPO}/data" refgen/sources.txt
fi
mkdir -p build_native/gen
./build_native/refgen ${REFGEN_ARGS[@]+"${REFGEN_ARGS[@]}"} "${REPO}/data/sandia.decay.xml" \
                      "${REPO}/data/sandia.reactiongamma.xml" refgen/sources.txt build_native/gen/ref_lines.json

emcmake cmake -S . -B build_wasm -G Ninja -DCMAKE_BUILD_TYPE=MinSizeRel -DLIGHT_PEAK_FILE_IO=${PEAK_FILE_IO} > /dev/null
cmake --build build_wasm --target light_wasm -- -j${NCORES}

python3 assemble_html.py --wasm-js build_wasm/light_wasm.js --ref-lines build_native/gen/ref_lines.json --output-dir dist
