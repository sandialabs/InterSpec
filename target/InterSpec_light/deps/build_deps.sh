#!/usr/bin/env bash
# Downloads the WASM build's dependencies at the versions below, and cross-compiles them with
#  Emscripten into deps/prefix: Eigen (headers), Abseil, and Ceres.
#
# Needs an activated emsdk (`source <emsdk>/emsdk_env.sh`), git, CMake and Ninja.
# Finished steps are skipped; delete deps/build and deps/prefix to rebuild from scratch.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${HERE}/src"
BUILD="${HERE}/build"
PREFIX="${HERE}/prefix"
NCORES=4

# The versions the page is built and tested with
EIGEN_URL="https://gitlab.com/libeigen/eigen.git"
EIGEN_COMMIT="bc3b39870ecb690a623a3f49149a358b95c5781d"    # 5.0.1 (2025-11-08)
ABSEIL_URL="https://github.com/abseil/abseil-cpp.git"
ABSEIL_COMMIT="255c84dadd029fd8ad25c5efb5933e47beaa00c7"   # LTS 2026-01, patch 1
CERES_URL="https://github.com/ceres-solver/ceres-solver.git"
CERES_COMMIT="2f946a582ae4a9e7ee0492030ec12d9b1f3dbade"    # master, 2026-03-22 (2.3 pre-release)

if ! command -v emcmake > /dev/null; then
  if [ -z "${EMSDK:-}" ]; then
    echo "ERROR: activate emsdk first: source <emsdk>/emsdk_env.sh" >&2
    exit 1
  fi
  # shellcheck disable=SC1091
  source "${EMSDK}/emsdk_env.sh" > /dev/null 2>&1
fi

# Must match the flags the app is compiled with (see ../CMakeLists.txt).  The prefix map keeps
#  this machine's paths out of the WASM (e.g., in Ceres' log messages).
WASM_FLAGS="-fwasm-exceptions -sWASM_LEGACY_EXCEPTIONS=1 -ffile-prefix-map=${HERE}/=deps/"

mkdir -p "${SRC}" "${BUILD}" "${PREFIX}/include"

# Shallow-fetches one commit (the commit hash also checks the download)
fetch_commit(){
  local name="$1" url="$2" commit="$3"
  local dir="${SRC}/${name}"
  if [ -d "${dir}" ]; then
    if [ "$(git -C "${dir}" rev-parse HEAD)" != "${commit}" ]; then
      echo "ERROR: ${dir} is not at ${commit}; delete it to re-download" >&2
      exit 1
    fi
    return
  fi
  git init -q "${dir}"
  git -C "${dir}" fetch -q --depth 1 "${url}" "${commit}"
  git -C "${dir}" -c advice.detachedHead=false checkout -q FETCH_HEAD
}

# Builds and installs one CMake project, unless already done at this commit
build_dep(){
  local name="$1" commit="$2"
  shift 2
  local stamp="${BUILD}/${name}.stamp"
  if [ -f "${stamp}" ] && [ "$(cat "${stamp}")" = "${commit}" ]; then
    return
  fi
  emcmake cmake -S "${SRC}/${name}" -B "${BUILD}/${name}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_CXX_STANDARD=20 \
    -DCMAKE_C_FLAGS="${WASM_FLAGS}" \
    -DCMAKE_CXX_FLAGS="${WASM_FLAGS}" \
    -DCMAKE_FIND_ROOT_PATH="${PREFIX}" \
    -DBUILD_SHARED_LIBS=OFF \
    "$@"
  cmake --build "${BUILD}/${name}" --target install -j ${NCORES}
  echo "${commit}" > "${stamp}"
}

fetch_commit eigen "${EIGEN_URL}" "${EIGEN_COMMIT}"
build_dep eigen "${EIGEN_COMMIT}" \
  -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF -DBUILD_TESTING=OFF \
  -DEIGEN_BUILD_PKGCONFIG=OFF -DEIGEN_BUILD_BLAS=OFF -DEIGEN_BUILD_LAPACK=OFF

fetch_commit abseil "${ABSEIL_URL}" "${ABSEIL_COMMIT}"
build_dep abseil "${ABSEIL_COMMIT}" \
  -DABSL_BUILD_TESTING=OFF -DABSL_USE_GOOGLETEST_HEAD=OFF -DABSL_PROPAGATE_CXX_STD=ON

EIGEN_CMAKE_DIR="$(dirname "$(find "${PREFIX}" -name Eigen3Config.cmake | head -1)")"

fetch_commit ceres "${CERES_URL}" "${CERES_COMMIT}"
build_dep ceres "${CERES_COMMIT}" \
  -Dabsl_DIR="${PREFIX}/lib/cmake/absl" -DEigen3_DIR="${EIGEN_CMAKE_DIR}" \
  -DLAPACK=OFF -DSUITESPARSE=OFF -DEIGENSPARSE=OFF -DACCELERATESPARSE=OFF -DUSE_CUDA=OFF \
  -DSCHUR_SPECIALIZATIONS=OFF -DCUSTOM_BLAS=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
  -DBUILD_BENCHMARKS=OFF -DBUILD_DOCUMENTATION=OFF -DEXPORT_BUILD_DIR=OFF \
  -DPROVIDE_UNINSTALL_TARGET=OFF

# A -pthread build produces shared-memory WASM, which cannot run from file://
if grep -rl -- "-pthread" "${BUILD}"/*/build.ninja "${PREFIX}/lib/cmake" 2> /dev/null | grep -q .; then
  echo "ERROR: -pthread found in dependency build flags" >&2
  exit 1
fi

echo "Dependencies installed to ${PREFIX}"
