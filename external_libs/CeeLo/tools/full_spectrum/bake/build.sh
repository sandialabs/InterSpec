#!/usr/bin/env bash
# CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
# and X-ray detector efficiency - developed as part of InterSpec.
#
# Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Stage 0 of every bake: configure and build the Release tree, then run the C++
# unit tests, the Python tests and the C++/Python lockstep checks.  Every step
# runs every time (an up-to-date build is quick); stop if anything fails.
#
#   source <geant4 install>/bin/geant4.sh      # so CMake finds GEANT4
#   tools/full_spectrum/bake/build.sh
#
# Dependencies CMake cannot find on its own go in CMAKE_ARGS, e.g.
#   CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/opt/deps" tools/full_spectrum/bake/build.sh
. "$(dirname "$0")/common.sh"
precision_tiers   # nothing to scale
family build
REDO=all

prerequisites() {
    command -v cmake > /dev/null || die "cmake is not on PATH"
    "$PY" - <<'EOF'
import sys
if sys.version_info < (3, 10):
    sys.exit(f"Python >= 3.10 is required; {sys.executable} is {sys.version.split()[0]}")
missing = []
for module in ("numpy", "scipy", "h5py", "matplotlib"):
    try:
        __import__(module)
    except ImportError:
        missing.append(module)
if missing:
    sys.exit("missing Python modules: " + " ".join(missing)
             + "  (python3 -m pip install -r tools/full_spectrum/requirements.txt)")
print(f"Python {sys.version.split()[0]} at {sys.executable}: numpy, scipy, h5py, matplotlib present")
EOF
    if [ -z "$(geant4_data_dir)" ]; then
        echo "Note: the GEANT4 environment is not set, so the GEANT4 programs will not"
        echo "      be built. Only the source-escape bake's optional PHOTONUCLEAR=1 and"
        echo "      GEANT4_CHECK=1 parts need them."
    fi
}

configure() {
    # Release, portable code generation: -march=native changes floating-point
    # results between machines and is unsupported under emulation.
    # shellcheck disable=SC2086
    cmake -S . -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DCEELO_NATIVE_ARCH=OFF \
        -DCEELO_BUILD_TESTS=ON ${CMAKE_ARGS:-} 2>&1 | tee "$LOGS/cmake.out"
    # Without Boost.Test CMake quietly builds no unit tests, and ctest would
    # then "pass" having tested nothing.
    if grep -q 'Boost.Test not found' "$LOGS/cmake.out"; then
        die "Boost.Test not found: install Boost (brew install boost) or point CMAKE_ARGS=-DCMAKE_PREFIX_PATH=... at it"
    fi
}

compile() {
    cmake --build "$BUILD" -j "$JOBS"
    if [ -x "$BIN/g4_source_escape_campaign" ]; then
        echo "GEANT4 programs built: $BIN/g4_source_escape_campaign, $BIN/g4_photonuclear_xs"
    else
        echo "Note: GEANT4 programs not built (GEANT4 not found at configure time)."
        echo "      Only the source-escape bake's optional PHOTONUCLEAR=1 and"
        echo "      GEANT4_CHECK=1 parts need them; to build them, source geant4.sh,"
        echo "      delete $BUILD/CMakeCache.txt and run this script again."
    fi
}

unit_tests() {
    (cd "$BUILD" && ctest -j "$JOBS" --output-on-failure)
}

python_tests() {
    "$PY" -m unittest discover -s $TOOLS/env_scatter_fit/tests -t .
    "$PY" -m unittest discover -s $TOOLS/tests -t .
}

# C++/Python parity of the shared bases and features (the B-spline check is
# pure Python against the C++ normalization constants).
lockstep() {
    "$PY" $TOOLS/check_basis_lockstep.py
    "$PY" $TOOLS/check_feature_lockstep.py --build "$BUILD"
    "$PY" $TOOLS/check_base_lockstep.py --build "$BUILD"
}

step prerequisites prerequisites
step configure configure
step compile compile
step unit_tests unit_tests
step python_tests python_tests
step lockstep lockstep
stamp "BUILD_DONE: $BUILD is built and verified"
