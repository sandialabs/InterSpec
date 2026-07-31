#!/bin/bash

handle_error() {
    echo "An error occurred building InterSpec, on line $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# set -e # Cause script to exit as soon as there is an error

# takes three or four arguments: [path to InterSpec code] [Path to CMake build Directory]
#                                [Working directory path] [optional: dependency prefix]
#  - InterSpec code is mapped to `/interspec`
#  - build directory mapped to `/build_app`
#  - build working directory mapped to /build_working_dir
#  - dependency prefix, if given, mapped to /deps_prefix
# For example: ./build_app.sh /interspec /build_app /build_working_dir /deps_prefix
#
# With four arguments the dependencies (boost, zlib, Wt, Eigen, Ceres) are taken from a prefix
# built by target/electron/build_linux_deps_from_docker.sh, which is what CI does - it caches that
# prefix, so a warm run compiles none of them.  With three arguments CMake's FetchContent downloads
# and builds them into the build directory instead: ~25 minutes longer, but needs no prior setup.
#
# This script *compiles* InterSpec and leaves the app payload in <CMake build dir>/app.  It does
# not package anything: turning that directory into a .deb/.rpm/AppImage/tarball is done by
# electron-builder, run outside this container (see .github/workflows/build_app.yml).  The split
# is deliberate - this container exists only because we need to compile against an old glibc, and
# the Debian/RPM tooling electron-builder wants is far easier to get on a normal Ubuntu host.
#
# This script is intended to be run from within the `quay.io/pypa/manylinux_2_28_x86_64` Docker image (Electron 28+ requires glibc 2.28+), for example, to build InterSpec you could
#  cd /tmp
#  mkdir build_interspec
#  cd build_interspec
#  git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_code
#  mkdir build_electron
#  mkdir build_working_dir
#  docker run --rm -it -v `pwd`/InterSpec_code:/interspec -v `pwd`/build_electron/:/build_app -v `pwd`/build_working_dir:/build_working_dir quay.io/pypa/manylinux_2_28_x86_64:latest /interspec/target/electron/build_linux_app_from_docker.sh /interspec /build_app /build_working_dir
#  #Then the app payload will be in build_electron/app, ready for:
#  #  cd InterSpec_code/target/electron && npm install \
#  #    && npx electron-builder --linux -c.directories.app=/tmp/build_interspec/build_electron/app

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
  echo "Error: This script requires three or four arguments."
  echo "    [path to InterSpec code] [Path to CMake build Directory] [Working directory path] [optional: dependency prefix]"
  exit 1  # Exit with an error code
fi

InterSpecCodePath=$1
CmakeBuildDir=$2
WorkingDir=$3
DepsPrefix=${4:-}

# Released Linux builds are x86_64, but the aarch64 manylinux image is a useful way to shake this
#  script out at native speed on an arm64 machine (and an arm64 Linux build is a long-standing
#  wish - see README.md).  Follow the machine we are actually on rather than hard-coding x64.
case "$(uname -m)" in
  x86_64)  NODE_ARCH="x64" ;;
  aarch64) NODE_ARCH="arm64" ;;
  *)       echo "Error: unsupported architecture $(uname -m)"; exit 1 ;;
esac
echo "Building for ${NODE_ARCH} ($(uname -m))"

echo "cd'ing into /build_working_dir"
cd ${WorkingDir}

cp "${InterSpecCodePath}/target/electron/package.json" .

git config --global --add safe.directory ${InterSpecCodePath}
export GIT_HASH=$(git -C ${InterSpecCodePath} rev-parse HEAD)
echo "GIT_HASH = ${GIT_HASH}"

echo "Will install npm and global packages"
yum update -y
yum install -y npm
# Use the 'n' package to install a fairly modern version of npm, by default we are on like version 6 or 8, which is too old to run some of the packages we need
npm install -g n
n 22.11.0
# Just in case, make sure to pick up latest node version
hash -r
# The default version of npm in manylinux_2_28_x86_64 is 6.14.11, which doesnt support the `exec` command, so we'll locally install a newer version of npm
npm install -g npm@10.9.0
hash -r
# uglify-js/uglifycss are used by cmake/DeployJsAndCss.cmake to minify our JS and CSS
npm install -g uglify-js
npm install -g uglifycss
npm install -g cmake-js

echo "Will install local npm packages"
# node-addon-api supplies napi.h; the Electron runtime itself is fetched later by
#  electron-builder, on the host, so it is not needed here.
npm install --save-dev node-addon-api --arch=${NODE_ARCH}

echo "CWD"
pwd
echo "ls node_modules"
ls node_modules

# We need to help find where node_modules directory is located in the Docker container
export NODE_MODULES_PARENT_DIR=$PWD


# -------------------------------------------------------------------------------------------------
# Where the dependencies come from.
#
# CMAKE_PREFIX_PATH is enough on its own - it outranks CMAKE_SYSTEM_PREFIX_PATH in every find_*()
#  search order, so the prefix beats /usr.  ZLIB_ROOT is the one belt-and-braces worth adding: this
#  image has /usr/include/zlib.h and /usr/lib64/libz.so, and picking those up would put
#  `libz.so.1` into the addon's DT_NEEDED, which check_elf_compat.py rejects on purpose.
#
# We deliberately do NOT pass -DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=OFF the way the Windows build
#  does.  On Linux that takes /usr/lib64 out of find_library()'s reach, and the top-level
#  CMakeLists.txt has `find_library( RT_LIB rt REQUIRES )` - the keyword is misspelt, so it parses
#  as a search path rather than REQUIRED and a miss silently yields RT_LIB-NOTFOUND, which then
#  gets handed straight to target_link_libraries().  find_program() for uglifyjs/terser and
#  FindThreads' -lpthread probe want system paths too.  (-DNO_CMAKE_SYSTEM_PATH=ON, also on the
#  Windows line, is a no-op: NO_CMAKE_SYSTEM_PATH is a find_*() keyword, not a variable.)
# -------------------------------------------------------------------------------------------------
DEP_ARGS=()
if [ -n "${DepsPrefix}" ]; then
  if [ ! -f "${DepsPrefix}/.interspec_deps_complete" ]; then
    echo "Error: '${DepsPrefix}' is not a completed InterSpec dependency prefix (no"
    echo "  .interspec_deps_complete sentinel).  Build one with"
    echo "  target/electron/build_linux_deps_from_docker.sh, or omit the fourth argument to have"
    echo "  CMake FetchContent build the dependencies instead."
    exit 1
  fi
  echo "Using dependency prefix ${DepsPrefix}:"
  cat "${DepsPrefix}/.interspec_deps_complete"
  DEP_ARGS=( --CDInterSpec_FETCH_DEPENDENCIES=OFF
             --CDCMAKE_PREFIX_PATH="${DepsPrefix}"
             --CDZLIB_ROOT="${DepsPrefix}" )
else
  echo "No dependency prefix given - CMake will fetch and build the dependencies."
  DEP_ARGS=( --CDInterSpec_FETCH_DEPENDENCIES=ON )
fi


# -------------------------------------------------------------------------------------------------
# ccache, for InterSpec's own ~330 objects.  (With a dependency prefix that is all that is left to
#  compile.)  Unlike a cached build directory, ccache is content-addressed, so it is immune to
#  `actions/checkout` rewriting every source mtime - which is the reason caching the build tree
#  never saved CI any compilation at all.
#
# ccache is not in the AlmaLinux 8 repos this image uses (el8 moved it to EPEL, which manylinux
#  does not enable), so take the official release binary: the musl-static build, which is fully
#  static and so has no glibc floor of its own, as a .tar.gz so we need no `xz` either.
#
# Never fail the build over this - ccache is an optimization, not a correctness requirement.
#
# NOTE ON ERROR HANDLING: `set -e` is commented out at the top of this file and we rely on a
#  `trap ... ERR`, which does not fire for commands inside an `if` condition or an && chain.  That
#  is why the install is written as one `if ... && ... ; then` chain: a failed download must not
#  abort the whole build.
# -------------------------------------------------------------------------------------------------
CCACHE_VERSION="4.13.6"
CCACHE_ARGS=()
if ! command -v ccache > /dev/null 2>&1; then
  case "$(uname -m)" in
    x86_64)  _cc_arch="x86_64";  _cc_sha256="09e0547a0c3b250a76675c33130366f1399f3580842fb360c052520d56214ead" ;;
    aarch64) _cc_arch="aarch64"; _cc_sha256="bff0e0c19165db8627c85c36b0885b3b180659eda67a028298d26565aad52f56" ;;
    *)       _cc_arch=""; _cc_sha256="" ;;
  esac

  if [ -z "${_cc_arch}" ]; then
    echo "WARNING: no pinned ccache build for $(uname -m) - compiling without ccache"
  else
    _cc_dir="ccache-${CCACHE_VERSION}-linux-${_cc_arch}-musl-static"
    _cc_tar="/tmp/${_cc_dir}.tar.gz"
    _cc_url="https://github.com/ccache/ccache/releases/download/v${CCACHE_VERSION}/${_cc_dir}.tar.gz"
    echo "Installing ccache ${CCACHE_VERSION} (${_cc_arch})"
    if curl -fsSL "${_cc_url}" -o "${_cc_tar}" \
       && echo "${_cc_sha256}  ${_cc_tar}" | sha256sum --check --status \
       && tar -xzf "${_cc_tar}" -C /tmp \
       && install -m 0755 "/tmp/${_cc_dir}/ccache" /usr/local/bin/ccache ; then
      echo "Installed ccache from the official release tarball"
    elif yum install -y epel-release > /dev/null 2>&1 && yum install -y ccache > /dev/null 2>&1 ; then
      echo "Installed ccache from EPEL"
    else
      echo "WARNING: could not install ccache - compiling without it"
    fi
  fi
fi

if command -v ccache > /dev/null 2>&1; then
  ccache --version | head -n 1
  # Zero the counters so the stats printed at the end of this script describe THIS build.
  ccache --zero-stats > /dev/null 2>&1 || true
  CCACHE_ARGS=( --CDCMAKE_C_COMPILER_LAUNCHER=ccache --CDCMAKE_CXX_COMPILER_LAUNCHER=ccache )
  # Note: we do not set CCACHE_SLOPPINESS=time_macros.  src/AppUtils.cpp and src/BatchInfoLog.cpp
  #  use __DATE__ (the latter stamps the compile date into batch-output N42 files), and ccache
  #  folds the current date into the hash for such translation units - so they hit normally within
  #  a day and re-compile after the date rolls over, which is exactly what we want.  (Only
  #  __TIME__ makes a file outright uncacheable, and nothing here uses it.)  Measured: a rebuild
  #  in a fresh build directory against a warm cache was 333/333 hits.
fi


echo "Will build InterSpec code"
# Note: libstdc++/libgcc are linked *dynamically* - see cmake/CxxRuntimePolicy.cmake.  Passing
#  `-static-libstdc++` here is what caused issue #51 (segfault in dl_init before main()); the
#  gcc-toolset compiler in this image already caps the required GLIBCXX at the AlmaLinux 8
#  base, so a dynamic link keeps the glibc-2.28 floor.
#  `--CDCMAKE_SHARED_LINKER_FLAGS=""` is passed explicitly rather than just omitted: `--CD` sets
#  a *cache* entry, so a build directory configured back when the flag was passed keeps it
#  forever.  CI no longer caches its build directory, so the hazard is now confined to the
#  reusable-build-directory workflow documented at the top of this file - but that workflow is
#  real, so the flags stay.  (cmake/CxxRuntimePolicy.cmake also strips these flags defensively;
#  this is the cheaper first line of defence.)
CMAKE_BUILD_PARALLEL_LEVEL=`nproc` cmake-js --directory ${InterSpecCodePath}/target/electron \
  --architecture "${NODE_ARCH}" --arch="${NODE_ARCH}" \
  --CDCMAKE_BUILD_TYPE="Release" \
  --CDCMAKE_SHARED_LINKER_FLAGS="" --CDCMAKE_EXE_LINKER_FLAGS="" --CDCMAKE_MODULE_LINKER_FLAGS="" \
  --CDInterSpec_NODE_MODULE_DIR="${WorkingDir}/node_modules" \
  "${DEP_ARGS[@]}" ${CCACHE_ARGS[@]+"${CCACHE_ARGS[@]}"} \
  --CDBUILD_AS_LOCAL_SERVER=OFF --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="${LEAFLET_KEY}" \
  --CDUSE_REL_ACT_TOOL=ON \
  --out="${CmakeBuildDir}" --target install

echo "ls ${CmakeBuildDir} (CmakeBuildDir)"
ls ${CmakeBuildDir}

if [ ! -f "${CmakeBuildDir}/app/InterSpecAddOn.node" ]; then
  echo "Error: ${CmakeBuildDir}/app/InterSpecAddOn.node was not produced."
  exit 1
fi

echo "This is a build of InterSpec for Linux, using InterSpec code git hash ${GIT_HASH}." > "${CmakeBuildDir}/app/build_information.txt"
echo "This build is untested - please contact InterSpec@sandia.gov for support." >> "${CmakeBuildDir}/app/build_information.txt"

# Fail here rather than shipping a binary that will not load on a supported distro; this is the
#  guard for https://github.com/sandialabs/InterSpec/issues/51.
echo "Checking the built binaries against the runtime floor we advertise"
python3 "${InterSpecCodePath}/target/electron/linux/check_elf_compat.py" "${CmakeBuildDir}/app"


# What actually got compiled.  This is the check that the dependency prefix is doing its job: with
#  a prefix there must be no boost_*/wt/wthttp/wtdbo*/zlib/ceres rows here at all, because those
#  targets do not exist in this build tree.  (Caching the *build directory* used to leave all 875
#  of them being recompiled every run, which is how we found out that cache was worthless.)
echo "=== object files produced, by target ==="
find "${CmakeBuildDir}" -name '*.o' | sed 's|.*/CMakeFiles/||; s|\.dir/.*||' | sort | uniq -c | sort -rn
echo "total objects: $(find "${CmakeBuildDir}" -name '*.o' | wc -l)"

if command -v ccache > /dev/null 2>&1; then
  echo "=== ccache statistics for this build ==="
  ccache --show-stats --verbose 2>/dev/null || ccache --show-stats || true
fi

echo "Done - app payload is in ${CmakeBuildDir}/app"
