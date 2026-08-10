#!/bin/sh

# Exit immediately if a command exits with a non-zero status
set -e

# Check the number of arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <interspec_src_dir> <working_directory> [<install_directory>]"
    exit 1
fi

# Set proxy, if you need it
#export http_proxy=http://proxy.sandia.gov:80
#export https_proxy=http://proxy.sandia.gov:80

# Everything this script needs is either a POSIX basic (sh, curl, tar, gzip, awk, sha256sum) or a
#  build tool (a C/C++ compiler, make, cmake, git).  Every archive downloaded is a gzipped tarball
#  and Eigen/Ceres come from git, so nothing here needs `unzip` - which is what lets the manylinux
#  CI container run this without installing a single package.
#
# On Alpine linux, you need at least the following packages
# apk add --no-cache alpine-sdk cmake linux-headers curl git

# If you need to insecurely bypass proxy or something, you can have curl not check the ssl certificate (not recommended)
curl_extra_args=""
# curl_extra_args="--insecure"


_working_dir_arg=$2
_install_dir_arg=$3
# Default to however many cores the build machine has, but stay overridable so it can be dialed
#  back (or up, where nproc under-reports), e.g. `_ncore=16 ./dep_build_linux.sh ...`
_ncore=${_ncore:-$(nproc 2>/dev/null || echo 4)}

if [ ! -d "$1" ]; then
  echo "The first argument (InterSpec source code directory '$1') is not a valid directory."
  exit 1
fi

# Get absolute path for interspec_src_dir
# Assumes $1 is an existing directory due to the check above
interspec_src=$(cd "$1" && pwd)
if [ $? -ne 0 ] || [ -z "$interspec_src" ]; then
    echo "Error: Failed to resolve absolute path for InterSpec source directory '$1'."
    exit 1
fi

PATCH_DIR="${interspec_src}/target/patches/"
if [ ! -d "${PATCH_DIR}" ]; then
  echo "The InterSpec patch directory at '${PATCH_DIR}' is invalid."
  exit 1
fi

# Use working_directory argument for initial path
working_directory="$_working_dir_arg"

# Create the working directory if it doesn't exist
if [ ! -d "$working_directory" ]; then
    mkdir -p "$working_directory"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create working directory '$working_directory'."
        exit 1
    fi
    echo "Created working directory: $working_directory"
else
    echo "Working directory already exists: $working_directory"
fi

# Handle install_directory: if not provided, use working_directory/install
if [ -z "$_install_dir_arg" ]; then
    install_directory="${working_directory}/install"
    echo "Install directory not specified, using default: '$install_directory'"
else
    install_directory="$_install_dir_arg"
fi

# Create the install directory if it doesn't exist
if [ ! -d "$install_directory" ]; then
    mkdir -p "$install_directory"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create install directory '$install_directory'."
        exit 1
    fi
    echo "Created install directory: $install_directory"
else
    echo "Install directory already exists: $install_directory"
fi

# Get absolute paths for working_directory and install_directory
# These directories should exist at this point
_orig_wd="$working_directory"
working_directory=$(cd "$working_directory" && pwd)
if [ $? -ne 0 ] || [ -z "$working_directory" ]; then
    echo "Error: Failed to resolve absolute path for working directory '$_orig_wd'."
    exit 1
fi

_orig_id="$install_directory"
install_directory=$(cd "$install_directory" && pwd)
if [ $? -ne 0 ] || [ -z "$install_directory" ]; then
    echo "Error: Failed to resolve absolute path for install directory '$_orig_id'."
    exit 1
fi


# Output the assigned variables
echo "InterSpec Source Directory: $interspec_src"
echo "Working Directory: $working_directory"
echo "Install Directory: $install_directory"


# git patching will fail if our working directory is a sub-directory of the InterSpec source code
case "$working_directory" in
  "$interspec_src"*) # Check if working_directory starts with interspec_src
    echo "The working directory (second argument) must not be a subdirectory of the InterSpec source directory (first argument)."
    exit 1
    ;;
esac

export MY_WT_PREFIX="$install_directory"

# CMake 4.0 made `cmake_minimum_required` below 3.5 a hard error.  The Wt patch we apply below
#  raises Wt itself, but zlib and Ceres are not covered, so raise the floor for every sub-build
#  here rather than patching each dependency.  The variable is new in CMake 3.31 and ignored
#  before that.  The top-level InterSpec CMakeLists.txt does the same thing for its own build.
export CMAKE_POLICY_VERSION_MINIMUM=3.5

# Build a static OpenSSL into the prefix and enable Wt's SSL support, so Wt::Http::Client can do
#  https:// - which is what the USE_NATIVE_HTTP_CLIENT backend needs on Linux, where the native
#  transport IS Wt::Http::Client.
#
# Default OFF so no existing build changes behaviour, and nobody pays the OpenSSL build cost
#  until they ask for it.  Turn on with `BUILD_OPENSSL=ON ./dep_build_linux.sh ...`.
#
# Must be OpenSSL from the prefix, never the distribution's: linking the system one would tie the
#  resulting binary to that distribution's exact OpenSSL soname.
BUILD_OPENSSL=${BUILD_OPENSSL:-OFF}
OPENSSL_VERSION=3.5.4
echo "BUILD_OPENSSL: ${BUILD_OPENSSL}"

# Guard against two stale-sentinel traps that would silently produce a broken prefix:
#
#  - `wt.installed` from an earlier BUILD_OPENSSL=OFF run makes the Wt block below a no-op, so Wt
#    keeps ENABLE_SSL=OFF while the prefix now *does* contain libssl.a - which looks like success
#    to InterSpec's cmake, but leaves https dead at runtime.
#  - `wt-3.7.1/wt.patched` from before this script gained the WtFindSsl.txt hunk means that hunk is
#    never applied, so Wt links a bare "-lcrypto" (the system one) against our static libssl.a.
#
# Both are silent, so refuse rather than guess.
if [ "${BUILD_OPENSSL}" = "ON" ] && [ -f "${working_directory}/wt.installed" ]; then
  if [ ! -f "${MY_WT_PREFIX}/include/Wt/WConfig.h" ] \
     || ! grep -q "^#define WT_WITH_SSL" "${MY_WT_PREFIX}/include/Wt/WConfig.h"; then
    echo ""
    echo "ERROR: BUILD_OPENSSL=ON, but this working directory already has a Wt built WITHOUT SSL."
    echo "       Wt will not be rebuilt while its sentinel exists, so the result would silently"
    echo "       have no https support.  Remove these and re-run:"
    echo "         rm ${working_directory}/wt.installed"
    echo "         rm -rf ${working_directory}/wt-3.7.1"
    echo "       (deleting the extracted Wt source also drops its 'wt.patched' marker, so the"
    echo "        WtFindSsl.txt fix gets applied)"
    echo ""
    exit 1
  fi
fi


# Define a function to download a file and check its hash
download_file() {
  _file_url="$1"
  _file_name="$2"
  _expected_sha256="$3"
  _downloaded_sha256="" # Initialize

  if [ -f "$_file_name" ]; then
    echo "File '$_file_name' already exists. Skipping download."
  else
    # --retry/--retry-connrefused because these are the only network steps in the whole build and
    #  a transient failure otherwise fails CI.  --fail so an HTTP error page is not silently saved
    #  as the "archive" (the sha256 check would catch it, but the error message would be useless).
    curl ${curl_extra_args} -fL --retry 5 --retry-delay 5 --retry-connrefused \
         --connect-timeout 30 "$_file_url" --output "$_file_name"
  fi

  # Ensure sha256sum is available (external dependency)
  if ! command -v sha256sum > /dev/null; then
    echo "Error: sha256sum command not found. Please install it."
    exit 1
  fi
  _downloaded_sha256=$(sha256sum "$_file_name" | awk '{ print $1 }')

  if [ "$_downloaded_sha256" != "$_expected_sha256" ]; then
    echo "Error: SHA-256 hash mismatch for '$_file_name'."
    echo "Expected: $_expected_sha256"
    echo "Got: $_downloaded_sha256"
    exit 1
  else
    echo "SHA-256 hash verified for '$_file_name'."
  fi
}


# Change to the working directory
cd "$working_directory"
if [ $? -ne 0 ]; then
    echo "Error: Failed to change to working directory '$working_directory'."
    exit 1
fi


## Build boost
if [ -f "${working_directory}/boost.installed" ]; then
    echo "Boost already installed (as indicated by existence of boost.installed file) - skipping."
else
  # Use the .tar.gz rather than the .zip: it needs no `unzip`, is a smaller download, and is the
  #  Unix-flavoured archive (the .zip is the Windows one, with CRLF line endings - which
  #  ./bootstrap.sh below does not appreciate).  Hashes are the ones published alongside the
  #  release, e.g. https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.gz.json
  #  archives.boost.io rather than sourceforge: sourceforge bounces through two redirects to a
  #  randomly chosen mirror, and when that mirror is unhealthy curl just stalls and exits 7 -
  #  observed while testing this script.  archives.boost.io is Boost's own archive and serves the
  #  file directly, with no redirect and no mirror roulette.
  _file_url="https://archives.boost.io/release/1.84.0/source/boost_1_84_0.tar.gz"
  _file_name="boost_1_84_0.tar.gz"
  _expected_sha256="a5800f405508f5df8114558ca9855d2640a2de8f0445f051fa1c7c3383045724"
  _src_dir="boost_1_84_0"
  #
  # Wt fails to compile against boost 1.85, but you just need to modify:
  #  - wt-3.7.1/src/web/FileUtils.C to include boost/filesystem.hpp
  #  - wt-3.7.1/src/http/Configuration.h to change `bool hasSslPasswordCallback()` to be { return !sslPasswordCallback_.empty(); }
  #_file_url="https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.gz"
  #_file_name="boost_1_85_0.tar.gz"
  #_expected_sha256="be0d91732d5b0cc6fbb275c7939974457e79b54d6f07ce2e3dfdd68bef883b0b"
  #_src_dir="boost_1_85_0"

  download_file "${_file_url}" "${_file_name}" "${_expected_sha256}"

  # untar and change into resulting directory
  if [ -d "${_src_dir}" ]; then
    echo "Boost already un-tarred, not doing again."
  else
    tar -xzf "${_file_name}" # tar is an external dependency
  fi

  cd "${_src_dir}"

  if [ -e "${working_directory}/boost.built" ]; then
    echo "Boost already built - not rebuilding."
  else
    if [ -e "b2" ]; then
      echo "Boost bootstrap looks to be done - not doing again"
    else
      # build the b2 executable
      ./bootstrap.sh --prefix="${MY_WT_PREFIX}"
    fi # if b2 already built / else

    ./b2 -j${_ncore} --without-python -q cxxflags="-std=c++17 -fPIC -fvisibility=default" linkflags="-std=c++17 -fPIC -fvisibility=default" cflags=-fPIC link=static variant=release threading=multi --build-dir=linux_build --prefix="${MY_WT_PREFIX}" -a install

    touch "${working_directory}/boost.built"
  fi

  touch "${working_directory}/boost.installed"
fi #if boost.installed exists / else


cd "${working_directory}"


## Build zlib
if [ -f "${working_directory}/zlib.installed" ]; then
    echo "zlib already installed (as indicated by existence of zlib.installed file) - skipping."
else
  _file_url="https://zlib.net/fossils/zlib-1.3.1.tar.gz"
  _file_name="zlib-1.3.1.tar.gz"
  _expected_sha256="9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23"
  _src_dir="zlib-1.3.1"

  download_file "${_file_url}" "${_file_name}" "${_expected_sha256}"

  if [ -d "${_src_dir}" ]; then
    echo "zlib already unzipped, not doing again."
  else
    tar -xzvf "${_file_name}" # tar is an external dependency
  fi

  cd "${_src_dir}"

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous zlib build directory."
  fi

  mkdir build
  cd build
  
  # CMAKE_INSTALL_LIBDIR=lib is set on every CMake sub-build here so the whole prefix keeps one
  #  layout.  boost's b2 and Wt always install into `lib`, but any dependency using GNUInstallDirs
  #  would pick `lib64` on a RHEL-family distro, and two places assume `lib`: the
  #  `rm -f .../lib/libz.so*` below, and cmake/FindWt.cmake's `HINTS ${Wt_INCLUDE_DIR}/../lib`.
  #  (zlib 1.3.1 has its own INSTALL_LIB_DIR that already defaults to lib, so here it is belt and
  #  braces; for Ceres it is doing real work.)
  cmake -DCMAKE_BUILD_TYPE=Release -DZLIB_BUILD_EXAMPLES=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" .. # cmake is an external dependency
  make -j${_ncore} install # make is an external dependency
  rm -rf ./*
  # Delete the shared zlib so everything downstream links the static libz.a.  A shared libz would
  #  become `DT_NEEDED libz.so.1` on the Electron addon, which target/electron/linux/check_elf_compat.py
  #  rejects: it is not on the allowlist, because we cannot assume end-user systems have it.
  rm -f "${MY_WT_PREFIX}"/lib/libz.so*

  touch "${working_directory}/zlib.installed"
fi #if zlib.installed exists / else


cd "${working_directory}"


## Build a static OpenSSL, if asked for (see BUILD_OPENSSL near the top of this file)
if [ "${BUILD_OPENSSL}" != "ON" ]; then
  echo "BUILD_OPENSSL is not ON - skipping OpenSSL, and Wt will be built without SSL support."
elif [ -f "${working_directory}/openssl.installed" ]; then
  echo "OpenSSL already installed (as indicated by existance of openssl.installed file) - skipping."
else
  file_url="https://github.com/openssl/openssl/releases/download/openssl-${OPENSSL_VERSION}/openssl-${OPENSSL_VERSION}.tar.gz"
  file_name="openssl-${OPENSSL_VERSION}.tar.gz"
  expected_sha256="967311f84955316969bdb1d8d4b983718ef42338639c621ec4c34fddef355e99"
  src_dir="openssl-${OPENSSL_VERSION}"

  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  if [ -d "${src_dir}" ]; then
    echo "OpenSSL already unzipped, not doing again."
  else
    tar -xzf "${file_name}"
  fi

  cd "${src_dir}"

  # -fPIC because everything else in this prefix is built that way and PIE executables are the
  #  default on Linux.  no-shared keeps it static, so nothing extra ships alongside InterSpec.
  # --openssldir gets compiled into the library and points into this build prefix, which will not
  #  exist on a user's machine - which is why the Wt backend calls setSslVerifyPath() explicitly
  #  rather than relying on OpenSSL's defaults or on SSL_CERT_DIR.
  # OpenSSL needs an explicit target triple; pick it from the build machine rather than assuming
  #  x86_64, since Linux ARM is common now (an aarch64 container on an Apple-silicon host, or an
  #  ARM server).  Getting this wrong produces a confusing compile failure, not a clear message.
  case "$(uname -m)" in
    x86_64)          ossl_target="linux-x86_64" ;;
    aarch64 | arm64) ossl_target="linux-aarch64" ;;
    armv7l | armv7)  ossl_target="linux-armv4" ;;
    ppc64le)         ossl_target="linux-ppc64le" ;;
    s390x)           ossl_target="linux64-s390x" ;;
    *)
      echo "Unrecognized machine type '$(uname -m)' - cannot pick an OpenSSL target."
      exit 1
      ;;
  esac
  echo "Building OpenSSL for target ${ossl_target}"

  ./Configure "${ossl_target}" no-shared no-tests no-docs no-legacy -fPIC \
    --prefix="${MY_WT_PREFIX}" \
    --openssldir="${MY_WT_PREFIX}/ssl" \
    --libdir=lib

  make -j${_ncore}
  make install_sw

  cd "${working_directory}"
  touch "${working_directory}/openssl.installed"
fi #if BUILD_OPENSSL / openssl.installed exists / else


cd "${working_directory}"



## Build Wt 3.7.1
if [ -f "${working_directory}/wt.installed" ]; then
    echo "Wt already installed (as indicated by existence of wt.installed file) - skipping."
else
  _file_url="https://github.com/emweb/wt/archive/3.7.1.tar.gz"
  _file_name="wt-3.7.1.tar.gz"
  _expected_sha256="232a2f0a2f3c4174c96872ab15ef7bc0a544d9401486a4c496a6b37f27cc10e7"
  _src_dir="wt-3.7.1"

  download_file "${_file_url}" "${_file_name}" "${_expected_sha256}"

  if [ -d "${_src_dir}" ]; then
    echo "Wt already unzipped, not doing again."
  else
    tar -xzvf "${_file_name}"
  fi

  cd "${_src_dir}"

  if [ -f "wt.patched" ]; then
    echo "Wt already patched, not doing again."
  else
    # git is an external dependency
    if ! command -v git > /dev/null; then
      echo "Error: git command not found. Please install it to apply patches."
      exit 1
    fi
    git apply --reject --ignore-space-change --ignore-whitespace "${PATCH_DIR}/wt/3.7.1/NormalBuild/wt_3.7.1_git.patch"
    touch wt.patched
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Wt build directory."
  fi

  mkdir build
  cd build

  # Wt 3.7.1 finds OpenSSL with its own cmake/WtFindSsl.txt, keyed on SSL_PREFIX - NOT with
  #  find_package(OpenSSL), so OPENSSL_ROOT_DIR would be silently ignored here.  The patch applied
  #  above also fixes that file to link libcrypto.a out of the prefix instead of emitting a bare
  #  "-lcrypto", which would otherwise pull in the *system* libcrypto alongside our static libssl.
  if [ "${BUILD_OPENSSL}" = "ON" ]; then
    WT_SSL_ARGS="-DENABLE_SSL=ON -DSSL_PREFIX=${MY_WT_PREFIX}"
  else
    WT_SSL_ARGS="-DENABLE_SSL=OFF"
  fi

  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DWT_CMAKE_FINDER_INSTALL_DIR="${MY_WT_PREFIX}/share" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBoost_INCLUDE_DIR="${MY_WT_PREFIX}/include" -DBOOST_PREFIX="${MY_WT_PREFIX}" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DHARU_PREFIX="${MY_WT_PREFIX}" ${WT_SSL_ARGS} -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++17" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCONFIGURATION=data/config/wt_config_web.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="${MY_WT_PREFIX}/etc/wt" -S ..
  #make -j${_ncore} install
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/wt.installed"
fi #if wt.installed exists / else


cd "${working_directory}"


## Install Eigen
if [ -f "${working_directory}/Eigen.installed" ]; then
    echo "Eigen already installed (as indicated by existance of Eigen.installed file) - skipping."
else
  # Build Eigen, which is required by ceres-solver, and used a few other places
  # in InterSpec if its available

  src_dir="eigen-3.x"
  # Get trunk version as of 20250114 to pickup some compile issues for c++20 (minimizing how much history we download)
  git_hash="2e76277bd049f7bec36b0f908c69734a42c5234f"
  
  if [ -d "${src_dir}" ]; then
    echo "Eigen cloned - not doing it again."
    cd "${src_dir}"
  else
    git clone --recursive https://gitlab.com/libeigen/eigen.git --branch master --single-branch --depth 1 "${src_dir}"
    cd "${src_dir}"
    git fetch --depth 1 origin ${git_hash}
    git checkout ${git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Eigen build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_LIBDIR=lib -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/Eigen.installed"
fi #if Eigen.installed exists / else


cd "${working_directory}"


## Build Ceres-Solver
if [ -f "${working_directory}/Ceres.installed" ]; then
    echo "Ceres-Solver already installed (as indicated by existence of Ceres.installed file) - skipping."
else
  # Build ceres-solver; this is the optimizer used for the relative efficiency
  # tool, and a small amount of the peak fitting.
  _src_dir="ceres-solver"
  # Get version 2.2.0, Oct 12, 2023 (minimizing how much history we download)
  _git_hash="85331393dc0dff09f6fb9903ab0c4bfa3e134b01"

  if [ -d "${_src_dir}" ]; then
    echo "Ceres-Solver already cloned - not doing it again."
    cd "${_src_dir}"
  else
    git clone --recursive https://github.com/ceres-solver/ceres-solver.git --branch master --single-branch --depth 1 "${_src_dir}"
    cd "${_src_dir}"
    git fetch --depth 1 origin ${_git_hash}
    git checkout ${_git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build_linux ]; then
    rm -r build_linux
    echo "Deleted previous Ceres-Solver build_linux directory."
  fi

  mkdir build_linux
  cd build_linux

  # BUILD_TESTING=OFF because `--target install` depends on `all`, and Ceres puts its whole test
  #  suite in `all` - that is dozens of test binaries compiled and then thrown away.
  # EXPORT_BUILD_DIR=OFF because ON registers this build tree in ~/.cmake/packages/Ceres, so a
  #  later find_package(Ceres) could resolve to a scratch directory instead of the prefix.
  #  (cmake/FetchInterSpecDeps.cmake sets both OFF for the same reasons.)
  cmake -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_LIBDIR=lib -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DUSE_CUDA=OFF -DEXPORT_BUILD_DIR=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF ..
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/Ceres.installed"
fi #if Ceres.installed exists / else

cd "${working_directory}"

echo "Have successfully installed all libraries to ${install_directory}"