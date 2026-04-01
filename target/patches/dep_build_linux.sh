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

# On Alpine linux, you need at least the following packages
# apk add --no-cache alpine-sdk cmake linux-headers curl git

# If you need to insecurely bypass proxy or something, you can have curl not check the ssl certificate (not recommended)
curl_extra_args=""
# curl_extra_args="--insecure"


_working_dir_arg=$2
_install_dir_arg=$3
_ncore=4

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


# Define a function to download a file and check its hash
download_file() {
  _file_url="$1"
  _file_name="$2"
  _expected_sha256="$3"
  _downloaded_sha256="" # Initialize

  if [ -f "$_file_name" ]; then
    echo "File '$_file_name' already exists. Skipping download."
  else
    curl ${curl_extra_args} -L "$_file_url" --output "$_file_name"
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
  _file_url="https://sourceforge.net/projects/boost/files/boost/1.84.0/boost_1_84_0.zip/download"
  _file_name="boost_1_84_0.zip"
  _expected_sha256="cc77eb8ed25da4d596b25e77e4dbb6c5afaac9cddd00dc9ca947b6b268cc76a4"
  #
  # Wt fails to compile against boost 1.85, but you just need to modify:
  #  - wt-3.7.1/src/web/FileUtils.C to include boost/filesystem.hpp
  #  - wt-3.7.1/src/http/Configuration.h to change `bool hasSslPasswordCallback()` to be { return !sslPasswordCallback_.empty(); }
  #_file_url="https://sourceforge.net/projects/boost/files/boost/1.85.0/boost_1_85_0.zip/download"
  #_file_name="boost_1_85_0.zip"
  #_expected_sha256="e712fe7eb1b9ec37ac25102525412fb4d74e638996443944025791f48f29408a"

  _src_dir="${_file_name%.*}" # POSIX compliant variable expansion

  download_file "${_file_url}" "${_file_name}" "${_expected_sha256}"

  # unzip and change into resulting directory
  if [ -d "${_src_dir}" ]; then
    echo "Boost already unzipped, not doing again."
  else
    unzip "${_file_name}" # unzip is an external dependency
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


## Build zlib 1.3.2
if [ -f "${working_directory}/zlib.installed" ]; then
    echo "zlib already installed (as indicated by existence of zlib.installed file) - skipping."
else
  _src_dir="zlib-1.3.2"
  # zlib version 1.3.2
  _git_hash="da607da739fa6047df13e66a2af6b8bec7c2a498"

  if [ -d "${_src_dir}" ]; then
    echo "zlib cloned - not doing it again."
    cd "${_src_dir}"
  else
    git clone --recursive https://github.com/madler/zlib.git --branch master --single-branch --depth 1 "${_src_dir}"
    cd "${_src_dir}"
    git fetch --depth 1 origin ${_git_hash}
    git checkout ${_git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous zlib build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_BUILD_TYPE=Release -DZLIB_BUILD_EXAMPLES=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" .. # cmake is an external dependency
  make -j${_ncore} install # make is an external dependency
  rm -rf ./*
  rm -f "${MY_WT_PREFIX}/lib/libz.so"
  rm -f "${MY_WT_PREFIX}/lib/libz.so.1"
  rm -f "${MY_WT_PREFIX}/lib/libz.so.1.3.2"

  touch "${working_directory}/zlib.installed"
fi #if zlib.installed exists / else


cd "${working_directory}"


## Build Wt 4.12.6
if [ -f "${working_directory}/wt.installed" ]; then
    echo "Wt already installed (as indicated by existence of wt.installed file) - skipping."
else
  _src_dir="wt-4.12.6"
  _git_hash="4d158458b06554192c63389e098ad03e5e35a713"

  if [ -d "${_src_dir}" ]; then
    echo "Wt cloned - not doing it again."
    cd "${_src_dir}"
  else
    git clone --recursive https://github.com/emweb/wt.git --branch master --single-branch --depth 1 "${_src_dir}"
    cd "${_src_dir}"
    git fetch --depth 1 origin ${_git_hash}
    git checkout ${_git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Wt build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBoost_INCLUDE_DIR="${MY_WT_PREFIX}/include" -DBOOST_PREFIX="${MY_WT_PREFIX}" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DHARU_PREFIX="${MY_WT_PREFIX}" -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DENABLE_FIREBIRD=OFF -DENABLE_MSSQLSERVER=OFF -DENABLE_OPENGL=ON -DENABLE_QT4=OFF -DENABLE_QT5=OFF -DENABLE_QT6=OFF -DENABLE_LIBWTTEST=ON -DENABLE_SAML=OFF -DENABLE_UNWIND=OFF -DHTTP_WITH_ZLIB=OFF -DCMAKE_CXX_STANDARD=17 -DCONFIGURATION=data/config/wt_config_web.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="${MY_WT_PREFIX}/etc/wt" -S ..
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

  _src_dir="eigen-5.x"
  # Get version 5.0.1 (minimizing how much history we download)
  _git_hash="bc3b39870ecb690a623a3f49149a358b95c5781d"

  if [ -d "${_src_dir}" ]; then
    echo "Eigen cloned - not doing it again."
    cd "${_src_dir}"
  else
    git clone --recursive https://gitlab.com/libeigen/eigen.git --branch master --single-branch --depth 1 "${_src_dir}"
    cd "${_src_dir}"
    git fetch --depth 1 origin ${_git_hash}
    git checkout ${_git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Eigen build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/Eigen.installed"
fi #if Eigen.installed exists / else


cd "${working_directory}"

## Build Abseil
if [ -f "${working_directory}/abseil.installed" ]; then
    echo "Abseil already installed (as indicated by existence of abseil.installed file) - skipping."
else
  # Abseil is required by ceres-solver.
  # Using LTS 20260107.1
  _src_dir="abseil-cpp"
  _git_hash="255c84dadd029fd8ad25c5efb5933e47beaa00c7"

  if [ -d "${_src_dir}" ]; then
    echo "Abseil cloned - not doing it again."
    cd "${_src_dir}"
  else
    git clone --recursive https://github.com/abseil/abseil-cpp.git --branch master --single-branch --depth 1 "${_src_dir}"
    cd "${_src_dir}"
    git fetch --depth 1 origin ${_git_hash}
    git checkout ${_git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Abseil build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DABSL_BUILD_TESTING=OFF -DABSL_USE_GOOGLETEST_HEAD=OFF -DCMAKE_CXX_STANDARD=17 -DBUILD_SHARED_LIBS=OFF ..
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/abseil.installed"
fi #if abseil.installed exists / else


cd "${working_directory}"


## Build Ceres-Solver
if [ -f "${working_directory}/Ceres.installed" ]; then
    echo "Ceres-Solver already installed (as indicated by existence of Ceres.installed file) - skipping."
else
  # Build ceres-solver; this is the optimizer used for the relative efficiency
  # tool, and a small amount of the peak fitting.
  git clone --recursive https://github.com/ceres-solver/ceres-solver.git --branch master --single-branch --depth 1
  cd ceres-solver
  # Get HEAD version from Mar 22, 2026 (minimizing how much history we download)
  git fetch --depth 1 origin 2f946a582ae4a9e7ee0492030ec12d9b1f3dbade
  git checkout 2f946a582ae4a9e7ee0492030ec12d9b1f3dbade
  git submodule update --init --recursive

  # Remove bundled abseil so Ceres uses the system-installed one (avoids CMake
  # export-set errors when bundled abseil targets arent in the CeresExport set)
  rm -rf third_party/abseil-cpp/*

  if [ -d build_linux ]; then
    rm -r build_linux
    echo "Deleted previous Ceres-Solver build_linux directory."
  fi

  mkdir build_linux
  cd build_linux

  cmake -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DACCELERATESPARSE=OFF -DUSE_CUDA=OFF -DEXPORT_BUILD_DIR=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_STANDARD=17 ..
  cmake --build . --config Release --target install --parallel ${_ncore}

  touch "${working_directory}/Ceres.installed"
fi #if Ceres.installed exists / else

cd "${working_directory}"

echo "Have successfully installed all libraries to ${install_directory}"