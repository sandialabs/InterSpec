#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check the number of arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <interspec_src_dir> <working_directory> <install_directory>"
    exit 1
fi

# Set proxy, if you need it
#export http_proxy=http://proxy.sandia.gov:80
#export https_proxy=http://proxy.sandia.gov:80

working_directory=$2
install_directory=$3

if [ ! -d "$1" ]; then
  echo "The first argument to the InterSpec source code directory is not valid"
  exit 1
fi

interspec_src=$(realpath "$1")
PATCH_DIR="${interspec_src}/target/patches/"
if [ ! -d "${PATCH_DIR}" ]; then
  echo "The InterSpec patch directory at '${PATCH_DIR}' is invalid."
  exit 1
fi

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

# Assign arguments to variables
working_directory=$(realpath "$working_directory")
install_directory=$(realpath "$install_directory")

# Output the assigned variables
echo "Working Directory: $working_directory"
echo "Install Directory: $install_directory"


# git patching will fail if our working directory is a sub-directory of the InterSpec source code
if [[ "$working_directory" == "$interspec_src"* ]]; then
  echo "The working directory (second argument) must not be a subdirectory of the InterSpec source directory (first argument)."
  exit 1
fi

export MACOSX_DEPLOYMENT_TARGET=10.13 # MacOS High Sierra (2017, supported through 2020)
export MY_WT_PREFIX=$install_directory


# Define a function to download a file and check its hash
download_file() {
  local file_url="$1"
  local file_name="$2"
  local expected_sha256="$3"

  if [ -f "$file_name" ]; then
    echo "File '$file_name' already exists. Skipping download."
  else
    curl -L "${file_url}" --output "${file_name}"
  fi

  downloaded_sha256=$(shasum -a 256 "$file_name" | awk '{ print $1 }')
  if [ "$downloaded_sha256" != "$expected_sha256" ]; then
    echo "Error: SHA-256 hash mismatch for '$file_name'."
    echo "Expected: $expected_sha256"
    echo "Got: $downloaded_sha256"
    exit 1
  else
    echo "SHA-256 hash verified for '$file_name'."
  fi
}


# Change to the working directory
cd "$working_directory"
if [ $? -ne 0 ]; then
    echo "Error: Failed to change to working directory '$working_directory'."
    exit 1
fi


## Build boost 1.78
if [ -f "${working_directory}/boost.installed" ]; then
    echo "Boost already installed (as indicated by existance of boost.installed file) - skipping."
else
  # Download the file using curl
  # Boost 1.78 seems to fail compile with Xcode 16.1 
  #file_url="https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download"
  #file_name="boost_1_78_0.zip"
  #expected_sha256="f22143b5528e081123c3c5ed437e92f648fe69748e95fa6e2bd41484e2986cc3"
  #
  # Boost 1.80 seems to fail compile with Xcode 16.1
  #file_url="https://sourceforge.net/projects/boost/files/boost/1.80.0/boost_1_80_0.zip/download"
  #file_name="boost_1_80_0.zip"
  #expected_sha256="e34756f63abe8ac34b35352743f17d061fcc825969a2dd8458264edb38781782"
  #
  # Boost 1.81 gives "duplicate symbol 'boost::phoenix::placeholders::uarg1" error when compiling InterSpec code
  #file_url="https://sourceforge.net/projects/boost/files/boost/1.81.0/boost_1_81_0.zip/download"
  #file_name="boost_1_81_0.zip"
  #expected_sha256="6e689b266b27d4db57f648b1e5c905c8acd6716b46716a12f6fc73fc80af842e"
  #
  # Boost 1.82 gives "duplicate symbol 'boost::phoenix::placeholders::uarg1" error when compiling InterSpec code
  #file_url="https://sourceforge.net/projects/boost/files/boost/1.82.0/boost_1_82_0.zip/download"
  #file_name="boost_1_82_0.zip"
  #expected_sha256="f7c9e28d242abcd7a2c1b962039fcdd463ca149d1883c3a950bbcc0ce6f7c6d9"
  #
  # The linking issue seems to be fixed in boost 1.84: https://lists.boost.org/Archives/boost/2023/11/255233.php
  file_url="https://sourceforge.net/projects/boost/files/boost/1.84.0/boost_1_84_0.zip/download"
  file_name="boost_1_84_0.zip"
  expected_sha256="cc77eb8ed25da4d596b25e77e4dbb6c5afaac9cddd00dc9ca947b6b268cc76a4"
  #
  # Wt fails to compile against boost 1.85, but you just need to modify:
  #  - wt-3.7.1/src/web/FileUtils.C to include boost/filesystem.hpp
  #  - wt-3.7.1/src/http/Configuration.h to change `bool hasSslPasswordCallback()` to be { return !sslPasswordCallback_.empty(); } 
  #file_url="https://sourceforge.net/projects/boost/files/boost/1.85.0/boost_1_85_0.zip/download"
  #file_name="boost_1_85_0.zip"
  #expected_sha256="e712fe7eb1b9ec37ac25102525412fb4d74e638996443944025791f48f29408a"
  


  src_dir="${file_name%.*}"
  
  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  # unzip and change into resulting directory
  if [ -d "${src_dir}" ]; then
    echo "Boost already unzipped, not doing again."
  else
    unzip "${file_name}"
  fi

  cd "${src_dir}"

  if [ -e "${working_directory}/boost.built" ]; then
    echo "Boost already built - not rebuilding."
  else
    if [ -e "b2" ]; then
      echo "Boost bootstrap looks to be done - not doing again"
    else
      # build the b2 executable
      ./bootstrap.sh cxxflags="-arch x86_64 -arch arm64 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" cflags="-arch x86_64 -arch arm64 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" linkflags="-arch x86_64 -arch arm64 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" --prefix=${MY_WT_PREFIX}
    fi # if b2 already built / else

    # build and stage boost for arm64
    ./b2 toolset=clang-darwin target-os=darwin architecture=arm abi=aapcs cxxflags="-stdlib=libc++ -arch arm64 -std=c++14 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" cflags="-arch arm64  -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" linkflags="-stdlib=libc++ -arch arm64 -std=c++14 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" link=static variant=release threading=multi --build-dir=macOS_arm64_build --prefix=${MY_WT_PREFIX} -a stage

    # copy arm libraries to a separate directory
    mkdir -p arm64 && cp stage/lib/libboost_* arm64/

    # build boost for x86_64 and install it (we'll copy over the libraries later)
    ./b2 toolset=clang-darwin target-os=darwin architecture=x86 cxxflags="-stdlib=libc++ -arch x86_64 -std=c++14 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" cflags="-arch x86_64 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" linkflags="-stdlib=libc++ -arch x86_64 -std=c++14 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" abi=sysv binary-format=mach-o link=static variant=release threading=multi --build-dir=macOS_x64_build --prefix=${MY_WT_PREFIX} -a install

    # move x86 libraries to a seperate directory
    mkdir x86_64 && mv ${MY_WT_PREFIX}/lib/libboost_* x86_64/

    # Now lipo libraries together
    mkdir -p universal
    for dylib in arm64/*; do 
      lipo -create -arch arm64 $dylib -arch x86_64 x86_64/$(basename $dylib) -output universal/$(basename $dylib); 
    done

    touch "${working_directory}/boost.built"
  fi

  # And finally copy the fat libraries to the prefix location
  cp universal/libboost_* ${MY_WT_PREFIX}/lib/

  touch "${working_directory}/boost.installed"
fi #if boost.installed exists / else


cd "${working_directory}"


## Build libpng
if [ -f "${working_directory}/libpng.installed" ]; then
    echo "libpng already installed (as indicated by existance of libpng.installed file) - skipping."
else
  # libpng is only necessary if you are going to build the macOS packaged app, with the Quick Look utility that provides previews of spectrum files in the Finder and various places throughout the OS.
  # We will build once for arm64, and then x86_64 and lipo the libraries together to create a fat library.
  file_url="http://prdownloads.sourceforge.net/libpng/libpng-1.6.37.tar.gz"
  file_name="libpng-1.6.37.tar.gz"
  expected_sha256="daeb2620d829575513e35fecc83f0d3791a620b9b93d800b763542ece9390fb4"
  src_dir="libpng-1.6.37"

  download_file "${file_url}" "${file_name}" "${expected_sha256}"


  if [ -d "${src_dir}" ]; then
    echo "libpng already unzipped, not doing again."
  else
    tar -xzvf "${file_name}"
  fi
  
  cd "${src_dir}"

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous libpng build directory."
  fi

  mkdir build
  cd build

  # First build for arm64 (-DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" doesnt seem to work)
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_OSX_ARCHITECTURES="arm64" -DPNG_ARM_NEON=on ..
  make -j10 install
  rm -rf ./*

  # Then build for x86_64
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_OSX_ARCHITECTURES="x86_64" ..
  make -j10

  # And now lipo the libraries together
  mkdir universal
  lipo -create -arch x86_64 libpng16.a -arch arm64 "${MY_WT_PREFIX}/lib/libpng16.a" -output universal/libpng16.a
  cp universal/libpng16.a "${MY_WT_PREFIX}/lib/"

  touch "${working_directory}/libpng.installed"
fi #if libpng.installe exists / else

cd "${working_directory}"


## Build libharu
if [ -f "${working_directory}/libharu.installed" ]; then
    echo "libharu already installed (as indicated by existance of libharu.installed file) - skipping."
else
  # libharu is only necessary if you are going to build the macOS packaged app, with the Quick Look utility.
  file_url="https://github.com/libharu/libharu/archive/RELEASE_2_3_0.zip"
  file_name="libharu_2.3.0.zip"
  expected_sha256="34de0ed4e994d9c704339292740b8019c33e8538814e104d882508a33517d1a8"
  src_dir="libharu-RELEASE_2_3_0"
   
  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  if [ -d "${src_dir}" ]; then
    echo "libharu already unzipped, not doing again."
  else
    unzip "${file_name}"
  fi
  
  cd "${src_dir}"

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous libharu build directory."
  fi

  mkdir build
  cd build

  # CMake will take care of building for x86_64 and arm64 at the same time
  cmake -DLIBHPDF_SHARED=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DPNG_LIBRARY_RELEASE= -DPNG_LIBRARY_RELEASE=${MY_WT_PREFIX}/lib/libpng.a -DPNG_PNG_INCLUDE_DIR=${MY_WT_PREFIX}/include -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
  make -j10 install

  touch "${working_directory}/libharu.installed"
fi #if libharu.installe exists / else


cd "${working_directory}"


## Build Wt 3.7.1
if [ -f "${working_directory}/wt.installed" ]; then
    echo "Wt already installed (as indicated by existance of wt.installed file) - skipping."
else
  # InterSpec code will only compile with exactly version 3.3.4 or 3.7.1.  
  # For 3.3.4, see version of this file from before Dec 26th 2021 for pathc and build instructions.

  file_url="https://github.com/emweb/wt/archive/3.7.1.tar.gz"
  file_name="wt-3.7.1.tar.gz"
  expected_sha256="232a2f0a2f3c4174c96872ab15ef7bc0a544d9401486a4c496a6b37f27cc10e7"
  src_dir="wt-3.7.1"
  
  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  if [ -d "${src_dir}" ]; then
    echo "Wt already unzipped, not doing again."
  else
    tar -xzvf "${file_name}"
  fi
  
  cd "${src_dir}"

  if [ -f "wt.patched" ]; then
    echo "Wt already pated, not doing again."
  else
    git apply --reject --ignore-space-change --ignore-whitespace "${PATCH_DIR}/wt/3.7.1/NormalBuild/wt_3.7.1_git.patch"
    touch wt.patched
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous libharu build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DBoost_INCLUDE_DIR="${MY_WT_PREFIX}/include" -DBOOST_PREFIX="${MY_WT_PREFIX}" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DHARU_PREFIX="${MY_WT_PREFIX}" -DHARU_LIB="${MY_WT_PREFIX}/lib/libhpdfs.a" -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++14" -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="${MY_WT_PREFIX}/etc/wt" -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -S ..
  make -j10 install
  touch "${working_directory}/wt.installed"
fi #if wt.installed exists / else


cd "${working_directory}"


## Install Eigen
if [ -f "${working_directory}/Eigen.installed" ]; then
    echo "Eigen already installed (as indicated by existance of Eigen.installed file) - skipping."
else
  # Build Eigen, which is required by ceres-solver, and used a few other places
  # in InterSpec if its available
  file_url="https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
  file_name="eigen-3.4.0.tar.gz"
  expected_sha256="8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72"
  src_dir="eigen-3.4.0"
  
  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  if [ -d "${src_dir}" ]; then
    echo "Eigen already unzipped, not doing again."
  else
    tar -xzvf "${file_name}"
  fi

  cd "${src_dir}"

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous libharu build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
  cmake --build . --config Release --target install

  touch "${working_directory}/Eigen.installed"
fi #if Eigen.installed exists / else


cd "${working_directory}"

## Build Ceres-Solver
if [ -f "${working_directory}/Ceres.installed" ]; then
    echo "Ceres-Solver already installed (as indicated by existance of Ceres.installed file) - skipping."
else
  # Build ceres-solver; this is the optimizer used for the relative efficiency
  # tool, and a small amount of the peak fitting.
  file_url="http://ceres-solver.org/ceres-solver-2.1.0.tar.gz"
  file_name="ceres-solver-2.1.0.tar.gz"
  expected_sha256="f7d74eecde0aed75bfc51ec48c91d01fe16a6bf16bce1987a7073286701e2fc6"
  src_dir="ceres-solver-2.1.0"
  
  download_file "${file_url}" "${file_name}" "${expected_sha256}"

  if [ -d "${src_dir}" ]; then
    echo "Ceres-Solver already unzipped, not doing again."
  else
    tar -xzvf "${file_name}"
  fi

  cd "${src_dir}"

  if [ -d build_macos ]; then
    rm -r build_macos
    echo "Deleted previous Ceres-Solver build_macos directory."
  fi

  mkdir build_macos
  cd build_macos

  cmake -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
  cmake --build . --config Release --target install -j 16

  touch "${working_directory}/Ceres.installed"
fi #if Ceres.installed exists / else

cd "${working_directory}"

echo "Have successfully installed all libraries to ${install_directory}"
