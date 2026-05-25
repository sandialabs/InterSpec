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
PATCH_DIR="${interspec_src}/target/dep_build/patches/"
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

export MACOSX_DEPLOYMENT_TARGET=13.3 # macOS Ventura (2023); required for C++20 libc++ features (std::format, std::ranges, std::jthread, <chrono> calendaring)
export MY_WT_PREFIX=$install_directory
export CMAKE_POLICY_VERSION_MINIMUM=3.5

# Optional: also build wxWidgets (off by default).  wxWidgets is only needed when
# building the target/wxwidgets desktop-app form; the local web-server and Electron
# targets do not require it.  Override by editing this line, or run with:
#   BUILD_WXWIDGETS=1 ./dep_build_macOS.sh <src> <build> <install>
BUILD_WXWIDGETS=${BUILD_WXWIDGETS:-0}
echo "BUILD_WXWIDGETS=${BUILD_WXWIDGETS} (set to 1 to also build wxWidgets)"


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
    ./b2 toolset=clang-darwin target-os=darwin architecture=arm abi=aapcs cxxflags="-stdlib=libc++ -arch arm64 -std=c++20 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" cflags="-arch arm64  -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" linkflags="-stdlib=libc++ -arch arm64 -std=c++20 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" link=static variant=release threading=multi --build-dir=macOS_arm64_build --prefix=${MY_WT_PREFIX} -a stage

    # copy arm libraries to a separate directory
    mkdir -p arm64 && cp stage/lib/libboost_* arm64/

    # build boost for x86_64 and install it (we'll copy over the libraries later)
    ./b2 toolset=clang-darwin target-os=darwin architecture=x86 cxxflags="-stdlib=libc++ -arch x86_64 -std=c++20 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" cflags="-arch x86_64 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" linkflags="-stdlib=libc++ -arch x86_64 -std=c++20 -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}" abi=sysv binary-format=mach-o link=static variant=release threading=multi --build-dir=macOS_x64_build --prefix=${MY_WT_PREFIX} -a install

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


## Build Wt 4.12.6
if [ -f "${working_directory}/wt.installed" ]; then
    echo "Wt already installed (as indicated by existance of wt.installed file) - skipping."
else
  src_dir="wt-4.12.6"
  git_hash="4d158458b06554192c63389e098ad03e5e35a713"

  if [ -d "${src_dir}" ]; then
    echo "Wt cloned - not doing it again."
    cd "${src_dir}"
  else
    git clone --recursive https://github.com/emweb/wt.git --branch master --single-branch --depth 1 "${src_dir}"
    cd "${src_dir}"
    git fetch --depth 1 origin ${git_hash}
    git checkout ${git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Wt build directory."
  fi

  mkdir build
  cd build

  # Force Wt to use its bundled Howard Hinnant date library for time-zone support.
  # Apples libc++ does not yet ship a complete C++20 <chrono> tzdb (no std::chrono::time_zone),
  # so the default WT_CPP20_DATE_TZ_IMPLEMENTATION=std (auto-selected when CMAKE_CXX_STANDARD>=20) fails to build.
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DBoost_INCLUDE_DIR="${MY_WT_PREFIX}/include" -DBOOST_PREFIX="${MY_WT_PREFIX}" -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DENABLE_HARU=OFF -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DENABLE_FIREBIRD=OFF -DENABLE_MSSQLSERVER=OFF -DENABLE_OPENGL=ON -DENABLE_QT4=OFF -DENABLE_QT5=OFF -DENABLE_QT6=OFF -DENABLE_LIBWTTEST=ON -DENABLE_SAML=OFF -DENABLE_UNWIND=OFF -DHTTP_WITH_ZLIB=OFF -DCMAKE_CXX_STANDARD=20 -DWT_CPP20_DATE_TZ_IMPLEMENTATION=date -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR="${MY_WT_PREFIX}/etc/wt" -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -S ..
  make -j$(sysctl -n hw.ncpu) install
  touch "${working_directory}/wt.installed"
fi #if wt.installed exists / else


cd "${working_directory}"


## Install Eigen
if [ -f "${working_directory}/Eigen.installed" ]; then
    echo "Eigen already installed (as indicated by existance of Eigen.installed file) - skipping."
else
  # Build Eigen, which is required by ceres-solver, and used a few other places
  # in InterSpec if its available

  src_dir="eigen-5.x"
  # Get version 5.0.1 (minimizing how much history we download)
  git_hash="bc3b39870ecb690a623a3f49149a358b95c5781d"
  
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

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
  cmake --build . --config Release --target install -j $(sysctl -n hw.ncpu)

  touch "${working_directory}/Eigen.installed"
fi #if Eigen.installed exists / else


cd "${working_directory}"

## Build Abseil
if [ -f "${working_directory}/abseil.installed" ]; then
    echo "Abseil already installed (as indicated by existence of abseil.installed file) - skipping."
else
  # Abseil is required by ceres-solver.
  # Using LTS 20260107.1
  src_dir="abseil-cpp"
  git_hash="255c84dadd029fd8ad25c5efb5933e47beaa00c7"

  if [ -d "${src_dir}" ]; then
    echo "Abseil cloned - not doing it again."
    cd "${src_dir}"
  else
    git clone --recursive https://github.com/abseil/abseil-cpp.git --branch master --single-branch --depth 1 "${src_dir}"
    cd "${src_dir}"
    git fetch --depth 1 origin ${git_hash}
    git checkout ${git_hash}
    git submodule update --init --recursive
  fi

  if [ -d build ]; then
    rm -r build
    echo "Deleted previous Abseil build directory."
  fi

  mkdir build
  cd build

  cmake -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DABSL_BUILD_TESTING=OFF -DABSL_USE_GOOGLETEST_HEAD=OFF -DCMAKE_CXX_STANDARD=20 -DBUILD_SHARED_LIBS=OFF -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
  cmake --build . --config Release --target install -j $(sysctl -n hw.ncpu)

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

  if [ -d build_macos ]; then
    rm -r build_macos
    echo "Deleted previous Ceres-Solver build_macos directory."
  fi

  mkdir build_macos
  cd build_macos

  cmake -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" -DACCELERATESPARSE=OFF -DUSE_CUDA=OFF -DEXPORT_BUILD_DIR=OFF -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_STANDARD=20 -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
  # Build single-threaded to avoid lipo race conditions with universal binary builds
  cmake --build . --config Release --target install -j $(sysctl -n hw.ncpu)

  touch "${working_directory}/Ceres.installed"
fi #if Ceres.installed exists / else

cd "${working_directory}"


## Build wxWidgets 3.2.1 (optional; needed only by the target/wxwidgets desktop-app form).
## Matches the version used by target/dep_build/dep_build_msvc2022.bat so the API surface stays
## consistent between platforms.
if [ "${BUILD_WXWIDGETS}" != "1" ]; then
    echo "Skipping wxWidgets build (set BUILD_WXWIDGETS=1 to enable)."
elif [ -f "${working_directory}/wxWidgets.installed" ]; then
    echo "wxWidgets already installed (as indicated by existence of wxWidgets.installed file) - skipping."
else
  file_url="https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.1/wxWidgets-3.2.1.zip"
  file_name="wxWidgets-3.2.1.zip"
  expected_sha256="cc8868c3c8ec4eddaf659a8b81589a3d83126d5afde012350f61031a607a56d8"
  src_dir="wxWidgets-3.2.1"
  # The wxWidgets source tree contains its own top-level "build/" directory
  # (build/cmake/main.cmake, build/bakefiles/version.bkl, ...) so we must not
  # use that name as our CMake binary directory or we'd clobber required sources.
  build_dir="wxWidgets-3.2.1-build"

  # The release zip does not contain a top-level directory, so create one and extract into it
  # (mirrors how target/dep_build/dep_build_msvc2022.bat lays things out).
  if [ -d "${src_dir}" ]; then
    echo "wxWidgets already extracted - not doing again."
    cd "${src_dir}"
  else
    mkdir -p "${src_dir}"
    cd "${src_dir}"
    download_file "${file_url}" "${file_name}" "${expected_sha256}"
    unzip -q "${file_name}"

    # Patch bundled libpng 1.6.x: it #include's the Classic Mac OS <fp.h> when
    # TARGET_OS_MAC is defined, which fails on modern macOS.  Upstream libpng
    # dropped this branch in 1.6.40 -- do the equivalent fix in-place.
    /usr/bin/sed -i '' 's|defined(__SC__) || defined(TARGET_OS_MAC)|defined(__SC__)|' src/png/pngpriv.h
  fi

  cd "${working_directory}"
  if [ -d "${build_dir}" ]; then
    rm -r "${build_dir}"
    echo "Deleted previous wxWidgets build directory."
  fi

  mkdir "${build_dir}"
  cd "${build_dir}"

  # Universal binary, static link, minimal extras.  On macOS wxWidgets' webview backend is WebKit
  # (system framework); wxUSE_WEBVIEW and wxUSE_STC are needed by target/wxwidgets'
  # find_package(wxWidgets REQUIRED core base net webview stc).
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET}" \
        -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" \
        -DCMAKE_PREFIX_PATH="${MY_WT_PREFIX}" \
        -DCMAKE_INSTALL_PREFIX="${MY_WT_PREFIX}" \
        -DwxBUILD_SHARED=OFF \
        -DwxBUILD_TESTS=OFF \
        -DwxBUILD_SAMPLES=OFF \
        -DwxBUILD_DEMOS=OFF \
        -DwxUSE_WEBVIEW=ON \
        -DwxUSE_STC=ON \
        "${working_directory}/${src_dir}"
  cmake --build . --config Release --target install -j $(sysctl -n hw.ncpu)

  touch "${working_directory}/wxWidgets.installed"
fi #if wxWidgets.installed exists / else

cd "${working_directory}"


echo "Have successfully installed all libraries to ${install_directory}"
