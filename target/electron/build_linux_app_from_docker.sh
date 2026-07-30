#!/bin/bash

handle_error() {
    echo "An error occurred building InterSpec, on line $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# set -e # Cause script to exit as soon as there is an error

# takes three argument: [path to InterSpec code] [Path to CMake build Directory] [Working directory path]
#  - InterSpec code is mapped to `/interspec`
#  - build directory mapped to `/build_app`
#  - build working directory mapped to /build_working_dir
# For example: ./build_app.sh /interspec /build_app /build_working_dir
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

if [ $# -ne 3 ]; then
  echo "Error: This script requires exactly three arguments."
  echo "    [path to InterSpec code] [Path to CMake build Directory] [Working directory path]"
  exit 1  # Exit with an error code
fi

InterSpecCodePath=$1
CmakeBuildDir=$2
WorkingDir=$3

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
npm install --save-dev node-addon-api --arch=x64

echo "CWD"
pwd
echo "ls node_modules"
ls node_modules

# We need to help find where node_modules directory is located in the Docker container
export NODE_MODULES_PARENT_DIR=$PWD


echo "Will build InterSpec code"
# Note: libstdc++/libgcc are linked *dynamically* - see cmake/CxxRuntimePolicy.cmake.  Passing
#  `-static-libstdc++` here is what caused issue #51 (segfault in dl_init before main()); the
#  gcc-toolset compiler in this image already caps the required GLIBCXX at the AlmaLinux 8
#  base, so a dynamic link keeps the glibc-2.28 floor.
CMAKE_BUILD_PARALLEL_LEVEL=`nproc` cmake-js --directory ${InterSpecCodePath}/target/electron  --architecture x64 --arch=x64 --CDCMAKE_BUILD_TYPE="Release" --CDInterSpec_NODE_MODULE_DIR="${WorkingDir}/node_modules" --CDInterSpec_FETCH_DEPENDENCIES=ON --CDBUILD_AS_LOCAL_SERVER=OFF --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="${LEAFLET_KEY}" --CDUSE_REL_ACT_TOOL=ON --out="${CmakeBuildDir}" --target install

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

echo "Done - app payload is in ${CmakeBuildDir}/app"
