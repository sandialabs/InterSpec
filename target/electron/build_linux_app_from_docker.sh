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
# This script is intended to be run from within the `quay.io/pypa/manylinux2014_x86_64` Docker image, for example, to build InterSpec you could
#  cd /tmp
#  mkdir build_interspec
#  cd build_interspec
#  git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_code
#  mkdir build_electron
#  mkdir build_working_dir
#  docker run --rm -it -v `pwd`/InterSpec_code:/interspec -v `pwd`/build_electron/:/build_app -v `pwd`/build_working_dir:/build_working_dir quay.io/pypa/manylinux_2_28_x86_64:latest /interspec/target/electron/Docker/build_linux_app_from_docker.sh /interspec /build_app /build_working_dir
#  #Then results will then be in build_working_dir/InterSpec-linux-x64

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
yum update
yum install -y npm zip
# Use the 'n' package to install a fairly modern version of npm, by default we are on like version 6 or 8, which is too old to run some of the packages we need
npm install -g n
n 20.18.0
# Just in case, make sure to pick up latest node version
hash -r
# The default version of npm in manylinux_2_28_x86_64 is 6.14.11, which doesnt support the `exec` command, so we'll locally install a newer version of npm
npm install -g npm@10.9.0
hash -r
npm install -g uglify-js
npm install -g uglifycss
npm install -g cmake-js

echo "Will install local npm packages"
npm install --save-dev node-addon-api --arch=x64
npm install electron --arch=x64
npm install --save-dev electron-packager

echo "CWD"
pwd
echo "ls node_modules"
ls node_modules

# We need to help find where node_modules directory is located in the Docker container
export NODE_MODULES_PARENT_DIR=$PWD


echo "Will build InterSpec code"
CMAKE_BUILD_PARALLEL_LEVEL=`nproc` cmake-js --directory ${InterSpecCodePath}/target/electron  --architecture x64 --arch=x64 --CDCMAKE_BUILD_TYPE="Release" --CDInterSpec_NODE_MODULE_DIR="${WorkingDir}/node_modules" --CDInterSpec_FETCH_DEPENDENCIES=ON --CDBUILD_AS_LOCAL_SERVER=OFF --CDCMAKE_SHARED_LINKER_FLAGS="-static-libgcc -static-libstdc++" --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="${LEAFLET_KEY}" --CDUSE_REL_ACT_TOOL=ON --out="${CmakeBuildDir}" --target install

echo "Will package InterSpec code"
if [ -d "${WorkingDir}/InterSpec-linux-x64" ]; then
  rm -rf "${WorkingDir}/InterSpec-linux-x64"
fi

echo "ls"
ls

echo "ls ${CmakeBuildDir} (CmakeBuildDir)"
ls ${CmakeBuildDir}

echo "ls ${WorkingDir} (WorkingDir)"
ls ${WorkingDir}
cd ${WorkingDir}
echo "Changed to: $(pwd)"

if [ -d "${WorkingDir}/app" ]; then
  echo "app was a valid working directory"
  rm -rf "${WorkingDir}/app"
else
  echo "There was no app directory in cwd"
fi

echo "Will try to copy app dir from ${CmakeBuildDir} to $(pwd)."
cp -r "${CmakeBuildDir}/app" "${WorkingDir}/"
echo "Have copied app dir from ${CmakeBuildDir} to ${WorkingDir}/."

echo "About to package"

npm exec electron-packager app InterSpec --overwrite=true --platform=linux --arch=x64 --protocol=interspec --protocol-name="InterSpec" --icon="${InterSpecCodePath}/target/electron/linux/InterSpec_desktop_icon_256x256.png" --prune=true --ignore=LICENSE.md --ignore=README.md

cp "${InterSpecCodePath}/NOTICE.html" "${WorkingDir}/InterSpec-linux-x64/"
echo "This is a build of InterSpec for Linux, using InterSpec code git hash ${GIT_HASH}." > "${WorkingDir}/InterSpec-linux-x64/build_information.txt"
echo "This build is untested - please contact InterSpec@sandia.gov for support." >> "${WorkingDir}/InterSpec-linux-x64/build_information.txt"

echo "About to zip"

zip -r InterSpec_app_Linux_Electron_latest_$(date +%Y-%m-%d).zip InterSpec-linux-x64

echo "Done"