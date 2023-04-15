# Creating a Distributable Electron Package

In order to create the [Electron]([https://electronjs.org/](https://electronjs.org/) packaged version of `InterSpec`, we create a node add-on that contains all the C++ code, and then request to start the `InterSpec` server from [main.js](app/main.js).  To create the add-on we use [cmake-js](https://www.npmjs.com/package/cmake-js) to build things, and [node-addon-api](https://www.npmjs.com/package/node-addon-api) to actually interface the C++ to JS.


You can either manually build the InterSpec dependencies (boost, Wt, zlib, and for macOS libpng and libharu), or you can have CMake fetch and build these dependencies.  Having CMake fetch and build the dependencies takes maybe an hour to build things the first time (however, only a couple minutes of this requires your attention) and then subsequent builds are a couple minutes; if you plan to make a substantial number of changes and re-compilations, then the manually built dependencies will yield a slightly faster compile time after each change.

## Building with manually compiled dependencies
To compile the InterSpec code, and package the Electron app, with the manually compiled dependencies (see [patches/README.md](../patches/README.md) for instructions), the following commands are a good base to start with:

```bash
npm install -g cmake-js 

# For macOS only, you may want to define a deployment target
export MACOSX_DEPLOYMENT_TARGET=10.12

cd /path/to/InterSpec/target/electron

# Install dependency for compiling a node.js add-on
npm install --save-dev node-addon-api --arch=x64

# If boost and Wt are in standard locations, you can just run
cmake-js

# Or to have a little more control over things
cmake-js --generator "Visual Studio 16 2019" \
         --architecture x64 --arch=x64 \
         --CDCMAKE_PREFIX_PATH=/path/to/boost/and/wt/and/zlib \
         --CDBoost_USE_STATIC_RUNTIME=ON \
         --CDCMAKE_BUILD_TYPE="Release" \
         --CDGOOGLE_MAPS_KEY="..." \
         --out="build_dir" \
         --target install

# If you make changes and want to recompile
CMAKE_BUILD_PARALLEL_LEVEL=12 cmake-js build --out="build_dir"
# Or you can use CMake to run the `make` command, which can be useful when
# the CMake generator isnt a command-line based system like (ex Xcode, MSVC)
cmake --build build_dir --config Release
# Or directly use the `make` command like:
ninja -C build_dir

# Make the "app" directory into the build directory
ninja -C build_dir install

# And also copy all the InterSpec resources to the 
# 'app' sub-directory of your build dir
cmake-js build --out=build_dir --target install
# Or
cmake --build build_dir --target install --config Release --clean-first
# Or
ninja -C build_dir install


# To run InterSpec without packaging everything, you
# need to install the Electron package
npm install electron --arch=x64

# Then to actually run things
npm start


# To create a packaged version of InterSpec, you
# need to install electron-packager
npm install electron-packager


# Then to actually create the distributable package, run on of the following
# (assuming your build directory is 'build_macos', 'build_win', or 'build_linux', 
#  otherwise you need to edit path in package.json)
npm run package-mac
npm run package-win
npm run package-linux

cp ../../NOTICE.html ./release-builds/InterSpec-PLATFORM-x64/
```

The resulting Electron package seems to run fine on a number of Linux distributions, but to run on the (now out of date) Ubuntu 12.04, you may need to install a couple packages from the package manager on the end-user system:

```bash
apt-get install libnss3
apt-get install libxss1
```

If you are cross-compiling, you can, for example, build the Linux package from macOS using a command like 

```bash
npm --target_arch=x64 --target_platform=linux run package-linux
```


## Building with CMake fetched and compiled dependencies
The following commands will compile and package the InterSpec code, starting from the Fedora 35 Docker image; there is nothing special about Fedora, and any of the distributions compatible with npm, node.js, and Electron should work.

Using the CMake FetchContent will fail on macOS because it currently doesnt fetch/build libpng and libharu, but otherwise it seems to work well on Windows, and various flavors of Linux.

```bash
# From host OS terminal - grab InterSpec source code
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec
cd InterSpec

# Start Fedora Docker container
docker run --rm -it -v `pwd`:/interspec fedora:35 sh

# Install the packages we need
dnf install -y make automake gcc gcc-c++ kernel-devel cmake git patch npm

# Change to directory we need
cd /interspec/target/electron/

# Install the NPM dependencies we need
npm install -g cmake-js
npm install --save-dev node-addon-api --arch=x64
npm install electron --arch=x64
npm install electron-packager

# Configure and build the node-add-on shared library
#  This step may take half and hour or more to clone into the Boost and Wt github repos, and 
#  then another half hour or hour to build
CMAKE_BUILD_PARALLEL_LEVEL=8 cmake-js --architecture x64 --arch=x64 --CDInterSpec_FETCH_DEPENDENCIES=ON --out=build_linux --target install

# Package the code up
npm run package-linux

# If you are using Docker, you will need to install some additional libraries to run things, although the GUI still probably wont run
#   dnf install nss atk at-spi2-atk cups libdrm gtk3

# Run the executable 
cd release-builds/InterSpec-linux-x64/
./InterSpec
```



## Using Docker to build Electron-based InterSpec package
These are the instructions for building the Electron-based InterSpec package for Linux, using the Python Many Linux docker image, and the `FetchContent` option of the CMake build system to compile boost and Wt.

### Building using a manylinux container
```bash
# From your host OS terminal, run the following commands
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_linux_electron_build
cd InterSpec_linux_electron_build/

# Grab the the oldest (currently) supported manylinux image to your machine
docker pull quay.io/pypa/manylinux2014_x86_64:latest

# Start a shell session within the image, mapping the InterSpec source 
#  directory to /interspec.  We'll also map port 8081 for testing.
docker run --rm -it -v `pwd`:/interspec quay.io/pypa/manylinux2014_x86_64:latest sh

# Get the dependancies we need to build InterSpec
yum update
yum install -y npm

# Make and cd into build directory - note this is in host filesystem incase we  
#  want to come back to things, but dont want to rebuild everything from scratch
cd /interspec/target/electron/

# Install the NPM dependancies we need
npm install uglify-js -g
npm install uglifycss -g
npm install cmake-js -g
npm install --save-dev node-addon-api --arch=x64
npm install electron --arch=x64
npm install electron-packager


# This next command will take like 10 or 20 minutes to clone into the boost and Wt repositories
#  Since we are are building a shared library Electron will load, we will compile shared libraries
#  to link against the static gcc libraries.
#  This command may take a long time to run.
CMAKE_BUILD_PARALLEL_LEVEL=6 cmake-js --architecture x64 --arch=x64 --CDCMAKE_BUILD_TYPE="Release" --CDInterSpec_FETCH_DEPENDENCIES=ON --CDBUILD_AS_LOCAL_SERVER=OFF --CDCMAKE_SHARED_LINKER_FLAGS="-static-libgcc -static-libstdc++" --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="..." --CDUSE_REL_ACT_TOOL=ON --out=build_manylinux_electron --target install

npm run package-manylinux
cp ../../NOTICE.html ./release-builds/InterSpec-linux-x64/

# The "InterSpec-linux-x64" directory contains the InterSpec package.
#  To transfer to another machine just zip this directory up and copy over.
# Or to run the program, run:
./release-builds/InterSpec-linux-x64/InterSpec
```


## Linux Considerations
The `InterSpec` module is really a shared library that node.js loads, therefore you need the `-fPIC` C/C++ compiler flag enabled not just for the `InterSpec` code, but for all of the static libraries you link it against, including boost, Wt, and zlib - which isnt the default when compiling static libraries for any of them, so when building them you may want to add `-fPIC -std=c++11` flags to the flags.


## Future Work
- In the future the build process may be improved to do the final packaging through CMake.  
- launch_options.json is intended to allow users to customize behaviors (like allow multiple windows, or not), but is not actually implemented.
