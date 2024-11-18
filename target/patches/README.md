# Building dependencies for InterSpec
When building InterSpec, you can either have CMake download and build all dependencies through the `InterSpec_FETCH_DEPENDENCIES` CMake option (e.g., `cmake -DInterSpec_FETCH_DEPENDENCIES=ON ..`), or you can build all dependencies and place them into a prefix location.  Having CMake fetch and build the dependencies is required for iOS and Android, while creating a prefix location with the dependencies pre-built is highly recommended for desktop OSes, especially when doing heavy development, as it is slightly faster to re-compile the InterSpec code against since the build system has much less files to track and check if they have been changed.  

Using package managers like `brew`, `pacman`, or `nuget` is very unlikely to work - and not suggested.  


Jump to:
- [macOS](#Building-dependencies-on-macOS-Catalina)
- [Windows](#Building-dependencies-on-Windows-10)
- [Linux](#Building-dependencies-on-Linux) 

For iOS and Android, only using `InterSpec_FETCH_DEPENDENCIES=ON` is supported.


# Building dependencies on macOS Catalina
First, lets set the directory we will install all the pre-requisites to:
```bash
export MACOSX_DEPLOYMENT_TARGET=10.13 # MacOS High Sierra (2017, supported through 2020)
export MY_WT_PREFIX=/path/to/install/prefix/to/macOS_wt3.7.1_prefix
export PATCH_DIR=/path/to/InterSpec/target/patches/
```

## Build boost 1.78
For instructions on building boost 1.65.1 for Wt 3.3.4, see versions of this file before Dec 26th 2021.

```bash
curl -L https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download --output boost_1_78_0.zip
unzip boost_1_78_0.zip
cd boost_1_78_0

# build the b2 executable
./bootstrap.sh cxxflags="-arch x86_64 -arch arm64 -mmacosx-version-min=10.13" cflags="-arch x86_64 -arch arm64 -mmacosx-version-min=10.13" linkflags="-arch x86_64 -arch arm64 -mmacosx-version-min=10.13" --prefix=${MY_WT_PREFIX}

# build and stage boost for arm64
./b2 toolset=clang-darwin target-os=darwin architecture=arm abi=aapcs cxxflags="-stdlib=libc++ -arch arm64 -std=c++14 -mmacosx-version-min=10.13" cflags="-arch arm64  -mmacosx-version-min=10.13" linkflags="-stdlib=libc++ -arch arm64 -std=c++14 -mmacosx-version-min=10.13" link=static variant=release threading=multi --build-dir=macOS_arm64_build --prefix=${MY_WT_PREFIX} -a stage

# copy arm libraries to a separate directory
mkdir -p arm64 && cp stage/lib/libboost_* arm64/

# build boost for x86_64 and install it (we'll copy over the libraries later)
./b2 toolset=clang-darwin target-os=darwin architecture=x86 cxxflags="-stdlib=libc++ -arch x86_64 -std=c++14 -mmacosx-version-min=10.13" cflags="-arch x86_64 -mmacosx-version-min=10.13" linkflags="-stdlib=libc++ -arch x86_64 -std=c++14 -mmacosx-version-min=10.13" abi=sysv binary-format=mach-o link=static variant=release threading=multi --build-dir=macOS_x64_build --prefix=${MY_WT_PREFIX} -a install

# move x86 libraries to a seperate directory
mkdir x86_64 && mv ${MY_WT_PREFIX}/lib/libboost_* x86_64/

# Now lipo libraries together
mkdir universal
for dylib in arm64/*; do 
  lipo -create -arch arm64 $dylib -arch x86_64 x86_64/$(basename $dylib) -output universal/$(basename $dylib); 
done

# And finally copy the fat libraries to the prefix location
cp universal/libboost_* ${MY_WT_PREFIX}/lib/

cd ..
```

## Build libpng
libpng is only necessary if you are going to build the macOS packaged app, with the Quick Look utility that provides previews of spectrum files in the Finder and various places throughout the OS.
We will build once for arm64, and then x86_64 and lipo the libraries together to create a fat library.
```bash
curl -L http://prdownloads.sourceforge.net/libpng/libpng-1.6.37.tar.gz --output libpng-1.6.37.tar.gz
tar -xzvf libpng-1.6.37.tar.gz
cd libpng-1.6.37
mkdir build
cd build

# First build for arm64 (-DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" doesnt seem to work)
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="arm64" -DPNG_ARM_NEON=on ..
make -j10 install
rm -rf ./*

# Then build for x86_64
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="x86_64" ..
make -j10

# And now lipo the libraries together
mkdir universal
lipo -create -arch x86_64 libpng16.a -arch arm64 ${MY_WT_PREFIX}/lib/libpng16.a -output universal/libpng16.a
cp universal/libpng16.a ${MY_WT_PREFIX}/lib/
cd ../..
```

## Build libharu
libharu is only necessary if you are going to build the macOS packaged app, with the Quick Look utility.
```bash
curl -L https://github.com/libharu/libharu/archive/RELEASE_2_3_0.zip --output libharu_2.3.0.zip
unzip libharu_2.3.0.zip
cd libharu-RELEASE_2_3_0
mkdir build
cd build

# CMake will take care of building for x86_64 and arm64 at the same time
cmake -DLIBHPDF_SHARED=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DPNG_LIBRARY_RELEASE= -DPNG_LIBRARY_RELEASE=${MY_WT_PREFIX}/lib/libpng.a -DPNG_PNG_INCLUDE_DIR=${MY_WT_PREFIX}/include -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
make -j10 install
cd ../..
```

## Build Wt 3.7.1
InterSpec code will only compile with exactly version 3.3.4 or 3.7.1.  
For 3.3.4, see version of this file from before Dec 26th 2021 for pathc and build instructions.

```bash
curl -L https://github.com/emweb/wt/archive/3.7.1.tar.gz --output wt-3.7.1.tar.gz
tar -xzvf wt-3.7.1.tar.gz

cd wt-3.7.1
git apply --reject --ignore-space-change --ignore-whitespace "${PATCH_DIR}/wt/3.7.1/NormalBuild/wt_3.7.1_git.patch"
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DBoost_INCLUDE_DIR=${MY_WT_PREFIX}/include -DBOOST_PREFIX=${MY_WT_PREFIX} -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DHARU_PREFIX=${MY_WT_PREFIX} -DHARU_LIB=${MY_WT_PREFIX}/lib/libhpdfs.a -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++14" -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=${MY_WT_PREFIX}/etc/wt -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
make -j10 install
```

## Install Eigen
```bash
# Build Eigen, which is required by ceres-solver, and used a few other places
# in InterSpec if its available
curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz --output eigen-3.4.0.tar.gz
tar -xzvf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
cmake --build . --config Release --target install
cd ../..
```

## Build Ceres-Solver
```bash
# Build ceres-solver; this is the optimizer used for the relative efficiency
# tool, and a small amount of the peak fitting.
curl -L http://ceres-solver.org/ceres-solver-2.1.0.tar.gz --output ceres-solver-2.1.0.tar.gz

# Double check the hash of the tar file
ceres_hash=`openssl SHA512 ceres-solver-2.1.0.tar.gz | cut -d' ' -f2`
wanted_hash='db777c4a4b85a9972628c714413fe9ca1a92cc6c8b14d80401842a02c4d7e7df4e9d236775e231e8337e47f21e47cd82c1a0fbcfeffbd9b654e46a66262996f2'
if [ "$ceres_hash" = "$wanted_hash" ]; then
  echo "Ceres hash validated."
else
    echo "Ceres hash did not validate - you should not procede."
fi

tar -xzvf ceres-solver-2.1.0.tar.gz
cd ceres-solver-2.1.0
mkdir build_macos
cd build_macos

cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
cmake --build . --config Release --target install -j 16
cd ../..
```

## Build wxWidgets (very optional)
wxWidgets for macOS is only for development purposes
```bash
mkdir wxWidgets-3.2.1
cd wxWidgets-3.2.1
curl -L https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.1/wxWidgets-3.2.1.zip --output wxWidgets-3.2.1.zip
tar -xzvf wxWidgets-3.2.1.zip


mkdir build_macOS && cd build_macOS
# TODO: we could/should turn off a lot of wxWidgets components
#  Untested of both arm and x64 architectures are build
cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DwxBUILD_SHARED=OFF -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
cmake --build . --config Release --target install -j10
```

## Build InterSpec to run on localhost
This is most useful for development.
```bash
cd ${PATCH_DIR}/../..
mkdir build_xcode
cd build_xcode
cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" -DTRY_TO_STATIC_LINK=ON -DUSE_SPECRUM_FILE_QUERY_WIDGET=ON -DUSE_TERMINAL_WIDGET=ON -G Xcode ..

# Debug builds reference '/external_libs/SpecUtils/d3_resources/d3.min.js'
# and similar for SpectrumChartD3.js/.css, so we can edit the files under
# version control (rather than copies), so to do this we need to create a 
# symbolic link
cd external_libs/SpecUtils
ln -s ../../../external_libs/SpecUtils/d3_resources .
cd ../..

# Make the directory user settings and such will be saved to
mkdir user_data

# Open the generated project in Xcode
open InterSpec.xcodeproj

# Then after compiling, point your browser to http://localhost:8080

## Build InterSpec as macOS app
Note, not currently building as fat binary (have not tried)
```bash
cd ${PATCH_DIR}/../osx/
mkdir build_macOS_app
cd build_macOS_app
cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -G Xcode ..
open InterSpecApp.xcodeproj
```


# Building dependencies on Windows 10 
For Windows, it is recommended to use the [dep_build_msvc2022.bat](dep_build_msvc2022.bat) script to build the prerequisites. However, instructions for manually building the prerequisites are also provided below, which are also compatible with Visual Studio 2019, or just the command compiler installed.

## Building dependencies using dep_build_msvc2022.bat
You will need to install Visual Studio 2022.

If you plan on building debug builds of InterSpec, edit the `builddebug` variable at the top of [dep_build_msvc2022.bat](dep_build_msvc2022.bat) to be equal to something - this will cause both Release and Debug builds of the prerequisite libraries to be built and installed.

The [dep_build_msvc2022.bat](dep_build_msvc2022.bat) script takes two argument, the first is the directory to use to compile the libraries, and the second is the path to install the libraries to.  This script is meant to be run from the "x64 Native Tools Command Prompt for VS 2022" terminal, that is in the "Visual Studio 2022" folder in the Start menu.



From the "x64 Native Tools Command Prompt for VS 2022" terminal, run:
```batch
# Checkout the InterSpec code
cd C:\Path\To\Put\InterSpec\Code
git clone --recursive https://github.com/sandialabs/interspec/

# Change CWD to one with compile script
cd interspec\target\patches

# Optionally set environment variables for the proxy (only if you are behind one)
set http_proxy=http://your-proxy:port
set https_proxy=http://your-proxy:port

# Set where we want all the libraries installed to
set MY_PREFIX=C:\install\msvc2022\x64\wt_3.7.1_prefix

# Run script to compile all the libraries
dep_build_msvc2022.bat C:\temp\interspec_prereq_build %MY_PREFIX%
```

You can now compile InterSpec, either by using CMake to explicitly create a VS project, or by opening the InteSpec directory in Visual Studio.  In either case, placing a file, `CMakeSettings.json`, with contents similar to the following, into the interspec code directory will help configure CMake.


Example `CMakeSettings.json` contents:
```
{
  "configurations": [
    {
      "name": "x64-Debug",
      "generator": "Ninja",
      "configurationType": "Debug",
      "inheritEnvironments": [ "msvc_x64_x64" ],
      "buildRoot": "${projectDir}\\out\\build\\${name}",
      "installRoot": "${projectDir}\\out\\install\\${name}",
      "cmakeCommandArgs": "-DCMAKE_PREFIX_PATH=\"C:\\install\\msvc2022\\x64\\wt_3.7.1_prefix\"",
      "buildCommandArgs": "",
      "ctestCommandArgs": "",
      "variables": [
        {
          "name": "BUILD_PARALLEL",
          "value": "/MP",
          "type": "STRING"
        }
      ]
    },
    {
      "name": "x64-Release",
      "configurationType": "Release",
      ...
    }
  ]
}
```
You will need to adjust the `prefix` path to suite where you chose to install the libraries to.

In MSVC, you should be able to then compile and run the code.  If you opened from the base directory of the InterSpec code, you can access the GUI at http://localhost:8080, by default.  You could instead open the `target\wxWidgets` directory, in which case lets you build the app form of the application.


## Manually building dependencies for Windows

The Windows release build of InterSpec uses either [wxWidgets](https://www.wxwidgets.org) (default) or [Electron](https://electronjs.org/) to render the app contents.


For Electron (https://electronjs.org/) to render the application; the InterSpec build scripts require using Node.js and NPM to build things; these can downloaded from https://nodejs.org/en/download/.
During the installation of Node.js, you can choose to install the [Chocolatey](https://chocolatey.org/) package manager; this will install the necessary Visual Studio command line compiler tools to install things (i.e., you dont need a full install of MSVC).
You also need to install CMake, most easily from https://cmake.org/download/.

The InterSpec build files are setup to use the MSVC static runtime - it is highly suggested to compile boost, Wt, zlib and optionally Ceres Solver and Eigen, from from source, following these instructions.

If you wish to only build the `wxWidgets` version of the app, you can choose to not use or install NPM or Electron in the below instructions.  Similarly, if you are only targeting `Electron`, you dont need wxWidgets.


From the Visual Studio 2022 "x64 Native Tools Command Prompt":
```bash
cd C:\temp
mkdir build
cd build

curl -L https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download --output boost_1_78_0.tar.gz
tar -xzvf boost_1_78_0.tar.gz
cd boost_1_78_0
bootstrap.bat
```

It shouldnt be necessary, but if you have trouble building boost, you may need to edit project-config.jam to change
```bash
using msvc ;
```
To something like
```bash
using msvc : 14.2 : "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64\cl.exe";
```


Then compile and install:

```bash
set MY_PREFIX=C:\install\msvc2022\x64\wt_3.7.1_prefix

.\b2.exe runtime-link=static link=static threading=multi variant=release address-model=64 architecture=x86 --prefix=%MY_PREFIX% --build-dir=win_build -j8 install
```


## Build zlib
`zlib` is used by InterSpec to allow you read-in  (usually by dragging and dropping) zip files that contain spectrum files.

```bash
git clone https://github.com/madler/zlib.git zlib-1.2.13
# Checkout version 1.2.13
git checkout 04f42ceca40f73e2978b50e93806c2a18c1281fc

mkdir build

cmake -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" ..
cmake --build . --config Debug --target install
cmake --build . --config Release --target install

# Optional: Get rid of zlib dll's, to avoid CMake possible always preffering them
del "%MY_PREFIX%\bin\zlib.dll"
del "%MY_PREFIX%\bin\zlibd.dll"
del "%MY_PREFIX%\lib\zlib.lib"
del "%MY_PREFIX%\lib\zlibd.lib"
```

To build Wt, you must patch the Wt source code, and then build it
```bash
curl -L https://github.com/emweb/wt/archive/3.7.1.tar.gz --output wt-3.7.1.tar.gz
tar -xzvf wt-3.7.1.tar.gz
cd wt-3.7.1

set WT_PATCH_FILE="%PATCH_DIR%\wt\3.7.1\NormalBuild\wt_3.7.1_git.patch"
git apply --ignore-space-change --ignore-whitespace %WT_PATCH_FILE%


mkdir build_msvc2022
cd build_msvc2022

cmake -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DBoost_INCLUDE_DIR=%MY_PREFIX%/include -DBOOST_PREFIX=%MY_PREFIX% -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DINSTALL_FINDWT_CMAKE_FILE=OFF -DCONFIGURATION=data/config/wt_config_electron.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=%MY_PREFIX%/etc/wt -DBoost_USE_STATIC_RUNTIME=ON ..

cmake --build . --config Release --target install
cmake --build . --config Debug --target install
```

If you plan to package InterSpec as an Electron application (e.g., normal desktop app), see the instructions in [patches](/target/electron/) for building the InterSpec code and packaging the application.


## Build Eigen and Ceres
If you plan to compile InterSpec with the relative activity tools enabled (CMake option `USE_REL_ACT_TOOL` - default is on), then you will also need to compile/install http://ceres-solver.org/installation.html#section-windows

```bash
# Build Eigen, which is required by ceres-solver, and used a few other places
# in InterSpec if its avaiable
curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz --output eigen-3.4.0.tar.gz
tar -xzvf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
cmake --build . --config Release --target install


# Build ceres-solver; this is the optimizer used for the relative efficiency
# tool, and a small amount of the peak fitting.
curl -L http://ceres-solver.org/ceres-solver-2.1.0.tar.gz --output ceres-solver-2.1.0.tar.gz
tar -xzvf ceres-solver-2.1.0.tar.gz
cd ceres-solver-2.1.0
mkdir build_msvc
cd build_msvc

cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDebug -DMSVC_USE_STATIC_CRT=ON ..
cmake --build . --config Debug --target install -j 16

cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded -DMSVC_USE_STATIC_CRT=ON ..
cmake --build . --config Release --target install -j 16
```

If you wish to build the wxWidgets target of InterSpec, then you will need wxWidgets:
```bash
mkdir wxWidgets-3.2.1
cd wxWidgets-3.2.1
curl -L https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.1/wxWidgets-3.2.1.zip --output wxWidgets-3.2.1.zip
tar -xzvf wxWidgets-3.2.1.zip

mkdir build_msvc
# TODO: we could/should turn off a lot of wxWidgets components
cmake -DCMAKE_PREFIX_PATH=%MY_PREFIX% -DCMAKE_INSTALL_PREFIX=%MY_PREFIX% -DwxUSE_WEBVIEW_EDGE=ON -DwxUSE_WEBVIEW_EDGE_STATIC=ON -DwxBUILD_USE_STATIC_RUNTIME=ON -DwxBUILD_SHARED=OFF ..
cmake --build . --config Debug --target install -j 16
cmake --build . --config Release --target install -j 16

# cmake doesnt seem to install WebView2Loader - so we'll manually copy its stuff
copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\x64\WebView2LoaderStatic.lib "%MY_PREFIX%\lib"
copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\x64\WebView2Guid.lib "%MY_PREFIX%\lib"
copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\include\WebView2.h "%MY_PREFIX%\include"
copy packages\Microsoft.Web.WebView2.1.0.705.50\build\native\include\WebView2EnvironmentOptions.h "%MY_PREFIX%\include"
```




# Building dependencies on Linux 
The commands to build both the prerequisites and InterSpec on Ubuntu 18.04 are below.

These commands build boost, Wt, and zlib and install them to a non-system directory prefix so they wont interfere with any other builds, or use the system package managers packages that may not be the exact version of Wt needed.


```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential cmake cmake-curses-gui git

# First, lets set the directory we will install all the pre-requisites to:
export MY_WT_PREFIX=/home/wcjohns/install/wt3.7.1_prefix
export PATCH_DIR=/mnt/hgfs/wcjohns/rad_ana/InterSpec_master/target/patches

# Move to a temporary location to build boost/wt
cd /tmp
mkdir build-wt
cd build-wt

# Download and build boost 1.78
curl -L https://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip/download --output boost_1_78_0.zip
unzip boost_1_78_0.zip
cd boost_1_78_0

# build the b2 executable
./bootstrap.sh --prefix=${MY_WT_PREFIX}

# build and install boost
./b2 cxxflags="-std=c++14 -fPIC" cflags="-fPIC" linkflags="-std=c++14" link=static variant=release threading=multi --build-dir=build --prefix=${MY_WT_PREFIX} -a install
cd ..


## Download and build Wt 3.7.1
curl -L https://github.com/emweb/wt/archive/3.7.1.tar.gz --output wt-3.7.1.tar.gz
tar -xzvf wt-3.7.1.tar.gz

cd wt-3.7.1
patch -u src/Wt/Render/CssParser.C -i ${PATCH_DIR}/wt/3.7.1/CssParser.C.patch
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DBOOST_PREFIX=${MY_WT_PREFIX} -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DINSTALL_FINDWT_CMAKE_FILE=OFF -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++14" -DCONFIGURATION=data/config/wt_config_electron.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=${MY_WT_PREFIX}/etc/wt -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
make -j10 install
cd ../..

## Download and build zlib (for supporting zipped spectrum files)
curl -L  https://www.zlib.net/zlib-1.2.12.tar.gz --output zlib-1.2.12.tar.gz
tar -xzvf zlib-1.2.12.tar.gz
cd zlib-1.2.12

mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
cmake --build . --config Release --target install
cd ../..
```

## Build Eigen and Ceres
If you plan to compile InterSpec with the relative activity tools enabled (CMake option `USE_REL_ACT_TOOL` - default is on), then you will also need to compile/install http://ceres-solver.org/installation.html#linux

```bash
# Build Eigen, which is required by ceres-solver, and used a few other places
# in InterSpec if its avaiable.
curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz --output eigen-3.4.0.tar.gz
tar -xzvf eigen-3.4.0.tar.gz
cd eigen-3.4.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_BUILD_TYPE=Release -DEIGEN_MPL2_ONLY=1 -DEIGEN_BUILD_SHARED_LIBS=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_TESTING=OFF ..
cmake --build . --config Release --target install


# Build ceres-solver; this is the optimizer used for the relative efficiency
# tool, and a small amount of the peak fitting.
curl -L http://ceres-solver.org/ceres-solver-2.1.0.tar.gz --output ceres-solver-2.1.0.tar.gz
tar -xzvf ceres-solver-2.1.0.tar.gz
cd ceres-solver-2.1.0
mkdir build_msvc
cd build_msvc

cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DMINIGLOG=ON -DGFLAGS=OFF -DCXSPARSE=OFF -DACCELERATESPARSE=OFF -DCUDA=OFF -DEXPORT_BUILD_DIR=ON -DBUILD_TESTING=ON -DBUILD_EXAMPLES=OFF -DPROVIDE_UNINSTALL_TARGET=OFF -DBUILD_SHARED_LIBS=OFF ..
cmake --build . --config Release --target install -j 16
```



You can now build a local web server:
```bash
cd ~/development
git clone --recursive https://github.com/sandialabs/interspec/
cd interspec
mkdir build; cd build

cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} ..
make -j8

# If you compiled it in debug mode, you may need symlink the D3.js based resources to where the web-requests expect to find them, e.g.
ln -s ../external_libs/SpecUtils/d3_resources ./external_libs/
# (this is so you can edit d3_resources/SpectrumChartD3.js and see changes without doing a file copy)


# And if all went well you can run the executable using the command:
./bin/InterSpec.exe --docroot . --http-address 127.0.0.1 --http-port 8080 -c ./data/config/wt_config_localweb.xml

# And then point your browser to http://localhost:8080 to access the application
```

Or if you want the desktop version of the application; see InterSpec/target/electron/README.md for full build instructions, but a summary is:
```bash
cd /path/to/InterSpec/target/electron

sudo apt-get install nodejs npm
npm install -g cmake-js
npm install --save-dev node-addon-api
npm install electron
npm install electron-packager

cmake-js --CDCMAKE_PREFIX_PATH=${MY_WT_PREFIX} --CDCMAKE_BUILD_TYPE="Release" --out="build_linux" --target install

# To run InterSpec:
electron build_linux/app

# Or to create a self-contained package of the app
npm run package-linux

#InterSpec is now in release-builds/InterSpec-linux-x64

# (note: you may need to upgrade to a newer version of npm and/or nodejs to
#        use latest versions of electron)
```




# Process for patching files.

## To make a patch
```bash
cp someFile.txt someFile.txt.orig
```
Then edit someFile.txt, then

```bash
diff -Naur someFile.txt.orig  someFile.txt > someFile.txt.patch
```

## To apply a patch
```bash
patch -u someFile.txt -i someFile.txt.patch
```

