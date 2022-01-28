Process for patching files.

# To make a patch
```bash
cp someFile.txt someFile.txt.orig
```
Then edit someFile.txt, then

```bash
diff -Naur someFile.txt.orig  someFile.txt > someFile.txt.patch
```

# To apply a patch
```bash
patch -u someFile.txt -i someFile.txt.patch
```

# Command to apply patches for Wt 3.3.4
```bash
patch -u wt-3.3.4/CMakeLists.txt -i /path/to/InterSpec/target/patches/wt/3.3.4/CMakeLists.txt.patch
patch -u wt-3.3.4/src/Wt/Render/CssParser.C -i /path/to/InterSpec/target/patches/wt/3.3.4/CssParser.C.patch
patch -u wt-3.3.4/src/http/RequestParser.C -i /path/to/InterSpec/target/patches/wt/3.3.4/RequestParser.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr.C -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr_impl.h -i /path/to/InterSpec/target/patches/wt/3.3.4/ptr_impl.h.patch
patch -u wt-3.3.4/src/Wt/Dbo/Query.C -i /path/to/InterSpec/target/patches/wt/3.3.4/Query.C.patch
patch -u wt-3.3.4/src/Wt/WDllDefs.h -i /path/to/InterSpec/target/patches/wt/3.3.4/WDllDefs.h.patch
```


# Building dependencies on macOS Catalina
First, lets set the directory we will install all the pre-requisites to:
```bash
export MACOSX_DEPLOYMENT_TARGET=10.12
export MY_WT_PREFIX=/path/to/install/prefix/to/macOS_wt3.3.4_prefix
export PATCH_DIR=/path/to/InterSpec/target/patches/
```

## Build boost 1.65.1
```bash
curl -L https://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.zip/download --output boost_1_65_1.tar.gz
tar -xzvf boost_1_65_1.tar.gz
cd boost_1_65_1

# Build the "b2" executable, but for newer versions of clang we need to add the
#  "-Wno-error-implicit-function-declaration" compile flag, which is done by:
patch -u tools/build/src/engine/build.sh -i ${PATCH_DIR}/boost/1_65_1/macOS/build.sh.patch
patch -u tools/build/src/engine/build.jam -i ${PATCH_DIR}/boost/1_65_1/macOS/build.jam.patch

# Build "b2"
./bootstrap.sh --prefix=${MY_WT_PREFIX}

# Update the C++ compile flags that will be used for building boost, again for the newer clang
patch -u tools/build/src/tools/darwin.jam -i ${PATCH_DIR}/boost/1_65_1/macOS/darwin.jam.patch

# Build and install boost
./b2 toolset=clang cxxflags="-stdlib=libc++ -std=c++11" linkflags="-stdlib=libc++ -std=c++11" --prefix=${MY_WT_PREFIX} link=static variant=release threading=multi --build-dir=macOS_build -j12 install
cd ..
```

## Build libpng
libpng is only necessary if you are going to build the macOS packaged app, with the Quick Look utility that provides previews of spectrum files in the Finder and various places throughout the OS.
```bash
curl -L http://prdownloads.sourceforge.net/libpng/libpng-1.6.37.tar.gz --output libpng-1.6.37.tar.gz
tar -xzvf libpng-1.6.37.tar.gz
cd libpng-1.6.37
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} ..
make -j12 install
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
cmake -DLIBHPDF_SHARED=OFF -DCMAKE_BUILD_TYPE=Release -DPNG_LIBRARY_RELEASE= -DPNG_LIBRARY_RELEASE=${MY_WT_PREFIX}/lib/libpng.a -DPNG_PNG_INCLUDE_DIR=${MY_WT_PREFIX}/include -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} ..
make -j12 install
cd ../..
```

## Build Wt 3.3.4
InterSpec code will only compile with exactly version 3.3.4 of Wt, and it must be patched according to the following instructions.

```bash
curl -L https://github.com/emweb/wt/archive/3.3.4.tar.gz --output wt-3.3.4.tar.gz
tar -xzvf wt-3.3.4.tar.gz

cd wt-3.3.4
patch -u CMakeLists.txt -i ${PATCH_DIR}/wt/3.3.4/CMakeLists.txt.patch
patch -u src/Wt/Render/CssParser.C -i ${PATCH_DIR}/wt/3.3.4/CssParser.C.patch
patch -u src/http/RequestParser.C -i ${PATCH_DIR}/wt/3.3.4/RequestParser.C.patch
patch -u src/Wt/Dbo/ptr.C -i ${PATCH_DIR}/wt/3.3.4/ptr.C.patch
patch -u src/Wt/Dbo/ptr -i ${PATCH_DIR}/wt/3.3.4/ptr.patch
patch -u src/Wt/Dbo/ptr_impl.h -i ${PATCH_DIR}/wt/3.3.4/ptr_impl.h.patch
patch -u src/Wt/Dbo/Query.C -i ${PATCH_DIR}/wt/3.3.4/Query.C.patch
patch -u src/Wt/WDllDefs.h -i ${PATCH_DIR}/wt/3.3.4/WDllDefs.h.patch
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DBoost_INCLUDE_DIR=${MY_WT_PREFIX}/include -DBOOST_PREFIX=${MY_WT_PREFIX} -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DHARU_PREFIX=${MY_WT_PREFIX} -DHARU_LIB=${MY_WT_PREFIX}/lib/libhpdfs.a -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=${MY_WT_PREFIX}/etc/wt ..
make -j16 install
```

## Build InterSpec as macOS app
```bash
cd ${PATCH_DIR}/../..
cmake -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DBUILD_AS_OSX_APP=ON -DTRY_TO_STATIC_LINK=ON -DUSE_SPECRUM_FILE_QUERY_WIDGET=ON -DUSE_TERMINAL_WIDGET=ON -G Xcode ..
open InterSpec.xcodeproj
```


# Building dependencies on Windows 10 
The Windows release build of InterSpec uses Electron (https://electronjs.org/) to render the application; the InterSpec build scripts require using Node.js and NPM to build things; these can downloaded from https://nodejs.org/en/download/.
During the installation of Node.js, you can choose to install the [Chocolatey](https://chocolatey.org/) package manager; this will install the necessary Visual Studio command line compiler tools to install things (i.e., you dont need a full install of MSVC).
You also need to install CMake, most easily from https://cmake.org/download/.

The InterSpec build files are setup to use the MSVC static runtime - it is highly suggested to compile boost, Wt, zlib from from source, following these instructions.


From the Visual Studio 2019 "x64 Native Tools Command Prompt":
```bash
cd C:\temp
mkdir build
cd build
curl -L https://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.zip/download --output boost_1_65_1.tar.gz
tar -xzvf boost_1_65_1.tar.gz
cd boost_1_65_1
bootstrap.bat
```

Edit project-config.jam to change
```bash
using msvc ;
```
To something like
```bash
using msvc : 14.2 : "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64\cl.exe";
```

Then compile and install:

```bash
set MY_BOOST_PREFIX=C:\install\msvc2019\x64\boost_1_65_1

.\b2.exe runtime-link=static link=static threading=multi variant=release address-model=64 architecture=x86 msvcver=msvc-14.2 --prefix=%MY_BOOST_PREFIX% --build-dir=win_build -j8 install
```


Build zlib
```bash
curl -L https://zlib.net/zlib-1.2.11.tar.gz --output zlib-1.2.11.tar.gz
tar -xzvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
mkdir build

set MY_ZLIB_PREFIX=C:\install\msvc2019\x64\zlib
cmake -DCMAKE_INSTALL_PREFIX=%MY_ZLIB_PREFIX% -DCMAKE_BUILD_TYPE=Release ..
# NOTE - you need to manually change all "MD" compiler flags to "MT" - I used the CMake GUI for this
cmake --build . --config Debug --target install
cmake --build . --config Release --target install
```

To build Wt, you must patch the Wt source code as described above, and then you can use the cmake GUI to configure Wt, and Visual Studio to build and install it.
You can use the following command from the command line to initiate the CMake config, but you also need to change all "MD" compiler flags to "MT", which its easiest to do this in the CMake gui.
```bash
set MY_WT_PREFIX=C:\install\msvc2019\x64\wt3.3.4
mkdir build_msvc2019
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=%MY_WT_PREFIX% -DBoost_INCLUDE_DIR=%MY_BOOST_PREFIX%/include -DBOOST_PREFIX=%MY_BOOST_PREFIX% -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=%MY_WT_PREFIX% -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++11" -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=%MY_WT_PREFIX%/etc/wt ..

cmake --build . --config Release --target install

#Note if you run into permission problems while files are being copied to install location, you maye need to change permissions on CMakes "share" directory, at least temporarily
#Note also I had to manually add the following to Wt's CMakeLists.txt to fix an issue when linking InterSpec
# SET( Boost_USE_STATIC_LIBS ON )
# add_definitions( -DBOOST_ALL_NO_LIB)
```

If you plan to package InterSpec as an Electron application (e.g., normal desktop app), see the instructions in [patches](/target/electron/) for building the InterSpec code and packaging the application.