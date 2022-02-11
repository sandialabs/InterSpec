


# Building dependencies on macOS Catalina
First, lets set the directory we will install all the pre-requisites to:
```bash
export MACOSX_DEPLOYMENT_TARGET=10.12
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
./bootstrap.sh cxxflags="-arch x86_64 -arch arm64" cflags="-arch x86_64 -arch arm64" linkflags="-arch x86_64 -arch arm64" --prefix=${MY_WT_PREFIX}

# build and stage boost for arm64
./b2 toolset=clang-darwin target-os=darwin architecture=arm abi=aapcs cxxflags="-stdlib=libc++ -arch arm64 -std=c++14" cflags="-arch arm64" linkflags="-stdlib=libc++ -arch arm64 -std=c++14" link=static variant=release threading=multi --build-dir=macOS_arm64_build --prefix=${MY_WT_PREFIX} -a stage

# copy arm libraries to a seperate directory
mkdir -p arm64 && cp stage/lib/libboost_* arm64/

# build boost for x86_64 and install it (we'll copy over the libraries later)
./b2 toolset=clang-darwin target-os=darwin architecture=x86 cxxflags="-stdlib=libc++ -arch x86_64 -std=c++14" cflags="-arch x86_64" linkflags="-stdlib=libc++ -arch x86_64 -std=c++14" abi=sysv binary-format=mach-o link=static variant=release threading=multi --build-dir=macOS_x64_build --prefix=${MY_WT_PREFIX} -a install

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
cmake -DCMAKE_BUILD_TYPE=Release -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="arm64" -DPNG_ARM_NEON=on ..
make -j10 install
rm -rf ./*

# Then build for x86_64
cmake -DCMAKE_BUILD_TYPE=Release -DPNG_SHARED=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="x86_64" ..
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
cmake -DLIBHPDF_SHARED=OFF -DCMAKE_BUILD_TYPE=Release -DPNG_LIBRARY_RELEASE= -DPNG_LIBRARY_RELEASE=${MY_WT_PREFIX}/lib/libpng.a -DPNG_PNG_INCLUDE_DIR=${MY_WT_PREFIX}/include -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
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
patch -u src/Wt/Render/CssParser.C -i ${PATCH_DIR}/wt/3.7.1/CssParser.C.patch
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${MY_WT_PREFIX} -DBoost_INCLUDE_DIR=${MY_WT_PREFIX}/include -DBOOST_PREFIX=${MY_WT_PREFIX} -DSHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${MY_WT_PREFIX} -DHARU_PREFIX=${MY_WT_PREFIX} -DHARU_LIB=${MY_WT_PREFIX}/lib/libhpdfs.a -DENABLE_SSL=OFF -DCONNECTOR_FCGI=OFF -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DENABLE_MYSQL=OFF -DENABLE_POSTGRES=OFF -DENABLE_PANGO=OFF -DINSTALL_FINDWT_CMAKE_FILE=ON -DHTTP_WITH_ZLIB=OFF -DWT_CPP_11_MODE="-std=c++14" -DCONFIGURATION=data/config/wt_config_osx.xml -DWTHTTP_CONFIGURATION=data/config/wthttpd -DCONFIGDIR=${MY_WT_PREFIX}/etc/wt -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
make -j10 install
```

## Build InterSpec as macOS app
Note, not currently building as fat binary (have not tried)
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
