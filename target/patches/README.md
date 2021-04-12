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


# Building dependancies on macOS Catalina
First, lets set the directory we will install all the pre-requesits to:
```bash
export MACOSX_DEPLOYMENT_TARGET=10.12
export MY_WT_PREFIX=/path/to/install/prefix/to/macOS_wt3.3.4_prefix
export PATCH_DIR=/path/to/InterSpec/target/patches/
```

## Build boost 1.65.1
```bash
curl -L https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --output boost_1_65_1.tar.gz
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