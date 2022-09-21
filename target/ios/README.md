# Building for iOS
To build for iOS we will take advantage of the CMake `FetchContent` style build of InterSpec to build boost and Wt for us.
What we'll do is build the InterSpec library, and then the InterSpec/InterSpec.xcodeproj Xcode project will build the app and link to the InterSpec library.  

```bash
cd InterSpec

# First we'll build for actual devices
mkdir build-iphoneos; cd build-iphoneos

unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../3rd_party/ios-cmake/ios.toolchain.cmake -DPLATFORM=OS64 -DIOS_DEPLOYMENT_TARGET=11.0 -DDEPLOYMENT_TARGET=11.0 -DENABLE_BITCODE_INT=OFF -DENABLE_BITCODE=OFF ..
# Note that ceres-solver-src/CMakeLists.txt by default defines -mllvm, which is not compatible with bitcode, unless you set ENABLE_BITCODE=ON; I'm not sure if -mllvm is actually going through to the final executable abd being used (so if we do use it, maybe define it upfront in our CMakeLists.txt)

cmake --build . --config Release -j8
# Note: working on doing a local build folder install... see CMAKE_INSTALL_PREFIX
cmake --build . --config Release -j8 --target install
cp -R _deps/wt-src/resources ./WtsRsrcs


# Now (optionally) build for simulator
mkdir build-iphonesimulator; cd build-iphonesimulator
unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../3rd_party/ios-cmake/ios.toolchain.cmake -DPLATFORM=SIMULATOR64 -DIOS_DEPLOYMENT_TARGET=11.0 -DDEPLOYMENT_TARGET=11.0 -DENABLE_BITCODE_INT=OFF -DENABLE_BITCODE=OFF ..
cmake --build . --config Release -j8
cp -R _deps/wt-src/resources ./WtsRsrcs
```

Now open InterSpec/InterSpec.xcodeproj and build the app.  
Note that the Xcode project has a "Targets" -> "InterSpec" -> "Build Phases" -> "Run Script" defined that triggers a `cmake --build` command for the `build-iphoneos` or `build-iphonesimulator` directory, so if the InterSpec library code was updated, the library will be re-built when you hit the build button in Xcode for the app project.


TODO: need to modify Wts CMakeLists.txt to not set the install directory
TODO: try adding a "-G Xcode" to the above, and then try integrating into iOS project, like at https://blog.tomtasche.at/2019/05/how-to-include-cmake-project-in-xcode.html - maybe cleaner that what I did now.


# [depreciated] Manually building Boost and Wt for iOS
The instructions/information below here are depreciated, and will be removed in the future; please use the `FetchContent` method of building boost/Wt.

## Building Boost for iOS
```bash
cd InterSpec/target/ios
curl -L https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz --output boost_1_65_1.tar.gz
tar -xzvf boost_1_65_1.tar.bz2
cd boost_1_65_1
```

It looks like bootstrap no longer works on macOS 10.15. So you will need to do the bootstrap on an older macOS version, or something, to create the `b2` executable, and copy it here.

```bash
# Modify some build flags to so we use the right atomics and such:
patch -u tools/build/example/user-config.jam -i ../patch/boost_1_65_1/user-config.jam.patch
patch -u tools/build/src/tools/darwin.jam -i ../../patch/boost_1_65_1/darwin.jam.patch
```

Then ran the following commands to build both for the actual devices, as well the x86 simulator.  Note if you are using a non-standard iPhone SDK, you will need to modify the path in the below commands.
```bash
# build for device
unset MACOSX_DEPLOYMENT_TARGET

./b2 --build-dir=iphone-build -sBOOST_BUILD_USER_CONFIG=tools/build/example/user-config.jam --stagedir=iphone-build/stage --prefix=../build/ios/prefix abi=aapcs binary-format=mach-o address-model=64 toolset=darwin variant=release architecture=arm optimization=space target-os=iphone macosx-version=iphone-14.4 define=_LITTLE_ENDIAN link=static threading=multi --with-random --with-regex --with-random --with-chrono --with-serialization --with-thread --with-signals --with-filesystem --with-system --with-date_time --with-exception --with-program_options --with-timer --with-atomic --with-iostreams cxxflags="-std=c++11 -stdlib=libc++ -DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -mios-version-min=9.0 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS.sdk" linkflags="-stdlib=libc++" -j12 install

# build for simulator
./b2 --build-dir=iphonesim-build -sBOOST_BUILD_USER_CONFIG=tools/build/example/user-config.jam --stagedir=iphonesim-build/stage --prefix=../build/ios/prefix binary-format=mach-o address-model=64 toolset=darwin-14.4~iphonesim architecture=x86 target-os=iphone variant=release optimization=space target-os=iphone macosx-version=iphonesim-14.4 link=static threading=multi --with-random --with-regex --with-random --with-chrono --with-serialization --with-thread --with-signals --with-filesystem --with-system --with-date_time --with-exception --with-program_options --with-timer --with-atomic --with-iostreams cxxflags="-std=c++11 -stdlib=libc++ -mios-version-min=9.0 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/iPhoneSimulator.platform/Developer/SDKs/iPhoneSimulator.sdk" linkflags="-stdlib=libc++" -j12 stage

mkdir -p ../build/libs/boost/lib/arm64
mkdir -p ../build/libs/boost/lib/x86_64

for NAME in atomic chrono serialization date_time exception filesystem iostreams program_options random regex signals system thread timer; do
	xcrun --sdk iphoneos lipo ../build/ios/prefix/lib/libboost_$NAME.a -thin arm64 -o ../build/libs/boost/lib/arm64/libboost_$NAME.a
done

# For x86, we seem to fail, so well just copy
cp iphonesim-build/stage/lib/libboost_* ../build/libs/boost/lib/x86_64/
# xcrun --sdk iphoneos lipo iphonesim-build/stage/lib/libboost_atomic.a -thin x86_64 -o ../build/libs/boost/lib/x86_64/libboost_atomic.a


# Create fat libs
mkdir ../build/libs/boost/lib/allarch
for NAME in atomic chrono serialization date_time exception filesystem iostreams program_options random regex signals system thread timer; do
	xcrun -sdk iphoneos lipo -create "../build/ios/prefix/lib/libboost_$NAME.a" "iphonesim-build/stage/lib/libboost_$NAME.a" -output "../build/libs/boost/lib/allarch/libboost_$NAME.a"
done

rm ../build/ios/prefix/lib/libboost_*

# Create uber lib - this isnt really necassary - we dont use it, but it could be nice to just link to a single boost library
cd ../build/libs/boost/lib/arm64
mkdir obj; cd obj
for NAME in atomic chrono serialization date_time exception filesystem iostreams program_options random regex signals system thread timer; do
	ar -x ../libboost_$NAME.a
done
cd ../../x86_64; mkdir obj; cd obj
for NAME in atomic chrono serialization date_time exception filesystem iostreams program_options random regex signals system thread timer; do
	ar -x ../libboost_$NAME.a
done

cd ../../arm64
xcrun -sdk iphoneos ar crus libboost.a obj/*.o
cd ../x86_64
xcrun -sdk iphoneos ar crus libboost.a obj/*.o
cd ..
lipo -c arm64/libboost.a x86_64/libboost.a -output ./boost.a

cd InterSpec/target/ios/build/ios/prefix/lib
cp ../../../libs/boost/lib/allarch/libboost_* .

cd InterSpec/target/ios
cp -r ./build/ios/prefix .
```


## Building Wt for iOS
```bash
cd InterSpec/target/ios/
mkdir build_wt; cd build_wt
bash
unset MACOSX_DEPLOYMENT_TARGET

curl -L https://github.com/emweb/wt/archive/3.3.4.tar.gz --output wt-3.3.4.tar.gz
tar -xzvf wt-3.3.4.tar.gz
patch -u wt-3.3.4/CMakeLists.txt -i ../../patches/wt/3.3.4/CMakeLists.txt.patch
patch -u wt-3.3.4/src/Wt/Render/CssParser.C -i ../../patches/wt/3.3.4/CssParser.C.patch
patch -u wt-3.3.4/src/http/RequestParser.C -i ../../patches/wt/3.3.4/RequestParser.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr.C -i ../../patches/wt/3.3.4/ptr.C.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr -i ../../patches/wt/3.3.4/ptr.patch
patch -u wt-3.3.4/src/Wt/Dbo/ptr_impl.h -i ../../patches/wt/3.3.4/ptr_impl.h.patch
patch -u wt-3.3.4/src/Wt/Dbo/Query.C -i ../../patches/wt/3.3.4/Query.C.patch
patch -u wt-3.3.4/src/Wt/WDllDefs.h -i ../../patches/wt/3.3.4/WDllDefs.h.patch

# Modify the build script for iOS included with Wt to work with current Xcode (tested on Xcode 12.4)
cp ./wt-3.3.4/target/osx/build-framework.sh .
patch -u build-framework.sh ../patch
patch -u build-framework.sh ../patch/wt/3.3.4/build-framework.sh.patch

# Started from patched Wt 3.3.4
./build-framework.sh

# For some reason I didnt bother to debug, the commands fail, so just then manually ran as:
xcrun -sdk iphoneos lipo -create ./build-arm64/src/libwttest.a ./build-x86-64/src/libwttest.a -output ../prefix/lib/libwttest.a
xcrun -sdk iphoneos lipo -create ./build-arm64/src/Wt/Dbo/backend/libwtdbosqlite3.a ./build-x86-64/src/Wt/Dbo/backend/libwtdbosqlite3.a -output ../prefix/lib/libwtdbosqlite3.a
xcrun -sdk iphoneos lipo -create ./build-arm64/src/http/libwthttp.a ./build-x86-64/src/http/libwthttp.a -output ../prefix/lib/libwthttp.a
xcrun -sdk iphoneos lipo -create ./build-arm64/src/Wt/Dbo/libwtdbo.a ./build-x86-64/src/Wt/Dbo/libwtdbo.a -output ../prefix/lib/libwtdbo.a
xcrun -sdk iphoneos lipo -create ./build-arm64/src/libwt.a ./build-x86-64/src/libwt.a -output ../prefix/lib/libwt.a
```

And you should now have a directory `InterSpec/target/ios/prefix` that contains all the boost and Wt includes and static libraries you need.



## Setting up the InterSpec libraries for the iOS InterSpec project to link to
```bash
cd InterSpec
mkdir build-iphoneos; cd build-iphoneos
IPHONEOS_DEPLOYMENT_TARGET=10.2 cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/iOSToolchain.cmake -DIOS=ON -DCMAKE_PREFIX_PATH=../target/ios/prefix ..

cd ..
mkdir build-iphonesimulator; cd build-iphonesimulator
IPHONEOS_DEPLOYMENT_TARGET=10.2 cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cmake/iOSToolchain.cmake -DIOS=ON -DCMAKE_PREFIX_PATH=../target/ios/prefix -DIOS_PLATFORM=SIMULATOR ..
```

iOS projects dont seem to like including the Wt web-resources forlder as 'resources', so we need a different name, so do:
```bash
cd InterSpec/target/ios/InterSpec
cp -r ./target/ios/prefix/share/Wt/resources ./XcodeFiles/WtsRsrcs
```

Now you can open the `InterSpec/target/ios/InterSpec/InterSpec.xcodeproj` and when you compile, it will run the `make` command for either the `build-iphoneos` or `build-iphonesimulator` target, depending on your build target.





