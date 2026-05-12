# Building for iOS
To build for iOS we will take advantage of the CMake `FetchContent` style build of InterSpec to build boost and Wt for us.
What we'll do is build the InterSpec library, and then the InterSpec/InterSpec.xcodeproj Xcode project will build the app and link to the InterSpec library.  

```bash
cd InterSpec/target/ios/InterSpec

# First we'll build for actual devices
mkdir build-iphoneos; cd build-iphoneos

unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../3rd_party/ios-cmake/ios.toolchain.cmake -DPLATFORM=OS64 -DIOS_DEPLOYMENT_TARGET=16.4 -DDEPLOYMENT_TARGET=16.4 -DENABLE_BITCODE_INT=OFF -DENABLE_BITCODE=OFF ..
# Note that ceres-solver-src/CMakeLists.txt by default defines -mllvm, which is not compatible with bitcode, unless you set ENABLE_BITCODE=ON; I'm not sure if -mllvm is actually going through to the final executable abd being used (so if we do use it, maybe define it upfront in our CMakeLists.txt)

cmake --build . --config Release -j8
# Note: working on doing a local build folder install... see CMAKE_INSTALL_PREFIX
cmake --build . --config Release -j8 --target install
cp -R _deps/wt-src/resources ./WtsRsrcs


# Now (optionally) build for simulator
mkdir build-iphonesimulator; cd build-iphonesimulator
unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../3rd_party/ios-cmake/ios.toolchain.cmake -DPLATFORM=SIMULATOR64 -DIOS_DEPLOYMENT_TARGET=16.4 -DDEPLOYMENT_TARGET=16.4 -DENABLE_BITCODE_INT=OFF -DENABLE_BITCODE=OFF ..
cmake --build . --config Release -j8
cp -R _deps/wt-src/resources ./WtsRsrcs
```

Now open InterSpec/InterSpec.xcodeproj and build the app.  
Note that the Xcode project has a "Targets" -> "InterSpec" -> "Build Phases" -> "Run Script" defined that triggers a `cmake --build` command for the `build-iphoneos` or `build-iphonesimulator` directory, so if the InterSpec library code was updated, the library will be re-built when you hit the build button in Xcode for the app project.


TODO: need to modify Wts CMakeLists.txt to not set the install directory
TODO: try adding a "-G Xcode" to the above, and then try integrating into iOS project, like at https://blog.tomtasche.at/2019/05/how-to-include-cmake-project-in-xcode.html - maybe cleaner that what I did now.







