# Building for iOS
To build for iOS we will take advantage of the CMake `FetchContent` style build of InterSpec to build boost and Wt for us.
What we'll do is build the InterSpec library, and then the InterSpec/InterSpec.xcodeproj Xcode project will build the app and link to the InterSpec library.  

The minimum supported iOS version is 13.0 (required because SpecUtils' bundled inja templating header uses `std::filesystem`, which Apple marks unavailable on iOS < 13.0).  The Xcode project's `IPHONEOS_DEPLOYMENT_TARGET` (currently 16.4) and the cmake `CMAKE_OSX_DEPLOYMENT_TARGET` must match; if you bump one, bump the other.

This uses CMake's built-in iOS support (`CMAKE_SYSTEM_NAME=iOS`), so no toolchain file is needed.

```bash
# Paths below assume the repo root is on disk; if your top-level directory is
# named something other than "InterSpec" (e.g. "InterSpec_master"), substitute
# accordingly.
cd InterSpec/target/ios/InterSpec

# First we'll build for actual devices
mkdir build-iphoneos; cd build-iphoneos

unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_NAME=iOS -DCMAKE_SYSTEM_PROCESSOR=arm64 -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=16.4 ..

# Initial configure + build pulls boost, Wt, ceres-solver, etc. via FetchContent and takes ~15-30 min.
cmake --build . --config Release -j8

# Wt's web assets: the root CMake leaves a copy at build-iphoneos/LibInterSpec/resources
# (with PRUNE_WT_RESOURCES=ON it is the slimmed set, see cmake/PruneWtResources.cmake), and
# target/ios/InterSpec/CMakeLists.txt mirrors it to build-iphoneos/WtsRsrcs.  A "Run Script"
# build phase copies WtsRsrcs into the .app.  It must NOT be named "resources": that collides
# case-insensitively with a bundle's reserved "Resources" directory and makes codesign fail
# with "bundle format unrecognized, invalid, or unsuitable".  data/config/wt_config_ios.xml
# sets resourcesURL=WtsRsrcs/ to match (Android keeps resources/ via wt_config_mobile.xml).


# Now (optionally) build for simulator
cd ..
mkdir build-iphonesimulator; cd build-iphonesimulator
unset MACOSX_DEPLOYMENT_TARGET
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_SYSTEM_NAME=iOS -DCMAKE_OSX_SYSROOT=iphonesimulator -DCMAKE_SYSTEM_PROCESSOR=arm64 -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=16.4 ..
cmake --build . --config Release -j8
```

Now open `InterSpec/target/ios/InterSpec/InterSpec.xcodeproj` and build the app.
The Xcode project has a "Targets" -> "InterSpec" -> "Build Phases" -> "Run Script" build phase that runs `cmake --build` against the matching `build-iphoneos` or `build-iphonesimulator` directory (selected via `${EFFECTIVE_PLATFORM_NAME}`), so subsequent edits to the InterSpec library sources will be rebuilt when you hit Build in Xcode.

If you change the iOS minimum deployment target, you need to update it in three places to keep them in sync:
  - `CMAKE_OSX_DEPLOYMENT_TARGET` in the cmake invocations above (and re-configure / rebuild from a clean build dir).
  - `IPHONEOS_DEPLOYMENT_TARGET` in `InterSpec/InterSpec.xcodeproj/project.pbxproj` (currently four occurrences: two project-level Debug/Release, two target-level Debug/Release).
  - Any signing/profile constraints in Xcode that may be tied to a minimum OS.


TODO: need to modify Wts CMakeLists.txt to not set the install directory
TODO: try adding a "-G Xcode" to the above, and then try integrating into iOS project, like at https://blog.tomtasche.at/2019/05/how-to-include-cmake-project-in-xcode.html - maybe cleaner that what I did now.







