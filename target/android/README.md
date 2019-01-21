Building for Android is a pain; you have to make sure the same compiler, flags,
runtime, etc are used for boost, Wt, and InterSpec (or at least that everything
is abi compatible) - if not, things might still compile, but randomly crash.  
And of course the boost build system is bjam, while Wt and InterSpec is CMake, 
and as it is now, InterSpec is being built by Android Studio, while Wt is from 
the command line.  Further, for the life of me I couldnt get CMake to actually
find boost and Wt libraries for android, even when I specified all the 
directories I could. I was unable to find an elegant way to do all this, so I 
just brute forced it.
These instructions are for May of 2018, and macOS was used.

-Download and install Android Studio.
-Then use the HelloJNI example project to determine which compile flags is used
 by cmake/Android-Studio.
  --First in build.gradle set options to what you want ('minSdkVersion 19',
    'targetSdkVersion 23', 'compileSdkVersion 25', and CMake arguments
    "-DANDROID_CPP_FEATURES=rtti exceptions", "-DANDROID_TOOLCHAIN=clang",
    "-DANDROID_PIE=ON", "-DANDROID_STL=c++_static").
  --Then build each architecture you might want, and then in like
    /Users/wcjohns/AndroidStudioProjects/HelloJNI/app/.externalNativeBuild/cmake/x86_64Release/x86_64
    open android_gradle_build.json and get the actual compile flags.
-Then use https://github.com/moritz-wundke/Boost-for-Android to build boost
  --Boost-for-Android is great, but you really need to check that it uses
    the same build settings as CMake will.  To do this edit the
    configs/user-config-boost-1_65_1-common.jam and
    configs/user-config-boost-1_65_1-{target}.jam target files that come with
    Boost-for-Androidso to the compile flags match the ones extracted from
    HelloJNI.  You need to do this for each architecture.
  --Also, you may need to edit build-android.sh to allow using NDK version 17
    (using same settings as 16)
  --run `./build-android.sh --boost=1.65.1 --toolchain=llvm --arch=armeabi-v7a,arm64-v8a,x86,x86_64 --prefix=/Users/wcjohns/install/android/ /Users/wcjohns/Library/Android/sdk/ndk-bundle`
-Then compile Wt 3.3.4.
  --Replace a few Wt files: target/android/wt_patch/wt-3.3.4/Android.C
    target/android/wt_patch/wt-3.3.4/CMakeLists.txt
    target/android/wt_patch/wt-3.3.4/base64.cpp
  --Manually build Wt for each architecture
    --```bash
      mkdir build_android_armeabi-v7a;cd build_android_armeabi-v7a`
      
      export MY_ANDROID_ABI="armeabi-v7a"
      export MY_ANDROID_PLATFORM="android-19"
      export MY_BOOST_DIR="/Users/wcjohns/install/android/${MY_ANDROID_ABI}"
      cmake -DANDROID_ABI=${MY_ANDROID_ABI} -DANDROID_PLATFORM=${MY_ANDROID_PLATFORM} \
            -DCMAKE_BUILD_TYPE=Release -DANDROID_NDK=/Users/wcjohns/Library/Android/sdk/ndk-bundle \
            -DCMAKE_CXX_FLAGS=-std=c++11 \
            -DANDROID_STL=c++_static \
            -DCMAKE_TOOLCHAIN_FILE=/Users/wcjohns/Library/Android/sdk/ndk-bundle/build/cmake/android.toolchain.cmake \
            -DANDROID_CPP_FEATURES=rtti exceptions -DANDROID_TOOLCHAIN=clang -DANDROID_PIE=ON \
            -DSHARED_LIBS=OFF -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF -DWT_CPP_11_MODE='-std=c++11' \
            -DENABLE_SSL=OFF \
            -DBoost_INCLUDE_DIR=${MY_BOOST_DIR}/include/boost-1_65_1 \
            -DBoost_DATE_TIME_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_date_time-clang-darwin-mt-1_65_1.a \
            -DBoost_FILESYSTEM_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_filesystem-clang-darwin-mt-1_65_1.a \
            -DBoost_PROGRAM_OPTIONS_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_program_options-clang-darwin-mt-1_65_1.a \
            -DBoost_RANDOM_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_random-clang-darwin-mt-1_65_1.a \
            -DBoost_REGEX_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_regex-clang-darwin-mt-1_65_1.a \
            -DBoost_SIGNALS_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_signals-clang-darwin-mt-1_65_1.a \
            -DBoost_SYSTEM_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_system-clang-darwin-mt-1_65_1.a \
            -DBoost_THREAD_LIBRARY_RELEASE=${MY_BOOST_DIR}/lib/libboost_thread-clang-darwin-mt-1_65_1.a \
            -DCMAKE_INSTALL_PREFIX=${MY_BOOST_DIR} \
            ..
      make -j8
      make install
      ```
      And of course you probably want to do this for the other architectures:
      ```bash
      export MY_ANDROID_ABI="arm64-v8a"
      export MY_ANDROID_PLATFORM="android-21"
      export MY_BOOST_DIR="/Users/wcjohns/install/android/${MY_ANDROID_ABI}"
      cmake ...
      ```

      ```
      export MY_ANDROID_ABI="x86"
      export MY_ANDROID_PLATFORM="android-19"
      export MY_BOOST_DIR="/Users/wcjohns/install/android/${MY_ANDROID_ABI}"
      cmake ...
      ```

      ```
      export MY_ANDROID_ABI="x86_64"
      export MY_ANDROID_PLATFORM="android-21"
      export MY_BOOST_DIR="/Users/wcjohns/install/android/${MY_ANDROID_ABI}"
      cmake ...
      ```

    --However at this point libwt.a is like 700 MB, and doesnt have the correct 
      headers, so you have to run:
      ```bash
      ${NDK_ROOT}/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-strip --strip-debug libwt.a`
      ```
      and 
      ```
      ${NDK_ROOT}/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-ranlib libwt.a
      ```
      But of course you will have to do this for each Wt library, and architecture 
      (make sure to use the correct ranlib and strip for each architecture though).
      You might also have to run ranlib over the boost libraries.
  -Open `InterSpec/target/android/InterSpec` in Android Studio, and adjust paths
   in the gradle build files; I wouldnt be suprisedif you had to adjust the main 
   InterSpec CMakeLists.txt too.
  -Update the files such as `data`, `InterSpec_resources`, `resources`, and 
  `example_data` by putting them into a zip file 
  `InterSpec/target/android/InterSpec/app/src/main/assets/interspec-assets.zip`
  -Good luck building things...









