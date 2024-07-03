To build for Android we will take advantage of the CMake `FetchContent` style build of InterSpec to build boost and Wt for us.

So you should just need to:
- Open `InterSpec/target/android/InterSpec` in Android Studio, and adjust paths
   in the gradle build files; I wouldnt be suprised if you had to adjust the main 
   InterSpec CMakeLists.txt too.
- The first build will take ~30 minutes to clone into boost, and Wt, and other dependancies, as well as patch and build them.
- Android Studio tends to be extremely slow and unresponsive (like hours to do an operation - maybe because of the amount of NDK building?), so
   building and loading to device from the command line tends to be much quicker.
   - ```bash
     cd target/android/InterSpec
     export JAVA_HOME="/Applications/Android Studio.app/Contents/jbr/Contents/Home/"
     ./gradlew assembleUniversal
     ./gradlew installUniversalRelease
     #or to make iterating a little faster (min of 1.5 minutes instead of 5 minutes), use just the ARM release
     ./gradlew assembleArm7 && ./gradlew installArm7Release
     ```
- Good luck building things...








