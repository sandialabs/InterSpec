# Building dependencies for InterSpec
When building InterSpec, you can either have CMake download and build all dependencies through the `InterSpec_FETCH_DEPENDENCIES` CMake option (e.g., `cmake -DInterSpec_FETCH_DEPENDENCIES=ON ..`), or you can build all dependencies and place them into a prefix location.  Having CMake fetch and build the dependencies is required for iOS and Android, while creating a prefix location with the dependencies pre-built is highly recommended for desktop OSes, especially when doing heavy development, as it is slightly faster to re-compile the InterSpec code against since the build system has much less files to track and check if they have been changed.  

Using package managers like `brew`, `pacman`, or `nuget` to install required is very unlikely to work - and not suggested.  


## Building with CMake Fetch 
You only need to add the `-DInterSpec_FETCH_DEPENDENCIES=ON` command-line option when configuring InterSpec code.
It may take 30 minutes or more to fetch and build the dependencies, but you dont need to manually build a prefix.
An example of doing this is:
```bash
git clone --recursive https://github.com/sandialabs/interspec/ InterSpec
cd InterSpec
mdkir build_fetch
cd build_fetch
cmake -DInterSpec_FETCH_DEPENDENCIES=ON -DBUILD_AS_LOCAL_SERVER=ON -DUSE_BATCH_CLI_TOOLS=OFF -DUSE_BATCH_GUI_TOOLS=OFF ..
cmake --build . -j4
./bin/InterSpec.exe --docroot . --http-address 127.0.0.1 --http-port 8080 -c ./data/config/wt_config_localweb.xml
```



## Creating a re-usable "prefix" directory for the dependencies
If you will be developing InterSpec, or iterating on builds, this is the recommended way to build the prerequisite.  When using this method, a `prefix` directory will be created that will contain all needed libraries to link InterSpec against, so when InterSpec is being built, it will not look in your systems default locations for libraries, and since the `prefix` is not your systems default path (e.g., `/usr/local/`), it wont affect anything else you might compile.

You can compile the dependencies using the following scripts for each of the operating systems:
- [dep_build_linux.sh](dep_build_linux.sh)
- [dep_build_macOS.sh](dep_build_macOS.sh)
- [dep_build_msvc2022.bat](dep_build_msvc2022.bat)  (requires MSVC 2019 or 2022)

Please open the scripts and read the directions at the top, but hte basic usuage is
```bash
dep_build_macOS.sh /path/to/source/code/for/InterSpec/ /some/temp/path/to/build/in/ /prefix/location/"
```

So to build for Linux from scratch, you could
```bash
cd /path/to/code
git clone --recursive https://github.com/sandialabs/interspec/ InterSpec
./InterSpec/target/patches/dep_build_linux.sh ./InterSpec /tmp/interspec_prefix_build /home/user/code/InterSpec_prefix

cd InterSpec
mkdir build
cd build

cmake -DCMAKE_PREFIX_PATH="/home/user/code/InterSpec_prefix" ..
cmake --build . -j4  #or `make -j4 `

# If you compiled it in debug mode, you may need symlink the D3.js based resources to where the web-requests expect to find them, e.g.
ln -s ../external_libs/SpecUtils/d3_resources ./external_libs/
# (this is so you can edit d3_resources/SpectrumChartD3.js and see changes without doing a file copy)

./bin/InterSpec.exe --docroot . --http-address 127.0.0.1 --http-port 8080 -c ./data/config/wt_config_localweb.xml

# You can now point your browser to http://127.0.0.1:8080
```



# Building on Windows
You will need to install Visual Studio 2022, or 2019 may work to.

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



# Building on Linux
You can now build a local web server:
```bash
cd ~/development
git clone --recursive https://github.com/sandialabs/interspec/
cd interspec
./InterSpec/target/patches/dep_build_linux.sh . /tmp/interspec_prefix_build ~/development/InterSpec_prefix
mkdir build; cd build

cmake -DCMAKE_PREFIX_PATH=~/development/InterSpec_prefix -DBUILD_AS_LOCAL_SERVER=ON -DUSE_BATCH_CLI_TOOLS=OFF -DUSE_BATCH_GUI_TOOLS=OFF ..
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
