# Creating a Distributable Electron Package

In order to create the [Electron]([https://electronjs.org/](https://electronjs.org/) packaged version of `InterSpec`, we create a node add-on that contains all the C++ code, and then request to start the `InterSpec` server from [main.js](app/main.js).  To create the add-on we use [cmake-js](https://www.npmjs.com/package/cmake-js) to build things, and [node-addon-api](https://www.npmjs.com/package/node-addon-api) to actually interface the C++ to JS.


```bash
npm install -g cmake-js 

# For macOS only, you may want to define a deployment target
export MACOSX_DEPLOYMENT_TARGET=10.12

cd /path/to/InterSpec/target/electron

# Install dependency for compiling a node.js add-on
npm install --save-dev node-addon-api

# If boost and Wt are in standard locations, you can just run
cmake-js

# Or to have a little more control over things
cmake-js --generator "Visual Studio 15 2017 Win64" \
         --CDBOOST_ROOT=/path/to/boost \
         --CDWt_INCLUDE_DIR=/path/to/wt/include \
         --CDCMAKE_BUILD_TYPE="Release" \
         --CDGOOGLE_MAPS_KEY="..." \
         --out="build_dir" \
         --target install

# If you make changes and want to recompile
cmake-js build --out="build_dir"
# Or you can use CMake to run the `make` command, which can be useful when
# the CMake generator isnt a command-line based system like (ex Xcode, MSVC)
cmake --build build_dir --config Release
# Or directly use the `make` command like:
ninja -C build_dir

# Make the "app" directory into the build directory
ninja -C build_dir install

# And also copy all the InterSpec resources to the 
# 'app' sub-directory of your build dir
cmake-js build --out=build_dir --target install
# Or
cmake --build build_dir --target install --config Release --clean-first
# Or
ninja -C build_dir install


# To run InterSpec without packaging everything, you
# need to install the Electron package
npm install electron

# Then to actually run things
npm start


# To create a packaged version of InterSpec, you
# need to install electron-packager
npm install electron-packager


# Then to actually create the distributable package, run on of the following
# (assuming your build directory is 'build_macos', 'build_win', or 'build_linux', 
#  otherwise you need to edit path in package.json)
npm run package-mac
npm run package-win
npm run package-linux
```

The resulting Electron package seems to run fine on a number of Linux distributions, but to run on the (now out of date) Ubuntu 12.04, you may need to install a couple packages from the package manager on the end-user system:

```bash
apt-get install libnss3
apt-get install libxss1
```

If you are cross-compiling, you can, for example, build the Linux package from macOS using a command like 

```bash
npm i --target_arch=x64 --target_platform=linux ws
```

## Linux Considerations
The `InterSpec` module is really a shared library that node.js loads, therefore you need the `-fPIC` C/C++ compiler flag enabled not just for the `InterSpec` code, but for all of the static libraries you link it against, including boost, Wt, and zlib - which isnt the default when compiling static libraries for any of them, so when building them you may want to add `-fPIC -std=c++11` flags to the flags.


## Future Work
- In the future the build process may be improved to do the final packaging through CMake.  
- launch_options.json is intended to allow users to customize behaviors (like allow multiple windows, or not), but is not actually implemented.
- Need to look into having only a single instance of the application running when user wants multiple windows open.