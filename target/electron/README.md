# Creating a Distributable Electron Package

In order to create the [Electron]([https://electronjs.org/](https://electronjs.org/) packaged version of `InterSpec`, we create a node add-on that contains all the C++ code, and then request to start the `InterSpec` server from [main.js](app/main.js).  To create the add-on we use [cmake-js](https://www.npmjs.com/package/cmake-js) to build things, and [node-addon-api](https://www.npmjs.com/package/node-addon-api) to actually interface the C++ to JS.


```bash
npm install -g cmake-js

# For macOS only, you may want to define a deployment target
export MACOSX_DEPLOYMENT_TARGET=10.10

cd /path/to/InterSpec/target/electron

# Install dependency for compiling a node.js add-on
npm install --save-dev node-addon-api

# If boost and Wt are in standard locations, you can just run
cmake-js

# Or to have a little more control over things
cmake-js --CDBOOST_ROOT=/path/to/boost \
         --CDWt_INCLUDE_DIR=/path/to/wt/include \
         --CDCMAKE_BUILD_TYPE="Release" \
         --out="build_dir"

# If you make changes and want to recompile
cmake-js build --out="build_dir"
# Or
ninja -C build_dir

# To run InterSpec without packaging everything, you
# need to install the Electron package
npm install electron

# Then to actually run things
npm start


# To create a packaged version of InterSpec, you
# need to install electron-packager
npm install electron-packager

# And also copy all the InterSpec resources to the 
# 'app' sub-directory of your build dir
ninja -C build_dir install

# Then to actually create the distributable package, run on of the  following
npm run package-mac
npm run package-win
npm run package-linux
```

The resulting Electron package seems to run fine on a number of Linux distributions, but to run on the (now out of date) Ubuntu 12.04, you may need to install a couple packages from the package manager on the end-user system:

```bash
apt-get install libnss3
apt-get install libxss1
```

If you are cross-compiling, you can build the Linux package from macOS using a command like 

```bash
npm i --target_arch=x64 --target_platform=linux ws
```

## Future Work

In the future the build process may be improved to reduce the number of steps.

In the future it may also be possible to build the `InterSpec` server code as a native node.js module for a tighter, cleaner, integration between the GUI window and the server code.
