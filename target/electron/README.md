# Creating a Distributable Electron Package

Currently, the [Electron]([https://electronjs.org/](https://electronjs.org/) packaged version of `InterSpec` starts the `InterSpec` server as a sub-process from [main.js](app/main.js).  So to create a Electron packaged version of `InterSpec` you must first build the server code using the CMake based build system, then package that into the Electron package, using the standard [node.js]([https://nodejs.org/en/](https://nodejs.org/) tools.

The below shows how to do this on macOS and Linux; Windows is essentually the same process, but you would use the CMake GUI and Visual Studio to build things (building using MingGW, Clang, or other build tools on Windows is not tested), and the npm terminal to package things.

```bash
#The CMake build directory must be named 'build_electron' for
#  the packaging scripts to work
cd /path/to/InterSpec
mkdir build_electron
cd build_electron

#For macOS only, you may want to define a deployment target
export MACOSX_DEPLOYMENT_TARGET=10.10

#Run CMake, letting it know you want to build the electron app.  
#  Specifying Boost and Wt paths is only necessary if they aren't 
#  in the standard locations
cmake -DBUILD_AS_ELECTRON_APP:BOOL="ON" \
      -DINCLUDE_ANALYSIS_TEST_SUITE:BOOL="OFF" \
      -DBOOST_ROOT=/path/to/boost \
      -DWt_INCLUDE_DIR=/path/to/wt/include \
      -DCMAKE_BUILD_TYPE="Release" \
      ..

#Build the InterSpec server executable; this will get packed 
#  into an Electron app
make -j8

#Change to electron target directory
cd ../target/electron/app

#(optional) Install electron and electron-packager globally on your system for development
npm install electron -g
npm install electron-packager -g

#Install the packages needed to make a distrbutable
#  You may need to edit package.json to get dependancies to work out
#  or manually locally install some of them 
npm install package.json

#If all has gone well, open an electron window and run InterSpec
npm start

#To create a distributable electron package, run one of the following
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
