# Instructions 20250422 for building InterSpec using FetchContent for serving over the web
Please note that security is not consided at all (e.g., InterSpec is currently running as root in container), and no testing has been done for the web-deployment build.

The [alpine_web_container.containerfile](alpine_web_container.containerfile) builds and packages a statically linked version of the InterSpec executable and app resources, then creates a very minimal container to run this from.  No testing has been done with this `musl` libc-based build, beyond it seems to work briefly checking a few spectra.


The [debian_web_container.dockerfile](debian_web_container.dockerfile) is likely a more trustworthy build, as it uses `glibc` for the c library, which is what the Electron based Linux builds of InterSpec use.  However, again, these containers have not been tested beyond "it seems to work".

## Building the Alpine or Debian containers

You can either use the InterSpec code from your local machine, or have the container clone the code from GitHub, and build from that.

To build using your local code, you just need to change your CWD to the base of your InterSpec code, then run the build script; ex.:
```bash
cd /path/to/InterSpec
docker build -t interspec -f ./target/docker/alpine_web_container.containerfile .

# Or for Debian
docker build -t interspec_debian -f ./target/docker/debian_web_container.containerfile .
```
If you dont start the build from the InterSpec codes base directory, the InterSpec GitHub repository will be cloned inside the build container.


There are two options you can specify to docker:
- `hostdepdir` (true/false, default false): a `build` directory will be created in your how CWD, and everything in `./build/_deps/*-src` will be reused.  This can be useful if you are making changes, and dont want to download the dependencies for each build.
- `hostbuilddir` (true/false, default false): if specified `true`, and `hostdepdir=true`, then will also keep intermediate build files on host filesystem.


After building the container, you can then run the container with a command like:
```bash
docker run --rm -v /path/to/save/user/data:/data -p 8078:8078/tcp interspec

# But also, you should probably limit things a little:
docker run --rm -v /path/to/save/user/data:/data --read-only --cap-drop=all --memory=8g --cpus="4" --security-opt=no-new-privileges -p 8078:8078/tcp interspec
```
You can then connect your host web browser to http://localhost:8078 to use InterSpec.



## Debugging the build
```bash
# Start the docker container, with debug permissions, mapping InterSpec code directory from host to /interspec:
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -p 8078:8078 -i -t -v `pwd`:/interspec interspec bash

# Run the executable
./bin/InterSpec.exe -c ./share/interspec/data/config/wt_config_web.xml --userdatadir=/data --http-port=8078 --http-address=0.0.0.0 --docroot=./share/interspec
```


## Manually building things 
If you want to manually build the InterSpec executable inside a Docker container, you can do something like:
```bash
git clone --recursive --depth 1 --branch https://github.com/sandialabs/InterSpec.git master ./InterSpec_code
cd InterSpec_code
docker run -it -v `pwd`:/work/src -p 127.0.0.1:8078:8078/tcp --platform linux/arm64 alpine:latest sh
# if you dont have your proxy certificates setup right, and you're brave, you may need to add '--no-check-certificate' to the next couple lines.
apk update
apk add --no-cache alpine-sdk cmake patch linux-headers suitesparse-dev patch curl uglify-js uglifycss git

# If you want to CMake Fetch to build the dependencies, you can:
cmake -B ./build_alpine -DInterSpec_FETCH_DEPENDENCIES=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_FOR_WEB_DEPLOYMENT=ON -DUSE_REL_ACT_TOOL=ON -DBUILD_AS_LOCAL_SERVER=OFF -DBoost_INCLUDE_DIR=./build/_deps/boost-src/libs -DUSE_SEARCH_MODE_3D_CHART=ON -DUSE_QR_CODES=ON -DUSE_DETECTION_LIMIT_TOOL=ON -DUSE_BATCH_CLI_TOOLS=OFF -DUSE_BATCH_GUI_TOOLS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ./src/

# Or if you prefer, you can build a "prefix" you can use in later Docker sessions, via:
#./target/patches/dep_build_linux.sh /work/src/ /tmp/ /work/src/alpine_prefix
# cmake -B ./build_alpine -DInterSpec_FETCH_DEPENDENCIES=OFF -DCMAKE_PREFIX_PATH=/work/src/alpine_prefix .....

cmake --build build_alpine -j8


cmake --install ./build_alpine --prefix ./InterSpec
./InterSpec/bin/InterSpec --config src/data/config/wt_config_web.xml --http-address 0.0.0.0 --http-port 8078 --docroot InterSpec/share/interspec/ --userdatadir /tmp
```


## If you want to rebuild from scatch
```bash
docker build --no-cache -t interspec -f alpine_web_container.dockerfile .

# or to clean things up
docker rmi $(docker images -f "dangling=true" -q)
docker rmi build_interspec
docker rmi interspec

# or to be really agressive:
docker system prune
```