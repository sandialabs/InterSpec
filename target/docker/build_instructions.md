# Instructions 20250422 for building InterSpec using FetchContent for serving over the web
Please note that security is not consided at all (e.g., InterSpec is currently running as root in container), and no testing has been done for the web-deployment build.

The [alpine_build_interspec.dockerfile](alpine_build_interspec.dockerfile) builds and packages the executable, then create a container to run it from.


## Building the container
```bash
docker build -t interspec -f alpine_web_container.dockerfile .
```
Inside the container, the InterSpec code will be cloned from github, then it and all its dependencies will be built, and a container made with the results.

## Running the container
The default build will assume you have mapped a directory to store user data (preferences, incomplete work, etc)
into `/data` in the container, so you will need to map this:
```bash
docker run -v /path/to/save/data/on/your/fs:/data  --rm -v /tmp -p 8078:8078/tcp interspec

#Or if you dont care about keeping user preferences and stuff around, you can just just map /data to a temp ephemeral dir
docker run --rm -v /data  --rm -v /tmp -p 8078:8078/tcp interspec

#But also, you should probably limit things a little:
docker run --read-only --cap-drop=all --memory=8g --cpus="4" --security-opt=no-new-privileges --rm -v /data  --rm -v /tmp -p 8078:8078/tcp interspec
```

You should now be able to point the browser on your host machine to http://localhost:8078/ and use InterSpec.


## Debugging the build
If you want to manually run the executable to debug, you can:
```bash
# Start the docker container, with debug permissions:
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -p 8080:8080 -i -t -v /Users/wcjohns/rad_ana/InterSpec:/interspec interspec bash

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

cmake -B ./build -DCMAKE_BUILD_TYPE=Release -DBUILD_FOR_WEB_DEPLOYMENT=ON -DUSE_REL_ACT_TOOL=ON -DBUILD_AS_LOCAL_SERVER=OFF -DInterSpec_FETCH_DEPENDENCIES=ON -DBoost_INCLUDE_DIR=./build/_deps/boost-src/libs -DUSE_SEARCH_MODE_3D_CHART=ON -DUSE_QR_CODES=ON -DUSE_DETECTION_LIMIT_TOOL=ON -DUSE_BATCH_TOOLS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ./src/

cmake --build build -j8


cmake --install ./build --prefix ./InterSpec
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