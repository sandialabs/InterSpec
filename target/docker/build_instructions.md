# Instructions 20220422 for building InterSpec using FetchContent for serving over the web
Please note that security is not consided at all (e.g., InterSpec is currently running as root in container), and no testing has been done for the web-deployment build.

## Building using a Alpine container
```bash
# From your host OS terminal, run the following commands
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_alpine_build
cd InterSpec_alpine_build

# Grab the Alpine image to your machine
docker pull alpine:latest

# Start a shell session within the image, mapping the InterSpec source 
#  directory to /interspec.  We'll also map port 8082 for testing.
docker run --rm -it -v `pwd`:/interspec -p 127.0.0.1:8082:8082/tcp alpine:latest sh

# Get the dependancies we need to build InterSpec
apk add --update alpine-sdk cmake patch linux-headers

# Optional to minimize JS/CSS
# apk add npm
# npm install uglify-js -g
# npm install uglifycss -g

# Make and cd into build directory - note this is in host filesystem incase we  
#  want to come back to things, but dont want to rebuild everything from scratch
cd interspec/
mkdir build_docker
cd build_docker

# This next command will take like 10 or 20 minutes to clone into the boost and Wt repositories
#  Also note that it will build a static executable, that *should* run on nearly any linux (but unteted)
cmake -DUSE_SQLITE3_DB=ON -DBUILD_FOR_WEB_DEPLOYMENT=ON -DInterSpec_FETCH_DEPENDENCIES=ON -DBUILD_AS_LOCAL_SERVER=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++ -static" ..

# And actually building everything will take maybe an hour the first time; after
#  that each rebuild should only take a minute or less.
cmake --build . --config Release -j6

# Grab html resources needed by Wt
cp -r _deps/wt-src/resources .

# Test running the newly compiled InterSpec executable (connect on your host OS browser at http://127.0.0.1:8082)
./InterSpec --docroot=. --http-address=0.0.0.0 --http-port=8082 --config=data/config/wt_config_localweb.xml --accesslog=- 


# Now create a directory with all the InterSpec resources you need
mkdir interspec_install
mkdir interspec_install/html_root
cp InterSpec interspec_install
cp -r ../data ../InterSpec_resources ../example_spectra _deps/wt-src/resources interspec_install/html_root/
```

Now to make a self-contained container - back in your host OS terminal:
```bash
# Get rid of old containers:
docker rmi interspec_alpine_web

# Build the container run the following from the InterSpec base-dir:
docker build -t 'interspec_alpine_web' -f target/docker/alpine_web_container.dockerfile .

# Run the InterSpec docker container, we need to supply a mount point (maps to /mnt/interspec_data in the container).
#  We can either mount our local filesystem, to save things into a directory, or we can use a docker volume
#  (you may want to adjust memory and cpu limits to fit your system)
#  To save InterSpec data to our local fileystem, run:
docker run --read-only --cap-drop=all --memory=8g --cpus="4" --security-opt=no-new-privileges --rm -v /tmp -v /path/on/your/fs/to/persist/user/data:/mnt/interspec_data  -p 127.0.0.1:8078:8078/tcp interspec_alpine_web

# To save data to a docker volume, named for example "interspec_data", run the command:
#  (although just niavely I ran into trouble with volume being mounted read-only...didnt debug)
docker run --read-only --cap-drop=all --memory=8g --cpus="4" --security-opt=no-new-privileges --rm -v /tmp -it --mount source=interspec_data,target=/mnt/interspec_data  -p 127.0.0.1:8078:8078/tcp interspec_alpine_web

# Or to use a efermeral temp drive
docker run --rm -v /mnt/interspec_data -p 127.0.0.1:8078:8078/tcp interspec_alpine_web


# To save the container to a file, run:
docker save -o interspec_alpine_web.20220424.docker interspec_alpine_web

# To load the container from a  file, run:
docker load -i interspec_alpine_web.20220424.docker
```





## Building using a Fedora container

```bash
# From your host OS terminal, run the following commands
git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_fedora_build
cd InterSpec_fedora_build

docker pull fedora:35

# Start a shell session within the image
docker run --rm -it -v `pwd`:/interspec -p 127.0.0.1:8082:8082/tcp fedora:35 sh

# Install dependancies needed to compile InterSpec
dnf install make automake gcc gcc-c++ kernel-devel cmake git patch

cd interspec/
mkdir build_docker
cd build_docker

# This next command will take like 10 or 20 minutes to clone into the boost and Wt repositories
cmake -DUSE_SQLITE3_DB=ON -DBUILD_FOR_WEB_DEPLOYMENT=ON -DInterSpec_FETCH_DEPENDENCIES=ON -DBUILD_AS_LOCAL_SERVER=OFF -DCMAKE_BUILD_TYPE=Release ..

# Build everything; may take while first time
cmake --build . --config Release -j6

# Grab html resources needed by Wt
cp -r _deps/wt-src/resources .

# Test running the newly compiled InterSpec executable
./InterSpec --docroot=. --http-address=0.0.0.0 --http-port=8082 --config=data/config/wt_config_localweb.xml --accesslog=- --no-compression

# Now create a directory with all the InterSpec resources you need
mkdir interspec_install
mkdir interspec_install/html_root
cp InterSpec interspec_install
cp -r ../data ../InterSpec_resources ../example_spectra _deps/wt-src/resources interspec_install/html_root/
```


Now to make a self-contained container - back in your host OS terminal:
```bash
# Get rid of old containers:
docker rmi interspec_web

# Build the container run the following from the InterSpec base-dir:
docker build -t 'interspec_web' -f target/docker/fedora_web_container.dockerfile .

# Run the InterSpec docker container:
docker run -p 127.0.0.1:8078:8078/tcp interspec_web

# To save the container to a file, run:
docker save -o interspec_web.20220422.docker interspec_web

# To load the container from a  file, run:
docker load -i interspec_web.20220422.docker
```

# Outdated instrcutions for building under Ubuntu, with manually compiled prerequisits

```bash
# Step 1) build and install boost 1.65.1 and Wt 3.3.4 into an Ubuntu 18.0.4 base image
docker build -t interspec_ubuntu_build_env:v1.0 -f "setup_build_env_ubuntu.dockerfile" ../..

# Step 2) Compile InterSpec.
docker build -t interspec_ubuntu_build:v1.0 -f "build_interspec_ubuntu.dockerfile" ../..

# Step 3) Enter docker image with compiled InterSpec.  The below command will let you use gdb to debug if you want and it maps the images port 8080 to your machines port 8080.
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -p 8080:8080 -i -t -v /Users/wcjohns/rad_ana/InterSpec:/interspec interspec_ubuntu_build:v1.0 bash

# Step 4) From within docker, run InterSpec
cd /tmp/build_interspec/
./bin/InterSpec.exe --docroot . --http-address 0.0.0.0 --http-port 8080 -c ./data/config/wt_config_localweb.xml

# Step 5) Point your browser to http://localhost:8080/ and use InterSpec.
```