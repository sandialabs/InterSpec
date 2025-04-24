# Build Deployable Image:
#   docker build --tag 'interspec_alpine' -f alpine_web_container.dockerfile .

# Run Prebuilt Image:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078/tcp interspec_alpine
# Or if you dont care about keeping user preferences and stuff around, you can just just map /data to a temp ephemeral dir
#   docker run --rf -v /data -p 8078:8078/tcp interspec_alpine

# Run From Dockerfile Directly:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078 -f alpine_web_container.dockerfile .

#  Optional: reuse existing build and/or source directories: Ensure ./build or ./src directory are populated and uncomment the COPY line

# If build fails, add the --no-cache option to the docker build command to force a fresh build
#  docker build --no-cache --tag 'interspec_alpine' -f alpine_web_container.dockerfile .

FROM alpine:3 AS build_interspec
ARG  repo=https://github.com/sandialabs/InterSpec.git
ARG  tag=master
WORKDIR /work

# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build

ARG  repo=https://github.com/sandialabs/InterSpec.git
ARG  tag=master
# RUN statements are broken up to allow loading cached images for debugging
RUN  apk add --no-cache \
        alpine-sdk \
        cmake \
        patch \
        linux-headers \
        suitesparse-dev \
        curl \
        uglify-js \
        uglifycss \
        git && \
    if [ ! -d ./src ]; then \
        git clone --recursive --branch $tag --depth=1 $repo ./src; \
    fi
RUN  cmake \
        -B ./build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FOR_WEB_DEPLOYMENT=ON \
        -DUSE_REL_ACT_TOOL=ON \
        -DBUILD_AS_LOCAL_SERVER=OFF \
        -DInterSpec_FETCH_DEPENDENCIES=ON \
        -DBoost_INCLUDE_DIR=./build/_deps/boost-src/libs \
        -DUSE_SEARCH_MODE_3D_CHART=ON \
        -DUSE_QR_CODES=ON \
        -DUSE_DETECTION_LIMIT_TOOL=ON \
        -DUSE_BATCH_TOOLS=OFF \
        -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" \
        -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
        -DZLIB_INCLUDE_DIR="_deps/zlib-build" \
        ./src
RUN  mkdir -p /InterSpec && \
     cmake --build build -j4
RUN  cmake --install ./build --prefix ./InterSpec && \
        rm -rf ./InterSpec/lib ./InterSpec/include ./InterSpec/share/eigen3/

#Web Server
FROM scratch
LABEL app="InterSpec"
COPY --from=build_interspec /work/InterSpec /interspec/
WORKDIR /interspec
EXPOSE 8078
ENTRYPOINT ["./bin/InterSpec", "--config=./share/interspec/data/config/wt_config_web.xml", "--userdatadir=/data", "--http-port=8078", "--http-address=0.0.0.0", "--docroot", "./share/interspec"]





# Then numeric group/user value of 280 was chosen randomly; it doesnt conflict with existing groups/users on dev or public server, and is below 1000 (e.g., a system user without a home directory or default shell)
#RUN groupadd --gid 280 interspec && useradd --uid 280 --gid interspec interspec
#RUN addgroup -S interspec && adduser --disabled-password --no-create-home -S interspec -G interspec
#USER interspec
# Or just use user guest
#USER guest
# You could keep the access log by chenging the entrypoint to: "--accesslog=/mnt/interspec_data/wt_access_log.txt"
# You could also edit the <log-file></log-file> element of data/config/wt_config_web.xml to save the stdout/stderr of InterSpec to a log file at /mnt/interspec_data/interspec_log.txt.