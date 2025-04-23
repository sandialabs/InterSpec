FROM debian:latest AS build
LABEL app="interspec_debian"
ARG BRANCH="master"
ARG REPO=https://github.com/sandialabs/InterSpec.git
WORKDIR /work
EXPOSE 8078
RUN apt update && apt upgrade -y && apt install -y --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    patch \
    libsuitesparse-dev \
    libjs-uglify-js \
    git
RUN git clone --branch $BRANCH --depth=1 $REPO ./src
RUN cmake \
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
        -DSpecUtils_PYTHON_BINDINGS=ON \
        -S ./src
RUN cmake --build build -j4; \
    mkdir -p release && \
    cmake --build build -j4 && \
    cmake --install ./build --prefix release && \
    rm -rf release/lib/cmake && \
    cd ./release &&\ 
    chmod -R a+r * &&\
    chmod a+x bin/InterSpec &&\
    chmod 777 /InterSpec &&\
    mkdir -p data &&\
    chmod 777 data &&\
    echo Build complete. 
#Web Server
FROM debian:latest
LABEL app="InterSpec"
COPY --from=build /work/release /interspec
WORKDIR /interspec
EXPOSE 8078
RUN apt update && apt upgrade -y && apt install -y --no-install-recommends \
        openblas \
        libstdc++ \
        libgcc 
SHELL ["sh", "-c"]
ENTRYPOINT ["./bin/InterSpec", "--config ./share/interspec/data/config/wt_config_web.xml", "--userdatadir=/data", "--http-port=8078", "--http-address=0.0.0.0", "--docroot", "./share/interspec"]



#ls -al
#docker build --tag 'interspec_debian' -f src/target/docker/debian_web_container.dockerfile .
#docker run --rm -v "$PWD":/work -v interspec_debian

# To use a pre-existing build directory and skip the configure step use:
# docker build --tag 'interspec_alpine_compileonly' -f target/docker/alpine_web_container_compileonly.dockerfile .
# STARTSERVER=True BUILDONLY=TRUE docker run --rm -v "$PWD":/work -v ~/appdir:/InterSpec interspec_alpine_compileonly

