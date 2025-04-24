FROM debian:latest AS build
LABEL app="interspec_build_debug"
ENV CONFIGURE=False
ENV BUILD=FALSE
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
CMD ["sh", "-c", \
                "if [ \"$CONFIGURE\" ]; then \
                    echo CONFIGURE is True, configuring InterSpec && \
                    cmake \
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
                        -S ./src ;\
                else \
                    echo CONFIGURE not set or False, skipping configure;\
                fi &&\
                if [ \"$BUILD\" ]; then \
                    echo BUILD is True, building InterSpec &&\
                    mkdir -p release && \
                    cmake --build build -j4 && \
                    cmake --install ./build --prefix release && \
                    rm -rf release/lib/cmake && \
                    cd ./release &&\ 
                    chmod -R a+r * &&\
                    chmod a+x bin/InterSpec &&\
                    mkdir -p data &&\
                    chmod 777 data &&\
                    echo Build complete. ;\
                else \
                    echo BUILD not set or False, skipping build; \
                fi &&\
                ./bin/InterSpec --config ./share/interspec/data/config/wt_config_web.xml --userdata=./ --http-port=8078 --http-address=0.0.0.0 --docroot ./share/interspec" \
    ]
#Web Server

#ls -al
#docker build --tag 'interspec_build_debug' -f src/target/docker/debian.dockerfile .

#Put the InterSpec repo in .src/, and the following container will put the compiled version in ./release
#docker run -e CONFIGURE=True -e BUILD=True --rm -v "$PWD":/work -it -p 8078:8078 interspec_build_debug

# To use a pre-existing build directory and skip the configure step use:
# docker build --tag 'interspec_alpine_compileonly' -f target/docker/alpine_web_container_compileonly.dockerfile .
# STARTSERVER=True BUILDONLY=TRUE docker run --rm -v "$PWD":/work -v ~/appdir:/InterSpec interspec_alpine_compileonly

