FROM alpine:latest
LABEL app="interspec_alpine_compileonly"
ENV BUILDONLY=null
ENV STARTSERVER=null
WORKDIR /work
EXPOSE 8078
RUN  apk add --no-cache \
        alpine-sdk \
        cmake \
        patch \
        linux-headers \
        suitesparse-dev patch \
        curl \
        uglify-js \
        uglifycss \
        git
CMD ["/bin/sh", "-c", \
                "if [ -z \"$BUILDONLY\" ]; then         \
                    echo BUILDONLY not set, configuring InterSpec;         \
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
                        -DUSE_SPECRUM_FILE_QUERY_WIDGET=ON \
                        -DSpecUtils_PYTHON_BINDINGS=ON \
                        . ; \
                else \
                    echo Building InterSpec with build only option; \
                fi && \
                mkdir -p /InterSpec && \
                cmake --build build -j4 && \
                cmake --install ./build --prefix /InterSpec && \
                rm -rf /InterSpec/lib/cmake && \
                echo Build complete. && \
                if [ -n \"$STARTSERVER\" ]; then         \
                    cd /InterSpec &&\ 
                    chmod -R a+r * &&\
                    chmod a+x bin/InterSpec &&\
                    chmod 777 /InterSpec &&\
                    mkdir -p data &&\
                    chmod 777 data &&\
                    ./bin/InterSpec --config ./share/interspec/data/config/wt_config_web.xml --userdatadir=./data --http-port=8078 --http-address=0.0.0.0;\
                fi" \
            ]



# To use a pre-existing build directory and skip the configure step use:
# docker build --tag 'interspec_alpine_compileonly' -f target/docker/alpine_web_container_compileonly.dockerfile .
# STARTSERVER=True BUILDONLY=TRUE docker run --rm -v "$PWD":/work -v ~/appdir:/InterSpec interspec_alpine_compileonly


