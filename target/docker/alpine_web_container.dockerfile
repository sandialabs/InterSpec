# Build:
#  git clone --recursive --depth 1 --branch https://github.com/sandialabs/InterSpec.git master ./InterSpec_code
#  cd InterSpec_code
#  docker build -t interspec -v `pwd`:/work/src -p 127.0.0.1:8078:8078/tcp -f alpine_web_container.dockerfile .
#To run, you can use
#  docker run -v /path/on/your/fs/to/persist/user/data:/data -p 8078:8078/tcp interspec
#Or if you dont care about keeping user preferences and stuff around, you can just just map /data to a temp ephemeral dir
#  docker run --rf -v /data -p 8078:8078/tcp interspec


FROM alpine:latest AS build
WORKDIR /work
# RUN statements are broken up to allow loading cached images for debugging
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
        ./src
RUN  mkdir -p /InterSpec && \
     cmake --build build -j4
RUN  cmake --install ./build --prefix ./InterSpec && \
        rm -rf ./InterSpec/lib/cmake

#Web Server
FROM alpine:latest
LABEL app="InterSpec"
COPY --from=build /work/InterSpec /interspec/
WORKDIR /interspec
EXPOSE 8078
RUN apk --no-cache add \
        openblas \
        libstdc++ \
        libgcc && \
        chmod -R a+r * && \
        chmod a+x bin/InterSpec &&  \
        chmod 777 /interspec
SHELL ["/bin/sh", "-c"]
ENTRYPOINT ["./bin/InterSpec", "-c ./share/interspec/data/config/wt_config_web.xml", "--userdatadir=/data", "--http-port=8078", "--http-address=0.0.0.0", "--docroot", "./share/interspec"]


# Then numeric group/user value of 280 was chosen randomly; it doesnt conflict with existing groups/users on dev or public server, and is below 1000 (e.g., a system user without a home directory or default shell)
#RUN groupadd --gid 280 interspec && useradd --uid 280 --gid interspec interspec
#RUN addgroup -S interspec && adduser --disabled-password --no-create-home -S interspec -G interspec
#USER interspec
# Or just use user guest
#USER guest

# If build fails, issue this command
# docker rmi $(docker images -f “dangling=true” -q)
# docker system prune
# You could keep the access log by chenging the entrypoint to: "--accesslog=/mnt/interspec_data/wt_access_log.txt"
# You could also edit the <log-file></log-file> element of data/config/wt_config_web.xml to save the stdout/stderr of InterSpec to a log file at /mnt/interspec_data/interspec_log.txt.
