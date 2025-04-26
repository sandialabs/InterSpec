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
ARG repo=https://github.com/sandialabs/InterSpec.git
ARG tag=master
ARG savebuild=false
WORKDIR /work

# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build

ARG  repo=https://github.com/sandialabs/InterSpec.git
ARG  tag=master
# RUN statements are broken up to allow loading cached images for debugging
RUN  if [[ ${savebuild} != "false"  ]]; then --mount type=bind,source=${PWD},target=/work,rw; fi &&\
        apk add --no-cache \
        alpine-sdk \
        cmake \
        patch \
        linux-headers \
        suitesparse-dev \
        curl \
        uglify-js \
        uglifycss \
        git && \
    if [ ! -d ./InterSpec_resources ]; then \ 
        git clone --recursive --branch $tag --depth=1 $repo .; \
    elif [[ ${savebuild} != "false"  ]]; then \
        git clone --recursive --branch $tag --depth=1 $repo .; \
    else \
        mkdir -p /tmp/src && \
        git clone clone --recursive --branch $tag --depth=1 $repo . && \
        cp -r /tmp/src/* ./ ; \
    fi && \
    cmake \
        -B ./build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FOR_CONTAINER_LIBC=ON \
        -S . && \
    cmake --build build -j4; \
    mkdir -p release && \
    cmake --build build -j4 && \
    cmake --install ./build --prefix /release && \
    rm -rf release/lib/cmake && \
    cd /release && \ 
    chmod -R a+r * && \
    chmod a+x bin/InterSpec && \
    mkdir -p data && \
    chmod 777 data && \
    if [[ "$savebuild" != "false"  ]]; then \
        cp -r /release /work/release; \
    fi && \
    echo Build complete. 

#Web Server
FROM scratch
LABEL app="InterSpec"
COPY --from=build_interspec /release /interspec/
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