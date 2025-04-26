# Build Deployable Image:
#  docker build --tag 'interspec_alpine' -f debian_web_container.dockerfile .

# Run Prebuilt Image:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078 interspec_debian

# Run From Dockerfile Directly:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078 -f alpine_web_container.dockerfile .

#  Optional: reuse existing build and/or source directories: Ensure ./build or ./src directory are populated and uncomment the COPY line

# If build fails, add the --no-cache option to the docker build command to force a fresh build
#  docker build --no-cache --build-arg true --tag 'interspec_debian' -f target/docker/debian_web_container.dockerfile .

FROM debian:stable-slim AS build
ARG repo=https://github.com/sandialabs/InterSpec.git
ARG tag=master
ARG savebuild=false
WORKDIR /work

# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build 

# RUN statements are broken up to allow loading cached images for debugging
RUN if [[ ${savebuild} != "false"  ]]; then --mount type=bind,source="${PWD}",target=/work,rw; fi &&\
apt update && apt upgrade -y && apt install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        cmake \
        patch \
        libsuitesparse-dev \
        libjs-uglify-js \
        zlib1g-dev \
        libssl-dev \
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
    if [ ! -d ./src ]; then \
        git clone --recursive --branch $BRANCH --depth=1 $REPO ./src; \
    fi
RUN cmake \
    -B ./build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_FOR_CONTAINER_LIBC=ON \
        -S ./src
RUN cmake --build build -j4; \
    mkdir -p release && \
    cmake --build build -j4 && \
    cmake --install ./build --prefix /release && \
    rm -rf release/lib/cmake && \
    cd /release && \ 
    chmod -R a+r * && \
    chmod a+x bin/InterSpec && \
    mkdir -p data && \
    chmod 777 data && \
    if [[ ${savebuild} != "false"  ]]; then \
        cp -r /release /work/release; \
    fi && \
    echo Build complete. 
#Web Server
FROM debian:stable-slim
LABEL app="InterSpec_Debian"
COPY --from=build /release /interspec/
WORKDIR /interspec
EXPOSE 8078
RUN apt update && apt upgrade -y && apt install -y --no-install-recommends \
        libblas3 \
        libstdc++6 \
        libc6 \
        libspqr2 && \
        mkdir -p /data && \
        rm -rf /var/cache/apt/archives /var/lib/apt/lists/*
CMD ["/interspec/bin/InterSpec", "--config", "/interspec/share/interspec/data/config/wt_config_web.xml", "--userdatadir=/data", "--http-port=8078", "--http-address=0.0.0.0", "--docroot", "/interspec/share/interspec"]


# Then numeric group/user value of 280 was chosen randomly; it doesnt conflict with existing groups/users on dev or public server, and is below 1000 (e.g., a system user without a home directory or default shell)
#RUN groupadd --gid 280 interspec && useradd --uid 280 --gid interspec interspec
#RUN addgroup -S interspec && adduser --disabled-password --no-create-home -S interspec -G interspec
#USER interspec
# Or just use user guest
#USER guest
# You could keep the access log by chenging the entrypoint to: "--accesslog=/mnt/interspec_data/wt_access_log.txt"
# You could also edit the <log-file></log-file> element of data/config/wt_config_web.xml to save the stdout/stderr of InterSpec to a log file at /mnt/interspec_data/interspec_log.txt.