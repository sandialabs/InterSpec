# Build Deployable Image:
#   docker build --tag 'interspec_alpine' -f alpine_web_container.containerfile .

# Run Prebuilt Image:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078/tcp interspec_alpine
# Or if you dont care about keeping user preferences and stuff around, you can just just map /data to a temp ephemeral dir
#   docker run --rf -v /data -p 8078:8078/tcp interspec_alpine

# Run From Dockerfile Directly:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078 -f alpine_web_container.dockerfile .

#  Optional: reuse existing source, dependency, build directories:
    # To reuse sources, use --build-args hostdepdir=true, place the sources in the build/deps directory and optionally the Interspec sources in ./ .
    # To reuse the complete build directory set --build-args hostbuilddir=true.

# ex: docker build --build-arg hostdepdir=true --build-arg hostbuilddir=true -t interspec -f alpine_web_container.containerfile .
# When hostdepdir=true, the container will reuse everything in ./build/_deps/*-src
# When hostbuilddir=true, the container will reuse everything in the ./build directory and will look for the source repo in the ./ directory.

# If build fails, add the --no-cache option to the docker build command to force a fresh build
#  docker build --no-cache --tag 'interspec_alpine' -f alpine_web_container.dockerfile .

FROM alpine:3 AS build_interspec
ARG repo=https://github.com/sandialabs/InterSpec.git \
ARG tag=master \
ARG hostdepdir=true \
ARG hostbuilddir=false \
WORKDIR /work

# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build

RUN --mount=type=bind,target=/mnt/host,rw    \
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
    if [ "$hostbuilddir" = "true" || "$hostdepdir" = "true" ]; then \
        echo "Using host build directory" && \
        srcdir=/mnt/host && \
        builddir=/mnt/host/build && \
        if [ "$hostdepdir" = true ] && [ "$hostbuilddir" != true ]; then \
            echo "Using host dependency directory" && \
            mkdir -p /mnt/host/build/_deps && \
            mkdir -p /mnt/host/_deps && \
            ls -al && \
            mv /mnt/host/build/_deps/*-src /mnt/host/_deps/ 2>/dev/null && \
            rm -rf /mnt/host/build && \
            ls -al /mnt/host/_deps && \
            mkdir -p /mnt/host/build && \
            mv /mnt/host/_deps /mnt/host/build/; \
        fi; \
    else \
        mkdir -p /work/src && \
        mkdir -p /work/build && \
        srcdir=/work/src && \
        builddir=/work/build; \
    fi && \
    if [ ! -d "$srcdir/InterSpec_resources" ]; then \ 
        echo "source repo missing, cloning" && \
        git clone --recursive --branch $tag --depth=1 $repo $srcdir; \
    fi && \
    cmake \
        -B "$builddir" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FOR_CONTAINER_MUSL=ON \
        -S "$srcdir" && \
    cmake --build "$builddir" -j4 && \
    mkdir -p release && \
    cmake --build "$builddir" -j4 && \
    cmake --install "$builddir" --prefix /release && \
    rm -rf release/lib/cmake && \
    cd /release && \ 
    chmod -R a+r * && \
    chmod a+x bin/InterSpec && \
    mkdir -p data && \
    chmod 777 data && \
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