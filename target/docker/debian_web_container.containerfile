# 
# Build Deployable Image:
#   docker build --tag 'interspec_debian' -f debian_web_container.containerfile .

# Run Prebuilt Image:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078/tcp interspec_debian
# Or if you dont care about keeping user preferences and stuff around, you can just just map /data to a temp ephemeral dir
#   docker run --rf -v /data -p 8078:8078/tcp interspec_debian

# Run From Dockerfile Directly:
#   docker run --rm -v path/to/your/datadir:/data -p 8078:8078 -f debian_web_container.dockerfile .

#  Optional: reuse existing source, dependency, build directories:
    # To reuse sources, use --build-args hostdepdir=true, place the sources in the build/deps directory and optionally the Interspec sources in ./ .
    # To reuse the complete build directory set --build-args hostbuilddir=true.

# ex: docker build --build-arg hostdepdir=true --build-arg hostbuilddir=true -t interspec -f debian_web_container.containerfile .
# When hostdepdir=true, the container will reuse everything in ./build/_deps/*-src
# When hostbuilddir=true, the container will reuse everything in the ./build directory and will look for the source repo in the ./ directory.

# If build fails, add the --no-cache option to the docker build command to force a fresh build
#  docker build --no-cache --tag 'interspec_debian' -f debian_web_container.dockerfile .

FROM debian:trixie-slim AS build
ARG repo=https://github.com/sandialabs/InterSpec.git \
ARG tag=master \
ARG hostdepdir=true \
ARG hostbuilddir=true \
WORKDIR /work
#libhpdf-dev libfcgi-dev freeglut3-dev libpango1.0-dev libglew-dev libopenblas-dev liblapack-dev libgmp-dev libmpfr-dev
# RUN statements are broken up to allow loading cached images for debugging
RUN --mount=type=bind,target=/mnt/host,rw    \
    apt update && apt upgrade -y && apt install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        libsuitesparse-dev \
        libjs-uglify-js \
        zlib1g-dev \
        libssl-dev \
        cmake \
        patch \
        git  && \
    if [ "$hostbuilddir" = true ] || [ "$hostdepdir" = true ]; then \
        echo "Using host build directory" && \
        srcdir=/mnt/host && \
        builddir=/mnt/host/build && \
        if [ "$hostdepdir" = true ] && [ "$hostbuilddir" != true ]; then \
            echo "Using host dependency directory" && \
            mkdir -p /mnt/host/build/_deps && \
            mkdir -p /mnt/host/_deps && \
            mv /mnt/host/build/_deps/*-src /mnt/host/_deps/ 2>/dev/null && \
            rm -rf /mnt/host/build && \
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
    -DCONTAINER_LIBC=ON \
    -DInterSpec_FETCH_DEPENDENCIES=ON \
    -S "$srcdir" && \
    cmake --build "$builddir" -j4 && \
    cmake --install "$builddir" --prefix /release && \
    cd /release && \ 
    chmod -R a+r * && \
    chmod a+x bin/InterSpec && \
    mkdir -p data && \
    chmod 777 data && \
    echo Build complete. 

FROM debian:trixie-slim
LABEL app="InterSpec_Debian"
COPY --from=build /release /interspec/
WORKDIR /interspec
EXPOSE 8078
RUN apt update && apt upgrade -y && apt install -y --no-install-recommends \
        libblas3 \
        libstdc++6 \
        libc6 \
        libsuitesparse-dev && \
        mkdir -p /data && \
        rm -rf /var/cache/apt/archives /var/lib/apt/lists/*
CMD ["/interspec/bin/InterSpec", "--config", "/interspec/share/interspec/data/config/wt_config_web.xml", "--userdatadir=/data", "--http-port=8078", "--http-address=0.0.0.0", "--docroot", "/interspec/share/interspec"]
