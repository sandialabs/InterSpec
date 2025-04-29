FROM debian:stable-slim AS build
ARG repo=https://github.com/sandialabs/InterSpec.git
ARG tag=master
ARG useHost=true
WORKDIR /work
# COPY --exclude src/build --exclude src* /work/src
# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build 

# RUN statements are broken up to allow loading cached images for debugging
RUN --mount=type=bind,target=/mnt/host,rw apt update && apt upgrade -y && apt install -y --no-install-recommends \
    ca-certificates \
    build-essential \
    libhpdf-dev \
    libfcgi-dev \
    freeglut3-dev \
    libpango1.0-dev \
    libgraphicsmagick++1-dev \
    libglew-dev \
    cmake \
    patch \
    libsuitesparse-dev \
    libjs-uglify-js \
    zlib1g-dev \
    libssl-dev \
    git && \  
    if [ "$useHost" = true ]; then \
        srcdir=/mnt/host/src && \
        builddir=/mnt/host/build && \
        mkdir -p mnt/host/src && \
        mkdir -p mnt/host/build; \
    else \
        mkdir -p /work/src && \
        mkdir -p /work/build && \
        srcdir=/work/src && \
        builddir=/work/build; \
    fi && \
    if [ ! -d "$srcdir/InterSpec_resources" ]; then \ 
        git clone --recursive --branch $tag --depth=1 $repo $srcdir; \
    fi && \
    cmake \
    -B "$builddir" \
    -DBUILD_FOR_CONTAINER_LIBC=ON \
    -DInterSpec_FETCH_DEPENDENCIES=ON \
    -S "$srcdir" && \
    cmake -B "$builddir" --build build -j4 && \
    cmake --install "$builddir" --prefix /release

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
