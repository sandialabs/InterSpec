# Builds InterSpec on Linux with the native HTTP transport, exercising the dependency-prefix
# path (target/patches/dep_build_linux.sh) with a static OpenSSL, so that Wt::Http::Client can
# actually do https.
#
# This exists to test the USE_NATIVE_HTTP_CLIENT Linux backend from a macOS development machine.
# The other containerfiles here use CMake FetchContent for dependencies and leave Wt's SSL off;
# this one deliberately uses the prefix scripts, because those are what a Linux release build
# uses and what the BUILD_OPENSSL flag was added to.
#
# Build (from the repository root):
#   docker build -f target/docker/linux_native_http_test.containerfile -t interspec_native_http .
#
# Run, then point a browser on the host at http://localhost:8078
#   docker run --rm -p 8078:8078 interspec_native_http
#
# To reach a proxy or an LLM endpoint running on the macOS host, use host.docker.internal.
# To supply a real llm_config.xml:
#   docker run --rm -p 8078:8078 -v /path/to/userdata:/data interspec_native_http

FROM debian:trixie-slim AS build

# `ca-certificates` matters for more than apt here: it is what populates /etc/ssl/certs, which is
#  the trust store the Wt backend probes at runtime.  A slim image without it would make every
#  https request fail certificate verification.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        cmake \
        ninja-build \
        curl \
        git \
        patch \
        perl \
        unzip \
        socat \
        libsuitesparse-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY . /work/InterSpec

# Build the dependency prefix, with a static OpenSSL and Wt SSL support turned on.  Everything
#  lands in /work/prefix; nothing links against the distribution's OpenSSL.
RUN _ncore=$(nproc) BUILD_OPENSSL=ON \
      bash /work/InterSpec/target/patches/dep_build_linux.sh \
        /work/InterSpec /work/deps /work/prefix

# Fail early and loudly if Wt did not actually come out with SSL support - otherwise the failure
#  would not surface until a request was made at runtime.
RUN grep -q "^#define WT_WITH_SSL" /work/prefix/include/Wt/WConfig.h \
    && echo "Wt built WITH SSL, as intended"

# Overridable so the same image can be built without the native transport, which is how you
# check whether a build failure is attributable to it.
ARG native_http=ON

# Note: this build is also what surfaced a pre-existing Linux link failure - ~22 "undefined
# reference to boost::thread::..." out of Wt's SocketNotifier.C, reproducible with
# --build-arg native_http=OFF, so not attributable to the native transport.  CMake emitted the
# Boost archives only *ahead* of the Wt ones, and GNU ld scans each archive once in order.  Fixed
# by repeating the Boost archives after Wt in CMakeLists.txt; see the comment there.

RUN cmake -G Ninja -S /work/InterSpec -B /work/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=/work/prefix \
      -DBOOST_ROOT=/work/prefix \
      -DUSE_NATIVE_HTTP_CLIENT=${native_http} \
      -DUSE_LLM_INTERFACE=ON \
      -DUSE_BATCH_CLI_TOOLS=OFF \
      -DUSE_BATCH_GUI_TOOLS=OFF \
    && cmake --build /work/build -j$(nproc)

# Confirm the binary is not linked against the distribution's OpenSSL - the whole point of
#  building it into the prefix statically.
RUN ldd /work/build/InterSpec | grep -E "libssl|libcrypto" \
    && (echo "ERROR: linked against a shared OpenSSL" && exit 1) \
    || echo "No shared OpenSSL linkage, as intended"

EXPOSE 8078
VOLUME /data

WORKDIR /work/build
# InterSpec only offers --http-address in BUILD_FOR_WEB_DEPLOYMENT builds, and that option is
# mutually exclusive with USE_LLM_INTERFACE (CMakeLists FATAL_ERRORs on the combination).  So the
# server necessarily binds container-localhost, which the host cannot reach - forward it.
CMD ["/bin/sh", "-c", \
     "./InterSpec --docroot . --http-port 8079 --userdatadir /data --config ./data/config/wt_config_localweb.xml & \
      exec socat TCP-LISTEN:8078,fork,reuseaddr TCP:127.0.0.1:8079"]
