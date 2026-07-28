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

# KNOWN PRE-EXISTING FAILURE (not caused by the native HTTP transport):
#
# The final link of InterSpec fails here with ~22 "undefined reference to boost::thread::..."
# from Wt's SocketNotifier.C.  Building this same image with --build-arg native_http=OFF fails
# identically, so it is not attributable to USE_NATIVE_HTTP_CLIENT.
#
# Cause: CMakeLists.txt has two byte-identical target_link_libraries() blocks for InterSpecLib
# (currently ~836-857 and ~933-953), both listing Boost::thread.  CMake de-duplicates them and
# keeps the first, which sits *before* the Wt libraries are added (~928).  GNU ld resolves static
# archives in one pass, so Wt's references to boost::thread go unresolved.  macOS does not hit
# this because ld64 is not order-sensitive in the same way.
#
# Everything up to the link is verified working on Linux: the static OpenSSL builds, Wt builds
# with SSL support (the WT_WITH_SSL assertion above passes), and both NativeHttpClient.cpp and
# NativeHttpClient_wt.cpp compile.  Fixing the link ordering is a separate change.

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
CMD ["./InterSpec", "--docroot", ".", "--http-address", "0.0.0.0", "--http-port", "8078", \
     "--userdatadir", "/data", "-c", "./data/config/wt_config_localweb.xml"]
