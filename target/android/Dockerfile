FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# System dependencies.
#
# On arm64 hosts we still need x86_64 glibc available so that the Android
# build-tools (aapt2 etc.), which Google ships as linux-x86_64 binaries
# only, can be transparently run by Rosetta/qemu/OrbStack.  Without
# libc6:amd64 the dynamic loader /lib64/ld-linux-x86-64.so.2 is missing
# and Gradle's aapt2 daemon fails with "Dynamic loader not found".
RUN dpkg --add-architecture amd64 \
    # arm64 Ubuntu's default sources point at ports.ubuntu.com, which
    # only carries arm64/ppc64el/s390x.  Rewrite Ubuntu 24.04's deb822
    # sources file to constrain ports.ubuntu.com to arm64 and add
    # archive.ubuntu.com for amd64/i386 so libc6:amd64 (the dynamic
    # loader Rosetta/OrbStack needs to run x86_64 binaries) can install.
    && rm -f /etc/apt/sources.list.d/ubuntu.sources \
    && printf 'Types: deb\nURIs: http://ports.ubuntu.com/ubuntu-ports/\nSuites: noble noble-updates noble-security noble-backports\nComponents: main restricted universe multiverse\nArchitectures: arm64\nSigned-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg\n\nTypes: deb\nURIs: http://archive.ubuntu.com/ubuntu/\nSuites: noble noble-updates noble-security\nComponents: main universe\nArchitectures: amd64 i386\nSigned-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg\n' > /etc/apt/sources.list.d/ubuntu.sources \
    && apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jdk-headless \
    git \
    cmake \
    ninja-build \
    zip \
    unzip \
    wget \
    curl \
    ca-certificates \
    libc6:amd64 \
    libstdc++6:amd64 \
    zlib1g:amd64 \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/lib/jvm/java-17-openjdk-* /usr/lib/jvm/java-17-openjdk

# Android SDK + NDK
ENV ANDROID_HOME=/opt/android-sdk
ENV PATH="${ANDROID_HOME}/cmdline-tools/latest/bin:${ANDROID_HOME}/platform-tools:${PATH}"

RUN mkdir -p "${ANDROID_HOME}/cmdline-tools" && \
    wget -q https://dl.google.com/android/repository/commandlinetools-linux-11076708_latest.zip \
         -O /tmp/cmdtools.zip && \
    unzip -q /tmp/cmdtools.zip -d "${ANDROID_HOME}/cmdline-tools" && \
    mv "${ANDROID_HOME}/cmdline-tools/cmdline-tools" "${ANDROID_HOME}/cmdline-tools/latest" && \
    rm /tmp/cmdtools.zip && \
    yes | sdkmanager --licenses > /dev/null 2>&1 && \
    sdkmanager \
        "platform-tools" \
        "platforms;android-35" \
        "build-tools;35.0.0" \
        "ndk;28.0.13004108"

# openjdk-17-jdk-headless installs to /usr/lib/jvm/java-17-openjdk-<arch>
# (amd64 or arm64). The symlink above resolves both to a stable path.
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk

WORKDIR /work
