FROM mcr.microsoft.com/devcontainers/base:alpine-3.21
LABEL app="interspec_alpine_base"
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