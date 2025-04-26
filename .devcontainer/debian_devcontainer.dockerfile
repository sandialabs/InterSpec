FROM debian:stable-slim AS build
ARG repo=https://github.com/sandialabs/InterSpec.git
ARG tag=master
ARG savebuild=false
WORKDIR /work

# Optional uncomment and populate directories
# COPY src ./src
# COPY build ./build 

# RUN statements are broken up to allow loading cached images for debugging
RUN --mount type=bind,source=data,target=/data ls -alh /data