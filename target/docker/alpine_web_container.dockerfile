# Then numeric group/user value of 280 was chosen randomly; it doesnt conflict with existing groups/users on dev or public server, and is below 1000 (e.g., a system user without a home directory or default shell)
#RUN groupadd --gid 280 interspec && useradd --uid 280 --gid interspec interspec
#RUN addgroup -S interspec && adduser --disabled-password --no-create-home -S interspec -G interspec
#USER interspec
# Or just use user guest
#USER guest

# If build fails, issue this command
# docker rmi $(docker images -f “dangling=true” -q)
# docker system prune
# You could keep the access log by chenging the entrypoint to: "--accesslog=/mnt/interspec_data/wt_access_log.txt"
# You could also edit the <log-file></log-file> element of data/config/wt_config_web.xml to save the stdout/stderr of InterSpec to a log file at /mnt/interspec_data/interspec_log.txt.

#docker run --name=interspec
#Builder
FROM alpine:latest AS build
WORKDIR /work
USER guest
ENV tag=v1.0.13
RUN apk add --no-cache \
    alpine-sdk \
    cmake \
    patch \
    linux-headers \
    suitesparse-dev patch \
    curl && \
    curl --proto '=https' --tlsv1.2 -sSf https://raw.githubusercontent.com/devmatteini/dra/refs/heads/main/install.sh | bash -s -- && \
    ./dra download -s InterSpec-{tag}-source-code.tar.gz sandialabs/InterSpec && \
    mkdir src build && \
    tar -xf *.tar.gz -C ./src/ --strip-components=1 && \
    cmake \
        -B ./build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FOR_WEB_DEPLOYMENT=ON \
        -DUSE_REL_ACT_TOOL=ON \
        -DBUILD_AS_LOCAL_SERVER=OFF \
        -DInterSpec_FETCH_DEPENDENCIES=ON \
        -DBoost_INCLUDE_DIR=./build/_deps/boost-src/libs \
        -DUSE_SEARCH_MODE_3D_CHART=ON \
        -DUSE_QR_CODES=ON \
        -DUSE_DETECTION_LIMIT_TOOL=ON \
        -DUSE_SPECRUM_FILE_QUERY_WIDGET=ON \
        -DSpecUtils_PYTHON_BINDINGS=ON  \
        ./src && \
        mkdir -p /InterSpec && \
        cmake --build ./build --target install --config Release --prefix ./InterSpec ./src && \
#Web Server
FROM alpine:latest
LABEL app="InterSpec"
COPY --from=build /InterSpec/bin /InterSpec/lib /InterSpec/example_spectra /InterSpec/Interspec_resources /InterSpec/resources /InterSpec/data /var/opt/interspec/
WORKDIR /var/opt/interspec
USER guest
EXPOSE 8078
SHELL ["/bin/sh", "-c"]
RUN chmod -R a+r /var/opt/interspec \
     && chmod a+x /var/opt/interspec/bin/InterSpec \
     && chmod -R uga-w /var/opt/interspec
ENTRYPOINT ["./bin/InterSpec", "--userdatadir ../userdata", "-c ../wt_config_web.xml", "--http-port=8078", "--config=/var/opt/interspec/html_root/data/config/wt_config_web.xml","--static-data-dir", "/var/opt/interspec/html_root/data/"]