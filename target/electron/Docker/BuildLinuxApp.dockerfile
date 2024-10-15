FROM quay.io/pypa/manylinux2014_x86_64:latest

WORKDIR /mnt/code/target/electron

RUN CMAKE_BUILD_PARALLEL_LEVEL=4 cmake-js --architecture x64 --arch=x64 --CDCMAKE_BUILD_TYPE="Release" --CDInterSpec_FETCH_DEPENDENCIES=ON --CDBUILD_AS_LOCAL_SERVER=OFF --CDCMAKE_SHARED_LINKER_FLAGS="-static-libgcc -static-libstdc++" --CDUSE_LEAFLET_MAP=ON --CDLEAFLET_MAPS_KEY="..." --CDUSE_REL_ACT_TOOL=ON --out=build_manylinux_electron --target install \
    && npm run package-manylinux \
    && cp ../../NOTICE.html ./release-builds/InterSpec-linux-x64/