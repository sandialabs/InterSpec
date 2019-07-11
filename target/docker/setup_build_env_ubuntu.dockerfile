#Dockerfile to compile InterSpec on (installs correct versions of Wt and Boost)

FROM ubuntu:18.04

# RUN export http_proxy=http://wwwproxy.ca.sandia.gov:80/
# RUN export https_proxy=https://wwwproxy.ca.sandia.gov:80/

ADD . /InterSpec

WORKDIR /tmp


#since this is just to build things, we'll use lots of run commands to take advantage of
# docker caching
RUN apt update && apt -y install bash build-essential binutils gcc cmake wget zlib1g-dev

#To Try and link to static runtime, add 'runtime-link=static' to below b2 command
RUN mkdir build; cd build \
&& wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz \
&& tar -xzvf boost_1_65_1.tar.gz \
&& cd boost_1_65_1 \
&& ./bootstrap.sh --without-icu --without-libraries=python,mpi,container,type_erasure,graph,graph_parallel,log,coroutine,stacktrace,test,wave --prefix=/usr/local \
&& ./b2 variant=release link=static cxxflags="-std=c++11" threading=multi -j4 install  \
&& cd .. \
&& rm -rf boost_1_65_1 boost_1_65_1.tar.gz

#if you want to try and static link, add '-DCMAKE_EXE_LINKER_FLAGS="-static"' to the below cmake command
RUN wget https://github.com/emweb/wt/archive/3.3.4.tar.gz \
&& tar -xzvf 3.3.4.tar.gz && cd wt-3.3.4 \
&& cp src/Wt/Dbo/backend/amalgamation/sqlite3.h /usr/local/include/ \
&& patch src/Wt/Render/CssParser.C /InterSpec/target/patches/wt/3.3.4/CssParser.C.patch \
&& patch src/http/RequestParser.C /InterSpec/target/patches/wt/3.3.4/RequestParser.C.patch \
&& mkdir build && cd build \
&& cmake -DWT_CPP_11_MODE="-std=c++11" -DWT_SIGNALS_IMPLEMENTATION="boost.signals2" \
      -DENABLE_SSL=OFF -DSHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release \
      -DCONNECTOR_FCGI=OFF -DCONNECTOR_HTTP=ON -DBUILD_EXAMPLES=OFF \
      -DINSTALL_RESOURCES=ON -DENABLE_LIBWTTEST=ON -DENABLE_LIBWTDBO=ON \
      -DENABLE_QT4=OFF -DENABLE_SQLITE=ON -DENABLE_EXT=OFF -DENABLE_OPENGL=ON \
      -DBOOST_PREFIX=/usr/local \
      .. \
&& make -j4 install \
&& cd ../../.. \
&& rm -rf build

# && patch -b src/Wt/Dbo/ptr_impl.h /InterSpec/target/patches/wt/3.3.4/ptr_impl.h.patch \
#Other Wt options we'll accept defaults on: ENABLE_HARU(ON), ENABLE_PANGO(ON), ENABLE_POSTGRES(ON), 
#   ENABLE_FIREBIRD(ON), ENABLE_MYSQL(ON), WT_NO_STD_LOCALE(OFF), WT_NO_STD_WSTRING(OFF)


#RUN mkdir build_interspec \
#&& cd build_interspec \
#&& cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_FOR_WEB_DEPLOYMENT=ON -DENABLE_RESOURCE_UPDATES=OFF \
#         -DBUILD_AS_LOCAL_SERVER=OFF -DUSE_SPECRUM_FILE_QUERY_WIDGET=OFF -DTRY_TO_STATIC_LINK=ON \
####         -DCMAKE_EXE_LINKER_FLAGS="-static -static-libstdc++ -static-libgcc"
#         /interspec \
#&& make -j4
#Doing the static linking causes the EXE to crash immediately.  Maybe this is to to PIE stuff?
#cd .. && mkdir InterSpec && cd InterSpec \
#&& cp ../build_interspec/bin/InterSpec.exe . \
#&& cp -r ../build_interspec/resources . \
#&& cp -r /interspec/InterSpec_resources . \
#&& cp -r /interspec/example_spectra . \
#&& cp -r /interspec/data .
#

#CMD ["/bin/bash"]
