#Dockerfile to compile InterSpec on (installs correct versions of Wt and Boost)

FROM alpine:3.9 as InterSpecBuilder

#ENV http_proxy wwwproxy.ca.sandia.gov:80
#ENV https_proxy wwwproxy.ca.sandia.gov:80


ADD . /InterSpec

WORKDIR /tmp

#since this is just to build things, we'll use lots of run commands to d=take advantage of
# docker caching

#gcc version 8.2.0
#cmake version 3.13.0
RUN apk add build-base binutils gcc abuild cmake linux-headers
RUN mkdir build; cd build \
&& wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz \
&& tar -xzvf boost_1_65_1.tar.gz \
&& cd boost_1_65_1
RUN ./bootstrap.sh --without-icu --without-libraries=python,mpi,container,type_erasure,graph,graph_parallel,log,coroutine,stacktrace,test,wave --prefix=/usr/local
RUN ./b2 variant=release link=static threading=multi runtime-link=static -j8 install 
WORKDIR /tmp
RUN wget https://github.com/emweb/wt/archive/3.3.4.tar.gz \
&& tar -xzvf 3.3.4.tar.gz && cd wt-3.3.4
RUN patch -b src/Wt/Render/CssParser.C /InterSpec/target/patches/wt/3.3.4/CssParser.C.patch
RUN patch -b src/http/RequestParser.C /InterSpec/target/patches/wt/3.3.4/RequestParser.C.patch
RUN patch -b src/Wt/Dbo/ptr_impl.h /InterSpec/target/patches/wt/3.3.4/ptr_impl.h.patch
RUN mkdir build && cd build
RUN cmake -DWT_CPP_11_MODE="-std=c++11" -DWT_SIGNALS_IMPLEMENTATION="boost.signals2" \
      -DENABLE_SSL=OFF -DSHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release \
      -DCONNECTOR_FCGI=OFF -DCONNECTOR_HTTP=ON -DBUILD_EXAMPLES=OFF \
      -DINSTALL_RESOURCES=ON -DENABLE_LIBWTTEST=ON -DENABLE_LIBWTDBO=ON \
      -DENABLE_QT4=OFF -DENABLE_SQLITE=ON -DENABLE_EXT=OFF -DENABLE_OPENGL=ON \
      -DBOOST_PREFIX=/usr/local \
      ..
RUN make -j8 install
WORKDIR /tmp
#RUN rm -rf build

#Other Wt options we'll accept defaults on: ENABLE_HARU(ON), ENABLE_PANGO(ON), ENABLE_POSTGRES(ON), 
#   ENABLE_FIREBIRD(ON), ENABLE_MYSQL(ON), WT_NO_STD_LOCALE(OFF), WT_NO_STD_WSTRING(OFF)


