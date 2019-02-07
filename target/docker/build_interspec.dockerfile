FROM interspec_dev_env:1.0

ADD . /interspec

#TODO: We are currently using the sqlite3 database to store things, we should setup a 
#      container to run mysql too.

RUN mkdir build_interspec \
&& cd build_interspec \
&& cmake -DCMAKE_BUILD_TYPE="Debug" -DBUILD_FOR_WEB_DEPLOYMENT=ON -DENABLE_RESOURCE_UPDATES=OFF \
         -DBUILD_AS_LOCAL_SERVER=OFF -DUSE_SPECRUM_FILE_QUERY_WIDGET=OFF -DTRY_TO_STATIC_LINK=ON \
         -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
         /interspec \
&& make -j4

####         -DCMAKE_EXE_LINKER_FLAGS="-static -static-libstdc++ -static-libgcc" \
#Adding the "-static" flag to the above causes the EXE to crash immediately.  Maybe this is to to PIE stuff?  GDB didnt help in figuring it out

RUN cd .. && mkdir InterSpec && cd InterSpec \
&& cp ../build_interspec/bin/InterSpec.exe . \
&& cp -r ../build_interspec/resources . \
&& cp -r /interspec/InterSpec_resources . \
&& cp -r /interspec/example_spectra . \
&& cp -r /interspec/data .