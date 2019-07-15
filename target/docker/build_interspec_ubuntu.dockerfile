FROM interspec_ubuntu_build_env:v1.0

ADD . /interspec

#TODO: We are currently using the sqlite3 database to store things, we should setup a 
#      container to run mysql too.

RUN cd /tmp && mkdir build_interspec && cd build_interspec \
&& cmake -DCMAKE_BUILD_TYPE="Debug" -DBUILD_FOR_WEB_DEPLOYMENT=ON -DENABLE_RESOURCE_UPDATES=OFF \
         -DBUILD_AS_LOCAL_SERVER=OFF -DUSE_SPECRUM_FILE_QUERY_WIDGET=OFF -DTRY_TO_STATIC_LINK=ON \
         /interspec \
&& make -j4

# This next step is just what you would want to do to make a stand-alone release of InterSpec.
# If you are going to use this, you probably want to change the "Debug" above, to "Release"
RUN mkdir -p /tmp/InterSpec && cd /tmp/InterSpec \
&& cp /tmp/build_interspec/bin/InterSpec.exe . \
&& cp -r /tmp/build_interspec/resources . \
&& cp -r /interspec/InterSpec_resources . \
&& cp -r /interspec/example_spectra . \
&& cp -r /interspec/data .