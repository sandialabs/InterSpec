
# Step 1) Create base image, with boost 1.65.1 and Wt 3.3.4 installed
docker build -t interspec_ubuntu_build_env:v1.0 -f "setup_build_env_ubuntu.dockerfile" ../..

# Step 2) Create and image with a compiled InterSpec 
docker build -t interspec_ubuntu_build:v1.0 -f "build_interspec_ubuntu.dockerfile" ../..

# Step 3) Enter image with compiled InterSpec.  This next line allows you to use gdb, and it maps the images port 8080 to your local machines port 8080
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -p 8080:8080 -i -t -v /Users/wcjohns/rad_ana/InterSpec:/interspec interspec_ubuntu_build:v1.0 bash

# Step 4) From within docker now
cd /tmp/build_interspec
./bin/InterSpec.exe --docroot . --http-address 0.0.0.0 --http-port 8080 -c ./data/config/wt_config_localweb.xml

# Step 5) from your normal browser, visit localhost:8080 to use InterSpec