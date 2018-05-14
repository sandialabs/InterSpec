echo "running buildscript.sh"
source ~/.profile
cd ../../../build_ios
make -j8
cd ../target/ios/InterSpec
