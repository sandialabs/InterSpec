#!/bin/sh

rm -f app/src/main/assets/interspec-assets.zip
cd ../../../ 
zip -9 -r --exclude=*OUO* --exclude=*ouo* --exclude=*.svn* --exclude=*.git* --exclude=*.DS_Store target/android/InterSpec/app/src/main/assets/interspec-assets.zip InterSpec_resources data example_spectra
cd build_xcode
zip -ur -9 -r --exclude=*OUO* --exclude=*ouo* --exclude=*.svn* --exclude=*.git* --exclude=*.DS_Store ../target/android/InterSpec/app/src/main/assets/interspec-assets.zip resources
cd ../target/android/InterSpec/
