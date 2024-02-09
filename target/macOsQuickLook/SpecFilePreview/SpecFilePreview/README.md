As of 20240208, previewing peaks exported from InterSpec isnt working.


You can build and test this QuickLook utility on its own. e.g.,

```bash
# You need to edit CMakeLists.txt to uncomment out `set( InterSpec_BASE_DIR ../../../../ )`

cd InterSpec/target/macOsQuickLook/SpecFilePreview/SpecFilePreview/
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/Users/wcjohns/install/macOS_wt3.7.1_prefix ..
make

rm -rf ./SpecFilePreview.qlgenerator
qlmanage -r
qlmanage -g ./SpecFilePreview.qlgenerator -c "gov.sandia.InterSpec.gamma-spectrum" -p '/path/to/somefile.n42'
make
qlmanage -g ./SpecFilePreview.qlgenerator -c "gov.sandia.InterSpec.gamma-spectrum" -p '/path/to/somefile.n42'
```

But note that the `qlmanage -r` and dummy call to qlmanage might be necessary, as it seems sticky for some reason.