#
# Builds Wt as a framework for iPhone/iPad devices.
# Copyright (C) 2011 Emweb bvba, Kessel-Lo, Belgium.
#
# Inspired by Pete Goodliffe's boost.sh script
#

: ${BOOST_FRAMEWORK_PATH:=`pwd`/framework}
: ${BOOST_PREFIX:=`pwd`/../prefix}
: ${IOS_SDK:=/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer}
: ${ISIM_SDK:=/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneSimulator.platform/Developer}
: ${SDK_VER:=`xcodebuild -showsdks | grep iphoneos | egrep "[[:digit:]]+\.[[:digit:]]+" -o | tail -1`}}
: ${XCODE_ROOT:=`xcode-select -print-path`}
: ${COMPILER:="clang++"}

: ${CMAKE:=cmake}

: ${BUILD_DIR:=`pwd`}
: ${WT_DIR:=`pwd`/wt-3.3.4}
: ${WT_VERSION:=3_3_4}
: ${FRAMEWORKDIR:=`pwd`}

BUILD_ARMV7_DIR=$BUILD_DIR/build-armv7
BUILD_ARMV7S_DIR=$BUILD_DIR/build-armv7s
BUILD_ARM64_DIR=$BUILD_DIR/build-arm64
BUILD_I386_DIR=$BUILD_DIR/build-i386
BUILD_X86_64_DIR=$BUILD_DIR/build-x86-64


#STAGE_DIR=$BUILD_DIR/stage
STAGE_DIR=`pwd`/../prefix
TMP_DIR=$BUILD_DIR/tmp

# PH: added -DBOOST_SIGNALS_NO_DEPRECATION_WARNING because Boost.Signals is deprecitated in 1.54
: ${COMMON_CMAKE_FLAGS=-DUSE_BOOST_FRAMEWORK=ON \
       -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_FRAMEWORK_PATH:FILEPATH=$BOOST_FRAMEWORK_PATH \
       -DBOOST_PREFIX=$BOOST_PREFIX \
       -DCMAKE_INSTALL_PREFIX=$STAGE_DIR \
       -DSHARED_LIBS=OFF \
       -DCONNECTOR_FCGI=OFF \
       -DBUILD_TESTS=OFF \
       -DBUILD_EXAMPLES=OFF \
       -DENABLE_GM=OFF \
       -DENABLE_HARU=OFF \
       -DENABLE_QT4=OFF \
       -DENABLE_SSL=OFF \
       -DENABLE_POSTGRES=OFF \
       -DENABLE_FIREBIRD=OFF \
       -DHTTP_WITH_ZLIB=OFF \
       -DWT_CPP_11_MODE="-std=c++11" \
       -DENABLE_SSL=OFF \
       -DWT_BOOST_DISCOVERY=OFF \
       -DBoost_INCLUDE_DIR:FILEPATH=$BOOST_PREFIX/include
       -DWT_SIGNALS_IMPLEMENTATION="boost.signals2" }


abort()
{
    echo
    echo "Aborted: $@"
    exit 1
}

clean()
{
    rm -rf $BUILD_ARMV7_DIR
    rm -rf $BUILD_ARMV7S_DIR
    rm -rf $BUILD_ARM64_DIR
    rm -rf $BUILD_I386_DIR
    rm -rf $BUILD_X86_64_DIR
    rm -rf $TMP_DIR
}


build-armv7()
{
   echo '======================================================='
   echo '= Building armv7'
   echo '======================================================='
   sleep 2
   [ -d $BUILD_ARMV7_DIR ] || mkdir -p $BUILD_ARMV7_DIR
   ( cd $BUILD_ARMV7_DIR; ${CMAKE} \
       $COMMON_CMAKE_FLAGS \
       -DCMAKE_C_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang \
       -DCMAKE_CXX_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/$COMPILER \
       -DCMAKE_OSX_SYSROOT:FILEPATH=$IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk \
       -DCMAKE_C_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -fvisibility=hidden -arch armv7 -pipe" \
       -DCMAKE_CXX_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -DBOOST_SIGNALS_NO_DEPRECATION_WARNING -std=c++11 -stdlib=libc++ -fvisibility=hidden -fvisibility-inlines-hidden -isysroot $IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk -arch armv7 -pipe" \
       $WT_DIR )
       #cd /Users/wcjohns/rad_ana/InterSpec/target/ios/wt/wt-3.3.4
       cd $BUILD_ARMV7_DIR
       make -j8 || abort "Failed building for arm7 architecture"
}

build-armv7s()
{
   echo '======================================================='
   echo '= Building armv7s'
   echo '======================================================='
   sleep 2
   [ -d $BUILD_ARMV7S_DIR ] || mkdir -p $BUILD_ARMV7S_DIR
   ( cd $BUILD_ARMV7S_DIR; ${CMAKE} \
       $COMMON_CMAKE_FLAGS \
       -DCMAKE_C_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang \
       -DCMAKE_CXX_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/$COMPILER \
       -DCMAKE_OSX_SYSROOT:FILEPATH=$IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk \
       -DCMAKE_C_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -fvisibility=hidden -isysroot $IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk -arch armv7s -pipe" \
       -DCMAKE_CXX_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -DBOOST_SIGNALS_NO_DEPRECATION_WARNING -std=c++11 -stdlib=libc++ -fvisibility=hidden -fvisibility-inlines-hidden -isysroot $IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk -arch armv7s -pipe" \
       $WT_DIR )
       #cd /Users/ckmorte/InterSpec/target/ios/wt/wt-3.3.4
       cd $BUILD_ARMV7S_DIR
       make -j8 || abort "Failed building for arm7s architecture"
}


build-arm64()
{
   echo '======================================================='
   echo '= Building arm64'
   echo '======================================================='
   sleep 2
   [ -d $BUILD_ARM64_DIR ] || mkdir -p $BUILD_ARM64_DIR
   ( cd $BUILD_ARM64_DIR; ${CMAKE} \
       $COMMON_CMAKE_FLAGS \
       -DCMAKE_C_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang \
       -DCMAKE_CXX_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/$COMPILER \
       -DCMAKE_OSX_SYSROOT:FILEPATH=$IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk \
       -DCMAKE_C_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -fvisibility=hidden -isysroot $IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk -arch arm64 -pipe" \
       -DCMAKE_CXX_FLAGS:STRING="-DBOOST_AC_USE_PTHREADS -DBOOST_SP_USE_PTHREADS -DBOOST_SIGNALS_NO_DEPRECATION_WARNING -std=c++11 -stdlib=libc++ -fvisibility=hidden -fvisibility-inlines-hidden -isysroot $IOS_SDK/SDKs/iPhoneOS${SDK_VER}.sdk -arch arm64 -pipe" \
       $WT_DIR )
       #cd /Users/ckmorte/InterSpec/target/ios/wt/wt-3.3.4
       cd $BUILD_ARM64_DIR
       make -j8 || abort "Failed building for arm64 architecture"
}

build-i386()
{
   echo '======================================================='
   echo '= Building i386'
   echo '======================================================='
   sleep 2
   [ -d $BUILD_I386_DIR ] || mkdir -p $BUILD_I386_DIR
   ( cd $BUILD_I386_DIR; ${CMAKE} \
       $COMMON_CMAKE_FLAGS \
       -DCMAKE_C_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang \
       -DCMAKE_CXX_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/$COMPILER \
       -DCMAKE_OSX_SYSROOT:FILEPATH=$ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk \
       -DCMAKE_C_FLAGS:STRING="-arch i386 -fvisibility=hidden -miphoneos-version-min=7.0 -isysroot $ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk" \
       -DCMAKE_CXX_FLAGS:STRING="-arch i386 -fvisibility=hidden -miphoneos-version-min=7.0 -fvisibility-inlines-hidden -isysroot $ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk -DBOOST_SIGNALS_NO_DEPRECATION_WARNING" \
       $WT_DIR )
       #cd /Users/ckmorte/InterSpec/target/ios/wt/wt-3.3.4
       cd $BUILD_I386_DIR
       make -j8 || abort "Failed building for i386 architecture"
}


build-X86-64()
{
   echo '======================================================='
   echo '= Building x86_64'
   echo '======================================================='
   sleep 2
   [ -d $BUILD_X86_64_DIR ] || mkdir -p $BUILD_X86_64_DIR
   ( cd $BUILD_X86_64_DIR; ${CMAKE} \
       $COMMON_CMAKE_FLAGS \
       -DCMAKE_C_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang \
       -DCMAKE_CXX_COMPILER:FILEPATH=$XCODE_ROOT/Toolchains/XcodeDefault.xctoolchain/usr/bin/$COMPILER \
       -DCMAKE_OSX_SYSROOT:FILEPATH=$ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk \
       -DCMAKE_C_FLAGS:STRING="-arch x86_64 -fvisibility=hidden -miphoneos-version-min=7.0 -isysroot $ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk" \
       -DCMAKE_CXX_FLAGS:STRING="-arch x86_64 -fvisibility=hidden -miphoneos-version-min=7.0 -fvisibility-inlines-hidden -isysroot $ISIM_SDK/SDKs/iPhoneSimulator${SDK_VER}.sdk -DBOOST_SIGNALS_NO_DEPRECATION_WARNING" \
       $WT_DIR )
       #cd /Users/ckmorte/InterSpec/target/ios/wt/wt-3.3.4
       cd $BUILD_X86_64_DIR
       make VERBOSE=1 -j8 || abort "Failed building for x86_64 architecture"
}

combineLibs()
{
    : ${1:?}
    BUILD_DIR=$1
    [ -d $TMP_DIR ] || mkdir -p $TMP_DIR
    ( mkdir -p $TMP_DIR/http &&
	cd $TMP_DIR/http &&
	ar x $BUILD_DIR/src/http/libwthttp.a ) || abort "Extract libwthttp.a failed"
    ( mkdir -p $TMP_DIR/dbosqlite3 &&
	cd $TMP_DIR/dbosqlite3 &&
        ar x $BUILD_DIR/src/Wt/Dbo/backend/libwtdbosqlite3.a sqlite3.o &&
        mv sqlite3.o db_sqlite3.o &&
        ar x $BUILD_DIR/src/Wt/Dbo/backend/libwtdbosqlite3.a Sqlite3.o) || abort "Extract libwtdbosqlite3.a failed"
    ( mkdir -p $TMP_DIR/dbo &&
	cd $TMP_DIR/dbo &&
        ar x $BUILD_DIR/src/Wt/Dbo/libwtdbo.a) || abort "Extract libwtdbo.a failed"
    ( cd $TMP_DIR
	ar x $BUILD_DIR/src/libwt.a &&
	ar -rcs $BUILD_DIR/Wt.a *.o */*.o)  || abort "Combine libs $1 failed."
    rm -rf $TMP_DIR
}

                    VERSION_TYPE=Alpha
                  FRAMEWORK_NAME=Wt
               FRAMEWORK_VERSION=A

       FRAMEWORK_CURRENT_VERSION=$WT_VERSION
 FRAMEWORK_COMPATIBILITY_VERSION=$WT_VERSION

createFramework()
{
    FRAMEWORK_BUNDLE=$FRAMEWORKDIR/$FRAMEWORK_NAME.framework

    rm -rf $FRAMEWORK_BUNDLE

    mkdir -p $FRAMEWORK_BUNDLE
    mkdir -p $FRAMEWORK_BUNDLE/Versions
    mkdir -p $FRAMEWORK_BUNDLE/Versions/$FRAMEWORK_VERSION
    mkdir -p $FRAMEWORK_BUNDLE/Versions/$FRAMEWORK_VERSION/Resources
    mkdir -p $FRAMEWORK_BUNDLE/Versions/$FRAMEWORK_VERSION/Headers
    mkdir -p $FRAMEWORK_BUNDLE/Versions/$FRAMEWORK_VERSION/Documentation

    ln -s $FRAMEWORK_VERSION               $FRAMEWORK_BUNDLE/Versions/Current
    ln -s Versions/Current/Headers         $FRAMEWORK_BUNDLE/Headers
    ln -s Versions/Current/Resources       $FRAMEWORK_BUNDLE/Resources
    ln -s Versions/Current/Documentation   $FRAMEWORK_BUNDLE/Documentation
    ln -s Versions/Current/$FRAMEWORK_NAME $FRAMEWORK_BUNDLE/$FRAMEWORK_NAME

    FRAMEWORK_INSTALL_NAME=$FRAMEWORK_BUNDLE/Versions/$FRAMEWORK_VERSION/$FRAMEWORK_NAME

    echo "Framework: Creating final library..."
    xcrun -sdk iphoneos lipo -create \
	-arch armv7 build-armv7/Wt.a \
	-arch armv7s build-armv7s/Wt.a \
	-arch i386 build-i386/Wt.a \
	-o $FRAMEWORK_INSTALL_NAME || abort "lipo failed"

    echo "Framework: Copying headers..."
cp -r ${STAGE_DIR}/include/Wt/* $FRAMEWORK_BUNDLE/Headers/

    echo "Framework: Creating plist..."
    cat > $FRAMEWORK_BUNDLE/Resources/Info.plist <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
	<key>CFBundleDevelopmentRegion</key>
	<string>English</string>
	<key>CFBundleExecutable</key>
	<string>${FRAMEWORK_NAME}</string>
	<key>CFBundleIdentifier</key>
	<string>org.boost</string>
	<key>CFBundleInfoDictionaryVersion</key>
	<string>6.0</string>
	<key>CFBundlePackageType</key>
	<string>FMWK</string>
	<key>CFBundleSignature</key>
	<string>????</string>
	<key>CFBundleVersion</key>
	<string>${FRAMEWORK_CURRENT_VERSION}</string>
</dict>
</plist>
EOF
}

#clean
#build-armv7
#build-armv7s
#build-arm64
#build-i386
#build-X86-64
#combineLibs $BUILD_ARMV7_DIR
#combineLibs $BUILD_ARMV7S_DIR
#combineLibs $BUILD_ARM64_DIR
#combineLibs $BUILD_I386_DIR
#combineLibs $BUILD_X86_64_DIR
#createFramework

[ -d ${STAGE_DIR}/lib ] || mkdir -p ${STAGE_DIR}/lib

xcrun -sdk iphoneos lipo -create ./build-armv7/src/http/libwthttp.a ./build-armv7s/src/http/libwthttp.a ./build-arm64/src/http/libwthttp.a ./build-i386/src/http/libwthttp.a ./build-x86-64/src/http/libwthttp.a -output ${STAGE_DIR}/lib/libwthttp.a
xcrun -sdk iphoneos lipo -create ./build-armv7/src/Wt/Dbo/backend/libwtdbosqlite3.a ./build-armv7s/src/Wt/Dbo/backend/libwtdbosqlite3.a ./build-arm64/src/Wt/Dbo/backend/libwtdbosqlite3.a ./build-i386/src/Wt/Dbo/backend/libwtdbosqlite3.a ./build-x86-64/src/Wt/Dbo/backend/libwtdbosqlite3.a -output ${STAGE_DIR}/lib/libwtdbosqlite3.a
xcrun -sdk iphoneos lipo -create ./build-armv7/src/Wt/Dbo/libwtdbo.a ./build-armv7s/src/Wt/Dbo/libwtdbo.a ./build-arm64/src/Wt/Dbo/libwtdbo.a ./build-i386/src/Wt/Dbo/libwtdbo.a ./build-x86-64/src/Wt/Dbo/libwtdbo.a -output ${STAGE_DIR}/lib/libwtdbo.a
xcrun -sdk iphoneos lipo -create ./build-armv7/src/libwt.a ./build-armv7s/src/libwt.a ./build-arm64/src/libwt.a ./build-i386/src/libwt.a ./build-x86-64/src/libwt.a -output ${STAGE_DIR}/lib/libwt.a
