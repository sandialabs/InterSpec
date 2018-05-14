# Adapt the pathnames in this file to suit your needs

# The path to the directory that contains your cross-compiled dependency
# libraries (boost, and optionally others such as graphicsmagick, ssl, ...)
#SET(ANDROID_STAGING_DIR /Users/wcjohns/install/android/arm-linux-androideabi-4.6 )
SET(ANDROID_STAGING_DIR /Users/wcjohns/install/android/armstandalone19abi4.8/sysroot/usr)

# The path to your Android SDK installation
SET(ANDROID_SDK_DIR /Users/wcjohns/install/android/android-sdk-macosx )

# The ID of the target must be compatible with the NDK used to build Wt.
# Target ID's are specific to your installation; to discover your target ID,
# use 'android list target'. The default value of 99 is chosen to ensure
# that you will verify this.
SET(ANDROID_SDK_TARGET_ID 6) 

# The path to your standalone toolchain directory.
# See docs/STANDALONE-TOOLCHAIN.html in your ndk's installation dir
#SET(ANDROID_NDK_TOOLS_DIR /Users/wcjohns/install/android/arm-linux-androideabi-4.6)
SET(ANDROID_NDK_TOOLS_DIR /Users/wcjohns/install/android/armstandalone19abi4.8)

# Below is the normal contents of a cmake toolchain file, you'll
# probably not need to modify it.

SET(TARGET_CC ${ANDROID_NDK_TOOLS_DIR}/bin/arm-linux-androideabi-gcc)
SET(TARGET_CXX ${ANDROID_NDK_TOOLS_DIR}/bin/arm-linux-androideabi-g++)
SET(ANDROID_STRIP ${ANDROID_NDK_TOOLS_DIR}/bin/arm-linux-androideabi-strip)

#SET(CMAKE_CXX_FLAGS "-fpic -ffunction-sections -funwind-tables -D__ARM_ARCH_5__ -D__ARM_ARCH_5T__ -D__ARM_ARCH_5E__ -D__ARM_ARCH_5TE__  -Wno-psabi -march=armv5te -mtune=xscale -msoft-float -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300 -DANDROID  -Wa,--noexecstack")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -finline-limit=300 -DANDROID" )


#To try next
#SET(CMAKE_CXX_FLAGS "-DGL_GLEXT_PROTOTYPES -fPIC -DANDROID -mthumb -Wno-psabi")
#SET(CMAKE_C_FLAGS "-DGL_GLEXT_PROTOTYPES -fPIC -DANDROID -mthumb -Wno-psabi")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a -mfloat-abi=softfp")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv7-a -mfloat-abi=softfp")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpu=neon")

#Flags extracted from androideabi-4.9 from the armv7-a section of make file
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DANDROID -fpic -ffunction-sections -funwind-tables -fstack-protector -no-canonical-prefixes" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv7-a -mfpu=vfpv3-d16 -mfloat-abi=softfp" )
#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon" )

SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -g -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -funswitch-loops -finline-limit=300" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -UNDEBUG -fno-omit-frame-pointer -fno-strict-aliasing" )

SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -no-canonical-prefixes -march=armv7-a -Wl,--fix-cortex-a8" )
SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-canonical-prefixes -march=armv7-a -Wl,--fix-cortex-a8" )
#SET( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -no-canonical-prefixes -march=armv7-a -Wl,--fix-cortex-a8" )
#SET( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -no-canonical-prefixes -march=armv7-a -Wl,--fix-cortex-a8" )


# this one is important
SET(CMAKE_SYSTEM_NAME Linux)

#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
SET(CMAKE_C_COMPILER ${TARGET_CC})

SET(CMAKE_CXX_COMPILER ${TARGET_CXX})

# where is the target environment
SET(CMAKE_FIND_ROOT_PATH ${ANDROID_NDK_TOOLS_DIR} ${ANDROID_STAGING_DIR})

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

