project(Cuba-3.0)
#right now this file only compiles "cuhre" version of the integration routines

#-----------------------------------------------------------------------------
# CMake administrivia
cmake_minimum_required( VERSION 2.8 )
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND cmake_policy )


set( LIBRARY_OUTPUT_PATH lib CACHE PATH "Path to libraries" FORCE )
set( UPDATE_TYPE "svn" )
set( CMAKE_COLOR_MAKEFILE ON CACHE BOOL "Color" FORCE )



SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3 -fomit-frame-pointer -ffast-math -Wall" )

#The following would normally be setup and put into config.h when using autoconf
#  however for right now I'll just hardcode the values
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_CONFIG_H=0" )
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_ERF=1" )
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_FORK=1" )
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_GETLOADAVG=1" )
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_LONG_DOUBLE=1" )
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_LONG_DOUBLE_WIDER=1" )
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_POWL=1" )
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DHAVE_SHMGET=1" )
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNOUNDERSCORE" )
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLONGLONGINT" )


INCLUDE_DIRECTORIES(
  .
  ../common
  ../..
)


set ( CUBA_SUAVE_SRCS
Suave.c
../common/WorkerIni.c
)

SET_SOURCE_FILES_PROPERTIES( ${CUBA_SUAVE_SRCS} PROPERTIES LANGUAGE C)
add_library(CubaSuave3.0 STATIC ${CUBA_SUAVE_SRCS})






