project(Cuba-3.0)
#right now this file only compiles "cuhre" version of the integration routines
cmake_minimum_required( VERSION 3.1...3.20 )



#The following would normally be setup and put into config.h when using autoconf
#  however for right now I'll just hardcode the values
IF(WIN32)
    set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_LONG_DOUBLE=1" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_LONG_DOUBLE_WIDER=1" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DNOUNDERSCORE" )
else(WIN32)
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_CONFIG_H=0" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_ERF=1" )
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_FORK=1" )
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_GETLOADAVG=1" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_LONG_DOUBLE=1" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_LONG_DOUBLE_WIDER=1" )
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_POWL=1" )
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DHAVE_SHMGET=1" )
    set( BUILD_FLAGS "${BUILD_FLAGS} -DNOUNDERSCORE" )
#   set( BUILD_FLAGS "${BUILD_FLAGS} -DLONGLONGINT" )
endif(WIN32)


#in order to compile with clang, we have to do a kinda shady workaround, see
#  Rule.c and decl.h
set( BUILD_FLAGS "${BUILD_FLAGS} -DNOUNDERSCORE -DMAX_NDIM=3 -DMAX_NCOMP=1" )

#SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3 -fomit-frame-pointer -ffast-math -Wall" )

set( CUBA_CUHRE_SRCS src/cuhre/Cuhre.c src/common/WorkerIni.c )

SET_SOURCE_FILES_PROPERTIES( ${CUBA_CUHRE_SRCS} PROPERTIES LANGUAGE C)
add_library(Cuba-3.0 STATIC ${CUBA_CUHRE_SRCS})
set_target_properties( Cuba-3.0 PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
target_include_directories( Cuba-3.0 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE src/common src/cuhre )




