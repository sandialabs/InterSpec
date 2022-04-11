cmake_minimum_required(VERSION 3.15...3.23)

#include( ../../cmake/FetchInterSpecDeps.txt )
 
set( ANDROID_CPP_FEATURES "rtti exceptions" CACHE INTERNAL "" )
set( ANDROID_PIE ON CACHE INTERNAL "" )
set( ANDROID_STL "c++_static" CACHE INTERNAL "" )
set( BOOST_EXCLUDE_LIBRARIES "locale" CACHE INTERNAL "" )


set( BOOST_LOCALE_ENABLE_ICONV OFF CACHE INTERNAL "" )
set( BOOST_LOCALE_ENABLE_ICU OFF CACHE INTERNAL "" )


set( InterSpec_FETCH_DEPENDENCIES ON CACHE INTERNAL "" )

message( "ANDROID_NDK_SYSROOT=${ANDROID_NDK_SYSROOT}")

project( InterSpecAndroid VERSION 1.10 )


add_definitions( -DANDROID )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )


add_library( InterSpecAppLib SHARED main.cpp )
target_link_libraries( InterSpecAppLib PUBLIC Boost::filesystem wt InterSpecLib log )
target_include_directories( InterSpecAppLib PUBLIC ${wt_BINARY_DIR} )
