cmake_minimum_required(VERSION 3.15...3.23)
 
set( SpecUtils_JAVA_SWIG OFF CACHE INTERNAL "" )
set( SpecUtils_USE_WT_THREADPOOL ON CACHE INTERNAL "" )
set( SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE INTERNAL "" )

set( ANDROID_CPP_FEATURES "rtti exceptions" CACHE INTERNAL "" )
set( ANDROID_PIE ON CACHE INTERNAL "" )
set( ANDROID_STL "c++_static" CACHE INTERNAL "" )
set( BOOST_EXCLUDE_LIBRARIES "locale" CACHE INTERNAL "" )


set( BOOST_LOCALE_ENABLE_ICONV OFF CACHE INTERNAL "" )
set( BOOST_LOCALE_ENABLE_ICU OFF CACHE INTERNAL "" )

set( BUILD_TESTING OFF CACHE INTERNAL "" )

set( InterSpec_FETCH_DEPENDENCIES ON CACHE INTERNAL "" )


set( BUILD_AS_LOCAL_SERVER OFF CACHE INTERNAL "" )
set( TRY_TO_STATIC_LINK ON CACHE INTERNAL "" )
set( MAX_SPECTRUM_MEMMORY_SIZE_MB 32 CACHE INTERNAL "" )
set( USE_DB_TO_STORE_SPECTRA ON CACHE INTERNAL "" )
set( USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE INTERNAL "" )
# Building the RelAct tool currently fails for Android because of erroring finding Eigen (which should be fetched) - I didnt spend time trying to fix. 
set( USE_REL_ACT_TOOL ON CACHE INTERNAL "" )
set( USE_LEAFLET_MAP ON CACHE INTERNAL "" )


message( "ANDROID_NDK_SYSROOT=${ANDROID_NDK_SYSROOT}")

project( InterSpecAndroid VERSION 1.0.12 )

add_definitions( -DANDROID )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

# TODO: test if ${log-lib} can replace log below
add_library( InterSpecAppLib SHARED android.cpp )
target_link_libraries( InterSpecAppLib PUBLIC Boost::filesystem wt InterSpecLib log )
target_include_directories( InterSpecAppLib PUBLIC ${wt_BINARY_DIR} )
