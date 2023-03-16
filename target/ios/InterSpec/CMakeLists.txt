cmake_minimum_required(VERSION 3.15...3.23)


set( BOOST_EXCLUDE_LIBRARIES "locale" CACHE INTERNAL "" )
set( BOOST_LOCALE_ENABLE_ICONV OFF CACHE INTERNAL "" )
set( BOOST_LOCALE_ENABLE_ICU OFF CACHE INTERNAL "" )
set( USE_BOOST_FRAMEWORK OFF CACHE INTERNAL "" )
#set( HAVE_LIBRT 0 CACHE INTERNAL "" )
#set( RT_required OFF CACHE INTERNAL "" )

set( InterSpec_IOS ON CACHE INTERNAL "" )
set( InterSpec_FETCH_DEPENDENCIES ON CACHE INTERNAL "" )


#set( CMAKE_XCODE_ATTRIBUTE_IPHONEOS_DEPLOYMENT_TARGET "11.0"
#        CACHE STRING "Set the minimum deployment target value."
#        FORCE
#    )

#For iOS we have to use libtool instead of ar so that the platform (armv7, armv7s, arm64, etc) are kept track of in the archive
#SET(CMAKE_AR "libtool" )

#CMake adds a cq after the archive command by default, we need to change this to "-static -o"
#SET(CMAKE_CXX_ARCHIVE_CREATE   "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")
#SET(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")

set(BUILD_AS_LOCAL_SERVER OFF CACHE INTERNAL "" )
set(PERFORM_DEVELOPER_CHECKS OFF CACHE INTERNAL "" )
set(TRY_TO_STATIC_LINK OFF CACHE INTERNAL "" )
set(USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE INTERNAL "" )
set(USE_DETECTION_LIMIT_TOOL OFF CACHE INTERNAL "" )
set(USE_QR_CODES ON CACHE INTERNAL "" )
set(USE_GOOGLE_MAP OFF CACHE INTERNAL "" )
set(USE_REMOTE_RID OFF CACHE INTERNAL "" )
set(USE_REL_ACT_TOOL OFF CACHE INTERNAL "" )
set(MAX_SPECTRUM_MEMMORY_SIZE_MB 32 CACHE INTERNAL "" )
set(USE_REL_ACT_TOOL ON CACHE INTERNAL "")
    


project( InterSpecIOS VERSION 1.10 )

add_definitions( -DIOS )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

add_library( InterSpecAppLib STATIC FileHandling.mm FileHandling.h )
target_link_libraries( InterSpecAppLib PUBLIC Boost::filesystem wt InterSpecLib )
target_include_directories( InterSpecAppLib PUBLIC ${wt_BINARY_DIR} )
