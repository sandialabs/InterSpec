#adapted from https://github.com/Kitware/CMake/tree/master/Tests/CFBundleTest
cmake_minimum_required(VERSION 3.1.0)

#project( SpecFilePreview )


set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( CMAKE_OSX_DEPLOYMENT_TARGET "10.10" )
set( Boost_USE_STATIC_LIBS ON )
SET( CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES} )

set( TRY_TO_STATIC_LINK ON CACHE BOOL "Try to static link" FORCE )

INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake/FindWt.cmake )
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake/WtFindBoost-cmake.txt )
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake/WtFindZlib.txt )

set( SOURCES
    main.c
    GenerateThumbnailForURL.m
    GeneratePreviewForURL.m
    	
    QLSpectrumChart.h
    QLSpectrumChart.cpp
    SpecPreviewCommon.h
    SpecPreviewCommon.cpp
      
    QLPeakDef.h
    QLPeakDef.cpp
    QLSpecMeas.h
    QLSpecMeas.cpp

    QLSpectrumDataModel.h
    QLSpectrumDataModel.cpp
    
    Info.plist
)

add_library( SpecFilePreview MODULE ${SOURCES} )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../../../../external_libs/SpecUtils ${CMAKE_CURRENT_BINARY_DIR}/SpecUtils )

set(sysroot)
if(CMAKE_OSX_SYSROOT)
  set(sysroot -isysroot ${CMAKE_OSX_SYSROOT})
endif()


FIND_LIBRARY( HPDF_LIBRARY NAME hpdfs HINTS ${Wt_INCLUDE_DIR}/../lib )
if( NOT HPDF_LIBRARY )
  message( FATAL "Could not find libharu - if Wt was built with PDF support, you will get link errors" )
endif()
  
FIND_LIBRARY( PNG_LIBRARY NAME png HINTS ${Wt_INCLUDE_DIR}/../lib )
if( NOT PNG_LIBRARY )
  message( FATAL "Could not find libpng - requred for libharu" )
endif()
  

# note that for some reason, the makefile and xcode generators use a different
# property to indicate where the Info.plist file is :-/ For that reason, we
# specify it twice so it will work both places
set_target_properties( SpecFilePreview PROPERTIES
    BUNDLE 1
    BUNDLE_EXTENSION qlgenerator
    XCODE_ATTRIBUTE_MACH_O_TYPE mh_bundle
    XCODE_ATTRIBUTE_INFOPLIST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist
    )
    #LINK_FLAGS "-Wl,-exported_symbols_list,\"${CMAKE_CURRENT_SOURCE_DIR}/ExportList_plugin.txt\"")

#export(TARGETS SpecFilePreview FILE ${CMAKE_CURRENT_BINARY_DIR}/exp.cmake)

target_link_libraries( SpecFilePreview PUBLIC 
  Threads::Threads
  SpecUtils
  ${HPDF_LIBRARY}
  ${PNG_LIBRARY}
  ${LIBZ_LIBRARY}
  ${ZLIB_LIBRARIES}
  ${Wt_LIBRARIES}
  ${Boost_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_RANDOM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  /System/Library/Frameworks/QuickLook.framework/QuickLook
  /System/Library/Frameworks/CoreServices.framework/CoreServices
)


target_include_directories( SpecFilePreview PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../../../external_libs 
                            ${Wt_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )
message( "SpecFilePreview: ${CMAKE_CURRENT_SOURCE_DIR}/../../../../external_libs" )


#install(TARGETS SpecFilePreview DESTINATION InterSpec.app/Contents/Library/QuickLook/ )
