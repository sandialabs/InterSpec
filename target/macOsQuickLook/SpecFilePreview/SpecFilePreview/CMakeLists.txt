#adapted from https://github.com/Kitware/CMake/tree/master/Tests/CFBundleTest
cmake_minimum_required(VERSION 3.15)

project( SpecFilePreview )

set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")

find_package( ZLIB REQUIRED )
find_package( Threads REQUIRED )

if( NOT DEFINED InterSpec_BASE_DIR )
  find_file( wt_find_file 
          NAMES "FindWt.cmake" 
          PATHS "${InterSpec_BASE_DIR}/cmake" 
                "../../../cmake" "../../../../cmake" 
                NO_DEFAULT_PATH 
                REQUIRED 
  )
  get_filename_component( InterSpec_BASE_DIR ${wt_find_file} DIRECTORY )
  get_filename_component( InterSpec_BASE_DIR ${InterSpec_BASE_DIR} DIRECTORY )
endif( NOT DEFINED InterSpec_BASE_DIR )


message( "InterSpec_BASE_DIR=${InterSpec_BASE_DIR}" )
message( "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}" )

INCLUDE( "${InterSpec_BASE_DIR}/cmake/FindWt.cmake" )



if( InterSpec_FETCH_DEPENDENCIES )
  set( Wt_INCLUDE_DIR "${wt_SOURCE_DIR}/src" )
else( InterSpec_FETCH_DEPENDENCIES )
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_MULTITHREADED ON)
  find_package(Boost 1.68 REQUIRED COMPONENTS thread date_time system filesystem program_options regex random )
endif( InterSpec_FETCH_DEPENDENCIES )


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

add_subdirectory( ${InterSpec_BASE_DIR}/external_libs/SpecUtils ${CMAKE_CURRENT_BINARY_DIR}/SpecUtils )

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
  
find_library(QUICK_LOOK QuickLook)
if (NOT QUICK_LOOK)
    message(FATAL_ERROR "QuickLook not found")
endif()

find_library(CORE_SERVICES CoreServices)
if (NOT CORE_SERVICES)
    message(FATAL_ERROR "CoreServices not found")
endif()

MARK_AS_ADVANCED(QUICK_LOOK CORE_SERVICES)


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
  ZLIB::ZLIB
  #wt wthttp wtdbosqlite3 wtdbo
  ${Wt_LIBRARIES}
  ${Wt_LIBRARY}
  ${Wt_HTTP_LIBRARY}
  ${Wt_DBO_LIBRARY}
  ${Wt_DBOSQLITE3_LIBRARY}
  Boost::system
  Boost::regex
  Boost::thread
  Boost::filesystem
  Boost::random
  Boost::program_options
  Boost::date_time
  ${QUICK_LOOK}
  ${CORE_SERVICES}
)


target_include_directories( SpecFilePreview 
  PUBLIC  
    ${Wt_INCLUDE_DIR}
    ${InterSpec_BASE_DIR}/external_libs
)


#install(TARGETS SpecFilePreview DESTINATION InterSpec.app/Contents/Library/QuickLook/ )
