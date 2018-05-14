cmake_minimum_required(VERSION 3.0)

#Note 20180508: References to "D3" or Java capabilities should be kept to off
#  as this has not been included in the public release of InterSpec.  It is 
#  expected to be included shortly.

#This file defines two useful variables in its parent scope:
#  SpecUtils_INCLUDE_DIRS: The list of directories that should be in the include
#                          path of the parent project to this library (SpecUtils
#                          directory and boost directory)
#  SpecUtils_LINK_TO_LIBS: The list of libraries the parent project needs to
#                          link to (boost libs, libSpecUtils, etc)

project(SpecUtils)

if(NOT ANDROID)
set( CMAKE_CXX_STANDARD 11 )
#set( CMAKE_CXX_STANDARD_REQUIRED ON )
endif(NOT ANDROID)


set( LIBRARY_OUTPUT_PATH lib CACHE PATH "Path to libraries" FORCE )
set( UPDATE_TYPE "svn" )
set( CMAKE_COLOR_MAKEFILE ON CACHE BOOL "Color" FORCE )

set( UPACKAGE "SpecUtils" )
set( PACKAGE "SpecUtils" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_FILESYSTEM_VERSION=3" )

option( ENABLE_D3_CHART_EXPORTING "Enables exporting of the D3 chart format" OFF )
option( SpecUtils_D3_SUPPORT_FILE_STATIC "Compiles the JS and CSS files into memory, rather than reading files on disk; disable to allow changing JS/CSS without having to recompile" OFF )

option( PERFORM_DEVELOPER_CHECKS "Performs additional computationally expensive tests during execution" OFF )
option( BUILD_AS_UNIT_TEST_SUITE "Builds as the test sweet target" OFF )
option( SPECTRUM_DATA_STRUCTS_NO_BOOST_LIB "Use boost header-only libraries in SpectrumDataStructs.h/.cpp.  Enabling this removes mutlithreadded file opening." OFF )

option( SpecUtils_USE_WT_THREADPOOL "Whether or not to use the Wt thread pool (only for use within InterSpec)" OFF )

option( REBIN_FILES_TO_SINGLE_BINNING "Rebins all spectra in spectrum file to the same energy calibration" OFF )
option( TRY_TO_STATIC_LINK "Try to link to static libs whenever possible" OFF )
option( SpecUtils_BUILD_TEST "Creates executuable to perform interactive regression test" OFF )
option( SpecUtils_PYTHON_BINDINGS "Creates python bindings to the c++ code" OFF )
option( SpecUtils_JAVA_SWIG "Creates swig/java bindings to the c++ code" OFF )

#This is a temporary (20170629) hack to deal with error messages... should be handled a lot better...
option( SpecUtils_EXTERNALLY_DEFINED_LOG_MESSAGE "If true, you must have a externally defined log_error_message(string,string,int) function defined when linking" OFF )

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if( TRY_TO_STATIC_LINK )
  set( Boost_USE_STATIC_LIBS ON )
  if( WIN32 )
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else( WIN32 )
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif( WIN32 )
else( TRY_TO_STATIC_LINK )
  set( Boost_USE_STATIC_LIBS OFF )
  set( CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES} )
endif( TRY_TO_STATIC_LINK )


set( boost_libs_required date_time system filesystem )


if( SpecUtils_BUILD_TEST OR BUILD_AS_UNIT_TEST_SUITE )
  set( boost_libs_required ${boost_libs_required} program_options )
endif( SpecUtils_BUILD_TEST OR BUILD_AS_UNIT_TEST_SUITE )

if( SpecUtils_PYTHON_BINDINGS )
  set( boost_libs_required ${boost_libs_required} python )
endif( SpecUtils_PYTHON_BINDINGS )


if( SpecUtils_JAVA_SWIG )
  FIND_PACKAGE( SWIG REQUIRED )
  FIND_PACKAGE( JNI REQUIRED )
  INCLUDE( ${SWIG_USE_FILE} )
  INCLUDE_DIRECTORIES( ${JAVA_INCLUDE_PATH} )
  INCLUDE_DIRECTORIES( ${JAVA_INCLUDE_PATH}/linux )
  INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/swig )
  SET( CMAKE_SWIG_FLAGS -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_BINARY_DIR} -DSpecUtils_JAVA_SWIG=1 )
  SET_SOURCE_FILES_PROPERTIES(swig/SpecUtils.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(swig/SpecUtils.i PROPERTIES SWIG_FLAGS "")
  SWIG_ADD_MODULE( SpecUtilsJni java swig/SpecUtils.i )
  SWIG_LINK_LIBRARIES( SpecUtilsJni ${JAVA_LIBRARIES} SpecUtils )
endif( SpecUtils_JAVA_SWIG )


find_package( Boost REQUIRED COMPONENTS ${boost_libs_required} )

if( NOT Boost_FOUND )
  message(FATAL_ERROR "Couldnt Find Boost")
endif( NOT Boost_FOUND )


find_path( RAPID_XML_DIR rapidxml.hpp ${CMAKE_CURRENT_SOURCE_DIR}/../rapidxml ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml )

#RAPID_XML_DIR isnt found when compiling for iOS...
#IF( ${RAPID_XML_DIR} MATCHES "RAPID_XML_DIR-NOTFOUND" )
#  message( "Current CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}" )
#  message( FATAL_ERROR "Couldnt find RAPID_XML_DIR" )
#ENDIF( ${RAPID_XML_DIR} MATCHES "RAPID_XML_DIR-NOTFOUND" )

#IF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
#  SET( PYTHON_LIBRARY /System/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib )
#  SET( PYTHON_INCLUDE_DIR /System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 )
#ENDIF()

if( SpecUtils_PYTHON_BINDINGS )
  FIND_PACKAGE(PythonLibs 2.7)
  IF(NOT PYTHONLIBS_FOUND)
    MESSAGE(FATAL_ERROR "Unable to find PythonLibs.")
  ENDIF()
endif( SpecUtils_PYTHON_BINDINGS )

INCLUDE_DIRECTORIES(
  SpecUtils
  ${RAPID_XML_DIR}/..
  ${CMAKE_CURRENT_BINARY_DIR}
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${PYTHON_INCLUDE_DIRS}
)


set( sources
     src/SpectrumDataStructs.cpp
     src/UtilityFunctions.cpp
     src/SpecUtilsAsync.cpp
)

set( headers
     SpecUtils/SpectrumDataStructs.h
     SpecUtils/UtilityFunctions.h
     SpecUtils/SpecUtilsAsync.h
     SpecUtils/SpecUtils_config.h.in
     SpecUtils/DetectiveModelFromSerial.hpp
)


if( ENABLE_D3_CHART_EXPORTING )
  list( APPEND sources src/D3SpectrumExport.cpp )
  list( APPEND headers SpecUtils/D3SpectrumExport.h )
  list( APPEND headers SpecUtils/D3SpectrumExportResources.h.in )


# IMPORTANT: These 2 lists below must be of same length, corresponding with each index
#   For outputting necessary library in D3 HTML files
  set( D3_SUPPORT_FILES
       d3_resources/d3.v3.min.js
       d3_resources/c.min.js
       d3_resources/numeric-1.2.6.min.js
       d3_resources/PeakFit.js
       d3_resources/SpectrumChartD3.js
       d3_resources/SpectrumChartD3.css
  )

  set( D3_SUPPORT_FILE_MACROS
       D3_MIN_JS
       CASSOWARY_MIN_JS
       NUMERIC_MIN_JS
       PEAK_FIT_JS
       SPECTRUM_CHART_D3_JS
       SPECTRUM_CHART_D3_CSS
  )

# https://stackoverflow.com/questions/7932205/parallel-iteration-over-lists-in-makefile-or-cmake-file
  list(LENGTH D3_SUPPORT_FILES D3_SUPPORT_FILES_LENGTH)
  math(EXPR D3_SUPPORT_FILES_LAST_INDEX "${D3_SUPPORT_FILES_LENGTH} - 1")

# Set the corresponding file content for each library file
  foreach(index RANGE ${D3_SUPPORT_FILES_LAST_INDEX})
    list(GET D3_SUPPORT_FILES ${index} D3_SUPPORT_FILE_NAME)
    list(GET D3_SUPPORT_FILE_MACROS ${index} D3_SUPPORT_FILE_MACRO_NAME)

    if( SpecUtils_D3_SUPPORT_FILE_STATIC )
      #Take from https://gist.github.com/sivachandran/3a0de157dccef822a230 to read file into hex array
      # reads source file contents as hex string
      file(READ ${D3_SUPPORT_FILE_NAME} hexString HEX)
      string(LENGTH ${hexString} hexStringLength)
      set(hexString "${hexString}00")
      math(EXPR arraySize "${hexStringLength} / 2")

      # adds '0x' prefix and comma suffix before and after every byte respectively
      string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1, " arrayValues ${hexString})
      # removes trailing comma
      string(REGEX REPLACE ", $" "" arrayValues ${arrayValues})

      # declares byte array and the length variables
      set( ${D3_SUPPORT_FILE_MACRO_NAME}_ARR "const unsigned char ${D3_SUPPORT_FILE_MACRO_NAME}[] = { ${arrayValues} };")
      set( ${D3_SUPPORT_FILE_MACRO_NAME}_LEN "const size_t ${D3_SUPPORT_FILE_MACRO_NAME}_SIZE = ${arraySize};")
    else( SpecUtils_D3_SUPPORT_FILE_STATIC )
      set( ${D3_SUPPORT_FILE_MACRO_NAME}_FILE "const char * const ${D3_SUPPORT_FILE_MACRO_NAME}_FILENAME = \"${CMAKE_CURRENT_SOURCE_DIR}/${D3_SUPPORT_FILE_NAME}\";" )
    endif( SpecUtils_D3_SUPPORT_FILE_STATIC )
  endforeach()

#We will create D3SpectrumExportResources.h that holds all the js/css, but this
#  is only done at CMake generation time, not build time, so can get out of date
#  See https://cmake.org/pipermail/cmake/2012-May/050221.html for a potential fix
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/SpecUtils/D3SpectrumExportResources.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/D3SpectrumExportResources.h
  )

  list( APPEND OTHER_SUPPORT_FILES
        d3_resources/PeakFit.js
        d3_resources/SpectrumChartD3.js
        d3_resources/SpectrumChartD3.min.index.js
        d3_resources/SpectrumChartD3.css
  )
endif( ENABLE_D3_CHART_EXPORTING )




set( LIBRARIES_TO_LINK_TO
  ${LIBRARIES_TO_LINK_TO}
  ${Boost_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
)


if( SpecUtils_BUILD_TEST OR BUILD_AS_UNIT_TEST_SUITE )
set( LIBRARIES_TO_LINK_TO ${LIBRARIES_TO_LINK_TO} ${Boost_PROGRAM_OPTIONS_LIBRARY} )
endif( SpecUtils_BUILD_TEST OR BUILD_AS_UNIT_TEST_SUITE )

if( SpecUtils_PYTHON_BINDINGS )
  set( sources ${sources} src/SpectrumDataStructs_py.cpp )
  set( LIBRARIES_TO_LINK_TO ${LIBRARIES_TO_LINK_TO} ${PYTHON_LIBRARIES} ${Boost_PYTHON_LIBRARY} )
endif( SpecUtils_PYTHON_BINDINGS )

if( SpecUtils_JAVA_SWIG )
  set( sources ${sources} swig/SpecUtils.i )
  set( LIBRARIES_TO_LINK_TO ${LIBRARIES_TO_LINK_TO} ${JAVA_LIBRARIES} )
endif( SpecUtils_JAVA_SWIG )


if( SpecUtils_PYTHON_BINDINGS OR SpecUtils_JAVA_SWIG )
  add_library( SpecUtils SHARED ${sources} ${headers} ${OTHER_SUPPORT_FILES} )
else( SpecUtils_PYTHON_BINDINGS OR SpecUtils_JAVA_SWIG )
  add_library( SpecUtils STATIC ${sources} ${headers} ${OTHER_SUPPORT_FILES} )
  set_target_properties( SpecUtils PROPERTIES PREFIX "lib" OUTPUT_NAME "SpecUtils" )
endif( SpecUtils_PYTHON_BINDINGS OR SpecUtils_JAVA_SWIG )


target_link_libraries( SpecUtils ${LIBRARIES_TO_LINK_TO} )

#if( NOT SpecUtils_PYTHON_BINDINGS AND NOT SpecUtils_JAVA_SWIG )
  # If we arent swiging or linking to python, then lets build a static version
  # of the library as well.
#  ADD_LIBRARY( SpecUtils-static STATIC ${sources} ${headers} )
#  SET_TARGET_PROPERTIES( SpecUtils-static PROPERTIES OUTPUT_NAME "SpecUtils")
#  SET_TARGET_PROPERTIES( SpecUtils-static PROPERTIES PREFIX "lib" )
#  target_link_libraries( SpecUtils-static ${LIBRARIES_TO_LINK_TO} )
#endif()


set( SpecUtils_LIBS SpecUtils ${LIBRARIES_TO_LINK_TO} CACHE INTERNAL "Libraries necassary to link with SpecUtils" FORCE )
set( SpecUtils_INC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
     "${CMAKE_CURRENT_SOURCE_DIR}/SpecUtils"
     ${CMAKE_CURRENT_BINARY_DIR}
     ${Boost_INCLUDE_DIRS}
     ${RAPID_XML_DIR}/.. CACHE INTERNAL "Include directories for SpecUtils library" FORCE )



if( SpecUtils_BUILD_TEST )
  if( NOT PERFORM_DEVELOPER_CHECKS )
    message( FATAL_ERROR "PERFORM_DEVELOPER_CHECKS must be enabled when building testing executable" )
  endif( NOT PERFORM_DEVELOPER_CHECKS )

  add_executable( regression_test testing/regression_test.cpp )
  target_link_libraries( regression_test SpecUtils ${LIBRARIES_TO_LINK_TO} )
endif( SpecUtils_BUILD_TEST )


configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/SpecUtils/SpecUtils_config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/SpecUtils_config.h
)


get_directory_property(hasParent PARENT_DIRECTORY)
if( hasParent  )
  list( APPEND INCLUDE_DIRS ${CMAKE_BINARY_DIR} ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} )
  set( SpecUtils_INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE )
  set( SpecUtils_LINK_TO_LIBS ${LIBRARIES_TO_LINK_TO} SpecUtils PARENT_SCOPE )
endif( hasParent )


