cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project( PeakFitImprove )

set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

if(MSVC)
  option(${PROJECT_NAME}_USE_MSVC_MultiThreadDLL "Use dynamically-link runtime library." OFF)

  if( ${PROJECT_NAME}_USE_MSVC_MultiThreadDLL)
    set(Boost_USE_STATIC_RUNTIME OFF)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  else()
    set(Boost_USE_STATIC_RUNTIME ON)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
  message(STATUS "Using CMAKE_MSVC_RUNTIME_LIBRARY ${CMAKE_MSVC_RUNTIME_LIBRARY}")
endif()


set( INCLUDE_ANALYSIS_TEST_SUITE ON CACHE BOOL "" )


set(BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "")
set(USE_DB_TO_STORE_SPECTRA ON CACHE BOOL "")
set(PERFORM_DEVELOPER_CHECKS ON CACHE BOOL "")
set(SpecUtils_ENABLE_EQUALITY_CHECKS ON CACHE BOOL "")
set(USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "")

set(TEST_SUITE_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test_data"
    CACHE STRING
    "Path to directory that contains the \"analysis_tests\" directory for saving N42 test states."
)

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

include(FetchContent)
FetchContent_Declare(
  openGA
  GIT_REPOSITORY git@github.com:Arash-codedev/openGA.git
  GIT_HASH        b2f292eaccfda77363bcc8f30c001392f3b71241 # latest as of 20240718
  GIT_SHALLOW    ON
)
FetchContent_MakeAvailable( openGA )


add_executable( PeakFitImprove PeakFitImprove.cpp )
target_link_libraries( PeakFitImprove PRIVATE InterSpecLib )

# For some reason ${openGA_SOURCE_DIR} doesnt look to be set - oh well, I'll just hard code things
target_include_directories( PeakFitImprove PUBLIC "${openGA_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}/_deps/openga-src/src/" )


if( CMAKE_GENERATOR STREQUAL "Xcode" )
  set_target_properties( PeakFitImprove PROPERTIES
    XCODE_GENERATE_SCHEME TRUE
    XCODE_SCHEME_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
endif( CMAKE_GENERATOR STREQUAL "Xcode" )
