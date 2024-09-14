cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project( InterSpec_batch )


# If we are the top-level project, we'll add libInterSpec to the build
get_directory_property(hasParent PARENT_DIRECTORY)


if( NOT hasParent )
  # These options have not been tested
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  set( BUILD_AS_WX_WIDGETS_APP ON CACHE BOOL "Set to build wxWidget target" )
  set( TRY_TO_STATIC_LINK ON CACHE BOOL "Lets statically link everything" )
  set( BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "N/A" )
  set( USE_OSX_NATIVE_MENU OFF CACHE BOOL "N/A" )
  set( USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE BOOL "Spectrum File Query Widget" )

  if(MSVC)
    add_definitions(-DBOOST_ALL_NO_LIB) #Prevent boost auto-linking, which seems to call in vc141 boost libs instead of vc142
    add_definitions(-D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS ) #0x0601==Win7
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

  add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )
endif( NOT hasParent )

set( sources batch_main.cpp )

add_executable( InterSpec_batch ${sources} )

target_link_libraries( InterSpec_batch PUBLIC InterSpecLib )
set_target_properties( InterSpec_batch PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )
set_target_properties( InterSpec_batch PROPERTIES OUTPUT_NAME "InterSpec_batch" )

# Turn off "note: The practice of declaring the Bind placeholders (_1, _2, ...) in the global namespace is deprecated...."
# warnings (these come from Wt I think, not our code)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -D BOOST_BIND_GLOBAL_PLACEHOLDERS")

if(WIN32)
  # Allow using wildcards from the command-line on Windows.
  set_target_properties( InterSpec_batch PROPERTIES LINK_FLAGS "wsetargv.obj")
endif(WIN32)

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

if( APPLE )
  install( TARGETS InterSpec_batch
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
  )
else( APPLE )
  set( RESOURCE_DEST_DIR "bin/.." )
  
  set(CPACK_BINARY_ZIP ON)
  set(CPACK_BINARY_NSIS OFF)
  set(CPACK_SOURCE_ZIP ON)
  set(CPACK_SOURCE_TGZ ON) #linux

  # We get an error if the destination is empty, or just ".", but the following seems to work
  install( TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "bin/.." )
endif( APPLE )

if( NOT hasParent )
  install(
    DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources
      ${CMAKE_CURRENT_SOURCE_DIR}/../../data
      ${CMAKE_CURRENT_SOURCE_DIR}/../../example_spectra
      ${WT_RESOURCES_DIRECTORY}
    DESTINATION
      ${RESOURCE_DEST_DIR}
  )

  install(
    FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/../../NOTICE.html
    DESTINATION
      ${RESOURCE_DEST_DIR}
  )

  install(
    FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/d3.v3.min.js
      ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
      ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
    DESTINATION
      "${RESOURCE_DEST_DIR}/InterSpec_resources"
  )

  install(
    FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SandiaDecay/sandia.decay.min.xml
    DESTINATION
      "${RESOURCE_DEST_DIR}/data"
    RENAME
      "sandia.decay.xml"
  )


  set(CPACK_PACKAGE_NAME "InterSpec")
  set(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
      "InterSpec assists in analyzing spectral nuclear radiation data, using a peak-based methodology.")
  set(CPACK_PACKAGE_VERSION "${INTERSPEC_VERSION}")
  #set(CPACK_PACKAGE_VERSION_MAJOR "1")
  #set(CPACK_PACKAGE_VERSION_MINOR "0")
  #set(CPACK_PACKAGE_VERSION_PATCH "0")
  #CPACK_INSTALL_PREFIX
  #set(CPACK_PACKAGE_INSTALL_DIRECTORY "release_builds/InterSpec-${CMAKE_SYSTEM_NAME}-x64")

  include(CPack)
endif( NOT hasParent )
