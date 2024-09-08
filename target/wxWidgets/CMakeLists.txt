cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project( wxInterSpec )

set( BUILD_AS_WX_WIDGETS_APP ON CACHE BOOL "Set to build wxWidget target" )
set( TRY_TO_STATIC_LINK ON CACHE BOOL "Lets statically link everything" )
set( BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "N/A" )
set( USE_OSX_NATIVE_MENU OFF CACHE BOOL "N/A" )
set( USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "Enable Spectrum File Query Widget" )
set( USE_BATCH_TOOLS ON CACHE BOOL "Enables (command line for now) batch processing." )
#option( USE_BATCH_TOOLS "Enables (command line for now) batch processing." ON )

IF(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) #Prevent boost auto-linking, which seems to call in vc141 boost libs instead of vc142
  add_definitions(-D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS ) #0x0601==Win7
ENDIF(WIN32)
 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if( APPLE )
  set( EXE_TYPE "MACOSX_BUNDLE" )
elseif(WIN32)
  set( EXE_TYPE "WIN32" )
endif( APPLE )

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

  add_definitions( -DwxUSE_RC_MANIFEST=1 -DwxUSE_DPI_AWARE_MANIFEST=2 )
endif()


if( NOT InterSpec_FETCH_DEPENDENCIES )
  # I havent gotten the above to find_package(...) call to actually work yet
  # with wxWidgets I compile with CMake on Windows, so we'll keep compiling
  # as a subdirectory for now... 
  set( wxWidgets_USE_STATIC ON CACHE BOOL "" )
  find_package( wxWidgets REQUIRED core base webview stc )

#  set(wxBUILD_SHARED OFF CACHE BOOL "")
#  if(WIN32)
#    set(wxUSE_WEBVIEW_EDGE ON CACHE BOOL "")
#    set(wxUSE_WEBVIEW_EDGE_STATIC ON CACHE BOOL "")
#    # TODO: static runtime should be set off of wxInterSpec_USE_MSVC_MultiThreadDLL
#    set(wxBUILD_USE_STATIC_RUNTIME ON CACHE BOOL "")
#  endif(WIN32)

  # We could pre-download the source code add add it, like:
  #add_subdirectory(3rd_party/wxWidgets-3.2.1 EXCLUDE_FROM_ALL)

#  # But for the moment we'll use CMake FetchContent to simplify other build scripts
#  message( "Using FetchContent to retrieve wxWidgets" )
#  include(FetchContent)
#  FetchContent_Declare( wxWidgets
#    GIT_REPOSITORY https://github.com/wxWidgets/wxWidgets.git
#    GIT_TAG        97e99707c5d2271a70cb686720b48dbf34ced496 # release-3.2.1
#    GIT_SHALLOW    ON
#  )

#  FetchContent_GetProperties(wxWidgets)
#  if(NOT wxWidgets_POPULATED)
#    FetchContent_Populate(wxWidgets)
#    message( "wxWidgets_SOURCE_DIR=${wxWidgets_SOURCE_DIR}")
#    message( "wxWidgets_BINARY_DIR=${wxWidgets_BINARY_DIR}")
#    add_subdirectory( ${wxWidgets_SOURCE_DIR} ${wxWidgets_BINARY_DIR} EXCLUDE_FROM_ALL)
#  endif(NOT wxWidgets_POPULATED)
endif( NOT InterSpec_FETCH_DEPENDENCIES )


add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )
 
if(WIN32)
  list(APPEND sources
    ../electron/windows/InterSpec.rc
    windows/InterSpec.manifest
  )
endif(WIN32)
 
list(APPEND sources
  InterSpecWxUtils.h
  InterSpecWxUtils.cpp
  InterSpecWxApp.h
  InterSpecWxApp.cpp
  InterSpecWebFrame.h
  InterSpecWebFrame.cpp
  wxMain.cpp
)

add_executable(${PROJECT_NAME} ${EXE_TYPE} ${sources})


if( USE_BATCH_TOOLS )
  add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../batch InterSpec_batch )
  add_dependencies( ${PROJECT_NAME} InterSpec_batch )

  # Post-build command to copy the DLL to the executable directory
  # (we build LibInterSpec as a static library if NOT using batch
  #  tools, but a shared library if including batch tools)
  if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      $<TARGET_FILE:InterSpecLib>
      $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  endif(WIN32)
endif( USE_BATCH_TOOLS )

if( InterSpec_FETCH_DEPENDENCIES )
  target_link_libraries(${PROJECT_NAME} PUBLIC wx::core wx::base wx::webview wx::stc ${wxWidgets_LIBRARIES} )
else( InterSpec_FETCH_DEPENDENCIES )
  # Not sure what wxWidgets_USE_FILE actually buys us 
  include( ${wxWidgets_USE_FILE} )

  # Not sure wx::core, wx::base, etc wont work in this case
  if( APPLE )
    # macOS needs wx::net, I am guessing for wxSingleInstanceChecker
    target_link_libraries(${PROJECT_NAME} PUBLIC ${wxWidgets_LIBRARIES} "-lwx_baseu_net-3.2")
  else( APPLE )
      target_link_libraries(${PROJECT_NAME} PUBLIC ${wxWidgets_LIBRARIES})
  endif( APPLE )
  target_include_directories( ${PROJECT_NAME} PUBLIC ${wxWidgets_INCLUDE_DIRS} )

  # We need to let the linker know where to find WebView2LoaderStatic.lib
  # TODO: use find file to find this better - for the moment we'll just assume we copied this file to the root lib directory
  target_link_directories( ${PROJECT_NAME} PUBLIC "${wxWidgets_ROOT_DIR}/lib" )
endif( InterSpec_FETCH_DEPENDENCIES )

target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_JS_LIB} InterSpecLib )
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "InterSpec")


set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)


if( APPLE )
  # TODO: for macOS
  #  - Create proper Info.plist
  #  - App icon doesnt seem to actually be set correctly
  #  - Need to have CMake "package" target be dependant on the ${PROJECT_NAME} executable
  #  - Remove the top os titlebar
  
  
  set( CPACK_GENERATOR "ZIP;DragNDrop" )
  set( CPACK_BINARY_ZIP ON )
  set( CPACK_BINARY_DRAGNDROP, ON )
  
  set( CPACK_BUNDLE_NAME "InterSpec.app" )
  set( CPACK_BUNDLE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/../osx/InterSpec.icns" )
  set( CPACK_BUNDLE_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/../osx/build_xcode//CMakeFiles/InterSpecApp.dir/Info.plist" )

  set( RESOURCE_DEST_DIR "InterSpec.app/Contents/Resources" )

  install(TARGETS ${PROJECT_NAME}
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

  if( USE_BATCH_TOOLS AND WIN32 )
    # The CMakeLists.txt for InterSpec_batch will take care of installing that exe,
    #  but we still need to copy over the InterSpec shared lib, since no install target 
    #  is defined for it.
    install( TARGETS InterSpecLib RUNTIME DESTINATION "bin/.." )
  endif( USE_BATCH_TOOLS AND WIN32 )
endif( APPLE )


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
