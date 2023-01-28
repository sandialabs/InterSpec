cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project( wxInterSpec )

set( BUILD_AS_WX_WIDGETS_APP ON CACHE BOOL "Set to build wxWidget target" )
set( TRY_TO_STATIC_LINK ON CACHE BOOL "Lets statically link everything" )
set( BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "N/A" )
set( USE_OSX_NATIVE_MENU OFF CACHE BOOL "N/A" )
set( USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "Enable Spectrum File Query Widget" )


IF(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) #Prevent boost auto-linking, which seems to call in vc141 boost libs instead of vc142
  add_definitions(-D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS ) #0x0601==Win7
ENDIF(WIN32)
 

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
 
add_executable(${PROJECT_NAME}  
  WIN32  
  ../electron/windows/InterSpec.rc 
  windows/InterSpec.manifest
  InterSpecWxUtils.h 
  InterSpecWxUtils.cpp
  InterSpecWxApp.h 
  InterSpecWxApp.cpp
  InterSpecWebFrame.h 
  InterSpecWebFrame.cpp
  wxMain.cpp
)


if( InterSpec_FETCH_DEPENDENCIES )
  target_link_libraries(${PROJECT_NAME} PUBLIC wx::core wx::base wx::webview wx::stc ${wxWidgets_LIBRARIES} )
else( InterSpec_FETCH_DEPENDENCIES )
  # Not sure what wxWidgets_USE_FILE actually buys us 
  include( ${wxWidgets_USE_FILE} )

  # Not sure wx::core, wx::base, etc wont work in this case
  target_link_libraries(${PROJECT_NAME} PUBLIC ${wxWidgets_LIBRARIES} )
  target_include_directories( ${PROJECT_NAME} PUBLIC ${wxWidgets_INCLUDE_DIRS} )

  # We need to let the linker know where to find WebView2LoaderStatic.lib
  # TODO: use find file to find this better - for the moment we'll just assume we copied this file to the root lib directory
  target_link_directories( ${PROJECT_NAME} PUBLIC "${wxWidgets_ROOT_DIR}/lib" )
endif( InterSpec_FETCH_DEPENDENCIES )

target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_JS_LIB} InterSpecLib )
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "InterSpec")


set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# We get an error if the destination is empty, or just ".", but the following seems to work
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION "bin/.." )


install(DIRECTORY 
  ${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources
  ${CMAKE_CURRENT_SOURCE_DIR}/../../data
  ${WT_RESOURCES_DIRECTORY}
DESTINATION 
  "bin/.." 
)


install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/../../NOTICE.html
DESTINATION 
  "bin/.." 
)

install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/d3.v3.min.js 
  ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
  ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
  ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SpecUtils/d3_resources/SpectrumChartD3StandAlone.css
DESTINATION 
  "InterSpec_resources"
)

install(FILES 
  ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml
DESTINATION 
  "data"
RENAME 
  "sandia.decay.xml"
)


set(CPACK_BINARY_ZIP ON)
set(CPACK_BINARY_NSIS OFF)
set(CPACK_SOURCE_ZIP ON)
set(CPACK_SOURCE_TGZ ON) #linux

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