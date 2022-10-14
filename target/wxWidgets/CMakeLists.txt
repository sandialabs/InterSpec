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
endif()

set(wxUSE_WEBVIEW_EDGE ON CACHE BOOL "")
set(wxUSE_WEBVIEW_EDGE_STATIC ON CACHE BOOL "")
set(wxBUILD_USE_STATIC_RUNTIME ON CACHE BOOL "")
set(wxBUILD_SHARED OFF CACHE BOOL "")

add_subdirectory(3rd_party/wxWidgets-3.2.1)

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

add_definitions( -DwxUSE_RC_MANIFEST=1 -DwxUSE_DPI_AWARE_MANIFEST=2 )
 
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


target_link_libraries(${PROJECT_NAME} PUBLIC wx::core wx::base wx::webview wx::stc ${wxWidgets_LIBRARIES} )
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_JS_LIB} InterSpecLib )
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "InterSpec")


