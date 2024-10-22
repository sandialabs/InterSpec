# Require CMake 3.15 to use MSVC_RUNTIME_LIBRARY
cmake_minimum_required(VERSION 3.1...3.20 FATAL_ERROR)
cmake_policy(SET CMP0091 NEW)
cmake_policy(SET CMP0042 NEW)

project( InterSpecAddOn )

add_definitions(-DNAPI_VERSION=3) #node.js v6.14.2*, 8.11.2+, v9.11.0+*, 10.0.0 and all later versions

include_directories(${CMAKE_JS_INC})

set( BUILD_AS_ELECTRON_APP ON CACHE BOOL "Set to build Electron target" )
set( TRY_TO_STATIC_LINK ON CACHE BOOL "Lets statically link everything" )
set( BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "N/A" )
set( USE_OSX_NATIVE_MENU OFF CACHE BOOL "N/A" )
set( USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "Enable Spectrum File Query Widget" )
set( USE_REMOTE_RID ON CACHE BOOL "Enables using remote RID tool" )
set( USE_BATCH_TOOLS ON CACHE BOOL "Enables using batch tool" )

IF(WIN32)
  add_definitions(-DBOOST_ALL_NO_LIB) #Prevent boost auto-linking, which seems to call in vc141 boost libs instead of vc142
  add_definitions(-D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS ) #0x0601==Win7
ENDIF(WIN32)

if(MSVC)
  option(${PROJECT_NAME}_USE_MSVC_MultiThreadDLL "Use dynamically-link runtime library." OFF)

  if( ${PROJECT_NAME}_USE_MSVC_MultiThreadDLL)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
  message(STATUS "Using CMAKE_MSVC_RUNTIME_LIBRARY ${CMAKE_MSVC_RUNTIME_LIBRARY}")
endif()

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )


include_directories( ${CMAKE_JS_INC})
set( SOURCE_FILES InterSpecAddOn.cpp InterSpecAddOn.h ElectronUtils.cpp ElectronUtils.h ../../js/AppHtmlMenu.js )
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${CMAKE_JS_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "" SUFFIX ".node")
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_JS_LIB} InterSpecLib )
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )


# If we are are running cmake from a random directory, and not the same place `node_modules` is at,
#   we should search aroung a bit.
#  We'll allow user to specify `InterSpec_NODE_MODULE_DIR` or set the environment variable 
#  `NODE_MODULES_PARENT_DIR` to help with this
#IF(WIN32)
#  execute_process( COMMAND cmd /C cd OUTPUT_VARIABLE CWD OUTPUT_STRIP_TRAILING_WHITESPACE )
#ELSE(WIN32)
#  execute_process( COMMAND pwd OUTPUT_VARIABLE CWD OUTPUT_STRIP_TRAILING_WHITESPACE )
#ENDIF(WIN32)

#find_file( NODE_ADDON_API_DIR napi.h 
#  PATHS 
#    node_modules/node-addon-api
#    ${CMAKE_CURRENT_SOURCE_DIR}/node_modules/node-addon-api
#    ${InterSpec_NODE_MODULE_DIR}/node_modules/node-addon-api
#    ${InterSpec_NODE_MODULE_DIR}/node-addon-api
#    ${InterSpec_NODE_MODULE_DIR}/..
#    ./node_modules/node-addon-api
#    ../node_modules/node-addon-api
#    ${CWD}
#    ${CWD}/node_modules/node-addon-api
#    $ENV{NODE_MODULES_PARENT_DIR}/node_modules/node-addon-api
#    /build_working_dir/node_modules/node-addon-api
#  PATH_SUFFIXES include
#)

#message( "Using node_modules directory at ${NODE_ADDON_API_DIR}")

#if(NOT NODE_ADDON_API_DIR)
#  message(FATAL_ERROR "Failed to find node_modules directory.")
#endif()

#get_filename_component(NODE_ADDON_API_DIR "${NODE_ADDON_API_DIR}" DIRECTORY)
#get_filename_component(NODE_MODULES_DIR "${NODE_ADDON_API_DIR}" DIRECTORY)

#message( "Using node_modules directory at NODE_ADDON_API_DIR=${NODE_ADDON_API_DIR} and NODE_MODULES_DIR=${NODE_MODULES_DIR}" )

#target_include_directories( ${PROJECT_NAME} PUBLIC 
#                            ${NODE_ADDON_API_DIR}
#                            ${NODE_MODULES_DIR}/node-api-headers/
#                            ${NODE_MODULES_DIR}/node-api-headers/include
#                            ${CMAKE_CURRENT_SOURCE_DIR}/../.. 
#)

if(MSVC AND CMAKE_JS_NODELIB_DEF AND CMAKE_JS_NODELIB_TARGET)
  # Generate node.lib
  execute_process(COMMAND ${CMAKE_AR} /def:${CMAKE_JS_NODELIB_DEF} /out:${CMAKE_JS_NODELIB_TARGET} ${CMAKE_STATIC_LINKER_FLAGS})
endif()

set( ELECTRON_SUPPORT_FILES
  app/main.js
  app/loading.html
  README.md
  package.json
)

add_custom_target( electron_resources SOURCES ${ELECTRON_SUPPORT_FILES} )


install( DIRECTORY ${WT_RESOURCES_DIRECTORY} DESTINATION ${PROJECT_BINARY_DIR}/app/ )
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../data 
         DESTINATION ${PROJECT_BINARY_DIR}/app/ 
         PATTERN "*OUO*" EXCLUDE
         PATTERN "*ouo*" EXCLUDE
)
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources DESTINATION ${PROJECT_BINARY_DIR}/app/ )

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec/InterSpec_resources/d3.v3.min.js
        ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec/InterSpec_resources/SpectrumChartD3.css
        ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec/InterSpec_resources/SpectrumChartD3.js
         DESTINATION ${PROJECT_BINARY_DIR}/app/InterSpec_resources/ )
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../example_spectra DESTINATION ${PROJECT_BINARY_DIR}/app/ )
install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/app/loading.html
        ${CMAKE_CURRENT_SOURCE_DIR}/app/main.js
        DESTINATION ${PROJECT_BINARY_DIR}/app/ 
)

#Need to get {d3.v3.min.js, SpectrumChartD3.js, SpectrumChartD3.css} and related into InterSpec_resources as well.
#  SHould be workin, but on Windows it doesnt look like it at the moment.

install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SandiaDecay/sandia.decay.min.xml
        # PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
        DESTINATION ${PROJECT_BINARY_DIR}/app/data/
        RENAME sandia.decay.xml )
install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/package_for_running.json
        DESTINATION ${PROJECT_BINARY_DIR}/app/
        RENAME package.json )
install( TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_BINARY_DIR}/app/ )



