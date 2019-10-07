cmake_policy(SET CMP0048 NEW)
project( InterSpecAddOn )

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

if( BOOST_ROOT AND Wt_INCLUDE_DIR )
elseif( APPLE )
  set( BOOST_ROOT "/Users/wcjohns/install/osx_wt_prefix/" CACHE BOOL "Set boost root dir to use" )
  set( Wt_INCLUDE_DIR "/Users/wcjohns/install/osx_wt_prefix/include/" CACHE BOOL "Set Wt include dir to use" )
elseif(WIN32)
endif()

set( BUILD_AS_ELECTRON_APP ON CACHE BOOL "Set to build Electron target" )
set( TRY_TO_STATIC_LINK ON CACHE BOOL "Lets statically link everything" )


IF(WIN32)
  set( MSVC_RUNTIME "static" )
  include( ../../cmake/ConfigureMsvc.txt )
  configure_msvc_runtime() 
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj -D_SCL_SECURE_NO_WARNINGS /MP /wd4996 /wd4267 /DWINVER=0x0601 /D_WIN32_WINNT=0x0601" )
    #0x0601==Win7, 0x0501==WinXP
ENDIF(WIN32)


add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )


include_directories( ${CMAKE_JS_INC})
set( SOURCE_FILES InterSpecAddOn.cpp )
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${CMAKE_JS_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "" SUFFIX ".node")
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_JS_LIB} InterSpecLib )

target_include_directories( ${PROJECT_NAME} PUBLIC 
                            node_modules/node-addon-api/ 
                            ${CMAKE_CURRENT_SOURCE_DIR}/../.. 
)


# Need to put all this stuff below into the "install" stage.
# Also need to copy InterSpecAddOn.node to the app directory
# Also need to fix up FindWt.cmake in SpecUtils.
# Then need to test out running InterSpec and packaging it, and updating the README.md


#CMAKE_INSTALL_PREFIX

install( DIRECTORY ${WT_RESOURCES_DIRECTORY} DESTINATION ${PROJECT_BINARY_DIR}/app/ )

install( DIRECTORY ${WT_RESOURCES_DIRECTORY} DESTINATION ${PROJECT_BINARY_DIR}/app/ )
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../data 
         DESTINATION ${PROJECT_BINARY_DIR}/app/ 
         PATTERN "*OUO*" EXCLUDE
         PATTERN "*ouo*" EXCLUDE
)
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources DESTINATION ${PROJECT_BINARY_DIR}/app/ )
install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../example_spectra DESTINATION ${PROJECT_BINARY_DIR}/app/ )
install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/app/loading.html
        ${CMAKE_CURRENT_SOURCE_DIR}/app/main.js
        ${CMAKE_CURRENT_SOURCE_DIR}/app/launch_options.json
        DESTINATION ${PROJECT_BINARY_DIR}/app/ 
)
#Need to get d3.min.js and related into InterSpec_resources as well!

install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
        DESTINATION ${PROJECT_BINARY_DIR}/app/data/
        RENAME sandia.decay.xml )
install( FILES ${CMAKE_CURRENT_SOURCE_DIR}/package_for_running.json
        DESTINATION ${PROJECT_BINARY_DIR}/app/
        RENAME package.json )
install( TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_BINARY_DIR}/app/ )



