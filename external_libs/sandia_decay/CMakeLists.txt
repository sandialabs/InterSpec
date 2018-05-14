cmake_minimum_required(VERSION 3.0)

#This file defines one useful variables in its parent scope:
#  SandiaDecay_INCLUDE_DIRS: The list of directories that should be in the include
#                          path of the parent project to this library (e.g., path to SandiaDecay.h)


project(SandiaDecay)

set( LIBRARY_OUTPUT_PATH lib CACHE PATH "Path to libraries" FORCE )

if(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_path( RAPID_XML_DIR rapidxml.hpp ${CMAKE_CURRENT_SOURCE_DIR}/../rapidxml ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${RAPID_XML_DIR} )

set( sources ${CMAKE_CURRENT_SOURCE_DIR}/SandiaDecay.cpp )
set( headers ${CMAKE_CURRENT_SOURCE_DIR}/SandiaDecay.h )

add_library( SandiaDecay STATIC ${sources} ${headers} ${OTHER_SUPPORT_FILES} )

get_directory_property(hasParent PARENT_DIRECTORY)
if( hasParent  )
  set( SandiaDecay_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE )
else( hasParent )
  add_subdirectory( tools )
  add_subdirectory( examples )
endif( hasParent )


