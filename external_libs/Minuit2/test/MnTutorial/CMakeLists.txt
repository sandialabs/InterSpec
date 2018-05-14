project(minuit2)


#-----------------------------------------------------------------------------
# CMake administrivia
cmake_minimum_required( VERSION 2.8 )
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND cmake_policy )


set( LIBRARY_OUTPUT_PATH lib CACHE PATH "Path to libraries" FORCE )
set( EXECUTABLE_OUTPUT_PATH bin CACHE PATH "Path to executables" FORCE )

set( UPDATE_TYPE "svn" )
set( CMAKE_COLOR_MAKEFILE ON CACHE BOOL "Color" FORCE )


#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_FILESYSTEM_VERSION=3")

#-----------------------------------------------------------------------------

#SET(CMAKE_MODULE_PATH "${WT_HISTOGRAM_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
#MESSAGE( ${WT_HISTOGRAM_SOURCE_DIR} )



INCLUDE_DIRECTORIES(
  ../../inc
  .
  ..
)

find_library( MINUIT2_LIBRARY minuit2 ../../build/lib )


add_executable( Quad12FMain Quad12FMain.cxx )
add_executable( Quad1FMain Quad1FMain.cxx )
add_executable( Quad4FMain Quad4FMain.cxx )
add_executable( Quad8FMain Quad8FMain.cxx )

target_link_libraries ( Quad12FMain ${MINUIT2_LIBRARY} )
target_link_libraries ( Quad1FMain ${MINUIT2_LIBRARY} )
target_link_libraries ( Quad4FMain ${MINUIT2_LIBRARY} )
target_link_libraries ( Quad8FMain ${MINUIT2_LIBRARY} )

