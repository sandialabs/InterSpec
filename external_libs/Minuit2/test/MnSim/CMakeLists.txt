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

set ( GAUS_SRCS
  GaussDataGen.cxx
  GaussFcn.cxx
  GaussFcn2.cxx
)

add_library ( gaus
  ${GAUS_SRCS}
)




find_library( MINUIT2_LIBRARY minuit2 ../../build/lib )


add_executable( PaulTest PaulTest.cxx )

target_link_libraries ( PaulTest ${MINUIT2_LIBRARY} gaus )

