cmake_minimum_required(VERSION 3.15...3.21)

include( ../../cmake/FetchInterSpecDeps.txt )

project( InterSpecAndroid VERSION 1.10 )

add_executable( example main.cpp)
target_link_libraries( example PUBLIC Boost::filesystem wt )
target_include_directories( example PUBLIC ${wt_BINARY_DIR} )
