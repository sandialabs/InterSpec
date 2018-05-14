cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayTools )

find_path( RAPID_XML_DIR rapidxml.hpp ${CMAKE_CURRENT_SOURCE_DIR}/../rapidxml ${CMAKE_CURRENT_SOURCE_DIR}/../../rapidxml ${CMAKE_CURRENT_SOURCE_DIR}/rapidxml )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${RAPID_XML_DIR} )

add_executable( remove_edit_history remove_edit_history.cpp )

include( cmake/FindSQLite3.cmake )
include_directories( ${SQLITE3_INCLUDE_DIRS} )
add_executable( add_tori_xray_to_xml_and_test add_tori_xray_to_xml_and_test.cpp )
target_link_libraries( add_tori_xray_to_xml_and_test SandiaDecay ${SQLITE3_LIBRARIES} )

add_executable( decay_database_benchmark decay_database_benchmark.cpp )
target_link_libraries( decay_database_benchmark SandiaDecay )

