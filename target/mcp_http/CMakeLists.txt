cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project( InterSpecMcpHttp )

# The multi-session MCP analysis server: MCP clients create sessions, exchange files, and drive
#  InterSpec (and its LLM assistant) over HTTP; each session's GUI lives in a tab of a headless
#  Chromium.  See README.md.
#
# It is a local-server build (the sessions are served to that browser over loopback), so the MCP
#  server is compiled in; anything that reaches out of the machine, or lets a session see the
#  server's file system, is left out.

set( BUILD_AS_MCP_HTTP_SERVER ON CACHE BOOL "Build the multi-session MCP HTTP analysis server" FORCE )
set( BUILD_AS_LOCAL_SERVER ON CACHE BOOL "N/A" FORCE )
set( BUILD_FOR_WEB_DEPLOYMENT OFF CACHE BOOL "N/A" FORCE )
set( BUILD_AS_ELECTRON_APP OFF CACHE BOOL "N/A" FORCE )
set( BUILD_AS_WX_WIDGETS_APP OFF CACHE BOOL "N/A" FORCE )
set( USE_LLM_INTERFACE ON CACHE BOOL "Required by the MCP server" FORCE )
set( USE_NATIVE_HTTP_CLIENT ON CACHE BOOL "LLM calls are made from C++, never the browser" FORCE )
set( USE_DB_TO_STORE_SPECTRA OFF CACHE BOOL "Sessions keep their state in exported N42 files" FORCE )
set( USE_REMOTE_RID OFF CACHE BOOL "No outbound calls from tools" FORCE )
set( USE_BATCH_CLI_TOOLS OFF CACHE BOOL "N/A" FORCE )
set( USE_BATCH_GUI_TOOLS OFF CACHE BOOL "N/A" FORCE )
set( USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE BOOL "Would browse the server's file system" FORCE )
set( USE_GOOGLE_MAP OFF CACHE BOOL "No outbound calls" FORCE )
set( USE_LEAFLET_MAP OFF CACHE BOOL "No outbound calls" FORCE )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING "Fitting is 10-100x slower in debug builds" FORCE )
endif()

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

# The agent-server pieces InterSpecLib calls into (it is built with BUILD_AS_MCP_HTTP_SERVER); like
#  target/electron, they are compiled here, and found as "target/mcp_http/<file>.h".
add_executable( InterSpecMcpServer
  mcp_http_main.cpp
  McpSessionManager.cpp McpSessionManager.h
  McpAgentTools.cpp McpAgentTools.h
)
target_link_libraries( InterSpecMcpServer PUBLIC InterSpecLib )
set_target_properties( InterSpecMcpServer PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO )

install( TARGETS InterSpecMcpServer DESTINATION bin )
install( FILES wt_config_mcp_http.xml DESTINATION share/interspec/data/config )
