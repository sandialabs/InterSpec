cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# As the top-level CMakeLists.txt does, but that is only read after project() below, by when the
#  default (the build machine's macOS version) is set; and an executable given to others must run
#  on older macOS.  The dependency prefix (target/dep_build) is built for the same version.
if( APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET )
  set( CMAKE_OSX_DEPLOYMENT_TARGET "13.3" CACHE STRING "macOS deployment target" FORCE )
endif()

project( InterSpecMcpStdio )

# The stdio MCP server: an MCP client (Claude Code, Claude Desktop, ...) launches it, and talks to
#  it over stdin/stdout.  It hosts a single InterSpec session with no browser and no web server,
#  using Wt's libwttest instead of its HTTP connector - so it can not share the InterSpecLib of
#  the desktop app, and builds its own, statically linked into one self-contained executable.
#
# It is a local-server build in all other respects (that is what gives it a writable user-data
#  directory), with no spectrum database (its user database is in memory), and nothing that
#  reaches out of the machine except the LLM API calls of the assistant.

set( BUILD_AS_MCP_STDIO_SERVER ON CACHE BOOL "Build the stdio MCP server" FORCE )
set( BUILD_AS_LOCAL_SERVER ON CACHE BOOL "N/A" FORCE )
set( BUILD_FOR_WEB_DEPLOYMENT OFF CACHE BOOL "N/A" FORCE )
set( BUILD_AS_ELECTRON_APP OFF CACHE BOOL "N/A" FORCE )
set( BUILD_AS_WX_WIDGETS_APP OFF CACHE BOOL "N/A" FORCE )
set( BUILD_AS_MCP_HTTP_SERVER OFF CACHE BOOL "N/A" FORCE )
set( USE_LLM_INTERFACE ON CACHE BOOL "Required by the MCP server" FORCE )
set( USE_NATIVE_HTTP_CLIENT ON CACHE BOOL "There is no browser to make LLM calls from" FORCE )
set( USE_DB_TO_STORE_SPECTRA OFF CACHE BOOL "The user database is in memory" FORCE )
set( USE_REMOTE_RID OFF CACHE BOOL "No outbound calls from tools" FORCE )
set( USE_BATCH_CLI_TOOLS OFF CACHE BOOL "N/A" FORCE )
set( USE_BATCH_GUI_TOOLS OFF CACHE BOOL "N/A" FORCE )
set( USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE BOOL "N/A" FORCE )
set( USE_GOOGLE_MAP OFF CACHE BOOL "No outbound calls" FORCE )
set( USE_LEAFLET_MAP OFF CACHE BOOL "No outbound calls" FORCE )
set( TRY_TO_STATIC_LINK ON CACHE BOOL "One self-contained executable" FORCE )
# Spectrum files no longer displayed are kept in memory up to this size, then written to temporary
#  files; an MCP client rarely goes back to an earlier file, so keep the process small (as on mobile).
set( MAX_SPECTRUM_MEMMORY_SIZE_MB 32 CACHE STRING "In-memory cache of spectrum files" FORCE )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING "Fitting is 10-100x slower in debug builds" FORCE )
endif()

if( MSVC )
  # QuickJS (quickjs-ng) uses C11 atomics, which MSVC has from Visual Studio 2022 17.5.
  if( MSVC_VERSION LESS 1935 )
    message( FATAL_ERROR "The stdio MCP server needs Visual Studio 2022 17.5 or newer, for QuickJS." )
  endif()

  # Static MSVC runtime, so the executable needs no redistributable.
  set( Boost_USE_STATIC_RUNTIME ON )
  set( CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>" )
endif()

if( UNIX AND NOT APPLE )
  # glibc stays dynamic (see below), and so must its libm: the static libm.a needs glibc's static
  #  libc.  TRY_TO_STATIC_LINK makes find_library() prefer `.a` files, so name libm up front for
  #  muparserx's find_library( M_LIBRARY m ).
  set( M_LIBRARY "m" CACHE STRING "The C math library; linked dynamically, as glibc is" FORCE )
endif()

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

# Raster images of the charts (PNG/JPEG, for the image tools) are drawn with plutovg (MIT; plain C,
#  no dependencies), using DejaVu Sans (Bitstream Vera license) compiled into the executable.
include( FetchContent )
FetchContent_Declare( plutovg
  GIT_REPOSITORY https://github.com/sammycage/plutovg.git
  GIT_TAG        bbd91f0d06a71491691b36330f29dffa4af87ccf  # v1.3.3
)
FetchContent_Declare( dejavu_fonts
  URL      https://github.com/dejavu-fonts/dejavu-fonts/releases/download/version_2_37/dejavu-fonts-ttf-2.37.zip
  URL_HASH SHA256=7576310b219e04159d35ff61dd4a4ec4cdba4f35c00e002a136f00e96a908b0a
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# The assistant's run_javascript tool runs in QuickJS (quickjs-ng; MIT): just the engine - not its
#  std/os modules, so the code has no file, network, or process access.
FetchContent_Declare( quickjs
  GIT_REPOSITORY https://github.com/quickjs-ng/quickjs.git
  GIT_TAG        6d46d07d04041b40f4f49eaa7fdebe44c314c699  # v0.17.0
)
set( BUILD_SHARED_LIBS OFF )
set( QJS_BUILD_LIBC OFF CACHE BOOL "" FORCE )
set( QJS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE )
set( QJS_ENABLE_INSTALL OFF CACHE BOOL "" FORCE )
FetchContent_MakeAvailable( plutovg dejavu_fonts quickjs )

set( EMBEDDED_FONT_CPP "${CMAKE_CURRENT_BINARY_DIR}/EmbeddedFont.cpp" )
add_custom_command(
  OUTPUT  "${EMBEDDED_FONT_CPP}"
  COMMAND ${CMAKE_COMMAND} -DINPUT_FILE=${dejavu_fonts_SOURCE_DIR}/ttf/DejaVuSans.ttf
                           -DOUTPUT_FILE=${EMBEDDED_FONT_CPP} -DSYMBOL=g_embedded_font_ttf
                           -P ${CMAKE_CURRENT_SOURCE_DIR}/embed_file.cmake
  DEPENDS "${dejavu_fonts_SOURCE_DIR}/ttf/DejaVuSans.ttf" "${CMAKE_CURRENT_SOURCE_DIR}/embed_file.cmake"
  COMMENT "Embedding DejaVuSans.ttf"
)

add_executable( InterSpecMcpStdio
  mcp_stdio_main.cpp
  PlutoPaintDevice.cpp PlutoPaintDevice.h
  QuickJsSandbox.cpp QuickJsSandbox.h
  "${EMBEDDED_FONT_CPP}"
)
target_link_libraries( InterSpecMcpStdio PUBLIC InterSpecLib plutovg::plutovg qjs )
set_target_properties( InterSpecMcpStdio PROPERTIES
  CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO
  OUTPUT_NAME "InterSpec-mcp-stdio"
)

if( UNIX AND NOT APPLE )
  # glibc stays dynamic (manylinux), but not the C++ runtime.
  target_link_options( InterSpecMcpStdio PRIVATE -static-libstdc++ -static-libgcc )
endif()

# Tests of the QuickJS sandbox: `ctest` in the build directory.
enable_testing()
add_executable( test_quickjs_sandbox tests/test_quickjs_sandbox.cpp QuickJsSandbox.cpp QuickJsSandbox.h )
target_link_libraries( test_quickjs_sandbox PRIVATE InterSpecLib qjs )
target_include_directories( test_quickjs_sandbox PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties( test_quickjs_sandbox PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES )
add_test( NAME quickjs_sandbox COMMAND test_quickjs_sandbox )

# Drives one request with the native HTTP client; run by tests/native_http_test.py.
add_executable( test_native_http tests/test_native_http.cpp )
target_link_libraries( test_native_http PRIVATE InterSpecLib )
set_target_properties( test_native_http PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES )

# `cmake --install` gives a self-contained directory to hand to others: the executable finds the
#  data in ../share/interspec/data, and InterSpec's text next to it (so the app root is
#  share/interspec).  Of InterSpec_resources only the text is needed, as there is no browser.
set( INTERSPEC_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." )
set( SHARE_DEST "share/interspec" )

install( TARGETS InterSpecMcpStdio DESTINATION bin )

# Local KnowledgeBase skills (data/llm_knowledge/*, git-ignored) are the developer's own, so not
#  distributed; sandia.decay.xml is a build-time copy, installed from SandiaDecay below.
install( DIRECTORY "${INTERSPEC_SRC_DIR}/data" DESTINATION "${SHARE_DEST}"
  PATTERN ".DS_Store" EXCLUDE
  REGEX "/data/llm_knowledge/.+" EXCLUDE
  REGEX "/data/gamma_spec_skills" EXCLUDE
  REGEX "/data/sandia\\.decay\\.xml$" EXCLUDE
)
install( FILES "${INTERSPEC_SRC_DIR}/data/llm_knowledge/README.md" DESTINATION "${SHARE_DEST}/data/llm_knowledge" )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/SandiaDecay/sandia.decay.min.xml"
         DESTINATION "${SHARE_DEST}/data" RENAME "sandia.decay.xml" )
install( DIRECTORY "${INTERSPEC_SRC_DIR}/InterSpec_resources/app_text" DESTINATION "${SHARE_DEST}/InterSpec_resources"
         PATTERN ".DS_Store" EXCLUDE )

install( FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md" "${INTERSPEC_SRC_DIR}/LICENSE.txt" "${INTERSPEC_SRC_DIR}/NOTICE.html"
         DESTINATION "${SHARE_DEST}" )

# Licenses of the code compiled into the executable that NOTICE.html does not reproduce in full.
set( LICENSE_DEST "${SHARE_DEST}/licenses" )
install( FILES "${dejavu_fonts_SOURCE_DIR}/LICENSE" DESTINATION "${LICENSE_DEST}" RENAME DejaVuFonts-LICENSE.txt )
install( FILES "${plutovg_SOURCE_DIR}/LICENSE" DESTINATION "${LICENSE_DEST}" RENAME plutovg-LICENSE.txt )
install( FILES "${quickjs_SOURCE_DIR}/LICENSE" DESTINATION "${LICENSE_DEST}" RENAME quickjs-ng-LICENSE.txt )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/SandiaDecay/LICENSE.txt" DESTINATION "${LICENSE_DEST}" RENAME SandiaDecay-LICENSE.txt )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/SpecUtils/LICENSE.txt" DESTINATION "${LICENSE_DEST}" RENAME SpecUtils-LICENSE.txt )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/CeeLo/LICENSE.txt" DESTINATION "${LICENSE_DEST}" RENAME CeeLo-LICENSE.txt )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/CeeLo/THIRD-PARTY-NOTICES.md" DESTINATION "${LICENSE_DEST}" RENAME CeeLo-THIRD-PARTY-NOTICES.md )
install( FILES "${INTERSPEC_SRC_DIR}/external_libs/muparserx-4.0.7/License.txt" DESTINATION "${LICENSE_DEST}" RENAME muparserx-LICENSE.txt )
set( SPECUTILS_3RDPARTY "${INTERSPEC_SRC_DIR}/external_libs/SpecUtils/3rdparty" )
install( FILES "${SPECUTILS_3RDPARTY}/nlohmann/LICENSE.MIT" DESTINATION "${LICENSE_DEST}" RENAME nlohmann_json-LICENSE.txt )
install( FILES "${SPECUTILS_3RDPARTY}/inja/LICENSE.txt" DESTINATION "${LICENSE_DEST}" RENAME inja-LICENSE.txt )
install( FILES "${SPECUTILS_3RDPARTY}/date/LICENSE.txt" DESTINATION "${LICENSE_DEST}" RENAME date-LICENSE.txt )
install( FILES "${SPECUTILS_3RDPARTY}/rapidxml/license.txt" DESTINATION "${LICENSE_DEST}" RENAME rapidxml-LICENSE.txt )
install( FILES "${SPECUTILS_3RDPARTY}/code_from_boost/LICENSE_1_0.txt" DESTINATION "${LICENSE_DEST}" RENAME Boost-LICENSE_1_0.txt )
