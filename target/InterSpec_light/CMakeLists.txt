cmake_minimum_required( VERSION 3.20 )
project( InterSpecLight CXX C )

# One project, two configurations:
#  - Emscripten (emcmake): builds the `light_wasm` module (.js + .wasm), which assemble_html.py
#    embeds in the single-file HTML.  Dependencies come from deps/build_deps.sh.
#  - Native: builds `refgen` (the reference-line generator; needs only this repository), and, if
#    Ceres is found (e.g., CMAKE_PREFIX_PATH set to InterSpec's dependency prefix), `light_cli` (the
#    test and debug driver).

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

option( LIGHT_PEAK_FILE_IO "Read and write peaks in N42-2012 and IAEA SPE files, compatibly with InterSpec" ON )

set( REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." )
get_filename_component( REPO_ROOT "${REPO_ROOT}" ABSOLUTE )

if( EMSCRIPTEN )
  set( LIGHT_DEPS_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/deps/prefix" CACHE PATH "Emscripten-built deps (deps/build_deps.sh)" )
  # Must match deps/build_deps.sh
  set( LIGHT_WASM_FLAGS "-fwasm-exceptions -sWASM_LEGACY_EXCEPTIONS=1" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIGHT_WASM_FLAGS}" )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LIGHT_WASM_FLAGS}" )
  # The page embeds the module, so build type MinSizeRel optimizes for size (-Oz: ~30% smaller
  #  than -O3, with no measurable slowdown in peak fitting).
  set( CMAKE_CXX_FLAGS_MINSIZEREL "-Oz -DNDEBUG" CACHE STRING "" FORCE )
  set( CMAKE_C_FLAGS_MINSIZEREL "-Oz -DNDEBUG" CACHE STRING "" FORCE )
  set( CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "-Oz" CACHE STRING "" FORCE )
  list( APPEND CMAKE_FIND_ROOT_PATH "${LIGHT_DEPS_PREFIX}" )
  set( absl_DIR "${LIGHT_DEPS_PREFIX}/lib/cmake/absl" CACHE PATH "" )
  file( GLOB_RECURSE _eigen_cfg "${LIGHT_DEPS_PREFIX}/share/*/Eigen3Config.cmake" )
  list( GET _eigen_cfg 0 _eigen_cfg )
  get_filename_component( _eigen_dir "${_eigen_cfg}" DIRECTORY )
  set( Eigen3_DIR "${_eigen_dir}" CACHE PATH "" )
  set( Ceres_DIR "${LIGHT_DEPS_PREFIX}/lib/cmake/Ceres" CACHE PATH "" )
  # Keep this machine's paths out of the WASM (e.g., __FILE__ in log and exception messages)
  add_compile_options( "-ffile-prefix-map=${REPO_ROOT}/=InterSpec/" )
  find_package( Ceres REQUIRED )
  set( LIGHT_BUILD_APP ON )
else()
  find_package( Ceres QUIET )
  if( Ceres_FOUND )
    set( LIGHT_BUILD_APP ON )
  else()
    message( STATUS "Ceres not found (set CMAKE_PREFIX_PATH to InterSpec's dependency prefix):"
                    " light_cli, the test driver, will not be built" )
    set( LIGHT_BUILD_APP OFF )
  endif()
endif()

# SpecUtils, configured like its own wasm binding (no threads, no zlib), plus the D3 JSON export.
set( SpecUtils_USING_NO_THREADING ON CACHE BOOL "" FORCE )
set( SpecUtils_C_BINDINGS OFF CACHE BOOL "" FORCE )
set( SpecUtils_ENABLE_D3_CHART ON CACHE BOOL "" FORCE )
set( SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE BOOL "" FORCE )
# Only used by SpecUtils' D3 HTML export (which the page doesnt offer); its default, SpecUtils' source
#  directory, would put this machine's path into the WASM
set( SpecUtils_D3_SCRIPTS_RUNTIME_DIR "." CACHE STRING "" FORCE )
set( SpecUtils_ENABLE_URI_SPECTRA OFF CACHE BOOL "" FORCE )
set( SpecUtils_FLT_PARSE_METHOD "strtod" CACHE STRING "" FORCE )
set( SpecUtils_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE )
set( SpecUtils_BUILD_FUZZING_TESTS OFF CACHE BOOL "" FORCE )
set( SpecUtils_BUILD_REGRESSION_TEST OFF CACHE BOOL "" FORCE )
set( SpecUtils_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE )
set( SpecUtils_JAVA_SWIG OFF CACHE BOOL "" FORCE )
set( SpecUtils_INJA_TEMPLATES OFF CACHE BOOL "" FORCE )
set( SpecUtils_USE_SIMD OFF CACHE BOOL "" FORCE )
set( SpecUtils_ENABLE_EQUALITY_CHECKS OFF CACHE BOOL "" FORCE )
set( PERFORM_DEVELOPER_CHECKS OFF CACHE BOOL "" FORCE )
add_subdirectory( "${REPO_ROOT}/external_libs/SpecUtils" "${CMAKE_BINARY_DIR}/SpecUtils" )


# Hard-forked InterSpec peak fitting / energy calibration code.
set( FORK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cpp/fork" )
set( APP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cpp/app" )

if( LIGHT_BUILD_APP )
  add_library( light_fork STATIC
    ${FORK_DIR}/src/PeakDef.cpp
    ${FORK_DIR}/src/PeakDists.cpp
    ${FORK_DIR}/src/PeakFit.cpp
    ${FORK_DIR}/src/PeakFitLM.cpp
    ${FORK_DIR}/src/PeakFitSpecImp.cpp
    ${FORK_DIR}/src/PeakFitUtils.cpp
    ${FORK_DIR}/src/PeakFitDetPrefs.cpp
    ${FORK_DIR}/src/EnergyCal.cpp
    ${FORK_DIR}/src/LightMath.cpp
    ${FORK_DIR}/src/PhysicalUnits.cpp
    $<$<BOOL:${LIGHT_PEAK_FILE_IO}>:${FORK_DIR}/src/PeakDefXml.cpp>
  )
  # The repo root is deliberately NOT an include directory, so "InterSpec/..." always
  # resolves to the fork.  external_libs/ provides VoigtDistribution/ and SandiaDecay/.
  target_include_directories( light_fork PUBLIC
    "${FORK_DIR}"
    "${REPO_ROOT}/external_libs"
    "${REPO_ROOT}/external_libs/SpecUtils/3rdparty"
  )
  target_compile_definitions( light_fork PUBLIC LIGHT_PEAK_FILE_IO=$<BOOL:${LIGHT_PEAK_FILE_IO}> )
  target_link_libraries( light_fork PUBLIC SpecUtils Ceres::ceres )
  if( EMSCRIPTEN )
    # The translation units with the hot fitting code (cost functions, peak shapes) stay at -O3 even
    #  when the rest is built for size; all -Oz is ~3x slower to fit.
    set_source_files_properties( ${FORK_DIR}/src/PeakFitLM.cpp ${FORK_DIR}/src/PeakDists.cpp ${FORK_DIR}/src/PeakFit.cpp
                                 PROPERTIES COMPILE_OPTIONS "-O3" )
  endif()


  # The app: session state + JSON dispatcher around the forked code
  add_library( light_app STATIC
    ${APP_DIR}/RefLib.cpp
    ${APP_DIR}/PeakCsv.cpp
    ${APP_DIR}/Session.cpp
    ${APP_DIR}/SessionPeaks.cpp
    ${APP_DIR}/SessionEnergyCal.cpp
    $<$<BOOL:${LIGHT_PEAK_FILE_IO}>:${APP_DIR}/PeakFileIo.cpp>
  )
  target_include_directories( light_app PUBLIC "${APP_DIR}" )
  target_link_libraries( light_app PUBLIC light_fork )


  if( EMSCRIPTEN )
    # LightApi.cpp is compiled into each executable (not light_app), so light_call is always linked
    add_executable( light_wasm ${APP_DIR}/wasm_main.cpp ${APP_DIR}/LightApi.cpp )
    target_link_libraries( light_wasm PRIVATE light_app )
    set_target_properties( light_wasm PROPERTIES SUFFIX ".js" )
    target_link_options( light_wasm PRIVATE
      -fwasm-exceptions
      -sWASM_LEGACY_EXCEPTIONS=1
      -sSINGLE_FILE=0
      -sMODULARIZE=1
      -sEXPORT_NAME=InterSpecLightModule
      -sALLOW_MEMORY_GROWTH=1
      -sINITIAL_MEMORY=67108864
      -sSTACK_SIZE=8388608
      -sFORCE_FILESYSTEM=1
      -sENVIRONMENT=web,node
      -sEXPORTED_FUNCTIONS=_light_call,_main
      -sEXPORTED_RUNTIME_METHODS=cwrap,FS,UTF8ToString
    )
  else()
    add_executable( light_cli cpp/tools/light_cli.cpp ${APP_DIR}/LightApi.cpp )
    target_link_libraries( light_cli PRIVATE light_app )
  endif()
endif( LIGHT_BUILD_APP )


if( NOT EMSCRIPTEN )
  # Reference-line library generator (needs the full decay data, so only built natively)
  add_executable( refgen
    refgen/make_ref_lines.cpp
    "${REPO_ROOT}/external_libs/SandiaDecay/SandiaDecay.cpp"
    ${FORK_DIR}/src/ReactionGamma.cpp
  )
  target_include_directories( refgen PRIVATE "${FORK_DIR}" "${REPO_ROOT}/external_libs"
                              "${REPO_ROOT}/external_libs/SpecUtils/3rdparty" "${CMAKE_CURRENT_SOURCE_DIR}/refgen" )
  target_link_libraries( refgen PRIVATE SpecUtils )
endif()
