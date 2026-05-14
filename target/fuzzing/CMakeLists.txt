# Standalone CMake project for InterSpec libFuzzer harnesses.  Builds a copy of
# InterSpecLib with ASan+UBSan+coverage instrumentation, then links each fuzz
# binary against it with the libFuzzer driver.
#
# Why standalone (not part of target/testing/): the sanitizer flags applied to
# InterSpecLib here would prevent the unit-test executables from linking, since
# they don't get rebuilt with the same flags.  Keeping this as its own project
# avoids the conflict.
#
# See README.md for build & run instructions.

cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project( InterSpecFuzzers )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

if( NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  message( FATAL_ERROR
    "Fuzzers require Clang/LLVM (Xcode's clang lacks the fuzzer runtime). "
    "Re-configure with:\n"
    "  cmake -DCMAKE_C_COMPILER=$(brew --prefix llvm)/bin/clang"
    " -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ ..."
  )
endif()


# Match the test build's expectations - InterSpecLib needs the same set of CMake
# options the unit-test build sets, otherwise feature flags drift.
set( BUILD_AS_UNIT_TEST_SUITE ON CACHE BOOL "" )
set( BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "" )
set( USE_DB_TO_STORE_SPECTRA ON CACHE BOOL "" )
set( PERFORM_DEVELOPER_CHECKS OFF CACHE BOOL "" )
set( SpecUtils_ENABLE_EQUALITY_CHECKS OFF CACHE BOOL "" )
set( USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "" )
set( USE_BATCH_CLI_TOOLS OFF CACHE INTERNAL "" )
set( USE_BATCH_GUI_TOOLS ON CACHE INTERNAL "" )

# Pull the InterSpec sources in as a subproject - this builds InterSpecLib here.
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )


# Sanitizer flags.  The fuzz binary links libFuzzer's driver via -fsanitize=fuzzer;
# the InterSpecLib library only needs coverage instrumentation (-fsanitize=fuzzer-no-link)
# so the symbol "main" comes from libFuzzer alone.
#
# -DNDEBUG: fuzz in release-asserts mode.  We want to surface real bugs that
# affect end users (ASan/UBSan still fire), not the many `assert()` calls used
# as developer round-trip-correctness checks - those compile out of release
# builds and aren't real safety issues.
set( FUZZ_INSTRUMENT_FLAGS
  -g -O1
  -DNDEBUG
  -fsanitize=fuzzer-no-link,address,undefined
  -fno-sanitize-recover=all
  -fno-omit-frame-pointer
)
set( FUZZ_DRIVER_LINK_FLAGS
  -fsanitize=fuzzer,address,undefined
)

target_compile_options( InterSpecLib PRIVATE ${FUZZ_INSTRUMENT_FLAGS} )
target_link_options( InterSpecLib PRIVATE
  -fsanitize=fuzzer-no-link,address,undefined
)


function( add_fuzz_target name )
  add_executable( ${name} ${name}.cpp )
  target_compile_options( ${name} PRIVATE ${FUZZ_INSTRUMENT_FLAGS} )
  target_link_options( ${name} PRIVATE ${FUZZ_DRIVER_LINK_FLAGS} )
  target_link_libraries( ${name} PRIVATE InterSpecLib )
  set_target_properties( ${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endfunction()


add_fuzz_target( fuzz_app_url )
add_fuzz_target( fuzz_color_theme )
add_fuzz_target( fuzz_gadras_detector_dat )
add_fuzz_target( fuzz_gadras_efficiency_csv )
add_fuzz_target( fuzz_rel_act_csv )
add_fuzz_target( fuzz_rel_eff_drf_csv )
add_fuzz_target( fuzz_state_xml )
