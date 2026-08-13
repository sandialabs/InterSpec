cmake_minimum_required(VERSION 3.15)
project(CeeLo VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimized build with debug info when no type is specified (single-
# config generators only).  This keeps timing/validation runs from accidentally
# building -O0 -- pass -DCMAKE_BUILD_TYPE=Debug explicitly for assert-checked
# correctness work.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Build type (Debug/Release/RelWithDebInfo/MinSizeRel)" FORCE)
endif()

# Options
option(CEELO_BUILD_TESTS "Build unit tests" ON)
option(CEELO_BUILD_EXAMPLES "Build examples" ON)
option(CEELO_BUILD_GEANT4_VALIDATION "Build GEANT4 validation (requires GEANT4)" OFF)
option(CEELO_WITH_SANDIADECAY "Build the SandiaDecay cascade adapter" OFF)

# Storable MC-parameterized detector-response payload (io/DetectorResponse,
# io/ResponseGenerator, ...) + its tests and the golden-fixture tool.  Needs
# rapidxml (see CEELO_RAPIDXML_DIR below), which is deliberately NOT vendored --
# so this defaults OFF and a plain `cmake -B build -S .` works with nothing but a
# compiler, CMake, Boost and Eigen.  Opt in with -DCEELO_BUILD_RESPONSE=ON.
option(CEELO_BUILD_RESPONSE
    "Build the storable detector-response payload, generator, and tools" OFF)

# rapidxml (header-only) is NOT vendored: for ODR safety we use the exact copy the
# host application already carries (InterSpec's SpecUtils/3rdparty).  The cache var
# must point at the directory CONTAINING the rapidxml/ folder; when unset we search
# the usual sibling locations for both the vendored and the side-by-side-dev layouts.
set(CEELO_RAPIDXML_DIR "" CACHE PATH
    "Directory containing rapidxml/rapidxml.hpp (e.g. SpecUtils/3rdparty)")
if(CEELO_BUILD_RESPONSE AND NOT CEELO_RAPIDXML_DIR)
    foreach(_cand
        "${CEELO_SANDIADECAY_DIR}/3rdparty"                         # SandiaDecay bundles rapidxml (works from a worktree when -DCEELO_SANDIADECAY_DIR is set)
        "${CMAKE_CURRENT_SOURCE_DIR}/scratch/SandiaDecay/3rdparty"  # dev-only SandiaDecay copy
        "${CMAKE_CURRENT_SOURCE_DIR}/../SpecUtils/3rdparty"
        "${CMAKE_CURRENT_SOURCE_DIR}/../InterSpec_wt4/external_libs/SpecUtils/3rdparty"
        "${CMAKE_CURRENT_SOURCE_DIR}/../InterSpec/external_libs/SpecUtils/3rdparty"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_wt4/external_libs/SpecUtils/3rdparty"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec/external_libs/SpecUtils/3rdparty")
        if(EXISTS "${_cand}/rapidxml/rapidxml.hpp")
            get_filename_component(_cand "${_cand}" ABSOLUTE)
            set(CEELO_RAPIDXML_DIR "${_cand}")
            break()
        endif()
    endforeach()
    if(NOT CEELO_RAPIDXML_DIR)
        message(FATAL_ERROR "CEELO_BUILD_RESPONSE requires rapidxml; could not find "
            "rapidxml/rapidxml.hpp near ${CMAKE_CURRENT_SOURCE_DIR}. Set "
            "CEELO_RAPIDXML_DIR to a directory containing rapidxml/ (e.g. "
            "SpecUtils/3rdparty).")
    endif()
    message(STATUS "CeeLo response payload using rapidxml from ${CEELO_RAPIDXML_DIR}")
endif()

# Where the SandiaDecay source tree lives (cascade-adapter dependency + example decay
# XML).  The core CeeLo library has no SandiaDecay dependency; only CeeLoCascade
# (CEELO_WITH_SANDIADECAY) uses this.  When unset we search the dev-only scratch/
# copy first (the "enriched" tree with cascade a2/a4 + level links), then a
# side-by-side InterSpec checkout; a host that already provides the SandiaDecay
# target (e.g. InterSpec when vendored) skips this entirely.
set(CEELO_SANDIADECAY_DIR "" CACHE PATH
    "Path to the SandiaDecay source tree used by the cascade adapter")
if(CEELO_WITH_SANDIADECAY AND NOT TARGET SandiaDecay AND NOT CEELO_SANDIADECAY_DIR)
    foreach(_cand
        "${CMAKE_CURRENT_SOURCE_DIR}/scratch/SandiaDecay"
        "${CMAKE_CURRENT_SOURCE_DIR}/../InterSpec_wt4/external_libs/SandiaDecay"
        "${CMAKE_CURRENT_SOURCE_DIR}/../InterSpec/external_libs/SandiaDecay"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_wt4/external_libs/SandiaDecay"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec/external_libs/SandiaDecay")
        if(EXISTS "${_cand}/SandiaDecay.h")
            get_filename_component(_cand "${_cand}" ABSOLUTE)
            set(CEELO_SANDIADECAY_DIR "${_cand}")
            break()
        endif()
    endforeach()
    if(NOT CEELO_SANDIADECAY_DIR)
        message(FATAL_ERROR "CEELO_WITH_SANDIADECAY requires a SandiaDecay source "
            "tree; none found near ${CMAKE_CURRENT_SOURCE_DIR}. Set "
            "CEELO_SANDIADECAY_DIR.")
    endif()
    message(STATUS "CeeLo cascade adapter using SandiaDecay from ${CEELO_SANDIADECAY_DIR}")
endif()

# Source-escape multiple-scattering / escape model (compile-time select, for A/B
# comparison of accuracy + speed across the three implementations):
#   gate  = empirical exit-state (1-albedo) skin-escape gate  [default, validated]
#   tails = first-principles: straggling + Highland-B2 + screened-Rutherford tail
#           + double-counting variance subtraction, no empirical gate
#   gs    = Goudsmit-Saunderson: soft-moment Gaussian core (theta0 from the soft
#           transport coefficient of the screened cross-section, replacing Highland)
#           + screened-Rutherford hard tail + straggling, no empirical gate
# Selected at configure time; maps to the CEELO_ESCAPE_MODEL macro (0|1|2).
set(CEELO_SOURCE_ESCAPE_MODEL "gate" CACHE STRING
    "Source-escape model: gate | tails | gs")
set_property(CACHE CEELO_SOURCE_ESCAPE_MODEL PROPERTY STRINGS gate tails gs)
if(CEELO_SOURCE_ESCAPE_MODEL STREQUAL "gate")
    set(_CEELO_ESCAPE_MODEL_NUM 0)
elseif(CEELO_SOURCE_ESCAPE_MODEL STREQUAL "tails")
    set(_CEELO_ESCAPE_MODEL_NUM 1)
elseif(CEELO_SOURCE_ESCAPE_MODEL STREQUAL "gs")
    set(_CEELO_ESCAPE_MODEL_NUM 2)
else()
    message(FATAL_ERROR
        "CEELO_SOURCE_ESCAPE_MODEL must be one of: gate | tails | gs")
endif()
message(STATUS "Source-escape model: ${CEELO_SOURCE_ESCAPE_MODEL} "
               "(CEELO_ESCAPE_MODEL=${_CEELO_ESCAPE_MODEL_NUM})")
add_compile_definitions(CEELO_ESCAPE_MODEL=${_CEELO_ESCAPE_MODEL_NUM})

# Eigen3 is the library's ONLY hard dependency.  Boost is needed solely for
# Boost.Test and is looked up in tests/CMakeLists.txt, which disables the test
# targets rather than failing the configure when Boost is absent.
find_package(Eigen3 REQUIRED)

# Enable native CPU optimizations for SIMD when not cross-compiling.
# Defaults OFF when CeeLo is built as a subproject (e.g. vendored in InterSpec):
# -march=native breaks distributable/universal binaries the host may produce.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_ceelo_native_default ON)
else()
    set(_ceelo_native_default OFF)
endif()
option(CEELO_NATIVE_ARCH "Build with -march=native" ${_ceelo_native_default})
if(CEELO_NATIVE_ARCH AND NOT CMAKE_CROSSCOMPILING)
    if(NOT MSVC)
        add_compile_options(-march=native)
    endif()
endif()

# Optional SandiaDecay adapter: pull in the SandiaDecay static library so the
# adapter (src/cascade/SandiaDecayCascade.cpp) can be built. Must precede src/.
# Guarded: a host project (e.g. InterSpec) may already provide the target.
if(CEELO_WITH_SANDIADECAY AND NOT TARGET SandiaDecay)
    add_subdirectory(${CEELO_SANDIADECAY_DIR} ${CMAKE_BINARY_DIR}/SandiaDecay)
endif()

add_subdirectory(src)

if(CEELO_BUILD_TESTS)
    enable_testing()
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME test_cross_section_provenance
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/tools/prepare_cross_sections/check_provenance.py
                    --require-manual-approval)
    endif()
    add_subdirectory(tests)
endif()

if(CEELO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Golden-fixture generator for host projects' regression tests (e.g.
# InterSpec's regen_ceelo_goldens target).
if(CEELO_BUILD_RESPONSE)
    add_executable(make_golden_response tools/make_golden_response.cpp)
    target_link_libraries(make_golden_response PRIVATE CeeLo)

    # EFFTRAN-style efficiency-transfer demo / validation harness.
    add_executable(efftran_transfer_demo tools/efftran_transfer_demo.cpp)
    target_link_libraries(efftran_transfer_demo PRIVATE CeeLo)
endif()

# Dev-only research/studies (gitignored; absent from the vendored archive, so this block
# is a no-op in a vendored checkout).
if(CEELO_BUILD_EXAMPLES AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/studies/CMakeLists.txt")
    add_subdirectory(studies)
endif()

if(CEELO_BUILD_GEANT4_VALIDATION)
    find_package(Geant4 REQUIRED)
    add_subdirectory(tools/geant4_validation)
endif()
