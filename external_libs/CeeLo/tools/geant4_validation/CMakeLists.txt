cmake_minimum_required(VERSION 3.15)
project(CeeLoG4Validation)

# This CMakeLists.txt is built only when the top-level option
#   CEELO_BUILD_GEANT4_VALIDATION=ON
# is set.  It requires a GEANT4 installation.
#
# To build:
#   cmake -B build -S /path/to/mc_det_eff_plan \
#         -DCEELO_BUILD_GEANT4_VALIDATION=ON \
#         -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4
#   cmake --build build/tools/geant4_validation
#
# Usage:
#   ./ceelo_g4val detector.gdml run.mac output.csv
#
# Workflow for generating reference data (run once, save to repository):
#   1.  Configure a geometry with EfficiencyCalculator, call export_geant4_gdml()
#       and export_geant4_macro() to write detector.gdml and run.mac.
#   2.  Build and run this application to produce output.csv.
#   3.  Copy output.csv to tests/data/geant4_reference/ and commit it.
#   4.  profiling/compare_validation.py then compares CeeLo's results against the
#       stored reference, without needing GEANT4 on CI/CD runners.

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(G4VAL_SOURCES
    src/main.cc
    src/DetectorConstruction.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/EventAction.cc
    src/SteppingAction.cc
)

add_executable(ceelo_g4val ${G4VAL_SOURCES})

target_include_directories(ceelo_g4val PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(ceelo_g4val PRIVATE
    ${Geant4_LIBRARIES}
)

target_compile_features(ceelo_g4val PRIVATE cxx_std_17)

# Copy macros directory to build directory.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/macros
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
