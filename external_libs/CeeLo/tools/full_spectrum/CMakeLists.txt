add_executable(ceelo_source_escape_campaign source_escape_campaign.cpp)
target_link_libraries(ceelo_source_escape_campaign PRIVATE CeeLo)

add_executable(ceelo_detector_response_campaign detector_response_campaign.cpp)
target_link_libraries(ceelo_detector_response_campaign PRIVATE CeeLo)

add_executable(ceelo_source_model_evaluate source_model_evaluate.cpp)
target_link_libraries(ceelo_source_model_evaluate PRIVATE CeeLo)

add_executable(
    ceelo_source_cross_section_features source_cross_section_features.cpp)
target_link_libraries(ceelo_source_cross_section_features PRIVATE CeeLo)

add_executable(ceelo_source_electron_campaign source_electron_campaign.cpp)
target_link_libraries(ceelo_source_electron_campaign PRIVATE CeeLo)

add_executable(
    ceelo_charged_particle_model_evaluate
    charged_particle_model_evaluate.cpp)
target_link_libraries(
    ceelo_charged_particle_model_evaluate PRIVATE CeeLo)

add_executable(
    ceelo_detector_model_evaluate
    detector_model_evaluate.cpp)
target_link_libraries(
    ceelo_detector_model_evaluate PRIVATE CeeLo)

# The detector payload's finite-distance depths, by point-source Monte Carlo.
add_executable(
    ceelo_detector_finite_distance
    detector_finite_distance.cpp)
target_link_libraries(
    ceelo_detector_finite_distance PRIVATE CeeLo)

add_executable(
    bench_detector_response_model
    bench_detector_response_model.cpp)
target_link_libraries(
    bench_detector_response_model PRIVATE CeeLo)

add_executable(
    bench_charged_particle_source
    bench_charged_particle_source.cpp)
target_link_libraries(
    bench_charged_particle_source PRIVATE CeeLo)

add_executable(
    ceelo_rate_surface_anchor
    rate_surface_anchor.cpp)
target_link_libraries(
    ceelo_rate_surface_anchor PRIVATE CeeLo)

# GEANT4 programs for the optional GEANT4 parts of the source-escape bake: the
# natural-element photonuclear attenuation table, and the photonuclear and
# acceptance-grid reference cells. Built only when GEANT4 is found (source
# geant4.sh before configuring).
find_package(Geant4 QUIET)
if(Geant4_FOUND)
    include(${Geant4_USE_FILE})
    add_executable(g4_source_escape_campaign g4_source_escape_campaign.cc)
    target_link_libraries(g4_source_escape_campaign PRIVATE
        CeeLo ${Geant4_LIBRARIES})
    target_compile_features(g4_source_escape_campaign PRIVATE cxx_std_17)

    add_executable(g4_photonuclear_xs g4_photonuclear_xs.cc)
    target_link_libraries(g4_photonuclear_xs PRIVATE
        CeeLo ${Geant4_LIBRARIES})
    target_compile_features(g4_photonuclear_xs PRIVATE cxx_std_17)
else()
    message(STATUS "GEANT4 not found: g4_source_escape_campaign and "
                   "g4_photonuclear_xs (optional GEANT4 bake parts) not built")
endif()

# Environment-scatter (air + ground + room) bake. The Python campaign driver
# is env_scatter_fit/ in this directory; it expects these binaries in
# <build>/tools/full_spectrum/.
#
# Offline truth generator for the next-event fluence estimator and producer of
# the CEEFIT-BLOCK-v1 fit-worker blocks; NOT part of the per-event compute()
# path.
find_package(Threads REQUIRED)
add_executable(env_kernel_bake env_kernel_bake.cpp)
target_link_libraries(env_kernel_bake PRIVATE CeeLo Threads::Threads)

# Reassembles a distributed bake from its per-machine shards, refusing
# grid/material mismatches, duplicate nodes, and silently incomplete tables.
add_executable(env_kernel_merge env_kernel_merge.cpp)
target_link_libraries(env_kernel_merge PRIVATE CeeLo)

# Verifies checksums and merge/continuation invariants for fit-worker blocks.
add_executable(env_fit_block_tool env_fit_block_tool.cpp)
target_link_libraries(env_fit_block_tool PRIVATE CeeLo)

# Dependency-free serialized-candidate evaluator used by the Python/C++
# publication-parity gate.  Its strict token protocol avoids a JSON dependency.
add_executable(env_fit_payload_eval env_fit_payload_eval.cpp)
target_link_libraries(env_fit_payload_eval PRIVATE CeeLo)

# Held-out registered-E0 gate: the mandatory registered-energy validation
# producer run by `env_scatter_fit registered-e0-report`.
add_executable(env_registered_e0_gate env_registered_e0_gate.cpp)
target_link_libraries(env_registered_e0_gate PRIVATE CeeLo)

# GEANT4 input generator (GDML + GPS macro) for environment-scatter validation.
add_executable(env_scatter_export env_scatter_export.cpp)
target_link_libraries(env_scatter_export PRIVATE CeeLo)
