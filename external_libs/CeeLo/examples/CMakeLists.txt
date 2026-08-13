function(add_ceelo_example EXE_NAME SRC_FILE)
    add_executable(${EXE_NAME} ${SRC_FILE})
    target_link_libraries(${EXE_NAME} PRIVATE CeeLo)
endfunction()

add_ceelo_example(simple_point_source    simple_point_source.cpp)
add_ceelo_example(extended_source        extended_source.cpp)
add_ceelo_example(cascade_correction     cascade_correction.cpp)
add_ceelo_example(generate_g4_benchmarks generate_g4_benchmarks.cpp)
add_ceelo_example(benchmark_mc_configs   benchmark_mc_configs.cpp)
add_ceelo_example(vpd_generate            vpd_generate.cpp)
add_ceelo_example(export_cascade_g4        export_cascade_g4.cpp)
add_ceelo_example(sphere_g4_validation     sphere_g4_validation.cpp)

# NOTE: one-off research studies (study0/7/8/9, study_grids, efficiency_grid_study,
# cascade_shielding_study) and perf micro-benchmarks (bench_compton_cap,
# bench_perf_verify) live in the dev-only studies/ tree, which is not part of the
# distributed source (studies/CMakeLists.txt),
# not here.

# Cascade-summing validation study (needs the SandiaDecay adapter).
if(CEELO_WITH_SANDIADECAY)
    add_executable(cascade_validation cascade_validation.cpp)
    target_link_libraries(cascade_validation PRIVATE CeeLoCascade)
    target_compile_definitions(cascade_validation PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml")

    add_executable(cascade_wtheta_demo cascade_wtheta_demo.cpp)
    target_link_libraries(cascade_wtheta_demo PRIVATE CeeLoCascade)
    target_compile_definitions(cascade_wtheta_demo PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml")

    add_executable(cascade_g4_compare cascade_g4_compare.cpp)
    target_link_libraries(cascade_g4_compare PRIVATE CeeLoCascade)
    target_compile_definitions(cascade_g4_compare PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml")

    add_executable(cascade_g4_spectrum cascade_g4_spectrum.cpp)
    target_link_libraries(cascade_g4_spectrum PRIVATE CeeLoCascade)
    target_compile_definitions(cascade_g4_spectrum PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml")

    # Cascade-summing MC observable producer for profiling/compare_validation.py.
    # Shares tests/cascade_ref_common.h (the alcyl run + reference reader) with the
    # ctest gate, so the geometry and observable definitions cannot drift.
    add_executable(cascade_observables cascade_observables.cpp)
    target_link_libraries(cascade_observables PRIVATE CeeLoCascade)
    target_include_directories(cascade_observables PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_compile_definitions(cascade_observables PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml"
        CASCADE_SUMMING_REF="${CMAKE_SOURCE_DIR}/tests/data/geant4_reference/cascade_summing_multi.csv")
endif()

# Comprehensive spectrum generator for all validation configs.  Lives under
# tools/geant4_validation/ (it belongs with the GEANT4 harness) but is built here
# because it links the CeeLo library, not GEANT4.
add_executable(generate_all_spectra
    ${CMAKE_SOURCE_DIR}/tools/geant4_validation/generate_all_spectra.cpp)
target_link_libraries(generate_all_spectra PRIVATE CeeLo)
