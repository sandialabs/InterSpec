# We use Boost.Test for unit testing.  It is the ONLY thing in this project that
# needs Boost, so a missing Boost disables the tests instead of failing the whole
# configure -- `cmake -B build -S .` must work with just a compiler and Eigen3.
# CONFIG mode first (CMake >= 4 removed the FindBoost module), then module mode
# for older Boost installs that ship no CMake package config.
find_package(Boost QUIET CONFIG COMPONENTS unit_test_framework)
if(NOT Boost_FOUND)
    find_package(Boost QUIET COMPONENTS unit_test_framework)
endif()

if(NOT Boost_FOUND OR NOT TARGET Boost::unit_test_framework)
    message(WARNING
        "Boost.Test not found -- unit tests will NOT be built. Install Boost "
        "(e.g. `brew install boost`, `apt install libboost-test-dev`) or pass "
        "-DBoost_DIR=<prefix>/lib/cmake/Boost-<version> to enable them. The "
        "CeeLo library itself does not need Boost.")
    return()
endif()

# Helper function to add a test executable
function(add_ceelo_test TEST_NAME TEST_SOURCE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE
        CeeLo
        Boost::unit_test_framework
    )
    # Define BOOST_TEST_DYN_LINK only when a shared unit_test_framework library is
    # available. The wt-4.13.2 prefix ships a static-only build, for which dynamic
    # linking would leave boost::unit_test::unit_test_main undefined at link time.
    get_target_property(_btf_type Boost::unit_test_framework TYPE)
    if(_btf_type STREQUAL "SHARED_LIBRARY")
        target_compile_definitions(${TEST_NAME} PRIVATE BOOST_TEST_DYN_LINK)
    endif()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction()

# Unit tests
add_ceelo_test(test_geometry test_geometry.cpp)
add_ceelo_test(test_geant4_export test_geant4_export.cpp)
add_ceelo_test(test_materials test_materials.cpp)
add_ceelo_test(test_cross_sections test_cross_sections.cpp)
add_ceelo_test(test_compton test_compton.cpp)
add_ceelo_test(test_efficiency test_efficiency.cpp)
add_ceelo_test(test_transport_physics test_transport_physics.cpp)
add_ceelo_test(test_extended_source test_extended_source.cpp)
add_ceelo_test(test_spherical_source test_spherical_source.cpp)
add_ceelo_test(test_cascade_correction test_cascade_correction.cpp)
add_ceelo_test(test_cascade_angular test_cascade_angular.cpp)
add_ceelo_test(test_analytic_cascade test_analytic_cascade.cpp)
add_ceelo_test(test_analytic_cascade_regression test_analytic_cascade_regression.cpp)
add_ceelo_test(test_csda test_csda.cpp)
add_ceelo_test(test_transmission test_transmission.cpp)
add_ceelo_test(test_solid_angle test_solid_angle.cpp)
# Detector-response payload tests (gated with the payload itself).
if(CEELO_BUILD_RESPONSE)
    add_ceelo_test(test_response_kernel test_response_kernel.cpp)
    add_ceelo_test(test_detector_response_io test_detector_response_io.cpp)
    add_ceelo_test(test_efficiency_transfer test_efficiency_transfer.cpp)
    # D-b structured probe family planner: pure (no MC), always registered.
    add_ceelo_test(test_structured_probe_bank test_structured_probe_bank.cpp)

    # End-to-end generation smoke test runs real (coarse) MC -- opt-in only.
    option(CEELO_RUN_MC_TESTS "Register MC-driven tests (slow) with ctest" OFF)
    if(CEELO_RUN_MC_TESTS)
        # The MC-vs-transfer consistency case is guarded by this define.
        target_compile_definitions(test_efficiency_transfer PRIVATE CEELO_RUN_MC_TESTS)
        add_ceelo_test(test_response_generator_smoke test_response_generator_smoke.cpp)
    else()
        # Still build it (compile coverage); just not registered with ctest.
        add_executable(test_response_generator_smoke test_response_generator_smoke.cpp)
        target_link_libraries(test_response_generator_smoke PRIVATE
            CeeLo Boost::unit_test_framework)
        get_target_property(_btf_type2 Boost::unit_test_framework TYPE)
        if(_btf_type2 STREQUAL "SHARED_LIBRARY")
            target_compile_definitions(test_response_generator_smoke PRIVATE BOOST_TEST_DYN_LINK)
        endif()
    endif()
endif()
add_ceelo_test(test_termination test_termination.cpp)
add_ceelo_test(test_fep_only test_fep_only.cpp)
add_ceelo_test(test_source_shielding test_source_shielding.cpp)
add_ceelo_test(test_asymmetric_shield test_asymmetric_shield.cpp)
add_ceelo_test(test_exponential_source test_exponential_source.cpp)
add_ceelo_test(test_dimension_estimator test_dimension_estimator.cpp)
add_ceelo_test(test_marinelli_beaker test_marinelli_beaker.cpp)
add_ceelo_test(test_position_bias test_position_bias.cpp)
add_ceelo_test(test_biasing_consistency test_biasing_consistency.cpp)
add_ceelo_test(test_virtual_depth test_virtual_depth.cpp)
add_ceelo_test(test_beta_spectrum test_beta_spectrum.cpp)
# Hollow-box (crate shell) and pipe source geometries: analytic path +
# sampling cases are MC-free; the expensive MC-consistency cases are guarded
# by CEELO_RUN_MC_TESTS.
add_ceelo_test(test_hollow_box_source test_hollow_box_source.cpp)
add_ceelo_test(test_pipe_source test_pipe_source.cpp)
if(CEELO_RUN_MC_TESTS)
    target_compile_definitions(test_hollow_box_source PRIVATE CEELO_RUN_MC_TESTS)
    target_compile_definitions(test_pipe_source PRIVATE CEELO_RUN_MC_TESTS)
endif()

# Optional SandiaDecay adapter test (gated on the adapter being built).
if(CEELO_WITH_SANDIADECAY)
    add_executable(test_cascade_adapter test_cascade_adapter.cpp)
    target_link_libraries(test_cascade_adapter PRIVATE
        CeeLoCascade
        Boost::unit_test_framework
    )
    get_target_property(_btf_type Boost::unit_test_framework TYPE)
    if(_btf_type STREQUAL "SHARED_LIBRARY")
        target_compile_definitions(test_cascade_adapter PRIVATE BOOST_TEST_DYN_LINK)
    endif()
    target_compile_definitions(test_cascade_adapter PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml")
    add_test(NAME test_cascade_adapter COMMAND test_cascade_adapter)

    # Cascade summing regression gate vs baked GEANT4 references (low-stats MC).
    add_executable(test_cascade_summing test_cascade_summing.cpp)
    target_link_libraries(test_cascade_summing PRIVATE
        CeeLoCascade
        Boost::unit_test_framework
    )
    if(_btf_type STREQUAL "SHARED_LIBRARY")
        target_compile_definitions(test_cascade_summing PRIVATE BOOST_TEST_DYN_LINK)
    endif()
    target_compile_definitions(test_cascade_summing PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml"
        CASCADE_SUMMING_REF="${CMAKE_CURRENT_SOURCE_DIR}/data/geant4_reference/cascade_summing_multi.csv")
    add_test(NAME test_cascade_summing COMMAND test_cascade_summing)

    # Analytic cascade-summing vs FullRealization MC on real nuclides (exercises
    # the vacancy K/L x-ray + EC-feed + 511 summing channels).
    add_executable(test_analytic_cascade_nuclides test_analytic_cascade_nuclides.cpp)
    target_link_libraries(test_analytic_cascade_nuclides PRIVATE
        CeeLoCascade
        Boost::unit_test_framework
    )
    if(_btf_type STREQUAL "SHARED_LIBRARY")
        target_compile_definitions(test_analytic_cascade_nuclides PRIVATE BOOST_TEST_DYN_LINK)
    endif()
    target_compile_definitions(test_analytic_cascade_nuclides PRIVATE
        SANDIA_DECAY_XML_PATH="${CEELO_SANDIADECAY_DIR}/sandia.decay.xml"
        ANALYTIC_MC_REF="${CMAKE_CURRENT_SOURCE_DIR}/data/ceelo_reference/analytic_cascade_mc.csv")
    add_test(NAME test_analytic_cascade_nuclides COMMAND test_analytic_cascade_nuclides)
endif()

# Benchmark (not a test, just a utility)
add_executable(bench_precision bench_precision.cpp)
target_link_libraries(bench_precision PRIVATE CeeLo)

# FEP-only stage-3 unbiasedness/FOM bench (not a ctest; run manually --
# real MC, minutes). Gate for enabling fep_only on the near-field scan.
if(CEELO_BUILD_RESPONSE)
    add_executable(bench_fep_only_stage3 bench_fep_only_stage3.cpp)
    target_link_libraries(bench_fep_only_stage3 PRIVATE CeeLo)
    target_include_directories(bench_fep_only_stage3 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../tools)
endif()

