cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

cmake_policy(SET CMP0167 NEW)

project( InterSpecTests )

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

if(MSVC)
  option(${PROJECT_NAME}_USE_MSVC_MultiThreadDLL "Use dynamically-link runtime library." OFF)

  if( ${PROJECT_NAME}_USE_MSVC_MultiThreadDLL)
    set(Boost_USE_STATIC_RUNTIME OFF)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  else()
    set(Boost_USE_STATIC_RUNTIME ON)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
  message(STATUS "Using CMAKE_MSVC_RUNTIME_LIBRARY ${CMAKE_MSVC_RUNTIME_LIBRARY}")

  # <windows.h> defines function-like min(a,b)/max(a,b) macros unless NOMINMAX is defined first,
  #  which turns any later `std::max( a, b )`, `std::numeric_limits<T>::max()` or even `obj.max()`
  #  into MSVC error C2589.  This seems to affect the unit tests mostly, so we will take the 
  #  heavy-handed approach of adding it to the compiler command line.
  add_compile_definitions( NOMINMAX )
endif()

enable_testing()


#if (WIN32)
#  add_definitions( -DBOOST_ALL_NO_LIB )
#  add_definitions( -DBOOST_ALL_DYN_LINK )
#endif()
#if(BUILD_AS_UNIT_TEST_SUITE)
#    target_link_libraries(
#        InterSpecLib
#        PUBLIC
#        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
#    )
#endif(BUILD_AS_UNIT_TEST_SUITE)


set( BUILD_AS_UNIT_TEST_SUITE ON CACHE BOOL "" )
set( INCLUDE_ANALYSIS_TEST_SUITE ON CACHE BOOL "" )
set( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE ON CACHE BOOL "")


set(BUILD_AS_LOCAL_SERVER OFF CACHE BOOL "")
set(USE_DB_TO_STORE_SPECTRA ON CACHE BOOL "")
set(PERFORM_DEVELOPER_CHECKS ON CACHE BOOL "")
set(SpecUtils_ENABLE_EQUALITY_CHECKS ON CACHE BOOL "")
set(USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE BOOL "")
set(USE_BATCH_CLI_TOOLS OFF CACHE INTERNAL "")
set(USE_BATCH_GUI_TOOLS ON CACHE INTERNAL "")

set(TEST_SUITE_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test_data"
    CACHE STRING
    "Path to directory that contains the \"analysis_tests\" directory for saving N42 test states."
)

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

# Set the common arguments that are interpreted by boost unit_test
#  Note: not using quotes in next line; if we do, then all arguments will be passed as a quoted
#        single argument, which isnt what we want.
set( BOOST_TEST_CL_ARGS --log_level=test_suite --catch_system_error=yes --report_level=detailed )

# The arguments passed to our code, so the code can know where static nuclear data is, as well spectra and such
set( DATA_DIR_ARGS --datadir=${CMAKE_CURRENT_SOURCE_DIR}/../../data --testfiledir=${CMAKE_CURRENT_SOURCE_DIR}/test_data )


# Drift gate for the Wt resource prune (target/tools/wt_resource_audit/wt_resource_audit.py): fails
# if the committed manifest no longer matches the Wt resources + InterSpec source
# (e.g. after a Wt upgrade adds/renames a resource or a widget's usage changes).
# Pure-Python, so it runs without building InterSpecLib (ctest -R WtResourcePruneManifest).
find_package( Python3 COMPONENTS Interpreter QUIET )
if( Python3_Interpreter_FOUND AND WT_RESOURCES_ORIG_DIR )
  add_test( NAME WtResourcePruneManifest
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/../tools/wt_resource_audit/wt_resource_audit.py" --check
            --wt-resources "${WT_RESOURCES_ORIG_DIR}"
            --interspec-root "${CMAKE_CURRENT_SOURCE_DIR}/../.."
            --manifest "${CMAKE_CURRENT_SOURCE_DIR}/../tools/wt_resource_audit/wt_resources_manifest.json"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif()

add_executable( test_PhysicalUnits test_PhysicalUnits.cpp )
target_link_libraries( test_PhysicalUnits PRIVATE InterSpecLib )
add_test( NAME PhysicalUnits
  COMMAND $<TARGET_FILE:test_PhysicalUnits> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PhysicalUnitsLocalized test_PhysicalUnitsLocalized.cpp )
target_link_libraries( test_PhysicalUnitsLocalized PRIVATE InterSpecLib )
add_test( NAME test_PhysicalUnitsLocalized
  COMMAND $<TARGET_FILE:test_PhysicalUnitsLocalized> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_findCharacteristics test_findCharacteristics.cpp )
target_link_libraries( test_findCharacteristics PRIVATE InterSpecLib )
add_test( NAME FindCharacteristicNuclides
  COMMAND $<TARGET_FILE:test_findCharacteristics> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_nuclideSearch test_nuclideSearch.cpp )
target_link_libraries( test_nuclideSearch PRIVATE InterSpecLib )
add_test( NAME NuclideSearch
  COMMAND $<TARGET_FILE:test_nuclideSearch> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_RefLines test_RefLines.cpp )
target_link_libraries( test_RefLines PRIVATE InterSpecLib )
add_test( NAME TesRefLines
  COMMAND $<TARGET_FILE:test_RefLines> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_CascadeRefLineSumming test_CascadeRefLineSumming.cpp )
target_link_libraries( test_CascadeRefLineSumming PRIVATE InterSpecLib )
add_test( NAME TCascadeRefLineSumming
  COMMAND $<TARGET_FILE:test_CascadeRefLineSumming> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_DoseCalc test_DoseCalc.cpp )
target_link_libraries( test_DoseCalc PRIVATE InterSpecLib )
add_test( NAME TDoseCalc
  COMMAND $<TARGET_FILE:test_DoseCalc> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_RelEffManual test_RelEffManual.cpp )
target_link_libraries( test_RelEffManual PRIVATE InterSpecLib )
add_test( NAME TestRelEffManual
  COMMAND $<TARGET_FILE:test_RelEffManual> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelEffAuto test_RelEffAuto.cpp )
target_link_libraries( test_RelEffAuto PRIVATE InterSpecLib )
add_test( NAME TestRelEffAuto
  COMMAND $<TARGET_FILE:test_RelEffAuto> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActAutoTruthManifest test_RelActAutoTruthManifest.cpp )
target_link_libraries( test_RelActAutoTruthManifest PRIVATE InterSpecLib )
add_test( NAME TRelActAutoTruthManifest
  COMMAND $<TARGET_FILE:test_RelActAutoTruthManifest> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActAutoReport test_RelActAutoReport.cpp )
target_link_libraries( test_RelActAutoReport PRIVATE InterSpecLib )
add_test( NAME TRelActAutoReport
  COMMAND $<TARGET_FILE:test_RelActAutoReport> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_ProfileApi test_RelActCalcAuto_ProfileApi.cpp )
target_link_libraries( test_RelActCalcAuto_ProfileApi PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_ProfileApi
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_ProfileApi> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_ProfileFit test_RelActCalcAuto_ProfileFit.cpp )
target_link_libraries( test_RelActCalcAuto_ProfileFit PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_ProfileFit
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_ProfileFit> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_ProfileScan test_RelActCalcAuto_ProfileScan.cpp )
target_link_libraries( test_RelActCalcAuto_ProfileScan PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_ProfileScan
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_ProfileScan> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Developer tool, deliberately not a ctest: it sweeps the local IDB/JRC corpora (about 1 GB of
# licensed scratch data that cannot live in the repository) for non-graceful RelActCalcAuto behavior
# - NaN/Inf, out-of-range mass fractions, unusable statuses with no explanation, throws, hangs - and
# dumps per-profile-sample delta-chi2 so scan convergence can be studied over a population.
add_executable( rel_act_corpus_sweep rel_act_corpus_sweep.cpp )
target_link_libraries( rel_act_corpus_sweep PRIVATE InterSpecLib )

# The JRC known-truth plutonium grid always builds, so its sources are compile-checked with every
# other test, but it is not registered as a ctest by default: each of the three cases runs a free-age
# solve plus bounded profile-likelihood scans over every Pu isotope, which makes the whole target a
# multi-hour run.  Enable the option to register it, and run it with `-L long`.
# Follows the INTERSPEC_ENABLE_DECON_COVERAGE_TEST precedent further down this file.
add_executable( test_RelActCalcAuto_Pu610775 test_RelActCalcAuto_Pu610775.cpp )
target_link_libraries( test_RelActCalcAuto_Pu610775 PRIVATE InterSpecLib )
option( INTERSPEC_ENABLE_JRC_PU_GRID_TEST
        "Register the licensed-fixture JRC plutonium profile grid as a ctest (runs for hours)" OFF )
if( INTERSPEC_ENABLE_JRC_PU_GRID_TEST )
  add_test( NAME TRelActCalcAuto_Pu610775
    COMMAND $<TARGET_FILE:test_RelActCalcAuto_Pu610775> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
  set_tests_properties( TRelActCalcAuto_Pu610775 PROPERTIES
    TIMEOUT 7200
    LABELS "long" )
endif()

add_executable( test_BatchSampleSelect test_BatchSampleSelect.cpp )
target_link_libraries( test_BatchSampleSelect PRIVATE InterSpecLib )
add_test( NAME TBatchSampleSelect
  COMMAND $<TARGET_FILE:test_BatchSampleSelect> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_BatchPeakMda test_BatchPeakMda.cpp )
target_link_libraries( test_BatchPeakMda PRIVATE InterSpecLib )
add_test( NAME TBatchPeakMda
  COMMAND $<TARGET_FILE:test_BatchPeakMda> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_BatchRelActAuto test_BatchRelActAuto.cpp )
target_link_libraries( test_BatchRelActAuto PRIVATE InterSpecLib )
add_test( NAME TBatchRelActAuto
  COMMAND $<TARGET_FILE:test_BatchRelActAuto> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PeakDists test_PeakDists.cpp )
target_link_libraries( test_PeakDists PRIVATE InterSpecLib )
add_test( NAME TPeakDists
  COMMAND $<TARGET_FILE:test_PeakDists> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_CrystalBall_CZT_fit test_CrystalBall_CZT_fit.cpp )
target_link_libraries( test_CrystalBall_CZT_fit PRIVATE InterSpecLib )
add_test( NAME TCrystalBall_CZT_fit
  COMMAND $<TARGET_FILE:test_CrystalBall_CZT_fit> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PeakFitLM_MultiPeak test_PeakFitLM_MultiPeak.cpp )
target_link_libraries( test_PeakFitLM_MultiPeak PRIVATE InterSpecLib )
add_test( NAME TPeakFitLM_MultiPeak
  COMMAND $<TARGET_FILE:test_PeakFitLM_MultiPeak> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_EnergyCal_imp test_RelActCalcAuto_EnergyCal_imp.cpp )
target_link_libraries( test_RelActCalcAuto_EnergyCal_imp PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_EnergyCal_imp
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_EnergyCal_imp> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_AgeDerivative test_RelActCalcAuto_AgeDerivative.cpp )
target_link_libraries( test_RelActCalcAuto_AgeDerivative PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_AgeDerivative
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_AgeDerivative> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_ConvergenceAudit test_RelActCalcAuto_ConvergenceAudit.cpp )
target_link_libraries( test_RelActCalcAuto_ConvergenceAudit PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_ConvergenceAudit
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_ConvergenceAudit> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_Search test_RelActCalcAuto_Search.cpp )
target_link_libraries( test_RelActCalcAuto_Search PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_Search
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_Search> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_Conditioning test_RelActCalcAuto_Conditioning.cpp )
target_link_libraries( test_RelActCalcAuto_Conditioning PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_Conditioning
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_Conditioning> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_CovarianceQuality test_RelActCalcAuto_CovarianceQuality.cpp )
target_link_libraries( test_RelActCalcAuto_CovarianceQuality PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_CovarianceQuality
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_CovarianceQuality> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_MassFracConstraint test_RelActCalcAuto_MassFracConstraint.cpp )
target_link_libraries( test_RelActCalcAuto_MassFracConstraint PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_MassFracConstraint
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_MassFracConstraint> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_FloatingPeakFwhm test_RelActCalcAuto_FloatingPeakFwhm.cpp )
target_link_libraries( test_RelActCalcAuto_FloatingPeakFwhm PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_FloatingPeakFwhm
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_FloatingPeakFwhm> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_ObsEffFilter test_RelActCalcAuto_ObsEffFilter.cpp )
target_link_libraries( test_RelActCalcAuto_ObsEffFilter PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_ObsEffFilter
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_ObsEffFilter> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_RelActCalcAuto_MultiCurve test_RelActCalcAuto_MultiCurve.cpp )
target_link_libraries( test_RelActCalcAuto_MultiCurve PRIVATE InterSpecLib )
add_test( NAME TRelActCalcAuto_MultiCurve
  COMMAND $<TARGET_FILE:test_RelActCalcAuto_MultiCurve> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_XRayWidthServer test_XRayWidthServer.cpp )
target_link_libraries( test_XRayWidthServer PRIVATE InterSpecLib )
add_test( NAME TXRayWidthServer
  COMMAND $<TARGET_FILE:test_XRayWidthServer> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)


add_executable( test_AppUtils test_AppUtils.cpp )
target_link_libraries( test_AppUtils PRIVATE InterSpecLib )
add_test( NAME TAppUtils
  COMMAND $<TARGET_FILE:test_AppUtils> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ShieldingSourceFitCalc test_ShieldingSourceFitCalc.cpp )
target_link_libraries( test_ShieldingSourceFitCalc PRIVATE InterSpecLib )
add_test( NAME TShieldingSourceFitCalc
  COMMAND $<TARGET_FILE:test_ShieldingSourceFitCalc> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_VolumetricIntegration test_VolumetricIntegration.cpp VolumetricIntegrationFixtures.h )
target_link_libraries( test_VolumetricIntegration PRIVATE InterSpecLib )
add_test( NAME TVolumetricIntegration
  COMMAND $<TARGET_FILE:test_VolumetricIntegration> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ShieldingDimLimit test_ShieldingDimLimit.cpp VolumetricIntegrationFixtures.h )
target_link_libraries( test_ShieldingDimLimit PRIVATE InterSpecLib )
add_test( NAME TShieldingDimLimit
  COMMAND $<TARGET_FILE:test_ShieldingDimLimit> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Act/Shield-fit truth tests vs CeeLo FullRealization cascade MC (Detective-X
# model grounded to the LANL curve).  Fixtures live in test_data/cascade_truth/
# (committed); the test never runs MC unless invoked with `--regenerate`.
add_executable( test_CascadeSummingFit test_CascadeSummingFit.cpp )
target_link_libraries( test_CascadeSummingFit PRIVATE InterSpecLib )
add_test( NAME TCascadeSummingFit
  COMMAND $<TARGET_FILE:test_CascadeSummingFit> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# CeeLo MC-parameterized detector-response validation: golden response XMLs +
# fresh MC probe banks live in test_data/ceelo_drf/ (committed fixtures); the
# test never runs Monte Carlo.
add_executable( test_CeeLoDrfIntegration test_CeeLoDrfIntegration.cpp )
target_link_libraries( test_CeeLoDrfIntegration PRIVATE InterSpecLib )
add_test( NAME TCeeLoDrfIntegration
  COMMAND $<TARGET_FILE:test_CeeLoDrfIntegration> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Regenerating the CeeLo golden fixtures runs REAL Monte Carlo (several
# minutes); opt in with -DINTERSPEC_REGEN_CEELO_GOLDENS=ON, then build the
# `regen_ceelo_goldens` target.  Deterministic seeds: regenerated probe
# truths reproduce within their recorded MC sigma.
option( INTERSPEC_REGEN_CEELO_GOLDENS "Add a target to regenerate CeeLo golden fixtures (runs minutes of MC)" OFF )
if( INTERSPEC_REGEN_CEELO_GOLDENS )
  set( _ceelo_fixture_dir "${CMAKE_CURRENT_SOURCE_DIR}/test_data/ceelo_drf" )
  add_custom_target( regen_ceelo_goldens
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_ceelo_fixture_dir}"
    # High-statistics fixtures: tight per-node precision with raised caps
    #  (max_events / max_cpu_s) so it is actually reached. czt_box is a box
    #  (4x the angular nodes via its phi axis), so a lighter target keeps it
    #  practical. Args: <preset> <dir> <profile> <precision> <n_probes> <seed>
    #  <max_events> <max_cpu_s> (CPU-s summed over worker threads; ~10x the
    #  old wall-s caps at ~10 threads, so the same budget, load-invariant).
    COMMAND $<TARGET_FILE:make_golden_response> nai3x3 "${_ceelo_fixture_dir}" general 0.001 120 1 200000000 600
    COMMAND $<TARGET_FILE:make_golden_response> hpge_coax "${_ceelo_fixture_dir}" general 0.001 120 1 200000000 600
    COMMAND $<TARGET_FILE:make_golden_response> czt_box "${_ceelo_fixture_dir}" far-field 0.0015 90 1 200000000 250
    # Detective-X-like coax HPGe (public specs), the Act/Shield cascade-truth
    #  detector; grounded to the ORTEC Detective-X_LANL_100cm curve by
    #  test_CascadeSummingFit's --regenerate mode (see test_data/cascade_truth/).
    COMMAND $<TARGET_FILE:make_golden_response> detective_x "${_ceelo_fixture_dir}" general 0.001 120 1 200000000 600
    COMMENT "Regenerating CeeLo golden response fixtures (real MC; ~hours at high stats)"
    VERBATIM
  )
  add_dependencies( regen_ceelo_goldens make_golden_response )
endif( INTERSPEC_REGEN_CEELO_GOLDENS )

option( RECORD_GEOM_GOLDENS "Print (instead of check) geometry golden values" OFF )
add_executable( test_ShieldingGeomGolden test_ShieldingGeomGolden.cpp VolumetricIntegrationFixtures.h )
target_link_libraries( test_ShieldingGeomGolden PRIVATE InterSpecLib )
if( RECORD_GEOM_GOLDENS )
  target_compile_definitions( test_ShieldingGeomGolden PRIVATE RECORD_GEOM_GOLDENS=1 )
endif()
add_test( NAME TShieldingGeomGolden
  COMMAND $<TARGET_FILE:test_ShieldingGeomGolden> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_DetectionLimit test_DetectionLimit.cpp )
target_link_libraries( test_DetectionLimit PRIVATE InterSpecLib )
add_test( NAME test_DetectionLimit
  COMMAND $<TARGET_FILE:test_DetectionLimit> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_DetectionLimitState test_DetectionLimitState.cpp )
# Needs Wt::Test as well: the state round-trip cases drive the real tool widgets, which is the
#  only level at which a control missing from `encodeStateToUrl()` is visible.
target_link_libraries( test_DetectionLimitState PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
add_test( NAME TDetectionLimitState
  COMMAND $<TARGET_FILE:test_DetectionLimitState> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# The deconvolution coverage study is a report, not a gate: it is opt-in, it is held to a
#  wall-clock budget rather than a fixed trial count, and `DeconCoverageGate` inside
#  test_DetectionLimit is what actually guards against regressions on every build.
option( INTERSPEC_ENABLE_DECON_COVERAGE_TEST
        "Register the opt-in deconvolution coverage study as a ctest (runs about an hour)" OFF )
if( INTERSPEC_ENABLE_DECON_COVERAGE_TEST )
  add_test( NAME TDeconCoverageStudy
    COMMAND $<TARGET_FILE:test_DetectionLimit> --run_test=DeconCoverageStudy
            ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
  set_tests_properties( TDeconCoverageStudy PROPERTIES
    TIMEOUT 7200
    LABELS "long"
    ENVIRONMENT "INTERSPEC_RUN_DECON_COVERAGE=1" )
endif()

add_executable( test_AppTextBundles test_AppTextBundles.cpp )
target_link_libraries( test_AppTextBundles PRIVATE InterSpecLib )
add_test( NAME TAppTextBundles
  COMMAND $<TARGET_FILE:test_AppTextBundles> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ShieldSourcePullTrend test_ShieldSourcePullTrend.cpp )
target_link_libraries( test_ShieldSourcePullTrend PRIVATE InterSpecLib )
add_test( NAME TShieldSourcePullTrend
  COMMAND $<TARGET_FILE:test_ShieldSourcePullTrend> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_MakeDrfSrcDef test_MakeDrfSrcDef.cpp )
target_link_libraries( test_MakeDrfSrcDef PRIVATE InterSpecLib )
add_test( NAME TMakeDrfSrcDef
  COMMAND $<TARGET_FILE:test_MakeDrfSrcDef> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_GammaInteractionCalc test_GammaInteractionCalc.cpp GeometryReference.h )
target_link_libraries( test_GammaInteractionCalc PRIVATE InterSpecLib )
add_test( NAME TGammaInteractionCalc
  COMMAND $<TARGET_FILE:test_GammaInteractionCalc> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Proves the double ray-tracing primitives, now thin `_imp<double>` wrappers, still agree with the
#  frozen pre-refactor bodies (LegacyGeometryRef.h) everywhere except the three documented bug
#  fixes, and that those fixes agree with independent references (GeometryReference.h).
#  Pure arithmetic - needs no static data - so it is fast enough to use as an inner loop.
#  Retires together with LegacyGeometryRef.h once the migration has shipped a release cycle.
add_executable( test_GeometryLegacyParity test_GeometryLegacyParity.cpp
                LegacyGeometryRef.h GeometryReference.h )
target_link_libraries( test_GeometryLegacyParity PRIVATE InterSpecLib )
add_test( NAME TGeometryLegacyParity
  COMMAND $<TARGET_FILE:test_GeometryLegacyParity> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# VoigtDistribution unit tests (requires Eigen3 and Ceres)
find_package(Eigen3 3.4...6.0 QUIET)
find_package(Ceres QUIET)
find_package(Boost COMPONENTS unit_test_framework QUIET)

if(Eigen3_FOUND AND Ceres_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "Building VoigtDistribution unit tests (Eigen3, Ceres, and Boost unit_test_framework found)")

    add_executable( test_VoigtDistribution
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/test_pseudo_voigt_exp_tail.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/test_faddeeva_templated.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/test_voigt_against_mit_imp.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/voigt_cdf_reference.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/voigt_distribution_tests.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/3rd_party/Faddeeva_MIT/Faddeeva.cc
    )

    target_include_directories( test_VoigtDistribution PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution
        ${CMAKE_CURRENT_SOURCE_DIR}/../../external_libs/VoigtDistribution/unit_tests/3rd_party
        ${EIGEN3_INCLUDE_DIRS}
    )

    target_link_libraries( test_VoigtDistribution PRIVATE
        Eigen3::Eigen
        Ceres::ceres
        Boost::unit_test_framework
    )

    add_test( NAME TVoigtDistribution
        COMMAND $<TARGET_FILE:test_VoigtDistribution> ${BOOST_TEST_CL_ARGS}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
else()
    message(STATUS "Skipping VoigtDistribution unit tests (requires Eigen3, Ceres, and Boost unit_test_framework)")
endif()

add_executable( test_SimpleActivityCalc test_SimpleActivityCalc.cpp )
target_link_libraries( test_SimpleActivityCalc PRIVATE InterSpecLib )
add_test( NAME TSimpleActivityCalc
  COMMAND $<TARGET_FILE:test_SimpleActivityCalc> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_BersteinPolynomial test_BersteinPolynomial.cpp )
target_link_libraries( test_BersteinPolynomial PRIVATE InterSpecLib )
add_test( NAME TBersteinPolynomial
  COMMAND $<TARGET_FILE:test_BersteinPolynomial> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_MakeDrfFit test_MakeDrfFit.cpp )
target_link_libraries( test_MakeDrfFit PRIVATE InterSpecLib )
add_test( NAME TMakeDrfFit
  COMMAND $<TARGET_FILE:test_MakeDrfFit> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_EnergyCal test_EnergyCal.cpp )
target_link_libraries( test_EnergyCal PRIVATE InterSpecLib )
add_test( NAME TEnergyCal
  COMMAND $<TARGET_FILE:test_EnergyCal> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_EnergyCalGainMatch test_EnergyCalGainMatch.cpp )
target_link_libraries( test_EnergyCalGainMatch PRIVATE InterSpecLib )
add_test( NAME TEnergyCalGainMatch
  COMMAND $<TARGET_FILE:test_EnergyCalGainMatch> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_EnergyCalGainGuess test_EnergyCalGainGuess.cpp )
target_link_libraries( test_EnergyCalGainGuess PRIVATE InterSpecLib )
add_test( NAME TEnergyCalGainGuess
  COMMAND $<TARGET_FILE:test_EnergyCalGainGuess> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PeakModel test_PeakModel.cpp )
target_link_libraries( test_PeakModel PRIVATE InterSpecLib )
add_test( NAME TPeakModel
  COMMAND $<TARGET_FILE:test_PeakModel> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_fitPeaksForSources test_fitPeaksForSources.cpp )
target_link_libraries( test_fitPeaksForSources PRIVATE InterSpecLib )
add_test( NAME TFitPeaksForSources
  COMMAND $<TARGET_FILE:test_fitPeaksForSources> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_BackgroundPeakRecovery test_BackgroundPeakRecovery.cpp )
target_link_libraries( test_BackgroundPeakRecovery PRIVATE InterSpecLib )
add_test( NAME TBackgroundPeakRecovery
  COMMAND $<TARGET_FILE:test_BackgroundPeakRecovery> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

if( USE_LLM_INTERFACE )
  add_executable( test_LlmToolRegistry test_LlmToolRegistry.cpp )
  target_link_libraries( test_LlmToolRegistry PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmToolRegistry
    COMMAND $<TARGET_FILE:test_LlmToolRegistry> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmPromptTemplate test_LlmPromptTemplate.cpp )
  target_link_libraries( test_LlmPromptTemplate PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmPromptTemplate
    COMMAND $<TARGET_FILE:test_LlmPromptTemplate> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmIsotopicsTool test_LlmIsotopicsTool.cpp )
  target_link_libraries( test_LlmIsotopicsTool PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmIsotopicsTool
    COMMAND $<TARGET_FILE:test_LlmIsotopicsTool> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmActivityFitTool test_LlmActivityFitTool.cpp )
  target_link_libraries( test_LlmActivityFitTool PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmActivityFitTool
    COMMAND $<TARGET_FILE:test_LlmActivityFitTool> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmRelActManualTool test_LlmRelActManualTool.cpp )
  target_link_libraries( test_LlmRelActManualTool PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmRelActManualTool
    COMMAND $<TARGET_FILE:test_LlmRelActManualTool> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmInterface test_LlmInterface.cpp )
  target_link_libraries( test_LlmInterface PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmInterface
    COMMAND $<TARGET_FILE:test_LlmInterface> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_AnalystChecks test_AnalystChecks.cpp )
  target_link_libraries( test_AnalystChecks PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TAnalystChecks
    COMMAND $<TARGET_FILE:test_AnalystChecks> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_executable( test_LlmTimeHistoryTools test_LlmTimeHistoryTools.cpp )
  target_link_libraries( test_LlmTimeHistoryTools PRIVATE InterSpecLib ${Wt_TEST_LIBRARY} )
  add_test( NAME TLlmTimeHistoryTools
    COMMAND $<TARGET_FILE:test_LlmTimeHistoryTools> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endif( USE_LLM_INTERFACE )

add_executable( test_DetectorPeakResponse test_DetectorPeakResponse.cpp )
target_link_libraries( test_DetectorPeakResponse PRIVATE InterSpecLib )
add_test( NAME TDetectorPeakResponse
  COMMAND $<TARGET_FILE:test_DetectorPeakResponse> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_DetectorEfficiency test_DetectorEfficiency.cpp )
target_link_libraries( test_DetectorEfficiency PRIVATE InterSpecLib )
add_test( NAME TDetectorEfficiency
  COMMAND $<TARGET_FILE:test_DetectorEfficiency> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ZipArchive test_ZipArchive.cpp )
target_link_libraries( test_ZipArchive PRIVATE InterSpecLib )
add_test( NAME TZipArchive
  COMMAND $<TARGET_FILE:test_ZipArchive> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_SecurityHardening test_SecurityHardening.cpp )
target_link_libraries( test_SecurityHardening PRIVATE InterSpecLib )
add_test( NAME TSecurityHardening
  COMMAND $<TARGET_FILE:test_SecurityHardening> ${BOOST_TEST_CL_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PeakFitLM test_PeakFitLM.cpp )
target_link_libraries( test_PeakFitLM PRIVATE InterSpecLib )
add_test( NAME TPeakFitLM
  COMMAND $<TARGET_FILE:test_PeakFitLM> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_PeakFitDetPrefs test_PeakFitDetPrefs.cpp )
target_link_libraries( test_PeakFitDetPrefs PRIVATE InterSpecLib )
add_test( NAME TPeakFitDetPrefs
  COMMAND $<TARGET_FILE:test_PeakFitDetPrefs> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_GadrasShieldScatter test_GadrasShieldScatter.cpp )
target_link_libraries( test_GadrasShieldScatter PRIVATE InterSpecLib )
add_test( NAME TGadrasShieldScatter
  COMMAND $<TARGET_FILE:test_GadrasShieldScatter> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ExportSpecFile test_ExportSpecFile.cpp )
target_link_libraries( test_ExportSpecFile PRIVATE InterSpecLib )
add_test( NAME TExportSpecFile
  COMMAND $<TARGET_FILE:test_ExportSpecFile> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable( test_ExportCAM test_ExportCAM.cpp )
target_link_libraries( test_ExportCAM PRIVATE InterSpecLib )
add_test( NAME TExportCAM
  COMMAND $<TARGET_FILE:test_ExportCAM> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# ./test_offlineAnalysis.cpp
#target_link_libraries(test_offlineAnalysis PUBLIC 
#            debug "${Wt_TEST_DEBUG_LIBRARY}"
#            optimized "${Wt_TEST_LIBRARY}"
#)


# For some reason in Xcode, the RUN_TESTS target doesnt cause the executables to be built if there
# is code changes - it will just run the already compiled tests.  To get around this we will define
# a target that depends on ALL_BUILD and RUN_TESTS.
# Note however, that Xcode wont show result of tests - it only will display an error if a test fails
if( CMAKE_GENERATOR STREQUAL "Xcode" )
  add_custom_target( build_then_run_tests DEPENDS ALL_BUILD RUN_TESTS )
endif( CMAKE_GENERATOR STREQUAL "Xcode" )

# libFuzzer-based fuzz harnesses live in a separate, standalone CMake project at
# target/fuzzing/ - they cannot be mixed with the unit-test build because the
# sanitizer-instrumented InterSpecLib does not link cleanly with un-instrumented
# test binaries.  See target/fuzzing/README.md.
