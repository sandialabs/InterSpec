# Tools for improving FitPeaksForNuclides::fit_peaks_for_nuclides against peak-fit corpora (see README.md).
# Added from target/testing/CMakeLists.txt, which provides the InterSpecLib target and the test arguments;
# the executables are written to the top of that build tree (e.g. target/testing/build_ninja/).

if( NOT TARGET InterSpecLib )
  message( FATAL_ERROR "peak_fit_improve_ai is built as part of target/testing (add_subdirectory from there)" )
endif()

# Fit-peaks-for-nuclides corpus evaluation (developer tool, deliberately not a ctest): scores
# FitPeaksForNuclides::fit_peaks_for_nuclides against hand-fit reference N42 files or GADRAS-injected
# spectra with truth, kept outside the repository, with parameter sweeps, determinism repeats and an
# HTML review gallery.  See harness/fit_peaks_corpus_eval.cpp (--help) for usage.
add_executable( fit_peaks_corpus_eval harness/fit_peaks_corpus_eval.cpp
                harness/FitPeaksCorpusScore.h harness/FitPeaksCorpusScore.cpp
                harness/FitPeaksCorpusReport.h harness/FitPeaksCorpusReport.cpp
                harness/PeakSearchEval.h harness/PeakSearchEval.cpp )
target_link_libraries( fit_peaks_corpus_eval PRIVATE InterSpecLib )

# Unit tests of the corpus scorer (synthetic peak sets; no data files needed).
add_executable( test_fitPeaksCorpusScore harness/test_fitPeaksCorpusScore.cpp
                harness/FitPeaksCorpusScore.h harness/FitPeaksCorpusScore.cpp
                harness/FitPeaksCorpusReport.h harness/FitPeaksCorpusReport.cpp )
target_link_libraries( test_fitPeaksCorpusScore PRIVATE InterSpecLib )
add_test( NAME TFitPeaksCorpusScore
  COMMAND $<TARGET_FILE:test_fitPeaksCorpusScore> ${BOOST_TEST_CL_ARGS} -- ${DATA_DIR_ARGS}
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# PeakFitLM fit-objective comparison (developer tool, not a ctest): bias, pulls, failure tail and CPU
# of the chi2 vs Poisson-likelihood objectives on Poisson replicas of synthetic and GADRAS-inject
# spectra.  See harness/peak_fit_objective_eval.cpp (--help) for usage.
add_executable( peak_fit_objective_eval harness/peak_fit_objective_eval.cpp
                harness/FitPeaksCorpusScore.h harness/FitPeaksCorpusScore.cpp )
target_link_libraries( peak_fit_objective_eval PRIVATE InterSpecLib )

set_target_properties( fit_peaks_corpus_eval test_fitPeaksCorpusScore peak_fit_objective_eval
                       PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )
