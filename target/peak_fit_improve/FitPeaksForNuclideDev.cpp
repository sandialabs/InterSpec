/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <chrono>
#include <random>
#include <string>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp> // For boost::asio::post

#include <Wt/WColor>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "FinalFit_GA.h"
#include "PeakFitImprove.h"
#include "InitialFit_GA.h"
#include "CandidatePeak_GA.h"

#include "InterSpec/FitPeaksForNuclides.h"

using namespace std;

// Types from FitPeaksForNuclides namespace - using type aliases for convenience
using PeakFitResult = FitPeaksForNuclides::PeakFitResult;
using PeakFitForNuclideConfig = FitPeaksForNuclides::PeakFitForNuclideConfig;
using GammaClusteringSettings = FitPeaksForNuclides::GammaClusteringSettings;

namespace FitPeaksForNuclideDev
{

// Struct to store cluster info including gamma line energies and amplitudes
struct ClusteredGammaInfo {
  double lower;
  double upper;
  vector<double> gamma_energies;    // energies of gamma lines in this cluster
  vector<double> gamma_amplitudes;  // expected peak areas/amplitudes
};


// Helper functions get_source_photons, get_source_age, should_use_desperation_shielding,
// and create_desperation_shielding have been moved to FitPeaksForNuclides namespace


// Forward declaration
// cluster_gammas_to_rois forward declaration removed - function moved to FitPeaksForNuclides namespace


// PeakFitResult is now defined in FitPeaksForNuclides namespace
// Using type alias for convenience
using PeakFitResult = FitPeaksForNuclides::PeakFitResult;

// GammaClusteringSettings is now defined in FitPeaksForNuclides namespace
// Using type alias for convenience
using GammaClusteringSettings = FitPeaksForNuclides::GammaClusteringSettings;


/** Combined scoring struct that evaluates both final fit quality and peak finding completeness.

 This struct provides a comprehensive assessment by combining:
 - FinalFitScore: Area, width, and position accuracy metrics
 - PeakFindAndFitWeights: Peak detection completeness and area accuracy
 - final_weight: Combined score (simple sum of find_weight + total_weight)
 */
struct CombinedPeakFitScore
{
  FinalFit_GA::FinalFitScore final_fit_score;
  InitialFit_GA::PeakFindAndFitWeights initial_fit_weights;
  CandidatePeak_GA::CandidatePeakScore candidate_peak_score;

  /** Combined weight score (find_weight + total_weight).
   * find_weight: higher is better (more correct peaks found)
   * total_weight: lower is better (more accurate areas)
   */
  double final_weight = 0.0;

  std::string print( const std::string &varname ) const
  {
    std::string answer;

    answer += "=== " + varname + " - Final Fit Metrics ===\n";
    answer += final_fit_score.print( varname + ".final_fit_score" );

    answer += "\n=== " + varname + " - Peak Finding Metrics ===\n";
    answer += varname + ".initial_fit_weights.find_weight =                     " + std::to_string(initial_fit_weights.find_weight) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_weight =          " + std::to_string(initial_fit_weights.def_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_weight =        " + std::to_string(initial_fit_weights.maybe_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_weight =          " + std::to_string(initial_fit_weights.not_wanted_area_weight) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_median_weight =   " + std::to_string(initial_fit_weights.def_wanted_area_median_weight) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_median_weight = " + std::to_string(initial_fit_weights.maybe_wanted_area_median_weight) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_median_weight =   " + std::to_string(initial_fit_weights.not_wanted_area_median_weight) + "\n";

    answer += "\n=== " + varname + " - Candidate Peak Metrics ===\n";
    answer += varname + ".candidate_peak_score.score =                           " + std::to_string(candidate_peak_score.score) + "\n";
    answer += varname + ".candidate_peak_score.num_peaks_found =                 " + std::to_string(candidate_peak_score.num_peaks_found) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_not_found =         " + std::to_string(candidate_peak_score.num_def_wanted_not_found) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_peaks_found =      " + std::to_string(candidate_peak_score.num_def_wanted_peaks_found) + "\n";
    answer += varname + ".candidate_peak_score.num_possibly_accepted_peaks_not_found = " + std::to_string(candidate_peak_score.num_possibly_accepted_peaks_not_found) + "\n";
    answer += varname + ".candidate_peak_score.num_extra_peaks =                " + std::to_string(candidate_peak_score.num_extra_peaks) + "\n";

    answer += "\n=== " + varname + " - Combined Metrics ===\n";
    answer += varname + ".final_weight = " + std::to_string(final_weight) + "\n";

    return answer;
  }//print(...)
};//struct CombinedPeakFitScore


// PeakFitForNuclideConfig is now defined in FitPeaksForNuclides namespace
// Using type alias for convenience
using PeakFitForNuclideConfig = FitPeaksForNuclides::PeakFitForNuclideConfig;


// Helper functions should_combine_peaks, combine_peaks, combine_overlapping_peaks_in_rois,
// and compute_observable_peaks have been moved to FitPeaksForNuclides namespace

// Helper functions find_valid_energy_range, rois_are_similar, compute_roi_chi2_significance,
// compute_filtered_chi2_dof, and fit_peaks_for_nuclide_relactauto have been moved to FitPeaksForNuclides namespace


// estimate_initial_rois_without_peaks, estimate_initial_rois_fallback, and estimate_initial_rois_using_relactmanual
// have been moved to FitPeaksForNuclides namespace
// Functions removed - implementations are now in src/FitPeaksForNuclides.cpp


/** Test function for estimate_initial_rois_without_peaks.

 Tests the function with Pu238, Pu239, and Pu241 sources and validates:
 - ROIs are created for each source
 - No ROIs overlap
 - All ROIs are at least 1 FWHM wide
 - All ROIs contain their source gamma energies

 @return true if all tests pass, false otherwise
 */
bool test_estimate_initial_rois_without_peaks()
{
  cout << "\n=== Testing estimate_initial_rois_without_peaks ===" << endl;

  bool all_tests_passed = true;

  try
  {
    // Setup test sources: Pu238, Pu239, Pu241
    std::vector<RelActCalcAuto::NucInputInfo> test_sources;

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
    {
      cerr << "TEST FAILED: Could not get decay database" << endl;
      return false;
    }

    const SandiaDecay::Nuclide * const pu238 = db->nuclide( "Pu238" );
    const SandiaDecay::Nuclide * const pu239 = db->nuclide( "Pu239" );
    const SandiaDecay::Nuclide * const pu241 = db->nuclide( "Pu241" );

    if( !pu238 || !pu239 || !pu241 )
    {
      cerr << "TEST FAILED: Could not find Pu238, Pu239, or Pu241 in database" << endl;
      return false;
    }

    RelActCalcAuto::NucInputInfo pu238_input;
    pu238_input.source = pu238;
    pu238_input.age = 20.0*PhysicalUnits::year;
    pu238_input.fit_age = false;

    RelActCalcAuto::NucInputInfo pu239_input;
    pu239_input.source = pu239;
    pu239_input.age = 20.0*PhysicalUnits::year;
    pu239_input.fit_age = false;

    RelActCalcAuto::NucInputInfo pu241_input;
    pu241_input.source = pu241;
    pu241_input.age = 20.0*PhysicalUnits::year;
    pu241_input.fit_age = false;


    test_sources.push_back( pu238_input );
    test_sources.push_back( pu239_input );
    test_sources.push_back( pu241_input );

    // Setup test parameters
    const bool isHPGe = true;
    const DetectorPeakResponse::ResolutionFnctForm fwhmFnctnlForm
        = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;

    // Typical HPGe FWHM coefficients: FWHM = sqrt(a + b*E + c*E^2) where E is in keV
    const std::vector<float> fwhm_coefficients = {0.5f, 0.0f, 0.0f};  // ~0.7 keV FWHM constant
    const double lower_fwhm_energy = 0.0;
    const double upper_fwhm_energy = 3000.0;
    const double min_valid_energy = 50.0;
    const double max_valid_energy = 2000.0;

    // Create test config
    PeakFitForNuclideConfig config;
    config.auto_rel_eff_sol_max_fwhm = 12.0;

    // Get generic HPGe DRF
    std::shared_ptr<const DetectorPeakResponse> drf = DetectorPeakResponse::getGenericHPGeDetector();
    if( !drf || !drf->isValid() )
    {
      cerr << "TEST FAILED: Could not get generic HPGe detector" << endl;
      return false;
    }

    // Call the function under test
    const std::vector<RelActCalcAuto::RoiRange> rois = FitPeaksForNuclides::estimate_initial_rois_without_peaks(
      test_sources, drf, isHPGe,
      fwhmFnctnlForm, fwhm_coefficients, lower_fwhm_energy, upper_fwhm_energy,
      min_valid_energy, max_valid_energy, config );

    cout << "Created " << rois.size() << " ROIs" << endl;

    // Test 1: Check that we have ROIs
    if( rois.empty() )
    {
      cerr << "TEST FAILED: No ROIs created" << endl;
      all_tests_passed = false;
    }
    else
    {
      cout << "PASS: ROIs were created" << endl;
    }

    // Test 2: Check no overlaps
    for( size_t i = 1; i < rois.size(); ++i )
    {
      if( rois[i].lower_energy < rois[i-1].upper_energy )
      {
        cerr << "TEST FAILED: ROIs " << (i-1) << " and " << i << " overlap: "
             << "[" << rois[i-1].lower_energy << ", " << rois[i-1].upper_energy << "] "
             << "and [" << rois[i].lower_energy << ", " << rois[i].upper_energy << "]" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: No overlapping ROIs" << endl;

    // Test 3: Check that each ROI is at least 1 FWHM wide
    for( size_t i = 0; i < rois.size(); ++i )
    {
      const double roi_width = rois[i].upper_energy - rois[i].lower_energy;
      const double mid_energy = 0.5 * (rois[i].lower_energy + rois[i].upper_energy);
      const double fwhm = DetectorPeakResponse::peakResolutionFWHM(
        static_cast<float>(mid_energy), fwhmFnctnlForm, fwhm_coefficients );

      if( roi_width < fwhm )
      {
        cerr << "TEST FAILED: ROI " << i << " is narrower than 1 FWHM: "
             << roi_width << " keV < " << fwhm << " keV FWHM at " << mid_energy << " keV" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs are at least 1 FWHM wide" << endl;

    // Test 4: Check that we have ROIs for each source
    // Collect expected gamma energies for each source (top 4 by BR*eff)
    std::map<const SandiaDecay::Nuclide*, std::vector<double>> expected_gammas;

    for( const RelActCalcAuto::NucInputInfo &src : test_sources )
    {
      // get_source_age and get_source_photons are now in anonymous namespace in FitPeaksForNuclides.cpp
      // For test purposes, use a simple approach
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( src.source );
      const double age = nuc ? src.age : 0.0;
      SandiaDecay::NuclideMixture mix;
      if( nuc )
        mix.addAgedNuclideByActivity( nuc, 1.0, age );
      const std::vector<SandiaDecay::EnergyRatePair> photons = mix.photons( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy );

      struct GammaScore
      {
        double energy;
        double score;
      };
      std::vector<GammaScore> candidates;

      for( const SandiaDecay::EnergyRatePair &photon : photons )
      {
        if( photon.energy < min_valid_energy || photon.energy > max_valid_energy )
          continue;

        const double br = photon.numPerSecond;
        const double eff = drf->intrinsicEfficiency( static_cast<float>(photon.energy) );
        const double score = br * eff;

        if( score > 0.0 )
          candidates.push_back( {photon.energy, score} );
      }

      // Sort and take top 4
      std::sort( candidates.begin(), candidates.end(),
        [](const GammaScore &a, const GammaScore &b) { return a.score > b.score; } );

      const size_t num_to_take = std::min( candidates.size(), size_t(4) );

      for( size_t i = 0; i < num_to_take; ++i )
        expected_gammas[nuc].push_back( candidates[i].energy );
    }

    // Check that ROIs cover expected gammas
    for( const auto &[nuc, gammas] : expected_gammas )
    {
      cout << "Checking " << nuc->symbol << " gammas:" << endl;

      for( const double gamma_energy : gammas )
      {
        bool found = false;
        for( const RelActCalcAuto::RoiRange &roi : rois )
        {
          if( gamma_energy >= roi.lower_energy && gamma_energy <= roi.upper_energy )
          {
            found = true;
            cout << "  " << gamma_energy << " keV: FOUND in ROI [" << roi.lower_energy
                 << ", " << roi.upper_energy << "]" << endl;
            break;
          }
        }

        if( !found )
        {
          cerr << "TEST FAILED: " << nuc->symbol << " gamma at " << gamma_energy
               << " keV not found in any ROI" << endl;
          all_tests_passed = false;
        }
      }
    }

    // Test 5: Verify all ROIs have Linear continuum type
    for( size_t i = 0; i < rois.size(); ++i )
    {
      if( rois[i].continuum_type != PeakContinuum::OffsetType::Linear )
      {
        cerr << "TEST FAILED: ROI " << i << " does not have Linear continuum type" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs have Linear continuum type" << endl;

    // Test 6: Verify all ROIs have Fixed range limits
    for( size_t i = 0; i < rois.size(); ++i )
    {
      if( rois[i].range_limits_type != RelActCalcAuto::RoiRange::RangeLimitsType::Fixed )
      {
        cerr << "TEST FAILED: ROI " << i << " does not have Fixed range limits" << endl;
        all_tests_passed = false;
      }
    }

    if( all_tests_passed )
      cout << "PASS: All ROIs have Fixed range limits" << endl;

  }
  catch( const std::exception &e )
  {
    cerr << "TEST FAILED: Exception thrown: " << e.what() << endl;
    all_tests_passed = false;
  }

  if( all_tests_passed )
    cout << "=== ALL TESTS PASSED ===" << endl;
  else
    cout << "=== SOME TESTS FAILED ===" << endl;

  cout << endl;

  return all_tests_passed;
}//test_estimate_initial_rois_without_peaks


PeakFitResult fit_peaks_for_nuclides(
  const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &foreground,
  const std::vector<RelActCalcAuto::SrcVariant> &sources,
  const std::shared_ptr<const SpecUtils::Measurement> &long_background,
  const std::shared_ptr<const DetectorPeakResponse> &drf_input,
  const PeakFitForNuclideConfig &config,
  const bool isHPGe )
{
  // Implementation moved to FitPeaksForNuclides namespace
  // This is now a wrapper that calls the new namespace version
  return FitPeaksForNuclides::fit_peaks_for_nuclides(
    auto_search_peaks, foreground, sources, long_background, drf_input, config, isHPGe );
}//fit_peaks_for_nuclides


// Struct to hold local minimum info for ROI breakpoint selection
struct LocalMinimum {
  size_t channel;                  // Absolute channel number of minimum
  double synthetic_value;
  double depth_score;              // For tiebreaking (higher = better)
  double statistical_significance; // Primary criterion (lower = better breakpoint)
};


// build_synthetic_spectrum has been moved to FitPeaksForNuclides namespace
// Removing duplicate implementation

// compute_significance_in_region has been moved to FitPeaksForNuclides namespace
// Removing duplicate implementation

// has_significant_peak_between has been moved to FitPeaksForNuclides namespace
// Removing duplicate implementation

// find_synthetic_minima has been moved to FitPeaksForNuclides namespace
// Removing duplicate implementation


// should_use_step_continuum and cluster_gammas_to_rois have been moved to FitPeaksForNuclides namespace
// Functions removed - implementations are now in src/FitPeaksForNuclides.cpp




/** Comprehensive peak fitting function for nuclides.

 This function performs the complete peak fitting workflow:
 1. Determines FWHM functional form from auto-search peaks or DRF
 2. Matches auto-search peaks to source nuclides
 3. Performs RelActManual fit to get initial relative efficiency estimate
 4. Clusters gamma lines into ROIs based on manual rel-eff
 5. Calls fit_peaks_for_nuclide_relactauto with three different configurations
 6. Returns the best result based on chi2

 @param auto_search_peaks Initial peaks found by automatic search
 @param foreground Foreground spectrum to fit
 @param source_nuclides Vector of source nuclides to fit
 @param long_background Long background spectrum (can be nullptr)
 @param drf Detector response function (can be nullptr, will use generic if needed)
 @param config Configuration for peak fitting parameters
 @param isHPGe Whether this is an HPGe detector
 @return PeakFitResult with status, error message, fit peaks, and solution
 */
// fit_peaks_for_nuclides is now defined in FitPeaksForNuclides namespace
// (forward declaration removed - implementation is in namespace)


// This is development playground for `AnalystChecks::fit_peaks_for_nuclides(...)`
void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  // Configuration for peak fitting - these values will eventually be optimized via GA
  PeakFitForNuclideConfig config;
  // TODO: lets better assign FWHM form that should be used based on detector type
  // TODO: add an enum for how to handle existing ROIs here {ignore, replace_source, no_overlap}

  // Get clustering settings from config
  const GammaClusteringSettings manual_settings = config.get_manual_clustering_settings();
  const GammaClusteringSettings auto_settings = config.get_auto_clustering_settings();


  double total_score = 0.0;

  const bool write_n42 = true;

  // Initialize HTML output file (only if we're writing N42 files)
  std::unique_ptr<ofstream> html_output;
  size_t chart_counter = 0;  // For unique div IDs

  if( write_n42 )
  {
    html_output = std::make_unique<ofstream>( "relact_auto_peakfit_dev.html" );
    try
    {
      D3SpectrumExport::write_html_page_header( *html_output, "RelActAuto Peak Fits", "InterSpec_resources" );
    }catch( std::exception &e )
    {
      cerr << "\n\nError: " << e.what() << endl;
      cerr << "You probably need to symbolically link InterSpec_resources to your CWD." << endl;
      exit(1);
    }

    *html_output << "<body>" << endl;
    *html_output << "<style>"
      << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
      << "table, th, td{ border: 1px solid black; }" << endl
      << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
      << "</style>" << endl;
  }

  for( const DataSrcInfo &src_info : srcs_info )
  {
    const InjectSourceInfo &src = src_info.src_info;
    assert( src.src_spectra.size() >= 12 );
    if( src.src_spectra.size() < 12 )
      throw runtime_error( "Unexpected number of measurements." );

    string src_name = src.src_name;
    const auto underscore_pos = src_name.find( '_' );
    if( underscore_pos != string::npos )
    {
      const string last_part = src_name.substr(underscore_pos);
      assert( (last_part == "_Sh")  || (last_part == "_Sh-Point") || (last_part == "_Sh_BPE") || (last_part == "_Sh_PE")
             || (last_part == "_Sh_100g") || (last_part == "_ShHeavy") || (last_part == "_Unsh")
             || (last_part == "_Unsh-Point") || (last_part == "_Phantom")
             || (last_part == "_Unsh_0200") || (last_part == "_Unsh_0300") || (last_part == "_Unsh_0400")
             || (last_part == "_Unsh_0500") || (last_part == "_Unsh_1000") || (last_part == "_Unsh_2000")
             || (last_part == "_Unsh_4000") || (last_part == "_Unsh_5000") || (last_part == "_Unsh_6000")
             || (last_part == "_Unsh_8000") || (last_part == "_Unsh_9000") || (last_part == "_Unsh_9330")
      );
      src_name = src_name.substr(0,underscore_pos);
    }//if( underscore_pos != string::npos )

    vector<RelActCalcAuto::SrcVariant> sources;
    if( src_name == "Tl201woTl202" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202")); //I think it does actually have Tl202 in it
    }else if( src_name == "Tl201wTl202" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202"));
    }else if( src_name == "Tl201" )
    {
      sources.push_back( db->nuclide("Tl201"));
      sources.push_back( db->nuclide("Tl202"));
    }else if( src_name == "I125" )
    {
      sources.push_back( db->nuclide("I125"));
      sources.push_back( db->nuclide("I126"));
    }else if( src_name == "U233" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U232"));
      sources.push_back( db->nuclide("U233"));
    }else if( src_name == "Pu238" )
    {
      sources.push_back( db->element( "Pu" ) );
      sources.push_back( db->nuclide("Pu238"));
      sources.push_back( db->nuclide("Pu239"));
      sources.push_back( db->nuclide("Pu241"));
    }else if( src_name == "Pu239" )
    {
      sources.push_back( db->element( "Pu" ) );
      sources.push_back( db->nuclide("Pu239"));
      sources.push_back( db->nuclide("Pu241"));
    }else if( src_name == "Uore" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U235") );
      sources.push_back( db->nuclide("U238") );
      sources.push_back( db->nuclide("Ra226") );
    }else if( src_name == "U235" )
    {
      sources.push_back( db->element( "U" ) );
      sources.push_back( db->nuclide("U235") );
      sources.push_back( db->nuclide("U238") );
    }else if( src_name == "Np237" )
    {
      sources.push_back( db->element( "Np" ) );
      sources.push_back( db->nuclide("Np237") );
    }else if( src_name == "Xe133" )
    {
      sources.push_back( db->nuclide("Xe133") );
      sources.push_back( db->nuclide("Xe133m") );
    }else if( src_name == "Cf252" )
    {
      sources.push_back( db->nuclide("Cf252") );
      
      const ReactionGamma * const rctn_db = ReactionGammaServer::database();
      vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      if( rctn_db )
        rctn_db->gammas( "H(n,g)", possible_rctns );  //May throw exception if reaction name is invalid
      assert( !possible_rctns.empty() );
      if( !possible_rctns.empty() )
        sources.push_back( possible_rctns[0].reaction );
      
      continue;
    }else if( src_name == "Am241Li" )
    {
      continue;
    }else
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( src_name );
      if( !nuc )
      {
        cout << "Unable to get nuclide for '" << src.src_name << "' - so will skip" << endl;
        continue;
      }
      sources.push_back( nuc );
      assert( nuc->symbol == src_name );
    }

    assert( !sources.empty() && !RelActCalcAuto::is_null( sources.front() ) );
    if( sources.empty() || RelActCalcAuto::is_null( sources.front() ) )
      throw runtime_error( "Failed to get a source" );

    const shared_ptr<SpecUtils::SpecFile> &spec_file = src.spec_file;
    const shared_ptr<const SpecUtils::Measurement> &foreground = src.src_spectra.front(); // TODO: we cold loop over all 16 of these histograms
    const shared_ptr<const SpecUtils::Measurement> &long_background = src.long_background;

    const bool isHPGe = true;
    const bool singleThreaded = false;
    shared_ptr<const DetectorPeakResponse> drf = nullptr; //Probably

    try
    {
      shared_ptr<const deque<shared_ptr<const PeakDef>>> dummy_origpeaks;

      const vector<shared_ptr<const PeakDef> > auto_search_peaks
          = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, drf, dummy_origpeaks, singleThreaded, isHPGe );

      const PeakFitResult curve_results = FitPeaksForNuclideDev::fit_peaks_for_nuclides(
        auto_search_peaks, foreground, sources, long_background, drf, config, isHPGe );

      if( curve_results.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
        throw std::runtime_error( "fit_peaks_for_nuclides failed for all RelEff curve types: " + curve_results.error_message );

      const RelActCalcAuto::RelActAutoSolution &solution = curve_results.solution;
      const std::vector<PeakDef> &fit_peaks = curve_results.observable_peaks;

#if( PERFORM_DEVELOPER_CHECKS )
      {// Begin check for overlapping ROIs in fit_peaks
        // Collect unique continuums (ROIs) and their bounds
        std::map<std::shared_ptr<const PeakContinuum>, std::vector<double>> continuum_to_peak_means;
        for( const PeakDef &peak : fit_peaks )
          continuum_to_peak_means[peak.continuum()].push_back( peak.mean() );

        struct RoiCheckInfo
        {
          double lower_energy;
          double upper_energy;
          std::vector<double> peak_means;
        };

        std::vector<RoiCheckInfo> rois_to_check;
        for( const auto &entry : continuum_to_peak_means )
        {
          RoiCheckInfo info;
          info.lower_energy = entry.first->lowerEnergy();
          info.upper_energy = entry.first->upperEnergy();
          info.peak_means = entry.second;
          rois_to_check.push_back( info );
        }

        // Sort ROIs by lower energy
        std::sort( begin(rois_to_check), end(rois_to_check),
          []( const RoiCheckInfo &a, const RoiCheckInfo &b ) { return a.lower_energy < b.lower_energy; } );

        // Check for overlaps
        bool found_overlap = false;
        for( size_t i = 1; i < rois_to_check.size(); ++i )
        {
          const RoiCheckInfo &prev_roi = rois_to_check[i - 1];
          const RoiCheckInfo &curr_roi = rois_to_check[i];

          if( curr_roi.lower_energy < prev_roi.upper_energy )
          {
            if( !found_overlap )
              cerr << "\n*** OVERLAPPING ROIs DETECTED in fit_peaks (observable_peaks) ***" << endl;
            found_overlap = true;

            cerr << "  ROI[" << (i-1) << "]: [" << prev_roi.lower_energy << ", " << prev_roi.upper_energy << "] keV, peaks at: ";
            for( double mean : prev_roi.peak_means )
              cerr << mean << " ";
            cerr << "keV" << endl;

            cerr << "  ROI[" << i << "]: [" << curr_roi.lower_energy << ", " << curr_roi.upper_energy << "] keV, peaks at: ";
            for( double mean : curr_roi.peak_means )
              cerr << mean << " ";
            cerr << "keV" << endl;
            cerr << "  OVERLAP: " << (prev_roi.upper_energy - curr_roi.lower_energy) << " keV" << endl;
          }
        }

        if( found_overlap )
          cerr << "*** END OVERLAPPING ROIs ***\n" << endl;
      }// End check for overlapping ROIs in fit_peaks
#endif

      cout << "Best RelEff curve type solution (" << solution.m_options.rel_eff_curves.front().name << ") selected with chi2/dof="
           << solution.m_chi2 << "/" << solution.m_dof << endl;

      // Score the fit results using only signal photopeaks
      CombinedPeakFitScore combined_score;

      // Calculate FinalFitScore using refactored helper
      combined_score.final_fit_score = FinalFit_GA::calculate_final_fit_score(
        fit_peaks,
        src_info.expected_signal_photopeaks,
        config.num_sigma_contribution
      );

      // Calculate PeakFindAndFitWeights using refactored helper
      combined_score.initial_fit_weights = InitialFit_GA::calculate_peak_find_weights(
        fit_peaks,
        src_info.expected_signal_photopeaks,
        config.num_sigma_contribution
      );

      // Calculate CandidatePeakScore using refactored helper
      combined_score.candidate_peak_score = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
        fit_peaks,
        src_info.expected_signal_photopeaks
      );

      // Correct score to exclude escape peaks from penalty counts
      CandidatePeak_GA::correct_score_for_escape_peaks(
        combined_score.candidate_peak_score,
        src_info.expected_signal_photopeaks
      );

      // Calculate combined final weight (simple sum)
      combined_score.final_weight = combined_score.initial_fit_weights.find_weight
                                    + combined_score.final_fit_score.total_weight 
                                    + combined_score.candidate_peak_score.score;

      // Add to cumulative total score
      total_score += combined_score.final_weight;

      if( PeakFitImprove::debug_printout )
      {
        cout << "Combined score for " << src.src_name << ":" << endl;
        cout << combined_score.print( "combined_score" ) << endl;
      }

      // Write N42 file with foreground, background, and fit peaks (only if fit succeeded)
      if( curve_results.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
      {
        if( write_n42 )
        {
        string outdir = "output_n42_relactauto";
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.detector_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.location_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        outdir = SpecUtils::append_path( outdir, src_info.live_time_name );
        if( !SpecUtils::is_directory( outdir ) && !SpecUtils::create_directory( outdir ) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;

        // Use the source name for the output file
        const string out_n42 = SpecUtils::append_path( outdir, src.src_name ) + "_relactauto_fit.n42";

        SpecMeas output;
        output.set_manufacturer( src.spec_file->manufacturer() );
        output.set_detector_type( src.spec_file->detector_type() );
        output.set_instrument_model( src.spec_file->instrument_model() );
        output.set_instrument_type( src.spec_file->instrument_type() );
        output.set_filename( src.spec_file->filename() );

        output.remove_measurements( output.measurements() );

        // Generate a unique UUID based on current time + random component
        {
          const auto now_time = std::chrono::system_clock::now();
          const auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
              now_time.time_since_epoch() ).count();
          std::random_device rd;
          std::mt19937 gen( rd() );
          std::uniform_int_distribution<uint32_t> dis( 0, 0xFFFFFFFF );
          std::ostringstream uuid_oss;
          uuid_oss << std::hex << std::setfill('0')
                   << std::setw(12) << now_ms << "-"
                   << std::setw(8) << dis(gen) << "-"
                   << std::setw(8) << dis(gen);
          output.set_uuid( uuid_oss.str() );
        }

        output.add_remark( combined_score.print( "combined_score" ) );
        output.add_remark( "RelActAuto solution chi2/dof=" + std::to_string( solution.m_chi2 ) + "/" + std::to_string( solution.m_dof ) );

        output.set_instrument_model( "RelActAuto Fit Test" );

        // Add foreground spectrum
        shared_ptr<SpecUtils::Measurement> out_foreground = make_shared<SpecUtils::Measurement>( *foreground );
        out_foreground->set_sample_number( 1 );
        out_foreground->set_title( src_name + " - Foreground" );
        const auto now = std::chrono::time_point_cast<std::chrono::microseconds>( std::chrono::system_clock::now() );
        out_foreground->set_start_time( now );
        output.add_measurement( out_foreground, false );

        // Add background spectrum if present
        if( long_background )
        {
          shared_ptr<SpecUtils::Measurement> out_background = make_shared<SpecUtils::Measurement>( *long_background );
          out_background->set_sample_number( 2 );
          out_background->set_title( src_name + " - Background" );
          out_background->set_start_time( now );
          output.add_measurement( out_background, false );
        }

        // Add fit peaks associated with foreground
        deque<shared_ptr<const PeakDef>> peaks;
        for( const PeakDef &p : fit_peaks )
          peaks.push_back( make_shared<PeakDef>( p ) );

        output.setPeaks( peaks, {1} );

        // Add RelActAuto GUI state so the fit can be loaded in InterSpec
        RelActCalcAuto::RelActAutoGuiState gui_state;
        gui_state.options = solution.m_options;
        gui_state.background_subtract = (long_background != nullptr);
        gui_state.show_ref_lines = true;

        // Set display energy range from ROIs or foreground spectrum
        if( !solution.m_final_roi_ranges.empty() )
        {
          gui_state.lower_display_energy = solution.m_final_roi_ranges.front().lower_energy;
          gui_state.upper_display_energy = solution.m_final_roi_ranges.back().upper_energy;
        }
        else
        {
          gui_state.lower_display_energy = foreground->gamma_energy_min();
          gui_state.upper_display_energy = foreground->gamma_energy_max();
        }

        gui_state.note = "RelActAuto fit from PeakFitImprove";
        gui_state.description = "Source: " + src.src_name;

        output.setRelActAutoGuiState( &gui_state );

        output.save2012N42File( out_n42, [=](){
          cerr << "Failed to write '" << out_n42 << "'" << endl;
        });

        if( PeakFitImprove::debug_printout )
          cout << "Wrote N42 file: " << out_n42 << endl;

        // Generate reference lines for the nuclides that were fit
        std::map<std::string, std::string> reference_lines_json;

        // Define colors to rotate through for different nuclides
        const std::vector<Wt::WColor> colors = {
          Wt::WColor(0, 0, 139),      // darkBlue
          Wt::WColor(0, 139, 139),    // darkCyan
          Wt::WColor(0, 100, 0),      // darkGreen
          Wt::WColor(139, 0, 139),    // darkMagenta
          Wt::WColor(255, 140, 0),    // darkOrange
          Wt::WColor(139, 0, 0),      // darkRed
        };

        size_t color_index = 0;

        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          const std::vector<RelActCalcAuto::NuclideRelAct> &nuclides = solution.m_rel_activities[rel_eff_index];

          for( const RelActCalcAuto::NuclideRelAct &nuc_info : nuclides )
          {
            const std::string nuc_name = nuc_info.name();

            // Create RefLineInput for this nuclide
            RefLineInput input;
            input.m_input_txt = nuc_name;

            // Set age if it was fit
            if( nuc_info.age_was_fit && nuc_info.age > 0 )
            {
              input.m_age = std::to_string(nuc_info.age) + " s";
            }

            // Rotate through colors
            input.m_color = colors[color_index % colors.size()];
            color_index++;

            input.m_showGammas = true;
            input.m_showXrays = true;
            input.m_showAlphas = false;
            input.m_showBetas = false;
            input.m_promptLinesOnly = false;
            input.m_lower_br_cutt_off = 0.0; // Only show lines with BR > 1%

            // Generate reference line info
            std::shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );

            if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string ref_json;
              ref_info->toJson( ref_json );
              reference_lines_json[nuc_name] = ref_json;
            }
          }
        }

        {// Begin add reference lines for valid energy range
          const pair<double,double> valid_range = FitPeaksForNuclides::find_valid_energy_range( foreground );
          const double min_valid_energy = valid_range.first;
          const double max_valid_energy = valid_range.second;

          if( (min_valid_energy > 0.0) && (max_valid_energy > min_valid_energy) )
          {
            // Create RefLineInput for minimum energy
            RefLineInput min_input;
            min_input.m_input_txt = std::to_string(min_valid_energy) + " keV";
            min_input.m_color = Wt::WColor(255, 165, 0); // Orange
            min_input.m_showGammas = true;
            min_input.m_showXrays = false;
            min_input.m_showAlphas = false;
            min_input.m_showBetas = false;
            min_input.m_promptLinesOnly = false;
            min_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> min_ref_info = ReferenceLineInfo::generateRefLineInfo( min_input );
            if( min_ref_info && min_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string min_ref_json;
              min_ref_info->toJson( min_ref_json );
              reference_lines_json["Valid Range Min"] = min_ref_json;
            }

            // Create RefLineInput for maximum energy
            RefLineInput max_input;
            max_input.m_input_txt = std::to_string(max_valid_energy) + " keV";
            max_input.m_color = Wt::WColor(255, 165, 0); // Orange
            max_input.m_showGammas = true;
            max_input.m_showXrays = false;
            max_input.m_showAlphas = false;
            max_input.m_showBetas = false;
            max_input.m_promptLinesOnly = false;
            max_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> max_ref_info = ReferenceLineInfo::generateRefLineInfo( max_input );
            if( max_ref_info && max_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string max_ref_json;
              max_ref_info->toJson( max_ref_json );
              reference_lines_json["Valid Range Max"] = max_ref_json;
            }
          }
        }// End add reference lines for valid energy range

        {// Begin add red reference lines for each energy in `combined_score.candidate_peak_score.def_expected_but_not_detected`
          const std::vector<ExpectedPhotopeakInfo> &not_detected = combined_score.candidate_peak_score.def_expected_but_not_detected;
          for( size_t i = 0; i < not_detected.size(); ++i )
          {
            const double energy = not_detected[i].effective_energy;

            RefLineInput red_line_input;
            red_line_input.m_input_txt = std::to_string(energy) + " keV";
            red_line_input.m_color = Wt::WColor(255, 0, 0); // Red
            red_line_input.m_showGammas = true;
            red_line_input.m_showXrays = false;
            red_line_input.m_showAlphas = false;
            red_line_input.m_showBetas = false;
            red_line_input.m_promptLinesOnly = false;
            red_line_input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> red_ref_info = ReferenceLineInfo::generateRefLineInfo( red_line_input );
            if( red_ref_info && red_ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string red_ref_json;
              red_ref_info->toJson( red_ref_json );
              reference_lines_json["Missing Peak " + std::to_string(i+1)] = red_ref_json;
            }
          }
        }// End add red reference lines for each energy in `combined_score.candidate_peak_score.def_expected_but_not_detected`

        // Setup chart options
        std::string title = src_name + " - RelActAuto Fit";
        std::string dataTitle = "";
        bool useLogYAxis = true;
        bool showVerticalGridLines = false;
        bool showHorizontalGridLines = false;
        bool legendEnabled = true;
        bool compactXAxis = true;
        bool showPeakUserLabels = false;
        bool showPeakEnergyLabels = false;
        bool showPeakNuclideLabels = false;
        bool showPeakNuclideEnergyLabels = false;
        bool showEscapePeakMarker = false;
        bool showComptonPeakMarker = false;
        bool showComptonEdgeMarker = false;
        bool showSumPeakMarker = false;
        bool backgroundSubtract = false;

        // Set energy range
        float xMin = 0.0f;
        float xMax = 3000.0f;
        if( !solution.m_final_roi_ranges.empty() )
        {
          xMin = static_cast<float>( solution.m_final_roi_ranges.front().lower_energy );
          xMax = static_cast<float>( solution.m_final_roi_ranges.back().upper_energy );
        }
        else if( foreground )
        {
          xMin = static_cast<float>( foreground->gamma_energy_min() );
          xMax = static_cast<float>( foreground->gamma_energy_max() );
        }

        D3SpectrumExport::D3SpectrumChartOptions options(
          title, "Energy (keV)", "Counts/Channel",
          dataTitle, useLogYAxis,
          showVerticalGridLines, showHorizontalGridLines,
          legendEnabled, compactXAxis,
          showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels,
          showPeakNuclideEnergyLabels, showEscapePeakMarker, showComptonPeakMarker,
          showComptonEdgeMarker, showSumPeakMarker, backgroundSubtract,
          xMin, xMax, reference_lines_json
        );

        // Setup spectrum data with peaks
        D3SpectrumExport::D3SpectrumOptions foreground_opts;
        foreground_opts.line_color = "black";
        foreground_opts.title = src_name;
        foreground_opts.display_scale_factor = 1.0;
        foreground_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

        // Convert fit_peaks to shared_ptr vector for peak_json
        vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
        for( const PeakDef &p : fit_peaks )
          fit_peaks_ptrs.push_back( make_shared<PeakDef>( p ) );

        foreground_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, foreground, Wt::WColor(), false );

        // Create unique div ID
        const string div_id = "chart_" + std::to_string(chart_counter);
        chart_counter++;

        // Write HTML fieldset and div
        *html_output << "<fieldset style=\"\">" << endl
          << "<legend>" << src.src_name << " (chi2/dof=" << solution.m_chi2 << "/" << solution.m_dof << ")</legend>" << endl;
        *html_output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;
        *html_output << "<script>" << endl;

        // Initialize chart
        D3SpectrumExport::write_js_for_chart( *html_output, div_id, options.m_dataTitle, options.m_xAxisTitle, options.m_yAxisTitle );

        // Set spectrum data
        std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
        measurements.emplace_back( foreground.get(), foreground_opts );

        // Add background if present
        if( long_background )
        {
          D3SpectrumExport::D3SpectrumOptions background_opts;
          background_opts.line_color = "steelblue";
          background_opts.title = "Background";
          background_opts.display_scale_factor = foreground->live_time() / long_background->live_time();
          background_opts.spectrum_type = SpecUtils::SpectrumType::Background;
          measurements.emplace_back( long_background.get(), background_opts );
        }

        D3SpectrumExport::write_and_set_data_for_chart( *html_output, div_id, measurements );

        // Write resize handler
        *html_output << R"delim(
const resizeChart)delim" << chart_counter-1 << R"delim( = function(){
  let height = window.innerHeight;
  let width = document.documentElement.clientWidth;
  let el = spec_chart_)delim" << div_id << R"delim(.chart;
  el.style.width = 0.8*width + "px";
  el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
  el.style.marginLeft = 0.05*width + "px";
  el.style.marginRight = 0.05*width + "px";
  )delim"
          << "  spec_chart_" << div_id << R"delim(.handleResize();
};

window.addEventListener('resize', resizeChart)delim" << chart_counter-1 << R"delim();
)delim" << endl;

        // Set chart options (this creates reference_lines_<div_id> variable)
        D3SpectrumExport::write_set_options_for_chart( *html_output, div_id, options );

        // Apply reference lines
        *html_output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;

        // Trigger initial resize
        *html_output << "resizeChart" << chart_counter-1 << "();" << endl;
        *html_output << "</script>" << endl;

        // Add fit score information with three sections
        *html_output << "<div style=\"margin: 10px auto; max-width: 800px;\">" << endl;
        *html_output << "<h4>Fit Quality Metrics</h4>" << endl;

        // === Combined Score (prominent display) ===
        *html_output << "<div style=\"background: #f0f8ff; padding: 10px; margin-bottom: 15px; border: 2px solid #4682b4; border-radius: 5px;\">" << endl;
        *html_output << "<h5 style=\"margin: 0 0 10px 0;\">Combined Score</h5>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\"><strong>Final Weight (Combined)</strong></td>"
                     << "<td style=\"text-align: right; padding: 5px;\"><strong>"
                     << combined_score.final_weight << "</strong></td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 0: Peak Find Score</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.candidate_peak_score.score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 1: Find Weight</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.initial_fit_weights.find_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px;\">Component 2: Total Weight (Avrg nsigma off)</td>"
                     << "<td style=\"text-align: right; padding: 5px;\">"
                     << combined_score.final_fit_score.total_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</div>" << endl;

        // === Final Fit Metrics (collapsible) ===
        *html_output << "<details>" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Final Fit Quality Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Area Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.area_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Width Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.width_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Position Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.position_score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Ignored Unexpected Peaks</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.ignored_unexpected_peaks << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Unexpected Peaks Sum Significance</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.unexpected_peaks_sum_significance << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Number of Peaks Used</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.num_peaks_used << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Total Weight (Area)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.final_fit_score.total_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;

        // === Peak Finding Metrics (collapsible) ===
        *html_output << "<details style=\"margin-top: 10px;\">" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Peak Finding & Area Accuracy Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Find Weight (Score)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.find_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Def Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.def_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Def Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.def_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Maybe Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.maybe_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Maybe Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.maybe_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Not Wanted Area Weight (Mean)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.not_wanted_area_weight << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Not Wanted Area Weight (Median)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.initial_fit_weights.not_wanted_area_median_weight << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;

        // === Candidate Peak Metrics (collapsible) ===
        *html_output << "<details style=\"margin-top: 10px;\">" << endl;
        *html_output << "<summary style=\"cursor: pointer; font-weight: bold; padding: 5px; background: #f5f5f5; border: 1px solid #ddd;\">"
                     << "Candidate Peak Detection Metrics</summary>" << endl;
        *html_output << "<table style=\"width: 100%; border-collapse: collapse; font-size: 0.9em; margin-top: 5px;\">" << endl;
        *html_output << "<tr><th style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">Metric</th>"
                     << "<th style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">Value</th></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Score</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.score << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Number of Peaks Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_peaks_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Not Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_def_wanted_not_found << "</td></tr>" << endl;
        if( !combined_score.candidate_peak_score.def_expected_but_not_detected.empty() )
        {
          *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Not Found (Energies)</td>"
                       << "<td style=\"text-align: left; padding: 5px; border: 1px solid #ddd;\">";
          for( size_t i = 0; i < combined_score.candidate_peak_score.def_expected_but_not_detected.size(); ++i )
          {
            if( i > 0 )
              *html_output << ", ";
            *html_output << std::fixed << std::setprecision(2) 
                         << combined_score.candidate_peak_score.def_expected_but_not_detected[i].effective_energy << " keV";
          }
          *html_output << "</td></tr>" << endl;
        }
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Definitely Wanted Peaks Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_def_wanted_peaks_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Possibly Accepted Peaks Not Found</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_possibly_accepted_peaks_not_found << "</td></tr>" << endl;
        *html_output << "<tr><td style=\"padding: 5px; border: 1px solid #ddd;\">Extra Peaks (Not Expected)</td>"
                     << "<td style=\"text-align: right; padding: 5px; border: 1px solid #ddd;\">"
                     << combined_score.candidate_peak_score.num_extra_peaks << "</td></tr>" << endl;
        *html_output << "</table>" << endl;
        *html_output << "</details>" << endl;
        
        *html_output << "<p>Definitely Wanted Peaks Not Found: "
                     << combined_score.candidate_peak_score.num_def_wanted_not_found << "</p>" << endl;

        *html_output << "</div>" << endl;

        *html_output << "</fieldset>" << endl;
        }//if( write_n42 )
      }//if( curve_results.status == Success )

    }catch( const std::exception &e )
    {
      //throw std::runtime_error("Error in fit_peaks_for_nuclides: " + string(e.what()));
      cerr << e.what() << endl;
    }
  }//for( const DataSrcInfo &src_info : srcs_info )

  // Print total score summary
  cout << endl << "========================================" << endl;
  cout << "Total score across all " << srcs_info.size() << " sources: " << total_score << endl;
  cout << "========================================" << endl << endl;

  // Close HTML file
  if( html_output )
  {
    *html_output << "</body>" << endl;
    *html_output << "</html>" << endl;
    html_output->close();

    if( PeakFitImprove::debug_printout )
      cout << "Wrote HTML file: relact_auto_peakfit_dev.html with " << chart_counter << " charts" << endl;
  }
}//void eval_peaks_for_nuclide( const DataSrcInfo &src_info )


}//namespace FitPeaksForNuclideDev
