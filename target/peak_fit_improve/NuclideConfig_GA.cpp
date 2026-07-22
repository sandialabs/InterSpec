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

#include <deque>
#include <mutex>
#include <atomic>
#include <chrono>
#include <ctime>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <thread>
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PeakFit.h"

#include "NuclideConfig_GA.h"
#include "InitialFit_GA.h"
#include "FinalFit_GA.h"
#include "CandidatePeak_GA.h"
#include "ClassifyDetType_GA.h"
#include "PeakFitImproveData.h"
#include "FitPeaksForNuclideDev.h"

using namespace std;


namespace NuclideConfig_GA
{

// Runtime flags set from CLI in PeakFitImprove main().  Default values
// preserve the pre-CLI behavior of these features.
bool sm_do_background_fit_trial = false;
bool sm_disable_auto_interferer_fit = false;
RelEffChi2CapMode sm_rel_eff_chi2_cap_mode = RelEffChi2CapMode::Fixed;
double sm_rel_eff_chi2_cap_fixed_value = 25.0;
PeakFitUtils::CoarseResolutionType sm_base_det_type = PeakFitUtils::CoarseResolutionType::Low;

// Checkpoint/resume state, set from the CLI in main() (all opt-in).  See NuclideConfig_GA.h.
std::string sm_checkpoint_name;
std::string sm_resume_path;
std::string sm_checkpoint_options_summary;


namespace ReportDetail
{
namespace
{
std::pair<double,double> spectroscopic_extent( const PrecomputedNuclideData &pd )
{
  return pd.foreground ? FitPeaksForNuclides::find_valid_energy_range( pd.foreground )
                       : std::make_pair( 0.0, 0.0 );
}


std::vector<ExpectedPhotopeakInfo> filter_truth_for_extent(
  const std::vector<ExpectedPhotopeakInfo> &truth,
  const std::pair<double,double> &extent )
{
  if( !(extent.second > extent.first) )
    return truth;

  std::vector<ExpectedPhotopeakInfo> filtered;
  filtered.reserve( truth.size() );
  for( const ExpectedPhotopeakInfo &peak : truth )
    if( (peak.effective_energy >= extent.first) && (peak.effective_energy <= extent.second) )
      filtered.push_back( peak );
  return filtered;
}
}//namespace


std::vector<ExpectedPhotopeakInfo> filter_truth_for_spectroscopic_extent(
  const std::vector<ExpectedPhotopeakInfo> &truth,
  const PrecomputedNuclideData &pd )
{
  return filter_truth_for_extent( truth, spectroscopic_extent(pd) );
}

std::string canonical_spectrum_key( const DataSrcInfo &info )
{
  return info.detector_name + "/" + info.location_name + "/" + info.live_time_name
      + "/" + info.src_info.src_name + "|file=" + info.src_info.file_base_path;
}


std::string canonical_spectrum_key( const PrecomputedNuclideData &pd )
{
  if( !pd.src_info )
    return "missing-spectrum-metadata";

  return canonical_spectrum_key( *pd.src_info );
}


std::string stable_spectrum_id( const PrecomputedNuclideData &pd )
{
  const std::string key = canonical_spectrum_key( pd );
  uint64_t hash = 14695981039346656037ULL;
  for( const unsigned char c : key )
  {
    hash ^= static_cast<uint64_t>( c );
    hash *= 1099511628211ULL;
  }

  std::string slug;
  slug.reserve( key.size() );
  bool last_was_dash = false;
  for( const unsigned char c : key )
  {
    const bool allowed = std::isalnum( c );
    if( allowed )
    {
      slug.push_back( static_cast<char>(std::tolower(c)) );
      last_was_dash = false;
    }else if( !last_was_dash )
    {
      slug.push_back( '-' );
      last_was_dash = true;
    }
  }
  while( !slug.empty() && slug.back() == '-' )
    slug.pop_back();
  if( slug.size() > 52 )
    slug.resize( 52 );

  std::ostringstream answer;
  answer << "spectrum-" << std::hex << std::setw(16) << std::setfill('0') << hash
         << "-" << slug;
  return answer.str();
}


std::string html_escape( const std::string &text )
{
  std::string answer;
  answer.reserve( text.size() + text.size()/8 );
  for( const char c : text )
  {
    switch( c )
    {
      case '&':  answer += "&amp;";  break;
      case '<':  answer += "&lt;";   break;
      case '>':  answer += "&gt;";   break;
      case '"': answer += "&quot;"; break;
      case '\'': answer += "&#39;";  break;
      default:   answer.push_back( c ); break;
    }
  }
  return answer;
}


size_t roi_channel_count( const SpecUtils::Measurement &measurement,
                          const double lower_energy, const double upper_energy )
{
  if( (upper_energy <= lower_energy) || (measurement.num_gamma_channels() == 0) )
    return 0;
  const size_t lower_channel = measurement.find_gamma_channel( static_cast<float>(lower_energy) );
  const size_t upper_channel = measurement.find_gamma_channel( static_cast<float>(upper_energy) );
  return (upper_channel >= lower_channel) ? (1 + upper_channel - lower_channel) : 0;
}


double roi_width_in_fwhm( const double lower_energy, const double upper_energy,
                          const double representative_fwhm )
{
  return ((upper_energy > lower_energy) && (representative_fwhm > 0.0))
    ? ((upper_energy - lower_energy) / representative_fwhm) : 0.0;
}


void validate_evaluation_coverage( const size_t selected_spectra,
                                   const ConfigEvaluation &evaluation )
{
  if( evaluation.spectra.size() != selected_spectra )
    throw std::runtime_error( "Evaluation does not contain one record per selected spectrum" );
  if( (evaluation.successes + evaluation.legitimate_empties
       + evaluation.mechanical_failures) != selected_spectra )
    throw std::runtime_error( "Evaluation status counts do not cover every selected spectrum" );
  std::set<std::string> spectrum_ids;
  std::set<std::string> anchor_ids;
  for( const SpectrumEvaluation &spectrum : evaluation.spectra )
  {
    if( spectrum.spectrum_id.empty() || spectrum.anchor_id.empty() )
      throw std::runtime_error( "Evaluation contains a spectrum without stable identifiers" );
    if( !spectrum_ids.insert(spectrum.spectrum_id).second
        || !anchor_ids.insert(spectrum.anchor_id).second )
      throw std::runtime_error( "Evaluation contains duplicate spectrum identifiers" );
  }
}


FitAccuracyBreakdown score_observable_peaks(
  const PrecomputedNuclideData &pd,
  const std::vector<PeakDef> &observable_peaks )
{
  constexpr double num_sigma_contribution = 1.5;
  const PeakFitUtils::CoarseResolutionType det_type = pd.src_info->det_type;
  const std::vector<ExpectedPhotopeakInfo> detector_scoring_peaks
    = PeakFitImproveData::filter_photopeaks_for_scoring(
        pd.src_info->expected_signal_photopeaks, det_type );
  const std::vector<ExpectedPhotopeakInfo> scoring_peaks
    = filter_truth_for_spectroscopic_extent( detector_scoring_peaks, pd );

  CombinedPeakFitScore combined_score;
  combined_score.final_fit_score = FinalFit_GA::calculate_final_fit_score(
    observable_peaks, scoring_peaks, num_sigma_contribution );
  combined_score.initial_fit_weights = InitialFit_GA::calculate_peak_find_weights(
    observable_peaks, scoring_peaks, num_sigma_contribution, det_type );
  combined_score.candidate_peak_score
    = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
        observable_peaks, scoring_peaks, det_type );
  CandidatePeak_GA::correct_score_for_escape_peaks(
    combined_score.candidate_peak_score, scoring_peaks );

  FitAccuracyBreakdown answer;
  answer.area_cost = combined_score.final_fit_score.total_weight;
  answer.find_reward = combined_score.initial_fit_weights.find_weight;
  answer.candidate_reward = combined_score.candidate_peak_score.score;
  answer.miss_fraction = PeakFitImproveData::missed_def_wanted_area_fraction(
    scoring_peaks, combined_score.candidate_peak_score.def_expected_but_not_detected,
    det_type );
  answer.missed_definitely_wanted
    = combined_score.candidate_peak_score.num_def_wanted_not_found;
  answer.extra_peaks = combined_score.candidate_peak_score.num_extra_peaks;
  answer.cost = answer.area_cost - answer.find_reward - answer.candidate_reward
      + sm_miss_penalty_weight*answer.miss_fraction;
  return answer;
}


void run_self_tests()
{
  const std::string escaped = html_escape( "a<&>\"'b" );
  if( escaped != "a&lt;&amp;&gt;&quot;&#39;b" )
    throw std::runtime_error( "HTML escaping regression" );

  DataSrcInfo info;
  info.detector_name = "Detector <A>";
  info.location_name = "Location";
  info.live_time_name = "300_seconds";
  info.src_info.src_name = "Cs137_Sh";
  PrecomputedNuclideData pd;
  pd.src_info = &info;
  const std::string first = stable_spectrum_id( pd );
  const std::string second = stable_spectrum_id( pd );
  if( first != second || first.find("spectrum-") != 0 )
    throw std::runtime_error( "stable spectrum ID regression" );

  DataSrcInfo other = info;
  other.location_name = "Other";
  pd.src_info = &other;
  if( stable_spectrum_id(pd) == first )
    throw std::runtime_error( "stable spectrum ID collision in fixture" );

  if( std::fabs(roi_width_in_fwhm(100.0, 112.0, 3.0) - 4.0) > 1.0e-12
      || roi_width_in_fwhm(112.0, 100.0, 3.0) != 0.0 )
    throw std::runtime_error( "ROI FWHM measurement regression" );

  // The no-peaks path is the authoritative metric used for both failures and valid empties.
  pd.src_info = &info;
  info.det_type = PeakFitUtils::CoarseResolutionType::High;
  const FitAccuracyBreakdown empty = score_observable_peaks( pd, {} );
  if( empty.cost != empty.area_cost - empty.find_reward - empty.candidate_reward
        + sm_miss_penalty_weight*empty.miss_fraction )
    throw std::runtime_error( "objective component consistency regression" );

  ConfigEvaluation coverage;
  coverage.spectra.resize( 3 );
  coverage.spectra[0].spectrum_id = "success";
  coverage.spectra[0].anchor_id = "success-anchor";
  coverage.spectra[1].spectrum_id = "empty";
  coverage.spectra[1].anchor_id = "empty-anchor";
  coverage.spectra[1].legitimate_empty = true;
  coverage.spectra[2].spectrum_id = "failure";
  coverage.spectra[2].anchor_id = "failure-anchor";
  coverage.spectra[2].mechanical_failure = true;
  coverage.successes = 1;
  coverage.legitimate_empties = 1;
  coverage.mechanical_failures = 1;
  validate_evaluation_coverage( 3, coverage );
  ConfigEvaluation duplicate = coverage;
  duplicate.spectra[2].anchor_id = duplicate.spectra[1].anchor_id;
  bool rejected_duplicate = false;
  try
  {
    validate_evaluation_coverage( 3, duplicate );
  }catch( const std::exception & )
  {
    rejected_duplicate = true;
  }
  if( !rejected_duplicate )
    throw std::runtime_error( "duplicate report anchor was not rejected" );
  bool rejected_missing = false;
  try
  {
    validate_evaluation_coverage( 4, coverage );
  }catch( const std::exception & )
  {
    rejected_missing = true;
  }
  if( !rejected_missing )
    throw std::runtime_error( "missing report spectrum was not rejected" );
}
}//namespace ReportDetail

// Module-level state for the GA (following InitialFit_GA pattern).
// The GA evaluator is wrapped to also return the foreground-only and raw
// background-penalty components, so per-individual breakdowns can be
// surfaced in reporting when sm_do_background_fit_trial is enabled.
static std::function<double( const PeakFitForNuclideConfig &, double *, double * )> ns_ga_eval_fcn;
static std::atomic<size_t> ns_num_evals_this_generation{0};
static bool sm_has_been_called = false;
static bool sm_set_best_genes = false;
static NuclideConfigSolution sm_best_genes;
static double sm_best_total_cost = std::numeric_limits<double>::max();
// Cached breakdown of the best-ever individual.  Only populated when
// sm_do_background_fit_trial is true (otherwise best_fg == best_total_cost).
static double sm_best_fg_cost = 0.0;
static double sm_best_bg_cost = 0.0;
static std::ofstream sm_output_file;

// Store precomputed data pointer for use by SO_report_generation
static const std::vector<PrecomputedNuclideData> *sm_precomputed_ptr = nullptr;


std::vector<RelActCalcAuto::SrcVariant> resolve_sources( const std::string &src_name )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Failed to open SandiaDecayDataBase" );

  // Strip shielding/config suffix from name
  string name = src_name;
  const auto underscore_pos = name.find( '_' );
  if( underscore_pos != string::npos )
    name = name.substr( 0, underscore_pos );

  vector<RelActCalcAuto::SrcVariant> sources;

  if( name == "Tl201woTl202" || name == "Tl201wTl202" || name == "Tl201" )
  {
    sources.push_back( db->nuclide( "Tl201" ) );
    sources.push_back( db->nuclide( "Tl202" ) );
  }
  else if( name == "I125" )
  {
    sources.push_back( db->nuclide( "I125" ) );
    sources.push_back( db->nuclide( "I126" ) );
  }
  else if( name == "U233" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U232" ) );
    sources.push_back( db->nuclide( "U233" ) );
  }
  else if( name == "Pu238" )
  {
    sources.push_back( db->element( "Pu" ) );
    sources.push_back( db->nuclide( "Pu238" ) );
    sources.push_back( db->nuclide( "Pu239" ) );
    sources.push_back( db->nuclide( "Pu241" ) );
  }
  else if( name == "Pu239" )
  {
    sources.push_back( db->element( "Pu" ) );
    sources.push_back( db->nuclide( "Pu239" ) );
    sources.push_back( db->nuclide( "Pu241" ) );
  }
  else if( name == "Uore" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U235" ) );
    sources.push_back( db->nuclide( "U238" ) );
    sources.push_back( db->nuclide( "Ra226" ) );
  }
  else if( name == "U235" )
  {
    sources.push_back( db->element( "U" ) );
    sources.push_back( db->nuclide( "U235" ) );
    sources.push_back( db->nuclide( "U238" ) );
  }
  else if( name == "Np237" )
  {
    sources.push_back( db->element( "Np" ) );
    sources.push_back( db->nuclide( "Np237" ) );
  }
  else if( name == "Xe133" )
  {
    sources.push_back( db->nuclide( "Xe133" ) );
    sources.push_back( db->nuclide( "Xe133m" ) );
  }
  else if( name == "Cf252" || name == "Am241Li" )
  {
    // Skip these sources
    return sources;
  }
  else
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
    if( !nuc )
    {
      cout << "Unable to get nuclide for '" << src_name << "' - so will skip" << endl;
      return sources;
    }
    sources.push_back( nuc );
  }

  return sources;
}//resolve_sources


std::vector<PrecomputedNuclideData> precompute_nuclide_data(
  const std::vector<DataSrcInfo> &srcs_info,
  const BackgroundMode bg_mode )
{
  std::vector<PrecomputedNuclideData> result;
  result.reserve( srcs_info.size() );

  cout << "Precomputing auto-search peaks for " << srcs_info.size() << " spectra..." << endl;

  // First, cheaply build up the entries (resolve sources, detector type, etc.),
  // skipping any that cannot be resolved.  The expensive search_for_peaks calls
  // are deferred and run in parallel below.  `result` is fully populated here and
  // never resized afterwards, so the worker tasks can safely write into their
  // own slot by index.
  for( size_t i = 0; i < srcs_info.size(); ++i )
  {
    const DataSrcInfo &info = srcs_info[i];
    const InjectSourceInfo &src = info.src_info;

    if( src.src_spectra.empty() )
      continue;

    // Resolve sources
    std::vector<RelActCalcAuto::SrcVariant> sources = resolve_sources( src.src_name );
    if( sources.empty() || RelActCalcAuto::is_null( sources.front() ) )
      continue;

    PrecomputedNuclideData pd;
    pd.src_info = &info;
    pd.sources = std::move( sources );
    pd.foreground = src.src_spectra.front();
    pd.drf = nullptr;

    // Set background based on mode
    switch( bg_mode )
    {
      case BackgroundMode::BackgroundSubtracted:
        pd.background = src.long_background;
        break;
      case BackgroundMode::NoBackground:
      case BackgroundMode::NoBackgroundFitNorm:
        pd.background = nullptr;
        break;
    }//switch( bg_mode )

    // Use det_type and skew_type from the DataSrcInfo (set during data loading)
    pd.det_type = info.det_type;

    // Create PeakFitDetPrefs with correct type
    pd.peak_fit_prefs = std::make_shared<PeakFitDetPrefs>();
    pd.peak_fit_prefs->m_det_type = pd.det_type;

    result.push_back( std::move( pd ) );
  }//for( size_t i = 0; i < srcs_info.size(); ++i )

  // Run the expensive search_for_peaks calls in parallel - one task per spectrum,
  // using the same worker count as the GA optimization.  Each task searches
  // single-threaded so we don't oversubscribe the CPU (num_threads * internal
  // search threads).
  boost::asio::thread_pool pool( PeakFitImprove::sm_num_optimization_threads );
  std::atomic<size_t> num_completed{ 0 };
  const size_t num_total = result.size();

  for( size_t idx = 0; idx < num_total; ++idx )
  {
    boost::asio::post( pool, [idx, num_total, &result, &num_completed](){
      PrecomputedNuclideData &pd = result[idx];

      const bool singleThreaded = true;
      std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> dummy_origpeaks;

      // The expensive call - done once
      pd.auto_search_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks(
        pd.foreground, pd.drf, dummy_origpeaks, singleThreaded, pd.peak_fit_prefs );

      // Background auto-search peaks for the GA background-false-positive penalty.
      // Gated on sm_do_background_fit_trial (CLI flag) - default off, so we skip
      // the expensive bg auto-search entirely unless the user opted in.  When
      // enabled we still only compute for entries that have a long_background
      // and at least one non-NORM-like source (otherwise the penalty path is a
      // no-op).
      const std::shared_ptr<const SpecUtils::Measurement> long_background
        = pd.src_info->src_info.long_background;
      if( sm_do_background_fit_trial && long_background )
      {
        bool any_non_norm = false;
        for( const RelActCalcAuto::SrcVariant &s : pd.sources )
        {
          if( !FitPeaksForNuclides::is_norm_like_for_ga( s ) )
          {
            any_non_norm = true;
            break;
          }
        }

        if( any_non_norm )
        {
          pd.background_auto_search_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks(
            long_background, pd.drf, dummy_origpeaks, singleThreaded, pd.peak_fit_prefs );
        }
      }//if( sm_do_background_fit_trial && long_background )

      const size_t done = num_completed.fetch_add( 1, std::memory_order_relaxed ) + 1;
      if( (done % 10) == 0 || done == num_total )
        cout << "  Precomputed " << done << " of " << num_total << " spectra" << endl;
    } );
  }//for( size_t idx = 0; idx < num_total; ++idx )

  pool.join();

  cout << "Precomputation complete: " << result.size() << " spectra ready for GA evaluation." << endl;

  return result;
}//precompute_nuclide_data


double compute_background_fit_penalty(
    const PrecomputedNuclideData &pd,
    const PeakFitForNuclideConfig &config,
    const BackgroundMode bg_mode,
    BackgroundFitDetail *detail_out )
{
  // Reset out-detail (caller may reuse).
  if( detail_out )
    *detail_out = BackgroundFitDetail{};

  // Skip if there's nothing to evaluate.
  if( pd.background_auto_search_peaks.empty() )
    return 0.0;

  const std::shared_ptr<const SpecUtils::Measurement> &long_bg
    = pd.src_info ? pd.src_info->src_info.long_background : nullptr;
  if( !long_bg )
    return 0.0;

  // The background fit treats long_background AS the foreground.  We
  // cannot subtract long_background from itself, so we pass nullptr for
  // the bg argument regardless of the global BackgroundMode.  We do mirror
  // NoBackgroundFitNorm by also fitting NORM peaks when that mode is set,
  // so the config is exercised consistently with the foreground pass.
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;
  if( bg_mode == BackgroundMode::NoBackgroundFitNorm )
    options |= FitPeaksForNuclides::FitNormBkgrndPeaks;
  if( sm_disable_auto_interferer_fit )
    options |= FitPeaksForNuclides::DisableAutoInterfererFit;

  const std::vector<std::shared_ptr<const PeakDef>> user_peaks;

  FitPeaksForNuclides::PeakFitResult result;
  try
  {
    result = FitPeaksForNuclides::fit_peaks_for_nuclides(
      pd.background_auto_search_peaks, long_bg, pd.sources, user_peaks,
      /*long_background=*/nullptr, pd.drf, options, config, pd.peak_fit_prefs );
  }catch( const std::exception &e )
  {
    if( detail_out )
    {
      detail_out->error_message = std::string("exception: ") + e.what();
      detail_out->penalty = sm_background_fit_penalty_per_spectrum_cap;
    }
    return sm_background_fit_penalty_per_spectrum_cap;
  }

  if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
  {
    if( detail_out )
    {
      detail_out->error_message = "fit status != Success: " + result.error_message;
      detail_out->penalty = sm_background_fit_penalty_per_spectrum_cap;
    }
    return sm_background_fit_penalty_per_spectrum_cap;
  }

  // Livetime normalization: a 3000s long_background gives ~sqrt(10)x the
  // significance of the same false-positive in a 300s foreground.  Divide
  // each peak's sig by sqrt(bg_lt / fg_lt) BEFORE the 7.5 cap.
  const double fg_lt = pd.foreground ? pd.foreground->live_time() : 1.0;
  const double bg_lt = long_bg->live_time();
  const double lt_ratio = (fg_lt > 0.0 && bg_lt > 0.0) ? (bg_lt / fg_lt) : 1.0;
  const double lt_norm = std::sqrt( std::max( lt_ratio, 1.0e-9 ) );

  // Per-peak filter thresholds:
  //   - peaks below low_energy_cutoff_kev are skipped (xrays / scatter are
  //     too easily generated in any spectrum to be a reliable signal)
  //   - peaks within annihilation_window_kev of 511.0 are skipped (ambient
  //     annihilation; sources with 511 in their gammas mis-attribute it)
  //   - peaks within ~1 FWHM of a strong NORM gamma are skipped (the bg
  //     peak is likely a real NORM peak the fit attributed to the source)
  constexpr double low_energy_cutoff_kev   = 30.0;
  constexpr double annihilation_kev        = 510.9989;
  constexpr double annihilation_window_kev = 2.0;

  double raw_penalty = 0.0;

  for( const PeakDef &peak : result.observable_peaks )
  {
    const SandiaDecay::Nuclide * const nuc = peak.parentNuclide();
    if( !nuc )
      continue;  // Floating / unattributed peaks (e.g. background 511) don't count.

    const double area = peak.amplitude();
    const double area_uncert = peak.amplitudeUncert();
    const double sig = (area_uncert > 0.0)
      ? (area / area_uncert)
      : ((area > 0.0) ? std::sqrt(area) : 0.0);
    const double sig_norm = sig / lt_norm;
    const double capped_sig = std::min( 7.5, std::max( 0.0, sig_norm ) );

    const double mean = peak.mean();
    const double fwhm = peak.fwhm();

    // Build suppression reason in priority order.
    std::string reason;
    if( FitPeaksForNuclides::is_norm_like_for_ga( RelActCalcAuto::SrcVariant( nuc ) ) )
      reason = "norm-like";
    else if( mean < low_energy_cutoff_kev )
      reason = "below 30 keV";
    else if( std::fabs( mean - annihilation_kev ) < annihilation_window_kev )
      reason = "near 511 keV";
    else if( FitPeaksForNuclides::is_near_strong_norm_gamma( mean, std::max(0.5, fwhm) ) )
      reason = "near NORM line";

    if( detail_out )
    {
      detail_out->source_attributed_peaks.push_back( peak );
      detail_out->normalized_significances.push_back( sig_norm );
      detail_out->suppression_reasons.push_back( reason );
    }

    if( reason.empty() )
      raw_penalty += capped_sig;
  }

  const double penalty = std::min( raw_penalty, sm_background_fit_penalty_per_spectrum_cap );

  if( detail_out )
    detail_out->penalty = penalty;

  return penalty;
}//compute_background_fit_penalty


std::string NuclideConfigSolution::to_string( const std::string &separator ) const
{
  std::ostringstream oss;
  // Full round-trip precision: the reformulated pipeline is decision-cascaded (an adaptive-extent
  // block accept/reject flip propagates through merging and the fit), so genes rounded to 6
  // significant digits could re-score several percent off their in-memory cost - which would make
  // checkpoint resume and NuclideConfigEval cross-checks unfaithful.
  oss << std::setprecision( std::numeric_limits<double>::max_digits10 );
  oss << "fwhm_functional_form=" << fwhm_functional_form << separator
      << "rel_eff_manual_base_rel_eff_uncert=" << rel_eff_manual_base_rel_eff_uncert << separator
      << "initial_nuc_match_cluster_num_sigma=" << initial_nuc_match_cluster_num_sigma << separator
      << "manual_eff_cluster_num_sigma=" << manual_eff_cluster_num_sigma << separator
      << "manual_releff_aicc_penalty=" << manual_releff_aicc_penalty << separator
      << "cont_order_aicc_penalty=" << cont_order_aicc_penalty << separator
      << "manual_keep_significance_z=" << manual_keep_significance_z << separator
      << "manual_rel_eff_sol_min_fwhm_roi=" << manual_rel_eff_sol_min_fwhm_roi << separator
      << "manual_rel_eff_sol_max_fwhm=" << manual_rel_eff_sol_max_fwhm << separator
      << "merge_tail_z=" << merge_tail_z << separator
      << "merge_clean_gap_fwhm=" << merge_clean_gap_fwhm << separator
      << "manual_roi_core_num_fwhm=" << manual_roi_core_num_fwhm << separator
      << "fwhm_form=" << fwhm_form << separator
      << "rel_eff_auto_base_rel_eff_uncert=" << rel_eff_auto_base_rel_eff_uncert << separator
      << "auto_rel_eff_cluster_num_sigma=" << auto_rel_eff_cluster_num_sigma << separator
      << "auto_keep_significance_z=" << auto_keep_significance_z << separator
      << "auto_roi_core_num_fwhm=" << auto_roi_core_num_fwhm << separator
      << "roi_extend_z=" << roi_extend_z << separator
      << "roi_max_num_fwhm=" << roi_max_num_fwhm << separator
      << "auto_rel_eff_sol_max_fwhm=" << auto_rel_eff_sol_max_fwhm << separator
      << "auto_rel_eff_sol_min_fwhm_roi=" << auto_rel_eff_sol_min_fwhm_roi << separator
      << "rel_eff_eqn_type=" << rel_eff_eqn_type << separator
      << "rel_eff_eqn_order=" << rel_eff_eqn_order << separator
      << "desperation_phys_model_atomic_number=" << desperation_phys_model_atomic_number << separator
      << "desperation_phys_model_areal_density_g_per_cm2=" << desperation_phys_model_areal_density_g_per_cm2 << separator
      << "nucs_of_el_same_age=" << nucs_of_el_same_age << separator
      << "phys_model_use_hoerl=" << phys_model_use_hoerl << separator
      << "fit_energy_cal=" << fit_energy_cal << separator
      << "roi_significance_z=" << roi_significance_z << separator
      << "observable_peak_initial_significance_threshold=" << observable_peak_initial_significance_threshold << separator
      << "observable_peak_final_significance_threshold=" << observable_peak_final_significance_threshold << separator
      << "step_cont_min_peak_significance=" << step_cont_min_peak_significance << separator
      << "step_trial_chi2_margin=" << step_trial_chi2_margin << separator
      << "initial_manual_rel_eff_max_chi2_dof=" << initial_manual_rel_eff_max_chi2_dof;
  return oss.str();
}//NuclideConfigSolution::to_string


bool NuclideConfigSolution::from_string( const std::string &line, const std::string &separator,
                                         NuclideConfigSolution &out )
{
  // Split "name=value<sep>name=value..." into a name->value map (order-independent; unknown keys are
  //  ignored so older checkpoints stay loadable as genes are added).
  std::map<std::string,std::string> kv;
  size_t pos = 0;
  while( true )
  {
    const size_t next = line.find( separator, pos );
    const size_t end = (next == std::string::npos) ? line.size() : next;
    const std::string token = line.substr( pos, end - pos );
    const size_t eq = token.find( '=' );
    if( eq != std::string::npos )
      kv[ token.substr( 0, eq ) ] = token.substr( eq + 1 );
    if( next == std::string::npos )
      break;
    pos = end + separator.size();
  }

  bool ok = true;
  const auto getI = [&kv,&ok]( const char *key, int &dst ){
    const auto it = kv.find( key );
    if( it == std::end(kv) ){ ok = false; return; }
    try{ dst = std::stoi( it->second ); }catch( ... ){ ok = false; }
  };
  const auto getD = [&kv,&ok]( const char *key, double &dst ){
    const auto it = kv.find( key );
    if( it == std::end(kv) ){ ok = false; return; }
    try{ dst = std::stod( it->second ); }catch( ... ){ ok = false; }
  };

  getI( "fwhm_functional_form", out.fwhm_functional_form );
  getD( "rel_eff_manual_base_rel_eff_uncert", out.rel_eff_manual_base_rel_eff_uncert );
  getD( "initial_nuc_match_cluster_num_sigma", out.initial_nuc_match_cluster_num_sigma );
  getD( "manual_eff_cluster_num_sigma", out.manual_eff_cluster_num_sigma );
  getD( "manual_releff_aicc_penalty", out.manual_releff_aicc_penalty );
  getD( "cont_order_aicc_penalty", out.cont_order_aicc_penalty );
  getD( "manual_keep_significance_z", out.manual_keep_significance_z );
  getD( "manual_rel_eff_sol_min_fwhm_roi", out.manual_rel_eff_sol_min_fwhm_roi );
  getD( "manual_rel_eff_sol_max_fwhm", out.manual_rel_eff_sol_max_fwhm );
  getD( "merge_tail_z", out.merge_tail_z );
  getD( "merge_clean_gap_fwhm", out.merge_clean_gap_fwhm );
  getD( "manual_roi_core_num_fwhm", out.manual_roi_core_num_fwhm );
  getI( "fwhm_form", out.fwhm_form );
  getD( "rel_eff_auto_base_rel_eff_uncert", out.rel_eff_auto_base_rel_eff_uncert );
  getD( "auto_rel_eff_cluster_num_sigma", out.auto_rel_eff_cluster_num_sigma );
  getD( "auto_keep_significance_z", out.auto_keep_significance_z );
  getD( "auto_roi_core_num_fwhm", out.auto_roi_core_num_fwhm );
  getD( "roi_extend_z", out.roi_extend_z );
  getD( "roi_max_num_fwhm", out.roi_max_num_fwhm );
  getD( "auto_rel_eff_sol_max_fwhm", out.auto_rel_eff_sol_max_fwhm );
  getD( "auto_rel_eff_sol_min_fwhm_roi", out.auto_rel_eff_sol_min_fwhm_roi );
  getI( "rel_eff_eqn_type", out.rel_eff_eqn_type );
  getI( "rel_eff_eqn_order", out.rel_eff_eqn_order );
  getD( "desperation_phys_model_atomic_number", out.desperation_phys_model_atomic_number );
  getD( "desperation_phys_model_areal_density_g_per_cm2", out.desperation_phys_model_areal_density_g_per_cm2 );
  getI( "nucs_of_el_same_age", out.nucs_of_el_same_age );
  getI( "phys_model_use_hoerl", out.phys_model_use_hoerl );
  getI( "fit_energy_cal", out.fit_energy_cal );
  getD( "roi_significance_z", out.roi_significance_z );
  getD( "observable_peak_initial_significance_threshold", out.observable_peak_initial_significance_threshold );
  getD( "observable_peak_final_significance_threshold", out.observable_peak_final_significance_threshold );
  getD( "step_cont_min_peak_significance", out.step_cont_min_peak_significance );
  getD( "step_trial_chi2_margin", out.step_trial_chi2_margin );
  getD( "initial_manual_rel_eff_max_chi2_dof", out.initial_manual_rel_eff_max_chi2_dof );

  return ok;
}//NuclideConfigSolution::from_string

// Resolves the effective initial_manual_rel_eff_max_chi2_dof based on the
// CLI mode.  In Disabled mode the cap is set high enough that the throw at
// src/FitPeaksForNuclides.cpp:5594 effectively never fires.  In Fixed mode
// the CLI-provided value is used.  In GAOptimized mode the gene's value is
// used directly.
static double resolve_chi2_dof_cap( const NuclideConfigSolution &p )
{
  switch( sm_rel_eff_chi2_cap_mode )
  {
    case RelEffChi2CapMode::Disabled:    return 1.0e6;
    case RelEffChi2CapMode::Fixed:       return sm_rel_eff_chi2_cap_fixed_value;
    case RelEffChi2CapMode::GAOptimized: return p.initial_manual_rel_eff_max_chi2_dof;
  }
  return sm_rel_eff_chi2_cap_fixed_value;
}//resolve_chi2_dof_cap


PeakFitForNuclideConfig genes_to_settings( const NuclideConfigSolution &p )
{
  // Start from the production default_config for the detector type this GA run targets
  // (sm_base_det_type, set in main()).  All GA-optimized fields are overwritten below; non-optimized
  // fields (e.g., skew_type, norm_css_color, shielding vectors) inherit that detector type's defaults
  // so the genome behaves the same way it will in production for this detector class.
  PeakFitForNuclideConfig config = PeakFitForNuclideConfig::default_config( sm_base_det_type );

  config.fwhm_functional_form = static_cast<DetectorPeakResponse::ResolutionFnctForm>(
    std::clamp( p.fwhm_functional_form, 0, static_cast<int>(DetectorPeakResponse::kNumResolutionFnctForm) - 1 ) );

  config.rel_eff_manual_base_rel_eff_uncert = p.rel_eff_manual_base_rel_eff_uncert;
  config.initial_nuc_match_cluster_num_sigma = p.initial_nuc_match_cluster_num_sigma;
  config.manual_eff_cluster_num_sigma = p.manual_eff_cluster_num_sigma;

  config.manual_releff_aicc_penalty = p.manual_releff_aicc_penalty;
  config.cont_order_aicc_penalty = p.cont_order_aicc_penalty;

  config.manual_keep_significance_z = p.manual_keep_significance_z;
  config.manual_rel_eff_sol_min_fwhm_roi = p.manual_rel_eff_sol_min_fwhm_roi;
  config.manual_rel_eff_sol_max_fwhm = p.manual_rel_eff_sol_max_fwhm;
  config.merge_tail_z = p.merge_tail_z;
  config.merge_clean_gap_fwhm = p.merge_clean_gap_fwhm;
  config.manual_roi_core_num_fwhm = p.manual_roi_core_num_fwhm;

  // FwhmForm: limited to Berstein_2 through Berstein_5
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
  config.fwhm_form = static_cast<RelActCalcAuto::FwhmForm>(
    std::clamp( p.fwhm_form,
                static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2),
                static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) ) );

  config.rel_eff_auto_base_rel_eff_uncert = p.rel_eff_auto_base_rel_eff_uncert;
  config.auto_rel_eff_cluster_num_sigma = p.auto_rel_eff_cluster_num_sigma;
  config.auto_keep_significance_z = p.auto_keep_significance_z;
  config.auto_roi_core_num_fwhm = p.auto_roi_core_num_fwhm;
  config.roi_extend_z = p.roi_extend_z;
  config.roi_max_num_fwhm = p.roi_max_num_fwhm;
  config.auto_rel_eff_sol_max_fwhm = p.auto_rel_eff_sol_max_fwhm;
  config.auto_rel_eff_sol_min_fwhm_roi = p.auto_rel_eff_sol_min_fwhm_roi;

  // RelEffEqnForm: 0=LnX, 1=LnY, 2=LnXLnY, 3=FramEmpirical, 4=FramPhysicalModel
  config.rel_eff_eqn_type = static_cast<RelActCalc::RelEffEqnForm>(
    std::clamp( p.rel_eff_eqn_type, 0, static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel) ) );
  config.rel_eff_eqn_order = static_cast<size_t>( std::clamp( p.rel_eff_eqn_order, 0, 6 ) );

  // FramPhysicalModel requires order 0; couple the genes here so the GA never samples an
  // (order>0, FramPhysicalModel) combination that would trip the engine's debug assert (and is
  // fitness-equivalent to order 0 in release, since the engine self-corrects).
  if( config.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    config.rel_eff_eqn_order = 0;

  config.desperation_phys_model_atomic_number = p.desperation_phys_model_atomic_number;
  config.desperation_phys_model_areal_density_g_per_cm2 = p.desperation_phys_model_areal_density_g_per_cm2;

  config.nucs_of_el_same_age = (p.nucs_of_el_same_age != 0);
  config.phys_model_use_hoerl = (p.phys_model_use_hoerl != 0);
  config.fit_energy_cal = (p.fit_energy_cal != 0);

  config.roi_significance_z = p.roi_significance_z;
  config.observable_peak_initial_significance_threshold = p.observable_peak_initial_significance_threshold;
  config.observable_peak_final_significance_threshold = p.observable_peak_final_significance_threshold;

  config.step_cont_min_peak_significance = p.step_cont_min_peak_significance;
  config.step_trial_chi2_margin = p.step_trial_chi2_margin;

  // Manual rel-eff chi2/dof cap: resolved from the CLI mode + (when
  // GA-optimized) the gene value.
  config.initial_manual_rel_eff_max_chi2_dof = resolve_chi2_dof_cap( p );

  // skew_type is not GA-optimized - leave at default (NoSkew)
  // Caller should override per detector type (e.g., DoubleSidedCrystalBall for CZT)

  return config;
}//genes_to_settings


NuclideConfigSolution settings_to_genes( const PeakFitForNuclideConfig &config )
{
  NuclideConfigSolution p;

  p.fwhm_functional_form = static_cast<int>( config.fwhm_functional_form );

  p.rel_eff_manual_base_rel_eff_uncert = config.rel_eff_manual_base_rel_eff_uncert;
  p.initial_nuc_match_cluster_num_sigma = config.initial_nuc_match_cluster_num_sigma;
  p.manual_eff_cluster_num_sigma = config.manual_eff_cluster_num_sigma;

  p.manual_releff_aicc_penalty = config.manual_releff_aicc_penalty;
  p.cont_order_aicc_penalty = config.cont_order_aicc_penalty;

  p.manual_keep_significance_z = config.manual_keep_significance_z;
  p.manual_rel_eff_sol_min_fwhm_roi = config.manual_rel_eff_sol_min_fwhm_roi;
  p.manual_rel_eff_sol_max_fwhm = config.manual_rel_eff_sol_max_fwhm;
  p.manual_roi_core_num_fwhm = config.manual_roi_core_num_fwhm;

  p.fwhm_form = static_cast<int>( config.fwhm_form );
  p.rel_eff_auto_base_rel_eff_uncert = config.rel_eff_auto_base_rel_eff_uncert;
  p.auto_rel_eff_cluster_num_sigma = config.auto_rel_eff_cluster_num_sigma;
  p.auto_keep_significance_z = config.auto_keep_significance_z;
  p.auto_roi_core_num_fwhm = config.auto_roi_core_num_fwhm;
  p.roi_extend_z = config.roi_extend_z;
  p.roi_max_num_fwhm = config.roi_max_num_fwhm;
  p.auto_rel_eff_sol_max_fwhm = config.auto_rel_eff_sol_max_fwhm;
  p.auto_rel_eff_sol_min_fwhm_roi = config.auto_rel_eff_sol_min_fwhm_roi;

  p.rel_eff_eqn_type = static_cast<int>( config.rel_eff_eqn_type );
  p.rel_eff_eqn_order = static_cast<int>( config.rel_eff_eqn_order );

  p.desperation_phys_model_atomic_number = config.desperation_phys_model_atomic_number;
  p.desperation_phys_model_areal_density_g_per_cm2 = config.desperation_phys_model_areal_density_g_per_cm2;

  p.nucs_of_el_same_age = config.nucs_of_el_same_age ? 1 : 0;
  p.phys_model_use_hoerl = config.phys_model_use_hoerl ? 1 : 0;
  p.fit_energy_cal = config.fit_energy_cal ? 1 : 0;

  p.merge_tail_z = config.merge_tail_z;
  p.merge_clean_gap_fwhm = config.merge_clean_gap_fwhm;

  p.roi_significance_z = config.roi_significance_z;
  p.observable_peak_initial_significance_threshold = config.observable_peak_initial_significance_threshold;
  p.observable_peak_final_significance_threshold = config.observable_peak_final_significance_threshold;

  p.step_cont_min_peak_significance = config.step_cont_min_peak_significance;
  p.step_trial_chi2_margin = config.step_trial_chi2_margin;

  p.initial_manual_rel_eff_max_chi2_dof = config.initial_manual_rel_eff_max_chi2_dof;

  return p;
}//settings_to_genes


void init_genes( NuclideConfigSolution &p, const std::function<double(void)> &rnd01 )
{
  // FWHM functional form [0, 3]
  p.fwhm_functional_form = static_cast<int>( 4 * rnd01() );
  if( p.fwhm_functional_form > 3 ) p.fwhm_functional_form = 3;

  // Manual RelEff parameters
  p.rel_eff_manual_base_rel_eff_uncert = 0.0 + 0.5 * rnd01();
  p.initial_nuc_match_cluster_num_sigma = 0.5 + 3.5 * rnd01();
  p.manual_eff_cluster_num_sigma = 0.5 + 3.5 * rnd01();

  // Manual RelEff equation forms/orders
  // Manual rel-eff form/order are selected per spectrum by AICc; only the penalty scales are genes.
  p.manual_releff_aicc_penalty = 0.5 + 7.5 * rnd01();
  p.cont_order_aicc_penalty = 0.5 + 7.5 * rnd01();

  // Manual clustering thresholds
  p.manual_keep_significance_z = 0.5 + 5.5 * rnd01();
  p.manual_rel_eff_sol_min_fwhm_roi = 0.5 + 2.5 * rnd01();
  p.manual_rel_eff_sol_max_fwhm = 5.0 + 25.0 * rnd01();
  p.merge_tail_z = 0.5 + 4.5 * rnd01();
  p.merge_clean_gap_fwhm = 0.25 + 2.75 * rnd01();
  p.manual_roi_core_num_fwhm = 0.75 + 1.75 * rnd01();

  // Auto RelEff parameters
  // FwhmForm: limited to Berstein_2(8) through Berstein_5(11)
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
  static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
  p.fwhm_form = 8 + static_cast<int>( 4 * rnd01() );
  if( p.fwhm_form > 11 ) p.fwhm_form = 11;

  p.rel_eff_auto_base_rel_eff_uncert = 0.0 + 0.5 * rnd01();
  p.auto_rel_eff_cluster_num_sigma = 0.5 + 4.5 * rnd01();
  p.auto_keep_significance_z = 0.5 + 7.5 * rnd01();
  p.auto_roi_core_num_fwhm = 0.75 + 1.75 * rnd01();
  p.roi_extend_z = 0.75 + 3.25 * rnd01();
  p.roi_max_num_fwhm = 2.0 + 6.0 * rnd01();
  p.auto_rel_eff_sol_max_fwhm = 4.0 + 21.0 * rnd01();
  p.auto_rel_eff_sol_min_fwhm_roi = 0.5 + 2.5 * rnd01();

  // RelActAuto model
  // RelEffEqnForm: 0..4 (LnX, LnY, LnXLnY, FramEmpirical, FramPhysicalModel)
  p.rel_eff_eqn_type = static_cast<int>( 5 * rnd01() );
  if( p.rel_eff_eqn_type > 4 ) p.rel_eff_eqn_type = 4;
  p.rel_eff_eqn_order = static_cast<int>( 7 * rnd01() );
  if( p.rel_eff_eqn_order > 6 ) p.rel_eff_eqn_order = 6;

  // Desperation physical model
  p.desperation_phys_model_atomic_number = 6.0 + 76.0 * rnd01();
  p.desperation_phys_model_areal_density_g_per_cm2 = 0.1 + 19.9 * rnd01();

  // Booleans
  p.nucs_of_el_same_age = rnd01() < 0.5 ? 0 : 1;
  p.phys_model_use_hoerl = rnd01() < 0.5 ? 0 : 1;
  p.fit_energy_cal = rnd01() < 0.5 ? 0 : 1;

  // ROI significance and observable peak thresholds
  p.roi_significance_z = 1.5 + 6.5 * rnd01();
  p.observable_peak_initial_significance_threshold = 1.0 + 4.0 * rnd01();
  p.observable_peak_final_significance_threshold = 0.5 + 4.5 * rnd01();

  // Step continuum
  p.step_cont_min_peak_significance = 5.0 + 95.0 * rnd01();
  p.step_trial_chi2_margin = 15.0 * rnd01();

  // Manual rel-eff chi2/dof cap.  Only randomize when GAOptimized; otherwise
  // genes_to_settings will override this with the CLI-resolved value.
  // Log-uniform in [5, 1000] gives the GA the freedom to effectively disable
  // the cut by picking a large value.
  if( sm_rel_eff_chi2_cap_mode == RelEffChi2CapMode::GAOptimized )
  {
    const double log_min = std::log( 5.0 );
    const double log_max = std::log( 1000.0 );
    p.initial_manual_rel_eff_max_chi2_dof = std::exp( log_min + (log_max - log_min) * rnd01() );
  }
  else
  {
    p.initial_manual_rel_eff_max_chi2_dof = sm_rel_eff_chi2_cap_fixed_value;
  }

  // skew_type is not GA-optimized
}//init_genes


bool eval_solution( const NuclideConfigSolution &p, NuclideConfigCost &c )
{
  const PeakFitForNuclideConfig config = genes_to_settings( p );

  double fg = 0.0, bg = 0.0;
  c.objective1 = ns_ga_eval_fcn( config, &fg, &bg );
  c.objective_fg = fg;
  c.objective_bg = bg;

  ns_num_evals_this_generation += 1;
  if( PeakFitImprove::debug_printout && ((ns_num_evals_this_generation % 10) == 0) )
    cout << "Have evaluated " << ns_num_evals_this_generation.load() << " individuals this generation." << endl;

  if( std::isnan( c.objective1 ) || std::isinf( c.objective1 ) )
  {
    cerr << "Got an objective of " << c.objective1 << " for " << p.to_string( ", " ) << endl;
    return false;
  }

  return true;
}//eval_solution


// Helper macro to reduce repetition in mutate for double genes
// Each gene: if random < threshold, apply mutation, check range
#define MUTATE_DOUBLE( field, lo, hi ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field += mu * (rnd01() - rnd01()); \
    in_range = in_range && (X_new.field >= (lo) && X_new.field < (hi)); \
  }

// For integer genes: randomly reassign with small probability
#define MUTATE_INT( field, lo, hi ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field = static_cast<int>( (lo) + ((hi) - (lo) + 1) * rnd01() ); \
    if( X_new.field > (hi) ) X_new.field = (hi); \
    in_range = in_range && (X_new.field >= (lo) && X_new.field <= (hi)); \
  }

// For boolean genes: toggle with small probability
#define MUTATE_BOOL( field ) \
  if( rnd01() < mutate_threshold ) { \
    X_new.field = (X_new.field == 0) ? 1 : 0; \
  }


NuclideConfigSolution mutate( const NuclideConfigSolution &X_base,
                              const std::function<double(void)> &rnd01,
                              double shrink_scale )
{
  NuclideConfigSolution X_new;
  const double mu = 0.2 * shrink_scale;
  const double mutate_threshold = PeakFitImprove::sm_ga_mutate_threshold;
  bool in_range;

  size_t num_tries = 0;

  do
  {
    num_tries += 1;
    if( num_tries > 1000 )
    {
      cerr << "NuclideConfig_GA::mutate: took over " << num_tries << " tries - reinitializing." << endl;
      std::random_device rng;
      std::uniform_real_distribution<double> unif_dist( 0.0, 1.0 );
      init_genes( X_new, [&](){ return unif_dist( rng ); } );
      return X_new;
    }

    in_range = true;
    X_new = X_base;

    // FWHM functional form
    MUTATE_INT( fwhm_functional_form, 0, 3 )

    // Manual RelEff doubles
    MUTATE_DOUBLE( rel_eff_manual_base_rel_eff_uncert, 0.0, 0.5 )
    MUTATE_DOUBLE( initial_nuc_match_cluster_num_sigma, 0.5, 4.0 )
    MUTATE_DOUBLE( manual_eff_cluster_num_sigma, 0.5, 4.0 )

    // AICc complexity-penalty scales (rel-eff form/order ladder + continuum-order selection)
    MUTATE_DOUBLE( manual_releff_aicc_penalty, 0.5, 8.0 )
    MUTATE_DOUBLE( cont_order_aicc_penalty, 0.5, 8.0 )

    // Manual clustering doubles
    MUTATE_DOUBLE( manual_keep_significance_z, 0.5, 6.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_min_fwhm_roi, 0.5, 3.0 )
    MUTATE_DOUBLE( manual_rel_eff_sol_max_fwhm, 5.0, 30.0 )
    MUTATE_DOUBLE( merge_tail_z, 0.5, 5.0 )
    MUTATE_DOUBLE( merge_clean_gap_fwhm, 0.25, 3.0 )
    MUTATE_DOUBLE( manual_roi_core_num_fwhm, 0.75, 2.5 )

    // Auto RelEff - fwhm_form limited to Berstein_2 through Berstein_5
    static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_2) == 8, "Berstein_2 enum value changed" );
    static_assert( static_cast<int>(RelActCalcAuto::FwhmForm::Berstein_5) == 11, "Berstein_5 enum value changed" );
    MUTATE_INT( fwhm_form, 8, 11 )
    MUTATE_DOUBLE( rel_eff_auto_base_rel_eff_uncert, 0.0, 0.5 )
    MUTATE_DOUBLE( auto_rel_eff_cluster_num_sigma, 0.5, 5.0 )
    MUTATE_DOUBLE( auto_keep_significance_z, 0.5, 8.0 )
    MUTATE_DOUBLE( auto_roi_core_num_fwhm, 0.75, 2.5 )
    MUTATE_DOUBLE( roi_extend_z, 0.75, 4.0 )
    MUTATE_DOUBLE( roi_max_num_fwhm, 2.0, 8.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_max_fwhm, 4.0, 25.0 )
    MUTATE_DOUBLE( auto_rel_eff_sol_min_fwhm_roi, 0.5, 3.0 )

    // RelActAuto model
    MUTATE_INT( rel_eff_eqn_type, 0, 4 )
    MUTATE_INT( rel_eff_eqn_order, 0, 6 )

    // Desperation physical model
    MUTATE_DOUBLE( desperation_phys_model_atomic_number, 6.0, 82.0 )
    MUTATE_DOUBLE( desperation_phys_model_areal_density_g_per_cm2, 0.1, 20.0 )

    // Booleans
    MUTATE_BOOL( nucs_of_el_same_age )
    MUTATE_BOOL( phys_model_use_hoerl )
    MUTATE_BOOL( fit_energy_cal )

    // ROI significance
    MUTATE_DOUBLE( roi_significance_z, 1.5, 8.0 )
    MUTATE_DOUBLE( observable_peak_initial_significance_threshold, 1.0, 5.0 )
    MUTATE_DOUBLE( observable_peak_final_significance_threshold, 0.5, 5.0 )

    // Step continuum
    MUTATE_DOUBLE( step_cont_min_peak_significance, 5.0, 100.0 )
    MUTATE_DOUBLE( step_trial_chi2_margin, 0.0, 15.0 )

    // Manual rel-eff chi2/dof cap - only meaningful when GAOptimized.  When
    // Fixed/Disabled, the gene value is ignored by resolve_chi2_dof_cap().
    if( sm_rel_eff_chi2_cap_mode == RelEffChi2CapMode::GAOptimized )
    {
      MUTATE_DOUBLE( initial_manual_rel_eff_max_chi2_dof, 5.0, 1000.0 )
    }

  } while( !in_range );

  return X_new;
}//mutate

#undef MUTATE_DOUBLE
#undef MUTATE_INT
#undef MUTATE_BOOL


// Helper for crossover of double genes
#define CROSSOVER_DOUBLE( field ) \
  if( rnd01() < crossover_threshold ) { \
    const double r = rnd01(); \
    X_new.field = r * X1.field + (1.0 - r) * X2.field; \
  } else if( rnd01() < 0.5 ) { \
    X_new.field = X2.field; \
  }

// For int genes: pick from one parent
#define CROSSOVER_INT( field ) \
  if( rnd01() < 0.5 ) { \
    X_new.field = X2.field; \
  }


NuclideConfigSolution crossover( const NuclideConfigSolution &X1,
                                 const NuclideConfigSolution &X2,
                                 const std::function<double(void)> &rnd01 )
{
  const double crossover_threshold = PeakFitImprove::sm_ga_crossover_threshold;

  NuclideConfigSolution X_new = X1;

  // FWHM functional form
  CROSSOVER_INT( fwhm_functional_form )

  // Manual RelEff doubles
  CROSSOVER_DOUBLE( rel_eff_manual_base_rel_eff_uncert )
  CROSSOVER_DOUBLE( initial_nuc_match_cluster_num_sigma )
  CROSSOVER_DOUBLE( manual_eff_cluster_num_sigma )

  // Manual equation forms/orders
  CROSSOVER_DOUBLE( manual_releff_aicc_penalty )
  CROSSOVER_DOUBLE( cont_order_aicc_penalty )

  // Manual clustering doubles
  CROSSOVER_DOUBLE( manual_keep_significance_z )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_min_fwhm_roi )
  CROSSOVER_DOUBLE( manual_rel_eff_sol_max_fwhm )
  CROSSOVER_DOUBLE( merge_tail_z )
  CROSSOVER_DOUBLE( merge_clean_gap_fwhm )
  CROSSOVER_DOUBLE( manual_roi_core_num_fwhm )

  // Auto RelEff
  CROSSOVER_INT( fwhm_form )
  CROSSOVER_DOUBLE( rel_eff_auto_base_rel_eff_uncert )
  CROSSOVER_DOUBLE( auto_rel_eff_cluster_num_sigma )
  CROSSOVER_DOUBLE( auto_keep_significance_z )
  CROSSOVER_DOUBLE( auto_roi_core_num_fwhm )
  CROSSOVER_DOUBLE( roi_extend_z )
  CROSSOVER_DOUBLE( roi_max_num_fwhm )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_max_fwhm )
  CROSSOVER_DOUBLE( auto_rel_eff_sol_min_fwhm_roi )

  // RelActAuto model
  CROSSOVER_INT( rel_eff_eqn_type )
  CROSSOVER_INT( rel_eff_eqn_order )

  // Desperation physical model
  CROSSOVER_DOUBLE( desperation_phys_model_atomic_number )
  CROSSOVER_DOUBLE( desperation_phys_model_areal_density_g_per_cm2 )

  // Booleans
  CROSSOVER_INT( nucs_of_el_same_age )
  CROSSOVER_INT( phys_model_use_hoerl )
  CROSSOVER_INT( fit_energy_cal )

  // ROI significance
  CROSSOVER_DOUBLE( roi_significance_z )
  CROSSOVER_DOUBLE( observable_peak_initial_significance_threshold )
  CROSSOVER_DOUBLE( observable_peak_final_significance_threshold )

  // Step continuum
  CROSSOVER_DOUBLE( step_cont_min_peak_significance )
  CROSSOVER_DOUBLE( step_trial_chi2_margin )

  // Manual rel-eff chi2/dof cap - only crossed when GAOptimized.
  if( sm_rel_eff_chi2_cap_mode == RelEffChi2CapMode::GAOptimized )
  {
    CROSSOVER_DOUBLE( initial_manual_rel_eff_max_chi2_dof )
  }

  return X_new;
}//crossover

#undef CROSSOVER_DOUBLE
#undef CROSSOVER_INT


double calculate_SO_total_fitness( const GA_Type::thisChromosomeType &X )
{
  double final_cost = 0.0;
  final_cost += X.middle_costs.objective1;
  return final_cost;
}//calculate_SO_total_fitness


void SO_report_generation( int generation_number,
                           const EA::GenerationType<NuclideConfigSolution,NuclideConfigCost> &last_generation,
                           const NuclideConfigSolution &best_genes )
{
  bool best_yet = false;
  if( !sm_set_best_genes || (last_generation.best_total_cost < sm_best_total_cost) )
  {
    best_yet = true;
    sm_set_best_genes = true;
    sm_best_genes = best_genes;
    sm_best_total_cost = last_generation.best_total_cost;
  }

  // When the bg trial is enabled, report the best individual's fg/bg breakdown.  The components are
  // already cached in the best chromosome's middle_costs (set by eval_solution), and openGA derives
  // best_genes from chromosomes[best_chromosome_index], so we read them directly instead of
  // re-evaluating: this reconciles exactly with best_total_cost and avoids a (potentially
  // non-deterministic) extra full-dataset eval per new-best generation.
  if( sm_do_background_fit_trial && best_yet )
  {
    const int best_idx = last_generation.best_chromosome_index;
    if( (best_idx >= 0) && (best_idx < static_cast<int>(last_generation.chromosomes.size())) )
    {
      sm_best_fg_cost = last_generation.chromosomes[best_idx].middle_costs.objective_fg;
      sm_best_bg_cost = last_generation.chromosomes[best_idx].middle_costs.objective_bg;
    }
  }
  const double best_fg = sm_do_background_fit_trial ? sm_best_fg_cost : last_generation.best_total_cost;
  const double best_bg = sm_do_background_fit_trial ? sm_best_bg_cost : 0.0;

  if( sm_do_background_fit_trial )
  {
    sm_output_file
      << generation_number << "\t"
      << last_generation.average_cost << "\t"
      << last_generation.best_total_cost << "\t"
      << best_fg << "\t"
      << best_bg << "\t"
      << "{" << best_genes.to_string( ", " ) << "}"
      << "\t" << (best_yet ? "BestYet" : "NoImprovement")
      << "\n\n";
  }
  else
  {
    sm_output_file
      << generation_number << "\t"
      << last_generation.average_cost << "\t"
      << last_generation.best_total_cost << "\t"
      << "{" << best_genes.to_string( ", " ) << "}"
      << "\t" << (best_yet ? "BestYet" : "NoImprovement")
      << "\n\n";
  }
  sm_output_file.flush();  // persist each generation so an interrupted run keeps its full history

  cout
    << "Generation [" << generation_number << "], "
    << "Best=" << last_generation.best_total_cost;
  if( sm_do_background_fit_trial )
  {
    const double weighted_bg = sm_background_fit_penalty_weight * best_bg;
    cout << " (fg=" << best_fg
         << ", bg_raw=" << best_bg
         << ", weighted_bg=" << weighted_bg
         << ", combined=" << (best_fg + weighted_bg) << ")";
  }
  cout << ", Average=" << last_generation.average_cost << ", "
    << (best_yet ? "Best generation yet" : "no improvement")
    << "\n  Best genes: {\n\t" << best_genes.to_string( "\n\t" ) << "\n}\n"
    << "Exe_time=" << last_generation.exe_time
    << endl;

  // Write HTML output for the best solution of this generation
  if( best_yet && sm_precomputed_ptr && !sm_precomputed_ptr->empty() )
  {
    try
    {
      const PeakFitForNuclideConfig config = genes_to_settings( sm_best_genes );
      // Gen-numbered filename so each generation's best-genes HTML is preserved (a visual history of
      //  how the best fit improves), rather than overwriting a single file.  n42 output is skipped for
      //  these per-generation writes (empty dir) to avoid the disk churn; the final run writes n42.
      const string html_filename = PeakFitImprove::sm_output_file_prefix
        + "nuclide_config_ga_best_gen" + std::to_string( generation_number ) + ".html";
      write_results_html_and_n42( *sm_precomputed_ptr, config, sm_best_genes, sm_background_mode, html_filename, /*n42_dir=*/"" );
      cout << "Wrote HTML results to " << html_filename << endl;
    }
    catch( const std::exception &e )
    {
      cerr << "Warning: failed to write HTML results for generation " << generation_number
           << ": " << e.what() << endl;
    }
  }//if( best_yet && sm_precomputed_ptr )

  // Checkpoint the full population's genes this generation (opt-in via --checkpoint) so an
  //  interrupted run can be continued with --resume.  Written atomically (temp + rename).  Genes
  //  only - the RNG/stall-counter are not persisted; a resume continues from the same population.
  if( !sm_checkpoint_name.empty() )
  {
    const std::string path = PeakFitImprove::sm_output_file_prefix + sm_checkpoint_name
                             + "_nuclide_config_ga_checkpoint.tsv";
    const std::string tmp_path = path + ".tmp";
    try
    {
      std::ofstream ckpt( tmp_path.c_str(), std::ios::out | std::ios::trunc );
      if( !ckpt )
        throw std::runtime_error( "could not open '" + tmp_path + "' for writing" );

      // Header ('#'-prefixed, skipped on load): the run-options summary lets --resume refuse an
      //  incompatible continuation; the rest is human-readable progress metadata.
      ckpt << "# options: " << sm_checkpoint_options_summary << "\n"
           << "# generation=" << generation_number
           << " population=" << last_generation.chromosomes.size()
           << " best_total_cost=" << last_generation.best_total_cost << "\n";
      for( const auto &chromo : last_generation.chromosomes )
        ckpt << chromo.genes.to_string( "\t" ) << "\n";

      ckpt.close();
      if( !ckpt )
        throw std::runtime_error( "error while writing '" + tmp_path + "'" );

      if( std::rename( tmp_path.c_str(), path.c_str() ) != 0 )
        throw std::runtime_error( "could not rename '" + tmp_path + "' to '" + path + "'" );
    }catch( const std::exception &e )
    {
      cerr << "Warning: failed to write checkpoint for generation " << generation_number
           << ": " << e.what() << endl;
    }
  }//if( !sm_checkpoint_name.empty() )

  ns_num_evals_this_generation = 0;
}//SO_report_generation


namespace
{
struct RoiReview
{
  const PeakContinuum *continuum = nullptr;
  std::vector<PeakDef> peaks;
  double lower = 0.0;
  double upper = 0.0;
  double representative_fwhm = 0.0;
  double width_fwhm = 0.0;
  size_t channels = 0;
  double model_pearson_rms = 0.0;
  double snip_continuum_rms = 0.0;
};


std::vector<RoiReview> review_rois( const PrecomputedNuclideData &pd,
                                    const SpectrumEvaluation &evaluation )
{
  std::vector<RoiReview> answer;
  if( !evaluation.has_fit_result || !pd.foreground )
    return answer;

  const PeakFitResult &result = evaluation.fit_result;
  const std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> groups
    = group_peaks_by_roi( result.fit_peaks );

  FitPeaksForNuclides::detail::GlobalContinuumEstimate snip;
  const std::shared_ptr<const DetectorPeakResponse> drf = result.solution.m_drf;
  if( drf )
  {
    const std::function<double(double)> fwhm = [drf]( const double energy ) -> double {
      return drf->peakResolutionFWHM( static_cast<float>(energy) );
    };
    snip = FitPeaksForNuclides::detail::make_global_continuum(
      pd.foreground, fwhm, pd.det_type, pd.foreground->gamma_energy_min(),
      pd.foreground->gamma_energy_max() );
  }

  for( const std::pair<const PeakContinuum *, std::vector<PeakDef>> &group : groups )
  {
    if( !group.first )
      continue;
    RoiReview roi;
    roi.continuum = group.first;
    roi.peaks = group.second;
    roi.lower = group.first->lowerEnergy();
    roi.upper = group.first->upperEnergy();
    roi.channels = ReportDetail::roi_channel_count( *pd.foreground, roi.lower, roi.upper );

    std::vector<double> fwhms;
    for( const PeakDef &peak : roi.peaks )
      if( peak.fwhm() > 0.0 )
        fwhms.push_back( peak.fwhm() );
    if( !fwhms.empty() )
    {
      std::sort( fwhms.begin(), fwhms.end() );
      roi.representative_fwhm = fwhms[fwhms.size()/2];
    }
    roi.width_fwhm = ReportDetail::roi_width_in_fwhm(
      roi.lower, roi.upper, roi.representative_fwhm );

    const size_t first = pd.foreground->find_gamma_channel( static_cast<float>(roi.lower) );
    const size_t last = pd.foreground->find_gamma_channel( static_cast<float>(roi.upper) );
    std::vector<const PeakDef *> peak_ptrs;
    for( const PeakDef &peak : roi.peaks )
      peak_ptrs.push_back( &peak );
    double model_sq = 0.0;
    double snip_sq = 0.0;
    size_t used = 0;
    for( size_t channel = first;
         channel <= last && channel < pd.foreground->num_gamma_channels(); ++channel )
    {
      const double lo = pd.foreground->gamma_channel_lower( channel );
      const double hi = pd.foreground->gamma_channel_upper( channel );
      const double data = pd.foreground->gamma_channel_content( channel );
      const double cont = group.first->offset_integral(
        lo, hi, pd.foreground, peak_ptrs.data(), peak_ptrs.size() );
      double model = cont;
      for( const PeakDef &peak : roi.peaks )
        model += peak.gauss_integral( lo, hi );
      const double variance = std::max( 1.0, model );
      model_sq += (data - model)*(data - model) / variance;
      if( snip.valid() && snip.snip )
      {
        const double snip_count = snip.snip->gamma_integral(
          static_cast<float>(lo), static_cast<float>(hi) );
        snip_sq += (cont - snip_count)*(cont - snip_count)
            / std::max( 1.0, data );
      }
      used += 1;
    }
    if( used )
    {
      roi.model_pearson_rms = std::sqrt( model_sq / used );
      if( snip.valid() )
        roi.snip_continuum_rms = std::sqrt( snip_sq / used );
    }
    answer.push_back( std::move(roi) );
  }
  return answer;
}


bool objective_truth_match( const PeakDef &peak, const ExpectedPhotopeakInfo &truth )
{
  const double mean = peak.mean();
  return ((mean > truth.roi_lower) && (mean < truth.roi_upper))
      || (((mean + peak.sigma()) > truth.roi_lower)
          && ((mean - peak.sigma()) < truth.roi_upper));
}


std::string background_mode_name( const BackgroundMode mode )
{
  switch( mode )
  {
    case BackgroundMode::BackgroundSubtracted: return "Supplied long background subtraction";
    case BackgroundMode::NoBackground: return "Raw foreground; no supplied background";
    case BackgroundMode::NoBackgroundFitNorm: return "Raw foreground with explicit NORM co-fit";
  }
  return "Unknown";
}


std::string resolution_class_name( const PeakFitUtils::CoarseResolutionType type )
{
  using PeakFitUtils::CoarseResolutionType;
  switch( type )
  {
    case CoarseResolutionType::High: return "HPGe/High";
    case CoarseResolutionType::Low: return "NaI/Low";
    case CoarseResolutionType::LaBr: return "LaBr";
    case CoarseResolutionType::CZT: return "CZT";
    case CoarseResolutionType::MedRes: return "Medium resolution";
    case CoarseResolutionType::LowOrMedRes: return "Low or medium resolution";
    case CoarseResolutionType::Unknown: return "Unknown";
  }
  return "Unknown";
}


std::string utc_timestamp()
{
  const std::time_t now = std::time( nullptr );
  const std::tm * const utc = std::gmtime( &now );
  std::ostringstream out;
  if( utc )
    out << std::put_time( utc, "%Y-%m-%dT%H:%M:%SZ" );
  return out.str();
}


void append_reference_line( ReferenceLineInfo &lines, const double energy,
                            const Wt::WColor &color, const std::string &label )
{
  ReferenceLineInfo::RefLine line;
  line.m_energy = energy;
  line.m_normalized_intensity = 1.0;
  line.m_drf_factor = 1.0;
  line.m_shield_atten = 1.0;
  line.m_particle_sf_applied = 1.0;
  line.m_color = color;
  line.m_decay_intensity = 1.0;
  line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
  line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
  line.m_attenuation_applies = false;
  line.m_decaystr = label;
  lines.m_ref_lines.push_back( line );
}


void write_cached_manual_review_report(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const NuclideConfigSolution &genes,
  const BackgroundMode bg_mode,
  const std::string &html_filename,
  const std::string &n42_output_dir,
  const ConfigEvaluation &evaluation,
  const std::string &run_metadata )
{
  ReportDetail::validate_evaluation_coverage( precomputed.size(), evaluation );
  for( size_t i = 0; i < precomputed.size(); ++i )
  {
    if( evaluation.spectra[i].spectrum_id
          != ReportDetail::canonical_spectrum_key(precomputed[i])
        || evaluation.spectra[i].anchor_id
          != ReportDetail::stable_spectrum_id(precomputed[i]) )
      throw std::runtime_error( "Evaluation record is not associated with its selected spectrum" );
  }

  std::ofstream html( html_filename );
  if( !html.good() )
    throw std::runtime_error( "Could not open HTML output: " + html_filename );
  std::ostringstream page_header;
  D3SpectrumExport::write_html_page_header(
    page_header, "FitPeaksForNuclides manual-review gallery", "InterSpec_resources" );
  std::string page_header_text = page_header.str();
  const size_t head_end = page_header_text.find( "<head>" );
  if( head_end == std::string::npos )
    throw std::runtime_error( "D3 HTML page header did not contain a head element" );
  page_header_text.insert( head_end + 6, "\r\n<meta charset=\"UTF-8\">" );
  html << page_header_text;

  struct Summary
  {
    size_t index = 0;
    std::string id;
    double cost = 0.0;
    double largest_area_error = 0.0;
    size_t missed = 0;
    size_t spurious = 0;
    double widest_fwhm = 0.0;
    size_t widest_channels = 0;
    size_t max_peaks_roi = 0;
    size_t stepped = 0;
    double continuum_mismatch = 0.0;
    std::string state;
    std::vector<RoiReview> rois;
  };
  std::vector<Summary> summaries( precomputed.size() );
  std::set<std::string> detector_models;

  for( size_t i = 0; i < precomputed.size(); ++i )
  {
    const PrecomputedNuclideData &pd = precomputed[i];
    const SpectrumEvaluation &ev = evaluation.spectra[i];
    detector_models.insert( pd.src_info->detector_name );
    Summary &summary = summaries[i];
    summary.index = i;
    summary.id = ev.anchor_id;
    summary.cost = ev.accuracy.cost
        + sm_background_fit_penalty_weight*ev.background_penalty;
    summary.missed = ev.accuracy.missed_definitely_wanted;
    summary.spurious = ev.accuracy.extra_peaks;
    summary.state = ev.mechanical_failure ? "failure"
      : (ev.legitimate_empty ? "empty" : "success");
    summary.rois = review_rois( pd, ev );
    for( const RoiReview &roi : summary.rois )
    {
      summary.widest_fwhm = std::max( summary.widest_fwhm, roi.width_fwhm );
      summary.widest_channels = std::max( summary.widest_channels, roi.channels );
      summary.max_peaks_roi = std::max( summary.max_peaks_roi, roi.peaks.size() );
      summary.continuum_mismatch = std::max(
        summary.continuum_mismatch, roi.snip_continuum_rms );
      if( roi.continuum && PeakContinuum::is_step_continuum(roi.continuum->type()) )
        summary.stepped += 1;
    }
    if( ev.has_fit_result )
    {
      const std::vector<ExpectedPhotopeakInfo> truth
        = PeakFitImproveData::filter_photopeaks_for_scoring(
            pd.src_info->expected_signal_photopeaks, pd.det_type );
      for( const ExpectedPhotopeakInfo &expected : truth )
      {
        const PeakDef *best = nullptr;
        for( const PeakDef &peak : ev.fit_result.observable_peaks )
        {
          if( objective_truth_match(peak, expected)
              && (!best || std::fabs(peak.mean() - expected.effective_energy)
                    < std::fabs(best->mean() - expected.effective_energy)) )
            best = &peak;
        }
        if( best )
        {
          constexpr double f_rel = 0.05;
          const double area = std::max( 1.0, expected.peak_area );
          const double sigma = std::sqrt( area + (f_rel*area)*(f_rel*area) );
          summary.largest_area_error = std::max( summary.largest_area_error,
            std::fabs(best->amplitude() - expected.peak_area) / sigma );
        }
      }
    }
  }

  html << "<body><style>"
       << "body{font-family:-apple-system,BlinkMacSystemFont,Arial,sans-serif;margin:0;background:#f6f7f9;color:#17202a;}"
       << "header,.report-block{max-width:1500px;margin:14px auto;padding:14px 18px;background:white;border:1px solid #ccd3da;border-radius:7px;}"
       << "table{border-collapse:collapse;width:100%;font-size:13px;}th,td{border:1px solid #c8ced5;padding:5px 7px;text-align:right;}th{background:#e9edf2;position:sticky;top:0;}td:first-child,th:first-child{text-align:left;}"
       << ".summary-table-wrap{max-height:65vh;overflow:auto}.failure{background:#ffe8e8}.empty{background:#fff7d6}.miss{color:#a40000;font-weight:bold}.spurious{color:#8a3d00;font-weight:bold}"
       << ".roi-review-row{transition:background-color .15s ease,box-shadow .15s ease}.roi-review-row.roi-selected{background:#dcecff;box-shadow:inset 4px 0 #2468b4;font-weight:600}"
       << ".spectrum{max-width:1500px;margin:12px auto;background:white;border:1px solid #bcc5cf;border-radius:7px;padding:0 12px;}"
       << ".spectrum>summary{cursor:pointer;padding:12px 4px;font-weight:600}.chart{min-height:300px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:10px}.meta{line-height:1.45}"
       << ".notice{border-left:5px solid #d17b00;padding:8px 12px;background:#fff7e6}.error{white-space:pre-wrap;color:#8b0000;background:#fff0f0;padding:8px}.warnings{white-space:pre-wrap;background:#fff8dd;padding:8px}"
       << ".noteable{cursor:pointer}.noteable:hover{outline:2px solid #6b8fd6}.note-button{margin:4px;padding:5px 9px}.mono{font-family:Menlo,monospace;word-break:break-word}"
       << ".modal{position:fixed;z-index:10000;top:0;right:0;bottom:0;left:0;background:rgba(0,0,0,.45);display:flex;align-items:center;justify-content:center}.modal[hidden]{display:none}.modal-card{background:white;border:1px solid #667;border-radius:6px;padding:16px;width:min(650px,90vw);max-height:85vh;overflow:auto;box-shadow:0 8px 30px rgba(0,0,0,.35)}#remark-text,#remark-export-text{box-sizing:border-box;width:100%}#remark-text{height:130px}.controls{display:flex;flex-wrap:wrap;align-items:center}.controls>*{margin:3px 7px 3px 0}"
       << "</style>\n";

  html << "<header><h1>FitPeaksForNuclides manual-review gallery</h1>"
       << "<p class=\"notice\"><b>Production vs diagnostics:</b> Plots and tables use the exact cached result returned during command-line evaluation; non-Success partial objects are explicitly labeled diagnostic rather than production. "
       << "R4 shadow output is not rendered as a fit. The public fitted-total overlay is distinct from objective-observable peaks; both are marked separately. Private solver/nuisance counts are disclosed in the representation audit but are not plotted because they use the solver's fitted-calibration frame. SNIP continuum mismatch is a clearly labeled report-only diagnostic, not an objective term.</p>"
       << "<div class=\"grid meta\"><div><b>Generated:</b> " << ReportDetail::html_escape(utc_timestamp())
       << "<br><b>Background mode:</b> " << ReportDetail::html_escape(background_mode_name(bg_mode))
       << "<br><b>Automatic interferer fitting:</b> "
       << (sm_disable_auto_interferer_fit ? "DISABLED (source-only baseline)" : "enabled")
       << "<br><b>Selected spectra:</b> " << evaluation.spectra.size()
       << "</div><div><b>Successes:</b> " << evaluation.successes
       << "<br><b>Legitimate empties:</b> " << evaluation.legitimate_empties
       << "<br><b>Mechanical failures:</b> " << evaluation.mechanical_failures
       << "<br><b>Objective:</b> " << evaluation.total_cost
       << " (foreground " << evaluation.total_fg << ", raw background "
       << evaluation.total_bg_raw << " × " << sm_background_fit_penalty_weight
       << " = " << (evaluation.total_bg_raw*sm_background_fit_penalty_weight)
       << ")</div></div>"
       << "<p><b>Detector models represented:</b> ";
  for( std::set<std::string>::const_iterator iter = detector_models.begin();
       iter != detector_models.end(); ++iter )
  {
    if( iter != detector_models.begin() )
      html << ", ";
    html << ReportDetail::html_escape( *iter );
  }
  html << "</p>"
       << "<p><b>Exact filters / reproduction metadata:</b> <span class=\"mono\">"
       << ReportDetail::html_escape(run_metadata) << "</span></p>"
       << "<details><summary>Configuration / gene values</summary><table><tr><th>Gene</th><th>Value</th></tr>";
  std::vector<std::string> gene_lines;
  SpecUtils::split( gene_lines, genes.to_string("\n"), "\n" );
  for( const std::string &line : gene_lines )
  {
    const size_t eq = line.find( '=' );
    if( eq != std::string::npos )
      html << "<tr><td class=\"mono\">" << ReportDetail::html_escape(line.substr(0,eq))
           << "</td><td class=\"mono\">" << ReportDetail::html_escape(line.substr(eq+1))
           << "</td></tr>";
  }
  html << "</table></details><details><summary>Aggregate objective components</summary><table>"
       << "<tr><th>Area cost</th><th>Find reward</th><th>Candidate reward</th><th>Sum miss fractions</th><th>Missed definitely wanted</th><th>Spurious</th></tr><tr>"
       << "<td>" << evaluation.accuracy_totals.area_cost << "</td><td>"
       << evaluation.accuracy_totals.find_reward << "</td><td>"
       << evaluation.accuracy_totals.candidate_reward << "</td><td>"
       << evaluation.accuracy_totals.miss_fraction << "</td><td>"
       << evaluation.accuracy_totals.missed_definitely_wanted << "</td><td>"
       << evaluation.accuracy_totals.extra_peaks << "</td></tr></table></details></header>\n";

  html << "<section class=\"report-block\"><h2>Spectrum index</h2>"
       << "<p>Diagnostic rankings help find review targets; they are not necessarily optimization-objective terms. Click any row to attach a remark.</p>"
       << "<div class=\"controls\"><label>Filter <input id=\"summary-filter\" type=\"search\" placeholder=\"source, detector, status\"></label>"
       << "<label>State <select id=\"state-filter\"><option value=\"\">all</option><option>success</option><option>empty</option><option>failure</option></select></label>"
       << "<label>Sort <select id=\"summary-sort\"><option value=\"cost\">worst scalar cost</option><option value=\"area\">largest report-only area pull</option><option value=\"missed\">missed peaks</option><option value=\"spurious\">spurious peaks</option><option value=\"fwhm\">widest ROI (FWHM)</option><option value=\"channels\">widest ROI (channels)</option><option value=\"shared\">most peaks sharing ROI</option><option value=\"stepped\">stepped continua</option><option value=\"mismatch\">continuum mismatch</option><option value=\"state\">failures / empties</option></select></label>"
       << "<button id=\"export-remarks\" class=\"note-button\">Copy/export remarks</button></div>"
       << "<div class=\"summary-table-wrap\"><table id=\"summary-table\"><thead><tr><th>Spectrum</th><th>Detector</th><th>Location</th><th>Live time</th><th>Status</th><th>Cost</th><th>Report-only max |area pull|</th><th>Missed</th><th>Spurious</th><th>Max ROI FWHM</th><th>Max channels</th><th>Peaks/ROI</th><th>Steps</th><th>SNIP mismatch</th></tr></thead><tbody>";
  for( const Summary &summary : summaries )
  {
    const PrecomputedNuclideData &pd = precomputed[summary.index];
    const DataSrcInfo &info = *pd.src_info;
    const SpectrumEvaluation &ev = evaluation.spectra[summary.index];
    html << "<tr class=\"noteable " << summary.state << "\" data-state=\"" << summary.state
         << "\" data-cost=\"" << summary.cost << "\" data-area=\"" << summary.largest_area_error
         << "\" data-missed=\"" << summary.missed << "\" data-spurious=\"" << summary.spurious
         << "\" data-fwhm=\"" << summary.widest_fwhm << "\" data-channels=\"" << summary.widest_channels
         << "\" data-shared=\"" << summary.max_peaks_roi << "\" data-stepped=\"" << summary.stepped
         << "\" data-mismatch=\"" << summary.continuum_mismatch << "\" data-note-context=\""
         << ReportDetail::html_escape(ev.spectrum_id + " [" + summary.id + "] summary row")
         << "\"><td><a href=\"#"
         << summary.id << "\">" << ReportDetail::html_escape(info.src_info.src_name)
         << "</a><br><small class=\"mono\">" << ReportDetail::html_escape(summary.id)
         << "</small></td><td>" << ReportDetail::html_escape(info.detector_name)
         << "</td><td>" << ReportDetail::html_escape(info.location_name)
         << "</td><td>" << ReportDetail::html_escape(info.live_time_name)
         << "</td><td>" << summary.state << "</td><td>" << summary.cost
         << "</td><td>" << summary.largest_area_error << "</td><td>" << summary.missed
         << "</td><td>" << summary.spurious << "</td><td>" << summary.widest_fwhm
         << "</td><td>" << summary.widest_channels << "</td><td>" << summary.max_peaks_roi
         << "</td><td>" << summary.stepped << "</td><td>" << summary.continuum_mismatch
         << "</td></tr>\n";
  }
  html << "</tbody></table></div></section>\n";

  // Every selected spectrum gets a section, including exceptions and legitimate empties.
  for( const Summary &summary : summaries )
  {
    const size_t i = summary.index;
    const PrecomputedNuclideData &pd = precomputed[i];
    const DataSrcInfo &info = *pd.src_info;
    const InjectSourceInfo &src = info.src_info;
    const SpectrumEvaluation &ev = evaluation.spectra[i];
    const PeakFitResult * const fit = ev.has_fit_result ? &ev.fit_result : nullptr;
    const std::string note_context = ev.spectrum_id + " [" + summary.id + "]";
    const std::string details_id = "spectrum_details_" + std::to_string(i);
    const std::string chart_id = "spectrum_chart_" + std::to_string(i);
    const std::string init_fn = "init_spectrum_chart_" + std::to_string(i);

    html << "<details class=\"spectrum " << summary.state << "\" id=\"" << summary.id
         << "\" data-spectrum-id=\"" << summary.id << "\"><summary>"
         << ReportDetail::html_escape(src.src_name) << " — "
         << ReportDetail::html_escape(info.detector_name) << " — "
         << ReportDetail::html_escape(info.location_name) << " — "
         << ReportDetail::html_escape(info.live_time_name) << " — "
         << ReportDetail::html_escape(ev.status) << " — cost " << summary.cost
         << "</summary><div id=\"" << details_id << "\">"
         << "<button class=\"note-button add-note\" data-note-context=\""
         << ReportDetail::html_escape(note_context + " spectrum section") << "\">Add remark</button>"
         << "<div class=\"grid meta\"><div><b>Stable spectrum ID:</b> <span class=\"mono\">"
         << ReportDetail::html_escape(summary.id) << "</span><br><b>Corpus key:</b> <span class=\"mono\">"
         << ReportDetail::html_escape(ev.spectrum_id) << "</span><br><b>Source truth:</b> "
         << ReportDetail::html_escape(src.src_name) << "<br><b>Detector/class:</b> "
         << ReportDetail::html_escape(info.detector_name) << " / "
         << ReportDetail::html_escape(resolution_class_name(info.det_type))
         << "<br><b>Location:</b> " << ReportDetail::html_escape(info.location_name)
         << "</div><div><b>Fit status:</b> " << ReportDetail::html_escape(ev.status)
         << "<br><b>Foreground live time:</b> " << (pd.foreground ? pd.foreground->live_time() : 0.0)
         << " s<br><b>Background live time:</b> " << (pd.background ? pd.background->live_time() : 0.0)
         << " s<br><b>Background plot scale:</b> "
         << ((pd.foreground && pd.background && pd.background->live_time() > 0.0)
             ? (pd.foreground->live_time()/pd.background->live_time()) : 0.0)
         << "<br><b>Automatic interferer fit:</b> "
         << (sm_disable_auto_interferer_fit ? "disabled" : "enabled");
    const std::pair<double,double> spectrum_extent = pd.foreground
      ? FitPeaksForNuclides::find_valid_energy_range( pd.foreground )
      : std::make_pair( 0.0, 0.0 );
    html << "<br><b>Spectroscopic extent:</b> "
         << ((spectrum_extent.second > spectrum_extent.first)
               ? (std::to_string(spectrum_extent.first) + "–" + std::to_string(spectrum_extent.second) + " keV")
               : "unavailable") << "</div></div>";

    if( !ev.error_message.empty() )
      html << "<h3>Complete error</h3><div class=\"error noteable\" data-note-context=\""
           << ReportDetail::html_escape(note_context + " error") << "\">"
           << ReportDetail::html_escape(ev.error_message) << "</div>";
    if( fit && !fit->warnings.empty() )
    {
      html << "<h3>Warnings</h3><div class=\"warnings noteable\" data-note-context=\""
           << ReportDetail::html_escape(note_context + " warnings") << "\">";
      for( const std::string &warning : fit->warnings )
        html << ReportDetail::html_escape(warning) << "<br>";
      html << "</div>";
    }

    html << "<h3>Spectrum and " << (ev.mechanical_failure ? "partial failed result" : "production fit")
         << "</h3><p>Black: measured foreground. Steel blue: supplied background scaled by live time. The derived background-subtracted foreground is intentionally omitted because the production fit peaks and continua are defined against the original foreground. "
         << (ev.mechanical_failure
              ? "The peak overlay is a partial object returned with the failed status; it is diagnostic and is not labeled as a production result. "
              : "The peak overlay is the combined public production result. ")
         << "It uses <code>PeakFitResult::fit_peaks</code> and draws that representation's fitted total model and per-ROI fitted continuum. Magenta mean markers show the separately refitted <code>observable_peaks</code> used by the objective; cyan markers show combined public fitted means. Private solver/nuisance means are intentionally not overlaid because they are expressed in the solver's fitted-calibration frame rather than the displayed foreground calibration. "
         << "Truth markers are green when matched, red only for objective-counted definite misses, amber for other scored-but-unmatched truth, and gray when excluded by the scoring filter. Click a fitted peak shape to highlight its ROI in the diagnostics table.</p>"
         << "<div id=\"" << chart_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>";

    // Chart setup is lazy: hundreds of embedded spectra do not instantiate hundreds of D3 charts
    // until their collapsible spectrum section is opened.
    std::map<std::string,std::string> references;
    ReferenceLineInfo matched_lines, missed_lines, other_unmatched_lines,
                      excluded_lines, public_peak_lines,
                      observable_peak_lines, extent_lines;
    for( ReferenceLineInfo *lines : {&matched_lines, &missed_lines,
          &other_unmatched_lines, &excluded_lines,
          &public_peak_lines, &observable_peak_lines, &extent_lines} )
    {
      lines->m_validity = ReferenceLineInfo::InputValidity::Valid;
      lines->m_source_type = ReferenceLineInfo::SourceType::OneOffSrcLines;
      lines->m_has_coincidences = false;
    }
    const std::vector<ExpectedPhotopeakInfo> &all_truth
      = info.expected_signal_photopeaks;
    const std::vector<ExpectedPhotopeakInfo> detector_scoring_truth
      = PeakFitImproveData::filter_photopeaks_for_scoring(
          info.expected_signal_photopeaks, info.det_type );
    const std::vector<ExpectedPhotopeakInfo> scoring_truth
      = ReportDetail::filter_truth_for_spectroscopic_extent( detector_scoring_truth, pd );
    const std::vector<PeakDef> empty_observable_peaks;
    const std::vector<PeakDef> &observable_peaks
      = fit ? fit->observable_peaks : empty_observable_peaks;
    const CandidatePeak_GA::CandidatePeakScore raw_candidate_score
      = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
          observable_peaks, scoring_truth, info.det_type );
    CandidatePeak_GA::CandidatePeakScore corrected_candidate_score
      = raw_candidate_score;
    CandidatePeak_GA::correct_score_for_escape_peaks(
      corrected_candidate_score, scoring_truth );
    const auto same_truth = []( const ExpectedPhotopeakInfo &lhs,
                                const ExpectedPhotopeakInfo &rhs ) -> bool {
      return std::fabs(lhs.effective_energy - rhs.effective_energy) < 1.0e-6
        && std::fabs(lhs.roi_lower - rhs.roi_lower) < 1.0e-6
        && std::fabs(lhs.roi_upper - rhs.roi_upper) < 1.0e-6;
    };
    const auto contains_truth = [&same_truth](
      const std::vector<ExpectedPhotopeakInfo> &haystack,
      const ExpectedPhotopeakInfo &needle ) -> bool {
      return std::any_of( haystack.begin(), haystack.end(),
        [&same_truth, &needle]( const ExpectedPhotopeakInfo &candidate ) {
          return same_truth( candidate, needle );
        } );
    };
    const auto is_scored_truth = [&scoring_truth]( const ExpectedPhotopeakInfo &truth ) -> bool {
      return std::any_of( scoring_truth.begin(), scoring_truth.end(),
        [&truth]( const ExpectedPhotopeakInfo &candidate ) {
          return std::fabs(candidate.effective_energy - truth.effective_energy) < 1.0e-6
            && std::fabs(candidate.roi_lower - truth.roi_lower) < 1.0e-6
            && std::fabs(candidate.roi_upper - truth.roi_upper) < 1.0e-6;
        } );
    };
    for( const ExpectedPhotopeakInfo &truth : all_truth )
    {
      if( !is_scored_truth(truth) )
      {
        append_reference_line( excluded_lines, truth.effective_energy,
          Wt::WColor(110,110,110), "Truth excluded by objective filter" );
        continue;
      }
      bool matched = false;
      if( fit )
        for( const PeakDef &peak : fit->observable_peaks )
          matched = matched || objective_truth_match( peak, truth );
      const bool objective_miss = !matched && contains_truth(
        corrected_candidate_score.def_expected_but_not_detected, truth );
      ReferenceLineInfo &lines = matched ? matched_lines
        : (objective_miss ? missed_lines : other_unmatched_lines);
      append_reference_line( lines, truth.effective_energy,
        matched ? Wt::WColor(0,128,0)
                : (objective_miss ? Wt::WColor(220,0,0) : Wt::WColor(190,120,0)),
        matched ? "Matched truth"
                : (objective_miss ? "Objective missed truth"
                                  : "Unmatched truth not counted as definite miss") );
    }
    if( spectrum_extent.second > spectrum_extent.first )
    {
      append_reference_line( extent_lines, spectrum_extent.first, Wt::WColor(220,120,0),
                             "Lower spectroscopic extent" );
      append_reference_line( extent_lines, spectrum_extent.second, Wt::WColor(220,120,0),
                             "Upper spectroscopic extent" );
    }
    if( fit )
    {
      for( const PeakDef &peak : fit->fit_peaks )
        append_reference_line( public_peak_lines, peak.mean(), Wt::WColor(0,145,170),
                               "Combined public fitted mean" );
      for( const PeakDef &peak : fit->observable_peaks )
        append_reference_line( observable_peak_lines, peak.mean(), Wt::WColor(180,0,180),
                               "Objective observable fitted mean" );
    }
    for( const std::pair<std::string,ReferenceLineInfo *> entry :
         std::vector<std::pair<std::string,ReferenceLineInfo *>>{
           {"Matched truth",&matched_lines}, {"Missed truth",&missed_lines},
           {"Other unmatched truth",&other_unmatched_lines},
           {"Truth not scored",&excluded_lines},
           {"Spectroscopic extent",&extent_lines},
           {"Combined public fitted means",&public_peak_lines},
           {"Objective observable fitted means",&observable_peak_lines} } )
    {
      if( !entry.second->m_ref_lines.empty() )
      {
        std::string json;
        entry.second->toJson( json );
        references[entry.first] = json;
      }
    }

    float x_min = pd.foreground ? pd.foreground->gamma_energy_min() : 0.0f;
    float x_max = pd.foreground ? pd.foreground->gamma_energy_max() : 3000.0f;
    const D3SpectrumExport::D3SpectrumChartOptions chart_options(
      "FitPeaksForNuclides review result", "Energy (keV)", "Counts/channel", "",
      true, false, false, true, true, false, false, false, false, false, false,
      false, false, false, x_min, x_max, references );

    D3SpectrumExport::D3SpectrumOptions foreground_options;
    foreground_options.line_color = "black";
    foreground_options.title = "Measured foreground";
    foreground_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
    if( fit )
    {
      std::vector<std::shared_ptr<const PeakDef>> ptrs;
      for( const PeakDef &peak : fit->fit_peaks )
        ptrs.push_back( std::make_shared<PeakDef>(peak) );
      foreground_options.peaks_json = PeakDef::peak_json(
        ptrs, pd.foreground, Wt::WColor(47,111,221), 254 );
    }

    html << "<script>function " << init_fn << "(){if(window['" << init_fn
         << "_done'])return;window['" << init_fn << "_done']=true;\n";
    D3SpectrumExport::write_js_for_chart(
      html, chart_id, chart_options.m_dataTitle,
      chart_options.m_xAxisTitle, chart_options.m_yAxisTitle );
    std::vector<std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions>> measurements;
    if( pd.foreground )
      measurements.emplace_back( pd.foreground.get(), foreground_options );
    if( pd.background )
    {
      D3SpectrumExport::D3SpectrumOptions options;
      options.line_color = "steelblue";
      options.title = "Supplied background (live-time scaled)";
      options.display_scale_factor = (pd.background->live_time() > 0.0)
        ? pd.foreground->live_time()/pd.background->live_time() : 1.0;
      options.spectrum_type = SpecUtils::SpectrumType::Background;
      measurements.emplace_back( pd.background.get(), options );
    }
    D3SpectrumExport::write_and_set_data_for_chart( html, chart_id, measurements );
    D3SpectrumExport::write_set_options_for_chart( html, chart_id, chart_options );
    html << "spec_chart_" << chart_id << ".setReferenceLines(reference_lines_" << chart_id
         << ");spec_chart_" << chart_id << ".chart.style.width='96%';spec_chart_"
         << chart_id << ".chart.style.height='430px';spec_chart_" << chart_id
         << ".handleResize();var chartEl=document.getElementById('" << chart_id
         << "');var section=document.getElementById('" << summary.id
         << "');chartEl.addEventListener('click',function(event){var path=event.target;if(!path||!path.classList||(!path.classList.contains('peakFill')&&!path.classList.contains('peakOutline')))return;var energy=parseFloat(path.getAttribute('data-energy'));if(!isFinite(energy))return;var rows=section.querySelectorAll('.roi-review-row');var selected=null;for(var i=0;i<rows.length;i++){rows[i].classList.remove('roi-selected');var lower=parseFloat(rows[i].getAttribute('data-lower'));var upper=parseFloat(rows[i].getAttribute('data-upper'));if(energy>=lower&&energy<=upper)selected=rows[i];}if(selected){selected.classList.add('roi-selected');selected.scrollIntoView(false);}},true);}\n"
         << "document.getElementById('" << summary.id
         << "').addEventListener('toggle',function(){if(this.open)" << init_fn << "();});</script>";

    html << "<h3>Objective components</h3><table><tr><th>Scalar cost</th><th>Area cost</th><th>Find reward</th><th>Candidate reward</th><th>Miss fraction</th><th>Missed definitely wanted</th><th>Spurious</th><th>Raw bg penalty</th><th>Weighted bg contribution</th></tr><tr><td>"
         << summary.cost << "</td><td>" << ev.accuracy.area_cost << "</td><td>"
         << ev.accuracy.find_reward << "</td><td>" << ev.accuracy.candidate_reward
         << "</td><td>" << ev.accuracy.miss_fraction << "</td><td>"
         << ev.accuracy.missed_definitely_wanted << "</td><td>" << ev.accuracy.extra_peaks
         << "</td><td>" << ev.background_penalty << "</td><td>"
         << (sm_background_fit_penalty_weight*ev.background_penalty)
         << "</td></tr></table>";

    if( sm_do_background_fit_trial )
    {
      html << "<h3>Background false-positive trial</h3><p>Normalized significance is live-time adjusted before the per-peak cap. A non-empty suppression reason means the peak did not contribute to the raw penalty.</p>";
      if( !ev.background_detail.error_message.empty() )
        html << "<div class=\"error noteable\" data-note-context=\""
             << ReportDetail::html_escape(note_context + " background trial error")
             << "\">" << ReportDetail::html_escape(ev.background_detail.error_message)
             << "</div>";
      html << "<table><tr><th>Peak</th><th>Source</th><th>Mean</th><th>Area</th><th>Normalized significance</th><th>Suppression / contribution</th></tr>";
      for( size_t bg_index = 0;
           bg_index < ev.background_detail.source_attributed_peaks.size(); ++bg_index )
      {
        const PeakDef &peak = ev.background_detail.source_attributed_peaks[bg_index];
        const double significance
          = (bg_index < ev.background_detail.normalized_significances.size())
              ? ev.background_detail.normalized_significances[bg_index] : 0.0;
        const std::string reason
          = (bg_index < ev.background_detail.suppression_reasons.size())
              ? ev.background_detail.suppression_reasons[bg_index] : std::string();
        html << "<tr class=\"noteable\" data-note-context=\""
             << ReportDetail::html_escape(note_context + " background peak "
                  + std::to_string(bg_index+1)) << "\"><td>" << (bg_index+1)
             << "</td><td>" << ReportDetail::html_escape(peak.sourceName())
             << "</td><td>" << peak.mean() << "</td><td>" << peak.amplitude()
             << "</td><td>" << significance << "</td><td>"
             << ReportDetail::html_escape(reason.empty() ? "contributes" : reason)
             << "</td></tr>";
      }
      if( ev.background_detail.source_attributed_peaks.empty() )
        html << "<tr><td colspan=\"6\">No source-attributed background peaks.</td></tr>";
      html << "</table>";
    }

    html << "<h3>Final public ROI diagnostics</h3><p>Widths and fitted-object bounds come from the shared <code>PeakContinuum</code> objects in <code>fit_peaks</code>. "
         << "Model Pearson RMS and fitted-continuum-vs-global-SNIP RMS are report-only diagnostics; neither changes the objective.</p>"
         << "<table><tr><th>ROI</th><th>Lower</th><th>Upper</th><th>Width keV</th><th>Channels</th><th>Width/FWHM</th><th>Peaks</th><th>OffsetType</th><th>Step/reference position</th><th>Model Pearson RMS</th><th>SNIP continuum RMS</th></tr>";
    for( size_t roi_index = 0; roi_index < summary.rois.size(); ++roi_index )
    {
      const RoiReview &roi = summary.rois[roi_index];
      html << "<tr class=\"noteable roi-review-row\" data-lower=\"" << roi.lower
           << "\" data-upper=\"" << roi.upper << "\" data-note-context=\""
           << ReportDetail::html_escape(note_context + " ROI " + std::to_string(roi_index+1))
           << "\"><td>ROI " << (roi_index+1) << "</td><td>" << roi.lower
           << "</td><td>" << roi.upper << "</td><td>" << (roi.upper-roi.lower)
           << "</td><td>" << roi.channels << "</td><td>" << roi.width_fwhm
           << "</td><td>" << roi.peaks.size() << "</td><td>"
           << ReportDetail::html_escape(PeakContinuum::offset_type_str(roi.continuum->type()))
           << "</td><td>" << (PeakContinuum::is_step_continuum(roi.continuum->type())
                ? std::to_string(roi.continuum->referenceEnergy()) : std::string("—"))
           << "</td><td>" << roi.model_pearson_rms << "</td><td>"
           << roi.snip_continuum_rms << "</td></tr>";
    }
    if( summary.rois.empty() )
      html << "<tr><td colspan=\"11\">No public fitted ROIs. The measured chart and truth table remain available for review.</td></tr>";
    html << "</table>";

    html << "<h3>Fitted peaks</h3><p>This table uses combined public <code>fit_peaks</code>. Observable association is a report-only nearest-position association to <code>observable_peaks</code>; the objective itself uses the observable vector directly.</p>"
         << "<table><tr><th>Peak</th><th>Source/provenance</th><th>Mean</th><th>Mean uncert</th><th>Area</th><th>Area uncert</th><th>Significance</th><th>ROI</th><th>Observable</th></tr>";
    if( fit )
    {
      for( size_t peak_index = 0; peak_index < fit->fit_peaks.size(); ++peak_index )
      {
        const PeakDef &peak = fit->fit_peaks[peak_index];
        size_t roi_number = 0;
        for( size_t roi_index = 0; roi_index < summary.rois.size(); ++roi_index )
          if( summary.rois[roi_index].continuum == peak.continuum().get() )
            roi_number = roi_index + 1;
        bool observable = false;
        for( const PeakDef &obs : fit->observable_peaks )
          observable = observable || (std::fabs(obs.mean()-peak.mean())
            <= std::max(0.01,0.25*std::min(obs.fwhm(),peak.fwhm())));
        const double uncert = peak.amplitudeUncert();
        const double significance = (uncert > 0.0) ? peak.amplitude()/uncert
          : ((peak.amplitude() > 0.0) ? std::sqrt(peak.amplitude()) : 0.0);
        std::string provenance = peak.sourceName();
        if( provenance.empty() ) provenance = peak.userLabel();
        if( provenance.empty() ) provenance = "unattributed/floating";
        html << "<tr class=\"noteable\" data-note-context=\""
             << ReportDetail::html_escape(note_context + " fitted peak " + std::to_string(peak_index+1))
             << "\"><td>" << (peak_index+1) << "</td><td>"
             << ReportDetail::html_escape(provenance) << "</td><td>" << peak.mean()
             << "</td><td>" << peak.meanUncert() << "</td><td>" << peak.amplitude()
             << "</td><td>" << peak.amplitudeUncert() << "</td><td>" << significance
             << "</td><td>" << (roi_number ? std::to_string(roi_number) : "—")
             << "</td><td>" << (observable ? "yes" : "no") << "</td></tr>";
      }
    }
    if( !fit || fit->fit_peaks.empty() )
      html << "<tr><td colspan=\"9\">No fitted public peaks.</td></tr>";
    html << "</table>";

    html << "<h3>Objective observable peaks</h3><p>This is the exact peak vector used by all displayed scalar objective components and Candidate miss/spurious classifications. It may differ from the combined public representation above.</p>"
         << "<table><tr><th>Peak</th><th>Source/provenance</th><th>Mean</th><th>Mean uncert</th><th>Area</th><th>Area uncert</th><th>Significance</th></tr>";
    if( fit )
    {
      for( size_t peak_index = 0; peak_index < fit->observable_peaks.size(); ++peak_index )
      {
        const PeakDef &peak = fit->observable_peaks[peak_index];
        const double uncert = peak.amplitudeUncert();
        const double significance = (uncert > 0.0) ? peak.amplitude()/uncert
          : ((peak.amplitude() > 0.0) ? std::sqrt(peak.amplitude()) : 0.0);
        std::string provenance = peak.sourceName();
        if( provenance.empty() ) provenance = peak.userLabel();
        if( provenance.empty() ) provenance = "unattributed/floating";
        html << "<tr class=\"noteable\" data-note-context=\""
             << ReportDetail::html_escape(note_context + " objective observable peak "
                  + std::to_string(peak_index+1)) << "\"><td>" << (peak_index+1)
             << "</td><td>" << ReportDetail::html_escape(provenance) << "</td><td>"
             << peak.mean() << "</td><td>" << peak.meanUncert() << "</td><td>"
             << peak.amplitude() << "</td><td>" << peak.amplitudeUncert()
             << "</td><td>" << significance << "</td></tr>";
      }
    }
    if( !fit || fit->observable_peaks.empty() )
      html << "<tr><td colspan=\"7\">No objective observable peaks.</td></tr>";
    html << "</table>";

    html << "<h3>Truth comparison and explicit spurious peaks</h3><p>Truth detection and spurious classification use the CandidatePeak objective's expected-ROI / fitted ±1σ overlap rule. The displayed best fitted match is the nearest among objective-compatible peaks. Area pull uses the same 5% truth-relative floor as the area objective.</p>"
         << "<table><tr><th>Truth energy</th><th>Truth area</th><th>Truth σ over background</th><th>Matched fitted peak</th><th>Fitted area</th><th>Area error (%)</th><th>Energy error</th><th>Normalized area error</th><th>Classification / reason</th></tr>";
    for( size_t truth_index = 0; truth_index < all_truth.size(); ++truth_index )
    {
      const ExpectedPhotopeakInfo &truth = all_truth[truth_index];
      const bool scored = is_scored_truth( truth );
      const PeakDef *best = nullptr;
      if( fit && scored )
      {
        for( const PeakDef &peak : fit->observable_peaks )
          if( objective_truth_match(peak,truth)
              && (!best || std::fabs(peak.mean()-truth.effective_energy)
                    < std::fabs(best->mean()-truth.effective_energy)) )
            best = &peak;
      }
      double area_pull = 0.0;
      if( best )
      {
        constexpr double f_rel = 0.05;
        const double area = std::max(1.0,truth.peak_area);
        const double sigma = std::sqrt(area + (f_rel*area)*(f_rel*area));
        area_pull = (best->amplitude()-truth.peak_area)/sigma;
      }
      const bool objective_miss = scored && !best && contains_truth(
        corrected_candidate_score.def_expected_but_not_detected, truth );
      const bool escape_corrected = scored && !best
        && contains_truth( raw_candidate_score.def_expected_but_not_detected, truth )
        && !objective_miss;
      html << "<tr class=\"noteable" << (objective_miss ? " miss" : "")
           << "\" data-note-context=\"" << ReportDetail::html_escape(
                note_context + " truth " + std::to_string(truth.effective_energy) + " keV")
           << "\"><td>" << truth.effective_energy << "</td><td>" << truth.peak_area
           << "</td><td>" << truth.nsigma_over_background << "</td><td>"
           << (best ? std::to_string(best->mean()) : std::string("—")) << "</td><td>"
           << (best ? std::to_string(best->amplitude()) : std::string("—")) << "</td><td>"
           << (best && (truth.peak_area != 0.0)
                 ? std::to_string(100.0*(best->amplitude()-truth.peak_area)/truth.peak_area)
                 : std::string("—")) << "</td><td>"
           << (best ? std::to_string(best->mean()-truth.effective_energy) : std::string("—"))
           << "</td><td>" << (best ? std::to_string(area_pull) : std::string("—"))
           << "</td><td>" << (!scored ? "not scored: detector-specific low-energy filter"
                : (best ? "matched by objective rule"
                : (objective_miss ? "MISSED definitely-wanted truth peak"
                : (escape_corrected ? "not penalized: objective escape-peak correction"
                                    : "not detected; below definitely-wanted threshold"))))
           << "</td></tr>";
    }
    if( all_truth.empty() )
      html << "<tr><td colspan=\"9\">No truth peaks are defined for this spectrum.</td></tr>";

    bool used_511_exemption = false;
    if( fit )
    {
      for( const PeakDef &peak : fit->observable_peaks )
      {
        bool matched = false;
        for( const ExpectedPhotopeakInfo &truth : scoring_truth )
          matched = matched || objective_truth_match(peak,truth);
        if( matched )
          continue;
        const bool exempt_511 = !used_511_exemption && (peak.mean() > 508.0) && (peak.mean() < 514.0);
        used_511_exemption = used_511_exemption || exempt_511;
        html << "<tr class=\"noteable " << (exempt_511 ? "" : "spurious")
             << "\" data-note-context=\"" << ReportDetail::html_escape(
                  note_context + " spurious fitted peak " + std::to_string(peak.mean()) + " keV")
             << "\"><td>—</td><td>—</td><td>—</td><td>" << peak.mean()
             << "</td><td>—</td><td>—</td><td>—</td><td>—</td><td>"
             << (exempt_511 ? "unexpected 511-keV peak; objective exempts the first one"
                            : "SPURIOUS by objective CandidatePeak rule")
             << "</td></tr>";
      }
    }
    html << "</table>";

    if( fit )
    {
      html << "<details><summary>Private solver/result representation audit</summary><p>Public combined fit peaks: "
           << fit->fit_peaks.size() << "; public uncombined peaks: "
           << fit->uncombined_fit_peaks.size() << "; objective observable peaks: "
           << fit->observable_peaks.size() << "; private solution peaks (may include hidden R6/NORM nuisance peaks): "
           << fit->solution.m_peaks_without_back_sub.size() << "; final solver ROI ranges: "
           << fit->solution.m_final_roi_ranges_in_spectrum_cal.size()
           << ". Private peaks are disclosed here but are not presented as production output.</p></details>";
    }
    html << "</div></details>\n";

    if( !n42_output_dir.empty() && fit )
    {
      const std::string detector_dir = SpecUtils::append_path(
        n42_output_dir, info.detector_name );
      if( !SpecUtils::is_directory(n42_output_dir) )
        SpecUtils::create_directory( n42_output_dir );
      if( !SpecUtils::is_directory(detector_dir) )
        SpecUtils::create_directory( detector_dir );
      SpecMeas output;
      output.remove_measurements( output.measurements() );
      std::shared_ptr<SpecUtils::Measurement> foreground
        = std::make_shared<SpecUtils::Measurement>( *pd.foreground );
      foreground->set_sample_number( 1 );
      output.add_measurement( foreground, false );
      std::deque<std::shared_ptr<const PeakDef>> peaks;
      for( const PeakDef &peak : fit->observable_peaks )
        peaks.push_back( std::make_shared<PeakDef>(peak) );
      output.setPeaks( peaks, {1} );
      const std::string output_path = SpecUtils::append_path(
        detector_dir, src.src_name + "_nuclide_config_ga.n42" );
      output.save2012N42File( output_path, [&output_path](){
        std::cerr << "Failed to write '" << output_path << "'" << std::endl;
      } );
    }
  }

  html << R"html(
<div id="remark-dialog" class="modal" role="dialog" aria-modal="true" aria-labelledby="remark-title" hidden><div class="modal-card"><h3 id="remark-title">Add review remark</h3><p id="remark-context"></p><textarea id="remark-text"></textarea><div><button type="button" id="remark-cancel">Cancel</button><button type="button" id="remark-accept">Accept</button></div></div></div>
<div id="remark-export-dialog" class="modal" role="dialog" aria-modal="true" aria-labelledby="remark-export-title" hidden><div class="modal-card"><h3 id="remark-export-title">Aggregated review remarks</h3><p>Copy this text back into the Codex session.</p><textarea id="remark-export-text" style="height:50vh"></textarea><div><button type="button" id="remark-copy">Copy</button><button type="button" id="remark-export-close">Close</button></div></div></div>
<script>
(function(){
  'use strict';
  var storageKey='interspec-fit-review:' + location.pathname;
  var activeContext='';
  var activeSelection='';
  var dialog=document.getElementById('remark-dialog');
  function selectedText(){var s=window.getSelection ? String(window.getSelection()) : '';return s.replace(/^\s+|\s+$/g,'');}
  function showModal(el){el.hidden=false;}function hideModal(el){el.hidden=true;}
  function openRemark(el){activeContext=el.getAttribute('data-note-context')||'report';activeSelection=selectedText();document.getElementById('remark-context').textContent=activeContext+(activeSelection?' — selected: '+activeSelection:'');document.getElementById('remark-text').value='';showModal(dialog);document.getElementById('remark-text').focus();}
  Array.prototype.forEach.call(document.querySelectorAll('.noteable,.add-note'),function(el){el.addEventListener('click',function(event){if(event.target.tagName==='A')return;openRemark(el);});});
  document.getElementById('remark-accept').addEventListener('click',function(){var text=document.getElementById('remark-text').value.replace(/^\s+|\s+$/g,'');if(!text)return;var notes=[];try{notes=JSON.parse(localStorage.getItem(storageKey)||'[]');}catch(e){}notes.push({context:activeContext,selection:activeSelection,remark:text,recorded:new Date().toISOString()});localStorage.setItem(storageKey,JSON.stringify(notes));hideModal(dialog);});
  document.getElementById('remark-cancel').addEventListener('click',function(){hideModal(dialog);});
  function exportNotes(){var notes=[];try{notes=JSON.parse(localStorage.getItem(storageKey)||'[]');}catch(e){}var lines=['InterSpec FitPeaksForNuclides review remarks'];for(var i=0;i<notes.length;i++){lines.push('\n['+(i+1)+'] '+notes[i].context+(notes[i].selection?'\nSelected: '+notes[i].selection:'')+'\nRemark: '+notes[i].remark+'\nRecorded: '+notes[i].recorded);}var text=lines.join('\n');document.getElementById('remark-export-text').value=text;showModal(document.getElementById('remark-export-dialog'));}
  document.getElementById('export-remarks').addEventListener('click',exportNotes);
  document.getElementById('remark-copy').addEventListener('click',function(){var area=document.getElementById('remark-export-text');area.select();document.execCommand('copy');});
  document.getElementById('remark-export-close').addEventListener('click',function(){hideModal(document.getElementById('remark-export-dialog'));});
  document.addEventListener('keydown',function(event){if(event.key==='Escape'){hideModal(dialog);hideModal(document.getElementById('remark-export-dialog'));}});
  function updateSummary(){var body=document.querySelector('#summary-table tbody');var rows=Array.prototype.slice.call(body.rows);var key=document.getElementById('summary-sort').value;var query=document.getElementById('summary-filter').value.toLowerCase();var state=document.getElementById('state-filter').value;rows.sort(function(a,b){if(key==='state'){var order={failure:2,empty:1,success:0};return order[b.dataset.state]-order[a.dataset.state];}return parseFloat(b.dataset[key]||0)-parseFloat(a.dataset[key]||0);});for(var i=0;i<rows.length;i++){body.appendChild(rows[i]);var text=rows[i].textContent.toLowerCase();rows[i].style.display=((!query||text.indexOf(query)>=0)&&(!state||rows[i].dataset.state===state))?'':'none';}}
  document.getElementById('summary-sort').addEventListener('change',updateSummary);document.getElementById('summary-filter').addEventListener('input',updateSummary);document.getElementById('state-filter').addEventListener('change',updateSummary);updateSummary();
  function openHashSpectrum(){if(!location.hash)return;var id;try{id=decodeURIComponent(location.hash.slice(1));}catch(e){return;}var target=document.getElementById(id);if(target&&target.tagName.toLowerCase()==='details'&&target.classList.contains('spectrum')){target.open=true;window.setTimeout(function(){target.scrollIntoView();},0);}}
  window.addEventListener('hashchange',openHashSpectrum);openHashSpectrum();
})();
</script></body></html>)html";
}


ConfigEvaluation evaluate_for_report(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const PeakFitForNuclideConfig &config,
  const BackgroundMode bg_mode )
{
  ConfigEvaluation evaluation;
  evaluation.spectra.resize( precomputed.size() );
  std::vector<std::pair<double,double>> costs( precomputed.size(), {0.0,0.0} );
  boost::asio::thread_pool pool(
    std::max<size_t>(1, PeakFitImprove::sm_num_optimization_threads) );

  for( size_t index = 0; index < precomputed.size(); ++index )
  {
    boost::asio::post( pool,
      [index, &precomputed, &config, bg_mode, &evaluation, &costs]() {
      const PrecomputedNuclideData &pd = precomputed[index];
      SpectrumEvaluation &record = evaluation.spectra[index];
      record.spectrum_id = ReportDetail::canonical_spectrum_key( pd );
      record.anchor_id = ReportDetail::stable_spectrum_id( pd );
      Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;
      if( bg_mode == BackgroundMode::NoBackgroundFitNorm )
        options |= FitPeaksForNuclides::FitNormBkgrndPeaks;
      if( sm_disable_auto_interferer_fit )
        options |= FitPeaksForNuclides::DisableAutoInterfererFit;

      try
      {
        const std::vector<std::shared_ptr<const PeakDef>> user_peaks;
        PeakFitResult result = FitPeaksForNuclides::fit_peaks_for_nuclides(
          pd.auto_search_peaks, pd.foreground, pd.sources, user_peaks,
          pd.background, pd.drf, options, config, pd.peak_fit_prefs );
        FitPeaksForNuclides::detail::take_roi_boundary_shadow_diagnostics();
        const bool success
          = (result.status == RelActCalcAuto::RelActAutoSolution::Status::Success);
        record.has_fit_result = true;
        record.accuracy = ReportDetail::score_observable_peaks(
          pd, result.observable_peaks );
        record.legitimate_empty = success && result.observable_peaks.empty()
          && (record.accuracy.missed_definitely_wanted == 0);
        record.mechanical_failure = !success;
        record.error_message = result.error_message;
        using Status = RelActCalcAuto::RelActAutoSolution::Status;
        switch( result.status )
        {
          case Status::Success: record.status = "Success"; break;
          case Status::NotInitiated: record.status = "NotInitiated"; break;
          case Status::FailedToSetupProblem: record.status = "FailedToSetupProblem"; break;
          case Status::FailToSolveProblem: record.status = "FailToSolveProblem"; break;
          case Status::UserCanceled: record.status = "UserCanceled"; break;
        }
        if( success )
          record.background_penalty = compute_background_fit_penalty(
            pd, config, bg_mode, &record.background_detail );
        FitPeaksForNuclides::detail::take_roi_boundary_shadow_diagnostics();
        costs[index] = {record.accuracy.cost, record.background_penalty};
        record.fit_result = std::move( result );
      }catch( const std::exception &e )
      {
        FitPeaksForNuclides::detail::take_roi_boundary_shadow_diagnostics();
        record.exception = true;
        record.mechanical_failure = true;
        record.status = "EXCEPTION";
        record.error_message = e.what();
        record.accuracy = ReportDetail::score_observable_peaks( pd, {} );
        costs[index] = {record.accuracy.cost,0.0};
      }
    } );
  }
  pool.join();

  for( size_t index = 0; index < evaluation.spectra.size(); ++index )
  {
    const SpectrumEvaluation &record = evaluation.spectra[index];
    evaluation.total_fg += costs[index].first;
    evaluation.total_bg_raw += costs[index].second;
    evaluation.accuracy_totals.area_cost += record.accuracy.area_cost;
    evaluation.accuracy_totals.find_reward += record.accuracy.find_reward;
    evaluation.accuracy_totals.candidate_reward += record.accuracy.candidate_reward;
    evaluation.accuracy_totals.miss_fraction += record.accuracy.miss_fraction;
    evaluation.accuracy_totals.missed_definitely_wanted
      += record.accuracy.missed_definitely_wanted;
    evaluation.accuracy_totals.extra_peaks += record.accuracy.extra_peaks;
    if( record.mechanical_failure ) evaluation.mechanical_failures += 1;
    else if( record.legitimate_empty ) evaluation.legitimate_empties += 1;
    else evaluation.successes += 1;
  }
  evaluation.total_cost = evaluation.total_fg
    + sm_background_fit_penalty_weight*evaluation.total_bg_raw;
  return evaluation;
}
}//anonymous namespace


void write_config_evaluation_tsv(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const ConfigEvaluation &evaluation,
  const std::string &filename )
{
  ReportDetail::validate_evaluation_coverage( precomputed.size(), evaluation );
  std::ifstream existing( filename );
  if( existing.good() )
    throw std::runtime_error( "Refusing to overwrite evaluation TSV: " + filename );

  std::ofstream out( filename );
  if( !out.good() )
    throw std::runtime_error( "Could not open evaluation TSV output: " + filename );
  out << std::setprecision(17);

  const auto tsv_text = []( std::string value ) {
    std::replace( value.begin(), value.end(), '\t', ' ' );
    std::replace( value.begin(), value.end(), '\r', ' ' );
    std::replace( value.begin(), value.end(), '\n', ' ' );
    return value;
  };
  const auto join_warnings = [&tsv_text]( const std::vector<std::string> &warnings ) {
    std::ostringstream value;
    for( size_t index = 0; index < warnings.size(); ++index )
    {
      if( index )
        value << " | ";
      value << tsv_text( warnings[index] );
    }
    return value.str();
  };

  out << "spectrum_key\tstable_id\tdetector\tcity\tlive_time\tsource\tfile_base_path"
      << "\tstatus\tlegitimate_empty\tmechanical_failure\tscalar_cost\tarea_cost"
      << "\tfind_reward\tcandidate_reward\tmiss_fraction\tdefinite_miss_count\tspurious_count"
      << "\tfitted_peak_count\tobservable_peak_count\troi_count\tmax_roi_width_kev"
      << "\tmax_roi_width_channels\tmax_roi_width_fwhm\tmax_peaks_per_roi\tstep_count"
      << "\tcontinuum_types\tmax_pearson_rms\tmax_snip_mismatch\tmax_report_area_pull"
      << "\tbackground_contribution\twarnings\terror\n";

  for( size_t index = 0; index < precomputed.size(); ++index )
  {
    const PrecomputedNuclideData &pd = precomputed[index];
    const SpectrumEvaluation &ev = evaluation.spectra[index];
    if( !pd.src_info || ev.spectrum_id != ReportDetail::canonical_spectrum_key(pd) )
      throw std::runtime_error( "Evaluation TSV record is not associated with its spectrum" );

    double max_width_kev = 0.0;
    size_t max_width_channels = 0;
    double max_width_fwhm = 0.0;
    size_t max_peaks_per_roi = 0;
    size_t step_count = 0;
    double max_pearson_rms = 0.0;
    double max_snip_mismatch = 0.0;
    std::set<std::string> continuum_types;
    const std::vector<RoiReview> rois = review_rois( pd, ev );
    for( const RoiReview &roi : rois )
    {
      max_width_kev = std::max( max_width_kev, roi.upper - roi.lower );
      max_width_channels = std::max( max_width_channels, roi.channels );
      max_width_fwhm = std::max( max_width_fwhm, roi.width_fwhm );
      max_peaks_per_roi = std::max( max_peaks_per_roi, roi.peaks.size() );
      max_pearson_rms = std::max( max_pearson_rms, roi.model_pearson_rms );
      max_snip_mismatch = std::max( max_snip_mismatch, roi.snip_continuum_rms );
      if( roi.continuum )
      {
        continuum_types.insert( PeakContinuum::offset_type_str(roi.continuum->type()) );
        if( PeakContinuum::is_step_continuum(roi.continuum->type()) )
          step_count += 1;
      }
    }

    double max_area_pull = 0.0;
    if( ev.has_fit_result )
    {
      const std::vector<ExpectedPhotopeakInfo> truth
        = PeakFitImproveData::filter_photopeaks_for_scoring(
            pd.src_info->expected_signal_photopeaks, pd.det_type );
      for( const ExpectedPhotopeakInfo &expected : truth )
      {
        const PeakDef *best = nullptr;
        for( const PeakDef &peak : ev.fit_result.observable_peaks )
        {
          if( objective_truth_match( peak, expected )
              && (!best || std::fabs(peak.mean() - expected.effective_energy)
                    < std::fabs(best->mean() - expected.effective_energy)) )
            best = &peak;
        }
        if( best )
        {
          constexpr double f_rel = 0.05;
          const double area = std::max( 1.0, expected.peak_area );
          const double sigma = std::sqrt( area + (f_rel*area)*(f_rel*area) );
          max_area_pull = std::max( max_area_pull,
            std::fabs(best->amplitude() - expected.peak_area) / sigma );
        }
      }
    }

    std::ostringstream continuum_text;
    for( std::set<std::string>::const_iterator iter = continuum_types.begin();
         iter != continuum_types.end(); ++iter )
    {
      if( iter != continuum_types.begin() )
        continuum_text << ',';
      continuum_text << *iter;
    }
    const DataSrcInfo &info = *pd.src_info;
    const double scalar_cost = ev.accuracy.cost
      + sm_background_fit_penalty_weight*ev.background_penalty;
    out << tsv_text(ev.spectrum_id) << '\t' << tsv_text(ev.anchor_id) << '\t'
        << tsv_text(info.detector_name) << '\t' << tsv_text(info.location_name) << '\t'
        << tsv_text(info.live_time_name) << '\t' << tsv_text(info.src_info.src_name) << '\t'
        << tsv_text(info.src_info.file_base_path) << '\t' << tsv_text(ev.status) << '\t'
        << ev.legitimate_empty << '\t' << ev.mechanical_failure << '\t' << scalar_cost << '\t'
        << ev.accuracy.area_cost << '\t' << ev.accuracy.find_reward << '\t'
        << ev.accuracy.candidate_reward << '\t' << ev.accuracy.miss_fraction << '\t'
        << ev.accuracy.missed_definitely_wanted << '\t' << ev.accuracy.extra_peaks << '\t'
        << (ev.has_fit_result ? ev.fit_result.fit_peaks.size() : 0) << '\t'
        << (ev.has_fit_result ? ev.fit_result.observable_peaks.size() : 0) << '\t'
        << rois.size() << '\t' << max_width_kev << '\t' << max_width_channels << '\t'
        << max_width_fwhm << '\t' << max_peaks_per_roi << '\t' << step_count << '\t'
        << continuum_text.str() << '\t' << max_pearson_rms << '\t' << max_snip_mismatch << '\t'
        << max_area_pull << '\t'
        << (sm_background_fit_penalty_weight*ev.background_penalty) << '\t'
        << join_warnings( ev.has_fit_result ? ev.fit_result.warnings : std::vector<std::string>() )
        << '\t' << tsv_text(ev.error_message) << '\n';
  }
}


void write_results_html_and_n42(
  const std::vector<PrecomputedNuclideData> &precomputed,
  const PeakFitForNuclideConfig &config,
  const NuclideConfigSolution &genes,
  const BackgroundMode bg_mode,
  const std::string &html_filename,
  const std::string &n42_output_dir,
  const ConfigEvaluation *evaluation,
  const std::string &run_metadata )
{
  if( evaluation )
  {
    write_cached_manual_review_report( precomputed, genes, bg_mode, html_filename,
                                       n42_output_dir, *evaluation, run_metadata );
    return;
  }
  const ConfigEvaluation generated_evaluation
    = evaluate_for_report( precomputed, config, bg_mode );
  write_cached_manual_review_report( precomputed, genes, bg_mode, html_filename,
                                     n42_output_dir, generated_evaluation, run_metadata );
  return;

#if 0
  // Superseded reporter retained temporarily for historical comparison.  This is intentionally
  // excluded from compilation; all callers return through the cached manual-review writer above.
  // Open HTML file
  std::ofstream html_output( html_filename );

  try
  {
    D3SpectrumExport::write_html_page_header( html_output, "NuclideConfig GA Results", "InterSpec_resources" );
  }
  catch( std::exception &e )
  {
    cerr << "Error writing HTML header: " << e.what() << endl;
    cerr << "You probably need to symbolically link InterSpec_resources to your CWD." << endl;
    return;
  }

  html_output << "<body>" << endl;
  html_output << "<p><strong>Automatic interferer fit:</strong> "
              << (sm_disable_auto_interferer_fit ? "disabled" : "enabled")
              << "</p>" << endl;
  html_output << "<style>"
    << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
    << "table, th, td{ border: 1px solid black; }" << endl
    << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
    << "</style>" << endl;

  // GA genes that produced every fit below.  Rendered from the authoritative NuclideConfigSolution
  //  serialization (key=value pairs) so the table never drifts from what the GA actually optimizes.
  //  Collapsed by default (<details>) since it is ~55 rows.
  {
    html_output << "<details style=\"width:90vw;margin:15px auto;\">" << endl
                << "<summary style=\"cursor:pointer;font-weight:bold;font-size:110%;\">"
                << "GA genes used to produce these fits</summary>" << endl;
    html_output << "<table class=\"TopLinesTable\"><tr><th>Gene</th><th>Value</th></tr>" << endl;
    std::vector<std::string> gene_lines;
    SpecUtils::split( gene_lines, genes.to_string( "\n" ), "\n" );
    for( const std::string &gl : gene_lines )
    {
      const size_t eq = gl.find( '=' );
      if( eq == std::string::npos )
        continue;
      html_output << "<tr><td style=\"text-align:left;font-family:monospace;\">" << gl.substr( 0, eq )
                  << "</td><td style=\"text-align:right;font-family:monospace;\">" << gl.substr( eq + 1 )
                  << "</td></tr>" << endl;
    }
    html_output << "</table></details>" << endl;
  }

  double total_score = 0.0;
  double total_fg = 0.0;
  double total_bg_raw = 0.0;

  const double num_sigma_contribution = 1.5;

  // Per-spectrum fit-quality, collected so we can list the worst-fitting spectra (the ones most
  //  likely to be catastrophic nuclides) in a summary table.  Ranked by the per-spectrum cost the GA
  //  actually minimizes; failed fits (which would otherwise be skipped) are recorded as the worst.
  struct SpectrumQuality
  {
    std::string src_name;
    std::string detector;
    double cost = 0.0;          // per-spectrum cost (combined_score.final_weight); higher = worse
    std::string status;         // "ok", "FIT FAILED", or "EXCEPTION"
    size_t n_missing = 0;       // expected ("definitely wanted") peaks not found
    size_t n_extra = 0;         // spurious / unexpected peaks
    double frac_area_missing = 0.0;  // fraction of expected signal peak-area not recovered
  };
  std::vector<SpectrumQuality> spectrum_quals;
  spectrum_quals.reserve( precomputed.size() );

  // The per-spectrum re-fit (fit_peaks_for_nuclides) dominates this function's runtime, so we run one
  //  task per spectrum in parallel - the same concurrent calls the GA already makes when scoring
  //  (score_one_spectrum is posted to a thread pool), so they are known thread-safe.  Each task builds
  //  its HTML <fieldset> fragment into a private stream and writes its own N42 file; the fragments are
  //  concatenated in spectrum order afterwards, so the report stays byte-for-byte deterministic.
  struct PerSpectrumResult
  {
    std::string html;            // the per-spectrum HTML fragment (empty for failed/excepted fits)
    double fg = 0.0;             // contribution to total_fg
    double bg_raw = 0.0;         // contribution to total_bg_raw (raw background penalty)
    double score = 0.0;          // contribution to total_score
    SpectrumQuality qual;        // recorded for every spectrum (ok / FIT FAILED / EXCEPTION)
  };
  std::vector<PerSpectrumResult> results( precomputed.size() );

  // Pre-create the N42 output directories serially - concurrent create_directory on a shared detector
  //  subdir would race.
  if( !n42_output_dir.empty() )
  {
    if( !SpecUtils::is_directory( n42_output_dir ) )
      SpecUtils::create_directory( n42_output_dir );
    for( const PrecomputedNuclideData &pd : precomputed )
    {
      const string outdir = SpecUtils::append_path( n42_output_dir, pd.src_info->detector_name );
      if( !SpecUtils::is_directory( outdir ) )
        SpecUtils::create_directory( outdir );
    }
  }//if( !n42_output_dir.empty() )

  boost::asio::thread_pool pool( std::max<size_t>( 1, PeakFitImprove::sm_num_optimization_threads ) );

  for( size_t spec_index = 0; spec_index < precomputed.size(); ++spec_index )
  {
   boost::asio::post( pool, [spec_index, &results, &precomputed, &config, bg_mode, &n42_output_dir, num_sigma_contribution]()
   {
    const PrecomputedNuclideData &pd = precomputed[spec_index];
    PerSpectrumResult &res = results[spec_index];

    // All existing `html_output << ...` lines below write into this private per-spectrum fragment;
    //  shadowing the name lets the body stay unchanged.  Merged into the file in order after join().
    std::ostringstream html_output;
    size_t chart_counter = 0;

    const InjectSourceInfo &src = pd.src_info->src_info;
    const string src_name = src.src_name;
    const string detector_name = pd.src_info->detector_name;

    // Build options based on background mode
    Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options;
    if( bg_mode == BackgroundMode::NoBackgroundFitNorm )
      options |= FitPeaksForNuclides::FitNormBkgrndPeaks;
    if( sm_disable_auto_interferer_fit )
      options |= FitPeaksForNuclides::DisableAutoInterfererFit;

    try
    {
      const std::vector<std::shared_ptr<const PeakDef>> user_peaks;

      const FitPeaksForNuclides::PeakFitResult result = FitPeaksForNuclides::fit_peaks_for_nuclides(
        pd.auto_search_peaks, pd.foreground, pd.sources, user_peaks,
        pd.background, pd.drf, options, config, pd.peak_fit_prefs );

      if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
      {
        // A failed foreground fit is the GA's worst case (scored at sm_fit_failure_penalty); record
        //  it so it shows up in the worst-spectra table rather than silently vanishing.
        SpectrumQuality q;
        q.src_name = src_name;
        q.detector = detector_name;
        q.cost = sm_fit_failure_penalty;
        q.status = "FIT FAILED";
        q.frac_area_missing = 1.0;
        res.qual = q;
        return;
      }

      const RelActCalcAuto::RelActAutoSolution &solution = result.solution;
      const std::vector<PeakDef> &fit_peaks = result.observable_peaks;

      // Score against the det-type-appropriate, low-energy-filtered expected peaks - identical to the
      //  GA objective (score_one_spectrum in PeakFitImprove.cpp), so reported totals match what the GA
      //  minimized.
      const PeakFitUtils::CoarseResolutionType det_type = pd.src_info->det_type;
      const std::vector<ExpectedPhotopeakInfo> scoring_peaks
        = PeakFitImproveData::filter_photopeaks_for_scoring( pd.src_info->expected_signal_photopeaks, det_type );

      // Score the results
      CombinedPeakFitScore combined_score;
      combined_score.final_fit_score = FinalFit_GA::calculate_final_fit_score(
        fit_peaks, scoring_peaks, num_sigma_contribution );
      combined_score.initial_fit_weights = InitialFit_GA::calculate_peak_find_weights(
        fit_peaks, scoring_peaks, num_sigma_contribution, det_type );
      combined_score.candidate_peak_score = CandidatePeak_GA::calculate_candidate_peak_score_for_source(
        fit_peaks, scoring_peaks, det_type );
      CandidatePeak_GA::correct_score_for_escape_peaks(
        combined_score.candidate_peak_score, scoring_peaks );

      // Mirror the GA objective (see score_one_spectrum in PeakFitImprove.cpp): lower-is-better, reward
      // terms subtracted, area-mismatch added, plus the missed-expected-area penalty.
      const double miss_fraction = PeakFitImproveData::missed_def_wanted_area_fraction(
          scoring_peaks, combined_score.candidate_peak_score.def_expected_but_not_detected, det_type );
      combined_score.final_weight = combined_score.final_fit_score.total_weight
                                  - ( combined_score.initial_fit_weights.find_weight
                                      + combined_score.candidate_peak_score.score )
                                  + sm_miss_penalty_weight * miss_fraction;
      const double fg_only_for_source = combined_score.final_weight;
      res.fg = fg_only_for_source;

      // Background false-positive penalty (with diagnostic detail for HTML).
      // compute_background_fit_penalty short-circuits to 0 when
      // sm_do_background_fit_trial is false.
      BackgroundFitDetail bg_detail;
      combined_score.background_fit_penalty
        = compute_background_fit_penalty( pd, config, bg_mode, &bg_detail );
      res.bg_raw = combined_score.background_fit_penalty;
      combined_score.final_weight
        += sm_background_fit_penalty_weight * combined_score.background_fit_penalty;

      res.score = combined_score.final_weight;

      // Record per-spectrum quality for the worst-spectra summary.  frac_area_missing is the same
      //  def-wanted missed-area fraction used by the objective's miss term, so the report column and
      //  the GA cost agree.
      {
        SpectrumQuality q;
        q.src_name = src_name;
        q.detector = detector_name;
        q.cost = combined_score.final_weight;
        q.status = "ok";
        q.n_missing = combined_score.candidate_peak_score.num_def_wanted_not_found;
        q.n_extra = combined_score.candidate_peak_score.num_extra_peaks;
        q.frac_area_missing = miss_fraction;
        res.qual = q;
      }

      // Write N42 file
      if( !n42_output_dir.empty() )
      {
        string outdir = n42_output_dir;
        if( !SpecUtils::is_directory( outdir ) )
          SpecUtils::create_directory( outdir );

        outdir = SpecUtils::append_path( outdir, pd.src_info->detector_name );
        if( !SpecUtils::is_directory( outdir ) )
          SpecUtils::create_directory( outdir );

        const string out_n42 = SpecUtils::append_path( outdir, src_name ) + "_nuclide_config_ga.n42";

        SpecMeas output;
        output.remove_measurements( output.measurements() );

        std::shared_ptr<SpecUtils::Measurement> out_fg = std::make_shared<SpecUtils::Measurement>( *pd.foreground );
        out_fg->set_sample_number( 1 );
        out_fg->set_title( src_name + " - Foreground" );
        output.add_measurement( out_fg, false );

        if( pd.background )
        {
          std::shared_ptr<SpecUtils::Measurement> out_bg = std::make_shared<SpecUtils::Measurement>( *pd.background );
          out_bg->set_sample_number( 2 );
          out_bg->set_title( src_name + " - Background" );
          output.add_measurement( out_bg, false );
        }

        std::deque<std::shared_ptr<const PeakDef>> peaks;
        for( const PeakDef &p : fit_peaks )
          peaks.push_back( std::make_shared<PeakDef>( p ) );
        output.setPeaks( peaks, {1} );

        output.save2012N42File( out_n42, [&](){ cerr << "Failed to write '" << out_n42 << "'" << endl; } );
      }//if( !n42_output_dir.empty() )

      // Write HTML chart
      {
        // Generate reference lines for nuclides
        std::map<std::string, std::string> reference_lines_json;

        const std::vector<Wt::WColor> colors = {
          Wt::WColor( 0, 0, 139 ),
          Wt::WColor( 0, 139, 139 ),
          Wt::WColor( 0, 100, 0 ),
          Wt::WColor( 139, 0, 139 ),
          Wt::WColor( 255, 140, 0 ),
          Wt::WColor( 139, 0, 0 ),
        };

        size_t color_index = 0;
        for( size_t rel_eff_index = 0; rel_eff_index < solution.m_rel_activities.size(); ++rel_eff_index )
        {
          for( const RelActCalcAuto::NuclideRelAct &nuc_info : solution.m_rel_activities[rel_eff_index] )
          {
            RefLineInput input;
            input.m_input_txt = nuc_info.name();
            if( nuc_info.age_was_fit && nuc_info.age > 0 )
              input.m_age = std::to_string( nuc_info.age ) + " s";
            input.m_color = colors[color_index % colors.size()];
            color_index++;
            input.m_showGammas = true;
            input.m_showXrays = true;
            input.m_showAlphas = false;
            input.m_showBetas = false;
            input.m_promptLinesOnly = false;
            input.m_lower_br_cutt_off = 0.0;

            std::shared_ptr<ReferenceLineInfo> ref_info = ReferenceLineInfo::generateRefLineInfo( input );
            if( ref_info && ref_info->m_validity == ReferenceLineInfo::InputValidity::Valid )
            {
              std::string ref_json;
              ref_info->toJson( ref_json );
              reference_lines_json[nuc_info.name()] = ref_json;
            }
          }//for( nuc_info )
        }//for( rel_eff_index )

        // Add red reference lines for missing peaks
        const std::vector<ExpectedPhotopeakInfo> &not_detected = combined_score.candidate_peak_score.def_expected_but_not_detected;
        for( size_t idx = 0; idx < not_detected.size(); ++idx )
        {
          RefLineInput red_input;
          red_input.m_input_txt = std::to_string( not_detected[idx].effective_energy ) + " keV";
          red_input.m_color = Wt::WColor( 255, 0, 0 );
          red_input.m_showGammas = true;
          red_input.m_showXrays = false;
          red_input.m_showAlphas = false;
          red_input.m_showBetas = false;
          red_input.m_promptLinesOnly = false;
          red_input.m_lower_br_cutt_off = 0.0;

          std::shared_ptr<ReferenceLineInfo> red_ref = ReferenceLineInfo::generateRefLineInfo( red_input );
          if( red_ref && red_ref->m_validity == ReferenceLineInfo::InputValidity::Valid )
          {
            std::string json;
            red_ref->toJson( json );
            reference_lines_json["Missing Peak " + std::to_string( idx + 1 )] = json;
          }
        }//for( not_detected )

        // Chart options
        float xMin = 0.0f, xMax = 3000.0f;
        if( !solution.m_final_roi_ranges.empty() )
        {
          xMin = static_cast<float>( solution.m_final_roi_ranges.front().lower_energy );
          xMax = static_cast<float>( solution.m_final_roi_ranges.back().upper_energy );
        }
        else if( pd.foreground )
        {
          xMin = static_cast<float>( pd.foreground->gamma_energy_min() );
          xMax = static_cast<float>( pd.foreground->gamma_energy_max() );
        }

        const string title = src_name + " - NuclideConfig GA Fit";
        D3SpectrumExport::D3SpectrumChartOptions chart_options(
          title, "Energy (keV)", "Counts/Channel",
          "", true, false, false, true, true,
          false, false, false, false, false, false,
          false, false, false, xMin, xMax, reference_lines_json );

        D3SpectrumExport::D3SpectrumOptions fg_opts;
        fg_opts.line_color = "black";
        fg_opts.title = src_name;
        fg_opts.display_scale_factor = 1.0;
        fg_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

        vector<shared_ptr<const PeakDef>> fit_peaks_ptrs;
        for( const PeakDef &p : fit_peaks )
          fit_peaks_ptrs.push_back( make_shared<PeakDef>( p ) );
        fg_opts.peaks_json = PeakDef::peak_json( fit_peaks_ptrs, pd.foreground, Wt::WColor(), false );

        // div / JS-function names must be unique across the whole file, since fragments are
        //  concatenated; key them on the spectrum index plus a per-fragment chart counter.
        const string div_id = "chart_" + std::to_string( spec_index ) + "_" + std::to_string( chart_counter );
        const string resize_fn = "resizeChart_" + std::to_string( spec_index ) + "_" + std::to_string( chart_counter );
        chart_counter++;

        html_output << "<fieldset>" << endl
          << "<legend>" << src_name << " &mdash; " << detector_name
          << " (chi2/dof=" << solution.m_chi2 << "/" << solution.m_dof << ")</legend>" << endl;
        html_output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;
        html_output << "<script>" << endl;

        D3SpectrumExport::write_js_for_chart( html_output, div_id, chart_options.m_dataTitle, chart_options.m_xAxisTitle, chart_options.m_yAxisTitle );

        std::vector<std::pair<const SpecUtils::Measurement *, D3SpectrumExport::D3SpectrumOptions>> measurements;
        measurements.emplace_back( pd.foreground.get(), fg_opts );

        if( pd.background )
        {
          D3SpectrumExport::D3SpectrumOptions bg_opts;
          bg_opts.line_color = "steelblue";
          bg_opts.title = "Background";
          bg_opts.display_scale_factor = pd.foreground->live_time() / pd.background->live_time();
          bg_opts.spectrum_type = SpecUtils::SpectrumType::Background;
          measurements.emplace_back( pd.background.get(), bg_opts );
        }

        D3SpectrumExport::write_and_set_data_for_chart( html_output, div_id, measurements );

        html_output << R"delim(
const )delim" << resize_fn << R"delim( = function(){
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

window.addEventListener('resize', )delim" << resize_fn << R"delim();
)delim" << endl;

        D3SpectrumExport::write_set_options_for_chart( html_output, div_id, chart_options );
        html_output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;
        html_output << resize_fn << "();" << endl;
        html_output << "</script>" << endl;

        // Score summary
        html_output << "<div style=\"margin: 10px auto; max-width: 800px;\">" << endl;
        html_output << "<h4>Score: " << combined_score.final_weight << "</h4>" << endl;
        html_output << "<p>Find weight: " << combined_score.initial_fit_weights.find_weight
                   << ", Final fit: " << combined_score.final_fit_score.total_weight
                   << ", Candidate: " << combined_score.candidate_peak_score.score << "</p>" << endl;
        html_output << "<p>Missing peaks: " << combined_score.candidate_peak_score.num_def_wanted_not_found
                   << ", Extra peaks: " << combined_score.candidate_peak_score.num_extra_peaks << "</p>" << endl;

        // Background-fit penalty section (only shown if computed; suppressed
        // when all sources for this entry are NORM-like or no long_background).
        if( !pd.background_auto_search_peaks.empty() )
        {
          const double weighted = sm_background_fit_penalty_weight * combined_score.background_fit_penalty;
          html_output << "<p><b>Background-fit penalty:</b> "
                     << combined_score.background_fit_penalty
                     << " (weighted contribution to score: " << weighted << ")";
          if( !bg_detail.error_message.empty() )
            html_output << " &mdash; " << bg_detail.error_message;
          html_output << "</p>" << endl;

          if( !bg_detail.source_attributed_peaks.empty() )
          {
            html_output << "<table style=\"border-collapse: collapse; margin: 5px auto;\">"
                       << "<tr><th>Energy (keV)</th><th>Area</th><th>AreaUncert</th>"
                       << "<th>Sig (norm)</th><th>Source</th>"
                       << "<th>Suppressed?</th></tr>" << endl;
            for( size_t i = 0; i < bg_detail.source_attributed_peaks.size(); ++i )
            {
              const PeakDef &bp = bg_detail.source_attributed_peaks[i];
              const SandiaDecay::Nuclide * const nuc = bp.parentNuclide();
              const std::string nuc_sym = nuc ? nuc->symbol : std::string("(none)");
              const double sn = (i < bg_detail.normalized_significances.size())
                                ? bg_detail.normalized_significances[i] : 0.0;
              const std::string &reason = (i < bg_detail.suppression_reasons.size())
                                ? bg_detail.suppression_reasons[i] : std::string();
              const std::string reason_html
                = reason.empty() ? std::string("<b>counted</b>") : reason;
              html_output << "<tr><td>" << bp.mean()
                         << "</td><td>" << bp.amplitude()
                         << "</td><td>" << bp.amplitudeUncert()
                         << "</td><td>" << sn
                         << "</td><td>" << nuc_sym
                         << "</td><td>" << reason_html << "</td></tr>" << endl;
            }
            html_output << "</table>" << endl;
          }

          // ---- Background-fit diagnostic chart -----------------------------
          // Render the long_background spectrum with the source-attributed
          // peaks overlaid, so the user can visually verify the false
          // positives.  Only emit when the bg trial is enabled AND we have
          // peaks to show.
          const std::shared_ptr<const SpecUtils::Measurement> &long_bg
            = pd.src_info ? pd.src_info->src_info.long_background : nullptr;
          if( sm_do_background_fit_trial
              && long_bg
              && !bg_detail.source_attributed_peaks.empty() )
          {
            html_output
              << "<h4 style=\"color:darkred; margin-top: 15px;\">"
              << "Background-fit diagnostic: peaks the source (" << src_name
              << ") attributed to the long-background spectrum (red overlay = false-positive candidates)"
              << "</h4>" << endl;

            // Determine an energy window around the bg-fit peaks for the
            // chart's x-axis.  Pad by 30 keV on either side, clamped to the
            // bg spectrum's extent.
            float bg_xmin = static_cast<float>( long_bg->gamma_energy_min() );
            float bg_xmax = static_cast<float>( long_bg->gamma_energy_max() );
            {
              double e_lo = std::numeric_limits<double>::max();
              double e_hi = std::numeric_limits<double>::lowest();
              for( const PeakDef &bp : bg_detail.source_attributed_peaks )
              {
                e_lo = std::min( e_lo, bp.mean() );
                e_hi = std::max( e_hi, bp.mean() );
              }
              if( e_lo <= e_hi )
              {
                bg_xmin = std::max( bg_xmin, static_cast<float>( e_lo - 30.0 ) );
                bg_xmax = std::min( bg_xmax, static_cast<float>( e_hi + 30.0 ) );
              }
            }

            const string bg_div_id = "chart_" + std::to_string( spec_index ) + "_" + std::to_string( chart_counter );
            const string bg_resize_fn = "resizeChart_" + std::to_string( spec_index ) + "_" + std::to_string( chart_counter );
            chart_counter++;

            const string bg_title = src_name + " — long_background (false positives)";
            D3SpectrumExport::D3SpectrumChartOptions bg_chart_options(
              bg_title, "Energy (keV)", "Counts/Channel",
              "", true, false, false, true, true,
              false, false, false, false, false, false,
              false, false, false, bg_xmin, bg_xmax,
              std::map<std::string,std::string>{} );

            D3SpectrumExport::D3SpectrumOptions bg_only_opts;
            bg_only_opts.line_color = "steelblue";
            bg_only_opts.title = "long_background";
            bg_only_opts.display_scale_factor = 1.0;
            bg_only_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;

            std::vector<std::shared_ptr<const PeakDef>> bg_fit_peaks_ptrs;
            bg_fit_peaks_ptrs.reserve( bg_detail.source_attributed_peaks.size() );
            for( const PeakDef &bp : bg_detail.source_attributed_peaks )
              bg_fit_peaks_ptrs.push_back( std::make_shared<PeakDef>( bp ) );
            bg_only_opts.peaks_json = PeakDef::peak_json(
              bg_fit_peaks_ptrs, long_bg, Wt::WColor("red"), false );

            html_output << "<div id=\"" << bg_div_id
                        << "\" class=\"chart\" oncontextmenu=\"return false;\""
                        << " style=\"border: 1px solid darkred;\"></div>" << endl;
            html_output << "<script>" << endl;

            D3SpectrumExport::write_js_for_chart(
              html_output, bg_div_id,
              bg_chart_options.m_dataTitle,
              bg_chart_options.m_xAxisTitle,
              bg_chart_options.m_yAxisTitle );

            std::vector<std::pair<const SpecUtils::Measurement *, D3SpectrumExport::D3SpectrumOptions>> bg_measurements;
            bg_measurements.emplace_back( long_bg.get(), bg_only_opts );

            D3SpectrumExport::write_and_set_data_for_chart( html_output, bg_div_id, bg_measurements );

            html_output << R"delim(
const )delim" << bg_resize_fn << R"delim( = function(){
  let height = window.innerHeight;
  let width = document.documentElement.clientWidth;
  let el = spec_chart_)delim" << bg_div_id << R"delim(.chart;
  el.style.width = 0.8*width + "px";
  el.style.height = Math.min(400,Math.max(200, Math.min(0.3*width,height-175))) + "px";
  el.style.marginLeft = 0.05*width + "px";
  el.style.marginRight = 0.05*width + "px";
  )delim"
              << "  spec_chart_" << bg_div_id << R"delim(.handleResize();
};

window.addEventListener('resize', )delim" << bg_resize_fn << R"delim();
)delim" << endl;

            D3SpectrumExport::write_set_options_for_chart( html_output, bg_div_id, bg_chart_options );
            html_output << bg_resize_fn << "();" << endl;
            html_output << "</script>" << endl;
          }//bg diagnostic chart
        }//if( background_auto_search_peaks present )

        html_output << "</div>" << endl;

        html_output << "</fieldset>" << endl;
      }// HTML chart block

      // Only fully-built fragments are emitted; an exception above leaves res.html empty (the
      //  spectrum still shows up in the worst-spectra table via res.qual).
      res.html = html_output.str();
    }
    catch( const std::exception &e )
    {
      cerr << "Error processing " << src_name << ": " << e.what() << endl;
      SpectrumQuality q;
      q.src_name = src_name;
      q.detector = detector_name;
      q.cost = sm_fit_failure_penalty;
      q.status = "EXCEPTION";
      q.frac_area_missing = 1.0;
      res.qual = q;
    }
   } );//boost::asio::post( pool, [&]() - one task per spectrum
  }//for( size_t spec_index = 0; spec_index < precomputed.size(); ++spec_index )

  pool.join();

  // Merge the per-spectrum fragments in spectrum order so the HTML is deterministic, and sum the
  //  per-spectrum totals.  Failed/excepted spectra have an empty fragment and zero contributions.
  for( const PerSpectrumResult &res : results )
  {
    html_output << res.html;
    total_fg += res.fg;
    total_bg_raw += res.bg_raw;
    total_score += res.score;
    spectrum_quals.push_back( res.qual );
  }

  // Worst-fitting spectra summary - the likely-catastrophic nuclides worth inspecting.  Two views,
  //  because no single rank catches both failure modes: (1) by the GA's per-spectrum cost, which
  //  surfaces area-accuracy errors on peaks that WERE found, and (2) by fraction of expected
  //  peak-area not fit, which surfaces total-misses (a fit that finds none of its expected peaks
  //  scores ~0 cost - no reward, but also no penalty - so it is invisible to view 1).
  if( !spectrum_quals.empty() )
  {
    size_t num_failed = 0;
    for( const SpectrumQuality &q : spectrum_quals )
      num_failed += (q.status != "ok");

    // Emit one worst-spectra table.  `quals` must already be sorted worst-first.
    auto write_worst_table = [&html_output,num_failed]( const std::vector<SpectrumQuality> &quals,
                                                        const std::string &heading,
                                                        const std::string &metric_desc )
    {
      // Show the worst ~10% (at least 10), but never fewer than all the outright failures.
      size_t num_show = std::max<size_t>( 10, static_cast<size_t>( std::ceil( 0.10 * quals.size() ) ) );
      num_show = std::max( num_show, num_failed );
      num_show = std::min( num_show, quals.size() );

      html_output << "<fieldset><legend>" << heading << "</legend>" << endl;
      html_output << "<p style=\"max-width:900px;margin:5px auto;\">Showing the " << num_show
                  << " worst of " << quals.size() << " spectra, ranked by " << metric_desc << ". "
                  << num_failed << " spectra failed to fit.</p>" << endl;
      html_output << "<table class=\"TopLinesTable\"><tr>"
                  << "<th>#</th><th>Source</th><th>Detector</th><th>Cost</th><th>Status</th>"
                  << "<th>Missing peaks</th><th>Extra peaks</th><th>Frac. area missing</th></tr>" << endl;
      const std::ios_base::fmtflags saved_flags = html_output.flags();
      const std::streamsize saved_prec = html_output.precision();
      html_output << std::fixed << std::setprecision(2);
      for( size_t i = 0; i < num_show; ++i )
      {
        const SpectrumQuality &q = quals[i];
        const bool ok = (q.status == "ok");
        html_output << "<tr" << (ok ? "" : " style=\"color:darkred;\"") << ">"
                    << "<td>" << (i+1) << "</td>"
                    << "<td>" << q.src_name << "</td>"
                    << "<td>" << q.detector << "</td>"
                    << "<td>" << q.cost << "</td>"
                    << "<td>" << q.status << "</td>"
                    << "<td>" << (ok ? std::to_string(q.n_missing) : std::string("-")) << "</td>"
                    << "<td>" << (ok ? std::to_string(q.n_extra) : std::string("-")) << "</td>"
                    << "<td>" << q.frac_area_missing << "</td>"
                    << "</tr>" << endl;
      }
      html_output.flags( saved_flags );
      html_output.precision( saved_prec );
      html_output << "</table></fieldset>" << endl;
    };//write_worst_table

    // View 1: ranked by per-spectrum cost (the value the GA minimizes).
    std::stable_sort( begin(spectrum_quals), end(spectrum_quals),
      []( const SpectrumQuality &a, const SpectrumQuality &b ){ return a.cost > b.cost; } );
    write_worst_table( spectrum_quals, "Worst-fitting spectra - by GA cost",
                       "per-spectrum cost (the value the GA minimizes; higher = worse)" );

    // View 2: ranked by fraction of expected peak-area not fit (catches total-misses).
    std::stable_sort( begin(spectrum_quals), end(spectrum_quals),
      []( const SpectrumQuality &a, const SpectrumQuality &b ){ return a.frac_area_missing > b.frac_area_missing; } );
    write_worst_table( spectrum_quals, "Worst-fitting spectra - by expected-area missed",
                       "fraction of expected peak-area not fit (higher = worse; catches total-misses)" );
  }//if( !spectrum_quals.empty() )

  if( sm_do_background_fit_trial )
  {
    const double weighted_bg = sm_background_fit_penalty_weight * total_bg_raw;
    html_output << "<h2 style=\"text-align: center;\">Total Score: " << total_score
                << "</h2>" << endl;
    html_output << "<p style=\"text-align: center;\">"
                << "foreground only = " << total_fg
                << "  &nbsp;|&nbsp;  background penalty (raw) = " << total_bg_raw
                << "  &nbsp;|&nbsp;  weighted bg contribution = " << weighted_bg
                << "  &nbsp;|&nbsp;  combined = fg + 0.25 * bg = " << total_score
                << "</p>" << endl;
  }
  else
  {
    html_output << "<h2 style=\"text-align: center;\">Total Score: " << total_score << "</h2>" << endl;
  }
  html_output << "</body>" << endl;
  html_output << "</html>" << endl;
  html_output.close();
#endif
}//write_results_html_and_n42


PeakFitForNuclideConfig do_nuclide_config_ga(
  const std::vector<PrecomputedNuclideData> &precomputed,
  std::function<double( const PeakFitForNuclideConfig &, double *, double * )> ga_eval_fcn )
{
  assert( !sm_has_been_called );
  if( sm_has_been_called )
  {
    cerr << "You may only call NuclideConfig_GA::do_nuclide_config_ga(...) once per program execution." << endl;
    exit( 1 );
  }

  sm_has_been_called = true;

  assert( !!ga_eval_fcn );
  if( !ga_eval_fcn )
    throw runtime_error( "Invalid eval function passed in." );

  ns_ga_eval_fcn = ga_eval_fcn;
  sm_precomputed_ptr = &precomputed;

  sm_output_file.open( PeakFitImprove::sm_output_file_prefix + "nuclide_config_ga_results.txt" );
  if( sm_do_background_fit_trial )
  {
    sm_output_file << "step\tcost_avg\tcost_best\tfg_best\tbg_best_raw"
                   << "\tsolution_best\tstatus\n";
  }
  else
  {
    sm_output_file << "step\tcost_avg\tcost_best\tsolution_best\tstatus\n";
  }

  EA::Chronometer timer;
  timer.tic();

  GA_Type ga_obj;
  ga_obj.problem_mode = EA::GA_MODE::SOGA;
  ga_obj.multi_threading = true;
  ga_obj.idle_delay_us = 1;
  ga_obj.dynamic_threading = (PeakFitImprove::sm_ga_population > PeakFitImprove::sm_num_optimization_threads);
  ga_obj.verbose = false;
  ga_obj.population = static_cast<unsigned int>( PeakFitImprove::sm_ga_population );
  ga_obj.N_threads = static_cast<int>( PeakFitImprove::sm_num_optimization_threads );
  ga_obj.generation_max = static_cast<int>( PeakFitImprove::sm_ga_generation_max );
  ga_obj.calculate_SO_total_fitness = calculate_SO_total_fitness;
  ga_obj.init_genes = init_genes;
  ga_obj.eval_solution = eval_solution;
  ga_obj.mutate = mutate;
  ga_obj.crossover = crossover;
  ga_obj.SO_report_generation = SO_report_generation;
  ga_obj.crossover_fraction = PeakFitImprove::sm_ga_crossover_fraction;
  ga_obj.mutation_rate = PeakFitImprove::sm_ga_mutation_rate;
  ga_obj.best_stall_max = static_cast<int>( PeakFitImprove::sm_ga_best_stall_max );
  ga_obj.elite_count = static_cast<int>( PeakFitImprove::sm_ga_elite_count );

  // Resume: seed the initial population from a checkpoint (--resume).  Genes only; openGA re-evaluates
  //  each via eval_solution and fills any remainder up to `population` with fresh random genes.
  if( !sm_resume_path.empty() )
  {
    std::ifstream resume_in( sm_resume_path.c_str() );
    if( !resume_in )
    {
      cerr << "Error: could not open --resume checkpoint '" << sm_resume_path << "'" << endl;
      exit( 1 );
    }

    std::vector<NuclideConfigSolution> seeded;
    size_t bad_lines = 0;
    std::string line;
    while( std::getline( resume_in, line ) )
    {
      if( !line.empty() && (line[0] == '#') )
      {
        // Refuse to resume into an incompatible run (different det-type, dataset, objective, ...).
        const std::string tag = "# options: ";
        if( (line.size() >= tag.size()) && (line.compare( 0, tag.size(), tag ) == 0)
            && (line.substr( tag.size() ) != sm_checkpoint_options_summary) )
        {
          cerr << "Error: --resume checkpoint is incompatible with this run.\n"
               << "  checkpoint options: " << line.substr( tag.size() ) << "\n"
               << "  this run's options:  " << sm_checkpoint_options_summary << "\n"
               << "Resume with matching options, or start a fresh run." << endl;
          exit( 1 );
        }
        continue;
      }//if( comment / header line )

      if( line.empty() )
        continue;

      NuclideConfigSolution s;
      if( NuclideConfigSolution::from_string( line, "\t", s ) )
        seeded.push_back( s );
      else
        ++bad_lines;
    }//while( read each checkpoint line )

    if( seeded.empty() )
    {
      cerr << "Error: no usable solutions parsed from --resume checkpoint '" << sm_resume_path << "'" << endl;
      exit( 1 );
    }

    cout << "Resuming from " << sm_resume_path << ": seeded " << seeded.size() << " solutions"
         << (bad_lines ? (" (" + std::to_string(bad_lines) + " unparseable lines skipped)") : "")
         << "." << endl;
    ga_obj.user_initial_solutions = std::move( seeded );
  }//if( !sm_resume_path.empty() )

  // Print a VERY rough estimate of per-generation wall time before the (long) solve starts, so the
  //  user can sanity-check the run scale.  Model (validated against past runs to ~1%):
  //     wall/gen ≈ evals × spectra × per_fit_seconds / effective_parallelism
  //  where evals ≈ population (gen 0) or population - elite_count (steady state), and
  //  effective_parallelism = min(outer_threads × inner_threads, hardware threads).  per_fit is a
  //  rough single-thread per-spectrum fit time keyed on detector resolution (HPGe fits are slower).
  {
    const size_t spectra = precomputed.size();
    const size_t pop = PeakFitImprove::sm_ga_population;
    const size_t elite = PeakFitImprove::sm_ga_elite_count;
    const size_t outer = std::max<size_t>( 1, PeakFitImprove::sm_num_optimization_threads );
    const size_t inner = std::max<size_t>( 1, PeakFitImprove::sm_num_threads_per_individual );
    const unsigned hw = std::max( 1u, std::thread::hardware_concurrency() );
    const size_t effective = std::min<size_t>( outer * inner, hw );
    const size_t evals_gen0 = pop;
    const size_t evals_steady = (pop > elite) ? (pop - elite) : pop;
    // Rough per-spectrum single-thread fit time (seconds); HPGe ~ measured 7.3 s, others ~ 4.4 s.
    const bool is_hpge = (sm_base_det_type == PeakFitUtils::CoarseResolutionType::High);
    const double per_fit = is_hpge ? 8.0 : 5.0;
    const double cold_factor = 1.6;  // gen 0/1 are slower: random configs hit the Ceres time cap

    auto hms = []( double s ) -> std::string {
      const long t = static_cast<long>( s + 0.5 );
      const long h = t / 3600, m = (t % 3600) / 60, sec = t % 60;
      std::ostringstream o;
      if( h ) o << h << "h ";
      if( h || m ) o << m << "m ";
      o << sec << "s";
      return o.str();
    };

    const double wall_steady = (evals_steady * (double)spectra * per_fit) / effective;
    const double wall_gen0 = cold_factor * (evals_gen0 * (double)spectra * per_fit) / effective;

    cout << "\n=== Rough runtime estimate (very approximate) ===\n"
         << "  spectra/individual: " << spectra
         << "   population: " << pop << " (elite " << elite << ")"
         << "   threads: " << outer << " outer x " << inner << " inner"
         << "   hw-threads: " << hw << "\n"
         << "  effective parallelism: " << effective << " concurrent fits"
         << (((outer*inner) > hw) ? "  [WARNING: outer x inner exceeds hw-threads -> oversubscribed]" : "")
         << "\n"
         << "  assumed per-spectrum fit: ~" << per_fit << " s (" << (is_hpge ? "HPGe" : "low/med-res") << ", rough)\n"
         << "  est. wall/generation: ~" << hms(wall_steady) << " steady-state"
         << "   (~" << hms(wall_gen0) << " for gen 0-1, slower)\n"
         << "  est. CPU-time/generation: ~" << std::fixed << std::setprecision(1)
         << (wall_steady * effective / 3600.0) << " CPU-hours\n"
         << "  NOTE: refined automatically once real per-generation Exe_time is reported below.\n"
         << "=================================================\n" << endl;
  }

  const EA::StopReason stop_reason = ga_obj.solve();

  cout << "The problem is optimized in " << timer.toc() << " seconds." << endl;
  cout << "Stop reason was: " << ga_obj.stop_reason_to_string( stop_reason ) << endl;

  sm_output_file.close();
  sm_precomputed_ptr = nullptr;

  // Print the best config in copy-pasteable format
  const PeakFitForNuclideConfig best_config = genes_to_settings( sm_best_genes );
  cout << "\n========================================" << endl;
  cout << "Best NuclideConfig GA result (cost=" << sm_best_total_cost << "):" << endl;
  cout << "Genes:\n\t" << sm_best_genes.to_string( "\n\t" ) << endl;
  cout << "========================================\n" << endl;

  // Write final HTML with best solution
  try
  {
    const string final_html = PeakFitImprove::sm_output_file_prefix + "nuclide_config_ga_final.html";
    const string final_n42_dir = PeakFitImprove::sm_output_file_prefix + "output_n42_nuclide_config_ga_final";
    write_results_html_and_n42( precomputed, best_config, sm_best_genes, sm_background_mode, final_html, final_n42_dir );
    cout << "Wrote final HTML results to " << final_html << endl;
  }
  catch( const std::exception &e )
  {
    cerr << "Warning: failed to write final HTML results: " << e.what() << endl;
  }

  return best_config;
}//do_nuclide_config_ga

}//namespace NuclideConfig_GA
