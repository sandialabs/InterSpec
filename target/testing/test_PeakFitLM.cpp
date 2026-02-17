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

#include <map>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>

#if( defined(WIN32) )
#undef min
#undef max
#endif

#define BOOST_TEST_MODULE PeakFitLM_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

std::string g_test_file_dir;

void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;

  s_have_set = true;

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  string datadir;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );

    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }

  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
      {
        datadir = d;
        break;
      }
    }
  }

  if( g_test_file_dir.empty() )
  {
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "../../../test_data",
                            "../../target/testing/test_data", "../../../target/testing/test_data" } )
    {
      if( SpecUtils::is_directory( d ) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }

  const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
    "Could not find 'sandia.decay.xml' at '" << sandia_decay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide( "U238" ), "SandiaDecayDataBase empty?" );
}//void set_data_dir()


/// Describes a test spectrum file and its expected characteristics.
struct TestFileConfig
{
  /// Filename relative to the `analysis_tests` directory (sibling of test_data)
  string filename;

  /// Expected detector resolution type
  PeakFitUtils::CoarseResolutionType resolution_type;

  /// Minimum number of peaks expected in the file
  size_t min_expected_peaks;

  /// Energy ranges to exclude from tight tolerance checks.  These are typically multi-peak ROIs
  /// where peaks are very close together and refitting can cause large deviations or peak loss.
  /// Peaks whose reference mean falls in any of these ranges are still printed but not included
  /// in the deviation statistics.
  vector<pair<double,double>> exclude_energy_ranges;
};


/// Returns true if `energy` falls within any of the exclude ranges.
bool is_excluded_energy( const double energy, const vector<pair<double,double>> &ranges )
{
  for( const pair<double,double> &r : ranges )
  {
    if( energy >= r.first && energy <= r.second )
      return true;
  }
  return false;
}


/// The list of test files to run all tests over.
static const vector<TestFileConfig> g_test_files = {
  //{ "u238_u235_1kg.n42",
  //  PeakFitUtils::CoarseResolutionType::High, 5,
  //  { {178.0, 190.0} } // 2-peak ROI at ~183/185 keV: refitting makes one peak vanish
  //},
  { "AEGIS_Eu152_surface_contamination.n42_20230622T113239.276178.n42",
    PeakFitUtils::CoarseResolutionType::High, 5,
    {}
    //{ {44.0, 48.0},      // 2-peak ROI at ~45/46 keV: very close peaks
    //  {1083.0, 1092.0} } // 2-peak ROI at ~1085/1089 keV
  },
};


/// Loaded spectrum data for a test file.  Cached to avoid reloading.
struct SpectrumTestData
{
  shared_ptr<SpecMeas> meas;
  shared_ptr<const SpecUtils::Measurement> foreground;
  vector<shared_ptr<const PeakDef>> original_peaks;
  string label; // short name for log output

  void load( const string &filename )
  {
    set_data_dir();

    const string spec_path = SpecUtils::append_path( g_test_file_dir,
      "../analysis_tests/" + filename );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ),
      "Could not find test file: '" << spec_path << "'" );

    meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( meas->load_N42_file( spec_path ),
      "Failed to open required file: '" << spec_path << "'" );

    BOOST_REQUIRE( meas->num_measurements() >= 1 );

    // Get the foreground measurement
    foreground = meas->measurements()[0];
    BOOST_REQUIRE( foreground );
    BOOST_REQUIRE( foreground->num_gamma_channels() > 64 );

    // Get peaks
    const set<int> sample_nums = { foreground->sample_number() };
    shared_ptr<const deque<shared_ptr<const PeakDef>>> allpeaks = meas->peaks( sample_nums );
    BOOST_REQUIRE_MESSAGE( allpeaks && !allpeaks->empty(),
      "No peaks in test spectrum '" << filename << "'" );

    original_peaks.assign( allpeaks->begin(), allpeaks->end() );

    // Extract a short label from the filename (up to first '.' or 30 chars)
    const size_t dot = filename.find( '.' );
    label = filename.substr( 0, std::min( dot, size_t(30) ) );
  }
};//struct SpectrumTestData


/// Cache of loaded test data, keyed by filename.
static map<string, SpectrumTestData> g_loaded_data;

SpectrumTestData &get_test_data( const TestFileConfig &config )
{
  auto it = g_loaded_data.find( config.filename );
  if( it != g_loaded_data.end() )
    return it->second;

  SpectrumTestData &data = g_loaded_data[config.filename];
  data.load( config.filename );

  BOOST_REQUIRE_MESSAGE( data.original_peaks.size() >= config.min_expected_peaks,
    "File '" << config.filename << "' has " << data.original_peaks.size()
    << " peaks, expected at least " << config.min_expected_peaks );

  cout << "Loaded " << data.original_peaks.size() << " peaks from '"
       << config.filename << "'." << endl;

  return data;
}//SpectrumTestData &get_test_data(...)


/// Helper to group peaks by shared continuum (ROI).
map<shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>>
group_peaks_by_roi( const vector<shared_ptr<const PeakDef>> &peaks )
{
  map<shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>> result;
  for( const shared_ptr<const PeakDef> &p : peaks )
    result[p->continuum()].push_back( p );
  return result;
}


/// 1-to-1 matching of reference peaks to test peaks using greedy nearest-neighbor.
/// Returns a vector of pairs (ref_index, test_index) for each matched pair.
/// Unmatched reference peaks are returned in `unmatched_ref`.
vector<pair<size_t,size_t>> match_peaks_1to1(
  const vector<shared_ptr<const PeakDef>> &ref_peaks,
  const vector<shared_ptr<const PeakDef>> &test_peaks,
  vector<size_t> &unmatched_ref,
  const double max_sigma_dist = 5.0 )
{
  vector<pair<size_t,size_t>> matches;
  vector<bool> test_used( test_peaks.size(), false );
  unmatched_ref.clear();

  // Build list of ref Gaussian peak indices sorted by amplitude (match strongest first)
  vector<size_t> ref_indices;
  for( size_t i = 0; i < ref_peaks.size(); ++i )
  {
    if( ref_peaks[i]->gausPeak() )
      ref_indices.push_back( i );
  }
  std::sort( ref_indices.begin(), ref_indices.end(), [&]( size_t a, size_t b ){
    return ref_peaks[a]->amplitude() > ref_peaks[b]->amplitude();
  } );

  for( const size_t ri : ref_indices )
  {
    const shared_ptr<const PeakDef> &ref = ref_peaks[ri];
    double best_diff = std::numeric_limits<double>::max();
    size_t best_ti = test_peaks.size();

    for( size_t ti = 0; ti < test_peaks.size(); ++ti )
    {
      if( test_used[ti] )
        continue;
      const double diff = fabs( test_peaks[ti]->mean() - ref->mean() );
      if( diff < best_diff )
      {
        best_diff = diff;
        best_ti = ti;
      }
    }

    if( (best_ti < test_peaks.size()) && (best_diff <= max_sigma_dist * ref->sigma()) )
    {
      matches.push_back( { ri, best_ti } );
      test_used[best_ti] = true;
    }
    else
    {
      unmatched_ref.push_back( ri );
    }
  }

  return matches;
}


struct DeviationStats
{
  double max_area_pct = 0.0;
  double max_fwhm_pct = 0.0;
  double max_mean_pct = 0.0;
  double max_area_peak_mean = 0.0;
  double max_fwhm_peak_mean = 0.0;
  size_t num_matched = 0;
  size_t num_unmatched = 0;

  /// All per-peak area deviations, sorted ascending
  vector<double> sorted_area_pcts;
  /// All per-peak FWHM deviations, sorted ascending
  vector<double> sorted_fwhm_pcts;

  /// Returns the Nth-from-worst area deviation (0 = worst, 1 = second worst, etc.)
  /// Returns 0.0 if n >= num_matched.
  double nth_worst_area( const size_t n ) const
  {
    if( n >= sorted_area_pcts.size() )
      return 0.0;
    return sorted_area_pcts[sorted_area_pcts.size() - 1 - n];
  }
  double nth_worst_fwhm( const size_t n ) const
  {
    if( n >= sorted_fwhm_pcts.size() )
      return 0.0;
    return sorted_fwhm_pcts[sorted_fwhm_pcts.size() - 1 - n];
  }
};

DeviationStats compute_deviations(
  const vector<shared_ptr<const PeakDef>> &ref_peaks,
  const vector<shared_ptr<const PeakDef>> &test_peaks,
  const string &label,
  const vector<pair<double,double>> &exclude_ranges = {} )
{
  DeviationStats stats;
  vector<size_t> unmatched;
  const vector<pair<size_t,size_t>> matches = match_peaks_1to1( ref_peaks, test_peaks, unmatched );

  size_t num_excluded = 0;

  cout << label << " per-peak deviations:" << endl;
  for( const auto &m : matches )
  {
    const shared_ptr<const PeakDef> &ref = ref_peaks[m.first];
    const shared_ptr<const PeakDef> &fit = test_peaks[m.second];

    const double area_pct = 100.0 * fabs( fit->amplitude() - ref->amplitude() )
                            / std::max( 1.0, fabs( ref->amplitude() ) );
    const double fwhm_pct = 100.0 * fabs( fit->fwhm() - ref->fwhm() )
                            / std::max( 0.01, ref->fwhm() );
    const double mean_pct = 100.0 * fabs( fit->mean() - ref->mean() )
                            / std::max( 0.01, ref->mean() );

    const bool excluded = is_excluded_energy( ref->mean(), exclude_ranges );

    cout << "  " << ref->mean() << " keV: area=" << area_pct
         << "%, fwhm=" << fwhm_pct << "%, mean=" << mean_pct << "%"
         << ", Chi2/dof(ref->fit): " << ref->chi2dof() << " --> " << fit->chi2dof()
         << (excluded ? "  [EXCLUDED from stats]" : "")
         << endl;

    if( excluded )
    {
      ++num_excluded;
      continue;
    }

    stats.sorted_area_pcts.push_back( area_pct );
    stats.sorted_fwhm_pcts.push_back( fwhm_pct );

    if( area_pct > stats.max_area_pct )
    {
      stats.max_area_pct = area_pct;
      stats.max_area_peak_mean = ref->mean();
    }
    if( fwhm_pct > stats.max_fwhm_pct )
    {
      stats.max_fwhm_pct = fwhm_pct;
      stats.max_fwhm_peak_mean = ref->mean();
    }
    if( mean_pct > stats.max_mean_pct )
      stats.max_mean_pct = mean_pct;
  }//for( const auto &m : matches )

  for( const size_t ui : unmatched )
  {
    const bool excluded = is_excluded_energy( ref_peaks[ui]->mean(), exclude_ranges );
    if( excluded )
      ++num_excluded;
    cout << "  " << ref_peaks[ui]->mean() << " keV: UNMATCHED"
         << (excluded ? "  [EXCLUDED]" : "") << endl;
  }

  stats.num_matched = stats.sorted_area_pcts.size();
  stats.num_unmatched = unmatched.size();

  std::sort( stats.sorted_area_pcts.begin(), stats.sorted_area_pcts.end() );
  std::sort( stats.sorted_fwhm_pcts.begin(), stats.sorted_fwhm_pcts.end() );

  cout << label << " summary (matched=" << stats.num_matched
       << ", unmatched=" << stats.num_unmatched
       << ", excluded=" << num_excluded << "):"
       << "\n  max area%=" << stats.max_area_pct << " (at " << stats.max_area_peak_mean << " keV)"
       << "\n  max fwhm%=" << stats.max_fwhm_pct << " (at " << stats.max_fwhm_peak_mean << " keV)"
       << "\n  max mean%=" << stats.max_mean_pct;
  if( stats.sorted_area_pcts.size() >= 3 )
    cout << "\n  3rd worst area%=" << stats.nth_worst_area( 2 )
         << ", 3rd worst fwhm%=" << stats.nth_worst_fwhm( 2 );
  cout << endl;

  return stats;
}


BOOST_AUTO_TEST_CASE( test_refit_all_peaks_no_skew )
{
  for( const TestFileConfig &config : g_test_files )
  {
    SpectrumTestData &td = get_test_data( config );
    const string prefix = "[" + td.label + " NoSkew] ";

    cout << prefix << "Refitting " << td.original_peaks.size() << " peaks..." << endl;

    const PeakFitLM::FitPeaksResults result = PeakFitLM::fit_peaks_in_spectrum_LM(
      td.original_peaks, td.foreground, 0.0, 0.0,
      config.resolution_type,
      std::nullopt, 0 );

    BOOST_REQUIRE_MESSAGE( result.status == PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success,
      prefix << "fit failed: " << result.error_message );

    const size_t num_gauss_orig = static_cast<size_t>(
      count_if( td.original_peaks.begin(), td.original_peaks.end(),
                []( const shared_ptr<const PeakDef> &p ){ return p->gausPeak(); } ) );

    cout << prefix << "Original Gaussian peaks: " << num_gauss_orig
         << ", fit peaks: " << result.fit_peaks.size()
         << ", lost peaks: " << result.lost_peaks.size() << endl;

    // Allow losing up to 2 peaks (multi-peak ROI significance culling)
    BOOST_CHECK_MESSAGE( result.fit_peaks.size() >= num_gauss_orig - 2,
      prefix << "Expected >= " << (num_gauss_orig - 2) << " fit peaks, got " << result.fit_peaks.size() );

    const DeviationStats stats = compute_deviations( td.original_peaks, result.fit_peaks,
      prefix + "refit", config.exclude_energy_ranges );

    // Peaks in multi-peak ROIs can shift significantly (esp. when one peak is borderline-significant),
    // so we check the 3rd-worst deviation against a tight tolerance, and the overall max against a
    // wider one.  Most peaks (single-peak ROIs) should refit almost exactly.
    BOOST_CHECK_MESSAGE( stats.nth_worst_area( 2 ) < 1.0,
      prefix << "3rd-worst area deviation " << stats.nth_worst_area( 2 ) << "% exceeds 1%" );
    BOOST_CHECK_MESSAGE( stats.nth_worst_fwhm( 2 ) < 5.0,
      prefix << "3rd-worst FWHM deviation " << stats.nth_worst_fwhm( 2 ) << "% exceeds 5%" );
    BOOST_CHECK_MESSAGE( stats.max_area_pct < 75.0,
      prefix << "Max area deviation " << stats.max_area_pct << "% exceeds 75% (even for outlier peaks)" );
    BOOST_CHECK_MESSAGE( stats.num_unmatched <= 2,
      prefix << "More than 2 peaks unmatched (" << stats.num_unmatched << ")" );

    // Verify values are not exactly the same (fitting actually ran)
    BOOST_CHECK_MESSAGE( stats.num_matched > 0 && stats.max_area_pct > 0.001,
      prefix << "All peak values are identical to input - fitting may not have run" );

    // No skew relation expected for nullopt skew
    BOOST_CHECK( !result.skew_relation.has_value() );
  }
}//test_refit_all_peaks_no_skew


BOOST_AUTO_TEST_CASE( test_refit_per_roi_matches_all_at_once )
{
  for( const TestFileConfig &config : g_test_files )
  {
    SpectrumTestData &td = get_test_data( config );
    const string prefix = "[" + td.label + " PerROI] ";

    // First: fit all peaks at once (baseline)
    const PeakFitLM::FitPeaksResults all_at_once = PeakFitLM::fit_peaks_in_spectrum_LM(
      td.original_peaks, td.foreground, 0.0, 0.0,
      config.resolution_type,
      std::nullopt, 0 );

    BOOST_REQUIRE( all_at_once.status == PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success );

    // Second: fit each ROI independently
    const auto roi_groups = group_peaks_by_roi( td.original_peaks );
    vector<shared_ptr<const PeakDef>> per_roi_peaks;

    for( const auto &kv : roi_groups )
    {
      const PeakFitLM::FitPeaksResults roi_result = PeakFitLM::fit_peaks_in_spectrum_LM(
        kv.second, td.foreground, 0.0, 0.0,
        config.resolution_type,
        std::nullopt, 0 );

      BOOST_REQUIRE_MESSAGE( roi_result.status == PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success,
        prefix << "Per-ROI fit failed: " << roi_result.error_message );

      per_roi_peaks.insert( end(per_roi_peaks), begin(roi_result.fit_peaks), end(roi_result.fit_peaks) );
    }

    BOOST_CHECK_MESSAGE( per_roi_peaks.size() == all_at_once.fit_peaks.size(),
      prefix << "Per-ROI fit returned " << per_roi_peaks.size() << " peaks vs "
      << all_at_once.fit_peaks.size() << " from all-at-once" );

    // Compare: per-ROI results should be very close to all-at-once (since with null skew,
    //  fit_peaks_in_spectrum_LM fits each ROI independently anyway)
    const DeviationStats stats = compute_deviations( all_at_once.fit_peaks, per_roi_peaks,
      prefix + "vs all-at-once" );

    const double tolerance_pct = 0.1;
    BOOST_CHECK_MESSAGE( stats.max_area_pct < tolerance_pct,
      prefix << "Max area deviation " << stats.max_area_pct << "% exceeds " << tolerance_pct << "%" );
    BOOST_CHECK_MESSAGE( stats.max_fwhm_pct < tolerance_pct,
      prefix << "Max FWHM deviation " << stats.max_fwhm_pct << "% exceeds " << tolerance_pct << "%" );
    BOOST_CHECK_MESSAGE( stats.max_mean_pct < tolerance_pct,
      prefix << "Max mean deviation " << stats.max_mean_pct << "% exceeds " << tolerance_pct << "%" );
  }
}//test_refit_per_roi_matches_all_at_once


void test_skew_type_helper( const PeakDef::SkewType skew_type,
                            const string &skew_name,
                            const double area_tolerance_pct,
                            const double fwhm_tolerance_pct,
                            const Wt::WFlags<PeakFitLM::PeakFitLMOptions> extra_options = 0 )
{
  for( const TestFileConfig &config : g_test_files )
  {
    SpectrumTestData &td = get_test_data( config );
    const string prefix = "[" + td.label + " " + skew_name + "] ";

    cout << prefix << "Fitting " << td.original_peaks.size() << " peaks..." << endl;

    PeakFitLM::FitPeaksResults result;
    BOOST_REQUIRE_NO_THROW(
      result = PeakFitLM::fit_peaks_in_spectrum_LM(
        td.original_peaks, td.foreground, 0.0, 0.0,
        config.resolution_type,
        skew_type, extra_options )
    );

    cout << prefix << "fit completed." << endl;

    BOOST_REQUIRE_MESSAGE( result.status == PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success,
      prefix << "fit failed: " << result.error_message );

    const size_t num_gauss_orig = static_cast<size_t>(
      count_if( td.original_peaks.begin(), td.original_peaks.end(),
                []( const shared_ptr<const PeakDef> &p ){ return p->gausPeak(); } ) );

    cout << prefix << "original Gaussian peaks: " << num_gauss_orig
         << ", fit peaks: " << result.fit_peaks.size()
         << ", lost peaks: " << result.lost_peaks.size() << endl;

    // All (or nearly all) original peaks should be in the output
    BOOST_CHECK_MESSAGE( result.fit_peaks.size() >= num_gauss_orig - 2,
      prefix << "Expected >= " << (num_gauss_orig - 2) << " fit peaks, got " << result.fit_peaks.size() );

    const DeviationStats stats = compute_deviations( td.original_peaks, result.fit_peaks,
      prefix, config.exclude_energy_ranges );

    // Multi-peak ROI peaks can shift significantly when changing skew type, so check 3rd-worst
    // deviation against the tight tolerance, and overall max against a wider one.
    BOOST_CHECK_MESSAGE( stats.nth_worst_area( 2 ) < area_tolerance_pct,
      prefix << "3rd-worst area deviation " << stats.nth_worst_area( 2 )
      << "% exceeds " << area_tolerance_pct << "%" );
    BOOST_CHECK_MESSAGE( stats.nth_worst_fwhm( 2 ) < fwhm_tolerance_pct,
      prefix << "3rd-worst FWHM deviation " << stats.nth_worst_fwhm( 2 )
      << "% exceeds " << fwhm_tolerance_pct << "%" );
    BOOST_CHECK_MESSAGE( stats.max_area_pct < 75.0,
      prefix << "Max area deviation " << stats.max_area_pct
      << "% exceeds 75% (even for outlier peaks)" );
    BOOST_CHECK_MESSAGE( stats.num_unmatched <= 2,
      prefix << "More than 2 peaks unmatched (" << stats.num_unmatched << ")" );

    // Check all fit peaks have the expected skew type (or VoigtPlusBortel exception)
    for( const shared_ptr<const PeakDef> &p : result.fit_peaks )
    {
      if( !p->gausPeak() )
        continue;
      const bool ok = (p->skewType() == skew_type)
        || (p->skewType() == PeakDef::SkewType::VoigtPlusBortel
            && (skew_type == PeakDef::SkewType::Bortel || skew_type == PeakDef::SkewType::GaussPlusBortel));
      BOOST_CHECK_MESSAGE( ok,
        prefix << "Peak at " << p->mean() << " keV has unexpected skew type "
        << static_cast<int>(p->skewType()) );
    }

    // If multi-ROI and not IndependentSkewValues, check skew_relation
    const auto roi_groups = group_peaks_by_roi( result.fit_peaks );
    const bool is_independent = extra_options.testFlag( PeakFitLM::PeakFitLMOptions::IndependentSkewValues );

    if( roi_groups.size() > 1 && !is_independent
       && (PeakDef::num_skew_parameters( skew_type ) > 0) )
    {
      BOOST_CHECK_MESSAGE( result.skew_relation.has_value(),
        prefix << "Expected skew_relation for multi-ROI fit" );

      if( result.skew_relation.has_value() )
      {
        BOOST_CHECK( result.skew_relation->skew_type == skew_type );
        BOOST_CHECK( result.skew_relation->energy_range.has_value() );

        if( result.skew_relation->energy_range.has_value() )
        {
          cout << prefix << "skew_relation energy_range: ["
               << result.skew_relation->energy_range->first << ", "
               << result.skew_relation->energy_range->second << "]" << endl;
        }
      }
    }
  }//for( const TestFileConfig &config : g_test_files )
}//test_skew_type_helper


BOOST_AUTO_TEST_CASE( test_fit_bortel_skew )
{
  // Bortel is a mild single-param skew; 3rd-worst peaks should stay < 5% area, < 7% FWHM
  test_skew_type_helper( PeakDef::SkewType::Bortel, "Bortel", 5.0, 7.0 );
}

BOOST_AUTO_TEST_CASE( test_fit_gaussexp_skew )
{
  // GaussExp adds an exponential tail - a more aggressive model change, so wider tolerances
  test_skew_type_helper( PeakDef::SkewType::GaussExp, "GaussExp", 15.0, 15.0 );
}

BOOST_AUTO_TEST_CASE( test_fit_double_bortel_skew )
{
  // DoubleBortel is similar to Bortel in effect
  test_skew_type_helper( PeakDef::SkewType::DoubleBortel, "DoubleBortel", 5.0, 5.0 );
}

BOOST_AUTO_TEST_CASE( test_bortel_independent_skew )
{
  // Independent skew per ROI allows each ROI to converge to different solutions, so wider tolerance
  test_skew_type_helper( PeakDef::SkewType::Bortel, "Bortel+Independent",
    10.0, 20.0,
    Wt::WFlags<PeakFitLM::PeakFitLMOptions>( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) );
}
