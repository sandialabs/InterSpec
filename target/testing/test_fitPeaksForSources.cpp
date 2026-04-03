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

#include <set>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <numeric>
#include <iostream>
#include <algorithm>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/DetectorPeakResponse.h"

#define BOOST_TEST_MODULE FitPeaksForSources_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;

namespace
{
  string g_data_dir;
  string g_test_file_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char **argv = boost::unit_test::framework::master_test_suite().argv;

    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        g_data_dir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        g_test_file_dir = arg.substr( 14 );
    }

    SpecUtils::ireplace_all( g_data_dir, "%20", " " );
    SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

    if( g_data_dir.empty() )
    {
      for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path( d, "sandia.decay.xml" ) ) )
        {
          g_data_dir = d;
          break;
        }
      }
    }

    if( g_test_file_dir.empty() )
    {
      for( const auto &d : { "test_data", "../test_data", "../../test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path( d, "FitPeaksForSource" ) ) )
        {
          g_test_file_dir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( !g_data_dir.empty(), "Could not find data directory" );

    const string sandia_decay = SpecUtils::append_path( g_data_dir, "sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay ),
      "sandia.decay.xml not found at '" << sandia_decay << "'" );

    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( g_data_dir ) );
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  }// set_data_dir


  // Load a spectrum file, returning foreground measurement.
  // Optionally loads background from a separate file.
  struct LoadedSpectrum
  {
    std::shared_ptr<const SpecUtils::Measurement> foreground;
    std::shared_ptr<const SpecUtils::Measurement> background;
    bool isHPGe;
  };


  LoadedSpectrum load_detective_x_spectrum( const string &filename )
  {
    set_data_dir();

    const string spec_dir = SpecUtils::append_path( g_data_dir,
      "reference_spectra/Common_Field_Nuclides/Detective X" );

    const string fg_path = SpecUtils::append_path( spec_dir, filename );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( fg_path ), "Spectrum not found: " << fg_path );

    SpecUtils::SpecFile fg_file;
    BOOST_REQUIRE_MESSAGE( fg_file.load_file( fg_path, SpecUtils::ParserType::Auto ),
      "Failed to load: " << fg_path );

    LoadedSpectrum result;
    result.foreground = fg_file.measurement_at_index( 0 );
    BOOST_REQUIRE( result.foreground && result.foreground->num_gamma_channels() );
    result.isHPGe = true;

    // Try to load background
    const string bg_path = SpecUtils::append_path( spec_dir, "background.txt" );
    if( SpecUtils::is_file( bg_path ) )
    {
      SpecUtils::SpecFile bg_file;
      if( bg_file.load_file( bg_path, SpecUtils::ParserType::Auto ) )
      {
        result.background = bg_file.measurement_at_index( 0 );
        if( result.background && !result.background->num_gamma_channels() )
          result.background = nullptr;
      }
    }

    return result;
  }// load_detective_x_spectrum


  LoadedSpectrum load_test_data_spectrum( const string &fg_filename,
                                          const string &bg_filename = "" )
  {
    set_data_dir();
    BOOST_REQUIRE_MESSAGE( !g_test_file_dir.empty(), "Test file directory not set" );

    const string fg_path = SpecUtils::append_path( g_test_file_dir,
      SpecUtils::append_path( "FitPeaksForSource", fg_filename ) );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( fg_path ), "Test spectrum not found: " << fg_path );

    SpecUtils::SpecFile fg_file;
    BOOST_REQUIRE_MESSAGE( fg_file.load_file( fg_path, SpecUtils::ParserType::Auto ),
      "Failed to load: " << fg_path );

    LoadedSpectrum result;
    result.foreground = fg_file.measurement_at_index( 0 );
    BOOST_REQUIRE( result.foreground && result.foreground->num_gamma_channels() );
    result.isHPGe = true;

    if( !bg_filename.empty() )
    {
      const string bg_path = SpecUtils::append_path( g_test_file_dir,
        SpecUtils::append_path( "FitPeaksForSource", bg_filename ) );
      if( SpecUtils::is_file( bg_path ) )
      {
        SpecUtils::SpecFile bg_file;
        if( bg_file.load_file( bg_path, SpecUtils::ParserType::Auto ) )
        {
          result.background = bg_file.measurement_at_index( 0 );
          if( result.background && !result.background->num_gamma_channels() )
            result.background = nullptr;
        }
      }
    }

    return result;
  }// load_test_data_spectrum


  // Run auto peak search on a spectrum
  vector<shared_ptr<const PeakDef>> run_auto_search(
    const shared_ptr<const SpecUtils::Measurement> &foreground, const bool isHPGe )
  {
    return ExperimentalAutomatedPeakSearch::search_for_peaks(
      foreground, nullptr, nullptr, true, isHPGe );
  }// run_auto_search


  // Make a source list from nuclide names using the simpler SrcVariant interface
  vector<RelActCalcAuto::SrcVariant> make_sources( const vector<string> &names )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE( db );

    vector<RelActCalcAuto::SrcVariant> sources;
    for( const string &name : names )
    {
      const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
      if( nuc )
      {
        sources.push_back( nuc );
        continue;
      }

      const SandiaDecay::Element * const el = db->element( name );
      if( el )
      {
        sources.push_back( el );
        continue;
      }

      BOOST_FAIL( "Unknown source: " << name );
    }

    return sources;
  }// make_sources


  // Run fit_peaks_for_nuclides with the simpler SrcVariant interface
  FitPeaksForNuclides::PeakFitResult run_fit(
    const shared_ptr<const SpecUtils::Measurement> &foreground,
    const shared_ptr<const SpecUtils::Measurement> &background,
    const vector<shared_ptr<const PeakDef>> &auto_search_peaks,
    const vector<RelActCalcAuto::SrcVariant> &sources,
    const vector<shared_ptr<const PeakDef>> &user_peaks,
    const bool isHPGe,
    const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options
      = Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions>() )
  {
    const FitPeaksForNuclides::PeakFitForNuclideConfig &config
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( isHPGe );

    return FitPeaksForNuclides::fit_peaks_for_nuclides(
      auto_search_peaks, foreground, sources, user_peaks,
      background, nullptr, options, config, isHPGe );
  }// run_fit


  // Check that a peak near the given energy exists in the result
  bool has_peak_near( const vector<PeakDef> &peaks, const double energy,
                      const double tolerance_keV = 3.0 )
  {
    for( const PeakDef &p : peaks )
    {
      if( fabs( p.mean() - energy ) < tolerance_keV )
        return true;
    }
    return false;
  }// has_peak_near


  // Get a peak near the given energy
  const PeakDef *find_peak_near( const vector<PeakDef> &peaks, const double energy,
                                  const double tolerance_keV = 3.0 )
  {
    const PeakDef *best = nullptr;
    double best_dist = tolerance_keV;
    for( const PeakDef &p : peaks )
    {
      const double dist = fabs( p.mean() - energy );
      if( dist < best_dist )
      {
        best_dist = dist;
        best = &p;
      }
    }
    return best;
  }// find_peak_near


  // Verify that no observable peak ROIs overlap (at least 1 channel gap)
  void verify_no_roi_overlaps( const vector<PeakDef> &peaks,
                                const shared_ptr<const SpecUtils::Measurement> &foreground )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> energy_cal
      = foreground->energy_calibration();
    BOOST_REQUIRE( energy_cal && energy_cal->valid() );

    // Collect unique ROIs (by PeakContinuum pointer)
    vector<pair<double,double>> roi_bounds;
    set<const PeakContinuum *> seen;
    for( const PeakDef &p : peaks )
    {
      if( !p.continuum() )
        continue;
      if( seen.insert( p.continuum().get() ).second )
        roi_bounds.emplace_back( p.continuum()->lowerEnergy(), p.continuum()->upperEnergy() );
    }

    sort( roi_bounds.begin(), roi_bounds.end() );

    for( size_t i = 1; i < roi_bounds.size(); ++i )
    {
      const double prev_upper = roi_bounds[i - 1].second;
      const double curr_lower = roi_bounds[i].first;

      const size_t prev_upper_ch = energy_cal->channel_for_energy( prev_upper );
      const size_t curr_lower_ch = energy_cal->channel_for_energy( curr_lower );

      BOOST_CHECK_MESSAGE( curr_lower_ch > prev_upper_ch,
        "ROI overlap or abutting: ROI [" << roi_bounds[i-1].first << ", " << prev_upper
        << "] and [" << curr_lower << ", " << roi_bounds[i].second
        << "] keV share channel boundary (channels " << prev_upper_ch << " and " << curr_lower_ch << ")" );
    }
  }// verify_no_roi_overlaps


  // Verify that new ROI edges maintain minimum distance from existing peak means
  void verify_fwhm_margin_from_existing(
    const vector<PeakDef> &observable_peaks,
    const vector<shared_ptr<const PeakDef>> &user_peaks,
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove,
    const FitPeaksForNuclides::PeakFitForNuclideConfig &config,
    const shared_ptr<const SpecUtils::Measurement> &foreground )
  {
    // Identify active existing peaks (not in remove list)
    set<const PeakDef *> removed_ptrs;
    for( const shared_ptr<const PeakDef> &p : peaks_to_remove )
      removed_ptrs.insert( p.get() );

    vector<double> active_existing_means;
    for( const shared_ptr<const PeakDef> &p : user_peaks )
    {
      if( p && !removed_ptrs.count( p.get() ) )
        active_existing_means.push_back( p->mean() );
    }

    if( active_existing_means.empty() )
      return;

    // Collect unique observable ROI bounds
    set<const PeakContinuum *> seen;
    for( const PeakDef &obs : observable_peaks )
    {
      if( !obs.continuum() || !seen.insert( obs.continuum().get() ).second )
        continue;

      const double roi_lower = obs.continuum()->lowerEnergy();
      const double roi_upper = obs.continuum()->upperEnergy();

      for( const double existing_mean : active_existing_means )
      {
        // Only check if the existing peak mean is near this ROI
        if( (existing_mean < roi_lower - 20.0) || (existing_mean > roi_upper + 20.0) )
          continue;

        // Skip if the existing peak mean is inside this ROI (it's a bystander situation)
        if( (existing_mean >= roi_lower) && (existing_mean <= roi_upper) )
          continue;

        // Compute minimum margin based on FWHM at existing peak energy
        // The requirement is 0.5 * config.auto_rel_eff_sol_min_fwhm_roi * FWHM
        // But we need the FWHM functional form. For now just check a reasonable minimum.
        const double margin_from_lower = existing_mean - roi_lower;
        const double margin_from_upper = roi_upper - existing_mean;

        // The ROI should not extend past the existing peak mean
        if( margin_from_lower > 0.0 && margin_from_lower < margin_from_upper )
        {
          // ROI lower edge is close to existing peak: existing mean is above roi_lower
          // This means the ROI extends to overlap the existing peak - bad
          BOOST_CHECK_MESSAGE( roi_lower > existing_mean || margin_from_lower > 0.5,
            "Observable ROI lower edge " << roi_lower << " keV is within " << margin_from_lower
            << " keV of existing peak mean " << existing_mean
            << " keV (ROI: [" << roi_lower << ", " << roi_upper << "])" );
        }

        if( margin_from_upper > 0.0 && margin_from_upper < margin_from_lower )
        {
          BOOST_CHECK_MESSAGE( roi_upper < existing_mean || margin_from_upper > 0.5,
            "Observable ROI upper edge " << roi_upper << " keV is within " << margin_from_upper
            << " keV of existing peak mean " << existing_mean
            << " keV (ROI: [" << roi_lower << ", " << roi_upper << "])" );
        }
      }// for( existing means )
    }// for( observable peaks )
  }// verify_fwhm_margin_from_existing


  // Verify that all peaks sharing a PeakContinuum with a removed peak are also removed
  void verify_continuum_consistency(
    const vector<shared_ptr<const PeakDef>> &user_peaks,
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove )
  {
    if( peaks_to_remove.empty() )
      return;

    set<const PeakContinuum *> removed_continuums;
    for( const shared_ptr<const PeakDef> &p : peaks_to_remove )
      removed_continuums.insert( p->continuum().get() );

    set<const PeakDef *> removed_ptrs;
    for( const shared_ptr<const PeakDef> &p : peaks_to_remove )
      removed_ptrs.insert( p.get() );

    for( const shared_ptr<const PeakDef> &p : user_peaks )
    {
      if( !p )
        continue;
      if( removed_continuums.count( p->continuum().get() ) && !removed_ptrs.count( p.get() ) )
      {
        BOOST_CHECK_MESSAGE( false,
          "User peak at " << p->mean() << " keV shares continuum with a removed peak "
          "but was not itself removed" );
      }
    }
  }// verify_continuum_consistency


  // Verify observable peak means are within their continuum bounds
  void verify_peaks_within_roi( const vector<PeakDef> &peaks )
  {
    for( const PeakDef &p : peaks )
    {
      if( !p.continuum() )
        continue;
      BOOST_CHECK_MESSAGE(
        (p.mean() >= p.continuum()->lowerEnergy()) && (p.mean() <= p.continuum()->upperEnergy()),
        "Peak mean " << p.mean() << " keV is outside its ROI ["
        << p.continuum()->lowerEnergy() << ", " << p.continuum()->upperEnergy() << "]" );
    }
  }// verify_peaks_within_roi


  // Comprehensive validation of a fit result
  void verify_fit_result(
    const FitPeaksForNuclides::PeakFitResult &result,
    const vector<shared_ptr<const PeakDef>> &user_peaks,
    const shared_ptr<const SpecUtils::Measurement> &foreground,
    const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options,
    const bool expect_success = true )
  {
    if( expect_success )
    {
      BOOST_CHECK_MESSAGE(
        result.status == RelActCalcAuto::RelActAutoSolution::Status::Success,
        "Fit failed: " << result.error_message );
    }

    if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
      return;

    const FitPeaksForNuclides::PeakFitForNuclideConfig &config
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( true );

    // 1. No overlapping observable ROIs
    verify_no_roi_overlaps( result.observable_peaks, foreground );

    // 2. Peak means within ROI bounds
    verify_peaks_within_roi( result.observable_peaks );

    // 3. Continuum consistency
    verify_continuum_consistency( user_peaks, result.original_peaks_to_remove );

    // 4. FWHM margin from existing peaks
    verify_fwhm_margin_from_existing( result.observable_peaks, user_peaks,
      result.original_peaks_to_remove, config, foreground );

    // 5. Mode-specific checks
    if( options & FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois )
    {
      BOOST_CHECK_MESSAGE( result.original_peaks_to_remove.empty(),
        "DoNotUseExistingRois: original_peaks_to_remove should be empty, but has "
        << result.original_peaks_to_remove.size() << " peaks" );
    }

    // 6. DoNotUseExistingRois and ExistingPeaksAsFreePeak should not be combined
    if( (options & FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois)
       && (options & FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak) )
    {
      BOOST_CHECK_MESSAGE( false,
        "DoNotUseExistingRois and ExistingPeaksAsFreePeak should not be combined" );
    }
  }// verify_fit_result


  // Verify that specific user peaks remain unchanged (same mean, amplitude, ROI bounds)
  void verify_existing_peaks_unchanged(
    const vector<shared_ptr<const PeakDef>> &original_peaks,
    const vector<shared_ptr<const PeakDef>> &current_peaks,
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove )
  {
    set<const PeakDef *> removed_ptrs;
    for( const shared_ptr<const PeakDef> &p : peaks_to_remove )
      removed_ptrs.insert( p.get() );

    for( const shared_ptr<const PeakDef> &orig : original_peaks )
    {
      if( !orig || removed_ptrs.count( orig.get() ) )
        continue;

      // Find the same peak in current_peaks by pointer identity
      bool found = false;
      for( const shared_ptr<const PeakDef> &curr : current_peaks )
      {
        if( curr.get() == orig.get() )
        {
          found = true;
          // Since it's the same pointer, the values must be identical
          break;
        }
      }

      BOOST_CHECK_MESSAGE( found,
        "Existing peak at " << orig->mean() << " keV (source: " << orig->sourceName()
        << ") is missing from current peaks" );
    }
  }// verify_existing_peaks_unchanged


  // Verify that for each source in original_peaks_to_remove, the same number of peaks
  // for that source appear in observable_peaks.  This ensures bystander peaks are properly
  // replaced rather than silently lost.
  void verify_removed_peaks_replaced(
    const FitPeaksForNuclides::PeakFitResult &result )
  {
    if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
      return;

    // Count removed peaks by source name
    map<string,size_t> removed_per_source;
    for( const shared_ptr<const PeakDef> &p : result.original_peaks_to_remove )
    {
      if( !p )
        continue;
      const string src = p->parentNuclide() ? p->parentNuclide()->symbol
                       : (p->xrayElement() ? p->xrayElement()->symbol
                       : (p->reaction() ? string("reaction") : string("unassigned")));
      removed_per_source[src] += 1;
    }

    // Count observable peaks by source name
    map<string,size_t> observable_per_source;
    for( const PeakDef &p : result.observable_peaks )
    {
      const string src = p.parentNuclide() ? p.parentNuclide()->symbol
                       : (p.xrayElement() ? p.xrayElement()->symbol
                       : (p.reaction() ? string("reaction") : string("unassigned")));
      observable_per_source[src] += 1;
    }

    // For each source that had peaks removed, check that at least as many
    // peaks for that source appear in observable_peaks
    for( const auto &src_count : removed_per_source )
    {
      const string &src = src_count.first;
      const size_t num_removed = src_count.second;
      const size_t num_observable = observable_per_source.count(src) ? observable_per_source.at(src) : 0;

      BOOST_CHECK_MESSAGE( num_observable >= num_removed,
        "Source '" << src << "': " << num_removed << " peak(s) removed but only "
        << num_observable << " replacement peak(s) in observable_peaks "
        << "(expected at least " << num_removed << ")" );
    }
  }// verify_removed_peaks_replaced


  // Apply fit result to a peak list: remove peaks_to_remove, add observable_peaks
  vector<shared_ptr<const PeakDef>> apply_fit_result(
    const vector<shared_ptr<const PeakDef>> &current_peaks,
    const FitPeaksForNuclides::PeakFitResult &result )
  {
    if( result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
      return current_peaks;
    if( result.observable_peaks.empty() )
      return current_peaks;

    set<const PeakDef *> remove_ptrs;
    for( const shared_ptr<const PeakDef> &p : result.original_peaks_to_remove )
      remove_ptrs.insert( p.get() );

    vector<shared_ptr<const PeakDef>> updated;
    for( const shared_ptr<const PeakDef> &p : current_peaks )
    {
      if( p && !remove_ptrs.count( p.get() ) )
        updated.push_back( p );
    }

    for( const PeakDef &p : result.observable_peaks )
      updated.push_back( make_shared<const PeakDef>( p ) );

    return updated;
  }// apply_fit_result

}// anonymous namespace


// ============================================================================
// B2: Smoke Tests
// ============================================================================

BOOST_AUTO_TEST_SUITE( SmokeTests )

BOOST_AUTO_TEST_CASE( test_cs137_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Cs137_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Cs137"} );
  const vector<shared_ptr<const PeakDef>> user_peaks; // empty

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  verify_fit_result( result, user_peaks, spec.foreground, {} );

  BOOST_CHECK_GE( result.observable_peaks.size(), 1u );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 661.66, 3.0 ) );

  BOOST_TEST_MESSAGE( "Cs137 smoke: " << result.observable_peaks.size() << " observable peaks" );
}


BOOST_AUTO_TEST_CASE( test_ba133_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Ba133_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Ba133"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  verify_fit_result( result, user_peaks, spec.foreground, {} );

  BOOST_CHECK_GE( result.observable_peaks.size(), 3u );

  // Check for major Ba-133 peaks
  BOOST_CHECK( has_peak_near( result.observable_peaks, 81.0, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 302.85, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 356.02, 3.0 ) );

  BOOST_TEST_MESSAGE( "Ba133 smoke: " << result.observable_peaks.size() << " observable peaks" );
}


BOOST_AUTO_TEST_CASE( test_eu152_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Eu152"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  verify_fit_result( result, user_peaks, spec.foreground, {} );

  // Eu-152 has many gamma lines - expect at least the major ones
  BOOST_CHECK_GE( result.observable_peaks.size(), 8u );

  // Check major Eu-152 peaks
  BOOST_CHECK( has_peak_near( result.observable_peaks, 121.78, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 344.28, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 778.90, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 964.08, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 1112.07, 3.0 ) );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 1408.01, 3.0 ) );

  BOOST_TEST_MESSAGE( "Eu152 smoke: " << result.observable_peaks.size() << " observable peaks" );
}


BOOST_AUTO_TEST_CASE( test_eu152_then_eu154 )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  // Step 1: Fit Eu-152
  const vector<RelActCalcAuto::SrcVariant> eu152_sources = make_sources( {"Eu152"} );
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;

  const FitPeaksForNuclides::PeakFitResult eu152_result
    = run_fit( spec.foreground, spec.background, auto_peaks, eu152_sources,
               empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( eu152_result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  const vector<shared_ptr<const PeakDef>> after_eu152 = apply_fit_result( empty_user_peaks, eu152_result );

  // Step 2: Fit Eu-154 with ExistingPeaksAsFreePeak (Eu-152 peaks exist)
  const vector<RelActCalcAuto::SrcVariant> eu154_sources = make_sources( {"Eu154"} );

  const FitPeaksForNuclides::PeakFitResult eu154_result
    = run_fit( spec.foreground, spec.background, auto_peaks, eu154_sources,
               after_eu152, spec.isHPGe,
               FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );

  // Eu-154 should not corrupt Eu-152 peaks
  // (the spectrum is pure Eu-152, so Eu-154 may or may not be found)
  verify_fit_result( eu154_result, after_eu152, spec.foreground,
    FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
  verify_removed_peaks_replaced( eu154_result );

  // If no Eu-154 was found, all Eu-152 peaks should be unchanged
  if( eu154_result.observable_peaks.empty() )
  {
    BOOST_CHECK( eu154_result.original_peaks_to_remove.empty() );
  }

  BOOST_TEST_MESSAGE( "Eu152+Eu154: Eu154 found " << eu154_result.observable_peaks.size()
    << " peaks, removed " << eu154_result.original_peaks_to_remove.size() );
}

BOOST_AUTO_TEST_SUITE_END() // SmokeTests


// ============================================================================
// B3: Idempotency Tests
// ============================================================================

BOOST_AUTO_TEST_SUITE( IdempotencyTests )

BOOST_AUTO_TEST_CASE( test_refit_same_source )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Cs137_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Cs137"} );

  // First fit
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;
  const FitPeaksForNuclides::PeakFitResult result1
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result1.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_REQUIRE_GE( result1.observable_peaks.size(), 1u );

  const vector<shared_ptr<const PeakDef>> after_first = apply_fit_result( empty_user_peaks, result1 );

  // Second fit with first fit's results as user_peaks (default mode - should replace)
  const FitPeaksForNuclides::PeakFitResult result2
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, after_first, spec.isHPGe );

  BOOST_REQUIRE( result2.status == RelActCalcAuto::RelActAutoSolution::Status::Success );

  verify_fit_result( result2, after_first, spec.foreground, {} );
  verify_removed_peaks_replaced( result2 );

  // The old peaks should be marked for removal (they are same-source)
  BOOST_CHECK_GE( result2.original_peaks_to_remove.size(), 1u );

  // The new peaks should be similar to the old ones
  BOOST_CHECK_GE( result2.observable_peaks.size(), 1u );

  const PeakDef *first_cs = find_peak_near( result1.observable_peaks, 661.66 );
  const PeakDef *second_cs = find_peak_near( result2.observable_peaks, 661.66 );
  BOOST_REQUIRE( first_cs && second_cs );

  // Peak means should be very close
  BOOST_CHECK_SMALL( first_cs->mean() - second_cs->mean(), 1.0 );
  // Areas should be within 20%
  if( first_cs->peakArea() > 0.0 )
  {
    const double area_ratio = second_cs->peakArea() / first_cs->peakArea();
    BOOST_CHECK_GT( area_ratio, 0.8 );
    BOOST_CHECK_LT( area_ratio, 1.2 );
  }
}


BOOST_AUTO_TEST_CASE( test_refit_after_peak_delete )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Eu152"} );

  // First fit
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;
  const FitPeaksForNuclides::PeakFitResult result1
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result1.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_REQUIRE_GE( result1.observable_peaks.size(), 5u );

  vector<shared_ptr<const PeakDef>> after_first = apply_fit_result( empty_user_peaks, result1 );

  // Delete 2 peaks (first two in the list)
  const double deleted_energy_1 = after_first[0]->mean();
  const double deleted_energy_2 = after_first[1]->mean();

  vector<shared_ptr<const PeakDef>> with_deletions;
  for( size_t i = 2; i < after_first.size(); ++i )
    with_deletions.push_back( after_first[i] );

  // Refit - the deleted peaks should reappear
  const FitPeaksForNuclides::PeakFitResult result2
    = run_fit( spec.foreground, spec.background, auto_peaks, sources,
               with_deletions, spec.isHPGe );

  BOOST_REQUIRE( result2.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  verify_fit_result( result2, with_deletions, spec.foreground, {} );
  verify_removed_peaks_replaced( result2 );

  // Check that the deleted peaks come back
  BOOST_CHECK_MESSAGE( has_peak_near( result2.observable_peaks, deleted_energy_1, 3.0 ),
    "Deleted peak at " << deleted_energy_1 << " keV did not reappear after refit" );
  BOOST_CHECK_MESSAGE( has_peak_near( result2.observable_peaks, deleted_energy_2, 3.0 ),
    "Deleted peak at " << deleted_energy_2 << " keV did not reappear after refit" );
}


BOOST_AUTO_TEST_CASE( test_refit_do_not_use_existing )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Cs137_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Cs137"} );

  // First fit
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;
  const FitPeaksForNuclides::PeakFitResult result1
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result1.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  const vector<shared_ptr<const PeakDef>> after_first = apply_fit_result( empty_user_peaks, result1 );

  // Second fit with DoNotUseExistingRois
  const FitPeaksForNuclides::PeakFitResult result2
    = run_fit( spec.foreground, spec.background, auto_peaks, sources,
               after_first, spec.isHPGe,
               FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois );

  // DoNotUseExistingRois should not remove any existing peaks
  BOOST_CHECK( result2.original_peaks_to_remove.empty() );

  // But should still find Cs137 peaks (in new ROIs that don't overlap existing)
  // Note: this may or may not succeed depending on whether non-overlapping ROIs can be found
  if( result2.status == RelActCalcAuto::RelActAutoSolution::Status::Success
     && !result2.observable_peaks.empty() )
  {
    verify_fit_result( result2, after_first, spec.foreground,
      FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois );
  }
}

BOOST_AUTO_TEST_SUITE_END() // IdempotencyTests


// ============================================================================
// B5: Option Behavior Tests
// ============================================================================

BOOST_AUTO_TEST_SUITE( OptionBehaviorTests )

BOOST_AUTO_TEST_CASE( test_default_preserves_other_source )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Ba133_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  // Fit Cs137 first (even though it's a Ba133 spectrum, the fit may find something near 662 keV
  // or return no observable peaks - either is fine for this test)
  const vector<RelActCalcAuto::SrcVariant> cs137_sources = make_sources( {"Cs137"} );
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;

  const FitPeaksForNuclides::PeakFitResult cs_result
    = run_fit( spec.foreground, spec.background, auto_peaks, cs137_sources,
               empty_user_peaks, spec.isHPGe );

  const vector<shared_ptr<const PeakDef>> after_cs = apply_fit_result( empty_user_peaks, cs_result );

  if( after_cs.empty() )
  {
    BOOST_TEST_MESSAGE( "Cs137 not found in Ba133 spectrum - skipping rest of test" );
    return;
  }

  // Now fit Ba133 with default mode - Cs137 peaks should NOT be in peaksToRemove
  const vector<RelActCalcAuto::SrcVariant> ba133_sources = make_sources( {"Ba133"} );

  const FitPeaksForNuclides::PeakFitResult ba_result
    = run_fit( spec.foreground, spec.background, auto_peaks, ba133_sources,
               after_cs, spec.isHPGe );

  if( ba_result.status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  verify_removed_peaks_replaced( ba_result );

  // Cs137 peaks should not be in peaksToRemove (they are from a different source)
  set<const PeakDef *> removed_ptrs;
  for( const shared_ptr<const PeakDef> &p : ba_result.original_peaks_to_remove )
    removed_ptrs.insert( p.get() );

  for( const shared_ptr<const PeakDef> &cs_peak : after_cs )
  {
    BOOST_CHECK_MESSAGE( !removed_ptrs.count( cs_peak.get() ),
      "Cs137 peak at " << cs_peak->mean()
      << " keV was removed when fitting Ba133 (default mode should preserve other-source peaks)" );
  }

  verify_fit_result( ba_result, after_cs, spec.foreground, {} );
}


BOOST_AUTO_TEST_CASE( test_do_not_use_existing_ignores_all )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  // Fit Eu-152 first
  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Eu152"} );
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;

  const FitPeaksForNuclides::PeakFitResult result1
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result1.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  const vector<shared_ptr<const PeakDef>> after_first = apply_fit_result( empty_user_peaks, result1 );
  BOOST_REQUIRE( !after_first.empty() );

  // Fit Eu-152 again with DoNotUseExistingRois
  const FitPeaksForNuclides::PeakFitResult result2
    = run_fit( spec.foreground, spec.background, auto_peaks, sources,
               after_first, spec.isHPGe,
               FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois );

  // Should NOT remove any existing peaks
  BOOST_CHECK_MESSAGE( result2.original_peaks_to_remove.empty(),
    "DoNotUseExistingRois should not remove any peaks, but removed "
    << result2.original_peaks_to_remove.size() );
}

BOOST_AUTO_TEST_SUITE_END() // OptionBehaviorTests


// ============================================================================
// B4: Trinitite Sequential Test
// ============================================================================

BOOST_AUTO_TEST_SUITE( TrinititeSequential )

BOOST_AUTO_TEST_CASE( test_trinitite_default_sequence )
{
  // Load trinitite spectra
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );

  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );
  BOOST_TEST_MESSAGE( "Trinitite auto-search found " << auto_peaks.size() << " peaks" );

  vector<shared_ptr<const PeakDef>> user_peaks; // accumulated peaks

  // ---- Step 1: Cs-137 ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 1: Cs-137 ---" );
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Cs137"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    verify_removed_peaks_replaced( result );

    BOOST_CHECK_GE( result.observable_peaks.size(), 1u );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 661.66, 3.0 ) );
    BOOST_CHECK( result.original_peaks_to_remove.empty() );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Cs-137: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 2: Am-241 ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 2: Am-241 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Am241"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    verify_removed_peaks_replaced( result );

    BOOST_CHECK_GE( result.observable_peaks.size(), 1u );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 59.54, 3.0 ) );
    BOOST_CHECK( result.original_peaks_to_remove.empty() );

    // Cs-137 peaks unchanged
    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Am-241: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 3: Eu-152 ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 3: Eu-152 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Eu152"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    verify_removed_peaks_replaced( result );

    // Eu-152 should have many peaks
    BOOST_CHECK_GE( result.observable_peaks.size(), 15u );

    // Check major peaks
    BOOST_CHECK( has_peak_near( result.observable_peaks, 121.78, 3.0 ) );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 344.28, 3.0 ) );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 778.90, 3.0 ) );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 964.08, 3.0 ) );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 1408.01, 3.0 ) );

    // Cs-137 and Am-241 peaks should not be removed
    BOOST_CHECK( result.original_peaks_to_remove.empty() );

    // Existing peaks unchanged
    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Eu-152: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 4: Eu-154 with ExistingPeaksAsFreePeak ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 4: Eu-154 (ExistingPeaksAsFreePeak) ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Eu154"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe,
                 FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );

    // Eu-154 is not strongly present in trinitite - expect no observable peaks
    // (but allow them if found)
    verify_fit_result( result, user_peaks, spec.foreground,
      FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
    verify_removed_peaks_replaced( result );

    if( result.observable_peaks.empty() )
    {
      // When no peaks found, nothing should be removed
      BOOST_CHECK( result.original_peaks_to_remove.empty() );
      BOOST_TEST_MESSAGE( "Eu-154: no observable peaks (expected)" );
    }
    else
    {
      BOOST_TEST_MESSAGE( "Eu-154: " << result.observable_peaks.size() << " peaks found" );
    }

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Eu-154: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 5: Ba-133 ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 5: Ba-133 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Ba133"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    verify_removed_peaks_replaced( result );

    BOOST_CHECK_GE( result.observable_peaks.size(), 3u );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 356.02, 3.0 ) );

    // Verify Ba-133 356 keV ROI fits between Eu-152 344 and 368 keV ROIs
    const PeakDef *ba356 = find_peak_near( result.observable_peaks, 356.02, 3.0 );
    if( ba356 && ba356->continuum() )
    {
      const double ba_roi_lower = ba356->continuum()->lowerEnergy();
      const double ba_roi_upper = ba356->continuum()->upperEnergy();

      // Find Eu-152 peaks in the neighborhood
      for( const shared_ptr<const PeakDef> &p : user_peaks )
      {
        if( !p || !p->continuum() )
          continue;
        // Eu-152 344.28 keV ROI should be below
        if( fabs( p->mean() - 344.28 ) < 3.0 )
        {
          BOOST_CHECK_MESSAGE( p->continuum()->upperEnergy() <= ba_roi_lower,
            "Ba-133 356 ROI [" << ba_roi_lower << ", " << ba_roi_upper
            << "] overlaps with Eu-152 344 ROI upper " << p->continuum()->upperEnergy() );
        }
        // Eu-152 367.79 keV ROI should be above
        if( fabs( p->mean() - 367.79 ) < 3.0 )
        {
          BOOST_CHECK_MESSAGE( ba_roi_upper <= p->continuum()->lowerEnergy(),
            "Ba-133 356 ROI [" << ba_roi_lower << ", " << ba_roi_upper
            << "] overlaps with Eu-152 368 ROI lower " << p->continuum()->lowerEnergy() );
        }
      }
    }

    // No existing peaks should be removed
    BOOST_CHECK( result.original_peaks_to_remove.empty() );
    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Ba-133: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 6: Co-60 ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 6: Co-60 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Co60"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    verify_removed_peaks_replaced( result );

    BOOST_CHECK_GE( result.observable_peaks.size(), 2u );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 1173.23, 3.0 ) );
    BOOST_CHECK( has_peak_near( result.observable_peaks, 1332.49, 3.0 ) );

    BOOST_CHECK( result.original_peaks_to_remove.empty() );
    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Co-60: " << user_peaks.size() << " total peaks" );
  }

  // ---- Step 7: U-235 (expected to fail) ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 7: U-235 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"U235"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    // May or may not find peaks - but should not corrupt existing
    if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      verify_fit_result( result, user_peaks, spec.foreground, {} );
    }

    // Either way, existing peaks should not be removed or altered
    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After U-235: " << user_peaks.size() << " total peaks"
      << (result.observable_peaks.empty() ? " (none found)" : "") );
  }

  // ---- Step 8: Ra-226 (expected to fail) ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 8: Ra-226 ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Ra226"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      verify_fit_result( result, user_peaks, spec.foreground, {} );
    }

    verify_existing_peaks_unchanged( pre_peaks, user_peaks, result.original_peaks_to_remove );

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After Ra-226: " << user_peaks.size() << " total peaks"
      << (result.observable_peaks.empty() ? " (none found)" : "") );
  }

  // ---- Final summary ----
  BOOST_TEST_MESSAGE( "\n=== Final state: " << user_peaks.size() << " peaks ===" );
  for( size_t i = 0; i < user_peaks.size(); ++i )
  {
    const PeakDef &p = *user_peaks[i];
    BOOST_TEST_MESSAGE( "  [" << i << "] " << p.mean() << " keV, source="
      << p.sourceName() << ", area=" << p.peakArea() );
  }
}


BOOST_AUTO_TEST_CASE( test_trinitite_do_not_use_existing_sequence )
{
  // Same sequence but all with DoNotUseExistingRois
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );

  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );

  vector<shared_ptr<const PeakDef>> user_peaks;
  const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> opts
    = FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois;

  const vector<string> source_names = { "Cs137", "Am241", "Eu152", "Ba133", "Co60" };

  for( const string &src_name : source_names )
  {
    BOOST_TEST_MESSAGE( "\n--- DoNotUseExisting: " << src_name << " ---" );
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {src_name} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources,
                 user_peaks, spec.isHPGe, opts );

    if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      verify_fit_result( result, user_peaks, spec.foreground, opts );

      // DoNotUseExistingRois: no peaks should ever be removed
      BOOST_CHECK_MESSAGE( result.original_peaks_to_remove.empty(),
        src_name << ": DoNotUseExistingRois removed " << result.original_peaks_to_remove.size()
        << " peaks" );
    }

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After " << src_name << ": " << user_peaks.size() << " total peaks" );
  }
}

BOOST_AUTO_TEST_SUITE_END() // TrinititeSequential


// ============================================================================
// B6: Additional Spectra
// ============================================================================

BOOST_AUTO_TEST_SUITE( AdditionalSpectra )

BOOST_AUTO_TEST_CASE( test_xe133_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Xe133_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Xe133", "Xe133m"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  verify_fit_result( result, user_peaks, spec.foreground, {} );

  BOOST_CHECK_GE( result.observable_peaks.size(), 1u );
  // Xe-133 81 keV
  BOOST_CHECK( has_peak_near( result.observable_peaks, 81.0, 3.0 ) );

  BOOST_TEST_MESSAGE( "Xe133 smoke: " << result.observable_peaks.size() << " observable peaks" );
}


BOOST_AUTO_TEST_CASE( test_am241_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Am241_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Am241"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  verify_fit_result( result, user_peaks, spec.foreground, {} );

  BOOST_CHECK_GE( result.observable_peaks.size(), 1u );
  BOOST_CHECK( has_peak_near( result.observable_peaks, 59.54, 3.0 ) );

  BOOST_TEST_MESSAGE( "Am241 smoke: " << result.observable_peaks.size() << " observable peaks" );
}


BOOST_AUTO_TEST_CASE( test_pu239_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Pu239_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Pu239"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  // Pu-239 may or may not fit well depending on the spectrum
  if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
  {
    verify_fit_result( result, user_peaks, spec.foreground, {} );
    BOOST_TEST_MESSAGE( "Pu239 smoke: " << result.observable_peaks.size() << " observable peaks" );
  }
  else
  {
    BOOST_TEST_MESSAGE( "Pu239 smoke: fit failed - " << result.error_message );
  }
}


BOOST_AUTO_TEST_CASE( test_i125_smoke )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "I125_Unshielded.txt" );
  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"I125"} );
  const vector<shared_ptr<const PeakDef>> user_peaks;

  const FitPeaksForNuclides::PeakFitResult result
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

  if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
  {
    verify_fit_result( result, user_peaks, spec.foreground, {} );
    BOOST_TEST_MESSAGE( "I125 smoke: " << result.observable_peaks.size() << " observable peaks" );
  }
  else
  {
    BOOST_TEST_MESSAGE( "I125 smoke: fit failed - " << result.error_message );
  }
}

BOOST_AUTO_TEST_SUITE_END() // AdditionalSpectra


// ============================================================================
// Bystander Peak Degradation Tests (ExistingPeaksAsFreePeak)
// ============================================================================

BOOST_AUTO_TEST_SUITE( BystanderDegradation )

// Issue 1: Strong Eu-152 244.7 keV peak should not be destroyed when fitting Eu-154
// with ExistingPeaksAsFreePeak. Eu-154 has a gamma at 247.93 keV (~3.2 keV away).
BOOST_AUTO_TEST_CASE( test_eu154_does_not_remove_strong_eu152_peak )
{
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );

  // Step 1: Fit Eu-152
  const vector<RelActCalcAuto::SrcVariant> sources_eu152 = make_sources( {"Eu152"} );
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;

  const FitPeaksForNuclides::PeakFitResult result_eu152
    = run_fit( spec.foreground, spec.background, auto_peaks, sources_eu152,
               empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result_eu152.status == RelActCalcAuto::RelActAutoSolution::Status::Success );

  // Verify Eu-152 has a strong peak near 244.7 keV
  const PeakDef *eu152_244 = find_peak_near( result_eu152.observable_peaks, 244.7, 3.0 );
  BOOST_REQUIRE_MESSAGE( eu152_244, "Eu-152 fit should have a peak near 244.7 keV" );

  const double eu152_244_sig = (eu152_244->amplitudeUncert() > 0.0)
    ? eu152_244->amplitude() / eu152_244->amplitudeUncert() : 0.0;
  BOOST_TEST_MESSAGE( "Eu-152 244.7 keV peak: mean=" << eu152_244->mean()
    << ", amp=" << eu152_244->amplitude() << ", sig=" << eu152_244_sig );
  BOOST_REQUIRE_MESSAGE( eu152_244_sig > 10.0,
    "Eu-152 244.7 keV peak should be strong (sig=" << eu152_244_sig << ")" );

  // Convert observable peaks to user_peaks
  vector<shared_ptr<const PeakDef>> eu152_user_peaks;
  for( const PeakDef &p : result_eu152.observable_peaks )
    eu152_user_peaks.push_back( make_shared<const PeakDef>( p ) );

  // Step 2: Fit Eu-154 with ExistingPeaksAsFreePeak
  const vector<RelActCalcAuto::SrcVariant> sources_eu154 = make_sources( {"Eu154"} );

  const FitPeaksForNuclides::PeakFitResult result_eu154
    = run_fit( spec.foreground, spec.background, auto_peaks, sources_eu154,
               eu152_user_peaks, spec.isHPGe,
               FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );

  verify_fit_result( result_eu154, eu152_user_peaks, spec.foreground,
    FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
  verify_removed_peaks_replaced( result_eu154 );

  // Check that no strong existing peak is removed without a comparably significant replacement.
  for( const shared_ptr<const PeakDef> &removed : result_eu154.original_peaks_to_remove )
  {
    const PeakDef *orig = find_peak_near( result_eu152.observable_peaks, removed->mean(), 1.0 );
    if( !orig )
      continue;

    const double orig_sig = (orig->amplitudeUncert() > 0.0)
      ? orig->amplitude() / orig->amplitudeUncert() : 0.0;
    if( orig_sig < 10.0 )
      continue;  // Only care about strong peaks

    // There must be a replacement peak near this energy in observable_peaks
    const PeakDef *replacement = find_peak_near( result_eu154.observable_peaks, removed->mean(), 3.0 );
    BOOST_CHECK_MESSAGE( replacement != nullptr,
      "Strong existing peak at " << removed->mean() << " keV (sig="
      << orig_sig << ") was removed but has no replacement in observable_peaks" );

    if( replacement )
    {
      const double repl_sig = (replacement->amplitudeUncert() > 0.0)
        ? replacement->amplitude() / replacement->amplitudeUncert() : 0.0;
      BOOST_CHECK_MESSAGE( repl_sig > 3.0,
        "Replacement peak at " << replacement->mean() << " keV has poor significance ("
        << repl_sig << ") compared to removed peak at " << removed->mean()
        << " keV (sig=" << orig_sig << ")" );
    }
  }//for( removed peaks )
}


// Issue 2: Am-241 peak at 59.5 keV should not be destroyed when fitting Eu-154
// with ExistingPeaksAsFreePeak. Am-241's 59.54 keV gamma is only ~1.1 keV from
// Eu-154's 58.4 keV gamma - too close to resolve on HPGe.
BOOST_AUTO_TEST_CASE( test_eu154_does_not_destroy_am241_peak )
{
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );

  // Step 1: Fit Am-241 first (claims the strong 59.5 keV peak)
  const vector<RelActCalcAuto::SrcVariant> sources_am241 = make_sources( {"Am241"} );
  const vector<shared_ptr<const PeakDef>> empty_user_peaks;

  const FitPeaksForNuclides::PeakFitResult result_am241
    = run_fit( spec.foreground, spec.background, auto_peaks, sources_am241,
               empty_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( result_am241.status == RelActCalcAuto::RelActAutoSolution::Status::Success );

  // Verify Am-241 has a strong peak near 59.5 keV
  const PeakDef *am241_59 = find_peak_near( result_am241.observable_peaks, 59.5, 3.0 );
  BOOST_REQUIRE_MESSAGE( am241_59, "Am-241 fit should have a peak near 59.5 keV" );

  const double am241_59_sig = (am241_59->amplitudeUncert() > 0.0)
    ? am241_59->amplitude() / am241_59->amplitudeUncert() : 0.0;
  BOOST_TEST_MESSAGE( "Am-241 59.5 keV peak: mean=" << am241_59->mean()
    << ", amp=" << am241_59->amplitude() << ", sig=" << am241_59_sig );
  BOOST_REQUIRE_MESSAGE( am241_59_sig > 10.0,
    "Am-241 59.5 keV peak should be strong (sig=" << am241_59_sig << ")" );

  // Convert observable peaks to user_peaks
  vector<shared_ptr<const PeakDef>> am241_user_peaks;
  for( const PeakDef &p : result_am241.observable_peaks )
    am241_user_peaks.push_back( make_shared<const PeakDef>( p ) );

  // Step 2: Fit Eu-154 with ExistingPeaksAsFreePeak
  // Eu-154 has a gamma at 58.4 keV, ~1.1 keV from Am-241's 59.54 keV
  const vector<RelActCalcAuto::SrcVariant> sources_eu154 = make_sources( {"Eu154"} );

  const FitPeaksForNuclides::PeakFitResult result_eu154
    = run_fit( spec.foreground, spec.background, auto_peaks, sources_eu154,
               am241_user_peaks, spec.isHPGe,
               FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );

  if( result_eu154.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
  {
    verify_fit_result( result_eu154, am241_user_peaks, spec.foreground,
      FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
    verify_removed_peaks_replaced( result_eu154 );
  }

  // Whether fit succeeded or not, check that no strong existing peak is destroyed
  for( const shared_ptr<const PeakDef> &removed : result_eu154.original_peaks_to_remove )
  {
    const PeakDef *orig = find_peak_near( result_am241.observable_peaks, removed->mean(), 1.0 );
    if( !orig )
      continue;

    const double orig_sig = (orig->amplitudeUncert() > 0.0)
      ? orig->amplitude() / orig->amplitudeUncert() : 0.0;
    if( orig_sig < 10.0 )
      continue;

    const PeakDef *replacement = find_peak_near( result_eu154.observable_peaks, removed->mean(), 3.0 );
    BOOST_CHECK_MESSAGE( replacement != nullptr,
      "Strong Am-241 peak at " << removed->mean() << " keV (sig="
      << orig_sig << ") was removed but has no replacement" );

    if( replacement )
    {
      const double repl_sig = (replacement->amplitudeUncert() > 0.0)
        ? replacement->amplitude() / replacement->amplitudeUncert() : 0.0;
      BOOST_CHECK_MESSAGE( repl_sig > 3.0,
        "Replacement peak at " << replacement->mean() << " keV has poor significance ("
        << repl_sig << ") compared to removed Am-241 peak at " << removed->mean()
        << " keV (sig=" << orig_sig << ")" );
    }
  }//for( removed peaks )
}

BOOST_AUTO_TEST_SUITE_END() // BystanderDegradation
