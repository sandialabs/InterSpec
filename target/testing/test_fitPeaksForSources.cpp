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
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <numeric>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <functional>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
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
    auto prefs = make_shared<PeakFitDetPrefs>();
    prefs->m_det_type = isHPGe ? PeakFitUtils::CoarseResolutionType::High
                               : PeakFitUtils::CoarseResolutionType::Low;
    return ExperimentalAutomatedPeakSearch::search_for_peaks(
      foreground, nullptr, nullptr, true, prefs );
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
    const PeakFitUtils::CoarseResolutionType det_type = isHPGe
        ? PeakFitUtils::CoarseResolutionType::High
        : PeakFitUtils::CoarseResolutionType::Low;
    auto peak_fit_prefs = make_shared<PeakFitDetPrefs>();
    peak_fit_prefs->m_det_type = det_type;

    const FitPeaksForNuclides::PeakFitForNuclideConfig &config
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( det_type );

    return FitPeaksForNuclides::fit_peaks_for_nuclides(
      auto_search_peaks, foreground, sources, user_peaks,
      background, nullptr, options, config, peak_fit_prefs );
  }// run_fit


  FitPeaksForNuclides::PeakFitResult run_fit_with_config(
    const shared_ptr<const SpecUtils::Measurement> &foreground,
    const shared_ptr<const SpecUtils::Measurement> &background,
    const vector<shared_ptr<const PeakDef>> &auto_search_peaks,
    const vector<RelActCalcAuto::SrcVariant> &sources,
    const vector<shared_ptr<const PeakDef>> &user_peaks,
    const bool isHPGe,
    const FitPeaksForNuclides::PeakFitForNuclideConfig &config,
    const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> options
      = Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions>() )
  {
    const PeakFitUtils::CoarseResolutionType det_type = isHPGe
        ? PeakFitUtils::CoarseResolutionType::High
        : PeakFitUtils::CoarseResolutionType::Low;
    auto peak_fit_prefs = make_shared<PeakFitDetPrefs>();
    peak_fit_prefs->m_det_type = det_type;
    return FitPeaksForNuclides::fit_peaks_for_nuclides(
      auto_search_peaks, foreground, sources, user_peaks,
      background, nullptr, options, config, peak_fit_prefs );
  }//run_fit_with_config(...)


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


  const PeakDef *find_source_gamma( const vector<PeakDef> &peaks, const string &symbol,
                                    const double gamma_energy,
                                    const double tolerance_keV = 3.0 )
  {
    const PeakDef *best = nullptr;
    double best_distance = tolerance_keV;
    for( const PeakDef &peak : peaks )
    {
      if( !peak.parentNuclide() || (peak.parentNuclide()->symbol != symbol)
          || !peak.hasSourceGammaAssigned() )
        continue;
      const double distance = std::fabs( peak.gammaParticleEnergy() - gamma_energy );
      if( distance < best_distance )
      {
        best = &peak;
        best_distance = distance;
      }
    }
    return best;
  }//find_source_gamma


  bool peak_areas_agree( const PeakDef &lhs, const PeakDef &rhs )
  {
    const double combined_uncert = std::hypot(
        std::max( 0.0, lhs.peakAreaUncert() ), std::max( 0.0, rhs.peakAreaUncert() ) );
    const double tolerance = std::max( 3.0*combined_uncert,
                                      0.20*std::fabs(rhs.peakArea()) );
    return std::fabs( lhs.peakArea() - rhs.peakArea() ) <= tolerance;
  }//peak_areas_agree


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
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( PeakFitUtils::CoarseResolutionType::High );

    // 1. No overlapping observable ROIs
    verify_no_roi_overlaps( result.observable_peaks, foreground );

    // 2. Peak means within ROI bounds
    verify_peaks_within_roi( result.observable_peaks );

    // 3. Continuum consistency
    verify_continuum_consistency( user_peaks, result.original_peaks_to_remove );

    // 4. FWHM margin from existing peaks
    verify_fwhm_margin_from_existing( result.observable_peaks, user_peaks,
      result.original_peaks_to_remove, config, foreground );

    // FitPeaksForNuclides finalizes its own automatic geometry.  Every range presented to
    // RelActAuto must therefore be Fixed and channel-disjoint, closing the solver's generic
    // CanBeBrokenUp significant-range recombination path for late found/edge/escape additions.
    const bool policy_enabled
      = (options & FitPeaksForNuclides::FitSrcPeaksOptions::DisableAutoInterfererFit)
        || (options & FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks)
        || (options & FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse);
    const shared_ptr<const SpecUtils::Measurement> fitted_foreground
      = result.solution.m_foreground ? result.solution.m_foreground : foreground;
    for( size_t roi_index = 0;
         policy_enabled && (roi_index < result.solution.m_options.rois.size()); ++roi_index )
    {
      const RelActCalcAuto::RoiRange &roi = result.solution.m_options.rois[roi_index];
      BOOST_CHECK( roi.range_limits_type
          == RelActCalcAuto::RoiRange::RangeLimitsType::Fixed );
      if( roi_index && fitted_foreground )
      {
        const RelActCalcAuto::RoiRange &previous
          = result.solution.m_options.rois[roi_index - 1];
        BOOST_CHECK_GT( fitted_foreground->find_gamma_channel(roi.lower_energy),
                        fitted_foreground->find_gamma_channel(previous.upper_energy) );
      }
    }

    // Every atom-safe partition that fell back to retaining incumbent geometry must document why.
    if( policy_enabled )
    {
      for( const FitPeaksForNuclides::AutomaticRoiDecisionDiagnostic &diag
             : result.automatic_roi_diagnostics )
      {
        if( diag.partition_infeasible )
          BOOST_CHECK_MESSAGE( !diag.reason.empty(),
            "partition_infeasible diagnostic at stage '" << diag.stage
            << "' has no reason recorded" );
      }
    }

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


BOOST_AUTO_TEST_CASE( test_rescue_recovers_marginal_line )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );

  // Remove one moderate Eu line from the data-confirmed seeding path.  The initial manual solve
  // still has the source's predicted ROIs, while the deliberately high fitted keep threshold
  // creates a genuine marginal-reject band for the one bounded R2 pass.
  const double withheld_energy = 1249.93;
  auto_peaks.erase( std::remove_if( std::begin(auto_peaks), std::end(auto_peaks),
    [withheld_energy]( const shared_ptr<const PeakDef> &peak ) {
      return peak && (std::fabs(peak->mean() - withheld_energy) < 2.0);
    } ), std::end(auto_peaks) );

  FitPeaksForNuclides::PeakFitForNuclideConfig config
    = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config(
        PeakFitUtils::CoarseResolutionType::High );
  config.auto_keep_significance_z = 5.5;
  config.manual_keep_significance_z = 5.5;

  const vector<shared_ptr<const PeakDef>> no_user_peaks;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_ecal_options;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;
  FitPeaksForNuclides::PeakFitResult result = run_fit_with_config(
      spec.foreground, spec.background, auto_peaks, make_sources({"Eu152"}),
      no_user_peaks, spec.isHPGe, config, no_ecal_options );
  BOOST_REQUIRE( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  for( const string &warning : result.warnings )
    BOOST_TEST_MESSAGE( "R2 recovery: " << warning );

  const bool admitted_rescue = std::any_of( std::begin(result.warnings),
    std::end(result.warnings), []( const string &warning ) {
      return warning.find("bounded fit-then-prune rescue") != string::npos;
    } );
  BOOST_CHECK_MESSAGE( admitted_rescue,
                       "High keep-z fit did not exercise the bounded rescue path" );
  BOOST_CHECK_MESSAGE( find_source_gamma(result.observable_peaks, "Eu152", withheld_energy, 0.5),
                       "The real withheld Eu152 marginal line was not restored" );

#if( PERFORM_DEVELOPER_CHECKS )
  struct RestoreRescue
  {
    ~RestoreRescue()
    {
      FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( true );
    }
  } restore_rescue;
  FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( false );
  FitPeaksForNuclides::PeakFitResult control = run_fit_with_config(
      spec.foreground, spec.background, auto_peaks, make_sources({"Eu152"}),
      no_user_peaks, spec.isHPGe, config, no_ecal_options );
  BOOST_REQUIRE( control.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  vector<double> rescued_only_energies;
  const auto find_rescued_only = [&]() {
    rescued_only_energies.clear();
    for( const PeakDef &peak : result.observable_peaks )
    {
      if( !peak.parentNuclide() || (peak.parentNuclide()->symbol != "Eu152")
          || !peak.hasSourceGammaAssigned() )
        continue;
      if( !find_source_gamma(control.observable_peaks, "Eu152",
                             peak.gammaParticleEnergy(), 0.2) )
        rescued_only_energies.push_back( peak.gammaParticleEnergy() );
    }
  };
  find_rescued_only();
  for( const double keep_z : { 6.5, 8.0, 10.0, 12.0, 15.0 } )
  {
    if( !rescued_only_energies.empty() )
      break;
    config.auto_keep_significance_z = keep_z;
    FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( true );
    result = run_fit_with_config( spec.foreground, spec.background, auto_peaks,
        make_sources({"Eu152"}), no_user_peaks, spec.isHPGe, config, no_ecal_options );
    FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( false );
    control = run_fit_with_config( spec.foreground, spec.background, auto_peaks,
        make_sources({"Eu152"}), no_user_peaks, spec.isHPGe, config, no_ecal_options );
    BOOST_REQUIRE( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
    BOOST_REQUIRE( control.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
    find_rescued_only();
  }
  for( const double energy : rescued_only_energies )
    BOOST_TEST_MESSAGE( "R2-only recovered Eu152 line at " << energy << " keV" );
  BOOST_CHECK_MESSAGE( !rescued_only_energies.empty(),
                       "R2-enabled and R2-disabled controls returned the same source lines" );
  const PeakDef * const rescued_peak = rescued_only_energies.empty() ? nullptr
    : find_source_gamma( result.observable_peaks, "Eu152", rescued_only_energies.front(), 0.2 );
  BOOST_REQUIRE( rescued_peak );
  BOOST_CHECK_GT( rescued_peak->peakArea(), 0.0 );
#endif
}


BOOST_AUTO_TEST_CASE( test_rescue_exception_retains_successful_incumbent )
{
#if( PERFORM_DEVELOPER_CHECKS )
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const double rescued_energy = 1249.93;
  auto_peaks.erase( std::remove_if( std::begin(auto_peaks), std::end(auto_peaks),
    [rescued_energy]( const shared_ptr<const PeakDef> &peak ) {
      return peak && (std::fabs(peak->mean() - rescued_energy) < 2.0);
    } ), std::end(auto_peaks) );
  FitPeaksForNuclides::PeakFitForNuclideConfig config
    = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config(
        PeakFitUtils::CoarseResolutionType::High );
  config.manual_keep_significance_z = 5.5;
  config.auto_keep_significance_z = 5.5;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_ecal_options;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;

  FitPeaksForNuclides::detail::force_next_bounded_rescue_evaluation_failure_for_test();
  const FitPeaksForNuclides::PeakFitResult result = run_fit_with_config(
      spec.foreground, spec.background, auto_peaks, make_sources({"Eu152"}), {},
      spec.isHPGe, config, no_ecal_options );
  BOOST_REQUIRE( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_CHECK( find_source_gamma(result.observable_peaks, "Eu152", 344.28, 0.5) );
  BOOST_CHECK( !result.observable_peaks.empty() );
  BOOST_CHECK( std::any_of(std::begin(result.warnings), std::end(result.warnings),
    []( const string &warning ) {
      return warning.find("rescue challenger threw") != string::npos;
    }) );
#endif
}


BOOST_AUTO_TEST_CASE( test_rescue_is_source_order_and_background_invariant )
{
  const LoadedSpectrum spec = load_detective_x_spectrum( "Eu152_Unshielded.txt" );
  vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const double rescued_energy = 1249.93;
  auto_peaks.erase( std::remove_if( std::begin(auto_peaks), std::end(auto_peaks),
    [rescued_energy]( const shared_ptr<const PeakDef> &peak ) {
      return peak && (std::fabs(peak->mean() - rescued_energy) < 2.0);
    } ), std::end(auto_peaks) );

  FitPeaksForNuclides::PeakFitForNuclideConfig config
    = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config(
        PeakFitUtils::CoarseResolutionType::High );
  config.manual_keep_significance_z = 5.5;
  config.auto_keep_significance_z = 5.5;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_ecal_options;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;

  const auto fit = [&]( const shared_ptr<const SpecUtils::Measurement> &background,
                        const vector<string> &source_names ) {
    return run_fit_with_config( spec.foreground, background, auto_peaks,
        make_sources(source_names), {}, spec.isHPGe, config, no_ecal_options );
  };
  const FitPeaksForNuclides::PeakFitResult results[] = {
    fit( nullptr, {"Eu152", "Cs137"} ),
    fit( nullptr, {"Cs137", "Eu152"} ),
    fit( spec.background, {"Eu152", "Cs137"} ),
    fit( spec.background, {"Cs137", "Eu152"} )
  };
  const char * const labels[] = {
    "raw Eu152,Cs137", "raw Cs137,Eu152",
    "supplied-background Eu152,Cs137", "supplied-background Cs137,Eu152"
  };

  const PeakDef *reference = nullptr;
  bool observable_presence[4] = { false, false, false, false };
  for( size_t result_index = 0; result_index < 4; ++result_index )
  {
    const FitPeaksForNuclides::PeakFitResult &result = results[result_index];
    BOOST_TEST_CONTEXT( labels[result_index] )
    {
    BOOST_REQUIRE( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
    const PeakDef * const rescued
      = find_source_gamma( result.uncombined_fit_peaks, "Eu152", rescued_energy, 0.5 );
    observable_presence[result_index] = static_cast<bool>(
      find_source_gamma(result.observable_peaks, "Eu152", rescued_energy, 0.5) );
    if( !rescued )
    {
      for( const string &warning : result.warnings )
        BOOST_TEST_MESSAGE( labels[result_index] << " warning: " << warning );
      const auto report_candidate = [&]( const char *name, const vector<PeakDef> &peaks ) {
        const PeakDef * const candidate
          = find_source_gamma( peaks, "Eu152", rescued_energy, 0.5 );
        BOOST_TEST_MESSAGE( labels[result_index] << ' ' << name << " 1249.93-keV peak: "
          << (candidate ? (std::to_string(candidate->peakArea()) + " +/- "
              + std::to_string(candidate->peakAreaUncert())) : string("absent")) );
      };
      report_candidate( "solution", result.solution.m_peaks_without_back_sub );
      report_candidate( "uncombined", result.uncombined_fit_peaks );
      report_candidate( "combined", result.fit_peaks );
      for( const PeakDef &peak : result.observable_peaks )
      {
        if( (peak.mean() < 1240.0) || (peak.mean() > 1260.0) )
          continue;
        BOOST_TEST_MESSAGE( labels[result_index] << " observable peak near rescue: mean="
          << peak.mean() << ", area=" << peak.peakArea() << " +/- " << peak.peakAreaUncert()
          << ", source=" << peak.sourceName() );
      }
    }
    BOOST_REQUIRE( rescued );
    BOOST_CHECK( std::any_of(std::begin(result.warnings), std::end(result.warnings),
      []( const string &warning ) {
        return warning.find("bounded fit-then-prune rescue") != string::npos;
      }) );
    if( reference )
      BOOST_CHECK_MESSAGE( peak_areas_agree(*reference, *rescued),
                           "R2 rescued area changed with source order/background mode" );
    else
      reference = rescued;
    }
  }
  BOOST_CHECK_EQUAL( observable_presence[0], observable_presence[1] );
  BOOST_CHECK_EQUAL( observable_presence[2], observable_presence[3] );
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

  // ---- Step 3: Eu-152 with a supplied background (the foreground-only R6 path is disabled) ----
  {
    BOOST_TEST_MESSAGE( "\n--- Step 3: Eu-152 (supplied-background path) ---" );
    const vector<shared_ptr<const PeakDef>> pre_peaks = user_peaks;
    // Eu-152's weak 1457 keV line sits on K40's strong 1460 keV line.  This sequence deliberately
    // supplies the measured background, so it verifies the legacy background-aware path rather than
    // the foreground-only R6 nuisance search exercised by test_r6_raw_interferer_transaction.
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

    // R6: the co-fit K40 interferer must NOT appear as a returned peak (dropped after co-fit); and
    // no spurious peak should be attributed to Eu-152 right on K40's 1460.8 keV line.
    for( const PeakDef &p : result.observable_peaks )
    {
      const bool is_k40 = p.parentNuclide() && (p.parentNuclide()->symbol == "K40");
      BOOST_CHECK_MESSAGE( !is_k40, "Eu-152-alone returned a K40-attributed peak at " << p.mean() << " keV" );
    }

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

    // Eu-154 is absent at useful strength in this spectrum.  This is the established bystander
    // regression: the bounded rescue pass must not resurrect Eu-154 from already-modeled peaks.
    verify_fit_result( result, user_peaks, spec.foreground,
      FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
    verify_removed_peaks_replaced( result );
    BOOST_CHECK_MESSAGE( result.observable_peaks.empty(),
                         "Bounded rescue resurrected absent Eu154 in trinitite" );
    BOOST_CHECK( result.original_peaks_to_remove.empty() );
    BOOST_TEST_MESSAGE( "Eu-154: no observable peaks (expected)" );

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

  // Fit each source independently (DoNotUseExistingRois) with a supplied background.  Automatic R6
  // discovery is intentionally inactive here; the raw-spectrum transaction tests cover that path.
  const vector<vector<string>> source_groups = {
    {"Cs137"}, {"Am241"}, {"Eu152"}, {"Ba133"}, {"Co60"}
  };

  for( const vector<string> &src_group : source_groups )
  {
    string group_label;
    for( const string &n : src_group )
      group_label += (group_label.empty() ? "" : "+") + n;

    BOOST_TEST_MESSAGE( "\n--- DoNotUseExisting: " << group_label << " ---" );
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( src_group );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources,
                 user_peaks, spec.isHPGe, opts );

    if( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      verify_fit_result( result, user_peaks, spec.foreground, opts );

      // DoNotUseExistingRois: no peaks should ever be removed
      BOOST_CHECK_MESSAGE( result.original_peaks_to_remove.empty(),
        group_label << ": DoNotUseExistingRois removed " << result.original_peaks_to_remove.size()
        << " peaks" );
    }

    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "After " << group_label << ": " << user_peaks.size() << " total peaks" );
  }
}

// Regression (2026-07 review, P1): the FitNormBkgrndPeaks path used to rebuild options.rois from
// the raw input ROIs, discarding the existing-ROI trimming and mixed-ROI setup - so NORM/source
// ROIs could cover existing other-source user peaks.  Fit Ba-133 first, then Cs-137 with NORM
// background peaks, and require that no new observable ROI covers the mean of a retained
// (not-removed) existing user peak.  Note: NORM fits disable background subtraction, so this case
// is on the slow side.
BOOST_AUTO_TEST_CASE( test_norm_fit_preserves_existing_rois )
{
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );

  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );

  // ---- Step 1: Ba-133, default mode, to establish existing user peaks/ROIs ----
  vector<shared_ptr<const PeakDef>> user_peaks;
  {
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Ba133"} );
    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe );

    verify_fit_result( result, user_peaks, spec.foreground, {} );
    BOOST_REQUIRE_GE( result.observable_peaks.size(), 1u );
    user_peaks = apply_fit_result( user_peaks, result );
    BOOST_TEST_MESSAGE( "Ba-133 established " << user_peaks.size() << " existing peaks" );
  }

  // ---- Step 2: Cs-137 with NORM background peaks, against the existing Ba-133 peaks ----
  {
    const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> opts
      = FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks;
    const vector<RelActCalcAuto::SrcVariant> sources = make_sources( {"Cs137"} );

    const FitPeaksForNuclides::PeakFitResult result
      = run_fit( spec.foreground, spec.background, auto_peaks, sources, user_peaks, spec.isHPGe, opts );

    verify_fit_result( result, user_peaks, spec.foreground, opts );
    verify_removed_peaks_replaced( result );

    BOOST_CHECK( has_peak_near( result.observable_peaks, 661.66, 3.0 ) );

    // The P1 regression check: a retained (not-removed) existing user peak's mean must not be
    // covered by any new observable ROI (per the FitSrcPeaksOptions default-mode contract).
    // A Ba-133 peak MAY legitimately be removed+replaced (mixed-ROI bystander flow, e.g. the
    // 356 keV ROI vs Pb-214 352 keV from the Ra-226 NORM chain) - those are excluded here.
    set<const PeakDef *> removed_ptrs;
    for( const shared_ptr<const PeakDef> &p : result.original_peaks_to_remove )
      removed_ptrs.insert( p.get() );

    for( const shared_ptr<const PeakDef> &user_peak : user_peaks )
    {
      if( !user_peak || removed_ptrs.count( user_peak.get() ) )
        continue;

      const double mean = user_peak->mean();
      set<const PeakContinuum *> seen_conts;
      for( const PeakDef &obs : result.observable_peaks )
      {
        if( !obs.continuum() || !seen_conts.insert( obs.continuum().get() ).second )
          continue;

        const bool covers = (mean >= obs.continuum()->lowerEnergy())
                            && (mean <= obs.continuum()->upperEnergy());
        BOOST_CHECK_MESSAGE( !covers,
          "NORM-fit observable ROI [" << obs.continuum()->lowerEnergy() << ", "
          << obs.continuum()->upperEnergy() << "] keV covers retained existing "
          << user_peak->sourceName() << " peak mean at " << mean << " keV" );
      }
    }
  }
}


BOOST_AUTO_TEST_CASE( test_eu152_interferer_matches_joint )
{
  // With a supplied background, the Eu-152-only fit should agree with an explicit {K40, Eu152}
  // reference and must not return a K40-attributed public peak.  The active foreground-only R6
  // comparison is covered more strictly by test_multisource_strong_norm_interferer_is_stable.
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<shared_ptr<const PeakDef>> no_user_peaks;

  const FitPeaksForNuclides::PeakFitResult r_auto
    = run_fit( spec.foreground, spec.background, auto_peaks,
               make_sources( {"Eu152"} ), no_user_peaks, spec.isHPGe );
  const FitPeaksForNuclides::PeakFitResult r_joint
    = run_fit( spec.foreground, spec.background, auto_peaks,
               make_sources( {"K40", "Eu152"} ), no_user_peaks, spec.isHPGe );

  BOOST_REQUIRE( r_auto.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_REQUIRE( r_joint.status == RelActCalcAuto::RelActAutoSolution::Status::Success );

  auto count_src = []( const vector<PeakDef> &peaks, const char *sym ) -> size_t {
    size_t n = 0;
    for( const PeakDef &p : peaks )
      if( p.parentNuclide() && (p.parentNuclide()->symbol == sym) )
        ++n;
    return n;
  };

  BOOST_TEST_MESSAGE( "auto: " << r_auto.observable_peaks.size() << " peaks ("
    << count_src(r_auto.observable_peaks,"Eu152") << " Eu152, "
    << count_src(r_auto.observable_peaks,"K40") << " K40); joint: "
    << r_joint.observable_peaks.size() << " peaks ("
    << count_src(r_joint.observable_peaks,"Eu152") << " Eu152, "
    << count_src(r_joint.observable_peaks,"K40") << " K40)" );

  // Diagnostic: dump the 1453-1465 keV region for both fits.
  for( const PeakDef &p : r_auto.observable_peaks )
    if( (p.mean() > 1453.0) && (p.mean() < 1465.0) )
      BOOST_TEST_MESSAGE( "  auto  1460-region peak@" << p.mean() << " area=" << p.amplitude()
        << " src=" << (p.parentNuclide()?p.parentNuclide()->symbol:string("none")) );
  for( const PeakDef &p : r_joint.observable_peaks )
    if( (p.mean() > 1453.0) && (p.mean() < 1465.0) )
      BOOST_TEST_MESSAGE( "  joint 1460-region peak@" << p.mean() << " area=" << p.amplitude()
        << " src=" << (p.parentNuclide()?p.parentNuclide()->symbol:string("none")) );

  // Eu-152-alone must not return a K40 peak.
  BOOST_CHECK_EQUAL( count_src( r_auto.observable_peaks, "K40" ), 0u );

  // Similar Eu-152 peak count between the two fits.
  const int dcount = (int)count_src(r_auto.observable_peaks,"Eu152")
                   - (int)count_src(r_joint.observable_peaks,"Eu152");
  BOOST_CHECK_MESSAGE( std::abs(dcount) <= 3, "Eu152 peak count auto-joint differs by " << dcount );

  // Major Eu-152 lines are present in both with similar area.
  const double key_lines[] = { 121.78, 344.28, 778.90, 964.08, 1408.01 };
  for( const double e : key_lines )
  {
    const PeakDef * const aa = find_peak_near( r_auto.observable_peaks, e, 3.0 );
    const PeakDef * const ja = find_peak_near( r_joint.observable_peaks, e, 3.0 );
    BOOST_CHECK_MESSAGE( aa, "auto missing Eu152 line " << e << " keV" );
    BOOST_CHECK_MESSAGE( ja, "joint missing Eu152 line " << e << " keV" );
    if( aa && ja && (ja->amplitude() > 0.0) )
    {
      const double ratio = aa->amplitude() / ja->amplitude();
      BOOST_CHECK_MESSAGE( (ratio > 0.75) && (ratio < 1.34),
        "Eu152 " << e << " keV area ratio auto/joint = " << ratio );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_r6_raw_interferer_transaction )
{
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  BOOST_REQUIRE( spec.foreground );
  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<shared_ptr<const PeakDef>> no_user_peaks;
  const FitPeaksForNuclides::PeakFitResult result = run_fit(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}),
      no_user_peaks, spec.isHPGe );
  BOOST_REQUIRE( result.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  for( const string &warning : result.warnings )
    BOOST_TEST_MESSAGE( warning );
  BOOST_REQUIRE( find_source_gamma(
      result.solution.m_peaks_without_back_sub, "K40", 1460.82, 0.5 ) );
  BOOST_CHECK( find_source_gamma(result.observable_peaks, "Eu152", 1408.01, 0.5) );
  BOOST_CHECK( find_source_gamma(result.observable_peaks, "Cs137", 661.657, 0.5) );
  BOOST_CHECK( !find_source_gamma(result.observable_peaks, "K40", 1460.82, 0.5) );

  std::set<string> nuisance_parents;
  for( const PeakDef &peak : result.solution.m_peaks_without_back_sub )
  {
    if( peak.parentNuclide() && (peak.parentNuclide()->symbol != "Eu152")
        && (peak.parentNuclide()->symbol != "Cs137") )
      nuisance_parents.insert( peak.parentNuclide()->symbol );
  }
  BOOST_CHECK_EQUAL( nuisance_parents.size(), 2u );

  // The caller must be able to turn the complete automatic R6 path off for controlled tuning and
  // for supplied-background workflows.  The default result above proves R6 is active; with the
  // opt-out, the same raw requested-source fit must remain source-only.
  const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_interferer_options
    = FitPeaksForNuclides::FitSrcPeaksOptions::DisableAutoInterfererFit;
  const FitPeaksForNuclides::PeakFitResult source_only = run_fit(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}),
      no_user_peaks, spec.isHPGe, no_interferer_options );
  BOOST_REQUIRE( source_only.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  verify_fit_result( source_only, no_user_peaks, spec.foreground, no_interferer_options );
  BOOST_CHECK( !find_source_gamma(
      source_only.solution.m_peaks_without_back_sub, "K40", 1460.82, 0.5 ) );
  BOOST_CHECK( find_source_gamma(source_only.observable_peaks, "Eu152", 1408.01, 0.5) );
  BOOST_CHECK( find_source_gamma(source_only.observable_peaks, "Cs137", 661.657, 0.5) );
  for( const PeakDef &peak : source_only.solution.m_peaks_without_back_sub )
  {
    BOOST_CHECK_MESSAGE( !peak.parentNuclide()
                         || (peak.parentNuclide()->symbol == "Eu152")
                         || (peak.parentNuclide()->symbol == "Cs137"),
                         "DisableAutoInterfererFit admitted unexpected nuisance "
                           << peak.parentNuclide()->symbol );
  }

  // An existing bystander represented as a floating peak at K40 must suppress the K40 nuisance,
  // avoiding two nearly-identical Gaussian components in the augmented solve.
  const shared_ptr<const PeakDef> k40_bystander = [&]() -> shared_ptr<const PeakDef> {
    for( const shared_ptr<const PeakDef> &peak : auto_peaks )
      if( peak && (std::fabs(peak->mean() - 1460.82) < 0.5) )
        return peak;
    return nullptr;
  }();
  BOOST_REQUIRE( k40_bystander );
  const FitPeaksForNuclides::PeakFitResult with_float = run_fit(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}),
      { k40_bystander }, spec.isHPGe,
      FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
  BOOST_REQUIRE( with_float.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_CHECK( !find_source_gamma(
      with_float.solution.m_peaks_without_back_sub, "K40", 1460.82, 0.5 ) );
  BOOST_CHECK( find_source_gamma(with_float.observable_peaks, "Eu152", 1408.01, 0.5) );

}


BOOST_AUTO_TEST_CASE( test_multisource_strong_norm_interferer_is_stable )
{
  // Raw trinitite has a strong K40 1460.82-keV line next to Eu152's weak 1457.64-keV line.
  // Without supplied background this exercises the active R6 nuisance path; with background it
  // verifies that the already-modeled line does not destabilize either requested-source order.
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  BOOST_REQUIRE( spec.foreground );
  BOOST_REQUIRE( spec.background );

  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );
  const vector<shared_ptr<const PeakDef>> no_user_peaks;

  const auto fit = [&]( const shared_ptr<const SpecUtils::Measurement> &background,
                        const vector<string> &source_names ) {
    return run_fit( spec.foreground, background, auto_peaks,
                    make_sources(source_names), no_user_peaks, spec.isHPGe );
  };

  const FitPeaksForNuclides::PeakFitResult raw_auto_ec
    = fit( nullptr, {"Eu152", "Cs137"} );
  const FitPeaksForNuclides::PeakFitResult raw_auto_ce
    = fit( nullptr, {"Cs137", "Eu152"} );
  const FitPeaksForNuclides::PeakFitResult raw_joint
    = fit( nullptr, {"K40", "Eu152", "Cs137"} );
  const FitPeaksForNuclides::PeakFitResult bg_auto_ec
    = fit( spec.background, {"Eu152", "Cs137"} );
  const FitPeaksForNuclides::PeakFitResult bg_auto_ce
    = fit( spec.background, {"Cs137", "Eu152"} );
  const FitPeaksForNuclides::PeakFitResult bg_joint
    = fit( spec.background, {"K40", "Eu152", "Cs137"} );

  const FitPeaksForNuclides::PeakFitResult * const results[] = {
    &raw_auto_ec, &raw_auto_ce, &raw_joint, &bg_auto_ec, &bg_auto_ce, &bg_joint
  };
  for( const FitPeaksForNuclides::PeakFitResult * const result : results )
  {
    BOOST_REQUIRE( result->status == RelActCalcAuto::RelActAutoSolution::Status::Success );
    BOOST_CHECK( !result->observable_peaks.empty() );
  }

  for( const string &warning : raw_auto_ec.warnings )
    BOOST_TEST_MESSAGE( "raw auto {Eu,Cs}: " << warning );

  const auto count_source = []( const vector<PeakDef> &peaks, const string &symbol ) {
    size_t count = 0;
    for( const PeakDef &peak : peaks )
      count += (peak.parentNuclide() && (peak.parentNuclide()->symbol == symbol)) ? 1u : 0u;
    return count;
  };

  const FitPeaksForNuclides::PeakFitResult * const automatic_results[] = {
    &raw_auto_ec, &raw_auto_ce, &bg_auto_ec, &bg_auto_ce
  };
  const double eu_anchors[] = { 344.28, 778.90, 1408.01 };
  for( const FitPeaksForNuclides::PeakFitResult * const result : automatic_results )
  {
    BOOST_CHECK_EQUAL( count_source(result->uncombined_fit_peaks, "K40"), 0u );
    BOOST_CHECK_EQUAL( count_source(result->fit_peaks, "K40"), 0u );
    BOOST_CHECK_EQUAL( count_source(result->observable_peaks, "K40"), 0u );
    BOOST_CHECK( find_source_gamma(result->observable_peaks, "Cs137", 661.657, 0.5) );
    for( const double energy : eu_anchors )
      BOOST_CHECK_MESSAGE( find_source_gamma(result->observable_peaks, "Eu152", energy, 0.5),
                           "Missing Eu152 anchor at " << energy << " keV" );
  }

  const auto compare_anchors = [&]( const FitPeaksForNuclides::PeakFitResult &lhs,
                                    const FitPeaksForNuclides::PeakFitResult &rhs ) {
    const PeakDef * const lhs_cs
      = find_source_gamma( lhs.uncombined_fit_peaks, "Cs137", 661.657, 0.5 );
    const PeakDef * const rhs_cs
      = find_source_gamma( rhs.uncombined_fit_peaks, "Cs137", 661.657, 0.5 );
    BOOST_REQUIRE( lhs_cs && rhs_cs );
    BOOST_CHECK( peak_areas_agree(*lhs_cs, *rhs_cs) );
    for( const double energy : eu_anchors )
    {
      const PeakDef * const lhs_peak
        = find_source_gamma( lhs.uncombined_fit_peaks, "Eu152", energy, 0.5 );
      const PeakDef * const rhs_peak
        = find_source_gamma( rhs.uncombined_fit_peaks, "Eu152", energy, 0.5 );
      BOOST_REQUIRE( lhs_peak && rhs_peak );
      BOOST_CHECK_MESSAGE( peak_areas_agree(*lhs_peak, *rhs_peak),
                           "Area mismatch at Eu152 " << energy << " keV" );
    }
  };

  compare_anchors( raw_auto_ec, raw_auto_ce );
  compare_anchors( raw_auto_ec, raw_joint );
  compare_anchors( bg_auto_ec, bg_auto_ce );
  compare_anchors( bg_auto_ec, bg_joint );

  // Prove the foreground-only arm actually retained a strong hidden K40 nuisance in the model.
  const PeakDef * const fitted_k40 = find_source_gamma(
      raw_auto_ec.solution.m_peaks_without_back_sub, "K40", 1460.82, 0.5 );
  const PeakDef * const fitted_k40_reversed = find_source_gamma(
      raw_auto_ce.solution.m_peaks_without_back_sub, "K40", 1460.82, 0.5 );
  const PeakDef * const fitted_eu1457 = find_source_gamma(
      raw_auto_ec.solution.m_peaks_without_back_sub, "Eu152", 1457.64, 0.5 );
  BOOST_REQUIRE( fitted_k40 );
  BOOST_REQUIRE( fitted_k40_reversed );
  BOOST_REQUIRE( fitted_eu1457 );
  BOOST_CHECK_GT( fitted_k40->peakAreaUncert(), 0.0 );
  BOOST_CHECK_GT( fitted_k40->peakArea() / fitted_k40->peakAreaUncert(), 10.0 );
  BOOST_CHECK_GT( fitted_k40->peakArea(), 5.0*fitted_eu1457->peakArea() );

  // The observable refit must not inflate the weak Eu152 line after its K40 nuisance is hidden.
  const PeakDef * const auto_eu1457
    = find_source_gamma( raw_auto_ec.observable_peaks, "Eu152", 1457.64, 0.5 );
  const PeakDef * const joint_eu1457
    = find_source_gamma( raw_joint.uncombined_fit_peaks, "Eu152", 1457.64, 0.5 );
  BOOST_REQUIRE( joint_eu1457 );
  BOOST_CHECK_MESSAGE( peak_areas_agree(*fitted_eu1457, *joint_eu1457),
                       "Solve-stage Eu152 1457-keV areas must agree before testing observable refit" );
  if( auto_eu1457 )
  {
    const double combined_uncert = std::hypot(
        std::max(0.0, auto_eu1457->peakAreaUncert()),
        std::max(0.0, joint_eu1457->peakAreaUncert()) );
    const double tolerance = std::max( 3.0*combined_uncert,
                                      0.20*std::fabs(joint_eu1457->peakArea()) );
    BOOST_CHECK_LE( auto_eu1457->peakArea() - joint_eu1457->peakArea(), tolerance );
  }
}


BOOST_AUTO_TEST_CASE( test_rescue_does_not_overfit_on_interferer )
{
#if( PERFORM_DEVELOPER_CHECKS )
  const LoadedSpectrum spec = load_test_data_spectrum(
    "trinitite_sample_b.n42", "trinitite_sample_b_background.n42" );
  const vector<shared_ptr<const PeakDef>> auto_peaks
    = run_auto_search( spec.foreground, spec.isHPGe );
  FitPeaksForNuclides::PeakFitForNuclideConfig config
    = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config(
        PeakFitUtils::CoarseResolutionType::High );
  config.manual_keep_significance_z = 15.0;
  config.auto_keep_significance_z = 15.0;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_ecal_options;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal;
  no_ecal_options |= FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal;

  struct RestoreRescue
  {
    ~RestoreRescue()
    {
      FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( true );
    }
  } restore_rescue;
  FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( true );
  const FitPeaksForNuclides::PeakFitResult enabled = run_fit_with_config(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}), {},
      spec.isHPGe, config, no_ecal_options );
  FitPeaksForNuclides::detail::set_bounded_rescue_enabled_for_test( false );
  const FitPeaksForNuclides::PeakFitResult disabled = run_fit_with_config(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}), {},
      spec.isHPGe, config, no_ecal_options );

  const auto did_rescue = []( const FitPeaksForNuclides::PeakFitResult &result ) {
    return std::any_of(std::begin(result.warnings), std::end(result.warnings),
      []( const string &warning ) {
        return warning.find("bounded fit-then-prune rescue") != string::npos;
      });
  };
  BOOST_REQUIRE( enabled.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_REQUIRE( disabled.status == RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_CHECK_MESSAGE( !did_rescue(enabled),
                       "A guarded marginal beside K40 was unexpectedly admitted" );
  BOOST_CHECK( find_source_gamma(enabled.observable_peaks, "Eu152", 1408.01, 0.5) );
  BOOST_CHECK( find_source_gamma(enabled.observable_peaks, "Cs137", 661.657, 0.5) );
  BOOST_CHECK( !find_source_gamma(enabled.observable_peaks, "K40", 1460.82, 0.5) );

  const PeakDef * const enabled_weak = find_source_gamma(
      enabled.solution.m_peaks_without_back_sub, "Eu152", 1457.64, 0.5 );
  const PeakDef * const disabled_weak = find_source_gamma(
      disabled.solution.m_peaks_without_back_sub, "Eu152", 1457.64, 0.5 );
  if( enabled_weak && disabled_weak )
  {
    const double combined_uncert = std::hypot(
        std::max(0.0, enabled_weak->peakAreaUncert()),
        std::max(0.0, disabled_weak->peakAreaUncert()) );
    const double tolerance = std::max( 3.0*combined_uncert,
                                      0.20*std::fabs(disabled_weak->peakArea()) );
    BOOST_CHECK_MESSAGE( enabled_weak->peakArea() - disabled_weak->peakArea() <= tolerance,
                         "R2 inflated Eu152 1457.6-keV area beside strong K40" );
  }
#endif
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

  // Step 1: Fit Eu-152 (with K40 - Eu-152 1457 keV overlaps K40's strong 1460 keV line)
  const vector<RelActCalcAuto::SrcVariant> sources_eu152 = make_sources( {"K40", "Eu152"} );
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


// ---------------------------------------------------------------------------
// Unit tests for the statistical detail:: helpers introduced by the
// dimensionless-parameter reformulation (local continuum estimate, adaptive
// ROI extent, clean-gap merge test, continuum-order selection).
// All use exact synthetic spectra (no Poisson noise) so assertions test the
// logic rather than noise robustness.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( StatisticalDetailHelpers )

namespace
{
  // Synthetic spectrum with per-keV continuum density given by `density`, plus optional
  // Gaussians described by (mean, sigma, total counts).
  std::shared_ptr<const SpecUtils::Measurement> make_synthetic_spectrum(
    const size_t nchannel,
    const float lower_energy,
    const float channel_width,
    const std::function<double(double)> &density,
    const std::vector<std::array<double,3>> &gaussians = {} )
  {
    auto cal = std::make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( nchannel, { lower_energy, channel_width }, {} );

    auto counts = std::make_shared<std::vector<float>>( nchannel, 0.0f );
    for( size_t i = 0; i < nchannel; ++i )
    {
      const double lo = lower_energy + i*channel_width;
      const double hi = lo + channel_width;
      const double mid = 0.5*(lo + hi);
      double val = density( mid ) * channel_width;
      for( const std::array<double,3> &g : gaussians )
      {
        const double t0 = (lo - g[0]) / (std::sqrt(2.0) * g[1]);
        const double t1 = (hi - g[0]) / (std::sqrt(2.0) * g[1]);
        val += g[2] * 0.5 * (std::erf(t1) - std::erf(t0));
      }
      (*counts)[i] = static_cast<float>( val );
    }

    auto meas = std::make_shared<SpecUtils::Measurement>();
    meas->set_gamma_counts( counts, 100.0f, 100.0f );
    meas->set_energy_calibration( cal );
    return meas;
  }//make_synthetic_spectrum
}//namespace


BOOST_AUTO_TEST_CASE( test_snip_joint_roi_boundary_shadow )
{
  using FitPeaksForNuclides::detail::GlobalContinuumEstimate;
  using FitPeaksForNuclides::detail::RoiBoundaryShadowGroup;
  using FitPeaksForNuclides::detail::RoiBoundaryShadowResult;
  using FitPeaksForNuclides::detail::optimize_roi_boundaries_shadow;

  const auto make_global = []( const shared_ptr<const SpecUtils::Measurement> &spectrum ) {
    GlobalContinuumEstimate global;
    global.snip = spectrum;
    global.foreground = spectrum;
    global.built = true;
    return global;
  };
  const vector<shared_ptr<const PeakDef>> no_unfit_peaks;

  // A wide, strongly non-polynomial baseline should be split into locally plausible intervals.
  const auto curved = make_synthetic_spectrum( 240, 0.0f, 1.0f,
    []( const double energy ) {
      return 8.0 + 5.0*std::sin(energy/13.0) + 0.00002*std::pow(energy - 120.0, 4.0);
    } );
  const auto fwhm4 = []( const double ){ return 4.0; };
  const vector<RoiBoundaryShadowGroup> wide_groups = {
    {40.0, 190.0, {60.0}}, {40.0, 190.0, {110.0}}, {40.0, 190.0, {165.0}}
  };
  const RoiBoundaryShadowResult curved_result = optimize_roi_boundaries_shadow(
      wide_groups, curved, make_global(curved), fwhm4, no_unfit_peaks, 50.0 );
  BOOST_REQUIRE( curved_result.valid );
  BOOST_CHECK_GE( curved_result.intervals.size(), 2u );
  BOOST_CHECK( std::isfinite(curved_result.legacy_total_score) );
  BOOST_CHECK_LT( curved_result.proposed_total_score, curved_result.legacy_total_score );
  for( size_t i = 0; i < curved_result.intervals.size(); ++i )
  {
    const auto &interval = curved_result.intervals[i];
    BOOST_CHECK_LE( interval.lower, interval.upper );
    BOOST_CHECK_EQUAL( interval.unmodeled_peak_conflicts, 0u );
    if( i )
      BOOST_CHECK_LE( curved_result.intervals[i-1].upper, interval.lower );
  }
  for( const RoiBoundaryShadowGroup &group : wide_groups )
  {
    const double energy = group.gamma_energies.front();
    const size_t covering = std::count_if( std::begin(curved_result.intervals),
      std::end(curved_result.intervals), [energy]( const auto &interval ) {
        return (energy >= interval.lower) && (energy <= interval.upper);
      } );
    BOOST_CHECK_EQUAL( covering, 1u );
  }

  // Overlapping source cores have no feasible boundary between them and must remain joint.
  const auto flat = make_synthetic_spectrum( 220, 0.0f, 1.0f,
      []( const double ){ return 10.0; } );
  const auto fwhm6 = []( const double ){ return 6.0; };
  const vector<RoiBoundaryShadowGroup> overlapping = {
    {85.0, 120.0, {100.0}}, {85.0, 120.0, {104.0}}
  };
  const RoiBoundaryShadowResult overlap_result = optimize_roi_boundaries_shadow(
      overlapping, flat, make_global(flat), fwhm6, no_unfit_peaks, 50.0 );
  BOOST_REQUIRE( overlap_result.valid );
  BOOST_CHECK_EQUAL( overlap_result.intervals.size(), 1u );

  // A real baseline discontinuity should select one of the production step families.
  const auto stepped = make_synthetic_spectrum( 220, 0.0f, 1.0f,
      []( const double energy ){ return (energy < 100.0) ? 5.0 : 22.0; } );
  const vector<RoiBoundaryShadowGroup> step_group = { {70.0, 130.0, {100.0}} };
  const auto raw_step_foreground = make_synthetic_spectrum( 220, 0.0f, 1.0f,
      []( const double energy ){
        return ((energy < 100.0) ? 8.0 : 25.0) + 0.005*energy;
      } );
  GlobalContinuumEstimate step_global = make_global( stepped );
  step_global.foreground = raw_step_foreground;
  const RoiBoundaryShadowResult step_result = optimize_roi_boundaries_shadow(
      step_group, raw_step_foreground, step_global, fwhm4, no_unfit_peaks, 50.0 );
  BOOST_REQUIRE( step_result.valid );
  BOOST_REQUIRE_EQUAL( step_result.intervals.size(), 1u );
  BOOST_CHECK( (step_result.intervals[0].continuum_type == PeakContinuum::OffsetType::FlatStep)
               || (step_result.intervals[0].continuum_type
                    == PeakContinuum::OffsetType::LinearStep) );

  GlobalContinuumEstimate invalid_global;
  const RoiBoundaryShadowResult fallback = optimize_roi_boundaries_shadow(
      step_group, stepped, invalid_global, fwhm4, no_unfit_peaks, 50.0 );
  BOOST_CHECK( !fallback.valid );
  BOOST_CHECK( !fallback.fallback_reason.empty() );

  // An unmodeled peak creates a one-FWHM exclusion gap.  The proposed intervals may end/start at
  // its edges but must never contain or bisect it.
  const shared_ptr<const PeakDef> unmodeled
    = make_shared<const PeakDef>( 100.0, 2.0, 1000.0 );
  const vector<RoiBoundaryShadowGroup> separated = {
    {45.0, 155.0, {60.0}}, {45.0, 155.0, {140.0}}
  };
  const RoiBoundaryShadowResult excluded = optimize_roi_boundaries_shadow(
      separated, flat, make_global(flat), fwhm4, {unmodeled}, 50.0 );
  BOOST_REQUIRE( excluded.valid );
  for( const auto &interval : excluded.intervals )
  {
    BOOST_CHECK( !((interval.lower < 104.0) && (interval.upper > 96.0)) );
    BOOST_CHECK_GE( interval.unmodeled_peak_conflicts, 1u );
    BOOST_CHECK( std::find_if( std::begin(interval.unmodeled_peak_energies),
      std::end(interval.unmodeled_peak_energies), []( const double energy ) {
        return std::fabs(energy - 100.0) < 0.1;
      } ) != std::end(interval.unmodeled_peak_energies) );
  }

  // Cores are clipped at the spectrum edge, while the configured maximum remains a hard fallback.
  const vector<RoiBoundaryShadowGroup> edge_group = { {0.0, 20.0, {1.0}} };
  const RoiBoundaryShadowResult edge_result = optimize_roi_boundaries_shadow(
      edge_group, flat, make_global(flat), fwhm4, no_unfit_peaks, 50.0 );
  BOOST_REQUIRE( edge_result.valid );
  BOOST_CHECK_GE( edge_result.intervals.front().lower, 0.0 );
  const RoiBoundaryShadowResult width_fallback = optimize_roi_boundaries_shadow(
      step_group, stepped, make_global(stepped), fwhm4, no_unfit_peaks, 1.0 );
  BOOST_CHECK( !width_fallback.valid );

  vector<RoiBoundaryShadowGroup> permuted = wide_groups;
  std::reverse( std::begin(permuted), std::end(permuted) );
  const RoiBoundaryShadowResult permuted_result = optimize_roi_boundaries_shadow(
      permuted, curved, make_global(curved), fwhm4, no_unfit_peaks, 50.0 );
  BOOST_REQUIRE( permuted_result.valid );
  BOOST_REQUIRE_EQUAL( permuted_result.intervals.size(), curved_result.intervals.size() );
  for( size_t i = 0; i < curved_result.intervals.size(); ++i )
  {
    BOOST_CHECK_SMALL( curved_result.intervals[i].lower
                       - permuted_result.intervals[i].lower, 0.05 );
    BOOST_CHECK_SMALL( curved_result.intervals[i].upper
                       - permuted_result.intervals[i].upper, 0.05 );
    BOOST_CHECK_EQUAL( curved_result.intervals[i].continuum_type,
                       permuted_result.intervals[i].continuum_type );
    BOOST_CHECK_SMALL( curved_result.intervals[i].normalized_continuum_mismatch
                       - permuted_result.intervals[i].normalized_continuum_mismatch, 1.0e-8 );
  }
}


BOOST_AUTO_TEST_CASE( test_estimate_local_continuum )
{
  using FitPeaksForNuclides::detail::LocalContinuumEstimate;
  using FitPeaksForNuclides::detail::estimate_local_continuum;

  // Flat continuum: 5 counts/keV over [0, 3000] keV, 1 keV channels
  const auto flat = make_synthetic_spectrum( 3000, 0.0f, 1.0f, []( double ){ return 5.0; } );

  const LocalContinuumEstimate flat_est = estimate_local_continuum( flat, 590.0, 610.0, 2.0, 0.5 );
  BOOST_REQUIRE( flat_est.valid );
  BOOST_CHECK_CLOSE( flat_est.integral( 595.0, 605.0 ), 50.0, 5.0 );

  // Sloped continuum: density falls linearly from 20 at 0 keV to 0 at 2000 keV
  const auto sloped = make_synthetic_spectrum( 2000, 0.0f, 1.0f,
      []( double e ){ return std::max( 0.0, 20.0 * (1.0 - e/2000.0) ); } );

  const LocalContinuumEstimate slope_est = estimate_local_continuum( sloped, 980.0, 1020.0, 2.0, 0.5 );
  BOOST_REQUIRE( slope_est.valid );
  // At 1000 keV density is 10/keV; over [990, 1010] expect ~200 counts
  BOOST_CHECK_CLOSE( slope_est.integral( 990.0, 1010.0 ), 200.0, 5.0 );

  // Degenerate inputs are flagged invalid rather than crashing
  const LocalContinuumEstimate bad = estimate_local_continuum( flat, 700.0, 650.0, 2.0, 0.5 );
  BOOST_CHECK( !bad.valid );
}//test_estimate_local_continuum


BOOST_AUTO_TEST_CASE( test_extend_roi_by_sidebands )
{
  using FitPeaksForNuclides::detail::AdaptiveExtentResult;
  using FitPeaksForNuclides::detail::extend_roi_by_sidebands;

  const double fwhm = 2.0;
  const auto fwhm_at = [fwhm]( double ){ return fwhm; };
  const std::vector<double> energies( 1, 600.0 );
  const std::vector<double> amps( 1, 1000.0 );

  // Case 1: flat continuum - extension should run to the cap on both sides
  const auto flat = make_synthetic_spectrum( 3000, 0.0f, 1.0f, []( double ){ return 5.0; },
                                             { {600.0, fwhm/2.355, 1000.0} } );
  const AdaptiveExtentResult full = extend_roi_by_sidebands(
      energies, amps, fwhm, flat, fwhm_at, {}, 1.5, 2.0, 5.0,
      PeakDef::SkewType::NoSkew, 0.0, 3000.0 );

  // Cap is 5 FWHM = 10 keV each side; block quantization can leave one block un-taken
  BOOST_CHECK_LT( full.lower, 600.0 - 0.8*5.0*fwhm );
  BOOST_CHECK_GT( full.upper, 600.0 + 0.8*5.0*fwhm );
  BOOST_CHECK_GT( full.sideband_lower_kev, 0.0 );
  BOOST_CHECK_GT( full.sideband_upper_kev, 0.0 );

  // Case 2: a large un-modeled structure at 610 keV must stop the high-side extension short,
  // while the clean low side still extends further out
  const auto bumped = make_synthetic_spectrum( 3000, 0.0f, 1.0f, []( double ){ return 5.0; },
      { {600.0, fwhm/2.355, 1000.0}, {610.0, fwhm/2.355, 5000.0} } );
  const AdaptiveExtentResult stopped = extend_roi_by_sidebands(
      energies, amps, fwhm, bumped, fwhm_at, {}, 1.5, 2.0, 8.0,
      PeakDef::SkewType::NoSkew, 0.0, 3000.0 );

  BOOST_CHECK_LT( stopped.upper, 609.0 );
  BOOST_CHECK_LT( stopped.lower, 600.0 - 0.8*8.0*fwhm );

  // Case 3: no usable spectrum - falls back to the core extent
  const AdaptiveExtentResult core_only = extend_roi_by_sidebands(
      energies, amps, fwhm, nullptr, fwhm_at, {}, 1.5, 2.0, 5.0,
      PeakDef::SkewType::NoSkew, 0.0, 3000.0 );
  BOOST_CHECK_CLOSE( core_only.lower, 600.0 - 1.5*fwhm, 1.0e-6 );
  BOOST_CHECK_CLOSE( core_only.upper, 600.0 + 1.5*fwhm, 1.0e-6 );
}//test_extend_roi_by_sidebands


BOOST_AUTO_TEST_CASE( test_find_clean_gap_between )
{
  using FitPeaksForNuclides::detail::find_clean_gap_between;

  const double fwhm = 2.0;
  const auto fwhm_at = [fwhm]( double ){ return fwhm; };
  const auto flat = make_synthetic_spectrum( 3000, 0.0f, 1.0f, []( double ){ return 5.0; } );

  double win_lo = 0.0, win_hi = 0.0;

  // Well-separated small peaks (10 FWHM apart): clean gap exists
  const std::vector<double> left_e( 1, 600.0 ), right_e( 1, 620.0 );
  const std::vector<double> small_amp( 1, 100.0 );
  BOOST_CHECK( find_clean_gap_between( left_e, small_amp, right_e, small_amp,
      600.0, 620.0, flat, fwhm_at, 2.0, 1.0, &win_lo, &win_hi ) );
  BOOST_CHECK_GT( win_lo, 600.0 - 1.0e-9 );
  BOOST_CHECK_LT( win_hi, 620.0 + 1.0e-9 );

  // Anchors closer than the required gap width: must merge (no room to anchor a continuum)
  const std::vector<double> close_right( 1, 601.5 );
  BOOST_CHECK( !find_clean_gap_between( left_e, small_amp, close_right, small_amp,
      600.0, 601.5, flat, fwhm_at, 2.0, 1.0, nullptr, nullptr ) );

  // At 3-FWHM separation the answer depends on amplitude vs continuum noise: small peaks leave
  // a clean anchoring window between them, but 1e7-count peaks put >> sqrt(continuum) of tail
  // into every candidate block (even the midpoint sits at only ~3.5 sigma) - must merge.
  const std::vector<double> mid_right( 1, 606.0 );
  BOOST_CHECK( find_clean_gap_between( left_e, small_amp, mid_right, small_amp,
      600.0, 606.0, flat, fwhm_at, 2.0, 1.0, nullptr, nullptr ) );
  const std::vector<double> huge_amp( 1, 1.0e7 );
  BOOST_CHECK( !find_clean_gap_between( left_e, huge_amp, mid_right, huge_amp,
      600.0, 606.0, flat, fwhm_at, 2.0, 1.0, nullptr, nullptr ) );

  // Zero/unknown amplitudes degrade to a pure gap-width test
  const std::vector<double> zero_amp( 1, 0.0 );
  BOOST_CHECK( find_clean_gap_between( left_e, zero_amp, right_e, zero_amp,
      600.0, 620.0, flat, fwhm_at, 2.0, 1.0, nullptr, nullptr ) );
}//test_find_clean_gap_between


BOOST_AUTO_TEST_CASE( test_central_automatic_roi_boundary_policy )
{
  using FitPeaksForNuclides::AutomaticRoiDecision;
  using FitPeaksForNuclides::detail::AutomaticRoiGroup;
  using FitPeaksForNuclides::detail::AutomaticRoiPolicyResult;
  using FitPeaksForNuclides::detail::AutomaticRoiPolicySettings;
  using FitPeaksForNuclides::detail::GlobalContinuumEstimate;
  using FitPeaksForNuclides::detail::evaluate_automatic_roi_boundary;

  const auto make_global = []( const shared_ptr<const SpecUtils::Measurement> &foreground,
                               const shared_ptr<const SpecUtils::Measurement> &snip ) {
    GlobalContinuumEstimate global;
    global.foreground = foreground;
    global.snip = snip;
    global.built = true;
    return global;
  };
  AutomaticRoiPolicySettings settings;
  settings.merge_tail_z = 2.0;
  settings.merge_clean_gap_fwhm = 1.0;
  settings.continuum_aicc_penalty = 2.0;
  settings.peak_core_num_fwhm = 1.0;
  settings.max_width_fwhm = 30.0;
  settings.stage = "deterministic unit test";

  // Piecewise-linear SNIP is exactly supportable by two child continua but not by any one
  // production family.  Strong, well-separated peaks therefore retain a boundary with no
  // statistically significant peak bridge.
  const double fwhm = 4.0;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto v_density = []( const double energy ) {
    return 8.0 + 0.20*std::fabs(energy - 650.0);
  };
  const auto clean_foreground = make_synthetic_spectrum( 900, 200.0f, 1.0f, v_density,
      { {600.0, fwhm/2.35482, 3000.0}, {700.0, fwhm/2.35482, 3000.0} } );
  const auto clean_snip = make_synthetic_spectrum( 900, 200.0f, 1.0f, v_density );
  GlobalContinuumEstimate clean_global = make_global( clean_foreground, clean_snip );
  AutomaticRoiGroup left{550.0, 660.0, {600.0}, {3000.0}, 1, false};
  AutomaticRoiGroup right{640.0, 750.0, {700.0}, {3000.0}, 1, false};
  const AutomaticRoiPolicyResult clean = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, &clean_global, fwhm_at, {}, settings );
  BOOST_CHECK( clean.decision == AutomaticRoiDecision::KeepSeparate );
  BOOST_CHECK( clean.diagnostic.sidebands_adequate );
  BOOST_CHECK_LE( clean.diagnostic.two_roi_aicc, clean.diagnostic.one_roi_aicc );
  BOOST_CHECK( clean.diagnostic.used_global_continuum );

  // High-statistics overlapping tails leave no defensible continuum window.
  const double overlap_fwhm = 6.0;
  const auto overlap_fwhm_at = [overlap_fwhm]( const double ){ return overlap_fwhm; };
  const auto overlap_foreground = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 10.0; },
      { {100.0, overlap_fwhm/2.35482, 1.0e7},
        {106.0, overlap_fwhm/2.35482, 1.0e7} } );
  const auto overlap_snip = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 10.0; } );
  GlobalContinuumEstimate overlap_global = make_global( overlap_foreground, overlap_snip );
  AutomaticRoiGroup overlap_left{80.0, 105.0, {100.0}, {1.0e7}, 1, false};
  AutomaticRoiGroup overlap_right{101.0, 130.0, {106.0}, {1.0e7}, 1, false};
  const AutomaticRoiPolicyResult overlap = evaluate_automatic_roi_boundary(
      overlap_left, overlap_right, overlap_foreground, &overlap_global,
      overlap_fwhm_at, {}, settings );
  BOOST_CHECK( overlap.decision == AutomaticRoiDecision::MergeInseparable );

  // An unmodeled core in the gap is neither a split point nor permission to merge through it.
  const shared_ptr<const PeakDef> unmodeled
    = make_shared<const PeakDef>( 650.0, fwhm/2.35482, 1000.0 );
  const AutomaticRoiPolicyResult blocked = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, &clean_global, fwhm_at, {unmodeled}, settings );
  BOOST_CHECK( blocked.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked );
  BOOST_CHECK( blocked.diagnostic.unmodeled_core_blocked );
  BOOST_CHECK_LE( blocked.exclusion_lower, 650.0 - fwhm );
  BOOST_CHECK_GE( blocked.exclusion_upper, 650.0 + fwhm );

  // Core occupancy is geometric, not a mean-only test: this peak mean is outside the anchor
  // interval, but its +1-FWHM core crosses the left edge of the proposed gap.
  const shared_ptr<const PeakDef> crossing_core
    = make_shared<const PeakDef>( 597.0, fwhm/2.35482, 1000.0 );
  const AutomaticRoiPolicyResult crossing_blocked = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, &clean_global, fwhm_at, {crossing_core}, settings );
  BOOST_CHECK( crossing_blocked.decision == AutomaticRoiDecision::UnmodeledFeatureBlocked );

  // A clean-looking window cannot authorize two ROIs without a valid common-domain continuum
  // comparison.  In particular, two unavailable sentinel scores must not compare as a tie.
  const AutomaticRoiPolicyResult no_aicc = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, nullptr, fwhm_at, {}, settings );
  BOOST_CHECK( (no_aicc.decision == AutomaticRoiDecision::MergeInseparable)
               || (no_aicc.decision == AutomaticRoiDecision::MergeInseparableWide) );
  BOOST_CHECK_EQUAL( no_aicc.diagnostic.two_roi_aicc,
                     std::numeric_limits<double>::max() );

  const auto invalid_fwhm = []( const double ) {
    return std::numeric_limits<double>::quiet_NaN();
  };
  const AutomaticRoiPolicyResult invalid_resolution = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, &clean_global, invalid_fwhm, {}, settings );
  BOOST_CHECK( invalid_resolution.decision == AutomaticRoiDecision::KeepSeparate );

  // Relabeling the same channels with a shifted current calibration preserves the decision and
  // boundary channel, proving the helper consumes one coherent working frame.
  const auto shifted_density = []( const double energy ) {
    return 8.0 + 0.20*std::fabs(energy - 660.0);
  };
  const auto shifted_foreground = make_synthetic_spectrum( 900, 210.0f, 1.0f,
      shifted_density,
      { {610.0, fwhm/2.35482, 3000.0}, {710.0, fwhm/2.35482, 3000.0} } );
  const auto shifted_snip = make_synthetic_spectrum( 900, 210.0f, 1.0f, shifted_density );
  GlobalContinuumEstimate shifted_global = make_global( shifted_foreground, shifted_snip );
  AutomaticRoiGroup shifted_left{560.0, 670.0, {610.0}, {3000.0}, 1, false};
  AutomaticRoiGroup shifted_right{650.0, 760.0, {710.0}, {3000.0}, 1, false};
  const AutomaticRoiPolicyResult shifted = evaluate_automatic_roi_boundary(
      shifted_left, shifted_right, shifted_foreground, &shifted_global,
      fwhm_at, {}, settings );
  BOOST_CHECK( shifted.decision == clean.decision );
  BOOST_CHECK_EQUAL(
      shifted_foreground->find_gamma_channel( shifted.boundary_energy ),
      clean_foreground->find_gamma_channel( clean.boundary_energy ) );
  BOOST_CHECK_EQUAL( shifted.diagnostic.boundary_channel,
                     clean.diagnostic.boundary_channel );
  BOOST_CHECK_EQUAL( shifted.diagnostic.calibration_num_channels,
                     clean.diagnostic.calibration_num_channels );
  BOOST_CHECK( shifted.diagnostic.used_global_continuum );

  // Protected mixed/fixed geometry is reported without mutating bounds or continuum metadata.
  left.protected_geometry = true;
  const AutomaticRoiGroup protected_before = left;
  const AutomaticRoiPolicyResult protected_result = evaluate_automatic_roi_boundary(
      left, right, clean_foreground, &clean_global, fwhm_at, {}, settings );
  BOOST_CHECK( protected_result.decision == AutomaticRoiDecision::ProtectedGeometry );
  BOOST_CHECK_EQUAL( left.lower, protected_before.lower );
  BOOST_CHECK_EQUAL( left.upper, protected_before.upper );
  BOOST_CHECK_EQUAL_COLLECTIONS( std::begin(left.peak_energies), std::end(left.peak_energies),
                                 std::begin(protected_before.peak_energies),
                                 std::end(protected_before.peak_energies) );

  // A causal eleven-group chain accumulates lineage and can exceed the soft onset only through
  // an explicit wide-inseparable outcome.
  settings.max_width_fwhm = 8.0;
  std::vector<std::array<double,3>> chain_gaussians;
  for( size_t group_index = 0; group_index < 11; ++group_index )
    chain_gaussians.push_back( {100.0 + 6.0*group_index,
                                overlap_fwhm/2.35482, 1.0e7} );
  const auto chain_foreground = make_synthetic_spectrum( 260, 0.0f, 1.0f,
      []( const double ){ return 10.0; }, chain_gaussians );
  const auto chain_snip = make_synthetic_spectrum( 260, 0.0f, 1.0f,
      []( const double ){ return 10.0; } );
  GlobalContinuumEstimate chain_global = make_global( chain_foreground, chain_snip );
  AutomaticRoiGroup chain{80.0, 105.0, {100.0}, {1.0e7}, 1, false};
  AutomaticRoiPolicyResult chain_decision;
  for( size_t group_index = 1; group_index < 11; ++group_index )
  {
    const double energy = 100.0 + 6.0*group_index;
    AutomaticRoiGroup next{energy - 5.0, energy + 5.0,
                           {energy}, {1.0e7}, 1, false};
    chain_decision = evaluate_automatic_roi_boundary(
        chain, next, chain_foreground, &chain_global, overlap_fwhm_at, {}, settings );
    BOOST_REQUIRE( (chain_decision.decision == AutomaticRoiDecision::MergeInseparable)
                   || (chain_decision.decision == AutomaticRoiDecision::MergeInseparableWide) );
    chain.upper = next.upper;
    chain.peak_energies.push_back( energy );
    chain.peak_areas.push_back( 1.0e7 );
    chain.joined_groups += 1;
  }
  BOOST_CHECK( chain_decision.decision == AutomaticRoiDecision::MergeInseparableWide );
  BOOST_CHECK_GT( chain_decision.diagnostic.width_pressure, 0.0 );
  BOOST_CHECK_EQUAL( std::string(FitPeaksForNuclides::automatic_roi_decision_name(
                         chain_decision.decision)),
                     "MergeInseparableWide" );
}//test_central_automatic_roi_boundary_policy


BOOST_AUTO_TEST_CASE( test_source_clean_challenger_model_selection_gates )
{
  using FitPeaksForNuclides::detail::should_accept_source_clean_challenger;
  using FitPeaksForNuclides::detail::should_try_source_clean_recovery;
  using FitPeaksForNuclides::detail::data_only_aicc;

  // A stable solution, or a single questionable lost line, never starts the transactional
  // challenger.  Two independent lost anchors do.
  BOOST_CHECK( !should_try_source_clean_recovery( 5, 5 ) );
  BOOST_CHECK( !should_try_source_clean_recovery( 5, 4 ) );
  BOOST_CHECK( should_try_source_clean_recovery( 5, 3 ) );

  BOOST_CHECK( should_accept_source_clean_challenger(
      true, 3, 5, 4, 4, 2.0, 1.5 ) );
  BOOST_CHECK( should_accept_source_clean_challenger(
      true, 0, 5, 0, 4, std::numeric_limits<double>::max(), 1.5 ) );
  BOOST_CHECK( !should_accept_source_clean_challenger(
      false, 3, 5, 4, 4, 2.0, 1.5 ) );  // failed solve rolls back
  BOOST_CHECK( !should_accept_source_clean_challenger(
      true, 3, 3, 4, 4, 2.0, 1.5 ) );   // no predicted-anchor recovery
  BOOST_CHECK( !should_accept_source_clean_challenger(
      true, 3, 5, 4, 3, 2.0, 1.5 ) );   // lost fitted source evidence
  BOOST_CHECK( !should_accept_source_clean_challenger(
      true, 3, 5, 4, 4, 2.0, 2.1 ) );   // data score did not improve
  BOOST_CHECK( !should_accept_source_clean_challenger( true, 3, 5, 4, 4, 2.0,
      std::numeric_limits<double>::quiet_NaN() ) );
  BOOST_CHECK( !should_accept_source_clean_challenger( true, 3, 5, 4, 4,
      std::numeric_limits<double>::max(), std::numeric_limits<double>::max() ) );

  BOOST_CHECK_CLOSE( data_only_aicc( 20.0, 20, 3, 2.0 ), 27.5, 1.0e-8 );
  BOOST_CHECK_EQUAL( data_only_aicc( 20.0, 4, 3, 2.0 ),
                     std::numeric_limits<double>::max() );
  BOOST_CHECK_EQUAL( data_only_aicc( std::numeric_limits<double>::quiet_NaN(),
                                    20, 3, 2.0 ),
                     std::numeric_limits<double>::max() );
}//test_source_clean_challenger_model_selection_gates


BOOST_AUTO_TEST_CASE( test_source_cluster_evidence_models )
{
  using FitPeaksForNuclides::detail::SourceClusterEvidenceDecision;
  using FitPeaksForNuclides::detail::SourceClusterEvidenceResult;
  using FitPeaksForNuclides::detail::evaluate_source_cluster_evidence;

  const double fwhm = 2.0;
  const double sigma = fwhm / 2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto flat = make_synthetic_spectrum(
      240, 0.0f, 1.0f, []( const double ){ return 20.0; } );

  // A predicted bridge with no measured peak is rejected by H0, despite its high predicted area.
  const SourceClusterEvidenceResult absent = evaluate_source_cluster_evidence(
      {100.0}, {1000.0}, 90.0, 112.0, flat, fwhm_at, {}, 3.0, 1.0, 2.0 );
  BOOST_CHECK( absent.decision == SourceClusterEvidenceDecision::RejectContinuumOnly );
  BOOST_CHECK( std::isfinite(absent.null_aicc) );

  // A measured source-shaped peak strongly favors the source-tied explanation.
  const auto source_data = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 20.0; }, { {100.0, sigma, 1200.0} } );
  std::shared_ptr<PeakDef> same_core = std::make_shared<PeakDef>( 100.0, sigma, 1200.0 );
  same_core->setAmplitudeUncert( 30.0 );
  const SourceClusterEvidenceResult present = evaluate_source_cluster_evidence(
      {100.0}, {1000.0}, 90.0, 112.0, source_data, fwhm_at, {same_core},
      3.0, 1.0, 2.0 );
  BOOST_CHECK( present.decision == SourceClusterEvidenceDecision::RetainSource );
  BOOST_CHECK_GT( present.source_likelihood_z, 3.0 );
  // A found peak in the same FWHM core confirms Hs; it is not a free contaminant challenger.
  BOOST_CHECK( !std::isfinite(present.free_feature_aicc) );

  // The local AICc comparison already charges Hs for its scale parameter.  Do not apply the
  // downstream observable-z gate a second time: a blended shoulder can be worth retaining even
  // when its isolated likelihood z is just below that later threshold.
  bool retained_subthreshold_aicc_winner = false;
  for( double area = 20.0; area <= 300.0; area += 10.0 )
  {
    const auto shoulder_data = make_synthetic_spectrum( 240, 0.0f, 1.0f,
        []( const double ){ return 20.0; }, { {100.0, sigma, area} } );
    const SourceClusterEvidenceResult shoulder = evaluate_source_cluster_evidence(
        {100.0}, {1000.0}, 90.0, 112.0, shoulder_data, fwhm_at, {},
        3.0, 1.0, 2.0 );
    retained_subthreshold_aicc_winner = retained_subthreshold_aicc_winner
        || ((shoulder.decision == SourceClusterEvidenceDecision::RetainSource)
            && (shoulder.source_likelihood_z < 3.0));
  }
  BOOST_CHECK( retained_subthreshold_aicc_winner );

  // Hs is an exact fixed-ratio mixture, not one moment-matched Gaussian.  A found peak on either
  // requested line is part of Hs provenance and must not be offered back as a free contaminant.
  const auto mixture_data = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 20.0; },
      { {99.0, sigma, 400.0}, {103.0, sigma, 1200.0} } );
  std::shared_ptr<PeakDef> peripheral_source
    = std::make_shared<PeakDef>( 99.0, sigma, 400.0 );
  peripheral_source->setAmplitudeUncert( 20.0 );
  const SourceClusterEvidenceResult mixture = evaluate_source_cluster_evidence(
      {99.0, 103.0}, {1.0, 3.0}, 92.0, 110.0, mixture_data, fwhm_at,
      {peripheral_source}, 3.0, 1.0, 2.0 );
  BOOST_CHECK( mixture.decision == SourceClusterEvidenceDecision::RetainSource );
  BOOST_CHECK( !std::isfinite(mixture.free_feature_aicc) );

  // When a distinct strong feature explains the local data better than the requested-source
  // shape, Hf wins.  The smaller source peak still makes Hs significant relative to H0.
  const auto contaminated = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 20.0; },
      { {100.0, sigma, 1500.0}, {108.0, sigma, 2200.0} } );
  std::shared_ptr<PeakDef> free_peak = std::make_shared<PeakDef>( 108.0, sigma, 2200.0 );
  free_peak->setAmplitudeUncert( 35.0 );
  const SourceClusterEvidenceResult free_wins = evaluate_source_cluster_evidence(
      {100.0}, {1500.0}, 90.0, 115.0, contaminated, fwhm_at, {free_peak},
      3.0, 1.0, 2.0 );
  BOOST_CHECK( free_wins.decision == SourceClusterEvidenceDecision::RejectFreeFeature );
  BOOST_CHECK_LT( free_wins.free_feature_aicc, free_wins.source_aicc );

  // Multiple distinct found features compete jointly and pay for every fitted amplitude.
  const auto two_contaminants = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 20.0; },
      { {100.0, sigma, 300.0}, {106.0, sigma, 1300.0}, {111.0, sigma, 1100.0} } );
  std::shared_ptr<PeakDef> free_one = std::make_shared<PeakDef>( 106.0, sigma, 1300.0 );
  std::shared_ptr<PeakDef> free_two = std::make_shared<PeakDef>( 111.0, sigma, 1100.0 );
  free_one->setAmplitudeUncert( 30.0 );
  free_two->setAmplitudeUncert( 30.0 );
  const SourceClusterEvidenceResult joint_free = evaluate_source_cluster_evidence(
      {100.0}, {300.0}, 92.0, 116.0, two_contaminants, fwhm_at,
      {free_one, free_two}, 3.0, 1.0, 2.0 );
  BOOST_CHECK( joint_free.decision == SourceClusterEvidenceDecision::RejectFreeFeature );
  BOOST_CHECK_LT( joint_free.free_feature_aicc, joint_free.source_aicc );

  const SourceClusterEvidenceResult invalid = evaluate_source_cluster_evidence(
      {}, {}, 90.0, 112.0, flat, fwhm_at, {}, 3.0, 1.0, 2.0 );
  BOOST_CHECK( invalid.decision == SourceClusterEvidenceDecision::InsufficientEvidence );
}//test_source_cluster_evidence_models


namespace
{
  using FitPeaksForNuclides::detail::RoiAtom;
  using FitPeaksForNuclides::detail::RoiAtomKind;
  using FitPeaksForNuclides::detail::AutomaticRoiComponent;
  using FitPeaksForNuclides::detail::AutomaticRoiPartitionResult;
  using FitPeaksForNuclides::detail::AutomaticRoiPartitionOutcome;
  using FitPeaksForNuclides::detail::AutomaticRoiReconcileResult;
  using FitPeaksForNuclides::detail::AutomaticRoiPolicySettings;
  using FitPeaksForNuclides::detail::AutomaticRoiPartitionConstraints;
  using FitPeaksForNuclides::detail::GlobalContinuumEstimate;

  RoiAtom mk_atom( const double energy, const double area,
                   const RoiAtomKind kind = RoiAtomKind::ModeledGamma )
  {
    RoiAtom a;
    a.id = FitPeaksForNuclides::detail::next_roi_atom_id();
    a.energy = energy;
    a.area = area;
    a.kind = kind;
    return a;
  }

  AutomaticRoiComponent mk_component( const double lower, const double upper,
      const std::shared_ptr<const SpecUtils::Measurement> &fg,
      std::vector<RoiAtom> atoms, const bool protected_geometry = false )
  {
    AutomaticRoiComponent c;
    c.lower = lower;
    c.upper = upper;
    c.first_channel = fg->find_gamma_channel( static_cast<float>(lower) );
    c.last_channel = fg->find_gamma_channel( static_cast<float>(upper) );
    std::sort( std::begin(atoms), std::end(atoms),
      []( const RoiAtom &a, const RoiAtom &b ){ return a.energy < b.energy; } );
    c.atoms = std::move( atoms );
    c.protected_geometry = protected_geometry;
    return c;
  }

  std::multiset<uint64_t> atom_ids( const std::vector<AutomaticRoiComponent> &comps,
                                    const std::vector<RoiAtom> &orphans = {} )
  {
    std::multiset<uint64_t> ids;
    for( const AutomaticRoiComponent &c : comps )
      for( const RoiAtom &a : c.atoms )
        ids.insert( a.id );
    for( const RoiAtom &a : orphans )
      ids.insert( a.id );
    return ids;
  }

  // The central invariant: the atom-ID multiset is preserved exactly-once, components are
  // channel-disjoint, and every atom energy lies within its owning component.
  void check_exact_once( const std::vector<AutomaticRoiComponent> &before,
      const std::vector<AutomaticRoiComponent> &after, const std::vector<RoiAtom> &orphans )
  {
    BOOST_CHECK( atom_ids(before) == atom_ids(after, orphans) );
    std::set<uint64_t> seen;
    for( const AutomaticRoiComponent &c : after )
      for( const RoiAtom &a : c.atoms )
        BOOST_CHECK_MESSAGE( seen.insert(a.id).second, "atom owned by two components" );
    for( size_t i = 1; i < after.size(); ++i )
      BOOST_CHECK_GT( after[i].first_channel, after[i-1].last_channel );
    for( const AutomaticRoiComponent &c : after )
      for( const RoiAtom &a : c.atoms )
      {
        BOOST_CHECK_GE( a.energy, c.lower - 1.0e-6 );
        BOOST_CHECK_LE( a.energy, c.upper + 1.0e-6 );
      }
  }

  GlobalContinuumEstimate mk_global( const std::shared_ptr<const SpecUtils::Measurement> &fg,
                                     const std::shared_ptr<const SpecUtils::Measurement> &snip )
  {
    GlobalContinuumEstimate g;
    g.foreground = fg;
    g.snip = snip;
    g.built = true;
    return g;
  }

  AutomaticRoiPolicySettings default_policy_settings()
  {
    AutomaticRoiPolicySettings s;
    s.merge_tail_z = 2.0;
    s.merge_clean_gap_fwhm = 1.0;
    s.continuum_aicc_penalty = 2.0;
    s.peak_core_num_fwhm = 1.0;
    s.max_width_fwhm = 30.0;
    s.stage = "partition unit test";
    return s;
  }
}//namespace


BOOST_AUTO_TEST_CASE( test_whole_component_measured_partition )
{
  using FitPeaksForNuclides::AutomaticRoiDecision;
  using FitPeaksForNuclides::detail::AutomaticRoiComponentPartitionResult;
  using FitPeaksForNuclides::detail::partition_overwide_automatic_component;

  const double fwhm = 2.0;
  const double sigma = fwhm / 2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto foreground = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 15.0; },
      { {100.0, sigma, 1500.0}, {120.0, sigma, 1800.0} } );
  const AutomaticRoiComponent wide = mk_component( 90.0, 130.0, foreground,
      { mk_atom(100.0, 1500.0), mk_atom(120.0, 1800.0) } );

  AutomaticRoiPolicySettings settings = default_policy_settings();
  settings.max_width_fwhm = 8.0;
  AutomaticRoiPartitionConstraints constraints;
  constraints.lowest_energy = 0.0;
  constraints.highest_energy = 240.0;
  constraints.left_barrier = -std::numeric_limits<double>::infinity();
  constraints.min_width_fwhm = 0.0;
  constraints.peak_core_num_fwhm = 1.0;

  const AutomaticRoiComponentPartitionResult split
    = partition_overwide_automatic_component( {wide}, foreground, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( split.valid );
  BOOST_REQUIRE( split.changed );
  BOOST_CHECK( split.diagnostic.decision == AutomaticRoiDecision::KeepSeparate );
  BOOST_CHECK_LT( split.diagnostic.two_roi_aicc, split.diagnostic.one_roi_aicc );
  BOOST_REQUIRE_EQUAL( split.components.size(), 2u );
  check_exact_once( {wide}, split.components, {} );
  BOOST_CHECK_GT( split.components[1].first_channel,
                  split.components[0].last_channel );

  // The minimum-gap rail applies before AICc selection: it can keep a dense or only modestly
  // resolved group intact even when separate continua would win the flexible measured-data fit.
  settings.minimum_partition_gap_fwhm = 12.0;
  const AutomaticRoiComponentPartitionResult minimum_gap_blocked
    = partition_overwide_automatic_component( {wide}, foreground, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( minimum_gap_blocked.valid );
  BOOST_CHECK( !minimum_gap_blocked.changed );
  BOOST_CHECK( minimum_gap_blocked.diagnostic.decision
               == AutomaticRoiDecision::MergeInseparableWide );
  check_exact_once( {wide}, minimum_gap_blocked.components, {} );
  // A visually clean continuum window is distinct evidence from modeled-core spacing.  The
  // override is explicit: it may admit this sparse valley without weakening the hard rail for
  // dense multiplets where the clean-window test fails.
  settings.allow_clean_gap_partition_override = true;
  const AutomaticRoiComponentPartitionResult clean_gap_override
    = partition_overwide_automatic_component( {wide}, foreground, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( clean_gap_override.valid );
  BOOST_CHECK( clean_gap_override.changed );
  check_exact_once( {wide}, clean_gap_override.components, {} );
  settings.allow_clean_gap_partition_override = false;
  settings.minimum_partition_gap_fwhm = 0.0;

  // A clean core gap is an explicit, separately tunable escape hatch for a final fitted ROI.
  // A deliberately high continuum-complexity penalty makes the otherwise valid two-ROI model
  // lose AICc; the default must retain that union.  Enabling the configured gap admits the same
  // atom-safe partition and labels it truthfully, rather than silently presenting it as an AICc
  // win.
  const auto close_peak_data = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 15.0; },
      { {100.0, sigma, 50.0}, {108.0, sigma, 60.0} } );
  const AutomaticRoiComponent unsupported_wide = mk_component( 94.0, 114.0,
      close_peak_data, { mk_atom(100.0, 50.0), mk_atom(108.0, 60.0) } );
  settings.continuum_aicc_penalty = 25.0;
  const AutomaticRoiComponentPartitionResult unsupported_default
    = partition_overwide_automatic_component( {unsupported_wide}, close_peak_data, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( unsupported_default.valid );
  BOOST_CHECK( !unsupported_default.changed );
  BOOST_CHECK_GE( unsupported_default.diagnostic.two_roi_aicc,
                  unsupported_default.diagnostic.one_roi_aicc );
  settings.force_partition_gap_fwhm = 2.0;
  const AutomaticRoiComponentPartitionResult unsupported_forced
    = partition_overwide_automatic_component( {unsupported_wide}, close_peak_data, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( unsupported_forced.valid );
  BOOST_CHECK( unsupported_forced.changed );
  BOOST_CHECK( unsupported_forced.diagnostic.decision == AutomaticRoiDecision::KeepSeparate );
  BOOST_CHECK( unsupported_forced.diagnostic.reason.find( "configured clean core gap" )
               != std::string::npos );
  check_exact_once( {unsupported_wide}, unsupported_forced.components, {} );
  settings.force_partition_gap_fwhm = 0.0;
  settings.continuum_aicc_penalty = 2.0;

  // A three-anchor component still produces the explicitly-scored two-continuum challenger, not
  // an unreported multi-segment DP result, and owns every atom exactly once.
  const auto three_peak_data = make_synthetic_spectrum( 240, 0.0f, 1.0f,
      []( const double ){ return 15.0; },
      { {96.0, sigma, 900.0}, {108.0, sigma, 700.0}, {124.0, sigma, 1400.0} } );
  const AutomaticRoiComponent three_anchor = mk_component( 88.0, 132.0, three_peak_data,
      { mk_atom(96.0, 900.0), mk_atom(108.0, 700.0), mk_atom(124.0, 1400.0) } );
  const AutomaticRoiComponentPartitionResult three_split
    = partition_overwide_automatic_component( {three_anchor}, three_peak_data, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( three_split.valid );
  BOOST_REQUIRE( three_split.changed );
  BOOST_REQUIRE_EQUAL( three_split.components.size(), 2u );
  check_exact_once( {three_anchor}, three_split.components, {} );

  // Below the same configured soft-width onset, the incumbent passes through unchanged.
  settings.max_width_fwhm = 30.0;
  const AutomaticRoiComponentPartitionResult nonwide
    = partition_overwide_automatic_component( {wide}, foreground, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( nonwide.valid );
  BOOST_CHECK( !nonwide.changed );
  BOOST_REQUIRE_EQUAL( nonwide.components.size(), 1u );
  check_exact_once( {wide}, nonwide.components, {} );

  // FWHM-connected cores are explicitly inseparable; retaining the union is not a silent
  // fallback and the atom ledger remains exact-once.
  settings.max_width_fwhm = 1.0;
  const AutomaticRoiComponent connected = mk_component( 94.0, 108.0, foreground,
      { mk_atom(100.0, 900.0), mk_atom(102.0, 800.0) } );
  const AutomaticRoiComponentPartitionResult inseparable
    = partition_overwide_automatic_component( {connected}, foreground, fwhm_at,
        {}, settings, constraints );
  BOOST_REQUIRE( inseparable.valid );
  BOOST_CHECK( !inseparable.changed );
  BOOST_CHECK( inseparable.diagnostic.decision
               == AutomaticRoiDecision::MergeInseparableWide );
  check_exact_once( {connected}, inseparable.components, {} );

  // Invalid scoring input cannot partially materialize children.
  const AutomaticRoiComponentPartitionResult invalid
    = partition_overwide_automatic_component( {wide}, nullptr, fwhm_at,
        {}, settings, constraints );
  BOOST_CHECK( !invalid.valid );
  BOOST_CHECK( !invalid.changed );
}//test_whole_component_measured_partition


// Test A: a non-dominant outer gamma in a multi-line left group (the Am241 failure shape) is kept
// in the child covering it - never clipped-and-dropped - and the split still happens.
BOOST_AUTO_TEST_CASE( test_partition_spatial_atom_reassignment_multiline )
{
  using FitPeaksForNuclides::detail::partition_automatic_roi_pair;
  const double fwhm = 4.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  // A piecewise-linear (V) continuum is better fit by two child continua than one, so the oracle
  // returns KeepSeparate (a flat baseline would be cheaper as a single continuum -> Merge).  Two
  // well-separated groups mirror the validated clean-valley geometry.
  const auto density = []( const double e ){ return 8.0 + 0.20*std::fabs(e-650.0); };
  const auto fg = make_synthetic_spectrum( 900, 200.0f, 1.0f, density,
      { {600.0, sigma, 3000.0}, {620.0, sigma, 1200.0}, {700.0, sigma, 3000.0} } );
  const auto snip = make_synthetic_spectrum( 900, 200.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  // Left group carries a NON-DOMINANT outer line at 620 keV (the Am241-style failure shape).
  const RoiAtom a600 = mk_atom( 600.0, 3000.0 );
  const RoiAtom a620 = mk_atom( 620.0, 1200.0 );
  const RoiAtom a700 = mk_atom( 700.0, 3000.0 );
  const AutomaticRoiComponent left = mk_component( 588.0, 632.0, fg, { a600, a620 } );
  const AutomaticRoiComponent right = mk_component( 688.0, 752.0, fg, { a700 } );

  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 200.0; cons.highest_energy = 1100.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  const AutomaticRoiPartitionResult res = partition_automatic_roi_pair( left, right, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons );

  BOOST_CHECK( res.outcome == AutomaticRoiPartitionOutcome::KeptSeparate );
  BOOST_REQUIRE_EQUAL( res.components.size(), 2u );
  check_exact_once( { left, right }, res.components, res.orphaned_atoms );
  // The non-dominant 620 keV line stays with 600 in the left child (boundary sits above 620),
  // and its core is fully covered - it is neither clipped nor dropped.
  const bool a620_in_left = std::any_of( std::begin(res.components[0].atoms),
      std::end(res.components[0].atoms), [&]( const RoiAtom &a ){ return a.id == a620.id; } );
  BOOST_CHECK( a620_in_left );
  BOOST_CHECK_GE( res.components[0].upper, 620.0 + fwhm );
}


// Test B: the chosen boundary lands on a channel edge with a one-channel gap - children never
// share a channel and their bounds equal exact channel edges.
BOOST_AUTO_TEST_CASE( test_partition_channel_rounding_boundary )
{
  using FitPeaksForNuclides::detail::partition_automatic_roi_pair;
  const double fwhm = 3.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double e ){ return 10.0 + 0.30*std::fabs(e-120.0); };
  const auto fg = make_synthetic_spectrum( 300, 0.0f, 1.0f, density,
      { {100.0, sigma, 6000.0}, {140.0, sigma, 6000.0} } );
  const auto snip = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  const AutomaticRoiComponent left = mk_component( 88.0, 118.0, fg, { mk_atom(100.0,6000.0) } );
  const AutomaticRoiComponent right = mk_component( 118.0, 154.0, fg, { mk_atom(140.0,6000.0) } );
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 300.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  const AutomaticRoiPartitionResult res = partition_automatic_roi_pair( left, right, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons );
  if( res.outcome == AutomaticRoiPartitionOutcome::KeptSeparate )
  {
    BOOST_REQUIRE_EQUAL( res.components.size(), 2u );
    check_exact_once( { left, right }, res.components, res.orphaned_atoms );
    // exactly one excluded gap channel between the children
    BOOST_CHECK_EQUAL( res.components[1].first_channel, res.components[0].last_channel + 2 );
    BOOST_CHECK_CLOSE( res.components[0].upper,
        fg->gamma_channel_upper( res.components[0].last_channel ), 1.0e-4 );
    BOOST_CHECK_CLOSE( res.components[1].lower,
        fg->gamma_channel_lower( res.components[1].first_channel ), 1.0e-4 );
  }
}


// Test C: when the atom cores cannot be separated by any channel, the pair MERGES - it never
// drops a side.  Exercised by making the partition core wider than the oracle's separation core.
BOOST_AUTO_TEST_CASE( test_partition_no_core_safe_channel_merges )
{
  using FitPeaksForNuclides::detail::partition_automatic_roi_pair;
  const double fwhm = 2.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double e ){ return 8.0 + 0.30*std::fabs(e-102.0); };
  const auto fg = make_synthetic_spectrum( 200, 50.0f, 1.0f, density,
      { {100.0, sigma, 5000.0}, {104.0, sigma, 5000.0} } );
  const auto snip = make_synthetic_spectrum( 200, 50.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  AutomaticRoiPolicySettings settings = default_policy_settings();
  settings.peak_core_num_fwhm = 0.25;   // oracle sees no significant bridge between the small cores
  const AutomaticRoiComponent left = mk_component( 94.0, 102.0, fg, { mk_atom(100.0,5000.0) } );
  const AutomaticRoiComponent right = mk_component( 102.0, 110.0, fg, { mk_atom(104.0,5000.0) } );
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 50.0; cons.highest_energy = 250.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.5;   // partition cores overlap -> no gap

  const AutomaticRoiPartitionResult res = partition_automatic_roi_pair( left, right, fg, &global,
      fwhm_at, {}, settings, cons );
  BOOST_CHECK( res.outcome == AutomaticRoiPartitionOutcome::Merged );
  BOOST_REQUIRE_EQUAL( res.components.size(), 1u );
  check_exact_once( { left, right }, res.components, res.orphaned_atoms );
  BOOST_CHECK_EQUAL( res.components[0].atoms.size(), 2u );
}


// Test D: a min-width child pinned against the spectrum edge cannot widen; the pair merges (or
// finds an alternate) but never drops an atom, and never exceeds the spectrum extent.
BOOST_AUTO_TEST_CASE( test_partition_underwidth_child_at_spectrum_edge )
{
  using FitPeaksForNuclides::detail::partition_automatic_roi_pair;
  const double fwhm = 1.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 8.0; };
  const auto fg = make_synthetic_spectrum( 152, 0.0f, 1.0f, density,
      { {146.0, sigma, 4000.0}, {150.0, sigma, 4000.0} } );
  const auto snip = make_synthetic_spectrum( 152, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  const double highest = fg->gamma_channel_upper( fg->num_gamma_channels() - 1 );
  const AutomaticRoiComponent left = mk_component( 142.0, 148.0, fg, { mk_atom(146.0,4000.0) } );
  const AutomaticRoiComponent right = mk_component( 148.0, 151.0, fg, { mk_atom(150.0,4000.0) } );
  AutomaticRoiPolicySettings settings = default_policy_settings();
  settings.peak_core_num_fwhm = 0.5;
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = highest;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 6.0; cons.peak_core_num_fwhm = 0.5;

  const AutomaticRoiPartitionResult res = partition_automatic_roi_pair( left, right, fg, &global,
      fwhm_at, {}, settings, cons );
  BOOST_CHECK( res.outcome != AutomaticRoiPartitionOutcome::Infeasible );
  check_exact_once( { left, right }, res.components, res.orphaned_atoms );
  for( const AutomaticRoiComponent &c : res.components )
    BOOST_CHECK_LE( c.upper, highest + 1.0e-6 );
}


// Test E: an unmodeled-feature exclusion band that would cut an admitted atom core is not carved
// through; the partition either finds a core-safe boundary or merges, never splitting the core.
BOOST_AUTO_TEST_CASE( test_partition_exclusion_band_clipping_core_declined )
{
  using FitPeaksForNuclides::detail::partition_automatic_roi_pair;
  const double fwhm = 4.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double e ){ return 8.0 + 0.20*std::fabs(e-650.0); };
  const auto fg = make_synthetic_spectrum( 900, 200.0f, 1.0f, density,
      { {600.0, sigma, 3000.0}, {648.0, sigma, 1200.0}, {700.0, sigma, 3000.0} } );
  const auto snip = make_synthetic_spectrum( 900, 200.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  const RoiAtom a648 = mk_atom( 648.0, 1200.0 );
  const AutomaticRoiComponent left = mk_component( 552.0, 660.0, fg,
      { mk_atom(600.0,3000.0), a648 } );
  const AutomaticRoiComponent right = mk_component( 640.0, 752.0, fg, { mk_atom(700.0,3000.0) } );
  // Unmodeled peak at 650 keV sits in the anchor gap; its core overlaps 648's core.
  const std::shared_ptr<const PeakDef> unmodeled
      = std::make_shared<const PeakDef>( 650.0, sigma, 900.0 );

  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 200.0; cons.highest_energy = 1100.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  const AutomaticRoiPartitionResult res = partition_automatic_roi_pair( left, right, fg, &global,
      fwhm_at, { unmodeled }, default_policy_settings(), cons );
  BOOST_CHECK( res.outcome != AutomaticRoiPartitionOutcome::Infeasible );
  check_exact_once( { left, right }, res.components, res.orphaned_atoms );
  // 648's core must lie wholly inside whichever component owns it (validator guarantees this).
  for( const AutomaticRoiComponent &c : res.components )
    for( const RoiAtom &a : c.atoms )
      if( a.id == a648.id )
      {
        BOOST_CHECK_LE( c.lower, 648.0 - fwhm + 1.0e-6 );
        BOOST_CHECK_GE( c.upper, 648.0 + fwhm - 1.0e-6 );
      }
}


// Test F: protected geometry is pinned - its bounds/metadata are bit-identical afterward, an atom
// whose core falls inside it is booked to it, and an atom straddling its edge is orphaned (never
// silently dropped from the ledger).
BOOST_AUTO_TEST_CASE( test_reconcile_protected_geometry_pins )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  const double fwhm = 2.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 300, 0.0f, 1.0f, density,
      { {100.0, sigma, 4000.0}, {150.0, sigma, 4000.0}, {158.0, sigma, 2000.0} } );
  const auto snip = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  // Protected ROI [148,156]; an automatic neighbor overlaps it with one atom inside (150) and one
  // atom straddling its upper edge (155.5 with core reaching past 156).
  AutomaticRoiComponent prot = mk_component( 148.0, 156.0, fg,
      { mk_atom(150.0, 4000.0) }, /*protected*/true );
  const RoiAtom a150b = mk_atom( 151.0, 500.0 );      // core [150,152] inside protected
  const RoiAtom a_straddle = mk_atom( 155.5, 600.0 ); // core [154.5,156.5] crosses the pin
  AutomaticRoiComponent autom = mk_component( 149.0, 162.0, fg, { a150b, a_straddle } );
  const AutomaticRoiComponent left = mk_component( 92.0, 108.0, fg, { mk_atom(100.0,4000.0) } );

  const double prot_lower = prot.lower, prot_upper = prot.upper;
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 300.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  std::vector<AutomaticRoiComponent> components = { left, prot, autom };
  const std::vector<AutomaticRoiComponent> before = components;
  const AutomaticRoiReconcileResult res = reconcile_automatic_components( components, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons, nullptr );

  BOOST_CHECK( res.valid );
  check_exact_once( before, res.components, res.orphaned_atoms );
  // The protected component survives with unchanged bounds.
  const auto prot_it = std::find_if( std::begin(res.components), std::end(res.components),
      []( const AutomaticRoiComponent &c ){ return c.protected_geometry; } );
  BOOST_REQUIRE( prot_it != std::end(res.components) );
  BOOST_CHECK_CLOSE( prot_it->lower, prot_lower, 1.0e-6 );
  BOOST_CHECK_CLOSE( prot_it->upper, prot_upper, 1.0e-6 );
  // The straddling atom is accounted for as an orphan.
  const bool straddle_orphaned = std::any_of( std::begin(res.orphaned_atoms),
      std::end(res.orphaned_atoms), [&]( const RoiAtom &a ){ return a.id == a_straddle.id; } );
  BOOST_CHECK( straddle_orphaned );
}


// Test G: overlapping input ROIs whose atoms share the overlap band are reconciled to exactly-once
// ownership - the direct regression for the flat-list double-claim path.
BOOST_AUTO_TEST_CASE( test_reconcile_overlapping_ownership_exact_once )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  using FitPeaksForNuclides::detail::assign_atoms_to_disjoint_rois;
  const double fwhm = 2.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 300, 0.0f, 1.0f, density,
      { {100.0, sigma, 5000.0}, {112.0, sigma, 4000.0}, {124.0, sigma, 5000.0} } );
  const auto snip = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  AutomaticRoiComponent left = mk_component( 92.0, 118.0, fg,
      { mk_atom(100.0,5000.0), mk_atom(112.0,4000.0) } );
  AutomaticRoiComponent right = mk_component( 110.0, 132.0, fg, { mk_atom(124.0,5000.0) } );
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 300.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  std::vector<AutomaticRoiComponent> components = { left, right };
  const std::vector<AutomaticRoiComponent> before = components;
  const AutomaticRoiReconcileResult res = reconcile_automatic_components( components, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons, nullptr );
  BOOST_CHECK( res.valid );
  check_exact_once( before, res.components, res.orphaned_atoms );

  // assign_atoms_to_disjoint_rois places each atom in at most one ROI.
  std::vector<RelActCalcAuto::RoiRange> rois( 2 );
  rois[0].lower_energy = 90.0; rois[0].upper_energy = 116.0;
  rois[1].lower_energy = 118.0; rois[1].upper_energy = 132.0;
  const std::vector<RoiAtom> universe = { mk_atom(100.0,1), mk_atom(112.0,1),
                                          mk_atom(124.0,1), mk_atom(200.0,1) };
  std::vector<std::vector<RoiAtom>> per_roi;
  std::vector<RoiAtom> unowned;
  assign_atoms_to_disjoint_rois( universe, rois, per_roi, unowned );
  size_t placed = unowned.size();
  for( const std::vector<RoiAtom> &v : per_roi ) placed += v.size();
  BOOST_CHECK_EQUAL( placed, universe.size() );   // exactly-once, none duplicated
  BOOST_CHECK_EQUAL( unowned.size(), 1u );         // the 200 keV atom is outside both ROIs
}

// Test G2: adaptive-style bound changes can leave an out-of-order, transitively overlapping list.
// The whole-list reconciliation must restore order, honor nonzero child-width expansion, preserve
// protected geometry, and retain every modeled atom exactly once.
BOOST_AUTO_TEST_CASE( test_reconcile_out_of_order_transitive_overlap )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  const double fwhm = 2.0;
  const double sigma = fwhm / 2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 240, 0.0f, 1.0f, density,
      { {50.0, sigma, 3000.0}, {100.0, sigma, 4000.0}, {140.0, sigma, 3500.0},
        {160.0, sigma, 4500.0}, {130.0, sigma, 1200.0}, {150.0, sigma, 1200.0} } );
  const auto snip = make_synthetic_spectrum( 240, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  AutomaticRoiComponent protected_component = mk_component(
      40.0, 60.0, fg, { mk_atom(50.0,3000.0) }, true );
  protected_component.continuum_type = PeakContinuum::OffsetType::Quadratic;
  protected_component.range_limits_type
      = RelActCalcAuto::RoiRange::RangeLimitsType::CanExpandForFwhm;
  const AutomaticRoiComponent left = mk_component(
      92.0, 150.0, fg, { mk_atom(100.0,4000.0) } );
  const AutomaticRoiComponent high = mk_component(
      110.0, 180.0, fg, { mk_atom(160.0,4500.0) } );
  const AutomaticRoiComponent middle = mk_component(
      120.0, 170.0, fg, { mk_atom(140.0,3500.0) } );
  const std::shared_ptr<const PeakDef> unmodeled_left
      = std::make_shared<const PeakDef>( 130.0, sigma, 1200.0 );
  const std::shared_ptr<const PeakDef> unmodeled_right
      = std::make_shared<const PeakDef>( 150.0, sigma, 1200.0 );

  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0;
  cons.highest_energy = 240.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 6.0;
  cons.peak_core_num_fwhm = 1.0;

  // Deliberately place the high-energy component first.  After the first pair is split, its right
  // child starts above the still-unprocessed middle component's lower bound; the former left fold
  // then merged those two in reverse anchor order and enlarged them back across the first child.
  std::vector<AutomaticRoiComponent> components
      = { high, protected_component, left, middle };
  const std::vector<AutomaticRoiComponent> before = components;
  const AutomaticRoiReconcileResult res = reconcile_automatic_components(
      components, fg, &global, fwhm_at, { unmodeled_left, unmodeled_right },
      default_policy_settings(), cons, nullptr );

  BOOST_REQUIRE_MESSAGE( res.valid, res.failure_reason );
  BOOST_CHECK( res.orphaned_atoms.empty() );
  check_exact_once( before, res.components, res.orphaned_atoms );
  BOOST_CHECK_EQUAL( atom_ids(before).size(), 4u );
  BOOST_CHECK_EQUAL( atom_ids(res.components).size(), 4u );

  for( size_t index = 1; index < res.components.size(); ++index )
  {
    BOOST_CHECK_GE( res.components[index].lower, res.components[index - 1].upper );
    BOOST_CHECK_GT( res.components[index].first_channel,
                    res.components[index - 1].last_channel );
  }
  for( const AutomaticRoiComponent &component : res.components )
  {
    if( !component.protected_geometry )
      BOOST_CHECK_GE( component.upper - component.lower, cons.min_width_fwhm*fwhm - 1.0e-6 );
    for( const RoiAtom &atom : component.atoms )
    {
      BOOST_CHECK_LE( component.lower, atom.energy - cons.peak_core_num_fwhm*fwhm + 1.0e-6 );
      BOOST_CHECK_GE( component.upper, atom.energy + cons.peak_core_num_fwhm*fwhm - 1.0e-6 );
    }
  }

  const auto protected_after = std::find_if(
      std::begin(res.components), std::end(res.components),
      []( const AutomaticRoiComponent &component ){ return component.protected_geometry; } );
  BOOST_REQUIRE( protected_after != std::end(res.components) );
  BOOST_CHECK_EQUAL( protected_after->lower, protected_component.lower );
  BOOST_CHECK_EQUAL( protected_after->upper, protected_component.upper );
  BOOST_CHECK( protected_after->continuum_type == protected_component.continuum_type );
  BOOST_CHECK( protected_after->range_limits_type == protected_component.range_limits_type );
}


// Test G3: a wide ROI whose only atom is nearer a narrow overlapping ROI's midpoint is fully
// "starved" by exact-once assignment, leaving it atom-empty.  The reconciler must keep the atom
// (in the narrow ROI) and never drop it via the zero-atom rejection.  (Regression: the zero-atom
// branch previously continue-dropped the non-empty side.)
BOOST_AUTO_TEST_CASE( test_reconcile_starved_wide_roi_keeps_atom )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  const double fwhm = 2.0;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  const auto snip = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  // Wide A=[100,140] (midpoint 120) and narrow B=[105,115] (midpoint 110) overlap; the single atom
  // at 110 is nearer B's midpoint, so exact-once assignment starves A of atoms.
  const RoiAtom a110 = mk_atom( 110.0, 3000.0 );
  AutomaticRoiComponent wide = mk_component( 100.0, 140.0, fg, {} );          // starved (empty)
  AutomaticRoiComponent narrow = mk_component( 105.0, 115.0, fg, { a110 } );  // owns the atom
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 300.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  std::vector<AutomaticRoiComponent> components = { wide, narrow };
  const std::vector<AutomaticRoiComponent> before = components;
  const AutomaticRoiReconcileResult res = reconcile_automatic_components( components, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons, nullptr );
  BOOST_CHECK( res.valid );
  check_exact_once( before, res.components, res.orphaned_atoms );
  // The atom survives somewhere.
  const bool present = std::any_of( std::begin(res.components), std::end(res.components),
      [&]( const AutomaticRoiComponent &c ){
        return std::any_of( std::begin(c.atoms), std::end(c.atoms),
            [&]( const RoiAtom &a ){ return a.id == a110.id; } ); } );
  BOOST_CHECK( present );
}


// Test H: evidence-only components (found-seed / floating features, no modeled gammas) survive
// reconciliation; the zero-atom rejection fires only for genuinely atom-empty ROIs.
BOOST_AUTO_TEST_CASE( test_reconcile_evidence_only_components )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  const double fwhm = 2.0, sigma = fwhm/2.35482;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 300, 0.0f, 1.0f, density,
      { {100.0, sigma, 5000.0}, {103.0, sigma, 4000.0} } );
  const auto snip = make_synthetic_spectrum( 300, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  AutomaticRoiComponent modeled = mk_component( 92.0, 108.0, fg, { mk_atom(100.0,5000.0) } );
  AutomaticRoiComponent evidence = mk_component( 100.0, 112.0, fg,
      { mk_atom(103.0, 4000.0, RoiAtomKind::FoundPeakEvidence) } );
  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 300.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  std::vector<AutomaticRoiComponent> components = { modeled, evidence };
  const std::vector<AutomaticRoiComponent> before = components;
  const AutomaticRoiReconcileResult res = reconcile_automatic_components( components, fg, &global,
      fwhm_at, {}, default_policy_settings(), cons, nullptr );
  BOOST_CHECK( res.valid );
  check_exact_once( before, res.components, res.orphaned_atoms );  // evidence atom survives
}


// Test I: randomized exact-once preservation.  Random small atom sets over random (possibly
// overlapping) component bounds and random protected flags must always validate, with orphans
// only ever arising from protected-boundary conflicts.
BOOST_AUTO_TEST_CASE( test_reconcile_randomized_exact_once )
{
  using FitPeaksForNuclides::detail::reconcile_automatic_components;
  const double fwhm = 2.0;
  const auto fwhm_at = [fwhm]( const double ){ return fwhm; };
  const auto density = []( const double ){ return 10.0; };
  const auto fg = make_synthetic_spectrum( 400, 0.0f, 1.0f, density );
  const auto snip = make_synthetic_spectrum( 400, 0.0f, 1.0f, density );
  GlobalContinuumEstimate global = mk_global( fg, snip );

  AutomaticRoiPartitionConstraints cons;
  cons.lowest_energy = 0.0; cons.highest_energy = 400.0;
  cons.left_barrier = -std::numeric_limits<double>::infinity();
  cons.min_width_fwhm = 0.0; cons.peak_core_num_fwhm = 1.0;

  const double core_hw = 1.0*fwhm;  // matches cons.peak_core_num_fwhm below
  std::mt19937 rng( 20260720u );
  std::uniform_int_distribution<int> ncomp_dist( 1, 5 );
  std::uniform_int_distribution<int> natom_dist( 1, 4 );
  std::uniform_real_distribution<double> center_dist( 40.0, 360.0 );
  std::uniform_real_distribution<double> half_dist( 6.0, 20.0 );
  std::uniform_real_distribution<double> prob( 0.0, 1.0 );

  for( int iter = 0; iter < 300; ++iter )
  {
    std::vector<AutomaticRoiComponent> components;
    const int ncomp = ncomp_dist( rng );
    bool used_protected = false;  // real protected ROIs are distinct; at most one per iteration
    for( int ci = 0; ci < ncomp; ++ci )
    {
      const double center = center_dist( rng );
      const double half = half_dist( rng );
      std::vector<RoiAtom> atoms;
      const int natom = natom_dist( rng );
      double min_e = center, max_e = center;
      for( int ai = 0; ai < natom; ++ai )
      {
        // Keep atoms inside a margin so their cores stay within the (core-extended) bounds below.
        const double e = std::min( 390.0, std::max( 10.0,
            center + (prob(rng) - 0.5)*2.0*(half - core_hw - 1.0) ) );
        atoms.push_back( mk_atom( e, 1000.0 ) );
        min_e = std::min( min_e, e );
        max_e = std::max( max_e, e );
      }
      // Real ROIs always include the core extent of their outermost atoms; model that here.
      const double lower = std::max( 2.0, std::min( center - half, min_e - core_hw - 0.5 ) );
      const double upper = std::min( 398.0, std::max( center + half, max_e + core_hw + 0.5 ) );
      const bool is_protected = !used_protected && (prob(rng) < 0.20);
      used_protected = used_protected || is_protected;
      components.push_back( mk_component( lower, upper, fg, atoms, is_protected ) );
    }
    const std::vector<AutomaticRoiComponent> before = components;
    const AutomaticRoiReconcileResult res = reconcile_automatic_components( components, fg, &global,
        fwhm_at, {}, default_policy_settings(), cons, nullptr );
    BOOST_REQUIRE_MESSAGE( res.valid, "randomized reconcile invalid at iter " << iter
        << ": " << res.failure_reason );
    check_exact_once( before, res.components, res.orphaned_atoms );
  }
}


BOOST_AUTO_TEST_CASE( test_select_continuum_order_by_sidebands )
{
  using FitPeaksForNuclides::detail::select_continuum_order_by_sidebands;

  // Linear continuum: sidebands are a straight line - Linear must win
  const auto lin = make_synthetic_spectrum( 3000, 0.0f, 1.0f,
      []( double e ){ return 20.0 - 0.005*e; } );
  BOOST_CHECK( select_continuum_order_by_sidebands( lin, 580.0, 620.0, 595.0, 605.0, 2.0 )
               == PeakContinuum::OffsetType::Linear );

  // Strongly curved continuum (quadratic in energy): Quadratic must win
  const auto quad = make_synthetic_spectrum( 3000, 0.0f, 1.0f,
      []( double e ){ const double d = (e - 600.0); return 50.0 + 0.05*d*d; } );
  BOOST_CHECK( select_continuum_order_by_sidebands( quad, 560.0, 640.0, 595.0, 605.0, 2.0 )
               == PeakContinuum::OffsetType::Quadratic );

  // Too few sideband channels: falls back to Linear
  BOOST_CHECK( select_continuum_order_by_sidebands( quad, 598.0, 602.0, 599.0, 601.0, 2.0 )
               == PeakContinuum::OffsetType::Linear );
}//test_select_continuum_order_by_sidebands


BOOST_AUTO_TEST_CASE( test_estimate_continuum_snip )
{
  // Flat 10 counts/keV continuum + a strong Gaussian (FWHM 10 keV), 1 keV channels.
  const double fwhm_kev = 10.0;
  const double sigma = fwhm_kev / 2.35482;
  const auto spec = make_synthetic_spectrum( 1024, 0.0f, 1.0f,
      []( double ){ return 10.0; },
      { { 512.0, sigma, 5000.0 } } );

  const std::function<double(double)> fwhm_at = [fwhm_kev]( double ){ return fwhm_kev; };

  // Order-2 filter, window 1.5*FWHM: taps land at +/-3.5 sigma, so the peak is erased cleanly.
  const auto cont = estimateContinuum( spec, fwhm_at, 1.5, 2, false, false );
  BOOST_REQUIRE( cont );
  BOOST_REQUIRE_EQUAL( cont->num_gamma_channels(), spec->num_gamma_channels() );

  // Min-filter property: never above the data
  for( size_t i = 0; i < 1024; ++i )
    BOOST_CHECK_LE( cont->gamma_channel_content(i), spec->gamma_channel_content(i) + 1.0e-3f );

  // Flat regions untouched; peak fully erased (continuum under the peak center ~ true 10)
  BOOST_CHECK_CLOSE( static_cast<double>(cont->gamma_channel_content(100)), 10.0, 1.0 );
  BOOST_CHECK( std::fabs( cont->gamma_channel_content(512) - 10.0 ) < 2.0 );

  // Order 6 at the same small window leaves an under-peak residual (the E4/E6 fractional taps
  // sample inside the peak) - it should sit clearly ABOVE the true continuum at the center.
  const auto cont_o6 = estimateContinuum( spec, fwhm_at, 1.5, 6, false, false );
  BOOST_REQUIRE( cont_o6 );
  BOOST_CHECK( cont_o6->gamma_channel_content(512) > cont->gamma_channel_content(512) + 5.0 );

  // A fwhm function returning a constant 125 keV (= 125 channels here), order 6, reproduces the
  // legacy fixed-window wrapper.
  const std::function<double(double)> legacy_win = []( double ){ return 125.0; };
  const auto cont_a = estimateContinuum( spec, legacy_win, 1.0, 6, false, false );
  const auto cont_b = estimateContinuum( spec );
  BOOST_REQUIRE( cont_a && cont_b );
  for( size_t i = 0; i < 1024; ++i )
    BOOST_CHECK_SMALL( cont_a->gamma_channel_content(i) - cont_b->gamma_channel_content(i), 1.0e-3f );

  // LLS-space and presmoothed variants (order 2): finite, and still recover the flat continuum
  // away from (and under) the peak on this noiseless spectrum
  const auto cont_lls = estimateContinuum( spec, fwhm_at, 1.5, 2, false, true );
  const auto cont_sm = estimateContinuum( spec, fwhm_at, 1.5, 2, true, false );
  BOOST_REQUIRE( cont_lls && cont_sm );
  for( const size_t i : { size_t(100), size_t(512), size_t(900) } )
  {
    BOOST_CHECK( std::isfinite( cont_lls->gamma_channel_content(i) ) );
    BOOST_CHECK( std::isfinite( cont_sm->gamma_channel_content(i) ) );
    BOOST_CHECK( std::fabs( cont_lls->gamma_channel_content(i) - 10.0 ) < 2.5 );
    BOOST_CHECK( std::fabs( cont_sm->gamma_channel_content(i) - 10.0 ) < 2.5 );
  }

  // Energy restriction: channels outside [200, 800] keV are left equal to the data (so a caller
  // gating on data-minus-continuum sees zero excess there), while the in-range continuum still
  // erases the 512 keV peak.
  const auto cont_r = estimateContinuum( spec, fwhm_at, 1.5, 2, false, false, 200.0, 800.0 );
  BOOST_REQUIRE( cont_r );
  BOOST_CHECK_EQUAL( cont_r->gamma_channel_content(50), spec->gamma_channel_content(50) );   // <200
  BOOST_CHECK_EQUAL( cont_r->gamma_channel_content(900), spec->gamma_channel_content(900) );  // >800
  BOOST_CHECK( cont_r->gamma_channel_content(512) < spec->gamma_channel_content(512) );        // in range, peak clipped
  BOOST_CHECK( std::fabs( cont_r->gamma_channel_content(512) - 10.0 ) < 2.0 );

  // Invalid inputs throw
  BOOST_CHECK_THROW( estimateContinuum( nullptr, fwhm_at, 1.5, 2, false, false ), std::exception );
  BOOST_CHECK_THROW( estimateContinuum( spec, fwhm_at, -1.0, 2, false, false ), std::exception );
  BOOST_CHECK_THROW( estimateContinuum( spec, fwhm_at, 1.5, 3, false, false ), std::exception );
  BOOST_CHECK_THROW( estimateContinuum( spec, std::function<double(double)>(), 1.5, 2, false, false ),
                     std::exception );
}//test_estimate_continuum_snip


BOOST_AUTO_TEST_CASE( test_interferer_detection_unit )
{
  using FitPeaksForNuclides::detail::InterfererCandidate;
  using FitPeaksForNuclides::detail::RequestedSourceGammas;
  using FitPeaksForNuclides::detail::find_strong_unmodeled_interferers;

  set_data_dir();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  const SandiaDecay::Nuclide * const k40   = db->nuclide( "K40" );
  const SandiaDecay::Nuclide * const eu152 = db->nuclide( "Eu152" );
  BOOST_REQUIRE( k40 && eu152 );

  const auto fwhm_at = []( double ){ return 2.0; };  // constant 2 keV FWHM (HPGe-ish near 1460)
  const double min_e = 50.0, max_e = 3000.0;

  // Multi-line Eu-152 with a weak 1457.6 keV line sitting near K-40's strong 1460.8 keV NORM line.
  RequestedSourceGammas eu;
  eu.source   = eu152;
  eu.energies = { 121.78, 344.28, 1457.64 };
  eu.yields   = {   0.28,   0.27,   0.005  };

  // Build a synthetic confirming auto-search peak (mean, area, area-uncert).
  const auto make_peak = []( double mean, double area, double uncert )
      -> std::shared_ptr<const PeakDef> {
    auto p = std::make_shared<PeakDef>( mean, 0.85, area );  // sigma ~ FWHM 2.0 keV
    p->setAmplitudeUncert( uncert );
    return p;
  };

  // Case A: K-40 1460.8 confirmed at z=50 -> exactly one K-40 nuclide candidate.
  {
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 5000.0, 100.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu }, peaks, fwhm_at, /*fit_norm_peaks=*/false, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_REQUIRE_EQUAL( c.size(), 1u );
    BOOST_CHECK_EQUAL( c[0].nuclide, k40 );
    BOOST_CHECK( !c[0].from_background_search );
    BOOST_CHECK_CLOSE( c[0].energy, 1460.82, 0.1 );
    BOOST_CHECK_GT( c[0].detection_z, 5.0 );
}

  // Case B: K-40 is itself a requested source -> already modeled, no candidate.
  {
    RequestedSourceGammas k40src;
    k40src.source   = k40;
    k40src.energies = { 1460.82 };
    k40src.yields   = { 0.1 };
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 5000.0, 100.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu, k40src }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_CHECK( c.empty() );
  }

  // Case C: 1460.8 present but not data-confirmed (z = 2 < 5) -> no candidate.
  {
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 40.0, 20.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_CHECK( c.empty() );
  }

  // Case C2: no confirming peak near the line at all -> no candidate.
  {
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 121.78, 5000.0, 100.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_CHECK( c.empty() );
  }

  // Case D: the source's own chain has a line on 1460.8 (source owns it) -> no candidate.
  {
    RequestedSourceGammas eu_owns = eu;
    eu_owns.energies.push_back( 1460.8 );
    eu_owns.yields.push_back( 0.004 );
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 5000.0, 100.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu_owns }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_CHECK( c.empty() );
  }

  // Case E: fitting NORM peaks -> K-40 already on the NORM curve, no candidate.
  {
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 5000.0, 100.0 ) };
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { eu }, peaks, fwhm_at, /*fit_norm_peaks=*/true, min_e, max_e, nullptr, nullptr, nullptr );
    BOOST_CHECK( c.empty() );
  }

  // Case F: doublet guard - a single-line source whose only line is < 1 FWHM from single-line K-40
  // is an unresolvable blend: skip and warn (only source.energies matters for the single-line test).
  {
    RequestedSourceGammas single;
    single.source   = eu152;
    single.energies = { 1460.3 };   // 0.5 keV from K-40 1460.82, well within 1 FWHM (2.0 keV)
    single.yields   = { 1.0 };
    const std::vector<std::shared_ptr<const PeakDef>> peaks = { make_peak( 1460.8, 5000.0, 100.0 ) };
    std::vector<std::string> warns;
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { single }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr, &warns );
    BOOST_CHECK( c.empty() );
    BOOST_CHECK( !warns.empty() );
  }

  // Case F2: geometric overlap alone is not a strong-interferer warning.  Without a confirming
  // foreground peak, the candidate and the structured R2 guard list must both stay empty.
  {
    RequestedSourceGammas single;
    single.source   = eu152;
    single.energies = { 1460.3 };
    single.yields   = { 1.0 };
    const std::vector<std::shared_ptr<const PeakDef>> peaks
      = { make_peak( 121.78, 5000.0, 100.0 ) };
    std::vector<std::string> warns;
    std::vector<double> guard_energies;
    const std::vector<InterfererCandidate> c = find_strong_unmodeled_interferers(
      { single }, peaks, fwhm_at, false, min_e, max_e, nullptr, nullptr, nullptr,
      &warns, nullptr, &guard_energies );
    BOOST_CHECK( c.empty() );
    BOOST_CHECK( warns.empty() );
    BOOST_CHECK( guard_energies.empty() );
  }

  // (The ambient-line sweep - Cs137/Co60 - is currently DISABLED in the helper because co-fitting an
  // ambient interferer destabilized the {K40,Eu152} joint fit; its unit test was removed with it.)

  // Regression: the {energy,parent} table refactor must not change is_near_strong_norm_gamma.
  BOOST_CHECK(  FitPeaksForNuclides::is_near_strong_norm_gamma( 1460.8, 1.0 ) );
  BOOST_CHECK( !FitPeaksForNuclides::is_near_strong_norm_gamma( 1457.6, 1.0 ) );
  BOOST_CHECK(  FitPeaksForNuclides::is_near_strong_norm_gamma( 609.31, 0.5 ) );
  BOOST_CHECK( !FitPeaksForNuclides::is_near_strong_norm_gamma( 500.0,  1.0 ) );
}//test_interferer_detection_unit


BOOST_AUTO_TEST_CASE( test_marginal_keep_classification_preserves_normal_accept_set )
{
  using FitPeaksForNuclides::detail::is_marginal_keep_reject;

  const double keep_z = 5.0;
  BOOST_CHECK( is_marginal_keep_reject( 15.1, 3.5, keep_z ) );
  BOOST_CHECK( is_marginal_keep_reject( 100.0, 5.0, keep_z ) );
  BOOST_CHECK( !is_marginal_keep_reject( 15.0, 4.0, keep_z ) );
  BOOST_CHECK( !is_marginal_keep_reject( 100.0, 3.499, keep_z ) );
  BOOST_CHECK( !is_marginal_keep_reject( 100.0, 5.001, keep_z ) );

  // The production accept predicate remains the pre-existing strict counts-and-z gate.  No point
  // can be both normally accepted and classified as marginal.
  for( const double z : { 0.0, 3.49, 3.5, 4.9, 5.0, 5.01, 8.0 } )
  {
    const bool normally_accepted = (100.0 > 15.0) && (z > keep_z);
    BOOST_CHECK( !(normally_accepted && is_marginal_keep_reject(100.0, z, keep_z)) );
  }
}

BOOST_AUTO_TEST_SUITE_END() // StatisticalDetailHelpers
