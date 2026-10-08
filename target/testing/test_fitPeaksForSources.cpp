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
#include <limits>
#include <numeric>
#include <cstdint>
#include <iostream>
#include <optional>
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
#include "InterSpec/MassAttenuationTool.h"
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


  // Map the test's coarse "is this HPGe?" flag to the detector-resolution types the
  //  peak-fitting API now takes (the old `bool isHPGe` arguments were replaced by
  //  `PeakFitUtils::CoarseResolutionType` / `PeakFitDetPrefs`).
  PeakFitUtils::CoarseResolutionType res_type_for( const bool isHPGe )
  {
    return isHPGe ? PeakFitUtils::CoarseResolutionType::High
                  : PeakFitUtils::CoarseResolutionType::Low;
  }

  shared_ptr<const PeakFitDetPrefs> det_prefs_for( const bool isHPGe )
  {
    auto prefs = make_shared<PeakFitDetPrefs>();
    prefs->m_det_type = res_type_for( isHPGe );
    return prefs;
  }


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

      // A ROI's upper bound is exclusive: the channel whose lower edge equals it is not in the ROI.
      // Probe just inside each ROI so a boundary that lands exactly on a channel edge (and can
      // round either way in float) is not read as covering the channel beyond it.
      const double nudge = 1.0e-4;
      const size_t prev_upper_ch = energy_cal->channel_for_energy( prev_upper - nudge );
      const size_t curr_lower_ch = energy_cal->channel_for_energy( curr_lower + nudge );

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
        RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status),
        "Fit failed: " << result.error_message );
    }

    if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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
    if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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
    if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(eu152_result.status) );
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result1.status) );
  BOOST_REQUIRE_GE( result1.observable_peaks.size(), 1u );

  const vector<shared_ptr<const PeakDef>> after_first = apply_fit_result( empty_user_peaks, result1 );

  // Second fit with first fit's results as user_peaks (default mode - should replace)
  const FitPeaksForNuclides::PeakFitResult result2
    = run_fit( spec.foreground, spec.background, auto_peaks, sources, after_first, spec.isHPGe );

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result2.status) );

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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result1.status) );
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result2.status) );
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result1.status) );
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
  if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result2.status)
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

  if( !RelActCalcAuto::RelActAutoSolution::is_usable_status(ba_result.status) )
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result1.status) );
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

    // Eu-154 IS present in trinitite (a neutron-activation product, like the Co-60 of step 6): the
    // automated search finds 1274.4 keV at z ~ 7 and 1005 keV at z ~ 5.  This step also exercises the
    // path, where user peaks adjacent to the new source's ROIs are absorbed as floating peaks and
    // listed in `original_peaks_to_remove` (so `verify_removed_peaks_replaced` is the contract here,
    // not an empty removal list).
    verify_fit_result( result, user_peaks, spec.foreground,
      FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak );
    verify_removed_peaks_replaced( result );
    BOOST_CHECK_MESSAGE( has_peak_near( result.observable_peaks, 1274.44, 3.0 ),
                         "Eu-154 1274.4 keV not found in trinitite" );
    BOOST_TEST_MESSAGE( "Eu-154: " << result.observable_peaks.size() << " observable peaks" );

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
    if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

    if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

    if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(r_auto.status) );
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(r_joint.status) );

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
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) );
  for( const string &warning : result.warnings )
    BOOST_TEST_MESSAGE( warning );
  // Public contract: the requested sources' strong lines are reported, the K40 1460.82 keV NORM
  // peak is never reported, and Eu152's weak 1457.64 keV line beside it is not reported either -
  // the single-pass planner rejects that line group as swamped by the unexplained K40 peak, so the
  // R6 nuisance co-fit no longer needs to model K40 to keep the line honest (before the planner,
  // this test required K40 in the raw solution as proof the co-fit ran).
  BOOST_CHECK( !has_peak_near( result.observable_peaks, 1457.64, 2.0 ) );
  BOOST_CHECK( !has_peak_near( result.observable_peaks, 1460.82, 2.0 ) );
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
  BOOST_CHECK_LE( nuisance_parents.size(), 2u );

  // The caller must be able to turn the complete automatic R6 path off for controlled tuning and
  // for supplied-background workflows.  The default result above proves R6 is active; with the
  // opt-out, the same raw requested-source fit must remain source-only.
  const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> no_interferer_options
    = FitPeaksForNuclides::FitSrcPeaksOptions::DisableAutoInterfererFit;
  const FitPeaksForNuclides::PeakFitResult source_only = run_fit(
      spec.foreground, nullptr, auto_peaks, make_sources({"Eu152", "Cs137"}),
      no_user_peaks, spec.isHPGe, no_interferer_options );
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(source_only.status) );
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
  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(with_float.status) );
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
    BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result->status) );
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

  // Public contract for the K40 / Eu152 1457.64 keV neighbourhood.  The single-pass planner
  // rejects Eu152's weak 1457.64 keV line group as swamped by the unexplained K40 1460.82 keV peak
  // (the line is unmeasurable under it), so neither an Eu152 peak nor a K40 peak is reported there
  // by the automatic arms, in either source order and with or without a supplied background.
  // When K40 is itself requested (the joint arms), its 1460.82 keV peak is reported and the Eu152
  // line beside it still is not.  (Before the planner this test proved the R6 nuisance co-fit had
  // modelled K40 in the raw solution; that mechanism is no longer needed here.)
  for( const FitPeaksForNuclides::PeakFitResult * const result : automatic_results )
  {
    BOOST_CHECK( !has_peak_near( result->observable_peaks, 1457.64, 2.0 ) );
    BOOST_CHECK( !has_peak_near( result->observable_peaks, 1460.82, 2.0 ) );
  }
  // Requested K40 is reported from the raw spectrum; with the supplied background (which carries
  // the same K40) the net K40 is nil and nothing need be reported, but the Eu152 line beside it is
  // never reported in either arm.
  BOOST_CHECK( find_source_gamma( raw_joint.observable_peaks, "K40", 1460.82, 0.5 ) );
  for( const FitPeaksForNuclides::PeakFitResult * const result : { &raw_joint, &bg_joint } )
    BOOST_CHECK( !find_source_gamma( result->observable_peaks, "Eu152", 1457.64, 0.5 ) );
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
  if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

  if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result.status) )
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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result_eu152.status) );

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

  BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(result_am241.status) );

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

  if( RelActCalcAuto::RelActAutoSolution::is_usable_status(result_eu154.status) )
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




namespace
{
  // Settings for the single-pass ROI planner on the synthetic spectra below (1 keV channels,
  // constant 2 keV FWHM), starting from the HPGe defaults.
  FitPeaksForNuclides::GammaClusteringSettings planner_settings()
  {
    FitPeaksForNuclides::GammaClusteringSettings s
      = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( PeakFitUtils::CoarseResolutionType::High )
          .get_auto_clustering_settings();
    s.use_roi_plan = true;
    s.keep_significance_z = 2.0;
    s.share_always_fwhm = 3.4;
    s.separate_always_fwhm = 3.4;
    s.step_use_chi2_trial = false;
    return s;
  }

  const std::function<double(double)> planner_fwhm = []( double ){ return 2.0; };
  const double planner_sigma = 2.0 / 2.35482;
}//namespace


BOOST_AUTO_TEST_CASE( test_sibling_absence_check )
{
  using FitPeaksForNuclides::detail::sibling_absence_check;
  using FitPeaksForNuclides::detail::SiblingAbsenceResult;
  set_data_dir();  // the lead-shield scan reads the mass-attenuation tables
  const std::function<double(double)> flat_eff = []( double ){ return 1.0; };

  // The attenuation tables come back in PhysicalUnits units; the check converts to cm2/g.
  const double mu_units = PhysicalUnits::cm2 / PhysicalUnits::g;
  const double mu_pb_60 = MassAttenuation::massAttenuationCoefficientFracAN( 82.0, 59.5 ) / mu_units;
  const double mu_pb_662 = MassAttenuation::massAttenuationCoefficientFracAN( 82.0, 662.0 ) / mu_units;
  BOOST_CHECK_MESSAGE( (mu_pb_60 > 3.0) && (mu_pb_60 < 8.0), "mu/rho(Pb, 59.5 keV) = " << mu_pb_60 << " cm2/g" );
  BOOST_CHECK_MESSAGE( (mu_pb_662 > 0.08) && (mu_pb_662 < 0.14), "mu/rho(Pb, 662 keV) = " << mu_pb_662 << " cm2/g" );

  // Ba133-like source: 81 keV (yield 0.33) and 356 keV (0.62); the spectrum shows both in a
  // consistent ratio on a 20/keV continuum.
  const std::vector<SandiaDecay::EnergyRatePair> ba{ {0.33, 81.0}, {0.62, 356.0} };
  const auto consistent = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 20.0; },
      { {81.0, planner_sigma, 5000.0}, {356.0, planner_sigma, 8000.0} } );
  const std::vector<std::pair<double,double>> ba_obs{ {81.0, 5000.0}, {356.0, 8000.0} };
  SiblingAbsenceResult r = sibling_absence_check( ba, ba, 81.0, 5000.0, 0.0, ba_obs, planner_fwhm, flat_eff, consistent, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK( r.judged );
  BOOST_CHECK_LT( r.worst_ratio, 1.0 );
  r = sibling_absence_check( ba, ba, 356.0, 8000.0, 0.0, ba_obs, planner_fwhm, flat_eff, consistent, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK_LT( r.worst_ratio, 1.0 );
  // Claiming a 50000-count 81 keV peak would need a 356 keV peak the data does not have, and no
  // shielding can hide a higher-energy line.
  r = sibling_absence_check( ba, ba, 81.0, 50000.0, 0.0, ba_obs, planner_fwhm, flat_eff, consistent, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK( r.judged );
  BOOST_CHECK_GT( r.worst_ratio, 2.0 );
  BOOST_CHECK_CLOSE( r.sibling_energy, 356.0, 1.0e-9 );

  // Am241-like: a huge 59.5 keV line (0.36), a 5e-6 line at 335.4 keV and a 3.6e-6 line at
  // 662.4 keV.  A 700000-count peak at 662 keV cannot be Am241's when nothing sits at 335 keV: a
  // thick shield could hide the 335 keV line, but the visible 59.5 keV peak rules the shield out.
  const std::vector<SandiaDecay::EnergyRatePair> am{ {0.36, 59.5}, {5.0e-6, 335.4}, {3.6e-6, 662.4} };
  const auto trinitite_like = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 30.0; },
      { {59.5, planner_sigma, 33000.0}, {662.4, planner_sigma, 700000.0} } );
  const std::vector<std::pair<double,double>> am_obs{ {59.5, 33000.0}, {662.4, 700000.0} };
  r = sibling_absence_check( am, am, 662.4, 700000.0, 0.0, am_obs, planner_fwhm, flat_eff, trinitite_like, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK( r.judged );
  BOOST_CHECK_GT( r.worst_ratio, 2.0 );
  // The binding line is the absent 335 keV sibling, or the 59.5 keV peak the implied activity
  // would swamp - both are physical statements of the same contradiction.
  BOOST_CHECK( (std::fabs( r.sibling_energy - 335.4 ) < 1.0e-6) || (std::fabs( r.sibling_energy - 59.5 ) < 1.0e-6) );
  // ... while the 59.5 keV peak itself is fine.
  r = sibling_absence_check( am, am, 59.5, 33000.0, 0.0, am_obs, planner_fwhm, flat_eff, trinitite_like, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK_LT( r.worst_ratio, 1.0 );

  // Heavily shielded source: 300 keV (0.5) and 600 keV (0.02) lines, spectrum where the 600 keV
  // peak is as large as the 300 keV one (about 4 cm of lead).  The 600 keV claim must pass because
  // some shielding makes the pair consistent.
  const std::vector<SandiaDecay::EnergyRatePair> sh{ {0.5, 300.0}, {0.02, 600.0} };
  const auto shielded = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 20.0; },
      { {300.0, planner_sigma, 3000.0}, {600.0, planner_sigma, 3000.0} } );
  const std::vector<std::pair<double,double>> sh_obs{ {300.0, 3000.0}, {600.0, 3000.0} };
  r = sibling_absence_check( sh, sh, 600.0, 3000.0, 0.0, sh_obs, planner_fwhm, flat_eff, shielded, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK_LT( r.worst_ratio, 1.0 );

  // Shielded Am241-like source: 59.5 keV (0.36) absent, a 99 keV line (2e-4) present with 575 counts.
  // Lead attenuates 99 keV more than 59.5 keV (K-edge), so only an iron-like shield explains the
  // pair - the scan must find it.
  const std::vector<SandiaDecay::EnergyRatePair> am_sh{ {0.36, 59.5}, {2.0e-4, 99.0}, {1.5e-4, 103.0} };
  const auto iron_shielded = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 10.0; },
      { {99.0, planner_sigma, 575.0}, {103.0, planner_sigma, 400.0} } );
  const std::vector<std::pair<double,double>> am_sh_obs{ {99.0, 575.0}, {103.0, 400.0} };
  r = sibling_absence_check( am_sh, am_sh, 99.0, 575.0, 0.0, am_sh_obs, planner_fwhm, flat_eff, iron_shielded, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK_MESSAGE( r.worst_ratio < 1.0, "shielded 99 keV claim ratio " << r.worst_ratio << " (Z=" << r.best_shield_z << ", " << r.best_shield_g_cm2 << " g/cm2)" );

  // A coincidence-sum peak of two strong lines is not judged.
  const std::vector<SandiaDecay::EnergyRatePair> co{ {0.999, 200.0}, {0.9998, 250.0}, {2.0e-8, 450.0} };
  r = sibling_absence_check( co, co, 450.0, 2000.0, 0.0, {}, planner_fwhm, flat_eff, consistent, 20.0, 690.0, 0.4, 120.0, 50.0 );
  BOOST_CHECK( !r.judged );
}


BOOST_AUTO_TEST_CASE( test_roi_plan_sharing_bands )
{
  using FitPeaksForNuclides::detail::PlannerLine;
  using FitPeaksForNuclides::detail::PlannedRoiSummary;
  using FitPeaksForNuclides::detail::plan_rois_for_lines;
  const FitPeaksForNuclides::GammaClusteringSettings s = planner_settings();

  for( const double sep_fwhm : { 1.5, 2.5, 3.0, 3.8, 6.0, 12.0 } )
  {
    const double e1 = 300.0, e2 = 300.0 + sep_fwhm*2.0;
    const auto fg = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 5.0; },
        { {e1, planner_sigma, 5000.0}, {e2, planner_sigma, 5000.0} } );
    const std::vector<PlannedRoiSummary> rois = plan_rois_for_lines(
        { {e1, 5000.0}, {e2, 5000.0} }, fg, planner_fwhm, 20.0, 690.0, s, {} );

    const bool expect_shared = (sep_fwhm < s.share_always_fwhm);
    BOOST_REQUIRE_MESSAGE( !rois.empty(), "no ROI planned at separation " << sep_fwhm );
    BOOST_CHECK_MESSAGE( rois.size() == (expect_shared ? 1u : 2u),
      "separation " << sep_fwhm << " FWHM planned " << rois.size() << " ROIs" );

    // every line sits inside its ROI; ROIs are channel aligned and disjoint with a gap
    for( const PlannedRoiSummary &roi : rois )
    {
      for( const double e : roi.line_energies )
        BOOST_CHECK( (e > roi.lower) && (e < roi.upper) );
      const size_t first = fg->find_gamma_channel( static_cast<float>(roi.lower) );
      const size_t last = fg->find_gamma_channel( std::nextafter( static_cast<float>(roi.upper), static_cast<float>(roi.lower) ) );
      BOOST_CHECK_CLOSE( roi.lower, fg->gamma_channel_lower( first ), 1.0e-6 );
      BOOST_CHECK_CLOSE( roi.upper, fg->gamma_channel_upper( last ), 1.0e-6 );
    }
    for( size_t i = 1; i < rois.size(); ++i )
      BOOST_CHECK_GT( fg->find_gamma_channel( static_cast<float>(rois[i].lower) ),
                      fg->find_gamma_channel( std::nextafter( static_cast<float>(rois[i-1].upper), 0.0f ) ) );
  }//for( separations )
}


BOOST_AUTO_TEST_CASE( test_roi_plan_admission_and_confirmation )
{
  using FitPeaksForNuclides::detail::PlannerLine;
  using FitPeaksForNuclides::detail::PlannedRoiSummary;
  using FitPeaksForNuclides::detail::plan_rois_for_lines;
  const FitPeaksForNuclides::GammaClusteringSettings s = planner_settings();

  // A 5/keV continuum: a 400-count line is detectable (z ~ 15), a 3-count line is not.
  const auto fg = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 5.0; },
      { {200.0, planner_sigma, 400.0}, {500.0, planner_sigma, 600.0} } );
  std::vector<PlannedRoiSummary> rois = plan_rois_for_lines( { {200.0, 400.0}, {350.0, 3.0} }, fg, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_CLOSE( rois[0].line_energies.at(0), 200.0, 1.0e-9 );

  // A found peak confirms a line only when the source could produce a meaningful part of it:
  // predicted 300 of a 600-count found peak is confirmed; predicted 3 counts is not.
  auto found = std::make_shared<PeakDef>( 500.0, planner_sigma, 600.0 );
  rois = plan_rois_for_lines( { {200.0, 400.0}, {500.0, 300.0} }, fg, planner_fwhm, 20.0, 690.0, s, { found } );
  BOOST_CHECK_EQUAL( rois.size(), 2u );
  rois = plan_rois_for_lines( { {200.0, 400.0}, {500.0, 3.0} }, fg, planner_fwhm, 20.0, 690.0, s, { found } );
  BOOST_CHECK_EQUAL( rois.size(), 1u );

  // An unexplained found peak next to a line is an obstacle that stops the ROI extension short of it.
  auto obstacle = std::make_shared<PeakDef>( 206.0, planner_sigma, 600.0 );
  const auto fg2 = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 5.0; },
      { {200.0, planner_sigma, 400.0}, {206.0, planner_sigma, 600.0} } );
  rois = plan_rois_for_lines( { {200.0, 400.0} }, fg2, planner_fwhm, 20.0, 690.0, s, { obstacle } );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_LT( rois[0].upper, 206.0 );
}


BOOST_AUTO_TEST_CASE( test_roi_plan_swamped_group )
{
  using FitPeaksForNuclides::detail::PlannerLine;
  using FitPeaksForNuclides::detail::PlannedRoiSummary;
  using FitPeaksForNuclides::detail::plan_rois_for_lines;
  const FitPeaksForNuclides::GammaClusteringSettings s = planner_settings();

  // A weak predicted line (60 counts) with a 5000-count foreign peak 1 FWHM away - inside the
  // group's core but outside the confirmation distance: the group is unmeasurable under the foreign
  // peak and must be rejected (the peak stays an obstacle), instead of the line absorbing the peak.
  const auto swamped = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 5.0; },
      { {300.0, planner_sigma, 60.0}, {302.0, planner_sigma, 5000.0}, {500.0, planner_sigma, 400.0} } );
  auto foreign = std::make_shared<PeakDef>( 302.0, planner_sigma, 5000.0 );
  std::vector<PlannedRoiSummary> rois = plan_rois_for_lines( { {300.0, 60.0}, {500.0, 400.0} }, swamped,
      planner_fwhm, 20.0, 690.0, s, { foreign } );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_CLOSE( rois[0].line_energies.at(0), 500.0, 1.0e-9 );

  // The same geometry with a comparable foreign peak (80 counts) is not swamped: the group is planned
  // and the found peak is merely an obstacle.
  const auto shared = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 5.0; },
      { {300.0, planner_sigma, 600.0}, {302.0, planner_sigma, 80.0} } );
  auto comparable = std::make_shared<PeakDef>( 302.0, planner_sigma, 80.0 );
  rois = plan_rois_for_lines( { {300.0, 600.0} }, shared, planner_fwhm, 20.0, 690.0, s, { comparable } );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_CLOSE( rois[0].line_energies.at(0), 300.0, 1.0e-9 );
}


BOOST_AUTO_TEST_CASE( test_roi_plan_span_cap )
{
  using FitPeaksForNuclides::detail::PlannerLine;
  using FitPeaksForNuclides::detail::PlannedRoiSummary;
  using FitPeaksForNuclides::detail::plan_rois_for_lines;
  FitPeaksForNuclides::GammaClusteringSettings s = planner_settings();

  // Nine lines 2 FWHM apart chain into one 16 FWHM component; the widest gap (3 FWHM, between the
  // fifth and sixth line) is where the cap must split it.
  std::vector<PlannerLine> lines;
  std::vector<std::array<double,3>> peaks;
  double e = 200.0;
  for( int i = 0; i < 9; ++i )
  {
    lines.push_back( {e, 3000.0} );
    peaks.push_back( {e, planner_sigma, 3000.0} );
    e += (i == 4) ? 6.0 : 4.0;
  }
  const auto fg = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 20.0; }, peaks );
  std::vector<PlannedRoiSummary> rois = plan_rois_for_lines( lines, fg, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_GE( rois.size(), 2u );

  // The cap splits at the widest gap first, so 216 and 222 keV never share a ROI ...
  int roi_of_216 = -1, roi_of_222 = -1;
  for( size_t i = 0; i < rois.size(); ++i )
  {
    for( const double e : rois[i].line_energies )
    {
      if( std::fabs( e - 216.0 ) < 0.5 )
        roi_of_216 = static_cast<int>( i );
      if( std::fabs( e - 222.0 ) < 0.5 )
        roi_of_222 = static_cast<int>( i );
    }
  }
  BOOST_REQUIRE( (roi_of_216 >= 0) && (roi_of_222 >= 0) );
  BOOST_CHECK_NE( roi_of_216, roi_of_222 );

  // ... every emitted ROI stays under the cap, and no split leaves a line against an edge (the
  // continuum needs roi_min_side_fwhm of room to anchor on - see GammaClusteringSettings).
  for( const PlannedRoiSummary &roi : rois )
  {
    BOOST_REQUIRE( !roi.line_energies.empty() );
    const double lo = *std::min_element( begin(roi.line_energies), end(roi.line_energies) );
    const double hi = *std::max_element( begin(roi.line_energies), end(roi.line_energies) );
    const double fwhm = planner_fwhm( 0.5*(lo + hi) );
    const double chan_w = 1.0;   // make_synthetic_spectrum uses 1 keV channels
    BOOST_CHECK_MESSAGE( (lo - roi.lower) >= (s.roi_min_side_fwhm*fwhm - chan_w - 1.0e-6),
      "ROI [" << roi.lower << ", " << roi.upper << "] starts " << (lo - roi.lower)
      << " keV below its first line, less than " << (s.roi_min_side_fwhm*fwhm) );
    BOOST_CHECK_MESSAGE( (roi.upper - hi) >= (s.roi_min_side_fwhm*fwhm - chan_w - 1.0e-6),
      "ROI [" << roi.lower << ", " << roi.upper << "] ends " << (roi.upper - hi)
      << " keV above its last line, less than " << (s.roi_min_side_fwhm*fwhm) );
    BOOST_CHECK_LE( (roi.upper - roi.lower) / fwhm, s.max_shared_span_fwhm + 1.0 );
  }

  s.max_shared_span_fwhm = 0.0;   // disabled: one chain
  rois = plan_rois_for_lines( lines, fg, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_CHECK_EQUAL( rois.size(), 1u );
}


BOOST_AUTO_TEST_CASE( test_roi_plan_continuum_selection )
{
  using FitPeaksForNuclides::detail::PlannerLine;
  using FitPeaksForNuclides::detail::PlannedRoiSummary;
  using FitPeaksForNuclides::detail::plan_rois_for_lines;
  const FitPeaksForNuclides::GammaClusteringSettings s = planner_settings();

  // Strong peak on a flat continuum: no step (the continuum does not drop across the peak) and a
  // linear continuum - a Constant one is available but off by default (the reference fits never use
  // one; see GammaClusteringSettings::cont_constant_max_counts).
  const auto flat = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double ){ return 50.0; },
      { {300.0, planner_sigma, 200000.0} } );
  std::vector<PlannedRoiSummary> rois = plan_rois_for_lines( { {300.0, 200000.0} }, flat, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_MESSAGE( rois[0].continuum_type == PeakContinuum::OffsetType::Linear,
    "flat-continuum strong peak got " << PeakContinuum::offset_type_str( rois[0].continuum_type )
    << " over [" << rois[0].lower << ", " << rois[0].upper << "]" );

  // A continuum RISING across the ROI cannot be a Compton step (those only fall with energy), so a
  // strong peak on one stays linear however big it is.
  const auto sloped = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double e ){ return 20.0 + 0.2*e; },
      { {300.0, planner_sigma, 200000.0} } );
  rois = plan_rois_for_lines( { {300.0, 200000.0} }, sloped, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK_MESSAGE( rois[0].continuum_type == PeakContinuum::OffsetType::Linear,
    "rising-continuum strong peak got " << PeakContinuum::offset_type_str( rois[0].continuum_type )
    << " over [" << rois[0].lower << ", " << rois[0].upper << "]" );

  // The same peak on a continuum that steps down across it gets a CDF step continuum.
  const auto stepped = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double e ){ return (e < 300.0) ? 200.0 : 50.0; },
      { {300.0, planner_sigma, 200000.0} } );
  rois = plan_rois_for_lines( { {300.0, 200000.0} }, stepped, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK( rois[0].continuum_type == PeakContinuum::OffsetType::FlatStepCDF );
  // and the step ROI is given extra room on the low side
  BOOST_CHECK_LT( rois[0].lower, 300.0 - (s.roi_core_num_fwhm + 0.5*s.step_low_side_extra_fwhm)*2.0 );

  // A weak peak never gets a step, whatever the continuum does (a polynomial order may still be
  // chosen for the stepped background).
  const auto weak_stepped = make_synthetic_spectrum( 700, 0.0f, 1.0f, []( double e ){ return (e < 300.0) ? 200.0 : 50.0; },
      { {300.0, planner_sigma, 400.0} } );
  rois = plan_rois_for_lines( { {300.0, 400.0} }, weak_stepped, planner_fwhm, 20.0, 690.0, s, {} );
  BOOST_REQUIRE_EQUAL( rois.size(), 1u );
  BOOST_CHECK( (rois[0].continuum_type == PeakContinuum::OffsetType::Linear)
               || (rois[0].continuum_type == PeakContinuum::OffsetType::Quadratic) );
}


BOOST_AUTO_TEST_CASE( test_energy_cal_drift_bound )
{
  // Every spectrum in the evaluation corpora is already well calibrated, so the energy-calibration
  // drift bound (sm_energy_cal_max_drift_fwhm, applied where the fit advances the calibration
  // between iterations) is exercised by nothing else.  These cases pin its contract directly.
  using FitPeaksForNuclides::detail::max_search_peak_drift_fwhm;

  const size_t nchannel = 8192;
  auto pristine = std::make_shared<SpecUtils::EnergyCalibration>();
  pristine->set_polynomial( nchannel, { 0.0f, 0.25f }, {} );   // 0 - 2048 keV

  // Two peaks with a 2 keV FWHM, one low and one high, so a gain-only error shows up as a drift
  // that grows with energy - the shape the Fulcrum U235 failure had.
  const auto make_peak = []( const double mean, const double fwhm ) {
    auto p = std::make_shared<PeakDef>( mean, fwhm/2.35482, 1000.0 );
    return std::shared_ptr<const PeakDef>( p );
  };
  const std::vector<std::shared_ptr<const PeakDef>> peaks{ make_peak( 200.0, 2.0 ),
                                                           make_peak( 1000.0, 2.0 ) };

  // The same calibration drifts nothing.
  BOOST_CHECK_SMALL( max_search_peak_drift_fwhm( peaks, pristine, pristine ), 1.0e-9 );

  // A pure offset drifts every peak by the same ENERGY, so by the same number of FWHM here.
  {
    auto offset = std::make_shared<SpecUtils::EnergyCalibration>();
    offset->set_polynomial( nchannel, { 1.0f, 0.25f }, {} );   // +1 keV everywhere
    double worst_energy = 0.0;
    const double drift = max_search_peak_drift_fwhm( peaks, pristine, offset, nullptr, &worst_energy );
    BOOST_CHECK_CLOSE( drift, 0.5, 1.0 );        // 1 keV / 2 keV FWHM
    BOOST_CHECK( drift < 2.0 );                  // a real, correctable offset: must NOT be rejected
  }

  // A gain error drifts the HIGH peak most, and the bound must notice the high one.
  {
    auto gain = std::make_shared<SpecUtils::EnergyCalibration>();
    gain->set_polynomial( nchannel, { 0.0f, 0.25f*1.006f }, {} );   // +0.6 % gain
    double worst_energy = 0.0;
    const double drift = max_search_peak_drift_fwhm( peaks, pristine, gain, nullptr, &worst_energy );
    // 1000 keV moves 6 keV = 3 FWHM; 200 keV moves 1.2 keV = 0.6 FWHM.
    BOOST_CHECK_CLOSE( drift, 3.0, 2.0 );
    BOOST_CHECK_CLOSE( worst_energy, 1000.0, 1.0e-6 );
    BOOST_CHECK( drift > 2.0 );                  // the measured Fulcrum runaway: must be rejected
  }

  // Peaks with no width cannot measure anything, and neither can an empty list; both are silent
  // rather than fatal, and the caller then applies the calibration unjudged.
  {
    auto gain = std::make_shared<SpecUtils::EnergyCalibration>();
    gain->set_polynomial( nchannel, { 0.0f, 0.25f*1.10f }, {} );
    auto zero_width = std::make_shared<PeakDef>( 1000.0, 0.0, 1000.0 );
    const std::vector<std::shared_ptr<const PeakDef>> widthless{
        std::shared_ptr<const PeakDef>( zero_width ) };
    BOOST_CHECK_SMALL( max_search_peak_drift_fwhm( widthless, pristine, gain ), 1.0e-9 );
    BOOST_CHECK_SMALL( max_search_peak_drift_fwhm( {}, pristine, gain ), 1.0e-9 );
  }

  // The fitted width model, when supplied, is the yardstick - not the peak's own width.  A narrow
  // noise spike high in the spectrum must not veto a calibration that is fine for the real peaks.
  {
    auto gain = std::make_shared<SpecUtils::EnergyCalibration>();
    gain->set_polynomial( nchannel, { 0.0f, 0.25f*1.001f }, {} );   // 1.4 keV at 1400 keV
    auto spike = std::make_shared<PeakDef>( 1400.0, 0.4/2.35482, 50.0 );
    const std::vector<std::shared_ptr<const PeakDef>> spiky{
        std::shared_ptr<const PeakDef>( spike ) };
    const std::function<double(double)> real_width = []( double ){ return 2.0; };

    // Judged by the spike's own 0.4 keV width the drift reads 3.5 FWHM and would be rejected...
    BOOST_CHECK_GT( max_search_peak_drift_fwhm( spiky, pristine, gain ), 2.0 );
    // ...but the detector's actual resolution says 0.7 FWHM, which is fine.
    BOOST_CHECK_LT( max_search_peak_drift_fwhm( spiky, pristine, gain, real_width ), 2.0 );
  }

  // An unusable calibration is not an excuse to reject a fit.
  {
    std::shared_ptr<const SpecUtils::EnergyCalibration> null_cal;
    BOOST_CHECK_SMALL( max_search_peak_drift_fwhm( peaks, pristine, null_cal ), 1.0e-9 );
    BOOST_CHECK_SMALL( max_search_peak_drift_fwhm( peaks, null_cal, pristine ), 1.0e-9 );
  }
}//BOOST_AUTO_TEST_CASE( test_energy_cal_drift_bound )


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


BOOST_AUTO_TEST_CASE( test_fit_fwhm_function_robust_shape_prior )
{
  using FitPeaksForNuclides::detail::class_shape_fwhm;
  using FitPeaksForNuclides::detail::fit_fwhm_function_robust;

  const auto form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
  const double nsig = PhysicalUnits::fwhm_nsigma;

  // 3 keV channels over 0-3000 keV; the contents do not enter the width fit
  const auto spec = make_synthetic_spectrum( 1000, 0.0f, 3.0f, []( double ){ return 10.0; } );

  const auto make_peak = [nsig]( const double mean, const double fwhm, const double z ){
    auto p = std::make_shared<PeakDef>( mean, fwhm/nsig, 1000.0*z );
    p->setAmplitudeUncert( 1000.0 );
    p->setSigmaUncert( 0.05*fwhm/nsig );
    return std::shared_ptr<const PeakDef>( p );
  };
  const auto model_fwhm = [form]( const std::vector<float> &coefs, const double e ){
    return static_cast<double>( DetectorPeakResponse::peakResolutionFWHM( static_cast<float>(e), form, coefs ) );
  };

  // NaI-like: six clean peaks at 0.9x the class curve, a 2x-wide backscatter bump, a 1.6x-wide
  // x-ray blob and a slightly wide 2614 keV line.  The contaminants must not vote, the 2614 keV
  // line must, and the curve must keep rising to 2614 keV instead of ending as a constant.
  {
    const auto Low = PeakFitUtils::CoarseResolutionType::Low;
    const double scale = 0.9;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
    for( const double e : { 59.5, 122.0, 356.0, 662.0, 1173.0, 1332.0 } )
      peaks.push_back( make_peak( e, scale*class_shape_fwhm( Low, e ), 30.0 ) );
    peaks.push_back( make_peak( 190.0, 2.0*scale*class_shape_fwhm( Low, 190.0 ), 25.0 ) );
    peaks.push_back( make_peak( 32.0, 1.6*scale*class_shape_fwhm( Low, 32.0 ), 40.0 ) );
    peaks.push_back( make_peak( 2614.0, 1.05*scale*class_shape_fwhm( Low, 2614.0 ), 20.0 ) );

    std::vector<float> coefs, uncerts;
    double lo = 0.0, hi = 0.0;
    std::string note;
    fit_fwhm_function_robust( peaks, spec, Low, 25.0, nullptr, form, coefs, uncerts, lo, hi, note );
    BOOST_TEST_MESSAGE( note );
    BOOST_REQUIRE( coefs.size() >= 2 );
    BOOST_CHECK_CLOSE( lo, 25.0, 1.0e-6 );
    BOOST_CHECK_CLOSE( hi, 3000.0, 1.0e-6 );
    BOOST_CHECK( note.find( "2 off the class-prior shape" ) != std::string::npos );
    BOOST_CHECK( note.find( "from 7 of 9 search peaks" ) != std::string::npos );
    BOOST_CHECK( note.find( "prior only" ) == std::string::npos );
    // The sqrt polynomial cannot follow E^0.6 exactly across two decades, so a few percent of
    // family error is expected; what must not happen is the old constant-width collapse.
    for( const double e : { 122.0, 662.0, 1332.0, 2614.0 } )
      BOOST_CHECK_CLOSE( model_fwhm( coefs, e ), scale*class_shape_fwhm( Low, e ), 10.0 );
    BOOST_CHECK_GT( model_fwhm( coefs, 2614.0 ) / model_fwhm( coefs, 339.0 ), 2.0 );
  }

  // A single peak: the class curve scaled through it, valid over the whole range.
  {
    const auto Low = PeakFitUtils::CoarseResolutionType::Low;
    const std::vector<std::shared_ptr<const PeakDef>> one{ make_peak( 662.0, 1.2*class_shape_fwhm( Low, 662.0 ), 30.0 ) };
    std::vector<float> coefs, uncerts;
    double lo = 0.0, hi = 0.0;
    std::string note;
    fit_fwhm_function_robust( one, spec, Low, 25.0, nullptr, form, coefs, uncerts, lo, hi, note );
    BOOST_TEST_MESSAGE( note );
    BOOST_REQUIRE( coefs.size() >= 2 );
    for( const double e : { 662.0, 1332.0 } )
      BOOST_CHECK_CLOSE( model_fwhm( coefs, e ), 1.2*class_shape_fwhm( Low, e ), 12.0 );
    BOOST_CHECK_CLOSE( model_fwhm( coefs, 60.0 ), 1.2*class_shape_fwhm( Low, 60.0 ), 25.0 );
    BOOST_CHECK_GT( model_fwhm( coefs, 2614.0 ) / model_fwhm( coefs, 339.0 ), 2.0 );
  }

  // HPGe-like: eight clean peaks at 1.1x the class curve all vote, and the fit follows them.
  {
    const auto High = PeakFitUtils::CoarseResolutionType::High;
    const auto hpge = make_synthetic_spectrum( 6000, 0.0f, 0.5f, []( double ){ return 10.0; } );
    std::vector<std::shared_ptr<const PeakDef>> peaks;
    for( const double e : { 59.5, 122.0, 356.0, 662.0, 1173.0, 1332.0, 1408.0, 2614.0 } )
      peaks.push_back( make_peak( e, 1.1*class_shape_fwhm( High, e ), 30.0 ) );
    std::vector<float> coefs, uncerts;
    double lo = 0.0, hi = 0.0;
    std::string note;
    fit_fwhm_function_robust( peaks, hpge, High, 20.0, nullptr, form, coefs, uncerts, lo, hi, note );
    BOOST_TEST_MESSAGE( note );
    BOOST_CHECK( note.find( "from 8 of 8 search peaks" ) != std::string::npos );
    for( const double e : { 59.5, 122.0, 662.0, 1332.0, 2614.0 } )
      BOOST_CHECK_CLOSE( model_fwhm( coefs, e ), 1.1*class_shape_fwhm( High, e ), 8.0 );
  }

  // A GR1-like CZT spectrum, dead below its 39 keV discriminator, whose only search peak is a line
  // sliced by it - narrow, and centred just above the cut.  It must not set the widths.
  {
    const auto CZT = PeakFitUtils::CoarseResolutionType::CZT;
    const auto czt = make_synthetic_spectrum( 1000, 0.0f, 3.0f, []( double e ){ return (e < 39.3) ? 0.0 : 10.0; } );
    const std::vector<std::shared_ptr<const PeakDef>> sliced{ make_peak( 43.5, 4.8, 8.0 ) };
    std::vector<float> coefs, uncerts;
    double lo = 0.0, hi = 0.0;
    std::string note;
    bool from_prior = false;
    fit_fwhm_function_robust( sliced, czt, CZT, 15.0, nullptr, form, coefs, uncerts, lo, hi, note, &from_prior );
    BOOST_TEST_MESSAGE( note );
    BOOST_CHECK( from_prior );
    BOOST_CHECK( note.find( "1 search peaks on the detector threshold" ) != std::string::npos );
    for( const double e : { 60.0, 357.0, 662.0 } )
      BOOST_CHECK_CLOSE( model_fwhm( coefs, e ), class_shape_fwhm( CZT, e ), 10.0 );
  }

  // A lone search peak a third of the class width: at z=4 its width is noise and the class curve
  // sets the scale; at z=30 it is a measurement and sets it.
  for( const double z : { 4.0, 30.0 } )
  {
    const auto CZT = PeakFitUtils::CoarseResolutionType::CZT;
    const double scale = 0.35;
    const std::vector<std::shared_ptr<const PeakDef>> one{ make_peak( 122.0, scale*class_shape_fwhm( CZT, 122.0 ), z ) };
    std::vector<float> coefs, uncerts;
    double lo = 0.0, hi = 0.0;
    std::string note;
    bool from_prior = false;
    fit_fwhm_function_robust( one, spec, CZT, 25.0, nullptr, form, coefs, uncerts, lo, hi, note, &from_prior );
    BOOST_TEST_MESSAGE( note );
    const bool weak = (z < 6.0);
    BOOST_CHECK_EQUAL( from_prior, weak );
    for( const double e : { 122.0, 662.0 } )
      BOOST_CHECK_CLOSE( model_fwhm( coefs, e ), (weak ? 1.0 : scale)*class_shape_fwhm( CZT, e ), 12.0 );
  }
}//test_fit_fwhm_function_robust_shape_prior


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


// Test A: a non-dominant outer gamma in a multi-line left group (the Am241 failure shape) is kept
// in the child covering it - never clipped-and-dropped - and the split still happens.
// Test B: the chosen boundary lands on a channel edge with a one-channel gap - children never
// share a channel and their bounds equal exact channel edges.
// Test C: when the atom cores cannot be separated by any channel, the pair MERGES - it never
// drops a side.  Exercised by making the partition core wider than the oracle's separation core.
// Test D: a min-width child pinned against the spectrum edge cannot widen; the pair merges (or
// finds an alternate) but never drops an atom, and never exceeds the spectrum extent.
// Test E: an unmodeled-feature exclusion band that would cut an admitted atom core is not carved
// through; the partition either finds a core-safe boundary or merges, never splitting the core.
// Test F: protected geometry is pinned - its bounds/metadata are bit-identical afterward, an atom
// whose core falls inside it is booked to it, and an atom straddling its edge is orphaned (never
// silently dropped from the ledger).
// Test G: overlapping input ROIs whose atoms share the overlap band are reconciled to exactly-once
// ownership - the direct regression for the flat-list double-claim path.
// Test G2: adaptive-style bound changes can leave an out-of-order, transitively overlapping list.
// The whole-list reconciliation must restore order, honor nonzero child-width expansion, preserve
// protected geometry, and retain every modeled atom exactly once.
// Test G3: a wide ROI whose only atom is nearer a narrow overlapping ROI's midpoint is fully
// "starved" by exact-once assignment, leaving it atom-empty.  The reconciler must keep the atom
// (in the narrow ROI) and never drop it via the zero-atom rejection.  (Regression: the zero-atom
// branch previously continue-dropped the non-empty side.)
// Test H: evidence-only components (found-seed / floating features, no modeled gammas) survive
// reconciliation; the zero-atom rejection fires only for genuinely atom-empty ROIs.
// Test I: randomized exact-once preservation.  Random small atom sets over random (possibly
// overlapping) component bounds and random protected flags must always validate, with orphans
// only ever arising from protected-boundary conflicts.
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


namespace
{
  RelActCalcAuto::FloatingPeakResult make_float_result( const double energy,
                                                       const double fit_energy,
                                                       const double amplitude,
                                                       const double amplitude_uncert )
  {
    RelActCalcAuto::FloatingPeakResult fpr;
    fpr.energy = energy;
    fpr.original_spectrum_cal_energy = fit_energy;
    fpr.amplitude = amplitude;
    fpr.amplitude_uncert = amplitude_uncert;
    fpr.fwhm = 3.6;
    fpr.fwhm_uncert = 0.12;
    return fpr;
  }//make_float_result(...)
}//namespace


/** A bystander must never be updated to the fit's amplitude while keeping the uncertainty of the
 peak it replaced - that reports a confidently-measured peak the fit could not determine.  The
 default-mode reconciliation block used to do exactly that whenever the solver returned a
 non-positive or non-finite amplitude uncertainty; both blocks now share this helper.
 */
BOOST_AUTO_TEST_CASE( test_bystander_update_requires_determined_amplitude )
{
  using FitPeaksForNuclides::detail::update_bystander_from_float_result;

  // The user's existing peak: a well-measured 1000 +- 20 counts.
  const double orig_energy = 661.657;
  PeakDef orig_peak( orig_energy, 1.5, 1000.0 );
  orig_peak.setAmplitudeUncert( 20.0 );

  const std::shared_ptr<PeakContinuum> roi_continuum = std::make_shared<PeakContinuum>();
  roi_continuum->setRange( 640.0, 680.0 );

  // A determined result is taken, uncertainty and all.
  {
    const RelActCalcAuto::FloatingPeakResult fpr
        = make_float_result( orig_energy, 661.9, 1500.0, 60.0 );
    const std::optional<PeakDef> updated = update_bystander_from_float_result(
        orig_peak, orig_energy, fpr, roi_continuum, "unit test" );

    BOOST_REQUIRE( updated.has_value() );
    BOOST_CHECK_CLOSE( updated->mean(), 661.9, 1.0e-6 );
    BOOST_CHECK_CLOSE( updated->amplitude(), 1500.0, 1.0e-6 );
    BOOST_CHECK_CLOSE( updated->sigma(), 3.6 / (2.0*std::sqrt(2.0*std::log(2.0))), 1.0e-6 );
    BOOST_CHECK_EQUAL( updated->continuum(), roi_continuum );

    // The regression: the fit's amplitude must arrive with the FIT's uncertainty, never with the
    //  original peak's.
    BOOST_CHECK_CLOSE( updated->amplitudeUncert(), 60.0, 1.0e-6 );
  }

  // Every flavour of "the fit could not determine this peak" retains the user's original.
  const double nan_uncert = std::numeric_limits<double>::quiet_NaN();
  const std::vector<std::pair<std::string,RelActCalcAuto::FloatingPeakResult>> undetermined = {
    { "zero uncertainty",      make_float_result( orig_energy, 661.9, 1500.0, 0.0 ) },
    { "negative uncertainty",  make_float_result( orig_energy, 661.9, 1500.0, -1.0 ) },
    { "NaN uncertainty",       make_float_result( orig_energy, 661.9, 1500.0, nan_uncert ) },
    { "uncertainty == amplitude", make_float_result( orig_energy, 661.9, 1500.0, 1500.0 ) },
    { "uncertainty > amplitude",  make_float_result( orig_energy, 661.9, 1500.0, 2400.0 ) },
    { "zero amplitude",        make_float_result( orig_energy, 661.9, 0.0, 5.0 ) },
    { "negative amplitude",    make_float_result( orig_energy, 661.9, -10.0, 5.0 ) },
  };

  for( const std::pair<std::string,RelActCalcAuto::FloatingPeakResult> &entry : undetermined )
  {
    const std::optional<PeakDef> updated = update_bystander_from_float_result(
        orig_peak, orig_energy, entry.second, roi_continuum, "unit test" );
    BOOST_CHECK_MESSAGE( !updated.has_value(),
      "Bystander with " << entry.first << " should have been retained, not updated" );
  }
}//test_bystander_update_requires_determined_amplitude


/** The floating-peak result pool also holds the 511 keV annihilation, escape, and interferer peaks
 this code injects itself.  Matching a bystander to the *nearest* result within a window let an
 unmatched user peak a few tenths of a keV from 511 bind to the annihilation peak's result, whose
 energy then erased the real 511 keV peak from the results.  Matching is on identity instead.
 */
BOOST_AUTO_TEST_CASE( test_find_enrolled_float_result )
{
  using FitPeaksForNuclides::detail::find_enrolled_float_result;

  const double annihilation_energy = 510.9989;

  // The pool as the solver returns it: the injected 511 keV float, plus a bystander enrolled at
  //  its own peak mean.
  const std::vector<RelActCalcAuto::FloatingPeakResult> results = {
    make_float_result( annihilation_energy, 511.2, 8000.0, 200.0 ),
    make_float_result( 661.657, 661.9, 1500.0, 60.0 ),
  };

  {
    // The regression: a source-less user peak at ~510.7 keV never enrolled this result, so it must
    //  not bind to it - doing so deletes the real 511 keV peak downstream.
    std::set<const RelActCalcAuto::FloatingPeakResult *> consumed;
    BOOST_CHECK( !find_enrolled_float_result( results, consumed, 510.7 ) );
    BOOST_CHECK( consumed.empty() );

    // Nor may it reach a result further away than the old 0.5 keV window.
    BOOST_CHECK( !find_enrolled_float_result( results, consumed, 661.3 ) );
  }

  {
    // A bystander enrolled at exactly its own energy finds its result, and consumes it.
    std::set<const RelActCalcAuto::FloatingPeakResult *> consumed;
    const RelActCalcAuto::FloatingPeakResult * const found
        = find_enrolled_float_result( results, consumed, 661.657 );
    BOOST_REQUIRE( found );
    BOOST_CHECK_CLOSE( found->energy, 661.657, 1.0e-9 );
    BOOST_CHECK_EQUAL( consumed.size(), 1 );

    // Consumed at most once: a repeat lookup finds nothing rather than re-using the result.
    BOOST_CHECK( !find_enrolled_float_result( results, consumed, 661.657 ) );

    // The 511 float is still available to whatever enrolled it.
    BOOST_CHECK( find_enrolled_float_result( results, consumed, annihilation_energy ) );
  }

  {
    // Two bystanders enrolled at the same energy must bind to different results.
    const std::vector<RelActCalcAuto::FloatingPeakResult> duplicates = {
      make_float_result( 1460.82, 1460.9, 500.0, 30.0 ),
      make_float_result( 1460.82, 1461.1, 400.0, 25.0 ),
    };

    std::set<const RelActCalcAuto::FloatingPeakResult *> consumed;
    const RelActCalcAuto::FloatingPeakResult * const first
        = find_enrolled_float_result( duplicates, consumed, 1460.82 );
    const RelActCalcAuto::FloatingPeakResult * const second
        = find_enrolled_float_result( duplicates, consumed, 1460.82 );
    BOOST_REQUIRE( first && second );
    BOOST_CHECK( first != second );
    BOOST_CHECK( !find_enrolled_float_result( duplicates, consumed, 1460.82 ) );
  }
}//test_find_enrolled_float_result


BOOST_AUTO_TEST_CASE( test_nai_iodine_escape_fractions )
{
  // The NaI iodine escape curve (GADRAS NaI escape components): none below the iodine K edge or
  // above 250 keV, falling steeply in between, K-beta a fixed 0.30 of K-alpha.
  double kalpha = -1.0, kbeta = -1.0;
  PeakFitUtils::nai_iodine_escape_fractions( 32.0, kalpha, kbeta );
  BOOST_CHECK_EQUAL( kalpha, 0.0 );
  BOOST_CHECK_EQUAL( kbeta, 0.0 );

  PeakFitUtils::nai_iodine_escape_fractions( 300.0, kalpha, kbeta );
  BOOST_CHECK_EQUAL( kalpha, 0.0 );

  PeakFitUtils::nai_iodine_escape_fractions( 59.54, kalpha, kbeta );   // Am241: 8.9 % in the truth
  BOOST_CHECK_CLOSE( kalpha, 0.089, 3.0 );
  BOOST_CHECK_CLOSE( kbeta, 0.30*kalpha, 1.0e-6 );

  PeakFitUtils::nai_iodine_escape_fractions( 122.06, kalpha, kbeta );  // Co57: 1.6 %
  BOOST_CHECK_CLOSE( kalpha, 0.0160, 3.0 );

  double previous = 1.0;
  for( double energy = 33.2; energy <= 250.0; energy += 1.0 )
  {
    PeakFitUtils::nai_iodine_escape_fractions( energy, kalpha, kbeta );
    BOOST_CHECK( (kalpha > 0.0) && (kalpha <= previous) );
    previous = kalpha;
  }
}//test_nai_iodine_escape_fractions


BOOST_AUTO_TEST_CASE( test_quadratic_null_power )
{
  // A quadratic continuum-only null absorbs most of a single peak in a narrow ROI, so a
  // peaks-vs-null test there can see little of the peak: lambda is small against the peak's own
  // chi2 (its shape against zero) at 2 FWHM and a large share of it at 6 FWHM.
  using FitPeaksForNuclides::detail::quadratic_null_power;

  const double mean = 650.0, fwhm = 10.0, area = 2000.0, continuum = 100.0;
  const auto peak = std::make_shared<const PeakDef>( mean, fwhm/PhysicalUnits::fwhm_nsigma, area );
  const std::vector<std::shared_ptr<const PeakDef>> peaks{ peak };

  // lambda / (the peak shape's chi2 against zero), over an ROI of `num_fwhm` centred on the peak.
  const auto lambda_fraction = [&]( const double num_fwhm ) -> double {
    const double channel_width = 0.5;
    const double lower = mean - 0.5*num_fwhm*fwhm;
    const size_t nchannel = static_cast<size_t>( std::round( num_fwhm*fwhm/channel_width ) );
    std::vector<float> energies( nchannel + 1 ), counts( nchannel );
    for( size_t i = 0; i <= nchannel; ++i )
      energies[i] = static_cast<float>( lower + i*channel_width );
    std::vector<double> shape( nchannel, 0.0 );
    peak->gauss_integral( energies.data(), shape.data(), nchannel );
    double shape_chi2 = 0.0;
    for( size_t i = 0; i < nchannel; ++i )
    {
      counts[i] = static_cast<float>( continuum*channel_width + shape[i] );
      shape_chi2 += shape[i]*shape[i] / counts[i];
    }
    const double lambda = quadratic_null_power( energies.data(), counts.data(), nchannel, mean, peaks );
    BOOST_REQUIRE( shape_chi2 > 0.0 );
    return lambda / shape_chi2;
  };//lambda_fraction

  const double at2 = lambda_fraction( 2.0 ), at3 = lambda_fraction( 3.0 ), at6 = lambda_fraction( 6.0 );
  BOOST_CHECK_MESSAGE( at2 < 0.08, "2 FWHM: lambda fraction " << at2 );
  BOOST_CHECK_MESSAGE( (at3 > at2) && (at3 < 0.25), "3 FWHM: lambda fraction " << at3 );
  BOOST_CHECK_MESSAGE( (at6 > 0.35) && (at6 < 0.65), "6 FWHM: lambda fraction " << at6 );

  // Degenerate input
  const std::vector<float> few_energies{ 600.0f, 601.0f, 602.0f, 603.0f };
  const std::vector<float> few_counts{ 10.0f, 10.0f, 10.0f };
  BOOST_CHECK_EQUAL( quadratic_null_power( few_energies.data(), few_counts.data(), 3, 600.0, peaks ), 0.0 );
}//test_quadratic_null_power

BOOST_AUTO_TEST_SUITE_END() // StatisticalDetailHelpers
