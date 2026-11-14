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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "InterSpec_config.h"

#include <cmath>
#include <limits>
#include <string>
#include <memory>
#include <vector>
#include <algorithm>

// See test_RelActCalcAuto_MassFracConstraint.cpp: winsock2.h must come before <windows.h> on MSVC.
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#define BOOST_TEST_MODULE RelActCalcAuto_Skew_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for the Rel. Act. tool is required for this test." );

namespace
{
  string g_test_file_dir;

  void set_data_dir()
  {
    static bool s_have_set = false;
    if( s_have_set )
      return;
    s_have_set = true;

    const int argc = framework::master_test_suite().argc;
    char ** const argv = framework::master_test_suite().argv;

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
      for( const char *d : { "data", "../data", "../../data", "../../../data" } )
      {
        if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
        {
          datadir = d;
          break;
        }
      }
    }

    if( g_test_file_dir.empty() )
    {
      for( const char *d : { "test_data", "../test_data", "../../test_data",
                             "../../target/testing/test_data" } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path(d, "RelActAutoReport") ) )
        {
          g_test_file_dir = d;
          break;
        }
      }
    }

    const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ),
                           "sandia.decay.xml not at '" << sandia_decay_file << "'" );

    BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    BOOST_REQUIRE_MESSAGE( db, "Error initializing SandiaDecayDataBase" );
  }//set_data_dir()


  /** The uranium HPGe test spectrum, its detector response, and the RelActAuto options embedded in it. */
  struct UraniumProblem
  {
    shared_ptr<const SpecUtils::Measurement> foreground, background;
    shared_ptr<const DetectorPeakResponse> drf;
    RelActCalcAuto::Options options;
  };

  UraniumProblem load_uranium_problem()
  {
    set_data_dir();

    const string spec_path = SpecUtils::append_path(
        SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), "U235_Unshielded_6000.n42" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ), "Test spectrum not at '" << spec_path << "'" );

    auto meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ),
                           "Failed to load '" << spec_path << "'" );

    UraniumProblem problem;
    for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
    {
      if( !m )
        continue;
      if( m->source_type() == SpecUtils::SourceType::Background )
        problem.background = m;
      else if( !problem.foreground )
        problem.foreground = m;
    }
    BOOST_REQUIRE_MESSAGE( problem.foreground, "No foreground measurement in '" << spec_path << "'" );

    problem.drf = meas->detector();
    BOOST_REQUIRE_MESSAGE( problem.drf, "Test N42 has no embedded detector response." );

    const unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
    BOOST_REQUIRE_MESSAGE( state, "Test N42 missing embedded RelActAuto state." );
    BOOST_REQUIRE_MESSAGE( !state->options.rel_eff_curves.empty(), "Embedded state has no rel-eff curve." );
    problem.options = state->options;

    return problem;
  }//load_uranium_problem()


  /** Solves `problem` with every skew coefficient of `type` fixed - `fixed_lower` at the spectrum's first
   channel, `fixed_upper` (energy-dependent coefficients only; NaN for none) at its last - and checks every
   peak gets the straight line through them, and that the reported values are the fixed ones. */
  void check_fixed_skew_respected( UraniumProblem problem, const PeakDef::SkewType type,
                                   const vector<double> &fixed_lower, const vector<double> &fixed_upper )
  {
    RelActCalcAuto::Options &options = problem.options;
    options.set_skew_type( type );
    options.lorentzian_xrays = false;   // needs NoSkew/GaussPlusBortel
    options.auto_simplify_model = false;
    options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;   // peak means = fit means

    const size_t num_skew = PeakDef::num_skew_parameters( type );
    BOOST_REQUIRE( (fixed_lower.size() == num_skew) && (fixed_upper.size() == num_skew) );
    for( size_t i = 0; i < num_skew; ++i )
    {
      options.fixed_lower_skew[i] = fixed_lower[i];
      if( !std::isnan( fixed_upper[i] ) )
        options.fixed_upper_skew[i] = fixed_upper[i];
    }

    RelActCalcAuto::RelActAutoSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, problem.foreground, problem.background,
                                    problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                                    nullptr ) );
    BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                           "Solve failed: status=" << static_cast<int>(sol.m_status)
                           << ", error='" << sol.m_error_message << "'" );
    BOOST_REQUIRE( !sol.m_fit_peaks_in_spectrums_cal.empty() );

    const size_t nchannel = problem.foreground->num_gamma_channels();
    const double lower_energy = problem.foreground->gamma_channel_lower( 0 );
    const double upper_energy = problem.foreground->gamma_channel_upper( nchannel - 1 );

    for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )
    {
      BOOST_REQUIRE( peak.skewType() == type );

      const double frac = (peak.mean() - lower_energy) / (upper_energy - lower_energy);
      for( size_t i = 0; i < num_skew; ++i )
      {
        const double expected = std::isnan( fixed_upper[i] ) ? fixed_lower[i]
                                         : (fixed_lower[i] + frac*(fixed_upper[i] - fixed_lower[i]));
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
        BOOST_CHECK_MESSAGE( fabs( peak.coefficient(ct) - expected ) < 1.0E-6*std::max( 1.0, fabs(expected) ),
                             "Peak at " << peak.mean() << " keV: skew coefficient " << i << " is "
                             << peak.coefficient(ct) << ", expected the fixed value " << expected );
      }
    }//for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )

    // The reported skew values are the fixed ones, and not fit.
    const vector<RelActCalcAuto::RelActAutoSolution::SkewParResult> skew = sol.skew_parameters();
    BOOST_REQUIRE_EQUAL( skew.size(), num_skew );
    for( size_t i = 0; i < num_skew; ++i )
    {
      BOOST_CHECK( !skew[i].was_fit );
      BOOST_CHECK( !skew[i].clamped );
      BOOST_CHECK_CLOSE( skew[i].lower, fixed_lower[i], 1.0E-6 );
      BOOST_CHECK_CLOSE( skew[i].upper, (std::isnan(fixed_upper[i]) ? fixed_lower[i] : fixed_upper[i]), 1.0E-6 );
    }
  }//check_fixed_skew_respected(...)


  /** Checks that every fitted peak's skew lies on one straight line in energy per coefficient (a constant,
   for coefficients without energy dependence), and that the solution reports that line at the spectrum's
   first-channel lower and last-channel upper energies - limited to the coefficient's allowed range, and
   flagged `clamped`, exactly when the line leaves it there.  Returns whether any coefficient was clamped. */
  bool check_reported_skew_is_the_fitted_line( const RelActCalcAuto::RelActAutoSolution &sol,
                                               const SpecUtils::Measurement &foreground )
  {
    const PeakDef::SkewType type = sol.m_options.skew_type;
    const size_t num_skew = PeakDef::num_skew_parameters( type );
    const vector<RelActCalcAuto::RelActAutoSolution::SkewParResult> skew = sol.skew_parameters();
    BOOST_REQUIRE_EQUAL( skew.size(), num_skew );

    const vector<PeakDef> &peaks = sol.m_fit_peaks_in_spectrums_cal;
    BOOST_REQUIRE_GE( peaks.size(), 2u );
    const auto by_energy = []( const PeakDef &a, const PeakDef &b ){ return a.mean() < b.mean(); };
    const PeakDef &lowest = *std::min_element( begin(peaks), end(peaks), by_energy );
    const PeakDef &highest = *std::max_element( begin(peaks), end(peaks), by_energy );
    BOOST_REQUIRE_GT( highest.mean() - lowest.mean(), 10.0 );

    const size_t nchannel = foreground.num_gamma_channels();
    const double spectrum_lower = foreground.gamma_channel_lower( 0 );
    const double spectrum_upper = foreground.gamma_channel_upper( nchannel - 1 );

    bool any_clamped = false;
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
      const double c_low = lowest.coefficient( ct ), c_high = highest.coefficient( ct );
      const auto line = [&]( const double energy ) -> double {
        return c_low + (energy - lowest.mean()) * (c_high - c_low) / (highest.mean() - lowest.mean());
      };

      for( const PeakDef &peak : peaks )
        BOOST_CHECK_MESSAGE( fabs( peak.coefficient(ct) - line(peak.mean()) ) < 1.0E-7*std::max( 1.0, fabs(c_low) ),
                             "Peak at " << peak.mean() << " keV: skew coefficient " << i << " is "
                             << peak.coefficient(ct) << ", off the line through the other peaks ("
                             << line(peak.mean()) << ")" );

      double range_lower = 0.0, range_upper = 0.0, start = 0.0, step = 0.0;
      BOOST_REQUIRE( PeakDef::skew_parameter_range( type, ct, range_lower, range_upper, start, step ) );
      const double at_lower = line( spectrum_lower ), at_upper = line( spectrum_upper );
      const bool leaves_range = (at_lower < range_lower) || (at_lower > range_upper)
                                || (at_upper < range_lower) || (at_upper > range_upper);
      BOOST_CHECK_MESSAGE( skew[i].clamped == leaves_range, "Skew coefficient " << i << ": the fitted line is "
                           << at_lower << " / " << at_upper << " at the spectrum ends (range " << range_lower
                           << " to " << range_upper << "), but clamped=" << skew[i].clamped );
      const double tol = 1.0E-6*std::max( 1.0, fabs(at_lower) + fabs(at_upper) );
      BOOST_CHECK_MESSAGE( fabs( skew[i].lower - std::clamp(at_lower, range_lower, range_upper) ) < tol,
                           "Skew coefficient " << i << ": reported low-energy value " << skew[i].lower
                           << ", the fitted line gives " << at_lower );
      BOOST_CHECK_MESSAGE( fabs( skew[i].upper - std::clamp(at_upper, range_lower, range_upper) ) < tol,
                           "Skew coefficient " << i << ": reported high-energy value " << skew[i].upper
                           << ", the fitted line gives " << at_upper );
      any_clamped = any_clamped || skew[i].clamped;
    }//for( loop over skew coefficients )

    const bool warned = std::any_of( begin(sol.m_warnings), end(sol.m_warnings), []( const string &w ){
      return SpecUtils::icontains( w, "would not reproduce this fit" );
    } );
    BOOST_CHECK_MESSAGE( warned == any_clamped, "Clamped skew warned about: " << warned
                         << ", any coefficient clamped: " << any_clamped );

    return any_clamped;
  }//check_reported_skew_is_the_fitted_line(...)


  /** `options` with every ROI made `LineAnchored` from its current limits, so it is re-sized from the fitted
   peak shape and the problem re-fit (with the first fit's model selection replayed). */
  RelActCalcAuto::Options with_line_anchored_rois( RelActCalcAuto::Options options )
  {
    for( RelActCalcAuto::RoiRange &roi : options.rois )
      roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::LineAnchored;
    return options;
  }
}//namespace


// Skew values the caller fixes are the skew every peak gets: an energy-dependent coefficient is the
//  straight line through its fixed values at the spectrum's first-channel lower and last-channel upper
//  energies (the convention PeakFitDetPrefs, PeakFitLM and the saved RelActAuto options share), and
//  any other coefficient is its fixed value.  Double Sided Crystal Ball is the type that used to be
//  special-cased (fit Gaussian first from hard-coded values, which replaced the fixed ones).
BOOST_AUTO_TEST_CASE( fixed_skew_is_respected )
{
  const UraniumProblem problem = load_uranium_problem();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // {alpha_low, n_low, alpha_high, n_high}; the alphas are energy dependent, the n's are not.
  check_fixed_skew_respected( problem, PeakDef::SkewType::DoubleSidedCrystalBall,
                              { 1.2, 3.0, 2.8, 6.0 }, { 2.4, nan, 1.6, nan } );

  // Both ExpGaussExp tails are energy dependent; the second is fixed without one.
  check_fixed_skew_respected( problem, PeakDef::SkewType::ExpGaussExp, { 0.9, 2.5 }, { 2.0, nan } );
}//BOOST_AUTO_TEST_CASE( fixed_skew_is_respected )


// The reported skew values (at the spectrum's ends) are the fitted skew's line there, limited to the allowed
//  range (and flagged and warned about) only when that line leaves it; unclamped, fixed at them, a re-solve
//  gives every peak the same skew.
BOOST_AUTO_TEST_CASE( reported_skew_reproduces_the_fit )
{
  // Fits `problem` with `type`, checks the reported values against the fitted line, and - if they were
  //  not clamped - that fixing the skew at them reproduces it; returns whether they were clamped.
  const auto check = []( const UraniumProblem &problem, const PeakDef::SkewType type ) -> bool {
    RelActCalcAuto::Options options = problem.options;
    options.set_skew_type( type );
    options.lorentzian_xrays = false;
    options.auto_simplify_model = false;
    options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;

    const auto solve = [&problem]( const RelActCalcAuto::Options &opts ) -> RelActCalcAuto::RelActAutoSolution {
      RelActCalcAuto::RelActAutoSolution sol;
      BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( opts, problem.foreground, problem.background,
                                        problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                                        nullptr ) );
      BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                             "Solve failed: status=" << static_cast<int>(sol.m_status)
                             << ", error='" << sol.m_error_message << "'" );
      return sol;
    };

    const RelActCalcAuto::RelActAutoSolution fit = solve( options );
    const bool any_clamped = check_reported_skew_is_the_fitted_line( fit, *problem.foreground );
    BOOST_TEST_MESSAGE( PeakDef::to_string(type) << ": reported skew clamped: " << any_clamped );

    // Clamped values can not reproduce the fit (the check above has confirmed that is all they are).
    if( any_clamped )
      return true;

    const vector<RelActCalcAuto::RelActAutoSolution::SkewParResult> skew = fit.skew_parameters();
    RelActCalcAuto::Options fixed = options;
    for( size_t i = 0; i < skew.size(); ++i )
    {
      fixed.fixed_lower_skew[i] = skew[i].lower;
      fixed.fixed_upper_skew[i] = skew[i].upper;
    }

    // Fixing the skew can change the ROIs (they are sized with it), so peaks are matched by energy (with
    //  the energy calibration not fit, a peak's mean is its line's energy).
    const RelActCalcAuto::RelActAutoSolution refit = solve( fixed );
    size_t num_matched = 0;
    for( const PeakDef &original : fit.m_fit_peaks_in_spectrums_cal )
    {
      const auto pos = std::find_if( begin(refit.m_fit_peaks_in_spectrums_cal), end(refit.m_fit_peaks_in_spectrums_cal),
                                     [&original]( const PeakDef &p ){ return fabs(p.mean() - original.mean()) < 1.0E-6; } );
      if( pos == end(refit.m_fit_peaks_in_spectrums_cal) )
        continue;
      num_matched += 1;
      const PeakDef &again = *pos;
      for( size_t i = 0; i < skew.size(); ++i )
      {
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
        BOOST_CHECK_MESSAGE( fabs( original.coefficient(ct) - again.coefficient(ct) ) < 1.0E-6,
                             "Peak at " << original.mean() << " keV: skew coefficient " << i << " fit "
                             << original.coefficient(ct) << ", but fixed at the reported values gives "
                             << again.coefficient(ct) );
      }
    }//for( const PeakDef &original : fit.m_fit_peaks_in_spectrums_cal )

    BOOST_CHECK_GT( num_matched, fit.m_fit_peaks_in_spectrums_cal.size() / 2 );
    return false;
  };//check lambda

  const UraniumProblem problem = load_uranium_problem();

  // The full spectrum: the HPGe's high-energy tail is at its no-skew limit, and its line, extrapolated
  //  to 0 and 3000 keV, leaves the allowed range.
  check( problem, PeakDef::SkewType::ExpGaussExp );

  // Trimmed to just past the fit range, with only the low-energy tail, the line stays in range - so this
  //  is the case that exercises fixing the skew at the reported values.
  UraniumProblem trimmed = problem;
  {
    auto foreground = make_shared<SpecUtils::Measurement>( *problem.foreground );
    double lowest_roi = std::numeric_limits<double>::max(), highest_roi = 0.0;
    for( const RelActCalcAuto::RoiRange &roi : problem.options.rois )
    {
      lowest_roi = std::min( lowest_roi, roi.lower_energy );
      highest_roi = std::max( highest_roi, roi.upper_energy );
    }
    const size_t first = foreground->find_gamma_channel( static_cast<float>(lowest_roi - 10.0) );
    const size_t last = foreground->find_gamma_channel( static_cast<float>(highest_roi + 10.0) );
    BOOST_REQUIRE( first < last );
    foreground->truncate_gamma_channels( first, last, false );
    trimmed.foreground = foreground;
    trimmed.background = nullptr;
  }
  BOOST_CHECK_MESSAGE( !check( trimmed, PeakDef::SkewType::GaussExp ),
                       "The trimmed-spectrum case, meant to exercise the re-solve, was clamped" );
}//BOOST_AUTO_TEST_CASE( reported_skew_reproduces_the_fit )


// Auto-simplify removes peak skew in nested steps - an energy dependence first, then the component (a tail)
//  together with the coefficients that then no longer matter (e.g., a Crystal Ball power law).  With every
//  removal forced through, the skew ends at its no-skew values at every energy, and replays (ROI refits)
//  do not fail on the -999.9 placeholder that a removed energy dependence holds.
BOOST_AUTO_TEST_CASE( nested_skew_simplification )
{
  const UraniumProblem problem = load_uranium_problem();

  for( const PeakDef::SkewType type : { PeakDef::SkewType::ExpGaussExp, PeakDef::SkewType::DoubleSidedCrystalBall } )
  {
    RelActCalcAuto::Options options = problem.options;
    options.set_skew_type( type );
    options.lorentzian_xrays = false;
    options.auto_simplify_model = true;
    options.auto_simplify_max_dchi2 = 1.0e9;  // accept every converged removal

    RelActCalcAuto::RelActAutoSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, problem.foreground, problem.background,
                                      problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                                      nullptr ) );
    BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                           PeakDef::to_string(type) << " solve failed: status=" << static_cast<int>(sol.m_status)
                           << ", error='" << sol.m_error_message << "'" );

    for( const string &warning : sol.m_warnings )
      BOOST_CHECK_MESSAGE( !SpecUtils::icontains( warning, "could not replay" ), "Unexpected warning: " << warning );

    const size_t num_skew = PeakDef::num_skew_parameters( type );
    vector<double> no_skew( num_skew, 0.0 );
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
      BOOST_REQUIRE( PeakDef::skew_no_skew_value( type, ct, no_skew[i] ) );
    }

    for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )
    {
      for( size_t i = 0; i < num_skew; ++i )
      {
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
        BOOST_CHECK_MESSAGE( fabs( peak.coefficient(ct) - no_skew[i] ) < 1.0E-9,
                             PeakDef::to_string(type) << " peak at " << peak.mean() << " keV: skew coefficient "
                             << i << " is " << peak.coefficient(ct) << ", not its no-skew value " << no_skew[i] );
      }
    }//for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )

    const vector<RelActCalcAuto::RelActAutoSolution::SkewParResult> skew = sol.skew_parameters();
    BOOST_REQUIRE_EQUAL( skew.size(), num_skew );
    for( size_t i = 0; i < num_skew; ++i )
    {
      BOOST_CHECK_CLOSE( skew[i].lower, no_skew[i], 1.0E-9 );
      BOOST_CHECK_CLOSE( skew[i].upper, no_skew[i], 1.0E-9 );
      BOOST_CHECK( skew[i].simplified );
      BOOST_CHECK( !skew[i].clamped );
    }

    BOOST_CHECK_NE( sol.m_frozen_model_policy_hash, UINT64_C(0) );
  }//for( both skew types )
}//BOOST_AUTO_TEST_CASE( nested_skew_simplification )


// ROIs re-sized from the fitted peak shape move the fit edges the skew is anchored at, so the refit must
//  re-express the first fit's skew line at its own edges, and replay the first fit's skew simplification
//  (the -999.9 placeholder and the dropped tail) on the re-sized problem.
BOOST_AUTO_TEST_CASE( skew_survives_roi_refits )
{
  const UraniumProblem problem = load_uranium_problem();

  const auto solve = [&problem]( const RelActCalcAuto::Options &opts ) -> RelActCalcAuto::RelActAutoSolution {
    RelActCalcAuto::RelActAutoSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( opts, problem.foreground, problem.background,
                                      problem.drf, {}, PeakFitUtils::coarse_det_type( problem.foreground, nullptr ),
                                      nullptr ) );
    BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status( sol.m_status ),
                           "Solve failed: status=" << static_cast<int>(sol.m_status)
                           << ", error='" << sol.m_error_message << "'" );
    // The refit's failure modes (see `RelActCalcAuto::solve`), which leave the first fit's results.
    for( const string &warning : sol.m_warnings )
      BOOST_CHECK_MESSAGE( !SpecUtils::icontains( warning, "replay" )
                           && !SpecUtils::icontains( warning, "ROIs sized from the fit peak shape" )
                           && !SpecUtils::icontains( warning, "Could not re-size the ROIs" ),
                           "Unexpected warning: " << warning );
    return sol;
  };

  RelActCalcAuto::Options options = with_line_anchored_rois( problem.options );
  options.set_skew_type( PeakDef::SkewType::ExpGaussExp );
  options.lorentzian_xrays = false;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;

  // An energy-dependent skew, carried through the refit.
  {
    options.auto_simplify_model = false;
    const RelActCalcAuto::RelActAutoSolution sol = solve( options );
    BOOST_TEST_MESSAGE( "Energy-dependent skew: " << sol.m_num_roi_refits << " ROI refit(s)" );
    BOOST_REQUIRE_GT( sol.m_num_roi_refits, size_t(0) );
    check_reported_skew_is_the_fitted_line( sol, *problem.foreground );
    for( const RelActCalcAuto::RelActAutoSolution::SkewParResult &par : sol.skew_parameters() )
      BOOST_CHECK( par.was_fit && !par.simplified && !par.energy_dependence_simplified );
  }

  // Every skew removal forced through on the first fit, and replayed on the refit.
  {
    options.auto_simplify_model = true;
    options.auto_simplify_max_dchi2 = 1.0e9;
    const RelActCalcAuto::RelActAutoSolution sol = solve( options );
    BOOST_TEST_MESSAGE( "Simplified skew: " << sol.m_num_roi_refits << " ROI refit(s)" );
    BOOST_REQUIRE_GT( sol.m_num_roi_refits, size_t(0) );

    const PeakDef::SkewType type = options.skew_type;
    for( size_t i = 0; i < PeakDef::num_skew_parameters( type ); ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
      double no_skew = 0.0;
      BOOST_REQUIRE( PeakDef::skew_no_skew_value( type, ct, no_skew ) );
      for( const PeakDef &peak : sol.m_fit_peaks_in_spectrums_cal )
        BOOST_CHECK_MESSAGE( fabs( peak.coefficient(ct) - no_skew ) < 1.0E-9, "Peak at " << peak.mean()
                             << " keV: skew coefficient " << i << " is " << peak.coefficient(ct)
                             << " after the refit, not its no-skew value " << no_skew );
    }
    for( const RelActCalcAuto::RelActAutoSolution::SkewParResult &par : sol.skew_parameters() )
      BOOST_CHECK( par.simplified && !par.clamped );
  }
}//BOOST_AUTO_TEST_CASE( skew_survives_roi_refits )
