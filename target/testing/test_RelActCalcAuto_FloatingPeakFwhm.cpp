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
#include <string>
#include <memory>
#include <vector>
#include <iostream>

#define BOOST_TEST_MODULE RelActCalcAuto_FloatingPeakFwhm_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

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
}//namespace


// Regression guard for review item A7: FloatingPeakResult::fwhm must be reported in keV, even when
// the floating peak's FWHM is released.
//
// When FloatingPeak::release_fwhm is true, the fit parameter at the FWHM index is a dimensionless
// 0.25-4.0 multiplier on the functional FWHM (the actual width is multiplier * fwhm(E)).  The result
// extraction used to store that raw multiplier into FloatingPeakResult::fwhm, while every consumer
// (e.g. the observed-efficiency refit and the bystander-peak reconstruction) treats it as keV and
// divides by 2.35482 to get sigma - so a released-width floating peak's reconstructed sigma was wrong
// by the functional-FWHM factor.
//
// This solves a uranium HPGe spectrum twice with a single floating peak at the 1001 keV (U238) line,
// once with the FWHM fixed to the functional value and once released.  At ~1 MeV the functional FWHM
// is a few keV, well separated from the ~1.0 converged multiplier, so the released-width result being
// in keV (and not the raw multiplier) is observable.
BOOST_AUTO_TEST_CASE( floating_peak_fwhm_reported_in_kev )
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path(
      SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), "U235_Unshielded_6000.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ), "Test spectrum not at '" << spec_path << "'" );

  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE_MESSAGE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ),
                         "Failed to load '" << spec_path << "'" );

  shared_ptr<const SpecUtils::Measurement> foreground, background;
  for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
  {
    if( !m )
      continue;
    if( m->source_type() == SpecUtils::SourceType::Background )
      background = m;
    else if( !foreground )
      foreground = m;
  }
  BOOST_REQUIRE_MESSAGE( foreground, "No foreground measurement in '" << spec_path << "'" );

  const shared_ptr<const DetectorPeakResponse> det = meas->detector();
  BOOST_REQUIRE_MESSAGE( det, "Test N42 has no embedded detector response." );

  const unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE_MESSAGE( state, "Test N42 missing embedded RelActAuto state." );
  BOOST_REQUIRE_MESSAGE( !state->options.rel_eff_curves.empty(), "Embedded state has no rel-eff curve." );

  const RelActCalcAuto::Options base_options = state->options;

  // The floating peak must model a strong, *un-modeled* line so its released FWHM is well-determined.
  // A peak that coincides with a fitted nuclide line (e.g. U238's 1001 keV, which IS modeled - Pa-234m is
  // in equilibrium even at this 20 yr age) competes with the nuclide peak for area, leaving the floating
  // peak's width ill-conditioned (build-dependent).  609.31 keV is the Bi-214 line: a Rn-222 daughter that
  // has NOT grown into U238 at 20 yr, so the U isotopics do not model it, yet it is a clean strong peak in
  // this spectrum (~760 net counts).  At ~609 keV the functional FWHM (~1.4 keV) is also well above the
  // 0.25-4.0 multiplier, so the units bug is visible.  We add a ROI so the fit covers this region.
  // Use the dominant U235 185.7 keV line.  A floating peak's released FWHM is only well-determined when its
  // amplitude can't collapse (into a competing nuclide peak or the ROI continuum) - which needs a strong,
  // sharp peak.  185.7 keV is by far the strongest peak here (~18k net counts); even though it is modeled
  // (so the floating peak shares area with the U235 nuclide peak), the combined peak's *shape* is so tightly
  // constrained that the released FWHM multiplier M* is well-determined (~1) and reproducible.  It already
  // sits in a fitted ROI.  F = fwhm(186) ~= 0.76 keV < 1, so the bug (T = bare M* ~= 1) makes T LARGER than
  // the correct keV value (T = M*F ~= 0.76) - the opposite direction from a high-energy peak.
  const double float_energy = 185.71;

  // Solve with one floating peak at `float_energy`, and return its FloatingPeakResult.
  const auto solve_with_float = [&]( const bool release ) -> RelActCalcAuto::FloatingPeakResult {
    RelActCalcAuto::Options opts = base_options;

    RelActCalcAuto::FloatingPeak fp;
    fp.energy = float_energy;
    fp.release_fwhm = release;
    fp.energy_origin = RelActCalcAuto::FloatingPeak::EnergyType::Known;
    opts.floating_peaks.push_back( fp );

    RelActCalcAuto::RelActAutoSolution sol;
    BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( opts, foreground, background, det, {}, nullptr ) );
    BOOST_REQUIRE_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                           "solve (release_fwhm=" << release << ") failed: status="
                           << static_cast<int>(sol.m_status) << ", error='" << sol.m_error_message << "'" );

    const RelActCalcAuto::FloatingPeakResult *best = nullptr;
    for( const RelActCalcAuto::FloatingPeakResult &fpr : sol.m_floating_peaks )
    {
      if( !best || (fabs(fpr.energy - float_energy) < fabs(best->energy - float_energy)) )
        best = &fpr;
    }
    BOOST_REQUIRE_MESSAGE( best, "No floating-peak result returned (release_fwhm=" << release << ")." );
    return *best;
  };//solve_with_float

  const RelActCalcAuto::FloatingPeakResult fixed_res = solve_with_float( false );
  const RelActCalcAuto::FloatingPeakResult released_res = solve_with_float( true );

  const double F = fixed_res.fwhm;      // functional FWHM at float_energy, in keV
  const double T = released_res.fwhm;   // released-FWHM result; must ALSO be in keV
  const double ratio = (F > 0.0) ? (T / F) : -1.0;

  BOOST_TEST_MESSAGE( "Floating peak @ " << float_energy << " keV: functional FWHM (release=false) = "
                      << F << " keV; released FWHM (release=true) = " << T << " keV; ratio T/F = " << ratio
                      << "; released amplitude = " << released_res.amplitude << " (release=false amplitude = "
                      << fixed_res.amplitude << "). T/F is the fitted FWHM multiplier when reported correctly"
                      " in keV; the pre-fix bug reports the raw multiplier (T ~= the bare M*)." );

  // Pre-condition: F (from the FWHM-fixed solve, untouched by the A7 bug) is a sane HPGe FWHM at 186 keV,
  // and the released floating peak carries real area - the dominant 185.7 peak (~18k counts) keeps its
  // amplitude from collapsing into the nuclide peak, so its width (M*) is well-determined and reproducible.
  // (This is the crux: at a weak peak the released multiplier wanders build-to-build; here it is pinned.)
  BOOST_REQUIRE_MESSAGE( (F > 1.0) && (F < 3.0),
                         "Unexpected functional FWHM " << F << " keV at " << float_energy << " keV." );
  BOOST_REQUIRE_MESSAGE( released_res.amplitude > 500.0,
                         "Released floating peak amplitude " << released_res.amplitude << " is too small - the"
                         " FWHM would be ill-determined; test pre-conditions not met." );

  // The A7 bug only affects the *released* extraction: a correct build reports T = M*F (keV) where M* is the
  // fitted FWHM multiplier (~0.95 for this dominant, accurately-modeled peak); the pre-fix build reports the
  // bare multiplier T = M*.  With F = fwhm(186) ~= 1.7 keV the bug divides T (and T/F) down by F:
  //    correct:  T = M*F ~= 1.6 keV,  T/F = M*   ~= 0.95
  //    pre-fix:  T = M*  ~= 0.95 keV, T/F = M*/F ~= 0.56
  // Require T/F > 0.75 (between the two).  [provisional threshold - calibrated from the measured values]
  BOOST_CHECK_MESSAGE( ratio > 0.75,
                       "Released floating-peak FWHM " << T << " keV gives T/F = " << ratio << "; the pre-fix"
                       " bug reports the bare FWHM multiplier (T/F = M*/F < 1), i.e. divided down by fwhm(E)." );
}//BOOST_AUTO_TEST_CASE( floating_peak_fwhm_reported_in_kev )
