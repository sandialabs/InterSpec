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

#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

#define BOOST_TEST_MODULE test_RelActCalcAuto_GradientCheck
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

static string g_test_file_dir;


static void set_data_dir()
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

  // Best-effort fallback if the args weren't passed (running outside ctest).
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
    for( const char *d : { "test_data", "../test_data", "../../test_data", "../../target/testing/test_data" } )
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
  BOOST_REQUIRE_MESSAGE( db && db->nuclide("U238"), "Error initing SandiaDecayDataBase" );
}//void set_data_dir()


/** Solves `n42_basename` (from `test_data/RelActAutoReport/`, with its embedded GUI state and DRF)
 using `cal_type`, with solve_ceres's derivative oracle (RELACT_GRADIENT_CHECK) turned on.  The
 oracle throws on a significant auto-diff vs finite-difference disagreement, at the start point and
 at convergence.

 A non-zero `cubic_term` re-expresses the file's linear polynomial calibration as {c0, c1, 0,
 cubic_term}: practically the same energies, but four terms, so channels are found numerically
 (`bracketed_newton_root`) rather than analytically.
 */
static void solve_with_gradient_check( const string &n42_basename,
                                       const RelActCalcAuto::EnergyCalFitType cal_type,
                                       const float cubic_term = 0.0f )
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path(
      SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), n42_basename );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ), "Test spectrum not at '" << spec_path << "'" );

  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ) );

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
  BOOST_REQUIRE( foreground );
  BOOST_REQUIRE( meas->detector() );

  if( cubic_term != 0.0f )
  {
    const shared_ptr<const SpecUtils::EnergyCalibration> orig = foreground->energy_calibration();
    BOOST_REQUIRE( orig && (orig->type() == SpecUtils::EnergyCalType::Polynomial) );
    vector<float> coefs = orig->coefficients();
    BOOST_REQUIRE_EQUAL( coefs.size(), 2 );
    coefs.push_back( 0.0f );
    coefs.push_back( cubic_term );

    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_polynomial( orig->num_channels(), coefs, orig->deviation_pairs() );
    BOOST_REQUIRE_EQUAL( cal->coefficients().size(), 4 );

    auto fg = make_shared<SpecUtils::Measurement>( *foreground );
    fg->set_energy_calibration( cal );
    foreground = fg;
    if( background )
    {
      auto bg = make_shared<SpecUtils::Measurement>( *background );
      bg->set_energy_calibration( cal );
      background = bg;
    }
  }//if( cubic_term != 0.0f )

  unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE( state );
  state->options.energy_cal_type = cal_type;

  const PeakFitUtils::CoarseResolutionType det_type
      = PeakFitUtils::effective_det_type( meas->peakFitDetPrefs(), foreground, meas );

#if( defined(_WIN32) )
  _putenv_s( "RELACT_GRADIENT_CHECK", "1" );
#else
  setenv( "RELACT_GRADIENT_CHECK", "1", 1 );
#endif

  RelActCalcAuto::RelActAutoSolution sol;
  string error_msg;
  try
  {
    sol = RelActCalcAuto::solve( state->options, foreground, background, meas->detector(), {}, det_type, nullptr );
  }catch( std::exception &e )
  {
    error_msg = e.what();
  }

#if( defined(_WIN32) )
  _putenv_s( "RELACT_GRADIENT_CHECK", "" );
#else
  unsetenv( "RELACT_GRADIENT_CHECK" );
#endif

  BOOST_CHECK_MESSAGE( error_msg.empty(), "solve threw: " << error_msg );
  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Solve did not succeed (status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "')" );
}//solve_with_gradient_check(...)


// The energy-offset column used to fail the oracle at the start point of this HPGe U fit (up to
//  60% "disagreement" on the 25% coax IDB spectra): `PeakFit::fit_continuum`'s 0.25-sigma
//  continuum-target hinge bends the residuals over ~1E-3 keV of offset, far below the oracle's
//  finite-difference step.  The auto-diff derivative was right; the oracle now refines its step
//  before judging.  Without PERFORM_DEVELOPER_CHECKS the oracle is compiled out, and this is just a
//  fit with the energy calibration fitted.
BOOST_AUTO_TEST_CASE( FittedEnergyCal_NonLinear )
{
  solve_with_gradient_check( "U235_Unshielded_6000.n42", RelActCalcAuto::EnergyCalFitType::NonLinearFit );
}


BOOST_AUTO_TEST_CASE( FittedEnergyCal_Linear )
{
  solve_with_gradient_check( "U235_Unshielded_6000.n42", RelActCalcAuto::EnergyCalFitType::LinearFit );
}


// A four-term calibration is inverted numerically.  That was a bisection stopping 0.001 keV short of
//  the root, which made the calibration mapping a staircase in the energy-calibration parameters: the
//  oracle saw ~0 where auto-diff gave ~50 (1312 bad entries), and with LinearFit the scalar and
//  auto-diff residuals differed.  The 1E-16 cubic term moves the top of this spectrum by ~5E-5 keV.
BOOST_AUTO_TEST_CASE( FittedEnergyCal_FourTermCalibration )
{
  solve_with_gradient_check( "U235_Unshielded_6000.n42", RelActCalcAuto::EnergyCalFitType::NonLinearFit, 1.0e-16f );
  solve_with_gradient_check( "U235_Unshielded_6000.n42", RelActCalcAuto::EnergyCalFitType::LinearFit, 1.0e-16f );
}
