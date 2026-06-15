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

#include <set>
#include <cmath>
#include <string>
#include <memory>
#include <utility>
#include <iostream>

#define BOOST_TEST_MODULE RelActCalcAuto_MassFracConstraint_suite
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


// Regression guard for the mass-fraction-constraint setup path in RelActCalcAuto::solve (review item A5).
//
// Before this test there was NO automated coverage of mass-fraction constraints, and the path was in fact
// broken in developer (assert-enabled) builds: the sum>1 feasibility rescale decoded the activity parameter
// with `- 0.5` instead of `- sm_activity_par_offset`, so any range constraint tripped the consistency assert
// in solve_ceres.  This test loads a uranium HPGe spectrum, adds a U235 [0,1] mass-fraction constraint
// (which routes through that rescale), solves, and checks the fit succeeds and decodes a sane, in-window
// mass fraction.  (The A5 "0.5*frac" start hack that prompted this work was found to be dead code - a second
// setup pass overwrites the constrained start before the solve - and was removed; the cross-corpus evidence
// lives in scratch/a5/.)
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_solves )
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

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u235 );

  // The embedded curve fits U232/U234/U235/U238, so U235 is present with other free U isotopes - a
  // mass-fraction constraint on it is valid (>=2 isotopes of the element, >=1 unconstrained).  A full
  // [0,1] window pushes the post-setup mass fraction to 1.0, exercising the sum>1 feasibility rescale.
  RelActCalcAuto::Options options = state->options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
  constraint.nuclide = u235;
  constraint.lower_mass_fraction = 0.0;
  constraint.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( constraint );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {}, nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Constrained solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 50.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );

  // The constrained nuclide must decode to a sane mass fraction inside (a hair outside, for round-off)
  // its window.  U235 in this sample is well above the lower bound, so also require it strictly positive.
  const pair<double,optional<double>> mf = sol.mass_enrichment_fraction( u235, 0 );
  BOOST_CHECK_MESSAGE( (mf.first > 0.0) && (mf.first <= constraint.upper_mass_fraction + 1.0e-6),
                       "Decoded U235 mass fraction " << mf.first << " is outside the constraint window ["
                       << constraint.lower_mass_fraction << ", " << constraint.upper_mass_fraction << "]." );
}//BOOST_AUTO_TEST_CASE( mass_fraction_constraint_solves )


// Companion: an element carrying TWO range mass-fraction constraints (U234 and U235).  Exercises the
// multi-constraint setup across both nuclide passes (each constrained nuclide set up exactly once) and the
// per-curve decode/extraction for more than one constrained nuclide.  (A fixed lower==upper constraint is
// intentionally NOT used here: RelActCalcManual's fixed-mass-fraction handling has contradictory debug
// asserts and is a separate pre-existing issue - see the A5 entry of RelActCalcAuto_review_2026-06.md.)
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_two_ranges )
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

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u234 = db->nuclide( "U234" );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u234 && u235 );

  // Range-constrain both U234 (to [0,0.05], comfortably bracketing its true ~1.5 wt%) and U235 (to [0,1]);
  // U238 and U232 stay free (the element keeps >=1 unconstrained nuclide).  Lower sum 0 < 1, so feasible.
  RelActCalcAuto::Options options = state->options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint range_u234;
  range_u234.nuclide = u234;
  range_u234.lower_mass_fraction = 0.0;
  range_u234.upper_mass_fraction = 0.05;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint range_u235;
  range_u235.nuclide = u235;
  range_u235.lower_mass_fraction = 0.0;
  range_u235.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( range_u234 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( range_u235 );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {}, nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Two-range constrained solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 50.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );

  // Both constrained nuclides must decode to a value inside (a hair outside, for round-off) their windows.
  const pair<double,optional<double>> mf234 = sol.mass_enrichment_fraction( u234, 0 );
  BOOST_CHECK_MESSAGE( (mf234.first >= -1.0e-6) && (mf234.first <= 0.05 + 1.0e-6),
                       "U234 mass fraction " << mf234.first << " is outside its window [0, 0.05]." );
  const pair<double,optional<double>> mf235 = sol.mass_enrichment_fraction( u235, 0 );
  BOOST_CHECK_MESSAGE( (mf235.first > 0.0) && (mf235.first <= 1.0 + 1.0e-6),
                       "U235 mass fraction " << mf235.first << " is outside [0,1]." );
}//BOOST_AUTO_TEST_CASE( mass_fraction_constraint_two_ranges )


// A FIXED (lower==upper) mass-fraction constraint alongside a range one on the same element.  This is the
// regression test for the RelActCalcManual fix: a fixed mass-fraction constraint removes its nuclide from the
// manual pre-solver's free fit, and the manual solver previously left it with the act-ratio -1.0 norm sentinel
// while its own asserts expected 1.0 - contradictory invariants that aborted assert-enabled builds.  Fixed
// mass-fraction nuclides now take norm 1.0 (matching range ones and release behavior), so this solves.
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_fixed_plus_range )
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

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u234 = db->nuclide( "U234" );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u234 && u235 );

  // Fix U234 at 1.5 wt% (lower==upper) and range-constrain U235 to [0,1]; U238 and U232 stay free.
  RelActCalcAuto::Options options = state->options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint fixed_u234;
  fixed_u234.nuclide = u234;
  fixed_u234.lower_mass_fraction = 0.015;
  fixed_u234.upper_mass_fraction = 0.015;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint range_u235;
  range_u235.nuclide = u235;
  range_u235.lower_mass_fraction = 0.0;
  range_u235.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( fixed_u234 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( range_u235 );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {}, nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Fixed+range constrained solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 50.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );

  // The fixed nuclide must decode back to its pinned mass fraction; the range nuclide must stay in-window.
  const pair<double,optional<double>> mf234 = sol.mass_enrichment_fraction( u234, 0 );
  BOOST_CHECK_MESSAGE( fabs(mf234.first - 0.015) < 1.0e-3,
                       "Fixed U234 mass fraction " << mf234.first << " is not pinned near 0.015." );
  const pair<double,optional<double>> mf235 = sol.mass_enrichment_fraction( u235, 0 );
  BOOST_CHECK_MESSAGE( (mf235.first > 0.0) && (mf235.first <= 1.0 + 1.0e-6),
                       "Range U235 mass fraction " << mf235.first << " is outside [0,1]." );
}//BOOST_AUTO_TEST_CASE( mass_fraction_constraint_fixed_plus_range )
