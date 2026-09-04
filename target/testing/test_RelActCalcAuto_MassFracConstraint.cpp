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
#include <random>
#include <string>
#include <memory>
#include <utility>
#include <iostream>
#include <algorithm>

#include <ceres/jet.h>

// InterSpec/RelActCalc_imp.hpp (included below) transitively pulls in InterSpec/GammaInteractionCalc.h
//  -> <boost/asio>, i.e. winsock2.h; while <boost/test/included/unit_test.hpp> includes <windows.h>,
//  which drags in the legacy WinSock.h.  Pulling winsock2.h in first keeps that later <windows.h>
//  from including WinSock.h -- otherwise MSVC fails with C1189 "WinSock.h has already been included".
//  Same approach as test_DetectorPeakResponse.cpp and test_ShieldingSourceFitCalc.cpp.
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

#define BOOST_TEST_MODULE RelActCalcAuto_MassFracConstraint_suite
#include <boost/test/included/unit_test.hpp>

// Pulls in the sigma-block mass-fraction machinery (RelActCalc::qmax_hinge / MassFracBlockSpec /
//  decode_mass_frac_block / invert_mass_frac_block) for the pure-math tests below; included after
//  <ceres/jet.h> and the std headers (per the header's own requirement) so it can be instantiated
//  as both double and ceres::Jet.
#include "InterSpec/RelActCalc_imp.hpp"

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

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


// Regression guard for an auto-simplify + peak-skew interaction in RelActCalcAuto::solve.
//
// The auto-simplify pass records every parameter it fixes into `as_fixed`, then builds the final Ceres
// manifold from `constant_parameters + as_fixed`.  The accept-branch used to add a fixed parameter even
// when it was already constant from problem setup, so the two lists overlapped and the concatenation
// contained duplicate indices - `ceres::SubsetManifold` then aborts ("The set of constant parameters
// cannot contain duplicates").
//
// Double Sided Crystal Ball is the reliable trigger: its two `n` exponents are never energy dependent
// (PeakDef::is_energy_dependent), so their second-energy-set copies are fixed unconditionally at setup;
// the accepted skew-removal candidate then re-added those already-constant indices.  This test loads a
// uranium HPGe spectrum, switches to DSCB skew with auto-simplify on, and (by forcing every converged
// removal to be accepted) requires the solve to complete - pre-fix it aborts the process at SubsetManifold.
BOOST_AUTO_TEST_CASE( auto_simplify_double_crystal_ball_no_duplicate_const )
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

  RelActCalcAuto::Options options = state->options;
  options.skew_type = PeakDef::SkewType::DoubleSidedCrystalBall;
  options.auto_simplify_model = true;
  // Force every converged removal to be accepted (regardless of chi2 cost) so the skew-removal candidate
  //  is guaranteed to run its accept-branch - that is the bookkeeping path that produced the duplicate
  //  constant index.  Whether real data would prefer the tail is irrelevant to the index bug under test;
  //  this just makes the regression deterministic instead of depending on the spectrum's skew preference.
  options.auto_simplify_max_dchi2 = 1.0e9;
  // `lorentzian_xrays` requires NoSkew/GaussPlusBortel; clear it so DSCB doesn't trip the setup check
  //  (an unrelated exception) before we reach the auto-simplify path under test.
  options.lorentzian_xrays = false;

  RelActCalcAuto::RelActAutoSolution sol;
  // Pre-fix this call aborts the process inside ceres::SubsetManifold; post-fix it returns normally.
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {},
                                    PeakFitUtils::coarse_det_type( foreground, nullptr ), nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "DSCB + auto-simplify solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );

  // Confirm the buggy code path was actually exercised: with removals force-accepted, auto-simplify must
  //  have dropped the skew and recorded it in the warnings (this is the path that produced the duplicate
  //  constant crash).
  const bool removed_skew = std::any_of( begin(sol.m_warnings), end(sol.m_warnings),
                                         []( const string &w ){ return SpecUtils::icontains( w, "peak skew" ); } );
  BOOST_CHECK_MESSAGE( removed_skew,
                       "Expected an auto-simplify 'peak skew' removal warning (the path that produced the"
                       " duplicate-constant crash); none found - the regression path may not be exercised." );

  BOOST_REQUIRE_EQUAL( sol.m_parameter_fixed_by_model_selection.size(),
                       sol.m_final_parameters.size() );
  BOOST_REQUIRE_EQUAL( sol.m_parameter_were_fit.size(),sol.m_final_parameters.size() );
  size_t num_model_fixed = 0;
  for( size_t index = 0; index < sol.m_final_parameters.size(); ++index )
  {
    if( !sol.m_parameter_fixed_by_model_selection[index] )
      continue;
    ++num_model_fixed;
    BOOST_CHECK( !sol.m_parameter_were_fit[index] );
    BOOST_CHECK( std::isfinite(sol.m_final_parameters[index]) );
  }
  BOOST_CHECK_GT( num_model_fixed,0U );
  BOOST_CHECK_NE( sol.m_frozen_model_policy_hash,UINT64_C(0) );
}//BOOST_AUTO_TEST_CASE( auto_simplify_double_crystal_ball_no_duplicate_const )


// =====================================================================================================
// Sigma-block (RelActCalc::MassFracBlockSpec / decode_mass_frac_block) - pure-math unit tests
//  (fast; no spectra/solve).  These replaced the former soft-cap tests when the soft-cap decode -
//  whose 0.95-of-budget knee made fractions above ~0.98 unreachable even for a single [0,1]
//  constraint - was replaced by the exact sigma-block reparameterization.
// =====================================================================================================
namespace
{
  // Decodes a block built from `windows` at the given sigma/g values (doubles).
  vector<double> decode_block( const RelActCalc::MassFracBlockSpec &spec,
                               const double sigma, const vector<double> &gs )
  {
    const size_t num_range = spec.lower.size();
    BOOST_REQUIRE( gs.size() == ((num_range > 1) ? (num_range - 1) : size_t(0)) );
    vector<double> fractions( num_range, 0.0 );
    RelActCalc::decode_mass_frac_block( spec, sigma, gs.data(), fractions.data() );
    return fractions;
  }
}//namespace


// B.1 - The quadratic hinges: one-sided bound properties (qmax >= max with excess <= r/4; qmin the
//  mirror), exactness outside the blend zone, and value+slope continuity at the seams (Jet vs FD).
BOOST_AUTO_TEST_CASE( sigma_block_hinge_properties )
{
  typedef ceres::Jet<double,1> Jet;

  const double a = 0.3, r = 1.0e-3;
  for( int k = -300; k <= 300; ++k )
  {
    const double x = a + k*1.0e-5;  //spans well past both seams
    const double qx = RelActCalc::qmax_hinge( x, a, r );
    const double true_max = (std::max)( x, a );

    BOOST_CHECK_MESSAGE( qx >= true_max - 1.0e-15, "qmax below true max at x=" << x );
    BOOST_CHECK_MESSAGE( qx <= true_max + 0.25*r + 1.0e-15, "qmax excess above r/4 at x=" << x );
    if( fabs(x - a) >= r )
      BOOST_CHECK_SMALL( qx - true_max, 1.0e-15 );  //exact outside the blend zone

    const double qn = RelActCalc::qmin_hinge( x, a, r );
    const double true_min = (std::min)( x, a );
    BOOST_CHECK_MESSAGE( qn <= true_min + 1.0e-15, "qmin above true min at x=" << x );
    BOOST_CHECK_MESSAGE( qn >= true_min - 0.25*r - 1.0e-15, "qmin deficit below -r/4 at x=" << x );

    // AD slope matches central finite differences (incl. across the seams)
    const Jet qj = RelActCalc::qmax_hinge( Jet(x, 0), a, r );
    const double h = 1.0e-9;
    const double fd = (RelActCalc::qmax_hinge( x + h, a, r ) - RelActCalc::qmax_hinge( x - h, a, r ))/(2.0*h);
    BOOST_CHECK_MESSAGE( fabs(qj.v[0] - fd) < 1.0e-5, "qmax AD/FD mismatch at x=" << x );
  }
}


// B.2 - Single range constraint (the dominant real usage): NO smoothing at all - `f == sigma` exactly
//  over the whole box.  Includes the reachability regression: a [0,1] window must decode to
//  1 - delta ~ 0.999999 at the top of the box (the old soft-cap capped out at ~0.9816).
BOOST_AUTO_TEST_CASE( sigma_block_m1_identity_and_reachability )
{
  // A [0,1] window with an unconstrained partner (mixed element)
  const RelActCalc::MassFracBlockSpec spec
                       = RelActCalc::make_mass_frac_block_spec( { {0.0, 1.0} }, 0.0, false );

  BOOST_CHECK_SMALL( spec.sig_lo, 1.0e-15 );
  BOOST_CHECK_CLOSE( spec.sig_hi, 1.0 - RelActCalc::ns_mass_frac_min_remainder_frac, 1.0e-9 );

  for( const double t : { 0.0, 1.0e-6, 0.25, 0.5, 0.75, 0.95, 0.999, 1.0 } )
  {
    const double sigma = spec.sig_lo + t*(spec.sig_hi - spec.sig_lo);
    const vector<double> f = decode_block( spec, sigma, {} );
    BOOST_REQUIRE_EQUAL( f.size(), size_t(1) );
    BOOST_CHECK_SMALL( f[0] - sigma, 1.0e-15 );  //identity - zero distortion anywhere in the box
  }

  // The D1 regression: the top of the window is reachable to within delta
  const vector<double> f_top = decode_block( spec, spec.sig_hi, {} );
  BOOST_CHECK_MESSAGE( f_top[0] >= 0.999,
                       "Single [0,1] constraint only reaches " << f_top[0]
                       << " - the reachability defect (D1) is back." );

  // A narrower window is exact at BOTH ends (no truncation when Sum upper < 1 - delta)
  const RelActCalc::MassFracBlockSpec spec2
                       = RelActCalc::make_mass_frac_block_spec( { {0.2, 0.7} }, 0.0, false );
  BOOST_CHECK_CLOSE( spec2.sig_lo, 0.2, 1.0e-12 );
  BOOST_CHECK_CLOSE( spec2.sig_hi, 0.7, 1.0e-12 );
}


// B.3 - Exactness fuzz over random polytopes: Sum fractions == sigma exactly, every lower bound exact,
//  every upper bound exact to ~width_radius/4 (~2.5e-13), for random in-box parameters.
BOOST_AUTO_TEST_CASE( sigma_block_exactness_fuzz )
{
  std::mt19937 rng( 20260701u );  //deterministic
  std::uniform_real_distribution<double> unit( 0.0, 1.0 );

  for( int trial = 0; trial < 500; ++trial )
  {
    const size_t num_range = 2 + (rng() % 4);            //2..5 windows
    const double fixed_sum = (trial % 3) ? 0.0 : 0.15;   //sometimes a fixed constraint too

    // Random windows, scaled so their lower sum stays well below the budget
    vector<pair<double,double>> windows( num_range );
    double lower_sum = 0.0;
    for( size_t k = 0; k < num_range; ++k )
    {
      const double l = 0.5*unit(rng)/num_range;
      const double u = l + (1.0e-4 + unit(rng))/(1.0 + 2.0*unit(rng)*num_range);
      windows[k] = { l, (std::min)(u, 1.0) };
      lower_sum += l;
    }
    if( (fixed_sum + lower_sum) >= 0.95 )
      continue;  //keep clearly feasible

    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( windows, fixed_sum, false );

    const double sigma = spec.sig_lo + unit(rng)*(spec.sig_hi - spec.sig_lo);
    vector<double> gs( num_range - 1 );
    for( double &g : gs )
      g = unit(rng);

    const vector<double> f = decode_block( spec, sigma, gs );

    double f_sum = 0.0;
    for( size_t k = 0; k < num_range; ++k )
    {
      f_sum += f[k];
      BOOST_CHECK_MESSAGE( f[k] >= (windows[k].first - 1.0e-14),
                           "trial " << trial << ": f[" << k << "]=" << f[k]
                           << " below lower bound " << windows[k].first );
      BOOST_CHECK_MESSAGE( f[k] <= (windows[k].second + 1.0e-11),
                           "trial " << trial << ": f[" << k << "]=" << f[k]
                           << " above upper bound " << windows[k].second );
    }
    BOOST_CHECK_MESSAGE( fabs(f_sum - sigma) < 1.0e-12,
                         "trial " << trial << ": Sum f=" << f_sum << " != sigma=" << sigma );
    BOOST_CHECK( (spec.fixed_sum + spec.sig_hi) <= (1.0 - 0.5*spec.delta) );
  }//for( int trial = 0; trial < 500; ++trial )
}


// B.4 - Corners: sigma at either end of its box, g values at 0/1 - fractions stay in-window and the
//  sum stays exact, including for the minimal m=2 block and a pinched conditional interval.
BOOST_AUTO_TEST_CASE( sigma_block_corners )
{
  const vector<vector<pair<double,double>>> cases = {
    { {0.0, 1.0}, {0.0, 1.0} },                    //two wide-open windows
    { {0.1, 0.4}, {0.2, 0.3} },                    //offset windows
    { {0.0, 0.05}, {0.0, 1.0}, {0.05, 0.2} },      //narrow + wide + medium
  };

  for( size_t c = 0; c < cases.size(); ++c )
  {
    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( cases[c], 0.0, false );
    const size_t num_range = cases[c].size();

    for( const double t : { 0.0, 1.0 } )
    {
      const double sigma = spec.sig_lo + t*(spec.sig_hi - spec.sig_lo);
      for( int g_pattern = 0; g_pattern < (1 << (num_range-1)); ++g_pattern )
      {
        vector<double> gs( num_range - 1 );
        for( size_t k = 0; k+1 < num_range; ++k )
          gs[k] = ((g_pattern >> k) & 1) ? 1.0 : 0.0;

        const vector<double> f = decode_block( spec, sigma, gs );
        double f_sum = 0.0;
        for( size_t k = 0; k < num_range; ++k )
        {
          f_sum += f[k];
          BOOST_CHECK_MESSAGE( (f[k] >= cases[c][k].first - 1.0e-12)
                                && (f[k] <= cases[c][k].second + 1.0e-11),
                               "case " << c << " t=" << t << " pattern=" << g_pattern
                               << ": f[" << k << "]=" << f[k] << " outside window" );
        }
        BOOST_CHECK_SMALL( f_sum - sigma, 1.0e-12 );
      }//for( g_pattern )
    }//for( t in {0,1} )
  }//for( case )
}


// B.5 - Inversion round-trip: interior target fractions -> (sigma, g) -> decode reproduces the targets.
BOOST_AUTO_TEST_CASE( sigma_block_inversion_round_trip )
{
  std::mt19937 rng( 20260702u );
  std::uniform_real_distribution<double> unit( 0.0, 1.0 );

  for( int trial = 0; trial < 200; ++trial )
  {
    const size_t num_range = 2 + (rng() % 3);   //2..4
    vector<pair<double,double>> windows( num_range );
    for( size_t k = 0; k < num_range; ++k )
    {
      const double l = 0.4*unit(rng)/num_range;
      const double u = l + (0.05 + 0.5*unit(rng))/num_range;
      windows[k] = { l, u };
    }

    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( windows, 0.0, false );

    // Interior targets: decode a random interior point, use the result as the target
    const double sigma_in = spec.sig_lo + (0.1 + 0.8*unit(rng))*(spec.sig_hi - spec.sig_lo);
    vector<double> gs_in( num_range - 1 );
    for( double &g : gs_in )
      g = 0.1 + 0.8*unit(rng);
    const vector<double> targets = decode_block( spec, sigma_in, gs_in );

    double sigma_out = 0.0;
    vector<double> gs_out( num_range - 1, 0.0 );
    RelActCalc::invert_mass_frac_block( spec, targets.data(), sigma_out, gs_out.data() );
    const vector<double> f = decode_block( spec, sigma_out, gs_out );

    for( size_t k = 0; k < num_range; ++k )
      BOOST_CHECK_MESSAGE( fabs(f[k] - targets[k]) < 1.0e-9,
                           "trial " << trial << ": round-trip f[" << k << "]=" << f[k]
                           << " vs target " << targets[k] );
  }//for( int trial = 0; trial < 200; ++trial )
}


// B.6 - Tiny-window regression (the U232 case): a window like [0, 0.9e-9] is a legitimate LIVE
//  parameter - the box must not collapse, the decode must span the whole window with a nonzero
//  derivative, alone (m=1) and alongside a wide [0,1] window (m=2, either order).
BOOST_AUTO_TEST_CASE( sigma_block_tiny_window_stays_live )
{
  typedef ceres::Jet<double,2> Jet;

  const double tiny_hi = 0.9e-9;

  {// m == 1: box is exactly the window (no delta truncation: Sum upper << 1 - delta)
    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( { {0.0, tiny_hi} }, 0.0, false );
    BOOST_CHECK_SMALL( spec.sig_lo, 1.0e-30 );
    BOOST_CHECK_CLOSE( spec.sig_hi, tiny_hi, 1.0e-9 );
    // The box is far wider than machine-epsilon-relative, so it must be treated as live:
    BOOST_CHECK( (spec.sig_hi - spec.sig_lo)
                  > 4.0*std::numeric_limits<double>::epsilon()*(std::max)(spec.sig_lo, fabs(spec.sig_hi)) );

    for( const double t : { 0.0, 0.5, 1.0 } )
    {
      const double sigma = spec.sig_lo + t*(spec.sig_hi - spec.sig_lo);
      const vector<double> f = decode_block( spec, sigma, {} );
      BOOST_CHECK_SMALL( f[0] - sigma, 1.0e-24 );  //identity, even at this scale
    }
  }

  {// m == 2: wide carrier [0,1] first, tiny window second (g-slot)
    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( { {0.0, 1.0}, {0.0, tiny_hi} }, 0.0, false );

    const double sigma = 0.3;
    Jet sigma_j( sigma, 0 ), g_j( 0.75, 1 );
    Jet fractions_j[2];
    RelActCalc::decode_mass_frac_block( spec, sigma_j, &g_j, fractions_j );

    // The tiny fraction spans its window and responds to its parameter (nonzero derivative)
    BOOST_CHECK( (fractions_j[1].a >= 0.0) && (fractions_j[1].a <= tiny_hi*(1.0 + 1.0e-6)) );
    BOOST_CHECK_MESSAGE( fabs(fractions_j[1].v[1]) > 0.1*tiny_hi,
                         "tiny-window fraction has (near-)zero derivative w.r.t. its own parameter: "
                         << fractions_j[1].v[1] );
    // ... and the sum stays exact
    BOOST_CHECK_SMALL( (fractions_j[0].a + fractions_j[1].a) - sigma, 1.0e-12 );
  }
}


// B.7 - All-constrained element: sigma is the constant leftover budget, fractions sum to exactly
//  (1 - fixed_sum) - i.e., to exactly 1 across the whole element.
BOOST_AUTO_TEST_CASE( sigma_block_all_constrained_sums_to_one )
{
  std::mt19937 rng( 20260703u );
  std::uniform_real_distribution<double> unit( 0.0, 1.0 );

  for( const double fixed_sum : { 0.0, 0.1 } )
  {
    // Windows that can bracket the 1-fixed_sum total: generous uppers
    const vector<pair<double,double>> windows = { {0.0, 1.0}, {0.0, 0.6}, {0.0, 0.8} };
    const RelActCalc::MassFracBlockSpec spec
                        = RelActCalc::make_mass_frac_block_spec( windows, fixed_sum, true );

    BOOST_CHECK_CLOSE( spec.sig_lo, 1.0 - fixed_sum, 1.0e-12 );
    BOOST_CHECK_CLOSE( spec.sig_hi, 1.0 - fixed_sum, 1.0e-12 );
    BOOST_CHECK_SMALL( spec.delta, 1.0e-30 );

    for( int trial = 0; trial < 50; ++trial )
    {
      vector<double> gs( windows.size() - 1 );
      for( double &g : gs )
        g = unit(rng);

      const vector<double> f = decode_block( spec, spec.sig_hi, gs );
      double f_sum = fixed_sum;
      for( size_t k = 0; k < f.size(); ++k )
      {
        f_sum += f[k];
        BOOST_CHECK( (f[k] >= windows[k].first - 1.0e-12) && (f[k] <= windows[k].second + 1.0e-11) );
      }
      BOOST_CHECK_MESSAGE( fabs(f_sum - 1.0) < 1.0e-12,
                           "all-constrained fractions sum to " << f_sum << " != 1" );
    }
  }//for( fixed_sum )
}


// B.8 - Jet gradient continuity: the decode's AD derivatives match central finite differences on a
//  sweep that crosses the hinge blend zones (the seams are C1 by construction).
BOOST_AUTO_TEST_CASE( sigma_block_jet_matches_finite_difference )
{
  typedef ceres::Jet<double,3> Jet;

  const vector<pair<double,double>> windows = { {0.0, 0.5}, {0.1, 0.4}, {0.0, 0.3} };
  const RelActCalc::MassFracBlockSpec spec
                      = RelActCalc::make_mass_frac_block_spec( windows, 0.0, false );

  const int nsteps = 400;
  for( int i = 0; i <= nsteps; ++i )
  {
    const double t = static_cast<double>(i)/nsteps;
    const double params[3] = { spec.sig_lo + t*(spec.sig_hi - spec.sig_lo),  //sigma sweeps its box
                               0.98,                                         //g_1 near its corner
                               0.5*t };                                      //g_2 sweeps

    Jet sigma_j( params[0], 0 ), gs_j[2] = { Jet(params[1], 1), Jet(params[2], 2) }, f_j[3];
    RelActCalc::decode_mass_frac_block( spec, sigma_j, gs_j, f_j );

    for( int p = 0; p < 3; ++p )
    {
      const double h = 1.0e-7;
      double up[3] = { params[0], params[1], params[2] }, dn[3] = { params[0], params[1], params[2] };
      up[p] += h;  dn[p] -= h;
      double f_up[3], f_dn[3];
      {
        const double g_up[2] = { up[1], up[2] }, g_dn[2] = { dn[1], dn[2] };
        RelActCalc::decode_mass_frac_block( spec, up[0], g_up, f_up );
        RelActCalc::decode_mass_frac_block( spec, dn[0], g_dn, f_dn );
      }

      for( int k = 0; k < 3; ++k )
      {
        const double fd = (f_up[k] - f_dn[k])/(2.0*h);
        BOOST_CHECK_MESSAGE( fabs(f_j[k].v[p] - fd) < 1.0e-4*(std::max)(1.0, fabs(fd)),
                             "AD/FD mismatch at step " << i << " for df[" << k << "]/dp[" << p
                             << "]: ad=" << f_j[k].v[p] << " fd=" << fd );
      }
    }//for( int p = 0; p < 3; ++p )
  }//for( int i = 0; i <= nsteps; ++i )
}


// =====================================================================================================
// A16 integration tests (full solve)
// =====================================================================================================
namespace
{
  // Loads the embedded-state uranium HPGe test case shared by the integration tests below.
  struct ULoadResult
  {
    shared_ptr<SpecMeas> meas;
    shared_ptr<const SpecUtils::Measurement> foreground, background;
    shared_ptr<const DetectorPeakResponse> det;
    RelActCalcAuto::Options options;   // copy of the embedded state's options, ready to tweak
  };

  bool load_u_test_case( ULoadResult &out )
  {
    const string spec_path = SpecUtils::append_path(
        SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), "U235_Unshielded_6000.n42" );
    if( !SpecUtils::is_file( spec_path ) )
      return false;
    out.meas = make_shared<SpecMeas>();
    if( !out.meas->load_file( spec_path, SpecUtils::ParserType::Auto ) )
      return false;
    for( const shared_ptr<const SpecUtils::Measurement> &m : out.meas->measurements() )
    {
      if( !m )
        continue;
      if( m->source_type() == SpecUtils::SourceType::Background )
        out.background = m;
      else if( !out.foreground )
        out.foreground = m;
    }
    out.det = out.meas->detector();
    const unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = out.meas->getRelActAutoGuiState();
    if( !out.foreground || !out.det || !state || state->options.rel_eff_curves.empty() )
      return false;
    out.options = state->options;
    return true;
  }


  /** A small live U problem used by deterministic repeat/permutation tests.  It retains the
   physical-model response and a nonlinear U-235 mass-fraction coordinate, but fixes nuisance
   families and keeps only the U-235/U-238 evidence windows so repeated solves remain inexpensive. */
  RelActCalcAuto::Options small_canonical_u_options(
                                      const RelActCalcAuto::Options &input,
                                      const SandiaDecay::Nuclide * const u235,
                                      const SandiaDecay::Nuclide * const u238 )
  {
    BOOST_REQUIRE( u235 && u238 );
    RelActCalcAuto::Options options = input;
    options.auto_profile_weak_mass_fractions = false;
    options.auto_simplify_model = false;
    options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
    options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
    options.fwhm_estimation_method
        = RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum;
    options.skew_type = PeakDef::SkewType::NoSkew;
    options.additional_br_uncert = 0.0;

    for( RelActCalcAuto::RoiRange &roi : options.rois )
      roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    options.rois.erase( std::remove_if(options.rois.begin(),options.rois.end(),
        []( const RelActCalcAuto::RoiRange &roi ) {
          const bool low_energy_u235 = (roi.lower_energy >= 140.0)
                                         && (roi.upper_energy <= 210.0);
          const bool high_energy_u238 = (roi.lower_energy >= 990.0)
                                         && (roi.upper_energy <= 1010.0);
          return !low_energy_u235 && !high_energy_u238;
        }),options.rois.end() );
    BOOST_REQUIRE_GE( options.rois.size(),3u );

    BOOST_REQUIRE_EQUAL( options.rel_eff_curves.size(),1u );
    RelActCalcAuto::RelEffCurveInput &curve = options.rel_eff_curves[0];
    curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
    curve.rel_eff_eqn_order = 0;
    curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
    curve.shielded_by_other_phys_model_curve_shieldings.clear();
    curve.mass_fraction_constraints.clear();
    curve.act_ratio_constraints.clear();

    const auto freeze_shield = [](
        const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &input_shield ) {
      if( !input_shield )
        return input_shield;
      auto fixed = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input_shield);
      fixed->fit_atomic_number = false;
      fixed->fit_areal_density = false;
      fixed->lower_fit_areal_density = fixed->areal_density;
      fixed->upper_fit_areal_density = fixed->areal_density;
      return std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>(fixed);
    };
    curve.phys_model_self_atten = freeze_shield(curve.phys_model_self_atten);
    for( std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &shield
           : curve.phys_model_external_atten )
      shield = freeze_shield(shield);

    curve.nuclides.erase( std::remove_if(curve.nuclides.begin(),curve.nuclides.end(),
        [u235,u238]( const RelActCalcAuto::NucInputInfo &source ) {
          const SandiaDecay::Nuclide * const nuclide
                                            = RelActCalcAuto::nuclide(source.source);
          return (nuclide != u235) && (nuclide != u238);
        }),curve.nuclides.end() );
    BOOST_REQUIRE_EQUAL( curve.nuclides.size(),2u );
    for( RelActCalcAuto::NucInputInfo &source : curve.nuclides )
    {
      source.fit_age = false;
      source.fit_age_min.reset();
      source.fit_age_max.reset();
      source.force_profile_mass_fraction = false;
    }

    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint u235_window;
    u235_window.nuclide = u235;
    u235_window.lower_mass_fraction = 0.10;
    u235_window.upper_mass_fraction = 0.50;
    curve.mass_fraction_constraints.push_back(u235_window);
    return options;
  }


  /** Restore the process-wide solve-thread limit even when a BOOST_REQUIRE aborts a test case.
   `max_solve_threads()` returns the resolved previous behavior; restoring that explicit value is
   behaviorally equivalent for the lifetime of this test process. */
  class SolveThreadLimitRestorer
  {
  public:
    SolveThreadLimitRestorer()
      : m_previous(RelActCalc::max_solve_threads())
    {}

    ~SolveThreadLimitRestorer()
    {
      RelActCalc::set_max_solve_threads(m_previous);
    }

    SolveThreadLimitRestorer( const SolveThreadLimitRestorer & ) = delete;
    SolveThreadLimitRestorer &operator=( const SolveThreadLimitRestorer & ) = delete;

  private:
    const unsigned m_previous;
  };
}//namespace


// MODEL-03: U-234 is tied to U-235 here, so U-235 controls a same-element isotope and no slot can
// scan its reported mass fraction exactly (deriving U-235's activity from its fraction while
// U-234's activity derives from U-235's would be circular).  The carrier engine must refuse with a
// structured `Failed` carrying that reason - never quote a confidently narrow inexact interval,
// and never silently drop the explicit request.
BOOST_AUTO_TEST_CASE( forced_ratio_constrained_u235_profile )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u234 = db->nuclide("U234");
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  BOOST_REQUIRE( u234 && u235 );

  RelActCalcAuto::Options options = tc.options;
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method
      = RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum;
  options.skew_type = PeakDef::SkewType::NoSkew;
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  options.rois.erase( std::remove_if(options.rois.begin(),options.rois.end(),
      []( const RelActCalcAuto::RoiRange &roi ){
        const bool low_energy_u235 = (roi.lower_energy >= 140.0) && (roi.upper_energy <= 210.0);
        const bool high_energy_u238 = (roi.lower_energy >= 990.0) && (roi.upper_energy <= 1010.0);
        return !low_energy_u235 && !high_energy_u238;
      }),options.rois.end() );
  BOOST_REQUIRE_GE( options.rois.size(),3u );

  RelActCalcAuto::RelEffCurveInput &curve = options.rel_eff_curves.at(0);
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnY;
  curve.rel_eff_eqn_order = 1;
  curve.phys_model_self_atten.reset();
  curve.phys_model_external_atten.clear();
  curve.nuclides.erase( std::remove_if(curve.nuclides.begin(),curve.nuclides.end(),
      []( const RelActCalcAuto::NucInputInfo &input ){
        const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
        return !nuclide || (nuclide->massNumber == 232);
      }),curve.nuclides.end() );
  bool found_u235 = false;
  for( RelActCalcAuto::NucInputInfo &input : curve.nuclides )
  {
    const bool is_u235 = RelActCalcAuto::nuclide(input.source) == u235;
    input.force_profile_mass_fraction = is_u235;
    found_u235 = found_u235 || is_u235;
  }
  BOOST_REQUIRE( found_u235 );

  RelActCalcAuto::RelEffCurveInput::ActRatioConstraint ratio;
  ratio.controlling_source = RelActCalcAuto::SrcVariant(u235);
  ratio.constrained_source = RelActCalcAuto::SrcVariant(u234);
  // Approximate the fixture's 1.5/25 wt% U-234/U-235 composition in activity coordinates.
  ratio.constrained_to_controlled_activity_ratio
      = (0.015*u234->activityPerGram()) / (0.25*u235->activityPerGram());
  curve.act_ratio_constraints.push_back(ratio);

  vector<shared_ptr<const PeakDef>> input_peaks;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> stored_peaks
      = static_cast<const SpecMeas &>(*tc.meas).peaks({tc.foreground->sample_number()});
  if( stored_peaks )
    input_peaks.assign(stored_peaks->begin(),stored_peaks->end());

  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         solution.m_error_message );

  const auto result = solution.mass_enrichment_result(u235,0);
  BOOST_REQUIRE_MESSAGE( result.profile.has_value(),
                         "Explicit U-235 profile request produced no structured result" );
  BOOST_CHECK( result.profile->reason
      == RelActCalcAuto::RelActAutoSolution::MassFractionProfileReason::Forced );
  BOOST_CHECK_MESSAGE( result.profile->status
      == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::Failed,
      "expected a structured refusal, got status "
        << static_cast<int>(result.profile->status) << ": " << result.profile->message );
  BOOST_CHECK_MESSAGE( !result.profile->message.empty(),
                       "a refusal must say why the fraction cannot be scanned exactly" );
  BOOST_CHECK( result.profile->intervals.empty() );
}


// MODEL-03: an explicit request profiles U-235 even when its quantity-specific local covariance is
// already usable.  This is intentionally a small two-isotope empirical problem so the production
// conditional solve (rather than fixture complexity) is what the regression measures.
/** The frozen, cheap U-235/U-238 physical-model configuration the forced-profile cases share.

 Everything that could move between two solves of the "same" problem is pinned (fixed ROIs, no
 energy-cal fit, fixed FWHM, no skew, frozen shields, no Hoerl correction, no auto-simplify), so a
 difference between two solves is attributable to the one thing under test.
 */
RelActCalcAuto::Options u235_forced_profile_options( const ULoadResult &tc,
                                                     const SandiaDecay::Nuclide * const u235,
                                                     const bool force_profile,
                                                     const bool add_constraint = true )
{
  RelActCalcAuto::Options options = tc.options;
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method
      = RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum;
  options.skew_type = PeakDef::SkewType::NoSkew;
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  options.rois.erase( std::remove_if(options.rois.begin(),options.rois.end(),
      []( const RelActCalcAuto::RoiRange &roi ){
        const bool low_energy_u235 = (roi.lower_energy >= 140.0) && (roi.upper_energy <= 210.0);
        const bool high_energy_u238 = (roi.lower_energy >= 990.0) && (roi.upper_energy <= 1010.0);
        return !low_energy_u235 && !high_energy_u238;
      }),options.rois.end() );

  RelActCalcAuto::RelEffCurveInput &curve = options.rel_eff_curves.at(0);
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  curve.rel_eff_eqn_order = 0;
  const auto freeze_shield = []( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &input ) {
    if( !input )
      return input;
    auto fixed = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input);
    fixed->fit_atomic_number = false;
    fixed->fit_areal_density = false;
    fixed->lower_fit_areal_density = fixed->areal_density;
    fixed->upper_fit_areal_density = fixed->areal_density;
    return std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>(fixed);
  };
  curve.phys_model_self_atten = freeze_shield(curve.phys_model_self_atten);
  for( std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &shield
       : curve.phys_model_external_atten )
    shield = freeze_shield(shield);
  curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
  curve.nuclides.erase( std::remove_if(curve.nuclides.begin(),curve.nuclides.end(),
      []( const RelActCalcAuto::NucInputInfo &input ){
        const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
        return !nuclide || ((nuclide->massNumber != 235) && (nuclide->massNumber != 238));
      }),curve.nuclides.end() );
  for( RelActCalcAuto::NucInputInfo &input : curve.nuclides )
    input.force_profile_mass_fraction
        = force_profile && (RelActCalcAuto::nuclide(input.source) == u235);

  if( add_constraint )
  {
    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint healthy_window;
    healthy_window.nuclide = u235;
    healthy_window.lower_mass_fraction = 0.10;
    healthy_window.upper_mass_fraction = 0.50;
    curve.mass_fraction_constraints.push_back(healthy_window);
  }

  return options;
}//u235_forced_profile_options(...)


/** Shared load boilerplate for the U-235 profile cases below. */
struct UProfileFixture
{
  ULoadResult tc;
  const SandiaDecay::Nuclide *u235 = nullptr;
  vector<shared_ptr<const PeakDef>> input_peaks;
};

bool load_u_profile_fixture( UProfileFixture &fixture )
{
  set_data_dir();
  if( !load_u_test_case(fixture.tc) )
    return false;
  fixture.u235 = DecayDataBaseServer::database()->nuclide("U235");
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> stored_peaks
      = static_cast<const SpecMeas &>(*fixture.tc.meas).peaks(
                                              {fixture.tc.foreground->sample_number()});
  if( stored_peaks )
    fixture.input_peaks.assign(stored_peaks->begin(),stored_peaks->end());
  return fixture.u235 != nullptr;
}//load_u_profile_fixture(...)


/** Invariants 4 and 7 of the profile design, checked through the public API.

 A profile pins a parameter with a `ceres::SubsetManifold` on the retained converged
 `ceres::Problem`; nothing is added to the objective.  Two obligations fail silently if broken:

  4. Requesting a profile must not move the fit itself: the solve is identical whether or not the
     converged problem is retained for profiling afterwards.
  7. After profiling, the shared parameter buffer is back at the unconstrained optimum and the
     production manifold is restored, so anything that reads the retained cost functor afterwards
     (uncertainties, reports, the rel-eff chart) describes exactly the point the solution claims.
 */
BOOST_AUTO_TEST_CASE( profile_request_is_inert_and_restores_the_optimum )
{
  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  const auto run = [&]( const bool force_profile ){
    RelActCalcAuto::RelActAutoSolution solution;
    BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
        u235_forced_profile_options(tc,u235,force_profile),
        tc.foreground,tc.background,tc.det,input_peaks,
        PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
    BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                           solution.m_error_message );
    return solution;
  };

  const RelActCalcAuto::RelActAutoSolution without_profile = run( false );
  const RelActCalcAuto::RelActAutoSolution with_profile = run( true );

  // The profile really ran only in the second solve.
  BOOST_REQUIRE( !without_profile.mass_enrichment_result(u235,0).profile.has_value() );
  BOOST_REQUIRE( with_profile.mass_enrichment_result(u235,0).profile.has_value() );

  // Invariant 4.  Nothing about the profile request participates in the minimization, so it must
  // not change the ANSWER.
  BOOST_REQUIRE( std::isfinite(without_profile.m_chi2) && std::isfinite(with_profile.m_chi2) );
  const double chi2_scale = std::max(1.0,std::fabs(without_profile.m_chi2));
  BOOST_CHECK_MESSAGE( std::fabs(with_profile.m_chi2 - without_profile.m_chi2) <= 1.0e-6*chi2_scale,
      "A profile request moved the objective: " << without_profile.m_chi2
        << " -> " << with_profile.m_chi2 );

  const double fraction_without = without_profile.mass_enrichment_result(u235,0).fraction;
  const double fraction_with = with_profile.mass_enrichment_result(u235,0).fraction;
  BOOST_REQUIRE( std::isfinite(fraction_without) && std::isfinite(fraction_with) );

  // Invariant 7.  `mass_enrichment_result` reads the retained cost functor, so a scan that left
  // the shared parameter buffer at its last conditional point would show up right here as a
  // nominal that no longer matches the unprofiled solve.
  BOOST_CHECK_MESSAGE( std::fabs(fraction_with - fraction_without) <= 1.0e-6,
      "Profiling did not restore the reported nominal: " << fraction_without
        << " -> " << fraction_with );

  // A returned solution never carries a live optimizer problem.
  for( const RelActCalcAuto::RelActAutoSolution *solution : { &without_profile,&with_profile } )
    BOOST_CHECK( !solution->m_profile_host );

  // Every sampled point of the scan is a conditional minimum of the SAME physical objective, so
  // none may sit below the unconstrained optimum: a delta-chi2 that came out negative would mean
  // the reported baseline was not the optimum, which the engine must route to the
  // deferred-rebaseline flow instead of reporting.
  const auto profile = *with_profile.mass_enrichment_result(u235,0).profile;
  BOOST_CHECK_GT( profile.sampled_delta_chi2.size(),0u );
  for( const std::pair<double,double> &sample : profile.sampled_delta_chi2 )
  {
    BOOST_CHECK_MESSAGE( std::isfinite(sample.first) && std::isfinite(sample.second),
        "Non-finite profile sample (" << sample.first << "," << sample.second << ")" );
    BOOST_CHECK_MESSAGE( sample.second >= -1.0e-6,
        "Profile sample at fraction " << sample.first << " reported delta-chi2 "
          << sample.second << ", below the unconstrained optimum" );
  }
}


BOOST_AUTO_TEST_CASE( forced_healthy_u235_profile )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  BOOST_REQUIRE( u235 );

  RelActCalcAuto::Options options = tc.options;
  options.auto_profile_weak_mass_fractions = false;
  options.auto_simplify_model = false;
  options.energy_cal_type = RelActCalcAuto::EnergyCalFitType::NoFit;
  options.fwhm_form = RelActCalcAuto::FwhmForm::Polynomial_2;
  options.fwhm_estimation_method
      = RelActCalcAuto::FwhmEstimationMethod::FixedToAllPeaksInSpectrum;
  options.skew_type = PeakDef::SkewType::NoSkew;
  for( RelActCalcAuto::RoiRange &roi : options.rois )
    roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  options.rois.erase( std::remove_if(options.rois.begin(),options.rois.end(),
      []( const RelActCalcAuto::RoiRange &roi ){
        const bool low_energy_u235 = (roi.lower_energy >= 140.0) && (roi.upper_energy <= 210.0);
        const bool high_energy_u238 = (roi.lower_energy >= 990.0) && (roi.upper_energy <= 1010.0);
        return !low_energy_u235 && !high_energy_u238;
      }),options.rois.end() );
  BOOST_REQUIRE_GE( options.rois.size(),3u );

  RelActCalcAuto::RelEffCurveInput &curve = options.rel_eff_curves.at(0);
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
  curve.rel_eff_eqn_order = 0;
  const auto freeze_shield = []( const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &input ) {
    if( !input )
      return input;
    auto fixed = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input);
    fixed->fit_atomic_number = false;
    fixed->fit_areal_density = false;
    fixed->lower_fit_areal_density = fixed->areal_density;
    fixed->upper_fit_areal_density = fixed->areal_density;
    return std::shared_ptr<const RelActCalc::PhysicalModelShieldInput>(fixed);
  };
  curve.phys_model_self_atten = freeze_shield(curve.phys_model_self_atten);
  for( std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &shield
       : curve.phys_model_external_atten )
    shield = freeze_shield(shield);
  curve.phys_model_corr.corr_fcn = RelActCalc::PhysModelCorrFcn::None;
  curve.nuclides.erase( std::remove_if(curve.nuclides.begin(),curve.nuclides.end(),
      []( const RelActCalcAuto::NucInputInfo &input ){
        const SandiaDecay::Nuclide * const nuclide = RelActCalcAuto::nuclide(input.source);
        return !nuclide || ((nuclide->massNumber != 235) && (nuclide->massNumber != 238));
      }),curve.nuclides.end() );
  BOOST_REQUIRE_EQUAL( curve.nuclides.size(),2u );
  bool found_u235 = false;
  for( RelActCalcAuto::NucInputInfo &input : curve.nuclides )
  {
    const bool is_u235 = RelActCalcAuto::nuclide(input.source) == u235;
    input.force_profile_mass_fraction = is_u235;
    found_u235 = found_u235 || is_u235;
  }
  BOOST_REQUIRE( found_u235 );
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint healthy_window;
  healthy_window.nuclide = u235;
  healthy_window.lower_mass_fraction = 0.10;
  healthy_window.upper_mass_fraction = 0.50;
  curve.mass_fraction_constraints.push_back(healthy_window);

  vector<shared_ptr<const PeakDef>> input_peaks;
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> stored_peaks
      = static_cast<const SpecMeas &>(*tc.meas).peaks({tc.foreground->sample_number()});
  if( stored_peaks )
    input_peaks.assign(stored_peaks->begin(),stored_peaks->end());

  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         solution.m_error_message );

  const auto result = solution.mass_enrichment_result(u235,0);
  BOOST_REQUIRE_MESSAGE(
      result.covariance_quality
        == RelActCalcAuto::RelActAutoSolution::MassFractionCovarianceQuality::Usable,
      "Expected a healthy local U-235 covariance before honoring the explicit profile request: "
        << result.diagnostic << " (fraction=" << result.fraction << ")" );
  BOOST_REQUIRE( result.profile.has_value() );
  BOOST_CHECK( result.profile->reason
      == RelActCalcAuto::RelActAutoSolution::MassFractionProfileReason::Forced );
  const auto profile_status = result.profile->status;
  BOOST_REQUIRE_MESSAGE(
      profile_status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::Complete
      || profile_status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::BoundaryLimited
      || profile_status == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::NonIdentifiable,
      result.profile->message );
  BOOST_CHECK_GT( result.profile->num_fits,0u );
  BOOST_CHECK_LE( result.profile->num_fits,32u );
  BOOST_REQUIRE_EQUAL( result.profile->intervals.size(),2u );
}


// A.8 - Overlap regression: an element carrying TWO wide overlapping constraints (U234 [0,1] and
//  U235 [0,1]).  The pre-A16 eval threw "Sum of constrained mass fractions ... > 1.0" whenever a
//  trial step drove the two fractions past the simplex face (a hard Ceres invalid-step wall); the
//  sigma-block makes Σ <= 1 - delta a hard box bound on ONE parameter, so the solve must simply
//  succeed with both fractions feasible.
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_two_wide_overlap )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case( tc ), "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u234 = db->nuclide( "U234" );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u234 && u235 );

  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint wide_u234, wide_u235;
  wide_u234.nuclide = u234;  wide_u234.lower_mass_fraction = 0.0;  wide_u234.upper_mass_fraction = 1.0;
  wide_u235.nuclide = u235;  wide_u235.lower_mass_fraction = 0.0;  wide_u235.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( wide_u234 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( wide_u235 );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Two-wide-overlap solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 50.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );

  const pair<double,optional<double>> mf234 = sol.mass_enrichment_fraction( u234, 0 );
  const pair<double,optional<double>> mf235 = sol.mass_enrichment_fraction( u235, 0 );
  BOOST_CHECK( (mf234.first >= -1.0e-6) && (mf234.first <= 1.0 + 1.0e-6) );
  BOOST_CHECK( (mf235.first >  0.0)     && (mf235.first <= 1.0 + 1.0e-6) );
  BOOST_CHECK_MESSAGE( (mf234.first + mf235.first) < 1.0,
                       "Sum of constrained U234+U235 fractions " << (mf234.first + mf235.first)
                       << " must be < 1." );
}


// A.9 - The constrained reparameterization is benign on a feasible config.  Adding a single feasible
//  U235 [0,1] constraint is the *same model in different coordinates* (the m=1 sigma-block decode is
//  the identity - see sigma_block_m1_identity_and_reachability), so it must reach a fit quality
//  comparable to an *unconstrained* solve of the same spectrum.  We deliberately do NOT assert
//  bit-identical enrichment: this spectrum's U235 fraction is only weakly determined (its cost
//  surface is shallow - the 185.7 keV peak alone sits ~10 sigma off the DRF-limited model), so the
//  two parametrizations settle on different near-degenerate minima (~25 vs ~28 wt%) at essentially
//  equal cost - a gap that long predates the constraint-decode changes.  The strict before/after
//  "feasible fits don't move" guard is the idb_enrichment_check over the HPGe corpus, not this
//  single shallow case.
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_feasible_invariance )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case( tc ), "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u235 );

  RelActCalcAuto::RelActAutoSolution sol_unc;
  BOOST_REQUIRE_NO_THROW( sol_unc = RelActCalcAuto::solve( tc.options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );
  BOOST_REQUIRE_MESSAGE( sol_unc.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                         "Unconstrained reference solve failed (status=" << static_cast<int>(sol_unc.m_status) << ")." );

  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint c;
  c.nuclide = u235;  c.lower_mass_fraction = 0.0;  c.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c );

  RelActCalcAuto::RelActAutoSolution sol_con;
  BOOST_REQUIRE_NO_THROW( sol_con = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );
  BOOST_REQUIRE_MESSAGE( sol_con.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                         "Single-constraint solve failed (status=" << static_cast<int>(sol_con.m_status) << ")." );

  // Same model, different coordinates: the soft-cap must not materially degrade the achievable fit.
  BOOST_REQUIRE( (sol_unc.m_dof > 0) && (sol_con.m_dof > 0) );
  const double chi2dof_unc = sol_unc.m_chi2 / sol_unc.m_dof;
  const double chi2dof_con = sol_con.m_chi2 / sol_con.m_dof;
  BOOST_CHECK_MESSAGE( chi2dof_con <= 2.0*chi2dof_unc + 1.0e-9,
                       "Constrained chi2/dof " << chi2dof_con << " is materially worse than unconstrained "
                       << chi2dof_unc << " - soft-cap reparametrization should not degrade a feasible fit." );

  // ... and the reported enrichment must be a feasible fraction.
  const double enr_con = sol_con.mass_enrichment_fraction( u235, 0 ).first;
  BOOST_CHECK_MESSAGE( (enr_con > 0.0) && (enr_con < 1.0),
                       "Constrained U235 enrichment " << enr_con << " is not a feasible fraction." );
}


// C.1 - ALL nuclides of the element mass-fraction constrained (previously a hard validation error):
//  the sigma-block frees the carrier slot to hold the element's total rel-mass scale, so this must
//  solve, with the element's fractions summing to exactly 1.
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_all_constrained_element )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case( tc ), "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u232 = db->nuclide( "U232" );
  const SandiaDecay::Nuclide * const u234 = db->nuclide( "U234" );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u232 && u234 && u235 && u238 );

  // Constrain every uranium isotope of the embedded curve (U232/U234/U235/U238): windows must be
  //  able to sum to exactly 1 (Sum lower = 0 <= 1 <= Sum upper).
  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint c232, c234, c235, c238;
  c232.nuclide = u232;  c232.lower_mass_fraction = 0.0;  c232.upper_mass_fraction = 1.0e-6;
  c234.nuclide = u234;  c234.lower_mass_fraction = 0.0;  c234.upper_mass_fraction = 0.05;
  c235.nuclide = u235;  c235.lower_mass_fraction = 0.0;  c235.upper_mass_fraction = 1.0;
  c238.nuclide = u238;  c238.lower_mass_fraction = 0.0;  c238.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c232 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c234 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c235 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c238 );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "All-constrained solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 50.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );

  // Every fraction in-window, and the element's fractions summing to exactly 1.
  const double f232 = sol.mass_enrichment_fraction( u232, 0 ).first;
  const double f234 = sol.mass_enrichment_fraction( u234, 0 ).first;
  const double f235 = sol.mass_enrichment_fraction( u235, 0 ).first;
  const double f238 = sol.mass_enrichment_fraction( u238, 0 ).first;

  BOOST_CHECK( (f232 >= -1.0e-9) && (f232 <= 1.0e-6 + 1.0e-9) );
  BOOST_CHECK( (f234 >= -1.0e-9) && (f234 <= 0.05 + 1.0e-6) );
  BOOST_CHECK( (f235 >  0.0)     && (f235 <= 1.0 + 1.0e-6) );
  BOOST_CHECK( (f238 >  0.0)     && (f238 <= 1.0 + 1.0e-6) );
  BOOST_CHECK_MESSAGE( fabs( (f232 + f234 + f235 + f238) - 1.0 ) < 1.0e-6,
                       "All-constrained U fractions sum to " << (f232 + f234 + f235 + f238)
                       << " instead of 1." );

  // This spectrum's enrichment is ~25-28 wt%; require the all-constrained fit lands in a sane band.
  BOOST_CHECK_MESSAGE( (f235 > 0.05) && (f235 < 0.6),
                       "All-constrained U235 fraction " << f235 << " is far from the expected ~0.25-0.28." );
}


// C.2 - A narrow high window that the data disagrees with: U235 constrained to [0.9, 1.0] when the
//  truth is ~0.26.  The hard sigma-block bounds must hold the reported fraction inside the window
//  (pinned at/near its lower edge) - hard constraints stay hard even against the data - and the
//  solve must still complete (no invalid-step wall).
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_narrow_high_window )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case( tc ), "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  BOOST_REQUIRE( u235 );

  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint c;
  c.nuclide = u235;  c.lower_mass_fraction = 0.9;  c.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Narrow-high-window solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  const double f235 = sol.mass_enrichment_fraction( u235, 0 ).first;
  BOOST_CHECK_MESSAGE( (f235 >= 0.9 - 1.0e-6) && (f235 <= 1.0 + 1.0e-6),
                       "U235 fraction " << f235 << " escaped its hard [0.9, 1.0] window." );
  // The data wants ~0.26, so the fit should be pressed against the lower edge of the window.
  BOOST_CHECK_MESSAGE( f235 < 0.95,
                       "U235 fraction " << f235 << " not pinned toward the windows lower edge - "
                       "unexpected for data whose enrichment is ~0.26." );
}


// C.3 - Fully specified isotopics: every uranium isotope FIXED (lower == upper, summing to 1); only
//  the element's total (the freed carrier slot) is fit.  Fractions must decode back pinned.
BOOST_AUTO_TEST_CASE( mass_fraction_constraint_all_fixed_element )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case( tc ), "Failed to load U235_Unshielded_6000.n42 test case." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u232 = db->nuclide( "U232" );
  const SandiaDecay::Nuclide * const u234 = db->nuclide( "U234" );
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u232 && u234 && u235 && u238 );

  // Pin the isotopics near this spectrum's expected composition; values must sum to exactly 1.
  const double f232 = 1.0e-8, f234 = 0.015, f235 = 0.26;
  const double f238 = 1.0 - f232 - f234 - f235;

  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint c232, c234, c235, c238;
  c232.nuclide = u232;  c232.lower_mass_fraction = c232.upper_mass_fraction = f232;
  c234.nuclide = u234;  c234.lower_mass_fraction = c234.upper_mass_fraction = f234;
  c235.nuclide = u235;  c235.lower_mass_fraction = c235.upper_mass_fraction = f235;
  c238.nuclide = u238;  c238.lower_mass_fraction = c238.upper_mass_fraction = f238;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c232 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c234 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c235 );
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c238 );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {},
                                    PeakFitUtils::coarse_det_type( tc.foreground, nullptr ), nullptr ) );

  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "All-fixed solve failed: status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "'" );
  if( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    return;

  // Fractions decode back to exactly their pinned values...
  BOOST_CHECK_SMALL( sol.mass_enrichment_fraction( u234, 0 ).first - f234, 1.0e-9 );
  BOOST_CHECK_SMALL( sol.mass_enrichment_fraction( u235, 0 ).first - f235, 1.0e-9 );
  BOOST_CHECK_SMALL( sol.mass_enrichment_fraction( u238, 0 ).first - f238, 1.0e-9 );

  // ... and the fitted element scale gives positive relative activities.
  BOOST_CHECK_MESSAGE( sol.m_dof > 0 && (sol.m_chi2 / sol.m_dof) < 100.0,
                       "chi2/dof = " << (sol.m_dof > 0 ? sol.m_chi2/sol.m_dof : -1.0) << " is unreasonable." );
}


// CONV-00/10: caller source order is presentation state, not optimization state.  Exercise both
// permutations of this deliberately small two-source mass-fraction problem; the five-repeat
// multi-curve permutation matrix lives in test_RelActCalcAuto_MultiCurve.cpp.
BOOST_AUTO_TEST_CASE( source_order_is_semantically_canonical )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u235 && u238 );
  const RelActCalcAuto::Options canonical_options
                           = small_canonical_u_options(tc.options,u235,u238);
  BOOST_REQUIRE_EQUAL( canonical_options.rel_eff_curves[0].nuclides.size(),2u );

  RelActCalcAuto::Options semantic_orders[2]{canonical_options,canonical_options};
  std::reverse( semantic_orders[1].rel_eff_curves[0].nuclides.begin(),
                semantic_orders[1].rel_eff_curves[0].nuclides.end() );

  const PeakFitUtils::CoarseResolutionType det_type
      = PeakFitUtils::coarse_det_type( tc.foreground, nullptr );
  std::uint64_t reference_gamma_hash = 0, reference_layout_hash = 0;
  double reference_objective = 0.0, reference_data_objective = 0.0;
  double reference_u235 = 0.0, reference_u238 = 0.0;
  bool exercised_orders[2]{false,false};

  for( size_t repeat = 0; repeat < 2; ++repeat )
  {
    const size_t order = repeat % 2;
    exercised_orders[order] = true;
    RelActCalcAuto::RelActAutoSolution solution;
    BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
        semantic_orders[order],tc.foreground,tc.background,tc.det,{},det_type,nullptr) );
    BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                           "Repeat " << repeat << " (source order " << order
                                     << ") failed: " << solution.m_error_message );
    BOOST_REQUIRE_NE( solution.m_frozen_gamma_membership_hash,UINT64_C(0) );
    BOOST_REQUIRE_NE( solution.m_frozen_layout_hash,UINT64_C(0) );
    BOOST_REQUIRE_EQUAL( solution.m_options.rel_eff_curves[0].nuclides.size(),
                         solution.m_rel_activities[0].size() );
    for( size_t row = 0; row < solution.m_rel_activities[0].size(); ++row )
      BOOST_CHECK( solution.m_rel_activities[0][row].source
                   == semantic_orders[order].rel_eff_curves[0].nuclides[row].source );

    const double fraction_u235 = solution.mass_enrichment_fraction(u235,0).first;
    const double fraction_u238 = solution.mass_enrichment_fraction(u238,0).first;
    BOOST_REQUIRE( std::isfinite(fraction_u235) && std::isfinite(fraction_u238) );
    if( repeat == 0 )
    {
      reference_gamma_hash = solution.m_frozen_gamma_membership_hash;
      reference_layout_hash = solution.m_frozen_layout_hash;
      reference_objective = solution.m_chi2;
      reference_data_objective = solution.m_chi2_data;
      reference_u235 = fraction_u235;
      reference_u238 = fraction_u238;
      continue;
    }

    BOOST_CHECK_EQUAL( solution.m_frozen_gamma_membership_hash,reference_gamma_hash );
    BOOST_CHECK_EQUAL( solution.m_frozen_layout_hash,reference_layout_hash );
    const double full_scale = (std::max)(1.0,std::fabs(reference_objective));
    const double data_scale = (std::max)(1.0,std::fabs(reference_data_objective));
    BOOST_CHECK_SMALL( solution.m_chi2-reference_objective,1.0e-8*full_scale );
    BOOST_CHECK_SMALL( solution.m_chi2_data-reference_data_objective,1.0e-8*data_scale );
    BOOST_CHECK_SMALL( fraction_u235-reference_u235,1.0e-7 );
    BOOST_CHECK_SMALL( fraction_u238-reference_u238,1.0e-7 );
  }
  BOOST_CHECK( exercised_orders[0] && exercised_orders[1] );
}


// AD-13/CONV-00: the production solver must freeze and select the same physical problem when its
// Ceres/ROI work is forced serial or allowed to use multiple threads.  Keep the nuisance families
// fixed so this specifically detects shared cache, summation-order, or semantic-layout drift.
BOOST_AUTO_TEST_CASE( production_thread_count_is_semantically_invariant )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );

  const RelActCalcAuto::Options options
      = small_canonical_u_options(tc.options,u235,u238);
  const PeakFitUtils::CoarseResolutionType det_type
      = PeakFitUtils::coarse_det_type(tc.foreground,nullptr);

  SolveThreadLimitRestorer restore_thread_limit;
  RelActCalc::set_max_solve_threads(0);
  const unsigned available_threads = RelActCalc::max_solve_threads();
  if( available_threads < 2 )
  {
    BOOST_TEST_MESSAGE( "Thread-count invariance skipped: only one hardware thread is available." );
    return;
  }

  RelActCalc::set_max_solve_threads(1);
  BOOST_REQUIRE_EQUAL( RelActCalc::max_solve_threads(),1u );
  RelActCalcAuto::RelActAutoSolution serial;
  BOOST_REQUIRE_NO_THROW( serial = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,{},det_type,nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(serial.m_status),
                         "Serial solve failed: " << serial.m_error_message );

  const unsigned parallel_limit = (std::min)(4u,available_threads);
  RelActCalc::set_max_solve_threads(parallel_limit);
  BOOST_REQUIRE_GT( RelActCalc::max_solve_threads(),1u );
  RelActCalcAuto::RelActAutoSolution parallel;
  BOOST_REQUIRE_NO_THROW( parallel = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,{},det_type,nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(parallel.m_status),
                         "Parallel solve failed: " << parallel.m_error_message );

  BOOST_REQUIRE_NE( serial.m_frozen_gamma_membership_hash,UINT64_C(0) );
  BOOST_REQUIRE_NE( serial.m_frozen_layout_hash,UINT64_C(0) );
  BOOST_CHECK_EQUAL( serial.m_frozen_gamma_membership_hash,
                     parallel.m_frozen_gamma_membership_hash );
  BOOST_CHECK_EQUAL( serial.m_frozen_layout_hash,parallel.m_frozen_layout_hash );
  BOOST_CHECK_EQUAL( serial.m_frozen_model_policy_hash,parallel.m_frozen_model_policy_hash );
  BOOST_CHECK_SMALL( serial.m_chi2-parallel.m_chi2,
                     1.0e-8*(std::max)(1.0,std::fabs(serial.m_chi2)) );
  BOOST_CHECK_SMALL( serial.m_chi2_data-parallel.m_chi2_data,
                     1.0e-8*(std::max)(1.0,std::fabs(serial.m_chi2_data)) );

  for( const SandiaDecay::Nuclide * const nuclide : {u235,u238} )
  {
    const double serial_fraction = serial.mass_enrichment_fraction(nuclide,0).first;
    const double parallel_fraction = parallel.mass_enrichment_fraction(nuclide,0).first;
    BOOST_REQUIRE( std::isfinite(serial_fraction) && std::isfinite(parallel_fraction) );
    BOOST_CHECK_SMALL( serial_fraction-parallel_fraction,1.0e-7 );
  }
  const double serial_u235 = serial.mass_enrichment_fraction(u235,0).first;
  const double parallel_u235 = parallel.mass_enrichment_fraction(u235,0).first;
  BOOST_CHECK( (serial_u235 >= 0.10) && (serial_u235 <= 0.50) );
  BOOST_CHECK( (parallel_u235 >= 0.10) && (parallel_u235 <= 0.50) );
}


// CONV-01/03: two curves with disjoint source identities have no EM/shared-line rescue candidate.
// They must nevertheless select a stable basin, report physical sole-isotope fractions, and avoid
// inventing cross-curve source-attribution rows.
BOOST_AUTO_TEST_CASE( no_shared_source_multicurve_is_deterministic )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );

  RelActCalcAuto::Options options = small_canonical_u_options(tc.options,u235,u238);
  BOOST_REQUIRE_EQUAL( options.rel_eff_curves.size(),1u );
  RelActCalcAuto::RelEffCurveInput u235_curve = options.rel_eff_curves.front();
  RelActCalcAuto::RelEffCurveInput u238_curve = options.rel_eff_curves.front();
  u235_curve.name = "U-235 only";
  u238_curve.name = "U-238 only";
  u235_curve.mass_fraction_constraints.clear();
  u238_curve.mass_fraction_constraints.clear();
  u235_curve.nuclides.erase( std::remove_if(u235_curve.nuclides.begin(),u235_curve.nuclides.end(),
      [u235]( const RelActCalcAuto::NucInputInfo &source ) {
        return RelActCalcAuto::nuclide(source.source) != u235;
      }),u235_curve.nuclides.end() );
  u238_curve.nuclides.erase( std::remove_if(u238_curve.nuclides.begin(),u238_curve.nuclides.end(),
      [u238]( const RelActCalcAuto::NucInputInfo &source ) {
        return RelActCalcAuto::nuclide(source.source) != u238;
      }),u238_curve.nuclides.end() );
  BOOST_REQUIRE_EQUAL( u235_curve.nuclides.size(),1u );
  BOOST_REQUIRE_EQUAL( u238_curve.nuclides.size(),1u );
  options.rel_eff_curves = {u235_curve,u238_curve};
  options.same_corr_fcn_for_all_rel_eff_curves = false;
  options.same_external_shielding_for_all_rel_eff_curves = false;

  const PeakFitUtils::CoarseResolutionType det_type
      = PeakFitUtils::coarse_det_type(tc.foreground,nullptr);
  RelActCalcAuto::RelActAutoSolution solutions[2];
  for( size_t repeat = 0; repeat < 2; ++repeat )
  {
    BOOST_REQUIRE_NO_THROW( solutions[repeat] = RelActCalcAuto::solve(
        options,tc.foreground,tc.background,tc.det,{},det_type,nullptr) );
    BOOST_REQUIRE_MESSAGE(
        RelActCalcAuto::RelActAutoSolution::is_usable_status(solutions[repeat].m_status),
        "No-shared-source repeat " << repeat << " failed: "
                                    << solutions[repeat].m_error_message );
    BOOST_REQUIRE_EQUAL( solutions[repeat].m_rel_activities.size(),2u );
    BOOST_REQUIRE_EQUAL( solutions[repeat].m_rel_activities[0].size(),1u );
    BOOST_REQUIRE_EQUAL( solutions[repeat].m_rel_activities[1].size(),1u );
    BOOST_CHECK_GT( solutions[repeat].m_rel_activities[0][0].rel_activity,0.0 );
    BOOST_CHECK_GT( solutions[repeat].m_rel_activities[1][0].rel_activity,0.0 );
    BOOST_CHECK_SMALL( solutions[repeat].mass_enrichment_fraction(u235,0).first-1.0,1.0e-12 );
    BOOST_CHECK_SMALL( solutions[repeat].mass_enrichment_fraction(u238,1).first-1.0,1.0e-12 );
    BOOST_CHECK( solutions[repeat].source_count_attributions().empty() );
    BOOST_CHECK( solutions[repeat].m_enrichment_diff_z.empty() );
  }

  BOOST_REQUIRE_NE( solutions[0].m_frozen_gamma_membership_hash,UINT64_C(0) );
  BOOST_REQUIRE_NE( solutions[0].m_frozen_layout_hash,UINT64_C(0) );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_gamma_membership_hash,
                     solutions[1].m_frozen_gamma_membership_hash );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_layout_hash,solutions[1].m_frozen_layout_hash );
  BOOST_CHECK_EQUAL( solutions[0].m_frozen_model_policy_hash,
                     solutions[1].m_frozen_model_policy_hash );
  BOOST_CHECK_SMALL( solutions[0].m_chi2-solutions[1].m_chi2,
                     1.0e-8*(std::max)(1.0,std::fabs(solutions[0].m_chi2)) );
  BOOST_CHECK_SMALL( solutions[0].m_chi2_data-solutions[1].m_chi2_data,
                     1.0e-8*(std::max)(1.0,std::fabs(solutions[0].m_chi2_data)) );
}


// CONV-01/03: a robust solve always runs the complete applicable candidate matrix (candidates are
// now warm in-frame re-solves ranked on the unsimplified model, before backward elimination).  In
// the small fixed-nuisance problem above exactly four named seeds are applicable: default, the two
// semantic activity splits, and the weakest-direction checkerboard.  The production diagnostic is
// deliberately asserted here so a successful early candidate cannot silently truncate the matrix.
BOOST_AUTO_TEST_CASE( forced_search_polishes_complete_applicable_candidate_matrix )
{
  set_data_dir();
  ULoadResult tc;
  BOOST_REQUIRE_MESSAGE( load_u_test_case(tc),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide("U235");
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );

  RelActCalcAuto::Options options = small_canonical_u_options(tc.options,u235,u238);
  options.auto_simplify_model = true; //candidates rank pre-elimination; elimination runs on the winner
  options.auto_simplify_max_dchi2 = 0.0;
  options.robust_solve = true;        //the multi-start candidate search is opt-in

  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,{},
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         solution.m_error_message );

  bool saw_forced_complete_matrix = false;
  for( const string &warning : solution.m_warnings )
  {
    if( SpecUtils::icontains(warning,"Deterministic candidate search was triggered")
        && SpecUtils::icontains(warning,"robust solve requested") )
    {
      BOOST_CHECK_MESSAGE( SpecUtils::icontains(warning,"polished 4 named candidates"),warning );
      saw_forced_complete_matrix = SpecUtils::icontains(warning,"polished 4 named candidates");
    }
  }
  BOOST_CHECK_MESSAGE( saw_forced_complete_matrix,
                       "No diagnostic confirmed all four applicable semantic candidates" );
}




/** What `check_u235_endpoint_oracle` measured. */
struct EndpointOracleOutcome
{
  size_t endpoints_checked = 0;
  double worst_shortfall = -std::numeric_limits<double>::infinity();
};

/** The refit oracle shared by the two endpoint gates below: for every likelihood-crossing endpoint
 of the baseline's U-235 profile, refit with `MassFractionConstraint(lower == upper)` at that
 endpoint and check the objective lands ON the threshold.

 On uranium the oracle is exact: with no Pu-242 correlation active the reported and uncorrected
 coordinates coincide, so constraining the mass fraction constrains precisely the coordinate the
 interval is quoted in.  (On plutonium it would NOT be - the correlation renormalizes the reported
 fraction away from anything a parameter constraint can reach - which is why the Pu arm of the
 harness secants the pin onto the reported coordinate instead; see
 `Pu610775::pinned_pu_endpoints_reproduce_the_threshold_on_a_fixed_fraction_refit`.)

 The gates are deliberately two-sided.  A refit BELOW the threshold (a positive shortfall) means
 the scan's conditional points sat above their true minima and the interval is too NARROW; a refit
 ABOVE it means the fitted crossing model was too flat and the interval too WIDE - the direction no
 uncertainty tool may fail in silently.  Note the oracle validates against the basin the scan
 itself explored: it measures local exactness, not global optimality.

 The refit sees exactly the objective the baseline minimized: the SPECTRUM peaks, not the fitted
 ones (substituting fitted peaks would change FWHM initialization and the nonlinear calibration
 anchors, corrupting the very delta-chi2 being measured).
 */
EndpointOracleOutcome check_u235_endpoint_oracle( const UProfileFixture &fixture,
                                                  const RelActCalcAuto::Options &options,
                                                  const RelActCalcAuto::RelActAutoSolution &baseline,
                                                  const char * const tag )
{
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const auto nominal = baseline.mass_enrichment_result( u235,0 );
  BOOST_REQUIRE_MESSAGE( nominal.profile.has_value(), "U-235 was not profiled" );
  BOOST_REQUIRE_EQUAL( nominal.profile->intervals.size(),2u );
  BOOST_REQUIRE( std::isfinite(baseline.m_chi2) );
  BOOST_REQUIRE_GE( baseline.m_cov_scale,0.0 );

  RelActCalcAuto::Options refit_options = options;
  refit_options.auto_profile_weak_mass_fractions = false;
  refit_options.auto_simplify_model = false;
  for( RelActCalcAuto::NucInputInfo &input : refit_options.rel_eff_curves.at(0).nuclides )
    input.force_profile_mass_fraction = false;

  const double cov_scale = (std::max)( 1.0,baseline.m_cov_scale );
  EndpointOracleOutcome outcome;

  for( const auto &interval : nominal.profile->intervals )
  {
    const double threshold = interval.delta_chi2;
    BOOST_REQUIRE( std::isfinite(threshold) && (threshold > 0.0) );

    const std::pair<double,RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind>
        endpoints[2] = { {interval.lower,interval.lower_kind},
                         {interval.upper,interval.upper_kind} };

    for( const auto &endpoint : endpoints )
    {
      // Only a likelihood crossing makes a claim about the objective.  A physical or input bound
      // says the scan ran out of feasible room, which is a different statement entirely.
      if( endpoint.second
          != RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind::LikelihoodCrossing )
        continue;

      const double requested = endpoint.first;
      BOOST_REQUIRE( std::isfinite(requested) && (requested > 0.0) && (requested < 1.0) );

      // REPLACE, never append: `check_nuclide_constraints()` throws on two constraints for one
      // nuclide, and the user-chart fixture already carries U-235's [0.10,0.50] window.
      RelActCalcAuto::Options trial = refit_options;
      RelActCalcAuto::RelEffCurveInput &curve = trial.rel_eff_curves.at(0);
      curve.mass_fraction_constraints.erase(
          std::remove_if( curve.mass_fraction_constraints.begin(),
                          curve.mass_fraction_constraints.end(),
              [u235]( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint ){
                return constraint.nuclide == u235;
              } ), curve.mass_fraction_constraints.end() );

      RelActCalcAuto::RelEffCurveInput::MassFractionConstraint pinned;
      pinned.nuclide = u235;
      pinned.lower_mass_fraction = requested;
      pinned.upper_mass_fraction = requested;
      curve.mass_fraction_constraints.push_back( pinned );

      RelActCalcAuto::RelActAutoSolution refit;
      BOOST_REQUIRE_NO_THROW( refit = RelActCalcAuto::solve(
          trial,fixture.tc.foreground,fixture.tc.background,baseline.m_drf,
          baseline.m_spectrum_peaks,
          PeakFitUtils::coarse_det_type(fixture.tc.foreground,nullptr),nullptr) );
      BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(refit.m_status),
                             refit.m_error_message );

      // The oracle must actually have reached the coordinate it was asked for, or it is measuring
      // something else entirely.
      const auto achieved = refit.mass_enrichment_result( u235,0 );
      BOOST_CHECK_MESSAGE( std::fabs(achieved.fraction - requested) <= 1.0e-5,
                           "oracle refit landed at " << achieved.fraction << " rather than "
                                                     << requested );

      const double delta = refit.m_chi2 - baseline.m_chi2;
      const double shortfall = threshold - delta;
      outcome.worst_shortfall = (std::max)( outcome.worst_shortfall,shortfall );
      ++outcome.endpoints_checked;

      BOOST_TEST_MESSAGE( tag << " q=" << requested << " threshold=" << threshold
                          << " refit_delta=" << delta << " shortfall=" << shortfall
                          << " cov_scale=" << baseline.m_cov_scale );

      // Both gates carry the same cov_scale factor, because the thresholds themselves are
      // cov_scale*{1,3.841} - an earlier draft paired a relative loose gate with an absolute hard
      // one, which differed by 17x on a real Pu fixture.
      const double loose = cov_scale*(std::max)( 0.25,0.15*threshold/cov_scale );
      BOOST_CHECK_MESSAGE( std::fabs(delta - threshold) <= loose,
          "refit at the predicted endpoint missed the threshold: delta=" << delta
          << " threshold=" << threshold << " tolerance=" << loose );

      BOOST_CHECK_MESSAGE( delta <= (threshold + 0.25*cov_scale),
          "the refit objective exceeded the threshold, so the profile model was too flat: delta="
          << delta << " threshold=" << threshold );
    }
  }

  BOOST_CHECK_MESSAGE( outcome.endpoints_checked > 0,
                       "no likelihood-crossing endpoint was available to validate" );

  // The loop above SKIPS every endpoint that is not a likelihood crossing, so an endpoint wrongly
  // demoted to a bound would be silently un-validated rather than failed - and the endpoint kind
  // is itself reported (as `lower_endpoint`/`upper_endpoint` in the report and LLM JSON), so a
  // misclassification asserts a feasibility bound that does not exist.  Pin both the count and
  // the classification: this U-235 fixture is well determined on both sides, so every endpoint of
  // both intervals must be a genuine crossing, and none may claim the scan's synthetic range cap
  // (which only an open activity/ratio box can reach - never the carrier chart).
  BOOST_CHECK_EQUAL( outcome.endpoints_checked, 4u );
  for( const auto &interval : baseline.mass_enrichment_result(fixture.u235,0).profile->intervals )
  {
    BOOST_CHECK( interval.lower_kind
                 != RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind::ScanRangeLimit );
    BOOST_CHECK( interval.upper_kind
                 != RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind::ScanRangeLimit );
  }

  return outcome;
}//check_u235_endpoint_oracle(...)


/** The endpoint gate, USER-chart arm: `u235_forced_profile_options`'s [0.10,0.50] window makes
 U-235 the sole range-constrained nuclide of its element, so the pinned carrier slot IS the
 reported fraction and the scan is exact by construction.  The shortfall must come out at
 essentially zero - which makes this case the harness's own self-test: if it does not, the harness
 is wrong and no measurement downstream can be trusted.

 Deviation from the plan, recorded deliberately: this lives beside the fixtures rather than in a
 separate `test_RelActCalcAuto_ProfilePinning.cpp`, because the uranium fixtures the harness needs
 and their frozen setup already live here, and duplicating them would create two things that must
 be kept identical.
 */
BOOST_AUTO_TEST_CASE( pinned_profile_endpoints_reproduce_the_threshold_on_refit )
{
  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );

  const RelActCalcAuto::Options options
      = u235_forced_profile_options( fixture.tc,fixture.u235,true );

  RelActCalcAuto::RelActAutoSolution baseline;
  BOOST_REQUIRE_NO_THROW( baseline = RelActCalcAuto::solve(
      options,fixture.tc.foreground,fixture.tc.background,fixture.tc.det,fixture.input_peaks,
      PeakFitUtils::coarse_det_type(fixture.tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(baseline.m_status),
                         baseline.m_error_message );

  const EndpointOracleOutcome outcome
      = check_u235_endpoint_oracle( fixture,options,baseline,"pinning-shortfall" );

  const double cov_scale = (std::max)( 1.0,baseline.m_cov_scale );
  BOOST_CHECK_MESSAGE( outcome.worst_shortfall <= 0.25*cov_scale,
      "the exactly-pinned fixture shows a shortfall of " << outcome.worst_shortfall
      << "; on an exact pin it must be ~0, so the harness is measuring the wrong thing" );
}


/** The endpoint gate, SYNTHETIC-chart arm: with NO user mass-fraction constraint, the profile
 machinery installs its own wide, non-binding [0,1] constraint per target so the pinned slot IS
 the reported coordinate, and the refit oracle must land on the threshold with ~zero shortfall -
 the same exactness the user-chart arm gets from its own window.  (Under the removed activity-pin
 engine this unconstrained fixture is precisely where the interval was measurably narrow -
 coverage tracking the squared pin/reported correlation - so a ~zero shortfall here is the whole point of the carrier
 reparameterization.)
 */
BOOST_AUTO_TEST_CASE( carrier_profile_endpoints_reproduce_the_threshold_on_refit )
{
  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );

  const RelActCalcAuto::Options options
      = u235_forced_profile_options( fixture.tc,fixture.u235,true,/*add_constraint=*/false );

  RelActCalcAuto::RelActAutoSolution baseline;
  BOOST_REQUIRE_NO_THROW( baseline = RelActCalcAuto::solve(
      options,fixture.tc.foreground,fixture.tc.background,fixture.tc.det,fixture.input_peaks,
      PeakFitUtils::coarse_det_type(fixture.tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(baseline.m_status),
                         baseline.m_error_message );

  const EndpointOracleOutcome outcome
      = check_u235_endpoint_oracle( fixture,options,baseline,"carrier-shortfall" );

  const double cov_scale = (std::max)( 1.0,baseline.m_cov_scale );
  BOOST_CHECK_MESSAGE( outcome.worst_shortfall <= 0.25*cov_scale,
      "the carrier engine shows a shortfall of " << outcome.worst_shortfall
      << "; the pinned slot IS the reported coordinate, so it must be ~0" );
}


/** The carrier reparameterization is installed and unwound per target on the retained problem, so
 the nominal answer must be bit-for-bit indistinguishable between a profiled and an unprofiled
 solve of the same problem - the public-surface form of the install's nominal-invariance check.
 */
BOOST_AUTO_TEST_CASE( carrier_reparam_leaves_nominal_unchanged )
{
  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  const auto run = [&]( const bool force_profile ){
    RelActCalcAuto::RelActAutoSolution solution;
    BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
        u235_forced_profile_options(tc,u235,force_profile,/*add_constraint=*/false),
        tc.foreground,tc.background,tc.det,input_peaks,
        PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
    BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                           solution.m_error_message );
    return solution;
  };

  const RelActCalcAuto::RelActAutoSolution without_profile = run( false );
  const RelActCalcAuto::RelActAutoSolution with_profile = run( true );

  const auto profiled = with_profile.mass_enrichment_result(u235,0);
  BOOST_REQUIRE( profiled.profile.has_value() );

  BOOST_REQUIRE( std::isfinite(without_profile.m_chi2) && std::isfinite(with_profile.m_chi2) );
  const double chi2_scale = (std::max)( 1.0,std::fabs(without_profile.m_chi2) );
  BOOST_CHECK_MESSAGE(
      std::fabs(with_profile.m_chi2 - without_profile.m_chi2) <= 1.0e-6*chi2_scale,
      "The carrier profile moved the objective: " << without_profile.m_chi2
        << " -> " << with_profile.m_chi2 );

  const double fraction_without = without_profile.mass_enrichment_result(u235,0).fraction;
  BOOST_REQUIRE( std::isfinite(fraction_without) && std::isfinite(profiled.fraction) );
  BOOST_CHECK_MESSAGE( std::fabs(profiled.fraction - fraction_without) <= 1.0e-9,
      "The carrier profile moved the reported nominal: " << fraction_without
        << " -> " << profiled.fraction );

  BOOST_REQUIRE_EQUAL( without_profile.m_final_parameters.size(),
                       with_profile.m_final_parameters.size() );
  for( size_t i = 0; i < without_profile.m_final_parameters.size(); ++i )
    BOOST_CHECK_MESSAGE(
        std::fabs(without_profile.m_final_parameters[i] - with_profile.m_final_parameters[i])
          <= 1.0e-12*(1.0 + std::fabs(without_profile.m_final_parameters[i])),
        "parameter " << i << " differs: " << without_profile.m_final_parameters[i]
          << " vs " << with_profile.m_final_parameters[i] );
}


/** A target the carrier route must refuse - activity bounds live on the very slot the install
 would re-purpose - reports a structured Failed carrying the reason.  There is deliberately no
 inexact engine to fall back to: an activity scan of this fixture was measured at
 a squared pin/reported correlation of ~0.004, i.e. an interval narrow by ~16x quoted with full confidence.
 */
BOOST_AUTO_TEST_CASE( carrier_refuses_bounded_slot_with_structured_failure )
{
  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  RelActCalcAuto::Options options
      = u235_forced_profile_options( tc,u235,true,/*add_constraint=*/false );
  for( RelActCalcAuto::NucInputInfo &input : options.rel_eff_curves.at(0).nuclides )
  {
    if( RelActCalcAuto::nuclide(input.source) == u235 )
    {
      // A 0.0 lower limit produces bit-identical Ceres bounds to the default, so the FIT is
      // unchanged; merely being present is what forces the carrier route to refuse.
      input.min_rel_act = 0.0;
    }
  }

  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         solution.m_error_message );

  const auto result = solution.mass_enrichment_result(u235,0);
  BOOST_REQUIRE_MESSAGE( result.profile.has_value(),
                         "the explicit request must still produce a structured result" );
  BOOST_CHECK_MESSAGE( result.profile->status
      == RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus::Failed,
      "expected a structured refusal, got status "
        << static_cast<int>(result.profile->status) << ": " << result.profile->message );
  BOOST_CHECK_MESSAGE( !result.profile->message.empty(),
                       "a refusal must say why the fraction cannot be scanned exactly" );
  BOOST_CHECK( result.profile->intervals.empty() );
}


/** The refit oracle for the `RelativeActivity` executor: fixing the activity at each claimed
 crossing endpoint must land the objective on the threshold - the same exactness discipline as
 the mass-fraction carrier (`check_u235_endpoint_oracle`), through the identity chart (the
 activity slot IS the reported coordinate, so there is no reparameterization to distrust). */
BOOST_AUTO_TEST_CASE( relative_activity_profile_endpoints_reproduce_the_threshold_on_refit )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;
  using ProfileStatus = RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus;
  using EndpointKind = RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind;

  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  RelActCalcAuto::Options options
      = u235_forced_profile_options( tc, u235, false, /*add_constraint=*/false );
  Target target;
  target.kind = Target::Kind::RelativeActivity;
  target.source = RelActCalcAuto::SrcVariant(u235);
  target.rel_eff_curve_index = 0;
  options.profile_targets = { target };

  RelActCalcAuto::RelActAutoSolution baseline;
  BOOST_REQUIRE_NO_THROW( baseline = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(baseline.m_status),
                         baseline.m_error_message );

  const RelActCalcAuto::RelActAutoSolution::ProfileResultEntry * const entry
      = baseline.profile_result( target );
  BOOST_REQUIRE_MESSAGE( entry, "the explicit RelativeActivity request produced no result" );
  BOOST_REQUIRE_MESSAGE( (entry->profile.status == ProfileStatus::Complete)
                          || (entry->profile.status == ProfileStatus::BoundaryLimited),
                         "unexpected profile status " << static_cast<int>(entry->profile.status)
                          << ": " << entry->profile.message );
  BOOST_REQUIRE_EQUAL( entry->profile.intervals.size(), 2u );
  BOOST_REQUIRE( std::isfinite(baseline.m_chi2) );
  const double cov_scale = (std::max)( 1.0, baseline.m_cov_scale );

  size_t endpoints_checked = 0;
  double worst_shortfall = -std::numeric_limits<double>::infinity();
  for( const RelActCalcAuto::RelActAutoSolution::MassFractionProfileInterval &interval
            : entry->profile.intervals )
  {
    const double threshold = interval.delta_chi2;
    const std::pair<double,EndpointKind> sides[2] = {
      { interval.lower, interval.lower_kind }, { interval.upper, interval.upper_kind } };
    for( const std::pair<double,EndpointKind> &side : sides )
    {
      // A physical or input bound says the scan ran out of feasible room - a different statement.
      if( side.second != EndpointKind::LikelihoodCrossing )
        continue;
      const double requested = side.first;
      if( !std::isfinite(requested) || (requested <= 0.0) )
        continue;

      RelActCalcAuto::Options trial
          = u235_forced_profile_options( tc, u235, false, /*add_constraint=*/false );
      for( RelActCalcAuto::NucInputInfo &input : trial.rel_eff_curves.at(0).nuclides )
      {
        if( RelActCalcAuto::nuclide(input.source) == u235 )
        {
          input.min_rel_act = requested;
          input.max_rel_act = requested;
        }
      }

      RelActCalcAuto::RelActAutoSolution refit;
      BOOST_REQUIRE_NO_THROW( refit = RelActCalcAuto::solve(
          trial,tc.foreground,tc.background,baseline.m_drf,baseline.m_spectrum_peaks,
          PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
      BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(refit.m_status),
                             refit.m_error_message );

      const double achieved = refit.rel_activity( target.source, 0 );
      BOOST_REQUIRE_MESSAGE(
          std::fabs(achieved - requested) <= 1.0e-5*(std::max)(1.0,std::fabs(requested)),
          "oracle refit landed at " << achieved << " rather than " << requested );

      const double delta = refit.m_chi2 - baseline.m_chi2;
      const double shortfall = threshold - delta;
      worst_shortfall = (std::max)( worst_shortfall, shortfall );
      BOOST_TEST_MESSAGE( "relact-shortfall q=" << requested << " threshold=" << threshold
                          << " refit_delta=" << delta << " shortfall=" << shortfall
                          << " cov_scale=" << cov_scale );

      const double loose = cov_scale*(std::max)( 0.25, 0.15*threshold/cov_scale );
      BOOST_CHECK_MESSAGE( std::fabs(delta - threshold) <= loose,
          "endpoint " << requested << " missed the threshold: delta=" << delta
           << " vs " << threshold );
      BOOST_CHECK_MESSAGE( delta <= (threshold + 0.25*cov_scale),
          "endpoint " << requested << " refit sits above threshold: interval too WIDE" );
      ++endpoints_checked;
    }//for( each side )
  }//for( each interval )

  BOOST_CHECK_MESSAGE( endpoints_checked > 0,
                       "no likelihood-crossing endpoints existed to check" );
  BOOST_CHECK_MESSAGE( worst_shortfall <= 0.25*cov_scale,
                       "worst shortfall " << worst_shortfall << " vs " << 0.25*cov_scale );
}


/** The refit oracle for the `ActivityRatio` executor: refitting with a fixed
 `ActRatioConstraint` at each claimed crossing endpoint must land the objective on the threshold.
 This is the acceptance gate for the slot-driven ratio reparameterization
 (`install_ratio_reparam`): a pinned bare numerator would under-cover because the denominator
 re-optimizes, and this oracle would catch exactly that. */
BOOST_AUTO_TEST_CASE( activity_ratio_profile_endpoints_reproduce_the_threshold_on_refit )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;
  using ProfileStatus = RelActCalcAuto::RelActAutoSolution::MassFractionProfileStatus;
  using EndpointKind = RelActCalcAuto::RelActAutoSolution::MassFractionProfileEndpointKind;

  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const SandiaDecay::Nuclide * const u238 = DecayDataBaseServer::database()->nuclide("U238");
  BOOST_REQUIRE( u235 && u238 );
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  RelActCalcAuto::Options options
      = u235_forced_profile_options( tc, u235, false, /*add_constraint=*/false );
  Target target;
  target.kind = Target::Kind::ActivityRatio;
  target.source = RelActCalcAuto::SrcVariant(u235);
  target.denominator = RelActCalcAuto::SrcVariant(u238);
  target.rel_eff_curve_index = 0;
  target.denominator_curve_index = 0;
  options.profile_targets = { target };

  RelActCalcAuto::RelActAutoSolution baseline;
  BOOST_REQUIRE_NO_THROW( baseline = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(baseline.m_status),
                         baseline.m_error_message );

  const RelActCalcAuto::RelActAutoSolution::ProfileResultEntry * const entry
      = baseline.profile_result( target );
  BOOST_REQUIRE_MESSAGE( entry, "the explicit ActivityRatio request produced no result" );
  BOOST_REQUIRE_MESSAGE( (entry->profile.status == ProfileStatus::Complete)
                          || (entry->profile.status == ProfileStatus::BoundaryLimited),
                         "unexpected profile status " << static_cast<int>(entry->profile.status)
                          << ": " << entry->profile.message );
  BOOST_REQUIRE_EQUAL( entry->profile.intervals.size(), 2u );
  const double cov_scale = (std::max)( 1.0, baseline.m_cov_scale );

  size_t endpoints_checked = 0;
  double worst_shortfall = -std::numeric_limits<double>::infinity();
  for( const RelActCalcAuto::RelActAutoSolution::MassFractionProfileInterval &interval
            : entry->profile.intervals )
  {
    const double threshold = interval.delta_chi2;
    const std::pair<double,EndpointKind> sides[2] = {
      { interval.lower, interval.lower_kind }, { interval.upper, interval.upper_kind } };
    for( const std::pair<double,EndpointKind> &side : sides )
    {
      if( side.second != EndpointKind::LikelihoodCrossing )
        continue;
      const double requested = side.first;
      if( !std::isfinite(requested) || (requested <= 0.0) )
        continue;

      RelActCalcAuto::Options trial
          = u235_forced_profile_options( tc, u235, false, /*add_constraint=*/false );
      RelActCalcAuto::RelEffCurveInput::ActRatioConstraint pin;
      pin.controlling_source = RelActCalcAuto::SrcVariant(u238);
      pin.constrained_source = RelActCalcAuto::SrcVariant(u235);
      pin.constrained_to_controlled_activity_ratio = requested;
      trial.rel_eff_curves.at(0).act_ratio_constraints.push_back( pin );

      RelActCalcAuto::RelActAutoSolution refit;
      BOOST_REQUIRE_NO_THROW( refit = RelActCalcAuto::solve(
          trial,tc.foreground,tc.background,baseline.m_drf,baseline.m_spectrum_peaks,
          PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
      BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(refit.m_status),
                             refit.m_error_message );

      const double achieved_num = refit.rel_activity( target.source, 0 );
      const double achieved_den = refit.rel_activity( target.denominator, 0 );
      BOOST_REQUIRE( std::isfinite(achieved_num) && std::isfinite(achieved_den)
                     && (achieved_den > 0.0) );
      const double achieved = achieved_num/achieved_den;
      BOOST_REQUIRE_MESSAGE(
          std::fabs(achieved - requested) <= 1.0e-5*(std::max)(1.0,std::fabs(requested)),
          "oracle refit landed at ratio " << achieved << " rather than " << requested );

      const double delta = refit.m_chi2 - baseline.m_chi2;
      const double shortfall = threshold - delta;
      worst_shortfall = (std::max)( worst_shortfall, shortfall );
      BOOST_TEST_MESSAGE( "ratio-shortfall q=" << requested << " threshold=" << threshold
                          << " refit_delta=" << delta << " shortfall=" << shortfall
                          << " cov_scale=" << cov_scale );

      const double loose = cov_scale*(std::max)( 0.25, 0.15*threshold/cov_scale );
      BOOST_CHECK_MESSAGE( std::fabs(delta - threshold) <= loose,
          "endpoint " << requested << " missed the threshold: delta=" << delta
           << " vs " << threshold );
      BOOST_CHECK_MESSAGE( delta <= (threshold + 0.25*cov_scale),
          "endpoint " << requested << " refit sits above threshold: interval too WIDE" );
      ++endpoints_checked;
    }//for( each side )
  }//for( each interval )

  BOOST_CHECK_MESSAGE( endpoints_checked > 0,
                       "no likelihood-crossing endpoints existed to check" );
  BOOST_CHECK_MESSAGE( worst_shortfall <= 0.25*cov_scale,
                       "worst shortfall " << worst_shortfall << " vs " << 0.25*cov_scale );
}


/** Schema smoke for the `Age` identity-chart executor: an explicit age target on a fitted age
 must produce a structured result (intervals or a reason), whatever the age's identifiability in
 this small frozen fixture. */
BOOST_AUTO_TEST_CASE( age_profile_target_produces_structured_result )
{
  using Target = RelActCalcAuto::Options::ProfileTarget;

  UProfileFixture fixture;
  BOOST_REQUIRE_MESSAGE( load_u_profile_fixture(fixture),
                         "Failed to load U235_Unshielded_6000.n42 test case." );
  const ULoadResult &tc = fixture.tc;
  const SandiaDecay::Nuclide * const u235 = fixture.u235;
  const vector<shared_ptr<const PeakDef>> &input_peaks = fixture.input_peaks;

  RelActCalcAuto::Options options
      = u235_forced_profile_options( tc, u235, false, /*add_constraint=*/false );
  options.rel_eff_curves.at(0).nucs_of_el_same_age = false;
  for( RelActCalcAuto::NucInputInfo &input : options.rel_eff_curves.at(0).nuclides )
  {
    if( RelActCalcAuto::nuclide(input.source) == u235 )
    {
      input.age = 20.0*PhysicalUnits::year;
      input.fit_age = true;
    }
  }
  Target target;
  target.kind = Target::Kind::Age;
  target.source = RelActCalcAuto::SrcVariant(u235);
  target.rel_eff_curve_index = 0;
  options.profile_targets = { target };
  BOOST_REQUIRE_EQUAL( target.why_not_usable(options), string() );

  RelActCalcAuto::RelActAutoSolution solution;
  BOOST_REQUIRE_NO_THROW( solution = RelActCalcAuto::solve(
      options,tc.foreground,tc.background,tc.det,input_peaks,
      PeakFitUtils::coarse_det_type(tc.foreground,nullptr),nullptr) );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(solution.m_status),
                         solution.m_error_message );

  const RelActCalcAuto::RelActAutoSolution::ProfileResultEntry * const entry
      = solution.profile_result( target );
  BOOST_REQUIRE_MESSAGE( entry, "the explicit Age request produced no result" );
  BOOST_CHECK_MESSAGE( !entry->profile.intervals.empty() || !entry->profile.message.empty(),
                       "age profile carried neither intervals nor a reason" );
}
