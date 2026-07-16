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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, background, det, {}, nullptr ) );

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
    const double true_max = std::max( x, a );

    BOOST_CHECK_MESSAGE( qx >= true_max - 1.0e-15, "qmax below true max at x=" << x );
    BOOST_CHECK_MESSAGE( qx <= true_max + 0.25*r + 1.0e-15, "qmax excess above r/4 at x=" << x );
    if( fabs(x - a) >= r )
      BOOST_CHECK_SMALL( qx - true_max, 1.0e-15 );  //exact outside the blend zone

    const double qn = RelActCalc::qmin_hinge( x, a, r );
    const double true_min = std::min( x, a );
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
      windows[k] = { l, std::min(u, 1.0) };
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
}//namespace


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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {}, nullptr ) );

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
  BOOST_REQUIRE_NO_THROW( sol_unc = RelActCalcAuto::solve( tc.options, tc.foreground, tc.background, tc.det, {}, nullptr ) );
  BOOST_REQUIRE_MESSAGE( sol_unc.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                         "Unconstrained reference solve failed (status=" << static_cast<int>(sol_unc.m_status) << ")." );

  RelActCalcAuto::Options options = tc.options;
  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint c;
  c.nuclide = u235;  c.lower_mass_fraction = 0.0;  c.upper_mass_fraction = 1.0;
  options.rel_eff_curves[0].mass_fraction_constraints.push_back( c );

  RelActCalcAuto::RelActAutoSolution sol_con;
  BOOST_REQUIRE_NO_THROW( sol_con = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {}, nullptr ) );
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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {}, nullptr ) );

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {}, nullptr ) );

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
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, tc.foreground, tc.background, tc.det, {}, nullptr ) );

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
