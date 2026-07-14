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

// Pulls in RelActCalc::mass_fraction_softcap_factor for the pure-math tests below; included after
//  <ceres/jet.h> and the std headers (per the header's own requirement) so it can be instantiated as
//  both double and ceres::Jet.
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
// A16 soft-cap (RelActCalc::mass_fraction_softcap_factor) - pure-math unit tests (fast; no spectra/solve)
// =====================================================================================================
namespace
{
  // Reconstructs the per-element mass-fraction decode the cost functor performs from the soft-cap
  //  factor, given each constrained nuclide's (lower, upper, g=rel_dist in [0,1]).
  struct SoftcapEval
  {
    double S = 0.0, B = 0.0, factor = 0.0, sum_constrained = 0.0;
    vector<double> v, f, f_des;
  };

  SoftcapEval eval_softcap( const vector<double> &lower, const vector<double> &upper,
                            const vector<double> &g,
                            const double eps = RelActCalc::ns_mass_frac_softcap_eps )
  {
    SoftcapEval r;
    const size_t n = lower.size();
    r.v.resize( n ); r.f.resize( n ); r.f_des.resize( n );
    double L = 0.0;
    for( size_t i = 0; i < n; ++i )
    {
      r.v[i] = g[i]*(upper[i] - lower[i]);
      r.f_des[i] = lower[i] + r.v[i];
      L += lower[i];
      r.S += r.v[i];
    }
    r.B = 1.0 - L;
    r.factor = RelActCalc::mass_fraction_softcap_factor( r.S, r.B, eps );
    for( size_t i = 0; i < n; ++i )
      r.f[i] = lower[i] + r.factor*r.v[i];
    r.sum_constrained = L + r.factor*r.S;
    return r;
  }
}//namespace


// A.1 - Feasible configs (variable demand comfortably below the knee) pass through unchanged: factor is
//  exactly 1, each decoded fraction equals the desired lerp, and stays inside its window.
BOOST_AUTO_TEST_CASE( softcap_feasible_is_identity )
{
  const vector<vector<double>> lowers = { {0.0, 0.0}, {0.1, 0.2}, {0.0, 0.05, 0.0} };
  const vector<vector<double>> uppers = { {0.5, 0.5}, {0.4, 0.5}, {0.3, 0.10, 0.2} };
  const vector<vector<double>> gs     = { {0.2, 0.3}, {0.5, 0.5}, {0.4, 0.50, 0.3} };

  for( size_t c = 0; c < lowers.size(); ++c )
  {
    const SoftcapEval r = eval_softcap( lowers[c], uppers[c], gs[c] );
    const double knee = RelActCalc::ns_mass_frac_softcap_knee*(1.0 - RelActCalc::ns_mass_frac_softcap_eps)*r.B;
    BOOST_REQUIRE_MESSAGE( r.S <= knee, "test case " << c << " is not in the identity region (S=" << r.S << ")" );

    BOOST_CHECK_SMALL( r.factor - 1.0, 1.0e-12 );
    for( size_t i = 0; i < r.f.size(); ++i )
    {
      BOOST_CHECK_SMALL( r.f[i] - r.f_des[i], 1.0e-12 );
      BOOST_CHECK( (r.f[i] >= lowers[c][i] - 1.0e-12) && (r.f[i] <= uppers[c][i] + 1.0e-12) );
    }
  }
}


// A.2 - Infeasible configs (Σ desired > 1) are capped: Σ actual < 1 strictly, every fraction stays in
//  window and is <= its desired value (monotone shrink), and the variable parts stay proportional.
BOOST_AUTO_TEST_CASE( softcap_infeasible_is_capped_and_proportional )
{
  const vector<double> lower = {0.0, 0.0}, upper = {1.0, 1.0}, g = {0.8, 0.7};
  const SoftcapEval r = eval_softcap( lower, upper, g );

  BOOST_CHECK( r.S > 1.0 );                       // genuinely infeasible at factor == 1
  BOOST_CHECK( r.factor < 1.0 );
  BOOST_CHECK( r.sum_constrained < 1.0 );         // structurally feasible
  for( size_t i = 0; i < r.f.size(); ++i )
  {
    BOOST_CHECK( (r.f[i] >= lower[i] - 1.0e-12) && (r.f[i] <= upper[i] + 1.0e-12) );
    BOOST_CHECK( r.f[i] <= r.f_des[i] + 1.0e-12 );      // monotone shrink
  }
  // proportional shrink: (f0-lower0)/(f1-lower1) == v0/v1  <=>  (f0-lower0)*v1 == (f1-lower1)*v0
  BOOST_CHECK_CLOSE( (r.f[0]-lower[0])*r.v[1], (r.f[1]-lower[1])*r.v[0], 1.0e-6 );
}


// A.3 - The decoded fraction never falls below the validated lower bound, for any demand; and the
//  variable sum is always held below (1-eps)*B (so the unconstrained remainder stays positive).
BOOST_AUTO_TEST_CASE( softcap_respects_lower_floor )
{
  const double eps = RelActCalc::ns_mass_frac_softcap_eps;
  for( const double B : { 1.0, 0.7, 0.2 } )
  {
    for( const double S : { 0.0, 0.1, 1.0, 5.0, 100.0, 1.0e6 } )
    {
      const double factor = RelActCalc::mass_fraction_softcap_factor( S, B, eps );
      BOOST_CHECK( (factor > 0.0) && (factor <= 1.0 + 1.0e-12) );
      BOOST_CHECK( (S*factor) <= (1.0 - eps)*B + 1.0e-9 );          // cap < (1-eps)*B
      const double lower_i = 0.123, v_i = 0.4;                       // representative nuclide
      BOOST_CHECK( (lower_i + factor*v_i) >= lower_i - 1.0e-12 );    // f_i >= lower_i
    }
  }
}


// A.4 - Sweeping a knob finely across the feasibility transition gives a continuous decode (no cliff),
//  and Σ never exceeds the (1 - eps*B) cap.  Direct regression against the old throw wall.
BOOST_AUTO_TEST_CASE( softcap_is_continuous_across_transition )
{
  const vector<double> lower = {0.0, 0.0}, upper = {1.0, 1.0};
  const int nsteps = 2000;
  double prev_factor = 0.0, prev_sum = 0.0;
  for( int k = 0; k <= nsteps; ++k )
  {
    const double g0 = static_cast<double>(k) / nsteps;             // 0 -> 1 ; S = g0 + 0.6 crosses the knee
    const SoftcapEval r = eval_softcap( lower, upper, { g0, 0.6 } );
    BOOST_CHECK( r.sum_constrained <= (1.0 - RelActCalc::ns_mass_frac_softcap_eps)*r.B + 1.0e-12 );
    if( k > 0 )
    {
      BOOST_CHECK_SMALL( r.factor - prev_factor, 5.0e-3 );
      BOOST_CHECK_SMALL( r.sum_constrained - prev_sum, 5.0e-3 );
    }
    prev_factor = r.factor;
    prev_sum = r.sum_constrained;
  }
}


// A.5 - Automatic differentiation (ceres::Jet) through the soft-cap matches central finite differences
//  in the feasible region, deep in the capped region, and at the knee.  Validates AD-safety (Option A).
BOOST_AUTO_TEST_CASE( softcap_ad_matches_finite_difference )
{
  typedef ceres::Jet<double,1> Jet;
  const double eps = RelActCalc::ns_mass_frac_softcap_eps;
  const double B = 1.0;
  const double knee = RelActCalc::ns_mass_frac_softcap_knee*(1.0 - eps)*B;

  const vector<double> test_S = { 0.0, 0.3, knee, 1.0, 1.5, 3.0 };
  for( const double S : test_S )
  {
    const Jet S_jet( S, 0 );                                        // value S, derivative seed at index 0
    const Jet f_jet = RelActCalc::mass_fraction_softcap_factor( S_jet, Jet(B), eps );
    const double ad = f_jet.v[0];
    BOOST_CHECK( std::isfinite( f_jet.a ) && std::isfinite( ad ) );

    const double h = 1.0e-6;
    const double fp = RelActCalc::mass_fraction_softcap_factor( S + h, B, eps );
    const double fm = RelActCalc::mass_fraction_softcap_factor( S - h, B, eps );
    const double fd = (fp - fm) / (2.0*h);
    BOOST_CHECK_MESSAGE( fabs(ad - fd) < 1.0e-5,
                         "AD/FD mismatch at S=" << S << ": ad=" << ad << " fd=" << fd );
  }
}


// A.6 - With all knobs at their lower bound (S == 0) the factor and its Jet derivative are finite (no 0/0).
BOOST_AUTO_TEST_CASE( softcap_s_zero_is_finite )
{
  typedef ceres::Jet<double,1> Jet;
  const double eps = RelActCalc::ns_mass_frac_softcap_eps;
  for( const double B : { 1.0, 0.5, 0.05 } )
  {
    const double factor = RelActCalc::mass_fraction_softcap_factor( 0.0, B, eps );
    BOOST_CHECK_SMALL( factor - 1.0, 1.0e-12 );
    const Jet f_jet = RelActCalc::mass_fraction_softcap_factor( Jet(0.0, 0), Jet(B), eps );
    BOOST_CHECK( std::isfinite( f_jet.a ) && std::isfinite( f_jet.v[0] ) );
    BOOST_CHECK_SMALL( f_jet.v[0], 1.0e-12 );                       // flat region -> zero slope
  }
}


// A.7 - A fixed (lower==upper) constraint contributes only a constant to L and is never scaled, while a
//  range nuclide on the same element stays feasible.
BOOST_AUTO_TEST_CASE( softcap_fixed_plus_range )
{
  const vector<double> lower = {0.3, 0.0}, upper = {0.3, 0.5}, g = {0.5, 0.4};   // nuc 0 fixed at 0.3
  const SoftcapEval r = eval_softcap( lower, upper, g );
  BOOST_CHECK_SMALL( r.v[0], 1.0e-15 );                            // fixed nuclide: zero variable amount
  BOOST_CHECK_SMALL( r.f[0] - 0.3, 1.0e-12 );                       // ... decodes to exactly its fixed fraction
  BOOST_CHECK( (r.f[1] >= 0.0) && (r.f[1] <= 0.5 + 1.0e-12) );
  BOOST_CHECK( r.sum_constrained < 1.0 );
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


// A.8 - A16 regression: an element carrying TWO wide overlapping constraints (U234 [0,1] and U235 [0,1]).
//  The old eval threw "Sum of constrained mass fractions ... > 1.0" whenever a trial step drove the two
//  fractions past the simplex face, which Ceres treats as a hard wall (invalid step) and could abort on.
//  The soft-cap makes Σ < 1 structural, so the solve must simply succeed with both fractions feasible.
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


// A.9 - The soft-cap reparametrization is benign on a feasible config.  Adding a single feasible
//  U235 [0,1] constraint is the *same model in different coordinates* (its decode is factor==1 here -
//  see softcap_feasible_is_identity), so it must reach a fit quality comparable to an *unconstrained*
//  solve of the same spectrum.  We deliberately do NOT assert bit-identical enrichment: this spectrum's
//  U235 fraction is only weakly determined (its cost surface is shallow - the 185.7 keV peak alone sits
//  ~10 sigma off the DRF-limited model), so the two parametrizations settle on different near-degenerate
//  minima (~25 vs ~28 wt%) at essentially equal cost - and that gap predates A16 (factor==1 makes the
//  decode identical to the old lerp).  The strict before/after "feasible fits don't move" guard is the
//  idb_enrichment_check over the HPGe corpus (verification step 4), not this single shallow case.
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
