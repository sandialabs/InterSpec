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
#include <string>
#include <memory>
#include <iostream>

#define BOOST_TEST_MODULE RelActAutoReport_suite
#include <boost/test/included/unit_test.hpp>

// `<windows.h>` (pulled in transitively on MSVC) defines `min`/`max` as macros
// that break `std::min`/`std::max` calls inside the vendored inja headers.
#if ( defined( WIN32 ) )
#undef min
#undef max
#endif

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

static string g_test_file_dir;

// Set the static data directory and locate the test_data directory; mirrors the helper in
// test_RelEffAuto.cpp.
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
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
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
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase appears empty" );

  const string rel_act_test_dir = SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( rel_act_test_dir ),
                         "RelActAutoReport test data dir not at '" << rel_act_test_dir << "'" );
}//void set_data_dir()


/** `RelActAutoSolution::print_html_report()` pulls its HTML template, plus the spectrum and
 rel-eff chart JS/CSS, out of `InterSpec_resources/` - relative to the Wt docroot, which with no
 `WApplication` is the CWD.  Non-Windows builds get a symlink to that directory in the build dir
 from CMake, but Windows builds do not (and the deployed `d3.v3.min.js`/`SpectrumChartD3.*` are
 only ever written into the binary dir), so there the HTML report simply cannot be generated.
 Returns whether every file `print_html_report()` needs is reachable, so that surface can be
 skipped rather than fail on an unrelated missing file.
 */
static bool html_report_resources_available()
{
  const char * const needed[] = {
    "static_text/auto_rel_act_report.tmplt.html",
    "d3.v3.min.js", "SpectrumChartD3.js", "SpectrumChartD3.css",
    "RelEffPlot.js", "RelEffPlot.css"
  };

  for( const char * const filename : needed )
  {
    if( !SpecUtils::is_file( SpecUtils::append_path( "InterSpec_resources", filename ) ) )
      return false;
  }

  return true;
}//bool html_report_resources_available()


// Returns true if any nuclide in any rel-eff curve of `sol` has a name containing `needle`.
static bool solution_has_nuclide_named( const RelActCalcAuto::RelActAutoSolution &sol,
                                        const string &needle )
{
  for( const vector<RelActCalcAuto::NuclideRelAct> &curve_acts : sol.m_rel_activities )
  {
    for( const RelActCalcAuto::NuclideRelAct &nra : curve_acts )
    {
      if( SpecUtils::icontains( nra.name(), needle ) && nra.rel_activity > 0.0 )
        return true;
    }
  }
  return false;
}


/** Loads `<n42_basename>` from `target/testing/test_data/RelActAutoReport/`, runs a full
 RelActCalcAuto solve using the embedded GUI state + DRF, and exercises the entire
 RelActAutoReport rendering pipeline (`solution_to_json` plus six templates: built-in
 `html`/`txt`/`json`, plus the three bundled README example templates).

 Asserts:
   - solve converged (`Status::Success`, `chi2/dof < 50`),
   - at least one fitted nuclide whose name contains `expected_nuclide_substring`,
   - all six templates render without throwing and produce non-empty output,
   - the JSON dump round-trips through `nlohmann::json::parse`,
   - the `expected_nuclide_substring` shows up in the rendered HTML / CSV / Markdown.
 */
static void do_one_fit_and_check( const string &n42_basename,
                                  const string &expected_nuclide_substring )
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path(
      SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" ), n42_basename );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spec_path ),
                         "Test spectrum not at '" << spec_path << "'" );

  // 1. Load the SpecMeas.  N42 carries fg/bg, embedded DRF, and the RelActAuto GUI state.
  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE_MESSAGE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ),
                         "Failed to load '" << spec_path << "'" );

  // 2. Pull out fg / bg.  Background is optional; many of these fits are foreground-only.
  shared_ptr<const SpecUtils::Measurement> foreground;
  shared_ptr<const SpecUtils::Measurement> background;
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

  // 3. Detector (DRF must be embedded in the N42).
  shared_ptr<const DetectorPeakResponse> det = meas->detector();
  BOOST_REQUIRE_MESSAGE( det, "N42 has no embedded detector" );

  // 4. RelActAuto GUI state must also be embedded.
  unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE_MESSAGE( state, "N42 missing embedded RelActAuto state" );

  // 5. Stamp filename + sample-numbers metadata onto Options so reports come out with the
  //    right "Filename, Samples: {…}" annotations.
  {
    const string fg_basename = SpecUtils::filename( spec_path );
    const set<int> fg_total = meas->sample_numbers();
    set<int> fg_samples, bg_samples;
    for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
    {
      if( !m )
        continue;
      if( m->source_type() == SpecUtils::SourceType::Background )
        bg_samples.insert( m->sample_number() );
      else if( m == foreground )
        fg_samples.insert( m->sample_number() );
    }
    const string bg_basename = background ? fg_basename : string();
    const set<int> bg_total   = background ? fg_total   : set<int>();
    RelActCalcAuto::set_input_spec_info( state->options,
                                         fg_basename, fg_total, fg_samples,
                                         bg_basename, bg_total, bg_samples );
  }

  // 6. Solve.
  const RelActCalcAuto::RelActAutoSolution sol
      = RelActCalcAuto::solve( state->options, foreground, background, det, {}, nullptr );

  // 7. Numerical sanity.
  BOOST_CHECK_MESSAGE( sol.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success,
                       "Solve did not succeed (status=" << static_cast<int>(sol.m_status)
                       << ", error='" << sol.m_error_message << "') for " << n42_basename );
  BOOST_CHECK( sol.m_dof > 0 );
  if( sol.m_dof > 0 )
    BOOST_CHECK_MESSAGE( sol.m_chi2 / sol.m_dof < 50.0,
                         "chi2/dof = " << (sol.m_chi2 / sol.m_dof)
                         << " is unreasonably large for " << n42_basename );

  // 7a. New goodness-of-fit stats: weighted R² (<= 1) and Jacobian condition number κ(J) (>= 1).
  BOOST_TEST_MESSAGE( "  " << n42_basename << ": R2 = " << sol.m_r2
                      << ", condition number kappa(J) = " << sol.m_jacobian_condition_number
                      << " (chi2/dof_data = "
                      << (sol.m_dof_data > 0 ? sol.m_chi2_data/static_cast<double>(sol.m_dof_data) : 0.0)
                      << ")" );
  BOOST_CHECK_MESSAGE( std::isnan(sol.m_r2) || (sol.m_r2 <= 1.0),
                       "R2 = " << sol.m_r2 << " exceeds 1 for " << n42_basename );
  BOOST_CHECK_MESSAGE( (sol.m_jacobian_condition_number < 0.0)
                         || (sol.m_jacobian_condition_number >= 1.0),
                       "condition number = " << sol.m_jacobian_condition_number
                       << " (< 1) for " << n42_basename );

  BOOST_CHECK_MESSAGE( solution_has_nuclide_named( sol, expected_nuclide_substring ),
                       "No nuclide with name containing '" << expected_nuclide_substring
                       << "' (and rel_activity > 0) in solution for " << n42_basename );

  // 7b. Uncertainty-reliability invariants: the data-only covariance rescale never deflates, the
  //     data-only DOF never exceeds the full DOF, and no reported enrichment uncertainty exhibits the
  //     rank-deficiency "variance collapse" (a sub-1e-4 relative uncertainty, which is now floored).
  BOOST_CHECK_MESSAGE( sol.m_cov_scale >= 1.0,
                       "m_cov_scale = " << sol.m_cov_scale << " (< 1) for " << n42_basename );
  BOOST_CHECK_MESSAGE( sol.m_dof_data <= sol.m_dof,
                       "m_dof_data (" << sol.m_dof_data << ") > m_dof (" << sol.m_dof << ") for " << n42_basename );

  if( !sol.m_covariance.empty() )
  {
    for( size_t curve = 0; curve < sol.m_rel_activities.size(); ++curve )
    {
      for( const RelActCalcAuto::NuclideRelAct &nra : sol.m_rel_activities[curve] )
      {
        const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( nra.source );
        if( !nuc || (nra.rel_activity <= 0.0) )
          continue;
        const pair<double,optional<double>> enr = sol.mass_enrichment_fraction( nuc, curve );
        if( !enr.second.has_value() || (enr.first <= 0.0) )
          continue;
        const double sigma = *enr.second;
        BOOST_CHECK_MESSAGE( std::isfinite(sigma) && (sigma >= 0.0),
                             "Enrichment uncertainty for " << nra.name() << " not finite/non-negative ("
                             << sigma << ") for " << n42_basename );
        // Skip single-isotope elements, whose fraction is legitimately ~1.0 with ~0 uncertainty.
        if( enr.first < 0.999 )
          BOOST_CHECK_MESSAGE( sigma >= 1.0e-4*enr.first,
                               "Enrichment uncertainty for " << nra.name() << " collapsed: sigma/value = "
                               << (sigma/enr.first) << " for " << n42_basename );
      }
    }
  }//if( !sol.m_covariance.empty() )

  // 8. JSON conversion.
  nlohmann::json data;
  BOOST_REQUIRE_NO_THROW( data = RelActAutoReport::solution_to_json( sol ) );
  BOOST_REQUIRE( data.contains("status") );
  BOOST_CHECK( data["status"].value("success", false) );

  // New goodness-of-fit fields must be present for the templates that reference them.
  BOOST_CHECK( data.contains("r2") );
  BOOST_CHECK( data.contains("condition_number") );
  BOOST_CHECK( data.contains("r2_str") );
  BOOST_CHECK( data.contains("condition_number_str") );

  // 9. Inja env construction.  (`inja::Environment` is non-copyable, so direct-construct.)
  inja::Environment env = RelActAutoReport::get_default_inja_env( "" );

  // 10. Smoke-test every template the README documents.  The first three are built-in shorthands
  //     resolved by name; the remaining three are concrete files in this test's data directory
  //     and are passed with their `include_dir` so `render_template` resolves them by filename.
  const string examples_dir = SpecUtils::append_path( g_test_file_dir, "RelActAutoReport" );
  struct TmpltCase { string tmplt; string include_dir; };
  const vector<TmpltCase> templates = {
    { "html", "" },
    { "txt",  "" },
    { "json", "" },
    { "example_csv.tmplt.csv",   examples_dir },
    { "example_md.tmplt.md",     examples_dir },
    { "example_peaks.tmplt.txt", examples_dir },
  };
  for( const TmpltCase &tc : templates )
  {
    string rendered;
    BOOST_REQUIRE_NO_THROW( rendered = RelActAutoReport::render_template( env, data, tc.tmplt, tc.include_dir ) );
    BOOST_CHECK_MESSAGE( !rendered.empty(),
                         "Template '" << tc.tmplt << "' rendered empty for " << n42_basename );

    if( tc.tmplt == "json" )
    {
      // The JSON dump should round-trip.
      BOOST_CHECK_NO_THROW( (void)nlohmann::json::parse( rendered ) );
    }else
    {
      // Every text-bearing template should mention the primary nuclide somewhere.
      BOOST_CHECK_MESSAGE( SpecUtils::icontains( rendered, expected_nuclide_substring ),
                           "Template '" << tc.tmplt << "' did not contain '"
                           << expected_nuclide_substring << "' for " << n42_basename );
    }
  }
}//do_one_fit_and_check


BOOST_AUTO_TEST_CASE( U235_unshielded )
{
  do_one_fit_and_check( "U235_Unshielded_6000.n42", "U235" );
}


BOOST_AUTO_TEST_CASE( Lu177m_shielded )
{
  do_one_fit_and_check( "ex_for_RelActAuto_report_Lu177m_Shielded.n42", "Lu177m" );
}


/** A solve that fails during setup returns a solution with `m_rel_eff_forms` filled in (it is
 populated from the options before validation runs) but `m_rel_eff_coefficients` still empty.
 Reporting on such a solution used to index one-past-the-end of `m_rel_eff_coefficients`, and the
 resulting garbage `vector` had a bogus `size()` that sent `RelActCalc::rel_eff_eqn_text` into an
 effectively endless string-append loop - the batch GUI just sat on "Performing Work" forever.
 */
BOOST_AUTO_TEST_CASE( report_of_failed_setup_solution )
{
  set_data_dir();

  RelActCalcAuto::RelActAutoSolution sol;
  sol.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
  sol.m_error_message = "Error initializing problem: No energy ranges defined.";

  RelActCalcAuto::RelEffCurveInput curve;
  curve.name = "Curve 0";
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 3;
  sol.m_options.rel_eff_curves.push_back( curve );

  sol.m_rel_eff_forms.push_back( RelActCalc::RelEffEqnForm::LnX );
  BOOST_REQUIRE( sol.m_rel_eff_coefficients.empty() );

  // The accessors must refuse, rather than read out of bounds.
  BOOST_CHECK_THROW( sol.rel_eff_txt( false, 0 ), std::exception );
  BOOST_CHECK_THROW( sol.rel_eff_eqn_js_function( 0 ), std::exception );
  BOOST_CHECK_THROW( sol.relative_efficiency( 186.0, 0 ), std::exception );

  nlohmann::json data;
  BOOST_REQUIRE_NO_THROW( data = RelActAutoReport::solution_to_json( sol ) );
  BOOST_REQUIRE( data.contains("status") );
  BOOST_CHECK( !data["status"].value("success", true) );

  BOOST_REQUIRE( data.contains("rel_eff_curves") );
  BOOST_REQUIRE_EQUAL( data["rel_eff_curves"].size(), 1 );
  const nlohmann::json &curve_json = data["rel_eff_curves"][0];
  BOOST_CHECK( curve_json.contains("equation_error") );
  BOOST_CHECK( !curve_json.contains("coefficients") );
  // The equation keys must still exist (empty), so templates can reference them unconditionally.
  BOOST_REQUIRE( curve_json.contains("equation_text") );
  BOOST_REQUIRE( curve_json.contains("equation_html") );
  BOOST_CHECK( curve_json["equation_text"].get<string>().empty() );
  BOOST_CHECK( curve_json["equation_html"].get<string>().empty() );

  // A FramPhysicalModel curve builds its equation from the cost functor / phys-model results, NOT
  //  from `m_rel_eff_coefficients` (a physical model with no correction function has nothing to put
  //  there), so the missing-coefficients gate must not apply to it.
  {
    RelActCalcAuto::RelActAutoSolution phys;
    phys.m_status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;

    RelActCalcAuto::RelEffCurveInput phys_curve;
    phys_curve.name = "Phys";
    phys_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;
    phys.m_options.rel_eff_curves.push_back( phys_curve );
    phys.m_rel_eff_forms.push_back( RelActCalc::RelEffEqnForm::FramPhysicalModel );
    BOOST_REQUIRE( phys.m_rel_eff_coefficients.empty() );

    // With no fit results there is no equation - but the complaint must be about the missing
    //  results, NOT about coefficients the physical model never uses.
    for( int is_js = 0; is_js < 2; ++is_js )
    {
      try
      {
        if( is_js )
          phys.rel_eff_eqn_js_function( 0 );
        else
          phys.rel_eff_txt( false, 0 );
        BOOST_CHECK_MESSAGE( false, "Physical model with no fit results should have thrown"
                                    " (is_js=" << is_js << ")" );
      }catch( std::exception &e )
      {
        const string msg = e.what();
        BOOST_CHECK_MESSAGE( msg.find("coefficients") == string::npos,
                             "Physical model should not be rejected for missing rel-eff"
                             " coefficients, but got: " << msg );
      }
    }//for( int is_js = 0; is_js < 2; ++is_js )

    // ...and the report must report that, not print a sentinel string into the user's equation.
    nlohmann::json phys_data;
    BOOST_REQUIRE_NO_THROW( phys_data = RelActAutoReport::solution_to_json( phys ) );
    BOOST_REQUIRE_EQUAL( phys_data["rel_eff_curves"].size(), 1 );
    BOOST_CHECK( phys_data["rel_eff_curves"][0]["equation_text"].get<string>().empty() );
    BOOST_CHECK( phys_data["rel_eff_curves"][0].contains("equation_error") );
  }

  // The report templates must still render for a failed solution.
  inja::Environment env = RelActAutoReport::get_default_inja_env( "" );
  for( const string &tmplt : vector<string>{ "html", "txt", "json" } )
  {
    string rendered, err;
    try
    {
      rendered = RelActAutoReport::render_template( env, data, tmplt, "" );
    }catch( std::exception &e )
    {
      err = e.what();
    }
    BOOST_CHECK_MESSAGE( err.empty(), "Template '" << tmplt << "' threw: " << err );
    BOOST_CHECK_MESSAGE( !rendered.empty(), "Template '" << tmplt << "' rendered empty." );
  }
}//BOOST_AUTO_TEST_CASE( report_of_failed_setup_solution )


/** End-to-end version of the above: run the real `solve(...)` with options that dont define a
 problem, then push the returned solution through every reporting entry point.  This covers the
 solution shape `solve` actually produces on an early setup failure (which differs from the
 hand-built one above), and would previously hang or read out of bounds.

 Note: against the unfixed code these checks HANG (an unbounded string append) rather than fail,
 so a CI timeout here is a real regression, not flake.
 */
BOOST_AUTO_TEST_CASE( reporting_on_unusable_solve )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );

  // A curve with a nuclide, but no energy ranges anywhere - i.e. the state of an "Isotopics by
  //  nuclides" tool that was opened and given a nuclide, but never given an energy range.
  RelActCalcAuto::Options options;
  RelActCalcAuto::RelEffCurveInput curve;
  curve.name = "Curve 0";
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 3;
  RelActCalcAuto::NucInputInfo nuc;
  nuc.source = db->nuclide( "U235" );
  BOOST_REQUIRE( RelActCalcAuto::nuclide(nuc.source) );
  curve.nuclides.push_back( nuc );
  options.rel_eff_curves.push_back( curve );
  BOOST_REQUIRE( !options.why_not_usable().empty() );

  // Something to hand `solve` as a foreground; its contents dont matter, we never get that far.
  std::shared_ptr<SpecUtils::Measurement> foreground = std::make_shared<SpecUtils::Measurement>();
  foreground->set_gamma_counts( std::make_shared<vector<float>>( 1024, 1.0f ), 300.0f, 300.0f );

  RelActCalcAuto::RelActAutoSolution sol;
  BOOST_REQUIRE_NO_THROW( sol = RelActCalcAuto::solve( options, foreground, nullptr, nullptr,
                                                       /*all_peaks=*/{}, /*cancel_calc=*/nullptr ) );
  BOOST_CHECK( sol.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success );
  BOOST_CHECK( !sol.m_error_message.empty() );

  // Every reporting surface must cope with the failed solution.
  std::stringstream summary_strm;
  BOOST_CHECK_NO_THROW( sol.print_summary( summary_strm ) );

  if( html_report_resources_available() )
  {
    std::stringstream html_strm;
    BOOST_CHECK_NO_THROW( sol.print_html_report( html_strm ) );
  }else
  {
    BOOST_TEST_MESSAGE( "Skipping print_html_report(): 'InterSpec_resources' is not in the CWD." );
  }

  nlohmann::json data;
  BOOST_REQUIRE_NO_THROW( data = RelActAutoReport::solution_to_json( sol ) );
  BOOST_CHECK( !data["status"].value("success", true) );

  inja::Environment env = RelActAutoReport::get_default_inja_env( "" );
  for( const string &tmplt : vector<string>{ "html", "txt", "json" } )
  {
    string rendered, err;
    try
    {
      rendered = RelActAutoReport::render_template( env, data, tmplt, "" );
    }catch( std::exception &e )
    {
      err = e.what();
    }
    BOOST_CHECK_MESSAGE( err.empty(), "Template '" << tmplt << "' threw: " << err );
  }

  // The above fails before any per-curve results are filled in, so every reporting loop is simply
  //  empty.  The more dangerous shape is a setup failure that happens AFTER `m_rel_eff_forms` is
  //  populated - then `m_rel_eff_forms` is non-empty while `m_rel_activities` /
  //  `m_rel_eff_coefficients` are still empty, and any loop bounded on the wrong vector reads out
  //  of bounds.  A zero live-time foreground fails at exactly that point.
  {
    RelActCalcAuto::Options ok_options = options;
    RelActCalcAuto::RoiRange roi;
    roi.lower_energy = 120.0;
    roi.upper_energy = 220.0;
    ok_options.rois.push_back( roi );
    BOOST_REQUIRE( ok_options.why_not_usable().empty() );

    std::shared_ptr<SpecUtils::Measurement> zero_lt = std::make_shared<SpecUtils::Measurement>();
    zero_lt->set_gamma_counts( std::make_shared<vector<float>>( 1024, 1.0f ), 0.0f, 0.0f );

    RelActCalcAuto::RelActAutoSolution late_fail;
    BOOST_REQUIRE_NO_THROW( late_fail = RelActCalcAuto::solve( ok_options, zero_lt, nullptr, nullptr,
                                                               /*all_peaks=*/{}, /*cancel_calc=*/nullptr ) );
    BOOST_CHECK( late_fail.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success );
    BOOST_CHECK_MESSAGE( !late_fail.m_rel_eff_forms.empty(),
                         "Expected this failure to occur after m_rel_eff_forms was populated;"
                         " if solve() changed, this test no longer covers what it means to." );
    BOOST_CHECK( late_fail.m_rel_activities.empty() );
    BOOST_CHECK( late_fail.m_rel_eff_coefficients.empty() );

    std::stringstream late_summary, late_html;
    BOOST_CHECK_NO_THROW( late_fail.print_summary( late_summary ) );
    if( html_report_resources_available() )
      BOOST_CHECK_NO_THROW( late_fail.print_html_report( late_html ) );

    nlohmann::json late_data;
    BOOST_REQUIRE_NO_THROW( late_data = RelActAutoReport::solution_to_json( late_fail ) );
    BOOST_REQUIRE_EQUAL( late_data["rel_eff_curves"].size(), 1 );
    BOOST_CHECK( late_data["rel_eff_curves"][0].contains("equation_error") );

    for( const string &tmplt : vector<string>{ "html", "txt", "json" } )
    {
      string err;
      try
      {
        RelActAutoReport::render_template( env, late_data, tmplt, "" );
      }catch( std::exception &e )
      {
        err = e.what();
      }
      BOOST_CHECK_MESSAGE( err.empty(), "Template '" << tmplt << "' threw on late failure: " << err );
    }
  }
}//BOOST_AUTO_TEST_CASE( reporting_on_unusable_solve )


/** `Options::why_not_usable()` is what keeps an unconfigured "Isotopics by nuclides" state from
 being handed to `solve(...)` (and to the batch tools).
 */
BOOST_AUTO_TEST_CASE( options_why_not_usable )
{
  set_data_dir();

  RelActCalcAuto::Options options;
  BOOST_CHECK( !options.why_not_usable().empty() );  //no rel eff curves

  RelActCalcAuto::RelEffCurveInput curve;
  curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
  curve.rel_eff_eqn_order = 3;
  options.rel_eff_curves.push_back( curve );
  BOOST_CHECK( !options.why_not_usable().empty() );  //no nuclides

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  RelActCalcAuto::NucInputInfo nuc;
  nuc.source = db->nuclide( "U235" );
  BOOST_REQUIRE( RelActCalcAuto::nuclide(nuc.source) );
  options.rel_eff_curves[0].nuclides.push_back( nuc );
  BOOST_CHECK( !options.why_not_usable().empty() );  //no energy ranges

  RelActCalcAuto::RoiRange roi;
  roi.lower_energy = 120.0;
  roi.upper_energy = 220.0;
  options.rois.push_back( roi );
  BOOST_CHECK( options.why_not_usable().empty() );
}//BOOST_AUTO_TEST_CASE( options_why_not_usable )


/** The bundled `std_rel_eff_summary` templates' multi-curve separation block is skipped for every
 solution the file-based cases above produce (single-curve fits have `curve_separation_status`
 "NotApplicable"), so its inja was previously never executed by tests - and this inja errors on
 unknown keys.  Render it against a synthetic multi-curve solution: the separation metrics are
 public members, so no fit is needed (numbers replicate the 2026-08 enriched-U-inside-lower-
 enriched-U case: marginal z's ~2 corroborated by a decisive merged-fit rejection).
 */
BOOST_AUTO_TEST_CASE( multi_curve_separation_block_renders )
{
  set_data_dir();

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const u235 = db->nuclide( "U235" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide( "U238" );
  BOOST_REQUIRE( u235 && u238 );

  RelActCalcAuto::RelActAutoSolution sol;
  sol.m_status = RelActCalcAuto::RelActAutoSolution::Status::Success;
  sol.m_options.rel_eff_curves.resize( 2 );
  sol.m_options.rel_eff_curves[0].name = "Inner";
  sol.m_options.rel_eff_curves[1].name = "Outer";
  sol.m_chi2_data = 571.0;
  sol.m_dof_data = 383;
  sol.m_r2 = 0.999;
  sol.m_jacobian_condition_number = 148.0;

  sol.m_evidence_purity.assign( 2, {} );
  sol.m_evidence_purity[0][u235] = 0.011;
  sol.m_evidence_purity[0][u238] = 0.46;
  sol.m_evidence_purity[1][u235] = 0.94;
  sol.m_evidence_purity[1][u238] = 0.54;

  sol.m_source_model_counts.assign( 2, {} );
  sol.m_source_model_counts[0][u235] = 500.0;
  sol.m_source_model_counts[0][u238] = 30000.0;
  sol.m_source_model_counts[1][u235] = 45000.0;
  sol.m_source_model_counts[1][u238] = 35000.0;

  RelActCalcAuto::RelActAutoSolution::CrossCurveCorrelation corr;
  corr.curve_a = 0;
  corr.curve_b = 1;
  corr.param_a = corr.param_b = "Act(U235)";
  corr.correlation = 0.5;
  sol.m_cross_curve_correlations.push_back( corr );
  sol.m_cross_curve_max_corr = corr;

  RelActCalcAuto::RelActAutoSolution::EnrichmentDiffZ diff;
  diff.curve_a = 0;
  diff.curve_b = 1;
  diff.nuclide = RelActCalcAuto::SrcVariant( u235 );
  diff.enrichment_a = 0.775;
  diff.sigma_a = 0.18;
  diff.enrichment_b = 0.00896;
  diff.sigma_b = 0.33;
  diff.z = 2.05;
  diff.reliable = true;
  sol.m_enrichment_diff_z.push_back( diff );

  RelActCalcAuto::RelActAutoSolution::MergedCurveComparison merged;
  merged.valid = true;
  merged.multi_chi2_data = 571.0;
  merged.multi_dof_data = 383;
  merged.merged_chi2_data = 586.17;
  merged.merged_dof_data = 385;
  merged.delta_chi2 = 15.17;
  merged.extra_dof_of_multi = 2;
  sol.m_merged_single_curve_comparison = merged;

  sol.finalize_curve_separation_status();
  BOOST_REQUIRE( sol.curves_detected_distinct() );

  const nlohmann::json data = RelActAutoReport::solution_to_json( sol );

  // The new additive keys.
  BOOST_REQUIRE( data.contains( "curves_distinct_basis" ) );
  BOOST_CHECK_EQUAL( data["curves_distinct_basis"].get<string>(), string("z_plus_merged") );
  BOOST_REQUIRE( data.contains( "source_count_attribution" ) );
  BOOST_CHECK_EQUAL( data["source_count_attribution"].size(), size_t(2) );
  BOOST_REQUIRE( data.contains( "evidence_purity" ) );
  BOOST_REQUIRE_EQUAL( data["evidence_purity"].size(), size_t(2) );
  BOOST_REQUIRE( !data["evidence_purity"][0].empty() );
  BOOST_CHECK( data["evidence_purity"][0][0].contains( "purity_str" ) );   //pre-existing key kept
  BOOST_CHECK( data["evidence_purity"][0][0].contains( "curve_label" ) );
  BOOST_CHECK( data["evidence_purity"][0][0].contains( "model_counts" ) );

  // Locate the bundled templates directory relative to the test-data dir (repo layout), with CWD
  //  fallbacks for running by hand from the build directory.
  string tmplt_dir;
  const string candidates[] = {
    SpecUtils::append_path( g_test_file_dir, "../../../InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts" ),
    "InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts",
    "../../InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts",
  };
  for( const string &candidate : candidates )
  {
    if( SpecUtils::is_file( SpecUtils::append_path( candidate, "std_rel_eff_summary.tmplt.txt" ) ) )
    {
      tmplt_dir = candidate;
      break;
    }
  }
  BOOST_REQUIRE_MESSAGE( !tmplt_dir.empty(), "Could not locate IsotopicsByNuclidesReportTmplts dir" );

  inja::Environment env = RelActAutoReport::get_default_inja_env( "" );
  for( const string tmplt : { string("std_rel_eff_summary.tmplt.txt"),
                              string("std_rel_eff_summary.tmplt.html") } )
  {
    string rendered, err;
    try
    {
      rendered = RelActAutoReport::render_template( env, data, tmplt, tmplt_dir );
    }catch( std::exception &e )
    {
      err = e.what();
    }
    BOOST_REQUIRE_MESSAGE( err.empty(), "Template '" << tmplt << "' threw: " << err );

    BOOST_CHECK_MESSAGE( rendered.find( "Peak-count attribution per nuclide" ) != string::npos,
                         "Template '" << tmplt << "' missing attribution table" );
    BOOST_CHECK_MESSAGE( rendered.find( "'Inner'" ) != string::npos,
                         "Template '" << tmplt << "' missing curve label" );
    BOOST_CHECK_MESSAGE( rendered.find( "Attributed share" ) != string::npos,
                         "Template '" << tmplt << "' missing attributed-share legend" );
    BOOST_CHECK_MESSAGE( rendered.find( "Evidence purity" ) == string::npos,
                         "Template '" << tmplt << "' still says 'Evidence purity'" );
  }
}//BOOST_AUTO_TEST_CASE( multi_curve_separation_block_renders )
