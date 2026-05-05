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

  BOOST_CHECK_MESSAGE( solution_has_nuclide_named( sol, expected_nuclide_substring ),
                       "No nuclide with name containing '" << expected_nuclide_substring
                       << "' (and rel_activity > 0) in solution for " << n42_basename );

  // 8. JSON conversion.
  nlohmann::json data;
  BOOST_REQUIRE_NO_THROW( data = RelActAutoReport::solution_to_json( sol ) );
  BOOST_REQUIRE( data.contains("status") );
  BOOST_CHECK( data["status"].value("success", false) );

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
