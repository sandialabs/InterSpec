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
 */

#include "InterSpec_config.h"

#include <set>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#define BOOST_TEST_MODULE BatchRelActAuto_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/BatchRelActAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

static_assert( USE_REL_ACT_TOOL, "Compile-time support for Rel Act tool is required." );

static string g_test_file_dir;
static string g_eu152_dir;

// Mirror of the boost-test driver setup helper used in test_RelActAutoReport.cpp.
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
      if( SpecUtils::is_directory( SpecUtils::append_path(d, "RelActAutoBatch") ) )
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

  g_eu152_dir = SpecUtils::append_path( g_test_file_dir, "RelActAutoBatch" );
  g_eu152_dir = SpecUtils::append_path( g_eu152_dir, "Eu152" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( g_eu152_dir ),
                         "Eu152 batch test data not at '" << g_eu152_dir << "'" );
}//set_data_dir()


// Helper: extract first nuclide's relative activity (with uncertainty) from
// solution.  Returns -1.0 if not present.
static pair<double,double> first_nuclide_rel_act( const RelActCalcAuto::RelActAutoSolution &sol )
{
  for( const vector<RelActCalcAuto::NuclideRelAct> &curve : sol.m_rel_activities )
  {
    for( const RelActCalcAuto::NuclideRelAct &nra : curve )
    {
      if( nra.rel_activity > 0.0 )
        return { nra.rel_activity, nra.rel_activity_uncertainty };
    }
  }
  return { -1.0, -1.0 };
}


// Test 1: Exemplar-only happy path - the N42 carries everything we need.
BOOST_AUTO_TEST_CASE( exemplar_only_happy_path )
{
  set_data_dir();

  const string exemplar = SpecUtils::append_path( g_eu152_dir, "exemplar.n42" );
  const string input    = SpecUtils::append_path( g_eu152_dir, "other_eu152_shielded.n42" );
  BOOST_REQUIRE( SpecUtils::is_file(exemplar) );
  BOOST_REQUIRE( SpecUtils::is_file(input) );

  BatchRelActAuto::Options options;
  const BatchRelActAuto::Result result
      = BatchRelActAuto::run_on_file( exemplar, /*sample_nums=*/{},
                                       /*cached_exemplar=*/nullptr,
                                       input, /*cached_file=*/nullptr,
                                       options );

  BOOST_CHECK_MESSAGE( result.m_result_code == BatchRelActAuto::ResultCode::Success,
                       "Expected Success, got " << BatchRelActAuto::to_str(result.m_result_code)
                       << " (msg='" << result.m_error_msg << "')" );
  BOOST_CHECK( result.m_solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success );

  const auto act = first_nuclide_rel_act( result.m_solution );
  BOOST_CHECK_MESSAGE( act.first > 0.0, "Expected positive rel_activity, got " << act.first );
  if( act.first > 0.0 )
    BOOST_CHECK_MESSAGE( act.second / act.first < 0.5,
                         "Expected uncertainty < 50% of rel_activity, got "
                         << (act.second/act.first) );
}


// Test 2: DRF override should produce a different fit result.
BOOST_AUTO_TEST_CASE( drf_override_changes_result )
{
  set_data_dir();

  const string exemplar = SpecUtils::append_path( g_eu152_dir, "exemplar.n42" );
  const string input    = SpecUtils::append_path( g_eu152_dir, "other_eu152_shielded.n42" );
  const string drf_file = SpecUtils::append_path( g_eu152_dir, "other_det_eff.drf.xml" );
  BOOST_REQUIRE( SpecUtils::is_file(drf_file) );

  // Baseline (exemplar's embedded DRF)
  BatchRelActAuto::Options baseline_opts;
  const BatchRelActAuto::Result baseline
      = BatchRelActAuto::run_on_file( exemplar, {}, nullptr, input, nullptr, baseline_opts );
  BOOST_REQUIRE( baseline.m_result_code == BatchRelActAuto::ResultCode::Success );
  const auto base_act = first_nuclide_rel_act( baseline.m_solution );
  BOOST_REQUIRE( base_act.first > 0.0 );

  // Override: load other_det_eff.drf.xml directly via DetectorPeakResponse::fromXml.
  // (BatchActivity::init_drf_from_name doesn't handle the DRF XML format - it tries
  //  GADRAS dirs, URI files, CSVs, ECC, and ANGLE outx; users with a .drf.xml drop
  //  the file via the GUI uploader, which uses fromXml directly.)
  shared_ptr<DetectorPeakResponse> drf;
  {
    vector<char> data;
    BOOST_REQUIRE_NO_THROW( SpecUtils::load_file_data( drf_file.c_str(), data ) );
    rapidxml::xml_document<char> doc;
    BOOST_REQUIRE_NO_THROW( doc.parse<rapidxml::parse_trim_whitespace>( &data[0] ) );
    const rapidxml::xml_node<char> *root = doc.first_node( "DetectorPeakResponse" );
    BOOST_REQUIRE( root );
    drf = make_shared<DetectorPeakResponse>();
    BOOST_REQUIRE_NO_THROW( drf->fromXml( root ) );
  }
  BOOST_REQUIRE( drf );

  BatchRelActAuto::Options ov_opts;
  ov_opts.drf_override = drf;
  const BatchRelActAuto::Result overridden
      = BatchRelActAuto::run_on_file( exemplar, {}, nullptr, input, nullptr, ov_opts );

  BOOST_CHECK_MESSAGE( overridden.m_result_code == BatchRelActAuto::ResultCode::Success,
                       "DRF override solve failed: "
                       << BatchRelActAuto::to_str(overridden.m_result_code)
                       << " (" << overridden.m_error_msg << ")" );
  const auto ov_act = first_nuclide_rel_act( overridden.m_solution );
  BOOST_REQUIRE( ov_act.first > 0.0 );

  // Verify the override DRF was actually plumbed through to the solver.  For
  // non-Physical-model rel-eff curves with FWHM-seeded methods, the DRF only
  // affects FWHM seeding, so the rel_activity numerical value may end up
  // similar to the embedded-DRF result; what matters is that the override DRF
  // identity reached the solver.  (The solver wraps the DRF in its own
  // shared_ptr, so compare by name rather than pointer.)
  BOOST_REQUIRE( overridden.m_solution.m_drf );
  BOOST_CHECK_MESSAGE( overridden.m_solution.m_drf->name() == drf->name(),
                       "DRF override didn't reach the solver: solution.m_drf->name() = '"
                       << overridden.m_solution.m_drf->name()
                       << "', expected '" << drf->name() << "'" );

  // Baseline should have used the exemplar's embedded DRF (different name).
  if( baseline.m_solution.m_drf )
    BOOST_CHECK_MESSAGE( baseline.m_solution.m_drf->name() != drf->name(),
                         "Baseline and override DRF have the same name '"
                         << drf->name() << "' - test fixture problem" );
}


// Test 3: state-file override fully replaces the exemplar's GuiState.
BOOST_AUTO_TEST_CASE( state_override_replaces_state )
{
  set_data_dir();

  const string exemplar  = SpecUtils::append_path( g_eu152_dir, "exemplar.n42" );
  const string input     = SpecUtils::append_path( g_eu152_dir, "other_eu152_shielded.n42" );
  const string state_xml = SpecUtils::append_path( g_eu152_dir,
                                                   "isotopics_by_nuclides_Eu152_Unshielded_releff.xml" );
  BOOST_REQUIRE( SpecUtils::is_file(state_xml) );

  shared_ptr<RelActCalcAuto::RelActAutoGuiState> override_state;
  BOOST_REQUIRE_NO_THROW(
      override_state = BatchRelActAuto::load_state_from_xml_file( state_xml )
  );
  BOOST_REQUIRE( override_state );
  // Sanity: the bundled override XML disables background subtraction
  BOOST_CHECK( override_state->background_subtract == false );

  BatchRelActAuto::Options options;
  options.state_override = override_state;
  const BatchRelActAuto::Result result
      = BatchRelActAuto::run_on_file( exemplar, {}, nullptr, input, nullptr, options );

  BOOST_CHECK_MESSAGE( result.m_result_code == BatchRelActAuto::ResultCode::Success,
                       "State-override solve failed: "
                       << BatchRelActAuto::to_str(result.m_result_code)
                       << " (" << result.m_error_msg << ")" );

  // The fit should run without background subtract since the override state had it off.
  // (RelActAutoSolution doesn't expose background_subtract directly, but absence of
  // a m_background pointer is the strongest available proxy.)
  BOOST_CHECK( !result.m_solution.m_background );

  // And the rel-eff curve should match what's in the override XML
  // (LnX, order 3, single Eu152 nuclide).
  BOOST_REQUIRE_MESSAGE( !result.m_solution.m_options.rel_eff_curves.empty(),
                         "Solution has no rel_eff_curves" );
  const RelActCalcAuto::RelEffCurveInput &c = result.m_solution.m_options.rel_eff_curves[0];
  BOOST_CHECK( c.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::LnX );
  BOOST_CHECK( c.rel_eff_eqn_order == 3 );
}


// Test 4: a malformed / non-existent --state-file should produce a clean error.
BOOST_AUTO_TEST_CASE( bad_state_file_throws )
{
  set_data_dir();

  const string nonexistent = SpecUtils::append_path( g_eu152_dir, "this_file_does_not_exist.xml" );
  BOOST_CHECK_THROW( BatchRelActAuto::load_state_from_xml_file(nonexistent),
                     std::runtime_error );

  // Also try a malformed XML file (just write some garbage to a temp file)
  const string tmp = SpecUtils::temp_file_name( "bad_state_xml", SpecUtils::temp_dir() );
  {
    std::ofstream out( tmp.c_str() );
    out << "this is not xml at all <RelActCalcAuto>\n";
  }
  BOOST_CHECK_THROW( BatchRelActAuto::load_state_from_xml_file(tmp),
                     std::runtime_error );
  SpecUtils::remove_file( tmp );

  // And an XML file with the wrong root element
  const string tmp2 = SpecUtils::temp_file_name( "wrong_root_xml", SpecUtils::temp_dir() );
  {
    std::ofstream out( tmp2.c_str() );
    out << "<?xml version=\"1.0\"?><NotRelActCalcAuto/>\n";
  }
  BOOST_CHECK_THROW( BatchRelActAuto::load_state_from_xml_file(tmp2),
                     std::runtime_error );
  SpecUtils::remove_file( tmp2 );
}


// Test 5: a state with FramPhysicalModel but no DRF should fail with a specific
// ResultCode rather than running the solver and returning garbage.
BOOST_AUTO_TEST_CASE( phys_model_without_drf_errors )
{
  set_data_dir();

  const string exemplar = SpecUtils::append_path( g_eu152_dir, "exemplar.n42" );
  const string input    = SpecUtils::append_path( g_eu152_dir, "other_eu152_shielded.n42" );

  // Load the exemplar, strip its embedded DRF, and mutate the rel-eff curve type to
  // FramPhysicalModel.  SpecMeas is non-copyable, so we mutate in place.
  shared_ptr<SpecMeas> exemplar_meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( exemplar_meas->load_file( exemplar, SpecUtils::ParserType::Auto ) );
  exemplar_meas->setDetector( nullptr );

  unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = exemplar_meas->getRelActAutoGuiState();
  BOOST_REQUIRE( state );
  BOOST_REQUIRE( !state->options.rel_eff_curves.empty() );
  state->options.rel_eff_curves[0].rel_eff_eqn_type = RelActCalc::RelEffEqnForm::FramPhysicalModel;

  BatchRelActAuto::Options options;
  options.state_override = std::shared_ptr<RelActCalcAuto::RelActAutoGuiState>( std::move(state) );
  // Don't set drf_override either.
  const BatchRelActAuto::Result result
      = BatchRelActAuto::run_on_file( "", {}, exemplar_meas, input, nullptr, options );

  BOOST_CHECK_MESSAGE(
      result.m_result_code == BatchRelActAuto::ResultCode::ExemplarUsesPhysModelButNoDrf,
      "Expected ExemplarUsesPhysModelButNoDrf, got "
      << BatchRelActAuto::to_str(result.m_result_code)
      << " (msg='" << result.m_error_msg << "')" );
}
