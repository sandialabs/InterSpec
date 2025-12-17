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

#include <string>
#include <iostream>

#include <Wt/Utils>
#include <Wt/WApplication>
#include <Wt/Test/WTestEnvironment>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE LlmToolRegistry_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using namespace boost::unit_test;
using json = nlohmann::json;

// TODO: as of 20251013, need to add tests for: `currie_mda_calc`, `add_analysis_peaks_for_source`, `automated_source_id_results`, and `nuclides_with_primary_gammas_in_energy_range` 


std::string g_test_file_dir;

// We need to set the static data directory, so the code knows where
//  like sandia.decay.xml is located.
void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;

  s_have_set = true;

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  string datadir;

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );

    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )

  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )

  const string sandia_deacay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_deacay_file ), "sandia.decay.xml not at '" << sandia_deacay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  // Set the writable data directory for tests - use the same as static data directory for simplicity
  // This is needed for tests that use InterSpec::writableDataDirectory() like list_isotopics_presets
  BOOST_REQUIRE_NO_THROW( InterSpec::setWritableDataDirectory( datadir ) );

  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
}//void set_data_dir()


// Helper class to manage InterSpec instance for tests
class InterSpecTestFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
  InterSpec *m_interspec;
  std::shared_ptr<LlmTools::ToolRegistry> m_tool_registry;


  InterSpecTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr ),
    m_tool_registry( nullptr )
  {
    set_data_dir();

    string wt_app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), "..");
    wt_app_root = SpecUtils::lexically_normalize_path(wt_app_root);

    // Create a test environment
    const std::string applicationPath = "";
    const std::string configurationFile = ""; //Add a XML config file that makes it so WLogger wont printout anything
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) ); //so wApp is avaiable everywhere

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    // Load LLM configuration
    std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
    BOOST_REQUIRE( llmConfig );

    // Create LlmTool for old tests that use llm_gui
    m_interspec->createLlmTool();

    // Also create a direct tool registry for new isotopics tests
    // This allows testing tools even when the LLM API is not enabled
    m_tool_registry = std::make_shared<LlmTools::ToolRegistry>( *llmConfig );
    BOOST_REQUIRE( m_tool_registry );

    // Load a test spectrum file for detector-related tests
    try
    {
      const string datadir = InterSpec::staticDataDirectory();

      {//Begin load foreground
        const string spectrum_file = SpecUtils::append_path( datadir, "reference_spectra/Common_Field_Nuclides/Detective X/Br82_Unshielded.txt" );
        BOOST_REQUIRE( SpecUtils::is_file(spectrum_file) );

        std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
        const bool loaded = meas->load_file( spectrum_file, SpecUtils::ParserType::Auto, spectrum_file );

        if( loaded && (meas->num_measurements() > 0) )
        {
          const shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index(0);
          BOOST_REQUIRE( !!m );
          m_interspec->setSpectrum( meas, {m->sample_number()}, SpecUtils::SpectrumType::Foreground, 0 );
        }
        else
        {
          BOOST_TEST_MESSAGE( "Failed to load test spectrum: " << spectrum_file );
        }
      }//End load foreground


      {//Begin load background
        const string background_file = SpecUtils::append_path( datadir, "reference_spectra/Common_Field_Nuclides/Detective X/background.txt" );
        BOOST_REQUIRE( SpecUtils::is_file(background_file) );

        std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
        const bool loaded = meas->load_file( background_file, SpecUtils::ParserType::Auto, background_file );

        if( loaded && (meas->num_measurements() > 0) )
        {
          const shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index(0);
          BOOST_REQUIRE( !!m );
          m_interspec->setSpectrum( meas, {m->sample_number()}, SpecUtils::SpectrumType::Background, 0 );
        }
        else
        {
          BOOST_TEST_MESSAGE( "Failed to load test background spectrum: " << background_file );
        }
      }//End load background
    }catch( const std::exception &e )
    {
      BOOST_TEST_MESSAGE( "Exception loading test spectrum: " << e.what() );
    }
  }

  ~InterSpecTestFixture()
  {
    m_tool_registry.reset();
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr; //deleted when m_env is deleted
  }


  const LlmTools::ToolRegistry &llmToolRegistry()
  {
    BOOST_REQUIRE( m_tool_registry );
    return *m_tool_registry;
  }

};




BOOST_AUTO_TEST_CASE( test_executeListIsotopicsPresets )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // List presets - no parameters needed
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  // Should return an object with "presets" array
  BOOST_CHECK( result.contains( "presets" ) );
  BOOST_CHECK( result["presets"].is_array() );

  // There should be at least some presets available
  BOOST_CHECK_MESSAGE( result["presets"].size() > 0, "Should have at least some isotopics presets" );

  // Check that presets have expected structure
  bool found_pu_preset = false;
  bool found_u_preset = false;

  for( const auto &preset : result["presets"] )
  {
    // Each preset should have a name
    BOOST_CHECK( preset.contains( "name" ) );
    BOOST_CHECK( preset["name"].is_string() );

    const string preset_name = preset["name"].get<string>();

    if( SpecUtils::icontains( preset_name, "Pu" ) || SpecUtils::icontains( preset_name, "plutonium" ) )
      found_pu_preset = true;
    if( SpecUtils::icontains( preset_name, "U" ) || SpecUtils::icontains( preset_name, "uranium" ) )
      found_u_preset = true;
  }

  BOOST_CHECK_MESSAGE( found_pu_preset || found_u_preset, "Should have at least one Pu or U preset" );

  cout << "Found " << result["presets"].size() << " isotopics presets" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeLoadIsotopicsPreset )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First list presets to get a valid preset name
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );
  BOOST_REQUIRE( list_result["presets"].is_array() );
  BOOST_REQUIRE( list_result["presets"].size() > 0 );

  // Find a preset that isn't the current configuration
  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;

    preset_name = preset["name"].get<string>();
    break;
  }

  BOOST_REQUIRE_MESSAGE( !preset_name.empty(), "Should have at least one non-current preset to load" );

  // Load the preset
  json params;
  params["preset"] = preset_name;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "load_isotopics_preset", params, fixture.m_interspec ) );

  // Check success
  BOOST_CHECK( result.contains( "success" ) );
  BOOST_CHECK( result["success"].get<bool>() == true );

  // Should have configuration details
  BOOST_CHECK( result.contains( "configuration" ) );

  cout << "Successfully loaded preset: " << preset_name << endl;
}


BOOST_AUTO_TEST_CASE( test_executeGetIsotopicsConfig )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First load a preset to ensure we have a configuration
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;
    preset_name = preset["name"].get<string>();
    break;
  }

  if( !preset_name.empty() )
  {
    json load_params;
    load_params["preset"] = preset_name;
    BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );
  }

  // Now get the configuration
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  // If we loaded a preset, should have a config
  if( !preset_name.empty() )
  {
    BOOST_CHECK( result.contains( "has_config" ) );
    BOOST_CHECK( result["has_config"].get<bool>() == true );

    // Check for expected fields
    BOOST_CHECK( result.contains( "nuclides" ) );
    BOOST_CHECK( result["nuclides"].is_array() );

    BOOST_CHECK( result.contains( "rois" ) );
    BOOST_CHECK( result["rois"].is_array() );

    cout << "Config has " << result["nuclides"].size() << " nuclides and "
         << result["rois"].size() << " ROIs" << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_executeGetIsotopicsConfigSchema )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_isotopics_config_schema", json::object(), fixture.m_interspec ) );

  // Check for expected sections
  BOOST_CHECK( result.contains( "structure" ) );
  BOOST_CHECK( result.contains( "rel_eff_curve" ) );
  BOOST_CHECK( result.contains( "nuclides" ) );
  BOOST_CHECK( result.contains( "rois" ) );
  BOOST_CHECK( result.contains( "global_options" ) );
  BOOST_CHECK( result.contains( "constraints" ) );
  BOOST_CHECK( result.contains( "shielding" ) );

  // Check for enum values
  BOOST_CHECK( result.contains( "rel_eff_eqn_types" ) );
  BOOST_CHECK( result["rel_eff_eqn_types"].is_array() );
  BOOST_CHECK( result["rel_eff_eqn_types"].size() >= 4 ); // LnX, LnY, LnXLnY, FramPhysicalModel

  BOOST_CHECK( result.contains( "fwhm_forms" ) );
  BOOST_CHECK( result["fwhm_forms"].is_array() );

  BOOST_CHECK( result.contains( "skew_types" ) );
  BOOST_CHECK( result["skew_types"].is_array() );

  BOOST_CHECK( result.contains( "continuum_types" ) );
  BOOST_CHECK( result["continuum_types"].is_array() );

  BOOST_CHECK( result.contains( "success" ) );
  BOOST_CHECK( result["success"].get<bool>() == true );

  cout << "Schema includes " << result["rel_eff_eqn_types"].size() << " equation types, "
       << result["fwhm_forms"].size() << " FWHM forms, "
       << result["skew_types"].size() << " skew types" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_Add )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First reset any existing config
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Add nuclides to a fresh configuration
  json params;
  params["action"] = "add";
  params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  nuc1["age"] = "15 years";
  params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "Pu240";
  nuc2["age"] = "15 years";
  params["nuclides"].push_back( nuc2 );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_nuclides", params, fixture.m_interspec ) );

  BOOST_CHECK( result.contains( "success" ) );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result.contains( "nuclide_count" ) );
  BOOST_CHECK( result["nuclide_count"].get<int>() == 2 );

  // Verify by getting config
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );
  BOOST_CHECK( config["nuclides"].size() == 2 );

  cout << "Added 2 nuclides successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_AddDuplicate )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add a nuclide
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json params;
  params["action"] = "add";
  params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  params["nuclides"].push_back( nuc1 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", params, fixture.m_interspec ) );

  // Try to add the same nuclide again - should fail
  BOOST_CHECK_THROW( registry.executeTool( "modify_isotopics_nuclides", params, fixture.m_interspec ), std::runtime_error );

  cout << "Correctly rejected duplicate nuclide" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_Update )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add a nuclide
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_params;
  add_params["action"] = "add";
  add_params["nuclides"] = json::array();

  json nuc;
  nuc["name"] = "Pu239";
  nuc["age"] = "10 years";
  add_params["nuclides"].push_back( nuc );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_params, fixture.m_interspec ) );

  // Update the nuclide's age
  json update_params;
  update_params["action"] = "update";
  update_params["nuclides"] = json::array();

  json update_nuc;
  update_nuc["name"] = "Pu239";
  update_nuc["age"] = "20 years";
  update_nuc["fit_age"] = true;
  update_nuc["age_min"] = "5 years";
  update_nuc["age_max"] = "30 years";
  update_params["nuclides"].push_back( update_nuc );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_nuclides", update_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );

  // Verify the update
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  BOOST_REQUIRE( config["nuclides"].size() == 1 );
  BOOST_CHECK( config["nuclides"][0]["fit_age"].get<bool>() == true );

  cout << "Updated nuclide age and fit settings successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_Remove )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add nuclides
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_params;
  add_params["action"] = "add";
  add_params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  add_params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "Pu240";
  add_params["nuclides"].push_back( nuc2 );

  json nuc3;
  nuc3["name"] = "Pu241";
  add_params["nuclides"].push_back( nuc3 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_params, fixture.m_interspec ) );

  // Remove one nuclide
  json remove_params;
  remove_params["action"] = "remove";
  remove_params["nuclides"] = json::array();

  json rem_nuc;
  rem_nuc["name"] = "Pu240";
  remove_params["nuclides"].push_back( rem_nuc );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_nuclides", remove_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["nuclide_count"].get<int>() == 2 );

  // Verify
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );
  BOOST_CHECK( config["nuclides"].size() == 2 );

  // Make sure Pu240 is gone
  bool found_pu240 = false;
  for( const auto &n : config["nuclides"] )
  {
    if( n["name"].get<string>() == "Pu240" )
      found_pu240 = true;
  }
  BOOST_CHECK( !found_pu240 );

  cout << "Removed nuclide successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_InvalidNuclide )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Try to add an invalid nuclide
  json params;
  params["action"] = "add";
  params["nuclides"] = json::array();

  json nuc;
  nuc["name"] = "NotARealNuclide123";
  params["nuclides"].push_back( nuc );

  BOOST_CHECK_THROW( registry.executeTool( "modify_isotopics_nuclides", params, fixture.m_interspec ), std::runtime_error );

  cout << "Correctly rejected invalid nuclide name" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_Add )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset config
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Add ROIs
  json params;
  params["action"] = "add";
  params["rois"] = json::array();

  json roi1;
  roi1["lower_energy"] = 120.0;
  roi1["upper_energy"] = 200.0;
  roi1["continuum_type"] = "Linear";
  params["rois"].push_back( roi1 );

  json roi2;
  roi2["lower_energy"] = 300.0;
  roi2["upper_energy"] = 450.0;
  roi2["continuum_type"] = "Quadratic";
  params["rois"].push_back( roi2 );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_rois", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["roi_count"].get<int>() == 2 );

  // Verify
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );
  BOOST_CHECK( config["rois"].size() == 2 );

  cout << "Added 2 ROIs successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_OverlapRejected )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset config
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Add first ROI
  json params1;
  params1["action"] = "add";
  params1["rois"] = json::array();

  json roi1;
  roi1["lower_energy"] = 100.0;
  roi1["upper_energy"] = 200.0;
  params1["rois"].push_back( roi1 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_rois", params1, fixture.m_interspec ) );

  // Try to add overlapping ROI
  json params2;
  params2["action"] = "add";
  params2["rois"] = json::array();

  json roi2;
  roi2["lower_energy"] = 150.0;  // Overlaps with first ROI
  roi2["upper_energy"] = 250.0;
  params2["rois"].push_back( roi2 );

  BOOST_CHECK_THROW( registry.executeTool( "modify_isotopics_rois", params2, fixture.m_interspec ), std::runtime_error );

  cout << "Correctly rejected overlapping ROI" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_Update )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add ROI
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_params;
  add_params["action"] = "add";
  add_params["rois"] = json::array();

  json roi;
  roi["lower_energy"] = 100.0;
  roi["upper_energy"] = 200.0;
  roi["continuum_type"] = "Linear";
  add_params["rois"].push_back( roi );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_rois", add_params, fixture.m_interspec ) );

  // Update the ROI - change bounds and continuum type
  json update_params;
  update_params["action"] = "update";
  update_params["rois"] = json::array();

  json update_roi;
  update_roi["lower_energy"] = 100.0;
  update_roi["upper_energy"] = 200.0;
  update_roi["new_lower_energy"] = 95.0;
  update_roi["new_upper_energy"] = 210.0;
  update_roi["continuum_type"] = "Quadratic";
  update_params["rois"].push_back( update_roi );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_rois", update_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );

  // Verify the update
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  BOOST_REQUIRE( config["rois"].size() == 1 );
  BOOST_CHECK( config["rois"][0]["lower_energy"].get<double>() == 95.0 );
  BOOST_CHECK( config["rois"][0]["upper_energy"].get<double>() == 210.0 );

  cout << "Updated ROI bounds successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_Remove )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add ROIs
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_params;
  add_params["action"] = "add";
  add_params["rois"] = json::array();

  json roi1;
  roi1["lower_energy"] = 100.0;
  roi1["upper_energy"] = 200.0;
  add_params["rois"].push_back( roi1 );

  json roi2;
  roi2["lower_energy"] = 300.0;
  roi2["upper_energy"] = 400.0;
  add_params["rois"].push_back( roi2 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_rois", add_params, fixture.m_interspec ) );

  // Remove one ROI
  json remove_params;
  remove_params["action"] = "remove";
  remove_params["rois"] = json::array();

  json rem_roi;
  rem_roi["lower_energy"] = 100.0;
  rem_roi["upper_energy"] = 200.0;
  remove_params["rois"].push_back( rem_roi );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_rois", remove_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["roi_count"].get<int>() == 1 );

  cout << "Removed ROI successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_ReplaceAll )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add ROIs
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_params;
  add_params["action"] = "add";
  add_params["rois"] = json::array();

  json roi1;
  roi1["lower_energy"] = 100.0;
  roi1["upper_energy"] = 200.0;
  add_params["rois"].push_back( roi1 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_rois", add_params, fixture.m_interspec ) );

  // Replace all ROIs with new set
  json replace_params;
  replace_params["action"] = "replace_all";
  replace_params["rois"] = json::array();

  json new_roi1;
  new_roi1["lower_energy"] = 50.0;
  new_roi1["upper_energy"] = 100.0;
  replace_params["rois"].push_back( new_roi1 );

  json new_roi2;
  new_roi2["lower_energy"] = 150.0;
  new_roi2["upper_energy"] = 250.0;
  replace_params["rois"].push_back( new_roi2 );

  json new_roi3;
  new_roi3["lower_energy"] = 300.0;
  new_roi3["upper_energy"] = 400.0;
  replace_params["rois"].push_back( new_roi3 );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_rois", replace_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["roi_count"].get<int>() == 3 );

  cout << "Replaced all ROIs successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsCurveSettings_EqnType )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset first
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;
    preset_name = preset["name"].get<string>();
    break;
  }

  if( preset_name.empty() )
  {
    cout << "Skipping test - no presets available" << endl;
    return;
  }

  json load_params;
  load_params["preset"] = preset_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );

  // Change equation type
  json params;
  params["rel_eff_eqn_type"] = "LnXLnY";
  params["rel_eff_eqn_order"] = 4;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_curve_settings", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["rel_eff_eqn_type"].get<string>() == "LnXLnY" );

  cout << "Changed equation type successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsCurveSettings_PhysicalModel )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset that already uses physical model to avoid complex state transitions
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  // Find a Pu preset - these typically use physical model
  string pu_preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;

    const string name = preset["name"].get<string>();
    if( SpecUtils::icontains( name, "Pu" ) || SpecUtils::icontains( name, "plutonium" ) )
    {
      pu_preset_name = name;
      break;
    }
  }

  if( pu_preset_name.empty() )
  {
    cout << "Skipping test - no Pu preset available" << endl;
    return;
  }

  json load_params;
  load_params["preset"] = pu_preset_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );

  // Get current config to verify it uses physical model
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  bool already_physical = false;
  if( config.contains( "rel_eff_curve" ) && config["rel_eff_curve"].contains( "rel_eff_eqn_type" ) )
  {
    already_physical = ( config["rel_eff_curve"]["rel_eff_eqn_type"].get<string>() == "FramPhysicalModel" );
  }

  if( !already_physical )
  {
    cout << "Loaded preset does not use physical model - skipping test" << endl;
    return;
  }

  // Modify physical model settings - just toggle Hoerl
  json params;
  params["use_hoerl"] = !config["rel_eff_curve"].value( "use_hoerl", false );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_curve_settings", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["rel_eff_eqn_type"].get<string>() == "FramPhysicalModel" );

  cout << "Modified physical model settings successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsOptions )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset first
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;
    preset_name = preset["name"].get<string>();
    break;
  }

  if( preset_name.empty() )
  {
    cout << "Skipping test - no presets available" << endl;
    return;
  }

  json load_params;
  load_params["preset"] = preset_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );

  // Modify options
  json params;
  params["fit_energy_cal"] = true;
  params["fwhm_form"] = "Gadras";
  params["skew_type"] = "NoSkew";
  params["note"] = "Test configuration";
  params["description"] = "Modified by unit test";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_options", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["changes"].is_array() );
  BOOST_CHECK( result["changes"].size() >= 3 );

  // Verify by getting config
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  BOOST_CHECK( config["fit_energy_cal"].get<bool>() == true );
  BOOST_CHECK( config["skew_type"].get<string>() == "NoSkew" );

  cout << "Modified options successfully: " << result["changes"].size() << " changes" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsOptions_BackgroundSubtract )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset first
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;
    preset_name = preset["name"].get<string>();
    break;
  }

  if( preset_name.empty() )
  {
    cout << "Skipping test - no presets available" << endl;
    return;
  }

  json load_params;
  load_params["preset"] = preset_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );

  // Try to enable background subtraction (should work since fixture loads a background)
  json params;
  params["background_subtract"] = true;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_options", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );

  cout << "Enabled background subtraction successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsConstraints_ActivityRatio )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add nuclides
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_nuc_params;
  add_nuc_params["action"] = "add";
  add_nuc_params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  add_nuc_params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "Pu240";
  add_nuc_params["nuclides"].push_back( nuc2 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_nuc_params, fixture.m_interspec ) );

  // Add activity ratio constraint
  json params;
  params["action"] = "add";
  params["constraint_type"] = "activity_ratio";
  params["constraints"] = json::array();

  json constraint;
  constraint["nuclide"] = "Pu240";
  constraint["relative_to"] = "Pu239";
  constraint["ratio"] = 0.0612;
  params["constraints"].push_back( constraint );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_constraints", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["modified"].size() == 1 );

  cout << "Added activity ratio constraint successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsConstraints_MassFraction )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add nuclides
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_nuc_params;
  add_nuc_params["action"] = "add";
  add_nuc_params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "U235";
  add_nuc_params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "U238";
  add_nuc_params["nuclides"].push_back( nuc2 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_nuc_params, fixture.m_interspec ) );

  // Add mass fraction constraint
  json params;
  params["action"] = "add";
  params["constraint_type"] = "mass_fraction";
  params["constraints"] = json::array();

  json constraint;
  constraint["nuclide"] = "U235";
  constraint["lower_fraction"] = 0.03;
  constraint["upper_fraction"] = 0.20;
  params["constraints"].push_back( constraint );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_constraints", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["modified"].size() == 1 );

  cout << "Added mass fraction constraint successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsConstraints_Remove )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add nuclides
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_nuc_params;
  add_nuc_params["action"] = "add";
  add_nuc_params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  add_nuc_params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "Pu240";
  add_nuc_params["nuclides"].push_back( nuc2 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_nuc_params, fixture.m_interspec ) );

  // Add constraint
  json add_params;
  add_params["action"] = "add";
  add_params["constraint_type"] = "activity_ratio";
  add_params["constraints"] = json::array();

  json constraint;
  constraint["nuclide"] = "Pu240";
  constraint["relative_to"] = "Pu239";
  constraint["ratio"] = 0.0612;
  add_params["constraints"].push_back( constraint );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_constraints", add_params, fixture.m_interspec ) );

  // Remove constraint
  json remove_params;
  remove_params["action"] = "remove";
  remove_params["constraint_type"] = "activity_ratio";
  remove_params["constraints"] = json::array();

  json rem_constraint;
  rem_constraint["nuclide"] = "Pu240";
  remove_params["constraints"].push_back( rem_constraint );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_constraints", remove_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );

  cout << "Removed constraint successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsConstraints_ClearAll )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset and add nuclides
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  json add_nuc_params;
  add_nuc_params["action"] = "add";
  add_nuc_params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  add_nuc_params["nuclides"].push_back( nuc1 );

  json nuc2;
  nuc2["name"] = "Pu240";
  add_nuc_params["nuclides"].push_back( nuc2 );

  json nuc3;
  nuc3["name"] = "Pu241";
  add_nuc_params["nuclides"].push_back( nuc3 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_nuc_params, fixture.m_interspec ) );

  // Add multiple constraints
  json add_params;
  add_params["action"] = "add";
  add_params["constraint_type"] = "activity_ratio";
  add_params["constraints"] = json::array();

  json constraint1;
  constraint1["nuclide"] = "Pu240";
  constraint1["relative_to"] = "Pu239";
  constraint1["ratio"] = 0.0612;
  add_params["constraints"].push_back( constraint1 );

  json constraint2;
  constraint2["nuclide"] = "Pu241";
  constraint2["relative_to"] = "Pu239";
  constraint2["ratio"] = 0.002;
  add_params["constraints"].push_back( constraint2 );

  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_constraints", add_params, fixture.m_interspec ) );

  // Clear all activity ratio constraints
  json clear_params;
  clear_params["action"] = "clear_all";
  clear_params["constraint_type"] = "activity_ratio";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_constraints", clear_params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["cleared"].get<int>() == 2 );

  cout << "Cleared all constraints successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeResetIsotopicsConfig )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset first
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string preset_name;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;
    preset_name = preset["name"].get<string>();
    break;
  }

  if( !preset_name.empty() )
  {
    json load_params;
    load_params["preset"] = preset_name;
    BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );
  }

  // Reset the config
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );

  // Verify config is cleared
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  BOOST_CHECK( config["has_config"].get<bool>() == false );

  cout << "Reset isotopics config successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsNuclides_WithXray )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset config
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Add x-ray source (element without mass number)
  json params;
  params["action"] = "add";
  params["nuclides"] = json::array();

  json nuc1;
  nuc1["name"] = "Pu239";
  params["nuclides"].push_back( nuc1 );

  json xray;
  xray["name"] = "Pu";  // X-ray source
  params["nuclides"].push_back( xray );

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "modify_isotopics_nuclides", params, fixture.m_interspec ) );

  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["nuclide_count"].get<int>() == 2 );

  cout << "Added nuclide and x-ray source successfully" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsRois_InvalidBounds )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Reset config
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Try to add ROI with invalid bounds (lower > upper)
  json params;
  params["action"] = "add";
  params["rois"] = json::array();

  json roi;
  roi["lower_energy"] = 200.0;
  roi["upper_energy"] = 100.0;  // Invalid: lower > upper
  params["rois"].push_back( roi );

  BOOST_CHECK_THROW( registry.executeTool( "modify_isotopics_rois", params, fixture.m_interspec ), std::runtime_error );

  cout << "Correctly rejected ROI with invalid bounds" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeModifyIsotopicsCurveSettings_InvalidOrderForPhysical )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a preset first - look for one that uses physical model
  json list_result;
  BOOST_REQUIRE_NO_THROW( list_result = registry.executeTool( "list_isotopics_presets", json::object(), fixture.m_interspec ) );

  string phys_model_preset;
  for( const auto &preset : list_result["presets"] )
  {
    if( preset.contains( "is_current" ) && preset["is_current"].get<bool>() )
      continue;

    // Look for Pu presets as they likely use physical model
    const string name = preset["name"].get<string>();
    if( SpecUtils::icontains( name, "Pu" ) || SpecUtils::icontains( name, "plutonium" ) )
    {
      phys_model_preset = name;
      break;
    }
  }

  if( phys_model_preset.empty() )
  {
    cout << "Skipping test - no Pu preset available for physical model test" << endl;
    return;
  }

  json load_params;
  load_params["preset"] = phys_model_preset;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "load_isotopics_preset", load_params, fixture.m_interspec ) );

  // Get current config to check if it's using physical model
  json config;
  BOOST_REQUIRE_NO_THROW( config = registry.executeTool( "get_isotopics_config", json::object(), fixture.m_interspec ) );

  // Check if the loaded preset uses physical model - if not, skip
  bool uses_physical = false;
  if( config.contains( "rel_eff_curve" ) && config["rel_eff_curve"].contains( "rel_eff_eqn_type" ) )
  {
    uses_physical = ( config["rel_eff_curve"]["rel_eff_eqn_type"].get<string>() == "FramPhysicalModel" );
  }

  if( !uses_physical )
  {
    cout << "Loaded preset does not use physical model - skipping test" << endl;
    return;
  }

  // Now try to set order on physical model - should fail
  json order_params;
  order_params["rel_eff_eqn_order"] = 5;

  BOOST_CHECK_THROW( registry.executeTool( "modify_isotopics_curve_settings", order_params, fixture.m_interspec ), std::runtime_error );

  cout << "Correctly rejected equation order for physical model" << endl;
}


BOOST_AUTO_TEST_CASE( test_executeIsotopics_Br82_EndToEnd )
{
  cout << "\n\nStarting test_executeIsotopics_Br82_EndToEnd" << endl;

  InterSpecTestFixture fixture;
  LlmTools::ToolRegistry registry( *LlmConfig::load() );

  // Reset any existing configuration
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "reset_isotopics_config", json::object(), fixture.m_interspec ) );

  // Step 1: Add Br82 nuclide
  json add_nuc_params;
  add_nuc_params["action"] = "add";
  add_nuc_params["nuclides"] = json::array({
    {
      {"name", "Br82"}
    }
  });
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_nuclides", add_nuc_params, fixture.m_interspec ) );
  cout << "Added Br82 nuclide" << endl;

  // Step 2: Add 10 ROIs with specified energy ranges
  json roi_params;
  roi_params["action"] = "replace_all";
  roi_params["rois"] = json::array({
    {{"lower_energy", 214.7}, {"upper_energy", 228.0}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 267.5}, {"upper_energy", 279.8}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 548.1}, {"upper_energy", 562.3}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 612.7}, {"upper_energy", 625.9}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 691.86}, {"upper_energy", 706.0}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 769.6}, {"upper_energy", 783.8}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 820.9}, {"upper_energy", 835.1}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 1037.0}, {"upper_energy", 1051.15}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 1309.25}, {"upper_energy", 1325.4}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}},
    {{"lower_energy", 1466.16}, {"upper_energy", 1480.34}, {"continuum_type", "Linear"}, {"range_type", "Fixed"}}
  });
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_rois", roi_params, fixture.m_interspec ) );
  cout << "Added 10 ROIs" << endl;

  // Step 3: Configure curve settings - empirical LnXLnY order 4
  json curve_params;
  curve_params["rel_eff_eqn_type"] = "LnXLnY";
  curve_params["rel_eff_eqn_order"] = 4;
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_curve_settings", curve_params, fixture.m_interspec ) );
  cout << "Set rel eff curve to LnXLnY order 4" << endl;

  // Step 4: Configure global options
  json options_params;
  options_params["fit_energy_cal"] = true;
  options_params["fwhm_form"] = "Polynomial_4";
  options_params["skew_type"] = "NoSkew";
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "modify_isotopics_options", options_params, fixture.m_interspec ) );
  cout << "Set options: fit_energy_cal=true, fwhm_form=Polynomial_4, skew_type=NoSkew" << endl;

  // Step 5: Perform the isotopics calculation
  json calc_result;
  BOOST_REQUIRE_NO_THROW( calc_result = registry.executeTool( "perform_isotopics_calculation", json::object(), fixture.m_interspec ) );

  // Check that calculation succeeded
  BOOST_REQUIRE( calc_result.contains("success") );
  BOOST_REQUIRE( calc_result["success"].get<bool>() );

  // Extract and verify the results
  cout << "\nIsotopics calculation complete!" << endl;

  // Check quality metrics
  BOOST_REQUIRE( calc_result.contains("quality") );
  const json &quality = calc_result["quality"];

  BOOST_REQUIRE( quality.contains("chi2") );
  BOOST_REQUIRE( quality.contains("dof") );

  const double chi2 = quality["chi2"].get<double>();
  const int dof = quality["dof"].get<int>();
  const double chi2_per_dof = chi2 / dof;

  cout << "Chi2/dof = " << chi2 << " / " << dof << " = " << chi2_per_dof << endl;

  // Verify chi2 and dof are in expected ranges (based on actual run results)
  // Actual values from run: chi2=4093.1, dof=287, chi2/dof=14.26
  BOOST_CHECK_CLOSE( chi2, 4093.1, 5.0 );  // within 5%
  BOOST_CHECK_EQUAL( dof, 287 );
  BOOST_CHECK_CLOSE( chi2_per_dof, 14.26, 5.0 );  // within 5%

  // Find and verify Br82 in the results
  BOOST_REQUIRE( calc_result.contains("isotopics") );
  const json &isotopics = calc_result["isotopics"];
  BOOST_REQUIRE( isotopics.is_array() );
  BOOST_REQUIRE( isotopics.size() == 1 );  // Should only have Br82

  bool found_br82 = false;
  for( const auto &nuc : isotopics )
  {
    BOOST_REQUIRE( nuc.contains("nuclide") );
    if( nuc["nuclide"].get<string>() == "Br82" )
    {
      found_br82 = true;

      BOOST_REQUIRE( nuc.contains("rel_activity") );
      BOOST_REQUIRE( nuc.contains("rel_activity_uncert") );

      const double rel_act = nuc["rel_activity"].get<double>();
      const double rel_act_uncert = nuc["rel_activity_uncert"].get<double>();

      cout << "Br82 relative activity = " << rel_act << " +- " << rel_act_uncert << endl;

      // Verify relative activity is in expected range
      // Actual value from run: 23.85 +- 0.377
      BOOST_CHECK_CLOSE( rel_act, 23.85, 5.0 );  // within 5%
      BOOST_CHECK_CLOSE( rel_act_uncert, 0.377, 10.0 );  // within 10%

      // Check that age was not fit (we didn't enable it)
      BOOST_REQUIRE( nuc.contains("age_was_fit") );
      BOOST_CHECK_EQUAL( nuc["age_was_fit"].get<bool>(), false );

      // Check mass fraction is 1.0 (only nuclide)
      BOOST_REQUIRE( nuc.contains("mass_fraction") );
      BOOST_CHECK_CLOSE( nuc["mass_fraction"].get<double>(), 1.0, 0.1 );

      break;
    }
  }

  BOOST_REQUIRE( found_br82 );

  cout << "Successfully validated Br82 isotopics calculation!" << endl;
}
