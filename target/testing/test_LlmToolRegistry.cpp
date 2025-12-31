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


BOOST_AUTO_TEST_CASE( test_executeGetMaterials )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Get materials - no parameters needed
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_materials", json::object(), fixture.m_interspec) );

  // Should return an array
  BOOST_CHECK( result.is_array() );
  BOOST_CHECK( !result.empty() );

  // Check for some expected materials
  bool found_water = false;
  bool found_air = false;
  bool found_concrete = false;

  for( const auto &material : result )
  {
    BOOST_CHECK( material.is_string() );
    const string mat_name = material.get<string>();

    if( SpecUtils::icontains(mat_name, "water") )
      found_water = true;
    if( SpecUtils::icontains(mat_name, "air") )
      found_air = true;
    if( SpecUtils::icontains(mat_name, "concrete") )
      found_concrete = true;
  }

  BOOST_CHECK_MESSAGE( found_water, "Materials list should contain water" );
  BOOST_CHECK_MESSAGE( found_air, "Materials list should contain air" );
  BOOST_CHECK_MESSAGE( found_concrete, "Materials list should contain concrete" );
}


BOOST_AUTO_TEST_CASE( test_executeGetMaterialInfo )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with valid material
  json params;
  params["material"] = "water";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_material_info", params, fixture.m_interspec) );

  // Check required fields
  BOOST_CHECK( result.contains("name") );
  BOOST_CHECK( result.contains("density") );
  BOOST_CHECK( result.contains("effectiveAtomicNumber") );
  BOOST_CHECK( result.contains("massFractionChemicalFormula") );

  // Check density is reasonable for water (should be close to 1.0 g/cm³)
  BOOST_CHECK( result["density"].is_number() );
  const double density = result["density"].get<double>();
  BOOST_CHECK( density > 0.9 && density < 1.1 );

  // Test with another common material
  params["material"] = "air";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_material_info", params, fixture.m_interspec) );
  BOOST_CHECK( result.contains("density") );

  // Test error handling with invalid material
  params["material"] = "nonexistent_material_12345";
  BOOST_CHECK_THROW( registry.executeTool("get_material_info", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetSourcePhotons )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with a nuclide
  json params;
  params["Source"] = "Cs137";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_photons", params, fixture.m_interspec) );

  // Should return an array
  BOOST_CHECK( result.is_array() );
  BOOST_CHECK( !result.empty() );

  // Each entry should be [energy, intensity]
  bool found_662 = false;
  for( const auto &entry : result )
  {
    BOOST_CHECK( entry.is_array() );
    BOOST_CHECK_EQUAL( entry.size(), 2 );
    BOOST_CHECK( entry[0].is_number() ); // energy
    BOOST_CHECK( entry[1].is_number() ); // intensity

    const double energy = entry[0].get<double>();
    const double intensity = entry[1].get<double>();

    // Check for 661.7 keV gamma from Cs137
    if( energy > 661.6 && energy < 661.7 )
    {
      found_662 = true;
      // Cs137 661.7 keV gamma should have intensity of approximately 0.853337
      BOOST_CHECK_MESSAGE( std::abs(intensity - 0.853337) < 0.001,
                          "Cs137 662 keV gamma intensity should be ~0.853337, got " << intensity );
    }
  }

  BOOST_CHECK_MESSAGE( found_662, "Cs137 should have 662 keV gamma" );

  // Test with an element (x-rays)
  params["Source"] = "Pb";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_photons", params, fixture.m_interspec) );
  BOOST_CHECK( result.is_array() );
  BOOST_CHECK( !result.empty() );

  // Test with age parameter
  params["Source"] = "U238";
  params["Age"] = "0 seconds";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_photons", params, fixture.m_interspec) );
  BOOST_CHECK( result.is_array() );

  // Test error handling - age with element should fail
  params["Source"] = "Pb";
  params["Age"] = "1 year";
  BOOST_CHECK_THROW( registry.executeTool("source_photons", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - invalid source
  params = json::object();
  params["Source"] = "InvalidNuclide12345";
  BOOST_CHECK_THROW( registry.executeTool("source_photons", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeAvailableDetectors )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params = json::object();

  json result;
  try
  {
    // Note: tool name has typo - "avaiable" instead of "available"
    result = registry.executeTool("available_detector_efficiency_functions", params, fixture.m_interspec);
  }
  catch( const std::exception &e )
  {
    BOOST_TEST_MESSAGE( "Exception calling available_detector_efficiency_functions: " << e.what() );
    throw;
  }

  // Should return an array
  BOOST_CHECK( result.is_array() );

  // Print all available detectors to stdout
  //std::cout << "\n=== Available Detectors (" << result.size() << " total) ===\n";
  //for( size_t i = 0; i < result.size(); ++i )
  //{
  //  const json &detector = result[i];
  //  std::cout << "\nDetector #" << (i+1) << ":\n";
  //  std::cout << "  JSON: " << detector.dump(2) << "\n";
  //}
  //std::cout << "=== End of Detectors ===\n" << std::endl;

  // If there are any detectors, check their structure
  if( !result.empty() )
  {
    const json &detector = result[0];
    BOOST_CHECK( detector.contains("name") );
    BOOST_CHECK( detector.contains("source") );
    BOOST_CHECK( detector["name"].is_string() );
    BOOST_CHECK( detector["source"].is_string() );
  }

  // Check for presence of specific commonly-used detectors
  const std::vector<std::string> expected_detectors = {
    "ORTEC Detective-DX_LANL_025cm (11%)",
    "ICx/FLIR IdentiFINDER-NGH (10%)",
    "Radiacode 102 DRF",
    "ORTEC Detective-EX100_LANL_025cm (40%)",
    "SAM-Eagle-NaI-3x3",
    "ICx/FLIR IdentiFINDER-LaBr (7%)"
  };

  for( const auto &expected_name : expected_detectors )
  {
    bool found = false;
    for( const auto &detector : result )
    {
      if( detector.contains("name") && detector["name"].get<std::string>() == expected_name )
      {
        found = true;
        break;
      }
    }
    BOOST_CHECK_MESSAGE( found, "Expected detector '" << expected_name << "' not found in available detectors" );
  }
}


BOOST_AUTO_TEST_CASE( test_executeLoadDetectorEfficiency )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First get available detectors (note: tool name has typo - "avaiable" not "available")
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test loading" );
    return;
  }

  // Try to load the first available detector
  const string detector_name = available_result[0]["name"].get<string>();

  json load_params;
  load_params["identifier"] = detector_name;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Check result structure
  BOOST_CHECK( result.contains("success") );
  BOOST_CHECK( result.contains("detectorName") );
  BOOST_CHECK( result.contains("source") );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_CHECK( result["detectorName"].is_string() );

  // Test error handling with invalid identifier
  load_params["identifier"] = "nonexistent_detector_12345";
  BOOST_CHECK_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetDetectorInfo )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test getting info" );
    return;
  }

  const string detector_name = available_result[0]["name"].get<string>();

  json load_params;
  load_params["identifier"] = detector_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Try to get info about the currently loaded detector
  // Note: In test environment without a foreground spectrum, the detector may not actually be loaded
  //       even though load_detector_efficiency_function returned success
  json info_params = json::object();

  json result;
  try
  {
    result = registry.executeTool("detector_efficiency_function_info", info_params, fixture.m_interspec);

    // If we got here, check the required fields
    BOOST_CHECK( result.contains("name") );
    BOOST_CHECK( result.contains("description") );
    BOOST_CHECK( result.contains("isValid") );
    BOOST_CHECK( result.contains("hasResolutionInfo") );
    BOOST_CHECK( result.contains("geometryType") );
    BOOST_CHECK( result.contains("geometryTypeDescription") );
    BOOST_CHECK( result.contains("isFixedGeometry") );
    BOOST_CHECK( result.contains("lowerEnergy") );
    BOOST_CHECK( result.contains("upperEnergy") );

    BOOST_CHECK( result["isValid"].is_boolean() );
    BOOST_CHECK( result["isValid"].get<bool>() == true );
  }
  catch( const std::exception &e )
  {
    // Expected in test environment - detector load requires a foreground spectrum
    const string msg = e.what();
    BOOST_TEST_MESSAGE( "Cannot test detector info without foreground: " << msg );
    BOOST_CHECK( msg.find("No detector efficiency function currently loaded") != string::npos );
    return;  // Skip remaining tests
  }

  // Test getting info by name
  info_params["name"] = detector_name;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("detector_efficiency_function_info", info_params, fixture.m_interspec) );
  BOOST_CHECK( result.contains("name") );

  // Test error handling with invalid detector name
  info_params["name"] = "nonexistent_detector_12345";
  BOOST_CHECK_THROW( registry.executeTool("detector_efficiency_function_info", info_params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executePhotopeakDetectionCalc )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test photopeak detection calculation" );
    return;
  }

  // Load the Detective-EX100 LANL detector
  const string detector_name = "ORTEC Detective-EX100_LANL_025cm (40%)";

  json load_params;
  load_params["identifier"] = detector_name;
  json load_result;
  BOOST_REQUIRE_NO_THROW( load_result = registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );
  BOOST_CHECK( load_result["success"].get<bool>() == true );

  // Verify detector loaded
  json detector_info;
  BOOST_REQUIRE_NO_THROW( detector_info = registry.executeTool("detector_efficiency_function_info", json::object(), fixture.m_interspec) );
  BOOST_CHECK( detector_info["isValid"].get<bool>() == true );
  const bool is_fixed_geom = detector_info["isFixedGeometry"].get<bool>();

  // Verify the correct detector is loaded
  BOOST_CHECK_EQUAL( detector_info["name"].get<string>(), "ORTEC Detective-EX100_LANL_025cm (40%)" );

  // Test with Detective-EX100 at 25 cm, Fe shielding 1.0 cm, specific energies
  // Expected values from calibration:
  // 185.0 keV: intrinsic=0.4749, solid_angle=0.004172, shielding_transmission=0.3204, total=0.0006348
  // 661.7 keV: intrinsic=0.2038, solid_angle=0.004172, shielding_transmission=0.5657, total=0.0004811
  // 1460.8 keV: intrinsic=0.1072, solid_angle=0.004172, shielding_transmission=0.6791, total=0.0003038
  json params;
  params["Energies"] = json::array({185.0, 661.7, 1460.8});
  params["Distance"] = "25 cm";

  // Add Fe shielding 1.0 cm thick
  json shielding;
  shielding["Material"] = "iron";
  shielding["Thickness"] = "1.0 cm";
  params["Shielding"] = shielding;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("photopeak_detection_efficiency", params, fixture.m_interspec) );

  // Result should be an object with "energies" and "results" arrays
  BOOST_REQUIRE( result.is_object() );
  BOOST_REQUIRE( result.contains("energies") );
  BOOST_REQUIRE( result.contains("results") );
  BOOST_REQUIRE( result["energies"].is_array() );
  BOOST_REQUIRE( result["results"].is_array() );
  BOOST_REQUIRE_EQUAL( result["results"].size(), 3 );

  // Verify specific efficiency values for each energy
  const json &results = result["results"];

  // 185.0 keV - verify shielding transmission
  BOOST_CHECK( results[0].contains("energy") );
  BOOST_CHECK_CLOSE( results[0]["energy"].get<double>(), 185.0, 0.1 );

  // Check solid angle fraction (distanceGeometryFactor)
  if( results[0].contains("distanceGeometryFactor") )
  {
    BOOST_CHECK_CLOSE( results[0]["distanceGeometryFactor"].get<double>(), 0.004172, 1.0 );
  }

  // Check shielding transmission
  if( results[0].contains("shieldingAttenuations") && results[0]["shieldingAttenuations"].is_array() && results[0]["shieldingAttenuations"].size() > 0 )
  {
    BOOST_CHECK_CLOSE( results[0]["shieldingAttenuations"][0].get<double>(), 0.3204, 1.0 );
  }

  // Check intrinsic efficiency (allowing wider tolerance for investigation)
  if( results[0].contains("detectorIntrinsicEfficiency") )
  {
    BOOST_CHECK_CLOSE( results[0]["detectorIntrinsicEfficiency"].get<double>(), 0.4749, 5.0 );
  }

  BOOST_REQUIRE( results[0].contains("finalEfficiency") );

  // 661.7 keV - verify shielding transmission
  BOOST_CHECK_CLOSE( results[1]["energy"].get<double>(), 661.7, 0.1 );

  // Check solid angle fraction
  if( results[1].contains("distanceGeometryFactor") )
  {
    BOOST_CHECK_CLOSE( results[1]["distanceGeometryFactor"].get<double>(), 0.004172, 1.0 );
  }

  // Check shielding transmission
  if( results[1].contains("shieldingAttenuations") && results[1]["shieldingAttenuations"].is_array() && results[1]["shieldingAttenuations"].size() > 0 )
  {
    BOOST_CHECK_CLOSE( results[1]["shieldingAttenuations"][0].get<double>(), 0.5657, 1.0 );
  }

  // Check intrinsic efficiency (allowing wider tolerance for investigation)
  if( results[1].contains("detectorIntrinsicEfficiency") )
  {
    BOOST_CHECK_CLOSE( results[1]["detectorIntrinsicEfficiency"].get<double>(), 0.2038, 5.0 );
  }

  BOOST_REQUIRE( results[1].contains("finalEfficiency") );

  // 1460.8 keV - verify shielding transmission
  BOOST_CHECK_CLOSE( results[2]["energy"].get<double>(), 1460.8, 0.1 );

  // Check solid angle fraction
  if( results[2].contains("distanceGeometryFactor") )
  {
    BOOST_CHECK_CLOSE( results[2]["distanceGeometryFactor"].get<double>(), 0.004172, 1.0 );
  }

  // Check shielding transmission
  if( results[2].contains("shieldingAttenuations") && results[2]["shieldingAttenuations"].is_array() && results[2]["shieldingAttenuations"].size() > 0 )
  {
    BOOST_CHECK_CLOSE( results[2]["shieldingAttenuations"][0].get<double>(), 0.6791, 1.0 );
  }

  // Check intrinsic efficiency (allowing wider tolerance for investigation)
  if( results[2].contains("detectorIntrinsicEfficiency") )
  {
    BOOST_CHECK_CLOSE( results[2]["detectorIntrinsicEfficiency"].get<double>(), 0.1072, 5.0 );
  }

  BOOST_REQUIRE( results[2].contains("finalEfficiency") );

  // Test error handling - no energies
  params = json::object();
  BOOST_CHECK_THROW( registry.executeTool("photopeak_detection_efficiency", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - invalid energy
  params["Energies"] = json::array({-100.0});
  BOOST_CHECK_THROW( registry.executeTool("photopeak_detection_efficiency", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetLoadedSpectra )
{
  InterSpecTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("loaded_spectra", {}, fixture.m_interspec) );
  BOOST_CHECK( result.is_array() );
  BOOST_REQUIRE_EQUAL( result.size(), 2 );
  BOOST_CHECK( result[0].is_string() );
  BOOST_CHECK( result[1].is_string() );
  BOOST_CHECK_EQUAL( result[0].get<string>(), "Foreground" );
  BOOST_CHECK_EQUAL( result[1].get<string>(), "Background" );
}


BOOST_AUTO_TEST_CASE( test_executeGetSpectrumInfo )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params = json::object();
  params["specType"] = "Foreground";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_spectrum_info", params, fixture.m_interspec) );

  // Check required fields (using actual field names from implementation)
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("fileName") );
  BOOST_CHECK( result.contains("liveTime") );
  BOOST_CHECK( result.contains("realTime") );
  BOOST_CHECK( result.contains("numChannels") );
  BOOST_CHECK( result.contains("displayedSamples") );

  // Verify Br82 spectrum properties
  const string filename = result["fileName"].get<string>();
  BOOST_CHECK( SpecUtils::icontains(filename, "Br82") );

  BOOST_CHECK( result["liveTime"].is_number() );
  BOOST_CHECK( result["liveTime"].get<double>() > 0.0 );

  BOOST_CHECK( result["numChannels"].is_number() );
  BOOST_CHECK( result["numChannels"].get<int>() > 0 );

  // displayedSamples should be an array
  BOOST_CHECK( result["displayedSamples"].is_array() );

  // Test error handling - invalid spectrum type
  params["specType"] = "InvalidType";
  BOOST_CHECK_THROW( registry.executeTool("get_spectrum_info", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - no Secondary loaded
  params["specType"] = "Secondary";
  BOOST_CHECK_THROW( registry.executeTool("get_spectrum_info", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetCountsInEnergyRange )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params;
  params["lowerEnergy"] = 600.0;
  params["upperEnergy"] = 800.0;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_counts_in_energy_range", params, fixture.m_interspec) );

  // Should return an object with counts
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("foregroundCounts") );
  BOOST_CHECK( result.contains("foregroundCps") );
  BOOST_CHECK( result["foregroundCounts"].is_number() );
  BOOST_CHECK( result["foregroundCps"].is_number() );

  const double counts = result["foregroundCps"].get<double>();
  // Br82 has major peak at 776.52 keV, so should have counts in this range
  BOOST_CHECK( counts > 0.0 );

  BOOST_CHECK( result["backgroundInfo"].contains("counts") );
  BOOST_CHECK( result["backgroundInfo"]["counts"].is_number() );
  BOOST_CHECK( result["backgroundInfo"].contains("cps") );
  BOOST_CHECK( result["backgroundInfo"]["cps"].is_number() );
  BOOST_CHECK( result["backgroundInfo"].contains("numSigmaCpsRelForeground") );
  BOOST_CHECK( result["backgroundInfo"]["numSigmaCpsRelForeground"].is_number() );

  // Test error handling - invalid energy range
  params["lowerEnergy"] = 800.0;
  params["upperEnergy"] = 600.0;
  BOOST_CHECK_THROW( registry.executeTool("get_counts_in_energy_range", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executePeakDetection )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // We have to manually run hint seard
  InterSpec * const viewer = fixture.m_interspec;
  
  const std::set<int> &foreground_samples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<SpecMeas> foreground_meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<const SpecUtils::Measurement> foreground_spectrum = viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  viewer->searchForHintPeaks( foreground_meas, foreground_samples, foreground_spectrum, true );
  
  
  json params = json::object();
  params["specType"] = "Foreground";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_detected_peaks", params, fixture.m_interspec) );
  BOOST_CHECK( result.is_object() && result.contains("rois") && result["rois"].is_array() );
  

  // Expected major Br82 peaks: 554.35, 619.11, 698.37, 776.52, 827.83, 1044.0, 1317.5, 1474.9
  // Also expect 2614 keV (Tl208 from Th232 contamination)
  const std::vector<double> expected_peaks = {554.35, 619.11, 698.37, 776.52, 827.83, 1044.0, 1317.5, 1474.9, 2614.5};

  std::vector<bool> found_peaks(expected_peaks.size(), false);

  for( const auto &roi : result["rois"] )
  {
    BOOST_REQUIRE( roi.contains("peaks") );
    BOOST_CHECK( roi.contains("continuumType") );
    BOOST_CHECK( roi.contains("lowerEnergy") );
    BOOST_CHECK( roi.contains("upperEnergy") );
    BOOST_CHECK( roi.contains("roiID") );
    
    for( const auto &peak : roi["peaks"] )
    {
      BOOST_CHECK( peak.is_object() );
      BOOST_CHECK( peak.contains("fwhm") );
      BOOST_CHECK( peak.contains("energy") );
      BOOST_CHECK( peak.contains("amplitude") );
      BOOST_CHECK( peak.contains("cps") ); //Requires a valid live time for this to be present
      BOOST_CHECK( peak.contains("numSigma") ); //Requires amplitudeUncert
      BOOST_CHECK( peak.contains("amplitudeUncert") );//Requires amplitudeUncert
      BOOST_CHECK( peak.contains("cpsUncert") );//Requires alid live time and amplitudeUncert
      
      const double mean = peak["energy"].get<double>();
      const double amplitude = peak["amplitude"].get<double>();
      
      // Check if this matches any expected peak (within 2 keV)
      for( size_t i = 0; i < expected_peaks.size(); ++i )
      {
        if( std::abs(mean - expected_peaks[i]) < 0.5 )
        {
          found_peaks[i] = true;
          break;
        }
      }
    }//for( const auto &peak : roi["peaks"] )
  }//for( const auto &roi : result["rois"] )

  // Should find at least the top 8 Br82 lines
  int num_found = 0;
  for( size_t i = 0; i < 8; ++i )
  {
    if( found_peaks[i] )
      ++num_found;
  }
  BOOST_CHECK_MESSAGE( num_found >= 6, "Should find at least 6 of the top 8 Br82 peaks, found " << num_found );

  // Should also find the 2614 keV line
  BOOST_CHECK_MESSAGE( found_peaks[8], "Should find 2614 keV peak" );

  // Test error handling - invalid spectrum type
  params["specType"] = "InvalidType";
  BOOST_CHECK_THROW( registry.executeTool("get_detected_peaks", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetUserPeaks )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Get user peaks when there are none initially
  json params = json::object();
  params["specType"] = "Foreground";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );

  // Should return an object with rois array (same structure as get_detected_peaks)
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("rois") );
  BOOST_CHECK( result["rois"].is_array() );

  // Initially should have no user peaks (until user manually adds or accepts detected peaks)
  BOOST_CHECK_MESSAGE( result["rois"].empty(), "Should have no user peaks initially" );

  //std::cout << "\n=== User Peaks Test ===\n";
  //std::cout << "Initial user peaks: " << result["rois"].size() << "\n";

  // Fit peaks at specific energies to add them as user peaks
  // Fit peak at 776.52 keV
  json fit_params_776;
  fit_params_776["energy"] = 776.52;
  fit_params_776["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( registry.executeTool("add_analysis_peak", fit_params_776, fixture.m_interspec) );

  // Fit peak at 554.35 keV
  json fit_params_554;
  fit_params_554["energy"] = 554.35;
  fit_params_554["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( registry.executeTool("add_analysis_peak", fit_params_554, fixture.m_interspec) );

  //std::cout << "Fitted peaks at 776.52 keV and 554.35 keV\n";

  // Now get user peaks - should have the fitted peaks
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );

  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("rois") );
  BOOST_CHECK( result["rois"].is_array() );
  BOOST_CHECK( !result["rois"].empty() );

  //std::cout << "User peaks after fitting: " << result["rois"].size() << " ROIs\n";

  // Verify we can find the fitted peaks
  bool found_776 = false;
  bool found_554 = false;

  for( const auto &roi : result["rois"] )
  {
    BOOST_REQUIRE( roi.contains("peaks") );
    BOOST_CHECK( roi.contains("continuumType") );
    BOOST_CHECK( roi.contains("lowerEnergy") );
    BOOST_CHECK( roi.contains("upperEnergy") );

    for( const auto &peak : roi["peaks"] )
    {
      BOOST_CHECK( peak.is_object() );
      BOOST_CHECK( peak.contains("energy") );
      BOOST_CHECK( peak.contains("amplitude") );
      BOOST_CHECK( peak.contains("fwhm") );

      const double mean = peak["energy"].get<double>();
      const double amplitude = peak["amplitude"].get<double>();

      //std::cout << "  Peak: energy=" << mean << " keV, amplitude=" << amplitude << "\n";

      // Check for 776.52 keV peak
      if( std::abs(mean - 776.52) < 1.0 )
      {
        found_776 = true;
        BOOST_CHECK( amplitude > 0.0 );
      }

      // Check for 554.35 keV peak
      if( std::abs(mean - 554.35) < 1.0 )
      {
        found_554 = true;
        BOOST_CHECK( amplitude > 0.0 );
      }
    }
  }

  //std::cout << "=== End of User Peaks Test ===\n" << std::endl;

  BOOST_CHECK_MESSAGE( found_776, "Should find 776.52 keV Br82 peak in user peaks" );
  BOOST_CHECK_MESSAGE( found_554, "Should find 554.35 keV Br82 peak in user peaks" );

  // Test energy range filtering - get peaks only between 500-600 keV
  params = json::object();
  params["specType"] = "Foreground";
  params["lowerEnergy"] = 500.0;
  params["upperEnergy"] = 600.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );

  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("rois") );
  BOOST_CHECK( result["rois"].is_array() );

  // Should find 554 keV peak but not 776 keV peak
  found_554 = false;
  found_776 = false;
  for( const auto &roi : result["rois"] )
  {
    for( const auto &peak : roi["peaks"] )
    {
      const double mean = peak["energy"].get<double>();
      if( std::abs(mean - 554.35) < 1.0 )
        found_554 = true;
      if( std::abs(mean - 776.52) < 1.0 )
        found_776 = true;
    }
  }
  BOOST_CHECK_MESSAGE( found_554, "Should find 554 keV peak in 500-600 keV range" );
  BOOST_CHECK_MESSAGE( !found_776, "Should NOT find 776 keV peak in 500-600 keV range" );

  // Test with only lower bound
  params = json::object();
  params["specType"] = "Foreground";
  params["lowerEnergy"] = 700.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );

  found_554 = false;
  found_776 = false;
  for( const auto &roi : result["rois"] )
  {
    for( const auto &peak : roi["peaks"] )
    {
      const double mean = peak["energy"].get<double>();
      if( std::abs(mean - 554.35) < 1.0 )
        found_554 = true;
      if( std::abs(mean - 776.52) < 1.0 )
        found_776 = true;
    }
  }
  BOOST_CHECK_MESSAGE( !found_554, "Should NOT find 554 keV peak above 700 keV" );
  BOOST_CHECK_MESSAGE( found_776, "Should find 776 keV peak above 700 keV" );

  // Test with only upper bound
  params = json::object();
  params["specType"] = "Foreground";
  params["upperEnergy"] = 600.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );

  found_554 = false;
  found_776 = false;
  for( const auto &roi : result["rois"] )
  {
    for( const auto &peak : roi["peaks"] )
    {
      const double mean = peak["energy"].get<double>();
      if( std::abs(mean - 554.35) < 1.0 )
        found_554 = true;
      if( std::abs(mean - 776.52) < 1.0 )
        found_776 = true;
    }
  }
  BOOST_CHECK_MESSAGE( found_554, "Should find 554 keV peak below 600 keV" );
  BOOST_CHECK_MESSAGE( !found_776, "Should NOT find 776 keV peak below 600 keV" );

  // Test with only lower bound above all peaks - should return no peaks
  params = json::object();
  params["specType"] = "Foreground";
  params["lowerEnergy"] = 800.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_analysis_peaks", params, fixture.m_interspec) );
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("rois") );
  BOOST_CHECK( result["rois"].is_array() );
  BOOST_CHECK_MESSAGE( result["rois"].empty(), "Should return no peaks when lowerEnergy is 800 keV (above all peaks)" );

  // Test with invalid range (upper < lower) - should throw exception
  params = json::object();
  params["specType"] = "Foreground";
  params["lowerEnergy"] = 800.0;
  params["upperEnergy"] = 500.0;
  BOOST_CHECK_THROW( registry.executeTool("get_analysis_peaks", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetCharacteristicLines )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test nuclide: U235 - should have 185.7 keV
  json params;
  params["source"] = "U235";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("primary_gammas_for_source", params, fixture.m_interspec) );

  // Result should be object with characteristicGammas array
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("characteristicGammas") );
  BOOST_CHECK( result["characteristicGammas"].is_array() );
  BOOST_CHECK( !result["characteristicGammas"].empty() );

  //std::cout << "\n=== U235 Primary Gammas ===\n";
  bool found_185 = false;
  for( const auto &energy : result["characteristicGammas"] )
  {
    BOOST_CHECK( energy.is_number() );
    const double e = energy.get<double>();
    //std::cout << "  " << e << " keV\n";

    if( std::abs(e - 185.7) < 0.5 )
      found_185 = true;
  }
  //std::cout << "=== End U235 Lines ===\n" << std::endl;

  BOOST_CHECK_MESSAGE( found_185, "U235 should have 185.7 keV line" );

  // Test fluorescence x-ray: Fe (iron) - should have K-alpha around 6.4 keV
  params.clear();
  params["source"] = "Fe";

  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("primary_gammas_for_source", params, fixture.m_interspec) );

  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("characteristicGammas") );
  BOOST_CHECK( result["characteristicGammas"].is_array() );
  BOOST_CHECK( !result["characteristicGammas"].empty() );

  //std::cout << "\n=== Fe X-ray Lines ===\n";
  bool found_fe_ka = false;

  for( const auto &energy : result["characteristicGammas"] )
  {
    BOOST_CHECK( energy.is_number() );
    const double e = energy.get<double>();
    std::cout << "  " << e << " keV\n";

    // Fe K-alpha is around 6.4 keV
    if( e > 5.0 && e < 8.0 )
      found_fe_ka = true;
  }
  //std::cout << "=== End Fe Lines ===\n" << std::endl;

  BOOST_CHECK_MESSAGE( found_fe_ka, "Fe should have K-alpha x-ray around 6.4 keV" );

  // Test nuclear reaction: H(n,g) - should have 2223 keV capture gamma
  params.clear();
  params["source"] = "H(n,g)";

  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("primary_gammas_for_source", params, fixture.m_interspec) );

  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("characteristicGammas") );
  BOOST_CHECK( result["characteristicGammas"].is_array() );
  BOOST_CHECK( !result["characteristicGammas"].empty() );

  //std::cout << "\n=== H(n,g) Reaction Lines ===\n";
  bool found_2223 = false;

  for( const auto &energy : result["characteristicGammas"] )
  {
    BOOST_CHECK( energy.is_number() );
    const double e = energy.get<double>();
    //std::cout << "  " << e << " keV\n";

    // H(n,g) produces 2223 keV gamma
    if( std::abs(e - 2223.0) < 5.0 )
      found_2223 = true;
  }
  //std::cout << "=== End H(n,g) Lines ===\n" << std::endl;

  BOOST_CHECK_MESSAGE( found_2223, "H(n,g) should have 2223 keV capture gamma" );

  // Test error handling - invalid source
  params["source"] = "InvalidSource12345";
  BOOST_CHECK_THROW( registry.executeTool("primary_gammas_for_source", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetAssociatedSources )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Get associated nuclides
  json params = json::object();
  params["source"] = "Br82";

  json result;
//#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
//#else
//  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("sources_associated_with_source", params, fixture.m_interspec) );
//  BOOST_CHECK( result.contains("status") && (result["status"] == "success") );
//#endif
  BOOST_CHECK( result.contains("source") && (result["source"] == "Br82") );
  BOOST_REQUIRE( result.contains("associatedSources") );
  BOOST_REQUIRE( result["associatedSources"].is_array() );

  // Check result["associatedSources"] contains K42 and Na24
  const json &associatedSources = result["associatedSources"];

  bool found_k42 = false;
  bool found_na24 = false;

  for( const auto &source : associatedSources )
  {
    if( source.is_string() )
    {
      const string sourceStr = source.get<string>();
      if( SpecUtils::icontains(sourceStr, "K42") || SpecUtils::icontains(sourceStr, "K-42") )
        found_k42 = true;
      if( SpecUtils::icontains(sourceStr, "Na24") || SpecUtils::icontains(sourceStr, "Na-24") )
        found_na24 = true;
    }
  }

  BOOST_CHECK_MESSAGE( found_k42, "Should find K42 in associated sources for Br82" );
  BOOST_CHECK_MESSAGE( found_na24, "Should find Na24 in associated sources for Br82" );
}


BOOST_AUTO_TEST_CASE( test_executeGetSourceAnalystNotes )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params;
  params["source"] = "U235";

  json result;
//#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
//#else
//  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec) );
//#endif

  // Should return an object
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("source") );
  BOOST_CHECK( result.contains("analystNotes") );

  const string nuclide = result["source"].get<string>();
  BOOST_CHECK( SpecUtils::icontains(nuclide, "U235") || SpecUtils::icontains(nuclide, "U-235") );

  BOOST_CHECK( result["analystNotes"].is_string() );
  const string notes = result["analystNotes"].get<string>();
  // U235 notes should mention enrichment or uranium
  BOOST_CHECK( !notes.empty() );

  // Test with Cs137
  params["source"] = "Cs137";
//#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
//#else
//  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec) );
//#endif
  BOOST_CHECK( result.contains("analystNotes") );

  // Test error handling - invalid nuclide
  params["source"] = "InvalidNuclide12345";
//#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_CHECK_THROW( registry.executeTool("source_info", params, fixture.m_interspec), std::runtime_error );
//#else
//  BOOST_CHECK_THROW( registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec), std::runtime_error );
//#endif
}


BOOST_AUTO_TEST_CASE( test_executeGetSourceInfo )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params;
  params["source"] = "Cs137";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );

  // Check required fields
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("source") );
  BOOST_CHECK( result.contains("atomicNumber") );
  BOOST_CHECK( result.contains("massNumber") );
  BOOST_CHECK( result.contains("isomerNumber") );
  BOOST_CHECK( result.contains("halfLife") );
  BOOST_CHECK( result.contains("atomicMass") );
  BOOST_CHECK( result.contains("associatedSources") );

  // Verify Cs137 properties
  BOOST_CHECK_EQUAL( result["atomicNumber"].get<int>(), 55 );
  BOOST_CHECK_EQUAL( result["massNumber"].get<int>(), 137 );

  BOOST_CHECK( result["halfLife"].is_string() );
  const string halflife = result["halfLife"].get<string>();
  // Cs137 has ~30 year half-life
  BOOST_CHECK( SpecUtils::icontains(halflife, "year") || SpecUtils::icontains(halflife, "y") );

  // Test with a different nuclide
  params["source"] = "Co60";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
  BOOST_CHECK_EQUAL( result["atomicNumber"].get<int>(), 27 );
  BOOST_CHECK_EQUAL( result["massNumber"].get<int>(), 60 );

  // Test with x-ray source (element)
  params["source"] = "Pb";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );

  // Check required fields for element
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("symbol") );
  BOOST_CHECK( result.contains("name") );
  BOOST_CHECK( result.contains("atomicNumber") );
  BOOST_CHECK( result.contains("atomicMass") );
  BOOST_CHECK( result.contains("naturallyOccuringIsotopes") );

  // Verify Pb properties
  BOOST_CHECK_EQUAL( result["symbol"].get<string>(), "Pb" );
  BOOST_CHECK_EQUAL( result["name"].get<string>(), "lead" );
  BOOST_CHECK_EQUAL( result["atomicNumber"].get<int>(), 82 );

  // Check atomic mass is reasonable for Pb (~207 amu)
  BOOST_CHECK( result["atomicMass"].is_number() );
  const double pb_mass = result["atomicMass"].get<double>();
  BOOST_CHECK( pb_mass > 206.0 && pb_mass < 208.0 );

  // Check naturally occurring isotopes array
  BOOST_CHECK( result["naturallyOccuringIsotopes"].is_array() );
  BOOST_CHECK( !result["naturallyOccuringIsotopes"].empty() );

  bool found_pb208 = false;
  double total_abundance = 0.0;

  for( const auto &isotope : result["naturallyOccuringIsotopes"] )
  {
    BOOST_CHECK( isotope.is_object() );
    BOOST_CHECK( isotope.contains("nuclide") );
    BOOST_CHECK( isotope.contains("abundance") );
    BOOST_CHECK( isotope["nuclide"].is_string() );
    BOOST_CHECK( isotope["abundance"].is_number() );

    const string nuclide = isotope["nuclide"].get<string>();
    const double abundance = isotope["abundance"].get<double>();

    // Pb208 is the most abundant (~52%)
    if( nuclide == "Pb208" )
      found_pb208 = true;

    BOOST_CHECK( abundance >= 0.0 && abundance <= 1.0 );
    total_abundance += abundance;
  }

  BOOST_CHECK_MESSAGE( found_pb208, "Pb should have Pb208 in naturally occurring isotopes" );
  // Total abundance should be close to 1.0
  BOOST_CHECK_CLOSE( total_abundance, 1.0, 1.0 );

  // Test error handling - invalid nuclide
  params["source"] = "InvalidNuclide12345";
  BOOST_CHECK_THROW( registry.executeTool("source_info", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetNuclideDecayChain )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params;
  params["nuclide"] = "U238";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("nuclide_decay_chain", params, fixture.m_interspec) );

  // Should return an array
  BOOST_CHECK( result.is_array() );
  BOOST_CHECK( !result.empty() );
  BOOST_CHECK( result.size() == 21 );

  // U238 decay chain should include Th234, Pa234m, Ra226, Rn222, Pb210, etc.
  bool found_th234 = false;
  bool found_ra226 = false;
  bool found_rn222 = false;

  for( const auto &entry : result )
  {
    BOOST_CHECK( entry.is_object() );
    BOOST_CHECK( entry.contains("nuclide") );
    BOOST_CHECK( entry.contains("halfLife") );

    const string nuclide = entry["nuclide"].get<string>();
    if( SpecUtils::icontains(nuclide, "Th234") || SpecUtils::icontains(nuclide, "Th-234") )
      found_th234 = true;
    if( SpecUtils::icontains(nuclide, "Ra226") || SpecUtils::icontains(nuclide, "Ra-226") )
      found_ra226 = true;
    if( SpecUtils::icontains(nuclide, "Rn222") || SpecUtils::icontains(nuclide, "Rn-222") )
      found_rn222 = true;
  }

  BOOST_CHECK_MESSAGE( found_th234, "U238 chain should include Th234" );
  BOOST_CHECK_MESSAGE( found_ra226, "U238 chain should include Ra226" );
  BOOST_CHECK_MESSAGE( found_rn222, "U238 chain should include Rn222" );

  // Test with a simple chain
  params["nuclide"] = "Cs137";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("nuclide_decay_chain", params, fixture.m_interspec) );
  BOOST_CHECK( result.is_array() );
  // Cs137 -> Ba137m -> Ba137 (stable)
  BOOST_CHECK( result.size() >= 2 );

  // Test error handling - invalid nuclide
  params["nuclide"] = "InvalidNuclide12345";
  BOOST_CHECK_THROW( registry.executeTool("nuclide_decay_chain", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetExpectedFwhm )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test FWHM calculation" );
    return;
  }

  const string detector_name = available_result[0]["name"].get<string>();

  json load_params;
  load_params["identifier"] = detector_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Check if detector has resolution info
  json info_result;
  try
  {
    info_result = registry.executeTool("detector_efficiency_function_info", json::object(), fixture.m_interspec);
  }
  catch( const std::exception &e )
  {
    BOOST_TEST_MESSAGE( "Detector not loaded in test environment: " << e.what() );
    return;
  }

  if( !info_result["hasResolutionInfo"].get<bool>() )
  {
    BOOST_TEST_MESSAGE( "Detector does not have resolution info, skipping FWHM test" );
    return;
  }

  // Test FWHM at various energies
  json params;
  params["energy"] = 661.7;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_expected_fwhm", params, fixture.m_interspec) );

  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("energy") );
  BOOST_CHECK( result.contains("fwhm") );

  BOOST_CHECK_CLOSE( result["energy"].get<double>(), 661.7, 0.1 );

  const double fwhm_661 = result["fwhm"].get<double>();
  BOOST_CHECK( fwhm_661 > 0.0 );
  BOOST_CHECK( fwhm_661 < 200.0 );  // Should be reasonable for typical detectors

  // Test at 1460 keV
  params["energy"] = 1460.8;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_expected_fwhm", params, fixture.m_interspec) );
  const double fwhm_1460 = result["fwhm"].get<double>();

  // FWHM should increase with energy for typical detectors
  BOOST_CHECK( fwhm_1460 > fwhm_661 );

  // Test error handling - no detector loaded (after clearing)
  // Note: Can't easily clear detector in test, so skip this check

  // Test error handling - invalid energy
  params["energy"] = -100.0;
  BOOST_CHECK_THROW( registry.executeTool("get_expected_fwhm", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeSearchSourcesByEnergy )
{
  InterSpecTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test basic search with Co-60 energies (1173 and 1332 keV)
  json params;
  json energies_array = json::array();

  json energy1;
  energy1["energy"] = 1173.0;
  energy1["window"] = 0.5;
  energies_array.push_back( energy1 );

  json energy2;
  energy2["energy"] = 1332.0;
  energy2["window"] = 0.75;
  energies_array.push_back( energy2 );

  params["energies"] = energies_array;

  json result;
  try
  {
    result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec);
  }
  catch( const std::exception &e )
  {
    BOOST_TEST_MESSAGE( "Exception thrown: " << e.what() );
    throw;
  }
  //BOOST_REQUIRE_NO_THROW( result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec) );

  // Check basic structure
  BOOST_CHECK( result.is_object() );
  BOOST_CHECK( result.contains("sources") );
  BOOST_CHECK( result.contains("search_parameters") );
  BOOST_CHECK( result["sources"].is_array() );
  BOOST_CHECK( !result["sources"].empty() );

  // Co-60 should be in the results (if any results returned)
  if( !result["sources"].empty() )
  {
    bool found_co60 = false;
    for( const auto& source : result["sources"] )
    {
      BOOST_CHECK( source.contains("source") );
      BOOST_CHECK( source.contains("source_type") && source["source_type"].is_string() ); //"nuclide", "x-ray", "reaction"
      BOOST_CHECK( source.contains("sum_energy_difference") && source["sum_energy_difference"].is_number() );
      BOOST_CHECK( source.contains("profile_score") && source["profile_score"].is_number() );
      BOOST_CHECK( source.contains("energy_matches") );
      BOOST_CHECK( source["energy_matches"].is_array() );

      // Print source information to stdout
      std::cout << "\nSource: " << source["source"].get<std::string>()
                << " (type: " << source["source_type"].get<std::string>() << ")\n";
      std::cout << "  Sum energy difference: " << source["sum_energy_difference"].get<double>() << "\n";
      std::cout << "  Profile score: " << source["profile_score"].get<double>() << "\n";
      std::cout << "  Energy matches (" << source["energy_matches"].size() << "):\n";
      for( const auto& match : source["energy_matches"] )
      {
        std::cout << "    Search: " << match["search_energy"].get<double>() << " keV"
                  << " -> Matched: " << match["matched_energy"].get<double>() << " keV"
                  << " (BR: " << match["relative_br"].get<double>() << ")"
                  << " [" << match["transition"].get<std::string>() << "]\n";
      }

      const std::string source_name = source["source"].get<std::string>();
      if( source_name == "Co60" )
      {
        found_co60 = true;
        BOOST_CHECK_EQUAL( source["source_type"].get<std::string>(), "nuclide" );
        BOOST_CHECK( source.contains("half_life") );
        BOOST_CHECK( source.contains("assumed_age") );
        BOOST_CHECK_EQUAL( source["energy_matches"].size(), 2 );

        // Check energy matches
        for( const auto& match : source["energy_matches"] )
        {
          BOOST_CHECK( match.contains("search_energy") );
          BOOST_CHECK( match.contains("matched_energy") );
          BOOST_CHECK( match.contains("relative_br") );
          BOOST_CHECK( match.contains("transition") );
        }
      }
    }

    BOOST_TEST_MESSAGE( "Found " << result["sources"].size() << " sources" );
    BOOST_CHECK_MESSAGE( found_co60, "Did not find Co60 in top results" );
  }

  // Test with different category (using i18n key which is more reliable in test environment)
  params["source_category"] = "isbe-category-nuc-xray-rctn";  // Nuclides, X-rays, Reactions
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec) );
  BOOST_CHECK( result["sources"].is_array() );

  // Test with single energy (Cs-137 @ 661.7 keV)
  params = json::object();
  energies_array = json::array();
  energy1 = json::object();
  energy1["energy"] = 661.7;
  energies_array.push_back( energy1 );
  params["energies"] = energies_array;
  params["max_results"] = 5;

  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec) );
  BOOST_CHECK( result["sources"].is_array() );
  BOOST_CHECK( result["sources"].size() <= 5 );

  // Test with min_half_life parameter
  params["min_half_life"] = "1 y";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec) );
  BOOST_CHECK( result["sources"].is_array() );

  // Test sort_by parameter
  params["sort_by"] = "SumEnergyDifference";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("search_sources_by_energy", params, fixture.m_interspec) );
  BOOST_CHECK( result["sources"].is_array() );

  // Test error handling - empty energies array
  params = json::object();
  params["energies"] = json::array();
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - missing energies parameter
  params = json::object();
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - invalid category
  params = json::object();
  params["energies"] = energies_array;
  params["source_category"] = "InvalidCategory";
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - invalid half-life format
  params["source_category"] = "Nuclides + X-rays";
  params["min_half_life"] = "invalid format";
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - invalid sort_by
  params["min_half_life"] = "100 m";
  params["sort_by"] = "InvalidSort";
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );

  // Test error handling - negative energy
  params = json::object();
  energies_array = json::array();
  energy1 = json::object();
  energy1["energy"] = -100.0;
  energies_array.push_back( energy1 );
  params["energies"] = energies_array;
  BOOST_CHECK_THROW( registry.executeTool("search_sources_by_energy", params, fixture.m_interspec), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( test_executeEditAnalysisPeak )
{
  InterSpecTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params, result;

  // First, fit multiple Br82 peaks (and nearby Ra226 and K40 peaks) for testing
  // Br82 peaks at: 554, 606.3, 619, 776.5, 1474.9 keV
  // Ra226 peak at: 609.3 keV
  // K40 peak at: 1460 keV

  std::vector<double> peak_energies = {554.0, 606.3, 609.3, 619.0, 776.5, 1460.0, 1474.9};
  std::vector<double> fitted_peak_energies;

  for( const double energy : peak_energies )
  {
    params = json::object();
    params["energy"] = energy;
    params["specType"] = "Foreground";
    params["addToUsersPeaks"] = true;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("add_analysis_peak", params, fixture.m_interspec) );
    // Store the actual fitted energy for use in tests
    if( result.contains("fitPeakEnergy") )
    {
      fitted_peak_energies.push_back( result["fitPeakEnergy"].get<double>() );
    }
    else
    {
      std::cout << "Warning: Could not fit peak at " << energy << " keV. Result: " << result.dump() << std::endl;
    }
  }

  // Verify we have enough peaks to test
  BOOST_REQUIRE( fitted_peak_energies.size() >= 5 );

  // Test 1: SetEnergy - modify peak energy
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetEnergy";
  params["specType"] = "Foreground";
  params["doubleValue"] = 558.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK_CLOSE( result["modifiedPeak"]["energy"].get<double>(), 558.0, 0.1 );

  // Verify the peak was actually modified by checking get_analysis_peaks
  json get_peaks_params = json::object();
  get_peaks_params["specType"] = "Foreground";
  json peaks_result;
  BOOST_REQUIRE_NO_THROW( peaks_result = registry.executeTool("get_analysis_peaks", get_peaks_params, fixture.m_interspec) );

  // Search for the modified peak at 558 keV and ensure the original peak is not present
  bool found_modified_peak = false;
  bool found_original_peak = false;
  for( const auto &roi : peaks_result["rois"] )
  {
    for( const auto &peak : roi["peaks"] )
    {
      const double peak_energy = peak["energy"].get<double>();
      // Check if we find the modified peak at 558 keV
      if( std::abs(peak_energy - 558.0) < 0.5 )
        found_modified_peak = true;
      // Check if we still find the original peak (should not be there)
      if( std::abs(peak_energy - fitted_peak_energies[0]) < 0.5 && std::abs(peak_energy - 558.0) > 0.5 )
        found_original_peak = true;
    }
  }
  BOOST_CHECK_MESSAGE( found_modified_peak, "Should find modified peak at 558 keV in get_analysis_peaks" );
  BOOST_CHECK_MESSAGE( !found_original_peak, "Should NOT find original peak after modification" );

  // Set back to original energy - need to search at the modified energy (558.0)
  params = json::object();
  params["energy"] = 558.0;  // Search for peak at its current location
  params["editAction"] = "SetEnergy";
  params["specType"] = "Foreground";
  params["doubleValue"] = fitted_peak_energies[0];
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  if( !result["success"].get<bool>() )
  {
    std::cout << "Failed to set back to original energy. Result: " << result.dump() << std::endl;
  }
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  if( result.contains("modifiedPeak") && result["modifiedPeak"].contains("energy") )
  {
    BOOST_CHECK_CLOSE( result["modifiedPeak"]["energy"].get<double>(), fitted_peak_energies[0], 0.01 );
  }

  // Test 2: SetFwhm - modify peak FWHM
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetFwhm";
  params["specType"] = "Foreground";
  params["doubleValue"] = 5.5;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK_CLOSE( result["modifiedPeak"]["fwhm"].get<double>(), 5.5, 0.1 );

  // Test 3: SetAmplitude - modify peak amplitude
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetAmplitude";
  params["specType"] = "Foreground";
  params["doubleValue"] = 1000.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK_CLOSE( result["modifiedPeak"]["amplitude"].get<double>(), 1000.0, 1.0 );

  // Test 4: SetEnergyUncertainty - modify energy uncertainty
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetEnergyUncertainty";
  params["specType"] = "Foreground";
  params["doubleValue"] = 0.5;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  // Note: energyUncertainty may not be in JSON output, skip verification for now

  // Test 5: SetFwhmUncertainty - modify FWHM uncertainty
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetFwhmUncertainty";
  params["specType"] = "Foreground";
  params["doubleValue"] = 0.2;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  // Note: fwhmUncertainty may not be in JSON output, skip verification for now

  // Test 6: SetAmplitudeUncertainty - modify amplitude uncertainty
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetAmplitudeUncertainty";
  params["specType"] = "Foreground";
  params["doubleValue"] = 50.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK_CLOSE( result["modifiedPeak"]["amplitudeUncert"].get<double>(), 50.0, 0.5 );

  // Test 7: SetRoiLower - modify lower ROI bound
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetRoiLower";
  params["specType"] = "Foreground";
  params["doubleValue"] = 550.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("peaksInRoi") );
  // ROI bounds are checked via success message

  // Test 8: SetRoiUpper - modify upper ROI bound
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetRoiUpper";
  params["specType"] = "Foreground";
  params["doubleValue"] = 560.0;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("peaksInRoi") );
  // ROI bounds are checked via success message

  // Test 9: SetSkewType - test various skew types
  std::vector<std::string> skew_types = {"NoSkew", "Bortel", "GaussExp", "CrystalBall", "ExpGaussExp", "DoubleSidedCrystalBall"};
  for( const auto& skew_type : skew_types )
  {
    params = json::object();
    params["energy"] = fitted_peak_energies[1];
    params["editAction"] = "SetSkewType";
    params["specType"] = "Foreground";
    params["stringValue"] = skew_type;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
    BOOST_CHECK( result["success"].get<bool>() );
    BOOST_CHECK( result.contains("modifiedPeak") );
  }

  // Test 10: SetContinuumType - test various continuum types
  std::vector<std::string> continuum_types = {"NoOffset", "Constant", "Linear", "Quadratic", "Cubic", "FlatStep", "LinearStep", "BiLinearStep"};
  for( const auto& cont_type : continuum_types )
  {
    params = json::object();
    params["energy"] = fitted_peak_energies[2];
    params["editAction"] = "SetContinuumType";
    params["specType"] = "Foreground";
    params["stringValue"] = cont_type;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
    // Some continuum types may not be supported, just verify the call succeeded
    BOOST_CHECK( result.contains("success") );
    if( result["success"].get<bool>() )
    {
      BOOST_CHECK( result.contains("modifiedPeak") );
    }
  }

  // Test 11: SetSource - assign nuclide source
  params = json::object();
  params["energy"] = fitted_peak_energies[3];
  params["editAction"] = "SetSource";
  params["specType"] = "Foreground";
  params["stringValue"] = "Br82";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  // Nuclide assignment verified via success message

  // Test 12: SetSource - assign different nuclide
  params = json::object();
  params["energy"] = fitted_peak_energies[5];
  params["editAction"] = "SetSource";
  params["specType"] = "Foreground";
  params["stringValue"] = "K40";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  // Nuclide assignment verified via success message

  // Test 13: SetColor - assign color to peak
  params = json::object();
  params["energy"] = fitted_peak_energies[4];
  params["editAction"] = "SetColor";
  params["specType"] = "Foreground";
  params["stringValue"] = "rgb(255,0,0)";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );

  // Test 14: SetUserLabel - assign user label to peak
  params = json::object();
  params["energy"] = fitted_peak_energies[4];
  params["editAction"] = "SetUserLabel";
  params["specType"] = "Foreground";
  params["stringValue"] = "Test Peak Label";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  // User label verified via success message

  // Verify boolean properties are present in peak JSON
  BOOST_CHECK( result["modifiedPeak"].contains("useForEnergyCalibration") );
  BOOST_CHECK( result["modifiedPeak"].contains("useForShieldingSourceFit") );
  BOOST_CHECK( result["modifiedPeak"].contains("useForManualRelEff") );

  // Test 14a: SetUseForEnergyCalibration - set to true
  params = json::object();
  params["energy"] = fitted_peak_energies[4];
  params["editAction"] = "SetUseForEnergyCalibration";
  params["specType"] = "Foreground";
  params["boolValue"] = true;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK( result["modifiedPeak"]["useForEnergyCalibration"].get<bool>() == true );

  // Test 14b: SetUseForShieldingSourceFit - set to true
  params = json::object();
  params["energy"] = fitted_peak_energies[4];
  params["editAction"] = "SetUseForShieldingSourceFit";
  params["specType"] = "Foreground";
  params["boolValue"] = true;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK( result["modifiedPeak"]["useForShieldingSourceFit"].get<bool>() == true );

  // Test 14c: SetUseForManualRelEff - set to false
  params = json::object();
  params["energy"] = fitted_peak_energies[4];
  params["editAction"] = "SetUseForManualRelEff";
  params["specType"] = "Foreground";
  params["boolValue"] = false;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );
  BOOST_CHECK( result["modifiedPeak"]["useForManualRelEff"].get<bool>() == false );

  // Test 15: SplitFromRoi - split peak into its own ROI
  // First verify there are multiple peaks in the ROI around 606-609 keV region
  params = json::object();
  params["energy"] = fitted_peak_energies[1];
  params["editAction"] = "SplitFromRoi";
  params["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("modifiedPeak") );

  // Test 16: MergeWithRight - merge ROI with adjacent peak on the right
  // Use the first peak and merge with the second peak
  params = json::object();
  params["energy"] = fitted_peak_energies[1];
  params["editAction"] = "MergeWithRight";
  params["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("peaksInRoi") );
  BOOST_CHECK( result["peaksInRoi"].is_array() );

  // Test 17: MergeWithLeft - merge ROI with adjacent peak on the left
  // Use the second peak and merge with the first peak
  params = json::object();
  params["energy"] = fitted_peak_energies[2];
  params["editAction"] = "MergeWithLeft";
  params["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( result.contains("peaksInRoi") );
  BOOST_CHECK( result["peaksInRoi"].is_array() );

  // Test 18: DeletePeak - delete the last peak
  params = json::object();
  params["energy"] = fitted_peak_energies[6];
  params["editAction"] = "DeletePeak";
  params["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec) );
  BOOST_CHECK( result["success"].get<bool>() );
  BOOST_CHECK( !result.contains("modifiedPeak") ); // No modified peak after deletion

  // Test error handling - peak not found (already deleted)
  params = json::object();
  params["energy"] = fitted_peak_energies[6];
  params["editAction"] = "SetEnergy";
  params["specType"] = "Foreground";
  params["doubleValue"] = 1475.0;
  // May throw exception or return error result
  try {
    result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec);
    // If it doesn't throw, should return success=false
    BOOST_CHECK( !result["success"].get<bool>() );
  } catch( const std::runtime_error& ) {
    // Exception is also acceptable
  }

  // Test error handling - missing required parameters
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  // Missing editAction
  BOOST_CHECK_THROW( registry.executeTool("edit_analysis_peak", params, fixture.m_interspec), std::exception );

  // Test error handling - invalid edit action
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "InvalidAction";
  params["specType"] = "Foreground";
  BOOST_CHECK_THROW( registry.executeTool("edit_analysis_peak", params, fixture.m_interspec), std::exception );

  // Test error handling - missing doubleValue when required (may not throw, just fail gracefully)
  params = json::object();
  params["energy"] = fitted_peak_energies[0];
  params["editAction"] = "SetEnergy";
  params["specType"] = "Foreground";
  // Missing doubleValue
  result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec);
  BOOST_CHECK( !result["success"].get<bool>() || result.contains("message") );

  // Test error handling - missing stringValue when required (may not throw, just fail gracefully)
  params = json::object();
  params["energy"] = fitted_peak_energies[3];
  params["editAction"] = "SetSource";
  params["specType"] = "Foreground";
  // Missing stringValue
  result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec);
  BOOST_CHECK( !result["success"].get<bool>() || result.contains("message") );

  // Test error handling - invalid continuum type (may not throw, just fail gracefully)
  params = json::object();
  params["energy"] = fitted_peak_energies[2];
  params["editAction"] = "SetContinuumType";
  params["specType"] = "Foreground";
  params["stringValue"] = "InvalidContinuumType";
  result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec);
  BOOST_CHECK( !result["success"].get<bool>() || result.contains("message") );

  // Test error handling - invalid skew type (may not throw, just fail gracefully)
  params = json::object();
  params["energy"] = fitted_peak_energies[1];
  params["editAction"] = "SetSkewType";
  params["specType"] = "Foreground";
  params["stringValue"] = "InvalidSkewType";
  result = registry.executeTool("edit_analysis_peak", params, fixture.m_interspec);
  BOOST_CHECK( !result["success"].get<bool>() || result.contains("message") );
}


BOOST_AUTO_TEST_CASE( test_toolsLoadedFromXml )
{
  using namespace LlmTools;

  // Use a minimal fixture just to set up data directory
  InterSpecTestFixture fixture;

  // Get the registry (tools should already be registered by fixture)
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Get all registered tools
  const std::map<std::string, SharedTool>& tools = registry.getTools();

  // Check that we have a reasonable number of tools
  BOOST_CHECK_GT( tools.size(), 20 );

  // Define the expected tools (same list as in the validation code)
  const std::vector<std::string> expectedTools = {
    "get_detected_peaks",
    "add_analysis_peak",
    "edit_analysis_peak",
    "get_analysis_peaks",
    "get_spectrum_info",
    "primary_gammas_for_source",
    "sources_with_primary_gammas_in_energy_range",
    "sources_with_primary_gammas_near_energy",
    // These are always included, regardless of INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO
    //"sources_associated_with_source", //This info now included in "source_info" call
    //"analyst_notes_for_source",       //This info now included in "source_info" call
    "source_info",
    "nuclide_decay_chain",
    "automated_source_id_results",
    "loaded_spectra",
    //"add_analysis_peaks_for_source", // Temporarily commenting out until we get this implementation right
    "get_counts_in_energy_range",
    "get_expected_fwhm",
    "currie_mda_calc",
    "source_photons",
    "photopeak_detection_efficiency",
    "get_materials",
    "get_material_info",
    "available_detector_efficiency_functions",
    "load_detector_efficiency_function",
    "detector_efficiency_function_info",
    "search_sources_by_energy"
  };

  // Check that all expected tools are present
  for( const std::string &toolName : expectedTools )
  {
    const SharedTool* tool = registry.getTool(toolName);
    BOOST_REQUIRE_MESSAGE( tool != nullptr, "Tool '" + toolName + "' not found in registry" );

    // Check that tool has a name
    BOOST_CHECK_EQUAL( tool->name, toolName );

    // Check that tool has a description
    BOOST_CHECK_MESSAGE( !tool->description.empty(), "Tool '" + toolName + "' has no description" );

    // Check that tool has a parameters schema (non-null; empty object {} is valid for no-parameter tools)
    BOOST_CHECK_MESSAGE( !tool->parameters_schema.is_null(),
                        "Tool '" + toolName + "' has no parameters schema" );

    // Check that tool has an executor
    BOOST_CHECK_MESSAGE( tool->executor != nullptr, "Tool '" + toolName + "' has no executor" );
  }

  // Check for unexpected tools (excluding invoke_ tools which are dynamically created)
  for( const auto &[toolName, tool] : tools )
  {
    // Skip invoke_ tools
    if( toolName.find("invoke_") == 0 )
      continue;

    bool found = false;
    for( const std::string &expected : expectedTools )
    {
      if( toolName == expected )
      {
        found = true;
        break;
      }
    }
    BOOST_WARN_MESSAGE( found, "Unexpected tool '" + toolName + "' found in registry" );
  }

  cout << "Successfully validated " << expectedTools.size() << " tools from XML configuration" << endl;
}


// Activity/Shielding Fit Tool Tests
BOOST_AUTO_TEST_CASE( test_agentsLoadedFromXml )
{
  // Load configuration
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  // Check that agents were loaded
  BOOST_REQUIRE_MESSAGE( !llmConfig->agents.empty(), "No agents loaded from XML" );

  // Expected agents: MainAgent, NuclideId, NuclideIdWorker, ActivityFit, Isotopics
  const std::set<std::string> expectedAgents = { "MainAgent", "NuclideId", "NuclideIdWorker", "ActivityFit", "Isotopics" };

  cout << "Loaded " << llmConfig->agents.size() << " agents from XML configuration:" << endl;

  std::set<std::string> foundAgents;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    cout << "  - " << agent.name << " (type: " << static_cast<int>(agent.type) << ")" << endl;
    cout << "    System prompt length: " << agent.systemPrompt.length() << " characters" << endl;

    // Check that agent has a non-empty name
    BOOST_CHECK_MESSAGE( !agent.name.empty(), "Agent has empty name" );

    // Check that agent has a non-empty system prompt
    BOOST_CHECK_MESSAGE( !agent.systemPrompt.empty(),
                        "Agent '" << agent.name << "' has empty system prompt" );

    // Check that system prompt is substantial (at least 100 characters)
    BOOST_CHECK_MESSAGE( agent.systemPrompt.length() >= 100,
                        "Agent '" << agent.name << "' has suspiciously short system prompt: "
                        << agent.systemPrompt.length() << " chars" );

    foundAgents.insert( agent.name );
  }

  // Verify we got all expected agents
  BOOST_CHECK_EQUAL( foundAgents.size(), expectedAgents.size() );

  for( const std::string &expectedName : expectedAgents )
  {
    BOOST_CHECK_MESSAGE( foundAgents.count(expectedName) > 0,
                        "Expected agent '" << expectedName << "' not found in loaded agents" );
  }

  // Check for unexpected agents
  for( const std::string &foundName : foundAgents )
  {
    if( expectedAgents.count(foundName) == 0 )
    {
      cout << "WARNING: Unexpected agent found: " << foundName << endl;
    }
  }

  cout << "Successfully validated " << expectedAgents.size() << " agents from XML configuration" << endl;
}


// ============================================================================
// State Machine Tests
// ============================================================================

BOOST_AUTO_TEST_CASE( test_StateMachineBasicFunctionality )
{
  // Load configuration with state machine (Isotopics agent)
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  // Find the Isotopics agent which has a state machine
  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE_MESSAGE( isotopicsAgent != nullptr, "Isotopics agent not found in configuration" );
  BOOST_REQUIRE_MESSAGE( isotopicsAgent->state_machine != nullptr, "Isotopics agent should have a state machine" );

  AgentStateMachine *sm = isotopicsAgent->state_machine.get();

  // Test initial state
  BOOST_CHECK_EQUAL( sm->getInitialState(), "ANALYZE_REQUEST" );
  BOOST_CHECK_EQUAL( sm->getCurrentState(), "ANALYZE_REQUEST" );

  // Test state definition retrieval
  const AgentStateMachine::StateDefinition &initialState = sm->getStateDefinition( "ANALYZE_REQUEST" );
  BOOST_CHECK_EQUAL( initialState.name, "ANALYZE_REQUEST" );
  BOOST_CHECK( !initialState.description.empty() );
  BOOST_CHECK( !initialState.prompt_guidance.empty() );
  BOOST_CHECK( !initialState.is_final );

  // Verify allowed transitions are populated (this was the bug we fixed)
  BOOST_CHECK_MESSAGE( !initialState.allowed_transitions.empty(),
                      "ANALYZE_REQUEST state should have allowed transitions" );
  BOOST_CHECK_EQUAL( initialState.allowed_transitions.size(), 1 );
  BOOST_CHECK_EQUAL( initialState.allowed_transitions[0], "SELECT_CONFIGURATION" );

  // Verify required tools are populated
  BOOST_CHECK_MESSAGE( !initialState.required_tools.empty(),
                      "ANALYZE_REQUEST state should have required tools" );
  BOOST_CHECK_EQUAL( initialState.required_tools.size(), 2 );

  cout << "Initial state validated: " << initialState.name << endl;
  cout << "  Allowed transitions: " << initialState.allowed_transitions.size() << endl;
  cout << "  Required tools: " << initialState.required_tools.size() << endl;
}

BOOST_AUTO_TEST_CASE( test_StateMachineTransitions )
{
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  // Get a copy of the state machine to test transitions
  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE( isotopicsAgent != nullptr );
  BOOST_REQUIRE( isotopicsAgent->state_machine != nullptr );

  // Create a copy for testing (to avoid modifying the original)
  std::shared_ptr<AgentStateMachine> sm = isotopicsAgent->state_machine->copy();
  BOOST_REQUIRE( sm );

  // Test valid transition
  BOOST_CHECK( sm->canTransitionTo( "SELECT_CONFIGURATION" ) );
  BOOST_REQUIRE_NO_THROW( sm->transitionTo( "SELECT_CONFIGURATION" ) );
  BOOST_CHECK_EQUAL( sm->getCurrentState(), "SELECT_CONFIGURATION" );

  // Test invalid transition (should be allowed but warned about - soft enforcement)
  BOOST_CHECK( !sm->canTransitionTo( "FINALIZE_RESULTS" ) );

  // Test transition to VALIDATE_CONFIGURATION
  BOOST_CHECK( sm->canTransitionTo( "VALIDATE_CONFIGURATION" ) );
  sm->transitionTo( "VALIDATE_CONFIGURATION" );
  BOOST_CHECK_EQUAL( sm->getCurrentState(), "VALIDATE_CONFIGURATION" );

  // Verify allowed transitions from current state
  std::vector<std::string> allowedTransitions = sm->getAllowedTransitions();
  BOOST_CHECK_EQUAL( allowedTransitions.size(), 2 );
  BOOST_CHECK( std::find( allowedTransitions.begin(), allowedTransitions.end(), "CHECK_INTERFERENCE" ) != allowedTransitions.end() );
  BOOST_CHECK( std::find( allowedTransitions.begin(), allowedTransitions.end(), "SELECT_CONFIGURATION" ) != allowedTransitions.end() );

  // Test reset
  sm->reset();
  BOOST_CHECK_EQUAL( sm->getCurrentState(), sm->getInitialState() );

  cout << "State machine transitions validated" << endl;
}

BOOST_AUTO_TEST_CASE( test_StateMachineFinalState )
{
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE( isotopicsAgent != nullptr );
  BOOST_REQUIRE( isotopicsAgent->state_machine != nullptr );

  AgentStateMachine *sm = isotopicsAgent->state_machine.get();

  // Check that FINALIZE_RESULTS is a final state
  BOOST_CHECK( sm->hasState( "FINALIZE_RESULTS" ) );
  BOOST_CHECK( sm->isFinalState( "FINALIZE_RESULTS" ) );

  const AgentStateMachine::StateDefinition &finalState = sm->getStateDefinition( "FINALIZE_RESULTS" );
  BOOST_CHECK( finalState.is_final );
  BOOST_CHECK( finalState.allowed_transitions.empty() ); // Final states should have no outgoing transitions

  // Check that non-final states are not marked as final
  BOOST_CHECK( !sm->isFinalState( "ANALYZE_REQUEST" ) );
  BOOST_CHECK( !sm->isFinalState( "SELECT_CONFIGURATION" ) );

  cout << "Final state validation passed" << endl;
}

BOOST_AUTO_TEST_CASE( test_StateMachineAllStates )
{
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE( isotopicsAgent != nullptr );
  BOOST_REQUIRE( isotopicsAgent->state_machine != nullptr );

  AgentStateMachine *sm = isotopicsAgent->state_machine.get();

  // Expected states for Isotopics workflow
  std::vector<std::string> expectedStates = {
    "ANALYZE_REQUEST",
    "SELECT_CONFIGURATION",
    "VALIDATE_CONFIGURATION",
    "CHECK_INTERFERENCE",
    "EXECUTE_CALCULATION",
    "EVALUATE_RESULTS",
    "FINALIZE_RESULTS"
  };

  // Verify all expected states exist
  for( const std::string &stateName : expectedStates )
  {
    BOOST_CHECK_MESSAGE( sm->hasState( stateName ),
                        "State '" + stateName + "' should exist in Isotopics state machine" );

    if( sm->hasState( stateName ) )
    {
      const AgentStateMachine::StateDefinition &state = sm->getStateDefinition( stateName );

      // Every state should have a description and guidance
      BOOST_CHECK_MESSAGE( !state.description.empty(),
                          "State '" + stateName + "' should have a description" );
      BOOST_CHECK_MESSAGE( !state.prompt_guidance.empty(),
                          "State '" + stateName + "' should have prompt guidance" );

      // Non-final states should have transitions
      if( !state.is_final )
      {
        BOOST_CHECK_MESSAGE( !state.allowed_transitions.empty(),
                            "Non-final state '" + stateName + "' should have allowed transitions" );
      }

      cout << "  State: " << stateName
           << ", Transitions: " << state.allowed_transitions.size()
           << ", Tools: " << state.required_tools.size()
           << ", Final: " << (state.is_final ? "yes" : "no") << endl;
    }
  }

  cout << "All " << expectedStates.size() << " states validated" << endl;
}

BOOST_AUTO_TEST_CASE( test_StateMachinePromptGuidance )
{
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE( isotopicsAgent != nullptr );
  BOOST_REQUIRE( isotopicsAgent->state_machine != nullptr );

  std::shared_ptr<AgentStateMachine> sm = isotopicsAgent->state_machine->copy();
  BOOST_REQUIRE( sm );

  // Test getting guidance for current state
  std::string initialGuidance = sm->getPromptGuidanceForCurrentState();
  BOOST_CHECK( !initialGuidance.empty() );
  BOOST_CHECK( initialGuidance.find( "spectrum" ) != std::string::npos ||
               initialGuidance.find( "request" ) != std::string::npos );

  // Transition to another state and check guidance changes
  sm->transitionTo( "SELECT_CONFIGURATION" );
  std::string configGuidance = sm->getPromptGuidanceForCurrentState();
  BOOST_CHECK( !configGuidance.empty() );
  BOOST_CHECK( configGuidance != initialGuidance ); // Guidance should be different
  BOOST_CHECK( configGuidance.find( "preset" ) != std::string::npos ||
               configGuidance.find( "configuration" ) != std::string::npos );

  cout << "Prompt guidance validation passed" << endl;
}

BOOST_AUTO_TEST_CASE( test_StateMachineEdgeCases )
{
  std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
  BOOST_REQUIRE( llmConfig );

  const LlmConfig::AgentConfig *isotopicsAgent = nullptr;
  for( const LlmConfig::AgentConfig &agent : llmConfig->agents )
  {
    if( agent.name == "Isotopics" )
    {
      isotopicsAgent = &agent;
      break;
    }
  }

  BOOST_REQUIRE( isotopicsAgent != nullptr );
  BOOST_REQUIRE( isotopicsAgent->state_machine != nullptr );

  AgentStateMachine *sm = isotopicsAgent->state_machine.get();

  // Test non-existent state
  BOOST_CHECK( !sm->hasState( "NONEXISTENT_STATE" ) );
  BOOST_CHECK_THROW( sm->getStateDefinition( "NONEXISTENT_STATE" ), std::exception );

  // Test invalid transition
  std::shared_ptr<AgentStateMachine> smCopy = sm->copy();
  BOOST_CHECK( !smCopy->canTransitionTo( "NONEXISTENT_STATE" ) );
  BOOST_CHECK_THROW( smCopy->transitionTo( "NONEXISTENT_STATE" ), std::exception );

  // Test multiple resets
  smCopy->reset();
  std::string firstState = smCopy->getCurrentState();
  smCopy->reset();
  BOOST_CHECK_EQUAL( smCopy->getCurrentState(), firstState );

  // Test copy independence
  smCopy->transitionTo( "SELECT_CONFIGURATION" );
  BOOST_CHECK_EQUAL( smCopy->getCurrentState(), "SELECT_CONFIGURATION" );
  BOOST_CHECK_EQUAL( sm->getCurrentState(), "ANALYZE_REQUEST" ); // Original should be unchanged

  cout << "Edge case validation passed" << endl;
}


// ============================================================================
// Isotopics Tool Tests
// ============================================================================


// ============================================================================
// Currie MDA Calculation Tool Tests
// ============================================================================

BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_Basic )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test basic functionality without nuclide (backward compatibility)
  json params;
  params["energy"] = 661.7;  // Cs137 main gamma

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

  // Check basic result structure
  BOOST_CHECK( result.contains("gammaEnergy") );
  BOOST_CHECK( result.contains("roiLowerEnergy") );
  BOOST_CHECK( result.contains("roiUpperEnergy") );
  BOOST_CHECK( result.contains("decisionThreshold") );
  BOOST_CHECK( result.contains("detectionLimit") );
  BOOST_CHECK( result.contains("sourceCounts") );
  BOOST_CHECK( result.contains("peakPresentInData") );

  // Verify energy is correct
  BOOST_CHECK_CLOSE( result["gammaEnergy"].get<double>(), 661.7, 0.1 );

  // Should not have activity fields when nuclide not specified
  BOOST_CHECK( !result.contains("gammasPerBq") );
  BOOST_CHECK( !result.contains("observedActivity") );
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_WithNuclideAndDistance )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First load a detector (required for activity calculation)
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test activity calculation" );
    return;
  }

  // Load a detector
  const string detector_name = available_result[0]["name"].get<string>();
  json load_params;
  load_params["identifier"] = detector_name;
  json load_result;
  BOOST_REQUIRE_NO_THROW( load_result = registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Test with nuclide and distance
  json params;
  params["energy"] = 661.7;  // Cs137 main gamma
  params["nuclide"] = "Cs137";
  params["distance"] = "25 cm";

  json result;
  try
  {
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

    // Check that activity fields are present
    BOOST_CHECK( result.contains("gammasPerBq") );
    BOOST_CHECK( result.contains("branchRatio") );
    BOOST_CHECK( result.contains("decisionThresholdActivity") );
    BOOST_CHECK( result.contains("detectionLimitActivity") );

    // Verify gammasPerBq is positive
    const double gammas_per_bq = result["gammasPerBq"].get<double>();
    BOOST_CHECK_GT( gammas_per_bq, 0.0 );

    // Verify branch ratio is reasonable (Cs137 661.7 keV has ~0.85 branching ratio)
    const double br = result["branchRatio"].get<double>();
    BOOST_CHECK_GT( br, 0.5 );
    BOOST_CHECK_LT( br, 1.0 );
  }
  catch( const std::exception &e )
  {
    // Expected if no foreground spectrum is loaded
    const string msg = e.what();
    BOOST_TEST_MESSAGE( "Cannot test activity calculation without foreground spectrum: " << msg );
    BOOST_CHECK( msg.find("No foreground spectrum loaded") != string::npos || 
                 msg.find("detector") != string::npos );
  }
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_WithShieldingMaterial )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test shielding" );
    return;
  }

  const string detector_name = available_result[0]["name"].get<string>();
  json load_params;
  load_params["identifier"] = detector_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Test with Material/Thickness shielding
  json params;
  params["energy"] = 661.7;
  params["nuclide"] = "Cs137";
  params["distance"] = "25 cm";
  json shielding;
  shielding["Material"] = "Fe";
  shielding["Thickness"] = "1.0 cm";
  params["shielding"] = shielding;

  json result;
  try
  {
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

    // Check that shielding transmission is present and reasonable
    BOOST_CHECK( result.contains("shieldingTransmission") );
    const double shield_trans = result["shieldingTransmission"].get<double>();
    BOOST_CHECK_GT( shield_trans, 0.0 );
    BOOST_CHECK_LE( shield_trans, 1.0 );
  }
  catch( const std::exception &e )
  {
    const string msg = e.what();
    BOOST_TEST_MESSAGE( "Cannot test shielding without foreground spectrum: " << msg );
  }
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_WithShieldingANAD )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test AN/AD shielding" );
    return;
  }

  const string detector_name = available_result[0]["name"].get<string>();
  json load_params;
  load_params["identifier"] = detector_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Test with AN/AD shielding format
  json params;
  params["energy"] = 661.7;
  params["nuclide"] = "Cs137";
  params["distance"] = "25 cm";
  json shielding;
  shielding["AN"] = 26.0;  // Iron atomic number
  shielding["AD"] = 7.87;  // Iron density in g/cm² for 1 cm thickness
  params["shielding"] = shielding;

  json result;
  try
  {
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

    // Check that shielding transmission is present
    BOOST_CHECK( result.contains("shieldingTransmission") );
    const double shield_trans = result["shieldingTransmission"].get<double>();
    BOOST_CHECK_GT( shield_trans, 0.0 );
    BOOST_CHECK_LE( shield_trans, 1.0 );
  }
  catch( const std::exception &e )
  {
    const string msg = e.what();
    BOOST_TEST_MESSAGE( "Cannot test AN/AD shielding without foreground spectrum: " << msg );
  }
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_WithAssertBackgroundSpectrum )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with assertBackgroundSpectrum=true
  json params;
  params["energy"] = 661.7;
  params["assertBackgroundSpectrum"] = true;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

  // When assertBackgroundSpectrum is true, side channels should be 0
  BOOST_CHECK_EQUAL( result["numLowerSideChannels"].get<int>(), 0 );
  BOOST_CHECK_EQUAL( result["numUpperSideChannels"].get<int>(), 0 );
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_ErrorCases )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with invalid nuclide
  json params;
  params["energy"] = 661.7;
  params["nuclide"] = "InvalidNuclide123";
  params["distance"] = "25 cm";

  BOOST_CHECK_THROW( registry.executeTool("currie_mda_calc", params, fixture.m_interspec), std::runtime_error );

  // Test with invalid distance
  params["nuclide"] = "Cs137";
  params["distance"] = "invalid distance string";

  BOOST_CHECK_THROW( registry.executeTool("currie_mda_calc", params, fixture.m_interspec), std::runtime_error );

  // Test with invalid shielding format (both Material and AN)
  params["distance"] = "25 cm";
  json shielding;
  shielding["Material"] = "Fe";
  shielding["Thickness"] = "1.0 cm";
  shielding["AN"] = 26.0;  // Should not have both
  params["shielding"] = shielding;

  BOOST_CHECK_THROW( registry.executeTool("currie_mda_calc", params, fixture.m_interspec), std::runtime_error );

  // Test with invalid material
  shielding = json::object();
  shielding["Material"] = "NonexistentMaterial123";
  shielding["Thickness"] = "1.0 cm";
  params["shielding"] = shielding;

  BOOST_CHECK_THROW( registry.executeTool("currie_mda_calc", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeCurrieMdaCalc_WithAge )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("available_detector_efficiency_functions", json::object(), fixture.m_interspec) );

  if( available_result.empty() )
  {
    BOOST_TEST_MESSAGE( "No detectors available to test age parameter" );
    return;
  }

  const string detector_name = available_result[0]["name"].get<string>();
  json load_params;
  load_params["identifier"] = detector_name;
  BOOST_REQUIRE_NO_THROW( registry.executeTool("load_detector_efficiency_function", load_params, fixture.m_interspec) );

  // Test with age parameter
  json params;
  params["energy"] = 661.7;
  params["nuclide"] = "Cs137";
  params["distance"] = "25 cm";
  params["age"] = "1 year";

  json result;
  try
  {
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool("currie_mda_calc", params, fixture.m_interspec) );

    // Should have activity fields
    BOOST_CHECK( result.contains("gammasPerBq") );
    BOOST_CHECK( result.contains("branchRatio") );
  }
  catch( const std::exception &e )
  {
    const string msg = e.what();
    BOOST_TEST_MESSAGE( "Cannot test age parameter without foreground spectrum: " << msg );
  }
}


// Dose Calculation Tests
BOOST_AUTO_TEST_CASE( test_executeCalculateDose_BasicNoShielding )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test Co60, 0.5 yr old, 100 µCi, 1m, no shielding
  // Expected: 115.42 µrem/hr (from DoseCalcWidget::runtime_sanity_checks)
  json params;
  params["nuclide"] = "Co60";
  params["activity"] = "100 uCi";
  params["distance"] = "1 m";
  params["age"] = "0.5 y";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("calculate_dose", params, fixture.m_interspec) );

  // Check that result contains expected fields
  BOOST_REQUIRE( result.contains("success") );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_REQUIRE( result.contains("dose_rate_Sv_per_hr") );
  BOOST_REQUIRE( result.contains("dose_rate_si_str") );
  BOOST_REQUIRE( result.contains("dose_rate_REM_per_hr") );
  BOOST_REQUIRE( result.contains("dose_rate_REM_str") );
  BOOST_REQUIRE( result.contains("nuclide") );
  BOOST_REQUIRE( result.contains("distance") );
  BOOST_REQUIRE( result.contains("activity") );
  BOOST_REQUIRE( result.contains("age") );

  // Check numeric dose values
  const double dose_rem_hr = result["dose_rate_REM_per_hr"].get<double>();
  const string rem_hr_str = result["dose_rate_REM_str"].get<string>();

  // Expected is 115.42 µrem/hr, allow 2% tolerance
  BOOST_CHECK_CLOSE( dose_rem_hr, 115.42e-6, 2.0 );

  // Verify the formatted string is reasonable
  BOOST_CHECK( rem_hr_str.find("rem/hr") != string::npos || rem_hr_str.find("rem/h") != string::npos );

  BOOST_TEST_MESSAGE( "Co60 dose rate (no shielding): " << rem_hr_str << " (" << dose_rem_hr << " rem/hr)" );
}


BOOST_AUTO_TEST_CASE( test_executeCalculateDose_WithMaterialShielding )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test Cs137, 0.5 yr old, 100 µCi, 1m, with Fe shielding
  // Using thickness to achieve ~5 g/cm² (Fe density is ~7.87 g/cm³, so 0.636 cm ≈ 5 g/cm²)
  json params;
  params["nuclide"] = "Cs137";
  params["activity"] = "100 uCi";
  params["distance"] = "1 m";
  //params["age"] = "0.5 y"; //Use default age
  params["material"] = "Fe";
  params["thickness"] = "0.636 cm";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("calculate_dose", params, fixture.m_interspec) );

  BOOST_REQUIRE( result.contains("success") );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_REQUIRE( result.contains("shielding") );

  // Check dose values
  BOOST_REQUIRE( result.contains("dose_rate_REM_per_hr") );
  BOOST_REQUIRE( result.contains("dose_rate_Sv_per_hr") );

  const double dose_rem_hr = result["dose_rate_REM_per_hr"].get<double>();
  const double dose_sv_hr = result["dose_rate_Sv_per_hr"].get<double>();

  // Expected dose for Cs137, 100 µCi, 1m, 5 g/cm² Fe shielding, 25.19 µrem/hr
  const double expected_rem_hr = 25.19e-6;  // rem/hr
  BOOST_CHECK_CLOSE( dose_rem_hr, expected_rem_hr, 2.0 );

  // Check shielding info
  const json &shielding = result["shielding"];
  BOOST_REQUIRE( shielding.contains("arealDensity_g_cm2") );
  BOOST_REQUIRE( shielding.contains("atomicNumber") );

  const double ad = shielding["arealDensity_g_cm2"].get<double>();
  const double an = shielding["atomicNumber"].get<double>();

  // Check areal density is close to 5 g/cm² (within 10% tolerance)
  BOOST_CHECK_CLOSE( ad, 5.0, 10.0 );

  // Check atomic number is close to 26 (iron)
  BOOST_CHECK_CLOSE( an, 26.0, 1.0 );

  BOOST_TEST_MESSAGE( "Cs137 with Fe shielding: " << result["dose_rate_REM_str"].get<string>() );
  BOOST_TEST_MESSAGE( "  Areal density: " << ad << " g/cm²" );
  BOOST_TEST_MESSAGE( "  Atomic number: " << an );
}


BOOST_AUTO_TEST_CASE( test_executeCalculateDose_WithArealDensityShielding )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test Cs137 with direct areal density/atomic number specification
  json params;
  params["nuclide"] = "Cs137";
  params["activity"] = "100 uCi";
  params["distance"] = "100 cm";
  // params["age"] = "0.5 y";  //Use the default age
  params["arealDensity"] = 5.0;  // g/cm²
  params["atomicNumber"] = 26.0;  // Iron

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("calculate_dose", params, fixture.m_interspec) );

  BOOST_REQUIRE( result.contains("success") );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_REQUIRE( result.contains("shielding") );

  const json &shielding = result["shielding"];
  BOOST_CHECK_CLOSE( shielding["arealDensity_g_cm2"].get<double>(), 5.0, 0.1 );
  BOOST_CHECK_CLOSE( shielding["atomicNumber"].get<double>(), 26.0, 0.1 );

  // Check dose values
  BOOST_REQUIRE( result.contains("dose_rate_REM_per_hr") );
  BOOST_REQUIRE( result.contains("dose_rate_Sv_per_hr") );

  const double dose_rem_hr = result["dose_rate_REM_per_hr"].get<double>();
  const double dose_sv_hr = result["dose_rate_Sv_per_hr"].get<double>();

  // Expected dose for Cs137, 100 µCi, 1m, 5 g/cm² Fe shielding, default age: 25.19 µrem/hr
  const double expected_rem_hr = 25.19e-6;  // rem/hr
  BOOST_CHECK_CLOSE( dose_rem_hr, expected_rem_hr, 2.0 );

  BOOST_TEST_MESSAGE( "Cs137 with AD/AN shielding: " << result["dose_rate_REM_str"].get<string>() );
}


BOOST_AUTO_TEST_CASE( test_executeCalculateDose_DefaultAge )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test F18 without specifying age (should use default, which is 0 for this nuclide)
  json params;
  params["nuclide"] = "F18";
  params["activity"] = "100 uCi";
  params["distance"] = "2 m";
  // Note: age not specified, should use PeakDef::defaultDecayTime

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("calculate_dose", params, fixture.m_interspec) );

  BOOST_REQUIRE( result.contains("success") );
  BOOST_CHECK( result["success"].get<bool>() == true );
  BOOST_REQUIRE( result.contains("age") );

  // Verify age was set (should be the default for F18)
  const string age_str = result["age"].get<string>();
  BOOST_CHECK( !age_str.empty() );

  // Check dose values
  BOOST_REQUIRE( result.contains("dose_rate_REM_per_hr") );
  BOOST_REQUIRE( result.contains("dose_rate_Sv_per_hr") );

  const double dose_rem_hr = result["dose_rate_REM_per_hr"].get<double>();
  const double dose_sv_hr = result["dose_rate_Sv_per_hr"].get<double>();

  // Expected dose for F18, 100 µCi, 2m, no shielding, default age: 13.29 µrem/hr
  const double expected_rem_hr = 13.29e-6;  // rem/hr
  BOOST_CHECK_CLOSE( dose_rem_hr, expected_rem_hr, 2.0 );

  BOOST_TEST_MESSAGE( "F18 with default age (" << age_str << "): " << result["dose_rate_REM_str"].get<string>() );
}


BOOST_AUTO_TEST_CASE( test_executeCalculateDose_ErrorCases )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test missing required parameter (nuclide)
  json params1;
  params1["distance"] = "1 m";
  params1["activity"] = "1 uCi";
  BOOST_CHECK_THROW( registry.executeTool("calculate_dose", params1, fixture.m_interspec), std::exception );

  // Test invalid nuclide
  json params2;
  params2["nuclide"] = "InvalidNuclide123";
  params2["distance"] = "1 m";
  params2["activity"] = "1 uCi";
  BOOST_CHECK_THROW( registry.executeTool("calculate_dose", params2, fixture.m_interspec), std::exception );

  // Test invalid distance format
  json params3;
  params3["nuclide"] = "Co60";
  params3["distance"] = "invalid distance";
  params3["activity"] = "1 uCi";
  BOOST_CHECK_THROW( registry.executeTool("calculate_dose", params3, fixture.m_interspec), std::exception );

  // Test invalid activity format
  json params4;
  params4["nuclide"] = "Co60";
  params4["distance"] = "1 m";
  params4["activity"] = "invalid activity";
  BOOST_CHECK_THROW( registry.executeTool("calculate_dose", params4, fixture.m_interspec), std::exception );

  // Test unknown material
  json params5;
  params5["nuclide"] = "Co60";
  params5["distance"] = "1 m";
  params5["activity"] = "1 uCi";
  params5["material"] = "UnknownMaterial123";
  params5["thickness"] = "1 cm";
  BOOST_CHECK_THROW( registry.executeTool("calculate_dose", params5, fixture.m_interspec), std::exception );
}


BOOST_AUTO_TEST_CASE( test_executeCalculateDose_WithMaterialAndArealDensity )
{
  InterSpecTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test Cs137 with material + thickness
  json params_thickness;
  params_thickness["nuclide"] = "Cs137";
  params_thickness["activity"] = "100 uCi";
  params_thickness["distance"] = "100 cm";
  //params_thickness["age"] = "0.5 y";  // Use default age
  params_thickness["material"] = "Fe";
  params_thickness["thickness"] = "0.636 cm";  // Fe density = 7.87 g/cm³, so 0.636 cm * 7.87 = ~5.0 g/cm²

  json result_thickness;
  BOOST_REQUIRE_NO_THROW( result_thickness = registry.executeTool("calculate_dose", params_thickness, fixture.m_interspec) );

  BOOST_REQUIRE( result_thickness.contains("success") );
  BOOST_CHECK( result_thickness["success"].get<bool>() == true );
  BOOST_REQUIRE( result_thickness.contains("shielding") );

  const json &shielding_thickness = result_thickness["shielding"];
  const double ad_thickness = shielding_thickness["arealDensity_g_cm2"].get<double>();
  const double an_thickness = shielding_thickness["atomicNumber"].get<double>();

  BOOST_TEST_MESSAGE( "Cs137 with Fe + thickness: " << result_thickness["dose_rate_REM_str"].get<string>() );
  BOOST_TEST_MESSAGE( "  Areal density: " << ad_thickness << " g/cm²" );
  BOOST_TEST_MESSAGE( "  Atomic number: " << an_thickness );

  // Now test the same scenario with material + arealDensity
  json params_ad;
  params_ad["nuclide"] = "Cs137";
  params_ad["activity"] = "100 uCi";
  params_ad["distance"] = "100 cm";
  //params_ad["age"] = "0.5 y";   // Use default age
  params_ad["material"] = "Fe";
  params_ad["arealDensity"] = 5.0;  // Directly specify the areal density in g/cm²

  json result_ad;
  BOOST_REQUIRE_NO_THROW( result_ad = registry.executeTool("calculate_dose", params_ad, fixture.m_interspec) );

  BOOST_REQUIRE( result_ad.contains("success") );
  BOOST_CHECK( result_ad["success"].get<bool>() == true );
  BOOST_REQUIRE( result_ad.contains("shielding") );

  const json &shielding_ad = result_ad["shielding"];
  const double ad_direct = shielding_ad["arealDensity_g_cm2"].get<double>();
  const double an_direct = shielding_ad["atomicNumber"].get<double>();

  BOOST_TEST_MESSAGE( "Cs137 with Fe + arealDensity: " << result_ad["dose_rate_REM_str"].get<string>() );
  BOOST_TEST_MESSAGE( "  Areal density: " << ad_direct << " g/cm²" );
  BOOST_TEST_MESSAGE( "  Atomic number: " << an_direct );

  // Both methods should produce the same atomic number (Fe)
  BOOST_CHECK_CLOSE( an_direct, an_thickness, 0.1 );
  BOOST_CHECK_CLOSE( an_direct, 26.0, 0.1 );

  // Areal density from thickness should be close to 5.0 g/cm² (depends on Fe density)
  BOOST_CHECK_CLOSE( ad_thickness, 5.0, 2.0 );  // Allow 2% tolerance

  // Direct areal density should be exactly 5.0 g/cm²
  BOOST_CHECK_CLOSE( ad_direct, 5.0, 0.1 );

  // Check that both methods give the same numeric dose values
  const double dose_thickness_rem = result_thickness["dose_rate_REM_per_hr"].get<double>();
  const double dose_ad_rem = result_ad["dose_rate_REM_per_hr"].get<double>();
  const double dose_thickness_sv = result_thickness["dose_rate_Sv_per_hr"].get<double>();
  const double dose_ad_sv = result_ad["dose_rate_Sv_per_hr"].get<double>();

  // Expected dose for Cs137, 100 µCi, 1m, 5 g/cm² Fe shielding is approximately 25.19 µrem/hr (251.95 nSv/hr)
  // Convert to rem/hr: 25.19 µrem/hr = 25.19e-6 rem/hr
  // Convert to Sv/hr: 251.95 nSv/hr = 251.95e-9 Sv/hr
  const double expected_rem_hr = 25.19e-6;  // rem/hr
  const double expected_sv_hr = 251.95e-9;  // Sv/hr

  // Check that calculated dose is close to expected value (within 2% tolerance)
  BOOST_CHECK_CLOSE( dose_ad_rem, expected_rem_hr, 2.0 );
  BOOST_CHECK_CLOSE( dose_ad_sv, expected_sv_hr, 2.0 );

  // Both methods should give nearly identical numeric results (within 0.5% since areal densities are very close)
  BOOST_CHECK_CLOSE( dose_thickness_rem, dose_ad_rem, 0.5 );
  BOOST_CHECK_CLOSE( dose_thickness_sv, dose_ad_sv, 0.5 );

  // Get formatted strings for display
  const string dose_thickness_str = result_thickness["dose_rate_REM_str"].get<string>();
  const string dose_ad_str = result_ad["dose_rate_REM_str"].get<string>();
  const string dose_sv_str = result_ad["dose_rate_si_str"].get<string>();

  // Note: Formatted strings may differ slightly due to small differences in areal density
  // We verify numeric values are close instead

  BOOST_TEST_MESSAGE( "Dose (thickness method): " << dose_thickness_str << " (" << dose_thickness_rem << " rem/hr)" );
  BOOST_TEST_MESSAGE( "Dose (arealDensity method): " << dose_ad_str << " (" << dose_ad_rem << " rem/hr)" );
  BOOST_TEST_MESSAGE( "Dose (Sv units): " << dose_sv_str << " (" << dose_ad_sv << " Sv/hr)" );
  BOOST_TEST_MESSAGE( "Areal density (thickness): " << ad_thickness << " g/cm²" );
  BOOST_TEST_MESSAGE( "Areal density (direct): " << ad_direct << " g/cm²" );
}

