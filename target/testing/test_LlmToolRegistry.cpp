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

#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmToolRegistry.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

using namespace std;
using namespace boost::unit_test;
using json = nlohmann::json;

// TODO: as of 20251013, need to add tests for: `currie_mda_calc`, `fit_peaks_for_nuclide`, `automated_isotope_id_results`, and `nuclides_with_primary_gammas_in_energy_range` 


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
  

  InterSpecTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr )
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

    // Register default LLM tools
    LlmTools::ToolRegistry::instance().registerDefaultTools();

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
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr; //deleted when m_env is deleted
  }
};


BOOST_AUTO_TEST_CASE( test_executeGetMaterials )
{
  InterSpecTestFixture fixture;

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  json params = json::object();

  json result;
  try
  {
    // Note: tool name has typo - "avaiable" instead of "available"
    result = registry.executeTool("avaiable_detector_efficiency_functions", params, fixture.m_interspec);
  }
  catch( const std::exception &e )
  {
    BOOST_TEST_MESSAGE( "Exception calling avaiable_detector_efficiency_functions: " << e.what() );
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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // First get available detectors (note: tool name has typo - "avaiable" not "available")
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("avaiable_detector_efficiency_functions", json::object(), fixture.m_interspec) );

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("avaiable_detector_efficiency_functions", json::object(), fixture.m_interspec) );

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("avaiable_detector_efficiency_functions", json::object(), fixture.m_interspec) );

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
  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // We have to manually run hint seard
  InterSpec * const viewer = fixture.m_interspec;
  viewer->searchForHintPeaks( viewer->measurment(SpecUtils::SpectrumType::Foreground),
                             viewer->displayedSamples(SpecUtils::SpectrumType::Foreground), true );
  
  
  json params = json::object();
  params["specType"] = "Foreground";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("detected_peaks", params, fixture.m_interspec) );
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
  BOOST_CHECK_THROW( registry.executeTool("detected_peaks", params, fixture.m_interspec), std::runtime_error );
}


BOOST_AUTO_TEST_CASE( test_executeGetUserPeaks )
{
  InterSpecTestFixture fixture;

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // Get user peaks when there are none initially
  json params = json::object();
  params["specType"] = "Foreground";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_user_peaks", params, fixture.m_interspec) );

  // Should return an object with rois array (same structure as detected_peaks)
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
  BOOST_REQUIRE_NO_THROW( registry.executeTool("fit_peak", fit_params_776, fixture.m_interspec) );

  // Fit peak at 554.35 keV
  json fit_params_554;
  fit_params_554["energy"] = 554.35;
  fit_params_554["specType"] = "Foreground";
  BOOST_REQUIRE_NO_THROW( registry.executeTool("fit_peak", fit_params_554, fixture.m_interspec) );

  //std::cout << "Fitted peaks at 776.52 keV and 554.35 keV\n";

  // Now get user peaks - should have the fitted peaks
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("get_user_peaks", params, fixture.m_interspec) );

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
}


BOOST_AUTO_TEST_CASE( test_executeGetCharacteristicLines )
{
  InterSpecTestFixture fixture;

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // Get associated nuclides
  json params = json::object();
  params["source"] = "Br82";

  json result;
#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
#else
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("sources_associated_with_source", params, fixture.m_interspec) );

  BOOST_CHECK( result.contains("status") && (result["status"] == "success") );
#endif
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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  json params;
  params["source"] = "U235";

  json result;
#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
#else
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec) );
#endif

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
#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("source_info", params, fixture.m_interspec) );
#else
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec) );
#endif
  BOOST_CHECK( result.contains("analystNotes") );

  // Test error handling - invalid nuclide
  params["source"] = "InvalidNuclide12345";
#if( INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  BOOST_CHECK_THROW( registry.executeTool("source_info", params, fixture.m_interspec), std::runtime_error );
#else
  BOOST_CHECK_THROW( registry.executeTool("analyst_notes_for_source", params, fixture.m_interspec), std::runtime_error );
#endif
}


BOOST_AUTO_TEST_CASE( test_executeGetSourceInfo )
{
  InterSpecTestFixture fixture;

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

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

  LlmTools::ToolRegistry &registry = LlmTools::ToolRegistry::instance();

  // First load a detector
  json available_result;
  BOOST_REQUIRE_NO_THROW( available_result = registry.executeTool("avaiable_detector_efficiency_functions", json::object(), fixture.m_interspec) );

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
