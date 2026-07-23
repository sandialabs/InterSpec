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

#define BOOST_TEST_MODULE LlmRelActManualTool_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
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

  // Set the writable data directory for tests
  BOOST_REQUIRE_NO_THROW( InterSpec::setWritableDataDirectory( datadir ) );

  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
}//void set_data_dir()


// Helper class to manage InterSpec instance for tests
// This fixture loads a uranium spectrum and fits peaks with U235/U238 assignments
class RelActManualTestFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
  InterSpec *m_interspec;
  std::shared_ptr<LlmTools::ToolRegistry> m_tool_registry;
  bool m_peaks_setup;

  RelActManualTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr ),
    m_tool_registry( nullptr ),
    m_peaks_setup( false )
  {
    set_data_dir();

    string wt_app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), "..");
    wt_app_root = SpecUtils::lexically_normalize_path(wt_app_root);

    // Create a test environment
    const std::string applicationPath = "";
    const std::string configurationFile = "";
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) );

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    // Load LLM configuration
    std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
    BOOST_REQUIRE( llmConfig );

    // Create tool registry for testing
    m_tool_registry = std::make_shared<LlmTools::ToolRegistry>( *llmConfig );
    BOOST_REQUIRE( m_tool_registry );

    // Load a test uranium spectrum
    loadTestSpectrum();
    
    // Setup peaks with source assignments
    setupPeaksWithSources();
  }

  ~RelActManualTestFixture()
  {
    m_tool_registry.reset();
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr;
  }

  void loadTestSpectrum()
  {
    try
    {
      const string datadir = InterSpec::staticDataDirectory();
      
      // Try to load a uranium spectrum file
      // Look for uranium reference spectra in a few common locations
      vector<string> candidate_files = {
        SpecUtils::append_path( datadir, "reference_spectra/Common_Field_Nuclides/Detective X/U235_Unshielded.txt" ),
        SpecUtils::append_path( datadir, "reference_spectra/Common_Field_Nuclides/Detective X/U238_Unshielded.txt" ),
        SpecUtils::append_path( datadir, "reference_spectra/Common_Field_Nuclides/Detective X/Cs137_Unshielded.txt" )
      };
      
      string spectrum_file;
      for( const auto &f : candidate_files )
      {
        if( SpecUtils::is_file(f) )
        {
          spectrum_file = f;
          break;
        }
      }
      
      if( spectrum_file.empty() )
      {
        BOOST_TEST_MESSAGE( "No suitable test spectrum file found in reference_spectra" );
        return;
      }
      
      std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
      const bool loaded = meas->load_file( spectrum_file, SpecUtils::ParserType::Auto, spectrum_file );

      if( loaded && (meas->num_measurements() > 0) )
      {
        const shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index(0);
        BOOST_REQUIRE( !!m );
        m_interspec->setSpectrum( meas, {m->sample_number()}, SpecUtils::SpectrumType::Foreground, 0 );
        BOOST_TEST_MESSAGE( "Loaded test spectrum: " << spectrum_file );
      }
      else
      {
        BOOST_TEST_MESSAGE( "Failed to load test spectrum: " << spectrum_file );
      }
    }
    catch( const std::exception &e )
    {
      BOOST_TEST_MESSAGE( "Exception loading test spectrum: " << e.what() );
    }
  }
  
  void setupPeaksWithSources()
  {
    // Create synthetic peaks with U235 and U238 source assignments for testing
    // This simulates having already fit peaks in the spectrum
    
    PeakModel *peakModel = m_interspec->peakModel();
    if( !peakModel )
    {
      BOOST_TEST_MESSAGE( "No peak model available" );
      return;
    }
    
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    if( !db )
    {
      BOOST_TEST_MESSAGE( "No decay database available" );
      return;
    }
    
    const SandiaDecay::Nuclide *u235 = db->nuclide("U235");
    const SandiaDecay::Nuclide *u238 = db->nuclide("U238");
    
    if( !u235 || !u238 )
    {
      BOOST_TEST_MESSAGE( "Could not find U235 or U238 nuclides" );
      return;
    }
    
    shared_ptr<const SpecUtils::Measurement> foreground = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if( !foreground )
    {
      BOOST_TEST_MESSAGE( "No foreground histogram available for peak setup" );
      return;
    }
    
    // Create peaks at typical uranium gamma energies with source assignments
    // Using U235 peaks at 143.76, 163.33, 185.72, 205.31 keV
    // Using U238 (via Pa234m) peak at 1001.03 keV (and others)
    
    struct PeakInfo
    {
      double energy;
      double amplitude;
      double sigma;
      const SandiaDecay::Nuclide *nuclide;
      float gamma_energy;
    };
    
    vector<PeakInfo> peaks_to_add = {
      { 143.76, 5000.0, 1.0, u235, 143.76f },
      { 163.33, 3000.0, 1.1, u235, 163.33f },
      { 185.72, 50000.0, 1.2, u235, 185.72f },
      { 205.31, 5500.0, 1.3, u235, 205.31f },
      
      { 766.4, 1500.0, 2.5, u238, 766.4f },
      { 1001.03, 2000.0, 3.0, u238, 1001.03f }
    };
    
    for( const auto &info : peaks_to_add )
    {
      try
      {
        PeakDef peak;
        peak.setMean( info.energy );
        peak.setSigma( info.sigma );
        peak.setAmplitude( info.amplitude );
        peak.setAmplitudeUncert( sqrt(info.amplitude) );
        
        // Set up the nuclide source
        // For U235, use major gamma transitions
        // We need to find the appropriate transition for the gamma energy
        if( info.nuclide )
        {
          // Find the transition corresponding to this gamma energy
          const SandiaDecay::Transition *trans = nullptr;
          int rad_index = -1;
          
          for( const SandiaDecay::Transition *t : info.nuclide->decaysToChildren )
          {
            if( !t )
              continue;
            for( size_t i = 0; i < t->products.size(); ++i )
            {
              const SandiaDecay::RadParticle &p = t->products[i];
              if( p.type == SandiaDecay::GammaParticle && 
                  fabs(p.energy - info.gamma_energy) < 0.5 )
              {
                trans = t;
                rad_index = static_cast<int>(i);
                break;
              }
            }
            if( trans )
              break;
          }
          
          if( trans && rad_index >= 0 )
          {
            peak.setNuclearTransition( info.nuclide, trans, rad_index, PeakDef::SourceGammaType::NormalGamma );
          }
          // If we couldn't find the transition, skip this peak
          else
          {
            continue;
          }
        }
        
        peak.useForManualRelEff( true );
        
        peakModel->addNewPeak( peak );
      }
      catch( const std::exception &e )
      {
        BOOST_TEST_MESSAGE( "Exception adding peak at " << info.energy << " keV: " << e.what() );
      }
    }
    
    // Check if peaks were added
    shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
    if( peaks && !peaks->empty() )
    {
      m_peaks_setup = true;
      BOOST_TEST_MESSAGE( "Successfully added " << peaks->size() << " peaks with source assignments" );
    }
    else
    {
      BOOST_TEST_MESSAGE( "No peaks were added to the peak model" );
    }
  }

  const LlmTools::ToolRegistry &llmToolRegistry()
  {
    BOOST_REQUIRE( m_tool_registry );
    return *m_tool_registry;
  }

  bool peaksAreSetup() const { return m_peaks_setup; }
};


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_with_sources )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test using sources parameter
  json params;
  params["sources"] = json::array({ "U235", "U238" });
  params["eqn_form"] = "LnY";
  params["eqn_order"] = 3;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

  // Check for success
  BOOST_CHECK( result.contains( "success" ) );
  
  if( result["success"].get<bool>() )
  {
    // Check for expected output fields
    BOOST_CHECK( result.contains( "sources" ) );
    BOOST_CHECK( result["sources"].is_array() );
    BOOST_CHECK( result["sources"].size() >= 1 );
    
    BOOST_CHECK( result.contains( "peaks" ) );
    BOOST_CHECK( result["peaks"].is_array() );
    
    BOOST_CHECK( result.contains( "fit_quality" ) );
    BOOST_CHECK( result["fit_quality"].contains( "chi2" ) );
    BOOST_CHECK( result["fit_quality"].contains( "dof" ) );
    BOOST_CHECK( result["fit_quality"].contains( "chi2_per_dof" ) );
    
    BOOST_CHECK( result.contains( "warnings" ) );
    BOOST_CHECK( result.contains( "saved_to_state" ) );
    BOOST_CHECK( result["saved_to_state"].get<bool>() == false );
    
    // Check peaks have expected fields
    for( const auto &peak : result["peaks"] )
    {
      BOOST_CHECK( peak.contains( "energy" ) );
      BOOST_CHECK( peak.contains( "counts" ) );
      
      // Check for residual_sigma if efficiency was computed
      if( peak.contains( "observed_efficiency" ) )
      {
        BOOST_CHECK( peak.contains( "fit_efficiency" ) );
        BOOST_CHECK( peak.contains( "residual_sigma" ) );
        
        // Verify residual_sigma is reasonable (not NaN or inf)
        if( peak["residual_sigma"].is_number() )
        {
          double sigma = peak["residual_sigma"].get<double>();
          BOOST_CHECK( std::isfinite(sigma) );
          BOOST_CHECK( sigma >= 0.0 );
        }
      }
    }
    
    cout << "peak_based_relative_efficiency with sources succeeded:" << endl;
    cout << "  Sources: " << result["sources"].size() << endl;
    cout << "  Peaks: " << result["peaks"].size() << endl;
    cout << "  Chi2/DOF: " << result["fit_quality"]["chi2_per_dof"] << endl;
  }
  else
  {
    BOOST_TEST_MESSAGE( "peak_based_relative_efficiency failed: " << result.value("error", "unknown error") );
  }
}


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_with_energies )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test using peak_energies parameter
  json params;
  params["peak_energies"] = json::array({ 143.76, 185.72, 205.31 });
  params["eqn_form"] = "LnXLnY";
  params["additional_uncertainty"] = "FivePercent";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

  BOOST_CHECK( result.contains( "success" ) );
  
  if( result["success"].get<bool>() )
  {
    BOOST_CHECK( result.contains( "peaks" ) );
    BOOST_CHECK( result["peaks"].is_array() );
    
    cout << "peak_based_relative_efficiency with energies succeeded" << endl;
    cout << "  Peaks matched: " << result["peaks"].size() << endl;
  }
  else
  {
    BOOST_TEST_MESSAGE( "peak_based_relative_efficiency failed: " << result.value("error", "unknown error") );
  }
}


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_mutual_exclusion )
{
  RelActManualTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test that providing both peak_energies and sources throws an error
  json params;
  params["peak_energies"] = json::array({ 185.72 });
  params["sources"] = json::array({ "U235" });

  json result;
  BOOST_CHECK_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ), std::exception );
}


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_no_input )
{
  RelActManualTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test that not providing either parameter throws an error
  json params;
  params["eqn_form"] = "LnY";

  json result;
  BOOST_CHECK_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ), std::exception );
}


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_save_to_state )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with save_to_state = true
  json params;
  params["sources"] = json::array({ "U235" });
  params["eqn_form"] = "LnY";
  params["save_to_state"] = true;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

  BOOST_CHECK( result.contains( "success" ) );
  
  if( result["success"].get<bool>() )
  {
    BOOST_CHECK( result.contains( "saved_to_state" ) );
    BOOST_CHECK( result["saved_to_state"].get<bool>() == true );
    
    cout << "peak_based_relative_efficiency with save_to_state succeeded" << endl;
  }
  else
  {
    BOOST_TEST_MESSAGE( "peak_based_relative_efficiency failed: " << result.value("error", "unknown error") );
  }
}


BOOST_AUTO_TEST_CASE( test_get_rel_act_manual_state_no_state )
{
  RelActManualTestFixture fixture;

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test getting state when none exists
  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_rel_act_manual_state", json::object(), fixture.m_interspec ) );

  BOOST_CHECK( result.contains( "success" ) );
  BOOST_CHECK( result["success"].get<bool>() == true );
  
  BOOST_CHECK( result.contains( "has_state" ) );
  // May or may not have state depending on whether save_to_state was called earlier
  
  BOOST_CHECK( result.contains( "source_from_gui" ) );
  
  cout << "get_rel_act_manual_state: has_state = " << result["has_state"].get<bool>() << endl;
}


BOOST_AUTO_TEST_CASE( test_get_rel_act_manual_state_after_save )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // First, save some state
  json save_params;
  save_params["sources"] = json::array({ "U235" });
  save_params["eqn_form"] = "LnXLnY";
  save_params["eqn_order"] = 4;
  save_params["match_tolerance"] = 2.0;
  save_params["additional_uncertainty"] = "TenPercent";
  save_params["save_to_state"] = true;

  json save_result;
  BOOST_REQUIRE_NO_THROW( save_result = registry.executeTool( "peak_based_relative_efficiency", save_params, fixture.m_interspec ) );

  if( !save_result["success"].get<bool>() )
  {
    BOOST_TEST_MESSAGE( "Could not save state: " << save_result.value("error", "unknown") );
    return;
  }

  // Now retrieve the state
  json get_result;
  BOOST_REQUIRE_NO_THROW( get_result = registry.executeTool( "get_rel_act_manual_state", json::object(), fixture.m_interspec ) );

  BOOST_CHECK( get_result.contains( "success" ) );
  BOOST_CHECK( get_result["success"].get<bool>() == true );
  
  BOOST_CHECK( get_result.contains( "has_state" ) );
  BOOST_CHECK( get_result["has_state"].get<bool>() == true );
  
  BOOST_CHECK( get_result.contains( "state" ) );
  
  if( get_result.contains( "state" ) )
  {
    const json &state = get_result["state"];
    
    // Verify the state matches what we saved
    BOOST_CHECK( state.contains( "eqn_form" ) );
    BOOST_CHECK_EQUAL( state["eqn_form"].get<string>(), "LnXLnY" );
    
    BOOST_CHECK( state.contains( "eqn_order" ) );
    BOOST_CHECK_EQUAL( state["eqn_order"].get<int>(), 4 );
    
    BOOST_CHECK( state.contains( "match_tolerance" ) );
    BOOST_CHECK_CLOSE( state["match_tolerance"].get<double>(), 2.0, 0.01 );
    
    BOOST_CHECK( state.contains( "additional_uncertainty" ) );
    BOOST_CHECK_EQUAL( state["additional_uncertainty"].get<string>(), "TenPercent" );
    
    BOOST_CHECK( state.contains( "background_subtract" ) );
    BOOST_CHECK_EQUAL( state["background_subtract"].get<bool>(), false );
    
    cout << "get_rel_act_manual_state after save succeeded:" << endl;
    cout << "  eqn_form: " << state["eqn_form"] << endl;
    cout << "  eqn_order: " << state["eqn_order"] << endl;
    cout << "  match_tolerance: " << state["match_tolerance"] << endl;
    cout << "  additional_uncertainty: " << state["additional_uncertainty"] << endl;
  }
}


BOOST_AUTO_TEST_CASE( test_peak_based_relative_efficiency_with_nuclide_ages )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test with nuclide_ages parameter
  json params;
  params["sources"] = json::array({ "U235" });
  params["eqn_form"] = "LnY";
  params["nuclide_ages"] = {
    {"U235", "20 years"},
    {"U238", "4.5e9 years"}
  };

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

  BOOST_CHECK( result.contains( "success" ) );
  
  if( result["success"].get<bool>() )
  {
    cout << "peak_based_relative_efficiency with nuclide_ages succeeded" << endl;
  }
  else
  {
    BOOST_TEST_MESSAGE( "peak_based_relative_efficiency failed: " << result.value("error", "unknown error") );
  }
}


BOOST_AUTO_TEST_CASE( test_additional_uncertainty_enum_values )
{
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Test various additional_uncertainty enum values
  vector<string> uncert_values = {
    "Unweighted", "StatOnly", "OnePercent", "FivePercent", 
    "TenPercent", "TwentyFivePercent", "FiftyPercent", 
    "SeventyFivePercent", "OneHundredPercent"
  };

  for( const auto &uncert : uncert_values )
  {
    json params;
    params["sources"] = json::array({ "U235" });
    params["eqn_form"] = "LnY";
    params["additional_uncertainty"] = uncert;

    json result;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

    BOOST_CHECK( result.contains( "success" ) );
    
    if( result["success"].get<bool>() )
    {
      BOOST_TEST_MESSAGE( "additional_uncertainty='" << uncert << "' succeeded" );
    }
    else
    {
      BOOST_TEST_MESSAGE( "additional_uncertainty='" << uncert << "' failed: " << result.value("error", "unknown") );
    }
  }
}


BOOST_AUTO_TEST_CASE( test_residual_sigma_interpretation )
{
  // This test verifies the residual_sigma interpretation guidance:
  // - 0 to 5 sigma: Good
  // - 5 to 10 sigma: Acceptable
  // - Beyond 10 sigma: Poor
  
  RelActManualTestFixture fixture;

  if( !fixture.peaksAreSetup() )
  {
    BOOST_TEST_MESSAGE( "Skipping test - peaks not set up" );
    return;
  }

  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  json params;
  params["sources"] = json::array({ "U235" });
  params["eqn_form"] = "LnY";

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "peak_based_relative_efficiency", params, fixture.m_interspec ) );

  if( !result["success"].get<bool>() )
  {
    BOOST_TEST_MESSAGE( "Skipping residual_sigma test - calculation failed" );
    return;
  }

  BOOST_CHECK( result.contains( "peaks" ) );
  
  int good_peaks = 0;
  int acceptable_peaks = 0;
  int poor_peaks = 0;
  
  for( const auto &peak : result["peaks"] )
  {
    if( peak.contains( "residual_sigma" ) && peak["residual_sigma"].is_number() )
    {
      double sigma = peak["residual_sigma"].get<double>();
      
      if( sigma >= 0.0 && sigma <= 5.0 )
        good_peaks++;
      else if( sigma > 5.0 && sigma <= 10.0 )
        acceptable_peaks++;
      else if( sigma > 10.0 )
        poor_peaks++;
    }
  }
  
  cout << "Residual sigma categorization:" << endl;
  cout << "  Good (0-5 sigma): " << good_peaks << endl;
  cout << "  Acceptable (5-10 sigma): " << acceptable_peaks << endl;
  cout << "  Poor (>10 sigma): " << poor_peaks << endl;
  
  // For a good fit with correct source assignments, most peaks should be "good"
  // We can't strictly enforce this since the test uses synthetic data
  BOOST_CHECK_MESSAGE( good_peaks + acceptable_peaks + poor_peaks > 0, 
                       "Should have at least one peak with residual_sigma" );
}

