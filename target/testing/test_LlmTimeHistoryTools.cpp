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
#include <vector>
#include <iostream>

#include <Wt/Utils>
#include <Wt/WApplication>
#include <Wt/Test/WTestEnvironment>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE LlmTimeHistoryTools_suite
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

  const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ), "sandia.decay.xml not at '" << sandia_decay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  BOOST_REQUIRE_NO_THROW( InterSpec::setWritableDataDirectory( datadir ) );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide( "U238" ), "SandiaDecayDataBase empty?" );
}//void set_data_dir()


// Helper class to manage InterSpec instance for time history tests
class TimeHistoryTestFixture
{
public:
  std::unique_ptr<Wt::Test::WTestEnvironment> m_env;
  InterSpecApp *m_app;
  std::unique_ptr<Wt::WApplication::UpdateLock> m_update_lock;
  InterSpec *m_interspec;
  std::shared_ptr<LlmTools::ToolRegistry> m_tool_registry;


  TimeHistoryTestFixture()
  : m_env( nullptr ),
    m_app( nullptr ),
    m_interspec( nullptr ),
    m_tool_registry( nullptr )
  {
    set_data_dir();

    string wt_app_root = SpecUtils::append_path( InterSpec::staticDataDirectory(), ".." );
    wt_app_root = SpecUtils::lexically_normalize_path( wt_app_root );

    const std::string applicationPath = "";
    const std::string configurationFile = "";
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::Application ) );
    m_env->setAppRoot( wt_app_root );

    m_app = new InterSpecApp( *m_env );
    m_update_lock.reset( new Wt::WApplication::UpdateLock( m_app ) );

    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    std::shared_ptr<const LlmConfig> llmConfig = LlmConfig::load();
    BOOST_REQUIRE( llmConfig );

    m_interspec->createLlmTool();

    m_tool_registry = std::make_shared<LlmTools::ToolRegistry>( *llmConfig );
    BOOST_REQUIRE( m_tool_registry );
  }

  ~TimeHistoryTestFixture()
  {
    m_tool_registry.reset();
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
    m_app = nullptr;
  }


  const LlmTools::ToolRegistry &llmToolRegistry()
  {
    BOOST_REQUIRE( m_tool_registry );
    return *m_tool_registry;
  }


  // Load a passthrough (search-mode) N42 file and set it as foreground with all samples displayed
  void loadPassthroughFile( const std::string &filename )
  {
    const string filepath = SpecUtils::append_path(
      SpecUtils::append_path( g_test_file_dir, "llm_time_history" ), filename );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( filepath ), "Test file not found: " << filepath );

    std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
    const bool loaded = meas->load_file( filepath, SpecUtils::ParserType::Auto, filepath );
    BOOST_REQUIRE_MESSAGE( loaded, "Failed to load: " << filepath );
    BOOST_REQUIRE_MESSAGE( meas->num_measurements() > 0, "No measurements in: " << filepath );
    BOOST_REQUIRE_MESSAGE( meas->passthrough(), "File is not passthrough: " << filepath );

    // Set all sample numbers as foreground
    const std::set<int> &allSamples = meas->sample_numbers();
    BOOST_REQUIRE( !allSamples.empty() );
    m_interspec->setSpectrum( meas, allSamples, SpecUtils::SpectrumType::Foreground, 0 );

    // Verify it loaded correctly
    const std::shared_ptr<SpecMeas> loadedMeas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    BOOST_REQUIRE( loadedMeas );
    BOOST_REQUIRE( loadedMeas->passthrough() );
  }//loadPassthroughFile(...)


  // Load a non-passthrough spectrum file
  void loadNonPassthroughFile()
  {
    const string datadir = InterSpec::staticDataDirectory();
    const string spectrum_file = SpecUtils::append_path( datadir,
      "reference_spectra/Common_Field_Nuclides/Detective X/Br82_Unshielded.txt" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( spectrum_file ), "Non-passthrough test file not found: " << spectrum_file );

    std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
    const bool loaded = meas->load_file( spectrum_file, SpecUtils::ParserType::Auto, spectrum_file );
    BOOST_REQUIRE_MESSAGE( loaded, "Failed to load non-passthrough file" );
    BOOST_REQUIRE_MESSAGE( meas->num_measurements() > 0, "No measurements in non-passthrough file" );

    const std::shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index( 0 );
    BOOST_REQUIRE( m );
    m_interspec->setSpectrum( meas, {m->sample_number()}, SpecUtils::SpectrumType::Foreground, 0 );
  }//loadNonPassthroughFile()
};// class TimeHistoryTestFixture


BOOST_AUTO_TEST_SUITE( LlmTimeHistoryTools )


// The three passthrough test files
static const std::vector<std::string> s_passthrough_files = {
  "Q345_Racetrack.n42",
  "Q67_Mobile.N42",
  "Q8910_Helicopter.N42"
};


BOOST_AUTO_TEST_CASE( test_get_spectrum_info_searchMode )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a passthrough file and verify searchMode is true
  fixture.loadPassthroughFile( "Q345_Racetrack.n42" );

  json params;
  params["specType"] = "Foreground";
  params["fileInfo"] = true;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_spectrum_info", params, fixture.m_interspec ) );
  BOOST_REQUIRE( result.is_object() );
  BOOST_REQUIRE( result.contains( "fileInfo" ) );
  BOOST_REQUIRE( result["fileInfo"].contains( "searchMode" ) );
  BOOST_CHECK_EQUAL( result["fileInfo"]["searchMode"].get<bool>(), true );

  // Now load a non-passthrough file and verify searchMode is false
  fixture.loadNonPassthroughFile();

  json result2;
  BOOST_REQUIRE_NO_THROW( result2 = registry.executeTool( "get_spectrum_info", params, fixture.m_interspec ) );
  BOOST_REQUIRE( result2.is_object() );
  BOOST_REQUIRE( result2.contains( "fileInfo" ) );
  BOOST_REQUIRE( result2["fileInfo"].contains( "searchMode" ) );
  BOOST_CHECK_EQUAL( result2["fileInfo"]["searchMode"].get<bool>(), false );
}//test_get_spectrum_info_searchMode


BOOST_AUTO_TEST_CASE( test_get_time_history_count_rates_basic )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  for( const std::string &filename : s_passthrough_files )
  {
    BOOST_TEST_MESSAGE( "Testing get_time_history_count_rates with file: " << filename );

    fixture.loadPassthroughFile( filename );

    json params = json::object();
    json result;
    BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_time_history_count_rates", params, fixture.m_interspec ) );

    // Verify result structure - now uses timeSlices array of objects
    BOOST_REQUIRE( result.is_object() );
    BOOST_REQUIRE_MESSAGE( result.contains( "timeSlices" ), "Missing timeSlices for " << filename );
    BOOST_REQUIRE_MESSAGE( result.contains( "numEntries" ), "Missing numEntries for " << filename );

    const json &timeSlices = result["timeSlices"];
    const int numEntries = result["numEntries"].get<int>();

    BOOST_REQUIRE( timeSlices.is_array() );
    BOOST_CHECK( numEntries > 0 );
    BOOST_CHECK_EQUAL( timeSlices.size(), static_cast<size_t>( numEntries ) );

    // Each entry should have startTime, duration, gammaCps
    for( size_t i = 0; i < timeSlices.size(); ++i )
    {
      const json &slice = timeSlices[i];
      BOOST_REQUIRE_MESSAGE( slice.contains( "startTime" ), "Missing startTime at index " << i << " in " << filename );
      BOOST_REQUIRE_MESSAGE( slice.contains( "duration" ), "Missing duration at index " << i << " in " << filename );
      BOOST_REQUIRE_MESSAGE( slice.contains( "gammaCps" ), "Missing gammaCps at index " << i << " in " << filename );

      BOOST_CHECK_MESSAGE( slice["gammaCps"].get<double>() >= 0.0,
        "Negative gammaCps at index " << i << " in " << filename );
      BOOST_CHECK_MESSAGE( slice["duration"].get<double>() > 0.0,
        "Non-positive duration at index " << i << " in " << filename );
    }

    // startTimes should be monotonically non-decreasing
    for( size_t i = 1; i < timeSlices.size(); ++i )
    {
      const double prev = timeSlices[i - 1]["startTime"].get<double>();
      const double curr = timeSlices[i]["startTime"].get<double>();
      BOOST_CHECK_MESSAGE( curr >= prev,
        "startTime not monotonic at index " << i << " in " << filename
        << ": " << prev << " > " << curr );
    }
  }//for( each passthrough file )
}//test_get_time_history_count_rates_basic


BOOST_AUTO_TEST_CASE( test_get_time_history_count_rates_maxEntries )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q345_Racetrack.n42" );

  // First, get the full count to know the total number of entries
  json fullParams = json::object();
  json fullResult;
  BOOST_REQUIRE_NO_THROW( fullResult = registry.executeTool( "get_time_history_count_rates", fullParams, fixture.m_interspec ) );
  const int fullNumEntries = fullResult["numEntries"].get<int>();
  BOOST_REQUIRE( fullNumEntries > 10 ); // Need more than 10 samples for this test to be meaningful

  // Now request at most 10 entries
  json params;
  params["maxEntries"] = 10;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_time_history_count_rates", params, fixture.m_interspec ) );

  const int numEntries = result["numEntries"].get<int>();
  BOOST_CHECK_MESSAGE( numEntries <= 10, "Expected at most 10 entries, got " << numEntries );
  BOOST_CHECK_EQUAL( result["timeSlices"].size(), static_cast<size_t>( numEntries ) );

  // Even with maxEntries, the startTimes should still be monotonically non-decreasing
  const json &slices = result["timeSlices"];
  for( size_t i = 1; i < slices.size(); ++i )
  {
    BOOST_CHECK( slices[i]["startTime"].get<double>() >= slices[i - 1]["startTime"].get<double>() );
  }
}//test_get_time_history_count_rates_maxEntries


BOOST_AUTO_TEST_CASE( test_get_time_history_count_rates_time_range )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q67_Mobile.N42" );

  // Get full data to know the time span
  json fullParams = json::object();
  json fullResult;
  BOOST_REQUIRE_NO_THROW( fullResult = registry.executeTool( "get_time_history_count_rates", fullParams, fixture.m_interspec ) );

  const json &fullSlices = fullResult["timeSlices"];
  BOOST_REQUIRE( fullSlices.size() >= 3 );

  const double firstTime = fullSlices[0]["startTime"].get<double>();
  const double lastTime = fullSlices[fullSlices.size() - 1]["startTime"].get<double>()
                         + fullSlices[fullSlices.size() - 1]["duration"].get<double>();
  const double totalSpan = lastTime - firstTime;
  BOOST_REQUIRE( totalSpan > 0.0 );

  // Request a sub-range in the middle third
  const double rangeStart = firstTime + (totalSpan / 3.0);
  const double rangeEnd = firstTime + (2.0 * totalSpan / 3.0);

  json params;
  params["startTime"] = rangeStart;
  params["endTime"] = rangeEnd;

  json result;
  BOOST_REQUIRE_NO_THROW( result = registry.executeTool( "get_time_history_count_rates", params, fixture.m_interspec ) );

  const int numEntries = result["numEntries"].get<int>();
  BOOST_CHECK( numEntries > 0 );
  BOOST_CHECK_MESSAGE( numEntries < fullResult["numEntries"].get<int>(),
    "Filtered result should have fewer entries than full result" );

  // All returned slices should overlap with requested range
  const json &filtSlices = result["timeSlices"];
  for( size_t i = 0; i < filtSlices.size(); ++i )
  {
    const double sampleStart = filtSlices[i]["startTime"].get<double>();
    const double sampleEnd = sampleStart + filtSlices[i]["duration"].get<double>();
    BOOST_CHECK_MESSAGE( sampleEnd >= rangeStart && sampleStart <= rangeEnd,
      "Sample at index " << i << " [" << sampleStart << ", " << sampleEnd
      << "] does not overlap with range [" << rangeStart << ", " << rangeEnd << "]" );
  }
}//test_get_time_history_count_rates_time_range


BOOST_AUTO_TEST_CASE( test_get_time_history_count_rates_energy_filter )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q8910_Helicopter.N42" );

  // Get full-range (unfiltered) count rates
  json fullParams = json::object();
  json fullResult;
  BOOST_REQUIRE_NO_THROW( fullResult = registry.executeTool( "get_time_history_count_rates", fullParams, fixture.m_interspec ) );

  // Get energy-filtered count rates (100 to 200 keV window)
  json filteredParams;
  filteredParams["energyMin"] = 100.0;
  filteredParams["energyMax"] = 200.0;

  json filteredResult;
  BOOST_REQUIRE_NO_THROW( filteredResult = registry.executeTool( "get_time_history_count_rates", filteredParams, fixture.m_interspec ) );

  // Both should have the same number of entries (energy filter does not remove time entries)
  BOOST_CHECK_EQUAL( fullResult["numEntries"].get<int>(), filteredResult["numEntries"].get<int>() );

  // Filtered CPS should be <= unfiltered CPS for each entry
  const json &fullSlicesE = fullResult["timeSlices"];
  const json &filtSlicesE = filteredResult["timeSlices"];
  BOOST_REQUIRE_EQUAL( fullSlicesE.size(), filtSlicesE.size() );

  for( size_t i = 0; i < fullSlicesE.size(); ++i )
  {
    const double fullVal = fullSlicesE[i]["gammaCps"].get<double>();
    const double filtVal = filtSlicesE[i]["gammaCps"].get<double>();
    BOOST_CHECK_MESSAGE( filtVal <= (fullVal + 1.0),
      "Filtered CPS (" << filtVal << ") > unfiltered CPS (" << fullVal
      << ") at index " << i );
  }
}//test_get_time_history_count_rates_energy_filter


BOOST_AUTO_TEST_CASE( test_set_displayed_time_ranges )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q345_Racetrack.n42" );

  // Get the full time span first
  json countRatesParams = json::object();
  json countRatesResult;
  BOOST_REQUIRE_NO_THROW( countRatesResult = registry.executeTool( "get_time_history_count_rates", countRatesParams, fixture.m_interspec ) );

  const json &tsSlices = countRatesResult["timeSlices"];
  BOOST_REQUIRE( tsSlices.size() >= 4 );

  const double firstTime = tsSlices[0]["startTime"].get<double>();
  const double lastTime = tsSlices[tsSlices.size() - 1]["startTime"].get<double>()
                         + tsSlices[tsSlices.size() - 1]["duration"].get<double>();
  const double totalSpan = lastTime - firstTime;
  BOOST_REQUIRE( totalSpan > 0.0 );

  // Set foreground to the first third and background to the last third
  const double fgEnd = firstTime + (totalSpan / 3.0);
  const double bgStart = firstTime + (2.0 * totalSpan / 3.0);

  const string fgRange = to_string( firstTime ) + "-" + to_string( fgEnd );
  const string bgRange = to_string( bgStart ) + "-" + to_string( lastTime );

  json setParams;
  setParams["foreground"] = fgRange;
  setParams["background"] = bgRange;

  json setResult;
  BOOST_REQUIRE_NO_THROW( setResult = registry.executeTool( "set_displayed_time_ranges", setParams, fixture.m_interspec ) );

  BOOST_REQUIRE( setResult.is_object() );
  BOOST_CHECK( setResult.contains( "status" ) );
  BOOST_CHECK_EQUAL( setResult["status"].get<string>(), "success" );

  // Should have foreground and background info
  BOOST_CHECK( setResult.contains( "foreground" ) );
  BOOST_CHECK( setResult.contains( "foregroundNumSamples" ) );
  BOOST_CHECK( setResult.contains( "background" ) );
  BOOST_CHECK( setResult.contains( "backgroundNumSamples" ) );

  // Foreground and background sample counts should be > 0
  BOOST_CHECK( setResult["foregroundNumSamples"].get<int>() > 0 );
  BOOST_CHECK( setResult["backgroundNumSamples"].get<int>() > 0 );

  // Verify that displayed samples actually changed
  const std::set<int> &fgSamples = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  BOOST_CHECK( !fgSamples.empty() );
  const std::set<int> &bgSamples = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Background );
  BOOST_CHECK( !bgSamples.empty() );

  // Foreground and background samples should not overlap
  for( const int sample : fgSamples )
  {
    BOOST_CHECK_MESSAGE( bgSamples.find( sample ) == bgSamples.end(),
      "Foreground sample " << sample << " found in background samples" );
  }
}//test_set_displayed_time_ranges


BOOST_AUTO_TEST_CASE( test_get_displayed_time_ranges )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q67_Mobile.N42" );

  // Get the full time span
  json countRatesParams = json::object();
  json countRatesResult;
  BOOST_REQUIRE_NO_THROW( countRatesResult = registry.executeTool( "get_time_history_count_rates", countRatesParams, fixture.m_interspec ) );

  const json &slicesGD = countRatesResult["timeSlices"];
  BOOST_REQUIRE( slicesGD.size() >= 4 );

  const double firstTime = slicesGD[0]["startTime"].get<double>();
  const double lastTime = slicesGD[slicesGD.size() - 1]["startTime"].get<double>()
                         + slicesGD[slicesGD.size() - 1]["duration"].get<double>();
  const double totalSpan = lastTime - firstTime;

  // Set foreground and background ranges
  const double fgEnd = firstTime + (totalSpan / 3.0);
  const double bgStart = firstTime + (2.0 * totalSpan / 3.0);

  json setParams;
  setParams["foreground"] = to_string( firstTime ) + "-" + to_string( fgEnd );
  setParams["background"] = to_string( bgStart ) + "-" + to_string( lastTime );
  BOOST_REQUIRE_NO_THROW( registry.executeTool( "set_displayed_time_ranges", setParams, fixture.m_interspec ) );

  // Now get the displayed ranges
  json getParams = json::object();
  json getResult;
  BOOST_REQUIRE_NO_THROW( getResult = registry.executeTool( "get_displayed_time_ranges", getParams, fixture.m_interspec ) );

  BOOST_REQUIRE( getResult.is_object() );

  // Foreground ranges should be present and non-null
  BOOST_REQUIRE( getResult.contains( "foreground" ) );
  BOOST_CHECK( !getResult["foreground"].is_null() );

  // Background ranges should be present and non-null (since we set them)
  BOOST_REQUIRE( getResult.contains( "background" ) );
  BOOST_CHECK( !getResult["background"].is_null() );

  // Total time range info should be present
  BOOST_CHECK( getResult.contains( "total_time_range_start" ) );
  BOOST_CHECK( getResult.contains( "total_time_range_end" ) );
  BOOST_CHECK( getResult.contains( "total_num_samples" ) );
  BOOST_CHECK( getResult["total_num_samples"].get<int>() > 0 );
}//test_get_displayed_time_ranges


BOOST_AUTO_TEST_CASE( test_time_range_roundtrip )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  fixture.loadPassthroughFile( "Q8910_Helicopter.N42" );

  // Get the full time span
  json countRatesParams = json::object();
  json countRatesResult;
  BOOST_REQUIRE_NO_THROW( countRatesResult = registry.executeTool( "get_time_history_count_rates", countRatesParams, fixture.m_interspec ) );

  const json &slicesRT = countRatesResult["timeSlices"];
  BOOST_REQUIRE( slicesRT.size() >= 4 );

  const double firstTime = slicesRT[0]["startTime"].get<double>();
  const double lastTime = slicesRT[slicesRT.size() - 1]["startTime"].get<double>()
                         + slicesRT[slicesRT.size() - 1]["duration"].get<double>();
  const double totalSpan = lastTime - firstTime;

  // Step 1: Set ranges
  const double fgEnd = firstTime + (totalSpan / 4.0);
  const double bgStart = firstTime + (3.0 * totalSpan / 4.0);

  json setParams1;
  setParams1["foreground"] = to_string( firstTime ) + "-" + to_string( fgEnd );
  setParams1["background"] = to_string( bgStart ) + "-" + to_string( lastTime );

  json setResult1;
  BOOST_REQUIRE_NO_THROW( setResult1 = registry.executeTool( "set_displayed_time_ranges", setParams1, fixture.m_interspec ) );

  // Record the sample numbers after first set
  const std::set<int> fgSamples1 = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  const std::set<int> bgSamples1 = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Background );
  BOOST_REQUIRE( !fgSamples1.empty() );
  BOOST_REQUIRE( !bgSamples1.empty() );

  // Step 2: Get the ranges back
  json getParams = json::object();
  json getResult;
  BOOST_REQUIRE_NO_THROW( getResult = registry.executeTool( "get_displayed_time_ranges", getParams, fixture.m_interspec ) );

  BOOST_REQUIRE( getResult.contains( "foreground" ) );
  BOOST_REQUIRE( getResult.contains( "background" ) );
  BOOST_REQUIRE( !getResult["foreground"].is_null() );
  BOOST_REQUIRE( !getResult["background"].is_null() );

  const string fgRangesStr = getResult["foreground"].get<string>();
  const string bgRangesStr = getResult["background"].get<string>();
  BOOST_REQUIRE( !fgRangesStr.empty() );
  BOOST_REQUIRE( !bgRangesStr.empty() );

  // Step 3: Set ranges again using the returned values
  json setParams2;
  setParams2["foreground"] = fgRangesStr;
  setParams2["background"] = bgRangesStr;

  json setResult2;
  BOOST_REQUIRE_NO_THROW( setResult2 = registry.executeTool( "set_displayed_time_ranges", setParams2, fixture.m_interspec ) );

  // Step 4: Verify sample numbers are very similar (may differ by at most 1 sample due to
  //  floating-point rounding in the time range string representation)
  const std::set<int> fgSamples2 = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
  const std::set<int> bgSamples2 = fixture.m_interspec->displayedSamples( SpecUtils::SpectrumType::Background );

  // Allow at most 1 sample difference due to float rounding at boundaries
  const int fgSizeDiff = static_cast<int>( fgSamples2.size() ) - static_cast<int>( fgSamples1.size() );
  const int bgSizeDiff = static_cast<int>( bgSamples2.size() ) - static_cast<int>( bgSamples1.size() );
  BOOST_CHECK_LE( std::abs( fgSizeDiff ), 1 );
  BOOST_CHECK_LE( std::abs( bgSizeDiff ), 1 );

  // All original samples should still be present in the second set (roundtrip should not lose samples)
  for( const int s : fgSamples1 )
    BOOST_CHECK( fgSamples2.count( s ) > 0 );
  for( const int s : bgSamples1 )
    BOOST_CHECK( bgSamples2.count( s ) > 0 );
}//test_time_range_roundtrip


BOOST_AUTO_TEST_CASE( test_non_passthrough_error )
{
  TimeHistoryTestFixture fixture;
  const LlmTools::ToolRegistry &registry = fixture.llmToolRegistry();

  // Load a non-passthrough file
  fixture.loadNonPassthroughFile();

  // get_time_history_count_rates should throw for non-passthrough data
  json params = json::object();
  BOOST_CHECK_THROW(
    registry.executeTool( "get_time_history_count_rates", params, fixture.m_interspec ),
    std::runtime_error
  );

  // set_displayed_time_ranges should also throw for non-passthrough data
  json setParams;
  setParams["foreground"] = "0-5";
  BOOST_CHECK_THROW(
    registry.executeTool( "set_displayed_time_ranges", setParams, fixture.m_interspec ),
    std::runtime_error
  );

  // get_displayed_time_ranges should also throw for non-passthrough data
  json getParams = json::object();
  BOOST_CHECK_THROW(
    registry.executeTool( "get_displayed_time_ranges", getParams, fixture.m_interspec ),
    std::runtime_error
  );
}//test_non_passthrough_error


BOOST_AUTO_TEST_SUITE_END()
