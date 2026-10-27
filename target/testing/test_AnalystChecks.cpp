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

#include <map>
#include <tuple>
#include <string>
#include <fstream>
#include <variant>
#include <iostream>
#include <algorithm>

#include <Wt/Utils.h>
#include <Wt/WApplication.h>
#include <Wt/Test/WTestEnvironment.h>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE AnalystChecks_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace boost::unit_test;

std::string g_test_file_dir;

// Set the static data directory
void set_data_dir()
{
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
  }

  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( g_test_file_dir, "%20", " " );

  // Search for the data directory if not specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }
  }

  // Search for test file directory if not specified
  if( g_test_file_dir.empty() )
  {
    const vector<string> candidate_test_dirs{
      "test_data",
      "../test_data",
      "target/testing/test_data",
      "../target/testing/test_data",
      "../../target/testing/test_data"
      , "../../../target/testing/test_data"
    };
    
    for( const auto &d : candidate_test_dirs )
    {
      if( SpecUtils::is_directory(d) )
      {
        g_test_file_dir = d;
        break;
      }
    }
  }

  if( datadir.empty() )
    throw runtime_error( "Could not find data directory" );

  if( g_test_file_dir.empty() )
    throw runtime_error( "Could not find test file directory" );

  // Initialize the decay database
  DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path(datadir, "sandia.decay.xml") );
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Failed to load decay database" );

  // Set other data directories
  InterSpec::setStaticDataDirectory( datadir );
}


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
    const std::string configurationFile = "";
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::EntryPointType::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) );

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    // Each session auto-saves its state into the user database (one per test process) and the next
    //  restores it, so without this a fixture can start with a previous session's spectrum
    //  already displayed - which makes `InterSpec::userOpenFile` treat the file we want as a
    //  possible background and only show a dialog no headless test answers.
    UserPreferences::setPreferenceValue<bool>( "AutoSaveSpectraToDb", false, m_interspec );

    // Load the Co56 test spectrum for escape peak tests
    const string test_file = SpecUtils::append_path( g_test_file_dir, "AnalystTests/escape_peak_check_Co56_Shielded_Fulcrum40h.n42" );
    BOOST_REQUIRE( SpecUtils::is_file(test_file) );

    m_interspec->userOpenFileFromFilesystem( test_file, test_file );
    Wt::WApplication::instance()->processEvents();

    shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    BOOST_REQUIRE( meas );
    BOOST_REQUIRE_MESSAGE( meas->filename() == test_file,
                          "Fixture spectrum was not loaded as the foreground: expected '"
                          << test_file << "', got '" << meas->filename() << "'" );
  }

  ~InterSpecTestFixture()
  {
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
  }
};


// A DRF the user picks once the spectrum is loaded brings its peak-fit prefs (e.g., a GADRAS
//  detector's peak shape) - unless the user set the spectrum's prefs in the GUI.
BOOST_FIXTURE_TEST_CASE( test_drf_change_adopts_peak_fit_prefs, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  BOOST_REQUIRE( meas );

  auto drf = make_shared<DetectorPeakResponse>();
  BOOST_REQUIRE_NO_THROW( drf->fromGadrasDirectory(
                            SpecUtils::append_path( g_test_file_dir, "gadras_detectors/CZT_1.5x2x2_text" ) ) );
  BOOST_REQUIRE( drf->peakFitDetPrefs() );

  m_interspec->detectorChanged().emit( drf );
  BOOST_REQUIRE( meas->peakFitDetPrefs() );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_peak_skew_type == PeakDef::SkewType::GadrasCZT );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_source == PeakFitDetPrefs::LoadingSource::FromDetectorPeakResponse );
  for( size_t i = 0; i < 6; ++i )
    BOOST_CHECK( meas->peakFitDetPrefs()->m_lower_energy_skew[i] == drf->peakFitDetPrefs()->m_lower_energy_skew[i] );

  // A DRF without prefs (or no DRF) doesnt leave the previous DRF's peak shape with the spectrum
  auto plain_drf = make_shared<DetectorPeakResponse>( *drf );
  plain_drf->setPeakFitDetPrefs( nullptr );
  m_interspec->detectorChanged().emit( plain_drf );
  BOOST_REQUIRE( meas->peakFitDetPrefs() );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_source != PeakFitDetPrefs::LoadingSource::FromDetectorPeakResponse );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_peak_skew_type != PeakDef::SkewType::GadrasCZT );

  m_interspec->detectorChanged().emit( drf );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_peak_skew_type == PeakDef::SkewType::GadrasCZT );
  m_interspec->detectorChanged().emit( nullptr );
  BOOST_CHECK( meas->peakFitDetPrefs()->m_peak_skew_type != PeakDef::SkewType::GadrasCZT );

  // Prefs the user set themselves are kept
  auto user_prefs = make_shared<PeakFitDetPrefs>();
  user_prefs->m_peak_skew_type = PeakDef::SkewType::Bortel;
  user_prefs->m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;
  meas->setPeakFitDetPrefs( user_prefs );
  m_interspec->detectorChanged().emit( drf );
  BOOST_CHECK( meas->peakFitDetPrefs() == user_prefs );
  m_interspec->detectorChanged().emit( plain_drf );
  BOOST_CHECK( meas->peakFitDetPrefs() == user_prefs );
}//test_drf_change_adopts_peak_fit_prefs


BOOST_FIXTURE_TEST_CASE( test_escape_peak_check_in_range_peaks, InterSpecTestFixture )
{
  // Test 1: Check 2087.24 keV is single escape of 2598.41 keV
  {
    AnalystChecks::EscapePeakCheckOptions options;
    options.energy = 2087.24;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

    BOOST_CHECK( result.parentPeak.has_value() );
    if( result.parentPeak.has_value() )
    {
      BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.parentPeak->escapeType), AnalystChecks::to_string(AnalystChecks::EscapePeakType::SingleEscape) );
      BOOST_CHECK_CLOSE( result.parentPeak->parentPeakEnergy, 2598.41, 1.0 ); // Within 1%
      BOOST_CHECK( result.parentPeak->sourceLabel.has_value() || result.parentPeak->userLabel.has_value() );
    }

    // Check potential energies are calculated correctly
    BOOST_REQUIRE( result.potentialSingleEscapePeakEnergy.has_value() );
    BOOST_REQUIRE( result.potentialDoubleEscapePeakEnergy.has_value() );
    BOOST_CHECK_CLOSE( *result.potentialSingleEscapePeakEnergy, 2087.24 - 510.9989, 0.1 );
    BOOST_CHECK_CLOSE( *result.potentialDoubleEscapePeakEnergy, 2087.24 - 2*510.9989, 0.1 );
    BOOST_CHECK( !result.belowPairProductionThreshold );
    BOOST_CHECK_CLOSE( result.potentialParentPeakSingleEscape, 2087.24 + 510.9989, 0.1 );
    BOOST_CHECK_CLOSE( result.potentialParentPeakDoubleEscape, 2087.24 + 2*510.9989, 0.1 );
  }

  // Test 2: Check 1576.38 keV is double escape of 2598.41 keV (NOT single escape of 2087.24)
  {
    AnalystChecks::EscapePeakCheckOptions options;
    options.energy = 1576.38;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

    BOOST_CHECK( result.parentPeak.has_value() );
    if( result.parentPeak.has_value() )
    {
      BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.parentPeak->escapeType), AnalystChecks::to_string(AnalystChecks::EscapePeakType::DoubleEscape) );
      BOOST_CHECK_CLOSE( result.parentPeak->parentPeakEnergy, 2598.41, 1.0 ); // Within 1%

      // Make sure it's NOT identified as single escape of 2087.24
      BOOST_CHECK( fabs(result.parentPeak->parentPeakEnergy - 2087.24) > 10.0 );
    }
  }
}


BOOST_FIXTURE_TEST_CASE( test_escape_peak_check_above_range_peaks, InterSpecTestFixture )
{
  // Test single escape peaks with parent above detector range
  const vector<double> single_escape_energies = {2690.87, 2742.35, 2761.96};
  const vector<string> expected_nuclide_labels = { "Co56 S.E. 3201.96 keV", "Co56 S.E. 3253.42 keV", "Co56 S.E. 3272.99 keV"};
  assert( single_escape_energies.size() == expected_nuclide_labels.size() );
  BOOST_REQUIRE_EQUAL( single_escape_energies.size(), expected_nuclide_labels.size() );

  for( size_t i = 0; i < single_escape_energies.size(); ++i )
  {
    const double energy = single_escape_energies[i];
      
    AnalystChecks::EscapePeakCheckOptions options;
    options.energy = energy;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

    BOOST_CHECK_MESSAGE( result.parentPeak.has_value(),
                         "Single escape peak at " << energy << " keV should have parent identified" );

    if( result.parentPeak.has_value() )
    {
      BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.parentPeak->escapeType), AnalystChecks::to_string(AnalystChecks::EscapePeakType::SingleEscape) );
      BOOST_CHECK_MESSAGE( result.parentPeak->parentPeakEnergy > 3000.0,
                           "Parent should be above 3000 keV for " << energy << " keV peak" );
      BOOST_CHECK( result.parentPeak->sourceLabel.has_value() );
      if( result.parentPeak->sourceLabel.has_value() )
      {
        BOOST_CHECK( !result.parentPeak->userLabel.has_value() );
        BOOST_CHECK_EQUAL( result.parentPeak->sourceLabel.value(), expected_nuclide_labels[i] );
      }
    }
  }
}


BOOST_FIXTURE_TEST_CASE( test_escape_peak_check_normal_peaks, InterSpecTestFixture )
{
  // Test a peak that is NOT an escape peak
  AnalystChecks::EscapePeakCheckOptions options;
  options.energy = 846.77;  // Co56 full-energy peak
  options.specType = SpecUtils::SpectrumType::Foreground;

  AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

  // Should not find a parent peak
  BOOST_CHECK( !result.parentPeak.has_value() );

  // But should still calculate potential energies
  BOOST_CHECK_CLOSE( result.potentialParentPeakSingleEscape, 846.77 + 510.9989, 0.1 );
  BOOST_CHECK_CLOSE( result.potentialParentPeakDoubleEscape, 846.77 + 2*510.9989, 0.1 );
}


BOOST_FIXTURE_TEST_CASE( test_escape_peak_check_search_windows, InterSpecTestFixture )
{
  AnalystChecks::EscapePeakCheckOptions options;
  options.energy = 1500.0;
  options.specType = SpecUtils::SpectrumType::Foreground;

  AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

  BOOST_CHECK( !result.parentPeak.has_value() );  // There is a 2015 keV peak we want to make sure we arent catching
    
  // Check that search windows are reasonable
  // Window should be max(1.0, min(0.5*fwhm, 15.0))
  BOOST_CHECK_GE( result.singleEscapeSearchWindow, 1.0 );
  BOOST_CHECK_LE( result.singleEscapeSearchWindow, 15.0 );

  BOOST_CHECK_GE( result.doubleEscapeSearchWindow, 1.0 );
  BOOST_CHECK_LE( result.doubleEscapeSearchWindow, 15.0 );

  BOOST_CHECK_GE( result.singleEscapeParentSearchWindow, 1.0 );
  BOOST_CHECK_LE( result.singleEscapeParentSearchWindow, 15.0 );

  BOOST_CHECK_GE( result.doubleEscapeParentSearchWindow, 1.0 );
  BOOST_CHECK_LE( result.doubleEscapeParentSearchWindow, 15.0 );
}


BOOST_FIXTURE_TEST_CASE( test_escape_peak_potential_energies, InterSpecTestFixture )
{
  // Test with arbitrary energy
  const double test_energy = 2000.0;
  const double electron_rest_mass = 510.9989;

  AnalystChecks::EscapePeakCheckOptions options;
  options.energy = test_energy;
  options.specType = SpecUtils::SpectrumType::Foreground;

  AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );

  // Verify potential energy calculations
  BOOST_REQUIRE( result.potentialSingleEscapePeakEnergy.has_value() );
  BOOST_REQUIRE( result.potentialDoubleEscapePeakEnergy.has_value() );
  BOOST_CHECK_CLOSE( *result.potentialSingleEscapePeakEnergy, test_energy - electron_rest_mass, 0.01 );
  BOOST_CHECK_CLOSE( *result.potentialDoubleEscapePeakEnergy, test_energy - 2*electron_rest_mass, 0.01 );
  BOOST_CHECK_CLOSE( result.potentialParentPeakSingleEscape, test_energy + electron_rest_mass, 0.01 );
  BOOST_CHECK_CLOSE( result.potentialParentPeakDoubleEscape, test_energy + 2*electron_rest_mass, 0.01 );
}


BOOST_FIXTURE_TEST_CASE( test_sum_peak_check, InterSpecTestFixture )
{
  // Load spectrum with summing peaks - not that this spectrum was incorrectly simulated and only has random-summing
  //  NOT cascade summing - however, we'll use it for the moment.
  const string test_file = SpecUtils::append_path( g_test_file_dir, "AnalystTests/Ba133_Cs137_RandomSummingOnly.n42" );

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(test_file),
                        "Test file not found: " << test_file );

  m_interspec->userOpenFileFromFilesystem( test_file, test_file );

  // Give the app a moment to load the file
  const shared_ptr<SpecMeas> sum_meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  BOOST_REQUIRE( sum_meas );
  // The fixture already loaded a different spectrum, so this open can be turned into a
  //  "is this a background?" prompt and silently skipped, leaving the fixture's spectrum on
  //  display - which would test the wrong data rather than fail.
  BOOST_REQUIRE_MESSAGE( sum_meas->filename() == test_file,
                        "Sum-peak spectrum was not loaded as the foreground: expected '"
                        << test_file << "', got '" << sum_meas->filename() << "'" );

  // Helper lambda to check sum peak results
  auto check_sum_peak = []( const double energy,
                            const AnalystChecks::SumPeakType expected_type,
                            const double expected_peak1_energy,
                            const double expected_peak2_energy,
                            InterSpec *interspec,
                            const string &description ) {
    AnalystChecks::SumPeakCheckOptions options;
    options.energy = energy;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, interspec );

    BOOST_CHECK_MESSAGE( result.sumPeakInfo.has_value(),
                        description << ": Peak at " << energy << " keV should be identified as a sum peak" );

    if( result.sumPeakInfo.has_value() )
    {
      const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;

      BOOST_CHECK_MESSAGE( info.sumType == expected_type,
                          description << ": Expected sum type " << AnalystChecks::to_string(expected_type)
                          << " but got " << AnalystChecks::to_string(info.sumType) << " (" << info.reason << ")" );

      BOOST_REQUIRE_MESSAGE( info.firstPeak && info.secondPeak,
                            description << ": Contributing peaks should be present" );

      if( info.firstPeak && info.secondPeak )
      {
        const double e1 = info.firstPeak->mean();
        const double e2 = info.secondPeak->mean();

        // Check if the energies match (in either order)
        const bool match1 = (fabs(e1 - expected_peak1_energy) < 2.0) && (fabs(e2 - expected_peak2_energy) < 2.0);
        const bool match2 = (fabs(e1 - expected_peak2_energy) < 2.0) && (fabs(e2 - expected_peak1_energy) < 2.0);

        BOOST_CHECK_MESSAGE( match1 || match2,
                            description << ": Expected peaks at " << expected_peak1_energy
                            << " and " << expected_peak2_energy << " keV, but got "
                            << e1 << " and " << e2 << " keV" );

        // Check that the sum is close to the target energy
        BOOST_CHECK_MESSAGE( fabs((e1 + e2) - energy) < 3.0,
                            description << ": Sum of contributing peaks (" << (e1 + e2)
                            << " keV) should match target energy (" << energy << " keV)" );
      }

      // Check that appropriate labels are generated
      const bool has_label = info.userLabel.has_value();
      BOOST_CHECK_MESSAGE( has_label,
                          description << ": Should have either userLabel" );
    }
  };

  // Test 1: 1323.2 keV peak should be random sum of 661.6 + 661.6 keV (Cs-137 self-sum)
  check_sum_peak( 1323.2, AnalystChecks::SumPeakType::RandomSum, 661.6, 661.6, m_interspec,
                  "Test 1: Cs-137 self-sum" );

  // Test 2: 712 keV peak is random sum of 356 + 356 keV (Ba-133 self-sum)
  check_sum_peak( 712.0, AnalystChecks::SumPeakType::RandomSum, 356.0, 356.0, m_interspec,
                  "Test 2: Ba-133 356 keV self-sum" );

  // Test 3: 1017 keV peak is random sum of 661.6 + 356 keV
  check_sum_peak( 1017.0, AnalystChecks::SumPeakType::RandomSum, 661.6, 356.0, m_interspec,
                  "Test 3: Cs-137 + Ba-133 sum" );

  // Test 4: 964 keV peak is random sum of 302.8 + 661.6 keV
  check_sum_peak( 964.0, AnalystChecks::SumPeakType::RandomSum, 302.8, 661.6, m_interspec,
                  "Test 4: Ba-133 302.8 + Cs-137 sum" );

  // Test 5: 437.1 keV peak is cascade-sum or unknown of 356 + 81 keV
  // Test with different distances to verify cascade detection is disabled at large distances
  // Note: This might be cascade-sum if Ba-133 has these as coincident gammas, or unknown
  {
    // Test 5a: No distance specified - cascade detection enabled
    {
      AnalystChecks::SumPeakCheckOptions options;
      options.energy = 437.1;
      options.specType = SpecUtils::SpectrumType::Foreground;
      // distance not specified

      AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

      BOOST_CHECK_MESSAGE( result.sumPeakInfo.has_value(),
                          "Test 5a: Peak at 437.1 keV should be identified as a sum peak (no distance)" );

      if( result.sumPeakInfo.has_value() )
      {
        const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;

        // 356 + 81 keV is a Ba-133 cascade, but this spectrum was simulated with random summing only
        //  (at 53% dead time), which accounts for nearly all of the 437 keV rate.
        const bool valid_type = (info.sumType == AnalystChecks::SumPeakType::RandomSum) ||
                                (info.sumType == AnalystChecks::SumPeakType::CascadeAndRandomSum);
        BOOST_CHECK_MESSAGE( valid_type,
                            "Test 5a: Expected RandomSum or CascadeAndRandomSum (no distance), got "
                            << AnalystChecks::to_string(info.sumType) << " - " << info.reason );

        if( info.firstPeak && info.secondPeak )
        {
          const double e1 = info.firstPeak->mean();
          const double e2 = info.secondPeak->mean();

          const bool match1 = (fabs(e1 - 356.0) < 2.0) && (fabs(e2 - 81.0) < 2.0);
          const bool match2 = (fabs(e1 - 81.0) < 2.0) && (fabs(e2 - 356.0) < 2.0);

          BOOST_CHECK_MESSAGE( match1 || match2,
                              "Test 5a: Expected peaks at 356 and 81 keV, got "
                              << e1 << " and " << e2 << " keV" );
        }
      }
    }

    // Test 5b: Small distance (5 cm) - cascade detection enabled
    {
      AnalystChecks::SumPeakCheckOptions options;
      options.energy = 437.1;
      options.specType = SpecUtils::SpectrumType::Foreground;
      options.distance = 5.0 * PhysicalUnits::cm;

      AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

      BOOST_CHECK_MESSAGE( result.sumPeakInfo.has_value(),
                          "Test 5b: Peak at 437.1 keV should be identified as a sum peak (5 cm)" );

      if( result.sumPeakInfo.has_value() )
      {
        const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;

        // Should get same result as no distance
        const bool valid_type = (info.sumType == AnalystChecks::SumPeakType::RandomSum) ||
                                (info.sumType == AnalystChecks::SumPeakType::CascadeAndRandomSum);
        BOOST_CHECK_MESSAGE( valid_type,
                            "Test 5b: Expected RandomSum or CascadeAndRandomSum (5 cm), got "
                            << AnalystChecks::to_string(info.sumType) << " - " << info.reason );

        if( info.firstPeak && info.secondPeak )
        {
          const double e1 = info.firstPeak->mean();
          const double e2 = info.secondPeak->mean();

          const bool match1 = (fabs(e1 - 356.0) < 2.0) && (fabs(e2 - 81.0) < 2.0);
          const bool match2 = (fabs(e1 - 81.0) < 2.0) && (fabs(e2 - 356.0) < 2.0);

          BOOST_CHECK_MESSAGE( match1 || match2,
                              "Test 5b: Expected peaks at 356 and 81 keV, got "
                              << e1 << " and " << e2 << " keV" );
        }
      }
    }

    // Test 5c: Large distance (100 cm) - cascade detection disabled, should be RandomSum if detected
    {
      AnalystChecks::SumPeakCheckOptions options;
      options.energy = 437.1;
      options.specType = SpecUtils::SpectrumType::Foreground;
      options.distance = 100.0 * PhysicalUnits::cm;

      AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

      // At 100 cm, cascade detection is disabled
      // The peak might still be detected as a sum, but should be RandomSum (not CascadeSum)
      if( result.sumPeakInfo.has_value() )
      {
        const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;

        // Should NOT be CascadeSum at this distance
        BOOST_CHECK_MESSAGE( (info.sumType != AnalystChecks::SumPeakType::CascadeSum)
                             && (info.sumType != AnalystChecks::SumPeakType::CascadeAndRandomSum),
                            "Test 5c: At 100 cm, should NOT be a cascade sum, got "
                            << AnalystChecks::to_string(info.sumType) );

        // If detected as a sum peak, should be RandomSum or Unknown
        const bool valid_type = (info.sumType == AnalystChecks::SumPeakType::RandomSum) ||
                               (info.sumType == AnalystChecks::SumPeakType::Unknown);
        BOOST_CHECK_MESSAGE( valid_type,
                            "Test 5c: At 100 cm, expected RandomSum or Unknown, got "
                            << AnalystChecks::to_string(info.sumType) );
      }
    }
  }

  // Test 6: 658.15 keV peak is random sum of 302.85 + 356.02 keV
  check_sum_peak( 658.15, AnalystChecks::SumPeakType::RandomSum, 302.85, 356.02, m_interspec,
                  "Test 6: Ba-133 302.85 + 356.02 sum" );

  // Test 7: 302.8 keV is NOT a sum peak
  {
    AnalystChecks::SumPeakCheckOptions options;
    options.energy = 302.8;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

    BOOST_CHECK_MESSAGE( !result.sumPeakInfo.has_value(),
                        "Test 7: Peak at 302.8 keV should NOT be identified as a sum peak" );
  }

  // Test 8: 661.64 keV is NOT a sum peak (Cs-137 full-energy peak)
  {
    AnalystChecks::SumPeakCheckOptions options;
    options.energy = 661.64;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

    BOOST_CHECK_MESSAGE( !result.sumPeakInfo.has_value(),
                        "Test 8: Peak at 661.64 keV (Cs-137 full-energy) should NOT be identified as a sum peak" );
  }

  // Test 9: Check search window is reasonable
  {
    AnalystChecks::SumPeakCheckOptions options;
    options.energy = 1000.0;
    options.specType = SpecUtils::SpectrumType::Foreground;

    AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );

    // Search window should be max(1.0, min(0.5*fwhm, 15.0))
    BOOST_CHECK_GE( result.searchWindow, 1.0 );
    BOOST_CHECK_LE( result.searchWindow, 15.0 );
  }
}


// ---------------------------------------------------------------------------------------------
// beta_continuum_check tests.
//
// Test PCFs (copied from the peak_fit_accuracy_inject_compact benchmark) each contain the
// source foreground record, a matched-live-time "Background" record, and a 10x-longer
// "Background" record.
// ---------------------------------------------------------------------------------------------

namespace
{
  /** Loads `filename` (relative to AnalystTests/) as both displayed Foreground and Background: the
   non-"Background" record becomes the foreground (or the short background record, when
   `background_as_foreground`), and the longest-live-time "Background" record the background. */
  void load_pcf_fore_and_back( InterSpec * const interspec, const std::string &filename,
                               const bool background_as_foreground = false )
  {
    const string path = SpecUtils::append_path( g_test_file_dir, "AnalystTests/" + filename );
    BOOST_REQUIRE( SpecUtils::is_file(path) );

    shared_ptr<SpecMeas> fore_meas = make_shared<SpecMeas>();
    BOOST_REQUIRE( fore_meas->load_file( path, SpecUtils::ParserType::Auto, path ) );
    shared_ptr<SpecMeas> back_meas = make_shared<SpecMeas>();
    BOOST_REQUIRE( back_meas->load_file( path, SpecUtils::ParserType::Auto, path ) );

    int fore_sample = -999, short_bkg_sample = -999, long_bkg_sample = -999;
    float short_bkg_lt = (std::numeric_limits<float>::max)(), long_bkg_lt = -1.0f;
    for( size_t i = 0; i < fore_meas->num_measurements(); ++i )
    {
      const shared_ptr<const SpecUtils::Measurement> m = fore_meas->measurement_at_index( static_cast<int>(i) );
      BOOST_REQUIRE( !!m );
      if( SpecUtils::icontains( m->title(), "Background" ) )
      {
        if( m->live_time() > long_bkg_lt )
        {
          long_bkg_lt = m->live_time();
          long_bkg_sample = m->sample_number();
        }
        if( m->live_time() < short_bkg_lt )
        {
          short_bkg_lt = m->live_time();
          short_bkg_sample = m->sample_number();
        }
      }else if( fore_sample == -999 )
      {
        fore_sample = m->sample_number();
      }
    }//for( measurements )

    BOOST_REQUIRE( fore_sample != -999 );
    BOOST_REQUIRE( long_bkg_sample != -999 );

    const int use_fore_sample = background_as_foreground ? short_bkg_sample : fore_sample;
    BOOST_REQUIRE( use_fore_sample != -999 );

    interspec->setSpectrum( fore_meas, {use_fore_sample}, SpecUtils::SpectrumType::Foreground, {} );
    Wt::WApplication::instance()->processEvents();
    interspec->setSpectrum( back_meas, {long_bkg_sample}, SpecUtils::SpectrumType::Background, {} );
    Wt::WApplication::instance()->processEvents();
  }//load_pcf_fore_and_back(...)
}//namespace


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_sr90_hpge, InterSpecTestFixture )
{
  load_pcf_fore_and_back( m_interspec, "BetaCheck/detx_sr90_unshielded.pcf" );

  const AnalystChecks::BetaContinuumCheckOptions options;
  const AnalystChecks::BetaContinuumCheckStatus result
      = AnalystChecks::beta_continuum_check( options, m_interspec );

  BOOST_TEST_MESSAGE( "Sr90 HPGe: verdict=" << AnalystChecks::to_string(result.verdict)
                      << ", elevSigma=" << result.continuumElevationSigma
                      << ", eTerm=" << (result.terminationEnergy ? *result.terminationEnergy : -1.0)
                      << ", eChar=" << (result.characteristicEnergy ? *result.characteristicEnergy : -1.0)
                      << ", nGammaPeaks=" << result.elevatedGammaPeaks.size() );

  BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.verdict),
                     AnalystChecks::to_string(AnalystChecks::BetaContinuumVerdict::BremLike) );
  BOOST_CHECK( result.elevatedGammaPeaks.empty() );
  BOOST_REQUIRE( result.terminationEnergy.has_value() );
  BOOST_CHECK_GT( *result.terminationEnergy, 500.0 );
  BOOST_CHECK_LT( *result.terminationEnergy, 2450.0 );

  bool have_sr90 = false;
  for( const AnalystChecks::BetaContinuumCandidateNuclide &candidate : result.candidateNuclides )
    have_sr90 = (have_sr90 || (candidate.nuclide == "Sr90"));
  BOOST_CHECK_MESSAGE( have_sr90, "Sr90 should be among the candidate nuclides" );
}//test_beta_continuum_sr90_hpge


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_sr90_nai, InterSpecTestFixture )
{
  load_pcf_fore_and_back( m_interspec, "BetaCheck/nai3x3_sr90_unshielded.pcf" );

  const AnalystChecks::BetaContinuumCheckOptions options;
  const AnalystChecks::BetaContinuumCheckStatus result
      = AnalystChecks::beta_continuum_check( options, m_interspec );

  BOOST_TEST_MESSAGE( "Sr90 NaI: verdict=" << AnalystChecks::to_string(result.verdict)
                      << ", elevSigma=" << result.continuumElevationSigma
                      << ", eTerm=" << (result.terminationEnergy ? *result.terminationEnergy : -1.0)
                      << ", eChar=" << (result.characteristicEnergy ? *result.characteristicEnergy : -1.0)
                      << ", nGammaPeaks=" << result.elevatedGammaPeaks.size() );

  BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.verdict),
                     AnalystChecks::to_string(AnalystChecks::BetaContinuumVerdict::BremLike) );
  BOOST_CHECK( result.elevatedGammaPeaks.empty() );
  BOOST_REQUIRE( result.terminationEnergy.has_value() );

  bool have_sr90 = false;
  for( const AnalystChecks::BetaContinuumCandidateNuclide &candidate : result.candidateNuclides )
    have_sr90 = (have_sr90 || (candidate.nuclide == "Sr90"));
  BOOST_CHECK_MESSAGE( have_sr90, "Sr90 should be among the candidate nuclides" );
}//test_beta_continuum_sr90_nai


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_cs137_not_brem, InterSpecTestFixture )
{
  load_pcf_fore_and_back( m_interspec, "BetaCheck/detx_cs137_shielded.pcf" );

  const AnalystChecks::BetaContinuumCheckOptions options;
  const AnalystChecks::BetaContinuumCheckStatus result
      = AnalystChecks::beta_continuum_check( options, m_interspec );

  BOOST_TEST_MESSAGE( "Cs137 shielded HPGe: verdict=" << AnalystChecks::to_string(result.verdict)
                      << ", nGammaPeaks=" << result.elevatedGammaPeaks.size() );

  BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.verdict),
                     AnalystChecks::to_string(AnalystChecks::BetaContinuumVerdict::NotBremLike) );

  bool have_662 = false;
  for( const AnalystChecks::BetaContinuumPeakInfo &peak : result.elevatedGammaPeaks )
    have_662 = (have_662 || (fabs(peak.energy - 661.66) < 3.0));
  BOOST_CHECK_MESSAGE( have_662, "661.66 keV should be among the elevated gamma peaks" );
}//test_beta_continuum_cs137_not_brem


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_background_as_foreground, InterSpecTestFixture )
{
  // The short "Background" record as foreground, the long one as background: nothing is elevated.
  load_pcf_fore_and_back( m_interspec, "BetaCheck/detx_cs137_shielded.pcf", true );

  const AnalystChecks::BetaContinuumCheckOptions options;
  const AnalystChecks::BetaContinuumCheckStatus result
      = AnalystChecks::beta_continuum_check( options, m_interspec );

  BOOST_TEST_MESSAGE( "Background-as-foreground: verdict=" << AnalystChecks::to_string(result.verdict)
                      << ", elevSigma=" << result.continuumElevationSigma );

  BOOST_CHECK_EQUAL( AnalystChecks::to_string(result.verdict),
                     AnalystChecks::to_string(AnalystChecks::BetaContinuumVerdict::NoContinuumElevation) );
  BOOST_CHECK_LT( result.continuumElevationSigma, 8.0 );
}//test_beta_continuum_background_as_foreground


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_tl204_xrays, InterSpecTestFixture )
{
  // Tl-204's EC branch gives Hg K x-rays (~69-83 keV); these are reported but must not veto.
  load_pcf_fore_and_back( m_interspec, "BetaCheck/nai3x3_tl204_unshielded.pcf" );

  const AnalystChecks::BetaContinuumCheckOptions options;
  const AnalystChecks::BetaContinuumCheckStatus result
      = AnalystChecks::beta_continuum_check( options, m_interspec );

  BOOST_TEST_MESSAGE( "Tl204 NaI: verdict=" << AnalystChecks::to_string(result.verdict)
                      << ", eTerm=" << (result.terminationEnergy ? *result.terminationEnergy : -1.0)
                      << ", nXrayPeaks=" << result.xrayRegionPeaks.size()
                      << ", nGammaPeaks=" << result.elevatedGammaPeaks.size() );

  const bool brem_or_ambiguous = (result.verdict == AnalystChecks::BetaContinuumVerdict::BremLike)
                                 || (result.verdict == AnalystChecks::BetaContinuumVerdict::Ambiguous);
  BOOST_CHECK_MESSAGE( brem_or_ambiguous, "Tl204 should grade BremLike or Ambiguous, got "
                       << AnalystChecks::to_string(result.verdict) );

  BOOST_REQUIRE( result.terminationEnergy.has_value() );
  BOOST_CHECK_LT( *result.terminationEnergy, 950.0 ); // Tl-204 endpoint is 763.7 keV
}//test_beta_continuum_tl204_xrays


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_no_background, InterSpecTestFixture )
{
  // Direct test of the worker-safe overload: null background must yield BackgroundNotLoaded,
  // and the check must be deterministic (identical results for identical inputs).
  const string path = SpecUtils::append_path( g_test_file_dir,
                                    "AnalystTests/BetaCheck/detx_sr90_unshielded.pcf" );
  BOOST_REQUIRE( SpecUtils::is_file(path) );

  shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( path, SpecUtils::ParserType::Auto, path ) );
  BOOST_REQUIRE( meas->num_measurements() > 0 );

  AnalystChecks::BetaContinuumCheckInput input;
  input.foreground = meas->measurement_at_index( 0 );
  input.expectedFwhm = []( const double energy ) -> float {
    return static_cast<float>( 2.0 + 0.03*std::sqrt(energy) );
  };

  const AnalystChecks::BetaContinuumCheckStatus first = AnalystChecks::beta_continuum_check( input );
  BOOST_CHECK_EQUAL( AnalystChecks::to_string(first.verdict),
                     AnalystChecks::to_string(AnalystChecks::BetaContinuumVerdict::BackgroundNotLoaded) );
  BOOST_CHECK( !first.warnings.empty() );

  const AnalystChecks::BetaContinuumCheckStatus second = AnalystChecks::beta_continuum_check( input );
  BOOST_CHECK_EQUAL( AnalystChecks::to_string(first.verdict), AnalystChecks::to_string(second.verdict) );
}//test_beta_continuum_no_background


// ---------------------------------------------------------------------------------------------
// Data-driven beta_continuum_check harness.
//
// Reads test_data/AnalystTests/beta_check_cases.csv (see that file for the column format) and
// runs the check on every row, so new brem spectra / shielded-gamma negatives / x-ray-machine
// spectra can be added without touching C++.
// ---------------------------------------------------------------------------------------------

namespace
{
  /** Resolve a manifest sample selector to a sample number of `meas`.
   'auto' = first record not titled "Background"; 'auto-bkg' = the "Background" record with the
   longest (prefer_long) or shortest live time, falling back to longest/shortest-live-time record
   overall; otherwise the integer sample number. */
  int resolve_manifest_sample( const std::shared_ptr<SpecMeas> &meas, const std::string &selector,
                               const bool prefer_long )
  {
    if( (selector != "auto") && (selector != "auto-bkg") )
      return std::stoi( selector );

    int best_sample = -999, best_any_sample = -999;
    float best_lt = prefer_long ? -1.0f : (std::numeric_limits<float>::max)();
    float best_any_lt = best_lt;

    for( size_t i = 0; i < meas->num_measurements(); ++i )
    {
      const shared_ptr<const SpecUtils::Measurement> m = meas->measurement_at_index( static_cast<int>(i) );
      if( !m )
        continue;
      const bool is_bkg = SpecUtils::icontains( m->title(), "Background" );
      const bool better_any = prefer_long ? (m->live_time() > best_any_lt) : (m->live_time() < best_any_lt);

      if( selector == "auto" )
      {
        if( !is_bkg )
          return m->sample_number();
      }else if( is_bkg )
      {
        const bool better = prefer_long ? (m->live_time() > best_lt) : (m->live_time() < best_lt);
        if( better )
        {
          best_lt = m->live_time();
          best_sample = m->sample_number();
        }
      }

      if( better_any )
      {
        best_any_lt = m->live_time();
        best_any_sample = m->sample_number();
      }
    }//for( measurements )

    if( best_sample != -999 )
      return best_sample;
    if( best_any_sample != -999 )
      return best_any_sample;  //'auto'/'auto-bkg' with no matching titles: any record
    throw std::runtime_error( "No usable record for selector '" + selector + "'" );
  }//resolve_manifest_sample(...)


  std::vector<std::string> split_csv_row( const std::string &line )
  {
    std::vector<std::string> fields;
    std::string field;
    for( const char c : line )
    {
      if( c == ',' )
      {
        fields.push_back( field );
        field.clear();
      }else
      {
        field += c;
      }
    }
    fields.push_back( field );
    for( std::string &f : fields )
      SpecUtils::trim( f );
    return fields;
  }//split_csv_row(...)
}//namespace


BOOST_FIXTURE_TEST_CASE( test_beta_continuum_case_manifest, InterSpecTestFixture )
{
  const string manifest_path = SpecUtils::append_path( g_test_file_dir,
                                                       "AnalystTests/beta_check_cases.csv" );
  BOOST_REQUIRE( SpecUtils::is_file(manifest_path) );

  std::ifstream manifest( manifest_path.c_str() );
  BOOST_REQUIRE( manifest.is_open() );

  size_t num_cases = 0, num_skipped = 0;
  string line;
  while( std::getline( manifest, line ) )
  {
    SpecUtils::trim( line );
    if( line.empty() || (line[0] == '#') )
      continue;

    const vector<string> fields = split_csv_row( line );
    BOOST_REQUIRE_MESSAGE( fields.size() >= 8, "Manifest row needs 8 columns: " << line );

    const string &fore_file = fields[0], &fore_sel = fields[1];
    const string &back_file = fields[2], &back_sel = fields[3];
    const string &allowed = fields[4], &min_term = fields[5], &max_term = fields[6];
    const string &note = fields[7];

    const auto resolve_path = [&]( const string &f ) -> string {
      if( SpecUtils::is_file(f) )
        return f;
      return SpecUtils::append_path( g_test_file_dir, "AnalystTests/" + f );
    };

    const string fore_path = resolve_path( fore_file );
    const string back_path = resolve_path( back_file );

    // Rows may reference spectra outside this repository (e.g. the benchmark dataset repo), which
    //  are not available on CI or a fresh checkout - skip those rows instead of failing.
    if( !SpecUtils::is_file(fore_path) || !SpecUtils::is_file(back_path) )
    {
      num_skipped += 1;
      BOOST_TEST_MESSAGE( "[beta manifest] skipping '" << note << "': spectrum file not available ("
                          << (SpecUtils::is_file(fore_path) ? back_file : fore_file) << ")" );
      continue;
    }

    num_cases += 1;

    shared_ptr<SpecMeas> fore_meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( fore_meas->load_file( fore_path, SpecUtils::ParserType::Auto, fore_file ),
                           "Failed to load " << fore_file );
    shared_ptr<SpecMeas> back_meas = make_shared<SpecMeas>();
    BOOST_REQUIRE_MESSAGE( back_meas->load_file( back_path, SpecUtils::ParserType::Auto, back_file ),
                           "Failed to load " << back_file );

    // Foreground 'auto-bkg' means the shortest-live-time background record (null tests).
    const int fore_sample = resolve_manifest_sample( fore_meas, fore_sel, false );
    const int back_sample = resolve_manifest_sample( back_meas, back_sel, true );

    m_interspec->setSpectrum( fore_meas, {fore_sample}, SpecUtils::SpectrumType::Foreground, {} );
    Wt::WApplication::instance()->processEvents();
    m_interspec->setSpectrum( back_meas, {back_sample}, SpecUtils::SpectrumType::Background, {} );
    Wt::WApplication::instance()->processEvents();

    const AnalystChecks::BetaContinuumCheckOptions options;
    AnalystChecks::BetaContinuumCheckStatus result;
    BOOST_REQUIRE_NO_THROW( result = AnalystChecks::beta_continuum_check( options, m_interspec ) );

    const string verdict = AnalystChecks::to_string( result.verdict );
    BOOST_TEST_MESSAGE( "[beta manifest] " << note << " (" << fore_file << "): verdict=" << verdict
        << " elevSigma=" << result.continuumElevationSigma
        << " eTerm=" << (result.terminationEnergy ? *result.terminationEnergy : -1.0)
        << " eChar=" << (result.characteristicEnergy ? *result.characteristicEnergy : -1.0)
        << " nGammaPk=" << result.elevatedGammaPeaks.size()
        << " nXrayPk=" << result.xrayRegionPeaks.size() );

    // Verdict must be one of the '|'-separated allowed values.
    vector<string> allowed_verdicts;
    SpecUtils::split( allowed_verdicts, allowed, "|" );
    const bool verdict_ok = (std::find( begin(allowed_verdicts), end(allowed_verdicts), verdict )
                             != end(allowed_verdicts));
    BOOST_CHECK_MESSAGE( verdict_ok, note << ": verdict '" << verdict
                         << "' not in allowed set '" << allowed << "'" );

    if( (min_term != "-") && !min_term.empty() )
    {
      BOOST_CHECK_MESSAGE( result.terminationEnergy && (*result.terminationEnergy >= std::stod(min_term)),
                           note << ": terminationEnergy below " << min_term << " keV" );
    }
    if( (max_term != "-") && !max_term.empty() )
    {
      BOOST_CHECK_MESSAGE( result.terminationEnergy && (*result.terminationEnergy <= std::stod(max_term)),
                           note << ": terminationEnergy above " << max_term << " keV" );
    }

    // Optional 9th column: a nuclide that must appear among the consistency candidates.
    if( (fields.size() >= 9) && !fields[8].empty() && (fields[8] != "-") )
    {
      bool have_candidate = false;
      for( const AnalystChecks::BetaContinuumCandidateNuclide &candidate : result.candidateNuclides )
        have_candidate = (have_candidate || (candidate.nuclide == fields[8]));
      BOOST_CHECK_MESSAGE( have_candidate, note << ": '" << fields[8]
                           << "' missing from ConsistentPureBetaNuclides" );
    }
  }//while( manifest rows )

  BOOST_TEST_MESSAGE( "[beta manifest] ran " << num_cases << " cases (" << num_skipped
                      << " skipped for missing spectrum files)" );
  BOOST_CHECK_GT( num_cases, 0 );
}//test_beta_continuum_case_manifest


namespace
{
  /** Loads a spectrum from `test_data/AnalystTests/` as the foreground, and optionally another as
   the background.  (Not `userOpenFileFromFilesystem`, which may ask if a second file is a background.) */
  shared_ptr<SpecMeas> load_analyst_test_spectrum( InterSpec * const interspec, const std::string &filename,
                                                   const std::string &background_filename = "" )
  {
    const auto load = [interspec]( const std::string &name, const SpecUtils::SpectrumType type ) -> shared_ptr<SpecMeas> {
      const string path = SpecUtils::append_path( g_test_file_dir, "AnalystTests/" + name );
      BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Missing test file " << path );

      shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
      BOOST_REQUIRE( meas->load_file( path, SpecUtils::ParserType::Auto, path ) );
      BOOST_REQUIRE( !meas->sample_numbers().empty() );

      interspec->setSpectrum( meas, { *meas->sample_numbers().begin() }, type, {} );
      Wt::WApplication::instance()->processEvents();
      return meas;
    };

    shared_ptr<SpecMeas> foreground = load( filename, SpecUtils::SpectrumType::Foreground );
    if( !background_filename.empty() )
      load( background_filename, SpecUtils::SpectrumType::Background );
    return foreground;
  }//load_analyst_test_spectrum(...)


  // NSDD "Flavor of the Month" January 2025 drill (LA-UR-24-33222): an unshielded Pb-203 HPGe
  //  spectrum (8192 channels, 6.6% dead time), with Tl K x-rays, 279.2, 401.3 and 680.5 keV lines.
  const char * const sm_q567_file = "NSDD_2025Jan_Q567_Pb203_HPGe.CHN";
}//namespace


BOOST_AUTO_TEST_CASE( test_characteristic_gammas_of_gammaless_nuclides )
{
  set_data_dir();

  // SandiaDecay lists a zero-rate 511 keV entry for every nuclide; for one with no photons at all
  //  that used to come out as its "characteristic" line.
  for( const char * const nuc : { "H3", "P32", "C14", "S35" } )
  {
    const vector<float> gammas = AnalystChecks::get_characteristic_gammas( nuc );
    BOOST_CHECK_MESSAGE( gammas.empty(), nuc << " should have no characteristic gammas, got " << gammas.size()
                         << (gammas.empty() ? string() : (" starting with " + std::to_string(gammas[0]))) );
  }

  const vector<float> cs137 = AnalystChecks::get_characteristic_gammas( "Cs137" );
  BOOST_REQUIRE( !cs137.empty() );
  BOOST_CHECK( std::any_of( begin(cs137), end(cs137), []( const float e ){ return fabs(e - 661.66f) < 0.5f; } ) );

  // A source's only gamma is characteristic of it, however weak (Po210: 803 keV in 1E-5 of decays)
  const vector<float> po210 = AnalystChecks::get_characteristic_gammas( "Po210" );
  BOOST_CHECK_MESSAGE( std::any_of( begin(po210), end(po210), []( const float e ){ return fabs(e - 803.1f) < 0.5f; } ),
                       "Po210 should have its 803 keV line as characteristic" );

  // The importance helper on photons that are never emitted
  BOOST_CHECK( AnalystChecks::cacl_estimated_gamma_importance( { {511.0f, 0.0f} } ).empty() );
  const vector<tuple<float,float,float>> one = AnalystChecks::cacl_estimated_gamma_importance( { {511.0f, 0.0f}, {100.0f, 2.0f} } );
  BOOST_REQUIRE_EQUAL( one.size(), 1 );
  BOOST_CHECK_CLOSE( std::get<2>(one[0]), 1.0f, 1.0E-4 );
}//test_characteristic_gammas_of_gammaless_nuclides


BOOST_FIXTURE_TEST_CASE( test_positron_emitters_for_511, InterSpecTestFixture )
{
  const vector<AnalystChecks::AnnihilationSource> sources = AnalystChecks::positron_emitters_for_annihilation( m_interspec );
  BOOST_REQUIRE( !sources.empty() );

  map<string,double> yields;
  for( size_t i = 0; i < sources.size(); ++i )
  {
    BOOST_REQUIRE( sources[i].nuclide );
    BOOST_CHECK_GE( sources[i].photons511PerDecay, 0.1 );
    if( i )
      BOOST_CHECK_GE( sources[i-1].photons511PerDecay, sources[i].photons511PerDecay );
    yields[sources[i].nuclide->symbol] = sources[i].photons511PerDecay;
  }

  BOOST_REQUIRE( yields.count( "Na22" ) );
  BOOST_CHECK_CLOSE( yields["Na22"], 1.8, 5.0 );
  BOOST_REQUIRE( yields.count( "F18" ) );
  BOOST_CHECK_CLOSE( yields["F18"], 1.93, 5.0 );
  BOOST_CHECK( yields.count( "Ga68" ) || yields.count( "Ge68" ) );
  BOOST_CHECK( !yields.count( "H3" ) );
  BOOST_CHECK( !yields.count( "Cs137" ) );
  BOOST_CHECK( !yields.count( "Cl36" ) );  // beta+ in ~1E-5 of decays

  // A nuclide may have 511.0 keV (annihilation) as a characteristic line only if it is a real
  //  positron emitter - not, e.g., Cl36, whose trace annihilation photons are its only gamma-ray-energy
  //  photons.  (Rh106's 511.86 keV is a real gamma, so is not in this window.)
  const vector<AnalystChecks::AnnihilationSource> any_positrons
                                 = AnalystChecks::positron_emitters_for_annihilation( m_interspec, 1.0E-12 );
  map<string,double> any_yields;
  for( const AnalystChecks::AnnihilationSource &src : any_positrons )
    any_yields[src.nuclide->symbol] = src.photons511PerDecay;

  size_t num_checked = 0;
  for( const auto &item : AnalystChecks::get_nuclides_with_characteristics_in_energy_range( 510.95, 511.05, m_interspec ) )
  {
    if( !std::holds_alternative<const SandiaDecay::Nuclide *>(item) )
      continue;
    const string symbol = std::get<const SandiaDecay::Nuclide *>(item)->symbol;
    num_checked += 1;
    BOOST_CHECK_MESSAGE( any_yields.count( symbol ) && (any_yields[symbol] >= AnalystChecks::sm_min_annihilation_photons_per_decay),
                         symbol << " has 511 keV as a characteristic line, with only "
                         << (any_yields.count( symbol ) ? any_yields[symbol] : 0.0) << " annihilation photons per decay" );
  }
  BOOST_CHECK_GT( num_checked, 2 );

  const vector<float> cl36 = AnalystChecks::get_characteristic_gammas( "Cl36" );
  BOOST_CHECK( std::none_of( begin(cl36), end(cl36), []( const float e ){ return fabs(e - 511.0f) < 1.0f; } ) );
}//test_positron_emitters_for_511


BOOST_FIXTURE_TEST_CASE( test_expected_fwhm_checks_drf, InterSpecTestFixture )
{
  const shared_ptr<SpecMeas> meas = load_analyst_test_spectrum( m_interspec, sm_q567_file );
  meas->setDetector( nullptr );

  // The spectrum's own peaks: about 0.55 keV at 70 keV, and 0.8 keV at 279 keV.
  string source;
  const float fwhm_no_drf = AnalystChecks::get_expected_fwhm( 143.7, m_interspec, &source );
  BOOST_TEST_MESSAGE( "Q567 FWHM at 143.7 keV, no DRF: " << fwhm_no_drf << " (" << source << ")" );
  BOOST_CHECK_GT( fwhm_no_drf, 0.5f );
  BOOST_CHECK_LT( fwhm_no_drf, 0.85f );
  BOOST_CHECK( SpecUtils::icontains( source, "fit to" ) );

  // A generic DRF, whose ~1.6 keV FWHM is twice what this detector does, is overruled by the peaks.
  auto generic = make_shared<DetectorPeakResponse>( "Generic HPGe", "" );
  generic->setIntrinsicEfficiencyFormula( "1.0", 5.0*PhysicalUnits::cm, PhysicalUnits::keV, 0.0f, 0.0f,
                                          DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  generic->setFwhmCoefficients( { 1.54f, 0.264f, 0.33f }, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn );
  meas->setDetector( generic );

  const float fwhm_generic = AnalystChecks::get_expected_fwhm( 143.7, m_interspec, &source );
  BOOST_TEST_MESSAGE( "Q567 FWHM at 143.7 keV, generic DRF: " << fwhm_generic << " (" << source << ")" );
  BOOST_CHECK_LT( fwhm_generic, 0.9f );
  BOOST_CHECK( SpecUtils::icontains( source, "not for this detector" ) );

  // A DRF that agrees with the peaks is used as is.
  auto matching = make_shared<DetectorPeakResponse>( "Matching HPGe", "" );
  matching->setIntrinsicEfficiencyFormula( "1.0", 5.0*PhysicalUnits::cm, PhysicalUnits::keV, 0.0f, 0.0f,
                                           DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
  matching->setFwhmCoefficients( { 0.25f, 1.6f }, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
  meas->setDetector( matching );

  const float fwhm_matching = AnalystChecks::get_expected_fwhm( 143.7, m_interspec, &source );
  BOOST_TEST_MESSAGE( "Q567 FWHM at 143.7 keV, matching DRF: " << fwhm_matching << " (" << source << ")" );
  BOOST_CHECK_EQUAL( source, "detector response function" );
  BOOST_CHECK_CLOSE( fwhm_matching, matching->peakResolutionFWHM( 143.7f ), 1.0E-3 );

  BOOST_CHECK_THROW( AnalystChecks::get_expected_fwhm( -100.0, m_interspec ), std::exception );
}//test_expected_fwhm_checks_drf


BOOST_FIXTURE_TEST_CASE( test_compton_kinematics_without_candidate, InterSpecTestFixture )
{
  load_analyst_test_spectrum( m_interspec, sm_q567_file );

  // 133-146 keV is the backscatter / Compton-edge plateau of the 279.2 keV line (E_bs = 133.4,
  //  E_edge = 145.8 keV), a broad feature the peak search does not fit.
  for( const double energy : { 133.0, 140.0, 146.0 } )
  {
    AnalystChecks::ComptonPeakCheckOptions options;
    options.energy = energy;
    const AnalystChecks::ComptonPeakCheckStatus result = AnalystChecks::compton_peak_check( options, m_interspec );

    BOOST_TEST_MESSAGE( "Compton " << energy << " keV: " << AnalystChecks::to_string(result.type) << ", "
                        << AnalystChecks::to_string(result.confidence) << " - " << result.reason );
    BOOST_CHECK( !result.reason.empty() );
    BOOST_REQUIRE_MESSAGE( !result.kinematicParents.empty(), "No Compton parent for " << energy << " keV" );

    const AnalystChecks::ComptonKinematicParent &best = result.kinematicParents.front();
    BOOST_CHECK_SMALL( best.parentEnergy - 279.15, 0.5 );
    BOOST_CHECK_CLOSE( best.backscatterEnergy, 133.4, 0.2 );
    BOOST_CHECK_CLOSE( best.comptonEdgeEnergy, 145.8, 0.2 );
    BOOST_CHECK( best.feature != AnalystChecks::ComptonScatterType::None );
    BOOST_CHECK_GE( static_cast<int>(best.confidence), static_cast<int>(AnalystChecks::ComptonScatterConfidence::Likely) );
    if( result.noCandidatePeak )
    {
      BOOST_CHECK( result.type != AnalystChecks::ComptonScatterType::None );
      BOOST_CHECK_GE( static_cast<int>(result.confidence), static_cast<int>(AnalystChecks::ComptonScatterConfidence::Likely) );
    }
  }//for( energies )

  // Nothing scatters to 600 keV: it is above every backscatter region, and at no Compton edge.
  AnalystChecks::ComptonPeakCheckOptions options;
  options.energy = 600.0;
  const AnalystChecks::ComptonPeakCheckStatus none = AnalystChecks::compton_peak_check( options, m_interspec );
  BOOST_CHECK( none.kinematicParents.empty() );
  BOOST_CHECK( none.noCandidatePeak );
  BOOST_CHECK( !none.reason.empty() );
}//test_compton_kinematics_without_candidate


BOOST_FIXTURE_TEST_CASE( test_escape_check_thresholds_and_xrays, InterSpecTestFixture )
{
  load_analyst_test_spectrum( m_interspec, sm_q567_file );

  AnalystChecks::EscapePeakCheckOptions options;

  // 170 keV is below the pair-production threshold; its would-be S.E. parent, 681 keV, is a real
  //  Pb-203 line, but is also below the threshold.
  {
    options.energy = 170.0;
    const AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );
    BOOST_CHECK( result.belowPairProductionThreshold );
    BOOST_CHECK( !result.potentialSingleEscapePeakEnergy.has_value() );
    BOOST_CHECK( !result.potentialDoubleEscapePeakEnergy.has_value() );
    BOOST_CHECK( !result.parentPeak.has_value() );
    BOOST_REQUIRE( result.peakAtPotentialSingleEscapeParent.has_value() );
    BOOST_CHECK_SMALL( result.peakAtPotentialSingleEscapeParent->peak->mean() - 680.6, 1.5 );
    BOOST_REQUIRE( !result.notes.empty() );
    bool mentions_threshold = false;
    for( const string &note : result.notes )
      mentions_threshold = mentions_threshold || SpecUtils::icontains( note, "1022" );
    BOOST_CHECK( mentions_threshold );
    BOOST_CHECK( SpecUtils::icontains( result.detectorMaterialAssumed, "Ge" ) );
  }

  // 391.4 keV is where a Ge K-alpha escape of the 401.3 keV line would be - at 401 keV the escape is
  //  far too rare to make a peak, so it is listed but not plausible.
  {
    options.energy = 391.44;
    const AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );
    bool found_ge_ka = false;
    for( const AnalystChecks::XrayEscapeInfo &x : result.xrayEscapes )
    {
      if( (x.element == "Ge") && (x.line == "Ka") )
      {
        found_ge_ka = true;
        BOOST_REQUIRE( x.parentPeak.has_value() );
        BOOST_CHECK_SMALL( x.parentPeak->peak->mean() - 401.32, 1.0 );
        BOOST_CHECK( !x.plausible );
        BOOST_CHECK_LT( x.maxEscapeFraction, 0.01 );
      }
    }//for( x-ray escapes )
    BOOST_CHECK( found_ge_ka );
  }

  // 350.02 keV (a Tl K-alpha + 279 keV sum) is 10.98 keV - a Ge K-beta escape - below the 361.6 keV
  //  sum peak, but is bigger than it, so it is not that peak's escape.
  {
    options.energy = 350.02;
    const AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, m_interspec );
    BOOST_CHECK( result.belowPairProductionThreshold );
    BOOST_REQUIRE( !result.xrayEscapes.empty() );
    for( const AnalystChecks::XrayEscapeInfo &x : result.xrayEscapes )
      BOOST_CHECK_MESSAGE( !x.plausible, "350.02 keV should not be a plausible " << x.element << " " << x.line << " escape" );
  }
}//test_escape_check_thresholds_and_xrays


BOOST_FIXTURE_TEST_CASE( test_sum_check_xray_cascade_and_pileup, InterSpecTestFixture )
{
  load_analyst_test_spectrum( m_interspec, sm_q567_file );

  // The sum checks need the peaks attributed to Pb203.
  AnalystChecks::FitPeaksForNuclideOptions fit_options;
  fit_options.sources = { "Pb203" };
  const AnalystChecks::FitPeaksForNuclideStatus fit_status = AnalystChecks::fit_peaks_for_nuclides( fit_options, m_interspec );
  BOOST_REQUIRE( !fit_status.fitPeaks.empty() );

  const auto check = [this]( const double energy ) -> AnalystChecks::SumPeakCheckStatus {
    AnalystChecks::SumPeakCheckOptions options;
    options.energy = energy;
    const AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );
    if( result.sumPeakInfo )
      BOOST_TEST_MESSAGE( "Sum " << energy << " keV: " << AnalystChecks::to_string(result.sumPeakInfo->sumType)
                          << " - " << result.sumPeakInfo->reason << "; tau=" << (result.pileupTau ? *result.pileupTau : -1.0)
                          << " (" << result.pileupTauSource << ")" );
    return result;
  };

  // A line can only sum with itself randomly: one Pb-203 decay emits at most one 279 keV gamma.
  {
    const AnalystChecks::SumPeakCheckStatus result = check( 558.34 );
    BOOST_REQUIRE( result.sumPeakInfo.has_value() );
    BOOST_CHECK( result.sumPeakInfo->sumType == AnalystChecks::SumPeakType::RandomSum );
    BOOST_CHECK( !result.sumPeakInfo->reason.empty() );
  }

  // Tl K x-rays (from the electron capture) are emitted with the 279 keV gamma.
  for( const double energy : { 350.02, 352.06, 361.61 } )
  {
    const AnalystChecks::SumPeakCheckStatus result = check( energy );
    BOOST_REQUIRE_MESSAGE( result.sumPeakInfo.has_value(), "No sum found at " << energy << " keV" );
    const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;
    BOOST_CHECK_MESSAGE( (info.sumType == AnalystChecks::SumPeakType::CascadeSum)
                         || (info.sumType == AnalystChecks::SumPeakType::CascadeAndRandomSum),
                         energy << " keV: expected a cascade sum, got " << AnalystChecks::to_string(info.sumType)
                         << " - " << info.reason );
    BOOST_CHECK( info.cascadeNuclide.has_value() && (*info.cascadeNuclide == "Pb203") );
    BOOST_CHECK( info.coincidenceType.has_value() && (*info.coincidenceType == "x-ray-gamma") );
    BOOST_CHECK( info.coincidenceFraction.has_value() && (*info.coincidenceFraction > 0.01) );
  }//for( Tl K x-ray + 279 keV sums )

  // Two K x-rays in one decay - the EC K-capture's, and the 279 keV transition's K-conversion's -
  //  sum to 143.7 keV (Ka2 70.83 + Ka1 72.87).  There is no fitted peak there, so no observed rate.
  {
    const AnalystChecks::SumPeakCheckStatus result = check( 143.7 );
    BOOST_REQUIRE( result.sumPeakInfo.has_value() );
    const AnalystChecks::SumPeakInfo &info = *result.sumPeakInfo;
    BOOST_CHECK_MESSAGE( info.sumType == AnalystChecks::SumPeakType::CascadeSum,
                         "143.7 keV: expected CascadeSum, got " << AnalystChecks::to_string(info.sumType) << " - " << info.reason );
    BOOST_CHECK( info.coincidenceType.has_value() && (*info.coincidenceType == "x-ray-x-ray") );
    BOOST_REQUIRE( info.coincidenceFraction.has_value() );
    // ~0.78 EC-K x 0.13 K-conversions, x 2 x 0.49 x 0.29 (Ka1, Ka2 per vacancy)
    BOOST_CHECK_GT( *info.coincidenceFraction, 0.015 );
    BOOST_CHECK_LT( *info.coincidenceFraction, 0.05 );
    BOOST_CHECK( !info.observedSumCps.has_value() );
    BOOST_CHECK( info.expectedRandomSumCps.has_value() );
    BOOST_CHECK( SpecUtils::icontains( info.reason, "no fitted peak" ) );
    BOOST_CHECK( !SpecUtils::icontains( info.reason, "could not be" ) );
  }

  // The pile-up constant from the 558 keV self-sum: 0.46 cps / (909 cps)^2 ~ 5.6E-7 s.
  {
    const AnalystChecks::SumPeakCheckStatus result = check( 352.06 );
    BOOST_REQUIRE( result.pileupTau.has_value() );
    BOOST_CHECK( !result.pileupTauIsUpperBound );
    BOOST_CHECK_GT( *result.pileupTau, 0.5*5.6E-7 );
    BOOST_CHECK_LT( *result.pileupTau, 2.0*5.6E-7 );
    BOOST_REQUIRE( result.sumPeakInfo && result.sumPeakInfo->expectedRandomSumCps );
    BOOST_CHECK_GT( *result.sumPeakInfo->expectedRandomSumCps, 0.3 );
    BOOST_CHECK_LT( *result.sumPeakInfo->expectedRandomSumCps, 3.0 );
  }
}//test_sum_check_xray_cascade_and_pileup


BOOST_FIXTURE_TEST_CASE( test_fit_user_peak_probe_reports_why, InterSpecTestFixture )
{
  load_analyst_test_spectrum( m_interspec, sm_q567_file );

  // On the 279 keV backscatter plateau there is no peak the automated thresholds accept; a probe
  //  says why, and gives what the data show there.
  AnalystChecks::FitPeakOptions options;
  options.energy = 143.7;
  options.doNotAddToAnalysisPeaks = true;
  const AnalystChecks::FitPeakStatus probe = AnalystChecks::fit_user_peak( options, m_interspec );
  BOOST_CHECK( !probe.significant );
  BOOST_REQUIRE( probe.diagnostic.has_value() );
  BOOST_TEST_MESSAGE( "Probe at 143.7 keV: " << probe.diagnostic->reason << "; fixed fit "
                      << probe.diagnostic->fixedFitArea << " +- " << probe.diagnostic->fixedFitAreaUncert );
  BOOST_CHECK( !probe.diagnostic->reason.empty() );
  BOOST_CHECK( probe.diagnostic->fixedFitValid );
  BOOST_CHECK_GT( probe.diagnostic->fixedFitAreaUncert, 0.0 );
  BOOST_CHECK_GT( probe.diagnostic->expectedFwhm, 0.4 );
  BOOST_CHECK_LT( probe.diagnostic->expectedFwhm, 1.0 );
  BOOST_CHECK( m_interspec->measurment( SpecUtils::SpectrumType::Foreground )->peaks(
                 m_interspec->displayedSamples( SpecUtils::SpectrumType::Foreground ) )->empty() );

  // Without the probe flag it is an error, whose message says why.
  options.doNotAddToAnalysisPeaks = false;
  try
  {
    AnalystChecks::fit_user_peak( options, m_interspec );
    BOOST_CHECK_MESSAGE( false, "Expected adding a peak at 143.7 keV to fail" );
  }catch( std::exception &e )
  {
    const string msg = e.what();
    BOOST_CHECK( SpecUtils::icontains( msg, "Could not fit peak" ) );
    BOOST_CHECK( SpecUtils::icontains( msg, "DoNotAddToAnalysisPeaks" ) );
  }

  // A real peak is still just fit.
  options.energy = 279.15;
  options.doNotAddToAnalysisPeaks = true;
  const AnalystChecks::FitPeakStatus real = AnalystChecks::fit_user_peak( options, m_interspec );
  BOOST_CHECK( real.significant );
  BOOST_REQUIRE( real.fitPeak );
  BOOST_CHECK_SMALL( real.fitPeak->mean() - 279.15, 0.3 );
}//test_fit_user_peak_probe_reports_why


// A candidate peak near a strong line's backscatter peak, when a weaker (here spurious) peak matches
//  its energy even better, is judged against the parent that could have produced it.  From the NSDD
//  January 2025 drill: Pb-203 behind ~3 cm Pb on NaI, where add-peaks-for-source fits a 5-sigma
//  "279 keV" peak (transmission ~4E-7, so not real), and 190 keV is the 681 keV line's backscatter
//  peak (E_bs = 185.8 keV).
BOOST_FIXTURE_TEST_CASE( test_compton_candidate_uses_plausible_parent, InterSpecTestFixture )
{
  load_analyst_test_spectrum( m_interspec, "NSDD_2025Jan_Q4_Pb203_NaI.CHN", "NSDD_2025Jan_Q4_Pb203_NaI_Background.CHN" );

  AnalystChecks::FitPeaksForNuclideOptions fit_options;
  fit_options.sources = { "Pb203" };
  const AnalystChecks::FitPeaksForNuclideStatus fit_status = AnalystChecks::fit_peaks_for_nuclides( fit_options, m_interspec );
  BOOST_REQUIRE( !fit_status.fitPeaks.empty() );

  AnalystChecks::ComptonPeakCheckOptions options;
  options.energy = 190.0;
  const AnalystChecks::ComptonPeakCheckStatus result = AnalystChecks::compton_peak_check( options, m_interspec );

  BOOST_CHECK( !result.reason.empty() );
  BOOST_TEST_MESSAGE( "Compton 190 keV (Q4 NaI): " << AnalystChecks::to_string(result.type) << ", "
                      << AnalystChecks::to_string(result.confidence) << ", noCandidate=" << result.noCandidatePeak
                      << ", parent=" << (result.parentInfo ? result.parentInfo->parentEnergy : -1.0)
                      << " - " << result.reason );

  BOOST_REQUIRE( !result.kinematicParents.empty() );
  BOOST_CHECK_SMALL( result.kinematicParents.front().parentEnergy - 680.5, 5.0 );

  // The verdict is about the 681 keV line, not the spurious 279 keV peak the candidate outweighs
  BOOST_REQUIRE_MESSAGE( !result.noCandidatePeak, "Expected a fitted candidate peak near 190 keV" );
  BOOST_REQUIRE( result.parentInfo.has_value() );
  BOOST_CHECK_SMALL( result.parentInfo->parentEnergy - 680.5, 5.0 );
  BOOST_CHECK( result.type != AnalystChecks::ComptonScatterType::None );
  BOOST_CHECK( result.confidence != AnalystChecks::ComptonScatterConfidence::No );
  BOOST_CHECK( !SpecUtils::icontains( result.reason, "larger than its kinematic parent" ) );

  // A verdict below VeryLikely names the criteria that held it there
  if( result.confidence != AnalystChecks::ComptonScatterConfidence::VeryLikely )
    BOOST_CHECK_MESSAGE( SpecUtils::icontains( result.reason, "rather than VeryLikely because of:" ), result.reason );
}//test_compton_candidate_uses_plausible_parent


// A real line at twice a line's energy is not that line's random sum, and must not be taken as a
//  self-sum peak when measuring the pile-up rate: Eu-152 has a real 688.67 keV line at 2 x 344.28 keV.
BOOST_FIXTURE_TEST_CASE( test_sum_check_real_line_at_twice_energy, InterSpecTestFixture )
{
  const string path = SpecUtils::append_path( g_test_file_dir,
                          "../analysis_tests/AEGIS_Eu152_surface_contamination.n42_20230622T113239.276178.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(path), "Missing test file " << path );
  shared_ptr<SpecMeas> meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( path, SpecUtils::ParserType::Auto, path ) );
  BOOST_REQUIRE( !meas->sample_numbers().empty() );
  m_interspec->setSpectrum( meas, meas->displayedSampleNumbers().empty() ? set<int>{ *meas->sample_numbers().begin() }
                                                                          : meas->displayedSampleNumbers(),
                            SpecUtils::SpectrumType::Foreground, {} );
  Wt::WApplication::instance()->processEvents();

  AnalystChecks::FitPeaksForNuclideOptions fit_options;
  fit_options.sources = { "Eu152" };
  AnalystChecks::fit_peaks_for_nuclides( fit_options, m_interspec );

  AnalystChecks::SumPeakCheckOptions options;
  options.energy = 688.67;
  const AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, m_interspec );
  if( result.sumPeakInfo )
  {
    BOOST_TEST_MESSAGE( "Sum 688.67 keV (Eu152): " << AnalystChecks::to_string(result.sumPeakInfo->sumType)
                        << " - " << result.sumPeakInfo->reason );
    BOOST_CHECK_MESSAGE( result.sumPeakInfo->sumType != AnalystChecks::SumPeakType::RandomSum,
                         "The real Eu152 688.67 keV line was called a random sum: " << result.sumPeakInfo->reason );
  }

  // Whatever tau is, it is not "measured" from 688.67 = 2 x 344.28, so not above the dead-time bound
  const shared_ptr<const SpecUtils::Measurement> spec = m_interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  BOOST_REQUIRE( spec );
  if( result.pileupTau && !result.pileupTauIsUpperBound
      && (spec->real_time() > spec->live_time()) && (spec->gamma_count_sum() > 0.0) )
  {
    BOOST_TEST_MESSAGE( "tau = " << *result.pileupTau << " (" << result.pileupTauSource << ")" );
    BOOST_CHECK_LE( *result.pileupTau, (spec->real_time() - spec->live_time()) / spec->gamma_count_sum() );
    BOOST_CHECK( !SpecUtils::icontains( result.pileupTauSource, "2 x 344" ) );
  }
}//test_sum_check_real_line_at_twice_energy
