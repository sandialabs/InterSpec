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

#define BOOST_TEST_MODULE AnalystChecks_suite
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/PhysicalUnits.h"
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
    m_env.reset( new Wt::Test::WTestEnvironment( applicationPath, configurationFile, Wt::Application ) );
    m_env->setAppRoot( wt_app_root );

    // Create the app
    m_app = new InterSpecApp( *m_env );

    m_update_lock.reset( new Wt::WApplication::UpdateLock(m_app) );

    // Get the InterSpec viewer instance
    m_interspec = m_app->viewer();
    BOOST_REQUIRE( m_interspec );

    // Load the Co56 test spectrum for escape peak tests
    const string test_file = SpecUtils::append_path( g_test_file_dir, "AnalystTests/escape_peak_check_Co56_Shielded_Fulcrum40h.n42" );
    BOOST_REQUIRE( SpecUtils::is_file(test_file) );

    m_interspec->userOpenFileFromFilesystem( test_file, test_file );
    Wt::WApplication::instance()->processEvents();

    shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    BOOST_REQUIRE( meas );
  }

  ~InterSpecTestFixture()
  {
    m_update_lock.reset();
    m_env.reset();
    m_interspec = nullptr;
  }
};


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
    BOOST_CHECK_CLOSE( result.potentialSingleEscapePeakEnergy, 2087.24 - 510.9989, 0.1 );
    BOOST_CHECK_CLOSE( result.potentialDoubleEscapePeakEnergy, 2087.24 - 2*510.9989, 0.1 );
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
  BOOST_CHECK_CLOSE( result.potentialSingleEscapePeakEnergy, test_energy - electron_rest_mass, 0.01 );
  BOOST_CHECK_CLOSE( result.potentialDoubleEscapePeakEnergy, test_energy - 2*electron_rest_mass, 0.01 );
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
  BOOST_REQUIRE( m_interspec->measurment(SpecUtils::SpectrumType::Foreground) );

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
                          << " but got " << AnalystChecks::to_string(info.sumType) );

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

        // Could be CascadeSum or Unknown depending on Ba-133 coincidence data
        const bool valid_type = (info.sumType == AnalystChecks::SumPeakType::CascadeSum) ||
                               (info.sumType == AnalystChecks::SumPeakType::Unknown);
        BOOST_CHECK_MESSAGE( valid_type,
                            "Test 5a: Expected CascadeSum or Unknown (no distance), got "
                            << AnalystChecks::to_string(info.sumType) );

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

        // Should get same result as no distance - CascadeSum or Unknown
        const bool valid_type = (info.sumType == AnalystChecks::SumPeakType::CascadeSum) ||
                               (info.sumType == AnalystChecks::SumPeakType::Unknown);
        BOOST_CHECK_MESSAGE( valid_type,
                            "Test 5b: Expected CascadeSum or Unknown (5 cm), got "
                            << AnalystChecks::to_string(info.sumType) );

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
        BOOST_CHECK_MESSAGE( info.sumType != AnalystChecks::SumPeakType::CascadeSum,
                            "Test 5c: At 100 cm, should NOT be CascadeSum, got "
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
