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
