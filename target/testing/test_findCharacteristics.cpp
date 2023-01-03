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
#include <map>
#include <deque>
#include <tuple>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <cstdlib>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE findCharacteristics_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


using namespace std;
using namespace boost::unit_test;

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
    cout << "Arg " << i << ": '" << argv[i] << "'" << endl;
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
  
  const string required_data_file = "findCharacteristics/202204_example_problem_1.n42";
  if( g_test_file_dir.empty() )
  {
    for( const auto &d : { "test_data", "../test_data", "../../test_data", "/Users/wcjohns/rad_ana/InterSpec/target/testing/test_data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, required_data_file) ) )
      {
        g_test_file_dir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }
  
  const string sandia_deacay_file = SpecUtils::append_path(datadir, "sandia.decay.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_deacay_file ), "sandia.decay.xml not at '" << sandia_deacay_file << "'" );
  
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );
  
  const string required_data_path = SpecUtils::append_path(g_test_file_dir, required_data_file);
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( required_data_path ), "'" << required_data_file << "' not at '" << required_data_path << "'" );
}//void set_data_dir()



BOOST_AUTO_TEST_CASE( testFindCharacteristics )
{
  set_data_dir();
  
  auto peak = make_shared<const PeakDef>( 661.66, 10, 1.0E5 );
  vector<string> char_nucs;
  IsotopeId::findCharacteristics( char_nucs, peak );
  BOOST_CHECK_MESSAGE( (char_nucs.size() > 1) && (char_nucs[0] == "Cs137 661.66 keV"),
                     "Error Finding Cs137 characteristics" );
  char_nucs.clear();
  
  
  peak = make_shared<const PeakDef>( 1001.0, 10, 1.0E5 );
  IsotopeId::findCharacteristics( char_nucs, peak );
  BOOST_CHECK_MESSAGE( (char_nucs.size() > 1) && (char_nucs[0] == "U238 1000.99 keV"),
                      "Error Finding U238 characteristics" );
  char_nucs.clear();
  
  
  peak = make_shared<const PeakDef>( 2103, 10, 1.0E5 );
  IsotopeId::findCharacteristics( char_nucs, peak );
  
  auto check_results_contain = []( const std::string &val, const vector<string> &character ) {
    auto pos = std::find( begin(character), end(character), val );
    BOOST_CHECK_MESSAGE( pos != end(character),
                        "IsotopeId::findCharacteristics results didnt include '" << val << "'" );
  };
  
  vector<string> required_results = { "I132 S.E. 2614.50 keV", "Th232 S.E. 2614.53 keV",
    "Th228 S.E. 2614.53 keV", "U232 S.E. 2614.53 keV", "Tl208 S.E. 2614.53 keV",
    "Pb212 S.E. 2614.53 keV", "Pb(n,g) S.E. 2614.53 keV", "Pb(n,n) S.E. 2614.68 keV",
    "U238 2101.90 keV", "Co62m 2104.60 keV", "Y92 2105.60 keV"
  };
  for( const string &required : required_results )
    check_results_contain( required, char_nucs );
  
  char_nucs.clear();
}//BOOST_AUTO_TEST_CASE( testFindCharacteristics ) 


BOOST_AUTO_TEST_CASE( testFindCandidates )
{
  // This is only a sanity check of the IsotopeId::findCandidates(...) function.
  //  We should include S.E., D.E., reactions, x-rays, and low-resolution data.
  set_data_dir();
  
  const string data_dir = InterSpec::staticDataDirectory();
  
  const string spec_file = "findCharacteristics/202204_example_problem_1.n42";
  const string spec_path = SpecUtils::append_path(g_test_file_dir, spec_file);
  BOOST_REQUIRE( SpecUtils::is_file(spec_path) );
  
  SpecMeas meas;
  BOOST_REQUIRE_MESSAGE( meas.load_N42_file(spec_path), "Failed to open required file: '" << spec_path << "'" );
  
  BOOST_REQUIRE( meas.num_measurements() == 2 );
  shared_ptr<const SpecUtils::Measurement> foreground = meas.measurements()[0];
  shared_ptr<const SpecUtils::Measurement> background = meas.measurements()[1];
  
  BOOST_REQUIRE( foreground && background );
  BOOST_REQUIRE( foreground->source_type() == SpecUtils::SourceType::Foreground );
  BOOST_REQUIRE( background->source_type() == SpecUtils::SourceType::Background );
  
  shared_ptr<DetectorPeakResponse> detector = meas.detector();
  BOOST_REQUIRE( detector && detector->isValid() );
  shared_ptr<deque<shared_ptr<const PeakDef>>> allpeaks = meas.peaks( {foreground->sample_number()} );
  BOOST_REQUIRE( allpeaks && !allpeaks->empty() );
  
  // We will strip all peaks of analyst assigned nuclide IDs
  shared_ptr<const PeakDef> cs137_peak, co60_peak, eu152_peak;
  auto allpeaks_no_id = make_shared<deque<shared_ptr<const PeakDef>>>();
  for( auto peak : *allpeaks )
  {
    auto new_peak = make_shared<PeakDef>( *peak );
    new_peak->clearSources();
    new_peak->setCandidateNuclides( {} );
    allpeaks_no_id->push_back( new_peak );
    if( fabs(new_peak->mean() - 344.3) < 0.5 )
      eu152_peak = new_peak;
    else if( fabs(new_peak->mean() - 661.4) < 0.5 )
      cs137_peak = new_peak;
    else if( fabs(new_peak->mean() - 1173.4) < 0.5 )
      co60_peak = new_peak;
  }//
  
  BOOST_REQUIRE( co60_peak && eu152_peak && cs137_peak );
  
  auto to_csv = []( const vector<string> &input ) -> string {
    string val = "{";
    for( size_t i = 0; i < input.size(); ++i )
      val += (i ? ", '" : "'") + input[i] + "'";
    return val + "}";
  };
  
  vector<string> suggestednucs;
  IsotopeId::findCandidates( suggestednucs, eu152_peak, allpeaks_no_id, detector, foreground );
  BOOST_CHECK_MESSAGE( !suggestednucs.empty() && (suggestednucs[0] == "Eu152 344.28 keV"),
                      "The top recommendation for the 344.3 keV should have been 'Eu152 344.28 keV'; suggestions were "
                      << to_csv(suggestednucs) );
  
  suggestednucs.clear();
  IsotopeId::findCandidates( suggestednucs, cs137_peak, allpeaks_no_id, detector, foreground );
  BOOST_CHECK_MESSAGE( std::count( begin(suggestednucs), end(suggestednucs), "Cs137 661.66 keV"),
                      "The 661 keV peak didnt have the expected 'Cs137 661.66 keV'; suggestions were "
                      << to_csv(suggestednucs) );
  
  suggestednucs.clear();
  IsotopeId::findCandidates( suggestednucs, co60_peak, allpeaks_no_id, detector, foreground );
  BOOST_CHECK_MESSAGE( !suggestednucs.empty() && (suggestednucs[0] == "Co60 1173.23 keV"),
                      "The top recommendation for the 1173 keV should have been 'Co60 1173.23 keV'; suggestions were "
                      << to_csv(suggestednucs) );
}//BOOST_AUTO_TEST_CASE( testFindCandidates )
