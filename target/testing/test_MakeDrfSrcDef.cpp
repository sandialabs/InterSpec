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
#include <string>
#include <vector>
#include <sstream>


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE testPhysicalUnits_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace boost::unit_test;

/* This file currently just tests reading `SrcLibLineInfo` info from files. */

namespace SetDataDirs
{
// The "InterSpec/data" directory that contains all the cross-sections, reactions, etc
std::string sm_static_data_dir = "";

// The sandia.decay.xml file in InterSpec/data is usually the minimized (~6MB) version of the file
//  without coincidences and stuff - we want to use the full ~30MB version of the file that is in
//  the SandiaDecay repository.
std::string sm_sandia_decay_file = "";

void set_data_dir()
{
  // We only need to initialize things once
  static bool s_have_set = false;
  if( s_have_set )
    return;
  
  s_have_set = true;
  
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  std::string datadir, test_file_dir;
  
  for( int i = 1; i < argc; ++i )
  {
    const std::string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      test_file_dir = arg.substr( 14 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  SpecUtils::ireplace_all( test_file_dir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data", "/Users/wcjohns/rad_ana/InterSpec/data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.reactiongamma.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  const std::string sandia_reaction_file = SpecUtils::append_path(datadir, "sandia.reactiongamma.xml");
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_reaction_file ), "sandia.reactiongamma.xml not at '" << sandia_reaction_file << "'" );
  
  // Set the static data directory so we have cross-sections, reactions, and all that
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  
  // Now find the full-version of sandia.decay.xml
  
  // What we actually want is the InterSpec code base-path - we'll search around a little for it
  const string decay_xml = "external_libs/SandiaDecay/sandia.decay.xml";
  string potential_code_dirs[] = {
    SpecUtils::append_path( test_file_dir, "../../../" ),
    "..", "../..", "../../..", "/Users/wcjohns/rad_ana/InterSpec/"
  };
  
  for( const auto &d : potential_code_dirs )
  {
    const string candidate = SpecUtils::append_path( d, decay_xml );
    if( SpecUtils::is_file(candidate) )
    {
      sm_sandia_decay_file = candidate;
      break;
    }
  }
  
  BOOST_REQUIRE_MESSAGE( !sm_sandia_decay_file.empty() && SpecUtils::is_file(sm_sandia_decay_file),
                        "Error finding the full version of sandia.decay.xml" );
  
  BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( sm_sandia_decay_file ) );
  
  sm_static_data_dir = datadir;
  
  // Make sure we can actually init the decay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing full SandiaDecayDataBase" );
  const SandiaDecay::Nuclide * const u238 = db->nuclide("U238");
  BOOST_REQUIRE_MESSAGE( u238, "Full SandiaDecayDataBase empty?" );
  
  // Now make sure we have the full database with x-rays and coincidences
  bool has_coincidences = false, has_xray = false;
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( u238, 1.0E-3 * SandiaDecay::curie );
  
  for( const SandiaDecay::NuclideActivityPair &nap : mixture.activity( 20*SandiaDecay::year ) )
  {
    for( const SandiaDecay::Transition *transition : nap.nuclide->decaysToChildren )
    {
      for( const SandiaDecay::RadParticle &particle : transition->products )
      {
        has_xray |= (particle.type == SandiaDecay::ProductType::XrayParticle);
        
        if( particle.type == SandiaDecay::GammaParticle )
          has_coincidences |= !particle.coincidences.empty();
      }//for( const SandiaDecay::RadParticle &particle : transition->products )
    }//for( const SandiaDecay::Transition *transition : nap.nuclide->decaysToChildren )
  }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
  
  BOOST_REQUIRE_MESSAGE( has_coincidences, "U238 decay in SandiaDecay didnt have coincidences" );
  BOOST_REQUIRE_MESSAGE( has_xray, "U238 decay in SandiaDecay didnt have x-rays" );
}//void set_data_dir()

}//namespace


BOOST_AUTO_TEST_CASE( sources_in_lib ) {
  SetDataDirs::set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( db );
  
  auto test_line = []( const string src_line, const SrcLibLineInfo &expected ){
    
    stringstream strm(src_line);
    vector<SrcLibLineInfo> info = SrcLibLineInfo::sources_in_lib( strm );
    BOOST_CHECK( info.size() == 1 );
    if( info.size() != 1 )
      return;
    
    const SrcLibLineInfo &src = info.at(0);
    
    BOOST_CHECK_CLOSE_FRACTION( src.m_activity, expected.m_activity, 1.0E-6 );
    BOOST_CHECK_EQUAL( src.m_nuclide, expected.m_nuclide );
    BOOST_CHECK_EQUAL( src.m_activity_date, expected.m_activity_date );
    BOOST_CHECK_EQUAL( src.m_source_name, expected.m_source_name );
    BOOST_CHECK_EQUAL( src.m_comments, expected.m_comments );
    BOOST_CHECK_EQUAL( src.m_line, src_line );
    
    BOOST_CHECK_CLOSE_FRACTION( src.m_activity_uncert, expected.m_activity_uncert, 1.0E-6 );
    BOOST_CHECK_CLOSE_FRACTION( src.m_distance, expected.m_distance, 1.0E-6 );
  };
  
  string src_line = "60Co_12916\xe2\x80\x82+5.03E+03   01-Jan-2024 ActivityUncert=1.0e3, Distance=5cm";
  SrcLibLineInfo expected;
  expected.m_activity = 5.03E+03 * PhysicalUnits::bq;
  expected.m_nuclide = db->nuclide( "Co60" );
  expected.m_activity_date = boost::posix_time::time_from_string( "2024-01-01 00:00:00.000" );
  expected.m_source_name = "60Co_12916";
  expected.m_comments = "ActivityUncert=1.0e3, Distance=5cm";
  expected.m_activity_uncert = 1.0e3 * PhysicalUnits::bq;
  expected.m_distance = 5.0 * PhysicalUnits::cm;
  test_line( src_line, expected );
  
  src_line = "60Co_12916 5.03E+03   01-Jan-2024 ActivityUncert=1.0e3";
  expected.m_distance = -1.0;
  expected.m_comments = "ActivityUncert=1.0e3";
  test_line( src_line, expected );
  
  src_line = "60Co_12916 5.03E+03   01-Jan-2024 Distance=5cm";
  expected.m_distance = 5.0 * PhysicalUnits::cm;
  expected.m_activity_uncert = -1.0;
  expected.m_comments = "Distance=5cm";
  test_line( src_line, expected );
  
  src_line = "60Co_12916 5.03E+03   01-Jan-2024";
  expected.m_distance = -1.0;
  expected.m_activity_uncert = -1.0;
  expected.m_comments = "";
  test_line( src_line, expected );
  
  
  src_line = "22NA_042713  3.600E+05  25-Apr-2020  SomeNote";
  expected.m_activity = 3.600E+05 * PhysicalUnits::bq;
  expected.m_nuclide = db->nuclide( "Na22" );
  expected.m_activity_date = boost::posix_time::time_from_string( "2020-04-25 00:00:00.000" );
  expected.m_source_name = "22NA_042713";
  expected.m_comments = "SomeNote";
  expected.m_activity_uncert = -1.0;
  expected.m_distance = -1.0;
  test_line( src_line, expected );
  
  
  src_line = "54MN_041711  3.511E+05  25-Apr-2020  AStorageLocation, Dist=24 cm, ActivityUncertainty=3.582E+04";
  expected.m_activity = 3.511E+05 * PhysicalUnits::bq;
  expected.m_nuclide = db->nuclide( "Mn54" );
  expected.m_activity_date = boost::posix_time::time_from_string( "2020-04-25 00:00:00.000" );
  expected.m_source_name = "54MN_041711";
  expected.m_comments = "AStorageLocation, Dist=24 cm, ActivityUncertainty=3.582E+04";
  expected.m_activity_uncert = 3.582E+04 * PhysicalUnits::bq;
  expected.m_distance = 24 * PhysicalUnits::cm;
  test_line( src_line, expected );
  
  
  src_line = "22NA_08251  4.107E+04  1-Aug-2019  SomePlace";
  expected.m_activity = 4.107E+04 * PhysicalUnits::bq;
  expected.m_nuclide = db->nuclide( "Na22" );
  expected.m_activity_date = boost::posix_time::time_from_string( "2019-08-01 00:00:00.000" );
  expected.m_source_name = "22NA_08251";
  expected.m_comments = "SomePlace";
  expected.m_activity_uncert = -1.0;
  expected.m_distance = -1.0;
  test_line( src_line, expected );
  
  
  // Test some invalid cases
  {
    stringstream strm("Invalid Line");
    vector<SrcLibLineInfo> info = SrcLibLineInfo::sources_in_lib( strm );
    BOOST_CHECK( info.size() == 0 );
  }
  
  {
    stringstream strm("22NA_08251\xe2\x80+4.107E+04  1-Aug-2019 \n"); //\xe2\x80 is missing \x82
    vector<SrcLibLineInfo> info = SrcLibLineInfo::sources_in_lib( strm );
    BOOST_CHECK( info.size() == 0 );
  }
  
  {
    stringstream strm("22NotANuc_08251  4.107E+04  1-Aug-2019  SomePlace\n");
    vector<SrcLibLineInfo> info = SrcLibLineInfo::sources_in_lib( strm );
    BOOST_CHECK( info.size() == 0 );
  }
  
  // Test multiple lines
  {
    stringstream strm( "   22NA_08251  4.107E+04  1-Aug-2019  SomePlace  \n"
                      "\n"
                      "  \n"
                      "54MN_041711  3.511E+05  25-Apr-2020  AStorageLocation, Dist=24 cm, ActivityUncertainty=3.582E+04\n"
                      "\n"
                      "  60Co_12916 5.03E+03   01-Jan-2024\n"
                      );
    vector<SrcLibLineInfo> info = SrcLibLineInfo::sources_in_lib( strm );
    BOOST_CHECK( info.size() == 3 );
  }
}//BOOST_AUTO_TEST_CASE( sources_in_lib )
