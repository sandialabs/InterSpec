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


//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DoseCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/DoseCalc.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DoseCalcWidget.h"
#include "InterSpec/GadrasSpecFunc.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;
using namespace boost::unit_test;


string find_data_dir( const string file )
{
  // The "InterSpec/data" directory that contains all the cross-sections, reactions, etc
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;
  
  std::string datadir;
  
  for( int i = 1; i < argc; ++i )
  {
    const std::string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }//for( int arg = 1; arg < argc; ++ arg )
  
  SpecUtils::ireplace_all( datadir, "%20", " " );
  
  // Search around a little for the data directory, if it wasnt specified
  if( datadir.empty() )
  {
    const string trial_paths[] = {
      "data",
      "../data",
      "../../data",
      "../../../data",
      "external_libs/SandiaDecay/",
      "../external_libs/SandiaDecay/",
      "../../external_libs/SandiaDecay/",
      "../../../external_libs/SandiaDecay/",
      "/Users/wcjohns/rad_ana/InterSpec/data/external_libs/SandiaDecay/",
      "/Users/wcjohns/rad_ana/InterSpec/data"
    };
    
    for( const auto &d : trial_paths )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, file) ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )
  
  return datadir;
}//void find_data_dir()


BOOST_AUTO_TEST_CASE( RuntimeSanityChecks )
{
  const string data_dir = find_data_dir( "GadrasContinuum.lib" );
  const string scatter_file = SpecUtils::append_path( data_dir, "GadrasContinuum.lib" );
  
  
  const string decay_xml = find_data_dir( "sandia.decay.xml" );
  const string sandia_decay_file = SpecUtils::append_path( decay_xml, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( !sandia_decay_file.empty() && SpecUtils::is_file(sandia_decay_file),
                        "Error finding sandia.decay.xml" );
  
  BOOST_REQUIRE_NO_THROW( DecayDataBaseServer::setDecayXmlFile( sandia_decay_file ) );
  
  
  // Set the static data directory so we have cross-sections, reactions, and all that
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( data_dir ) );
  
  
  BOOST_REQUIRE( SpecUtils::is_file(scatter_file) );
  
  std::unique_ptr<GadrasScatterTable> scatter_table;
  BOOST_REQUIRE_NO_THROW( scatter_table.reset( new GadrasScatterTable(scatter_file) ) );
  
  try
  {
    DoseCalcWidget::runtime_sanity_checks( scatter_table.get() );
  }catch( std::exception &e )
  {
    BOOST_CHECK_MESSAGE( false, e.what() );
  }
  
}//BOOST_AUTO_TEST_CASE( RuntimeSanityChecks )

