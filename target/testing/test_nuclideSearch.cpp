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
#define BOOST_TEST_MODULE nuclideSearch_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"


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



BOOST_AUTO_TEST_CASE( testNuclideSearch )
{
  set_data_dir();

  shared_ptr<IsotopeSearchByEnergyModel::SearchWorkingSpace> orig_workingspace
      = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>();

  orig_workingspace->energies = { 121.886, 244.983, 344.534 };
  orig_workingspace->windows = { 2.03, 2.04, 2.1 };
  //orig_workingspace->detector_response_function = std::shared_ptr<const DetectorPeakResponse>;
  //orig_workingspace->user_peaks = std::vector<std::shared_ptr<const PeakDef>>;
  //orig_workingspace->automated_search_peaks = std::vector<std::shared_ptr<const PeakDef>>;
  //orig_workingspace->displayed_measurement  =  std::shared_ptr<const SpecUtils::Measurement>;
  //orig_workingspace->foreground = std::shared_ptr<SpecMeas>;
  //orig_workingspace->foreground_samplenums = std::set<int>;
  orig_workingspace->isHPGe = true;
  //orig_workingspace->error_msg = string;
  //orig_workingspace->matches = std::vector< std::vector<IsotopeMatch> > ;
  orig_workingspace->sortColumn = IsotopeSearchByEnergyModel::Column::ParentIsotope;
  orig_workingspace->sortOrder = Wt::SortOrder::AscendingOrder;
  //orig_workingspace->undoSentry = std::shared_ptr<void>
  //orig_workingspace->searchdoneCallback = searchdoneCallback

  //The 244 keV of Eu152 has rel. BR of 0.27, so if we set min BR to just below that, we should get one result
  double min_rel_br = 0.26;
  const double minHalfLife = 0.0;
  Wt::WFlags<IsotopeSearchByEnergyModel::RadSource> radiation;
  radiation |= IsotopeSearchByEnergyModel::RadSource::NuclideGammaOrXray;

  const std::vector<const SandiaDecay::Element *> elements;
  const std::vector<const SandiaDecay::Nuclide *> nuclides;
  const std::vector<const ReactionGamma::Reaction *> reactions;

  boost::function< void(void) > updatefcn;

  shared_ptr<IsotopeSearchByEnergyModel::SearchWorkingSpace> workingspace
      = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>( *orig_workingspace );

  IsotopeSearchByEnergyModel::setSearchEnergies( workingspace, min_rel_br, minHalfLife, radiation,
                                elements, nuclides, reactions, "", updatefcn );

  BOOST_CHECK( workingspace->matches.size() == 1 );
  if( workingspace->matches.size() > 0 )
  {
    const std::vector<IsotopeSearchByEnergyModel::IsotopeMatch> &match = workingspace->matches[0];
    BOOST_CHECK_EQUAL( match.size(), 3 );
    BOOST_REQUIRE( match.size() > 0 && match[0].m_nuclide );
    BOOST_CHECK( match[0].m_nuclide->symbol == "Eu152" );
  }

  // Check if we set Rel Br to just above the value for the 244 keV gamma, we shouldnt get any result
  min_rel_br = 0.28;
  workingspace = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>( *orig_workingspace );
  IsotopeSearchByEnergyModel::setSearchEnergies( workingspace, min_rel_br, minHalfLife, radiation,
                                elements, nuclides, reactions, "", updatefcn );

  BOOST_CHECK( workingspace->matches.empty() );


  // We'll check a little more complicated decay: the 1001 keV line of U238
  orig_workingspace->energies = { 1000.99 };
  orig_workingspace->windows = { 1.0 };

  min_rel_br = 0.16; //actual value is 0.17
  workingspace = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>( *orig_workingspace );

  const auto has_nuc = []( const vector<vector<IsotopeSearchByEnergyModel::IsotopeMatch>> &matches, const string &nuc) -> bool {
    for( const auto &match : matches )
    {
      for( const IsotopeSearchByEnergyModel::IsotopeMatch &iso : match )
      {
        if( iso.m_nuclide && iso.m_nuclide->symbol == nuc )
          return true;
      }
    }
    return false;
  };


  IsotopeSearchByEnergyModel::setSearchEnergies( workingspace, min_rel_br, minHalfLife, radiation,
                                elements, nuclides, reactions, "", updatefcn );
  BOOST_CHECK( has_nuc( workingspace->matches, "U238") );

  min_rel_br = 0.18;
  workingspace = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>( *orig_workingspace );
  IsotopeSearchByEnergyModel::setSearchEnergies( workingspace, min_rel_br, minHalfLife, radiation,
                                elements, nuclides, reactions, "", updatefcn );
  BOOST_CHECK( !has_nuc( workingspace->matches, "U238") );
}//BOOST_AUTO_TEST_CASE( testNuclideSearch )
