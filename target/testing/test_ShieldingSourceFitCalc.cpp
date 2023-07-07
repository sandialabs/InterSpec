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

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ShieldingSourceFitCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"


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


BOOST_AUTO_TEST_CASE( ShieldingInfoUri )
{
  using namespace ShieldingSourceFitCalc;
  
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  const Material * const iron = matdb.material("Fe (iron)");
  BOOST_REQUIRE( iron );
  
  const Material * const uranium = matdb.material("U (uranium)");
  BOOST_REQUIRE( uranium );
  
  {// Begin test Trace source serialization
    TraceSourceInfo trace;
    trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
    trace.m_fitActivity = true;
    BOOST_REQUIRE( trace.m_nuclide = db->nuclide( "U238" ) );
    trace.m_activity = 2.0*PhysicalUnits::microCi / PhysicalUnits::cm2;
    trace.m_relaxationDistance = 0;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
    
    trace.m_type = GammaInteractionCalc::TraceActivityType::ExponentialDistribution;
    trace.m_relaxationDistance = 5*PhysicalUnits::cm;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
  }// End test Trace source serialization
  
  
  {// Begin test Trace source serialization, with empty src
    TraceSourceInfo trace;
    trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
    trace.m_fitActivity = false;
    trace.m_nuclide = nullptr;
    trace.m_activity = 0.0f;
    trace.m_relaxationDistance = 0.0f;
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( trace.serialize( &doc ) );
      
      TraceSourceInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node() ) );
      BOOST_REQUIRE_NO_THROW( TraceSourceInfo::equalEnough( trace, from_xml ) );
    }
  }// End test Trace source serialization, with empty src
  

  
  
  {// Begin test ShieldingInfo encode/decode from URI/XML
    // URI currently does not encode self-attenuating source, trace-source information, or "truth" values.
    ShieldingInfo info;
    
    info.m_geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
    info.m_isGenericMaterial = false;
    info.m_forFitting = true;
    info.m_material = make_shared<Material>( *iron );
    info.m_dimensions[0] = 1.2*PhysicalUnits::meter;
    info.m_dimensions[1] = 1.1*PhysicalUnits::meter;
    info.m_dimensions[2] = 1.0E-6*PhysicalUnits::meter;
    info.m_fitDimensions[0] = false;
    info.m_fitDimensions[1] = true;
    info.m_fitDimensions[2] = false;
    
    
    {// Begin test URI - doesnt include truth values or trace/self-atten sources, so test before we set those
      const string uri = info.encodeStateToUrl();
      
      ShieldingInfo from_uri;
      BOOST_REQUIRE_NO_THROW( from_uri.handleAppUrl( uri, &matdb ) );
      BOOST_CHECK_NO_THROW( ShieldingInfo::equalEnough( info, from_uri ) );
    }// End test URI
    
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    // Truth info only tested for XML serialization
    info.m_truthDimensions[0] = 1.11E-6f;
    info.m_truthDimensions[1] = boost::optional<float>();
    info.m_truthDimensions[2] = 9.99E-3f;
    
    info.m_truthDimensionsTolerances[0] = 1.11E-7f;
    info.m_truthDimensionsTolerances[1] = 100.0f;
    info.m_truthDimensionsTolerances[2] = boost::optional<float>();
#endif
    
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( info.serialize( &doc ) );
      
      ShieldingInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node(), &matdb ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
    
    
    // Now modify `info` to include trace and self-attenuating sources
    info.m_material = make_shared<Material>( *uranium );
    
    
    // Self-atten source stuff
    info.m_fitMassFrac = true;
    
    {// Begin add self-atten source
      map<const SandiaDecay::Nuclide *,double> &el_to_nuc_frac = info.m_nuclideFractions;
      
      const SandiaDecay::Nuclide *nuc = db->nuclide( "U-238" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 0.8;
      
      nuc = db->nuclide( "U-235" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 0.15;
      
      nuc = db->nuclide( "U-232" );
      BOOST_REQUIRE( nuc );
      el_to_nuc_frac[nuc] = 2.5E-8;
    }// End add self-atten source
    
    // Trace-source stuff
    {// Begin add trace source source
      
      TraceSourceInfo trace;
      trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
      trace.m_fitActivity = false;
      trace.m_nuclide = db->nuclide( "Cs137" );
      trace.m_activity = 2.0*PhysicalUnits::microCi / PhysicalUnits::cm2;
      trace.m_relaxationDistance = 0.0f;
      info.m_traceSources.push_back( trace );
      
      
      trace.m_type = GammaInteractionCalc::TraceActivityType::ActivityPerGram;
      trace.m_fitActivity = true;
      trace.m_nuclide = db->nuclide( "Ba133" );
      trace.m_activity = 5.0*PhysicalUnits::microCi / PhysicalUnits::gram;
      trace.m_relaxationDistance = 0.0f;
      info.m_traceSources.push_back( trace );
      
      trace.m_type = GammaInteractionCalc::TraceActivityType::ExponentialDistribution;
      trace.m_fitActivity = true;
      trace.m_nuclide = db->nuclide( "I131" );
      trace.m_activity = 1.0E5*PhysicalUnits::bq / PhysicalUnits::m2;
      trace.m_relaxationDistance = 1.23*PhysicalUnits::cm;
      info.m_traceSources.push_back( trace );
    }// End add trace source source
    
    {
      rapidxml::xml_document<char> doc;
      BOOST_REQUIRE_NO_THROW( info.serialize( &doc ) );
      
      ShieldingInfo from_xml;
      BOOST_REQUIRE_NO_THROW( from_xml.deSerialize( doc.first_node(), &matdb ) );
      BOOST_REQUIRE_NO_THROW( ShieldingInfo::equalEnough( info, from_xml ) );
    }
  }// End test ShieldingInfo encode/decode from URI
  
}//BOOST_AUTO_TEST_CASE( ShieldingInfoUri )


