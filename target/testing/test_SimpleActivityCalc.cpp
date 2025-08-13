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

#include <cmath>
#include <string>
#include <cstdio>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>


#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SimpleActivityCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GammaInteractionCalc_suite
//#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace boost::unit_test;
using namespace GammaInteractionCalc;


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



BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_Serialization )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo deserialization
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test basic XML serialization and deserialization
  {
    // Create a test state with various data
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 661.657; // Cs-137 peak
    original_state.nuclideName = "Cs137";
    original_state.nuclideAgeStr = "30 y";
    original_state.distanceStr = "1 m";
    original_state.geometryType = SimpleActivityGeometryType::Point;
    
    // Set up a valid shielding using optional
    original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
    original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::Spherical;
    original_state.shielding->m_isGenericMaterial = false;
    original_state.shielding->m_forFitting = false;
    original_state.shielding->m_dimensions[0] = 2.0; // radius in cm 
    original_state.shielding->m_dimensions[1] = 0.0;
    original_state.shielding->m_dimensions[2] = 0.0;
    original_state.shielding->m_fitDimensions[0] = false;
    original_state.shielding->m_fitDimensions[1] = false;
    original_state.shielding->m_fitDimensions[2] = false;
    
    // Create XML document and serialize
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
    doc.append_node( root );
    
    BOOST_REQUIRE_NO_THROW( original_state.serialize( root ) );
    
    // Check that the XML contains expected nodes
    rapidxml::xml_node<char> *state_node = XML_FIRST_NODE( root, "SimpleActivityCalcState" );
    BOOST_REQUIRE_MESSAGE( state_node, "SimpleActivityCalcState node not created" );
    
    // Check version attribute
    rapidxml::xml_attribute<char> *version_attr = state_node->first_attribute( "version", 7 );
    BOOST_REQUIRE_MESSAGE( version_attr, "Version attribute not found" );
    BOOST_CHECK_EQUAL( std::string(version_attr->value()), "0" );
    
    // Check individual field nodes exist
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "PeakEnergy" ), "PeakEnergy node missing" );
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "NuclideName" ), "NuclideName node missing" );
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "NuclideAgeStr" ), "NuclideAgeStr node missing" );
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "DistanceStr" ), "DistanceStr node missing" );
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "GeometryType" ), "GeometryType node missing" );
    BOOST_CHECK_MESSAGE( XML_FIRST_NODE( state_node, "Shielding" ), "Shielding node missing" );
    
    // Test deserialization
    SimpleActivityCalcState deserialized_state;
    BOOST_REQUIRE_NO_THROW( deserialized_state.deSerialize( state_node, &matdb ) );
    
    // Verify all fields were restored correctly
    BOOST_CHECK_CLOSE( deserialized_state.peakEnergy, original_state.peakEnergy, 1e-6 );
    BOOST_CHECK_EQUAL( deserialized_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( deserialized_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( deserialized_state.distanceStr, original_state.distanceStr );
    BOOST_CHECK( deserialized_state.geometryType == original_state.geometryType );
    
    // Verify shielding data  
    BOOST_CHECK( deserialized_state.shielding.has_value() );
    BOOST_CHECK( original_state.shielding.has_value() );
    BOOST_REQUIRE_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( *deserialized_state.shielding, *original_state.shielding ) );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_Serialization )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_GeometryTypes )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo deserialization
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test serialization/deserialization of all geometry types
  const std::vector<SimpleActivityGeometryType> geometry_types = {
    SimpleActivityGeometryType::Point,
    SimpleActivityGeometryType::Plane,
    SimpleActivityGeometryType::SelfAttenuating,
    SimpleActivityGeometryType::TraceSrc
  };
  
  for( const auto geom_type : geometry_types )
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1332.5; // Co-60 peak
    original_state.nuclideName = "Co60";
    original_state.nuclideAgeStr = "5.27 y";
    original_state.distanceStr = "50 cm";
    original_state.geometryType = geom_type;
    
    // Set up shielding using optional
    original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
    original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
    original_state.shielding->m_isGenericMaterial = true;
    original_state.shielding->m_fitDimensions[0] = false;
    original_state.shielding->m_fitDimensions[1] = false;
    original_state.shielding->m_fitDimensions[2] = false;
    original_state.shielding->m_forFitting = false;
    original_state.shielding->m_dimensions[0] = 82; //atomic number
    original_state.shielding->m_dimensions[1] = 10 * PhysicalUnits::g_per_cm2; //areal density
    
    // Serialize
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
    doc.append_node( root );
    
    BOOST_REQUIRE_NO_THROW( original_state.serialize( root ) );
    
    // Deserialize
    rapidxml::xml_node<char> *state_node = XML_FIRST_NODE( root, "SimpleActivityCalcState" );
    SimpleActivityCalcState deserialized_state;
    BOOST_REQUIRE_NO_THROW( deserialized_state.deSerialize( state_node, &matdb ) );
    
    // Verify geometry type is preserved
    BOOST_CHECK_MESSAGE( deserialized_state.geometryType == geom_type, 
                        "Geometry type mismatch for type " << static_cast<int>(geom_type) );
    BOOST_CHECK_MESSAGE( deserialized_state == original_state,
                        "SimpleActivityCalcState didnt make it around" );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_GeometryTypes )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_ErrorHandling )
{
  set_data_dir();
  
  // Test error handling for invalid XML
  {
    // Test with nullptr node
    SimpleActivityCalcState state;
    BOOST_CHECK_THROW( state.deSerialize( nullptr, nullptr ), std::exception );
  }
  
  // Test with wrong node name
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *wrong_node = doc.allocate_node( rapidxml::node_element, "WrongNodeName" );
    doc.append_node( wrong_node );
    
    SimpleActivityCalcState state;
    BOOST_CHECK_THROW( state.deSerialize( wrong_node, nullptr ), std::exception );
  }
  
  // Test with invalid geometry type value
  {
    rapidxml::xml_document<char> doc;
    rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "SimpleActivityCalcState" );
    doc.append_node( root );
    
    // Add version attribute
    rapidxml::xml_attribute<char> *version_attr = doc.allocate_attribute( "version", "0" );
    root->append_attribute( version_attr );
    
    // Add required fields
    rapidxml::xml_node<char> *energy_node = doc.allocate_node( rapidxml::node_element, "PeakEnergy", "661.657" );
    root->append_node( energy_node );
    
    // Add invalid geometry type
    rapidxml::xml_node<char> *geom_node = doc.allocate_node( rapidxml::node_element, "GeometryType", "999" );
    root->append_node( geom_node );
    
    SimpleActivityCalcState state;
    BOOST_CHECK_THROW( state.deSerialize( root, nullptr ), std::exception );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_ErrorHandling )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_RoundTrip )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo deserialization
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test round-trip serialization with complex data
  SimpleActivityCalcState original_state;
  original_state.peakEnergy = 1173.228; // Co-60 first peak
  original_state.nuclideName = "Co60";
  original_state.nuclideAgeStr = "2.5 years";
  original_state.distanceStr = "75 cm";
  original_state.geometryType = SimpleActivityGeometryType::SelfAttenuating;
  
  // Set up complex physical shielding using optional
  original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
  original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::Spherical;
  original_state.shielding->m_isGenericMaterial = false;
  original_state.shielding->m_forFitting = true;
  original_state.shielding->m_dimensions[0] = 1.5; // radius in cm
  original_state.shielding->m_dimensions[1] = 0.0; // length in cm
  original_state.shielding->m_dimensions[2] = 0.0; // unused
  original_state.shielding->m_fitDimensions[0] = true;
  original_state.shielding->m_fitDimensions[1] = false;
  original_state.shielding->m_fitDimensions[2] = false;
  
  // First serialization
  rapidxml::xml_document<char> doc1;
  rapidxml::xml_node<char> *root1 = doc1.allocate_node( rapidxml::node_element, "root" );
  doc1.append_node( root1 );
  
  BOOST_REQUIRE_NO_THROW( original_state.serialize( root1 ) );
  
  
  // First deserialization
  rapidxml::xml_node<char> *state_node1 = XML_FIRST_NODE( root1, "SimpleActivityCalcState" );
  SimpleActivityCalcState intermediate_state;
  BOOST_REQUIRE_NO_THROW( intermediate_state.deSerialize( state_node1, &matdb ) );
  
  BOOST_CHECK( intermediate_state == original_state );
  
  // Second serialization
  rapidxml::xml_document<char> doc2;
  rapidxml::xml_node<char> *root2 = doc2.allocate_node( rapidxml::node_element, "root" );
  doc2.append_node( root2 );
  
  BOOST_REQUIRE_NO_THROW( intermediate_state.serialize( root2 ) );
  
  // Second deserialization
  rapidxml::xml_node<char> *state_node2 = XML_FIRST_NODE( root2, "SimpleActivityCalcState" );
  SimpleActivityCalcState final_state;
  BOOST_REQUIRE_NO_THROW( final_state.deSerialize( state_node2, &matdb ) );
  
  // Verify all data survived round-trip exactly
  BOOST_CHECK_CLOSE( final_state.peakEnergy, original_state.peakEnergy, 1e-10 );
  BOOST_CHECK_EQUAL( final_state.nuclideName, original_state.nuclideName );
  BOOST_CHECK_EQUAL( final_state.nuclideAgeStr, original_state.nuclideAgeStr );
  BOOST_CHECK_EQUAL( final_state.distanceStr, original_state.distanceStr );
  BOOST_CHECK( final_state.geometryType == original_state.geometryType );
  
  // Verify complex shielding data
  BOOST_CHECK( final_state.shielding.has_value() );
  BOOST_CHECK( original_state.shielding.has_value() );
  BOOST_REQUIRE_NO_THROW( ShieldingSourceFitCalc::ShieldingInfo::equalEnough( *final_state.shielding, *original_state.shielding ) );
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_RoundTrip )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_EmptyFields )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo deserialization
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test serialization/deserialization with empty string fields
  SimpleActivityCalcState original_state;
  original_state.peakEnergy = 511.0; // Annihilation peak
  original_state.nuclideName = ""; // Empty nuclide name
  original_state.nuclideAgeStr = ""; // Empty age
  original_state.distanceStr = ""; // Empty distance
  original_state.geometryType = SimpleActivityGeometryType::Plane;
  
  // No shielding for this test (leave as empty optional)
  // original_state.shielding is already default-constructed as empty optional
  
  // Serialize
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *root = doc.allocate_node( rapidxml::node_element, "root" );
  doc.append_node( root );
  
  BOOST_REQUIRE_NO_THROW( original_state.serialize( root ) );
  
  // Deserialize
  rapidxml::xml_node<char> *state_node = XML_FIRST_NODE( root, "SimpleActivityCalcState" );
  SimpleActivityCalcState deserialized_state;
  BOOST_REQUIRE_NO_THROW( deserialized_state.deSerialize( state_node, &matdb ) );
  
  // Verify empty fields are preserved
  BOOST_CHECK_CLOSE( deserialized_state.peakEnergy, original_state.peakEnergy, 1e-6 );
  BOOST_CHECK_EQUAL( deserialized_state.nuclideName, original_state.nuclideName );
  BOOST_CHECK_EQUAL( deserialized_state.nuclideAgeStr, original_state.nuclideAgeStr );
  BOOST_CHECK_EQUAL( deserialized_state.distanceStr, original_state.distanceStr );
  BOOST_CHECK( deserialized_state.geometryType == original_state.geometryType );
  
  // Verify shielding is empty in both states
  BOOST_CHECK( !deserialized_state.shielding.has_value() );
  BOOST_CHECK( !original_state.shielding.has_value() );
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_XML_EmptyFields )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_EqualEnough )
{
  set_data_dir();
  
  // Test that equalEnough throws exceptions with meaningful messages
  SimpleActivityCalcState state1, state2;
  
  // Test identical states (should not throw)
  BOOST_REQUIRE_NO_THROW( SimpleActivityCalcState::equalEnough( state1, state2 ) );
  
  // Test different peak energies
  state2.peakEnergy = 123.456;
  BOOST_CHECK_THROW( SimpleActivityCalcState::equalEnough( state1, state2 ), std::runtime_error );
  
  // Test different nuclide names
  state2.peakEnergy = state1.peakEnergy; // Reset
  state2.nuclideName = "Cs137";
  BOOST_CHECK_THROW( SimpleActivityCalcState::equalEnough( state1, state2 ), std::runtime_error );
  
  // Test different geometry types
  state2.nuclideName = state1.nuclideName; // Reset
  state2.geometryType = SimpleActivityGeometryType::Plane;
  BOOST_CHECK_THROW( SimpleActivityCalcState::equalEnough( state1, state2 ), std::runtime_error );
  
  // Test different shielding presence
  state2.geometryType = state1.geometryType; // Reset
  state2.shielding = ShieldingSourceFitCalc::ShieldingInfo();
  BOOST_CHECK_THROW( SimpleActivityCalcState::equalEnough( state1, state2 ), std::runtime_error );
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_EqualEnough )
