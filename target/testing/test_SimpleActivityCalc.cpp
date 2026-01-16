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
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"

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


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_BasicEncoding )
{
  set_data_dir();
  
  // Test basic URL encoding without shielding
  {
    SimpleActivityCalcState state;
    state.peakEnergy = 661.657;
    state.nuclideName = "Cs137";
    state.nuclideAgeStr = "30 y";
    state.distanceStr = "1 m";
    state.geometryType = SimpleActivityGeometryType::Point;
    
    const std::string url = state.encodeToUrl();
    
    BOOST_CHECK( !url.empty() );
    BOOST_CHECK( url.find("V=1") != std::string::npos );
    BOOST_CHECK( url.find("E=661.6") != std::string::npos );
    BOOST_CHECK( url.find("N=Cs137") != std::string::npos );
    BOOST_CHECK( url.find("AGE=30 y") != std::string::npos );
    BOOST_CHECK( url.find("DIST=1 m") != std::string::npos );
    BOOST_CHECK( url.find("GEOM=Point") != std::string::npos );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_BasicEncoding )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_BasicDecoding )
{
  set_data_dir();
  
  // Test basic URL decoding without shielding
  {
    const std::string url = "V=1&E=661.657&N=Cs137&AGE=30 y&DIST=1 m&GEOM=Point";
    
    SimpleActivityCalcState state;
    BOOST_REQUIRE_NO_THROW( state.decodeFromUrl( url ) );
    
    BOOST_CHECK_CLOSE( state.peakEnergy, 661.657, 1e-6 );
    BOOST_CHECK_EQUAL( state.nuclideName, "Cs137" );
    BOOST_CHECK_EQUAL( state.nuclideAgeStr, "30 y" );
    BOOST_CHECK_EQUAL( state.distanceStr, "1 m" );
    BOOST_CHECK( state.geometryType == SimpleActivityGeometryType::Point );
    BOOST_CHECK( !state.shielding.has_value() );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_BasicDecoding )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_RoundTrip )
{
  set_data_dir();
  
  // Test round-trip URL encoding/decoding
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1173.228;
    original_state.nuclideName = "Co60";
    original_state.nuclideAgeStr = "2.5 years";
    original_state.distanceStr = "75 cm";
    original_state.geometryType = SimpleActivityGeometryType::Plane;
    
    const std::string url = original_state.encodeToUrl();
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    BOOST_CHECK_LT( fabs(decoded_state.peakEnergy - original_state.peakEnergy), 0.01 );
    BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( decoded_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( decoded_state.distanceStr, original_state.distanceStr );
    BOOST_CHECK( decoded_state.geometryType == original_state.geometryType );
    BOOST_CHECK( !decoded_state.shielding.has_value() );
    BOOST_CHECK( !original_state.shielding.has_value() );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_RoundTrip )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_AllGeometryTypes )
{
  set_data_dir();
  
  // Test URL encoding/decoding for all geometry types
  const std::vector<SimpleActivityGeometryType> geometry_types = {
    SimpleActivityGeometryType::Point,
    SimpleActivityGeometryType::Plane,
    SimpleActivityGeometryType::SelfAttenuating,
    SimpleActivityGeometryType::TraceSrc
  };
  
  for( const auto geom_type : geometry_types )
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1332.5;
    original_state.nuclideName = "Co60";
    original_state.nuclideAgeStr = "5.27 y";
    original_state.distanceStr = "50 cm";
    original_state.geometryType = geom_type;
    
    const std::string url = original_state.encodeToUrl();
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    BOOST_CHECK_MESSAGE( decoded_state.geometryType == geom_type, 
                        "Geometry type mismatch for type " << static_cast<int>(geom_type) );
    BOOST_CHECK_CLOSE( decoded_state.peakEnergy, original_state.peakEnergy, 1e-6 );
    BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_AllGeometryTypes )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_EnergyPrecision )
{
  set_data_dir();
  
  // Test that peak energy matches to nearest 0.01 as specified in requirements
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 661.65789; // More precision than 0.01
    original_state.nuclideName = "Cs137";
    original_state.geometryType = SimpleActivityGeometryType::Point;
    
    const std::string url = original_state.encodeToUrl();
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    // Energy should match within 0.01 keV as specified
    BOOST_CHECK( std::abs(decoded_state.peakEnergy - original_state.peakEnergy) <= 0.01 );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_EnergyPrecision )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_SpecialCharacters )
{
  set_data_dir();
  
  // Test URL encoding/decoding with special characters
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 511.0;
    original_state.nuclideName = "F18"; // Contains numbers
    original_state.nuclideAgeStr = "109.8 m"; // Contains period
    original_state.distanceStr = "1.5 m"; // Contains period
    original_state.geometryType = SimpleActivityGeometryType::Point;
    
    const std::string url = original_state.encodeToUrl();
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( decoded_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( decoded_state.distanceStr, original_state.distanceStr );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_SpecialCharacters )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_EmptyFields )
{
  set_data_dir();
  
  // Test URL encoding/decoding with empty fields
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1173.2;
    original_state.nuclideName = "";
    original_state.nuclideAgeStr = "";
    original_state.distanceStr = "";
    original_state.geometryType = SimpleActivityGeometryType::Plane;
    
    const std::string url = original_state.encodeToUrl();
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    BOOST_CHECK_CLOSE( decoded_state.peakEnergy, original_state.peakEnergy, 1e-6 );
    BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( decoded_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( decoded_state.distanceStr, original_state.distanceStr );
    BOOST_CHECK( decoded_state.geometryType == original_state.geometryType );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_EmptyFields )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_ErrorHandling )
{
  set_data_dir();
  
  // Test error handling for invalid URLs
  {
    SimpleActivityCalcState state;
    
    // Test with empty string
    BOOST_CHECK_THROW( state.decodeFromUrl( "" ), std::exception );
    
    // Test with invalid version
    BOOST_CHECK_THROW( state.decodeFromUrl( "V=2&E=661" ), std::exception );
    
    // Test with missing version
    BOOST_CHECK_THROW( state.decodeFromUrl( "E=661&N=Cs137" ), std::exception );
    
    // Test with invalid energy
    BOOST_CHECK_THROW( state.decodeFromUrl( "V=1&E=invalid" ), std::exception );
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_ErrorHandling )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_WithGenericShielding )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo handling
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test URL round-trip with generic shielding
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1173.228;
    original_state.nuclideName = "Co60";
    original_state.nuclideAgeStr = "5.27 y";
    original_state.distanceStr = "1 m";
    original_state.geometryType = SimpleActivityGeometryType::Point;
    
    // Set up generic shielding
    original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
    original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::Spherical;
    original_state.shielding->m_isGenericMaterial = true;
    original_state.shielding->m_forFitting = false;
    original_state.shielding->m_dimensions[0] = 82; // atomic number (Lead)
    original_state.shielding->m_dimensions[1] = 5.0 * PhysicalUnits::g_per_cm2; // areal density
    original_state.shielding->m_dimensions[2] = 0.0;
    original_state.shielding->m_fitDimensions[0] = false;
    original_state.shielding->m_fitDimensions[1] = false;
    original_state.shielding->m_fitDimensions[2] = false;
    
    const std::string url = original_state.encodeToUrl();
    BOOST_CHECK( !url.empty() );
    BOOST_CHECK( url.find("V=1") != std::string::npos );
    
    SimpleActivityCalcState decoded_state;
    BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url ) );
    
    // Verify basic fields
    BOOST_CHECK_LT( fabs(decoded_state.peakEnergy - original_state.peakEnergy), 0.01 );
    BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( decoded_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( decoded_state.distanceStr, original_state.distanceStr );
    BOOST_CHECK( decoded_state.geometryType == original_state.geometryType );
    
    // Verify shielding was preserved
    BOOST_CHECK( decoded_state.shielding.has_value() );
    BOOST_CHECK( original_state.shielding.has_value() );
    
    if( decoded_state.shielding.has_value() && original_state.shielding.has_value() )
    {
      BOOST_CHECK( decoded_state.shielding->m_isGenericMaterial );
      BOOST_CHECK( decoded_state.shielding->m_geometry == original_state.shielding->m_geometry );
      BOOST_CHECK_CLOSE( decoded_state.shielding->m_dimensions[0], original_state.shielding->m_dimensions[0], 1e-6 );
      BOOST_CHECK_CLOSE( decoded_state.shielding->m_dimensions[1], original_state.shielding->m_dimensions[1], 1e-6 );
    }
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_WithGenericShielding )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_WithMaterialShielding )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo handling
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test URL round-trip with material-based shielding
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 661.657;
    original_state.nuclideName = "Cs137";
    original_state.nuclideAgeStr = "30 y";
    original_state.distanceStr = "50 cm";
    original_state.geometryType = SimpleActivityGeometryType::Point;
    
    // Set up material-based shielding (Iron)
    const Material *iron = matdb.material( "Iron" );
    if( iron ) // Only run this test if Iron material is available
    {
      original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
      original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::Spherical;
      original_state.shielding->m_isGenericMaterial = false;
      original_state.shielding->m_forFitting = false;
      original_state.shielding->m_material = make_shared<Material>( *iron );
      original_state.shielding->m_dimensions[0] = 2.5; // thickness in cm
      original_state.shielding->m_dimensions[1] = 0.0;
      original_state.shielding->m_dimensions[2] = 0.0;
      original_state.shielding->m_fitDimensions[0] = false;
      original_state.shielding->m_fitDimensions[1] = false;
      original_state.shielding->m_fitDimensions[2] = false;
      
      const std::string url = original_state.encodeToUrl();
      BOOST_CHECK( !url.empty() );
      BOOST_CHECK( url.find("V=1") != std::string::npos );
      
      SimpleActivityCalcState decoded_state;
      BOOST_REQUIRE_NO_THROW( decoded_state.decodeFromUrl( url, &matdb ) );
      
      // Verify basic fields
      BOOST_CHECK_LT( fabs(decoded_state.peakEnergy - original_state.peakEnergy), 0.01 );
      BOOST_CHECK_EQUAL( decoded_state.nuclideName, original_state.nuclideName );
      BOOST_CHECK_EQUAL( decoded_state.nuclideAgeStr, original_state.nuclideAgeStr );
      BOOST_CHECK_EQUAL( decoded_state.distanceStr, original_state.distanceStr );
      BOOST_CHECK( decoded_state.geometryType == original_state.geometryType );
      
      // Verify shielding was preserved
      BOOST_CHECK( decoded_state.shielding.has_value() );
      BOOST_CHECK( original_state.shielding.has_value() );
      
      if( decoded_state.shielding.has_value() && original_state.shielding.has_value() )
      {
        BOOST_CHECK( !decoded_state.shielding->m_isGenericMaterial );
        BOOST_CHECK( decoded_state.shielding->m_geometry == original_state.shielding->m_geometry );
        BOOST_CHECK_CLOSE( decoded_state.shielding->m_dimensions[0], original_state.shielding->m_dimensions[0], 1e-6 );
        
        // Verify material name matches
        if( decoded_state.shielding->m_material && original_state.shielding->m_material )
        {
          BOOST_CHECK_EQUAL( decoded_state.shielding->m_material->name, original_state.shielding->m_material->name );
        }
      }
    }
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_WithMaterialShielding )


BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_ShieldingRoundTripMultiple )
{
  set_data_dir();
  
  // Set up MaterialDB for ShieldingInfo handling
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  // Test multiple encode/decode cycles to ensure stability
  {
    SimpleActivityCalcState original_state;
    original_state.peakEnergy = 1332.5;
    original_state.nuclideName = "Co60";
    original_state.nuclideAgeStr = "5.27 y";
    original_state.distanceStr = "75 cm";
    original_state.geometryType = SimpleActivityGeometryType::Plane;
    
    // Set up shielding with parameters compatible with SimpleActivityCalc constraints
    original_state.shielding = ShieldingSourceFitCalc::ShieldingInfo();
    original_state.shielding->m_geometry = GammaInteractionCalc::GeometryType::Spherical; // SimpleActivityCalc always uses Spherical
    original_state.shielding->m_isGenericMaterial = true;
    original_state.shielding->m_forFitting = false; // SimpleActivityCalc sets this to false
    original_state.shielding->m_dimensions[0] = 26; // atomic number (Iron)
    original_state.shielding->m_dimensions[1] = 7.87 * PhysicalUnits::g_per_cm2; // areal density
    original_state.shielding->m_dimensions[2] = 0.0;
    original_state.shielding->m_fitDimensions[0] = false; // SimpleActivityCalc sets all to false
    original_state.shielding->m_fitDimensions[1] = false;
    original_state.shielding->m_fitDimensions[2] = false;
    
    // First round-trip
    const std::string url1 = original_state.encodeToUrl();
    SimpleActivityCalcState state1;
    BOOST_REQUIRE_NO_THROW( state1.decodeFromUrl( url1 ) );
    
    // Second round-trip
    const std::string url2 = state1.encodeToUrl();
    SimpleActivityCalcState state2;
    BOOST_REQUIRE_NO_THROW( state2.decodeFromUrl( url2 ) );
    
    // Third round-trip
    const std::string url3 = state2.encodeToUrl();
    SimpleActivityCalcState final_state;
    BOOST_REQUIRE_NO_THROW( final_state.decodeFromUrl( url3 ) );
    
    // Verify all data survived multiple round-trips
    BOOST_CHECK_LT( fabs(final_state.peakEnergy - original_state.peakEnergy), 0.01 );
    BOOST_CHECK_EQUAL( final_state.nuclideName, original_state.nuclideName );
    BOOST_CHECK_EQUAL( final_state.nuclideAgeStr, original_state.nuclideAgeStr );
    BOOST_CHECK_EQUAL( final_state.distanceStr, original_state.distanceStr );
    BOOST_CHECK( final_state.geometryType == original_state.geometryType );
    
    // Verify shielding data survived multiple round-trips
    BOOST_CHECK( final_state.shielding.has_value() );
    BOOST_CHECK( original_state.shielding.has_value() );
    
    if( final_state.shielding.has_value() && original_state.shielding.has_value() )
    {
      BOOST_CHECK( final_state.shielding->m_isGenericMaterial == original_state.shielding->m_isGenericMaterial );
      BOOST_CHECK( final_state.shielding->m_forFitting == original_state.shielding->m_forFitting );
      BOOST_CHECK( final_state.shielding->m_geometry == original_state.shielding->m_geometry );
      BOOST_CHECK_CLOSE( final_state.shielding->m_dimensions[0], original_state.shielding->m_dimensions[0], 1e-6 );
      BOOST_CHECK_CLOSE( final_state.shielding->m_dimensions[1], original_state.shielding->m_dimensions[1], 1e-6 );
      BOOST_CHECK( final_state.shielding->m_fitDimensions[0] == original_state.shielding->m_fitDimensions[0] );
      BOOST_CHECK( final_state.shielding->m_fitDimensions[1] == original_state.shielding->m_fitDimensions[1] );
    }
  }
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalcState_URL_ShieldingRoundTripMultiple )


struct SimpleActivityTestConfig
{
  std::string description;
  std::string spectrumFile;
  double targetPeakEnergy;
  SimpleActivityGeometryType geometryType;
  std::optional<ShieldingSourceFitCalc::ShieldingInfo> shielding;
  bool backgroundSubtract;
  double distance;
  double expectedActivity;
  double expectedActivityUncertainty;
  double tolerancePercent;
  
  SimpleActivityTestConfig(const std::string& desc, const std::string& file, double energy, 
                          SimpleActivityGeometryType geom, bool bgSub = false, double dist = 100.0)
    : description(desc), spectrumFile(file), targetPeakEnergy(energy), geometryType(geom),
      backgroundSubtract(bgSub), distance(dist), expectedActivity(-1.0), 
      expectedActivityUncertainty(-1.0), tolerancePercent(1.0) {}
};


BOOST_AUTO_TEST_CASE( SimpleActivityCalc_Ra226_609keV_Tests )
{
  set_data_dir();
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  MaterialDB matdb;
  
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  BOOST_REQUIRE_NO_THROW( matdb.parseGadrasMaterialFile( materialfile, db, false ) );
  
  auto setup_generic_shielding = []() -> ShieldingSourceFitCalc::ShieldingInfo {
    ShieldingSourceFitCalc::ShieldingInfo shielding;
    shielding.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
    shielding.m_isGenericMaterial = true;
    shielding.m_forFitting = false;
    shielding.m_dimensions[0] = 26;
    shielding.m_dimensions[1] = 2.0 * PhysicalUnits::g_per_cm2;
    shielding.m_dimensions[2] = 0.0;
    shielding.m_fitDimensions[0] = false;
    shielding.m_fitDimensions[1] = false;
    shielding.m_fitDimensions[2] = false;
    return shielding;
  };
  
  auto setup_material_shielding = [&matdb]() -> ShieldingSourceFitCalc::ShieldingInfo {
    ShieldingSourceFitCalc::ShieldingInfo shielding;
    shielding.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
    shielding.m_isGenericMaterial = false;
    shielding.m_forFitting = false;
    
    const Material *iron = matdb.material( "Fe (iron)" );
    BOOST_REQUIRE( iron );
    shielding.m_material = std::make_shared<Material>( *iron );
    
    shielding.m_dimensions[0] = 0.5 * PhysicalUnits::cm;
    shielding.m_dimensions[1] = 0.0;
    shielding.m_dimensions[2] = 0.0;
    shielding.m_fitDimensions[0] = false;
    shielding.m_fitDimensions[1] = false;
    shielding.m_fitDimensions[2] = false;
    return shielding;
  };
  
  // Its not a test from first principles, but we'll test computation against some "accepted" answers (e.g., checked
  //  using the Activity/Shielding Fit tool.
  std::vector<SimpleActivityTestConfig> test_configs = {
    SimpleActivityTestConfig("Point source, no shielding, no background subtraction", 
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Point,
                             false,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Point source, generic shielding, no background subtraction",
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Point,
                             false,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Point source, generic shielding, with background subtraction",
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Point,
                             true,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Point source, material shielding",
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Point,
                             true,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Trace source",
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::TraceSrc,
                             true,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Infinite plane",
                           "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Plane,
                             true,
                             100.0*PhysicalUnits::cm),
    SimpleActivityTestConfig("Point source, no shielding, no background subtraction - 100m",
                             "test_data/SimpleActivityCalc/Ra226 Shielded.n42", 609.3,
                             SimpleActivityGeometryType::Point,
                             false,
                             100.0*PhysicalUnits::m)
  };
  
  
  // Point source, no shielding, no background subtraction
  test_configs[0].expectedActivity = 4.51*PhysicalUnits::microCi;
  test_configs[0].expectedActivityUncertainty = 0.04584*PhysicalUnits::microCi;
  
  // Point source, generic shielding, no background subtraction
  test_configs[1].shielding = setup_generic_shielding();
  test_configs[1].expectedActivity = 5.24*PhysicalUnits::microCi;
  test_configs[1].expectedActivityUncertainty = 0.05328*PhysicalUnits::microCi;
  
  // Point source, generic shielding, with background subtraction
  test_configs[2].shielding = setup_generic_shielding();
  test_configs[2].expectedActivity = 4.85*PhysicalUnits::microCi;
  test_configs[2].expectedActivityUncertainty = 0.05525*PhysicalUnits::microCi;
  
  // Point source, material shielding, with background subtraction
  test_configs[3].shielding = setup_material_shielding();
  test_configs[3].expectedActivity = 5.61*PhysicalUnits::microCi;
  test_configs[3].expectedActivityUncertainty = 0.0639*PhysicalUnits::microCi;
  
  // Trace source, with background subtraction
  test_configs[4].shielding = setup_material_shielding();
  ShieldingSourceFitCalc::TraceSourceInfo trace;
  trace.m_type = GammaInteractionCalc::TraceActivityType::TotalActivity;
  trace.m_fitActivity = true;
  //trace.m_nuclide = nuc;
  trace.m_activity = 0.0; //units are according to #m_type
  trace.m_relaxationDistance = 0.0; //only applicable to #TraceActivityType::ExponentialDistribution
  test_configs[4].shielding->m_traceSources.push_back( trace );
  test_configs[4].expectedActivity = 5.15*PhysicalUnits::microCi;
  test_configs[4].expectedActivityUncertainty = 0.05875*PhysicalUnits::microCi;
  
  // Infinite plane, with background subtraction
  test_configs[5].expectedActivity = 18.93*PhysicalUnits::curie*1.0E-12 / PhysicalUnits::cm2;
  test_configs[5].expectedActivityUncertainty = 0.21584*PhysicalUnits::curie*1.0E-12 / PhysicalUnits::cm2;
  
  //"Point source, no shielding, no background subtraction - 100m"
  test_configs[6].expectedActivity = 124.75*PhysicalUnits::curie*1.0E-3;
  test_configs[6].expectedActivityUncertainty = 1.27*PhysicalUnits::curie*1.0E-3;
  
  for( size_t i = 0; i < test_configs.size(); ++i )
  {
    const auto &config = test_configs[i];
    
    //cout << "Running test " << (i+1) << ": " << config.description << endl;
    
    const string analysis_tests = SpecUtils::append_path(g_test_file_dir, ".." );
    const string spec_path = SpecUtils::append_path(analysis_tests, config.spectrumFile);
    
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(spec_path),
                          "Spectrum file not found: '" << spec_path << "'" );
    
    SpecMeas meas;
    BOOST_REQUIRE_MESSAGE( meas.load_N42_file(spec_path), 
                          "Failed to load spectrum from '" << spec_path << "'" );
    
    BOOST_REQUIRE_MESSAGE( meas.num_measurements() >= 1, 
                          "Expected at least one measurement in spectrum file" );
    
    std::shared_ptr<const SpecUtils::Measurement> foreground, background;
    
    for( size_t j = 0; j < meas.num_measurements(); ++j )
    {
      auto measurement = meas.measurements()[j];
      if( measurement->source_type() == SpecUtils::SourceType::Foreground )
        foreground = measurement;
      else if( measurement->source_type() == SpecUtils::SourceType::Background && config.backgroundSubtract )
        background = measurement;
    }
    
    BOOST_REQUIRE_MESSAGE( foreground, "No foreground measurement found in spectrum file" );
    
    if( config.backgroundSubtract )
    {
      BOOST_REQUIRE_MESSAGE( background, "Background subtraction requested but no background measurement found" );
    }
    
    auto foreground_peaks = meas.peaks( {foreground->sample_number()} );
    BOOST_REQUIRE_MESSAGE( foreground_peaks && !foreground_peaks->empty(), 
                          "No peaks found in foreground measurement" );
    
    std::shared_ptr<const PeakDef> target_peak, bg_peak;
    
    for( auto peak : *foreground_peaks )
    {
      if( std::abs(peak->mean() - config.targetPeakEnergy) < 1.0 )
      {
        target_peak = peak;
        break;
      }
    }
    
    BOOST_REQUIRE_MESSAGE( target_peak, 
                          "Could not find " << config.targetPeakEnergy << " keV peak in foreground spectrum" );
    
    if( config.backgroundSubtract && background )
    {
      auto background_peaks = meas.peaks( {background->sample_number()} );
      if( background_peaks && !background_peaks->empty() )
      {
        for( auto peak : *background_peaks )
        {
          if( std::abs(peak->mean() - config.targetPeakEnergy) < 2.0 )
          {
            bg_peak = peak;
            break;
          }
        }
      }
    }
    
    SimpleActivityCalcInput input;
    input.peak = target_peak;
    input.background_peak = bg_peak;
    input.detector = meas.detector();
    input.foreground = foreground;
    input.background = background;
    input.distance = config.distance;
    input.geometryType = config.geometryType;
    input.age = PeakDef::defaultDecayTime( target_peak->parentNuclide() );
    
    if( config.shielding.has_value() )
    {
      input.shielding = { config.shielding.value() };
      if( !input.shielding.front().m_traceSources.empty() )
        input.shielding.front().m_traceSources.front().m_nuclide = target_peak->parentNuclide();
    }
    
    SimpleActivityCalcResult result = SimpleActivityCalc::performCalculation( input );
    
    BOOST_CHECK_MESSAGE( result.successful, "Calculation failed for test '" << config.description 
                        << "': " << result.errorMessage );
    
    if( result.successful )
    {
      //cout << "  Activity: " << result.activity << " ± " << result.activityUncertainty << endl;
      
      BOOST_REQUIRE( config.expectedActivity > 0.0 );
      
      double percent_diff = std::abs(result.activity - config.expectedActivity) / config.expectedActivity * 100.0;
      BOOST_CHECK_MESSAGE( percent_diff <= config.tolerancePercent,
                          "Activity " << result.activity << " differs from expected "
                          << config.expectedActivity << " by " << percent_diff << "% (tolerance: "
                          << config.tolerancePercent << "%)" );
    }//if( result.successful )
  }//for( size_t i = 0; i < test_configs.size(); ++i )
}//BOOST_AUTO_TEST_CASE( SimpleActivityCalc_Ra226_609keV_Tests )

