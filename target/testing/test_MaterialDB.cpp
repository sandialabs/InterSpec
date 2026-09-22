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
#include <vector>
#include <iostream>

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#define BOOST_TEST_MODULE MaterialDB_suite
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;
using namespace boost::unit_test;

namespace
{
const double g_cm3 = PhysicalUnits::g / PhysicalUnits::cm3;

// We need to set the static data directory, so the code knows where sandia.decay.xml and
//  MaterialDataBase.txt are located.
void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  const int argc = framework::master_test_suite().argc;
  char ** const argv = framework::master_test_suite().argv;

  string datadir;
  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
  }

  SpecUtils::ireplace_all( datadir, "%20", " " );

  if( datadir.empty() )
  {
    for( const auto &d : { "data", "../data", "../../data", "../../../data" } )
    {
      if( SpecUtils::is_file( SpecUtils::append_path(d, "sandia.decay.xml") ) )
      {
        datadir = d;
        break;
      }
    }//for( loop over candidate dirs )
  }//if( datadir.empty() )

  const string sandia_decay_file = SpecUtils::append_path( datadir, "sandia.decay.xml" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( sandia_decay_file ), "sandia.decay.xml not at '" << sandia_decay_file << "'" );

  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE_MESSAGE( db, "Error initing SandiaDecayDataBase" );
  BOOST_REQUIRE_MESSAGE( db->nuclide("U238"), "SandiaDecayDataBase empty?" );

  BOOST_REQUIRE_NO_THROW( MaterialDB::instance() );
  BOOST_REQUIRE( MaterialDB::instance() );
}//void set_data_dir()


/** Serializes `mat` with Material::toXml, prints the document to text, re-parses it, and returns
 the material Material::fromXml gives back - i.e., a full round-trip through what would be on disk.
 */
shared_ptr<const Material> xml_round_trip( const Material &mat )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();

  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> * const root = doc.allocate_node( rapidxml::node_element, "Test" );
  doc.append_node( root );

  rapidxml::xml_node<char> * const node = mat.toXml( root );
  BOOST_REQUIRE( node );
  BOOST_CHECK_EQUAL( SpecUtils::xml_name_str(node), "MaterialDefinition" );
  BOOST_REQUIRE( XML_FIRST_NODE(root, "MaterialDefinition") == node );

  string xml;
  rapidxml::print( back_inserter(xml), doc, 0 );
  BOOST_REQUIRE( !xml.empty() );

  vector<char> buffer( begin(xml), end(xml) );
  buffer.push_back( '\0' );

  rapidxml::xml_document<char> parsed;
  BOOST_REQUIRE_NO_THROW( parsed.parse<rapidxml::parse_default>( &buffer[0] ) );

  const rapidxml::xml_node<char> * const parsed_root = parsed.first_node();
  BOOST_REQUIRE( parsed_root );
  const rapidxml::xml_node<char> * const parsed_node = XML_FIRST_NODE( parsed_root, "MaterialDefinition" );
  BOOST_REQUIRE( parsed_node );

  shared_ptr<const Material> answer;
  BOOST_REQUIRE_NO_THROW( answer = Material::fromXml( parsed_node, db ) );
  BOOST_REQUIRE( answer );

  return answer;
}//xml_round_trip(...)

}//namespace


BOOST_AUTO_TEST_CASE( MaterialXmlRoundTrip )
{
  set_data_dir();

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( matdb && db );

  // An element-only material
  const shared_ptr<const Material> iron = matdb->material( "Fe" );
  BOOST_REQUIRE( iron );
  BOOST_REQUIRE( !iron->elements.empty() );
  {
    const shared_ptr<const Material> from_xml = xml_round_trip( *iron );
    BOOST_CHECK_NO_THROW( Material::equalEnough( *iron, *from_xml ) );
    BOOST_CHECK( Material::sameComposition( *iron, *from_xml ) );
    BOOST_CHECK_EQUAL( from_xml->name, iron->name );
    BOOST_CHECK_EQUAL( from_xml->description, iron->description );
    BOOST_CHECK_EQUAL( static_cast<int>(from_xml->source), static_cast<int>(iron->source) );
    BOOST_CHECK_CLOSE( from_xml->density / g_cm3, iron->density / g_cm3, 1.0E-4 );
  }

  // A material with nuclide components
  const shared_ptr<const Material> spent_fuel = matdb->material( "Spent Fuel" );
  BOOST_REQUIRE( spent_fuel );
  BOOST_REQUIRE( !spent_fuel->nuclides.empty() );
  BOOST_REQUIRE( !spent_fuel->elements.empty() );
  {
    const shared_ptr<const Material> from_xml = xml_round_trip( *spent_fuel );
    BOOST_CHECK_NO_THROW( Material::equalEnough( *spent_fuel, *from_xml ) );
    BOOST_CHECK_EQUAL( from_xml->nuclides.size(), spent_fuel->nuclides.size() );
    BOOST_CHECK_EQUAL( from_xml->elements.size(), spent_fuel->elements.size() );
  }

  // A user-entered chemical formula, with an explicit density
  const shared_ptr<const Material> formula = MaterialDB::materialFromChemicalFormula( "C0.5H0.2Ni0.6 d=2.2", db );
  BOOST_REQUIRE( formula );
  BOOST_CHECK_CLOSE( formula->density / g_cm3, 2.2, 1.0E-4 );
  {
    const shared_ptr<const Material> from_xml = xml_round_trip( *formula );
    BOOST_CHECK_NO_THROW( Material::equalEnough( *formula, *from_xml ) );
    BOOST_CHECK_EQUAL( static_cast<int>(from_xml->source), static_cast<int>(Material::kUser) );
  }

  // A density-modified copy of a database material - the case the whole mechanism exists for
  auto modified = make_shared<Material>( *iron );
  modified->density = static_cast<float>( 7.5 * g_cm3 );
  BOOST_CHECK( Material::sameComposition( *iron, *modified ) );
  BOOST_CHECK( !(*iron == *modified) );
  BOOST_CHECK( *iron != *modified );
  BOOST_CHECK_THROW( Material::equalEnough( *iron, *modified ), std::exception );
  {
    const shared_ptr<const Material> from_xml = xml_round_trip( *modified );
    BOOST_CHECK_NO_THROW( Material::equalEnough( *modified, *from_xml ) );
    BOOST_CHECK_THROW( Material::equalEnough( *iron, *from_xml ), std::exception );
    BOOST_CHECK_CLOSE( from_xml->density / g_cm3, 7.5, 1.0E-4 );
  }
}//BOOST_AUTO_TEST_CASE( MaterialXmlRoundTrip )


BOOST_AUTO_TEST_CASE( MaterialFromXmlWithoutDatabaseEntry )
{
  set_data_dir();

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( matdb && db );

  // A definition whose name is not in the database must still load - this is what keeps old
  //  files working if a material is removed from MaterialDataBase.txt
  const shared_ptr<const Material> lead = matdb->material( "Pb" );
  BOOST_REQUIRE( lead );

  auto renamed = make_shared<Material>( *lead );
  renamed->name = "Unobtainium";
  renamed->description = "Not in any database";
  BOOST_CHECK_THROW( matdb->material( renamed->name ), std::exception );

  const shared_ptr<const Material> from_xml = xml_round_trip( *renamed );
  BOOST_CHECK_NO_THROW( Material::equalEnough( *renamed, *from_xml ) );
  BOOST_CHECK_EQUAL( from_xml->name, "Unobtainium" );

  // Invalid definitions must throw, not produce garbage
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> * const root = doc.allocate_node( rapidxml::node_element, "Test" );
  doc.append_node( root );
  rapidxml::xml_node<char> * const node = lead->toXml( root );
  BOOST_REQUIRE( node );

  rapidxml::xml_node<char> * const element_node = XML_FIRST_NODE( node, "Element" );
  BOOST_REQUIRE( element_node );
  rapidxml::xml_attribute<char> * const symbol = XML_FIRST_ATTRIB( element_node, "Symbol" );
  BOOST_REQUIRE( symbol );
  symbol->value( "Xx", 2 );
  BOOST_CHECK_THROW( Material::fromXml( node, db ), std::exception );

  BOOST_CHECK_THROW( Material::fromXml( root, db ), std::exception );      //wrong node name
  BOOST_CHECK_THROW( Material::fromXml( nullptr, db ), std::exception );
}//BOOST_AUTO_TEST_CASE( MaterialFromXmlWithoutDatabaseEntry )


BOOST_AUTO_TEST_CASE( MaterialLookupHelpers )
{
  set_data_dir();

  const shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  BOOST_REQUIRE( matdb && db );

  const shared_ptr<const Material> iron = matdb->material( "Fe" );
  BOOST_REQUIRE( iron );

  // materialFromNameOrFormula: database names (in any of the accepted spellings), then formulas,
  //  and null (never a throw) for anything else
  BOOST_CHECK( MaterialDB::materialFromNameOrFormula( "Fe", db ) == iron );
  BOOST_CHECK( MaterialDB::materialFromNameOrFormula( iron->name, db ) == iron );
  BOOST_CHECK( MaterialDB::materialFromNameOrFormula( "iron", db ) == iron );

  const shared_ptr<const Material> formula = MaterialDB::materialFromNameOrFormula( "C0.5H0.2Ni0.6", db );
  BOOST_REQUIRE( formula );
  BOOST_CHECK_CLOSE( formula->density / g_cm3, 1.3, 1.0E-3 );

  BOOST_CHECK( !MaterialDB::materialFromNameOrFormula( "definitely not a material", db ) );
  BOOST_CHECK( !MaterialDB::materialFromNameOrFormula( "", db ) );

  // materialFromDefinitionOrName: the databases own instance when the definition matches it,
  //  the definition when it differs, and the name when there is no definition
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> * const root = doc.allocate_node( rapidxml::node_element, "Test" );
  doc.append_node( root );

  rapidxml::xml_node<char> * const iron_node = iron->toXml( root );
  BOOST_CHECK( MaterialDB::materialFromDefinitionOrName( iron_node, iron->name, db ) == iron );
  BOOST_CHECK( MaterialDB::materialFromDefinitionOrName( iron_node, "", db ) == iron );

  auto modified = make_shared<Material>( *iron );
  modified->density = static_cast<float>( 7.5 * g_cm3 );
  rapidxml::xml_node<char> * const modified_node = modified->toXml( root );
  const shared_ptr<const Material> resolved = MaterialDB::materialFromDefinitionOrName( modified_node, iron->name, db );
  BOOST_REQUIRE( resolved );
  BOOST_CHECK( resolved != iron );
  BOOST_CHECK_NO_THROW( Material::equalEnough( *modified, *resolved ) );

  auto renamed = make_shared<Material>( *iron );
  renamed->name = "Unobtainium";
  rapidxml::xml_node<char> * const renamed_node = renamed->toXml( root );
  const shared_ptr<const Material> from_def = MaterialDB::materialFromDefinitionOrName( renamed_node, "Unobtainium", db );
  BOOST_REQUIRE( from_def );
  BOOST_CHECK_EQUAL( from_def->name, "Unobtainium" );

  BOOST_CHECK( MaterialDB::materialFromDefinitionOrName( nullptr, "Fe", db ) == iron );
  BOOST_CHECK( !MaterialDB::materialFromDefinitionOrName( nullptr, "Unobtainium", db ) );

  // An unparsable definition falls back to the name
  rapidxml::xml_node<char> * const broken_node = iron->toXml( root );
  rapidxml::xml_node<char> * const density_node = XML_FIRST_NODE( broken_node, "Density" );
  BOOST_REQUIRE( density_node );
  density_node->value( "not a number", 12 );
  BOOST_CHECK( MaterialDB::materialFromDefinitionOrName( broken_node, "Fe", db ) == iron );
  BOOST_CHECK( !MaterialDB::materialFromDefinitionOrName( broken_node, "Unobtainium", db ) );
}//BOOST_AUTO_TEST_CASE( MaterialLookupHelpers )
