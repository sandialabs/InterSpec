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

#include <array>
#include <limits>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <sstream>
#include <functional>
#include <iostream>

#define BOOST_TEST_MODULE test_DrfAppearance
#include <boost/test/included/unit_test.hpp>

#include <Wt/Utils.h>

#include <nlohmann/json.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DrfAppearance.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;


namespace
{
  // A side view in the DRF-appearance frame: mm, face centre at (0,0), face on the left.
  const char * const ns_icon =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<!-- test icon -->\n"
    "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"-5 -60 210 120\">\n"
    "  <g fill=\"#000\"><path d=\"M0,-40 H200 V40 H0 Z\"/><path d=\"M20,-60 H60 V-40 H20 Z\" fill-opacity=\"0.5\"/></g>\n"
    "</svg>\n";

  string svg_with( const string &inner, const string &root_attribs = "viewBox=\"0 0 10 10\"" )
  {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" " + root_attribs + ">" + inner + "</svg>";
  }

  // One triangle: x in [0, 0.1], y in [-0.05, 0.05], z in [-0.2, 0] (metres); with `glb_bin`, the
  //  bytes a .glb would carry instead of the data URI.
  nlohmann::json triangle_gltf( string *glb_bin = nullptr )
  {
    const float pos[9] = { 0.0f, -0.05f, 0.0f,  0.1f, 0.0f, -0.2f,  0.0f, 0.05f, -0.1f };
    const uint16_t idx[4] = { 0, 1, 2, 0 };  //the last is padding
    string bin( reinterpret_cast<const char *>( pos ), sizeof(pos) );
    bin.append( reinterpret_cast<const char *>( idx ), sizeof(idx) );

    nlohmann::json doc;
    doc["asset"] = { {"version", "2.0"} };
    doc["scene"] = 0;
    doc["scenes"] = { { {"nodes", {0}} } };
    doc["nodes"] = { { {"mesh", 0}, {"name", "body"} } };
    doc["meshes"] = { { {"primitives", { { {"attributes", { {"POSITION", 0} }}, {"indices", 1}, {"material", 0} } }} } };
    doc["materials"] = { { {"pbrMetallicRoughness", { {"baseColorFactor", {0.5, 0.5, 0.5, 1.0}} }} } };
    doc["accessors"] = {
      { {"bufferView", 0}, {"componentType", 5126}, {"count", 3}, {"type", "VEC3"},
        {"min", {0.0, -0.05, -0.2}}, {"max", {0.1, 0.05, 0.0}} },
      { {"bufferView", 1}, {"componentType", 5123}, {"count", 3}, {"type", "SCALAR"} } };
    doc["bufferViews"] = { { {"buffer", 0}, {"byteOffset", 0}, {"byteLength", 36} },
                           { {"buffer", 0}, {"byteOffset", 36}, {"byteLength", 6} } };
    doc["buffers"] = { { {"byteLength", bin.size()} } };
    if( glb_bin )
      *glb_bin = bin;
    else
      doc["buffers"][0]["uri"] = "data:application/octet-stream;base64," + Wt::Utils::base64Encode( bin, false );
    return doc;
  }//triangle_gltf()


  void append_le_uint32( string &out, const uint32_t v )
  {
    for( int i = 0; i < 4; ++i )
      out.push_back( static_cast<char>( (v >> (8*i)) & 0xFF ) );
  }

  string to_glb( const nlohmann::json &doc, const string &bin )
  {
    string json = doc.dump();
    while( json.size() % 4 )
      json.push_back( ' ' );

    string out = "glTF";
    append_le_uint32( out, 2 );
    append_le_uint32( out, static_cast<uint32_t>( 12 + 8 + json.size() + 8 + bin.size() ) );
    append_le_uint32( out, static_cast<uint32_t>( json.size() ) );
    append_le_uint32( out, 0x4E4F534A );
    out += json;
    append_le_uint32( out, static_cast<uint32_t>( bin.size() ) );
    append_le_uint32( out, 0x004E4942 );
    out += bin;
    return out;
  }//to_glb(...)


  shared_ptr<DrfAppearance> test_appearance()
  {
    auto look = make_shared<DrfAppearance>();
    look->icon_svg = DrfAppearance::validated_svg( ns_icon );
    look->model_gltf = DrfAppearance::validated_gltf( triangle_gltf().dump() );
    return look;
  }


  shared_ptr<DetectorPeakResponse> test_drf()
  {
    auto drf = make_shared<DetectorPeakResponse>( "look test", "a DRF with a picture" );
    drf->setIntrinsicEfficiencyFormula( "0.3*exp(-0.001*x)", 5.0*PhysicalUnits::cm,
                                        PhysicalUnits::keV, 0.0f, 0.0f,
                                        DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );
    return drf;
  }


  /** The DRF as XML, and read back. */
  shared_ptr<DetectorPeakResponse> xml_round_trip( const DetectorPeakResponse &drf, string *xml_out = nullptr )
  {
    rapidxml::xml_document<char> doc;
    drf.toXml( &doc, &doc );
    string xml;
    rapidxml::print( std::back_inserter( xml ), doc, 0 );
    if( xml_out )
      *xml_out = xml;

    vector<char> buffer( begin(xml), end(xml) );
    buffer.push_back( '\0' );
    rapidxml::xml_document<char> in_doc;
    in_doc.parse<rapidxml::parse_trim_whitespace>( buffer.data() );

    auto answer = make_shared<DetectorPeakResponse>();
    answer->fromXml( in_doc.first_node( "DetectorPeakResponse" ) );
    return answer;
  }


  string test_file_dir()
  {
    const int argc = boost::unit_test::framework::master_test_suite().argc;
    char ** const argv = boost::unit_test::framework::master_test_suite().argv;
    string data_dir, file_dir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        data_dir = arg.substr( 10 );
      if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
        file_dir = arg.substr( 14 );
    }

    if( !data_dir.empty() )
    {
      DecayDataBaseServer::setDecayXmlFile( SpecUtils::append_path( data_dir, "sandia.decay.xml" ) );
      InterSpec::setStaticDataDirectory( data_dir );
    }
    return file_dir;
  }//test_file_dir()
}//namespace


BOOST_AUTO_TEST_CASE( svg_validation )
{
  const string svg = DrfAppearance::validated_svg( ns_icon );
  BOOST_CHECK( svg.find( "<?xml" ) == string::npos );       //declaration dropped
  BOOST_CHECK( svg.find( "test icon" ) == string::npos );   //and comments
  BOOST_CHECK( svg.find( "viewBox=\"-5 -60 210 120\"" ) != string::npos );
  BOOST_CHECK_EQUAL( DrfAppearance::validated_svg( svg ), svg );  //stable, so a re-read keeps it

  DrfAppearance look;
  look.icon_svg = svg;
  const array<double,4> vb = look.icon_view_box();
  BOOST_CHECK_EQUAL( vb[0], -5.0 );
  BOOST_CHECK_EQUAL( vb[1], -60.0 );
  BOOST_CHECK_EQUAL( vb[2], 210.0 );
  BOOST_CHECK_EQUAL( vb[3], 120.0 );

  // The data URI decodes back to the stored SVG
  const string uri = look.icon_data_uri();
  BOOST_REQUIRE( SpecUtils::starts_with( uri, "data:image/svg+xml;base64," ) );
  BOOST_CHECK_EQUAL( Wt::Utils::base64Decode( uri.substr( 26 ) ), svg );

  // Width follows the viewBox, capped
  const string html = look.icon_html( 16.0, 36.0 );
  BOOST_CHECK( html.find( "width:28.0px;height:16.0px;" ) != string::npos );
  BOOST_CHECK( look.icon_html( 16.0, 20.0 ).find( "width:20.0px;" ) != string::npos );

  // A missing namespace is added, as an image needs it
  BOOST_CHECK( DrfAppearance::validated_svg( "<svg viewBox=\"0 0 1 1\"/>" ).find( "http://www.w3.org/2000/svg" ) != string::npos );

  // Nothing that could run or fetch
  const vector<string> bad = {
    svg_with( "<script>alert(1)</script>" ),
    svg_with( "<foreignObject><div/></foreignObject>" ),
    svg_with( "<rect onclick=\"alert(1)\"/>" ),
    svg_with( "<rect/>", "viewBox=\"0 0 10 10\" onload=\"alert(1)\"" ),
    svg_with( "<image href=\"https://example.com/a.png\"/>" ),
    svg_with( "<use xlink:href=\"other.svg#a\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"/>" ),
    svg_with( "<a href=\"javascript:alert(1)\"><rect/></a>" ),
    svg_with( "<rect style=\"fill:url(https://example.com/x)\"/>" ),
    svg_with( "<style>@import url(x.css);</style>" ),
    svg_with( "<rect/>", "" ),                                    //no viewBox
    svg_with( "<rect/>", "viewBox=\"0 0 0 10\"" ),                //zero width
    svg_with( "<rect/>", "viewBox=\"0 0 10\"" ),                  //three numbers
    "<!DOCTYPE svg [<!ENTITY x \"y\">]>" + svg_with( "<rect/>" ),
    "<html/>",
    "not xml at all <",
    svg_with( "<rect/>" + string( 70*1024, ' ' ) )                 //too big
  };
  for( const string &s : bad )
    BOOST_CHECK_THROW( DrfAppearance::validated_svg( s ), std::exception );

  // ...but in-document references and embedded rasters are fine
  BOOST_CHECK_NO_THROW( DrfAppearance::validated_svg( svg_with(
      "<defs><linearGradient id=\"g\"/></defs><rect fill=\"url(#g)\"/>"
      "<image href=\"data:image/png;base64,AAAA\"/>" ) ) );
}//svg_validation


BOOST_AUTO_TEST_CASE( gltf_validation )
{
  const string gltf = DrfAppearance::validated_gltf( triangle_gltf().dump( 2 ) );
  BOOST_CHECK( gltf.find( '\n' ) == string::npos );  //stored compact

  DrfAppearance look;
  look.model_gltf = gltf;
  array<double,6> bounds;
  BOOST_REQUIRE( look.model_bounds_m( bounds ) );
  BOOST_CHECK_SMALL( bounds[0], 1.0E-9 );
  BOOST_CHECK_CLOSE( bounds[3], 0.1, 1.0E-6 );
  BOOST_CHECK_CLOSE( bounds[1], -0.05, 1.0E-6 );
  BOOST_CHECK_CLOSE( bounds[4], 0.05, 1.0E-6 );
  BOOST_CHECK_CLOSE( bounds[2], -0.2, 1.0E-6 );
  BOOST_CHECK_SMALL( bounds[5], 1.0E-9 );

  // Node transforms move the bounds: a 1 m translation along x
  {
    nlohmann::json moved = nlohmann::json::parse( gltf );
    moved["nodes"][0]["translation"] = { 1.0, 0.0, 0.0 };
    DrfAppearance moved_look;
    moved_look.model_gltf = DrfAppearance::validated_gltf( moved.dump() );
    BOOST_REQUIRE( moved_look.model_bounds_m( bounds ) );
    BOOST_CHECK_CLOSE( bounds[0], 1.0, 1.0E-6 );
    BOOST_CHECK_CLOSE( bounds[3], 1.1, 1.0E-6 );
  }

  // A .glb is converted to embedded JSON, with the same geometry
  {
    string bin;
    const nlohmann::json glb_doc = triangle_gltf( &bin );
    const string converted = DrfAppearance::validated_gltf( to_glb( glb_doc, bin ) );
    const nlohmann::json j = nlohmann::json::parse( converted );
    BOOST_REQUIRE( j["buffers"][0].contains( "uri" ) );
    BOOST_CHECK( SpecUtils::starts_with( j["buffers"][0]["uri"].get<string>(), "data:" ) );
    DrfAppearance glb_look;
    glb_look.model_gltf = converted;
    BOOST_REQUIRE( glb_look.model_bounds_m( bounds ) );
    BOOST_CHECK_CLOSE( bounds[3], 0.1, 1.0E-6 );
  }

  // What can not travel inside a DRF, or be drawn
  const auto modified = []( std::function<void(nlohmann::json &)> change ) -> string {
    nlohmann::json doc = triangle_gltf();
    change( doc );
    return doc.dump();
  };
  const vector<string> bad = {
    modified( []( nlohmann::json &d ){ d["buffers"][0]["uri"] = "model.bin"; } ),
    modified( []( nlohmann::json &d ){ d["asset"]["version"] = "1.0"; } ),
    modified( []( nlohmann::json &d ){ d["extensionsRequired"] = { "KHR_draco_mesh_compression" }; } ),
    modified( []( nlohmann::json &d ){ d["meshes"][0]["primitives"][0]["mode"] = 1; } ),       //lines only
    modified( []( nlohmann::json &d ){ d["accessors"][0]["count"] = 30; } ),                    //past its view
    modified( []( nlohmann::json &d ){ d["bufferViews"][1]["byteLength"] = 600; } ),            //past the buffer
    modified( []( nlohmann::json &d ){ d["bufferViews"][0]["byteStride"] = 24; } ),             //interleaved
    modified( []( nlohmann::json &d ){ d["accessors"][0]["componentType"] = 5123; } ),          //quantized
    modified( []( nlohmann::json &d ){ d["meshes"][0]["primitives"][0]["indices"] = 7; } ),     //no accessor
    "{ not json",
    "glTF\x01\x00\x00\x00 truncated"
  };
  for( const string &s : bad )
    BOOST_CHECK_THROW( DrfAppearance::validated_gltf( s ), std::exception );
}//gltf_validation


BOOST_AUTO_TEST_CASE( drf_xml_round_trip )
{
  const shared_ptr<DetectorPeakResponse> plain = test_drf();
  string plain_xml;
  xml_round_trip( *plain, &plain_xml );

  shared_ptr<DetectorPeakResponse> drf = test_drf();
  const uint64_t hash = drf->hashValue();
  drf->setAppearance( test_appearance() );

  // Cosmetic: no change to the hash or to the version written, so older readers still read it
  BOOST_CHECK_EQUAL( drf->hashValue(), hash );
  string xml;
  const shared_ptr<DetectorPeakResponse> back = xml_round_trip( *drf, &xml );
  const auto version_of = []( const string &s ){ return s.substr( s.find( "version=" ), 11 ); };
  BOOST_CHECK_EQUAL( version_of( xml ), version_of( plain_xml ) );
  BOOST_CHECK( xml.find( "<Appearance" ) != string::npos );

  BOOST_REQUIRE( back->appearance() );
  BOOST_CHECK( *back->appearance() == *drf->appearance() );
  BOOST_CHECK_EQUAL( back->hashValue(), hash );
  BOOST_CHECK( *back == *drf );
#if( PERFORM_DEVELOPER_CHECKS )
  BOOST_CHECK_NO_THROW( DetectorPeakResponse::equalEnough( *drf, *back ) );

  // ...and a dropped appearance is caught, since the hash can not catch it
  BOOST_CHECK_THROW( DetectorPeakResponse::equalEnough( *drf, *plain ), std::exception );
#endif
  BOOST_CHECK( !(*drf == *plain) );

  // An empty appearance is no appearance
  drf->setAppearance( make_shared<DrfAppearance>() );
  BOOST_CHECK( !drf->appearance() );
}//drf_xml_round_trip


BOOST_AUTO_TEST_CASE( drf_extra_round_trip )
{
  // The database keeps it in the extras column; a DRF with only a picture must still write one.
  shared_ptr<DetectorPeakResponse> drf = test_drf();
  BOOST_CHECK( drf->drfExtraToXmlString().empty() );
  drf->setAppearance( test_appearance() );
  const string extra = drf->drfExtraToXmlString();
  BOOST_REQUIRE( !extra.empty() );

  auto back = test_drf();
  back->setDrfExtraFromXmlString( extra );
  BOOST_REQUIRE( back->appearance() );
  BOOST_CHECK( *back->appearance() == *drf->appearance() );

  // Re-reading an extras column without one clears it
  back->setDrfExtraFromXmlString( "" );
  BOOST_CHECK( !back->appearance() );
}//drf_extra_round_trip


BOOST_AUTO_TEST_CASE( bad_appearance_keeps_drf )
{
  // An icon that fails the check (here, with a script) is dropped on read; the model and the DRF
  //  itself are kept.
  shared_ptr<DetectorPeakResponse> drf = test_drf();
  auto look = test_appearance();
  look->icon_svg = svg_with( "<script>alert(1)</script>" );  //as a hand-edited file might have
  drf->setAppearance( look );

  shared_ptr<DetectorPeakResponse> back;
  BOOST_REQUIRE_NO_THROW( back = xml_round_trip( *drf ) );
  BOOST_CHECK( back->isValid() );
  BOOST_REQUIRE( back->appearance() );
  BOOST_CHECK( back->appearance()->icon_svg.empty() );
  BOOST_CHECK( !back->appearance()->model_gltf.empty() );

  // A format this build does not know is skipped, not an error
  const char *future =
    "<DrfExtra><Appearance version=\"3\"><Icon format=\"png\">AAAA</Icon></Appearance></DrfExtra>";
  auto other = test_drf();
  BOOST_REQUIRE_NO_THROW( other->setDrfExtraFromXmlString( future ) );
  BOOST_CHECK( !other->appearance() );
  BOOST_CHECK( other->isValid() );
}//bad_appearance_keeps_drf


BOOST_AUTO_TEST_CASE( n42_round_trip )
{
  // Saved with a spectrum file, the picture comes back with the DRF
  const string dir = test_file_dir();
  const string file = SpecUtils::append_path( dir, "uranium_40pct_selfatten_HPGe_15cm.n42" );
  if( dir.empty() || !SpecUtils::is_file( file ) )
  {
    BOOST_TEST_MESSAGE( "Skipping N42 round trip: pass --testfiledir=..." );
    return;
  }

  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( file, SpecUtils::ParserType::Auto ) );
  shared_ptr<DetectorPeakResponse> drf = test_drf();
  drf->setAppearance( test_appearance() );
  meas->setDetector( drf );

  stringstream strm;
  BOOST_REQUIRE( meas->write_2012_N42( strm ) );
  string data = strm.str();

  auto back = make_shared<SpecMeas>();
  BOOST_REQUIRE( back->load_N42_from_data( &data[0], &data[0] + data.size() ) );
  BOOST_REQUIRE( back->detector() );
  BOOST_REQUIRE( back->detector()->appearance() );
  BOOST_CHECK( *back->detector()->appearance() == *drf->appearance() );
}//n42_round_trip


BOOST_AUTO_TEST_CASE( hostile_inputs )
{
  // A node tree that lists the same children over and over (glTF allows each node one parent) is
  //  refused; and walking one anyway - as a hand-made DRF file could hold - visits each node once,
  //  rather than 50^19 times.
  {
    nlohmann::json doc = triangle_gltf();
    nlohmann::json nodes = nlohmann::json::array();
    for( int i = 0; i < 20; ++i )
    {
      nlohmann::json n = { {"mesh", 0} };
      if( i < 19 )
        n["children"] = vector<int>( 50, i + 1 );
      nodes.push_back( n );
    }
    doc["nodes"] = nodes;
    BOOST_CHECK_THROW( DrfAppearance::validated_gltf( doc.dump() ), std::exception );

    DrfAppearance look;
    look.model_gltf = doc.dump();
    array<double,6> bounds;
    BOOST_CHECK( look.model_bounds_m( bounds ) );
  }

  // What InterSpec does not use is not kept: neither a deep nesting (which serializing would
  //  recurse into) nor text (which would reach the browser).
  {
    nlohmann::json doc = triangle_gltf();
    doc["nodes"][0]["name"] = "<script>alert(1)</script>";
    doc["nodes"][0]["extras"] = { {"note", "<img src=x onerror=alert(1)>"} };
    string text = doc.dump();
    text.insert( text.size() - 1, ",\"extensions\":" + string( 100000, '[' ) + string( 100000, ']' ) );

    const string clean = DrfAppearance::validated_gltf( text );
    BOOST_CHECK( clean.find( "script" ) == string::npos );
    BOOST_CHECK( clean.find( "onerror" ) == string::npos );
    BOOST_CHECK( clean.find( "extensions" ) == string::npos );
    BOOST_CHECK_EQUAL( DrfAppearance::validated_gltf( clean ), clean );  //stable, so a re-read keeps it
  }

  // Accessor ranges can not wrap around
  {
    nlohmann::json doc = triangle_gltf();
    doc["accessors"][0]["count"] = uint64_t(1) << 62;
    BOOST_CHECK_THROW( DrfAppearance::validated_gltf( doc.dump() ), std::exception );
    doc = triangle_gltf();
    doc["accessors"][0]["byteOffset"] = std::numeric_limits<uint64_t>::max() - 4;
    BOOST_CHECK_THROW( DrfAppearance::validated_gltf( doc.dump() ), std::exception );
  }

  // Deep SVG nesting is refused before it is parsed - including when hidden behind a quoted ">"
  {
    string inner, hidden;
    for( int i = 0; i < 100; ++i )
    {
      inner += "<g>";
      hidden += "<g a=\"/>\">";
    }
    for( int i = 0; i < 100; ++i )
    {
      inner += "</g>";
      hidden += "</g>";
    }
    BOOST_CHECK_THROW( DrfAppearance::validated_svg( svg_with( inner ) ), std::exception );
    BOOST_CHECK_THROW( DrfAppearance::validated_svg( svg_with( hidden ) ), std::exception );
    BOOST_CHECK_NO_THROW( DrfAppearance::validated_svg( svg_with( "<g><g><g><rect/></g></g></g>" ) ) );
  }

  // What drawing editors store - author, local file paths - does not travel with the icon
  {
    const string svg = DrfAppearance::validated_svg(
      "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\""
      " xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\""
      " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" viewBox=\"0 0 10 10\""
      " sodipodi:docname=\"jdoe-icon.svg\" inkscape:export-filename=\"/Users/jdoe/Desktop/icon.png\">"
      "<metadata><rdf:RDF><rdf:Description>Jane Doe</rdf:Description></rdf:RDF></metadata>"
      "<sodipodi:namedview inkscape:window-width=\"1920\"/>"
      "<rect width=\"10\" height=\"10\" inkscape:label=\"body\"/></svg>" );
    BOOST_CHECK( svg.find( "jdoe" ) == string::npos );
    BOOST_CHECK( svg.find( "Jane" ) == string::npos );
    BOOST_CHECK( svg.find( "metadata" ) == string::npos );
    BOOST_CHECK( svg.find( "inkscape" ) == string::npos );
    BOOST_CHECK( svg.find( "<rect" ) != string::npos );
  }

  // A model that fails the check when a DRF is read is dropped; the icon and the DRF are kept
  {
    shared_ptr<DetectorPeakResponse> drf = test_drf();
    auto look = test_appearance();
    look->model_gltf = "{\"asset\":{\"version\":\"2.0\"},\"buffers\":[{\"uri\":\"file:///etc/passwd\",\"byteLength\":1}]}";
    drf->setAppearance( look );

    shared_ptr<DetectorPeakResponse> back;
    BOOST_REQUIRE_NO_THROW( back = xml_round_trip( *drf ) );
    BOOST_REQUIRE( back->appearance() );
    BOOST_CHECK( back->appearance()->model_gltf.empty() );
    BOOST_CHECK( !back->appearance()->icon_svg.empty() );
  }
}//hostile_inputs
