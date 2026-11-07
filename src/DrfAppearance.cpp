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

#include <set>
#include <cmath>
#include <array>
#include <limits>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <utility>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <Wt/Utils.h>

#include <nlohmann/json.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/DrfAppearance.h"

using namespace std;


namespace
{
  const size_t ns_max_svg_bytes = 64*1024;
  const size_t ns_max_gltf_bytes = 2*1024*1024;

  // Deeper SVG nesting than any drawing needs; rapidxml parses recursively, so a crafted file with
  //  thousands of levels could otherwise overflow a (small, non-main thread) stack.
  const size_t ns_max_svg_depth = 64;

  /** Lower-cased local name (any namespace prefix removed) of an XML element or attribute. */
  string local_name( const char *name, const size_t len )
  {
    string answer( name, len );
    const size_t colon = answer.find( ':' );
    if( colon != string::npos )
      answer = answer.substr( colon + 1 );
    SpecUtils::to_lower_ascii( answer );
    return answer;
  }//local_name(...)


  /** Throws unless every `url(...)` in CSS (or an attribute like `fill="url(#g)"`) is a fragment
   within the document; also refuses `@import`. */
  void check_css_urls( string css )
  {
    SpecUtils::to_lower_ascii( css );
    if( css.find( "@import" ) != string::npos )
      throw runtime_error( "SVG icons can not import style sheets." );

    size_t pos = 0;
    while( (pos = css.find( "url(", pos )) != string::npos )
    {
      pos += 4;
      while( (pos < css.size()) && (isspace( static_cast<unsigned char>(css[pos]) )
                                    || (css[pos] == '"') || (css[pos] == '\'')) )
        ++pos;
      if( (pos >= css.size()) || (css[pos] != '#') )
        throw runtime_error( "SVG icons can not refer to anything outside themselves." );
    }
  }//check_css_urls(...)


  /** The deepest element nesting in an XML text, found without parsing it (quoted attribute values,
   comments, and <!...> / <?...?> sections are skipped). */
  size_t max_element_depth( const string &xml )
  {
    size_t depth = 0, max_depth = 0, pos = 0;
    while( (pos = xml.find( '<', pos )) != string::npos )
    {
      const char next = ((pos + 1) < xml.size()) ? xml[pos + 1] : '\0';
      if( (next == '!') || (next == '?') )
      {
        const bool comment = (xml.compare( pos, 4, "<!--" ) == 0);
        const size_t end = comment ? xml.find( "-->", pos ) : xml.find( '>', pos );
        if( end == string::npos )
          break;
        pos = end + 1;
        continue;
      }

      if( next == '/' )
      {
        depth = depth ? (depth - 1) : 0;
        ++pos;
        continue;
      }

      // An opening tag: find its end, skipping quoted values, to see whether it closes itself.
      char quote = '\0';
      size_t i = pos + 1;
      for( ; i < xml.size(); ++i )
      {
        if( quote )
          quote = (xml[i] == quote) ? '\0' : quote;
        else if( (xml[i] == '"') || (xml[i] == '\'') )
          quote = xml[i];
        else if( xml[i] == '>' )
          break;
      }
      if( i >= xml.size() )
        break;

      if( xml[i - 1] != '/' )
        max_depth = std::max( max_depth, ++depth );
      pos = i + 1;
    }//while( more tags )

    return max_depth;
  }//max_element_depth(...)


  /** The namespace prefix of an XML name ("" if none), lower-cased. */
  string name_prefix( const char *name, const size_t len )
  {
    const string full( name, len );
    const size_t colon = full.find( ':' );
    string answer = (colon == string::npos) ? string() : full.substr( 0, colon );
    SpecUtils::to_lower_ascii( answer );
    return answer;
  }


  /** Removes what drawing editors store but nothing draws - `<metadata>` (author, licence), and
   Inkscape/Sodipodi elements and attributes (window sizes, local file paths) - since the icon
   travels inside shared DRF and spectrum files. */
  void strip_editor_data( rapidxml::xml_node<char> *el )
  {
    for( rapidxml::xml_attribute<char> *attr = el->first_attribute(); attr; )
    {
      rapidxml::xml_attribute<char> *next = attr->next_attribute();
      const string prefix = name_prefix( attr->name(), attr->name_size() );
      const string name( attr->name(), attr->name_size() );
      const bool keep = prefix.empty() || (prefix == "xlink") || (prefix == "xml")
                        || (name == "xmlns:xlink") || (name == "xmlns:svg");
      if( !keep )
        el->remove_attribute( attr );
      attr = next;
    }

    for( rapidxml::xml_node<char> *child = el->first_node(); child; )
    {
      rapidxml::xml_node<char> *next = child->next_sibling();
      if( child->type() == rapidxml::node_element )
      {
        const string prefix = name_prefix( child->name(), child->name_size() );
        if( (local_name( child->name(), child->name_size() ) == "metadata")
            || (!prefix.empty() && (prefix != "svg")) )
          el->remove_node( child );
        else
          strip_editor_data( child );
      }
      child = next;
    }
  }//strip_editor_data(...)


  /** Recursively checks an SVG element, and its children, for anything that could run or fetch. */
  void check_svg_element( const rapidxml::xml_node<char> *el )
  {
    static const set<string> forbidden{ "script", "foreignobject", "iframe", "object", "embed",
                                        "handler", "listener", "audio", "video" };
    const string name = local_name( el->name(), el->name_size() );
    if( forbidden.count( name ) )
      throw runtime_error( "SVG icons can not contain '" + name + "' elements." );

    for( const rapidxml::xml_attribute<char> *attr = el->first_attribute(); attr;
        attr = attr->next_attribute() )
    {
      const string aname = local_name( attr->name(), attr->name_size() );
      const string value( attr->value(), attr->value_size() );

      if( SpecUtils::istarts_with( aname, "on" ) )
        throw runtime_error( "SVG icons can not have event handlers (" + aname + ")." );

      if( (aname == "href") || (aname == "src") )
      {
        string lvalue = value;
        SpecUtils::trim( lvalue );
        SpecUtils::to_lower_ascii( lvalue );
        const bool local = SpecUtils::starts_with( lvalue, "#" );
        const bool raster = SpecUtils::starts_with( lvalue, "data:image/png" )
                            || SpecUtils::starts_with( lvalue, "data:image/jpeg" )
                            || SpecUtils::starts_with( lvalue, "data:image/gif" );
        if( !local && !raster )
          throw runtime_error( "SVG icons can not refer to anything outside themselves." );
      }//if( a link )

      if( value.find( "url(" ) != string::npos )
        check_css_urls( value );
    }//for( loop over attributes )

    for( const rapidxml::xml_node<char> *child = el->first_node(); child;
        child = child->next_sibling() )
    {
      switch( child->type() )
      {
        case rapidxml::node_element:
          check_svg_element( child );
          break;

        case rapidxml::node_data:
        case rapidxml::node_cdata:
          if( name == "style" )
            check_css_urls( string( child->value(), child->value_size() ) );
          break;

        default:
          break;
      }//switch( child->type() )
    }//for( loop over child nodes )
  }//check_svg_element(...)


  /** Parses an SVG viewBox ("x y w h", space or comma separated); throws if it is not usable. */
  array<double,4> parse_view_box( string text )
  {
    for( char &c : text )
    {
      if( c == ',' )
        c = ' ';
    }

    array<double,4> answer{ 0.0, 0.0, 0.0, 0.0 };
    istringstream strm( text );
    for( double &v : answer )
    {
      if( !(strm >> v) || !std::isfinite( v ) )
        throw runtime_error( "The SVG's viewBox is not four numbers." );
    }

    string extra;
    if( (strm >> extra) || (answer[2] <= 0.0) || (answer[3] <= 0.0) )
      throw runtime_error( "The SVG's viewBox must be four numbers, with a positive width and height." );

    return answer;
  }//parse_view_box(...)


  /** The text of an XML element: its value, or the content of its first CDATA section. */
  string element_text( const rapidxml::xml_node<char> *node )
  {
    if( node->value_size() )
      return string( node->value(), node->value_size() );

    for( const rapidxml::xml_node<char> *child = node->first_node(); child;
        child = child->next_sibling() )
    {
      if( child->type() == rapidxml::node_cdata )
        return string( child->value(), child->value_size() );
    }

    return string();
  }//element_text(...)


  uint32_t read_le_uint32( const string &data, const size_t pos )
  {
    const unsigned char *p = reinterpret_cast<const unsigned char *>( data.data() + pos );
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
  }


  /** A binary glTF (.glb) as glTF JSON, its BIN chunk embedded as buffer 0's data URI. */
  nlohmann::json glb_to_json( const string &data )
  {
    if( (data.size() < 20) || (read_le_uint32( data, 4 ) != 2) )
      throw runtime_error( "Only version 2 binary glTF (.glb) files are supported." );

    const size_t length = read_le_uint32( data, 8 );
    if( length > data.size() )
      throw runtime_error( "The .glb file is truncated." );

    string json_chunk, bin_chunk;
    size_t pos = 12;
    while( (pos + 8) <= length )
    {
      const size_t chunk_length = read_le_uint32( data, pos );
      const uint32_t chunk_type = read_le_uint32( data, pos + 4 );
      pos += 8;
      if( (pos + chunk_length) > length )
        throw runtime_error( "The .glb file is truncated." );

      if( chunk_type == 0x4E4F534A )  // "JSON"
        json_chunk = data.substr( pos, chunk_length );
      else if( (chunk_type == 0x004E4942) && bin_chunk.empty() )  // "BIN"
        bin_chunk = data.substr( pos, chunk_length );

      pos += chunk_length;
    }//while( more chunks )

    nlohmann::json doc = nlohmann::json::parse( json_chunk );

    if( !bin_chunk.empty() )
    {
      if( !doc.contains( "buffers" ) || !doc["buffers"].is_array() || doc["buffers"].empty()
          || doc["buffers"][0].contains( "uri" ) )
        throw runtime_error( "The .glb file's binary chunk does not match its buffers." );

      doc["buffers"][0]["uri"] = "data:application/octet-stream;base64,"
                                 + Wt::Utils::base64Encode( bin_chunk, false );
    }//if( !bin_chunk.empty() )

    return doc;
  }//glb_to_json(...)


  size_t gltf_component_size( const int component_type )
  {
    switch( component_type )
    {
      case 5120: case 5121: return 1;  // (unsigned) byte
      case 5122: case 5123: return 2;  // (unsigned) short
      case 5125: case 5126: return 4;  // unsigned int, float
    }
    throw runtime_error( "The glTF model has an unknown accessor component type." );
  }//gltf_component_size(...)


  size_t gltf_type_count( const string &type )
  {
    if( type == "SCALAR" ) return 1;
    if( type == "VEC2" ) return 2;
    if( type == "VEC3" ) return 3;
    if( type == "VEC4" ) return 4;
    if( type == "MAT4" ) return 16;
    throw runtime_error( "The glTF model has an unsupported accessor type (" + type + ")." );
  }//gltf_type_count(...)


  /** Checks an accessor a triangle mesh uses is in range of its data, and of the expected kind. */
  void check_gltf_accessor( const nlohmann::json &doc, const vector<size_t> &buffer_sizes,
                            const size_t index, const bool is_indices )
  {
    const nlohmann::json &accessors = doc.at( "accessors" );
    if( index >= accessors.size() )
      throw runtime_error( "The glTF model refers to a missing accessor." );

    const nlohmann::json &acc = accessors[index];
    if( acc.contains( "sparse" ) || !acc.contains( "bufferView" ) )
      throw runtime_error( "Sparse glTF accessors are not supported." );

    const int component_type = acc.at( "componentType" ).get<int>();
    const string type = acc.at( "type" ).get<string>();
    if( is_indices && ((type != "SCALAR")
                       || ((component_type != 5121) && (component_type != 5123) && (component_type != 5125))) )
      throw runtime_error( "The glTF model's triangle indices must be unsigned integers." );

    if( !is_indices && ((type != "VEC3") || (component_type != 5126)) )
      throw runtime_error( "The glTF model's vertex positions and normals must be (non-quantized) floats." );

    const size_t view_index = acc.at( "bufferView" ).get<size_t>();
    const nlohmann::json &views = doc.at( "bufferViews" );
    if( view_index >= views.size() )
      throw runtime_error( "The glTF model refers to a missing buffer view." );

    const nlohmann::json &view = views[view_index];
    const size_t elem_size = gltf_component_size( component_type ) * gltf_type_count( type );
    if( view.contains( "byteStride" ) && (view["byteStride"].get<size_t>() != elem_size) )
      throw runtime_error( "Interleaved glTF vertex data is not supported; please export the model"
                           " without interleaving." );

    const size_t buffer_index = view.at( "buffer" ).get<size_t>();
    if( buffer_index >= buffer_sizes.size() )
      throw runtime_error( "The glTF model refers to a missing buffer." );

    // Written so that no sum or product can wrap around.
    const size_t buffer_size = buffer_sizes[buffer_index];
    const size_t view_offset = view.value( "byteOffset", size_t(0) );
    const size_t view_length = view.at( "byteLength" ).get<size_t>();
    const size_t acc_offset = acc.value( "byteOffset", size_t(0) );
    const size_t count = acc.at( "count" ).get<size_t>();
    if( (view_offset > buffer_size) || (view_length > (buffer_size - view_offset))
        || (acc_offset > view_length) || (count > ((view_length - acc_offset) / elem_size)) )
      throw runtime_error( "The glTF model's data is shorter than its accessors say." );
  }//check_gltf_accessor(...)


  /** A glTF index (or count): a non-negative integer, else throws. */
  size_t gltf_index( const nlohmann::json &v )
  {
    if( !v.is_number_integer() || (v.get<int64_t>() < 0) )
      throw runtime_error( "The glTF model has an invalid index or count." );
    return v.get<size_t>();
  }


  /** `n` finite numbers (or, with n == 0, one to sixteen), else throws. */
  nlohmann::json gltf_numbers( const nlohmann::json &v, const size_t n )
  {
    if( !v.is_array() || (n && (v.size() != n)) || v.empty() || (v.size() > 16) )
      throw runtime_error( "The glTF model has a malformed number list." );
    nlohmann::json answer = nlohmann::json::array();
    for( const nlohmann::json &x : v )
    {
      if( !x.is_number() || !std::isfinite( x.get<double>() ) )
        throw runtime_error( "The glTF model has a malformed number list." );
      answer.push_back( x.get<double>() );
    }
    return answer;
  }


  /** The glTF rebuilt from only what InterSpec draws, each value type-checked: no names, extras,
   extensions, or textures, so neither a deep nesting nor arbitrary text in a file from elsewhere
   reaches the serializer or the browser.  Buffers are decoded and re-encoded. */
  nlohmann::json rebuilt_gltf( const nlohmann::json &doc )
  {
    const auto list = []( const nlohmann::json &parent, const char *key ) -> const nlohmann::json & {
      static const nlohmann::json empty = nlohmann::json::array();
      if( !parent.is_object() || !parent.contains( key ) )
        return empty;
      if( !parent[key].is_array() )
        throw runtime_error( string( "The glTF model's " ) + key + " is not a list." );
      return parent[key];
    };
    const auto copy_index = []( const nlohmann::json &from, nlohmann::json &to, const char *key ) {
      if( from.contains( key ) )
        to[key] = gltf_index( from[key] );
    };

    nlohmann::json clean;
    clean["asset"] = { {"version", "2.0"} };
    if( doc.contains( "scene" ) )
      clean["scene"] = gltf_index( doc["scene"] );

    clean["scenes"] = nlohmann::json::array();
    for( const nlohmann::json &scene : list( doc, "scenes" ) )
    {
      nlohmann::json roots = nlohmann::json::array();
      for( const nlohmann::json &n : list( scene, "nodes" ) )
        roots.push_back( gltf_index( n ) );
      clean["scenes"].push_back( { {"nodes", roots} } );
    }

    clean["nodes"] = nlohmann::json::array();
    for( const nlohmann::json &node : list( doc, "nodes" ) )
    {
      nlohmann::json c = nlohmann::json::object();
      copy_index( node, c, "mesh" );
      if( node.contains( "children" ) )
      {
        c["children"] = nlohmann::json::array();
        for( const nlohmann::json &child : list( node, "children" ) )
          c["children"].push_back( gltf_index( child ) );
      }
      static const pair<const char *, size_t> transforms[] = { { "matrix", 16 }, { "translation", 3 },
                                                               { "rotation", 4 }, { "scale", 3 } };
      for( const pair<const char *, size_t> &t : transforms )
      {
        if( node.contains( t.first ) )
          c[t.first] = gltf_numbers( node[t.first], t.second );
      }
      clean["nodes"].push_back( c );
    }//for( each node )

    clean["meshes"] = nlohmann::json::array();
    for( const nlohmann::json &mesh : list( doc, "meshes" ) )
    {
      nlohmann::json prims = nlohmann::json::array();
      for( const nlohmann::json &prim : list( mesh, "primitives" ) )
      {
        nlohmann::json c = nlohmann::json::object(), attributes = nlohmann::json::object();
        if( prim.contains( "attributes" ) && prim["attributes"].is_object() )
        {
          copy_index( prim["attributes"], attributes, "POSITION" );
          copy_index( prim["attributes"], attributes, "NORMAL" );
        }
        c["attributes"] = attributes;
        copy_index( prim, c, "indices" );
        copy_index( prim, c, "mode" );
        copy_index( prim, c, "material" );
        prims.push_back( c );
      }
      clean["meshes"].push_back( { {"primitives", prims} } );
    }//for( each mesh )

    clean["materials"] = nlohmann::json::array();
    for( const nlohmann::json &mat : list( doc, "materials" ) )
    {
      nlohmann::json pbr = nlohmann::json::object();
      if( mat.contains( "pbrMetallicRoughness" ) && mat["pbrMetallicRoughness"].is_object()
          && mat["pbrMetallicRoughness"].contains( "baseColorFactor" ) )
        pbr["baseColorFactor"] = gltf_numbers( mat["pbrMetallicRoughness"]["baseColorFactor"], 4 );
      clean["materials"].push_back( { {"pbrMetallicRoughness", pbr} } );
    }

    clean["accessors"] = nlohmann::json::array();
    for( const nlohmann::json &acc : list( doc, "accessors" ) )
    {
      if( !acc.is_object() || acc.contains( "sparse" ) )
        throw runtime_error( "Sparse glTF accessors are not supported." );
      nlohmann::json c = nlohmann::json::object();
      copy_index( acc, c, "bufferView" );
      copy_index( acc, c, "byteOffset" );
      c["componentType"] = gltf_index( acc.at( "componentType" ) );
      c["count"] = gltf_index( acc.at( "count" ) );
      c["type"] = acc.at( "type" ).get<string>();
      gltf_type_count( c["type"].get<string>() );  //throws if not a type we know
      for( const char *key : { "min", "max" } )
      {
        if( acc.contains( key ) )
          c[key] = gltf_numbers( acc[key], 0 );
      }
      clean["accessors"].push_back( c );
    }//for( each accessor )

    clean["bufferViews"] = nlohmann::json::array();
    for( const nlohmann::json &view : list( doc, "bufferViews" ) )
    {
      nlohmann::json c = nlohmann::json::object();
      c["buffer"] = gltf_index( view.at( "buffer" ) );
      c["byteLength"] = gltf_index( view.at( "byteLength" ) );
      copy_index( view, c, "byteOffset" );
      copy_index( view, c, "byteStride" );
      clean["bufferViews"].push_back( c );
    }

    // Every buffer must be embedded, so the model can travel inside the DRF.
    clean["buffers"] = nlohmann::json::array();
    for( const nlohmann::json &buffer : list( doc, "buffers" ) )
    {
      const string uri = (buffer.is_object() && buffer.contains( "uri" ) && buffer["uri"].is_string())
                         ? buffer["uri"].get<string>() : string();
      const size_t comma = uri.find( ',' );
      if( !SpecUtils::starts_with( uri, "data:" ) || (comma == string::npos)
          || (uri.substr( 0, comma ).find( ";base64" ) == string::npos) )
        throw runtime_error( "The model's data is in a separate file; please export it as a single"
                             " .glb, or as glTF with embedded buffers." );

      const string data = Wt::Utils::base64Decode( uri.substr( comma + 1 ) );
      if( data.size() < gltf_index( buffer.at( "byteLength" ) ) )
        throw runtime_error( "The model's embedded data is shorter than it says." );
      clean["buffers"].push_back( { {"byteLength", data.size()},
                                    {"uri", "data:application/octet-stream;base64," + Wt::Utils::base64Encode( data, false )} } );
    }//for( each buffer )

    // A node may have only one parent (as glTF requires): then walking the tree visits each node once.
    vector<int> parents( clean["nodes"].size(), 0 );
    for( const nlohmann::json &node : clean["nodes"] )
    {
      for( const nlohmann::json &child : node.value( "children", nlohmann::json::array() ) )
      {
        const size_t i = child.get<size_t>();
        if( (i >= parents.size()) || (++parents[i] > 1) )
          throw runtime_error( "The glTF model's node tree is malformed." );
      }
    }

    return clean;
  }//rebuilt_gltf(...)


  typedef array<double,16> Mat4;  // column-major, as glTF stores them

  Mat4 mat_mult( const Mat4 &a, const Mat4 &b )
  {
    Mat4 answer{};
    for( size_t col = 0; col < 4; ++col )
    {
      for( size_t row = 0; row < 4; ++row )
      {
        double sum = 0.0;
        for( size_t k = 0; k < 4; ++k )
          sum += a[k*4 + row] * b[col*4 + k];
        answer[col*4 + row] = sum;
      }
    }
    return answer;
  }//mat_mult(...)


  /** A glTF node's local transform, from its `matrix` or its translation/rotation/scale. */
  Mat4 node_matrix( const nlohmann::json &node )
  {
    Mat4 answer{ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    if( node.contains( "matrix" ) )
    {
      const vector<double> m = node["matrix"].get<vector<double>>();
      if( m.size() == 16 )
        std::copy( begin(m), end(m), begin(answer) );
      return answer;
    }

    const vector<double> t = node.value( "translation", vector<double>{ 0.0, 0.0, 0.0 } );
    const vector<double> q = node.value( "rotation", vector<double>{ 0.0, 0.0, 0.0, 1.0 } );
    const vector<double> s = node.value( "scale", vector<double>{ 1.0, 1.0, 1.0 } );
    if( (t.size() != 3) || (q.size() != 4) || (s.size() != 3) )
      return answer;

    const double x = q[0], y = q[1], z = q[2], w = q[3];
    answer = Mat4{ (1 - 2*(y*y + z*z))*s[0], (2*(x*y + z*w))*s[0],     (2*(x*z - y*w))*s[0],     0,
                   (2*(x*y - z*w))*s[1],     (1 - 2*(x*x + z*z))*s[1], (2*(y*z + x*w))*s[1],     0,
                   (2*(x*z + y*w))*s[2],     (2*(y*z - x*w))*s[2],     (1 - 2*(x*x + y*y))*s[2], 0,
                   t[0],                     t[1],                     t[2],                     1 };
    return answer;
  }//node_matrix(...)


  /** Grows `bounds` by the transformed POSITION boxes of a node's meshes, and its children's. */
  void add_node_bounds( const nlohmann::json &doc, const size_t node_index, const Mat4 &parent,
                        array<double,6> &bounds, vector<bool> &visited )
  {
    // Each node once: a crafted file listing the same children over and over (or a cycle) could
    //  otherwise make this walk take forever.
    const nlohmann::json &nodes = doc.at( "nodes" );
    if( (node_index >= nodes.size()) || visited.at( node_index ) )
      return;
    visited[node_index] = true;

    const nlohmann::json &node = nodes[node_index];
    const Mat4 world = mat_mult( parent, node_matrix( node ) );

    if( node.contains( "mesh" ) )
    {
      const nlohmann::json &mesh = doc.at( "meshes" ).at( node["mesh"].get<size_t>() );
      for( const nlohmann::json &prim : mesh.at( "primitives" ) )
      {
        if( !prim.contains( "attributes" ) || !prim["attributes"].contains( "POSITION" ) )
          continue;

        const nlohmann::json &acc = doc.at( "accessors" ).at( prim["attributes"]["POSITION"].get<size_t>() );
        if( !acc.contains( "min" ) || !acc.contains( "max" ) )
          continue;

        const vector<double> lo = acc["min"].get<vector<double>>(), hi = acc["max"].get<vector<double>>();
        if( (lo.size() != 3) || (hi.size() != 3) )
          continue;

        for( int corner = 0; corner < 8; ++corner )
        {
          const double p[3] = { (corner & 1) ? hi[0] : lo[0], (corner & 2) ? hi[1] : lo[1],
                                (corner & 4) ? hi[2] : lo[2] };
          for( size_t i = 0; i < 3; ++i )
          {
            const double v = world[i] * p[0] + world[4 + i] * p[1] + world[8 + i] * p[2] + world[12 + i];
            bounds[i] = std::min( bounds[i], v );
            bounds[3 + i] = std::max( bounds[3 + i], v );
          }
        }//for( each corner of the box )
      }//for( each primitive )
    }//if( node.contains( "mesh" ) )

    if( node.contains( "children" ) )
    {
      for( const nlohmann::json &child : node["children"] )
        add_node_bounds( doc, child.get<size_t>(), world, bounds, visited );
    }
  }//add_node_bounds(...)
}//namespace


bool DrfAppearance::empty() const
{
  return icon_svg.empty() && model_gltf.empty();
}


bool DrfAppearance::operator==( const DrfAppearance &rhs ) const
{
  return (icon_svg == rhs.icon_svg) && (model_gltf == rhs.model_gltf);
}


string DrfAppearance::icon_data_uri() const
{
  if( icon_svg.empty() )
    return string();
  return "data:image/svg+xml;base64," + Wt::Utils::base64Encode( icon_svg, false );
}//icon_data_uri()


string DrfAppearance::icon_html( const double height_px, const double max_width_px ) const
{
  if( icon_svg.empty() )
    return string();

  double width_px = height_px;
  try
  {
    const array<double,4> vb = icon_view_box();
    width_px = std::min( max_width_px, height_px * vb[2] / vb[3] );
  }catch( std::exception & )
  {
  }

  const string uri = icon_data_uri();
  char size[64] = { '\0' };
  snprintf( size, sizeof(size), "width:%.1fpx;height:%.1fpx;", width_px, height_px );

  return "<span class=\"DrfIcon\" style=\"" + string(size) + "-webkit-mask-image:url('" + uri
         + "');mask-image:url('" + uri + "');\"></span>";
}//icon_html(...)


array<double,4> DrfAppearance::icon_view_box() const
{
  if( icon_svg.empty() )
    throw runtime_error( "No icon" );

  vector<char> buffer( begin(icon_svg), end(icon_svg) );
  buffer.push_back( '\0' );
  rapidxml::xml_document<char> doc;
  doc.parse<0>( buffer.data() );

  const rapidxml::xml_node<char> *root = doc.first_node();
  const rapidxml::xml_attribute<char> *vb = root ? root->first_attribute( "viewBox", 7 ) : nullptr;
  if( !vb )
    throw runtime_error( "The icon has no viewBox." );

  return parse_view_box( string( vb->value(), vb->value_size() ) );
}//icon_view_box()


bool DrfAppearance::model_bounds_m( array<double,6> &bounds ) const
{
  if( model_gltf.empty() )
    return false;

  try
  {
    const nlohmann::json doc = nlohmann::json::parse( model_gltf );
    const double inf = std::numeric_limits<double>::infinity();
    array<double,6> answer{ inf, inf, inf, -inf, -inf, -inf };
    const Mat4 identity{ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    vector<bool> visited( doc.value( "nodes", nlohmann::json::array() ).size(), false );

    if( doc.contains( "scenes" ) && !doc["scenes"].empty() )
    {
      const size_t scene = doc.value( "scene", size_t(0) );
      for( const nlohmann::json &node : doc["scenes"].at( scene ).value( "nodes", nlohmann::json::array() ) )
        add_node_bounds( doc, node.get<size_t>(), identity, answer, visited );
    }else if( doc.contains( "nodes" ) )
    {
      // No scene: every node that is nobody's child is a root.
      set<size_t> children;
      for( const nlohmann::json &node : doc["nodes"] )
      {
        for( const nlohmann::json &child : node.value( "children", nlohmann::json::array() ) )
          children.insert( child.get<size_t>() );
      }
      for( size_t i = 0; i < doc["nodes"].size(); ++i )
      {
        if( !children.count( i ) )
          add_node_bounds( doc, i, identity, answer, visited );
      }
    }//if( scenes ) / else

    if( !(answer[0] <= answer[3]) )
      return false;

    bounds = answer;
    return true;
  }catch( std::exception &e )
  {
    cerr << "DrfAppearance::model_bounds_m: " << e.what() << endl;
  }

  return false;
}//model_bounds_m(...)


string DrfAppearance::validated_svg( const string &svg )
{
  if( svg.size() > ns_max_svg_bytes )
    throw runtime_error( "The SVG is larger than 64 kB; please simplify it, or save it as a plain or"
                         " optimized SVG." );

  string lower = svg;
  SpecUtils::to_lower_ascii( lower );
  if( lower.find( "<!entity" ) != string::npos )
    throw runtime_error( "SVG icons can not define XML entities." );

  if( max_element_depth( svg ) > ns_max_svg_depth )
    throw runtime_error( "The SVG nests its elements too deeply." );

  vector<char> buffer( begin(svg), end(svg) );
  buffer.push_back( '\0' );

  // Default flags skip the XML declaration, DOCTYPE, comments, and processing instructions, so they
  //  are not written back out below.
  rapidxml::xml_document<char> doc;
  try
  {
    doc.parse<0>( buffer.data() );
  }catch( std::exception &e )
  {
    throw runtime_error( "The file is not valid XML (" + string( e.what() ) + ")." );
  }

  rapidxml::xml_node<char> *root = doc.first_node();
  if( !root || (local_name( root->name(), root->name_size() ) != "svg") || root->next_sibling() )
    throw runtime_error( "The file is not an SVG image." );

  const rapidxml::xml_attribute<char> *vb = root->first_attribute( "viewBox", 7 );
  if( !vb )
    throw runtime_error( "The SVG has no viewBox; InterSpec needs one to know its size (in mm)"
                         " and where the detector face is." );
  parse_view_box( string( vb->value(), vb->value_size() ) );  //throws if not usable

  strip_editor_data( root );
  check_svg_element( root );

  // An SVG shown as an image is only rendered with the SVG namespace declared.
  if( !root->first_attribute( "xmlns", 5 ) && (root->name_size() == 3) )
    root->append_attribute( doc.allocate_attribute( "xmlns", "http://www.w3.org/2000/svg" ) );

  string answer;
  rapidxml::print( std::back_inserter( answer ), doc, rapidxml::print_no_indenting );

  return answer;
}//validated_svg(...)


string DrfAppearance::validated_gltf( const string &gltf_or_glb )
{
  if( gltf_or_glb.size() > ns_max_gltf_bytes )
    throw runtime_error( "The model is larger than 2 MB; please reduce its number of triangles." );

  nlohmann::json doc;
  try
  {
    const bool is_glb = (gltf_or_glb.size() >= 4) && (gltf_or_glb.compare( 0, 4, "glTF" ) == 0);
    doc = is_glb ? glb_to_json( gltf_or_glb ) : nlohmann::json::parse( gltf_or_glb );
  }catch( nlohmann::json::exception &e )
  {
    throw runtime_error( "The file is not a glTF model (" + string( e.what() ) + ")." );
  }

  nlohmann::json clean;
  try
  {
    if( !doc.is_object() || !doc.contains( "asset" ) || !doc["asset"].is_object()
        || !doc["asset"].contains( "version" ) || !doc["asset"]["version"].is_string() )
      throw runtime_error( "The file is not a glTF model." );

    const string version = doc["asset"]["version"].get<string>();
    if( !SpecUtils::starts_with( version, "2." ) )
      throw runtime_error( "Only glTF 2.0 models are supported (this one is version "
                           + version.substr( 0, 16 ) + ")." );

    if( doc.contains( "extensionsRequired" ) && doc["extensionsRequired"].is_array()
        && !doc["extensionsRequired"].empty() )
    {
      const nlohmann::json &ext = doc["extensionsRequired"][0];
      throw runtime_error( "The model needs the glTF extension "
                           + (ext.is_string() ? ext.get<string>().substr( 0, 64 ) : string( "?" ))
                           + ", which InterSpec can not read; please export it without compression." );
    }

    clean = rebuilt_gltf( doc );
    doc = nlohmann::json();   //frees the original, which may be deeply nested, iteratively

    vector<size_t> buffer_sizes;
    for( const nlohmann::json &buffer : clean["buffers"] )
      buffer_sizes.push_back( buffer["byteLength"].get<size_t>() );

    size_t num_triangle_primitives = 0;
    for( const nlohmann::json &mesh : clean["meshes"] )
    {
      for( const nlohmann::json &prim : mesh["primitives"] )
      {
        if( prim.value( "mode", size_t(4) ) != 4 )
          continue;  //points and lines are not drawn

        const nlohmann::json &attributes = prim["attributes"];
        if( !attributes.contains( "POSITION" ) )
          continue;

        check_gltf_accessor( clean, buffer_sizes, attributes["POSITION"].get<size_t>(), false );
        if( attributes.contains( "NORMAL" ) )
          check_gltf_accessor( clean, buffer_sizes, attributes["NORMAL"].get<size_t>(), false );
        if( prim.contains( "indices" ) )
          check_gltf_accessor( clean, buffer_sizes, prim["indices"].get<size_t>(), true );

        ++num_triangle_primitives;
      }//for( each primitive )
    }//for( each mesh )

    if( !num_triangle_primitives )
      throw runtime_error( "The model has no triangle meshes." );
  }catch( nlohmann::json::exception &e )
  {
    throw runtime_error( "The glTF model is malformed (" + string( e.what() ) + ")." );
  }

  return clean.dump();
}//validated_gltf(...)


void DrfAppearance::toXml( ::rapidxml::xml_node<char> *parent, ::rapidxml::xml_document<char> *doc ) const
{
  if( empty() )
    return;

  rapidxml::xml_node<char> *base_node = doc->allocate_node( rapidxml::node_element, "Appearance" );
  base_node->append_attribute( doc->allocate_attribute( "version", "0" ) );
  parent->append_node( base_node );

  const auto add_text = [doc,base_node]( const char *name, const char *format, const string &text ){
    const char *val = doc->allocate_string( text.c_str(), text.size() + 1 );
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, name, val,
                                                         strlen( name ), text.size() );
    node->append_attribute( doc->allocate_attribute( "format", format ) );
    base_node->append_node( node );
  };

  if( !icon_svg.empty() )
    add_text( "Icon", "svg", icon_svg );

  if( !model_gltf.empty() )
    add_text( "ExteriorModel", "gltf", model_gltf );
}//toXml(...)


shared_ptr<const DrfAppearance> DrfAppearance::fromXml( const ::rapidxml::xml_node<char> *node )
{
  if( !node )
    return nullptr;

  // A format this build doesnt know (e.g., a later addition) is skipped, not an error.
  const auto has_format = []( const rapidxml::xml_node<char> *el, const char *format ) -> bool {
    const rapidxml::xml_attribute<char> *attr = el->first_attribute( "format", 6 );
    return !attr || (string( attr->value(), attr->value_size() ) == format);
  };

  auto answer = make_shared<DrfAppearance>();

  const rapidxml::xml_node<char> *icon = node->first_node( "Icon", 4 );
  if( icon && has_format( icon, "svg" ) )
  {
    // Re-checked on every read: the XML may have come from anyone.
    try
    {
      answer->icon_svg = validated_svg( element_text( icon ) );
    }catch( std::exception &e )
    {
      cerr << "DrfAppearance::fromXml: ignoring the DRF's icon: " << e.what() << endl;
    }
  }//if( icon )

  const rapidxml::xml_node<char> *model = node->first_node( "ExteriorModel", 13 );
  if( model && has_format( model, "gltf" ) )
  {
    try
    {
      answer->model_gltf = validated_gltf( element_text( model ) );
    }catch( std::exception &e )
    {
      cerr << "DrfAppearance::fromXml: ignoring the DRF's exterior model: " << e.what() << endl;
    }
  }//if( model )

  if( answer->empty() )
    return nullptr;

  return answer;
}//fromXml(...)
