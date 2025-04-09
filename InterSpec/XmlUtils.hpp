#ifndef XmlUtils_hpp
#define XmlUtils_hpp
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

#include <cstdio>
#include <string>
#include <assert.h>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"


/** Some XML helper functions for serializing to/from XML.
 */
namespace XmlUtils
{
  
  inline void append_float_node( rapidxml::xml_node<char> *base_node, const char *node_name, const double value )
  {
    assert( base_node && base_node->document() );
    rapidxml::xml_document<char> *doc = base_node->document();
    
    char buffer[128];
    snprintf( buffer, sizeof(buffer), "%1.8e", value );
    const char *strvalue = doc->allocate_string( buffer );
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, node_name, strvalue );
    base_node->append_node( node );
  }//append_float_node(...)


  inline void append_bool_node( rapidxml::xml_node<char> *base_node, const char *node_name, const bool value )
  {
    assert( base_node && base_node->document() );
    rapidxml::xml_document<char> *doc = base_node->document();
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, node_name, (value ? "true" : "false") );
    base_node->append_node( node );
  }//append_bool_node(...)


  template <class T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr >
  inline void append_int_node( rapidxml::xml_node<char> *base_node, const char *node_name, const T value )
  {
    assert( base_node && base_node->document() );
    rapidxml::xml_document<char> *doc = base_node->document();
    const char *strvalue = doc->allocate_string( std::to_string(value).c_str() );
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, node_name, strvalue );
    base_node->append_node( node );
  }//append_bool_node(...)


  inline rapidxml::xml_node<char> *append_string_node( rapidxml::xml_node<char> *base_node, const char *node_name, const std::string &value )
  {
    assert( base_node && base_node->document() );
    rapidxml::xml_document<char> *doc = base_node->document();
    const char *strvalue = doc->allocate_string( value.c_str(), value.size() + 1 );
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, node_name, strvalue );
    base_node->append_node( node );
    
    return node;
  }//append_string_node(...)


  inline void append_attrib( rapidxml::xml_node<char> *base_node, const std::string &name, const std::string &value )
  {
    assert( base_node && base_node->document() );
    rapidxml::xml_document<char> *doc = base_node->document();
    
    const char *name_str = doc->allocate_string( name.c_str() );
    const char *value_str = doc->allocate_string( value.c_str() );
    
    rapidxml::xml_attribute<char> *attrib = doc->allocate_attribute( name_str, value_str );
    base_node->append_attribute( attrib );
  }//void append_attrib(...)


  inline void append_version_attrib( rapidxml::xml_node<char> *base_node, const int version )
  {
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%i", version );
    append_attrib( base_node, "version", buffer );
  }//append_version_attrib(...)


  template<size_t n>
  inline double get_float_node_value( const rapidxml::xml_node<char> * const parent_node, const char (&name)[n] )
  {
    assert( parent_node );
    assert( name );
    
    if( !parent_node )
      throw std::runtime_error( "null parent node." );
    
    const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent_node, name);
    if( !node )
      throw std::runtime_error( "Missing node '" + std::string(name) + "'" );
    
    double answer;
    if( !SpecUtils::parse_double(node->value(), node->value_size(), answer) )
      throw std::runtime_error( "Value ('" + SpecUtils::xml_value_str(node) + "') of node '"
                          + std::string(name) + "' was a valid float." );
    
    return answer;
  }//double get_float_node_value(...)


  template<size_t n>
  inline bool get_bool_node_value( const rapidxml::xml_node<char> * const parent_node, const char (&name)[n] )
  {
    assert( parent_node );
    assert( name );
    
    if( !parent_node )
      throw std::runtime_error( "null parent node." );
    
    const rapidxml::xml_node<char> *node = XML_FIRST_NODE(parent_node, name);
    if( !node )
      throw std::runtime_error( "Missing node '" + std::string(name) + "'" );
    
    if( XML_VALUE_ICOMPARE(node,"yes") || XML_VALUE_ICOMPARE(node,"true") || XML_VALUE_ICOMPARE(node,"1") )
      return true;
    
    if( !XML_VALUE_ICOMPARE(node,"no") && !XML_VALUE_ICOMPARE(node,"false") && !XML_VALUE_ICOMPARE(node,"0") )
      throw std::runtime_error( "Invalid boolean value in node '" + std::string(name) + "' with value '"
                          + SpecUtils::xml_value_str(node) );
    
    return false;
  }//bool get_bool_node_value(...)

  template<size_t n>
  inline int get_int_attribute( const rapidxml::xml_node<char> * const node, const char (&name)[n] )
  {
    assert( node );
    assert( name );
    
    const rapidxml::xml_attribute<char> *att = XML_FIRST_ATTRIB(node, name);
    if( !att )
      throw std::runtime_error( "Missing attribute '" + std::string(name) + "'" );
    
    int answer;
    if( !SpecUtils::parse_int(att->value(), att->value_size(), answer) )
      throw std::runtime_error( "Invalid integer value in attribute '" + std::string(name) + "' with value '"
                          + SpecUtils::xml_value_str(att) + "'" );
    
    return answer;
  }//int get_int_attribute(...)

  inline void check_xml_version( const rapidxml::xml_node<char> * const node, const int required_version )
  {
    assert( node );
    const rapidxml::xml_attribute<char> *att = XML_FIRST_ATTRIB(node, "version");
    
    const int version = get_int_attribute( node, "version" );
    
    if( (version < 0) || (version > required_version) )
      throw std::runtime_error( "Invalid version: " + std::to_string(version) + ".  "
                          + "Only up to version " + std::to_string(required_version)
                          + " supported." );
  }//check_xml_version(...)

  template<size_t n>
  inline const rapidxml::xml_node<char> *get_required_node( const rapidxml::xml_node<char> *parent, const char (&name)[n] )
  {
    assert( parent );
    const auto child_node = XML_FIRST_INODE(parent, name);
    if( !child_node )
      throw std::runtime_error( "No <" + std::string(name) + "> node" );
    
    return child_node;
  }//get_required_node(...)
   
}//namespace XmlUtils

#endif //XmlUtils_hpp
