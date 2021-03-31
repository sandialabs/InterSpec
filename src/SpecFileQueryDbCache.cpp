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

#include <Wt/Utils>
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>

#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/Dbo/backend/Sqlite3>

#include "pugixml.hpp"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/SpecFileQuery.h"
//#include "InterSpec/InterSpecApp.h" //for passMessage debugging
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/SpecFileQueryDbCache.h"


#include <boost/config.hpp>


#if( defined(BOOST_NO_CXX11_HDR_CODECVT) )
//gcc 4.8 doesnt have enable_if_t or codecvt...
#define HAS_STD_ENABLE_IF 0
#include <boost/core/enable_if.hpp>
#else
#define HAS_STD_ENABLE_IF 1
#endif

using namespace std;
using namespace Wt;

#if(WT_VERSION<=0x3030400)
//Arrrg: Wt 3.3.4 at least has bug (fixed by 4.0.3) where when you bind a blob
// field, SQLITE_TRANSIENT option is not used (SQLITE_STATIC instead) so your
// vector<unsigned char> would have to stick around until after writing to
// database, which we cant do here.
// So instead we have to write text to the database (and actually we cant just
// put blob data into text, because then when reading back from database, we
// will only get to the first zero byte), which is probably way slower, and was
// a pain to discover this issue and implement around it.
#include <sstream>
#include <type_traits>
#include <Wt/WStringStream>
#include <boost/tokenizer.hpp>
#include <boost/io/detail/quoted_manip.hpp>
#define USE_TEXT_FOR_BLOB 1
#else
#define USE_TEXT_FOR_BLOB 0
static_assert( 0, "Serialization to database hasnt been checked (but is probably okay, assuming this version of Wt doesnt have the bug 3.3.4 does)" );
#endif

namespace
{
  /*
  struct xml_string_writer: pugi::xml_writer
  {
    std::string result;
    
    virtual void write(const void* data, size_t size)
    {
      result.append(static_cast<const char*>(data), size);
    }
  };//struct xml_string_writer:
  */
  
  bool xml_files_small_enough( const std::string &filename, void *userdata )
  {
    if( !SpecUtils::iends_with(filename, ".xml") )
      return false;
    
    const size_t fsize = SpecUtils::file_size(filename);
    if( (fsize < SpecFileQuery::EventXmlTest::sm_min_event_xml_file_size)
       || (fsize > SpecFileQuery::EventXmlTest::sm_max_event_xml_file_size) )
      return false;
    
    return true;
  }//xml_files_small_enough(...)
}//namespace


namespace Wt {
  namespace Dbo {
    template<>
    struct dbo_traits<SpecFileInfoToQuery> : dbo_default_traits {
      static const char *surrogateIdField() { return nullptr; }
    };
    // Necessary if you want to use ptr<const SpecFileInfoToQuery>
    template<> struct dbo_traits<const SpecFileInfoToQuery> : dbo_traits<SpecFileInfoToQuery> {};

//How to specialize types for storing in the database.  If we want to get rid of
//  bitely and store each field in the table manually, we need to specialize:
//  bool, set<float>, size_t, set<size_t>
//  set<SpecUtils::EnergyCalType>, set<string>.
//DetectorAnalysis can probably be serialized to XML to store in a text field.

    std::string sql_value_traits<size_t>::type(SqlConnection *conn, int size)
    {
      return conn->longLongType() + " not null";
    }
    void sql_value_traits<size_t>::bind(size_t v, SqlStatement *statement, int column, int size)
    {
      statement->bind(column, static_cast<long long>(v) );
    }
    bool sql_value_traits<size_t>::read(size_t& v, SqlStatement *statement, int column, int size)
    {
      long long intValue;
      bool result = statement->getResult( column, &intValue );
      if( result )
        v = intValue;
      return result;
    }
    
    //We dont have to worry about portability, so
    template<class T>
    void bind_pod_set(const std::set<T> &v, SqlStatement *statement, int column, int size)
    {
      static_assert( std::is_pod<T>::value, "T must be POD" );
      
#if( USE_TEXT_FOR_BLOB )
      static_assert( std::is_integral<T>::value || std::is_enum<T>::value || std::is_floating_point<T>::value, "T must be int or float" );
      
      WStringStream strm;
      bool first = true;
      for( const auto &val : v )
      {
        if( std::is_integral<T>::value || std::is_enum<T>::value )
          strm << string(first ? "" : ",") << static_cast<long long>(val);
        else
          strm << string(first ? "" : ",") << static_cast<double>(val);
        first = false;
      }
      statement->bind( column, strm.str() );
#else
      std::vector<unsigned char> bindata( sizeof(T)*v.size() );

      size_t pos = 0;
      for( const T &f : v )
        memcpy( &(bindata[sizeof(T)*pos++]), &f, sizeof(T) );
      statement->bind( column, bindata );
#endif
    }
    
#if( USE_TEXT_FOR_BLOB )
    bool do_splitting( const std::string &data, set<float> &v )
    {
      vector<float> valfloats;
      if( !SpecUtils::split_to_floats( data.c_str(), data.size(), valfloats ) )
      {
        cerr << "do_splitting: Failed to split to floats" << endl;
        return false;
      }
      for( const auto &i : valfloats )
        v.insert( i );
      return true;
    }//do_splitting( float )


#if( HAS_STD_ENABLE_IF )
    template<class T,
      typename = std::enable_if_t<std::is_integral<T>::value || std::is_enum<T>::value> >
#else
    template<class T,
    typename = boost::enable_if_c<boost::is_integral<T>::value || boost::is_enum<T>::value> >
#endif
    bool do_splitting( const std::string &data, set<T> &v )
    {
      vector<long long> valints;
      if( !SpecUtils::split_to_long_longs( data.c_str(), data.size(), valints ) )
      {
        cerr << "read_pod_set: Failed to split to long long" << endl;
        return false;
      }
      for( const auto &i : valints )
        v.insert( static_cast<T>(i) );
      return true;
    }//do_splitting( integral type )

#endif //#if( USE_TEXT_FOR_BLOB )

    template<class T>
    bool read_pod_set(std::set<T> &v, SqlStatement *statement, int column, int size)
    {
      static_assert( std::is_pod<T>::value, "T must be POD" );

      v.clear();
#if( USE_TEXT_FOR_BLOB )
      string data;
      if( !statement->getResult( column, &data, size ) )
        return false;
      return do_splitting( data, v );
#else
      std::vector<unsigned char> bindata;
      
      bool result = statement->getResult( column, &bindata, size );
      result = (result && ((bindata.size() % sizeof(T))==0));
      
      //ToDo: remove below assert after testing.
      if( result && (((bindata.size() % sizeof(T))!=0)) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Database re-read in of a set<POD> did not have proper size!  Programming logic error." );
#endif
        result = false;
      }
      
      if( result )
      {
        T val;
        for( size_t i = 0; i < bindata.size(); i += sizeof(T) )
        {
          memcpy( &val, &(bindata[i]), sizeof(T));
          v.insert( val );
        }
      }
      return result;
#endif
    }
  
    
    std::string sql_value_traits<std::set<float>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      return conn->blobType();
#endif
    }
    void sql_value_traits<std::set<float>>::bind(const std::set<float> &v, SqlStatement *statement, int column, int size)
    {
      bind_pod_set( v, statement, column, size );
    }
    bool sql_value_traits<std::set<float>>::read(std::set<float> &v, SqlStatement *statement, int column, int size)
    {
      return read_pod_set( v, statement, column, size );
    }
    
    
    std::string sql_value_traits<std::set<size_t>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      return conn->blobType();
#endif
    }
    void sql_value_traits<std::set<size_t>>::bind(const std::set<size_t> &v, SqlStatement *statement, int column, int size)
    {
      bind_pod_set( v, statement, column, size );
    }
    bool sql_value_traits<std::set<size_t>>::read(std::set<size_t> &v, SqlStatement *statement, int column, int size)
    {
      return read_pod_set( v, statement, column, size );
    }
    
    
    std::string sql_value_traits<std::set<SpecUtils::EnergyCalType>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      return conn->blobType();
#endif
    }
    void sql_value_traits<std::set<SpecUtils::EnergyCalType>>::bind(const std::set<SpecUtils::EnergyCalType> &v, SqlStatement *statement, int column, int size)
    {
      bind_pod_set( v, statement, column, size );
    }
    bool sql_value_traits<std::set<SpecUtils::EnergyCalType>>::read(std::set<SpecUtils::EnergyCalType> &v, SqlStatement *statement, int column, int size)
    {
      return read_pod_set( v, statement, column, size );
    }

    
    std::string sql_value_traits<std::set<std::time_t>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      return conn->blobType();
#endif
    }
    void sql_value_traits<std::set<std::time_t>>::bind(const std::set<std::time_t> &v, SqlStatement *statement, int column, int size)
    {
      bind_pod_set( v, statement, column, size );
    }
    bool sql_value_traits<std::set<std::time_t>>::read(std::set<std::time_t> &v, SqlStatement *statement, int column, int size)
    {
      return read_pod_set( v, statement, column, size );
    }
    
    
    std::string sql_value_traits<std::set<std::string>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      return conn->blobType();
#endif
    }
    void sql_value_traits<std::set<std::string>>::bind(const std::set<std::string> &v, SqlStatement *statement, int column, int size)
    {
      size_t totallen = sizeof(size_t)*v.size();
      for( const string &str : v )
        totallen += str.size();

#if( USE_TEXT_FOR_BLOB )
      //In principal it looks like spectrum file almost never contain quotes,
      //  but JIC we will go through the hassle of quoting things
      bool first = true;
      stringstream strm;
      for( const auto &val : v )
      {
        if( val.find_last_of("\",") != string::npos )
          strm << (first ? "" : ",") << boost::io::quoted(val);
        else
          strm << (first ? "" : ",") << val;
        first = false;
      }
      statement->bind( column, strm.str() );
#else
      std::vector<unsigned char> bindata( totallen );
      
      size_t pos = 0;
      for( const string &str : v )
      {
        const size_t len = str.size();
        memcpy( &(bindata[pos]), &len, sizeof(size_t) );
        pos += sizeof(size_t);
        if( len )
          memcpy( &(bindata[pos]), &(str[0]), len );
        pos += len;
      }
      
      statement->bind( column, bindata );
#endif
    }
    bool sql_value_traits<std::set<std::string>>::read(std::set<std::string> &v, SqlStatement *statement, int column, int size)
    {
      v.clear();
#if( USE_TEXT_FOR_BLOB )
      string line;
      bool result = statement->getResult( column, &line, size );
      
      if( !result )
        return result;
      
      //ToDo: I havent checked if boost::tokenizer will throw exceptiuon or not!
      //      I dont believe it will, but I'll wrap in a try/catch for the moment.
      try
      {
        typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
        boost::escaped_list_separator<char> separator("\\",",", "\"");
        Tokeniser t( line, separator );
        for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
        {
          if( it->find_last_of("\",") != string::npos )
          {
            string unqotval;
            stringstream unqotstrm(*it);
            unqotstrm >> boost::io::quoted(unqotval); //Also unsure if boost::io::quoted could throw
            v.emplace( std::move(unqotval) );
          }else
          {
            v.insert(*it);
          }
      }
      }catch(...)
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Boost tokenizer does throw exception!" );
#endif
      }
      return result;
#else
      std::vector<unsigned char> bindata;
      
      bool result = statement->getResult( column, &bindata, size );
      
      if( !result )
        return result;
      
      if( bindata.empty() )
        return result;
      
      size_t position = 0;
      const size_t nbytes = bindata.size();
      while( position < nbytes )
      {
        if( (position + sizeof(size_t)) > nbytes )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "Database re-read in of a set<string> did not have proper size!  Programming logic error." );
#endif
          v.clear();
          return false;
        }//if( something unexpected is up )
        
        size_t nstrbytes;
        memcpy( &nstrbytes, &(bindata[position]), sizeof(size_t) );
        position += sizeof(size_t);
        
        if( (position + nstrbytes) > nbytes )
        {
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, "Database re-read in of a set<string> did not have proper size to read a str!  Programming logic error." );
#endif
          v.clear();
          return false;
        }//
        
        string str;
        str.resize(nstrbytes);
        if( nstrbytes )
          memcpy( &(str[0]), &(bindata[position]), sizeof(nstrbytes) );
        position += nstrbytes;
        v.emplace_hint(v.end(), std::move(str) );
      }//while( position < nbytes )
      
      if( position != nbytes )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Database re-read in of a set<string> had left over bytes!  Programming logic error." );
#endif
      }
      
      return result;
#endif
    }//read( std::set<std::string> )
  
    
    std::string sql_value_traits<SpecUtils::DetectorAnalysis>::type(SqlConnection *conn, int size)
    {
      return conn->textType(-1);
    }
    void sql_value_traits<SpecUtils::DetectorAnalysis>::bind( const SpecUtils::DetectorAnalysis &v, SqlStatement *statement, int column, int size)
    {
      rapidxml::xml_document<char> xmldoc;
      std::mutex xmldocmutex;
      add_analysis_results_to_2012_N42(v, &xmldoc, xmldocmutex );
        
      string xmlstring;
      rapidxml::print(std::back_inserter(xmlstring), xmldoc, 0);
        
      statement->bind(column, xmlstring );
    }
    bool sql_value_traits<SpecUtils::DetectorAnalysis>::read( SpecUtils::DetectorAnalysis &v, SqlStatement *statement, int column, int size )
    {
      v.reset();
      string xmlstring;
      bool result = statement->getResult( column, &xmlstring, size );
      if( !result || xmlstring.size()<20 )
        return result;
        
      try
      {
        rapidxml::xml_document<char> doc;
        doc.parse<rapidxml::parse_default>( &(xmlstring[0]) );
          
        set_analysis_info_from_n42( doc.first_node(), v );
      }catch(...)
      {
        cerr << "Failed to parse DetectorAnalysis XML: " << xmlstring << endl;
        return false;
      }
        
      return result;
    }
    
  
    std::string sql_value_traits<std::map<std::string,std::vector<std::string>>>::type(SqlConnection *conn, int size)
    {
#if( USE_TEXT_FOR_BLOB )
      return conn->textType(-1);
#else
      static_assert( 0, "You should specialize sql_value_traits::type for map<string,vector<string>> to store as binary data (or just use the slower XML version)" );
#endif
    }
    void sql_value_traits< std::map<std::string,std::vector<std::string>> >::bind(
                           const std::map<std::string,std::vector<std::string>> &v,
                           SqlStatement *statement, int column, int size)
    {
#if( USE_TEXT_FOR_BLOB )
/*  //For some reason I couldnt figure out in 10 minutes, the below constructing of an XML document using pugixml was not working
 // The "value" nodes werent getting created, and instead their values was being put into the "filter" nodes - weird!
 // Also, maybe the text values were garbled.
      pugi::xml_document xmldoc;
      
      for( const auto &np : v )
      {
        if( np.second.empty() )
          continue;
        
        pugi::xml_node filter_node = xmldoc.append_child("filter");
        pugi::xml_attribute attrib = filter_node.append_attribute("name");
        attrib.set_value( np.first.c_str() );
        
        for( const string &val : np.second )
        {
          pugi::xml_node child = filter_node.append_child(pugi::node_pcdata);
          child.set_name( "value" );
          child.set_value( val.c_str() );
        }
      }//for( const auto &np : v )
      
      xml_string_writer output;
      xmldoc.save( output );
 
      cout << "event_xml_filter_values->\"" << output.result << "\"" << endl << endl;
 
      statement->bind(column, output.result );
*/
      //Lets run with scisors: we will not allocate any rapidxml string, but
      //  instead rely on the values in 'v' - this should be fine, but also also
      //  using this level of indirection (map->pair->{string->c_str(),vector->string->c_str())
      //  seems a little precarious
      rapidxml::xml_document<char> xmldoc;
      
      for( const auto &np : v )
      {
        if( np.first.empty() || np.second.empty() )
          continue;
        
        auto filter_node = xmldoc.allocate_node( rapidxml::node_element, "filter", nullptr, 6, 0 );
        filter_node->append_attribute( xmldoc.allocate_attribute("name", np.first.c_str(),4,np.first.size()) );
        xmldoc.append_node( filter_node );
        
        for( const string &val : np.second )
        {
          auto child = xmldoc.allocate_node( rapidxml::node_element, "value", val.c_str(), 5, val.size() );
          filter_node->append_node( child );
        }
      }//for( const auto &np : v )
      
      string xml_data;
      rapidxml::print( std::back_inserter(xml_data), xmldoc, 0 );
      
      statement->bind(column, xml_data );
#else
      static_assert( 0, "You should specialize sql_value_traits::bind for map<string,vector<string>> to store as binary data (or just use the slower XML version)" );
#endif
    }
    bool sql_value_traits< std::map<std::string,std::vector<std::string>> >::read(
                           std::map<std::string,std::vector<std::string>> &v,
                           SqlStatement *statement, int column, int size )
    {
      v.clear();
#if( USE_TEXT_FOR_BLOB )
      string xmlstring;
      bool result = statement->getResult( column, &xmlstring, size );
      if( !result || xmlstring.size() < 20 )
        return result;
      
      pugi::xml_document doc;
      pugi::xml_parse_result xml_result = doc.load_string( xmlstring.c_str() );
      
      if( !xml_result )
      {
        //Shouldnt ever happen
        stringstream msg;
        msg << "Failed to parse event_xml_filter_values: \n"
        "Error description: " << xml_result.description() << "\n"
        "Error offset: " << xml_result.offset << ", size(xml)=" << xmlstring.size();
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, msg.str().c_str() );
#else
        cerr << msg.str().c_str() << endl;
#endif
        return false;
      }//if( failed to parse XML )
      
      for( const auto &filter_node : doc.children("filter") )
      {
        const auto name = filter_node.attribute("name").as_string();
        vector<string> values;
        for( const auto &value_node : filter_node.children("value") )
          values.push_back( value_node.text().as_string() );
        if( values.size() )
          v[name] = std::move(values);
      }//
      
      return result;
#else
      static_assert( 0, "You should specialize sql_value_traits::read for map<string,vector<string>> to store as binary data (or just use the slower XML version)" );
#endif
    }
  }//namespace Dbo
}//namespace Wt


std::vector<EventXmlFilterInfo> EventXmlFilterInfo::parseJsonString( const std::string &sourcestr )
{
  vector<EventXmlFilterInfo> eventinfos;
  
  try
  {
    Json::Value q;
    Wt::Json::parse( sourcestr, q );
    
    const Json::Array &basearray = q;
    
    
    int arraypos = 0;
    string base_node_test;
    for( const Json::Value &val : basearray )
    {
      const Json::Object &valobj = val;
      ++arraypos;
      if( arraypos == 1 )
      {
        //The very first entry in the array is different than all the rest; it has
        // four properties "author", "comments", "last_modified", "base_node_test"
        // We only care about "base_node_test".
        const auto &node = valobj.get("base_node_test");
        if( !node.isNull() && node.type()==Json::StringType )
          base_node_test = static_cast<const string &>(node);
        continue;
      }//if( arraypos == 0 )
      
      
      EventXmlFilterInfo info;
      info.m_base_node_test = base_node_test;
      
      const auto &label = valobj.get("label");
      if( label.isNull() || label.type()!=Json::StringType )
        throw runtime_error( "Each Event XML JSON object must have a string label" );
      
      info.m_label = Wt::Utils::htmlEncode( static_cast<const string &>(label) );
      
      //Make sure no duplicate labels
      bool alreadyExists = false;
      for( size_t i = 0 ; !alreadyExists && i < eventinfos.size(); ++i )
        alreadyExists = (eventinfos[i].m_label==info.m_label);
      if( alreadyExists )
        throw runtime_error( "Event XML JSON has objects with duplicate label values ('" + info.m_label + "')" );
      
      //Make sure the JSON wont stomp over a spectrum file query field label.
      for( SpecFileQuery::FileDataField f = SpecFileQuery::FileDataField(0);
          f != SpecFileQuery::FileDataField::NumFileDataFields;
          f = SpecFileQuery::FileDataField(f+1) )
        if( info.m_label == to_string(f) )
          throw runtime_error( "Event XML JSON has objects un-allowed label value ('" + info.m_label + "')" );
      
      //"show_result_in_gui": true,
      const auto &showInGui = valobj.get("show_result_in_gui");
      if( showInGui.isNull() || showInGui.type()!=Json::BoolType )
        info.show_in_gui_table = true;
      else
        info.show_in_gui_table = static_cast<bool>( showInGui );
        
        
      const size_t len = SpecUtils::utf8_str_len(info.m_label.c_str(), info.m_label.size());
      if( len==0 || len > 32 )
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' has an invalid label; must be between 1 and 32 characters." );
      
      string typestr;
      const auto &type = valobj.get("type");
      if( !type.isNull() && type.type()==Json::StringType )
        typestr = static_cast<const string &>(type);
      
      if( SpecUtils::iequals_ascii(typestr, "select") )
        info.m_type = EventXmlFilterInfo::InputType::Select;
      else if( SpecUtils::iequals_ascii(typestr, "text") )
        info.m_type = EventXmlFilterInfo::InputType::Text;
      else if( SpecUtils::iequals_ascii(typestr, "date") )
        info.m_type = EventXmlFilterInfo::InputType::Date;
      else
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' does not have a valid type property (can have string values 'select' or 'text')." );
      
      const auto &xpath = valobj.get("xpath");
      if( xpath.isNull() || xpath.type()!=Json::StringType || (static_cast<const string &>(xpath).empty()) )
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' does not have a string xpath property specified." );
      info.m_xpath = static_cast<const string &>(xpath);
      
      try
      {
        pugi::xpath_query q(info.m_xpath.c_str(), nullptr); //should throw if invalid xpath
      }catch( std::exception & )
      {
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' has an invalid xpath given." );
      }
      
      const auto &placeholder = valobj.get("placeholder");
      if( !placeholder.isNull() && placeholder.type()==Json::StringType )
      {
        info.m_placeholder = Wt::Utils::htmlEncode( static_cast<const string &>(placeholder) );
        const size_t len = SpecUtils::utf8_str_len(info.m_placeholder.c_str(), info.m_placeholder.size());
        if( len > 64 )
          throw runtime_error( "The Event XML query condition '" + info.m_label
                              + "' has aa placeholder text larger than 64 characters (not allowed)." );
      }
      
      const auto &operators = valobj.get("operators");
      if( operators.isNull() || operators.type()!=Json::ArrayType )
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' does not have the property 'operators' defined as an array of strings." );
      
      const string str_allowed_ops[] = { "contains", "equal", "not equal",
        "does not contain", "begins with", "does not begin with", "ends with",
        "does not end with", "regex"
      };
      const string discrete_allowed_ops[] = { "equal", "not equal" };
      const string date_allowed_ops[] = { "greater","less","equal","not_equal" };
      
      for( const Json::Value &op : static_cast<const Json::Array &>(operators) )
      {
        if( op.type() != Json::StringType )
          throw runtime_error( "The Event XML query condition '" + info.m_label
                              + "' has an array entry to 'operators' that is not a string." );
        string strval = static_cast<const string &>(op);
        SpecUtils::trim( strval );
        SpecUtils::to_lower_ascii( strval );
        
        switch( info.m_type )
        {
          case EventXmlFilterInfo::InputType::Text:
            if( std::find(begin(str_allowed_ops), end(str_allowed_ops), strval) == end(str_allowed_ops) )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has an invalid 'operators' entry for a string field ('" + strval + "')." );
            break;
          case EventXmlFilterInfo::InputType::Select:
            if( std::find(begin(discrete_allowed_ops), end(discrete_allowed_ops), strval) == end(discrete_allowed_ops) )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has an invalid 'operators' entry for a select ('" + strval + "')." );
            
          case EventXmlFilterInfo::InputType::Date:
            if( std::find(begin(date_allowed_ops), end(date_allowed_ops), strval) == end(date_allowed_ops) )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has an invalid 'operators' entry for a date ('" + strval + "')." );
        }//switch( info.m_type )
        
        info.m_operators.push_back( strval );
      }//for( const Json::Value &op : static_cast<const Json::Array &>(operators) )
      
      if( info.m_operators.empty() )
        throw runtime_error( "The Event XML query condition '" + info.m_label
                            + "' does not have at least one string entry in its 'operators' property." );
      
      
      switch( info.m_type )
      {
        case EventXmlFilterInfo::InputType::Text:
          break;
          
        case EventXmlFilterInfo::InputType::Select:
        {
          const auto &discreetops = valobj.get("discreet_options");
          if( discreetops.isNull() || discreetops.type()!=Json::ArrayType )
            throw runtime_error( "The Event XML query condition '" + info.m_label
                                + "' does not have the property 'discreet_options' defined as an array of strings "
                                "(required since a type=='select' specified)." );
          
          for( const Json::Value &op : static_cast<const Json::Array &>(discreetops) )
          {
            if( op.type() != Json::StringType )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has an array entry to 'discreet_options' that is not a string." );
            string strval = Wt::Utils::htmlEncode( static_cast<const string &>(op) );
            SpecUtils::trim( strval );
            const size_t len = SpecUtils::utf8_str_len(strval.c_str(), strval.size());
            if( len==0 )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has an empty discreet option (not allowed)." );
            
            if( len > 64 )
              throw runtime_error( "The Event XML query condition '" + info.m_label
                                  + "' has a discreet option longer than the allowed 64 characters." );
            
            info.m_discreet_options.push_back( strval );
          }//for( const Json::Value &op : static_cast<const Json::Array &>(discreetops) )
          
          if( info.m_discreet_options.size() < 2 )
            throw runtime_error( "The Event XML query condition '" + info.m_label
                                + "' does not have at least two discreet option as required." );
          
          break;
        }//case EventXmlFilterInfo::InputType::Select:
          
        case EventXmlFilterInfo::InputType::Date:
          break;
      }//switch( info.m_type )
      
      eventinfos.push_back( info );
    }//for( const auto Json::Value &val : basearray )
  }catch( Wt::Json::ParseError &e )
  {
    cerr << "Error: " << e.what() << endl << endl << endl;
    throw runtime_error( "The Event XML JSON defininting conditions is not valid JSON - not adding in Event XML search terms." );
  }catch( std::runtime_error & )
  {
    throw;
  }//try /catch
  
  return eventinfos;
}//vector<EventXmlFilterInfo> parseJsonString(const string &)




SpecFileInfoToQuery::SpecFileInfoToQuery()
{
  reset();
}

void SpecFileInfoToQuery::reset()
{
  file_path.clear();
  file_size = file_path_hash = 0;
  is_file = is_spectrum_file = is_event_xml_file = false;
  
  filename.clear();
  detector_names.clear();
  serial_number.clear();
  manufacturer.clear();
  model.clear();
  uuid.clear();
  file_remarks.clear();
  record_remarks.clear();
  location_name.clear();
  
  has_riid_analysis = false;
  riid_ana = SpecUtils::DetectorAnalysis();
  detector_type = SpecUtils::DetectorType::Unknown;
  passthrough = false;
  total_livetime = total_realtime = 0.0f;
  contained_neutron = contained_dev_pairs = contained_gps = false;
  
  energy_cal_types.clear();
  individual_spectrum_live_time.clear();
  individual_spectrum_real_time.clear();
  number_of_samples = number_of_records = 0;
  number_of_gamma_channels.clear();
  max_gamma_energy.clear();
  mean_latitude = mean_longitude = -999.9f;
  
  neutron_count_rate.clear();
  gamma_count_rate.clear();
  start_times.clear();
  
  event_xml_filter_values.clear();
}//void SpecFileInfoToQuery::reset()


void SpecFileInfoToQuery::fill_info_from_file( const std::string filepath )
{
  reset();
  
  file_path = filepath;
  is_file = SpecUtils::is_file(filepath);
  if( !is_file )
    return;
  
  file_size = SpecUtils::file_size(filepath);
  file_path_hash = std::hash<std::string>()(filepath);
  
  SpecUtils::SpecFile meas;
  const bool loaded = meas.load_file(filepath, SpecUtils::ParserType::Auto, filepath);
  if( !loaded )
  {
    //see if "<event>" is in the first 128 bytes.
#ifdef _WIN32
    const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16( filepath );
    ifstream datstrm( wfilepath.c_str(), ios_base::binary|ios_base::in );
#else
    ifstream datstrm( filepath.c_str(), ios_base::binary|ios_base::in );
#endif
    
    if( datstrm.is_open() )
    {
      string startdata( 129, '\0' );
      if( datstrm.read(&startdata[0], 128) )
        is_event_xml_file = (startdata.find( "<event>" ) != string::npos);
    }
    return;
  }//if( not a spectrum file )
  
  is_spectrum_file = true;
  
  meas.set_filename( filepath );
  
  filename = meas.filename();
  
  detector_names.insert( meas.detector_names().begin(), meas.detector_names().end() );
  serial_number = meas.instrument_id();
  manufacturer = meas.manufacturer();
  model = meas.instrument_model();
  uuid = meas.uuid();
  file_remarks.insert( meas.remarks().begin(), meas.remarks().end() );
  location_name = meas.measurement_location_name();
  
  has_riid_analysis = !!meas.detectors_analysis();
  if( meas.detectors_analysis() )
    riid_ana = *meas.detectors_analysis();
  
  detector_type = meas.detector_type();
  passthrough = meas.passthrough();
  
  number_of_samples = meas.sample_numbers().size();
  number_of_records = meas.num_measurements();
  
  contained_gps = meas.has_gps_info();
  if( contained_gps )
  {
    mean_latitude = meas.mean_latitude();
    mean_longitude = meas.mean_longitude();
  }
  
  total_livetime = total_realtime = 0.0f;
  contained_neutron = contained_dev_pairs = false;
  
  for( const auto m : meas.measurements() )
  {
    total_livetime += m->live_time();
    total_realtime += m->real_time();
    
    record_remarks.insert( m->remarks().begin(), m->remarks().end() );
    if( m->num_gamma_channels() > 6 )//skip over GMTubes and other gross-count gamma detectors
    {
      energy_cal_types.insert( m->energy_calibration_model() );
      
      if( m->live_time() > 0.0f )
        individual_spectrum_live_time.insert( m->live_time() );
      if( m->real_time() > 0.0f )
        individual_spectrum_real_time.insert( m->real_time() );
      number_of_gamma_channels.insert( m->num_gamma_channels() );
      const float lt = ((m->live_time() > 0.001f) ? m->live_time() : m->real_time());
      if( lt > 0.001f )
        gamma_count_rate.insert( m->gamma_count_sum() / lt );
      max_gamma_energy.insert( m->gamma_energy_max() );
    }
    
    if( m->contained_neutron() )
    {
      contained_neutron = true;
      const float rt = ((m->real_time() > 0.001f) ? m->real_time() : m->live_time());
      if( rt > 0.001f )
        neutron_count_rate.insert( m->neutron_counts_sum() / rt );
    }
    
    contained_dev_pairs = (contained_dev_pairs || !m->deviation_pairs().empty());
    
    if( !m->start_time().is_special() )
    {
      const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
      boost::posix_time::time_duration::sec_type x = (m->start_time() - epoch).total_seconds();
      start_times.insert( static_cast<time_t>(x) );
    }
  }
}//void fill_info_from_file( const std::string filepath )


void SpecFileInfoToQuery::fill_event_xml_filter_values( const std::string &filepath,
                                                       const std::vector<EventXmlFilterInfo> &xmlfilters )
{
  event_xml_filter_values.clear();
  
  if( xmlfilters.empty() )
    return;
  
  const auto path = SpecUtils::parent_path( filepath );
  auto candidates = SpecUtils::ls_files_in_directory( path, &xml_files_small_enough, nullptr );
  
  //Assume all xmlfilters have the same filter.
  //   ToDo: add development check for this
  const string &base_node_test = xmlfilters[0].m_base_node_test;
  
  for( const auto &xmlfilename : candidates )
  {
#ifdef _WIN32
    const std::wstring wxmlfilename = SpecUtils::convert_from_utf8_to_utf16( xmlfilename );
    ifstream f( wxmlfilename.c_str() );
#else
    ifstream f( xmlfilename.c_str() );
#endif
    
    if( !f )
      continue;
    
    if( !base_node_test.empty() )
    {
      string start( 128, ' ');
      f.read( &(start[0]), 127 );
      if( !SpecUtils::icontains(start, base_node_test) )
        continue;
    }//if( !base_node_test.empty() )
    
    pugi::xml_document doc;
#ifdef _WIN32
    if( !doc.load_file( wxmlfilename.c_str() ) )
      continue;
#else
    if( !doc.load_file( xmlfilename.c_str() ) )
      continue;
#endif
    
    for( const EventXmlFilterInfo &test : xmlfilters )
    {
      try
      {
        pugi::xpath_node_set tools = doc.select_nodes( test.m_xpath.c_str() );
        
        for( pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it )
        {
          const string strval = it->node().text().as_string();
          if( !strval.empty() ) //probably shouldnt ever be empty, but could if the xpath is wonky
            event_xml_filter_values[test.m_label].push_back( strval );
        }//for( pugi::xpath_node_set::const_iterator it = tools.begin(); it != tools.end(); ++it )
      }
#if( PERFORM_DEVELOPER_CHECKS )
      catch( std::exception &e ){
        log_developer_error( __func__, ("Unexpected exception performing xpath query: " + string(e.what())).c_str() );
      }
#else
      catch( std::exception & ){
        //We should have already checked the xpath - so we really shouldnt ever get here.
      }
#endif
    }//for( const auto &xmlfilename : candidates )
    
    //In principle if we made it here we found a Event XML file - so lets not
    //  keep wasting time looking in other files (on the other hand, would it be
    //  reasonable to allow using more than one Event XML file?)
    break;
  }//for( const EventXmlFilterInfo &test : xmlfilters )
  
}//void fill_event_xml_filter_values( const std::string filepath )



SpecFileQueryDbCache::SpecFileQueryDbCache( const bool use_db_caching,
                                            const std::string &path,
                                            const std::vector<EventXmlFilterInfo> &xmlfilters )
  : m_use_db_caching( use_db_caching ),
    m_using_persist_caching( false ),
    m_fs_path( path ),
    m_xmlfilters( xmlfilters )
{
  m_stop_caching = false;
  m_doing_caching = false;
  
  if( m_use_db_caching )
    m_using_persist_caching = init_existing_persisted_db();
  
  if( m_use_db_caching && !m_using_persist_caching )
  {
    string db_location = SpecUtils::temp_file_name( "interspec_file_query", SpecUtils::temp_dir() );
    m_use_db_caching = open_db( db_location, true );
  }//if( m_use_db_caching )
}//SpecFileQueryDbCache( constructor )


SpecFileQueryDbCache::~SpecFileQueryDbCache()
{
  if( m_use_db_caching )
  {
    std::unique_lock<std::mutex> lock( m_cv_mutex );
    m_stop_caching = true;
    
    if( m_doing_caching )
      m_cv.wait( lock, [this]() -> bool { return !m_doing_caching;} );
  }//if( m_use_db_caching )
  
  if( m_use_db_caching )
  {
    std::lock_guard<std::mutex> lock( m_db_mutex );
    m_db_session.reset();
    m_db.reset();
    
    if( !m_using_persist_caching )
      SpecUtils::remove_file( m_db_location );
  }//if( m_use_db_caching )
}//~SpecFileQueryDbCache()


std::string SpecFileQueryDbCache::construct_persisted_db_filename( std::string persisted_path )
{
  SpecUtils::make_canonical_path(persisted_path);
  
  const auto path_hash = std::hash<std::string>()( persisted_path );
  return SpecUtils::append_path( persisted_path, "InterSpec_file_query_cache_" + std::to_string(path_hash) + ".sqlite3" );
}//std::string construct_persisted_db_filename( std::string basepath )


bool SpecFileQueryDbCache::open_db( const std::string &path, const bool create_tables )
{
  try
  {
    std::unique_ptr<Wt::Dbo::backend::Sqlite3> db( new Wt::Dbo::backend::Sqlite3(path) );
    std::unique_ptr<Wt::Dbo::Session> db_session( new Wt::Dbo::Session() );
    
    db->setProperty( "show-queries", "false" );
    db_session->setConnection( *db );
    
    db_session->mapClass<SpecFileInfoToQuery>( "SpecFileInfoToQuery" );
    if( create_tables )
      db_session->createTables();
    
    m_use_db_caching = true;
    m_db_location = path;
    m_db_session = std::move( db_session );
    m_db = std::move( db );
  }catch( Wt::Dbo::Exception &e )
  {
    cerr << "Failed to create SpecFileQueryDbCache session, Dbo::Exception: " << e.what() << endl;
    cerr << endl;
    return false;
  }catch( std::exception &e )
  {
    cerr << "Failed to create SpecFileQueryDbCache session, std::exception: " << e.what() << endl;
    cerr << endl;
    return false;
  }
  return true;
}//bool open_db( const std::string &path, const bool map_classes );


bool SpecFileQueryDbCache::init_existing_persisted_db()
{
  if( !SpecUtils::is_directory( m_fs_path ) )
  {
    std::cerr << "m_fs_path='" << m_fs_path << "', is not a path" << std::endl;
    return false;
  }

  // The read-only atribute of a Windows directory doesnt mean you cant read/write in that dir
  //  so we'll just always skip checking if you can RW in a directory, even on POSIX, and 
  //  instead wait for openeing of the database to fail.
  //if( !SpecUtils::can_rw_in_directory( m_fs_path ) )
  //{
  //  std::cerr << "Can not RW in directory '" << m_fs_path << "'." << std::endl;
  //  return false;
  //}

  const string persisted_path = construct_persisted_db_filename(m_fs_path);
  
  if( !SpecUtils::is_file( persisted_path ) )
  {
    std::cerr << "Persisted db path '" << persisted_path << "' is not a file" << std::endl; 
    return false;
  }

  try
  {
    std::unique_ptr<Wt::Dbo::backend::Sqlite3> db( new Wt::Dbo::backend::Sqlite3(persisted_path) );
    std::unique_ptr<Wt::Dbo::Session> db_session( new Wt::Dbo::Session() );
    
    db->setProperty( "show-queries", "false" );
    db_session->setConnection( *db );
    db_session->mapClass<SpecFileInfoToQuery>( "SpecFileInfoToQuery" );
    
    
    //If there are any entries in the database, check the schema is okay by
    //  grabbing a result.  If no entries, try to insert something and then
    //  remove it.
    Wt::Dbo::Transaction trans( *db_session );
    auto query = db_session->find<SpecFileInfoToQuery>();
    auto results = query.resultList();
    
    //Note, not calling results.size() because it can take like 10 seconds for
    //  a large database, much faster just to iterate over results.
    bool gotentry = false;
    for( Dbo::collection<Dbo::ptr<SpecFileInfoToQuery>>::const_iterator iter = results.begin();
         !gotentry && (iter != results.end()); ++iter )
    {
      gotentry = true;
      
      //Force a read from the database (I think reading from DB is lazy, but not sure)
      cout << "First entry is for '" << (*iter)->file_path << "'" << endl;
    }
    trans.commit();
    
    if( !gotentry )
    {
      //This happens when the cache file in the temp location 
      // doesnt have any entries yet because the recursive_ls is
      //  still running.
      
      /*
       if( SpecUtils::file_size( persisted_path ) < 16 * 1024 )
       {
         //Delete the file and create a new one.
       } else
       {
       //Something may be up - pass a message to the user?
         throw runtime_error( "No results in database" );
       }
       */
      
      Dbo::ptr<SpecFileInfoToQuery> dummy;
      {
        Wt::Dbo::Transaction innertrans( *db_session );
        SpecFileInfoToQuery *dummyptr = new SpecFileInfoToQuery();
        dummyptr->file_path = "Dummy";
        dummy = db_session->add( dummyptr );
        innertrans.commit();
      }
      
      {
        Wt::Dbo::Transaction innertrans( *db_session );
        dummy.remove();
        innertrans.commit();
      }
    }//if( !gotentry )
  
    m_db = std::move( db );
    m_db_session = std::move( db_session );
    m_db_location = persisted_path;
    m_using_persist_caching = true;
    return true;
  }catch( std::exception &e )
  {
    m_db_session.reset();
    m_db.reset();
    //passMessage( "Caught exception initind persisted DB; " + string(e.what()), "", 3 );
    cerr << "\n\nIn persisting cache, caught: " << e.what() << endl << endl;
  }
  return false;
}//bool init_existing_persisted_db()


#if( PERFORM_DEVELOPER_CHECKS )
bool operator==( const SpecFileInfoToQuery &lhs, const SpecFileInfoToQuery &rhs )
{
  if( rhs.file_path != lhs.file_path )
    return false;
  if( rhs.file_size != lhs.file_size )
    return false;
  if( rhs.file_path_hash != lhs.file_path_hash )
    return false;
  if( rhs.is_file != lhs.is_file )
    return false;
  if( rhs.is_spectrum_file != lhs.is_spectrum_file )
    return false;
  if( rhs.is_event_xml_file != lhs.is_event_xml_file )
    return false;
  if( rhs.filename != lhs.filename )
    return false;
  if( rhs.detector_names != lhs.detector_names )
    return false;
  if( rhs.serial_number != lhs.serial_number )
    return false;
  if( rhs.manufacturer != lhs.manufacturer )
    return false;
  if( rhs.model != lhs.model )
    return false;
  if( rhs.uuid != lhs.uuid )
    return false;
  if( rhs.file_remarks != lhs.file_remarks )
    return false;
  if( rhs.record_remarks != lhs.record_remarks )
    return false;
  if( rhs.location_name != lhs.location_name )
    return false;
  if( rhs.has_riid_analysis != lhs.has_riid_analysis )
    return false;
  if( rhs.has_riid_analysis )
  {
    if( rhs.riid_ana.remarks_ != rhs.riid_ana.remarks_ )
      return false;
    if( rhs.riid_ana.algorithm_name_ != rhs.riid_ana.algorithm_name_ )
      return false;
    if( rhs.riid_ana.algorithm_component_versions_ != rhs.riid_ana.algorithm_component_versions_ )
      return false;
    if( rhs.riid_ana.algorithm_creator_ != rhs.riid_ana.algorithm_creator_ )
      return false;
    if( rhs.riid_ana.algorithm_description_ != rhs.riid_ana.algorithm_description_ )
      return false;
    if( rhs.riid_ana.algorithm_result_description_ != rhs.riid_ana.algorithm_result_description_ )
      return false;
    
    //rhs.riid_ana.results_ ==
  }
  
  if( rhs.detector_type != lhs.detector_type )
    return false;
  if( rhs.passthrough != lhs.passthrough )
    return false;
  if( rhs.total_livetime != lhs.total_livetime )
    return false;
  if( rhs.total_realtime != lhs.total_realtime )
    return false;
  if( rhs.contained_neutron != lhs.contained_neutron )
    return false;
  if( rhs.contained_dev_pairs != lhs.contained_dev_pairs )
    return false;
  if( rhs.contained_gps != lhs.contained_gps )
    return false;
  if( rhs.energy_cal_types != lhs.energy_cal_types )
    return false;
  if( rhs.individual_spectrum_live_time != lhs.individual_spectrum_live_time )
    return false;
  if( rhs.individual_spectrum_real_time != lhs.individual_spectrum_real_time )
    return false;
  if( rhs.number_of_samples != lhs.number_of_samples )
    return false;
  if( rhs.number_of_records != lhs.number_of_records )
    return false;
  if( rhs.number_of_gamma_channels != lhs.number_of_gamma_channels )
    return false;
  if( rhs.max_gamma_energy != lhs.max_gamma_energy )
    return false;
  if( rhs.mean_latitude != lhs.mean_latitude )
    return false;
  if( rhs.mean_longitude != lhs.mean_longitude )
    return false;
  if( rhs.neutron_count_rate != lhs.neutron_count_rate )
    return false;
  if( rhs.gamma_count_rate != lhs.gamma_count_rate )
    return false;
  if( rhs.start_times != lhs.start_times )
    return false;
  return true;
}
#endif

void SpecFileQueryDbCache::cache_results( const std::vector<std::string> &&files )
{
  if( !m_use_db_caching )
    return;
  
  {//begin scope lock on m_cv_mutex
    std::unique_lock<std::mutex> lock( m_cv_mutex );
    if( m_stop_caching )
      return;
    
    if( m_doing_caching )
    {
      //Shouldnt ever happen
      cerr << "SpecFileQueryDbCache::cache_results can not be called multiple times" << endl;
      return;
    }
    
    {
      std::lock_guard<std::mutex> lock( m_db_mutex ); //always locked inside of m_cv_mutex
      m_doing_caching = (!!m_db && !!m_db_session);
    }
    
    if( !m_doing_caching )
      return;
  }//end scope lock on m_cv_mutex
  
  
  for( const string filename : files )
  {
    {//begin lock on m_cv_mutex
      std::unique_lock<std::mutex> lock( m_cv_mutex );
      if( m_stop_caching )
      {
        m_doing_caching = false;
        m_cv.notify_all();
        return;
      }
    }//end lock on m_cv_mutex
    
    try
    {
      {//begin lock on m_db_mutex
        std::lock_guard<std::mutex> lock( m_db_mutex );
        
        const long long filenamehash = static_cast<long long>( std::hash<std::string>()(filename) );
        const size_t filesize = SpecUtils::file_size(filename);
        
        Wt::Dbo::Transaction trans( *m_db_session );
        auto results = m_db_session->find<SpecFileInfoToQuery>().where( "file_path_hash = ?" ).bind(filenamehash).resultList();
        
        bool have_in_db = false;
        if( results.size() == 1 )
        {
          have_in_db = (results.front()->file_size == static_cast<long long>(filesize));
          if( !have_in_db )
            results.front().remove();
        }if( results.size() > 1 )
        {
          for( auto p : results )
            p.remove();
        }
        
        if( have_in_db )
          continue;
        trans.commit();
      }//end lock on m_db_mutex
      
      auto dbinforaw = new SpecFileInfoToQuery();;
      Wt::Dbo::ptr<SpecFileInfoToQuery> dbinfo( dbinforaw );
      dbinforaw->fill_info_from_file(filename);
      dbinforaw->fill_event_xml_filter_values(filename,m_xmlfilters);
      
      {//begin lock on m_db_mutex
        std::lock_guard<std::mutex> lock( m_db_mutex );
        Wt::Dbo::Transaction trans( *m_db_session );
        
        //We could be adding this file info uncessarily to the database, but I think end-logic will be fine...
        m_db_session->add( dbinfo );
        
        trans.commit();
      }//end lock on m_db_mutex
      
#if( PERFORM_DEVELOPER_CHECKS )
      {//Begin check we can read back in the identical object from the database
        std::lock_guard<std::mutex> lock( m_db_mutex );
        
        Wt::Dbo::Transaction trans( *m_db_session );
        auto results = m_db_session->find<SpecFileInfoToQuery>().where( "file_path_hash = ?" ).bind(dbinfo->file_path_hash).resultList();
        if( !results.size() )
        {
          log_developer_error( __func__, "Failed to find SpecFileInfoToQuery I just saved!!  Programming logic error." );
        }else
        {
          auto fromdb = results.front();
          if( !((*fromdb) == (*dbinfo)) )
          {
            log_developer_error( __func__, "The SpecFileInfoToQuery from database is not equal to the one saved to database!  Programming logic error." );
          }
        }
        trans.commit();
      }//End check we can read back in the identical object from the database
#endif
    }catch( Wt::Dbo::Exception &e )
    {
      //I think we get here mostly when an entry for a particular hash is already
      //  in the database.
      cerr << "Dbo::Exception caching spec files to database - oh well: " << e.what() << endl;
    }catch( std::exception &e )
    {
      cerr << "std::exception caching spec files to database - oh well: " << e.what() << endl;
    }
  }//for( const string filename : files )
  
  
  {//begin lock on m_cv_mutex
    std::lock_guard<std::mutex> lock( m_cv_mutex );
    m_doing_caching = false;
    m_cv.notify_all();
  }//end lock on m_cv_mutex
}//void cache_results()


void SpecFileQueryDbCache::stop_caching()
{
  if( m_use_db_caching )
  {
    std::unique_lock<std::mutex> lock( m_cv_mutex );
    m_stop_caching = true;
  
    if( m_doing_caching )
      m_cv.wait( lock, [this]() -> bool { return !m_doing_caching;} );
  }//if( m_use_db_caching )
}//void stop_caching();


void SpecFileQueryDbCache::allow_start_caching()
{
  if( m_use_db_caching )
  {
    std::unique_lock<std::mutex> lock( m_cv_mutex );
    m_stop_caching = false;
  }
}


bool SpecFileQueryDbCache::caching_enabled()
{
  return m_use_db_caching;
}//bool caching_enabled()


bool SpecFileQueryDbCache::is_persistand() const
{
  return (m_use_db_caching && m_using_persist_caching);
}//bool is_persistand() const;


void SpecFileQueryDbCache::set_persist( const bool persist )
{
  if( !m_use_db_caching )
    throw runtime_error( "Could not persist caching because caching is not enabled." );
  
  if( persist == m_using_persist_caching )
    return;
 
  std::lock_guard<std::mutex> lock( m_db_mutex );
  
  m_using_persist_caching = persist;
  
  if( persist )
  {
    const string old_db = m_db_location;
    
    //If the directory already has a persisted store - use if (note: this wont
    //  be the case normally (ever?) as we would already be using it on class
    //  initialization)
    if( init_existing_persisted_db() )
    {
      if( old_db != m_db_location )
        SpecUtils::remove_file(old_db);
      return;
    }
    
    const string newfilename = construct_persisted_db_filename( m_fs_path );
    if( SpecUtils::iends_with( newfilename, ".sqlite3" ) //JIC
        && SpecUtils::is_file(newfilename) )
      SpecUtils::remove_file(newfilename);
    
    if( m_db && m_db_session )
    {
      if( newfilename == m_db_location )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Database new and old locations same - Programming logic error." );
#endif
        return;
      }
      
      m_db_session.reset();
      m_db.reset();

      if( SpecUtils::is_file(m_db_location) ) //This is being a little over-cautious, but whatever.
      {
        if( !SpecUtils::rename_file( m_db_location, newfilename ) )
        {
          m_db_location = "";
          m_use_db_caching = false;
          m_using_persist_caching = false;

          throw runtime_error( "Failed to move existing cache database from '" + m_db_location + "' to '" + newfilename + "'.  Disabling caching." );
        }//

        const string oldloc = m_db_location;
        m_db_location = newfilename;
       
        if( init_existing_persisted_db() )
            return;
       
        m_db_location = "";
        m_use_db_caching = false;
        m_using_persist_caching = false;

        throw runtime_error( "Failed to open cache database after moving from '" 
                             + oldloc + "' to '" + newfilename + "'.  Disabling caching." );
      }else
      {
        //We can actually fall through here, and a new data base will attempt to be made
      }//if( SpecUtils::is_file(m_db_location) )
    }//if( had a cache database )
    
    if( open_db( newfilename, true ) )
      return;
    
    m_db_location = "";
    m_use_db_caching = false;
    m_using_persist_caching = false;
    
    throw runtime_error( "Failed to open cache database in '" + m_fs_path + "'.  Disabling caching." );
  }//if( persist )

  
  //if we're here, we want to stop using a persisted database
  if( SpecUtils::is_file(m_db_location) )
  {
    bool create_tables = false;
    string newname = SpecUtils::temp_file_name( "interspec_file_query", SpecUtils::temp_dir() );
    
    if( !SpecUtils::rename_file(m_db_location, newname) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      char buffer[2048];
      snprintf( buffer, sizeof(buffer), "Failed to rename persisted file from '%s' to '%s'.",
               m_db_location.c_str(), newname.c_str() );
      log_developer_error( __func__, "Failed to rename persisted file from." );
#endif
      create_tables = true;
      newname = SpecUtils::temp_file_name( "interspec_file_query", SpecUtils::temp_dir() );
    }
    
    if( open_db(newname, create_tables) )
      return;
    
    m_db_location = "";
    m_use_db_caching = false;
    m_using_persist_caching = false;
    
    throw runtime_error( "Failed to open database after moving it to a temp folder - diabling caching" );
  }//if( is_file(m_db_location) )
  
  
  const string newname = SpecUtils::temp_file_name( "interspec_file_query", SpecUtils::temp_dir() );
  if( open_db(newname, true) )
    return;
  
  m_db_location = "";
  m_use_db_caching = false;
  m_using_persist_caching = false;
  
  throw runtime_error( "Failed to open database in a temp folder - diabling caching" );
}//void set_persist( const bool persist, const std::string basepath );


const std::vector<EventXmlFilterInfo> &SpecFileQueryDbCache::xml_filters() const
{
  return m_xmlfilters;
}


std::unique_ptr<SpecFileInfoToQuery> SpecFileQueryDbCache::spec_file_info( const std::string &filepath )
{
  //ToDo: A copy or two of the results could probably be eliminated in this function.
  std::unique_ptr<SpecFileInfoToQuery> info( new SpecFileInfoToQuery() );
  
  if( !m_use_db_caching )
  {
    info->fill_info_from_file( filepath );
    info->fill_event_xml_filter_values(filepath,m_xmlfilters);
    return info;
  }
  
  const long long filenamehash = static_cast<long long>( std::hash<std::string>()(filepath) );
  const size_t filesize = SpecUtils::file_size(filepath);
  
  try
  {
    {//begin check in DB
      std::lock_guard<std::mutex> lock( m_db_mutex );
      
      if( !m_db || !m_db_session )
      {
        //Shouldnt ever get here!
        info->fill_info_from_file( filepath );
        info->fill_event_xml_filter_values(filepath,m_xmlfilters);
        return info;
      }
      
      Wt::Dbo::Transaction trans( *m_db_session );
      auto results = m_db_session->find<SpecFileInfoToQuery>().where( "file_path_hash = ?" ).bind(filenamehash).resultList();
      if( results.size() == 1 )
      {
        Wt::Dbo::ptr<SpecFileInfoToQuery> c = results.front();
        if( c->file_size == filesize )
        {
          *info = *c;
          return info;
        }
      }else
      {
        for( auto p : results )
          p.remove();
      }
      trans.commit();
    }//end check in DB
    
    //Wasnt in the database
    info->fill_info_from_file( filepath );
    info->fill_event_xml_filter_values( filepath, m_xmlfilters );
    
    auto dbinforaw = new SpecFileInfoToQuery();
    Wt::Dbo::ptr<SpecFileInfoToQuery> dbinfo( dbinforaw );
    *dbinforaw = *info;
    
    {//begin check in DB
      std::lock_guard<std::mutex> lock( m_db_mutex );
      Wt::Dbo::Transaction trans( *m_db_session );
      
      m_db_session->add( dbinfo );
      trans.commit();
    }//end check in DB
  }catch( Wt::Dbo::Exception &e )
  {
    cerr << "Caught Dbo::Exception in spec_file_info: " << e.what() << endl;
    info->fill_info_from_file( filepath );
    info->fill_event_xml_filter_values( filepath, m_xmlfilters );
  }catch( std::exception &e )
  {
    cerr << "Caught std::Exception in spec_file_info: " << e.what() << endl;
    info->fill_info_from_file( filepath );
    info->fill_event_xml_filter_values( filepath, m_xmlfilters );
  }
  
  return info;
}//SpecFileInfoToQuery spec_file_info( const std::string &filepath )



