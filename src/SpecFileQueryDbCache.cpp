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

#include <numeric>

#include <Wt/Utils>
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/WStringStream>

#include <Wt/Dbo/Dbo>
#include <Wt/Dbo/WtSqlTraits>
#include <Wt/Dbo/backend/Sqlite3>

#include "pugixml.hpp"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/SpecFileQuery.h"
#include "InterSpec/FarmOptions.h"
#include "InterSpec/FarmAnalysis.h"
#include "InterSpec/EnrichmentResults.h"
#include "InterSpec/SpecFileQueryDbCache.h"
#include "SpecUtils/EnergyCalibration.h"

#include <nlohmann/json.hpp>

#include <boost/config.hpp>
#include <boost/io/quoted.hpp>


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
//#define USE_TEXT_FOR_BLOB 0
//static_assert( 0, "Serialization to database hasnt been checked (but is probably okay, assuming this version of Wt doesnt have the bug 3.3.4 does)" );
#define USE_TEXT_FOR_BLOB 1

#if( WT_VERSION >= 0x4000000)
#ifdef _MSC_VER
#pragma message("You should check if USE_TEXT_FOR_BLOB can be set to zero, for better serialization to database, since Wt is newer than 3.3.4")
#else
#warning "You should check if USE_TEXT_FOR_BLOB can be set to zero, for better serialization to database, since Wt is newer than 3.3.4"
#endif
#endif //#if( WT_VERSION >= 0x4000000)

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
  
  
  /** Function taken from full-spectrum 20231002 */
  bool potentially_analyze_derived_data( const SpecUtils::SpecFile &spec )
  {
    // Right now we will only use derived data from Verifinder detectors, since they will show
    //  up as searchmode data, but their derived data is what we would sum anyway
    
    if( !spec.num_measurements() )
      return false;
    
    bool potentially_use = false;
    switch( spec.detector_type() )
    {
      case SpecUtils::DetectorType::VerifinderNaI:
      case SpecUtils::DetectorType::VerifinderLaBr:
        //We'll use derived data for the Verifinder, if we have it
        potentially_use = spec.contains_derived_data();
        break;
        
      default:
        // For all other systems, we will only consider using derived data, if thats the only data
        //  we have.  The meaning of derived data is not well-specified in N42 files, so we should
        //  probably manually inspect contents of systems before using derived data from them.
        potentially_use = (spec.contains_derived_data() && !spec.contains_non_derived_data());
        break;
    }//switch( spec->detector_type() )
    
    return potentially_use;
  }//potentially_analyze_derived_data(...)
  
  /** Function taken from full-spectrum 20231002 */
  void get_derived_measurements( const SpecUtils::SpecFile &spec,
                                set<shared_ptr<const SpecUtils::Measurement>> &foreground,
                                set<shared_ptr<const SpecUtils::Measurement>> &background )
  {
    foreground.clear();
    background.clear();
    
    try
    {
      for( const auto &m : spec.measurements() )
      {
        if( !m || !m->derived_data_properties() || m->num_gamma_channels() < 32 )
          continue;
        
        const uint32_t properties = m->derived_data_properties();
        assert( properties & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::IsDerived) );
        const bool ioi_sum = (properties & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ItemOfInterestSum));
        const bool for_ana = (properties & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::UsedForAnalysis));
        const bool processed = (properties & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::ProcessedFurther));
        const bool back_sub = (properties & static_cast<uint32_t>(SpecUtils::Measurement::DerivedDataProperties::BackgroundSubtracted));
        
        if( back_sub || processed )
          continue;
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::Foreground:
            foreground.insert( m );
            break;
            
          case SpecUtils::SourceType::Background:
            background.insert( m );
            break;
            
          case SpecUtils::SourceType::Unknown:
            //This makes it so the order of seeing Foreground marked record an a IOI sum matters
            //  ... whatever for now
            if( ioi_sum && foreground.empty() )
              foreground.insert( m );
            break;
            
          case SpecUtils::SourceType::IntrinsicActivity:
          case SpecUtils::SourceType::Calibration:
            break;
        }//switch( m->source_type() )
      }//for( const auto &m : derived->measurements() )
      
      if( foreground.size() > 1 )
        throw runtime_error( "Multiple foreground" );
      
      if( background.size() > 1 )
        throw runtime_error( "Multiple background" );
      
      if( foreground.empty() )
        throw runtime_error( "No foreground in derived data" );
      
      if( background.empty() )
        throw runtime_error( "No background in derived data" );
    }catch( std::exception &e )
    {
      foreground.clear();
      background.clear();
    }//try / catch to get derived data spectra to use
  }//get_derived_measurements
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
  contains_background = contains_calibration = contains_intrinsic = false;
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
  
  start_time_ioi = std::time_t(0);

  event_xml_filter_values.clear();

  // FARM fields
  farm_peaks_json.clear();
  farm_peaks_full_json.clear();
  gadras_rid_json.clear();
  isotopics_result_json.clear();
  spectrum_mean = 0.0;
  spectrum_variance = 0.0;
  spectrum_skewness = 0.0;
  spectrum_kurtosis = 0.0;
  farm_min_channel_with_data = -1;
  farm_max_channel_with_data = -1;
  farm_foreground_total_gamma_counts = 0.0;
  farm_foreground_num_gamma_channels = 0;
  farm_foreground_has_neutrons = false;
  farm_foreground_neutron_count = 0.0;
  farm_foreground_neutron_live_time = 0.0f;
  farm_foreground_min_gamma_count = 0.0f;
  farm_foreground_max_gamma_count = 0.0f;
  farm_energy_cal_json.clear();
}//void SpecFileInfoToQuery::reset()


void SpecFileInfoToQuery::fill_info_from_file( const std::string filepath, const Farm::FarmOptions &farm_options )
{
  reset();
  
  file_path = filepath;
  is_file = SpecUtils::is_file(filepath);
  if( !is_file )
    return;
  
  file_size = SpecUtils::file_size(filepath);
  file_path_hash = std::hash<std::string>()(filepath);
  
  // TODO: For files that are not N42 or XML files (as determined by filename), and less than a number of MB, read the file into memory, and then parse.  This will minimize hard-drive seeks (for spinning drives), and allow multithreaded searches to work better - hopefully.
  
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
  
  meas.set_filename( SpecUtils::filename(filepath) );
  
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
  
  set<std::time_t> fore_start_times, back_start_times, cal_start_times;
  set<std::time_t> intrinsic_start_times, unknown_start_times;
  
  // Usually if you think of measurement time, you think of it as IOI time, and you dont want
  //  to include the included background (that you may or may not want to use), cal, or intrinsic
  //  radiation spectra.
  // Also, you dont want to double count between derived and non-derived data.
  // However, you dont always get a spectrum marked as foreground, or whatever, so we will go
  //  through some trouble to do our best to choose the most reasonable time from an explosion
  //  of options :(
  double all_live_time = 0, all_real_time = 0;
  double fore_live_time = 0, fore_real_time = 0;
  double back_live_time = 0, back_real_time = 0;
  double cal_live_time = 0, cal_real_time = 0;
  double intrinsic_live_time = 0, intrinsic_real_time = 0;
  double unknown_live_time = 0, unknown_real_time = 0;
  
  double derived_fore_live_time = 0, derived_fore_real_time = 0;
  double derived_back_live_time = 0, derived_back_real_time = 0;
  double derived_cal_live_time = 0, derived_cal_real_time = 0;
  double derived_intrinsic_live_time = 0, derived_intrinsic_real_time = 0;
  double derived_unknown_live_time = 0, derived_unknown_real_time = 0;
  
  
  map<int,vector<shared_ptr<const SpecUtils::Measurement>>> samples_to_meas;
  
  // If this is a system we want derived data from, we will go all-in and
  //  only look at derived data; most fields will then only show derived data,
  //  but some, like number of samples/measurements will show all..
  const bool use_derived = potentially_analyze_derived_data( meas );
  if( use_derived )
  {
    set<shared_ptr<const SpecUtils::Measurement>> foreground, background;
    get_derived_measurements( meas, foreground, background );
    for( const auto &m : foreground )
      samples_to_meas[m->sample_number()].push_back( m );
    for( const auto &m : background )
      samples_to_meas[m->sample_number()].push_back( m );
    
    
  }//if( use_derived )
  
  if( !use_derived || (samples_to_meas.empty()) )
  {
    for( const int sample_num : meas.sample_numbers() )
    {
      for( const string &det_name : meas.detector_names() )
      {
        const auto m = meas.measurement( sample_num, det_name );
        if( m && !m->derived_data_properties() )
          samples_to_meas[m->sample_number()].push_back( m );
      }//for( const string &det_name : meas.detector_names() )
    }//for( const int sample_num : meas.sample_numbers() )
  }//if( !use_derived || (samples_to_meas.empty()) )
  
  for( const pair<int,vector<shared_ptr<const SpecUtils::Measurement>>> &sample_meas : samples_to_meas )
  {
    vector<float> live_times, real_times;
    bool is_fore = false, is_back = false, is_cal = false, is_intrinsic = false, is_unknown = false;
    size_t num_meas_this_sample = 0;
    
    // In principle as we move through the detectors for a sample number, is_fore, is_back, is_cal,
    //  and is_intrinsic _could_ change, but _shouldn't_, but we'll loop through and check anyway
    bool is_derived_sample = false;
    for( const shared_ptr<const SpecUtils::Measurement> m : sample_meas.second )
    {
      if( m && m->derived_data_properties() )
        is_derived_sample = true;
    }
    
    for( const shared_ptr<const SpecUtils::Measurement> m : sample_meas.second )
    {
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity: is_intrinsic = true; break;
        case SpecUtils::SourceType::Calibration:       is_cal = true;       break;
        case SpecUtils::SourceType::Background:        is_back = true;      break;
        case SpecUtils::SourceType::Foreground:        is_fore = true;      break;
        case SpecUtils::SourceType::Unknown:           is_unknown = true;   break;
      }//switch( m->source_type() )
      
      record_remarks.insert( m->remarks().begin(), m->remarks().end() );
      if( m->num_gamma_channels() > 6 )//skip over GMTubes and other gross-count gamma detectors
      {
        num_meas_this_sample += 1;
        energy_cal_types.insert( m->energy_calibration_model() );
        
        if( m->live_time() > 0.0f )
        {
          live_times.push_back( m->live_time() );
          individual_spectrum_live_time.insert( m->live_time() );
        }
        
        if( m->real_time() > 0.0f )
        {
          real_times.push_back( m->real_time() );
          individual_spectrum_real_time.insert( m->real_time() );
        }
        
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
      
      if( !SpecUtils::is_special(m->start_time()) )
      {
        const std::time_t this_time = chrono::system_clock::to_time_t( m->start_time() );
        start_times.insert( this_time );
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::IntrinsicActivity:
            intrinsic_start_times.insert( this_time );
            break;
            
          case SpecUtils::SourceType::Calibration:
            cal_start_times.insert( this_time );
            break;
            
          case SpecUtils::SourceType::Background:
            back_start_times.insert( this_time );
            break;
            
          case SpecUtils::SourceType::Foreground:
            fore_start_times.insert( this_time );
            break;
            
          case SpecUtils::SourceType::Unknown:
            unknown_start_times.insert( this_time );
            break;
        }//switch( m->source_type() )
      }//if( !SpecUtils::is_special(m->start_time()) )
    }//for( const string &det_name : meas.detector_names() )
    
    contains_background  |= is_back;
    contains_calibration |= is_cal;
    contains_intrinsic   |= is_intrinsic;
      
    float avrg_live_time = std::accumulate( begin(live_times), end(live_times), 0.0f );
    if( avrg_live_time > 0.0f )
      avrg_live_time /= live_times.size();
    
    float avrg_real_time = std::accumulate( begin(real_times), end(real_times), 0.0f );
    if( avrg_real_time > 0.0f )
      avrg_real_time /= real_times.size();
    
    all_live_time += avrg_live_time;
    all_real_time += avrg_real_time;
    
    if( is_fore )
    {
      if( is_derived_sample )
      {
        derived_fore_live_time += avrg_live_time;
        derived_fore_real_time += avrg_real_time;
      }else
      {
        fore_live_time += avrg_live_time;
        fore_real_time += avrg_real_time;
      }
    }else if( is_unknown && !is_intrinsic && !is_cal && !is_back )
    {
      if( is_derived_sample )
      {
        derived_unknown_live_time += avrg_live_time;
        derived_unknown_real_time += avrg_real_time;
      }else
      {
        unknown_live_time += avrg_live_time;
        unknown_real_time += avrg_real_time;
      }
    }else if( is_back )
    {
      if( is_derived_sample )
      {
        derived_back_live_time += avrg_live_time;
        derived_back_real_time += avrg_real_time;
      }else
      {
        back_live_time += avrg_live_time;
        back_real_time += avrg_real_time;
      }
    }else if( is_cal )
    {
      if( is_derived_sample )
      {
        derived_cal_live_time += avrg_live_time;
        derived_cal_real_time += avrg_real_time;
      }else
      {
        cal_live_time += avrg_live_time;
        cal_real_time += avrg_real_time;
      }
    }else if( is_intrinsic )
    {
      if( is_derived_sample )
      {
        derived_intrinsic_live_time += avrg_live_time;
        derived_intrinsic_real_time += avrg_real_time;
      }else
      {
        intrinsic_live_time += avrg_live_time;
        intrinsic_real_time += avrg_real_time;
      }
    }else
    {
      assert( 0 );
    }
  }//for( const int sample_num : meas.sample_numbers() )
  
  if( fore_live_time > 0 )
    total_livetime += fore_live_time;
  else if( derived_fore_live_time > 0 )
    total_livetime += derived_fore_live_time;
  else if( unknown_live_time > 0 )
    total_livetime += unknown_live_time;
  else if( derived_unknown_live_time > 0 )
    total_livetime += derived_unknown_live_time;
  else if( back_live_time > 0 )
    total_livetime += back_live_time;
  else if( cal_live_time > 0 )
    total_livetime += cal_live_time;
  else if( derived_cal_live_time > 0 )
    total_livetime += derived_cal_live_time;
  else if( intrinsic_live_time > 0 )
    total_livetime += intrinsic_live_time;
  else if( derived_intrinsic_live_time > 0 )
    total_livetime += derived_intrinsic_live_time;
  else
  {
    assert( all_live_time == 0 );
    total_livetime += all_live_time;
  }
  
  
  if( fore_real_time > 0 )
    total_realtime += fore_real_time;
  else if( derived_fore_real_time > 0 )
    total_realtime += derived_fore_real_time;
  else if( unknown_real_time > 0 )
    total_realtime += unknown_real_time;
  else if( derived_unknown_real_time > 0 )
    total_realtime += derived_unknown_real_time;
  else if( back_real_time > 0 )
    total_realtime += back_real_time;
  else if( cal_real_time > 0 )
    total_realtime += cal_real_time;
  else if( derived_cal_real_time > 0 )
    total_realtime += derived_cal_real_time;
  else if( intrinsic_real_time > 0 )
    total_realtime += intrinsic_real_time;
  else if( derived_intrinsic_real_time > 0 )
    total_realtime += derived_intrinsic_real_time;
  else
  {
    assert( all_real_time == 0 );
    total_realtime += all_real_time;
  }
  

  if( !fore_start_times.empty() )
    start_time_ioi = *begin(fore_start_times);
  else if( !unknown_start_times.empty() )
    start_time_ioi = *begin(unknown_start_times);
  else if( !back_start_times.empty() )
    start_time_ioi = *begin(back_start_times);
  else if( !cal_start_times.empty() )
    start_time_ioi = *begin(cal_start_times);
  else if( !intrinsic_start_times.empty() )
    start_time_ioi = *begin(intrinsic_start_times);
  else if( !start_times.empty() )
    start_time_ioi = *begin(start_times);


  // ============= Create Combined Foreground/Background for FARM Analysis =============
  // These are needed for peak search, GADRAS, isotopics, and statistical moments

  std::shared_ptr<SpecUtils::Measurement> farm_foreground;
  std::shared_ptr<SpecUtils::Measurement> farm_background;

  if( farm_options.enable_farm_analysis )
  {
    // Collect measurements by source type
    std::set<int> foreground_sample_nums, background_sample_nums, unknown_sample_nums;

    for( const int sample_num : meas.sample_numbers() )
    {
      // Get source type for this sample - check first measurement with this sample
      SpecUtils::SourceType src_type = SpecUtils::SourceType::Unknown;
      for( const std::shared_ptr<const SpecUtils::Measurement> &m : meas.sample_measurements(sample_num) )
      {
        if( m && m->num_gamma_channels() > 6 )
        {
          src_type = m->source_type();
          break;
        }
      }

      switch( src_type )
      {
        case SpecUtils::SourceType::Foreground:
          foreground_sample_nums.insert( sample_num );
          break;
        case SpecUtils::SourceType::Background:
          background_sample_nums.insert( sample_num );
          break;
        case SpecUtils::SourceType::Unknown:
          unknown_sample_nums.insert( sample_num );
          break;
        default:
          break; // Skip calibration, intrinsic, etc.
      }
    }//for( sample_num )

    // Use unknown samples as foreground if no foreground samples exist
    if( foreground_sample_nums.empty() )
      foreground_sample_nums = unknown_sample_nums;

    // If still no foreground, or more than one foreground sample, take just the first, for right now
    if( foreground_sample_nums.size() > 1 )
    {
      const int first_sample = *foreground_sample_nums.begin();
      foreground_sample_nums.clear();
      foreground_sample_nums.insert( first_sample );
    }

    // Only use background if exactly one sample
    if( background_sample_nums.size() != 1 )
      background_sample_nums.clear();

    // Sum all detectors for the selected foreground sample(s)
    if( !foreground_sample_nums.empty() )
    {
      std::vector<std::shared_ptr<const SpecUtils::Measurement>> fore_meass;
      for( const int sample_num : foreground_sample_nums )
      {
        for( const std::shared_ptr<const SpecUtils::Measurement> &m : meas.sample_measurements(sample_num) )
        {
          if( m && m->num_gamma_channels() > 6 )
            fore_meass.push_back( m );
        }
      }

      if( !fore_meass.empty() )
      {
        farm_foreground = meas.sum_measurements( foreground_sample_nums,
                                                  meas.detector_names(),
                                                  nullptr );
        if( farm_foreground )
          farm_foreground->set_source_type( SpecUtils::SourceType::Foreground );
      }
    }//if( !foreground_sample_nums.empty() )

    // Sum all detectors for the selected background sample
    if( !background_sample_nums.empty() )
    {
      farm_background = meas.sum_measurements( background_sample_nums,
                                                meas.detector_names(),
                                                nullptr );
      if( farm_background )
        farm_background->set_source_type( SpecUtils::SourceType::Background );
    }
  }//if( farm_options.enable_farm_analysis )

  // Load a default DRF for peak search and isotopics (background-thread safe)
  std::shared_ptr<DetectorPeakResponse> farm_drf;
  if( farm_options.enable_farm_analysis && farm_foreground )
  {
    SpecUtils::DetectorType det_type = detector_type;
    if( det_type == SpecUtils::DetectorType::Unknown )
      det_type = SpecMeas::guessDetectorTypeFromFileName( filename );

    if( !farm_drf )
    {
      try
      {
        farm_drf = DrfSelect::initARelEffDetector( det_type, manufacturer, model );
      }catch( std::exception & )
      {
      }
    }//if( !farm_drf )
    
    if( !farm_drf )
    {
      try
      {
        farm_drf = DrfSelect::initAGadrasDetector( det_type, nullptr ); //TODO: Doesnt have access to user "GadrasDRFPath" preference, so may not be a great search
      }catch( std::exception & )
      {
      }
    }//if( !farm_drf )
  }//if( farm_options.enable_farm_analysis && farm_foreground )


  // ============= Calculate Statistical Moments =============
  if( farm_options.enable_farm_analysis
      && farm_foreground && farm_foreground->num_gamma_channels() > 0 )
  {
    const std::shared_ptr<const std::vector<float>> gamma_counts = farm_foreground->gamma_counts();
    const std::shared_ptr<const SpecUtils::EnergyCalibration> cal =
        farm_foreground->energy_calibration();

    if( cal && cal->valid() && gamma_counts && !gamma_counts->empty() )
    {
      const size_t num_channel = gamma_counts->size();

      // First pass: calculate mean and total weight
      double sum_w = 0.0;
      double sum_wx = 0.0;

      for( size_t i = 0; i < num_channel; ++i )
      {
        const double y = static_cast<double>( (*gamma_counts)[i] );
        if( y <= 0.0 )
          continue;

        const double x = 0.5 * (cal->energy_for_channel(i) + cal->energy_for_channel(i + 1));
        sum_w += y;
        sum_wx += y * x;
      }

      if( sum_w > 0.0 )
      {
        spectrum_mean = sum_wx / sum_w;

        // Second pass: calculate variance (and accumulate for higher moments)
        double sum_diff2 = 0.0;
        double sum_diff3 = 0.0;
        double sum_diff4 = 0.0;

        for( size_t i = 0; i < num_channel; ++i )
        {
          const double y = static_cast<double>( (*gamma_counts)[i] );
          if( y <= 0.0 )
            continue;

          const double x = 0.5 * (cal->energy_for_channel(i) + cal->energy_for_channel(i + 1));
          const double diff = x - spectrum_mean;
          const double diff2 = diff * diff;
          const double diff3 = diff2 * diff;
          const double diff4 = diff2 * diff2;

          sum_diff2 += y * diff2;
          sum_diff3 += y * diff3;
          sum_diff4 += y * diff4;
        }

        spectrum_variance = sum_diff2 / sum_w;

        if( spectrum_variance > 0.0 )
        {
          const double std_dev = std::sqrt( spectrum_variance );
          const double std_dev3 = std_dev * std_dev * std_dev;
          const double std_dev4 = std_dev3 * std_dev;

          spectrum_skewness = (sum_diff3 / sum_w) / std_dev3;
          spectrum_kurtosis = (sum_diff4 / sum_w) / std_dev4 - 3.0;  // Excess kurtosis
        }
      }//if( sum_w > 0.0 )
    }//if( valid calibration and counts )
  }//if( farm_foreground )


  // ============= FARM: Foreground Summary + Energy Cal =============
  if( farm_options.enable_farm_analysis && farm_foreground )
  {
    const std::shared_ptr<const std::vector<float>> gamma_counts = farm_foreground->gamma_counts();
    farm_foreground_num_gamma_channels = static_cast<int>( farm_foreground->num_gamma_channels() );
    farm_foreground_total_gamma_counts = farm_foreground->gamma_count_sum();

    // Channel range and min/max counts
    if( gamma_counts && !gamma_counts->empty() )
    {
      float min_nonzero = std::numeric_limits<float>::max();
      float max_count   = 0.0f;
      int   first_ch    = -1;
      int   last_ch     = -1;

      for( size_t i = 0; i < gamma_counts->size(); ++i )
      {
        const float c = (*gamma_counts)[i];
        if( c > 0.0f )
        {
          if( first_ch < 0 )
            first_ch = static_cast<int>( i );
          last_ch = static_cast<int>( i );
          if( c < min_nonzero )
            min_nonzero = c;
          if( c > max_count )
            max_count = c;
        }
      }

      farm_min_channel_with_data  = first_ch;
      farm_max_channel_with_data  = last_ch;
      farm_foreground_min_gamma_count = (first_ch >= 0) ? min_nonzero : 0.0f;
      farm_foreground_max_gamma_count = max_count;
    }

    // Neutron info
    farm_foreground_has_neutrons = farm_foreground->contained_neutron();
    if( farm_foreground_has_neutrons )
    {
      const std::vector<float> &ncounts = farm_foreground->neutron_counts();
      double nsum = 0.0;
      for( float v : ncounts )
        nsum += v;
      farm_foreground_neutron_count    = nsum;
      farm_foreground_neutron_live_time = farm_foreground->neutron_live_time();
    }

    // Energy calibration JSON
    const std::shared_ptr<const SpecUtils::EnergyCalibration> cal =
        farm_foreground->energy_calibration();

    if( cal && cal->valid() )
    {
      nlohmann::json cal_json;
      cal_json["min_energy"] = cal->lower_energy();
      cal_json["max_energy"] = cal->upper_energy();

      const SpecUtils::EnergyCalType cal_type = cal->type();

      if( cal_type == SpecUtils::EnergyCalType::Polynomial
          || cal_type == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
      {
        const std::vector<float> &coeffs = cal->coefficients();
        cal_json["polynomial_energy_coeffs"] = std::vector<double>( coeffs.begin(), coeffs.end() );
      }
      else if( cal_type == SpecUtils::EnergyCalType::FullRangeFraction )
      {
        const std::vector<float> poly_coeffs =
            SpecUtils::fullrangefraction_coef_to_polynomial( cal->coefficients(),
                                                            cal->num_channels() );
        cal_json["polynomial_energy_coeffs"] = std::vector<double>( poly_coeffs.begin(),
                                                                    poly_coeffs.end() );
      }
      else if( cal_type == SpecUtils::EnergyCalType::LowerChannelEdge )
      {
        const size_t nch = cal->num_channels();
        if( nch > 0 )
          cal_json["average_amp_gain"] = static_cast<double>(
              (cal->upper_energy() - cal->lower_energy()) / static_cast<float>( nch ) );
      }

      // Deviation pairs (present on Polynomial, FRF, or DefaultPolynomial)
      const std::vector<std::pair<float,float>> &dev_pairs = cal->deviation_pairs();
      if( !dev_pairs.empty() )
      {
        nlohmann::json dp_array = nlohmann::json::array();
        for( const std::pair<float,float> &dp : dev_pairs )
          dp_array.push_back( nlohmann::json::array( { dp.first, dp.second } ) );
        cal_json["nonlinear_deviation_pairs"] = dp_array;
      }

      farm_energy_cal_json = cal_json.dump();
    }//if( cal && cal->valid() )
  }//if( farm_options.enable_farm_analysis && farm_foreground )


  // ============= FARM: Peak Search =============
  if( farm_options.enable_farm_analysis && farm_foreground )
  {
    const bool is_hpge = (detector_type == SpecUtils::DetectorType::Fulcrum
                          || detector_type == SpecUtils::DetectorType::Fulcrum40h);
    Farm::perform_peak_search( farm_foreground, farm_drf, is_hpge, *this );
  }

  // ============= FARM: GADRAS Full Spectrum Isotope ID =============
  if( farm_options.enable_farm_analysis && farm_options.enable_gadras_rid
      && !farm_options.gadras_exe_path.empty() && farm_foreground )
  {
    // Build a minimal SpecFile with just foreground (and optionally background) for GADRAS
    std::shared_ptr<SpecUtils::SpecFile> gadras_spec = std::make_shared<SpecUtils::SpecFile>( meas );
    gadras_spec->remove_measurements( gadras_spec->measurements() );
    gadras_spec->add_measurement( farm_foreground, !farm_background );
    if( farm_background )
      gadras_spec->add_measurement( farm_background, true );

    gadras_rid_json = Farm::run_gadras_full_spectrum_id_analysis(
        farm_options.gadras_exe_path,
        gadras_spec,
        detector_type,
        farm_options.synthesize_background_if_missing );
  }

  
  nlohmann::json isotopics_json = nlohmann::json::array();
  
  // ============= FARM: Isotopics via RelActCalcAuto =============
  const bool do_u_enrich = Farm::should_do_uranium_isotopics( farm_peaks_json );
  const bool do_pu_enrich = Farm::should_do_plutonium_isotopics( farm_peaks_json );
  
  
  if( farm_options.enable_farm_analysis && farm_options.enable_relact_isotopics
     && farm_foreground && (do_u_enrich || do_pu_enrich) )
  {
    // Uranium isotopics: requires 185.7 keV + (205.3 or 143.8 keV) peaks
    if( do_u_enrich && !do_pu_enrich )
    {
      const std::string u_opts_path = SpecUtils::append_path( InterSpec::staticDataDirectory(), "rel_act/HPGe U (120-1001 keV).xml" );
      const Farm::EnrichmentResults u_result = Farm::run_relact_isotopics(
          farm_foreground, farm_background,
          std::vector<std::shared_ptr<const PeakDef>>{}, u_opts_path, farm_drf );
      isotopics_json.push_back( u_result.toJson() );
    }else if( do_pu_enrich && !do_u_enrich )
    {
      const std::string pu_opts_path = SpecUtils::append_path( InterSpec::staticDataDirectory(), "rel_act/HPGe Pu (120-780 keV).xml" );
      const Farm::EnrichmentResults pu_result = Farm::run_relact_isotopics(
          farm_foreground, farm_background,
          std::vector<std::shared_ptr<const PeakDef>>{}, pu_opts_path, farm_drf );
      isotopics_json.push_back( pu_result.toJson() );
    }else
    {
      cerr << "TODO: Do not have MOX implemented for RelAct" << endl;
    }
  }
  
  
  // ============= FARM: Isotopics via FRAM =============
  if( farm_options.enable_farm_analysis && farm_options.enable_fram_isotopics
      && !farm_options.fram_exe_path.empty() && farm_foreground
     && (do_u_enrich || do_pu_enrich) )
  {
    
    // Write a temp N42 for FRAM input
    const std::string fram_tmp = SpecUtils::temp_file_name( "farm_fram_", SpecUtils::temp_dir() );
    {
      std::ofstream fram_out( fram_tmp, std::ios::binary );
      if( fram_out.is_open() )
      {
        SpecUtils::SpecFile fram_spec = meas;
        fram_spec.remove_measurements( fram_spec.measurements() );
        fram_spec.add_measurement( farm_foreground, !farm_background );
        if( farm_background )
          fram_spec.add_measurement( farm_background, true );
        fram_spec.write_2012_N42( fram_out );
      }
    }

    assert( SpecUtils::is_file( fram_tmp ) );
    
    if( SpecUtils::is_file( fram_tmp ) )
    {
      const Farm::EnrichmentResults fram_result = Farm::run_fram_isotopics(
            farm_options.fram_exe_path, farm_options.fram_output_path, fram_tmp, do_u_enrich, do_pu_enrich );
      isotopics_json.push_back( fram_result.toJson() );
    }//if( SpecUtils::is_file( fram_tmp ) )
    
    SpecUtils::remove_file( fram_tmp );
  }
  
  if( !isotopics_json.empty() )
    isotopics_result_json = isotopics_json.dump();

  // ============= FARM: Write Fertilized N42 =============
  if( farm_options.enable_farm_analysis && farm_options.write_fertilized_n42 && farm_foreground )
  {
    SpecUtils::SpecFile output_spec = meas;

    Farm::write_fertilized_n42( filepath, output_spec, *this );
  }

}//void fill_info_from_file( const std::string filepath, ... )


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
    const string prefix = m_farm_options.enable_farm_analysis
                          ? "interspec_file_query_FARM" : "interspec_file_query";
    string db_location = SpecUtils::temp_file_name( prefix, SpecUtils::temp_dir() );
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


std::string SpecFileQueryDbCache::construct_persisted_db_filename() const
{
  string canonical_path = m_fs_path;
  SpecUtils::make_canonical_path( canonical_path );

  const size_t path_hash = std::hash<std::string>()( canonical_path );
  string fname = "InterSpec_file_query_cache_" + std::to_string( path_hash );

  if( m_farm_options.enable_farm_analysis )
    fname += "_FARM";

  fname += "_v3.sqlite3";

  return SpecUtils::append_path( canonical_path, fname );
}//std::string construct_persisted_db_filename()


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

  const string persisted_path = construct_persisted_db_filename();
  
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
  if( rhs.contains_background != lhs.contains_background )
    return false;
  if( rhs.contains_calibration != lhs.contains_calibration )
    return false;
  if( rhs.contains_intrinsic != lhs.contains_intrinsic )
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
  if( rhs.start_time_ioi != lhs.start_time_ioi )
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
      dbinforaw->fill_info_from_file( filename, m_farm_options );
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
    
    const string newfilename = construct_persisted_db_filename();
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


void SpecFileQueryDbCache::set_farm_options( const Farm::FarmOptions &opts )
{
  if( m_farm_options == opts )
    return;

  m_farm_options = opts;

  if( m_use_db_caching )
    resetCache();
}//void set_farm_options()


void SpecFileQueryDbCache::resetCache()
{
  // Wait for any in-progress caching to finish, then prevent new caching
  {
    std::unique_lock<std::mutex> lock( m_cv_mutex );
    m_stop_caching = true;
    if( m_doing_caching )
      m_cv.wait( lock, [this]() -> bool { return !m_doing_caching; } );
  }

  {
    std::lock_guard<std::mutex> lock( m_db_mutex );

    const string old_location = m_db_location;
    m_db_session.reset();
    m_db.reset();

    // Remove the old temp DB file; persisted DBs are kept (keyed by _FARM suffix)
    if( !m_using_persist_caching && !old_location.empty() )
      SpecUtils::remove_file( old_location );

    // Try to open an existing persisted DB at the new name
    if( m_using_persist_caching && init_existing_persisted_db() )
    {
      // Re-allow caching so cache_results can proceed
      std::lock_guard<std::mutex> cv_lock( m_cv_mutex );
      m_stop_caching = false;
      return;
    }

    // Otherwise create a fresh DB (persisted or temp)
    string new_location;
    if( m_using_persist_caching )
      new_location = construct_persisted_db_filename();
    else
    {
      const string prefix = m_farm_options.enable_farm_analysis
                            ? "interspec_file_query_FARM" : "interspec_file_query";
      new_location = SpecUtils::temp_file_name( prefix, SpecUtils::temp_dir() );
    }

    open_db( new_location, true );
  }

  // Re-allow caching so cache_results can proceed
  {
    std::lock_guard<std::mutex> lock( m_cv_mutex );
    m_stop_caching = false;
  }
}//void resetCache()


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
    info->fill_info_from_file( filepath, m_farm_options );
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
        info->fill_info_from_file( filepath, m_farm_options );
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
    info->fill_info_from_file( filepath, m_farm_options );
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
    cerr << "Caught Dbo::Exception in spec_file_info: '" << e.what() <<"', backend code: '"
         << e.code() << "'" << endl;
    info->fill_info_from_file( filepath, m_farm_options );
    info->fill_event_xml_filter_values( filepath, m_xmlfilters );
  }catch( std::exception &e )
  {
    cerr << "Caught std::Exception in spec_file_info: " << e.what() << endl;
    info->fill_info_from_file( filepath, m_farm_options );
    info->fill_event_xml_filter_values( filepath, m_xmlfilters );
  }
  
  return info;
}//SpecFileInfoToQuery spec_file_info( const std::string &filepath )



