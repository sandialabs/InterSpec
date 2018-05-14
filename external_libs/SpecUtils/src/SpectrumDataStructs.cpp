/**
 SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 Copyright (C) 2016 William Johnson
 
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


#include "SpecUtils_config.h"

#if( USE_D3_EXPORTING )
#include "D3SupportFiles.h"
#endif //#if( USE_D3_EXPORTING )

#include <ctime>
#include <vector>
#include <memory>
#include <string>
#include <cctype>
#include <locale>
#include <limits>
#include <numeric>
#include <fstream>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <algorithm>
#include <functional>
//#define __STDC_FORMAT_MACROS
//#include <inttypes.h>
#include <sys/stat.h>

#include <boost/functional/hash.hpp>

#if( BUILD_AS_UNIT_TEST_SUITE )
#include <boost/test/unit_test.hpp>
#endif

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/UtilityFunctions.h"
#include "SpecUtils/SpectrumDataStructs.h"

#if( ENABLE_D3_CHART_EXPORTING )
#include "SpecUtils/D3SpectrumExport.h"
#endif

using namespace std;


using UtilityFunctions::trim;
using UtilityFunctions::split;
using UtilityFunctions::iequals;
using UtilityFunctions::to_lower;
using UtilityFunctions::to_upper;
using UtilityFunctions::contains;
using UtilityFunctions::icontains;
using UtilityFunctions::starts_with;
using UtilityFunctions::istarts_with;
using UtilityFunctions::ireplace_all;

using UtilityFunctions::time_from_string;

//For references in comments similar to 'refDXMA3HSRA6', see
//  documentation/comment_refernce_to_ouo_source.txt for source file information

#if( !SpecUtils_EXTERNALLY_DEFINED_LOG_MESSAGE )
#define passMessage(m,s,p) std::cerr << s << ": " << m << std::endl;

#else
#define passMessage(m,s,p) log_error_message(m,s,p)

extern void log_error_message( const std::string &message, const std::string &source, const int priority );
#endif


//Including DetectiveModelFromSerial.hpp is done last since it relies on having
//  all the 'using' statments and other code included already 
#include "SpecUtils/DetectiveModelFromSerial.hpp"


namespace
{
  //Sometimes a detector wont have a name, but we still need to refer to it in
  //  2012 XML files in various places, so we'll use s_unnamed_det_placeholder
  //  to do this.  A side effect to this is we have to be careful to use it
  //  everywhere, and account for it when reading it back in.
  const std::string s_unnamed_det_placeholder = "unamed";
  
  
  template <class T>
  std::istream &readBinaryData( std::istream &input, T &val )
  {
    input.read( (char *)&val, sizeof(T) );
    return input;
  }
  
  //writeBinaryData(): returns number of bytes written
  template <class T>
  size_t writeBinaryData( std::ostream &input, T &val )
  {
    input.write( (const char *)&val, sizeof(T) );
    return sizeof(T);
  }
  
  bool toInt( const std::string &str, int &f )
  {
    const int nconvert = sscanf( str.c_str(), "%i", &f );
    return (nconvert == 1);
  }
  
  bool toFloat( const std::string &str, float &f )
  {
    //ToDO: should probably use UtilityFunctions::parse_float(...) for consistency/speed
    const int nconvert = sscanf( str.c_str(), "%f", &f );
    return (nconvert == 1);
  }
  
  bool toDouble( const std::string &str, double &f )
  {
    const int nconvert = sscanf( str.c_str(), "%lf", &f );
    return (nconvert == 1);
  }
  
  /** Adds the vector of 'input' float vectors, to results.
      The result will be resized to larges input float vector size.
   */
  void add_to( vector<float> &results,
              const vector< std::shared_ptr<const vector<float> > > &input )
  {
    results.clear();
    if( input.empty() )
      return;
    
    size_t max_size = 0;
    for( const auto &currarayptr : input )
      max_size = std::max( max_size, currarayptr->size() );
    
    results.resize( max_size, 0.0f );
    
    for( const auto &currarayptr : input )
    {
      const vector<float> &curraray = (*currarayptr);
      
      const size_t size = curraray.size();
      if( size > results.size() )  //should never evaluate to true, but JIC
        results.resize( size, 0.0f );
      
      for( size_t i = 0; i < size; ++i )
        results[i] += curraray[i];
    }
  }//add_to(...)
  
  
  void sum_with_rebin( vector<float> &results,
                      const std::shared_ptr<const Measurement> &binning,
                      const vector<std::shared_ptr<const Measurement> > &datas )
  {
    assert( !!binning );
    
    const ShrdConstFVecPtr &wantedenergies = binning->channel_energies();
    
    const size_t nbin = wantedenergies->size();
    if( results.size() < nbin )
      results.resize( nbin, 0.0f );
    
    assert( !!wantedenergies );
    
    for( size_t i = 0; i < datas.size(); ++i )
    {
      const std::shared_ptr<const Measurement> &d = datas[i];
      const ShrdConstFVecPtr &dataenergies = d->channel_energies();
      const ShrdConstFVecPtr &channel_counts = d->gamma_counts();
      
      if( !dataenergies || !channel_counts )
      {
        cerr << "sum_with_rebin(...): found spectrum with no bin" << endl;
        continue;
      }//if( !dataenergies )
      
      if( dataenergies == wantedenergies )
      {
        for( size_t j = 0; j < nbin; ++j )
          results[j] += (*channel_counts)[j];
      }else if( channel_counts->size() > 3 )
      {
        vector<float> resulting_counts;
        rebin_by_lower_edge( *dataenergies, *channel_counts,
                            *wantedenergies, resulting_counts );
        
        assert( resulting_counts.size() == nbin );
        
        for( size_t j = 0; j < nbin; ++j )
          results[j] += resulting_counts[j];
      }//if( dataenergies == wantedenergies )
      
    }//for( size_t i = 0; i < datas.size(); ++i )
  }//void sum_with_rebin(...)

  
  
  bool not_alpha_numeric(char c)
  {
    return !(std::isalnum(c) || c==' ');
  }
  
  
  double ortecLatOrLongStrToFlt( std::string input )
  {
    input.erase( std::remove_if(input.begin(), input.end(), not_alpha_numeric),
                 input.end());
    trim( input );
    
    char dir;
    float degress, minutes, seconds;
    const size_t npar = sscanf( input.c_str(),"%f %f %f %c",
                                &degress, &minutes, &seconds, &dir );
    
    if( npar == 4 )
    {
      const double sign = ((dir=='N'||dir=='E') ? 1.0 : -1.0);
      return sign * (degress + (minutes/60.0) + (seconds/3600.0));
    }
    
    return -999.9;
  };//ortecLatOrLongStrToFlt(...)
  
  bool dev_pair_less_than( const std::pair<float,float> &lhs,
                           const std::pair<float,float> &rhs )
  {
    return (lhs.first < rhs.first);
  }
  
  //parse_deg_min_sec_lat_lon(...): not real well implemented, only tested on
  //  MicroRaider files
  bool parse_deg_min_sec_lat_lon( const char *str, const size_t len,
                                  double &lat, double &lon )
  {
//    "25°47\"17.820' N / 80°19\"25.500' W"
    lat = lon = -999.9;
    
    if( !str || !len )
      return false;
    
    const char *end = str + len;
    const char *pos = std::find( str, end, '/' );
    if( pos != end )
    {
      string latstr(str,pos), lonstr(pos+1,end);
      for( size_t i = 0; i < latstr.size(); ++i )
        if( !isalnum(latstr[i]) )
          latstr[i] = ' ';
      for( size_t i = 0; i < lonstr.size(); ++i )
        if( !isalnum(lonstr[i]) )
          lonstr[i] = ' ';
      ireplace_all( latstr, "  ", " " );
      ireplace_all( lonstr, "  ", " " );
      
      lat = ortecLatOrLongStrToFlt( latstr );
      lon = ortecLatOrLongStrToFlt( lonstr );
      
      //cerr << "latstr='" << latstr << "'-->" << lat << endl;
      //cerr << "lonstr='" << lonstr << "'-->" << lon << endl;
      
      return (Measurement::valid_longitude(lon)
              && Measurement::valid_latitude(lat));
    }//if( pos != end )
    
    return false;
  }//bool parse_deg_min_sec_lat_lon(
  
  boost::posix_time::ptime datetime_ole_to_posix(double ole_dt)
  {
    static const boost::gregorian::date ole_zero(1899,12,30);
    
    boost::gregorian::days d(static_cast<int64_t>(ole_dt) );
    boost::posix_time::ptime pt(ole_zero + d);
    
    ole_dt -= d.days();
    ole_dt *= 24 * 60 * 60 * 1000;
    
    return pt + boost::posix_time::milliseconds( std::abs( static_cast<int64_t>(ole_dt) ) );
    
    /*
     typedef typename time_type::date_type date_type;
     typedef typename time_type::date_duration_type date_duration_type;
     typedef typename time_type::time_duration_type time_duration_type;
     using boost::math::modf;
     static const date_type base_date(1899, Dec, 30);
     static const time_type base_time(base_date, time_duration_type(0,0,0));
     int dayOffset, hourOffset, minuteOffset, secondOffset;
     double fraction = fabs(modf(oa_date, &dayOffset)) * 24; // fraction = hours
     fraction = modf(fraction, &hourOffset) * 60; // fraction = minutes
     fraction = modf(fraction, &minuteOffset) * 60; // fraction = seconds
     modf(fraction, &secondOffset);
     time_type t(base_time);
     t += time_duration_type(hourOffset, minuteOffset, secondOffset);
     t += date_duration_type(dayOffset);
     return t;
     */
  }
  
  
  /* IAEA block labels that represent items to put into the remarks_ variable of
     MeasurementInfo.
   */
  const char * const ns_iaea_comment_labels[] =
  {
    "Comment", "AcquisitionMode", "CrystalType", "Confidence",
    "MinDoseRate", "MaxDoseRate", "AvgDoseRate", "MinNeutrons", "MaxNeutrons",
    "DetectorLength", "DetectorDiameter", "BuiltInSrcType", "BuiltInSrcActivity",
    "HousingType", "GMType", "He3Pressure", "He3Length", "He3Diameter",
    "ModMaterial", "ModVolume", "ModThickness", "LastSourceStabTime",
    "LastSourceStabFG", "LastCalibTime", "LastCalibSource", "LastCalibFG",
    "LastCalibFWHM", "LastCalibTemp", "StabilType", "StartupStatus",
    "TemperatureBoard", "TemperatureBoardRange", "BatteryVoltage", "Uptime",
    "DoseRate", "DoseRateMax20min", "BackgroundSubtraction",
    "FWHMCCoeff", "ROI", "CalibPoint", "NeutronAlarm",
    "GammaDetector", "NeutronDetector", "SurveyId", "EventNumber",
    "Configuration"
  };//const char * const ns_iaea_comment_labels = {...}
  
  
  /* IAEA block labels that represent information to be put into
     component_versions_ member variable of MeasurementInfo.
   */
  const char * const ns_iaea_version_labels[] =
  {
    "Hardware", "TemplateLibraryVersion", "NativeAlgorithmVersion",
    "ApiVersion", "Firmware", "Operating System", "Application",
    "SoftwareVersion"
  };
  
  string pad_iaea_prefix( string label )
  {
    label.resize( 22, ' ' );
    return label + ": ";
  }
  
  string print_to_iaea_datetime( const boost::posix_time::ptime &t )
  {
    char buffer[256];
    const int day = static_cast<int>( t.date().day() );
    const int month = static_cast<int>( t.date().month() );
    const int year = static_cast<int>( t.date().year() );
    const int hour = static_cast<int>( t.time_of_day().hours() );
    const int minutes = static_cast<int>( t.time_of_day().minutes() );
    const int seconds = static_cast<int>( t.time_of_day().seconds() );
    
    snprintf( buffer, sizeof(buffer), "%02d.%02d.%04d %02d:%02d:%02d",
              day, month, year, hour, minutes, seconds );
    
    return buffer;
  }//print_iaea_datetime(...)
  
  
  //During parsing we abuse the remarks to hold PCF specific information, so
  //  lets extract that back out now - soryy for this horriblness.
  string findPcfRemark( const char *start, const vector<string> &remarks )
  {
    const size_t start_len = strlen(start);
    for( size_t i = 0; i < remarks.size(); ++i )
      if( UtilityFunctions::istarts_with(remarks[i], start) )
      {
        string val = remarks[i].substr( start_len );
        const size_t pos = val.find_first_not_of( " :\t\n\r=");
        if( pos != string::npos )
          val = val.substr(pos);
        return val;
      }
    return "";
  }//string findPcfRemark( const char *start, const vector<string> &remarks )
  
  
  string parse_pcf_field( const string &header, size_t offset, size_t len )
  {
#if(PERFORM_DEVELOPER_CHECKS)
    if( offset+len > header.size() )
    {
      log_developer_error( BOOST_CURRENT_FUNCTION, "Logic error in parse_pcf_field" );
      throw runtime_error( "Logic error in parse_pcf_field" );
    }
#endif
    string field( header.begin() + offset, header.begin() + offset + len );
    const size_t zeropos = field.find_first_of( '\0' );
    if( zeropos != string::npos )
      field = field.substr(0,zeropos);
    UtilityFunctions::trim( field );
    
    return field;
  };//parse_pcf_field
}//anaomous namespace



namespace
{
  //The functions and macros in this namespace (I know, macros dont obey
  //  namespaces, but whatever for now) are a small scale first try towards
  //  using non-destructive XML parsing so memory mapped files can be used as
  //  input to the file reader or whatever.  The goal is to enable reasonably
  //  mistake free parsing, while cutting down the verbosity of the code, and
  //  avoiding unecassary allocations.
  //
  //The general convention is that macros are totally unsafe and dont check for
  //  null pointers or anything, while functions are safe to call with nullptrs.
  //  The functions with the postfix "_nso" are name space optional functions
  //  meant to help with cases where some XML files use a namespace, and others
  //  do not (rapidxml is name-space unaware, so we have to deal with them).
  //
  //TODO: Finish converting older sections of code to use these functions/macros
  
  //Function to determine static c-string length at compile time.
  template<size_t N>
  size_t lengthof(const char (&)[N])
  {
    return N - 1;
  }
  
#define XML_VALUE_COMPARE( node, cstr ) (rapidxml::internal::compare((node)->value(), (node)->value_size(), cstr, lengthof(cstr), true))
#define XML_VALUE_ICOMPARE( node, cstr ) (rapidxml::internal::compare((node)->value(), (node)->value_size(), cstr, lengthof(cstr), false))
#define XML_NAME_COMPARE( node, cstr ) (rapidxml::internal::compare((node)->name(), (node)->name_size(), cstr, lengthof(cstr), true))
#define XML_NAME_ICOMPARE( node, cstr ) (rapidxml::internal::compare((node)->name(), (node)->name_size(), cstr, lengthof(cstr), false))
  
#define XML_FIRST_NODE(node,name)((node)->first_node(name,lengthof(name),true))
#define XML_FIRST_ATTRIB(node,name)((node)->first_attribute(name,lengthof(name),true))

#define XML_FIRST_NODE_CHECKED(node,name)((node) ? (node)->first_node(name,lengthof(name),true) : (rapidxml::xml_node<char> *)0)
#define XML_FIRST_ATTRIB_CHECKED(node,name)((node) ? (node)->first_attribute(name,lengthof(name),true) : (rapidxml::xml_node<char> *)0)
  
#define XML_NEXT_TWIN(node)((node)->next_sibling((node)->name(), (node)->name_size()))
#define XML_NEXT_TWIN_CHECKED(node)((node) ? (node)->next_sibling((node)->name(), (node)->name_size()): (rapidxml::xml_node<char> *)0)
  
  
  template<class Ch>
  std::string xml_value_str( const rapidxml::xml_base<Ch> *node )
  {
    if( !node || !node->value_size() )
      return "";
    return std::string( node->value(), node->value() + node->value_size() );
  }
  
  template<class Ch>
  std::string xml_name_str( const rapidxml::xml_base<Ch> *node )
  {
    if( !node || !node->name_size() )
      return "";
    return std::string( node->name(), node->name() + node->name_size() );
  }
  
  template<size_t n, size_t m>
  const rapidxml::xml_node<char> *xml_first_node_nso( const rapidxml::xml_node<char> *parent, const char (&name)[n], const char (&ns)[m] )
  {
    if( m < 2 )
    {
      return parent ? parent->first_node(name, n-1) : (const rapidxml::xml_node<char> *)0;
    }else
    {
      if( !parent )
        return 0;
      const rapidxml::xml_node<char> *answer = parent->first_node(name, n-1);
      if( !answer )
      {
        char newname[n+m-1];
        memcpy(newname, ns, m-1);
        memcpy(newname + m - 1, name, n);
        answer = parent->first_node(newname, m+n-2);
      }
      
      return answer;
    }
  }
  
  template<size_t n>
  const rapidxml::xml_node<char> *xml_first_node_nso( const rapidxml::xml_node<char> *parent, const char (&name)[n], const std::string &ns )
  {
    if( ns.size() < 2 )
    {
      return parent ? parent->first_node(name, n-1) : (const rapidxml::xml_node<char> *)0;
    }else
    {
      if( !parent )
        return 0;
      const rapidxml::xml_node<char> *answer = parent->first_node(name, n-1);
      if( !answer )
      {
        const std::string newname = ns + name;
        answer = parent->first_node(newname.c_str(), newname.size() );
      }
      
      return answer;
    }
  }
  
  bool xml_value_to_flt( const rapidxml::xml_node<char> *node, float &val )
  {
    val = 0.0f;
    if( !node )
      return false;
    return UtilityFunctions::parse_float( node->value(), node->value_size(), val );
  }
  
  std::string get_n42_xmlns( const rapidxml::xml_node<char> *node )
  {
    const string node_name = xml_name_str(node);
    const size_t colon_pos = node_name.find(':');
    if( colon_pos != string::npos && UtilityFunctions::icontains(node_name, "n42") )
      return node_name.substr( 0, colon_pos + 1 );
    return "";
  }//std::string get_n42_xmlns( const rapidxml::xml_node<char> *node )
}//anonomous namespace for XML utilities


namespace
{
  //anaonomous namespace for functions to help parse N42 files, that wont be
  //  usefull outside of this file
  
  //is_gamma_spectrum(): Trys to determine  if the spectrum_node cooresponds
  //  to a gamma or neutron node in the XML.
  //  Will throw std::runtime_exception if if cant unambiguosly tell if its a
  //  gamma spectrum or not
  bool is_gamma_spectrum( const rapidxml::xml_attribute<char> *detector_attrib,
                         const rapidxml::xml_attribute<char> *type_attrib,
                         const rapidxml::xml_node<char> *det_type_node,
                         const rapidxml::xml_node<char> *spectrum_node )
  {
    bool is_gamma = false, is_nuetron = false;
    
    //The ICD1 Spec it shall be the case that <DetectorType>
    //  will say 'Neutron' for neutron detectors, so if we find this, we
    //  wont bother to mess with names
    if( det_type_node && det_type_node->value_size() )
    {
      const string det_type = xml_value_str(det_type_node);
      if( icontains(det_type,"neutron")
         || icontains(det_type,"GMTube") )
        return false;
      if( icontains(det_type,"Gamma") )
        return true;
    }//if( det_type_node && det_type_node->value_size() )
    
    if( detector_attrib )
    {
      string name = xml_value_str(detector_attrib);
      to_lower( name );
      
      if( contains(name, "neutron") )
        is_nuetron = true;
      
      if( icontains(name, "GMTube") )
        is_nuetron = true;
      
      if( contains(name, "pha") )
        is_gamma = true;
      
      if( contains(name, "gamma") )
        is_gamma = true;
      
      if( name == "tungsten" ) //some FLIR identiFINDER
        return true;
      
      if( !is_nuetron && !is_gamma ) //try to match the name
      {
        const size_t len = name.length();
        bool matches_convention = (len >= 2);
        if( len >= 1 )
        {
          const char c = name[0];
          matches_convention |= (c=='a' || c=='b' || c=='c' || c=='d' );
        }
        if( len >= 2 )
        {
          const char c = name[1];
          matches_convention |= (isdigit(c) || c=='a' || c=='b' || c=='c' || c=='d' );
        }
        if( len >= 3 )
        {
          const char c = name[2];
          matches_convention |= (isdigit(c) || c=='n');
        }
        
        if( matches_convention )
          matches_convention = !UtilityFunctions::icontains( name, "Unknown" );
        
        if( matches_convention )
        {
          const char c = name[name.size()-1];
          is_nuetron = (c == 'n');
          is_gamma = (isdigit(c) != 0);
        }//if( matches_convention )
      }//if( !is_nuetron && !is_gamma ) //try to match the name
    }//if( detector_attrib )
    
    
    if( type_attrib )
    {
      const string name = xml_value_str(type_attrib);
      if( icontains(name, "pha") || icontains(name, "Gamma") )
        is_gamma = true;
    }//if( type_attrib )
    
    
    if( is_nuetron == is_gamma )
    {
      for( const rapidxml::xml_node<char> *node = spectrum_node; node; node = node->parent() )
      {
        const rapidxml::xml_attribute<char> *attrib = node->first_attribute( "DetectorType", 12 );
        const string textstr = xml_value_str(attrib);
        if( textstr.length() )
        {
          is_gamma = (icontains( textstr, "gamma" ) || icontains( textstr, "LaBr" ) || icontains( textstr, "NaI" ) );
          is_nuetron = icontains( textstr, "neutron" );
          break;
        }//if( textstr.length() )
      }//while( (parent = det_type_node->parent()) )
    }//if( is_nuetron == is_gamma )
    
    
    if( is_nuetron == is_gamma )
    {
      //We should probably just assume its a gamma detector here....
      stringstream msg;
      msg << SRC_LOCATION << "\n\tFound spectrum thats ";
      
      if( is_nuetron )
        msg << "a neutron and a gamma spectrum Detector=";
      else
        msg << "neither neutron or gamma spectrum Detector=";
      
      if( detector_attrib && detector_attrib->value_size() )
        msg << xml_value_str(detector_attrib);
      else
        msg << "NULL";
      
      msg << ", Type=";
      
      if( type_attrib && type_attrib->value_size() )
        msg << xml_value_str(type_attrib);
      else
        msg << "NULL";
      
      throw std::runtime_error( msg.str() );
    }//if( is_nuetron && is_gamma )
    
    return is_gamma;
  }//bool is_gamma_spectrum()

  
  //is_occupied(): throws exception on error
  bool is_occupied( const rapidxml::xml_node<char> *uccupied_node )
  {
    if( uccupied_node && uccupied_node->value_size() )
    {
      bool occupied = false;
      
      if( uccupied_node->value()[0] == '0' )
        occupied = false;
      else if( uccupied_node->value()[0] == '1' )
        occupied = true;
      else if( XML_VALUE_ICOMPARE(uccupied_node, "true") )
        occupied = true;
      else if( XML_VALUE_ICOMPARE(uccupied_node, "false") )
        occupied = false;
      else
      {
        stringstream msg;
        msg << SRC_LOCATION << "\n\tUnknown Occupied node value: '"
            << xml_value_str(uccupied_node) << "'";
        cerr << msg.str() << endl;
        throw std::runtime_error( msg.str() );
      }
      
      return occupied;
    }//if( uccupied_node is valid )
    
    throw std::runtime_error( "NULL <Occupied> node" );
  }//bool is_occupied( rapidxml::xml_node<char> *uccupied_node )

  
  const rapidxml::xml_attribute<char> *find_detector_attribute( const rapidxml::xml_node<char> *spectrum )
  {
    rapidxml::xml_attribute<char> *attribute = spectrum->first_attribute( "Detector", 8 );
    if( attribute )
      return attribute;
    
    using rapidxml::internal::compare;
    
    for( const rapidxml::xml_node<char> *node = spectrum->parent();
        (node && !XML_VALUE_ICOMPARE(node, "DetectorData"));
        node = node->parent() )
    {
      attribute = node->first_attribute( "Detector", 8 );
      if( attribute )
        return attribute;
    }//for( search parents incase it was out in the wrong node )
    
    //Avid N42 files contain a "Sensor" attribute in the <Spectrum> node that
    //  I think should be the detector name (unconfirmed for multiple detector
    //  systems as of 20150528)
    attribute = spectrum->first_attribute( "Sensor", 6 );
    
    return attribute;
  }//const rapidxml::xml_attribute<char> *find_detector_attribute( const rapidxml::xml_node<char> *spectrum )
  
  //speed_from_node(): returns speed in m/s; throws exception on error
  float speed_from_node( const rapidxml::xml_node<char> *speed_node )
  {
    if( !speed_node || !speed_node->value_size() )
      throw std::runtime_error( "speed_from_node(...): NULL <Speed> node" );
    
    float speed = 0.0f;
    if( !UtilityFunctions::parse_float( speed_node->value(), speed_node->value_size(), speed ) )
    {
      stringstream msg;
      msg << SRC_LOCATION << "\n\tUnable to convert '" << xml_value_str(speed_node)
      << "' to a float";
      cerr << msg.str() << endl;
      throw runtime_error( msg.str() );
    }//if( couldnt convert to float
    
    if( speed < 0.00000001f )
      return 0.0f;
    
    const rapidxml::xml_attribute<char> *unit_attrib = XML_FIRST_ATTRIB( speed_node, "Units" );
    if( !unit_attrib || !unit_attrib->value_size() )
    {
      cerr << SRC_LOCATION << "\n\t:Warning no units attribut avaliable in "
      << "<Speed> node, assuming m/s" << endl;
      return speed;
    }//if( no unit attribute" )
    
    string units = xml_value_str( unit_attrib );
    trim( units );
    to_lower( units );
    if( units == "mph" )
      return 0.44704f * speed;
    if( units == "m/s" )
      return speed;
    
    stringstream msg;
    msg << SRC_LOCATION << "\n\tUnknown speed units: '" << speed
    << "' - please fix";
    cerr << msg.str() << endl;
    throw std::runtime_error( msg.str() );
    
    return speed;
  }//float speed_from_node( rapidxml::xml_node<char> *speed_node )

}//namespace


#if(PERFORM_DEVELOPER_CHECKS)
void log_developer_error( const char *location, const char *error )
{
  static std::recursive_mutex s_dev_error_log_mutex;
  static ofstream s_dev_error_log( "developer_errors.log", ios::app | ios::out );
  
  std::unique_lock<std::recursive_mutex> loc( s_dev_error_log_mutex );
  
  boost::posix_time::ptime time = boost::posix_time::second_clock::local_time();
  
  const string timestr = boost::posix_time::to_iso_extended_string( time );
//  const string timestr = UtilityFunctions::to_iso_string( time );
  
  s_dev_error_log << timestr << ": " << location << endl << error << "\n\n" << endl;
  cerr << timestr << ": " << location << endl << error << "\n\n" << endl;
}//void log_developer_error( const char *location, const char *error )
#endif //#if(PERFORM_DEVELOPER_CHECKS)


//Analogous to Measurement::compare_by_sample_det_time; compares by
// sample_number, and then detector_number_, but NOT by start_time_
struct MeasurementInfoLessThan
{
  const int sample_number, detector_number;
  
  MeasurementInfoLessThan( int sample_num, int det_num )
  : sample_number( sample_num ), detector_number( det_num ) {}
  bool operator()( const MeasurementShrdPtr& lhs, const MeasurementShrdPtr &dummy )
  {
    assert( !dummy );
    if( !lhs )
      return false;
    
    if( lhs->sample_number() == sample_number )
      return (lhs->detector_number() < detector_number);
    return (lhs->sample_number() < sample_number);
  }//operator()
};//struct MeasurementInfoLessThan



struct MeasurementCalibInfo
{
  Measurement::EquationType equation_type;

  size_t nbin;
  vector<float> coefficients;
  vector< pair<float,float> > deviation_pairs_;
  ShrdConstFVecPtr original_binning;
  mutable ShrdConstFVecPtr binning;

  string calib_id; //optional
  
  MeasurementCalibInfo( MeasurementShrdPtr meas )
  {
    equation_type = meas->energy_calibration_model();
    nbin = meas->gamma_counts()->size();
    coefficients = meas->calibration_coeffs();
    deviation_pairs_ = meas->deviation_pairs();
    original_binning = meas->channel_energies();
    if( meas->channel_energies() && !meas->channel_energies()->empty())
      binning = meas->channel_energies();
  }//MeasurementCalibInfo constructor

  MeasurementCalibInfo()
  {
    nbin = 0;
    equation_type = Measurement::UnknownEquationType;
  }//MeasurementCalibInfo constructor

  void strip_end_zero_coeff()
  {
    for( int i = static_cast<int>(coefficients.size()) - 1;
         (i>=0) && (fabs(coefficients[i]) < 0.00000000001); --i )
      coefficients.pop_back();
  }//void strip_end_zero_coeff()


  void fill_binning() const
  {
    if( binning && !binning->empty())
      return;

    switch( equation_type )
    {
      case Measurement::Polynomial:
        binning = polynomial_binning( coefficients, nbin, deviation_pairs_ );
      break;

      case Measurement::FullRangeFraction:
        binning = fullrangefraction_binning( coefficients, nbin, deviation_pairs_ );
      break;

      case Measurement::LowerChannelEdge:
      {
        std::shared_ptr< vector<float> > energies
                                        = std::make_shared< vector<float> >();
        
        *energies  = coefficients;
        
        if( (nbin > 2) && (energies->size() == nbin) )
        {
          energies->push_back(2.0f*((*energies)[nbin-1]) - (*energies)[nbin-2]);
        }else if( energies->size() < nbin )
        {
          //Deal with error
          cerr << "fill_binning(): warning, not enough cahnnel energies to be"
                  " valid; expected at least " << nbin << ", but only got "
               << energies->size() << endl;
          energies.reset();
        }else if( energies->size() > (nbin+1) )
        {
          //Make it so
          cerr << "fill_binning(): warning, removing channel energy values!" << endl;
          energies->resize( nbin + 1 );
        }
          
        
        binning = energies;
        
        break;
      }//case Measurement::LowerChannelEdge:

      case Measurement::UnknownEquationType:
      break;
    }//switch( meas->energy_calibration_model_ )


    if( !binning )
    {
      cerr << SRC_LOCATION << "\n\tWarning, failed to sett binning" << endl;
      binning.reset(); //because I never make sure biining is valid...
    }

  }//void fill_binning()

  bool operator<( const MeasurementCalibInfo &rhs ) const
  {
    if( nbin != rhs.nbin )
      return (nbin < rhs.nbin);
    
    if( equation_type != rhs.equation_type )
      return (equation_type < rhs.equation_type);
    
    if( coefficients.size() != rhs.coefficients.size() )
      return (coefficients.size() < rhs.coefficients.size());
    
    
    for( size_t i = 0; i < coefficients.size(); ++i )
    {
      const float leftcoef = coefficients[i];
      const float rightcoef = rhs.coefficients[i];
      const float maxcoef = std::max( fabs(leftcoef), fabs(rightcoef) );
      if( fabs(leftcoef - rightcoef) > (1.0E-5 * maxcoef) )
        return coefficients[i] < rhs.coefficients[i];
    }//for( size_t i = 0; i < coefficients.size(); ++i )
    
    if( deviation_pairs_.size() != rhs.deviation_pairs_.size() )
      return (deviation_pairs_.size() < rhs.deviation_pairs_.size());
    
    for( size_t i = 0; i < deviation_pairs_.size(); ++i )
    {
      const pair<float,float> &l = deviation_pairs_[i];
      const pair<float,float> &r = rhs.deviation_pairs_[i];
      const float maxenergy = std::max( fabs(l.first), fabs(r.first) );
      
      if( fabs(l.first - r.first) > (1.0E-5 * maxenergy) )
        return l.first < r.first;
      
      const float maxdeviation = std::max( fabs(l.second), fabs(r.second) );
      
      if( fabs(l.second - r.second) > (1.0E-5 * maxdeviation) )
        return l.second < r.second;
    }//for( size_t i = 0; i < deviation_pairs_.size(); ++i )


    if( original_binning && !original_binning->empty() )
    {
      if( !rhs.original_binning || rhs.original_binning->empty() )
        return false;

      //a full lexicographical_compare is too much if the binnings are same, but
      //  different objects, so well just look at the objects address for all
      //  cases
      return (original_binning.get() < rhs.original_binning.get());

//      if( binning == rhs.binning )
//        return false;
//      return lexicographical_compare( binning->begin(), binning->end(),
//                                      rhs.binning->begin(), rhs.binning->end(), std::less<float>() );
    }else if( rhs.original_binning && !rhs.original_binning->empty() )
      return true;

    return false;
  }//bool operatr<(...)

  bool operator==( const MeasurementCalibInfo &rhs ) const
  {
    const bool rhsLt = operator<(rhs);
    const bool lhsLt = rhs.operator<(*this);
    return !lhsLt && !rhsLt;
  }
};//struct MeasurementCalibInfo

bool MeasurementInfo::calibration_is_valid( const Measurement::EquationType type,
                         const std::vector<float> &ineqn,
                         const std::vector< std::pair<float,float> > &devpairs,
                         size_t nbin )
{
  std::unique_ptr< std::vector<float> > frfeqn;
  const std::vector<float> *eqn = &ineqn;

  if( type == Measurement::Polynomial )
  {
    frfeqn.reset( new vector<float>() );
    *frfeqn = polynomial_coef_to_fullrangefraction( ineqn, nbin );
    eqn = frfeqn.get();
  }//if( type == Measurement::Polynomial )
  
  switch( type )
  {
    case Measurement::Polynomial:
    case Measurement::FullRangeFraction:
    {
      for( const float &val : (*eqn) )
      {
        if( IsInf(val) || IsNan(val) )
          return false;
      }
      
      const float nearend   = fullrangefraction_energy( float(nbin-2), *eqn, nbin, devpairs );
      const float end       = fullrangefraction_energy( float(nbin-1), *eqn, nbin, devpairs );
      const float begin     = fullrangefraction_energy( 0.0f,      *eqn, nbin, devpairs );
      const float nearbegin = fullrangefraction_energy( 1.0f,      *eqn, nbin, devpairs );
 
      const bool valid = !( (nearend >= end) || (begin >= nearbegin) );
      
#if( PERFORM_DEVELOPER_CHECKS )
      if( valid )
      {
        ShrdConstFVecPtr binning = fullrangefraction_binning( *eqn, nbin, devpairs );
        for( size_t i = 1; i < binning->size(); ++i )
        {
          if( (*binning)[i] <= (*binning)[i-1] )
          {
            char buffer[512];
            snprintf( buffer, sizeof(buffer),
                      "Energy of bin %i is lower or equal to the one before it",
                      int(i) );
            stringstream msg;
            msg << "Energy of bin " << i << " is lower or equal to the one before it"
                << " for coefficients {";
            for( size_t j = 0; j < eqn->size(); ++j )
              msg << (j ? ", " : "") << (*eqn)[j];
            msg << "} and deviation pairs {";
            
            for( size_t j = 0; j < devpairs.size(); ++j )
              msg << (j ? ", {" : "{") << devpairs[j].first << ", " << devpairs[j].second << "}";
            msg << "}";
            
            log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
          }
        }//for( size_t i = 1; i < binning->size(); ++i )
      }//if( valid )
#endif
      return valid;
    }//case Measurement::FullRangeFraction:
      
    case Measurement::LowerChannelEdge:
    {
      for( size_t i = 1; i < ineqn.size(); ++i )
        if( ineqn[i-1] >= ineqn[i] )
          return false;
      return (!ineqn.empty() && (ineqn.size()>=nbin));
    }//case Measurement::LowerChannelEdge:
      
    case Measurement::UnknownEquationType:
    break;
  }//switch( type )
  
  return false;
}//checkFullRangeFractionCoeffsValid



double gamma_integral( const std::shared_ptr<const Measurement> &hist,
                 const float minEnergy, const float maxEnergy )
{
  if( !hist )
    return 0.0;
  
  const double gamma_sum = hist->gamma_integral( minEnergy, maxEnergy );

#if( PERFORM_DEVELOPER_CHECKS )
  double check_sum = 0.0;
  
  const int lowBin = hist->FindFixBin( minEnergy );
  int highBin = hist->FindFixBin( maxEnergy );
  if( highBin > hist->GetNbinsX() )
    highBin = hist->GetNbinsX();
  
  if( lowBin == highBin )
  {
    const float binWidth = hist->GetBinWidth( lowBin );
    const double frac = (maxEnergy - minEnergy) / binWidth;
    check_sum = frac * hist->GetBinContent( lowBin );
  }else
  {
    if( lowBin > 0 )
    {
      const float lowerLowEdge = hist->GetBinLowEdge( lowBin );
      const float lowerBinWidth = hist->GetBinWidth( lowBin );
      const float lowerUpEdge = lowerLowEdge + lowerBinWidth;
    
      double fracLowBin = 1.0;
      if( minEnergy > lowerLowEdge )
        fracLowBin = (lowerUpEdge - minEnergy) / lowerBinWidth;
    
      check_sum += fracLowBin * hist->GetBinContent( lowBin );
    }//if( lowBin > 0 )
  
    if( highBin > 0 )
    {
      const float upperLowEdge = hist->GetBinLowEdge( highBin );
      const float upperBinWidth = hist->GetBinWidth( highBin );
      const float upperUpEdge = upperLowEdge + upperBinWidth;
    
      double fracUpBin = 1.0;
      if( maxEnergy < upperUpEdge )
        fracUpBin  = (maxEnergy - upperLowEdge) / upperBinWidth;
    
      check_sum += fracUpBin * hist->GetBinContent( highBin );
    }//if( highBin > 0 && lowBin!=highBin )
  
    for( int bin = lowBin + 1; bin < highBin; ++bin )
      check_sum += hist->GetBinContent( bin );
  }//if( lowBin == highBin ) / else
  
  
  if( check_sum != gamma_sum && !IsNan(check_sum) && !IsInf(check_sum)
      && lowBin != highBin )
  {
    char buffer[512];
    snprintf( buffer, sizeof(buffer),
              "Gamma Integral using new method varied from old methodology;"
              " %f (old) vs %f (new); lowBin=%i, highBin=%i" ,
              check_sum, gamma_sum, lowBin, highBin );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }//if( check_sum != gamma_sum )
#endif  //#if( PERFORM_DEVELOPER_CHECKS )

  return gamma_sum;
}//float integral(...)



double Measurement::gamma_integral( float lowerx, float upperx ) const
{
  double sum = 0.0;
  
  if( !channel_energies_ || !gamma_counts_ || channel_energies_->size() < 2
      || gamma_counts_->size() < 2 )
    return sum;

  const vector<float> &x = *channel_energies_;
  const vector<float> &y = *gamma_counts_;
  
  const size_t nchannel = x.size();
  const float maxX = 2.0f*x[nchannel-1] - x[nchannel-2];
  
  lowerx = std::min( maxX, std::max( lowerx, x[0] ) );
  upperx = std::max( x[0], std::min( upperx, maxX ) );
  
  if( lowerx == upperx )
    return sum;
  
  if( lowerx > upperx )
    std::swap( lowerx, upperx );
 
  //need to account for edgecase of incase x.size() != y.size()
  const size_t maxchannel = gamma_counts_->size() - 1;
  const size_t lowerChannel = min( find_gamma_channel( lowerx ), maxchannel );
  const size_t upperChannel = min( find_gamma_channel( upperx ), maxchannel );
  
  const float lowerLowEdge = x[lowerChannel];
  const float lowerBinWidth = (lowerChannel < (nchannel-1))
                                      ? x[lowerChannel+1] - x[lowerChannel]
                                      : x[lowerChannel]   - x[lowerChannel-1];
  const float lowerUpEdge = lowerLowEdge + lowerBinWidth;
  
  if( lowerChannel == upperChannel )
  {
    const double frac = (upperx - lowerx) / lowerBinWidth;
    return frac * y[lowerChannel];
  }
  
  const double fracLowBin = (lowerUpEdge - lowerx) / lowerBinWidth;
  sum += fracLowBin * y[lowerChannel];
  
  const float upperLowEdge = x[upperChannel];
  const float upperBinWidth = (upperChannel < (nchannel-1))
                                      ? x[upperChannel+1] - x[upperChannel]
                                      : x[upperChannel]   - x[upperChannel-1];
//  const float upperUpEdge = upperLowEdge + upperBinWidth;
  const double fracUpBin  = (upperx - upperLowEdge) / upperBinWidth;
  sum += fracUpBin * y[upperChannel];

  for( size_t channel = lowerChannel + 1; channel < upperChannel; ++channel )
    sum += y[channel];
  
  return sum;
}//double gamma_integral( const float lowEnergy, const float upEnergy ) const;


double Measurement::gamma_channels_sum( size_t startbin, size_t endbin ) const
{
  double sum = 0.0;
  if( !gamma_counts_ )
    return sum;
  
  const size_t nchannels = gamma_counts_->size();
  
  if( startbin >= nchannels )
    return sum;
  
  endbin = std::min( endbin, nchannels-1 );
  
  if( startbin > endbin )
    std::swap( startbin, endbin );
  
  for( size_t channel = startbin; channel <= endbin; ++channel )
    sum += (*gamma_counts_)[channel];
  
  return sum;
}//double gamma_channels_sum( size_t startbin, size_t endbin ) const;



const char *descriptionText( const SpectrumType type )
{
  switch( type )
  {
    case kForeground:       return "Foreground";
    case kSecondForeground: return "Secondary";
    case kBackground:       return "Background";
  }//switch( type )
  
  return "";
}//const char *descriptionText( const SpectrumType type )


const char *suggestedNameEnding( const SaveSpectrumAsType type )
{
  switch( type )
  {
    case kTxtSpectrumFile:                return "txt";
    case kCsvSpectrumFile:                return "csv";
    case kPcfSpectrumFile:                return "pcf";
    case kXmlSpectrumFile:                return "n42";
    case k2012N42SpectrumFile:            return "n42";
    case kChnSpectrumFile:                return "chn";
    case kBinaryIntSpcSpectrumFile:       return "spc";
    case kBinaryFloatSpcSpectrumFile:     return "spc";
    case kAsciiSpcSpectrumFile:           return "spc";
    case kExploraniumGr130v0SpectrumFile: return "dat";
    case kExploraniumGr135v2SpectrumFile: return "dat";
    case kIaeaSpeSpectrumFile:            return "spe";
#if( ENABLE_D3_CHART_EXPORTING )
    case kD3HtmlSpectrumFile:             return "html";
#endif
    case kNumSaveSpectrumAsType:          break;
  }//switch( m_format )
  
  return "";
}//const char *suggestedNameEnding( const SaveSpectrumAsType type )


SpectrumType spectrumTypeFromDescription( const char *descrip )
{
  if( strcmp(descrip,descriptionText(kForeground)) == 0 )
    return kForeground;
  if( strcmp(descrip,descriptionText(kSecondForeground)) == 0 )
    return kSecondForeground;
  if( strcmp(descrip,descriptionText(kBackground)) == 0 )
    return kBackground;

  throw runtime_error( "spectrumTypeFromDescription(...): invalid descrip: "
                        + string(descrip) );
  
  return kForeground;
}//SpectrumType spectrumTypeFromDescription( const char *descrip )


const char *descriptionText( const SaveSpectrumAsType type )
{
  switch( type )
  {
    case kTxtSpectrumFile:                return "TXT";
    case kCsvSpectrumFile:                return "CSV";
    case kPcfSpectrumFile:                return "PCF";
    case kXmlSpectrumFile:                return "2006 N42";
    case k2012N42SpectrumFile:            return "2012 N42";
    case kChnSpectrumFile:                return "CHN";
    case kBinaryIntSpcSpectrumFile:       return "Integer SPC";
    case kBinaryFloatSpcSpectrumFile:     return "Float SPC";
    case kAsciiSpcSpectrumFile:           return "ASCII SPC";
    case kExploraniumGr130v0SpectrumFile: return "GR130 DAT";
    case kExploraniumGr135v2SpectrumFile: return "GR135v2 DAT";
    case kIaeaSpeSpectrumFile:            return "IAEA SPE";
#if( ENABLE_D3_CHART_EXPORTING )
    case kD3HtmlSpectrumFile:             return "D3 HTML";
#endif
    case kNumSaveSpectrumAsType:          return "";
  }
  return "";
}//const char *descriptionText( const SaveSpectrumAsType type )


int sample_num_from_remark( const std::string &remark )
{
  size_t pos = remark.find( "Survey" );

  if( pos == string::npos )
  {
//    cerr << "Remark '" << remark << "'' didnt contain a sample num" << endl;
    return -1;
  }

  pos = remark.find_first_not_of( " \t=", pos+6 );
  if( pos == string::npos )
  {
    cerr << "Remark '" << remark << "'' didnt have a integer sample num" << endl;
    return -1;
  }

  int num = -1;
  if( !(stringstream(remark.c_str()+pos) >> num) )
  {
     cerr << "sample_num_from_remark(...): Error converting '"
          << remark.c_str()+pos << "' to float" << endl;
    return -1;
  }//if( cant convert result to int )

  return num;
}//int sample_num_from_remark( const std::string &remark )

//Takes a line like "Speed = 5 mph" and returns the speed in m/s.
//  Returns 0 upon failure.
float speed_from_remark( std::string remark )
{
  to_lower( remark );
  size_t pos = remark.find( "speed" );

  if( pos == string::npos )
    return 0.0;

  pos = remark.find_first_not_of( "= \t", pos+5 );
  if( pos == string::npos )
    return 0.0;

  const string speedstr = remark.substr( pos );

  float speed = 0.0;
  if( !toFloat( speedstr, speed) )
  {
    cerr << "speed_from_remark(...): couldn conver to number: '"
         << speedstr << "'" << endl;
    return 0.0;
  }//if( !(stringstream(speedstr) >> speed) )


  for( size_t i = 0; i < speedstr.size(); ++i )
  {
    if( (!isdigit(speedstr[i])) && (speedstr[i]!=' ') && (speedstr[i]!='\t') )
    {
      float convertion = 0.0f;

      const string unitstr = speedstr.substr( i );
      const size_t unitstrlen = unitstr.size();

      if( unitstrlen>=3 && unitstr.substr(0,3) == "m/s" )
        convertion = 1.0f;
      else if( unitstrlen>=3 && unitstr.substr(0,3) == "mph" )
        convertion = 0.44704f;
      else
        cerr << "speed_from_remark(...): Unknown speed unit: '"
             << unitstrlen << "'" << endl;

      return convertion*speed;
    }//if( we found the start of the units )
  }//for( size_t i = 0; i < speedstr.size(); ++i )

  return 0.0;
}//float speed_from_remark( const std::string &remark )


//Looks for GADRAS style detector names in remarks, or something from the N42
//  conventions of 'Aa1', 'Aa2', etc.
//  Returns empty string on failure.
std::string detector_name_from_remark( const std::string &remark )
{
  //Check for the Gadras convention similar to "Det=Aa1"
  if( UtilityFunctions::icontains(remark, "det") )
  {
    //Could use a regex here... maybe someday I'll upgrade
    string remarkcopy = remark;
    
    string remarkcopylowercase = remarkcopy;
    UtilityFunctions::to_lower( remarkcopylowercase );
    
    size_t pos = remarkcopylowercase.find( "det" );
    if( pos != string::npos )
    {
      remarkcopy = remarkcopy.substr(pos);
      pos = remarkcopy.find_first_of( " =" );
      if( pos != string::npos )
      {
        string det_identifier = remarkcopy.substr(0,pos);
        UtilityFunctions::to_lower( det_identifier );
        UtilityFunctions::trim( det_identifier ); //I dont htink this is necassarry
        if( det_identifier=="det" || det_identifier=="detector" )
        {
          remarkcopy = remarkcopy.substr(pos);  //get rid of the "det="
          while( remarkcopy.length() && (remarkcopy[0]==' ' || remarkcopy[0]=='=') )
            remarkcopy = remarkcopy.substr(1);
          pos = remarkcopy.find_first_of( ' ' );
          if( pos != string::npos && pos > 1 && pos < 35 ) //1 and 25 arbitrary
            return remarkcopy.substr(0,pos);
        }
      }
    }
    
  }//if( UtilityFunctions::icontains(remark, "det") )
  
  
  vector<string> split_contents;
  split( split_contents, remark, ", \t\r\n" );

  for( const string &field : split_contents )
  {
    if( (field.length() < 3) ||  !isdigit(field[field.size()-1])
        || (field.length() > 4) ||  (field[1] != 'a') )
      continue;

    if( field[0]!='A' && field[0]!='B' && field[0]!='C' && field[0]!='D' )
      continue;

    return field;
  }//for( size_t i = 0; i < split_contents.size(); ++i )

  return "";
}//std::string detector_name_from_remark( const std::string &remark )



float time_duration_in_seconds( const std::string &duration )
{
  return time_duration_in_seconds( duration.c_str(), duration.size() );
}


float time_duration_in_seconds( const char *duration_str, const size_t len )
{
  float durration = 0.0f;

  if( !duration_str || !len )
    return durration;

  const char *orig = duration_str;
  const char *end = duration_str + len;
  
  while( duration_str < end )
  {
    while( (duration_str < end) && !isdigit(*duration_str) )
      ++duration_str;

    if( duration_str >= end )
      break;

    const char *num_start = duration_str;
    while( (duration_str < end) && (isdigit(*duration_str) || duration_str[0]=='.') )
      ++duration_str;
    const char *num_end = duration_str;

    float unit = 1.0;
    char unitchar = (num_end<end) ? (*num_end) : 's';

    if( unitchar=='s' || unitchar=='S' )
      unit = 1.0;
    else if( unitchar=='m' || unitchar=='M' )
      unit = 60.0;
    else if( unitchar=='h' || unitchar=='H' )
      unit = 3600.0;

    float num_val = 0.0;
    const string numberstr( num_start, num_end );
    if( toFloat( numberstr, num_val ) )
    {
      durration += num_val * unit;
    }else
    {
      cerr << "Error parsing live time from '" << string(orig,orig+len) << "'" << endl;
      return durration;
    }
  }//while( *duration_str )

  return durration;
}//float time_duration_in_seconds( const char *duration_str )


float dose_units_usvPerH( const char *str, const size_t str_length )
{
  if( !str )
    return 0.0f;
  
  if( icontains( str, str_length, "uSv", 3 )
     || icontains( str, str_length, "\xc2\xb5Sv", 4) )
    return 1.0f;
  
  //One sievert equals 100 rem.
  if( icontains( str, str_length, "&#xB5;Rem/h", 11 ) ) //micro
    return 0.01f;

  return 0.0f;
}//float dose_units_usvPerH( const char *str )


bool likely_not_spec_file( const std::string &fullpath )
{
  const std::string extension = UtilityFunctions::file_extension( fullpath );
  const std::string filename = UtilityFunctions::filename( fullpath );
  const std::string dir = UtilityFunctions::parent_path( fullpath );
  
  const char * const bad_exts[] = { ".jpg", ".jpeg", ".zip", ".docx", ".png",
      ".pdf", ".html", ".xtk3d", ".xtk", ".doc", /*".txt",*/ ".An1", ".rpt",
      ".ufo", ".bmp", ".IIF", ".xls", ".xlsx", ".ds_store", ".kmz", ".msg",
      ".exe", ".jp2", ".wmv", /*".gam",*/ ".pptx", ".htm", ".ppt", ".mht",
      ".ldb", ".lis", ".zep", ".ana", ".eft", ".clb", ".lib", ".wav", ".gif",
      ".wmf", ".phd", ".log", ".vi", ".incident", ".tiff", ".cab", ".ANS",
      ".menc", ".tif", ".psd", ".mdb", ".drill", ".lnk", ".mov", ".rtf", ".shx",
      ".dbf", ".prj", ".sbn", ".shb", ".inp1", ".bat", ".xps", ".svy", ".ini",
      ".2", ".mp4", ".sql", ".gz", ".url", ".zipx", ".001", ".002", ".003" };
  const size_t num_bad_exts = sizeof(bad_exts) / sizeof(bad_exts[0]);
  
  for( size_t i = 0; i < num_bad_exts; ++i )
    if( UtilityFunctions::iequals(extension, bad_exts[i]) )
      return true;
    
  if( //(filename.find("_AA_")!=std::string::npos && UtilityFunctions::iequals(extension, ".n42") )
     //|| (filename.find("_AN_")!=std::string::npos && UtilityFunctions::iequals(extension, ".n42") )
     filename.find("Neutron.n42") != std::string::npos
     || filename.find(".xml.XML") != std::string::npos
     || filename.find("results.xml") != std::string::npos
     || filename.find("Rebin.dat") != std::string::npos
     || filename.find("Detector.dat") != std::string::npos
     || UtilityFunctions::iends_with( fullpath, ".html")
     || filename.empty()
     || filename[0] == '.'
     || extension.empty()
     || UtilityFunctions::file_size(fullpath) < 100
     )
    return true;
  
  return false;
}//bool likely_not_spec_file( const boost::filesystem::path &file )



void setAnalysisInformation( const rapidxml::xml_node<char> *analysis_node,
                             std::shared_ptr<DetectorAnalysis> analysis )
{
  if( !analysis_node || !analysis )
    return;
 
//  boost::posix_time::ptime ana_start_time;
//  const rapidxml::xml_node<char> *AnalysisStartDateTime = XML_FIRST_NODE( analysis_node, "AnalysisStartDateTime" );
//  if( AnalysisStartDateTime )
//    ana_start_time = time_from_string( xml_value_str(AnalysisStartDateTime).c_str() );
  
  const rapidxml::xml_node<char> *nuc_ana_node = XML_FIRST_NODE( analysis_node, "NuclideAnalysis" );
  if( !nuc_ana_node )
    nuc_ana_node = XML_FIRST_NODE(analysis_node, "NuclideAnalysisResults");
  if( !nuc_ana_node )
  {
    //Special case for RadSeeker
    if( XML_FIRST_NODE(analysis_node, "Nuclide") )
      nuc_ana_node = analysis_node;
    
    //<Algorithm>
    //<AlgorithmVendor>SomeBrand</AlgorithmVendor>
    //<AlgorithmName>SomeName</AlgorithmName>
    //<AlgorithmVersion>2.9.3.701</AlgorithmVersion>
    //<FirmwareVersion>Not Set</FirmwareVersion>
    //<SoftwareVersion>2.24.144</SoftwareVersion>
    //<Parameters ParameterName="IsotopeLibrary" ParameterValue="11215" />
    //<Parameters ParameterName="DSP" ParameterValue="1.1.3.2" />
    //</Algorithm>
  }
  
//  if( !nuc_ana_node )
//    return;
  
  for( const rapidxml::xml_node<char> *remark_node = XML_FIRST_NODE(analysis_node, "Remark");
      remark_node;
      remark_node = XML_NEXT_TWIN(remark_node) )
  {
    string remark = xml_value_str( remark_node );
    UtilityFunctions::trim(remark);
    if(!remark.empty())
      analysis->remarks_.push_back( remark );
  }
  
  const rapidxml::xml_node<char> *algo_info_node = XML_FIRST_NODE_CHECKED( analysis_node->parent(), "Algorithm" );
  
  const rapidxml::xml_node<char> *algo_name_node = XML_FIRST_NODE( analysis_node, "AnalysisAlgorithmName" );
  if( !algo_name_node && algo_info_node )
    algo_name_node = XML_FIRST_NODE( algo_info_node, "AlgorithmName" );
  analysis->algorithm_name_ = xml_value_str(algo_name_node);
  
  if( analysis->algorithm_name_.empty() && nuc_ana_node )
  {
    const rapidxml::xml_attribute<char> *att = XML_FIRST_ATTRIB( nuc_ana_node, "AlgorithmName" );
    analysis->algorithm_name_ = xml_value_str(att);
  }
    
  const rapidxml::xml_node<char> *algo_version_node = XML_FIRST_NODE( analysis_node, "AnalysisAlgorithmVersion" );
  if( !algo_version_node && algo_info_node )
    algo_version_node = XML_FIRST_NODE( algo_info_node, "AlgorithmVersion" );
  analysis->algorithm_version_ = xml_value_str(algo_version_node);
  
  if( analysis->algorithm_version_.empty() && nuc_ana_node )
  {
    const rapidxml::xml_attribute<char> *algorithm_version_att = nuc_ana_node->first_attribute( "AlgorithmVersion", 16 );
    analysis->algorithm_version_ = xml_value_str(algorithm_version_att);
  }
  
  trim( analysis->algorithm_version_ ); //RadEagle has a new line in it...
  
  const rapidxml::xml_node<char> *algo_creator_node = XML_FIRST_NODE( analysis_node, "AnalysisAlgorithmCreatorName" );
  if( !algo_creator_node && algo_info_node )
    algo_creator_node = XML_FIRST_NODE( algo_info_node, "AlgorithmVendor" );
  analysis->algorithm_creator_ = xml_value_str(algo_creator_node);
  
  const rapidxml::xml_node<char> *algo_desc_node = XML_FIRST_NODE( analysis_node, "AnalysisAlgorithmDescription" );
  analysis->algorithm_description_ = xml_value_str(algo_desc_node);
  
  
  if( analysis->algorithm_description_.empty() && algo_info_node )
  {
    //RadSeeker specific
    const rapidxml::xml_node<char> *FirmwareVersion = XML_FIRST_NODE( algo_info_node, "FirmwareVersion" );
    const rapidxml::xml_node<char> *SoftwareVersion = XML_FIRST_NODE( algo_info_node, "SoftwareVersion" );
    
    string desc;
    if( FirmwareVersion && FirmwareVersion->value_size() )
      desc += "FirmwareVersion: " + xml_value_str(FirmwareVersion);
    if( SoftwareVersion && SoftwareVersion->value_size() )
      desc += string(desc.empty()?"":", ") + "SoftwareVersion: " + xml_value_str(SoftwareVersion);
    
    for( const rapidxml::xml_node<char> *Parameters = XML_FIRST_NODE( algo_info_node, "Parameters" );
        Parameters;
        Parameters = XML_NEXT_TWIN(Parameters) )
    {
      const rapidxml::xml_attribute<char> *name = XML_FIRST_ATTRIB(Parameters, "ParameterName");
      const rapidxml::xml_attribute<char> *value XML_FIRST_ATTRIB(Parameters, "ParameterValue");
      if( name && value && name->value_size() && value->value_size() )
        desc += string(desc.empty()?"":", ") + xml_value_str(name) + ": " + xml_value_str(value);
    }//for( loop over parameters )
    
    analysis->algorithm_description_.swap( desc );
  }//if( analysis->algorithm_description_.empty() && algo_info_node )
  
  
  const rapidxml::xml_node<char> *result_desc_node = XML_FIRST_NODE( analysis_node, "AnalysisResultDescription" );
  if( !result_desc_node )
    result_desc_node = XML_FIRST_NODE( analysis_node, "ThreatDescription" );
  analysis->algorithm_result_description_ = xml_value_str(result_desc_node);
  
  
  //Should loop over NuclideAnalysis/ NuclideAnalysisResults ...
  for( ; nuc_ana_node; nuc_ana_node = XML_NEXT_TWIN(nuc_ana_node) )
  {
    for( const rapidxml::xml_node<char> *nuclide_node = XML_FIRST_NODE( nuc_ana_node, "Nuclide" );
         nuclide_node;
         nuclide_node = XML_NEXT_TWIN(nuclide_node) )
    {
      const rapidxml::xml_node<char> *remark_node = nuclide_node->first_node( "Remark", 6 );
      const rapidxml::xml_node<char> *nuclide_name_node = nuclide_node->first_node( "NuclideName", 11 );
      const rapidxml::xml_node<char> *nuclide_type_node = nuclide_node->first_node( "NuclideType", 11 );
      const rapidxml::xml_node<char> *confidence_node = nuclide_node->first_node( "NuclideIDConfidenceIndication", 29 );
      if( !confidence_node )
        confidence_node = XML_FIRST_NODE(nuclide_node, "NuclideIDConfidence");  //RadSeeker
      const rapidxml::xml_node<char> *id_desc_node = nuclide_node->first_node( "NuclideIDConfidenceDescription", 30 );
      const rapidxml::xml_node<char> *position_node = nuclide_node->first_node( "SourcePosition", 14 );
      const rapidxml::xml_node<char> *id_indicator_node = nuclide_node->first_node( "NuclideIdentifiedIndicator", 26 ); //says 'true' or 'false', seen in refZ077SD6DVZ
      const rapidxml::xml_node<char> *confidence_value_node = nuclide_node->first_node( "NuclideIDConfidenceValue", 24 ); //seen in refJHQO7X3XFQ (FLIR radHUNTER UL-LGH)
      const rapidxml::xml_node<char> *catagory_desc_node = nuclide_node->first_node( "NuclideCategoryDescription", 26 );  //seen in refQ7M2PV6MVJ
    
      //Some files list true/false if a nuclide in its "trigger" list is present
      //  If a nuclide is not there, lets not include it in results.
      if( id_indicator_node && XML_VALUE_ICOMPARE(id_indicator_node, "false") )
      {
        //Could also check if NuclideIDConfidenceValue is there and if so make sure its value is equal to "0" 
        continue;
      }
    
      DetectorAnalysisResult result;
      result.remark_ = xml_value_str(remark_node);
      result.nuclide_ = xml_value_str(nuclide_name_node);
      
      if( nuclide_type_node && nuclide_type_node->value_size() )
        result.nuclide_type_ = xml_value_str(nuclide_type_node);
      else if( catagory_desc_node && catagory_desc_node->value_size() )
        result.nuclide_type_ = xml_value_str(catagory_desc_node);
    
      if( confidence_node && confidence_node->value_size() )
        result.id_confidence_ = xml_value_str(confidence_node);
      else if( confidence_value_node && confidence_value_node->value_size() )
        result.id_confidence_ = xml_value_str(confidence_value_node);
    
      const rapidxml::xml_node<char> *nuc_activity_node = nuclide_node->first_node( "NuclideActivityValue", 20 );
      if( nuc_activity_node && nuc_activity_node->value_size() )
      {
       const rapidxml::xml_attribute<char> *activity_units_att = nuc_activity_node->first_attribute( "units", 5 );

        double activity_units = 1.0e+3;
        if( activity_units_att && activity_units_att->value_size() )
        {
          //This is a mini implimentation of PhysicalUnits::stringToActivity(...)
          string letters = xml_value_str(activity_units_att);
        
          if( istarts_with(letters, "n" ) )
            activity_units = 1.0E-9;
          else if( istarts_with(letters, "u" )
                  || istarts_with(letters, "micro" )
                  || istarts_with(letters, "\xc2\xb5" ) )
            activity_units = 1.0E-6;
          else if( starts_with(letters, "m" )
                  || istarts_with(letters, "milli" ) )
            activity_units = 1.0E-3;
          else if( istarts_with(letters, "b" )
                  || istarts_with(letters, "c" ) )
            activity_units = 1.0;
          else if( istarts_with(letters, "k" ) )
            activity_units = 1.0E+3;
          else if( starts_with(letters, "M" )
                  || istarts_with(letters, "mega" ) )
            activity_units = 1.0E+6;
          else
            activity_units = 0.0;
        
          const bool hasb = icontains( letters, "b" );
          const bool hasc = icontains( letters, "c" );
        
          if( hasc && !hasb )
            activity_units *= 3.7e+10;
          else
            activity_units = 1.0e+3; //0.0;
        }//if( activity_units_att && activity_units_att->value() )
      
        xml_value_to_flt(nuc_activity_node, result.activity_);
        result.activity_ *= activity_units;
      }//if( nuc_activity_node && nuc_activity_node->value_size() )
    
      if( position_node )
      {
        const rapidxml::xml_node<char> *location_node = position_node->first_node( "RelativeLocation", 16 );
        if( location_node )
        {
          const rapidxml::xml_node<char> *dist_node = location_node->first_node( "DistanceValue", 13 );
          //Could check 'units' attribute, but only valid value is "m"
          if( xml_value_to_flt(dist_node, result.distance_ ) )
            result.distance_ *= 1000.0;
        }//if( position_node )
      }//if( position_node )

    
      const rapidxml::xml_node<char> *extention_node = nuclide_node->first_node( "NuclideExtension", 16 );
      if( extention_node )
      {
        const rapidxml::xml_node<char> *SampleRealTime = extention_node->first_node( "SampleRealTime", 14 );
        const rapidxml::xml_node<char> *Detector = extention_node->first_node( "Detector", 8 );
      
        if( SampleRealTime && SampleRealTime->value_size() )
          result.real_time_ = time_duration_in_seconds( SampleRealTime->value(), SampleRealTime->value_size() );

        result.detector_ = xml_value_str(Detector);
      }//if( extention_node )
    

      //Some N42 files include analysis results for all nuclides, even if they
      //  arent present, let not keep those.
      if( id_desc_node && XML_VALUE_ICOMPARE(id_desc_node, "Not present") )
      {
        continue;
      }
    
      if( result.remark_.empty() && result.nuclide_.empty()
          && result.nuclide_type_.empty()
          && result.dose_rate_ <= 0.0 )
        continue;
    
//    cerr << "result.remark_=" << result.remark_
//         << ", result.nuclide_=" << result.nuclide_
//         << ", result.nuclide_type_=" << result.nuclide_type_
//         << ", result.id_confidence_=" << result.id_confidence_ << endl;
    
      analysis->results_.push_back( result );
    }//For( loop over nuclide analysis )
  }//for( ; nuc_ana_node; nuc_ana_node = XML_NEXT_TWIN(nuc_ana_node) )

  for( const rapidxml::xml_node<char> *dose_node = analysis_node->first_node( "DoseAnalysisResults", 19 );
      dose_node;
      dose_node = XML_NEXT_TWIN(dose_node) )
  {
    const rapidxml::xml_node<char> *remark_node = dose_node->first_node( "Remark", 6 );
    //const rapidxml::xml_node<char> *start_time_node = dose_node->first_node( "StartTime", 9 );
    const rapidxml::xml_node<char> *avrg_dose_node = dose_node->first_node( "AverageDoseRateValue", 20 );
    const rapidxml::xml_node<char> *total_dose_node = dose_node->first_node( "TotalDoseValue", 14 );
    const rapidxml::xml_node<char> *position_node = dose_node->first_node( "SourcePosition", 14 );
    
    DetectorAnalysisResult result;
    result.remark_ = xml_value_str(remark_node);

    //if( start_time_node && start_time_node->value() )
      //result.start_time_ = time_from_string( xml_value_str(start_time_node).c_str() );
    
    //assuming in units of micro-sievert per hour (could check 'units' attribute
    //  but only valid value is uSv/h anyway
    xml_value_to_flt(avrg_dose_node, result.dose_rate_ );
    
    if( total_dose_node )
    {
      float total_dose;
      xml_value_to_flt(total_dose_node, total_dose );
      if( (result.dose_rate_> 0.0f) && (total_dose > 0.0f) )
      {
        result.real_time_ = total_dose / result.dose_rate_;
      }else if( total_dose > 0.0f )
      {
        //Uhhhg, this isnt correct, but maybye better than nothing?
        result.dose_rate_ = total_dose;
      }
    }
    
    if( position_node )
    {
      const rapidxml::xml_node<char> *location_node = position_node->first_node( "RelativeLocation", 16 );
      if( location_node )
      {
        const rapidxml::xml_node<char> *dist_node = location_node->first_node( "DistanceValue", 13 );
        
        //Could check 'units' attribute, but only val;id value is "m"
        if( xml_value_to_flt( dist_node, result.distance_ ) )
          result.distance_ *= 1000.0;
      }//if( position_node )
    }//if( position_node )
    
    analysis->results_.push_back( result );
  }//for( loop over DoseAnalysisResults nodes )
}//void setAnalysisInformation(...)



const std::string &detectorTypeToString( const DetectorType type )
{
  static const string GR135Detector             = "GR135";
  static const string IdentiFinderDetector      = "IdentiFINDER";
  static const string IdentiFinderNGDetector    = "IdentiFINDER-NG";
  static const string IdentiFinderLaBr3Detector = "IdentiFINDER-LaBr3";
  static const string DetectiveDetector         = "Detective";
  static const string DetectiveExDetector       = "Detective-EX";
  static const string DetectiveEx100Detector    = "Detective-EX100";
  static const string OrtecIDMPortalDetector    = "Ortec IDM Portal";
  static const string SAIC8Detector             = "SAIC8";
  static const string Falcon5kDetector          = "Falcon 5000";
  static const string UnknownDetector           = "Unknown";
  static const string MicroDetectiveDetector    = "MicroDetective";
  static const string MicroRaiderDetector       = "MicroRaider";
  static const string Sam940Detector            = "SAM940";
  static const string Sam940Labr3Detector       = "SAM940LaBr3";
  static const string Sam945Detector            = "SAM945";
  static const string Rsi701Detector            = "RS-701";
  static const string Rsi705Detector            = "RS-705";
  static const string RadHunterNaIDetector      = "RadHunterNaI";
  static const string RadHunterLaBr3Detector    = "RadHunterLaBr3";
  static const string AvidRsiDetector           = "RSI-Unspecified";
  static const string RadEagleNaiDetector       = "RadEagle NaI 3x1";
  static const string RadEagleCeBr2InDetector   = "RadEagle CeBr3 2x1";
  static const string RadEagleCeBr3InDetector   = "RadEagle CeBr3 3x0.8";
  static const string RadEagleLaBrDetector      = "RadEagle LaBr3 2x1";
  
  
//  GN3, InSpector 1000 LaBr3, Pager-S, SAM-Eagle-LaBr, GR130, SAM-Eagle-NaI-3x3
//  InSpector 1000 NaI, RadPack, SpiR-ID LaBr3, Interceptor, Radseeker, SpiR-ID NaI
//  GR135Plus, LRM, Raider, HRM, LaBr3PNNL, Transpec, Falcon 5000, Ranger
//  MicroDetective, FieldSpec, IdentiFINDER-NG, SAM-935, NaI 3x3, SAM-Eagle-LaBr3

  switch( type )
  {
    case kGR135Detector:
      return GR135Detector;
    case kIdentiFinderNGDetector:
      return IdentiFinderNGDetector;
//  kIdentiFinderNGDetector,   //I dont have any examples of this
    case kIdentiFinderDetector:
      return IdentiFinderDetector;
    case kIdentiFinderLaBr3Detector:
      return IdentiFinderLaBr3Detector;
    case kDetectiveDetector:
      return DetectiveDetector;
    case kDetectiveExDetector:
      return DetectiveExDetector;
    case kDetectiveEx100Detector:
      return DetectiveEx100Detector;
    case kOrtecIDMPortalDetector:
      return OrtecIDMPortalDetector;
    case kSAIC8Detector:
      return SAIC8Detector;
    case kFalcon5000:
      return Falcon5kDetector;
    case kUnknownDetector:
      return UnknownDetector;
    case kMicroDetectiveDetector:
      return MicroDetectiveDetector;
    case kMicroRaiderDetector:
      return MicroRaiderDetector;
    case kSam940:
      return Sam940Detector;
    case kSam945:
      return Sam945Detector;
    case kSam940LaBr3:
      return Sam940Labr3Detector;
    case kRsi701:
      return Rsi701Detector;
    case kRadHunterNaI:
      return RadHunterNaIDetector;
    case kRadHunterLaBr3:
      return RadHunterLaBr3Detector;
    case kRsi705:
      return Rsi705Detector;
    case kAvidRsi:
      return AvidRsiDetector;
    case kOrtecRadEagleNai:
      return RadEagleNaiDetector;
    case kOrtecRadEagleCeBr2Inch:
      return RadEagleCeBr2InDetector;
    case kOrtecRadEagleCeBr3Inch:
      return RadEagleCeBr3InDetector;
    case kOrtecRadEagleLaBr:
      return RadEagleLaBrDetector;
  }//switch( type )

  return UnknownDetector;
}//const std::string &detectorTypeToString( const DetectorType type )


Measurement::Measurement()
{
  reset();
}//Measurement()


bool Measurement::valid_longitude( const double longitude )
{
  return (fabs(static_cast<double>(longitude))<=180.0 && !IsInf(longitude) );
}

bool Measurement::valid_latitude( const double latitude )
{
  return (fabs(static_cast<double>(latitude))<=90.0 && !IsInf(latitude) );
}


size_t Measurement::memmorysize() const
{
  size_t size = sizeof(*this);

  //now we need to take care of all the non-simple objects
  size += detector_name_.capacity()*sizeof(string::value_type);
  size += detector_type_.capacity()*sizeof(string::value_type);

  for( const string &r : remarks_ )
    size += r.capacity()*sizeof(string::value_type);

  size += title_.capacity()*sizeof(string::value_type);

  if( gamma_counts_ )
    size += sizeof(*(gamma_counts_.get())) + gamma_counts_->capacity()*sizeof(float);
  if( channel_energies_ )
    size += sizeof(*(channel_energies_.get())) + channel_energies_->capacity()*sizeof(float);

  size += calibration_coeffs_.capacity()*sizeof(float);
  size += neutron_counts_.capacity()*sizeof(float);

  size += deviation_pairs_.capacity()*sizeof(std::pair<float,float>);

  return size;
}//size_t Measurement::memmorysize() const


bool Measurement::compare_by_sample_det_time( const MeasurementConstShrdPtr &lhs,
                           const MeasurementConstShrdPtr &rhs )
{
  if( !lhs || !rhs )
    return false;

  if( lhs->sample_number_ != rhs->sample_number_ )
    return (lhs->sample_number_ < rhs->sample_number_);

  if( lhs->detector_number_ != rhs->detector_number_ )
    return (lhs->detector_number_ < rhs->detector_number_);

  if( lhs->start_time_ != rhs->start_time_ )
    return (lhs->start_time_ < rhs->start_time_);
  
  return (lhs->source_type() < rhs->source_type());
}//lesThan(...)


void Measurement::reset()
{
  live_time_ = 0.0f;
  real_time_ = 0.0f;

  sample_number_ = 1;
  occupied_ = UnknownOccupancyStatus;
  gamma_count_sum_ = 0.0;
  neutron_counts_sum_ = 0.0;
  speed_ = 0.0;
  detector_name_ = "";
  detector_number_ = -1;
  detector_type_ = "";
  quality_status_ = Missing;

  source_type_       = UnknownSourceType;
  energy_calibration_model_ = UnknownEquationType;

  contained_neutron_ = false;

  latitude_ = longitude_ = -999.9;
  position_time_ = boost::posix_time::not_a_date_time;

  remarks_.clear();
  start_time_ = boost::posix_time::not_a_date_time;
  calibration_coeffs_.clear();
  deviation_pairs_.clear();
  channel_energies_ = std::make_shared<vector<float> >();//XXX I should test not bothering to place an empty vector in this pointer
  gamma_counts_ = std::make_shared<vector<float> >();  //XXX I should test not bothering to place an empty vector in this pointer
  neutron_counts_.clear();
}//void reset()



void Measurement::combine_gamma_channels( const size_t ncombine )
{
  if( !gamma_counts_ )
    return;
  
  const size_t nchannelorig = gamma_counts_->size();
  const size_t nnewchann = nchannelorig / ncombine;
  
  if( !nchannelorig || ncombine==1 )
    return;
  
  if( !ncombine || ((nchannelorig % ncombine) != 0) || ncombine > nchannelorig )
    throw runtime_error( "combine_gamma_channels(): invalid input." );

#if( PERFORM_DEVELOPER_CHECKS )
  const double pre_gammasum = std::accumulate( gamma_counts_->begin(),
                                            gamma_counts_->end(), double(0.0) );
  const float pre_lower_e = (*channel_energies_)[0];
  const float pre_upper_e = (*channel_energies_)[nchannelorig - ncombine];
#endif  //#if( PERFORM_DEVELOPER_CHECKS )
  
  std::shared_ptr< vector<float> > newchanneldata
                        = std::make_shared<vector<float> >( nnewchann, 0.0f );
  
  for( size_t i = 0; i < nchannelorig; ++i )
    (*newchanneldata)[i/ncombine] += (*gamma_counts_)[i];
  
  //nchannelorig is std::min( gamma_counts_->size(), channel_energies_->size() )
  //  which practically gamma_counts_->size()==channel_energies_->size(), but
  //  jic
  //for( size_t i = nchannelorig; i < gamma_counts_->size(); ++i )
    //(*newchanneldata)[nnewchann-1] += (*gamma_counts_)[i];
  
  gamma_counts_ = newchanneldata;
  
  switch( energy_calibration_model_ )
  {
    case Polynomial:
      for( size_t i = 1; i < calibration_coeffs_.size(); ++i )
        calibration_coeffs_[i] *= std::pow( float(ncombine), float(i) );
      channel_energies_ = polynomial_binning( calibration_coeffs_, nnewchann,
                                              deviation_pairs_ );
      break;
      
    case FullRangeFraction:
      channel_energies_ = fullrangefraction_binning( calibration_coeffs_,
                                                   nnewchann, deviation_pairs_ );
      break;
      
    case LowerChannelEdge:
    case UnknownEquationType:
      if( !!channel_energies_ )
      {
        std::shared_ptr<vector<float> > newbinning
                               = std::make_shared<vector<float> >(nnewchann);
        
        for( size_t i = 0; i < nchannelorig; i += ncombine )
          (*newbinning)[i/ncombine] = (*channel_energies_)[i];
        
        channel_energies_ = newbinning;
      }//if( !!channel_energies_ )
      break;
  }//switch( energy_calibration_model_ )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  const double post_gammasum = std::accumulate( gamma_counts_->begin(),
                                            gamma_counts_->end(), double(0.0) );
  const float post_lower_e = channel_energies_->front();
  const float post_upper_e = channel_energies_->back();
  
  if( fabs(post_gammasum - pre_gammasum) > (0.00001*std::max(fabs(post_gammasum),fabs(pre_gammasum))) )
  {
    char buffer[512];
    snprintf( buffer, sizeof(buffer),
             "Gamma sum changed from %f to %f while combining channels.",
              pre_gammasum, post_gammasum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }//if( gamma sum changed )
  
  if( fabs(post_lower_e - pre_lower_e) > 0.0001 )
  {
    char buffer[512];
    snprintf( buffer, sizeof(buffer), "Lower energy of spectrum changed from "
              "%f to %f while combining channels.", pre_lower_e, post_lower_e );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }//if( lower energy changed )
  
  
  if( fabs(post_upper_e - pre_upper_e) > 0.0001 )
  {
    char buffer[512];
    snprintf( buffer, sizeof(buffer), "Upper energy of spectrum changed from %f"
             " to %f while combining channels.", pre_upper_e, post_upper_e );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }//if( lower energy changed )
  
#endif  //#if( PERFORM_DEVELOPER_CHECKS )
}//void combine_gamma_channels( const size_t nchann )


size_t MeasurementInfo::do_channel_data_xform( const size_t nchannels,
                std::function< void(std::shared_ptr<Measurement>) > xform )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  size_t nchanged = 0;
  std::shared_ptr<Measurement> firstcombined;
  
  set<size_t> nchannelset, othernchannel;
  set<MeasurementCalibInfo> othercalibs;
  map<MeasurementCalibInfo, vector<std::shared_ptr<Measurement> > > calibs;
  
  for( size_t i = 0; i < measurements_.size(); ++i )
  {
    std::shared_ptr<Measurement> &m = measurements_[i];
    //    const size_t nchannel = m->num_gamma_channels();
    
    if( !m->gamma_channel_contents()
       || m->gamma_channel_contents()->size() != nchannels )
    {
      if( !!m->gamma_channel_contents() && !m->gamma_channel_contents()->empty())
      {
        othernchannel.insert( m->gamma_channel_contents()->size() );
        MeasurementCalibInfo info(m);
        info.binning.reset();
        info.original_binning.reset();
        othercalibs.insert( info );
      }
      
      continue;
    }//if( not a gamma measurment )
    
    xform( m );
    
    MeasurementCalibInfo info( m );
    info.binning.reset();
    info.original_binning.reset();
    
    if(!calibs[info].empty())
    {
      cerr << "Making it so measurment shared channel_energies_ with " << calibs[info][0].get() << endl;
      m->channel_energies_ = calibs[info][0]->channel_energies_;
    }
    
    calibs[info].push_back( m );
    if( !!m->channel_energies_ )
      nchannelset.insert( m->channel_energies_->size() );
    
    ++nchanged;
  }//for( size_t i = 0; i < measurements_.size(); ++i )
  
  //ToDo: better test this function.
  cerr << "There are calibs.size()=" << calibs.size() << endl;
  cerr << "There are othercalibs.size()=" << othercalibs.size() << endl;
  
  if( nchanged )
  {
    if( calibs.size() > 1 || othercalibs.size() > 1
       || (calibs.size()==1 && othercalibs.size()==1
           && !((calibs.begin()->first) == (*othercalibs.begin()))) )
    {
      cerr << "Un-setting properties_flags_ kHasCommonBinning bit" << endl;
      properties_flags_ &= (~kHasCommonBinning);
    }else
    {
      cerr << "Setting properties_flags_ kHasCommonBinning bit" << endl;
      properties_flags_ |= kHasCommonBinning;
    }
  
    if( (nchannelset.size() > 1) || (othernchannel.size() > 1)
       || (nchannelset.size()==1 && othernchannel.size()==1
           && (*nchannelset.begin())!=(*othernchannel.begin())) )
    {
      cerr << "Un-setting properties_flags_ kAllSpectraSameNumberChannels bit" << endl;
      properties_flags_ &= (~kAllSpectraSameNumberChannels);
    }else
    {
      cerr << "Setting properties_flags_ kAllSpectraSameNumberChannels bit" << endl;
      properties_flags_ |= kAllSpectraSameNumberChannels;
    }
  
    modifiedSinceDecode_ = modified_ = true;
  }//if( nchanged )
  
  return nchanged;
}//size_t do_channel_data_xform( const size_t nchannels, const std::function<std::shared_ptr<Measurement> > &xform )


size_t MeasurementInfo::combine_gamma_channels( const size_t ncombine,
                                                const size_t nchannels )
{
  if( ((nchannels % ncombine) != 0) || !nchannels || !ncombine )
    throw runtime_error( "MeasurementInfo::combine_gamma_channels(): invalid input" );
  
  try
  {
    auto xform = [ncombine]( std::shared_ptr<Measurement> m ){ m->combine_gamma_channels(ncombine); };
    return do_channel_data_xform( nchannels, xform );
  }catch( std::exception &e )
  {
    throw runtime_error( "MeasurementInfo::combine_gamma_channels():" + string(e.what()) );
  }
  
  return 0;
}//size_t combine_gamma_channels( const size_t, const size_t )



void MeasurementInfo::combine_gamma_channels( const size_t ncombine,
                             const std::shared_ptr<const Measurement> &meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  MeasurementShrdPtr m = measurment( meas );
  if( !m )
    throw runtime_error( "MeasurementInfo::combine_gamma_channels(): measurment"
                         " passed in is not owned by this MeasurementInfo." );
  
  m->combine_gamma_channels( ncombine );
  
  //Could actually check for a common binning, or try to share channel_energies_
  //  here, but whatever for now.
  if( measurements_.size() > 1 )
  {
    properties_flags_ &= (~kHasCommonBinning);
    properties_flags_ &= (~kAllSpectraSameNumberChannels);
  }//if( measurements_.size() > 1 )
  
  modifiedSinceDecode_ = modified_ = true;
}//void combine_gamma_channels( const int nchann, &meas )



std::vector<float> polynomial_cal_remove_first_channels( const int nchannel,
                                                  const std::vector<float> &a )
{
  const float n = static_cast<float>( nchannel );  //to avoid integer overflow
  std::vector<float> answer( a.size() );
  
  if( a.size() == 2 )
  {
    answer[0] = a[0] + n*a[1];
    answer[1] = a[1];
  }else if( a.size() == 3 )
  {
    answer[0] = a[0] + n*a[1] + n*n*a[2];
    answer[1] = a[1] + 2.0f*n*a[2];
    answer[2] = a[2];
  }else if( a.size() == 4 )
  {
    answer[0] = n*n*n*a[3] + n*n*a[2] + n*a[1] + a[0];
    answer[1] = 3.0f*n*n*a[3] + 2.0f*n*a[2] + a[1];
    answer[2] = 3.0f*n*a[3] + a[2];
    answer[3] = a[3];
  }else if( a.size() == 5 )
  {
    answer[0] = n*n*n*n*a[4] + n*n*n*a[3] + n*n*a[2] + n*a[1] + a[0];
    answer[1] = 4.0f*n*n*n*a[4] + 3.0f*n*n*a[3] + 2.0f*n*a[2] + a[1];
    answer[2] = 6.0f*n*n*a[4] + 3.0f*n*a[3] + a[2];
    answer[3] = 4.0f*n*a[4] + a[3];
    answer[4] = a[4];
  }else if( a.size() >= 6 )
  {
    //Ignore anything past the 6th coefficient (its questionalble if we even
    //  need to go up this high).
    /*
     //Sage code to solve for the coefficients.
     %var S, x, n, a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5
     S = b0 + b1*(x-n) + b2*(x-n)^2 + b3*(x-n)^3 + b4*(x-n)^4 + b5*(x-n)^5
     
     r5 = (a5 == S.expand().coefficients(x)[5][0]).solve(b5)
     r4 = (a4 == S.expand().coefficients(x)[4][0]).subs_expr(r5[0]).solve(b4)
     r3 = (a3 == S.expand().coefficients(x)[3][0]).subs_expr(r4[0],r5[0]).solve(b3)
     r2 = (a2 == S.expand().coefficients(x)[2][0]).subs_expr(r3[0],r4[0],r5[0]).solve(b2)
     r1 = (a1 == S.expand().coefficients(x)[1][0]).subs_expr(r2[0],r3[0],r4[0],r5[0]).solve(b1)
     r0 = (a0 == S.expand().coefficients(x)[0][0]).subs_expr(r1[0],r2[0],r3[0],r4[0],r5[0]).solve(b0)
     */
    answer.resize( 6, 0.0f );
    answer[0] = pow(n,5.0f)*a[5] + pow(n,4.0f)*a[4] + pow(n,3.0f)*a[3] + pow(n,2.0f)*a[2] + n*a[1] + a[0];
    answer[1] = 5*pow(n,4.0f)*a[5] + 4*pow(n,3.0f)*a[4] + 3*pow(n,2.0f)*a[3] + 2.0f*n*a[2] + a[1];
    answer[2] = 10*pow(n,3.0f)*a[5] + 6*pow(n,2.0f)*a[4] + 3*n*a[3] + a[2];
    answer[3] = 10*pow(n,2.0f)*a[5] + 4*n*a[4] + a[3];
    answer[4] = 5*n*a[5] + a[4];
    answer[5] = a[5];
  }
  
  return answer;
}//polynomial_cal_remove_first_channels(...)



void Measurement::truncate_gamma_channels( const size_t keep_first_channel,
                                          const size_t keep_last_channel,
                                          const bool keep_under_over_flow )
{
  if( !gamma_counts_ || gamma_counts_->empty() )
    return;
  
  const size_t nprevchannel = gamma_counts_->size();
  
  if( keep_last_channel >= nprevchannel )
    throw runtime_error( "Measurement::truncate_gamma_channels(): invalid upper channel." );
  
  if( keep_first_channel > keep_last_channel )
    throw runtime_error( "Measurement::truncate_gamma_channels(): invalid channel range." );
  
  double underflow = 0.0, overflow = 0.0;
  if( keep_under_over_flow )
  {
    for( size_t i = 0; i < keep_first_channel; ++i )
      underflow += (*gamma_counts_)[i];
    for( size_t i = keep_last_channel + 1; i < nprevchannel; ++i )
      overflow += (*gamma_counts_)[i];
  }//if( keep_under_over_flow )
  
  const size_t nnewchannel = 1 + keep_last_channel - keep_first_channel;
  
  std::shared_ptr<vector<float> > newchannelcounts
                            = std::make_shared<vector<float> >(nnewchannel);
  
  for( size_t i = keep_first_channel; i <= keep_last_channel; ++i )
    (*newchannelcounts)[i-keep_first_channel] = (*gamma_counts_)[i];
  
  newchannelcounts->front() += static_cast<float>( underflow );
  newchannelcounts->back()  += static_cast<float> ( overflow );

  
#if( PERFORM_DEVELOPER_CHECKS )
  if( keep_under_over_flow )
  {
    const double newsum = std::accumulate( newchannelcounts->begin(),
                                        newchannelcounts->end(), double(0.0) );
    if( fabs(newsum - gamma_count_sum_) > 0.001 )
    {
      char buffer[512];
      snprintf( buffer, sizeof(buffer),
               "Cropping channel counts resulted gamma sum disagreement, "
                "expected new sum to equal old sum, but instead got %f (new) vs"
                " %f (old).", newsum, gamma_count_sum_ );
      log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
    }
  }//if( keep_under_over_flow )
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  
  gamma_counts_ = newchannelcounts;
  
  if( !keep_under_over_flow )
  {
    gamma_count_sum_ = 0.0;
    for( const float f : *gamma_counts_ )
      gamma_count_sum_ += f;
  }//if( !keep_under_over_flow )

  
#if( PERFORM_DEVELOPER_CHECKS )
  ShrdConstFVecPtr old_binning = channel_energies_;
#endif  //#if( PERFORM_DEVELOPER_CHECKS )
  
  const int n = static_cast<int>( keep_first_channel );
  const vector<float> a = calibration_coeffs_;

  
  switch( energy_calibration_model_ )
  {
    case Polynomial:
    {
      calibration_coeffs_ = polynomial_cal_remove_first_channels( n, a );      
      channel_energies_ = polynomial_binning( calibration_coeffs_, nnewchannel,
                                              deviation_pairs_ );
      break;
    }//case Polynomial:
      
    case FullRangeFraction:
    {
      const vector<float> oldpoly = fullrangefraction_coef_to_polynomial( a,
                                                              nprevchannel );
      
      const vector<float> newpoly = polynomial_cal_remove_first_channels( n,
                                                                    oldpoly );
      
      vector<float> newfwf = polynomial_coef_to_fullrangefraction(
                                                        newpoly, nnewchannel );
      
      if( a.size() > 4 )
      {
//        const float x = static_cast<float>(i)/static_cast<float>(nbin);
//        const float low_e_coef = (a.size() > 4) ? a[4] : 0.0f;
//        val += low_e_coef / (1.0f+60.0f*x);

        //XXX - fix up here, if doable.  I dont think this can be exactly done,
        //  but maybye make it match up at the bottom binning, since this term
        //  only effects low energy.
      }//if( a.size() > 4 )
      
      calibration_coeffs_ = newfwf;
      channel_energies_ = fullrangefraction_binning( calibration_coeffs_,
                                            nnewchannel, deviation_pairs_ );
      break;
    }//case FullRangeFraction:
      
    case LowerChannelEdge:
    case UnknownEquationType:
      if( !!channel_energies_ )
      {
        std::shared_ptr< vector<float> > newbinning
                             = std::make_shared<vector<float> >(nnewchannel);
        
        for( size_t i = keep_first_channel; i <= keep_last_channel; ++i )
          (*newbinning)[i-keep_first_channel] = (*channel_energies_)[i];
        
        channel_energies_ = newbinning;
      }//if( !!channel_energies_ )
      break;
  }//switch( energy_calibration_model_ )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  if( !!old_binning && !!channel_energies_ )
  {
    for( size_t i = keep_first_channel; i <= keep_last_channel; ++i )
    {
      const float newval = (*channel_energies_)[i-keep_first_channel];
      const float oldval = (*old_binning)[i];
      if( fabs(newval-oldval) > 0.001f )
      {
        char buffer[512];
        snprintf( buffer, sizeof(buffer),
                  "Cropping channel counts resulted in disagreement of channel"
                  " energies by old channel %i which had energy %f but now has"
                  " energy %f (new channel number %i)",
                  int(i), oldval, newval, int(i-keep_first_channel) );
        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
        break;
      }//if( fabs(newval-oldval) > 0.001f )
    }//for( size_t i = keep_first_channel; i <= keep_last_channel; ++i )
  }//if( !!old_binning && !!channel_energies_ )
#endif  //#if( PERFORM_DEVELOPER_CHECKS )
}//size_t Measurement::truncate_gamma_channels(...)


size_t MeasurementInfo::truncate_gamma_channels( const size_t keep_first_channel,
                               const size_t keep_last_channel,
                               const size_t nchannels,
                               const bool keep_under_over_flow )
{
  try
  {
    auto xform = [=]( std::shared_ptr<Measurement> m ){
      m->truncate_gamma_channels(keep_first_channel, keep_last_channel, keep_under_over_flow);
    };
    
    return do_channel_data_xform( nchannels, xform );
  }catch( std::exception &e )
  {
    throw runtime_error( "MeasurementInfo::truncate_gamma_channels():"
                          + string(e.what()) );
  }
  
  return 0;
}//size_t truncate_gamma_channels(...)


void MeasurementInfo::truncate_gamma_channels( const size_t keep_first_channel,
                             const size_t keep_last_channel,
                             const bool keep_under_over_flow,
                             const std::shared_ptr<const Measurement> &meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  
  MeasurementShrdPtr m = measurment( meas );
  if( !m )
    throw runtime_error( "MeasurementInfo::truncate_gamma_channels(): measurment"
                        " passed in is not owned by this MeasurementInfo." );
  
  m->truncate_gamma_channels( keep_first_channel, keep_last_channel,
                              keep_under_over_flow );
  
  //Could actually check for a common binning, or try to share channel_energies_
  //  here, but whatever for now.
  if( measurements_.size() > 1 )
  {
    properties_flags_ &= (~kHasCommonBinning);
    properties_flags_ &= (~kAllSpectraSameNumberChannels);
  }//if( measurements_.size() > 1 )
  
  modifiedSinceDecode_ = modified_ = true;
}//void MeasurementInfo::truncate_gamma_channels(...)


void Measurement::set_2006_N42_spectrum_node_info( const rapidxml::xml_node<char> *spectrum )
{

  if( !spectrum )
    throw runtime_error( "set_2006_N42_spectrum_node_info: Recieved NULL 'Spectrum' node" );

  const string xmlns = get_n42_xmlns( spectrum );
  
  for( const rapidxml::xml_node<char> *remark_node = xml_first_node_nso( spectrum, "Remark", xmlns );
       remark_node;
       remark_node = XML_NEXT_TWIN(remark_node) )
  {
    string remark_from_node = xml_value_str( remark_node );
    
    vector<string> remark_lines;
    split( remark_lines, remark_from_node, "\r\n" );
    
    for( string &remark : remark_lines )
    {
      trim( remark );
      if( remark.empty() )
        continue;
    
      if( UtilityFunctions::istarts_with( remark, "Title:") )
      {
        remark = remark.substr(6);
        trim( remark );
        title_ = remark;
        continue;
      }
    
      remarks_.push_back( remark );
      
      if( sample_number_ < 0 )
      {
        sample_number_ = sample_num_from_remark( remarks_.back() );
      }else
      {
        const int samplen = sample_num_from_remark( remarks_.back() );
        if( samplen!=sample_number_ && samplen>=0 )
          cerr << "Got multiple sample_nums: " << sample_number_
               << " vs: " << samplen << " from " << remarks_.back() << endl;
        
      //marking it intrinsic activity will happen further down from the 'ID'
      //  attribute, so we wont wast cpu time here checking the remark for itww
//      if( UtilityFunctions::icontains( remark, "intrinsic activity") )
//        source_type_ = Measurement::IntrinsicActivity;
      }

      const float thisspeed = speed_from_remark( remark );
      if( thisspeed > 0.0f )
        speed_ = thisspeed;

      const string found_detector_name = detector_name_from_remark( remarks_.back() );
      if( !found_detector_name.empty() && detector_name_.empty() )
        detector_name_ = found_detector_name;
      else if( detector_name_ != found_detector_name )
        cerr << "Got multiple detector names: " << detector_name_
             << " vs " << found_detector_name << endl;
    }//for( string remark, remark_lines )
  }//for( loop over remark_nodes )

  const rapidxml::xml_attribute<char> *sample_num_att = spectrum->first_attribute( "SampleNumber", 12 );
  if( sample_num_att )
  {
    const string strvalue = xml_value_str( sample_num_att );

    if( sample_number_ >= 2 )
      cerr << SRC_LOCATION << "\n\tWarning: replacing sample_number_="
           << sample_number_ << " with whatever will come from "
           << strvalue << endl;

    
    if( !toInt( strvalue, sample_number_ ) && !strvalue.empty() )
      cerr << SRC_LOCATION << "\n\tWarning: couldnt convert '" << strvalue
           << "' to an int" << endl;
    else sample_number_ = 1;
  }//if( sample_num_att )
  
  const rapidxml::xml_node<char> *src_type_node = xml_first_node_nso( spectrum, "SourceType", xmlns );
  
  if( src_type_node )
  {
    if( XML_VALUE_ICOMPARE(src_type_node, "Item") )
      source_type_ = Foreground;
    else if( XML_VALUE_ICOMPARE(src_type_node, "Background") )
      source_type_ = Background;
    else if( XML_VALUE_ICOMPARE(src_type_node, "Calibration") )
      source_type_ = Calibration;
    else if( XML_VALUE_ICOMPARE(src_type_node, "Stabilization") ) //RadSeeker HPRDS files have the "Stabilization" source type, which looks like an intrinsic source
      source_type_ = IntrinsicActivity;
    else if( XML_VALUE_ICOMPARE(src_type_node, "IntrinsicActivity") )
      source_type_ = IntrinsicActivity;
    else
      source_type_ = UnknownSourceType;
    
    
  }//if( src_type_node )
  
  const rapidxml::xml_attribute<char> *id_att = spectrum->first_attribute( "ID", 2, false );
  if( id_att )
  {
    if( XML_VALUE_ICOMPARE( id_att, "intrinsicActivity" ) )
      source_type_ = IntrinsicActivity;
  }// if( id_att )

  const rapidxml::xml_node<char> *uccupied_node = xml_first_node_nso( spectrum, "Occupied", xmlns );

  try
  {
    if( is_occupied( uccupied_node ) ) occupied_ = Occupied;
    else                               occupied_ = NotOccupied;
  }catch(...){                         occupied_ = UnknownOccupancyStatus; }

  const rapidxml::xml_node<char> *det_type_node = xml_first_node_nso( spectrum, "DetectorType", xmlns );
  if( det_type_node && det_type_node->value_size() )
    detector_type_ = xml_value_str( det_type_node );

  const rapidxml::xml_attribute<char> *quality_attrib = spectrum->first_attribute( "Quality", 7 );
  if( quality_attrib && quality_attrib->value_size() )
  {
    if( XML_VALUE_ICOMPARE( quality_attrib, "Good" ) )
      quality_status_ = Good;
    else if( XML_VALUE_ICOMPARE( quality_attrib, "Suspect" ) )
      quality_status_ = Suspect;
    else if( XML_VALUE_ICOMPARE( quality_attrib, "Bad" ) )
      quality_status_ = Bad;
    else if( XML_VALUE_ICOMPARE( quality_attrib, "Missing" )
             || XML_VALUE_ICOMPARE( quality_attrib, "Unknown" ) )
      quality_status_ = Missing;
    else
    {
      cerr << SRC_LOCATION << "\n\tWarning: unknow quality status: '"
           << quality_attrib->value() << "' setting to Missing." << endl;
      quality_status_ = Missing;
    }//if(0.../else/...
  }//if( quality_attrib is valid )

  const rapidxml::xml_attribute<char> *detector_attrib = find_detector_attribute( spectrum );

  if( detector_attrib && detector_attrib->value_size() )
  {
    if(!detector_name_.empty())
      cerr << SRC_LOCATION << "\n\tWarning: replacing detector name '"
           << detector_name_ << "'' with '" << xml_value_str(detector_attrib) << "'"
            << endl;
    detector_name_ = xml_value_str(detector_attrib);
  }//if( detector_attrib && detector_attrib->value() )

  const rapidxml::xml_node<char> *live_time_node  = xml_first_node_nso( spectrum, "LiveTime", xmlns );
  const rapidxml::xml_node<char> *real_time_node  = xml_first_node_nso( spectrum, "RealTime", xmlns );
  const rapidxml::xml_node<char> *start_time_node = xml_first_node_nso( spectrum, "StartTime", xmlns );

  if( live_time_node )
    live_time_ = time_duration_in_seconds( live_time_node->value(), live_time_node->value_size() );
  if( real_time_node )
    real_time_ = time_duration_in_seconds( real_time_node->value(), real_time_node->value_size() );
  
  if( !start_time_node && spectrum->parent() )
    start_time_node = xml_first_node_nso( spectrum->parent(), "StartTime", xmlns );
  
  if( start_time_node )
    start_time_ = time_from_string( xml_value_str(start_time_node).c_str() );

  
  //XXX Things we should look for!
  //Need to handle case <Calibration Type="FWHM" FWHMUnits="Channels"> instead of right now only handling <Calibration Type="Energy" EnergyUnits="keV">

  const rapidxml::xml_node<char> *channel_data_node = xml_first_node_nso( spectrum, "ChannelData", xmlns );  //can have attribute Compression, Start(The channel number (one-based) of the first value in this element), ListMode(string)
  
  if( !channel_data_node )
  {
    //The N42 analysis result file refU35CG8VWRM get here a lot (its not a valid
    //  spectrum file)
    const char *msg = "Error, didnt find <ChannelData> under <Spectrum>";
    throw runtime_error( msg );
  }//if( !channel_data_node )
  
  const rapidxml::xml_attribute<char> *compress_attrib = channel_data_node->first_attribute(
                                                          "Compression", 11 );
  const string compress_type = xml_value_str( compress_attrib );
  ShrdFVecPtr contents = std::make_shared< vector<float> >();
  
  //Some variants have a <Data> tag under the <ChannelData> node.
  const rapidxml::xml_node<char> *datanode = xml_first_node_nso( channel_data_node, "Data", xmlns );
  if( datanode && datanode->value_size() )
    channel_data_node = datanode;

  
  const bool compressed_zeros = iequals(compress_type, "CountedZeroes");
  
  //XXX - this next call to split_to_floats(...) is not safe for non-destructively parsed XML!!!  Should fix.
  UtilityFunctions::split_to_floats( channel_data_node->value(), *contents, " ,\r\n\t", compressed_zeros );
//  UtilityFunctions::split_to_floats( channel_data_node->value(), channel_data_node->value_size(), *contents );

  if( compressed_zeros )
    expand_counted_zeros( *contents, *contents );
  else if( (compress_type!="") && (contents->size()>2) && !UtilityFunctions::iequals(compress_type, "None" ) )
  {
    stringstream msg;
    msg << SRC_LOCATION << "\n\tUnknown spectrum compression type: '"
        << compress_type << "', Compression atribute value='"
        << xml_value_str(compress_attrib) << "'";

    cerr << msg.str() << endl;
    throw runtime_error( msg.str() );
  }//if( counted zeros ) / else some other compression

  //Fix cambio zero compression
  if( compressed_zeros )
  {
    for( float &val : *contents )
    {
      if( val > 0.0f && val <= 2.0f*FLT_MIN )
        val = 0.0f;
    }
  }//if( compressed_zeros )
  
  const rapidxml::xml_attribute<char> *type_attrib = spectrum->first_attribute( "Type", 4 );
  
  if( !type_attrib )
    type_attrib = spectrum->first_attribute( "DetectorType", 12 );
  
  if( !type_attrib && spectrum->parent() )
    type_attrib = spectrum->parent()->first_attribute( "DetectorType", 12 );            //<SpectrumMeasurement>
  if( !type_attrib && spectrum->parent() && spectrum->parent()->parent() )
    type_attrib = spectrum->parent()->parent()->first_attribute( "DetectorType", 12 );  //<DetectorMeasurement> node
  
  bool is_gamma = true;

  try
  {
    is_gamma = is_gamma_spectrum( detector_attrib, type_attrib,
                                  det_type_node, spectrum );
  }catch( std::exception &e )
  {
    if( !channel_data_node || (channel_data_node->value_size() < 10) )
      cerr << SRC_LOCATION << "\n\t: Coudlnt determine detector type: "
           << e.what() << endl << "\tAssuming is a gamma detector" << endl;
  }

  if( is_gamma )
  {
    //The bellow handles a special case for Raytheon-Variant L-1 (see refSJHFSW1DZ4)
    const rapidxml::xml_node<char> *specsize = spectrum->first_node( "ray:SpectrumSize", 16 );
    if( !!contents && contents->size() && specsize && specsize->value_size() )
    {
      vector<int> sizes;
      const char *str = specsize->value();
      const size_t strsize = specsize->value_size();
      if( UtilityFunctions::split_to_ints( str, strsize, sizes )
         && sizes.size() == 1 )
      {
        const size_t origlen = gamma_counts_->size();
        const size_t newlen = static_cast<size_t>(sizes[0]);
        if( newlen >= 64
            && newlen != origlen
            && newlen < origlen
            && (origlen % newlen)==0 )
        {
          contents->resize( newlen );
          
#if( PERFORM_DEVELOPER_CHECKS )
          char buffer[512];
          snprintf( buffer, sizeof(buffer),
                    "Reducing channel data from %i to %i channels on advice of"
                    " <ray:SpectrumSize>; note that this is throwing away %i"
                    " channels", int(origlen), int(newlen), int(origlen-newlen) );
          log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
#endif  //#if( PERFORM_DEVELOPER_CHECKS )
        }
      }//if( UtilityFunctions::split_to_ints( str, strsize, sizes ) )
    }//if( specsize_node && specsize_node->value_size() )
    
    
    contained_neutron_ = false;
    gamma_counts_ = contents;

    
    for( const rapidxml::xml_node<char> *calibration_node = xml_first_node_nso( spectrum, "Calibration", xmlns );
        calibration_node;
        calibration_node = XML_NEXT_TWIN(calibration_node) )
    {
      try
      {
        decode_n42_2006_binning( calibration_node,
                                calibration_coeffs_, energy_calibration_model_ );
        
        break;
      } catch ( std::exception &e )
      {
        calibration_coeffs_.clear();
        energy_calibration_model_ = UnknownEquationType;
      }
    }
    
//    const rapidxml::xml_node<char> *nonlinarity_node = xml_first_node_nso( spectrum, "NonlinearityCorrection", xmlns );
    for( const float x : *(gamma_counts_) )
      gamma_count_sum_ += x;
  }else
  {
    contained_neutron_ = true;
    if( neutron_counts_.size() < contents->size() )
      neutron_counts_.resize( contents->size(), 0.0 );

    for( size_t i = 0; i < contents->size(); ++i )
    {
      neutron_counts_[i] += contents->operator[](i);
      neutron_counts_sum_ += contents->operator[](i);
    }//for( loop over neutron counts )
  }//if( is_gamma ) / else
}//void set_2006_N42_spectrum_node_info( rapidxml::xml_node<char> *measurementNode )



void Measurement::set_n42_2006_spectrum_calibration_from_id( const rapidxml::xml_node<char> *doc_node,
                                                const rapidxml::xml_node<char> *spectrum_node )
{
  if( !doc_node || !spectrum_node )
    return;

  const string xmlns = get_n42_xmlns( spectrum_node );
  
  const rapidxml::xml_attribute<char> *cal_IDs_att = XML_FIRST_ATTRIB( spectrum_node, "CalibrationIDs" );
  
  vector<string> cal_ids;
  split( cal_ids, xml_value_str(cal_IDs_att), " \t" );

  
  //If there is only calibration node, but the ID we want doesnt match the only
  //  calibration node, then let the match work off of the first two charcters
  //  of the ID.
  //This is to allow the N42 files in C:\GADRAS\Detector\HPRDS\SmithsNaI
  //  to decode.
  
  int ncalnodes = 0;
  for( const rapidxml::xml_node<char> *node = xml_first_node_nso( doc_node, "Calibration", xmlns );
      node; node = XML_NEXT_TWIN(node) )
  {
    ++ncalnodes;
  }
  
  if( cal_ids.empty() && ncalnodes != 1 )
    return;
  
  for( const rapidxml::xml_node<char> *cal_node = xml_first_node_nso( doc_node, "Calibration", xmlns );
       cal_node;
       cal_node = XML_NEXT_TWIN(cal_node) )
  {
    const rapidxml::xml_attribute<char> *id_att = cal_node->first_attribute( "ID", 2, false );

//    if( id_att && id_att->value_size()
//        && compare( cal_IDs_att->value(), cal_IDs_att->value_size(),
//                    id_att->value(), id_att->value_size(), case_sensitive) )
    const string id = xml_value_str( id_att );
    bool id_match = (!id.empty() && (std::find(cal_ids.begin(), cal_ids.end(), id) != cal_ids.end()));
    
    if( !id_match && ncalnodes==1 && id.empty() )
      id_match = true;
    
    if( !id_match && ncalnodes==1 && cal_ids.empty() )
      id_match = true;
    
    if( !id_match && ncalnodes==1 && cal_ids.size()==1 )
    {
      string calid = *cal_ids.begin();
      calid = calid.size()>1 ? calid.substr(0,2) : calid;
      string idid = id.size()>1 ? id.substr(0,2) : id;
      id_match = (calid == idid);
    }
    
    if( id_match )
    {
      const rapidxml::xml_attribute<char> *type_att = XML_FIRST_ATTRIB( cal_node, "Type" );
      const rapidxml::xml_attribute<char> *unit_att = XML_FIRST_ATTRIB( cal_node, "EnergyUnits" );

      if( type_att && !XML_VALUE_COMPARE(type_att, "Energy") )
      {
        continue;
      }//if( not energy calibration node )

      float units = 1.0f;
      if( unit_att )
      {
        if( XML_VALUE_ICOMPARE(unit_att, "eV") )
          units = 0.001f;
        else if( XML_VALUE_ICOMPARE(unit_att, "keV") )
          units = 1.0f;
        else if( XML_VALUE_ICOMPARE(unit_att, "MeV") )
          units = 1000.0f;
      }//if( unit_att )

      const rapidxml::xml_node<char> *array_node = xml_first_node_nso( cal_node, "ArrayXY", xmlns );
      const rapidxml::xml_node<char> *eqn_node = xml_first_node_nso( cal_node, "Equation", xmlns ); //hmm, was this a mistake

      if( array_node && !eqn_node )
      {
//This function has been added in for FLIR identiFINDER N42 files
//        <Calibration Type="Energy" ID="calibration" EnergyUnits="keV">
//          <ArrayXY X="Channel" Y="Energy">
//            <PointXY>
//              <X>1 0</X>
//              <Y>3 0</Y>
//            </PointXY>
//          </ArrayXY>
//        </Calibration>

        const rapidxml::xml_node<char> *point_node = xml_first_node_nso( array_node, "PointXY", xmlns );
        if( !point_node )
          continue;

        const rapidxml::xml_node<char> *x_node = xml_first_node_nso( point_node, "X", xmlns );
        const rapidxml::xml_node<char> *y_node = xml_first_node_nso( point_node, "Y", xmlns );

        if( x_node && x_node->value_size() && y_node && y_node->value_size() )
        {
          //y_node ex: "3 0"
          float coeffval = 0.0;
          if( xml_value_to_flt(y_node, coeffval) )
          {
            energy_calibration_model_ = FullRangeFraction;

            calibration_coeffs_.clear();
            calibration_coeffs_.push_back( 0.0 );
            if( gamma_counts_ && gamma_counts_->size() )
            {
              calibration_coeffs_.push_back( gamma_counts_->size() * coeffval );
            }else
            {
              cerr << SRC_LOCATION << "\n\tI expected gamma counts to be filled out by now" << endl;
              calibration_coeffs_.push_back( 3000.0  );
            }
          }else
          {
            cerr << SRC_LOCATION << "\nExpected XML float list" << endl;
  //          calibration_coeffs_.push_back( 3000.0 / gamma_counts_->size() );
          }//Faile
        }//if( x and y nodes are valid )
        
        if( units != 1.0f )
          for( size_t i = 0; i < calibration_coeffs_.size(); ++i )
            calibration_coeffs_[i] *= units;
        
        return;
      }else if( eqn_node )
      {
        try
        {
          decode_n42_2006_binning( cal_node, calibration_coeffs_, energy_calibration_model_ );
          return;
        }catch( std::exception & )
        {
          //
        }
        
      }//if( array_node ) else if( eqn_node )
    }//if( this is the calibration we want )
  }//for( loop over calibrations )
}//set_n42_2006_spectrum_calibration_from_id


MeasurementShrdPtr MeasurementInfo::measurment( MeasurementConstShrdPtr meas )
{
  for( const auto &m : measurements_ )
  {
    if( m == meas )
      return m;
  }//for( const auto &m : measurements_ )
  
  return MeasurementShrdPtr();
}//measurment(...)

//set_live_time(...) and set_real_time(...) update both the measurment
//  you pass in, as well as *this.  If measurment is invalid, or not in
//  measurments_, than an exception is thrown.
void MeasurementInfo::set_live_time( const float lt,
                                     MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_live_time(...): measurment"
                         " passed in didnt belong to this MeasurementInfo" );

  const float oldLifeTime = meas->live_time();
  ptr->live_time_ = lt;
  gamma_live_time_ += (lt-oldLifeTime);
  modified_ = modifiedSinceDecode_ = true;
}//set_live_time(...)

void MeasurementInfo::set_real_time( const float rt, MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_real_time(...): measurment"
                         " passed in didnt belong to this MeasurementInfo" );

  const float oldRealTime = ptr->live_time();
  ptr->real_time_ = rt;
  gamma_real_time_ += (rt - oldRealTime);
  modified_ = modifiedSinceDecode_ = true;
}//set_real_time(...)


void MeasurementInfo::add_measurment( MeasurementShrdPtr meas,
                                      const bool doCleanup )
{
  if( !meas )
    return;
  
  vector< MeasurementShrdPtr >::iterator meas_pos;
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  meas_pos = lower_bound( measurements_.begin(), measurements_.end(),
                          meas, &Measurement::compare_by_sample_det_time );
 
  if( (meas_pos!=measurements_.end()) && ((*meas_pos)==meas) )
    throw runtime_error( "MeasurementInfo::add_measurment: duplicate meas" );
  
  //Making sure detector names/numbers are kept track of here instead of in
  //  cleanup_after_load() makes sure to preserve sample and detector numbers
  //  of the Measurments already in this MeasurementInfo object
  const string &detname = meas->detector_name_;
  vector<std::string>::const_iterator namepos
         = std::find( detector_names_.begin(), detector_names_.end(), detname );
  
  if( namepos == detector_names_.end() )
  {
    detector_names_.push_back( detname );
    int detnum = -1;
    for( const int d : detector_numbers_ )
      detnum = std::max( d, detnum );
    detnum += 1;
    meas->detector_number_ = detnum;
    detector_numbers_.push_back( detnum );
    if( meas->contained_neutron_ )
      neutron_detector_names_.push_back( detname );
  }else
  {
    const size_t index = namepos - detector_names_.begin();
    meas->detector_number_ = detector_numbers_.at(index);
  }//if( we dont already have a detector with this name ) / else
  
  const int detnum = meas->detector_number_;
  int samplenum = meas->sample_number();
  
  meas_pos = lower_bound( measurements_.begin(), measurements_.end(),
                         MeasurementShrdPtr(),
                         MeasurementInfoLessThan(samplenum, detnum) );
  
  if( meas_pos != measurements_.end()
     && (*meas_pos)->sample_number()==samplenum
     && (*meas_pos)->detector_number()==detnum )
  {
    const int last_sample = (*sample_numbers_.rbegin());
    meas_pos = lower_bound( measurements_.begin(), measurements_.end(),
                            MeasurementShrdPtr(),
                            MeasurementInfoLessThan(last_sample, detnum) );
    if( meas_pos == measurements_.end()
       || (*meas_pos)->sample_number()!=last_sample
       || (*meas_pos)->detector_number()!=detnum  )
    {
      samplenum = last_sample;
    }else
    {
      samplenum = last_sample + 1;
    }
    
    meas->sample_number_ = samplenum;
  }//if( we need to modify the sample number )
  
  sample_numbers_.insert( meas->sample_number_ );
  
  meas_pos = upper_bound( measurements_.begin(), measurements_.end(),
                          meas, &Measurement::compare_by_sample_det_time );
  
  measurements_.insert( meas_pos, meas );
  
  if( doCleanup )
  {
    cleanup_after_load();
  }else
  {
    gamma_count_sum_    += meas->gamma_count_sum_;
    neutron_counts_sum_ += meas->neutron_counts_sum_;
    gamma_live_time_    += meas->live_time_;
    gamma_real_time_    += meas->real_time_;
    
    sample_to_measurments_.clear();
    for( size_t measn = 0; measn < measurements_.size(); ++measn )
    {
      MeasurementShrdPtr &meas = measurements_[measn];
      sample_to_measurments_[meas->sample_number_].push_back( measn );
    }
  }//if( doCleanup ) / else
  
  modified_ = modifiedSinceDecode_ = true;
}//void add_measurment( MeasurementShrdPtr meas, bool doCleanup )


void MeasurementInfo::remove_measurments(
                                         const vector<MeasurementConstShrdPtr> &meas )
{
  if( meas.empty() )
    return;
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const size_t norigmeas = measurements_.size();
  if( meas.size() > norigmeas )
    throw runtime_error( "MeasurementInfo::remove_measurments:"
                        " to many input measurments to remove" );
  
  //This bellow implementation is targeted for MeasurementInfo's with lots of
  //  measurements_, and empiracally is much faster than commented out
  //  implementation bellow (which is left in because its a bit 'cleaner')
  vector<bool> keep( norigmeas, true );
  
  for( size_t i = 0; i < meas.size(); ++i )
  {
    const MeasurementConstShrdPtr &m = meas[i];
    
    map<int, std::vector<size_t> >::const_iterator iter
    = sample_to_measurments_.find( m->sample_number_ );
    
    if( iter != sample_to_measurments_.end() )
    {
      const vector<size_t> &indexs = iter->second;
      size_t i;
      for( i = 0; i < indexs.size(); ++i )
      {
        if( measurements_[indexs[i]] == m )
        {
          keep[indexs[i]] = false;
          break;
        }
      }//for( size_t i = 0; i < indexs.size(); ++i )
      
      if( i == indexs.size() )
        throw runtime_error( "MeasurementInfo::remove_measurments: invalid meas" );
    }//if( iter != sample_to_measurments_.end() )
  }//for( size_t i = 0; i < meas.size(); ++i )
  
  vector< MeasurementShrdPtr > surviving;
  surviving.reserve( norigmeas - meas.size() );
  for( size_t i = 0; i < norigmeas; ++i )
  {
    if( keep[i] )
      surviving.push_back( measurements_[i] );
  }
  
  measurements_.swap( surviving );
  
  /*
   for( size_t i = 0; i < meas.size(); ++i )
   {
   const MeasurementConstShrdPtr &m = meas[i];
   vector< MeasurementShrdPtr >::iterator pos;
   if( measurements_.size() > 100 )
   {
   pos = std::lower_bound( measurements_.begin(), measurements_.end(),
   m, &Measurement::compare_by_sample_det_time );
   }else
   {
   pos = std::find( measurements_.begin(), measurements_.end(), m );
   }
   
   if( pos == measurements_.end() || ((*pos)!=m) )
   throw runtime_error( "MeasurementInfo::remove_measurments: invalid meas" );
   
   measurements_.erase( pos );
   }
   */
  
  cleanup_after_load();
  modified_ = modifiedSinceDecode_ = true;
}//void remove_measurments( const vector<MeasurementConstShrdPtr> meas )


void MeasurementInfo::remove_measurment( MeasurementConstShrdPtr meas,
                                         const bool doCleanup )
{
  if( !meas )
    return;
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  vector< MeasurementShrdPtr >::iterator pos
                = std::find( measurements_.begin(), measurements_.end(), meas );
  
  if( pos == measurements_.end() )
    throw runtime_error( "MeasurementInfo::remove_measurment: invalid meas" );
  
  measurements_.erase( pos );
  
  //Could actually fix up detector_names_, detector_numbers_,
  //  neutron_detector_names_,
  
  if( doCleanup )
  {
    cleanup_after_load();
  }else
  {
    gamma_count_sum_    -= meas->gamma_count_sum_;
    neutron_counts_sum_ -= meas->neutron_counts_sum_;
    gamma_live_time_    -= meas->live_time_;
    gamma_real_time_    -= meas->real_time_;
    
    sample_numbers_.clear();
    sample_to_measurments_.clear();
    
    //should update detector names and numbers too
    set<string> detectornames;
    
    const size_t nmeas = measurements_.size();
    for( size_t measn = 0; measn < nmeas; ++measn )
    {
      const int samplenum = measurements_[measn]->sample_number_;
      sample_numbers_.insert( samplenum );
      sample_to_measurments_[samplenum].push_back( measn );
      detectornames.insert( measurements_[measn]->detector_name_ );
    }
    
    //check that detectornames matches
    for( size_t i = 0; i < detector_names_.size();  )
    {
      if( !detectornames.count(detector_names_[i]) )
      {
        detector_names_.erase( detector_names_.begin() + i );
        detector_numbers_.erase( detector_numbers_.begin() + i );
        continue;
      }
      
      ++i;
    }//for( size_t i = 0; i < detector_names_.size();  )
  }//if( doCleanup ) / else
  
  modified_ = modifiedSinceDecode_ = true;
}//void remove_measurment( MeasurementShrdPtr meas, bool doCleanup );


void MeasurementInfo::set_start_time( const boost::posix_time::ptime &timestamp,
                    const MeasurementConstShrdPtr meas  )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_start_time(...): measurment"
                        " passed in didnt belong to this MeasurementInfo" );
  
  ptr->set_start_time( timestamp );
  modified_ = modifiedSinceDecode_ = true;
}//set_start_time(...)

void MeasurementInfo::set_remarks( const std::vector<std::string> &remarks,
                 const MeasurementConstShrdPtr meas  )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_remarks(...): measurment"
                        " passed in didnt belong to this MeasurementInfo" );
  
  ptr->set_remarks( remarks );
  modified_ = modifiedSinceDecode_ = true;
}//set_remarks(...)

void MeasurementInfo::set_source_type( const Measurement::SourceType type,
                                    const MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_source_type(...): measurment"
                        " passed in didnt belong to this MeasurementInfo" );
  
  ptr->set_source_type( type );
  modified_ = modifiedSinceDecode_ = true;
}//set_source_type(...)


void MeasurementInfo::set_position( double longitude, double latitude,
                                    boost::posix_time::ptime position_time,
                                    const MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_position(...): measurment"
                        " passed in didnt belong to this MeasurementInfo" );
  
  ptr->longitude_ = longitude;
  ptr->latitude_ = latitude;
  ptr->position_time_ = position_time;
  
  
  int nGpsCoords = 0;
  mean_latitude_ = mean_longitude_ = 0.0;
  for( auto &meas : measurements_ )
  {
    if( meas->has_gps_info() )
    {
      ++nGpsCoords;
      mean_latitude_ += meas->latitude();
      mean_longitude_ += meas->longitude();
    }
  }//for( auto &meas : measurements_ )
  
  if( nGpsCoords == 0 )
  {
    mean_latitude_ = mean_longitude_ = -999.9;
  }else
  {
    mean_latitude_ /= nGpsCoords;
    mean_longitude_ /= nGpsCoords;
  }//if( !nGpsCoords ) / else
  
  
  modified_ = modifiedSinceDecode_ = true;
}//set_position(...)_


void MeasurementInfo::set_title( const std::string &title,
                                 const MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_title(...): measurment"
                        " passed in didnt belong to this MeasurementInfo" );
  
  ptr->set_title( title );
  
  modified_ = modifiedSinceDecode_ = true;
}//void MeasurementInfo::set_title(...)


void MeasurementInfo::set_contained_neutrons( const bool contained,
                                             const float counts,
                                             const MeasurementConstShrdPtr meas )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  MeasurementShrdPtr ptr = measurment( meas );
  if( !ptr )
    throw runtime_error( "MeasurementInfo::set_containtained_neutrons(...): "
                        "measurment passed in didnt belong to this "
                        "MeasurementInfo" );
  
  ptr->contained_neutron_ = contained;
  if( contained )
  {
    ptr->neutron_counts_.resize( 1 );
    ptr->neutron_counts_[0] = counts;
    ptr->neutron_counts_sum_ = counts;
  }else
  {
    ptr->neutron_counts_.resize( 0 );
    ptr->neutron_counts_sum_ = 0.0;
  }
  
  modified_ = modifiedSinceDecode_ = true;
}//void set_containtained_neutrons(...)


int MeasurementInfo::occupancy_number_from_remarks() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  for( const string &str : remarks_ )
  {
    if( istarts_with( str, "Occupancy number = " ) )
    {
      const string valstr = str.substr( 19 );
      int val;
      if( toInt(valstr,val) )
        return val;
    }else if( istarts_with( str, "OccupancyNumber" ) )
    {
      string valstr = str.substr( 15 );
      size_t pos = valstr.find_first_not_of( " :\t\n\r=" );
      if( pos != string::npos )
      {
        valstr = valstr.substr( pos );
        int val;
        if( toInt(valstr,val) )
          return val;
      }
    }
    
    //if( has "Occupancy number = " )
    
  }//for( const string &str : remarks_ )

  return -1;
}//int occupancy_number_from_remarks() const


MeasurementConstShrdPtr MeasurementInfo::measurement( const int sample_number,
                                             const std::string &det_name ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  vector<string>::const_iterator det_name_iter;
  det_name_iter = find( detector_names_.begin(), detector_names_.end(), det_name );

  if( det_name_iter == detector_names_.end() )
  {
    cerr << "Didnt find detector named '" << det_name
         << "' in detector_names_" << endl;
    return MeasurementConstShrdPtr();
  }//if( det_name_iter == detector_names_.end() )

  const size_t det_index = det_name_iter - detector_names_.begin();
  const int detector_number = detector_numbers_[det_index];

  return measurement( sample_number, detector_number );
}//MeasurementConstShrdPtr measurement( const int, const string & )



MeasurementConstShrdPtr MeasurementInfo::measurement( const int sample_number,
                                               const int detector_number ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  if( properties_flags_ & kNotSampleDetectorTimeSorted )
  {
    std::map<int, std::vector<size_t> >::const_iterator pos;
    pos = sample_to_measurments_.find( sample_number );
    
    if( pos != sample_to_measurments_.end() )
    {
      const vector<size_t> &indices = pos->second;
      for( const size_t ind : indices )
        if( measurements_[ind]->detector_number_ == detector_number )
          return measurements_[ind];
    }//if( pos != sample_to_measurments_.end() )
    
    return MeasurementConstShrdPtr();
  }//if( properties_flags_ & kNotSampleDetectorTimeSorted )
  
  
  std::vector< MeasurementShrdPtr >::const_iterator meas_pos;
  meas_pos = lower_bound( measurements_.begin(), measurements_.end(),
                          MeasurementShrdPtr(),
                          MeasurementInfoLessThan(sample_number, detector_number) );
  if( meas_pos == measurements_.end()
     || (*meas_pos)->sample_number()!=sample_number
     || (*meas_pos)->detector_number()!=detector_number )
  {
    return MeasurementConstShrdPtr();
  }//if( meas_pos == measurements_.end() )

  return *meas_pos;

  //Below kept around just incase we want to check the above
//
  //could rely on measurements_ being sorted here
//  for( const auto &meas : measurements_ )
//    if( meas->sample_number_==sample_number && meas->detector_number==detector_number )
//    {
//      assert( meas_pos != measurements_.end() );
//      MeasurementShrdPtr lb = *meas_pos;
//      if( lb.get() != meas.get() )
//      {
//        for( const auto &meas : measurements_ )
//          cout << "\t" << meas->sample_number_ << "\t" << meas->detector_number << endl;
//      }//if( lb.get() != meas.get() )
//      assert( lb.get() == meas.get() );
//      return meas;
//    }
//  cerr << "Didnt find detector " << detector_number << " for sample number "
//       << sample_number_ << endl;
//  return MeasurementConstShrdPtr();
}//MeasurementConstShrdPtr measurement( const int , onst int )



vector<MeasurementConstShrdPtr> MeasurementInfo::sample_measurements( const int sample ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  vector< MeasurementConstShrdPtr > answer;

  std::map<int, std::vector<size_t> >::const_iterator pos;
  pos = sample_to_measurments_.find( sample );

  if( pos != sample_to_measurments_.end() )
  {
    const vector<size_t> &indices = pos->second;
    for( const size_t ind : indices )
      answer.push_back( measurements_.at(ind) );
  }//if( pos != sample_to_measurments_.end() )

  return answer;

  /*
  //get all measurements_ cooresponding to 'sample_number', where
  //  sample_number may not be a zero based index (eg may start at one)
  vector< MeasurementConstShrdPtr > answer;

  if( measurements_.empty() )
    return answer;

  //we could relies on measurements being sorted and use lower_bound/upper_bound
  for( const auto &meas : measurements_ )
    if( meas->sample_number_ == sample )
      answer.push_back( meas );

  return answer;
  */
}//vector<const MeasurementConstShrdPtr> measurements( const int sample ) const







void Measurement::recalibrate_by_eqn( const std::vector<float> &eqn,
                                      const DeviationPairVec &dev_pairs,
                                      Measurement::EquationType type,
                                      ShrdConstFVecPtr binning )
{
  // Note: this will run multiple times per calibration

  if( !binning && gamma_counts_ && gamma_counts_->size() )
  {
    switch( type )
    {
      case Polynomial:
        binning = polynomial_binning( eqn, gamma_counts_->size(), dev_pairs );
      break;

      case FullRangeFraction:
        binning = fullrangefraction_binning( eqn, gamma_counts_->size(), dev_pairs );
      break;

      case LowerChannelEdge:
        cerr << SRC_LOCATION << "\n\tWarning, you probabhouldnt be calling "
             << "this function for EquationType==LowerChannelEdge, as its probably "
             << "not as efficient as the other version of this function"
             << endl;
        binning.reset( new vector<float>(eqn) );
      break;

      case UnknownEquationType:
        throw runtime_error( "Measurement::recalibrate_by_eqn(...): Can not call with EquationType==UnknownEquationType" );
    }//switch( type )
  }//if( !binning )


  if( !binning || ((binning->size() != gamma_counts_->size())
      && (type!=LowerChannelEdge || gamma_counts_->size()>binning->size())) )
  {
    stringstream msg;
    msg << "Measurement::recalibrate_by_eqn(...): "
           "You can not use this funtion to change the number"
           " of bins. The spectrum has " << gamma_counts_->size()
        << " bins, you tried to give it " << binning->size() << " bins."
        << endl;
    throw runtime_error( msg.str() );
  }//if( channel_energies_->size() != gamma_counts_->size() )


  channel_energies_   = binning;
  if( type != Measurement::LowerChannelEdge )
    calibration_coeffs_ = eqn;
  else
    calibration_coeffs_.clear();
  energy_calibration_model_  = type;
  deviation_pairs_ = dev_pairs;
  
  if( type == LowerChannelEdge )
    deviation_pairs_.clear();
  
}//void recalibrate_by_eqn( const std::vector<float> &eqn )


void Measurement::rebin_by_eqn( const std::vector<float> &eqn,
                                const DeviationPairVec &dev_pairs,
                                Measurement::EquationType type )
{
  if( !gamma_counts_ || gamma_counts_->empty() )
    throw std::runtime_error( "Measurement::rebin_by_eqn(...): "
                              "gamma_counts_ must be defined");

  ShrdConstFVecPtr binning;
  const size_t nbin = gamma_counts_->size(); //channel_energies_->size();

  switch( type )
  {
    case Measurement::Polynomial:
      binning = polynomial_binning( eqn, nbin, dev_pairs );
    break;

    case Measurement::FullRangeFraction:
      binning = fullrangefraction_binning( eqn, nbin, dev_pairs );
    break;

    case Measurement::LowerChannelEdge:
      binning.reset( new vector<float>(eqn) );
    break;

    case Measurement::UnknownEquationType:
      throw std::runtime_error( "Measurement::rebin_by_eqn(...): "
                                "can not be called with type==UnknownEquationType");
    break;
  }//switch( type )

  rebin_by_eqn( eqn, dev_pairs, binning, type );
}//void rebin_by_eqn( const std::vector<float> eqn )


void Measurement::rebin_by_eqn( const std::vector<float> &eqn,
                                const DeviationPairVec &dev_pairs,
                                ShrdConstFVecPtr binning,
                                Measurement::EquationType type )
{
  rebin_by_lower_edge( binning );

  if( type == UnknownEquationType )
    throw std::runtime_error( "Measurement::rebin_by_eqn(...): "
                              "can not be called with type==UnknownEquationType");

  if( type != LowerChannelEdge )
    calibration_coeffs_ = eqn;
  else
    calibration_coeffs_.clear();

  deviation_pairs_    = dev_pairs;
  energy_calibration_model_  = type;
}//void rebin_by_eqn(...)



void rebin_by_lower_edge( const std::vector<float> &original_energies,
                         const std::vector<float> &original_counts,
                         const std::vector<float> &new_energies,
                         std::vector<float> &resulting_counts )
{
  const size_t old_nbin = min( original_energies.size(), original_counts.size());
  const size_t new_nbin = new_energies.size();
  
  if( old_nbin < 4 )
    throw runtime_error( "rebin_by_lower_edge: input must have more than 3 bins" );
  
  if( original_energies.size() < original_counts.size() )
    throw runtime_error( "rebin_by_lower_edge: input energies and gamma counts"
                        " have mismatched number of channels" );
  
  if( new_nbin < 4 )
    throw runtime_error( "rebin_by_lower_edge: output have more than 3 bins" );
  
  resulting_counts.resize( new_nbin, 0.0f );
  
  size_t newbinnum = 0;
  while( new_energies[newbinnum] < original_energies[0] && newbinnum < (new_nbin-1) )
    resulting_counts[newbinnum++] = 0.0f;
  
  //new_energies[newbinnum] is now >= original_energies[0]
  if( newbinnum && (new_energies[newbinnum] > original_energies[0]) )
  {
    if(new_energies[newbinnum] >= original_energies[1])
    {
      resulting_counts[newbinnum-1] = original_counts[0];
      resulting_counts[newbinnum-1] += static_cast<float>( original_counts[1]
           * (double(new_energies[newbinnum]) - double(original_energies[1]))
             / (double(original_energies[2]) - double(original_energies[1])) );
    }else
    {
      resulting_counts[newbinnum-1] = static_cast<float>( original_counts[0]
              * (double(new_energies[newbinnum]) - double(original_energies[0]))
              / (double(original_energies[1]) - double(original_energies[0])) );
    }
  }
  
  size_t oldbinlow = 0, oldbinhigh = 0;
  for( ; newbinnum < new_nbin; ++newbinnum )
  {
    const double newbin_lower = new_energies[newbinnum];
    const double newbin_upper = ( ((newbinnum+1)<new_nbin)
                                 ? new_energies[newbinnum+1]
                                 : (2.0*new_energies[new_nbin-1] - new_energies[new_nbin-2]));
    
    double old_lower_low, old_lower_up, old_upper_low, old_upper_up;
    
    for( ; oldbinlow < old_nbin; ++oldbinlow )
    {
      old_lower_low = original_energies[oldbinlow];
      if( (oldbinlow+1) < old_nbin )
        old_lower_up = original_energies[oldbinlow+1];
      else
        old_lower_up = 2.0*original_energies[oldbinlow]-original_energies[oldbinlow-1];
      
      if( newbin_lower >= old_lower_low && newbin_lower < old_lower_up )
        break;
    }
    
    double sum_lower_to_upper = 0.0;
    for( oldbinhigh = oldbinlow; oldbinhigh < old_nbin; ++oldbinhigh )
    {
      old_upper_low = original_energies[oldbinhigh];
      if( (oldbinhigh+1) < old_nbin )
        old_upper_up = original_energies[oldbinhigh+1];
      else
        old_upper_up = 2.0*original_energies[oldbinhigh]-original_energies[oldbinhigh-1];
      
      if( newbin_upper >= old_upper_low && newbin_upper < old_upper_up )
        break;
      sum_lower_to_upper += original_counts[oldbinhigh];
    }
    
    //new binning goes higher than old binning, lets take care of the last
    //  fraction of a bin, and zero everything else out
    if( oldbinhigh == old_nbin )
    {
      old_upper_low = original_energies[old_nbin-1];
      old_upper_up = 2.0*original_energies[old_nbin-1]-original_energies[old_nbin-2];
      
      //maybe more numerically accurate if we remove the next line, and uncomment
      //  the next two commented lines.  Didnt do because I didnt want to test
      resulting_counts[newbinnum] = static_cast<float>( sum_lower_to_upper );
      
      if( oldbinlow != old_nbin )
      {
        const double lower_old_width = old_lower_up - old_lower_low;
        const double lower_bin_delta_counts = original_counts[oldbinlow];
        const double lower_delta_energy = double(newbin_lower)-double(original_energies[oldbinlow]);
        const double lower_frac_energy = lower_delta_energy/lower_old_width;
        
//        sum_lower_to_upper -= lower_bin_delta_counts*lower_frac_energy;
        resulting_counts[newbinnum] -= static_cast<float>(lower_bin_delta_counts*lower_frac_energy);
      }
      
//      resulting_counts[newbinnum] = static_cast<float>( sum_lower_to_upper );
      
      newbinnum++;
      while( newbinnum < new_nbin )
        resulting_counts[newbinnum++] = 0.0f;
      break;
    }
    
    const double lower_old_width = old_lower_up - old_lower_low;
    const double lower_bin_delta_counts = original_counts[oldbinlow];
    const double lower_delta_energy = double(newbin_lower)-double(original_energies[oldbinlow]);
    const double lower_frac_energy = lower_delta_energy/lower_old_width;
    
    const double upper_old_width = old_upper_up - old_upper_low;
    const double upper_bin_delta_counts = original_counts[oldbinhigh];
    const double upper_delta_energy = double(newbin_upper)-double(original_energies[oldbinhigh]);
    const double upper_frac_energy = upper_delta_energy/upper_old_width;
    
    //interpolate sum height at newbin_lower
    resulting_counts[newbinnum] = static_cast<float>( sum_lower_to_upper + upper_bin_delta_counts*upper_frac_energy-lower_bin_delta_counts*lower_frac_energy );
  }//for( ; newbinnum < (new_nbin-1); ++newbinnum )
  
  
  //capture case where new energies starts higher than the original energies, so
  //  we should put the contents of the lower energy bins into the first bin
  //  of the new counts
  if( original_energies[0] < new_energies[0] )
  {
    size_t i = 0;
    while( (original_energies[i+1]<new_energies[0]) && (i<(old_nbin-1)) )
      resulting_counts[0] += original_counts[i++];
    
    //original_energies[i] is now >= new_energies[0]
    
    if( i < old_nbin )
      resulting_counts[0] += static_cast<float>( original_counts[i]*(double(new_energies[0])-double(original_energies[i]))/(double(original_energies[i+1])-original_energies[i]) );
  }//if( original_energies[0] < new_energies[0] )
  
  //Now capture the case where the old binning extends further than new binning
  const float upper_old_energy = 2.0f*original_energies[old_nbin-1] - original_energies[old_nbin-2];
  const float upper_new_energy = 2.0f*new_energies[new_nbin-1] - new_energies[new_nbin-2];
  if( upper_old_energy > upper_new_energy )
  {
    if( oldbinhigh < (old_nbin-1) )
      resulting_counts[new_nbin-1] += original_counts[oldbinhigh] * (original_energies[oldbinhigh]-upper_new_energy)/(original_energies[oldbinhigh+1]-original_energies[oldbinhigh]);
    else
      resulting_counts[new_nbin-1] += original_counts[oldbinhigh] * (original_energies[oldbinhigh]-upper_new_energy)/(original_energies[oldbinhigh]-original_energies[oldbinhigh-1]);
    
    for( ; oldbinhigh < old_nbin; ++oldbinhigh )
      resulting_counts[new_nbin-1] += original_counts[oldbinhigh];
  }//if( original_energies.back() > new_energies.back() )
  
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  double oldsum = 0.0, newsum = 0.0;
  for( const float f : original_counts )
    oldsum += f;
  for( const float f : resulting_counts )
    newsum += f;
  
  if( (fabs(oldsum - newsum) > (old_nbin*0.1/16384.0)) && (fabs(oldsum - newsum) > 1.0E-7*oldsum) )
  {
    static std::mutex m;
    
    {
      std::lock_guard<std::mutex> loc(m);
      ofstream output( "rebin_by_lower_edge_error.txt", ios::app );
      output << "  std::vector<float> original_energies{" << std::setprecision(9);
      for( size_t i = 0; i < original_energies.size(); ++i )
          output << (i?", ":"") << original_energies[i];
      output << "};" << endl;
    
      output << "  std::vector<float> original_counts{" << std::setprecision(9);
      for( size_t i = 0; i < original_counts.size(); ++i )
          output << (i?", ":"") << original_counts[i];
      output << "};" << endl;
    
      output << "  std::vector<float> new_energies{" << std::setprecision(9);
      for( size_t i = 0; i < new_energies.size(); ++i )
          output << (i?", ":"") << new_energies[i];
      output << "};" << endl << endl;
      
      output << "  std::vector<float> new_counts{" << std::setprecision(9);
      for( size_t i = 0; i < resulting_counts.size(); ++i )
        output << (i?", ":"") << resulting_counts[i];
      output << "};" << endl << endl;
    }

  
    char buffer[1024];
    snprintf( buffer, sizeof(buffer),
              "rebin_by_lower_edge gives %1.8E vs pre rebin of %1.8E, an unacceptable error.\n",
              newsum, oldsum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
#endif
}//void rebin_by_lower_edge(...)


void Measurement::rebin_by_lower_edge( ShrdConstFVecPtr new_binning_shrd )
{
  if( new_binning_shrd.get() == channel_energies_.get() )
    return;

  if( !gamma_counts_ )
    return;
  
  {//being codeblock to do work
    if( !new_binning_shrd )
      throw runtime_error( "rebin_by_lower_edge: invalid new binning" );
    
    if( !channel_energies_ || channel_energies_->empty() )
      throw runtime_error( "rebin_by_lower_edge: no initial binning" );
    
    const size_t new_nbin = new_binning_shrd->size();
    ShrdFVecPtr rebinned_gamma_counts( new vector<float>(new_nbin) );
    
    ::rebin_by_lower_edge( *channel_energies_, *gamma_counts_,
                            *new_binning_shrd, *rebinned_gamma_counts );
    
    channel_energies_ = new_binning_shrd;
    gamma_counts_ = rebinned_gamma_counts;
    
    deviation_pairs_.clear();
    calibration_coeffs_.clear();
    energy_calibration_model_ = LowerChannelEdge;
  }//end codeblock to do works
}//void rebin_by_lower_edge(...)



void Measurement::decode_n42_2006_binning(  const rapidxml::xml_node<char> *calibration_node,
                                    vector<float> &coeffs,
                                    Measurement::EquationType &eqnmodel )
{
  coeffs.clear();

  if( !calibration_node )
  {
    const string msg = "decode_n42_2006_binning(...): Couldnt find node 'Calibration'";
    throw std::runtime_error( msg );
  }//if( !calibration_node )

  string xmlns = get_n42_xmlns(calibration_node);
  if( xmlns.empty() && calibration_node->parent() )
    xmlns = get_n42_xmlns(calibration_node->parent());
  
  const rapidxml::xml_attribute<char> *type = calibration_node->first_attribute( "Type", 4 );

  if( type && type->value_size() )
  {
    if( XML_VALUE_ICOMPARE(type, "FWHM") )
      throw runtime_error( "decode_n42_2006_binning(...): passed in FWHM cal node" );
    
    //20160601: Not adding in an explicit comparison for energy, but probably should...
    //if( !XML_VALUE_ICOMPARE(type, "Energy") )
    //  throw runtime_error( "decode_n42_2006_binning(...): passed in non-energy cal node - " + xml_value_str(type) );
  }//if( type && type->value_size() )
  
  const rapidxml::xml_attribute<char> *units = calibration_node->first_attribute( "EnergyUnits", 11 );

  string unitstr = xml_value_str(units);

  const rapidxml::xml_node<char> *equation_node = xml_first_node_nso( calibration_node, "Equation", xmlns );
  
  if( equation_node )
  {
    const rapidxml::xml_node<char> *coeff_node = xml_first_node_nso( equation_node, "Coefficients", xmlns );

    if( coeff_node )
    {
      if( coeff_node->value_size() )
      {
        coeffs.clear();
        UtilityFunctions::split_to_floats( coeff_node->value(),
                                         coeff_node->value_size(), coeffs );
      }else
      {
        //SmithsNaI HPRDS
        rapidxml::xml_attribute<char> *subeqn = coeff_node->first_attribute( "Subequation", 11 );
        if( subeqn && subeqn->value_size() )
        {
          coeffs.clear();
          UtilityFunctions::split_to_floats( subeqn->value(),
                                             subeqn->value_size(), coeffs );
        }
      }//if( coeff_node->value_size() ) / else
      
      while( coeffs.size() && coeffs.back()==0.0f )
        coeffs.erase( coeffs.end()-1 );
      
      float units = 1.0f;
      if( unitstr == "eV" )
        units = 0.001f;
      else if( unitstr == "MeV" )
        units = 1000.0f;
      
      if( units != 1.0f )
      {
        for( float &f : coeffs )
          f *= units;
      }
      
//      static std::recursive_mutex maa;
//      std::unique_lock<std::recursive_mutex> scoped_lock( maa );
//      cerr << "coeffs.size()=" << coeffs.size() << "from '"
//           << coeff_node->value() << "'"
//           << " which has length " << value_size
//           << " and strlen of " << strlen(coeff_node->value()) << endl;
    }else
    {
      const string msg = "Couldnt find node 'Coefficients'";
      throw runtime_error( msg );
    }//if( coeff_node ) / else
    
    
    const rapidxml::xml_attribute<char> *model = equation_node->first_attribute( "Model", 5 );
    //    rapidxml::xml_attribute<char> *units = equation_node->first_attribute( "Units", 5 );
    
    eqnmodel = UnknownEquationType;
    
    string modelstr = xml_value_str(model);
    if( modelstr == "Polynomial" )
      eqnmodel = Polynomial;
    else if( modelstr == "FullRangeFraction" )
      eqnmodel = FullRangeFraction;
    else if( modelstr == "LowerChannelEdge" || modelstr == "LowerBinEdge")
      eqnmodel = LowerChannelEdge;
    else if( modelstr == "Other" )
    {
      rapidxml::xml_attribute<char> *form = equation_node->first_attribute( "Form", 4 );
      const string formstr = xml_value_str(form);
      if( icontains(formstr, "Lower edge") )
        eqnmodel = LowerChannelEdge;
    }//if( modelstr == ...) / if / else
    
    //Lets try to guess the equation type for Polynomial or FullRangeFraction
    if( eqnmodel == UnknownEquationType )
    {
      if( (coeffs.size() < 5) && (coeffs.size() > 1) )
      {
        if( coeffs[1] < 10.0 )
          eqnmodel = Polynomial;
        else if( coeffs[1] > 1000.0 )
          eqnmodel = FullRangeFraction;
      }//if( coeffs.size() < 5 )
    }//if( eqnmodel == UnknownEquationType )
    
    
    if( eqnmodel == UnknownEquationType )
    {
      eqnmodel = UnknownEquationType;
      coeffs.clear();
      cerr << "Equation model is not polynomial" << endl;
      const string msg = "Equation model is not polynomial or FullRangeFraction, but is " +
      (modelstr.size()?modelstr:string("NULL"));
      cerr << msg << endl;
      throw runtime_error( msg );
    }//if( model isnt defined )
  }else
  {
    const string msg = "Couldnt find node 'Equation'";
    throw runtime_error( msg );
  }//if( equation_node ) / else

}//void decode_n42_2006_binning(...)



std::vector<float> fullrangefraction_coef_to_polynomial(
                                            const std::vector<float> &coeffs,
                                            const size_t nbinint )
{
  //nbin needs to be a float, not an integer, or else integer overflow 
  //  can happen on 32 bit systems when its cubed.
  const float nbin = static_cast<float>( nbinint );
  const size_t npars = coeffs.size();
  float a0 = 0.0f, a1 = 0.0f, a2 = 0.0f, a3 = 0.0f;
  float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;

  if( npars > 0 )
    a0 = coeffs[0];
  if( npars > 1 )
    a1 = coeffs[1];
  if( npars > 2 )
    a2 = coeffs[2];
  if( npars > 3 )
    a3 = coeffs[3];

  //Logic from from "Energy Calibration Conventions and Parameter Conversions"
  //  Dean J. Mitchell, George Lasche and John Mattingly, 6418 / MS-0791
  //The bellow "- (0.5f*a1/nbin)" should actually be "+ (0.5f*a1/nbin)"
  //  according to the paper above, I changed it to compensate for the
  //  similar hack in polynomial_coef_to_fullrangefraction.
//  c0 = a0 + (0.5f*a1/nbin) + 0.25f*a2/(nbin*nbin) - (a3/(8.0f*nbin*nbin*nbin));
//  c1 = (a1/nbin) - (a2/(nbin*nbin)) + 0.75f*(a3/(nbin*nbin*nbin));
//  c2 = (a2/(nbin*nbin)) - 1.5f*(a3/(nbin*nbin*nbin));
//  c3 = a3/(nbin*nbin*nbin);
  
  //See notes in polynomial_coef_to_fullrangefraction(...) for justification
  //  of using the bellow, instead of the above.
  c0 = a0;
  c1 = a1 / nbin;
  c2 = a2 / (nbin*nbin);
  c3 = a3 / (nbin*nbin*nbin);

  vector<float> answer;
  if( c0!=0.0f || c1!=0.0f || c2!=0.0f || c3!=0.0f )
  {
    answer.push_back( c0 );
    if( c1!=0.0f || c2!=0.0f || c3!=0.0f )
    {
      answer.push_back( c1 );
      if( c2 != 0.0f || c3!=0.0f )
      {
        answer.push_back( c2 );
        if( c3!=0.0f )
          answer.push_back( c3 );
      }
    }//if( p1!=0.0 || p2!=0.0 )
  }//if( p0!=0.0 || p1!=0.0 || p2!=0.0 )

  return answer;
}//std::vector<float> fullrangefraction_coef_to_polynomial( const std::vector<float> &coeffs, const size_t nbin )


vector<float> mid_channel_polynomial_to_fullrangeFraction(
                               const vector<float> &coeffs, const size_t nbini )
{
  const float nbin = static_cast<float>( nbini );
  float a0 = 0.0f, a1 = 0.0f, a2 = 0.0f, a3 = 0.0f;
  float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;
  
  if( coeffs.size() > 0 )
    c0 = coeffs[0];
  if( coeffs.size() > 1 )
    c1 = coeffs[1];
  if( coeffs.size() > 2 )
    c2 = coeffs[2];
  if( coeffs.size() > 3 )
    c3 = coeffs[3];
  
  //Logic from from "Energy Calibration Conventions and Parameter Conversions"
  //  Dean J. Mitchell, George Lasche and John Mattingly, 6418 / MS-0791
  //The bellow "- 0.5f*c1" should actually be "+ 0.5f*c1" according to the
  //  paper above, but to agree with PeakEasy I changed this.
  a0 = c0 - 0.5f*c1 + 0.25f*c2 + (1.0f/8.0f)*c3;
  a1 = nbin*(c1 + c2 + 0.75f*c3);
  a2 = nbin*nbin*(c2 + 1.5f*c3);
  a3 = nbin*nbin*nbin*c3;
  
  vector<float> answer;
  answer.push_back( a0 );
  answer.push_back( a1 );
  if( a2 != 0.0f || a3 != 0.0f )
    answer.push_back( a2 );
  if( a3 != 0.0f )
    answer.push_back( a3 );
  
  return answer;
}//mid_channel_polynomial_to_fullrangeFraction()


vector<float> polynomial_coef_to_fullrangefraction( const vector<float> &coeffs,
                                                    const size_t nbin )
{
  float a0 = 0.0f, a1 = 0.0f, a2 = 0.0f, a3 = 0.0f;
  float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;

  if( coeffs.size() > 0 )
    c0 = coeffs[0];
  if( coeffs.size() > 1 )
    c1 = coeffs[1];
  if( coeffs.size() > 2 )
    c2 = coeffs[2];
  if( coeffs.size() > 3 )
    c3 = coeffs[3];

  //Logic from from "Energy Calibration Conventions and Parameter Conversions"
  //  Dean J. Mitchell, George Lasche and John Mattingly, 6418 / MS-0791
  //The bellow "- 0.5f*c1" should actually be "+ 0.5f*c1" according to the
  //  paper above, but to agree with PeakEasy I changed this.
//  a0 = c0 - 0.5f*c1 + 0.25f*c2 + (1.0f/8.0f)*c3;
//  a1 = nbin*(c1 + c2 + 0.75f*c3);
//  a2 = nbin*nbin*(c2 + 1.5f*c3);
//  a3 = nbin*nbin*nbin*c3;
  
  //In talking to Lee Harding, the bellow conversion is what should be used.
  //In talking to Dean, he said (in reference to the note referenced above):
  //  "The memo that you reference is confusing because neither Canberra nor
  //   Raytheon applied the definition correctly with their ASP systems (despite
  //   the fact that Bob Huckins of Canberra wrote the relevant part of the 2006
  //   N42 standard)."
  a0 = c0;
  a1 = nbin*c1;
  a2 = nbin*nbin*c2;
  a3 = nbin*nbin*nbin*c3;
  
  
  vector<float> answer;
  answer.push_back( a0 );
  answer.push_back( a1 );
  if( a2 != 0.0f || a3 != 0.0f )
    answer.push_back( a2 );
  if( a3 != 0.0f )
    answer.push_back( a3 );

  return answer;
}//vector<float> polynomial_coeef_to_fullrangefraction( const vector<float> &coeffs )


ShrdConstFVecPtr polynomial_binning( const vector<float> &coeffs,
                                     const size_t nbin,
                                     const DeviationPairVec &dev_pairs )
{
  ShrdFVecPtr answer( new vector<float>(nbin, 0.0) );
  const size_t ncoeffs = coeffs.size();
  
  for( size_t i = 0; i < nbin; i++ )
  {
    float val = 0.0;
    for( size_t c = 0; c < ncoeffs; ++c )
      val += coeffs[c] * pow( static_cast<float>(i), static_cast<float>(c) );
    answer->operator[](i) = val;
  }//for( loop over bins, i )
  
  return apply_deviation_pair( answer, dev_pairs );
  //  const vector<float> fff_coef = polynomial_coef_to_fullrangefraction( coeffs, nbin );
  //  return fullrangefraction_binning( fff_coef, nbin, dev_pairs );
}//ShrdConstFVecPtr polynomial_binning( const vector<float> &coefficients, size_t nbin )


ShrdConstFVecPtr fullrangefraction_binning( const vector<float> &coeffs,
                                            const size_t nbin,
                                            const DeviationPairVec &dev_pairs  )
{
  ShrdFVecPtr answer( new vector<float>(nbin, 0.0f) );
  const size_t ncoeffs = std::min( coeffs.size(), size_t(4) );
  const float low_e_coef = (coeffs.size() > 4) ? coeffs[4] : 0.0f;
  
  for( size_t i = 0; i < nbin; i++ )
  {
    const float x = static_cast<float>(i)/static_cast<float>(nbin);
    float &val = answer->operator[](i);
    for( size_t c = 0; c < ncoeffs; ++c )
      val += coeffs[c] * pow(x,static_cast<float>(c) );
    val += low_e_coef / (1.0f+60.0f*x);
  }//for( loop over bins, i )

  return apply_deviation_pair( answer, dev_pairs );
}//ShrdConstFVecPtr fullrangefraction_binning(...)


float find_bin_fullrangefraction( const double energy,
                                  const std::vector<float> &coeffs,
                                  const size_t nbin,
                                  const DeviationPairVec &devpair,
                                  const double accuracy )
{
  size_t ncoefs = 0;
  for( size_t i = 0; i < coeffs.size(); ++i )
    if( fabs(coeffs[i]) > std::numeric_limits<float>::epsilon() )
      ncoefs = i+1;

  if( ncoefs < 2  )
    throw std::runtime_error( "find_bin_fullrangefraction(...): must pass"
                              " in at least two coefficients" );
  
  if( ncoefs < 4 && devpair.empty() )
  {
    if( ncoefs == 2  )
    {
      //  energy =  coeffs[0] + coeffs[1]*(bin/nbins)
      return static_cast<float>( nbin * (energy - coeffs[0]) / coeffs[1] );
    }//if( coeffs.size() == 2  )
    
    //Note purposeful use of double precision
    
    const double a = double(coeffs[0]) - double(energy);
    const double b = coeffs[1];
    const double c = coeffs[2];
    
    //energy = coeffs[0] + coeffs[1]*(bin/nbin) + coeffs[2]*(bin/nbin)*(bin/nbin)
    //--> 0 = a + b*(bin/nbin) + c*(bin/nbin)*(bin/nbin)
    //roots at (-b +- sqrt(b*b-4*a*c))/(2c)
    
    const double sqrtarg = b*b-4.0f*a*c;
    
    if( sqrtarg >= 0.0 )
    {
      const double root_1 = (-b + sqrt(sqrtarg))/(2.0f*c);
      const double root_2 = (-b - sqrt(sqrtarg))/(2.0f*c);
    
      vector<double> roots;
      if( root_1 > 0.0 && root_1 < static_cast<double>(nbin) )
        roots.push_back( root_1 );
      if( root_2 > 0.0 && root_2 < static_cast<double>(nbin) )
        roots.push_back( root_2 );

      if( roots.size() == 1 )
        return static_cast<float>( nbin * roots[0] );
    
      if( roots.size() == 2 )
      {
        const double e1 = coeffs[0] + coeffs[1]*root_1 + coeffs[2]*root_1*root_1;
        const double e2 = coeffs[0] + coeffs[1]*root_2 + coeffs[2]*root_2*root_2;
        if( fabs(e1-e2) < static_cast<double>(accuracy) )
          return static_cast<float>( nbin * root_1 );
      }//if( roots.size() == 2 )

#if( BUILD_AS_UNIT_TEST_SUITE )
      //This is kinda a hack...
      const string testname = boost::unit_test::framework::current_test_case().p_name;
      if( testname.find( "FindEnergy") != string::npos )
        throw runtime_error( "Failed to find exact solution for bin number" );
#endif
    
#if( PERFORM_DEVELOPER_CHECKS )
      stringstream msg;
      msg << "find_bin_fullrangefraction(): found " << roots.size()
          << " energy solutions, shouldnt have happened: "
          << root_1 << " and " << root_2 << " for " << energy << " keV"
          << " so root coorespond to "
          << coeffs[0] + coeffs[1]*root_1 + coeffs[2]*root_1*root_1
          << " and "
          << coeffs[0] + coeffs[1]*root_2 + coeffs[2]*root_2*root_2
          << ". I will attempt to recover, but please check results of operation.";
//    passMessage( msg.str(), "", 3 );
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
#endif
      cerr << SRC_LOCATION << "\n\tWarning, couldnt algebraicly find bin number\n"
           << "\tcoeffs[0]=" << coeffs[0] << ", coeffs[1]=" << coeffs[1]
           << ", coeffs[2]=" << coeffs[2] << ", energy=" << energy << endl;
    }//if( sqrtarg >= 0.0f )
  }//if( ncoefs < 4 && devpair.empty() )
    
  
  float lowbin = 0.0;
  float highbin = static_cast<float>( nbin );
  float testenergy = fullrangefraction_energy( highbin, coeffs, nbin, devpair );
  while( testenergy < energy )
  {
    highbin *= 2.0f;
    testenergy = fullrangefraction_energy( highbin, coeffs, nbin, devpair );
  }//while( testenergy < energy )
  
  testenergy = fullrangefraction_energy( lowbin, coeffs, nbin, devpair );
  while( testenergy > energy )
  {
    lowbin -= nbin;
    testenergy = fullrangefraction_energy( lowbin, coeffs, nbin, devpair );
  }//while( testenergy < energy )
  
  
  float bin = lowbin + ((highbin-lowbin)/2.0f);
  testenergy = fullrangefraction_energy( bin, coeffs, nbin, devpair );
  float dx = static_cast<float>( fabs(testenergy-energy) );
  
  while( dx > accuracy )
  {
    if( highbin == lowbin )
    {
      cerr << "Possible error in find_bin_fullrangefraction... check out" << endl;
      throw runtime_error( "find_bin_fullrangefraction(...): error finding bin coorespongin to deired energy (this shouldnt happen)" );
    }

    if( testenergy == energy )
      return bin;
    if( testenergy > energy )
      highbin = bin;
    else
      lowbin = bin;

    bin = lowbin + ((highbin-lowbin)/2.0f);
    testenergy = fullrangefraction_energy( bin, coeffs, nbin, devpair );
    dx = static_cast<float>( fabs(testenergy-energy) );
  }//while( dx > accuracy )
  
  return bin;
}//float find_bin_fullrangefraction(...)


float fullrangefraction_energy( float bin_number,
                                const std::vector<float> &coeffs,
                                const size_t nbin,
                                const DeviationPairVec &deviation_pairs )
{
  const float x = bin_number/static_cast<float>(nbin);
  float val = 0.0;
  for( size_t c = 0; c < coeffs.size(); ++c )
    val += coeffs[c] * pow(x,static_cast<float>(c) );
  return val + offset_due_to_deviation_pairs( val, deviation_pairs );
}//float fullrangefraction_energy(...)


float offset_due_to_deviation_pairs( float energy, const DeviationPairVec &dps )
{
  if( dps.empty() )
    return 0.0;
  
  const size_t ndev_pair = dps.size() - 1;

  float answer = 0.0f;
  for( size_t dev_pairn = 0; dev_pairn < ndev_pair; ++dev_pairn )
  {
    const float lower_e = dps[dev_pairn].first;
    const float upper_e = dps[dev_pairn+1].first;
    const float lower_shift = dps[dev_pairn].second;
    const float upper_shift = dps[dev_pairn+1].second;
    const float residual_shift = (upper_shift -lower_shift)/(upper_e - lower_e);

    if( (energy >= lower_e) && (energy < upper_e) )
      answer += (energy - lower_e)*residual_shift;
    else if( energy >= upper_e )
      answer += (upper_shift - lower_shift);
  }//for( size_t dev_pairn = 0; dev_pairn < ndev_pair; ++dev_pairn )
  
  return answer;
}//float offset_due_to_deviation_pairs(...)


ShrdConstFVecPtr apply_deviation_pair( ShrdConstFVecPtr binning,
                             const std::vector< std::pair<float,float> > &dps )
{
  //This function adapted by wcjohns from a previous TRB sub-project, and has
  //  not been tested.  Further, it linearly interpolates between deviation
  //  pairs, causing bad artifacts in data; a spline interpolation should be
  //  used.
  //ToDo - XXX - convert to using a spline interpolation of deviation pairs

  if( dps.empty() )
    return binning;
  
  ShrdFVecPtr answer( new vector<float>( *binning ) );
  vector<float> &Ex = *answer;
  const vector<float> &Ex_org = *binning;
  const size_t nchannel = binning->size();
  const size_t ndev_pair = dps.size() - 1;

  // Offset prior to the first Deviation pair
  size_t ch = upper_bound( Ex.begin(), Ex.end(), dps[0].first ) - Ex.begin();

  for( size_t dev_pairn = 0; dev_pairn < ndev_pair; ++dev_pairn )
  {
    const float lower_e = dps[dev_pairn].first;
    const float upper_e = dps[dev_pairn+1].first;
    const float lower_shift = dps[dev_pairn].second;
    const float upper_shift = dps[dev_pairn+1].second;
    const float residual_shift = (upper_shift -lower_shift)/(upper_e - lower_e);

    for( size_t channel = ch; channel < nchannel; ++channel )
    {
      const float this_orig_e = Ex_org[channel];

      if( (this_orig_e >= lower_e) && (this_orig_e < upper_e) )
      {
        const float shift = (this_orig_e - lower_e)*residual_shift;
        Ex[channel] += shift;
        
        ch = channel;
      }//if( (this_orig_e >= lower_e) && (this_orig_e < upper_e) )
    }//for( size_t f = ch; f < nchannel; ++f )

    for( size_t channel = ch+1; channel < nchannel; ++channel )
      Ex[channel] += (upper_shift - lower_shift);
  }//for( size_t dev_pairn = 0; dev_pairn < ndev_pair; ++dev_pairn )
  
  return answer;
}//std::vector<float> apply_deviation_pair(...)



float bin_number_to_energy_polynomial( const float bin,
                                       const std::vector<float> &coeffs,
                                       const size_t nbin )
{
  if( coeffs.empty() )
    return 0.0;

  const int npars = static_cast<int>( coeffs.size() );

  if( npars > 3 )
    cerr << "Warning: bin_number_to_energy_polynomial(...) only uses first three "
            "coefficients right now, ignoring higher order terms" << endl;

  double a0 = 0.0, a1 = 0.0, a2 = 0.0;

  if( npars > 0)
  {
    a0 += coeffs[0];

    if( npars > 1 )
    {
      a0 += coeffs[1]/2.0 + coeffs[2]/4.0;
      a1 += nbin*coeffs[1];

      if( npars > 2 )
      {
        a0 += coeffs[2]/4.0;
        a1 += nbin*coeffs[2];
        a2 += nbin*nbin*coeffs[2];
      }//if( npars > 2 )
    }//if( npars > 1 )
  }//if( npars > 0)

  const float x = bin / static_cast<float>(nbin);
  return float( a0 + a1*x + a2*x*x );
}//float bin_number_to_energy_polynomial(...)


void expand_counted_zeros( const vector<float> &data, vector<float> &return_answer )
{
  vector<float> answer;
  answer.reserve( 1024 );
  vector<float>::const_iterator iter;
  for( iter = data.begin(); iter != data.end(); iter++)
  {
    if( (*iter != 0.0f) || (iter+1==data.end()) || (*(iter+1)==0.0f) )
      answer.push_back(*iter);
    else
    {
      iter++;
      const size_t nZeroes = ((iter==data.end()) ? 0u : static_cast<size_t>(floor(*iter + 0.5f)) );

      if( ((*iter) <= 0.5f) || ((answer.size() + nZeroes) > 131072) )
        throw runtime_error( "Invalid counted zeros: too many total elements, or negative number of zeros" );
      
      for( size_t k = 0; k < nZeroes; ++k )
        answer.push_back( 0.0f );
    }//if( at a non-zero value, the last value, or the next value is zero) / else
  }//for( iterate over data, iter )

  answer.swap( return_answer );
}//vector<float> expand_counted_zeros( const vector<float> &data )


void compress_to_counted_zeros( const std::vector<float> &input, std::vector<float> &results )
{
  results.clear();

  const size_t nBin = input.size();

  for( size_t bin = 0; bin < nBin; ++bin )
  {
      const bool isZero = (fabs(input[bin]) < 1E-8);

      if( !isZero ) results.push_back( input[bin] );
      else          results.push_back( 0.0f );

      if( isZero )
      {
          size_t nBinZeroes = 0;
          while( ( bin < nBin ) && ( fabs( input[bin] ) < 1E-8) )
          {
            ++nBinZeroes;
            ++bin;
          }//while more zero bins

          results.push_back( static_cast<float>(nBinZeroes) );

          if( bin != nBin )
            --bin;
      }//if( input[bin] == 0.0 )
  }//for( size_t bin = 0; bin < input.size(); ++bin )
}//void compress_to_counted_zeros(...)

#if( PERFORM_DEVELOPER_CHECKS )

namespace
{
#define compare_field(l,r,f) if( (l.f) != (r.f) ) return (l.f) < (r.f);
  bool compare_DetectorAnalysisResult( const DetectorAnalysisResult &lhs, const DetectorAnalysisResult &rhs )
  {
    compare_field(lhs,rhs,remark_);
    compare_field(lhs,rhs,nuclide_);
    compare_field(lhs,rhs,activity_);
    compare_field(lhs,rhs,nuclide_type_);
    compare_field(lhs,rhs,id_confidence_);
    compare_field(lhs,rhs,distance_);
    compare_field(lhs,rhs,dose_rate_);
    compare_field(lhs,rhs,real_time_);
    //compare_field(lhs,rhs,start_time_);
    compare_field(lhs,rhs,detector_);
    return false;
  };
}//namespace

void DetectorAnalysisResult::equalEnough( const DetectorAnalysisResult &lhs,
                                          const DetectorAnalysisResult &rhs )
{
  if( lhs.remark_ != rhs.remark_ )
    throw runtime_error( "DetectorAnalysisResult remark for LHS ('"
                         + lhs.remark_ + "') doesnt match RHS ('"
                         + rhs.remark_ + "')" );
  if( lhs.nuclide_ != rhs.nuclide_ )
    throw runtime_error( "DetectorAnalysisResult nuclide for LHS ('"
                         + lhs.nuclide_ + "') doesnt match RHS ('"
                         + rhs.nuclide_ + "')" );
  if( lhs.nuclide_type_ != rhs.nuclide_type_ )
    throw runtime_error( "DetectorAnalysisResult nuclide type for LHS ('"
                         + lhs.nuclide_type_ + "') doesnt match RHS ('"
                         + rhs.nuclide_type_ + "')" );
  if( lhs.id_confidence_ != rhs.id_confidence_ )
    throw runtime_error( "DetectorAnalysisResult ID confifence for LHS ('"
                         + lhs.id_confidence_ + "') doesnt match RHS ('"
                         + rhs.id_confidence_ + "')" );
  if( lhs.detector_ != rhs.detector_ )
    throw runtime_error( "DetectorAnalysisResult detector for LHS ('"
                         + lhs.detector_ + "') doesnt match RHS ('"
                         + rhs.detector_ + "')" );
  //if( lhs.start_time_ != rhs.start_time_ )
  //  throw runtime_error( "DetectorAnalysisResult start time for LHS ("
  //                      + UtilityFunctions::to_iso_string(lhs.start_time_) + ") doesnt match RHS ("
  //                    + UtilityFunctions::to_iso_string(rhs.start_time_) + ")" );
  
  char buffer[1024];
  if( fabs(lhs.activity_ - rhs.activity_) > 0.00001 )
  {
    snprintf( buffer, sizeof(buffer),
             "DetectorAnalysisResult activity of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.activity_, rhs.activity_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.distance_ - rhs.distance_) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer),
             "DetectorAnalysisResult distance of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.activity_, rhs.activity_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.dose_rate_ - rhs.dose_rate_) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer),
             "DetectorAnalysisResult dose rate of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.dose_rate_, rhs.dose_rate_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.real_time_ - rhs.real_time_) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer),
             "DetectorAnalysisResult dose rate of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.real_time_, rhs.real_time_ );
    throw runtime_error( buffer );
  }
}//void DetectorAnalysisResult::equalEnough(...)


void DetectorAnalysis::equalEnough( const DetectorAnalysis &lhs,
                                    const DetectorAnalysis &rhs )
{
  char buffer[1024];
  
  vector<string> lhsremarks, rhsremarks;
  
  for( string remark : lhs.remarks_ )
  {
//    remark.erase( std::remove_if(remark.begin(), remark.end(), not_alpha_numeric), remark.end());
    trim( remark );
    ireplace_all( remark, "  ", " " );
    if( remark.size() )
      lhsremarks.push_back( remark );
  }
  
  for( string remark : rhs.remarks_ )
  {
    //    remark.erase( std::remove_if(remark.begin(), remark.end(), not_alpha_numeric), remark.end());
    trim( remark );
    ireplace_all( remark, "  ", " " );
    if( remark.size() )
      rhsremarks.push_back( remark );
  }
  
  stable_sort( lhsremarks.begin(), lhsremarks.end() );
  stable_sort( rhsremarks.begin(), rhsremarks.end() );
  
  
  if( lhsremarks.size() != rhsremarks.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "Number of Analysis remarks for LHS (%i) doesnt match RHS %i",
             int(lhsremarks.size()), int(rhsremarks.size()) );
    throw runtime_error( buffer );
  }
  
  for( size_t i = 0; i < rhsremarks.size(); ++i )
  {
    if( lhsremarks[i] != rhsremarks[i] )
    {
      snprintf( buffer, sizeof(buffer),
               "Analysis remark %i for LHS ('%s') doesnt match RHS ('%s')",
               int(i), lhsremarks[i].c_str(), rhsremarks[i].c_str() );
      throw runtime_error( buffer );
    }
  }//for( size_t i = 0; i < rhsremarks.size(); ++i )
  
  if( lhs.algorithm_version_ != rhs.algorithm_version_ )
    throw runtime_error( "Analysis algorithm version for LHS ('"
                        + lhs.algorithm_version_ + "') doesnt match RHS ('"
                        + rhs.algorithm_version_ + "')" );
  
  if( lhs.algorithm_name_ != rhs.algorithm_name_ )
    throw runtime_error( "Analysis algorithm name for LHS ('"
                        + lhs.algorithm_name_ + "') doesnt match RHS ('"
                        + rhs.algorithm_name_ + "')" );
  
  if( lhs.algorithm_creator_ != rhs.algorithm_creator_ )
    throw runtime_error( "Analysis algorithm creator for LHS ('"
                        + lhs.algorithm_creator_ + "') doesnt match RHS ('"
                        + rhs.algorithm_creator_ + "')" );
  
  if( lhs.algorithm_description_ != rhs.algorithm_description_ )
    throw runtime_error( "Analysis algorithm description for LHS ('"
                        + lhs.algorithm_description_ + "') doesnt match RHS ('"
                        + rhs.algorithm_description_ + "')" );

  if( lhs.algorithm_result_description_ != rhs.algorithm_result_description_ )
    throw runtime_error( "Analysis algorithm result description for LHS ('"
                        + lhs.algorithm_result_description_ + "') doesnt match RHS ('"
                        + rhs.algorithm_result_description_ + "')" );
  
  
  
  if( lhs.results_.size() != rhs.results_.size() )
  {
    stringstream msg;
    msg << "Differnt number of analysis results for LHS ("
    << lhs.results_.size() << ") vs RHS (" << rhs.results_.size() << "):"
    << endl;
    
    for( DetectorAnalysisResult l : lhs.results_ )
    {
      msg << "\tLHS: remark='" << l.remark_ << "', nuclide='"
           << l.nuclide_ << "', doserate="
           << l.dose_rate_ << ", activity=" << l.activity_
           << ", id confidence='" << l.id_confidence_ << "'"
           << ", distance=" << l.distance_ << endl;
    }
    
    for( DetectorAnalysisResult l : rhs.results_ )
    {
      msg << "\t RHS: remark='" << l.remark_ << "', nuclide='"
          << l.nuclide_ << "', doserate="
          << l.dose_rate_ << ", activity=" << l.activity_
          << ", id confidence='" << l.id_confidence_ << "'"
          << ", distance=" << l.distance_ << endl;
    }
    
    throw runtime_error( msg.str() );
  }//if( lhs.results_.size() != rhs.results_.size() )

  //Ordering of the results is not garunteed in a round trip, since we use a
  //  single analysis result type, but N42 2012 defines a few different ones.
  vector<DetectorAnalysisResult> rhsres = rhs.results_;
  vector<DetectorAnalysisResult> lhsres = lhs.results_;
  std::sort( lhsres.begin(), lhsres.end(), &compare_DetectorAnalysisResult );
  std::sort( rhsres.begin(), rhsres.end(), &compare_DetectorAnalysisResult );
  
  for( size_t i = 0; i < rhsres.size(); ++i )
    DetectorAnalysisResult::equalEnough( lhsres[i], rhsres[i] );
}//void DetectorAnalysis::equalEnough(...)



void Measurement::equalEnough( const Measurement &lhs, const Measurement &rhs )
{
  char buffer[1024];
  if( fabs(lhs.live_time_ - rhs.live_time_) > 0.00001 )
  {
    snprintf( buffer, sizeof(buffer),
             "Live time of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.live_time_, rhs.live_time_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.real_time_ - rhs.real_time_) > 0.00001 )
  {
    snprintf( buffer, sizeof(buffer),
             "Real time of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.real_time_, rhs.real_time_ );
    throw runtime_error( buffer );
  }
  
  if( lhs.contained_neutron_ != rhs.contained_neutron_ )
  {
    snprintf( buffer, sizeof(buffer),
             "LHS %s contain neutrons while RHS %s",
             (lhs.contained_neutron_?"did":"did not"),
             (rhs.contained_neutron_?"did":"did not") );
    throw runtime_error( buffer );
  }

  if( lhs.sample_number_ != rhs.sample_number_ )
    throw runtime_error( "LHS sample number some how didnt equal RHS sample number" );
  
  if( lhs.occupied_ != rhs.occupied_ )
  {
    snprintf( buffer, sizeof(buffer),
             "Ocupied of LHS (%i) differend form RHS (%i)",
             int(lhs.occupied_), int(rhs.occupied_) );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.gamma_count_sum_ - rhs.gamma_count_sum_) > (0.00001*std::max(fabs(lhs.gamma_count_sum_),fabs(rhs.gamma_count_sum_))) )
  {
    snprintf( buffer, sizeof(buffer),
             "Gamma count sum of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.gamma_count_sum_, rhs.gamma_count_sum_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.neutron_counts_sum_ - rhs.neutron_counts_sum_) > (0.00001*std::max(fabs(lhs.neutron_counts_sum_),fabs(rhs.neutron_counts_sum_))) )
  {
    snprintf( buffer, sizeof(buffer),
             "Neutron count sum of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.neutron_counts_sum_, rhs.neutron_counts_sum_ );
    throw runtime_error( buffer );
  }

  if( fabs(lhs.speed_ - rhs.speed_) > 0.01 )
  {
    snprintf( buffer, sizeof(buffer),
             "Speed of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.speed_, rhs.speed_ );
    throw runtime_error( buffer );
  }

  if( lhs.detector_name_ != rhs.detector_name_ )
    throw runtime_error( "Detector name for LHS ('" + lhs.detector_name_
                        + "') doesnt match RHS ('" + rhs.detector_name_ + "')" );

  if( lhs.detector_number_ != rhs.detector_number_ )
  {
    snprintf( buffer, sizeof(buffer),
             "Detector number of LHS (%i) doesnt match RHS (%i)",
             lhs.detector_number_, rhs.detector_number_ );
    throw runtime_error( buffer );
  }

  if( lhs.detector_type_ != rhs.detector_type_
     && rhs.detector_type_!="Gamma and Neutron"
     && ("NaI, " + lhs.detector_type_) != rhs.detector_type_
     && ("LaBr3, " + lhs.detector_type_) != rhs.detector_type_
     && ("unknown, " + lhs.detector_type_) != rhs.detector_type_
     )
    throw runtime_error( "Detector type for LHS ('" + lhs.detector_type_
                        + "') doesnt match RHS ('" + rhs.detector_type_ + "')" );

  if( lhs.quality_status_ != rhs.quality_status_ )
  {
    snprintf( buffer, sizeof(buffer),
             "Quality status of LHS (%i) different from RHS (%i)",
             int(lhs.quality_status_), int(rhs.quality_status_) );
    throw runtime_error( buffer );
  }
  
  if( lhs.source_type_ != rhs.source_type_ )
  {
    snprintf( buffer, sizeof(buffer),
             "Source type of LHS (%i) different from RHS (%i)",
             int(lhs.source_type_), int(rhs.source_type_) );
    throw runtime_error( buffer );
  }


  if( lhs.energy_calibration_model_ != rhs.energy_calibration_model_
      && !!lhs.channel_energies_ )
  {
    if( (lhs.energy_calibration_model_ == Polynomial
         && rhs.energy_calibration_model_ == FullRangeFraction)
       || (lhs.energy_calibration_model_ == FullRangeFraction
           && rhs.energy_calibration_model_ == Polynomial) )
    {
      const size_t nbin = lhs.channel_energies_->size();
      vector<float> lhscoefs = lhs.calibration_coeffs_;
      vector<float> rhscoefs = rhs.calibration_coeffs_;
      
      if( rhs.energy_calibration_model_ == Polynomial )
        rhscoefs = polynomial_coef_to_fullrangefraction( rhscoefs, nbin );
      if( lhs.energy_calibration_model_ == Polynomial )
        lhscoefs = polynomial_coef_to_fullrangefraction( lhscoefs, nbin );
      
      if( rhscoefs.size() != lhscoefs.size() )
        throw runtime_error( "Calibration coefficients LHS and RHS do not match size after converting to be same type" );
      
      for( size_t i = 0; i < rhscoefs.size(); ++i )
      {
        const float a = lhscoefs[i];
        const float b = rhscoefs[i];
        
        if( fabs(a - b) > (1.0E-5*std::max(fabs(a),fabs(b))) )
        {
          snprintf( buffer, sizeof(buffer),
                   "Calibration coefficient %i of LHS (%1.8E) doesnt match RHS (%1.8E) (note, concerted calib type)",
                   int(i), a, b );
          throw runtime_error( buffer );
        }
      }
    }else if( lhs.gamma_count_sum_>0.0 || rhs.gamma_count_sum_>0.0 )
    {
      snprintf( buffer, sizeof(buffer),
               "Calibration model of LHS (%i) different from RHS (%i)",
               int(lhs.energy_calibration_model_),
               int(rhs.energy_calibration_model_) );
      throw runtime_error( buffer );
    }
  }
  
  
  const set<string> nlhsremarkss( lhs.remarks_.begin(), lhs.remarks_.end() );
  const set<string> nrhsremarkss( rhs.remarks_.begin(), rhs.remarks_.end() );
  
  const vector<string> nlhsremarks( nlhsremarkss.begin(), nlhsremarkss.end() );
  const vector<string> nrhsremarks( nrhsremarkss.begin(), nrhsremarkss.end() );
  
  if( nlhsremarks.size() != nrhsremarks.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "Number of remarks in LHS (%i) doesnt match RHS (%i)",
             int(nlhsremarks.size()), int(nrhsremarks.size()) );
    
//    for( size_t i = 0; i < nlhsremarks.size(); ++i )
//      cerr << "LHS: '" << nlhsremarks[i] << "'" << endl;
//    for( size_t i = 0; i < nrhsremarks.size(); ++i )
//      cerr << "RHS: '" << nrhsremarks[i] << "'" << endl;
    
    throw runtime_error( buffer );
  }
  
  for( size_t i = 0; i < nlhsremarks.size(); ++i )
  {
    if( nlhsremarks[i] != nrhsremarks[i] )
    {
      snprintf( buffer, sizeof(buffer),
               "Remark %i in LHS ('%s') doesnt match RHS ('%s')",
               int(i), nlhsremarks[i].c_str(), nrhsremarks[i].c_str() );
      throw runtime_error( buffer );
    }
  }

  if( lhs.start_time_ != rhs.start_time_ )
    throw runtime_error( "Start time for LHS ("
                        + UtilityFunctions::to_iso_string(lhs.start_time_) + ") doesnt match RHS ("
                        + UtilityFunctions::to_iso_string(rhs.start_time_) + ")" );
  
  if( lhs.energy_calibration_model_ == rhs.energy_calibration_model_ )
  {
    vector<float> lhscalcoef = lhs.calibration_coeffs_;
    vector<float> rhscalcoef = rhs.calibration_coeffs_;
    
    while( lhscalcoef.size() && lhscalcoef.back() == 0.0f )
      lhscalcoef.erase( lhscalcoef.end() - 1 );
    while( rhscalcoef.size() && rhscalcoef.back() == 0.0f )
      rhscalcoef.erase( rhscalcoef.end() - 1 );

    if( lhscalcoef.size() != rhscalcoef.size() )
    {
      snprintf( buffer, sizeof(buffer),
               "Number of calibration coefficients of LHS (%i) doesnt match RHS (%i)",
               int(lhscalcoef.size()),
               int(rhscalcoef.size()) );
      throw runtime_error( buffer );
    }
  
    for( size_t i = 0; i < rhscalcoef.size(); ++i )
    {
      const float a = lhscalcoef[i];
      const float b = rhscalcoef[i];
    
      if( fabs(a - b) > (1.0E-5*std::max(fabs(a),fabs(b))) )
      {
        snprintf( buffer, sizeof(buffer),
                 "Calibration coefficient %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
                 int(i), a, b );
        throw runtime_error( buffer );
      }
    }
  }//if( energy_calibration_model_ == rhs.energy_calibration_model_ )
  
  if( lhs.deviation_pairs_.size() != rhs.deviation_pairs_.size() )
    throw runtime_error( "Number of deviation pairs of LHS and RHS dont match" );
  
  for( size_t i = 0; i < lhs.deviation_pairs_.size(); ++i )
  {
    if( fabs(lhs.deviation_pairs_[i].first - rhs.deviation_pairs_[i].first) > 0.001 )
    {
      snprintf( buffer, sizeof(buffer),
               "Energy of deviation pair %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
               int(i), lhs.deviation_pairs_[i].first,
               rhs.deviation_pairs_[i].first );
      throw runtime_error( buffer );
    }
    
    if( fabs(lhs.deviation_pairs_[i].second - rhs.deviation_pairs_[i].second) > 0.001 )
    {
      snprintf( buffer, sizeof(buffer),
               "Offset of deviation pair %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
               int(i), lhs.deviation_pairs_[i].second,
               rhs.deviation_pairs_[i].second );
      throw runtime_error( buffer );
    }
  }//for( size_t i = 0; i < lhs.deviation_pairs_.size(); ++i )
  
  if( (!lhs.channel_energies_) != (!rhs.channel_energies_) )
  {
    snprintf( buffer, sizeof(buffer), "Channel energies avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.channel_energies_?"missing":"available"),
             (!rhs.channel_energies_?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  if( !!lhs.channel_energies_ )
  {
    if( lhs.channel_energies_->size() != rhs.channel_energies_->size() )
    {
      snprintf( buffer, sizeof(buffer),
               "Number of channel energies of LHS (%i) doesnt match RHS (%i)",
               int(lhs.channel_energies_->size()),
               int(rhs.channel_energies_->size()) );
      throw runtime_error( buffer );
    }
    
    for( size_t i = 0; i < lhs.channel_energies_->size(); ++i )
    {
      const float lhsenergy = lhs.channel_energies_->at(i);
      const float rhsenergy = rhs.channel_energies_->at(i);
      const float diff = fabs(lhsenergy - rhsenergy );
      if( diff > 0.00001*std::max(fabs(lhsenergy), fabs(rhsenergy)) && diff > 0.00001 )
      {
        snprintf( buffer, sizeof(buffer),
                 "Energy of channel %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
                 int(i), lhs.channel_energies_->at(i),
                 rhs.channel_energies_->at(i) );
        throw runtime_error( buffer );
      }
    }
  }//if( !!channel_energies_ )
  
  if( (!lhs.gamma_counts_) != (!rhs.gamma_counts_) )
  {
    snprintf( buffer, sizeof(buffer), "Gamma counts avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.gamma_counts_?"missing":"available"),
             (!rhs.gamma_counts_?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  if( !!lhs.gamma_counts_ )
  {
    if( lhs.gamma_counts_->size() != rhs.gamma_counts_->size() )
    {
      snprintf( buffer, sizeof(buffer),
               "Number of gamma channels of LHS (%i) doesnt match RHS (%i)",
               int(lhs.gamma_counts_->size()),
               int(rhs.gamma_counts_->size()) );
      throw runtime_error( buffer );
    }
    
    for( size_t i = 0; i < lhs.gamma_counts_->size(); ++i )
    {
      const float a = lhs.gamma_counts_->at(i);
      const float b = rhs.gamma_counts_->at(i);
      if( fabs(a-b) > 1.0E-6*std::max( fabs(a), fabs(b) ) )
      {
        cerr << "LHS:";
        for( size_t j = (i>4) ? i-4: 0; j < lhs.gamma_counts_->size() && j < (i+5); ++j )
        {
          if( i==j ) cerr << "__";
          cerr << lhs.gamma_counts_->at(j);
          if( i==j ) cerr << "__, ";
          else cerr << ", ";
        }
        cerr << endl;

        cerr << "RHS:";
        for( size_t j = (i>4) ? i-4: 0; j < rhs.gamma_counts_->size() && j < (i+5); ++j )
        {
          if( i==j ) cerr << "__";
          cerr << rhs.gamma_counts_->at(j);
          if( i==j ) cerr << "__, ";
          else cerr << ", ";
        }
        cerr << endl;
        
        snprintf( buffer, sizeof(buffer),
                 "Counts in gamma channel %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
                 int(i), lhs.gamma_counts_->at(i),
                 rhs.gamma_counts_->at(i) );
        throw runtime_error( buffer );
      }
    }
  }//if( !!channel_energies_ )
  
  if( lhs.neutron_counts_.size() != rhs.neutron_counts_.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "Number of neutron channels of LHS (%i) doesnt match RHS (%i)",
             int(lhs.neutron_counts_.size()),
             int(rhs.neutron_counts_.size()) );
    throw runtime_error( buffer );
  }
  
  for( size_t i = 0; i < lhs.neutron_counts_.size(); ++i )
  {
    if( fabs(lhs.neutron_counts_[i] - rhs.neutron_counts_[i] ) > 0.0001 )
    {
      snprintf( buffer, sizeof(buffer),
               "Counts in neutron channel %i of LHS (%1.8E) doesnt match RHS (%1.8E)",
               int(i), lhs.neutron_counts_[i],
               rhs.neutron_counts_[i] );
      throw runtime_error( buffer );
    }
  }
  
  if( fabs(lhs.latitude_ - rhs.latitude_) > 0.00001 )
  {
    snprintf( buffer, sizeof(buffer),
             "Latitude of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.latitude_, rhs.latitude_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.longitude_ - rhs.longitude_) > 0.00001 )
  {
    snprintf( buffer, sizeof(buffer),
             "Longitude of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.longitude_, rhs.longitude_ );
    throw runtime_error( buffer );
  }

  if( lhs.position_time_ != rhs.position_time_ )
    throw runtime_error( "Position time for LHS ("
                        + UtilityFunctions::to_iso_string(lhs.position_time_) + ") doesnt match RHS ("
                        + UtilityFunctions::to_iso_string(rhs.position_time_) + ")" );

  if( lhs.title_ != rhs.title_ )
    throw runtime_error( "Title for LHS ('" + lhs.title_
                        + "') doesnt match RHS ('" + rhs.title_ + "')" );
}//void equalEnough( const Measurement &lhs, const Measurement &rhs )


void MeasurementInfo::equalEnough( const MeasurementInfo &lhs,
                                   const MeasurementInfo &rhs )
{
  std::lock( lhs.mutex_, rhs.mutex_ );
  std::unique_lock<std::recursive_mutex> lhs_lock( lhs.mutex_, std::adopt_lock_t() );
  std::unique_lock<std::recursive_mutex> rhs_lock( rhs.mutex_, std::adopt_lock_t() );

  char buffer[1024];
  if( fabs(lhs.gamma_live_time_ - rhs.gamma_live_time_) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer),
              "MeasurementInfo: Live time of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.gamma_live_time_, rhs.gamma_live_time_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.gamma_real_time_ - rhs.gamma_real_time_) > 0.001 )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Real time of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.gamma_real_time_, rhs.gamma_real_time_ );
    throw runtime_error( buffer );
  }
  
  const double gamma_sum_diff = fabs(lhs.gamma_count_sum_ - rhs.gamma_count_sum_);
  const double max_gamma_sum = std::max(fabs(lhs.gamma_count_sum_), fabs(rhs.gamma_count_sum_));
  if( gamma_sum_diff > 0.1 || gamma_sum_diff > 1.0E-6*max_gamma_sum )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Gamma sum of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.gamma_count_sum_, rhs.gamma_count_sum_ );
    throw runtime_error( buffer );
  }
  
  if( fabs(lhs.neutron_counts_sum_ - rhs.neutron_counts_sum_) > 0.01 )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Neutron sum of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.neutron_counts_sum_, rhs.neutron_counts_sum_ );
    throw runtime_error( buffer );
  }
  
  if( lhs.filename_ != rhs.filename_ )
    throw runtime_error( "MeasurementInfo: Filename of LHS (" + lhs.filename_
                         + ") doenst match RHS (" + rhs.filename_ + ")" );
  
  if( lhs.detector_names_.size() != rhs.detector_names_.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Number of detector names of LHS (%i) doesnt match RHS (%i)",
             int(lhs.detector_names_.size()), int(rhs.detector_names_.size()) );
    throw runtime_error( buffer );
  }
  
  const set<string> lhsnames( lhs.detector_names_.begin(),
                              lhs.detector_names_.end() );
  const set<string> rhsnames( rhs.detector_names_.begin(),
                              rhs.detector_names_.end() );
  
  if( lhsnames != rhsnames )
    throw runtime_error( "MeasurementInfo: Detector names do not match for LHS and RHS" );
 
  if( lhs.detector_numbers_.size() != rhs.detector_numbers_.size()
      || lhs.detector_numbers_.size() != lhs.detector_names_.size() )
    throw runtime_error( "MeasurementInfo: Inproper number of detector numbers - wtf" );
  
  for( size_t i = 0; i < lhs.detector_names_.size(); ++i )
  {
    const size_t pos = std::find( rhs.detector_names_.begin(),
                                  rhs.detector_names_.end(),
                                  lhs.detector_names_[i] )
                       - rhs.detector_names_.begin();
    if( lhs.detector_numbers_[i] != rhs.detector_numbers_[pos] )
      throw runtime_error( "MeasurementInfo: Detector number for detector '"
                           + lhs.detector_names_[i] + "' dont match" );
  }//for( size_t i = 0; i < lhs.detector_names_.size(); ++i )
  
  if( lhs.neutron_detector_names_.size() != rhs.neutron_detector_names_.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Number of neutron detector names of LHS (%i) doesnt match RHS (%i)",
             int(lhs.neutron_detector_names_.size()),
             int(rhs.neutron_detector_names_.size()) );
    throw runtime_error( buffer );
  }
  
  const set<string> nlhsnames( lhs.neutron_detector_names_.begin(),
                               lhs.neutron_detector_names_.end() );
  const set<string> nrhsnames( rhs.neutron_detector_names_.begin(),
                               rhs.neutron_detector_names_.end() );
  
  if( nlhsnames != nrhsnames )
    throw runtime_error( "MeasurementInfo: Neutron detector names dont match for LHS and RHS" );

  
  if( lhs.lane_number_ != rhs.lane_number_ )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Lane number of LHS (%i) doesnt match RHS (%i)",
             lhs.lane_number_, rhs.lane_number_ );
    throw runtime_error( buffer );
  }
  
  if( lhs.measurement_location_name_ != rhs.measurement_location_name_ )
    throw runtime_error( "MeasurementInfo: Measurment location name of LHS ('"
                         + lhs.measurement_location_name_
                         + "') doesnt match RHS ('"
                         + rhs.measurement_location_name_ + "')" );
  
  if( lhs.inspection_ != rhs.inspection_ )
    throw runtime_error( "MeasurementInfo: Inspection of LHS ('" + lhs.inspection_
                        + "') doesnt match RHS ('" + rhs.inspection_ + "')" );

  string leftoperator = lhs.measurment_operator_;
  string rightoperator = rhs.measurment_operator_;
  ireplace_all( leftoperator, "\t", " " );
  ireplace_all( leftoperator, "  ", " " );
  trim( leftoperator );
  ireplace_all( rightoperator, "\t", " " );
  ireplace_all( rightoperator, "  ", " " );
  trim( rightoperator );
  
  if( leftoperator != rightoperator )
    throw runtime_error( "MeasurementInfo: Measurment operator of LHS ('"
                         + lhs.measurment_operator_ + "') doesnt match RHS ('"
                         + rhs.measurment_operator_ + ")" );

  if( lhs.sample_numbers_.size() != rhs.sample_numbers_.size() )
  {
    cout << "lhs.measurements_.size()=" << lhs.measurements_.size() << endl;
    cout << "rhs.measurements_.size()=" << rhs.measurements_.size() << endl;
    
    for( size_t i = 0; i < lhs.measurements_.size(); ++i )
    {
      cout << "LHS: DetName=" << lhs.measurements_[i]->detector_name_
          << " {" << lhs.measurements_[i]->sample_number_
          << "," << lhs.measurements_[i]->detector_number_ << "}"
           << ", SumGamma=" << lhs.measurements_[i]->gamma_count_sum_
           << ", SumNeutron=" << lhs.measurements_[i]->neutron_counts_sum_
           << ", LiveTime=" << lhs.measurements_[i]->live_time_
           << ", RealTime=" << lhs.measurements_[i]->real_time_
           << ", StarTime=" << lhs.measurements_[i]->start_time_
      << endl;
    }
    
    for( size_t i = 0; i < rhs.measurements_.size(); ++i )
    {
      cout << "RHS: DetName=" << rhs.measurements_[i]->detector_name_
      << " {" << rhs.measurements_[i]->sample_number_
               << "," << rhs.measurements_[i]->detector_number_ << "}"
      << ", SumGamma=" << rhs.measurements_[i]->gamma_count_sum_
      << ", SumNeutron=" << rhs.measurements_[i]->neutron_counts_sum_
      << ", LiveTime=" << rhs.measurements_[i]->live_time_
      << ", RealTime=" << rhs.measurements_[i]->real_time_
      << ", StarTime=" << rhs.measurements_[i]->start_time_
      << endl;
    }
    
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Number of sample numbers in LHS (%i) doesnt match RHS (%i)",
             int(lhs.sample_numbers_.size()), int(rhs.sample_numbers_.size()) );
    throw runtime_error( buffer );
  }
  
  if( lhs.sample_numbers_ != rhs.sample_numbers_ )
  {
    stringstream lhssamples, rhssamples;
    for( auto sample : lhs.sample_numbers_ )
      lhssamples << (sample==(*lhs.sample_numbers_.begin()) ? "":",") << sample;
    for( auto sample : rhs.sample_numbers_ )
      rhssamples << (sample== (*rhs.sample_numbers_.begin()) ? "":",") << sample;
    throw runtime_error( "MeasurementInfo: Sample numbers of RHS (" + rhssamples.str()
                         + ") and LHS (" + lhssamples.str() + ") doent match" );
  }
  
  if( lhs.detector_type_ != rhs.detector_type_ )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: LHS detector type (%i) doesnt match RHS (%i)",
             int(lhs.detector_type_), int(rhs.detector_type_) );
    throw runtime_error( buffer );
  }
  
  string lhsinst = convert_n42_instrument_type_from_2006_to_2012( lhs.instrument_type_ );
  string rhsinst = convert_n42_instrument_type_from_2006_to_2012( rhs.instrument_type_ );
  if( lhsinst != rhsinst )
  {
    throw runtime_error( "MeasurementInfo: Instrument type of LHS ('" + lhs.instrument_type_
                    + "') doesnt match RHS ('" + rhs.instrument_type_ + "')" );
  }
  
  if( lhs.manufacturer_ != rhs.manufacturer_ )
    throw runtime_error( "MeasurementInfo: Manufacturer of LHS ('" + lhs.manufacturer_
                        + "') doesnt match RHS ('" + rhs.manufacturer_ + "')" );
  
  if( lhs.instrument_model_ != rhs.instrument_model_ )
    throw runtime_error( "MeasurementInfo: Instrument model of LHS ('" + lhs.instrument_model_
                    + "') doesnt match RHS ('" + rhs.instrument_model_ + "')" );
  
  if( lhs.instrument_id_ != rhs.instrument_id_ )
    throw runtime_error( "MeasurementInfo: Instrument ID model of LHS ('" + lhs.instrument_id_
                       + "') doesnt match RHS ('" + rhs.instrument_id_ + "')" );
  
  
  if( fabs(lhs.mean_latitude_ - rhs.mean_latitude_) > 0.000001 )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Mean latitude of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.mean_latitude_, rhs.mean_latitude_ );
    throw runtime_error( buffer );
  }

  if( fabs(lhs.mean_longitude_ - rhs.mean_longitude_) > 0.000001 )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Mean longitude of LHS (%1.8E) doesnt match RHS (%1.8E)",
             lhs.mean_longitude_, rhs.mean_longitude_ );
    throw runtime_error( buffer );
  }

  if( lhs.properties_flags_ != rhs.properties_flags_ )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Properties flags of LHS (%x) doesnt match RHS (%x)",
             static_cast<unsigned int>(lhs.properties_flags_),
             static_cast<unsigned int>(rhs.properties_flags_) );
    throw runtime_error( buffer );
  
  }
  
  for( const int sample : lhs.sample_numbers_ )
  {
    for( const int detnum : lhs.detector_numbers_ )
    {
      MeasurementConstShrdPtr lhsptr = lhs.measurement( sample, detnum );
      MeasurementConstShrdPtr rhsptr = rhs.measurement( sample, detnum );
      
      if( (!lhsptr) != (!rhsptr) )
      {
        snprintf( buffer, sizeof(buffer), "MeasurementInfo: Measurment avaialblity for LHS (%s)"
                 " doesnt match RHS (%s) for sample %i and detector number %i",
              (!lhsptr?"missing":"available"), (!rhsptr?"missing":"available"),
              sample, detnum  );
        throw runtime_error(buffer);
      }
      
      if( !lhsptr )
        continue;
      
      try
      {
        Measurement::equalEnough( *lhsptr, *rhsptr );
      }catch( std::exception &e )
      {
        snprintf( buffer, sizeof(buffer), "MeasurementInfo: Sample %i, Detector num %i: %s",
                  sample, detnum, e.what() );
        throw runtime_error( buffer );
      }
    }//for( const int detnum : lhs.detector_numbers_ )
  }//for( const int sample : lhs.sample_numbers_ )
  
  if( (!lhs.detectors_analysis_) != (!rhs.detectors_analysis_) )
  {
    snprintf( buffer, sizeof(buffer), "MeasurementInfo: Detector analysis avaialblity for LHS (%s)"
             " doesnt match RHS (%s)",
             (!lhs.detectors_analysis_?"missing":"available"),
             (!rhs.detectors_analysis_?"missing":"available") );
    throw runtime_error(buffer);
  }
  
  vector<string> nlhsremarkss, nrhsremarkss;
  for( string r : lhs.remarks_ )
  {
    while( r.find("  ") != string::npos )
      ireplace_all( r, "  ", " " );
    
    if( !starts_with(r,"N42 file created by")
        && !starts_with(r,"N42 file created by") )
      nlhsremarkss.push_back( r );
  }
  
  for( string r : rhs.remarks_ )
  {
    while( r.find("  ") != string::npos )
      ireplace_all( r, "  ", " " );
    
    if( !starts_with(r,"N42 file created by") )
      nrhsremarkss.push_back( r );
  }
  
  stable_sort( nlhsremarkss.begin(), nlhsremarkss.end() );
  stable_sort( nrhsremarkss.begin(), nrhsremarkss.end() );
  
  if( nlhsremarkss.size() != nrhsremarkss.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Number of remarks in LHS (%i) doesnt match RHS (%i)",
             int(nlhsremarkss.size()), int(nrhsremarkss.size()) );
    
    for( string a : nlhsremarkss )
    cout << "\tLHS: " << a << endl;
    for( string a : nrhsremarkss )
    cout << "\tRHS: " << a << endl;
    
    throw runtime_error( buffer );
  }
  
  for( size_t i = 0; i < nlhsremarkss.size(); ++i )
  {
    string lhsremark = nlhsremarkss[i];
    string rhsremark = nrhsremarkss[i];
    UtilityFunctions::trim( lhsremark );
    UtilityFunctions::trim( rhsremark );
    
    if( lhsremark != rhsremark )
    {
      snprintf( buffer, sizeof(buffer),
               "MeasurementInfo: Remark %i in LHS ('%s') doesnt match RHS ('%s')",
               int(i), nlhsremarkss[i].c_str(), nrhsremarkss[i].c_str() );
      throw runtime_error( buffer );
    }
  }

  vector<pair<string,string> > lhscompvsn, rhscompvsn;
  
  
  for( size_t i = 0; i < lhs.component_versions_.size(); ++i )
  {
    const string &n = lhs.component_versions_[i].first;
    if( n != "InterSpec" && n!= "InterSpecN42Serialization" && n != "Software"
        && !UtilityFunctions::istarts_with(n, "Original Software") )
      lhscompvsn.push_back( lhs.component_versions_[i] );
  }
  
  for( size_t i = 0; i < rhs.component_versions_.size(); ++i )
  {
    const string &n = rhs.component_versions_[i].first;
    if( n != "InterSpec" && n!= "InterSpecN42Serialization" && n != "Software"
        && !UtilityFunctions::istarts_with(n, "Original Software") )
      rhscompvsn.push_back( rhs.component_versions_[i] );
  }
  
  
  if( lhscompvsn.size() != rhscompvsn.size() )
  {
    snprintf( buffer, sizeof(buffer),
             "MeasurementInfo: Number of component versions in LHS (%i) doesnt match RHS (%i)",
             int(lhscompvsn.size()), int(rhscompvsn.size()) );
    
    for( size_t i = 0; i < lhscompvsn.size(); ++i )
      cout << "\tLHS: " << lhscompvsn[i].first << ": " << lhscompvsn[i].second << endl;
    for( size_t i = 0; i < rhscompvsn.size(); ++i )
      cout << "\tRHS: " << rhscompvsn[i].first << ": " << rhscompvsn[i].second << endl;
    throw runtime_error( buffer );
  }
  
  stable_sort( lhscompvsn.begin(), lhscompvsn.end() );
  stable_sort( rhscompvsn.begin(), rhscompvsn.end() );
  
  for( size_t i = 0; i < lhscompvsn.size(); ++i )
  {
    pair<string,string> lhsp = lhscompvsn[i];
    pair<string,string> rhsp = rhscompvsn[i];
    
    UtilityFunctions::trim( lhsp.first );
    UtilityFunctions::trim( lhsp.second );
    UtilityFunctions::trim( rhsp.first );
    UtilityFunctions::trim( rhsp.second );
    
    if( lhsp.first != rhsp.first )
    {
      snprintf( buffer, sizeof(buffer),
               "MeasurementInfo: Component Version %i name in LHS ('%s') doesnt match RHS ('%s')",
               int(i), lhsp.first.c_str(), rhsp.first.c_str() );
      throw runtime_error( buffer );
    }
    
    if( lhsp.second != rhsp.second )
    {
      snprintf( buffer, sizeof(buffer),
               "MeasurementInfo: Component Version %i valiue in LHS ('%s') doesnt match RHS ('%s')",
               int(i), lhsp.second.c_str(), rhsp.second.c_str() );
      throw runtime_error( buffer );
    }
  }//for( size_t i = 0; i < lhscompvsn.size(); ++i )
  
  
  //Make UUID last, since its sensitive to the other variables changing, and we
  //  want to report on those first to make fixing easier.
  if( !!lhs.detectors_analysis_ )
    DetectorAnalysis::equalEnough( *lhs.detectors_analysis_,
                                  *rhs.detectors_analysis_ );
  
  
  if( lhs.uuid_ != rhs.uuid_ )
    throw runtime_error( "MeasurementInfo: UUID of LHS (" + lhs.uuid_
                        + ") doesnt match RHS (" + rhs.uuid_ + ")" );

//  bool modified_;
//  bool modifiedSinceDecode_;
}//void equalEnough( const Measurement &lhs, const Measurement &rhs )
#endif //#if( PERFORM_DEVELOPER_CHECKS )


MeasurementInfo::MeasurementInfo()
{
  reset();
}


MeasurementInfo::~MeasurementInfo()
{
}


MeasurementInfo::MeasurementInfo( const MeasurementInfo &rhs )
{
  *this = rhs;
}


const MeasurementInfo &MeasurementInfo::operator=( const MeasurementInfo &rhs )
{
  if( this == &rhs )
    return *this;
//  std::unique_lock<std::recursive_mutex> lhs_lock( mutex_, std::defer_lock );
//  std::unique_lock<std::recursive_mutex> rhs_lock( rhs.mutex_, std::defer_lock );

  //XXX - this next section of code seems to really cause troubles!
  //      need to investigate it, and whatever code is calling it
  std::lock<std::recursive_mutex,std::recursive_mutex>( mutex_, rhs.mutex_ );
  std::unique_lock<std::recursive_mutex> lhs_lock( mutex_, std::adopt_lock_t() );
  std::unique_lock<std::recursive_mutex> rhs_lock( rhs.mutex_, std::adopt_lock_t() );

  reset();
  
  gamma_live_time_        = rhs.gamma_live_time_;
  gamma_real_time_        = rhs.gamma_real_time_;
  gamma_count_sum_        = rhs.gamma_count_sum_;
  neutron_counts_sum_     = rhs.neutron_counts_sum_;

  filename_               = rhs.filename_;
  detector_names_         = rhs.detector_names_;
  detector_numbers_       = rhs.detector_numbers_;
  neutron_detector_names_ = rhs.neutron_detector_names_;

  uuid_                   = rhs.uuid_;
  remarks_                = rhs.remarks_;
  lane_number_             = rhs.lane_number_;
  measurement_location_name_ = rhs.measurement_location_name_;
  inspection_             = rhs.inspection_;
  measurment_operator_    = rhs.measurment_operator_;
  sample_numbers_         = rhs.sample_numbers_;
  sample_to_measurments_  = rhs.sample_to_measurments_;
  detector_type_          = rhs.detector_type_;
  instrument_type_        = rhs.instrument_type_;
  manufacturer_           = rhs.manufacturer_;
  instrument_model_       = rhs.instrument_model_;
  instrument_id_          = rhs.instrument_id_;
  component_versions_     = rhs.component_versions_;
  mean_latitude_          = rhs.mean_latitude_;
  mean_longitude_         = rhs.mean_longitude_;
  properties_flags_       = rhs.properties_flags_;
  
  modified_               = rhs.modified_;
  modifiedSinceDecode_    = rhs.modifiedSinceDecode_;
  
  measurements_.clear();
  for( size_t i = 0; i < rhs.measurements_.size(); ++i )
  {
    MeasurementShrdPtr ptr( new Measurement( *rhs.measurements_[i] ) );
    measurements_.push_back( ptr );
  }

/*
  //As long as we've done everything right above, we shouldnt need to call
  //  either of the following, right?
  cleanup_after_load();
*/

  return *this;
}//operator=





const Measurement &Measurement::operator=( const Measurement &rhs )
{
  if( this == &rhs )
    return *this;

  live_time_ = rhs.live_time_;
  real_time_ = rhs.real_time_;

  contained_neutron_ = rhs.contained_neutron_;

  sample_number_ = rhs.sample_number_;
  occupied_ = rhs.occupied_;
  gamma_count_sum_ = rhs.gamma_count_sum_;
  neutron_counts_sum_ = rhs.neutron_counts_sum_;
  speed_ = rhs.speed_;
  detector_name_ = rhs.detector_name_;
  detector_number_ = rhs.detector_number_;
  detector_type_ = rhs.detector_type_;
  quality_status_ = rhs.quality_status_;

  source_type_ = rhs.source_type_;
  energy_calibration_model_ = rhs.energy_calibration_model_;

  remarks_ = rhs.remarks_;
  start_time_ = rhs.start_time_;

  calibration_coeffs_ = rhs.calibration_coeffs_;
  deviation_pairs_ = rhs.deviation_pairs_;

  channel_energies_ = rhs.channel_energies_;
  gamma_counts_ = rhs.gamma_counts_;

//  if( rhs.gamma_counts_ )
//    gamma_counts_.reset( new vector<float>( *rhs.gamma_counts_ ) );
//  else
//    gamma_counts_.reset();

  neutron_counts_ = rhs.neutron_counts_;
  
  latitude_       = rhs.latitude_;
  longitude_      = rhs.longitude_;
  position_time_  = rhs.position_time_;

  title_ = rhs.title_;

  return *this;
}//const Measurement &operator=( const Measurement &rhs )





bool MeasurementInfo::load_file( const std::string &filename,
               ParserType parser_type,
               std::string orig_file_ending )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  bool success = false;
  switch( parser_type )
  {
    case k2006Icd1Parser:
    case K2012ICD1Parser:
      success = load_N42_file( filename );
    break;

    case kSpcParser:
      success = load_spc_file( filename );
    break;

    case kGR135Parser:
      success = load_binary_exploranium_file( filename );
    break;

    case kPcfParser:
      success = load_pcf_file( filename );
    break;

    case kChnParser:
      success = load_chn_file( filename );
    break;

    case kIaeaParser:
      success = load_iaea_file( filename );
    break;

    case kTxtOrCsvParser:
      success = load_txt_or_csv_file( filename );
    break;

    case kCanberraCnfParser:
      success = load_camberra_cnf_file( filename );
    break;
      
    case kTracsMpsParser:
      success = load_tracs_mps_file( filename );
    break;
      
    case kAramParser:
      success = load_aram_file( filename );
    break;
      
    case kSPMDailyFile:
      success = load_spectroscopic_daily_file( filename );
    break;
      
    case kAmptekMca:
      success = load_amptek_file( filename );
    break;
      
    case kOrtecListMode:
      success = load_ortec_listmode_file( filename );
    break;
      
    case kMicroRaider:
      success = load_micro_raider_file( filename );
    break;
      
    case kAutoParser:
    {
      bool triedPcf = false, triedSpc = false,
          triedNativeIcd1 = false, triedTxt = false, triedGR135 = false,
          triedChn = false, triedIaea = false, triedCnf = false,
          triedMps = false, triedSPM = false, triedMCA = false,
          triedOrtecLM = false, triedMicroRaider = false, triedAram = false;
      if( !orig_file_ending.empty() )
      {
        const size_t period_pos = orig_file_ending.find_last_of( '.' );
        if( period_pos != string::npos )
          orig_file_ending = orig_file_ending.substr( period_pos+1 );
        to_lower( orig_file_ending );

        if( orig_file_ending=="pcf")
        {
          triedPcf = true;
          success = load_pcf_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="pcf")

        if( orig_file_ending=="dat")
        {
          triedGR135 = true;
          success = load_binary_exploranium_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="dat")

        if( orig_file_ending=="spc")
        {
          triedSpc = true;
          success = load_spc_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="dat")

        if( orig_file_ending=="n42" || orig_file_ending=="xml"
            || orig_file_ending=="icd1" || orig_file_ending=="icd")
        {
          triedNativeIcd1 = true;
          success = load_N42_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="n42")
        
        if( orig_file_ending=="chn" )
        {
          triedChn = true;
          success = load_chn_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="chn" )

        if( orig_file_ending=="spe" )
        {
          triedIaea = true;
          success = load_iaea_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="chn" )
        
        if( orig_file_ending=="txt" )
        {
          triedSPM = true;
          success = load_spectroscopic_daily_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" )
        
        if( orig_file_ending=="txt" || orig_file_ending=="csv" )
        {
          triedTxt = true;
          success = load_txt_or_csv_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )

        if( orig_file_ending=="cnf" )
        {
          triedCnf = true;
          success = load_camberra_cnf_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )
        
        if( orig_file_ending=="mps" )
        {
          triedMps = true;
          success = load_tracs_mps_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )
        
        if( orig_file_ending=="gam" /* || orig_file_ending=="neu" */ )
        {
          triedAram = true;
          success = load_aram_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )
        
        if( orig_file_ending=="mca" )
        {
          triedMCA = true;
          success = load_amptek_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )

        if( orig_file_ending=="lis" )
        {
          triedOrtecLM = true;
          success = load_ortec_listmode_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="txt" || orig_file_ending=="csv" )
        
        
        if( orig_file_ending=="xml" )
        {
          triedMicroRaider = true;
          success = load_micro_raider_file( filename );
          if( success ) break;
        }//if( orig_file_ending=="xml" )
      }//if( !orig_file_ending.empty() ) / else

      if( !success && !triedSpc )
        success = load_spc_file( filename );

      if( !success && !triedGR135 )
        success = load_binary_exploranium_file( filename );

      if( !success && !triedNativeIcd1 )
        success = load_N42_file( filename );

      if( !success && !triedPcf )
        success = load_pcf_file( filename );

      if( !success && !triedChn )
        success = load_chn_file( filename );

      if( !success && !triedIaea )
        success = load_iaea_file( filename );

      if( !success && !triedSPM )
        success = load_spectroscopic_daily_file( filename );
      
      if( !success && !triedTxt )
        success = load_txt_or_csv_file( filename );

      if( !success && !triedCnf )
        success = load_camberra_cnf_file( filename );
      
      if( !success && !triedMps )
        success = load_tracs_mps_file( filename );
      
      if( !success && !triedMCA )
        success = load_amptek_file( filename );
      
      if( !success && !triedMicroRaider )
        success = load_micro_raider_file( filename );
      
      if( !success && !triedAram )
        success = load_aram_file( filename );
      
      if( !success && !triedOrtecLM )
        success = load_ortec_listmode_file( filename );
      
       break;
    }//case kAutoParser
  };//switch( parser_type )

  set_filename( filename );

  if( num_measurements() == 0 )
    reset();

  return (success && num_measurements());
}//bool load_file(...)



bool comp_by_start_time_source( const MeasurementShrdPtr &lhs,
                           const MeasurementShrdPtr &rhs )
{
  if( !lhs || !rhs)
    return (!rhs) < (!lhs);
  
  const boost::posix_time::ptime &left = lhs->start_time();
  const boost::posix_time::ptime &right = rhs->start_time();
  
  if( left == right )
    return (lhs->source_type() < rhs->source_type());
  
  if( left.is_special() && !right.is_special() )
    return true;
  
  return (left < right);
}

void  MeasurementInfo::set_sample_numbers_by_time_stamp()
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  if( measurements_.empty() )
    return;

  //If we're here, we need to (re) assign some sample numbers
  
  //This function can be really slow, so I'm experimenting with a faster way
  //  For a certain large file in debug mode, took from 9.058687s to 0.232338s
  //  (in release mode this was from 1.808691s to 0.053467s)
  //  Right now this faster method is only enabled for really large files with
  //  more than 500 measurments - this is because the faster method does not
  //  preserve existing sample numbers

  if( measurements_.size() > 500 )
  {
    vector<MeasurementShrdPtr> sorted_meas;

    sorted_meas.reserve( measurements_.size() );

    for( auto &m : measurements_ )
    {
      if( m )
        sorted_meas.push_back( m );
    }//for( auto &m : measurements_ )

    stable_sort( sorted_meas.begin(), sorted_meas.end(), &comp_by_start_time_source );

    int sample_num = 1;
    vector<MeasurementShrdPtr>::iterator start, end, iter;
    for( start=end=sorted_meas.begin(); start != sorted_meas.end(); start=end )
    {
      while( (end != sorted_meas.end())
             && ((*end)->start_time_ == (*start)->start_time_)  )
        ++end;

      typedef map<string,int> StrIntMap;
      StrIntMap detectors;

      for( iter = start; iter != end; ++iter )
      {
        MeasurementShrdPtr &meas = *iter;

        if( detectors.count(meas->detector_name_) == 0 )
          detectors[meas->detector_name_] = 0;
        else
          ++detectors[meas->detector_name_];

        meas->sample_number_ = sample_num + detectors[meas->detector_name_];
      }//for( iter = start; iter != end; ++iter )

      int largest_delta = 0;
      for( const StrIntMap::value_type &val : detectors )
        largest_delta = max( largest_delta, val.second );

      sample_num = sample_num + largest_delta + 1;
    }//for( loop over time ranges )

  }else
  {
    typedef std::map<int, vector<MeasurementShrdPtr > > SampleToMeasMap;
    typedef map<boost::posix_time::ptime, SampleToMeasMap > TimeToSamplesMeasMap;
    
    TimeToSamplesMeasMap time_meas_map;

    for( const auto &m : measurements_ )
    {
      if( !m )
        continue;
      
      const int detnum = m->detector_number_;
      
      //If the time is invalid, we'll put this measurment after all the others.
      //If its an IntrinsicActivity, we'll put it before any of the others.
      if( m->source_type() == Measurement::IntrinsicActivity )
        time_meas_map[boost::posix_time::neg_infin][detnum].push_back( m );
      else if( m->start_time_.is_special() )
        time_meas_map[boost::posix_time::pos_infin][detnum].push_back( m );
      else
        time_meas_map[m->start_time_][detnum].push_back( m );
    }//for( auto &m : measurements_ )
    
    int sample = 1;
    
    for( TimeToSamplesMeasMap::value_type &t : time_meas_map )
    {
      SampleToMeasMap &measmap = t.second;
      
      //Make an attempt to sort the measurments in a reproducable, unique way
      //  because measurments wont be in the same order due to the decoding
      //  being multithreaded for some of the parsers.
//20150609: I think the multithreaded parsing has been fixed to yeild a
//  deterministic ordering always.  This only really matters for spectra where
//  the same start time is given for all records, or no start time at all is
//  given.  I think in these cases we want to assume the order that the
//  measurments were taken is the same as the order in the file.

      
      size_t nsamples = 0;
      for( const SampleToMeasMap::value_type &s : measmap )
        nsamples = std::max( nsamples, s.second.size() );
      
      for( size_t i = 0; i < nsamples; ++i )
      {
        for( SampleToMeasMap::value_type &s : measmap )
        {
          if( i < s.second.size() )
            s.second[i]->sample_number_ = sample;
        }
        
        ++sample;
      }//for( size_t i = 0; i < nsamples; ++i )
    }//for( const TimeToSamplesMeasMap::value_type &t : time_meas_map )
    
  }//if( measurements_.size() > 500 ) / else
  
  stable_sort( measurements_.begin(), measurements_.end(), &Measurement::compare_by_sample_det_time );
}//void  set_sample_numbers_by_time_stamp()


bool MeasurementInfo::has_unique_sample_and_detector_numbers() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  //First we will check that the comnination of sample numbers and detector
  //  numbers is unique, and also that all measurments with the same sample
  //  number have the same time stamp.  If we pass both of thers conditions,
  //  then we'll return since there is no need re-assign sample numbers.
  std::map<int, vector<int> > sampleNumsToSamples;
  std::map<int,std::set<boost::posix_time::ptime> > sampleToTimes;
  const size_t nmeas = measurements_.size();
  for( size_t i = 0; i < nmeas; ++i )
  {
    const MeasurementShrdPtr &m = measurements_[i];
    vector<int> &meass = sampleNumsToSamples[m->sample_number_];
    
    if( std::find(meass.begin(),meass.end(),m->detector_number_) != meass.end() )
      return false;
    
    meass.push_back( m->detector_number_ );
    
    std::set<boost::posix_time::ptime> &timesSet = sampleToTimes[m->sample_number_];
    if( !m->start_time_.is_special() )
      timesSet.insert( m->start_time_ );
    if( timesSet.size() > 1 )
      return false;
  }//for( MeasurementShrdPtr &m : measurements_ )
  
  return true;
}//bool has_unique_sample_and_detector_numbers() const


void  MeasurementInfo::ensure_unique_sample_numbers()
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( has_unique_sample_and_detector_numbers() )
  {
    stable_sort( measurements_.begin(), measurements_.end(), &Measurement::compare_by_sample_det_time );
  }else
  {
    set_sample_numbers_by_time_stamp();
  }
  
  
  //XXX - TODO should validate this a little further and check performance
  //      impact!  Also, should be proactive about not needing to hit the
  //      expensive "fix" bellow.
  //
  //Here we will check the first two sample number, and if they are '1' and '2'
  //  repectively, we will not do anything.  If first sample is not 1, but
  //  second sample is 2, we will change first sample to 1.  Otherwise if
  //  first two sample numbers isnt {1,2}, we will change all sample numbers to
  //  start at 1 and increase continuosly by 1.
  //  (note: this mess of logic is "inspired" by heuristics, and that actually
  //   looping through all measurements of a large file is expensive)
  set<int> sample_numbers;
  for( size_t i = 0; (sample_numbers.size() < 3) && (i < measurements_.size()); ++i )
    sample_numbers.insert( measurements_[i]->sample_number_ );
  
  if( sample_numbers.empty() )
    return;
  
  if( sample_numbers.size() == 1 )
  {
    for( auto &m : measurements_ )
      if( m->sample_number_ <= 0 )
        m->sample_number_ = 1;
    return;
  }
  
  const auto first_val = sample_numbers.begin();
  const auto second_val = next(first_val);
  int first_sample_val = *first_val;
  const int second_sample_val = *second_val;
  if( (first_sample_val + 1) != second_sample_val )
    first_sample_val = second_sample_val - 1;  //garunteed first_sample_val will be unique
  
  if( second_sample_val != 2 )
  {
    for( auto &m : measurements_ )
      sample_numbers.insert( m->sample_number_ );
    
    const vector<int> sample_numbers_vec( sample_numbers.begin(), sample_numbers.end() );
    const auto sn_begin = begin(sample_numbers_vec);
    const auto sn_end = end(sample_numbers_vec);
    
    for( auto &m : measurements_ )
    {
      auto pos = lower_bound( sn_begin, sn_end, m->sample_number_ );
      m->sample_number_ = static_cast<int>( (pos - sn_begin) + 1 );
    }
    
    return;
  }//if( (*first_val) != 1 )
  
  if( first_sample_val != (*first_val) )
  {
    const int old_first_sample = *first_val;
    for( size_t i = 0; (measurements_[i]->sample_number_!=second_sample_val) && (i < measurements_.size()); ++i )
      if( measurements_[i]->sample_number_ == old_first_sample )
        measurements_[i]->sample_number_ = first_sample_val;
  }
}//void ensure_unique_sample_numbers()


std::set<std::string> MeasurementInfo::find_detector_names() const
{
//  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  set<string> det_names;

  for( const auto &meas : measurements_ )
    det_names.insert( meas->detector_name_ );

  return det_names;
}//set<string> find_detector_names() const


void MeasurementInfo::cleanup_after_load( const unsigned int flags )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const bool rebinToCommonBinning = (flags & RebinToCommonBinning);
  
  
  //When loading the example passthrough N42 file, this function
  //  take about 60% of the parse time - due almost entirely to
  //  MeasurementInfo::rebin_by_polunomial_eqn
  try
  {
    recalc_total_counts();
    
    set<string> gamma_detector_names; //can be gamma+nutron
    const set<string> det_names = find_detector_names();
    
    int nGpsCoords = 0;
    mean_latitude_ = mean_longitude_ = 0.0;
    
    properties_flags_ = 0x0;
    
    vector<string> names( det_names.begin(), det_names.end() );
    typedef map<int,string> IntStrMap;
    IntStrMap num_to_name_map;
    set<string> neut_det_names;
    vector<string>::const_iterator namepos;
    
    for( size_t meas_index = 0; meas_index < measurements_.size(); ++meas_index )
    {
      MeasurementShrdPtr &meas = measurements_[meas_index];
      
      namepos = find( names.begin(), names.end(), meas->detector_name_ );
      if( namepos != names.end() )
      {
        meas->detector_number_ = static_cast<int>( namepos - names.begin() );
        num_to_name_map[meas->detector_number_] = meas->detector_name_;
      }else cerr << "Couldnt find detector '" << meas->detector_name_ << "'!" << endl;
      
      if( meas->gamma_counts_ && meas->gamma_counts_->size() )
        gamma_detector_names.insert( meas->detector_name_ );
      
#if( PERFORM_DEVELOPER_CHECKS )
      if( meas->neutron_counts_sum_ > 0.00001 && !meas->contained_neutron_ )
      {
        char buffer[1024];
        snprintf( buffer, sizeof(buffer),
                 "Spectrum contained %f neutrons, but neutron_counts_sum_ was not set. File=\"%s\"",
                 meas->neutron_counts_sum_, filename_.c_str() );
        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
      }
#endif
      
      if( meas->contained_neutron_ )
        neut_det_names.insert( meas->detector_name_ );
      
      //make sure the file didnt just have all zeros for the linear and
      //  higher calibration coefficients
      bool allZeros = true, invalid_calib = false;
      for( size_t i = 1; i < meas->calibration_coeffs_.size(); ++i )
        allZeros &= (fabs(meas->calibration_coeffs_[i])<1.0E-14);
      
      //Do a basic sanity check of is the calibration is reasonable.
      if( !allZeros && meas->gamma_counts_ && meas->gamma_counts_->size() )
      {
        switch( meas->energy_calibration_model_ )
        {
          case Measurement::Polynomial:
            if( meas->calibration_coeffs_.size() < 2
                || fabs(meas->calibration_coeffs_[0])>1000.0 )
            {
              allZeros = true;
              invalid_calib = true;
            }
          case Measurement::FullRangeFraction:
            if( meas->calibration_coeffs_.size() < 2
               || fabs(meas->calibration_coeffs_[0]) > 1000.0f
               || meas->calibration_coeffs_[1] < (FLT_EPSILON*meas->gamma_counts_->size()) )
            {
              allZeros = true;
              invalid_calib = true;
            }
            break;
            
          case Measurement::LowerChannelEdge:
            //should check that at least as many energies where provided as
            //  meas->gamma_counts_->size()
            break;
            
          case Measurement::UnknownEquationType:
            //Nothing to check here.
            break;
        }//switch( meas->energy_calibration_model_ )
      }//if( !allZeros && meas->gamma_counts_ && meas->gamma_counts_->size() )
      
      if( allZeros && (meas->energy_calibration_model_==Measurement::Polynomial
           || meas->energy_calibration_model_==Measurement::FullRangeFraction)
          && meas->gamma_counts_ && meas->gamma_counts_->size() )
      {
        if( invalid_calib )
        {
          stringstream remark;
          remark << "Energy calibration provided by detector was invalid {";
          for( size_t i = 0; i < meas->calibration_coeffs_.size(); ++i )
            remark << (i?", ":"") << meas->calibration_coeffs_[i];
          remark << "}";
          meas->remarks_.push_back( remark.str() );
        }//if( invalid_calib )
        
        meas->energy_calibration_model_ = Measurement::UnknownEquationType;
        meas->calibration_coeffs_.clear();
        meas->deviation_pairs_.clear();
        meas->channel_energies_.reset();
      }//if( we dont have a calibration...)
      
      if( allZeros && meas->gamma_counts_ && meas->gamma_counts_->size() )
      {
        //look to see if we can grab a calibration for this same detector from
        // another sample.
        
        //MeasurementShrdPtr &meas = measurements_[meas_index];
        
        //invalid_calib
      }//if( allZeros )
      
      if( meas->has_gps_info() )
      {
        //        a lat long of (0 0) probably isnt a valid GPS coordinate
        if( (fabs(meas->latitude_) < 1.0E-6)
           && (fabs(meas->longitude_) < 1.0E-6) )
        {
          meas->latitude_ = meas->longitude_ = -999.9;
          meas->position_time_ = boost::posix_time::not_a_date_time;
        }else
        {
          ++nGpsCoords;
          mean_latitude_ += meas->latitude();
          mean_longitude_ += meas->longitude();
        }
      }else if( !meas->position_time_.is_special() )
      {
        meas->position_time_ = boost::posix_time::not_a_date_time;
      }
      
      meas->contained_neutron_ |= (meas->neutron_counts_sum_>0.0
                                   || !meas->neutron_counts_.empty());
    }//for( auto &meas : measurements_ )
    
    
    if( nGpsCoords == 0
       || (fabs(mean_latitude_)<1.0E-6 && fabs(mean_longitude_)<1.0E-6) )
    {
      mean_latitude_ = mean_longitude_ = -999.9;
    }else
    {
      mean_latitude_ /= nGpsCoords;
      mean_longitude_ /= nGpsCoords;
    }//if( !nGpsCoords ) / else
    
    if( !Measurement::valid_longitude(mean_longitude_)
       || !Measurement::valid_latitude(mean_latitude_) )
      mean_latitude_ = mean_longitude_ = -999.9;
    
    if( flags & DontChangeOrReorderSamples )
    {
      if( !has_unique_sample_and_detector_numbers() )
        properties_flags_ |= kNotUniqueSampleDetectorNumbers;
      
      for( size_t i = 1; i < measurements_.size(); ++i )
      {
        if( measurements_[i-1]->start_time_.is_special()
           || measurements_[i]->start_time_.is_special() )
          continue;
        
        if( measurements_[i-1]->start_time_ > measurements_[i]->start_time_ )
          properties_flags_ |= kNotTimeSortedOrder;
        
        if( !Measurement::compare_by_sample_det_time(measurements_[i-1],measurements_[i]) )
          properties_flags_ |= kNotSampleDetectorTimeSorted;
      }//for( size_t i = 1; i < measurements_.size(); ++i )
    }else
    {
      ensure_unique_sample_numbers();
    }
    
    detector_numbers_.clear();
    detector_names_.clear();
    neutron_detector_names_.clear();
    
    for( const IntStrMap::value_type &t : num_to_name_map )
    {
      detector_numbers_.push_back( t.first );
      detector_names_.push_back( t.second );
    }//for( const IntStrMap::value_type &t : num_to_name_map )
    
    neutron_detector_names_.insert( neutron_detector_names_.end(), neut_det_names.begin(), neut_det_names.end() );
    
    //Make a set of all unique calibrations, so we can later allow Measurements
    //  with the same calibration share a vector<float> of bin energies.
    set<MeasurementCalibInfo> calib_infos_set;
    
    //Build up a mapping from detector names, to _a_ calibration for it, so if
    //  any spectra for that detector doesnt have a calibration associated with
    //  it, we can assign it this one.  Note: the MeasurementCalibInfo wont be
    //  used to assign calibration from, it will be used to index
    //  calib_infos_set.
    map<string,MeasurementCalibInfo> detname_to_calib_map;
    
    //Need to go through and set binning on spectrums here
    for( auto &meas : measurements_ )
    {
      if( !meas->gamma_counts_ || meas->gamma_counts_->empty() )
        continue;
      
      MeasurementCalibInfo info( meas );
      
      if( info.coefficients.size()
         || (info.binning && info.binning->size()) )
      {
        calib_infos_set.insert( info );
        
        
        const string &name = meas->detector_name_;
        //        if( detname_to_calib_map.find(name) != detname_to_calib_map.end()
        //            && !(detname_to_calib_map[name] == info) )
        //          cerr << SRC_LOCATION << "\n\tWarning same detector " << name
        //               << " has multiple calibrations" << endl;
        detname_to_calib_map[name] = info;
      }//if( a non-empty info ) / else
    }//for( auto &meas : measurements_ )
    
    
    int n_times_guess_cal = 0;
    size_t nbins = 0;
    bool all_same_num_bins = true;
    
    
    //In order to figure out if a measurment is passthrough, we'll use these
    //  two variables.
    int pt_num_items = 0;
    float pt_averageRealTime = 0;
    
    //ensure_unique_sample_numbers(); has already been called a
    sample_numbers_.clear();
    sample_to_measurments_.clear();
    
    typedef std::pair<boost::posix_time::ptime,float> StartAndRealTime;
    typedef map<int, StartAndRealTime > SampleToTimeInfoMap;
    SampleToTimeInfoMap samplenum_to_starttime;
    
    const size_t nmeas = measurements_.size();
    for( size_t measn = 0; measn < nmeas; ++measn )
    {
      MeasurementShrdPtr &meas = measurements_[measn];
      sample_numbers_.insert( meas->sample_number_ );
      sample_to_measurments_[meas->sample_number_].push_back( measn );
      
      if( !meas->gamma_counts_ || meas->gamma_counts_->empty() )
        continue;
      
      if( nbins==0 )
        nbins = meas->gamma_counts_->size();
      if( nbins != meas->gamma_counts_->size() )
        all_same_num_bins = false;
      
      //20180221: Removed check on measurement type because some ROSA portal
      //  occupancies have less than 4 non-background samples.
      if( //meas->source_type_ != Measurement::Background
         //&& meas->source_type_ != Measurement::Calibration &&
        meas->source_type_ != Measurement::IntrinsicActivity
         && meas->sample_number() >= 0
         && meas->live_time() > 0.00000001
         && meas->real_time() > 0.00000001
         && meas->real_time() < 5.0 )
      {
        ++pt_num_items;
        pt_averageRealTime += meas->real_time_;
        
        if( !meas->start_time().is_special() )
        {
          const int samplenum = meas->sample_number();
          const boost::posix_time::ptime &st = meas->start_time();
          const float rt = meas->real_time();
          
          SampleToTimeInfoMap::iterator pos = samplenum_to_starttime.find( samplenum );
          if( pos == samplenum_to_starttime.end() )
            pos = samplenum_to_starttime.insert( make_pair(samplenum, make_pair(st,rt)) ).first;
          pos->second.second = max( rt, pos->second.second );
        }
        
        
      }//if( a candidate for a passthrough spectrum )
      
      MeasurementCalibInfo thisinfo( meas );
      set<MeasurementCalibInfo>::const_iterator pos;
      pos = calib_infos_set.find(thisinfo); //
      
      //If this Measurement doesnt have a calibration, see if another one for
      //  the same detector exists.
      if( pos == calib_infos_set.end() )
        pos = calib_infos_set.find(detname_to_calib_map[meas->detector_name_]);
      
      
      try
      {
        if( pos == calib_infos_set.end() )
          throw runtime_error( "" );
      
        pos->fill_binning();
        meas->channel_energies_ = pos->binning;
        meas->calibration_coeffs_ = pos->coefficients;
        meas->deviation_pairs_ = pos->deviation_pairs_;
        meas->energy_calibration_model_ = pos->equation_type;
        
        if( !meas->channel_energies_ )
          throw runtime_error( "" );
          
        //Force calibration_coeffs_ to free the memory for LowerChannelEdge
        if( pos->equation_type == Measurement::LowerChannelEdge )
          vector<float>().swap( meas->calibration_coeffs_ );
      }catch( std::exception & )
      {
        if( meas->gamma_counts_ && !meas->gamma_counts_->empty())
        {
          if( !n_times_guess_cal && manufacturer_!="ICx Technologies" )
          {
            const char * msg = "Warning, could not find a calibration, so "
            "assuming energy range 0 to 3MeV, with equal sized bins";
            passMessage( msg, "", 1 );
          }
        
          ++n_times_guess_cal;
          meas->energy_calibration_model_ = Measurement::FullRangeFraction;
          meas->calibration_coeffs_.resize( 2 );
          meas->calibration_coeffs_[0] = 0.0;
          meas->calibration_coeffs_[1] = 3000.0;
          MeasurementCalibInfo info( meas );
        
          pos = calib_infos_set.find(info);
        
          if( pos == calib_infos_set.end() )
          {
            calib_infos_set.insert( info );
            pos = calib_infos_set.find(info);
          }//if( pos == calib_infos_set.end() )
        
          if( pos != calib_infos_set.end() )
          {
            pos->fill_binning();
            meas->channel_energies_ = pos->binning;
          }else
          {
            //shouldnt ever get here! - but JIC
            info.fill_binning();
            meas->channel_energies_ = info.binning;
          }//if( pos != calib_infos_set.end() )  e;se
        }//if( meas->gamma_counts_ && meas->gamma_counts_->size() )
      }//try( to fill binning ) / catch
    }//for( auto &meas : measurements_ )
    
    
    bool is_passthrough = true;
    
    if( sample_numbers_.size() < 5 || detector_numbers_.empty() )
      is_passthrough = false;
    
    if( pt_averageRealTime <= 0.00000001 )
      is_passthrough = false;
    
    pt_averageRealTime /= pt_num_items;
    
    //In principle should check that measurements were taken sequentially as well
    //is_passthrough = is_passthrough && ( (pt_num_items>5) && (pt_averageRealTime < 2.5) );
    is_passthrough = is_passthrough && ( (pt_num_items>5) && pt_num_items > static_cast<size_t>(0.75*nmeas) );
    
    //Go through and verify the measurments are mostly sequential (one
    //  measurement right after the next), for at least the occupied/item/unknown
    //  portions of the file.  The background portion, as well as the transition
    //  between background/foreground/etc is a little harder to deal with since
    //  background may have gaps, or may have been taken a while before the item
    if( !is_passthrough && samplenum_to_starttime.size() > 20 )
    {
      int nnotadjacent = 0, nadjacent = 0;
      
      SampleToTimeInfoMap::const_iterator next = samplenum_to_starttime.begin(), iter;
      for( iter = next++; next != samplenum_to_starttime.end(); ++iter, ++next )
      {
        const boost::posix_time::ptime &st = iter->second.first;
        const boost::posix_time::ptime &next_st = next->second.first;
        
        const float rt = iter->second.second;
        const boost::posix_time::time_duration duration = boost::posix_time::microseconds( rt*1.0E6 );
        
        boost::posix_time::time_duration diff = ((st + duration) - next_st);
        if( diff.is_negative() )
          diff = -diff;
        
        if( diff < (duration/100) )
          ++nadjacent;
        else
          ++nnotadjacent;
      }
      
      is_passthrough = (10*nnotadjacent < nadjacent);
    }//if( !is_passthrough && pt_num_items>20 )
    
  
    if( all_same_num_bins )
      properties_flags_ |= kAllSpectraSameNumberChannels;
    
    if( is_passthrough )
      properties_flags_ |= kPassthroughOrSearchMode;
    
    //XXX - as of 20121130, the bellow code which makes sure all samples
    //      of passthrough data has bining information (as apposed to just first
    //      detector that some of the detector manufactureers do), has not been
    //      well tested - just kinda seems to work for SAIC detectors.
    if( is_passthrough )
    {
      const int first_sample_number = *sample_numbers_.begin();
      map<string,MeasurementShrdPtr> det_meas_map;
      
      for( const auto &meas : measurements_ )
      {
        if( meas->sample_number() == first_sample_number
           && meas->gamma_counts_ && !meas->gamma_counts_->empty()
                                     && (!meas->calibration_coeffs_.empty()
                                         || (meas->channel_energies_ && !meas->channel_energies_->empty())) )
          det_meas_map[meas->detector_name_] = meas;
      }//for( const auto &meas : measurements_ )
      
      for( const auto &meas : measurements_ )
      {
        if( meas->sample_number() == first_sample_number )
          continue;
        
        const string &name = meas->detector_name_;
        if( meas->gamma_counts_ && !meas->gamma_counts_->empty()
                                   && det_meas_map.count( name )
           && meas->calibration_coeffs_.empty()
              && ( !meas->channel_energies_ || meas->channel_energies_->empty()
                                               || (det_meas_map[name]->energy_calibration_model_ == Measurement::LowerChannelEdge
                   && meas->energy_calibration_model_==Measurement::UnknownEquationType) ) )
        {
          if( !meas->channel_energies_ || meas->channel_energies_->empty())
          {
            meas->channel_energies_   = det_meas_map[name]->channel_energies_;
            meas->deviation_pairs_    = det_meas_map[name]->deviation_pairs_;
          }//if( !meas->channel_energies_ || !meas->channel_energies_->size() )
          
          meas->energy_calibration_model_  = det_meas_map[name]->energy_calibration_model_;
          meas->calibration_coeffs_ = det_meas_map[name]->calibration_coeffs_;
        }
      }//for( const MeasurementShrdPtr &meas : measurements_ )
    }//if( is_passthrough )
    
    if( rebinToCommonBinning && all_same_num_bins && !measurements_.empty()
       && ((gamma_detector_names.size() > 1) || is_passthrough) )
    {
      properties_flags_ |= kHasCommonBinning;
      
      if( (calib_infos_set.size() > 1) )
      {
        properties_flags_ |= kRebinnedToCommonBinning;
        
        //If we have more than one gamma detector, than we have to move them
        //  to have a common energy binning, to display properly
        size_t nbin = 0;
        float min_energy = 99999.9f, max_energy = -99999.9f;
        for( const auto &meas : measurements_ )
        {
          nbin = max( nbin , meas->channel_energies_->size() );
          if(!meas->channel_energies_->empty())
          {
            min_energy = min( min_energy, meas->channel_energies_->front() );
            max_energy = max( max_energy, meas->channel_energies_->back() );
          }//if( meas->channel_energies_.size() )
        }//for( const auto &meas : measurements_ )
        
        size_t nbinShift = nbin - 1;
        const float channel_width = (max_energy - min_energy) /
        static_cast<float>( nbinShift );
        vector<float> poly_eqn( 2, 0.0 );
        poly_eqn[0] = min_energy;  // - 0.5*channel_width; // min_energy;
        poly_eqn[1] = channel_width;
        
        for( const auto &meas : measurements_ )
        {
          if( meas->gamma_counts_->size() > 16 )//dont rebin SAIC or LUDLUMS
          {
            rebin_by_eqn( poly_eqn, DeviationPairVec(), Measurement::Polynomial );
            break;
          }//if( meas->gamma_counts_->size() > 16 )
        }//for( const auto &meas : measurements_ )
      }//if( calib_infos_set.size() > 1 )
    }else if( all_same_num_bins && !measurements_.empty() && calib_infos_set.size()==1 )
    {
      properties_flags_ |= kHasCommonBinning;
    }else if( !all_same_num_bins )
    {
    }//if( !measurements_.empty() )
    
    if( uuid_.empty() )
      uuid_ = generate_psuedo_uuid();
    
    set_detector_type_from_other_info();

    //Lets get rid of duplicate component_versions
    set< pair<string,string> > component_versions_s;
    std::vector<std::pair<std::string,std::string> > component_versions_nondup;
    for( size_t i = 0; i < component_versions_.size(); ++i )
    {
      if( component_versions_s.count( component_versions_[i] ) )
        continue;
      component_versions_s.insert( component_versions_[i] );
      component_versions_nondup.push_back( component_versions_[i] );
    }
    component_versions_.swap( component_versions_nondup );
    
    
#if( PERFORM_DEVELOPER_CHECKS )
    //Check to make sure all neutron detector names can be found in detector names
    {
      const vector<string>::const_iterator begindet = detector_names_.begin();
      const vector<string>::const_iterator enddet = detector_names_.end();
      for( const std::string &ndet : neutron_detector_names_ )
      {
        if( std::find(begindet,enddet,ndet) == enddet )
        {
          char buffer[1024];
          snprintf( buffer, sizeof(buffer),
                   "Found a neutron detector name not in the list of all detector names: %s\n",
                   ndet.c_str() );
          log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
        }
      }
    }
#endif
    
#if( PERFORM_DEVELOPER_CHECKS )
    static int ntest = 0;
    if( ntest++ < 10 )
      cerr << "Warning, testing rebingin ish" << endl;
    
    const double prev_gamma_count_sum_ = gamma_count_sum_;
    const double prev_neutron_counts_sum_ = neutron_counts_sum_;
    
    recalc_total_counts();
    
    if( fabs(gamma_count_sum_ - prev_gamma_count_sum_) > 0.01 )
    {
      char buffer[1024];
      snprintf( buffer, sizeof(buffer),
               "Before rebinning and gamma coutn sum=%10f and afterwards its %10f\n",
               prev_gamma_count_sum_, gamma_count_sum_ );
      log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
    }
    
    if( fabs(neutron_counts_sum_ - prev_neutron_counts_sum_) > 0.01 )
    {
      char buffer[1024];
      snprintf( buffer, sizeof(buffer),
               "Before rebinning and gamma count sum=%10f and afterwards its %10f\n",
               prev_neutron_counts_sum_, neutron_counts_sum_ );
      log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
    }
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  }catch( std::exception &e )
  {
    stringstream msg;
    msg << "From " << SRC_LOCATION << " caught error:\n\t" << e.what();
    throw runtime_error( msg.str() );
  }//try / catch
  
  modified_ = modifiedSinceDecode_ = false;
}//void cleanup_after_load()


void MeasurementInfo::set_detector_type_from_other_info()
{
  using UtilityFunctions::contains;
  using UtilityFunctions::icontains;
  
  
  if( detector_type_ != kUnknownDetector )
    return;
  
  const string &model = instrument_model_;
//  const string &id = instrument_id_;
  
  if( icontains(model,"SAM")
      && (contains(model,"940") || icontains(model,"Eagle+")) )
  {
    if( icontains(model,"LaBr") )
      detector_type_ = kSam940LaBr3;
    else
      detector_type_ = kSam940;
    
    cerr << "ASAm940 model=" << model << endl;
    
    return;
  }
  
  if( icontains(model,"SAM") && contains(model,"945") )
  {
    //if( icontains(model,"LaBr") )
      //detector_type_ = kSam945LaBr3;
    //else
    detector_type_ = kSam945;
    return;
  }
  
  //Dont know what the 'ULCS' that some models have in their name is
  if( icontains(model,"identiFINDER") && icontains(model,"NG") )
  {
    detector_type_ = kIdentiFinderNGDetector;
    return;
  }
  
  if( icontains(model,"identiFINDER") && icontains(model,"LG") )
  {
    detector_type_ = kIdentiFinderLaBr3Detector;
    return;
  }
  
  if( icontains(model,"RS-701") )
  {
    detector_type_ = kRsi701;
    return;
  }
  
  if( icontains(model,"RS-705") )
  {
    detector_type_ = kRsi705;
    return;
  }
  
  if( icontains(model,"RS???") /*&& icontains(id,"Avid")*/ )
  {
    detector_type_ = kAvidRsi;
    return;
  }
  
  if( icontains(model,"radHUNTER") )
  {
    if( icontains(model,"UL-LGH") )
      detector_type_ = kRadHunterLaBr3;
    else
      detector_type_ = kRadHunterNaI;
    return;
  }
  
  
  if( icontains(model,"RadEagle") || icontains(model,"RE 3") || icontains(model,"RE 2") )
  {
    if( UtilityFunctions::icontains(model,"3SG") ) //RADEAGLE NaI(Tl) 3x1, GM Handheld RIID
    {
      detector_type_ = kOrtecRadEagleNai;
    }else if( UtilityFunctions::icontains(model,"2CG") ) //RADEAGLE CeBr3 2x1, GM Handheld RIID.
    {
      detector_type_ = kOrtecRadEagleCeBr2Inch;
    }else if( UtilityFunctions::icontains(model,"3CG") ) //RADEAGLE CeBr3 3x0.8, GM Handheld RIID
    {
      detector_type_ = kOrtecRadEagleCeBr3Inch;
    }else if( UtilityFunctions::icontains(model,"2LG") ) //RADEAGLE LaBr3(Ce) 2x1, GM Handheld RIID
    {
      detector_type_ = kOrtecRadEagleLaBr;
    }else
    {
#if(PERFORM_DEVELOPER_CHECKS)
      log_developer_error( BOOST_CURRENT_FUNCTION, ("Unrecognized RadEagle Model: " + model).c_str() );
#endif
    }
    
    //Make the modle human readable
    if( istarts_with(model, "RE ") )
      instrument_model_ = "RadEagle " + instrument_model_.substr(3);
    
    return;
  }//if( a rad eagle )
}//void set_detector_type_from_other_info()


#if( PERFORM_DEVELOPER_CHECKS )
double MeasurementInfo::deep_gamma_count_sum() const
{
  double deep_gamma_sum = 0.0;
  for( const auto &meas : measurements_ )
  {
    if( !meas )
      continue;
    if( !!meas->gamma_counts_ )
    {
      for( const float f : *(meas->gamma_counts_) )
        deep_gamma_sum += f;
    }
  }//for( const auto &meas : measurements_ )
  
  return deep_gamma_sum;
}//double deep_gamma_count_sum() const

double MeasurementInfo::deep_neutron_count_sum() const
{
  double deep_sum = 0.0;
  for( const auto &meas : measurements_ )
  {
    if( !meas )
      continue;
    for( const float f : meas->neutron_counts_ )
      deep_sum += f;
  }//for( const auto &meas : measurements_ )
  
  return deep_sum;
}//double deep_neutron_count_sum() const;

#endif //#if( PERFORM_DEVELOPER_CHECKS )


void MeasurementInfo::recalc_total_counts()
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  gamma_live_time_ = 0.0f;
  gamma_real_time_ = 0.0f;
  gamma_count_sum_ = 0.0;
  neutron_counts_sum_ = 0.0;
  
  //We have an issue that the same data can be repeated, due to multiple energy
  //  calibrations, doubling neutron or gamma counts...  Right now we are just
  //  ignoring this.
  
  for( const auto &meas : measurements_ )
  {
    if( meas )
    {
      if( meas->gamma_counts_ && !meas->gamma_counts_->empty())
      {
        gamma_live_time_ += meas->live_time_;
        gamma_real_time_ += meas->real_time_;
      }

      gamma_count_sum_    += meas->gamma_count_sum_;
      neutron_counts_sum_ += meas->neutron_counts_sum_;
    }//if( meas )
  }//for( const auto &meas : measurements_ )
  
#if( PERFORM_DEVELOPER_CHECKS )
  const double deep_gamma_sum = deep_gamma_count_sum();
  const double deep_neutron_sum = deep_neutron_count_sum();
  
  if( fabs(deep_gamma_sum - gamma_count_sum_) > 0.1
      && fabs(deep_gamma_sum - gamma_count_sum_) > 1.0E-7*max(deep_gamma_sum, gamma_count_sum_) )
  {
    char buffer[1024];
    snprintf( buffer, sizeof(buffer),
             "recalc_total_counts() found a discrepance for sum gammas depending"
            " on if a shallow or deep count was done: %9f for shallow, %9f for"
            " deep\n", gamma_count_sum_, deep_gamma_sum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
  
  if( fabs(deep_neutron_sum - neutron_counts_sum_) > 0.1 )
  {
    char buffer[1024];
    snprintf( buffer, sizeof(buffer),
              "recalc_total_counts() found a discrepance for sum nuetrons depending"
              " on if a shallow or deep count was done: %9f for shallow, %9f for"
              " deep\n", neutron_counts_sum_, deep_neutron_sum );
    log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
  }
  
#endif //#if( PERFORM_DEVELOPER_CHECKS )
}//void recalc_total_counts()




std::string MeasurementInfo::generate_psuedo_uuid() const
{
  std::size_t seed = 0;
  
  
  boost::hash_combine( seed, gamma_live_time_ );
  boost::hash_combine( seed, gamma_real_time_ );
  
  boost::hash_combine( seed, gamma_count_sum_ );
  boost::hash_combine( seed, neutron_counts_sum_ );
//  boost::hash_combine( seed, filename_ );
  boost::hash_combine( seed, detector_names_ );
//  boost::hash_combine( seed, detector_numbers_ );
  boost::hash_combine( seed, neutron_detector_names_ );
  if( !remarks_.empty() )
    boost::hash_combine( seed, remarks_ );
  boost::hash_combine( seed, lane_number_ );
  if( !measurement_location_name_.empty() )
  boost::hash_combine( seed, measurement_location_name_ );
  if( !inspection_.empty() )
    boost::hash_combine( seed, inspection_ );
  
  boost::hash_combine( seed, instrument_type_ );
  boost::hash_combine( seed, manufacturer_ );
  boost::hash_combine( seed, instrument_model_ );
  
  if( Measurement::valid_latitude(mean_latitude_)
     && Measurement::valid_longitude(mean_longitude_) )
  {
    boost::hash_combine( seed, mean_latitude_ );
    boost::hash_combine( seed, mean_longitude_ );
  }
  
  //Note, not including properties_flags_
  
  boost::hash_combine( seed, instrument_id_ );
  boost::hash_combine( seed, measurements_.size() );
//  boost::hash_combine( seed, detectors_analysis_ );
  boost::hash_combine( seed, int(detector_type_) );
  boost::hash_combine( seed, measurment_operator_ );
  
  for( const std::shared_ptr<const Measurement> meas : measurements_ )
  {
    boost::hash_combine( seed, meas->live_time() );
    boost::hash_combine( seed, meas->real_time() );
    boost::hash_combine( seed, meas->gamma_count_sum() );
    boost::hash_combine( seed, meas->neutron_counts_sum() );
    
    if( Measurement::valid_latitude(meas->latitude_) )
      boost::hash_combine( seed, meas->latitude_ );
    if( Measurement::valid_longitude(meas->longitude_) )
      boost::hash_combine( seed, meas->longitude_ );
    //  boost::hash_combine( seed, position_time_ );
  }//for( const std::shared_ptr<const Measurement> meas : measurements_ )

  
//  std::set<int> sample_numbers_;
//  std::shared_ptr<const DetectorAnalysis> detectors_analysis_;
  
  string uuid;
  
  if(!measurements_.empty() && measurements_[0]
      && !measurements_[0]->start_time().is_special() )
    uuid = UtilityFunctions::to_iso_string( measurements_[0]->start_time() );
  else
    uuid = UtilityFunctions::to_iso_string( time_from_string( "1982-07-28 23:59:59" ) );
  
  //uuid something like: "20020131T100001,123456789"
  if( uuid.size() >= 15 )
    uuid = uuid.substr(2,6) + uuid.substr(9,2) + "-" + uuid.substr(11,4) + "-4"
           + ((uuid.size()>=18) ? uuid.substr(16,2) : string("00"));
  
  const uint64_t seed64 = seed;
//  char buffer[64];
//  snprintf( buffer, sizeof(buffer), "%.16llu", seed64 ); //PRIu64
//  const string seedstr = buffer;
  stringstream seedstrm;
  seedstrm  << setw(16) << setfill('0') << seed64;
  const string seedstr = seedstrm.str();
  
  if( seedstr.size() >= 16 )
    uuid += seedstr.substr(0,1) + "-a" + seedstr.substr(1,3) + "-"
            + seedstr.substr(4,12);
  
  //Now in form (using 1982 data above) of 82072823-5959-4xxx-a-xxxxxxxxxxxx
  //  and where the x are hexadecimal digits
  
  return uuid;
}//std::string generate_psuedo_uuid() const


bool MeasurementInfo::load_spc_file( const std::string &filename )
{
  reset();
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;

  uint8_t firstbyte;
  file.read( (char *) (&firstbyte), 1 );
  file.seekg( 0, ios::beg );
  const bool isbinary = (firstbyte == 0x1);
//  const bool istext = (char(firstbyte) == 'S');

  if( !isbinary && /*!istext*/ !isalpha(firstbyte) )
  {
//    cerr << "SPC file '" << filename << "'is not binary or text firstbyte="
//         << (char)firstbyte << endl;
    return false;
  }//if( !isbinary && !istext )

  bool loaded = false;
  if( isbinary )
    loaded = load_from_binary_spc( file );
  else
    loaded =  load_from_iaea_spc( file );

  if( loaded )
    filename_ = filename;

  return loaded;
}//bool load_spc_file( const std::string &filename )



bool MeasurementInfo::load_chn_file( const std::string &filename )
{
  reset();
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;

  uint8_t firstbyte;
  file.read( (char *) (&firstbyte), 1 );
  file.seekg( 0, ios::beg );

  if( firstbyte != uint8_t(255) )
  {
//    cerr << "CHN file '" << filename << "'does not have expected first byte of"
//         << " 255 firstbyte=" << int(firstbyte) << endl;
    return false;
  }//if( wrong first byte )

  const bool loaded = load_from_chn( file );

  if( loaded )
    filename_ = filename;

  return loaded;
}//bool load_chn_file( const std::string &filename )


bool MeasurementInfo::load_iaea_file( const std::string &filename )
{
  reset();
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;

  uint8_t firstbyte;
  file.read( (char *) (&firstbyte), 1 );
  file.seekg( 0, ios::beg );

  if( firstbyte != '$' )
  {
//    cerr << "IAEA file '" << filename << "'does not have expected first chacter"
//         << " of '$', firstbyte=" << int(firstbyte)
//         << " (" << char(firstbyte) << ")" << endl;
    return false;
  }//if( wrong first byte )

  const bool loaded = load_from_iaea( file );

  if( loaded )
    filename_ = filename;

  return loaded;
}//bool load_iaea_file(...)


bool MeasurementInfo::load_from_iaea_spc( std::istream &input )
{
  //Function is currently not very robust to line ending changes, or unexpected
  //  whitespaces.  Aslo parsing of channel counts coult be sped up probably.

  reset();
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  std::shared_ptr<DetectorAnalysis> analysis;
  
  MeasurementShrdPtr meas = std::make_shared<Measurement>();

  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();

//There are quite a number of fileds that Measurment or MeasurmentInfo class
//  does not yet implement, so for now we will just put them into the remarks
  
  string detector_type = "";

  try
  {
  string line;
  int Length = -1;


  //going through and making sure this is an ASCII file wont work, because
  //  there is often a subscript 3 (ascii code 179) or infinity symbols...
  //  So instead we'll insist the first non-empty line of file must start with
  //  three different alphanumneric characters.  We could probably tighten this
  //  up to apply to all non-empty lines in the file.
  int linenum = 0, nnotrecognized = 0;
  bool tested_first_line = false;

  while( input.good() )
  {
    const istream::pos_type sol_pos = input.tellg();
    
//    getline( input, line, '\r' );
    
    const size_t max_len = 1024*1024;  //allows a line length of 64k fields, each of 16 characters, which his more than any spectrum file should get
    UtilityFunctions::safe_get_line( input, line, max_len );
    
    if( line.size() >= (max_len-1) )
      throw runtime_error( "Line greater than 1MB" );
    
    trim( line );

    if( line.empty() )
      continue;

    if( !tested_first_line )
    {
      tested_first_line = true;
      if( line.size() < 3
          || !isalnum(line[0]) || !isalnum(line[1]) || !isalnum(line[2])
          || !(line[0]!=line[1] || line[0]!=line[2] || line[2]!=line[3]) )
        throw runtime_error( "File failed constraint that first three charcters"
                             " on the first non-empty line must be alphanumeric"
                             " and not be equal" );
    }//if( !tested_first_line )


    const size_t colonpos = line.find(':');
    const size_t info_pos = line.find_first_not_of(": ", colonpos);

    bool is_remark = false;
    //Check if its a remark field
    for( const char * const label : ns_iaea_comment_labels )
    {
      if( istarts_with( line, label ) )
      {
        is_remark = true;
        if( info_pos != string::npos )
        {
          string remark = label;
          remark += " : " + line.substr(info_pos);
          trim( remark );
          remarks_.push_back( remark );
        }
        break;
      }//if( istarts_with( line, label ) )
    }//for( const char * const label : ns_iaea_comment_labels )

    bool is_version = false;
    for( const char * const label : ns_iaea_version_labels )
    {
      if( istarts_with( line, label ) )
      {
        is_version = true;
        if( info_pos != string::npos )
          component_versions_.push_back( make_pair(label,line.substr(info_pos)) );
        break;
      }//if( istarts_with( line, label ) )
    }//for( const char * const label : ns_iaea_comment_labels )
    
    
    if( is_version )
    {
      //nothing to to here
    }else if( is_remark )
    {
      //Go through and look for warning signs...
      if( istarts_with( line, "BackgroundSubtraction" )
          && info_pos != string::npos
          && !icontains( line.substr(info_pos), "No" ) )
        passMessage( "Instrument may have been in background subtract mode.",
                     "load_from_iaea_spc(...)", 3 );
    }else if( istarts_with( line, "SpectrumName" ) )
    {//SpectrumName        : ident903558-21_2012-07-26_07-10-55-003.spc
      if( info_pos != string::npos )
      {
        if( UtilityFunctions::icontains( line.substr(info_pos), "ident") )
        {
          detector_type_ = kIdentiFinderNGDetector;
                           //TODO: kIdentiFinderLaBr3Detector
          manufacturer_ = "FLIR";
          instrument_model_ = "identiFINDER";
        }else if( UtilityFunctions::icontains(line, "Raider") )
        {
          detector_type_ = kMicroRaiderDetector;
          instrument_model_ = "MicroRaider";
          manufacturer_ = "FLIR";
        }
      }//if( info_pos != string::npos )
    }else if( istarts_with( line, "DetectorType" ) )
    {//DetectorType        : NaI
      if( info_pos != string::npos )
         detector_type = line.substr(info_pos);
    }
    //"GammaDetector" now included in ns_iaea_comment_labels[].
    //else if( istarts_with( line, "GammaDetector" ) )
    //{
      //ex: "GammaDetector           : NaI 35x51"
      //if( info_pos != string::npos )
        //remarks_.push_back( "Gamma Detector: " + line.substr(info_pos) );
        //meas->detector_type_ = line.substr(info_pos);
    //}
    //"NeutronDetector" now included in ns_iaea_comment_labels[].
    //else if( istarts_with( line, "NeutronDetector" ) )
    //{
      //if( info_pos != string::npos )
        //remarks_.push_back( "Neutron Detector: " + line.substr(info_pos) );
    //}
    else if( istarts_with( line, "XUnit" ) )
    {//XUnit        : keV
      if( info_pos != string::npos
          && !istarts_with( line.substr(info_pos), "keV") )
        cerr << SRC_LOCATION << "\n\t" << "Unexpected x-unit: "
             << line.substr(info_pos) << endl;
    }else if( istarts_with( line, "YUnit" ) ) //        :
    {
    }else if( istarts_with( line, "Length" ) )
    {//Length       : 1024
      if( info_pos != string::npos )
        Length = atoi( line.c_str() + info_pos );
    }else if( istarts_with( line, "SubSpcNum" ) )
    {//SubSpcNum    : 1
      int SubSpcNum = 1;
      if( info_pos != string::npos )
        SubSpcNum = atoi( line.c_str() + info_pos );
      if( SubSpcNum > 1 )
      {
        const string msg = SRC_LOCATION + "\n\tASCII Spc files only support "
                           "reading files with one spectrum right now";
        throw std::runtime_error( msg );
      }
    }else if( istarts_with( line, "StartSubSpc" ) )
    {//StartSubSpc  : 0
    }else if( istarts_with( line, "StopSubSpc" ) )
    {//StopSubSpc   : 0
    }else if( istarts_with( line, "Realtime" ) )
    {//Realtime     : 300.000
      if( info_pos != string::npos )
        meas->real_time_ = static_cast<float>( atof( line.c_str() + info_pos ) );
    }else if( istarts_with( line, "Livetime" )
              || istarts_with( line, "Liveime" )
              || istarts_with( line, "Lifetime" ) )
    {//Livetime     : 300.000
      if( info_pos != string::npos )
        meas->live_time_ = static_cast<float>( atof( line.c_str() + info_pos ) );
    }else if( istarts_with( line, "Deadtime" ) )
    {//Deadtime     : 0.000
    }else if( istarts_with( line, "FastChannel" ) )
    {//FastChannel  : 69008
    }else if( istarts_with( line, "Starttime" ) )
    {//Starttime    : 28.08.2012 16:12:26
      if( info_pos != string::npos )
        meas->start_time_ = time_from_string( line.substr( info_pos ).c_str() );
    }else if( istarts_with( line, "Stoptime" ) )
    {//Stoptime     : 28.08.2012 16:17:25
//      if( info_pos != string::npos )
//         meas-> = time_from_string( line.substr( info_pos ).c_str() );
    }else if( istarts_with( line, "NeutronCounts" )
              || istarts_with( line, "SumNeutrons" ) )
    {//NeutronCounts         : 0
      const float num_neut = static_cast<float>( atof( line.substr( info_pos ).c_str() ) );
      if( info_pos != string::npos && meas->neutron_counts_.empty() )
         meas->neutron_counts_.push_back( num_neut );
      else if( info_pos != string::npos )
        meas->neutron_counts_[0] += num_neut;
      meas->neutron_counts_sum_ += num_neut;
      meas->contained_neutron_ = true;
    }
//    FWHMCCoeff            : a=0.000000000E+000 b=0.000000000E+000 c=0.000000000E+000 d=0.000000000E+000'
    else if( starts_with( line, "CalibCoeff" ) )
    {//CalibCoeff   : a=0.000000000E+000 b=0.000000000E+000 c=3.000000000E+000 d=0.000000000E+000
      float a = 0.0f, b = 0.0f, c = 0.0f, d = 0.0f;
      const size_t apos = line.find( "a=" );
      const size_t bpos = line.find( "b=" );
      const size_t cpos = line.find( "c=" );
      const size_t dpos = line.find( "d=" );
      if( apos < (line.size()-2) )
        a = static_cast<float>( atof( line.c_str() + apos + 2 ) );
      if( bpos < (line.size()-2) )
        b = static_cast<float>( atof( line.c_str() + bpos + 2 ) );
      if( cpos < (line.size()-2) )
        c = static_cast<float>( atof( line.c_str() + cpos + 2 ) );
      if( dpos < (line.size()-2) )
        d = static_cast<float>( atof( line.c_str() + dpos + 2 ) );
//      if( a!=0.0 || b!=0.0 || d!=0.0 )
//        cerr << SRC_LOCATION << "\n\tUnknown calibration coefficient meaning -"
//             << " wcjohns should check into this although he thinks it should "
//             << "be okay" << endl;

      meas->energy_calibration_model_ = Measurement::Polynomial;
      meas->calibration_coeffs_.resize( 2 );
      meas->calibration_coeffs_[0] = b;
      meas->calibration_coeffs_[1] = c;
    }else if( istarts_with( line, "NuclideID1" )
             || istarts_with( line, "NuclideID2" )
             || istarts_with( line, "NuclideID3" )
             || istarts_with( line, "NuclideID4" ) )
    {
      //"8 Annih. Rad."
      //"- Nuc. U-233"
      //"5 NORM K-40"
      //"- Ind.Ir-192s"
      if( info_pos != string::npos )
      {
        if( !analysis )
          analysis = std::make_shared<DetectorAnalysis>();
        
        DetectorAnalysisResult result;
        
        string info = line.substr(info_pos);
        string::size_type delim = info.find_first_of( ' ' );
        if( delim == 1 && (std::isdigit(info[0]) || info[0]=='-') )
        {
          result.id_confidence_ = info.substr(0,delim);
          info = info.substr(delim);
          UtilityFunctions::trim( info );
          delim = info.find_first_of( " ." );
          
          string nuctype = info.substr( 0, delim );
          
          if( UtilityFunctions::istarts_with(nuctype, "Ann")
             || UtilityFunctions::istarts_with(nuctype, "Nuc")
             || UtilityFunctions::istarts_with(nuctype, "NORM")
             || UtilityFunctions::istarts_with(nuctype, "Ind")
             || UtilityFunctions::istarts_with(nuctype, "Cal")
             || UtilityFunctions::istarts_with(nuctype, "x")
             || UtilityFunctions::istarts_with(nuctype, "med")
             || UtilityFunctions::istarts_with(nuctype, "cos")
             || UtilityFunctions::istarts_with(nuctype, "bac")
             || UtilityFunctions::istarts_with(nuctype, "TENORM")
             || UtilityFunctions::istarts_with(nuctype, "bre")
             //any others? this list so far was just winged
             )
          {
            result.nuclide_type_ = nuctype;
            result.nuclide_ = info.substr( delim );
            UtilityFunctions::trim( result.nuclide_ );
            
            if(!result.nuclide_.empty() && result.nuclide_[0]=='.' )
            {
              result.nuclide_ = result.nuclide_.substr(1);
              UtilityFunctions::trim( result.nuclide_ );
              result.nuclide_type_ += ".";
            }
            
            result.remark_ = line.substr(info_pos); //just in case
          }else
          {
            //Leaving bellow line in because I only tested above parsing on a handfull of files (20161010).
            cerr << "Unknown radiation type in ana  result: '" << result.nuclide_type_ << "'" << endl;
            result.nuclide_ = line.substr(info_pos);
          }
        }else
        {
          result.nuclide_ = line.substr(info_pos);
        }
        
        
        analysis->results_.push_back( result );
      }//if( info_pos != string::npos )
    }else if( istarts_with( line, "IDLibrary" ) )
    {
      //Comes in files with "NuclideID1" and "NuclideID2" lines, after all the
      //  nuclides.
      if( !analysis )
        analysis = std::make_shared<DetectorAnalysis>();
      analysis->remarks_.push_back( "Library: " + line.substr(info_pos) );
    }else if( istarts_with( line, "SpectrumText" ) )
    {//SpectrumText : 0

    }else if( istarts_with( line, "SerialNumber" ) )
    {
      if( info_pos != string::npos )
        instrument_id_ = line.substr(info_pos);
    }else if( istarts_with( line, "UUID" ) )
    {
      if( info_pos != string::npos )
        uuid_ = line.substr(info_pos);
    }else if( istarts_with( line, "Manufacturer" ) )
    {
      if( info_pos != string::npos )
        manufacturer_ = line.substr(info_pos);
    }else if( istarts_with( line, "ModelNumber" ) )
    {
      if( info_pos != string::npos )
        instrument_model_ = line.substr(info_pos);
    }else if( istarts_with( line, "OperatorInformation" ) )
    {
      measurment_operator_ = line.substr(info_pos);
    }else if( istarts_with( line, "GPSValid" ) )
    {
      if( UtilityFunctions::icontains(line, "no") )
      {
        meas->latitude_ = -999.9;
        meas->longitude_= -999.9;
      }
    }else if( istarts_with( line, "GPS" ) )
    {
      line = line.substr(info_pos);
      const string::size_type pos = line.find( '/' );
      if( pos != string::npos )
      {
        for( size_t i = 0; i < line.size(); ++i )
          if( !isalnum(line[i]) )
            line[i] = ' ';
        
        string latstr = line.substr(0,pos);
        string lonstr = line.substr(pos+1);
        trim( latstr );
        trim( lonstr );
        
        meas->latitude_ = ortecLatOrLongStrToFlt( latstr );
        meas->longitude_ = ortecLatOrLongStrToFlt( lonstr );
      }else
      {
        cerr << SRC_LOCATION << ": couldnt split lat lon" << endl;
      }
    }else if( istarts_with( line, "DeviceId" ) )
    {
      instrument_id_ = line.substr(info_pos);
      trim( instrument_id_ );
    }else if( istarts_with( line, "Nuclide0" )
              || istarts_with( line, "Nuclide1" )
              || istarts_with( line, "Nuclide2" )
              || istarts_with( line, "Nuclide3" ) )
    {
      //some identiFINDER 2 LGH detectors makes it here.
      
      const istream::pos_type currentpos = input.tellg();
      
      //"Nuclide0" line is sometimes followed by "Strength0", "Class0",
      //  and "Confidence0" lines, so lets try and grab them.
      string strengthline, classline, confidenceline;
      try
      {
        UtilityFunctions::safe_get_line( input, strengthline, max_len );
        UtilityFunctions::safe_get_line( input, classline, max_len );
        UtilityFunctions::safe_get_line( input, confidenceline, max_len );
        
        
        const size_t strength_colonpos = strengthline.find(':');
        const size_t strength_info_pos = strengthline.find_first_not_of(": ", strength_colonpos);
        
        const size_t class_colonpos = classline.find(':');
        const size_t class_info_pos = classline.find_first_not_of(": ", class_colonpos);
        
        const size_t conf_colonpos = confidenceline.find(':');
        const size_t conf_info_pos = confidenceline.find_first_not_of(": ", conf_colonpos);
        
        if( !UtilityFunctions::istarts_with( strengthline, "Strength" )
            || !UtilityFunctions::istarts_with( classline, "Class" )
            || !UtilityFunctions::istarts_with( confidenceline, "Confidence" )
            || class_info_pos == string::npos
            || strength_info_pos == string::npos
            || conf_info_pos == string::npos )
          throw runtime_error( "" );
        
        strengthline = strengthline.substr( strength_info_pos );
        classline = classline.substr( class_info_pos );
        confidenceline = confidenceline.substr( conf_info_pos );
      }catch(...)
      {
        input.seekg( currentpos );
        strengthline.clear();
        classline.clear();
        confidenceline.clear();
      }
      
      
      if( !analysis )
        analysis = std::make_shared<DetectorAnalysis>();
      
      DetectorAnalysisResult result;
      result.nuclide_ = line.substr(info_pos);
      result.nuclide_type_ = classline;
      result.id_confidence_ = confidenceline;
      if(!strengthline.empty())
        result.remark_ = "Strength " + strengthline;
      
      analysis->results_.push_back( result );
    }else if( line.length() && isdigit(line[0]) )
    {
      ShrdFVecPtr channel_data( new vector<float> );
      
      input.seekg( sol_pos, ios::beg );
      while( UtilityFunctions::safe_get_line( input, line ) )
      {
        trim(line);
        
        if( line.empty() && Length>=0 && channel_data->size()==Length )
        {
          //ref8MLQDKLR3E seems to have a bunch of extra zeros at the end of the
          //  file (after a line break), so lets deal with this in a way that we
          //  can still try to enforce Length==channel_data->size() at the end
          if( !input.eof() )
          {
            istream::pos_type pos;
            do
            {
              pos = input.tellg();
              trim( line );
              if(!line.empty() && (line[0]<'0' || line[0]>'9') )
              {
                input.seekg( pos, ios::beg );
                break;
              }
            }while( UtilityFunctions::safe_get_line( input, line ) );
          }//if( not at the end of the file )
          
          break;
        }//if( we hit an empty line, and weve read the expected number of channels )
        
        if(!line.empty() && (line[0]<'0' || line[0]>'9') )
          break;
        
        vector<float> linefloats;
        UtilityFunctions::split_to_floats( line.c_str(), line.length(), linefloats );
        channel_data->insert( channel_data->end(), linefloats.begin(), linefloats.end() );
      }//while( UtilityFunctions::safe_get_line( input, line ) )
      
      if( size_t(Length) != channel_data->size() )
      {
        bool isPowerOfTwo = ((Length != 0) && !(Length & (Length - 1)));
        if( isPowerOfTwo && Length >= 1024 && size_t(Length) < channel_data->size() )
        {
          channel_data->resize( Length );
        }else if( Length > 0 )
        {
          stringstream msg;
          msg << SRC_LOCATION << "\n\tExpected to read " << Length
              << " channel datas, but instead read " << channel_data->size();
          throw std::runtime_error( msg.str() );
        }//if( Length > 0 && size_t(Length) != channel_data->size() )
      }//if( size_t(Length) != channel_data->size() )
      
      
      meas->gamma_counts_ = channel_data;
      for( const float a : *channel_data )
        meas->gamma_count_sum_ += a;
    }
    else
    {
      if( !linenum && line.length() )
      {
        for( size_t i = 0; i < line.size(); ++i )
          if( (line[i] & 0x80) )
            throw runtime_error( "Unknown tag and non-ascii character in first non-empty line" );
      }
      
      if( UtilityFunctions::istarts_with(line, "TSA,") )
        throw runtime_error( "This is probably a TSA file, not a Ascii Spc" );
      
      ++nnotrecognized;
      if( nnotrecognized > 15 && nnotrecognized >= linenum )
        throw runtime_error( "To many unregognized begining lines" );
     
#if(PERFORM_DEVELOPER_CHECKS)
      cerr << "Warning: MeasurementInfo::load_from_iaea_spc(...):  I didnt recognize line: '"
           << line << "'" << endl;
#endif
    }//if / else to figure out what this line cooresponds to
    
    ++linenum;
  }//while( input.good() )
  }catch( std::exception & )
  {
//    cerr  << SRC_LOCATION << "caught: " << e.what() << endl;    
    
    reset();
    input.clear();
    input.seekg( orig_pos, ios::beg );
    return false;
  }

  
  //identiFINDER 2 NGH spectrum files will have spectrum number as their UUID,
  //  so to create a bit more unique UUID, lets add in the serial number to the
  //  UUID, like in the other identiFINDER formats.
  if(!uuid_.empty() && uuid_.size() < 5 && !instrument_id_.empty())
    uuid_ = instrument_id_ + "/" + uuid_;
  
  if( meas->gamma_counts_->size() < 9 )
  {
    reset();
//    cerr << "MeasurementInfo::load_from_iaea_spc(...): did not read any spectrum info"
//         << endl;
    return false;
  }//if( meas->gamma_counts_->empty() )

  measurements_.push_back( meas );

  detectors_analysis_ = analysis;
  
  //A temporary message untile I debug detector_type_ a little more
  if( icontains(instrument_model_,"identiFINDER")
      && ( (icontains(instrument_model_,"2") && !icontains(instrument_model_,"LG")) || icontains(instrument_model_,"NG")))
     detector_type_ = kIdentiFinderNGDetector;
  else if( icontains(detector_type,"La") && !detector_type.empty())
  {
    cerr << "Has " << detector_type << " is this a LaBr3? Cause I'm assuming it is" << endl;
    //XXX - this doesnt actually catch all LaBr3 detectors
    detector_type_ = kIdentiFinderLaBr3Detector;
  }else if( icontains(instrument_model_,"identiFINDER") && icontains(instrument_model_,"LG") )
  {
    cout << "Untested kIdentiFinderLaBr3Detector association!" << endl;
    detector_type_ = kIdentiFinderLaBr3Detector;
  }else if( icontains(instrument_model_,"identiFINDER") )
  {
    detector_type_ = kIdentiFinderDetector;
  }

//  if( detector_type_ == kUnknownDetector )
//  {
//    cerr << "I couldnt find detector type for ASCII SPC file" << endl;
//  }

  cleanup_after_load();

  return true;
}//bool load_from_iaea_spc( std::istream &input )



bool MeasurementInfo::write_ascii_spc( std::ostream &output,
                                     std::set<int> sample_nums,
                                     const std::set<int> &det_nums ) const
{
  try
  {
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( sample_nums.empty() )
    sample_nums = sample_numbers_;
  
  const size_t ndet = detector_numbers_.size();
  vector<bool> detectors( ndet, true );
  if( !det_nums.empty() )
  {
    for( size_t i = 0; i < ndet; ++i )
      detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
  }//if( det_nums.empty() )
  
  
  MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
  
  if( !summed || !summed->gamma_counts() || summed->gamma_counts()->empty() )
    return false;
  
  if(!summed->title().empty())
    output << pad_iaea_prefix( "SpectrumName" ) << summed->title() << "\r\n";
  else
    output << pad_iaea_prefix( "SpectrumName" ) << filename_ << "\r\n";
  
  output << pad_iaea_prefix( "XUnit" ) << "keV\r\n";
  output << pad_iaea_prefix( "YUnit" ) << "\r\n";
  output << pad_iaea_prefix( "Length" ) << summed->gamma_counts_->size() << "\r\n";
  output << pad_iaea_prefix( "SubSpcNum" )   << "1\r\n";
  output << pad_iaea_prefix( "StartSubSpc" ) << "0\r\n";
  output << pad_iaea_prefix( "StopSubSpc" )  << "0\r\n";
  
  int ncomment = 0;
  bool printedFWHMCCoeff = false;
  for( const string &remark : remarks_ )
  {
    bool used = false;
    for( const char * const label : ns_iaea_comment_labels )
    {
      const string prefix = label + string(" : ");
      if( UtilityFunctions::istarts_with(remark, prefix) )
      {
        output << pad_iaea_prefix(label) << remark.substr( prefix.size() ) << "\r\n";
        printedFWHMCCoeff |= UtilityFunctions::iequals(label,"FWHMCCoeff");
        used = true;
        break;
      }
    }//for( const char * const label : ns_iaea_comment_labels )
    
    if( !used )
    {
      ++ncomment;
      output << pad_iaea_prefix("Comment") << remark << "\r\n";
    }
  }//for( const string &remark : remarks_ )
  
  if( !ncomment )
    output << pad_iaea_prefix("Comment") << "\r\n";
  
  char buffer[256];
  if( summed->real_time_ > 0.0f )
  {
    snprintf( buffer, sizeof(buffer), "%.3f", summed->real_time_ );
    output << pad_iaea_prefix( "Realtime" ) << buffer << "\r\n";
  }
  
  if( summed->live_time_ > 0.0f )
  {
    snprintf( buffer, sizeof(buffer), "%.3f", summed->live_time_ );
    output << pad_iaea_prefix( "Livetime" ) << buffer << "\r\n";
  }
  
  if( (summed->real_time_ > 0.0f) && (summed->live_time_ > 0.0f) )
  {
    snprintf( buffer, sizeof(buffer), "%.3f", (summed->real_time_ - summed->live_time_) );
    output << pad_iaea_prefix( "Deadtime" ) << buffer << "\r\n";
  }
  
  //I dont know what FastChannel is
  //output << pad_iaea_prefix( "FastChannel" ) << "3229677" << "\r\n";
  
  for( const pair<std::string,std::string> &cmpnt : component_versions_ )
  {
    for( const char * const label : ns_iaea_version_labels )
    {
      if( cmpnt.first == label )
      {
        output << pad_iaea_prefix(cmpnt.first) << cmpnt.second << "\r\n";
        break;
      }//if( cmpnt.first == label )
    }
  }//for( const pair<std::string,std::string> &cmpnt : component_versions_ )
  
  
  if( !summed->start_time_.is_special() )
  {
    output << pad_iaea_prefix( "Starttime" ) << print_to_iaea_datetime(summed->start_time_) << "\r\n";

    //Add stop time if we
    if( (sample_nums.size()==1 && det_nums.size()==1) )
    {
      float intsec, fracsec;
      fracsec = std::modf( summed->real_time_, &intsec );
    
      boost::posix_time::ptime endtime = summed->start_time_
                                          + boost::posix_time::seconds( static_cast<int>(intsec) )
                                          + boost::posix_time::microseconds( static_cast<int>( floor((1.0e6f * fracsec) + 0.5f) ) );
    
      output << pad_iaea_prefix( "StopTime" ) << print_to_iaea_datetime( endtime ) << "\r\n";
    }//
  }//if( !summed->start_time_.is_special() )
  
  
  if( summed->contained_neutron_ )
    output << pad_iaea_prefix( "NeutronCounts" ) << static_cast<int>(floor(summed->neutron_counts_sum_ + 0.5)) << "\r\n";
  
  vector<float> calcoefs;
  
  if( summed->energy_calibration_model_ == Measurement::Polynomial )
    calcoefs = summed->calibration_coeffs_;
  else if( summed->energy_calibration_model_ == Measurement::FullRangeFraction )
    calcoefs = fullrangefraction_coef_to_polynomial( summed->calibration_coeffs_, summed->gamma_counts_->size() );
  
  const size_t ncoef = calcoefs.size();
  const float a = 0.0f, b = (ncoef ? calcoefs[0] : 0.0f),
              c = (ncoef>1 ? calcoefs[1] : 0.0f), d = (ncoef>2 ? calcoefs[2] : 0.0f);
  
  snprintf( buffer, sizeof(buffer), "a=%.9e b=%.9e c=%.9e d=%.9e", a, b, c, d );
  output << pad_iaea_prefix( "CalibCoeff" ) << buffer << "\r\n";
  
  //Not sure why this line always appears in files, but whatever
  if( !printedFWHMCCoeff )
    output << pad_iaea_prefix( "FWHMCCoeff" ) << "a=0.000000000E+000 b=0.000000000E+000 c=0.000000000E+000 d=0.000000000E+000\r\n";
  
  
  if(!instrument_id_.empty())
  {
    //We see two variants of how the serial number is specified, so lets put
    //  both into the file incase an analysis program only looks for one of them
    output << pad_iaea_prefix( "SerialNumber" ) << instrument_id_ << "\r\n";
    output << pad_iaea_prefix( "DeviceId" ) << instrument_id_ << "\r\n";
  }
  
  string uuid = uuid_;
  if(!instrument_id_.empty() && istarts_with( uuid, (instrument_id_+"/") ) )
    uuid = uuid.substr(instrument_id_.size()+1);
  if(!uuid.empty())
    output << pad_iaea_prefix( "UUID" ) << uuid.substr(instrument_id_.size()+1) << "\r\n";
  
  if(!manufacturer_.empty())
    output << pad_iaea_prefix( "Manufacturer" ) << manufacturer_ << "\r\n";
  
  if(!instrument_model_.empty())
    output << pad_iaea_prefix( "ModelNumber" ) << instrument_model_ << "\r\n";
  
  if(!measurment_operator_.empty())
    output << pad_iaea_prefix( "OperatorInformation" ) << measurment_operator_ << "\r\n";
  
  if( summed->has_gps_info() )
  {
    output << pad_iaea_prefix( "GPSValid" ) << "yes\r\n";
    //Should probably put into degree, minute, second notation, but some other time...
    output << pad_iaea_prefix( "GPS" ) << summed->latitude_ << "," << summed->longitude_ << "\r\n";
  }//if( summed->has_gps_info() )
  
  
  if( detectors_analysis_ && !detectors_analysis_->results_.empty())
  {
    for( size_t i = 0; i < detectors_analysis_->results_.size(); ++i )
    {
      const DetectorAnalysisResult &res = detectors_analysis_->results_[i];
      
      //We see two ways analysis results are conveyed in SPC files; lets make an
      //  attempt at having the output be consistent with the input, in terms
      //  of SPC files.  This of course will be a inconsistent for converting
      //  other file formats to SPC, but such is life
      if(!res.nuclide_.empty() && !res.nuclide_type_.empty())
      {
        const string postfix = string("") + char('0' + i);
        output << pad_iaea_prefix( "Nuclide" + postfix ) << res.nuclide_ << "\r\n";
        if( UtilityFunctions::istarts_with(res.remark_, "Strength ") )
          output << pad_iaea_prefix( "Strength" + postfix ) << res.nuclide_ << "\r\n";
        else
          output << pad_iaea_prefix( "Strength" + postfix ) << "\r\n";
        output << pad_iaea_prefix( "Class" + postfix ) << res.nuclide_type_ << "\r\n";
        output << pad_iaea_prefix( "Confidence" + postfix ) << res.id_confidence_ << "\r\n";
      }else if(!res.nuclide_.empty())
      {
        const string postfix = string("") + char('1' + i);
        output << pad_iaea_prefix( "NuclideID" + postfix ) << res.nuclide_ << "\r\n";
      }
    }//
  }//if( have analsysi results to output )
  
  
  //I've only seen the SpectrumText line either empty, having only a '0'...
  output << pad_iaea_prefix( "SpectrumText" ) << "\r\n";
  
  const vector<float> &counts = *summed->gamma_counts_;
  output << counts[0];
  for( size_t i = 1; i < counts.size(); ++i )
    output << (((i%8)==0) ? "\r\n" : ",") << counts[i];
  output << "\r\n";
  
    /*
     //Other tags found in identiFINDER SPE files:
     $RT:
     12
     $DT:
     9
     $SPEC_INTEGRAL:
     2769
     $FLIR_DATASET_NUMBER:
     284
     $FLIR_GAMMA_DETECTOR_INFORMATION:
     SodiumIodide
     Cesium137
     NaI 35x51
     
     $FLIR_NEUTRON_DETECTOR_INFORMATION:
     ?He tube 3He3/608/15NS
     $FLIR_SPECTRUM_TYPE:
     Measurement
     UserControlled
     $FLIR_DOSE_RATE_SWMM:
     12
     1.301
     0.062
     0.141
     $FLIR_NEUTRON_SWMM:
     12
     0.000
     0.000
     0.000
     $FLIR_REACHBACK:
     
     */
    
  //Still need to fix up DetectorType, and deal with instrument model coorectly
  /*
      }else if( istarts_with( line, "DetectorType" ) )
      {//DetectorType        : NaI
        if( info_pos != string::npos )
          detector_type = line.substr(info_pos);
          }
  //A temporary message untile I debug detector_type_ a little more
  if( icontains(instrument_model_,"identiFINDER")
     && ( (icontains(instrument_model_,"2") && !icontains(instrument_model_,"LG")) || icontains(instrument_model_,"NG")))
    detector_type_ = kIdentiFinderNGDetector;
  else if( icontains(detector_type,"La") && detector_type.size() )
  {
    cerr << "Has " << detector_type << " is this a LaBr3? Cause I'm assuming it is" << endl;
    //XXX - this doesnt actually catch all LaBr3 detectors
    detector_type_ = kIdentiFinderLaBr3Detector;
  }else if( icontains(instrument_model_,"identiFINDER") && icontains(instrument_model_,"LG") )
  {
    cout << "Untested kIdentiFinderLaBr3Detector association!" << endl;
    detector_type_ = kIdentiFinderLaBr3Detector;
  }else if( icontains(instrument_model_,"identiFINDER") )
  {
    detector_type_ = kIdentiFinderDetector;
  }
  */
  }catch( std::exception & )
  {
    return false;
  }

  return true;
}//bool write_ascii_spc(...)


bool MeasurementInfo::write_binary_spc( std::ostream &output,
                                    const MeasurementInfo::SpcBinaryType type,
                                    std::set<int> sample_nums,
                                    const std::set<int> &det_nums ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( sample_nums.empty() )
    sample_nums = sample_numbers_;
  
  const size_t ndet = detector_numbers_.size();
  vector<bool> detectors( ndet, true );
  if( !det_nums.empty() )
  {
    for( size_t i = 0; i < ndet; ++i )
      detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
  }//if( det_nums.empty() )
  
  //const size_t initialpos = output.tellp();
  
  MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
  
  if( !summed || !summed->gamma_counts() )
    return false;
  
  //  if( summed->energy_calibration_model() != Measurement::FullRangeFraction
  //     && summed->energy_calibration_model() != Measurement::Polynomial )
  //    return false;
  
  //require the number of channels to be a power of two.
  const uint16_t n_channel = static_cast<uint16_t>( summed->gamma_counts()->size() ); // Number of Channels;
  const bool isPowerOfTwo = ((n_channel != 0) && !(n_channel & (n_channel - 1)));
  if( !isPowerOfTwo )
    return false;
  
  //see http://www.ortec-online.com/download/ortec-software-file-structure-manual.pdf
  const int16_t wINFTYP = 1; // Must be 1
  const int16_t wFILTYP = (type==IntegerSpcType ? 1 : 5);
  const int16_t wSkip1[2] = { 0, 0 } ;
  const int16_t wACQIRP = 3; // Acquisition information record pointer
  const int16_t wSAMDRP = 4; // Sample description record pointer
  const int16_t wDETDRP = 5; // Detector description record pointer
  const int16_t wSKIP2[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  const int16_t wCALDES = 0; //8; // Calibration description record pointer
  const int16_t wCALRP1 = 6; // First calibration data record pointer
  const int16_t wCALRP2 = 0; //7; // Second calibration data record pointer
  const int16_t wEFFPRP = 0; // Efficiency pairs record pointer (first record)
  const int16_t wROIRP1 = 0; // Record number of the first of two ROI recs
  const int16_t wEPRP = 0;   // Energy pairs record pointer
  const int16_t wEPN = 0;    // Number of energy pairs records
  const int16_t wSkip3[6] = { 0, 0, 0, 0, 0, 0 };
  const int16_t wEFFPNM = 0; // Number of efficiency pairs records  256
  const int16_t wSPCTRP = 9; //was 25 for example file //could probably be wCALDES // pointer to spectrum
  //n_channel
  //We can fit 32 floats per 128 byte record
  const int16_t wSPCRCN = (n_channel / 32) + (((n_channel%32)!=0) ? 1 : 0); // Number records for the spectrum
  
  const int16_t wABSTCHN = 0;// Physical start channel for data
  
  float sACQTIM = 0.0; //10616.5 is 2014-Sep-19 12:14:57  // Date and time acquisition start in DECDAY format
  double dACQTI8 = 0.0; //10616.5 // Date and time as double precision DECDAY
  
  if( !summed->start_time().is_special() )
  {
    const boost::posix_time::ptime &startime = summed->start_time();
    const boost::gregorian::date epic_date( 1979, boost::gregorian::Jan, 1 );
    const boost::gregorian::days daydiff = startime.date() - epic_date;
    const double dayfrac = startime.time_of_day().total_microseconds() / (24.0*60.0*60.0*1.0E6);
    dACQTI8 = daydiff.days() + dayfrac;
    sACQTIM = static_cast<float>( dACQTI8 );
  }//if( !summed->start_time().is_special() )
  
  const int16_t wSkip4[4] = { 0, 0, 0, 0 };
  const int16_t wCHNSRT = 0; // Start channel number
  const float sRLTMDT = summed->real_time(); // Real time in seconds
  const float sLVTMDT = summed->live_time(); // Live time in seconds
  const int16_t wSkip50 = 0;
  const int16_t framRecords = 0;//     Pointer to FRAM records
  const int16_t TRIFID = 0; // I*2     Pointer to TRIFID records
  const int16_t NaI = 0; //I*2     Pointer to NaI records
  const int16_t Location = 0;// I*2     Pointer to Location records
  const int16_t MCSdata = 0; //I*2     Number of channels of MCS data appended to the histogram data in this file
  const int16_t expansionHeader = 2;//     Pointer to expansion header record (i.e. second header record)
  const int16_t reserved[5] = { 0, 0, 0, 0, 0 }; // 57-62                 Reserved (must be 0)
  const float RRSFCT = 0.0;//  R*4     Total random summing factor
  const int16_t zeroword = 0;
  //word number
  writeBinaryData( output, wINFTYP );            //1
  writeBinaryData( output, wFILTYP );            //2
  output.write( (const char *)&wSkip1[0], 2*2 ); //4
  writeBinaryData( output, wACQIRP );            //5
  writeBinaryData( output, wSAMDRP );            //6
  writeBinaryData( output, wDETDRP );            //7
  output.write( (const char *)&wSKIP2[0], 2*9 ); //8
  writeBinaryData( output, wCALDES );            //17
  writeBinaryData( output, wCALRP1 );            //18
  writeBinaryData( output, wCALRP2 );            //19
  writeBinaryData( output, wEFFPRP );            //20
  writeBinaryData( output, wROIRP1 );            //21
  writeBinaryData( output, wEPRP );              //22
  writeBinaryData( output, wEPN );               //23
  output.write( (const char *)&wSkip3[0], 2*6 ); //24
  writeBinaryData( output, wEFFPNM );            //30
  writeBinaryData( output, wSPCTRP );            //31
  writeBinaryData( output, wSPCRCN );            //32
  writeBinaryData( output, n_channel );          //33
  writeBinaryData( output, wABSTCHN );           //34
  writeBinaryData( output, sACQTIM );            //35
  writeBinaryData( output, dACQTI8 );            //37
  
  output.write( (const char *)&wSkip4[0], 2*4 ); //41
  writeBinaryData( output, wCHNSRT );            //45
  writeBinaryData( output, sRLTMDT );            //46
  writeBinaryData( output, sLVTMDT );            //48
  writeBinaryData( output, wSkip50 );            //50
  writeBinaryData( output, framRecords );        //51
  
  //writeBinaryData( output, zeroword );           //52
  
  writeBinaryData( output, TRIFID );             //53
  writeBinaryData( output, NaI );                //54
  
  writeBinaryData( output, Location );           //
  writeBinaryData( output, MCSdata );            //
  writeBinaryData( output, expansionHeader );    //
  output.write( (const char *)&reserved[0], 2*5 ); //55
  
  //output.write( (const char *)&reserved[0], 2*3 ); //60
  
  //  output.write( (const char *)&reserved[0], 2*5 ); //
  writeBinaryData( output, RRSFCT ); //63
  
  //20160915: we're actually at pos 126 right now, so lets write in
  writeBinaryData( output, zeroword );
  
  const uint8_t zero_byte = 0;
  
  //Write expansion header information
  size_t pos = 128;
  size_t poswanted = (expansionHeader-1)*128;
  while( pos < poswanted )
  {
    ++pos;
    output.write( (const char *)&zero_byte, 1 );
    pos += writeBinaryData( output, zero_byte );
  }
  
  int16_t firstReportPtr = 0;
  
  {
    const int16_t recordID = 111;
    const int16_t gpsPointer = 0;  //I havent been able to reliable decode files with a GPS record, so we wont write this record
    
    if( summed->contained_neutron() || summed->has_gps_info() )
      firstReportPtr = wSPCTRP + wSPCRCN;
    
    pos += writeBinaryData( output, recordID );
    pos += writeBinaryData( output, gpsPointer );
    pos += writeBinaryData( output, firstReportPtr );
  }
  
  //write Acquisition information record pointer
  //1      Default spectrum file name (stored as 16 ASCII characters)
  //17     Date in the form DD-MMM-YY* (stored as 12 ASCII characters).
  //       The * character should be ignored if it is not a "1". If it is a "1",
  //       it indicates the data is after the year 2000.
  //29     Time in the form HH:MM:SS (stored as 10 ASCII characters)
  //39     Live Time rounded to nearest second (stored as 10 ASCII characters)
  //49     Real Time rounded to nearest second (stored as 10 ASCII characters)
  //59–90: Reserved
  //91     Start date of sample collection (10 ASCII characters)
  //103    Start time of sample collection (8 ASCII characters)
  //111    Stop date of sample collection (10 ASCII characters
  //121    Stop time of sample collection (8 ASCII characters)
  
  poswanted = (wACQIRP-1)*128;
  assert( (expansionHeader == 0) || (wACQIRP > expansionHeader) );
  while( pos < poswanted )
    pos += writeBinaryData( output, zero_byte );
  
  const char *defaultname = 0;
  switch( detector_type_ )
  {
    case kGR135Detector:          case kIdentiFinderDetector:
    case kIdentiFinderNGDetector: case kIdentiFinderLaBr3Detector:
    case kSAIC8Detector:          case kFalcon5000:
    case kUnknownDetector:        case kMicroRaiderDetector:
    case kRsi701: case kRsi705: case kAvidRsi: case kSam940LaBr3: case kSam940:
    case kOrtecRadEagleNai: case kOrtecRadEagleCeBr2Inch:
    case kOrtecRadEagleCeBr3Inch: case kOrtecRadEagleLaBr:
    case kSam945:
    case kRadHunterNaI: case kRadHunterLaBr3:
      defaultname = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
      break;
      
    case kDetectiveDetector:      case kDetectiveExDetector:
    case kDetectiveEx100Detector: case kOrtecIDMPortalDetector:
    case kMicroDetectiveDetector:
      defaultname = "DetectiveEX.SPC";
      break;
  }//switch( detector_type_ )
  
  pos += 16;
  output.write( defaultname, 16 );
  
  
  string datestr;
  if( summed->start_time().date().is_special() )
  {
    datestr = "01-Jan-001\0\0\0";
  }else
  {
    const int daynum = summed->start_time().date().day();
    if( daynum < 10 )
      datestr += "0";
    datestr += std::to_string(daynum);
    datestr += "-";
    switch( summed->start_time().date().month() )
    {
      case boost::gregorian::Jan: datestr += "Jan"; break;
      case boost::gregorian::Feb: datestr += "Feb"; break;
      case boost::gregorian::Mar: datestr += "Mar"; break;
      case boost::gregorian::Apr: datestr += "Apr"; break;
      case boost::gregorian::May: datestr += "May"; break;
      case boost::gregorian::Jun: datestr += "Jun"; break;
      case boost::gregorian::Jul: datestr += "Jul"; break;
      case boost::gregorian::Aug: datestr += "Aug"; break;
      case boost::gregorian::Sep: datestr += "Sep"; break;
      case boost::gregorian::Oct: datestr += "Oct"; break;
      case boost::gregorian::Nov: datestr += "Nov"; break;
      case boost::gregorian::Dec: datestr += "Dec"; break;
      case boost::gregorian::NotAMonth: case boost::gregorian::NumMonths:
        datestr += "\0\0\0";
        break;
    }//switch( summed->start_time().date().month() )
    
    datestr += "-";
    const int yearnum = summed->start_time().date().year() % 100;
    if( yearnum < 10 )
      datestr += "0";
    datestr += std::to_string(yearnum);
    datestr += (summed->start_time().date().year() > 1999 ? "1" : "0");
    datestr.resize( 13, '\0' );
  }
  
  pos += 12;
  output.write( &datestr[0], 12 );
  
  char timestr[11] = { '\0' };
  const int hournum = summed->start_time().time_of_day().hours();
  const int minutenum = summed->start_time().time_of_day().minutes();
  const int secondnum = summed->start_time().time_of_day().seconds();
  snprintf( timestr, sizeof(timestr), "%02d:%02d:%02d",
           hournum, minutenum, secondnum );
  
  pos += 10;
  output.write( timestr, 10 );
  
  const int rtint = int(floor(summed->live_time()+0.5f));
  const int ltint = int(floor(summed->real_time()+0.5f));
  string ltIntStr = std::to_string(rtint );
  string rtIntStr = std::to_string(ltint);
  ltIntStr.resize( 11, '\0' );
  rtIntStr.resize( 11, '\0' );
  
  pos += 10;
  output.write( &ltIntStr[0], 10 );
  pos += 10;
  output.write( &rtIntStr[0], 10 );
  
  for( int i = 0; i < 32; ++i )
    pos += writeBinaryData( output, zero_byte );
  
  //The start date/time and end date/time of sample collection doesnt seem to
  //  be correct in detector generated SPC file, so we'll just put something
  //  here that is not correct
  const char start_date_of_sample_collection[13] = "\0\0\0\0\0\0\0\0\0\0\0\0";
  pos += 12;
  output.write( start_date_of_sample_collection, 12 );
  const char start_time_of_sample_collection[9] = "10:59:03";
  pos += 8;
  output.write( start_time_of_sample_collection, 8 );
  const char stop_date_of_sample_collection[11] = "25-JAN-081";
  pos += 10;
  output.write( stop_date_of_sample_collection, 10 );
  const char stop_time_of_sample_collection[9] = "10:59:03";
  pos += 8;
  output.write( stop_time_of_sample_collection, 8 );
  
  
  
  //Write Sample description record (only works if input file format was SPC)
  poswanted = (wSAMDRP-1)*128;
  while( pos < poswanted )
    pos += writeBinaryData( output, zero_byte );
  string sampledescrip = summed->title_;
  for( const string &s : remarks_ )
  {
    if( UtilityFunctions::starts_with(s, "Sample Description: ") )
      sampledescrip = " " + s.substr( 20 );
  }

  trim( sampledescrip );
  
  sampledescrip.resize( 128, '\0' );
  pos += 128;
  output.write( &sampledescrip[0], 128 );
  
  
  //Write Detector description record pointer
  poswanted = (wDETDRP-1)*128;
  assert( wDETDRP > wSAMDRP );
  while( pos < poswanted )
    pos += writeBinaryData( output, zero_byte );
  string descrip = instrument_id_;
  descrip.resize( 128, '\0' );
  pos += 128;
  output.write( &descrip[0], 128 );
  
  
  //First calibration data record
  poswanted = (wCALRP1-1)*128;
  assert( wCALRP1 > wDETDRP );
  while( pos < poswanted )
    pos += writeBinaryData( output, zero_byte );
  
  {//begin codeblock to write energy calibration information
    const int16_t wAFIT = 0; // 2 Above knee efficiency calibration type
    const int16_t wBFIT = 0; // 4 Below knee efficiency calibration type
    const int16_t wEFFPRS = 0; // 6 Number of efficiency pairs
    const int16_t wNCH = 0; // 8 number of channels in spectrum
    const float sKNEE = 0.0f; // 12 Detector knee (keV)
    const float sASIG = 0.0f; // 16 2-sigma uncertainty above knee
    const float sBSIG = 0.0f; // 20 2-sigma uncertainty below knee
    float sEC1 = 0.0f; // 24 Energy vs. channel coefficient A
    float sEC2 = 0.0f; // 28 Energy vs. channel coefficient B
    float sEC3 = 0.0f; // 32 Energy vs. channel coefficient C
    const float sFC1 = 0.0f; //FWHM vs. channel coefficient A (actually in most file)
    const float sFC2 = 0.0f; //FWHM vs. channel coefficient B (actually in most file)
    const float sFC3 = 0.0f; //FWHM vs. channel coefficient C (actually in most file)
    const float sPE1 = 0.0f; //Above knee eff. vs. energy coeff A or poly coeff 1
    const float sPE2 = 0.0f; //Above knee eff. vs. energy coeff B or poly coeff 2
    const float sPE3 = 0.0f; //Above knee eff. vs. energy coeff C or poly coeff 3
    const float sSE1 = 0.0f; //Below knee eff. vs. energy coeff A or poly coeff 4
    const float sSE2 = 0.0f; //Below knee eff. vs. energy coeff B or poly coeff 5
    const float sSE3 = 0.0f; //Below knee eff. vs. energy coeff C or poly coeff 6
    const int16_t wFWHTYP = 0; // FWHM type
    const int16_t wRES1 = 0; // reserved
    const int16_t wRES2 = 3; // reserved
    const int16_t wENGPRS = 0; // Number of energy pairs
    const int16_t wDETNUM = 0; //Detector Number
    const int16_t wNBKNEE = 0; // Number of calibration points below knee
    const float sENA2 = 0.0f; // Temp energy calibration
    const float sENB2 = 0.0f; // Temp energy calibration
    const float sENC2 = 0.0f; // Temp energy calibration
    const float sCALUNC = 0.0f; //Calibration source uncertainty
    const float sCALDIF = 0.0f; // Energy calibration difference
    const float sR7 = 0.0f; // Polynomial coefficient 7
    const float sR8 = 0.0f; // Polynomial coefficient 8
    const float sR9 = 0.0f; // Polynomial coefficient 9
    const float sR10 = 0.0f; // Polynomial coefficient 10
    
    
    vector<float> calib_coef = summed->calibration_coeffs();
    if( summed->energy_calibration_model() == Measurement::FullRangeFraction )
      calib_coef = fullrangefraction_coef_to_polynomial( calib_coef, n_channel );
    else if( summed->energy_calibration_model() != Measurement::Polynomial )
      calib_coef.clear();
    
    if( calib_coef.size() > 0 )
      sEC1 = calib_coef[0];
    if( calib_coef.size() > 1 )
      sEC2 = calib_coef[1];
    if( calib_coef.size() > 2 )
      sEC3 = calib_coef[2];
    
    pos += writeBinaryData( output, wAFIT );
    pos += writeBinaryData( output, wBFIT ); // 4 Below knee efficiency calibration type
    pos += writeBinaryData( output, wEFFPRS ); // 6 Number of efficiency pairs
    pos += writeBinaryData( output, wNCH ); // 8 number of channels in spectrum
    pos += writeBinaryData( output, sKNEE ); // 12 Detector knee (keV)
    pos += writeBinaryData( output, sASIG ); // 16 2-sigma uncertainty above knee
    pos += writeBinaryData( output, sBSIG ); // 20 2-sigma uncertainty below knee
    pos += writeBinaryData( output, sEC1 ); // 24 Energy vs. channel coefficient A
    pos += writeBinaryData( output, sEC2 ); // 28 Energy vs. channel coefficient B
    pos += writeBinaryData( output, sEC3 ); // 32 Energy vs. channel coefficient C
    pos += writeBinaryData( output, sFC1 ); //FWHM vs. channel coefficient A (actually in most file)
    pos += writeBinaryData( output, sFC2 ); //FWHM vs. channel coefficient B (actually in most file)
    pos += writeBinaryData( output, sFC3 ); //FWHM vs. channel coefficient C (actually in most file)
    pos += writeBinaryData( output, sPE1 ); //Above knee eff. vs. energy coeff A or poly coeff 1
    pos += writeBinaryData( output, sPE2 ); //Above knee eff. vs. energy coeff B or poly coeff 2
    pos += writeBinaryData( output, sPE3 ); //Above knee eff. vs. energy coeff C or poly coeff 3
    pos += writeBinaryData( output, sSE1 ); //Below knee eff. vs. energy coeff A or poly coeff 4
    pos += writeBinaryData( output, sSE2 ); //Below knee eff. vs. energy coeff B or poly coeff 5
    pos += writeBinaryData( output, sSE3 ); //Below knee eff. vs. energy coeff C or poly coeff 6
    pos += writeBinaryData( output, wFWHTYP ); // FWHM type
    pos += writeBinaryData( output, wRES1 ); // reserved
    pos += writeBinaryData( output, wRES2 ); // reserved
    pos += writeBinaryData( output, wENGPRS ); // Number of energy pairs
    pos += writeBinaryData( output, wDETNUM ); //Detector Number
    pos += writeBinaryData( output, wNBKNEE ); // Number of calibration points below knee
    pos += writeBinaryData( output, sENA2 ); // Temp energy calibration
    pos += writeBinaryData( output, sENB2 ); // Temp energy calibration
    pos += writeBinaryData( output, sENC2 ); // Temp energy calibration
    pos += writeBinaryData( output, sCALUNC ); //Calibration source uncertainty
    pos += writeBinaryData( output, sCALDIF ); // Energy calibration difference
    pos += writeBinaryData( output, sR7 ); // Polynomial coefficient 7
    pos += writeBinaryData( output, sR8 ); // Polynomial coefficient 8
    pos += writeBinaryData( output, sR9 ); // Polynomial coefficient 9
    pos += writeBinaryData( output, sR10 ); // Polynomial coefficient 10
  }//end codeblock to write energy calibration information
  
  
  
  //Second calibration data record
  if( wCALRP2 > 0 )
  {
    //static_assert( wCALRP2 > wCALRP1, "");
    poswanted = (wCALRP2-1)*128;
    while( pos < poswanted )
      pos += writeBinaryData( output, zero_byte );
  }//if( wCALRP2 > 0 )
  
  
  //Calibration description record
  if( wCALDES > 0 )
  {
      static_assert( (wCALRP2 == 0) || (wCALDES > wCALRP2), "");
    
    poswanted = (wCALDES-1)*128;
    while( pos < poswanted )
      pos += writeBinaryData( output, zero_byte );
  }//if( wCALDES > 0 )
  
  
  //Write spectrum information
  poswanted = (wSPCTRP-1)*128;
  while( pos < poswanted )
    pos += writeBinaryData( output, zero_byte );
  
  pos += 4*n_channel;
  const vector<float> &channel_data = *summed->gamma_counts();
  if( type == IntegerSpcType )
  {
    vector<uint32_t> int_channel_data( n_channel );
    for( int16_t i = 0; i < n_channel; ++i )
      int_channel_data[i] = static_cast<uint32_t>( channel_data[i] );  //should we round instead of truncate?
    output.write( (const char *)&int_channel_data[0], 4*n_channel );
  }else
  {
    output.write( (const char *)&channel_data[0], 4*n_channel );
  }//if( file is integer channel data ) / else float data
  
  if( firstReportPtr > 0 )
  {
    assert( firstReportPtr>0 && static_cast<size_t>(128*(firstReportPtr-1)) >= pos );
    while( pos < poswanted )
      pos += writeBinaryData( output, zero_byte );
    
    string information;
    if( summed->has_gps_info() )
    {
      int degrees, minutes, seconds;
      double val = summed->latitude();
      degrees = static_cast<int>( floor( fabs(val) ) );
      val = 60.0*(val - degrees);
      minutes = static_cast<int>( floor( val ) );
      val = 60.0*(val - minutes);
      seconds = static_cast<int>( floor( val + 0.5 ) );
      
      char buffer[128];
      snprintf( buffer, sizeof(buffer), "Latitude %d %d %d %s\n",
               degrees, minutes, seconds, (summed->latitude()>0 ? "N" : "S") );
      information += buffer;
      
      val = summed->longitude();
      degrees = static_cast<int>( floor( fabs(val) ) );
      val = 60.0*(val - degrees);
      minutes = static_cast<int>( floor( val ) );
      val = 60.0*(val - minutes);
      seconds = static_cast<int>( floor( val + 0.5 ) );
      
      snprintf( buffer, sizeof(buffer), "Longitude %d %d %d %s\n",
               degrees, minutes, seconds, (summed->longitude()>0 ? "E" : "W") );
      information += buffer;
    }//if( summed->has_gps_info() )
    
    
    if( summed->contained_neutron() )
    {
      char buffer[256];
      const int nneut = static_cast<int>( floor(summed->neutron_counts_sum()+0.5) );
      snprintf( buffer, sizeof(buffer), "Total neutron counts = %d\n", nneut );
      information += buffer;
      
      for( const string &s : remarks_ )
      {
        if( s.find("Total neutron count time = ") != string::npos )
          information += s + "\n";
      }
    }//if( summed->contained_neutron() )
    
    //Should consider adding: "Total neutron count time = ..."
    
    
    if( information.size() >= 2048 )
      information.resize( 2047 );
    
    const uint16_t ntxtbytes = static_cast<uint16_t>( information.size() + 1 );
    const uint16_t sourcecode = 0;  //I dont actually know what this is for
    pos += writeBinaryData( output, ntxtbytes );
    pos += writeBinaryData( output, sourcecode );
    
    pos += ntxtbytes;
    output.write( information.c_str(), ntxtbytes );
    
    //Advance the file position to the next record start, to keep file size a
    //  mutliple of 128 bytes (not sure if this is needed)
    while( (pos % 128) )
      pos += writeBinaryData( output, zero_byte );
    
    //Should upgrade to include analysis results if they are available
    //    if( !!detectors_analysis_ ){ ... }
    
    //"Found Nuclides"
    //"Suspect Nuclides"
    //"Top Lines"
    //"GPS"  or something like "GPS Location not determined"
    //"Gamma Dose Rate"
    //"Version"
    //"ID Report"
  }//if( firstReportPtr > 0 )
  
  
  return true;
}//bool write_binary_spc(...)



bool MeasurementInfo::load_from_binary_spc( std::istream &input )
{
/*
 This function was implemented by hand-decoding binary SPC files by wcjohns.
 However, I modified it a little bit when I found
 http://www.ortec-online.com/download/ortec-software-file-structure-manual.pdf
 (which I originally wasnt aware of until 20141107)
 
 //smallint is signed 2-byte integer
 //word is unsigned 32 bit integer
 //singe is a 4 byte float
 
 TSpcHdr = record
 1  1    wINFTYP: smallint; // Must be 1
 2  1    wFILTYP: smallint; // Must be 1 (integer SPC) or 5 (real SPC)
 4  2    wSkip1: array[1..2] of word;
 5  1    wACQIRP: smallint; // Acquisition information record pointer
 6  1    wSAMDRP: smallint; // Sample description record pointer
 7  1    wDETDRP: smallint; // Detector description record pointer
 16 9    wSKIP2: array[1..9] of word;
 17 1    wCALDES: smallint; // Calibration description record pointer
 18 1    wCALRP1: smallint; // First calibration data record pointer
 19 1    wCALRP2: smallint; // Second calibration data record pointer
 20 1    wEFFPRP: smallint; // Efficiency pairs record pointer (first record)
 21 1    wROIRP1: smallint; // Record number of the first of two ROI recs
 22 1    wEPRP: smallint; // Energy pairs record pointer
 23 1    wEPN: smallint; // Number of energy pairs records
 29 6    wSkip3: array[1..6] of word;
 30 1    wEFFPNM: smallint; // Number of efficiency pairs records
 31 1    wSPCTRP: smallint; // pointer to spectrum
 32 1    wSPCRCN: smallint; // number of records in sp
 33 1    wSPCCHN: word; // Number of Channels;
 34 1    wABSTCHN: smallint; // Physical start channel for data
 36 2    sACQTIM: single; // Date and time acquisition start in DECDAY format
 40 4    dACQTI8: double; // Date and time as double precision DECDAY
 44 4    wSkip4: array[1..4] of word;
 45 1    wCHNSRT: smallint; // Start channel number
 47 2    sRLTMDT: single; // Real time in seconds
 49 2    sLVTMDT: single; // Live time in seconds
 50 1    wSkip50: word;
 51 1           I*2     Pointer to FRAM records
 52 1           I*2     Pointer to TRIFID records
 53 1           I*2     Pointer to NaI records
 54 1           I*2     Pointer to Location records
 55 1           I*2     Number of channels of MCS data appended to the histogram data in this file
 56 1           I*2     Pointer to expansion header record (i.e. second header record)
 57-62                 Reserved (must be 0)
 63-64 RRSFCT  R*4     Total random summing factor
 end;
 This data is referenced from the second header record. You first must determine
 whether you have a second header record by reading word 56 of the first header record.
 If this word is non-zero then it points to the second header record which is just
 an extension of the first header record. The current File structures manual has
 an error starting with the description of word 53 of Record 1 in section 4.11.1.
 The current and correct definition for the last few words in record 1 is:
 
 WORD  NAME    Type   Description
 51            I*2     Pointer to FRAM records
 52            I*2     Pointer to TRIFID records
 53            I*2     Pointer to NaI records
 54            I*2     Pointer to Location records
 55            I*2     Number of channels of MCS data appended to the histogram data in this file
 56            I*2     Pointer to expansion header record (i.e. second header record)
 57-62                 Reserved (must be 0)
 63-64 RRSFCT  R*4     Total random summing factor
 
 The following is a description of the expansion header (second header) record
 
 WORD  Type   Description
 1           I*2      Record Identifier (should be decimal 111 for second header record)
 2           I*2      Pointer to GPS data record
 3           I*2      Pointer to 1st record of App specific report text. 1st 16-bits
 are the integer number of text bytes and the 2nd 16 bits are
 the integer report source code (bit 0=Detective-EX). Immediately
 following are the specified number of text bytes. This report
 may span several 128 byte records.
 4-64                Reserved (must be 0)
 
 The App specific report text record(s) are where the Detective-EX stores the
 identification report. The report starts with the fourth byte in the first record
 and continues for as many bytes as specified in the first 16 bit integer in the
 first record. The text is formatted with carriage returns and line feeds so that
 the block could be sent directly to a printer. The last line contains the neutron
 count rate from Detective units that have been updated to add that information
 to their report.
 
 The following is a description of the new entries in the hardware parameters records.
 
 To read the hardware records you must first read the first record to get the total
 number of records (Word 64 of the first hardware parameters record). Next you read
 records 2 through 4 assuming the total number of records is greater than or equal 4.
 If the total number of records is less than 4 then read only the specified number of records.
 
 Next read the 5 Multi-nuclide MDA Preset records if the total number of records
 is greater than or equal 9. The Multi-nuclide MDA Preset records are always complete
 so stop reading records after the first 4 if the total number of records is less
 than 9 (in other words, don't try to read a partial set of MDA Preset records).
 
 Finally, read the Monitors records if the total number of records is greater
 than 9. All remaining bytes in the record block are either stored monitor values
 or padding to round up to the next record. Monitors are stored as pairs of ASCII
 encoded, null terminated strings. The first string in a pair is the monitor label
 and the second string is the monitor value. All bytes after the last string pair
 should be set to zero to indicate zero length string pairs. This prevents garbage
 strings from being read from the file.
 
 The new stuff starts with hardware parameters record 4. The following describes
 parameters record 4:
 
 WORD   Type   Description
 1-4     bit    Validity flags
 5       I*2    Start delay in seconds (DART)
 6       I*2    Conserve delay in seconds (DART)
 7       I*2    Off delay in seconds (DART)
 8       I*2    Power mode (zero=On, nonzero=Conserve)
 9-10    I*4    Hardware status flags from SHOW_STATUS command
 11      I*2    LFR Enabled (zero=disabled, nonzero=enabled)
 12      I*2    ETP Enabled (zero=disabled, nonzero=enabled)
 13-14   R*4    ETP mode protection time in microseconds
 15-16   R*4    Actual HV in volts
 17      I*2    Resolution Enhancer Enabled (zero=disabled, nonzero=enabled)
 18-19   R*4    Manual PZ Setting in microseconds
 20-32          Unused (must be zero)
 33-64   C*16   Data view name strings (4 x 16-character strings, padded with spaces)
 
 If you wanted to include the GV analysis output, then you could read the UFO file.
 It has the same name as the spectrum file with the extension of UFO. The format is
 in the file manual. It isnít complicated, but is another table-driven content file.
 */
  
  reset();
  
  if( !input.good() )
    return false;
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  const istream::pos_type orig_pos = input.tellg();

  try
  {
    input.seekg( 0, ios::end );
    const istream::pos_type eof_pos = input.tellg();
    input.seekg( orig_pos, ios::beg );

    const size_t size = static_cast<size_t>( 0 + eof_pos - orig_pos );
    
    bool foundNeutronDet = false;
    string latitudeStr, longitudeStr;
    
    //As painful as this is, I will read in each variable bellow individually
    //  rather than using a struct, because packing may be an issue
    //  see http://stackoverflow.com/questions/7789668/why-would-the-size-of-a-packed-structure-be-different-on-linux-and-windows-when
    int16_t wINFTYP; // Must be 1
    int16_t wFILTYP; // Must be 1 (integer SPC) or 5 (real SPC)
    int16_t wSkip1[2];
    int16_t wACQIRP; // Acquisition information record pointer
    int16_t wSAMDRP; // Sample description record pointer
    int16_t wDETDRP; // Detector description record pointer
    int16_t wSKIP2[9];
    int16_t wCALDES; // Calibration description record pointer
    int16_t wCALRP1; // First calibration data record pointer
    int16_t wCALRP2; // Second calibration data record pointer
    int16_t wEFFPRP; // Efficiency pairs record pointer (first record)
    int16_t wROIRP1; // Record number of the first of two ROI recs
    int16_t wEPRP; // Energy pairs record pointer
    int16_t wEPN; // Number of energy pairs records
    int16_t wSkip3[6];
    int16_t wEFFPNM; // Number of efficiency pairs records
    int16_t wSPCTRP; // pointer to spectrum
    int16_t wSPCRCN; // number of records in sp
    uint16_t n_channel; // Number of Channels;
    int16_t wABSTCHN; // Physical start channel for data
    float sACQTIM; // Date and time acquisition start in DECDAY format
    double dACQTI8; // Date and time as double precision DECDAY
    int16_t wSkip4[4];
    int16_t wCHNSRT; // Start channel number
    float sRLTMDT; // Real time in seconds
    float sLVTMDT; // Live time in seconds
    int16_t wSkip50;
    int16_t framRecords;//     Pointer to FRAM records
    int16_t TRIFID; // I*2     Pointer to TRIFID records
    int16_t NaI; //I*2     Pointer to NaI records
    int16_t Location;// I*2     Pointer to Location records
    int16_t MCSdata; //I*2     Number of channels of MCS data appended to the histogram data in this file
    int16_t expansionHeader;//     Pointer to expansion header record (i.e. second header record)
    int16_t reserved[5]; // 57-62                 Reserved (must be 0)
    int16_t RRSFCT[2];//  R*4     Total random summing factor

    readBinaryData( input, wINFTYP ); // Must be 1
    if( wINFTYP != 1 )
      throw runtime_error( "First byte indicates not a binary SPC file" );
    
    readBinaryData( input, wFILTYP ); // Must be 1 (integer SPC) or 5 (real SPC)
    if( wFILTYP != 1 && wFILTYP != 5 )
      throw runtime_error( "Second byte indicates not a binary SPC file" );
    
    input.read( (char *)&wSkip1[0], 2*2 );
    readBinaryData( input, wACQIRP ); // Acquisition information record pointer
    readBinaryData( input, wSAMDRP ); // Sample description record pointer
    readBinaryData( input, wDETDRP ); // Detector description record pointer
    input.read( (char *)&wSKIP2[0], 2*9 );
    readBinaryData( input, wCALDES ); // Calibration description record pointer
    readBinaryData( input, wCALRP1 ); // First calibration data record pointer
    readBinaryData( input, wCALRP2 ); // Second calibration data record pointer
    readBinaryData( input, wEFFPRP ); // Efficiency pairs record pointer (first record)
    readBinaryData( input, wROIRP1 ); // Record number of the first of two ROI recs
    readBinaryData( input, wEPRP ); // Energy pairs record pointer
    readBinaryData( input, wEPN ); // Number of energy pairs records
    input.read( (char *)&wSkip3[0], 2*6 );
    readBinaryData( input, wEFFPNM ); // Number of efficiency pairs records
    readBinaryData( input, wSPCTRP ); // pointer to spectrum
    readBinaryData( input, wSPCRCN ); // number of records in sp
    readBinaryData( input, n_channel ); // Number of Channels;
    
    bool isPowerOfTwo = ((n_channel != 0) && !(n_channel & (n_channel - 1)));
    
    //check if it is one channel less than a power of two - Simon Labov emailed
    //  wcjohns an example file like this on 20131216
    if( !isPowerOfTwo && n_channel>=255 )
      isPowerOfTwo = !((n_channel+1) & n_channel);
    
    if( !isPowerOfTwo )
      throw runtime_error( "loadBinarySpcFile(...) - unexpected number of data channels" );
    
    readBinaryData( input, wABSTCHN ); // Physical start channel for data
    readBinaryData( input, sACQTIM ); // Date and time acquisition start in DECDAY format
    readBinaryData( input, dACQTI8 ); // Date and time as double precision DECDAY
    input.read( (char *)&wSkip4[0], 2*4 );
    readBinaryData( input, wCHNSRT ); // Start channel number
    readBinaryData( input, sRLTMDT ); // Real time in seconds
    readBinaryData( input, sLVTMDT ); // Live time in seconds
    readBinaryData( input, wSkip50 );
    readBinaryData( input, framRecords );//     Pointer to FRAM records
    readBinaryData( input, TRIFID ); // I*2     Pointer to TRIFID records
    readBinaryData( input, NaI ); //I*2     Pointer to NaI records
    readBinaryData( input, Location );// I*2     Pointer to Location records
    readBinaryData( input, MCSdata ); //I*2     Number of channels of MCS data appended to the histogram data in this file
    readBinaryData( input, expansionHeader );//     Pointer to expansion header record (i.e. second header record)
    input.read( (char *)&reserved[0], 2*5 );
    input.read( (char *)&RRSFCT[0], 2*2 );
    
    if( !input.good() )
      throw runtime_error( "Error reading header data" );
  
/*
    cout << "wINFTYP=" <<wINFTYP << endl;
    cout << "wINFTYP=" <<wINFTYP << endl;
    cout << "wACQIRP=" <<wACQIRP << endl;
    cout << "wSAMDRP=" << wSAMDRP<< endl;
    cout << "wDETDRP=" <<wDETDRP << endl;
    cout << "wCALDES=" << wCALDES<< endl;
    cout << "wCALRP1=" << wCALRP1<< endl;
    cout << "wCALRP2=" <<wCALRP2 << endl;
    cout << "wEFFPRP=" << wEFFPRP<< endl;
    cout << "wROIRP1=" << wROIRP1<< endl;
    cout << "wEPRP=" << wEPRP<< endl;
    cout << "wEPN=" <<wEPN << endl;
    cout << "wEFFPNM=" <<wEFFPNM << endl;
    cout << "wSPCTRP=" <<wSPCTRP << endl;
    cout << "wSPCRCN=" <<wSPCRCN << endl;
    cout << "n_channel=" << n_channel<< endl;
    cout << "wABSTCHN=" << wABSTCHN<< endl;
    cout << "sACQTIM=" << sACQTIM<< endl;
    cout << "dACQTI8=" << dACQTI8<< endl;
    cout << "wCHNSRT=" << wCHNSRT<< endl;
    cout << "sRLTMDT=" <<sRLTMDT << endl;
    cout << "sLVTMDT=" << sLVTMDT<< endl;
    cout << "wSkip50=" << wSkip50<< endl;
    cout << "framRecords=" <<framRecords << endl;
    cout << "TRIFID=" << TRIFID<< endl;
    cout << "NaI=" << NaI<< endl;
    cout << "Location=" << Location<< endl;
    cout << "MCSdata=" <<MCSdata << endl;
    cout << "expansionHeader=" << expansionHeader<< endl;
    cout << "RRSFCT[0]=" << RRSFCT[0] << endl;
    cout << "RRSFCT[1]=" << RRSFCT[1] << endl;
*/
    
    //calibration data
    vector<float> calib_coefs;
    if( wCALRP1 > 0 )
    {
      input.seekg( (wCALRP1-1)*128 + orig_pos, ios::beg );
      int16_t wAFIT; // 2 Above knee efficiency calibration type
      int16_t wBFIT; // 4 Below knee efficiency calibration type
      int16_t wEFFPRS; // 6 Number of efficiency pairs
      int16_t wNCH; // 8 number of channels in spectrum
      float sKNEE; // 12 Detector knee (keV)
      float sASIG; // 16 2-sigma uncertainty above knee
      float sBSIG; // 20 2-sigma uncertainty below knee
      float sEC1; // 24 Energy vs. channel coefficient A
      float sEC2; // 28 Energy vs. channel coefficient B
      float sEC3; // 32 Energy vs. channel coefficient C
      float sFC1; //FWHM vs. channel coefficient A (actually in most file)
      float sFC2; //FWHM vs. channel coefficient B (actually in most file)
      float sFC3; //FWHM vs. channel coefficient C (actually in most file)
      float sPE1; //Above knee eff. vs. energy coeff A or poly coeff 1
      float sPE2; //Above knee eff. vs. energy coeff B or poly coeff 2
      float sPE3; //Above knee eff. vs. energy coeff C or poly coeff 3
      float sSE1; //Below knee eff. vs. energy coeff A or poly coeff 4
      float sSE2; //Below knee eff. vs. energy coeff B or poly coeff 5
      float sSE3; //Below knee eff. vs. energy coeff C or poly coeff 6
      int16_t wFWHTYP; // FWHM type
      int16_t wRES1; // reserved
      int16_t wRES2; // reserved
      int16_t wENGPRS; // Number of energy pairs
      int16_t wDETNUM; //Detector Number
      int16_t wNBKNEE; // Number of calibration points below knee
      float sENA2; // Temp energy calibration
      float sENB2; // Temp energy calibration
      float sENC2; // Temp energy calibration
      float sCALUNC; //Calibration source uncertainty
      float sCALDIF; // Energy calibration difference
      float sR7; // Polynomial coefficient 7
      float sR8; // Polynomial coefficient 8
      float sR9; // Polynomial coefficient 9
      float sR10; // Polynomial coefficient 10
    
      
      readBinaryData( input, wAFIT ); // 2 Above knee efficiency calibration type
      readBinaryData( input, wBFIT ); // 4 Below knee efficiency calibration type
      readBinaryData( input, wEFFPRS ); // 6 Number of efficiency pairs
      readBinaryData( input, wNCH ); // 8 number of channels in spectrum
      readBinaryData( input, sKNEE ); // 12 Detector knee (keV)
      readBinaryData( input, sASIG ); // 16 2-sigma uncertainty above knee
      readBinaryData( input, sBSIG ); // 20 2-sigma uncertainty below knee
      readBinaryData( input, sEC1 ); // 24 Energy vs. channel coefficient A
      readBinaryData( input, sEC2 ); // 28 Energy vs. channel coefficient B
      readBinaryData( input, sEC3 ); // 32 Energy vs. channel coefficient C
      readBinaryData( input, sFC1 ); //FWHM vs. channel coefficient A (actually in most file)
      readBinaryData( input, sFC2 ); //FWHM vs. channel coefficient B (actually in most file)
      readBinaryData( input, sFC3 ); //FWHM vs. channel coefficient C (actually in most file)
      readBinaryData( input, sPE1 ); //Above knee eff. vs. energy coeff A or poly coeff 1
      readBinaryData( input, sPE2 ); //Above knee eff. vs. energy coeff B or poly coeff 2
      readBinaryData( input, sPE3 ); //Above knee eff. vs. energy coeff C or poly coeff 3
      readBinaryData( input, sSE1 ); //Below knee eff. vs. energy coeff A or poly coeff 4
      readBinaryData( input, sSE2 ); //Below knee eff. vs. energy coeff B or poly coeff 5
      readBinaryData( input, sSE3 ); //Below knee eff. vs. energy coeff C or poly coeff 6
      readBinaryData( input, wFWHTYP ); // FWHM type
      readBinaryData( input, wRES1 ); // reserved
      readBinaryData( input, wRES2 ); // reserved
      readBinaryData( input, wENGPRS ); // Number of energy pairs
      readBinaryData( input, wDETNUM ); //Detector Number
      readBinaryData( input, wNBKNEE ); // Number of calibration points below knee
      readBinaryData( input, sENA2 ); // Temp energy calibration
      readBinaryData( input, sENB2 ); // Temp energy calibration
      readBinaryData( input, sENC2 ); // Temp energy calibration
      readBinaryData( input, sCALUNC ); //Calibration source uncertainty
      readBinaryData( input, sCALDIF ); // Energy calibration difference
      readBinaryData( input, sR7 ); // Polynomial coefficient 7
      readBinaryData( input, sR8 ); // Polynomial coefficient 8
      readBinaryData( input, sR9 ); // Polynomial coefficient 9
      readBinaryData( input, sR10 ); // Polynomial coefficient 10
    
      calib_coefs.push_back( sEC1 );
      calib_coefs.push_back( sEC2 );
      calib_coefs.push_back( sEC3 );

/*
      cout << "wAFIT=" <<wAFIT << endl;
      cout << "wBFIT=" <<wBFIT << endl;
      cout << "wEFFPRS=" <<wEFFPRS << endl;
      cout << "wNCH=" << wNCH<< endl;
      cout << "sKNEE=" <<sKNEE << endl;
      cout << "sASIG=" <<sASIG << endl;
      cout << "sBSIG=" <<sBSIG << endl;
      cout << "sEC1=" <<sEC1 << endl;
      cout << "sEC2=" <<sEC2 << endl;
      cout << "sEC3=" <<sEC3 << endl;
      cout << "sFC1=" <<sFC1 << endl;
      cout << "sFC2=" <<sFC2 << endl;
      cout << "sFC3=" <<sFC3 << endl;
      cout << "sPE1=" <<sPE1 << endl;
      cout << "sPE2=" <<sPE2 << endl;
      cout << "sPE3=" <<sPE3 << endl;
      cout << "sSE1=" <<sSE1 << endl;
      cout << "sSE2=" <<sSE2 << endl;
      cout << "sSE3=" <<sSE3 << endl;
      cout << "wFWHTYP=" <<wFWHTYP << endl;
      cout << "wRES1=" <<wRES1 << endl;
      cout << "wRES2=" <<wRES2 << endl;
      cout << "wENGPRS=" <<wENGPRS << endl;
      cout << "wDETNUM=" <<wDETNUM << endl;
      cout << "wNBKNEE=" <<wNBKNEE << endl;
      cout << "sENA2=" <<sENA2 << endl;
      cout << "sENB2=" <<sENB2 << endl;
      cout << "sENC2=" <<sENC2 << endl;
      cout << "sCALUNC=" <<sCALUNC << endl;
      cout << "sCALDIF=" <<sCALDIF << endl;
      cout << "sR7=" <<sR7 << endl;
      cout << "sR8=" <<sR8 << endl;
      cout << "sR9=" <<sR9 << endl;
      cout << "sR10=" << sR10<< endl;
*/
    }//if( wCALRP1 > 0 )
    
    if( !input.good() )
      throw runtime_error( "Error reading calibration data" );
    
    string instrument_id;
    double sum_gamma = 0.0;
    double total_neutrons = 0.0;
    float total_neutron_count_time = 0.0;
    std::shared_ptr<DetectorAnalysis> analysis;
    std::shared_ptr< vector<float> > channel_data( new vector<float>( n_channel ) );
    
  
    boost::posix_time::ptime meas_time;
    string title;
    string manufacturer = "Ortec";
    string inst_model = "Detective";
    string type_instrument = "RadionuclideIdentifier";
    DetectorType type_detector = kUnknownDetector;
      
    
    {//begin codeblock to get acquisition information
      char namedata[17], datedata[19];
      datedata[18] = namedata[16] = '\0';
      datedata[9] = ' ';
      
      
      //Not working correctly - f-it for right now
      if( wACQIRP > 0 )
      {
        input.seekg( 128*(wACQIRP-1) + orig_pos, ios::beg );
        input.read( namedata, 16 );
        input.read( datedata, 9 );
        input.read( datedata+10, 3 );  //just burning off 3 bytes
        input.read( datedata+10, 8 );
      
        if( !input.good() )
          throw runtime_error( "Didnt successfully read date data" );
      
        string name( namedata, namedata+16 );
        trim( name );
  
      
        //name seems to always be 'DetectiveEX.SPC'
        if( istarts_with( name, "Detective" ) )
        {
          type_instrument = "Radionuclide Identifier";
          manufacturer = "Ortec";
          type_detector = kDetectiveDetector;
        }//if( istarts_with( name, "DetectiveEX" ) )
      
        try
        {
          meas_time = time_from_string( datedata );
          //cout << "meas_time=" << UtilityFunctions::to_iso_string( meas_time ) << endl;
        }catch(...)
        {
          cerr << "MeasurementInfo::loadBinarySpcFile(...): invalid date string: "
               << datedata << endl;
        }
      
/*
        cout << "name=" << name << endl;
        cout << "datedata='" << datedata << "'" << endl;
*/
      }//if( wACQIRP > 0 )
      
      if( wACQIRP > 0 )
      {
        input.seekg( 128*(wACQIRP-1) + orig_pos + 90, ios::beg );
        char start_date_of_sample_collection[11] = { '\0' };
        char start_time_of_sample_collection[9] = { '\0' };
        char stop_date_of_sample_collection[11] = { '\0' };
        char stop_time_of_sample_collection[9] = { '\0' };
        input.read( start_date_of_sample_collection, 10 );
        char dummy[2];
        input.read( dummy, 2 );
        input.read( start_time_of_sample_collection, 8 );
        input.read( stop_date_of_sample_collection, 10 );
        input.read( stop_time_of_sample_collection, 8 );
      
/*
        cout << "start_date_of_sample_collection='" << start_date_of_sample_collection << "'" << endl;
        cout << "start_time_of_sample_collection='" << start_time_of_sample_collection << "'" << endl;
        cout << "stop_date_of_sample_collection='" << stop_date_of_sample_collection << "'" << endl;
        cout << "stop_time_of_sample_collection='" << stop_time_of_sample_collection << "'" << endl;
*/
      }//if( wACQIRP > 0 )
    }//end codeblock to get acquisition information
  
    //TODO: look for the following strings and include in results
    //State Of Health^@OK
    //Gamma Dose Rate^@0.07 ?Sv/h
    //  [Ge]^@Detector Temperature^@OK
    //  Battery Voltage^@15.37 Volts
    //  Battery Time Remaining^@125 Min
    //  Cooler Body Temperature^@OK
    //  Cooler Drive Voltage^@OK
    //  Cold-Tip Temperature^@OK
    //  HV Bias^@-3509 V
    
    if( expansionHeader )
    {
      input.seekg( 128*(expansionHeader-1) + orig_pos, ios::beg );
      
      if( !input.good() )
      {
        stringstream msg;
        msg << "Unable to read expansion header in file, possible pointer "
            << expansionHeader << " (location " << 128*(expansionHeader-1)
            << " of size=" << size << ")" << endl;
        throw runtime_error( msg.str() );
      }//if( !input.good() )
      
      int16_t recordID, gpsPointer, firstReportPtr;
      readBinaryData( input, recordID );
      readBinaryData( input, gpsPointer );
      readBinaryData( input, firstReportPtr );
      
      if( recordID != 111 )
      {
        gpsPointer = firstReportPtr = 0;
        cerr << "Binary SPC file has invalid expansion header" << endl;
      }//if( recordID != 111 )
      
      if( gpsPointer )
      {
        //See refD3BAVOI7JG
//        input.seekg( 128*(gpsPointer-1) + orig_pos, ios::beg );
        
//        if( input.good() )
//        {
//          uint16_t ntxtbytes;
//          readBinaryData( input, ntxtbytes );
//        
//          cerr << "ntxtbytes=" << ntxtbytes << endl;
//          ntxtbytes = std::min(ntxtbytes, uint16_t(120));
//          vector<char> data(ntxtbytes+1);
//          data[ntxtbytes] = '\0';
//          input.read( &data[0], ntxtbytes );
//          for( size_t i = 0; i < ntxtbytes; ++i )
//            cout << i << ": " << int(data[i]) << ", '" << data[i] << "'" << endl;
//        }else
//        {
//          cerr << "Failed to be able to read GPS REcord" << endl;
//        }
        cerr << "Binary SPC file has not yet implemented GPS coordinates decoding" << endl;
      }

      
      if( firstReportPtr > 0 )
      {
        input.seekg( 128*(firstReportPtr-1) + orig_pos, ios::beg );
        if( !input.good() )
        {
          stringstream msg;
          msg << "Unable to read report in file, possible bad report pointer "
              << firstReportPtr << " (location " << 128*(firstReportPtr-1)
              << " of size=" << size << ")" << endl;
          throw runtime_error( msg.str() );
        }//if( !input.good() )
        
        uint16_t ntxtbytes, sourcecode;
        readBinaryData( input, ntxtbytes );
        readBinaryData( input, sourcecode );

        if( ntxtbytes > 2048 )
          ntxtbytes = 0;
          
        if( ntxtbytes > 0 )
        {
          vector<char> data(ntxtbytes+1);
          data[ntxtbytes] = '\0';
          input.read( &data[0], ntxtbytes );

          string term;
          
          {//being codeblock to look for neutrons
            //Apparently capitalization isnt consistent, so will convert to
            //  lowercase; I didnt test yet for longitude/latitude, or nuclide
            //  IDs, so I dont want to convert to lower case for those yet.
            string datastr( data.begin(), data.end() );
            UtilityFunctions::to_lower( datastr );
            
            term = "total neutron counts = ";
            string::const_iterator positer, enditer;
            positer = std::search( datastr.begin(), datastr.end(),
                                   term.begin(), term.end() );
            if( positer != datastr.end() )
            {
              term = "total neutron counts = ";
              positer = std::search( datastr.begin(), datastr.end(),
                                    term.begin(), term.end() );
              foundNeutronDet = true;
              total_neutrons = atof( &(*(positer+term.size())) );
            }else
            {
              term = "neutron counts";
              positer = std::search( datastr.begin(), datastr.end(),
                                    term.begin(), term.end() );
              if( positer != datastr.end() )
              {
                foundNeutronDet = true;
//                positer += 17;
                total_neutrons = atof( &(*(positer+term.size())) ); //atof( &(*positer) );
              }
            }
            
          
            term = "total neutron count time = ";
            positer = std::search( datastr.begin(), datastr.end(),
                                  term.begin(), term.end() );
            if( positer != datastr.end() )
            {
              foundNeutronDet = true;
              total_neutron_count_time = static_cast<float>( atof( &(*positer) + 27 ) );
            }
          }//end codeblock to look for neutrons
          
          
          
          //Other strings to look for:
          //"Found Nuclides"
          //"Suspect Nuclides"
          //"Top Lines"
          //"GPS"  or something like "GPS Location not determined"
          //"Gamma Dose Rate"
          //"Version"
          //"ID Report"
          //"Total neutron count time = "
          //"Average neutron count rate = "
          //
          //Maybye look for "ID Report"
        
          vector<char>::iterator positer, enditer;
          term = "Latitude";
          positer = std::search( data.begin(), data.end(),
                                term.begin(), term.end() );
          if( positer != data.end() )
          {
            positer += 8;
            while( (positer!=data.end()) && !isdigit(*positer) )
              ++positer;
            enditer = positer;
            while( (enditer!=data.end()) && (*enditer)!='\n' )
              ++enditer;
            latitudeStr.insert(latitudeStr.end(), positer, enditer );
            latitudeStr.erase( std::remove_if(latitudeStr.begin(), latitudeStr.end(), not_alpha_numeric), latitudeStr.end());
          }//if( there is latitude info )
          
          term = "Longitude";
          positer = std::search( data.begin(), data.end(),
                                term.begin(), term.end() );
          if( positer != data.end() )
          {
            positer += 9;
            while( (positer!=data.end()) && !isdigit(*positer) )
              ++positer;
            enditer = positer;
            while( (enditer!=data.end()) && (*enditer)!='\n' )
              ++enditer;
            longitudeStr.insert(longitudeStr.end(), positer, enditer );
            longitudeStr.erase(std::remove_if(longitudeStr.begin(), longitudeStr.end(), not_alpha_numeric), longitudeStr.end());
          }//if( there is longitude info )
        
          
          string found_term = "Found Nuclides";
          vector<char>::iterator nucpos = std::search( data.begin(), data.end(),
                                                       found_term.begin(), found_term.end() );
          if( nucpos == data.end() )
          {
            found_term = "Found:";
            nucpos = std::search( data.begin(), data.end(), found_term.begin(), found_term.end() );
            
            //Found: actually looks a little different, ex:
            //  Found: Industrial(2)  NORM(1) Other(2)^M
            //  Co60    Co57    Mn54    K40     Co56
            //XXX TODO: Which right now were putting Industrial/NORM/OTHER as
            //          their own analysis results, in addition to the actual
            //          nuclides. Should fix.
          }
          
          if( nucpos != data.end() )
          {
            if( !analysis )
              analysis = std::make_shared<DetectorAnalysis>();
            
            //Should reformat this list seperated by newlines to csv or somehting
            //  and also I dont know if "Suspect Nuclides" is garunteed to be there
            string suspect_term = "Suspect Nuclides";
            vector<char>::iterator suspectpos = std::search( nucpos, data.end(),
                                                        suspect_term.begin(), suspect_term.end() );
            if( suspectpos == data.end() )
            {
              suspect_term = "Suspect:";
              suspectpos = std::search( nucpos, data.end(), suspect_term.begin(), suspect_term.end() );
            }
            
            
            string found_nucs_str( nucpos+found_term.size(), suspectpos );
            vector<string> found_nucs, suspect_nucs;
            split( found_nucs, found_nucs_str, "\t,\n\r\0");
              
            for( string &nuc : found_nucs )
            {
              nuc.erase(std::remove_if(nuc.begin(), nuc.end(), not_alpha_numeric), nuc.end());
              trim( nuc );
              ireplace_all( nuc, "  ", " " );
              
              if( icontains( nuc, "keep counting" ) )
              {
                remarks_.push_back( nuc );
              }else if( nuc.size() )
              {
                DetectorAnalysisResult result;
                result.remark_ = "Found";
                //result.id_confidence_ = "Found";
                result.nuclide_ = nuc;
                analysis->results_.push_back( result );
              }
            }//for( string &nuc : found_nucs )
            
            const string lines_term = "Top Lines";
            vector<char>::iterator linesiter = std::search( suspectpos, data.end(),
                                                           lines_term.begin(), lines_term.end() );
            
            if( suspectpos != data.end() )
            {
              string suspect_nucs_str = string( suspectpos+suspect_term.size(), linesiter );
              string::size_type endpos = suspect_nucs_str.find_first_of("\0");
              if( endpos != string::npos )
                suspect_nucs_str.substr(0,endpos);
              
              split( suspect_nucs, suspect_nucs_str, "\t,\n\r\0" );
              for( string &nuc : suspect_nucs )
              {
                nuc.erase( std::remove_if(nuc.begin(), nuc.end(), not_alpha_numeric), nuc.end());
                trim( nuc );
                ireplace_all( nuc, "  ", " " );
              
                if( icontains( nuc, "keep counting" ) )
                {
                  remarks_.push_back( nuc );
                }else if( nuc.size() )
                {
                  DetectorAnalysisResult result;
                  result.remark_ = "Suspect";
                  //result.id_confidence_ = "Suspect";
                  result.nuclide_ = nuc;
                  analysis->results_.push_back( result );
                }
              }//for( string &nuc : found_nucs )
            }//if( suspectpos != data.end() )
            
//               const char *dose_rate = strstr( &data[0], "Gamma Dose Rate" );
//               if( dose_rate )
//               cerr << "Found dose rate: " << dose_rate << endl << endl;
              
            if( (data.end() - (linesiter+lines_term.size())) > 0 )
            {
              string toplines( linesiter+10, data.end() );
              vector<string> lines;
              split( lines, toplines, "\r\n" );
              
              for( size_t i = 0; i < lines.size(); ++i )
              {
//                 lines[i].erase( std::remove_if(lines[i].begin(), lines[i].end(), not_alpha_numeric), lines[i].end());
                trim( lines[i] );
                ireplace_all( lines[i], "  ", " " );
                ireplace_all( lines[i], "\t", "&#009;" );  //replace tab characters with tab character code
                
                if( lines[i].empty()
                     || istarts_with(lines[i], "Longitude")
                     || istarts_with(lines[i], "GPS") )
                  break;

                char buffer[256];
                snprintf( buffer, sizeof(buffer), "Top Line %i: %s", int(i), lines[i].c_str() );
                analysis->remarks_.push_back( buffer );
              }//for( size_t i = 0; i < lines.size(); ++i )
            }//if( there are top lines )
          }//if( nucpos != data.end() )
        }//if( ntxtbytes )
        
      }//if( firstReportPtr )
    }//if( expansionHeader )
      

    if( wSAMDRP > 0 )
    {
      input.seekg( 128*(wSAMDRP-1) + orig_pos, ios::beg );
      vector<char> data(128);
      input.read( &data[0], 128 );
      
      data.erase(std::remove_if(data.begin(), data.end(), not_alpha_numeric), data.end());
      
      string remark( data.begin(), data.end() );
      trim( remark );
      if( remark.size() )
        remarks_.push_back( "Sample Description: " + remark );
    }//if( wSAMDRP )
      
    if( wDETDRP > 0 && input.seekg(128*(wDETDRP-1) + orig_pos, ios::beg) )
    {
      input.seekg( 128*(wDETDRP-1) + orig_pos, ios::beg );
      vector<char> data(128);
      input.read( &data[0], 128 );
      data.erase(std::remove_if(data.begin(), data.end(), not_alpha_numeric), data.end());
      
      if( data.size() )
      {
        data.push_back( '\0' );
        instrument_id = &data[0];
        size_t len = instrument_id.find_last_not_of( " \0\t" );
        if( len != string::npos )
          instrument_id = instrument_id.substr( 0, len + 1 );
        else
          instrument_id = "";
      }//if( data.size() )
      
      trim( instrument_id );
      ireplace_all( instrument_id, "\n", " " );
      ireplace_all( instrument_id, "\r", " " );
      ireplace_all( instrument_id, "  ", " " );
      
  
      //Some detective EX100s, like in ref49KB84PGM4, have the serial number in
      //  not in the standard position (but may have the model in the standard
      //  position), so we will attempt to check this out and add the SN.
      vector<int16_t> calibrationpos;
      if( wCALRP1 > 0 )
        calibrationpos.push_back( wCALRP1 );
      if( wCALRP2 > 0 )
        calibrationpos.push_back( wCALRP2 );
      if( wCALDES > 0 )
        calibrationpos.push_back( wCALDES );
      
      for( size_t i = 0; i < calibrationpos.size(); ++i )
      {
        vector<char> data(128);
        
        if( !input.seekg(128*(calibrationpos[i]-1) + orig_pos, ios::beg) )
          continue;
        if( !input.read( &data[0], 128 ) )
          continue;
        
        data.push_back( '\0' );
        string strdata = &data[0];
        size_t pos = strdata.find( "SN:" );
        if( pos != string::npos )
        {
          strdata = strdata.substr( pos+3 );
          pos = strdata.find_last_not_of( " \0\t" );
          if( pos != string::npos )
            strdata = strdata.substr( 0, pos + 1 );
          strdata.erase(std::remove_if(strdata.begin(), strdata.end(), not_alpha_numeric), strdata.end());
          ireplace_all( strdata, "\n", " " );
          ireplace_all( strdata, "\r", " " );
          ireplace_all( strdata, "  ", " " );
          trim( strdata );
          if( strdata.length() && (instrument_id.find(strdata)==string::npos) )
            instrument_id += (instrument_id.size() ? " " : "") + strdata;
        }//if( pos != string::npos )
      }//for( size_t i = 0; i < calibrationpos.size(); ++i )
      
      type_detector = detective_model_from_serial( instrument_id );
      switch( type_detector )
      {
        case kDetectiveExDetector:
          inst_model = foundNeutronDet ? "DetectiveEX" : "DetectiveDX";
          break;
          
        case kDetectiveEx100Detector:
          inst_model = foundNeutronDet ? "DetectiveEX100" : "DetectiveDX100";
          break;
          
        case kMicroDetectiveDetector:
          inst_model = "MicroDetective";
          break;
          
        default:
          type_detector = kDetectiveDetector;
          inst_model = "Detective";
          break;
      }//switch( new_type )
      
    }//if( wDETDRP > 0 && input.seekg(128*(wDETDRP-1) + orig_pos, ios::beg) )
    
//    cout << instrument_id << ": " << sFC1 << ", " << sFC2 << ", " << sFC3 << endl;
  
//      
//    cerr << "instrument_id=" << instrument_id << endl;
    /*
      const char *note = &(data[0])+384;
      const char *end_note = strstr( note, "\r" );
      if( (end_note - note) > 256 )
      end_note = note + 256;
      string remark = string( note, end_note );
      trim( remark );
    */
      
    //read in channel data
    vector<float> &counts_ref = *channel_data;
    input.seekg( 128*(wSPCTRP-1) + orig_pos, ios::beg );
    
    if( !input.good() )
      throw runtime_error( "Unable to read channel data" );
    
    const size_t last_expected = static_cast<size_t>( 4*n_channel + 128*(wSPCTRP-1) + orig_pos );
    if( last_expected > size_t(12+eof_pos) )  //12 is a just in case...
      throw runtime_error( "File not expected size" );
    
    if( wFILTYP == 1 )
    {
      vector<uint32_t> int_channel_data( n_channel );
      input.read( (char *)&int_channel_data[0], 4*n_channel );
        
      for( size_t i = 0; i < n_channel; ++i )
        counts_ref[i] = static_cast<float>( int_channel_data[i] );
    }else //if( wFILTYP == 5 )
    {
      input.read( (char *) &(counts_ref[0]), 4*n_channel );
    }//if( file is integer channel data ) / else float data
      
    counts_ref[0] = 0;
    counts_ref[n_channel-1] = 0;
      
    for( size_t i = 0; i < n_channel; ++i )
      sum_gamma += counts_ref[i];
      
    manufacturer_       = manufacturer;
    instrument_type_    = type_instrument;
    instrument_model_   = inst_model;
    detector_type_      = type_detector;
    detectors_analysis_ = analysis;
    instrument_id_      = instrument_id;
      
    MeasurementShrdPtr meas = std::make_shared<Measurement>();
    
    if( sLVTMDT < 0.0 || sRLTMDT < 0.0 )
      throw runtime_error( "Invalid real or live time" );
    
    meas->live_time_ = sLVTMDT;
    meas->real_time_ = sRLTMDT;
    meas->start_time_ = meas_time;
    meas->gamma_counts_ = channel_data;
    meas->gamma_count_sum_ = sum_gamma;

    meas->calibration_coeffs_ = calib_coefs;
    meas->energy_calibration_model_ = Measurement::Polynomial; //Measurement::FullRangeFraction
    
    //File ref9HTGHJ9SXR has the neutron information in it, but
    //  the serial number claims this is a micro-DX (no neutron detector), and
    //  a found nuclide is neutron on hydrogen, but yet no neutrons are actually
    //  detected.  This implies we should check if the serial number contains
    //  "DX" and no neutrons are detected, then set meas->contained_neutron_
    //  to false.
    //  Not doing this now due to ambiguity.
    meas->contained_neutron_ = foundNeutronDet;
    meas->neutron_counts_sum_ = total_neutrons;
    if( foundNeutronDet || (total_neutrons>0.0) )
      meas->neutron_counts_.push_back( static_cast<float>(total_neutrons) );

    if( total_neutron_count_time > 0.0 )
    {
      char buffer[128];
      snprintf( buffer, sizeof(buffer),
          "Total neutron count time = %f seconds", total_neutron_count_time );
      meas->remarks_.push_back( buffer );
    }
    
    
    if( longitudeStr.size() && latitudeStr.size() )
    {
      meas->latitude_ = ortecLatOrLongStrToFlt(latitudeStr);
      meas->longitude_ = ortecLatOrLongStrToFlt(longitudeStr);
    }//if( longitudeStr.size() && latitudeStr.size() )
    
    measurements_.push_back( meas );

    cleanup_after_load();
  }catch( std::exception &e )
  {
    cerr  << SRC_LOCATION << "caught:\n" << e.what() << endl;
    reset();
    input.clear();
    input.seekg( orig_pos, ios::beg );
    return false;
  }//try / catch

  return true;
}//bool loadBinarySpcFile( std::istream &input )



bool MeasurementInfo::load_binary_exploranium_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;
  
  const bool loaded = load_from_binary_exploranium( file );
  
  if( loaded )
    filename_ = filename;
  
  return loaded;
}//bool load_binary_exploranium_file( const std::string &file_name )





bool MeasurementInfo::load_from_binary_exploranium( std::istream &input )
{
  //Currently doesnt:
  //  -CHSUM	Checksum not checked
  //  -Dose rate not extracted
  //  -Software version not checked/extracted
  //  -Roi/Line results not checked
  //  -Check for library being used
  //  -Spectrum, Survey, CZT, and Dose not checked for or used
  //  -Channel number hasnt been verified
  //  -GR135 v1 and GR130 not verified
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  set<const char *> warnings;
  
  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  input.seekg( 0, ios::end );
  const istream::pos_type eof_pos = input.tellg();
  input.seekg( orig_pos, ios::beg );
  
  const size_t size = static_cast<size_t>( 0 + eof_pos - orig_pos );
  
  if( size < 513 )
    return false;
  
  try
  {
    char charbuff[128];
    vector<char> buffer( size );
    if( !input.read( &buffer[0], size) )
      throw runtime_error( "Failed Read" );
    
    const char *zs = "ZZZZ";
    const char *start = &buffer[0];
    const char *end = start + size;
  
    const string asc = "ASC";
    const string asd = "ASD";
    
    //We have to first identify the byte each record starts at before decoding
    //  any of them, since we need to know the record size while decoding the
    //  record (could probably be clever about avoiding this, but whatever).
    vector<size_t> recordstarts;
    for( size_t i = 0; i < (size-1-2*256);
         i = (std::search(start+i+1, end, zs, zs+4) - start) )
    {
      const char *data = start + i;
      
      if( data[0]!='Z' || data[1]!='Z' || data[2]!='Z' || data[3]!='Z' )
        continue;
      
      //Record type Type: A=Spectrum, S=Survey, C=CZT, D=Dose
      const bool is135v2 = (asc.find(data[8])!=string::npos);
      const bool is135v1 = (asc.find(data[4])!=string::npos);
      const bool is130v0 = (asd.find(data[6])!=string::npos);
    
      const bool dtSpectrum = (is135v2 && (data[8] == 'A'))
                              || (is135v1 && (data[4] == 'A'))
                              || (is130v0 && (data[6] == 'A'));
      const bool dtCzt = (is135v2 && (data[8] == 'C'))
                          || (is135v1 && (data[4] == 'C'));
    
      const size_t record_size = recordstarts.size() ? i - recordstarts.back() : 0;
    
      const size_t num_channels = (is130v0 ? 256 : 1024);
      const bool valid = (is135v2 || is135v1 || is130v0)
                         && (dtSpectrum || dtCzt)
                         && ((record_size > num_channels*2) || recordstarts.empty());
      if( valid )
      {
        recordstarts.push_back( i );
        i += 512;
      }//if( valid )
    }//for( loop over file and identify starts of records )
    
    for( size_t j = 0; j < recordstarts.size(); ++j )
    {
      const char *data = start + recordstarts[j];
      
      const bool is135v2 = (asc.find(data[8])!=string::npos);
      const bool is135v1 = (asc.find(data[4])!=string::npos);
      const bool is130v0 = (asd.find(data[6])!=string::npos);
      
      const bool dtCzt = (is135v2 && (data[8] == 'C'))
                         || (is135v1 && (data[4] == 'C'));
      
      const bool found1350 = (data[4] == '1') && (data[5] == '3') && (data[6] == '5');
      const bool found1024 = (data[5] == 0x0) && (data[6] == 0x4);
        
      if( is135v2 && !found1350 )
        warnings.insert( "The header is missing the string \"135\" at offsets 5, 6, and 7" );
        
      if( is135v1 && !found1024 )
        warnings.insert( "The header is missing the 16-bit integer \"1024\" at offsets 6 and 7'" );

      const size_t record_size = ( (j==(recordstarts.size()-1))
                                  ? (size - recordstarts[j])
                                  : (recordstarts[j+1] - recordstarts[j]));
      
      std::shared_ptr<Measurement> meas( new Measurement() );
      
      //try to get the date/time
      const char *datepos = data + 7;
      if( is130v0 )
        datepos = data + 7;
      else if( is135v1 )
      {
        if( record_size == 2099 )
          datepos = data + 7;
        else if( record_size == 2124 )
          datepos = data + 13;
        else if( record_size == 2127 )
          datepos = data + 13;
      }else //if( is135v2 )
        datepos = data + 13;
      
      uint8_t timeinfo[6]; //= {year, month, day, hour, minutes, seconds}
      for( size_t i = 0; i < 6; ++i )
        timeinfo[i] = 10*(((datepos[i]) & 0xF0) >> 4) + ((datepos[i]) & 0xF);
      
      try
      {
        const boost::gregorian::date the_date( 2000 + timeinfo[0], timeinfo[1], timeinfo[2]);
        const boost::posix_time::time_duration the_time( timeinfo[3], timeinfo[4], timeinfo[5]);
        meas->start_time_ = boost::posix_time::ptime( the_date, the_time );
        
        //if( is130v0 ){
        //  float battery_voltage = (10*(((datepos[7]) & 0xF0) >> 4) + ((datepos[7]) & 0xF)) / 10.0f;
        //}
      }catch(...)
      {
      }
      
      if( meas->start_time_.is_special() )
      {
        //We're desperate here - lets search the entire header area.
        for( size_t i = 0; i < 75; ++i )
        {
          datepos = data + 4 + i;
          for( size_t i = 0; i < 6; ++i )
            timeinfo[i] = 10*(((datepos[i]) & 0xF0) >> 4) + ((datepos[i]) & 0xF);
          try
          {
            const boost::gregorian::date the_date( 2000 + timeinfo[0], timeinfo[1], timeinfo[2]);
            const boost::posix_time::time_duration the_time( timeinfo[3], timeinfo[4], timeinfo[5]);
            meas->start_time_ = boost::posix_time::ptime( the_date, the_time );
            break;
          }catch(...)
          {
          }
        }//for( size_t i = 0; i < 75; ++i )
      }//if( meas->start_time_.is_special() )
      
      //neutrn info
      meas->contained_neutron_ = is135v2;
      if( is135v2 )
      {
        uint16_t neutsum;
        memcpy( &neutsum, data+36, 2 );
        meas->neutron_counts_sum_ = neutsum;
        meas->neutron_counts_.resize( 1, static_cast<float>(neutsum) );
      }//if( is135v2 )
      
      //Serial number
      uint16_t serialnum = 0, version = 0;
      if( is135v2 )
      {
        memcpy( &serialnum, data+40, 2 );
        memcpy( &version, data+42, 2 );
      }else if( is130v0 )
      {
        memcpy( &serialnum, data+27, 2 );
        memcpy( &version, data+29, 2 );
      }else if( is135v1 )
      {
        memcpy( &serialnum, data+28, 2 );
        memcpy( &version, data+30, 2 );
      }

      
      
      if( serialnum != 0 )
      {
        snprintf( charbuff, sizeof(charbuff), "%i", static_cast<int>(serialnum) );
        instrument_id_ = charbuff;
      }//if( serialnum != 0 )
      
      
      size_t calpos = 0;
      if( is130v0 || (is135v1 && record_size==2099) )
        calpos = 31;
      else if( is135v2 || (is135v1 && (record_size==2124)) )
        calpos = 44;
      
      if( calpos )
      {
        meas->calibration_coeffs_.resize( 3 );
        memcpy( &(meas->calibration_coeffs_[0]), data + calpos, 3*4 );
        meas->energy_calibration_model_ = Measurement::Polynomial;
      }//if( calpos )
      
      size_t datapos;
      uint16_t nchannels;
      
      if( is130v0 )
      {
        uint16_t real_time;
        uint32_t live_time_thousanths;
        nchannels = 256;
//        memcpy( &nchannels, data+5, 2 );
        memcpy( &real_time, data+14, 2 );
        memcpy( &live_time_thousanths, data+47, 4 );  //channels 1,2
        datapos = 47;
        
        const float lt = live_time_thousanths/1000.0f;
        const float rt = static_cast<float>( real_time );
        meas->real_time_ = std::max( rt, lt );
        meas->live_time_ = std::min( rt, lt );
        
        meas->contained_neutron_ = false;
      }else if( is135v1 )
      {
        memcpy( &nchannels, data+5, 2 );
        
        uint32_t live_time_thousanths;
        if( record_size == 2099 )
        {
          memcpy( &live_time_thousanths, data+50, 4 );
          datapos = 50;
        }else  // if( record_size == 2124 )
        {
          memcpy( &live_time_thousanths, data+75, 4 );
          datapos = 75;
        }
        
        meas->contained_neutron_ = false;
        meas->live_time_ = meas->real_time_ = live_time_thousanths / 1000.0f;
      }else //if( is135v2 )
      {
        /*
        char dose_unit = data[31];
        uint32_t dose;
        memcpy( &dose, data + 32, sizeof(dose) );
        
        switch( dose_unit )
        {
          case 'R': case 0:   //not sure how works...
          case 'G': case 1:
          case 'S': case 2:
          default:
            break;
        }
         */
        
        //uint16_t pilup_pulses;
        //memcpy( &pilup_pulses, data + 38, sizeof(pilup_pulses) );
        
        uint16_t nneutrons;
        uint32_t real_time_thousanths, live_time_thousanths;
        memcpy( &nchannels, data + 19, sizeof(uint16_t) );
        memcpy( &real_time_thousanths, data+21, sizeof(uint32_t) );
        memcpy( &nneutrons, data+36, sizeof(uint16_t) );
        memcpy( &live_time_thousanths, data+75, sizeof(uint32_t) ); // live time in mSec (channels 1,2)
        datapos = 75;
        
        
 
        
        const float lt = live_time_thousanths / 1000.0f;
        const float rt = real_time_thousanths / 1000.0f;
        meas->real_time_ = std::max( rt, lt );
        meas->live_time_ = std::min( rt, lt );
        
        meas->contained_neutron_ = true;
        meas->neutron_counts_sum_ = static_cast<double>( nneutrons );
        meas->neutron_counts_.resize( 1, static_cast<float>(nneutrons) );
      }//if( is130v0 ) / ....
      
      const float month = 30.0f*24.0f*60.0f*60.0f;
      if( meas->live_time_ < 0.1f || meas->live_time_ > month )
        meas->live_time_ = 0.0f;
      if( meas->real_time_ < 0.1f || meas->real_time_ > month )
        meas->real_time_ = 0.0f;
      
      if( is130v0 )
        meas->detector_type_ = "";
      else if( is135v1 )
        meas->detector_type_ = "";
      
      if( is130v0 )
         meas->detector_type_ = "GR-130";
      else if( is135v1 )
         meas->detector_type_ = "GR-135 v1";  // {Serial #'+SerialNumber+', '}  {Fnot invariant for same file}
      else if( is135v2 )
         meas->detector_type_ = "GR-135 v2";
      if( dtCzt )
        meas->detector_type_ += ", CZT";
    
      snprintf( charbuff, sizeof(charbuff),
                ", RecordSize: %i bytes", int(record_size) );
      meas->title_ += meas->detector_type_ + charbuff;
      
      meas->detector_number_ = dtCzt;
      meas->sample_number_ = static_cast<int>( j+1 );
      
      const size_t expected_num_channels = (is130v0 ? 256 : 1024);
      if( expected_num_channels != nchannels )
      {
        if( nchannels != 256 )
          nchannels = expected_num_channels;
        warnings.insert( "The expected and read number of channels didnt agree" );
      }
      
      std::shared_ptr< vector<float> > gamma_counts = std::make_shared<vector<float> >( nchannels, 0.0f );
      meas->gamma_counts_ = gamma_counts;
      vector<float> &channel_data = *gamma_counts.get();
      channel_data[0] = channel_data[1] = 0.0f;
      channel_data[channel_data.size()-1] = 0.0f;
      
	  if( nchannels >= 1 )
	  {
        for( uint16_t i = 2; i < (nchannels-1); ++i )
        {
          uint16_t counts;
          memcpy( &counts, data + datapos + 2*i, 2 );
          channel_data[i] = static_cast<float>( counts );
          meas->gamma_count_sum_ += counts;
        }//for( size_t i = 0; i < nchannels; ++i )
	  }

      if( is135v1 && meas->calibration_coeffs_.empty() )
      {
        const float nbin = static_cast<float>( nchannels );
        
        for( size_t i = 0; i < (record_size - 1024*2); ++i )
        {
          float cal[3];
          memcpy( cal, data+1, 12 );
          const float x = cal[0] + cal[1]*1.0f + cal[2]*1.0f*1.0f;
          const float y = cal[0] + cal[1]*nbin + cal[2]*nbin*nbin;
          const bool valid = (x > -100.0f) && (x < 100.0f) && (y > 400.0f) && (y < 4000.0f);
          if( valid )
          {
            meas->energy_calibration_model_ = Measurement::Polynomial;
            meas->calibration_coeffs_.push_back( cal[0] );
            meas->calibration_coeffs_.push_back( cal[1] );
            meas->calibration_coeffs_.push_back( cal[2] );
            warnings.insert( "Warning. Irregular GR energy calibration apparently found" );
            break;
          }//if( valid )
        }//for( size_t i = 0; i < (record_size - 1024*2); ++i )
      }//if( is135v1 && meas->calibration_coeffs_.empty() )
      
      if( meas->calibration_coeffs_.empty() )
      {
        float EgyCal[3] = { 0.0f };
        const float nbin = static_cast<float>( nchannels );
        
        if( is135v1 || is135v2 )
        {
          if( dtCzt )
          {
            EgyCal[1] = (122.06f - 14.4f)/(126.0f - 10.0f);
            EgyCal[0] = 14.4f - EgyCal[1]*10.0f;
            EgyCal[2] = 0.0f;
            warnings.insert( "Warning. Default GR135 energy calibration for CZT"
                             " has been assumed" );
          }else
          {
            EgyCal[0] = 0.11533801f;
            EgyCal[1] = 2.8760445f;
            EgyCal[2] = 0.0006023737f;
            warnings.insert( "Warning. Default GR135 energy calibration for NaI"
                             " has been assumed." );
          }
        }else if( is130v0 )
        {
          EgyCal[0] = -21.84f;
          EgyCal[1] = 3111.04f/(nbin+1.0f);
          EgyCal[2] = 432.84f/(nbin+1.0f)/(nbin+1.0f);
          warnings.insert( "Warning. Default GR130 energy calibration for NaI"
                           " has been assumed." );
        }
        
        meas->calibration_coeffs_.push_back( EgyCal[0] );
        meas->calibration_coeffs_.push_back( EgyCal[1] );
        meas->calibration_coeffs_.push_back( EgyCal[2] );
        meas->energy_calibration_model_ = Measurement::Polynomial;
      }//if( meas->calibration_coeffs_.empty() )
      
      if( j == 0 )
      {
        manufacturer_ = "Exploranium";
        instrument_model_ = is130v0 ? "GR130" : "GR135";
        instrument_type_ = "Radionuclide Identifier";
        if( !is130v0 )
          detector_type_ = kGR135Detector;
      }//if( j == 0 )
      
      measurements_.push_back( meas );
    }//for( size_t j = 0; j < recordstarts.size(); ++j )
    
    cleanup_after_load();
  }catch( std::exception & )
  {
    input.clear();
    input.seekg( orig_pos, ios::beg );

    reset();
    return false;
  }//try / catch
  
  const bool success = !measurements_.empty();
  
  if( success && warnings.size() )
  {
    for( const char *msg : warnings )
    {
      passMessage( msg, "load_binary_exploranium_file(...)", 1 );
#if(PERFORM_DEVELOPER_CHECKS)
      log_developer_error( BOOST_CURRENT_FUNCTION, msg );
#endif
    }
  }//if( success )
  
  return success;
}//void load_from_binary_exploranium()


bool MeasurementInfo::write_binary_exploranium_gr130v0( std::ostream &output ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const size_t nmeas = measurements_.size();
  
  try
  {
    int nwrote = 0;
  
    for( size_t measn = 0; measn < nmeas; ++measn )
    {
      const std::shared_ptr<const Measurement> origmeas = measurements_[measn];
      if( !origmeas || !origmeas->gamma_counts_ || origmeas->gamma_counts_->size() < 5 )
        continue;
      
      Measurement meas( *origmeas );
      
      const size_t noutchannel = 256;
      const size_t ninputchannel = meas.gamma_counts_->size();
      
      if( (ninputchannel > noutchannel) && ((ninputchannel % noutchannel)!=0) )
      {
        const size_t ncombine = ninputchannel / noutchannel;
        meas.combine_gamma_channels( ncombine );
      }else if( ninputchannel > noutchannel )
      {
        const float min_e = meas.gamma_energy_min();
        const float delta_e = meas.gamma_energy_max() - min_e;
        vector<float> coeffs;
        coeffs.push_back( min_e );
        coeffs.push_back( delta_e / noutchannel );
        
        std::shared_ptr<const vector<float> > new_binning =
        polynomial_binning( coeffs, noutchannel, DeviationPairVec() );
        
        meas.rebin_by_lower_edge( new_binning );
        meas.recalibrate_by_eqn( coeffs, DeviationPairVec(),
                                Measurement::Polynomial, new_binning );
      }//
      
      //Lets write to a buffer before writing to the stream to make things a
      //  little easier in terms of keeping track of position.
      char buffer[561];
      memset( buffer, 0, sizeof(buffer) );
      
      memcpy( buffer + 0, "ZZZZ", 4 );

      const uint16_t record_length = 560; //total length of the record
      memcpy( buffer + 4, &record_length, 2 );
      buffer[6] = 'A';

      if( !meas.start_time_.is_special() && meas.start_time_.date().year() > 2000 )
      {
        buffer[7] = static_cast<uint8_t>(meas.start_time_.date().year() - 2000);
        buffer[8] = static_cast<uint8_t>(meas.start_time_.date().month());
        buffer[9] = static_cast<uint8_t>(meas.start_time_.date().day());
        buffer[10] = static_cast<uint8_t>(meas.start_time_.time_of_day().hours());
        buffer[11] = static_cast<uint8_t>(meas.start_time_.time_of_day().minutes());
        buffer[12] = static_cast<uint8_t>(meas.start_time_.time_of_day().seconds());
        
        for( size_t i = 7; i < 13; ++i )
          buffer[i] = (((buffer[i] / 10) << 4) | (buffer[i] % 10));
      }//if( !meas.start_time_.is_special() )
      
      //battery volt*10 in BCD
      //buffer[13] = uint8_t(3.3*10);
      //buffer[13] = (((buffer[13] / 10) << 4) | (buffer[13] % 10));
      
      const uint16_t real_time = static_cast<uint16_t>( floor(meas.real_time_ + 0.5f) );
      memcpy( buffer + 14, &real_time, 2 );
      buffer[16] = (meas.gamma_energy_max() > 2000.0f ? 1 : 0); //coarse gain (0=1.5MeV, 1=3.0MeV)
      buffer[17] = 0; //Fgain: fine gain (0-255)
      //buffer[18] = buffer[19] = 0; //unsigned int		Peak		stab. peak position in channels *10
      //buffer[20] = buffer[21] = 0; //unsigned int		FW		stab. peak resolution *10 %
      buffer[22] = 'R'; //dose meter mode (R, G, S) {rem, grey, sievert}?
      //buffer[23] = buffer[24] = buffer[25] = buffer[26] = 0;//geiger		GM tube accumulated dose
      
      uint16_t serialint;
      if( !(stringstream(instrument_id_) >> serialint) )
        serialint = 0;
      memcpy( buffer + 27, &serialint, 2 );
      
      uint16_t softwareversion = 301;
      memcpy( buffer + 29, &softwareversion, 2 );
      
      buffer[31] = 'C'; //modification 	"C" Customs, "G" Geological
      
      //buffer[32] throughw buffer[46] 15 bytes spare
      
      const uint32_t ltime = static_cast<uint32_t>( floor( 1000.0f * meas.live_time() + 0.5f) );
      memcpy( buffer + 47, &ltime, 4 );
      
      const vector<float> &gammcounts = *meas.gamma_counts_;
      vector<uint16_t> channelcounts( gammcounts.size(), 0 );
      for( size_t i = 0; i < gammcounts.size(); ++i )
        channelcounts[i] = gammcounts[i];
      
      if( gammcounts.size() <= 256 )
        memcpy( buffer + 51, &channelcounts[2], (2*gammcounts.size() - 6) );
      else
        memcpy( buffer + 51, &channelcounts[2], 506 );
      
      const uint16_t cosmicchanel = 0; //cosmic		cosmic channel (channel 256)
      memcpy( buffer + 557, &cosmicchanel, 2 );
      
      buffer[559] = 0; //CHSUM		check sum
      
      output.write( buffer, 560 );
      
      ++nwrote;
    }//for( size_t measn = 0; measn < nmeas; ++measn )
    
    if( !nwrote )
      throw runtime_error("Failed to write any spectrums");
  }catch( std::exception &e )
  {
    return false;
  }
  
  return true;
}//bool write_binary_exploranium_gr130v0(...)



bool MeasurementInfo::write_binary_exploranium_gr135v2( std::ostream &output ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const size_t nmeas = measurements_.size();
  
  try
  {
    int nwrote = 0;
    
    for( size_t measn = 0; measn < nmeas; ++measn )
    {
      const std::shared_ptr<const Measurement> origmeas = measurements_[measn];
      if( !origmeas || !origmeas->gamma_counts_ || origmeas->gamma_counts_->size() < 5 )
        continue;
      
      Measurement meas( *origmeas );
      
      const uint16_t noutchannel = 1024;
      const size_t ninputchannel = meas.gamma_counts_->size();
      
      if( (ninputchannel > noutchannel) && ((ninputchannel % noutchannel)!=0) )
      {
        const size_t ncombine = meas.gamma_counts_->size() / noutchannel;
        meas.combine_gamma_channels( ncombine );
      }else if( ninputchannel > noutchannel )
      {
        const float min_e = meas.gamma_energy_min();
        const float delta_e = meas.gamma_energy_max() - min_e;
        vector<float> coeffs;
        coeffs.push_back( min_e );
        coeffs.push_back( delta_e / noutchannel );
        
        std::shared_ptr<const vector<float> > new_binning =
             polynomial_binning( coeffs, noutchannel, DeviationPairVec() );
        
        meas.rebin_by_lower_edge( new_binning );
        meas.recalibrate_by_eqn( coeffs, DeviationPairVec(),
                                Measurement::Polynomial, new_binning );
      }//
      
      
      //Lets write to a buffer before writing to the stream to make things a
      //  little easier in terms of keeping track of position.
      char buffer[76+2*noutchannel];
      memset( buffer, 0, sizeof(buffer) );
      
      memcpy( buffer + 0, "ZZZZ", 4 );
      memcpy( buffer + 4, "1350", 4 );
      
      const bool is_czt = UtilityFunctions::icontains(meas.detector_name_, "CZT");
      buffer[8] = (is_czt ? 'C' : 'A');
      
    
      //10-13		unsigned long		sequence number	number of spectra ever measured
      //  *this is not true to the input file*
      const uint32_t samplnum = meas.sample_number_;
      memcpy( buffer + 9, &samplnum, sizeof(samplnum) );
      
      if( !meas.start_time_.is_special() && meas.start_time_.date().year() > 2000 )
      {
        buffer[13] = static_cast<uint8_t>(meas.start_time_.date().year() - 2000);
        buffer[14] = static_cast<uint8_t>(meas.start_time_.date().month());
        buffer[15] = static_cast<uint8_t>(meas.start_time_.date().day());
        buffer[16] = static_cast<uint8_t>(meas.start_time_.time_of_day().hours());
        buffer[17] = static_cast<uint8_t>(meas.start_time_.time_of_day().minutes());
        buffer[18] = static_cast<uint8_t>(meas.start_time_.time_of_day().seconds());
        
        for( size_t i = 13; i < 19; ++i )
          buffer[i] = (((buffer[i] / 10) << 4) | (buffer[i] % 10));
      }//if( !meas.start_time_.is_special() )
      
      memcpy( buffer + 19, &noutchannel, sizeof(noutchannel) );
  
      const uint32_t real_time_thousanths = static_cast<uint32_t>( floor(1000*meas.real_time_ + 0.5f) );
      memcpy( buffer + 21, &real_time_thousanths, 4 );
      
      //26,27		unsigned int		gain			gain ( 0 - 1023)
      //28,29		unsigned int		Peak			stab. peak position in channels * 10
      //30,31		unsigned int		FW			stab. peak resolution * 10 %
      //32		unsigned char		unit			dose meter unit ('R', 'G', 'S')
      //33-36		unsigned long		geiger			GM tube accumulated dose
      
      const uint16_t nneutron = static_cast<uint16_t>( floor(meas.neutron_counts_sum() + 0.5) );
      memcpy( buffer + 36, &nneutron, 2 );
      
      //39,40		unsigned int		pileup			pileup pulses
      
      uint16_t serialint;
      if( !(stringstream(instrument_id_) >> serialint) )
        serialint = 0;
      memcpy( buffer + 40, &serialint, 2 );
      
      uint16_t softwareversion = 201;
      memcpy( buffer + 42, &softwareversion, 2 );
      
      if( meas.calibration_coeffs_.size() )
        memcpy( buffer + 44, &(meas.calibration_coeffs_[0]), 4*meas.calibration_coeffs_.size() );
      
      //57		char			temp[0]			display temperature
      //58		char			temp[1]			battery temperature
      //59		char			temp[2]			detector temperature
      //60-75		char			spare			16 bytes

      //live time isnt in the spec, but it was found here in the parsing code.
      const uint32_t live_time_thousanths = static_cast<uint32_t>( floor(1000*meas.live_time_ + 0.5f) );
      memcpy( buffer + 75, &live_time_thousanths, sizeof(uint32_t) ); // live time in mSec (channels 1,2)
      
      const vector<float> &counts = *meas.gamma_counts_;
      vector<uint16_t> channelcounts( counts.size(), 0 );
      for( size_t i = 0; i < counts.size(); ++i )
        channelcounts[i] = static_cast<uint16_t>( min( counts[i], 65535.0f ) );
      
      if( counts.size() <= noutchannel )
        memcpy( buffer + 79, &channelcounts[2], 2*counts.size() - 4 );
      else
        memcpy( buffer + 79, &channelcounts[2], 2*(noutchannel-2) );
      
      //2*nch+76	unsigned char		CHSUM		check sum
      
      output.write( buffer, 76+2*noutchannel );
      
      ++nwrote;
    }//for( size_t measn = 0; measn < nmeas; ++measn )
    
    if( !nwrote )
      throw runtime_error("Failed to write any spectrums");
  }catch( std::exception &e )
  {
    return false;
  }
  
  return true;
}//bool write_binary_exploranium_gr135v2(...)


#if( ENABLE_D3_CHART_EXPORTING )
bool MeasurementInfo::write_d3_html( ostream &ostr,
                                     const D3SpectrumExport::D3SpectrumChartOptions &options,
                                     std::set<int> sample_nums,
                                     const std::set<int> &det_nums ) const
{
  try
  {
    std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
    if( sample_nums.empty() )
      sample_nums = sample_numbers_;
  
    const size_t ndet = detector_numbers_.size();
    vector<bool> detectors( ndet, true );
    if( !det_nums.empty() )
    {
      for( size_t i = 0; i < ndet; ++i )
        detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
    }//if( det_nums.empty() )
  
  
    MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
  
    if( !summed || !summed->gamma_counts() || summed->gamma_counts()->empty() )
      return false;
  
    
    vector< pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    D3SpectrumExport::D3SpectrumOptions spec_options;
    measurements.push_back( pair<const Measurement *,D3SpectrumExport::D3SpectrumOptions>(summed.get(),spec_options) );
    
    return D3SpectrumExport::write_d3_html( ostr, measurements, options );
  }catch( std::exception &e )
  {
     return false;
  }
  
  return true;
}
#endif

bool MeasurementInfo::write_iaea_spe( ostream &output,
                                      set<int> sample_nums,
                                      const set<int> &det_nums ) const
{
  //www.ortec-online.com/download/ortec-software-file-structure-manual.pdf
  
  try
  {
    std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
    
    if( sample_nums.empty() )
      sample_nums = sample_numbers_;
    
    const size_t ndet = detector_numbers_.size();
    vector<bool> detectors( ndet, true );
    if( !det_nums.empty() )
    {
      for( size_t i = 0; i < ndet; ++i )
        detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
    }//if( det_nums.empty() )
    
    
    MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
    
    if( !summed || !summed->gamma_counts() || summed->gamma_counts()->empty() )
      return false;
    
    char buffer[256];
    
    string title = summed->title();
    UtilityFunctions::ireplace_all( title, "\r\n", " " );
    UtilityFunctions::ireplace_all( title, "\r", " " );
    UtilityFunctions::ireplace_all( title, "\n", " " );
    
    if( title.size() )
      output << "$SPEC_ID:\r\n" << title << "\r\n";
    
    vector<string> remarks = remarks_;
    remarks.insert( remarks.end(), summed->remarks_.begin(), summed->remarks_.end() );
    
    if( remarks.size() || title.size() )
    {
      output << "$SPEC_REM:\r\n";
      
      for( string remark : remarks )
      {
        UtilityFunctions::ireplace_all( remark, "\r\n", " " );
        UtilityFunctions::ireplace_all( remark, "\r", " " );
        UtilityFunctions::ireplace_all( remark, "\n", " " );
        if( remark.size() )
          output << remark << "\r\n";
      }
    }//if( remarks.size() )
    
    if( !summed->start_time_.is_special() )
    {
      // mm/dd/yyyy hh:mm:ss "02/29/2016 14:31:47"
      const int year = summed->start_time_.date().year();
      const int month = summed->start_time_.date().month();
      const int day = summed->start_time_.date().day();
      const int hour = summed->start_time_.time_of_day().hours();
      const int mins = summed->start_time_.time_of_day().minutes();
      const int secs = summed->start_time_.time_of_day().seconds();
      //double frac = summed->start_time_.time_of_day().fractional_seconds()
        //             / double(boost::posix_time::time_duration::ticks_per_second());
      
      //snprintf( buffer, sizeof(buffer), "%.2i/%.2i/%.4i %.2i:%.2i:%09.6f",
        //        month, day, year, hour, mins, (secs+frac) );
      snprintf( buffer, sizeof(buffer), "%.2i/%.2i/%.4i %.2i:%.2i:%.2i",
                 month, day, year, hour, mins, secs );

      output << "$DATE_MEA:\r\n" << buffer << "\r\n";
    }
    
    if( summed->real_time_ > 0.0f && summed->live_time_ > 0.0f )
    {
      //output << "$MEAS_TIM:\r\n" << summed->live_time_ << " " << summed->real_time_ << "\r\n";
      snprintf( buffer, sizeof(buffer), "%.5f %.5f", summed->live_time_, summed->real_time_ );
      output << "$MEAS_TIM:\r\n" << buffer << "\r\n";
    }

    const vector<float> &counts = *summed->gamma_counts_;
    if( counts.size() )
    {
      output << "$DATA:\r\n0 " << (counts.size()-1) << "\r\n";
      for( size_t i = 0; i < counts.size(); ++i )
        output << counts[i] << "\r\n";
    }//if( counts.size() )
    
    vector<float> coefs;
    if( summed->energy_calibration_model_ == Measurement::Polynomial )
      coefs = summed->calibration_coeffs_;
    else if( summed->energy_calibration_model_ == Measurement::FullRangeFraction )
      coefs = fullrangefraction_coef_to_polynomial( summed->calibration_coeffs_, counts.size() );
    
    if( coefs.size() )
    {
      output << "$ENER_FIT:\r\n";
      for( size_t i = 0; i < coefs.size(); ++i )
        output << (i?" ":"") << coefs[i];
      output << "\r\n";
      output << "$MCA_CAL:\r\n" << coefs.size() << "\r\n";
      for( size_t i = 0; i < coefs.size(); ++i )
        output << (i?" ":"") << coefs[i];
      output << "\r\n";
    }//if( coefs.size() )
    
    output << "$ENDRECORD:\r\n";
  }catch( std::exception & )
  {
    return false;
  }
  
  return true;
}//write_iaea_spe(...)

bool MeasurementInfo::load_pcf_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();


  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;

//  const size_t orig_pos = file.tellg();
//  file.seekg( 0, ios::end );
//  const istream::pos_type eof_pos = file.tellg();
//  file.seekg( orig_pos, ios::beg );
//  const size_t filelen = 0 + eof_pos - orig_pos;
//  string data;
//  data.resize( filelen );
//  file.read( &(data[0]), filelen );
//  stringstream strm( data );
//  const bool loaded = load_from_pcf( strm );

  const bool loaded = load_from_pcf( file );


  if( loaded )
    filename_ = filename;

  return loaded;
}//bool load_pcf_file( const std::string &filename )


const string SpectrumNodeDecodeWorker_failed_decode_title = "AUniqueStringToMarkThatThisDecodingFailed";
struct SpectrumNodeDecodeWorker
{
  const rapidxml::xml_node<char> *m_spec_node;
  std::mutex *m_mutex;
  std::shared_ptr<Measurement> m_meas;
  std::shared_ptr<DetectorAnalysis> m_analysis_info;
  const rapidxml::xml_node<char> *m_dose_data_parent;
  const rapidxml::xml_node<char> *m_doc;

  SpectrumNodeDecodeWorker( const rapidxml::xml_node<char> *node_in,
                            std::mutex *mutex_ptr,
                            std::shared_ptr<Measurement> meas,
                            std::shared_ptr<DetectorAnalysis> analysis_info_ptr,
                            const rapidxml::xml_node<char> *dose_data_parent,
                            const rapidxml::xml_node<char> *doc_node )
    : m_spec_node( node_in ),
      m_mutex( mutex_ptr ),
      m_meas( meas ),
      m_analysis_info( analysis_info_ptr ),
      m_dose_data_parent( dose_data_parent ),
      m_doc( doc_node )
  {}
  
  static void filter_valid_measurments( vector< MeasurementShrdPtr > &meass )
  {
    vector< MeasurementShrdPtr > valid_meass;
    valid_meass.reserve( meass.size() );
    
    for( auto &meas : meass )
    {
//      const bool hasGammaData = (meas->gamma_counts() && meas->gamma_counts()->size());
//      const bool hasNeutData = (meas->neutron_counts().size() || meas->contained_neutron_);
//      if( hasGammaData || hasNeutData )

      if( meas->title_ != SpectrumNodeDecodeWorker_failed_decode_title )
        valid_meass.push_back( meas );
    }//for( MeasurementShrdPtr &meas : meass )
    
    meass.swap( valid_meass );
  }

  void operator()()
  {
    try
    {
      const string xmlns = get_n42_xmlns(m_spec_node);
      
      m_meas->set_2006_N42_spectrum_node_info( m_spec_node );
      
      if( m_meas->calibration_coeffs().empty()
          && (!m_meas->channel_energies() || m_meas->channel_energies()->empty()) )
      {
        m_meas->set_n42_2006_spectrum_calibration_from_id( m_doc, m_spec_node );
      }

      if( m_dose_data_parent )
      {
        //If m_spec_node has any immediate siblings, then we need to be careful in setting
        //  The count dose informiaton, so lets count them
        int nspectra = 0;
        if( m_spec_node->parent() )
        {
          for( const rapidxml::xml_node<char> *node = m_spec_node->parent()->first_node( m_spec_node->name(), m_spec_node->name_size() );
               node;
               node = XML_NEXT_TWIN(node) )
            ++nspectra;
        }
        
        for( const rapidxml::xml_node<char> *dose_data = xml_first_node_nso( m_dose_data_parent, "CountDoseData", xmlns );
            dose_data;
            (dose_data = XML_NEXT_TWIN(dose_data)))
        {
          if( nspectra < 2 )
          {
            m_meas->set_n42_2006_count_dose_data_info( dose_data, m_analysis_info, m_mutex );
          }else
          {
            const rapidxml::xml_node<char> *starttime = XML_FIRST_NODE(dose_data, "StartTime");
            //const rapidxml::xml_node<char> *realtime = XML_FIRST_NODE(dose_data, "SampleRealTime");
            if( /*!realtime ||*/ !starttime )
              continue;
            
            const boost::posix_time::ptime startptime = UtilityFunctions::time_from_string( xml_value_str(starttime).c_str() );
            if( startptime.is_special() )
              continue;
            
            boost::posix_time::time_duration thisdelta = startptime - m_meas->start_time_;
            if( thisdelta < boost::posix_time::time_duration(0,0,0) )
              thisdelta = -thisdelta;
            //const float realtime = time_duration_in_seconds(realtime->value(), realtime->value_size());
            
            if( thisdelta < boost::posix_time::time_duration(0,0,10) )
              m_meas->set_n42_2006_count_dose_data_info( dose_data, m_analysis_info, m_mutex );
          }
        }//for( loop over CountDoseData nodes, dos )
      }//if( measurement_node_for_cambio )

      
      //HPRDS (see refF4TD3P2VTG) have start time and remark as siblings to
      //  m_spec_node
      const rapidxml::xml_node<char> *parent = m_spec_node->parent();
      if( parent )
      {
        const rapidxml::xml_node<char> *remark = xml_first_node_nso( parent, "Remark", xmlns );
        string remarkstr = xml_value_str( remark );
        trim( remarkstr );
        if( remarkstr.size() )
          m_meas->remarks_.push_back( remarkstr );
        
        const rapidxml::xml_node<char> *start_time = xml_first_node_nso( parent, "StartTime", xmlns );
        if( start_time && start_time->value_size() && m_meas->start_time_.is_special()
            && m_meas->source_type_ != Measurement::IntrinsicActivity )
          m_meas->start_time_ = time_from_string( xml_value_str(start_time).c_str() );
      }//if( parent )
    }
#if(PERFORM_DEVELOPER_CHECKS)
    catch( std::exception &e )
    {
      m_meas->reset();
      m_meas->title_ = SpectrumNodeDecodeWorker_failed_decode_title;
      if( !UtilityFunctions::icontains( e.what(), "didnt find <ChannelData>" ) )
      {
        char buffer[256];
        snprintf( buffer, sizeof(buffer), "Caught: %s", e.what() );
        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
      }//if( !UtilityFunctions::icontains( e.what(), "didnt find <ChannelData>" ) )
    }//try / catch
#else
    catch( std::exception & )
    {
      m_meas->reset();
      m_meas->title_ = SpectrumNodeDecodeWorker_failed_decode_title;
    }
#endif
  }//void operator()()
};//struct SpectrumNodeDecodeWorker




struct GrossCountNodeDecodeWorker
{
  const rapidxml::xml_node<char> *m_node;
  Measurement *m_meas;

  GrossCountNodeDecodeWorker( const rapidxml::xml_node<char> *node,
                              Measurement *newmeas )
    : m_node( node ),
      m_meas( newmeas )
  {}

  void operator()()
  {
    try
    {
      m_meas->set_n42_2006_gross_count_node_info( m_node );
    }catch( std::exception & )
    {
//      cerr << SRC_LOCATION << "\n\tCaught: " << e.what() << endl;
      m_meas->reset();
    }
  }//void operator()()
};//struct GrossCountNodeDecodeWorker


void MeasurementInfo::set_n42_2006_deviation_pair_info( const rapidxml::xml_node<char> *info_node,
                                            std::vector<MeasurementShrdPtr> &measurs_to_update )
{
  if( !info_node )
    return;

  set<string> det_names_set;

  for( const rapidxml::xml_node<char> *nl_corr_node = xml_first_node_nso( info_node, "NonlinearityCorrection", "dndons:" );
       nl_corr_node;
       nl_corr_node = XML_NEXT_TWIN(nl_corr_node) )
  {
    if( det_names_set.empty() )
      det_names_set = find_detector_names();

    const rapidxml::xml_attribute<char> *det_attrib = XML_FIRST_ATTRIB( nl_corr_node, "Detector" );

    if( !det_attrib && det_names_set.size()>1 )
    {
      cerr << SRC_LOCATION << "\n\tWarning, no Detector attribute "
           << "in <dndons:NonlinearityCorrection> node; skipping" << endl;
      continue;
    }//if( !det_attrib )

    const string det_name = det_attrib ? xml_value_str(det_attrib) : *(det_names_set.begin());
    bool have_det_with_this_name = false;
    for( const string &name : det_names_set )
    {
      //We may have appended "_intercal_" + "..." to the detector name, lets check
      have_det_with_this_name = (name == det_name || (istarts_with(name, det_name) && icontains(name, "_intercal_")));
      if( have_det_with_this_name )
        break;
    }
    
    if( !have_det_with_this_name )
    {
      cerr << SRC_LOCATION << "\n\tWarning: could find nedetctor name '"
           << det_name << "' in Measurements loaded, skipping deviation "
           << "pair" << endl;
      continue;
    }//if( an invalid detector name )

    vector< pair<float,float> > deviatnpairs;
    for( const rapidxml::xml_node<char> *dev_node = xml_first_node_nso( nl_corr_node, "Deviation", "dndons:" );
         dev_node;
         dev_node = XML_NEXT_TWIN(dev_node) )
    {
      if( dev_node->value_size() )
      {
        vector<float> devpair;
        const bool success = UtilityFunctions::split_to_floats( dev_node->value(), dev_node->value_size(), devpair );
        
        if( success && devpair.size()>=2 )
          deviatnpairs.push_back( pair<float,float>(devpair[0],devpair[1]) );
        else
          cerr << "Could not put '" << xml_value_str(dev_node) << "' into deviation pair" << endl;
      }//if( dev_node->value_size() )
    }//for( loop over <dndons:Deviation> )

    for( auto &meas : measurs_to_update )
      if( meas->detector_name_ == det_name ||
         (istarts_with(meas->detector_name_, det_name) && icontains(meas->detector_name_, "_intercal_")) )
        meas->deviation_pairs_ = deviatnpairs;
  }//for( loop over <dndons:NonlinearityCorrection> )
}//void set_n42_2006_deviation_pair_info(...)


void MeasurementInfo::set_n42_2006_instrument_info_node_info( const rapidxml::xml_node<char> *info_node )
{
  if( !info_node )
    return;
  
  string xmlns = get_n42_xmlns(info_node);
  if( xmlns.empty() && info_node->parent() )
    xmlns = get_n42_xmlns(info_node->parent());
  
  const rapidxml::xml_node<char> *type_node = xml_first_node_nso( info_node, "InstrumentType", xmlns );
  if( type_node && !XML_VALUE_COMPARE(type_node, "unknown") && !XML_VALUE_ICOMPARE(type_node, "Other") )
    instrument_type_ = xml_value_str( type_node );
  
  const rapidxml::xml_node<char> *manufacturer_node = xml_first_node_nso( info_node, "Manufacturer", xmlns );
  if( manufacturer_node && !XML_VALUE_COMPARE(manufacturer_node, "unknown") )
    manufacturer_ = xml_value_str( manufacturer_node );

  const rapidxml::xml_node<char> *model_node = xml_first_node_nso( info_node, "InstrumentModel", xmlns );
  if( model_node && !XML_VALUE_COMPARE(model_node, "unknown") )
    instrument_model_ = xml_value_str( model_node );

  const rapidxml::xml_node<char> *id_node = xml_first_node_nso( info_node, "InstrumentID", xmlns );
  if( id_node && !XML_VALUE_COMPARE(id_node, "unknown") )
    instrument_id_ = xml_value_str( id_node );

  const rapidxml::xml_node<char> *probe_node = xml_first_node_nso( info_node, "ProbeType", xmlns );
  if( probe_node && !XML_VALUE_COMPARE(probe_node, "unknown") )
  {
    const string val = xml_value_str( probe_node );
    vector<string> fields;
    UtilityFunctions::split( fields, val, "," );
    
    if( fields.size() == 1 )
    {
      trim( fields[0] );
      const size_t gamma_pos = fields[0].find( "Gamma Detector:" );
      const size_t neut_pos = fields[0].find( "Neutron Detector:" );
      
      if( gamma_pos != string::npos && neut_pos != string::npos
          && neut_pos > gamma_pos )
      {
        //identiFINDER 2 NGH make it here
        //"Gamma Detector: NaI 35x51 Neutron Detector: ³He tube 3He3/608/15NS"
        string gamma = fields[0].substr( 0, neut_pos );
        string neutron = fields[0].substr( neut_pos );
        
        UtilityFunctions::trim( gamma );
        UtilityFunctions::trim( neutron );
        
        remarks_.push_back( gamma );
        remarks_.push_back( neutron );
//        contained_neutron_ = true;
      }else
      {
        //identiFINDER 2 NG makes it here. "Gamma Detector: NaI 35x51"
        remarks_.push_back( fields[0] );
      }
    }else
    {
      //Sam940s make it here
      for( string &field : fields )
      {
        trim( field );
        if( UtilityFunctions::istarts_with( field, "Serial") )
          instrument_id_ += ", Probe " + field;
        else if( UtilityFunctions::istarts_with( field, "Type") )
          instrument_model_ += "," + field.substr(4);
      }//for( string &field : fields )
    }
  }//if( probe_node && !XML_VALUE_COMPARE(probe_node, "unknown") )
  
  
  const rapidxml::xml_node<char> *lane_number_node = xml_first_node_nso( info_node, "LaneNumber", xmlns );
  if( lane_number_node && lane_number_node->value_size() )
  {
    const string lanestr = xml_value_str( lane_number_node );
    
    char *pEnd = NULL;
    const char *lanecstr = lanestr.c_str();
    
    const long int val = strtol( lanecstr, &pEnd, 10 );
    
    if( val || (lanecstr != pEnd) )
      lane_number_ = static_cast<int>( val );
  }//if( lane_number_node && lane_number_node->value_size() )

  const rapidxml::xml_node<char> *inst_version = xml_first_node_nso( info_node, "InstrumentVersion", xmlns );
  if( inst_version && inst_version->value_size() )
  {
    //"DETH-008"
    //"4822400HHB-D"
    //"Firmware DETD-306 Software 3.1.11.27229"
    //"Hardware: 4C	Firmware: 5.00.54	Operating System: 1.2.040	Application: 2.37"
    vector<string> fields;
    const string value = xml_value_str( inst_version );
    
    vector<pair<string,string> > subcomponents;
    const size_t ntab = std::count(value.begin(), value.end(), '\t');
    const size_t nsemi = std::count(value.begin(), value.end(), ':');
    
    bool hassub = false;
    if( nsemi == (ntab+1) )
    {
      hassub = true;
      //A rough, and easily breakable hack for:
      //"Hardware: 4C	Firmware: 5.00.54	Operating System: 1.2.040	Application: 2.37"
      UtilityFunctions::split( fields, value, "\t");
      for( size_t i = 0; hassub && i < fields.size(); ++i )
      {
        vector<string> subfields;
        UtilityFunctions::split( subfields, fields[i], ":");
        if( subfields.size() == 2 )
          subcomponents.push_back( make_pair(subfields[0], subfields[1]) );
        else
          hassub = false;
      }
    }
    
    if( !hassub )
    {
      UtilityFunctions::split( fields, value, " \t:");
    
      if( fields.size() )
      {
        if( (fields.size() % 2) != 0 )
        {
          component_versions_.push_back( make_pair("System", xml_value_str(inst_version)) );
        }else
        {
          for( size_t i = 0; i < fields.size(); i += 2 )
          {
            string &name = fields[i];
            string &value = fields[i+1];
            UtilityFunctions::trim(name);
            UtilityFunctions::trim(value);
            UtilityFunctions::ireplace_all( name, ":", "" );
            component_versions_.push_back( make_pair(name,value) );
          }
        }
      }//if( fields.size() )
    }
  }//if( inst_version && inst_version->value_size() )
  
  //<Canberra:Version>2.0.0.8</Canberra:Version>
  inst_version =  XML_FIRST_NODE(info_node,"Canberra:Version");
  if( inst_version && inst_version->value_size() )
    component_versions_.push_back( make_pair("CanberraVersion", xml_value_str(inst_version)) );

  
  //RadSeeker HPRDS.  Grab the detector type and append to model for now...
  //  This is a bit of a hack, and should be improved.
  const rapidxml::xml_node<char> *det_setup = xml_first_node_nso( info_node, "DetectorSetup", "sym:" );
  if( det_setup )
  {
    for( const rapidxml::xml_node<char> *det = XML_FIRST_NODE( det_setup, "Detector" );
         det; det = XML_NEXT_TWIN(det) )
    {
      const rapidxml::xml_attribute<char> *type_attrib = XML_FIRST_ATTRIB( det, "Type" );
      if( type_attrib && XML_VALUE_ICOMPARE(type_attrib,"MCA") )
      {
        const rapidxml::xml_node<char> *id_settings = XML_FIRST_NODE(det, "IdentificationSettings");
        if( id_settings )
        {
          const rapidxml::xml_attribute<char> *Material = XML_FIRST_ATTRIB( id_settings, "Material" );
          const rapidxml::xml_attribute<char> *Size = XML_FIRST_ATTRIB( id_settings, "Size" );
          const rapidxml::xml_attribute<char> *Name = XML_FIRST_ATTRIB( id_settings, "Name" );
          if( Material || Size || Name )
          {
            //string val = "Gamma Detector: " + xml_value_str(Material) + " " + xml_value_str(Size) + " " + xml_value_str(Name);
            //remarks_.push_back( val );
            
            string val = xml_value_str(Material) + " " + xml_value_str(Size) + " " + xml_value_str(Name);
            UtilityFunctions::trim( val );
            UtilityFunctions::ireplace_all( val, "  ", " " );
            if( instrument_model_.size() )
              val = " " + val;
            instrument_model_ += val;
          }
        }
      }
    }
  }//if( det_setup )
  
}//void set_n42_2006_instrument_info_node_info( const rapidxml::xml_node<char> *info_node )



void Measurement::set_n42_2006_count_dose_data_info( const rapidxml::xml_node<char> *dose_data,
                                        std::shared_ptr<DetectorAnalysis> analysis_info,
                                        std::mutex *analysis_mutex )
{
  if( !dose_data )
    return;

  string xmlns = get_n42_xmlns(dose_data);
  if( xmlns.empty() && dose_data->parent() )
    xmlns = get_n42_xmlns(dose_data->parent());
  
  Measurement *m = this;
  
  //Detective N42 2006 files have a CountRate node that you need to multiply
  //  by the live time
  
  const rapidxml::xml_node<char> *count_node = xml_first_node_nso( dose_data, "CountRate", xmlns );
  const rapidxml::xml_node<char> *realtime_node = xml_first_node_nso( dose_data, "SampleRealTime", xmlns );
  
  
  const rapidxml::xml_attribute<char> *det_attrib = XML_FIRST_ATTRIB(dose_data, "DetectorType");
  
  if( count_node && count_node->value_size()
     && (!det_attrib || XML_VALUE_ICOMPARE(det_attrib, "Neutron")) )
  {
    try
    {
      if( !realtime_node || !realtime_node->value_size() )
        throw runtime_error( "Couldnt find realtime for neutron count rate" );
      
      const float realtime = time_duration_in_seconds( realtime_node->value(),
                                                  realtime_node->value_size() );
      if( realtime <= 0.0f )
        throw runtime_error( "Couldnt read realtime" );
      
      rapidxml::xml_attribute<char> *units_attrib = count_node->first_attribute( "Units", 5 );
      if( units_attrib && !units_attrib->value_size() )
        units_attrib = 0;
      
      if( units_attrib && units_attrib->value_size()
         && !UtilityFunctions::icontains( xml_value_str(units_attrib), "CPS") )
        throw runtime_error( "Neutron count rate not in CPS" );
      
      float countrate;
      if( !xml_value_to_flt(count_node, countrate) )
        throw runtime_error( "Neutron count rate is non-numeric" );
      
      m->neutron_counts_sum_ = countrate * realtime;
      m->neutron_counts_.resize(1);
      m->neutron_counts_[0] = countrate * realtime;
      m->contained_neutron_ = true;
      m->remarks_.push_back( "Neutron Real Time: " + xml_value_str(realtime_node) );
      
      if( (m->real_time_ > FLT_EPSILON)
          && fabs(m->live_time_ - realtime) > 0.1f*m->live_time_ )
      {
        const string msg = "Warning: The neutron live time may not correspond to the gamma live time.";
        passMessage( msg.substr(9), "", 3 );
        if( std::find(m->remarks_.begin(),m->remarks_.end(),msg) == m->remarks_.end() )
          m->remarks_.push_back( msg );
      }
      
      rapidxml::xml_node<char> *starttime_node = dose_data->first_node( "StartTime", 9 );
      if( starttime_node && starttime_node->value_size() )
      {
        boost::posix_time::ptime starttime = time_from_string( xml_value_str(starttime_node).c_str() );
        if( !starttime.is_special() && !m->start_time_.is_special() )
        {
          if( (starttime-m->start_time_) > boost::posix_time::minutes(1) )
            m->remarks_.push_back( "Warning: neutron start time doesnt match gamma start time!" );
        }
      }//if( starttime_node && starttime_node->value_size() )
    }catch( std::exception &e )
    {
      m->remarks_.push_back( "Error decoding neutron count rate: " + string(e.what()) );
      //Could probably give more information here, but I doubt it will ever actually
      //  occur, so not bothering.
    }//try / catch
  }//if( we have count rate and sample real time nodes )
  
  
  if( !det_attrib )
    return;

  if( XML_VALUE_ICOMPARE(det_attrib, "Neutron") )
  {
    const rapidxml::xml_node<char> *counts = xml_first_node_nso( dose_data, "Counts", xmlns );
    if( counts && counts->value_size() )
    {
      float neut_counts = 0.0;
      if( xml_value_to_flt(counts,neut_counts) )
      {
        m->neutron_counts_sum_ += neut_counts;
        if( m->neutron_counts_.empty() )
          m->neutron_counts_.push_back( neut_counts );
        else if( m->neutron_counts_.size() == 1 )
          m->neutron_counts_[0] += neut_counts;
        else
          cerr << "Have both neutron spectrum and neutron dose count" << endl;

        //XXX
        //  We get this dose rate info from Cambio wether or not there was
        //  a neutron detector, so we will only mark the detector as a neutron
        //  detector if at least one sample has at least one neutron
        //   ... I'm not a huge fan of this...
        m->contained_neutron_ |= (m->neutron_counts_[0]>0.0);
      }else
        cerr << "Error converting neutron counts '" << xml_value_str(counts)
             << "' to float; ignoring" << endl;
    }//if( counts )
  }else if( XML_VALUE_ICOMPARE(det_attrib, "Gamma") )
  {
    const rapidxml::xml_node<char> *remark_node = xml_first_node_nso( dose_data, "Remark", xmlns );
    //const rapidxml::xml_node<char> *start_time_node = xml_first_node_nso( dose_data, "StartTime", xmlns );
    const rapidxml::xml_node<char> *real_time_node = xml_first_node_nso( dose_data, "SampleRealTime", xmlns );
    const rapidxml::xml_node<char> *dose_node = xml_first_node_nso( dose_data, "DoseRate", xmlns );
    
    DetectorAnalysisResult thisana;
    thisana.remark_ = xml_value_str( remark_node );
    //thisana.start_time_ = time_from_string( xml_value_str(start_time_node).c_str() );
    if( real_time_node && real_time_node->value_size() )
      thisana.real_time_ = time_duration_in_seconds( real_time_node->value(),
                                                real_time_node->value_size() );
    if( dose_node && dose_node->value_size() )
    {
      try
      {
        const rapidxml::xml_attribute<char> *units_attrib = XML_FIRST_ATTRIB( dose_node, "Units" );
        if( !units_attrib || !units_attrib->value_size() )
          throw runtime_error( "No units for dose" );
        
        xml_value_to_flt( dose_node, thisana.dose_rate_ );
        thisana.dose_rate_ *= dose_units_usvPerH( units_attrib->value(), units_attrib->value_size() );
      }catch(...){}
    }//if( dose_node && dose_node->value_size() )
    
    if( thisana.nuclide_.size()
       //|| !thisana.start_time_.is_special()
       || thisana.activity_ > 0.0
       //|| thisana.remark_.size()
       || thisana.nuclide_type_.size()
       || thisana.id_confidence_.size()
       //|| thisana.detector_.size()
       || thisana.distance_ > 0.0
       || thisana.dose_rate_ > 0.0 )
    {
      std::unique_ptr< std::lock_guard<std::mutex> > lock;
      if( analysis_mutex )
        lock.reset( new std::lock_guard<std::mutex>( *analysis_mutex ) );
      analysis_info->results_.push_back( thisana );
    }//if( 
  }//if( neutron detector ) / else( gamma detector )


}//void set_n42_2006_count_dose_data_info( const rapidxml::xml_node<char> *dose_data, MeasurementShrdPtr m )



void Measurement::set_n42_2006_gross_count_node_info( const rapidxml::xml_node<char> *gross_count_meas )
{
  if( !gross_count_meas )
    throw runtime_error( "!gross_count_measurement" );

  string xmlns = get_n42_xmlns(gross_count_meas);
  if( xmlns.empty() && gross_count_meas->parent() )
    xmlns = get_n42_xmlns(gross_count_meas->parent());
  
  Measurement *m = this;

  bool is_neuteron = m->contained_neutron_;

  if( !is_neuteron )
  {
    const rapidxml::xml_attribute<char> *det_type_attrib = XML_FIRST_ATTRIB( gross_count_meas, "DetectorType" );

    if( !det_type_attrib )
    {
      const rapidxml::xml_node<char> *parent = gross_count_meas->parent();
      if( parent )
        det_type_attrib = XML_FIRST_ATTRIB( parent, "DetectorType" );
    }//if( !det_type_attrib )
    if( det_type_attrib )
      is_neuteron = icontains( xml_value_str(det_type_attrib), "Neutron" );
  }//if( !is_neuteron )

  if( !is_neuteron )
    throw runtime_error( "!is_neuteron" );

  if( m->neutron_counts_sum_ > 0.0001 )
    throw runtime_error( "m->totalNeutronCounts > 0.0001" );

  float nprev = 0.0;
  for( size_t i = 0; i < m->neutron_counts_.size(); ++i )
    nprev +=  m->neutron_counts_[i];

  if( nprev > 0.0001 )
    throw runtime_error( "nprev > 0.0001" );

  m->contained_neutron_ = true;
  m->neutron_counts_.resize( 1 );
  m->neutron_counts_[0] = 0.0;

  const rapidxml::xml_node<char> *node = xml_first_node_nso( gross_count_meas, "GrossCounts", xmlns );

  if( node )
  {
    if( UtilityFunctions::split_to_floats( node->value(), node->value_size(), m->neutron_counts_ ) )
      m->neutron_counts_sum_ = std::accumulate( m->neutron_counts_.begin(), m->neutron_counts_.end(), 0.0f, std::plus<float>() );
    else
      m->neutron_counts_sum_ = 0.0;
  }//if( node )


  //attempt to set detctor name
  if( m->detector_name_.empty() )
  {
    const rapidxml::xml_attribute<char> *det_atrrib = NULL;
    for( node = gross_count_meas; !det_atrrib && node; node = node->parent() )
      det_atrrib = XML_FIRST_ATTRIB( node, "Detector" );
    if( det_atrrib )
      m->detector_name_ = xml_value_str( det_atrrib );
  }//if( m->detector_name_.empty() )

}//void set_n42_2006_gross_count_node_info(...)


void Measurement::set_n42_2006_detector_data_node_info( const rapidxml::xml_node<char> *det_data_node,
                              std::vector<MeasurementShrdPtr> &measurs_to_update )
{
  string xmlns = get_n42_xmlns(det_data_node);
  if( xmlns.empty() && det_data_node && det_data_node->parent() )
    xmlns = get_n42_xmlns(det_data_node->parent());
  
  const rapidxml::xml_node<char> *speed_node = xml_first_node_nso( det_data_node, "Speed", xmlns );
  const rapidxml::xml_node<char> *occupancy_node = xml_first_node_nso( det_data_node, "Occupied", xmlns );
  const rapidxml::xml_node<char> *start_time_node = xml_first_node_nso( det_data_node, "StartTime", xmlns );
  const rapidxml::xml_node<char> *sample_real_time_node = xml_first_node_nso( det_data_node, "SampleRealTime", xmlns );

  float real_time = 0.0, speed = 0.0;
  boost::posix_time::ptime start_time;
  Measurement::OccupancyStatus occupied = Measurement::UnknownOccupancyStatus;

  if( sample_real_time_node && sample_real_time_node->value_size() )
    real_time = time_duration_in_seconds( sample_real_time_node->value(), sample_real_time_node->value_size() );

  if( start_time_node )
     start_time = time_from_string( xml_value_str(start_time_node).c_str() );

  try{
    speed = speed_from_node( speed_node );
  }catch(...){}

  try{
    if( is_occupied( occupancy_node ) )
      occupied = Measurement::Occupied;
    else
      occupied = Measurement::NotOccupied;
  }catch(...){}


  for( auto &meas : measurs_to_update )
  {
    if( meas->occupied_ == Measurement::UnknownOccupancyStatus )
      meas->occupied_ = occupied;

    if( meas->speed_ < 0.00000001f )
      meas->speed_ = speed;

    //dont set the start time for identiFINDER IntrinsicActivity spectrum (which
    //  doesnt have the time in the file), based on the StartTime node under
    //  <DetectorData> (of which, this time disagrees with the actual spectrum
    //  StartTime, so I dont know what to do about this anyway)
    if( meas->start_time_.is_special()
        && (meas->source_type_ != Measurement::IntrinsicActivity) )
      meas->start_time_ = start_time;

    if( meas->real_time_ < 0.000001f )
      meas->real_time_ = real_time;
    
    //For neutron only detectors, a lot of times live time isnt given, so fill
    //  this out
    if( meas->contained_neutron_ && meas->live_time_ < 0.000001f
       && (!meas->gamma_counts_ || meas->gamma_counts_->empty()) )
      meas->live_time_ = meas->real_time_;
  }//foeach( MeasurementShrdPtr &meas, measurements_this_node )
}//void set_n42_2006_detector_data_node_info(  )


void MeasurementInfo::set_n42_2006_measurment_location_information(
                          const rapidxml::xml_node<char> *measured_item_info_node,
                          std::vector<MeasurementShrdPtr> added_measurements )
{
  if( !measured_item_info_node )
    return;
  
  string xmlns = get_n42_xmlns(measured_item_info_node);
  if( xmlns.empty() && measured_item_info_node->parent() )
    xmlns = get_n42_xmlns(measured_item_info_node->parent());
  
  vector<string> remarks;

  //XXX - get the remarks - but for right now just put into remarks_
  for( const rapidxml::xml_node<char> *remark_node = xml_first_node_nso( measured_item_info_node, "Remark", xmlns );
       remark_node;
       remark_node = XML_NEXT_TWIN(remark_node) )
  {
    string remark = xml_value_str(remark_node);
    trim( remark );
    if(!remark.empty())
      remarks_.push_back( remark );
  }
  
  bool read_in_gps = false;
  double latitude = -999.9;
  double longitude = -999.9;
  boost::posix_time::ptime position_time;
    
  const rapidxml::xml_node<char> *meas_loc_name = NULL;
  const rapidxml::xml_node<char> *meas_loc = xml_first_node_nso( measured_item_info_node, "MeasurementLocation", xmlns );
    
  if( !meas_loc )
    meas_loc = xml_first_node_nso( measured_item_info_node, "InstrumentLocation", xmlns );
  
  if( !meas_loc && XML_NAME_ICOMPARE(measured_item_info_node, "InstrumentLocation") )
    meas_loc = measured_item_info_node;
    
    
  if( meas_loc )
  {
    meas_loc_name = xml_first_node_nso( meas_loc, "MeasurementLocationName", xmlns );
    const rapidxml::xml_node<char> *coord_node = xml_first_node_nso( meas_loc, "Coordinates", xmlns );
      
    //We should actually loop over Coordinate nodes here and average values
    //  that occur while the measurment is being taken...
    if( coord_node && coord_node->value() )
    {
      if( (stringstream(xml_value_str(coord_node)) >> latitude >> longitude) )
      {
        //If there are a third coordinate is there, it is elevation in meters.
        //
        const rapidxml::xml_base<char> *time_element = XML_FIRST_ATTRIB( coord_node, "Time" );
        
        if( !time_element ) //Raytheon Portal
          time_element = xml_first_node_nso(meas_loc, "GPSDateTime", "ray:");
        
        if( time_element )
          position_time = time_from_string( xml_value_str(time_element).c_str() );
        
        //const rapidxml::xml_base<char> *delusion_element = xml_first_node_nso(meas_loc, "DelusionOfPrecision", "ray:");
      }
    }//if( coord_node )
  }//if( meas_loc )
    
  //I've notices some Detecive-EX100 have coordiantes specified similar to:
  //  <Coordinates Time="2017-08-18T10:33:52-04:00">3839.541600 -7714.840200 32</Coordinates>
  //  Which are actuall 38 degrees 39.541600 minutes North, 77 degrees 14.840200 minutes West, and an elevation of 32m
  //  which is 38.659028, -77.247333
  //  Unfortuanetly I havent yet noticed a solid tell in the file for this besides invalid coordinates
  if( !Measurement::valid_latitude(latitude)
      && !Measurement::valid_longitude(longitude)
      && (fabs(latitude)>999.99)
      && (fabs(longitude)>999.99) )
  {
    //As of 20170818 I have only tested this code for a handfull of points in the US
    //  (all positive latitude, and negative longitude).
    double lat_deg = floor( fabs(latitude) / 100.0);
    double lon_deg = floor( fabs(longitude) / 100.0);
    
    lat_deg += (fabs(latitude) - 100.0*lat_deg) / 60.0;
    lon_deg += (fabs(longitude) - 100.0*lon_deg) / 60.0;
    
    lat_deg *= (latitude > 0 ? 1.0 : -1.0);
    lon_deg *= (longitude > 0 ? 1.0 : -1.0);
    
    if( Measurement::valid_latitude(lat_deg) && Measurement::valid_longitude(lon_deg) )
    {
      latitude = lat_deg;
      longitude = lon_deg;
    }
  }//if( invalid lat/long, but maybe specified in a wierd format )
    
  if( Measurement::valid_latitude(latitude)
      && Measurement::valid_longitude(longitude) )
  {
    for( auto &meas : added_measurements )
    {
      read_in_gps          = true;
      meas->latitude_      = latitude;
      meas->longitude_     = longitude;
      meas->position_time_ = position_time;
    }//for( MeasurementShrdPtr &meas : measurements_this_node )
  }//if( !valid(latitude) &&  !valid(longitude) )
    
  if( !meas_loc_name )
    meas_loc_name = xml_first_node_nso( measured_item_info_node, "MeasurementLocationName", xmlns );
  measurement_location_name_ = xml_value_str( meas_loc_name );
    
  const rapidxml::xml_node<char> *operator_node = xml_first_node_nso( measured_item_info_node, "MeasurementOperator", xmlns );
  measurment_operator_ = xml_value_str( operator_node );
}//void set_n42_2006_measurment_location_information()



void MeasurementInfo::load_2006_N42_from_doc( const rapidxml::xml_node<char> *document_node )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  //CambioN42, ORTEC IDM, Thermo: document_node == "N42InstrumentData"
  
  //CambioN42: I think has multiple <Measurement> nodes
  //          The <Measurement> node then has daughter <Spectrum>, and possibly <CountDoseData>
  //            The <Spectrum> has daughters <StartTime>, <LiveTime>, <RealTime>, <Calibration>, and <ChannelData>, <Cambio:TagChar>, <Cambio:OriginalFileType>, <Cambio:Version>
  //          Note that cambio n42 files does necassarily preserve detector name, sample number, speed, etc
  //ORTEC IDM, Thermo, SAIC8 all have single <Measurement> node, which has daughters <InstrumentInformation>, <InstrumentType>, <Manufacturer>, <InstrumentModel>, <InstrumentID>, <LaneNumber>, <dndons:DetectorStatus>, <MeasuredItemInformation>, and <DetectorData>
  //           There is a <DetectorData> node for each timeslice, with daughtes <StartTime>, <SampleRealTime>, <Occupied>, <Speed>, and <DetectorMeasurement>
  //            The <DetectorMeasurement> node then has daughter <SpectrumMeasurement>
  //              The <SpectrumMeasurement> node then has daughters <SpectrumAvailable>, <Spectrum>
  //                There are then a <spectrum> node for each detector, which has the daughters <RealTime>, <LiveTime>, <SourceType>, <DetectorType>, <Calibration>, and <ChannelData>
  
  if( !document_node )
    throw runtime_error( "load_2006_N42_from_doc: Invalid document node" );
  
  //is_cambio actually indicates a "spectrometer" type file, see
  //  http://physics.nist.gov/Divisions/Div846/Gp4/ANSIN4242/2005/simple.html
  //  for example nucsafe g4 predators also go down this path
  const rapidxml::xml_node<char> *firstmeas = document_node->first_node( "Measurement", 11 );
  const bool is_cambio = (firstmeas && firstmeas->first_node( "Spectrum", 8 ));
  
  std::mutex queue_mutex;
  std::shared_ptr<DetectorAnalysis> analysis_info( new DetectorAnalysis() );
  
  SpecUtilsAsync::ThreadPool workerpool;
  
  string xmlns = get_n42_xmlns(document_node);
  if( xmlns.empty() )
  {
    for( auto attrib = document_node->first_attribute(); attrib; attrib = attrib->next_attribute() )
    {
      //Some files use the "n42ns" namespace, IDK
      const string name = xml_name_str(attrib);
      if( UtilityFunctions::starts_with(name, "xmlns:" ) && UtilityFunctions::icontains(name, "n42") )
        xmlns = name.substr(6) + ":";
    }//for( check for xmlns:n42ns="..."
  }//if( doc_node_name has a namespace in it ) / else
  
  if( xmlns.empty() )
    xmlns = "n42:";  //might not actually be needed, but JIC until I bother to test.
  
  
  if( is_cambio )
  {
    vector<const rapidxml::xml_node<char> *> countdoseNodes;
    vector<const rapidxml::xml_node<char> *> locationNodes, instInfoNodes;
    
    for( const rapidxml::xml_node<char> *measurement = xml_first_node_nso( document_node, "Measurement", xmlns );
        measurement;
        measurement = XML_NEXT_TWIN(measurement) )
    {
      
      for( const rapidxml::xml_node<char> *spectrum = xml_first_node_nso( measurement, "Spectrum", xmlns );
          spectrum;
          spectrum = XML_NEXT_TWIN(spectrum) )
      {
        std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
        measurements_.push_back( meas );
        workerpool.post( SpectrumNodeDecodeWorker( spectrum, &queue_mutex,
                                                  meas, analysis_info,
                                                  measurement, document_node ) );
      }
      
      for( const rapidxml::xml_node<char> *dose = xml_first_node_nso( measurement, "CountDoseData", xmlns );
          dose;
          dose = XML_NEXT_TWIN(dose) )
      {
        countdoseNodes.push_back( dose );
        //We dont want to set file level information from this <Measurment> node
        //  so we will continue;
        //continue;
      }//if( spectrum ) / else ( dose )
      
      //XML files from "Princeton Gamma-Tech Instruments, Inc." detectors
      //  will have the InstrumentInformation, AnalysisResults, and
      //  Calibration here
      const rapidxml::xml_node<char> *inst_info = xml_first_node_nso( measurement, "InstrumentInformation", xmlns );
      set_n42_2006_instrument_info_node_info( inst_info );
      instInfoNodes.push_back( inst_info );
      
      const rapidxml::xml_node<char> *analysis_node = xml_first_node_nso( measurement, "AnalysisResults", xmlns );
      if( analysis_node )
        setAnalysisInformation( analysis_node, analysis_info );
      
      //Try to find the MeasuredItemInformation node - this is all just barely hacked in
      const rapidxml::xml_node<char> *item_info_node = xml_first_node_nso( measurement, "MeasuredItemInformation", xmlns );
      if( !item_info_node && inst_info )
        item_info_node = xml_first_node_nso( inst_info, "MeasuredItemInformation", xmlns );
      
      if( !item_info_node )  //HPRDS files
        item_info_node = xml_first_node_nso( measurement, "InstrumentLocation", xmlns );
      
      locationNodes.push_back( item_info_node );
      
      //This should have already been taken care of by inst_info
//      const rapidxml::xml_node<char> *info_node = measurement->first_node( "InstrumentInformation", 21 );
//      set_n42_2006_instrument_info_node_info( info_node );
    }//for( loop over <Measurement> nodes )
    
    workerpool.join();
    
    SpectrumNodeDecodeWorker::filter_valid_measurments( measurements_ );
    
    //This loop may be taken care of by SpectrumNodeDecodeWorker anyway
    for( size_t i = 0; i < countdoseNodes.size(); ++i )
    {
      //Detectives N42 are what end up here.
      //Detectives list the neutron information in "CountDoseData" nodes
      //  If this is the case, we are ignoring <InstrumentInformation>,
      //  <MeasuredItemInformation> nodes.  In
      const rapidxml::xml_node<char> *dose = countdoseNodes[i];

      const rapidxml::xml_attribute<char> *dettype = XML_FIRST_ATTRIB(dose, "DetectorType");
      if( dettype && !XML_VALUE_ICOMPARE(dettype, "Neutron") )
        continue;
      
      bool addedCountDose = false;
      
      //find the nearest measurements and set intfo from it
      if( measurements_.size() == 1 )
      {
        //If only one measurment, set the data no matter whate
        measurements_[0]->set_n42_2006_count_dose_data_info( dose, analysis_info, &queue_mutex );
        addedCountDose = true;
      }else if( measurements_.size() )
      {
        const rapidxml::xml_node<char> *starttime_node = xml_first_node_nso( dose, "StartTime", xmlns );
        if( starttime_node && starttime_node->value_size() )
        {
          const boost::posix_time::ptime starttime = time_from_string( xml_value_str(starttime_node).c_str() );
          if( !starttime.is_special() )
          {
            int nearestindex = -1;
            boost::posix_time::time_duration smallestdelta = boost::posix_time::time_duration(10000,0,0);
            for( size_t j = 0; j < measurements_.size(); j++ )
            {
              const boost::posix_time::ptime &thisstartime = measurements_[j]->start_time_;
              boost::posix_time::time_duration thisdelta = thisstartime - starttime;
              if( thisdelta < boost::posix_time::time_duration(0,0,0) )
                thisdelta = -thisdelta;
              if( thisdelta < smallestdelta )
              {
                smallestdelta = thisdelta;
                nearestindex = (int)j;
              }
            }//for( size_t j = 1; j < measurements_.size(); j++ )
            
            //Make sure the start time maches within 1 minute (arbitrary), there
            //  is apparently a mode in detectives where the start time of
            //  neutrons can be well before the gamma, which causes confusion.
            if( nearestindex >= 0
                && smallestdelta < boost::posix_time::time_duration(0,1,0) )
            {
              if( !measurements_[nearestindex]->contained_neutron_ )
              {
                measurements_[nearestindex]->set_n42_2006_count_dose_data_info( dose, analysis_info, &queue_mutex );
                addedCountDose = true;
              }
              
              if( measurements_.size() == 2
                 && measurements_[nearestindex]->source_type() == Measurement::Foreground
                 && measurements_[nearestindex?0:1]->source_type() == Measurement::Background )
              {
                //For nucsafe g4 predator
                const rapidxml::xml_attribute<char> *det_attrib = XML_FIRST_ATTRIB(dose, "DetectorType");
                const rapidxml::xml_node<char> *backrate_node = xml_first_node_nso(dose, "BackgroundRate", "Nucsafe:");
                MeasurementShrdPtr back = measurements_[nearestindex?0:1];
                
                if( !back->contained_neutron_
                   && det_attrib && backrate_node && XML_VALUE_ICOMPARE(det_attrib, "Neutron") )
                {
                  MeasurementShrdPtr back = measurements_[nearestindex?0:1];
                  float rate;
                  if( xml_value_to_flt( backrate_node, rate) )
                  {
                    rate *= back->real_time_;
                    back->contained_neutron_ = true;
                    back->neutron_counts_.resize(1);
                    back->neutron_counts_[0] = rate;
                    back->neutron_counts_sum_ = rate;
                  }
                }//( if( the foregorund neutron info als ocontains background neutron rate )
              }//if( there is a foreground and background spectrum for this measurment, and we're assigning the neutron info to foreground )
            }//if( found an appropriate measurment to put this neutron info in )
          }//if( !starttime.is_special() )
        }//if( starttime_node && starttime_node->value_size() )
      }//if( we only have one other measurment ) / else, find nearest one
      
//#if(PERFORM_DEVELOPER_CHECKS)
//      if( !addedCountDose )
//        log_developer_error( BOOST_CURRENT_FUNCTION, "Failed to add count dose data!" );
//#endif
    }//for( size_t i = 0; i < countdoseNodes.size(); ++i )
    
    
    //Add in GPS and other location information.
    //Note that this is a hack to read in information saved by this same class
    //  in 2006 N42 files,
    //I've also seen this for refOAD5L6QOTM, where this construct is actually
    //  wrong and we should matchup the GPS coordinated using a time based
    //  approach
    for( size_t i = 0; i < locationNodes.size(); ++i )
    {
      if( locationNodes[i] && (i<measurements_.size()) && measurements_[i] )
      {
        vector<MeasurementShrdPtr> addedmeas( 1, measurements_[i] );
        set_n42_2006_measurment_location_information( locationNodes[i], addedmeas );
      }//if( locationNodes[i] && measurements_[i] )
    }//for( size_t i = 0; i < locationNodes.size(); ++i )
    
    //Add in deviation pairs for the detectors.
    //Note that this is a hack to read in deviation pairs saved from this code,
    //  and I havent seen any simple spectrometer style files actually use this.
    //Also, it is completely untested if it actually works.
    //For the ref8TINMQY7PF, the bellow index based matching up shows this
    //  approach may not always be correct...
    for( size_t i = 0; i < instInfoNodes.size(); ++i )
    {
      if( instInfoNodes[i] && (i<measurements_.size()) && measurements_[i] )
      {
        vector<MeasurementShrdPtr> addedmeas( 1, measurements_[i] );
        set_n42_2006_deviation_pair_info( instInfoNodes[i], addedmeas );
      }//if( instInfoNodes[i] && measurements_[i] )
    }//for( size_t i = 0; i < locationNodes.size(); ++i )
    
  }else
  {
    //This bellow loop could be parallized a bit more, in terms of all
    //  <Spectrum> tags of all <Measurement>'s being put in the queue to work
    //  on asyncronously.  In practice though this isnt necassarry since it
    //  looks like most passthrough files only have one <Measurement> tag anyway
    
    for( const rapidxml::xml_node<char> *measurement = xml_first_node_nso( document_node,"Measurement",xmlns);
        measurement;
        measurement = XML_NEXT_TWIN(measurement) )
    {
      const rapidxml::xml_attribute<char> *uuid_attrib = measurement->first_attribute( "UUID", 4 );
      if( uuid_attrib && uuid_attrib->value_size() )
      {
        const string thisuuid = xml_value_str( uuid_attrib );
        
        if( uuid_.empty() )
          uuid_ = thisuuid;
        else if( uuid_.length() < 32 )
          uuid_ += (uuid_.empty() ? "" : " ") + thisuuid;
      }//if( uuid_attrib && uuid_attrib->value_size() )
      
      for( const rapidxml::xml_node<char> *remark = xml_first_node_nso(measurement, "Remark", xmlns );
          remark;
          remark = XML_NEXT_TWIN(remark) )
      {
        string str = xml_value_str(remark);
        trim( str );
        if( str.size() )
          remarks_.push_back( str );
      }//for( loop over remarks )
      

      rapidxml::xml_attribute<char> *inspection_attrib = XML_FIRST_ATTRIB( measurement, "dndons:Inspection" );
      inspection_ = xml_value_str( inspection_attrib );
      
      vector<MeasurementShrdPtr> added_measurements;
      
      //<SpectrumMeasurement> nodes should be under <DetectorMeasurement> nodes,
      //  however "Avid Annotated Spectrum" files (see ref40568MWYJS)
      //  dont obey this.  I didnt test this addition very well.
      if( measurement->first_node( "SpectrumMeasurement", 19 ) )
      {
        vector< MeasurementShrdPtr > measurements_this_node;
        
        
        for( const rapidxml::xml_node<char> *spec_meas_node = xml_first_node_nso( measurement, "SpectrumMeasurement", xmlns );
             spec_meas_node;
             spec_meas_node = XML_NEXT_TWIN(spec_meas_node) )
        {
          //each detector will have a <Spectrum> node (for each timeslice)
          for( const rapidxml::xml_node<char> *spectrum = xml_first_node_nso( spec_meas_node, "Spectrum", xmlns );
               spectrum;
               spectrum = XML_NEXT_TWIN(spectrum) )
          {
            std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
            measurements_this_node.push_back( meas );
            workerpool.post( SpectrumNodeDecodeWorker( spectrum, &queue_mutex,
                                                       meas, analysis_info, NULL,
                                              document_node ) );
          }//for( rapidxml::xml_node<char> *spectrum = ... )
        }//for( const rapidxml::xml_node<char> *spec_meas_node = ... )
        
        workerpool.join();

        
        SpectrumNodeDecodeWorker::filter_valid_measurments( measurements_this_node );
        
        for( auto &meas : measurements_this_node )
        {
          measurements_.push_back( meas );
          added_measurements.push_back( meas );
        }//for( auto &meas : measurements_this_node )
      }//if( measurement->first_node( "SpectrumMeasurement", 19 ) )
      
      
      
      //each timeslice will have a <DetectorData> node
      for( const rapidxml::xml_node<char> *det_data_node = xml_first_node_nso( measurement, "DetectorData", xmlns );
           det_data_node;
           det_data_node = XML_NEXT_TWIN(det_data_node) )
      {
        vector< MeasurementShrdPtr > measurements_this_node, gross_counts_this_node;
        
        //For files that give multiple spectra, but with different energy
        //  calibration, for the same data, lets keep track of the <Spectrum>
        //  node each Measurement cooresponds to, so we can append the detector
        //  name with the calibration name (kinda a hack, see
        //  #energy_cal_variants and #keep_energy_cal_variant.
        vector< pair<std::shared_ptr<Measurement>, const rapidxml::xml_node<char> *> > meas_to_spec_node;
        
        //Raytheon portals do this wierd thing of putting two calibrations into
        //  one ChannelData node.. lets find and fix this.  The pointer (second)
        //  will be to the second set of calibration coefficients.
        vector< pair<std::shared_ptr<Measurement>, const rapidxml::xml_node<char> *> > multiple_cals;
        
        
        //Only N42 files that I saw in refMK252OLFFE use this next loop
        //  to initiate the SpectrumNodeDecodeWorkers
        for( const rapidxml::xml_node<char> *spectrum = xml_first_node_nso( det_data_node, "Spectrum", xmlns );
             spectrum;
             spectrum = XML_NEXT_TWIN(spectrum) )
        {
          std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
          measurements_this_node.push_back( meas );
          
          workerpool.post( SpectrumNodeDecodeWorker( spectrum, &queue_mutex,
                            meas, analysis_info, (rapidxml::xml_node<char> *)0,
                                  document_node ) );
          
          if( spectrum->first_attribute( "CalibrationIDs", 14 ) )
            meas_to_spec_node.push_back( make_pair( meas, spectrum ) );
        }//for( loop over spectrums )
        
        
        //Most N42 files use the following loop
        for( const rapidxml::xml_node<char> *det_meas_node = xml_first_node_nso( det_data_node, "DetectorMeasurement", xmlns );
            det_meas_node;
            det_meas_node = XML_NEXT_TWIN(det_meas_node) )
        {
          //Look for <SpectrumMeasurement> measurments
          for( const rapidxml::xml_node<char> *spec_meas_node = xml_first_node_nso( det_meas_node, "SpectrumMeasurement", xmlns );
              spec_meas_node;
              spec_meas_node = XML_NEXT_TWIN(spec_meas_node) )
          {
            //each detector will have a <Spectrum> node (for each timeslice)
            
            for( const rapidxml::xml_node<char> *spectrum = xml_first_node_nso( spec_meas_node, "Spectrum", xmlns );
                spectrum;
                spectrum = XML_NEXT_TWIN(spectrum) )
            {
              std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
              measurements_this_node.push_back( meas );
              
              workerpool.post( SpectrumNodeDecodeWorker( spectrum, &queue_mutex,
                                              meas, analysis_info, det_meas_node, document_node ) );
              
              if( XML_FIRST_ATTRIB(spectrum, "CalibrationIDs") )
                meas_to_spec_node.push_back( make_pair( meas, spectrum ) );
              
              //Raytheon portals do a funky thing where they have both a 2.5MeV and a
              //  8 MeV scale - providing two sets of calibration coefficients under
              //  the SpectrumMeasurement->Spectrum->Calibration->Equation->Coefficients tag, and then placing the
              //  spectra, one after the other, in the <ChannelData>
              {//begin weird Ratheyon check
                const rapidxml::xml_node<char> *firstcal = XML_FIRST_NODE_CHECKED(spectrum, "Calibration");
                const rapidxml::xml_node<char> *firsteqn = XML_FIRST_NODE_CHECKED(firstcal, "Equation");
                const rapidxml::xml_node<char> *firstcoef = XML_FIRST_NODE_CHECKED(firsteqn	, "Coefficients");
                const rapidxml::xml_node<char> *secondcoef = XML_NEXT_TWIN_CHECKED(firstcoef);
                if( secondcoef && secondcoef->value_size() )
                  multiple_cals.push_back( std::make_pair(meas,secondcoef) );
              }//end weird Ratheyon check
            }//for( rapidxml::xml_node<char> *spectrum = ... )
          }//for( const rapidxml::xml_node<char> *spec_meas_node ... )
          
          
          //now look for <GrossCountMeasurement> measurments
          for( const rapidxml::xml_node<char> *gross_count_meas = xml_first_node_nso( det_meas_node, "GrossCountMeasurement", xmlns );
              gross_count_meas;
              (gross_count_meas = XML_NEXT_TWIN(gross_count_meas)))
          {
            std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
            gross_counts_this_node.push_back( meas );
            workerpool.post( GrossCountNodeDecodeWorker( gross_count_meas, meas.get() ) );
          }//for( loop over CountDoseData nodes, dos )
          
      
          //This next section is probably not nesassary
          rapidxml::xml_attribute<char> *detector_attrib = det_meas_node->first_attribute( "Detector", 8 );
          if( detector_attrib && XML_VALUE_ICOMPARE(detector_attrib, "ORTEC Portal") )
          {
            detector_type_ = kOrtecIDMPortalDetector;
          }
        }//for( loop over <DetectorMeasurement> nodes under current <DetectorData> )
        
        workerpool.join();
        
        //Avid system in ref67CSUPJ531 has one Gamma, and one Nuetron
        //  DetectorMeasurment node for each <DetectorData> node, where the
        //  neutron detector isnt explicitly named, but it should clearly go
        //  with the gamma detector
        if( measurements_this_node.size() == gross_counts_this_node.size() )
        {
          for( size_t i = 0; i < measurements_this_node.size(); ++i )
          {
            std::shared_ptr<Measurement> &lhs = measurements_this_node[i];
            std::shared_ptr<Measurement> &rhs = gross_counts_this_node[i];
            if( rhs->detector_name_.empty()
                || lhs->detector_name_ == rhs->detector_name_ )
            {
              lhs->neutron_counts_     = rhs->neutron_counts_;
              lhs->contained_neutron_  = rhs->contained_neutron_;
              lhs->neutron_counts_sum_ = rhs->neutron_counts_sum_;
              rhs.reset();
            }
          }//for( size_t i = 0; i < measurements_this_node.size(); ++i )
        }//if( measurements_this_node.size() == gross_counts_this_node.size() )
        
        
        SpectrumNodeDecodeWorker::filter_valid_measurments( measurements_this_node );
        

        //Now go through and look for Measurments with neutron data, but no
        //  gamma data, and see if can combine with onces that have gamma data...
        bool combined_gam_neut_spec = false;
        for( auto &neut : measurements_this_node )
        {
          if( !neut || !neut->contained_neutron_ || (neut->gamma_counts_ && neut->gamma_counts_->size()) )
            continue;
          
          for( const auto &gam : measurements_this_node )
          {
            if( !gam || gam == neut || gam->contained_neutron_ )
              continue;
            
            string gamdetname = gam->detector_name_;
            string::size_type intercal_pos = gamdetname.find("_intercal_");
            if( intercal_pos != string::npos )
              gamdetname = gamdetname.substr(0,intercal_pos);
            
            bool gamma_matches_neutron = (gamdetname == neut->detector_name_);
            if( !gamma_matches_neutron )
              gamma_matches_neutron = ((gamdetname+"N") == neut->detector_name_);
            
            //Some detectors will have names like "ChannelAGamma", "ChannelANeutron", "ChannelBGamma", "ChannelBNeutron", etc
            if( !gamma_matches_neutron && UtilityFunctions::icontains(gamdetname, "Gamma")
                && UtilityFunctions::icontains(neut->detector_name_, "Neutron") )
            {
              //Note!: this next line alters gamdetname, so this should be the last test!
              UtilityFunctions::ireplace_all(gamdetname, "Gamma", "Neutron" );
              gamma_matches_neutron = UtilityFunctions::iequals(gamdetname, neut->detector_name_);
            }
            
            if( !gamma_matches_neutron )
              continue;

            //A couple basic checks
            
            //Note the real time for RadSeekers can easily be 0.45 seconds off
            if( neut->real_time_ > 0.0f && fabs(neut->real_time_ - gam->real_time_) > 1.0 )
              continue;
            
            //if( neut->live_time_ > 0.0f && fabs(neut->live_time_ - gam->live_time_) > 0.0001 )
              //continue;
            
            if( !neut->start_time_.is_special() && !gam->start_time_.is_special()
               && neut->start_time_ != gam->start_time_ )
              continue;
            
            combined_gam_neut_spec = true;
            gam->neutron_counts_     = neut->neutron_counts_;
            gam->contained_neutron_  = neut->contained_neutron_;
            gam->neutron_counts_sum_ = neut->neutron_counts_sum_;
            
            //We should have a proper neutron real/live time fields, but until I add that, lets hack it.
            if( neut->real_time_ > 0.0f )
            {
              char buffer[128];
              snprintf( buffer, sizeof(buffer), "Neutron Real Time: %.5f s", neut->real_time_ );
              gam->remarks_.push_back( buffer );
            }//
            
            if( neut->live_time_ > 0.0f )
            {
              char buffer[128];
              snprintf( buffer, sizeof(buffer), "Neutron Live Time: %.5f s", neut->live_time_ );
              gam->remarks_.push_back( buffer );
            }//
            
            neut.reset();
            break;
          }//for( const auto &gam : measurements_this_node )
        }//for( const auto &neut : measurements_this_node )
        
        if( combined_gam_neut_spec )
        {
          vector< MeasurementShrdPtr >::iterator begin, end, newend;
          begin = measurements_this_node.begin();
          end = measurements_this_node.end();
          newend = std::remove( begin, end, MeasurementShrdPtr() );
          measurements_this_node.resize( newend - begin );
        }//if( combined_gam_neut_spec )
      
        
        if( gross_counts_this_node.size() )
        {
          //The Measurments in gross_counts might be just neutrons from
          //  detectors with the same name as a gamma detector, if
          //  so we need to combine measurments so ensure_unique_sample_numbers()
          //  wont go wonky
          std::lock_guard<std::mutex> lock( queue_mutex );
          
          for( auto &gross : gross_counts_this_node )
          {
            if( !gross
               || !gross->contained_neutron_
               || gross->gamma_count_sum_ > 0.000001 )
              continue;
            
            bool gross_used = false;
            for( auto &spec : measurements_this_node )
            {
              
              bool gamma_matches_neutron = (spec->detector_name_ == gross->detector_name_);
              if( !gamma_matches_neutron )
                gamma_matches_neutron = ((spec->detector_name_+"N") == gross->detector_name_);
              
              //Some detectors will have names like "ChannelAGamma", "ChannelANeutron", "ChannelBGamma", "ChannelBNeutron", etc
              if( !gamma_matches_neutron && UtilityFunctions::icontains(spec->detector_name_, "Gamma")
                 && UtilityFunctions::icontains(gross->detector_name_, "Neutron") )
              {
                string gamdetname = spec->detector_name_;
                UtilityFunctions::ireplace_all(gamdetname, "Gamma", "Neutron" );
                gamma_matches_neutron = UtilityFunctions::iequals(gamdetname, gross->detector_name_);
                
                //TODO: should probably get rid of "Gamma" in the detector name
                //  now that we have combined them but I havent evaluated the
                //  effect of this yet...
                //if( gamma_matches_neutron )
                //  UtilityFunctions::ireplace_all( spec->detector_name_, "Gamma" );
              }
              
              if( !gamma_matches_neutron )
                continue;
              
              if( spec->contained_neutron_
                  && (spec->neutron_counts_ != gross->neutron_counts_) )
              {
                cerr << SRC_LOCATION << "\n\tWarning: confusing gross count"
                << " situation" << endl;
                continue;
              }//if( spec->neutron_counts_sum_ > 0.000001 )
              
              spec->neutron_counts_     = gross->neutron_counts_;
              spec->contained_neutron_  = gross->contained_neutron_;
              spec->neutron_counts_sum_ = gross->neutron_counts_sum_;
              gross_used = true;
              
              //Gamma data may be represented by multiple spectra, so we have to
              //  loop over all measurements_this_node and not just assume one
              //  match, so dont break once we found a match
            }//for( auto &spectru : measurements_this_node )
            
            if( gross_used )
              gross.reset();  //we dont want this measurment anymore
          }//for( auto &gross : gross_count_meas )
          
          for( const auto &gross : gross_counts_this_node )
          {
            if( gross )
              measurements_this_node.push_back( gross );
          }//for( const auto &gross : gross_counts_this_node )
        }//if( gross_counts_this_node.size() )
        
        Measurement::set_n42_2006_detector_data_node_info( det_data_node, measurements_this_node );
        
        for( size_t multical_index = 0; multical_index < multiple_cals.size(); ++multical_index )
        {
          MeasurementShrdPtr &meas = multiple_cals[multical_index].first;
          const rapidxml::xml_node<char> *second_coefs = multiple_cals[multical_index].second;
          
          if( find( measurements_this_node.begin(), measurements_this_node.end(), meas ) == measurements_this_node.end() )
            continue;
          
          //Sanity check to make sure doing this matches the example data I have
          //  where this bizare fix is needed.
          if( !meas->gamma_counts_ || meas->gamma_counts_->size()!=2048
             || meas->energy_calibration_model_!=Measurement::Polynomial )
            continue;
          
          if( meas->channel_energies_ && meas->channel_energies_->size() )
            meas->channel_energies_.reset();
          
          std::shared_ptr<const vector<float> > oldcounts = meas->gamma_counts_;
          std::shared_ptr<vector<float> > lowerbins
          = std::make_shared<vector<float> >( oldcounts->begin(), oldcounts->begin()+1024 );
          std::shared_ptr<vector<float> > upperbins
          = std::make_shared<vector<float> >( oldcounts->begin()+1024, oldcounts->end() );
          
          MeasurementShrdPtr newmeas = std::make_shared<Measurement>(*meas);
          
          meas->gamma_counts_ = lowerbins;
          meas->gamma_count_sum_ = 0;
          for( const float a : *lowerbins )
            meas->gamma_count_sum_ += a;
          
          newmeas->gamma_counts_ = upperbins;
          newmeas->gamma_count_sum_ = 0;
          for( const float a : *upperbins )
            newmeas->gamma_count_sum_ += a;
          
          vector<float> coeffs;
          if( UtilityFunctions::split_to_floats( second_coefs->value(),
                                                second_coefs->value_size(), coeffs )
             && coeffs != meas->calibration_coeffs_ )
          {
            newmeas->calibration_coeffs_ = coeffs;
            
            //Should come up with better way to decide which calibrarion is which
            //  or make more general or something
            meas->detector_name_ += "_intercal_9MeV";
            newmeas->detector_name_ += "_intercal_2.5MeV";
            
            measurements_this_node.push_back( newmeas );
          }else
          {
#if(PERFORM_DEVELOPER_CHECKS)
            log_developer_error( BOOST_CURRENT_FUNCTION, "Failed to split second energy calibration coefficents into floats" );
#endif
          }//if( second_coefs contained arrray of floats )
        }//for( multiple_cals )
        
        //Make sure if any of the spectrum had the <SourceType> tag, but some
        //  others didnt, we propogate this info to them.  This notably effects
        //  rad assist detectors
        Measurement::SourceType sourcetype = Measurement::UnknownSourceType;
        for( auto &m : measurements_this_node )
        {
          if( !m ) continue;
          if( sourcetype == Measurement::UnknownSourceType )
            sourcetype = m->source_type_;
          else if( m->source_type_ != Measurement::UnknownSourceType )
            sourcetype = max( sourcetype, m->source_type_ );
        }//for( auto &m : measurements_this_node )
        
        for( auto &m : measurements_this_node )
        {
          if( m && (m->source_type_ == Measurement::UnknownSourceType) )
            m->source_type_ = sourcetype;
        }//for( auto &m : measurements_this_node )
        
       //Look for multiple spectra representing the same data, but that actually
       // have different calibrations.
       // See comments for #energy_cal_variants and #keep_energy_cal_variant.
        const vector< MeasurementShrdPtr >::const_iterator datastart = measurements_this_node.begin();
        const vector< MeasurementShrdPtr >::const_iterator dataend = measurements_this_node.end();
        
        for( size_t i = 1; i < meas_to_spec_node.size(); ++i )
        {
          std::shared_ptr<Measurement> &meas = meas_to_spec_node[i].first;
          
          //This check to make sure the Measurement hasnt been discarded is
          //  probably most costly than worth it.
          if( std::find(datastart,dataend,meas) == dataend )
            continue;
          
          vector< pair<std::shared_ptr<Measurement>, const rapidxml::xml_node<char> *> > samenames;
          for( size_t j = 0; j < i; ++j )
          {
            std::shared_ptr<Measurement> &innermeas = meas_to_spec_node[j].first;
            
            if( innermeas->detector_name_ == meas->detector_name_
                && innermeas->start_time_ == meas->start_time_
                && fabs(innermeas->real_time_ - meas->real_time_) < 0.001
               && fabs(innermeas->live_time_ - meas->live_time_) < 0.001 )
            {
              samenames.push_back( meas_to_spec_node[j] );
            }
          }//for( size_t j = 0; j < i; ++j )
          
          if( samenames.size() )
          {
            const rapidxml::xml_node<char> *spec = meas_to_spec_node[i].second;
            const rapidxml::xml_attribute<char> *cal = spec->first_attribute( "CalibrationIDs", 14 );
            meas->detector_name_ += "_intercal_" + xml_value_str( cal );
            for( size_t j = 0; j < samenames.size(); ++j )
            {
              spec = samenames[j].second;
              cal = spec->first_attribute( "CalibrationIDs", 14 );
              samenames[j].first->detector_name_ += "_intercal_" + xml_value_str( cal );
            }
          }//if( samenames.size() )
        }//for( size_t i = 1; i < meas_to_spec_node.size(); ++i )
      
        for( auto &meas : measurements_this_node )
        {
          measurements_.push_back( meas );
          added_measurements.push_back( meas );
        }//for( auto &meas : measurements_this_node )
      }//for( loop over <DetectorData> nodes )
    
      
      //I'm a bad person and am allowing ICD2 files to get to here...  I guess I could
      //  claim code reuse or something, but really its laziness to allow parsing of
      //  ICD2 files for someone when I dont have enough time to properly implement
      //  parsing.
      //  I think this one loop is the only specialization of code to allow this
      //  IF were gonna keep this around, should parallelize the parsing, similar to the above loops
      for( const rapidxml::xml_node<char> *icd2_ana_res = xml_first_node_nso( measurement, "AnalysisResults", "dndoarns:" );
          icd2_ana_res;
          icd2_ana_res = XML_NEXT_TWIN(icd2_ana_res) )
      {
        for( const rapidxml::xml_node<char> *gamma_data = xml_first_node_nso( icd2_ana_res, "AnalyzedGammaData", "dndoarns:" );
            gamma_data;
            gamma_data = XML_NEXT_TWIN(gamma_data) )
        {
          vector< MeasurementShrdPtr > measurements_this_node;
          vector<const rapidxml::xml_node<char> *> spectrum_nodes;
          
          const rapidxml::xml_node<char> *node = xml_first_node_nso( gamma_data, "BackgroundSpectrum", "dndoarns:" );
          if( node )
          {
            std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
            spectrum_nodes.push_back( node );
            measurements_this_node.push_back( meas );
            workerpool.post( SpectrumNodeDecodeWorker( node, &queue_mutex, meas, analysis_info, NULL, document_node ) );
          }
          
          for( node = xml_first_node_nso( gamma_data, "SpectrumSummed", "dndoarns:" );
              node; node = XML_NEXT_TWIN(node) )
          {
            std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
            spectrum_nodes.push_back( node );
            measurements_this_node.push_back( meas );
            workerpool.post( SpectrumNodeDecodeWorker( node, &queue_mutex, meas, analysis_info, NULL, document_node ) );
          }//for( loop over <SpectrumSummed> )
          
          workerpool.join();
          
          for( size_t i = 0; i < measurements_this_node.size(); ++i )
          {
            //<dndoarns:SpectrumSummed dndoarns:DetectorSubset="Partial" dndoarns:NuclidesIdentified="Th-232-001 Ra-226-002">
            std::shared_ptr<Measurement> meas = measurements_this_node[i];
            vector<std::shared_ptr<Measurement> > measv( 1, meas );
            SpectrumNodeDecodeWorker::filter_valid_measurments( measv );
            
            if( measv.empty() )
              continue;
            
            string name = xml_name_str(spectrum_nodes[i]);
            if( UtilityFunctions::icontains(name, "BackgroundSpectrum" ) )
            {
              meas->title_ += " " + string("Background");
             
              //Check if calibration is valid, and if not, fill it in from the
              //  next spectrum... seems a little sketch but works on all the
              //  example files I have.
              if( Measurement::UnknownEquationType == meas->energy_calibration_model_
                  && meas->calibration_coeffs_.empty()
                 && (i < (measurements_this_node.size()-1))
                 && measurements_this_node[i+1]
                 && !measurements_this_node[i+1]->calibration_coeffs_.empty()
                 && measurements_this_node[i+1]->energy_calibration_model_ != Measurement::UnknownEquationType )
              {
                meas->energy_calibration_model_ = measurements_this_node[i+1]->energy_calibration_model_;
                meas->calibration_coeffs_ = measurements_this_node[i+1]->calibration_coeffs_;
                meas->deviation_pairs_ = measurements_this_node[i+1]->deviation_pairs_;
              }
            }//if( UtilityFunctions::icontains(name, "BackgroundSpectrum" ) )
              
            const rapidxml::xml_attribute<char> *nuclides_att = XML_FIRST_ATTRIB( spectrum_nodes[i], "dndoarns:NuclidesIdentified" );
            const string nucstr = xml_value_str( nuclides_att );  //ex "Th-232-001 Ra-226-002"
            //vector<string> nuclides;
            //UtilityFunctions::split( nuclides, nucstr, " \t\n," );
            if( nucstr.size() )
              meas->title_ +=  " Nuclides Reported: " + nucstr + ".";
              
            string detectors;
            map<string,string> det_to_sequence;
            for( const rapidxml::xml_node<char> *subset = xml_first_node_nso( spectrum_nodes[i], "SubsetSampleList", "dndoarns:" );
                subset;
                subset = XML_NEXT_TWIN(subset) )
            {
              //<dndoarns:SubsetSampleList Detector="Aa4">299 300 301 302 303 304 305</dndoarns:SubsetSampleList>
              const rapidxml::xml_attribute<char> *det_att = XML_FIRST_ATTRIB(subset, "Detector");
              const string detname = xml_value_str(det_att);
                
              if( detname.size() )
                detectors += (detectors.size() ? ", " : "") + detname;
                
              const string value = xml_value_str( subset );
              if( subset && subset->value_size() )
              {
                vector<int> samples;
                if( UtilityFunctions::split_to_ints( subset->value(), subset->value_size(), samples ) )
                {
                  const set<int> samplesset( samples.begin(), samples.end() );
                  const string sequence = UtilityFunctions::sequencesToBriefString( samplesset );
                  if( sequence.size() )
                    det_to_sequence[detname] = sequence;
                }else
                {
                  string value = xml_value_str( subset );
                  if( value.size() )
                    det_to_sequence[detname] = value;
                  
                }
              }//if( subset && subset->value_size() )
            }//for( loop over <SubsetSampleList> nodes )
              
            bool all_same = !det_to_sequence.empty();
            for( map<string,string>::const_iterator i = det_to_sequence.begin();
                 all_same && i != det_to_sequence.end(); ++i )
            {
              all_same = (i->second==det_to_sequence.begin()->second);
            }
              
            if( all_same )
            {
              meas->remarks_.push_back( "SampleNumbers: " + det_to_sequence.begin()->second );
            }else
            {
              for( map<string,string>::const_iterator i = det_to_sequence.begin();
                   i != det_to_sequence.end(); ++i )
                meas->remarks_.push_back( "Detector " + i->first + " SampleNumbers: " + i->second );
            }
              
            if( detectors.size() )
              meas->title_ += " Detectors: " + detectors + ". ";
              
            UtilityFunctions::trim( meas->title_ );
              
            measurements_.push_back( meas );
            added_measurements.push_back( meas );
          }//for( loop over <SubsetSampleList> )
        }//for( loop over <AnalyzedGammaData> )
      }//for( loop over <AnalysisResults> in ICD2 file )
      
      
      
      
      const rapidxml::xml_node<char> *info_node = xml_first_node_nso( measurement, "InstrumentInformation", xmlns );
      set_n42_2006_instrument_info_node_info( info_node );
      
      //Try to find the MeasuredItemInformation node - this is all just barely hacked in
      const rapidxml::xml_node<char> *item_info_node = xml_first_node_nso( measurement, "MeasuredItemInformation", xmlns );
      if( !item_info_node && info_node )
        item_info_node = xml_first_node_nso( info_node, "MeasuredItemInformation", xmlns );
      if( !item_info_node ) //HPRDS
        item_info_node = xml_first_node_nso( measurement, "InstrumentLocation", xmlns );
      
      set_n42_2006_measurment_location_information( item_info_node, added_measurements );
      set_n42_2006_deviation_pair_info( info_node, added_measurements );
      
      const rapidxml::xml_node<char> *analysis_node = xml_first_node_nso( measurement, "AnalysisResults", xmlns );
      
      if( !analysis_node )
      {
        // RadSeeker stores its results information in Event->AnalysisResults->RadiationDataAnalysis
        //  Which has a field <SpectrumID>...<SpectrumID> that should match the
        //  SpectrumID attribute of DetectorData->Spectrum.  Lets try to find this
        //  analysis information
        const rapidxml::xml_node<char> *DetectorData = xml_first_node_nso( measurement, "DetectorData", xmlns );
        const rapidxml::xml_node<char> *Spectrum = xml_first_node_nso( DetectorData, "Spectrum", xmlns );
        const rapidxml::xml_attribute<char> *SpectrumIDAttrib = Spectrum ? XML_FIRST_ATTRIB(Spectrum, "SpectrumID") : 0;
        
        if( SpectrumIDAttrib && SpectrumIDAttrib->value_size() )
        {
          const rapidxml::xml_node<char> *AnalysisResults = xml_first_node_nso( document_node->parent(), "AnalysisResults", xmlns );
          const rapidxml::xml_node<char> *RadiationDataAnalysis = xml_first_node_nso( AnalysisResults, "RadiationDataAnalysis", xmlns );
          for( const rapidxml::xml_node<char> *SpectrumID = xml_first_node_nso( RadiationDataAnalysis, "SpectrumID", xmlns );
              SpectrumID && !analysis_node;
              SpectrumID = XML_NEXT_TWIN(SpectrumID) )
          {
            if( rapidxml::internal::compare(SpectrumID->value(), SpectrumID->value_size(),
                                            SpectrumIDAttrib->value(), SpectrumIDAttrib->value_size(), false) )
              analysis_node = RadiationDataAnalysis;
          }
        }//if( SpectrumID && SpectrumID->value_size() )
      }//if( !analysis_node )
      
      setAnalysisInformation( analysis_node, analysis_info );
      
      
      //THe identiFINDER has its neutron info in a <CountDoseData> node under the <Measurement> node
      for( const rapidxml::xml_node<char> *count_dose_data_node = xml_first_node_nso( measurement, "CountDoseData", xmlns );
          count_dose_data_node;
          count_dose_data_node = XML_NEXT_TWIN(count_dose_data_node) )
      {
        const rapidxml::xml_attribute<char> *dettype = XML_FIRST_ATTRIB(count_dose_data_node, "DetectorType");
        if( !dettype || !XML_VALUE_ICOMPARE(dettype, "Neutron") )
          continue;
        
        const rapidxml::xml_node<char> *starttime_node  = xml_first_node_nso( count_dose_data_node, "StartTime", xmlns );
        const rapidxml::xml_node<char> *realtime_node   = xml_first_node_nso( count_dose_data_node, "SampleRealTime", xmlns );
        const rapidxml::xml_node<char> *counts_node     = xml_first_node_nso( count_dose_data_node, "Counts", xmlns );
        const rapidxml::xml_node<char> *count_rate_node = xml_first_node_nso( count_dose_data_node, "CountRate", xmlns ); //Only seen in RadEagle files so far
        const rapidxml::xml_node<char> *remark_node     = xml_first_node_nso( count_dose_data_node, "Remark", xmlns );
        
        const bool has_counts_node = (counts_node && counts_node->value_size());
        const bool has_count_rate_node = (count_rate_node && count_rate_node->value_size());
        
        if( !starttime_node || !starttime_node->value_size()
            || !realtime_node || !realtime_node->value_size()
           || !(has_counts_node || has_count_rate_node) )
          continue;

        
        if( remark_node )
        {
          //RadEagle files contain three <CountDoseData> elements that
          //  each have a single <remark> element with values of "Minimum",
          //  "Average", or "Maximum".  These all give CountRate (not
          //  dose, or total counts), so lets only take the "Average" one
          //  20171218: I have only seen files with zero neutron counts, so this
          //            is all likely not yet correct (example measurements that
          //            did contain neutrons in SPE file, didnt have them in the
          //            N42 file, so something may be odd).
           if( XML_VALUE_ICOMPARE(remark_node,"Minimum") || XML_VALUE_ICOMPARE(remark_node,"Maximum") )
             continue;
        }//if( remark_node )
        
        float counts;
        if( has_counts_node )
        {
          if( !xml_value_to_flt( counts_node, counts ) )
          continue;
        }else
        {
          if( !xml_value_to_flt( count_rate_node, counts ) )
            continue;
        }
        
        const boost::posix_time::ptime start_time = time_from_string( xml_value_str(starttime_node).c_str() );
        if( start_time.is_special() )
          continue;
        
        const float realtimesec = time_duration_in_seconds( realtime_node->value(), realtime_node->value_size() );
        
        if( !has_counts_node && realtimesec > 0.0 )
          counts *= realtimesec;
        
        for( auto &meas : added_measurements )
        {
          if( !meas || meas->contained_neutron_
              || meas->start_time_ != start_time )
            continue;
          
          if( fabs(realtimesec - meas->real_time_) >= 1.0 )
            continue;
          
          meas->contained_neutron_ = true;
          meas->neutron_counts_.resize( 1 );
          meas->neutron_counts_[0] = counts;
          meas->neutron_counts_sum_ = counts;
          break;
        }
      }//for( loop over <CountDoseData> nodes )
    
    }//for( loop over all measurements )
  }//if( is_cambio ) / else
 
  
  //Some HPRDS files will have Instrument information right bellow the document node.
  //  For example see refR96VUZSYYM
  //  (I dont think checking here will ever overide info found for a specific
  //  spectrum, but not totally sure, so will leave in a warning, which can be taken out after running the tests)
  {
    const rapidxml::xml_node<char> *info_node = xml_first_node_nso( document_node, "InstrumentInformation", xmlns );
    if( info_node )
    {
      if( instrument_type_.size() )
        cerr << "MeasurementInfo::load_2006_N42_from_doc(): may be overwriting InstrumentInformation already gathered from a specific spectrum" << endl;
      set_n42_2006_instrument_info_node_info( info_node );
    }
  }
  
  for( const rapidxml::xml_node<char> *remark = xml_first_node_nso( document_node, "Remark", xmlns );
      remark;
      remark = XML_NEXT_TWIN(remark) )
  {
    string str = xml_value_str( remark );
    trim( str );
    if( str.size() )
      remarks_.push_back( str );
  }//for( loop over remarks )
  
  
  if( analysis_info && analysis_info->results_.size() )
    detectors_analysis_ = analysis_info;
  
  if( measurements_.empty() )
  {
    stringstream msg;
    msg << SRC_LOCATION << "\n\rNo Measurments found inside ICD1/XML file";
    throw runtime_error( msg.str() );
  }//if( measurements_.empty() )
  
  
  //Lets try to figure out if we can fill out detector_type_
  if( iequals( manufacturer_,"ORTEC" ) )
  {
    if( iequals(instrument_model_,"OSASP") )
      detector_type_ = kOrtecIDMPortalDetector;
    else if( icontains(instrument_model_,"100") )
      detector_type_ = kDetectiveEx100Detector;
    else if( icontains(instrument_model_,"Detective-EX") )
      detector_type_ = kDetectiveExDetector;
    else if( icontains(instrument_model_,"Detective") && contains(instrument_model_,"100") )
      detector_type_ = kDetectiveEx100Detector;
    else if( icontains(instrument_model_,"Detective") && icontains(instrument_model_,"micro") )
      detector_type_ = kMicroDetectiveDetector;
    else if( icontains(instrument_model_,"Detective") )
      detector_type_ = kDetectiveDetector;
  }else if( iequals(instrument_type_,"PVT Portal")
           && iequals(manufacturer_,"SAIC") )
  {
    detector_type_ = kSAIC8Detector;
  }else if( icontains(instrument_model_,"identiFINDER")
           //&& icontains(manufacturer_,"FLIR")
           )
  {
    //WE could probably get more specific detector fine tuning here, cause we have:
    //      <InstrumentModel>identiFINDER 2 ULCS-TNG</InstrumentModel>
    //      <InstrumentVersion>Hardware: 4C	Firmware: 5.00.54	Operating System: 1.2.040	Application: 2.37</InstrumentVersion>
    
    if( icontains(instrument_model_,"LG") )
      detector_type_ = kIdentiFinderLaBr3Detector;
    else
      detector_type_ = kIdentiFinderNGDetector;
  }else if( icontains(manufacturer_,"FLIR") || icontains(instrument_model_,"Interceptor") )
  {
    
  }else if( icontains(instrument_model_,"SAM940") || icontains(instrument_model_,"SAM 940") || icontains(instrument_model_,"SAM Eagle") )
  {
    if( icontains(instrument_model_,"LaBr") )
      detector_type_ = kSam940LaBr3;
    else
      detector_type_ = kSam940;
  }else if( istarts_with(instrument_model_,"RE ") || icontains(instrument_model_,"RadEagle") || icontains(instrument_model_,"Rad Eagle" ) )
  {
    if( !manufacturer_.empty() && !icontains(manufacturer_, "ortec") )
      manufacturer_ += " (Ortec)";
    else if( !icontains(manufacturer_, "ortec") )
      manufacturer_ = "Ortec";
    //set_detector_type_from_other_info() will set detector_type_ ...
  }else if( icontains(instrument_model_,"SAM") && icontains(instrument_model_,"945") )
  {
    detector_type_ = kSam945;
  }else if( (icontains(manufacturer_,"ICx Radiation") || icontains(manufacturer_,"FLIR"))
            && icontains(instrument_model_,"Raider") )
  {
    detector_type_ = kMicroRaiderDetector;
  }else if( icontains(manufacturer_,"Canberra Industries, Inc.") )
  {
    //Check to see if detectors like "Aa1N+Aa2N", or "Aa1N+Aa2N+Ba1N+Ba2N+Ca1N+Ca2N+Da1N+Da2N"
    //  exist and if its made up of other detectors, get rid of those spectra
    
  }else if( icontains(instrument_type_,"SpecPortal")
            && icontains(manufacturer_,"SSC Pacific")
            && icontains(instrument_model_,"MPS Pod") )
  {
    //Gamma spectrum is in CPS, so multiply each spectrum by live time.
    //  Note that there is no indication of this in the file, other than
    //  fracitonal counts
    for( auto &m : measurements_ )
    {
      if( !m || (m->live_time_ < 1.0f) )  //1 second is arbitrary
        continue;
      
      //We could probably add some count rate sanity check here too.
      
      if( m->contained_neutron_ )
      {
        for( auto &f : m->neutron_counts_ )
          f *= m->live_time_;
        m->neutron_counts_sum_ *= m->live_time_;
      }
      
      if( m->gamma_counts_ )
      {
        for( auto &f : const_cast<vector<float>&>(*m->gamma_counts_) )  //We know this is safe right here.
          f *= m->live_time_;
        m->gamma_count_sum_ *= m->live_time_;
      }
      if( m->gamma_counts_ || m->contained_neutron_ )
        m->remarks_.push_back( "Gamma/Neutron counts have been mutliplied by live time, "
                               "to account for observed shortcommings of this detectors N42-2006 format." );
    }//for( auto &m : measurements_ )
  }else if( manufacturer_.size() || instrument_model_.size() )
  {
//    if( (manufacturer_=="ICx Technologies" && instrument_model_=="identiFINDER") )
//    {
//    }
    
    if( !(manufacturer_=="Princeton Gamma-Tech Instruments, Inc." && instrument_model_=="RIIDEye")
        && !(manufacturer_=="ICx Technologies" && instrument_model_=="")
        && !(manufacturer_=="Radiation Solutions Inc." /* && instrument_model_=="RS-701"*/)
        && !(manufacturer_=="Raytheon" && instrument_model_=="Variant L" )
       && !(manufacturer_=="Avid Annotated Spectrum" /* && instrument_model_==""*/)
       && !(manufacturer_=="Mirion Technologies" && instrument_model_=="model Pedestrian G")
       && !(manufacturer_=="Princeton Gamma-Tech Instruments, Inc." && instrument_model_=="")
       && !(manufacturer_=="Nucsafe" && instrument_model_=="G4_Predator")
       && !(manufacturer_=="Princeton Gamma-Tech Instruments, Inc." && instrument_model_=="Model 135")
       && !(manufacturer_=="" && instrument_model_=="Self-Occuluding Quad NaI Configuration")
       && !(manufacturer_=="" && instrument_model_=="3x3x12 inch NaI Side Ortec Digibase MCA")
       && !(manufacturer_=="Berkeley Nucleonics Corp." && instrument_model_=="SAM 945")
       && !(manufacturer_=="Canberra Industries, Inc." && instrument_model_=="ASP EDM")
       && !(manufacturer_=="Smiths Detection" && instrument_model_=="RadSeeker_DL")
       && !(manufacturer_=="" && instrument_model_=="")
       )
      cerr << "Unknown detector type: maufacturer=" << manufacturer_ << ", ins_model=" << instrument_model_ << endl;
    
    //Unknown detector type: maufacturer=Smiths Detection, ins_model=RadSeeker_CS
    //Unknown detector type: maufacturer=Smiths Detection, ins_model=RadSeeker_DL
    //Unknown detector type: maufacturer=RSL/SPAWAR, ins_model=Multipod MPS
  }
  
  cleanup_after_load();
}//bool load_2006_N42_from_doc( rapidxml::xml_node<char> *document_node )


void Measurement::popuplate_channel_energies_from_coeffs()
{
  if( !gamma_counts_ || gamma_counts_->empty() )
    throw runtime_error( "popuplate_channel_energies_from_coeffs():"
                         " no gamma counts" );
  
  if( calibration_coeffs_.size() < 2 )
    throw runtime_error( "popuplate_channel_energies_from_coeffs():"
                         " no calibration coefficients" );
  
  if( !!channel_energies_ && !channel_energies_->empty() )
    throw runtime_error( "popuplate_channel_energies_from_coeffs():"
                         " channel energies already populated" );
  
  const size_t nbin = gamma_counts_->size();
    
  switch( energy_calibration_model_ )
  {
    case Measurement::Polynomial:
      channel_energies_ = polynomial_binning( calibration_coeffs_, 
                                              nbin, deviation_pairs_ );
      break;
        
    case Measurement::FullRangeFraction:
      channel_energies_ = fullrangefraction_binning( calibration_coeffs_, 
                                                    nbin, deviation_pairs_ );
      break;
        
    case Measurement::LowerChannelEdge:
      channel_energies_.reset( new vector<float>(calibration_coeffs_) );
      break;
        
    case Measurement::UnknownEquationType:
      throw runtime_error( "popuplate_channel_energies_from_coeffs():"
                           " unknown equation type" );
      break;
  }//switch( meas->energy_calibration_model_ )
}//void popuplate_channel_energies_from_coeffs()


void Measurement::add_calibration_to_2012_N42_xml(
                                rapidxml::xml_node<char> *RadInstrumentData,
                                std::mutex &xmldocmutex,
                                const int i ) const
{
  using namespace rapidxml;
  const char *val = (const char *)0;
  char buffer[32];
  
  rapidxml::xml_document<char> *doc = RadInstrumentData->document();
  
  const char *coefname = 0;
  xml_node<char> *EnergyCalibration = 0, *node = 0;
  
  stringstream valuestrm;
  std::vector<float> coefs = calibration_coeffs();
  const size_t nbin = gamma_counts()->size();
  
  switch( energy_calibration_model() )
  {
    case Measurement::FullRangeFraction:
      coefs = fullrangefraction_coef_to_polynomial( coefs, nbin );
      //note intential fallthrough
    case Measurement::Polynomial:
    {
      coefname = "CoefficientValues";
      
      //Technically this node must have exactly three coeficients, but here we'll
      //  slip in some more if we have them...
      const size_t ncoef = std::max( size_t(3), coefs.size() );
      for( size_t j = 0; j < ncoef; ++j )
      {
        snprintf( buffer, sizeof(buffer), "%s%.9g", (j?" ":""), (j>=coefs.size() ? 0.0f : coefs[j]) );
        valuestrm << buffer;
      }
     
      break;
    }//case Measurement::Polynomial:
      
    case Measurement::LowerChannelEdge:
    case Measurement::UnknownEquationType:  //Taking a guess here (we should alway capture FRF, Poly, and LowrBinEnergy correctly,
      if( !!channel_energies_ || calibration_coeffs().size() )
      {
        coefname = "EnergyBoundaryValues";
        const std::vector<float> &b = (!!channel_energies_ ? *channel_energies_ : calibration_coeffs());
        
        //This next part should be formatted better!
        for( size_t j = 0; j < b.size(); ++j )
          valuestrm << (j?" ":"") << b[j];
        
        //According to N42 2012 standard, we must give energy of upper channel
        if( b.size() && b.size()<=nbin )
          valuestrm << " " << ((2.0f*b[b.size()-1])-b[b.size()-2]);
      }//if( we have some sort of coefficients )
    break;
  }//switch( energy_calibration_model() )
  
  
  if( coefname )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    EnergyCalibration = doc->allocate_node( node_element, "EnergyCalibration" );
    RadInstrumentData->append_node( EnergyCalibration );
    
    snprintf( buffer, sizeof(buffer), "EnergyCal%i", int(i) );
    val = doc->allocate_string( buffer );
    EnergyCalibration->append_attribute( doc->allocate_attribute( "id", val ) );
    
    val = doc->allocate_string( valuestrm.str().c_str() );
    node = doc->allocate_node( node_element, coefname, val );
    EnergyCalibration->append_node( node );
  }//if( coefname )
  
  const DeviationPairVec &devpairs = deviation_pairs();
  if( devpairs.size() )
  {
    stringstream EnergyValuesStrm, EnergyDeviationValuesStrm;
    for( size_t j = 0; j < devpairs.size(); ++j )
    {
      snprintf( buffer, sizeof(buffer), "%s%.9g", (j?" ":""), devpairs[j].first );
      EnergyValuesStrm << buffer;
      snprintf( buffer, sizeof(buffer), "%s%.9g", (j?" ":""), devpairs[j].second );
      EnergyDeviationValuesStrm << buffer;
    }
    
    std::lock_guard<std::mutex> lock( xmldocmutex );
    if( !EnergyCalibration )
    {
      EnergyCalibration = doc->allocate_node( node_element, "EnergyCalibration" );
      RadInstrumentData->append_node( EnergyCalibration );
    }
    
    val = doc->allocate_string( EnergyValuesStrm.str().c_str() );
    node = doc->allocate_node( node_element, "EnergyValues", val );
    EnergyCalibration->append_node( node );
    
    val = doc->allocate_string( EnergyDeviationValuesStrm.str().c_str() );
    node = doc->allocate_node( node_element, "EnergyDeviationValues", val );
    EnergyCalibration->append_node( node );
  }//if( devpairs.size() )
}//void add_calibration_to_2012_N42_xml(...)


//getCalibrationToSpectrumMap(...): builds map from the binning shared vector to
//  the the index of a MeasurementShrdPtr that has that binning
namespace
{
typedef std::map< std::shared_ptr< const std::vector<float> >, size_t > BinningToIndexMap;

void insert_N42_calibration_nodes( const std::vector< MeasurementShrdPtr > &measurements,
                                ::rapidxml::xml_node<char> *RadInstrumentData,
                                std::mutex &xmldocmutex,
                                BinningToIndexMap &calToSpecMap )
{
  using namespace rapidxml;
  calToSpecMap.clear();
  
  for( size_t i = 0; i < measurements.size(); ++i )
  {
    const MeasurementShrdPtr &meas = measurements[i];
    
    if( !meas || !meas->gamma_counts() || meas->gamma_counts()->empty() )
      continue;
    
    std::shared_ptr< const std::vector<float> > binning = meas->channel_energies();
    BinningToIndexMap::const_iterator iter = calToSpecMap.find( binning );
    if( iter == calToSpecMap.end() )
    {
      calToSpecMap[binning] = i;
      meas->add_calibration_to_2012_N42_xml( RadInstrumentData, xmldocmutex, int(i) );
    }//if( iter == calToSpecMap.end() )
  }//for( size_t i = 0; i < measurements.size(); ++i )
}//void insert_N42_calibration_nodes(...)
}//namespace


const std::string &convert_n42_instrument_type_from_2006_to_2012( const std::string &classcode )
{
  static const string PortalMonitor = "Portal Monitor";
  static const string SpecPortal = "Spectroscopic Portal Monitor";
  static const string RadionuclideIdentifier = "Radionuclide Identifier";
  static const string PersonalRadiationDetector = "Spectroscopic Personal Radiation Detector"; //hmm, prob not best
  static const string SurveyMeter = "Backpack or Personal Radiation Scanner";
  static const string Spectrometer = "Spectroscopic Personal Radiation Detector";
  
  if( iequals(classcode, "PortalMonitor") || iequals(classcode, "PVT Portal") )
    return PortalMonitor;
  else if( iequals(classcode, "SpecPortal") )
    return SpecPortal;
  else if( iequals(classcode, "RadionuclideIdentifier") )
    return RadionuclideIdentifier;
  else if( iequals(classcode, "PersonalRadiationDetector") )
    return PersonalRadiationDetector;
  else if( iequals(classcode, "SurveyMeter") )
    return SurveyMeter;
  else if( iequals(classcode, "Spectrometer") )
    return Spectrometer;
  
  return classcode;
}//const std::string &convert_n42_instrument_type_2006_to_2012( const std::string &input )


void MeasurementInfo::add_analysis_results_to_2012_N42(
                                 const DetectorAnalysis &ana,
                                 ::rapidxml::xml_node<char> *RadInstrumentData,
                                 std::mutex &xmldocmutex )
{
  /* The relevant part of the N42 XSD (https://www.nist.gov/file/42636, or
     https://www.nist.gov/programs-projects/ansiieee-n4242-standard)
   <xsd:complexType name="AnalysisResultsType">
   <xsd:annotation>
   <xsd:documentation>A data type to provide information on the results of a radiation data analysis.</xsd:documentation>
   <xsd:appinfo>
   </xsd:appinfo>
   </xsd:annotation>
   <xsd:complexContent>
   <xsd:extension base="n42:OptIdComplexObjectType">
   <xsd:sequence>
   <xsd:element ref="n42:AnalysisStartDateTime" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisComputationDuration" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisAlgorithmName" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisAlgorithmCreatorName" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisAlgorithmDescription" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisAlgorithmVersion" minOccurs="0" maxOccurs="unbounded"/>
   <xsd:element ref="n42:AnalysisAlgorithmSetting" minOccurs="0" maxOccurs="unbounded"/>
   <xsd:element ref="n42:AnalysisResultStatusCode" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisConfidenceValue" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:AnalysisResultDescription" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:RadAlarm" minOccurs="0" maxOccurs="unbounded"/>
   <xsd:element ref="n42:NuclideAnalysisResults" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:SpectrumPeakAnalysisResults" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:GrossCountAnalysisResults" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:DoseAnalysisResults" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:ExposureAnalysisResults" minOccurs="0" maxOccurs="1"/>
   <xsd:element ref="n42:Fault" minOccurs="0" maxOccurs="unbounded"/>
   <xsd:element ref="n42:AnalysisResultsExtension" minOccurs="0" maxOccurs="unbounded"/>
   </xsd:sequence>
   <xsd:attribute name="radMeasurementGroupReferences" type="xsd:IDREFS" use="optional">
   <xsd:annotation>
   <xsd:documentation>Identifies the RadMeasurementGroup element(s) within the N42 XML document that applies to this particular analysis. There shall be no duplicate IDREF values in the list.</xsd:documentation>
   </xsd:annotation>
   </xsd:attribute>
   <xsd:attribute name="derivedDataReferences" type="xsd:IDREFS" use="optional">
   <xsd:annotation>
   <xsd:documentation>Identifies the DerivedData element(s) within the N42 XML document that applies to this particular analysis. There shall be no duplicate IDREF values in the list.</xsd:documentation>
   </xsd:annotation>
   </xsd:attribute>
   <xsd:attribute name="radMeasurementReferences" type="xsd:IDREFS" use="optional">
   <xsd:annotation>
   <xsd:documentation>Identifies the RadMeasurement element(s) within the N42 XML document that applies to a particular analysis. There shall be no duplicate IDREF values in the list.</xsd:documentation>
   </xsd:annotation>
   </xsd:attribute>
   </xsd:extension>
   </xsd:complexContent>
   </xsd:complexType>
   */
  using namespace rapidxml;
  const char *val = 0;
  char buffer[512];
  xml_attribute<> *attr = 0;
  
  ::rapidxml::xml_document<char> *doc = RadInstrumentData->document();
  
  xml_node<char> *AnalysisResults = 0;
  
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    AnalysisResults = doc->allocate_node( node_element, "AnalysisResults" );
    RadInstrumentData->append_node( AnalysisResults );
  }
  
  for( size_t i = 0; i < ana.remarks_.size(); ++i )
  {
    const string &remark = ana.remarks_[i];
    if( remark.size() )
    {
      std::lock_guard<std::mutex> lock( xmldocmutex );
      val = doc->allocate_string( remark.c_str(), remark.size()+1 );
      xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
      AnalysisResults->append_node( remark );
    }
  }
  
  if( ana.algorithm_name_.size() )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( ana.algorithm_name_.c_str(), ana.algorithm_name_.size()+1 );
    xml_node<char> *name = doc->allocate_node( node_element, "AnalysisAlgorithmName", val );
    AnalysisResults->append_node( name );
  }
  
  if( ana.algorithm_creator_.size() )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( ana.algorithm_creator_.c_str(), ana.algorithm_creator_.size()+1 );
    xml_node<char> *creator = doc->allocate_node( node_element, "AnalysisAlgorithmCreatorName", val );
    AnalysisResults->append_node( creator );
  }
  
  if( ana.algorithm_description_.size() )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( ana.algorithm_description_.c_str(), ana.algorithm_description_.size()+1 );
    xml_node<char> *desc = doc->allocate_node( node_element, "AnalysisAlgorithmDescription", val );
    AnalysisResults->append_node( desc );
  }
  
  if( ana.algorithm_version_.size() )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( ana.algorithm_version_.c_str(), ana.algorithm_version_.size()+1 );
    xml_node<char> *version = doc->allocate_node( node_element, "AnalysisAlgorithmVersion", val );
    AnalysisResults->append_node( version );
  }
  
  //<AnalysisAlgorithmSetting>
  //<AnalysisResultStatusCode>
  //<AnalysisConfidenceValue>
  
  if( ana.algorithm_result_description_.size() )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( ana.algorithm_result_description_.c_str(), ana.algorithm_result_description_.size()+1 );
    xml_node<char> *desc = doc->allocate_node( node_element, "AnalysisResultDescription", val );
    AnalysisResults->append_node( desc );
  }
  
  //<RadAlarm>
  //<NuclideAnalysisResults>
  //<SpectrumPeakAnalysisResults>
  //<GrossCountAnalysisResults>
  //<DoseAnalysisResults>
  //<ExposureAnalysisResults>
  //<Fault>
  

  
  


  
  
  xml_node<char> *result_node = 0;
  
  for( const DetectorAnalysisResult &result : ana.results_ )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    if( result.nuclide_.size() )
    {
      if( !result_node )
      {
        result_node = doc->allocate_node( node_element, "NuclideAnalysisResults" );
        AnalysisResults->append_node( result_node );
      }
      
/*
      xml_node<char> *AnalysisStartDateTime = (xml_node<char> *)0;
      if( !AnalysisStartDateTime && !result.start_time_.is_special() )
      {
        const string timestr = UtilityFunctions::to_iso_string( result.start_time_ );
        val = doc->allocate_string( timestr.c_str(), timestr.size()+1 );
        AnalysisStartDateTime = doc->allocate_node( node_element, "StartTime", val );
        AnalysisResults->append_node( AnalysisStartDateTime );
      }
*/
      
      xml_node<char> *nuclide_node = doc->allocate_node( node_element, "Nuclide" );
      result_node->append_node( nuclide_node );
      
      val = doc->allocate_string( result.nuclide_.c_str(), result.nuclide_.size()+1 );
      xml_node<char> *NuclideName = doc->allocate_node( node_element, "NuclideName", val );
      nuclide_node->append_node( NuclideName );
      
      if( result.remark_.size() )
      {
        val = doc->allocate_string( result.remark_.c_str(), result.remark_.size()+1 );
        xml_node<char> *Remark = doc->allocate_node( node_element, "Remark", val );
        nuclide_node->append_node( Remark );
      }
      
      if( result.activity_ > 0.0f )
      {
        snprintf( buffer, sizeof(buffer), "%1.8E", (result.activity_/1000.0f) );
        val = doc->allocate_string( buffer );
        xml_node<char> *NuclideActivityValue = doc->allocate_node( node_element, "NuclideActivityValue", val );
        attr = doc->allocate_attribute( "units", "kBq" );
        NuclideActivityValue->append_attribute( attr );
        nuclide_node->append_node( NuclideActivityValue );
      }
      
      if( result.nuclide_type_.size() )
      {
        val = doc->allocate_string( result.nuclide_type_.c_str(), result.nuclide_type_.size()+1 );
        xml_node<char> *NuclideType = doc->allocate_node( node_element, "NuclideType", val );
        nuclide_node->append_node( NuclideType );
      }
      
      if( result.id_confidence_.size() )
      {
        val = doc->allocate_string( result.id_confidence_.c_str(), result.id_confidence_.size()+1 );
        xml_node<char> *NuclideIDConfidenceDescription = doc->allocate_node( node_element, "NuclideIDConfidenceIndication", val ); //NuclideIDConfidenceDescription
        nuclide_node->append_node( NuclideIDConfidenceDescription );
      }
      
      xml_node<char> *NuclideExtension = 0;
      
      if( result.real_time_ > 0.0f )
      {
        NuclideExtension = doc->allocate_node( node_element, "NuclideExtension" );
        nuclide_node->append_node( NuclideExtension );
        snprintf( buffer, sizeof(buffer), "PT%fS", result.real_time_ );
        val = doc->allocate_string( buffer );
        xml_node<char> *SampleRealTime = doc->allocate_node( node_element, "SampleRealTime", val );
        NuclideExtension->append_node( SampleRealTime );
      }//if( we should record some more info )
      
      if( result.distance_ > 0.0f )
      {
        xml_node<char> *SourcePosition = doc->allocate_node( node_element, "SourcePosition" );
        nuclide_node->append_node( SourcePosition );
        
        xml_node<char> *RelativeLocation = doc->allocate_node( node_element, "RelativeLocation" );
        SourcePosition->append_node( RelativeLocation );
        
        snprintf( buffer, sizeof(buffer), "%f", (result.distance_/1000.0f) );
        attr = doc->allocate_attribute( "units", "m" );
        val = doc->allocate_string( buffer );
        xml_node<char> *DistanceValue = doc->allocate_node( node_element, "DistanceValue", val );
        RelativeLocation->append_node( DistanceValue );
      }
      
      
      if( result.detector_.size() )
      {
        if( !NuclideExtension )
        {
          NuclideExtension = doc->allocate_node( node_element, "NuclideExtension" );
          nuclide_node->append_node( NuclideExtension );
        }
        
        val = doc->allocate_string( result.detector_.c_str(), result.detector_.size()+1 );
        xml_node<char> *Detector = doc->allocate_node( node_element, "Detector", val );
        NuclideExtension->append_node( Detector );
      }
    }//if( result.nuclide_.size() )
    
    if( result.dose_rate_ > 0.0f )
    {
      xml_node<char> *DoseAnalysisResults = doc->allocate_node( node_element, "DoseAnalysisResults" );
      AnalysisResults->append_node( DoseAnalysisResults );
      
      if( result.remark_.size() )
      {
        val = doc->allocate_string( result.remark_.c_str(), result.remark_.size()+1 );
        xml_node<char> *Remark = doc->allocate_node( node_element, "Remark", val );
        DoseAnalysisResults->append_node( Remark );
      }
  
/*
      if( !AnalysisStartDateTime && !result.start_time_.is_special() )
      {
        const string timestr = UtilityFunctions::to_iso_string( result.start_time_ );
        val = doc->allocate_string( timestr.c_str(), timestr.size()+1 );
        AnalysisStartDateTime = doc->allocate_node( node_element, "StartTime", val );
        AnalysisResults->append_node( AnalysisStartDateTime );
      }
*/
      
      snprintf( buffer, sizeof(buffer), "%1.8E", result.dose_rate_ );
      val = doc->allocate_string( buffer );
      xml_node<char> *AverageDoseRateValue = doc->allocate_node( node_element, "AverageDoseRateValue", val );
      attr = doc->allocate_attribute( "units", "\xc2\xb5Sv/h" ); //\xc2\xb5 --> micro
      AverageDoseRateValue->append_attribute( attr );
      DoseAnalysisResults->append_node( AverageDoseRateValue );
      
      if( result.real_time_ > 0.f )
      {
        snprintf( buffer, sizeof(buffer), "%1.8E", (result.dose_rate_*result.real_time_) );
        val = doc->allocate_string( buffer );
        xml_node<char> *TotalDoseValue = doc->allocate_node( node_element, "TotalDoseValue", val );
        attr = doc->allocate_attribute( "units", "\xc2\xb5Sv" ); //\xc2\xb5 --> micro
        TotalDoseValue->append_attribute( attr );
        DoseAnalysisResults->append_node( TotalDoseValue );
      }
      
      if( result.distance_ > 0.0f )
      {
        xml_node<char> *SourcePosition = doc->allocate_node( node_element, "SourcePosition" );
        DoseAnalysisResults->append_node( SourcePosition );
        
        xml_node<char> *RelativeLocation = doc->allocate_node( node_element, "RelativeLocation" );
        SourcePosition->append_node( RelativeLocation );
        
        snprintf( buffer, sizeof(buffer), "%f", (result.distance_/1000.0f) );
        attr = doc->allocate_attribute( "units", "m" );
        val = doc->allocate_string( buffer );
        xml_node<char> *DistanceValue = doc->allocate_node( node_element, "DistanceValue", val );
        RelativeLocation->append_node( DistanceValue );
        RelativeLocation->append_attribute( attr );
      }
    }//if( result.dose_rate_ > 0.0f )
    
    
    //Children of <Nuclide> that are not checked for:
    //<NuclideIdentifiedIndicator>,
    //<NuclideActivityUncertaintyValue>,
    //<NuclideMinimumDetectableActivityValue>,
    //<NuclideCategoryDescription>,
    //<NuclideSourceGeometryCode>,
    //<NuclideShieldingAtomicNumber>,
    //<NuclideShieldingArealDensityValue>,
    //<NuclideExtension>
  }//for( const DetectorAnalysisResult &result : ana.results_ )
}//void add_analysis_results_to_2012_N42(...)


std::string MeasurementInfo::determine_rad_detector_kind_code() const
{
  string det_kind = "Other";
  switch( detector_type_ )
  {
    case kDetectiveDetector: case kDetectiveExDetector:
    case kDetectiveEx100Detector: case kOrtecIDMPortalDetector:
    case kFalcon5000: case kMicroDetectiveDetector:
      det_kind = "HPGe";
      break;
      
    case kGR135Detector: case kIdentiFinderDetector:
    case kIdentiFinderNGDetector: case kRadHunterNaI:
    case kRsi701: case kRsi705: case kAvidRsi: case kOrtecRadEagleNai:
    case kSam940: case kSam945:
      det_kind = "NaI";
      break;
      
    case kIdentiFinderLaBr3Detector: case kRadHunterLaBr3: case kSam940LaBr3:
    case kOrtecRadEagleLaBr:
      det_kind = "LaBr3";
      break;
      
    case kOrtecRadEagleCeBr2Inch:
    case kOrtecRadEagleCeBr3Inch:
      det_kind = "CeBr3";
      break;
      
    case kSAIC8Detector:
      det_kind = "PVT";
      break;
      
    case kMicroRaiderDetector:
      det_kind = "CZT";
      break;
      
    case kUnknownDetector:
      if( num_gamma_channels() > 4100 )
        det_kind = "HPGe";
      else if( manufacturer_=="Raytheon" && instrument_model_=="Variant L" )
        det_kind = "NaI";
      else if( manufacturer_=="Mirion Technologies" && instrument_model_=="model Pedestrian G" )
        det_kind = "NaI";
      else if( manufacturer_=="Nucsafe" && instrument_model_=="G4_Predator" )
        det_kind = "PVT";
      break;
  }//switch( detector_type_ )
  
  return det_kind;
}//determine_rad_detector_kind_code()


std::shared_ptr< ::rapidxml::xml_document<char> > MeasurementInfo::create_2012_N42_xml() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  using namespace rapidxml;
  
  std::mutex xmldocmutex;
  BinningToIndexMap calToSpecMap;
  SpecUtilsAsync::ThreadPool workerpool;
  
  std::shared_ptr<xml_document<char> > doc = std::make_shared<xml_document<char> >();
  
  const char *val = (const char *)0;
  
  xml_node<char> *declaration = doc->allocate_node( node_declaration );
  doc->append_node( declaration );
  xml_attribute<> *attr = doc->allocate_attribute( "version", "1.0" );
  declaration->append_attribute( attr );
  attr = doc->allocate_attribute( "encoding", "utf-8" );
  declaration->append_attribute( attr );
  
  
  xml_node<char> *RadInstrumentData = doc->allocate_node( node_element, "RadInstrumentData" );
  doc->append_node( RadInstrumentData );
  
  val = doc->allocate_string( uuid_.c_str(), uuid_.size()+1 );
  attr = doc->allocate_attribute( "n42DocUUID", val );
  RadInstrumentData->append_attribute( attr );
  
  attr = doc->allocate_attribute( "xmlns", "http://physics.nist.gov/N42/2012/N42" );
  RadInstrumentData->append_attribute( attr );
  
  {
    const boost::posix_time::ptime t = boost::posix_time::second_clock::universal_time();
    const string datetime = UtilityFunctions::to_extended_iso_string(t) + "Z";
    val = doc->allocate_string( datetime.c_str(), datetime.size()+1 );
    attr = doc->allocate_attribute( "n42DocDateTime", val );
    RadInstrumentData->append_attribute( attr );
  }
  
  workerpool.post( std::bind( &insert_N42_calibration_nodes,
                               std::cref(measurements_),
                               RadInstrumentData,
                               std::ref(xmldocmutex),
                               std::ref(calToSpecMap) ) );
  
  string original_creator;
  for( size_t i = 0; i < remarks_.size(); ++i )
  {
    if( remarks_[i].empty()
        || UtilityFunctions::starts_with(remarks_[i], "InstrumentVersion:")
        || UtilityFunctions::starts_with(remarks_[i], "Instrument ")
       )
      continue;
    
    if( UtilityFunctions::starts_with(remarks_[i], "N42 file created by: " ) )
    {
      original_creator = remarks_[i].substr(21);
      continue;
    }
    
    std::lock_guard<std::mutex> lock( xmldocmutex );
    val = doc->allocate_string( remarks_[i].c_str(), remarks_[i].size()+1 );
    
    xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
    RadInstrumentData->append_node( remark );
  }//for( size_t i = 0; i < remarks_.size(); ++i )
  
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    if( original_creator.empty() )
    {
      xml_node<char> *RadInstrumentDataCreatorName = doc->allocate_node( node_element, "RadInstrumentDataCreatorName", "InterSpec" );
      RadInstrumentData->append_node( RadInstrumentDataCreatorName );
    }else
    {
      UtilityFunctions::ireplace_all( original_creator, "InterSpec", "" );
      UtilityFunctions::ireplace_all( original_creator, "  ", "" );
      original_creator = "InterSpec. Original file by " + original_creator;
      
      val = doc->allocate_string( original_creator.c_str() );
      xml_node<char> *RadInstrumentDataCreatorName = doc->allocate_node( node_element, "RadInstrumentDataCreatorName", "InterSpec" );
      RadInstrumentData->append_node( RadInstrumentDataCreatorName );
    }//if( original_creator.empty() ) / else
    
//    xml_node<char> *interspec_node = doc->allocate_node( node_element, "InterSpec" );
//    RadInstrumentData->append_node( interspec_node );
  
    //Should add original file name here
//    val = doc->allocate_string( filename_.c_str(), filename_.size()+1 );
//    xml_node<char> *node = doc->allocate_node( node_element, "FileName", val );
//    interspec_node->append_node( node );
  }
  
  
  //We have to convert from 2006 N42 to 2012 instrument typs
  string classcode = convert_n42_instrument_type_from_2006_to_2012( instrument_type_ );
  
  //Class code is required tag.
  if( classcode.empty() )
    classcode = "Other";
  
  string descrip;
  if( lane_number_ >= 0 )
    descrip += "Lane " + std::to_string(lane_number_);
  if( measurement_location_name_.size() )
    descrip += (descrip.size() ? " Location " : "Location ") + measurement_location_name_;
  if( inspection_.size() )
    descrip += (descrip.size() ? " Inspection: " : "Inspection: ") + inspection_;
  
  xml_node<char> *RadInstrumentInformation = 0;
  
  //Child elements of RadInstrumentInformation must be in the order:
  //RadInstrumentManufacturerName
  //RadInstrumentIdentifier
  //RadInstrumentModelName
  //RadInstrumentDescription
  //RadInstrumentClassCode
  //RadInstrumentVersion
  //RadInstrumentQualityControl
  //RadInstrumentCharacteristics
  //RadInstrumentInformationExtension
  
  
  {//Begin lock code block
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    RadInstrumentInformation = doc->allocate_node( node_element, "RadInstrumentInformation" );
    RadInstrumentData->append_node( RadInstrumentInformation );
    
    attr = doc->allocate_attribute( "id", "InstInfo1" );
    RadInstrumentInformation->append_attribute( attr );
  

    if( manufacturer_.size() )
    {
      val = doc->allocate_string( manufacturer_.c_str(), manufacturer_.size()+1 );
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentManufacturerName", val ) );
    }else
    {
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentManufacturerName", "unknown" ) );
    }
    
    if( instrument_id_.size() )
    {
      val = doc->allocate_string( instrument_id_.c_str(), instrument_id_.size()+1 );
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentIdentifier", val ) );
    }
    
    if( instrument_model_.size() )
    {
      val = doc->allocate_string( instrument_model_.c_str(), instrument_model_.size()+1 );
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentModelName", val ) );
    }else
    {
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentModelName", "unknown" ) );
    }
    
    if( descrip.size() )
    {
      val = doc->allocate_string( descrip.c_str(), descrip.size()+1 );
      RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentDescription", val ) );
    }
    
    val = doc->allocate_string( classcode.c_str(), classcode.size()+1 );
    RadInstrumentInformation->append_node( doc->allocate_node( node_element, "RadInstrumentClassCode", val ) );
  }//End lock code block
  
  
  for( size_t i = 0; i < component_versions_.size(); ++i )
  {
    const string &name = component_versions_[i].first;
    const string &version = component_versions_[i].second;
    
    if( UtilityFunctions::icontains( name, "Software") && version == "Unknown" )
      continue;
  
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    xml_node<char> *version_node = version_node = doc->allocate_node( node_element, "RadInstrumentVersion" );
    RadInstrumentInformation->append_node( version_node );
    
    if( UtilityFunctions::iequals( name, "Software" ) )
      val = doc->allocate_string( ("Original " +  name).c_str() );
    else
      val = doc->allocate_string( name.c_str(), name.size()+1 );
    version_node->append_node( doc->allocate_node( node_element, "RadInstrumentComponentName", val ) );
    val = doc->allocate_string( version.c_str(), version.size()+1 );
    version_node->append_node( doc->allocate_node( node_element, "RadInstrumentComponentVersion", val ) );
  }//for( size_t i = 0; i < component_versions_.size(); ++i )
  
  {//Begin lock code block
    /*
     At a minimum, there shall be an instance of this element with the
     component name Software that describes the version of the software and/or
     firmware that produced the current N42 XML document.
     ...
     To determine N42 file compatibility, the user shall always specify the
     “Software” component.
     */
    
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    xml_node<char> *instvrsn = doc->allocate_node( node_element, "RadInstrumentVersion" );
    RadInstrumentInformation->append_node( instvrsn );
    
    xml_node<char> *cmpntname = doc->allocate_node( node_element, "RadInstrumentComponentName", "Software" );
    instvrsn->append_node( cmpntname );
    
    //Could add InterSpec_VERSION
    xml_node<char> *cmpntvsn = doc->allocate_node( node_element, "RadInstrumentComponentVersion", "InterSpec" );
    instvrsn->append_node( cmpntvsn );

    //This next remark is invalid acording to the n42.xsd file, but valid according to the PDF spec
    //xml_node<char> *remark = doc->allocate_node( node_element, "Remark", "N42 file created by InterSpec" );
    //instvrsn->append_node( remark );
    
    
    //Could add InterSpec_VERSION
    
    char buff[8];
    snprintf( buff, sizeof(buff), "%d", MeasurementInfo_2012N42_VERSION );
    
    instvrsn = doc->allocate_node( node_element, "RadInstrumentVersion" );
    RadInstrumentInformation->append_node( instvrsn );
    
    cmpntname = doc->allocate_node( node_element, "RadInstrumentComponentName", "InterSpecN42Serialization" );
    instvrsn->append_node( cmpntname );
  
    val = doc->allocate_string( buff );
    cmpntvsn = doc->allocate_node( node_element, "RadInstrumentComponentVersion", val );
    instvrsn->append_node( cmpntvsn );
  }//End lock code block
  
  //The <RadInstrumentQualityControl> node would be inserted here
  
  xml_node<char> *RadInstrumentCharacteristics = 0;
  xml_node<char> *RadInstrumentInformationExtension = 0;
  
  if(!measurment_operator_.empty())
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    RadInstrumentCharacteristics = doc->allocate_node( node_element, "RadInstrumentCharacteristics" );
    RadInstrumentInformation->append_node( RadInstrumentCharacteristics );
    
    xml_node<char> *CharacteristicGroup = doc->allocate_node( node_element, "CharacteristicGroup" );
    RadInstrumentCharacteristics->append_node( CharacteristicGroup );
    
    xml_node<char> *Characteristic = doc->allocate_node( node_element, "Characteristic" );
    CharacteristicGroup->append_node( Characteristic );
    
    xml_node<char> *CharacteristicName = doc->allocate_node( node_element, "CharacteristicName", "Operator Name" );
    Characteristic->append_node( CharacteristicName );
    
    string operatorname = measurment_operator_;
//    ireplace_all( operatorname, "\t", "&#009;" );//replace tab characters with tab character code
//    ireplace_all( operatorname, "&", "&#38;" );
    
    val = doc->allocate_string( operatorname.c_str(), operatorname.size()+1 );
    xml_node<char> *CharacteristicValue = doc->allocate_node( node_element, "CharacteristicValue", val );
    Characteristic->append_node( CharacteristicValue );
  }//if( measurment_operator_.size() )
  
  
  if( detector_type_ != kUnknownDetector )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    if( !RadInstrumentInformationExtension )
      RadInstrumentInformationExtension = doc->allocate_node( node_element, "RadInstrumentInformationExtension" );
    
    const string &typestr = detectorTypeToString(detector_type_);
    xml_node<char> *type = doc->allocate_node( node_element, "InterSpec:DetectorType", typestr.c_str() );
    RadInstrumentInformationExtension->append_node( type );
  }
  
  if( RadInstrumentInformationExtension )
    RadInstrumentInformation->append_node( RadInstrumentInformationExtension );
  
  
  const size_t ndetectors = detector_names_.size();
  for( size_t i = 0; i < ndetectors; ++i )
  {
    std::lock_guard<std::mutex> lock( xmldocmutex );
    
    xml_node<char> *RadDetectorInformation = doc->allocate_node( node_element, "RadDetectorInformation" );
    RadInstrumentData->append_node( RadDetectorInformation );
    
    if( detector_names_[i].empty() )
      val = s_unnamed_det_placeholder.c_str();
    else
      val = doc->allocate_string( detector_names_[i].c_str(), detector_names_[i].size()+1 );
      
    xml_attribute<> *name = doc->allocate_attribute( "id", val );
    RadDetectorInformation->append_attribute( name );
    
    //<RadDetectorName>: free-form text name for the detector description, ex, "Gamma Front/Right"

    string rad_det_desc;
    for( size_t j = 0; j < measurements_.size(); ++j )
    {
      if( detector_names_[i] == measurements_[j]->detector_name_ )
      {
        rad_det_desc = measurements_[j]->detector_type_;
        break;
      }
    }
    
    if( std::count(neutron_detector_names_.begin(), neutron_detector_names_.end(), detector_names_[i] ) )
    {
      bool hasgamma = false;
      
      for( int sample : sample_numbers_ )
      {
        MeasurementConstShrdPtr m = measurement( sample, detector_numbers_[i] );
        if( m && m->gamma_counts_ && !m->gamma_counts_->empty())
        {
          hasgamma = true;
          break;
        }
      }//for( int sample : sample_numbers_ )
      
      if( hasgamma )
      {
        xml_node<char> *RadDetectorCategoryCode = doc->allocate_node( node_element, "RadDetectorCategoryCode", "Gamma" );
        RadDetectorInformation->append_node( RadDetectorCategoryCode );
        
        if( !rad_det_desc.empty() )
          rad_det_desc += ", ";
        rad_det_desc += "Gamma and Neutron";
      }else
      {
        xml_node<char> *RadDetectorCategoryCode = doc->allocate_node( node_element, "RadDetectorCategoryCode", "Neutron" );
        RadDetectorInformation->append_node( RadDetectorCategoryCode );
      }
    }else
    {
      xml_node<char> *RadDetectorCategoryCode = doc->allocate_node( node_element, "RadDetectorCategoryCode", "Gamma" );
      RadDetectorInformation->append_node( RadDetectorCategoryCode );
    }
    
    if(!rad_det_desc.empty())
    {
      val = doc->allocate_string( rad_det_desc.c_str(), rad_det_desc.size()+1 );
      xml_node<char> *RadDetectorDescription = doc->allocate_node( node_element, "RadDetectorDescription", val );
      RadDetectorInformation->append_node( RadDetectorDescription );
    }
    
    //"RadDetectorDescription"
    
    
    
    const string det_kind = determine_rad_detector_kind_code();
    val = doc->allocate_string( det_kind.c_str() );
    xml_node<char> *RadDetectorKindCode = doc->allocate_node( node_element, "RadDetectorKindCode", val );
    RadDetectorInformation->append_node( RadDetectorKindCode );
    
//    if( det_kind == "Other" )
//    {
//      val = "InterSpec could not determine detector type.";
//      xml_node<char> *RadDetectorDescription = doc->allocate_node( node_element, "RadDetectorDescription", val );
//      RadDetectorInformation->append_node( RadDetectorDescription );
//    }
    
//    xml_node<char> *RadDetectorInformationExtension = doc->allocate_node( node_element, "RadDetectorInformationExtension" );
//    RadInstrumentData->append_node( RadDetectorInformationExtension );
    
//detector_numbers_ are assigned in cleanup_after_load(), so dont need to be saved to the file
//    snprintf( buffer, sizeof(buffer), "%i", detector_numbers_[i] );
//    val = doc->allocate_string( buffer );
//    xml_node<char> *DetectorNumber = doc->allocate_node( node_element, "InterSpec:DetectorNumber", val );
//    RadDetectorInformationExtension->append_node( DetectorNumber );
  }//for( size_t i = 0; i < ndetectors; ++i )
  
  workerpool.join();  //ensures all calibrations have been written to the DOM and can be referenced using calToSpecMap
  
  //For the case of portal data, where the first sample is a long backgorund and
  //  the rest of the file is the occupancy, GADRAS has some "special"
  //  requirements on id attribute values...
  bool first_sample_was_back = false;
  const vector<int> sample_nums_vec( sample_numbers_.begin(), sample_numbers_.end() );
  
  for( set<int>::const_iterator sni = sample_numbers_.begin(); sni != sample_numbers_.end(); ++sni )
  {
    const int sample_num = *sni;
    const vector< std::shared_ptr<const Measurement> > smeas
                                            = sample_measurements( sample_num );
    if( smeas.empty() )  //probably shouldnt happen, but jic
      continue;
    
    vector<size_t> calid;
    for( size_t i = 0; i < smeas.size(); ++i )
    {
      const ShrdConstFVecPtr &binning = smeas[i]->channel_energies();
      calid.push_back( calToSpecMap[binning] );
    }
    
    
    //Check if all smeas have same start and real times, if so, write them under
    //  a single <RadMeasurement> tag, otherwise write under sepearte ones.
    //Note: this test is abit brittle, and maybe not quite fair.  For example on
    //  my portal the neutron
    bool allsame = true;
    boost::posix_time::ptime starttime = smeas[0]->start_time();
    float rtime = smeas[0]->real_time();

    //Bellow allows up to (an arbitrarily chosen) 50 ms difference in start and
    //  real times to allow rounding errors during parsing; There are also some
    //  multi-detector systems (ex Anthony NM HPGe)  where the various
    //  sub-detectors arnet perfectly synced, but they are close enough to
    //  capture the intent of the <RadMeasurement> node.  Note that this does
    //  mess with real time (and dead time) a little for some RPM systems where
    //  the real time for each subdetector isnt always the same for each
    //  timeslice; to make up for this, a <Remark>LiveTime: ...</Remark> live
    //  time is added in these situations...
    for( size_t i = 1; i < smeas.size(); ++i )
    {
      const boost::posix_time::ptime tst = smeas[i]->start_time();
      starttime = ((tst.is_special() || (starttime < tst)) ? starttime : tst);
      rtime = max( rtime, smeas[i]->real_time() );
    }
    
    for( size_t i = 1; allsame && (i < smeas.size()); ++i )
    {
      std::shared_ptr<const Measurement> m = smeas[i];
      
      //Allow the first record in the file, if it is a background, to not need
      //  to have the same start time or real time as the other measurments in
      //  this sample.  This is to accomidate some portals whos RSPs may
      //  accumulate backgrounds a little weird, but none-the-less should all be
      //  considered teh same sample.
      if( m->source_type() != Measurement::SourceType::Background
         || (sample_num != (*sample_numbers_.begin())) )
      {
        if( starttime != m->start_time() )
        {
          const boost::posix_time::time_duration diff = starttime - m->start_time();
          if( diff.is_special() || std::abs( diff.total_microseconds() ) > 50000 )
            allsame = false;
        }
        if( fabs(rtime - m->real_time()) > 0.05f )
          allsame = false;
      }//if( not background or not first sample )
      
      allsame = allsame && (smeas[i]->source_type() == smeas[0]->source_type());
    }//for( check if all measurments can go under a single <RadMeasurement> tag )
    
    
    //allsame = false;
    
    if( allsame )
    {
      char RadMeasurementId[32];
      xml_node<char> *RadMeasurement = 0;
      
      if( passthrough() )
      {
        //Catch the case where its a searchmode file, or a portal file where
        //  the first sample is a long background.  Apparently GADRAS relies
        //  on the "id" attribute of RadMeasurement for this...
        if( (sample_num == (*sample_numbers_.begin()))
            && (smeas[0]->source_type() == Measurement::SourceType::Background)
            && smeas[0]->live_time() > 10.0 )
        {
          first_sample_was_back = true;
          snprintf( RadMeasurementId, sizeof(RadMeasurementId), "Background" );
        }else
        {
          int sn = sample_num;
          if( first_sample_was_back )
          {
            auto pos = lower_bound(begin(sample_nums_vec), end(sample_nums_vec), sample_num );
            sn = (pos - begin(sample_nums_vec));
          }
          snprintf( RadMeasurementId, sizeof(RadMeasurementId), "Survey %i", sn );
        }
      }else
      {
        snprintf( RadMeasurementId, sizeof(RadMeasurementId), "Sample%i", sample_num );
      }

      {
        std::lock_guard<std::mutex> lock( xmldocmutex );
        RadMeasurement = doc->allocate_node( node_element, "RadMeasurement" );
        RadInstrumentData->append_node( RadMeasurement );
      
        char *val = doc->allocate_string( RadMeasurementId );
        xml_attribute<> *attr = doc->allocate_attribute( "id", val );
        RadMeasurement->append_attribute( attr );
      }
      
      workerpool.post( std::bind( &MeasurementInfo::add_spectra_to_measurment_node_in_2012_N42_xml,
                                    RadMeasurement, smeas, calid,
                                   std::ref(xmldocmutex) ) );
    }else //if( allsame )
    {
      for( size_t i = 0; i < smeas.size(); ++i )
      {
        const vector< std::shared_ptr<const Measurement> > thismeas(1,smeas[i]);
        const vector<size_t> thiscalid( 1, calid[i] );
        
        char RadMeasurementId[32];
        xml_node<char> *RadMeasurement = 0;
        snprintf( RadMeasurementId, sizeof(RadMeasurementId), "Sample%iDet%i", sample_num, smeas[i]->detector_number_ );
        
        {
          std::lock_guard<std::mutex> lock( xmldocmutex );
          RadMeasurement = doc->allocate_node( node_element, "RadMeasurement" );
          RadInstrumentData->append_node( RadMeasurement );
          
          char *val = doc->allocate_string( RadMeasurementId );
          xml_attribute<> *attr = doc->allocate_attribute( "id", val );
          RadMeasurement->append_attribute( attr );
        }
        
        workerpool.post( std::bind( &MeasurementInfo::add_spectra_to_measurment_node_in_2012_N42_xml,
                                     RadMeasurement, thismeas, thiscalid,
                                     std::ref(xmldocmutex) ) );
      }//for( loop over measuremtns for this sample number )
    }//if( allsame ) / else
  }//for( )

/*
  for( size_t i = 0; i < measurements_.size(); ++i )
  {
    const ShrdConstFVecPtr &binning = measurements_[i]->channel_energies();

    workerpool.post( std::bind( &Measurement::add_to_2012_N42_xml, measurements_[i],
                                  RadInstrumentData,
                                 calToSpecMap[binning],
                                 std::ref(xmldocmutex) ) );
  }//for( size_t i = 0; i < measurements_.size(); ++i )
*/
  
  if( !!detectors_analysis_ )
    add_analysis_results_to_2012_N42( *detectors_analysis_,
                                 RadInstrumentData, xmldocmutex );
  
  workerpool.join();
  
  
  return doc;
}//rapidxml::xml_node<char> *create_2012_N42_xml() const


void MeasurementInfo::add_spectra_to_measurment_node_in_2012_N42_xml(
                                   ::rapidxml::xml_node<char> *RadMeasurement,
                                   const std::vector< std::shared_ptr<const Measurement> > measurments,
                                   const std::vector<size_t> calibids,
                                   std::mutex &xmldocmutex )
{
  using namespace ::rapidxml;
  
  try
  {
    //Some checks that should never actualy trigger
    if( !RadMeasurement )
      throw runtime_error( "null RadMeasurement" );
    if( measurments.empty() )
      throw runtime_error( "with empty input" );
    if( measurments.size() != calibids.size() )
      throw runtime_error( "measurments.size != calibids.size" );
 
    string radMeasID;
    xml_document<char> *doc = 0;
    
    {
      std::lock_guard<std::mutex> lock( xmldocmutex );
      doc = RadMeasurement->document();
      radMeasID = xml_value_str( XML_FIRST_ATTRIB(RadMeasurement, "id") );
    }
    
    const char *val = 0;
    char buffer[256];
    
    //not dealing with radItemInformationReferences and radMeasurementGroupReferences attributes
    
    //Need child tags of <RadMeasurement> in following order
    //MeasurementClassCode, exactly once
    //StartDateTime, exactly once
    //RealTimeDuration, exactly once
    
    //Spectrum, 0 or more
    //GrossCounts, 0 or more
    //DoseRate, 0 or more
    //TotalDose, 0 or more
    //ExposureRate, 0 or more
    //TotalExposure, 0 or more
    
    //RadInstrumentState, 0 or more
    //RadDetectorState, 0 or more
    //RadItemState, 0 or more
    //OccupancyIndicator, 0 or more
    //RadMeasurementExtension, 0 or more
    
    //Since all samples might not have occupancy/speed/gps info, lets loop
    // through and grab it.  Not this is an artifact of this code not
    // originally being modeled after N42 2012.  In principle, this loop
    // shouldnt have an effect in the vast majority (maybe all I know of) of
    // situations, but jic
    float speed = measurments[0]->speed_;
    boost::posix_time::ptime starttime = measurments[0]->start_time();
    
    Measurement::OccupancyStatus occupancy = measurments[0]->occupied_;
    Measurement::SourceType source_type = measurments[0]->source_type();
    
    bool has_gps = false;
    string positiontime;
    char latitude[16], longitude[16];
    float realtime_used = measurments[0]->real_time_;
    
    for( size_t i = 0; i < measurments.size(); ++i )
    {
      realtime_used = max( measurments[i]->real_time_, realtime_used );
      const boost::posix_time::ptime tst = measurments[i]->start_time();
      starttime = ((tst.is_special() || (starttime < tst)) ? starttime : tst);
      
      speed = max( measurments[i]->speed_, speed );
        
      if( measurments[i]->occupied_ == Measurement::Occupied )
        occupancy = measurments[i]->occupied_;
      else if( occupancy == Measurement::UnknownOccupancyStatus )
        occupancy = measurments[i]->occupied_;
      else if( measurments[i]->occupied_ == Measurement::NotOccupied && occupancy == Measurement::UnknownOccupancyStatus )
        occupancy = measurments[i]->occupied_;
        
      if( !has_gps && measurments[i]->has_gps_info() )
      {
        has_gps = true;
        snprintf( latitude, sizeof(latitude), "%.12f", measurments[i]->latitude_ );
        snprintf( longitude, sizeof(longitude), "%.12f", measurments[i]->longitude_ );
        if( !measurments[i]->position_time_.is_special() )
          positiontime = UtilityFunctions::to_extended_iso_string(measurments[i]->position_time_) + "Z";
      }//if( !has_gps )
        
      if( measurments[i]->source_type_ != Measurement::UnknownSourceType )
        source_type = std::max( measurments[i]->source_type_, source_type );
    }//for( size_t i = 1; i < measurments.size(); ++i )
  
    
    char realtime[32], speedstr[32];

    snprintf( realtime, sizeof(realtime), "PT%fS", realtime_used );
    snprintf( speedstr, sizeof(speedstr), "%.8f", speed );
  
    const string startstr = UtilityFunctions::to_extended_iso_string(starttime) + "Z";
    
    const char *classcode = (const char *)0;
    const char *occupied = (const char *)0;
    switch( source_type )
    {
      case Measurement::Background:         classcode = "Background";        break;
      case Measurement::Calibration:        classcode = "Calibration";       break;
      case Measurement::Foreground:         classcode = "Foreground";        break;
      case Measurement::IntrinsicActivity:  classcode = "IntrinsicActivity"; break;
      case Measurement::UnknownSourceType:  classcode = "NotSpecified";      break;
    }//switch( source_type_ )
    
    switch( occupancy )
    {
      case Measurement::NotOccupied: occupied = "false"; break;
      case Measurement::Occupied:    occupied = "true";  break;
      case Measurement::UnknownOccupancyStatus:          break;
    }//switch( occupied_ )
    
    {
      std::lock_guard<std::mutex> lock( xmldocmutex );
      val = doc->allocate_string( classcode );
      RadMeasurement->append_node( doc->allocate_node( node_element, "MeasurementClassCode", val ) );
      
      if( !measurments[0]->start_time_.is_special() )
      {
        val = doc->allocate_string( startstr.c_str(), startstr.size()+1 );
        RadMeasurement->append_node( doc->allocate_node( node_element, "StartDateTime", val, 13, startstr.size() ) );
      }
      
      if( measurments[0]->real_time_ > 0.0f )
      {
        val = doc->allocate_string( realtime );
        RadMeasurement->append_node( doc->allocate_node( node_element, "RealTimeDuration", val ) );
      }
    }
    
    //Since gross count nodes have to come after
    vector< xml_node<char> *> spectrum_nodes, gross_count_nodes, det_states;
    
    for( size_t i = 0; i < measurments.size(); ++i )
    {
      const size_t calibid = calibids[i];
      const std::shared_ptr<const Measurement> m = measurments[i];
      
      char livetime[32], calibstr[16], spec_idstr[48];
    
      string neutcounts;
      snprintf( livetime, sizeof(livetime), "PT%fS", m->live_time_ );
      snprintf( calibstr, sizeof(calibstr), "EnergyCal%i", static_cast<int>(calibid) );

      if( UtilityFunctions::icontains(radMeasID, "Det") )
      {
        //This is case where all measurements of a sample number did not have a similar
        //  start time or background/foreground status so each sample/detector
        //  gets its own <RadMeasurement> element, with an id like "Sample3Det1"
        snprintf( spec_idstr, sizeof(spec_idstr), "%sSpectrum", radMeasID.c_str() );
      }else if( !radMeasID.empty() )
      {
        //radMeasID will be "Background", "Survey XXX" if passthrough() that
        //  starts with a long background, and "SampleXXX" otherwise.
        snprintf( spec_idstr, sizeof(spec_idstr), "%sDet%iSpectrum", radMeasID.c_str(), m->detector_number_ );
      }else
      {
        //Probably shouldnt ever make it here.
        snprintf( spec_idstr, sizeof(spec_idstr), "Sample%iDet%iSpectrum", m->sample_number_, m->detector_number_ );
      }
      

      const string detnam = !m->detector_name_.empty() ? m->detector_name_ : s_unnamed_det_placeholder;
      
      //Bellow choice of zero compressing if the gamma sum is less than 15 times the
      //  number of gamma channels is arbitrarily chosen, and has not been
      //  benchmarked or checked it is a reasonable value
      const bool zerocompressed = (!!m->gamma_counts_ && (m->gamma_count_sum_<15.0*m->gamma_counts_->size()));
      vector<float> compressedchannels;
  
      if( zerocompressed )
        compress_to_counted_zeros( *m->gamma_counts_, compressedchannels );
  
      const vector<float> &data = (zerocompressed || !m->gamma_counts_) ? compressedchannels : (*m->gamma_counts_);
  
      string channeldata;
      if( !zerocompressed )
        channeldata.reserve( 3*m->gamma_counts_->size() ); //3 has not been verified to be reasonalbe

      const size_t nchannel = data.size();
  
      //The hope is that writing 8 channels data at a time will be faster than one
      //  at a time - however I didnt check that it is, or check that doing somrthign
      //  like 16 or 32 would be faster.
      //"%.9G" specifies use exponential form (i.e. "1.234E5") if shorter than
      //  decimal (i.e 123450), printing up to 9 significant digits.  Also, it looks
      //  like the shortest expressible form of integers are used (e.g. 0.0f prints
      //  as "0", and 101.0f prints as "101").  The maximum sig figs is specified
      //  since floats get converted to doubles when given as arguments of printf.
      //  Also, 8 was chosen since we have integer acuracy of floats up to
      //  16,777,216 (after this floats have less precision than int).
      //  Also note that if we wanted to garuntee a round-trip of float-text-float
      //  we could use "%1.8e" or "%.9g".
      //For a lot of great float information, see:
      //  https://randomascii.wordpress.com/2013/02/07/float-precision-revisited-nine-digit-float-portability/
      if( (nchannel % 8) == 0 )
      {
        for( size_t i = 0; i < nchannel; i += 8 )
        {
          snprintf( buffer, sizeof(buffer),
                   (i?" %.8G %.8G %.8G %.8G %.8G %.8G %.8G %.8G"
                     :"%.8G %.8G %.8G %.8G %.8G %.8G %.8G %.8G"),
                   data[i], data[i+1], data[i+2], data[i+3],
                  data[i+4], data[i+5], data[i+6], data[i+7] );
          channeldata += buffer;
        }//for( size_t i = 0; i < nchannel; i += 8 )
      }else
      {
        for( size_t i = 0; i < nchannel; ++i )
        {
          snprintf( buffer, sizeof(buffer), (i?" %.8G":"%.8G"), data[i] );
          channeldata += buffer;
        }//for( size_t i = 0; i < nchannel; i += 8 )
      }//if( (nchannel % 8) == 0 )
  
      if( m->neutron_counts_.size() > 1 )
      {
        for( size_t i = 0; i < m->neutron_counts_.size(); ++i )
        {
          snprintf( buffer, sizeof(buffer), (i?" %.8G":"%.8G"), m->neutron_counts_[i] );
          neutcounts += buffer;
        }//for( size_t i = 0; i < nchannel; i += 8 )
      }else
      {
        snprintf( buffer, sizeof(buffer), "%.8G", m->neutron_counts_sum_ );
        neutcounts += buffer;
      }
  
    
      std::lock_guard<std::mutex> lock( xmldocmutex );
  
      if( !!m->gamma_counts_ && !m->gamma_counts_->empty())
      {
        xml_node<char> *Spectrum = doc->allocate_node( node_element, "Spectrum" );
        spectrum_nodes.push_back( Spectrum );
  
        //If there is a slight mismatch between the live times of this sample
        //  (~50 ms), we will still include all detectors in the same sample,
        //  but put in a remark notting a difference.  This is absolutely a
        //  hack, but some sort of comprimise is needed to cram stuff into N42
        //  2012 files from arbitrary sources.
        if( fabs(m->real_time_ - realtime_used) > 0.00001 )
        {
          char thisrealtime[64];
          snprintf( thisrealtime, sizeof(thisrealtime), "RealTime: PT%fS", m->real_time_ );
          val = doc->allocate_string( thisrealtime );
          xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
          Spectrum->append_node( remark );
        }
        
        if(!m->title_.empty())
        {
          const string title = "Title: " + m->title_;
          val = doc->allocate_string( title.c_str(), title.size()+1 );
          xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
          Spectrum->append_node( remark );
        }
      
        for( size_t i = 0; i < m->remarks_.size(); ++i )
        {
          if( m->remarks_[i].empty() )
            continue;
          const char *val = doc->allocate_string( m->remarks_[i].c_str(), m->remarks_[i].size()+1 );
          xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
          Spectrum->append_node( remark );
        }//for( size_t i = 0; i < remarks_.size(); ++i )
      
        val = doc->allocate_string( calibstr );
        Spectrum->append_attribute( doc->allocate_attribute( "energyCalibrationReference", val ) );

        val = doc->allocate_string( detnam.c_str(), detnam.size()+1 );
        Spectrum->append_attribute( doc->allocate_attribute( "radDetectorInformationReference", val, 31, detnam.size() ) );
  
        //Add required ID attribute
        val = doc->allocate_string( spec_idstr );
        Spectrum->append_attribute( doc->allocate_attribute( "id", val ) );
    
        if( m->live_time_ > 0.0f )
        {
          val = doc->allocate_string( livetime );
          xml_node<char> *LiveTimeDuration = doc->allocate_node( node_element, "LiveTimeDuration", val );
          Spectrum->append_node( LiveTimeDuration );
        }//if( live_time_ > 0.0f )
  
        if(!channeldata.empty())
        {
          val = doc->allocate_string( channeldata.c_str(), channeldata.size()+1 );
          xml_node<char> *ChannelData = doc->allocate_node( node_element, "ChannelData", val, 11, channeldata.size() );
          Spectrum->append_node( ChannelData );
   
          if( zerocompressed )
            ChannelData->append_attribute( doc->allocate_attribute( "compressionCode", "CountedZeroes" ) );
        }//if( channeldata.size() )
      }//if( !!gamma_counts_ && gamma_counts_->size() )
    
      if( m->contained_neutron_ )
      {
        xml_node<char> *GrossCounts = doc->allocate_node( node_element, "GrossCounts" );
        gross_count_nodes.push_back( GrossCounts );
    
        if( (!m->gamma_counts_ || m->gamma_counts_->empty())  )
        {
          if( fabs(m->real_time_ - realtime_used) > 0.00001 )
          {
            char thisrealtime[64];
            snprintf( thisrealtime, sizeof(thisrealtime), "RealTime: PT%fS", m->real_time_ );
            val = doc->allocate_string( thisrealtime );
            xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
            GrossCounts->append_node( remark );
          }
          
          if(!m->title_.empty())
          {
            const string title = "Title: " + m->title_;
            val = doc->allocate_string( title.c_str(), title.size()+1 );
            xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
            GrossCounts->append_node( remark );
          }//if( m->title_.size() )
          
          for( size_t i = 0; i < m->remarks_.size(); ++i )
          {
            if( m->remarks_[i].empty() )
              continue;
            const char *val = doc->allocate_string( m->remarks_[i].c_str(), m->remarks_[i].size()+1 );
            xml_node<char> *remark = doc->allocate_node( node_element, "Remark", val );
            GrossCounts->append_node( remark );
          }//for( size_t i = 0; i < remarks_.size(); ++i )
        }//if( (!m->gamma_counts_ || m->gamma_counts_->empty())  )
    
        char neutId[32];
        if( radMeasID.empty() )
          snprintf( neutId, sizeof(neutId), "Sample%iNeutron%i", m->sample_number_, m->detector_number_ );
        else
          snprintf( neutId, sizeof(neutId), "%sNeutron%i", radMeasID.c_str(), m->detector_number_ );
        
        val = doc->allocate_string( neutId );
        GrossCounts->append_attribute( doc->allocate_attribute( "id", val ) );
        val = doc->allocate_string( detnam.c_str(), detnam.size()+1 );
        GrossCounts->append_attribute( doc->allocate_attribute( "radDetectorInformationReference", val, 31, detnam.size() ) );
    
        val = doc->allocate_string( livetime );
        xml_node<char> *LiveTimeDuration = doc->allocate_node( node_element, "LiveTimeDuration", val );
        GrossCounts->append_node( LiveTimeDuration );
    
        val = doc->allocate_string( neutcounts.c_str(), neutcounts.size()+1 );
        xml_node<char> *CountData = doc->allocate_node( node_element, "CountData", val, 9, neutcounts.size() );
        GrossCounts->append_node( CountData );
      }//if( contained_neutron_ )
      
      
      switch( measurments[i]->quality_status_ )
      {
        case Measurement::Good:
          //When reading in the 2012 N42, we will assume good unless indicated otherwise
        break;
            
        case Measurement::Suspect: case Measurement::Bad:
        {
          xml_node<char> *RadDetectorState = doc->allocate_node( node_element, "RadDetectorState" );
          det_states.push_back( RadDetectorState );
          
          const char *val = doc->allocate_string( detnam.c_str() );
          xml_attribute<char> *att = doc->allocate_attribute( "radDetectorInformationReference", val );
          RadDetectorState->append_attribute( att );
          
          val = ((measurments[i]->quality_status_==Measurement::Suspect) ? "Warning" : "Fatal" ); //"Error" is also an option
          RadDetectorState->append_node( doc->allocate_node( node_element, "Fault", val ) );
          break;
        }//case Suspect: case Bad:
            
        case Measurement::Missing:
        {
          //This next line is InterSpec specific for round-tripping files
          xml_node<char> *RadDetectorState = doc->allocate_node( node_element, "RadDetectorState" );
          det_states.push_back( RadDetectorState );
          
          const char *val = doc->allocate_string( detnam.c_str() );
          xml_attribute<char> *att = doc->allocate_attribute( "radDetectorInformationReference", val );
          RadDetectorState->append_attribute( att );
            
          xml_node<char> *remark = doc->allocate_node( node_element, "Remark", "InterSpec could not determine detector state." );
          RadDetectorState->append_node( remark );
          break;
        }
      }//switch( quality_status_ )
    }//for( loop over input measurments )
  
    
    {//start put <Spectrum> and <GrossCount> nodes into tree
      std::lock_guard<std::mutex> lock( xmldocmutex );
      for( size_t i = 0; i < spectrum_nodes.size(); ++i )
        RadMeasurement->append_node( spectrum_nodes[i] );
      for( size_t i = 0; i < gross_count_nodes.size(); ++i )
        RadMeasurement->append_node( gross_count_nodes[i] );
    }//end put <Spectrum> and <GrossCount> nodes into tree
    
    
    {//begin add other information
      std::lock_guard<std::mutex> lock( xmldocmutex );
      
      if( has_gps )
      {
        xml_node<char> *RadInstrumentState = doc->allocate_node( node_element, "RadInstrumentState" );
        RadMeasurement->append_node( RadInstrumentState );
        
        xml_node<char> *StateVector = doc->allocate_node( node_element, "StateVector" );
        RadInstrumentState->append_node( StateVector );
        
        xml_node<char> *GeographicPoint = doc->allocate_node( node_element, "GeographicPoint" );
        StateVector->append_node( GeographicPoint );
        
        
        val = doc->allocate_string( latitude );
        xml_node<char> *LatitudeValue = doc->allocate_node( node_element, "LatitudeValue", val );
        GeographicPoint->append_node( LatitudeValue );
        
        val = doc->allocate_string( longitude );
        xml_node<char> *LongitudeValue = doc->allocate_node( node_element, "LongitudeValue", val );
        GeographicPoint->append_node( LongitudeValue );
        
        //<PositionTime> is an InterSpec addition since it didnt look like there wa a place for it in the spec
        if(!positiontime.empty())
        {
          val = doc->allocate_string( positiontime.c_str(), positiontime.size()+1 );
          xml_node<char> *PositionTime = doc->allocate_node( node_element, "PositionTime", val, 12, positiontime.size() );
          GeographicPoint->append_node( PositionTime );
        }
      }//if( has_gps_info() )
      
      for( size_t i = 0; i < det_states.size(); ++i )
        RadMeasurement->append_node( det_states[i] );
      
      if( speed > 0.0f )
      {
        xml_node<char> *RadItemState = doc->allocate_node( node_element, "RadItemState" );
        RadMeasurement->append_node( RadItemState );
        
        xml_node<char> *StateVector = doc->allocate_node( node_element, "StateVector" );
        RadItemState->append_node( StateVector );
        
        val = doc->allocate_string( speedstr );
        xml_node<char> *SpeedValue = doc->allocate_node( node_element, "SpeedValue", val );
        StateVector->append_node( SpeedValue );
      }//if( speed_ > 0 )
    
      if( occupied )
      {
        val = doc->allocate_string( occupied );
        RadMeasurement->append_node( doc->allocate_node( node_element, "OccupancyIndicator", val ) );
      }
    }//end add other information
    

//Potential child nodes of <RadMeasurement> we could
//<GrossCounts>, <DoseRate>, <TotalDose>, <ExposureRate>, <TotalExposure>,
//  <RadInstrumentState>, <RadDetectorState>, <RadItemState>, <RadMeasurementExtension>
    
  }catch( std::exception &e )
  {
    cerr << "Measurement::add_spectra_to_measurment_node_in_2012_N42_xml(...): something horrible happened!: " << e.what() << endl;
  }//try catch
}//void Measurement::add_to_2012_N42_xml(...)



bool MeasurementInfo::write_2012_N42( std::ostream& ostr ) const
{
  std::shared_ptr< rapidxml::xml_document<char> > doc = create_2012_N42_xml();
  
//  if( !!doc )
//    rapidxml::print( ostr, *doc, rapidxml::print_no_indenting );
  if( !!doc )
    rapidxml::print( ostr, *doc );

  
  return !!doc;
}//bool write_2012_N42( std::ostream& ostr ) const


std::string MeasurementInfo::concat_2012_N42_characteristic_node( const rapidxml::xml_node<char> *char_node )
{
  //      const rapidxml::xml_attribute<char> *char_id = char_node->first_attribute( "id", 2 );
  const rapidxml::xml_attribute<char> *date = char_node->first_attribute( "valueDateTime", 13 );
  const rapidxml::xml_attribute<char> *limits = char_node->first_attribute( "valueOutOfLimits", 16 );
  const rapidxml::xml_node<char> *remark_node = char_node->first_node( "Remark", 6 );
  const rapidxml::xml_node<char> *name_node = char_node->first_node( "CharacteristicName", 18 );
  const rapidxml::xml_node<char> *value_node = char_node->first_node( "CharacteristicValue", 19 );
  const rapidxml::xml_node<char> *unit_node = char_node->first_node( "CharacteristicValueUnits", 24 );
//  const rapidxml::xml_node<char> *class_node = char_node->first_node( "CharacteristicValueDataClassCode", 32 );
  
  string comment;
  if( name_node && name_node->value_size() )
    comment = xml_value_str( name_node );
  
  if( (date && date->value_size()) || (limits && limits->value_size())
     || (remark_node && remark_node->value_size()) )
  {
    comment += "(";
    comment += xml_value_str( date );
    if( limits && limits->value_size() )
    {
      if( comment[comment.size()-1] != '(' )
        comment += ", ";
      comment += "value out of limits: ";
      comment += xml_value_str(limits);
    }//if( limits && limits->value_size() )
    
    if( remark_node && remark_node->value_size() )
    {
      if( comment[comment.size()-1] != '(' )
        comment += ", ";
      comment += "remark: ";
      comment += xml_value_str(remark_node);
    }//if( limits && limits->value_size() )
    comment += ")";
  }//if( an attribute has info )

  if( value_node )
    comment += string(":") + xml_value_str(value_node);
  
  if( unit_node && !XML_VALUE_ICOMPARE(unit_node, "unit-less") )
    comment += " " + xml_value_str(unit_node);
  
  return comment;
}//std::string concat_2012_N42_characteristic_node( const rapidxml::xml_node<char> *node )


void MeasurementInfo::set_2012_N42_instrument_info( const rapidxml::xml_node<char> *info_node )
{
  if( !info_node )
    return;
  
//  const rapidxml::xml_attribute<char> *id = info_node->first_attribute( "id", 2 );
  
  const rapidxml::xml_node<char> *remark_node = info_node->first_node( "Remark", 6 );
  if( remark_node )
    remarks_.push_back( xml_value_str(remark_node) );
  
  
  const rapidxml::xml_node<char> *manufacturer_name_node = info_node->first_node( "RadInstrumentManufacturerName", 29 );
  if( manufacturer_name_node && !XML_VALUE_COMPARE(manufacturer_name_node,"unknown") )
    manufacturer_ = xml_value_str(manufacturer_name_node);  //ex. "FLIR Systems"
  
  const rapidxml::xml_node<char> *instr_id_node = info_node->first_node( "RadInstrumentIdentifier", 23, false );
  if( instr_id_node && instr_id_node->value_size() )
    instrument_id_ = xml_value_str(instr_id_node);          //ex. "932222-76"
  
  const rapidxml::xml_node<char> *model_node = info_node->first_node( "RadInstrumentModelName", 22 );
  
  if( !model_node )  //file_format_test_spectra/n42_2006/identiFINDER/20130228_184247Preliminary2010.n42
    model_node = info_node->first_node( "RadInstrumentModel", 18 );
  if( model_node && !XML_VALUE_COMPARE(model_node,"unknown") )
    instrument_model_ = xml_value_str( model_node );          //ex. "identiFINDER 2 ULCS-TNG"
  
  const rapidxml::xml_node<char> *desc_node = info_node->first_node( "RadInstrumentDescription", 24 );
  //Description: Free-form text describing the radiation measurement instrument.
  //Usage: This information can be a general description of the radiation measurement instrument or its use.
  if( desc_node && desc_node->value_size() )
  {
    string val = xml_value_str( desc_node );
    
    string::size_type lanepos = val.find( "Lane " );
    if( lanepos != string::npos )
    {
      if( !(stringstream(val.substr(lanepos+5) ) >> lane_number_) )
      {
        lanepos = string::npos;
        cerr << "Failed to read lane number from '" << val << "'" << endl;
      }
    }//if( lanepos != string::npos )
    
    
    string::size_type locationpos = val.find( "Location " );
    if( locationpos != string::npos )
      measurement_location_name_ = val.substr(locationpos+9);
    else if( (locationpos = val.find( " at " ))==string::npos )
      measurement_location_name_ = val.substr(locationpos+4);
    
    if( measurement_location_name_.find("Inspection:") != string::npos )
      measurement_location_name_ = measurement_location_name_.substr(0,measurement_location_name_.find("Inspection:"));
    
    string::size_type inspectionpos = val.find( "Inspection: " );
    if( inspectionpos != string::npos )
      inspection_ = val.substr(inspectionpos+12);
    if( inspection_.find("Location ") != string::npos )
      inspection_ = inspection_.substr(0,inspection_.find("Location "));
    
    UtilityFunctions::trim( inspection_ );
    UtilityFunctions::trim( measurement_location_name_ );
    
    if( lanepos==string::npos && (locationpos==string::npos && val.size()<8) )
      remarks_.push_back( string("Instrument Description: ") + xml_value_str( desc_node ) );
  }
  
  
  const rapidxml::xml_node<char> *infoextension_node = info_node->first_node( "RadInstrumentInformationExtension", 33 );
//  const rapidxml::xml_node<char> *operator_node = infoextension_node ? infoextension_node->first_node("InterSpec:MeasurmentOperator",18) : (const rapidxml::xml_node<char> *)0;
  
  //<InterSpec:Inspection> node is vestigual as of 20160607, and only kept in to read older files (with MeasurementInfo_2012N42_VERSION==1), of which, there are probably not many, if any around
  const rapidxml::xml_node<char> *inspection_node = infoextension_node ? infoextension_node->first_node("InterSpec:Inspection",20) : (const rapidxml::xml_node<char> *)0;
  
  const rapidxml::xml_node<char> *detector_type_node = infoextension_node ? infoextension_node->first_node("InterSpec:DetectorType",22) : (const rapidxml::xml_node<char> *)0;
  
//  
  
//  if( operator_node )
//    measurment_operator_ = xml_value_str( operator_node );
  if( inspection_node )
    inspection_ = xml_value_str( inspection_node );
  
  if( detector_type_node )
  {
    const string type = xml_value_str( detector_type_node );
    for( DetectorType i = DetectorType(0); i <= kSam945; i = DetectorType(i+1) )
    {
      if( type == detectorTypeToString(i) )
      {
        detector_type_ = i;
        break;
      }
    }
  }//if( detector_type_node )

  
  const rapidxml::xml_node<char> *class_code_node = info_node->first_node( "RadInstrumentClassCode", 22 );
  instrument_type_ = xml_value_str( class_code_node );
  
  if( UtilityFunctions::iequals( instrument_type_, "Other" ) )
    instrument_type_ = "";
  
  for( const rapidxml::xml_node<char> *version_node = XML_FIRST_NODE(info_node, "RadInstrumentVersion");
      version_node;
      version_node = XML_NEXT_TWIN(version_node) )
  {
    const rapidxml::xml_node<char> *name = version_node->first_node( "RadInstrumentComponentName", 26 );
    const rapidxml::xml_node<char> *version = version_node->first_node( "RadInstrumentComponentVersion", 29 );
    
    //RadInstrumentVersion is requiren for N42 2012, so should add it as another field
    string comment;
    if( version_node->value_size() && !XML_VALUE_COMPARE( version_node, "unknown" ) )
    {
      //This is actually invalid...
      component_versions_.push_back( make_pair("unknown",xml_value_str(version_node) ) );
    }else if( XML_VALUE_COMPARE( version_node, "Software" ) && XML_VALUE_COMPARE( version_node, "InterSpec" ) )
    {
      //Skip this, since it is written by InterSpec.
    }else if( name && version )
    {
      string namestr = xml_value_str(name);
      if( istarts_with(namestr, "Original Software") )
        namestr = namestr.substr( 9 );
      
      component_versions_.push_back( make_pair( namestr, xml_value_str(version) ) );
      //const rapidxml::xml_node<char> *remark = version_node->first_node( "Remark", 6 );
      //if( remark && remark->value_size() )
        //comment += ((name->value_size() || version->value_size()) ? ", " : " ")
                   //+ string("Remark: ") + xml_value_str(remark_node);
    }//if( unknown ) / else
  }//for( loop over "RadInstrumentVersion" nodes )
  
  
  const rapidxml::xml_node<char> *qc_node = info_node->first_node( "RadInstrumentQualityControl", 27 );
  if( qc_node )
  {
    const rapidxml::xml_attribute<char> *id = qc_node->first_attribute( "id", 2 );
    const rapidxml::xml_node<char> *remark_node = qc_node->first_node( "Remark", 6 );
    const rapidxml::xml_node<char> *date_node = qc_node->first_node( "InspectionDateTime", 18 );
    const rapidxml::xml_node<char> *indicator_node = qc_node->first_node( "InCalibrationIndicator", 22 );
    string comment = "Calibration Check";
    if( id && id->value_size() )
      comment = xml_value_str(id) + string(" ") + comment;
    if( date_node && date_node->value_size() )
      comment += string(" ") + xml_value_str(date_node);
    if( indicator_node && indicator_node->value_size() )
      comment += string(" pass=") + xml_value_str(indicator_node);
    if( remark_node && remark_node->value_size() )
      comment += string(", remark: ") + xml_value_str(remark_node);
  }//if( qc_node )
  
  for( const rapidxml::xml_node<char> *charac_node = info_node->first_node( "RadInstrumentCharacteristics", 28 );
      charac_node;
      charac_node = XML_NEXT_TWIN(charac_node) )
  {
//    const rapidxml::xml_attribute<char> *id = charac_node->first_attribute( "id", 2 );
    const rapidxml::xml_node<char> *remark_node = charac_node->first_node( "Remark", 6 );
    
    if( remark_node )
    {
      string remark = xml_value_str( remark_node );
      trim( remark );
      if(!remark.empty())
        remarks_.push_back( remark );
    }
    
    for( const rapidxml::xml_node<char> *char_node = charac_node->first_node( "Characteristic", 14 );
        char_node;
        char_node = XML_NEXT_TWIN(char_node) )
    {
      const string comment = concat_2012_N42_characteristic_node( char_node );
      if(!comment.empty())
        remarks_.push_back( comment );
    }//for( loop over "Characteristic" nodes )
    
    for( const rapidxml::xml_node<char> *group_node = charac_node->first_node( "CharacteristicGroup", 19 );
        group_node;
        group_node = XML_NEXT_TWIN(group_node) )
    {
//      const rapidxml::xml_attribute<char> *char_id = group_node->first_attribute( "id", 2 );
      const rapidxml::xml_node<char> *remark_node = group_node->first_node( "Remark", 6 );
      const rapidxml::xml_node<char> *name_node = group_node->first_node( "CharacteristicGroupName", 23 );
      
      string precursor;
      if( (name_node && name_node->value_size()) || (remark_node && remark_node->value_size()) )
      {
        precursor += "[";
        if( name_node && name_node->value_size() )
          precursor += xml_value_str(name_node);
        if( remark_node )
        {
          if( precursor.size() > 1 )
            precursor += " ";
          precursor += string("(remark: ") + xml_value_str(remark_node) + string(")");
        }
        precursor += "] ";
      }//if( name_node or remark_node )
      
      for( const rapidxml::xml_node<char> *char_node = group_node->first_node( "Characteristic", 14 );
          char_node;
          char_node = XML_NEXT_TWIN(char_node) )
      {
        
        using rapidxml::internal::compare;
        const rapidxml::xml_node<char> *name_node = char_node->first_node( "CharacteristicName", 18 );
        
        
        if( name_node && XML_VALUE_ICOMPARE(name_node, "Operator Name") )
        {
          const rapidxml::xml_node<char> *value_node = char_node->first_node( "CharacteristicValue", 19 );
          measurment_operator_ = xml_value_str(value_node);
        }else
        {
          const string comment = concat_2012_N42_characteristic_node( char_node );
          if(!comment.empty())
            remarks_.push_back( precursor + comment );
        }//
      }//for( loop over "Characteristic" nodes )
    }//for( loop over "CharacteristicGroup" nodes )
  }//for( loop over "RadInstrumentCharacteristics" nodes )
  

//  rapidxml::xml_node<char> *info_extension_node = info_node->first_node( "RadInstrumentInformationExtension", 33 );
//  if( info_extension_node )
//    cerr << "Warning, the RadInstrumentInformationExtension tag isnt "
//         << "implemented for N42 2012 files yet." << endl;
}//void set_2012_N42_instrument_info( rapidxml::xml_node<char> *inst_info_node )


void get_2012_N42_energy_calibrations( map<string,MeasurementCalibInfo> &calibrations,
                            const rapidxml::xml_node<char> *data_node,
                            vector<string> &remarks )
{
  using UtilityFunctions::split_to_floats;
  
  
  for( const rapidxml::xml_node<char> *cal_node = XML_FIRST_NODE(data_node, "EnergyCalibration");
      cal_node;
      cal_node = XML_NEXT_TWIN(cal_node) )
  {
    string id;
    rapidxml::xml_attribute<char> *id_att = cal_node->first_attribute( "id", 2, false );
    
    
    if( !id_att || !id_att->value_size() )
      id_att = cal_node->first_attribute( "Reference", 9, false );
    
    if( id_att && id_att->value_size() )
      id = xml_value_str(id_att);  //if attribute is _required_, so well rely on this
    
    rapidxml::xml_node<char> *remark_node = cal_node->first_node( "Remark", 6 );
    rapidxml::xml_node<char> *coef_val_node = cal_node->first_node( "CoefficientValues", 17 );
    rapidxml::xml_node<char> *energy_boundry_node = cal_node->first_node( "EnergyBoundaryValues", 20 );
    rapidxml::xml_node<char> *date_node = cal_node->first_node( "CalibrationDateTime", 19 );
    
    if( !coef_val_node )
      coef_val_node = cal_node->first_node( "Coefficients", 12 );
    
    if( remark_node && remark_node->value_size() )
      remarks.push_back( "Calibration for " + id + " remark: " + xml_value_str(remark_node) );
    
    if( date_node && date_node->value_size() )
      remarks.push_back( id + " calibrated " + xml_value_str(date_node) );

    MeasurementCalibInfo info;
    if( coef_val_node && coef_val_node->value_size() )
    {
      info.equation_type = Measurement::Polynomial;
      vector<string> fields;
      const char *data = coef_val_node->value();
      const size_t len = coef_val_node->value_size();
      if( !split_to_floats( data, len, info.coefficients ) )
        throw runtime_error( "Invalid calibration value: " + xml_value_str(coef_val_node) );
      
      while(!info.coefficients.empty() && info.coefficients.back()==0.0f )
        info.coefficients.erase( info.coefficients.end()-1 );
      
      if( info.coefficients.size() < 2 )
      {
        cerr << "Warning: found a EnergyCalibration CoefficientValues with "
             << info.coefficients.size() << " coefficients, which isnt enough, "
             << "skipping calibration" << endl;
        continue;
      }//if( info.coefficients.size() < 2 )
      
      info.calib_id = xml_value_str( coef_val_node->first_attribute( "id", 2 ) );
      
      //  XXX - deviation pairs not yet not, because I have no solid examples files.
      //  The EnergyDeviationValues and EnergyValues child elements provide a means
      //  to account for the difference in the energy predicted by the second-order
      //  polynomial equation and the true energy.
      rapidxml::xml_node<char> *energy_values_node = cal_node->first_node( "EnergyValues", 12 );
      rapidxml::xml_node<char> *energy_deviation_node = cal_node->first_node( "EnergyDeviationValues", 21 );
      if( energy_deviation_node && energy_values_node
         && energy_deviation_node->value_size() && energy_values_node->value_size() )
      {
        try
        {
          vector<float> energies, deviations;
          
          const char *energiesstr = energy_values_node->value();
          const size_t energystrsize = energy_values_node->value_size();
          const char *devstrs = energy_deviation_node->value();
          const size_t devstrsize = energy_deviation_node->value_size();
          
          if( !split_to_floats( energiesstr, energystrsize, energies ) )
            throw runtime_error( "" );
        
          if( !split_to_floats( devstrs, devstrsize, deviations ) )
            throw runtime_error( "" );
          
          if( energies.size() != deviations.size() )
            throw runtime_error( "" );
          
          vector< pair<float,float> > &devpairs = info.deviation_pairs_;
          
          for( size_t i = 0; i < energies.size(); ++i )
            devpairs.push_back( make_pair(energies[i],deviations[i]) );
          
          stable_sort( devpairs.begin(), devpairs.end(), &dev_pair_less_than );
        }catch( std::exception & )
        {
          passMessage( "Deviation pairs from file appear to be invalid, not using", "", 3 );
        }
      }//if( energy deviation information is available )
      
    }else if( energy_boundry_node )
    {
      info.equation_type = Measurement::LowerChannelEdge;
      
      const char *data = energy_boundry_node->value();
      const size_t len = energy_boundry_node->value_size();
      
      if( !split_to_floats( data, len, info.coefficients ) )
        throw runtime_error( "Failed to parse lower channel energies" );
    }else
    {
      cerr << "Warning, found an invalid EnergyCalibration node" <<endl;
      continue;
    }
    
    if( calibrations.count(id) != 0 )
      cerr << "Warning, overwriting calibration '" << id << "'" << endl;
    calibrations[id] = info;
  }//for( loop over "Characteristic" nodes )
  
}//get_2012_N42_energy_calibrations(...)




void MeasurementInfo::decode_2012_N42_detector_state_and_quality( MeasurementShrdPtr meas,
                                                           const rapidxml::xml_node<char> *meas_node )
{
  using rapidxml::internal::compare;
  
  if( !meas_node || !meas )
    return;
  
  meas->quality_status_ = Measurement::Good;  //2012 N42 defaults to good
  const rapidxml::xml_node<char> *detector_state_node = meas_node->first_node( "RadDetectorState", 16 );
  
  if( detector_state_node )
  {
    const rapidxml::xml_node<char> *fault = detector_state_node->first_node( "Fault", 5 );
    const rapidxml::xml_node<char> *remark = XML_FIRST_NODE( detector_state_node, "Remark" );
    
    if( fault && fault->value_size() )
    {
      if( XML_VALUE_ICOMPARE( fault, "Fatal" )
          || XML_VALUE_ICOMPARE( fault, "Error" ) )
        meas->quality_status_ = Measurement::Bad;
      else if( XML_VALUE_ICOMPARE( fault, "Warning" ) )
        meas->quality_status_ = Measurement::Suspect;
    }else if( !detector_state_node->first_node() ||
              (remark && UtilityFunctions::starts_with( xml_value_str(remark), "InterSpec could not")) )
    {
      meas->quality_status_ = Measurement::Missing; //InterSpec Specific
    }
  }//if( detector_state_node )


  const rapidxml::xml_node<char> *inst_state_node = meas_node->first_node( "RadInstrumentState", 18 );
  if( !inst_state_node )
    inst_state_node = meas_node->first_node( "RadItemState", 12 );  //bubbletech portal stotes gps info under this tag
  
  if( inst_state_node )
  {
    const rapidxml::xml_node<char> *StateVector = inst_state_node->first_node( "StateVector", 11 );
    const rapidxml::xml_node<char> *GeographicPoint = (StateVector ? StateVector->first_node( "GeographicPoint", 15 ) : (const rapidxml::xml_node<char> *)0);
  
    if( GeographicPoint )
    {
      const rapidxml::xml_node<char> *longitude = GeographicPoint->first_node( "LongitudeValue", 14 );
      const rapidxml::xml_node<char> *latitude = GeographicPoint->first_node( "LatitudeValue", 13 );
//      const rapidxml::xml_node<char> *elevation= GeographicPoint->first_node( "ElevationValue", 14 );
      const rapidxml::xml_node<char> *time = GeographicPoint->first_node( "PositionTime", 12 ); //An InterSpec Addition
    
      if( !longitude )
        longitude = GeographicPoint->first_node( "Longitude", 9 );
      if( !latitude )
        latitude = GeographicPoint->first_node( "Latitude", 8 );
//      if( !elevation )
//        elevation = GeographicPoint->first_node( "Elevation", 9 );
      
      const string longstr = xml_value_str( longitude );
      const string latstr = xml_value_str( latitude );
//      const string elevstr = xml_value_str( elevation );
      const string timestr = xml_value_str( time );
      
      if( longstr.size() )
        meas->longitude_ = atof( longstr.c_str() );
      if( latstr.size() )
        meas->latitude_ = atof( latstr.c_str() );
//      if( elevstr.size() )
//        meas->elevation_ = atof( elevstr.c_str() );
      
      if( timestr.size()
          && Measurement::valid_longitude(meas->longitude_)
          && Measurement::valid_latitude(meas->latitude_) )
        meas->position_time_ =  time_from_string( timestr.c_str() );
      
//      static std::mutex sm;
//      {
//        std::lock_guard<std::mutex> lock(sm);
//        cerr << "time=" << time << ", timestr=" << timestr << ", meas->position_time_=" << UtilityFunctions::to_iso_string(meas->position_time_) << endl;
//      }
    }//if( GeographicPoint )
  }//if( RadInstrumentState )
  
  rapidxml::xml_node<char> *extension_node = meas_node->first_node( "RadMeasurementExtension", 23 );
  
  if( extension_node )
  {
    //This is vestigial for MeasurementInfo_2012N42_VERSION==2
    rapidxml::xml_node<char> *title_node = extension_node->first_node( "InterSpec:Title", 15 );
    meas->title_ = xml_value_str( title_node );
    
    //This is vestigial for MeasurementInfo_2012N42_VERSION==1
    rapidxml::xml_node<char> *type_node = extension_node->first_node( "InterSpec:DetectorType", 22 );
    meas->detector_type_ = xml_value_str( type_node );
  }//if( detector_type_.size() || title_.size() )
}//void decode_2012_N42_detector_state_and_quality(...)


void MeasurementInfo::decode_2012_N42_rad_measurment_node(
                                     vector< MeasurementShrdPtr > &measurments,
                                     const rapidxml::xml_node<char> *meas_node,
                                     const IdToDetectorType *id_to_dettype_ptr,
                                     DetectorToCalibInfo *calibrations_ptr,
                                     std::mutex &meas_mutex,
                                     std::mutex &calib_mutex )
{
//  static std::mutex cerrmutex;
  
  try
  {
//  const IdToDetectorType *id_to_dettype_ptr
//               = boost::any_cast< const IdToDetectorType * >( id_to_dettypeany );
//  DetectorToCalibInfo *calibrations_ptr
//             = boost::any_cast< DetectorToCalibInfo *>( calibrations_any );
  
//  if( !id_to_dettype_ptr && !id_to_dettypeany.empty() )
//    throw logic_error( "decode_2012_N42_rad_measurment_node(...): invalid id_to_dettypeany" );
//    
//  if( !calibrations_ptr && !calibrations_any.empty() )
//    throw logic_error( "decode_2012_N42_rad_measurment_node(...): invalid calibrations_any" );
  
    vector<string> remarks;
    float real_time = 0.0;
    boost::posix_time::ptime start_time;
    Measurement::SourceType spectra_type = Measurement::UnknownSourceType;
    Measurement::OccupancyStatus occupied = Measurement::UnknownOccupancyStatus;
    
    rapidxml::xml_attribute<char> *meas_att = meas_node->first_attribute( "id", 2, false );
//    rapidxml::xml_attribute<char> *info_att = meas_node->first_attribute( "radItemInformationReferences", 28 );
//    rapidxml::xml_attribute<char> *group_att = meas_node->first_attribute( "radMeasurementGroupReferences", 29 );
  
    //Try to grab sample number from the 'id' attribute of <RadMeasurement>
    int sample_num_from_meas_node = -999;
    const string meas_id_att_str = xml_value_str( meas_att );
    if( meas_id_att_str.size() )
    {
      if( UtilityFunctions::icontains(meas_id_att_str,"background")
          && !UtilityFunctions::icontains(meas_id_att_str,"Survey")
          && !UtilityFunctions::icontains(meas_id_att_str,"Sample") )
      {
        sample_num_from_meas_node = 0;
      }else if( sscanf( meas_id_att_str.c_str(), "Survey %i", &(sample_num_from_meas_node)) == 1 )
      {
      }else if( sscanf( meas_id_att_str.c_str(), "Sample%i", &(sample_num_from_meas_node)) == 1 )
      {
      }//else ... another format I dont recall seeing.
    }//if( samp_det_str.size() )
    
    for( rapidxml::xml_node<char> *remark_node = meas_node->first_node( "Remark", 6 );
        remark_node;
        remark_node = XML_NEXT_TWIN(remark_node) )
    {
      string remark = xml_value_str( remark_node );
      trim( remark );
      if( remark.size() )
        remarks.push_back( remark );
    }//for( loop over remarks _
    
    rapidxml::xml_node<char> *class_code_node = meas_node->first_node( "MeasurementClassCode", 20 ); //XML_FIRST_NODE( meas_node, "MeasurementClassCode" )
    if( class_code_node && class_code_node->value_size() )
    {
      if( XML_VALUE_ICOMPARE(class_code_node, "Foreground") )
        spectra_type = Measurement::Foreground;
      else if( XML_VALUE_ICOMPARE(class_code_node, "Background") )
        spectra_type = Measurement::Background;
      else if( XML_VALUE_ICOMPARE(class_code_node, "Calibration") )
        spectra_type = Measurement::Calibration;
      else if( XML_VALUE_ICOMPARE(class_code_node, "IntrinsicActivity") )
        spectra_type = Measurement::IntrinsicActivity;
      else if( XML_VALUE_ICOMPARE(class_code_node, "NotSpecified") )
        spectra_type = Measurement::UnknownSourceType;
    }//if( class_code_node && class_code_node->value_size() )

    //Special check for RadSeeker.
    if( spectra_type == Measurement::UnknownSourceType
        && meas_att && XML_VALUE_ICOMPARE(meas_att, "Stabilization") )
      spectra_type = Measurement::IntrinsicActivity;
    
    rapidxml::xml_node<char> *time_node = meas_node->first_node( "StartDateTime", 13 );
    if( time_node && time_node->value_size() )
      start_time = time_from_string( xml_value_str(time_node).c_str() );
  
    rapidxml::xml_node<char> *real_time_node = meas_node->first_node( "RealTimeDuration", 16 );
    if( !real_time_node )
      real_time_node = meas_node->first_node( "RealTime", 8 );
    if( real_time_node && real_time_node->value_size() )
      real_time = time_duration_in_seconds( real_time_node->value(), real_time_node->value_size() );
  
    rapidxml::xml_node<char> *occupancy_node = meas_node->first_node( "OccupancyIndicator", 18 );
    if( occupancy_node && occupancy_node->value_size() )
    {
      if( XML_VALUE_ICOMPARE(occupancy_node, "true") )
        occupied = Measurement::Occupied;
      else if( XML_VALUE_ICOMPARE(occupancy_node, "false") )
        occupied = Measurement::NotOccupied;
    }//if( occupancy_node && occupancy_node->value_size() )

    vector< MeasurementShrdPtr > spectrum_meas, neutron_meas;
    
    //Lets track the Measurement to calibration id value incase multiple spectra
    //  are given for the same detector and <RadMeasurement>, but with different
    //  energy binnings.
    vector< pair<std::shared_ptr<Measurement>,string> > meas_to_cal_id;
    
    
    for( const rapidxml::xml_node<char> *spectrum_node = meas_node->first_node( "Spectrum", 8 );
         spectrum_node;
         spectrum_node = XML_NEXT_TWIN(spectrum_node) )
    {
      using rapidxml::internal::compare;
      rapidxml::xml_attribute<char> *id_att = spectrum_node->first_attribute( "id", 2, false );
      rapidxml::xml_attribute<char> *det_info_att = spectrum_node->first_attribute( "radDetectorInformationReference", 31, false );  //case sensitive false for refAO7WGOXDJ4
      rapidxml::xml_attribute<char> *calib_att = spectrum_node->first_attribute( "energyCalibrationReference", 26, false );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "fullEnergyPeakEfficiencyCalibrationReference", 44 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "FWHMCalibrationReference", 24 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "intrinsicDoubleEscapePeakEfficiencyCalibrationReference", 55 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "intrinsicFullEnergyPeakEfficiencyCalibrationReference", 53 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "intrinsicSingleEscapePeakEfficiencyCalibrationReference", 55 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "radRawSpectrumReferences", 24 );
    //      rapidxml::xml_attribute<char> *_att = spectrum_node->first_attribute( "totalEfficiencyCalibrationReference", 35 );

    
      std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
      DetectionType det_type = GammaDetection;
    
      //Get the detector name from the XML det_info_att if we have it, otherwise
      //  if there is only one detector description in the file, we will assume
      //  this spetrum is from that.
      if( det_info_att && det_info_att->value_size() )
        meas->detector_name_ = xml_value_str( det_info_att );
      else if( id_to_dettype_ptr && id_to_dettype_ptr->size()==1 )
        meas->detector_name_ = id_to_dettype_ptr->begin()->first;
    
      if( meas->detector_name_ == s_unnamed_det_placeholder )
        meas->detector_name_.clear();
    
      if( id_to_dettype_ptr )
      {
        map<string,pair<DetectionType,string> >::const_iterator det_iter
                              = id_to_dettype_ptr->find( meas->detector_name_ );
        if( det_iter != id_to_dettype_ptr->end() )
        {
          det_type = det_iter->second.first;
          meas->detector_type_ = det_iter->second.second; //e.x. "HPGe 50%"
        }//if( det_iter != id_to_dettype_ptr->end() )
      }//if( id_to_dettype_ptr )

      const rapidxml::xml_node<char> *live_time_node = spectrum_node->first_node( "LiveTimeDuration", 16 );
      if( !live_time_node )
        live_time_node = spectrum_node->first_node( "LiveTime", 8 );
      const rapidxml::xml_node<char> *channel_data_node = spectrum_node->first_node( "ChannelData", 11 );
    
      for( size_t i = 0; i < remarks.size(); ++i )
        meas->remarks_.push_back( remarks[i] );
    
      bool use_remark_real_time = false;
      
      for( rapidxml::xml_node<char> *remark_node = spectrum_node->first_node( "Remark", 6 );
           remark_node;
           remark_node = XML_NEXT_TWIN(remark_node) )
      {
        string remark = xml_value_str( remark_node );
        trim( remark );
        if( remark.empty() )
          continue;
          
        if( UtilityFunctions::istarts_with( remark, "RealTime:") )
        {
          //Starting with MeasurementInfo_2012N42_VERSION==3, a slightly more
          //  accurate RealTime may be recorded in the remark if necassary...
          //  see notes in create_2012_N42_xml() and add_spectra_to_measurment_node_in_2012_N42_xml()
          //snprintf( thisrealtime, sizeof(thisrealtime), "RealTime: PT%fS", realtime_used );
          remark = remark.substr( 9 );
          UtilityFunctions::trim( remark );
          meas->real_time_ = time_duration_in_seconds( remark );
          
          use_remark_real_time = (meas->real_time_ > 0.0);
        }else if( UtilityFunctions::istarts_with( remark, "Title:") )
        {
          //Starting with MeasurementInfo_2012N42_VERSION==3, title is encoded as a remark prepended with 'Title: '
          remark = remark.substr( 6 );
          UtilityFunctions::trim( remark );
          meas->title_ += remark;
        }else if( remark.size() )
        {
          meas->remarks_.push_back( remark );
        }
      }//for( loop over remarks )

    
      //This next line is specific to file written by InterSpec
      //const string samp_det_str = xml_value_str( meas_att ); //This was used pre 20180225, however I believe this was wrong due to it probably not containing DetXXX - we'll see.
      const string samp_det_str = xml_value_str( spectrum_node );
      if( samp_det_str.size() )
      {
        if( UtilityFunctions::istarts_with(samp_det_str, "background") )
        { 
          meas->sample_number_ = 0;
        }else if( sscanf( samp_det_str.c_str(), "Survey %i", &(meas->sample_number_)) == 1 )
        {
        }else if( sscanf( samp_det_str.c_str(), "Sample%i", &(meas->sample_number_)) == 1 )
        {
        }else if( sample_num_from_meas_node != -999 )
        {
          meas->sample_number_ = sample_num_from_meas_node;
        }
      }else if( sample_num_from_meas_node != -999 )
      {
        meas->sample_number_ = sample_num_from_meas_node;
      }//if( samp_det_str.size() )
    
#if(PERFORM_DEVELOPER_CHECKS)
      if( (sample_num_from_meas_node != -999) && (meas->sample_number_ != sample_num_from_meas_node) )
      {
        char buffer[512];
        snprintf( buffer, sizeof(buffer),
                  "Found a case where RadMeasurement id ('%s') gave a different"
                  " sample number than Spectrum id ('%s').",
                  meas_id_att_str.c_str(), samp_det_str.c_str() );
        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
      }
#endif
    
      if( !use_remark_real_time )
        meas->real_time_ = real_time;
      
      meas->start_time_ = start_time;
      meas->source_type_ = spectra_type;
      
      //For the sake of file_format_test_spectra/n42_2006/identiFINDER/20130228_184247Preliminary2010.n42
      if( meas->source_type_ == Measurement::UnknownSourceType
         && UtilityFunctions::iequals(meas->detector_name_, "intrinsicActivity")  )
        meas->source_type_ = Measurement::IntrinsicActivity;
      
      meas->occupied_ = occupied;
    
      if( live_time_node && live_time_node->value_size() )
        meas->live_time_ = time_duration_in_seconds( live_time_node->value(), live_time_node->value_size() );
    
      std::shared_ptr<vector<float> > gamma_counts = std::make_shared<vector<float> >();
    
      if( channel_data_node && channel_data_node->value_size() )
      {
        const char *char_data = channel_data_node->value();
        const size_t char_data_len = channel_data_node->value_size();
        UtilityFunctions::split_to_floats( char_data, char_data_len, *gamma_counts );
      
        rapidxml::xml_attribute<char> *comp_att = channel_data_node->first_attribute( "compressionCode", 15 );
        if( comp_att && XML_VALUE_ICOMPARE(comp_att, "CountedZeroes") )
          expand_counted_zeros( *gamma_counts, *gamma_counts );
      }//if( channel_data_node && channel_data_node->value() )

      //Well swapp meas->gamma_count_sum_ and meas->neutron_counts_sum_ in a bit if we need to
      meas->gamma_count_sum_ = 0.0;
      for( const float a : *gamma_counts )
        meas->gamma_count_sum_ += a;
    
      const rapidxml::xml_node<char> *RadItemState = meas_node->first_node( "RadItemState", 12 );
      const rapidxml::xml_node<char> *StateVector = (RadItemState ? RadItemState->first_node( "StateVector", 11 ) : (const rapidxml::xml_node<char> *)0);
      const rapidxml::xml_node<char> *SpeedValue = (StateVector ? StateVector->first_node( "SpeedValue", 10 ) : (const rapidxml::xml_node<char> *)0);
    
      if( SpeedValue && SpeedValue->value_size() )
      {
        const string val = xml_value_str( SpeedValue );
        if( !(stringstream(val) >> (meas->speed_)) )
          cerr << "Failed to convert '" << val << "' to a numeric speed" << endl;
      }//if( speed_ > 0 )
    
      if( det_type == OtherDetection )
        continue;

      bool is_gamma = (det_type == GammaDetection);
      if( det_type == GammaAndNeutronDetection )
      {
        const string att_val = xml_value_str( id_att );
        is_gamma = !icontains( att_val, "Neutron" );
        if( !calib_att || !calib_att->value_size() )
          is_gamma = false;
      }//if( det_type == GammaAndNeutronDetection )
    
      if( is_gamma )
      {
        meas->gamma_counts_ = gamma_counts;
      
        if( !calib_att || !calib_att->value_size() )
        {
          cerr << "Found a gamma spectrum without calibration "
               << "information, skipping it" << endl;
          continue;
        }//if( !calib_att || !calib_att->value_size() )
      
        if( !gamma_counts || gamma_counts->empty() )
        {
          cerr << "Found a gamma spectrum with no channel data, skipping" << endl;
          continue;
        }//if( !meas->gamma_counts_ || meas->gamma_counts_.empty() )

        if( calibrations_ptr )
        {//begin code block to deal with calibrations
          std::lock_guard<std::mutex> lock( calib_mutex );

          const string detnam = xml_value_str(calib_att);
        
          map<string,MeasurementCalibInfo>::iterator calib_iter
                                             = calibrations_ptr->find( detnam );
    
          if( calib_iter == calibrations_ptr->end() )
          {
            cerr << "No calibration information for gamma spcetrum claiming to "
                 << "have calibration data '" << detnam
                 << "' skipping" << endl;
            continue;
          }//if( calib_iter == calibrations_ptr->end() )
      
          MeasurementCalibInfo &calib = calib_iter->second;
      
          calib.nbin = meas->gamma_counts_->size();
          calib.fill_binning();
      
          if( !calib.binning )
          {
            cerr << "Calibration somehow invalid for '" << detnam
                 << "', skipping filling out." << endl;
            continue;
          }//if( !calib.binning )
          
          meas->calibration_coeffs_ = calib.coefficients;
          meas->deviation_pairs_    = calib.deviation_pairs_;
          meas->channel_energies_   = calib.binning;
          meas->energy_calibration_model_  = calib.equation_type;
          
          if( calib.calib_id.size() )
            meas_to_cal_id.push_back( make_pair(meas,calib.calib_id) );
        }//if( calibrations_ptr )

        meas->contained_neutron_ = false;
      }else
      {
        meas->neutron_counts_sum_ = meas->gamma_count_sum_;
        meas->gamma_count_sum_ = 0.0;
        meas->gamma_counts_.reset();
        meas->contained_neutron_ = true;
      //        if( gamma_counts )
      //          meas->neutron_counts_.swap( *gamma_counts );
      }//if( is_gamma ) / else
    //      const rapidxml::xml_node<char> *extension_node = meas_node->first_node( "SpectrumExtension", 17 );
    
      decode_2012_N42_detector_state_and_quality( meas, meas_node );
    
      spectrum_meas.push_back( meas );
    }//for( loop over "Spectrum" nodes )
    
    //flir radHUNTER N42 files is the inspiration for this next loop that
    // checks if there is a min, max, and total neutron <GrossCounts> node
    //  for this <RadMeasurement> node.
    bool min_neut = false, max_neut = false, total_neut = false, has_other = false;
    for( rapidxml::xml_node<char> *node = meas_node->first_node( "GrossCounts", 11 );
        node;
        node = XML_NEXT_TWIN(node) )
    {
      const rapidxml::xml_attribute<char> *att = node->first_attribute( "radDetectorInformationReference", 31, false );
      const bool is_min = XML_VALUE_ICOMPARE(att, "minimumNeutrons");
      const bool is_max = XML_VALUE_ICOMPARE(att, "maximumNeutrons");
      const bool is_total = XML_VALUE_ICOMPARE(att, "totalNeutrons");
      
      min_neut = (min_neut || is_min);
      max_neut = (max_neut || is_max);
      total_neut = (total_neut || is_total);
      has_other = (has_other || (!is_min && !is_max && !is_total));
    }//for( loop over GrossCounts nodes )
    
    const bool has_min_max_total_neutron = ((min_neut && max_neut && total_neut) && !has_other);
    
    for( rapidxml::xml_node<char> *gross_counts_node = meas_node->first_node( "GrossCounts", 11 );
         gross_counts_node;
         gross_counts_node = XML_NEXT_TWIN(gross_counts_node) )
    {
      const rapidxml::xml_node<char> *live_time_node = gross_counts_node->first_node( "LiveTimeDuration", 16 );
      const rapidxml::xml_node<char> *count_data_node = gross_counts_node->first_node( "CountData", 9 );
      const rapidxml::xml_attribute<char> *det_info_att = gross_counts_node->first_attribute( "radDetectorInformationReference", 31, false );
    
      std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
      DetectionType det_type = GammaDetection;
    
      const string det_info_ref = xml_value_str( det_info_att );
      
      if( det_info_ref.empty() )
      {
        cerr << "Found GrossCounts node with no radDetectorInformationReference"
             << endl;
        continue;
      }//if( det_info_ref.empty() )
    
      
      if( has_min_max_total_neutron )
      {
        if( !UtilityFunctions::iequals(det_info_ref, "totalNeutrons") )
          continue;
      }
      
      meas->detector_name_ = det_info_ref;
      if( meas->detector_name_ == s_unnamed_det_placeholder )
        meas->detector_name_.clear();
    
      if( id_to_dettype_ptr )
      {
        map<string,pair<DetectionType,string> >::const_iterator det_iter
                             = id_to_dettype_ptr->find( meas->detector_name_ );
        if( det_iter == id_to_dettype_ptr->end() )
        {
          cerr << "No detector information for '" << meas->detector_name_
               << "' so skipping" << endl;
          continue;
        }//if( !id_to_dettype_ptr->count( meas->detector_name_ ) )
      
        det_type = det_iter->second.first;
        meas->detector_type_ = det_iter->second.second; //e.x. "HPGe 50%"
      }//if( id_to_dettype_ptr )

      if( icontains( det_info_ref, "Neutrons" ) )
        det_type = NeutronDetection;
    
      if( det_type != NeutronDetection
          && (det_type != GammaAndNeutronDetection || !meas_node->first_node("Spectrum",8) ) )
      {
        cerr << "Found a non neutron GrossCount node, skipping" << endl;
        continue;
      }
    
      //This next line is specific to file written by InterSpec
      //const string sample_det_att = xml_value_str( meas_att ); //See notes above about pre 20180225,
      const string sample_det_att = xml_value_str( gross_counts_node );
      if( sample_det_att.size() )
      {
        if( UtilityFunctions::istarts_with(sample_det_att, "background") )
        {
          meas->sample_number_ = 0;
        }else if( sscanf( sample_det_att.c_str(), "Survey %i", &(meas->sample_number_) ) == 1 )
        {
        }else if( sscanf( sample_det_att.c_str(), "Sample%i", &(meas->sample_number_) ) == 1 )
        {
        }else if( sample_num_from_meas_node != -999 )
        {
          meas->sample_number_ = 0;
        }else
        {
#if(PERFORM_DEVELOPER_CHECKS)
          char buffer[256];
          snprintf( buffer, sizeof(buffer),
                   "Unrecognized 'id' attribute of Spectrum node: '%s'", sample_det_att.c_str() );
          log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
#endif
        }
      }//if( sample_det_att.size() )
    
      bool use_remark_real_time = false;
      
      for( size_t i = 0; i < remarks.size(); ++i )
        meas->remarks_.push_back( remarks[i] );
    
      for( rapidxml::xml_node<char> *remark_node = gross_counts_node->first_node( "Remark", 6 );
           remark_node;
           remark_node = XML_NEXT_TWIN(remark_node) )
      {
        string remark = xml_value_str( remark_node );
        trim( remark );
        if( remark.empty() )
          continue;
        
        if( UtilityFunctions::istarts_with( remark, "RealTime:") )
        {
          //See notes in equivalent portion of code for the <Spectrum> tag
          remark = remark.substr( 9 );
          UtilityFunctions::trim( remark );
          meas->real_time_ = time_duration_in_seconds( remark );
          use_remark_real_time = (meas->real_time_ > 0.0);
        }else if( UtilityFunctions::istarts_with( remark, "Title:") )
        {
          //See notes in equivalent portion of code for the <Spectrum> tag
          remark = remark.substr( 6 );
          UtilityFunctions::trim( remark );
          meas->title_ += remark;
        }else if( remark.size() )
        {
          meas->remarks_.push_back( remark );
        }
      }//for( loop over remark nodes )
    
      if( !use_remark_real_time )
        meas->real_time_ = real_time;
      meas->start_time_ = start_time;
      meas->source_type_ = spectra_type;
      meas->occupied_ = occupied;
    
      if( live_time_node && live_time_node->value_size() )
        meas->live_time_ = time_duration_in_seconds( live_time_node->value(), live_time_node->value_size() );
    
      
      const rapidxml::xml_node<char> *RadItemState = meas_node->first_node( "RadItemState", 12 );
      const rapidxml::xml_node<char> *StateVector = (RadItemState ? RadItemState->first_node( "StateVector", 11 ) : (const rapidxml::xml_node<char> *)0);
      const rapidxml::xml_node<char> *SpeedValue = (StateVector ? StateVector->first_node( "SpeedValue", 10 ) : (const rapidxml::xml_node<char> *)0);
      
      if( SpeedValue && SpeedValue->value_size() )
      {
        const string val = xml_value_str( SpeedValue );
        if( !(stringstream(val) >> (meas->speed_)) )
          cerr << "Failed to convert '" << val << "' to a numeric speed" << endl;
      }//if( speed_ > 0 )
      
      
      meas->contained_neutron_ = true;

      if( !count_data_node || !count_data_node->value_size() )
        count_data_node = gross_counts_node->first_node( "GrossCountData", 14 );
    
      if( !count_data_node || !count_data_node->value_size() )
      {
        cerr << "Found a GrossCount node without a CountData node, skipping" << endl;
        continue;
      }
    
      UtilityFunctions::split_to_floats( count_data_node->value(),
                                         count_data_node->value_size(),
                                         meas->neutron_counts_ );
      for( size_t i = 0; i < meas->neutron_counts_.size(); ++i )
        meas->neutron_counts_sum_ += meas->neutron_counts_[i];
    
      decode_2012_N42_detector_state_and_quality( meas, meas_node );
    
      neutron_meas.push_back( meas );
    }//for( loop over GrossCounts Node )
  
    vector<MeasurementShrdPtr> meas_to_add;
    if( spectrum_meas.size() == neutron_meas.size() )
    {
      for( size_t i = 0; i < spectrum_meas.size(); ++i )
      {
        MeasurementShrdPtr &gamma = spectrum_meas[i];
        MeasurementShrdPtr &neutron = neutron_meas[i];
      
        gamma->neutron_counts_ = neutron->neutron_counts_;
        gamma->neutron_counts_sum_ = neutron->neutron_counts_sum_;
        gamma->contained_neutron_ = neutron->contained_neutron_;
        for( const string &s : neutron->remarks_ )
        {
          if( std::find(gamma->remarks_.begin(), gamma->remarks_.end(), s)
              == gamma->remarks_.end() )
          gamma->remarks_.push_back( s );
        }
      
        meas_to_add.push_back( gamma );
      }//for( size_t i = 0; i < spectrum_meas.size(); ++i )
    }else
    {
      vector< pair<MeasurementShrdPtr,MeasurementShrdPtr> > gamma_and_neut_pairs;
      for( size_t i = 0; i < spectrum_meas.size(); ++i )
      {
        if( !spectrum_meas[i] )
          continue;
        
        for( size_t j = 0; j < neutron_meas.size(); ++j )
        {
          if( !neutron_meas[j] )
            continue;
          
          const string &gdetname = spectrum_meas[i]->detector_name_;
          const string &ndetname = neutron_meas[j]->detector_name_;
          if( gdetname.size() < 2 || ndetname.size() < 2 )
            continue;
          
          if( (gdetname == ndetname)
              || (istarts_with(ndetname, gdetname) && UtilityFunctions::icontains(ndetname, "neut")) )
          {
            gamma_and_neut_pairs.push_back( make_pair(spectrum_meas[i], neutron_meas[j]) );
            spectrum_meas[i].reset();
            neutron_meas[j].reset();
            break;
          }
        }//for( size_t j = 0; j < neutron_meas.size(); ++j )
      }//for( size_t i = 0; i < spectrum_meas.size(); ++i )
      
      for( size_t i = 0; i < gamma_and_neut_pairs.size(); ++i )
      {
        MeasurementShrdPtr &gamma = gamma_and_neut_pairs[i].first;
        MeasurementShrdPtr &neutron = gamma_and_neut_pairs[i].second;
        
        gamma->neutron_counts_ = neutron->neutron_counts_;
        gamma->neutron_counts_sum_ = neutron->neutron_counts_sum_;
        gamma->contained_neutron_ = neutron->contained_neutron_;
        for( const string &s : neutron->remarks_ )
        {
          if( std::find(gamma->remarks_.begin(), gamma->remarks_.end(), s)
             == gamma->remarks_.end() )
            gamma->remarks_.push_back( s );
        }
        
        meas_to_add.push_back( gamma );
      }//for( size_t i = 0; i < gamma_and_neut_pairs.size(); ++i )
      
      
      for( size_t i = 0; i < spectrum_meas.size(); ++i )
        if( spectrum_meas[i] )
          meas_to_add.push_back( spectrum_meas[i] );
    
      for( size_t i = 0; i < neutron_meas.size(); ++i )
        if( neutron_meas[i] )
          meas_to_add.push_back( neutron_meas[i] );
    }//
    
  //XXX - todo - should implement the bellow
  //    rapidxml::xml_node<char> *dose_rate_node = meas_node->first_node( "DoseRate", 8 );
  //    rapidxml::xml_node<char> *total_dose_node = meas_node->first_node( "TotalDose", 9 );
  //    rapidxml::xml_node<char> *exposure_rate_node = meas_node->first_node( "ExposureRate", 12 );
  //    rapidxml::xml_node<char> *total_exposure_node = meas_node->first_node( "TotalExposure", 13 );
  //    rapidxml::xml_node<char> *instrument_state_node = meas_node->first_node( "RadInstrumentState", 18 );
  //    rapidxml::xml_node<char> *item_state_node = meas_node->first_node( "RadItemState", 12 );
    
    //Check for duplicate spectrum in spectrum_meas for the same detector, but
    //  with different calibrations.
    //  See comments for #energy_cal_variants and #keep_energy_cal_variant.
    //Note: as of 20160531, this duplicate spectrum stuff is untested.
    const vector<MeasurementShrdPtr>::const_iterator beginmeas = meas_to_add.begin();
    const vector<MeasurementShrdPtr>::const_iterator endmeas = meas_to_add.end();
    for( size_t i = 1; i < meas_to_cal_id.size(); ++i )
    {
      MeasurementShrdPtr &meas = meas_to_cal_id[i].first;
      if( std::find(beginmeas,endmeas,meas) == endmeas )
        continue;
      
      vector< pair<std::shared_ptr<Measurement>,string> > samenames;
      for( size_t j = 0; j < i; ++j )
      {
        std::shared_ptr<Measurement> &innermeas = meas_to_cal_id[j].first;

        if( std::find(beginmeas,endmeas,innermeas) == endmeas )
          continue;
        
        if( innermeas->detector_name_ == meas->detector_name_
           && innermeas->start_time_ == meas->start_time_
           && fabs(innermeas->real_time_ - meas->real_time_) < 0.001
           && fabs(innermeas->live_time_ - meas->live_time_) < 0.001 )
        {
          samenames.push_back( make_pair(innermeas, meas_to_cal_id[j].second ) );
        }
      }//for( size_t j = 0; j < i; ++j )
      
      if( samenames.size() )
      {
        meas->detector_name_ += "_intercal_" + meas_to_cal_id[i].second;
        for( size_t j = 0; j < samenames.size(); ++j )
          samenames[j].first->detector_name_ += "_intercal_" + samenames[j].second;
      }//if( samenames.size() )
    }//for( size_t i = 1; i < meas_to_cal_id.size(); ++i )
    
    {
      std::lock_guard<std::mutex> lock( meas_mutex );
      measurments.insert( measurments.end(), meas_to_add.begin(), meas_to_add.end() );
    }
  }catch( std::exception &e )
  {
    std::lock_guard<std::mutex> lock( meas_mutex );
    cerr << "Error decoding MeasurementInfo::decode2012N42SpectrumNode(...): "
         << e.what() << endl;
  }//try / catch
}//void MeasurementInfo::decode2012N42SpectrumNode( const rapidxml::xml_node<char> *spectrum )



void MeasurementInfo::load_2012_N42_from_doc( const rapidxml::xml_node<char> *data_node )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !data_node )
    throw runtime_error( "load_2012_N42_from_doc: Invalid document node" );
  
  if( !XML_NAME_ICOMPARE(data_node, "RadInstrumentData") )
    throw runtime_error( "load_2012_N42_from_doc: Unable to get RadInstrumentData node" );
  
  rapidxml::xml_attribute<char> *uuid_att = data_node->first_attribute( "n42DocUUID", 10 );
  if( uuid_att && uuid_att->value_size() )
    uuid_ = xml_value_str(uuid_att);
  
  /*
   If the N42 XML document has been created by translation from an original data file,
   include the name of the software that was used and the name and format type of the original data file.
   
   This information is kinda duplicate information of <RadInstrumentVersion><RadInstrumentComponentName>Software</>...</>
   
   */
  rapidxml::xml_node<char> *creator_node = data_node->first_node( "RadInstrumentDataCreatorName", 28 );
  if( creator_node && creator_node->value_size() )
    remarks_.push_back( "N42 file created by: " + xml_value_str(creator_node) );
  
  for( rapidxml::xml_node<char> *remark_node = data_node->first_node( "Remark", 6 );
      remark_node;
      remark_node = XML_NEXT_TWIN(remark_node) )
  {
    string remark = xml_value_str( remark_node );
    trim( remark );
    if( remark.size() )
      remarks_.push_back( remark );
  }
  
  rapidxml::xml_node<char> *inst_info_node = data_node->first_node( "RadInstrumentInformation", 24 );
  set_2012_N42_instrument_info( inst_info_node );
  
  map<string,MeasurementCalibInfo> calibrations;
  get_2012_N42_energy_calibrations( calibrations, data_node, remarks_ );
  
//XXX - implement using RadItemInformation
//  for( const rapidxml::xml_node<char> *rad_item_node = data_node->first_node( "RadItemInformation", 18 );
//      rad_item_node;
//      rad_item_node = XML_NEXT_TWIN(rad_item_node) )
//  {
//    rapidxml::xml_attribute<char> *id_att = rad_item_node->first_attribute( "id", 2, false );
//    rapidxml::xml_node<char> *remark_node = rad_item_node->first_node( "Remark", 6 );
//    rapidxml::xml_node<char> *descrip_node = rad_item_node->first_node( "RadItemDescription", 18 );
//    rapidxml::xml_node<char> *quantity_node = rad_item_node->first_node( "RadItemQuantity", 15 );
//    rapidxml::xml_node<char> *geometry_node = rad_item_node->first_node( "RadItemMeasurementGeometryDescription", 37 );
//    rapidxml::xml_node<char> *characteristics_node = rad_item_node->first_node( "RadItemCharacteristics", 22 );
//    rapidxml::xml_node<char> *extension_node = rad_item_node->first_node( "RadItemInformationExtension", 27 );
//  }//for( loop over "RadItemInformation" nodes )
  
//  map<string,int> detname_to_num; //only used when reading in InterSpec generated files
  map<string,pair<DetectionType,string> > id_to_dettype;
  
  for( const rapidxml::xml_node<char> *info_node = data_node->first_node( "RadDetectorInformation", 22 );
      info_node;
      info_node = XML_NEXT_TWIN(info_node) )
  {
    rapidxml::xml_attribute<char> *id_att   = info_node->first_attribute( "id", 2, false );
//    rapidxml::xml_node<char> *remark_node   = info_node->first_node( "Remark", 6 );
    rapidxml::xml_node<char> *name_node     = info_node->first_node( "RadDetectorName", 15 );
    rapidxml::xml_node<char> *category_node = info_node->first_node( "RadDetectorCategoryCode", 23 );
    
    //<RadDetectorKindCode> returns "NaI", "HPGe", etc (see determine_rad_detector_kind_code())
    //  and should be utilized at some point.  But would require adding a field to MeasurmentInfo
    //  I think to kind of do it properly.
    
    //rapidxml::xml_node<char> *kind_node     = info_node->first_node( "RadDetectorKindCode", 19 );
    rapidxml::xml_node<char> *descrip_node  = info_node->first_node( "RadDetectorDescription", 22 );
    rapidxml::xml_node<char> *length_node   = info_node->first_node( "RadDetectorLengthValue", 22 );
    rapidxml::xml_node<char> *width_node    = info_node->first_node( "RadDetectorWidthValue", 21 );
    rapidxml::xml_node<char> *depth_node    = info_node->first_node( "RadDetectorDepthValue", 21 );
    rapidxml::xml_node<char> *diameter_node = info_node->first_node( "RadDetectorDiameterValue", 24 );
    rapidxml::xml_node<char> *volume_node   = info_node->first_node( "RadDetectorVolumeValue", 22 );
//    rapidxml::xml_node<char> *characteristics_node = info_node->first_node( "RadDetectorCharacteristics", 26 );
//    rapidxml::xml_node<char> *extension_node = info_node->first_node( "RadDetectorInformationExtension", 31 );
    
    string name = xml_value_str(id_att);
    if( name == s_unnamed_det_placeholder )
    {
      name.clear();
    }else
    {
      if( name.empty() )
        name = xml_value_str(name_node);
    
      if( name.empty() )
      {
        rapidxml::xml_attribute<char> *ref_att = info_node->first_attribute( "Reference", 9, false );
        name = xml_value_str(ref_att);
      }
    }//if( name == s_unnamed_det_placeholder ) / else
    
//    rapidxml::xml_node<char> *detinfo_extension_node = info_node->first_node( "RadDetectorInformationExtension", 31 );
//    rapidxml::xml_node<char> *det_num_node = detinfo_extension_node ? detinfo_extension_node->first_node( "InterSpec:DetectorNumber", 24 ) : (rapidxml::xml_node<char> *)0;
//    if( det_num_node && det_num_node->value() )
//    {
//      int detnum;
//      if( sscanf(det_num_node->value(), "%d", &detnum) == 1 )
//        detname_to_num[name] = detnum;
//    }
    
    DetectionType type = GammaDetection; //OtherDetection;
    if( category_node && category_node->value_size() )
    {
      if( XML_VALUE_ICOMPARE(category_node, "Gamma") )
        type = GammaDetection;
      else if( XML_VALUE_ICOMPARE(category_node, "Neutron") )
        type = NeutronDetection;
      else
        type = OtherDetection;
      
      const string desc = xml_value_str( descrip_node );
      if( icontains( desc, "Gamma" ) && icontains( desc, "Neutron" ) )
         type = GammaAndNeutronDetection;
    }//if( category_node && category_node->value_size() )
    
    string descrip; // = xml_value_str( kind_node );
    
    descrip = xml_value_str( descrip_node );
    UtilityFunctions::ireplace_all( descrip, ", Gamma and Neutron", "" );
    UtilityFunctions::ireplace_all( descrip, "Gamma and Neutron", "" );
  
    if( length_node && length_node->value_size() )
    {
      if( descrip.length() )
        descrip += ", ";
      descrip += string("Length: ") + xml_value_str(length_node) + string(" cm");
    }
    
    if( width_node && width_node->value_size() )
    {
      if( descrip.length() )
        descrip += ", ";
      descrip += string("Width: ") + xml_value_str(width_node) + string(" cm");
    }
    
    if( depth_node && depth_node->value_size() )
    {
      if( descrip.length() )
        descrip += ", ";
      descrip += string("Depth: ") + xml_value_str(depth_node) + string(" cm");
    }

    if( diameter_node && diameter_node->value_size() )
    {
      if( descrip.length() )
        descrip += ", ";
      descrip += string("Diameter: ") + xml_value_str(diameter_node) + string(" cm");
    }

    if( volume_node && volume_node->value_size() )
    {
      if( descrip.length() )
        descrip += ", ";
      descrip += string("Volume: ") + xml_value_str(volume_node) + string(" cc");
    }
    
    if( type==GammaDetection || type==NeutronDetection || type==GammaAndNeutronDetection )
      id_to_dettype[name] = pair<DetectionType,string>(type,descrip);
  }//for( loop over "RadDetectorInformation" nodes )

  
  
  rapidxml::xml_node<char> *analysis_node = data_node->first_node( "AnalysisResults", 15 );
  if( analysis_node )
  {
    std::shared_ptr<DetectorAnalysis> analysis_info = std::make_shared<DetectorAnalysis>();
    setAnalysisInformation( analysis_node, analysis_info );
//    if( analysis_info->results_.size() )
      detectors_analysis_ = analysis_info;
  }//if( analysis_node )
  
  SpecUtilsAsync::ThreadPool workerpool;
  
//  Need to make order in measurements_ reproducable
  vector<std::shared_ptr<std::mutex> > meas_mutexs;
  vector< std::shared_ptr< vector<std::shared_ptr<Measurement> > > > measurements_each_meas;
  
  std::mutex meas_mutex, calib_mutex;
  //vector< boost::function<void()> > rad_meas_queue;
  
  for( const rapidxml::xml_node<char> *meas_node = data_node->first_node( "RadMeasurement", 14 );
      meas_node;
      meas_node = XML_NEXT_TWIN(meas_node) )
  {
    std::shared_ptr<std::mutex> mutexptr = std::make_shared<std::mutex>();
    std::shared_ptr< vector<std::shared_ptr<Measurement> > > these_meas = std::make_shared< vector<std::shared_ptr<Measurement> > >();
    
    meas_mutexs.push_back( mutexptr );
    measurements_each_meas.push_back( these_meas );
    
    //auto worker = [these_meas,meas_node,&id_to_dettype,&calibrations,mutexptr,&calib_mutex](){
      //MeasurementInfo::decode_2012_N42_rad_measurment_node( &these_meas, meas_node, &id_to_dettype, &calibrations, *mutexptr, calib_mutex );
    //};
    //workerpool.post( worker );
    workerpool.post( std::bind( &MeasurementInfo::decode_2012_N42_rad_measurment_node,
                std::ref(*these_meas),
                meas_node,
                &id_to_dettype,
                &calibrations,
                std::ref(*mutexptr),
                std::ref(calib_mutex) ) );
  }//for( loop over "RadMeasurement" nodes )

  workerpool.join();

  for( size_t i = 0; i < measurements_each_meas.size(); ++i )
    for( size_t j = 0; j < measurements_each_meas[i]->size(); ++j )
      measurements_.push_back( (*measurements_each_meas[i])[j] );
  
  
  //test for files like "file_format_test_spectra/n42_2006/identiFINDER/20130228_184247Preliminary2010.n42"
  if( measurements_.size() == 2
     && inst_info_node && inst_info_node->first_node( "RadInstrumentModel", 18 ) )
  {
    bool has_spectra = (measurements_[0]->detector_name_ == "spectra");
    bool has_intrinsic = (measurements_[0]->detector_name_ == "intrinsicActivity");
    
    has_spectra |= (measurements_[1]->detector_name_ == "spectra");
    has_intrinsic |= (measurements_[1]->detector_name_ == "intrinsicActivity");
    
    if( has_spectra && has_intrinsic )
    {
      detector_names_.clear();
      neutron_detector_names_.clear();
      measurements_[0]->detector_name_ = measurements_[1]->detector_name_ = "";
    }
  }//if( measurements_.size() == 2 )
  
  
/*
    XXX - Elements still not addressed/implemented
    <RadMeasurementGroup>   //FLIR has
    <FWHMCalibration>
    <TotalEfficiencyCalibration>
    <FullEnergyPeakEfficiencyCalibration>
    <IntrinsicFullEnergyPeakEfficiencyCalibration>
    <IntrinsicSingleEscapePeakEfficiencyCalibration>
    <IntrinsicDoubleEscapePeakEfficiencyCalibration>
    <EnergyWindows>
    <DerivedData>
    <AnalysisResults>     //FLIR has
    <MultimediaData>
  <RadInstrumentDataExtension>
*/

  if( measurements_.empty() )
    throw runtime_error( "No valid measurments in 2012 N42 file." );
  
  cleanup_after_load();  //DontChangeOrReorderSamples
}//bool load_2012_N42_from_doc( rapidxml::xml_node<char> *document_node )


bool MeasurementInfo::load_from_N42_document( const rapidxml::xml_node<char> *document_node )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !document_node || !document_node->name_size() )
    throw runtime_error( "no first node" );
  
  const string doc_node_name = xml_name_str(document_node);
  if( doc_node_name == "RadInstrumentData" )
  {
    load_2012_N42_from_doc( document_node );
  }else if( doc_node_name == "Event" )
  {
    //HPRD files
    const rapidxml::xml_node<char> *daughter = document_node->first_node( "N42InstrumentData", 17 );
    
    if( !daughter )
      throw runtime_error( "Unrecognized N42 file structure" );
    
    load_2006_N42_from_doc( daughter );
    
    bool hprds = false;
    const rapidxml::xml_node<char> *dataformat = document_node->first_node( "ThisDataFormat", 14 );
    if( dataformat && dataformat->value_size() )
      hprds = UtilityFunctions::icontains( xml_value_str(dataformat), "HPRDS" );
    
    if( hprds )
    {
      const rapidxml::xml_node<char> *node = document_node->first_node( "OnsetDateTime", 13 );
      node = document_node->first_node( "EventCategory", 13 );
      if( node && node->value_size() )
        remarks_.push_back( string("Event Category ") + xml_value_str(node) );
      
      node = document_node->first_node( "EventType", 9 );
      if( node && node->value_size() )
        remarks_.push_back( string("Event Type ") + xml_value_str(node) );
      
      node = document_node->first_node( "EventCode", 9 );
      if( node && node->value_size() )
        remarks_.push_back( string("Event Code ") + xml_value_str(node) );
      
      node = document_node->first_node( "EventNumber", 11 );
      if( node && node->value_size() )
        remarks_.push_back( string("Event Number ") + xml_value_str(node) );
    
      if( measurements_.size() == 2 )
      {
        std::shared_ptr<Measurement> gamma, neutron;
        
        for( int i = 0; i < 2; ++i )
        {
          const string &dettype = measurements_[i]->detector_type_;
          if( UtilityFunctions::icontains( dettype, "Gamma" ) )
            gamma = measurements_[i];
          else if( UtilityFunctions::icontains( dettype, "Neutron" ) )
            neutron = measurements_[i];
        }//for( int i = 0; i < 2; ++i )
        
        if( gamma && neutron && gamma != neutron
            && neutron->num_gamma_channels() < 2 )
        {
          gamma->neutron_counts_     = neutron->neutron_counts_;
          gamma->neutron_counts_sum_ = neutron->neutron_counts_sum_;
          gamma->contained_neutron_  = neutron->contained_neutron_;
          
          measurements_.erase( find(measurements_.begin(),measurements_.end(),neutron) );
          
          //cleanup_after_load() wont refresh neutron_detector_names_ unless
          //  its empty
          neutron_detector_names_.clear();
          
          cleanup_after_load();
        }//if( gamma && neutron && gamma != neutron )
      }else if( measurements_.size() > 10 )  //10 chosen arbitrarily - meant to seperate from files that have hundred of samples we dont want
      {
        //We want to keep only
        set<int> keepersamples;
        std::vector< std::shared_ptr<Measurement> > keepers;
        
        for( const auto &m : measurements_ )
        {
          bool keep = false;
          if( m->source_type() == Measurement::Background )
            keep = true;
          for( const std::string &c : m->remarks_ )
            keep |= UtilityFunctions::icontains( c, "count" );
          if( keep )
          {
            keepersamples.insert( m->sample_number_ );
            keepers.push_back( m );
          }
        }//for( const auto &m : measurements_ )
      
      
        {//begin codeblock to filter remarks
          vector<string> remarks_to_keep;
          remarks_to_keep.reserve( remarks_.size() );
          for( const string &remark : remarks_ )
          {
            if( !UtilityFunctions::icontains(remark, "DNDORadiationMeasurement") )
              remarks_to_keep.push_back( remark );
          }
          remarks_.swap( remarks_to_keep );
        }//end codeblock to filter remarks
        
        const vector<int> samples( keepersamples.begin(), keepersamples.end() );
        const vector<int>::const_iterator sbegin = samples.begin();
        const vector<int>::const_iterator send = samples.begin();
        
        if( keepers.size() > 1 )
        {
          typedef map< boost::posix_time::ptime, vector<std::shared_ptr<Measurement> > > time_to_meas_t;
          time_to_meas_t time_to_meas;
          
          for( size_t i = 0; i < keepers.size(); ++i )
          {
            MeasurementShrdPtr &m = keepers[i];
            const int oldsn = m->sample_number_;
            const vector<int>::const_iterator pos = std::find(sbegin, send, oldsn);
            m->sample_number_ = 1 + static_cast<int>(pos - sbegin);
            time_to_meas[m->start_time_].push_back( m );
          }
          
          //ORTEC MicroDetective2 HPRDS files will have a gamma, nuetron, and
          //  GMTube spectrum under their <DetectorData> node, so lets try to
          //  combine gamma and neutron together.
          for( time_to_meas_t::value_type &vt : time_to_meas )
          {
            vector<std::shared_ptr<Measurement> > &meass = vt.second;
            
            size_t ngamma = 0, nneut = 0, ngm = 0;
            std::shared_ptr<Measurement> gm_meas, gamma_meas, neut_meas;
            for( size_t i = 0; i < meass.size(); ++i )
            {
              MeasurementShrdPtr &m = meass[i];
              
              if( UtilityFunctions::iequals(m->detector_name_, "gamma" ) )
              {
                ++ngamma;
                gamma_meas = m;
              }else if( UtilityFunctions::iequals(m->detector_name_, "neutron" ) )
              {
                ++nneut;
                neut_meas = m;
              }else if( UtilityFunctions::iequals(m->detector_name_, "GMTube" ) )
              {
                ++ngm;
                gm_meas = m;
              }
            }//for( size_t i = 0; i < meass.size(); ++i )
            
            if( ngamma == 1 && nneut == 1
               && neut_meas->gamma_count_sum_ < 1.0
               && (!gamma_meas->contained_neutron_ || gamma_meas->neutron_counts_sum_ < 1.0) )
            {
              gamma_meas->neutron_counts_     = neut_meas->neutron_counts_;
              gamma_meas->neutron_counts_sum_ = neut_meas->neutron_counts_sum_;
              gamma_meas->contained_neutron_  = neut_meas->contained_neutron_;
              
              //cleanup_after_load() wont refresh neutron_detector_names_ unless
              //  its empty
              neutron_detector_names_.clear();
              
              keepers.erase( find( keepers.begin(), keepers.end(), neut_meas ) );
              if( ngm == 1 )
              {
                //We could add a comment here about the GMTube counts...
                keepers.erase( find( keepers.begin(), keepers.end(), gm_meas ) );
              }
            }//if( we have one neutron and one gamma measurment with same start time )
          }//for( time_to_meas_t::value_type &vt : time_to_meas )
        
        
          measurements_ = keepers;
          cleanup_after_load();
        }//if( keepers.size() > 1 )
      }//if( measurements_.size() == 2 ) / else if( measurements_.size() > 2 )
      
    }//if( hprds )
  }else
  {
    rapidxml::xml_node<char> *daughter = document_node->first_node();
    
    if( daughter && daughter->first_node( "Measurement", 11 ) )
      document_node = daughter;
    
    load_2006_N42_from_doc( document_node );
  }//if( N42.42 2012 ) else if( variation on N42.42 2006 )
  
  return true;
}//bool load_from_N42_document( rapidxml::xml_node<char> *document_node )


bool MeasurementInfo::load_micro_raider_file( const std::string &filename )
{
  ifstream input( filename.c_str(), ios_base::binary|ios_base::in );
  
  if( !input.is_open() )
    return false;
  
  try
  {
    reset();
    ::rapidxml::file<char> input_file( input );
    const bool success = load_from_micro_raider_from_data( input_file.data() );
    
    if( success )
      filename_ = filename;
    
    return success;
  }catch( std::exception & )
  {
    reset();
    return false;
  }//try/catch
  
  return false;
}//bool load_micro_raider_file(...)


bool MeasurementInfo::load_from_micro_raider_from_data( const char *data )
{
  try
  {
    typedef char XmlChar;
    rapidxml::xml_document<XmlChar> doc;
    
    //Casting 'data' to non-const, BUT rapidxml::parse_non_destructive
    //  _should_ garuntee the source is no altered.  It's a little bit shaky
    //  but it allows us to potentially use mmap
    doc.parse<rapidxml::parse_non_destructive |rapidxml::allow_sloppy_parse>( (XmlChar *)data );
    const rapidxml::xml_node<XmlChar> *IdResult = XML_FIRST_NODE((&doc),"IdResult");
    
    if( !IdResult )
      throw runtime_error( "Invalid Micro Raider XML document" );
    
    const rapidxml::xml_node<XmlChar> *DeviceId = XML_FIRST_NODE(IdResult,"DeviceId");
    const rapidxml::xml_node<XmlChar> *SurveyId = XML_FIRST_NODE(IdResult,"SurveyId");
    const rapidxml::xml_node<XmlChar> *UUID = XML_FIRST_NODE(IdResult,"UUID");
    const rapidxml::xml_node<XmlChar> *EventNumber = XML_FIRST_NODE(IdResult,"EventNumber");
    const rapidxml::xml_node<XmlChar> *CrystalType = XML_FIRST_NODE(IdResult,"CrystalType");
    const rapidxml::xml_node<XmlChar> *UserMode = XML_FIRST_NODE(IdResult,"UserMode");
    const rapidxml::xml_node<XmlChar> *StartTime = XML_FIRST_NODE(IdResult,"StartTime");
//    const rapidxml::xml_node<XmlChar> *StopTime = XML_FIRST_NODE(IdResult,"StopTime");
    const rapidxml::xml_node<XmlChar> *GPS = XML_FIRST_NODE(IdResult,"GPS");
    const rapidxml::xml_node<XmlChar> *RealTime = XML_FIRST_NODE(IdResult,"RealTime");
    const rapidxml::xml_node<XmlChar> *LiveTime = XML_FIRST_NODE(IdResult,"LiveTime");
//    const rapidxml::xml_node<XmlChar> *DeadTime = XML_FIRST_NODE(IdResult,"DeadTime");
    const rapidxml::xml_node<XmlChar> *DoseRate = XML_FIRST_NODE(IdResult,"DoseRate");
//    const rapidxml::xml_node<XmlChar> *CountRate = XML_FIRST_NODE(IdResult,"CountRate");
    const rapidxml::xml_node<XmlChar> *NeutronCountRate = XML_FIRST_NODE(IdResult,"NeutronCountRate");
    
    
    const rapidxml::xml_node<XmlChar> *Nuclide = XML_FIRST_NODE(IdResult,"Nuclide");
    const rapidxml::xml_node<XmlChar> *Image = XML_FIRST_NODE(IdResult,"Image");
    const rapidxml::xml_node<XmlChar> *VoiceRecording = XML_FIRST_NODE(IdResult,"VoiceRecording");
    const rapidxml::xml_node<XmlChar> *Spectrum = XML_FIRST_NODE(IdResult,"Spectrum");
    
    if( !Spectrum || !Spectrum->value_size() )
      throw runtime_error( "No Spectrum Node" );
    
    std::shared_ptr< vector<float> > channel_counts
                                         = std::make_shared<vector<float> >();
    
    const bool validchannel
               = UtilityFunctions::split_to_floats( Spectrum->value(),
                                      Spectrum->value_size(), *channel_counts );
    if( !validchannel || channel_counts->empty() )
      throw runtime_error( "Couldnt parse channel counts" );
    
    MeasurementShrdPtr meas = std::make_shared<Measurement>();
    
    meas->gamma_counts_ = channel_counts;
    
    meas->gamma_count_sum_ = 0.0;
    const size_t nchannel = meas->gamma_counts_->size();
    for( size_t i = 0; i < nchannel; ++i )
      meas->gamma_count_sum_ += (*channel_counts)[i];
    
    instrument_id_ = xml_value_str( DeviceId );
    if( SurveyId )
      remarks_.push_back( "Survey ID: " + xml_value_str( SurveyId ) );
    uuid_ = xml_value_str( UUID );
    if( EventNumber )
      remarks_.push_back( "EventNumber: " + xml_value_str( EventNumber ) );
    if( CrystalType )
      remarks_.push_back( "CrystalType: " + xml_value_str( CrystalType ) );
    if( UserMode )
      remarks_.push_back( "CrystalType: " + xml_value_str( UserMode ) );
    
    //Unecasary allocation to get time.
    const string start_time = xml_value_str(StartTime);
    meas->start_time_ = UtilityFunctions::time_from_string( start_time.c_str() );

    if( GPS && GPS->value_size() )
    {
      rapidxml::xml_attribute<XmlChar> *att = GPS->first_attribute("Valid",5);
      
      if( !att || XML_VALUE_ICOMPARE(att,"True") )
        parse_deg_min_sec_lat_lon(GPS->value(), GPS->value_size(),
                                  meas->latitude_, meas->longitude_ );
    }//if( GPS )
    
    if( RealTime && RealTime->value_size() )
      meas->real_time_ = time_duration_in_seconds( RealTime->value(), RealTime->value_size() );
    if( LiveTime && LiveTime->value_size() )
      meas->live_time_ = time_duration_in_seconds( LiveTime->value(), LiveTime->value_size() );
    
    
    std::shared_ptr<DetectorAnalysis> detana;
    
    if( DoseRate && DoseRate->value_size() )
    {
      DetectorAnalysisResult detanares;
      const float doseunit = dose_units_usvPerH( DoseRate->value(),
                                                 DoseRate->value_size() );
      
      //Avoidable allocation here bellow
      float dose;
      if( (stringstream(xml_value_str(DoseRate)) >> dose) )
        detanares.dose_rate_ = dose * doseunit;
      else
        cerr << "Failed to turn '" << xml_value_str(DoseRate) << "' into a dose" << endl;
      
      //detanares.start_time_ = meas->start_time_;
      detanares.real_time_ = meas->real_time_;
      
      if( !detana )
        detana = std::make_shared<DetectorAnalysis>();
      detana->results_.push_back( detanares );
    }//if( DoseRate && DoseRate->value_size() )
    
    if( NeutronCountRate && NeutronCountRate->value_size() )
    {
      const string neutroncountstr = xml_value_str(NeutronCountRate);
      meas->neutron_counts_.resize( 1, 0.0f );
      
      float neutrons = 0.0f;
      if( toFloat(neutroncountstr,neutrons) )
      {
        if( meas->real_time_ > 0.0f )
          neutrons *= meas->real_time_;
        else if( meas->live_time_ > 0.0f )
          neutrons *= meas->live_time_;
        else
          meas->remarks_.push_back( "NeutronCountRate: " + neutroncountstr + " (error computing gross count)" ); //meh, should be fine...
        
        meas->neutron_counts_sum_ = neutrons;
        meas->neutron_counts_[0] = neutrons;
        meas->contained_neutron_ = true;
      }else
      {
        meas->remarks_.push_back( "NeutronCountRate: " + neutroncountstr );
        cerr << "Failed to read '" << neutroncountstr << "' as neutroncountstr"
             << endl;
      }
    }//if( NeutronCountRate && NeutronCountRate->value_size() )
    
    while( Nuclide )
    {
      DetectorAnalysisResult res;
      
      const rapidxml::xml_node<XmlChar> *NuclideName = XML_FIRST_NODE(Nuclide,"NuclideName");
      const rapidxml::xml_node<XmlChar> *NuclideType = XML_FIRST_NODE(Nuclide,"NuclideType");
      const rapidxml::xml_node<XmlChar> *NuclideIDConfidenceIndication = XML_FIRST_NODE(Nuclide,"NuclideIDConfidenceIndication");
      const rapidxml::xml_node<XmlChar> *NuclideIDStrengthIndication = XML_FIRST_NODE(Nuclide,"NuclideIDStrengthIndication");
      const rapidxml::xml_node<XmlChar> *NuclideDescription = XML_FIRST_NODE(Nuclide,"NuclideDescription");
//      const rapidxml::xml_node<XmlChar> *NuclideInstruction = XML_FIRST_NODE(Nuclide,"NuclideInstruction");
      const rapidxml::xml_node<XmlChar> *NuclideHPRDSType = XML_FIRST_NODE(Nuclide,"NuclideHPRDSType");
    
      res.nuclide_ = xml_value_str(NuclideName);
      res.nuclide_type_ = xml_value_str(NuclideType);
      res.id_confidence_ = xml_value_str(NuclideIDConfidenceIndication);
      
      const string strength = xml_value_str(NuclideIDStrengthIndication);
      if( strength.size() )
        res.remark_ += "strength: " + strength;
      if( NuclideHPRDSType && NuclideHPRDSType->value_size() )
        res.remark_ += (res.remark_.size() ? ". " : "") + xml_value_str(NuclideHPRDSType);
      if( NuclideDescription && NuclideDescription->value_size() )
        res.remark_ += (res.remark_.size() ? ". " : "") + xml_value_str(NuclideDescription);
      
      Nuclide = XML_NEXT_TWIN(Nuclide);
      
      if( !detana )
        detana = std::make_shared<DetectorAnalysis>();
      detana->results_.push_back( res );
    }//while( Nuclide )
    
    while( Image && Image->value_size() )
    {
      remarks_.push_back( "Image: " + xml_value_str(Image) );
      Image = XML_NEXT_TWIN(Image);
    }
    
    if( VoiceRecording && VoiceRecording->value_size() )
      remarks_.push_back( "VoiceRecording: " + xml_value_str(VoiceRecording) );
    
    detectors_analysis_ = detana;
    
    //Following values taken from a Micro Raider ICD1 N42 2006 file
    manufacturer_ = "ICx Radiation";
    instrument_model_ = "Raider";
    instrument_type_ = "Radionuclide Identifier";  //or PersonalRadiationDetector
    detector_type_ = kMicroRaiderDetector;
    
    measurements_.push_back( meas );
    
    cleanup_after_load();
    
    return true;
  }catch( std::exception & )
  {
    return false;
  }
  
  
//  "<IdResult"
  
  return false;
}//bool load_from_micro_raider_from_data( const char *data )


bool MeasurementInfo::load_from_N42( std::istream &input )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  
  try
  {
    ::rapidxml::file<char> input_file( input );
    return MeasurementInfo::load_N42_from_data( input_file.data() );
  }catch( std::exception & )
  {
    input.clear();
    input.seekg( orig_pos, ios::beg );
    reset();
    return false;
  }//try/catch
  
  return true;
}//bool load_from_N42( char *data )


bool MeasurementInfo::load_N42_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    rapidxml::file<char> input_file( filename.c_str() );  //throws runtime_error upon failure
    const bool loaded = MeasurementInfo::load_N42_from_data( input_file.data() );
    
    if( !loaded )
      throw runtime_error( "Failed to load" );
    
    filename_ = filename;
  }catch(...)
  {
    reset();
    return false;
  }
  
  return true;
}//bool load_N42_file( const std::string &filename );


bool MeasurementInfo::load_N42_from_data( char *data )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    reset();
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace|rapidxml::allow_sloppy_parse>( data );
    rapidxml::xml_node<char> *document_node = doc.first_node();
    
    const bool loaded = load_from_N42_document( document_node );
    
    if( !loaded )
      throw runtime_error( "Failed to load" );
  }catch(...)
  {
    reset();
    return false;
  }
  
  return true;
}//bool load_N42_from_data( char *data )




void MeasurementInfo::rebin_by_eqn( const std::vector<float> &eqn,
                                    const DeviationPairVec &dev_pairs,
                                    Measurement::EquationType type )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  ShrdConstFVecPtr new_binning; //to save memory we will only create one binning object

  SpecUtilsAsync::ThreadPool threadpool;

  for( auto &m : measurements_ )
  {
    if( m->gamma_counts_->empty() )
      continue;

    const bool noRecalNeeded = ( m->energy_calibration_model_ == type
                                 && ((type==Measurement::LowerChannelEdge && m->channel_energies_ && (*(m->channel_energies_)==eqn))
                                     || (type!=Measurement::LowerChannelEdge && m->calibration_coeffs_ == eqn) )
                                 && m->channel_energies_->size() > 0
                                 && m->deviation_pairs_ == dev_pairs );

    if( !new_binning )
    {
      if( !noRecalNeeded )
        m->rebin_by_eqn( eqn, dev_pairs, type );
      new_binning = m->channel_energies_;
    }//if( !new_binning )

    if( !noRecalNeeded )
    {
      threadpool.post( std::bind( &Measurement::rebin_by_lower_edge,
                                          m, std::cref(new_binning)  ) );
    }//if( !noRecalNeeded )
  }//for( auto &m : measurements_ )

  threadpool.join();

  for( auto &m : measurements_ )
  {
    m->calibration_coeffs_ = eqn;
    m->deviation_pairs_    = dev_pairs;
    m->energy_calibration_model_  = type;
  }//for( auto &m : measurements_ )
  
  properties_flags_ |= kHasCommonBinning;
  
  modified_ = modifiedSinceDecode_ = true;
}//void rebin_by_eqn( const std::vector<float> &eqn )


void MeasurementInfo::recalibrate_by_lower_edge( ShrdConstFVecPtr binning )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  for( auto &m : measurements_ )
    if( !(m->gamma_counts_->empty()) )
      m->recalibrate_by_eqn( vector<float>(), DeviationPairVec(), Measurement::LowerChannelEdge, binning );
  
  modified_ = modifiedSinceDecode_ = true;
}//void recalibrate_by_lower_edge( const std::vector<float> &binning )



void MeasurementInfo::recalibrate_by_eqn( const std::vector<float> &eqn,
                                          const DeviationPairVec &dev_pairs,
                                          Measurement::EquationType type )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );


  ShrdConstFVecPtr new_binning; //to save memory we will only create one binning object

  for( auto &m : measurements_ )
  {
    if( m->gamma_counts_->empty() )
      continue;

    if( !new_binning )
    {
      m->recalibrate_by_eqn( eqn, dev_pairs, type, ShrdConstFVecPtr() );
      new_binning = m->channel_energies_;
    }else
    {
      m->recalibrate_by_eqn( eqn, dev_pairs, type, new_binning );
    }
  }//for( auto &m : measurements_ )
  
  properties_flags_ |= kHasCommonBinning;
  
  modified_ = modifiedSinceDecode_ = true;
}//void recalibrate_by_eqn( const std::vector<float> &eqn )


//If only certain detectors are specified, then those detectors will be
//  recalibrated, and the other detectors will be rebinned.
void MeasurementInfo::recalibrate_by_eqn( const std::vector<float> &eqn,
                                          const DeviationPairVec &dev_pairs,
                                          Measurement::EquationType type,
                                          const vector<string> &detectors,
                                          const bool rebin_other_detectors )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( detectors.empty() )
    throw runtime_error( "recalibrate_by_eqn(...): an empty set of detectors"
                         " to recalibrate was passed in." );
  
  //If we are recalibrating all the detectors, lets just call the other version
  //  of this function, but it doenst really matter if we dont
  if( detectors == detector_names_ )
  {
    recalibrate_by_eqn( eqn, dev_pairs, type );
    return;
  }//if( detectors == detector_names_ )
  
  
  ShrdConstFVecPtr new_binning; //to save memory we will only create one binning object
  vector<MeasurementShrdPtr> rebinned;
  
  for( auto &m : measurements_ )
  {
    if( m->gamma_counts_->empty() )
      continue;
    
    vector<string>::const_iterator name_pos = std::find( detectors.begin(),
                                        detectors.end(), m->detector_name() );
    if( name_pos == detectors.end() )
    {
      rebinned.push_back( m );
    }else if( !new_binning )
    {
      m->recalibrate_by_eqn( eqn, dev_pairs, type, ShrdConstFVecPtr() );
      new_binning = m->channel_energies_;
    }else
    {
      m->recalibrate_by_eqn( eqn, dev_pairs, type, new_binning );
    }
  }//for( auto &m : measurements_ )
  
  if( !new_binning )
    throw runtime_error( "recalibrate_by_eqn(...): no valid detector"
                        " names were passed in to recalibrate." );
  
  if( rebin_other_detectors && rebinned.size() )
  {
    if( rebinned.size() > 4 )
    {
    SpecUtilsAsync::ThreadPool threadpool;
    for( auto &m : rebinned )
    {
      threadpool.post( std::bind( &Measurement::rebin_by_lower_edge, m,
                                    std::cref(new_binning) ) );
    }
    threadpool.join();
    }else
    {
      for( auto &m : rebinned )
        m->rebin_by_lower_edge( new_binning );
    }
  
    for( auto &m : rebinned )
    {
      m->calibration_coeffs_ = eqn;
      m->deviation_pairs_    = dev_pairs;
      m->energy_calibration_model_  = type;
    }//for( auto &m : rebinned )
    
    properties_flags_ |= kHasCommonBinning;
  }else
  {
    //check to see if the binning of all the measurments_ happens to be same
    bool allsame = true;
    for( size_t i = 0; allsame && i < rebinned.size(); ++i )
      allsame = ((rebinned[i]->energy_calibration_model_ == type)
                 && (rebinned[i]->calibration_coeffs_ == eqn)
                 && (rebinned[i]->deviation_pairs_ == dev_pairs));

    if( allsame )
      properties_flags_ |= kHasCommonBinning;
    else
      properties_flags_ &= (~kHasCommonBinning);
  }//if( rebin_other_detectors && rebinned.size() )

  modified_ = modifiedSinceDecode_ = true;
}//void recalibrate_by_eqn(...)




size_t MeasurementInfo::memmorysize() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  size_t size = sizeof(*this);

  size += filename_.capacity()*sizeof(string::value_type);
  for( const string &s : detector_names_ )
    size += s.capacity()*sizeof(string::value_type);
  size += detector_numbers_.capacity()*sizeof(int);
  for( const string &s : neutron_detector_names_ )
    size += s.capacity()*sizeof(string::value_type);

  size += uuid_.capacity()*sizeof(string::value_type);
  for( const string &str : remarks_ )
    size += str.capacity()*sizeof(string::value_type);
  size += measurement_location_name_.capacity()*sizeof(string::value_type);
  size += inspection_.capacity()*sizeof(string::value_type);
  size += sample_numbers_.size()*sizeof(int);

  size += sample_to_measurments_.size() * sizeof(vector<size_t>);
  typedef std::map<int, std::vector<size_t> > InIndVecMap;
  for( const InIndVecMap::value_type &t : sample_to_measurments_ )
    size += t.second.capacity() * sizeof(size_t);

  size += instrument_type_.capacity()*sizeof(string::value_type);
  size += manufacturer_.capacity()*sizeof(string::value_type);
  size += instrument_model_.capacity()*sizeof(string::value_type);
  size += instrument_id_.capacity()*sizeof(string::value_type);

  size += measurements_.capacity() * sizeof( MeasurementShrdPtr );

  map< const vector<float> *, int> binnings;
  for( const auto &m : measurements_ )
  {
    size += m->memmorysize();
    const vector<float> *bin_ptr = m->channel_energies_.get();
    if( binnings.find( bin_ptr ) == binnings.end() )
      binnings[ bin_ptr ] = 0;
    ++binnings[bin_ptr];
  }//for( const auto &m, measurements_ )

  //What about this->channel_energies_ ?  In principle we've already counted
  //  the size of this->channel_energies_
  for( const map< const vector<float> *, int>::value_type &entry : binnings )
    if( entry.first && entry.second > 1 )
      size -= ((entry.second-1)*( sizeof(*(entry.first)) + entry.first->capacity()*sizeof(float)));

  return size;
}//size_t memmorysize() const


bool MeasurementInfo::passthrough() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  return (properties_flags_ & kPassthroughOrSearchMode);
}//bool passthrough() const


size_t MeasurementInfo::suggested_gamma_binning_index(
                                            const std::set<int> &sample_numbers,
                                            const vector<bool> &det_to_use ) const
{
  //XXX - this function is a quickly implemented (temporary?) hack
  size_t index;
  ShrdConstFVecPtr binning_ptr;
    
  if( detector_numbers_.size() != det_to_use.size() )
    throw runtime_error( "MeasurementInfo::suggested_gamma_binning_index():"
                         " invalid det_to_use" );
  
  const bool same_nchannel = ((properties_flags_ & kAllSpectraSameNumberChannels) != 0);
  
  for( size_t i = 0; i < measurements_.size(); ++i )
  {
    const MeasurementShrdPtr &meas = measurements_[i];
    
    if( !sample_numbers.empty() && !sample_numbers.count(meas->sample_number_) )
      continue;
    
    const size_t det_index = std::find( detector_numbers_.begin(),
                                        detector_numbers_.end(),
                                        meas->detector_number_ )
                              - detector_numbers_.begin();
    
    const ShrdConstFVecPtr &thisbinning = meas->channel_energies();
    
    if( det_to_use.at(det_index) && !!thisbinning && !thisbinning->empty() )
    {
#if(!PERFORM_DEVELOPER_CHECKS)
      if( !binning_ptr || (binning_ptr->size() < thisbinning->size()) )
      {
        if( same_nchannel )
          return i;
        
        index = i;
        binning_ptr = thisbinning;
      }//if( !binning_ptr || (binning_ptr->size() < thisbinning->size()) )
#else
      if( !binning_ptr || (binning_ptr->size() < thisbinning->size()) )
      {
        if( !!binning_ptr
            && (binning_ptr->size() != thisbinning->size())
            && same_nchannel )
        {
          char buffer[512];
          snprintf( buffer, sizeof(buffer),
                    "Found instance of differening number of gamma channels,"
                    " when I shouldnt have; measurment %i had %i channels,"
                    " while measurment %i had %i channels.",
                   static_cast<int>(index),
                   static_cast<int>(binning_ptr->size()),
                   static_cast<int>(i),
                   static_cast<int>(thisbinning->size() ) );
          log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
        }
        
        index = i;
        binning_ptr = thisbinning;
      }//if( !binning_ptr || (binning_ptr->size() < thisbinning->size()) )
#endif
    }//if( use thisbinning )
  }//for( size_t i = 0; i < measurements_.size(); ++i )
  
  if( !binning_ptr )
    throw runtime_error( "MeasurementInfo::suggested_gamma_binning_index():"
                         " no valid measurments." );
  
  return index;
}//suggested_gamma_binning_index(...)


MeasurementShrdPtr MeasurementInfo::sum_measurements( const set<int> &sample_num,
                                            const vector<int> &det_nums ) const
{
  vector<bool> det_to_use( detector_numbers_.size(), false );
 
  for( const int num : det_nums )
  {
    vector<int>::const_iterator pos = std::find( detector_numbers_.begin(),
                                                 detector_numbers_.end(),
                                                 num );
    if( pos == detector_numbers_.end() )
      throw runtime_error( "MeasurementInfo::sum_measurements(): invalid detector number in the input" );
    
    const size_t index = pos - detector_numbers_.begin();
    det_to_use[index] = true;
  }//for( const int num : det_nums )
  
  return sum_measurements( sample_num, det_to_use );
}//sum_measurements(...)


MeasurementShrdPtr MeasurementInfo::sum_measurements( const set<int> &sample_num,
                                                     const vector<int> &det_nums,
                                                     const MeasurementConstShrdPtr binTo ) const
{
  vector<bool> det_to_use( detector_numbers_.size(), false );
  
  for( const int num : det_nums )
  {
    vector<int>::const_iterator pos = std::find( detector_numbers_.begin(),
                                                detector_numbers_.end(),
                                                num );
    if( pos == detector_numbers_.end() )
      throw runtime_error( "MeasurementInfo::sum_measurements(): invalid detector number in the input" );
    
    const size_t index = pos - detector_numbers_.begin();
    det_to_use[index] = true;
  }//for( const int num : det_nums )
  
  return sum_measurements( sample_num, det_to_use, binTo );
}//sum_measurements(...)



MeasurementShrdPtr MeasurementInfo::sum_measurements(
                                         const std::set<int> &sample_numbers,
                                         const vector<bool> &det_to_use ) const
{
  //For example example passthrough (using late 2011 mac book pro):
  //  Not using multithreading:          It took  0.135964s wall, 0.130000s user + 0.000000s system = 0.130000s CPU (95.6%)
  //  Using multithreading (4 threads):  It took  0.061113s wall, 0.220000s user + 0.000000s system = 0.220000s CPU (360.0%)
  //For example example HPGe Ba133:
  //  Not using multithreading:         It took  0.001237s wall, 0.000000s user
  //  Using multithreading (2 threads): It took  0.001248s wall, 0.000000s user

  size_t calIndex;
  try
  {
    calIndex = suggested_gamma_binning_index( sample_numbers, det_to_use );
  }catch(...)
  {
    return MeasurementShrdPtr();
  }
  
  return sum_measurements( sample_numbers, det_to_use, measurements_[calIndex] );
}//sum_measurements(...)


MeasurementShrdPtr MeasurementInfo::sum_measurements( const std::set<int> &sample_numbers,
                                      const vector<bool> &det_to_use,
                                      const MeasurementConstShrdPtr binto ) const
{
  if( !binto || measurements_.empty() )
    return MeasurementShrdPtr();
  
  MeasurementShrdPtr dataH = std::make_shared<Measurement>();
  
  dataH->energy_calibration_model_  = binto->energy_calibration_model_;
  dataH->calibration_coeffs_ = binto->calibration_coeffs_;
  dataH->deviation_pairs_    = binto->deviation_pairs_;
  dataH->channel_energies_   = binto->channel_energies();
  
  if( measurements_.size() == 1 )
    dataH->set_title( measurements_[0]->title_ );
  else
    dataH->set_title( filename_ );
  
  if( detector_names_.size() != det_to_use.size() )
    throw runtime_error( "MeasurementInfo::sum_measurements(...): "
                        "det_to_use.size() != sample_measurements.size()" );
  
  dataH->contained_neutron_ = false;
  dataH->real_time_ = 0.0f;
  dataH->live_time_ = 0.0f;
  dataH->gamma_count_sum_ = 0.0;
  dataH->neutron_counts_sum_ = 0.0;
  dataH->sample_number_ = -1;
  dataH->detector_name_ = "Summed";
  dataH->detector_number_ = -1;
  if( sample_numbers.size() == 1 )
    dataH->sample_number_ = *sample_numbers.begin();
  dataH->start_time_ = boost::posix_time::pos_infin;
  
  const size_t nenergies = dataH->channel_energies_ ? dataH->channel_energies_->size() : size_t(0);
  
  size_t ndet_to_use = 0;
  for( const bool i : det_to_use ) //could use std::count_if(...), for std::for_each(...) ...
    ndet_to_use += static_cast<size_t>( i );

  
  bool allBinningIsSame = ((properties_flags_ & kHasCommonBinning) != 0);
  
  if( allBinningIsSame )
  {
    const vector< std::shared_ptr<Measurement> >::const_iterator pos
               = std::find( measurements_.begin(), measurements_.end(), binto );
    const bool bintoInMeas = (pos != measurements_.end());
    if( !bintoInMeas )
      allBinningIsSame = (measurements_[0]->channel_energies() == binto->channel_energies());
  }//if( allBinningIsSame && measurements_.size() )

  //any less than 'min_per_thread' than the additional memorry allocation isnt
  //  worth it - briefly determined on my newer mbp using both example
  //  example passthrough (16k channel) and a 512 channel NaI portal passthrough
  const size_t min_per_thread = 8;

  size_t num_thread = static_cast<size_t>( SpecUtilsAsync::num_physical_cpu_cores() );
  const size_t num_potential_specta = ndet_to_use * sample_numbers.size();
  num_thread = min( num_thread, num_potential_specta/min_per_thread );
  num_thread = max( size_t(1), num_thread );
  
  vector< vector< std::shared_ptr<const Measurement> > > specs( num_thread );
  vector< vector< std::shared_ptr<const vector<float> > > > spectrums( num_thread );
  
  int current_total_sample_num = 0;
  set<string> remarks;
  for( const int sample_number : sample_numbers )
  {
    for( size_t index = 0; index < detector_names_.size(); ++index )
    {
      const std::string &det = detector_names_[index];
      MeasurementConstShrdPtr meas = measurement( sample_number, det );
      
      if( !meas )
        continue;
      
      std::shared_ptr<const vector<float> > spec = meas->gamma_counts();
      
      const size_t spec_size = (spec ? spec->size() : (size_t)0);
      
      //we'll allow measurement->channel_energies() to have more bins, since
      //  if meas->energy_calibration_model() == LowerChannelEdge then ther will be
      //  an extra bin to mark overflow
      if( allBinningIsSame && (spec_size > nenergies) )
      {
        stringstream msg;
        msg << SRC_LOCATION << "\n\tspec.size()=" << spec_size
            << "  measurement->channel_energies().size()=" << nenergies;
        cerr << msg.str() << endl;
        throw runtime_error( msg.str() );
      }//if( spec.size() > (binning_ptr->size()) )
      
      //Could add consistency check here to make sure all spectra are same size
      
      if( det_to_use[index] )
      {
        dataH->start_time_ = std::min( dataH->start_time_, meas->start_time_ );
        
        if( binto->gamma_counts_ && (binto->gamma_counts_->size() > 3) )
        {
          if( meas->gamma_counts_ && (meas->gamma_counts_->size() > 3) )
          {
            dataH->live_time_ += meas->live_time();
            dataH->real_time_ += meas->real_time();
          }
        }else
        {
          dataH->live_time_ += meas->live_time();
          dataH->real_time_ += meas->real_time();
        }
        
        dataH->neutron_counts_sum_ += meas->neutron_counts_sum();
        dataH->contained_neutron_ |= meas->contained_neutron_;
        
        if( dataH->neutron_counts_.size() < meas->neutron_counts_.size() )
          dataH->neutron_counts_.resize( meas->neutron_counts_.size(), 0.0f );
        const size_t nneutchannel = meas->neutron_counts_.size();
        for( size_t i = 0; i < nneutchannel; ++i )
          dataH->neutron_counts_[i] += meas->neutron_counts_[i];
        
        for( const std::string &remark : meas->remarks_ )
          remarks.insert( remark );
        
        if( spec_size )
        {
          dataH->gamma_count_sum_ += meas->gamma_count_sum();
          const size_t thread_num = current_total_sample_num % num_thread;
          specs[thread_num].push_back( meas );
          spectrums[thread_num].push_back( meas->gamma_counts() );
          ++current_total_sample_num;
        }
      }//if( det_to_use[index] )
    }//for( size_t index = 0; index < detector_names_.size(); ++index )
  }//for( const int sample_number : sample_numbers )
  
  if( !current_total_sample_num )
    return dataH->contained_neutron_ ? dataH : MeasurementShrdPtr();
  
  //If we are only summing one sample, we can preserve some additional
  //  information
  if( current_total_sample_num == 1 )
  {
    dataH->latitude_ = specs[0][0]->latitude_;
    dataH->longitude_ = specs[0][0]->longitude_;
    dataH->position_time_ = specs[0][0]->position_time_;
    dataH->sample_number_ = specs[0][0]->sample_number_;
    dataH->occupied_ = specs[0][0]->occupied_;
    dataH->speed_ = specs[0][0]->speed_;
    dataH->detector_name_ = specs[0][0]->detector_name_;
    dataH->detector_number_ = specs[0][0]->detector_number_;
    dataH->detector_type_ = specs[0][0]->detector_type_;
    dataH->quality_status_ = specs[0][0]->quality_status_;
  }//if( current_total_sample_num == 1 )
  
  if( allBinningIsSame )
  {
    if( num_thread > 1 )
    {
      //Should consider using calloc( )  and free...
      vector< vector<float> > results( num_thread );
    
      SpecUtilsAsync::ThreadPool threadpool;
      for( size_t i = 0; i < num_thread; ++i )
        threadpool.post( std::bind( &add_to, std::ref(results[i]), std::cref(spectrums[i]) ) );
      threadpool.join();
    
      //Note: in principle for a multicore machine (>16 physical cores), we
      //      could combine results using a few different threads down to less
      //      than 'min_per_thread'
      const size_t spec_size = results[0].size();
      std::shared_ptr<vector<float> > result_vec = std::make_shared< vector<float> >(spec_size, 0.0f);
      dataH->gamma_counts_ = result_vec;
    
      for( size_t i = 0; i < num_thread; ++i )
      {
        if( results[i].size() )
        {
          const float *spec_array = &(results[i][0]);
          for( size_t bin = 0; bin < spec_size; ++bin )
            result_vec->operator[](bin) += spec_array[bin];
        }
      }//for( size_t i = 0; i < num_thread; ++i )
    }else
    {
      if( spectrums.size()!=1 || spectrums[0].empty() )
        throw runtime_error( string(SRC_LOCATION) + "\n\tSerious programming logic error" );
    
      const vector< std::shared_ptr<const vector<float> > > &spectra = spectrums[0];
      const size_t num_spectra = spectra.size();
    
      const size_t spec_size = spectra[0]->size();
      vector<float> *result_vec = new vector<float>( spec_size, 0.0 );
      dataH->gamma_counts_.reset( result_vec );
      float *result_vec_raw = &((*result_vec)[0]);
    
      for( size_t i = 0; i < num_spectra; ++i )
      {
        const size_t len = spectra[i]->size();
        const float *spec_array = &(spectra[i]->operator[](0));
      
        //Using size_t to get around possible variable size issues.
        //  In principle I would expect this bellow loop to get vectorized by the
        //  compiler - but I havent actually checked for this.
        for( size_t bin = 0; bin < spec_size && bin < len; ++bin )
          result_vec_raw[bin] += spec_array[bin];
      }//for( size_t i = 0; i < num_thread; ++i )
    }//if( num_thread > 1 ) / else
  
  }else //if( allBinningIsSame )
  {
    if( sample_numbers.size() > 1 || det_to_use.size() > 1 )
      cerr << "sum_measurements for case with without a common binning not tested yet!" << endl;
    
    vector< vector<float> > results( num_thread );
    SpecUtilsAsync::ThreadPool threadpool;
    for( size_t i = 0; i < num_thread; ++i )
      threadpool.post( std::bind( &sum_with_rebin, std::ref(results[i]), std::cref(dataH), std::cref(specs[i]) ) );
    threadpool.join();
    
    const size_t spec_size = results[0].size();
    vector<float> *result_vec = new vector<float>( results[0] );
    dataH->gamma_counts_.reset( result_vec );
    float *result_vec_raw = &((*result_vec)[0]);
    
    for( size_t i = 1; i < num_thread; ++i )
    {
      const size_t len = results[i].size();
      const float *spec_array = &(results[i][0]);
      for( size_t bin = 0; bin < spec_size && bin < len; ++bin )
        result_vec_raw[bin] += spec_array[bin];
    }//for( size_t i = 0; i < num_thread; ++i )
  }//if( allBinningIsSame ) / else
  
  if( dataH->start_time_.is_infinity() )
    dataH->start_time_ = boost::posix_time::not_a_date_time;
  
  for( const std::string &remark : remarks )
    dataH->remarks_.push_back( remark );

  return dataH;
}//std::shared_ptr<Measurement> sum_measurements( int &, int &, const SpecMeas & )



set<size_t> MeasurementInfo::gamma_channel_counts() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  set<size_t> answer;

  for( const auto &meas : measurements_ )
    if( meas && meas->channel_energies_ && !meas->channel_energies_->empty() )
      answer.insert( meas->channel_energies_->size() );

  return answer;
}//std::set<size_t> gamma_channel_counts() const


size_t MeasurementInfo::num_gamma_channels() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  for( const auto &meas : measurements_ )
    if( meas && meas->channel_energies_ && !meas->channel_energies_->empty() )
      return meas->channel_energies_->size();

  return 0;
}//size_t MeasurementInfo::num_gamma_channels() const


struct KeepNBinSpectraStruct
{
  typedef std::vector< MeasurementShrdPtr >::const_iterator MeasVecIter;

  size_t m_nbin;
  MeasVecIter m_start, m_end;
  std::vector< MeasurementShrdPtr > *m_keepers;

  KeepNBinSpectraStruct( size_t nbin,
                         MeasVecIter start, MeasVecIter end,
                         std::vector< MeasurementShrdPtr > *keepers )
    : m_nbin( nbin ), m_start( start ), m_end( end ), m_keepers( keepers )
  {
  }

  void operator()()
  {
    m_keepers->reserve( m_end - m_start );
    for( MeasVecIter iter = m_start; iter != m_end; ++iter )
    {
      const MeasurementShrdPtr &m = *iter;
      const ShrdConstFVecPtr channel_energies = (m ? m->channel_energies()
                                                   : ShrdConstFVecPtr());

      if( !( ( !channel_energies || channel_energies->size()!=m_nbin )
          && (channel_energies || !channel_energies->empty()) ) )
      {
        m_keepers->push_back( m );
      }
    }//for( MeasVecIter iter = m_start; iter != m_end; ++iter )
  }//void operator()()
};//struct KeepNBinSpectraStruct


//return number of removed spectra
size_t MeasurementInfo::keep_n_bin_spectra_only( size_t nbin )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  size_t nremoved = 0;
  const size_t nstart = measurements_.size();

  std::vector< MeasurementShrdPtr > newmeass;
  newmeass.reserve( nstart );

  if( nstart < 100 )
  {
    KeepNBinSpectraStruct worker( nbin, measurements_.begin(),
                                  measurements_.end(), &newmeass );
    worker();
    nremoved = nstart - newmeass.size();
  }else
  {
    const int nthread = SpecUtilsAsync::num_logical_cpu_cores();
    const size_t meas_per_thread = nstart / nthread;

    SpecUtilsAsync::ThreadPool threadpool;
//    vector<KeepNBinSpectraStruct> workers;
    const KeepNBinSpectraStruct::MeasVecIter begin = measurements_.begin();

    size_t nsections = nstart / meas_per_thread;
    if( (nstart % meas_per_thread) != 0 )
      nsections += 1;

    vector< vector< MeasurementShrdPtr > > answers( nsections );

    size_t sec_num = 0;
    for( size_t pos = 0; pos < nstart; pos += meas_per_thread )
    {
      KeepNBinSpectraStruct::MeasVecIter start, end;
      start = begin + pos;
      if( (pos + meas_per_thread) <= nstart )
        end = start + meas_per_thread;
      else
        end = begin + nstart;

      threadpool.post( KeepNBinSpectraStruct(nbin, start, end, &(answers[sec_num]) ) );
      ++sec_num;
    }//for( size_t pos = 0; pos < nstart; pos += meas_per_thread )

    if( nsections != sec_num )
      throw runtime_error( SRC_LOCATION + "\n\tSerious logic error here!" );

    threadpool.join();

    for( size_t i = 0; i < nsections; ++i )
    {
      const vector< MeasurementShrdPtr > &thismeas = answers[i];
      newmeass.insert( newmeass.end(), thismeas.begin(), thismeas.end() );
    }//for( size_t i = 0; i < nsections; ++i )
    nremoved = nstart - newmeass.size();
  }//if( nstart < 100 )

  if( nremoved )
  {
    measurements_.swap( newmeass );
    cleanup_after_load();
  }//if( nremoved )

  return nremoved;
}//size_t keep_n_bin_spectra_only( size_t nbin )


size_t MeasurementInfo::remove_neutron_measurments()
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  size_t nremoved = 0;

  for( size_t i = 0; i < measurements_.size();  )
  {
    MeasurementShrdPtr &m = measurements_[i];
    if( m->contained_neutron_
        && (!m->channel_energies_ || m->channel_energies_->empty()) )
    {
      ++nremoved;
      measurements_.erase( measurements_.begin()+i );
    }else ++i;
  }//for( size_t i = 0; i < norig; ++i )

  if( nremoved )
  {
    cleanup_after_load();
  }//if( nremoved )

  if( nremoved )
    modified_ = modifiedSinceDecode_ = true;
  
  return nremoved;
}//size_t remove_neutron_measurments();


set<string> MeasurementInfo::energy_cal_variants() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  set<string> answer;
  
  for( const std::string &detnam : detector_names_ )
  {
    const size_t pos = detnam.find( "_intercal_" );
    if( pos != string::npos )
      answer.insert( detnam.substr(pos + 10) );
  }

  return answer;
}//set<string> energy_cal_variants() const


size_t MeasurementInfo::keep_energy_cal_variant( const std::string variant )
{
  const string ending = "_intercal_" + variant;
  std::vector< MeasurementShrdPtr > keepers;
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const set<string> origvaraints = energy_cal_variants();
  
  if( !origvaraints.count(variant) )
    throw runtime_error( "MeasurementInfo::keep_energy_cal_variant():"
                         " measurment did not contain an energy variant named '"
                         + variant + "'" );
  
  if( origvaraints.size() == 1 )
    return 0;
  
  keepers.reserve( measurements_.size() );
  
  for( auto &ptr : measurements_ )
  {
    const std::string &detname = ptr->detector_name_;
    const size_t pos = detname.find( "_intercal_" );
    if( pos == string::npos )
    {
      keepers.push_back( ptr );
    }else if( ((pos + ending.size()) == detname.size())
               && (strcmp(detname.c_str()+pos+10,variant.c_str())==0) )
    {
      ptr->detector_name_ = detname.substr( 0, pos );
      keepers.push_back( ptr );
    }
    //else
      //cout << "Getting rid of: " << detname << endl;
  }//for( auto &ptr : measurements_ )
  
  measurements_.swap( keepers );
  cleanup_after_load();
  
  modified_ = modifiedSinceDecode_ = true;
  
  return (keepers.size() - measurements_.size());
}//void keep_energy_cal_variant( const std::string variant )



int MeasurementInfo::background_sample_number() const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  //XXX - maybe could be sped up by using sample_numbers()
  //      and/or sample_measurements(..)
  for( const auto &meas : measurements_ )
    if( meas->source_type_ == Measurement::Background )
      return meas->sample_number_;

  return numeric_limits<int>::min();
}//int background_sample_number() const


void MeasurementInfo::reset()
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  gamma_live_time_ = 0.0;
  gamma_real_time_ = 0.0;
  gamma_count_sum_ = 0.0;
  neutron_counts_sum_ = 0.0;
  mean_latitude_ = mean_longitude_ = -999.9;
  properties_flags_ = 0x0;
  filename_ =  "";
  detector_names_.clear();
  neutron_detector_names_.clear();
  uuid_.clear();
  remarks_.clear();
  lane_number_ = -1;
  measurement_location_name_.clear();
  inspection_.clear();
  measurment_operator_.clear();
  sample_numbers_.clear();
  sample_to_measurments_.clear();
  detector_type_ = kUnknownDetector;
  instrument_type_.clear();
  manufacturer_.clear();
  instrument_model_.clear();
  instrument_id_.clear();
  measurements_.clear();
  detector_numbers_.clear();
  modified_ = modifiedSinceDecode_ = false;
}//void MeasurementInfo::reset()

  
DetectorAnalysisResult::DetectorAnalysisResult()
{
  reset();
}


void DetectorAnalysisResult::reset()
{
  remark_.clear();
  nuclide_.clear();
  activity_ = -1.0;            //in units of becquerel (eg: 1.0 == 1 bq)
  nuclide_type_.clear();
  id_confidence_.clear();
  distance_ = -1.0f;
  dose_rate_ = -1.0f;
  real_time_ = -1.0f;           //in units of secons (eg: 1.0 = 1 s)
  //start_time_ = boost::posix_time::ptime();
  detector_.clear();
}//void DetectorAnalysisResult::reset()


DetectorAnalysis::DetectorAnalysis()
{
  reset();
}

void DetectorAnalysis::reset()
{
  remarks_.clear();
  algorithm_name_.clear();
  algorithm_version_.clear();
  algorithm_creator_.clear();
  algorithm_description_.clear();
  algorithm_result_description_.clear();
  
  results_.clear();
}//void DetectorAnalysis::reset()



void MeasurementInfo::write_to_file( const std::string filename,
                                     const SaveSpectrumAsType format ) const
{
  set<int> samples, detectors;
  
  {
    std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
    samples = sample_numbers_;
    detectors = set<int>( detector_numbers_.begin(), detector_numbers_.end() );
  }
  
  write_to_file( filename, samples, detectors, format );
}//void write_to_file(...)



void MeasurementInfo::write_to_file( const std::string name,
                   const std::set<int> sample_nums,
                   const std::set<int> det_nums,
                   const SaveSpectrumAsType format ) const
{
  if( UtilityFunctions::is_file(name) || UtilityFunctions::is_directory(name) )
    throw runtime_error( "File (" + name + ") already exists, not overwriting" );
  
  std::ofstream output( name.c_str(), ios::out | ios::binary );

  if( !output )
    throw runtime_error( "Failed to open file (" + name + ") for writing" );
  
  write( output, sample_nums, det_nums, format );
}//void write_to_file(...)


void MeasurementInfo::write_to_file( const std::string name,
                   const std::vector<int> sample_nums_vector,
                   const std::vector<int> det_nums_vector,
                   const SaveSpectrumAsType format ) const
{
    //copy vectors into sets 
   const std::set<int> sample_nums_set( sample_nums_vector.begin(),
                                        sample_nums_vector.end() );
   const std::set<int> det_nums_set( det_nums_vector.begin(),
                                     det_nums_vector.end() );
    //write the file
    write_to_file( name, sample_nums_set, det_nums_set, format);
}//write_to_file(...)



void MeasurementInfo::write( std::ostream &strm,
           std::set<int> sample_nums,
           const std::set<int> det_nums,
           const SaveSpectrumAsType format ) const
{ 
  //Its a little heavy handed to lock mutex_ for this entire function, but
  //  techncically necassry since operator= operations are a shallow copy
  //  of the Measurements, meaning there could be an issue if the user modifies
  //  a Measurement in another thread.
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( sample_nums.empty() )
    throw runtime_error( "No sample numbers were specified to write out" );
  
  if( det_nums.empty() )
    throw runtime_error( "No detector numbers were specified to write out" );
  
  for( const int sample : sample_nums )
  {
    if( !sample_numbers_.count(sample ) )
      throw runtime_error( "Specified invalid sample number to write out" );
  }
  
  for( const int detnum : det_nums )
  {
    if( std::find(detector_numbers_.begin(), detector_numbers_.end(), detnum) == detector_numbers_.end() )
      throw runtime_error( "Specified invalid detector number to write out" );
  }
  
  MeasurementInfo info = *this;
  
  if( (sample_nums != sample_numbers_)
      || (det_nums.size() != detector_numbers_.size()) )
  {
    vector< std::shared_ptr<const Measurement> > toremove;
    for( std::shared_ptr<const Measurement> oldm : info.measurements() )
    {
      if( !sample_nums.count(oldm->sample_number())
         || !det_nums.count(oldm->detector_number()) )
      {
        toremove.push_back( oldm );
      }
    }//for( oldm : info.measurements() )
    
    info.remove_measurments( toremove );
  }//if( we dont want all the measurments )
  
  if( info.measurements_.empty() )
    throw runtime_error( "No Measurements to write out" );
  
  const std::set<int> &samples = info.sample_numbers_;
  const set<int> detectors( info.detector_numbers_.begin(), info.detector_numbers_.end() );
  
  bool success = false;
  switch( format )
  {
    case kTxtSpectrumFile:
      success = info.write_txt( strm );
      break;
      
    case kCsvSpectrumFile:
      success = info.write_csv( strm );
      break;
      
    case kPcfSpectrumFile:
      success = info.write_pcf( strm );
      break;
      
    case kXmlSpectrumFile:
      success = info.write_2006_N42( strm );
      break;
      
    case k2012N42SpectrumFile:
      success = info.write_2012_N42( strm );
      break;
      
    case kChnSpectrumFile:
      success = info.write_integer_chn( strm, samples, detectors );
      break;
      
    case kBinaryIntSpcSpectrumFile:
      success = info.write_binary_spc( strm, IntegerSpcType, samples, detectors );
      break;
      
    case kBinaryFloatSpcSpectrumFile:
      success = info.write_binary_spc( strm, FloatSpcType, samples, detectors );
      break;
      
    case kAsciiSpcSpectrumFile:
      success = info.write_ascii_spc( strm, samples, detectors );
      break;
      
    case kExploraniumGr130v0SpectrumFile:
      success = info.write_binary_exploranium_gr130v0( strm );
      break;
            
    case kExploraniumGr135v2SpectrumFile:
      success = info.write_binary_exploranium_gr135v2( strm );
      break;
      
    case kIaeaSpeSpectrumFile:
      success = info.write_iaea_spe( strm, samples, detectors );
      break;

#if( ENABLE_D3_CHART_EXPORTING )
    case kD3HtmlSpectrumFile:
    {
      D3SpectrumExport::D3SpectrumChartOptions options;
      success = info.write_d3_html( strm, options, samples, detectors );
      break;
    }
#endif
      
    case kNumSaveSpectrumAsType:
      throw runtime_error( "Invalid output format specified" );
      break;
  }//switch( format )
  
  if( !success )
    throw runtime_error( "Failed to write to output" );
  
}//write_to_file(...)



bool MeasurementInfo::write_pcf( std::ostream &outputstrm ) const
{
#if(PERFORM_DEVELOPER_CHECKS)
  //The input stream may not support tellp(), so for testing to make sure we
  //  have the positions correct, we'll use a stringstream, and then swap
  //  later on.
  stringstream ostr;
#else
  std::ostream &ostr = outputstrm;
#endif
  
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  try
  {
    size_t nchannel_file = num_gamma_channels();
    
    if( !(properties_flags_ & kAllSpectraSameNumberChannels) )
    {
      const set<size_t> n_counts = gamma_channel_counts();
      if( n_counts.size() )
        nchannel_file = *n_counts.rbegin();
    }//if( !(properties_flags_ & kAllSpectraSameNumberChannels) )
    
    if( !nchannel_file )
      throw runtime_error( "No measurments to write to PCF." );
    
    //We need to round nchannel_file up to the nearest 64 channels since each
    // record must be multiple of 256 bytes
    if( (nchannel_file % 64) != 0 )
      nchannel_file += (64 - (nchannel_file % 64));
    
    //We want to put the detector name in the "Title" of the PCF, but only if
    //  there is more than one detector.
    //  XXX - this is an expensive operation!
    set<string> gamma_det_names;
    
    for( const auto &m : measurements_ )
    {
      if( m && m->gamma_counts_ && m->gamma_counts_->size() )
        gamma_det_names.insert( m->detector_name_ );
    }
    
    const size_t num_gamma_detectors = gamma_det_names.size();
    
    
    std::basic_string<char> fileid;
    int16_t NRPS = 1 + static_cast<int16_t>( 4.0*nchannel_file/256.0 );
    fileid.resize( 2, '\0' );
    memcpy( &(fileid[0]), &NRPS, 2 );
    fileid += "DHS       " + uuid();
    fileid.resize( 48 , ' ' );
    fileid += inspection_;
    fileid.resize( 64, ' ' );
    int16_t lanenum = static_cast<int16_t>( lane_number_ );
    fileid.resize( 66, 0 );
    memcpy( &(fileid[64]), &lanenum, 2 );
    
    
    
    
    for( size_t i = 0; i < remarks_.size(); ++i )
    {
      string val = remarks_[i];
      UtilityFunctions::trim( val );
      if( val.empty()
         || UtilityFunctions::istarts_with(val, "ItemDescription")
         || UtilityFunctions::istarts_with(val, "CargoType")
         || UtilityFunctions::istarts_with(val, "ItemToDetectorDistance")
         || UtilityFunctions::istarts_with(val, "OccupancyNumber") )
        continue;
      fileid += (i ? "\r\n" : "") + val;
    }
    
    fileid.resize( 92, ' ' );
    
    fileid += instrument_type();
    fileid.resize( 120, ' ' );
    fileid += manufacturer();
    fileid.resize( 148, ' ' );
    fileid += instrument_model();
    fileid.resize( 166, ' ' );
    fileid += instrument_id();
    fileid.resize( 184, ' ' );
    
    string item_description = findPcfRemark("ItemDescription",remarks_);
    if( item_description.size() > 20 )
      item_description = item_description.substr(0,20);
    fileid += item_description;
    fileid.resize( 204, ' ' );
    
    fileid += measurement_location_name_;
    fileid.resize( 220, ' ' );
    
    if( has_gps_info() )
    {
      //we only have 16 bytes, here; we'll try printing to 7 decimal points, and
      //  if too long, print 5, and then 4, etc..
      size_t len = 63;
      char valbuffer[64] = { 0 };
      for( int ndecimals = 7; len > 16 && ndecimals > 2; --ndecimals )
      {
        char frmtstr[32];
        snprintf( frmtstr, sizeof(frmtstr), "%%.%if,%%.%if", ndecimals, ndecimals);
        snprintf( valbuffer, sizeof(valbuffer), frmtstr, mean_latitude(), mean_longitude() );
        len = strlen( valbuffer );
      }
      fileid += valbuffer;
    }//if( has_gps_info() )
    
    fileid.resize( 236, ' ' );
    
    fileid.resize( 238, '\0' ); //2-byte signed integer of Item to detector distance
    string item_dist_str = findPcfRemark("ItemToDetectorDistance",remarks_);
    int16_t itemdistance = static_cast<int16_t>( atoi( item_dist_str.c_str() ) );
    memcpy( &(fileid[236]), &itemdistance, 2 );
    
    fileid.resize( 240, '\0' ); //2-byte signed integer of Occupancy number
    int occnum = occupancy_number_from_remarks();
    if( occnum >= 0 )
    {
      const int16_t occ = static_cast<int16_t>( occnum );
      memcpy( &(fileid[238]), &occ, 2 );
    }
    
    string cargo_type = findPcfRemark("CargoType",remarks_);
    if( cargo_type.size() > 16 )
      cargo_type = cargo_type.substr(0,16);
    fileid += cargo_type;
    
    fileid.resize( 256, ' ' );
    ostr.write( &fileid[0], fileid.size() );
    
    //Find the deviation pairs to use in this file; we'll use the first onese
    //  we find.
    vector< pair<float,float> > dev_pairs;
    for( const auto &meas : measurements_ )
    {
      if( !!meas->gamma_counts_ && !meas->gamma_counts_->empty()
         && !!meas->channel_energies() && !meas->channel_energies()->empty() )
      {
        dev_pairs = meas->deviation_pairs_;
        break;
      }
    }//for( const auto &meas, measurements_ )
    
    std::basic_string<char> header; //has some structure ( eg says 'DeviationPairsInFile')
    header = "DeviationPairsInFile";
    header.resize( 256, ' ' );
    
    ostr.write( &header[0], header.size() );
    const size_t npairs = 80*256/sizeof(float);
    float dev_pair_data[npairs] = { 0.0f };
    for( size_t i = 0; i < dev_pairs.size() && i < npairs/2; ++i )
    {
      dev_pair_data[2*i] = dev_pairs[i].first;
      dev_pair_data[2*i+1] = dev_pairs[i].second;
    }//for( size_t i = 0; i < dev_pairs.size(); ++i )
    
    if( dev_pairs.size() )
      cerr << "MeasurementInfo::write_pcf: currently am not writing deviation pairs to PCF file coorectly" << endl;
    
    
    ostr.write( (char *)dev_pair_data, sizeof(dev_pair_data) );
    
    //Backgrounds dont count toward sample numbers for GADRAS; it also assumes
    //  samples start at 1 (like FORTRAN...).  So we will hack things a bit
    //  for passthrough()s; there is a little bit of checking to make sure
    //  sample numbers are kept in the same order as original, btu its not super
    //  robust or tested.
    vector<int> passthrough_samples;  //will be sorted with at most one of each value
    
    for( size_t i = 0; i < measurements_.size(); ++i )
    {
#if(PERFORM_DEVELOPER_CHECKS)
      istream::pos_type file_pos = ostr.tellp();
      
      if( (file_pos % 256) != 0 )
      {
        char buffer[256];
        snprintf( buffer, sizeof(buffer),
                 "When writing PCF file, at file position %i at start of spectrum %i when should be at a multiple of 256", int(file_pos), int(i) );
        log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
      }
#endif
      
      MeasurementConstShrdPtr meas = measurements_[i];
      
      if( !meas || !meas->gamma_counts_ || meas->gamma_counts_->empty() )
        continue;
      
      string spectrum_title;  //ex: 'Survey 1 Det=Aa1 Background @250cm'
      string additional_info; //ex '<T>HPGe 50%</T>' or '<T>Gamma</T>'
      string source_list;
      string collection_time;  //Formatted like: '2010-02-24T00:08:24.82Z'
      
      char character_tag;
      float live_time, true_time, halflife = 0.0, molecular_weight = 0.0,
      spectrum_multiplier = 0.0, offset, gain, quadratic, cubic, low_energy,
      neutron_counts;
      const int32_t num_channel = static_cast<int32_t>( meas->gamma_counts_->size() );
      
      live_time = meas->live_time_;
      true_time = meas->real_time_;
      
      char buffer[128];
      
      int sample_num = meas->sample_number_;
      if( passthrough() && (meas->source_type() != Measurement::Background) )
      {
        auto pos = std::lower_bound( begin(passthrough_samples), end(passthrough_samples), sample_num );
        sample_num = (pos - passthrough_samples.begin());
        
        //Will almost always be an insertion at the end - so its not the worst possible...
        if( (pos != end(passthrough_samples)) && ((*pos) != sample_num) )
          passthrough_samples.insert( pos, sample_num );
      }
      
      if( passthrough()
         && (meas->sample_number_ >= 0)
         && !UtilityFunctions::icontains( meas->title_, "sample" )
         && !UtilityFunctions::icontains( meas->title_, "survey" ) )
      {
        if( meas->source_type() == Measurement::Background )
          snprintf( buffer, sizeof(buffer), " Background" );
        else
          snprintf( buffer, sizeof(buffer), " Survey %i", sample_num );
        spectrum_title += buffer;
      }
      
      if( num_gamma_detectors > 1 )
        spectrum_title += (spectrum_title.empty() ? "Det=" : ": Det=") + meas->detector_name_;
      
      if( !passthrough()
         && !UtilityFunctions::icontains( meas->title_, "Background" )
         && !UtilityFunctions::icontains( meas->title_, "Foreground" ) )
      {
        if( meas->source_type_ == Measurement::Background )
          spectrum_title += " Background";
        else
          spectrum_title += " Foreground";
      }//if( not already labeled foreground or background )
      
      if( meas->speed_ > 0.0
         && !UtilityFunctions::icontains( meas->title_, "speed" ) )
      {
        snprintf( buffer, sizeof(buffer), " Speed %f m/s", meas->speed_ );
        spectrum_title += buffer;
      }
      
      
      //Lets keep the title from repeating background/foreground
      string old_title = meas->title_;
      const char * const check_paterns[] = { "background", "foreground" };
      for( const char *p : check_paterns )
        if( UtilityFunctions::icontains(old_title, p) && UtilityFunctions::icontains(spectrum_title, p) )
          UtilityFunctions::ireplace_all(old_title, p, "" );
      
      spectrum_title += " " + old_title;
      
      trim( spectrum_title );
      UtilityFunctions::ireplace_all( spectrum_title, "  ", " " );
      
      if( !meas->start_time_.is_special() )
      {
        stringstream stimestrm;
        stimestrm << meas->start_time_;
        collection_time = stimestrm.str();
      }//if( !meas->start_time_.is_special() )
      
      character_tag = '\0';//meas->cambio_tag_char();
      
      //From phone conversation with Dean 20170816:
      //  The meaning of the 'tag' character is highly overloaded, and can mean,
      //  among other usses:
      //    '-' not occupied, and anything else occupied - for RPM data
      //    '-' use a dashed line when plotting
      //    '<' Use filled region style when plotting
      //    'T' Calibration from thorium
      //    'K' Calibration from potasium
      
      if( passthrough() )
      {
        if( meas->occupied() == Measurement::NotOccupied )
          character_tag = '-';
        else if( meas->occupied() == Measurement::Occupied )
          character_tag = ' ';
      }
      
      vector<float> calib_coef = meas->calibration_coeffs_;
      if( meas->energy_calibration_model_ == Measurement::Polynomial )
        calib_coef = polynomial_coef_to_fullrangefraction( calib_coef, meas->gamma_counts_->size() );
      
      offset         = (calib_coef.size() > 0) ? calib_coef[0] : 0.0f;
      gain           = (calib_coef.size() > 1) ? calib_coef[1] : 0.0f;
      quadratic      = (calib_coef.size() > 2) ? calib_coef[2] : 0.0f;
      cubic          = (calib_coef.size() > 3) ? calib_coef[3] : 0.0f;
      low_energy     = 0.0f;
      if( meas->energy_calibration_model_ == Measurement::FullRangeFraction )
        low_energy   = (calib_coef.size() > 4) ? calib_coef[4] : 0.0f;
      
      
      const float occupied = (meas->occupied()==Measurement::Occupied ? 1.0f : 0.0f);
      
      //Isnt there more terms in newer PCF files?  Or are those only in the
      //  Detector.dat files
      
      neutron_counts = static_cast<float>( meas->neutron_counts_sum_ );
      
      spectrum_title.resize( 60, ' ' );
      additional_info.resize( 60, ' ' );
      source_list.resize( 60, ' ' );
      collection_time.resize( 23, ' ' );
      ostr.write( &(spectrum_title[0]), spectrum_title.size() );
      ostr.write( &(additional_info[0]), additional_info.size() );
      ostr.write( &(source_list[0]), source_list.size() );
      ostr.write( &(collection_time[0]), collection_time.size() );
      ostr.write( &character_tag, 1 );
      ostr.write( (char *)&live_time, 4 );
      ostr.write( (char *)&true_time, 4 );
      ostr.write( (char *)&halflife, 4 );
      ostr.write( (char *)&molecular_weight, 4 );
      ostr.write( (char *)&spectrum_multiplier, 4 );
      ostr.write( (char *)&offset, 4 );
      ostr.write( (char *)&gain, 4 );
      ostr.write( (char *)&quadratic, 4 );
      ostr.write( (char *)&cubic, 4 );
      ostr.write( (char *)&low_energy, 4 );
      ostr.write( (char *)&occupied, 4 ); //4-byte floating point, Occupancy flag (0 = unoccupied, 1 = occupied)
      ostr.write( (char *)&neutron_counts, 4 ); //
      ostr.write( (char *)&num_channel, 4 );
      
      if( meas->deviation_pairs_ == dev_pairs )
      {
        ostr.write( (char *)&(meas->gamma_counts_->operator[](0)), 4*num_channel );
      }else
      {
        //PCF files apply the deviation pairs to the whole file
        vector<float> new_counts;
        ShrdConstFVecPtr origbinning = meas->channel_energies_;
        if( !origbinning )
          origbinning = fullrangefraction_binning( calib_coef, num_channel,
                                                  meas->deviation_pairs_ );
        ShrdConstFVecPtr newbinning = fullrangefraction_binning( calib_coef,
                                                                num_channel, dev_pairs );
        rebin_by_lower_edge( *origbinning, *(meas->gamma_counts_),
                            *newbinning, new_counts );
        ostr.write( (char *)&(new_counts[0]), 4*num_channel );
      }
      
      //Incase this spectrum has less channels than 'nchannel_file'
      if( nchannel_file != num_channel )
      {
        //        assert( nchannel_file > num_channel );
        char dummies[4] = {'\0','\0','\0','\0'};
        for( size_t i = num_channel; i < nchannel_file; ++i )
          ostr.write( dummies, 4 );
      }//if( nchannel_file != num_channel )
      
    }//for( auto meas, measurements_ )
    
#if(PERFORM_DEVELOPER_CHECKS)
    if( !ostr.bad() )
    {
      outputstrm << ostr.rdbuf();
      return true;
    }
    return false;
#else
    return !ostr.bad();
#endif
  }catch( std::exception &e )
  {
    cerr << SRC_LOCATION << "\n\tCaught " << e.what() << endl;
  }
  
  return false;
}//bool write_pcf( std::ostream& ostr ) const


bool Measurement::write_2006_N42_xml( std::ostream& ostr ) const
{
  const char *endline = "\r\n";

  //ostr << "  <Measurement>" << endline;

  string detname = detector_name_;
  if( detname == "" )
    detname = s_unnamed_det_placeholder;

  
  if( contained_neutron_ )
  {
    ostr << "    <CountDoseData DetectorType=\"Neutron\">" << endline
         << "      <Counts>" << neutron_counts_sum_ << "</Counts>" << endline
         << "    </CountDoseData>" << endline;
  }//if( contained_neutron_ )

  //XXX - deviation pairs completeley untested, and probably not conformant to
  //      N42 2006 specification (or rather the dndons extention).
  //XXX - also, could put more detector information here!
/*
  if( deviation_pairs_.size() )
  {
    ostr << "    <InstrumentInformation>" << endline;
    ostr << "      <dndons:NonlinearityCorrection";
    if( detector_name_.size() )
      ostr << " Detector=\"" << detector_name_ << "\"";
    ostr << ">";
    
    for( size_t i = 0; i < deviation_pairs_.size(); ++i )
    {
      ostr << "        <dndons:Deviation>"
           << deviation_pairs_[i].first << " " << deviation_pairs_[i].second
           << "</dndons:Deviation>" << endline;
    }//for( size_t i = 0; i < deviation_pairs_.size(); ++i )
    
    ostr << "      </dndons:NonlinearityCorrection>";
    ostr << "    </InstrumentInformation>" << endline;
  }//if( deviation_pairs_.size() )
*/
  
//Further information we should write
//   std::string detector_name_;
//   int detector_number_;
//   QualityStatus quality_status_;
//   std::string title_;

/*
  if( Measurement::valid_latitude(latitude_)
      && Measurement::valid_longitude(longitude_) )
  {
    ostr << "    <MeasuredItemInformation>" << endline
         << "      <MeasurementLocation>" << endline;

//The MeasurementInfo class contains the member variable
//      measurement_location_name_, so we cant write the following
//  if( measurement_location_name_.size() )
//    ostr << "        <MeasurementLocationName>" << measurement_location_name_
//         << "</MeasurementLocationName>" << endline;
    ostr << "        <Coordinates";
    if( position_time_.is_special() )
      ostr << ">";
    else
      ostr << " Time=\"" << UtilityFunctions::to_extended_iso_string(position_time_) << "Z\">";
    
    ostr << latitude_ << " " << longitude_ << "</Coordinates>" << endline;
    ostr << "      </MeasurementLocation>" << endline;
    ostr << "    </MeasuredItemInformation>" << endline;
  }//if( valid(latitude_) && !valid(longitude) )
*/
  
  
  ostr << "    <Spectrum Type=\"PHA\"";
  
  ostr << " Detector=\"" << detname << "\"";
  if( sample_number_ > 0 )
    ostr << " SampleNumber=\"" << sample_number_ << "\"";
  
  switch( quality_status_ )
  {
    case Good:    ostr << " Quality=\"Good\""; break;
    case Suspect: ostr << " Quality=\"Suspect\""; break;
    case Bad:     ostr << " Quality=\"Bad\""; break;
    case Missing:
//      ostr << " Quality=\"Missing\"";
    break;
  }//switch( quality_status_ )
  
  ostr << ">" << endline;
  
  vector<string> remarks;
  
  if( title_.size() )
    remarks.push_back( "Title: " + title_ );
  
  bool wroteSurvey = false, wroteName = false, wroteSpeed = false;
  
  for( size_t i = 0; i < remarks_.size(); ++i )
  {
    const string &remark = remarks_[i];
    remarks.push_back( remark );
    
    if( i == 0 )
    {
      wroteSurvey = (remark.find( "Survey" ) != string::npos);
      wroteName   = (remark.find( detector_name_ ) != string::npos);
      wroteSpeed  = (remark.find( "Speed" ) != string::npos);
    }//if( i == 0 )
  }//for( size_t i = 0; i < remarks_.size(); ++i )

  if( remarks_.empty()
      && (sample_number_>=0 || !detector_name_.empty() || speed_>0.00000001) )
  {
    string thisremark;
    if( sample_number_>=0 && !wroteSurvey )
    {
      thisremark = "Survey " + std::to_string(sample_number_);
    }
    
    if( !detector_name_.empty() && !wroteName )
    {
      if(!thisremark.empty())
        thisremark += " ";
      thisremark += detector_name_;
    }
    
    if( (speed_ > 0.00000001) && !wroteSpeed )
    {
      if(!thisremark.empty())
       thisremark += " ";
      thisremark += "Speed " + std::to_string(speed_) + " m/s";
    }
    trim( thisremark );
    
    if(!thisremark.empty())
      remarks.push_back( thisremark );
  }//if( remarks_.empty() )
  
  //We're only allowed to have one Remark element here.
  if(!remarks.empty())
  {
    ostr << "      <Remark>";
    for( size_t i = 0; i < remarks.size(); ++i )
    {
      if( i )
        ostr << endline;
      ostr << remarks[i];
    }
    
    ostr << "</Remark>";
  }//if( remarks.size() )
  
  /*
  switch( occupied_ )
  {
    case UnknownOccupancyStatus:
    break;
    case Occupied:
    case NotOccupied:
      ostr << "      <Occupied>" << (occupied_==Occupied)
           << "</Occupied>" << endline;
    break;
  }//switch( occupied_ )
  */
  
  /*
  if( !start_time_.is_special() )
    ostr << "      <StartTime>"
         << UtilityFunctions::to_extended_iso_string(start_time_)
         << "Z</StartTime>" << endline;
//    ostr << "      <StartTime>" << start_time_ << "</StartTime>" << endline;
  */
  
  ostr << "      <RealTime>PT" << real_time_ << "S</RealTime>" << endline;
  ostr << "      <LiveTime>PT" << live_time_ << "S</LiveTime>" << endline;
  
  switch( source_type_ )
  {
    case IntrinsicActivity: ostr << "      <SourceType>Other</SourceType>" << endline; break; break;
    case Calibration:       ostr << "      <SourceType>Calibration</SourceType>" << endline; break;
    case Background:        ostr << "      <SourceType>Background</SourceType>" << endline; break;
    case Foreground:        ostr << "      <SourceType>Item</SourceType>" << endline; break;
    case UnknownSourceType: break;
  }//switch( source_type_ )
  
  if(!detector_type_.empty())
    ostr << "      <DetectorType>" << detector_type_ << "</DetectorType>" << endline;
  
  ostr << "      <Calibration Type=\"Energy\" EnergyUnits=\"keV\">" << endline
       << "        <Equation Model=\"";

  switch( energy_calibration_model_ )
  {
    case Polynomial:           ostr << "Polynomial";        break;
    case FullRangeFraction:    ostr << "FullRangeFraction"; break;
    case LowerChannelEdge:     ostr << "LowerChannelEdge";  break;
    case UnknownEquationType:  ostr << "Unknown";           break;
  }//switch( energy_calibration_model_ )

  ostr << "\">" << endline;

  ostr << "          <Coefficients>";
  for( size_t i = 0; i < calibration_coeffs_.size(); ++i )
    ostr << (i ? " " : "") << calibration_coeffs_[i];
 
  //Not certain we need this next loop, but JIC
  if( energy_calibration_model_ == LowerChannelEdge && calibration_coeffs_.empty()
      && !!channel_energies_ && !channel_energies_->empty())
  {
    for( size_t i = 0; i < channel_energies_->size(); ++i )
      ostr << (i ? " " : "") << (*channel_energies_)[i];
  }
  
  ostr << "</Coefficients>" << endline
       << "        </Equation>" << endline
       << "      </Calibration>" << endline;

  ostr << "      <ChannelData Compression=\"CountedZeroes\">";

  vector<float> compressed_counts;
  compress_to_counted_zeros( *gamma_counts_, compressed_counts );

  const size_t nCompressChannels = compressed_counts.size();
  for( size_t i = 0; i < nCompressChannels; ++i )
  {
    const size_t pos = (i%12);
    if( pos == 0 )
      ostr << endline;
    else
      ostr << " ";
    if( compressed_counts[i] == 0.0 )
      ostr << "0";
    else ostr << compressed_counts[i];
  }//for( size_t i = 0; i < compressed_counts.size(); ++i )

  ostr << "      </ChannelData>" << endline
//       << "      <DHS:Comment>Spectrum Saved From SrbSpectrumViewer</DHS:Comment>" << endline
//       << "      <DHS:InterSpecCompileDate>" << COMPILE_DATE_AS_INT << "</DHS:InterSpecCompileDate>" << endline
       << "    </Spectrum>" << endline
//       << "</Measurement>" << endline;
      ;
  
  return true;
}//bool write_2006_N42_xml( std::ostream& ostr ) const


bool MeasurementInfo::write_2006_N42( std::ostream& ostr ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  const char *endline = "\r\n";

  ostr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endline
       << "<N42InstrumentData xmlns=\"http://physics.nist.gov/Divisions/Div846/Gp4/ANSIN4242/2005/ANSIN4242\"" << endline
       << "xmlns:n42ns=\"http://physics.nist.gov/Divisions/Div846/Gp4/ANSIN4242/2005/ANSIN4242\"" << endline
       << "xmlns:dndons=\"http://www.DNDO.gov/N42Schema/2006/DNDOSchema\"" << endline
       << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endline
       << "xmlns:Cambio=\"Cambio\"" << endline
       << "xmlns:DHS=\"DHS\"" << endline
       << "xsi:schemaLocation=\"http://physics.nist.gov/Divisions/Div846/Gp4/ANSIN4242/2005/ANSIN4242" << endline
       << "http://physics.nist.gov/Divisions/Div846/Gp4/ANSIN4242/2005/ANSIN4242.xsd\">" << endline;
  
  ostr << "<Measurement UUID=\"" << uuid_ << "\">" << endline;

  ostr << "  <InstrumentInformation>" << endline;
  
  if(!instrument_type_.empty())
  {
    if( instrument_type_ == "PortalMonitor" || instrument_type_ == "SpecPortal"
        || instrument_type_ == "RadionuclideIdentifier" || instrument_type_ == "PersonalRadiationDetector"
        || instrument_type_ == "SurveyMeter" || instrument_type_ == "Spectrometer"
        || instrument_type_ == "Other" )
    {
      ostr << "    <InstrumentType>" << instrument_type_ << "</InstrumentType>" << endline;
    }else if( instrument_type_ == "Portal Monitor" )
      ostr << "    <InstrumentType>PortalMonitor</InstrumentType>" << endline;
    else if( instrument_type_ == "Radionuclide Identifier" )
      ostr << "    <InstrumentType>RadionuclideIdentifier</InstrumentType>" << endline;
    else if( instrument_type_ == "Spectroscopic Portal Monitor" )
      ostr << "    <InstrumentType>SpecPortal</InstrumentType>" << endline;
    else if( instrument_type_ == "Personal Radiation Detector" )
      ostr << "    <InstrumentType>PersonalRadiationDetector</InstrumentType>" << endline;
    else if( instrument_type_ == "Spectroscopic Personal Radiation Detector" )
      ostr << "    <InstrumentType>PersonalRadiationDetector</InstrumentType>" << endline;
    else if( instrument_type_ == "Transportable System" )
      ostr << "    <InstrumentType>Other</InstrumentType>" << endline;
    else if( instrument_type_ == "Gamma Handheld" )
      ostr << "    <InstrumentType>Spectrometer</InstrumentType>" << endline;
    else
      ostr << "<!-- <InstrumentType>" << instrument_type_ << "</InstrumentType> -->" << endline;
  }//if( instrument_type_.size() )
  
  if(!manufacturer_.empty())
    ostr << "    <Manufacturer>" << manufacturer_ << "</Manufacturer>" << endline;
  if(!instrument_model_.empty())
    ostr << "    <InstrumentModel>" << instrument_model_ << "</InstrumentModel>" << endline;
  if(!instrument_id_.empty())
    ostr << "    <InstrumentID>" << instrument_id_ << "</InstrumentID>" << endline;

  for( const string &detname : detector_names_ )
  {
    ostr << "    <dndons:DetectorStatus Detector=\""
         << (detname.empty() ? s_unnamed_det_placeholder : detname)
         << "\" Operational=\"true\"/>" << endline;
  }
  
  vector<string> unwrittendets = detector_names_;
  
  for( size_t i = 0; !unwrittendets.empty() && i < measurements_.size(); ++i )
  {
    const MeasurementShrdPtr &meas = measurements_[i];
    if( !meas )
      continue;
    
    vector<string>::iterator pos = std::find( unwrittendets.begin(), unwrittendets.end(), meas->detector_name_ );
    if( pos == unwrittendets.end() )
      continue;
    unwrittendets.erase( pos );

    if( !meas->gamma_counts() || meas->gamma_counts()->empty() )
      continue;

    string name = meas->detector_name_;
    if( name == "" )
      name = s_unnamed_det_placeholder;
    
    if( meas->deviation_pairs_.size() )
    {
      ostr << "    <dndons:NonlinearityCorrection Detector=\"" << name << "\">" << endline;
      for( size_t j = 0; j < meas->deviation_pairs_.size(); ++j )
      {
        ostr << "      <dndons:Deviation>" << meas->deviation_pairs_[j].first
             << " " << meas->deviation_pairs_[j].second << "</dndons:Deviation>" << endline;
      }
      ostr << "    </dndons:NonlinearityCorrection>" << endline;
    }
  }//for( size_t i = 0; i < measurements.size(); ++i )
  
  
  
  if( measurement_location_name_.size() || measurment_operator_.size() )
  {
    ostr << "    <MeasuredItemInformation>" << endline;
    if( measurement_location_name_.size() )
      ostr << "      <MeasurementLocationName>" << measurement_location_name_ << "</MeasurementLocationName>" << endline;
    //<Coordinates>40.12 10.67</Coordinates>
    if( measurment_operator_.size() )
      ostr << "      <MeasurementOperator>" << measurment_operator_ << "</MeasurementOperator>" << endline;
    ostr << "    </MeasuredItemInformation>" << endline;
  }
  
  ostr << "  </InstrumentInformation>" << endline;

  if( inspection_.size() )
    ostr << "  <dndons:Inspection>" << inspection_ << "</dndons:Inspection>" << endline;

  
  for( const int samplenum : sample_numbers_ )
  {
    vector<MeasurementConstShrdPtr> meass = sample_measurements( samplenum );
    if( meass.empty() )
      continue;
    
    //Look for min start time
    boost::posix_time::ptime starttime = meass[0]->start_time();
    float rtime = meass[0]->real_time_;
    float speed = meass[0]->speed_;
    Measurement::OccupancyStatus occstatus = meass[0]->occupied_;
    
    for( size_t i = 1; i < meass.size(); ++i )
    {
      const boost::posix_time::ptime tst = meass[i]->start_time();
      starttime = ((tst.is_special() || (starttime < tst)) ? starttime : tst);
      rtime = max( rtime, meass[i]->real_time_ );
      speed = max( speed, meass[i]->speed_ );
      if( occstatus == Measurement::UnknownOccupancyStatus )
        occstatus = meass[i]->occupied_;
      else if( meass[i]->occupied_ != Measurement::UnknownOccupancyStatus )
        occstatus = max( occstatus, meass[i]->occupied_ );
    }
    
    
    ostr << "  <DetectorData>" << endline;
    if( !starttime.is_special() )
      ostr << "    <StartTime>" << UtilityFunctions::to_extended_iso_string(starttime) << "Z</StartTime>" << endline;
    if( rtime > 0.0f )
     ostr << "    <SampleRealTime>PT" << rtime << "S</SampleRealTime>" << endline;
    if( occstatus != Measurement::UnknownOccupancyStatus )
      ostr << "    <Occupied>" << (occstatus==Measurement::NotOccupied ? "0" : "1") << "</Occupied>" << endline;
    if( speed > 0.0f )
      ostr << "    <Speed Units=\"m/s\">" << speed << "</Speed>" << endline;

    string detsysname = measurement_location_name_;
    if( lane_number_ >= 0 )
      detsysname += "Lane" + std::to_string(lane_number_);
    if( inspection_.size() )
      detsysname += inspection_;
    if( detsysname.empty() )
      detsysname = "detector";
    
    ostr << "    <DetectorMeasurement Detector=\"" << detsysname << "\" DetectorType=\"Other\">" << endline;
    ostr << "      <SpectrumMeasurement>" << endline;
    ostr << "        <SpectrumAvailable>1</SpectrumAvailable>" << endline;
    
    for( const std::shared_ptr<const Measurement> meas : meass )
      meas->write_2006_N42_xml( ostr );
    
    ostr << "      </SpectrumMeasurement>" << endline;
    ostr << "    </DetectorMeasurement>" << endline;
    ostr << "  </DetectorData>" << endline;
  }
  
  //for( const std::shared_ptr<const Measurement> meas : measurements_ )
    //meas->write_2006_N42_xml( ostr );

  
  ostr << "</Measurement>" << endline;
  ostr << "</N42InstrumentData>" << endline;

  return !ostr.bad();
}//bool write_2006_N42( std::ostream& ostr ) const

bool Measurement::write_csv( std::ostream& ostr ) const
{
  const char *endline = "\r\n";

  ostr << "Energy, Data" << endline;

  for( size_t i = 0; i < gamma_counts_->size(); ++i )
    ostr << channel_energies_->at(i) << "," << gamma_counts_->operator[](i) << endline;

  ostr << endline;

  return !ostr.bad();
}//bool Measurement::write_csv( std::ostream& ostr ) const


bool MeasurementInfo::write_csv( std::ostream& ostr ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  for( const std::shared_ptr<const Measurement> meas : measurements_ )
    meas->write_csv( ostr );

  return !ostr.bad();
}//bool write_csv( std::ostream& ostr ) const


bool MeasurementInfo::load_from_pcf( std::istream &input )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  
  try
  {
    input.exceptions( ios::failbit | ios::badbit );

    input.seekg( 0, ios::end );
    const istream::pos_type eof_pos = input.tellg();
    input.seekg( orig_pos, ios::beg );
    const size_t filelen = static_cast<size_t>( 0 + eof_pos - orig_pos );
    
    if( filelen && (filelen < 512) )
      throw runtime_error( "File to small" );

    
    std::basic_string<char> fileid;  //There looks to be some structure in this, ex 'DHS       84567b0f-baa9-4176-9b8c-a066426063a9Secondary       Measured Data             SpecPortal                  ORTEC                       OSASP             Serial #10000033                      Anthony New Mexi                ?'
    fileid.resize( 256 );
    std::basic_string<char> header; //has some structure ( eg says 'DeviationPairsInFile')
    header.resize( 256 );
    
    
    input.read( &(fileid[0]), fileid.size() ); //256
    
    uint16_t NRPS;  //number of 256 byte records per spectrum
    memcpy( &NRPS, &(fileid[0]), 2 );
    const size_t bytes_per_record = 256 * static_cast<uint32_t>(NRPS);
    
    if( NRPS == 0 || (filelen && (NRPS*size_t(256)) > filelen) )
    {
      char buff[256];
      snprintf( buff, sizeof(buff), "Invalid number 256 segments per records, NRPS=%i", int(NRPS) );
      throw runtime_error( buff );
      //      throw runtime_error( "Invalid number 256 segments per records" );
    }
    
    //Usually expect fileid[2]=='D', fileid[3]=='H' fileid[4]=='S'
    const bool is_dhs_version = (fileid[2]=='D' || fileid[3]=='H' || fileid[4]=='S');
    bool goodheader = is_dhs_version;
    if( !goodheader )
      goodheader = (fileid.substr(2,3) == "   ");
    
    if( !goodheader )
    {
      //something like '783 - 03/06/15 18:10:28'
      const size_t pos = fileid.find( " - " );
      goodheader = (pos != string::npos && pos < 20
                     && fileid[pos+5]=='/' && fileid[pos+7]=='/'
                     && fileid[pos+9]=='/' && fileid[pos+11]==' ');
    }//if( !goodheader )
    
    if( !goodheader )
      throw runtime_error( "Unexpected fileID: '" + string( &fileid[0], &fileid[3] ) + "'" );

    input.read( &(header[0]), header.size() );  //512
    std::vector< std::pair<float,float> > deviation_pairs;
    
    if( header.find("DeviationPairsInFile") != string::npos )
    {
      if( header.find("DeviationPairsInFileCompressed") != string::npos )
      {
        //loop over columns (4)
        //  loop over panels (8)
        //    loop over MCAs (8)
        //      loop over deviation pairs (20)
        //         parse energy (4-byte floating point)
        //         parse offset (4-byte floating point)

        int16_t dev_pair_data[2*80*256/sizeof(float)];
        input.read( (char *)dev_pair_data, sizeof(dev_pair_data) );
        size_t last_pair = 0;
        const size_t devpairsize = sizeof(dev_pair_data)/sizeof(dev_pair_data[0]);
        for( size_t i = 0; i < devpairsize; ++i )
          if( dev_pair_data[i] != 0 )
            last_pair = i;
        if( (last_pair%2) == 1 )
          ++last_pair;
        for( size_t i = 0; i < last_pair; i+=2 )
          deviation_pairs.push_back( make_pair( dev_pair_data[i], dev_pair_data[i+1] ) );
        
        if( last_pair )
        {
          remarks_.push_back( "MeasurementInfo::load_from_pcf: deviation pairs are currently not being parsed from PCF files correctly" );
#if(PERFORM_DEVELOPER_CHECKS)
          log_developer_error( BOOST_CURRENT_FUNCTION, "MeasurementInfo::load_from_pcf: deviation pairs are currently not being parsed from PCF files correctly" );
#endif
        }
      }else
      {
        //There can be 128 different detectors defined by their column index, panel index, and MCA index
        //Each detector can have 20 deviation pairs.
        //loop over columns (2)   {these are portal columns - see N42 standard for which one where}
        //  loop over panels (8)  {panels in each column - see N42 standard for which one where}
        //    loop over MCAs (8)  {individual detectors in each panel}
        //      loop over deviation pairs (20)
        //         parse energy (4-byte floating point)
        //         parse offset (4-byte floating point)
        // float dev_pair_data[2/*columns*/][8/*panels*/][8/*MCAs*/][20/*pairs*/][2]  //5120
        
        float dev_pair_data[80*256/sizeof(float)];  //Is this entire region for deviation pairs?
                                                  //This region only here when header contains 'DeviationPairsInFile'

        input.read( (char *)dev_pair_data, sizeof(dev_pair_data) );

        size_t last_pair = 0;
        const size_t devpairsize = sizeof(dev_pair_data)/sizeof(dev_pair_data[0]);
        for( size_t i = 0; i < devpairsize; ++i )
          if( dev_pair_data[i] != 0.0 )
            last_pair = i;
        if( (last_pair%2) == 1 )
          ++last_pair;
        for( size_t i = 0; i < last_pair; i+=2 )
          deviation_pairs.push_back( make_pair( dev_pair_data[i], dev_pair_data[i+1] ) );
        
        if( last_pair )
        {
          remarks_.push_back( "MeasurementInfo::load_from_pcf: deviation pairs are currently not being parsed from PCF files correctly" );
#if(PERFORM_DEVELOPER_CHECKS)
          log_developer_error( BOOST_CURRENT_FUNCTION, "MeasurementInfo::load_from_pcf: deviation pairs are currently not being parsed from PCF files correctly" );
#endif
        }
      }
    }else
    {
      //If this is not the "DHS" version of a PCF file with the extended header
      //  information, after the file header contents are:
      //  byte offset, data type, description
      //  0          , int16_t  , Number of records per spectrum (NRPS)
      //  2          , char[3]  , Version
      //  5          , char[4]  , Energy calibration label (unused)
      //  9          , float[5] , Energy calibration
      //  Then I guess a bunch of garbage to get up to 356 bytes.
      //  Note that each spectrum record usually has its own calibration, so
      //  this one in the header can usually be ignored.
      const size_t current_pos = static_cast<size_t>( 0 + input.tellg() );
      input.seekg( current_pos-256, ios::beg );
    }//if( header.find("DeviationPairsInFile") != string::npos ) / else


    double latitude = -999.9, longitude = -999.9;
    
    if( is_dhs_version )
    {
      int16_t lanenumber, item_dist, occ_num;
      //const string last_modified = parse_pcf_field(fileid,5,7);
      uuid_ = parse_pcf_field(fileid,12,36);
      inspection_ = parse_pcf_field(fileid,48,16);  //Secondary, Primary
      memcpy( &lanenumber, &(fileid[64]), 2 );
      if( lanenumber > 0 )
        lane_number_ = lanenumber;
      const string measremark = parse_pcf_field(fileid,66,26);
      if(!measremark.empty())
        remarks_.push_back( measremark );
      instrument_type_ = parse_pcf_field(fileid,92,28);
      manufacturer_ = parse_pcf_field(fileid,120,28);
      instrument_model_ = parse_pcf_field(fileid,148,18);
      instrument_id_ = parse_pcf_field(fileid,166,18);
      const string item_description = parse_pcf_field(fileid,184,20);
      if(!item_description.empty())
        remarks_.push_back( "ItemDescription: " + measremark );
      measurement_location_name_ = parse_pcf_field(fileid,204,16);
      const string meas_coords = parse_pcf_field(fileid,220,16);
      vector<string> meas_coords_components;
      UtilityFunctions::split(meas_coords_components, meas_coords, " ,\t\r\n");
      if( meas_coords_components.size() > 2 )
      {
        //Totally untested as of 20170811 - beacuase I have never seen a PCF file with coordinates...
        //parse_deg_min_sec_lat_lon(...)
        //ortecLatOrLongStrToFlt()
        if( !toDouble( meas_coords_components[0], latitude )
           || !toDouble( meas_coords_components[0], longitude )
           || !Measurement::valid_latitude(latitude)
           || !Measurement::valid_longitude(longitude) )
        {
          latitude = longitude = -999.9;
#if(PERFORM_DEVELOPER_CHECKS)
          char buffer[256];
          snprintf( buffer, sizeof(buffer),
                   "PCF file had non empty coordinates string '%s', but didnt return valid coordinates", meas_coords.c_str() );
          log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
#endif
        }
      }//if( meas_coords_components.size() > 2 )
      
      memcpy( &item_dist, &(fileid[236]), 2 );
      if( item_dist > 0 )
        remarks_.push_back( "ItemToDetectorDistance: " + std::to_string(item_dist) + " cm" );
      
      memcpy( &occ_num, &(fileid[238]), 2 );
      if( occ_num > 0 )
        remarks_.push_back( "OccupancyNumber: " + std::to_string(occ_num) );
      
      const string cargo_type = parse_pcf_field(fileid,240,16);
      if( cargo_type.size() )
        remarks_.push_back( "CargoType: " + cargo_type );
    }//if( is_dhs_version )
    
    map<string,int> sample_numbers;
    map<MeasurementShrdPtr,int> sample_numbers_from_title;
    
    /* XXX: TODO: we should try to parse detector name from the title, and
     set had neutrons based on that
    */
    bool any_contained_neutron = false, all_contained_neutron = true;
    
    //If deviation pairs in file, we are at: 512+(20*256)=5632
    // else we are at 256

    while( input.good() && static_cast<size_t>(input.tellg()) < (filelen-256) )
    {
      const istream::pos_type specstart = input.tellg();

      basic_string<char> spectrum_title;  //ex: 'Background Aa1  Distance=250 cm'
      basic_string<char> additional_info; //ex '<T>HPGe 50%</T>' or '<T>Gamma</T>'
      basic_string<char> source_list;
      basic_string<char> collection_time;  //Formatted like: '2010-02-24T00:08:24.82Z'

      spectrum_title.resize(60);
      additional_info.resize(60);
      source_list.resize(60);
      collection_time.resize(23);

      input.read( &(spectrum_title[0]), spectrum_title.size() );   //60
      input.read( &(additional_info[0]), additional_info.size() ); //120
      input.read( &(source_list[0]), source_list.size() );         //180
      input.read( &(collection_time[0]), collection_time.size() ); //203

      
      size_t lastpos = spectrum_title.find_last_not_of( '\0' );
      if( lastpos != string::npos )
        spectrum_title = spectrum_title.substr( 0, lastpos + 1 );
      
      //Look for Det=...
      //Look for "Survey" or "Sample"
      
      
      lastpos = additional_info.find_last_not_of( '\0' );
      if( lastpos != string::npos )
        additional_info = additional_info.substr( 0, lastpos + 1 );
      
//      cerr << "source_list=" << source_list << endl;
//      cerr << "spectrum_title='" << spectrum_title << "'" << endl;
//      static_assert( sizeof(float) == 4, "Float must be 4 bytes" );
//      static_assert(std::numeric_limits<float>::digits >= 32);

      int32_t num_channel;
      char character_tag;
      vector<float> energy_cal_terms( 5, 0.0f );
      float live_time, true_time, halflife, molecular_weight,
            spectrum_multiplier, occupancy_flag, neutron_counts;
      
      input.read( &character_tag, 1 );                            //204
      input.read( (char *)&live_time, 4 );                        //208
      input.read( (char *)&true_time, 4 );                        //212
      input.read( (char *)&halflife, 4 );                         //216
      input.read( (char *)&molecular_weight, 4 );                 //220
      input.read( (char *)&spectrum_multiplier, 4 );              //224
      input.read( (char *)&energy_cal_terms[0], 4 );              //228
      input.read( (char *)&energy_cal_terms[1], 4 );              //232
      input.read( (char *)&energy_cal_terms[2], 4 );              //236
      input.read( (char *)&energy_cal_terms[3], 4 );              //240
      input.read( (char *)&energy_cal_terms[4], 4 );              //244
      input.read( (char *)&occupancy_flag, 4 );                   //248
      input.read( (char *)&neutron_counts, 4 );                   //252
      input.read( (char *)&num_channel, 4 );                      //256

      //we have now read 256 bytes for this record

      if( num_channel == 0 )
      {
        //lets advance to the next expected spectrum
        if( input.seekg( 0 + specstart + bytes_per_record, ios::beg ) )
          continue;
        else
          break;
      }//if( num_channel == 0 )

      //XXX - bellow, we are assuming we shouldnt ever get zero channels
      if( num_channel < 0 || num_channel>65536 )
      {
        char buffer[64];
        snprintf( buffer, sizeof(buffer),
                  "Invaid number of channels: %i", int(num_channel) );
        throw runtime_error( buffer );
      }//if( num_channel < 0 || num_channel>65536 )

      std::shared_ptr< vector<float> > channel_data = std::make_shared<vector<float> >( num_channel );
      input.read( (char *)&(channel_data->operator[](0)), 4*num_channel );

      const istream::pos_type specend = input.tellg();
      const size_t speclen = static_cast<size_t>( 0 + specend - specstart );
      if( speclen != bytes_per_record )
      {
        if( speclen > bytes_per_record )
          cerr << SRC_LOCATION << "\n\tUnexpected record length, expected "
               << (256*NRPS) << " but got length " << speclen << ", warning, am "
               << " forcing correct position in file" << endl;

        //For the last spectrum in the file may extend beyond the end of the
        //  file since NRPS may be larger than necessary to capture all the
        //  spectral information, so the extra space is left out of the file
        const size_t nextpos = static_cast<size_t>( 0 + specstart + bytes_per_record );
        if( nextpos > filelen )
          input.seekg( filelen, ios::beg );
        else
          input.seekg( 0 + specstart + bytes_per_record, ios::beg );
      }//if( speclen != (4*NRPS) )

      if( spectrum_multiplier > 1.0 )
      {
        for( float &f : *channel_data )
          f *= spectrum_multiplier;
      }//if( spectrum_multiplier!=0.0 && spectrum_multiplier!=1.0 )


      MeasurementShrdPtr meas = std::make_shared<Measurement>();
      measurements_.push_back( meas );
      meas->reset();

      meas->live_time_ = live_time;
      meas->real_time_ = true_time;
      meas->latitude_ = latitude;
      meas->longitude_ = longitude;
      
      const bool has_neutrons = (neutron_counts > 0.00000001);
      meas->contained_neutron_ = has_neutrons;
      any_contained_neutron = (any_contained_neutron || has_neutrons);
      all_contained_neutron = (all_contained_neutron && has_neutrons);
 
      for( const float f : *channel_data )
        meas->gamma_count_sum_ += f;
      meas->neutron_counts_.resize( 1 );
      meas->neutron_counts_[0] = neutron_counts;
      meas->neutron_counts_sum_ = neutron_counts;
      meas->speed_ = speed_from_remark( spectrum_title );
      meas->detector_name_ = detector_name_from_remark( spectrum_title );

      if( !sample_numbers.count(meas->detector_name_) )
        sample_numbers[meas->detector_name_] = 0;
      meas->sample_number_ = ++sample_numbers[meas->detector_name_];
      
      if( occupancy_flag == 0.0f )  //should be zero or one, but lets be loose
        meas->occupied_ = Measurement::NotOccupied;
      else if( occupancy_flag == 1.0f )
        meas->occupied_ = Measurement::Occupied;
      else
        meas->occupied_ = Measurement::UnknownOccupancyStatus;  //shouldnt ever get here, but file values are whack sometimes.
      
      const int titlesamplenum = sample_num_from_remark( spectrum_title );
      if( titlesamplenum >= 0 )
        sample_numbers_from_title[meas] = titlesamplenum;
      
      sample_numbers[meas->detector_name_] = std::max( meas->sample_number_,
                                        sample_numbers[meas->detector_name_] );
      
      
      
      meas->start_time_ = time_from_string( collection_time.c_str() );
      
       //XXX test for Background below not tested
      if( UtilityFunctions::icontains( spectrum_title, "Background" ) )
        meas->source_type_ = Measurement::Background;
      else //if( spectrum_title.find("Foreground") != string::npos )
        meas->source_type_ = Measurement::Foreground;
       //else meas->source_type_ = Measurement::UnknownSourceType

      trim( spectrum_title );
//      if( spectrum_title.size() )
//        meas->remarks_.push_back( spectrum_title );
      meas->title_ = spectrum_title;

      trim( additional_info );
      if( additional_info.size() )
        meas->remarks_.push_back( additional_info );

//      meas->cambio_tag_char_ = character_tag;
      
      while( energy_cal_terms.size() && (energy_cal_terms.back()==0.0f) )
        energy_cal_terms.erase( energy_cal_terms.begin() + energy_cal_terms.size() - 1 );
      
      meas->energy_calibration_model_ = Measurement::FullRangeFraction;
      meas->calibration_coeffs_ = energy_cal_terms;

      meas->deviation_pairs_ = deviation_pairs; //this is not correct if there is more than one detector in the file!
      meas->gamma_counts_ = channel_data;
    }//while( stream.good() && stream.tellg() < size )

    if( any_contained_neutron && !all_contained_neutron )
    {
      for( auto &p : measurements_ )
        p->contained_neutron_ = true;
    }else if( !any_contained_neutron )
    {
      for( auto &p : measurements_ )
      {
        p->neutron_counts_.clear();
        p->neutron_counts_sum_ = 0.0;
      }
    }//if( we had some neutrons ) / else
    

    //Assign sample numbers if the titles all specified the sample numbers
    if( sample_numbers_from_title.size() == measurements_.size() )
    {
      for( auto &p : measurements_ )
        p->sample_number_ = sample_numbers_from_title[p];
    }
    
    //if( any_contained_neutron && !all_contained_neutron )

    if( measurements_.empty() )
      throw runtime_error( "Didnt read in any Measurments" );

    cleanup_after_load();
  }catch( std::exception & /*e*/ )
  {
    input.clear();
    input.seekg( orig_pos, ios::beg );
//    cerr << SRC_LOCATION << "\n\tCaught:" << e.what() << endl;
    reset();
    return false;
  }//try / catch

  return true;
}//bool load_from_pcf( std::istream& istr )



bool MeasurementInfo::load_from_chn( std::istream &input )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );

  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  input.seekg( 0, ios::end );
  const istream::pos_type eof_pos = input.tellg();
  input.seekg( orig_pos, ios::beg );
  const size_t size = static_cast<size_t>( 0 + eof_pos - orig_pos );

  try
  {
    if( size < 548 )
      throw runtime_error( "File to small to be a CHN file." );
    
    const size_t header_size = 32;

    vector<char> buffer( header_size );
    if( !input.read( &buffer[0], header_size ) )
      throw runtime_error( SRC_LOCATION + " Error reading from file stream" );

    int16_t firstval;
    memcpy( &firstval, &(buffer[0]), sizeof(uint16_t) );
    
    if( firstval != -1 )
      throw runtime_error( "Invalid first value" );

    uint16_t firstchannel, numchannels;
    memcpy( &firstchannel, &(buffer[28]), sizeof(uint16_t) );
    memcpy( &numchannels, &(buffer[30]), sizeof(uint16_t) );
    
//    if( numchannels )
//      numchannels = (size - 32 - 516) / 4;
    
//    const bool isPowerOfTwo = ((numchannels != 0) && !(numchannels & (numchannels - 1)));
//    const bool evenNumChnnel = ((numchannels%2)==0) || (numchannels==16383);
//    if( !evenNumChnnel || numchannels<2 ) // || !isPowerOfTwo
//      throw runtime_error( "CHN file with unexpected number of data channels: " + std::to_string(numchannels) );

#if(PERFORM_DEVELOPER_CHECKS)
    if( firstchannel != 0 )
    {
      //if( firstchannel==1 ), we should shift the bin contents over by 1 to the left
      char buffer[256];
      snprintf( buffer, sizeof(buffer),
                "Found a first channel offset of %i", int(firstchannel) );
      log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
      if( firstchannel != 1 )
        cerr << "firstchannel==" << firstchannel << endl;  //I've never ran into this, at least in the corpus of files I have
    }
#endif
    
    if( size < (header_size+4*numchannels) )
      throw runtime_error( "CHN Filesize smaller than expected" );

    //first non-zero spectrum counts at channel 120, so I assume first channel data
    //is in channel 32, 36, or 40 (PeakEasy has 22 leading zeros, Cambio 21)

    const string monthstr(  &(buffer[18]), &(buffer[18]) + 3 );
    const string daystr(    &(buffer[16]), &(buffer[16]) + 2 );
    const string yearstr(   &(buffer[21]), &(buffer[21]) + 2 );
    const bool isNot2k = (buffer[23]=='0');  //I think PeakEasy defaults to putting a '\0' in this spot, so well assume year > 2k unless actually properly noted
    const string hourstr(   &(buffer[24]), &(buffer[24]) + 2 );
    const string minutestr( &(buffer[26]), &(buffer[26]) + 2 );
    const string secondstr( &(buffer[6]), &(buffer[6]) + 2 );

    uint32_t realTimeTimesFifty, liveTimeTimesFifty;
    memcpy( &realTimeTimesFifty, &(buffer[8]), sizeof(uint32_t) );
    memcpy( &liveTimeTimesFifty, &(buffer[12]), sizeof(uint32_t) );
    
    const float realTime = realTimeTimesFifty / 50.0f;
    const float liveTime = liveTimeTimesFifty / 50.0f;

    double gamma_sum = 0.0;
    ShrdFVecPtr channel_data( new vector<float>(numchannels) );
    vector<float> &channel_data_ref = *channel_data;
    
    if( numchannels > 2 )
    {
      vector<uint32_t> int_channel_data( numchannels );
      if( !input.read( (char *)(&int_channel_data[0]), 4*numchannels ) )
        throw runtime_error( SRC_LOCATION + " Error reading from file stream" );

      int_channel_data[0] = 0;
      int_channel_data[1] = 0;
      int_channel_data[numchannels-2] = 0;
      int_channel_data[numchannels-1] = 0;
      
      for( uint16_t i = 0; i < numchannels; ++i )
      {
        gamma_sum += int_channel_data[i];
        channel_data_ref[i] = static_cast<float>( int_channel_data[i] );
      }//for( size_t i = 0; i < numchannels; ++i )
    }//if( numchannels > 2 )

    const istream::pos_type current_pos = input.tellg();
    const size_t bytes_left = static_cast<size_t>( 0 + eof_pos - current_pos );


    buffer.resize( bytes_left );
    if( !input.read( &buffer[0], bytes_left ) )
      throw runtime_error( SRC_LOCATION + " Error reading from file stream" );

    int16_t chntype;
    memcpy( &chntype, &(buffer[0]), sizeof(int16_t) );
    
#if(PERFORM_DEVELOPER_CHECKS)
    //Files such as ref985OS89O82 can have value chntype=-1
    if( chntype != -102 && chntype != -101 )
    {
      stringstream msg;
      msg << "Found a chntype with unexpected value: " << chntype;
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
    }
#endif
    
    float calibcoefs[3] = { 0.0f, 0.0f, 0.0f };
    if( bytes_left >= 12 )
    {
      if( chntype == -102 )
        memcpy( calibcoefs, &(buffer[4]), 3*sizeof(float) );
      else
        memcpy( calibcoefs, &(buffer[4]), 2*sizeof(float) );
//      float FWHM_Zero_32767 = *(float *)(&(buffer[16]));
//      float FWHM Slope = *(float *)(&(buffer[20]));
    }//if( we potentially have calibration info )

    //at 256 we have 1 byte for DetDescLength, and then 63 bytes for DetDesc
    string detdesc;
    uint8_t DetDescLength;
    memcpy( &DetDescLength, &(buffer[256]), 1 );
    const size_t det_desc_index = 257;
    const size_t enddet_desc = det_desc_index + DetDescLength;
    if( DetDescLength && enddet_desc < bytes_left && DetDescLength < 64 )
    {
      detdesc = string( &(buffer[det_desc_index]), &(buffer[enddet_desc]) );
      trim( detdesc );
    }//if( title_index < bytes_left )
    
    
    string title;
    uint8_t SampleDescLength;
    memcpy( &SampleDescLength, &(buffer[320]), 1 );
    const size_t title_index = 321;
    const size_t endtitle = title_index + SampleDescLength;
    if( SampleDescLength && endtitle < bytes_left && SampleDescLength < 64 )
    {
      title = string( &(buffer[title_index]), &(buffer[endtitle]) );
      trim( title );
    }//if( title_index < bytes_left )

    MeasurementShrdPtr meas = std::make_shared<Measurement>();
    meas->live_time_ = liveTime;
    meas->real_time_ = realTime;
    meas->gamma_count_sum_ = gamma_sum;

    if( (fabs(calibcoefs[0])<1.0E-12 && fabs(calibcoefs[1])<1.0E-12)
        || (fabs(calibcoefs[0])<1.0E-12 && fabs(calibcoefs[1]-1.0)<1.0E-8) )
    {
      calibcoefs[0] = calibcoefs[1] = calibcoefs[2] = 0.0;
      meas->calibration_coeffs_.clear();
//      meas->calibration_coeffs_.push_back( 0.0f );
//      meas->calibration_coeffs_.push_back( 1.0f );
    }else
    {
      meas->energy_calibration_model_ = Measurement::Polynomial;
      meas->calibration_coeffs_.push_back( calibcoefs[0] );
      meas->calibration_coeffs_.push_back( calibcoefs[1] );
      if( calibcoefs[2] != 0.0f )
        meas->calibration_coeffs_.push_back( calibcoefs[2] );
    }//if( calibcoefs[0]==0.0 && calibcoefs[1]==1.0 )


    if( channel_data && channel_data->size() )
      meas->gamma_counts_ = channel_data;

    const string datestr = daystr + "-" + monthstr
                           + (isNot2k ? "-19" : "-20") + yearstr
                           + " " + hourstr + ":" + minutestr + ":" + secondstr;

    meas->start_time_ = time_from_string( datestr.c_str() );

    if( title.size() )
      meas->title_ = title;

    if( detdesc.size() )
      meas->remarks_.push_back( "Detector Description: " + detdesc );
    
    measurements_.push_back( meas );

    cleanup_after_load();
  }catch( std::runtime_error &e )
  {
    input.clear();
    input.seekg( orig_pos, ios::beg );
    //cerr << SRC_LOCATION << "\n\tCaught:" << e.what() << endl;
    reset();
    return false;
  }//try / catch


  return true;
}//bool load_from_chn( std::istream &input )


bool MeasurementInfo::write_integer_chn( ostream &ostr, set<int> sample_nums,
                                         const set<int> &det_nums ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  if( sample_nums.empty() )
    sample_nums = sample_numbers_;
  
  const size_t ndet = detector_numbers_.size();
  vector<bool> detectors( ndet, true );
  if( !det_nums.empty() )
  {
    for( size_t i = 0; i < ndet; ++i )
      detectors[i] = (det_nums.count(detector_numbers_[i]) != 0);
  }//if( det_nums.empty() )
  
  MeasurementShrdPtr summed = sum_measurements( sample_nums, detectors );
  
  if( !summed || !summed->gamma_counts() )
    return false;
  
  const ShrdConstFVecPtr &fgammacounts = summed->gamma_counts();
  
  int16_t val16 = -1;
  ostr.write( (const char *)&val16, sizeof(int16_t) );
  //index=2

  val16 = 0; //McaNumber
  ostr.write( (const char *)&val16, sizeof(int16_t) );
  //index=4
  
  val16 = 1; //Segment, set to 1 in UMCBI
  ostr.write( (const char *)&val16, sizeof(int16_t) );
  //index=6
  
  char buffer[256];
  
  const boost::posix_time::ptime &starttime = summed->start_time_;
  if( starttime.is_special() )
    buffer[0] = buffer[1] = '0';
  else
    snprintf( buffer, sizeof(buffer), "%02d", int(starttime.time_of_day().seconds()) );
  ostr.write( buffer, 2 );
  //index=8
  
  const uint32_t realTimeTimesFifty = static_cast<uint32_t>(50.0f*summed->real_time_);
  const uint32_t liveTimeTimesFifty = static_cast<uint32_t>(50.0f*summed->live_time_);
  ostr.write( (const char *)&realTimeTimesFifty, 4 );
  ostr.write( (const char *)&liveTimeTimesFifty, 4 );
  //index=16
  
  if( starttime.is_special() )
    buffer[0] = buffer[1] = '0';
  else
    snprintf( buffer, sizeof(buffer), "%02d", int(starttime.date().day_number()) );
  ostr.write( buffer, 2 );
  //index=18

  const char *monthstr = "   ";
  
  try
  {
    switch( starttime.date().month().as_enum() )
    {
      case boost::gregorian::Jan: monthstr = "Jan"; break;
      case boost::gregorian::Feb: monthstr = "Feb"; break;
      case boost::gregorian::Mar: monthstr = "Mar"; break;
      case boost::gregorian::Apr: monthstr = "Apr"; break;
      case boost::gregorian::May: monthstr = "May"; break;
      case boost::gregorian::Jun: monthstr = "Jun"; break;
      case boost::gregorian::Jul: monthstr = "Jul"; break;
      case boost::gregorian::Aug: monthstr = "Aug"; break;
      case boost::gregorian::Sep: monthstr = "Sep"; break;
      case boost::gregorian::Oct: monthstr = "Oct"; break;
      case boost::gregorian::Nov: monthstr = "Nov"; break;
      case boost::gregorian::Dec: monthstr = "Dec"; break;
      case boost::gregorian::NotAMonth:
      case boost::gregorian::NumMonths:
       break;
    }//switch( starttime.date().month().as_enum() )
  }catch(...)
  {
    //Here when month is invalid...
  }
  ostr.write( monthstr, 3 );
  //index=21
  
  if( starttime.is_special() )
    buffer[0] = buffer[1] = '0';
  else
    snprintf( buffer, sizeof(buffer), "%02d", int(starttime.date().year()%100) );
  ostr.write( buffer, 2 );
  //index=23
  
  if( !starttime.is_special() && starttime.date().year() >= 2000 )
    ostr.write( "1", 1 );
  else
    ostr.write( "0", 1 );
  //index=24
  
  if( starttime.is_special() )
    buffer[0] = buffer[1] = '0';
  else
    snprintf( buffer, sizeof(buffer), "%02d", int(starttime.time_of_day().hours()) );
  ostr.write( buffer, 2 );
  //index=26

  
  if( starttime.is_special() )
    buffer[0] = buffer[1] = '0';
  else
    snprintf( buffer, sizeof(buffer), "%02d", int(starttime.time_of_day().minutes()) );
  ostr.write( buffer, 2 );
  //index=28

  uint16_t firstchannel = 0;
  ostr.write( (const char *)&firstchannel, 2 );
  //index=30
  
  const uint16_t numchannels = static_cast<uint16_t>( fgammacounts->size() );
  ostr.write( (const char *)&numchannels, 2 );
  //index=32
  
  //Not actually sure if we want to write the channel data at index 32, 34, 36, or 40...
  //  Also, there may be a need to shift the channels by one left or right.
  //  Also, not certain about values in first channel or two...
  vector<uint32_t> intcounts( numchannels );
  for( uint16_t i = 0; i < numchannels; ++i )
    intcounts[i] = static_cast<uint32_t>( (*fgammacounts)[i] + 0.5f );
  ostr.write( (const char *)&intcounts[0], numchannels*4 );
  
  vector<float> calibcoef = summed->calibration_coeffs();
  
  switch( summed->energy_calibration_model() )
  {
    case Measurement::Polynomial:
      break;
    case Measurement::FullRangeFraction:
      calibcoef = fullrangefraction_coef_to_polynomial( calibcoef, fgammacounts->size() );
      break;
    case Measurement::LowerChannelEdge:
    case Measurement::UnknownEquationType:
      calibcoef.clear();
      break;
  }//switch( summed->energy_calibration_model() )
  
  while( calibcoef.size() < 3 )
    calibcoef.push_back( 0.0f );
  
  int16_t chntype = -102;
  ostr.write( (const char *)&chntype, 2 );
  
  buffer[0] = buffer[1] = '0';
  ostr.write( buffer, 2 );
  ostr.write( (const char *)&calibcoef[0], 3*sizeof(float) );
  
  float fummyf = 0.0f;
  ostr.write( (const char *)&fummyf, sizeof(float) );
  ostr.write( (const char *)&fummyf, sizeof(float) );
  ostr.write( (const char *)&fummyf, sizeof(float) );
  
  for( size_t i = 0; i < 228; ++i )
    ostr.write( "\0", 1 );
  
  
  
  string detdesc = summed->title_;
  for( const string &remark : remarks_ )
  {
    if( UtilityFunctions::starts_with( remark, "Detector Description: " ) )
      detdesc = " " + remark.substr( 22 );
  }
  
  trim( detdesc );
  
  if( detdesc.length() > 63 )
    detdesc = detdesc.substr(0,63);
  
  uint8_t len = static_cast<uint8_t>( detdesc.size() );
  ostr.write( (const char *)&len, 1 );
  ostr.write( detdesc.c_str(), len+1 );
  for( size_t i = len+1; i < 63; ++i )
    ostr.write( "\0", 1 );
  
  string titlestr;
  if( measurements_.size()==1 )
    titlestr = measurements_[0]->title_;
  
  if( titlestr.length() > 63 )
    titlestr = titlestr.substr(0,63);
  len = static_cast<uint8_t>( titlestr.size() );
  ostr.write( (const char *)&len, 1 );
  ostr.write( titlestr.c_str(), len+1 );
  for( size_t i = len+1; i < 63; ++i )
    ostr.write( "\0", 1 );
  
  for( size_t i = 0; i < 128; ++i )
    ostr.write( "\0", 1 );
  
  return true;
}//bool write_integer_chn(...)


void Measurement::set_info_from_avid_mobile_txt( std::istream &istr )
{
  //There is a variant of refQQZGMTCC93, RSL mobile system ref8T2SZ11TQE
  
  using UtilityFunctions::safe_get_line;
  using UtilityFunctions::split_to_floats;
  
  const istream::pos_type orig_pos = istr.tellg();
  
  try
  {
    string line;
    if( !UtilityFunctions::safe_get_line(istr, line) )
      throw runtime_error(""); //"Failed getting first line"
  
    if( line.size() < 8 || line.size() > 100 )
      throw runtime_error(""); //"First line not reasonable length"
    
    const size_t first_invalid_char = line.substr(0,8).find_first_not_of( "0123456789 ,\r\n\t+-e." );
    
    if( first_invalid_char != string::npos )
      throw runtime_error( "" ); //"Invalid character in first 8 characters"
    
    vector<float> fline;
    if( !split_to_floats(line, fline) || fline.size()!=4 )
      throw runtime_error( "" ); //"First line not real time then calibration coefs"
    
    const vector<float> eqn( fline.begin() + 1, fline.end() );
    const float realtime = fline[0];
    
    if( realtime < -FLT_EPSILON )
      throw runtime_error( "" ); //"First coefficient not real time"
    
    if( !safe_get_line(istr, line) )
      throw runtime_error(""); //"Failed getting second line"
    
    if( !split_to_floats(line, fline) )
      throw runtime_error( "" ); //"Second line not floats"
    
    if( fline.size() < 127 && fline.size() != 2 )
      throw runtime_error( "" ); //"Invalid second line"
    
    //If we got here, this is probably a valid file
    std::shared_ptr< vector<float> > counts
                                       = std::make_shared< vector<float> >();
    
    if( fline.size() >= 127 )
    {
      //Second line is CSV of channel counts
      if( UtilityFunctions::safe_get_line(istr, line) && line.size() )
        throw runtime_error(""); //"Only expected two lines"
      
      counts->swap( fline );
    }else
    {
      //Rest of file is \t seperated column with two columns per line
      //  "channel\tcounts"
      float channelnum = fline[0];
      const float counts0 = fline[1];
    
      if( fabs(channelnum) > FLT_EPSILON && fabs(channelnum - 1.0) > FLT_EPSILON )
        throw runtime_error( "" ); //"First column doesnt refer to channel number"
    
      if( counts0 < -FLT_EPSILON )
        throw runtime_error( "" ); //"Second column doesnt refer to channel counts"

      channelnum = channelnum - 1.0f;
      istr.seekg( orig_pos, ios::beg );
      UtilityFunctions::safe_get_line( istr, line );
    
      while( safe_get_line( istr, line ) )
      {
        trim( line );
        if( line.empty() ) //Sometimes file will have a newline at the end of the file
          continue;
        
        if( !split_to_floats(line, fline) || fline.size() != 2 )
          throw runtime_error( "" ); //"Unexpected number of fields on a line"
        
        if( fabs(channelnum + 1.0f - fline[0]) > 0.9f /*FLT_EPSILON*/ )
          throw runtime_error( "" ); //"First column is not channel number"
        
        channelnum = fline[0];
        counts->push_back( fline[1] );
      }//while( UtilityFunctions::safe_get_line( istr, line ) )
    }//if( fline.size() >= 127 )
    
    const size_t nchannel = counts->size();
    if( nchannel < 127 )
      throw runtime_error(""); //"Not enought channels"
    
    const vector< pair<float,float> > devpairs;
    const bool validcalib
           = MeasurementInfo::calibration_is_valid( Polynomial, eqn, devpairs,
                                                    nchannel );
    if( !validcalib )
      throw runtime_error( "" ); //"Invalid calibration"
    
//    real_time_ = realtime;
    live_time_ = realtime;
    contained_neutron_ = false;
    deviation_pairs_.clear();
    calibration_coeffs_ = eqn;
    energy_calibration_model_ = Polynomial;
    neutron_counts_.clear();
    gamma_counts_ = counts;
    neutron_counts_sum_ = gamma_count_sum_ = 0.0;
    for( const float f : *counts )
      gamma_count_sum_ += f;
  }catch( std::exception &e )
  {
    istr.seekg( orig_pos, ios::beg );
    throw e;
  }
}//void set_info_from_avid_mobile_txt( std::istream& istr )


void Measurement::set_info_from_txt_or_csv( std::istream& istr )
{
  const istream::pos_type orig_pos = istr.tellg();
  
  errno = 0;
  
  try
  {
    set_info_from_avid_mobile_txt( istr );
    return;
  }catch(...)
  {
  }
  
  //I feel as though this function can be improved in terms of being more robust
  //  to reading input, as well as shortened or re-factored
  //Also, I hacked to make this function very quickly, so I'm sure the code is
  //  even more so crap that the rest in this file
  reset();
  const int kChannel = 0, kEnergy = 1, kCounts = 2;//, kSecondRecord = 3;
  float energy_units = 1.0f;
  
  map<size_t,int> column_map;
  vector<string>::const_iterator pos;

  string line;
  
  size_t nlines_used = 0, nlines_total = 0;
  
  const size_t maxlen = 1024*1024; //should be long enough for even the largest spectra
  while( UtilityFunctions::safe_get_line(istr, line, maxlen) )
  {
    if( line.size() > (maxlen-5) )
      throw runtime_error( "Found to long of line" );
    
    trim( line );
    to_lower( line );
    
    if( line.empty() )
      continue;

    ++nlines_total;
    
    vector<string> split_fields, fields;
    const char *delim = (line.find(',') != string::npos) ? "," : "\t, ";
    UtilityFunctions::split( split_fields, line, delim );

    fields.reserve( split_fields.size() );
    for( string s : split_fields )
    {
      trim( s );
      if( !s.empty() )
        fields.push_back( s );
    }//for( string s : split_fields )

    const size_t nfields = fields.size();

    if( !nfields )
      continue;

    if( isdigit(fields.at(0).at(0)) )
    {
      if( column_map.empty() )
      {
        if( nfields==1 )
        {
          column_map[0] = kCounts;
        }if( nfields==2 && isdigit(fields.at(1).at(0)) )
        {
          column_map[0] = kEnergy;
          column_map[1] = kCounts;
        }else if( nfields<9
                  && isdigit(fields.at(1).at(0)) && isdigit(fields.at(2).at(0)) )
        {
          column_map[0] = kChannel;
          column_map[1] = kEnergy;
          column_map[2] = kCounts;
        }else
        {
          throw runtime_error( string("unrecognized line that started with digit '")
                                + fields[0][0] + string("'") );
        }
      }//if( column_map.empty() )

      if( fields.size() == 4 )
      {
        vector<float> cals;
        if( UtilityFunctions::split_to_floats( line.c_str(), line.size(), cals ) )
        {
          //refY2EF53S0BD
          vector<float> eqn;
          if( cals.size() )
            eqn.insert( eqn.end(), cals.begin()+1, cals.end() );
          
          if( eqn.size()>=3 && fabs(eqn[0]) < 3000.0f && eqn[1] >= 0.0f )
          {
            //I think cals[0] might be real time
//            cerr << "cals[0]=" << cals[0] << endl;
            const vector< pair<float,float> > devpairs;
            
            const istream::pos_type current_pos = istr.tellg();
            
            string channeldata;
            if( UtilityFunctions::safe_get_line(istr, channeldata, maxlen) )
            {
              ++nlines_total;
              const istream::pos_type post_pos = istr.tellg();
              istr.seekg( 0, ios::end );
              const istream::pos_type eof_pos = istr.tellg();
              istr.seekg( post_pos, ios::beg );
              if( post_pos == eof_pos )
                istr.setstate( ios::eofbit );
              
              std::shared_ptr<vector<float> > channels( new vector<float>() );
              if( UtilityFunctions::split_to_floats( channeldata.c_str(), channeldata.size(), *channels ) )
              {
                if( post_pos == eof_pos )
                {
                  ++nlines_used;
                  const size_t nchan = channels->size();
                  const bool validCalib
                         = MeasurementInfo::calibration_is_valid( Polynomial, eqn,
                                                              devpairs, nchan );
                  
                  if( validCalib && nchan >= 128 )
                  {
                    gamma_counts_ = channels;
                    energy_calibration_model_ = Polynomial;
                    calibration_coeffs_ = eqn;
                    channel_energies_ = polynomial_binning( eqn, nchan, devpairs );
                    
                    break;
                  }//if( some reasonalbe number of channels )
                }else if( channels->size() == 2 )
                {
                  // refV6GHP7WTWX
                  channels->clear();
                  
                  string str;
                  while( UtilityFunctions::safe_get_line(istr, str, maxlen) )
                  {
                    ++nlines_total;
                    trim(str);
                    if( str.empty() )
                      continue;
                    
                    vector<float> vals;
                    if( !UtilityFunctions::split_to_floats( str.c_str(), str.size(), vals ) )
                    {
                      channels.reset();
                      break;
                    }
                    if( vals.size() != 2 )
                    {
                      channels.reset();
                      break;
                    }
                    
                    ++nlines_used;
                    channels->push_back( vals[1] );
                  }//
                  
                  if( !!channels )
                  {
                    const bool validCalib
                      = MeasurementInfo::calibration_is_valid( Polynomial, eqn, devpairs,
                                                              channels->size() );
                    
                    if( validCalib && channels->size() >= 64 )
                    {
                      ++nlines_used;
                      
                      live_time_ = cals[0];
                      gamma_counts_ = channels;
                      energy_calibration_model_ = Polynomial;
                      const size_t ncalcoef = size_t((cals[3]==0.0f) ? 2 : 3);
                      calibration_coeffs_.resize( ncalcoef );
                      for( size_t i = 0; i < ncalcoef; ++i )
                        calibration_coeffs_[i] = cals[1+i];
                      channel_energies_ = polynomial_binning( calibration_coeffs_,
                                                             channels->size(),
                                                             DeviationPairVec() );
                      
                      break;
                    }//if( some reasonalbe number of channels )
                  }//if( !!channels )
                }//if( there were exactly two lines ) / else
              }//if( could split the second line into floats )
            
            }//if( we could get a second line )
            
            //If we didnt 'break' above, we werent successful at reading the
            //  file
            istr.seekg( current_pos, ios::beg );
          }//if( potentially calibration data )
        }//if( the line was made of 4 numbers )
      }//if( fields.size() == 4 )
      
      ShrdFVecPtr energies( new vector<float>() ), counts( new vector<float>());
      std::shared_ptr<vector<int> > channels( new vector<int>() );

      //After we hit a line that no longer starts with numbers, we actually want
      // to leave istr at the beggining of that line so if another spectrum
      // comes aftwards, we wont lose its first line of information
      istream::pos_type position = istr.tellg();
      
      do
      {
        ++nlines_total;
        
        if( line.size() > (maxlen-5) )
          throw runtime_error( "Found to long of line" );
        
        vector<string> split_fields, fields;
        trim( line );
        split( split_fields, line, "\t, " );

        fields.reserve( split_fields.size() );
        for( const string &s : split_fields )
          if( !s.empty() )
            fields.push_back( s );


        if( fields.empty() )
          continue;

        if( !isdigit( fields.at(0).at(0) ) )
        {
          istr.seekg( position, ios::beg );
          break;
        }

        int channel = 0;
        float energy = 0.0f, count = 0.0f;//, count2 = 0.0f;
        for( size_t col = 0; col < fields.size(); ++col )
        {
          if( column_map.count(col) )
          {
            switch( column_map[col] )
            {
              case kChannel:       channel = atoi(fields[col].c_str()); break;
              case kEnergy:        energy  = static_cast<float>(atof(fields[col].c_str())); break;
              case kCounts:        count   = static_cast<float>(atof(fields[col].c_str())); break;
//              case kSecondRecord: count2  = static_cast<float>(atof(fields[col].c_str()));  break;
            }//switch( column_map[col] )
          }//if( column_map.count(col) )
        }//for( size_t col = 0; col < fields.size(); ++col )

        if( IsNan(energy) || IsInf(energy) )
          continue;
        if( IsNan(count) || IsInf(count) /*|| IsNan(count2) || IsInf(count2)*/ )
          continue;
        
//        if( errno )
//          throw runtime_error( "Error converting to float" );
        
        energy *= energy_units;
        
        if( (energies->size() && (energies->back() > energy) )
            || ( channels->size() && (channels->back() > channel) ) )
        {
          throw runtime_error( "Found decreasing energy" );
        }//if( energies->size() && (energies->back() > energy) )
        
        ++nlines_used;
        energies->push_back( energy );
        counts->push_back( count );
        channels->push_back( channel );
        position = istr.tellg();
      }while( UtilityFunctions::safe_get_line( istr, line, maxlen ) );

      if( counts->empty() )
        throw runtime_error( "Didnt find and channel counts" );
      
      gamma_counts_ = counts;

      if( energies->size() && (energies->back()!=0.0f) )
      {
        energy_calibration_model_ = LowerChannelEdge;
        calibration_coeffs_ = *energies;
//        channel_energies_ = energies;
      }else //if( channels->size() && (channels->back()!=0) )
      {
        if( (energy_calibration_model_ == Polynomial)
            && (calibration_coeffs_.size() > 1) && (calibration_coeffs_.size() < 10) )
        {
          //nothing to do here
        }else
        {
          energy_calibration_model_ = Polynomial;
          calibration_coeffs_.resize( 2 );
          calibration_coeffs_[0] = 0.0f;
          calibration_coeffs_[1] = 3000.0f / counts->size();
//        channel_energies_ = polynomial_binning( calibration_coeffs_, counts->size(),
//                                                           DeviationPairVec() );
        }
      }//if( energies ) / else channels

      break;
    }else if( column_map.empty()
        && (  istarts_with( fields[0], "channel" )
            || istarts_with( fields[0], "counts" )
            || istarts_with( fields[0], "data" )
            || istarts_with( fields[0], "energy" )
            || istarts_with( fields[0], "Ch" ) ) )
    {
      ++nlines_used;
      
      for( size_t i = 0; i < nfields; ++i )
      {
        if( starts_with( fields[i], "channel" )
           || starts_with( fields[i], "ch" ) )
        {
          column_map[i] = kChannel;
        }else if( starts_with( fields[i], "energy" )
                || starts_with( fields[i], "en" ) )
        {
          column_map[i] = kEnergy;
          if( UtilityFunctions::contains( fields[i], "mev" ) )
            energy_units = 1000.0f;
        }else if( starts_with( fields[i], "counts" )
                  || starts_with( fields[i], "data" )
                  || starts_with( fields[i], "selection" ) )
        {
//          bool hasRecordOne = false;
//          for( map<size_t,int>::const_iterator iter = column_map.begin();
//              iter != column_map.end(); ++iter )
//          {
//            hasRecordOne |= (iter->second == kCounts);
//          }
          
//          if( hasRecordOne )
//            column_map[i] = kSecondRecord;
//          else
            column_map[i] = kCounts;
        }
      }//for( size_t i = 0; i < nfields; ++i )
      
    }else if( starts_with( fields[0], "remark" ) )
    {
      ++nlines_used;
      bool used = false;
      pos = std::find( fields.begin(), fields.end(), "starttime" );
      if( (pos != fields.end()) && ((pos+1) != fields.end()) )
      {
        used = true;
        start_time_ = time_from_string( (pos+1)->c_str() );
      }
      pos = std::find( fields.begin(), fields.end(), "livetime" );
      if( (pos != fields.end()) && ((pos+1) != fields.end()) )
      {
        used = true;
        live_time_ = time_duration_in_seconds( *(pos+1) );
      }
      pos = std::find( fields.begin(), fields.end(), "realtime" );
      if( (pos != fields.end()) && ((pos+1) != fields.end()) )
      {
        used = true;
        real_time_ = time_duration_in_seconds( *(pos+1) );
      }

      try
      {
        if( sample_number_ < 0 )
        {
          sample_number_ = sample_num_from_remark( line );
          used |= (sample_number_ > -1);
        }
      }catch(...){}
      
      try
      {
        if( speed_ == 0.0 )
        {
         speed_ = speed_from_remark( line );
         used |= (speed_ != 0.0);
        }
      }catch(...){}
      
      try
      {
        if( detector_name_.empty() )
        {
          detector_name_ = detector_name_from_remark( line );
          used |= !detector_name_.empty();
        }
      }catch(...){}
      
      if( !used )
      {
        string::size_type pos = line.find_first_of( " :\t" );
        if( pos < line.size() )
        {
          string remark = line.substr( pos + 1 );
          pos = remark.find_first_not_of( " :\t" );
          remarks_.push_back( remark.substr( pos ) );
        }
      }
    }else if( starts_with( fields[0], "starttime" ) )
    {
      ++nlines_used;
      
      string timestr;
      if( nfields > 1 )
        timestr = fields[1];
      if( nfields > 2 )
        timestr += (" " + fields[2]);
      start_time_ = time_from_string( timestr.c_str() );
    }else if( starts_with( fields[0], "livetime" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
        live_time_ = time_duration_in_seconds( fields[1] );
    }else if( starts_with( fields[0], "realtime" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
        real_time_ = time_duration_in_seconds( fields[1] );
    }else if( starts_with( fields[0], "neutroncount" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
      {
        if( !(stringstream(fields[1]) >> neutron_counts_sum_) )
          throw runtime_error( "Invalid neutroncount: " + fields[1] );
        contained_neutron_ = true;
      }
    }else if( starts_with( fields[0], "samplenumber" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
      {
        if( !(stringstream(fields[1]) >> sample_number_) )
          throw runtime_error( "Invalid samplenumber: " + fields[1] );
      }
    }else if( starts_with( fields[0], "detectorname" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
        detector_name_ = fields[1];
    }else if( starts_with( fields[0], "detectortype" ) )
    {
      ++nlines_used;
      
      if( nfields > 1 )
        detector_type_ = fields[1];
    }else if( starts_with( fields[0], "title" ) )
    {
      ++nlines_used;
      
      string::size_type pos = line.find_first_of( " :\t" );
      if( pos < line.size() )
      {
        title_ = line.substr( pos + 1 );
        pos = title_.find_first_not_of( " :\t" );
        title_ = title_.substr( pos );
      }
    }else if( starts_with( fields[0], "calibcoeff" ) )
    {
      //CalibCoeff   : a=0.000000000E+000 b=0.000000000E+000 c=3.000000000E+000 d=0.000000000E+000
      ++nlines_used;
      
      float a = 0.0f, b = 0.0f, c = 0.0f, d = 0.0f;
      const size_t apos = line.find( "a=" );
      const size_t bpos = line.find( "b=" );
      const size_t cpos = line.find( "c=" );
      const size_t dpos = line.find( "d=" );
      if( apos < (line.size()-2) )
        a = static_cast<float>( atof( line.c_str() + apos + 2 ) );
      if( bpos < (line.size()-2) )
        b = static_cast<float>( atof( line.c_str() + bpos + 2 ) );
      if( cpos < (line.size()-2) )
        c = static_cast<float>( atof( line.c_str() + cpos + 2 ) );
      if( dpos < (line.size()-2) )
        d = static_cast<float>( atof( line.c_str() + dpos + 2 ) );
        
      if( c > 0 )
      {
        energy_calibration_model_ = Measurement::Polynomial;
        calibration_coeffs_.resize( 2 );
        calibration_coeffs_[0] = b;
        calibration_coeffs_[1] = c;
      }
    }
    
  }//while( getline( istr, line ) )

  if( nlines_total < 10 || nlines_used < static_cast<size_t>( ceil(0.25*nlines_total) ))
  {
    reset();
    istr.seekg( orig_pos, ios::beg );
    istr.clear( ios::failbit );
    throw runtime_error( "Not enough (useful) lines in the file." );
  }//
  
  if( !gamma_counts_ || gamma_counts_->size() < 5 || calibration_coeffs_.empty() )
  {
    reset();
    istr.seekg( orig_pos, ios::beg );
    istr.clear( ios::failbit );
    stringstream msg;
    msg << SRC_LOCATION << "\n\tI was unable to load the spectrum, probably"
        << " due to missing data or an invalid line somewhere";
    throw runtime_error( msg.str() );
  }//if( !gamma_counts_ || gamma_counts_->empty() )

  if( contained_neutron_ )
  {
    neutron_counts_.resize( 1 );
    neutron_counts_[0] = static_cast<float>( neutron_counts_sum_ );
  }//if( contained_neutron_ )

  for( const float f : *gamma_counts_ )
    gamma_count_sum_ += f;
}//void set_info_from_txt_or_csv( std::istream& istr )


bool MeasurementInfo::load_from_iaea( std::istream& istr )
{
  //channel data in $DATA:
  //live time, real time in $MEAS_TIM:
  //measurment datetime in $DATE_MEA:
  //Description in $SPEC_ID
  //Ploynomial calibration coefficients in $ENER_FIT: as well as $MCA_CAL:

  if( !istr.good() )
    return false;
  
  const istream::pos_type orig_pos = istr.tellg();
  
  reset();
  std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();

  //Each line should be terminated with "line feed" (0Dh)
  // and "carriage return" (0Ah), but for now I'm just using safe_get_line(...),
  // just to be safe.
  bool skip_getline = false;
  
  try
  {
    string line;
    while( UtilityFunctions::safe_get_line( istr, line ) )
    {
      trim(line);
      if( !line.empty() )
        break;
    }//while( UtilityFunctions::safe_get_line( istr, line ) )

    if( line.empty() || line[0]!='$' )
      throw runtime_error( "IAEA file first line must start with a '$'" );

    bool neutrons_were_cps = false;
    std::shared_ptr<DetectorAnalysis> anaresult;
    
    do
    {
      trim(line);
      to_upper(line);
      skip_getline = false;

      if( starts_with(line,"$DATA:") )
      {
        //RadEagle files contains a seemingly duplicate section of DATA: $TRANSFORMED_DATA:
        
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error( "Error reading DATA section of IAEA file" );

        trim(line);
        vector<string> channelstrs;
        split( channelstrs, line, " \t," );

        unsigned int firstchannel = 0, lastchannel = 0;
        if( channelstrs.size() == 2 )
        {
          try
          {
            firstchannel = atol( channelstrs[0].c_str() );
            lastchannel = atol( channelstrs[1].c_str() );
          }catch(...)
          {
          }
        }else
        {
          passMessage( "Error reading DATA section of IAEA file, "
                       "unexpected number of fields in first line.",
                       "MeasurementInfo::load_from_iaea()", 1 );
        }//if( firstlineparts.size() == 2 )

        double sum = 0.0;
        std::shared_ptr< vector<float> > channel_data( new vector<float>() );
        if( firstchannel < lastchannel )
          channel_data->reserve( lastchannel - firstchannel + 1 );

        //XXX - for some reason I think this next test condition is a little
        //      fragile...
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )

          vector<float> linevalues;
          const bool ok = UtilityFunctions::split_to_floats( line.c_str(), line.size(), linevalues );
        
          if( !ok )
          {
            char buffer[1024];
            snprintf( buffer, sizeof(buffer),
                      "Error converting channel data to counts for line: '%s'",
                      line.c_str() );
            cerr << buffer << endl;
#if(PERFORM_DEVELOPER_CHECKS)
            log_developer_error( BOOST_CURRENT_FUNCTION, buffer );
#endif
            continue;
          }//if( !ok )
          
//          sum += std::accumulate( linevalues.begin(), linevalues.end(), std::plus<float>() );
          for( size_t i = 0; i < linevalues.size(); ++i )
            sum += linevalues[i];
          channel_data->insert( channel_data->end(), linevalues.begin(), linevalues.end() );
        }//while( we havent reached the nex section )

        meas->gamma_counts_ = channel_data;
        meas->gamma_count_sum_ = sum;
      }else if( starts_with(line,"$MEAS_TIM:") )
      {
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error( "Error reading MEAS_TIM section of IAEA file" );
        vector<string> fields;
        split( fields, line, " \t," );
        if( fields.size() == 2 )
        {
          meas->live_time_ = static_cast<float>( atof( fields[0].c_str() ) );
          meas->real_time_ = static_cast<float>( atof( fields[1].c_str() ) );
        }else
        {
          passMessage( "Error reading MEAS_TIM section of IAEA file, "
                       "unexpected number of fields.",
                       "MeasurementInfo::load_from_iaea()", 1 );
        }//if( firstlineparts.size() == 2 )
      }else if( starts_with(line,"$DATE_MEA:") )
      {
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error( "Error reading MEAS_TIM section of IAEA file" );

        trim(line);

        try
        {
          //Nominally formated like: "mm/dd/yyyy hh:mm:ss" (ex "02/29/2016 14:31:47")
          //  which time_from_string should get right...
          meas->start_time_ = time_from_string( line.c_str() );
        }catch(...)
        {
          passMessage( "Unable to convert date/time '" + line +
                       "' to a valid posix time",
                       "MeasurementInfo::load_from_iaea()", 1 );
        }
      }else if( starts_with(line,"$SPEC_ID:") )
      {
        string remark;
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.empty() )
            continue;
          
          if( starts_with(line,"RE ") )
          {
            //Is a ORTEC RADEAGLE
            //See http://www.ortec-online.com/-/media/ametekortec/brochures/radeagle.pdf to decode model number/options
            //Example "RE 3SG-H-GPS"
            //const bool has_gps = UtilityFunctions::icontains(line,"-GPS");
            //const bool has_neutron = UtilityFunctions::icontains(line,"-H");
            const bool has_underwater = UtilityFunctions::icontains(line,"SGA");
            if( has_underwater )
              remarks_.push_back( "Detector has under water option" );
            
            manufacturer_ = "Ortec";
            instrument_model_ = "RadEagle " + line.substr(3);
            instrument_type_ = "RadionuclideIdentifier";
          }else if( starts_with(line,"SN#") )
          {
            instrument_id_ = line.substr(3);
            trim(instrument_id_);
          }else if( starts_with(line,"HW#") )
          { //ex. "HW# HW 2.1 SW 2.34"
            line = line.substr(3);
            string hw, sw;
            const size_t hw_pos = line.find("HW");
            if( hw_pos != string::npos )
            {
              hw = line.substr( hw_pos + 2 );
              const size_t pos = hw.find("SW");
              if( pos != string::npos )
                hw = hw.substr(0,pos);
            }//if( hw_pos != string::npos )
            
            const size_t sw_pos = line.find("SW");
            if( sw_pos != string::npos )
            {
              sw = line.substr( sw_pos + 2 );
              const size_t pos = sw.find("HW");
              if( pos != string::npos )
                sw = sw.substr(0,pos);
            }//if( sw_pos != string::npos )
            
            trim( hw );
            trim( sw );
            
            if( !hw.empty() )
              component_versions_.push_back( pair<string,string>("HardwareVersion",hw) );
            
            if( !sw.empty() )
              component_versions_.push_back( pair<string,string>("SoftwareVersion",sw) );
            
         // }else if( starts_with(line,"AW#") )  //Not sure what AW is
          //{
          }else
          {
            remark += (!remark.empty() ? " " : "") + line;
          }
        }//while( UtilityFunctions::safe_get_line( istr, line ) )

        remarks_.push_back( remark );
      }else if( starts_with(line,"$ENER_FIT:")
                || starts_with(line,"$GAIN_OFFSET_XIA:") )
      {
        if( !starts_with(line,"$GAIN_OFFSET_XIA:") || meas->calibration_coeffs_.empty() )
        {
          if( !UtilityFunctions::safe_get_line( istr, line ) )
            throw runtime_error( "Error reading ENER_FIT section of IAEA file" );

          trim(line);
          meas->energy_calibration_model_ = Measurement::Polynomial;

          vector<float> coefs;
          if( UtilityFunctions::split_to_floats( line.c_str(), line.size(), coefs ) )
          {
            bool allZeros = true;
            for( const float c : coefs )
              allZeros &= (fabs(c)<1.0E-08);
          
            if( !allZeros )
              meas->calibration_coeffs_ = coefs;
          }
        }else
        {
          UtilityFunctions::safe_get_line( istr, line );
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }
      }else if( starts_with(line,"$MCA_CAL:") )
      {
        try
        {
          if( !UtilityFunctions::safe_get_line( istr, line ) )
            throw runtime_error("Error reading MCA_CAL section of IAEA file");
          trim(line);
          
          const int npar = atoi( line.c_str() );
          
          if( line.empty() || npar < 1 )
            throw runtime_error("Invalid number of parameters");
          
          if( !UtilityFunctions::safe_get_line( istr, line ) )
            throw runtime_error("Error reading MCA_CAL section of IAEA file");
          trim(line);
          meas->energy_calibration_model_ = Measurement::Polynomial;

          vector<float> coefs;
          try
          {
            const bool success = UtilityFunctions::split_to_floats(
                                             line.c_str(), line.size(), coefs );

            if( !success )
              coefs.clear();
            
            //make sure the file didnt just have all zeros
            bool allZeros = true;
            for( const float c : coefs )
              allZeros &= (fabs(c)<1.0E-08);
            
            if( !allZeros && (coefs.size() == npar) )
            {
              meas->calibration_coeffs_ = coefs;
            }else
            {
              coefs.clear();
              if( allZeros )
              {
                passMessage( "Calibration in file is {0, 0, 0}; not using.",
                             "", 1 );
              }else
              {
                stringstream msg;
                msg << "Unexpected number of calibration parameters in "
                       "IAEA file, expected " << npar << " found "
                    << coefs.size();
                passMessage( msg.str(), "MeasurementInfo::load_from_iaea()", 1 );
//              throw runtime_error( "Unexpected number of energy calibration parameters" );
              }
            }
          }catch(...)
          {
          }
        }catch( exception &e )
        {
          cerr << SRC_LOCATION << "\n\terror in MCA_CAL section of IAEA file"
               << "\n\t" << e.what();
        }
      }else if( starts_with(line,"$GPS:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          string valuestr;
          const string::size_type pos = line.find( "=" );
          if( pos != string::npos )
            valuestr = line.substr( pos );
          trim( valuestr );
          
          if( starts_with( line, "Lon=") )
          {
            if( sscanf( valuestr.c_str(), "%lf", &(meas->longitude_) ) != 1 )
              meas->longitude_ = -999.9;
          }else if( starts_with( line, "Lat=") )
          {
            if( sscanf( valuestr.c_str(), "%lf", &(meas->latitude_) ) != 1 )
              meas->latitude_ = -999.9;
          }else if( starts_with( line, "Speed=") )
          {
            if( sscanf( valuestr.c_str(), "%f", &(meas->speed_) ) != 1 )
              meas->speed_ = 0.0;
          }else if( !line.empty() )
            remarks_.push_back( line ); //also can be Alt=, Dir=, Valid=
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$GPS_COORDINATES:") )
      {
        if( UtilityFunctions::safe_get_line( istr, line ) )
          remarks_.push_back( "GPS Coordinates: " + line );
      }else if( starts_with(line,"$NEUTRONS:") )
      { //ex "0.000000  (total)"
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error("Error reading NEUTRONS section of IAEA file");
        trim( line );
        float val;
        if( toFloat(line,val) )
        {
          meas->neutron_counts_.push_back( val );
          meas->neutron_counts_sum_ += val;
          meas->contained_neutron_ = true;
        }else
          passMessage( "Error parsing neutron counts from line: " + line,
                       "MeasurementInfo::load_from_iaea()", 1 );
      }else if( starts_with(line,"$NEUTRONS_LIVETIME:") )
      { //ex "267706.437500  (sec)"
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error("Error reading NEUTRONS_LIVETIME section of IAEA file");
        trim( line );
        meas->remarks_.push_back( "Neutron Live Time: " + line );
      }else if( starts_with(line,"$NEUTRON_CPS:") )
      { //found in RadEagle SPE files
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error("Error reading NEUTRON_CPS section of IAEA file");
        trim( line );
        float val;
        if( toFloat(line,val) )
        {
          neutrons_were_cps = true;
          meas->neutron_counts_.push_back( val );
          meas->neutron_counts_sum_ += val;
          meas->contained_neutron_ = true;
        }else
          passMessage( "Error parsing neutron cps from line: " + line,
                      "MeasurementInfo::load_from_iaea()", 1 );
      }else if( starts_with(line,"$SPEC_REM:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          if( !line.empty() )
            meas->remarks_.push_back( line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$ROI:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
//          if( !line.empty() )
//            meas->remarks_.push_back( line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$ROI_INFO:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          if( !line.empty() )
          {
            vector<float> parts;
            UtilityFunctions::split_to_floats( line.c_str(), line.size(), parts );
            if( parts.size() > 7 )
            {
              const int roinum     = static_cast<int>( parts[0] );
              const int startbin   = static_cast<int>( parts[1] );
              const int endbin     = static_cast<int>( parts[2] );
              const float meanbin  = parts[3];
              const float fwhmbins = parts[4];
              const int roiarea    = static_cast<int>( parts[5] );
              const int peakarea   = static_cast<int>( parts[6] );
              const int areauncert = static_cast<int>( parts[7] );
              
              char buffer[512];
              snprintf( buffer, sizeof(buffer),
                        "ROI in file: { \"roinum\": %i, \"startbin\": %i,"
                        " \"endbin\": %i, \"meanbin\": %.2f, \"fwhmbins\": %.2f,"
                        " \"roiarea\": %i, \"peakarea\": %i, "
                        "\"peakareauncert\": %i }",
                        roinum, startbin, endbin, meanbin, fwhmbins, roiarea,
                        peakarea, areauncert );
              
              meas->remarks_.push_back( buffer );
            }//if( parts.size() > 7 )
          }
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$ENER_DATA:") || starts_with(line,"$MCA_CAL_DATA:") )
      {
        if( UtilityFunctions::safe_get_line( istr, line ) )
        {
          const size_t num = static_cast<size_t>( atol( line.c_str() ) );
          vector<pair<int,float> > bintoenergy;
          while( UtilityFunctions::safe_get_line( istr, line ) )
          {
            trim(line);
            if( starts_with(line,"$") )
            {
              skip_getline = true;
              break;
            }//if( we have overrun the data section )
            
            if( !line.empty() )
            {
              vector<float> parts;
              UtilityFunctions::split_to_floats( line.c_str(), line.size(), parts );
              if( parts.size() == 2 )
                bintoenergy.push_back( make_pair(static_cast<int>(parts[0]),parts[1]) );
            }
          }//while( UtilityFunctions::safe_get_line( istr, line ) )
          
          if( num && bintoenergy.size()==num )
          {
            stringstream remarkstrm;
            remarkstrm << "Calibration in file from:";
            for( size_t i = 0; i < num; ++i )
              remarkstrm << (i?", ":" ") << "bin " << bintoenergy[i].first
                         << "->" << bintoenergy[i].second << " keV";
            meas->remarks_.push_back( remarkstrm.str() );
          }
        }//if( file says how manyy entries to expect )
      }else if( starts_with(line,"$SHAPE_CAL:") )
      {
        //I think this is FWHM calibration parameters - skipping this for now
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }//while( UtilityFunctions::safe_get_line( istr, line ) )

/*
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error("Error reading SHAPE_CA section of IAEA file");
        trim(line);
        unsigned int npar;
        if( !(stringstream(line) >> npar ) )
          throw runtime_error( "Invalid number of parameter: " + line );

        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error("Error reading SHAPE_CA section of IAEA file");
        trim(line);

        vector<float> coefs;
        vector<string> fields;
        split( fields, line, " \t," );
        try
        {
          for( const string &field : fields )
          {
            float val;
            if( !(stringstream(field) >> val) )
               throw runtime_error( "Invalid calibration parameter: " + field );
            coefs.push_back( val );
          }
 
          if( npar != coefs.size() )
            throw runtime_error( "Unexpected number of parameters in SHAPE_CA block of IAEA file." );
        }catch(...)
        {
        }
*/
      }else if( starts_with(line,"$PEAKLABELS:") )
      {
        //XXX - it would be nice to keep the peak lables, or at least search for
        //      these peaks and assign the isotopes..., but for now we'll ignore
        //      them :(
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$CAMBIO:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$APPLICATION_ID:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          //line something like "identiFINDER 2 NG, Application: 2.37, Operating System: 1.2.040"
          if( line.size() )
            remarks_.push_back( line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$DEVICE_ID:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( UtilityFunctions::icontains( line, "identiFINDER 2 LG" ) )
          {
            detector_type_ = kIdentiFinderLaBr3Detector;
            instrument_model_ = line;
            manufacturer_ = "FLIR";
            
            if( UtilityFunctions::icontains( line, "LGH" ) )
              meas->contained_neutron_ = true;
          }else if( UtilityFunctions::icontains( line, "identiFINDER 2 NG") ) //"nanoRaider ZH"
          {
            detector_type_ = kIdentiFinderNGDetector;
            instrument_model_ = line;
            manufacturer_ = "FLIR";
            
            if( UtilityFunctions::icontains( line, "NGH" ) )
              meas->contained_neutron_ = true;
          }else if( UtilityFunctions::istarts_with( line, "SN#") )
          {
            line = line.substr(3);
            UtilityFunctions::trim( line );
            if( line.size() )
              instrument_id_ = line;
          }else
            remarks_.push_back( line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_DATASET_NUMBER:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          //line something like "identiFINDER 2 NG, Application: 2.37, Operating System: 1.2.040"
          if( line.size() )
            remarks_.push_back( "FLIR DATSET NUMBER: " + line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_GAMMA_DETECTOR_INFORMATION:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.size() )
            remarks_.push_back( "GAMMA DETECTOR INFORMATION: " + line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_NEUTRON_DETECTOR_INFORMATION:") )
      {
        meas->contained_neutron_ = true;
        
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.size() )
            remarks_.push_back( "NEUTRON DETECTOR INFORMATION: " + line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_SPECTRUM_TYPE:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( UtilityFunctions::icontains( line, "IntrinsicActivity" ) )
            meas->source_type_ = Measurement::IntrinsicActivity;
          else if( UtilityFunctions::icontains( line, "Measurement" ) )
            meas->source_type_ = Measurement::Foreground;
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_REACHBACK:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.size() )
            remarks_.push_back( "Reachback url: " + line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_DOSE_RATE_SWMM:") )
      {
        //I dont understand the structure of this... so just putting in remarks.
        string remark;
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.size() )
            remark += (remark.size()?", ":"") + line;
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
        
        if( remark.size() )
          remarks_.push_back( "Dose information: " + remark );
      }else if( starts_with(line,"$FLIR_ANALYSIS_RESULTS:") )
      {
        vector<string> analines;
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          if( !line.empty() )
            analines.push_back( line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
        
        if( analines.empty() )
          continue;
        
        int numresults = 0;
        if( !toInt(analines[0], numresults) || numresults <= 0 )
          continue;
        
        if( (analines.size()-1) != numresults*4 )
        {
          string remark;
          for( size_t i = 0; i < analines.size(); ++i )
            remark += (remark.size()?", ":"") + analines[i];
          remarks_.push_back( "FLIR_ANALYSIS_RESULTS not in expected format: " + remark );
          continue;
        }
          
        if( !anaresult )
          anaresult = std::make_shared<DetectorAnalysis>();
        
        for( int ananum = 0; ananum < numresults; ++ananum )
        {
          DetectorAnalysisResult result;
          result.nuclide_       = analines.at(1+4*ananum); //"Th-232 or U-232"
          result.nuclide_type_  = analines.at(2+4*ananum); //"NORM"
          result.id_confidence_ = analines.at(4+4*ananum); //"5" - not really sure this is actually confidence...
          result.remark_        = analines.at(3+4*ananum); //"Innocent"
          
          anaresult->results_.push_back( result );
        }//for( int ananum = 0; ananum < numresults; ++ananum )
      }else if( starts_with(line,"$DOSE_RATE:") )
      {  //Dose rate in uSv.  Seen in RadEagle
        if( !UtilityFunctions::safe_get_line( istr, line ) )
          throw runtime_error( "Error reading DOSE_RATE section of IAEA file" );
        
        trim(line);
        skip_getline = starts_with(line,"$");
        
        DetectorAnalysisResult result;
        if( !toFloat(line, result.dose_rate_) )
        {
          remarks_.push_back( "Error reading DOSE_RATE, line: " + line );
        }else
        {
          if( !anaresult )
            anaresult = std::make_shared<DetectorAnalysis>();
          anaresult->results_.push_back( result );
        }
      }else if( starts_with(line,"$RADIONUCLIDES:") )
      { //Have only seen one file with this , and it only had a single nuclide
        //Cs137*[9.58755]
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          //Super unsure about the formatting so far (20171218)!
          const size_t nuc_end = line.find_first_of("*[");
          if( nuc_end != string::npos )
          {
            DetectorAnalysisResult result;
            result.nuclide_ = line.substr(0,nuc_end);
          
            //result.nuclide_type_  = "";
            const size_t conf_start = line.find('[');
            const size_t conf_end = line.find(']');
            if( conf_end > conf_start && conf_end != string::npos )
              result.id_confidence_ = line.substr(conf_start+1,conf_end-conf_start-1);
            result.remark_        = line;  //20171218: include everything for until I get better confirmation of format...
            
          
            if( !anaresult )
              anaresult = std::make_shared<DetectorAnalysis>();
            anaresult->results_.push_back( result );
          }//if( nuc_end != string::npos )
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$SPEC_INTEGRAL:") )
      {
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line.size() )
            remarks_.push_back( "SPEC_INTEGRAL: " + line );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$IDENTIFY_PARAMETER:") )
      {
        vector<float> calibcoefs;
        vector< pair<float,float> > fwhms;
        while( !skip_getline && UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          const bool slope = UtilityFunctions::icontains( line, "Energieeichung_Steigung");
          const bool offset = UtilityFunctions::icontains( line, "Energieeichung_Offset");
          const bool quad = UtilityFunctions::icontains( line, "Energieeichung_Quadrat");
          if( slope || offset || quad )
          {
            if( UtilityFunctions::safe_get_line( istr, line ) )
            {
              trim(line);
              if( starts_with(line,"$") )
              {
                skip_getline = true;
                break;
              }//if( we have overrun the data section )
              
              if( slope )
              {
                if( calibcoefs.size() < 2 )
                  calibcoefs.resize( 2, 0.0f );
                if( !toFloat( line, calibcoefs[1] ) )
                  throw runtime_error( "Couldnt convert to cal slope to flt: " + line );
              }else if( offset )
              {
                if( calibcoefs.size() < 1 )
                  calibcoefs.resize( 1, 0.0f );
                if( !toFloat( line, calibcoefs[0] ) )
                  throw runtime_error( "Couldnt convert to cal offset to flt: " + line );
              }else if( quad )
              {
                if( calibcoefs.size() < 3 )
                  calibcoefs.resize( 3, 0.0f );
                if( !toFloat( line, calibcoefs[2] ) )
                  throw runtime_error( "Couldnt convert to cal quad to flt: " + line );
              }
            }//if( UtilityFunctions::safe_get_line( istr, line ) )
          }//if( UtilityFunctions::icontains( line, "Energieeichung_Steigung") )
          
/*
            FWHM_FWHM1:
            18.5
            FWHM_Energie1:
            122
            FWHM_FWHM2:
            43.05
            FWHM_Energie2:
            662
            Detektortyp:
            2
 */
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
        
        if( calibcoefs.size() && !meas->calibration_coeffs_.size() )
        {
          meas->energy_calibration_model_ = Measurement::Polynomial;
          meas->calibration_coeffs_ = calibcoefs;
        }
      }else if( starts_with(line,"$NON_LINEAR_DEVIATIONS:") )
      {
        UtilityFunctions::safe_get_line( istr, line );
        trim(line);
        
        const size_t npairs = static_cast<size_t>( atoi(line.c_str()) );
        
        if( line.empty() || npairs < 1 )
          continue;
        
        DeviationPairVec pairs;
        
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          float energy, deviation;
          if( !(stringstream(line) >> energy >> deviation ) )
          {
            pairs.clear();
            break;
          }//if( !(stringstream(line) >> energy >> deviation ) )

          pairs.push_back( pair<float,float>(energy,deviation) );
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
        
        if( pairs.size() == npairs )
        {
          meas->deviation_pairs_ = pairs;
        }else if( npairs || pairs.size() )
        {
          char buffer[256];
          snprintf( buffer, sizeof(buffer),
                    "Error parsing deviation pairs, expected %i, read in %i; "
                    "not using", int(npairs), int(pairs.size()) );
          passMessage( buffer, "MeasurementInfo::load_from_iaea()", 1 );
        }//if( pairs.size() == npairs )
      }else if( starts_with(line,"$ENDRECORD:") )
      {
        measurements_.push_back( meas );
        meas.reset( new Measurement() );
      }else if( starts_with(line,"$RT:")
               || starts_with(line,"$DT:") )
      {
        //Burn off things we dont care about
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }
      }else if( starts_with(line,"$IDENTIFY_NUKLIDE:")
               || starts_with(line,"$IDENTIFY_PEAKS:")
               || starts_with(line,"$PRESETS:")
               || starts_with(line,"$ICD_TYPE:")
               || starts_with(line,"$TEMPERATURE:")
              )
      {
        //Burn off things we dont care about
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
          
          if( line == "Live Time" )
          {
            if( UtilityFunctions::safe_get_line( istr, line ) )
            {
              if( starts_with(line,"$") )
              {
                skip_getline = true;
                break;
              }//if( we have overrun the data section )
              remarks_.push_back( "A preset live time of " + line + " was used" );
            }//if( UtilityFunctions::safe_get_line( istr, line ) )
          }//if( line == "Live Time" )
//          if( line.size() )
//            cout << "Got Something" << endl;
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( starts_with(line,"$FLIR_NEUTRON_SWMM:")
               || starts_with(line,"$TRANSFORMED_DATA:") )
      {
        //we'll just burn through this section of the file since wcjohns doesnt
        //  understand the purpose or structure of these sections
        while( UtilityFunctions::safe_get_line( istr, line ) )
        {
          trim(line);
          if( starts_with(line,"$") )
          {
            skip_getline = true;
            break;
          }//if( we have overrun the data section )
        }//while( UtilityFunctions::safe_get_line( istr, line ) )
      }else if( !line.empty() && line != "END" )
      {
        cerr << "Unrecognized line '" << line << "'" << endl;
      }

  /*
   *GADRAS users manual recomends having $SAMPLE_DESCRIPTION:, $DATA_TYPE:, $DETECTOR_TYPE:, $RADIATION_TYPE:, $METHOD_TYPE:, $MCA_CAL_DATA:
   *IAEA pdfs say the following are also possible.
      $ADC:
      $ADD_HV:
      $BT:
      $COUNTS:
      $DT:
      $DTC:
      $ENER_DATA:
      $ENER_DATA_X:
      $ENER_FIT:
      $FAST_DISCR:
      $GAIN_VALUE:
      $HV:
      $INPUT:
      $INSP_INFO:
      $MCA_166_ID:
      $MCA_REPEAT:
      $MEAS_TIM:
      $MCS_ADD_DATA:
      $MCS_AMP_DATA:
      $MCS_AMP_ROI:
      $MCS_AMP_ROI_INFO:
      $MCS_CHANNELS:
      $MCS_INPUT:
      $MCS_SWEEPS:
      $MCS_TIME:
      $MODE:
      $PD_COUNTS:
      $POWER:
      $POWER_STATE:
      $PRESETS:
      $PUR:
      $PZC_VALUE:
      $REC_COUNTER:
      $REC_ERROR_COUNTER:
      $ROI:
      $ROI_INFO:
      $RT:
      $SCANDU:
      $SCANDU_RESULTS:
      $Sensors:
      $SINGLE_POINTS:
      $SLOW_DISCR:
      $SPEC_INTEGRAL:
      $SPEC_REM:
      $STAB:
      $STAB_COUNTER:
      $STAB_OFFSET:
      $STAB_OFFSET_MAX:
      $STAB_OFFSET_MIN:
      $STAB_PARAM:
      $TEMPERATURE:
      $TDF:
      $THR:
      $UF6_INSP_INFO:
      $WATCHMON_ROI_INFO:
      $WINSCAN_INFO:
      $WINSPEC_AUTOMATION:
      $WINSPEC_INFO:
      */
    }while( skip_getline || UtilityFunctions::safe_get_line( istr, line) );
    //  while( getline( istr, line ) )

    if( neutrons_were_cps )
    {
      if( meas->real_time_ > 0.0f )
      {
        meas->neutron_counts_sum_ *= meas->real_time_;
        for( size_t i = 0; i < meas->neutron_counts_.size(); ++i )
          meas->neutron_counts_[i] *= meas->real_time_;
      }else
      {
        meas->remarks_.push_back( "Neutron counts is in counts per second (real time was zero, so could not determine gross counts)" );
      }
    }//if( neutrons_were_cps )
    
    if( anaresult )
      detectors_analysis_ = anaresult;
    
    if( (measurements_.empty() || measurements_.back()!=meas)
        && meas && meas->gamma_counts_ && meas->gamma_counts_->size() )
      measurements_.push_back( meas );
  }catch(...)
  {
    istr.clear();
    istr.seekg( orig_pos, ios::beg );
    reset();
    return false;
  }

  cleanup_after_load();

  return true;
}//bool load_from_iaea( std::istream& ostr )


bool MeasurementInfo::load_from_Gr135_txt( std::istream &input )
{
  //See data file for refIED7MP8PYR for an example file
  reset();

  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  
  try
  {
    string line;
    if( !UtilityFunctions::safe_get_line( input, line ) )
      return false;
  
    vector<string> headers;
    split( headers, line, "\t" );

    if( headers.empty() )
      throw runtime_error( "Not the GR135 header expected" );
  
    //each header will look like:
    // "1899715091 Oct. 09 2013  12:29:38 T counts Live time (s) 279.4 neutron 1 gieger 194"
    std::vector< MeasurementShrdPtr > measurments;
    std::vector< ShrdFVecPtr > gammacounts;
    
    for( size_t i = 0; i < headers.size(); ++i )
    {
      const string &header = headers[i];
      if( header.empty() )
        continue;

      MeasurementShrdPtr meas( new Measurement() );
      
      string::size_type pos = header.find( ' ' );
      if( pos == string::npos )
        throw runtime_error( "Invalid GR135 measurment header" );
      
      const string measIDStr = header.substr( 0, pos );
      string timestampStr = header.substr( pos+1 );
      pos = timestampStr.find( " T " );
      if( pos == string::npos )
        pos = timestampStr.find( " R " );
      if( pos == string::npos )
        pos = timestampStr.find( " Q " );
      if( pos == string::npos )
        pos = timestampStr.find( " counts " );
      
      if( pos == string::npos )
        throw runtime_error( "Couldnt find end of GR135 timestamp string" );
      timestampStr = timestampStr.substr( 0, pos );  //"Oct. 09 2013  13:08:29"
      
      meas->start_time_ = time_from_string( timestampStr.c_str() );
#if(PERFORM_DEVELOPER_CHECKS)
      if( meas->start_time_.is_special() )
        log_developer_error( BOOST_CURRENT_FUNCTION, ("Failed to extract measurment start time from: '" + header  + "'").c_str() );
#endif
      
      pos = header.find( "Live time (s)" );
      if( pos == string::npos )
        throw runtime_error( "Couldnt find Live time" );
    
      string liveTimeStr = header.substr( pos + 13 );
      trim( liveTimeStr );
      pos = liveTimeStr.find( " " );
      if( pos != string::npos )
        liveTimeStr = liveTimeStr.substr( 0, pos );
      
      if( !liveTimeStr.empty() )
      {
        float val;
        if( !toFloat( liveTimeStr, val ) )
          throw runtime_error( "Error converting live time to float" );
        meas->live_time_ = static_cast<float>( val );
      }//if( !liveTimeStr.empty() )

      pos = header.find( "neutron" );
      if( pos != string::npos )
      {
        string neutronStr = header.substr( pos + 7 );
        trim( neutronStr );
        
        if( !neutronStr.empty() )
        {
          float val;
          if( !toFloat( neutronStr, val ) )
            throw runtime_error( "Error converting neutron counts to float" );
          meas->neutron_counts_.resize( 1, val );
          meas->neutron_counts_sum_ = val;
          meas->contained_neutron_ = true;
        }//if( !liveTimeStr.empty() )
//        cerr << "meas->neutron_counts_sum_=" << meas->neutron_counts_sum_ << endl;
      }//if( pos != string::npos )
      
      pos = header.find( "gieger" );
      if( pos != string::npos )
      {
        string remark = header.substr( pos );
        trim( remark );
        meas->remarks_.push_back( remark );
      }
      
      ShrdFVecPtr channelcounts( new vector<float>() );
      channelcounts->reserve( 1024 );
      meas->gamma_counts_ = channelcounts;
      meas->gamma_count_sum_ = 0.0;
      
      meas->sample_number_ = static_cast<int>(i) + 1;
      
      //Some default calibration coefficients that are kinda sorta close
//      meas->calibration_type_ = Measurement::GuesedCalibration;
      meas->energy_calibration_model_ = Measurement::Polynomial;
      meas->calibration_coeffs_.push_back( 0 );
      meas->calibration_coeffs_.push_back( 3.0 );
      
      measurments.push_back( meas );
      gammacounts.push_back( channelcounts );
      measurements_.push_back( meas );
    }//for( size_t i = 0; i < headers.size(); ++i )
  
    if( gammacounts.empty() )
      throw runtime_error( "No GR135 txt file header" );
    
    std::vector<float> counts;
    
    while( UtilityFunctions::safe_get_line( input, line ) )
    {
      if( line.empty() )
        continue;
      
      UtilityFunctions::split_to_floats( line.c_str(), line.size(), counts );
      if( counts.size() != measurments.size() )
        throw runtime_error( "Unexpected number of channel counts" );
    
      for( size_t i = 0; i < counts.size(); ++i )
      {
        gammacounts.at(i)->push_back( counts[i] );
        measurments[i]->gamma_count_sum_ += counts[i];
      }
    }//while( UtilityFunctions::safe_get_line( input, line ) )
    
    const uint32_t len = static_cast<uint32_t>( gammacounts[0]->size() );
    const bool isPowerOfTwo = ((len != 0) && !(len & (len - 1)));
    if( !isPowerOfTwo )
      throw runtime_error( "Invalid number of channels" );
    
    cleanup_after_load();
    
    return true;
  }catch( std::exception &e )
  {
    cerr << "load_from_Gr135_txt(...) caught: " << e.what() << endl;
  }//try / catch
  
  input.clear();
  input.seekg( orig_pos, ios::beg );
  
  reset();

  return false;
}//bool load_from_Gr135_txt( std::istream &istr )


bool MeasurementInfo::load_from_txt_or_csv( std::istream &istr )
{
  reset();

  if( !istr.good() )
    return false;
  
  const std::streampos startpos = istr.tellg();
  
  string firstdata;
  firstdata.resize( 20, '\0' );
  if( !istr.read(&(firstdata[0]), 19) )
    return false;

  //Non- exaustive list of formats that we might be able to extract a spectrum
  //  from, but we really shouldnt, because its N42
  const char *not_allowed_txt[] = { "<?xml", "<Event", "<N42InstrumentData" };
  for( const char *txt : not_allowed_txt )
  {
    if( UtilityFunctions::icontains(firstdata, txt) )
      return false;
  }
  
  istr.seekg( startpos, ios::beg );
  
  while( istr.good() )
  {
    try
    {
      MeasurementShrdPtr m = std::make_shared<Measurement>();
      m->set_info_from_txt_or_csv( istr );
      
      measurements_.push_back( m );
    }catch( exception & )
    {
//      cerr << SRC_LOCATION << "\n\tCaught: " << e.what() << endl;
      break;
    }
  }//while( istr.good() )

  
  if( measurements_.empty() )
  {
    reset();
    istr.clear();
    istr.seekg( startpos, ios::end );
    return false;
  }
  
  try
  {
    cleanup_after_load();
  }catch( std::exception &e )
  {
    cerr << SRC_LOCATION << "\n\tCaught: " << e.what() << endl;
    reset();
    istr.clear();
    istr.seekg( startpos, ios::end );
    return false;
  }//try / catch

  if( measurements_.empty() )
  {
    istr.clear();
    istr.seekg( startpos, ios::end );
    reset();
    return false;
  }
  
  return true;
}//bool load_from_txt_or_csv( std::ostream& ostr )


bool MeasurementInfo::load_spectroscopic_daily_file( const std::string &filename )
{
  ifstream input( filename.c_str(), ios_base::binary|ios_base::in );
  
  if( !input.is_open() )
    return false;
  
  char buffer[8];
  input.get( buffer, sizeof(buffer)-1 );
  buffer[sizeof(buffer)-1] = '\0'; //JIC
  
  string bufferstr = buffer;
  const bool isSDF = ((bufferstr.size() > 3 && bufferstr[2]==',')
                      && ( UtilityFunctions::starts_with( bufferstr, "GB" )
                          || UtilityFunctions::starts_with( bufferstr, "NB" )
                          || UtilityFunctions::starts_with( bufferstr, "S1" )
                          || UtilityFunctions::starts_with( bufferstr, "S2" )
                          || UtilityFunctions::starts_with( bufferstr, "GS" )
                          || UtilityFunctions::starts_with( bufferstr, "GS" )
                          || UtilityFunctions::starts_with( bufferstr, "NS" )
                          || UtilityFunctions::starts_with( bufferstr, "ID" )
                          || UtilityFunctions::starts_with( bufferstr, "AB" )));
  if( !isSDF )
    return false;
  
  input.seekg( 0, ios_base::beg );
  
  const bool success = load_from_spectroscopic_daily_file( input );
  
  if( !success )
    return false;
  
  filename_ = filename;
  //      Field 4, the equipment specifier, is as follows:
  //        -SPM-T for a Thermo ASP-C
  //        -SPM-C for a Canberra ASP-C
  //        -RDSC1 for the Radiation Detector Straddle Carrier in primary
  //        -RDSC2 for the Radiation Detector Straddle Carrier in secondary
  //        -MRDIS2 for the Mobile Radiation Detection and Identification System in secondary
  //        ex. refG8JBF6M229
  vector<string> fields;
  UtilityFunctions::split( fields, filename, "_" );
  if( fields.size() > 3 )
  {
    if( fields[3] == "SPM-T" )
    {
      manufacturer_ = "Thermo";
      instrument_model_ = "ASP";
    }else if( fields[3] == "SPM-C" )
    {
      manufacturer_ = "Canberra";
      instrument_model_ = "ASP";
    }else if( fields[3] == "RDSC1" )
    {
      inspection_ = "Primary";
      instrument_model_ = "Radiation Detector Straddle Carrier";
    }else if( fields[3] == "RDSC2" )
    {
      inspection_ = "Secondary";
      instrument_model_ = "Radiation Detector Straddle Carrier";
    }else if( fields[3] == "MRDIS2" )
    {
      inspection_ = "Secondary";
      instrument_model_ = "Mobile Radiation Detection and Identification System";
    }
  }//if( fields.size() > 3 )
      
  return true;
}//bool load_spectroscopic_daily_file( const std::string &filename )


bool MeasurementInfo::load_amptek_file( const std::string &filename )
{
  ifstream input( filename.c_str(), ios_base::binary|ios_base::in );
  
  if( !input.is_open() )
    return false;
  
  const bool success = load_from_amptek_mca( input );
  
  if( success )
    filename_ = filename;
  
  return success;
}//bool load_amptek_file( const std::string &filename )


bool MeasurementInfo::load_ortec_listmode_file( const std::string &filename )
{
  ifstream input( filename.c_str(), ios_base::binary|ios_base::in );
  
  if( !input.is_open() )
    return false;
  
  const bool success = load_from_ortec_listmode( input );
  
  if( success )
    filename_ = filename;
  
  return success;
}


bool MeasurementInfo::load_txt_or_csv_file( const std::string &filename )
{
  try
  {
    std::unique_ptr<ifstream> input( new ifstream( filename.c_str(), ios_base::binary|ios_base::in ) );

    if( !input->is_open() )
      return false;

    //lets make sure this is an ascii file
    //Really we should make sure its a UTF-8 file, so we can export from Excell
    while( input->good() )
    {
      const int c = input->get();
      if( input->good() && c>127 )
        return false;
    }//while( input.good() )

    //we have an ascii file if we've made it here
    input->clear();
    input->seekg( 0, ios_base::beg );

    
    //Check to see if this is a GR135 text file
    string firstline;
    UtilityFunctions::safe_get_line( *input, firstline );
    
    bool success = false;
    
    const bool isGR135File = contains( firstline, "counts Live time (s)" )
                             && contains( firstline, "gieger" );
    
    if( isGR135File )
    {
      input->seekg( 0, ios_base::beg );
      success = load_from_Gr135_txt( *input );
    }
    
    const bool isSDF = ((!success && firstline.size() > 3 && firstline[2]==',')
                         && ( UtilityFunctions::starts_with( firstline, "GB" )
                         || UtilityFunctions::starts_with( firstline, "NB" )
                         || UtilityFunctions::starts_with( firstline, "S1" )
                         || UtilityFunctions::starts_with( firstline, "S2" )
                         || UtilityFunctions::starts_with( firstline, "GS" )
                         || UtilityFunctions::starts_with( firstline, "GS" )
                         || UtilityFunctions::starts_with( firstline, "NS" )
                         || UtilityFunctions::starts_with( firstline, "ID" )
                         || UtilityFunctions::starts_with( firstline, "AB" )));
    if( isSDF )
    {
      input->close();
      input.reset();
      
      success = load_spectroscopic_daily_file( filename );
      
      if( success )
        return true;
      
      input.reset( new ifstream( filename.c_str(), ios_base::binary|ios_base::in ) );
    }//if( isSDF )
    
    if( !success )
    {
      input->clear();
      input->seekg( 0, ios_base::beg );
      success = load_from_txt_or_csv( *input );
    }
    
    if( success )
      filename_ = filename;
    else
      reset();

    return success;
  }catch(...)
  {}//try / catch

  reset();

  return false;
}//bool load_txt_or_csv_file( const std::string &filename )


namespace SpectroscopicDailyFile
{
  struct DailyFileS1Info
  {
    bool success;
    std::string detTypeStr;
    std::string appTypeStr;
    int nchannels;
    std::vector<float> calibcoefs;
    std::string algorithmVersion;
    //caputures max 1 attribute...
    struct params{ std::string name, value, attname, attval; };
    vector<DailyFileS1Info::params> parameters;
  };//struct DailyFileS1Info
  
  
  struct DailyFileEndRecord
  {
    bool success;
    std::string alarmColor; //Red, Yellow, Gree
    int occupancyNumber;
    boost::posix_time::ptime lastStartTime;
    std::string icd1FileName;
    float entrySpeed, exitSpeed;
  };//struct DailyFileEndRecord

  
  struct DailyFileAnalyzedBackground
  {
    enum BackgroundType{ Gamma, Neutrons };
    
    bool success;
    BackgroundType type;
    float realTime;
    std::shared_ptr< std::vector<float> > spectrum;
  };//struct DailyFileAnalyzedBackground
  
  struct DailyFileNeutronSignal
  {
    bool success;
    int numTimeSlicesAgregated;
    int timeChunkNumber;  //identical to one used for gamma
    vector<float> counts;  //Aa1, Aa2, Aa3, Aa4, Ba1, Ba2, Ba3, Ba4, Ca1, Ca2, Ca3, Ca4, Da1, Da2, Da3, Da4
  };//enum DailyFileNeutronSignal
  
  struct DailyFileGammaSignal
  {
    bool success;
    std::string detectorName;
    int timeChunkNumber;
    std::shared_ptr< std::vector<float> > spectrum;
  };//struct DailyFileGammaSignal
  
  
  struct DailyFileGammaBackground
  {
    bool success;
    std::string detectorName;
    std::shared_ptr< std::vector<float> > spectrum;
  };//struct DailyFileGammaBackground
  
  struct DailyFileNeutronBackground
  {
    bool success;
    float realTime;
    vector<float> counts;
  };//struct DailyFileNeutronBackground
  
  
  void parse_s1_info( const char * const data, const size_t datalen, DailyFileS1Info &info )
  {
    const string s1str( data, datalen );
    
    info.calibcoefs.clear();
    
    vector<string> s1fields;
    UtilityFunctions::split( s1fields, s1str, "," );
    if( s1fields.size() < 5 )
    {
      cerr << "parse_s1_info(): Invalid S1 line" << endl;
      info.success = false;
      return;
    }
    
    info.detTypeStr = s1fields[1];  //NaI or HPGe
    info.appTypeStr = s1fields[2];  //SPM, RDSC, MRDIS
    info.nchannels = atoi( s1fields[3].c_str() );  //typically 512 or 4096
    info.algorithmVersion = s1fields[4];
    
    if( info.nchannels <= 0 )
    {
      cerr << "parse_s1_info(): Invalid claimed number of channels" << endl;
      info.nchannels = 512;
    }
    
    for( size_t i = 5; i < (s1fields.size()-1); i += 2 )
    {
      DailyFileS1Info::params p;
      
      p.name = s1fields[i];
      p.value = s1fields[i+1];
      
      const string::size_type spacepos = p.name.find(' ');
      if( spacepos != string::npos )
      {
        const string::size_type equalpos = p.name.find('=',spacepos);
        if( equalpos != string::npos )
        {
          p.attval = p.name.substr( equalpos+1 );
          p.attname = p.name.substr( spacepos + 1, equalpos - spacepos - 1 );
          p.name = p.name.substr( 0, spacepos );
        }//if( equalpos != string::npos )
      }//if( spacepos != string::npos )
    }//for( size_t i = 5; i < (s1fields.size()-1); i += 2 )

    if( info.calibcoefs.empty() )
    {
      cerr << "parse_s1_info(): warning, couldnt find calibration coeffecicents"
           << endl;
      info.calibcoefs.push_back( 0.0f );
      info.calibcoefs.push_back( 3000.0f / info.nchannels );
    }//if( info.calibcoefs.empty() )

    info.success = true;
  }//bool parse_s1_info( const std::string &s1str, DailyFileS1Info &info )
  
  void parse_s2_info( const char * const data, const size_t datalen,
                      map<string,vector< pair<float,float> > > &answer )
  {
    const std::string s2str( data, datalen );
    answer.clear();
  
    vector<string> s2fields;
    UtilityFunctions::split( s2fields, s2str, "," );
    string detname;
    for( size_t i = 1; i < (s2fields.size()-1); )
    {
      const string &field = s2fields[i];
      const string &nextfield = s2fields[i+1];
      
      if( field.empty() || nextfield.empty() )
      {
        i += 2;
        continue;
      }
      
      if( isdigit(field[0]) )
      {
        const float energy = static_cast<float>( atof( field.c_str() ) );
        const float offset = static_cast<float>( atof( nextfield.c_str() ) );
        answer[detname].push_back( pair<float,float>(energy,offset) );
        i += 2;
      }else
      {
        detname = s2fields[i];
        ++i;
      }
    }//for( size_t i = 1; i < (s2fields.size()-1); )
  }//void parse_s2_info(...)
  

  void parse_end_record( const char * const data, const size_t datalen, DailyFileEndRecord &info )
  {
    const std::string str( data, datalen );
    
    vector<string> fields;
    UtilityFunctions::split( fields, str, "," );
    
    if( fields.size() < 5 )
    {
      info.success = false;
      return;
    }
    
    info.alarmColor = fields[1];
    info.occupancyNumber = atoi( fields[2].c_str() );
    info.lastStartTime = time_from_string( fields[3].c_str() );
    
//    cout << "'" << fields[3] << "'--->" << info.lastStartTime << endl;
    info.icd1FileName = fields[4];
    info.entrySpeed = (fields.size()>5) ? static_cast<float>(atof(fields[5].c_str())) : 0.0f;
    info.exitSpeed = (fields.size()>6) ? static_cast<float>(atof(fields[6].c_str())) : info.entrySpeed;
    
    info.success = true;
  }//bool parse_end_record(...)
  
  void parse_analyzed_background( const char * const data, const size_t datalen,
                                 DailyFileAnalyzedBackground &info )
  {
    const string line( data, datalen );
    std::string::size_type pos1 = line.find( ',' );
    assert( line.substr(0,pos1) == "AB" );
    std::string::size_type pos2 = line.find( ',', pos1+1 );
    if( pos2 == string::npos )
    {
      cerr << "parse_analyzed_background: unexpected EOL 1" << endl;
      info.success = false;
      return;
    }
    
    string type = line.substr( pos1+1, pos2-pos1-1 );
    if( UtilityFunctions::iequals( type, "Gamma" ) )
      info.type = DailyFileAnalyzedBackground::Gamma;
    else if( UtilityFunctions::iequals( type, "Neutron" ) )
      info.type = DailyFileAnalyzedBackground::Neutrons;
    else
    {
      cerr << "parse_analyzed_background: invalid type '" << type << "'" << endl;
      info.success = false;
      return;
    }
    
    pos1 = pos2;
    info.realTime = static_cast<float>( atof( line.c_str() + pos1 + 1 ) );
    pos1 = line.find( ',', pos1+1 );
    
    if( pos1 == string::npos )
    {
      cerr << "parse_analyzed_background: unexpected EOL 2" << endl;
      info.success = false;
      return;
    }
    
    info.spectrum.reset( new vector<float>() );
    
    if( info.type == DailyFileAnalyzedBackground::Neutrons )
    {
      const float nneut = static_cast<float>( atof( line.c_str() + pos1 + 1 ) );
      info.spectrum->resize( 1, nneut );
    }else
    {
      const char *start = line.c_str() + pos1 + 1;
      const size_t len = line.size() - pos1 - 2;
      const bool success
              = UtilityFunctions::split_to_floats( start, len, *info.spectrum );
      
      if( !success )
      {
        cerr << "parse_analyzed_background: did not decode spectrum" << endl;
        info.success = false;
        return;
      }
    }//if( neutron ) / ( gamma )
    
    info.success = true;
  }//void parse_analyzed_background(...)
  
  void parse_neutron_signal( const char * const data, const size_t datalen,
                                 DailyFileNeutronSignal &info )
  {
    const string line( data, datalen );
    std::string::size_type pos = line.find( ',' );
    if( pos == string::npos )
    {
      info.success = false;
      return;
    }
    
    const char *start = line.c_str() + pos + 1;
    const size_t len = line.size() - pos - 1;
    
    vector<float> vals;
    const bool success = UtilityFunctions::split_to_floats( start, len, vals );
    if( !success || vals.size() < 2 )
    {
      cerr << "parse_neutron_signal: did not decode spectrum" << endl;
      info.success = false;
      return;
    }
    
    info.numTimeSlicesAgregated = static_cast<int>( vals[0] );
    info.timeChunkNumber = static_cast<int>( vals[1] );
    info.counts.clear();
    info.counts.insert( info.counts.end(), vals.begin()+2, vals.end() );
    
    info.success = true;
  }//bool parse_analyzed_background
  
  void parse_gamma_signal( const char * const data, const size_t datalen,
                           DailyFileGammaSignal &info )
  {
    const string line( data, datalen );
    std::string::size_type pos1 = line.find( ',' );
    if( pos1 == string::npos )
    {
      info.success = false;
      return;
    }
    
    std::string::size_type pos2 = line.find( ',', pos1+1 );
    if( pos2 == string::npos )
    {
      info.success = false;
      return;
    }
    
    info.detectorName = line.substr( pos1+1, pos2-pos1-1 );
    pos1 = pos2;
    pos2 = line.find( ',', pos1+1 );
    if( pos2 == string::npos )
    {
      info.success = false;
      return;
    }
    
    info.timeChunkNumber = static_cast<int>( atoi( line.c_str() + pos1 + 1 ) );
    info.spectrum.reset( new vector<float>() );
    
    vector<float> vals;
    const char *start = line.c_str() + pos2 + 1;
    const size_t len = line.size() - pos2 - 1;
    const bool success = UtilityFunctions::split_to_floats( start, len, *info.spectrum );
    if( !success || info.spectrum->size() < 2 )
    {
      cerr << "parse_gamma_signal: did not decode spectrum" << endl;
      info.success = false;
      return;
    }
    
    info.success = true;
    return;
  }//void parse_gamma_signal()

  
  void parse_gamma_background( const char * const data, const size_t datalen,
                               DailyFileGammaBackground &info )
  {
    const string line( data, datalen );
    std::string::size_type pos1 = line.find( ',' );
    if( pos1 == string::npos )
    {
      info.success = false;
      return;
    }
    
    std::string::size_type pos2 = line.find( ',', pos1+1 );
    if( pos2 == string::npos )
    {
      info.success = false;
      return;
    }
    
    info.detectorName = line.substr( pos1+1, pos2-pos1-1 );
    info.spectrum.reset( new vector<float>() );
    
    const char *start = line.c_str() + pos2 + 1;
    const size_t len = line.size() - pos2 - 1;
    const bool success = UtilityFunctions::split_to_floats( start, len, *info.spectrum );
    if( !success || info.spectrum->size() < 2 )
    {
      cerr << "parse_gamma_background: did not decode spectrum" << endl;
      info.success = false;
      return;
    }
    
    info.success = true;
  }//bool parse_gamma_background(...)
    
  
  void parse_neutron_background( const char * const data, const size_t datalen,
                                 DailyFileNeutronBackground &info )
  {
    const string line( data, datalen );
    std::string::size_type pos1 = line.find( ',' );
    if( pos1 == string::npos )
    {
      info.success = false;
      return;
    }
    
    std::string::size_type pos2 = line.find( ',', pos1+1 );
    if( pos2 == string::npos )
    {
      info.success = false;
      return;
    }
    
    info.realTime = static_cast<float>( atof( line.c_str() + pos1 + 1 ) );
    
    const char *start = line.c_str() + pos2 + 1;
    const size_t len = line.size() - pos2 - 1;
    
    //Files like fail: ref1E5GQ2SW76
    const bool success = UtilityFunctions::split_to_floats( start, len, info.counts );
    if( !success || info.counts.size() < 2 )
    {
      cerr << "parse_neutron_background: did not decode counts" << endl;
      info.success = false;
      return;
    }

    info.success = true;
    return;
  }//bool parse_neutron_background(...)
  
}//namespace SpectroscopicDailyFile


bool MeasurementInfo::load_from_spectroscopic_daily_file( std::istream &input )
{
/* The daily file is a comma separated value file, with a carriage return and
   line feed denoting the end of each line.  The file is saved as a text (.txt) 
   file.  Spaces are not necessary after each comma, in an effort to minimize 
   the overall size of the file.
   In the future, all data provided in the Daily File will be energy calibrated, 
   according to the calibration parameters provided in the S1 and/or S2 lines.
   This is to ensure that calibration is being done correctly and consistently 
   between various institutions, laboratories, and individuals.  The calibration 
   parameters will be provided in case an individual wants to “unwrap” the data 
   back to the source information provided in the ICD-1 file.  This is not done 
   for GS line data as of Revision 3 of this data.
*/
  
  //This is a rough hack in; it would be nice to mmap things and read in that
  //  way, as well as handle potential errors better
  
#define DO_SDF_MULTITHREAD 0
  //Intitial GetLineSafe() for "dailyFile6.txt"
  //  To parse SDF took:  1.027030s wall, 1.000000s user + 0.020000s system = 1.020000s CPU (99.3%)
  //  Total File openeing time was:  1.164266s wall, 1.300000s user + 0.110000s system = 1.410000s CPU (121.1%)
  //
  //Adding an extra indirection of creating a copy of a string
  //  To parse SDF took:  1.162067s wall, 1.120000s user + 0.030000s system = 1.150000s CPU (99.0%)
  //  Total File openeing time was:  1.313293s wall, 1.420000s user + 0.130000s system = 1.550000s CPU (118.0%)
  //
  //Reading the file all at once
  //  To parse SDF took:  1.023828s wall, 0.980000s user + 0.030000s system = 1.010000s CPU (98.6%)
  //  Total File openeing time was:  1.191765s wall, 1.330000s user + 0.160000s system = 1.490000s CPU (125.0%)
  //
  //Making things niavly multithreaded
  //  To parse SDF took:  0.864120s wall, 1.140000s user + 0.110000s system = 1.250000s CPU (144.7%)
  //  Total File openeing time was:  0.995905s wall, 1.410000s user + 0.190000s system = 1.600000s CPU (160.7%)
  //
  //With error checking
  //  To parse SDF took:  0.855769s wall, 1.140000s user + 0.110000s system = 1.250000s CPU (146.1%)
  //
  //With current multithreaded implementation:
  //  To parse SDF took:  0.971223s wall, 0.950000s user + 0.020000s system = 0.970000s CPU (99.9%)
  //  Total File openeing time was:  1.102778s wall, 1.230000s user + 0.110000s system = 1.340000s CPU (121.5%)
  //
  //So I think I'll just stick to single threaded parsing for now since it only
  //  increases speed by ~15% to do mutliithreaded, at the cost of memory.
  //
  //TODO: Check whether creating the Measurment objects in a multithreaded
  //      fashion significantly helps things.
  //      Make it so the worker functions in the SpectroscopicDailyFile
  //      namespace dont re-copy all the strings passed in.
  
  
  using namespace SpectroscopicDailyFile;
  
  int occupancy_num = 0, background_num = 0, s1_num = 0, s2_num = 0;
  vector<DailyFileS1Info> s1infos;
  vector< map<string,vector< pair<float,float> > > > detname_to_devpairs;
  
  map<int,int> background_to_s1_num, background_to_s2_num;
  map<int,int> occupancy_num_to_s1_num, occupancy_num_to_s2_num;
  map<int,vector<std::shared_ptr<DailyFileGammaBackground> > > gamma_backgrounds;
  map<int,std::shared_ptr<DailyFileNeutronBackground> > neutron_backgrounds;  //*should* only have one per background number
  map<int, boost::posix_time::ptime > end_background;
  
  map<int,vector<std::shared_ptr<DailyFileGammaSignal> > > gamma_signal;
  map<int,vector<std::shared_ptr<DailyFileNeutronSignal> > > neutron_signal;
  
  map<int,DailyFileEndRecord> end_occupancy;
  
  //We *should* only have one analyzed background per occupancy, so we'll go
  //  with this assumption
  map<int,DailyFileAnalyzedBackground> analyzed_gamma_backgrounds;
  map<int,DailyFileAnalyzedBackground> analyzed_neutron_backgrounds;
  
  
  set<string> detectorNames;

  reset();
  const istream::pos_type orig_pos = input.tellg();
  
#if( DO_SDF_MULTITHREAD )
  
  if( !input.good() )
    return false;
  
  vector<char> filedata;
  input.seekg( 0, ios::end );
  const istream::pos_type eof_pos = input.tellg();
  input.seekg( orig_pos, ios::beg );
  
  const size_t filelength = 0 + eof_pos - orig_pos;
  filedata.resize( filelength + 1 );
  input.read( &filedata[0], filelength );
  if( !input.good() )
  {
    input.clear();
    input.seekg( orig_pos, ios::beg );
    return false;
  }
  
  filedata[filelength] = '\0';
  
  size_t pos = 0;
  const char * const data = &filedata[0];
  
  SpecUtilsAsync::ThreadPool pool;
#else
  string line;
#endif
  
  int nUnrecognizedLines = 0, nLines = 0;
  
#if( DO_SDF_MULTITHREAD )
  while( pos < filelength )
#else
  while( UtilityFunctions::safe_get_line( input, line ) )
#endif
  {
    
#if( DO_SDF_MULTITHREAD )
    const char * const linestart = data + pos;
    while( pos < filelength && data[pos] != '\r' && data[pos] != '\n' )
      ++pos;
    
    const char * const lineend = data + pos;
    const size_t linelen = lineend - linestart;
    
    ++pos;
    if( pos < filelength && data[pos]=='\n' )
      ++pos;
#else
    const size_t linelen = line.length();
    const char * const linestart = line.c_str();
#endif
    
    ++nLines;
    if( linelen < 3 )
      continue;
    
    const string linetype(linestart,2);
    
    //dates are writtn as yyyy-mm-ddThh:mm:ss.000Z-hh:mm
    
    if( linetype == "S1" )
    {
/*     First line of setup parameters
       “S1” is the first line of Setup Parameters.
       The second element on the S1 line is the type of detector, either “NaI” 
       or “HPGe”.
       The third element on the S1 line is the application, either “SPM”, “RDSC”
       , or “MRDIS”.
       The fourth element on the S1 line is the number of channels used per 
       detector.  This is typically 512 for NaI systems or 4096 for HPGe 
       systems.
       The fifth element on the S1 line is the algorithm version number.  This 
       is taken from the <AlgorithmVersion> tag in the ICD-2 file.
       The next series of elements are variable, and are additional setup 
       parameters taken from the ICD-2 file.  They are the children of the
       <AlgorithmParameter> tag in the ICD-2 file, and alternate <ParameterName>
       , <ParameterValue>.
       For example, if a ParameterName was “NSigma” and the ParameterValue was 
       “8”, and the next ParameterName was “Width” and the ParameterValue was 
       “14”, the S1 line would include the text:
       NSigma, 8, Width, 14
       This process continues through all <ParameterName> and <ParameterValue> 
       elements in the ICD-2 file.
*/
      DailyFileS1Info info;
      parse_s1_info( linestart, linelen, info );
      if( !info.success )
      {
        cerr << "load_from_spectroscopic_daily_file(): S1 line invalid" << endl;
        input.clear();
        input.seekg( orig_pos, ios::beg );
        return false;
      }
      
      s1infos.push_back( info );
      
      s1_num = static_cast<int>( s1infos.size() - 1 );
    }else if( linetype == "S2" )
    {
/*     “S2” is the second line of Setup Parameters.  
       This line provides any detector-specific calibration information that is 
       included in the ICD-1 file as a <NonlinearityCorrection> tag.
       The <NonlinearityCorrection> tags are listed in detector order.  So the 
       S2 line should read:
           S2, Aa1, 81, -5, 122, -6, …, Aa2, 81, -4, 122, -6, …, Aa3, 81, -5, …
       These elements are important for properly calibrating the detectors when 
       plotting spectra.
*/
      map<string,vector< pair<float,float> > > devpairs;
      parse_s2_info( linestart, linelen, devpairs );
      detname_to_devpairs.push_back( devpairs );
      
      s2_num = static_cast<int>( detname_to_devpairs.size() - 1 );
    }else if( linetype == "GB" )
    {
/*
      The monitors will produce a background on a regular, periodic interval.  
      The period is currently every 30 minutes; however, it is conceivable that 
      this could change in the future.  These periodic backgrounds are crucial 
      to evaluating the long-term health of the system as well for detecting and 
      troubleshooting intermittent failures.
*/
/*
      Each gamma detector has its own GB line.
      The specific gamma detector is denoted in the second element of the GB 
      line (e.g., Aa1, Ca3).
      The remaining elements are the channel counts in each channel, in order.  
      These are taken from the <ChannelData> child of the same 
      <BackgroundSpectrum> that the timestamp was derived from.  If 
      “zeros compression” is used, the data must be expanded such that each line 
      has 512 (or 4096) channels of data.  The data must also be energy calibrated.
*/
      std::shared_ptr<DailyFileGammaBackground> info( new DailyFileGammaBackground() );
      
#if( DO_SDF_MULTITHREAD )
      pool.post( std::bind( &parse_gamma_background, linestart, linelen, std::ref(*info) ) );
#else
      parse_gamma_background( linestart, linelen, *info );
      
      if( !info->success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing gamma background" << endl;
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( !info->success )
      
      detectorNames.insert( info->detectorName );
#endif
      
      gamma_backgrounds[background_num].push_back( info );
    }else if( linetype == "NB" )
    {
/*    All neutron detectors are combined on a single NB line.
      The second element of the NB line is the duration of the background, in 
      seconds.  To help align columns, the duration in seconds should always be 
      a three-digit number (e.g., 030 for a 30-second background).  The daily 
      file assumes that the duration of the neutron background is identical to 
      the duration of the gamma background, which is supported by operational 
      observations.
      The remaining elements of the NB line are the counts for each neutron 
      detector recorded over the background period.  The detectors are listed in 
      order – Aa1N, then Aa2N, then Ba1N, etc.
*/
      std::shared_ptr<DailyFileNeutronBackground> info( new DailyFileNeutronBackground() );
      
#if( DO_SDF_MULTITHREAD )
      pool.post( std::bind(&parse_neutron_background, linestart, linelen, std::ref(*info) ) );
#else
      parse_neutron_background( linestart, linelen, *info );
      
      if( !info->success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing neutron background" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }
#endif
      
      neutron_backgrounds[background_num] = info;
    }else if( linetype == "BX" )
    {
/*    Then end of the background is denoted by a BX line.
      The second element of the BX line is the timestamp of the background.  The 
      timestamp is the <StartTime> child to the <BackgroundSpectrum> tag in the 
      periodic background files.
*/
      if( linelen < 4 )
      {
        cerr << "load_from_spectroscopic_daily_file(): invalid BX line lenght"
             << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }
      
      const string line( linestart, linelen );
      end_background[background_num] = time_from_string( line.c_str() + 3 );
      background_to_s1_num[background_num] = s1_num;
      background_to_s2_num[background_num] = s2_num;
      
      ++background_num;
    }else if( linetype == "GS" )
    {
/*    Signals – GS and NS
      ICD-1 file contains large amounts of pre- and post-occupancy data, all 
      recorded in 0.1 second intervals.  These high-frequency, long-duration 
      measurements significantly increase the size of the ICD-1 file.  This 
      level of data fidelity is useful for detailed analysis; however, the cost 
      would be too great to use these as daily files, and is not warranted for 
      daily file analysis.  The signal lines in the daily file make two 
      important concessions to reduce overall file size:
        -The only data included in the daily file is occupied data.  Pre- and 
         post-occupancy data are discarded.
        -The data is only recorded in 1 second intervals, not 0.1 second 
         intervals.
      Future versions of this document may relax some of these concessions to 
      include data other than occupied, or at a higher frequency.
*/
/*    The second element of the GS line is the gamma detector name, for example 
      Aa1 or Da3.
      The third element of the GS line is the time chunk number, starting from 
      001.  In this version of the specification, the first ten time slices will 
      be aggregated into the first time chunk (001); the next ten time slices 
      will be aggregated into the second time chunk (002); and so on, resulting 
      in 1 second time chunks.  In the future, it is conceivable that a 
      different time chunk size (0.5 second, 0.2 second, or 2 seconds) may be 
      utilized.  The time chunk number serves as a timestamp within the occupancy.
      The remaining elements of the GS line are the counts in each channel for 
      that detector, aggregated over one second and energy calibrated per the 
      parameters provided in the S1 and S2 lines (and any other source as 
      necessary).  These are taken directly from the <ChannelData> elements.  
      Unfortunately, since these are taken directly from the ICD-1 file, GS line 
      data is not energy calibrated as of this version.
*/
      std::shared_ptr<DailyFileGammaSignal> info( new DailyFileGammaSignal );
      
#if( DO_SDF_MULTITHREAD )
      pool.post( std::bind(&parse_gamma_signal, linestart, linelen, std::ref(*info) ) );
#else
      parse_gamma_signal( linestart, linelen, *info );
      
      if( !info->success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing gamma signal" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( !info->success )
      
      detectorNames.insert( info->detectorName );
#endif
      
      gamma_signal[occupancy_num].push_back( info );
    }else if( linetype == "NS" )
    {
/*    Neutron Signal
      The second element of the NS line is the number of time slices used to
      form one time chunk, represented as a 3 digit number to help align 
      columns.  In this case, since ten time slices contribute to each chunk, 
      the second element of the NS line should read, “010”.  (Again, future 
      versions could change this to 005 or 002 or 020.)
      The third element of the NS line is the time chunk number.  This should be 
      identical to the time chunk number used in the gamma signal.
      The remaining elements of the NS line are the counts from each detector 
      for the one second interval.  These are taken directly from the 
      <ChannelData> elements.  The signals are listed in order of the detectors: 
      Aa1N, Aa2N, Ba1N, and so forth.
*/
      std::shared_ptr<DailyFileNeutronSignal> info( new DailyFileNeutronSignal() );
      
#if( DO_SDF_MULTITHREAD )
      pool.post( std::bind(&parse_neutron_signal,linestart, linelen, std::ref(*info) ) );
#else
      parse_neutron_signal( linestart, linelen, *info );
      
      if( !info->success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing neutron signal" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( !info->success )
#endif
      
      neutron_signal[occupancy_num].push_back( info );
    }else if( linetype == "ID" )
    {
/*    One line is provided for each radionuclide identification.  Even if 
      multiple identifications are made within the same detector subset and time 
      slices, a new line should be provided for each radionuclide identified.  
      If a radionuclide is identified multiple times within the same occupancy 
      (based on different time slices or different detector subsets), a separate
      line should be provided for each ID.
      The second element of the ID line is the radionuclide identified.  The 
      nuclide ID comes from the <NuclideName> tag in the ICD-2 file.
      The next elements of the ID line are the detectors that were used to make
      the identification.  These stem from the <SubsetSampleList> element in the 
      ICD-2 file.
      The next elements of the ID line are the time slices that were used in 
      making the identification.  These are taken directly from the 
      <SubsetSampleList> tag name in the ICD-2 file.
      The ID field will state “NONE” if no radionuclide was identified.  This is 
      made clear if the <AlarmDescription> tag reads “No Alarm”.
*/
      //ToDo, implement this
    }else if( linetype == "AB" )
    {
/*    When evaluating an alarm, the background used by the algorithm is 
      extremely important to capture.  This is provided in the ICD-2 file as a
      common background aggregated over all gamma detectors, over a long period 
      of time – typically 300 seconds.
*/
/*    The second element of the AB line is “Neutron” for the neutron background.
      The third element of the AB line is the duration of the background, in 
      seconds.  This is taken from the <RealTime> child of the
      <BackgroundSpectrum> element in the ICD-2 file, the same as the gamma 
      background.
      The fourth element of the AB line is the <BackgroundCounts> child of the 
      <GrossCountSummed> element.  This is the sum of the counts from all 
      neutron detectors over the background period.
*/
      DailyFileAnalyzedBackground info;
      parse_analyzed_background( linestart, linelen, info );
      
      if( !info.success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing analyzed background" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }
      
      if( info.type == DailyFileAnalyzedBackground::Gamma )
        analyzed_gamma_backgrounds[occupancy_num] = info;
      else
        analyzed_neutron_backgrounds[occupancy_num] = info;
    }else if( linetype == "GX" )
    {
/*    The end of the occupancy is denoted by a GX line.
      The second element on the GX line is the alarm light color, either Red, 
      Yellow, or Green, taken from the <AlarmLightColor> tag in the ICD-2 file.  
      This is useful for categorizing alarms in follow-up analysis.
      The third element on the GX line is the occupancy number, taken from the 
      <occupancyNumber> tag in the ICD-2 file.
      The fourth element on the GX line is the timestamp, taken from the last 
      time stamp in the ICD-1 file.  This is the <StartTime> child of the last 
      <DetectorData> element in the file.  This methodology should also work for 
      the RDSC, which does not record pre- and post-occupancy data.
      The fifth element of the GX line is the filename of the ICD-1 file that 
      provided data for this occupancy.  This information may be useful in case 
      the actual ICD-1 file is needed for additional analysis.
      The sixth element on the GX line is the entry speed, taken from the 
      <vehicleEntrySpeed> tag in the ICD-2 file.
      The seventh element on the GX line is the exit speed, taken from the 
      <vehicleExitSpeed> tag in the ICD-2 file.  This sixth element (exit speed) 
      may or may not exist, if the monitor records it.
*/
      
      DailyFileEndRecord info;
      parse_end_record( linestart, linelen, info );
      
      if( !info.success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing end of record line" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }
      
      end_occupancy[occupancy_num] = info;
      
      occupancy_num_to_s1_num[occupancy_num] = s1_num;
      occupancy_num_to_s2_num[occupancy_num] = s2_num;
      
      ++occupancy_num;
    }else
    {
      string line(linestart, linelen);
      UtilityFunctions::trim( line );
      
      if( !line.empty() )
      {
        ++nUnrecognizedLines;
        const double fracBad = double(nUnrecognizedLines) / nLines;
        if( (nUnrecognizedLines > 10) && (fracBad > 0.1) )
        {
          input.clear();
          input.seekg( orig_pos, ios::beg );
          return false;
        }
      
        cerr << "load_from_spectroscopic_daily_file, unrecognized line begining: "
             << linetype << endl;
      }//if( !line.empty() )
    }//if / else (determine what this line means)
  }//while( UtilityFunctions::safe_get_line( input, line ) )
  
#if( DO_SDF_MULTITHREAD )
  pool.join();
#endif
  
  //Probably not necassary, but JIC
  background_to_s1_num[background_num] = s1_num;
  background_to_s2_num[background_num] = s2_num;
  occupancy_num_to_s1_num[occupancy_num] = s1_num;
  occupancy_num_to_s2_num[occupancy_num] = s2_num;
  
  //TODO: convert so that we are sure we are using the correct setup, incase
  //      there are multiple setup lines
  //Probably just create a struct to hold the information, and parse all the
  //      setups.

  
  //Heres what we have to work with:
  if( s1infos.empty() )
  {
    cerr << "load_from_spectroscopic_daily_file(): Either S1 line missing"
         << endl;
    
    input.clear();
    input.seekg( orig_pos, ios::beg );
    
    return false;
  }//
  

#if( DO_SDF_MULTITHREAD )
  //The bellow two loops are probably quite wasteful, and not necassary
  for( map<int,vector<std::shared_ptr<DailyFileGammaBackground> > >::const_iterator i = gamma_backgrounds.begin();
      i != gamma_backgrounds.end(); ++i )
  {
    for( size_t j = 0; j < i->second.size(); ++j )
      detectorNames.insert( i->second[j]->detectorName );
  }

  for( map<int,vector<std::shared_ptr<DailyFileGammaSignal> > >::const_iterator i = gamma_signal.begin();
      i != gamma_signal.end(); ++i )
  {
    for( size_t j = 0; j < i->second.size(); ++j )
      detectorNames.insert( i->second[j]->detectorName );
  }
#endif  //#if( DO_SDF_MULTITHREAD )
  
  
  map<string,int> detNameToNum;
  int detnum = 0;
  for( const string &name : detectorNames )
    detNameToNum[name] = detnum++;

  vector< std::shared_ptr<Measurement> > backgroundMeasurments, signalMeasurments;
  
  int max_occupancie_num = 0;
  
  for( int occnum = 0; occnum < occupancy_num; ++occnum )
  {
    const map<int,int>::const_iterator s1pos = occupancy_num_to_s1_num.find(occnum);
    const map<int,int>::const_iterator s2pos = occupancy_num_to_s2_num.find(occnum);
    
    if( s1pos == occupancy_num_to_s1_num.end() || s1pos->second >= int(s1infos.size()) )
    {
      cerr << "Serious programing logic error in "
      "load_from_spectroscopic_daily_file(): 0" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }
    
    const DailyFileS1Info &sinfo = s1infos[s1pos->second];
    const map<string,vector< pair<float,float> > > *devpairs = 0;
    if( s2pos != occupancy_num_to_s2_num.end() && s2pos->second<int(detname_to_devpairs.size()))
      devpairs = &(detname_to_devpairs[s2pos->second]);

    const map<int,vector<std::shared_ptr<DailyFileGammaSignal> > >::const_iterator gammaiter
                                                  = gamma_signal.find( occnum );
    if( gammaiter == gamma_signal.end() )
    {
      cerr << "Serious programing logic error in "
      "load_from_spectroscopic_daily_file(): 0.1" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }//if( gammaiter == gamma_signal.end() )
    
    const map<int,vector<std::shared_ptr<DailyFileNeutronSignal> > >::const_iterator neutiter
                                                = neutron_signal.find( occnum );
    
    const map<int,DailyFileEndRecord>::const_iterator endrecorditer
                                                 = end_occupancy.find( occnum );
    if( endrecorditer == end_occupancy.end() )
    {
      cerr << "Serious programing logic error in "
      "load_from_spectroscopic_daily_file(): 0.2" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }//if( endrecorditer == end_occupancy.end() )
    
    const DailyFileEndRecord &endrecord = endrecorditer->second;
    const vector<std::shared_ptr<DailyFileGammaSignal> > &gammas = gammaiter->second;
    const vector<std::shared_ptr<DailyFileNeutronSignal> > *nutsignal = 0;
    if( neutiter != neutron_signal.end() )
      nutsignal = &neutiter->second;
    
    const DailyFileAnalyzedBackground *gammaback = 0, *neutback = 0;
    map<int,DailyFileAnalyzedBackground>::const_iterator backiter;
    backiter = analyzed_gamma_backgrounds.find( occnum );
    if( backiter != analyzed_gamma_backgrounds.end() )
      gammaback = &backiter->second;
    backiter = analyzed_neutron_backgrounds.find( occnum );
    if( backiter != analyzed_neutron_backgrounds.end() )
      neutback = &backiter->second;
    
    //Place the analyzed background into signalMeasurments
    if( gammaback )
    {
      std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
      
      meas->detector_number_    = static_cast<int>( detNameToNum.size() );
      meas->detector_name_      = "sum";
      meas->gamma_counts_       = gammaback->spectrum;
      meas->sample_number_      = 1000*endrecord.occupancyNumber;
      meas->source_type_        = Measurement::Background;
      meas->occupied_           = Measurement::NotOccupied;
      meas->energy_calibration_model_  = Measurement::Polynomial;
      meas->calibration_coeffs_ = sinfo.calibcoefs;
      meas->remarks_.push_back( "Analyzed Background (sum over all detectors" );
      meas->real_time_ = meas->live_time_ = 0.1f*detNameToNum.size()*gammaback->realTime;
      
/*
      meas->start_time_         = endrecord.lastStartTime;
      if( gammas.size() )
      {
        //This is a bit of a hack; I want the analyzed backgrounds to appear
        //  just before the analyzed spectra, so to keep this being the case
        //  we have to falsify the time a bit, because the measurments will get
        //  sorted later on according to start time
        const int totalChunks = gammas[gammas.size()-1].timeChunkNumber;
        
        const DailyFileNeutronSignal *neut = 0;
        if( nutsignal && nutsignal->size() )
          neut = &(*nutsignal)[0];
        
        const float realTime = neut ? 0.1f*neut->numTimeSlicesAgregated : 1.0f;
        const float timecor = realTime * (totalChunks - 0.5);
        const boost::posix_time::seconds wholesec( static_cast<int>(floor(timecor)) );
        const boost::posix_time::microseconds fracsec( static_cast<int>(1.0E6 * (timecor-floor(timecor))) );
        meas->start_time_ -= wholesec;
        meas->start_time_ -= fracsec;
        
        cout << "Background meas->sample_number_=" << meas->sample_number_ << " has time " << meas->start_time_ << endl;
      }//if( gammas.size() )
*/
      
      if( neutback && !neutback->spectrum )
      {
        meas->neutron_counts_ = *neutback->spectrum;
        meas->neutron_counts_sum_ = 0.0;
        for( const float f : meas->neutron_counts_ )
          meas->neutron_counts_sum_ += f;
        meas->contained_neutron_ = true;
      }//if( neutback )
      
      meas->gamma_count_sum_ = 0.0;
      if( !!meas->gamma_counts_ )
      {
        for( const float f : *meas->gamma_counts_ )
          meas->gamma_count_sum_ += f;
      }//if( !!meas->gamma_counts_ )
      
      signalMeasurments.push_back( meas );
//      backgroundMeasurments.push_back( meas );
    }//if( gammaback )
    
    for( size_t i = 0; i < gammas.size(); ++i )
    {
      const DailyFileNeutronSignal *neut = 0;
      const DailyFileGammaSignal &gamma = *gammas[i];
      
#if( DO_SDF_MULTITHREAD )
      if( !gamma.success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing gamma signal" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( !gamma.success )
#endif
      
      if( nutsignal )
      {
        for( size_t j = 0; j < nutsignal->size(); ++j )
        {
          if( (*nutsignal)[j]->timeChunkNumber == gamma.timeChunkNumber )
          {
            neut = ((*nutsignal)[j]).get();
            break;
          }
        }//for( size_t j = 0; j < nutsignal->size(); ++j )
      }//if( nutsignal )
      
#if( DO_SDF_MULTITHREAD )
      if( neut && !neut->success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing neutron signal" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( neut && !neut->success )
#endif
      
      std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
      
      meas->detector_number_    = detNameToNum[gamma.detectorName];
      meas->detector_name_      = gamma.detectorName;
      meas->gamma_counts_       = gamma.spectrum;
      meas->sample_number_      = 1000*endrecord.occupancyNumber + gamma.timeChunkNumber;
      meas->source_type_        = Measurement::Foreground;
      meas->occupied_           = Measurement::Occupied;
      meas->energy_calibration_model_  = Measurement::Polynomial;
      meas->calibration_coeffs_ = sinfo.calibcoefs;
      meas->speed_              = 0.5f*(endrecord.entrySpeed + endrecord.exitSpeed);
      meas->start_time_         = endrecord.lastStartTime;
      meas->remarks_.push_back( "ICD1 Filename: " + endrecord.icd1FileName );
      meas->remarks_.push_back( "Alarm Color: " + endrecord.alarmColor );
      meas->remarks_.push_back( "Occupancy Number: " + std::to_string(endrecord.occupancyNumber) );
      
      max_occupancie_num = std::max( endrecord.occupancyNumber, max_occupancie_num );
      
      meas->gamma_count_sum_ = 0.0;
      if( !!meas->gamma_counts_ )
      {
        for( const float f : *meas->gamma_counts_ )
          meas->gamma_count_sum_ += f;
      }
      
      meas->contained_neutron_ = false;
      meas->live_time_ = meas->real_time_ = 1.0f;
      if( neut )
      {
        meas->live_time_ = meas->real_time_ = 0.1f*neut->numTimeSlicesAgregated;
        
        if( meas->detector_number_ < static_cast<int>(neut->counts.size()) )
        {
          meas->neutron_counts_sum_ = neut->counts[meas->detector_number_];
          meas->neutron_counts_.resize( 1 );
          meas->neutron_counts_[0] = static_cast<float>( meas->neutron_counts_sum_ );
          meas->contained_neutron_ = true;
        }else
        {
          meas->neutron_counts_sum_ = 0.0;
        }
      }//if( neut )
      
      const int totalChunks = gammas[gammas.size()-1]->timeChunkNumber;
      const float dtMeasStart = meas->real_time_ * (totalChunks - 1);
      const float timecor = dtMeasStart * float(totalChunks-gamma.timeChunkNumber)/float(totalChunks);
      const boost::posix_time::seconds wholesec( static_cast<int>(floor(timecor)) );
      const boost::posix_time::microseconds fracsec( static_cast<int>(1.0E6 * (timecor-floor(timecor))) );
      meas->start_time_ -= wholesec;
      meas->start_time_ -= fracsec;
      
      signalMeasurments.push_back( meas );
    }//for( size_t i = 0; i < gammas.size(); ++i )
  }//for( int occnum = 0; occnum < occupancy_num; ++occnum )

  
  
  for( int backnum = 0; backnum < background_num; ++backnum )
  {
    const map<int,int>::const_iterator s1pos = background_to_s1_num.find(backnum);
    const map<int,int>::const_iterator s2pos = background_to_s2_num.find(backnum);
  
    if( s1pos == background_to_s1_num.end() || s1pos->second >= int(s1infos.size()) )
    {
      cerr << "Serious programing logic error in "
              "load_from_spectroscopic_daily_file(): 1" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }
    
    const DailyFileS1Info &sinfo = s1infos[s1pos->second];
    const map<string,vector< pair<float,float> > > *devpairs = 0;
	const int ndets = static_cast<int>( detname_to_devpairs.size() );
    if( s2pos != background_to_s2_num.end() && s2pos->second < ndets )
      devpairs = &(detname_to_devpairs[s2pos->second]);

    const map<int,vector<std::shared_ptr<DailyFileGammaBackground> > >::const_iterator gammaback
                                             = gamma_backgrounds.find(backnum);
    
    const map<int,std::shared_ptr<DailyFileNeutronBackground> >::const_iterator neutback
                                           = neutron_backgrounds.find(backnum);
    
    const map<int, boost::posix_time::ptime >::const_iterator backtimestamp
                                                = end_background.find(backnum);
    
    if( gammaback == gamma_backgrounds.end() )
    {
      cerr << "Serious programing logic error in "
              "load_from_spectroscopic_daily_file(): 1.1" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }
    
    if( backtimestamp == end_background.end() )
    {
      cerr << "Serious programing logic error in "
              "load_from_spectroscopic_daily_file(): 1.2" << endl;
      
      input.clear();
      input.seekg( orig_pos, ios::beg );
      
      return false;
    }
    
    const vector<std::shared_ptr<DailyFileGammaBackground> > &backgrounds = gammaback->second;
    const boost::posix_time::ptime &timestamp = backtimestamp->second;
    
    for( size_t i = 0; i < backgrounds.size(); ++i )
    {
      const DailyFileGammaBackground &back = *backgrounds[i];
      
#if( DO_SDF_MULTITHREAD )
      if( !back.success )
      {
        cerr << "load_from_spectroscopic_daily_file(): "
                "Error Parsing gamma background" << endl;
        
        input.clear();
        input.seekg( orig_pos, ios::beg );
        
        return false;
      }//if( !back.success )
#endif
      
      std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
      
      meas->source_type_        = Measurement::Background;
      meas->detector_name_      = back.detectorName;
      meas->detector_number_    = detNameToNum[back.detectorName];
      meas->gamma_counts_       = back.spectrum;
      meas->start_time_         = timestamp;
      meas->energy_calibration_model_  = Measurement::Polynomial;
      meas->calibration_coeffs_ = sinfo.calibcoefs;
      meas->occupied_           = Measurement::NotOccupied;
      
      meas->sample_number_ = 1000*(max_occupancie_num+1) + backnum;
      
      if( !meas->gamma_counts_ )
        cerr << "Warning, invalid gamma counts" << endl;
      else if( static_cast<int>(meas->gamma_counts_->size()) != sinfo.nchannels )
        cerr << "WArning, mismatch in spectrum size, got "
             << meas->gamma_counts_->size() << " expected " << sinfo.nchannels
             << endl;
      
      if( devpairs )
      {
        map<string,vector< pair<float,float> > >::const_iterator pos
                                      = devpairs->find( meas->detector_name_ );
        meas->deviation_pairs_ = pos->second;
      }//if( devpairs )
      
      meas->gamma_count_sum_ = 0.0;
      if( !!meas->gamma_counts_ )
      {
        for( const float f : *meas->gamma_counts_ )
          meas->gamma_count_sum_ += f;
      }//if( !!meas->gamma_counts_ )
      
      meas->contained_neutron_ = false;
      if( neutback != neutron_backgrounds.end() )
      {
        const DailyFileNeutronBackground &neutbackground = *neutback->second;
        
#if( DO_SDF_MULTITHREAD )
        if( !neutbackground.success )
        {
          cerr << "load_from_spectroscopic_daily_file(): "
                  "Error Parsing neutron background" << endl;
          
          input.clear();
          input.seekg( orig_pos, ios::beg );
          
          return false;
        }
#endif
        
        meas->live_time_ = meas->real_time_ = neutbackground.realTime;
        const int nneutdet = static_cast<int>( neutbackground.counts.size() );
        if( meas->detector_number_ < nneutdet )
        {
          const float counts = neutbackground.counts[meas->detector_number_];
          meas->neutron_counts_.resize( 1 );
          meas->neutron_counts_[0] = counts;
          meas->neutron_counts_sum_ = counts;
          meas->contained_neutron_ = true;
        }//if( meas->detector_number_ < neutbackground.counts.size() )
      }//if( neutback != neutron_backgrounds.end() )
      
      backgroundMeasurments.push_back( meas );
    }//for( size_t i = 0; i < backgrounds.size(); ++i )
  }//for( int backnum = 0; backnum < background_num; ++backnum )
  
  
  for( std::shared_ptr<Measurement> &m : signalMeasurments )
    measurements_.push_back( m );
  
  for( std::shared_ptr<Measurement> &m : backgroundMeasurments )
    measurements_.push_back( m );

  
  for( size_t i = 0; i < s1infos.size(); ++i )
  {
    const DailyFileS1Info &sinfo = s1infos[i];
    remarks_.push_back( "Algorithm Version: " + sinfo.algorithmVersion );
    remarks_.push_back( "Portal Type: " + sinfo.appTypeStr );  //SPM, RDSC, MRDIS
    instrument_type_ = sinfo.detTypeStr;
    
    if( sinfo.appTypeStr == "SPM" )
      instrument_model_ = "ASP";
    else if( sinfo.appTypeStr == "RDSC" )
      instrument_model_ = "Radiation Detector Straddle Carrier";
    else if( sinfo.appTypeStr == "MRDIS" )
      instrument_model_ = "Mobile Radiation Detection and Identification System";

    for( const DailyFileS1Info::params &p : sinfo.parameters )
    {
      string remark = p.name + " = " + p.value;
      if( p.attname.size() && p.attval.size() )
        remark += ", " + p.attname + " = " +p.attval;
      remarks_.push_back( remark );
    }//for( const DailyFileS1Info::params &p : sinfo.parameters )
  }//for( size_t i = 0; i < s1infos.size(); ++i )
  
  
#if( REBIN_FILES_TO_SINGLE_BINNING )
  cleanup_after_load( StandardCleanup | DontChangeOrReorderSamples );
#else
  cleanup_after_load();
#endif
  
  if( properties_flags_ & kNotSampleDetectorTimeSorted )
    cerr << "load_from_spectroscopic_daily_file: kNotSampleDetectorTimeSorted is set" << endl;
  
  if( properties_flags_ & kNotTimeSortedOrder )
    cerr << "load_from_spectroscopic_daily_file: kNotTimeSortedOrder is set" << endl;

  
  if( properties_flags_ & kNotUniqueSampleDetectorNumbers )
    cerr << "load_from_spectroscopic_daily_file: kNotUniqueSampleDetectorNumbers is set" << endl;

  
  if( properties_flags_ & kAllSpectraSameNumberChannels )
    cerr << "load_from_spectroscopic_daily_file: kAllSpectraSameNumberChannels is set" << endl;

  
  if( properties_flags_ & kHasCommonBinning )
    cerr << "load_from_spectroscopic_daily_file: kHasCommonBinning is set" << endl;


  return true;
}//bool load_from_spectroscopic_daily_file( std::istream &input )


namespace
{
  string getAmptekMcaLineInfo( const string &data, const string &heading )
  {
    const size_t pos = data.find( heading );
    if( pos == string::npos )
      return "";
    const size_t end = data.find_first_of( "\r\n", pos );
    if( end == string::npos )
      return "";
    return data.substr( pos + heading.size(), end - pos - heading.size() );
  }
}//anonomous namespace



bool MeasurementInfo::load_from_amptek_mca( std::istream &input )
{
  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  
  char firstline[18];
  input.read( firstline, 17 );
  firstline[17] = '\0';
  
  if( strcmp(firstline, "<<PMCA SPECTRUM>>") != 0 )
    return false;
  
  input.seekg( 0, ios::end );
  const istream::pos_type eof_pos = input.tellg();
  input.seekg( orig_pos, ios::beg );
  
  const size_t filesize = static_cast<size_t>( 0 + eof_pos - orig_pos );
  
  
  //Assume maximum file size of 2.5 MB - this is _way_ more than expected for
  //  even a 16k channel spectrum
  if( filesize > 2.5*1024*1024 )
    return false;
  
  try
  {
    MeasurementShrdPtr meas( new Measurement() );
  
    string filedata;
    filedata.resize( filesize );
    input.read( &filedata[0], filesize );
  
    string lineinfo = getAmptekMcaLineInfo( filedata, "TAG - " );
    if( !lineinfo.empty() )
      remarks_.push_back( "Tag: " + lineinfo );
  
    lineinfo = getAmptekMcaLineInfo( filedata, "DESCRIPTION - " );
    if( !lineinfo.empty() )
      remarks_.push_back( "Description: " + lineinfo );
  
    lineinfo = getAmptekMcaLineInfo( filedata, "GAIN - " );
    if( !lineinfo.empty() )
    {
      float gain;
      if( toFloat(lineinfo,gain) )
      {
        meas->calibration_coeffs_.push_back( 0.0f );
        meas->calibration_coeffs_.push_back( gain );
        meas->energy_calibration_model_ = Measurement::Polynomial;
      }
    }//if( !lineinfo.empty() )
  
  
    lineinfo = getAmptekMcaLineInfo( filedata, "LIVE_TIME - " );
    if( !lineinfo.empty() )
      meas->live_time_ = static_cast<float>( atof( lineinfo.c_str() ) );
  
    lineinfo = getAmptekMcaLineInfo( filedata, "REAL_TIME - " );
    if( !lineinfo.empty() )
      meas->real_time_ = static_cast<float>( atof( lineinfo.c_str() ) );
  
    lineinfo = getAmptekMcaLineInfo( filedata, "START_TIME - " );
    if( !lineinfo.empty() )
      meas->start_time_ = time_from_string( lineinfo.c_str() );
  
    lineinfo = getAmptekMcaLineInfo( filedata, "SERIAL_NUMBER - " );
    if( !lineinfo.empty() )
      instrument_id_ = lineinfo;
  
    size_t datastart = filedata.find( "<<DATA>>" );
    if( datastart == string::npos )
      throw runtime_error( "File doesnt contain <<DATA>> section" );
  
    datastart += 8;
    while( datastart < filedata.size() && !isdigit(filedata[datastart]) )
      ++datastart;
    
    const size_t dataend = filedata.find( "<<END>>", datastart );
    if( dataend == string::npos )
      throw runtime_error( "File doesnt contain <<END>> for data section" );
  
    const size_t datalen = dataend - datastart - 1;
    
    std::shared_ptr<vector<float> > counts( new vector<float>() );
    meas->gamma_counts_ = counts;
    
    const bool success = UtilityFunctions::split_to_floats(
                               filedata.c_str() + datastart, datalen, *counts );
    if( !success )
      throw runtime_error( "Couldnt parse channel data" );
    
    meas->gamma_count_sum_ = 0.0;
    for( const float f : *counts )
      meas->gamma_count_sum_ += f;
    
    const size_t dp5start = filedata.find( "<<DP5 CONFIGURATION>>" );
    if( dp5start != string::npos )
    {
      const size_t dp5end = filedata.find( "<<DP5 CONFIGURATION END>>", dp5start );
      
      if( dp5end != string::npos )
      {
        vector<string> lines;
        const string data = filedata.substr( dp5start, dp5end - dp5start );
        UtilityFunctions::split( lines, data, "\r\n" );
        for( size_t i = 1; i < lines.size(); ++i )
          meas->remarks_.push_back( lines[i] );
      }//if( dp5end != string::npos )
    }//if( dp5start == string::npos )
    
    
    const size_t dppstart = filedata.find( "<<DPP STATUS>>" );
    if( dppstart != string::npos )
    {
      const size_t dppend = filedata.find( "<<DPP STATUS END>>", dppstart );
      
      if( dppend != string::npos )
      {
        vector<string> lines;
        const string data = filedata.substr( dppstart, dppend - dppstart );
        UtilityFunctions::split( lines, data, "\r\n" );
        for( size_t i = 1; i < lines.size(); ++i )
        {
          if( UtilityFunctions::starts_with( lines[i], "Serial Number: " )
              && instrument_id_.size() < 3 )
            instrument_id_ = lines[i].substr( 15 );
          else if( UtilityFunctions::starts_with( lines[i], "Device Type: " ) )
            instrument_model_ = lines[i].substr( 13 );
          else
            remarks_.push_back( lines[i] );
        }
      }//if( dppend != string::npos )
    }//if( dp5start == string::npos )
  
    
    measurements_.push_back( meas );
   
    cleanup_after_load();
    return true;
  }catch( std::exception & )
  {
    reset();
    input.clear();
    input.seekg( orig_pos, ios::beg );
  }
  
  return false;
}//bool MeasurementInfo::load_from_amptek_mca( std::istream &input )



bool MeasurementInfo::load_from_ortec_listmode( std::istream &input )
{
  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  
  try
  {
    //http://www.ortec-online.com/download/List-Mode-File-Formats.pdf
    
    //For reference:
    //  2^21 = 2097152  (e.g. event us clock overflows every 2.097 seconds)
    //  2^31 = 2147483648
    //  2^32 = 4294967296
    //  2147.483647s = 35.79m (e.g. 31 bit clock us overflows every ~36 minutes)
    //  2^30 = 1073741824
    //  1073741.824 = 298 hours (e.g. digiBASE-E ms overflows every ~12 days, can ignore)
    const streampos orig_pos = input.tellg();
    
    double olestartdate;
    uint8_t energy_cal_valid, shape_cal_valid;
    float realtime, livetime;
    float offset, gain, quadratic;
    float shapeoffset, shapegain, shapequad;
    int32_t magicnum, lmstyle, conversiongain, detectoridint;
    char devaddress[81] = {'\0'}, mcb_type[10] = {'\0'}, energyunits[5] = {'\0'};
    char serialnumstr[17] = {'\0'}, txtdcreption[81] = {'\0'}, dummy[9];
    
    
    if( !input.read( (char *)&magicnum, 4) )
      throw runtime_error( "" );  //Failed to read from input stream
    
    if( magicnum != -13 )
      throw runtime_error( "Incorrect leading 4 bytes for .LIS file" );
    
    if( !input.read( (char *)&lmstyle, 4) )
      throw runtime_error( "" );  //Failed to read from input stream
    
    if( lmstyle != 1 && lmstyle != 2 && lmstyle != 4 )
      throw runtime_error( "Unrecognized listmode format" );
    
    if( lmstyle != 1 && lmstyle != 4 )
      throw runtime_error( "Listmode data not in digiBASE/digiBASE-E format (PRO List not supported yet)" );
    
    if( !input.read( (char *)&olestartdate, 8) )
      throw runtime_error( "" );
    
    if( !input.read( devaddress, 80) || !input.read( mcb_type, 9)
       || !input.read( serialnumstr, 16) || !input.read( txtdcreption, 80) )
      throw runtime_error( "" );
    
    if( !input.read( (char *)&energy_cal_valid, 1 ) || !input.read(energyunits,4)
       || !input.read( (char *)&offset, 4 ) || !input.read( (char *)&gain, 4 )
       || !input.read( (char *)&quadratic, 4 ) )
      throw runtime_error( "" );  //Failed reading energy calibration
    
    if( !input.read( (char *)&shape_cal_valid, 1 )
       || !input.read( (char *)&shapeoffset, 4 )
       || !input.read( (char *)&shapegain, 4 )
       || !input.read( (char *)&shapequad, 4 ) )
      throw runtime_error( "" ); //Failed reading shape calibration coefficents
    
    if( !input.read( (char *)&conversiongain, 4) )
      throw runtime_error( "" );
    
    if( !input.read( (char *)&detectoridint, 4) )
      throw runtime_error( "" );
    
    if( !input.read( (char *)&realtime, 4) || !input.read( (char *)&livetime, 4) )
      throw runtime_error( "" );
    
    if( !input.read(dummy, 9) )
      throw runtime_error( "" );
    
    assert( input.tellg() == (orig_pos + streampos(256)) );
  
    size_t ninitialbin = 1024;
    switch( lmstyle )
    {
      case 1: ninitialbin = 1024; break;
      case 2: ninitialbin = 1024; break;  //16k?
      case 4: ninitialbin = 2048; break;
    }
    
    std::shared_ptr< vector<float> > histogram = std::make_shared< vector<float> >( ninitialbin, 0.0f );
    
    uint32_t event;
    
    if( lmstyle == 1 )
    {
      //We need to track overflows in the 31bit microsecond counter, so we will
      //  check if the current 31bit timestamp is less than the previous, and if
      //  so know a overflow occured, and add 2^31 to timeepoch.
      uint32_t previous_time = 0;
      
      //Incase real time isnt excplicitly given in the header we will grab the
      //  first and last events timetampts
      uint64_t firsttimestamp = 0, lasttimestamp = 0;
      
      //To track measurments longer than 35.79 minutes we need to keep track of
      //  more than a 31 bit clock.
      uint64_t timeepoch = 0;
      
      //Bits 20 through 31 of the timestamp (e.g. the part not given with actual
      //  hits)
      uint32_t time_msb = 0;
      
      //It appears that events with a zero value of the 21 bit timestamp will be
      //  transmitted before the 31 bit clock update (it will be the next 32bit
      //  event), so we have to add 2^21 to these zero clock value hits - however
      //  since these events are rare to investigate, I'm adding in an additional
      //  check to make sure that the 31 bit clock wasnt just sent, since I cant be
      //  sure the ordering is always consistent.
      bool prev_was_timestamp = false;
      
      //First, look for timestamp.  If first two timestamps most significant bits
      //  (20 though 31) are different, then the data is starting with the time_msb
      //  first given in file.  If first two timestamps are the same, then begging
      //  timestamp is this value minus 2^21.
      uint32_t firsttimestamps[2] = { 0 };
      for( int i = 0; i < 2 && input.read((char *)&event, 4); )
      {
        if( event > 0x7fffffff )
          firsttimestamps[i++] = (event & 0x7fe00000);
      }
    
      if( firsttimestamps[0] && firsttimestamps[0] == firsttimestamps[1] )
        time_msb = firsttimestamps[0] - 2097152;
      else
        time_msb = firsttimestamps[0];
      
      input.seekg( orig_pos + streampos(256) );
      
      
      for( int64_t eventnum = 0; input.read((char *)&event, 4); ++eventnum )
      {
        if( event <= 0x7fffffff )
        {
          //Bits   Description
          //31     0 for event
          //30-21  Amplitude of pulse (10 bits)
          //20-0   Time in microseconds that the event occured (21 bits)
            
          //10-bit ADC value
          uint32_t amplitude = (uint32_t) ((event & 0x7fe00000) >> 21);
            
          //21-bit timestamp
          uint32_t time_lsb = (event & 0x001fffff);
          //Correct for time_lsb with value zero
          time_lsb = ((time_lsb || prev_was_timestamp) ? time_lsb : 2097152);
          const uint64_t timestamp = timeepoch + time_msb + time_lsb;
          
          if( amplitude > 16384 )
            throw runtime_error( "To high of a channel number" );
            
          amplitude = (amplitude > 0 ? amplitude-1 : amplitude);
          
          if( amplitude >= histogram->size() )
          {
            const size_t next_power_of_two = std::pow(2, std::ceil(log(amplitude)/log(2)));
            histogram->resize( next_power_of_two, 0.0f );
          }
            
          ++((*histogram)[amplitude]);
            
          firsttimestamp = (firsttimestamp ? firsttimestamp : timestamp);
          lasttimestamp = timestamp;
            
          prev_was_timestamp = false;
        }else
        {
          //	The	number	rolls	over	to	0	every 2.097152 seconds. In	order	to	track
          //  the rollovers, a “time only” event is sent from the digiBASE to the
          //  computer every 1.048576 seconds.
            
          //Bits   Description
          //31     1 for time
          //30-0   Current time in microseconds
          const uint32_t this_time = (uint32_t) (event & 0x7fffffff);
          
          if( this_time < previous_time )
            timeepoch += 2147483648u;
          previous_time = this_time;
          time_msb = (this_time & 0xffe00000);
          prev_was_timestamp = true;
        }//if( a hit ) / else( a timestamp )
      }//for( int64_t eventnum = 0; input.read((char *)&event, 4); ++eventnum )
      
      
      if( realtime == 0.0 )
        realtime = 1.0E-6*(lasttimestamp - firsttimestamp);
      if( livetime == 0.0 )
        livetime = 1.0E-6*(lasttimestamp - firsttimestamp);
    }else if( lmstyle == 4 )
    {
      uint32_t firstlivetimes[2] = { 0 }, firstrealtimes[2] = { 0 };
      
      for( int i = 0, j = 0, k = 0; (i < 2 || j < 2) && input.read((char *)&event, 4); ++k)
      {
        const bool msb = (event & 0x80000000);
        const bool ssb = (event & 0x40000000);
        
        if( msb && ssb )
        {
          //ADC Word
        }else if( msb )
        {
          //RT Word, 30-bit Real Time counter in 10 mS Ticks
          firstrealtimes[i++] = (event & 0x3FFFFFFF);
        }else if( ssb )
        {
          //LT Word, 30 Bit Live Time counter in 10 mS Ticks
          firstlivetimes[j++] = (event & 0x3FFFFFFF);
        }else
        {
          //Ext Sync, 13-bit Real Time counter in 10 mS Ticks, 17-bit time pre-scale in 80 nS ticks
        }
      }
      
      uint64_t firsttimestamp = 0, lasttimestamp = 0;
      uint64_t realtime_ns = 10*1000000 * uint64_t(firstrealtimes[0]);
      uint64_t livetime_ns = 10*1000000 * uint64_t(firstlivetimes[0]);
      
      //if( firsttimestamps[0] && firsttimestamps[0] == firsttimestamps[1] )
        //time_msb = firsttimestamps[0] - 2097152;
      //else
        //time_msb = firsttimestamps[0];
      
      
      input.seekg( orig_pos + streampos(256) );
      for( int64_t eventnum = 0; input.read((char *)&event, 4); ++eventnum )
      {
        //Untested!
        if( (event & 0xC0000000) == 0xC0000000 )
        {
          //From DIGIBASE-E-MNL.pdf:
          //  Data[31:30] = 3
          //  Data[29] = Memory Routing bit
          //  Data[28:17] = Conversion Data[11:0]
          //  Data[16:0] = RealTime PreScale [17:1]
          
          //ADC Word: 17-bit time pre-scale in 80 nS Ticks, 13-bit ADC channel number
          //          const uint32_t adc_value = (event & 0x1FFF);
          const uint32_t ticks = (event & 0x0001FFFF);
          const uint64_t timestamp_ns = realtime_ns + 80*ticks;
          
          firsttimestamp = (!firsttimestamp ? timestamp_ns : firsttimestamp);
          lasttimestamp = timestamp_ns;
          
          
          //Have 2048 channels
          uint32_t amplitude = (uint32_t(event & 0x0FFE0000) >> 17);
          
          if( amplitude > 16384 )
            throw runtime_error( "To high of a channel number" );
          
          if( amplitude >= histogram->size() )
          {
            const size_t next_power_of_two = std::pow(2, std::ceil(log(amplitude)/log(2)));
            histogram->resize( next_power_of_two, 0.0f );
          }
          
          ++((*histogram)[amplitude]);
        }else if( (event & 0x80000000) == 0x80000000 )
        {
          //RT Word: 30-bit Real Time counter in 10 mS Ticks
          realtime_ns = 10*1000000*uint64_t(event & 0x3FFFFFFF);
        }else if( (event & 0x40000000) == 0x40000000 )
        {
          //LT Word: 30 Bit Live Time counter in 10 mS Ticks
          livetime_ns = 10*1000000*uint64_t(event & 0x3FFFFFFF);
        }else if( ((event ^ 0xC0000000) >> 30) == 0x00000003 )
        {
          /*
           //Ext Sync: 13-bit Real Time counter in 10 mS Ticks; 17-bit time pre-scale in 80 nS ticks
           //The Ext Sync words contain the value of the external input pulse counters.
           // The external pulse counters count the positive pulses at the external
           // input of the spectrometer. The sync time is calculated by adding the
           //  real time counter to the time pre-scale value.
           const uint32_t tenms_ticks = (uint32_t(event & 0x3FFE0000) >> 17);
           const uint32_t ns_ticks = (event & 0x1FFFF);
           
           
           from: https://oaktrust.library.tamu.edu/bitstream/handle/1969.1/ETD-TAMU-2012-08-11450/CraneWowUserManual.pdf?sequence=3&isAllowed=y
           Care must be taken as the real time
           portion of the stamp, found at bits 12-0, will reset every 80 stamps. This means that the sync time
           stamps roll over every 8 seconds. Note that the time stamp is in units of 100ms, such that a reading of
           “10” would equal one second. The remaining bits are a prescale in the units of 80ns that can be used to
           pinpoint the sync pulse in relation to the digiBASE-Es own clock. See the Sync Gating section of
           Explanations of Terms for more details.
           */
        }
      }//for( int64_t eventnum = 0; input.read((char *)&event, 4); ++eventnum )
      
      //realtime=33.02, from hits->32.62
      //livetime=33.02, from hits->32.2
      //cout << "realtime=" << realtime << ", from hits->" << (1.0E-9*(lasttimestamp - firsttimestamp)) << endl;
      //cout << "livetime=" << realtime << ", from hits->" << (1.0E-9*(livetime_ns - (10*1000000 * uint64_t(firstlivetimes[0])))) << endl;
      
      if( realtime == 0.0 ) //not exact, but close enough
        realtime = 1.0E-9*(lasttimestamp - firsttimestamp);
      
      if( livetime == 0.0 ) //get us within ~20ms live time calc, close enough!
        livetime = 1.0E-9*(livetime_ns - (10*1000000 * uint64_t(firstlivetimes[0])));
    }else if( lmstyle == 2 )
    {
      assert( 0 );
      //This one is pretty complicated, so I would definetly need some example
      //  data to test against.
    }else
    {
      assert( 0 );
    }//if( lmstyle == 1 ) / else
    
    
    
    const float gammasum = std::accumulate( histogram->begin(), histogram->end(), 0.0f, std::plus<float>() );
    if( gammasum < 1.0 && realtime == 0.0f )
      throw runtime_error( "" ); //"Empty listmode file"
    
    std::shared_ptr<Measurement> meas = std::make_shared<Measurement>();
    
    meas->live_time_ = livetime;
    meas->real_time_ = realtime;
    meas->contained_neutron_ = false;
    meas->sample_number_ = 1;
    meas->occupied_ = Measurement::UnknownOccupancyStatus;
    meas->gamma_count_sum_ = gammasum;
    meas->neutron_counts_sum_ = 0.0;
    meas->speed_ = 0.0;  //in m/s
    meas->detector_name_ = ((lmstyle==1) ? "digiBASE" : "digiBASE-E");
    meas->detector_number_ = 0;
    meas->detector_type_ = meas->detector_name_ + " ListMode data";
    meas->quality_status_ = Measurement::Missing;
    meas->source_type_ = Measurement::UnknownSourceType;
    meas->energy_calibration_model_ = Measurement::Polynomial;
    
    //std::vector<std::string>  remarks_;
    
    if( olestartdate > 0 )
      meas->start_time_ = datetime_ole_to_posix( olestartdate );
    
    vector<float> energycoef;
    energycoef.push_back( offset );
    energycoef.push_back( gain );
    if( quadratic != 0.0 )
      energycoef.push_back( quadratic );
    const DeviationPairVec devpairs;
    
    
    //Should check energyunits=="keV"
    if( energy_cal_valid && (gain != 0 || quadratic != 0)
       && calibration_is_valid( Measurement::Polynomial, energycoef, devpairs, histogram->size() ) )
    {
      meas->calibration_coeffs_ = energycoef;
    }
    
    meas->deviation_pairs_ = devpairs;
    //meas->channel_energies_ = ...;
    meas->gamma_counts_ = histogram;
    //meas->neutron_counts_ = ...;
    //meas->latitude_;  //set to -999.9 if not specified
    //meas->longitude_; //set to -999.9 if not specified
    //meas->position_time_;
    

    meas->title_ = txtdcreption;
    
    instrument_type_ = "Spectroscopic Personal Radiation Detector";
    manufacturer_ = "ORTEC";
    instrument_model_ = ((lmstyle==1) ? "digiBASE" : "digiBASE-E");
    instrument_id_ = serialnumstr;
    if( instrument_id_.empty() && detectoridint )
    {
      char buffer[32];
      snprintf( buffer, sizeof(buffer), "%i", int(detectoridint) );
      instrument_id_ = buffer;
    }

    if( strlen(txtdcreption) )
      remarks_.push_back( "Description: " + string(txtdcreption) );
    if( strlen(devaddress) )
      remarks_.push_back( "Device Address: " + string(devaddress) );
    if( strlen(mcb_type) )
      remarks_.push_back( "MCB Type: " + string(mcb_type) );
    
    measurements_.push_back( meas );

    cleanup_after_load();
    
    if( measurements_.empty() )
      throw std::runtime_error( "no measurments" );
    
    return true;
  }catch( std::exception & )
  {
    reset();
    input.clear();
    input.seekg( orig_pos, ios::beg );
  }
  
  return false;
}//bool load_from_ortec_listmode( std::istream &input )


namespace
{
//Function to find a specific block (e.g., 512kb) of information in a CNF file.
bool findCamberraBlock( const uint8_t B, const size_t SBeg, size_t &pos,
                        std::istream &input, const size_t streamsize )
{
  for( pos = SBeg; (pos+512) < streamsize; pos += 512 )
  {
    input.seekg( pos, std::ios::beg );
    uint8_t bytes[2];
    input.read( (char *)bytes, 2 );
    if( (bytes[0] == B) && (bytes[1] == 0x20) )
      return true;
  }
  return false;
}//findCamberraBlock(...)

//Function to read a 32bit float (i.e., energy calibration coefficient) from a
//  Canbera CNF file
float readCamberaFloat( std::istream &input )
{
  float val;
  uint8_t Buf[4], TmpBuf[4];
  input.read( (char *)Buf, 4 );
  TmpBuf[0] = Buf[2];
  TmpBuf[1] = Buf[3];
  TmpBuf[2] = Buf[0];
  TmpBuf[3] = Buf[1];
  
  memcpy( &val, TmpBuf, sizeof(val) );
  
  return 0.25f*val;
}//float readCamberaFloat( std::istream &input )

}//namespace to help read Canberra files

bool MeasurementInfo::load_camberra_cnf_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;
  const bool loaded = load_from_camberra_cnf( file );
  if( loaded )
    filename_ = filename;
  
  return loaded;
}//bool load_camberra_cnf_file( const std::string &filename );


bool MeasurementInfo::load_from_camberra_cnf( std::istream &input )
{
  //Content removed 20180509 pending input from expert
  return false;
}//bool load_from_camberra_cnf( std::istream &input )


bool MeasurementInfo::load_tracs_mps_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;
  const bool loaded = load_from_tracs_mps( file );
  if( loaded )
    filename_ = filename;
  
  return loaded;
}//bool load_tracs_mps_file( const std::string &filename )


bool MeasurementInfo::load_from_tracs_mps( std::istream &input )
{
/*
 Cabin Data
 Byte offset	Size	Description
 0	8	Memory address
 8	4	Memory address
 12	4	Connect Status
 16	4	Event
 20	4	Neutron AlarmLevel
 24	4	Gamma Alarm Level
 28	4	Ratio Alarm Level
 32	8	Latitude
 40	8	Longitude
 48	4	GPS Time of Day
 52	4	#1 pod status
 56	4	#2 pod status
 60	4	#1 det status
 64	4	#2 det status
 68	4	#3 det status
 72	4	#4 det status
 76	4	Index Number
 80	4	Neutron GC
 84	4	Gamma GC
 88	2048	Sum Spectra
 2136	4	Pod1 Index Number
 2140	4	Pod1 deltaTau
 
 2144	4	Pod1 Det1 Neutron GC
 2148	4	Pod1 Det2 Neutron GC
 
 2152	4	Pod1 Det1 Gamma GC
 2156	4	Pod1 Det2 Gamma GC
 2160	4	Pod1 Det1 DAC
 2164	4	Pod1 Det2 DAC
 2168	4	Pod1 Det1 calibration Peak
 2172	4	Pod1 Det2 calibration Peak
 2176	4	Pod1 Det1 calibration peak found
 2180	4	Pod1 Det2 calibration peak found
 
 2184	2048	Pod1 Det1 spectra
 4232	2	Pod1 Det1 clock time
 4234	2	Pod1 Det1 dead time
 4236	2	Pod1 Det1 live time
 
 4238	2048	Pod1 Det2 spectra
 6286	2	Pod1 Det2 clock time
 6288	2	Pod1 Det2 dead time
 6290	2	Pod1 Det2 live time
 
 6292	4	Pod2 Index Number
 6296	4	Pod2 deltaTau
 
 6300	4	Pod2 Det1 Neutron GC
 6304	4	Pod2 Det2 Neutron GC
 
 6308	4	Pod2 Det1 Gamma GC
 6312	4	Pod2 Det2 Gamma GC
 6316	4	Pod2 Det1 DAC
 6320	4	Pod2 Det2 DAC
 6324	4	Pod2 Det1 calibration Peak
 6328	4	Pod2 Det2 calibration Peak
 6332	4	Pod2 Det1 calibration peak found
 6336	4	Pod2 Det2 calibration peak found
 
 6340	2048	Pod2 Det1 spectra
 8388	2	Pod2 Det1 clock time
 8390	2	Pod2 Det1 dead time
 8392	2	Pod2 Det1 live time
 
 8394	2048	Pod2 Det2 spectra
 10442	2	Pod2 Det2 clock time
 10444	2	Pod2 Det2 dead time
 10446	2	Pod2 Det2 live time
 
 10448	4	Radar Altimeter
 10452	128	GPS String
 10580	8	GPS Source
 10588	6	GPS Age
 10594	3	GPS Num SV
 10597
*/
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  const size_t samplesize = 10597;
  
  if( !input.good() )
    return false;
  
  const istream::pos_type orig_pos = input.tellg();
  input.seekg( 0, ios::end );
  const istream::pos_type eof_pos = input.tellg();
  input.seekg( orig_pos, ios::beg );

  const size_t filesize = static_cast<size_t>( 0 + eof_pos - orig_pos );
  const size_t numsamples = filesize / samplesize;
  const size_t leftoverbytes = (filesize % samplesize);
  
  if( leftoverbytes )
    return false;
  
  try
  {
    for( size_t sample = 0; sample < numsamples; ++sample )
    {
      double lat, lon;
      const size_t startpos = sample * samplesize;
      uint32_t gpsTOD, indexNum, overallGammaGC, overallNeutronGC, radarAltimiter;
    
      if( !input.seekg(startpos+32, ios::beg) )
        throw runtime_error( "Failed seek 1" );
      if( !input.read( (char *)&lat, sizeof(lat) ) )
        throw runtime_error( "Failed read 2" );
      if( !input.read( (char *)&lon, sizeof(lon) ) )
        throw runtime_error( "Failed read 3" );
      if( !input.read( (char *)&gpsTOD, sizeof(gpsTOD) ) )
        throw runtime_error( "Failed read 4" );
      if( !input.seekg(startpos+76, ios::beg) )
        throw runtime_error( "Failed read 5" );
      if( !input.read( (char *)&indexNum, sizeof(indexNum) ) )
        throw runtime_error( "Failed read 6" );
      if( !input.read( (char *)&overallNeutronGC, sizeof(overallNeutronGC) ) )
        throw runtime_error( "Failed read 7" );
      if( !input.read( (char *)&overallGammaGC, sizeof(overallGammaGC) ) )
        throw runtime_error( "Failed read 8" );
  
      char gpsstr[129];
      gpsstr[128] = '\0';
      if( !input.seekg(startpos+10448, ios::beg) )
        throw runtime_error( "Failed seek 9" );
      if( !input.read( (char *)&radarAltimiter, sizeof(radarAltimiter) ) )
        throw runtime_error( "Failed read 10" );
      if( !input.read( gpsstr, 128 ) )
        throw runtime_error( "Failed read 11" );
      
      for( size_t i = 0; i < 4; ++i )
      {
        const char *title;
        size_t datastart, neutrongc, gammagc, detstatus;
     
        switch( i )
        {
          case 0:
            detstatus = 60;
            datastart = 2184;
            gammagc   = 2152;
            neutrongc = 2144;
            title = "Pod 1, Det 1";
          break;
        
          case 1:
            detstatus = 64;
            datastart = 4238;
            gammagc   = 2156;
            neutrongc = 2148;
            title = "Pod 1, Det 2";
          break;
    
          case 2:
            detstatus = 68;
            datastart = 6340;
            gammagc   = 6308;
            neutrongc = 6300;
            title = "Pod 2, Det 1";
          break;
        
          case 3:
            detstatus = 72;
            datastart = 8394;
            gammagc   = 6312;
            neutrongc = 6304;
            title = "Pod 2, Det 2";
          break;
        }//switch( i )
    
        uint32_t neutroncount;
        uint16_t channeldata[1024];
        uint16_t realtime, livetime, deadtime;
        uint32_t gammaGC, detDAC, calPeak, calPeakFound, status, dummy;
      
        if( !input.seekg(startpos+detstatus, ios::beg) )
          throw runtime_error( "Failed seek 12" );
        if( !input.read( (char *)&status, sizeof(status) ) )
          throw runtime_error( "Failed read 13" );
      
        if( !input.seekg(startpos+gammagc, ios::beg) )
          throw runtime_error( "Failed seek 14" );
    
        if( !input.read( (char *)&gammaGC, sizeof(gammaGC) ) )
          throw runtime_error( "Failed read 15" );
        if( !input.read( (char *)&dummy, sizeof(dummy) ) )
          throw runtime_error( "Failed read 16" );
        
        if( !input.read( (char *)&detDAC, sizeof(detDAC) ) )
          throw runtime_error( "Failed read 17" );
        if( !input.read( (char *)&dummy, sizeof(dummy) ) )
          throw runtime_error( "Failed read 18" );
      
        if( !input.read( (char *)&calPeak, sizeof(calPeak) ) )
          throw runtime_error( "Failed read 19" );
        if( !input.read( (char *)&dummy, sizeof(dummy) ) )
          throw runtime_error( "Failed read 20" );
      
        if( !input.read( (char *)&calPeakFound, sizeof(calPeakFound) ) )
          throw runtime_error( "Failed read 21" );
        if( !input.read( (char *)&dummy, sizeof(dummy) ) )
          throw runtime_error( "Failed read 22" );
      
        if( !input.seekg(startpos+datastart, ios::beg) )
          throw runtime_error( "Failed seek 23" );
      
        if( !input.read( (char *)channeldata, sizeof(channeldata) ) )
          throw runtime_error( "Failed read 24" );
    
        if( !input.read( (char *)&realtime, sizeof(realtime) ) )
          throw runtime_error( "Failed read 25" );
    
      //if realtime == 6250, then its about 1 second... - I have no idea what these units means (25000/4?)
      
        if( !input.read( (char *)&deadtime, sizeof(deadtime) ) )
          throw runtime_error( "Failed read 26" );
    
        if( !input.read( (char *)&livetime, sizeof(livetime) ) )
          throw runtime_error( "Failed read 27" );
    
        if( !input.seekg(startpos+neutrongc, ios::beg) )
          throw runtime_error( "Failed seek 28" );
    
        if( !input.read( (char *)&neutroncount, sizeof(neutroncount) ) )
          throw runtime_error( "Failed read 29" );
      
      
        MeasurementShrdPtr m( new Measurement() );
        m->live_time_ = livetime / 6250.0f;
        m->real_time_ = realtime / 6250.0f;
        m->contained_neutron_ = (((i%2)!=1) || neutroncount);
        m->sample_number_ = static_cast<int>( sample + 1 );
        m->occupied_ = Measurement::UnknownOccupancyStatus;
        m->gamma_count_sum_ = 0.0;
        m->neutron_counts_sum_ = neutroncount;
//        m->speed_ = ;
        m->detector_name_ = title;
        m->detector_number_ = static_cast<int>( i );
//        m->detector_type_ = "";
        m->quality_status_ = (status==0 ? Measurement::Good : Measurement::Suspect);
        m->source_type_  = Measurement::UnknownSourceType;
        m->energy_calibration_model_ = Measurement::Polynomial;
//        m->start_time_ = ;
        m->calibration_coeffs_.push_back( 0.0f );
//        m->calibration_coeffs_.push_back( 3.0 );
        m->calibration_coeffs_.push_back( 1460.0f / calPeakFound );
//        m->channel_energies_  //dont need to fill out here
      
        vector<float> *gammacounts = new vector<float>( 1024 );
        m->gamma_counts_.reset( gammacounts );
        for( size_t i = 0; i < 1024; ++i )
        {
          const float val = static_cast<float>( channeldata[i] );
          m->gamma_count_sum_ += val;
          (*gammacounts)[i] = val;
        }
      
        if( m->contained_neutron_ )
          m->neutron_counts_.resize( 1, static_cast<float>(neutroncount) );

        m->latitude_ = lat;
        m->longitude_ = lon;
//        m->position_time_ = ;//
        m->title_ = title;
      
        measurements_.push_back( m );
      }//for( size_t i = 0; i < 4; ++i )
    }//for( size_t sample = 0; sample < numsamples; ++sample )
    
    cleanup_after_load();
    
    if( measurements_.empty() )
      throw std::runtime_error( "no measurments" );
  }catch( std::exception &e )
  {
    //cerr << SRC_LOCATION << "\n\tCaught: " << e.what() << endl;
    input.clear();
    input.seekg( orig_pos, ios::beg );
    reset();
    return false;
  }//try / catch
  
  return true;
}//bool load_from_tracs_mps( std::istream &input )


bool MeasurementInfo::load_aram_file( const std::string &filename )
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  reset();
  
  ifstream file( filename.c_str(), ios_base::binary|ios_base::in );
  if( !file.is_open() )
    return false;
  const bool loaded = load_from_aram( file );
  if( loaded )
    filename_ = filename;
  
  return loaded;
}//bool load_aram_file( const std::string &filename )


bool MeasurementInfo::load_from_aram( std::istream &input )
{
  if( !input )
    return false;
  
  const istream::iostate origexceptions = input.exceptions();
  input.exceptions( istream::goodbit );  //make so stream never throws
  const istream::pos_type start_pos = input.tellg();
  
  try
  {
    string line;
    if( !UtilityFunctions::safe_get_line( input, line, 1024 ) )
      throw std::runtime_error( "failed to get line" );
    
    
    
  }catch( std::exception &e )
  {
    reset();
    input.seekg( start_pos );
    input.clear();
  }//try / catch
  
  input.exceptions( origexceptions );
  
  
  return false;
}//bool load_from_aram( std::istream &input )



bool Measurement::write_txt( std::ostream& ostr ) const
{
  const char *endline = "\r\n";
  char buffer[128];
  
  ostr << endline << endline;
  
  for( size_t i = 0; i < remarks_.size(); ++i )
  {
    string remark = remarks_[i];
    if( i == 0 )
    {
      if( (remark.find( "Survey" ) == string::npos) && sample_number_>=0 )
      {
        snprintf( buffer, sizeof(buffer), " Survey %i ", sample_number_ );
        remark += buffer;
      }
      
      const string found_name = detector_name_from_remark( remark );
      if( found_name.empty() && !detector_name_.empty() )
        remark += " " + detector_name_ + " ";
      
      if( (remark.find( "Speed" ) == string::npos) && (speed_>0.000000001) )
      {
        snprintf( buffer, sizeof(buffer), " Speed %f m/s", speed_ );
        remark += buffer;
      }
    }//if( i == 0 )
    
    ostr << "Remark: " << remark << endline;
  }//for( size_t i = 0; i < remarks_.size(); ++i )
  
  if( !start_time_.is_special() )
    ostr << "StartTime " << start_time_ << "" << endline;
  ostr << "LiveTime " << live_time_ << " seconds" << endline;
  ostr << "RealTime " << real_time_ << " seconds" << endline;
  ostr << "SampleNumber " << sample_number_ << endline;
  ostr << "DetectorName " << detector_name_ << endline;
  ostr << "DetectorType " << detector_type_ << endline;
  
  if( has_gps_info() )
  {
    ostr << "Latitude: " << latitude_ << endline;
    ostr << "Longitude: " << longitude_ << endline;
    if( !position_time_.is_special() )
      ostr << "Position Time: "
      << UtilityFunctions::to_iso_string(position_time_) << endline;
  }
  
  ostr << "EquationType ";
  switch( energy_calibration_model_ )
  {
    case Polynomial:          ostr << "Polynomial"; break;
    case FullRangeFraction:   ostr << "FullRangeFraction"; break;
    case LowerChannelEdge:    ostr << "LowerChannelEdge"; break;
    case UnknownEquationType: ostr << "Unknown"; break;
  }//switch( energy_calibration_model_ )
  
  
  ostr << endline << "Coefficients ";
  for( size_t i = 0; i < calibration_coeffs_.size(); ++i )
    ostr << (i ? " " : "") << calibration_coeffs_[i];
  
  if( energy_calibration_model_ == LowerChannelEdge && calibration_coeffs_.empty()
     && !!channel_energies_ && channel_energies_->size() )
  {
    for( size_t i = 0; i < channel_energies_->size(); ++i )
      ostr << (i ? " " : "") << (*channel_energies_)[i];
  }
  
  ostr << endline;
  
  if( contained_neutron_ )
    ostr << "NeutronCount " << neutron_counts_sum_ << endline;
  
  //  ostr "Channel" << " "
  //       << setw(12) << ios::left << "Energy" << " "
  //       << setw(12) << ios::left << "Counts" << endline;
  ostr << "Channel" << " "
  "Energy" << " "
  "Counts" << endline;
  const size_t nChannel = gamma_counts_->size();
  
  for( size_t i = 0; i < nChannel; ++i )
  {
    //    ostr << setw(12) << ios::right << i
    //         << setw(12) << ios::right << channel_energies_->at(i)
    //         << setw(12) << ios::right << gamma_counts_->operator[](i)
    //         << endline;
    ostr << i << " " << channel_energies_->at(i)
    << " " << gamma_counts_->operator[](i)
    << endline;
  }//for( size_t i = 0; i < compressed_counts.size(); ++i )
  
  ostr << endline;
  
  return !ostr.bad();
}//bool write_txt( std::ostream& ostr ) const

bool MeasurementInfo::write_txt( std::ostream& ostr ) const
{
  std::unique_lock<std::recursive_mutex> scoped_lock( mutex_ );
  
  const char *endline = "\r\n";
  ostr << "Original File Name: " << filename_ << endline;
  ostr << "File save from SrbSpectrumViewer compiled on "
  << COMPILE_DATE_AS_INT << endline;
  ostr << "TotalGammaLiveTime: " << gamma_live_time_ << " seconds" << endline;
  ostr << "TotalRealTime: " << gamma_real_time_ << " seconds" << endline;
  ostr << "TotalGammaCounts: " << gamma_count_sum_ << " seconds" << endline;
  ostr << "TotalNeutron: " << neutron_counts_sum_ << " seconds" << endline;
  if( instrument_id_.size() )
    ostr << "Serial number " << instrument_id_ << endline;
  
  for( const string &remark : remarks_ )
    ostr << "Remark: " << remark << endline;
  
  for( const std::shared_ptr<const Measurement> m : measurements_ )
    m->write_txt( ostr );
  
  return !ostr.bad();
}//bool write_txt( std::ostream& ostr ) const



