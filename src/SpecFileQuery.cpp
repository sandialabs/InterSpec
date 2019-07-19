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

#include <mutex>
#include <atomic>
#include <cerrno>
#include <limits>
#include <sstream>
#include <functional>

#include <boost/regex.hpp>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/SpecFileQuery.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/SpecFileQueryDbCache.h"

using namespace std;

namespace SpecFileQuery
{
  const size_t EventXmlTest::sm_min_event_xml_file_size = 128;
  const size_t EventXmlTest::sm_max_event_xml_file_size = 64*1024;
}



namespace SpecFileQuery
{
#if( BUILD_AS_OSX_APP )
  static std::atomic<bool> sm_selecting_directory( false );
  static std::mutex sm_search_directory_mutex;
  static std::string sm_search_directory;

  bool isSelectingDirectory()
  {
    return sm_selecting_directory.load();
  }
  
  void setIsSelectingDirectory( bool isDir )
  {
    sm_selecting_directory = isDir;
  }
  
  void setSearchDirectory( std::string path )
  {
    std::lock_guard<std::mutex> lock( sm_search_directory_mutex );
    sm_search_directory = path;
  }
  
  std::string getSearchDirectory()
  {
    std::lock_guard<std::mutex> lock( sm_search_directory_mutex );
    return sm_search_directory;
  }
#endif
  
    
    
    
  
  SpecTest::SpecTest()
  : m_searchField( NumFileDataFields )
  {
  }
  
  //Will throw if field is not a string search filed
  void SpecTest::set_test( const FileDataField field, const std::string &searchstr, const TextFieldSearchType type )
  {
    switch( field )
    {
      //String type searches
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
        break;
        
      case DetectionSystemType:
      case SearchMode:
      case ContainedNuetronDetector:
      case ContainedDeviationPairs:
      case HasGps:
      case HasRIIDAnalysis:
      case EnergyCalibrationType:
      case TotalLiveTime:
      case TotalRealTime:
      case IndividualSpectrumLiveTime:
      case IndividualSpectrumRealTime:
      case NumberOfSamples:
      case NumberOfRecords:
      case NumberOfGammaChannels:
      case MaximumGammaEnergy:
      case Latitude:
      case Longitude:
      case NeutronCountRate:
      case GammaCountRate:
      case StartTime:
      case NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ") + to_string(field) + " not a string field" );
        break;
    }//switch( m_searchField )
    
    m_searchField = field;
    m_searchString = searchstr;
    m_stringSearchType = type;
  }//set_test string field )
  
  void SpecTest::set_discreet_test( const FileDataField field, const int val )
  {
    switch( field )
    {
      case DetectionSystemType: case SearchMode:
      case ContainedNuetronDetector:
      case ContainedDeviationPairs:
      case HasGps:
      case HasRIIDAnalysis:
      case EnergyCalibrationType:
        break;
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case NumberOfSamples: case NumberOfRecords: case NumberOfGammaChannels:
      case MaximumGammaEnergy: case Latitude: case Longitude:
      case NeutronCountRate: case GammaCountRate:
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
      case StartTime: case NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ") + to_string(field) + " not a discrete field" );
        break;
    }//switch( m_searchField )
    
    m_searchField = field;
    m_discreteOption = val;
  }
  
  void SpecTest::set_numeric_test( const FileDataField field, const double value, const NumericFieldMatchType type )
  {
    switch( field )
    {
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case NumberOfSamples: case NumberOfRecords: case NumberOfGammaChannels:
      case MaximumGammaEnergy: case Latitude: case Longitude:
      case NeutronCountRate: case GammaCountRate:
        break;
        
      case DetectionSystemType: case SearchMode: case ContainedNuetronDetector:
      case ContainedDeviationPairs: case HasGps: case EnergyCalibrationType:
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide: case HasRIIDAnalysis:
      case StartTime: case NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ") + to_string(field) + " not a numeric field" );
        break;
    }//switch( m_searchField )
    
    m_numeric = value;
    m_searchField = field;
    m_compareType = type;
  }
  
  
  //6 year old, 
  
  void SpecTest::set_time_test( const FileDataField field, boost::posix_time::ptime comptime, const NumericFieldMatchType type )
  {
    switch( field )
    {
      case StartTime:
        break;
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case NumberOfSamples: case NumberOfRecords: case NumberOfGammaChannels:
      case MaximumGammaEnergy: case Latitude: case Longitude:
      case NeutronCountRate: case GammaCountRate:
      case DetectionSystemType: case SearchMode:
      case ContainedNuetronDetector: case ContainedDeviationPairs: case HasGps:
      case EnergyCalibrationType:
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide: case HasRIIDAnalysis:
      case NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ") + to_string(field) + " not a time field" );
        break;
    }//switch( m_searchField )
    
    
    m_searchField = field;
    m_compareType = type;
    m_time = comptime;
  }//void set_time_test(...)
  
  
  bool SpecTest::test_string( const std::string &teststr, const TextFieldSearchType &t, const std::string &ss )
  {
    if( ss.empty() )
      return true;
    
    switch( t )
    {
      case TextIsExact:
        return UtilityFunctions::iequals( teststr, ss );
        
      case TextNotEqual:
        return !UtilityFunctions::iequals( teststr, ss );
        
      case TextIsContained:
        return UtilityFunctions::icontains( teststr, ss );
        
      case TextDoesNotContain:
        return !UtilityFunctions::icontains( teststr, ss );
        
      case TextStartsWith:
        return UtilityFunctions::istarts_with( teststr, ss );
        
      case TextDoesNotStartWith:
        return !UtilityFunctions::istarts_with( teststr, ss );
        
      case TextEndsWith:
        return UtilityFunctions::iends_with( teststr, ss );
      case TextDoesNotEndWith:
        return !UtilityFunctions::iends_with( teststr, ss );
        
      case TextRegex:
      {
        //Next line can throw, and I think is doing a lot of needless work inside
        //  the loop, so should maybe be optimized some day
        boost::regex expression( ss, boost::regex::icase );
        return boost::regex_match( teststr, expression );
      }//case TextRegex:
    }//switch( t )
    
    throw runtime_error( "SpecTest::test_string Needs updating" );
    
    return false;
  }//bool SpecTest::test_string(...)
  
  
  bool SpecTest::test_date( const std::string &teststr, const NumericFieldMatchType &t, const boost::posix_time::ptime &dt )
  {
    auto testdate = UtilityFunctions::time_from_string(teststr.c_str());
    if( testdate.is_special() )
      return false;
    
    if( dt.is_special() )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION, "Was passed in a posix_time::ptime that is invalid!");
#endif
      return false;
    }
    
    switch( t )
    {
      case ValueIsExact:
        return (testdate.date() == dt.date());
      case ValueIsGreaterThan:
        return (testdate > dt);
      case ValueIsLessThan:
        return (testdate < dt);
      case ValueIsNotEqual:
        return (testdate.date() != dt.date());
    }//switch( t )
    
    return false;
  }//bool SpecTest::test_date(...)
  
  
  bool SpecTest::test( const SpecFileInfoToQuery &meas ) const
  {
    if( !meas.is_spectrum_file )
      return false;
   
    /* Consider changing so that for "not equal" cases, returning false if any
     searched field matches.  For example: if a file has remarks {"a","b"} and
     user asks for Remark TextDoesNotContain "b", then we should return false,
     where right now we are returning true.
     
    bool invert_str_search = false;
    TextFieldSearchType str_search_type = m_stringSearchType;
    if( str_search_type == TextFieldSearchType::TextNotEqual )
    {
      str_search_type = TextFieldSearchType::TextIsExact;
      invert_str_search = true;
    }else if( str_search_type == TextFieldSearchType::TextDoesNotContain )
    {
      str_search_type = TextFieldSearchType::TextIsContained;
      invert_str_search = true;
    }else if( str_search_type == TextFieldSearchType::TextDoesNotStartWith )
    {
      str_search_type = TextFieldSearchType::TextStartsWith;
      invert_str_search = true;
    }else if( str_search_type == TextFieldSearchType::TextDoesNotEndWith )
    {
      str_search_type = TextFieldSearchType::TextEndsWith;
      invert_str_search = true;
    }
     */
  
    
    switch( m_searchField )
    {
      case ParentPath:
        return test_string( UtilityFunctions::parent_path(meas.filename), m_stringSearchType, m_searchString );
        
      case Filename:
        return test_string( meas.filename, m_stringSearchType, m_searchString );
        
      case DetectorName:
      {
        bool found = false;
        for( auto it = begin(meas.detector_names); !found && (it != end(meas.detector_names)); ++it )
          found |= test_string( *it, m_stringSearchType, m_searchString );
        return found;
      }
        
      case SerialNumber:
        return test_string( meas.serial_number, m_stringSearchType, m_searchString );
        
      case Manufacturer:
        return test_string( meas.manufacturer, m_stringSearchType, m_searchString );
        
      case Model:
        return test_string( meas.model, m_stringSearchType, m_searchString );
        
      case Uuid:
        return test_string( meas.uuid, m_stringSearchType, m_searchString );
        
      case Remark:
      {
        bool found = false;
        
        for( auto it = begin(meas.file_remarks); !found && (it != end(meas.file_remarks)); ++it )
          found = test_string( *it, m_stringSearchType, m_searchString );
        
        if( found )
          return found;
        
        for( auto it = begin(meas.record_remarks); !found && (it != end(meas.record_remarks)); ++it )
          found = test_string( *it, m_stringSearchType, m_searchString );

        return found;
      }
        
      case LocationName:
        return test_string( meas.location_name, m_stringSearchType, m_searchString );
        
      case HasRIIDAnalysis:
        return (m_discreteOption==1 ? meas.has_riid_analysis : !meas.has_riid_analysis);
        
      case AnalysisResultText:
      {
        if( !meas.has_riid_analysis )
          return false;
        
        const vector<string> &remarks = meas.riid_ana.remarks_;
        const vector<DetectorAnalysisResult> &res = meas.riid_ana.results_;
        
        for( const auto &nv : meas.riid_ana.algorithm_component_versions_ )
        {
          if( test_string( nv.first, m_stringSearchType, m_searchString ) )
            return true;
          if( test_string( nv.second, m_stringSearchType, m_searchString ) )
            return true;
        }
        
        if( test_string( meas.riid_ana.algorithm_name_, m_stringSearchType, m_searchString ) )
          return true;
        
        for( size_t i = 0; i < remarks.size(); ++i )
          if( test_string( remarks[i], m_stringSearchType, m_searchString ) )
            return true;
        
        for( size_t i = 0; i < res.size(); ++i )
        {
          const DetectorAnalysisResult &r = res[i];
          if( test_string( r.detector_, m_stringSearchType, m_searchString ) )
            return true;
          if( test_string( r.id_confidence_, m_stringSearchType, m_searchString ) )
            return true;
          if( test_string( r.nuclide_, m_stringSearchType, m_searchString ) )
            return true;
          if( test_string( r.nuclide_type_, m_stringSearchType, m_searchString ) )
            return true;
          if( test_string( r.remark_, m_stringSearchType, m_searchString ) )
            return true;
        }//for( size_t i = 0; i < res.size(); ++i )
        
        return false;
      }//case AnalysisResultText:
        
        
      case AnalysisResultNuclide:
      {
        if( !meas.has_riid_analysis )
          return false;
        
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        const SandiaDecay::Nuclide * const nuc = db->nuclide( m_searchString );
      
        //const vector<string> &remarks = ana->remarks_;
        const vector<DetectorAnalysisResult> &res = meas.riid_ana.results_;
        
        if( nuc )
        {
          bool found_nuclide = false;
          
          for( size_t i = 0; !found_nuclide && i < res.size(); ++i )
          {
            //Maybe there was somethign like "Am 241" that splitting them up
            //  would mess up.
            const SandiaDecay::Nuclide *testnuc = db->nuclide( res[i].nuclide_ );
            if( nuc == testnuc )
              return true;
            
            //Just in case nuclide has mutliple fields for some reason, we'll
            //  split them up and try; this isnt foolproof,
            //  ex "found Am-241, Th 232" would fail for looking for Th232, but
            //  its better than nothing for now
            vector<string> fields;
            UtilityFunctions::split( fields, res[i].nuclide_, " \t\n,;" );
            
            
            for( size_t j = 0; !found_nuclide && j < fields.size(); ++j )
            {
              testnuc = db->nuclide( fields[j] );
              if( !testnuc && ((j+1)<fields.size()) )
                testnuc = db->nuclide( fields[j] + fields[j+1] ); //Try to make uf for things like "Th 232"
              
              found_nuclide = (nuc == testnuc);
            }
          }//for( size_t i = 0; i < res.size(); ++i )
          
          //assert( (m_stringSearchType==TextDoesNotContain) || (m_stringSearchType==TextIsContained) );
          
          return (m_stringSearchType==TextDoesNotContain) ? !found_nuclide : found_nuclide;
        }else
        {
          //Go through here and check if the specified text is contained within nuclide_, nuclide_type_, or remark_
          for( size_t i = 0; i < res.size(); ++i )
          {
            const DetectorAnalysisResult &r = res[i];
            if( test_string( r.nuclide_, m_stringSearchType, m_searchString ) )
              return true;
            if( test_string( r.nuclide_type_, m_stringSearchType, m_searchString ) )
              return true;
            if( test_string( r.remark_, m_stringSearchType, m_searchString ) )
              return true;
          }//for( size_t i = 0; i < res.size(); ++i )
        }
        
        return false;
      }//case AnalysisResultNuclide:
        
      case DetectionSystemType:
      {
        const DetectorType type = DetectorType(m_discreteOption);
        return (meas.detector_type == type);
        
        break;
      }//case DetectionSystemType:
        
      case SearchMode:
        return ((m_discreteOption==0 && !meas.passthrough) || (m_discreteOption && meas.passthrough));
        
      case TotalLiveTime:
      {
        switch( m_compareType )
        {
          case ValueIsExact:       return (fabs(meas.total_livetime - m_numeric) < 0.001);
          case ValueIsNotEqual:    return (fabs(meas.total_livetime - m_numeric) >= 0.001);
          case ValueIsLessThan:    return (meas.total_livetime < m_numeric);
          case ValueIsGreaterThan: return (meas.total_livetime > m_numeric);
        }
        return false;
      }//case TotalLiveTime:
        
      case TotalRealTime:
      {
        switch( m_compareType )
        {
          case ValueIsExact:       return (fabs(meas.total_realtime-m_numeric) < 0.001);
          case ValueIsNotEqual:    return (fabs(meas.total_realtime-m_numeric) >= 0.001);
          case ValueIsLessThan:    return (meas.total_realtime < m_numeric);
          case ValueIsGreaterThan: return (meas.total_realtime > m_numeric);
        }
        return false;
      }//case TotalRealTime:
        
      case ContainedNuetronDetector:
        return ((m_discreteOption==1) ? meas.contained_neutron : !meas.contained_neutron);
        
      case ContainedDeviationPairs:
        return ((m_discreteOption==1) ? meas.contained_dev_pairs : !meas.contained_dev_pairs);
        
      case HasGps:
        return (m_discreteOption==1 ? meas.contained_gps : !meas.contained_gps);
        
      case EnergyCalibrationType:
      {
        const auto type = Measurement::EquationType(m_discreteOption);
        for( const auto cal : meas.energy_cal_types )
          if( cal == type )
            return true;
        return false;
      }//case EnergyCalibrationType:
        
      case IndividualSpectrumLiveTime:
      {
        for( const auto lt : meas.individual_spectrum_live_time )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if( fabs(lt-m_numeric) < 0.001 ) return true; break;
            case ValueIsNotEqual:    if( fabs(lt-m_numeric) > 0.001 ) return true; break;
            case ValueIsLessThan:    if(lt < m_numeric) return true; break;
            case ValueIsGreaterThan: if(lt > m_numeric) return true; break;
          }
        }
        return false;
      }//case IndividualSpectrumLiveTime:
        
      case IndividualSpectrumRealTime:
      {
        for( const auto rt : meas.individual_spectrum_real_time )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if( fabs(rt-m_numeric) < 0.001 ) return true; break;
            case ValueIsNotEqual:    if( fabs(rt-m_numeric) > 0.001 ) return true; break;
            case ValueIsLessThan:    if(rt < m_numeric) return true; break;
            case ValueIsGreaterThan: if(rt > m_numeric) return true; break;
          }
        }
        return false;
      }//case IndividualSpectrumRealTime:
        
      case NumberOfSamples:
      {
        //The number of samples we will ever deal with will always be able to
        //  be represented exactly as a double.
        const double nsample = static_cast<double>( meas.number_of_samples );
        
        switch( m_compareType )
        {
          case ValueIsExact:       if(nsample == m_numeric) return true; break;
          case ValueIsNotEqual:    if(nsample != m_numeric) return true; break;
          case ValueIsLessThan:    if(nsample < m_numeric) return true;  break;
          case ValueIsGreaterThan: if(nsample > m_numeric) return true;  break;
        }
        return false;
      }//case NumberOfSamples:
        
      case NumberOfRecords:
      {
        //The number of measurements we will ever deal with will always be able
        //  to be represented exactly as a double.
        const double nsample = static_cast<double>( meas.number_of_records );
        
        switch( m_compareType )
        {
          case ValueIsExact:       if(nsample == m_numeric) return true; break;
          case ValueIsNotEqual:    if(nsample != m_numeric) return true; break;
          case ValueIsLessThan:    if(nsample < m_numeric) return true;  break;
          case ValueIsGreaterThan: if(nsample > m_numeric) return true;  break;
        }
        return false;
      }//case NumberOfRecords:
        
      case NumberOfGammaChannels:
      {
        for( const auto n : meas.number_of_gamma_channels )
        {
          if( !n )
            continue;
          
          const double nchann = static_cast<double>( n );
          
          switch( m_compareType )
          {
            case ValueIsExact:       if(nchann == m_numeric) return true; break;
            case ValueIsNotEqual:    if(nchann != m_numeric) return true; break;
            case ValueIsLessThan:    if(nchann < m_numeric) return true;  break;
            case ValueIsGreaterThan: if(nchann > m_numeric) return true;  break;
          }
        }//for( const auto &m : meas->measurements() )
        
        return false;
      }//case NumberOfGammaChannels:
        
      case MaximumGammaEnergy:
      {
        for( const auto energy : meas.max_gamma_energy )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if( fabs(energy-m_numeric) < 0.1 ) return true; break;
            case ValueIsNotEqual:    if( fabs(energy-m_numeric) > 0.1 ) return true; break;
            case ValueIsLessThan:    if(energy < m_numeric) return true; break;
            case ValueIsGreaterThan: if(energy > m_numeric) return true; break;
          }
        }//for( const auto &m : meas->measurements() )
        
        return false;
      }//case MaximumGammaEnergy:
        
      case Latitude:
      {
        if( !meas.contained_gps )
          return false;
        
        switch( m_compareType )
        {
          case ValueIsExact:       return (fabs(meas.mean_latitude-m_numeric) < 0.000001);
          case ValueIsNotEqual:    return (fabs(meas.mean_latitude-m_numeric) > 0.000001);
          case ValueIsLessThan:    return (meas.mean_latitude < m_numeric);
          case ValueIsGreaterThan: return (meas.mean_latitude > m_numeric);
        }
      }//case Latitude:
        
      case Longitude:
      {
        if( !meas.contained_gps )
          return false;
        
        switch( m_compareType )
        {
          case ValueIsExact:       return (fabs(meas.mean_longitude - m_numeric) < 0.000001);
          case ValueIsNotEqual:    return (fabs(meas.mean_longitude - m_numeric) > 0.000001);
          case ValueIsLessThan:    return (meas.mean_longitude < m_numeric);
          case ValueIsGreaterThan: return (meas.mean_longitude > m_numeric);
        }
      }//case Longitude:
        
      case NeutronCountRate:
      {
        for( const auto cps : meas.neutron_count_rate )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if( fabs(cps-m_numeric) < 1.0E-6 ) return true; break;
            case ValueIsNotEqual:    if( fabs(cps-m_numeric) > 1.0E-6 ) return true; break;
            case ValueIsLessThan:    if(cps < m_numeric) return true; break;
            case ValueIsGreaterThan: if(cps > m_numeric) return true; break;
          }
        }//for( const auto &m : meas->measurements() )
        
        return false;
      }//case NeutronCountRate:
        
      case GammaCountRate:
      {
        for( const auto cps : meas.gamma_count_rate )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if( fabs(cps-m_numeric) < 1.0E-6 ) return true; break;
            case ValueIsNotEqual:    if( fabs(cps-m_numeric) > 1.0E-6 ) return true; break;
            case ValueIsLessThan:    if(cps < m_numeric) return true; break;
            case ValueIsGreaterThan: if(cps > m_numeric) return true; break;
          }
        }//for( const auto &m : meas->measurements() )
        
        return false;
      }//case GammaCountRate:
        
        
      case StartTime:
      {
        const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
        const boost::posix_time::time_duration::sec_type test_time = (m_time - epoch).total_seconds();
        
        for( const auto starttime : meas.start_times )
        {
          switch( m_compareType )
          {
            case ValueIsExact:
            {
              if( abs(starttime - test_time) < 60 ) //note: we only have minute resoltion on the gui selector
                return true;
              break;
            }
              
            case ValueIsNotEqual:
            {
              if( abs(starttime - test_time) > 60 ) //We only have minute resoltion on the gui selector
                return true;
              break;
            }
            case ValueIsLessThan:    if(starttime < test_time) return true; break;
            case ValueIsGreaterThan: if(starttime > test_time) return true; break;
          }
        }
        return false;
      }//case StartTime:
        
      case NumFileDataFields:
        return false;
        break;
    }//switch( m_searchField )
    
    return false;
  }//bool test( const SpecFileInfoToQuery &meas ) const
  
  
  
  const char *to_string( const FileDataField field )
  {
    /* These values must exactly match those the 'id' field of the filters in
       SpecFileQueryWidget.js.
     */
    switch( field )
    {
      case ParentPath:                 return "ParentPath";
      case Filename:                   return "Filename";
      case DetectorName:               return "Detector name";
      case SerialNumber:               return "Serial number";
      case Manufacturer:               return "Manufacturer";
      case Model:                      return "Model";
      case Uuid:                       return "UUID";
      case Remark:                     return "Remark";
      case LocationName:               return "Location name";
      case AnalysisResultText:         return "RIID Ana result";
      case AnalysisResultNuclide:      return "RIID IDed nuclide";
      case HasRIIDAnalysis:            return "Has RIID Analysis";
      case DetectionSystemType:        return "Detector system";
      case ContainedNuetronDetector:   return "Has Neutron";
      case ContainedDeviationPairs:    return "Has Dev. Pairs";
      case HasGps:                     return "Has GPS Info";
      case EnergyCalibrationType:      return "Energy Cal Type";
      case SearchMode:                 return "Aquisition mode";
      case TotalLiveTime:              return "Sum live time";
      case TotalRealTime:              return "Sum real time";
      case IndividualSpectrumLiveTime: return "Spectrum live time";
      case IndividualSpectrumRealTime: return "Spectrum real time";
      case NumberOfSamples:            return "Num. Time Samples";
      case NumberOfRecords:            return "Num. Records";
      case NumberOfGammaChannels:      return "Num Gamma Channels";
      case MaximumGammaEnergy:         return "Max Gamma Energy";
      case Latitude:                   return "Latitude";
      case Longitude:                  return "Longitude";
      case NeutronCountRate:           return "Neutron CPS";
      case GammaCountRate:             return "Gamma CPS";
      case StartTime:                  return "Start Time";
        
      case NumFileDataFields:
        break;
    }//switch( m_searchField )
    
    return "";
  }//to_string( const FileDataField field )
  

  bool from_string( const std::string &fieldName, FileDataField &field )
  {
    for( FileDataField i = FileDataField(0); i < NumFileDataFields; i = FileDataField(i+1) )
    {
      if( fieldName == to_string(i) )
      {
        field = i;
        return true;
      }
    }
    
    field = NumFileDataFields;
    return false;
  }//bool from_string( const std::string &fieldName, FileDataField &field )

  
  
  const char *to_string( const TextFieldSearchType type )
  {
    switch( type )
    {
      case TextIsExact:          return "equal";
      case TextNotEqual:         return "not equal";
      case TextIsContained:      return "contains";
      case TextDoesNotContain:   return "does not contain";
      case TextStartsWith:       return "begins with";
      case TextDoesNotStartWith: return "does not begin with";
      case TextEndsWith:         return "ends with";
      case TextDoesNotEndWith:   return "does not end with";
      case TextRegex:            return "regex";
        break;
    }
    
    return "";
  }//const char *to_string( const TextFieldSearchType type );
  
  
  bool from_string( const std::string &val, TextFieldSearchType &type )
  {
    const TextFieldSearchType types[] = { TextIsExact, TextNotEqual,
      TextIsContained, TextDoesNotContain, TextStartsWith, TextDoesNotStartWith,
      TextEndsWith, TextDoesNotEndWith, TextRegex
    };
    
    for( const TextFieldSearchType &t : types )
    {
      if( val == to_string(t) )
      {
        type = t;
        return true;
      }
    }
    
    return false;
  }//bool from_string( const std::string &val, TextFieldSearchType &type )
  
  
  const char *to_string( const LogicType type )
  {
    switch( type )
    {
      case LogicalOr:         return "OR";
      case LogicalAnd:        return "AND";
      case LogicalNot:        return "NOT";
      case LogicalOpenParan:  return "Open Parenthesis";
      case LogicalCloseParan: return "close Parenthesis";
      case NumLogicType:      break;
    }
    return "";
  }
  
  const char *to_string( const NumericFieldMatchType type )
  {
    switch( type )
    {
      case ValueIsExact:       return "equal";
      case ValueIsLessThan:    return "less than";
      case ValueIsGreaterThan: return "greater than";
      case ValueIsNotEqual:    return "not equal";
    }
    return "";
  }
  
  bool from_string( const std::string &val, NumericFieldMatchType &type )
  {
    if( val == "equal" )
    {
      type = NumericFieldMatchType::ValueIsExact;
      return true;
    }else if( val == "less than" || val == "less" )
    {
      type = NumericFieldMatchType::ValueIsLessThan;
      return true;
    }else if( val == "greater than" || val == "greater" )
    {
      type = NumericFieldMatchType::ValueIsGreaterThan;
      return true;
    }else if( val == "not equal" || val == "not_equal" )
    {
      type = NumericFieldMatchType::ValueIsNotEqual;
      return true;
    }

    return false;
  }//bool from_string( const std::string &val, NumericFieldMatchType &type )
  
  
  std::string SpecTest::summary() const
  {
    string summary = to_string( m_searchField );
    
    switch( m_searchField )
    {
        //String type searches
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
        summary += string(" ") + to_string( m_stringSearchType )
        + string(" \"") + m_searchString + "\"";
        break;
        
      case HasRIIDAnalysis:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have RIID Analysis Results";
        break;
        
      case DetectionSystemType:
        summary += " is a " + detectorTypeToString( DetectorType(m_discreteOption) );
        break;
        
      case SearchMode:
        summary += " is " + string(m_discreteOption ? "search/passthrough" : "dwell");
        break;
        
      case ContainedNuetronDetector:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have a neutron detectoor";
        break;
        
      case ContainedDeviationPairs:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have non-linear deviation pairs";
        break;
        
      case HasGps:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have GPS information";
        break;
        
      case EnergyCalibrationType:
        switch( Measurement::EquationType(m_discreteOption) )
        {
          case Measurement::EquationType::Polynomial: summary += "polynomial"; break;
          case Measurement::EquationType::FullRangeFraction: summary += "full range fraction"; break;
          case Measurement::EquationType::LowerChannelEdge: summary += "lower channel energy"; break;
          case Measurement::EquationType::InvalidEquationType:
          case Measurement::EquationType::UnspecifiedUsingDefaultPolynomial:
            summary += "unknown";
            break;
        }
        summary += " energy calibration";
        break;
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case NumberOfSamples: case NumberOfRecords: case NumberOfGammaChannels:
      case MaximumGammaEnergy: case Latitude: case Longitude:
      case NeutronCountRate: case GammaCountRate:
        summary += " is " + string(to_string(m_compareType)) + " " + PhysicalUnits::printToBestTimeUnits(m_numeric);
        break;
        
      case StartTime:
        summary += " is " + string(to_string(m_compareType)) + " " + UtilityFunctions::to_iso_string(m_time);
        break;
        
      case NumFileDataFields:
        break;
    }//switch( m_searchField )
    
    return summary;
  }//string summary()
  
  
  void SpecTest::isvalid()
  {
    switch( m_searchField )
    {
      //String type searches
      case ParentPath: case Filename: case DetectorName: case SerialNumber:
      case Manufacturer: case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
      {
        switch( m_stringSearchType )
        {
          case TextNotEqual:
          case TextIsExact:
            break;
            
          case TextIsContained:
          case TextDoesNotContain:
          case TextStartsWith:
          case TextDoesNotStartWith:
          case TextEndsWith:
          case TextDoesNotEndWith:
            if( m_searchString.empty() )
              throw runtime_error( string("Search string for ") + to_string(m_stringSearchType) + " can not be empty" );
            break;
            
          case TextRegex:
            if( m_searchString.empty() )
              throw runtime_error( string("Search string for ") + to_string(m_stringSearchType) + " can not be empty" );
            
            try
            {
              boost::regex expression( m_searchString, boost::regex::icase );
              boost::regex_match( "Somedummystr", expression );
            }catch( std::exception &e )
            {
              cout << e.what() << endl;
              throw runtime_error( "\"" + m_searchString + "\" is an invalid regular expression" );
            }
            
            break;
        }//switch ( m_stringSearchType)
        
        break;
      }//case All string searches
        
      case DetectionSystemType:
      {
        //We could check m_discreteOption is
        bool gotit = false;
        switch( DetectorType(m_discreteOption) )
        {
          case kGR135Detector: case kIdentiFinderDetector:
          case kIdentiFinderNGDetector: case kIdentiFinderLaBr3Detector:
          case kDetectiveDetector: case kDetectiveExDetector:
          case kDetectiveEx100Detector: case kDetectiveEx200Detector:
          case kDetectiveX:
          case kSAIC8Detector: case kFalcon5000: case kUnknownDetector:
          case kMicroDetectiveDetector: case kMicroRaiderDetector:
          case kRadHunterNaI: case kRadHunterLaBr3: case kRsi701: case kRsi705:
          case kAvidRsi:
          case kOrtecRadEagleNai: case kOrtecRadEagleCeBr2Inch:
          case kOrtecRadEagleCeBr3Inch: case kOrtecRadEagleLaBr:
          case kSam940LaBr3: case kSam940: case kSam945: case kSrpm210:
            gotit = true;
            break;
        }
        
        if( !gotit )
          throw runtime_error( "Invalid DetectorType (programming error)" );
        break;
      }//case DetectionSystemType:
        
      case SearchMode:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid aquitision mode (programming error)" );
        break;
      }
      
      case HasRIIDAnalysis:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid RIID Analysis present option (programming error)" );
        break;
      }
        
      case ContainedNuetronDetector:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid neutron detector presence (programming error)" );
        break;
      }
        
      case ContainedDeviationPairs:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid deviation pair presence (programming error)" );
        break;
      }
        
      case HasGps:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid GPS info presence (programming error)" );
        break;
      }
        
      case EnergyCalibrationType:
      {
        switch( m_discreteOption )
        {
          case Measurement::EquationType::Polynomial:
          case Measurement::EquationType::FullRangeFraction:
          case Measurement::EquationType::LowerChannelEdge:
          case Measurement::EquationType::UnspecifiedUsingDefaultPolynomial:
          case Measurement::EquationType::InvalidEquationType:
            break;
          default:
            throw runtime_error( "Invalid energy calibration type (programming error)" );
        }
      }
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case NumberOfSamples: case NumberOfRecords: case NumberOfGammaChannels:
      case MaximumGammaEnergy: case Latitude: case Longitude:
      case NeutronCountRate: case GammaCountRate:
      {
        if( IsInf(m_numeric) || IsNan(m_numeric) )
          throw runtime_error( "Invalid time duration" );
        break;
      }
        
      case StartTime:
      {
        if( m_time.is_special() )
          throw runtime_error( "Invalid Time" );
      }
        
      case NumFileDataFields:
        break;
    }//switch( m_searchField )
  }//void isvalid()
  
  
 
  
  EventXmlTest::EventXmlTest()
  {
    m_testType = TestType::NotSet;
  }//
  
  
  void EventXmlTest::set_string_test_info( const std::string &test_label,
                                           const std::string &test_string,
                                           const TextFieldSearchType &comparison_type )
  {
    m_test_string = test_string;
    m_test_label = test_label;
    m_fieldTestType = comparison_type;
    m_testType = TestType::String;
    m_test_time = boost::posix_time::ptime();
  }//set_string_test_info(...)
  
  
  void EventXmlTest::set_date_test_info( const std::string &test_label,
                                        const boost::posix_time::ptime &test_time,
                                        const NumericFieldMatchType &comparison_type )
  {
    m_testType = TestType::Date;
    m_test_label = test_label;
    m_test_time = test_time;
    m_dateTestType = comparison_type;
    m_test_string = "";
    if( m_test_time.is_special() )
      throw runtime_error( "Time to test against not a valid" );
  }//set_date_test_info(...)

  
  bool EventXmlTest::test( const SpecFileInfoToQuery &meas ) const
  {
    const auto iter = meas.event_xml_filter_values.find(m_test_label);
    if( iter == end(meas.event_xml_filter_values) )
      return false;
    
    const vector<string> &field_values = iter->second;
    
    
    for( const string &strval : field_values )
    {
      switch( m_testType )
      {
        case TestType::Date:
          if( SpecTest::test_date( strval, m_dateTestType, m_test_time ) )
            return true;
          break;
          
        case TestType::String:
          if( SpecTest::test_string( strval, m_fieldTestType, m_test_string ) )
            return true;
          break;
          
        case TestType::NotSet:
          throw runtime_error( "TestType not set!" );
          break;
      }//switch( m_testType )
    }//for( const string &val : field_values )
    
    return false;
  }//EventXmlTest::test(...)
  
  
  
  std::string EventXmlTest::summary() const
  {
    switch( m_testType )
    {
      case TestType::Date:
        return "Event XML with comparing '" + m_test_label + "' " + to_string(m_dateTestType) + " '" + UtilityFunctions::to_common_string(m_test_time,true) + "'";
        break;
        
      case TestType::String:
        return "Event XML with comparing '" + m_test_label + "' " + to_string(m_fieldTestType) + " '" + m_test_string + "'";
        break;
        
      case TestType::NotSet:
        return "Event XML with comparing '" + m_test_label + "' is an undefined test";
        break;
    }//switch( m_testType )
    
    return "Error";
  }//std::string summary() const
  
  
  void EventXmlTest::isvalid()
  {
    if( m_test_label.empty() )
      throw runtime_error( "Test label is empty" );
    
    switch( m_testType )
    {
      case TestType::NotSet:
        throw runtime_error( "Test not set." );
        break;
        
      case TestType::Date:
      {
        if( m_test_time.is_special() )
          throw runtime_error( "Test datetime ('" + m_test_string + "') is invalid" );
        break;
      }
        
      case TestType::String:
      {
        if( m_test_string.empty() )
          throw runtime_error( "String to test against is empty." );
        
        bool validTestField = false;
        switch( m_fieldTestType )
        {
          case TextIsExact: case TextNotEqual: case TextIsContained:
          case TextDoesNotContain: case TextStartsWith:
          case TextDoesNotStartWith: case TextEndsWith: case TextDoesNotEndWith:
          case TextRegex:
            validTestField = true;
            break;
        }//switch( m_fieldTestType )
    
        if( !validTestField )
          throw runtime_error( "String test type somehow isnt valid." );
        
        break;
      }
    }//switch( m_testType )
  }//void isvalid()
  
  
  SpecLogicTest::SpecLogicTest()
  {
  }
  
  void SpecLogicTest::addCondition( const SpecTest &test )
  {
    m_fields.push_back( boost::any(test) );
  }
  
  void SpecLogicTest::addLogic( const LogicType test )
  {
    m_fields.push_back( boost::any(test) );
  }
  
  void SpecLogicTest::addEventXmlTest( const EventXmlTest &test )
  {
    m_fields.push_back( boost::any(test) );
  }
  
  std::string SpecLogicTest::summary() const
  {
    stringstream strm;
    print_equation( m_fields, strm );
    return strm.str();
  }
  
  std::ostream &SpecLogicTest::print_equation( vector<boost::any> fields, std::ostream &strm )
  {
    for( size_t i = 0; i < fields.size(); ++i )
    {
      try
      {
        const LogicType val = boost::any_cast<LogicType>( fields[i] );
        switch( val )
        {
          case LogicalOr:         strm << " || "; break;
          case LogicalAnd:        strm << " && "; break;
          case LogicalNot:        strm << "!";    break;
          case LogicalOpenParan:  strm << "(";    break;
          case LogicalCloseParan: strm << ")";    break;
          case NumLogicType: break;
        }
      }catch(...){}
      
      try
      {
        const bool val = boost::any_cast<bool>( fields[i] );
        strm << (val ? "True" : "False");
      }catch(...){}
      
      try
      {
        const SpecTest test = boost::any_cast<SpecTest>( fields[i] );
        strm << "[" << test.summary() << "]";
      }catch(...){}
      
      try
      {
        const EventXmlTest test = boost::any_cast<EventXmlTest>( fields[i] );
        strm << "[" << test.summary() << "]";
      }catch(...){}
      
    }//for( size_t i = 0; i < fields.size(); ++i )
    
    return strm;
  }//std::ostream &print_equation(...)
  
  
  bool SpecLogicTest::evaluate( std::vector<boost::any> fields, const SpecFileInfoToQuery &meas )
  {
    if( !meas.is_file )
      return false;
    
    if( fields.empty() )
      return true;
    
    //cerr << "Initial test:\n\t";
    //print_equation( fields, cerr ) << endl;
    
    //find open and close parenthesis, and evaluate them recursively
    for( size_t i = 0; i < fields.size(); ++i )
    {
      try
      {
        const LogicType startlogic = boost::any_cast<LogicType>( fields[i] );
        if( startlogic == LogicalOpenParan )
        {
          int nparen = 1;
          size_t closepos = i + 1;
          
          for( ; closepos < fields.size(); ++closepos )
          {
            try
            {
              const LogicType closelogic = boost::any_cast<LogicType>( fields[closepos] );
              if( closelogic == LogicalOpenParan )
                ++nparen;
              else if( closelogic == LogicalCloseParan )
                --nparen;
            }catch( ... ){}
            
            if( nparen == 0 )
              break;
          }//for( ; closepos < fields.size(); ++closepos )
          
          if( closepos >= fields.size() )
            throw runtime_error( "Failed to find closing parenthesis, invalid expression" );
          
          vector<boost::any> inside( fields.begin() + i + 1, fields.begin() + closepos );
          
          //cerr << "Had " << fields.size() << " elements\n\t";
          //print_equation( fields, cerr ) << endl;
          //cerr << "Evaluating " << inside.size() << " elements\n\t";
          //print_equation( inside, cerr ) << endl;
          
          
          fields.erase( fields.begin() + i, fields.begin() + closepos + 1 );
          //cerr << "Reduced down to " << fields.size() << " elements (will add one)\n\t";
          //print_equation( fields, cerr ) << endl;
          
          const bool answer = SpecLogicTest::evaluate( inside, meas );
          
          fields.insert( fields.begin() + i, boost::any(answer) );
          //cerr << "New eqn becomes:" << endl;
          //print_equation( fields, cerr ) << endl;
        }//if( startlogic == LogicalOpenParan )
        
        continue;
      }catch( ... ){}
      
      try
      {
        //At some point in the future should imlement doing things a little
        //  better so we can short-circuit some logic.  But since parsing the
        //  file is the lions share of the time for this test, this optimization
        //  can be put off.
        const SpecTest test = boost::any_cast<SpecTest>( fields[i] );
        const bool pass = test.test( meas );
        fields[i] = boost::any( pass );
        continue;
      }catch(... ){}
      
      
      try
      {
        const EventXmlTest test = boost::any_cast<EventXmlTest>( fields[i] );
        const bool pass = test.test( meas );
        fields[i] = boost::any( pass );
        continue;
      }catch(...){}
      
      //hmmm... should we ever get here?
    }//for( size_t i = 0; i < fields.size(); ++i )
    
    
    //cerr << "After getting rid of parenthesis and evaluating terms this eqn becomes:\n\t";
    //print_equation( fields, cerr ) << endl;
    
    //go through and toggle flip NOTs
    for( size_t i = 0; i < fields.size(); ++i )
    {
      try
      {
        const LogicType notlogic = boost::any_cast<LogicType>( fields[i] );
        if( notlogic == LogicalNot )
        {
          if( i == (fields.size()-1) )
            throw runtime_error( "! encountered at end of equation" );
          
          try
          {
            const bool nextval = boost::any_cast<bool>( fields[i+1] );
            fields[i+1] = boost::any( !nextval );
            fields.erase( fields.begin() + i );
          }catch(...)
          {
            throw runtime_error( "Expected boolean following !" );
          }
        }
      }catch(...){}
    }//for( size_t i = 0; i < fields.size(); ++i )
    
    //    cerr << "After getting rid of NOTS, eqn becomes:\n\t";
    //    print_equation( fields, cerr ) << endl;
    
    
    if( (fields.size()%2) == 0 )
      throw runtime_error( "Expect there to be an odd number of elements at this point" );
    
    for( size_t i = 0; i < fields.size(); ++i )
    {
      if( (i%2) == 0 )
      {
        try
        {
          boost::any_cast<bool>(fields[i]);
        }catch(...)
        {
          throw runtime_error( "Expect all boost::any's to be bools in even locations at this point" );
        }
      }else
      {
        try
        {
          const LogicType logic = boost::any_cast<LogicType>( fields[i] );
          if( logic != LogicalOr && logic != LogicalAnd )
            throw runtime_error( "Expect all logic to be AND or OR only at this point" );
        }catch(...)
        {
          throw runtime_error( "Expect all boost::any's to be LogicType's in odd locations at this point" );
        }
      }
    }//for( size_t i = 0; i < fields.size(); ++i )
    
    bool answer = boost::any_cast<bool>( fields[0] );
    
    for( size_t i = 1; i < (fields.size()-1); i += 2 )
    {
      const LogicType logic = boost::any_cast<LogicType>( fields[i] );
      const bool nextval = boost::any_cast<bool>( fields[i+1] );
      if( logic == LogicalOr )
        answer = (answer || nextval);
      else if( logic == LogicalAnd )
        answer = (answer && nextval);
    }
    
    //cerr << "Eqn evaluated to " << (answer ? "True" : "False") << endl << endl << endl;
    
    return answer;
  }//bool SpecLogicTest::evaluate( std::vector<boost::any> fields, const SpecFileInfoToQuery &meas )
  
  
  bool SpecLogicTest::test( const SpecFileInfoToQuery &meas ) const
  {
    return evaluate( m_fields, meas );
  }
  
  void SpecLogicTest::isvalid()
  {
    if( m_fields.empty() )
      return;
    
    //1) Same number of ( as ), and they must match
    //2) !, ||, and && must be followed by ( or a condition
    //3) except for last condition, all conditions and ) must be followed by || && or )
    //4) All regular exression matches should have valid regex
    //5) And && followed by a condition then an || is ambiguouse, and visa versa
    
    {//1) Check matching parenthesis
      int nparan = 0;
      for( size_t i = 0; nparan >= 0 && i < m_fields.size(); ++i )
      {
        try
        {
          const LogicType t = boost::any_cast<LogicType>( m_fields[i] );
          if( t == LogicalOpenParan )
            ++nparan;
          else if( t == LogicalCloseParan )
            --nparan;
        }catch( ... ){}
      }//for( size_t i = 0; i < m_fields.size(); ++i )
      
      if( nparan < 0 )
        throw runtime_error( "Closing parenthesis doesnt match" );
      else if( nparan > 0 )
        throw runtime_error( "Missing closing parenthesis" );
    }//end 1
    
    
    {//2 !, ||, or && must be followed by ( or a condition
      for( size_t i = 0; i < m_fields.size(); ++i )
      {
        LogicType t;
        try
        {
          t = boost::any_cast<LogicType>( m_fields[i] );
        }catch(...)
        {
          continue;
        }
        
        if( t == LogicalOr || t == LogicalAnd || t == LogicalNot )
        {
          string oper;
          if( t == LogicalOr )
            oper = "OR";
          else if( t == LogicalAnd )
            oper = "AND";
          else if( t == LogicalNot )
            oper = "NOT";
          
          if( (i+1) == m_fields.size() )
            throw runtime_error( oper + " (last item) is dangling" );
          
          bool nextislogic = false;
          LogicType nextt;
          
          try
          {
            nextt = boost::any_cast<LogicType>( m_fields[i+1] );
            nextislogic = true;
          }catch(...){ }
          
          /*
           bool nextistest = false;
           SpecTest nexttest;
           try
           {
           nexttest = boost::any_cast<SpecTest>( m_fields[i+1] );
           nextistest = true;
           }catch(...){}
           */
          
          if( nextislogic && nextt != LogicalOpenParan && nextt != LogicalNot )
            throw runtime_error( oper + " must be followed by a test, open parenthesis or a NOT" );
        }
      }//for( size_t i = 0; i < m_fields.size(); ++i )
    }//end 2
    
    
    {//BEGIN 3) except for last condition, all conditions must be followed by || && or )
      for( size_t i = 0; i < m_fields.size(); ++i )
      {
        bool istest = false, isxmltest = false;
        SpecTest test;
        EventXmlTest xmltest;
        LogicType closeparen;
        try
        {
          try
          {
          test = boost::any_cast<SpecTest>( m_fields[i] );
          istest = true;
          }catch(...){}
          
          try
          {
            xmltest = boost::any_cast<EventXmlTest>( m_fields[i] );
            isxmltest = true;
          }catch(...){}
          
          if( !istest && !isxmltest )
            throw runtime_error("");
        }catch(...)
        {
          try
          {
            closeparen = boost::any_cast<LogicType>( m_fields[i] );
            if( closeparen != LogicalCloseParan )
              continue;
          }catch(...)
          {
            continue;
          }
        }
        
        if( (i+1) < m_fields.size() )
        {
          try
          {
            const LogicType nextt = boost::any_cast<LogicType>( m_fields[i+1] );
            
            switch( nextt )
            {
              case LogicalOr: case LogicalAnd: case LogicalCloseParan:
                break;
              case LogicalNot: case LogicalOpenParan: case NumLogicType:
                if( istest )
                  throw runtime_error( test.summary() + " should not be followed " + to_string(nextt) );
                else if( isxmltest )
                  throw runtime_error( xmltest.summary() + " should not be followed " + to_string(nextt) );
                else
                  throw runtime_error( "Closing parenthesis should not be followed " + string(to_string(nextt)) );
                break;
            }//switch( nextt )
          }catch(...)
          {
            if( istest )
              throw runtime_error( test.summary() + " should be followed by an operation" );
            else if( isxmltest )
              throw runtime_error( xmltest.summary() + " should be followed by an operation" );
            else
              throw runtime_error( "Closing parenthesis should be followed by an operation" );
          }
        }//if( not the last condition )
      }//for( size_t i = 0; i < m_fields.size(); ++i )
    }//END 3) except for last condition, all conditions must be followed by || && or )
    
    {//BEGIN 4) All regular exression matches should have valid regex
      for( size_t i = 0; i < m_fields.size(); ++i )
      {
        bool istest = false, isxmltest = false;
        SpecTest test;
        EventXmlTest xmltest;
        try
        {
          test = boost::any_cast<SpecTest>( m_fields[i] );
          istest = true;
        }catch(...){}
        
        try
        {
          xmltest = boost::any_cast<EventXmlTest>( m_fields[i] );
          isxmltest = false;
        }catch(...){}
        
        if( istest )
          test.isvalid();
        if( isxmltest )
          xmltest.isvalid();
      }//for( size_t i = 0; i < m_fields.size(); ++i )
    }//END 4) All regular exression matches should have valid regex
    
    {//BEGIN 5) And && followed by a condition then an || is ambiguouse, and visa versa
      //cout << "Need to enforce no ambigous order of evaluation of ORs and ANDs" << endl;
    }//END 5) And && followed by a condition then an || is ambiguouse, and visa versa
    
  }//void isvalid()
}//namespace SpecFileQuery

