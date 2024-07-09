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

#include "3rdparty/date/include/date/date.h"

#include "InterSpec/SpecMeas.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/SpecFileQuery.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/SpecFileQueryDbCache.h"

using namespace std;


// We *could* include InterSpecUser.h to get the to_ptime and to_time_point functions, or we
//  could move these function to a common header.... but instead copy/pasting, for the moment.
namespace
{
boost::posix_time::ptime to_ptime( const std::chrono::system_clock::time_point &rhs )
{
  auto dp = ::date::floor<::date::days>(rhs);
  auto ymd = ::date::year_month_day{dp};
  auto time = ::date::make_time(rhs - dp);
  
  boost::posix_time::time_duration td( time.hours().count(),
                                      time.minutes().count(),
                                      time.seconds().count(),
                                      date::floor<chrono::microseconds>(time.subseconds()).count()
                                      );
  
  const unsigned month_num = static_cast<unsigned>(ymd.month());
  boost::gregorian::greg_month month = boost::date_time::months_of_year(month_num);
  
  const boost::gregorian::date d( static_cast<int>(ymd.year()), month, (unsigned)ymd.day() );
  
  return boost::posix_time::ptime( d, td );
}//to_ptime(...)

std::chrono::system_clock::time_point to_time_point( const boost::posix_time::ptime &rhs )
{
  if( rhs.is_special() )
    return {};
  
  unsigned short year = static_cast<unsigned short>( rhs.date().year() );
  const unsigned month = rhs.date().month().as_number();
  const unsigned day = rhs.date().day().as_number();
  
  const int64_t nmicro = rhs.time_of_day().total_microseconds();
  date::year_month_day ymd{ date::year(year), date::month(month), date::day(day) };
  
  date::sys_days days = ymd;
  std::chrono::system_clock::time_point tp = days;
  tp += chrono::microseconds(nmicro);
  
  return tp;
}
}//namespace

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
      case FileDataField::ParentPath:
      case FileDataField::Filename:
      case FileDataField::DetectorName:
      case FileDataField::SerialNumber:
      case FileDataField::Manufacturer:
      case FileDataField::Model:
      case FileDataField::Uuid:
      case FileDataField::Remark:
      case FileDataField::LocationName:
      case FileDataField::AnalysisResultText:
      case FileDataField::AnalysisResultNuclide:
        break;
        
      case FileDataField::DetectionSystemType:
      case FileDataField::SearchMode:
      case FileDataField::ContainedNuetronDetector:
      case FileDataField::ContainedDeviationPairs:
      case FileDataField::HasGps:
      case FileDataField::HasRIIDAnalysis:
      case FileDataField::EnergyCalibrationType:
      case FileDataField::TotalLiveTime:
      case FileDataField::TotalRealTime:
      case FileDataField::IndividualSpectrumLiveTime:
      case FileDataField::IndividualSpectrumRealTime:
      case FileDataField::NumberOfSamples:
      case FileDataField::NumberOfRecords:
      case FileDataField::NumberOfGammaChannels:
      case FileDataField::MaximumGammaEnergy:
      case FileDataField::Latitude:
      case FileDataField::Longitude:
      case FileDataField::NeutronCountRate:
      case FileDataField::GammaCountRate:
      case FileDataField::StartTimeIoI:
      case FileDataField::MeasurementsStartTimes:
      case FileDataField::NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ")
                             + to_string(field) + " not a string field" );
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
      case FileDataField::DetectionSystemType:
      case FileDataField::SearchMode:
      case FileDataField::ContainedNuetronDetector:
      case FileDataField::ContainedDeviationPairs:
      case FileDataField::HasGps:
      case FileDataField::HasRIIDAnalysis:
      case FileDataField::EnergyCalibrationType:
        break;
        
      case FileDataField::TotalLiveTime:
      case FileDataField::TotalRealTime:
      case FileDataField::IndividualSpectrumLiveTime:
      case FileDataField::IndividualSpectrumRealTime:
      case FileDataField::NumberOfSamples:
      case FileDataField::NumberOfRecords:
      case FileDataField::NumberOfGammaChannels:
      case FileDataField::MaximumGammaEnergy:
      case FileDataField::Latitude:
      case FileDataField::Longitude:
      case FileDataField::NeutronCountRate:
      case FileDataField::GammaCountRate:
      case FileDataField::ParentPath:
      case FileDataField::Filename:
      case FileDataField::DetectorName:
      case FileDataField::SerialNumber:
      case FileDataField::Manufacturer:
      case FileDataField::Model:
      case FileDataField::Uuid:
      case FileDataField::Remark:
      case FileDataField::LocationName:
      case FileDataField::AnalysisResultText:
      case FileDataField::AnalysisResultNuclide:
      case FileDataField::StartTimeIoI:
      case FileDataField::MeasurementsStartTimes:
      case FileDataField::NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ")
                            + to_string(field) + " not a discrete field" );
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
      case FileDataField::StartTimeIoI: case MeasurementsStartTimes:
      case NumFileDataFields:
        throw runtime_error( string("SpecTest::set_test: ") + to_string(field) + " not a numeric field" );
        break;
    }//switch( m_searchField )
    
    m_numeric = value;
    m_searchField = field;
    m_compareType = type;
  }
  
  
  void SpecTest::set_time_test( const FileDataField field, boost::posix_time::ptime comptime, const NumericFieldMatchType type )
  {
    switch( field )
    {
      case FileDataField::StartTimeIoI:
      case FileDataField::MeasurementsStartTimes:
        break;
        
      case FileDataField::TotalLiveTime: case FileDataField::TotalRealTime:
      case FileDataField::IndividualSpectrumLiveTime: case FileDataField::IndividualSpectrumRealTime:
      case FileDataField::NumberOfSamples: case FileDataField::NumberOfRecords:
      case FileDataField::NumberOfGammaChannels: case FileDataField::MaximumGammaEnergy:
      case FileDataField::Latitude: case FileDataField::Longitude:
      case FileDataField::NeutronCountRate: case FileDataField::GammaCountRate:
      case FileDataField::DetectionSystemType: case FileDataField::SearchMode:
      case FileDataField::ContainedNuetronDetector: case FileDataField::ContainedDeviationPairs:
      case FileDataField::HasGps: case FileDataField::EnergyCalibrationType:
      case FileDataField::ParentPath: case FileDataField::Filename:
      case FileDataField::DetectorName: case FileDataField::SerialNumber:
      case FileDataField::Manufacturer: case FileDataField::Model:
      case FileDataField::Uuid: case FileDataField::Remark:
      case FileDataField::LocationName: case FileDataField::AnalysisResultText:
      case FileDataField::AnalysisResultNuclide: case FileDataField::HasRIIDAnalysis:
      case FileDataField::NumFileDataFields:
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
        return SpecUtils::iequals_ascii( teststr, ss );
        
      case TextNotEqual:
        return !SpecUtils::iequals_ascii( teststr, ss );
        
      case TextIsContained:
        return SpecUtils::icontains( teststr, ss );
        
      case TextDoesNotContain:
        return !SpecUtils::icontains( teststr, ss );
        
      case TextStartsWith:
        return SpecUtils::istarts_with( teststr, ss );
        
      case TextDoesNotStartWith:
        return !SpecUtils::istarts_with( teststr, ss );
        
      case TextEndsWith:
        return SpecUtils::iends_with( teststr, ss );
      case TextDoesNotEndWith:
        return !SpecUtils::iends_with( teststr, ss );
        
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
    auto testdate = to_ptime( SpecUtils::time_from_string(teststr.c_str()) );
    if( testdate.is_special() )
      return false;
    
    if( dt.is_special() )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Was passed in a posix_time::ptime that is invalid!");
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
      case FileDataField::ParentPath:
        return test_string( SpecUtils::parent_path(meas.filename), m_stringSearchType, m_searchString );
        
      case FileDataField::Filename:
        return test_string( meas.filename, m_stringSearchType, m_searchString );
        
      case FileDataField::DetectorName:
      {
        bool found = false;
        for( auto it = begin(meas.detector_names); !found && (it != end(meas.detector_names)); ++it )
          found |= test_string( *it, m_stringSearchType, m_searchString );
        return found;
      }
        
      case FileDataField::SerialNumber:
        return test_string( meas.serial_number, m_stringSearchType, m_searchString );
        
      case FileDataField::Manufacturer:
        return test_string( meas.manufacturer, m_stringSearchType, m_searchString );
        
      case FileDataField::Model:
        return test_string( meas.model, m_stringSearchType, m_searchString );
        
      case FileDataField::Uuid:
        return test_string( meas.uuid, m_stringSearchType, m_searchString );
        
      case FileDataField::Remark:
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
        
      case FileDataField::LocationName:
        return test_string( meas.location_name, m_stringSearchType, m_searchString );
        
      case FileDataField::HasRIIDAnalysis:
        return (m_discreteOption==1 ? meas.has_riid_analysis : !meas.has_riid_analysis);
        
      case FileDataField::AnalysisResultText:
      {
        if( !meas.has_riid_analysis )
          return false;
        
        const vector<string> &remarks = meas.riid_ana.remarks_;
        const vector<SpecUtils::DetectorAnalysisResult> &res = meas.riid_ana.results_;
        
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
          const SpecUtils::DetectorAnalysisResult &r = res[i];
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
        
        
      case FileDataField::AnalysisResultNuclide:
      {
        if( !meas.has_riid_analysis )
          return false;
        
        bool found_nuclide = false;
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        
        // Sometimes user may want to search for a large number of isotopes (e.g., any spectrum
        //  file containing a medical isotope), so we'll let the user "OR" isotopes by separating
        //  them with a comma or semi-colon.
        vector<string> user_wanted_nuclides;
        SpecUtils::split( user_wanted_nuclides, m_searchString, ",;" );
        
        for( const auto &search_nuc : user_wanted_nuclides )
        {
          const SandiaDecay::Nuclide * const nuc = db->nuclide( search_nuc );
          
          //const vector<string> &remarks = ana->remarks_;
          const vector<SpecUtils::DetectorAnalysisResult> &res = meas.riid_ana.results_;
          
          if( nuc )
          {
            for( size_t i = 0; !found_nuclide && i < res.size(); ++i )
            {
              //Maybe there was something like "Am 241" that splitting them up
              //  would mess up.
              const SandiaDecay::Nuclide *testnuc = db->nuclide( res[i].nuclide_ );
              if( nuc == testnuc )
                found_nuclide = true;
              
              //Just in case nuclide has mutliple fields for some reason, we'll
              //  split them up and try; this isnt foolproof,
              //  ex "found Am-241, Th 232" would fail for looking for Th232, but
              //  its better than nothing for now
              vector<string> fields;
              if( !found_nuclide )
                SpecUtils::split( fields, res[i].nuclide_, " \t\n,;" );
              
              for( size_t j = 0; !found_nuclide && j < fields.size(); ++j )
              {
                testnuc = db->nuclide( fields[j] );
                if( !testnuc && ((j+1)<fields.size()) )
                  testnuc = db->nuclide( fields[j] + fields[j+1] ); //Try to make uf for things like "Th 232"
                
                found_nuclide = (nuc == testnuc);
              }
            }//for( size_t i = 0; i < res.size(); ++i )
          }else
          {
            //Go through here and check if the specified text is contained within nuclide_,
            //  nuclide_type_, or remark_
            //
            // Note we arent searching by `m_stringSearchType`, but instead always by
            //  TextIsContained, since the user may have entered multiple comma separated things,
            //  and we want to "OR" those together.
            for( size_t i = 0; !found_nuclide && (i < res.size()); ++i )
            {
              const SpecUtils::DetectorAnalysisResult &r = res[i];
              if( test_string( r.nuclide_, TextIsContained, m_searchString ) )
                found_nuclide = true;
              if( test_string( r.nuclide_type_, TextIsContained, m_searchString ) )
                found_nuclide = true;
              if( test_string( r.remark_, TextIsContained, m_searchString ) )
                found_nuclide = true;
            }//for( size_t i = 0; i < res.size(); ++i )
          }//if( user is searching for valid nuclide ) / else (user searching for like "HEU" )
          
          if( found_nuclide )
            break;
        }//for( const auto search_nuc : user_wanted_nuclides )
        
        assert( (m_stringSearchType == TextDoesNotContain) || (m_stringSearchType == TextIsContained) );
        
        return (m_stringSearchType == TextDoesNotContain) ? !found_nuclide : found_nuclide;
      }//case AnalysisResultNuclide:
        
      case FileDataField::DetectionSystemType:
      {
        const SpecUtils::DetectorType type = SpecUtils::DetectorType(m_discreteOption);
        return (meas.detector_type == type);
        
        break;
      }//case DetectionSystemType:
        
      case FileDataField::SearchMode:
        return ((m_discreteOption==0 && !meas.passthrough) || (m_discreteOption && meas.passthrough));
        
      case FileDataField::TotalLiveTime:
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
        
      case FileDataField::TotalRealTime:
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
        
      case FileDataField::ContainedNuetronDetector:
        return ((m_discreteOption==1) ? meas.contained_neutron : !meas.contained_neutron);
        
      case FileDataField::ContainedDeviationPairs:
        return ((m_discreteOption==1) ? meas.contained_dev_pairs : !meas.contained_dev_pairs);
        
      case FileDataField::HasGps:
        return (m_discreteOption==1 ? meas.contained_gps : !meas.contained_gps);
        
      case FileDataField::EnergyCalibrationType:
      {
        const auto type = SpecUtils::EnergyCalType(m_discreteOption);
        for( const auto cal : meas.energy_cal_types )
          if( cal == type )
            return true;
        return false;
      }//case EnergyCalibrationType:
        
      case FileDataField::IndividualSpectrumLiveTime:
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
        
      case FileDataField::IndividualSpectrumRealTime:
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
        
      case FileDataField::NumberOfSamples:
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
        
      case FileDataField::NumberOfRecords:
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
        
      case FileDataField::NumberOfGammaChannels:
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
        
      case FileDataField::MaximumGammaEnergy:
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
        
      case FileDataField::Latitude:
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
        
      case FileDataField::Longitude:
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
        
      case FileDataField::NeutronCountRate:
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
        
      case FileDataField::GammaCountRate:
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
        
      case FileDataField::StartTimeIoI:
      {
        const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
        const boost::posix_time::time_duration::sec_type test_time = (m_time - epoch).total_seconds();
        
        if( meas.start_time_ioi == 0 )
          return false;
        
        switch( m_compareType )
        {
          case ValueIsExact:
          {
            if( abs(meas.start_time_ioi - test_time) < 60 ) //note: we only have minute resoltion on the gui selector
              return true;
            break;
          }
            
          case ValueIsNotEqual:
          {
            if( abs(meas.start_time_ioi - test_time) > 60 ) //We only have minute resoltion on the gui selector
              return true;
            break;
          }
          case ValueIsLessThan:    if(meas.start_time_ioi < test_time) return true; break;
          case ValueIsGreaterThan: if(meas.start_time_ioi > test_time) return true; break;
        }//switch( m_compareType )
        
        return false;
      }//case FileDataField::StartTimeIoI:
        
      case FileDataField::MeasurementsStartTimes:
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
        
      case FileDataField::NumFileDataFields:
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
      case FileDataField::ParentPath:                 return "ParentPath";
      case FileDataField::Filename:                   return "Filename";
      case FileDataField::DetectorName:               return "Detector name";
      case FileDataField::SerialNumber:               return "Serial number";
      case FileDataField::Manufacturer:               return "Manufacturer";
      case FileDataField::Model:                      return "Model";
      case FileDataField::Uuid:                       return "UUID";
      case FileDataField::Remark:                     return "Remark";
      case FileDataField::LocationName:               return "Location name";
      case FileDataField::AnalysisResultText:         return "RIID Ana result";
      case FileDataField::AnalysisResultNuclide:      return "RIID IDed nuclide";
      case FileDataField::HasRIIDAnalysis:            return "Has RIID Analysis";
      case FileDataField::DetectionSystemType:        return "Detector system";
      case FileDataField::ContainedNuetronDetector:   return "Has Neutron";
      case FileDataField::ContainedDeviationPairs:    return "Has Dev. Pairs";
      case FileDataField::HasGps:                     return "Has GPS Info";
      case FileDataField::EnergyCalibrationType:      return "Energy Cal Type";
      case FileDataField::SearchMode:                 return "Acquisition mode";
      case FileDataField::TotalLiveTime:              return "Sum live time";
      case FileDataField::TotalRealTime:              return "Sum real time";
      case FileDataField::IndividualSpectrumLiveTime: return "Spectrum live time";
      case FileDataField::IndividualSpectrumRealTime: return "Spectrum real time";
      case FileDataField::NumberOfSamples:            return "Num. Time Samples";
      case FileDataField::NumberOfRecords:            return "Num. Records";
      case FileDataField::NumberOfGammaChannels:      return "Num Gamma Channels";
      case FileDataField::MaximumGammaEnergy:         return "Max Gamma Energy";
      case FileDataField::Latitude:                   return "Latitude";
      case FileDataField::Longitude:                  return "Longitude";
      case FileDataField::NeutronCountRate:           return "Neutron CPS";
      case FileDataField::GammaCountRate:             return "Gamma CPS";
      case FileDataField::StartTimeIoI:               return "Start Time";
      case FileDataField::MeasurementsStartTimes:     return "Meas. Start Times";
        
      case FileDataField::NumFileDataFields:
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
      case FileDataField::ParentPath: case FileDataField::Filename:
      case FileDataField::DetectorName: case FileDataField::SerialNumber:
      case FileDataField::Manufacturer: case FileDataField::Model:
      case FileDataField::Uuid: case FileDataField::Remark:
      case FileDataField::LocationName: case FileDataField::AnalysisResultText:
      case FileDataField::AnalysisResultNuclide:
        summary += string(" ") + to_string( m_stringSearchType )
        + string(" \"") + m_searchString + "\"";
        break;
        
      case FileDataField::HasRIIDAnalysis:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have RIID Analysis Results";
        break;
        
      case FileDataField::DetectionSystemType:
        summary += " is a " + SpecUtils::detectorTypeToString( SpecUtils::DetectorType(m_discreteOption) );
        break;
        
      case FileDataField::SearchMode:
        summary += " is " + string(m_discreteOption ? "search/passthrough" : "dwell");
        break;
        
      case FileDataField::ContainedNuetronDetector:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have a neutron detectoor";
        break;
        
      case FileDataField::ContainedDeviationPairs:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have non-linear deviation pairs";
        break;
        
      case FileDataField::HasGps:
        summary += "does " + string(m_discreteOption ? "" : "not ") + " have GPS information";
        break;
        
      case FileDataField::EnergyCalibrationType:
        switch( SpecUtils::EnergyCalType(m_discreteOption) )
        {
          case SpecUtils::EnergyCalType::Polynomial: summary += "polynomial"; break;
          case SpecUtils::EnergyCalType::FullRangeFraction: summary += "full range fraction"; break;
          case SpecUtils::EnergyCalType::LowerChannelEdge: summary += "lower channel energy"; break;
          case SpecUtils::EnergyCalType::InvalidEquationType:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
            summary += "unknown";
            break;
        }
        summary += " energy calibration";
        break;
        
      case FileDataField::TotalLiveTime: case FileDataField::TotalRealTime:
      case FileDataField::IndividualSpectrumLiveTime: case FileDataField::IndividualSpectrumRealTime:
      case FileDataField::NumberOfSamples: case FileDataField::NumberOfRecords:
      case FileDataField::NumberOfGammaChannels: case FileDataField::MaximumGammaEnergy:
      case Latitude: case FileDataField::Longitude:
      case FileDataField::NeutronCountRate: case FileDataField::GammaCountRate:
        summary += " is " + string(to_string(m_compareType)) + " " + PhysicalUnits::printToBestTimeUnits(m_numeric);
        break;
        
      case FileDataField::StartTimeIoI:
      case FileDataField::MeasurementsStartTimes:
        summary += " is " + string(to_string(m_compareType)) + " " + SpecUtils::to_iso_string( std::chrono::time_point_cast<std::chrono::microseconds>(to_time_point(m_time)));
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
      case FileDataField::ParentPath: case FileDataField::Filename:
      case FileDataField::DetectorName: case FileDataField::SerialNumber:
      case FileDataField::Manufacturer: case FileDataField::Model:
      case FileDataField::Uuid: case FileDataField::Remark:
      case FileDataField::LocationName: case AnalysisResultText:
      case FileDataField::AnalysisResultNuclide:
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
        
      case FileDataField::DetectionSystemType:
      {
        //We could check m_discreteOption is
        bool gotit = false;
        switch( SpecUtils::DetectorType(m_discreteOption) )
        {
          case SpecUtils::DetectorType::Exploranium:
          case SpecUtils::DetectorType::IdentiFinder:
          case SpecUtils::DetectorType::IdentiFinderNG:
          case SpecUtils::DetectorType::IdentiFinderLaBr3:
          case SpecUtils::DetectorType::DetectiveUnknown:
          case SpecUtils::DetectorType::DetectiveEx:
          case SpecUtils::DetectorType::DetectiveEx100:
          case SpecUtils::DetectorType::DetectiveEx200:
          case SpecUtils::DetectorType::DetectiveX:
          case SpecUtils::DetectorType::SAIC8:
          case SpecUtils::DetectorType::Falcon5000:
          case SpecUtils::DetectorType::Unknown:
          case SpecUtils::DetectorType::MicroDetective:
          case SpecUtils::DetectorType::MicroRaider:
          case SpecUtils::DetectorType::RadHunterNaI:
          case SpecUtils::DetectorType::RadHunterLaBr3:
          case SpecUtils::DetectorType::Rsi701:
          case SpecUtils::DetectorType::Rsi705:
          case SpecUtils::DetectorType::AvidRsi:
          case SpecUtils::DetectorType::OrtecRadEagleNai:
          case SpecUtils::DetectorType::OrtecRadEagleCeBr2Inch:
          case SpecUtils::DetectorType::OrtecRadEagleCeBr3Inch:
          case SpecUtils::DetectorType::OrtecRadEagleLaBr:
          case SpecUtils::DetectorType::Sam940LaBr3:
          case SpecUtils::DetectorType::Sam940:
          case SpecUtils::DetectorType::Sam945:
          case SpecUtils::DetectorType::Srpm210:
          case SpecUtils::DetectorType::IdentiFinderTungsten:
          case SpecUtils::DetectorType::IdentiFinderR425NaI:
          case SpecUtils::DetectorType::IdentiFinderR425LaBr:
          case SpecUtils::DetectorType::IdentiFinderR500NaI:
          case SpecUtils::DetectorType::IdentiFinderR500LaBr:
          case SpecUtils::DetectorType::IdentiFinderUnknown:
          case SpecUtils::DetectorType::Interceptor:
          case SpecUtils::DetectorType::RIIDEyeNaI:
          case SpecUtils::DetectorType::RIIDEyeLaBr:
          case SpecUtils::DetectorType::RadSeekerNaI:
          case SpecUtils::DetectorType::RadSeekerLaBr:
          case SpecUtils::DetectorType::VerifinderNaI:
          case SpecUtils::DetectorType::VerifinderLaBr:
          case SpecUtils::DetectorType::KromekD3S:
          case SpecUtils::DetectorType::RadiaCode:
          case SpecUtils::DetectorType::Fulcrum:
          case SpecUtils::DetectorType::Fulcrum40h:
          case SpecUtils::DetectorType::Sam950:
            
            gotit = true;
            break;
        }
        
        if( !gotit )
          throw runtime_error( "Invalid DetectorType (programming error)" );
        break;
      }//case DetectionSystemType:
        
      case FileDataField::SearchMode:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid aquitision mode (programming error)" );
        break;
      }
      
      case FileDataField::HasRIIDAnalysis:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid RIID Analysis present option (programming error)" );
        break;
      }
        
      case FileDataField::ContainedNuetronDetector:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid neutron detector presence (programming error)" );
        break;
      }
        
      case FileDataField::ContainedDeviationPairs:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid deviation pair presence (programming error)" );
        break;
      }
        
      case FileDataField::HasGps:
      {
        if( m_discreteOption!=0 && m_discreteOption!=1 )
          throw runtime_error( "Invalid GPS info presence (programming error)" );
        break;
      }
        
      case FileDataField::EnergyCalibrationType:
      {
        switch( static_cast<SpecUtils::EnergyCalType>(m_discreteOption) )
        {
          case SpecUtils::EnergyCalType::Polynomial:
          case SpecUtils::EnergyCalType::FullRangeFraction:
          case SpecUtils::EnergyCalType::LowerChannelEdge:
          case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          case SpecUtils::EnergyCalType::InvalidEquationType:
            break;
          default:
            throw runtime_error( "Invalid energy calibration type (programming error)" );
        }
      }
        
      case FileDataField::TotalLiveTime: case FileDataField::TotalRealTime:
      case FileDataField::IndividualSpectrumLiveTime: case FileDataField::IndividualSpectrumRealTime:
      case FileDataField::NumberOfSamples: case FileDataField::NumberOfRecords:
      case FileDataField::NumberOfGammaChannels: case FileDataField::MaximumGammaEnergy:
      case FileDataField::Latitude: case FileDataField::Longitude:
      case FileDataField::NeutronCountRate: case FileDataField::GammaCountRate:
      {
        if( IsInf(m_numeric) || IsNan(m_numeric) )
          throw runtime_error( "Invalid time duration" );
        break;
      }
        
      case FileDataField::StartTimeIoI:
      case FileDataField::MeasurementsStartTimes:
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
        return "Event XML with comparing '" + m_test_label + "' " + to_string(m_dateTestType) + " '" + SpecUtils::to_common_string(std::chrono::time_point_cast<std::chrono::microseconds>(to_time_point(m_test_time)),true) + "'";
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
  
  std::ostream &SpecLogicTest::print_equation( std::vector<boost::any> fields, std::ostream &strm )
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

