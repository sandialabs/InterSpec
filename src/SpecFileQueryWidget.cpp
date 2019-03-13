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

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLabel>
#include <Wt/WTheme>
#include <Wt/WAnchor>
#include <Wt/WServer>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WIOService>
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/WFileUpload>
#if( WT_VERSION >= 0x3030800 )
#include <Wt/WTimeEdit>
#endif
#include <Wt/WDatePicker>
#include <Wt/WTimePicker>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/Http/Request>
#include <Wt/WApplication>
#include <Wt/WStandardItem>
#include <Wt/Http/Response>
#include <Wt/Json/Serializer>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecFileQuery.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/SpecFileQueryWidget.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/SpecFileQueryDbCache.h"

#include "js/SpecFileQueryWidget.js"


#if( !BUILD_AS_ELECTRON_APP && defined(WIN32) && BUILD_AS_LOCAL_SERVER )
#include <shellapi.h>  //for ShellExecute
#endif

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


using namespace Wt;
using namespace std;


namespace
{
  using namespace SpecFileQuery;
  
  bool file_smaller_than( const std::string &path, void *maxsizeptr )
  {
    const size_t *maxsize = (const size_t *)maxsizeptr;
    if( UtilityFunctions::file_size(path) > (*maxsize) )
      return false;
    
    return true;
  }
  
  bool maybe_spec_file( const std::string &path, void *maxsizeptr )
  {
    if( likely_not_spec_file( path ) )
      return false;
    
    if( maxsizeptr )
      return file_smaller_than( path, maxsizeptr );
    
    return true;
  }
  
  
  void add_logic( SpecFileQuery::SpecLogicTest &test, const Json::Value &cond, const Json::Array &rules,
                  const std::vector<EventXmlFilterInfo> &eventXmlTests )
  {
    using namespace SpecFileQuery;
    
    if( cond.isNull() || rules==Json::Array::Empty )
      throw runtime_error( "Empty condition or logic array" );
    
    const std::string &condstr = cond;
    if( cond != "AND" && cond != "OR" )
      throw runtime_error( "Found condition that is not 'AND' or 'OR;" );
    
    const LogicType logic = ((condstr=="AND") ? LogicType::LogicalAnd : LogicType::LogicalOr);
    
    if( rules.size() > 1 )
      test.addLogic( LogicType::LogicalOpenParan );
    
    for( size_t i = 0; i < rules.size(); ++i )
    {
      if( i )
        test.addLogic( logic );
      
      const Json::Object &t = rules[i];
      if( t.contains("condition") && t.contains("rules") )
      {
        const Json::Value &subCond = t.get( "condition" );
        const Json::Array &subRules = t.get( "rules" );
        
        add_logic( test, subCond, subRules, eventXmlTests );
      }else if( t.contains("field") && t.contains("operator") && t.contains("value") )
      {
        const std::string &field_str = t.get("field");
        const std::string &operator_str = t.get("operator");
        
        FileDataField field;
        
        if( !from_string(field_str, field) )
        {
          const EventXmlFilterInfo *info = nullptr;
          
          for( size_t i = 0; !info && i < eventXmlTests.size(); ++i )
          {
            if( eventXmlTests[i].m_label == field_str )
              info = &(eventXmlTests[i]);
          }
          
          if( !info )
            throw runtime_error( "Could not convert field '" + Wt::Utils::htmlEncode(field_str) + "' to a testable field." );
          
          //If we are here it is an Event XML test
          const std::string &value_str = t.get("value");
          
          EventXmlTest testitem;
          
          switch( info->m_type )
          {
            case EventXmlFilterInfo::InputType::Text:
            case EventXmlFilterInfo::InputType::Select:
            {
              TextFieldSearchType type;
              if( !from_string(operator_str, type) )
                throw runtime_error( "Couldnt not map '" + Wt::Utils::htmlEncode(operator_str) + "' to a TextFieldSearchType." );
              
              testitem.set_string_test_info( info->m_label, value_str, type );
              break;
            }
              
            case EventXmlFilterInfo::InputType::Date:
            {
              NumericFieldMatchType type;
              if( !from_string(operator_str, type) )
                throw runtime_error( "Couldnt not map '" + Wt::Utils::htmlEncode(operator_str) + "' to a NumericFieldMatchType." );
              
              const boost::posix_time::ptime test_time = UtilityFunctions::time_from_string( value_str.c_str() );
              if( test_time.is_special() )
                throw runtime_error( "Could not convert string '" + Wt::Utils::htmlEncode(value_str) + "' to a time" );
              
              testitem.set_date_test_info( info->m_label, test_time, type );
              break;
            }
          }//switch( testinfo.m_type )
          
          test.addEventXmlTest( testitem );
          
          continue;
        }//if( !from_string(field_str, field) )
        
        SpecTest testitem;
        
        switch( field )
        {
          case FileDataField::SearchMode:
          case FileDataField::DetectionSystemType:
          case FileDataField::ContainedNuetronDetector:
          case FileDataField::ContainedDeviationPairs:
          case FileDataField::HasGps:
          case HasRIIDAnalysis:
          case FileDataField::EnergyCalibrationType:
          {
            int val = t.get("value");
          
            if( field == FileDataField::DetectionSystemType )
            {
              //The below switch must match whats in JS - a little brittle, bu oh well for now
              switch( val )
              {
                case 0: val = static_cast<int>(DetectorType::kGR135Detector); break;
                case 1: val = static_cast<int>(DetectorType::kIdentiFinderDetector); break;
                case 2: val = static_cast<int>(DetectorType::kIdentiFinderNGDetector); break;
                case 3: val = static_cast<int>(DetectorType::kIdentiFinderLaBr3Detector); break;
                case 4: val = static_cast<int>(DetectorType::kDetectiveDetector); break;
                case 5: val = static_cast<int>(DetectorType::kDetectiveExDetector); break;
                case 6: val = static_cast<int>(DetectorType::kDetectiveEx100Detector); break;
                case 7: val = static_cast<int>(DetectorType::kDetectiveEx200Detector); break;
                case 8: val = static_cast<int>(DetectorType::kSAIC8Detector); break;
                case 9: val = static_cast<int>(DetectorType::kFalcon5000); break;
                case 10: val = static_cast<int>(DetectorType::kMicroDetectiveDetector); break;
                case 11: val = static_cast<int>(DetectorType::kMicroRaiderDetector); break;
                case 12: val = static_cast<int>(DetectorType::kSam940); break;
                case 13: val = static_cast<int>(DetectorType::kSam940LaBr3); break;
                case 14: val = static_cast<int>(DetectorType::kSam945); break;
                case 15: val = static_cast<int>(DetectorType::kRsi701); break;
                case 16: val = static_cast<int>(DetectorType::kRsi705); break;
                case 17: val = static_cast<int>(DetectorType::kAvidRsi); break;
                case 18: val = static_cast<int>(DetectorType::kRadHunterNaI); break;
                case 19: val = static_cast<int>(DetectorType::kRadHunterLaBr3); break;
                case 20: val = static_cast<int>(DetectorType::kOrtecRadEagleNai); break;
                case 21: val = static_cast<int>(DetectorType::kOrtecRadEagleCeBr2Inch); break;
                case 22: val = static_cast<int>(DetectorType::kOrtecRadEagleCeBr3Inch); break;
                case 23: val = static_cast<int>(DetectorType::kOrtecRadEagleLaBr); break;
                case 24: val = static_cast<int>(DetectorType::kSrpm210); break;
                case 25: val = static_cast<int>(DetectorType::kUnknownDetector); break;
                default:
                  throw runtime_error( "Unknown DetectionSystemType value type" );
              }//switch( val )
            }else if( field == FileDataField::EnergyCalibrationType )
            {
              //The below switch must match whats in JS - a little brittle, bu oh well for now
              switch( val )
              {
                case 0: val = static_cast<int>(Measurement::EquationType::Polynomial); break;
                case 1: val = static_cast<int>(Measurement::EquationType::FullRangeFraction); break;
                case 2: val = static_cast<int>(Measurement::EquationType::LowerChannelEdge); break;
                case 3: val = static_cast<int>(Measurement::EquationType::UnspecifiedUsingDefaultPolynomial); break;
                default:
                  throw runtime_error( "Unknown EnergyCalibrationType value type" );
              }//switch( val )
            }//if( DetectionSystemType ) / elif (EnergyCalibrationType)
          
            testitem.set_discreet_test( field, val );
            break;
          }//case FileDataField::SearchMode:
            
          case FileDataField::TotalLiveTime:
          case FileDataField::TotalRealTime:
          case FileDataField::IndividualSpectrumLiveTime:
          case FileDataField::IndividualSpectrumRealTime:
          {
            const std::string &value_str = t.get("value");
            
            NumericFieldMatchType matchtype;
            
            if( !from_string(operator_str, matchtype) ) //shouldnt ever happen!
              throw runtime_error( "Could not map '" + operator_str + "' to a NumericFieldMatchType" );
            
            //if( !(stringstream(value_str) >> value) ) //Should have been valideted client-side, but lets be safe
            //  throw runtime_error( "Could not convert '" + value_str + "' to a floating point" );
            double value = std::numeric_limits<double>::infinity();
            try{ value = PhysicalUnits::stringToTimeDuration( value_str ); }catch(...){}
            
            testitem.set_numeric_test( field, value, matchtype );
            break;
          }//case( a float field )
            
          case FileDataField::NumberOfSamples:
          case FileDataField::NumberOfRecords:
          case FileDataField::NumberOfGammaChannels:
          case FileDataField::MaximumGammaEnergy:
          case FileDataField::Latitude:
          case FileDataField::Longitude:
          case FileDataField::NeutronCountRate:
          case FileDataField::GammaCountRate:
          {
            NumericFieldMatchType matchtype;
            if( !from_string(operator_str, matchtype) ) //shouldnt ever happen!
              throw runtime_error( "Could not map '" + operator_str + "' to a NumericFieldMatchType" );
            
            const double value = t.get("value");
            testitem.set_numeric_test( field, value, matchtype );
            break;
          }
            
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
          {
            const std::string &value_str = t.get("value");
            
            TextFieldSearchType type;
            if( !from_string(operator_str, type) )
              throw runtime_error( "Couldnt not map '" + operator_str + "' to a TextFieldSearchType." );
            
            testitem.set_test( field, value_str, type );
            break;
          }//case( a string field )
            
          case FileDataField::StartTime:
          {
            const std::string &value_str = t.get("value");
            
            NumericFieldMatchType matchtype;
            if( !from_string(operator_str, matchtype) ) //shouldnt ever happen!
              throw runtime_error( "Could not map '" + operator_str + "' to a NumericFieldMatchType for time comparison" );
            
            boost::posix_time::ptime comptime = UtilityFunctions::time_from_string( value_str.c_str() );
            if( comptime.is_special() ) //Can definetly happen!
              throw runtime_error( "The string '" + value_str + "' is not a valid date/time string." );
            
            testitem.set_time_test( field, comptime, matchtype );
            break;
          }//case StartTime
            
          case FileDataField::NumFileDataFields:
            assert( 0 );
            break;
        }//switch( field )
        
        test.addCondition( testitem );
      }else
      {
        throw runtime_error( "Unexpected query JSON of: " + Json::serialize(t) + "" );
      }
    }//for( size_t i = 0; i < rules.size(); ++i )
    
    if( rules.size() > 1 )
      test.addLogic( LogicType::LogicalCloseParan );
  }//void add_logic(...)
  
  
  void output_csv_field( std::ostream &out, std::string s)
  {
    //adapted from http://darrendev.blogspot.com/2009/11/escaping-csv-in-c.html 20160623
    //There are two escaping rules for each field in a comma-separated value row:
    //  1. Change each double quote to two double quotes.
    //  2. Surround with double quotes if the field contains a comma or double quote.
    
    //ToDo: see if this agrees with boost::io::quoted(...)
    
    if( s.find('"') != std::string::npos ) //Escape double-quotes
    {
      std::string::size_type pos = 0;
      while( 1 )
      {
        pos = s.find( '"', pos );
        if( pos == std::string::npos )
          break;
        s.replace(pos,1,"\"\"");
        pos += 2; //Need to skip over those two quotes, to avoid an infinite loop!
      }
      
      out << '"' << s << '"';
    }else if( s.find(',') != std::string::npos
              || s.find('\n') != std::string::npos
              || s.find('\r') != std::string::npos )  //Need to surround with "..."
    {
      out << '"' << s << '"';
    }else
      out << s;
  }//output_csv_field
  
  
  
  class ResultCsvResource : public Wt::WResource
  {
  private:
    WAbstractItemModel *m_model;
    
  public:
    ResultCsvResource( WAbstractItemModel *model, WPushButton *parent )
     : WResource( parent ),
       m_model( model )
    {
    }
    
    virtual ~ResultCsvResource()
    {
      beingDeleted();
    }
    
  
    virtual void handleRequest( const Wt::Http::Request &request,
                                Wt::Http::Response &response )
    {
      suggestFileName( "FileSearchResults.csv", WResource::Attachment ); //WResource::NoDisposition
      response.setMimeType( "text/csv" );
      
      if( !m_model )
        return;
      
      boost::any description = m_model->headerData(0, Horizontal, Wt::UserRole);
      response.out() << asString(description).toUTF8() << "\r\n";
      
      const int ncol = m_model->columnCount();
      for( int col = 0; col < ncol; ++col )
      {
        const boost::any title = m_model->headerData( col, Horizontal, DisplayRole );
        const WString titlestr = Wt::asString( title );
        response.out() << (col ? "," : "");
        output_csv_field(	response.out(), titlestr.toUTF8() );
      }
      response.out() << "\r\n";
      
      const int nrow = m_model->rowCount();
      for( int row = 0; row < nrow; ++row )
      {
        for( int col = 0; col < ncol; ++col )
        {
          const boost::any txt = m_model->data( row, col, Wt::DisplayRole /* (col==FileDataField::Filename ? UserRole : DisplayRole)*/ );
          const WString txtstr = Wt::asString( txt );
          const string utf8txt = txtstr.toUTF8();
          
          response.out() << (col ? "," : "");
          output_csv_field(	response.out(), utf8txt );
        }
        
        response.out() << "\r\n";
      }//for( int row = 0; row < nrow; ++row )
    }//handleRequest

  };//class ResultCsvResource : public Wt::WResource

  
  vector<string> get_result_fields( const SpecFileInfoToQuery &meas, const string &base_search_dir,
                                    const std::vector<EventXmlFilterInfo> &xmlfilters )
  {
    const size_t max_cell_size = 1024;
    vector<string> row( NumFileDataFields + xmlfilters.size(), "" );
    
    for( FileDataField f = FileDataField(0); f < NumFileDataFields; f = FileDataField(f+1) )
    {
      switch( f )
      {
        case FileDataField::ParentPath:
          try
          {
            row[f] = UtilityFunctions::fs_relative( base_search_dir, UtilityFunctions::parent_path(meas.filename) );
          }catch(...)
          {
#if( PERFORM_DEVELOPER_CHECKS )
            log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected exception getting relative path between fi" );
#endif
            row[f] = UtilityFunctions::parent_path(meas.filename);
            if( UtilityFunctions::starts_with(row[f], base_search_dir.c_str()) )
              row[f] = row[f].substr(base_search_dir.size());
          }
          break;
          
        case FileDataField::Filename:
          row[f] = meas.filename;
          break;
          
        case FileDataField::DetectorName:
          for( const auto &name : meas.detector_names )
            row[f] += ((row[f].empty() || name.empty()) ? "" : ";") + name;
          break;
          
        case FileDataField::SerialNumber:
          row[f] = meas.serial_number;
          break;
          
        case FileDataField::Manufacturer:
          row[f] = meas.manufacturer;
          break;
          
        case FileDataField::Model:
          row[f] = meas.model;
          break;
          
        case FileDataField::Uuid:
          row[f] = meas.uuid;
          break;
          
        case FileDataField::Remark:
          for( const auto &remark : meas.file_remarks )
          {
            if( row[f].size() < max_cell_size && !remark.empty() )
              row[f] += (row[f].size() ? "\n" : "") + remark;
          }
          
          for( const auto &remark : meas.record_remarks )
          {
            if( row[f].size() < max_cell_size && !remark.empty() )
              row[f] += (row[f].size() ? "\n" : "") + remark;
          }
          break;
          
        case FileDataField::LocationName:
          row[f] = meas.location_name;
          break;
          
        case HasRIIDAnalysis:
          row[f] = (meas.has_riid_analysis ? "Has RIID Ana" : "No RIID Ana");
          break;
          
        case FileDataField::AnalysisResultText:
        {
          if( meas.has_riid_analysis )
          {
            const DetectorAnalysis &anares = meas.riid_ana;
            for( size_t i = 0; i < anares.remarks_.size() && row[f].size() < max_cell_size; ++i )
            {
              if( anares.remarks_[i].size() )
                row[f] += (row[f].size() ? "\n" : "") + anares.remarks_[i];
            }
            
            for( size_t i = 0; i < anares.results_.size() && row[f].size() < max_cell_size; ++i )
            {
              string r;
              if( anares.results_[i].nuclide_type_.size() )
                r += anares.results_[i].nuclide_type_ + ";";
              if( anares.results_[i].id_confidence_.size() )
                r += (r.size() ? ";" : "") + anares.results_[i].id_confidence_ + ";";
              if( anares.results_[i].remark_.size() )
                r += (r.size() ? ";" : "") + anares.results_[i].remark_ + ";";
              if( r.size() )
                row[f] += (row[f].size() ? "\n" : "") + r;
            }
          }//if( have analysis results )
          break;
        }//case FileDataField::AnalysisResultText:
          
          
        case FileDataField::AnalysisResultNuclide:
        {
          if( meas.has_riid_analysis )
          {
            for( size_t i = 0; i < meas.riid_ana.results_.size() && row[f].size() < max_cell_size; ++i )
              if( meas.riid_ana.results_[i].nuclide_.size() )
                row[f] += (row[f].size() ? ";" : "") + meas.riid_ana.results_[i].nuclide_;
          }
          break;
        }//case FileDataField::AnalysisResultNuclide:
          
        case FileDataField::DetectionSystemType:
          row[f] = detectorTypeToString( meas.detector_type );
          break;
          
        case FileDataField::SearchMode:
          row[f] = (meas.passthrough ? "Search/Passthrough" : "Dwell");
          break;
          
        case FileDataField::ContainedNuetronDetector:
          row[f] = (meas.contained_neutron ? "Has Neutron" : "No Neutron");
          break;
          
        case FileDataField::ContainedDeviationPairs:
          row[f] = (meas.contained_dev_pairs ? "Has Dev Pair" : "No Dev Pair");
          break;
          
        case FileDataField::HasGps:
          row[f] = (meas.contained_gps ? "Has GPS" : "No GPS");
          break;
          
        case FileDataField::EnergyCalibrationType:
        {
          for( const auto &caltype : meas.energy_cal_types )
          {
            row[f] += (row[f].size() ? ";" : "");
            switch( caltype )
            {
              case Measurement::EquationType::Polynomial:          row[f] += "polynomial"; break;
              case Measurement::EquationType::FullRangeFraction:   row[f] += "full range fraction"; break;
              case Measurement::EquationType::LowerChannelEdge:    row[f] += "lower channel energy"; break;
              case Measurement::EquationType::InvalidEquationType:
              case Measurement::EquationType::UnspecifiedUsingDefaultPolynomial:
                row[f] += "unknown";
                break;
            }//switch( m->energy_calibration_model() )
          }//for( const auto &m : meas->measurements() )
          break;
        }//case FileDataField::EnergyCalibrationType:
          
        case FileDataField::TotalLiveTime:
          row[f] = PhysicalUnits::printToBestTimeUnits( meas.total_livetime );
          break;
          
        case FileDataField::TotalRealTime:
          row[f] = PhysicalUnits::printToBestTimeUnits( meas.total_realtime );
          break;
          
        case FileDataField::IndividualSpectrumLiveTime:
        {
          for( const auto lt : meas.individual_spectrum_live_time )
          {
            row[f] += (row[f].size()?";":"") + PhysicalUnits::printToBestTimeUnits(lt);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::IndividualSpectrumLiveTime:
          
        case FileDataField::IndividualSpectrumRealTime:
        {
          for( const auto lt : meas.individual_spectrum_real_time )
          {
            row[f] += (row[f].size()?";":"") + PhysicalUnits::printToBestTimeUnits(lt);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::IndividualSpectrumRealTime:
          
        case FileDataField::NumberOfSamples:
          row[f] = std::to_string( meas.number_of_samples );
          break;
          
        case FileDataField::NumberOfRecords:
          row[f] = std::to_string( meas.number_of_records );
          break;
          
        case FileDataField::NumberOfGammaChannels:
        {
          for( const auto nchan : meas.number_of_gamma_channels )
          {
            row[f] += (row[f].size()?";":"") + std::to_string(nchan);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::NumberOfGammaChannels:
          
        case FileDataField::MaximumGammaEnergy:
        {
          for( const auto ene : meas.max_gamma_energy )
          {
            char buff[64];
            snprintf( buff, sizeof(buff), "%.2f", ene);
            row[f] += (row[f].empty()?"":";") + string(buff);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::MaximumGammaEnergy:
          
        case FileDataField::Latitude:
          if( meas.contained_gps )
            row[f] = std::to_string( meas.mean_latitude );
          break;
          
        case FileDataField::Longitude:
          if( meas.contained_gps )
            row[f] = std::to_string( meas.mean_longitude );
          break;
          
        case FileDataField::NeutronCountRate:
        {
          for( const float cps : meas.neutron_count_rate )
          {
            char buff[64];
            snprintf( buff, sizeof(buff), "%.2g", cps);
            row[f] += (row[f].empty()?"":";") + string(buff);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::NeutronCountRate:
          
        case FileDataField::GammaCountRate:
        {
          for( const float cps : meas.gamma_count_rate )
          {
            char buff[64];
            snprintf( buff, sizeof(buff), "%.2g", cps);
            row[f] += (row[f].empty()?"":";") + string(buff);
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }//case FileDataField::NeutronCountRate:
          
        case FileDataField::StartTime:
        {
          for( const auto &st : meas.start_times )
          {
            const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
            row[f] += (row[f].size()?";":"") + UtilityFunctions::to_iso_string( epoch + boost::posix_time::seconds(st) );
            if( row[f].size() > max_cell_size )
              break;
          }
          break;
        }
          
        case FileDataField::NumFileDataFields:
          break;
      }//switch( f )
    }//for( FileDataField f = FileDataField(0); f < NumFileDataFields; f = FileDataField(f+1) )
    
    
    for( size_t i = 0; i < xmlfilters.size(); ++i )
    {
      const string &label = xmlfilters[i].m_label;
      const auto pos = meas.event_xml_filter_values.find(label);
      if( pos == end(meas.event_xml_filter_values) )
        continue;
      
      const int index = static_cast<int>(i) + NumFileDataFields;
      
      for( const string &res : pos->second )
      {
        row[index] += (row[index].size() ? "; " : "") + res;
      }//for( size_t result = 0; result < pos->second.size(); ++result )
    }//for( size_t i = 0; i < xmlfilters.size(); ++i )
    
    //ToDo: do we really need to replace commas and quotes?  I think we properly
    //      quote things when saving to CSV.  Oh well for now.
    for( string &col : row )
    {
      UtilityFunctions::ireplace_all(col,",", " ");
      UtilityFunctions::ireplace_all(col, "\"", "&quot;");
      
      if( col.size() > max_cell_size )
      {
        UtilityFunctions::utf8_limit_str_size( col, max_cell_size-5 ); //count '\n' as twoa characters for windows
        col += "\n...";
      }
    }//for( string &col : row )
    
    return row;
  }//vector<string> get_result_fields( const SpecFileInfoToQuery &meas )
  
  

  struct HaveSeenUuid
  {
    HaveSeenUuid( bool require )
      : m_require( require )
    {
    }
    
    const bool m_require;
    set<string> m_uuids;
    std::mutex m_mutex;
    
    bool have_seen( const string &uuid )
    {
      if( !m_require )
        return false;
      
      std::unique_lock<std::mutex> lock( m_mutex );
      const set<string>::iterator pos = m_uuids.find( uuid );
      
      if( pos == m_uuids.end() )
      {
        m_uuids.insert( uuid );
        return false;
      }
      
      return true;
    }
  };//struct HaveSeenUuid
}//namespace


/** Class made to hold results of the file query.  Had to specialized because
    WStandardItemModel requires you to add new data one row at a time, which 
    would cause a change in the RowStretchTreeView, so it would take over
    a minute just to insert new results for the Data file set that had
    ~20k rows.  This also allows us to be a little more memorry efficient as
    well.
 */
class ResultTableModel : public WAbstractItemModel
{
protected:
  std::shared_ptr< vector< vector<string> > > m_result;
  int m_sortcolumn;
  Wt::SortOrder m_sortorder;
  boost::any m_summary;
  std::vector<std::string> m_eventXmlColNames;
  
  
  struct row_less_than_t
  {
    const int m_column;
    const SortOrder m_order;
    
    row_less_than_t( const int column, const SortOrder order )
    : m_column( column ), m_order( order )
    {
    }
    
    bool operator()( const vector<string> &lhs, const vector<string> &rhs ) const
    {
      if( m_column < 0 )
        return false;
      
      const int nlhs = static_cast<int>(lhs.size());
      const int nrhs = static_cast<int>(rhs.size());
      
      if( m_column >= nlhs )
        return true;
      
      if( m_column >= nrhs )
        return false;
      
      bool lesthan = false;
      
      //Take care of when it is an Event XML column
      if( m_column >= NumFileDataFields )
      {
        lesthan = (lhs[m_column] < rhs[m_column]);
        return ((m_order==Wt::AscendingOrder) ? lesthan : !lesthan);
      }//if( m_column >= NumFileDataFields )
      
      
      switch( FileDataField(m_column) )
      {
        case FileDataField::ParentPath:  //be lazy and just use full file path
        case Filename:
          lesthan = (UtilityFunctions::filename(lhs[m_column]) < UtilityFunctions::filename(rhs[m_column]));
          break;
          
        case DetectorName: case SerialNumber: case Manufacturer:
        case Model: case Uuid: case Remark: case LocationName:
        case AnalysisResultText: case AnalysisResultNuclide:
        case HasRIIDAnalysis:
        case DetectionSystemType: case SearchMode:
        case ContainedNuetronDetector: case ContainedDeviationPairs:
        case HasGps: case EnergyCalibrationType:
          lesthan = (lhs[m_column] < rhs[m_column]);
          break;
        
        case TotalLiveTime:
        case TotalRealTime:
        case IndividualSpectrumLiveTime:
        case IndividualSpectrumRealTime:
        {
          const size_t lhssemi = lhs[m_column].find( ';' );
          const size_t rhssemi = lhs[m_column].find( ';' );
          const string lhsstr = ((lhssemi==string::npos) ? lhs[m_column] : lhs[m_column].substr(0,lhssemi));
          const string rhsstr = ((rhssemi==string::npos) ? rhs[m_column] : rhs[m_column].substr(0,lhssemi));
          
          bool pl = true, pr = true;
          double lhsval, rhsval;
          try
          {
            lhsval = PhysicalUnits::stringToTimeDuration( lhsstr );
          }catch(...)
          {
            pl = false;
          }
          
          try
          {
            rhsval = PhysicalUnits::stringToTimeDuration( rhsstr );
          }catch(...)
          {
            pr = false;
          }
          
          if( pl && pr )
            lesthan = (lhsval < rhsval);
          else
            lesthan = (pl < pr);
          break;
        }
          
        case NumberOfSamples:
        case NumberOfRecords:
        case NumberOfGammaChannels:
        case MaximumGammaEnergy:
        case Latitude:
        case Longitude:
        case NeutronCountRate:
        case GammaCountRate:
        {
          const size_t lhssemi = lhs[m_column].find( ';' );
          const size_t rhssemi = lhs[m_column].find( ';' );
          const string lhsstr = ((lhssemi==string::npos) ? lhs[m_column] : lhs[m_column].substr(0,lhssemi));
          const string rhsstr = ((rhssemi==string::npos) ? rhs[m_column] : rhs[m_column].substr(0,lhssemi));
          
          double lhsval, rhsval;
          const bool pl = !!(stringstream(lhsstr) >> lhsval);
          const bool pr = !!(stringstream(rhsstr) >> rhsval);
          
          if( pl && pr )
            lesthan = (lhsval < rhsval);
          else
            lesthan = (pl < pr);
          break;
        }//
          
        case StartTime:
        {
          const size_t lhssemi = lhs[m_column].find( ';' );
          const size_t rhssemi = lhs[m_column].find( ';' );
          const string lhsstr = ((lhssemi==string::npos) ? lhs[m_column] : lhs[m_column].substr(0,lhssemi));
          const string rhsstr = ((rhssemi==string::npos) ? rhs[m_column] : rhs[m_column].substr(0,lhssemi));
          
          const boost::posix_time::ptime l = UtilityFunctions::time_from_string( lhsstr.c_str() );
          const boost::posix_time::ptime r = UtilityFunctions::time_from_string( rhsstr.c_str() );
          
          lesthan = (l < r);
          
          break;
        }
          
        case NumFileDataFields:
          break;
      }//switch( FileDataField(column) )
      
      return ((m_order==Wt::AscendingOrder) ? lesthan : !lesthan);
    }
  };
  
  
  
  
public:
  ResultTableModel( WObject *parent = 0 )
  : WAbstractItemModel( parent ),
  m_sortcolumn( -1 ),
  m_sortorder( Wt::AscendingOrder )
  {
    
  }
  
  void setEventXmlColumns( const std::vector<std::string> &collnames )
  {
    if( collnames == m_eventXmlColNames )
      return;
    
    const int norig = static_cast<int>( m_eventXmlColNames.size() );
    const int nafter = static_cast<int>( collnames.size() );
    
    if( nafter > norig )
      beginInsertColumns( WModelIndex(), NumFileDataFields+norig, NumFileDataFields+nafter-1 );
    else if( norig > nafter )
      beginRemoveColumns( WModelIndex(), NumFileDataFields+nafter, NumFileDataFields+norig-1 );
    else
      layoutAboutToBeChanged().emit();
      
    m_eventXmlColNames = collnames;
    
    if( nafter > norig )
      endInsertColumns();
    else if( norig > nafter )
      endRemoveColumns();
    else
      layoutChanged().emit();
  }//void setEventXmlColumns( const std::vector<std::string> &collnames )
  
  
  void appendResults( std::shared_ptr< vector< vector<string> > > result )
  {
    if( !result || result->empty() )
      return;
    
    if( !m_result )
      m_result = make_shared< vector<vector<string>> >();
    
    if( result && result->size() )
    {
      const int nprevrow = static_cast<int>( m_result->size() );
      const int nrowadd = static_cast<int>( result->size() );
      
      beginInsertRows( WModelIndex(), nprevrow, nprevrow + nrowadd - 1 );
      m_result->insert( end(*m_result), begin(*result), end(*result) );
      endInsertRows();
      
      if( m_sortcolumn >= 0 )
      {
        layoutAboutToBeChanged().emit();
        std::stable_sort( m_result->begin(), m_result->end(), row_less_than_t(m_sortcolumn,m_sortorder) );
        layoutChanged().emit();
      }
    }
  }//void appendResults( std::shared_ptr< vector< vector<string> > > result )
  
  
  void setResult( std::shared_ptr< vector< vector<string> > > result )
  {
    if( m_result && m_result->size() )
    {
      beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_result->size() -1) );
      m_result.reset();
      endRemoveRows();
    }
    
    if( result && result->size() )
    {
      if( m_sortcolumn >= 0 )
        std::stable_sort( result->begin(), result->end(), row_less_than_t(m_sortcolumn,m_sortorder) );
      
      beginInsertRows( WModelIndex(), 0, static_cast<int>(result->size() - 1) );
      m_result = result;
      endInsertRows();
    }
  }
  
  virtual int columnCount( const WModelIndex &parent = WModelIndex() ) const
  {
    if( parent.isValid() )
      return 0;
    return NumFileDataFields + static_cast<int>( m_eventXmlColNames.size() );
  }
  
  virtual int rowCount( const WModelIndex &parent = WModelIndex() ) const
  {
    if( !m_result )
      return 0;
    return static_cast<int>(m_result->size());
  }
  
  
  virtual WModelIndex parent( const WModelIndex & ) const
  {
    return WModelIndex();
  }
  
  
  virtual boost::any data( const WModelIndex &index, int role = DisplayRole ) const
  {
    if( !index.isValid() || !m_result )
      return boost::any();
    
    const int row = index.row();
    const int col = index.column();
    if( row < 0 || col < 0 || row >= static_cast<int>(m_result->size()) )
      return boost::any();
    
    if( (col >= NumFileDataFields)
       && ((col - NumFileDataFields) >= static_cast<int>(m_eventXmlColNames.size())) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected (to large of) column requested (larger than ever expected)" );
#endif
      return boost::any();
    }
    
    const vector<string> &fields = (*m_result)[row];
    
    if( col == FileDataField::Filename )
    {
      if( role == DisplayRole )
        return WString(UtilityFunctions::filename(fields[col]));
      if( role == UserRole )
        return WString(fields[col]);
      return boost::any();
    }
    
    if( role != DisplayRole )
      return boost::any();
    
    if( col >= static_cast<int>(fields.size()) )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected (to large of) column requested we dont have data for" );
#endif
      return boost::any();
    }
    
    return WString( fields[col] );
  }//data(...)
  
  
  virtual bool setHeaderData( int section, Orientation orientation, const boost::any &value, int role = EditRole )
  {
    if( section != 0 || orientation != Horizontal || role != Wt::UserRole )
      return false;
    m_summary = value;
    return true;
  }
  
  virtual boost::any headerData( int section, Orientation orientation = Horizontal, int role = DisplayRole ) const
  {
    if( section < 0 || orientation != Horizontal )
      return boost::any();
    
    if( section == 0 && role == Wt::UserRole )
      return m_summary;
    
    if( role != DisplayRole )
      return boost::any();
    
    if( section >= NumFileDataFields )
    {
      const int index = section - NumFileDataFields;
      if( index >= static_cast<int>(m_eventXmlColNames.size()) )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected (to large of) column requested for header data" );
#endif
        return boost::any();
      }
      return WString(m_eventXmlColNames[index]);
    }//if( section >= NumFileDataFields )
    
    
    switch( SpecFileQuery::FileDataField(section) )
    {
      case FileDataField::ParentPath:  return WString("Parent Path");
      case Filename:                   return WString("Filename");
      case DetectorName:               return WString("Detectors Name");
      case SerialNumber:               return WString("Serial Number");
      case Manufacturer:               return WString("Manufacturer");
      case Model:                      return WString("Model");
      case Uuid:                       return WString("UUID");
      case Remark:                     return WString("Remarks");
      case LocationName:               return WString("Location Name");
      case AnalysisResultText:         return WString("Analysis Txt");
      case HasRIIDAnalysis:            return WString("RIID Results");
      case AnalysisResultNuclide:      return WString("Analysis Nuclide");
      case DetectionSystemType:        return WString("Detector Type");
      case SearchMode:                 return WString("Aquisition Mode");
      case ContainedNuetronDetector:   return WString("Has Neutron");
      case ContainedDeviationPairs:    return WString("Dev. Pairs");
      case HasGps:                     return WString("GPS Info");
      case EnergyCalibrationType:      return WString("Energy Cal. Type");
      case TotalLiveTime:              return WString("Sum Live Time");
      case TotalRealTime:              return WString("Sum Wall Time");
      case IndividualSpectrumLiveTime: return WString("Spec Live Times");
      case IndividualSpectrumRealTime: return WString("Spec Wall Times");
      case NumberOfSamples:            return WString("Num Samples");
      case NumberOfRecords:            return WString("Num Records");
      case NumberOfGammaChannels:      return WString("Num Gamma Channels");
      case MaximumGammaEnergy:         return WString("Max Gamma Energy");
      case Latitude:                   return WString("Latitude");
      case Longitude:                  return WString("Longitude");
      case NeutronCountRate:           return WString("Neutron CPS");
      case GammaCountRate:             return WString("Gamma CPS");
      case StartTime:                  return WString("Start Time");
        
      case SpecFileQuery::NumFileDataFields: break;
    }
      
    return boost::any();
  }//headerData
  
  virtual WModelIndex index(int row, int column, const WModelIndex &parent = WModelIndex() ) const
  {
    if( parent.isValid() )
      return WModelIndex();
    return WAbstractItemModel::createIndex( row, column, (void *)0 );
  }
  
  
  virtual void sort( int column, SortOrder order = AscendingOrder )
  {
    if( column == m_sortcolumn && order == m_sortorder )
      return;
    
    m_sortcolumn = column;
    m_sortorder = order;
    
    if( !m_result || (m_result->size() < 2) )
      return;
    
    layoutAboutToBeChanged().emit();
    
    std::stable_sort( m_result->begin(), m_result->end(), row_less_than_t(m_sortcolumn,m_sortorder) );
    
    layoutChanged().emit();
  }
};//class ResultTableModel





SpecFileQueryWidget::SpecFileQueryWidget( InterSpec *viewer, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_stopUpdate( std::make_shared< std::atomic<bool> >(false) ),
    m_widgetDeleted( std::make_shared< std::atomic<bool> >(false) ),
    m_app( wApp ),
    m_viewer( viewer ),
    m_conditions( nullptr ),
    m_update( nullptr ),
    m_cancelUpdate( nullptr ),
    m_loadSelectedFile( nullptr ),
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
    m_openSelectedDir( nullptr ),
#endif
    m_resultmodel( nullptr ),
    m_resultview( nullptr ),
#if( BUILD_AS_ELECTRON_APP )
    m_pathSelectedSignal( nullptr ),
#else
    m_baseLocation( nullptr ),
#endif
    m_recursive( nullptr ),
    m_filterByExtension( nullptr ),
    m_maxFileSize( nullptr ),
    m_filterUnique( nullptr ),
    m_cacheParseResults( nullptr ),
    m_persistCacheResults( nullptr ),
    m_optionsBtn( nullptr ),
    m_optionsMenu( nullptr ),
    m_numberFiles( nullptr ),
    m_numberResults( nullptr ),
    m_csv( nullptr ),
    m_queryChanged( this, "fileSearchQueryChanged", false ),
    m_searchRequested( this, "fileSearchRequested", false )
{
  init();
}//SpecFileQueryWidget constructor


SpecFileQueryWidget::~SpecFileQueryWidget()
{
  m_widgetDeleted->store( true );
  
  for( auto &c : m_path_caches )
  {
    if( c.second )  //should always be true
      c.second->stop_caching();
  }
  
#if( BUILD_AS_OSX_APP )
  SpecFileQuery::setIsSelectingDirectory( false );
  SpecFileQuery::setSearchDirectory( "" );
#endif
}//~SpecFileQueryWidget()


void SpecFileQueryWidget::init()
{
  wApp->useStyleSheet( "InterSpec_resources/SpecFileQueryWidget.css" );
  
  wApp->useStyleSheet( "InterSpec_resources/assets/js/QueryBuilder2.5.2/css/query-builder.default.min.css" );
  wApp->useStyleSheet( "InterSpec_resources/assets/js/QueryBuilder2.5.2/css/QueryBuilderFakeBootstrap.css" );
  
  wApp->require( "InterSpec_resources/assets/js/moment-2.22.2/moment.min.js" );
  
  wApp->require( "InterSpec_resources/assets/js/QueryBuilder2.5.2/js/query-builder.standalone.min.js" );
  wApp->require( "InterSpec_resources/assets/js/QueryBuilder2.5.2/i18n/query-builder.en.js" );
  
  LOAD_JAVASCRIPT(wApp, "js/SpecFileQueryWidget.js", "SpecFileQueryWidget", wtjsFileQueryInit);
  
  const string addfilters = prepareEventXmlFilters();
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 9, 9, 9, 0 );
  setLayout( layout );
  
  WContainerWidget *line = new WContainerWidget();
  WGridLayout *linelayout = new WGridLayout( line );
  
  string defpath;
  int maxsize = 32;
  bool dofilter = true, dorecursive = true, docache = true, instantToolTip = true;
  if( m_viewer )
  {
    try
    {
      dofilter = InterSpecUser::preferenceValue<bool>( "SpecFileQueryFilter", m_viewer );
      docache = InterSpecUser::preferenceValue<bool>( "SpecFileQueryCacheParse", m_viewer );
      dorecursive = InterSpecUser::preferenceValue<bool>( "SpecFileQueryRecursive", m_viewer );
      maxsize = InterSpecUser::preferenceValue<int>( "SpecFileQueryMaxSize", m_viewer );
      defpath = InterSpecUser::preferenceValue<string>( "SpecFileQueryPath", m_viewer );
      instantToolTip = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );

      
      if( defpath == "None" )
        defpath = "";
      
      if( !UtilityFunctions::is_directory( defpath ) )
        defpath = "";
      
      if( maxsize < 1 || maxsize > 2048 )
        maxsize = 32;
    }catch(...)
    {
      cerr << "Failed to get def path" << endl;
    }
  }
  
  WLabel *label = new WLabel( "Base Path:" );
  linelayout->addWidget( label, 0, 0 );
  label->setMargin( 3, Wt::Top );
  
#if( BUILD_AS_ELECTRON_APP )
  //Based see https://jaketrent.com/post/select-directory-in-electron/
  m_pathSelectedSignal.reset( new Wt::JSignal<std::string>( this, "BaseDirSelected", false ) );
  
  const string uploadname = id() + "PathPicker";
  const string uploadhtml = "<input id=\"" + uploadname + "\" type=\"file\" webkitdirectory=\"\" />";
  
  WText *uploadtext = new WText( uploadhtml, XHTMLUnsafeText );
  linelayout->addWidget( uploadtext, 0, 1 );
  
  //TODO: put in error handling!
  wApp->doJavaScript( "document.getElementById('" + uploadname + "').onchange = function(event){"
                         "var outputDir = document.getElementById('" + uploadname + "').files[0].path;"
                         "Wt.emit( \"" + id() + "\", { name: 'BaseDirSelected' }, outputDir );"
                     "};"
                     );
  m_pathSelectedSignal->connect( boost::bind( &SpecFileQueryWidget::newElectronPathSelected, this, _1 ) );
#elif( BUILD_AS_OSX_APP )
  SpecFileQuery::setIsSelectingDirectory( true );
  setSearchDirectory( "" );
  m_baseLocation = new WFileUpload();
  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::newMacOsPathSelected );
  linelayout->addWidget( m_baseLocation, 0, 1 );
#else
  m_baseLocation = new WLineEdit();
  
  if( !defpath.empty() )
    m_baseLocation->setText( defpath );
  linelayout->addWidget( m_baseLocation, 0, 1 );
  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_baseLocation->enterPressed().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_baseLocation->setEmptyText( "Path to base directory to search" );
#endif
  
  m_optionsBtn = new WPushButton( "Options" );
  m_optionsBtn->addStyleClass( "SpecFileQueryOptionsBtn" );
  
  m_optionsMenu = new PopupDivMenu( m_optionsBtn, PopupDivMenu::MenuType::TransientMenu );
  linelayout->addWidget( m_optionsBtn, 0, 2 );
  
  auto item = m_optionsMenu->addMenuItem( "recursive" );
  item->setCheckable( true );
  const char *tooltip = "Recursively searches sub-directories or the selected directory.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  
  m_recursive = item->checkBox();
  if( !m_recursive ) //This shouldnt ever happen, but JIC
  {
    m_recursive = new WCheckBox( "recursive" );
    item->addWidget(m_recursive);
  }
  
  m_recursive->setChecked( dorecursive );
  m_recursive->checked().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_recursive->unChecked().connect( this, &SpecFileQueryWidget::basePathChanged );
  
  linelayout->setColumnStretch( 1, 1 );
  
  layout->addWidget( line, 0, 0 );
  
  line = new WContainerWidget();
  line->setStyleClass( "SpecFileQueryConditions" );
  m_conditions = new WContainerWidget( line );
  
  layout->addWidget( line, 1, 0 );
  layout->setRowStretch( 1, 3 );
  
  line = new WContainerWidget();
  linelayout = new WGridLayout( line );
  
//  linelayout->addWidget( new WContainerWidget(), 0, linelayout->columnCount() );
//  linelayout->setColumnStretch( linelayout->columnCount()-1, 1 );

  auto maxFileSizeDiv = new WContainerWidget();
  label = new WLabel( "Max file size:", maxFileSizeDiv );
  label->setMargin( 4, Wt::Bottom );
  label->setMargin( 15, Wt::Left );
  
  m_maxFileSize = new WSpinBox( maxFileSizeDiv );
  m_maxFileSize->setMaximum( 2*1024 );
  m_maxFileSize->setMinimum( 0 );
  m_maxFileSize->setValue( maxsize );
  m_maxFileSize->setWidth( 55 );
  label->setBuddy( m_maxFileSize );
  m_maxFileSize->valueChanged().connect( this, &SpecFileQueryWidget::basePathChanged );
  
  label = new WLabel( "MB", maxFileSizeDiv );
  label->setMargin( 6, Wt::Right );
  label->setMargin( 4, Wt::Bottom );
  label->setBuddy( m_maxFileSize );

  item = m_optionsMenu->addWidget( maxFileSizeDiv );
  tooltip = "Maximum size of file to attempt to parse as a spectrum file.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  
  item = m_optionsMenu->addMenuItem( "pre-filter by extension" );
  item->setCheckable( true );
  tooltip = "Eliminates common file types (ex. zip, doc, avi, etc) from search to speed results up.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  m_filterByExtension = item->checkBox();
  if( !m_filterByExtension ) //shouldnt ever happen
  {
    m_filterByExtension = new WCheckBox( "pre-filter by extension" );
    item->addWidget( m_filterByExtension );
  }
  
  m_filterByExtension->setChecked( dofilter );
  m_filterByExtension->checked().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_filterByExtension->unChecked().connect( this, &SpecFileQueryWidget::basePathChanged );
  
  
  item = m_optionsMenu->addMenuItem( "filter duplicate files" );
  item->setCheckable( true );
  tooltip = "Filters duplicate spectrum files through use of a hash of spectral and meta-information.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  m_filterUnique = item->checkBox();
  if( !m_filterUnique ) //shouldnt ever happen
  {
    m_filterUnique = new WCheckBox( "filter duplicate files" );
    item->addWidget( m_filterUnique );
  }
  
  m_filterUnique->setChecked( dofilter );
  m_filterUnique->checked().connect( this, &SpecFileQueryWidget::setResultsStale );
  m_filterUnique->unChecked().connect( this, &SpecFileQueryWidget::setResultsStale );
  
  
  linelayout->addWidget( new WContainerWidget(), 0, linelayout->columnCount() );
  linelayout->setColumnStretch( linelayout->columnCount()-1, 6 );
  
  
  item = m_optionsMenu->addMenuItem( "cache parse result" );
  item->setCheckable( true );
  tooltip = "Allows caching spectrum file parsing results for much faster subsequent searches of the same directory if you change search criteria.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  m_cacheParseResults = item->checkBox();
  if( !m_cacheParseResults ) //shouldnt ever happen
  {
    m_cacheParseResults = new WCheckBox( "Cache parse result" );
    item->addWidget( m_cacheParseResults );
  }
  
  m_cacheParseResults->setChecked( docache );
  m_cacheParseResults->checked().connect( this, &SpecFileQueryWidget::doCacheChanged );
  m_cacheParseResults->unChecked().connect( this, &SpecFileQueryWidget::doCacheChanged );
  
  
  
  item = m_optionsMenu->addMenuItem( "persist cache result" );
  item->setCheckable( true );
  tooltip = "Saves the results of parse caching to a file in the directory selected to be searched."
            " This allows the caching to persist across usages of this tool and InterSpec.";
  HelpSystem::attachToolTipOn( item, tooltip, instantToolTip, HelpSystem::Left );
  m_persistCacheResults = item->checkBox();
  if( !m_persistCacheResults ) //shouldnt ever happen
  {
    m_persistCacheResults = new WCheckBox( "Cache parse result" );
    item->addWidget( m_persistCacheResults );
  }
  
  m_persistCacheResults->setChecked( false );
  m_persistCacheResults->checked().connect( this, &SpecFileQueryWidget::doPersistCacheChanged );
  m_persistCacheResults->unChecked().connect( this, &SpecFileQueryWidget::doPersistCacheChanged );
  
  
  
  m_cancelUpdate = new WPushButton( "Cancel" );
  m_cancelUpdate->setHidden( true );
  linelayout->addWidget( m_cancelUpdate, 0, linelayout->columnCount() );
  m_cancelUpdate->clicked().connect( this, &SpecFileQueryWidget::cancelUpdate );
  
  m_update = new WPushButton( "Update" );
  m_update->clicked().connect(
  "function(){"
    "var result = $('#" + m_conditions->id() + "').queryBuilder('getRules',{ allow_invalid: true });"
    "var resultjson = (!$.isEmptyObject(result)) ? JSON.stringify(result, null, 2) : 'empty';"
    "Wt.emit('" + id() + "', 'fileSearchRequested', resultjson);"
  "}" );
  
  m_queryChanged.connect( boost::bind( &SpecFileQueryWidget::queryChangedCallback, this, _1 ) );
  m_searchRequested.connect( boost::bind( &SpecFileQueryWidget::searchRequestedCallback, this, _1 ) );
  
  m_update->disable();
  linelayout->addWidget( m_update, 0, linelayout->columnCount() );
  
  layout->addWidget( line, 2, 0 );
  
  m_resultview = new RowStretchTreeView();
  m_resultview->addStyleClass( "SpecFileQueryResults" );
  m_resultview->setRootIsDecorated( false ); //makes the tree look like a table! :)
  m_resultview->setAlternatingRowColors( true );
  m_resultview->setColumnWidth( 0, WLength(20,WLength::FontEx) );
  layout->addWidget( m_resultview, 3, 0 );
  layout->setRowStretch( 3, 5 );

  vector<string> collnames;
  for( const auto &filter : m_eventXmlFilters )
    collnames.push_back( filter.m_label );
  
  m_resultmodel = new ResultTableModel( this );
  m_resultmodel->setEventXmlColumns( collnames );
  
  m_resultview->setModel( m_resultmodel );
  m_resultview->setColumnHidden( DetectorName, true );
  m_resultview->setColumnHidden( Remark, true );
  m_resultview->setColumnHidden( LocationName, true );
  m_resultview->setColumnHidden( AnalysisResultText, true );
  m_resultview->setColumnHidden( HasRIIDAnalysis, true );
  
  for( size_t i = 0; i < m_eventXmlFilters.size(); ++i )
  {
    m_resultview->setColumnHidden( static_cast<int>(NumFileDataFields+i),
                                   !m_eventXmlFilters[i].show_in_gui_table );
  }
  
  m_resultview->selectionChanged().connect( this, &SpecFileQueryWidget::selectionChanged );
  m_resultview->setSelectionBehavior( Wt::SelectRows );
  m_resultview->setSelectionMode( Wt::SingleSelection );
  line = new WContainerWidget();
  linelayout = new WGridLayout( line );
  layout->addWidget( line, 4, 0 );
  
  m_numberFiles = new WText( "Select path to search" );
  linelayout->addWidget( m_numberFiles, 0, 0 );
  
  m_numberResults = new WText();
  linelayout->addWidget( m_numberResults, 0, 1 );
  
  linelayout->setColumnStretch( 0, 1 );
  linelayout->setColumnStretch( 1, 2 );
  
  
  m_csv = new WPushButton();
  ResultCsvResource *csvresource = new ResultCsvResource( m_resultmodel, m_csv );
  m_csv->setIcon( "InterSpec_resources/images/download_small.png" );
  m_csv->setLink( WLink(csvresource) );
  m_csv->setLinkTarget( Wt::TargetNewWindow );
  m_csv->setText( "CSV" );
  m_csv->setStyleClass( "CsvLinkBtn" );
  m_csv->disable();
  
  m_loadSelectedFile = new WPushButton( "Load Selected" );
  linelayout->addWidget( m_loadSelectedFile, 0, linelayout->columnCount(), AlignRight );
  m_loadSelectedFile->clicked().connect( this, &SpecFileQueryWidget::loadSelected );
  m_loadSelectedFile->disable();

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
#if( __APPLE__ )
  const char *show_msg = "Show in Finder...";
#else
  const char *show_msg = "Show in Explorer...";
#endif
  m_openSelectedDir = new WPushButton( show_msg );
  m_openSelectedDir->setToolTip( "Opens the directory containing the file in your operating systems file explorer." );
  linelayout->addWidget( m_openSelectedDir, 0, linelayout->columnCount(), AlignRight );
  m_openSelectedDir->clicked().connect( this, &SpecFileQueryWidget::openSelectedFilesParentDir );
  m_openSelectedDir->disable();
#endif

  
  linelayout->addWidget( m_csv, 0, linelayout->columnCount(), AlignBottom );
  
  const string init_js_call = "Wt.WT.FileQueryInit('" + m_conditions->id() + "','" + id() + "', " + addfilters + ");";
  
  doJavaScript( init_js_call );
  
#if( !BUILD_AS_OSX_APP )
  if( defpath.size() )
  {
    //TODO: Should do this basePathChanged in a non-gui thread.  Not sure if all
    //  the protections against a non-null database is set up, or if there could
    //  maybe be a race condition, so not niavely doing now.
    //  The reason fo rthis is if there is a persisted database, it can take a
    //  few seconds on the main GUI thread to do the initial sanity query.
    basePathChanged();
  }
#endif
}//void init()


std::string SpecFileQueryWidget::prepareEventXmlFilters()
{
  //An example to check on unique values in event XMLs:
  //grep -o "<urgency>.*</urgency>" ~/Desktop/spectra/Data/ExportScriptDownloaded/TE-18-*/*.xml | cut -d '>' -f2 | cut -d '<' -f1 > tmp.txt
  //sort < tmp.txt | uniq
  //
  //To recursively remove a file type from set of directories:
  //  find . -iname '*.kml' -exec rm {} \;
  //To list all file extensions recursively in a directory
  //  find . -type f | perl -ne 'print $1 if m/\.([^.\/]+)$/' | sort -u
  //To find the largest directories in the current directory
  //  du -a . | sort -n -r | head -n 20
  //To unzip all zip files in current directory, and recursively in each subdirectory
  //  find . -iname '*.zip' -exec sh -c 'unzip -d "${1%.*}" "$1"' _ {} \;
  
  
  string filters = "[";
  
  vector<string> jsonfilenames;
  try
  {
    const string jsonfilename = InterSpecUser::preferenceValue<string>("FileQueryEventQmlyFieldsFile", m_viewer);
    UtilityFunctions::split( jsonfilenames, jsonfilename, ";" );
  }catch(...)
  {
  }
  
  const string default_file = UtilityFunctions::append_path( InterSpec::staticDataDirectory(), "file_query_event_xml_fields.json" );
  bool hasDefault = false;
  for( const auto &f : jsonfilenames )
    hasDefault = hasDefault || (UtilityFunctions::filename(f) == "file_query_event_xml_fields.json");
  if( !hasDefault )
    jsonfilenames.push_back( default_file );
    
  
  for( const string jsonfilename : jsonfilenames )
  {
    ifstream inputjsonfile( jsonfilename.c_str() );
  
    if( !inputjsonfile.is_open() )
    {
      if( !jsonfilename.empty() && !UtilityFunctions::icontains(jsonfilename,"ouo") )
        passMessage( "Unable to open JSON file '" + jsonfilename + "' defining fields to search in Event XML.",
                   "", WarningWidget::WarningMsgHigh );
      continue;
    }//if( couldnt open file )
  
    const istream::pos_type orig_pos = inputjsonfile.tellg();
    inputjsonfile.seekg( 0, ios::end );
    const istream::pos_type eof_pos = inputjsonfile.tellg();
    inputjsonfile.seekg( orig_pos, ios::beg );
    const size_t filelen = 0 + eof_pos - orig_pos;
  
    if( filelen < 2 )
    {
      passMessage( "JSON file '" + jsonfilename + "' defining fields to search in Event XML is empty.",
                  "", WarningWidget::WarningMsgHigh );
      continue;
    }
  
    if( filelen > 256*1024 )
    {
      passMessage( "JSON file '" + jsonfilename + "' defining fields to search in Event XML is larger than allowed 256 kb - not using.",
                  "", WarningWidget::WarningMsgHigh );
      continue;
    }
  
    string jsonsourcestr( filelen, ' ' );
    inputjsonfile.read( &jsonsourcestr[0], filelen );
  
    //Allow a comment section at the beggining of the file, e.g., /*....*/
    UtilityFunctions::trim( jsonsourcestr );
    if( jsonsourcestr.size() > 4 && jsonsourcestr[0]=='/' && jsonsourcestr[0]=='*' )
    {
      auto pos = jsonsourcestr.find("*/");
      if( pos != string::npos )
        jsonsourcestr = jsonsourcestr.substr(pos+2);
    }
  
    std::vector<EventXmlFilterInfo> these_filters;
    try
    {
      these_filters = EventXmlFilterInfo::parseJsonString( jsonsourcestr );
    }catch( std::exception &e )
    {
      passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
      continue;
    }
  

    for( const auto &et : these_filters )
    {
      m_eventXmlFilters.push_back( et );
      
      string thisfilter;
      const string label = Wt::Utils::htmlEncode( et.m_label );
    
      thisfilter += (m_eventXmlFilters.size()>1 ? ",{" : "{");
      thisfilter += "\n\tid: '" + label + "'"
                  ",\n\tlabel: '" + label + "'";
    
      switch( et.m_type )
      {
        case EventXmlFilterInfo::InputType::Select:
        {
          thisfilter += ",\n\ttype: 'string'";
          thisfilter += ",\n\tinput: 'select'";
          thisfilter += ",\n\tvalues: [";
          for( size_t i = 0; i < et.m_discreet_options.size(); ++i )
            thisfilter += string(i?",":"") + "'" + et.m_discreet_options[i] + "'";
          thisfilter += "]";
          break;
        }//case EventXmlFilterInfo::InputType::Select:
        
        case EventXmlFilterInfo::InputType::Text:
        {
          thisfilter += ",\n\ttype: 'string'";
          if( !et.m_placeholder.empty() )
            thisfilter += ",\n\tplaceholder: '" + et.m_placeholder + "'";
          thisfilter += ",\n\tsize: 65";
          break;
        }
    
        case EventXmlFilterInfo::InputType::Date:
        {
          thisfilter += ",\n\ttype: 'date'";
          thisfilter += ",\n\tvalidation: {format:'YYYY/MM/DD'}";
          thisfilter += ",\n\tplaceholder: 'YYYY/MM/DD'";
        }//
      }//switch( et.second.m_type )
    
      thisfilter += ",\n\toperators: [";
      for( size_t i = 0; i < et.m_operators.size(); ++i )
        thisfilter += (i?",'":"'") + et.m_operators[i] + "'";
    
      thisfilter += "]";
      thisfilter += "\n}";
    
      filters += thisfilter;
    }//for( const auto &et : m_eventXmlFilters )
  }//for( const string jsonfilename : jsonfilenames )
    
  filters += "]";
  
  return filters;
}//std::string prepareEventXmlFilters()


void SpecFileQueryWidget::setResultsStale()
{
  m_stopUpdate->store( true );
  
  if( m_resultview->selectedIndexes().size() )
    m_resultview->setSelectedIndexes( WModelIndexSet() );
  
  if( m_resultmodel->rowCount() )
    m_resultmodel->setResult( std::shared_ptr< vector< vector<string> > >() );
    //m_resultmodel->removeRows( 0, m_resultmodel->rowCount() );
  
  m_resultview->disable();
  m_csv->disable();
  
  m_optionsBtn->enable();
  m_optionsMenu->enable();
  
  if( m_openSelectedDir )
    m_openSelectedDir->disable();
  
  if( m_loadSelectedFile )
    m_loadSelectedFile->disable();
  
  //Validate logic, and give hint if something wrong, or else
  bool queryvalid = true;
  string errormsg;
  try
  {
    queryFromJson( m_queryJson, m_eventXmlFilters );
  }catch( std::exception &e )
  {
    errormsg = e.what();
    queryvalid = false;
  }
  
  if( queryvalid )
  {
    m_numberResults->removeStyleClass( "SpecFileQueryErrorTxt" );
    m_update->enable();
    m_numberResults->setText( "Tap Update to refresh results" );
  }else
  {
    m_numberResults->setText( errormsg );
    m_numberResults->addStyleClass( "SpecFileQueryErrorTxt" );
  }
}//void setResultsStale()



void SpecFileQueryWidget::doCacheChanged()
{
  const std::string basepath = baseDirectory();
  
  auto map_iter = m_path_caches.find( basepath );
  if( map_iter != end(m_path_caches) )
  {
    const bool docache = m_cacheParseResults->isChecked();
    
    m_persistCacheResults->setUnChecked();
    m_persistCacheResults->setDisabled( !docache );
    
    std::shared_ptr<SpecFileQueryDbCache> database = map_iter->second;
    if( database && (database->caching_enabled()==docache) )
      return;//probably shouldnt ever happen, but whatever.
    
    m_path_caches.erase(map_iter);
    if( !docache )
      database->stop_caching();
  }//if( map_iter != end(m_path_caches) )
  
  basePathChanged();
}//void doCacheChanged()


void SpecFileQueryWidget::doPersistCacheChanged()
{
  const bool persist = m_persistCacheResults->isChecked();
  const std::string basepath = baseDirectory();
  auto map_iter = m_path_caches.find( basepath );
  if( (map_iter != end(m_path_caches)) && map_iter->second )
  {
    try
    {
      map_iter->second->set_persist( persist );
    }catch( std::exception &e )
    {
      passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
      m_cacheParseResults->setChecked( false );
      m_persistCacheResults->setChecked( false );
      m_persistCacheResults->disable();
    }
  }//if( we can find the correct cache )
}//void doPersistCacheChanged();


#if( BUILD_AS_ELECTRON_APP )
void SpecFileQueryWidget::newElectronPathSelected( std::string path )
{
  m_basePath = path;
  basePathChanged();
}
#elif( BUILD_AS_OSX_APP )
void SpecFileQueryWidget::newMacOsPathSelected()
{
  m_basePath = SpecFileQuery::getSearchDirectory();
  basePathChanged();
}//void newMacOsPathSelected()
#endif


void SpecFileQueryWidget::basePathChanged()
{
  setResultsStale();
  
  m_numberFiles->setText( "" );
  m_numberResults->setText( "" );
  m_update->enable();


  const std::string basepath = baseDirectory();

  if( basepath.empty() || !UtilityFunctions::is_directory( basepath ) )
  {
    m_numberResults->setText( "Not a valid base directory" );
#if( !BUILD_AS_OSX_APP && !BUILD_AS_ELECTRON_APP )
    m_baseLocation->addStyleClass( "SpecFileQueryErrorTxt" );
#endif
    m_update->disable();
    return;
  }
  
  const bool recursive = m_recursive->isChecked();
  const bool docache = m_cacheParseResults->isChecked();
  const bool filter = m_filterByExtension->isChecked();
  const bool filterUnique = m_filterUnique->isChecked();
  const int maxsize = ((m_maxFileSize->value()==0) ? 2*1024 : m_maxFileSize->value());
  
  string prefpath;
  int prevmaxsize = 32;
  bool prevrecursive = true, prevfilter = true, prevFilterUnique = true, prevcache = true;
  try
  {
    prevmaxsize = InterSpecUser::preferenceValue<int>( "SpecFileQueryMaxSize", m_viewer );
    prevcache = InterSpecUser::preferenceValue<bool>( "SpecFileQueryCacheParse", m_viewer );
    prevfilter = InterSpecUser::preferenceValue<bool>( "SpecFileQueryFilter", m_viewer );
    prevrecursive = InterSpecUser::preferenceValue<bool>( "SpecFileQueryRecursive", m_viewer );
    prevFilterUnique = InterSpecUser::preferenceValue<bool>( "SpecFileQueryUnique", m_viewer );
#if( BUILD_AS_OSX_APP || BUILD_AS_ELECTRON_APP )
    prefpath = "";
#else
    prefpath = InterSpecUser::preferenceValue<string>( "SpecFileQueryPath", m_viewer );
    if( prefpath == "None" )
      prefpath = "";
#endif
  }catch(...)
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected error grabbing prefefence value" );
#endif
  }
  
  try
  {
    if( prevrecursive != recursive )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryRecursive", recursive, m_viewer );
    if( prevcache != docache )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryCacheParse", docache, m_viewer );
    if( prevfilter != filter )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryFilter", filter, m_viewer );
    if( maxsize != prevmaxsize )
      InterSpecUser::setPreferenceValue<int>( m_viewer->m_user, "SpecFileQueryFilter", maxsize, m_viewer );
    if( filterUnique != prevFilterUnique )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryUnique", filterUnique, m_viewer );
#if( !BUILD_AS_OSX_APP && !BUILD_AS_ELECTRON_APP )
    if( prefpath != basepath )
      InterSpecUser::setPreferenceValue<string>( m_viewer->m_user, "SpecFileQueryPath", basepath, m_viewer );
#endif
  }catch( ... )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( BOOST_CURRENT_FUNCTION, "Unexpected error setting prefefence value" );
#endif
  }
  
#if( !BUILD_AS_OSX_APP  && !BUILD_AS_ELECTRON_APP )
  m_baseLocation->removeStyleClass( "SpecFileQueryErrorTxt" );
#endif
  
  for( auto &i : m_path_caches )
  {
    if( i.second )
      i.second->stop_caching();
  }
  
  const bool cache_in_db = m_cacheParseResults->isChecked();
  
  std::shared_ptr<SpecFileQueryDbCache> database;
  auto map_iter = m_path_caches.find( basepath );
  if( map_iter != end(m_path_caches) )
    database = map_iter->second;
  if( !database )
    database = m_path_caches[basepath] = make_shared<SpecFileQueryDbCache>( cache_in_db, basepath, m_eventXmlFilters );
  
  if( cache_in_db )
  {
    database->allow_start_caching();
    m_persistCacheResults->enable();
    m_persistCacheResults->setChecked( database->is_persistand() );
  }else
  {
    m_persistCacheResults->disable();
    m_persistCacheResults->setChecked( false );
  }
  
  WServer *server = WServer::instance();  //can this ever be NULL?
  if( server )
  {
    m_numberFiles->setText( "Updating # of files" );
    const size_t maxsize_mb = static_cast<size_t>(maxsize*1024*1024);
    server->ioService().post( boost::bind( &SpecFileQueryWidget::updateNumberFiles,
                                           basepath, recursive, filter, maxsize_mb,
                                           this, wApp->sessionId(), m_widgetDeleted,
                                           database ) );
  }else
  {
    cerr << SRC_LOCATION << "\n\tWarning: couldnt get WServer to post to - not good" << endl << endl;
  }//if( server ) / else
}//void basePathChanged()



void SpecFileQueryWidget::updateNumberFiles( const string srcdir,
                                    const bool recursive, const bool extfilter,
                                    const size_t maxsize,
                                    SpecFileQueryWidget *querywidget,
                                    const std::string sessionid,
                                    std::shared_ptr< std::atomic<bool> > widgetdeleted,
                                    std::shared_ptr<SpecFileQueryDbCache> database )
{
  vector<string> files;
  
  UtilityFunctions::file_match_function_t filterfcn = extfilter ? &maybe_spec_file : &file_smaller_than;
  
#if( USE_DIRECTORY_ITERATOR_METHOD )
  double updatetime = UtilityFunctions::get_wall_time();
  const bool docache = (database && database->caching_enabled());
  
  size_t nfiles = 0;
  boost::filesystem::recursive_directory_iterator diriter( srcdir, boost::filesystem::symlink_option::recurse );
  const boost::filesystem::recursive_directory_iterator dirend;
  
  while( diriter != dirend )
  {
    if( *widgetdeleted )
      return;
    
    string filename = diriter->path().string<std::string>();
    
    const bool is_dir = boost::filesystem::is_directory( diriter.status() ); //folows symlinks to see if target of symlink is a directory
    bool is_file = (diriter.status().type() == boost::filesystem::file_type::regular_file);  //folows symlinks
    
    if( !recursive && is_dir )
    {
      is_file = false; //JIC
      diriter.no_push();  //Dont recurse down into directories if we arent doing a recursive search
    }
    
    bool is_simlink_dir = false;
    if( is_dir && recursive )
    {
      //If this is a directory, check if we are actually on a symlink to a
      //  directory, because if so, we need to check for cyclical links.
      boost::system::error_code symec;
      const auto symstat = boost::filesystem::symlink_status( diriter->path(), symec );
      is_simlink_dir = (!symec && (symstat.type()==boost::filesystem::file_type::symlink_file));
    }
    
    if( is_simlink_dir )
    {
      auto resvedpath = boost::filesystem::read_symlink( diriter->path() );
      if( resvedpath.is_relative() )
        resvedpath = diriter->path().parent_path() / resvedpath;
      resvedpath = boost::filesystem::canonical( resvedpath );
      auto pcanon = boost::filesystem::canonical( diriter->path().parent_path() );
      if( UtilityFunctions::starts_with( pcanon.string<string>(), resvedpath.string<string>().c_str() ) )
        diriter.no_push();  //Dont recurse down into directories
    }//if( is_simlink_dir && recursive )
    
    
    if( is_file && filterfcn( filename, (void *)&maxsize ) )
    {
      if( docache && (files.size() < 100000) )
        files.push_back( filename );
      
      ++nfiles;
      if( (nfiles % 500) == 0 )
      {
        const double nowtime = UtilityFunctions::get_wall_time();
        if( (nowtime - updatetime) > 1.0 )
        {
          updatetime = nowtime;
          if( !(*widgetdeleted) )
            WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateNumberFilesInGui,
                                                              nfiles, false, srcdir, recursive, extfilter, querywidget, widgetdeleted ) );
          else
            return;
        }
      }//if( (nfiles % 500) == 0 )
    }//if( this is a potential file we should check on )
    
    boost::system::error_code ec;
    diriter.increment(ec);
    while( ec && (diriter!=dirend) )
    {
      std::cerr << "Error While Accessing : " << diriter->path().string() << " :: " << ec.message() << '\n';
      diriter.increment(ec);
    }
  }//while( diriter != dirend )
  
  if( !(*widgetdeleted) )
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateNumberFilesInGui,
                                                      nfiles, true, srcdir, recursive, extfilter, querywidget, widgetdeleted ) );
#else
  if( recursive )
    files = UtilityFunctions::recursive_ls( srcdir, filterfcn, (void *)&maxsize );
  else
    files = UtilityFunctions::ls_files_in_directory( srcdir, filterfcn, (void *)&maxsize );
  
  const size_t nfiles = files.size();
  
  if( !(*widgetdeleted) )
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateNumberFilesInGui,
                                                    nfiles, true, srcdir, recursive, extfilter, querywidget, widgetdeleted ) );
  
#endif
  
  
  if( database && database->caching_enabled() )
    database->cache_results( std::move(files) );
}//updateNumberFiles


std::string SpecFileQueryWidget::baseDirectory()
{
#if( BUILD_AS_OSX_APP || BUILD_AS_ELECTRON_APP )
  return m_basePath;
#else
  return m_baseLocation->text().toUTF8();
#endif
}//std::string baseDirectory()


void SpecFileQueryWidget::updateNumberFilesInGui( const size_t nfiles,
                            const bool completed,
                            const std::string srcdir,
                            const bool recursive,
                            const bool extfilter,
                            SpecFileQueryWidget *querywidget,
                            std::shared_ptr< std::atomic<bool> > widgetdeleted )
{
  //We should be operating in our main WApplication thread here
  assert( wApp );  //lets double check that.
  
  if( (*widgetdeleted) )
    return;
  
  const string newdir = querywidget->baseDirectory();
  
  //Check if the user has changed pre-filter chriteria since this search was
  //  launched, and if so dont update the GUI (another search has been launched
  //  and will update things)
  if( newdir != srcdir
      || newdir.empty()
      || recursive != querywidget->m_recursive->isChecked()
      || extfilter != querywidget->m_filterByExtension->isChecked() )
      return;
  
#if( USE_DIRECTORY_ITERATOR_METHOD )
  if( completed )
  {
    querywidget->m_numberFiles->setText( (to_string(nfiles) + " initial files").c_str() );
  }else
  {
    string msg = to_string(nfiles) + "+ files...";
    querywidget->m_numberFiles->setText( msg.c_str() );
  }
#else
  const string nfilesstr = (nfiles<99000 ? to_string(nfiles) : string(">100k"));
  querywidget->m_numberFiles->setText( (nfilesstr + " initial files").c_str() );
#endif //USE_DIRECTORY_ITERATOR_METHOD
  
  wApp->triggerUpdate();
}//updateNumberFilesInGui(...)


void SpecFileQueryWidget::setResultFieldVisibility( const SpecFileQuery::FileDataField field, const bool visible )
{
  m_resultview->setColumnHidden( field, visible );
}//void setResultFieldVisibility(...)


void SpecFileQueryWidget::queryChangedCallback( const std::string &queryJson )
{
  //should have been checked client side that the query has actually changed,
  //  so we'll skip doing that here.
  m_queryJson = queryJson;
  
  setResultsStale();
}//void queryChangedCallback( const std::string &queryJson );


SpecFileQuery::SpecLogicTest SpecFileQueryWidget::queryFromJson( const std::string &queryJson,
                const std::vector<EventXmlFilterInfo> &eventXml )
{
  SpecFileQuery::SpecLogicTest test;
  
  try
  {
    if( UtilityFunctions::istarts_with(queryJson, "empty") )
      throw runtime_error( "<p>Unable to translate query into test logic.</p>"
                            "<p>Please take a screenshot of query and send to "
                            "<a href=\"mailto:interspec@sandia.gov\" target=\"_blank\">interspec@sandia.gov</a></p>" );
  
    Json::Object q;
    Wt::Json::parse( queryJson, q );
  
    const auto &valid = q.get("valid");
  
    if( valid.isNull() )
      throw std::runtime_error( "Query logic is not formated as expected." );
  
    if( !valid )
      throw std::runtime_error( "Query is currently not valid, please fix errors." );
  
    const Json::Value &baseCond = q.get( "condition" );
    const Json::Array &baseRules = q.get( "rules" );
  
    if( baseCond.isNull() || baseRules==Json::Array::Empty ) //shouldnt ever happen
      throw std::runtime_error( "Error converting query to search logic: case condition or rules not avaialble" );
  
    add_logic( test, baseCond, baseRules, eventXml );
  
    test.isvalid();  //throws exception giving description of issue
  }catch( Wt::Json::ParseError &e )
  {
    throw runtime_error( "Error parsing JSON from user interface: " + string(e.what()) );
  }catch( std::runtime_error & )
  {
    throw;
  }//try /catch
  
  return test;
}//SpecFileQuery::SpecLogicTest SpecFileQueryWidget::queryFromJson( const std::string &json )


void SpecFileQueryWidget::searchRequestedCallback( const std::string &queryJson )
{
  //cout << "Got query string: " << queryJson << endl;
  m_queryJson = queryJson;
  
  try
  {
    SpecFileQuery::SpecLogicTest test = queryFromJson( queryJson, m_eventXmlFilters );
    
    //cout << "Translated JSON to: " << test.summary() << endl;
    
    const std::string basepath = baseDirectory();
    const bool recursive = m_recursive->isChecked();
    const bool filter = m_filterByExtension->isChecked();
    const size_t maxsize = 1024*1024*(m_maxFileSize->value() ? m_maxFileSize->value() : 2*1024);
    const bool filterUnique = m_filterUnique->isChecked();
    
    m_cancelUpdate->show();
    m_update->hide();
#if( BUILD_AS_ELECTRON_APP )
#else
    m_baseLocation->disable();
#endif
    m_conditions->disable();
    
    //ToDo: need to make it so you cant modify any of the QUeryBuilder widgets until the search is ended.  Maybe something like:
    //var rule = $('#builder').queryBuilder('getRules', { get_flags: true });
    // Then set 'readonly' flag to true
    //$('#builder').queryBuilder('setRules', rule);
    //  but for right now we'll just take the heavyhanded approach and manually siable all the widget.
    doJavaScript( "$('#" + m_conditions->id() + "').find('select, input, button, .btn').prop('disabled',true).addClass('Wt-disabled');" );
    
    m_csv->disable();
    m_optionsBtn->disable();
    m_optionsMenu->disable();
    m_resultview->disable();
    m_numberResults->setText( "Updating results" );
    if( m_resultmodel->rowCount() )
      m_resultmodel->setResult( nullptr );
    
    unsigned long options = 0x0;
    if( recursive )
      options |= SearchRecursive;
    if( filter )
      options |= FilterByFilename;
    if( filterUnique )
      options |= FilterDuplicates;
    
    std::shared_ptr<SpecFileQueryDbCache> database;
    auto map_iter = m_path_caches.find( basepath );
    if( map_iter != end(m_path_caches) )
      database = map_iter->second;
    
    if( !database )
      database = std::make_shared<SpecFileQueryDbCache>( false, basepath, m_eventXmlFilters );
    
    m_stopUpdate->store( false );
    WServer::instance()->ioService().post(
                                          boost::bind( &SpecFileQueryWidget::doSearch, this, basepath, options,
                                                      maxsize, test, wApp->sessionId(), database, m_stopUpdate,
                                                      m_widgetDeleted ) );
  }catch( std::runtime_error &e )
  {
    passMessage( "Current query is invalid: " + string(e.what()), "", WarningWidget::WarningMsgHigh );
  }//try /catch
}//void searchRequestedCallback( std::string queryJson )


void SpecFileQueryWidget::finishUpdate( std::shared_ptr< std::vector< std::vector<std::string> > > result,
                                        const std::string description,
                                        const double clock_seconds,
                                        const bool wasCanceled,
                                        std::shared_ptr< std::atomic<bool> > widgetDeleted )
{
  if( widgetDeleted->load() )
    return;
  
  m_cancelUpdate->hide();
  m_update->show();
#if( BUILD_AS_ELECTRON_APP )
#else
  m_baseLocation->enable();
#endif
  m_conditions->enable();
  m_resultview->enable();
  m_csv->enable();
  m_optionsBtn->enable();
  m_optionsMenu->enable();
  doJavaScript( "$('#" + m_conditions->id() + "').find('select, input, button, .btn').prop('disabled',false).removeClass('Wt-disabled');" );
  
  if( wasCanceled )
  {
    m_csv->disable();
    m_numberResults->setText( "Search Canceled" );
    setResultsStale();
    return;
  }
  
  m_resultmodel->setHeaderData( 0, Horizontal, boost::any(WString(description)), Wt::UserRole );
  //m_resultmodel->setResult( result );
  m_resultmodel->appendResults( result );
  
  const string timestr = PhysicalUnits::printToBestTimeUnits( clock_seconds, 1, 1.0 );

  const int numRows = m_resultmodel->rowCount();
  m_numberResults->setText( std::to_string(numRows) + " matching files. Took " + timestr + "." );
  
  wApp->triggerUpdate();
}//finishUpdate(...)


void SpecFileQueryWidget::cancelUpdate()
{
  m_stopUpdate->store( true );
  setResultsStale();
}//cancelUpdate(...)


void SpecFileQueryWidget::updateSearchStatus( const size_t nfilestotal,
                                              const size_t nfileschecked,
                                              const std::string specialmsg,
                                              std::shared_ptr< std::vector< std::vector<std::string> > > results,
                                              std::shared_ptr< std::atomic<bool> > widgetDeleted )
{
  if( widgetDeleted->load() )
    return;
  
  if( results && results->size() )
    m_resultmodel->appendResults( results );
  
  const int nfilesaccepted = m_resultmodel->rowCount();
  
  char buffer[256];
  if( nfilestotal == 0 )
  {
    snprintf( buffer, sizeof(buffer), "Checked %i files with %i passing%s",
             static_cast<int>(nfileschecked), nfilesaccepted, specialmsg.c_str() );
  }else
  {
    snprintf( buffer, sizeof(buffer), "%i of %i checked with %i passing%s",
              static_cast<int>(nfileschecked), static_cast<int>(nfilestotal),
              nfilesaccepted, specialmsg.c_str() );
  }
  m_numberResults->setText( buffer );
  
  
  wApp->triggerUpdate();
}//void updateSearchStatus(...)



void testfile( const string &filename, vector<string> &result,
               const SpecFileQuery::SpecLogicTest &query,
               HaveSeenUuid &uniquecheck,
               std::shared_ptr< SpecFileQueryDbCache > database,
               const string &base_search_dir )
{
  try
  {
    result.clear();
    
    std::unique_ptr<SpecFileInfoToQuery> db_test_info = database->spec_file_info( filename );
    
    if( !db_test_info || !db_test_info->is_spectrum_file /*&& !db_test_info.is_event_xml_file*/ )
      return;
    
    if( uniquecheck.have_seen( db_test_info->uuid ) )
      return;
    
    const bool testresult = query.test( *db_test_info );
    
    if( testresult )
      result = get_result_fields( *db_test_info, base_search_dir, database->xml_filters() );
  }catch( ... )
  {
    cerr << "Caught exception testing file: " << filename << endl;
  }
}//testfile



void SpecFileQueryWidget::doSearch( const std::string basedir,
                                    unsigned long options,
                                    const size_t maxsize,
                                    const SpecFileQuery::SpecLogicTest query,
                                    const std::string sessionid,
                                    std::shared_ptr< SpecFileQueryDbCache > database,
                                    std::shared_ptr< std::atomic<bool> > stopUpdate,
                                    std::shared_ptr< std::atomic<bool> > widgetDeleted )
{
  std::shared_ptr< vector<vector<string> > > result;
 
  const bool recursive = (options & SearchRecursive);
  const bool extfilter = (options & FilterByFilename);
  const bool filterUnique = (options & FilterDuplicates);
  
  const double starttime = UtilityFunctions::get_wall_time();
  const double startcpu = UtilityFunctions::get_cpu_time();
  
  stringstream description;
  boost::posix_time::ptime utctime = boost::posix_time::second_clock::universal_time();
  boost::posix_time::ptime localtime = boost::posix_time::second_clock::local_time();
  description << "Search performed at " << UtilityFunctions::to_extended_iso_string( localtime )
              << " local (" << UtilityFunctions::to_extended_iso_string( utctime )
              << " UTC)\r\n";
  if( recursive )
    description << "Recursive" << " search of \"" << basedir << "\"\r\n";
  else
   description << "Non-recursive" << " search of \"" << basedir << "\"\r\n";

  if( extfilter )
    description << "Candidate files pre-filtered based on file extension and other niave properties to exclude common non-spectrum files\r\n";
  description << "Candidate files pre-filtered to have a size under " << (maxsize/1024/1024) << " MB\r\n";
  if( filterUnique )
    description << "Only considered unique spectrum files.\r\n";
  description << "Search Query: " << query.summary() << "\r\n\r\n";
  
  HaveSeenUuid uniqueCheck( filterUnique );
  
  //If we are currently caching files to the database, we should stop that
  //  so we dont end up stepping on eachothers toes when parsing the same file
  //  twice (which happens when the search catches up to the caching).
  if( database )
    database->stop_caching();
  
  size_t num_files_pass = 0;
  
  try
  {
    if( stopUpdate->load() )
      throw std::runtime_error( "" );
    
    result = std::make_shared< vector<vector<string> > >();
    
    
    typedef vector<string>(*ls_fcn_t)( const string &, UtilityFunctions::file_match_function_t, void * );
    UtilityFunctions::file_match_function_t filterfcn = extfilter ? &maybe_spec_file : &file_smaller_than;
    
    
#if( USE_DIRECTORY_ITERATOR_METHOD )
    if( stopUpdate->load() )
      throw std::runtime_error( "" );
    
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                      this, 0, 0, "", result, widgetDeleted ) );
#else
    ls_fcn_t lsfcn = &UtilityFunctions::recursive_ls;
    if( !recursive )
      lsfcn = &UtilityFunctions::ls_files_in_directory;
    
    const vector<string> files = lsfcn( basedir, filterfcn, (void *)&maxsize );
    const int nfiles = static_cast<int>( files.size() );
    
    description << "There were " << nfiles << " candidate files after pre-filtering\r\n";
    
    if( stopUpdate->load() )
      throw std::runtime_error( "" );
    
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                      this, nfiles, 0, "", result, widgetDeleted ) );
#endif //USE_DIRECTORY_ITERATOR_METHOD
    
    int nupdates_sent = 0;
    double lastupdate = UtilityFunctions::get_wall_time();
    
    //We could choose something like 100 for nfile_at_a_time with the only
    //  dowsides being N42 files would parse any faster due to not enough
    //  threads being available, and also the fact that when parsing multiple
    //  file, hard drive bandwidth is probably the main limiting factor.
#if( defined(WIN32) )
    //The penalty of multiple seeks is redicuolous on spinning drives - just use a single thread...
    const int nfile_at_a_time = 1;
#else
    const int nfile_at_a_time = SpecUtilsAsync::num_physical_cpu_cores();
#endif
    
    SpecUtilsAsync::ThreadPool pool;
    
#if( USE_DIRECTORY_ITERATOR_METHOD )
    size_t ncheckssubmitted = 0;
    std::mutex result_mutex;
    
    boost::filesystem::recursive_directory_iterator diriter( basedir, boost::filesystem::symlink_option::recurse );
    const boost::filesystem::recursive_directory_iterator dirend;
    
    while( diriter != dirend )
    {
      if( stopUpdate->load() )
        throw runtime_error("");
      
      string filename = diriter->path().string<std::string>();
      
      const bool is_dir = boost::filesystem::is_directory( diriter.status() ); //folows symlinks to see if target of symlink is a directory
      bool is_file = (diriter.status().type() == boost::filesystem::file_type::regular_file);  //folows symlinks
     
      if( !recursive && is_dir )
      {
        is_file = false; //JIC
        diriter.no_push();  //Dont recurse down into directories if we arent doing a recursive search
      }
      
      bool is_simlink_dir = false;
      if( is_dir && recursive )
      {
        //If this is a directory, check if we are actually on a symlink to a
        //  directory, because if so, we need to check for cyclical links.
        boost::system::error_code symec;
        const auto symstat = boost::filesystem::symlink_status( diriter->path(), symec );
        is_simlink_dir = (!symec && (symstat.type()==boost::filesystem::file_type::symlink_file));
      }
      
      if( is_simlink_dir )
      {
        auto resvedpath = boost::filesystem::read_symlink( diriter->path() );
        if( resvedpath.is_relative() )
          resvedpath = diriter->path().parent_path() / resvedpath;
        resvedpath = boost::filesystem::canonical( resvedpath );
        auto pcanon = boost::filesystem::canonical( diriter->path().parent_path() );
        if( UtilityFunctions::starts_with( pcanon.string<string>(), resvedpath.string<string>().c_str() ) )
          diriter.no_push();  //Dont recurse down into directories
      }//if( is_simlink_dir && recursive )
      
      
      if( is_file && filterfcn( filename, (void *)&maxsize ) )
      {
        if( (ncheckssubmitted % nfile_at_a_time) == 0 )
        {
          pool.join();
          
          const double now = UtilityFunctions::get_wall_time();
          if( now > (lastupdate + 1.0) || !nupdates_sent )
          {
            ++nupdates_sent;
            num_files_pass += result->size();
            
            std::unique_lock<std::mutex> lock( result_mutex );
            WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                              this, 0, ncheckssubmitted, "", result, widgetDeleted ) );
            result = std::make_shared< vector<vector<string> > >();
            lastupdate = now;
          }
        }
        
        pool.post( [filename,&query,&uniqueCheck,&database,&basedir,&result,&result_mutex](){
          vector<string> testres;
          testfile( filename, testres, query, uniqueCheck, database, basedir );
          if( !testres.empty() )
          {
            std::unique_lock<std::mutex> lock( result_mutex );
            result->push_back( testres );
          }
        } );
        
        ++ncheckssubmitted;
      }//if( this is a potential file we should check on )
      
      boost::system::error_code ec;
      diriter.increment(ec);
      while( ec && (diriter!=dirend) )
      {
        std::cerr << "Error While Accessing : " << diriter->path().string() << " :: " << ec.message() << '\n';
        diriter.increment(ec);
      }
    }//while( diriter != dirend )
    
    pool.join();
    
#else
    
    for( int i = 0; i < nfiles; i += nfile_at_a_time )
    {
      if( stopUpdate->load() )
        throw runtime_error("");
      
      const int nfilethisone = std::min( nfiles-i, nfile_at_a_time );
      
      vector< vector<string> > testres( nfilethisone );
      
      for( int j = 0; j < nfilethisone; ++j )
        pool.post( boost::bind( &testfile, boost::cref(files[i+j]), boost::ref(testres[j]),
                                boost::cref(query), boost::ref(uniqueCheck), database,
                                boost::cref(basedir) ) );
      pool.join();
      
      for( int j = 0; j < nfilethisone; ++j )
        if( testres[j].size() )
          result->push_back( testres[j] );
      
      const double now = UtilityFunctions::get_wall_time();
      if( now > (lastupdate + 2.0) || !nupdates_sent )
      {
        ++nupdates_sent;
        num_files_pass += result->size();
        WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                          this, nfiles, i+nfilethisone, "", result, widgetDeleted ) );
        result = std::make_shared< vector<vector<string> > >();
        lastupdate = now;
      }
    }//for( size_t i = 0; i < nfiles; ++i  )
#endif //USE_DIRECTORY_ITERATOR_METHOD
  }catch( ... )
  {
    const double total_clock_time = (UtilityFunctions::get_wall_time() - starttime);
    WServer::instance()->post( sessionid, boost::bind( &SpecFileQueryWidget::finishUpdate, this, result, description.str(), total_clock_time, true, widgetDeleted ) );
  }
  
  const double total_clock_time = (UtilityFunctions::get_wall_time() - starttime);
  const double finishtime = UtilityFunctions::get_wall_time();
  const double finishcpu = UtilityFunctions::get_cpu_time();
  
  description << "Search took " << (finishtime - starttime) << " wall seconds and " << (finishcpu - startcpu) << " cpu seconds \r\n";
  description << "There were " << (num_files_pass+result->size()) << " files that satisfied the test conditions\r\n";
  
  //ToDo: I guess there is, in principle, a slight chance that the previous post
  //      to updateSearchStatus() could happen _after_ this next post to
  //      finishUpdate().  But I guess life will go on for now (I dont think this
  //      will actually happen because although its not gaurnteed in the future,
  //      the underlying boost::asio scedules jobs serially - I think).
  WServer::instance()->post( sessionid, boost::bind( &SpecFileQueryWidget::finishUpdate, this, result, description.str(), total_clock_time,  false, widgetDeleted ) );
}//doSearch(...)
  


void SpecFileQueryWidget::selectionChanged()
{
  if( !m_viewer )
    return;
  
  const bool sel = (m_resultview->selectedIndexes().size() == 1);
  if( m_loadSelectedFile )
    m_loadSelectedFile->setEnabled( sel );
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  if( m_openSelectedDir )
    m_openSelectedDir->setEnabled( sel );
#endif
}//selectionChanged()



#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
void SpecFileQueryWidget::openSelectedFilesParentDir()
{
  const WModelIndexSet selected = m_resultview->selectedIndexes();
  if(selected.size() != 1 )
    return;
  
  const int row = selected.begin()->row();
  WModelIndex findex = m_resultmodel->index( row, SpecFileQuery::Filename );
  boost::any fnany = m_resultmodel->data( findex, Wt::UserRole );
  std::string filenameandy = asString(fnany).toUTF8();
  
  const string parentdir = UtilityFunctions::parent_path(filenameandy);
  if( !UtilityFunctions::is_directory(parentdir) )
    return;
  
#if( BUILD_AS_ELECTRON_APP )
  //Opens the directory in the operating system
  //doJavaScript( "const {shell} = require('electron'); shell.openItem('" + parentdir + "');" );  //Works on macOS.  Need to test on WIndows.
  
  //Opens the directory and highlights file in oeprating system file manager
  //For windows need to ecape '\'.

#if( defined(WIN32) )
  //UtilityFunctions::ireplace_all( filenameandy, "\\\\", "[%%%%%%%%]" );
  UtilityFunctions::ireplace_all( filenameandy, "\\", "\\\\" );
  //UtilityFunctions::ireplace_all( filenameandy, "[%%%%%%%%]", "\\\\\\\\" );
  //ToDo: Check into compilation optimization 
  //      Add time it took for a search to result text
#endif

  doJavaScript( "const {shell} = require('electron'); shell.showItemInFolder('" + filenameandy + "');" );  //Works on macOS.  Need to test on WIndows.
#elif( BUILD_AS_LOCAL_SERVER || __APPLE__ )
#if( __APPLE__ )
  const string command = "open '" + parentdir + "'";
  system( command.c_str() );  //Havent tested with macOS Sandboxing yet.
#elif( BUILD_AS_LOCAL_SERVER && defined(WIN32) )
  ShellExecute(NULL, NULL, parentdir.c_str(), NULL, NULL,  SW_SHOWNORMAL);
  //static_assert( 0, "You need to implement opening a Explorer window or something" );
#else
  static_assert( 0, "You need to implement opening a Finder window or something" );
#endif
#endif
}//void SpecFileQueryWidget::openSelectedFilesParentDir()
#endif  //#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )


void SpecFileQueryWidget::loadSelected()
{
  if( !m_viewer || !m_viewer->fileManager() )
    return;
  
  const WModelIndexSet selected = m_resultview->selectedIndexes();
  if(selected.size() != 1 )
    return;
  
  const int row = selected.begin()->row();
  WModelIndex findex = m_resultmodel->index( row, SpecFileQuery::Filename );
  boost::any fnany = m_resultmodel->data( findex, Wt::UserRole );
  const std::string filenameandy = asString(fnany).toUTF8();
  
  if( !UtilityFunctions::is_file(filenameandy) )
  {
    passMessage( "Sorry, '" + filenameandy + "' file appears to no longer be accessable", "", WarningWidget::WarningMsgHigh );
    m_resultview->setSelectedIndexes( WModelIndexSet() );
    return;
  }
  
  
  const bool loaded = m_viewer->fileManager()->loadFromFileSystem( filenameandy, kForeground, kAutoParser );
  
  if( !loaded )
    passMessage( "Sorry, '" + filenameandy + "' can not be displayed", "", WarningWidget::WarningMsgHigh );
}//loadSelected()
  
