/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <cerrno>
#include <limits>
#include <sstream>

#include <boost/regex.hpp>

#include <Wt/WText>
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
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/SpecFileQueryWidget.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace Wt;
using namespace std;

namespace SpecQuery
{
  
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
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
        break;
        
      case DetectionSystemType:
      case SearchMode:
      case TotalLiveTime:
      case TotalRealTime:
      case IndividualSpectrumLiveTime:
      case IndividualSpectrumRealTime:
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
        break;

      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
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
        break;
        
      case DetectionSystemType: case SearchMode:
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
      case StartTime: case NumFileDataFields:
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
      case StartTime:
        break;
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
      case DetectionSystemType: case SearchMode:
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
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
      
      case TextIsContained:
        return UtilityFunctions::icontains( teststr, ss );
        
      case TextStartsWith:
        return UtilityFunctions::istarts_with( teststr, ss );
        
      case TextEndsWith:
        return UtilityFunctions::iends_with( teststr, ss );
        
      case TextRegex:
      {
        //Next line can throw, and I think is doing a lot of needless work inside
        //  the loop, so should maybe be optimized some day
        boost::regex expression( ss, boost::regex::icase );
        return boost::regex_match( teststr, expression );
      }//case TextRegex:
      break;
    }//switch( t )
    
    throw runtime_error( "SpecTest::test_string Needs updating" );
    
    return false;
  }//SpecTest::test_string
  
  
  bool SpecTest::test( const std::shared_ptr<const MeasurementInfo> meas ) const
  {
    if( !meas )
      return false;

    switch( m_searchField )
    {
      case Filename:
        return test_string( meas->filename(), m_stringSearchType, m_searchString );
        
      case DetectorName:
      {
        const vector<string> &detnames = meas->detector_names();
        bool found = false;
        for( size_t i = 0; !found && i < detnames.size(); ++i )
          found |= test_string( detnames[i], m_stringSearchType, m_searchString );
        return found;
      }
        
      case SerialNumber:
        return test_string( meas->instrument_id(), m_stringSearchType, m_searchString );
        
      case Manufacturer:
        return test_string( meas->manufacturer(), m_stringSearchType, m_searchString );
        
      case Model:
        return test_string( meas->instrument_model(), m_stringSearchType, m_searchString );
        
      case Uuid:
        return test_string( meas->uuid(), m_stringSearchType, m_searchString );
        
      case Remark:
      {
        const vector<string> &remarks = meas->remarks();
        bool found = false;
        for( size_t i = 0; !found && i < remarks.size(); ++i )
          found |= test_string( remarks[i], m_stringSearchType, m_searchString );
        
        if( found )
          return found;
        
        std::vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
        for( size_t i = 0; !found && i < meass.size(); ++i )
        {
          const vector<string> &measremarks = meas->remarks();
          bool found = false;
          for( size_t i = 0; !found && i < measremarks.size(); ++i )
            found |= test_string( measremarks[i], m_stringSearchType, m_searchString );
        }
        
        return found;
      }
        
      case LocationName:
        return test_string( meas->measurement_location_name(), m_stringSearchType, m_searchString );
        
      case AnalysisResultText:
      {
        const std::shared_ptr<const DetectorAnalysis> ana = meas->detectors_analysis();
        
        if( !ana )
          return false;
        
        const vector<string> &remarks = ana->remarks_;
        const vector<DetectorAnalysisResult> &res = ana->results_;
        
        if( test_string( ana->algorithm_version_, m_stringSearchType, m_searchString ) )
          return true;
        
        if( test_string( ana->algorithm_name_, m_stringSearchType, m_searchString ) )
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
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        const SandiaDecay::Nuclide * const nuc = db->nuclide( m_searchString );
        
        const std::shared_ptr<const DetectorAnalysis> ana = meas->detectors_analysis();
        
        if( !ana )
          return false;
        
        //const vector<string> &remarks = ana->remarks_;
        const vector<DetectorAnalysisResult> &res = ana->results_;
        
        if( nuc )
        {
          for( size_t i = 0; i < res.size(); ++i )
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
            
            for( size_t j = 0; j < fields.size(); ++j )
            {
              testnuc = db->nuclide( fields[j] );
              if( nuc == testnuc )
                return true;
            }
          }//for( size_t i = 0; i < res.size(); ++i )
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
        return (meas->detector_type() == type);
        
        break;
      }//case DetectionSystemType:
        
      case SearchMode:
        return ((m_discreteOption==0 && !meas->passthrough()) || (m_discreteOption && meas->passthrough()));
        
      case TotalLiveTime:
      {
        switch( m_compareType )
        {
          case ValueIsExact:       return (meas->gamma_live_time() == m_numeric);
          case ValueIsLessThan:    return (meas->gamma_live_time() < m_numeric);
          case ValueIsGreaterThan: return (meas->gamma_live_time() > m_numeric);
        }
        
        cerr << "SHouldnt be here m_compareType" << endl;
        return false;
      }
        
      case TotalRealTime:
      {
        switch( m_compareType )
        {
          case ValueIsExact:       return (meas->gamma_real_time() == m_numeric);
          case ValueIsLessThan:    return (meas->gamma_real_time() < m_numeric);
          case ValueIsGreaterThan: return (meas->gamma_real_time() > m_numeric);
        }
        
        cerr << "Shouldnt be here m_compareType" << endl;
        return false;
      }
        
      case IndividualSpectrumLiveTime:
      {
        const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
        for( size_t i = 0; i < meass.size(); ++i )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if(meass[i]->live_time() == m_numeric) return true;
            case ValueIsLessThan:    if(meass[i]->live_time() < m_numeric) return true;
            case ValueIsGreaterThan: if(meass[i]->live_time() > m_numeric) return true;
          }
        }
        return false;
      }//case IndividualSpectrumLiveTime:
        
      case IndividualSpectrumRealTime:
      {
        const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
        for( size_t i = 0; i < meass.size(); ++i )
        {
          switch( m_compareType )
          {
            case ValueIsExact:       if(meass[i]->real_time() == m_numeric) return true;
            case ValueIsLessThan:    if(meass[i]->real_time() < m_numeric) return true;
            case ValueIsGreaterThan: if(meass[i]->real_time() > m_numeric) return true;
          }
        }
        return false;
      }//case IndividualSpectrumRealTime:
        

      case StartTime:
      {
        const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
        for( size_t i = 0; i < meass.size(); ++i )
        {
          if( meass[i]->start_time().is_special() )
            continue;
          
          switch( m_compareType )
          {
            case ValueIsExact:
            {
              //We only have minute resoltion on the gui selector
              boost::posix_time::time_duration dur = meass[i]->start_time() - m_time;
              if( dur.is_negative() )
                dur = -dur;
                
              if( dur < boost::posix_time::minutes(1) )
                return true;
            }
              
            case ValueIsLessThan:    if(meass[i]->start_time() < m_time) return true;
            case ValueIsGreaterThan: if(meass[i]->start_time() > m_time) return true;
          }
        }
        return false;
      }//case StartTime:
        
      case NumFileDataFields:
        return false;
        break;
    }//switch( m_searchField )
    
    return false;
  }
  
  const char *to_string( const FileDataField field )
  {
    switch( field )
    {
      case Filename:              return "Filename";
      case DetectorName:          return "Detector name";
      case SerialNumber:          return "Serial number";
      case Manufacturer:          return "Manufacturer";
      case Model:                 return "Model";
      case Uuid:                  return "UUID";
      case Remark:                return "Remark";
      case LocationName:          return "Location name";
      case AnalysisResultText:    return "Analysis result";
      case AnalysisResultNuclide: return "Identified nuclide";
      case DetectionSystemType:   return "Detector system";
      case SearchMode:            return "Aquisition mode";
      case TotalLiveTime:         return "Sum live time";
      case TotalRealTime:         return "Sum real time";
      case IndividualSpectrumLiveTime: return "Spectrum live time";
      case IndividualSpectrumRealTime: return "Spectrum real time";
      case StartTime:             return "Start Time";
        
      case NumFileDataFields:
        break;
    }//switch( m_searchField )
    
    return "";
  }//to_string( const FileDataField field )
  
  const char *to_string( const TextFieldSearchType type )
  {
    switch( type )
    {
      case TextIsExact:     return "exactly matches";
      case TextIsContained: return "contains";
      case TextStartsWith:  return "begins with";
      case TextEndsWith:    return "ends with";
      case TextRegex:       return "matches regular expression";
        break;
    }
    
    return "";
  }//const char *to_string( const TextFieldSearchType type );
  
  
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
    }
	return "";
  }
  
  std::string SpecTest::summary() const
  {
    string summary = to_string( m_searchField );
    
    switch( m_searchField )
    {
      //String type searches
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
        summary += string(" ") + to_string( m_stringSearchType )
                   + string(" \"") + m_searchString + "\"";
        break;
        
      case DetectionSystemType:
        summary += " is a " + detectorTypeToString( DetectorType(m_discreteOption) );
        break;
      
      case SearchMode:
        summary += " is " + string(m_discreteOption ? "search/passthrough" : "dwell");
        break;
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
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
      case Filename: case DetectorName: case SerialNumber: case Manufacturer:
      case Model: case Uuid: case Remark: case LocationName:
      case AnalysisResultText: case AnalysisResultNuclide:
      {
        switch ( m_stringSearchType)
        {
          case TextIsExact:
          case TextIsContained:
            break;
            
          case TextStartsWith:
          case TextEndsWith:
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
          case kDetectiveEx100Detector: case kOrtecIDMPortalDetector:
          case kSAIC8Detector: case kFalcon5000: case kUnknownDetector:
          case kMicroDetectiveDetector: case kMicroRaiderDetector:
          case kRadHunterNaI: case kRadHunterLaBr3: case kRsi701: case kRsi705:
          case kAvidRsi:
          case kOrtecRadEagleNai: case kOrtecRadEagleCeBr2Inch:
          case kOrtecRadEagleCeBr3Inch: case kOrtecRadEagleLaBr:
          case kSam940LaBr3: case kSam940: case kSam945:
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
        
      case TotalLiveTime: case TotalRealTime:
      case IndividualSpectrumLiveTime: case IndividualSpectrumRealTime:
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
  
  
  SpecFileQuery::SpecFileQuery()
  {
  }
  
  void SpecFileQuery::addCondition( const SpecTest &test )
  {
    m_fields.push_back( boost::any(test) );
  }
  
  void SpecFileQuery::addLogic( const LogicType test )
  {
    m_fields.push_back( boost::any(test) );
  }
  
  
  std::string SpecFileQuery::summary() const
  {
    stringstream strm;
    print_equation( m_fields, strm );
    return strm.str();
  }
  
  std::ostream &SpecFileQuery::print_equation( vector<boost::any> fields, std::ostream &strm )
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
    }//for( size_t i = 0; i < fields.size(); ++i )
    
    return strm;
  }//std::ostream &print_equation(...)
  
  
  bool SpecFileQuery::evaluate( vector<boost::any> fields,
                                const std::shared_ptr<const MeasurementInfo> &meas )
  {
    if( !meas )
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
          
          const bool answer = SpecFileQuery::evaluate( inside, meas );
          
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
      }catch(... ){}
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
  }//bool evaluate(...)
  
  
  bool SpecFileQuery::test( const std::shared_ptr<const MeasurementInfo> meas ) const
  {
    return evaluate( m_fields, meas );
  }
  
  
  void SpecFileQuery::isvalid()
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
        bool istest = false;
        SpecTest test;
        LogicType closeparen;
        try
        {
          test = boost::any_cast<SpecTest>( m_fields[i] );
          istest = true;
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
                else
                  throw runtime_error( "Closing parenthesis should not be followed " + string(to_string(nextt)) );
                break;
            }//switch( nextt )
          }catch(...)
          {
            if( istest )
              throw runtime_error( test.summary() + " should be followed by an operation" );
            else
               throw runtime_error( "Closing parenthesis should be followed by an operation" );
          }
        }//if( not the last condition )
      }//for( size_t i = 0; i < m_fields.size(); ++i )
    }//END 3) except for last condition, all conditions must be followed by || && or )
    
    {//BEGIN 4) All regular exression matches should have valid regex
      for( size_t i = 0; i < m_fields.size(); ++i )
      {
        SpecTest test;
        try
        {
          test = boost::any_cast<SpecTest>( m_fields[i] );
        }catch(...)
        {
          continue;
        }
        
        test.isvalid();
      }//for( size_t i = 0; i < m_fields.size(); ++i )
    }//END 4) All regular exression matches should have valid regex
    
    {//BEGIN 5) And && followed by a condition then an || is ambiguouse, and visa versa
      cout << "Need to enforce no ambigous order of evaluation of ORs and ANDs" << endl;
    }//END 5) And && followed by a condition then an || is ambiguouse, and visa versa
    
  }//void isvalid()
}//namespace SpecQuery


namespace
{
  using namespace SpecQuery;
  
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
  
  
  
  //GUI Widgets used
  class ConditionWidget : public WContainerWidget
  {
  public:
    ConditionWidget( WContainerWidget *parent = 0 )
      : WContainerWidget( parent )
    {
      addStyleClass( "SpecFileQueryCond" );
      
      WGridLayout *layout = new WGridLayout();
      layout->setContentsMargins( 3, 2, 3, 2 );
      layout->setVerticalSpacing( 0 );
      layout->setHorizontalSpacing( 0 );
      setLayout( layout );
      
      m_typeSelect = new WComboBox();
      layout->addWidget( m_typeSelect, 0, 0, AlignMiddle );
      
      
      for( SpecQuery::FileDataField t = SpecQuery::FileDataField(0);
          t < SpecQuery::NumFileDataFields; t = SpecQuery::FileDataField(t+1) )
      {
        switch( t )
        {
          case SpecQuery::Filename:        m_typeSelect->addItem( "Filename" );        break;
          case DetectorName:               m_typeSelect->addItem( "Detectors Name" );  break;
          case SerialNumber:               m_typeSelect->addItem( "Serial Number" );   break;
          case Manufacturer:               m_typeSelect->addItem( "Manufacturer" );    break;
          case Model:                      m_typeSelect->addItem( "Model" );           break;
          case Uuid:                       m_typeSelect->addItem( "UUID" );            break;
          case Remark:                     m_typeSelect->addItem( "Remarks" );         break;
          case LocationName:               m_typeSelect->addItem( "Location Name" );   break;
          case AnalysisResultText:         m_typeSelect->addItem( "Analysis Txt" );    break;
          case AnalysisResultNuclide:      m_typeSelect->addItem( "Analysis Nuclide" );break;
          case DetectionSystemType:        m_typeSelect->addItem( "Detector Type" );   break;
          case SearchMode:                 m_typeSelect->addItem( "Aquisition Mode" ); break;
          case TotalLiveTime:              m_typeSelect->addItem( "Sum Live Time" ); break;
          case TotalRealTime:              m_typeSelect->addItem( "Sum Real Time" ); break;
          case IndividualSpectrumLiveTime: m_typeSelect->addItem( "Spec Live Time" ); break;
          case IndividualSpectrumRealTime: m_typeSelect->addItem( "Spec Real Time" ); break;
          case StartTime:                  m_typeSelect->addItem( "Start Time" ); break;
        
          case SpecQuery::NumFileDataFields: break;
        }
      }//for( loop over FileDataField )
      
      
      for( SpecQuery::LogicType t = SpecQuery::LogicalOr;
          t < SpecQuery::NumLogicType; t = SpecQuery::LogicType(t+1) )
      {
        switch( t )
        {
          case SpecQuery::LogicalOr:         m_typeSelect->addItem( "Logical OR" );        break;
          case SpecQuery::LogicalAnd:        m_typeSelect->addItem( "Logical AND" );       break;
          case SpecQuery::LogicalNot:        m_typeSelect->addItem( "Logical NOT" );       break;
          case SpecQuery::LogicalOpenParan:  m_typeSelect->addItem( "Open Parenthesis" );  break;
          case SpecQuery::LogicalCloseParan: m_typeSelect->addItem( "Close Parenthesis" ); break;
          case SpecQuery::NumLogicType: break;
        }
      }//for( loop over LogicType )
        
      
      m_txtContent = new WContainerWidget();
      m_txtContent->setMargin( 2 );
      layout->addWidget( m_txtContent, 0, 1 );
      layout->setColumnStretch( 1, 1 );
      
      
      WGridLayout *txtLayout = new WGridLayout();
      txtLayout->setVerticalSpacing( 0 );
      txtLayout->setContentsMargins( 0, 0, 0, 0 );
      m_txtContent->setLayout( txtLayout );
      
      m_txtSearchType = new WComboBox();
      m_txtSearchType->setMargin( 2, Wt::Top );
      txtLayout->addWidget( m_txtSearchType, 0, 0 );
      
      
      for( SpecQuery::TextFieldSearchType t = SpecQuery::TextIsExact;
           t <= SpecQuery::TextRegex; t = SpecQuery::TextFieldSearchType(t+1) )
      {
        switch( t )
        {
          case TextIsExact:     m_txtSearchType->addItem( "Matches" );       break;
          case TextIsContained: m_txtSearchType->addItem( "Contains" );      break;
          case TextStartsWith:  m_txtSearchType->addItem( "Starts with" );   break;
          case TextEndsWith:    m_txtSearchType->addItem( "Ends with" );     break;
          case TextRegex:       m_txtSearchType->addItem( "Regex Matches" ); break;
        }
      }//for( loop over TextFieldSearchType )
      
      m_txtSearchType->setCurrentIndex( TextIsContained );
      m_txtSearchType->activated().connect( this, &ConditionWidget::searchTextTypeChanged );
      
      m_searchTxt = new WLineEdit();
      txtLayout->addWidget( m_searchTxt, 0, 1 );
      m_searchTxt->changed().connect( this, &ConditionWidget::searchTextChanged );
      m_searchTxt->enterPressed().connect( this, &ConditionWidget::searchTextChanged );
      txtLayout->setColumnStretch( 1, 1 );
      
      m_discreteContent = new WContainerWidget();
      m_discreteContent->setMargin( 2 );
      layout->addWidget( m_discreteContent, 0, 2 );
      layout->setColumnStretch( 2, 1 );
      m_discreteContent->hide();
      
      WGridLayout *discreteLayout = new WGridLayout();
      discreteLayout->setVerticalSpacing( 0 );
      discreteLayout->setContentsMargins( 0, 0, 0, 0 );
      m_discreteContent->setLayout( discreteLayout );
      
      m_discreteTxt = new WText();
      m_discreteTxt->setMargin( 15, Wt::Left );
      m_discreteTxt->setMargin( 2, Wt::Top );
      discreteLayout->addWidget( m_discreteTxt, 0, 0 );
      
      m_discreteCombo = new WComboBox();
      m_discreteCombo->setMargin( 5, Wt::Left );
      m_discreteCombo->setMargin( 2, Wt::Top );
      discreteLayout->addWidget( m_discreteCombo, 0, 1 );
      m_discreteCombo->activated().connect( this, &ConditionWidget::searchTextTypeChanged );
      
      discreteLayout->addWidget( new WContainerWidget(), 0, 2 );
      discreteLayout->setColumnStretch( 2, 1 );
      
      
      
      m_numericContent = new WContainerWidget();
      layout->addWidget( m_numericContent, 0, 3 );
      layout->setColumnStretch( 3, 1 );
      
      WGridLayout *numericLayout = new WGridLayout();
      numericLayout->setVerticalSpacing( 0 );
      numericLayout->setContentsMargins( 0, 0, 0, 0 );
      m_numericContent->setLayout( numericLayout );
      
      m_numericCombo = new WComboBox();
      m_numericCombo->setMargin( 5, Wt::Left );
      m_numericCombo->setMargin( 2, Wt::Top );
      numericLayout->addWidget( m_numericCombo, 0, 1 );
      m_numericCombo->addItem( "Matches" ); //ValueIsExact
      m_numericCombo->addItem( "Less than" ); //ValueIsLessThan,
      m_numericCombo->addItem( "Greater than" ); //ValueIsGreaterThan,
      m_numericCombo->setCurrentIndex( 1 );
      m_numericCombo->activated().connect( this, &ConditionWidget::searchTextTypeChanged );
      
      m_numericInput = new WLineEdit();
      m_numericInput->setMargin( 5, Wt::Left );
      numericLayout->addWidget( m_numericInput, 0, 2 );
      m_numericInput->changed().connect( this, &ConditionWidget::searchTextChanged );
      m_numericInput->enterPressed().connect( this, &ConditionWidget::searchTextChanged );
      numericLayout->setColumnStretch( 2, 1 );
      
      m_timeContent = new WContainerWidget();
      layout->addWidget( m_timeContent, 0, layout->columnCount() );
      layout->setColumnStretch( layout->columnCount()-1, 1 );
      m_timeContent->hide();
      WGridLayout *timeLayout = new WGridLayout();
      timeLayout->setVerticalSpacing( 0 );
      timeLayout->setContentsMargins( 0, 0, 0, 0 );
      m_timeContent->setLayout( timeLayout );

      m_timeCombo = new WComboBox();
      m_timeCombo->addItem( "Matches" ); //ValueIsExact
      m_timeCombo->addItem( "Less than" ); //ValueIsLessThan,
      m_timeCombo->addItem( "Greater than" ); //ValueIsGreaterThan,
      m_timeCombo->setCurrentIndex( 1 );
      m_timeCombo->activated().connect( this, &ConditionWidget::searchTextTypeChanged );
      timeLayout->addWidget( m_timeCombo, 0, 0 );
      
      m_dateInput = new WDatePicker();
      timeLayout->addWidget( m_dateInput, 0, 1 );
      m_dateInput->changed().connect( this, &ConditionWidget::searchTextChanged );
      
#if( WT_VERSION >= 0x3030800 )
      WTimeEdit *timeedit = new WTimeEdit();
      m_timeInput = new WTimePicker( timeedit );
      timeLayout->addWidget( m_timeInput, 0, 2 );
#else
      m_timeInput = new WTimePicker();
      timeLayout->addWidget( m_timeInput, 0, 2 );
#endif
      
      m_timeInput->selectionChanged().connect( this, &ConditionWidget::searchTextChanged );
      
      timeLayout->addWidget( new WContainerWidget(), 0, 3 );
      timeLayout->setColumnStretch( 3, 1 );

      
      
      //Add in some icons to add after, or delete this row
      WText *addIcon = new WText();
      addIcon->setStyleClass( "AddIcon" );
      layout->addWidget( addIcon, 0, layout->columnCount() );
      addIcon->clicked().connect( this, &ConditionWidget::addRequested );
      
      WText *closeIcon = new WText();
      closeIcon->setStyleClass( "DeleteIcon" );
      closeIcon->clicked().connect( this, &ConditionWidget::removeRequested );
      layout->addWidget( closeIcon, 0, layout->columnCount() );
      
      
      m_typeSelect->setCurrentIndex( SerialNumber );
      m_typeSelect->activated().connect( this, &ConditionWidget::typeChanged );
      
      typeChanged( SerialNumber );
    }//ConditionWidget constructor
    
    
    bool isCondition()
    {
      return (m_typeSelect->currentIndex() < NumFileDataFields);
    }
    
    bool isCloseParan()
    {
      return (m_typeSelect->currentIndex() == (int(NumFileDataFields) + int(LogicalCloseParan)));
    }
    
    void setIsAnd()
    {
      const int index = int(NumFileDataFields) + LogicalAnd;
      m_typeSelect->setCurrentIndex( index );
      typeChanged( index );
    }
    
    
    SpecQuery::SpecTest test()
    {
      const int index = m_typeSelect->currentIndex();
      if( index < 0 || (index >= NumFileDataFields) )
        throw runtime_error( "ConditionWidget can not provide SpecTest for a logical widget" );
      
      const FileDataField type = FileDataField( index );
      
      SpecQuery::SpecTest t;
    
      switch( type )
      {
        case Filename: case DetectorName: case SerialNumber: case Manufacturer:
        case Model: case Uuid: case Remark: case LocationName:
        case AnalysisResultText: case AnalysisResultNuclide:
        {
          const string searctxt = m_searchTxt->text().toUTF8();
          const TextFieldSearchType searchtype = TextFieldSearchType( m_txtSearchType->currentIndex() );
          t.set_test( type, searctxt, searchtype );
          break;
        }
          
        case DetectionSystemType:
        case SearchMode:
          t.set_discreet_test( type, m_discreteCombo->currentIndex() );
          break;
        
        case TotalLiveTime:
        case TotalRealTime:
        case IndividualSpectrumLiveTime:
        case IndividualSpectrumRealTime:
        {
          double value = 0.0;
          try
          {
            const string val = m_numericInput->text().toUTF8();
            value = PhysicalUnits::stringToTimeDuration( val );
          }catch(...)
          {
            value = std::numeric_limits<double>::infinity();
          }
          
          //const bool parsed = UtilityFunctions::parse_float( val.c_str(), val.size(), value );
          t.set_numeric_test( type, value, NumericFieldMatchType(m_numericCombo->currentIndex()) );
          break;
        }
          
        case StartTime:
        {
          const WDate wanteddate = m_dateInput->date();
          const WTime wantedtime = m_timeInput->time();
          const boost::gregorian::date date = wanteddate.toGregorianDate();
          const boost::posix_time::time_duration duration = wantedtime.toTimeDuration();
          
          const boost::posix_time::ptime time( date, duration );
          
          t.set_time_test( type, time, NumericFieldMatchType(m_timeCombo->currentIndex()) );
        }
          
        case NumFileDataFields:
          break;
      }//switch( type )
      
      return t;
    }//SpecQuery::SpecTest test()
    
    LogicType logic()
    {
      const int index = m_typeSelect->currentIndex();
      if( (index < NumFileDataFields) || ((index-NumFileDataFields) >= NumLogicType) )
        throw runtime_error( "ConditionWidget can not provide LogicType for a test widget" );
      return LogicType( index - NumFileDataFields );
    }//logic()
    
    
    Wt::Signal<> &changed(){ return m_changed; }
    Wt::Signal<WWebWidget *> &add(){ return m_add; }
    Wt::Signal<WWebWidget *> &remove(){ return m_remove; }
    
  protected:
    
    void typeChanged( int index )
    {
      setStyleClass( "SpecFileQueryCond" );
      
      if( index < NumFileDataFields )
      {
        m_searchTxt->enable();
        m_searchTxt->setText( "" );
        m_searchTxt->setPlaceholderText( "" );
        
        m_txtSearchType->show();
        m_searchTxt->setAttributeValue( "style", "text-align: left;" );
        
        m_txtContent->hide();
        m_discreteContent->hide();
        m_numericContent->hide();
        m_timeContent->hide();
        
        const SpecQuery::FileDataField field = SpecQuery::FileDataField(index);
        
        const char *styleclass = "FileQueryTextSearch";
        
        switch( field )
        {
          case Filename:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter filesystem path" );
            break;
            
          case DetectorName:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter a detector name in file (ex A1, BA2, ...)" );
            break;
            
          case SerialNumber:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter serial number given in file" );
            break;
            
          case Manufacturer:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter manufacturer wanted (searches from file or infered)" );
            break;
            
          case Model:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter model wanted (searches from file or infered)" );
            break;
            
          case Uuid:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter UUID search (in file or else computed)" );
            break;
            
          case Remark:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter remarks in file" );
            break;
            
          case LocationName:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter location text given in file" );
            break;
            
          case AnalysisResultText:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter text contained in analysis" );
            break;
            
          case AnalysisResultNuclide:
            m_txtContent->show();
            m_searchTxt->setPlaceholderText( "Enter nuclides given by analysis in file (smart, ex \"Co60\" matches \"co 60\" and \"CO-60\", etc.)" );
            break;
            
          case DetectionSystemType:
            m_discreteContent->show();
            m_discreteTxt->setText( "determined to be" );
            m_discreteCombo->clear();
            
            for( DetectorType type = DetectorType(0); type <= kSam945; type = DetectorType(type+1) )
            {
              switch( type )
              {
                case kGR135Detector:             m_discreteCombo->addItem( "GR135" ); break;
                case kIdentiFinderDetector:      m_discreteCombo->addItem( "identiFINDER (first gen)" ); break;
                case kIdentiFinderNGDetector:    m_discreteCombo->addItem( "identiFINDER NG/H" ); break;
                case kIdentiFinderLaBr3Detector: m_discreteCombo->addItem( "identiFINDER LaBr3" ); break;
                case kDetectiveDetector:         m_discreteCombo->addItem( "Detective (unknown model)" ); break;
                case kDetectiveExDetector:       m_discreteCombo->addItem( "Detective EX/DX" ); break;
                case kDetectiveEx100Detector:    m_discreteCombo->addItem( "Detective EX/DX-100" ); break;
                case kOrtecIDMPortalDetector:    m_discreteCombo->addItem( "Ortec IDM" ); break;
                case kSAIC8Detector:             m_discreteCombo->addItem( "SAIC8 Portal" ); break;
                case kFalcon5000:                m_discreteCombo->addItem( "Falcon 5000" ); break;
                case kUnknownDetector:           m_discreteCombo->addItem( "Unknown Detector" ); break;
                case kMicroDetectiveDetector:    m_discreteCombo->addItem( "Micro Detective" ); break;
                case kMicroRaiderDetector:       m_discreteCombo->addItem( "Micro Raider" ); break;
                case kRadHunterNaI:              m_discreteCombo->addItem( "Rad Hunter NaI" ); break;
                case kRadHunterLaBr3:            m_discreteCombo->addItem( "Rad Hunter LaBr3" ); break;
                case kRsi701:                    m_discreteCombo->addItem( "RSI 701" ); break;
                case kRsi705:                    m_discreteCombo->addItem( "RSI 705" ); break;
                case kAvidRsi:                   m_discreteCombo->addItem( "Avid RSI" ); break;
                case kOrtecRadEagleNai:          m_discreteCombo->addItem( "Ortec RadEagle NaI(Tl) 3x1" ); break;
                case kOrtecRadEagleCeBr2Inch:    m_discreteCombo->addItem( "Ortec RadEagle CeBr3 2x1" ); break;
                case kOrtecRadEagleCeBr3Inch:    m_discreteCombo->addItem( "Ortec RadEagle CeBr3 3x0.8" ); break;
                case kOrtecRadEagleLaBr:         m_discreteCombo->addItem( "Ortec RadEagle LaBr3(Ce) 2x1" ); break;
                case kSam940LaBr3:               m_discreteCombo->addItem( "Sam 940 LaBr3" ); break;
                case kSam940:                    m_discreteCombo->addItem( "Sam 940" ); break;
                case kSam945:                    m_discreteCombo->addItem( "Sam 945" ); break;
              }//switch( type )
              
              m_discreteCombo->setCurrentIndex( kUnknownDetector );
            }//for( loop over DetectorType )
            break;
            
          case SearchMode:
            m_discreteContent->show();
            m_discreteTxt->setText( "Detector type is" );
            m_discreteCombo->clear();
            
            m_discreteCombo->addItem( "Dwell" );
            m_discreteCombo->addItem( "Search/Passthrough" );
            m_discreteCombo->setCurrentIndex( 0 );
            break;
            
          case TotalLiveTime:
            m_numericContent->show();
            m_numericInput->setEmptyText( "Enter total live time of file (ex \"5s\" or \"3days 2h 15s\")" );
            break;
            
          case TotalRealTime:
            m_numericContent->show();
            m_numericInput->setEmptyText( "Enter total wall time of file (ex \"5s\" or \"3days 2h 15s\")" );
            break;
            
          case IndividualSpectrumLiveTime:
            m_numericContent->show();
            m_numericInput->setEmptyText( "Enter live time for any spectrum in file (ex \"5.4 seconds\")" );
            break;
            
          case IndividualSpectrumRealTime:
            m_numericContent->show();
            m_numericInput->setEmptyText( "Enter wall time for any spectrum in file (ex \"5.4 seconds\")" );
            break;
            
          case StartTime:
            m_timeContent->show();
            break;
            
          case NumFileDataFields:
            throw runtime_error( "Logic error for SpecQuery::LogicType 0" );
            break;
        }//switch( field )
        
        addStyleClass( styleclass );
      }else
      {
        const SpecQuery::LogicType type = SpecQuery::LogicType(index - NumFileDataFields);
        if( type >= NumLogicType )
          throw runtime_error( "Logic error for SpecQuery::LogicType" );
        
        m_txtContent->show();
        m_discreteContent->hide();
        
        m_searchTxt->disable();
        m_searchTxt->setPlaceholderText( "" );
        
        m_txtSearchType->hide();
        m_searchTxt->setAttributeValue( "style", "text-align: center;" );
        
        switch( type )
        {
          case LogicalOr:
            addStyleClass( "FileQueryOr" );
            m_searchTxt->setText( "OR" );
            break;
            
          case LogicalAnd:
            addStyleClass( "FileQueryAnd" );
            m_searchTxt->setText( "AND" );
            break;
            
          case LogicalNot:
            addStyleClass( "FileQueryNot" );
            m_searchTxt->setText( "NOT" );
            break;
            
          case LogicalOpenParan:
            addStyleClass( "FileQueryOpenParan" );
            m_searchTxt->setText( "Open Parenthesis" );
            break;
            
          case LogicalCloseParan:
            addStyleClass( "FileQueryCloseParan" );
            m_searchTxt->setText( "Close Parenthesis" );
            break;
            
          case NumLogicType:
            assert( 0 );
            break;
        }//switch( type )
        
      }//if( index < NumFileDataFields ) / else
      
      m_changed.emit();
    }//void typeChanged( int index )
    
    void searchTextTypeChanged()
    {
      m_changed.emit();
    }
    
    void searchTextChanged()
    {
      m_changed.emit();
    }
    
    void removeRequested()
    {
      m_remove.emit( this );
    }
    
    void addRequested()
    {
      m_add.emit( this );
    }
    
  protected:
    
    WComboBox *m_typeSelect;
    
    WContainerWidget *m_txtContent;
    WLineEdit *m_searchTxt;
    WComboBox *m_txtSearchType;
    
    WContainerWidget *m_discreteContent;
    WText *m_discreteTxt;
    WComboBox *m_discreteCombo;
    
    WContainerWidget *m_numericContent;
    WComboBox *m_numericCombo;
    WLineEdit *m_numericInput;
    
    WContainerWidget *m_timeContent;
    WComboBox *m_timeCombo;
    WDatePicker *m_dateInput;
    WTimePicker *m_timeInput;
    
    
    
    Wt::Signal<> m_changed;
    Wt::Signal<WWebWidget *> m_add, m_remove;
  };//class ConditionWidget
  
  
  void output_csv_field( std::ostream &out, std::string s)
  {
    //adapted from http://darrendev.blogspot.com/2009/11/escaping-csv-in-c.html 20160623
    //There are two escaping rules for each field in a comma-separated value row:
    //  1. Change each double quote to two double quotes.
    //  2. Surround with double quotes if the field contains a comma or double quote.
    
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
    ResultCsvResource( WAbstractItemModel *model, WAnchor *parent )
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
          const boost::any txt = m_model->data( row, col, (col==0 ? UserRole : DisplayRole) );
          const WString txtstr = Wt::asString( txt );
          response.out() << (col ? "," : "");
          output_csv_field(	response.out(), txtstr.toUTF8() );
        }
        
        response.out() << "\r\n";
      }//for( int row = 0; row < nrow; ++row )
    }//handleRequest

  };//class PeakCsvResource : public Wt::WResource

  
  
  vector<string> get_result_fields( std::shared_ptr<MeasurementInfo> &meas )
  {
    vector<string> row( NumFileDataFields, "" );
    std::shared_ptr<const DetectorAnalysis> anares = meas->detectors_analysis();
    
    for( FileDataField f = FileDataField(0); f < NumFileDataFields; f = FileDataField(f+1) )
    {
      switch( f )
      {
        case Filename:
          row[f] = meas->filename();
          break;
          
        case DetectorName:
          for( size_t i = 0; i < meas->detector_names().size(); ++i )
            row[f] += (i ? ";" : "") + meas->detector_names()[i];
          break;
          
        case SerialNumber:
          row[f] = meas->instrument_id();
          break;
          
        case Manufacturer:
          row[f] = meas->manufacturer();
          break;
          
        case Model:
          row[f] = meas->instrument_model();
          break;
          
        case Uuid:
          row[f] = meas->uuid();
          break;
          
        case Remark:
          for( size_t i = 0; i < meas->remarks().size() && row[f].size() < 1024; ++i )
            row[f] += (row[f].size() ? "\n" : "") + meas->remarks()[i];
          
          for( size_t i = 0; i < meas->num_measurements() && row[f].size() < 1024; ++i )
          {
            std::shared_ptr<const Measurement>  m = meas->measurement(i);
            for( size_t j = 0; j < m->remarks().size() && row[f].size() < 1024; ++j )
              row[f] += (row[f].size() ? "\n" : "") + m->remarks()[j];
          }
          
          if( row[f].size() >= 1024 )
            row[f] += "\n...";
          break;
          
        case LocationName:
          row[f] = meas->measurement_location_name();
          break;
          
        case AnalysisResultText:
          if( anares )
          {
            for( size_t i = 0; i < anares->remarks_.size() && row[f].size() < 1024; ++i )
              row[f] += (row[f].size() ? "\n" : "") + anares->remarks_[i];
            
            for( size_t i = 0; i < anares->results_.size() && row[f].size() < 1024; ++i )
            {
              string r;
              if( anares->results_[i].nuclide_type_.size() )
                r += anares->results_[i].nuclide_type_ + ";";
              if( anares->results_[i].id_confidence_.size() )
                r += (r.size() ? ";" : "") + anares->results_[i].id_confidence_ + ";";
              if( anares->results_[i].remark_.size() )
                r += (r.size() ? ";" : "") + anares->results_[i].remark_ + ";";
              if( r.size() )
                row[f] += (row[f].size() ? "\n" : "") + r;
            }
            
            if( row[f].size() >= 1024 )
              row[f] += "\n...";
          }//if( anares )
          break;
          
          
        case AnalysisResultNuclide:
          if( anares )
          {
            for( size_t i = 0; i < anares->results_.size(); ++i )
              row[f] += (row[f].size() ? ";" : "") + anares->results_[i].nuclide_;
          }
          break;
          
        case DetectionSystemType:
          row[f] = detectorTypeToString( meas->detector_type() );
          break;
          
        case SearchMode:
          row[f] = (meas->passthrough() ? "Search/Passthrough" : "Dwell");
          break;
          
        case TotalLiveTime:
          row[f] = PhysicalUnits::printToBestTimeUnits( meas->gamma_live_time() );
          break;
          
        case TotalRealTime:
          row[f] = PhysicalUnits::printToBestTimeUnits( meas->gamma_live_time() );
          break;
          
        case IndividualSpectrumLiveTime:
        {
          const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
          for( size_t i = 0; i < meass.size(); ++i )
          {
            row[f] += (i?";":"") + PhysicalUnits::printToBestTimeUnits( meass[i]->live_time() );
            if( row[f].size() > 1024 )
              break;
          }
          break;
        }
          
        case IndividualSpectrumRealTime:
        {
          const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
          for( size_t i = 0; i < meass.size(); ++i )
          {
            row[f] += (i?";":"") + PhysicalUnits::printToBestTimeUnits( meass[i]->real_time() );
            if( row[f].size() > 1024 )
              break;
          }
          break;
        }
        
        case StartTime:
        {
          const vector< std::shared_ptr<const Measurement> > meass = meas->measurements();
          for( size_t i = 0; i < meass.size(); ++i )
          {
            row[f] += (i?";":"") + UtilityFunctions::to_iso_string( meass[i]->start_time() );
            if( row[f].size() > 1024 )
              break;
          }
          
          break;
        }
          
        case NumFileDataFields:
          break;
      }//switch( f )
    }//for( FileDataField f = FileDataField(0); f < NumFileDataFields; f = FileDataField(f+1) )
    
    return row;
  }//vector<string> get_result_fields( std::shared_ptr<const MeasurementInfo> &meas )
  

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
    a minute just to insert new results for the Traige Data file set that had
    18.7k rows.  This also allows us to be a little more memorry efficient as 
    well.
 */
class ResultTableModel : public WAbstractItemModel
{
protected:
  std::shared_ptr< vector< vector<string> > > m_result;
  int m_sortcolumn;
  Wt::SortOrder m_sortorder;
  boost::any m_summary;
  
  
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
      
      switch( FileDataField(m_column) )
      {
        case Filename:
          lesthan = (UtilityFunctions::filename(lhs[m_column]) < UtilityFunctions::filename(rhs[m_column]));
          break;
          
        case DetectorName: case SerialNumber: case Manufacturer:
        case Model: case Uuid: case Remark: case LocationName:
        case AnalysisResultText: case AnalysisResultNuclide:
        case DetectionSystemType:
        case SearchMode:
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
    return NumFileDataFields;
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
    if( row < 0 || col < 0 || col >= NumFileDataFields
       || row >= static_cast<int>(m_result->size()) )
      return boost::any();
    
    const vector<string> &fields = (*m_result)[row];
    
    if( col == 0 )
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
      return boost::any();
    
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
    if( section < 0 || section >= NumFileDataFields || orientation != Horizontal )
      return boost::any();
    if( section == 0 && role == Wt::UserRole )
      return m_summary;
    
    if( role != DisplayRole )
      return boost::any();
    
    switch( SpecQuery::FileDataField(section) )
    {
      case Filename:     return WString("Filename");
      case DetectorName:          return WString("Detectors Name");
      case SerialNumber:          return WString("Serial Number");
      case Manufacturer:          return WString("Manufacturer");
      case Model:                 return WString("Model");
      case Uuid:                  return WString("UUID");
      case Remark:                return WString("Remarks");
      case LocationName:          return WString("Location Name");
      case AnalysisResultText:    return WString("Analysis Txt");
      case AnalysisResultNuclide: return WString("Analysis Nuclide");
      case DetectionSystemType:   return WString("Detector Type");
      case SearchMode:            return WString("Aquisition Mode");
      case TotalLiveTime:         return WString("Sum Live Time");
      case TotalRealTime:         return WString("Sum Wall Time");
      case IndividualSpectrumLiveTime: return WString("Spec Live Times");
      case IndividualSpectrumRealTime: return WString("Spec Wall Times");
      case StartTime:             return WString("Start Time");
        
      case SpecQuery::NumFileDataFields: break;
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
    m_conditions( 0 ),
    m_addCondition( 0 ),
    m_update( 0 ),
    m_cancelUpdate( 0 ),
    m_loadSelectedFile( 0 ),
    m_resultmodel( 0 ),
    m_resultview( 0 ),
    m_baseLocation( 0 ),
    m_recursive( 0 ),
    m_filterByExtension( 0 ),
    m_maxFileSize( 0 ),
    m_filterUnique( 0 ),
    m_numberFiles( 0 ),
    m_numberResults( 0 ),
    m_csv( 0 )
{
  init();
}//SpecFileQueryWidget constructor


SpecFileQueryWidget::~SpecFileQueryWidget()
{
  m_widgetDeleted->store( true );
}//~SpecFileQueryWidget()


void SpecFileQueryWidget::init()
{
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 9, 9, 9, 0 );
  setLayout( layout );
  
  WContainerWidget *line = new WContainerWidget();
  WGridLayout *linelayout = new WGridLayout( line );
  
  string defpath;
  int maxsize = 32;
  bool dofilter = true, dorecursive = true;
  if( m_viewer )
  {
    try
    {
      dofilter = InterSpecUser::preferenceValue<bool>( "SpecFileQueryFilter", m_viewer );
      dorecursive = InterSpecUser::preferenceValue<bool>( "SpecFileQueryRecursive", m_viewer );
      maxsize = InterSpecUser::preferenceValue<int>( "SpecFileQueryMaxSize", m_viewer );
      defpath = InterSpecUser::preferenceValue<string>( "SpecFileQueryPath", m_viewer );
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
  
  m_baseLocation = new WLineEdit();
  
  if( !defpath.empty() )
    m_baseLocation->setText( defpath );
  
  linelayout->addWidget( m_baseLocation, 0, 1 );
  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_baseLocation->enterPressed().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_baseLocation->setEmptyText( "Path to base directory to search" );
  
  m_recursive = new WCheckBox( "recursive" );
  m_recursive->setChecked( dorecursive );
  m_recursive->setMargin( 2, Wt::Top );
  linelayout->addWidget( m_recursive, 0, 2 );
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
  
  m_addCondition = new WText();
  m_addCondition->setStyleClass( "AddIcon" );
  m_addCondition->clicked().connect( this, &SpecFileQueryWidget::addCondition );
  linelayout->addWidget( m_addCondition, 0, linelayout->columnCount() );
  
//  linelayout->addWidget( new WContainerWidget(), 0, linelayout->columnCount() );
//  linelayout->setColumnStretch( linelayout->columnCount()-1, 1 );

  label = new WLabel( "Max file size:" );
  label->setMargin( 4, Wt::Bottom );
  label->setMargin( 15, Wt::Left );
  linelayout->addWidget( label, 0, linelayout->columnCount(), AlignBottom );
  
  m_maxFileSize = new WSpinBox();
  m_maxFileSize->setMaximum( 2*1024 );
  m_maxFileSize->setMinimum( 1 );
  m_maxFileSize->setValue( maxsize );
  m_maxFileSize->setWidth( 55 );
  linelayout->addWidget( m_maxFileSize, 0, linelayout->columnCount(), AlignCenter );
  label->setBuddy( m_maxFileSize );
  m_maxFileSize->setToolTip( "Maximum size of file to attempt to parse as a spectrum file." );
  m_maxFileSize->valueChanged().connect( this, &SpecFileQueryWidget::basePathChanged );
  
  
  label = new WLabel( "MB" );
  label->setMargin( 6, Wt::Right );
  label->setMargin( 4, Wt::Bottom );
  linelayout->addWidget( label, 0, linelayout->columnCount(), AlignBottom );
  label->setBuddy( m_maxFileSize );

  m_filterByExtension = new WCheckBox( "pre-filter by extension" );
  m_filterByExtension->setChecked( dofilter );
  m_filterByExtension->setMargin( 3, Wt::Bottom );
  linelayout->addWidget( m_filterByExtension, 0, linelayout->columnCount(), AlignBottom );
  m_filterByExtension->checked().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_filterByExtension->unChecked().connect( this, &SpecFileQueryWidget::basePathChanged );
  m_filterByExtension->setToolTip( "Eliminates common file types (ex. zip, doc, avi, etc) from search to speed results up." );
  m_filterByExtension->setMargin( 5, Wt::Right );
  
  
  m_filterUnique = new WCheckBox( "filter duplicate files" );
  m_filterUnique->setChecked( dofilter );
  m_filterUnique->setMargin( 3, Wt::Bottom );
  linelayout->addWidget( m_filterUnique, 0, linelayout->columnCount(), AlignBottom );
  m_filterUnique->checked().connect( this, &SpecFileQueryWidget::setResultsStale );
  m_filterUnique->unChecked().connect( this, &SpecFileQueryWidget::setResultsStale );
  m_filterUnique->setToolTip( "Filters duplicate spectrum files through use of a hash of spectral and meta-information." );
  m_filterUnique->setMargin( 5, Wt::Right );
  
  
  linelayout->addWidget( new WContainerWidget(), 0, linelayout->columnCount() );
  linelayout->setColumnStretch( linelayout->columnCount()-1, 6 );
  
  m_cancelUpdate = new WPushButton( "Cancel" );
  m_cancelUpdate->setHidden( true );
  linelayout->addWidget( m_cancelUpdate, 0, linelayout->columnCount() );
  m_cancelUpdate->clicked().connect( this, &SpecFileQueryWidget::cancelUpdate );
  
  m_update = new WPushButton( "Update" );
  m_update->clicked().connect( this, &SpecFileQueryWidget::startUpdate );
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
  
  //m_resultmodel = new WStandardItemModel( 0, NumFileDataFields, this );
  m_resultmodel = new ResultTableModel( this );
  
  m_resultview->setModel( m_resultmodel );
  m_resultview->setColumnHidden( DetectorName, true );
  m_resultview->setColumnHidden( Remark, true );
  m_resultview->setColumnHidden( LocationName, true );
  m_resultview->setColumnHidden( AnalysisResultText, true );
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
  linelayout->setColumnStretch( 1, 1 );
  
  if( m_viewer )
  {
    m_loadSelectedFile = new WPushButton( "Load Selected" );
    linelayout->addWidget( m_loadSelectedFile, 0, 2, AlignRight );
    m_loadSelectedFile->disable();
    m_loadSelectedFile->clicked().connect( this, &SpecFileQueryWidget::loadSelected );
  }
  
  m_csv = new WAnchor();
  ResultCsvResource *csvresource = new ResultCsvResource( m_resultmodel, m_csv );
  m_csv->setLink( WLink(csvresource) );
  m_csv->setTarget( TargetNewWindow );
  m_csv->setText( "csv" );
  m_csv->disable();
  
  linelayout->addWidget( m_csv, 0, 3, AlignBottom );
  
  addCondition();
  
  if( defpath.size() )
    basePathChanged();
}//void init()
  

SpecQuery::SpecFileQuery SpecFileQueryWidget::query()
{
  SpecQuery::SpecFileQuery q;
  
  const vector<WWidget *> kids = m_conditions->children();
  for( size_t i = 0; i < kids.size(); ++i )
  {
    ConditionWidget *cond = dynamic_cast<ConditionWidget *>( kids[i] );
    if( !cond )
      continue;
    
    if( cond->isCondition() )
      q.addCondition( cond->test() );
    else
      q.addLogic( cond->logic() );
  }//for( size_t i = 0; i < kids.size(); ++i )
  
  return q;
}//SpecFileQuery query()



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
  
  if( m_loadSelectedFile )
    m_loadSelectedFile->disable();
  
  //Validate logic, and give hint if something wrong, or else
  
  SpecFileQuery q = query();
  
  try
  {
    q.isvalid();
    m_numberResults->removeStyleClass( "SpecFileQueryErrorTxt" );
    m_update->enable();
    m_numberResults->setText( "Tap Update to refresh results" );
  }catch( std::exception &e )
  {
    m_numberResults->setText( e.what() );
    m_numberResults->addStyleClass( "SpecFileQueryErrorTxt" );
  }
}//void setResultsStale()


void SpecFileQueryWidget::basePathChanged()
{
  setResultsStale();
  
  m_numberFiles->setText( "" );
  m_numberResults->setText( "" );
  m_update->enable();
  
  const std::string basepath = m_baseLocation->text().toUTF8();
  if( basepath.empty() || !UtilityFunctions::is_directory( basepath ) )
  {
    m_numberResults->setText( "Not a valid base directory" );
    m_baseLocation->addStyleClass( "SpecFileQueryErrorTxt" );
    m_update->disable();
    return;
  }
  
  const bool recursive = m_recursive->isChecked();
  const bool filter = m_filterByExtension->isChecked();
  const bool filterUnique = m_filterUnique->isChecked();
  const int maxsize = m_maxFileSize->value();
  
  string prefpath;
  int prevmaxsize = 32;
  bool prevrecursive = true, prevfilter = true, prevFilterUnique = true;
  try
  {
    prevmaxsize = InterSpecUser::preferenceValue<int>( "SpecFileQueryMaxSize", m_viewer );
    prevfilter = InterSpecUser::preferenceValue<bool>( "SpecFileQueryFilter", m_viewer );
    prevrecursive = InterSpecUser::preferenceValue<bool>( "SpecFileQueryRecursive", m_viewer );
    prevFilterUnique = InterSpecUser::preferenceValue<bool>( "SpecFileQueryUnique", m_viewer );
    prefpath = InterSpecUser::preferenceValue<string>( "SpecFileQueryPath", m_viewer );
    if( prefpath == "None" )
      prefpath = "";
  }catch(...){}
  
  try
  {
    if( prevrecursive != recursive )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryRecursive", recursive, m_viewer );
    if( prevfilter != filter )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryFilter", filter, m_viewer );
    if( maxsize != prevmaxsize )
      InterSpecUser::setPreferenceValue<int>( m_viewer->m_user, "SpecFileQueryFilter", maxsize, m_viewer );
    if( filterUnique != prevFilterUnique )
      InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user, "SpecFileQueryUnique", filterUnique, m_viewer );
    if( prefpath != basepath )
      InterSpecUser::setPreferenceValue<string>( m_viewer->m_user, "SpecFileQueryPath", basepath, m_viewer );
  }catch( ... ){}
  
  m_baseLocation->removeStyleClass( "SpecFileQueryErrorTxt" );

  WServer *server = WServer::instance();  //can this ever be NULL?
  if( server )
  {
    m_numberFiles->setText( "Updating # of files" );
    const size_t maxsize_mb = static_cast<size_t>(maxsize*1024*1024);
    server->ioService().post( boost::bind( &SpecFileQueryWidget::updateNumberFiles,
                                           basepath, recursive, filter, maxsize_mb,
                                           this, wApp->sessionId(), m_widgetDeleted ) );
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
                                    std::shared_ptr< std::atomic<bool> > widgetdeleted )
{
  UtilityFunctions::file_match_function_t filterfcn = extfilter ? &maybe_spec_file : &file_smaller_than;
  
  vector<string> files;
  if( recursive )
    files = UtilityFunctions::recursive_ls( srcdir, filterfcn, (void *)&maxsize );
  else
    files = UtilityFunctions::ls_files_in_directory( srcdir, filterfcn, (void *)&maxsize );
  
  const size_t nfiles = files.size();
  
  if( !(*widgetdeleted) )
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateNumberFilesInGui,
                                                    nfiles, srcdir, recursive, extfilter, querywidget, widgetdeleted ) );
}//updateNumberFiles


void SpecFileQueryWidget::updateNumberFilesInGui( const size_t nfiles,
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
  
  const string newdir = querywidget->m_baseLocation->text().toUTF8();
  
  //Check if the user has changed pre-filter chriteria since this search was
  //  launched, and if so dont update the GUI (another search has been launched
  //  and will update things)
  if( newdir != srcdir
      || recursive != querywidget->m_recursive->isChecked()
      || extfilter != querywidget->m_filterByExtension->isChecked() )
      return;
    
  const string nfilesstr = (nfiles<99000 ? to_string(nfiles) : string(">100k"));
  querywidget->m_numberFiles->setText( (nfilesstr + " initial files").c_str() );
    
  wApp->triggerUpdate();
}//updateNumberFilesInGui(...)


void SpecFileQueryWidget::setResultFieldVisibility( const SpecQuery::FileDataField field, const bool visible )
{
  m_resultview->setColumnHidden( field, visible );
}//void setResultFieldVisibility(...)


void SpecFileQueryWidget::addConditionAt( int index )
{
  const vector<WWidget *> kids = m_conditions->children();
  
  if( index < 0 )
    index = 0;
  if( index > static_cast<int>(kids.size()) )
    index = static_cast<int>(kids.size());
  
  ConditionWidget *lastcond = kids.size() ? dynamic_cast<ConditionWidget *>( kids.back() ) : (ConditionWidget *)0;
  if( lastcond && (lastcond->isCondition() || lastcond->isCloseParan()) )
  {
    ConditionWidget *cond = new ConditionWidget();
    m_conditions->insertWidget( index, cond );
    cond->setIsAnd();
    cond->add().connect( this, &SpecFileQueryWidget::addConditionAfter );
    cond->remove().connect( this, &SpecFileQueryWidget::removeCondition );
    cond->changed().connect( this, &SpecFileQueryWidget::setResultsStale );
    index += 1;
  }//if( lastcond )
  
  ConditionWidget *cond = new ConditionWidget();
  m_conditions->insertWidget( index++, cond );
  cond->add().connect( this, &SpecFileQueryWidget::addConditionAfter );
  cond->remove().connect( this, &SpecFileQueryWidget::removeCondition );
  cond->changed().connect( this, &SpecFileQueryWidget::setResultsStale );
  
  setResultsStale();
}//void addConditionAt( int index )


void SpecFileQueryWidget::addCondition()
{
  addConditionAt( static_cast<int>( m_conditions->children().size() ) );
}//void addCondition()


void SpecFileQueryWidget::addConditionAfter( WWebWidget *ww )
{
  const int index = m_conditions->indexOf( ww );
  addConditionAt( index + 1 );
}//void addConditionAfter( WWebWidget *ww )


void SpecFileQueryWidget::removeCondition( WWebWidget *ww )
{
  if( ww && (m_conditions->indexOf(ww) >= 0) )
    delete ww;

  setResultsStale();
}//void SpecFileQueryWidget::removeCondition()


void SpecFileQueryWidget::startUpdate()
{
  const std::string basepath = m_baseLocation->text().toUTF8();
  const bool recursive = m_recursive->isChecked();
  const bool filter = m_filterByExtension->isChecked();
  const SpecQuery::SpecFileQuery q = query();
  const size_t maxsize = 1024*1024*m_maxFileSize->value();
  const bool filterUnique = m_filterUnique->isChecked();
  
  m_cancelUpdate->show();
  m_update->hide();
  m_baseLocation->disable();
  m_addCondition->disable();
  m_conditions->disable();
  m_maxFileSize->disable();
  m_filterUnique->disable();
  m_recursive->disable();
  m_filterByExtension->disable();
  m_csv->disable();
  m_resultview->disable();
  m_numberResults->setText( "Updating results" );
  if( m_resultmodel->rowCount() )
    m_resultmodel->setResult( std::shared_ptr< vector< vector<string> > >() );
    //m_resultmodel->removeRows( 0, m_resultmodel->rowCount() );
  
  unsigned long options = 0x0;
  if( recursive )
    options |= SearchRecursive;
  if( filter )
    options |= FilterByFilename;
  if( filterUnique )
    options |= FilterDuplicates;
  
  m_stopUpdate->store( false );
  WServer::instance()->ioService().post(
        boost::bind( &SpecFileQueryWidget::doSearch, this, basepath, options,
                     maxsize, q, wApp->sessionId(), m_stopUpdate,
                     m_widgetDeleted ) );
}//void startUpdate()


void SpecFileQueryWidget::finishUpdate( std::shared_ptr< std::vector< std::vector<std::string> > > result,
                                        const std::string description,
                                        const bool wasCanceled,
                                        std::shared_ptr< std::atomic<bool> > widgetDeleted )
{
  if( widgetDeleted->load() )
    return;
  
  m_cancelUpdate->hide();
  m_update->show();
  m_baseLocation->enable();
  m_maxFileSize->enable();
  m_filterUnique->enable();
  m_conditions->enable();
  m_addCondition->enable();
  m_recursive->enable();
  m_filterByExtension->enable();
  m_resultview->enable();
  m_csv->enable();
  
  if( wasCanceled )
  {
    m_csv->disable();
    m_numberResults->setText( "Search Canceled" );
    setResultsStale();
    return;
  }
  
  m_resultmodel->setHeaderData( 0, Horizontal, boost::any(WString(description)), Wt::UserRole );
  m_resultmodel->setResult( result );
  
  /*
  for( size_t i = 0; i < result->size(); ++i )
  {
    const vector<string> &row = (*result)[i];
    
    std::vector< WStandardItem * > rowitems( row.size() );
    for( size_t j = 0; j < row.size(); ++j )
    {
      if( j == SpecQuery::Filename )
      {
        rowitems[j] = new WStandardItem( WString( UtilityFunctions::filename(row[j]) ) );
        rowitems[j]->setData( WString(row[j]), Wt::UserRole );
      }else
      {
        rowitems[j] = new WStandardItem( WString(row[j]) );
      }
    }
    m_resultmodel->appendRow( rowitems );
  }
   */
  
  char buffer[256];
  snprintf( buffer, sizeof(buffer), "%i matching files", static_cast<int>(result->size()) );
  
  m_numberResults->setText( buffer );
  
  wApp->triggerUpdate();
}//finishUpdate(...)


void SpecFileQueryWidget::cancelUpdate()
{
  m_stopUpdate->store( true );
  setResultsStale();
}//cancelUpdate(...)


void SpecFileQueryWidget::updateSearchStatus( const size_t nfilestotal,
                                              const size_t nfileschecked,
                                              const size_t nfilesaccepted,
                                              const std::string specialmsg,
                                              std::shared_ptr< std::atomic<bool> > widgetDeleted )
{
  if( widgetDeleted->load() )
    return;
  
  char buffer[256];
  snprintf( buffer, sizeof(buffer), "%i of %i checked with %i passing%s",
            static_cast<int>(nfileschecked), static_cast<int>(nfilestotal),
            static_cast<int>(nfilesaccepted), specialmsg.c_str() );
  m_numberResults->setText( buffer );
  
  wApp->triggerUpdate();
}//void updateSearchStatus(...)



void testfile( const string &filename, vector<string> &result,
               const SpecQuery::SpecFileQuery &query,
               HaveSeenUuid &uniquecheck )
{
  try
  {
    result.clear();
    
    std::shared_ptr<MeasurementInfo> meas = std::make_shared<MeasurementInfo>();
    if( meas->load_file( filename, kAutoParser, filename ) )
    {
      meas->set_filename( filename );
      
      if( uniquecheck.have_seen( meas->uuid() ) )
        return;
      
      const bool testresult = query.test( meas );
      if( testresult )
        result = get_result_fields( meas );
    }
  }catch( ... )
  {
    cerr << "Caught exception testing file: " << filename << endl;
  }
}//testfile


void SpecFileQueryWidget::doSearch( const std::string basedir,
                                    unsigned long options,
                                    const size_t maxsize,
                                    const SpecQuery::SpecFileQuery query,
                                    const std::string sessionid,
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
  
  try
  {
    if( stopUpdate->load() )
      throw std::runtime_error( "" );
    
    result = std::make_shared< vector<vector<string> > >();
    
    
    typedef vector<string>(*ls_fcn_t)( const string &, UtilityFunctions::file_match_function_t, void * );
    UtilityFunctions::file_match_function_t filterfcn = extfilter ? &maybe_spec_file : &file_smaller_than;
    
    ls_fcn_t lsfcn = &UtilityFunctions::recursive_ls;
    if( !recursive )
      lsfcn = &UtilityFunctions::ls_files_in_directory;
    
    const vector<string> files = lsfcn( basedir, filterfcn, (void *)&maxsize );
    const int nfiles = static_cast<int>( files.size() );
    
    description << "There were " << nfiles << " candidate files after pre-filtering\r\n";
    
    if( stopUpdate->load() )
      throw std::runtime_error( "" );
    
    WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                      this, nfiles, 0, result->size(), "", widgetDeleted ) );
    
    double lastupdate = UtilityFunctions::get_wall_time();
    
    //We could choose something like 100 for nfile_at_a_time with the only
    //  dowsides being N42 files would parse any faster due to not enough
    //  threads being available, and also the fact that when parsing multiple
    //  file, hard drive bandwidth is probably the main limiting factor.
    const int nfile_at_a_time = SpecUtilsAsync::num_physical_cpu_cores();
    
    SpecUtilsAsync::ThreadPool pool;
  
    for( int i = 0; i < nfiles; i += nfile_at_a_time )
    {
      if( stopUpdate->load() )
        throw runtime_error("");
      
      const double now = UtilityFunctions::get_wall_time();
      if( now > (lastupdate + 1.0) )
      {
        WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                    this, nfiles, i, result->size(), "", widgetDeleted ) );
        lastupdate = now;
      }
      
      const int nfilethisone = std::min( nfiles-i, nfile_at_a_time );
      
      vector< vector<string> > testres( nfilethisone );
      
      for( int j = 0; j < nfilethisone; ++j )
        pool.post( boost::bind( &testfile, boost::cref(files[i+j]), boost::ref(testres[j]), boost::cref(query), boost::ref(uniqueCheck) ) );
      pool.join();
      
      for( int j = 0; j < nfilethisone; ++j )
        if( testres[j].size() )
          result->push_back( testres[j] );
      
      if( nfilethisone > 1000 )
        WServer::instance()->post( sessionid, boost::bind(&SpecFileQueryWidget::updateSearchStatus,
                                                          this, nfiles, i, result->size(), ", Updating GUI.", widgetDeleted ) );
    }//for( size_t i = 0; i < nfiles; ++i )
  }catch( ... )
  {
    WServer::instance()->post( sessionid, boost::bind( &SpecFileQueryWidget::finishUpdate, this, result, description.str(), true, widgetDeleted ) );
  }
  
  const double finishtime = UtilityFunctions::get_wall_time();
  const double finishcpu = UtilityFunctions::get_cpu_time();
  
  description << "Search took " << (finishtime - starttime) << " wall seconds and " << (finishcpu - startcpu) << " cpu seconds \r\n";
  description << "There were " << result->size() << " files that satisfied the test conditions\r\n";
  WServer::instance()->post( sessionid, boost::bind( &SpecFileQueryWidget::finishUpdate, this, result, description.str(), false, widgetDeleted ) );
}//doSearch(...)
  


void SpecFileQueryWidget::selectionChanged()
{
  if( !m_viewer )
    return;
  
  const bool sel = (m_resultview->selectedIndexes().size() == 1);
  m_loadSelectedFile->setEnabled( sel );
}//selectionChanged()


void SpecFileQueryWidget::loadSelected()
{
  if( !m_viewer || !m_viewer->fileManager() )
    return;
  
  const WModelIndexSet selected = m_resultview->selectedIndexes();
  if(selected.size() != 1 )
    return;
  
  const int row = selected.begin()->row();
  WModelIndex findex = m_resultmodel->index( row, SpecQuery::Filename );
  boost::any fnany = m_resultmodel->data( findex, Wt::UserRole );
  std::string filenameandy = asString( fnany ).toUTF8();
  
  const std::string filename = asString( selected.begin()->data(Wt::UserRole) ).toUTF8();
  
  if( !UtilityFunctions::is_file(filename) )
  {
    passMessage( "Sorry, '" + filename + "' file appears to no longer be accessable", "", WarningWidget::WarningMsgHigh );
    m_resultview->setSelectedIndexes( WModelIndexSet() );
    return;
  }
  
  
  const bool loaded = m_viewer->fileManager()->loadFromFileSystem( filename, kForeground, kAutoParser );
  
  if( !loaded )
    passMessage( "Sorry, '" + filename + "' can not be displayed", "", WarningWidget::WarningMsgHigh );
}//loadSelected()
  
