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

#include <map>
#include <regex>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string.hpp>

#include <Wt/WMenu>
#include <Wt/Utils>
#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WTable>
#include <Wt/WSignal>
#include <Wt/WServer>
#include <Wt/WRandom>
#include <Wt/WDateTime>
#include <Wt/WTextArea>
#include <Wt/WResource>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WIOService>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WFileUpload>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WEnvironment>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WStackedWidget>
#include <Wt/Dbo/QueryModel>
#include <Wt/WRegExpValidator>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>


#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"


#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/DrfChart.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DetectorPeakResponse.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#endif

// The regex in GCC 4.8.x does not have working regex
#if( defined(__GLIBCXX__) && (__cplusplus < 201402L) )
static_assert( defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) \
              || defined(_GLIBCXX_REGEX_STATE_LIMIT) \
              || (defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE > 4), "GCC 4.8 is not supported due to buggy regex implementation" );
#endif



using namespace std;
using namespace Wt;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


using SpecUtils::Measurement;
using SpecUtils::DetectorType;

namespace
{
  void right_select_item( WMenu *menu, WMenuItem *item )
  {
    menu->select( item );
    item->triggered().emit( item ); //
  }
  
  class UtcToLocalTimeDelegate : public Wt::WItemDelegate
  {
    int64_t m_now;
    int m_timeZoneOffset;
  public:
    UtcToLocalTimeDelegate( Wt::WObject *parent = 0 )
    : WItemDelegate( parent )
    {
      m_now = std::time(nullptr);
      if( wApp )
        m_timeZoneOffset = wApp->environment().timeZoneOffset();
      else
        m_timeZoneOffset = 0;
    }
    virtual ~UtcToLocalTimeDelegate(){}
    virtual Wt::WWidget *update( Wt::WWidget *widget,
                                const Wt::WModelIndex &index,
                                Wt::WFlags<Wt::ViewItemRenderFlag > flags )
    {
      if( flags & RenderEditing )
        throw runtime_error( "UtcToLocalTimeDelegate not for editing" );
      
      if( !(flags & RenderEditing) )
      {
        WText *text = dynamic_cast<WText *>( widget );
        
        if( !text )
          widget = text = new WText();
        
        int64_t val;
        try
        {
          val = boost::any_cast<int64_t>( index.data() );
          string valstr;
          if( val > 0 )
          {
            SpecUtils::time_point_t ptt = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::from_time_t( time_t(val) ) );
            ptt += std::chrono::seconds( 60*m_timeZoneOffset );
            valstr = SpecUtils::to_common_string(ptt, true);
          }
          
          text->setText( WString::fromUTF8( valstr ) );
        }catch( std::exception &e )
        {
          cerr << "UtcToLocalTimeDelegate caught: " << e.what() << endl;
          text->setText( "" );
          return widget;
        }//try / catch
      }//if( !(flags & RenderEditing) )
      
      widget->setStyleClass( asString( index.data(StyleClassRole) ) );
      if( flags & RenderSelected )
        widget->addStyleClass(  "Wt-selected" );

      return widget;
    }
  };//class UtcToLocalTimeDelegate
  
  
  /** A struct to hold the info a user may copy/paste into formula field. */
  struct CoefToFormulaInfo
  {
    bool m_valid_coefficients = false;   // Wether input is valid-ish user-specified coeficient list
    string m_fcn;                        // Will be empty if `m_valid_coefficients` is false
    vector<double> m_coefficients;
    string m_det_name;                   // Empty if not user specified
    double m_energy_units = 0.0;         // Will be zero if not specified, otherwise keV or MeV
    string m_detector_setback;           // Empty if not user specified
    string m_det_diam;                   // Empty if not user specified
    string m_source_to_crystal_distance; // Empty if not user specified
    float m_lower_energy = 0.0f;         // Will be zero if not user-specified
    float m_upper_energy = 0.0f;         // Will be zero if not user-specified
    bool m_fixed_geom = false;           // Currently not used, but set if known
  };//struct CoefToFormulaInfo
  
  
  /** Function to check if someone copy-pasted in coefficiecnts from a spreadsheet, instead
   of entering a mathematical formula.  If this is the case, we will assume the standard functional form
   f(x) = exp( C_0 + C_1*log(x)  + C_2*log(x)^2 + ... ).
   
   Returns true if coefficients were entered, and also changes the input string to have the proper
   mathematical function.
   */
  CoefToFormulaInfo check_if_coef_for_formula( const string &fcn )
  {
    CoefToFormulaInfo results;
    
    //  We will allow space tab, comma, semicolon, newline seperators
    vector<string> potential_coefs;
    SpecUtils::split( potential_coefs, fcn, ",\t; \n\r" );
    
    // There is an Excel sheet that users could copy-paste an entire column of
    //  data from, that would like like:
    //  { `visible?`, `Det. Name`, `cal. desc`, `Comment`, `far-field/fixed-geom`, `Det Rad(cm)`,
    //    `Source to det dist (cm)`, `Det Setback (cm)`, `Det Area`, `Source to crystal face`,
    //   `Lower Energy`, `Upper Energy`, `keV/MeV`, coefficients... }
    
    vector<string> add_info;
    for( size_t i = 0; i < potential_coefs.size(); ++i )
    {
      const bool is_keV = SpecUtils::iequals_ascii(potential_coefs[i], "keV");
      const bool is_MeV = SpecUtils::iequals_ascii(potential_coefs[i], "MeV");
      
      if( is_keV || is_MeV )
      {
        const auto first_coef_iter = begin(potential_coefs) + i;
        add_info.insert( end(add_info), begin(potential_coefs), first_coef_iter );
        potential_coefs.erase( begin(potential_coefs), first_coef_iter );
        
        // We *could* have the first entry as Yes/No, to indicate if DRF is visible or not;
        //  we will discard this entry.
        if( !add_info.empty()
            && (SpecUtils::iequals_ascii(add_info[0], "Yes")
                || SpecUtils::iequals_ascii(add_info[0], "No")) )
        {
          add_info.erase( begin(add_info) );
        }//
        
        // Now find the next "Yes" or "No" - which indicates if it is a far-field point source
        //   calibration, or if "No", then a fixed geometry DRF
        //  Before it, the name, source description, or comments can all have delimiters in them,
        //  or be empty, so we don't know how to use this info super easily.
        for( size_t i = 0; i < add_info.size(); ++i )
        {
          if( (SpecUtils::iequals_ascii(add_info[i], "Yes")
               || SpecUtils::iequals_ascii(add_info[i], "No")) )
          {
            float dummy;
            if( (i > 0) && !add_info[0].empty() &&
               !SpecUtils::parse_float(add_info[0].c_str(), add_info[0].size(), dummy) )
            {
              results.m_det_name = add_info[0];
            }
            
            add_info.erase( begin(add_info), begin(add_info) + i );
            break;
          }//
        }
        
        break;
      }//if( is_keV || is_MeV )
    }//for( int i = 0; i < static_cast<int>(potential_coefs.size()); ++i )
    
    
    if( ((add_info.size() == 9) || (add_info.size() == 10))
       && (SpecUtils::iequals_ascii( add_info[0], "No" )
           || SpecUtils::iequals_ascii( add_info[0], "Yes" )) )
    {
      results.m_fixed_geom = SpecUtils::iequals_ascii( add_info[0], "No" );
      
      double det_radius = 0, src_to_det_face = 0, det_setback = 0;
      double src_to_crystal = 0, low_energy = 0, up_energy = 0;
      
      if( SpecUtils::parse_double(add_info[1].c_str(), add_info[1].size(), det_radius) )
        results.m_det_diam = PhysicalUnits::printToBestLengthUnits( 2.0 * det_radius * PhysicalUnits::cm, 5 );
      
      if( SpecUtils::parse_double(add_info[2].c_str(), add_info[2].size(), src_to_det_face) )
        src_to_det_face *= PhysicalUnits::cm;
      
      // 3 - source radius
      
      if( SpecUtils::parse_double(add_info[4].c_str(), add_info[4].size(), det_setback) )
        results.m_detector_setback = add_info[4] + " cm";
      
      // 5 - Detector Area
      
      if( SpecUtils::parse_double(add_info[6].c_str(), add_info[6].size(), src_to_crystal) )
        results.m_source_to_crystal_distance = add_info[6] + " cm";
      
      // 7 - Relative Efficiency (%) - doesnt seem to
      
      // Second to last element, lower energy
      const string &low_energy_str = add_info[add_info.size()-2];
      if( SpecUtils::parse_double(low_energy_str.c_str(), low_energy_str.size(), low_energy) )
        results.m_lower_energy = low_energy * PhysicalUnits::keV;
      
      // Last element, upper energy
      const string &up_energy_str = add_info[add_info.size()-1];
      if( SpecUtils::parse_double(up_energy_str.c_str(), up_energy_str.size(), up_energy) )
        results.m_upper_energy = up_energy * PhysicalUnits::keV;
      
      if( low_energy >= up_energy )
        results.m_lower_energy = results.m_upper_energy = 0.0;
    }//if( add_info_fields.size() == 12 )
    
    size_t last_non_zero_coef = 0;
    vector<double> potential_coef_values;
    try
    {
      for( auto iter = begin(potential_coefs); iter != end(potential_coefs); ++iter )
      {
        const bool is_keV = SpecUtils::iequals_ascii(*iter, "keV");
        const bool is_MeV = SpecUtils::iequals_ascii(*iter, "MeV");
        if( is_keV || is_MeV )
        {
          results.m_energy_units = is_keV ? PhysicalUnits::keV : PhysicalUnits::MeV;
          potential_coefs.erase( iter );
          break;
        }
      }//for( size_t i = 0; i < potential_coefs.size(); ++i )
      
      
      for( size_t i = 0; i < potential_coefs.size(); ++i )
      {
        const string &txt = potential_coefs[i];
        size_t idx = 0;
        const double val = stod( txt, &idx ); //throws exception if invalid float
        
        // Dont accept inf, NaN, or allow trailing characters
        if( IsInf(val) || IsNan(val) || (idx != txt.size()) )
          throw exception();
        
        potential_coef_values.push_back( val );
        if( val != 0.0 )
          last_non_zero_coef = i;
      }//
    }catch( std::exception & )
    {
      // If we're here, its an invalid double, or there where characters after the number
      results.m_valid_coefficients = false;
      return results;
    }//try catch to parse entered text as exactly a list of coefficecnts
    
    
    if( (potential_coef_values.size() <= 1) || (last_non_zero_coef <= 1) )
    {
      results.m_valid_coefficients = false;
      return results;
    }
    
    // If we are here, we are quite confident the user entered a series of coeficients.
    assert( potential_coefs.size() == potential_coef_values.size() );
    assert( !potential_coefs.empty() && (last_non_zero_coef < potential_coefs.size()) );
      
    potential_coefs.resize( last_non_zero_coef + 1 );
    potential_coef_values.resize( last_non_zero_coef + 1 );
    
    results.m_coefficients = potential_coef_values;
    
    results.m_fcn = "";
    
    for( size_t i = 0; i < potential_coefs.size(); ++i )
    {
      if( potential_coef_values[i] == 0.0 )
        continue;
        
      const string mfunc = (i > 0) ? "*log(x)" : "";
      const string power = (i > 1) ? ("^" + to_string(i)) : string();
      results.m_fcn += (results.m_fcn.empty() ? "exp( " : " + ") + potential_coefs[i] + mfunc + power;
    }
    results.m_fcn += results.m_fcn.empty() ? "" : " )";
    
    results.m_valid_coefficients = true;
    
    return results;
  }//bool check_if_coef_for_formula( string &fcn )

  
  //Allow the file to be comma or tab delimited, but if an individual field
  //  contains a comma or tab, then the field must be quoted by a double
  //  quote.  Note that if you just copy cells from Microsoft Excel, that
  //  contain a comma, and then past into a text editor, fields with a comma
  //  are not quoted.
  void split_escaped_csv( vector<string> &fields, const string &line )
  {
    fields.clear();
    
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokeniser;
    boost::escaped_list_separator<char> separator("\\",",\t", "\"");
    Tokeniser t( line, separator );
    for( Tokeniser::iterator it = t.begin(); it != t.end(); ++it )
      fields.push_back(*it);
  }//void split_escaped_csv(...)
  
 
  //Eliminate any duplicate strings in the vector.
  //  In the future we may decide to remove duplicates based just on the `leaf` file name
  void remove_duplicate_paths( vector<string> &paths )
  {
    // Note that the following requires sorting paths first, (eg std::sort(begin(paths),end(paths)))
    //  because std::unique only makes neighboring elements unique
    //  paths.erase( std::unique( begin(paths), end(paths) ), end(paths) );
  
    const auto paths_begin = begin(paths);
    const auto paths_end = end(paths);
  
    auto end_pos = paths_begin;
  
    for( auto iter = paths_begin; iter != paths_end; ++iter )
    {
      if( std::find(paths_begin, end_pos, *iter) == end_pos )
      {
        if( end_pos != iter )
          *end_pos = *iter;
        ++end_pos;
      }
    }//for( auto iter = paths_begin; iter != paths_end; ++iter )
  
    paths.erase( end_pos, end(paths) );
  }//remove_duplicate_paths(...)


 /** Searches static and user-writable data directories for the given file.
  */
 string find_valid_path( string pathstr, const bool isdir )
 {
#if( BUILD_FOR_WEB_DEPLOYMENT )
   if( SpecUtils::is_absolute_path(pathstr) )
     return "";
   
   SpecUtils::ireplace_all(pathstr, "..", "");
   pathstr = SpecUtils::lexically_normalize_path(pathstr);
   return pathstr;
#endif
   
   if( SpecUtils::is_absolute_path(pathstr) )
     return pathstr;
   
   const auto file_check = isdir ? &SpecUtils::is_directory : &SpecUtils::is_file;
   
   if( file_check(pathstr) )
     return pathstr;
   
   vector<string> trial_paths;
   
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
   try
   {
     const string user_data_dir = InterSpec::writableDataDirectory();
     if( !user_data_dir.empty() )
     {
       trial_paths.push_back( user_data_dir );
       trial_paths.push_back( SpecUtils::append_path( user_data_dir, "drfs" ) );
     }
   }catch( std::exception & )
   {
   }
#endif
   
   string trialpath = InterSpec::staticDataDirectory();
   trial_paths.push_back( trialpath );
   trialpath = SpecUtils::lexically_normalize_path( SpecUtils::append_path( trialpath, "..") );
   trial_paths.push_back( trialpath );
   
   trialpath = SpecUtils::append_path( SpecUtils::append_path( trialpath, ".."), "data_ouo" );
   trialpath = SpecUtils::lexically_normalize_path( trialpath );
   trial_paths.push_back( trialpath );
   
   
   for( const string &p : trial_paths )
   {
     trialpath = SpecUtils::append_path(p, pathstr);
     if( file_check(trialpath) )
       return trialpath;
   }//for( const string &p : trial_paths )
   
   return pathstr;
  }//find_valid_path


  string complete_drf_path( string pathstr )
  {
    return find_valid_path( pathstr, false );
  }


  string complete_drf_dir_path( string pathstr )
  {
    return find_valid_path( pathstr, true );
  }
  
  void potential_rel_eff_files( const string &base_dir, vector<string> &potential_paths )
  {
    const vector<string> tsv_files = SpecUtils::recursive_ls( base_dir, ".tsv");
    for( const auto &p : tsv_files )
      potential_paths.push_back( p );
    
    const string drfsdir = SpecUtils::append_path( base_dir, "drfs" );
    
    const vector<string> drf_files = SpecUtils::recursive_ls( drfsdir );
    for( const auto &p : drf_files )
    {
      const string fname = SpecUtils::filename(p);
      if( SpecUtils::iequals_ascii( fname, "Efficiency.csv") )
        continue;
      
      const string type = SpecUtils::file_extension(fname);
      if( !SpecUtils::iequals_ascii( type, ".csv") && !SpecUtils::iequals_ascii( type, ".tsv") )
        continue;
      
      // Note that this next line will cause duplicate TSV entries - we'll remove them later
      potential_paths.push_back( p );
    }//for( const auto &p : drf_files )
  }//potential_rel_eff_files(..)

/*
#if( PERFORM_DEVELOPER_CHECKS )
void check_url_serialization( std::shared_ptr<const DetectorPeakResponse> drf )
{
  if( !drf || !drf->isValid() )
    return;
  
  string appUrl;
  try
  {
    appUrl = drf->toAppUrl();
    DetectorPeakResponse rhs;
    
    rhs.fromAppUrl( appUrl );
    
    DetectorPeakResponse::equalEnough( *drf, rhs );
  }catch( std::exception &e )
  {
    cerr << "Failed to decode DRF from URL: " << e.what() << endl;
    assert( 0 );
    throw runtime_error( "Failed to decode DRF from URL: " + string(e.what()) );
  }
  
  cout << "Detector '" << drf->name() <<  "' has initial DRF len=" << appUrl.size()
       << " and passes going to and then from URL" << endl;
}//void check_url_serialization( std::shared_ptr<const DetectorPeakResponse> drf )
#endif //#if( PERFORM_DEVELOPER_CHECKS )
*/

class DrfDownloadResource : public Wt::WResource
{
  DrfSelect * const m_drfSelect;
  Wt::WApplication *m_app;
  
public:
  DrfDownloadResource( DrfSelect *drfSelect )
  : WResource( drfSelect ),
    m_drfSelect( drfSelect ),
    m_app( WApplication::instance() )
  {
    assert( m_app );
    assert( m_drfSelect );
  }
  
  virtual ~DrfDownloadResource()
  {
    beingDeleted();
  }
  
  virtual void handleRequest( const Wt::Http::Request &request,
                             Wt::Http::Response &response )
  {
    WApplication::UpdateLock lock( m_app );
    
    if( !lock )
    {
      log("error") << "Failed to WApplication::UpdateLock in DrfDownloadResource.";
      
      response.out() << "Error grabbing application lock to form DrfDownloadResource resource;"
                        " please report to InterSpec@sandia.gov.";
      response.setStatus(500);
      assert( 0 );
      
      return;
    }//if( !lock )
    
    if( !m_drfSelect )
      return;
    
    std::shared_ptr<DetectorPeakResponse> det = m_drfSelect->detector();
    if( !det || !det->isValid() )
      return;
    
    string name = det->name();
    if( !name.empty() )
      name += ".";
    name += "drf.xml";
    
    //Remove bad filename characters
    const string notallowed = "\\/:?\"<>|*";
    for( auto it = begin(name) ; it < end(name) ; ++it )
    {
      if( notallowed.find(*it) != string::npos )
        *it = '_';
    }
    
    suggestFileName( name, WResource::Attachment );
    response.setMimeType( "application/octet-stream" );
    
    try
    {
      rapidxml::xml_document<char> doc;
      det->toXml( &doc, &doc );
      rapidxml::print( response.out(), doc, 0 );
    }catch( std::exception &e )
    {
      log("error") << "Failed to to writ XML for DRF: " << e.what();
      
      response.out() << "Error creating DRF XML.";
      response.setStatus(500);
      return;
    }
  }
};//class DrfDownloadResource

}//namespace


class RelEffFile;

enum class MatchDetectorStatus
{
  NotInited,
  Match,
  NoMatch
};//


/** This class represents potentially many CSV or TSV relative efficiency
    detector response files (that each may have multiple DRFs).
    Looks in the InterSpec::writableDataDirectory() and InterSpec::staticDataDirectory() directories for TSVs/CSVs
    This widget allows users to add or remove files.
 */
class RelEffDetSelect : public Wt::WContainerWidget
{
  friend class RelEffFile;
  
protected:
  InterSpec *m_interspec;
  DrfSelect *m_drfSelect;
  WContainerWidget *m_files; //holds the RelEffFile objects..
  
public:
  RelEffDetSelect( InterSpec *interspec, DrfSelect *detedit, WContainerWidget *parent = 0 );
  
  virtual ~RelEffDetSelect(){}
  
  void userSelectedRelEffDetSelect();
  
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  void addFile();
#endif
  
  MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  void detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det );
  
  void docreate();
  
  //Attempts to do a defered rednering, but it looks like thisisnt actually
  //  the case
  virtual void load();
};//class RelEffDetSelect


/** Represents a single CSV/TSV DRF file, which may have multiple DRFs.
 */
class RelEffFile : public Wt::WContainerWidget
{
  void init();
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  void handleUserAskedRemove();
#endif
  
  std::string m_existingFilePath; //Will be non-empty only if constructor with a path name is used, or user has asked to save the file for later
  WFileUpload *m_fileUpload; //Will be nullptr if constructor with a string is used

  RelEffDetSelect *m_relEffDetSelect;
  DrfSelect *m_drfSelect;
  
  Wt::WText *m_credits;
  Wt::WComboBox *m_detectorSelect;
  std::vector<std::shared_ptr<DetectorPeakResponse> > m_responses;
  
public:
  /** Constructor to point to an existing file on disk. */
  RelEffFile( std::string file,
             RelEffDetSelect *parentSelect,
             DrfSelect *drfSelect,
             WContainerWidget *parent );
  
  /** Constructor that will create a file upload area for user to upload a file */
  RelEffFile( RelEffDetSelect *parentSelect,
             DrfSelect *drfSelect,
             WContainerWidget *parent );
  
  
  ~RelEffFile(){}
  
  void selectNone();
  
  void detectorSelectCallback();
  
  MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  void initDetectors();
  void detectorSelected( const int index );
  
  void handleFileUpload();

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  void handleSaveFileForLater();
#endif
  
  const std::vector<std::shared_ptr<DetectorPeakResponse> > &availableDrfs();
};//class RelEffFile


class GadrasDirectory : public Wt::WContainerWidget
{
  friend class GadrasDetSelect;
  
protected:
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  std::string m_directory;
#else
  WLineEdit *m_directoryEdit;
  WPushButton *m_setDirectoryButton;
#endif
  GadrasDetSelect *m_parent;
  DrfSelect *m_drfSelect;
  
  Wt::WComboBox *m_detectorSelect;
  std::vector<std::shared_ptr<DetectorPeakResponse> > m_responses;

  Wt::WText *m_msg;
  Wt::WPushButton *m_deleteBtn;
  
  GadrasDirectory( std::string dorectory, GadrasDetSelect *parentSelect,
                   DrfSelect *drfSelect, WContainerWidget *parent );
  
  
  ~GadrasDirectory(){}
  
  std::string directory();
  
  void selectNone();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  void dirPathChanged();
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )
  
  void detectorSelectCallback();
  
  MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  //parseDetector() returns null on error
  static std::shared_ptr<DetectorPeakResponse> parseDetector( const string directory );
  static vector<string> recursive_list_gadras_drfs( const string &sourcedir );
  
  void initDetectors();
  void detectorSelected( const int index );
};//class GadrasDirectory


/** Handles selecting a GADRAS DRF
 */
class GadrasDetSelect : public Wt::WContainerWidget
{
protected:
  InterSpec *m_interspec;
  DrfSelect *m_drfSelect;
  
  WContainerWidget *m_directories; //holds the GadrasDirectory objects..
  
public:
  GadrasDetSelect( InterSpec *interspec, DrfSelect *detedit, WContainerWidget *parent = 0 );
  
  virtual ~GadrasDetSelect(){}
  
  MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det );
  
  std::shared_ptr<DetectorPeakResponse> selectedDetector();
  
  void addDirectory();
  
  void removeDirectory( GadrasDirectory *dirWidget );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  void saveFilePathToUserPreferences();
#endif
  
  void detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det );
  
  void docreate();
  
  //Attempts to do a defered rednering, (not actually verified)
  virtual void load();
};//class GadrasDetSelect


////////////////////////////////////////////////////////////////////////////////
//////////////////   Begin RelEffFile implementation   /////////////////////////
////////////////////////////////////////////////////////////////////////////////



RelEffFile::RelEffFile( std::string file,
                       RelEffDetSelect *parentSelect,
                       DrfSelect *drfSelect,
                       WContainerWidget *parent )
: WContainerWidget( parent ),
  m_existingFilePath( file ),
  m_fileUpload( nullptr ),
  m_relEffDetSelect( parentSelect ),
  m_drfSelect( drfSelect ),
  m_credits( nullptr )
{
  init();
}//RelEffFile


RelEffFile::RelEffFile( RelEffDetSelect *parentSelect,
                        DrfSelect *drfSelect,
                        WContainerWidget *parent )
: WContainerWidget( parent ),
m_existingFilePath(),
m_fileUpload( nullptr ),
m_relEffDetSelect( parentSelect ),
m_drfSelect( drfSelect ),
m_credits( nullptr )
{
  init();
}//RelEffFile


void RelEffFile::init()
{
  addStyleClass( "RelEffFile" );
  
  WContainerWidget *topdiv = new WContainerWidget( this );

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  WPushButton *closeIcon = new WPushButton( topdiv );
  closeIcon->addStyleClass( "closeicon-wtdefault" );
  closeIcon->setToolTip( WString::tr("ref-tt-remove-drf-file") );
  //closeIcon->setAttributeValue( "style", "position: relative; top: 3px; right: 3px;" + closeIcon->attributeValue("style") );
  closeIcon->clicked().connect( this, &RelEffFile::handleUserAskedRemove );
#endif
  
  if( m_existingFilePath.size() )
  {
    string user_data_dir;
    try
    {
      if( !SpecUtils::is_absolute_path(m_existingFilePath)
          && !SpecUtils::is_file(m_existingFilePath) )
      {
        m_existingFilePath = complete_drf_path( m_existingFilePath );
      }
    }catch( std::exception & )
    {
    }
    
    string filename = SpecUtils::filename( m_existingFilePath );
    WText *filenameDisplay = new WText( WString::fromUTF8(filename), topdiv );
    filenameDisplay->addStyleClass( "RelEffFixedPath" );
  }else
  {
    m_fileUpload = new WFileUpload( topdiv );
    m_fileUpload->setInline( false );
    m_fileUpload->uploaded().connect( this, &RelEffFile::handleFileUpload );
    m_fileUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge,
                                                      boost::placeholders::_1) );
    m_fileUpload->changed().connect( m_fileUpload, &WFileUpload::upload );
  }//if( m_existingFilePath.size() )
  
  
  
  WContainerWidget *bottomDiv = new WContainerWidget( this );
  
  m_detectorSelect = new WComboBox( bottomDiv );
  m_detectorSelect->setInline( false );
  m_detectorSelect->setMaximumSize( WLength(90,WLength::Percentage), WLength::Auto );
  m_detectorSelect->setNoSelectionEnabled( true );
  m_detectorSelect->activated().connect( this, &RelEffFile::detectorSelectCallback );
  
  
  m_credits = new WText( "", Wt::XHTMLText, bottomDiv );
  m_credits->setInline( false );
  m_credits->addStyleClass( "RelEffFileCredits" );
  
  initDetectors();
}//void RelEffFile::init()

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
void RelEffFile::handleUserAskedRemove()
{
  if( m_existingFilePath.empty() )
  {
    delete this;
    return;
  }
  
  SimpleDialog *dialog = new SimpleDialog( WString::tr("ref-remove-drf-window-title"),
                                          WString::tr("ref-remove-drf-window-txt") );
  
  Wt::WPushButton *yes = dialog->addButton( WString::tr("Yes") );
  dialog->addButton( WString::tr("No") );
  
  string filepath = m_existingFilePath;
  yes->clicked().connect( std::bind([filepath](){
    
    string pathToDel = SpecUtils::lexically_normalize_path(filepath);
    
    if( !SpecUtils::is_absolute_path(pathToDel) )
    {
      passMessage( "Error removing TSV or CSV relative-efficiency DRF file - not an absolute path.",
                  WarningWidget::WarningMsgHigh );
      return;
    }
    
    
    try
    {
      const string userDir = InterSpec::writableDataDirectory();
      if( !SpecUtils::starts_with(pathToDel, userDir.c_str()) )
        throw runtime_error( "Not in users data path." );
    }catch( std::exception &e )
    {
      passMessage( "Error removing TSV or CSV relative-efficiency DRF file - file path is in unexpected base.",
                  WarningWidget::WarningMsgHigh );
      return;
    }//
    
    
    try
    {
      boost::filesystem::remove(pathToDel);
      cout << "Removed '" << filepath << "'" << endl;
    }catch( std::exception &e )
    {
      cerr << "Error removing file '" << filepath << "': " << e.what() << endl;
      passMessage( "Error removing TSV or CSV relative-efficiency DRF file - invalid path.",
                  WarningWidget::WarningMsgHigh );
      return;
    }
    
    
    passMessage( WString::tr("ref-file-removed-toast").arg(SpecUtils::filename(pathToDel)),
                WarningWidget::WarningMsgInfo );
  }) );
  
  
  delete this;
}//void RelEffFile::handleUserAskedRemove()

void RelEffFile::handleSaveFileForLater()
{
  if( !m_fileUpload )
  {
    assert( 0 ); //shouldnt ever happen
    return;
  }
  
  if( m_fileUpload->empty() || m_responses.empty() )
  {
    passMessage( "No valid file available for saving.", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const string spoolPath = m_fileUpload->spoolFileName();
  const string displayName = m_fileUpload->clientFileName().toUTF8();
  
  try
  {
    std::string datadir = InterSpec::writableDataDirectory();
    if( datadir.empty() )
      throw runtime_error( "Writable data directory not set." );
    
    datadir = SpecUtils::append_path( datadir, "drfs" );
    
    if( SpecUtils::create_directory(datadir) == 0 ) //-1 means already existed, 1 means created
      throw runtime_error( "Could not create 'drfs' directory in app data directory." );
    
    //displayName
    string filename = SpecUtils::filename( displayName );
    const string orig_extension = SpecUtils::file_extension( filename );
    assert( orig_extension.size() <= filename.size() );
    
    if( orig_extension.size() )
      filename = filename.substr( 0, filename.size() - orig_extension.size() );
    
    const int offset = wApp->environment().timeZoneOffset();
    auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    now += std::chrono::seconds(60*offset);
    string timestr = SpecUtils::to_vax_string(now); //"2014-Sep-19 14:12:01.62"
    const string::size_type pos = timestr.find( ' ' );
    //std::string timestr = SpecUtils::to_extended_iso_string( now ); //"2014-04-14T14:12:01.621543"
    //string::size_type pos = timestr.find( 'T' );
    if( pos != string::npos )
      timestr = timestr.substr(0,pos);
    SpecUtils::ireplace_all( timestr, "-", "_" );
    
    filename += "_" + timestr + orig_extension;
    const string outputname = SpecUtils::append_path( datadir, filename );
    
    boost::filesystem::copy( spoolPath, outputname );
    
    m_existingFilePath = outputname;
    
    passMessage( WString::tr("ref-file-saved-for-later").arg(filename),
                WarningWidget::WarningMsgInfo );
  }catch( std::exception &e )
  {
    // Dont ever expect to really get here.
    cerr << "RelEffFile::handleSaveFileForLater(): Caught exception trying to save '"
         << displayName << "' from '" << spoolPath << "'." << endl;
    passMessage( WString::tr("ref-err-saving-drf-file"), WarningWidget::WarningMsgHigh );
  }//try / catch
}//void handleSaveFileForLater();
#endif  //#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )


void RelEffFile::handleFileUpload()
{
  selectNone();

  initDetectors();
  
  if( m_responses.empty() )
  {
    passMessage( WString::tr("ref-err-invalid-csv-drf"), WarningWidget::WarningMsgHigh );
    return;
  }//if( invalid file )
  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  SimpleDialog *dialog = new SimpleDialog( WString::tr("ref-save-drf-window-title"),
                                WString::tr("ref-save-drf-window-txt") );
  
  Wt::WPushButton *yes = dialog->addButton( WString::tr("Yes") );
  dialog->addButton( WString::tr("No") );
  yes->clicked().connect( this, &RelEffFile::handleSaveFileForLater );
#endif
}//void handleFileUpload()



void RelEffFile::selectNone()
{
  if( m_detectorSelect && m_detectorSelect->count() )
    m_detectorSelect->setCurrentIndex( 0 );
}//void selectNone()


void RelEffFile::detectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det;
 
  const int index = m_detectorSelect->currentIndex();
  
  //Item at index 0 is always a non-detector.
  if( index > 0 && (index-1) < m_responses.size() )
    det = m_responses[index-1];
  
  if( m_relEffDetSelect )
    m_relEffDetSelect->detectorSelected( this, det );
  
  if( m_drfSelect )
  {
    m_drfSelect->setDetector( det );
    m_drfSelect->emitChangedSignal();
  }//if( m_drfSelect )
}//void detectorSelectCallback()


MatchDetectorStatus RelEffFile::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( !det )
  {
    m_detectorSelect->setCurrentIndex( 0 );
    return MatchDetectorStatus::NoMatch;
  }//if( !det )
  
  for( size_t i = 0; i < m_responses.size(); ++i )
  {
    if( m_responses[i]->hashValue() == det->hashValue() )
    {
      m_detectorSelect->setCurrentIndex( static_cast<int>(i) + 1 );
      return MatchDetectorStatus::Match;
    }
  }//for( size_t i = 0; i < m_responses.size(); ++i )
    
  m_detectorSelect->setCurrentIndex( 0 );
  
  if( m_detectorSelect->count() <= 1 )
    return MatchDetectorStatus::NotInited;
  
  return MatchDetectorStatus::NoMatch;
}//void trySelectDetector(...)


void RelEffFile::initDetectors()
{
  m_responses.clear();
  m_credits->setText( "" );
  while( m_detectorSelect->count() )
    m_detectorSelect->removeItem( 0 );
  
  vector<string> credits;
  string pathstr;
  
  if( m_fileUpload )
  {
    if( m_fileUpload->empty() )
    {
      m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("ref-dont-use-drf-in-file") ) );
      m_detectorSelect->setCurrentIndex( 0 );
      m_detectorSelect->hide();
      m_credits->setText( WString::tr("ref-upload-csv-file") );
      
      return;
    }//if( m_fileUpload->empty() )
    
    pathstr = m_fileUpload->spoolFileName();
  }else
  {
    pathstr = m_existingFilePath;
  }

#ifdef _WIN32
  const std::wstring wpathstr = SpecUtils::convert_from_utf8_to_utf16(pathstr);
  std::ifstream input( wpathstr.c_str(), ios::in | ios::binary );
#else
  std::ifstream input( pathstr.c_str(), ios::in | ios::binary );
#endif
  
  const bool file_opened = input.is_open();
  
  if( file_opened )
  {
    DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, m_responses );
    
//#if( PERFORM_DEVELOPER_CHECKS )
//    for( auto drf : m_responses )
//      check_url_serialization( drf );
//#endif
    
    if( m_responses.empty() )
    {
      m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("ref-dont-use-drf-in-file") ) );
      m_detectorSelect->hide();
    }else
    {
      m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("ref-select-det-from-file") ) );
      m_detectorSelect->show();
    }
    
    if( m_responses.empty() )
      credits.push_back( WString("<span style=\"color:red;\">{1}</span>").arg( WString::tr("ref-no-drfs-in-file") ).toUTF8() );
    
    for( const auto &det : m_responses )
      m_detectorSelect->addItem( det->name() );
  }else
  {
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
    credits.push_back( WString("<span style=\"color:red;\">{1}</span>").arg( WString::tr("ref-no-rel-eff-drfs-in-file") ).toUTF8() );
#else
    credits.push_back( WString("<span style=\"color:red;\">{1}</span>").arg( WString::tr("ref-couldnt-open-drf-file") ).toUTF8() );
#endif
    m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("ref-no-drfs-available") ) );
    m_detectorSelect->hide();
  }//if( file_opened )
  
  if( file_opened && m_responses.empty() )
    credits.push_back( WString::tr("<span style=\"color:red;\">{1}</span>").arg( WString::tr("ref-file-has-no-drfs") ).toUTF8() );
  
  string creditHtml;
  for( string credit : credits )
  {
    WString c = WString::fromUTF8(credit);
    Wt::Utils::removeScript(c);
    creditHtml += "<div>" + c.toUTF8() + "</div>";
  }
  
  m_credits->setText( creditHtml );
  
  // Call DrfSelect::detectorsWereInited() since we may be updating from a worker thread;
  //  the main thread will call later anyway.
  m_drfSelect->detectorsWereInited();
  
  m_detectorSelect->setCurrentIndex( 0 );
}//initDetectors()



void RelEffFile::detectorSelected( const int index )
{
  std::shared_ptr<DetectorPeakResponse> det;
  if( index > 0 )
    det = m_responses.at( index - 1 );
  
  if( !m_drfSelect )
    return;
  
  m_drfSelect->setDetector( det );
  m_drfSelect->emitChangedSignal();
}//void detectorSelected( const int index )





const std::vector<std::shared_ptr<DetectorPeakResponse> > &RelEffFile::availableDrfs()
{
  return m_responses;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////   End RelEffFile implementation   //////////////////////////
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////   Begin RelEffDetSelect implementation   //////////////////////
////////////////////////////////////////////////////////////////////////////////
RelEffDetSelect::RelEffDetSelect( InterSpec *interspec, DrfSelect *detedit, WContainerWidget *parent )
    : WContainerWidget( parent ),
      m_interspec( interspec ),
      m_drfSelect( detedit ),
      m_files( nullptr )
{
}


#if( !BUILD_FOR_WEB_DEPLOYMENT )
void RelEffDetSelect::addFile()
{
  if( m_files )
    new RelEffFile( this, m_drfSelect, m_files );
}//void addFile()


#endif //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void RelEffDetSelect::detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det )
{
  if( !m_files || !file )
    return;
  
  for( WWidget *w : m_files->children() )
  {
    auto child = dynamic_cast<RelEffFile *>(w);
    if( child && (child != file) )
      child->selectNone();
  }//for( WWidget *w : m_files->children() )
}//void detectorSelected( RelEffFile *file, std::shared_ptr<DetectorPeakResponse> det );


void RelEffDetSelect::userSelectedRelEffDetSelect()
{
  
}//void userSelectedRelEffDetSelect()


MatchDetectorStatus RelEffDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  docreate();
  
  if( !m_files )
    return MatchDetectorStatus::NoMatch;
  
  auto children = m_files->children();
  
  bool anyNotInitied = false, found = false;
  for( auto w : children )
  {
    auto child = dynamic_cast<RelEffFile *>( w );
    if( !child )
      continue;
    
    if( found )
    {
      child->selectNone();
    }else
    {
      const MatchDetectorStatus status = child->trySelectDetector(det);
      found = (status == MatchDetectorStatus::Match);
      anyNotInitied |= (status == MatchDetectorStatus::NotInited);
    }
  }//for( size_t i = 0; i < children.size(); ++i )
  
  if( found )
    return MatchDetectorStatus::Match;
  
  return anyNotInitied ? MatchDetectorStatus::NotInited : MatchDetectorStatus::NoMatch;
}//void trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


void RelEffDetSelect::docreate()
{
  if( m_files )
    return;
  
  addStyleClass( "DrfSelectArea" );
  
  m_files = new WContainerWidget( this );
  m_files->addStyleClass( "DrfSelectAreaFiles" );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  WContainerWidget *holder = new WContainerWidget( this );
  holder->addStyleClass( "DrfSelectAreaFooter" );
  //holder->setToolTip( WString::tr("reds-tt-add-drf-btn") );
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  HelpSystem::attachToolTipOn( holder, WString::tr("reds-tt-add-drf-btn"),
                              showToolTips, HelpSystem::ToolTipPosition::Left );
  
  WPushButton *addIcon = new WPushButton( holder );
  addIcon->setStyleClass( "DrfSelectAddFile Wt-icon" );
  addIcon->clicked().connect( this, &RelEffDetSelect::addFile );
#endif
  
  vector<string> user_data_paths;

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  try
  {
    const string userDir = InterSpec::writableDataDirectory();
    potential_rel_eff_files( userDir, user_data_paths );
  }catch( std::exception & )
  {
    cerr << "Couldnt call into InterSpec::writableDataDirectory()" << endl;
  }
#endif
  
  try
  {
    const string dataDir = InterSpec::staticDataDirectory();
    potential_rel_eff_files( dataDir, user_data_paths );
  }catch( std::exception & )
  {
    cerr << "Couldnt call into InterSpec::staticDataDirectory()" << endl;
  }

  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || ANDROID || IOS )
  if( !InterSpecApp::isPrimaryWindowInstance() )
  {
    // TODO: decide if we want to do anything here
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif

  
  remove_duplicate_paths( user_data_paths );
  
  size_t num_user = 0;
  for( const string &path : user_data_paths )
  {
    RelEffFile *f = new RelEffFile( path, this, m_drfSelect, m_files );
    
    //We'll remove the widget if the DRF is invalid, or make it so path cant be changed for valid
    //  files (this then doesnt provide any feedback to the user that a DRF they placed in their
    //  data directory is invalid, but I'm not sure if dealing with that here would be a net
    //  improvement vs the alternative of possibly showing lots of files not intended to be DRF
    //  files)
    //
    const vector<shared_ptr<DetectorPeakResponse>> drfs = f->availableDrfs();
    if( drfs.empty() )
    {
      delete f;
    }else
    {
      num_user += 1;
    }
  }//for( const string &path : user_data_paths )
  
  if( !num_user )
  {
    new RelEffFile( this, m_drfSelect, m_files );
  }
}//void docreate()


// Attempts to do a defered rendering.  Note that if trySelectDetector(...) has
//  already been called, then so will have docreate().
void RelEffDetSelect::load()
{
  if( !loaded() )
    docreate();
  WContainerWidget::load();
}//void load()
////////////////////////////////////////////////////////////////////////////////
/////////////////   End RelEffDetSelect implementation   ///////////////////////
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
////////////////   Begin GadrasDetSelect implementation   //////////////////////
////////////////////////////////////////////////////////////////////////////////


///************************ Begin Need to implement for GADRAS **************/////
GadrasDetSelect::GadrasDetSelect( InterSpec *interspec, DrfSelect *detedit, WContainerWidget *parent )
  : WContainerWidget(parent),
    m_interspec( interspec ),
    m_drfSelect( detedit ),
    m_directories( nullptr )
{
  //GADRAS-DRF URL: https://rsicc.ornl.gov/codes/psr/psr6/psr-610.html
}


void GadrasDetSelect::load()
{
  if( !loaded() )
    docreate();
  WContainerWidget::load();
}//void load()


void GadrasDetSelect::docreate()
{
  if( m_directories )
    return;
  
  addStyleClass( "DrfSelectArea" );
  
  m_directories = new WContainerWidget( this );
  m_directories->addStyleClass( "DrfSelectAreaFiles" );
  
  //Need to point the GUI to the appropriate directory, and implement to an `ls` to find detctors with both Detcotr.dat and Efficy.csv.
  const string drfpaths = UserPreferences::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  
  WContainerWidget *footer = new WContainerWidget( this );
  footer->addStyleClass( "DrfSelectAreaFooter" );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  //footer->setToolTip( WString::tr("reds-tt-add-dir") );
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  HelpSystem::attachToolTipOn( footer, WString::tr("reds-tt-add-dir"), showToolTips, HelpSystem::ToolTipPosition::Left );
  
  
  WPushButton *addIcon = new WPushButton( footer );
  addIcon->setStyleClass( "DrfSelectAddFile Wt-icon" );
  addIcon->clicked().connect( this, &GadrasDetSelect::addDirectory );
#endif
  
  string pathstr;
  vector<string> paths;
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  try
  {
    if( m_interspec )
      pathstr = UserPreferences::preferenceValue<string>( "GadrasDRFPath", m_interspec );
  }catch( std::exception & )
  {
    passMessage( "Error retrieving 'GadrasDRFPath' preference.", WarningWidget::WarningMsgHigh );
  }
#endif

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  try
  {
    //ToDo: do more testing and use only the Android implementation
#if( ANDROID )
    const string basestr = InterSpec::writableDataDirectory();
    const vector<string> subdirs = SpecUtils::ls_directories_in_directory( basestr );
    for( const string &subdirpath : subdirs )
    {
      auto subsubdirs = GadrasDirectory::recursive_list_gadras_drfs(subdirpath);
      if( subsubdirs.size() )
        pathstr += (pathstr.empty() ? "" : ";") + subdirpath;
    }//for( const string &subdir : subdirs )
#else
    using namespace boost::filesystem;
    auto itr = directory_iterator( InterSpec::writableDataDirectory() );
    
    for( ; itr != directory_iterator(); ++itr )
    {
      const boost::filesystem::path &p = itr->path();
      const string pstr = p.string<string>();
      if( SpecUtils::is_directory( pstr ) )
      {
        auto subdirs = GadrasDirectory::recursive_list_gadras_drfs(pstr);
        if( subdirs.size() )
          pathstr += (pathstr.empty() ? "" : ";") + pstr;
      }
    }//for( loop over
#endif //if ANDROID / else
  }catch( std::exception &e )
  {
    cerr << "Got exception looking for GADRAS DRFs in user document dir: " << e.what() << endl;
  }//try / catch
#endif
  
  SpecUtils::split( paths, pathstr, "\r\n;" );
  
  //Make sure we always at least have the default generic detectors available.
  bool hasGeneric = false;
  for( size_t i = 0; !hasGeneric && (i < paths.size()); ++i )
    hasGeneric = (paths[i].find("GenericGadrasDetectors") != string::npos);
  
  if( !hasGeneric )
  {
    const string drfpaths = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GenericGadrasDetectors" );
    paths.push_back( drfpaths );
  }
  
    
  for( const string &path : paths )
  {
    auto dir = new GadrasDirectory( path, this, m_drfSelect, m_directories );
    
    if( path.find("GenericGadrasDetectors") != string::npos )
    {
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
      dir->m_directoryEdit->disable();
      dir->m_setDirectoryButton->hide();
      dir->m_deleteBtn->hide();
#endif
    }
    
#if( !defined(WIN32) )
    if( SpecUtils::istarts_with( path, "C:\\" ) )
      dir->hide();
#endif
  }

#if( IOS )
  const char *gadrasToolTip_key = "reds-gadras-drf-info-ios";
#else
  const char *gadrasToolTip_key = "reds-gadras-drf-info";
#endif
  
  WText *useInfo = new WText( WString::tr(gadrasToolTip_key), footer );
  useInfo->setStyleClass("DrfTypeDescrip");
}//void docreate()



MatchDetectorStatus GadrasDetSelect::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  docreate();  //We need to initialize the widget in order for m_directories to be true
  
  if( !m_directories )
    return MatchDetectorStatus::NoMatch;
  
  if( !det || det->name().empty() )
  {
    for( auto w : m_directories->children() )
    {
      GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
      if( !d )
        continue;  //shouldnt ever happen, but JIC
      d->m_detectorSelect->setCurrentIndex( -1 );
    }//for( auto w : m_directories->children() )
    
    return MatchDetectorStatus::NoMatch;
  }//if( !det || det->name().empty() )
  
  for( auto w : m_directories->children() )
  {
    GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
    if( !d )
      continue;  //shouldnt ever happen, but JIC
    
    for( size_t i = 0; i < d->m_responses.size(); ++i )
    {
      auto response = d->m_responses[i];
      
      if( (response == det) || (det->hashValue() == response->hashValue()) )
      {
        d->m_detectorSelect->setCurrentIndex( static_cast<int>(i) + 1 );
        return MatchDetectorStatus::Match;
      }
    }//for( size_t i = 0; i < d->m_responses.size(); ++i )
    
    d->m_detectorSelect->setCurrentIndex( 0 );
  }//for( auto w : m_directories->children() )
  
  return MatchDetectorStatus::NoMatch;
}//MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


std::shared_ptr<DetectorPeakResponse> GadrasDetSelect::selectedDetector()
{
  if( !m_directories )
    return nullptr;
  
  for( auto w : m_directories->children() )
  {
    GadrasDirectory *d = dynamic_cast<GadrasDirectory *>( w );
    const int currentIndex = d ? d->m_detectorSelect->currentIndex() : -1;
    if( currentIndex >= 1 )
    {
      try
      {
        WAbstractItemModel *m = d->m_detectorSelect->model();
        const string p = boost::any_cast<std::string>( m->data( currentIndex, 0, Wt::UserRole ) );
        auto answer = DrfSelect::initAGadrasDetectorFromDirectory( p );
        
        if( p.find("GenericGadrasDetectors") != string::npos )
          answer->setDrfSource( DetectorPeakResponse::DrfSource::DefaultGadrasDrf );
        else
          answer->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
        
        return answer;
      }catch( std::exception &e )
      {
        passMessage( WString::tr("reds-failed-parse-gadras-det").arg(e.what()),
                    WarningWidget::WarningMsgHigh );
      }//try / catch
    }//if( user has selected an index )
  }//for( auto w : m_directories->children() )
  
  return nullptr;
}//std::shared_ptr<DetectorPeakResponse> selectedDetector()


void GadrasDetSelect::addDirectory()
{
  if( m_directories )
    new GadrasDirectory( "", this, m_drfSelect, m_directories );
}//void addDirectory()


void GadrasDetSelect::removeDirectory( GadrasDirectory *dirWidget )
{
  if( !m_directories || !dirWidget )
    return;
  
  for( WWidget *w : m_directories->children() )
  {
    if( dynamic_cast<GadrasDirectory *>(w) == dirWidget )
    {
      delete w;
      break;
    }
  }//for( WWidget *w : m_files->children() )
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  saveFilePathToUserPreferences();
#endif
}//void removeDirectory( GadrasDirectory *dirWidget )


#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void GadrasDetSelect::saveFilePathToUserPreferences()
{
  if( !m_directories )
    return;
  
  string dirs;
  for( WWidget *w : m_directories->children() )
  {
    auto d = dynamic_cast<GadrasDirectory *>(w);
    if( !d )
      continue;
    const string p = d->m_directoryEdit->text().toUTF8();
    dirs += ((dirs.size() && !p.empty()) ? ";" : "") + p;
  }//for( WWidget *w : m_directories->children() )
  
  try
  {
    UserPreferences::setPreferenceValue( "GadrasDRFPath", dirs, m_interspec );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("reds-err-saving-gad-path").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//void saveFilePathToUserPreferences()
#endif


void GadrasDetSelect::detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det )
{
  if( !m_directories || !dir )
    return;
  
  for( WWidget *w : m_directories->children() )
  {
    auto child = dynamic_cast<GadrasDirectory *>(w);
    if( child && (child != dir) )
      child->selectNone();
  }//for( WWidget *w : m_files->children() )
}//void detectorSelected( GadrasDirectory *dir, std::shared_ptr<DetectorPeakResponse> det )


GadrasDirectory::GadrasDirectory( std::string directory, GadrasDetSelect *parentSelect,
                  DrfSelect *drfSelect, WContainerWidget *parent )
: WContainerWidget( parent ),
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  m_directory( directory ),
#else
  m_directoryEdit( new WLineEdit( Wt::WString::fromUTF8(directory) ) ),
  m_setDirectoryButton( new WPushButton( WString::tr("ds-set-btn") ) ),
#endif
  m_parent( parentSelect ),
  m_drfSelect( drfSelect ),
  m_detectorSelect( nullptr ),
  m_msg( nullptr ),
  m_deleteBtn( nullptr )
{
  setObjectName( "GadDir" + Wt::WRandom::generateId() );
  
  addStyleClass( "RelEffFile" );
  
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
#else
  m_directoryEdit->setAttributeValue( "ondragstart", "return false" );
  
  auto interspec = InterSpec::instance();
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec );
  
  
  WContainerWidget *topdiv = new WContainerWidget( this );
  m_deleteBtn = new WPushButton( topdiv );
  m_deleteBtn->addStyleClass( "closeicon-wtdefault" );
  //m_deleteBtn->setToolTip( WString::tr("reds-tt-remove-gad-dir") );
  HelpSystem::attachToolTipOn( m_deleteBtn, WString::tr("reds-tt-remove-gad-dir"),
                              showToolTips, HelpSystem::ToolTipPosition::Left );
  
  m_deleteBtn->clicked().connect( boost::bind( &GadrasDetSelect::removeDirectory, parentSelect, this ) );
  
  
//If we wanted to actually select directories, could do similar to file query widget... which webkitdirectory no longer seems to work, so see file query widget use of electrons dialog
//#if( BUILD_AS_ELECTRON_APP )
//  m_pathSelectedSignal.reset( new Wt::JSignal<std::string>( this, "BaseDirSelected", false ) );
//  const string uploadname = id() + "PathPicker";
//  const string uploadhtml = "<input id=\"" + uploadname + "\" type=\"file\" webkitdirectory=\"\" />";
  
//  WText *uploadtext = new WText( uploadhtml, XHTMLUnsafeText );
//  linelayout->addWidget( uploadtext, 0, 1 );
  
  //TODO: put in error handling!
//  wApp->doJavaScript( "document.getElementById('" + uploadname + "').onchange = function(event){"
//                     "var outputDir = document.getElementById('" + uploadname + "').files[0].path;"
//                     "Wt.emit( \"" + id() + "\", { name: 'BaseDirSelected' }, outputDir );"
//                     "};"
//                     );
  //m_pathSelectedSignal->connect( boost::bind( &SpecFileQueryWidget::newElectronPathSelected, this, boost::placeholders::_1 ) );
//#elif( BUILD_AS_OSX_APP )
//  SpecFileQuery::setIsSelectingDirectory( true );
//  setSearchDirectory( "" );
//  m_baseLocation = new WFileUpload();
//  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::newMacOsPathSelected );
//  linelayout->addWidget( m_baseLocation, 0, 1 );
//#else
  
  //TODO: Check if 'directory' is a subdirectory of InterSpec::writableDataDirectory()
  //and if so replace the beggingin of directory with {DataDir}/... and make it
  //  so the edit cant be changed and set button is hidden, and then deal with
  //  in GadrasDirectory::directory()
/*
  bool isInDataDir = false;
  string dirCononical = directory;
  string appDataDir = InterSpec::writableDataDirectory();
  if( SpecUtils::make_canonical_path(appDataDir)
     && SpecUtils::make_canonical_path(dirCononical)
     && appDataDir.size()>2 && dirCononical.size()>2 )
  {
    while( !isInDataDir && (dirCononical.size()+1) >= appDataDir.size() )
    {
      if( dirCononical == appDataDir )
      {
        isInDataDir = true;
        dirCononical = directory;
        SpecUtils::make_canonical_path(dirCononical);
        dirCononical = fs_relative( appDataDir, dirCononical );
        dirCononical = SpecUtils::append_path( "{AppDataDir}", dirCononical );
      }else
      {
        dirCononical = SpecUtils::parent_path(dirCononical);
        SpecUtils::make_canonical_path( dirCononical ); //to make sure we get ending / or \ right
      }
    }
  }
*/
  
  
/*
 ToDo: custonize path picking for electron...
#if( BUILD_AS_ELECTRON_APP )
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
  m_pathSelectedSignal->connect( boost::bind( &SpecFileQueryWidget::newElectronPathSelected, this, boost::placeholders::_1 ) );
#elif( BUILD_AS_OSX_APP )
 //For macOS dont saved picked directory to preferences in DB as sandboxing will mess this up.
  SpecFileQuery::setIsSelectingDirectory( true );
  setSearchDirectory( "" );
  m_baseLocation = new WFileUpload();
  m_baseLocation->changed().connect( this, &SpecFileQueryWidget::newMacOsPathSelected );
  linelayout->addWidget( m_baseLocation, 0, 1 );
#else
*/
  
  string user_data_dir;
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  try
  {
    user_data_dir = InterSpec::writableDataDirectory();
  }catch( std::exception & )
  {
  }
#endif
  
  if( !user_data_dir.empty() )
  {
    const string data_norm = SpecUtils::lexically_normalize_path( user_data_dir );
    const string file_norm = SpecUtils::lexically_normalize_path( directory );
    const string drf_norm = SpecUtils::append_path( data_norm, "drfs" );
    
    if( SpecUtils::istarts_with(file_norm, drf_norm) )
    {
      const string relpath = SpecUtils::fs_relative( drf_norm, directory );
      if( (relpath.size() >= 2) && !SpecUtils::istarts_with(relpath, "..") )
        directory = relpath;
    }else if( SpecUtils::istarts_with(file_norm, data_norm) )
    {
      const string relpath = SpecUtils::fs_relative( user_data_dir, directory );
      if( (relpath.size() >= 2) && !SpecUtils::istarts_with(relpath, "..") )
        directory = relpath;
    }
  }//if( !user_data_dir.empty() )
  
  
  topdiv->addWidget( m_directoryEdit );
  m_directoryEdit->setText( directory );
  m_directoryEdit->setTextSize( 48 );
  
#if( BUILD_AS_OSX_APP || IOS )
  m_directoryEdit->setAttributeValue( "autocorrect", "off" );
  m_directoryEdit->setAttributeValue( "spellcheck", "off" );
#endif

  
  topdiv->addWidget( m_setDirectoryButton );
  m_setDirectoryButton->setMargin( 5, Wt::Left );
  m_setDirectoryButton->disable();
  //m_setDirectoryButton->setToolTip( WString::tr("reds-tt-gad-set-dir") );
  HelpSystem::attachToolTipOn( m_setDirectoryButton, WString::tr("reds-tt-gad-set-dir"),
                              showToolTips, HelpSystem::ToolTipPosition::Left );
  
  m_setDirectoryButton->clicked().connect( this, &GadrasDirectory::dirPathChanged );
  
  //m_fileEdit->changed().connect( m_setFileButton, &WPushButton::enable );
  m_directoryEdit->textInput().connect( m_setDirectoryButton, &WPushButton::enable );
  m_directoryEdit->enterPressed().connect( this, &GadrasDirectory::dirPathChanged );
#endif
  
  WContainerWidget *bottomDiv = new WContainerWidget( this );
  
  m_detectorSelect = new WComboBox( bottomDiv );
  m_detectorSelect->setNoSelectionEnabled( true );
  m_detectorSelect->activated().connect( this, &GadrasDirectory::detectorSelectCallback );
  
#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  m_directoryEdit->changed().connect( m_detectorSelect, &WPushButton::disable );
#endif
  
  m_msg = new WText( bottomDiv );
  m_msg->setInline( false );
  
  initDetectors();
}//GadrasDirectory constructor


std::string GadrasDirectory::directory()
{
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) )
  return find_valid_path( m_directory, true );
#else
  return find_valid_path( m_directoryEdit->text().toUTF8(), true );
#endif
}


void GadrasDirectory::selectNone()
{
  if( m_detectorSelect->count() )
    m_detectorSelect->setCurrentIndex( 0 );
}


#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
void GadrasDirectory::dirPathChanged()
{
  m_parent->saveFilePathToUserPreferences();
  
  initDetectors();
  
  m_setDirectoryButton->disable();
  m_detectorSelect->enable();
}//void dirPathChanged()
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void GadrasDirectory::detectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det;
  
  const int index = m_detectorSelect->currentIndex();
  
  //Item at index 0 is always a non-detector.
  if( index > 0 && (index-1) < m_responses.size() )
    det = m_responses[index-1];
  
  if( m_parent )
    m_parent->detectorSelected( this, det );
  
  m_drfSelect->setDetector( det );
  m_drfSelect->emitChangedSignal();
}//void detectorSelectCallback()

  
MatchDetectorStatus GadrasDirectory::trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( !det )
  {
    m_detectorSelect->setCurrentIndex( 0 );
    return MatchDetectorStatus::NoMatch;
  }//if( !det )
  
  for( size_t i = 0; i < m_responses.size(); ++i )
  {
    if( m_responses[i]->hashValue() == det->hashValue() )
    {
      m_detectorSelect->setCurrentIndex( static_cast<int>(i) + 1 );
      return MatchDetectorStatus::Match;
    }
  }//for( size_t i = 0; i < m_responses.size(); ++i )
  
  m_detectorSelect->setCurrentIndex( 0 );
  
  if( m_detectorSelect->count() < 2 )
    return MatchDetectorStatus::NotInited;
  
  return MatchDetectorStatus::NoMatch;
}//MatchDetectorStatus trySelectDetector( std::shared_ptr<DetectorPeakResponse> det )


//parseDetector() returns null on error
std::shared_ptr<DetectorPeakResponse> GadrasDirectory::parseDetector( string path )
{
  try
  {
    const string name = SpecUtils::filename( path );
    auto drf = make_shared<DetectorPeakResponse>( name, "" );
    drf->fromGadrasDirectory( path );
    
//#if( PERFORM_DEVELOPER_CHECKS )
//    check_url_serialization( drf );
//#endif
    
    return drf;
  }catch( std::exception &e )
  {
    cerr << "GadrasDirectory::parseDetector() caught: " << e.what() << endl;
  }
 
  return nullptr;
}//shared_ptr<DetectorPeakResponse> parseDetector( const string directory )


vector<string> GadrasDirectory::recursive_list_gadras_drfs( const string &sourcedir )
{
  vector<string> files;
  if( !SpecUtils::is_directory( sourcedir ) )
    return files;
  
  const string csv_file = SpecUtils::append_path( sourcedir, "Efficiency.csv");
  const string dat_file = SpecUtils::append_path( sourcedir, "Detector.dat");
  
  if( SpecUtils::is_file(csv_file) && SpecUtils::is_file(dat_file) )
    files.push_back( sourcedir );
  
  //ToDo: Maybe we should use the Android implemenatation always.
#if( ANDROID )
  const vector<string> subdirs = SpecUtils::ls_directories_in_directory( sourcedir );
  for( const string &subdirpath : subdirs )
  {
    auto subsubdirs = recursive_list_gadras_drfs( subdirpath );
    files.insert( end(files), begin(subsubdirs), end(subsubdirs) );
  }//for( const string &subdir : subdirs )
#else
using namespace boost::filesystem;
  directory_iterator end_itr; // default construction yields past-the-end
  
  directory_iterator itr;
  try
  {
    itr = directory_iterator( sourcedir );
  }catch( std::exception & )
  {
    //ex: boost::filesystem::filesystem_error: boost::filesystem::directory_iterator::construct: Permission denied: "..."
    return files;
  }
  
  for( ; itr != end_itr; ++itr )
  {
    const boost::filesystem::path &p = itr->path();
    const string pstr = p.string<string>();
    if( SpecUtils::is_directory( pstr ) )
    {
      auto subdirs = recursive_list_gadras_drfs(pstr);
      files.insert( end(files), begin(subdirs), end(subdirs) );
    }
  }//for( loop over
#endif //if ANDROID / else
  
  return files;
}//vector<string> recursive_list_gadras_drfs( const string &sourcedir )


void GadrasDirectory::initDetectors()
{
  m_detectorSelect->clear();
  m_responses.clear();
  
  m_msg->setText( "" );
  m_msg->hide();
  
  string basedir = directory();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  basedir = complete_drf_dir_path( basedir );
#endif
  
  
  if( !SpecUtils::is_directory( basedir ) )
  {
    m_msg->setText( WString("<span style=\"color:red;\">{1}</span>").arg( WString::tr("reds-err-not-valid-dir") ) );
    m_msg->show();
    m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("reds-invalid-dir") ) );
    m_detectorSelect->setCurrentIndex( 0 );
    m_detectorSelect->disable();
    return;
  }//if( !SpecUtils::is_directory( basedir ) )
  
  m_detectorSelect->disable();
  this->disable();
  
  
  const std::string objname = objectName();
  const std::string sessid = wApp->sessionId();
  
  auto updategui = [this,objname]( std::vector<std::shared_ptr<DetectorPeakResponse> > drfs ){
    
    //Lets make sure this widget hasnt been deleted, by looking for it in the DOM
    WWidget *w = wApp->findWidget(objname);
    if( !w )
    {
      cerr << "Couldnt find widget '" << objname << "' in the DOM." << endl;
      return;
    }
    
    //cout << "Found widget '" << objname << "' in the DOM!" << endl;
    
    for( auto drf : drfs )
    {
      if( drf )
        m_responses.push_back( drf );
    }
    
    if( m_responses.empty() )
    {
      m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("reds-no-drf-in-dir") ) );
      m_detectorSelect->disable();
    }else
    {
      m_detectorSelect->addItem( WString("<{1}>").arg( WString::tr("ref-select-det-from-file") ) );
      m_detectorSelect->enable();
    }
    
    std::sort( begin(m_responses), end(m_responses),
               []( const shared_ptr<DetectorPeakResponse> &lhs,
                   const shared_ptr<DetectorPeakResponse> &rhs ) -> bool {
     if( !lhs || !rhs )
       return false;
      return lhs->name() < rhs->name();
    } );
    
    
    for( const auto drf : m_responses )
      m_detectorSelect->addItem( drf->name() );
    m_detectorSelect->setCurrentIndex( 0 );
    
    if( m_responses.empty() )
    {
      m_msg->setText( WString("<span style=\"color:red;\">{1}</span>").arg("reds-recursive-no-drfs-in-dir") );
      m_msg->show();
    }else
    {
      m_detectorSelect->enable();
    }
    
    this->enable();
    
    // Call DrfSelect::detectorsWereInited() since we are updating from a worker thread;
    //  the main thread will check, and fail if m_detector is a GADRAS detector,
    //  before we get here.
    m_drfSelect->detectorsWereInited();
    
    wApp->triggerUpdate();
  };//updategui lambda
  
  auto searchpaths = [basedir, objname, sessid, updategui](){
    const vector<string> dirs = recursive_list_gadras_drfs( basedir );
    std::vector<std::shared_ptr<DetectorPeakResponse> > dets( dirs.size(), nullptr );
    
    SpecUtilsAsync::ThreadPool pool;
    for( size_t i = 0; i < dirs.size(); ++i )
      pool.post( [i,&dets,&dirs](){ dets[i] = GadrasDirectory::parseDetector( dirs[i] ); } );
    pool.join();
    
    //for( const auto &d : dirs )
    //{
    //  auto drf = parseDetector( d );
    //  if( drf )
    //    m_responses.push_back( drf );
    //}
    //end section that should be done off the main thread
    
    Wt::WServer *server = Wt::WServer::instance();
    if( server )
      server->post(sessid, std::bind( [updategui,dets](){ updategui(dets); } ) );
  };//searchpaths lamda
  
  
  Wt::WServer *server = Wt::WServer::instance();
  Wt::WIOService &io = server->ioService();
  io.boost::asio::io_service::post( searchpaths );
}//void initDetectors()


void GadrasDirectory::detectorSelected( const int index )
{
  std::shared_ptr<DetectorPeakResponse> det;
  if( index > 0 )
    det = m_responses.at( index - 1 );
  
  m_drfSelect->setDetector( det );
  m_drfSelect->emitChangedSignal();
}//void detectorSelected( const int index )





DetectorDisplay::DetectorDisplay( InterSpec *specViewer,
                                  SpectraFileModel *fileModel,
                                  WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_text( NULL ),
    m_interspec( specViewer ),
    m_fileModel( fileModel )
{
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  HelpSystem::attachToolTipOn( this, WString::tr("app-tt-det-disp"), showToolTips ); //"app-tt-det-disp" defined in InterSpec.xml localisation, so we dont need to load DrfSelect.xml for this widget
  
  addStyleClass( "DetectorDisplay" );  //In InterSpec.css since this widget is loaded almost always at initial load time anyway

  new WImage( "InterSpec_resources/images/detector_small_white.png", this );
  new WText( WString("{1}:").arg(WString::tr("Detector") ), this );
  const bool isMobile = (m_interspec && m_interspec->isMobile());
  
  WString txt = WString("<font style=\"font-weight:100;color:#CFCFCF;\">&lt;{1}&gt;</font>")
                .arg( WString::tr(isMobile ? "app-det-select-txt-mobile" : "app-det-select-txt") );
  
  m_text = new WText( txt, XHTMLText, this );
  m_text->addStyleClass( "DetName" );

  std::shared_ptr<DetectorPeakResponse> detector;
  auto meas = specViewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( meas )
    detector = meas->detector();
  m_currentDetector = detector;

  if( detector && detector->isValid() )
    m_text->setText( detector->name() );

  specViewer->detectorChanged().connect( this, &DetectorDisplay::setDetector );
  specViewer->detectorModified().connect( this, &DetectorDisplay::setDetector );

  clicked().connect( this, &DetectorDisplay::editDetector );
}//DetectorDisplay constructor


DetectorDisplay::~DetectorDisplay()
{
} //~DetectorDisplay()


std::shared_ptr<DetectorPeakResponse> DetectorDisplay::detector()
{
  return m_currentDetector.lock();
}// std::shared_ptr<DetectorPeakResponse> DetectorDisplay::detector()


std::shared_ptr<const DetectorPeakResponse> DetectorDisplay::detector() const
{
  return m_currentDetector.lock();
}// std::shared_ptr<const DetectorPeakResponse> DetectorDisplay::detector() const


void DetectorDisplay::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  m_currentDetector = det;

  if( det && det->isValid() )
  {
    m_text->setText( det->name() );
  }else
  {
    const bool isMobile = (m_interspec && m_interspec->isMobile());
    WString txt = WString("<font style=\"font-weight:100;color:#CFCFCF;\">&lt;{1}&gt;</font>")
                  .arg( WString::tr(isMobile ? "app-det-select-txt-mobile" : "app-det-select-txt") );
    m_text->setText( txt );
  }
}//void setDetector( std::shared_ptr<DetectorPeakResponse> det )


void DetectorDisplay::editDetector()
{
  if( m_interspec )
    m_interspec->showDrfSelectWindow();
}//void editDetector()


std::shared_ptr<DetectorPeakResponse> DrfSelect::detector()
{
  return m_detector;
} //std::shared_ptr<DetectorPeakResponse> DrfSelect::detector()

DrfSelect::DrfSelect( std::shared_ptr<DetectorPeakResponse> currentDet,
                InterSpec *specViewer,
                SpectraFileModel *fileModel,
                AuxWindow* auxWindow )
  : WContainerWidget(),
    m_footer( nullptr ),
    m_interspec( specViewer ),
    m_fileModel( fileModel ),
    m_chart( 0 ),
    m_detector( currentDet ),
    m_gui_select_matches_det( false ),
    m_tabs( nullptr ),
    m_detectorDiameter( nullptr ),
    m_uploadedDetNameDiv( nullptr ),
    m_uploadedDetName( nullptr ),
    m_detectrDiameterDiv( nullptr ),
    m_efficiencyCsvUpload( nullptr ),
    m_detectrDotDatDiv( nullptr ),
    m_detectorDotDatUpload( nullptr ),
    m_fixedGeometryCb( nullptr ),
    m_acceptButton( nullptr ),
    m_cancelButton( nullptr ),
    m_noDrfButton( nullptr ),
    m_xmlDownload( nullptr ),
    m_detectorManualFunctionName( nullptr ),
    m_detectorManualFunctionText( nullptr ),
    m_detectorManualDescription( nullptr ),
    m_eqnEnergyGroup( nullptr ),
    m_drfType( nullptr ),
    m_detectorManualDiameterLabel( nullptr ),
    m_detectorManualDiameterText( nullptr ),
    m_detectorManualDistText( nullptr ),
    m_detectorManualDistLabel( nullptr ),
    m_detectorManualMinEnergy( nullptr ),
    m_detectorManualMaxEnergy( nullptr ),
    m_manualSetButton( nullptr ),
    m_gadrasDetSelect( nullptr ),
    m_relEffSelect( nullptr ),
    m_drfTypeMenu( nullptr ),
    m_drfTypeStack( nullptr ),
    m_deleteButton( nullptr ),
    m_DBtable( nullptr ),
    m_model( nullptr ),
    m_defaultForSerialNumber( nullptr ),
    m_defaultForDetectorModel( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/DrfSelect.css" );
  m_interspec->useMessageResourceBundle( "DrfSelect" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", specViewer );
  
  WGridLayout *mainLayout = new WGridLayout();
  
  setLayout( mainLayout );
  mainLayout->setContentsMargins(0, 0, 0, 0);
  mainLayout->setHorizontalSpacing( 0 );
  
#if( IOS )
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  if( app )
  {
    float safeAreas[4] = { 0.0f };
    InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
    app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
    if( safeAreas[3] > 1 )
      setPadding( WLength(safeAreas[3], WLength::Pixel), Side::Left );
  }//if( app )
#endif
  
  const bool is_mobile = (!specViewer || specViewer->isMobile());
  const bool is_phone = (!specViewer || specViewer->isMobile());
  int w = m_interspec->renderedWidth();
  int h = m_interspec->renderedHeight();
  if( is_mobile && (w < 100) )
  {
    w = wApp->environment().screenWidth();
    h = wApp->environment().screenHeight();
  }
  const bool narrow_layout = (is_phone && ((w > 100) && (w < 600)));
  if( narrow_layout )
    addStyleClass( "DrfSelectNarrow" );
  
  mainLayout->setVerticalSpacing( (is_mobile && !narrow_layout) ? 15 : 5 ); //so the chart resizer handle will show up
  
  specViewer->detectorChanged().connect( this, &DrfSelect::setDetector );
  specViewer->detectorModified().connect( this, &DrfSelect::setDetector );

  if( currentDet )
  {
    m_originalDetectorCopy.reset( new DetectorPeakResponse() );
    (*m_originalDetectorCopy) = (*currentDet);
    m_originalDetector = currentDet;
    m_previousDetectorDef = *currentDet;
  }//if( currentDet )

  m_previousEmmittedDetector = m_detector;
  
  m_chart = new DrfChart();
  
  if( is_phone )
  {
    const double chart_width = narrow_layout ? std::min( std::max( 100.0, w - 40.0 ), 450.0 ) : 325;
    const double chart_height = chart_width * 125.0/325.0;
    m_chart->resize( chart_width, chart_height );
    mainLayout->addWidget( m_chart, 0, 0, AlignCenter );
  }else
  {
    mainLayout->addWidget( m_chart, 0, 0 );
  }//if( we are on a phone ) / else
  
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  distValidator->setInvalidBlankText( "0.0 cm" );
  distValidator->setMandatory( true );
  
  WContainerWidget *lowerContent = new WContainerWidget();
  mainLayout->addWidget( lowerContent, 1, 0 );
  
  WGridLayout *lowerLayout = new WGridLayout();
  lowerContent->setLayout( lowerLayout );
  lowerLayout->setContentsMargins( 0, 2, 0, 0 );
  
  //Pre-size the lower content to accommodate the tallest DRF type ("Formula"),
  //  so everything wont change size when the user selects the different types.
  //  Although, as it stands now if you add a bunch of paths to GADRAS/Rel. Eff.
  //  you can make its content larger than the 190 px.
  //ToDo: improve the sizing of this layout to not have hardcoded sizes!
  lowerContent->resize( WLength::Auto, WLength(190,WLength::Pixel) );
  
  m_drfTypeStack = new Wt::WStackedWidget();
  m_drfTypeStack->addStyleClass( "UseInfoStack DetEditContent" );
  
  m_drfTypeMenu = new WMenu( m_drfTypeStack, Wt::Vertical );
  WContainerWidget *menuHolder = new WContainerWidget();
  menuHolder->addWidget( m_drfTypeMenu );
  
  if( narrow_layout )
  {
    m_drfTypeMenu->addStyleClass( "VerticalNavMenu HorizontalMenu HeavyNavMenu DetEditMenuHorizontal" );
    menuHolder->setOverflow( WContainerWidget::Overflow::OverflowAuto, Wt::Orientation::Horizontal );
    menuHolder->setOverflow( WContainerWidget::Overflow::OverflowHidden, Wt::Orientation::Vertical );
    
    lowerLayout->addWidget( menuHolder, 0, 0 );
    lowerLayout->addWidget( m_drfTypeStack, 1, 0 );
    lowerLayout->setRowStretch( 1, 1 );
  }else
  {
    m_drfTypeMenu->addStyleClass( "VerticalNavMenu SideMenu HeavyNavMenu DetEditMenu" );
    menuHolder->setOverflow( WContainerWidget::Overflow::OverflowAuto, Wt::Orientation::Vertical );
    menuHolder->setOverflow( WContainerWidget::Overflow::OverflowHidden, Wt::Orientation::Horizontal );
    
    lowerLayout->addWidget( menuHolder, 0, 0, 2, 1 );
    lowerLayout->addWidget( m_drfTypeStack, 0, 1 );
    lowerLayout->setColumnStretch( 1, 1 );
  }
  
  
  
  WContainerWidget *defaultOptions = new WContainerWidget();
  mainLayout->addWidget( defaultOptions, 2, 0 );
  //lowerLayout->addWidget( defaultOptions, 1, 1 );
  
  mainLayout->setRowResizable( 0, true, WLength(is_phone ? 125 : 250, WLength::Unit::Pixel) );
  
  mainLayout->setRowStretch( 1, 1 );
  
  
  //-------------------------------------
  //--- 1)  GADRAS
  //-------------------------------------
  
  m_gadrasDetSelect = new GadrasDetSelect( m_interspec, this );
  
  //-------------------------------------
  //--- 2)  Relative Efficiencies
  //-------------------------------------

  m_relEffSelect = new RelEffDetSelect( m_interspec, this );

  //-------------------------------------
  //--- 3)  Upload
  //-------------------------------------

  
  WContainerWidget *uploadDetTab = new WContainerWidget();
  
  WText *descrip = new WText( WString::tr("ds-csv-upload-desc"), uploadDetTab );
  descrip->setInline( false );
  descrip->setStyleClass("DetectorLabel");
  descrip->setInline( false );
  
  
  m_efficiencyCsvUpload = new WFileUpload( uploadDetTab );
  m_efficiencyCsvUpload->setInline( false );
  m_efficiencyCsvUpload->uploaded().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::EfficiencyCsvUploaded ) );
  m_efficiencyCsvUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge,
                                                             boost::placeholders::_1) );
  m_efficiencyCsvUpload->changed().connect( m_efficiencyCsvUpload, &WFileUpload::upload );
  m_efficiencyCsvUpload->setInline( false );
 
  m_detectrDiameterDiv = new WContainerWidget( uploadDetTab );
  m_detectrDiameterDiv->addStyleClass( "DetectorDiamDiv" );
  WLabel *label = new WLabel( WString::tr("ds-enter-gad-det-diam"), m_detectrDiameterDiv );
  label->setInline( false );

  label = new WLabel( WString::tr("ds-det-diam"), m_detectrDiameterDiv );
  m_detectorDiameter = new WLineEdit( "0 cm", m_detectrDiameterDiv );
  label->setBuddy( m_detectorDiameter );

  m_detectorDiameter->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detectorDiameter->setAttributeValue( "autocorrect", "off" );
  m_detectorDiameter->setAttributeValue( "spellcheck", "off" );
#endif
  
  m_detectorDiameter->setValidator( distValidator );
  m_detectorDiameter->setTextSize( 10 );
  m_detectorDiameter->changed().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );
  m_detectorDiameter->enterPressed().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );
  m_detectorDiameter->blurred().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::DetectorDiameterChanged ) );

  m_fixedGeometryCb = new WCheckBox( WString::tr("ds-fixed-geom"), uploadDetTab );
  m_fixedGeometryCb->setInline( false );
  m_fixedGeometryCb->addStyleClass( "FixedGeometry" );
  m_fixedGeometryCb->hide();
  m_fixedGeometryCb->checked().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::FixedGeometryChanged ) );
  m_fixedGeometryCb->unChecked().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::FixedGeometryChanged ) );
  
  
  m_uploadedDetNameDiv = new WContainerWidget( m_detectrDiameterDiv );
  label = new WLabel( WString("{1}:").arg( WString::tr("Name") ), m_uploadedDetNameDiv );
  m_uploadedDetName = new WLineEdit( m_uploadedDetNameDiv );
  label->setBuddy( m_uploadedDetName );
  
  m_uploadedDetName->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_uploadedDetName->setAttributeValue( "autocorrect", "off" );
  m_uploadedDetName->setAttributeValue( "spellcheck", "off" );
#endif
  
  m_uploadedDetName->textInput().connect( this, &DrfSelect::handleUserChangedUploadedDrfName );
  m_uploadedDetName->setTextSize( 30 );
  
  m_detectrDotDatDiv = new WContainerWidget( m_detectrDiameterDiv );
  m_detectrDotDatDiv->addStyleClass( "DetectorDotDatDiv" );

  label = new WLabel( WString::tr("ds-gad-det-file-label"), m_detectrDotDatDiv );
  label->setInline( false );
  m_detectorDotDatUpload = new WFileUpload( m_detectrDotDatDiv );
  m_detectorDotDatUpload->setInline( false );
  m_detectorDotDatUpload->uploaded().connect( boost::bind( &DrfSelect::fileUploadedCallback, this, UploadCallbackReason::DetectorDotDatUploaded ) );
  m_detectorDotDatUpload->fileTooLarge().connect( boost::bind(&SpecMeasManager::fileTooLarge,
                                                              boost::placeholders::_1) );
  m_detectorDotDatUpload->changed().connect( m_detectorDotDatUpload, &WFileUpload::upload );
  m_detectrDiameterDiv->hide();
  m_detectrDiameterDiv->setHiddenKeepsGeometry( true );

    
  //-------------------------------------
  //--- 4)  Manual
  //-------------------------------------

  const char *const diamtxt = "2.2 cm";
    
  WContainerWidget *formulaDiv = new WContainerWidget();
  WString nameLabel = WString::tr("ds-manual-det-desc");
  if( narrow_layout )
    nameLabel = WString("{1}:").arg(nameLabel);
  
  WText *selfDefineLabel = new WText( nameLabel, formulaDiv );
  selfDefineLabel->setInline( false );
  selfDefineLabel->setStyleClass("DetectorLabel");
  
  WTable *formulaTable = new WTable( formulaDiv );
  formulaTable->addStyleClass( "FormulaDrfTbl" );
  WTableCell *cell = formulaTable->elementAt( 0, 0 );
  label = new WLabel( WString::tr("ds-manual-det-name-label"), cell );
  if( narrow_layout )
  {
    cell = formulaTable->elementAt( 1, 0 );
    cell->setColumnSpan( 2 );
  }else
  {
    cell = formulaTable->elementAt( 0, 1 );
  }
  m_detectorManualFunctionName = new WLineEdit( cell );
  label->setBuddy( m_detectorManualFunctionName );

  m_detectorManualFunctionName->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detectorManualFunctionName->setAttributeValue( "autocorrect", "off" );
  m_detectorManualFunctionName->setAttributeValue( "spellcheck", "off" );
#endif
  
  m_detectorManualFunctionName->setStyleClass("DrfSelectFunctionalFormText");
  m_detectorManualFunctionName->setText( WString::tr("app-default-manual-det-name") );
  m_detectorManualFunctionName->setEmptyText( WString::tr("ds-manual-det-name-empty-txt") );
  m_detectorManualFunctionName->setInline(false);
  //m_detectorManualFunctionName->keyWentUp().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_detectorManualFunctionName->textInput().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));

  if( narrow_layout )
  {
    m_detectorManualFunctionName->setWidth( WLength(w - 65,WLength::Pixel) );
    m_detectorManualFunctionName->setMargin( WLength::Auto, Wt::Side::Left );
  }
  
  cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  label = new WLabel( WString::tr("ds-manual-function-label"), cell );
  if( narrow_layout )
  {
    cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
    cell->setColumnSpan( 2 );
  }else
  {
    cell = formulaTable->elementAt( cell->row(), 1 );
  }
  cell->setRowSpan( 2 );
  m_detectorManualFunctionText = new WTextArea( cell );
  //m_detectorManualFunctionText->setColumns(60);
  label->setBuddy( m_detectorManualFunctionText );
  
  m_detectorManualFunctionText->setRows(3);
  m_detectorManualFunctionText->setStyleClass("DrfSelectFunctionalFormText");
  //m_detectorManualFunctionText->setText( WString::tr("ds-manual-det-empty-fcn") );
  m_detectorManualFunctionText->setEmptyText( WString::tr("ds-manual-det-empty-fcn") );
  m_detectorManualFunctionText->setInline(false);
  //m_detectorManualFunctionText->keyWentUp().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_detectorManualFunctionText->textInput().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  if( narrow_layout )
  {
    m_detectorManualFunctionText->setWidth( WLength(w - 65,WLength::Pixel) );
    m_detectorManualFunctionText->setMargin( WLength::Auto, Wt::Side::Left );
  }
  
  if( narrow_layout )
    cell = formulaTable->elementAt( cell->row() + 2, 0 );
  else
    cell = formulaTable->elementAt( cell->row() + 1, 0 );
  
  Wt::WContainerWidget *energyContainer = new Wt::WContainerWidget( cell );
  energyContainer->setMargin( 6, Wt::Top );
  m_eqnEnergyGroup = new WButtonGroup( energyContainer );
  WRadioButton *button = new WRadioButton( "keV", energyContainer );
  m_eqnEnergyGroup->addButton( button, 0 );
  button = new WRadioButton( "MeV", energyContainer );
  button->setMargin( 5, Wt::Left );
  m_eqnEnergyGroup->addButton( button, 1 );
  m_eqnEnergyGroup->checkedChanged().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_eqnEnergyGroup->setSelectedButtonIndex( 0 );
    
  HelpSystem::attachToolTipOn( energyContainer, WString::tr("ds-tt-manual-energy-unit"),
                              showToolTips, HelpSystem::ToolTipPosition::Top );

  cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  
  WString descLabel = WString::tr("Description");
  if( narrow_layout )
    descLabel = WString("{1}:").arg(descLabel);
  
  label = new WLabel( descLabel, cell );
  if( narrow_layout )
  {
    cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
    cell->setColumnSpan( 2 );
  }else
  {
    cell = formulaTable->elementAt( cell->row(), 1 );
  }
  m_detectorManualDescription = new WLineEdit( cell );
  label->setBuddy( m_detectorManualDescription );
  if( narrow_layout )
  {
    m_detectorManualDescription->setWidth( WLength(w - 65,WLength::Pixel) );
    m_detectorManualDescription->setMargin( WLength::Auto, Wt::Side::Left );
    m_detectorManualDescription->setInline( false );
  }else
  {
    m_detectorManualDescription->setWidth( WLength(100,WLength::Percentage) );
  }
  m_detectorManualDescription->setEmptyText( WString::tr("ds-tt-manual-desc") );
  m_detectorManualDescription->changed().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  m_detectorManualDescription->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detectorManualDescription->setAttributeValue( "autocorrect", "off" );
  m_detectorManualDescription->setAttributeValue( "spellcheck", "off" );
#endif
  
  cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  m_detectorManualDiameterLabel = new WLabel( WString::tr("ds-manual-det-diam-label"), cell );
  cell = formulaTable->elementAt( cell->row(), 1 );
  m_detectorManualDiameterText = new WLineEdit( cell );
  m_detectorManualDiameterLabel->setBuddy( m_detectorManualDiameterText );
  m_detectorManualDiameterText->setValidator( distValidator );
  //m_detectorManualDiameterText->setText( diamtxt );
  m_detectorManualDiameterText->setEmptyText( diamtxt );
  m_detectorManualDiameterText->setValidator( distValidator );
  m_detectorManualDiameterText->blurred().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_detectorManualDiameterText->enterPressed().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  //m_detectorManualDiameterText->keyWentUp().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_detectorManualDiameterText->textInput().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  
  m_detectorManualDiameterText->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_detectorManualDiameterText->setAttributeValue( "autocorrect", "off" );
  m_detectorManualDiameterText->setAttributeValue( "spellcheck", "off" );
#endif
  
  
  cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  label = new WLabel( WString::tr("ds-energy-range"), cell );
  
  cell = formulaTable->elementAt( cell->row(), 1 );
  label = new WLabel( WString::tr("Min."), cell );
  m_detectorManualMinEnergy = new NativeFloatSpinBox( cell );
  m_detectorManualMinEnergy->setPlaceholderText( WString::tr("ds-(optional)") );
  m_detectorManualMinEnergy->setRange( 0, 15000 );
  m_detectorManualMinEnergy->setSpinnerHidden( true );
  m_detectorManualMinEnergy->setWidth( WLength(9,Wt::WLength::FontEx) );
  m_detectorManualMinEnergy->setText( "" );
  label->setBuddy( m_detectorManualMinEnergy );
  m_detectorManualMinEnergy->valueChanged().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  if( narrow_layout )
    cell = formulaTable->elementAt( formulaTable->rowCount(), 1 );
    
  label = new WLabel( WString::tr("Max."), cell );
  if( !narrow_layout )
    label->setMargin( 20, Wt::Side::Left );
  m_detectorManualMaxEnergy = new NativeFloatSpinBox( cell );
  m_detectorManualMaxEnergy->setPlaceholderText( WString::tr("ds-(optional)") );
  m_detectorManualMaxEnergy->setRange( 0, 15000 );
  m_detectorManualMaxEnergy->setSpinnerHidden( true );
  m_detectorManualMaxEnergy->setWidth( WLength(9,Wt::WLength::FontEx) );
  m_detectorManualMaxEnergy->setText( "" );
  label->setBuddy( m_detectorManualMaxEnergy );
  if( !narrow_layout )
  {
    label = new WLabel("(keV)", cell );
    label->setMargin( 10, Wt::Side::Left );
  }
  
  m_detectorManualMaxEnergy->valueChanged().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  
  cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  cell->setColumnSpan( 2 );
  m_drfType = new WComboBox( cell );
  m_drfType->addItem( WString::tr("ds-intrinsic-eff") );
  m_drfType->addItem( WString::tr("ds-abs-eff") );
  m_drfType->addItem( WString::tr("ds-fixed-geom-total") );
  m_drfType->addItem( WString::tr("ds-fixed-geom-cm2") );
  m_drfType->addItem( WString::tr("ds-fixed-geom-m2") );
  m_drfType->addItem( WString::tr("ds-fixed-geom-gram") );
  //m_drfType->setToolTip( WString::tr("ds-tt-manual-drf-type") );
  HelpSystem::attachToolTipOn( m_drfType, WString::tr("ds-tt-manual-drf-type"), showToolTips );
  
  
  m_drfType->setCurrentIndex( 0 );
  m_drfType->activated().connect( this, &DrfSelect::verifyManualDefinition );
  
  if( narrow_layout )
    cell = formulaTable->elementAt( formulaTable->rowCount(), 0 );
  m_detectorManualDistLabel = new WLabel( WString::tr("ds-dist-label"), cell );
  m_detectorManualDistLabel->setMargin( 10, Wt::Left );
  
  if( narrow_layout )
    cell = formulaTable->elementAt( formulaTable->rowCount(), 1 );
  m_detectorManualDistText = new Wt::WLineEdit( cell );
  m_detectorManualDistText->setValidator( distValidator );
  m_detectorManualDistText->setEmptyText("1 m");
  if( !narrow_layout )
    m_detectorManualDistText->setHiddenKeepsGeometry( true );
  m_detectorManualDistLabel->hide();
  m_detectorManualDistText->hide();
  //m_detectorManualDistText->keyWentUp().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  m_detectorManualDistText->textInput().connect(boost::bind(&DrfSelect::verifyManualDefinition, this));
  
  
  if( narrow_layout )
    cell = formulaTable->elementAt( formulaTable->rowCount(), 1 );
  m_manualSetButton = new WPushButton( WString::tr("ds-set-btn"), cell );
  //m_manualSetButton->setInline(false);
  m_manualSetButton->setFloatSide( Wt::Right );
  //selfDefineLayout->setColumnStretch( 1, 1 );
  m_manualSetButton->clicked().connect(boost::bind(&DrfSelect::setFormulaDefineDetector, this));
    
  //-------------------------------------
  //--- 5)  Recent
  //-------------------------------------
  WContainerWidget *recentDiv = new WContainerWidget( );
  WGridLayout* recentDivLayout = new WGridLayout( recentDiv );
  recentDivLayout->setContentsMargins( 0, 0, 0, 0 );

  const Dbo::ptr<InterSpecUser> &user = m_interspec->user();
  
  m_DBtable = new RowStretchTreeView();
  m_DBtable->setRootIsDecorated	(	false ); //makes the tree look like a table! :)
  
  m_DBtable->addStyleClass( "DbSpecFileSelectTable" );
  
  //We have to create a independant Dbo::Session for this class since the
  //  m_viewer->m_user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_sql.reset( new DataBaseUtils::DbSession( *m_interspec->sql() ) );

  m_model = new Dbo::QueryModel< Dbo::ptr<DetectorPeakResponse> >( m_DBtable );

  const auto userid = user.id();
  m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                           .where( "InterSpecUser_id = ? ").bind( userid )
                           .orderBy( "-1*m_lastUsedUtc" )
                    );
  
  /*  ToDo: - add column for when was last used, but give time in how long ago (remove diameter column)
            - Add bullets or dropdown to allow filtering by "All", "Uploaded", "Entered Formula", "Made From Data"
   */
  
  m_model->addColumn( "m_name" );
  m_model->addColumn( "m_description" );
  m_model->addColumn( "m_lastUsedUtc" );
  
  
  m_model->setHeaderData(  0, Horizontal, WString::tr("Name"), DisplayRole );
  m_DBtable->setColumnWidth( 0, 200 );
  m_model->setHeaderData(  1, Horizontal, WString::tr("Description"), DisplayRole );
  m_DBtable->setColumnWidth( 1, 150 );
  m_model->setHeaderData(  2, Horizontal, WString::tr("ds-last-used"), DisplayRole );
  m_DBtable->setColumnWidth( 2, 150 );
 
  
  WItemDelegate *delegate = new UtcToLocalTimeDelegate( m_DBtable );
  m_DBtable->setItemDelegateForColumn( 2, delegate );

  m_DBtable->setModel( m_model );
  m_DBtable->setAlternatingRowColors( true );
  m_DBtable->setSelectionMode( SingleSelection );
  
  for( int col = 0; col < m_model->columnCount(); ++col )
  {
    m_DBtable->setColumnHidden( col, false );
    m_DBtable->setSortingEnabled( col, true );
  }//for( int col = 0; col < model->columnCount(); ++col )
  
  m_DBtable->selectionChanged().connect( this, &DrfSelect::dbTableSelectionChanged );

  m_deleteButton = new WPushButton( WString::tr("Delete") );
  m_deleteButton->setIcon( "InterSpec_resources/images/minus_min_white.png" );
  m_deleteButton->setDefault(true);
  m_deleteButton->disable();
  m_deleteButton->clicked().connect( this, &DrfSelect::deleteDBTableSelected );

  WComboBox *filter = new WComboBox();
  filter->addItem( WString::tr("ds-det-type-filter-all") );
  filter->addItem( WString::tr("ds-det-type-filter-uploaded") );
  filter->addItem( WString::tr("ds-det-type-filter-formula") );
  filter->addItem( WString::tr("ds-det-type-filter-data") );
  filter->addItem( WString::tr("ds-det-type-filter-N42") );
  filter->setCurrentIndex( 0 );
  
  filter->changed().connect( std::bind( [filter,this,userid](){
    m_DBtable->setSelectedIndexes( WModelIndexSet() );
    m_deleteButton->disable();
    
    if( filter->currentIndex() == 0 ) //"All"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? ").bind( userid ), true );
    }else if( filter->currentIndex() == 1 ) //"Uploaded"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND (m_efficiencySource = ? OR m_efficiencySource = ?)")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserImportedIntrisicEfficiencyDrf )
                        .bind( DetectorPeakResponse::DrfSource::UserImportedGadrasDrf ) , true );
    }else if( filter->currentIndex() == 2 ) //"Entered Formula"
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf ) , true );
    }else if( filter->currentIndex() == 3 )
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::UserCreatedDrf ) , true );
    }else if( filter->currentIndex() == 4 )
    {
      m_model->setQuery( m_sql->session()->find< DetectorPeakResponse >()
                        .where( "InterSpecUser_id = ? AND m_efficiencySource = ?")
                        .bind( userid )
                        .bind( DetectorPeakResponse::DrfSource::FromSpectrumFileDrf ) , true );
    }else
    {
      return; //shouldnt happen
    }
    
    m_DBtable->sortByColumn( 2, Wt::SortOrder::DescendingOrder );
  } ) );
  
  recentDivLayout->addWidget( m_DBtable, 0, 0, 1, 3 );
  recentDivLayout->addWidget( m_deleteButton, 1, 0, AlignLeft);
  
  label = new WLabel( WString::tr("ds-show-label") );
  label->setBuddy( filter );
  recentDivLayout->addWidget( label, 1, 1, AlignRight | AlignMiddle );
  
  recentDivLayout->addWidget( filter, 1, 2, AlignMiddle );
  recentDivLayout->setRowStretch( 0, 1 );
  //recentDiv->setOverflow(Wt::WContainerWidget::OverflowHidden);
  //recentDiv->setMaximumSize( WLength::Auto, 180 );
  recentDiv->setHeight( WLength( 100, WLength::Unit::Percentage ) );
  m_DBtable->sortByColumn( 2, Wt::SortOrder::DescendingOrder );
  
  //--------------------------------------------------------------------------------
  
  WMenuItem *item = m_drfTypeMenu->addItem( WString::tr("ds-mi-gadras"), m_gadrasDetSelect );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( WString::tr("ds-mi-rel-eff"), m_relEffSelect );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( WString::tr("ds-mi-import"), uploadDetTab );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( WString::tr("ds-mi-formula"), formulaDiv );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  item = m_drfTypeMenu->addItem( WString::tr("ds-mi-previous"), recentDiv );
  item->clicked().connect( boost::bind(&right_select_item, m_drfTypeMenu, item) );
  
  //m_drfTypeStack->setHeight( WLength(185.0) );
  
  
  auto meas = specViewer->measurment( SpecUtils::SpectrumType::Foreground );
  
  if( meas && !meas->instrument_id().empty() )
  {
    // Hack to keep checkboxes visible
    if( !defaultOptions->hasStyleClass( "DrfDefaultOptions" ) )
      defaultOptions->addStyleClass( "DrfDefaultOptions" );
    
    m_defaultForSerialNumber = new WCheckBox( WString::tr("ds-use-for-serialnum-cb").arg( meas->instrument_id() ), defaultOptions );
    m_defaultForSerialNumber->setInline( false );
    //mainLayout->addWidget( m_defaultForSerialNumber, mainLayout->rowCount(), 0, 1, mainLayout->columnCount() );
  }//if( have serial number )
  
  if( meas && (meas->detector_type()!=DetectorType::Unknown || !meas->instrument_model().empty()) )
  {
    string model;
    if( meas->detector_type() == DetectorType::Unknown )
      model = meas->instrument_model();
    else
      model = detectorTypeToString( meas->detector_type() );
    
    // Hack to keep checkboxes visible
    if( !defaultOptions->hasStyleClass( "DrfDefaultOptions" ) )
      defaultOptions->addStyleClass( "DrfDefaultOptions" );
    
    m_defaultForDetectorModel = new WCheckBox( WString::tr("ds-use-for-model-cb").arg(model), defaultOptions );
    m_defaultForDetectorModel->setInline( false );
    //mainLayout->addWidget( m_defaultForDetectorModel, mainLayout->rowCount(), 0, 1, mainLayout->columnCount() );
  }//if( have
  
  WContainerWidget *secondFooter = nullptr;
  if( auxWindow )
  {
    m_footer = auxWindow->footer();
    if( narrow_layout )
    {
      secondFooter = new WContainerWidget();
      mainLayout->addWidget( secondFooter, mainLayout->rowCount(), 0, 1, 2 );
      secondFooter->addStyleClass( "DetEditNarrowSecondFooter" );
    }else
    {
      secondFooter = m_footer;
    }
  }else
  {
    m_footer = new WContainerWidget();
    secondFooter = m_footer;
    mainLayout->addWidget( m_footer, mainLayout->rowCount(), 0, 1, 2 );
  }
  
  AuxWindow::addHelpInFooter( m_footer, "detector-edit-dialog" );
  
  // Add cancel button first on phones, so it will stay all the way to the left
  if( auxWindow && auxWindow->isPhone() )
    m_cancelButton = auxWindow->addCloseButtonToFooter( WString::tr("Cancel") );
  
  DrfDownloadResource *xmlResource = new DrfDownloadResource( this );
#if( BUILD_AS_OSX_APP || IOS )
  m_xmlDownload = new WAnchor( WLink(xmlResource), secondFooter );
  m_xmlDownload->setTarget( AnchorTarget::TargetNewWindow );
  m_xmlDownload->setStyleClass( "LinkBtn DownloadLink DrfXmlDownload" );
#else
  m_xmlDownload = new WPushButton( secondFooter );
  m_xmlDownload->setIcon( "InterSpec_resources/images/download_small.svg" );
  m_xmlDownload->setLink( WLink(xmlResource) );
  m_xmlDownload->setLinkTarget( Wt::TargetNewWindow );
  m_xmlDownload->setStyleClass( "LinkBtn DownloadBtn DrfXmlDownload" );
    
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  m_xmlDownload->clicked().connect( std::bind([xmlResource](){
    android_download_workaround(xmlResource, "drf.xml");
  }) );
#endif //ANDROID

#endif
  
  m_xmlDownload->setText( "XML" );
  m_xmlDownload->setHidden( !m_detector || !m_detector->isValid() );
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( secondFooter );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DrfXmlDownload" );
  
  qr_btn->clicked().connect( std::bind( [this](){
  
    if( !m_detector || !m_detector->isValid() )
    {
      passMessage( WString::tr("ds-no-drf-loaded"), WarningWidget::WarningMsgHigh );
      return;
    }
  
    try
    {
      const string url = "interspec://drf/specify?" + Wt::Utils::urlEncode(m_detector->toAppUrl());
      QrCode::displayTxtAsQrCode( url, WString::fromUTF8(m_detector->name()),
                                 WString::fromUTF8(m_detector->description()) );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES
  
  m_noDrfButton = new WPushButton( WString::tr("ds-no-det-btn"), secondFooter );
  m_noDrfButton->clicked().connect( this, &DrfSelect::finishWithNoDetector );
  if( specViewer && !specViewer->isPhone() )
    m_noDrfButton->addStyleClass( "NoDrfBtn" );
  
  WPushButton *changeFwhm = new WPushButton( WString::tr("ds-fit-fwhm-btn"), secondFooter );
  if( specViewer && !specViewer->isPhone() )
    changeFwhm->addStyleClass( "NoDrfBtn" );
  changeFwhm->clicked().connect( this, &DrfSelect::handleFitFwhmRequested );
  
  if( auxWindow && !auxWindow->isPhone() )
  {
    m_cancelButton = auxWindow->addCloseButtonToFooter( WString::tr("Cancel") );
  }else if( !auxWindow )
  {
    m_cancelButton = new WPushButton( WString::tr("Close") );
    m_footer->insertWidget( m_footer->count(), m_cancelButton );
  }
  
  m_cancelButton->clicked().connect( this, &DrfSelect::cancelAndFinish );

  
  if( specViewer && !specViewer->isPhone() )
  {
    m_acceptButton = new WPushButton( WString::tr("Accept"), m_footer );
    m_acceptButton->setFloatSide(Wt::Right);
    
    m_cancelButton->setIcon( "InterSpec_resources/images/reject.png" );
    HelpSystem::attachToolTipOn( m_cancelButton, WString::tr("ds-tt-cancel"), showToolTips );
    
    m_acceptButton->setIcon( "InterSpec_resources/images/accept.png" );
    HelpSystem::attachToolTipOn( m_acceptButton, WString::tr("ds-tt-accept"), showToolTips );
  }else
  {
    m_acceptButton = new WPushButton( WString::tr("ds-use-det-btn"), m_footer );
    m_acceptButton->addStyleClass( "CenterBtnInMblAuxWindowHeader" );
  }//if( isMobile() ) / else
  
  m_acceptButton->clicked().connect( this, &DrfSelect::acceptAndFinish );
  if( !currentDet || !currentDet->isValid() )
    setAcceptButtonEnabled( false );
  
  m_drfTypeMenu->select( 0 );
  
  setGuiToCurrentDetector();
}//DrfSelect constructor


void DrfSelect::createChooseDrfDialog( vector<shared_ptr<DetectorPeakResponse>> inputdrfs,
                                       WString mainMsgHtml,
                                       string creditsHtml,
                                       std::function<void()> saveDrfsCallBack )
{
  wApp->useStyleSheet( "InterSpec_resources/DrfSelect.css" );
  
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  for( auto d : inputdrfs )
  {
    if( d && d->isValid() )
      drfs.push_back( d );
  }
  
  if( drfs.empty() )
    return;
  
  auto interspec = InterSpec::instance();
  auto meas = interspec ? interspec->measurment(SpecUtils::SpectrumType::Foreground) : nullptr;
  if( interspec )
    interspec->useMessageResourceBundle( "DrfSelect" );
  
  auto selected_drf = make_shared<shared_ptr<DetectorPeakResponse>>();
  
  const char *title_key = (drfs.size() == 1) ? "ds-use-drf" : "ds-select-drf";
  
  SimpleDialog *dialog = new SimpleDialog( WString::tr(title_key) );
  dialog->addStyleClass( "DrfFileSelectDialog" );
  
  WPushButton *cancel = dialog->addButton( WString::tr("Cancel") );
  WPushButton *accept = dialog->addButton( WString::tr("Accept") );
  
  //WGridLayout *layout = new WGridLayout( dialog->contents() );
  //layout->setContentsMargins( 0, 0, 0, 0 );
  
  dialog->contents()->setOverflow( Overflow::OverflowHidden, Wt::Orientation::Horizontal );
  
  DrfChart *chart = new DrfChart();
  chart->addStyleClass( "DrfFileSelectChart" );
  
  //layout->addWidget( chart, 0, 0 );
  dialog->contents()->addWidget( chart );
  
  
  if( !mainMsgHtml.empty() )
  {
    WText *gen_desc = new WText( mainMsgHtml );
    gen_desc->addStyleClass( "DrfFileSelectMainDesc" );
    
    //layout->addWidget( gen_desc, layout->rowCount(), 0 );
    gen_desc->setInline( false );
    dialog->contents()->addWidget( gen_desc );
  }//if( !descrip_html.empty() )
  
  
  // Create drfDescription so we can capture in lambda; we'll add to layout in a bit
  WText *drfDescription = new WText();
  drfDescription->addStyleClass( "DrfFileSelectDrfDesc" );

  auto showDrf = [drfs,chart,accept,drfDescription,selected_drf](const int index){
    if( (index < 0) || (index >= drfs.size()) )
    {
      // Shouldnt happen, but JIC
      accept->disable();
      drfDescription->setText( "" );
      chart->updateChart( nullptr );
      (*selected_drf) = nullptr;
      return;
    }
    
    auto drf = drfs[index];
    (*selected_drf) = drf;
    chart->updateChart( drf );
    
    WString drfdesc;
    if( !drf->description().empty() )
      drfdesc = WString("<b>{1}:</b>&nbsp;{2}").arg( WString::tr("ds-drf-desc") ).arg( drf->description() );
    
    drfDescription->setHidden( drfdesc.empty() );
    drfDescription->setText( drfdesc );
  };//showDrf lambda
  
  auto bestGuessDrf = []( std::shared_ptr<const SpecMeas> meas, const vector<shared_ptr<DetectorPeakResponse>> &drfs ) -> pair<shared_ptr<DetectorPeakResponse>,int> {
    if( !meas
       || ((meas->detector_type() == DetectorType::Unknown) && meas->instrument_model().empty()) )
      return {nullptr,-1};
    
    string model;
    if( meas->detector_type() == DetectorType::Unknown )
      model = meas->instrument_model();
    else
      model = detectorTypeToString( meas->detector_type() );
  
    // TODO: this bestGuessDrf function is not very good...
    //       instead of removing non-alhpanumeric, should match first word, or something...
    
    // Remove all non-alphanumeric
    auto removeNonAlpha = []( string &str ){
      str.erase( std::remove_if(begin(str), end(str), [](char a){return !isalnum(a);}), end(str) );
    };
    
    removeNonAlpha( model );
    if( model.size() < 2 )
      return {nullptr,-1};
    
    for( int i = 0; i < static_cast<int>(drfs.size()); ++i )
    {
      const shared_ptr<DetectorPeakResponse> &d = drfs[i];
      
      string name = d->name();
      removeNonAlpha( name );
      
      const size_t pos = name.find( model );
      if( pos != string::npos )
        return {d,i};
    }
    
    return {nullptr,-1};
  };//bestGuessDrf lamda
  
  int initialDrf = 0;
  
  if( drfs.size() > 1 )
  {
    // TODO: try to guess which DRF is closest, and select that one
    initialDrf = 1;
    
    const auto bestDrf = bestGuessDrf( meas, drfs );
    if( bestDrf.first && (bestDrf.second >= 0) )
      initialDrf = bestDrf.second;
  }//if( drfs.size() > 1 )
  
  showDrf( initialDrf );
  
  if( drfs.size() > 1 )
  {
    WComboBox *combo = new WComboBox();
    combo->addStyleClass( "DrfFileSelectCombo" );
    
    //layout->addWidget( combo, layout->rowCount(), 0, AlignCenter );
    combo->setInline( false );
    dialog->contents()->addWidget( combo );
    
    for( auto drf : drfs )
      combo->addItem( drf->name() );
    combo->setCurrentIndex( 0 );
    combo->activated().connect( std::bind(showDrf, std::placeholders::_1) );
  }//if( drfs.size() > 1 )
  
  
  //layout->addWidget( drfDescription, layout->rowCount(), 0 );
  drfDescription->setInline( false );
  dialog->contents()->addWidget( drfDescription );


  if( !creditsHtml.empty() )
  {
    WText *credits = new WText( creditsHtml );
    credits->addStyleClass( "DrfFileSelectCredits" );
    
    //layout->addWidget( credits, layout->rowCount(), 0 );
    dialog->contents()->addWidget( credits );
  }//if( !descrip_html.empty() )

  const bool makeDrfsSaveCb = !!saveDrfsCallBack;
  const bool makeSerialNumCb = (meas && !meas->instrument_id().empty());
  const bool makeModelCb = (meas
                            && ((meas->detector_type()!=DetectorType::Unknown)
                                 || !meas->instrument_model().empty()));
  
  WContainerWidget *saveCbDiv = nullptr;
  WCheckBox *saveCb = nullptr;
  WCheckBox *defaultForSerialNumber = nullptr, *defaultForDetectorModel = nullptr;
  
  if( makeDrfsSaveCb || makeSerialNumCb || makeModelCb )
  {
    saveCbDiv = new WContainerWidget( dialog->contents() );
    saveCbDiv->addStyleClass( "DrfFileSelectSaveCbs" );
  }
  
  if( makeDrfsSaveCb )
  {
    saveCb = new WCheckBox( WString::tr("ds-save-for-later-cb"), saveCbDiv );
    //layout->addWidget( saveCb, layout->rowCount(), 0, AlignLeft );
    saveCb->setInline( false );
    dialog->contents()->addWidget( saveCb );
  }//if( makeDrfsSaveCb )
  
  if( makeSerialNumCb )
  {
    WString msg = WString::tr("ds-use-for-serialnum-cb").arg( meas->instrument_id() );
    defaultForSerialNumber = new WCheckBox( msg, saveCbDiv );
    defaultForSerialNumber->setInline( false );
  }//if( have serial number )
  
  if( makeModelCb )
  {
    string model;
    if( meas->detector_type() == DetectorType::Unknown )
      model = meas->instrument_model();
    else
      model = detectorTypeToString( meas->detector_type() );
    
    
    WString msg = WString::tr("ds-use-for-model-cb").arg( model );
    defaultForDetectorModel = new WCheckBox( msg, saveCbDiv );
    defaultForDetectorModel->setInline( false );
  }//if( makeModelCb )
  
  
  accept->clicked().connect( std::bind( [=](){
    auto drf = *selected_drf;
    auto interspec = InterSpec::instance();
    if( !drf || !interspec )
      return;
      
    auto sql = interspec->sql();
    const Dbo::ptr<InterSpecUser> &user = interspec->user();
    DrfSelect::updateLastUsedTimeOrAddToDb( drf, user.id(), sql );
    interspec->detectorChanged().emit( drf ); //This loads it to the foreground spectrum file
      
    if( saveCb && saveCb->isChecked() && !!saveDrfsCallBack )
      saveDrfsCallBack();
    
    
    auto meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    
    
    if( meas && drf && defaultForSerialNumber && defaultForSerialNumber->isChecked() )
    {
      UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorSerialNumber;
      WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
        DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, meas );
      } ) );
    }//if( m_defaultForSerialNumber and is checked )
    
    if( meas && drf && defaultForDetectorModel && defaultForDetectorModel->isChecked() )
    {
      UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorModelName;
      WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
        DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, meas );
      } ) );
    }//if( m_defaultForDetectorModel and is checked )
  }) );
  
  if( interspec && (interspec->renderedWidth() > 100) && (interspec->renderedHeight() > 50) )
  {
    double dialogWidth = std::min(800.0, 0.5*interspec->renderedWidth());
    
    double chartHeight = 0.5*interspec->renderedHeight();
    chartHeight = std::min( std::min( chartHeight, 325.0 ), dialogWidth );
    
    chart->setWidth( dialogWidth - 25 );
    chart->setHeight( chartHeight );
    
    dialog->setWidth( dialogWidth );
    dialog->setMaximumSize( 0.5*interspec->renderedWidth(), 0.95*interspec->renderedHeight() );
  }else
  {
    // TODO: have some kind of a callback to resize the dialog, and/or call into wApp->environment().screenWidth() and wApp->environment().screenHeight(), for phones/tablets (untested if this works).
    
    // We get here like on iOS if opening a deep-link url is what opened the application.
    // A iPhone 13 has about 812 x 375 pixels root DOM (including past the notch), and max-width
    //  of a SimpleDialog is 50% width, so we'll resize the chart to a kinda small, but good-enough
    //  area of 350 x 200, which makes the SimpleDialog about 376 x 323 px
    chart->resize( 350, 200 );
    //dialog->setWidth( 640 );
  }
  
  dialog->rejectWhenEscapePressed();
}//createChooseDrfDialog(...)


void DrfSelect::handle_app_url_drf( const std::string &url_query )
{
  try
  {
    auto drf = make_shared<DetectorPeakResponse>();
    drf->fromAppUrl( url_query );
    
    assert( drf->isValid() );
    
    DrfSelect::createChooseDrfDialog( {drf}, "DRF received from external program.", "" );
  }catch( std::exception &e )
  {
    wApp->log( "error" ) << "App URL was invalid DRF: " << e.what();
    throw runtime_error( WString::tr("ds-err-drf-uri-invalid").toUTF8() );
  }
}//void DrfSelect::handle_app_url_drf( const std::string &url_query )


void DrfSelect::setAcceptButtonEnabled( const bool enable )
{
  m_acceptButton->setEnabled( enable );
  m_noDrfButton->setHidden( !m_detector );
  m_xmlDownload->setHidden( !m_detector );
  if( m_defaultForSerialNumber )
    m_defaultForSerialNumber->setEnabled( enable );
  if( m_defaultForDetectorModel )
    m_defaultForDetectorModel->setEnabled( enable );
}


void DrfSelect::handleUserChangedUploadedDrfName()
{
  if( !m_uploadedDetName || !m_detector )
    return;
  
  string value = m_uploadedDetName->text().toUTF8();
  if( value.empty() )
  {
    const int offset = wApp->environment().timeZoneOffset();
    auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    now += std::chrono::seconds(60*offset);
    value = SpecUtils::to_vax_string(now);
    m_uploadedDetName->setText( WString::fromUTF8(value) );
  }//if( value.empty() )
  
  m_detector->setName( value );
}//handleUserChangedUploadedDrfName()


void DrfSelect::handleFitFwhmRequested()
{
  const bool use_auto_fit_peaks_too = true;
  MakeFwhmForDrfWindow *tool = m_interspec->fwhmFromForegroundWindow( use_auto_fit_peaks_too );
}//void handleFitFwhmRequested()


void DrfSelect::handleFitFwhmFinished( std::shared_ptr<DetectorPeakResponse> drf )
{
  if( drf )
    updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql );
  done().emit();
}//void handleFitFwhmFinished( std::shared_ptr<DetectorPeakResponse> drf )


void DrfSelect::deleteDBTableSelected()
{
  WModelIndexSet indices = m_DBtable->selectedIndexes();
  if( !indices.size() )
  {
    m_deleteButton->disable();
    return;
  }//if( !indices.size() )
  
  WModelIndex index = *indices.begin();
  DataBaseUtils::DbTransaction transaction( *m_sql );
  Dbo::ptr<DetectorPeakResponse> dbfile = m_model->stableResultRow( index.row() );
  dbfile.remove();
  transaction.commit();
  m_model->reload();
} //deleteDBTableSelected

// Grabs the detector if selected in DB table
void DrfSelect::dbTableSelectionChanged()
{
  if( m_DBtable->selectedIndexes().size() )
    m_deleteButton->enable();
  
  bool failed = false;
  std::shared_ptr<DetectorPeakResponse> det;

  try
  {
    if( !m_model  )
      return;
    
    const WModelIndexSet indices = m_DBtable->selectedIndexes();
    
    if( indices.empty() )
    {
      m_deleteButton->disable();
      return;
    }//if( !indices.empty() )
    
    m_deleteButton->enable();
    const WModelIndex index = *indices.begin();
    Dbo::ptr<DetectorPeakResponse> dbfile = m_model->stableResultRow( index.row() );
    det = std::shared_ptr<DetectorPeakResponse>(new DetectorPeakResponse(*dbfile));
  }catch(...)
  {
    failed = true;
    passMessage( "Error getting from DetectorPeakResponse table", WarningWidget::WarningMsgHigh );
  }//try / catch
  
  m_detector = det;
  m_gui_select_matches_det = true;
  setAcceptButtonEnabled( !failed );
  emitChangedSignal();
}//void dbTableSelectionChanged()

void DrfSelect::setGuiToCurrentDetector()
{
  //----initialize
  
  if( !m_detector )
  {
    //m_drfTypeMenu->select( 0 );
  }else
  {
    switch( m_detector->drfSource() )
    {
      case DetectorPeakResponse::DefaultGadrasDrf:
      case DetectorPeakResponse::UserAddedGadrasDrf:
      {
        m_drfTypeMenu->select( 0 );
        
        const MatchDetectorStatus status = m_gadrasDetSelect->trySelectDetector( m_detector );
        
        m_gui_select_matches_det = (status == MatchDetectorStatus::Match);
        break;
      }
      
      case DetectorPeakResponse::UserAddedRelativeEfficiencyDrf:
      case DetectorPeakResponse::DefaultRelativeEfficiencyDrf:
      {
        m_drfTypeMenu->select( 1 );
        const MatchDetectorStatus status = m_relEffSelect->trySelectDetector( m_detector );
        
        m_gui_select_matches_det = (status == MatchDetectorStatus::Match);
        break;
      }
        
      case DetectorPeakResponse::UserSpecifiedFormulaDrf:
      {
        m_drfTypeMenu->select( 3 );
        m_detectorManualFunctionName->setText( m_detector->name() );
        m_detectorManualFunctionText->setText( m_detector->efficiencyFormula() );
        m_detectorManualDescription->setText( m_detector->description() );
        const float diam = m_detector->detectorDiameter();
        const string diamstr = PhysicalUnits::printToBestLengthUnits( diam );
        m_detectorManualDiameterText->setText( diamstr );
        const float units = m_detector->efficiencyEnergyUnits();
        int energygrp = 0;
        if( fabs(PhysicalUnits::MeV-units) < fabs(PhysicalUnits::keV-units) )
          energygrp = 1;
        m_eqnEnergyGroup->setSelectedButtonIndex( energygrp );
        m_drfType->setCurrentIndex( 0 );
        
        m_gui_select_matches_det = true;
        break;
      }//case kUserEfficiencyEquationSpecified:
        
      case DetectorPeakResponse::UserImportedIntrisicEfficiencyDrf:
      {
        m_drfTypeMenu->select( 2 );
        m_gui_select_matches_det = true;
        break;
      }
      
      case DetectorPeakResponse::UnknownDrfSource:
      case DetectorPeakResponse::UserImportedGadrasDrf:
      case DetectorPeakResponse::UserCreatedDrf:
      case DetectorPeakResponse::IsocsEcc:
      case DetectorPeakResponse::FromSpectrumFileDrf:
      {
        m_drfTypeMenu->select( 4 );
        try
        {
          const double start_time = SpecUtils::get_wall_time();
          for( int row = 0; row < m_model->rowCount(); ++row )
          {
            Wt::Dbo::ptr<DetectorPeakResponse> drf = m_model->resultRow(row);
            if( drf && (drf->hashValue() == m_detector->hashValue()) )
            {
              WModelIndexSet indexset;
              indexset.insert( m_model->index(row,0) );
              m_DBtable->setSelectedIndexes( indexset );
              m_gui_select_matches_det = true;
              break;
            }
            
            //On a quick test, it took about 25ms when the very first row
            //  matched the detector; even if we loop through all rows, it
            //  takes about the same amount of time
            const double now_time = SpecUtils::get_wall_time();
            if( (now_time-start_time) > 0.25 )
              break;
          }//for( loop over previous DRFs )
        }catch( std::exception &e )
        {
          cerr << "Caught exception trying to select previous DRF from database to select for the GUI: " << e.what() << endl;
        }//try / catch
        break;
      }
    };//switch( m_detector->drfSource() )
  }//if( m_detector == NULL ) / else
  
  selectButton( m_drfTypeStack, m_drfTypeMenu, false );
}// setGuiToCurrentDetector()


void DrfSelect::detectorsWereInited()
{
  // If we try to set the GUI state to match m_detector, it may cause the tab on the left
  //  (i.e., "GADRAS", "Rel. Eff.", "Import", etc) to be changed, which may be overiding
  //  what the user wants.
  if( !m_gui_select_matches_det )
    setGuiToCurrentDetector();
}//void detectorsWereInited();


DrfSelect::~DrfSelect()
{
} //~DrfSelect()


void DrfSelect::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  if( m_detector != det )
  {
    m_detector = det;
    m_gui_select_matches_det = false;
    setGuiToCurrentDetector();
  }
}//void setDetector( std::shared_ptr<DetectorPeakResponse> det )


void DrfSelect::verifyManualDefinition()
{ //if we change anything in the text fields, we should disable Accept
  setAcceptButtonEnabled( false );
  
  bool valid = true;
  float lowerEnergy = 0.0f, upperEnergy = 0.0f;
  
  if( !m_detectorManualMinEnergy->text().empty()
     && (m_detectorManualMinEnergy->validate() == WValidator::State::Valid) )
  {
    lowerEnergy = m_detectorManualMinEnergy->value();
  }else if( !m_detectorManualMinEnergy->text().empty() )
  {
    valid = false;
  }
  
  if( !m_detectorManualMaxEnergy->text().empty()
     && (m_detectorManualMaxEnergy->validate() == WValidator::State::Valid) )
  {
    upperEnergy = m_detectorManualMaxEnergy->value();
  }else if( !m_detectorManualMaxEnergy->text().empty() )
  {
    valid = false;
  }
  
  
  string fcn_txt = m_detectorManualFunctionText->text().toUTF8();
  
  const CoefToFormulaInfo coef_formula_info = check_if_coef_for_formula( fcn_txt );
  if( coef_formula_info.m_valid_coefficients )
  {
    fcn_txt = coef_formula_info.m_fcn;
    
    if( coef_formula_info.m_energy_units == PhysicalUnits::keV )
      m_eqnEnergyGroup->setSelectedButtonIndex( 0 );
    else if( coef_formula_info.m_energy_units == PhysicalUnits::MeV )
      m_eqnEnergyGroup->setSelectedButtonIndex( 1 );
    
    if( coef_formula_info.m_upper_energy > coef_formula_info.m_lower_energy )
    {
      if( m_detectorManualMinEnergy->text().empty() )
        m_detectorManualMinEnergy->setValue( coef_formula_info.m_lower_energy );
      
      if( m_detectorManualMaxEnergy->text().empty() )
        m_detectorManualMaxEnergy->setValue( coef_formula_info.m_upper_energy );
    }//if( coef_formula_info.m_upper_energy > coef_formula_info.m_lower_energy )
    
    if( !coef_formula_info.m_det_diam.empty()
       && (m_detectorManualDiameterText->text().empty()
           || (m_detectorManualDiameterText->text().toUTF8() == "2.2 cm") ) )
    {
      m_detectorManualDiameterText->setText( coef_formula_info.m_det_diam );
    }
    
    if( !coef_formula_info.m_source_to_crystal_distance.empty()
       && (m_detectorManualDistText->isHidden() || m_detectorManualDistText->text().empty() )
       && !coef_formula_info.m_fixed_geom )
    {
      m_drfType->setCurrentIndex( 1 );
      m_detectorManualDistText->setText( coef_formula_info.m_source_to_crystal_distance );
    }
    
    const WString &pre_det_wname = m_detectorManualFunctionName->text();
    const string prev_det_name = pre_det_wname.toUTF8();
    
    
    if( !coef_formula_info.m_det_name.empty()
       && (prev_det_name.empty() || (pre_det_wname.key() == "app-default-manual-det-name")) )
    {
      m_detectorManualFunctionName->setText( coef_formula_info.m_det_name );
    }
    
    if( coef_formula_info.m_fixed_geom )
      m_drfType->setCurrentIndex( 3 ); //Act/cm2
  }//if( user_entered_coefs )
  
  const bool is_intrinsic = (m_drfType->currentIndex() == 0);
  const bool is_absolute = (m_drfType->currentIndex() == 1);
  const bool is_fixed_geom = (m_drfType->currentIndex() >= 2);
  
  assert( (static_cast<int>(is_intrinsic)
           + static_cast<int>(is_absolute)
           + static_cast<int>(is_fixed_geom)) == 1 );
  DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarField;
  switch( m_drfType->currentIndex() )
  {
    case 0:  geom_type = DetectorPeakResponse::EffGeometryType::FarField; break;
    case 1:  geom_type = DetectorPeakResponse::EffGeometryType::FarField; break;
    case 2:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct; break;
    case 3:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2; break;
    case 4:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2; break;
    case 5:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default:
      assert( 0 );
      break;
  }//switch( m_drfType->currentIndex() )
  
  m_detectorManualDistLabel->setHidden( is_intrinsic || is_fixed_geom );
  m_detectorManualDistText->setHidden( is_intrinsic || is_fixed_geom );
  m_detectorManualDiameterText->setHidden( is_fixed_geom );
  m_detectorManualDiameterLabel->setHidden( is_fixed_geom );
  
  if( is_absolute )
  {
    if( m_detectorManualDistText->text().empty()
       || (m_detectorManualDistText->validate() != WValidator::State::Valid) )
      valid = false;
  }//if( is_absolute )
  
  if( !is_fixed_geom && (m_detectorManualDiameterText->validate() != WValidator::State::Valid) )
    valid = false;
  
  
  double det_diam = 5.08*PhysicalUnits::cm;
  
  try
  {
    const string diamtxt = m_detectorManualDiameterText->text().toUTF8();
    det_diam = PhysicalUnits::stringToDistance( diamtxt );
  }catch( std::exception & )
  {
    if( is_fixed_geom )
      det_diam = 0.0;
    else
      valid = false;
  }
  
  if( is_absolute )
  {
    double dist = 1*PhysicalUnits::meter;
    try
    {
      const string disttxt = m_detectorManualDistText->text().narrow();
      dist = static_cast<float>( PhysicalUnits::stringToDistance( disttxt ) );
      if( dist < 0.001*PhysicalUnits::cm )
        throw exception();
    }catch( std::exception & )
    {
      valid = false;
      if( !m_detectorManualDistText->hasStyleClass( "Wt-invalid" ) )
        m_detectorManualDistText->addStyleClass( "Wt-invalid" );
    }//try / catch
    
    const double gfactor = DetectorPeakResponse::fractionalSolidAngle( det_diam, dist );
    if( !fcn_txt.empty() )
      fcn_txt = std::to_string(1.0/gfactor) + "*(" + fcn_txt + ")";
  }//if( absolute eff )
  
  const float energyUnits = static_cast<float>( m_eqnEnergyGroup->checkedId()==0
                                               ? PhysicalUnits::keV : PhysicalUnits::MeV);
  
  if( upperEnergy <= lowerEnergy )
    lowerEnergy = upperEnergy = 0.0f;

  try
  {
    if( fcn_txt.empty() )
      throw exception();
    
    DetectorPeakResponse detec( "temp", "temp " );
    detec.setIntrinsicEfficiencyFormula( fcn_txt, static_cast<float>(det_diam),
                                        energyUnits, lowerEnergy, upperEnergy, geom_type );
    
    if( m_detectorManualFunctionText->hasStyleClass("Wt-invalid") )
      m_detectorManualFunctionText->removeStyleClass( "Wt-invalid" );
  }catch( std::exception & )
  {
    valid = false;
    if( m_detectorManualFunctionText->text().empty() )
    {
      if( m_detectorManualFunctionText->hasStyleClass("Wt-invalid") )
        m_detectorManualFunctionText->removeStyleClass( "Wt-invalid" );
    }else
    {
      if( !m_detectorManualFunctionText->hasStyleClass("Wt-invalid") )
        m_detectorManualFunctionText->addStyleClass( "Wt-invalid" );
    }
  }//try / catch

  m_manualSetButton->setEnabled( valid );
}//verifyManualDefinition()


void DrfSelect::setFormulaDefineDetector()
{
  string fcn = m_detectorManualFunctionText->text().narrow();
  
  float coef_units = 0.0;
  const CoefToFormulaInfo coef_formula_info = check_if_coef_for_formula( fcn );
  if( coef_formula_info.m_valid_coefficients )
  {
    fcn = coef_formula_info.m_fcn;
    m_detectorManualFunctionText->setText( WString::fromUTF8(fcn) );
  }//if( check_if_coef_for_formula(fcn, coef_units) )
  
  const string name = m_detectorManualFunctionName->text().toUTF8();
  string descr = m_detectorManualDescription->text().toUTF8();
  if( descr.empty() )
  {
    descr = "f(x)=" + fcn;
    if( descr.size() > 25 )
      descr = descr.substr( 0, 25 ) + "...";
    m_detectorManualDescription->setText( descr );
  }//if( descr.empty() )

  const bool is_intrinsic = (m_drfType->currentIndex() == 0);
  const bool is_absolute = (m_drfType->currentIndex() == 1);
  const bool is_fixed_geom = (m_drfType->currentIndex() >= 2);
  
  assert( (static_cast<int>(is_intrinsic)
           + static_cast<int>(is_absolute)
           + static_cast<int>(is_fixed_geom)) == 1 );
  DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarField;
  switch( m_drfType->currentIndex() )
  {
    case 0:  geom_type = DetectorPeakResponse::EffGeometryType::FarField; break;
    case 1:  geom_type = DetectorPeakResponse::EffGeometryType::FarField; break;
    case 2:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct; break;
    case 3:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2; break;
    case 4:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2; break;
    case 5:  geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram; break;
    default:
      assert( 0 );
      break;
  }//switch( m_drfType->currentIndex() )
  
  assert( (static_cast<int>(is_intrinsic)
           + static_cast<int>(is_absolute)
           + static_cast<int>(is_fixed_geom)) == 1 );
  
  
  shared_ptr<DetectorPeakResponse> detec = make_shared<DetectorPeakResponse>( name, descr );
  
  try
  {
    const string diamtxt = m_detectorManualDiameterText->text().narrow();
    double det_diam = 0.0;
    try
    {
      det_diam = PhysicalUnits::stringToDistance( diamtxt );
    }catch( std::exception &e )
    {
      if( !is_fixed_geom )
        throw;
    }
    
    if( !is_fixed_geom && (det_diam < (is_fixed_geom ? 0.0 : 0.0001)) )
      throw runtime_error( WString::tr("ds-err-det-diam-pos").toUTF8() );
    
    if( is_absolute )
    {
      const string disttxt = m_detectorManualDistText->text().narrow();
      const float dist = (float)PhysicalUnits::stringToDistance( disttxt );
      if( dist < 0.0001f )
        throw runtime_error( WString::tr("ds-err-det-dist-pos").toUTF8() );
      
      
      const double gfactor = DetectorPeakResponse::fractionalSolidAngle(
                                                              det_diam, dist );
        
      fcn = std::to_string(1.0/gfactor) + "*(" + fcn + ")";
    }//if( is_absolute )
    
    const float energyUnits
              = static_cast<float>(m_eqnEnergyGroup->checkedId()==0
                                    ? PhysicalUnits::keV : PhysicalUnits::MeV);
    
    float minEnergy = 0.0f, maxEnergy = 0.0f;
    if( !m_detectorManualMinEnergy->text().empty()
        && (m_detectorManualMinEnergy->validate() == WValidator::State::Valid) )
      minEnergy = m_detectorManualMinEnergy->value();
    
    if( !m_detectorManualMaxEnergy->text().empty()
        && (m_detectorManualMaxEnergy->validate() == WValidator::State::Valid) )
      maxEnergy = m_detectorManualMaxEnergy->value();
    
    if( maxEnergy < minEnergy )
    {
      minEnergy = maxEnergy = 0.0f;
      m_detectorManualMinEnergy->setValueText( "" );
      m_detectorManualMaxEnergy->setValueText( "" );
      passMessage( WString::tr("ds-warn-min-max-energy-swapped"), WarningWidget::WarningMsgHigh );
    }//if( maxEnergy < minEnergy )
    
    if( m_detectorManualDiameterText->hasStyleClass("Wt-invalid") )
      m_detectorManualDiameterText->removeStyleClass( "Wt-invalid" );
    
    try
    {
      detec->setIntrinsicEfficiencyFormula( fcn, det_diam, energyUnits,
                                           minEnergy, maxEnergy, geom_type );
    }catch( std::exception &e )
    {
      if( !fcn.empty() && !m_detectorManualDiameterText->hasStyleClass("Wt-invalid") )
        m_detectorManualDiameterText->addStyleClass( "Wt-invalid" );
      throw;
    }
    
    detec->setDrfSource( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf );
    
//#if( PERFORM_DEVELOPER_CHECKS )
//    check_url_serialization( detec );
//#endif
    
    updateChart();
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    setAcceptButtonEnabled( false );
    m_manualSetButton->enable();
    updateChart();
    return;
  }//try / catch
  
  m_manualSetButton->disable();
  m_detector = detec;
  m_gui_select_matches_det = true;
  setAcceptButtonEnabled( true );

  emitChangedSignal();
}//setFormulaDefineDetector(const std::string fcn)



// This calls the action to initialize the charts to whatever
// is selected.  Each stack page should have something to update
void DrfSelect::selectButton( WStackedWidget *stack,
                                 WMenu *menu,
                                 bool activateCallBack)
{
  //makes one of the stack visible
  stack->setCurrentIndex( menu->currentIndex() );  //I think this is now not needed
  
  if( activateCallBack )
  {
    //clear the chart when we change, regardless
    m_detector = nullptr;
    m_gui_select_matches_det = false;
    
    switch( menu->currentIndex() )
    {
      case 0:
        gadrasDetectorSelectCallback();
        break;
    
      case 1:
        relEffDetectorSelectCallback();
        break;
        
      case 2:
        fileUploadedCallback( UploadCallbackReason::ImportTabChosen );
        break;
        
      case 3:
        setFormulaDefineDetector();
        break;
        
      case 4:
        dbTableSelectionChanged();
        break;
    } //switch
  }//if( activateCallBack )
  
  updateChart();
} //DrfSelect::selectButton(Wt::WStackedWidget * stack, int index)


void DrfSelect::updateChart()
{
  m_chart->updateChart( m_detector );
} //DrfSelect::updateChart()


std::shared_ptr<DetectorPeakResponse> DrfSelect::parseRelEffCsvFile( const std::string filename )
{
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream csvfile( wfilename.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream csvfile( filename.c_str(), ios_base::binary | ios_base::in );
#endif
  
  if( !csvfile.is_open() )
    return nullptr;
  
  string line;
  int nlineschecked = 0;
  
  string drfname, drfdescrip;
  bool fixed_geometry = false;
  bool foundMeV = false, foundKeV = false;
  
  //ToDo: Need to implement getting lines safely where a quoted field may span
  //      several lines.
  while( SpecUtils::safe_get_line(csvfile, line, 2048) && (++nlineschecked < 100) )
  {
    foundKeV |= SpecUtils::icontains( line, "kev" );
    foundMeV |= SpecUtils::icontains( line, "mev" );
    fixed_geometry |= SpecUtils::icontains( line, "Fixed Geometry DRF" );
    
    vector<string> fields;
    split_escaped_csv( fields, line );
    
    if( fields.size() == 2 )
    {
      if( fields[0] == "# Name" )
        drfname = fields[1];
      if( fields[0] == "# Description" )
        drfdescrip = fields[1];
    }
    
    if( fields.size() < 16 )
      continue;
    
    if( !SpecUtils::iequals_ascii( fields[3], "c0")
       || !SpecUtils::iequals_ascii( fields[4], "c1")
       || !SpecUtils::iequals_ascii( fields[5], "c2")
       || !SpecUtils::iequals_ascii( fields[6], "c3")
       || !SpecUtils::icontains( fields[15], "radius") )
      continue;
    
    // Nominally the header is "FixedGeometry", and be in column 17, but we'll be a little loose
    int fixed_geom_col = -1;
    for( size_t i = 0; (fixed_geom_col < 0) && (i < fields.size()); ++i )
    {
      if( SpecUtils::icontains( fields[i], "Fixed") && SpecUtils::icontains( fields[i], "Geom") )
        fixed_geom_col = static_cast<int>( i );
    }
    
    //Okay, next line should be
    if( !SpecUtils::safe_get_line(csvfile, line, 2048) )
      throw runtime_error( "Couldnt get next line" );
    
    split_escaped_csv( fields, line );
    
    try
    {
      vector<float> coefs( 8, 0.0f ), coef_uncerts;
      for( int i = 0; i < 8; ++i )
      {
        string field = fields.at(3+i);
        coefs[i] = (field.empty() ? 0.0f : std::stof(field));
      }
      
      for( int i = 7; i >= 0; --i )
      {
        if( coefs[i] != 0.0f )
          break;
        else
          coefs.resize( coefs.size() - 1 );
      }
      
      if( coefs.size() < 1 )
        continue;
      
      const float dist = std::stof( fields.at(14) ) * PhysicalUnits::cm;
      const float radius = std::stof( fields.at(15) ) * PhysicalUnits::cm;
      
      if( (fixed_geom_col >= 0) && (static_cast<int>(fields.size()) >= fixed_geom_col) )
      {
        fixed_geometry |= (SpecUtils::icontains( fields[fixed_geom_col], "1")
                           || SpecUtils::icontains( fields[fixed_geom_col], "yes")
                           || SpecUtils::icontains( fields[fixed_geom_col], "true")
                           || SpecUtils::istarts_with(fields[fixed_geom_col], "y"));
      }//if( has_fixed_geom_col )
      
      
      const string name = (fields[0].empty() ? drfname : fields[0]);
      const float energUnits = ((foundKeV && !foundMeV) ? 1.0f : 1000.0f);
      float lowerEnergy = 0.0f, upperEnergy = 0.0f;
      
      // Lets try to get coefs uncertainties
      if( !SpecUtils::safe_get_line(csvfile, line, 2048) )
      {
        try
        {
          vector<string> uncert_strs;
          split_escaped_csv( uncert_strs, line );
          if( !uncert_strs.empty()
             && SpecUtils::istarts_with(uncert_strs[0], "# 1 sigma Uncert")
             && (uncert_strs.size() >= (3 + coefs.size())) )
          {
            coef_uncerts.resize( coefs.size(), 0.0f );
            for( int i = 0; i < coefs.size(); ++i )
            {
              string field = uncert_strs.at(3+i);
              coef_uncerts[i] = (field.empty() ? 0.0f : std::stof(field));
            }
          }//if( it looks like we found the uncertainty line )
        }catch( std::exception &e )
        {
          cerr << "Error caught parsing DRF eff uncertainties: " << e.what() << endl;
          coef_uncerts.clear();
        }
      }//if( we got another line we'll check if its the uncertainties )
      
      DetectorPeakResponse::EffGeometryType geom_type = DetectorPeakResponse::EffGeometryType::FarField;
      if( fixed_geometry )
        geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;
      
      auto det = std::make_shared<DetectorPeakResponse>( fields[0], drfdescrip );
      
      det->fromExpOfLogPowerSeriesAbsEff( coefs, coef_uncerts, dist, 2.0f*radius, energUnits,
                                          lowerEnergy, upperEnergy, geom_type );
      
      //Look for the line that gives the appropriate energy range.
#define POS_DECIMAL_REGEX "\\+?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?\\s*"
      const char * const rng_exprsn_txt = "Valid energy range:\\s*(" POS_DECIMAL_REGEX ")\\s*keV to\\s*(" POS_DECIMAL_REGEX ")\\s*keV.";
      
      std::regex range_expression( rng_exprsn_txt );
      
      while( SpecUtils::safe_get_line(csvfile, line, 2048) && (++nlineschecked < 100) )
      {
        if( SpecUtils::icontains( line, "Full width half maximum (FWHM) follows equation" ) )
        {
          const bool isConstPlusSqrt = SpecUtils::icontains( line, "A0 + A1*sqrt" );
          const bool isSqrt = !isConstPlusSqrt && SpecUtils::icontains( line, "sqrt(" );  //Else contains "GadrasEqn"
          auto form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
          if( isConstPlusSqrt )
            form = DetectorPeakResponse::ResolutionFnctForm::kConstantPlusSqrtEnergy;
          else if( isSqrt )
            form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
          
          int nlinecheck = 0;
          while( SpecUtils::safe_get_line(csvfile, line, 2048)
                && !SpecUtils::icontains(line, "Values")
                && (++nlinecheck < 15) )
          {
          }
        
          vector<string> fwhm_fields;
          split_escaped_csv( fwhm_fields, line );
          
          if( fwhm_fields.size() > 1 && SpecUtils::icontains(fwhm_fields[0], "Values") )
          {
            try
            {
              vector<float> coefs;
              for( size_t i = 1; i < fwhm_fields.size(); ++i )
                coefs.push_back( stof(fwhm_fields[i]) );
              if( !coefs.empty() )
                det->setFwhmCoefficients( coefs, form );
            }catch(...)
            {
            }
            break;
          }
        }//if( start of FWHM section of CSV file )
        
        std::smatch range_matches;
        if( std::regex_search( line, range_matches, range_expression ) )
        {
          lowerEnergy = std::stof( range_matches[1] );
          upperEnergy = std::stof( range_matches[6] );
          det->setEnergyRange( lowerEnergy, upperEnergy );
          break;
        }
      }//while( getline )
      
//#if( PERFORM_DEVELOPER_CHECKS )
//      check_url_serialization( det );
//#endif

      return det;
    }catch(...)
    {
      continue;
    }
  }//while( more lines )
  
  return nullptr;
}//parseRelEffCsvFile(...)



void DrfSelect::fileUploadedCallback( const UploadCallbackReason context )
{
  auto updateUserName = [this](){
    
    string userDrfFilename = m_efficiencyCsvUpload->clientFileName().toUTF8();
    if( SpecUtils::iends_with( userDrfFilename, ".csv" ) )
      userDrfFilename = userDrfFilename.substr(0, userDrfFilename.size()-4);
    if( userDrfFilename.empty() )
      userDrfFilename = "Uploaded";
    
    // Efficiency.csv is a really common name, as is a few other short names; lets
    //  add current date/time to these short names as *some* type of differentiator.
    //  Its something.
    if( userDrfFilename.size() < 15 )
    {
      const int offset = wApp->environment().timeZoneOffset();
      auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      now += chrono::seconds(60*offset);
      userDrfFilename += " " + SpecUtils::to_vax_string(now);
    }
    
    if( m_detector && m_uploadedDetName )
    {
      m_detector->setName( userDrfFilename );
      m_uploadedDetName->setText( WString::fromUTF8(userDrfFilename) );
    }
  };//updateUserName(...)
  
  bool is_fixed_geometry = false;
  m_uploadedDetName->setValueText( "" );
  if( !m_efficiencyCsvUpload->empty() )
  {
    m_detectrDiameterDiv->show();
  }else
  {
    m_detectrDiameterDiv->hide();
  }

  switch( context )
  {
    case UploadCallbackReason::ImportTabChosen:
      break;
      
    case UploadCallbackReason::DetectorDiameterChanged:
      break;
      
    case UploadCallbackReason::FixedGeometryChanged:
      m_detectorDiameter->setHidden( m_fixedGeometryCb->isChecked() );
      break;
      
    case UploadCallbackReason::DetectorDotDatUploaded:
      break;
      
    case UploadCallbackReason::EfficiencyCsvUploaded:
      if( !m_efficiencyCsvUpload->empty() )
      {
        const string filename = m_efficiencyCsvUpload->spoolFileName();
        auto det = DrfSelect::parseRelEffCsvFile( filename );
        
        if( det )
        {
          auto detDiamStr = PhysicalUnits::printToBestLengthUnits( det->detectorDiameter() );
          m_detectorDiameter->setText( detDiamStr );
          m_detectrDotDatDiv->hide();
          m_detectorDiameter->disable();
          m_detector = det;
          
          const bool fixed_geometry = det->isFixedGeometry();
          m_fixedGeometryCb->show();
          m_fixedGeometryCb->setChecked( fixed_geometry );
          m_fixedGeometryCb->disable();
          m_detectrDiameterDiv->setHidden( fixed_geometry );
          
          setAcceptButtonEnabled( true );
          updateUserName();
          emitChangedSignal();
          return;
        }else
        {
          m_detectrDotDatDiv->show();
          m_detectrDiameterDiv->show();
          m_fixedGeometryCb->hide();
          m_fixedGeometryCb->setChecked( false );
          m_detectorDiameter->enable();
        }//if( det )
      }//if( we have a file to test )
      break;
  }//switch( context )

  const bool isGadrasDet = (!m_efficiencyCsvUpload->empty()
                            && !m_detectorDotDatUpload->empty()
                            && m_efficiencyCsvUpload->spoolFileName().size()
                            && m_detectorDotDatUpload->spoolFileName().size());

  const bool fixed_geometry = (m_fixedGeometryCb->isVisible() && m_fixedGeometryCb->isChecked());
  
  float diameter = -1.0f;
  if( !fixed_geometry )
  {
    try
    {
      diameter = static_cast<float>( PhysicalUnits::stringToDistance( m_detectorDiameter->text().toUTF8() ) );
      if( diameter < float(0.001*PhysicalUnits::cm) )
        diameter = 0.0f;
    }catch(...)
    {
      m_detectorDiameter->setText( "" );
    }
  }//if( !fixed_geometry )
  
  const bool isDiamDet = (!m_efficiencyCsvUpload->empty()
                          && diameter>0.0
                          && !fixed_geometry
                          && m_efficiencyCsvUpload->spoolFileName().size() );

  if( !isGadrasDet && !isDiamDet && !fixed_geometry )
  {
    if( diameter > 0.0 )
      passMessage( WString::tr("ds-need-det-diam"), WarningWidget::WarningMsgHigh );
    return;
  }//if( !isGadrasDet && !isDiamDet )

  const string csvFileName = m_efficiencyCsvUpload->spoolFileName();

#ifdef _WIN32
  const std::wstring wcsvFileName = SpecUtils::convert_from_utf8_to_utf16(csvFileName);
  ifstream csvfile( wcsvFileName.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream csvfile( csvFileName.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !csvfile.is_open() )
    return;

  auto det = std::make_shared<DetectorPeakResponse>();
  try
  {
    if( isGadrasDet )
    {
      const string dotDatFileName = m_detectorDotDatUpload->spoolFileName();
#ifdef _WIN32
      const std::wstring wdotDatFileName = SpecUtils::convert_from_utf8_to_utf16(dotDatFileName);
      ifstream datfile( wdotDatFileName.c_str(), ios_base::binary|ios_base::in );
#else
      ifstream datfile( dotDatFileName.c_str(), ios_base::binary|ios_base::in );
#endif
      if( !datfile.is_open() )
        return;
      det->fromGadrasDefinition( csvfile, datfile );
      det->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedGadrasDrf );
      if( context == UploadCallbackReason::DetectorDotDatUploaded )
        m_detectorDiameter->setText( PhysicalUnits::printToBestLengthUnits(det->detectorDiameter()) );
    }else
    {
      det->fromEnergyEfficiencyCsv( csvfile, diameter, float(PhysicalUnits::keV), DetectorPeakResponse::EffGeometryType::FarField );
      det->setDrfSource( DetectorPeakResponse::DrfSource::UserImportedIntrisicEfficiencyDrf );
      
      m_fixedGeometryCb->show();
      m_fixedGeometryCb->enable();
    }//if( isGadrasDet ) / else
  }catch( std::exception &e )
  {
    setAcceptButtonEnabled( false );
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
//#if( PERFORM_DEVELOPER_CHECKS )
//  check_url_serialization( det );
//#endif

  
  m_detector = det;
  setAcceptButtonEnabled( true );
  updateUserName();
  emitChangedSignal();
}//void fileUploadedCallback();


void DrfSelect::relEffDetectorSelectCallback()
{
  m_relEffSelect->userSelectedRelEffDetSelect();
}//void relEffDetectorSelectCallback()


void DrfSelect::gadrasDetectorSelectCallback()
{
  std::shared_ptr<DetectorPeakResponse> det = m_gadrasDetSelect->selectedDetector();

  m_detector = det;
  m_gui_select_matches_det = true;
  setAcceptButtonEnabled( !!det );
  
  emitChangedSignal();
}//void gadrasDetectorSelectCallback()


void DrfSelect::updateLastUsedTimeOrAddToDb( std::shared_ptr<DetectorPeakResponse> drf,
                                                long long db_user_id,
                                                std::shared_ptr<DataBaseUtils::DbSession> sql )
{
  if( !drf || !sql )
    return;
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    uint64_t hash = drf->hashValue();
    auto result = sql->session()->find<DetectorPeakResponse>()
                      .where("InterSpecUser_id = ? AND Hash = ?")
                      .bind( db_user_id )
                      .bind( reinterpret_cast<int64_t&>(hash) )
                      .resultList();
    int nupdated = 0;
    for( auto iter = result.begin(); iter != result.end(); ++iter )
    {
      (*iter).modify()->updateLastUsedTimeToNow();
      ++nupdated;
    }//
    
//    cout << "Have seen Detector '" << drf->name() << "' " << nupdated << " times." << endl;
    
    if( !nupdated )
    {
      //Create a separate DetectorPeakResponse because shared_ptr and
      //  dbo::ptr don't work well together
      DetectorPeakResponse *tempDetector = new DetectorPeakResponse( *drf );
      tempDetector->updateLastUsedTimeToNow();
      
      tempDetector->m_user = db_user_id; //adds the current user to the detectorpeakresponse insertion into DB
      sql->session()->add( tempDetector );
    }
    
    transaction.commit();
  }catch( std::exception &e )
  {
    passMessage( "Error writing DetectorPeakResponse to DB: " + std::string(e.what()), WarningWidget::WarningMsgHigh );
  }//try / catch
}//void updateLastUsedTimeOrAddToDb(...)



std::shared_ptr<DetectorPeakResponse> DrfSelect::getUserPreferredDetector(
                                std::shared_ptr<DataBaseUtils::DbSession> sql,
                                Wt::Dbo::ptr<InterSpecUser> user,
                                const std::string &serial_number,
                                SpecUtils::DetectorType detType,
                                const std::string &detector_model )
{
  std::shared_ptr<DetectorPeakResponse> answer;
  if( !sql || !user )
    return answer;
  
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    Wt::Dbo::ptr<UseDrfPref> pref;
    
    if( !serial_number.empty() )
    {
      auto results = sql->session()->find<UseDrfPref>()
                     .where( "InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?" )
                     .bind( user.id() )
                     .bind( UseDrfPref::UseDrfType::UseDetectorSerialNumber )
                     .bind( serial_number )
                     .orderBy("id desc") //shouldnt have an effect because we should only get at most one result
                     .resultList();
      if( results.size() )
        pref = results.front();
    }
    
    if( !pref )
    {
      string model;
      if( detType == DetectorType::Unknown )
        model = detector_model;
      else
        model = detectorTypeToString( detType );
      if( !model.empty() )
      {
        auto results = sql->session()->find<UseDrfPref>()
                       .where( "InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?" )
                       .bind( user.id() )
                       .bind( UseDrfPref::UseDrfType::UseDetectorModelName )
                       .bind( model )
                       .orderBy("id desc")  //shouldnt have an effect because we should only get at most one result
                       .resultList();
        
        if( results.size() )
          pref = results.front();
      }
    }//if( !results.size() )
    
    if( !pref )
    {
      //cout << "Could not find default DRF preference in database" << endl;
      transaction.commit();
      return answer;
    }
    
    if( pref->m_drfIndex < 0 )
    {
      //cout << "Default DRF preference in database was invalid" << endl;
      pref.remove();
      transaction.commit();
      return answer;
    }
    
    auto drflist = sql->session()->find<DetectorPeakResponse>()
    .where("InterSpecUser_id = ? AND id = ?")
    .bind( user.id() )
    .bind( pref->m_drfIndex )
    .resultList();
    
    for( auto iter = drflist.begin(); !answer && iter != drflist.end(); ++iter )
      answer = std::make_shared<DetectorPeakResponse>( **iter );
    
    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "Got exception getting user preffered detector from DB: " << e.what() << endl;
  }//try catch see if the user has a preference
  
  //if( answer )
  //  cout << "Got user preffered default det" << endl;
  //else
  //  cout << "Did not get user preffered default det" << endl;
  
  return answer;
};//getUserPreferredDetector(...)


void DrfSelect::setUserPrefferedDetector( std::shared_ptr<DetectorPeakResponse> drf,
                                     std::shared_ptr<DataBaseUtils::DbSession> sql,
                                     Wt::Dbo::ptr<InterSpecUser> user,
                                     UseDrfPref::UseDrfType prefType,
                                     std::shared_ptr<const SpecUtils::SpecFile> meas )
{
  //Check if there is a UseDrfPref for the measurement.  If so, we will delete
  //  it.  We will then add one to the database.
  
  if( !meas || !drf || !sql || !user )
    return;
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    string criteria;
    switch( prefType )
    {
      case UseDrfPref::UseDrfType::UseDetectorModelName:
        if( meas->detector_type() == DetectorType::Unknown )
          criteria = meas->instrument_model();
        else
          criteria = detectorTypeToString( meas->detector_type() );
        break;
        
      case UseDrfPref::UseDrfType::UseDetectorSerialNumber:
        criteria = meas->instrument_id();
        break;
    }//switch( prefType )
    
    //Remove preferences (there should only be a max of one) from the DB
    auto prevpref = user->drfPrefs().find().where( "MatchField = ? AND Criteria = ?" )
                         .bind( prefType )
                         .bind( criteria )
                         .resultList();
    for( auto iter = prevpref.begin(); iter != prevpref.end(); ++iter )
      iter->remove();
    
    
    //Note: adding the preference to the database, even if criteria is empty
    uint64_t uhash = drf->hashValue();
    int64_t ihash = reinterpret_cast<int64_t&>(uhash);
    auto drflist = sql->session()->find<DetectorPeakResponse>()
                       .where("InterSpecUser_id = ? AND Hash = ?")
                       .bind( user.id() )
                       .bind( ihash )
                       .resultList();
    if( drflist.size() < 1 )
      throw runtime_error( "Could not find detector with hash '" + std::to_string(uhash) + " in DB." );
    
    UseDrfPref *newprefraw = new UseDrfPref();
    Wt::Dbo::ptr<UseDrfPref> newpref( newprefraw );
      
    newprefraw->m_user = user;
    newprefraw->m_field = prefType;
    newprefraw->m_criteria = criteria;
    newprefraw->m_flags = 0x00;
    newprefraw->m_drfIndex = drflist.front().id();
    
    newpref = sql->session()->add( newpref );
    
    auto newid = newpref.id();
    auto drfindex = newpref->m_drfIndex;
    
    transaction.commit();
    
    cout << "Added preferred detector to DB for criteria='" << criteria
         << "', DRF index=" << drfindex << ", as UseDrfPref.id="
         << newid << endl;
  }catch( std::exception &e )
  {
    cerr << "setUserPrefferedDetector: Got exception setting user preferred detector to DB: " << e.what() << endl;
  }//try catch see if the user has a preference
  
}//void setUserPrefferedDetector(...)


void DrfSelect::acceptAndFinish()
{
  updateLastUsedTimeOrAddToDb( m_detector, m_interspec->user().id(), m_sql );
  
  emitChangedSignal();
  //emitModifiedSignal();
  
  auto meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  auto sql = m_interspec->sql();
  const Dbo::ptr<InterSpecUser> &user = m_interspec->user();
  auto drf = m_detector;
  
  if( meas && drf && m_defaultForSerialNumber && m_defaultForSerialNumber->isChecked() )
  {
    UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorSerialNumber;
    WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
      DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, meas );
    } ) );
  }//if( m_defaultForSerialNumber and is checked )
  
  if( meas && drf && m_defaultForDetectorModel && m_defaultForDetectorModel->isChecked() )
  {
    UseDrfPref::UseDrfType preftype = UseDrfPref::UseDrfType::UseDetectorModelName;
    WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
      DrfSelect::setUserPrefferedDetector( drf, sql, user, preftype, meas );
    } ) );
  }//if( m_defaultForDetectorModel and is checked )
  
  
  done().emit();
}//void acceptAndFinish()


void DrfSelect::cancelAndFinish()
{
  m_detector = m_originalDetector;
  if( m_detector && m_originalDetectorCopy )
    (*m_detector) = (*m_originalDetectorCopy);

  emitChangedSignal();
  //emitModifiedSignal();
  done().emit();
}//void cancelAndFinish()


void DrfSelect::finishWithNoDetector()
{
  m_detector.reset();
  emitChangedSignal();
  //emitModifiedSignal();
  done().emit();
}//void finishWithNoDetector()


vector<pair<string,string>> DrfSelect::avaliableGadrasDetectors() const
{
  vector<pair<string,string>> answer;

  //Look through paths specified in user prefernce 'GadrasDRFPath' and find
  //  directoires containing both a 'Detector.dat' and 'Efficiency.csv' file.
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::staticDataDirectory();
  const string drfpaths = SpecUtils::append_path( datadir, "GenericGadrasDetectors" )
                          + ";" + SpecUtils::append_path( datadir, "OUO_GadrasDetectors" );
#else
  const string drfpaths = UserPreferences::preferenceValue<string>( "GadrasDRFPath", m_interspec );
#endif
  
  vector<string> paths;
  SpecUtils::split( paths, drfpaths, "\r\n;" );
  
  for( std::string path : paths )
  {
    SpecUtils::trim( path );
    const vector<string> csvs = SpecUtils::recursive_ls( path, "Efficiency.csv" ); //ending string is case-insensitive
    
    for( const std::string &csv : csvs )
    {
      const string parent = SpecUtils::parent_path(csv);
      vector<string> files = SpecUtils::ls_files_in_directory(parent);
      bool found_dat = false;
      for( size_t i = 0; !found_dat && i < files.size(); ++i )
        found_dat = SpecUtils::iends_with( files[i], "Detector.dat" );
      if( found_dat )
      {
        try
        {
          //fs_relative probably shouldnt throw, but JIC.
          //fs_relative( const std::string &from_path, const std::string &to_path )
          string displayname = SpecUtils::fs_relative( path, parent );
          if( displayname.empty() || displayname == "." )
            displayname = SpecUtils::filename(parent);
          
          answer.push_back( make_pair(parent, displayname) );
        }catch(...)
        {
        }
      }//if( found_dat )
    }
  }//for( const std::string &path : paths )

  return answer;
}//vector<string> avaliableGadrasDetectors() const


shared_ptr<DetectorPeakResponse> DrfSelect::initARelEffDetector( const SpecUtils::DetectorType type,
                                                                std::string manufacturer,
                                                                std::string model )
{

  // Get a list of potential model names for a given DetectorType.
  //  Note that depending on the file format, capitalization, or non alpha-numeric character
  //  seem to commonly not be consistent (e.g., using a space instead of dash, etc).
  // TODO: Need to go through data from the various model detectors and fill in the rest of the
  //       potential names - and then maybe make this below lambda a stand-alone function
  auto potentialDetModelNames = []( const SpecUtils::DetectorType type ) -> vector<string> {
    switch( type )
    {
      case SpecUtils::DetectorType::Exploranium:
        return { "GR135" };
        
      case SpecUtils::DetectorType::IdentiFinder:
        return { "IdentiFINDER" };
        
      case SpecUtils::DetectorType::IdentiFinderNG:
        return { "IdentiFINDER-NG", "IdentiFINDER 2 NG" };
        
      case SpecUtils::DetectorType::IdentiFinderLaBr3:
        return { "IdentiFINDER-LaBr", "IdentiFINDER-LG" };
        
      case SpecUtils::DetectorType::IdentiFinderR425NaI:
        return { "IdentiFINDER-R425" };
      
      case SpecUtils::DetectorType::IdentiFinderR425LaBr:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::IdentiFinderR500NaI:
        return { "IdentiFINDER-R500-NaI" };
        
      case SpecUtils::DetectorType::IdentiFinderTungsten:
      case SpecUtils::DetectorType::IdentiFinderR500LaBr:
      case SpecUtils::DetectorType::IdentiFinderUnknown:
        break;
        
      case SpecUtils::DetectorType::DetectiveUnknown:
        return { "Detective" };
        
      case SpecUtils::DetectorType::DetectiveEx:
        return { "Detective-DX" };
        
      case SpecUtils::DetectorType::MicroDetective:
        return { "MicroDetective", "uDetective" };
        
      case SpecUtils::DetectorType::DetectiveEx100:
        return { "Detective-EX100", "Detective-DX100", "EX100", "DX100" };
        
      case SpecUtils::DetectorType::DetectiveEx200:
        return { "Detective-EX200", "Detective-DX200" };
        
      case SpecUtils::DetectorType::DetectiveX:
        return { "Detective-X" };
        
      case SpecUtils::DetectorType::SAIC8:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::Falcon5000:
        return { "Falcon5000", "Falcon5k" };
        
      case SpecUtils::DetectorType::Unknown:
        break;
        
      case SpecUtils::DetectorType::MicroRaider:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::Interceptor:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::RadHunterNaI:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::RadHunterLaBr3:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::Rsi701:
        return { "Rsi701" };
        
      case SpecUtils::DetectorType::Rsi705:
        return { "Rsi705" };
        
      case SpecUtils::DetectorType::AvidRsi:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::OrtecRadEagleNai:
        return { "RadEagle", "ORTEC RadEagle" };
        break;
        
      case SpecUtils::DetectorType::OrtecRadEagleCeBr2Inch:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::OrtecRadEagleCeBr3Inch:
        //return  { }
        break;
        
      case SpecUtils::DetectorType::OrtecRadEagleLaBr:
        return { "RadEagleLaBr" };
        break;
        
      case SpecUtils::DetectorType::Sam940LaBr3:
        //return { };
        break;
        
      case SpecUtils::DetectorType::Sam940:
        return { "Sam-940" };
        
      case SpecUtils::DetectorType::Sam945:
        return { "Sam-945" };
        
      case SpecUtils::DetectorType::Sam950:
        return { "BNC SAM-950 3x3 NaI", "SAM-950" };
        
      case SpecUtils::DetectorType::Srpm210:
        //return { }
        break;
        
      case SpecUtils::DetectorType::RIIDEyeNaI:
        //return { }
        break;
        
      case SpecUtils::DetectorType::RIIDEyeLaBr:
        //return { }
        break;
        
      case SpecUtils::DetectorType::RadSeekerNaI:
        return { "RadSeeker-DL" };
        
      case SpecUtils::DetectorType::RadSeekerLaBr:
        return { "Radseeker-LaBr3" };
        break;
        
      case SpecUtils::DetectorType::VerifinderNaI:
        return { "Verifinder-SN23N", "Verifinder-NaI" };
        break;
        
      case SpecUtils::DetectorType::VerifinderLaBr:
        //return { }
        break;
        
      case SpecUtils::DetectorType::KromekD3S:
        return { "D3S", "Kromek D3S", "D3", "Kromek D3S", "Kromek D3" };
      
      case SpecUtils::DetectorType::RadiaCode:
        return { "Radiacode 102", "Radiacode", "Radiacode 101" };
        break;
        
      case SpecUtils::DetectorType::Fulcrum:
        //return {};
        break;
        
      case SpecUtils::DetectorType::Fulcrum40h:
        return { "Fulcrum40h" };
        
    }//switch( type )
    
    return {};
  };//potentialDetModelNames lambda
  
  // Remove all non-alphanumeric
  auto removeNonAlphaNumAndLower = []( string &str ){
    str.erase( std::remove_if(begin(str), end(str), [](char a){return !isalnum(a);}), end(str) );
    SpecUtils::to_lower_ascii( str );
  };
  
  
  vector<string> potential_names;
  
  for( string &s : potentialDetModelNames(type) )
  {
    removeNonAlphaNumAndLower( s );
    if( s.size() > 2 )
      potential_names.push_back( s );
  }//for( string &s : potentialDetModelNames(type) )
  
  removeNonAlphaNumAndLower( model );
  
  // If we have DetectorType, we will trust that more than the string model number, because, for
  //  example some IndentiFINDER-NG will list the model as "IndentiFINDER" in the N42 file, while
  //  SpecUtils, goes into the comments and uses detector size to infer detector model in this case.
  if( (type == SpecUtils::DetectorType::Unknown) && (model.size() > 2) )
    potential_names.push_back( model );
  
  removeNonAlphaNumAndLower( manufacturer );
  
  if( potential_names.empty() )
    throw std::runtime_error( "Could not find Rel. Eff. detector response functions" );
  
  
  vector<string> user_paths, default_paths;
  
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  try
  {
    const string userDir = InterSpec::writableDataDirectory();
    potential_rel_eff_files( userDir, user_paths );
  }catch( std::exception & )
  {
    cerr << "Couldnt call into InterSpec::writableDataDirectory()" << endl;
  }
#endif
  
  try
  {
    const string staticDataDir = InterSpec::staticDataDirectory();
    potential_rel_eff_files( staticDataDir, default_paths );
    
    // TODO: we could prioritize "common_drfs.tsz"
  }catch( std::exception & )
  {
    cerr << "Couldnt call into InterSpec::staticDataDirectory()" << endl;
  }
  
  vector<string> paths = user_paths;
  paths.insert( end(paths), begin(default_paths), end(default_paths) );
  
  for( const string &filename : paths )
  {
#ifdef _WIN32
    const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
    std::ifstream input( wfilename.c_str() );
#else
    std::ifstream input( filename.c_str() );
#endif
    
    if( !input.is_open() )
      continue;
  
    string line;
    while( SpecUtils::safe_get_line( input, line, 2048 ) )
    {
      SpecUtils::trim( line );
      if( line.empty() || line[0]=='#' )
        continue;
      
      auto det = DetectorPeakResponse::parseSingleCsvLineRelEffDrf( line );
      string det_name = det ? det->name() : string();

      removeNonAlphaNumAndLower( det_name );
      
      if( !det || !det->isValid() || det_name.empty() )
        continue;
      
      bool isMatch = false;
      
      for( const string &candidate_name : potential_names )
      {
        const size_t pos = det_name.find( candidate_name );
        if( pos == string::npos )
          continue;
        
        // For Identifinder (and probably RadEagle, RadSeeker, etc), we want the model to go all the
        //  way to the end of the DRF name - since we dont want to match like "IdentiFINDER-NGH",
        //  when the detector is "just" an IdentiFINDER.
        switch( type )
        {
          case SpecUtils::DetectorType::IdentiFinder:
          {
            if( (pos + candidate_name.size()) < det_name.size() )
              continue;
            
            break;
          }//case (we should )
            
          default:
            break;
        }//switch( type )
        
        isMatch = true;
        break;
      }//for( const string &candidate_name : potential_names )
      
      
      if( isMatch )
      {
        const auto filenampos = std::find(begin(default_paths), end(default_paths), filename );
        
        if( filenampos == end(default_paths) )
          det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf );
        else
          det->setDrfSource( DetectorPeakResponse::DrfSource::DefaultRelativeEfficiencyDrf );
        
        return det;
      }
    }//while( SpecUtils::safe_get_line( input, line ) )
  }//for( const string &filename : paths )
  
  throw runtime_error( "Coudlnt find detector in relative efficiency files" );
  
  return std::shared_ptr<DetectorPeakResponse>();
}//initARelEffDetector( int type )


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetector(
                                                            const SpecUtils::DetectorType type, InterSpec *interspec )

{
  using SpecUtils::DetectorType;
  
  string name;
  switch( type )
  {
    case DetectorType::Exploranium:       name = "GR135";              break;
    case DetectorType::IdentiFinder:      name = "identiFINDER-N";     break;
    case DetectorType::IdentiFinderNG:    name = "identiFINDER-NGH";   break;
    case DetectorType::IdentiFinderLaBr3: name = "identiFINDER-LaBr3"; break;
    case DetectorType::DetectiveUnknown:  name = "Detective";          break;
    case DetectorType::DetectiveEx:       name = "Detective-EX";       break;
    case DetectorType::DetectiveEx100:    name = "Detective-EX100";    break;
    case DetectorType::DetectiveEx200:    name = "Ortec IDM Portal";   break;
    case DetectorType::DetectiveX:        name = "Detective X";        break;
    case DetectorType::SAIC8:             name = "";                   break;
    case DetectorType::Falcon5000:        name = "Falcon 5000";        break;
    case DetectorType::Unknown:           name = "";                   break;
    case DetectorType::MicroDetective:    name = "MicroDetective";     break;
      
    //case DetectorType::RadHunterNaI:            name = ""; break;
    //case DetectorType::RadHunterLaBr3:          name = ""; break;
    case DetectorType::Rsi701:                    name = "NaI 2x4x16"; break;
    case DetectorType::Rsi705:                    name = "NaI 2x4x16"; break;
    case DetectorType::AvidRsi:                   name = "NaI 2x4x16"; break;
    case DetectorType::OrtecRadEagleNai:          name = "RadEagle"; break;
    //case DetectorType::OrtecRadEagleCeBr2Inch: name = ""; break;
    //case DetectorType::OrtecRadEagleCeBr3Inch: name = ""; break;
    //case DetectorType::OrtecRadEagleLaBr: name = ""; break;
    //case DetectorType::Sam940LaBr3: name = ""; break;
    case DetectorType::Sam940:                     name = "SAM-945"; break;
    case DetectorType::Sam945:                     name = "SAM-945"; break;
    case DetectorType::Srpm210:                    name = "SRPM-210"; break;
      
    case DetectorType::RIIDEyeNaI:                 name = "RIIDEyeX-GN1"; break;
    case DetectorType::Interceptor:                name = "Interceptor"; break;
      
    case DetectorType::MicroRaider:                name = "Raider"; break;
      
    case DetectorType::RadSeekerNaI:               name = "RadSeeker-NaI"; break;
    case DetectorType::RadSeekerLaBr:              name = "Radseeker-LaBr3"; break;
    
    case DetectorType::VerifinderNaI:              name = "Verifinder-NaI"; break;
    case DetectorType::VerifinderLaBr:             name = "Verifinder-LaBr3"; break;
      
    case DetectorType::IdentiFinderR500NaI:        name = "IdentiFINDER-R500-NaI";   break;
    case DetectorType::IdentiFinderR500LaBr:       name = "IdentiFINDER-R500-LaBr3"; break;
      
    case SpecUtils::DetectorType::KromekD3S:       name = "D3S"; break;
    case DetectorType::Fulcrum40h:                 name = "Fulcrum40h"; break;
    case DetectorType::IdentiFinderR425NaI: 	     name = "IdentiFINDER-R425"; break;
    case DetectorType::Sam950:                     name = "SAM-950GN-N30"; break;
    case DetectorType::Sam940LaBr3:                name = "SAM-Eagle-LaBr3"; break;
      
    case DetectorType::OrtecRadEagleCeBr2Inch:
    case DetectorType::OrtecRadEagleCeBr3Inch:
    case DetectorType::OrtecRadEagleLaBr:
    case DetectorType::RadHunterLaBr3:
    case DetectorType::RadHunterNaI:
    case DetectorType::IdentiFinderTungsten:
    case DetectorType::IdentiFinderUnknown:
    case DetectorType::RIIDEyeLaBr:
    case DetectorType::Fulcrum:
    case DetectorType::IdentiFinderR425LaBr:
    case SpecUtils::DetectorType::RadiaCode:
      // \TODO: fill in these names
      break;
  }//switch( type )
  
  if( name.empty() )
    throw runtime_error( "There is no GADRAS detector response function for a "
                        + SpecUtils::detectorTypeToString(type) );

  std::shared_ptr<DetectorPeakResponse> det;
  
  try
  {
    det = initAGadrasDetector( name, interspec );
    if( det )
      return det;
  }catch(...)
  {
  }
  
  //We couldnt find an exact match - lets be a little looser
  string dettype, deteff;
  switch( type )
  {
    case DetectorType::Exploranium:
      break;
    case DetectorType::IdentiFinder:
      break;
    case DetectorType::IdentiFinderNG:
      break;
    case DetectorType::OrtecRadEagleNai:
      break;
      
    case DetectorType::IdentiFinderLaBr3:
      break;
      
    case DetectorType::DetectiveUnknown:
    case DetectorType::DetectiveEx:
    case DetectorType::MicroDetective:
      break;
    case DetectorType::DetectiveEx100:
      break;
    case DetectorType::Falcon5000:
      break;
  
    default:
      break;
  }//switch( type )
  
  const string secondname = dettype + " " + deteff;
  
  try
  {
    det = initAGadrasDetector( secondname, interspec );
    if( det )
      return det;
  }catch(...)
  {
  }
  
  throw runtime_error( "Could not find GADRAS detector names " + name + " or " + secondname  );
  return det;
}//initAGadrasDetector


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetector( const std::string &currentDetName, InterSpec *interspec )
{
  //Grab "GadrasDRFPath" and split it by semicolon and newlines, then go through
  //  and look for sub-folders with the given name that contain Detector.data
  //  and Efficiency.csv.
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::staticDataDirectory();
  const string drfpaths = SpecUtils::append_path( datadir, "GenericGadrasDetectors" )
                          + ";" + SpecUtils::append_path( datadir, "OUO_GadrasDetectors" );
#else
  const string drfpaths = UserPreferences::preferenceValue<string>( "GadrasDRFPath", interspec );
#endif
  
  
  vector<string> paths;
  SpecUtils::split( paths, drfpaths, "\r\n;" );
  for( string basepath : paths )
  {
    SpecUtils::trim( basepath );
    const string path = SpecUtils::append_path(basepath, currentDetName);
    if( !SpecUtils::is_directory(path) )
      continue;
    
    string thiscsv, thisdat;
    const vector<string> files = SpecUtils::ls_files_in_directory( path );
    for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
    {
      //We will be loose with upper/lower case
      if( SpecUtils::iends_with(files[i], "Efficiency.csv") )
        thiscsv = files[i];
      else if( SpecUtils::iends_with(files[i], "Detector.dat") )
        thisdat = files[i];
    }//for( size_t i = 0; i < files.size(); ++i )
    
    if( thiscsv.size() && thisdat.size() )
    {
      auto det = initAGadrasDetectorFromDirectory( path );
      if( det )
      {
        if( SpecUtils::icontains( basepath, "GenericGadrasDetectors") )
          det->setDrfSource( DetectorPeakResponse::DrfSource::DefaultGadrasDrf );
        else
          det->setDrfSource( DetectorPeakResponse::DrfSource::UserAddedGadrasDrf );
      }
      
      return det;
    }
  }//for( const string &basepath : paths )
  
  throw runtime_error( "Could not find GADRAS detector names " + currentDetName );
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetector( const std::string &name );


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetectorFromDirectory( const std::string &path )
{
  auto det = std::make_shared<DetectorPeakResponse>();

  if( !SpecUtils::is_directory(path) )
    throw runtime_error( "'" + path + "' is not a directory" );
  
  //Lets be lienient on upper vs lower case.
  string thiscsv, thisdat;
  const vector<string> files = SpecUtils::ls_files_in_directory( path );
  for( size_t i = 0; (thiscsv.empty() || thisdat.empty()) && i < files.size(); ++i )
  {
    //We will be loose with upper/lower case
    if( SpecUtils::iends_with(files[i], "Efficiency.csv") )
      thiscsv = files[i];
    else if( SpecUtils::iends_with(files[i], "Detector.dat") )
      thisdat = files[i];
  }//for( size_t i = 0; i < files.size(); ++i )
  
  if( thiscsv.empty() )
    throw runtime_error( "Could not find a Efficiency.csv file in '" + path + "'" );
  
  if( thisdat.empty() )
    throw runtime_error( "Could not find a Detector.dat file in '" + path + "'" );
  
#ifdef _WIN32
  const std::wstring wthiscsv = SpecUtils::convert_from_utf8_to_utf16(thiscsv);
  ifstream effstrm( wthiscsv.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream effstrm( thiscsv.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !effstrm.is_open() )
    throw runtime_error( "Could not open '" + thiscsv + "'" );

#ifdef _WIN32
  const std::wstring wthisdat = SpecUtils::convert_from_utf8_to_utf16(thisdat);
  ifstream datstrm( wthisdat.c_str(), ios_base::binary|ios_base::in );
#else
  ifstream datstrm( thisdat.c_str(), ios_base::binary|ios_base::in );
#endif
  
  if( !datstrm.is_open() )
    throw runtime_error( "Could not open '" + thisdat + "'" );

  det->fromGadrasDefinition( effstrm, datstrm );
  det->setName( SpecUtils::filename( SpecUtils::parent_path(thiscsv) ) );

//#if( PERFORM_DEVELOPER_CHECKS )
//  check_url_serialization( det );
//#endif

  return det;
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetectorFromDirectory()


void DrfSelect::emitChangedSignal()
{
  //Make sure this is a necessary signal to emit
  if( m_previousEmmittedDetector == m_detector )
    return;
  
  if( !m_detector && !m_previousDetectorDef.isValid() )
    return;
  
  if( m_detector && ((*m_detector)==m_previousDetectorDef) )
    return;

  const auto prev = m_previousEmmittedDetector;
  const auto current = m_detector;
  
  //Now update m_previousDetectorDef and m_previousEmmittedDetector so we can
  //  appropriately filter out unecessary signal emits in the future
  if( m_detector )
    m_previousDetectorDef = *m_detector;
  else
    m_previousDetectorDef.reset();
  m_previousEmmittedDetector = m_detector;

  m_interspec->detectorChanged().emit( m_detector );
  setAcceptButtonEnabled( !!m_detector );
  
  updateChart();
  
  
  auto undo = [prev](){
    InterSpec *viewer = InterSpec::instance();
    
    // This widget *should* (but as of 20230428, this looks broke) hear the
    //  detectorChanged() signal, and update its state based on this; see
    //  #DrfSelect::setDetector.
    if( viewer )
      viewer->detectorChanged().emit( prev );
  };
  
  auto redo = [current](){
    InterSpec *viewer = InterSpec::instance();
    if( viewer )
      viewer->detectorChanged().emit( current );
  };
  
  UndoRedoManager *undoManager = m_interspec->undoRedoManager();
  if( undoManager )
    undoManager->addUndoRedoStep( undo, redo, "Change DRF" );
}//void emitChangedSignal()


/*
void DrfSelect::emitModifiedSignal()
{
  //Make sure this is a necessary signal to emit
  if( !m_detector && !m_previousDetectorDef.isValid() )
    return;

  if( m_detector && ((*m_detector)==m_previousDetectorDef) )
    return;

  //Now update m_previousDetectorDef and m_previousEmmittedDetector so we can
  //  appropriately filter out unecessary signal emits in the future
  if( m_detector )
    m_previousDetectorDef = *m_detector;
  else
    m_previousDetectorDef.reset();
  m_previousEmmittedDetector = m_detector;

  m_interspec->detectorChanged().emit( m_detector );
}//void emitModifiedSignal()
*/

Wt::Signal<> &DrfSelect::done()
{
  return m_doneSignal;
}


DrfSelectWindow::DrfSelectWindow( InterSpec *viewer )
  : AuxWindow( WString::tr("window-title-select-det-eff-fcn"),
              Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                   | AuxWindowProperties::TabletNotFullScreen
                   | AuxWindowProperties::DisableCollapse
                   | AuxWindowProperties::EnableResize
                   | AuxWindowProperties::SetCloseable ),
    m_edit( NULL ),
    m_interspec( viewer )
{
  assert( m_interspec );
  if( !m_interspec )
    m_interspec = InterSpec::instance(); //jic
  if( !m_interspec )
    throw runtime_error( "DrfSelectWindow: invalid ptr to InterSpec" );
  
  shared_ptr<SpecMeas> meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<DetectorPeakResponse> det = meas ? meas->detector() : nullptr;
  SpecMeasManager *fileManager = m_interspec->fileManager();
  SpectraFileModel *fileModel = fileManager ? fileManager->model() : nullptr;
  
  m_edit = new DrfSelect( det, m_interspec, fileModel, this );
  
  contents()->setPadding( 0 );
  contents()->setOverflow( WContainerWidget::Overflow::OverflowHidden );
  stretcher()->addWidget( m_edit, 0, 0 );
  stretcher()->setContentsMargins( 6, 2, 6, 0 );
  
  m_edit->done().connect( boost::bind( &DrfSelectWindow::acceptAndDelete, this ) );
  finished().connect( boost::bind( &DrfSelectWindow::acceptAndDelete, this ) );
  rejectWhenEscapePressed();
  
  show();
  
  if( !m_isPhone )
  {
    // None this would have effect for phone
    resize( WLength(650,WLength::Pixel), WLength(610,WLength::Pixel));
    centerWindow();
    
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }//
}//DrfSelectWindow constructor


DrfSelectWindow::~DrfSelectWindow()
{
}


DrfSelect *DrfSelectWindow::widget()
{
  return m_edit;
}

void DrfSelectWindow::acceptAndDelete( DrfSelectWindow *window )
{
  window->m_edit->emitChangedSignal();  //will only emit if necessary
  //window->m_edit->emitModifiedSignal(); //will only emit if necessary
  
  if( window->m_interspec )
    window->m_interspec->closeDrfSelectWindow();
  else
    AuxWindow::deleteAuxWindow( window );
}//void acceptAndDelete( DrfSelectWindow *window )
