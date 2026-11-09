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
#include <set>
#include <array>
#include <mutex>
#include <filesystem>
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

#include <Wt/WMenu.h>
#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WBreak.h>
#include <Wt/WLabel.h>
#include <Wt/WTimer.h>
#include <Wt/WTable.h>
#include <Wt/WSignal.h>
#include <Wt/WServer.h>
#include <Wt/WRandom.h>
#include <Wt/WDateTime.h>
#include <Wt/WTextArea.h>
#include <Wt/WResource.h>
#include <Wt/WLineEdit.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WIOService.h>
#include <Wt/WTableCell.h>
#include <Wt/WFileUpload.h>
#include <Wt/WPushButton.h>
#include <Wt/WGridLayout.h>
#include <Wt/WApplication.h>
#include <Wt/WRadioButton.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WEnvironment.h>
#include <Wt/Http/Request.h>
#include <Wt/Http/Response.h>
#include <Wt/WItemDelegate.h>
#include <Wt/WStackedWidget.h>
#include <Wt/Dbo/QueryModel.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WStandardItemModel.h>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/DrfChart.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/DrfLibrary.h"
#include "InterSpec/WidgetUtils.h"
#include "InterSpec/DrfLibraryWidgets.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/DrfAppearance.h"
#include "InterSpec/GadrasDetectorDat.h"
#include "InterSpec/MakeMcResponseForDrf.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/DrfModifyWidget.h"
#include "InterSpec/DrfImportWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

// CeeLo (external_libs/CeeLo/src): the MC-parameterized response.
#include "io/DetectorResponse.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#endif

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
#include "InterSpec/DirectorySelector.h"
#endif

#if( BUILD_AS_OSX_APP )
#include "target/macos/macOsUtils.h"
#endif

// The regex in GCC 4.8.x does not have working regex
#if( defined(__GLIBCXX__) && (__cplusplus < 201402L) )
static_assert( defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) \
              || defined(_GLIBCXX_REGEX_STATE_LIMIT) \
              || (defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE > 4), "GCC 4.8 is not supported due to buggy regex implementation" );
#endif



using namespace std;
using namespace Wt;


using SpecUtils::Measurement;
using SpecUtils::DetectorType;

namespace
{
  
  
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
    string m_source_to_face_distance;    // Empty if not user specified; from the detector face
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
      
      const bool have_face_dist = SpecUtils::parse_double( add_info[2].c_str(), add_info[2].size(),
                                                           src_to_det_face );
      
      // 3 - source radius
      
      if( SpecUtils::parse_double(add_info[4].c_str(), add_info[4].size(), det_setback)
         && (det_setback > 0.0) )
        results.m_detector_setback = add_info[4] + " cm";
      
      // 5 - Detector Area
      
      const bool have_crystal_dist = SpecUtils::parse_double( add_info[6].c_str(), add_info[6].size(),
                                                              src_to_crystal );
      
      // The Distance field is measured from the detector FACE, as every InterSpec distance is (the
      //  setback goes in its own field), so use the "Source to Detector Face" column; the "Source to
      //  Crystal Face" one is that plus the setback.  (Same as DetectorPeakResponse's relative
      //  efficiency CSV parser.)
      if( have_face_dist && (src_to_det_face > 0.0) )
        results.m_source_to_face_distance = add_info[2] + " cm";
      else if( have_crystal_dist && (src_to_crystal > 0.0) )
        results.m_source_to_face_distance
                   = SpecUtils::printCompact( std::max( 0.0, src_to_crystal - std::max( 0.0, det_setback ) ), 6 ) + " cm";
      
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
    
    if( SpecUtils::is_directory(drfsdir) )
    {
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
    }//if( SpecUtils::is_directory(drfsdir) )
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


}//namespace




DetectorDisplay::DetectorDisplay( InterSpec *specViewer )
  : WContainerWidget(),
    m_genericIcon( nullptr ),
    m_drfIcon( nullptr ),
    m_text( nullptr ),
    m_interspec( specViewer )
{
  assert( specViewer );

  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", specViewer );

  HelpSystem::attachToolTipOn( this, WString::tr("app-tt-det-disp"), showToolTips ); //"app-tt-det-disp" defined in InterSpec.xml localisation, so we dont need to load DrfSelect.xml for this widget

  addStyleClass( "DetectorDisplay" );  //In InterSpec.css since this widget is loaded almost always at initial load time anyway

  // The icon says "detector", so there is no "Detector:" label on screen - only for screen readers.
  //  A button for the keyboard too: focusable, and Enter / Space open the dialog.
  setAttributeValue( "role", "button" );
  setAttributeValue( "tabindex", "0" );
  keyWentDown().connect( "function(s,e){ if( (e.key === 'Enter') || (e.key === ' ') ){ e.preventDefault(); s.click(); } }" );

  m_genericIcon = addNew<WContainerWidget>();
  m_genericIcon->setInline( true );
  m_genericIcon->addStyleClass( "DrfIcon DrfIconGeneric" );
  m_drfIcon = addNew<WText>( "", Wt::TextFormat::UnsafeXHTML );  //see DrfAppearance::icon_html
  m_drfIcon->addStyleClass( "DrfIconHolder" );
  m_drfIcon->hide();

  m_text = addNew<WText>( "", Wt::TextFormat::Plain );  //the DRF's name is user text
  m_text->addStyleClass( "DetName" );

  std::shared_ptr<DetectorPeakResponse> detector;
  auto meas = specViewer->measurment(SpecUtils::SpectrumType::Foreground);
  if( meas )
    detector = meas->detector();
  setDetector( detector );

  specViewer->detectorChanged().connect( this, &DetectorDisplay::setDetector );
  specViewer->detectorModified().connect( this, &DetectorDisplay::setDetector );

  clicked().connect( this, &DetectorDisplay::editDetector );
}//DetectorDisplay constructor


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
  const bool have_det = (det && det->isValid());

  const shared_ptr<const DrfAppearance> look = have_det ? det->appearance() : nullptr;
  const string icon_html = look ? look->icon_html( 16.0, 36.0 ) : string();
  m_drfIcon->setText( WString::fromUTF8( icon_html ) );
  m_drfIcon->setHidden( icon_html.empty() );
  m_genericIcon->setHidden( !icon_html.empty() );

  const WString name = have_det ? WString::fromUTF8( det->name() ) : WString::tr("app-det-none");
  m_text->setText( name );
  m_text->toggleStyleClass( "DetNone", !have_det );

  WString aria_label = WString::tr("detector-label");
  aria_label += " ";
  aria_label += name;
  setAttributeValue( "aria-label", aria_label );
}//void setDetector( std::shared_ptr<DetectorPeakResponse> det )


void DetectorDisplay::editDetector()
{
  m_interspec->showDrfSelectWindow();
}//void editDetector()



void DrfSelect::createChooseDrfDialog( vector<shared_ptr<DetectorPeakResponse>> inputdrfs,
                                       WString mainMsgHtml,
                                       string creditsHtml,
                                       std::function<void()> saveDrfsCallBack,
                                       DrfSelect::ChooseDrfHooks hooks )
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
  
  SimpleDialog *dialog = SimpleDialog::make( WString::tr(title_key) );
  dialog->addStyleClass( "DrfFileSelectDialog" );
  
  WPushButton *cancel = dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
  WPushButton *accept = dialog->addButton( WString::tr("Accept"), WidgetUtils::ButtonRole::Affirm );
  
  //WGridLayout *layout = new WGridLayout( dialog->contents() );
  //layout->setContentsMargins( 0, 0, 0, 0 );
  
  dialog->contents()->setOverflow( Overflow::Hidden, Wt::Orientation::Horizontal );
  
  DrfChart *chart = dialog->contents()->addNew<DrfChart>();
  chart->addStyleClass( "DrfFileSelectChart" );

  if( !mainMsgHtml.empty() )
  {
    WText *gen_desc = dialog->contents()->addNew<WText>( mainMsgHtml );
    gen_desc->addStyleClass( "DrfFileSelectMainDesc" );
    gen_desc->setInline( false );
  }//if( !descrip_html.empty() )

  // Create drfDescription now for lambda capture; will be added to dialog->contents() after combo
  auto drfDescriptionOwned = std::make_unique<WText>();
  WText *drfDescription = drfDescriptionOwned.get();
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
    WComboBox *combo = dialog->contents()->addNew<WComboBox>();
    combo->addStyleClass( "DrfFileSelectCombo" );
    combo->setInline( false );

    for( auto drf : drfs )
      combo->addItem( drf->name() );
    combo->setCurrentIndex( 0 );
    combo->activated().connect( combo, [showDrf]( int arg ){ showDrf( arg ); } );
  }//if( drfs.size() > 1 )

  drfDescription->addStyleClass( "DrfFileSelectDrfDesc" );
  drfDescription->setInline( false );
  dialog->contents()->addWidget( std::move(drfDescriptionOwned) );

  if( !creditsHtml.empty() )
  {
    WText *credits = dialog->contents()->addNew<WText>( creditsHtml );
    credits->addStyleClass( "DrfFileSelectCredits" );
  }//if( !descrip_html.empty() )

  // Caller-specific controls sit between what the DRF *is* and what to do with it.
  if( hooks.addContent )
    hooks.addContent( dialog->contents() );
  
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
    saveCbDiv = dialog->contents()->addNew<WContainerWidget>();
    saveCbDiv->addStyleClass( "DrfFileSelectSaveCbs" );
  }
  
  if( makeDrfsSaveCb )
  {
    saveCb = saveCbDiv->addNew<WCheckBox>( WString::tr("ds-save-for-later-cb") );
    saveCb->setInline( false );
  }//if( makeDrfsSaveCb )

  if( makeSerialNumCb )
  {
    const WString msg = WString::tr("ds-use-for-serialnum-cb").arg( meas->instrument_id() );
    defaultForSerialNumber = saveCbDiv->addNew<WCheckBox>( msg );
    defaultForSerialNumber->setInline( false );
  }//if( have serial number )

  if( makeModelCb )
  {
    string model;
    if( meas->detector_type() == DetectorType::Unknown )
      model = meas->instrument_model();
    else
      model = detectorTypeToString( meas->detector_type() );

    const WString msg = WString::tr("ds-use-for-model-cb").arg( model );
    defaultForDetectorModel = saveCbDiv->addNew<WCheckBox>( msg );
    defaultForDetectorModel->setInline( false );
  }//if( makeModelCb )
  
  
  accept->clicked().connect( accept, [=](){
    auto drf = *selected_drf;
    auto interspec = InterSpec::instance();
    if( !drf || !interspec )
      return;
      
    // Before the DB write and before it is applied: a hook that modifies the DRF changes its hash,
    //  which is the key the "Previous" row is stored under, so it has to run first or the stored
    //  detector is not the one the user gets.
    if( hooks.beforeAccept )
      hooks.beforeAccept( drf );
    
    auto sql = interspec->sql();
    const Dbo::ptr<InterSpecUser> &user = interspec->user();
    DrfSelect::updateLastUsedTimeOrAddToDb( drf, user.id(), sql );

    const shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<DetectorPeakResponse> prev = meas ? meas->detector() : nullptr;
    interspec->detectorChanged().emit( drf ); //This loads it to the foreground spectrum file

    // Through InterSpec::instance(), as for "Use" in the Detector Efficiency dialog
    UndoRedoManager * const undoManager = interspec->undoRedoManager();
    if( (prev != drf) && undoManager && undoManager->canAddUndoRedoNow() )
    {
      undoManager->addUndoRedoStep( [prev](){
        InterSpec * const viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( prev );
      }, [drf](){
        InterSpec * const viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( drf );
      }, "Change DRF" );
    }//if( undoManager )
      
    if( saveCb && saveCb->isChecked() && !!saveDrfsCallBack )
      saveDrfsCallBack();

    // On this (the session) thread: copying the `Dbo::ptr` user is not thread safe.
    for( const pair<WCheckBox *,UseDrfPref::UseDrfType> &cb
          : { make_pair( defaultForSerialNumber, UseDrfPref::UseDrfType::UseDetectorSerialNumber ),
              make_pair( defaultForDetectorModel, UseDrfPref::UseDrfType::UseDetectorModelName ) } )
    {
      if( meas && cb.first && cb.first->isChecked() )
        DrfSelect::setUserPrefferedDetector( drf, sql, user, cb.second, meas );
    }//for( both "use as default" checkboxes )
    
    // Last, so any follow-up the caller offers is on top of a detector already in use and saved.
    if( hooks.afterAccept )
      hooks.afterAccept( drf );
  } );
  
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


namespace
{
/** The detector-modeling choices offered when an app-URL brings in a DRF that knows its physical
 shape.  A URL can never carry a Monte-Carlo response (~17 KB deflated, far past a QR code), so the
 shape arrives inert - these are the ways to make it count.  Ordered cheapest first; the index is
 the combo index.
 */
enum class UrlDrfModeling : int
{
  /** Leave the DRF as the classic far-field curve, accurate only at the reference location it was
   characterized at.  The geometry is still stored with the detector.
   */
  FlatDisk = 0,
  
  /** EFFTRAN-style transfer through the geometry, anchored on the DRF's own curve.  Instant and
   deterministic, and what every other geometry-bearing import does - see the
   CeeLoUtils::attachCurveTransferResponse call sites.
   */
  GeometryTransfer = 1,
  
  /** Monte Carlo of only the on-axis energy backbone, with the ray-traced kernel carrying the
   distance/angle transfer - roughly a tenth the work of a full characterization.
   */
  PartialMc = 2,
  
  /** A full characterization over the general profile. */
  StandardMc = 3,
  
  /** A full characterization with the extra near-field detail of the contact profile. */
  ExtendedMc = 4
};//enum class UrlDrfModeling


/** Whether a choice needs a Monte-Carlo run (and so the accuracy control, and the hand-off to the
 tool that can show progress) rather than being instant.
 */
bool url_modeling_needs_mc( const UrlDrfModeling choice )
{
  return (choice == UrlDrfModeling::PartialMc)
         || (choice == UrlDrfModeling::StandardMc)
         || (choice == UrlDrfModeling::ExtendedMc);
}//url_modeling_needs_mc(...)
}//namespace


void DrfSelect::handle_app_url_drf( const std::string &url_query )
{
  try
  {
    auto drf = make_shared<DetectorPeakResponse>();
    drf->fromAppUrl( url_query );
    
    assert( drf->isValid() );
    
    InterSpec *interspec = InterSpec::instance();
    if( interspec )
      interspec->useMessageResourceBundle( "DrfSelect" );
    
    // A URL never carries a Monte-Carlo response, so a DRF that arrives knowing its shape has a
    //  choice to make about what that shape is used for.  Ask it here, where the user already is,
    //  rather than leaving the detector a flat disk and hoping they find the Modify dialog.
    const bool geom_only = (drf->geometry() && !drf->ceeloResponse() && !drf->isFixedGeometry());
    
    const WString msg = WString::tr( geom_only ? "ds-url-drf-geom-only" : "ds-url-drf-received" );
    
    ChooseDrfHooks hooks;
    
    if( geom_only )
    {
      // Shared by the hooks below: the combos live in the dialog, so they are reached through
      //  observing_ptrs that null themselves if the dialog goes first.
      auto modeling = make_shared<Wt::Core::observing_ptr<WComboBox>>();
      auto accuracy = make_shared<Wt::Core::observing_ptr<WComboBox>>();
      
      hooks.addContent = [modeling,accuracy]( WContainerWidget *parent ){
        WContainerWidget *row = parent->addNew<WContainerWidget>();
        row->addStyleClass( "DrfUrlModelingRow" );
        
        WLabel *label = row->addNew<WLabel>( WString::tr("ds-url-modeling-label") );
        WComboBox *combo = row->addNew<WComboBox>();
        label->setBuddy( combo );
        
        combo->addItem( WString::tr("ds-url-modeling-flat") );
        combo->addItem( WString::tr("ds-url-modeling-transfer") );
        combo->addItem( WString::tr("ds-url-modeling-partial-mc") );
        combo->addItem( WString::tr("ds-url-modeling-standard-mc") );
        combo->addItem( WString::tr("ds-url-modeling-extended-mc") );
        combo->setCurrentIndex( static_cast<int>(UrlDrfModeling::GeometryTransfer) );
        
        // Its own row: it qualifies the choice above rather than sitting beside it, and only
        //  appears once a Monte Carlo is actually going to run.
        WContainerWidget *accRow = parent->addNew<WContainerWidget>();
        accRow->addStyleClass( "DrfUrlModelingRow" );
        
        WLabel *accLabel = accRow->addNew<WLabel>( WString::tr("ds-url-accuracy-label") );
        WComboBox *acc = accRow->addNew<WComboBox>();
        accLabel->setBuddy( acc );
        acc->addItem( WString::tr("ds-url-accuracy-fast") );
        acc->addItem( WString::tr("ds-url-accuracy-normal") );
        acc->addItem( WString::tr("ds-url-accuracy-thorough") );
        acc->setCurrentIndex( 1 );  //Normal
        accRow->hide();
        
        // Most people meeting this dialog will not know what any of the options mean, so say what
        //  the selected one does, right here.
        WText *desc = parent->addNew<WText>( WString::tr("ds-url-modeling-desc-transfer") );
        desc->addStyleClass( "DrfUrlModelingDesc" );
        desc->setInline( false );
        
        combo->changed().connect( combo, [combo,desc,accRow](){
          const UrlDrfModeling choice = static_cast<UrlDrfModeling>( combo->currentIndex() );
          
          const char *key = "ds-url-modeling-desc-transfer";
          switch( choice )
          {
            case UrlDrfModeling::FlatDisk:         key = "ds-url-modeling-desc-flat";        break;
            case UrlDrfModeling::GeometryTransfer: key = "ds-url-modeling-desc-transfer";    break;
            case UrlDrfModeling::PartialMc:        key = "ds-url-modeling-desc-partial-mc";  break;
            case UrlDrfModeling::StandardMc:       key = "ds-url-modeling-desc-standard-mc"; break;
            case UrlDrfModeling::ExtendedMc:       key = "ds-url-modeling-desc-extended-mc"; break;
          }
          desc->setText( WString::tr(key) );
          
          accRow->setHidden( !url_modeling_needs_mc(choice) );
        } );
        
        *modeling = combo;
        *accuracy = acc;
      };//hooks.addContent
      
      // Runs before the DB write, so what gets stored is the detector the user chose - attaching a
      //  response changes the hash the "Previous" row is keyed on.
      hooks.beforeAccept = [modeling]( shared_ptr<DetectorPeakResponse> accepted ){
        if( !accepted )
          return;
        
        // Default to the transfer if the combo is somehow gone: it is what every other
        //  geometry-bearing import gives, and costs nothing.
        const UrlDrfModeling choice = (*modeling)
                  ? static_cast<UrlDrfModeling>( (*modeling)->currentIndex() )
                  : UrlDrfModeling::GeometryTransfer;
        
        // A Monte-Carlo run happens after the detector is in use (see afterAccept); it needs a seed
        //  with no response attached, which is exactly what we have here.
        if( choice == UrlDrfModeling::GeometryTransfer )
          CeeLoUtils::attachCurveTransferResponse( *accepted );

        // "Saved with the detector, but not used" (ds-url-modeling-desc-flat).
        if( choice == UrlDrfModeling::FlatDisk )
          accepted->setGeometryDisabled( true );
      };//hooks.beforeAccept
      
      hooks.afterAccept = [modeling,accuracy]( shared_ptr<DetectorPeakResponse> accepted ){
        if( !accepted || !(*modeling) )
          return;
        
        const UrlDrfModeling choice = static_cast<UrlDrfModeling>( (*modeling)->currentIndex() );
        if( !url_modeling_needs_mc(choice) )
          return;
        
        // QuickMc's generation forces the far-field profile, so the profile only distinguishes the
        //  two full characterizations.
        const MakeMcResponseForDrf::Method method = (choice == UrlDrfModeling::PartialMc)
                        ? MakeMcResponseForDrf::Method::QuickMc
                        : MakeMcResponseForDrf::Method::FullMc;
        const ceelo::ResponseProfile profile = (choice == UrlDrfModeling::ExtendedMc)
                        ? ceelo::ResponseProfile::Contact
                        : ceelo::ResponseProfile::General;
        
        MakeMcResponseForDrf::Precision precision = MakeMcResponseForDrf::Precision::Normal;
        if( *accuracy )
        {
          switch( (*accuracy)->currentIndex() )
          {
            case 0:  precision = MakeMcResponseForDrf::Precision::Fast;     break;
            case 2:  precision = MakeMcResponseForDrf::Precision::Thorough; break;
            default: precision = MakeMcResponseForDrf::Precision::Normal;   break;
          }
        }//if( *accuracy )
        
        // Minutes of work, so it goes to the tool that already has progress, ETA and cancel -
        //  opened on the right tab with the run already started, so the user does not have to know
        //  the tool exists, let alone how to drive it.
        InterSpec *viewer = InterSpec::instance();
        if( !viewer )
          return;
        
        DrfModifyWindow *window = viewer->showDrfModifyWindow( accepted );
        if( !window || !window->tool() )
          return;
        
        if( !window->tool()->startMcCharacterization( method, profile, precision ) )
          passMessage( WString::tr("ds-url-mc-not-started"), WarningWidget::WarningMsgHigh );
      };//hooks.afterAccept
    }//if( geom_only )
    
    DrfSelect::createChooseDrfDialog( {drf}, msg, "", nullptr, hooks );
  }catch( std::exception &e )
  {
    wApp->log( "error" ) << "App URL was invalid DRF: " << e.what();
    throw runtime_error( WString::tr("ds-err-drf-uri-invalid").toUTF8() );
  }
}//void DrfSelect::handle_app_url_drf( const std::string &url_query )



namespace
{
/** Undo/redo helpers for the DRF-select dialogs "Modify..." window.

 Both resolve the dialog when the step runs, rather than capturing it: the DRF-select window itself
 is undoable, so by the time one of these executes the `DrfSelect` that recorded it may be a
 different instance (or gone).  See the comments on `UndoRedoManager::addUndoRedoStep`.
 */
void reopenDrfSelectModifyWindow()
{
  InterSpec *viewer = InterSpec::instance();
  DrfSelectWindow *window = viewer ? viewer->showDrfSelectWindow() : nullptr;
  if( window && window->widget() )
    window->widget()->openModifyWindow();
}//reopenDrfSelectModifyWindow()


void closeDrfSelectModifyWindow()
{
  InterSpec *viewer = InterSpec::instance();
  DrfSelectWindow *window = viewer ? viewer->drfSelectWindow() : nullptr;
  if( window && window->widget() )
    window->widget()->programmaticallyCloseModifyWindow();
}//closeDrfSelectModifyWindow()
}//namespace


DrfModifyWindow *DrfSelect::modifyWindow()
{
  return m_modifyWindow.get();
}//DrfModifyWindow *modifyWindow()


void DrfSelect::programmaticallyCloseModifyWindow()
{
  if( !m_modifyWindow )
    return;

  // Null the member first, so the `finished()` handler below early-returns rather than adding a
  //  second (recursive) undo step.
  DrfModifyWindow *window = m_modifyWindow.get();
  m_modifyWindow = nullptr;
  AuxWindow::deleteAuxWindow( window );
}//void programmaticallyCloseModifyWindow()



std::vector<std::string> DrfSelect::potential_rel_eff_det_files()
{
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

  return user_data_paths;
}//DrfSelect::potential_rel_eff_det_files()


std::vector<std::string> DrfSelect::recursive_list_gadras_drfs( const std::string &sourcedir )
{
  vector<string> files;
  if( !SpecUtils::is_directory( sourcedir ) )
    return files;

  const string dat_file = SpecUtils::append_path( sourcedir, "Detector.dat");

  // A Detector.dat is enough on its own: it defines the crystal, and a
  //  Monte-Carlo characterization of that geometry supplies the efficiency an
  //  Efficiency.csv would have.  Such a directory is listed as "geometry only"
  //  and opens the Modify dialog rather than installing a detector that has no
  //  sensitivity yet - see GadrasDirectory::initDetectors.
  if( SpecUtils::is_file(dat_file) )
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
}//DrfSelect::recursive_list_gadras_drfs()


std::vector<std::string> DrfSelect::potential_gadras_det_dirs( InterSpec *interspec )
{
  string pathstr;
  vector<string> paths;

#if( !BUILD_FOR_WEB_DEPLOYMENT && !defined(IOS) )
  try
  {
    if( interspec )
      pathstr = UserPreferences::preferenceValue<string>( "GadrasDRFPath", interspec );
  }catch( std::exception & )
  {
    if( interspec )
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
      auto subsubdirs = DrfSelect::recursive_list_gadras_drfs(subdirpath);
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
        auto subdirs = DrfSelect::recursive_list_gadras_drfs(pstr);
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

  return paths;
}//DrfSelect::potential_gadras_det_dirs()



std::shared_ptr<DetectorPeakResponse> DrfSelect::parseInterSpecRelEffCsvFile( const std::string filename )
{
#ifdef _WIN32
  const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(filename);
  ifstream csvfile( wfilename.c_str(), ios_base::binary | ios_base::in );
#else
  ifstream csvfile( filename.c_str(), ios_base::binary | ios_base::in );
#endif
  
  if( !csvfile.is_open() )
    return nullptr;
  
  return DetectorPeakResponse::parseInterSpecRelEffCsv( csvfile );
}//parseInterSpecRelEffCsvFile(...)



void DrfSelect::updateLastUsedTimeOrAddToDb( std::shared_ptr<DetectorPeakResponse> drf,
                                                long long db_user_id,
                                                std::shared_ptr<DataBaseUtils::DbSession> sql,
                                                const bool appearance_is_authoritative )
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
      // The hash doesnt include the peak-fit prefs, and rows stored before they were saved to the DB
      //  (see `DetectorPeakResponse::drfExtraToXmlString`) dont have them - so give them them now.
      if( !(*iter)->peakFitDetPrefs() && drf->peakFitDetPrefs() )
        (*iter).modify()->setPeakFitDetPrefs( drf->peakFitDetPrefs() );

      // Nor the appearance, which the user may have just attached, changed, or removed.
      const shared_ptr<const DrfAppearance> row_look = (*iter)->appearance(), new_look = drf->appearance();
      const bool same_look = (!row_look && !new_look) || (row_look && new_look && (*row_look == *new_look));
      if( !same_look && (new_look || appearance_is_authoritative) )
        (*iter).modify()->setAppearance( new_look );

      (*iter).modify()->updateLastUsedTimeToNow();
      ++nupdated;
    }//
    
//    cout << "Have seen Detector '" << drf->name() << "' " << nupdated << " times." << endl;
    
    if( !nupdated )
    {
      //Create a separate DetectorPeakResponse because shared_ptr and
      //  dbo::ptr don't work well together
      auto tempDetector = std::make_unique<DetectorPeakResponse>( *drf );
      tempDetector->updateLastUsedTimeToNow();

      tempDetector->m_user = static_cast<int>( db_user_id ); //adds the current user to the detectorpeakresponse insertion into DB
      sql->session()->add( std::move(tempDetector) );
    }
    
    transaction.commit();
  }catch( std::exception &e )
  {
    passMessage( "Error writing DetectorPeakResponse to DB: " + std::string(e.what()), WarningWidget::WarningMsgHigh );
  }//try / catch
}//void updateLastUsedTimeOrAddToDb(...)



std::shared_ptr<DetectorPeakResponse> DrfSelect::getUserPreferredDetector(
                                std::shared_ptr<DataBaseUtils::DbSession> sql,
                                const long long db_user_id,
                                const std::string &serial_number,
                                SpecUtils::DetectorType detType,
                                const std::string &detector_model )
{
  std::shared_ptr<DetectorPeakResponse> answer;
  if( !sql || (db_user_id < 0) )
    return answer;
  
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *sql );
    
    Wt::Dbo::ptr<UseDrfPref> pref;
    
    if( !serial_number.empty() )
    {
      auto results = sql->session()->find<UseDrfPref>()
                     .where( "InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?" )
                     .bind( db_user_id )
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
                       .bind( db_user_id )
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
    .bind( db_user_id )
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
    
    auto newprefOwned = std::make_unique<UseDrfPref>();
    newprefOwned->m_user = user;
    newprefOwned->m_field = prefType;
    newprefOwned->m_criteria = criteria;
    newprefOwned->m_flags = 0x00;
    newprefOwned->m_drfIndex = drflist.front().id();

    Wt::Dbo::ptr<UseDrfPref> newpref = sql->session()->add( std::move(newprefOwned) );
    
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

      case SpecUtils::DetectorType::H3D400:
        return { "H3D-400", "H3D 400", "H3D M400", "H3D A400", "H3D H400", "H3D S400" };

      case SpecUtils::DetectorType::KromekD3S:
        return { "D3S", "Kromek D3S", "D3", "Kromek D3S", "Kromek D3" };

      case SpecUtils::DetectorType::KromekD5:
        return { "D5", "KromekD5", "Kromek-D5", "Kromek D5" };

      case SpecUtils::DetectorType::KromekGR1:
        return { "GR1", "Kromek GR1", "Kromek-GR1" };

      case SpecUtils::DetectorType::RadiaCodeCsI10:
        return { "Radiacode 102", "Radiacode", "Radiacode 101", "Radiacode 103" };
        break;

      case SpecUtils::DetectorType::RadiaCodeCsI14:
        return { "Radiacode 110" };
        break;

      case SpecUtils::DetectorType::RadiaCodeGAGG10:
        return { "Radiacode 103G" };
        break;

      case SpecUtils::DetectorType::Raysid:
        return { "Raysid" };
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


#if( BUILD_AS_OSX_APP )
namespace
{
  void access_bookmarked_locations( const DrfLibrary::LocationSettings &settings );
}
#endif


std::string DrfSelect::gadrasDrfSearchPaths( InterSpec *interspec )
{
#if( BUILD_FOR_WEB_DEPLOYMENT )
  const string datadir = InterSpec::staticDataDirectory();
  return SpecUtils::append_path( datadir, "GenericGadrasDetectors" )
         + ";" + SpecUtils::append_path( datadir, "OUO_GadrasDetectors" );
#else
  string pref;
  try
  {
    pref = UserPreferences::preferenceValue<string>( "DrfLibraryLocations", interspec );
  }catch( std::exception & )
  {
  }

  // The Library's locations have not been set up yet (the Detector Efficiency dialog does that
  //  the first time it opens): the old preference.
  if( pref.empty() )
    return UserPreferences::preferenceValue<string>( "GadrasDRFPath", interspec );

  // The generic detectors always (they are the defaults for each detector type), and the other
  //  GADRAS folders the Library shows.
  const DrfLibrary::LocationSettings settings = DrfLibrary::parseLocationSettings( pref );
#if( BUILD_AS_OSX_APP )
  access_bookmarked_locations( settings );
#endif

  string paths;
  for( const DrfLibrary::Location &loc
                 : DrfLibrary::applyLocationSettings( DrfLibrary::standardLocations(), settings ) )
  {
    const bool generic = (loc.role == DrfLibrary::LocationRole::GenericGadras)
                         || (loc.role == DrfLibrary::LocationRole::ExtraGadras);
    const bool shown_gadras = loc.isFolder && !loc.hidden && (loc.formats & DrfLibrary::FormatGadras);
    if( generic || shown_gadras )
      paths += (paths.empty() ? "" : ";") + loc.path;
  }
  return paths;
#endif
}//std::string gadrasDrfSearchPaths( InterSpec *interspec )


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetector( const SpecUtils::DetectorType type,
                                                                      const std::string &searchPaths )

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
    case DetectorType::H3D400:                     name = "CZT_H3D_M400_ORNL_25cm"; break;
      
    case DetectorType::IdentiFinderR500NaI:        name = "IdentiFINDER-R500-NaI";   break;
    case DetectorType::IdentiFinderR500LaBr:       name = "IdentiFINDER-R500-LaBr3"; break;
      
    case SpecUtils::DetectorType::KromekD3S:       name = "D3S"; break;
    case SpecUtils::DetectorType::KromekD5:        name = "Kromek-D5"; break;
    case SpecUtils::DetectorType::KromekGR1:       name = "Kromek-GR1-CZT"; break;
    case SpecUtils::DetectorType::Raysid:          name = "Raysid"; break;
    case SpecUtils::DetectorType::RadiaCodeCsI10:  name = "Radiacode 102"; break;
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
    case SpecUtils::DetectorType::RadiaCodeCsI14:
    case SpecUtils::DetectorType::RadiaCodeGAGG10:
      // \TODO: fill in these names
      break;
  }//switch( type )
  
  if( name.empty() )
    throw runtime_error( "There is no GADRAS detector response function for a "
                        + SpecUtils::detectorTypeToString(type) );

  std::shared_ptr<DetectorPeakResponse> det;
  
  try
  {
    det = initAGadrasDetector( name, searchPaths );
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
    det = initAGadrasDetector( secondname, searchPaths );
    if( det )
      return det;
  }catch(...)
  {
  }
  
  throw runtime_error( "Could not find GADRAS detector names " + name + " or " + secondname  );
  return det;
}//initAGadrasDetector


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetector( const std::string &currentDetName,
                                                                      const std::string &searchPaths )
{
  //Split searchPaths by semicolon and newlines, then go through and look for sub-folders with the
  //  given name that contain Detector.data and Efficiency.csv.
  vector<string> paths;
  SpecUtils::split( paths, searchPaths, "\r\n;" );
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


std::shared_ptr<DetectorPeakResponse> DrfSelect::initAGadrasDetectorFromDirectory(
                                          const std::string &path,
                                          const bool allow_missing_efficiency_csv )
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
  
  if( thiscsv.empty() && !allow_missing_efficiency_csv )
    throw runtime_error( "Could not find a Efficiency.csv file in '" + path + "'" );
  
  if( thisdat.empty() )
    throw runtime_error( "Could not find a Detector.dat file in '" + path + "'" );
  
  if( thiscsv.empty() )
  {
    // Geometry only: everything the Detector.dat defines, and no efficiency.
    //  The caller is responsible for not presenting this as a finished detector.
    det->fromGadrasDirectory( path, true );
    return det;
  }//if( no Efficiency.csv )

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
  det->setName( SpecUtils::filename( path ) );

//#if( PERFORM_DEVELOPER_CHECKS )
//  check_url_serialization( det );
//#endif

  return det;
}//std::shared_ptr<DetectorPeakResponse> initAGadrasDetectorFromDirectory()



// ---------------------------------------------------------------------------------------------
//  DrfFormulaForm - Import › From a formula
// ---------------------------------------------------------------------------------------------

/** Defines a detector from an efficiency formula, or the coefficients of the standard
 exp( c0 + c1*ln(x) + c2*ln(x)^2 + ... ) form, plus the crystal diameter.  Less common choices (an
 absolute efficiency at a distance, fixed geometries, the energy range) are under "More options".
 */
class DrfFormulaForm : public Wt::WContainerWidget
{
public:
  DrfFormulaForm();

  /** The detector the form currently defines, or null if it does not define a valid one. */
  std::shared_ptr<DetectorPeakResponse> candidate() const { return m_candidate; }

  /** Emitted whenever #candidate may have changed. */
  Wt::Signal<> &changed() { return m_changed; }

protected:
  void setCoefficientMode( const bool coefs );
  void update();

  /** Builds the detector; throws with a message for the user if the inputs do not define one. */
  std::shared_ptr<DetectorPeakResponse> build() const;

  bool m_coefMode;
  WLineEdit *m_name;
  WPushButton *m_formulaBtn, *m_coefBtn;
  WTextArea *m_text;
  WText *m_hint, *m_check;
  WRadioButton *m_kev, *m_mev;
  std::shared_ptr<WButtonGroup> m_units;
  WLineEdit *m_diameter, *m_setback, *m_description;
  WPushButton *m_moreBtn;
  WContainerWidget *m_more;
  WComboBox *m_type;
  WContainerWidget *m_distanceRow;
  WLineEdit *m_distance, *m_minEnergy, *m_maxEnergy;

  std::shared_ptr<DetectorPeakResponse> m_candidate;
  Wt::Signal<> m_changed;
};//class DrfFormulaForm


namespace
{
  /** A labelled row of the formula form. */
  WContainerWidget *form_row( WContainerWidget *parent, const WString &label, WFormWidget *buddy = nullptr )
  {
    WContainerWidget *row = parent->addNew<WContainerWidget>();
    row->addStyleClass( "DrfFormRow" );
    WLabel *lbl = row->addNew<WLabel>( label );
    lbl->addStyleClass( "DrfFormLabel" );
    if( buddy )
      lbl->setBuddy( buddy );
    return row;
  }

  void no_autocorrect( WFormWidget *w )
  {
    w->setAttributeValue( "autocorrect", "off" );
    w->setAttributeValue( "spellcheck", "false" );
  }

  /** Parses an optional energy (keV) field: 0 if empty; throws if not a number. */
  float optional_energy( const WLineEdit *edit )
  {
    const string txt = SpecUtils::trim_copy( edit->text().toUTF8() );
    if( txt.empty() )
      return 0.0f;
    float value = 0.0f;
    if( !SpecUtils::parse_float( txt.c_str(), txt.size(), value ) || (value < 0.0f) )
      throw runtime_error( WString::tr("dsl-f-err-energy").toUTF8() );
    return value;
  }
}//namespace


DrfFormulaForm::DrfFormulaForm()
  : WContainerWidget(),
    m_coefMode( false ),
    m_name( nullptr ), m_formulaBtn( nullptr ), m_coefBtn( nullptr ), m_text( nullptr ),
    m_hint( nullptr ), m_check( nullptr ), m_kev( nullptr ), m_mev( nullptr ), m_units(),
    m_diameter( nullptr ), m_setback( nullptr ), m_description( nullptr ), m_moreBtn( nullptr ),
    m_more( nullptr ), m_type( nullptr ), m_distanceRow( nullptr ), m_distance( nullptr ),
    m_minEnergy( nullptr ), m_maxEnergy( nullptr ), m_candidate(), m_changed()
{
  addStyleClass( "DrfFormulaForm" );

  // Name
  m_name = new WLineEdit();  //owned once added below
  {
    std::unique_ptr<WLineEdit> owned( m_name );
    form_row( this, WString::tr("dsl-f-name"), m_name )->addWidget( std::move(owned) );
  }
  m_name->setText( WString::tr("app-default-manual-det-name") );
  no_autocorrect( m_name );

  // Efficiency: formula or coefficients
  WContainerWidget *eff_row = form_row( this, WString::tr("dsl-f-eff") );
  eff_row->addStyleClass( "Top" );
  WContainerWidget *eff = eff_row->addNew<WContainerWidget>();
  eff->addStyleClass( "DrfFormStack" );
  WContainerWidget *seg = eff->addNew<WContainerWidget>();
  seg->addStyleClass( "DrfSeg Small" );
  seg->setAttributeValue( "role", "group" );
  seg->setAttributeValue( "aria-label", WString::tr("dsl-f-eff-as") );
  m_formulaBtn = seg->addNew<WPushButton>( WString::tr("dsl-f-formula") );
  m_coefBtn = seg->addNew<WPushButton>( WString::tr("dsl-f-coefs") );
  m_formulaBtn->clicked().connect( [this](){ setCoefficientMode( false ); } );
  m_coefBtn->clicked().connect( [this](){ setCoefficientMode( true ); } );

  m_text = eff->addNew<WTextArea>();
  m_text->setRows( 3 );
  m_text->addStyleClass( "DrfMono" );
  m_text->setAttributeValue( "aria-label", WString::tr("dsl-f-eff") );
  no_autocorrect( m_text );
  m_hint = eff->addNew<WText>( "", TextFormat::Plain );
  m_hint->addStyleClass( "DrfFormHint" );
  m_hint->setInline( false );
  m_check = eff->addNew<WText>( "", TextFormat::Plain );
  m_check->addStyleClass( "DrfFormCheck" );
  m_check->setInline( false );

  // Energy units
  WContainerWidget *unit_row = form_row( this, WString::tr("dsl-f-units") );
  WContainerWidget *units = unit_row->addNew<WContainerWidget>();
  units->addStyleClass( "DrfFormInline" );
  m_units = std::make_shared<WButtonGroup>();
  m_kev = units->addNew<WRadioButton>( "keV" );
  m_mev = units->addNew<WRadioButton>( "MeV" );
  m_units->addButton( m_kev, 0 );
  m_units->addButton( m_mev, 1 );
  m_units->setSelectedButtonIndex( 0 );

  // Diameter and setback
  WContainerWidget *diam_row = form_row( this, WString::tr("dsl-f-diam") );
  WContainerWidget *diam = diam_row->addNew<WContainerWidget>();
  diam->addStyleClass( "DrfFormInline" );
  m_diameter = diam->addNew<WLineEdit>( "5.08 cm" );
  m_diameter->addStyleClass( "DrfFormShort" );
  m_diameter->setAttributeValue( "aria-label", WString::tr("dsl-f-diam") );
  WLabel *setback_label = diam->addNew<WLabel>( WString::tr("dsl-f-setback") );
  setback_label->addStyleClass( "DrfFormLabel" );
  m_setback = diam->addNew<WLineEdit>();
  m_setback->addStyleClass( "DrfFormShort" );
  m_setback->setPlaceholderText( "0 cm" );
  setback_label->setBuddy( m_setback );
  diam->addNew<WText>( WString::tr("dsl-f-setback-hint"), TextFormat::Plain )->addStyleClass( "DrfFormHint" );
  no_autocorrect( m_diameter );
  no_autocorrect( m_setback );

  // Description
  m_description = new WLineEdit();
  {
    std::unique_ptr<WLineEdit> owned( m_description );
    form_row( this, WString::tr("dsl-f-desc"), m_description )->addWidget( std::move(owned) );
  }
  m_description->setPlaceholderText( WString::tr("dsl-f-desc-placeholder") );

  // More options: what the efficiency means, and the energy range
  m_moreBtn = addNew<WPushButton>( WString::tr("dsl-f-more") );
  m_moreBtn->addStyleClass( "DrfLinkBtn DrfFormMore" );
  m_more = addNew<WContainerWidget>();
  m_more->addStyleClass( "DrfFormMoreBody" );
  m_more->hide();
  m_moreBtn->clicked().connect( [this](){
    m_more->setHidden( !m_more->isHidden() );
    m_moreBtn->setAttributeValue( "aria-expanded", m_more->isHidden() ? "false" : "true" );
  } );
  m_moreBtn->setAttributeValue( "aria-expanded", "false" );

  m_type = new WComboBox();
  {
    std::unique_ptr<WComboBox> owned( m_type );
    form_row( m_more, WString::tr("dsl-f-type"), m_type )->addWidget( std::move(owned) );
  }
  for( const char *key : { "ds-intrinsic-eff", "ds-abs-eff", "ds-fixed-geom-total",
                           "ds-fixed-geom-cm2", "ds-fixed-geom-m2", "ds-fixed-geom-gram" } )
    m_type->addItem( WString::tr(key) );
  m_type->setCurrentIndex( 0 );

  m_distance = new WLineEdit();
  {
    std::unique_ptr<WLineEdit> owned( m_distance );
    m_distanceRow = form_row( m_more, WString::tr("dsl-f-distance"), m_distance );
    m_distanceRow->addWidget( std::move(owned) );
  }
  m_distance->addStyleClass( "DrfFormShort" );
  m_distance->setPlaceholderText( "25 cm" );
  m_distanceRow->hide();

  WContainerWidget *range_row = form_row( m_more, WString::tr("dsl-f-range") );
  WContainerWidget *range = range_row->addNew<WContainerWidget>();
  range->addStyleClass( "DrfFormInline" );
  m_minEnergy = range->addNew<WLineEdit>();
  m_minEnergy->addStyleClass( "DrfFormShort" );
  m_minEnergy->setPlaceholderText( WString::tr("dsl-f-range-low") );
  m_minEnergy->setAttributeValue( "aria-label", WString::tr("dsl-f-range-low") );
  range->addNew<WText>( WString::fromUTF8("–"), TextFormat::Plain );
  m_maxEnergy = range->addNew<WLineEdit>();
  m_maxEnergy->addStyleClass( "DrfFormShort" );
  m_maxEnergy->setPlaceholderText( WString::tr("dsl-f-range-high") );
  m_maxEnergy->setAttributeValue( "aria-label", WString::tr("dsl-f-range-high") );
  range->addNew<WText>( "keV", TextFormat::Plain );

  // Every edit re-checks the detector
  for( WLineEdit *edit : { m_name, m_diameter, m_setback, m_description, m_distance, m_minEnergy, m_maxEnergy } )
    edit->textInput().connect( this, &DrfFormulaForm::update );
  m_text->textInput().connect( this, &DrfFormulaForm::update );
  m_units->checkedChanged().connect( this, &DrfFormulaForm::update );
  m_type->activated().connect( this, &DrfFormulaForm::update );

  setCoefficientMode( false );
}//DrfFormulaForm constructor


void DrfFormulaForm::setCoefficientMode( const bool coefs )
{
  m_coefMode = coefs;
  m_formulaBtn->toggleStyleClass( "Pressed", !coefs );
  m_formulaBtn->setAttributeValue( "aria-pressed", coefs ? "false" : "true" );
  m_coefBtn->toggleStyleClass( "Pressed", coefs );
  m_coefBtn->setAttributeValue( "aria-pressed", coefs ? "true" : "false" );
  m_hint->setText( WString::tr( coefs ? "dsl-f-coefs-hint" : "dsl-f-formula-hint" ) );
  m_text->setPlaceholderText( WString::tr( coefs ? "dsl-f-coefs-placeholder" : "dsl-f-formula-placeholder" ) );
  update();
}//setCoefficientMode(...)


std::shared_ptr<DetectorPeakResponse> DrfFormulaForm::build() const
{
  string fcn = SpecUtils::trim_copy( m_text->text().toUTF8() );
  if( fcn.empty() )
    throw runtime_error( "" );  //nothing entered yet; not an error to show

  // Coefficients (typed, or pasted from a spreadsheet, possibly with units and dimensions)
  const CoefToFormulaInfo coefs = check_if_coef_for_formula( fcn );
  float energy_units = static_cast<float>( (m_units->checkedId() == 1) ? PhysicalUnits::MeV : PhysicalUnits::keV );
  if( coefs.m_valid_coefficients )
  {
    fcn = coefs.m_fcn;
    if( coefs.m_energy_units > 0.0 )
      energy_units = static_cast<float>( coefs.m_energy_units );
  }else if( m_coefMode )
  {
    throw runtime_error( WString::tr("dsl-f-err-coefs").toUTF8() );
  }

  using GeomType = DetectorPeakResponse::EffGeometryType;
  const int type_index = m_type->currentIndex();
  const bool is_absolute = (type_index == 1);
  const bool is_fixed_geom = (type_index >= 2);
  GeomType geom_type = GeomType::FarFieldIntrinsic;
  switch( type_index )
  {
    case 2: geom_type = GeomType::FixedGeomTotalAct;   break;
    case 3: geom_type = GeomType::FixedGeomActPerCm2;  break;
    case 4: geom_type = GeomType::FixedGeomActPerM2;   break;
    case 5: geom_type = GeomType::FixedGeomActPerGram; break;
    default: break;
  }

  double diameter = 0.0;
  if( !is_fixed_geom )
  {
    try
    {
      diameter = PhysicalUnits::stringToDistance( m_diameter->text().toUTF8() );
    }catch( std::exception & )
    {
      throw runtime_error( WString::tr("ds-err-det-diam-pos").toUTF8() );
    }
    if( diameter < 0.0001 )
      throw runtime_error( WString::tr("ds-err-det-diam-pos").toUTF8() );
  }//if( !is_fixed_geom )

  double setback = 0.0;
  const string setback_txt = SpecUtils::trim_copy( m_setback->text().toUTF8() );
  if( !setback_txt.empty() && !is_fixed_geom )
  {
    try
    {
      setback = std::max( 0.0, PhysicalUnits::stringToDistance( setback_txt ) );
    }catch( std::exception & )
    {
      throw runtime_error( WString::tr("dsl-f-err-setback").toUTF8() );
    }
  }//if( a setback was entered )

  if( is_absolute )
  {
    // An absolute efficiency at a distance is turned into an intrinsic one.
    double distance = 0.0;
    try
    {
      distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
    }catch( std::exception & )
    {
    }
    if( distance < 0.0001 )
      throw runtime_error( WString::tr("ds-err-det-dist-pos").toUTF8() );
    const double gfactor = DetectorPeakResponse::fractionalSolidAngle( diameter, distance + setback );
    fcn = std::to_string( 1.0/gfactor ) + "*(" + fcn + ")";
  }//if( is_absolute )

  float min_energy = optional_energy( m_minEnergy ), max_energy = optional_energy( m_maxEnergy );
  if( (min_energy == 0.0f) && (max_energy == 0.0f) && coefs.m_valid_coefficients
      && (coefs.m_upper_energy > coefs.m_lower_energy) )
  {
    min_energy = coefs.m_lower_energy;
    max_energy = coefs.m_upper_energy;
  }
  if( max_energy <= min_energy )
    min_energy = max_energy = 0.0f;

  string name = SpecUtils::trim_copy( m_name->text().toUTF8() );
  if( name.empty() )
    name = coefs.m_det_name.empty() ? WString::tr("app-default-manual-det-name").toUTF8() : coefs.m_det_name;
  string description = SpecUtils::trim_copy( m_description->text().toUTF8() );
  if( description.empty() )
  {
    description = "f(x)=" + fcn;
    if( description.size() > 25 )
      description = description.substr( 0, 25 ) + "...";
  }

  auto drf = make_shared<DetectorPeakResponse>( name, description );
  drf->setIntrinsicEfficiencyFormula( fcn, static_cast<float>(diameter), energy_units,
                                      min_energy, max_energy, geom_type );
  drf->setDrfSource( DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf );
  if( setback > 0.0 )
    drf->setDetectorSetback( setback );

  // A formula that evaluates to nonsense at a typical energy is not a detector.
  const double eff = drf->farFieldIntrinsicEfficiency( 661.657f );
  if( !std::isfinite( eff ) || (eff <= 0.0) )
    throw runtime_error( WString::tr("dsl-f-err-eval").toUTF8() );

  return drf;
}//build()


void DrfFormulaForm::update()
{
  m_distanceRow->setHidden( m_type->currentIndex() != 1 );

  m_candidate.reset();
  try
  {
    m_candidate = build();
    m_text->removeStyleClass( "Wt-invalid" );

    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%.3g", m_candidate->farFieldIntrinsicEfficiency( 661.657f ) );
    m_check->setText( WString::tr("dsl-f-check-ok").arg( buffer ) );
    m_check->removeStyleClass( "Bad" );
    m_check->show();
  }catch( std::exception &e )
  {
    const string msg = e.what();
    const bool nothing_entered = SpecUtils::trim_copy( m_text->text().toUTF8() ).empty();
    m_text->toggleStyleClass( "Wt-invalid", !nothing_entered );
    m_check->setText( WString::fromUTF8( msg ) );
    m_check->addStyleClass( "Bad" );
    m_check->setHidden( nothing_entered || msg.empty() );
  }//try / catch

  m_changed.emit();
}//update()


// ---------------------------------------------------------------------------------------------
//  DrfSelect
// ---------------------------------------------------------------------------------------------

namespace
{
  const char * const ns_chip_keys[] = { "dsl-chip-all", "dsl-origin-builtin", "dsl-origin-gadras",
    "dsl-origin-isocs", "dsl-origin-angle", "dsl-chip-reltable", "dsl-origin-usermade" };
  const int ns_chip_origins[] = { -1, static_cast<int>(DrfLibrary::Origin::BuiltIn),
    static_cast<int>(DrfLibrary::Origin::Gadras), static_cast<int>(DrfLibrary::Origin::Isocs),
    static_cast<int>(DrfLibrary::Origin::Angle), static_cast<int>(DrfLibrary::Origin::RelEffTable),
    static_cast<int>(DrfLibrary::Origin::UserMade) };


  /** The detector a user chose, for one kind of match, to be picked automatically; null if none.
   Runs off the session thread (the session is mutex guarded).
   */
  shared_ptr<const DetectorPeakResponse> preferred_drf( shared_ptr<DataBaseUtils::DbSession> sql,
                                                        const long long user_id,
                                                        const UseDrfPref::UseDrfType type,
                                                        const string &criteria )
  {
    if( !sql || (user_id < 0) || criteria.empty() )
      return nullptr;

    try
    {
      // Only the column: loading UseDrfPref objects would touch the session's shared InterSpecUser
      //  object (its belongsTo pointer), whose reference count is not thread safe.
      DataBaseUtils::DbTransaction transaction( *sql );
      const Wt::Dbo::collection<long long> indices = sql->session()->query<long long>(
                           "SELECT DrfIndex FROM UseDrfPref"
                           " WHERE InterSpecUser_id = ? AND MatchField = ? AND Criteria = ?"
                           " ORDER BY id DESC" )
                           .bind( user_id ).bind( static_cast<int>(type) ).bind( criteria ).resultList();
      long long index = -1;
      for( const long long i : indices )
      {
        index = i;
        break;
      }

      shared_ptr<const DetectorPeakResponse> answer;
      if( index >= 0 )
      {
        const auto drfs = sql->session()->find<DetectorPeakResponse>()
                            .where( "InterSpecUser_id = ? AND id = ?" )
                            .bind( user_id ).bind( index ).resultList();
        if( !drfs.empty() )
          answer = make_shared<const DetectorPeakResponse>( *drfs.front() );
      }
      transaction.commit();
      return answer;
    }catch( std::exception &e )
    {
      cerr << "DrfSelect: failed reading the preferred detector: " << e.what() << endl;
    }
    return nullptr;
  }//preferred_drf(...)


  /** What the foreground spectrum says about its detector, gathered on the session thread for the
   Library scan's suggestions. */
  struct ForegroundInfo
  {
    string serial, model, manufacturer, modelCriteria;
    SpecUtils::DetectorType type = SpecUtils::DetectorType::Unknown;
    shared_ptr<const DetectorPeakResponse> fileDrf, inUse;
    string gadrasPaths;
    shared_ptr<DataBaseUtils::DbSession> sql;
    long long userId = -1;
  };//struct ForegroundInfo


  /** Detectors to suggest for the foreground spectrum (runs off the session thread). */
  vector<DrfLibrary::Suggestion> make_suggestions( const ForegroundInfo &info )
  {
    vector<DrfLibrary::Suggestion> answer;
    const auto add = [&answer]( shared_ptr<const DetectorPeakResponse> drf, const uint32_t reason ){
      if( !drf )
        return;
      for( DrfLibrary::Suggestion &s : answer )
      {
        if( s.drf->hashValue() == drf->hashValue() )
        {
          s.reasons |= reason;
          return;
        }
      }
      DrfLibrary::Suggestion s;
      s.drf = drf;
      s.reasons = reason;
      answer.push_back( s );
    };//add lambda

    add( preferred_drf( info.sql, info.userId, UseDrfPref::UseDrfType::UseDetectorSerialNumber, info.serial ),
         DrfLibrary::SuggestForSerialNumber );
    add( preferred_drf( info.sql, info.userId, UseDrfPref::UseDrfType::UseDetectorModelName, info.modelCriteria ),
         DrfLibrary::SuggestForModel );
    if( info.fileDrf && info.fileDrf->isValid() )
      add( info.fileDrf, DrfLibrary::SuggestInSpectrumFile );

    if( info.type != SpecUtils::DetectorType::Unknown )
    {
      shared_ptr<const DetectorPeakResponse> def;
      try
      {
        def = DrfSelect::initARelEffDetector( info.type, info.manufacturer, info.model );
      }catch( std::exception & )
      {
      }

      try
      {
        if( !def )
          def = DrfSelect::initAGadrasDetector( info.type, info.gadrasPaths );
      }catch( std::exception & )
      {
      }

      if( def && def->isValid() )
        add( def, DrfLibrary::SuggestDefaultForType );
    }//if( a known detector type )

    // The detector in use is always in the Library, whether or not it is saved anywhere.
    if( info.inUse && info.inUse->isValid() )
      add( info.inUse, DrfLibrary::SuggestNone );

    return answer;
  }//make_suggestions(...)


#if( BUILD_AS_OSX_APP )
  /** The key of the security-scoped bookmark that keeps a user-chosen location readable across
   launches of the (sandboxed) macOS app. */
  string bookmark_key( const string &path )
  {
    string key = path;
    try
    {
      key = SpecUtils::lexically_normalize_path( path );
    }catch( std::exception & )
    {
    }
    for( char &c : key )
    {
      if( (c == '/') || (c == '\\') || (c == ':') )
        c = '_';
    }
    return "SecBookmark_DrfLibrary_" + key;
  }//bookmark_key(...)

  /** Starts accessing the bookmarked user locations - once per run of the app. */
  void access_bookmarked_locations( const DrfLibrary::LocationSettings &settings )
  {
    static std::mutex s_mutex;
    static std::set<string> s_started;

    vector<string> paths;
    for( const DrfLibrary::Location &loc : settings.userLocations )
      paths.push_back( loc.path );
    for( const pair<const string,string> &p : settings.standardPaths )
      paths.push_back( p.second );

    std::lock_guard<std::mutex> lock( s_mutex );
    for( const string &path : paths )
    {
      if( s_started.count( path ) )
        continue;
      string resolved;
      if( macOsUtils::resolveAndStartAccessingSecurityScopedBookmark( bookmark_key( path ), resolved ) )
        s_started.insert( path );
    }
  }//access_bookmarked_locations(...)
#endif //BUILD_AS_OSX_APP


  /** The user's changes to the Library's locations.  The first time, the folders of the old
   "GadrasDRFPath" preference (beyond its defaults) become locations the user added.
   */
  DrfLibrary::LocationSettings load_location_settings( InterSpec *viewer )
  {
    DrfLibrary::LocationSettings settings;
#if( !BUILD_FOR_WEB_DEPLOYMENT )
    try
    {
      // Read before any write (see UserPreferences::setPreferenceValue on new preferences).
      const string pref = UserPreferences::preferenceValue<string>( "DrfLibraryLocations", viewer );
      settings = DrfLibrary::parseLocationSettings( pref );

      if( pref.empty() )
      {
        string old_pref;
        try
        {
          old_pref = UserPreferences::preferenceValue<string>( "GadrasDRFPath", viewer );
        }catch( std::exception & )
        {
        }

        set<string> standard_paths;
        for( const DrfLibrary::Location &loc : DrfLibrary::standardLocations() )
          standard_paths.insert( SpecUtils::lexically_normalize_path( loc.path ) );

        vector<string> paths;
        SpecUtils::split( paths, old_pref, "\r\n;" );
        for( string path : paths )
        {
          // A folder missing now (an unmounted drive, say) is kept - the Locations page shows it
          //  as not found - unless it is a relative path that could not be resolved.
          path = complete_drf_dir_path( SpecUtils::trim_copy( path ) );
          if( path.empty() || (!SpecUtils::is_directory( path ) && !SpecUtils::is_absolute_path( path ))
              || standard_paths.count( SpecUtils::lexically_normalize_path( path ) ) )
            continue;

          DrfLibrary::Location loc;
          loc.id = loc.path = path;
          loc.formats = DrfLibrary::FormatGadras;
          settings.userLocations.push_back( loc );
        }//for( string path : paths )

        UserPreferences::setPreferenceValue<string>( "DrfLibraryLocations",
                                     DrfLibrary::serializeLocationSettings( settings ), viewer );
      }//if( never saved )
    }catch( std::exception &e )
    {
      cerr << "DrfSelect: could not read the Library locations: " << e.what() << endl;
    }
#endif
    return settings;
  }//load_location_settings(...)


  /** Saves the Library location settings; false (with a message to the user) if they are too long. */
  bool save_location_settings( InterSpec *viewer, const DrfLibrary::LocationSettings &settings )
  {
#if( BUILD_FOR_WEB_DEPLOYMENT )
    return false;
#else
    const string value = DrfLibrary::serializeLocationSettings( settings );
    if( value.size() > UserOption::sm_max_value_str_len )
    {
      passMessage( WString::tr("dsl-loc-too-many"), WarningWidget::WarningMsgHigh );
      return false;
    }
    UserPreferences::setPreferenceValue<string>( "DrfLibraryLocations", value, viewer );
    return true;
#endif
  }//save_location_settings(...)


  /** The Library's locations: the standard ones, with the user's changes, then the user's own. */
  vector<DrfLibrary::Location> library_locations( InterSpec *viewer )
  {
    const DrfLibrary::LocationSettings settings = load_location_settings( viewer );
#if( BUILD_AS_OSX_APP )
    access_bookmarked_locations( settings );
#endif
    return DrfLibrary::applyLocationSettings( DrfLibrary::standardLocations(), settings );
  }//library_locations(...)


  /** A name shortened for a button. */
  string short_name( const string &name )
  {
    if( SpecUtils::utf8_str_len( name.c_str(), name.size() ) <= 30 )
      return name;
    string answer = name;
    SpecUtils::utf8_limit_str_size( answer, 29 );
    return answer + "…";
  }//short_name(...)


  /** A serial number without a leading "Serial" or "S/N" some files include (e.g. "Serial #1XXX"),
   since the text it goes in says what it is. */
  string without_serial_prefix( const string &serial )
  {
    for( const char *prefix : { "serial", "s/n", "sn:" } )
    {
      if( SpecUtils::istarts_with( serial, prefix ) )
      {
        const string rest = SpecUtils::trim_copy( serial.substr( strlen(prefix) ) );
        if( !rest.empty() )
          return (rest[0] == ':') ? SpecUtils::trim_copy( rest.substr(1) ) : rest;
      }
    }
    return serial;
  }//without_serial_prefix(...)


  /** Whether an entry IS the detector with this hash: its own hash, or the hash it has under its
   sibling source (see Summary::sourceAliasHashes) - not the hash of a database row folded into it,
   which may be a different detector of the same name. */
  bool has_hash( const DrfLibrary::Entry &entry, const uint64_t hash )
  {
    if( !entry.summary )
      return false;
    const vector<uint64_t> &aliases = entry.summary->sourceAliasHashes;
    return (entry.summary->hash == hash) || (std::find( begin(aliases), end(aliases), hash ) != end(aliases));
  }
}//namespace


DrfSelect::DrfSelect( std::shared_ptr<DetectorPeakResponse> currentDet,
                      InterSpec *specViewer,
                      AuxWindow *parentWindow )
  : WContainerWidget(),
    m_interspec( specViewer ),
    m_sql(),
    m_inUse( currentDet ),
    m_justChanged( false ),
    m_justUsedHash( 0 ),
    m_tab( Tab::Library ),
    m_narrow( false ),
    m_showDetails( false ),
    m_detailsHeader( nullptr ),
    m_detailsBack( nullptr ),
    m_actionBar( nullptr ),
    m_entries(),
    m_pickedId(),
    m_userPicked( false ),
    m_pickInUseAfterScan( false ),
    m_scrolledToPicked( false ),
    m_pickedDrf(),
    m_pickedError(),
    m_suggestedNote(),
    m_scanGeneration( 0 ),
    m_loadGeneration( 0 ),
    m_scanCancel(),
    m_scanning( false ),
    m_originChip( -1 ),
    m_withGeometry( false ),
    m_stripInUse( nullptr ), m_stripNone( nullptr ),
    m_stripIcon( nullptr ), m_stripJust( nullptr ), m_stripName( nullptr ), m_stripMeta( nullptr ),
    m_libTab( nullptr ), m_importTab( nullptr ), m_stack( nullptr ),
    m_list( nullptr ), m_search( nullptr ), m_searchTimer( nullptr ), m_checkingTimer( nullptr ),
    m_chips(), m_geoChip( nullptr ),
    m_fileModeBtn( nullptr ), m_formulaModeBtn( nullptr ), m_importStack( nullptr ),
    m_importWidget( nullptr ), m_formats( nullptr ), m_foundList( nullptr ), m_formula( nullptr ),
    m_batch(), m_found(), m_foundPick( std::string::npos ), m_foundFolder(), m_batchGeneration( 0 ),
    m_useAfterScanRoot(), m_useAfterScanName(), m_useAfterScanFallback(), m_useWhenLoadedId(),
    m_closeAfterUse( false ),
    m_preview( nullptr ),
    m_manageLocations( nullptr ),
    m_locationsPage( nullptr ),
    m_scanLocations(),
    m_locationStatuses(),
    m_autoPick( nullptr ), m_autoSerial( nullptr ), m_autoModel( nullptr ),
    m_cancelBtn( nullptr ), m_primaryBtn( nullptr ),
    m_notice( nullptr ), m_noticeText( nullptr ),
    m_modifyWindow( nullptr ),
    m_doneSignal()
{
  assert( m_interspec );
  if( !m_interspec )
    throw runtime_error( "DrfSelect: no InterSpec" );

  m_interspec->useMessageResourceBundle( "DrfSelect" );
  wApp->useStyleSheet( "InterSpec_resources/DrfSelect.css" );
  addStyleClass( "DrfSelect" );

  m_sql = m_interspec->sql();

  m_interspec->detectorChanged().connect( this, &DrfSelect::setDetector );
  m_interspec->detectorModified().connect( this, &DrfSelect::setDetector );

  // A phone, or a window too narrow for the list beside the preview: one screen at a time.
  {
    const int width = m_interspec->renderedWidth();
    m_narrow = m_interspec->isPhone() || ((width > 100) && (width < 860));
    toggleStyleClass( "DrfNarrow", m_narrow );
  }

  // ---- Details screen header (phone): back to the Library or Import
  m_detailsHeader = addNew<WContainerWidget>();
  m_detailsHeader->addStyleClass( "DrfDetailsHead" );
  m_detailsBack = m_detailsHeader->addNew<WPushButton>( WString::tr("dsl-back-library") );
  m_detailsBack->addStyleClass( "DrfLinkBtn DrfDetailsBack" );
  m_detailsBack->clicked().connect( [this](){ setDetailsShown( false ); } );
  m_detailsHeader->addNew<WText>( WString::tr("dsl-details"), TextFormat::Plain )->addStyleClass( "DrfDetailsTitle" );
  m_detailsHeader->hide();

  // ---- "In use now" strip
  WContainerWidget *strip = addNew<WContainerWidget>();
  strip->addStyleClass( "DrfStrip" );

  m_stripInUse = strip->addNew<WContainerWidget>();
  m_stripInUse->addStyleClass( "DrfStripRow" );
  m_stripIcon = m_stripInUse->addNew<WText>( "", TextFormat::UnsafeXHTML );  //see DrfAppearance::icon_html
  m_stripIcon->addStyleClass( "DrfStripIcon" );
  WContainerWidget *strip_text = m_stripInUse->addNew<WContainerWidget>();
  strip_text->addStyleClass( "DrfStripText" );
  WContainerWidget *strip_label = strip_text->addNew<WContainerWidget>();
  strip_label->addStyleClass( "DrfStripLabelLine" );
  strip_label->addNew<WText>( WString::tr("dsl-in-use-now"), TextFormat::Plain )->addStyleClass( "DrfStripLabel" );
  m_stripJust = strip_label->addNew<WText>( WString::tr("dsl-just-changed"), TextFormat::Plain );
  m_stripJust->addStyleClass( "DrfJustChanged" );
  m_stripName = strip_text->addNew<WText>( "", TextFormat::Plain );
  m_stripName->addStyleClass( "DrfStripName" );
  m_stripName->setInline( false );
  m_stripMeta = strip_text->addNew<WText>( "", TextFormat::Plain );
  m_stripMeta->addStyleClass( "DrfStripMeta" );
  m_stripMeta->setInline( false );
  WContainerWidget *strip_btns = m_stripInUse->addNew<WContainerWidget>();
  strip_btns->addStyleClass( "DrfStripBtns" );
  WPushButton *show_btn = strip_btns->addNew<WPushButton>( WString::tr("dsl-show-in-list") );
  show_btn->addStyleClass( "DrfBtn" );
  show_btn->clicked().connect( this, &DrfSelect::showInUseInList );
  WPushButton *none_btn = strip_btns->addNew<WPushButton>( WString::tr("dsl-use-none") );
  none_btn->addStyleClass( "DrfBtn" );
  none_btn->clicked().connect( this, &DrfSelect::useNoDetector );

  m_stripNone = strip->addNew<WContainerWidget>();
  m_stripNone->addStyleClass( "DrfStripRow None" );
  m_stripNone->addNew<WText>( "<span class=\"DrfIcon DrfIconGeneric\"></span>", TextFormat::UnsafeXHTML )
                             ->addStyleClass( "DrfStripIcon" );
  WContainerWidget *none_text = m_stripNone->addNew<WContainerWidget>();
  none_text->addStyleClass( "DrfStripText" );
  none_text->addNew<WText>( WString::tr("dsl-no-det-in-use"), TextFormat::Plain )->addStyleClass( "DrfStripLabel" );
  WText *none_desc = none_text->addNew<WText>( WString::tr("dsl-no-det-in-use-text"), TextFormat::Plain );
  none_desc->addStyleClass( "DrfStripNoneText" );
  none_desc->setInline( false );

  // ---- Tabs
  WContainerWidget *tabs = addNew<WContainerWidget>();
  tabs->addStyleClass( "DrfTabs" );
  tabs->setAttributeValue( "role", "tablist" );
  m_libTab = tabs->addNew<WPushButton>( WString::tr("dsl-tab-library"), TextFormat::XHTML );
  m_libTab->addStyleClass( "DrfTab" );
  m_libTab->setAttributeValue( "role", "tab" );
  m_libTab->clicked().connect( [this](){ showTab( Tab::Library ); } );
  m_importTab = tabs->addNew<WPushButton>( WString::tr("dsl-tab-import") );
  m_importTab->addStyleClass( "DrfTab DrfTabImport" );
  m_importTab->setAttributeValue( "role", "tab" );
  m_importTab->clicked().connect( [this](){ showTab( Tab::Import ); } );

  // ---- Body: the Library or Import, beside the preview
  WContainerWidget *body = addNew<WContainerWidget>();
  body->addStyleClass( "DrfBody" );
  m_stack = body->addNew<WStackedWidget>();
  m_stack->addStyleClass( "DrfMain" );

  // Library
  {
    auto lib_owned = std::make_unique<WContainerWidget>();
    WContainerWidget *lib = lib_owned.get();
    lib->addStyleClass( "DrfLibPane" );
    m_stack->addWidget( std::move(lib_owned) );

    m_list = lib->addNew<DrfLibraryList>();
    m_list->setNarrow( m_narrow );
    m_list->setExtraSearchText( []( const DrfLibrary::Entry &entry ) -> string {
      return DrfLibraryText::originLabel( entry.origin ).toUTF8();
    } );

    WContainerWidget *search_row = lib->addNew<WContainerWidget>();
    search_row->addStyleClass( "DrfSearchRow" );
    m_search = search_row->addNew<WLineEdit>();
    m_search->addStyleClass( "DrfSearch" );
    m_search->setPlaceholderText( WString::tr("dsl-search-placeholder") );
    m_search->setAttributeValue( "aria-label", WString::tr("dsl-search-label") );
    m_search->setAttributeValue( "type", "search" );
    no_autocorrect( m_search );

    // Adding folders the app reads needs a file system the user can see: not on the web or phones.
    m_manageLocations = search_row->addNew<WPushButton>( WString::tr("dsl-manage-locations") );
    m_manageLocations->addStyleClass( "DrfBtn DrfManageLocations" );
    m_manageLocations->clicked().connect( [this](){ showTab( Tab::Locations ); } );
#if( BUILD_FOR_WEB_DEPLOYMENT || IOS || ANDROID )
    m_manageLocations->hide();
#endif

    WContainerWidget *chips = lib->addNew<WContainerWidget>();
    chips->addStyleClass( "DrfChips" );
    chips->setAttributeValue( "role", "group" );
    chips->setAttributeValue( "aria-label", WString::tr("dsl-chips-label") );
    for( size_t i = 0; i < sizeof(ns_chip_origins)/sizeof(ns_chip_origins[0]); ++i )
    {
      WPushButton *chip = chips->addNew<WPushButton>( WString::tr(ns_chip_keys[i]), TextFormat::XHTML );
      chip->addStyleClass( "DrfChip" );
      const int origin = ns_chip_origins[i];
      chip->clicked().connect( [this, origin](){
        m_originChip = origin;
        updateChips();
        updateFilter();
      } );
      m_chips.push_back( { origin, chip } );
    }
    chips->addNew<WContainerWidget>()->addStyleClass( "DrfChipDivider" );
    m_geoChip = chips->addNew<WPushButton>( WString::tr("dsl-chip-geometry"), TextFormat::XHTML );
    m_geoChip->addStyleClass( "DrfChip DrfChipGeometry" );
    m_geoChip->clicked().connect( [this](){
      m_withGeometry = !m_withGeometry;
      updateChips();
      updateFilter();
    } );
  }

  // Import
  {
    auto imp_owned = std::make_unique<WContainerWidget>();
    WContainerWidget *imp = imp_owned.get();
    imp->addStyleClass( "DrfImportPane" );
    m_stack->addWidget( std::move(imp_owned) );

    WContainerWidget *top = imp->addNew<WContainerWidget>();
    top->addStyleClass( "DrfImportTop" );
    WContainerWidget *seg = top->addNew<WContainerWidget>();
    seg->addStyleClass( "DrfSeg Large" );
    seg->setAttributeValue( "role", "group" );
    seg->setAttributeValue( "aria-label", WString::tr("dsl-import-how") );
    m_fileModeBtn = seg->addNew<WPushButton>( WString::tr("dsl-from-file") );
    m_fileModeBtn->addStyleClass( "DrfSegFile" );
    m_formulaModeBtn = seg->addNew<WPushButton>( WString::tr("dsl-from-formula") );
    m_formulaModeBtn->addStyleClass( "DrfSegFormula" );

    WPushButton *make_drf = top->addNew<WPushButton>( WString::tr("dsl-make-from-cal") );
    make_drf->addStyleClass( "DrfLinkBtn DrfMakeDrf" );
    make_drf->clicked().connect( [this](){
      InterSpec * const viewer = m_interspec;
      viewer->showMakeDrfWindow();
      m_doneSignal.emit();  //closes (deletes) this dialog
    } );

    m_importStack = imp->addNew<WStackedWidget>();
    m_importStack->addStyleClass( "DrfImportStack" );
    // From a file: the drop area (DrfImportWidget, for one detector), the formats it reads, and -
    //  after a drop of several detectors' files - the detectors found.
    WContainerWidget *file_page = m_importStack->addNew<WContainerWidget>();
    file_page->addStyleClass( "DrfFilePage" );
    // Not given the current detector: this tab is for bringing in a new one.
    m_importWidget = file_page->addNew<DrfImportWidget>( DrfImportWidget::Host::ImportTab );
    m_foundList = file_page->addNew<DrfFoundList>();
    m_foundList->hide();

    // On a phone, files can also arrive from other apps
    if( m_interspec->isMobile() )
    {
      WText *share_note = file_page->addNew<WText>( WString::tr("dsl-share-note"), TextFormat::Plain );
      share_note->addStyleClass( "DrfFormHint DrfShareNote" );
      share_note->setInline( false );
    }

    m_formats = file_page->addNew<WContainerWidget>();
    m_formats->addStyleClass( "DrfFormats" );
    WText *formats_title = m_formats->addNew<WText>( WString::tr("dsl-fmt-title"), TextFormat::Plain );
    formats_title->addStyleClass( "DrfSectionTitle" );
    formats_title->setInline( false );
    WContainerWidget *tiles = m_formats->addNew<WContainerWidget>();
    tiles->addStyleClass( "DrfFormatTiles" );
    for( const char *key : { "gadras", "isocs", "angle", "table", "interspec", "csv" } )
    {
      WContainerWidget *tile = tiles->addNew<WContainerWidget>();
      tile->addStyleClass( "DrfFormatTile" );
      WText *name = tile->addNew<WText>( WString::tr( string("dsl-fmt-") + key ), TextFormat::Plain );
      name->addStyleClass( "DrfFormatName" );
      name->setInline( false );
      WText *desc = tile->addNew<WText>( WString::tr( string("dsl-fmt-") + key + "-desc" ), TextFormat::Plain );
      desc->addStyleClass( "DrfFormHint" );
      desc->setInline( false );
    }

    m_formula = m_importStack->addNew<DrfFormulaForm>();

    const auto set_mode = [this]( const bool formula ){
      m_importStack->setCurrentIndex( formula ? 1 : 0 );
      m_fileModeBtn->toggleStyleClass( "Pressed", !formula );
      m_fileModeBtn->setAttributeValue( "aria-pressed", formula ? "false" : "true" );
      m_formulaModeBtn->toggleStyleClass( "Pressed", formula );
      m_formulaModeBtn->setAttributeValue( "aria-pressed", formula ? "true" : "false" );
      handleImportChanged();
    };
    m_fileModeBtn->clicked().connect( [set_mode](){ set_mode( false ); } );
    m_formulaModeBtn->clicked().connect( [set_mode](){ set_mode( true ); } );
    set_mode( false );

    m_importWidget->changed().connect( this, &DrfSelect::handleImportChanged );
    m_importWidget->characterizeRequested().connect( this, &DrfSelect::openModifyWindow );
    m_importWidget->batchDropped().connect( this, &DrfSelect::handleBatchDropped );
    m_foundList->picked().connect( [this]( const size_t index ){
      pickFound( index );
      if( m_narrow )
        setDetailsShown( true );
    } );
    m_foundList->cleared().connect( this, &DrfSelect::clearFound );
    m_formula->changed().connect( this, &DrfSelect::handleImportChanged );
  }

  // Library › Locations
  {
    auto page_owned = std::make_unique<DrfLocationsPage>();
    m_locationsPage = page_owned.get();
    m_stack->addWidget( std::move(page_owned) );
    m_locationsPage->backRequested().connect( [this](){ showTab( Tab::Library ); } );
    m_locationsPage->shownChanged().connect( this, &DrfSelect::setLocationShown );
    m_locationsPage->rescanRequested().connect( this, &DrfSelect::rescanLocation );
    m_locationsPage->removeRequested().connect( this, &DrfSelect::confirmRemoveLocation );
    m_locationsPage->locateRequested().connect( this, &DrfSelect::locateLocation );
    m_locationsPage->addFolderRequested().connect( this, &DrfSelect::startAddLocationFolder );
    m_locationsPage->addFileRequested().connect( this, &DrfSelect::startAddLocationFile );
  }

  m_preview = body->addNew<DrfPreviewPanel>( DrfPreviewPanel::Context::Library );
  m_preview->setReference( m_inUse );
  m_preview->setNarrow( m_narrow );
  if( m_narrow )
  {
    m_manageLocations->hide();
    // The "In use now" card opens its details
    m_stripInUse->clicked().connect( [this](){
      const string id = entryIdFor( m_inUse );
      if( !id.empty() )
        pickEntry( id );
    } );
    // A card that acts as a button, for keyboards and screen readers too
    m_stripInUse->setAttributeValue( "role", "button" );
    m_stripInUse->setAttributeValue( "tabindex", "0" );
    m_stripInUse->keyWentDown().connect( "function(s,e){ if( (e.key === 'Enter') || (e.key === ' ') ){ e.preventDefault(); s.click(); } }" );
    m_preview->useNoneRequested().connect( [this](){
      useNoDetector();
      setDetailsShown( false );  //back to the Library, where the notice and "no detector" card show
    } );
  }

  // ---- A notice (e.g., after "Use none")
  m_notice = addNew<WContainerWidget>();
  m_notice->addStyleClass( "DrfNotice" );
  m_notice->setAttributeValue( "role", "status" );
  m_noticeText = m_notice->addNew<WText>( "", TextFormat::Plain );
  WPushButton *dismiss = m_notice->addNew<WPushButton>();
  dismiss->addStyleClass( "DrfNoticeClose" );
  dismiss->setAttributeValue( "aria-label", WString::tr("dsl-dismiss") );
  dismiss->clicked().connect( m_notice, &WWidget::hide );
  m_notice->hide();

  // ---- Footer: the auto-pick checkboxes, Cancel, and the main button
  WContainerWidget *footer = parentWindow ? parentWindow->footer() : addNew<WContainerWidget>();
  if( !parentWindow )
    footer->addStyleClass( "DrfFooter" );
  AuxWindow::addHelpInFooter( footer, "detector-edit-dialog" );

  // On a phone the window's footer is at the top of the screen; what acts on the detector shown
  //  goes in a bar at the bottom instead.
  WContainerWidget * const action_bar = m_narrow ? addNew<WContainerWidget>() : footer;
  if( m_narrow )
  {
    action_bar->addStyleClass( "DrfBottomBar" );
    m_actionBar = action_bar;
  }

  m_autoPick = action_bar->addNew<WContainerWidget>();
  m_autoPick->addStyleClass( "DrfAutoPick" );
  m_autoPick->addNew<WText>( WString::tr("dsl-auto-pick"), TextFormat::Plain )->addStyleClass( "DrfAutoPickLabel" );

  const shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  string serial, model;
  if( meas )
  {
    serial = meas->instrument_id();
    model = (meas->detector_type() != SpecUtils::DetectorType::Unknown)
            ? SpecUtils::detectorTypeToString( meas->detector_type() ) : meas->instrument_model();
  }
  if( !serial.empty() )
    m_autoSerial = m_autoPick->addNew<WCheckBox>( WString::tr("dsl-auto-serial").arg( without_serial_prefix(serial) ) );
  if( !model.empty() )
    m_autoModel = m_autoPick->addNew<WCheckBox>( WString::tr("dsl-auto-model").arg( model ) );

  if( !model.empty() && !serial.empty() )
    m_suggestedNote = model + " · " + WString::tr("dsl-note-serial").arg( without_serial_prefix(serial) ).toUTF8();
  else if( !model.empty() || !serial.empty() )
    m_suggestedNote = model.empty() ? WString::tr("dsl-note-serial").arg( without_serial_prefix(serial) ).toUTF8() : model;
  m_list->setSuggestedNote( WString::fromUTF8( m_suggestedNote ) );

  m_cancelBtn = footer->addNew<WPushButton>( WString::tr( m_narrow ? "dsl-close" : "Cancel" ) );
  WidgetUtils::applyButtonRole( m_cancelBtn, WidgetUtils::ButtonRole::Dismiss );
  m_cancelBtn->clicked().connect( [this](){ m_doneSignal.emit(); } );

  m_primaryBtn = action_bar->addNew<WPushButton>();
  WidgetUtils::applyButtonRole( m_primaryBtn, WidgetUtils::ButtonRole::Affirm );
  m_primaryBtn->addStyleClass( "DrfPrimary" );
  m_primaryBtn->clicked().connect( this, &DrfSelect::handlePrimary );

  // ---- Timers: search after a pause in typing; "Checking for changes" only if a scan is slow
  m_searchTimer = addChild( std::make_unique<WTimer>() );
  m_searchTimer->setSingleShot( true );
  m_searchTimer->setInterval( std::chrono::milliseconds(250) );
  m_searchTimer->timeout().connect( this, &DrfSelect::updateFilter );
  m_search->textInput().connect( [this](){
    m_searchTimer->stop();
    m_searchTimer->start();
  } );
  m_search->enterPressed().connect( this, &DrfSelect::updateFilter );

  m_checkingTimer = addChild( std::make_unique<WTimer>() );
  m_checkingTimer->setSingleShot( true );
  m_checkingTimer->setInterval( std::chrono::milliseconds(150) );
  m_checkingTimer->timeout().connect( [this](){
    if( m_scanning && !m_entries )
    {
      m_list->setChecking( true );
      m_list->rebuild();
    }
  } );

  // ---- Wiring the list and preview
  m_list->picked().connect( this, &DrfSelect::pickEntry );
  m_list->importRequested().connect( [this](){ showTab( Tab::Import ); } );
  m_preview->customizeRequested().connect( this, &DrfSelect::openModifyWindow );
  m_preview->qrRequested().connect( this, &DrfSelect::showQrCode );
  m_preview->deleteRequested().connect( this, &DrfSelect::confirmDeletePicked );
  m_preview->parentRequested().connect( [this](){
    const DrfLibrary::Entry * const entry = pickedEntry();
    if( !entry || !m_entries )
      return;
    const int parent = DrfLibrary::parentIndex( *m_entries, static_cast<size_t>(entry - &(*m_entries)[0]) );
    if( parent >= 0 )
    {
      pickEntry( (*m_entries)[parent].id );
      m_list->scrollToPicked();
    }
  } );

  // Until the Library is read, preview the detector in use - it is already in memory.
  updateInUseStrip();
  showTab( Tab::Library );
  updateChips();
  startScan();
}//DrfSelect constructor


DrfSelect::~DrfSelect()
{
  if( m_scanCancel )
    m_scanCancel->store( true );

  if( m_modifyWindow )
    AuxWindow::deleteAuxWindow( m_modifyWindow.get() );
}//~DrfSelect()


Wt::Signal<> &DrfSelect::done()
{
  return m_doneSignal;
}


std::shared_ptr<DetectorPeakResponse> DrfSelect::detector()
{
  return m_inUse;
}


std::shared_ptr<const DetectorPeakResponse> DrfSelect::previewedDetector() const
{
  return m_preview ? m_preview->detector() : nullptr;
}


void DrfSelect::startScan()
{
  if( m_scanCancel )
    m_scanCancel->store( true );
  m_scanCancel = std::make_shared<std::atomic<bool>>( false );
  const uint64_t generation = ++m_scanGeneration;
  m_scanning = true;
  m_checkingTimer->start();

  // Everything the worker needs, gathered here on the session thread.
  DrfLibrary::ScanRequest request;
  request.locations = library_locations( m_interspec );
  m_scanLocations = make_shared<const vector<DrfLibrary::Location>>( request.locations );
  m_locationStatuses.reset();
  updateLocationsPage();
  request.sql = m_sql;
  request.dbUserId = m_interspec->user() ? m_interspec->user().id() : -1;
  request.cancel = m_scanCancel;

  ForegroundInfo info;
  const shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( meas )
  {
    info.serial = meas->instrument_id();
    info.model = meas->instrument_model();
    info.manufacturer = meas->manufacturer();
    info.type = meas->detector_type();
    info.modelCriteria = (info.type == SpecUtils::DetectorType::Unknown)
                         ? info.model : SpecUtils::detectorTypeToString( info.type );
    info.fileDrf = meas->detectorFromFile();
  }
  info.inUse = m_inUse ? make_shared<const DetectorPeakResponse>( *m_inUse ) : nullptr;
  info.gadrasPaths = gadrasDrfSearchPaths( m_interspec );
  info.sql = m_sql;
  info.userId = request.dbUserId;

  const WidgetUtils::WidgetHandle self( this );
  const string session_id = wApp->sessionId();
  request.progress = [self, session_id, generation]( const DrfLibrary::Snapshot &partial ){
    const auto copy = make_shared<const DrfLibrary::Snapshot>( partial );
    WServer::instance()->post( session_id, [self, generation, copy](){
      DrfSelect * const widget = self.resolve_as<DrfSelect>();
      if( !widget || (widget->m_scanGeneration != generation) )
        return;
      widget->receiveSnapshot( copy, false );
      if( wApp )
        wApp->triggerUpdate();
    } );
  };//request.progress

  const auto result = make_shared<DrfLibrary::Snapshot>();
  AppUtils::run_on_ioservice( [request, info, result]() mutable {
      request.suggestions = make_suggestions( info );
      *result = DrfLibrary::scan( request );
    },
    [self, generation, result](){
      DrfSelect * const widget = self.resolve_as<DrfSelect>();
      if( !widget || (widget->m_scanGeneration != generation) )
        return;
      widget->receiveSnapshot( result, true );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//startScan()


void DrfSelect::receiveSnapshot( std::shared_ptr<const DrfLibrary::Snapshot> snapshot, const bool final )
{
  if( !snapshot )
    return;

  if( final )
  {
    m_scanning = false;
    m_checkingTimer->stop();
    m_locationStatuses = make_shared<const vector<DrfLibrary::LocationStatus>>( snapshot->locations );
  }

  m_entries = make_shared<const vector<DrfLibrary::Entry>>( snapshot->entries );

  // Until the user picks one, the preview follows the detector in use - whose entry can change
  //  between scans (e.g., from an in-memory one, to the file it turns out to be from).
  const string in_use_id = entryIdFor( m_inUse );
  if( m_pickInUseAfterScan && final )
  {
    m_pickInUseAfterScan = false;
    m_userPicked = false;
  }
  string pick = m_userPicked ? m_pickedId : in_use_id;

  // A picked detector can be listed under a new id after a scan (e.g. the spectrum's own detector,
  //  once used, is a database row): find it again by what it is.
  if( m_userPicked && !pick.empty() && m_pickedDrf
      && std::none_of( begin(*m_entries), end(*m_entries),
                       [&pick]( const DrfLibrary::Entry &e ){ return e.id == pick; } ) )
  {
    const string again = entryIdFor( m_pickedDrf );
    if( !again.empty() )
    {
      pick = again;
      m_pickedId = again;  //same detector: nothing to reload
    }
  }

  m_list->setChecking( false );
  m_list->setEntries( m_entries );
  m_list->setInUse( in_use_id );
  m_list->setPicked( pick );
  m_list->rebuild();

  // The picked detector's entry may be new (or newly read): load it if need be.
  const bool changed = (pick != m_pickedId);
  if( changed )
    m_pickedDrf.reset();
  const DrfLibrary::Entry *entry = nullptr;
  for( size_t i = 0; !pick.empty() && (i < m_entries->size()); ++i )
  {
    if( (*m_entries)[i].id == pick )
      entry = &(*m_entries)[i];
  }
  // Also reload when the entry's file changed under it (e.g. edited, then "Rescan").
  const bool stale = entry && entry->summary && m_pickedDrf && !has_hash( *entry, m_pickedDrf->hashValue() );
  if( stale )
    m_pickedDrf.reset();
  if( changed || (entry && entry->summary && !m_pickedDrf) )
    selectEntry( pick );

  if( !m_scrolledToPicked && !m_pickedId.empty() && final )
  {
    m_scrolledToPicked = true;
    m_list->scrollToPicked();
  }

  updateInUseStrip();  //its source may only now be known (e.g., a just-saved customized copy)
  updateLocationsPage();
  updateChips();
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();

  // Last: on a phone, using the kept detector closes (deletes) this dialog.
  if( final )
    useKeptEntry();
}//receiveSnapshot(...)


const DrfLibrary::Entry *DrfSelect::pickedEntry() const
{
  if( !m_entries || m_pickedId.empty() )
    return nullptr;
  for( const DrfLibrary::Entry &entry : *m_entries )
  {
    if( entry.id == m_pickedId )
      return &entry;
  }
  return nullptr;
}//pickedEntry()


std::string DrfSelect::entryIdFor( const std::shared_ptr<const DetectorPeakResponse> &drf ) const
{
  if( !drf || !m_entries )
    return string();
  const uint64_t hash = drf->hashValue();
  for( const DrfLibrary::Entry &entry : *m_entries )
  {
    if( !entry.hidden && has_hash( entry, hash ) )
      return entry.id;
  }
  return string();
}//entryIdFor(...)


void DrfSelect::pickEntry( const std::string &entryId )
{
  m_userPicked = true;
  m_useWhenLoadedId.clear();  //the user moved on from a kept folder's detector still loading
  if( m_tab != Tab::Library )
    showTab( Tab::Library );
  selectEntry( entryId );
  if( m_narrow )
    setDetailsShown( true );
}//pickEntry(...)


void DrfSelect::setDetailsShown( const bool show )
{
  m_showDetails = (m_narrow && show);
  toggleStyleClass( "ShowDetails", m_showDetails );
  m_detailsHeader->setHidden( !m_showDetails );
  m_detailsBack->setText( WString::tr( (m_tab == Tab::Import) ? "dsl-back-import" : "dsl-back-library" ) );
  if( !m_showDetails && (m_tab == Tab::Library) )
    m_list->scrollToPicked();
  updatePrimaryButton();
  updateAutoPick();
}//setDetailsShown(...)


void DrfSelect::selectEntry( const std::string &entryId )
{
  const bool same = (entryId == m_pickedId);
  m_pickedId = entryId;
  m_list->setPicked( entryId );

  const DrfLibrary::Entry * const entry = pickedEntry();
  if( !same || !m_pickedDrf )
  {
    m_pickedDrf.reset();
    m_pickedError = WString();
    ++m_loadGeneration;

    if( entry && entry->summary )
    {
      if( m_inUse && has_hash( *entry, m_inUse->hashValue() ) )
      {
        m_pickedDrf = m_inUse;
      }else if( entry->kind == DrfLibrary::EntryKind::InMemory )
      {
        m_pickedDrf = entry->inMemoryDrf;
      }else
      {
        // Reading a detector can take a moment (a CeeLo response, a .par grid): do it off the
        //  session thread, keeping only the result for the detector still picked.
        const uint64_t generation = m_loadGeneration;
        const DrfLibrary::Entry entry_copy = *entry;
        const shared_ptr<DataBaseUtils::DbSession> sql = m_sql;
        const WidgetUtils::WidgetHandle self( this );
        struct Loaded
        {
          shared_ptr<const DetectorPeakResponse> drf;
          string error;
        };
        const auto loaded = make_shared<Loaded>();

        AppUtils::run_on_ioservice( [entry_copy, sql, loaded](){
            try
            {
              loaded->drf = DrfLibrary::loadDrf( entry_copy, sql );
            }catch( std::exception &e )
            {
              loaded->error = e.what();
            }
          },
          [self, generation, loaded](){
            DrfSelect * const widget = self.resolve_as<DrfSelect>();
            if( !widget || (widget->m_loadGeneration != generation) )
              return;
            widget->m_pickedDrf = loaded->drf;
            if( !loaded->drf )
            {
              widget->m_pickedError = WString::tr("dsl-load-error").arg( loaded->error );
              widget->m_useWhenLoadedId.clear();
            }
            if( loaded->drf && !widget->m_useWhenLoadedId.empty()
                && (widget->m_useWhenLoadedId == widget->m_pickedId) )
            {
              widget->m_useWhenLoadedId.clear();
              widget->useDetector( make_shared<DetectorPeakResponse>( *loaded->drf ) );
              if( widget->m_closeAfterUse )
              {
                widget->m_doneSignal.emit();  //deletes the dialog
                if( wApp )
                  wApp->triggerUpdate();
                return;
              }
            }
            widget->updatePreview();
            widget->updatePrimaryButton();
            if( wApp )
              wApp->triggerUpdate();
          } );
      }//if( in use ) / in memory / else
    }//if( entry && entry->summary )
  }//if( a different detector )

  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//selectEntry(...)


void DrfSelect::showTab( const Tab tab )
{
  if( m_showDetails )
  {
    m_showDetails = false;
    removeStyleClass( "ShowDetails" );
    m_detailsHeader->hide();
  }
  m_tab = tab;
  m_stack->setCurrentIndex( (tab == Tab::Library) ? 0 : ((tab == Tab::Import) ? 1 : 2) );

  // Locations is a page of the Library
  for( WPushButton *btn : { m_libTab, m_importTab } )
  {
    const bool on = ((btn == m_libTab) == (tab != Tab::Import));
    btn->toggleStyleClass( "Selected", on );
    btn->setAttributeValue( "aria-selected", on ? "true" : "false" );
  }

  m_preview->setHidden( tab == Tab::Locations );
  if( tab == Tab::Locations )
    updateLocationsPage();
  m_preview->setContext( (tab == Tab::Import) ? DrfPreviewPanel::Context::Import
                                              : DrfPreviewPanel::Context::Library );
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//showTab(...)


void DrfSelect::updateInUseStrip()
{
  const bool have = (m_inUse && m_inUse->isValid());
  m_stripInUse->setHidden( !have );
  m_stripNone->setHidden( have );
  if( !have )
    return;

  const shared_ptr<const DrfLibrary::Summary> summary = DrfLibrary::Summary::fromDrf( *m_inUse );
  const DrfLibrary::Entry *entry = nullptr;
  const string id = entryIdFor( m_inUse );
  for( size_t i = 0; m_entries && !id.empty() && (i < m_entries->size()); ++i )
  {
    if( (*m_entries)[i].id == id )
      entry = &(*m_entries)[i];
  }
  const DrfLibrary::Origin origin = entry ? entry->origin : DrfLibrary::originForSource( m_inUse->drfSource() );

  m_stripIcon->setText( WString::fromUTF8( DrfLibraryText::iconHtml( summary.get(), 18.0, 32.0 ) ) );
  m_stripName->setText( WString::fromUTF8( m_inUse->name() ) );
  m_stripMeta->setText( DrfLibraryText::metaLine( *summary, origin ) );
  m_stripJust->setHidden( !m_justChanged );
}//updateInUseStrip()


std::shared_ptr<DetectorPeakResponse> DrfSelect::importCandidate() const
{
  if( m_importStack->currentIndex() == 1 )
    return m_formula->candidate();

  if( m_batch )
  {
    // A copy, so the found detectors stay as they were read
    if( m_found && (m_foundPick < m_found->size()) && (*m_found)[m_foundPick].drf )
      return make_shared<DetectorPeakResponse>( *(*m_found)[m_foundPick].drf );
    return nullptr;
  }

  return m_importWidget->candidate();
}//importCandidate()


std::shared_ptr<DetectorPeakResponse> DrfSelect::importSeed() const
{
  if( m_importStack->currentIndex() == 1 )
    return nullptr;

  if( m_batch )
  {
    if( m_found && (m_foundPick < m_found->size()) )
    {
      const DrfImport::FoundDetector &found = (*m_found)[m_foundPick];
      if( (found.source.status == DrfImport::Status::NeedsCharacterization) && found.source.base )
        return make_shared<DetectorPeakResponse>( *found.source.base );
    }
    return nullptr;
  }

  return m_importWidget->characterizationSeed();
}//importSeed()


namespace
{
  /** The top folder all of a drop's files are in, or empty if they are not all in one. */
  string drop_top_folder( const DrfDropBatch &batch )
  {
    string top;
    for( size_t i = 0; i < batch.files.size(); ++i )
    {
      const string &path = batch.files[i].first;
      const size_t pos = path.find( '/' );
      if( pos == string::npos )
        return string();
      const string first = path.substr( 0, pos );
      if( i == 0 )
        top = first;
      else if( first != top )
        return string();
    }
    return top;
  }//drop_top_folder(...)


  /** Where the Library will list a found detector once its folder is kept at `root`: its GADRAS
   folder, or its main file. */
  string kept_entry_path( const string &root, const string &top_folder, const DrfImport::FoundDetector &found )
  {
    if( found.files.empty() )
      return string();
    string relative = found.files.front();
    if( !top_folder.empty() && (relative.size() > top_folder.size()) )
      relative = relative.substr( top_folder.size() + 1 );
    const bool gadras_folder = (found.files.size() == 2)
                && ((found.kind == DrfImport::FileKind::GadrasEfficiencyCsv)
                    || (found.kind == DrfImport::FileKind::EfficiencyCsv)
                    || (found.kind == DrfImport::FileKind::GadrasDetectorDat));
    const string path = SpecUtils::append_path( root, relative );
    return gadras_folder ? SpecUtils::parent_path( path ) : path;
  }//kept_entry_path(...)


  std::filesystem::path to_fs_path_local( const string &utf8 )
  {
#ifdef _WIN32
    return std::filesystem::path( SpecUtils::convert_from_utf8_to_utf16( utf8 ) );
#else
    return std::filesystem::path( utf8 );
#endif
  }
}//namespace


void DrfSelect::handleBatchDropped( std::shared_ptr<const DrfDropBatch> batch )
{
  if( !batch )
    return;

  m_batch = batch;
  m_found.reset();
  m_foundPick = std::string::npos;
  m_foundFolder = drop_top_folder( *batch );
  const uint64_t generation = ++m_batchGeneration;

  m_importWidget->hide();
  m_formats->hide();
  m_foundList->show();
  m_foundList->setReading( WString::tr("dsl-found-reading").arg( static_cast<int>(batch->files.size()) ) );
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();

  // Reading and grouping the files can take a while (a .par grid takes seconds): off the session thread.
  struct Result
  {
    vector<DrfImport::FoundDetector> found;
    vector<string> skipped;
  };
  const auto result = make_shared<Result>();
  const WidgetUtils::WidgetHandle self( this );

  AppUtils::run_on_ioservice( [batch, result](){
      vector<DrfImport::InputFile> files;
      for( const pair<string,string> &f : batch->files )
      {
        DrfImport::InputFile input;
        input.relativePath = f.first;
        try
        {
          input.data = DrfImport::readFile( f.second );
        }catch( std::exception & )
        {
        }
        files.push_back( input );
      }
      result->found = DrfImport::groupFiles( files, &result->skipped );
    },
    [self, generation, result, batch](){
      DrfSelect * const widget = self.resolve_as<DrfSelect>();
      if( !widget || (widget->m_batchGeneration != generation) )
        return;
      widget->receiveFound( make_shared<const vector<DrfImport::FoundDetector>>( std::move(result->found) ),
                            result->skipped.size() + batch->numSkipped );
      if( wApp )
        wApp->triggerUpdate();
    } );
}//handleBatchDropped(...)


void DrfSelect::receiveFound( std::shared_ptr<const std::vector<DrfImport::FoundDetector>> found,
                              const size_t numSkipped )
{
  m_found = found;
  m_foundPick = std::string::npos;

  vector<DrfFoundList::Item> items;
  size_t first_ready = std::string::npos;
  for( size_t i = 0; found && (i < found->size()); ++i )
  {
    const DrfImport::FoundDetector &d = (*found)[i];
    DrfFoundList::Item item;
    item.name = d.name;

    // "GADRAS Detector.dat · Detector.dat + Efficiency.csv"
    WString format = DrfImportWidget::kindLabel( d.kind );
    string names;
    for( const string &f : d.files )
      names += (names.empty() ? "" : " + ") + SpecUtils::filename( f );
    format += WString::fromUTF8( " · " + names );
    item.format = format;

    const shared_ptr<const DetectorPeakResponse> drf = d.drf ? d.drf : d.source.base;
    shared_ptr<const DrfLibrary::Summary> summary;
    if( drf )
    {
      summary = DrfLibrary::Summary::fromDrf( *drf );
      item.crystal = DrfLibraryText::crystalLabel( *summary );
    }

    switch( d.source.status )
    {
      case DrfImport::Status::Ready:
        if( d.drf )
        {
          item.tag = WString::tr( (summary && (summary->range == DrfLibrary::Range::NearFieldOk))
                                  ? "dsl-found-ready-near" : "dsl-found-ready" );
          if( first_ready == std::string::npos )
            first_ready = i;
        }else
        {
          // Its file reads fine, but no detector could be made with the default options
          item.tag = WString::tr( "dsl-found-needs-input" );
          item.usable = false;
        }
        break;
      case DrfImport::Status::NeedsCharacterization:
        item.tag = WString::tr( "dsl-found-geometry" );
        break;
      case DrfImport::Status::Error:
        item.tag = WString::tr( "dsl-found-error" );
        item.usable = false;
        break;
      default:
        item.tag = WString::tr( "dsl-found-needs-input" );
        item.usable = false;
        break;
    }//switch( d.source.status )

    items.push_back( item );
  }//for( loop over found detectors )

  const size_t num = items.size();
  const WString title = m_foundFolder.empty()
                        ? WString::tr("dsl-found-files").arg( static_cast<int>(m_batch ? m_batch->files.size() : 0) )
                        : WString::fromUTF8( m_foundFolder );
  const WString summary = WString::tr( (num == 0) ? "dsl-found-none" : ((num == 1) ? "dsl-found-1" : "dsl-found-n") )
                            .arg( static_cast<int>(num) );
  WString subtitle;
  if( numSkipped )
    subtitle = WString::tr("dsl-found-skipped").arg( static_cast<int>(numSkipped) );
  if( m_batch && m_batch->numRefused )
  {
    const WString refused = WString::tr("dsl-found-refused").arg( static_cast<int>(m_batch->numRefused) );
    subtitle = subtitle.empty() ? refused : WString::fromUTF8( subtitle.toUTF8() + " · " + refused.toUTF8() );
  }

  bool can_keep = false;
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  can_keep = (num > 1);
#endif

  // In the macOS app a dropped folder can be listed where it is, rather than copied.
  if( m_batch && !m_batch->nativeFolder.empty() )
  {
    can_keep = (num > 0);
    m_foundList->setKeepDescription( WString::tr("dsl-found-keep-desc-link") );
  }else
  {
    // The phone layout has no Locations page to point to
    m_foundList->setKeepDescription( WString::tr( m_narrow ? "dsl-found-keep-desc-phone" : "dsl-found-keep-desc") );
  }

  m_foundList->setFound( title, subtitle, summary, items, can_keep );
  if( first_ready != std::string::npos )
  {
    pickFound( first_ready );
  }else
  {
    updatePreview();
    updatePrimaryButton();
    updateAutoPick();
  }
}//receiveFound(...)


void DrfSelect::pickFound( const size_t index )
{
  if( !m_found || (index >= m_found->size()) )
    return;
  m_foundPick = index;
  m_foundList->setPicked( index );
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//pickFound(...)


void DrfSelect::clearFound()
{
  ++m_batchGeneration;  //drops a read still running
  m_batch.reset();
  m_found.reset();
  m_foundPick = std::string::npos;
  m_foundFolder.clear();
  m_foundList->hide();
  m_importWidget->show();
  m_formats->show();
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//clearFound()


void DrfSelect::keepFoundAndUse()
{
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  if( !m_batch || !m_found || (m_foundPick >= m_found->size()) )
    return;
  const DrfImport::FoundDetector &picked = (*m_found)[m_foundPick];

  // A folder whose own path we know (dropped in the macOS app) is listed where it is.
  if( !m_batch->nativeFolder.empty() )
  {
    const string folder = m_batch->nativeFolder;
    m_useAfterScanRoot = folder;
    m_useAfterScanName = picked.name;
    m_useAfterScanPath = kept_entry_path( folder, m_foundFolder, picked );
    m_useAfterScanFallback = picked.drf ? make_shared<DetectorPeakResponse>( *picked.drf ) : nullptr;
    clearFound();
    m_originChip = -1;
    m_withGeometry = false;
    m_search->setText( "" );
    updateChips();
    updateFilter();
    showTab( Tab::Library );
    if( addUserLocation( folder ) )  //saves it, bookmarks it on macOS, and re-reads the Library
      return;

    // Could not be added: use the detector as dropped (saved under Made by you)
    const shared_ptr<DetectorPeakResponse> fallback = m_useAfterScanFallback;
    m_useAfterScanRoot.clear();
    m_useAfterScanFallback.reset();
    if( fallback )
      useDetector( fallback );
    if( m_closeAfterUse )
      m_doneSignal.emit();  //deletes this dialog
    return;
  }

  // Only the files the detectors were made from; into a new folder of the user-data `drfs` folder.
  set<string> needed;
  for( const DrfImport::FoundDetector &d : *m_found )
    needed.insert( begin(d.files), end(d.files) );

  try
  {
    const string drfs = SpecUtils::append_path( InterSpec::writableDataDirectory(), "drfs" );
    string name = m_foundFolder.empty() ? WString::tr("dsl-found-default-folder").toUTF8() : m_foundFolder;
    for( char &c : name )
    {
      if( string("\\/:?\"<>|*").find( c ) != string::npos )
        c = '_';
    }
    string dest = SpecUtils::append_path( drfs, name );
    for( int i = 2; SpecUtils::is_directory( dest ) || SpecUtils::is_file( dest ); ++i )
      dest = SpecUtils::append_path( drfs, name + " (" + std::to_string(i) + ")" );

    // Copied into a hidden folder (the Library skips names starting with '.'), renamed once
    //  complete - so a copy that fails half-way never shows up as broken detectors.
    const string partial = SpecUtils::append_path( drfs, "." + SpecUtils::filename( dest ) + ".partial" );
    std::error_code ec;
    std::filesystem::remove_all( to_fs_path_local( partial ), ec );
    try
    {
      for( const pair<string,string> &f : m_batch->files )
      {
        if( !needed.count( f.first ) )
          continue;
        string relative = f.first;
        if( !m_foundFolder.empty() && (relative.size() > m_foundFolder.size()) )
          relative = relative.substr( m_foundFolder.size() + 1 );
        const string to = SpecUtils::append_path( partial, relative );
        std::filesystem::create_directories( to_fs_path_local( SpecUtils::parent_path( to ) ) );
        std::filesystem::copy_file( to_fs_path_local( f.second ), to_fs_path_local( to ) );
      }
      std::filesystem::rename( to_fs_path_local( partial ), to_fs_path_local( dest ) );
    }catch( std::exception & )
    {
      std::filesystem::remove_all( to_fs_path_local( partial ), ec );
      throw;
    }

    m_useAfterScanRoot = dest;
    m_useAfterScanPath = kept_entry_path( dest, m_foundFolder, picked );
    m_useAfterScanName = picked.name;
    m_useAfterScanFallback = picked.drf ? make_shared<DetectorPeakResponse>( *picked.drf ) : nullptr;
    showNotice( WString::tr("dsl-found-kept").arg( SpecUtils::filename( dest ) ) );
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dsl-found-keep-error").arg( Wt::Utils::htmlEncode( string(e.what()) ) ), WarningWidget::WarningMsgHigh );
    m_useAfterScanRoot.clear();
    if( picked.drf )
      useDetector( make_shared<DetectorPeakResponse>( *picked.drf ) );
    return;
  }//try / catch

  clearFound();
  m_originChip = -1;
  m_withGeometry = false;
  m_search->setText( "" );
  updateChips();
  updateFilter();
  showTab( Tab::Library );
  startScan();
#endif
}//keepFoundAndUse()


void DrfSelect::useKeptEntry()
{
  if( m_useAfterScanRoot.empty() || !m_entries )
    return;

  const string root = m_useAfterScanRoot;
  const string name = m_useAfterScanName;
  const string expected = m_useAfterScanPath;
  const shared_ptr<DetectorPeakResponse> fallback = m_useAfterScanFallback;
  m_useAfterScanRoot.clear();
  m_useAfterScanName.clear();
  m_useAfterScanPath.clear();
  m_useAfterScanFallback.reset();

  // The entry from the kept folder: at the picked detector's file (and of its name, for a table),
  //  else of its name, else the only one there.
  const auto in_root = [&root]( const string &path ) -> bool {
    if( path.size() == root.size() )
      return path == root;
    return (path.size() > root.size()) && (path.compare( 0, root.size(), root ) == 0)
           && ((path[root.size()] == '/') || (path[root.size()] == '\\'));
  };

  const DrfLibrary::Entry *by_path = nullptr, *by_path_and_name = nullptr, *by_name = nullptr, *only = nullptr;
  size_t num_there = 0;
  for( const DrfLibrary::Entry &entry : *m_entries )
  {
    if( entry.hidden || !entry.summary || !in_root( entry.path ) )
      continue;
    num_there += 1;
    only = &entry;
    const bool same_name = SpecUtils::iequals_ascii( entry.summary->name, name );
    const bool same_path = !expected.empty()
                 && (SpecUtils::lexically_normalize_path( entry.path ) == SpecUtils::lexically_normalize_path( expected ));
    if( same_path && same_name )
      by_path_and_name = &entry;
    if( same_path && !by_path )
      by_path = &entry;
    if( same_name && !by_name )
      by_name = &entry;
  }
  const DrfLibrary::Entry *match = by_path_and_name ? by_path_and_name
                                   : (by_path ? by_path : (by_name ? by_name : ((num_there == 1) ? only : nullptr)));

  if( !match )
  {
    if( fallback )
      useDetector( fallback );
    if( m_closeAfterUse )
      m_doneSignal.emit();  //deletes this dialog
    return;
  }

  const string id = match->id;
  m_userPicked = true;
  m_useWhenLoadedId = id;
  selectEntry( id );
  m_list->scrollToPicked();
  if( m_pickedDrf && (m_useWhenLoadedId == id) )
  {
    m_useWhenLoadedId.clear();
    useDetector( make_shared<DetectorPeakResponse>( *m_pickedDrf ) );
    if( m_closeAfterUse )
      m_doneSignal.emit();  //deletes this dialog
  }
}//useKeptEntry()


void DrfSelect::updatePreview()
{
  if( m_tab == Tab::Import )
  {
    const shared_ptr<DetectorPeakResponse> candidate = importCandidate();
    if( !candidate || !candidate->isValid() )
    {
      m_preview->showNothing();
      return;
    }

    DrfPreviewPanel::Details details;
    details.origin = DrfLibrary::originForSource( candidate->drfSource() );
    details.where = (m_importStack->currentIndex() == 1) ? WString::tr("dsl-where-formula").toUTF8() : string();
    details.whereIsPath = false;
    m_preview->showSummary( DrfLibrary::Summary::fromDrf( *candidate ), details );
    m_preview->setDetector( candidate );
    return;
  }//if( Import tab )

  const DrfLibrary::Entry * const entry = pickedEntry();

  // Before the Library is read: the detector in use.
  if( !entry && !m_entries && m_inUse && m_inUse->isValid() )
  {
    DrfPreviewPanel::Details details;
    details.origin = DrfLibrary::originForSource( m_inUse->drfSource() );
    details.inUse = true;
    m_preview->showSummary( DrfLibrary::Summary::fromDrf( *m_inUse ), details );
    m_preview->setDetector( m_inUse );
    return;
  }

  if( !entry || !entry->summary )
  {
    m_preview->showNothing();
    return;
  }

  DrfPreviewPanel::Details details;
  details.origin = entry->origin;
  details.inUse = (m_inUse && has_hash( *entry, m_inUse->hashValue() ));
  // A database row listed on its own: made by the user, or the record of a detector whose file is
  //  gone - either way only deleting it removes it.
  details.deletable = (entry->kind == DrfLibrary::EntryKind::DbRow) && (entry->dbId >= 0);

  switch( entry->kind )
  {
    case DrfLibrary::EntryKind::DbRow:
      details.where = WString::tr( entry->summary->parentHash ? "dsl-where-copy" : "dsl-where-saved" ).toUTF8();
      details.whereIsPath = false;
      break;
    case DrfLibrary::EntryKind::InMemory:
      details.where = WString::tr( (entry->origin == DrfLibrary::Origin::ThisSpectrum)
                                   ? "dsl-where-spectrum" : "dsl-where-memory" ).toUTF8();
      details.whereIsPath = false;
      break;
    default:
    {
      // Paths inside InterSpec's own data are shown from "data/".
      string path = entry->path;
      try
      {
        const string data_dir = InterSpec::staticDataDirectory();
        if( !data_dir.empty() && SpecUtils::istarts_with( path, data_dir ) )
          path = SpecUtils::append_path( "data", path.substr( data_dir.size() + ((path.size() > data_dir.size()) ? 1 : 0) ) );
      }catch( std::exception & )
      {
      }
      details.where = path;
      break;
    }
  }//switch( entry->kind )

  const size_t index = static_cast<size_t>( entry - &(*m_entries)[0] );
  const int parent = DrfLibrary::parentIndex( *m_entries, index );
  if( (parent >= 0) && (*m_entries)[parent].summary )
    details.parentName = (*m_entries)[parent].summary->name;
  for( size_t i = 0; i < m_entries->size(); ++i )
  {
    if( (i != index) && !(*m_entries)[i].hidden && (DrfLibrary::parentIndex( *m_entries, i ) == static_cast<int>(index)) )
      details.numCopies += 1;
  }

  m_preview->showSummary( entry->summary, details );
  m_preview->setDetector( m_pickedDrf );
  if( !m_pickedError.empty() )
    m_preview->showError( m_pickedError );
}//updatePreview()


void DrfSelect::updatePrimaryButton()
{
  WString label;
  bool enabled = true;

  if( m_tab == Tab::Locations )
  {
    label = WString::tr("dsl-loc-back");
  }else if( m_tab == Tab::Import )
  {
    const shared_ptr<DetectorPeakResponse> candidate = importCandidate();
    if( candidate && candidate->isValid() )
      label = WString::tr("dsl-use-named").arg( short_name( candidate->name() ) );
    else if( importSeed() )
      label = WString::tr("dsl-compute-eff");
    else
    {
      label = WString::tr("dsl-use-detector");
      enabled = false;
    }
  }else
  {
    const DrfLibrary::Entry * const entry = pickedEntry();
    const bool in_use = entry && m_inUse && has_hash( *entry, m_inUse->hashValue() );
    if( entry && entry->summary && (entry->summary->range == DrfLibrary::Range::NeedsCharacterization) )
    {
      label = WString::tr("dsl-compute-eff");
      enabled = !!m_pickedDrf;
    }else if( entry && entry->summary && !in_use )
    {
      label = WString::tr("dsl-use-named").arg( short_name( entry->summary->name ) );
      enabled = entry->summary->hasEfficiency && m_pickedDrf;
    }else if( !entry && !m_inUse )
    {
      label = WString::tr("dsl-use-detector");
      enabled = false;
    }else
    {
      label = WString::tr("dsl-done");
    }
  }//if( Import ) / else

  // On a phone, the main button belongs to the Details screen - and on Import, previews first.
  if( m_narrow )
  {
    bool show = m_showDetails;
    if( !m_showDetails && (m_tab == Tab::Import) && (importCandidate() || importSeed()) )
    {
      label = WString::tr("dsl-preview");
      enabled = true;
      show = true;
    }else if( m_showDetails && (m_tab == Tab::Import) && enabled && importCandidate() )
    {
      label = WString::tr("dsl-import-and-use");
    }
    m_primaryBtn->setHidden( !show );
    if( m_actionBar )
      m_actionBar->setHidden( !show );
  }//if( m_narrow )

  m_primaryBtn->setText( label );
  m_primaryBtn->setEnabled( enabled );
}//updatePrimaryButton()


void DrfSelect::handlePrimary()
{
  if( m_tab == Tab::Locations )
  {
    showTab( Tab::Library );
    return;
  }

  // Phone: see the detector before importing it
  if( m_narrow && !m_showDetails && (m_tab == Tab::Import) )
  {
    if( importCandidate() || importSeed() )
      setDetailsShown( true );
    return;
  }

  if( m_tab == Tab::Import )
  {
    const shared_ptr<DetectorPeakResponse> candidate = importCandidate();
    if( candidate && candidate->isValid() )
    {
      if( m_batch && m_foundList->keep() )
      {
        m_closeAfterUse = m_narrow;
        keepFoundAndUse();
        return;
      }
      useDetector( candidate );
      if( m_narrow )
      {
        // A phone has no room to show it in the Library: done, as when using a Library detector
        m_doneSignal.emit();  //deletes this dialog: nothing after this
        return;
      }
      if( m_batch )
        clearFound();
      showUsedInLibrary();
    }else if( importSeed() )
    {
      openModifyWindow();
    }
    return;
  }//if( Import tab )

  const DrfLibrary::Entry * const entry = pickedEntry();
  const bool in_use = entry && m_inUse && has_hash( *entry, m_inUse->hashValue() );
  if( entry && entry->summary && (entry->summary->range == DrfLibrary::Range::NeedsCharacterization) )
  {
    openModifyWindow();
  }else if( entry && entry->summary && !in_use )
  {
    if( m_pickedDrf && entry->summary->hasEfficiency )
    {
      useDetector( make_shared<DetectorPeakResponse>( *m_pickedDrf ) );
      if( m_narrow )
        m_doneSignal.emit();  //on a phone, using a detector is the end of it; deletes this dialog
    }
  }else
  {
    // "Done": the auto-pick boxes still apply to the detector in use.
    const shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
    const shared_ptr<DetectorPeakResponse> drf = m_inUse;
    const shared_ptr<DataBaseUtils::DbSession> sql = m_sql;
    const Dbo::ptr<InterSpecUser> user = m_interspec->user();
    for( const pair<WCheckBox *,UseDrfPref::UseDrfType> &box
          : { make_pair( m_autoSerial, UseDrfPref::UseDrfType::UseDetectorSerialNumber ),
              make_pair( m_autoModel, UseDrfPref::UseDrfType::UseDetectorModelName ) } )
    {
      if( meas && drf && box.first && !box.first->isHidden() && box.first->isChecked() )
      {
        const UseDrfPref::UseDrfType type = box.second;
        DrfSelect::setUserPrefferedDetector( drf, sql, user, type, meas );  //a small write; session thread
      }
    }//for( the two boxes )

    m_doneSignal.emit();
  }
}//handlePrimary()


void DrfSelect::useDetector( std::shared_ptr<DetectorPeakResponse> drf )
{
  if( !drf || !drf->isValid() )
    return;

  // The detector belongs to the foreground spectrum; without one there is nothing to use it for.
  if( !m_interspec->measurment( SpecUtils::SpectrumType::Foreground ) )
  {
    showNotice( WString::tr("dsl-notice-no-spectrum") );
    return;
  }

  // Already in use (e.g. "Use" pressed twice): no second undo step.
  if( m_inUse && (m_inUse->hashValue() == drf->hashValue()) )
  {
    updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql );
    return;
  }

  // Recorded first, so the Library lists it (with its "Last used") on the rescan below.
  updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql );

  const shared_ptr<SpecMeas> meas = m_interspec->measurment( SpecUtils::SpectrumType::Foreground );
  const shared_ptr<DataBaseUtils::DbSession> sql = m_sql;
  const Dbo::ptr<InterSpecUser> user = m_interspec->user();
  for( const pair<WCheckBox *,UseDrfPref::UseDrfType> &box
        : { make_pair( m_autoSerial, UseDrfPref::UseDrfType::UseDetectorSerialNumber ),
            make_pair( m_autoModel, UseDrfPref::UseDrfType::UseDetectorModelName ) } )
  {
    if( meas && box.first && !box.first->isHidden() && box.first->isChecked() )
    {
      const UseDrfPref::UseDrfType type = box.second;
      DrfSelect::setUserPrefferedDetector( drf, sql, user, type, meas );  //a small write; session thread
    }
  }//for( the two boxes )

  const shared_ptr<DetectorPeakResponse> prev = m_inUse;
  m_justUsedHash = drf->hashValue();
  m_interspec->detectorChanged().emit( drf );  //calls setDetector(...)

  // The undo/redo steps go through InterSpec, not this dialog, which may be gone by then.
  UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
  if( undoManager )
  {
    undoManager->addUndoRedoStep( [prev](){
      InterSpec * const viewer = InterSpec::instance();
      if( viewer )
        viewer->detectorChanged().emit( prev );
    }, [drf](){
      InterSpec * const viewer = InterSpec::instance();
      if( viewer )
        viewer->detectorChanged().emit( drf );
    }, "Change DRF" );
  }//if( undoManager )

  m_notice->hide();
  startScan();
}//useDetector(...)


void DrfSelect::useNoDetector()
{
  if( !m_inUse )
    return;

  const shared_ptr<DetectorPeakResponse> prev = m_inUse;
  m_justUsedHash = 0;
  m_interspec->detectorChanged().emit( nullptr );

  UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
  if( undoManager )
  {
    undoManager->addUndoRedoStep( [prev](){
      InterSpec * const viewer = InterSpec::instance();
      if( viewer )
        viewer->detectorChanged().emit( prev );
    }, [](){
      InterSpec * const viewer = InterSpec::instance();
      if( viewer )
        viewer->detectorChanged().emit( nullptr );
    }, "Change DRF" );
  }//if( undoManager )

  showNotice( WString::tr("dsl-notice-none") );
}//useNoDetector()


void DrfSelect::setDetector( std::shared_ptr<DetectorPeakResponse> det )
{
  m_inUse = det;
  m_justChanged = (det && m_justUsedHash && (det->hashValue() == m_justUsedHash));

  m_preview->setReference( m_inUse );
  updateInUseStrip();

  // The Library may not have the detector in use yet (a new one): it shows up on the next scan.
  m_list->setInUse( entryIdFor( m_inUse ) );
  m_list->rebuild();
  if( m_inUse && entryIdFor( m_inUse ).empty() && !m_scanning )
    startScan();

  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//setDetector(...)


void DrfSelect::showInUseInList()
{
  if( m_tab != Tab::Library )
    showTab( Tab::Library );

  // Clear what could hide it
  m_originChip = -1;
  m_withGeometry = false;
  m_search->setText( "" );
  updateChips();
  updateFilter();

  const string id = entryIdFor( m_inUse );
  if( !id.empty() )
  {
    pickEntry( id );
    m_list->scrollToPicked();
  }
}//showInUseInList()


void DrfSelect::updateChips()
{
  size_t num_all = 0, num_geom = 0;
  map<int,size_t> num_by_origin;
  for( size_t i = 0; m_entries && (i < m_entries->size()); ++i )
  {
    const DrfLibrary::Entry &entry = (*m_entries)[i];
    if( entry.hidden || (!entry.summary && !entry.pending) )
      continue;
    num_all += 1;
    num_by_origin[static_cast<int>(entry.origin)] += 1;
    if( entry.summary && entry.summary->hasGeometry )
      num_geom += 1;
  }//for( loop over entries )

  const auto chip_text = []( const WString &label, const size_t count ) -> WString {
    return WString::fromUTF8( Wt::Utils::htmlEncode( label.toUTF8() )
                              + " <span class=\"DrfChipCount\">" + std::to_string(count) + "</span>" );
  };

  for( size_t i = 0; i < m_chips.size(); ++i )
  {
    const int origin = m_chips[i].first;
    WPushButton * const chip = m_chips[i].second;
    const size_t count = (origin < 0) ? num_all : num_by_origin[origin];
    const bool on = (origin == m_originChip);
    chip->setText( chip_text( WString::tr( ns_chip_keys[i] ), count ) );
    chip->toggleStyleClass( "Pressed", on );
    chip->setAttributeValue( "aria-pressed", on ? "true" : "false" );
    // A source with nothing in it (e.g., no ISOCS here) is not worth a chip, unless picked.
    chip->setHidden( (origin >= 0) && (count == 0) && !on && m_entries );
  }//for( loop over chips )

  m_geoChip->setText( chip_text( WString::tr("dsl-chip-geometry"), num_geom ) );
  m_geoChip->toggleStyleClass( "Pressed", m_withGeometry );
  m_geoChip->setAttributeValue( "aria-pressed", m_withGeometry ? "true" : "false" );

  m_libTab->setText( WString::fromUTF8( Wt::Utils::htmlEncode( WString::tr("dsl-tab-library").toUTF8() )
                     + (m_entries ? (" <span class=\"DrfTabCount\">" + std::to_string(num_all) + "</span>") : string()) ) );
}//updateChips()


void DrfSelect::updateFilter()
{
  DrfLibrary::Filter filter;
  filter.allOrigins = (m_originChip < 0);
  if( m_originChip >= 0 )
    filter.origin = static_cast<DrfLibrary::Origin>( m_originChip );
  filter.withGeometry = m_withGeometry;
  filter.query = m_search->text().toUTF8();

  m_list->setFilter( filter );
  m_list->rebuild();
}//updateFilter()


void DrfSelect::updateAutoPick()
{
  bool show = (m_autoSerial || m_autoModel) && (m_tab != Tab::Locations) && (!m_narrow || m_showDetails);
  if( show && (m_tab == Tab::Import) )
  {
    const shared_ptr<DetectorPeakResponse> candidate = importCandidate();
    show = (candidate && candidate->isValid());
  }else if( show )
  {
    const DrfLibrary::Entry * const entry = pickedEntry();
    show = entry ? (entry->summary && entry->summary->hasEfficiency) : (m_inUse && m_inUse->isValid());
  }

  m_autoPick->setHidden( !show );
}//updateAutoPick()


void DrfSelect::handleImportChanged()
{
  if( m_tab != Tab::Import )
    return;
  updatePreview();
  updatePrimaryButton();
  updateAutoPick();
}//handleImportChanged()


void DrfSelect::showNotice( const Wt::WString &message )
{
  m_noticeText->setText( message );
  m_notice->show();
}//showNotice(...)


void DrfSelect::updateLocationsPage()
{
  if( !m_locationsPage || !m_scanLocations )
    return;

  size_t num_made = 0;
  for( size_t i = 0; m_entries && (i < m_entries->size()); ++i )
  {
    const DrfLibrary::Entry &entry = (*m_entries)[i];
    if( !entry.hidden && (entry.origin == DrfLibrary::Origin::UserMade) )
      num_made += 1;
  }

  m_locationsPage->setLocations( *m_scanLocations,
                                 m_locationStatuses ? *m_locationStatuses : vector<DrfLibrary::LocationStatus>(),
                                 num_made );
}//updateLocationsPage()


void DrfSelect::setLocationShown( const std::string &id, const bool shown )
{
  DrfLibrary::LocationSettings settings = load_location_settings( m_interspec );

  bool is_user = false;
  for( DrfLibrary::Location &loc : settings.userLocations )
  {
    if( loc.id == id )
    {
      loc.hidden = !shown;
      is_user = true;
    }
  }

  if( !is_user )
  {
    if( shown )
      settings.hiddenStandard.erase( id );
    else
      settings.hiddenStandard.insert( id );
  }

  if( save_location_settings( m_interspec, settings ) )
    startScan();
}//setLocationShown(...)


void DrfSelect::rescanLocation( const std::string &id )
{
  for( size_t i = 0; m_scanLocations && (i < m_scanLocations->size()); ++i )
  {
    if( (*m_scanLocations)[i].id == id )
      DrfLibrary::invalidatePath( (*m_scanLocations)[i].path );
  }
  startScan();
}//rescanLocation(...)


void DrfSelect::confirmRemoveLocation( const std::string &id )
{
  const DrfLibrary::LocationSettings settings = load_location_settings( m_interspec );
  const auto pos = std::find_if( begin(settings.userLocations), end(settings.userLocations),
                                 [&id]( const DrfLibrary::Location &loc ){ return loc.id == id; } );
  if( pos == end(settings.userLocations) )
    return;

  const string name = SpecUtils::filename( pos->path );
  SimpleDialog * const dialog = SimpleDialog::make<SimpleDialog>(
                  WString::tr("dsl-loc-remove-title").arg( Wt::Utils::htmlEncode( name ) ),
                  WString::tr("dsl-loc-remove-text") );
  WPushButton * const remove_btn = dialog->addButton( WString::tr("dsl-loc-remove"), WidgetUtils::ButtonRole::Destructive );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );

  // The dialog can outlive this one: find it again when the button is clicked.
  const Wt::Core::observing_ptr<DrfSelect> self( this );
  remove_btn->clicked().connect( [self, id](){
    if( self )
      self->removeLocation( id );
  } );
}//confirmRemoveLocation(...)


void DrfSelect::removeLocation( const std::string &id )
{
  DrfLibrary::LocationSettings settings = load_location_settings( m_interspec );
  const auto pos = std::find_if( begin(settings.userLocations), end(settings.userLocations),
                                 [&id]( const DrfLibrary::Location &loc ){ return loc.id == id; } );
  if( pos == end(settings.userLocations) )
    return;

  const string path = pos->path;
  settings.userLocations.erase( pos );
  if( !save_location_settings( m_interspec, settings ) )
    return;

#if( BUILD_AS_OSX_APP )
  macOsUtils::removeSecurityScopedBookmark( bookmark_key( path ) );
#endif

  DrfLibrary::invalidatePath( path );
  showNotice( WString::tr("dsl-loc-removed").arg( SpecUtils::filename( path ) ) );
  startScan();
}//removeLocation(...)


bool DrfSelect::addUserLocation( const std::string &path )
{
  if( path.empty() || (!SpecUtils::is_directory( path ) && !SpecUtils::is_file( path )) )
    return false;

  DrfLibrary::LocationSettings settings = load_location_settings( m_interspec );
  const string normalized = SpecUtils::lexically_normalize_path( path );
  for( const DrfLibrary::Location &loc : settings.userLocations )
  {
    if( SpecUtils::lexically_normalize_path( loc.path ) == normalized )
    {
      showNotice( WString::tr("dsl-loc-already").arg( path ) );
      startScan();  //harmless, and what a caller waiting on the scan (keepFoundAndUse) needs
      return true;
    }
  }

  DrfLibrary::Location loc;
  loc.id = loc.path = path;
  loc.isFolder = SpecUtils::is_directory( path );
  settings.userLocations.push_back( loc );
  if( !save_location_settings( m_interspec, settings ) )
    return false;

#if( BUILD_AS_OSX_APP )
  // Without the bookmark the (sandboxed) app cannot read the folder after it restarts.
  if( !macOsUtils::createAndStoreSecurityScopedBookmark( path, bookmark_key( path ) ) )
    passMessage( WString::tr("dsl-loc-bookmark-error").arg( Wt::Utils::htmlEncode( path ) ), WarningWidget::WarningMsgMedium );
#endif

  showNotice( WString::tr("dsl-loc-added").arg( SpecUtils::filename( path ) ) );
  startScan();
  return true;
}//addUserLocation(...)


void DrfSelect::startAddLocationFolder()
{
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  SimpleDialog * const dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dsl-loc-add-folder-title") );
  WPushButton * const ok = dialog->addButton( WString::tr("dsl-loc-add"), WidgetUtils::ButtonRole::Affirm );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
  ok->disable();

  DirectorySelector * const selector = dialog->contents()->addNew<DirectorySelector>();
  selector->pathValidityChanged().connect( ok, [ok]( bool valid ){ ok->setDisabled( !valid ); } );

  const Wt::Core::observing_ptr<DrfSelect> self( this );
  ok->clicked().connect( [self, selector](){
    if( self )
      self->addUserLocation( selector->path() );
  } );
#endif
}//startAddLocationFolder()


void DrfSelect::startAddLocationFile()
{
#if( BUILD_AS_OSX_APP )
  // The macOS app is sandboxed: the file must come from the system's picker.
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    const WidgetUtils::WidgetHandle self( this );
    macOsUtils::showFilePicker( WString::tr("dsl-loc-add-file-title").toUTF8(), "", true, false, false,
      [self]( const vector<string> &paths ){
        // Called on the session thread (macOsUtils posts it there)
        DrfSelect * const widget = self.resolve_as<DrfSelect>();
        if( widget && !paths.empty() )
          widget->addUserLocation( paths.front() );
      } );
    return;
  }
#endif

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  SimpleDialog * const dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dsl-loc-add-file-title") );
  WPushButton * const ok = dialog->addButton( WString::tr("dsl-loc-add"), WidgetUtils::ButtonRole::Affirm );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
  ok->disable();

  WContainerWidget * const row = dialog->contents()->addNew<WContainerWidget>();
  row->addStyleClass( "DrfLocFileRow" );
  WLabel * const label = row->addNew<WLabel>( WString::tr("dsl-loc-file-path") );
  WLineEdit * const edit = row->addNew<WLineEdit>();
  edit->setPlaceholderText( WString::tr("dsl-loc-file-placeholder") );
  label->setBuddy( edit );
  WText * const hint = dialog->contents()->addNew<WText>( "", TextFormat::Plain );
  hint->addStyleClass( "DrfFormHint" );
  hint->setInline( false );

  edit->textInput().connect( [edit, ok, hint](){
    const string path = SpecUtils::trim_copy( edit->text().toUTF8() );
    const bool valid = !path.empty() && SpecUtils::is_file( path );
    ok->setDisabled( !valid );
    hint->setText( (valid || path.empty()) ? WString() : WString::tr("dsl-loc-file-not-found") );
  } );

  const Wt::Core::observing_ptr<DrfSelect> self( this );
  ok->clicked().connect( [self, edit](){
    if( self )
      self->addUserLocation( SpecUtils::trim_copy( edit->text().toUTF8() ) );
  } );
#endif
}//startAddLocationFile()


void DrfSelect::locateLocation( const std::string &id )
{
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  WString name = WString::fromUTF8( id );
  for( size_t i = 0; m_scanLocations && (i < m_scanLocations->size()); ++i )
  {
    const DrfLibrary::Location &loc = (*m_scanLocations)[i];
    if( loc.id != id )
      continue;
    switch( loc.role )
    {
      case DrfLibrary::LocationRole::GadrasInstall: name = WString::tr("dsl-loc-gadras"); break;
      case DrfLibrary::LocationRole::IsocsInstall:  name = WString::tr("dsl-loc-isocs");  break;
      case DrfLibrary::LocationRole::AngleInstall:  name = WString::tr("dsl-loc-angle");  break;
      default: break;
    }
  }

  SimpleDialog * const dialog = SimpleDialog::make<SimpleDialog>( WString::tr("dsl-loc-locate-title").arg( Wt::Utils::htmlEncode( name ) ) );
  WPushButton * const ok = dialog->addButton( WString::tr("dsl-loc-use-folder"), WidgetUtils::ButtonRole::Affirm );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );
  ok->disable();

  DirectorySelector * const selector = dialog->contents()->addNew<DirectorySelector>();
  selector->pathValidityChanged().connect( ok, [ok]( bool valid ){ ok->setDisabled( !valid ); } );

  const Wt::Core::observing_ptr<DrfSelect> self( this );
  ok->clicked().connect( [self, selector, id](){
    if( !self )
      return;
    const string path = selector->path();
    DrfLibrary::LocationSettings settings = load_location_settings( self->m_interspec );
    settings.standardPaths[id] = path;
    settings.hiddenStandard.erase( id );
    if( !save_location_settings( self->m_interspec, settings ) )
      return;
#if( BUILD_AS_OSX_APP )
    macOsUtils::createAndStoreSecurityScopedBookmark( path, bookmark_key( path ) );
#endif
    self->startScan();
  } );
#endif
}//locateLocation(...)


void DrfSelect::showQrCode()
{
#if( USE_QR_CODES )
  const shared_ptr<const DetectorPeakResponse> drf = m_preview->detector();
  if( !drf || !drf->isValid() )
  {
    passMessage( WString::tr("ds-no-drf-loaded"), WarningWidget::WarningMsgHigh );
    return;
  }

  try
  {
    const string url = drf->toAppUrlQr();
    SimpleDialog *dialog = QrCode::displayTxtAsQrCode( url, WString::fromUTF8(drf->name()),
                                                       WString::fromUTF8(drf->description()) );

    // The URL carries the detector's shape, but a Monte-Carlo response is far too large for a QR
    //  code - so say so, rather than let the recipient believe they received it.
    if( dialog && drf->ceeloResponse() )
    {
      WText *note = dialog->contents()->addNew<WText>(
                     WString::tr( drf->hasImportedGrid() ? "ds-qr-no-grid-note" : "ds-qr-no-mc-note" ) );
      note->addStyleClass( "DrfQrNoMcNote" );
      note->setInline( false );
    }
  }catch( std::exception &e )
  {
    passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
#endif //USE_QR_CODES
}//showQrCode()


void DrfSelect::confirmDeletePicked()
{
  const DrfLibrary::Entry * const entry = pickedEntry();
  if( !entry || !entry->summary || (entry->kind != DrfLibrary::EntryKind::DbRow) || (entry->dbId < 0) )
    return;

  const long long db_id = entry->dbId;
  const string name = entry->summary->name;

  SimpleDialog * const dialog = SimpleDialog::make<SimpleDialog>(
                     WString::tr("dsl-delete-title").arg( Wt::Utils::htmlEncode( name ) ), WString::tr("dsl-delete-text") );
  WPushButton * const delete_btn = dialog->addButton( WString::tr("dsl-delete"), WidgetUtils::ButtonRole::Destructive );
  dialog->addButton( WString::tr("Cancel"), WidgetUtils::ButtonRole::Dismiss );

  // The dialog can outlive this one: find it again when the button is clicked.
  const Wt::Core::observing_ptr<DrfSelect> self( this );
  delete_btn->clicked().connect( [self, db_id, name](){
    if( self )
      self->deleteEntry( db_id, name );
  } );
}//confirmDeletePicked()


void DrfSelect::deleteEntry( const long long dbId, const std::string &name )
{
  try
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    const long long user_id = m_interspec->user().id();

    const auto rows = m_sql->session()->find<DetectorPeakResponse>()
                        .where( "InterSpecUser_id = ? AND id = ?" ).bind( user_id ).bind( dbId ).resultList();
    for( Dbo::ptr<DetectorPeakResponse> row : rows )
      row.remove();

    // Nor pick it automatically for any spectrum any more.
    const auto prefs = m_sql->session()->find<UseDrfPref>()
                         .where( "InterSpecUser_id = ? AND DrfIndex = ?" ).bind( user_id ).bind( dbId ).resultList();
    for( Dbo::ptr<UseDrfPref> pref : prefs )
      pref.remove();

    transaction.commit();
  }catch( std::exception &e )
  {
    passMessage( WString::tr("dsl-delete-error").arg( Wt::Utils::htmlEncode( name ) ).arg( Wt::Utils::htmlEncode( string(e.what()) ) ), WarningWidget::WarningMsgHigh );
    return;
  }//try / catch

  DrfLibrary::invalidateDbRow( dbId );
  if( pickedEntry() && (pickedEntry()->dbId == dbId) )
  {
    m_pickedId.clear();
    m_pickedDrf.reset();
    m_list->setPicked( "" );
  }
  updatePreview();
  updatePrimaryButton();
  startScan();
}//deleteEntry(...)


void DrfSelect::handleFitFwhmFinished( std::shared_ptr<DetectorPeakResponse> drf )
{
  if( drf )
    updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql );
  m_doneSignal.emit();
}//handleFitFwhmFinished(...)


void DrfSelect::handleMcResponseFinished( std::shared_ptr<DetectorPeakResponse> drf )
{
  if( drf )
    updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql );
  m_doneSignal.emit();
}//handleMcResponseFinished(...)


void DrfSelect::openModifyWindow()
{
  if( m_modifyWindow )
  {
    m_modifyWindow->show();
    m_modifyWindow->centerWindowHeavyHanded();
    return;
  }

  // On the Import tab, Customize is for the file being imported - including a geometry-only one
  //  (e.g., a lone Detector.dat), which is characterized there.
  shared_ptr<DetectorPeakResponse> seed;
  if( m_tab == Tab::Import )
  {
    seed = importCandidate();
    if( !seed )
      seed = importSeed();
  }else if( m_preview->detector() )
  {
    seed = make_shared<DetectorPeakResponse>( *m_preview->detector() );
  }

  if( !seed )
    return;

  m_modifyWindow = AuxWindow::make<DrfModifyWindow>( m_interspec, seed );

  // Hide the Customize window before handling its result (the same order the standalone FWHM / MC
  //  tools use in InterSpec::fwhmFromForegroundWindow / showMcResponseWindow).
  m_modifyWindow->tool()->updatedDrf().connect( m_modifyWindow.get(), &AuxWindow::hide );
  m_modifyWindow->tool()->updatedDrf().connect( this, &DrfSelect::handleModifyFinished );
  m_modifyWindow->finished().connect( this, [this](){
    if( !m_modifyWindow )
      return;

    AuxWindow::deleteAuxWindow( m_modifyWindow.get() );
    assert( !m_modifyWindow );

    UndoRedoManager *undoRedo = m_interspec->undoRedoManager();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      undoRedo->addUndoRedoStep( [](){ reopenDrfSelectModifyWindow(); },
                                 [](){ closeDrfSelectModifyWindow(); },
                                 "Close modify-DRF tool" );
    }//if( undoRedo )
  } );

  UndoRedoManager *undoRedo = m_interspec->undoRedoManager();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    undoRedo->addUndoRedoStep( [](){ closeDrfSelectModifyWindow(); },
                               [](){ reopenDrfSelectModifyWindow(); },
                               "Show modify-DRF tool" );
  }//if( undoRedo )
}//void openModifyWindow()


void DrfSelect::handleModifyFinished( std::shared_ptr<DetectorPeakResponse> drf )
{
  if( !drf )
    return;

  // The customized detector is used straight away, and listed under the one it was made from.
  updateLastUsedTimeOrAddToDb( drf, m_interspec->user().id(), m_sql, true );  //so a removed icon stays removed
  useDetector( drf );
  showUsedInLibrary();
}//handleModifyFinished(...)


void DrfSelect::showUsedInLibrary()
{
  // Once the rescan useDetector(...) started lists it: picked, with nothing filtering it out.
  m_pickInUseAfterScan = true;
  m_scrolledToPicked = false;
  m_originChip = -1;
  m_withGeometry = false;
  m_search->setText( "" );
  updateChips();
  updateFilter();
  if( m_tab != Tab::Library )
    showTab( Tab::Library );
}//showUsedInLibrary()


// ---------------------------------------------------------------------------------------------
//  DrfSelectWindow
// ---------------------------------------------------------------------------------------------

DrfSelectWindow::DrfSelectWindow( InterSpec *viewer )
  : AuxWindow( WString::tr("dsl-window-title"),
              AuxWindowProperties::IsModal
                   | AuxWindowProperties::TabletNotFullScreen
                   | AuxWindowProperties::DisableCollapse
                   | AuxWindowProperties::EnableResize
                   | AuxWindowProperties::SetCloseable ),
    m_edit( nullptr ),
    m_interspec( viewer )
{
  assert( m_interspec );
  if( !m_interspec )
    m_interspec = InterSpec::instance(); //jic
  if( !m_interspec )
    throw runtime_error( "DrfSelectWindow: invalid ptr to InterSpec" );

  addStyleClass( "DrfSelectWindow" );

  shared_ptr<SpecMeas> meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<DetectorPeakResponse> det = meas ? meas->detector() : nullptr;
  {
    auto editOwned = std::make_unique<DrfSelect>( det, m_interspec, this );
    m_edit = editOwned.get();
    contents()->setPadding( 0 );
    contents()->setOverflow( Overflow::Hidden );
    stretcher()->addWidget( std::move(editOwned), 0, 0 );
  }
  stretcher()->setContentsMargins( 0, 0, 0, 0 );

  // Closing in any way applies nothing: a detector is only applied by "Use".
  const auto close = [this](){
    if( m_interspec )
      m_interspec->closeDrfSelectWindow();
    else
      AuxWindow::deleteAuxWindow( this );
  };
  m_edit->done().connect( this, close );
  finished().connect( this, [close]( Wt::DialogCode ){ close(); } );
  rejectWhenEscapePressed();

  show();

  if( !m_isPhone )
  {
    resize( WLength(1120,WLength::Unit::Pixel), WLength(880,WLength::Unit::Pixel) );
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }
}//DrfSelectWindow constructor


DrfSelectWindow::~DrfSelectWindow()
{
}


DrfSelect *DrfSelectWindow::widget()
{
  return m_edit;
}
