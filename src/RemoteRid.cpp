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
#include <string>
#include <sstream>
#include <iostream>

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
#include <cstdlib>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#elif( defined(_WIN32) )
#define WIN32_LEAN_AND_MEAN
#include <stdio.h>
#include <windows.h>
#include <libloaderapi.h>
#endif

#include <boost/process.hpp>
#include <boost/filesystem.hpp>
#if( defined(_WIN32) )
#include <boost/process/windows.hpp>
#endif
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )


#include <Wt/Utils>
#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStackedWidget>
#include <Wt/WRegExpValidator>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RemoteRid.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#endif

using namespace std;
using namespace Wt;


namespace
{
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
/** Runs an executable from the filesystem with the specified arguments.
 
 TODO: currently doesnt handle error super well - mostly because I dont fully trust the return codes of `full-spec` yet, but we can probably improve things a bit as that project matures.
 */
std::string run_external_command( const string &exe, const vector<string> &args )
{
  namespace bp = boost::process;
  
  const auto pp = boost::filesystem::path(exe).parent_path();
  
  bp::ipstream proc_stdout, proc_stderr; //reading pipe-stream
  
  
#ifdef _WIN32
  // Keep from poping up a terminal window (also maybe try `bp::windows::hide`)
  bp::child c( exe, bp::args(args), bp::start_dir(pp), bp::std_out > proc_stdout, bp::std_err > proc_stderr, bp::windows::create_no_window );
#else
  bp::child c( exe, bp::args(args), bp::start_dir(pp), bp::std_out > proc_stdout, bp::std_err > proc_stderr );
#endif
  
  c.wait();

  string output( std::istreambuf_iterator<char>(proc_stdout), {});
  string error( std::istreambuf_iterator<char>(proc_stderr), {});
  
  int result = c.exit_code();
  
  //cout << "Process output: " << output << endl;
  //cout << "Process stdout: " << error << endl;
  
  // Throw exception only if return code is not success, and either there is some error output,
  //  or no stdout
  if( (result != EXIT_SUCCESS) && (error.size() || output.empty()) )
    throw runtime_error( error );
  
  return output;
  
  /*
#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#define WEXITSTATUS
#endif
  
  FILE *stream = popen( cmd.c_str(), "r" );
  //FILE *stream = popen( "echo 'Hello'", "r" );
  
  if( !stream )
    throw runtime_error( "Failed to start command '" + cmd + "'" );
  
  string result;
  try
  {
    char buf[128];
    std::size_t bytesread;
    
    while( (bytesread = fread( buf, sizeof(buf[0]), sizeof(buf), stream)) )
    {
      result += string(buf,bytesread);
    }
  }catch(...)
  {
    pclose( stream );
    throw;
  }
  
  const int plcose_val = pclose( stream );
  const int exitcode = WEXITSTATUS( plcose_val );
  
  cout << "exitcode=" << exitcode << endl;
  
  SpecUtils::trim(result);
  
  return result;
   */
  /*
  vector<char> buffer( 512, '\0' );
  
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe( popen(cmd.c_str(), "r"), pclose );
  if( !pipe )
    throw std::runtime_error("Failed to execute command '" + cmd + "'" );
  
  while( fgets(buffer.data(), buffer.size(), pipe.get() ) != nullptr )
    result += buffer.data();
  
  return result;
   */
}

#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
}//namespace


namespace RestRidImp
{
class ExternalRidWidget;

// Some structs to represent the results from the Full-Spectrum analysis service
struct FullSpecIsotope
{
  string name;
  string type;
  string confidenceStr;
  
  double countRate;
  double confidence;
};//struct FullSpecIsotope


struct FullSpecResults
{
  int code = -1;
  int analysisError = 0;
  string errorMessage;
  vector<string> analysisWarnings;
  string drf;
  double chi2 = -1.0;
  //estimatedDose: in PhysicalUnits (i.e., 1.0E-6*PhysicalUnits::rem/PhysicalUnits::hour)
  double estimatedDose = -1.0;
  double stuffOfInterest = 1.0;
  double alarmBasisDuration = -1.0;
  
  string foregroundDescription;
  string foregroundTitle;
  string isotopeString;
  
  string analysisType;
  
  vector<FullSpecIsotope> isotopes;
};//struct FullSpecResults



FullSpecResults json_to_results( const string &input )
{
  /*
   // Input will look something like:
   {
     "alarmBasisDuration" : 5025.6,
     "analysisError" : 0,
     "analysisWarnings" : ["background synthesized ..."],
     "chi2" : 1.849143743515015,
     "code" : 0,
     "drf" : "Detective-EX200",
     "errorMessage" : null,
     "foregroundDescription" : "17-Apr-2010 17:33:19, 521.6 &gamma; cps, real time: 5119.3 s",
     "foregroundTitle" : "Foreground",
     "isotopeString" : "Ba133 (H)",
     "isotopes" : [{"confidence" : 9.4, "confidenceStr" : "H", "countRate" : 309.1, "name" : "Ba133", "type" : "Industrial"}],
     "stuffOfInterest" : 9.47948169708252
   }
   */
  
  Json::Object result;
  try
  {
    Json::parse( input, result );
  }catch( std::exception &e )
  {
    throw runtime_error( "Error parsing analysis results: <code>" + string(e.what()) + "</code>" );
  }
  
  Json::Value code = result.contains("code") ? result["code"].toNumber() : Json::Value::Null;
  if( code.isNull() )
    throw runtime_error( "Analysis response JSON did not include a 'code' field." );
  
  FullSpecResults answer;
  
  answer.code = code;
  answer.errorMessage = result["errorMessage"].orIfNull("");
  answer.analysisError = result["analysisError"].orIfNull(0);
  answer.drf = result["drf"].orIfNull("");
  answer.chi2 = result["chi2"].orIfNull( -1.0 );
  answer.estimatedDose = result["estimatedDose"].orIfNull( -1.0 ); //estimation from spectrum, in micro-rem/hour
  if( answer.estimatedDose > 0 )
    answer.estimatedDose *= 1.0E-6 * PhysicalUnits::rem / PhysicalUnits::hour;
  answer.stuffOfInterest = result["stuffOfInterest"].orIfNull( -1.0 );
  answer.alarmBasisDuration = result["alarmBasisDuration"].orIfNull( -1.0 );
  answer.foregroundDescription = result["foregroundDescription"].orIfNull("");
  answer.foregroundTitle = result["foregroundTitle"].orIfNull("");
  answer.isotopeString = result["isotopeString"].orIfNull("");
  answer.analysisType = result["analysisType"].orIfNull("AnalysisTypeNotSpecified");
  
  if( result.contains("analysisWarnings") )
  {
    const Json::Array &warnings = result["analysisWarnings"];
    for( const Json::Value &warning : warnings )
      answer.analysisWarnings.push_back( warning.orIfNull( "invalid-warning" ) );
  }//if( result.contains("analysisWarnings") )
  
  if( result.contains("isotopes") )
  {
    const Json::Array &isotopes = result["isotopes"];
    for( size_t i = 0; i < isotopes.size(); ++i )
    {
      const Json::Object &obj = isotopes[i];
      
      FullSpecIsotope iso;
      
      iso.confidence = obj.get("confidence").orIfNull( -1.0 );
      iso.countRate = obj.get("countRate").orIfNull( -1.0 );
      iso.name = obj.get("name").orIfNull( "null" );
      iso.type = obj.get("type").orIfNull( "null" );
      iso.confidenceStr = obj.get("confidenceStr").orIfNull( "null" );
      
      answer.isotopes.push_back( iso );
    }//for( size_t i = 0; i < isotopes.size(); ++i )
  }//if( result.contains("isotopes") )
  
  return answer;
}//FullSpecResults json_to_results(string)


/** A class to stream the current spectrum data and analysis options to the JavaScript FormData
 format expected by the FullSpectrum service.
 
 Having the client-side javascript make a request to this resource, and then use form data from
 its response is a little round-about, but I dont want to require linking Wt against OpenSSL (or
 another TLS library), so we'll have the client-side browser do the requests, and then send the
 results to this server for display.  Also I chose to use this WResource based method, rather than
 just making hundred of kilobytes (or even megabytes) of JavaScript and then calling
 wApp->doJavaScript(...), because I assume performance will be better, and there is probably some
 practical size limit to what wApp->doJavaScript(...) can handle.
 
 We are declaring up here, and implementing below ExternalRidWidget so everything is defined when its
 used.
 */
class RestRidInputResource : public Wt::WResource
{
  Wt::WFlags<RemoteRid::AnaFileOptions> m_flags;
  std::string m_drf;
  Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
  InterSpec *m_interspec;
  
public:
  RestRidInputResource( InterSpec *interspec, WObject* parent = 0 );
  
  void setAnaFlags( Wt::WFlags<RemoteRid::AnaFileOptions> flags );
  void setDrf( const std::string &drf );
  
  virtual ~RestRidInputResource();
  virtual void handleRequest( const Wt::Http::Request &request,
                             Wt::Http::Response &response );
};//class RestRidInputResource : public Wt::WResource





/** This class handles REST communication/facilitation.
 It is inserted into the DOM so signals are exposed, but it does not display any content, and
 its size is set to zero.
 
 We could probably make this a "lighter" widget than a WContainerWidget, but effort to check
 this out probably wasnt worth it.
 */
class RestRidInterface : public Wt::WContainerWidget
{
public:
  JSignal<std::string> m_info_response_signal;
  JSignal<std::string> m_info_error_signal;
  JSignal<std::string> m_analysis_results_signal;
  JSignal<std::string> m_analysis_error_signal;
  
  RestRidInputResource *m_resource;
  
public:
  RestRidInterface( InterSpec *viewer, Wt::WContainerWidget *parent );
  void setAnaFlags( Wt::WFlags<RemoteRid::AnaFileOptions> flags );
  void setDrf( const string &drf );
  void startRestAnalysis( const string ana_service_url );
  void requestRestServiceInfo( const string url );
  void deleteSelfForToast();
  void handleAutoAnaResultSuccess( string response );
  void handleAutoAnaResultFailure( string response );
};//class RestRidInterface





class ExternalRidWidget : public Wt::WContainerWidget
{
public:
  enum class ServiceType
  {
    Rest, Exe
  };
  
protected:
  InterSpec *m_interspec;
  RemoteRid *m_remote_rid;
  ServiceType m_type;
  
  /// m_url is either the REST API URL (minus the "/analysis" or "/info" part), or path to
  /// executable, depending on ServiceType.
  WString m_prev_url;
  WLineEdit *m_url;
  WStackedWidget *m_drf_stack;
  WPushButton *m_retrieve_drfs_btn;
  int m_current_drf_index; // For undo/redo
  WComboBox *m_drf_select;
  WCheckBox *m_onlyDisplayedSamples;
  WPushButton *m_submit;
  WStackedWidget *m_status_stack;
  WText *m_error;
  WText *m_status;
  WText *m_result;
  
  WCheckBox *m_alwaysDoAnalysisCb;
  WCheckBox *m_autoAnaInDialog;

  RestRidInterface *m_rest_interface;
  
  void setMostRecentUrl( const string &url )
  {
    try
    {
      const char *prefname = "";
      switch( m_type )
      {
        case ServiceType::Rest:
          prefname = "ExternalRidUrl";
          break;
          
        case ServiceType::Exe:
          prefname = "ExternalRidExe";
          break;
      }//switch( m_type )
      
      const string prev_urls = InterSpecUser::preferenceValue<string>( prefname, m_interspec );
      vector<string> prev;
      SpecUtils::split( prev, prev_urls, ";" );
      
      if( !prev.empty() && (prev.front() == url) )
        return;
      
      prev.erase( std::remove_if(begin(prev), end(prev), [&url](const string & a){return a == url;}), end(prev) );
      prev.insert( begin(prev), url );
      
      string new_url;
      for( size_t i = 0; i < prev.size(); ++i )
      {
        if( (new_url.size() + 1 + prev[i].size()) >= UserOption::sm_max_value_str_len )
          break;
        new_url += (i ? ";" : "") + prev[i];
      }
      
      InterSpecUser::setPreferenceValue( m_interspec->m_user, prefname, new_url, m_interspec );
    }catch( std::exception &e )
    {
      cerr << "Error in setMostRecentUrl( '" << url << "' ): " << e.what() << endl;
      assert( 0 );
    }// try / catch
  }//void setMostRecentUrl( const string &url )
  
  
public:
  static std::string mostRecentUserUrl( const ServiceType type, InterSpec *interspec )
  {
    const char *prefname = "";
    switch( type )
    {
      case ServiceType::Rest:
        prefname = "ExternalRidUrl";
        break;
        
      case ServiceType::Exe:
        prefname = "ExternalRidExe";
        break;
    }//switch( m_type )

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    // Make option to just copy FullSpectrum into InterSpecs directory for easy access
    auto check_for_default = [type]() -> string {
      if( type == ServiceType::Rest )
        return "";
      
#if( defined(_WIN32) )
      // Assumes Electron build
      string exe_name = "resources/app/FullSpectrum/full-spec.exe";
#elif( APPLE )
      // Assumes app build
      string exe_name = "../Resources/FullSpectrum/full-spec";
#else
      // Assumes Electron build
      string exe_name = "resources/app/FullSpectrum/full-spec";
#endif

      if( AppUtils::locate_file( exe_name, false, 0, true ) )
        return exe_name;
      
      return "";
    };//check_for_default() lamda
#endif // !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT
    
    try
    {
      const string prev_urls = InterSpecUser::preferenceValue<string>( prefname, interspec );
      vector<string> prev;
      SpecUtils::split( prev, prev_urls, ";" );
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
      if( prev.empty() )
        return check_for_default();
#else
      if( prev.empty() )
        return "";
#endif
      
      return prev.front();
    }catch( std::exception &e )
    {
      cerr << "Error in mostRecentUserUrl(): " << e.what() << endl;
      assert( 0 );
    }//try / catch
    
    return "";
  }//std::string mostRecentUserUrl()
  
public:
  ExternalRidWidget( const ServiceType type, RemoteRid *remoterid, InterSpec *interspec, WContainerWidget *parent = 0 )
  : WContainerWidget( parent ),
   m_interspec( interspec ),
   m_remote_rid( remoterid ),
   m_type( type ),
   m_prev_url(),
   m_url( nullptr ),
   m_drf_stack( nullptr ),
   m_retrieve_drfs_btn( nullptr ),
   m_current_drf_index( -1 ),
   m_drf_select( nullptr ),
   m_onlyDisplayedSamples( nullptr ),
   m_submit( nullptr ),
   m_status_stack( nullptr ),
   m_error( nullptr ),
   m_status( nullptr ),
   m_result( nullptr ),
   m_alwaysDoAnalysisCb( nullptr ),
   m_autoAnaInDialog( nullptr ),
   m_rest_interface( nullptr )
  {
    assert( m_remote_rid );
    assert( interspec );
    
    addStyleClass( "ExternalRidWidget" );
    
    WGridLayout *layout = new WGridLayout( this );
    layout->setContentsMargins( 0, 0, 0, 0 );
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    
    const char *url_label = "";
    switch( m_type )
    {
      case ServiceType::Rest:
        url_label = "URL";
        break;
        
      case ServiceType::Exe:
        url_label = "Exe Path";
        break;
    }//switch( m_type )
    
    WLabel *label = new WLabel( url_label );
    layout->addWidget( label, 0, 0, AlignMiddle );
    m_url = new WLineEdit();
    label->setBuddy( m_url );
    m_url->changed().connect( this, &ExternalRidWidget::urlChanged );
    m_url->enterPressed().connect( this, &ExternalRidWidget::urlChanged );
    
    switch( m_type )
    {
      case ServiceType::Rest:
      {
        const char *url_regex = "^(http(s)?:\\/\\/)[\\w.-]+(?:\\.[\\w\\.-]+)*[\\w\\-\\._~:/?#[\\]@!\\$&'\\(\\)\\*\\+,;=.]+$";
        WRegExpValidator *validator = new WRegExpValidator( url_regex, m_url );
        m_url->setValidator( validator );
        m_url->setPlaceholderText( "URL of FullSpectrum Service" );
        break;
      }//case ServiceType::Rest:
        
      case ServiceType::Exe:
      {
        // TODO: need to implement path validator, and actually for macOS and Electron builds implement allowing to select executable!
        m_url->setPlaceholderText( "Path to FullSpectrum executable" );
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )
    
    m_prev_url = WString::fromUTF8( mostRecentUserUrl(m_type,m_interspec) );
    m_url->setText( m_prev_url );
    
    layout->addWidget( m_url, 0, 1, AlignMiddle );
    layout->setColumnStretch( 1, 1 );
  

    m_drf_stack = new WStackedWidget();
    m_drf_stack->addStyleClass( "DrfStack" );
    layout->addWidget( m_drf_stack, 1, 0, 1, 2 );
    
    WContainerWidget *container = new WContainerWidget();
    WGridLayout *sub_layout = new WGridLayout( container );
    m_retrieve_drfs_btn = new WPushButton( "Retrieve DRFs" );
    sub_layout->addWidget( m_retrieve_drfs_btn, 0, 0, AlignCenter | AlignMiddle );
    m_retrieve_drfs_btn->clicked().connect( this, &ExternalRidWidget::requestServiceInfo );
    
    m_drf_stack->addWidget( container );
    
    container = new WContainerWidget();
    sub_layout = new WGridLayout( container );
    label = new WLabel( "DRF:" );
    sub_layout->addWidget( label, 0, 0 );
    m_drf_select = new WComboBox();
    sub_layout->addWidget( m_drf_select, 0, 1 );
    sub_layout->setColumnStretch( 1, 1 );
    m_current_drf_index = -1;
    m_drf_select->activated().connect( this, &ExternalRidWidget::handleUserChangedDrf );
    
    m_drf_stack->addWidget( container );
    m_drf_stack->setCurrentIndex( 0 );
    
    m_status_stack = new WStackedWidget();
    m_status_stack->addStyleClass( "RestRidInfoStack" );
    layout->addWidget( m_status_stack, 2, 0, 1, 2 );
    
    m_error = new WText();
    m_error->addStyleClass( "RestRidError" );
    m_error->setInline( false );
    m_status_stack->addWidget( m_error );
    
    m_status = new WText();
    m_status->addStyleClass( "RestRidStatus" );
    m_status->setInline( false );
    m_status_stack->addWidget( m_status );
    
    m_result = new WText( "After entering URL of web service, please (optionally) retrieve DRFs, and then press &quot;<em>Submit</em>&quot;" );
    m_result->addStyleClass( "RestRidResult" );
    m_result->setInline( false );
    m_status_stack->addWidget( m_result );
    
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_result) );
    
    
    m_onlyDisplayedSamples = new WCheckBox( "Only Displayed Samples" );
    m_onlyDisplayedSamples->setToolTip( "Normally search-mode and portal data are analyzed in their"
                                        " entirety, but selecting this option will cause only the"
                                        " displayed foreground/background to be analyzed." );
    m_onlyDisplayedSamples->addStyleClass( "AlwaysSubmitAna" );
    layout->addWidget( m_onlyDisplayedSamples, 3, 0, 1, 2 );
    auto meas = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    m_onlyDisplayedSamples->setHidden( !meas || !meas->passthrough() );
    
    
    m_submit = new WPushButton( "Submit Analysis" );
    layout->addWidget( m_submit, 4, 0, 1, 2, AlignCenter );
    m_submit->clicked().connect( this, &ExternalRidWidget::submitForAnalysis );
    m_submit->disable();
    
    WContainerWidget *cb_row = new WContainerWidget();
    cb_row->addStyleClass( "RidOptionCbRow" );
    layout->addWidget( cb_row, 5, 0, 1, 2 );
    
    m_alwaysDoAnalysisCb = new WCheckBox( "Always call on spectrum load", cb_row );
    m_alwaysDoAnalysisCb->addStyleClass( "AlwaysSubmitAna" );
    
    m_autoAnaInDialog = new WCheckBox( "Show dialog", cb_row );
    m_autoAnaInDialog->addStyleClass( "ShowDialogOpt" );
    m_autoAnaInDialog->setToolTip( "If checked, a dialog that must be dismissed by the user will be"
                                  " shown with RID results. If not checked, only a"
                                  " notification that will disappear after a short time, will be"
                                  " shown with RID results." );
    m_autoAnaInDialog->hide();
    
    const ExternalRidAuotCallPref pref = RemoteRid::external_rid_call_pref( m_interspec );
    
    switch( m_type )
    {
      case ServiceType::Rest:
      {
        switch( pref )
        {
          case ExternalRidAuotCallPref::DoNotCall:
          case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
          case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
            m_autoAnaInDialog->setChecked( true );
          case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
            m_alwaysDoAnalysisCb->setChecked( true );
            m_autoAnaInDialog->show();
            break;
        }//switch( pref )
        
        break;
      }//case ServiceType::Rest:
        
      case ServiceType::Exe:
      {
        switch( pref )
        {
          case ExternalRidAuotCallPref::DoNotCall:
          case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
          case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
            m_autoAnaInDialog->setChecked( true );
          case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
            m_alwaysDoAnalysisCb->setChecked( true );
            m_autoAnaInDialog->show();
            break;
        }//switch( pref )
        
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )
    
    m_alwaysDoAnalysisCb->checked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    m_alwaysDoAnalysisCb->unChecked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    m_autoAnaInDialog->checked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    m_autoAnaInDialog->unChecked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    
    m_interspec->displayedSpectrumChanged().connect( this, &ExternalRidWidget::displayedSpectrumChanged );
    
    if( m_type == ServiceType::Rest )
    {
      m_rest_interface = new RestRidInterface( m_interspec, nullptr );
      layout->addWidget( m_rest_interface, 6, 0, 1, 2 );
      
      m_rest_interface->m_info_response_signal.connect( this, &ExternalRidWidget::handleInfoResponse );
      m_rest_interface->m_info_error_signal.connect( this, &ExternalRidWidget::handleInfoResponseError );
      m_rest_interface->m_analysis_results_signal.connect( this, &ExternalRidWidget::handleResultResponse );
      m_rest_interface->m_analysis_error_signal.connect( this, &ExternalRidWidget::handleResultResponseError );
    }//if( m_type == ServiceType::Rest )
    
    urlChanged();
  }//ExternalRidWidget
  
  void uncheckAlwaysCallAnalysis()
  {
    m_alwaysDoAnalysisCb->setChecked( false );
    m_autoAnaInDialog->setHidden( true );
  }
  
  void handleAlwaysCallCheckChanged()
  {
    const auto prev_pref_value = RemoteRid::external_rid_call_pref( m_interspec );
    
    ExternalRidAuotCallPref pref_value = ExternalRidAuotCallPref::DoNotCall;
    
    m_autoAnaInDialog->setHidden( !m_alwaysDoAnalysisCb->isChecked() );
    
    if( m_alwaysDoAnalysisCb->isChecked() )
    {
      const bool dialog = m_autoAnaInDialog->isChecked();
      switch( m_type )
      {
        case ServiceType::Rest:
        {
          pref_value = dialog ? ExternalRidAuotCallPref::AlwaysUseRestWithDialog
                              : ExternalRidAuotCallPref::AlwaysUseRestWithToast;
          break;
        }//case ServiceType::Rest:
          
        case ServiceType::Exe:
        {
          pref_value = dialog ? ExternalRidAuotCallPref::AlwaysUseExeWithDialog
                              : ExternalRidAuotCallPref::AlwaysUseExeWithToast;
          break;
        }//case ServiceType::Exe:
      }//switch( m_type )
      
      
#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT )
      m_remote_rid->alwaysCallRestAnaChecked();
#else
      switch( m_type )
      {
        case ServiceType::Rest:
          m_remote_rid->alwaysCallRestAnaChecked();
          break;
        case ServiceType::Exe:
          m_remote_rid->alwaysCallExeAnaChecked();
          break;
      }//switch( m_type )
#endif
    }//if( m_alwaysDoAnalysisCb->isChecked() )
    
    try
    {
      InterSpecUser::setPreferenceValue( m_interspec->m_user, "AlwaysCallExternalRid",
                                        static_cast<int>(pref_value), m_interspec );
    }catch( std::exception & )
    {
      assert(0);
    }
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo_redo = []( const ExternalRidAuotCallPref value ){
        InterSpec *viewer = InterSpec::instance();
        RemoteRid *rid = viewer ? viewer->remoteRid() : nullptr;
        if( !rid )
          return;
      
        RestRidImp::ExternalRidWidget *tool = rid->restRidTool();
        
#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT )
#else
        switch( value )
        {
          case ExternalRidAuotCallPref::DoNotCall:
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
          case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
          case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
            tool = rid->exeRidTool();
            break;
        }//switch( value )
        
        WCheckBox *always_cb = tool ? tool->m_alwaysDoAnalysisCb : nullptr;
        WCheckBox *dialog_cb = tool ? tool->m_autoAnaInDialog : nullptr;
        
        switch( value )
        {
          case ExternalRidAuotCallPref::DoNotCall:
            if( always_cb )
              always_cb->setChecked( false );
            if( dialog_cb )
              dialog_cb->setHidden( true );
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
          case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
            if( always_cb )
              always_cb->setChecked( true );
            if( dialog_cb )
            {
              dialog_cb->setHidden( false );
              dialog_cb->setChecked( false );
            }
            break;
            
          case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
          case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
            if( always_cb )
              always_cb->setChecked( true );
            if( dialog_cb )
            {
              dialog_cb->setHidden( false );
              dialog_cb->setChecked( true );
            }
            break;
        }//switch( value )
        
        if( tool )
          tool->handleAlwaysCallCheckChanged();
#endif //if ANDROID/IOS/WEB / else
      };//undo_redo
      
      undoRedo->addUndoRedoStep( [=](){ undo_redo(prev_pref_value); },
                                [=](){ undo_redo(pref_value); },
                                "Toggle always do External Rid on spectrum load" );
    }//if( undo )
  }//handleAlwaysCallCheckChanged()
  
  
  Wt::WFlags<RemoteRid::AnaFileOptions> anaFlags() const
  {
    Wt::WFlags<RemoteRid::AnaFileOptions> flags;
    if( !m_onlyDisplayedSamples->isHidden() && m_onlyDisplayedSamples->isChecked() )
      flags |= RemoteRid::AnaFileOptions::OnlyDisplayedSearchSamples;
    
    return flags;
  }//anaFlags()
  
protected:
  
  void startRestAnalysis()
  {
    assert( m_rest_interface );
    if( m_rest_interface )
    {
      m_rest_interface->setDrf( drf() );
      m_rest_interface->setAnaFlags( anaFlags() );
      m_rest_interface->startRestAnalysis( m_url->text().toUTF8() );
    }
  }//void startRestAnalysis()
  
  void requestRestServiceInfo()
  {
    assert( m_rest_interface );
    if( m_rest_interface )
      m_rest_interface->requestRestServiceInfo( m_url->text().toUTF8() );
  }//void requestRestServiceInfo()

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  void receiveExeAnalysis( std::shared_ptr<int> rcode, std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )
  {
    assert( rcode && result && m );
    assert( WApplication::instance() );
    
    std::lock_guard<mutex> lock( *m );
    
    const int code = *rcode;
    const string &res = *result;
    if( code == 0 )
    {
      handleResultResponse( res );
    }else
    {
      handleResultResponseError( res );
    }
    
    wApp->triggerUpdate();
  }//void receiveExeAnalysis( std::shared_ptr<string> result )
#endif  //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )

public:
  static void displayAutoRidAnaResult( std::shared_ptr<int> rcode,
                                                std::shared_ptr<string> result,
                                                std::shared_ptr<std::mutex> m )
  {
    assert( rcode && result && m );
    WApplication *app = WApplication::instance();
    InterSpec *interspec = InterSpec::instance();
    assert( app && interspec );
    if( !app || !interspec )
      return;
    
    const ExternalRidAuotCallPref pref = RemoteRid::external_rid_call_pref( InterSpec::instance() );
    
    bool make_dialog = false;
    
    switch( pref )
    {
      case ExternalRidAuotCallPref::DoNotCall:
        assert( 0 );
        return;
        break;
        
      case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
      case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
        make_dialog = false;
        break;
        
      case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
      case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
        make_dialog = true;
        break;
    }//switch( pref )
    
    std::lock_guard<mutex> lock( *m );
    
    string message;
    const int code = *rcode;
    const string &res = *result;
        
    if( code == 0 )
    {
      // `code == 0` just means the http response was okay - there could still be GADRAS error
      vector<string> nuclides;
      vector<pair<string,string>> iso_descrips;
      std::unique_ptr<const FullSpecResults> result_ptr;
      try
      {
        const FullSpecResults result = json_to_results( res );
        result_ptr.reset( new FullSpecResults(result) );
        
        string iso_str = result.isotopeString;
        SpecUtils::ireplace_all( iso_str, "+", ", " );
        
        if( iso_str.empty() )
        {
          if( result.errorMessage.size() )
          {
            result_ptr.reset();
            message = "GADRAS RID Error: " + result.errorMessage;
          }else if( result.code )
          {
            result_ptr.reset();
            message = "GADRAS RID Error code " + to_string(result.code);
          }else
          {
            message = "GADRAS RID: nothing identified";
          }
        }else
        {
          message = "GADRAS RID: " + iso_str;
        }
        
        const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
        assert( db );
        for( const auto &iso : result.isotopes )
        {
          const SandiaDecay::Nuclide *nuc = db->nuclide(iso.name);
          if( nuc )
            nuclides.push_back( nuc->symbol );
          
          // ReferencePhotopeakDisplay::updateOtherNucsDisplay() uses following convention
          //  to decide if a user should be able to click on a result
          if( nuc )
            iso_descrips.emplace_back( nuc->symbol, iso.type );
          else
            iso_descrips.emplace_back( iso.name, "" );
        }
        
        // We should set the results to ReferencePhotopeakDisplay, and then in
        //  ReferencePhotopeakDisplay::handleSpectrumChange clear out the old results
        //  when the file changes.  Also, should modify the auto submital to submit
        //  the file only if all samples are displayed, and otherwise just the
        //  displayed spectra, and do this when the sample numbers change.
        //
        //  TODO: the results we just got might be for the previous file
        ReferencePhotopeakDisplay *reflines = interspec->referenceLinesWidget();
        if( reflines )
          reflines->setExternalRidResults( "GADRAS", iso_descrips );
      }catch( std::exception &e )
      {
        message = "Failed to parse GADRAS RID results.";
        
        cerr << "displayAutoRidAnaResult: Failed to parse results from GADRAS: "
        << e.what() << endl;
      }//try / catch
      
      
      if( make_dialog && result_ptr )
      {
        WStringStream msg_html_strm;
        generateResultHtml( msg_html_strm, *result_ptr );
        const string msg_html = msg_html_strm.str();
          
        wApp->useStyleSheet( "InterSpec_resources/RemoteRid.css" );
        SimpleDialog *dialog = new SimpleDialog( "External RID Results" );
        
        WText *contents = new WText( msg_html, dialog->contents() );
        contents->addStyleClass( "content RestRidResult" );
        contents->setInline( false );
        
        WPushButton *btn = dialog->addButton( "Close" );
        btn->clicked().connect( std::bind( [rcode,result,m](){
          UndoRedoManager *undoRedo = UndoRedoManager::instance();
          if( undoRedo && undoRedo->canAddUndoRedoNow() )
          {
            auto undo = [=](){ displayAutoRidAnaResult( rcode, result, m ); };
            auto redo = [=](){
              InterSpec *interspec = InterSpec::instance();
              if( interspec )
                interspec->programaticallyCloseAutoRemoteRidResultDialog();
            };
            undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close external RID dialog" );
          }
        } ) );//Close button clicked
          
        btn = dialog->addButton( "RID Tool" );
        btn->clicked().connect( std::bind([=](){
          InterSpec *interspec = InterSpec::instance();
          assert( interspec );
          if( !interspec )
            return;
          interspec->createRemoteRidWindow();
          
          UndoRedoManager *undoRedo = UndoRedoManager::instance();
          if( undoRedo && undoRedo->canAddUndoRedoNow() )
          {
            auto undo = [=](){
              InterSpec *interspec = InterSpec::instance();
              if( interspec )
                interspec->deleteRemoteRidWindow();
              displayAutoRidAnaResult( rcode, result, m );
            };
            auto redo = [=](){
              InterSpec *interspec = InterSpec::instance();
              if( interspec )
              {
                interspec->programaticallyCloseAutoRemoteRidResultDialog();
                interspec->createRemoteRidWindow();
              }
            };
            undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close external RID dialog" );
          }
        }));//RID Tool button clicked
        
        
        interspec->setAutoRemoteRidResultDialog( dialog );
        dialog->finished().connect( interspec, &InterSpec::handleAutoRemoteRidResultDialogClose );
        
        UndoRedoManager *undoRedo = UndoRedoManager::instance();
        if( undoRedo && undoRedo->canAddUndoRedoNow() )
        {
          auto undo = [](){
            InterSpec *interspec = InterSpec::instance();
            if( interspec )
              interspec->programaticallyCloseAutoRemoteRidResultDialog();
          };
          auto redo = [rcode,result,m](){
            displayAutoRidAnaResult( rcode, result, m );
          };
          undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Show External RID dialog" );
        }//if( make undo/redo )
      }else
      {
        WStringStream js;
        js << message;
        
        js << "<div class=\"RemoteRidToastButtons\">";
        if( !nuclides.empty() )
        {
          string nucstr;
          for( size_t i = 0; i < nuclides.size(); ++i )
            nucstr += (i ? "," : "") + nuclides[i];
          
          js <<
          "<div onclick=\"Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'showRemoteRidRefLines-" << nucstr << "');"
          "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
          "return false;\">Show Ref. Lines</div>";
        }//if( !nuclides.empty() )
        
        js <<
        "<div onclick=\"Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'openRemoteRidTool');"
        "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
        "return false;\">Open Remote RID</div>";
        
        js << "</div>";
        
        InterSpec *viewer = InterSpec::instance();
        if( viewer )
          viewer->warningWidget()->addMessageUnsafe( js.str(), WarningWidget::WarningMsgExternalRiid, 20000 );
      }//if( make_dialog ) / else
    }else // if( code != 0 )
    {
      // TODO: check for return code, etc., and customize message based on that
      
      string msg = res;
      if( SpecUtils::icontains(msg , "Failed to fetch") //chrome
         || SpecUtils::icontains(msg , "Load failed")   //Safari
         || SpecUtils::icontains(msg , "NetworkError when attempting") ) //Firefox
      {
        msg = "Error calling analysis: incorrect URL or no internet.";
      }
      
      WStringStream js;
      js << msg
         <<
      "<div class=\"RemoteRidToastButtons\">"
      "<div onclick=\"Wt.emit( $('.specviewer').attr('id'),{name:'miscSignal'}, 'openRemoteRidTool');"
      "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
      "return false;\">Open Remote RID</div>"
      "<div onclick=\"Wt.emit( $('.specviewer').attr('id'),{name:'miscSignal'}, 'disableRemoteRid');"
      "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
      "return false;\">Disable Remote RID</div>"
      "</div>";
      
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        viewer->warningWidget()->addMessageUnsafe( js.str(), WarningWidget::WarningMsgExternalRiid, 20000 );
    }//if( code == 0 ) / else
    
    wApp->triggerUpdate();
  }//void displayAutoRidAnaResult( std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )

public:

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  static void startExeAnalysis( string exe_path, ExternalRidWidget *parent, string drf, shared_ptr<SpecUtils::SpecFile> spec_file )
  {
    vector<string> arguments;
    arguments.push_back( "--mode=command-line" );
    arguments.push_back( "--out-format=json" );
    arguments.push_back( "--drf" );
    arguments.push_back( drf );
    
    if( spec_file->num_measurements() < 2 )
      arguments.push_back( "--synthesize-background=1" );
    
    if( !AppUtils::locate_file(exe_path, false, 0, true ) )
    {
      if( parent )
      {
        parent->m_status_stack->setCurrentIndex( parent->m_status_stack->indexOf(parent->m_error) );
        parent->m_error->setText( "Error locating executable: '" + exe_path + "'" );
      }else
      {
        cerr << __func__ << ": Error locating executable: '" + exe_path + "'" << endl;
      }
      
      return;
    }//if( !AppUtils::locate_file(exe_path, false, 0, true ) )
    
    
    const string tmpfilename = SpecUtils::temp_file_name( "interspec_ana_" + wApp->sessionId(),
                                                          SpecUtils::temp_dir() );
    
    if( SpecUtils::is_file(tmpfilename) || SpecUtils::is_directory(tmpfilename) )
    {
      // Shouldnt ever happen
      if( parent )
      {
        parent->m_status_stack->setCurrentIndex( parent->m_status_stack->indexOf(parent->m_error) );
        parent->m_error->setText( "Error creating temp file name." );
      }else
      {
        cerr << __func__ << ": Error creating temp file name." << endl;
      }
      
      return;
    }//if( SpecUtils::is_file(tmpfilename) || SpecUtils::is_directory(tmpfilename) )
    
    {//begin writing temp file
      ofstream tmpfile( tmpfilename.c_str(), ios::out | ios::binary );
      if( !tmpfile.is_open() )
      {
        if( parent )
        {
          parent->m_status_stack->setCurrentIndex( parent->m_status_stack->indexOf(parent->m_error) );
          parent->m_error->setText( "Error opening temp file." );
        }else
        {
          cerr << __func__ << ": Error opening temp file." << endl;
        }
        
        return;
      }//if( !tmpfile.is_open() )
      
      spec_file->write_2012_N42( tmpfile );
    }//end writing temp file
    
    arguments.push_back( tmpfilename );
    
    
    const string appsession = WApplication::instance()->sessionId();
    auto result = std::make_shared<string>();
    auto rcode = std::make_shared<int>(-1);
    auto m = make_shared<mutex>();
    boost::function<void()> doUpdateFcn;
    
    if( parent )
      doUpdateFcn = wApp->bind( boost::bind( &ExternalRidWidget::receiveExeAnalysis, parent, rcode, result, m ) );
    else
      doUpdateFcn = boost::bind( &ExternalRidWidget::displayAutoRidAnaResult, rcode, result, m );
    
    auto commandRunner = [tmpfilename,exe_path,arguments,doUpdateFcn,rcode,result,m,appsession](){
      try
      {
        string results = run_external_command( exe_path, arguments );
        SpecUtils::trim( results );
        
        if( results.empty() )
          throw runtime_error( "No output from running executable." );
        
        std::lock_guard<mutex> lock( *m );
        *rcode = 0;
        *result = results;
      }catch( std::exception &e )
      {
        std::lock_guard<mutex> lock( *m );
        *rcode = -1;
        *result = e.what();
      }
      
      if( !SpecUtils::remove_file( tmpfilename ) )
      {
        cerr << "\n\nFailed to remove temporary file: '" << tmpfilename << "'\n\n";
        assert( 0 );
      }
      
      WServer *server = WServer::instance();
      if( server )
        server->post( appsession, doUpdateFcn );
    };//commandRunner lamda
    
    WServer::instance()->ioService().boost::asio::io_service::post( commandRunner );
  }//void startExeAnalysis()
#endif
  
  void submitForAnalysis()
  {
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_status->setText( "Requesting analysis" );

#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT )
    startRestAnalysis();
#else
    switch( m_type )
    {
      case ServiceType::Rest:
        startRestAnalysis();
        break;
        
      case ServiceType::Exe:
      {
        const string drf_name = drf();
        const string exe_path = m_url->text().toUTF8();
        shared_ptr<SpecUtils::SpecFile> spec_file = RemoteRid::fileForAnalysis(m_interspec, anaFlags());
        assert( spec_file );
        
        if( !spec_file )
        {
          m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
          m_error->setText( "No displayed spectrum" );
          
          return;
        }//if( !spec_file )
        
        startExeAnalysis( exe_path, this, drf_name, spec_file );
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )
#endif //#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT ) / else
  }//void submitForAnalysis()


#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  void receiveExeDrfInfo( std::shared_ptr<int> success, std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )
  {
    assert( success && result && m );
    assert( WApplication::instance() );
    
    std::lock_guard<mutex> lock( *m );
    
    const int successval = *success;
    const string &res = *result;
    if( successval == 0 )
    {
      string json = "{\"Options\":[{\"name\":\"drf\","
                        "\"possibleValues\":" + res + "}]}";
      handleInfoResponse( json );
    }else
    {
      handleInfoResponseError( res );
    }
    
    wApp->triggerUpdate();
  }//void receiveExeDrfInfo( std::shared_ptr<string> result )
  
  
  void requestExeInfo()
  {
    const string exe_path = m_url->text().toUTF8();
    
    const string appsession = WApplication::instance()->sessionId();
    auto result = std::make_shared<string>();
    auto success = std::make_shared<int>(0);
    auto m = make_shared<mutex>();
    auto doUpdateFcn = wApp->bind( boost::bind( &ExternalRidWidget::receiveExeDrfInfo, this, success, result, m ) );
    
    auto commandRunner = [exe_path,doUpdateFcn,result,success,m,appsession](){
      try
      {
        vector<string> args{"--command-line", "--out-format", "json", "--drfs"};
        string drfs = run_external_command( exe_path, args );
        SpecUtils::trim( drfs );
        
        if( drfs.empty() )
          throw runtime_error( "No output from running executable." );
        
        std::lock_guard<mutex> lock( *m );
        *success = 0;
        *result = drfs;
      }catch( std::exception &e )
      {
        std::lock_guard<mutex> lock( *m );
        *success = -1;
        *result = e.what();
      }
      
      WServer *server = WServer::instance();
      if( server )
        server->post( appsession, doUpdateFcn );
    };//commandRunner lamda
    
    WServer::instance()->ioService().boost::asio::io_service::post( commandRunner );
  }//void requestExeInfo()
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  
  void requestServiceInfo()
  {
    m_retrieve_drfs_btn->disable();
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "Requesting service information." );
    m_current_drf_index = -1;
    m_drf_select->clear();

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    switch( m_type )
    {
      case ServiceType::Rest:
        requestRestServiceInfo();
        break;
        
      case ServiceType::Exe:
        requestExeInfo();
        break;
    }//switch( m_type )
#else
    requestRestServiceInfo();
#endif
  }//void requestServiceInfo()
  
  
  void handleInfoResponse( const std::string &msg )
  {
    m_current_drf_index = -1;
    m_drf_select->clear();
    
    try
    {
      /*
       // The message returned will look like:
      {"Options":
       [{"comment":"Optional name of the Detector Response Function.",
         "name":"drf",
         "possibleValues":["auto","Detective-EX","Detective-EX100",...],
         "required":false,
         "type":"Enumerated"
        },{
         "comment": "Optional array of integer...",
         "name":"foregroundSampleNumbers",
         "required":false,
         "type":"Array of integers"
        },{
         "comment":"Optional array of integer...",
         "name":"backgroundSampleNumbers",
         "required":false,
         "type":"Array of integers"
       },{
         "comment":"Synthesize the background...",
         "name":"synthesizeBackground",
         "required":false,
         "type":"Boolean"
       }],
       "comment":"To make an analysis request, you must POST to /v1/Analysis Using multipart/form-data...",
       "versions":{"ApiInterface":"v1","analysis":"GADRAS 19.2.3","compileDate":"Jul  7 2022"}}
*/
      Json::Object result;
      Json::parse( msg, result );
      
      if( !result.contains("Options") )
        throw runtime_error( "No 'Options' component to service info." );
      
      const Json::Array &options = result["Options"];
      
      vector<string> drfs;
      
      for( Wt::Json::Array::const_iterator iter = options.begin();
          iter != options.end(); ++iter )
      {
        const Wt::Json::Object &option = *iter;
        if( !option.contains("name") || !option.contains("possibleValues") )
        {
          cerr << "JSON option does not contain 'name' or 'possibleValues' keys???" << endl;
          continue;
        }
        
        const Json::Value &name = option.get("name");
        if( name != "drf" )
          continue;
        
        const Json::Array &possibleValues = option.get("possibleValues");
        for( size_t i = 0; i < possibleValues.size(); ++i )
        {
          const Json::Value &value = possibleValues[i];
          if( value.type() != Wt::Json::Type::StringType)
             continue;
             
          const WString val = value;
          const string val_utf8 = val.toUTF8();
          
          // We'll limit the DRF name to 64 characters (arbitrary length)
          if( val_utf8.size() < 64 )
            drfs.push_back( val_utf8 );
          else
            cerr << "A DRF was unexpectedly long: '" << val_utf8 << "' and wont be used." << endl;
        }//for( loop over DRF
        
        // We dont need any more values
        break;
      }//for( iterate over JSON entries );
      
      
      string infotxt;
      if( result.contains("versions") )
      {
        const Wt::Json::Object &versions = result.get("versions");
        for( const string &name : versions.names() )
        {
          const WString val = versions.get(name).toString();
          infotxt += name + ": " + val.toUTF8() + "<br />";
        }
      }//if( result.contains("versions") )
      
      // We'll limit the number of DRFs, to 75 (arbitrary)
      if( drfs.size() > 75 )
      {
        cerr << "There were " << drfs.size() << " DRFs - only using the firs 75." << endl;
        drfs.resize( 75 );
      }
      
      for( const string &v : drfs )
      {
        //cout << "DRF: " << v << endl;
        m_drf_select->addItem( v );
      }
      
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
      m_status->setText( infotxt );
      m_drf_stack->setCurrentIndex( 1 );
    }catch( std::exception &e )
    {
      m_current_drf_index = -1;
      m_drf_select->clear();
      m_retrieve_drfs_btn->disable();
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
      m_result->setText( "" );
      m_status->setText( "" );
      m_error->setText( "Error parsing service information results: <code>" + string(e.what()) + "</code>" );
    }//try / catch
  }//void handleInfoResponse()
  
  
  void handleInfoResponseError( std::string msg )
  {
    m_current_drf_index = -1;
    m_drf_select->clear();
    m_drf_stack->setCurrentIndex( 0 );
    m_retrieve_drfs_btn->disable();
    m_result->setText( "" );
    m_status->setText( "" );
    
    if( SpecUtils::icontains( msg , "Failed to fetch") //chrome
       || SpecUtils::icontains( msg , "Load failed")   //Safari
       || SpecUtils::icontains( msg , "NetworkError when attempting") ) //Firefox
    {
      msg = "incorrect URL or no internet.";
    }
    
    m_error->setText( "Error retrieving service information: <code>" + msg + "</code>" );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
  }//void handleInfoResponseError()

  
  static void generateResultHtml( WStringStream &rslttxt, const FullSpecResults &results )
  {
    const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    rslttxt << "<table class=\"ResultTable\"><tbody>\n"
    << "\t<tr>"
    << "\t\t<th>Nuclide</th>\n"
    << "\t\t<th>Confidence</th>\n"
    << "\t\t<th>Category</th>\n";
    
    if( SpecUtils::iequals_ascii(results.analysisType, "Simple")
       || SpecUtils::icontains(results.analysisType, "NotSpecified")
       || results.analysisType.empty() )
    {
      rslttxt << "\t\t<th>Count Rate</th>\n";
    }else
    {
      rslttxt << "\t\t<th>Max CPS</th>\n";
    }
    
    rslttxt << "\t</tr>";
    
    if( results.isotopes.empty() )
    {
      rslttxt << "\t<tr>"
      << "\t\t<td colspan=\"4\" style=\"text-align: center; vertical-align: middle;\">None Found</td>\n\t</tr>\n";
    }else
    {
      const size_t nres = results.isotopes.size();
      
      for( const FullSpecIsotope &res_iso : results.isotopes )
      {
        const string &iso = res_iso.name;
        string type = res_iso.type;
        string conf = res_iso.confidenceStr;
        const float count_rate = res_iso.countRate;
        
        if( conf == "H" )
          conf = "High";
        else if( conf == "F" )
          conf = "Fair";
        else if( conf == "L" )
          conf = "Low";
        else
          cerr << "Unknown confidence '" << conf << "' for nuclide " << iso << endl;
        
        char buffer[64] = { '\0' };
        if( count_rate < 1.0E-7 )  //FLT_EPSILON is usually 1.19209e-07
          snprintf( buffer, sizeof(buffer), "--" );
        else
          snprintf( buffer, sizeof(buffer), "%.4g", count_rate );
        
        rslttxt << "\t<tr>\n"
        <<"\t\t<td>" << Wt::Utils::htmlEncode(iso) << "</td>\n"
        <<"\t\t<td>" << Wt::Utils::htmlEncode(conf) << "</td>\n"
        <<"\t\t<td>" << Wt::Utils::htmlEncode(type) << "</td>\n"
        <<"\t\t<td>" << std::string(buffer) << "</td>\n"
        <<"\t</tr>\n";
      }//for( size_t i = 0; i < nres; ++i )
    }//if( output.isotope_names.empty() ) / else
    
    rslttxt << "</tbody></table>";
    
    if( results.code != 0 )
      rslttxt << "<div class=\"ResultChi2\">Result Code: " << results.code << "</div>";
    
    string dose_str;
    if( results.estimatedDose > 0 )
    {
      dose_str = "&gamma; dose="
                 + PhysicalUnits::printToBestEquivalentDoseRateUnits( results.estimatedDose, 1, useBq )
                 + ", ";
    }
    
    if( (results.chi2 > 0.0001f) && (results.isotopes.size() > 0) )
    {
      char buffer[64] = { '\0' };
      snprintf( buffer, sizeof(buffer), "%.2f", results.chi2 );
      rslttxt << "<div class=\"ResultChi2\">"
              << dose_str
              << "&chi;<sup>2</sup>=" << std::string(buffer)
              << ", DRF:&nbsp;" << Wt::Utils::htmlEncode(results.drf)
              << "</div>";
    }else
    {
      rslttxt << "<div class=\"ResultChi2\">"
              << dose_str
              << "DRF:&nbsp;" << Wt::Utils::htmlEncode(results.drf)
              << "</div>";
    }//if( output.chi_sqr > 0.0f )
    
    if( !results.analysisWarnings.empty() )
    {
      rslttxt << "<div class=\"AnaWarnings\">\n";
      for( const string &warning : results.analysisWarnings )
        rslttxt << "<div>" << Wt::Utils::htmlEncode(warning) << "</div>\n";
      rslttxt << "</div>\n";
    }//if( result.contains("analysisWarnings") )
  }//void generateResultHtml( WStringStream &rslttxt, const FullSpecResults &results )
  
  
  void handleResultResponse( const std::string &msg )
  {
    m_submit->enable();
    
    FullSpecResults results;
    try
    {
      results = json_to_results( msg );
    }catch( std::exception &e )
    {
      m_error->setText( "Error parsing analysis results: <code>" + string(e.what()) + "</code>" );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
    }
    
    try
    {
      if( !results.errorMessage.empty() )
        throw runtime_error( "Analysis service returned error message: " + results.errorMessage );
      
      WStringStream rslttxt;
      rslttxt << "<div class=\"ResultLabel\">Results:</div>";
      generateResultHtml( rslttxt, results );
      
      m_result->setText( rslttxt.str() );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_result) );
    }catch( std::exception &e )
    {
      m_error->setText( e.what() );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
    }
  }//void handleResultResponse()
  
  
  void handleResultResponseError( std::string msg )
  {
    m_submit->enable();
    m_result->setText( "" );
    m_status->setText( "" );
    
    if( SpecUtils::icontains( msg , "Failed to fetch") //chrome
       || SpecUtils::icontains( msg , "Load failed")   //Safari
       || SpecUtils::icontains( msg , "NetworkError when attempting") ) //Firefox
    {
      msg = "incorrect URL or no internet.";
    }
    
    m_error->setText( msg );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
  }//void handleResultResponseError()
  
  ServiceType serviceType() const
  {
    return m_type;
  }
  
  string url() const
  {
    return m_url->text().toUTF8();
  }
  
  bool alwaysAutoCallAnalysis() const
  {
    return m_alwaysDoAnalysisCb->isChecked();
  }
  
  bool showDialogForAutoCalledResults() const
  {
    if( !m_alwaysDoAnalysisCb->isChecked() )
      return false;
    return m_autoAnaInDialog->isChecked();
  }
  
  bool isValidUrl()
  {
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    switch( m_type )
    {
      case ServiceType::Rest:
      {
        switch( m_url->validate() )
        {
          case Wt::WValidator::Invalid:
          case Wt::WValidator::InvalidEmpty:
            return false;
            
          case Wt::WValidator::Valid:
            return true;
        }//switch( m_url->validate() )
        
        break;
      }//case ServiceType::Rest:
        
        
      case ServiceType::Exe:
      {
        string url = m_url->text().toUTF8();
        try
        {
          if( !AppUtils::locate_file(url, false, 0, true) )
            return false;
          
          const boost::filesystem::file_status status = boost::filesystem::status(url);
          if( !boost::filesystem::exists(status) )
          {
            assert( 0 );
            return false;
          }
          
          const boost::filesystem::perms p = status.permissions();
          
          // We'll check if _anyone_ can execute the file, and if so return true - we are not
          //  checking if _we_ can execute the file, which would require determining if _we_
          //  or our group can execute
          if( (p & boost::filesystem::perms::owner_exe)
             || (p & boost::filesystem::perms::group_exe)
             || (p & boost::filesystem::perms::others_all)
             )
          {
            return true;
          }
          
          return false;
        }catch( std::exception & )
        {
          return false;
        }

        assert( 0 );
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )
#else
    switch( m_url->validate() )
    {
      case Wt::WValidator::Invalid:
      case Wt::WValidator::InvalidEmpty:
        return false;

      case Wt::WValidator::Valid:
        return true;
    }//switch( m_url->validate() )
#endif

    assert( 0 );
    return false;
  }//bool isValidUrl()
  
  
  void urlChanged()
  {
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "" );
    m_current_drf_index = -1;
    m_drf_select->clear();
    m_drf_stack->setCurrentIndex( 0 );
      
    const WString new_url = m_url->text();
    const WString old_url = m_prev_url;
    m_prev_url = new_url;
    
    auto spec = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    
    const bool validUrl = isValidUrl();
    const char *errmsg = "";
    switch( m_type )
    {
      case ServiceType::Rest:
      {
        errmsg = "Invalid URL";
        break;
      }//case ServiceType::Rest:
        
      case ServiceType::Exe:
      {
        errmsg = "Invalid path to executable.";
        
        const char *invalid_style_class = "Wt-invalid";
        if( validUrl )
          m_url->removeStyleClass(invalid_style_class);
        else if( !m_url->hasStyleClass(invalid_style_class) )
          m_url->addStyleClass(invalid_style_class);
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )
    
    
    if( validUrl )
    {
      const string url = m_url->text().toUTF8();
      
      m_submit->setHidden( !spec || url.empty() );
      m_retrieve_drfs_btn->setHidden( !spec || url.empty() );
      m_retrieve_drfs_btn->setEnabled( spec && !url.empty() );
      
      if( !url.empty() )
        setMostRecentUrl( url );
      
      m_submit->enable();
    }else
    {
      m_submit->hide();
      m_retrieve_drfs_btn->hide();
      m_retrieve_drfs_btn->disable();
      m_status->setText( "Invalid URL" );
      m_submit->disable();
    }//if( !validUrl ) / else
    
    if( m_rest_interface )
    {
      m_rest_interface->setAnaFlags( anaFlags() );
      m_rest_interface->setDrf( drf() );
    }
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( (new_url != old_url) && undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      const bool isexe = (m_type == ServiceType::Exe);
      
      auto undo_redo = [new_url,old_url,isexe]( const bool is_undo ){
        InterSpec *viewer = InterSpec::instance();
        RemoteRid *rid = viewer ? viewer->remoteRid() : nullptr;
        
        RestRidImp::ExternalRidWidget *tool = rid->restRidTool();
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
        if( isexe )
          tool = rid->exeRidTool();
#endif
        if( !tool )
          return;
        
        const WString &url = is_undo ? old_url : new_url;
        tool->m_url->setText( url );
        tool->setMostRecentUrl( url.toUTF8() );
        tool->urlChanged();
      };//undo_redo
      
      auto undo = [undo_redo](){ undo_redo(true); };
      auto redo = [undo_redo](){ undo_redo(false); };
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Update External RID URL" );
    }//if( undo )
  }//urlChanged()
  
  
  void handleUserChangedDrf()
  {
    if( m_rest_interface )
    {
      m_rest_interface->setAnaFlags( anaFlags() );
      m_rest_interface->setDrf( drf() );
    }//if( m_rest_interface )
    
    const int new_drf_index = m_drf_select->currentIndex();
    const int prev_drf_index = m_current_drf_index;
    m_current_drf_index = new_drf_index;
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( (new_drf_index != prev_drf_index) && undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      const bool isexe = (m_type == ServiceType::Exe);
      
      auto undo_redo = [new_drf_index,prev_drf_index,isexe]( const bool is_undo ){
        InterSpec *viewer = InterSpec::instance();
        RemoteRid *rid = viewer ? viewer->remoteRid() : nullptr;
        
        RestRidImp::ExternalRidWidget *tool = rid->restRidTool();
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
        if( isexe )
          tool = rid->exeRidTool();
#endif
        if( !tool )
          return;
        
        const int nitems = tool->m_drf_select->count();
        const int index = is_undo ? prev_drf_index : new_drf_index;
        if( index < nitems )
          tool->m_drf_select->setCurrentIndex( index );
      };
      
      auto undo = [undo_redo](){ undo_redo(true); };
      auto redo = [undo_redo](){ undo_redo(false); };
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change External RID DRF" );
    }//if( undo )
  }//void handleUserChangedDrf()
  
  
  void displayedSpectrumChanged()
  {
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "" );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    
    string url = m_url->text().toUTF8();
    if( !isValidUrl() )
      url = "";
    
    auto spec = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !spec || url.empty() )
    {
      m_submit->setHidden( true );
      m_onlyDisplayedSamples->setHidden( true );
      
      if( !spec )
        m_status->setText( "No spectrum displayed to submit." );
      else if( url.empty() )
        m_status->setText( "Entered URL looks to be invalid." );
    }else
    {
      m_submit->setHidden( false );
      m_status->setText( "Click &quot;Submit Analysis&quot; to get RIID results." );
      
      auto m = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
      m_onlyDisplayedSamples->setHidden( !m || !m->passthrough() );
    }
    
    if( m_rest_interface )
    {
      m_rest_interface->setAnaFlags( anaFlags() );
      m_rest_interface->setDrf( drf() );
    }
  }//void displayedSpectrumChanged()
  
  
public:
  string drf() const
  {
    if( (m_drf_select->currentIndex() < 0) || !m_drf_select->count() )
      return "auto";
    
    return m_drf_select->currentText().toUTF8();
  }//string drf() const
  
};//ExternalRidWidget



RestRidInterface::RestRidInterface( InterSpec *viewer, Wt::WContainerWidget *parent )
: Wt::WContainerWidget( parent ),
m_info_response_signal( this, "info_response", false ),
m_info_error_signal( this, "info_error", false ),
m_analysis_results_signal( this, "analysis_results", false ),
m_analysis_error_signal( this, "analysis_error", false ),
m_resource( new RestRidInputResource(viewer, this) )
{
  setPositionScheme( PositionScheme::Absolute ); //or Fixed
  setWidth( 0 );
  setHeight( 0 );
}


void RestRidInterface::setAnaFlags( Wt::WFlags<RemoteRid::AnaFileOptions> flags )
{
  m_resource->setAnaFlags( flags );
}

void RestRidInterface::setDrf( const string &drf )
{
  m_resource->setDrf( drf );
}

void RestRidInterface::startRestAnalysis( const string ana_service_url )
{
  const string &resource_url = m_resource->url();
  
  // Note: the JS is defined in the C++ to avoid having a separate JS file that will end up being
  // packaged in builds that wont use it (i.e., most InterSpec builds); I dont like having various
  // `fetch` commands around that can call off of localhost; if nothing else this is a bad look.
  WStringStream js;
  js <<
  "\n(function(){"
  "let sendAnaRequest = function(formData){\n"
  "  fetch('" << ana_service_url << "/analysis', {method: 'POST', body: formData})\n"
  "    .then(response => {\n"
  "       if( !response.ok ) {\n"
  "         throw new Error('Server returned status ' + response.status);\n"
  "       }\n"
  // We could require server to declare content JSON, but we wont for the moment
  //"       const contentType = response.headers.get('content-type');\n"
  //"       if( !contentType || !contentType.includes('application/json') ){\n"
  //"         throw new TypeError('Server did not return JSON');\n"
  //"       }\n"
  "       return response.json();\n"
  "   })\n"
  "    .then(json_results => {\n"
  //"      console.log( 'Got JSON analysis response:', json_results );\n"
  "      const result_str = JSON.stringify(json_results);\n"
  "      " << m_analysis_results_signal.createCall("result_str") << "\n"
  "    })\n"
  "  .catch( error => {\n"
  "    console.error( 'Error retrieving results from server:', error);\n"
  "    const error_str = 'Error retrieving results from server: ' + error.toString();\n"
  "  " << m_analysis_error_signal.createCall( "error_str" ) << ";\n"
  "  });\n"
  "};\n";
  
  js <<
  "\n"
  "fetch('" << resource_url << "', {method: 'POST'})\n"
  "  .then( response => response.formData() )\n"
  "  .then( data => {\n"
  //"  console.log( 'Got FormData wt app!' );\n"
  "    sendAnaRequest(data);\n"
  "  })\n"
  "  .catch( error => {\n"
  "    console.error( 'Error getting analysis input:', error);\n"
  "    const error_str = 'Error getting analysis input:' + error.toString();\n"
  "    " << m_analysis_error_signal.createCall( "error_str" ) << ";\n"
  "});\n})();";
  
  doJavaScript( js.str() );
}//void startRestAnalysis()

void RestRidInterface::requestRestServiceInfo( const string url )
{
  WStringStream js;
  
  js <<
  "fetch('" << url << "/info', {method: 'POST'})\n"
  "  .then( response => response.json() )\n"
  "  .then( data => {\n"
//            << "console.log( 'Got service info data:', data );"
  "    const data_str = JSON.stringify(data);\n"
  "    " << m_info_response_signal.createCall( "data_str" ) << ";\n"
  "  })\n"
  "  .catch( error => {\n"
  "  console.error( 'Error retrieving info:', error);\n"
  "  const error_str = error.toString();\n"
  "  " << m_info_error_signal.createCall( "error_str" ) << ";\n"
  "});";
  
  doJavaScript( js.str() );
}//void requestRestServiceInfo()



void RestRidInterface::deleteSelfForToast()
{
  // Function called after a RestRidInterface object has served its purpose in getting a
  //  result back (or an error) for a "Toast" notification analysis (e.g., RID ana performed when
  //  spectrum is changed in InterSpec).
  //
  // We could leave the RestRidInterface in the DOM, and just update its ana flags and URL, but
  //  for the moment, we'll just delete the widget after its intended use... although I *think*
  //  we should always get a call-back for success or failure is the case, I'm not 100% certain
  //  ... need to do some more testing
  
  const string rest_id = id();
  RestRidInterface *instance = this;
  
  // I think its fine to just call `delete this`, but we'll do a little delay, "just because"
  WServer::instance()->schedule( 1000, wApp->sessionId(), [rest_id, instance](){
    assert( wApp );
    WWidget *w = wApp->domRoot()->findById( rest_id );
    RestRidInterface *ww = dynamic_cast<RestRidInterface *>( w );
    
    // IF the user did a "Clear Session..." the widget could already be deleted
    if( ww && (ww == instance) )
      delete instance;
    else
      cerr << "When deleting RestRidInterface, got w = " << w << ", ww = " << ww << ", instance = " << instance << endl;
  } );
}//deleteSelfForToast()

void RestRidInterface::handleAutoAnaResultSuccess( string response )
{
  auto result = std::make_shared<string>( response );
  auto rcode = std::make_shared<int>( 0 );
  auto m = make_shared<mutex>();
  
  RestRidImp::ExternalRidWidget::displayAutoRidAnaResult( rcode, result, m );
  
  deleteSelfForToast();
}//handleAutoAnaResultSuccess(...)

void RestRidInterface::handleAutoAnaResultFailure( string response )
{
  auto result = std::make_shared<string>( response );
  auto rcode = std::make_shared<int>( -1 );
  auto m = make_shared<mutex>();
  
  RestRidImp::ExternalRidWidget::displayAutoRidAnaResult( rcode, result, m );
  
  deleteSelfForToast();
}//handleAutoAnaResultFailure(...)



RestRidInputResource::RestRidInputResource( InterSpec *interspec, WObject* parent )
: WResource( parent ),
m_flags( 0 ),
m_drf( "auto" ),
m_app( WApplication::instance() ),
m_interspec( interspec )
{
  assert( m_app && m_interspec );
}

RestRidInputResource::~RestRidInputResource()
{
  beingDeleted();
}


void RestRidInputResource::setAnaFlags( Wt::WFlags<RemoteRid::AnaFileOptions> flags )
{
  m_flags = flags;
}//void setAnaFlags( Wt::WFlags<RemoteRid::AnaFileOptions> flags )


void RestRidInputResource::setDrf( const std::string &drf )
{
  m_drf = drf;
}

void RestRidInputResource::handleRequest( const Wt::Http::Request &request,
                                         Wt::Http::Response &response )
{
  WApplication::UpdateLock lock( m_app );
  
  if( !lock )
  {
    log("error") << "Failed to WApplication::UpdateLock in RestRidInputResource.";
    
    response.out() << "Error grabbing application lock to form RestRidInputResource resource; please report to InterSpec@sandia.gov.";
    response.setStatus(500);
    assert( 0 );
    
    return;
  }//if( !lock )
  
  // Generate some JSON then
  shared_ptr<SpecUtils::SpecFile> spec_file = RemoteRid::fileForAnalysis(m_interspec, m_flags);
  assert( spec_file );
  
  string options = "{\"drf\": \"" + m_drf + "\"";
  if( spec_file && (spec_file->num_measurements() < 2) )
    options += ", \"synthesizeBackground\": true";
  options += "}";
  
  const string boundary = "InterSpec_MultipartBoundary_InterSpec";
  
  string n42_content;
  if( spec_file )
  {
    stringstream n42strm;
    spec_file->write_2012_N42( n42strm );
    n42_content = n42strm.str();
  }else
  {
    cerr << "RestRidInputResource::handleRequest() logic error: invalid SpecUtils::SpecFile\n";
    n42_content = "No Spectrum";
  }
  
  string body =
  "--" + boundary
  + "\r\nContent-Disposition: form-data; name=\"options\""
  + "\r\nContent-type: application/octet-stream"
  + "\r\n\r\n" + options + "\r\n"
  + "--" + boundary
  
  + "\r\nContent-Disposition: form-data; name=\"foreground\"; filename=\"spectrum.n42\""
  + "\r\nContent-type: application/octet-stream"
  + "\r\n\r\n" + n42_content + "\r\n"
  += "--" + boundary + "--\r\n";
  
  response.setMimeType( "multipart/form-data; boundary=" + boundary );
  response.setContentLength( body.size() );
  
  //cout << "Body:\n\n" << body << "\n\n\n" << endl;
  
  response.out() << body;
}//handleRequest(...)

}//namespace



ExternalRidAuotCallPref RemoteRid::external_rid_call_pref( InterSpec *viewer )
{
  ExternalRidAuotCallPref answer = ExternalRidAuotCallPref::DoNotCall;
  try
  {
    if( !viewer )
      throw runtime_error( "Invalid InterSpec instance" );
    
    const int pref_value = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", viewer );
    const auto pref = static_cast<ExternalRidAuotCallPref>( pref_value );
    
    bool valid = false;
    switch( pref )
    {
      case ExternalRidAuotCallPref::DoNotCall:
      case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
      case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
      case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
      case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
        answer = pref;
        valid = true;
        break;
      
      //Not using default: to preserve compiler warning if we add a new value and dont get it here
    }//switch( pref )
    
    if( !valid )
      throw runtime_error( "Invalid AlwaysCallExternalRid value (" + std::to_string(pref_value) + ")" );
  }catch( std::exception &e )
  {
    Wt::log("error") << "RemoteRid::external_rid_call_pref: error " << e.what();
  }//try / catch
  
  return answer;
}//ExternalRidAuotCallPref external_rid_call_pref( InterSpec *viewer )


SimpleDialog *RemoteRid::startRemoteRidDialog( InterSpec *viewer,
                                      function<void(AuxWindow *, RemoteRid *)> callback )
{
  assert( viewer );
  if( !viewer )
    return nullptr;
  
  const bool showWarning = InterSpecUser::preferenceValue<bool>( "ExternalRidWarn", viewer );
  
  if( !showWarning )
  {
    auto res = RemoteRid::createDialog( viewer );
    if( callback )
      callback( res.first, res.second );
    
    return nullptr;
  }//if( !showWarning )
  
  wApp->useStyleSheet( "InterSpec_resources/RemoteRid.css" );
  
  SimpleDialog *dialog = new SimpleDialog( "Warning" );
  dialog->addStyleClass( "ExternalRidWarningDialog" );
  dialog->setWidth( 300 );
  
  string warntxt = "<p>Using this feature will send your data to"
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  " either another program on your computer that you will choose, or to"
#endif
  " a remote URL that you will enter.</p>"
  "<p>Are you sure you would like to continue?</p>"
  ;
  WText *txt = new WText( warntxt, dialog->contents() );
  txt->addStyleClass( "content" );
  txt->setInline( false );
  
  WCheckBox *cb = new WCheckBox( "Don't show again", dialog->contents() );
  cb->addStyleClass( "NoShowAgain" );
  cb->setInline( false );
  
  Wt::WPushButton *btn = dialog->addButton( "Cancel" );
  btn->clicked().connect( std::bind([callback](){
    if( callback )
      callback( nullptr, nullptr );
  }) );
  
  btn = dialog->addButton( "Continue" );
  btn->clicked().connect( std::bind([viewer,callback,cb](){
    if( cb->isChecked() )
      InterSpecUser::setPreferenceValue(viewer->m_user, "ExternalRidWarn", false, viewer );
      
    auto res = RemoteRid::createDialog( viewer );
    if( callback )
      callback( res.first, res.second );
  }) );
  
  return dialog;
}//void startRemoteRidDialog( viewer, callback )


pair<AuxWindow *, RemoteRid *> RemoteRid::createDialog( InterSpec *viewer )
{
  AuxWindow *window = new AuxWindow( WString::fromUTF8("External RID"),
                                    (Wt::WFlags<AuxWindowProperties> (AuxWindowProperties::DisableCollapse)
                                     | AuxWindowProperties::SetCloseable
                                     | AuxWindowProperties::TabletNotFullScreen)
                                    );
  
  WPushButton *close = window->addCloseButtonToFooter();
  close->clicked().connect( window, &AuxWindow::hide );
  window->addHelpInFooter( window->footer(), "external-rid" );
  
  RemoteRid *w = new RemoteRid( viewer, window->contents() );
  
  window->resizeToFitOnScreen();
  window->centerWindowHeavyHanded();
  
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( window->footer() );
  qr_btn->setText( "QR Code" );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [w](){
    try
    {
      const string url = "interspec://remoterid/?" + Wt::Utils::urlEncode(w->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, "Remote-RID State", "Current Remote-RID state." );
    }catch( std::exception &e )
    {
      passMessage( "Error creating QR code: " + std::string(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES
  
  return { window, w };
}//RemoteRid::createDialog(...)


RemoteRid::RemoteRid( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer ),
  m_rest_rid( nullptr )
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  , m_menu( nullptr )
  , m_exe_rid( nullptr )
#endif
{
  assert( viewer );
  if( !viewer )
    return;
 
  wApp->useStyleSheet( "InterSpec_resources/RemoteRid.css" );
  
  addStyleClass( "RemoteRid" );
  
  m_rest_rid = new RestRidImp::ExternalRidWidget( RestRidImp::ExternalRidWidget::ServiceType::Rest, this, viewer );
  
  WGridLayout *layout = new WGridLayout( this );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );
  
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  m_exe_rid = new RestRidImp::ExternalRidWidget( RestRidImp::ExternalRidWidget::ServiceType::Exe, this, viewer );
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  stack->setOverflow( WContainerWidget::OverflowAuto );
  
  stack->addWidget( m_rest_rid );
  stack->addWidget( m_exe_rid );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_menu = new WMenu( stack, Wt::Vertical );
  //const bool is_phone = m_interspec->isPhone();
  //const char *style_class = is_phone ? "VerticalNavMenuPhone HeavyNavMenuPhone"
  //                                   : "VerticalNavMenu HeavyNavMenu SideMenu";
  const char * const style_class = "VerticalNavMenu HeavyNavMenu SideMenu";
  m_menu->addStyleClass( style_class );
  
  
  WMenuItem *item = new WMenuItem( "Remote URL", m_rest_rid );
  m_menu->addItem( item );
  
  item->clicked().connect( std::bind([item,this](){
    m_menu->select( item );
    item->triggered().emit( item );
  }) );
  
  item = new WMenuItem( "Executable", m_exe_rid );
  m_menu->addItem( item );
  item->clicked().connect( std::bind([item,this](){
    m_menu->select( item );
    item->triggered().emit( item );
  }) );
  
  
  layout->setContentsMargins( 9, 0, 9, 0 );
  layout->addWidget( m_menu, 0, 0, AlignLeft );
  
  layout->addWidget( stack, 0, 1 );
  layout->setColumnStretch( 1, 1 );
  
  // If we have it setup to always call the EXE, then show this tab, otherwise show REST tab.
  const auto pref_value = RemoteRid::external_rid_call_pref( m_interspec );
  
  switch( pref_value )
  {
    case ExternalRidAuotCallPref::DoNotCall:
      break;
      
    case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
    case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
      m_menu->select( 0 );
      break;
      
    case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
    case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
      m_menu->select( 1 );
      break;
  }//switch( prev_value )
  
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo )
  {
    m_menu->itemSelected().connect( std::bind([](){
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( !undoRedo || !undoRedo->canAddUndoRedoNow() )
        return;
      auto undo_redo = [](){
        InterSpec *viewer = InterSpec::instance();
        RemoteRid *rid = viewer ? viewer->remoteRid() : nullptr;
        if( rid && rid->m_menu )
          rid->m_menu->select( rid->m_menu->currentIndex() ? 0 : 1 );
      };
      undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Change RID location menu" );
    }) );
  }//if( undoRedo )
  
#else
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setColumnStretch( 0, 1 );
  layout->addWidget( m_rest_rid, 0, 0 );
#endif
}//RemoteRid constructor



std::shared_ptr<SpecUtils::SpecFile> RemoteRid::fileForAnalysis( InterSpec *interspec,
                                                              const WFlags<AnaFileOptions> flags )
{
  try
  {
    if( !interspec )
      return nullptr;
    
    shared_ptr<const SpecMeas> foreground_file = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    if( !foreground_file )
      return nullptr;
    
    auto answer = make_shared<SpecUtils::SpecFile>( *foreground_file );
    
    // Get rid of some information the external service doesnt need
    answer->set_instrument_id( "" );
    answer->set_uuid( "" );
    answer->set_measurement_location_name( "" );
    answer->set_parse_warnings( {} );
    answer->set_filename( "" );
    answer->set_remarks( {} );
    answer->set_detectors_analysis( {} );
    answer->clear_multimedia_data();
    
    for( const auto &m : answer->measurements() )
    {
      answer->set_title( "", m );
      answer->set_remarks( {}, m );
      if( m->has_gps_info() )
        answer->set_position( -999.9, -999.9, SpecUtils::time_point_t{}, m );
    }
    
    
    if( foreground_file->passthrough() && !flags.testFlag(OnlyDisplayedSearchSamples) )
      return answer;
    
    answer->remove_measurements( answer->measurements() );
    
    auto disp_fore = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !disp_fore )
      return nullptr;
    
    auto fore = make_shared<SpecUtils::Measurement>( *disp_fore );
    fore->set_source_type(SpecUtils::SourceType::Foreground);
    
    auto disp_back = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
    
    answer->add_measurement( fore, !disp_back );
    if( disp_back )
    {
      auto back = make_shared<SpecUtils::Measurement>( *disp_back );
      back->set_source_type( SpecUtils::SourceType::Background );
      back->set_title( "" );
      back->set_remarks( {} );
      if( back->has_gps_info() )
        back->set_position( -999.9, -999.9, SpecUtils::time_point_t{} );
      
      answer->add_measurement( back, true );
    }//if( disp_back )
    
    return answer;
  }catch( std::exception &e )
  {
    cerr << "Unexpected exception in RemoteRid::fileForAnalysis(): " << e.what() << endl;
    assert( 0 );
  }//try / catch
  
  return nullptr;
}//std::shared_ptr<SpecUtils::SpecFile> fileForAnalysis()


void RemoteRid::alwaysCallRestAnaChecked()
{
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  m_exe_rid->uncheckAlwaysCallAnalysis();
#endif
}//void alwaysCallRestAnaChecked()

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
void RemoteRid::alwaysCallExeAnaChecked()
{
  m_rest_rid->uncheckAlwaysCallAnalysis();
}//void alwaysCallExeAnaChecked()
#endif


RestRidImp::ExternalRidWidget *RemoteRid::restRidTool()
{
  return m_rest_rid;
}

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
RestRidImp::ExternalRidWidget *RemoteRid::exeRidTool()
{
  return m_exe_rid;
}
#endif


void RemoteRid::startAutomatedOnLoadAnalysis( InterSpec *interspec,
                                              const Wt::WFlags<AnaFileOptions> flags )
{
  WApplication *app = WApplication::instance();
  if( !interspec || !app )
    return;
  
  const auto pref_value = RemoteRid::external_rid_call_pref( interspec );
  auto service_type = RestRidImp::ExternalRidWidget::ServiceType::Exe;
  switch( pref_value )
  {
    case ExternalRidAuotCallPref::DoNotCall:
      return;
      
    case ExternalRidAuotCallPref::AlwaysUseRestWithToast:
    case ExternalRidAuotCallPref::AlwaysUseRestWithDialog:
      service_type = RestRidImp::ExternalRidWidget::ServiceType::Rest;
      break;
    
    case ExternalRidAuotCallPref::AlwaysUseExeWithToast:
    case ExternalRidAuotCallPref::AlwaysUseExeWithDialog:
      service_type = RestRidImp::ExternalRidWidget::ServiceType::Exe;
      break;
  }//switch( pref_value )
  
  shared_ptr<SpecUtils::SpecFile> meas = fileForAnalysis( interspec, flags );
  if( !meas )
    return;
  
  const string uri = RestRidImp::ExternalRidWidget::mostRecentUserUrl( service_type, interspec );
  
  switch( service_type )
  {
    case RestRidImp::ExternalRidWidget::ServiceType::Rest:
    {
      RestRidImp::RestRidInterface *rest = new RestRidImp::RestRidInterface( interspec, wApp->domRoot() );
      rest->m_resource->setAnaFlags( flags );
      rest->m_resource->setDrf( "auto" );
      rest->m_analysis_results_signal.connect( rest, &RestRidImp::RestRidInterface::handleAutoAnaResultSuccess );
      rest->m_analysis_error_signal.connect( rest, &RestRidImp::RestRidInterface::handleAutoAnaResultFailure );
      rest->startRestAnalysis( uri );
    
      break;
    }//case RestRidImp::ExternalRidWidget::ServiceType::Rest:
      
    case RestRidImp::ExternalRidWidget::ServiceType::Exe:
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
      RestRidImp::ExternalRidWidget::startExeAnalysis( uri, nullptr, "auto", meas );
#else
      assert( 0 );
      throw runtime_error( "Cannot execute EXE analysis." );
#endif
      break;
  }//switch( service_type )
}//void startAutomatedOnLoadAnalysis( InterSpec *interspec )


void RemoteRid::handleOpeningRemoteRidTool( InterSpec *interspec )
{
  if( !interspec )
    return;
  
  interspec->createRemoteRidWindow();
  WApplication *app = WApplication::instance();
  if( app )
    app->triggerUpdate();
}//void handleOpeningRemoteRidTool( InterSpec *interspec )


void RemoteRid::disableAutoRemoteRid( InterSpec *interspec )
{
  if( !interspec )
    return;
  
  try
  {
    InterSpecUser::setPreferenceValue(interspec->m_user, "AlwaysCallExternalRid", 0, interspec);
  }catch(...)
  {
    assert( 0 );
  }
}//void disableAutoRemoteRid( InterSpec *interspec )


void RemoteRid::handleShowingRemoteRidRefLines( InterSpec *interspec, std::string signal )
{
  if( !interspec )
    return;
  
  SpecUtils::ireplace_all( signal, "showRemoteRidRefLines-", "" );
  vector<string> nuclides;
  SpecUtils::split( nuclides, signal, "," );
  
  ReferencePhotopeakDisplay *display = interspec->referenceLinesWidget();
  
  bool hadToCreateWindow = false;
  if( !display && !interspec->toolTabsVisible() )
  {
    hadToCreateWindow = true;
    interspec->showGammaLinesWindow();
    display = interspec->referenceLinesWidget();
  }
  
  if( !display )
  {
    cerr << "Failed to get referenceLinesWidget!" << endl;
    assert( 0 );
    return;
  }
  
  display->clearAllLines();
  display->setShieldingMaterialAndThickness( "", "0 cm" );
  
  size_t ndisp = 0;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  for( string &nuc : nuclides )
  {
    if( ndisp )
      display->persistCurentLines();
  
    if( !db->nuclide(nuc) )
      continue;
    
    display->setNuclideAndAge( nuc, "" );
    ++ndisp;
  }//for( string &nuc : nuclides )
  
  if( hadToCreateWindow && interspec->isPhone() )
    interspec->closeGammaLinesWindow();
  
  WApplication *app = WApplication::instance();
  if( app )
    app->triggerUpdate();
}//void handleShowingRemoteRidRefLines( InterSpec *interspec, std::string signal)


void RemoteRid::handleAppUrl( std::string query_str )
{
  {// Begin remove leading protocal/path from query_str
    const string path = "remoterid?";
    const char * const str_start = query_str.c_str();
    const char * const str_end = str_start + query_str.size();
    const auto pos = std::search( str_start, str_end, begin(path), end(path),
                                 [](unsigned char a, unsigned char b) -> bool {
      return (std::tolower(a) == std::tolower(b));
    } );
    
    if( pos != str_end )
      query_str = query_str.substr( (pos - str_start) + path.size() );
  }// End remove leading protocal/path from query_str
  
  map<string,string> parts = AppUtils::query_str_key_values( query_str );
  
  if( !parts.count("VER") || (parts["VER"] != "1") )
    throw runtime_error( "RemoteRid::handleAppUrl: missing or invalid 'VER'" );
  
  // This is kinda arbitrary, but we wont support setting both a EXE and URL; doing both might be confusing to the
  //  user, and we dont want to take any chances, since this involves potentually sending data off the
  //  users device.
  if( parts.count("PATH") && parts.count("URL") )
    throw runtime_error( "interspec://remoterid?... URI can not contain both PATH and URL parameters." );
  
  if( !parts.count("PATH") && !parts.count("URL") )
    throw runtime_error( "interspec://remoterid?... URI must contain either PATH or URL parameters." );
  
  const bool always_call = parts.count("ALWAYS")
              ? ((parts["ALWAYS"] == "1") || SpecUtils::iequals_ascii(parts["ALWAYS"], "true"))
              : false;
  const bool show_dialog = (always_call && parts.count("DIALOG"))
              ? ((parts["DIALOG"] == "1") || SpecUtils::iequals_ascii(parts["DIALOG"], "true"))
              : false;
  
  string exe_path = parts.count("PATH") ? Wt::Utils::urlDecode(parts["PATH"]) : string();
  const string url_path = parts.count("URL") ? Wt::Utils::urlDecode(parts["URL"]) : string();
  
  // No matter what value is given for "&none=...", or what other parameters are specified, we will try to reset things
  if( parts.count("NONE") || (url_path.empty() && exe_path.empty()) )
  {
    SimpleDialog *dialog = new SimpleDialog( "Stop using External-RID?",
                                            "This will reset you External-RID preferences.<br />"
                                            "Any URL"
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
                                            " or EXE's"
#endif
                                            " you have entered will be removed, and the External-RID"
                                            " service will no longer be used." );
    WPushButton *cont = dialog->addButton( "Continue" );
    cont->clicked().connect( std::bind([](){
      InterSpec *interspec = InterSpec::instance();
      assert( interspec );
      if( !interspec )
        return;
      InterSpecUser::setPreferenceValue( interspec->m_user, "AlwaysCallExternalRid", static_cast<int>(0), interspec );
      InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidWarn", true, interspec );
      InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidUrl", string(), interspec );
      InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidExe", string(), interspec );
      passMessage( "External-RID preferences have been reset.", WarningWidget::WarningMsgInfo );
    }) );
    dialog->addButton( "Cancel" );
    
    return;
  }//if( parts.count("NONE") )
  
#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT )
  if( !exe_path.empty() )
    throw runtime_error( "This build of InterSpec does not support PATH in a interspec://remoterid?... URI." );
#else
  if( !exe_path.empty() )
  {
    if( SpecUtils::icontains(exe_path,"http")
       || (exe_path.find(",;#<$+%>!`&*'\"|{?\"=}@") != string::npos) )
      throw runtime_error( "Filename contained invalid characters." );
    
    if( !AppUtils::locate_file(exe_path, false, 0, true) )
      throw runtime_error( "Could not find the '" + exe_path + "' executable asked for on your system." );
  }//if( !exe_path.empty() )
#endif
  
  const string title = "Use External-RID Service?";
  const string desc = "Are you sure you would like to use the '"
  + (exe_path.empty() ? url_path : exe_path) + "' service to call"
  + (always_call ? " whenever spectrum files are loaded" : " when the External-RID tool is used")
  + "?<br />"
  + (show_dialog ? "A modal-dialog with results will be shown." : "")
  + "Your spectrum data will be sent to this "
  + (exe_path.empty() ? "service." : "executable.")
  + "<br /><em>InterSpec</em> does not check the validity, trustworthiness, or anything else about this "
  + (exe_path.empty() ? "service" : "executable") + " - it just blindly sends it your spectroscopy data."
  "<br />If you are unsure, select <b>No</b>.";
  
  SimpleDialog *dialog = new SimpleDialog( title, desc );
  WPushButton *btn = dialog->addButton( "Yes" );
  
  const auto set_to_new_prefs = [exe_path,url_path,always_call,show_dialog](){
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    
    ExternalRidAuotCallPref pref_value = ExternalRidAuotCallPref::DoNotCall;
    if( always_call )
    {
      if( !url_path.empty() )
      {
        pref_value = show_dialog ? ExternalRidAuotCallPref::AlwaysUseRestWithDialog
                                 : ExternalRidAuotCallPref::AlwaysUseRestWithToast;
      }else
      {
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
        pref_value = show_dialog ? ExternalRidAuotCallPref::AlwaysUseExeWithDialog
                                 : ExternalRidAuotCallPref::AlwaysUseExeWithToast;
#else
        pref_value = ExternalRidAuotCallPref::DoNotCall;
        passMessage( "External-RID can-not call executable on this platform"
                     " - disabling auto-calling of external RID.", WarningWidget::WarningMsgHigh );
#endif
      }//if( !url_path.empty() ) / else
    }//if( always_call )
    
    InterSpecUser::setPreferenceValue( interspec->m_user, "AlwaysCallExternalRid", static_cast<int>(pref_value), interspec );
    //InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidWarn", false, interspec );
    if( !url_path.empty() )
      InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidUrl", url_path, interspec );
    if( !exe_path.empty() )
      InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidExe", exe_path, interspec );
    passMessage( "External-RID preferences have been updated.", WarningWidget::WarningMsgInfo );
  };//set_to_new_prefs
  
  
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  int prev_always_call = 0;
  string prev_url_path, prev_exe_path;
  
  try
  {
    prev_always_call = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", interspec );
  }catch( std::exception & ){ cerr << "Exception getting AlwaysCallExternalRid." << endl; }
  
  try
  {
    prev_url_path = InterSpecUser::preferenceValue<string>( "ExternalRidUrl", interspec );
  }catch( std::exception & ){ cerr << "Exception getting ExternalRidUrl" << endl; }
  
  try
  {
    prev_exe_path = InterSpecUser::preferenceValue<string>( "ExternalRidExe", interspec );
  }catch( std::exception & ){ cerr << "Exception getting ExternalRidExe." << endl; }
  
  const auto set_to_old_prefs = [prev_always_call, prev_url_path, prev_exe_path](){
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    
    InterSpecUser::setPreferenceValue( interspec->m_user, "AlwaysCallExternalRid", prev_always_call, interspec );
    //InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidWarn", false, interspec );
    InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidUrl", prev_url_path, interspec );
    InterSpecUser::setPreferenceValue( interspec->m_user, "ExternalRidExe", prev_exe_path, interspec );
    passMessage( "External-RID preferences have been reverted back to original.", WarningWidget::WarningMsgInfo );
  };//set_to_old_prefs
  
  
  btn->clicked().connect( std::bind([set_to_old_prefs,set_to_new_prefs](){
    set_to_new_prefs();
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( !undoRedo || undoRedo->isInUndoOrRedo() )
      return;
    undoRedo->addUndoRedoStep( set_to_old_prefs, set_to_new_prefs, "Set Remote-RID from URI." );
  }) );
  
  dialog->addButton( "No" );
}//void handleAppUrl( std::string query_str )


std::string RemoteRid::encodeStateToUrl() const
{
  string url = "ver=1";
  
  RestRidImp::ExternalRidWidget *w = nullptr;
#if( ANDROID || IOS || BUILD_FOR_WEB_DEPLOYMENT )
  w = m_rest_rid;
#else
  w = (m_menu->currentIndex() == 0) ? m_rest_rid : m_exe_rid;
#endif
  assert( w );
  if( !w || !w->isValidUrl() )
    return url + "&none=1";
  
  url += ((w->serviceType() == RestRidImp::ExternalRidWidget::ServiceType::Rest) ? "&url=" : "&path=");
  url += Wt::Utils::urlEncode(w->url());
  
  if( w->alwaysAutoCallAnalysis() )
    url += "&always=1";
  
  if( w->showDialogForAutoCalledResults() )
    url += "&dialog=1";
  
  return url;
}//std::string encodeStateToUrl() const
