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
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )

#include <boost/process.hpp>
#include <boost/filesystem.hpp>

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

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/RemoteRid.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"


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


/** Looks at the file path passed in to try and find the file wether it is a relative path, or
 an absolute path.
 
 Returns if the path could be found, and modifies the filename passed in such that you can open that
 path as an ifstream or something.
 */
bool locate_file( string &filename, const bool is_dir )
{
  auto check_exists = [is_dir]( const string &name ) -> bool {
    return is_dir ? SpecUtils::is_directory(name) : SpecUtils::is_file(name);
  };//auto check_exists
  
  if( SpecUtils::is_absolute_path(filename) )
    return check_exists(filename);
  
  // Check if path is there, relative to CWD
  if( check_exists(filename) )
    return true;
  
  // We'll look relative to the executables path, but note that if we started from a symlink, I
  //  think it will resolve relative to actual executable
  try
  {
#ifdef __APPLE__
    char path_buffer[PATH_MAX + 1] = { '\0' };
    uint32_t size = PATH_MAX + 1;
    
    if (_NSGetExecutablePath(path_buffer, &size) != 0) {
      return false;
    }
    
    path_buffer[PATH_MAX] = '\0'; // JIC
    const string exe_path = path_buffer;
#elif( defined(_WIN32) )
    //static_assert( 0, "Need to test this EXE path stuff on Windows..." );
    wchar_t wbuffer[2*MAX_PATH];
    const DWORD len = GetModuleFileNameW( NULL, wbuffer, 2*MAX_PATH );
    if( len <= 0 )
      throw runtime_error( "Call to GetModuleFileName falied" );
    
    const string exe_path = SpecUtils::convert_from_utf16_to_utf8( wbuffer );
#else // if __APPLE__ / Win32 / else
    
    char path_buffer[PATH_MAX + 1] = { '\0' };
    const ssize_t ret = readlink("/proc/self/exe", path_buffer, PATH_MAX);
    
    if( (ret == -1) || (ret > PATH_MAX) )
      throw runtime_error( "Failed to read line" );
    
    assert( ret < PATH_MAX );
    path_buffer[ret] = '\0';
    path_buffer[PATH_MAX] = '\0'; // JIC
    const string exe_path = path_buffer;
#endif // else not __APPLE__
    
    string trial_parent_path = exe_path;
    if( !SpecUtils::make_canonical_path(trial_parent_path) )
      throw runtime_error( "Failed to make trial_parent_path canonical" );
    
    trial_parent_path = SpecUtils::parent_path(trial_parent_path);
    
    string trialpath = SpecUtils::append_path( trial_parent_path, filename );
    
    if( check_exists(trialpath) )
    {
      filename = trialpath;
      return true;
    }
    
    if( boost::filesystem::is_symlink(trialpath) )
    {
      trialpath = boost::filesystem::read_symlink(trialpath).string<string>();
      if( check_exists(trialpath) )
      {
        filename = trialpath;
        return true;
      }
    }//if( is symlink )
    
    // We could try again, going up one more directory, incase executable is in "Debug" directory
    //  or something, but we'll skip this for the moment
    //trial_parent_path = SpecUtils::parent_path(trial_parent_path);
    //trialpath = SpecUtils::append_path( trial_parent_path, filename );
    //if( check_exists(trialpath) )
    //{
    //  filename = trialpath;
    //  return true;
    //}
  }catch( std::exception &e )
  {
    //cerr << "Caught exception: " << e.what() << endl;
  }
  
  
  // Check relative to environments "PATH"
  const char *env_path = std::getenv("PATH");
  if( !env_path )
    return false;
  
  // grab the PATH system variable, and check if the input path is relative to any of those
  vector<string> paths;
#if( defined(_WIN32) )
  const char *delims = ";";
#else
  const char *delims = ":";
#endif
  SpecUtils::split( paths, env_path, ";" );
  
  for( string base : paths )
  {
    const string trialpath = SpecUtils::append_path( base, filename );
    
    if( check_exists(filename) )
    {
      filename = trialpath;
      return true;
    }
  }//for( string base : paths )
  
  // Out of luck, we failed if we're here.

  return false;
}//bool locate_file( ... )
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
  void handleAnaResultSuccessForToast( string response );
  void handleAnaResultFailureForToast( string response );
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
      
      if( locate_file( exe_name, false ) )
        return exe_name;
      
      return "";
    };//check_for_default() lamda
    
    
    try
    {
      const string prev_urls = InterSpecUser::preferenceValue<string>( prefname, interspec );
      vector<string> prev;
      SpecUtils::split( prev, prev_urls, ";" );
      if( prev.empty() )
        return check_for_default();
      
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
        const char *url_regex = "^(http(s)?:\\/\\/)[\\w.-]+(?:\\.[\\w\\.-]+)+[\\w\\-\\._~:/?#[\\]@!\\$&'\\(\\)\\*\\+,;=.]+$";
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
    
    m_alwaysDoAnalysisCb = new WCheckBox( "Always submit analysis on spectrum load" );
    m_alwaysDoAnalysisCb->addStyleClass( "AlwaysSubmitAna" );
    layout->addWidget( m_alwaysDoAnalysisCb, 5, 0, 1, 2 );
    
    const int always_call_index = (m_type == ServiceType::Rest) ? 1 : 2;
    const int always_call = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", m_interspec );
    if( always_call == always_call_index )
      m_alwaysDoAnalysisCb->setChecked( true );
    
    m_alwaysDoAnalysisCb->checked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    m_alwaysDoAnalysisCb->unChecked().connect( this, &ExternalRidWidget::handleAlwaysCallCheckChanged );
    
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
  }
  
  void handleAlwaysCallCheckChanged()
  {
    const int always_call_index = (m_type == ServiceType::Rest) ? 1 : 2;
    
    if( m_alwaysDoAnalysisCb->isChecked() )
    {
      try
      {
        InterSpecUser::setPreferenceValue(m_interspec->m_user, "AlwaysCallExternalRid", always_call_index, m_interspec);
        switch( m_type )
        {
          case ServiceType::Rest:
            m_remote_rid->alwaysCallRestAnaChecked();
            break;
          case ServiceType::Exe:
            m_remote_rid->alwaysCallExeAnaChecked();
            break;
        }//switch( m_type )
      }catch( std::exception & )
      {
        assert(0);
      }
    }else
    {
      try
      {
        InterSpecUser::setPreferenceValue(m_interspec->m_user, "AlwaysCallExternalRid", 0, m_interspec);
      }catch( std::exception & )
      {
        assert(0);
      }
    }
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo_redo = [always_call_index](){
        InterSpec *viewer = InterSpec::instance();
        RemoteRid *rid = viewer ? viewer->remoteRid() : nullptr;
        if( !rid )
          return;
      
        RestRidImp::ExternalRidWidget *tool = rid->restRidTool();
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
        if( always_call_index == 2 )
          tool = rid->exeRidTool();
#endif
        WCheckBox *cb = tool ? tool->m_alwaysDoAnalysisCb : nullptr;
        
        if( cb )
          cb->setChecked( !cb->isChecked() );
        
        if( tool )
          tool->handleAlwaysCallCheckChanged();
      };//undo_redo
      
      undoRedo->addUndoRedoStep( undo_redo, undo_redo, "Toggle always do External Rid on spectrum load" );
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
  
public:
  static void makeToastNotificationForAnaResult( std::shared_ptr<int> rcode, std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )
  {
    assert( rcode && result && m );
    WApplication *app = WApplication::instance();
    assert( app );
    if( !app )
      return;
    
    std::lock_guard<mutex> lock( *m );
    
    string message;
    const int code = *rcode;
    const string &res = *result;
    if( code == 0 )
    {
      vector<string> nuclides;
      vector<pair<string,string>> iso_descrips;
      
      try
      {
        const FullSpecResults result = json_to_results( res );
        
        string iso_str = result.isotopeString;
        SpecUtils::ireplace_all( iso_str, "+", ", " );
        
        if( iso_str.empty() )
        {
          if( result.errorMessage.size() )
            message = "GADRAS RID Error: " + result.errorMessage;
          else if( result.code )
            message = "GADRAS RID Error code " + to_string(result.code);
          else
            message = "GADRAS RID: nothing identified";
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
        InterSpec *interspec = InterSpec::instance();
        ReferencePhotopeakDisplay *reflines = interspec ? interspec->referenceLinesWidget() : nullptr;
        if( reflines )
          reflines->setExternalRidResults( "GADRAS", iso_descrips );
      }catch( std::exception &e )
      {
        message = "Failed to parse GADRAS RID results.";
        
        cerr << "makeToastNotificationForAnaResult: Failed to parse results from GADRAS: "
        << e.what() << endl;
      }//try / catch
      
      WStringStream js;
      js << message;
      
      js << "<div class=\"RemoteRidToastButtons\">";
      if( !nuclides.empty() )
      {
        string nucstr;
        for( size_t i = 0; i < nuclides.size(); ++i )
          nucstr += (i ? "," : "") + nuclides[i];
        
        js <<
        "<div onclick=\"Wt.emit('" << app->root()->id() << "',{name:'miscSignal'}, 'showRemoteRidRefLines-" << nucstr << "');"
        "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
        "return false;\">Show Ref. Lines</div>";
      }//if( !nuclides.empty() )
      
      js <<
      "<div onclick=\"Wt.emit('" << app->root()->id() << "',{name:'miscSignal'}, 'openRemoteRidTool');"
      "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
      "return false;\">Open Remote RID</div>";
      
      js << "</div>";
      
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        viewer->warningWidget()->addMessageUnsafe( js.str(), WarningWidget::WarningMsgExternalRiid, 20000 );
    }else // if( code != 0 )
    {
      // TODO: check for return code, etc., and customize message based on that
      
      WStringStream js;
      js << res;
      js <<
      "<div class=\"RemoteRidToastButtons\">"
      "<div onclick=\"Wt.emit('" << app->root()->id() << "',{name:'miscSignal'}, 'openRemoteRidTool');"
      "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
      "return false;\">Open Remote RID</div>"
      "<div onclick=\"Wt.emit('" << app->root()->id() << "',{name:'miscSignal'}, 'disableRemoteRid');"
      "try{$(this.parentElement.parentElement.parentElement).remove();}catch(e){}"
      "return false;\">Disable Remote RID</div>"
      "</div>";
      
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        viewer->warningWidget()->addMessageUnsafe( js.str(), WarningWidget::WarningMsgExternalRiid, 20000 );
    }//if( code == 0 ) / else
    
    wApp->triggerUpdate();
  }//void makeToastNotificationForAnaResult( std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )

public:
  static void startExeAnalysis( string exe_path, ExternalRidWidget *parent, string drf, shared_ptr<SpecUtils::SpecFile> spec_file )
  {
    vector<string> arguments;
    arguments.push_back( "--mode=command-line" );
    arguments.push_back( "--out-format=json" );
    arguments.push_back( "--drf" );
    arguments.push_back( drf );
    
    if( spec_file->num_measurements() < 2 )
      arguments.push_back( "--synthesize-background=1" );
    
    if( !locate_file(exe_path, false) )
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
    }//if( !locate_file(exe_path, false) )
    
    
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
      doUpdateFcn = boost::bind( &ExternalRidWidget::makeToastNotificationForAnaResult, rcode, result, m );
    
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
  
  
  void submitForAnalysis()
  {
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_status->setText( "Requesting analysis" );
    
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
  }//void submitForAnalysis()
  
  
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
    
  
  void requestServiceInfo()
  {
    m_retrieve_drfs_btn->disable();
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "Requesting service information." );
    m_current_drf_index = -1;
    m_drf_select->clear();
    
    switch( m_type )
    {
      case ServiceType::Rest:
        requestRestServiceInfo();
        break;
        
      case ServiceType::Exe:
        requestExeInfo();
        break;
    }//switch( m_type )
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
  
  
  void handleInfoResponseError( const std::string &msg )
  {
    m_current_drf_index = -1;
    m_drf_select->clear();
    m_drf_stack->setCurrentIndex( 0 );
    m_retrieve_drfs_btn->disable();
    m_result->setText( "" );
    m_status->setText( "" );
    m_error->setText( "Error retrieving service information: <code>" + msg + "</code>" );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
  }//void handleInfoResponseError()

  
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
      
      if( (results.chi2 > 0.0001f) && (results.isotopes.size() > 0) )
      {
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.2f", results.chi2 );
        rslttxt << "<div class=\"ResultChi2\">&chi;<sup>2</sup>=" << std::string(buffer)
        << ", DRF:&nbsp;" << Wt::Utils::htmlEncode(results.drf) << "</div>";
      }else
      {
        rslttxt << "<div class=\"ResultChi2\">DRF:&nbsp;" << Wt::Utils::htmlEncode(results.drf)
        << "</div>";
      }//if( output.chi_sqr > 0.0f )
      
      if( !results.analysisWarnings.empty() )
      {
        rslttxt << "<div class=\"AnaWarnings\">\n";
        for( const string &warning : results.analysisWarnings )
          rslttxt << "<div>" << Wt::Utils::htmlEncode(warning) << "</div>\n";
        rslttxt << "</div>\n";
      }//if( result.contains("analysisWarnings") )
      
      m_result->setText( rslttxt.str() );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_result) );
    }catch( std::exception &e )
    {
      m_error->setText( e.what() );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
    }
  }//void handleResultResponse()
  
  
  void handleResultResponseError( const std::string &msg )
  {
    m_submit->enable();
    m_result->setText( "" );
    m_status->setText( "" );
    
    m_error->setText( msg );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
  }//void handleResultResponseError()
  
  
  bool isValidUrl()
  {
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
          if( !locate_file(url, false) )
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
        //blah blah blah
        break;
      }//case ServiceType::Exe:
    }//switch( m_type )

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
  "    .then(response => response.json() )\n"
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
    // cout << "When deleting RestRidInterface, got w = " << w << ", ww = " << ww << ", instance = " << instance << endl;
    assert( w );
    
    delete instance;
  } );
}//deleteSelfForToast()

void RestRidInterface::handleAnaResultSuccessForToast( string response )
{
  auto result = std::make_shared<string>( response );
  auto rcode = std::make_shared<int>( 0 );
  auto m = make_shared<mutex>();
  
  RestRidImp::ExternalRidWidget::makeToastNotificationForAnaResult( rcode, result, m );
  
  deleteSelfForToast();
}//handleAnaResultSuccessForToast(...)

void RestRidInterface::handleAnaResultFailureForToast( string response )
{
  auto result = std::make_shared<string>( response );
  auto rcode = std::make_shared<int>( -1 );
  auto m = make_shared<mutex>();
  
  RestRidImp::ExternalRidWidget::makeToastNotificationForAnaResult( rcode, result, m );
  
  deleteSelfForToast();
}//handleAnaResultFailureForToast(...)



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
  AuxWindow *window = new AuxWindow( WString::fromUTF8("External RIID"),
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
  const bool is_phone = m_interspec->isPhone();
  const char *style_class = is_phone ? "VerticalNavMenuPhone HeavyNavMenuPhone SideMenuPhone"
                                     : "VerticalNavMenu HeavyNavMenu SideMenu";
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
  layout->setColumnStretch( 0, 1 );
  
  // If we have it setup to always call the EXE, then show this tab, otherwise show REST tab.
  const int always_call = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", viewer );
  m_menu->select( (always_call == 2) ? 1 : 0 );
  
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
      back->set_source_type(SpecUtils::SourceType::Background);
      answer->add_measurement( back, true );
    }
    
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
  
  const int always_call = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", interspec );
  if( (always_call != 1) && (always_call != 2) )
    return;
  
  shared_ptr<SpecUtils::SpecFile> meas = fileForAnalysis( interspec, flags );
  if( !meas )
    return;
  
  
  const auto service_type = (always_call == 1) ? RestRidImp::ExternalRidWidget::ServiceType::Rest
                                               : RestRidImp::ExternalRidWidget::ServiceType::Exe;
  
  const string uri = RestRidImp::ExternalRidWidget::mostRecentUserUrl( service_type, interspec );
  
  switch( service_type )
  {
    case RestRidImp::ExternalRidWidget::ServiceType::Rest:
    {
      RestRidImp::RestRidInterface *rest = new RestRidImp::RestRidInterface( interspec, wApp->domRoot() );
      rest->m_resource->setAnaFlags( flags );
      rest->m_resource->setDrf( "auto" );
      rest->m_analysis_results_signal.connect( rest, &RestRidImp::RestRidInterface::handleAnaResultSuccessForToast );
      rest->m_analysis_error_signal.connect( rest, &RestRidImp::RestRidInterface::handleAnaResultFailureForToast );
      rest->startRestAnalysis( uri );
    
      break;
    }//case RestRidImp::ExternalRidWidget::ServiceType::Rest:
      
    case RestRidImp::ExternalRidWidget::ServiceType::Exe:
      RestRidImp::ExternalRidWidget::startExeAnalysis( uri, nullptr, "auto", meas );
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
