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


using namespace std;
using namespace Wt;


namespace
{
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
std::string run_external_command( const string &exe, const vector<string> &args )
{
  namespace bp = boost::process;
  
  const auto pp = boost::filesystem::path(exe).parent_path();
  
  bp::ipstream proc_stdout, proc_stderr; //reading pipe-stream
  bp::child c( exe, bp::args(args), bp::start_dir(pp), bp::std_out > proc_stdout, bp::std_err > proc_stderr );
  
  c.wait();

  std::string output(std::istreambuf_iterator<char>(proc_stdout), {});
  std::string error(std::istreambuf_iterator<char>(proc_stderr), {});
  
  cout << "Process output: " << output << endl;
  cout << "Process stdout: " << error << endl;
  
  int result = c.exit_code();
  
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
  RemoteRid *m_remote_rid;
  ExternalRidWidget *m_widget;
  Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
  InterSpec *m_interspec;
  
public:
  RestRidInputResource( RemoteRid *remote_rid, ExternalRidWidget *widget, InterSpec *interspec, WObject* parent = 0 );
  
  virtual ~RestRidInputResource();
  virtual void handleRequest( const Wt::Http::Request &request,
                             Wt::Http::Response &response );
};//class RestRidInputResource : public Wt::WResource



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
  WLineEdit *m_url;
  WStackedWidget *m_drf_stack;
  WPushButton *m_retrieve_drfs_btn;
  WComboBox *m_drf_select;
  WPushButton *m_submit;
  WStackedWidget *m_status_stack;
  WText *m_error;
  WText *m_status;
  WText *m_result;
  RestRidInputResource *m_resource;
  WCheckBox *m_alwaysDoAnalysisCb;
  
  
  std::unique_ptr<JSignal<std::string>> m_info_response_signal;
  std::unique_ptr<JSignal<std::string>> m_info_error_signal;
  std::unique_ptr<JSignal<std::string>> m_analysis_results_signal;
  std::unique_ptr<JSignal<std::string>> m_analysis_error_signal;
  
  
  std::string mostRecentUserUrl()
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
      if( prev.empty() )
        return "";
      
      return prev.front();
    }catch( std::exception &e )
    {
      cerr << "Error in mostRecentUserUrl(): " << e.what() << endl;
      assert( 0 );
    }//try / catch
    
    return "";
  }//std::string mostRecentUserUrl()
  
  
  
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
      
      if( !prev.empty() && (prev.back() == url) )
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
  ExternalRidWidget( const ServiceType type, RemoteRid *remoterid, InterSpec *interspec, WContainerWidget *parent = 0 )
  : WContainerWidget( parent ),
   m_interspec( interspec ),
   m_remote_rid( remoterid ),
   m_type( type ),
   m_url( nullptr ),
   m_drf_stack( nullptr ),
   m_retrieve_drfs_btn( nullptr ),
   m_drf_select( nullptr ),
   m_submit( nullptr ),
   m_status_stack( nullptr ),
   m_error( nullptr ),
   m_status( nullptr ),
   m_result( nullptr ),
   m_resource( nullptr ),
   m_alwaysDoAnalysisCb( nullptr ),
   m_info_response_signal{},
   m_info_error_signal{},
   m_analysis_results_signal{},
   m_analysis_error_signal{}
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
    
    m_url->setText( WString::fromUTF8( mostRecentUserUrl() ) );
    
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
    
    m_submit = new WPushButton( "Submit Analysis" );
    layout->addWidget( m_submit, 3, 0, 1, 2, AlignCenter );
    m_submit->clicked().connect( this, &ExternalRidWidget::submitForAnalysis );
    m_submit->disable();
    
    if( m_type == ServiceType::Rest )
      m_resource = new RestRidInputResource( m_remote_rid, this, m_interspec, this );
    
    m_alwaysDoAnalysisCb = new WCheckBox( "Always submit analysis on spectrum load" );
    m_alwaysDoAnalysisCb->addStyleClass( "AlwaysSubmitAna" );
    layout->addWidget( m_alwaysDoAnalysisCb, 4, 0, 1, 2 );
    
    const int always_call_index = (m_type == ServiceType::Rest) ? 1 : 2;
    const int always_call = InterSpecUser::preferenceValue<int>( "AlwaysCallExternalRid", m_interspec );
    if( always_call == always_call_index )
      m_alwaysDoAnalysisCb->setChecked( true );
    
    m_alwaysDoAnalysisCb->checked().connect( std::bind([this,always_call_index](){
      try
      {
        InterSpecUser::setPreferenceValue(m_interspec->m_user, "AlwaysCallExternalRid", always_call_index, m_interspec);
        m_remote_rid->alwaysCallRestAnaChecked();
      }catch( std::exception & )
      {
        assert(0);
      }
    }) );
    
    m_alwaysDoAnalysisCb->unChecked().connect( std::bind([this](){
      try
      {
        InterSpecUser::setPreferenceValue(m_interspec->m_user, "AlwaysCallExternalRid", 0, m_interspec);
      }catch( std::exception & )
      {
        assert(0);
      }
    }) );
    
    
    m_interspec->displayedSpectrumChanged().connect( this, &ExternalRidWidget::displayedSpectrumChanged );
    
    if( m_type == ServiceType::Rest )
    {
      m_info_response_signal.reset( new JSignal<std::string>( this, "info_response", false ) );
      m_info_error_signal.reset( new JSignal<std::string>( this, "info_error", false ) );
      m_analysis_results_signal.reset( new JSignal<std::string>( this, "analysis_results", false ) );
      m_analysis_error_signal.reset( new JSignal<std::string>( this, "analysis_error", false ) );
      
      m_info_response_signal->connect( boost::bind( &ExternalRidWidget::handleInfoResponse, this, boost::placeholders::_1 ) );
      m_info_error_signal->connect( boost::bind( &ExternalRidWidget::handleInfoResponseError, this, boost::placeholders::_1 ) );
      m_analysis_results_signal->connect( boost::bind( &ExternalRidWidget::handleResultResponse, this, boost::placeholders::_1 ) );
      m_analysis_error_signal->connect( boost::bind( &ExternalRidWidget::handleResultResponseError, this, boost::placeholders::_1 ) );
    }//if( m_type == ServiceType::Rest )
    
    urlChanged();
  }//ExternalRidWidget
  
  void uncheckAlwaysCallAnalysis()
  {
    m_alwaysDoAnalysisCb->setChecked( false );
  }
  
protected:
  
  void startRestAnalysis()
  {
    const string ana_service_url = m_url->text().toUTF8();
    const string &resource_url = m_resource->url();
    
    // Note: the JS is defined in the C++ to avoid having a separate JS file that will end up being
    // packaged in builds that wont use it (i.e., most InterSpec builds); I dont like having various
    // `fetch` commands around that can call off of localhost; if nothing else this is a bad look.
    
    WStringStream js;
    js <<
    "let sendAnaRequest = function(formData){\n"
    "  fetch('" << ana_service_url << "/analysis', {method: 'POST', body: formData})\n"
    "    .then(response => response.json() )\n"
    "    .then(json_results => {\n"
    "      console.log( 'Got JSON analysis response:', json_results );\n"
    "      const result_str = JSON.stringify(json_results);\n"
    << "      " << m_analysis_results_signal->createCall("result_str") << "\n"
    "    })\n"
    "  .catch( error => {\n"
    "    console.error( 'Error getting analysis input:', error);\n"
    "    const error_str = 'Error retrieving results from server: ' + error.toString();\n"
    "  " << m_analysis_error_signal->createCall( "error_str" ) << ";\n"
    "  });\n"
    "};\n";
    
    js << "fetch('" << resource_url << "', {method: 'POST'})\n"
    << ".then( response => response.formData() )\n"
    << ".then( data => {\n"
    //<< "  console.log( 'Got FormData wt app!' );\n"
    << "  sendAnaRequest(data);"
    //<< "  const data_str = JSON.stringify(data);\n"
    //<< "  " << m_info_response_signal.createCall( "data_str" ) << ";\n"
    << "})\n"
    << ".catch( error => {\n"
    << "  console.error( 'Error getting analysis input:', error);"
    << "  const error_str = 'Error getting analysis input:' + error.toString();\n"
    << "  " << m_analysis_error_signal->createCall( "error_str" ) << ";\n"
    << "});\n";
    
    doJavaScript( js.str() );
  }//void startRestAnalysis()
  
  
  void receiveExeAnalysis( std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )
  {
    assert( result && m );
    assert( WApplication::instance() );
    
    std::lock_guard<mutex> lock( *m );
    
    const string &res = *result;
    if( SpecUtils::starts_with(res, "Success:") )
    {
      const string json = SpecUtils::trim_copy(res.substr(8));
      handleResultResponse( json );
    }else if( SpecUtils::starts_with(res, "Fail:") )
    {
      const string result = SpecUtils::trim_copy(res.substr(5));
      handleResultResponseError( result );
    }else
    {
      assert( 0 );
    }
    
    wApp->triggerUpdate();
  }//void receiveExeAnalysis( std::shared_ptr<string> result )
  
  
  void startExeAnalysis()
  {
    shared_ptr<SpecUtils::SpecFile> spec_file = RemoteRid::fileForAnalysis(m_interspec);
    assert( spec_file );
    
    if( !spec_file )
    {
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
      m_error->setText( "No displayed spectrum" );
      
      return;
    }//if( !spec_file )
    
    vector<string> arguments;
    arguments.push_back( "--mode=command-line" );
    arguments.push_back( "--out-format=json" );
    arguments.push_back( "--drf" );
    arguments.push_back( drf() );
    
    if( spec_file->num_measurements() < 2 )
      arguments.push_back( "--synthesize-background=1" );
    
    const string tmpfilename = SpecUtils::temp_file_name( "interspec_ana_" + wApp->sessionId(),
                                                          SpecUtils::temp_dir() );
    
    if( SpecUtils::is_file(tmpfilename) || SpecUtils::is_directory(tmpfilename) )
    {
      // Shouldnt ever happen
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
      m_error->setText( "Error creating temp file name." );
      return;
    }
    
    {//begin writing temp file
      ofstream tmpfile( tmpfilename.c_str(), ios::out | ios::binary );
      if( !tmpfile.is_open() )
      {
        m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
        m_error->setText( "Error opening temp file." );
        return;
      }
      
      spec_file->write_2012_N42( tmpfile );
    }//end writing temp file
    
    arguments.push_back( tmpfilename );
    
    
    const string exe_path = m_url->text().toUTF8();
    
    const string appsession = WApplication::instance()->sessionId();
    auto result = std::make_shared<string>();
    auto m = make_shared<mutex>();
    auto doUpdateFcn = wApp->bind( boost::bind( &ExternalRidWidget::receiveExeAnalysis, this, result, m ) );
    
    auto commandRunner = [tmpfilename,exe_path,arguments,doUpdateFcn,result,m,appsession](){
      try
      {
        string results = run_external_command( exe_path, arguments );
        
        if( results.empty() )
          throw runtime_error( "No output from running executable." );
        
        std::lock_guard<mutex> lock( *m );
        *result = "Success:" + results;
      }catch( std::exception &e )
      {
        std::lock_guard<mutex> lock( *m );
        *result = "Fail:" + string(e.what());
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
        startExeAnalysis();
        break;
    }//switch( m_type )
  }//void submitForAnalysis()
  
  
  void requestRestServiceInfo()
  {
    const string url = m_url->text().toUTF8();
    
    WStringStream js;
    
    js << "fetch('" << url << "/info', {method: 'POST'})"
    << ".then( response => response.json() )"
    << ".then( data => {"
    //          << "console.log( 'Got data:', data );"
    << "const data_str = JSON.stringify(data);"
    << m_info_response_signal->createCall( "data_str" )
    << "})"
    << ".catch( error => {"
    << "console.error( 'Error retrieving info:', error);"
    << "const error_str = error.toString();"
    << m_info_error_signal->createCall( "error_str" )
    << "});";
    
    doJavaScript( js.str() );
  }//void requestRestServiceInfo()
  
  
  void receiveExeDrfInfo( std::shared_ptr<string> result, std::shared_ptr<std::mutex> m )
  {
    assert( result && m );
    assert( WApplication::instance() );
    
    std::lock_guard<mutex> lock( *m );
    
    const string &res = *result;
    if( SpecUtils::starts_with(res, "Success:") )
    {
      string json = "{\"Options\":[{\"name\":\"drf\","
                        "\"possibleValues\":" + SpecUtils::trim_copy(res.substr(8)) + "}]}";
      handleInfoResponse( json );
    }else if( SpecUtils::starts_with(res, "Fail:") )
    {
      handleInfoResponseError( res.substr(5) );
    }else
    {
      assert( 0 );
    }
    
    wApp->triggerUpdate();
  }//void receiveExeDrfInfo( std::shared_ptr<string> result )
  
  
  void requestExeInfo()
  {
    const string exe_path = m_url->text().toUTF8();
    
    const string appsession = WApplication::instance()->sessionId();
    auto result = std::make_shared<string>();
    auto m = make_shared<mutex>();
    auto doUpdateFcn = wApp->bind( boost::bind( &ExternalRidWidget::receiveExeDrfInfo, this, result, m ) );
    
    auto commandRunner = [exe_path,doUpdateFcn,result,m,appsession](){
      try
      {
        vector<string> args{"--command-line", "--out-format", "json", "--drfs"};
        string drfs = run_external_command( exe_path, args );
        
        if( drfs.empty() )
          throw runtime_error( "No output from running executable." );
        
        std::lock_guard<mutex> lock( *m );
        *result = "Success:" + drfs;
      }catch( std::exception &e )
      {
        std::lock_guard<mutex> lock( *m );
        *result = "Fail:" + string(e.what());
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
        rslttxt << "<div class=\"ResultChi2\">&chi;<sup>2</sup>=" << std::string(buffer) << "</div>";
      }//if( output.chi_sqr > 0.0f )
      
      if( !results.analysisWarnings.empty() )
      {
        rslttxt << "<div class=\"AnaWarnings\">\n";
        for( const string &warning : results.analysisWarnings )
          rslttxt << "<div>" << warning << "</div>\n";
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
        const string url = m_url->text().toUTF8();
        try
        {
          const boost::filesystem::file_status status = boost::filesystem::status(url);
          if( !boost::filesystem::exists(status) )
            return false;
          
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
    m_drf_select->clear();
    m_drf_stack->setCurrentIndex( 0 );
      
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
  }//urlChanged()
  
  
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
      if( !spec )
        m_status->setText( "No spectrum displayed to submit." );
      else if( url.empty() )
        m_status->setText( "Entered URL looks to be invalid." );
    }else
    {
      m_submit->setHidden( false );
      m_status->setText( "Click &quot;Submit Analysis&quot; to get RIID results." );
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




RestRidInputResource::RestRidInputResource( RemoteRid *remote_rid, ExternalRidWidget *widget, InterSpec *interspec, WObject* parent )
: WResource( parent ),
m_remote_rid( remote_rid ),
m_widget( widget ),
m_app( WApplication::instance() ),
m_interspec( interspec )
{
  assert( m_remote_rid && m_widget && m_app && m_interspec );
}

RestRidInputResource::~RestRidInputResource()
{
  beingDeleted();
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
  shared_ptr<SpecUtils::SpecFile> spec_file = RemoteRid::fileForAnalysis(m_interspec);
  assert( spec_file );
  
  const string drf = m_widget->drf();
  string options = "{\"drf\": \"" + drf + "\"";
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


void RemoteRid::startRemoteRidDialog( InterSpec *viewer,
                                      function<void(AuxWindow *, RemoteRid *)> callback )
{
  assert( viewer );
  if( !viewer )
    return;
  
  const bool showWarning = InterSpecUser::preferenceValue<bool>( "ExternalRidWarn", viewer );
  
  if( !showWarning )
  {
    auto res = RemoteRid::createDialog( viewer );
    if( callback )
      callback( res.first, res.second );
    
    return;
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
}//void startRemoteRidDialog( viewer, callback )


pair<AuxWindow *, RemoteRid *> RemoteRid::createDialog( InterSpec *viewer )
{
  AuxWindow *window = new AuxWindow( WString::fromUTF8("External RIID"),
                                    (Wt::WFlags<AuxWindowProperties> (AuxWindowProperties::DisableCollapse)
                                     | AuxWindowProperties::SetCloseable
                                     | AuxWindowProperties::IsModal
                                     | AuxWindowProperties::TabletNotFullScreen)
                                    );
  
  WPushButton *close = window->addCloseButtonToFooter();
  close->clicked().connect( window, &AuxWindow::hide );
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
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
  
  WMenu *menu = new WMenu( stack, Wt::Vertical );
  const bool is_phone = m_interspec->isPhone();
  const char *style_class = is_phone ? "VerticalNavMenuPhone HeavyNavMenuPhone SideMenuPhone"
                                     : "VerticalNavMenu HeavyNavMenu SideMenu";
  menu->addStyleClass( style_class );
  
  
  WMenuItem *item = new WMenuItem( "Remote URL", m_rest_rid );
  menu->addItem( item );
  
  item->clicked().connect( std::bind([item,menu](){
    menu->select( item );
    item->triggered().emit( item );
  }) );
  
  item = new WMenuItem( "Executable", m_exe_rid );
  menu->addItem( item );
  item->clicked().connect( std::bind([item,menu](){
    menu->select( item );
    item->triggered().emit( item );
  }) );
  
  
  layout->setContentsMargins( 9, 0, 9, 0 );
  layout->addWidget( menu, 0, 0, AlignLeft );
  
  layout->addWidget( stack, 0, 1 );
  layout->setColumnStretch( 0, 1 );
  
  menu->select( 0 );
  
#else
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setColumnStretch( 0, 1 );
  layout->addWidget( m_rest_rid, 0, 0 );
#endif
}//RemoteRid constructor



std::shared_ptr<SpecUtils::SpecFile> RemoteRid::fileForAnalysis( InterSpec *interspec )
{
  try
  {
    if( !interspec )
      return nullptr;
    
    shared_ptr<const SpecMeas> foreground_file = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    if( !foreground_file )
      return nullptr;
    
    auto answer = make_shared<SpecUtils::SpecFile>( *foreground_file );
    
    if( foreground_file->passthrough() )
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
