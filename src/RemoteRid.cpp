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

#include <string>
#include <sstream>
#include <iostream>

#include <Wt/Utils>
#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStackedWidget>
#include <Wt/WRegExpValidator>

#include "SpecUtils/SpecFile.h"
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
class RestRidWidget;

/** A class to stream the current spectrum data and analysis options to the JavaScript FormData
 format expected by the FullSpectrum service.
 
 Having the client-side javascript make a request to this resource, and then use form data from
 its response is a little round-about, but I dont want to require linking Wt against OpenSSL (or
 another TLS library), so we'll have the client-side browser do the requests, and then send the
 results to this server for display.  Also I chose to use this WResource based method, rather than
 just making hundred of kilobytes (or even megabytes) of JavaScript and then calling
 wApp->doJavaScript(...), because I assume performance will be better, and there is probably some
 practical size limit to what wApp->doJavaScript(...) can handle.
 
 We are declaring up here, and implementing below RestRidWidget so everything is defined when its
 used.
 */
class RestRidInputResource : public Wt::WResource
{
  RemoteRid *m_remote_rid;
  RestRidWidget *m_widget;
  Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
  InterSpec *m_interspec;
  
public:
  RestRidInputResource( RemoteRid *remote_rid, RestRidWidget *widget, InterSpec *interspec, WObject* parent = 0 );
  
  virtual ~RestRidInputResource();
  virtual void handleRequest( const Wt::Http::Request &request,
                             Wt::Http::Response &response );
};//class RestRidInputResource : public Wt::WResource



class RestRidWidget : public Wt::WContainerWidget
{
protected:
  InterSpec *m_interspec;
  RemoteRid *m_remote_rid;
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
  
  JSignal<std::string> m_info_response_signal;
  JSignal<std::string> m_info_error_signal;
  JSignal<std::string> m_analysis_results_signal;
  JSignal<std::string> m_analysis_error_signal;
  
  
  std::string mostRecentUserUrl()
  {
    try
    {
      const string prev_urls = InterSpecUser::preferenceValue<string>( "ExternalRidUrl", m_interspec );
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
      const string prev_urls = InterSpecUser::preferenceValue<string>( "ExternalRidUrl", m_interspec );
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
      
      InterSpecUser::setPreferenceValue( m_interspec->m_user, "ExternalRidUrl", new_url, m_interspec );
    }catch( std::exception &e )
    {
      cerr << "Error in setMostRecentUrl( '" << url << "' ): " << e.what() << endl;
      assert( 0 );
    }// try / catch
  }//void setMostRecentUrl( const string &url )
  
  
public:
  RestRidWidget( RemoteRid *remoterid, InterSpec *interspec, WContainerWidget *parent = 0 )
  : WContainerWidget( parent ),
   m_remote_rid( remoterid ),
   m_interspec( interspec ),
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
   m_info_response_signal( this, "info_response", false ),
   m_info_error_signal( this, "info_error", false ),
   m_analysis_results_signal( this, "analysis_results", false ),
   m_analysis_error_signal( this, "analysis_error", false )
  {
    assert( m_remote_rid );
    assert( interspec );
    
    addStyleClass( "RestRidWidget" );
    
    WGridLayout *layout = new WGridLayout( this );
    
    WLabel *label = new WLabel( "URL" );
    layout->addWidget( label, 0, 0, AlignMiddle );
    m_url = new WLineEdit();
    label->setBuddy( m_url );
    layout->addWidget( m_url, 0, 1, AlignMiddle );
    layout->setColumnStretch( 1, 1 );
    
    const char *url_regex = "^(http(s)?:\\/\\/)[\\w.-]+(?:\\.[\\w\\.-]+)+[\\w\\-\\._~:/?#[\\]@!\\$&'\\(\\)\\*\\+,;=.]+$";
    WRegExpValidator *validator = new WRegExpValidator( url_regex, m_url );
    m_url->setValidator( validator );
    m_url->setPlaceholderText( "URL of FullSpectrum Service" );
    
    m_url->setText( WString::fromUTF8( mostRecentUserUrl() ) );
    
    m_url->changed().connect( this, &RestRidWidget::urlChanged );
    m_url->enterPressed().connect( this, &RestRidWidget::urlChanged );
    

    m_drf_stack = new WStackedWidget();
    m_drf_stack->addStyleClass( "DrfStack" );
    layout->addWidget( m_drf_stack, 1, 0, 1, 2 );
    
    WContainerWidget *container = new WContainerWidget();
    WGridLayout *sub_layout = new WGridLayout( container );
    m_retrieve_drfs_btn = new WPushButton( "Retrieve DRFs" );
    sub_layout->addWidget( m_retrieve_drfs_btn, 0, 0, AlignCenter | AlignMiddle );
    m_retrieve_drfs_btn->clicked().connect( this, &RestRidWidget::requestServiceInfo );
    
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
    m_submit->clicked().connect( this, &RestRidWidget::submitForAnalysis );
    m_submit->disable();
    
    m_resource = new RestRidInputResource( m_remote_rid, this, m_interspec, this );
    
    m_interspec->displayedSpectrumChanged().connect( this, &RestRidWidget::displayedSpectrumChanged );
    
    m_info_response_signal.connect( boost::bind( &RestRidWidget::handleInfoResponse, this, boost::placeholders::_1 ) );
    m_info_error_signal.connect( boost::bind( &RestRidWidget::handleInfoResponseError, this, boost::placeholders::_1 ) );
    m_analysis_results_signal.connect( boost::bind( &RestRidWidget::handleResultResponse, this, boost::placeholders::_1 ) );
    m_analysis_error_signal.connect( boost::bind( &RestRidWidget::handleResultResponseError, this, boost::placeholders::_1 ) );
    
    urlChanged();
  }//RestRidWidget
  
protected:
  
  void submitForAnalysis()
  {
    const string ana_service_url = m_url->text().toUTF8();
    const string &resource_url = m_resource->url();

    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_status->setText( "Requesting analysis" );
    
    
    WStringStream js;
    js <<
    "let sendAnaRequest = function(formData){\n"
    "  fetch('" << ana_service_url << "/analysis', {method: 'POST', body: formData})\n"
    "    .then(response => response.json() )\n"
    "    .then(json_results => {\n"
    "      console.log( 'Got JSON analysis response:', json_results );\n"
    "      const result_str = JSON.stringify(json_results);\n"
    << "      " << m_analysis_results_signal.createCall("result_str") << "\n"
    "    })\n"
    "  .catch( error => {\n"
    "    console.error( 'Error getting analysis input:', error);\n"
    "    const error_str = 'Error retrieving results from server: ' + error.toString();\n"
    "  " << m_analysis_error_signal.createCall( "error_str" ) << ";\n"
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
    << "  " << m_analysis_error_signal.createCall( "error_str" ) << ";\n"
    << "});\n";
    
    doJavaScript( js.str() );
  }//void submitForAnalysis()
  
  
  void requestServiceInfo()
  {
    m_retrieve_drfs_btn->disable();
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "Requesting service information." );
    m_drf_select->clear();
    
    const string url = m_url->text().toUTF8();
    
    WStringStream js;
    
    js << "fetch('" << url << "/info', {method: 'POST'})"
       << ".then( response => response.json() )"
       << ".then( data => {"
//          << "console.log( 'Got data:', data );"
          << "const data_str = JSON.stringify(data);"
          << m_info_response_signal.createCall( "data_str" )
       << "})"
       << ".catch( error => {"
          << "console.error( 'Error retrieving info:', error);"
          << "const error_str = error.toString();"
          << m_info_error_signal.createCall( "error_str" )
       << "});";
    
    doJavaScript( js.str() );
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
    
    Json::Object result;
    try
    {
      Json::parse( msg, result );
    }catch( std::exception &e )
    {
      m_error->setText( "Error parsing analysis results: <code>" + string(e.what()) + "</code>" );
      m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_error) );
    }
    
    try
    {
      Json::Value code = result.contains("code") ? result["code"].toNumber() : Json::Value::Null;
      if( code.isNull() )
        throw runtime_error( "Analysis response JSON did not include a 'code' field." );
      
      string errorMessage = result["errorMessage"].orIfNull("");
      if( !errorMessage.empty() )
        throw runtime_error( "Analysis service returned error message: " +errorMessage );
      
      const double codeVal = code;
      
      /*
       {
       "alarmBasisDuration" : 5025.6201171875,
       "analysisError" : 0,
       "analysisWarnings" : ["background synthesized ..."],
       "chi2" : 1.849143743515015,
       "code" : 0,
       "drf" : "Detective-EX200",
       "errorMessage" : null,
       "foregroundDescription" : "17-Apr-2010 17:33:19, 521.6 &gamma; cps, real time: 5119.3 s",
       "foregroundTitle" : "Foreground",
       "isotopeString" : "Ba133 (H)",
       "isotopes" : [
       {
       "confidence" : 9.47948169708252,
       "confidenceStr" : "H",
       "countRate" : 309.19342041015625,
       "name" : "Ba133",
       "type" : "Industrial"
       }
       ],
       "stuffOfInterest" : 9.47948169708252
       }
       */
      
      cout << Json::serialize(result) << endl;
    
      
      WStringStream rslttxt;
      rslttxt << "<div>\n";
      rslttxt << "<div class=\"ResultLabel\">Results:</div>";
      
      double chi_sqr = -1.0;
      std::vector<std::string> isotope_names;
      std::vector<std::string> isotope_types;
      std::vector<float> isotope_count_rates; //!< If negative, ignore
      std::vector<float> isotope_confidences; //!< If negative, ignore
      std::vector<std::string> isotope_confidence_strs;
      
      //"drf" : "Detective-EX200",
      //"isotopes" : [{"confidence" : 9.4, "confidenceStr" : "H", "countRate" : 309.2, "name" : "Ba133", "type" : "Industrial"}]
      if( result.contains("chi2") )
        chi_sqr = result["chi2"];
      
      if( result.contains("isotopes") )
      {
        const Json::Array &isotopes = result["isotopes"];
        for( size_t i = 0; i < isotopes.size(); ++i )
        {
          const Json::Object &obj = isotopes[i];
          const double isotope_confidence = obj.get("confidence").orIfNull( -1.0 );
          const double isotope_count_rate = obj.get("countRate").orIfNull( -1.0 );
          const string isotope_name = obj.get("name").orIfNull( "null" );
          const string isotope_type = obj.get("type").orIfNull( "null" );
          const string isotope_confidence_str = obj.get("confidenceStr").orIfNull( "null" );
          
          isotope_names.push_back( isotope_name );
          isotope_types.push_back( isotope_type );
          isotope_confidence_strs.push_back( isotope_confidence_str );
          isotope_count_rates.push_back( static_cast<float>(isotope_count_rate) );
          isotope_confidences.push_back( static_cast<float>(isotope_confidence) );
        }//for( size_t i = 0; i < isotopes.size(); ++i )
      }//if( result.contains("isotopes") )
      
      if( codeVal != 0.0 )
        rslttxt << "<div class=\"ResultChi2\">Result Code: " << static_cast<int>(codeVal) << "</div>";
      
      if( (chi_sqr > 0.0001f) && (isotope_names.size() > 0) )
      {
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.2f", chi_sqr );
        rslttxt << "<div class=\"ResultChi2\">&chi;<sup>2</sup>=" << std::string(buffer) << "</div>";
      }//if( output.chi_sqr > 0.0f )
      rslttxt << "</div>\n";
      
      rslttxt << "<table class=\"ResultTable\"><tbody>\n"
      << "\t<tr>"
      << "\t\t<th>Nuclide</th>\n"
      << "\t\t<th>Confidence</th>\n"
      << "\t\t<th>Category</th>\n";
      
      //if( output.ana_type == Analysis::AnalysisType::Simple )
        rslttxt << "\t\t<th>Count Rate</th>\n";
      //else
      //  rslttxt << "\t\t<th>MaxCountRate</th>\n";
      rslttxt << "\t</tr>";
      
      if( isotope_names.empty() )
      {
        rslttxt << "\t<tr>"
        << "\t\t<td colspan=\"4\" style=\"text-align: center; vertical-align: middle;\">None Found</td>\n\t</tr>\n";
      }else
      {
        const size_t nres = isotope_names.size();
        assert( nres == isotope_confidences.size() );
        assert( nres == isotope_types.size() );
        assert( nres == isotope_count_rates.size() );
        assert( nres == isotope_confidence_strs.size() );
        
        for( size_t i = 0; i < nres; ++i )
        {
          const string &iso = isotope_names[i];
          string type = isotope_types[i];
          string conf = isotope_confidence_strs[i];
          const float count_rate = isotope_count_rates[i];
          
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
  
  
  void urlChanged()
  {
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "" );
    m_drf_select->clear();
    m_drf_stack->setCurrentIndex( 0 );
      
    auto spec = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    
    switch( m_url->validate() )
    {
      case Wt::WValidator::Invalid:
      case Wt::WValidator::InvalidEmpty:
        m_submit->hide();
        m_retrieve_drfs_btn->hide();
        m_retrieve_drfs_btn->disable();
        m_status->setText( "Invalid URL" );
        m_submit->disable();
        break;
        
      case Wt::WValidator::Valid:
      {
        const string url = m_url->text().toUTF8();
        
        m_submit->setHidden( !spec || url.empty() );
        m_retrieve_drfs_btn->setHidden( !spec || url.empty() );
        m_retrieve_drfs_btn->setEnabled( spec && !url.empty() );
        
        if( !url.empty() )
          setMostRecentUrl( url );
        
        m_submit->enable();
        break;
      }//case Wt::WValidator::Valid:
    }//switch( m_url->validate() )
  }//urlChanged()
  
  void displayedSpectrumChanged()
  {
    m_result->setText( "" );
    m_error->setText( "" );
    m_status->setText( "" );
    m_status_stack->setCurrentIndex( m_status_stack->indexOf(m_status) );
    
    string url = m_url->text().toUTF8();
    if( m_url->validate() != Wt::WValidator::Valid )
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
  
};//RestRidWidget




RestRidInputResource::RestRidInputResource( RemoteRid *remote_rid, RestRidWidget *widget, InterSpec *interspec, WObject* parent )
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
  shared_ptr<SpecUtils::SpecFile> spec_file = m_remote_rid->fileForAnalysis();
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



#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
class ExeRidWidget : public Wt::WContainerWidget
{
  InterSpec *m_interspec;
  
public:
  ExeRidWidget( InterSpec *interspec, WContainerWidget *parent = 0 )
  : WContainerWidget( parent ),
  m_interspec( interspec )
  {
    addStyleClass( "ExeRidWidget" );
    
    new WText( "This will eventually be the ExeRidWidget", this );
  }
};//class ExeRidWidget
#endif //#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )

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
  window->addHelpInFooter( window->footer(), "external-rid" );
  
  RemoteRid *w = new RemoteRid( viewer, window->contents() );
  
  window->resizeToFitOnScreen();
  window->centerWindowHeavyHanded();
  
  return { window, w };
}//RemoteRid::createDialog(...)


RemoteRid::RemoteRid( InterSpec *viewer, Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_interspec( viewer )
{
  assert( viewer );
  if( !viewer )
    return;
 
  wApp->useStyleSheet( "InterSpec_resources/RemoteRid.css" );
  
  addStyleClass( "RemoteRid" );
  
  RestRidWidget *remote = new RestRidWidget( this, viewer );
  
  WGridLayout *layout = new WGridLayout( this );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setRowStretch( 0, 1 );
  
  
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  ExeRidWidget *exe = new ExeRidWidget( viewer );
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  stack->setOverflow( WContainerWidget::OverflowAuto );
  
  stack->addWidget( remote );
  stack->addWidget( exe );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  WMenu *menu = new WMenu( stack, Wt::Vertical );
  const bool is_phone = m_interspec->isPhone();
  const char *style_class = is_phone ? "VerticalNavMenuPhone HeavyNavMenuPhone SideMenuPhone"
                                     : "VerticalNavMenu HeavyNavMenu SideMenu";
  menu->addStyleClass( style_class );
  
  
  WMenuItem *item = new WMenuItem( "Remote URL", remote );
  menu->addItem( item );
  
  item->clicked().connect( std::bind([item,menu](){
    menu->select( item );
    item->triggered().emit( item );
  }) );
  
  item = new WMenuItem( "Executable", exe );
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
  layout->addWidget( remote, 0, 0 );
#endif
  
  
  
  //RestRidWidget
  
  
  // "ExternalRidExe"
  
}//RemoteRid constructor



std::shared_ptr<SpecUtils::SpecFile> RemoteRid::fileForAnalysis()
{
  try
  {
    if( !m_interspec )
      return nullptr;
    
    shared_ptr<const SpecMeas> foreground_file = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    if( !foreground_file )
      return nullptr;
    
    auto answer = make_shared<SpecUtils::SpecFile>( *foreground_file );
    
    if( foreground_file->passthrough() )
      return answer;
    
    answer->remove_measurements( answer->measurements() );
    
    auto disp_fore = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !disp_fore )
      return nullptr;
    
    auto fore = make_shared<SpecUtils::Measurement>( *disp_fore );
    fore->set_source_type(SpecUtils::SourceType::Foreground);
    
    auto disp_back = m_interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
    
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
