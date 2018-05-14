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

#include <mutex>
#include <string>
#include <stdio.h>

#if( !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )
#include <pwd.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <Lmcons.h>
#endif

#if( defined(WIN32) )
#undef min
#undef max
#endif 

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WBreak>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WIOService>
#include <Wt/WFitLayout>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WProgressBar>
#include <Wt/WBorderLayout>
#include <Wt/WBootstrapTheme>
#include <Wt/WContainerWidget>
#include <Wt/WContainerWidget>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WCartesianChart>
#include <Wt/WMessageResourceBundle>

#if( USE_SRB_HEADER_FOOTER )
#include <boost/filesystem.hpp>
#endif

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  std::mutex AppInstancesMutex;
  std::set<InterSpecApp *> AppInstances;
  //note: could potentially use Wt::WServer::instance()->sessions() to retrieve
  //      sessionIds.
#endif
  
  /*
#if( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  void convert_datetime()
  {
    //A dummy function for windows to cause time_from_string_boost(...)
    //  to initialize locales, to avoid the ~1 second delay it takes to do this.
    try
    {
      UtilityFunctions::time_from_string_boost( "2012-11-23T09:28:36Z" );
    }catch( std::exception & )
    {
    }
  }
#endif
  */
  
#if( ALLOW_URL_TO_FILESYSTEM_MAP && (INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS) )
  std::string uri_decode( const std::string &sSrc )
  {
    //adapted from http://www.codeguru.com/cpp/cpp/algorithms/strings/article.php/c12759/URI-Encoding-and-Decoding.htm
    
    const char HEX2DEC[256] =
    {
      /*       0  1  2  3   4  5  6  7   8  9  A  B   C  D  E  F */
      /* 0 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 1 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 2 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 3 */  0, 1, 2, 3,  4, 5, 6, 7,  8, 9,-1,-1, -1,-1,-1,-1,
      
      /* 4 */ -1,10,11,12, 13,14,15,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 5 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 6 */ -1,10,11,12, 13,14,15,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 7 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      
      /* 8 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* 9 */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* A */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* B */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      
      /* C */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* D */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* E */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1,
      /* F */ -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1
    };
    
    const unsigned char *pSrc = (const unsigned char *)sSrc.c_str();
    const int SRC_LEN = sSrc.length();
    const unsigned char * const SRC_END = pSrc + SRC_LEN;
    // last decodable '%'
    const unsigned char * const SRC_LAST_DEC = SRC_END - 2;
    
    vector<char> newchars( SRC_LEN );
    char *pEnd = &newchars[0];
    
    while (pSrc < SRC_LAST_DEC)
    {
      if (*pSrc == '%')
      {
        char dec1, dec2;
        if (-1 != (dec1 = HEX2DEC[*(pSrc + 1)])
            && -1 != (dec2 = HEX2DEC[*(pSrc + 2)]))
        {
          *pEnd++ = (dec1 << 4) + dec2;
          pSrc += 3;
          continue;
        }
      }
      
      *pEnd++ = *pSrc++;
    }//while (pSrc < SRC_LAST_DEC)
    
    // the last 2- chars
    while (pSrc < SRC_END)
      *pEnd++ = *pSrc++;
    
    return std::string( &newchars[0], pEnd);
  }
#endif  //#if( ALLOW_URL_TO_FILESYSTEM_MAP )
  
#if( USE_SRB_HEADER_FOOTER )
  /** Class meant to show a PHP based site-header for when deployed on the TRB
      WebPage.
   */
  class SrbHeader: public Wt::WContainerWidget
  {
  public:
    SrbHeader::SrbHeader( const std::string &appTitle, const std::string &baseurl,
                         WContainerWidget *parentParam )
    : WContainerWidget(parentParam),
      m_headerText( NULL ),
      m_applicationTitle( NULL )
    {
      static_assert( !(defined(WIN32) || defined(_WIN32)), "PHP header/footers not supported on Windows." );
      
      setId( "SrbHeader" );
      //  setAttributeValue( "style", "" );
      //  setInline(true); //set to <span>
      
      string content = "";
      
      const string srb_common_dir = InterSpecApp::getSrbComonResourcePath();
      
      const string php_command = "php " + srb_common_dir + "/assets/includes/global-nav.php";
      FILE *phpheader = popen( php_command.c_str(), "r" );
      
      if( phpheader != NULL )
      {
        int character;
        while( (character = getc(phpheader)) != EOF ) content += (char)character;
        
        const int close_status = pclose( phpheader );
        
        if( close_status != EXIT_SUCCESS )
        {
          string errmsg = "The php--->c++ srb header status gave a exit status of "
          + std::to_string(close_status)
          + ", you should probably fix this - the php command was '"
          + php_command + "'";
          
          cerr << errmsg << endl; //so it will show up in apaches log
        }//if( close_status != EXIT_SUCCESS )
      } else {
        cerr << SRC_LOCATION << "\n\tFailed to execute '" << php_command << "'\n";
      }//if( opened FILE output ) / else
      
      if( content == "" ) //if for some reason the php command failed, lets default to something
      {
        string errmsg = "InterSpec is NOT using the php-header for some reason; this should be checked into";
        cerr << errmsg << endl;
        
        content =
        "\n  <ul>\n"\
        "    <li><a href=\"https://" + baseurl + "/srb/\" title=\"DHS TRB Webtools Homepage\" id=\"webtools-logo\"><div><span>DHS TRB Webtools Homepage</span></div></a></li> \n"\
        "    <li id=\"select-another-tool\">\n"\
        "      <div id=\"select-header\"><a id=\"select-header-link\">Select another tool</a></div>\n"\
        "      <ul id=\"select-list\" style=\"display:none;\">\n"\
        "        <li><a href=\"https://" + baseurl + "/srb/reports/cbp/\">Automated Reports</a></li>\n"\
        "        <li><a href=\"https://" + baseurl + "/Cambio/\">Cambio</a></li>\n"\
        "        <li><a href=\"https://" + baseurl + "/srb/dhs-server/\">DHS Data Management</a></li>\n"\
        "        <li><a href=\"https://" + baseurl + "/GADRAS\">GADRAS</a></li>\n"\
        "        <li><a href=\"https://" + baseurl + "/srb/viz/lsslogs/\">LSS logs: Google Earth Viz</a></li>\n"\
        "        <li><a href=\"https://" + baseurl + "/wiki/srb/siteinfo/\">POE Wiki</a></li>\n"\
        "      </ul>\n"\
        "    </li>\n"\
        "  </ul>\n\n";
      }//if( content == "" )
      
      m_headerText = new WText( content, Wt::XHTMLUnsafeText, this );
      m_headerText->setId( "global-nav" );
      m_headerText->setInline(false); //set to div
      
      if( appTitle != "" )
      {
        content =
        "\n  <span id=\"site-logo\"></span>\n"\
        "  <h1><a href=\"#\">" + appTitle + "</a></h1>\n\n";
        m_applicationTitle = new WText( content, Wt::XHTMLUnsafeText, this );
        m_applicationTitle->setId( "site-title" );
        m_applicationTitle->setInline(false); //set to div
      } // if( appTitle != "" )
    } // SrbHeader constructor
  protected:
    Wt::WText *m_headerText;
    Wt::WText *m_applicationTitle;
  };//class SrbHeader
  
  /** Class meant to show a PHP based site-footer for when deployed on the TRB
   WebPage.
   */
  class SrbFooter: public Wt::WContainerWidget
  {
  public:
    SrbFooter::SrbFooter( WContainerWidget *parentParam )
    : WContainerWidget(parentParam),
      m_logoText( NULL ),
      m_footerText( NULL ),
      m_footerCloseText( NULL )
    {
      //  setId( "SrbFooter" );
      setId( "footer" );
      setInline(false);
      
      // sponsors
      const string srb_common_url = InterSpecApp::getSrbComonResourceUrl();
      const string srb_common_dir = InterSpecApp::getSrbComonResourcePath();
      string php_dir = srb_common_dir + "/assets/includes/";
      
      string sponsors_html = "";
      string global_sponsors_html = "";
      InterSpecApp::exectutePhpScript( global_sponsors_html, php_dir, "global-sponsors.php", srb_common_url );
      sponsors_html += global_sponsors_html;
      WText *sponsors = new WText( sponsors_html, Wt::XHTMLUnsafeText, this );
      sponsors->setInline( false );
      sponsors->setId( "sponsors" );
      
      // footer
      string global_footer_html = "";
      const int rcode = InterSpecApp::exectutePhpScript( global_footer_html, php_dir, "global-footer.php", srb_common_url );
      if( rcode == EXIT_SUCCESS )
      {
        m_footerText = new WText( global_footer_html, Wt::XHTMLUnsafeText, this );
        m_footerText->setId( "footer-content" );
        m_footerText->setInline(false); //set to div
      }//if( rcode == EXIT_SUCCESS )
      
    }//SrbFooter( constructor )
  protected:
    Wt::WText *m_logoText;
    Wt::WText *m_footerText;
    Wt::WText *m_footerCloseText;
  };//class SrbHeader
#endif //USE_SRB_HEADER_FOOTER
}//namespace


InterSpecApp::InterSpecApp( const WEnvironment &env )
  :  WApplication( env )
    , m_viewer( 0 )
    , m_layout( NULL )
#if(USE_SRB_HEADER_FOOTER)
    , m_header( NULL )
    , m_footer( NULL )
    , m_showHeaderFooterButton( NULL )
    , m_hideHeaderFooterButton( NULL )
#endif //USE_SRB_HEADER_FOOTER
    , m_lastAccessTime( boost::posix_time::microsec_clock::local_time() )
    , m_activeTimeInSession( 0, 0, 0 )
{
  enableUpdates( true );
  
#if( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
  //UtilityFunctions::time_from_string_boost(...) takes an extra second or two
  //  to initialize the first time its called, so we'll do this in the
  //  background. time_from_string_boost() is only used on Windows.
//  if( !loadedSpecFile )
    //WServer::instance()->ioService().post( &convert_datetime );
#endif
  

  //Might as well initialize the DecayDataBaseServer, but in the background
  WServer::instance()->ioService().post( &DecayDataBaseServer::initialize );
 
  setupDomEnvironment();
  setupWidgets( true );
}//InterSpecApp constructor


InterSpecApp::~InterSpecApp()
{
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  {
    std::lock_guard<std::mutex> lock( AppInstancesMutex );
    AppInstances.erase( this );
  }
#endif

}//~InterSpecApp()


void InterSpecApp::setupDomEnvironment()
{
#if( BUILD_AS_ELECTRON_APP )
  //To make nodeIntegration=true work:
  //  https://stackoverflow.com/questions/32621988/electron-jquery-is-not-defined
  doJavaScript( "if (typeof module === 'object') {window.module = module; module = undefined;}", false );

#if( USE_ELECTRON_NATIVE_MENU )
  //Some support to use native menu...
  doJavaScript( "const {remote} = require('electron');const {Menu, MenuItem} = remote;", false );
#endif
  
#endif
  
  setTitle( "TRB InterSpec" );
  
  //Call tempDirectory() to set global variable holding path to the temp file
  //  directory, to ensure this will be avaialble at all points in the future;
  //  this is a workaround for FCGI deployment.
  tempDirectory();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  //Setting this threshold doesnt seem to make a huge difference, except
  //  for on Windows using websockets, it fixes an issue where pre-loaded
  //  tabs and menus wont otherwise work. (default is 5000 {see WebRenderer.C}
  //  despite what Wt documentation says)
  const int threshold_bytes = 20*1024*1024;
  setTwoPhaseRenderingThreshold( threshold_bytes );
#endif
  
  //XXX - it appears Firefox/Opera doesnt always load the wt.css file, lets force it
  //  if( env.agentIsGecko() || env.agentIsOpera() )
  //const string theme = "default"; //cssTheme();
  //useStyleSheet( WApplication::resourcesUrl() + "/themes/" + theme + "/wt.css");
  setCssTheme( "default" );
  
  //Ethan added this to test bootstrap looks
  //  WBootstrapTheme *theme = new WBootstrapTheme(this);
  //  theme->setVersion(  Wt::WBootstrapTheme(3) );
  //  setTheme( theme );
  
#if(USE_SRB_HEADER_FOOTER)
  const string srb_common_url = getSrbComonResourceUrl();
  useStyleSheet( srb_common_url + "/assets/css/styles.css" );
#endif
  
  //Adding css/js files for pNotify popup notifications
  //TODO: remove pNotify if we use qTip2
  requireJQuery("InterSpec_resources/assets/js/jquery-1.11.1.min.js");
  //  require("InterSpec_resources/assets/js/pnotify.custom.min.js");
  //  useStyleSheet( "InterSpec_resources/assets/css/pnotify.custom.min.css" );
  //  useStyleSheet( "InterSpec_resources/assets/css/bootstrap.css" );
  //  doJavaScript("PNotify.prototype.options.styling = \"bootstrap3\";");
  
  //for qTip2
  useStyleSheet( "InterSpec_resources/assets/css/jquery.qtip.min.css" );
  require("InterSpec_resources/assets/js/jquery.qtip.min.js");
  require("InterSpec_resources/assets/js/imagesloaded.pkg.min.js");
  
  
  useStyleSheet( "InterSpec_resources/InterSpec.css" );
  doJavaScript( "if(typeof console==='undefined'){console={log:function(){}};}" );
  
  
  styleSheet().addRule( "input[type=\"text\"]", "font-size:0.95em;" );
  styleSheet().addRule( "button[type=\"button\"]", "font-size:0.9em;" );
  
  //This is needed to fix keyboard hiding in iOS
#if(IOS || ANDROID)
  const char *fix_ios_js = INLINE_JAVASCRIPT(
    var t=document.createElement('meta');
    t.name = "viewport";
    t.content = "initial-scale=1.0, user-scalable=no, height=device-height, width=device-width";
    document.getElementsByTagName('head')[0].appendChild(t);
    $(document).on('blur', 'input, textarea', function() {
      setTimeout(function(){ window.scrollTo(document.body.scrollLeft, document.body.scrollTop); }, 0);
    });
  );
  
  doJavaScript( fix_ios_js );
#endif
  
  if( isMobile() )
  {
    //Prevent mobile from showing spinner
    const char *prevent_spinner_js = INLINE_JAVASCRIPT(
      $("<style>").prop("type", "text/css")
                  .html(".Wt-spinbox {background-image: inherit !important; \
                         background-repeat: inherit  !important; \
                         background-position: inherit  !important; \
                         padding-right: inherit  !important; }")
                  .appendTo("head");
    );
    
    cerr << "prevent_spinner_js='" << prevent_spinner_js << "'" << endl;
    
    //Prevent mobile from hovering white
    const char *prevent_hovering_white = INLINE_JAVASCRIPT(
      $("<style>").prop("type", "text/css")
                  .html("ul.Wt-popupmenu > li > a > span > span:hover {color: black !important;}")
                  .appendTo("head");
    );
    
    doJavaScript( prevent_spinner_js );
    doJavaScript( prevent_hovering_white );
  }//if( isMobile() )
  
#if( BUILD_AS_ELECTRON_APP )
  doJavaScript( "if (window.module) module = window.module;", false );
#endif
}//void setupDomEnvironment()


void InterSpecApp::setupWidgets( const bool attemptStateLoad  )
{
  if( m_viewer )
  {
    //for some reason calling 'delete m_layout' doesnt cause m_viewer to
    //  destruct?
    delete m_viewer;
    m_viewer = 0;
  }
  
#if(USE_SRB_HEADER_FOOTER)
  if( m_header )
  {
    delete m_header;
    m_header = 0;
  }
  
  if( m_footer )
  {
    delete m_footer;
    m_footer = 0;
  }
#endif
  
  if( m_layout )
  {
    delete m_layout;
    root()->clear();
  }
  
  try
  {
    m_viewer = new InterSpec();
  }catch( std::exception &e )
  {
#if( BUILD_AS_UNIT_TEST_SUITE )
    cerr << "There was an expection initializing a InterSpec: " << e.what() << endl;
    throw e;
#endif //#if( BUILD_AS_UNIT_TEST_SUITE )
    
    string msg = "There was a problem initializing a necessary resource"
    " for InterSpec";
    WText *errorText = new WText( msg, root() );
    errorText->setInline( false );
    errorText->setAttributeValue( "style",
                                 "width:100%;"
                                 "text-align:center;"
                                 "font-size:150%;"
                                 "font-family:\"Courier New\", Courier, monospace;" );
    msg = "error:";
    errorText = new WText( msg, root() );
    errorText->setAttributeValue( "style", "font-family:\"Courier New\", Courier, monospace;" );
    
    WContainerWidget *errorDiv = new WContainerWidget( root() );
    errorText = new WText( e.what(), errorDiv );
    errorText->setAttributeValue( "style", "font-family:\"Times New Roman\", Times, serif;" );
    
    new WBreak( root() );
    
    msg = "Please contact wcjohns@sandia.gov and/or srb@sandia.gov to fix this error.";
    errorText = new WText( msg, root() );
    errorText->setAttributeValue( "style", "font-family:Courier New; color: blue;" );
    WApplication::quit();
    return;
  }//try / catch to create a InterSpec object
  
  root()->addStyleClass( "specviewer" );
  
  m_layout = new WGridLayout();
  root()->setLayout( m_layout );
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setHorizontalSpacing( 0 );
  m_layout->setVerticalSpacing( 0 );
  
  
#if(!USE_SRB_HEADER_FOOTER)
  m_layout->addWidget( m_viewer, 0, 0, 1, 1 );
  m_layout->setRowStretch( 0, 10 );
#else
  PopupDivMenu *displayOptions = m_viewer->displayOptionsPopupDiv();
  
  displayOptions->addSeparator();
  m_showHeaderFooterButton = displayOptions->addMenuItem( "Show TRB Header and Footer" );
  m_hideHeaderFooterButton = displayOptions->addMenuItem( "Hide TRB Header and Footer" );
  m_showHeaderFooterButton->triggered().connect( boost::bind( &InterSpecApp::showHeaderFooter, this, true ) );
  m_hideHeaderFooterButton->triggered().connect( boost::bind( &InterSpecApp::hideHeaderFooter, this, true ) );
  m_hideHeaderFooterButton->hide();
  
  string urlLocal = "hekili.ca.sandia.gov";  //url(); ?? wApp->url(); ??
  m_header = new SrbHeader( "", urlLocal );
  m_footer = new SrbFooter();
  m_layout->addWidget( m_header, 0, 0, 1, 1 );
  m_layout->addWidget( m_viewer, 1, 0, 1, 1 );
  m_layout->addWidget( m_footer, 2, 0, 1, 1 );
  m_layout->setRowStretch( 0, 0 );
  m_layout->setRowStretch( 1, 10 );
  m_layout->setRowStretch( 2, 0 );
  
  if( m_viewer->m_user->preferenceValue<bool>( "ShowHeaderFooter" ) )
    showHeaderFooter( false );
  else
    hideHeaderFooter( false );
#endif //USE_SRB_HEADER_FOOTER
  
  
  m_clearSession.reset( new JSignal<>(this, "clearSession", false) );
  m_clearSession->connect( this, &InterSpecApp::clearSession );
  
  bool loadedSpecFile = false;
  const Http::ParameterMap &parmap = environment().getParameterMap();
  
#if( ALLOW_URL_TO_FILESYSTEM_MAP )
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS )
  //Allowing opening a file via the URL could potentially be a security issue
  //  (not entirely sure how, but JIC), so will restrict this feature to
  //  development builds only (outside of development we dont need this anyway).
  const Http::ParameterMap::const_iterator specfileiter
  = parmap.find( "specfilename" );
  if( (specfileiter != parmap.end()) && specfileiter->second.size() )
  {
    const string filename = uri_decode( specfileiter->second[0] );
    loadedSpecFile = m_viewer->userOpenFileFromFilesystem( filename );
    if( loadedSpecFile )
      cerr << "Openend file specified by URL '" << filename << "'" << endl;
    else
      cerr << "Invalid specfile specified in URL '" << filename << "'" << endl;
  }//if( speciter != parmap.end() )
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS )

  
  const Http::ParameterMap::const_iterator speciter = parmap.find( "specfile" );
  if( (speciter != parmap.end()) && speciter->second.size() )
  {
    try
    {
      const int id = std::stoi( speciter->second[0] );
      std::cerr << "SpecViewerApp: Will try to open file with ID="
      << id << " (" << speciter->second[0] << ")" << std::endl;
      loadedSpecFile = m_viewer->loadFileFromDbFilesystemLink( id, false );
    } catch( std::exception & )
    {
      cerr << "Invalid specfile" << endl;
    }//try / catch
  }//if( speciter != parmap.end() )
#endif  //#if( ALLOW_URL_TO_FILESYSTEM_MAP )


#if( USE_DB_TO_STORE_SPECTRA )
  //Check to see if we should load the apps last saved state
  try
  {
    const bool saveSpectra =  true; // InterSpecUser::preferenceValue<bool>("SaveSpectraToDb", m_viewer );
    bool saveState = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb" /* "SaveStateToDbOnExit"*/, m_viewer );
    saveState = (saveState && attemptStateLoad);
    
#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
    //Make it so the app will prompt user to load the state
    //  - Experimental!
    //  - maybe should make it so it always does this, except on iOS and Android
    const bool loadPrev = InterSpecUser::preferenceValue<bool>("LoadPrevStateOnStart", m_viewer );
    const bool promptLoad = InterSpecUser::preferenceValue<bool>("PromptStateLoadOnStart", m_viewer );
    if( !loadPrev && !promptLoad )
      saveState = false;
#endif
    
    //if the URL contains something like "&restore=no", then we wont load in
    //  stored state
    if( saveState && saveSpectra )
    {
      Http::ParameterMap::const_iterator iter = parmap.find( "restore" );
      saveState = !(iter != parmap.end() && iter->second.size()
                    && (iter->second[0] == "0" || iter->second[0] == "no"
                        || iter->second[0] == "false"));
    }//if( saveState && saveSpectra )
    
    if( !loadedSpecFile && saveSpectra && saveState )
    {
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
      DataBaseUtils::DbTransaction transaction( *sql );
      
      const string query = "InterSpecUser_id = "
      + std::to_string( m_viewer->m_user.id() )
      + " AND (StateType = "
      + std::to_string( int(UserState::kEndOfSessionTemp) )
      + " OR StateType = "
      + std::to_string( int(UserState::kEndOfSessionHEAD) ) + ")";
      Dbo::ptr<UserState> state = sql->session()->find<UserState>()
      .where( query ).orderBy( "id desc" ).limit( 1 );
      
      //Since we only read from the database, its probably fine if the commit
      //  fails (with SQLite3 this will happen if the database is locked because
      //  of another write to it is happening), so we'll allow this.
      try
      {
        transaction.commit();
      }catch( std::exception &e )
      {
        cerr << "Failed to commit transaction while loading state, caught: "
        << e.what() << endl;
      }//try / catch
      
      if( state )
      {
        loadedSpecFile = true;
#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
        if( promptLoad )
        {
          //Create a dialog asking if the user wants to pick up
          AuxWindow *loadStateDialog = new AuxWindow( "Load previous state?" ,true);
        
          WContainerWidget *dialogDiv = loadStateDialog->contents();
          
          WLabel *text = new WLabel( "Would you like to pick up where you left off last time?" );
          text->setInline( false );
          text->setMargin( 15 );
          dialogDiv->addWidget( text );
          
          WCheckBox *cb = new WCheckBox( "Dont ask again" );
          cb->setInline( false );
          cb->setMargin( 30, Wt::Left );
          cb->setMargin( 10, Wt::Bottom );
          dialogDiv->addWidget( cb );
          auto changePromptPref = [this](bool dont_ask ){
            try {
              InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user,
                              "PromptStateLoadOnStart", dont_ask, m_viewer );
            }catch( std::exception &e ) {
              cerr << "Failed setting 'PromptStateLoadOnStart' user prefrence" << endl;
            }
          };
        
          auto changeDoLoadPref = [this](bool load ){
            try {
              InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user,
                                    "LoadPrevStateOnStart", load, m_viewer );
            }catch( std::exception &e ) {
              cerr << "Failed setting 'LoadPrevStateOnStart' user prefrence" << endl;
            }
          };
        
          
          cb->unChecked().connect( std::bind([changePromptPref](){ changePromptPref(false); }) );
          cb->checked().connect( std::bind([changePromptPref](){ changePromptPref(true); }) );
          
          WPushButton *nobutton = loadStateDialog->addCloseButtonToFooter("No");
          
          nobutton->clicked().connect( std::bind([this,cb,loadStateDialog,changeDoLoadPref](){
            if( cb->isChecked() )
              changeDoLoadPref(false);
            AuxWindow::deleteAuxWindow(loadStateDialog);
          }) );
           
        
          WPushButton *yesbutton = loadStateDialog->addCloseButtonToFooter("Yes");
          
          yesbutton->clicked().connect( std::bind([this,cb,state,loadStateDialog,changeDoLoadPref](){
            if( cb->isChecked() )
              changeDoLoadPref(true);
            m_viewer->loadStateFromDb( state );
            AuxWindow::deleteAuxWindow(loadStateDialog);
          }));
           
        
          if( loadPrev )
            yesbutton->setFocus();
          else
            nobutton->setFocus();
          
          loadStateDialog->centerWindow();
          loadStateDialog->disableCollapse();
          loadStateDialog->rejectWhenEscapePressed();
          loadStateDialog->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, loadStateDialog ) );
          loadStateDialog->show();
        }else
        {
#endif
        m_viewer->loadStateFromDb( state );
        
        WStringStream js;
        js << "<div onclick="
        "\"Wt.emit('" << id() << "',{name:'clearSession'});"
        "$('.qtip.jgrowl:visible:last').remove();"
        "return false;\" "
        "class=\"clearsession\"><span class=\"clearsessiontxt\">Start Fresh Session</span></div>";
        
        passMessage( "Resuming where you left off on " + state->name + js.str(),
                    "", WarningWidget::WarningMsgInfo );
#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
        }
#endif
      }else
      {
        cerr << "Could not load previously saved state" << endl;
      }
    }//if( saveSpectra && saveState )
  }catch( std::exception &e )
  {
    cerr << "Failed to load app state, caught: " << e.what() << endl;
  }//try / catch
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  
#if( USE_DB_TO_STORE_SPECTRA && (ALLOW_URL_TO_FILESYSTEM_MAP || INCLUDE_ANALYSIS_TEST_SUITE ) )
  bool userstate = true;
  Http::ParameterMap::const_iterator stateiter = parmap.end();
  
#if( ALLOW_URL_TO_FILESYSTEM_MAP )
  stateiter = parmap.find( "appstate" );
#endif
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  userstate = (stateiter != parmap.end());
  if( !userstate )
    stateiter = parmap.find( "teststate" );
#endif
  
  if( !loadedSpecFile && stateiter!=parmap.end() && stateiter->second.size() )
  {
    const int id = std::stoi( stateiter->second[0] );
    const int type = int(userstate ? UserState::kUserSaved : UserState::kForTest);
    char query[128];
    snprintf( query, sizeof(query), "StateType=%i AND id=%i", type, id );
    
    
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    
    Dbo::ptr<UserState> state = sql->session()->find<UserState>().where( query );
    transaction.commit();
    std::cerr << "SpecViewerApp: Will try to open test state ID="
    << id << " (" << id << ")" << std::endl;
    
    if( state )
    {
      m_viewer->loadStateFromDb( state );
      loadedSpecFile = true;
    }else
    {
      cerr << "Could not load test state" << endl;
    }//if( state ) / else
  }//if( (stateiter != parmap.end()) && stateiter->second.size() )
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )

  
  auto showWelcomeCallback = [this,loadedSpecFile](){
    //If we already loaded a spectrum, then dont show welcome dialog (or IE
    //  warning dialog)
    if( loadedSpecFile )
      return;
    
    //If client is internet explorer, show a warning before the welcome dialog
    if( !environment().agentIsIE() )
    {
      m_viewer->showWelcomeDialog();
    }else
    {
      AuxWindow *dialog = m_viewer->showIEWarningDialog();
      if( dialog )
        dialog->finished().connect( boost::bind( &InterSpec::showWelcomeDialog, m_viewer, false ) );
      else
        m_viewer->showWelcomeDialog();
    }// if( not in IE ) / else
  };//auto showWelcomeCallback
  
  //Need to check if this is the latest version of InterSpec this user has used;
  //  if not, need to set "HasAgreedToUseTerms" to false.
  
  //  Check if user has agreed to terms, if not, show terms window
  const int previousAgreedVersion = InterSpecUser::preferenceValue<int>( "VersionAgreedToUseTerms" , m_viewer );
  if( previousAgreedVersion != COMPILE_DATE_AS_INT )
  {
    m_viewer->showLicenseAndDisclaimersWindow( true, showWelcomeCallback );
  }else if( !loadedSpecFile )
  {
    showWelcomeCallback();
  }//if( specfile.empty() )
  
  //Only display warning message on navigation away from current session if in
  //  a browser environment, and not this isnt a developers build
#if( (BUILD_FOR_WEB_DEPLOYMENT || BUILD_AS_LOCAL_SERVER) && !PERFORM_DEVELOPER_CHECKS  )
  
#if( !OPTIMISTICALLY_SAVE_USER_STATE )
  setConfirmCloseMessage( "This action will terminate your current session." );
#else
  //The bellow solution of catching onbeforeunload catches navigating to a new
  //  url, or closing a tab - but for refreshes, the ThinkOfLeave signal is
  //  never transmitted back to the server...
  m_thinkingLeaveSignal.reset( new Wt::JSignal<>( this, "ThinkOfLeave" ) );
  m_thinkingLeaveSignal->connect( this, &InterSpecApp::prepareForEndOfSession );
  
  const char *closejs = INLINE_JAVASCRIPT
  (
   window.onbeforeunload = function(event) {
     var e = event || window.event;
     var txt = "This action will terminate your current session.";
     if(e)
       e.returnValue = txt;
     //    console.log( 'User tried navigating, refreshing, or closing' );
     //       window.onunload();
     Wt.emit( $(document).data('AppId'),{name:'ThinkOfLeave'} );
     //       Wt._p_.update( $(document).data('AppId'), 'ThinkOfLeave', {name:'ThinkOfLeave'}, false ); //doesnt seem to help, but I really dont know how to call it
     return txt;
   }
   );
  doJavaScript( "$(document).data('AppId', '" + id() + "');" );
  doJavaScript( closejs );
#endif
#endif
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  {
    std::lock_guard<std::mutex> lock( AppInstancesMutex );
    AppInstances.insert( this );
  }
#endif
  
  if( m_viewer->m_user )
    cerr << "Have started session " << sessionId() << " for user "
    << m_viewer->m_user->userName() << endl;
}//void setupWidgets()



std::string InterSpecApp::getUserNameFromEnvironment() const
{
  string remoteUser = environment().getCgiValue( "REMOTE_USER" );
  UtilityFunctions::ireplace_all( remoteUser, " ", "" );
  UtilityFunctions::to_lower( remoteUser );
  return remoteUser;
}//std::string getUserNameFromEnvironment() const


std::string InterSpecApp::getSrbComonResourceUrl()
{
#if( !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )
  char hostname[256];
  const int namestatus = gethostname( hostname, 255 );
  if( namestatus != 0 )
    hostname[0] = '\0';
  
  if( strstr(hostname, "outrage") || strstr(hostname, "hekili") )
  {
    const string url = wApp->url();
    
    //get the leading directory deployment
    string userDir = "";
    vector< string > splitVec;
    UtilityFunctions::split( splitVec, url, "/" );
    
    if( splitVec.size() ) userDir = splitVec[0];
    //if url started with a '/', then splitVec[0]==""
    if( userDir.empty() && (splitVec.size() > 1) ) userDir = splitVec[1];
    
    return wApp->makeAbsoluteUrl( "/" + userDir + "/common" );
  }//if( on outrage )
#endif
  return "srb_common_for_not_on_outrage";
}//std::string getSrbComonResourceUrl()


std::string InterSpecApp::getSrbComonResourcePath()
{
#if( !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )
  char hostname[256];
  const int namestatus = gethostname( hostname, 255 );
  if( namestatus != 0 )
    hostname[0] = '\0';
  
  if( strstr(hostname, "outrage") || strstr(hostname, "hekili") )
  {
    //ex. '/~srb/fastcgi/wcjohns/app.fcgi?wtd=...'
    const string url = wApp->url();
    
    if( UtilityFunctions::starts_with(url, "/srb/" )     )
      return "/home/srb-prod/public_html/common";
    if( UtilityFunctions::starts_with(url, "/~srb/" )    )
      return "/home/srb/public_html/common";
    if( UtilityFunctions::starts_with(url, "/srb-dev/" )
       || UtilityFunctions::starts_with(url, "/~srb-dev/" ) )
      return "/home/srb-dev/public_html/common";
    if( UtilityFunctions::starts_with(url, "/srb-qc/" )
       || UtilityFunctions::starts_with(url, "/~srb-qc/" ) )
      return "/home/srb-qc/public_html/common";
    if( UtilityFunctions::starts_with(url, "/srb-prod/" )
       || UtilityFunctions::starts_with(url, "/~srb-prod/" ) )
      return "/home/srb-prod/public_html/common";
    
    string errmsg = "The InterSpecApp gui is being executed at an unknown url='"
    + wApp->url()
    + "'; this should be checked into";
    cerr << errmsg << endl; //so it makes it into the apache log files
  }//if( on outrage )
#endif
  return "srb_common_for_not_on_outrage";
}//std::string getSrbComonResourcePath()


std::string InterSpecApp::tempDirectory()
{
  static std::mutex s_tempDirMutex;
  static std::string s_tempDir;
  
  std::lock_guard<std::mutex> lock( s_tempDirMutex );
  
  if( s_tempDir.size() )
    return s_tempDir;
  
  WApplication *app = WApplication::instance();
  
  if( app )
  {
    s_tempDir = app->environment().getCgiValue("TMPDIR");
    if( s_tempDir.empty() )
      s_tempDir = app->environment().getCgiValue("TMP");
    if( s_tempDir.empty() )
      s_tempDir = app->environment().getCgiValue("TEMP");
    if( s_tempDir.empty() )
      s_tempDir = app->environment().getCgiValue("TEMPDIR");
    
    //Must not be spcified as a CGI value, lets just go with what the operating
    //  system says
    if( s_tempDir.empty() || !UtilityFunctions::is_directory(s_tempDir) )
      s_tempDir = UtilityFunctions::temp_dir();
    
    return s_tempDir;
  }//if( app )
  
  return UtilityFunctions::temp_dir();
}//void tempDirectory()



std::string InterSpecApp::userNameFromOS()
{
  const char *sysusername = getenv( "USER" );  //*NIX
  if( sysusername )
    return sysusername;
  
  if( !sysusername )
    sysusername = getenv( "USERNAME" );        //Windows
  if( sysusername )
    return sysusername;
  
#if( !(defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64)) )
  //I havent tested the following code
  const uid_t uid = geteuid();
  const struct passwd *pw = getpwuid( uid );
  if( pw )
    return pw->pw_name;
#else
  TCHAR username[UNLEN+1];
  DWORD size = UNLEN + 1;
  GetUserName( username, &size );
  return username;
#endif
  
  return "";
}//std::string userNameFromOS()


#if USE_SRB_HEADER_FOOTER
int InterSpecApp::exectutePhpScript( std::string &results,
                      const std::string &path_to_script,
                      const std::string &script_name,
                      const string &path_to_srb_common )
{
  results = "";
  
  boost::filesystem::path script = script_name;
  boost::filesystem::path script_location = path_to_script;
  
  if( path_to_script.empty() )
    script_location = ".";
  
  boost::system::error_code ec;
  script = boost::filesystem::path(script_location) / script;
  script = boost::filesystem::absolute( script );
  

  if( !boost::filesystem::is_regular_file( script, ec ) )
  {
    string errmsg = "Couln't access the file '" + script.generic_string()  + "'";
    cerr << errmsg << endl; //so it will show up in apaches log
    return EXIT_FAILURE;
  }//if( we can touch the file ).
  
  const string host = wApp->environment().hostName();
  const string user = static_cast<InterSpecApp *>(wApp)->getUserNameFromEnvironment();
  string uri = host;
  const size_t ind = uri.find_first_of( "/" );
  if( ind != string::npos ) uri = uri.substr( ind );
  else                      uri = "";
  
  //TODO: It doesnt appear the below is actually working!
  const string file_contents = "<?php $PATH_TO_COMMON=\""
  + path_to_srb_common +
  "/\"; "
  "$HTTP_HOST=\""     + host + "\"; "  //Contents of the Host: header from the current request, if there is one.
  "$PHP_AUTH_USER=\"" + user + "\"; "  //When doing HTTP authentication this variable is set to the username provided by the user.
  "$REQUEST_URI=\""   + uri  + "\"; "  //The URI which was given in order to access this page; for instance, '/index.html'.
  "$_SERVER['" + host + "']=\"" + host + "\"; "
  "$_SERVER['" + user + "']=\"" + user + "\"; "
  "$_SERVER['" + uri + "']=\"" + uri + "\"; "
  " include_once( '"
  + script.generic_string()
  + "' ); \?>";
  
  
  string tempFileName = UtilityFunctions::temp_file_name( "IntersSpecPhpScript", InterSpecApp::tempDirectory() );
  FILE *tempFile = fopen( tempFileName.c_str(), "w+" );
  
  if( tempFile == NULL )
  {
    string errmsg = "exectutePhpScript(...): couldnt open temp file";
    cerr << errmsg << endl;
    return EXIT_FAILURE;
  }//if( tempFile == NULL )
  
  fprintf( tempFile, "%s", file_contents.c_str() );
  
  fflush( tempFile );
  //  ::pclose( tempFile );
  fclose( tempFile );
  
  string php_command = "php " + tempFileName;
  
  FILE *phpheader = popen( php_command.c_str(), "r" );
  
  int close_status = EXIT_FAILURE;
  
  if( phpheader != NULL )
  {
    int character;
    while( (character = getc(phpheader)) != EOF ) results += (char)character;
    
    close_status = pclose( phpheader );
    
    if( close_status != EXIT_SUCCESS )
    {
      string errmsg = "The php--->c++/tex status gave a exit status of "
      + std::to_string(close_status)
      + ", you should probably fix this '"
      + php_command + "'";
      cerr << errmsg << endl; //so it will show up in apaches log
    }//if( close_status != EXIT_SUCCESS )
  }else
  {
    string errmsg = "The InterSpecApp failed to execute '" + php_command + "'";
    cerr << errmsg << endl;
  }//if( opened FILE output ) / else
  
  UtilityFunctions::remove_file( tempFileName );
  
  return close_status;
}//int exectutePhpScript( const std::string &command, std::string &result )
#endif // USE_SRB_HEADER_FOOTER

InterSpec *InterSpecApp::viewer()
{
  return m_viewer;
}//InterSpec* viewer()


#if( !BUILD_FOR_WEB_DEPLOYMENT )
std::string InterSpecApp::sessionUrlId()
{
  WApplication::UpdateLock lock( this );

  const Http::ParameterMap &parmap = environment().getParameterMap();
  for( const Http::ParameterMap::value_type &p : parmap )
  {
    if( UtilityFunctions::iequals(p.first, "externalid") && !p.second.empty() )
      return p.second.front();
  }//for( const Http::ParameterMap::value_type &p : parmap )
  
  return "";
}//const std::string &sessionUrlId() const


bool InterSpecApp::userOpenFromFileSystem( const std::string &path )
{
  if( !m_viewer )
    return false;
  return m_viewer->userOpenFileFromFilesystem( path );
}//bool userOpenFromFileSystem(...)


#if( ALLOW_URL_TO_FILESYSTEM_MAP )
bool InterSpecApp::openFileFromDbFileSystemLink( int index )
{
  return m_viewer->loadFileFromDbFilesystemLink( index, true );
}//bool openFileFromDbFileSystemLink( int index );
#endif

std::set<InterSpecApp *> InterSpecApp::runningInstances()
{
  std::lock_guard<std::mutex> lock( AppInstancesMutex );
  return AppInstances;
}//set<InterSpecApp *> runningInstances()


InterSpecApp *InterSpecApp::instanceFromExtenalIdString( const std::string &idstr )
{
  if( idstr.empty() )
    return (InterSpecApp *)0;
  
  std::lock_guard<std::mutex> lock( AppInstancesMutex );
  for( InterSpecApp *app : AppInstances )
  {
    Wt::WApplication::UpdateLock lock( app );
    if( app->sessionUrlId() == idstr )
      return app;
  }
  
  return (InterSpecApp *)0;
}//InterSpecApp *instanceFromExtenalIdString( const std::string &idstr )

#endif //#if( !BUILD_FOR_WEB_DEPLOYMENT )

#if( defined(WIN32) && BUILD_AS_ELECTRON_APP )
  //When users drag files from Outlook on windows into the app
  //  you can call the following functions
void InterSpecApp::dragEventWithFileContentsStarted()
{
  m_viewer->dragEventWithFileContentsStarted();
}

void InterSpecApp::dragEventWithFileContentsFinished()
{
  m_viewer->dragEventWithFileContentsFinished();
}
#endif

#if( BUILD_AS_ELECTRON_APP )
bool InterSpecApp::isElectronInstance()
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( !app )
	return false;
  WApplication::UpdateLock lock(app);
  if( !lock )
    return false;

  const string externalid = app->sessionUrlId();
  
  //XXX - should compare externalid to ElectronUtils::external_id()

  return !externalid.empty();
}
#endif

void InterSpecApp::notify( const Wt::WEvent& event )
{
//  cout << "InterSpecApp::notify starting" << endl;
#if( !BUILD_AS_UNIT_TEST_SUITE )
  try
  {
#endif
    const bool userEvent = (event.eventType() == Wt::UserEvent);
    if( userEvent )
    {
      using namespace boost::posix_time;
      const ptime thistime = microsec_clock::local_time();
      time_duration duration = thistime - m_lastAccessTime;
      if( duration < minutes(5) )
        m_activeTimeInSession += duration;
      m_lastAccessTime = thistime;
    }//if( userEvent )

     WApplication::notify( event );
    
    //Note that event.eventType() may have change (although I dont know how/why)
//    if( userEvent )
//    {
//      std::shared_ptr<const SpecMeas> meas = m_viewer->measurment(kForeground);
//      std::shared_ptr<const SpecMeas> back = m_viewer->measurment(kBackground);
//      std::shared_ptr<const SpecMeas> secn = m_viewer->measurment(kSecondForeground);
//      
//      if( (meas&&meas->modified())
//          || (back&&back->modified()) || (secn&&secn->modified())
//         //or the spectrum for any of these changed
//         )
//      {
//        //save the incremental difference to the database
//      }//
//    }//if( userEvent )
#if( !BUILD_AS_UNIT_TEST_SUITE )
  }catch( std::exception &e )
  {
    string msg = "There was an unexpected error, application state may be suspect: ";
    msg += e.what();
    svlog( msg, "", WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    char message[1024];
    snprintf(message, sizeof(message), "Uncaught exception in event loop: %s", e.what() );
    log_developer_error( BOOST_CURRENT_FUNCTION, message );
#endif
    
  }//try/catch
#endif //#if( !BUILD_AS_UNIT_TEST_SUITE )

//  cout << "\tending InterSpecApp::notify" << endl;
}//void notify( const Wt::WEvent& event )


void InterSpecApp::unload()
{
  cerr << "\n\nRecieved unload from " << sessionId() << endl;
  WApplication::unload();
}//void unload()


void InterSpecApp::prepareForEndOfSession()
{
  try
  {
    if( m_viewer && m_viewer->m_user )
    {
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
      
      {
        DataBaseUtils::DbTransaction transaction( *sql );
        
        //If the user has mutliple sessions going on, the bellow will quash
        //  eachother out.  Could probably be fixed by
        //  m_viewer->m_user.reread();
        
        m_viewer->m_user.modify()->addUsageTimeDuration( m_activeTimeInSession );
        cerr << "Added " << m_activeTimeInSession << " usage time for user "
        << m_viewer->m_user->userName() << endl;
        m_activeTimeInSession = boost::posix_time::time_duration(0,0,0,0);
        transaction.commit();
      }
      
#if( USE_DB_TO_STORE_SPECTRA )
      //Check to see if we should save the apps state
      const bool saveSpectra = true; //InterSpecUser::preferenceValue<bool>( "SaveSpectraToDb", m_viewer );
      const bool saveState = InterSpecUser::preferenceValue<bool>(
                                                                  "AutoSaveSpectraToDb" /*"SaveStateToDbOnExit"*/, m_viewer );
        
        //Clean up the kEndOfSessions from before
        bool cleanupStates = true;
        if( cleanupStates )
        {
            DataBaseUtils::DbTransaction transaction( *sql );
            const string query = "InterSpecUser_id = "
            + std::to_string( m_viewer->m_user.id() )
            + " AND (StateType = "
            + std::to_string( int(UserState::kEndOfSessionTemp))
            + " OR StateType = "
            + std::to_string( int(UserState::kEndOfSessionHEAD) )
            + ")";
            Dbo::collection<Dbo::ptr<UserState> > states
            = sql->session()->find<UserState>().where( query );
            for( Dbo::collection<Dbo::ptr<UserState> >::iterator iter = states.begin();
                iter != states.end(); ++iter )
            {
                if ((*iter)->stateType==UserState::kEndOfSessionTemp)
                {   //delete temporary kEndOfSession
                    (*iter).remove();
                } else
                {   //do not delete, but just change this HEAD to kUserSaved state (no longer kEndOfSession
                    (*iter).modify()->stateType = UserState::kUserSaved;
                }
            }
            
            transaction.commit();
        }//if( cleanupStates )

        
      //Make sure there is at least a spectra valid in order to save this state
      std::shared_ptr<const SpecMeas> foreground = m_viewer->measurment( kForeground );
      std::shared_ptr<const SpecMeas> second = m_viewer->measurment( kSecondForeground );
      std::shared_ptr<const SpecMeas> background = m_viewer->measurment( kBackground );
        
      if( saveSpectra && saveState && (foreground || second || background))
      {
        const int offset = environment().timeZoneOffset();
        WString desc = "End of Session";
        WString name = WDateTime::currentDateTime().addSecs(60*offset).toString( DATE_TIME_FORMAT_STR );
        
        Wt::Dbo::ptr<UserState> dbstate;
        std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
        DataBaseUtils::DbTransaction transaction( *sql );
          
        UserState *state = new UserState();
        state->user = m_viewer->m_user;
        state->stateType = UserState::kEndOfSessionTemp;
        state->creationTime = WDateTime::currentDateTime();
        state->name = name;
        state->description = desc;
          
          //check if Save menu needs to be updated
          if( m_viewer->getCurrentStateID() >= 0 )
          {
              try
              {
                  std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
                  DataBaseUtils::DbTransaction transaction( *sql );
                  Dbo::ptr<UserState> currentState = sql->session()->find<UserState>( "where id = ?" ).bind( m_viewer->getCurrentStateID() );
                  transaction.commit();
                  
                  if (currentState)
                  {
                      //Saves to the HEAD first
                      m_viewer->startStoreStateInDb( false, false, false, true ); //save snapshot
                      transaction.commit();

                      return;
                  } //UserState exists
                  else
                  { //no UserState
                      //do nothing, a temporary kEndOfState HEAD snapshot will be created.
                  } //no UserState
              }catch( std::exception &e )
              {
                  cerr << "Could not access current Snapshot State ID while EndOfSession " << e.what() << endl;
              }//try / catch
          }

          
          
          
        dbstate = sql->session()->add( state );
        
        try
        {
          m_viewer->saveStateToDb( dbstate );
        }catch( std::exception &e )
        {
          cerr << "InterSpecApp::prepareForEndOfSession(): " << e.what() << endl;
        }
        transaction.commit();
      }//if( saveSpectra && saveState )
#endif //#if( USE_DB_TO_STORE_SPECTRA )
    }//if( we can accumulate usage stats )
    
  }catch( std::exception &e )
  {
    cerr << "InterSpecApp::prepareForEndOfSession() caught: " << e.what() << endl;
  }//try / catch
  
  cerr << "Have prepared for end of session " << sessionId() << endl;
}//void InterSpecApp::prepareForEndOfSession()


void InterSpecApp::clearSession()
{
  setupWidgets( false );
}//void clearSession()



void InterSpecApp::finalize()
{
  prepareForEndOfSession();
  
  if( m_viewer && m_viewer->m_user )
    cerr << "Have finalized session " << sessionId() << " for user "
         << m_viewer->m_user->userName() << endl;
}//void InterSpecApp::finalize()


#if(USE_SRB_HEADER_FOOTER)
void InterSpecApp::showHeaderFooter( bool setUserOption )
{
  if( setUserOption )
    InterSpecUser::setPreferenceValue( m_viewer->m_user, "ShowHeaderFooter", true, m_viewer );
  
  m_header->show();
  m_footer->show();

  m_showHeaderFooterButton->setHidden( true, WAnimation() );
  m_hideHeaderFooterButton->setHidden( false, WAnimation() );
}//void showHeaderFooter()


void InterSpecApp::hideHeaderFooter( bool setUserOption )
{
  if( setUserOption )
    InterSpecUser::setPreferenceValue( m_viewer->m_user, "ShowHeaderFooter", false, m_viewer );

  m_header->hide();
  m_footer->hide();

  // Toggle the visibility of the show/hide buttons
  m_showHeaderFooterButton->setHidden( false, WAnimation() );
  m_hideHeaderFooterButton->setHidden( true, WAnimation() );
}//void hideHeaderFooter()
#endif


void InterSpecApp::svlog( const WString& message, const WString& source, int priority )
{
  if( m_viewer )
    m_viewer->logMessage( message, source, priority );
}//void svlog(WString&, WString&)


bool InterSpecApp::isMobile() const
{
  const WEnvironment &env = environment();
  const bool isMob = (   env.agentIsMobileWebKit()
                      || env.agentIsIEMobile()
                      || env.agentIsMobileWebKit()
                      || env.agent()==WEnvironment::MobileWebKitiPhone
                      || env.agent()==WEnvironment::MobileWebKitAndroid
                      || env.agent()==WEnvironment::MobileWebKit
                      || env.userAgent().find("Opera Mobi") != std::string::npos
                      || env.userAgent().find("Android") != std::string::npos
                      || env.userAgent().find("RIM ") != std::string::npos
                      || env.userAgent().find("iPad") != std::string::npos
                      || env.userAgent().find("iPhone") != std::string::npos
                      || env.userAgent().find("iPod") != std::string::npos
                      || env.userAgent().find("Mobile") != std::string::npos
                      );
  return isMob;
}//bool isMobile() const


bool InterSpecApp::isPhone() const
{
  //see: (iOS) http://www.enterpriseios.com/wiki/UserAgent
  //     (Android) http://www.gtrifonov.com/2011/04/15/google-android-user-agent-strings-2/
  
  const WEnvironment &env = environment();
  return ( env.userAgent().find("iPhone") != std::string::npos
           || env.userAgent().find("iPod") != std::string::npos
           || (env.userAgent().find("Android") != std::string::npos
               && env.userAgent().find("Mobile") != std::string::npos)
           || env.userAgent().find("RIM ") != std::string::npos);
}//bool InterSpecApp::isPhone()

bool InterSpecApp::isTablet() const
{
    const WEnvironment &env = environment();
    const string &agent = env.userAgent();
    
    return (agent.find("iPad") != string::npos
            || (agent.find("Android") != string::npos
                && agent.find("Mobile") == string::npos)
            );
}//bool isTablet() const
