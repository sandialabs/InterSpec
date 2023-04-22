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
#include <Wt/WTimer>
#include <Wt/WLabel>
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

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"

#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#endif

#if( BUILD_AS_OSX_APP )
#include "target/osx/macOsUtils.h"
#endif

#if( BUILD_AS_WX_WIDGETS_APP )
#include "target/wxWidgets/InterSpecWxUtils.h"
#endif

#if( !BUILD_FOR_WEB_DEPLOYMENT )
#include "InterSpec/InterSpecServer.h"
#endif

#if( USE_REMOTE_RID )
#include "InterSpec/RemoteRid.h"
#endif


using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

#if( BUILD_AS_ELECTRON_APP )
WT_DECLARE_WT_MEMBER
(IsElectronInstance, Wt::JavaScriptFunction, "IsElectronInstance",
  function() {
   return (typeof window !== 'undefined' && typeof window.process === 'object' && window.process.type === 'renderer');
} );
#endif


namespace
{
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  std::mutex AppInstancesMutex;
  std::set<InterSpecApp *> AppInstances;
  //note: could potentially use Wt::WServer::instance()->sessions() to retrieve
  //      sessionIds.
#endif
}//namespace


InterSpecApp::InterSpecApp( const WEnvironment &env )
  :  WApplication( env ),
     m_viewer( 0 ),
     m_layout( nullptr ),
     m_lastAccessTime( std::chrono::steady_clock::now() ),
     m_activeTimeInSession{ std::chrono::seconds(0) }
#if( IOS )
    , m_orientation( InterSpecApp::DeviceOrientation::Unknown )
    , m_safeAreas{ 0.0f }
#endif
{
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( !checkExternalTokenFromUrl() )
  {
    setTitle( "Error loading" );
    new WText( "Invalid 'apptoken' specified, not loading", root() );
    WApplication::quit( "InterSpec is configured to not allow arbitrary sessions" );
    return;
  }//if( !checkExternalTokenFromUrl() )
#endif
 
  enableUpdates( true );
  
  // Lets get a jump on initializing DecayDataBaseServer.
  // Note however, in InterSpecServer::startServer and InterSpecServer::startWebServer we have
  //  already called DecayDataBaseServer::initialize() as soon as the thread pool was available, so
  //  this call will be a waste.
  //  TODO: once the localhost server is converted to using InterSpecServer, remove this next line
  WServer::instance()->ioService().boost::asio::io_service::post( &DecayDataBaseServer::initialize );
   
  setupDomEnvironment();
  setupWidgets( true );

  Wt::log("debug") << "Done in InterSpecApp constructor.";
}//InterSpecApp constructor


InterSpecApp::~InterSpecApp()
{
  Wt::log("debug") << "Entering ~InterSpecApp() destructor.";

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( !m_externalToken.empty() )
    InterSpecServer::set_session_destructing( m_externalToken.c_str() );
#endif

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  {
    std::lock_guard<std::mutex> lock( AppInstancesMutex );
    AppInstances.erase( this );
  }
#endif
}//~InterSpecApp()


#if( !BUILD_FOR_WEB_DEPLOYMENT )
bool InterSpecApp::checkExternalTokenFromUrl()
{
  m_primaryApp = false;
  
  const bool allow_untokened = InterSpecServer::allow_untokened_sessions();
  const Http::ParameterMap &parmap = environment().getParameterMap();
  for( const Http::ParameterMap::value_type &p : parmap )
  {
    if( (SpecUtils::iequals_ascii(p.first, "externalid")
          || SpecUtils::iequals_ascii(p.first, "apptoken"))
       && !p.second.empty() )
    {
      m_externalToken = p.second.front();
      
      const int status = InterSpecServer::set_session_loaded( m_externalToken.c_str() );
      
      pair<bool,InterSpecServer::SessionType> session_type
                                         = InterSpecServer::session_type( m_externalToken.c_str() );
      
      m_primaryApp = false;
      if( session_type.first )
      {
        switch( session_type.second )
        {
          case InterSpecServer::SessionType::PrimaryAppInstance:
            m_primaryApp = true;
            break;
            
          case InterSpecServer::SessionType::ExternalBrowserInstance:
            m_primaryApp = false;
            break;
        }//switch( session_type.second )
      }//if( session_type.first )
      
      return (allow_untokened || ((status==0) && session_type.first));
    }
  }//for( const Http::ParameterMap::value_type &p : parmap )
  
  return allow_untokened;
}//bool checkExternalTokenFromUrl()
#endif  //#if( !BUILD_FOR_WEB_DEPLOYMENT )


void InterSpecApp::setupDomEnvironment()
{
#if( BUILD_AS_ELECTRON_APP )
  LOAD_JAVASCRIPT(wApp, "InterSpecApp.cpp", "InterSpecApp", wtjsIsElectronInstance);
  
  //To make nodeIntegration=true work:
  //  https://stackoverflow.com/questions/32621988/electron-jquery-is-not-defined
  const char *fixElectronJs =
  "if( Wt.WT.IsElectronInstance() )"
    "if (typeof module === 'object') {window.module = module; module = undefined;}";
  doJavaScript( fixElectronJs, false );
  
  if( isPrimaryWindowInstance() )
  {
#if( USING_ELECTRON_NATIVE_MENU )
    //Some support to use native menu
    doJavaScript( "const {Menu, MenuItem} = remote;", false );
#endif
  }//if( isPrimaryWindowInstance() )
  
#endif //BUILD_AS_ELECTRON_APP

  setTitle( "InterSpec" );
  
  //Call tempDirectory() to set global variable holding path to the temp file
  //  directory, to ensure this will be available at all points in the future;
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
  
  setCssTheme( "default" );  //"polished" is the other option
  
  // Use newer version of jQuery than Wt uses by default
  requireJQuery("InterSpec_resources/assets/js/jquery-3.6.0.min.js");
  
  //for qTip2
  useStyleSheet( "InterSpec_resources/assets/js/qTip2-3.0.3/jquery.qtip.min.css" );
  require("InterSpec_resources/assets/js/qTip2-3.0.3/jquery.qtip.min.js");
  
  //qTip2 has a dependancy on the imagesloaded plugin, but from a quick and niave test, it doesnt seem to be compulsory.
  //Leaving next line commented out as a reminder, for the moment, that we may need it
  //require("InterSpec_resources/assets/js/imagesloaded.pkg.min.js");
  
  // For older browsers (primarily the WkWebView for macOS Mojave (10.14) and older).
  require("InterSpec_resources/assets/js/resize-observer-polyfill/ResizeObserver.js", "ResizeObserver");
  
  useStyleSheet( "InterSpec_resources/InterSpec.css" );
  doJavaScript( "if(typeof console==='undefined'){console={log:function(){}};}" );
  
  styleSheet().addRule( "input[type=\"text\"]", "font-size:0.95em;" );
  styleSheet().addRule( "button[type=\"button\"]", "font-size:0.9em;" );
  
  //This is needed to fix keyboard hiding in iOS
#if(IOS || ANDROID)
  //TODO: this can be added to wt_config.xml instead
  const char *fix_ios_js = INLINE_JAVASCRIPT(
    var t=document.createElement('meta');
    t.name = "viewport";
    t.content = "initial-scale=1.0, maximum-scale=1.0, user-scalable=no, height=device-height, width=device-width, viewport-fit=cover";
    document.getElementsByTagName('head')[0].appendChild(t);
    $(document).on('blur', 'input, textarea', function() {
      setTimeout(function(){ window.scrollTo(document.body.scrollLeft, document.body.scrollTop); }, 0);
    });
  );
  
  doJavaScript( fix_ios_js );
#endif
  
#if( BUILD_AS_OSX_APP && !PERFORM_DEVELOPER_CHECKS )
  domRoot()->setAttributeValue( "oncontextmenu", "return false;" );
#endif
  
  // Define some javascript to artificially trigger a resize event; this is a hack used a few
  //  places to invoke the Wt layout stuff, when its needed, but not being triggered automatically.
  //  To call this code, use something like:
  //    wApp->doJavaScript(wApp->javaScriptClass() + ".TriggerResizeEvent();");
  declareJavaScriptFunction( "TriggerResizeEvent",
  "function(){"
    "const a = function(ms){"
      "setTimeout( function(){ "
        + wApp->javaScriptClass() + ".layouts2.scheduleAdjust();"
        "window.dispatchEvent(new Event('resize')); "
      "}, ms );"
    "};"
    // The number of calls, or when the calls are made to trigger resizes has not been investigated
    "a(0); a(50); a(500); a(2500);"
  "}"
  );
  

#if(  BUILD_AS_WX_WIDGETS_APP || (BUILD_AS_ELECTRON_APP && !USING_ELECTRON_NATIVE_MENU) )
  // A workaround to allow moving the app window around by the titlebar area, when a 
  //  AuxWindow or SimpleDialog is showing
  if (isPrimaryWindowInstance())
  {
#if(  BUILD_AS_WX_WIDGETS_APP )
    declareJavaScriptFunction("MouseDownOnDialogCover",
      "function(el, evt) {"
        "if( (evt.button!==0) || (el.id !== evt.target.id) || evt.pageY > 30) return;"
        "if(window.wx) window.wx.postMessage('MouseDownInTitleBar');"
        "evt.stopPropagation();"
        "evt.preventDefault();"
      "}"
    );
#endif

    // Setup a function to raise windows controls (minimize, maximize, close), 
    // to above the dialog-cover, so we dont block people from being able to 
    // close the app.  I dont think the setTimeout(...) is needed, but JIC
    //usage: wApp->doJavaScript(wApp->javaScriptClass() + ".RaiseWinCntrlsAboveCover();");
    declareJavaScriptFunction("RaiseWinCntrlsAboveCover",
      "function(){ "
        "setTimeout(function(){"
          "const z = $('.Wt-dialogcover').css('z-index');"
          "if(z) $('.window-controls-container').css('z-index',z+1);"
        "},250);"
      "}"
    );
  }//if (isPrimaryWindowInstance())
#endif
  
  //cout << "wApp->javaScriptClass()='" << wApp->javaScriptClass() << "'" << endl; //Prints "Wt"
  
  
  if( isMobile() )
  {
    // if isTablet(), then these css files may later get unloaded in setupWidgets().
    useStyleSheet( "InterSpec_resources/InterSpecSafeArea.css" );
    useStyleSheet( "InterSpec_resources/InterSpecMobileCommon.css" );
    
    if( isAndroid() )
      useStyleSheet( "InterSpec_resources/InterSpecMobileDroid.css" );
    else
      useStyleSheet( "InterSpec_resources/InterSpecMobileApple.css" );
    
    //Prevent mobile from showing spinner
    const char *prevent_spinner_js = INLINE_JAVASCRIPT(
      $("<style>").prop("type", "text/css")
                  .html(".Wt-spinbox {background-image: inherit !important; \
                         background-repeat: inherit  !important; \
                         background-position: inherit  !important; \
                         padding-right: inherit  !important; }")
                  .appendTo("head");
    );
    
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
  
  // Pre-load some CSS we will likely encounter anyway; avoids some glitching when loading new
  //  widgets.
  //  (The CSS wont be re-loaded later, so maybe we dont hurt anything doing it here too - maybe)
  WServer::instance()->schedule( 100, sessionId(), [](){
    wApp->useStyleSheet( "InterSpec_resources/SimpleDialog.css" );
    wApp->useStyleSheet( "InterSpec_resources/DrfSelect.css" );
    
    // anything else relevant?
    wApp->triggerUpdate();
  } );
  
  
#if( IOS )
  //Check if device orientation and SafeAreas are specified in the URL.
  const Http::ParameterMap &parmap = environment().getParameterMap();
  const Http::ParameterMap::const_iterator iter = parmap.find( "SafeAreas" );
  if( iter != parmap.end() && iter->second.size() )
  {
    //Try to parse out an integer followed by four floats
    int orientation;
    float top, right, bottom, left;
    const int nargscanned = sscanf( iter->second[0].c_str(), "%i,%f,%f,%f,%f",
                              &orientation, &top, &right, &bottom, &left );
    
    if( nargscanned == 5 )
    {
      const DeviceOrientation orient = DeviceOrientation(orientation);
      Wt::log("debug") << "Received initial orientation: " << orientation
           << " and SafeAreas={" << top << ", " << right << ", " << bottom << ", " << left << "}";
      if( ((orient != DeviceOrientation::LandscapeLeft)
           && (orient != DeviceOrientation::LandscapeRight))
         || (top < 0)    || (top > 50)
         || (right < 0)  || (right > 75)
         || (bottom < 0) || (bottom > 50)
         || (left < 0)   || (left > 75) )
      {
        Wt::log("error") << "Initial device orientations invalid!";
      }else
      {
        setSafeAreaInsets( orientation, top, right, bottom, left );
      }
    }else
    {
      Wt::log("error") << "Failed to parse SafeAreas='" << iter->second[0] << "'";
    }
  }
#endif

}//void setupDomEnvironment()


void InterSpecApp::setupWidgets( const bool attemptStateLoad  )
{
  if( m_viewer )
  {
    //for some reason calling 'delete m_layout' doesnt cause m_viewer to
    //  destruct?
    delete m_viewer;
    m_viewer = nullptr;
    
#if( USING_ELECTRON_NATIVE_MENU )
    //ToDo: need to clear out all the native menus...
#elif( USE_OSX_NATIVE_MENU )
#endif
  }
  
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
    cerr << "There was an exception initializing a InterSpec: " << e.what() << endl;
    throw e;
#endif //#if( BUILD_AS_UNIT_TEST_SUITE )
    
    string msg = "There was a problem initializing a necessary resource for InterSpec";
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
    
    msg = "<br />Please contact wcjohns@sandia.gov and/or interspec@sandia.gov to fix this error.";
    errorText = new WText( msg, root() );
    errorText->setAttributeValue( "style", "font-family:Courier New; color: blue;" );
    WApplication::quit();
    return;
  }//try / catch to create a InterSpec object
  
  root()->addStyleClass( "specviewer" );
  
  // TODO: we could add an explicit CSS class 
  // @media screen and (max-device-width: 640px) { ... }
  //if( isPhone() )
  //  root()->addStyleClass( "is-phone" );  //see also LandscapeRight and LandscapeLeft CSS classes
  //else if( isTablet() )
  //  root()->addStyleClass( "is-tablet" );
  //if( isMobile() )
  //  root()->addStyleClass( "is-mobile" );
  
  if( !m_miscSignal )
  {
    m_miscSignal.reset( new JSignal<std::string>(root(), "miscSignal", false) );
    m_miscSignal->connect( this, &InterSpecApp::miscSignalHandler );
  }
  
  m_layout = new WGridLayout();
  root()->setLayout( m_layout );
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setHorizontalSpacing( 0 );
  m_layout->setVerticalSpacing( 0 );
  
  
  m_layout->addWidget( m_viewer, 0, 0, 1, 1 );
  m_layout->setRowStretch( 0, 10 );
  
  // TODO: It seems sometimes when closing a AuxWindow that is modal, the grey background doesnt actually get hidden... not totally sure why, but we *could* maybe help this with the below - but I havent tested it doesnt have odd side-effects, so not currently enabling, but should consider
  //globalEscapePressed().connect( "function(){$('.Wt-dialogcover').hide();}" );
  
  
  bool loadedSpecFile = false;
  const Http::ParameterMap &parmap = environment().getParameterMap();
  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( !m_externalToken.empty() )
  {
    const string initial_file = InterSpecServer::file_to_open_on_load( m_externalToken );
    if( !initial_file.empty() )
    {
      
      if( SpecUtils::istarts_with(initial_file, "interspec://")
         || SpecUtils::istarts_with(initial_file, "RADDATA://G0/") )
      {
        try
        {
          loadedSpecFile = handleAppUrl( initial_file );
        }catch( std::exception &e )
        {
          passMessage( "Error handling deep-link (1): " + string(e.what()),
                      WarningWidget::WarningMsgHigh );
          wApp->log( "error" ) << "InterSpecApp::setupWidgets: invalid URL: " << e.what();
        }
      }else
      {
        loadedSpecFile = m_viewer->userOpenFileFromFilesystem( initial_file );
      }//if( SpecUtils::istarts_with(foregroundPath, "interspec://") )
     
      InterSpecServer::clear_file_to_open_on_load( m_externalToken );
      
      if( loadedSpecFile )
      {
        Wt::log("debug") << "Opened file/url from initial request: '" << initial_file << "'";
      }else
      {
        //userOpenFileFromFilesystem will call SpecMeasManager::displayInvalidFileMsg for us
        Wt::log("error") << "Invalid file/url from initial request: '" << initial_file << "'";
      }
    }//if( !initial_file.empty() )
  }//if( !m_externalToken.empty() )
#endif
  

#if( USE_DB_TO_STORE_SPECTRA )
  //Check to see if we should load the apps last saved state
  try
  {
    const bool saveSpectra =  true; // InterSpecUser::preferenceValue<bool>("SaveSpectraToDb", m_viewer );
    bool saveState = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
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
        Wt::log("error") << "Failed to commit transaction while loading state, caught: "
        << e.what();
      }//try / catch
      
      if( state )
      {
        loadedSpecFile = true;
#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
        if( promptLoad )
        {
          //Create a dialog asking if the user wants to pick up
          AuxWindow *loadStateDialog = new AuxWindow( "Load previous state?", (AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen) );
        
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
              Wt::log("error") << "Failed setting 'PromptStateLoadOnStart' user prefrence/";
            }
          };
        
          auto changeDoLoadPref = [this](bool load ){
            try {
              InterSpecUser::setPreferenceValue<bool>( m_viewer->m_user,
                                    "LoadPrevStateOnStart", load, m_viewer );
            }catch( std::exception &e ) {
              Wt::log("error") << "Failed setting 'LoadPrevStateOnStart' user prefrence.";
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
        js << "Resuming where you left off on " << state->name.toUTF8()
           << "<div onclick="
        "\"Wt.emit('" << root()->id() << "',{name:'miscSignal'}, 'clearSession');"
        //"$('.qtip.jgrowl:visible:last').remove();"
        "try{$(this.parentElement.parentElement).remove();}catch(e){}"
        "return false;\" "
        "class=\"clearsession\"><span class=\"clearsessiontxt\">Start Fresh Session</span></div>";
        
        m_viewer->warningWidget()->addMessageUnsafe( js.str(), WarningWidget::WarningMsgInfo, 5000 );
#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
        }
#endif
      }else
      {
        Wt::log("error") << "Could not load previously saved state.";
      }
    }//if( saveSpectra && saveState )
  }catch( std::exception &e )
  {
    Wt::log("error") << "Failed to load app state, caught: '" << e.what() << "'.";
  }//try / catch
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  
#if( USE_DB_TO_STORE_SPECTRA && INCLUDE_ANALYSIS_TEST_SUITE )
  bool userstate = true;
  Http::ParameterMap::const_iterator stateiter = parmap.end();
  
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
    Wt::log("info") << "SpecViewerApp: Will try to open test state ID="
    << id << " (" << id << ").";
    
    if( state )
    {
      m_viewer->loadStateFromDb( state );
      loadedSpecFile = true;
    }else
    {
      Wt::log("error") << "Could not load test state.";
    }//if( state ) / else
  }//if( (stateiter != parmap.end()) && stateiter->second.size() )
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )

  
  auto showWelcomeCallback = [this,loadedSpecFile](){
    //If we already loaded a spectrum, then dont show welcome dialog (or IE
    //  warning dialog)
    if (loadedSpecFile)
    {
      return;
    }

    //If client is internet explorer, show a warning before the welcome dialog
    if( !environment().agentIsIE() )
    {
      //Using WTimer as a workaround for iOS so screen size and safe-area and
      //  such can all get setup before creating a AuxWindow; otherwise size of
      //  window will be all messed up.
      WTimer::singleShot( 10, boost::bind( &InterSpec::showWelcomeDialog, m_viewer, false) );
      
#if( IOS || ANDROID )
      //For iPhoneX* devices we should trigger a resize once
      doJavaScript( javaScriptClass() + ".TriggerResizeEvent();" );
#endif
      
    }else
    {
      AuxWindow *dialog = m_viewer->showIEWarningDialog();
      if( dialog )
        dialog->finished().connect( boost::bind( &InterSpec::showWelcomeDialog, m_viewer, false ) );
      else
        m_viewer->showWelcomeDialog();
    }// if( not in IE ) / else
  };//auto showWelcomeCallback
  
  if( !loadedSpecFile )
  {
    showWelcomeCallback();
  }//if( specfile.empty() )
  
  //Only display warning message on navigation away from current session if in
  //  a browser environment, and not this isnt a developers build
#if( (BUILD_FOR_WEB_DEPLOYMENT || BUILD_AS_LOCAL_SERVER) && !PERFORM_DEVELOPER_CHECKS  )
  
#if( !OPTIMISTICALLY_SAVE_USER_STATE )
  setConfirmCloseMessage( "This action will terminate your current session." );
#else
  //The below solution of catching onbeforeunload catches navigating to a new
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
  
#if( BUILD_AS_OSX_APP || IOS || BUILD_AS_ELECTRON_APP  || BUILD_AS_WX_WIDGETS_APP  )
  auto themeiter = parmap.find( "colortheme" );
  if( themeiter != parmap.end() && themeiter->second.size() )
    m_viewer->osThemeChange( themeiter->second[0] );
#endif

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  {
    std::lock_guard<std::mutex> lock( AppInstancesMutex );
    AppInstances.insert( this );
  }
#endif
  
  if( m_viewer->m_user )
    Wt::log("debug") << "Have started session " << sessionId() << " for user "
    << m_viewer->m_user->userName() << ".";
  
#if( USING_ELECTRON_NATIVE_MENU )
  if( isPrimaryWindowInstance() )
  {
    //ToDo: time if WTimer::singleShot or post() is better
    //WServer::instance()->post( sessionId(), [](){ PopupDivMenu::triggerElectronMenuUpdate(); } );
    PopupDivMenu::triggerElectronMenuUpdate();
  }//if( !m_externalToken.empty() )
#endif
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || ANDROID || BUILD_AS_WX_WIDGETS_APP  )
  // TODO 20220405: in macOS app I see the error "Wt: decodeSignal(): signal 'of7g69f.SucessfullyLoadedConf' not exposed"
  //                Not sure how to fix this, atm
  m_sucessfullyLoadedSignal.reset( new Wt::JSignal<>( m_viewer, "SucessfullyLoadedConf", false ) );
  m_sucessfullyLoadedSignal->connect( this, &InterSpecApp::loadSuccesfullCallback );
  doJavaScript( "setTimeout(function(){" + m_sucessfullyLoadedSignal->createCall() + "}, 250);" );
#endif
  
  // The user can override the mobile setting for tablets, so we'll fix-up what CSS we load here.
  //  I think we could delay loading the CSS for mobile from setupDomEnvironment() to here, without
  //  any effect, but I havent tested yet - and I dont feel totally great about how things are
  //  structured w.r.t. the "TabletUseDesktopMenus" user preference yet.
  if( isTablet() )
  {
    // The user could have changed the "TabletUseDesktopMenus" preference, and done a
    //  "Clear Session...".
    //  InterSpecApp::isTablet() does not account for the user preference, while
    //  InterSpec::isTablet() does.
    //  Note that we are not unloading InterSpecSafeArea.css, even if user has requested "Desktop"
    //   version of app
    if( m_viewer->isTablet() )
    {
      // WApplication::useStyleSheet will essentially be a no-op if the file has already been loaded
      useStyleSheet( "InterSpec_resources/InterSpecMobileCommon.css" );
      if( isAndroid() )
        useStyleSheet( "InterSpec_resources/InterSpecMobileDroid.css" );
      else
        useStyleSheet( "InterSpec_resources/InterSpecMobileApple.css" );
    }else
    {
      // WApplication::useStyleSheet will essentially be a no-op if the file has not been loaded
      removeStyleSheet( "InterSpec_resources/InterSpecMobileCommon.css" );
      if( isAndroid() )
        removeStyleSheet( "InterSpec_resources/InterSpecMobileDroid.css" );
      else
        removeStyleSheet( "InterSpec_resources/InterSpecMobileApple.css" );
    }//
  }//if( isTablet() )
}//void setupWidgets()



std::string InterSpecApp::getUserNameFromEnvironment() const
{
  string remoteUser = environment().getCgiValue( "REMOTE_USER" );
  SpecUtils::ireplace_all( remoteUser, " ", "" );
  SpecUtils::to_lower_ascii( remoteUser );
  return remoteUser;
}//std::string getUserNameFromEnvironment() const


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
    if( s_tempDir.empty() || !SpecUtils::is_directory(s_tempDir) )
      s_tempDir = SpecUtils::temp_dir();
    
    return s_tempDir;
  }//if( app )
  
  return SpecUtils::temp_dir();
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


InterSpec *InterSpecApp::viewer()
{
  return m_viewer;
}//InterSpec* viewer()


std::chrono::steady_clock::time_point::duration InterSpecApp::activeTimeInCurrentSession() const
{
  return m_activeTimeInSession;
}

#if( !BUILD_FOR_WEB_DEPLOYMENT )
std::string InterSpecApp::externalToken()
{
  WApplication::UpdateLock lock( this );
  return m_externalToken;
}//const std::string &externalToken() const


bool InterSpecApp::isPrimaryWindowInstance()
{
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( !app )
    return false;
  
  WApplication::UpdateLock lock(app);
  if( !lock )
    return false;
  
  if( !app->m_primaryApp )
    return false;
  
  const string &externalid = app->m_externalToken;
  const int status = InterSpecServer::session_status( externalid.c_str() );
  
  const pair<bool,InterSpecServer::SessionType> session_type = InterSpecServer::session_type( externalid.c_str() );
  
  //app->m_primaryApp should only be true if its SessionType::PrimaryAppInstance, but we'll double
  //  check this, JIC
  assert( !session_type.first || (session_type.second == InterSpecServer::SessionType::PrimaryAppInstance) );
  
  return (((status == 1) || (status == 2))
           && (session_type.first && session_type.second == InterSpecServer::SessionType::PrimaryAppInstance));
}//bool isPrimaryWindowInstance()


InterSpecApp *InterSpecApp::instanceFromExtenalToken( const std::string &idstr )
{
  if( idstr.empty() )
    return (InterSpecApp *)0;
  
  std::lock_guard<std::mutex> lock( AppInstancesMutex );
  //Wt::log("debug") << "There are AppInstances=" << AppInstances.size() << "; we want session: '" << idstr << "'";
  for( InterSpecApp *app : AppInstances )
  {
    Wt::WApplication::UpdateLock lock( app );
    //Wt::log("debug") << "\tSession '" << idstr << "' An instance=" << app->externalToken();
    if( app->externalToken() == idstr )
      return app;
  }
  
  return nullptr;
}//InterSpecApp *instanceFromExtenalToken( const std::string &idstr )


bool InterSpecApp::userOpenFromFileSystem( const std::string &path )
{
  if( !m_viewer )
    return false;
  return m_viewer->userOpenFileFromFilesystem( path );
}//bool userOpenFromFileSystem(...)


std::set<InterSpecApp *> InterSpecApp::runningInstances()
{
  std::lock_guard<std::mutex> lock( AppInstancesMutex );
  return AppInstances;
}//set<InterSpecApp *> runningInstances()

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



#if( BUILD_AS_OSX_APP || IOS || BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
void InterSpecApp::osThemeChange( std::string name )
{
  auto server = WServer::instance();
  if( !server )
    return;
  
  server->postAll( [name](){
    auto app = dynamic_cast<InterSpecApp *>(wApp);
    if( !app )
      return;  //probably shouldnt ever happen
    app->m_viewer->osThemeChange(name);
    app->triggerUpdate();
  } );
}//osThemeChange(...)
#endif


#if( IOS )
void InterSpecApp::setSafeAreaInsets( const int orientation, const float top,
                       const float right, const float bottom,
                       const float left )
{
  m_orientation = static_cast<DeviceOrientation>( orientation );
  m_safeAreas[0] = top;
  m_safeAreas[1] = right;
  m_safeAreas[2] = bottom;
  m_safeAreas[3] = left;
  
  //ToDo: see if triggering a resize event is ever necessary
  //  doJavaScript( javaScriptClass() + ".TriggerResizeEvent();" );
  
  // Note that CSS takes care of insets, mostly by detecting the LandscapeLeft and LandscapeRight
  //  CSS classes, which are set by the `DoOrientationChange` javascript function
  
  Wt::log("debug") << "Set safe area insets: orientation=" << orientation
       << ", safeAreas={" << top << ", " << right << ", "
       << bottom << ", " << left << "}";
}//setSafeAreaInsets(...)

void InterSpecApp::getSafeAreaInsets( InterSpecApp::DeviceOrientation &orientation,
                                      float &top, float &right, float &bottom, float &left )
{
  orientation = m_orientation;
  top = m_safeAreas[0];
  right = m_safeAreas[1];
  bottom = m_safeAreas[2];
  left = m_safeAreas[3];
}//getSafeAreaInsets(...)
#endif

void InterSpecApp::notify( const Wt::WEvent& event )
{
#if( !BUILD_AS_UNIT_TEST_SUITE )
  try
  {
#endif
    const bool userEvent = (event.eventType() == Wt::UserEvent);
    if( userEvent )
    {
      const auto thistime = std::chrono::steady_clock::now();
      const auto duration = thistime - m_lastAccessTime;
      if( duration < std::chrono::seconds(300) )
        m_activeTimeInSession += duration;
      m_lastAccessTime = thistime;
    }//if( userEvent )

     WApplication::notify( event );
#if( !BUILD_AS_UNIT_TEST_SUITE )
  }catch( std::exception &e )
  {
    string msg = "There was an unexpected error, application state may be suspect: ";
    msg += e.what();
    svlog( msg, WarningWidget::WarningMsgHigh );
    
    char message[512];
    snprintf( message, sizeof(message), "Uncaught exception in event loop: '%s'.", e.what() );
    Wt::log("error") << message;
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, message );
#endif
  }//try/catch
#endif //#if( !BUILD_AS_UNIT_TEST_SUITE )
}//void notify( const Wt::WEvent& event )


void InterSpecApp::unload()
{
  Wt::log("info") << "Received unload from " << sessionId() << ".";
  WApplication::unload();
}//void unload()


void InterSpecApp::prepareForEndOfSession()
{
  if( !m_viewer || !m_viewer->m_user )
    return;

  try
  {
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
      
      {
        DataBaseUtils::DbTransaction transaction( *sql );
        
        //If the user has multiple sessions going on, the below will quash
        //  each other out.  Could probably be fixed by
        //  m_viewer->m_user.reread();
        
        m_viewer->m_user.modify()->addUsageTimeDuration( m_activeTimeInSession );
        
        const chrono::milliseconds nmilli = chrono::duration_cast<chrono::milliseconds>(m_activeTimeInSession);
        Wt::log("info") << "At session end, added " << nmilli.count() << " ms usage time for user "
             << m_viewer->m_user->userName();
        
        m_activeTimeInSession = std::chrono::seconds(0);
        transaction.commit();
      }
      
#if( USE_DB_TO_STORE_SPECTRA )
      //Check to see if we should save the apps state
      const bool saveState = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
        
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
                { 
                  //delete temporary kEndOfSession
                  iter->remove();
                } else
                {   
                  //do not delete, but just change this HEAD to kUserSaved state (no longer kEndOfSession
                   iter->modify()->stateType = UserState::kUserSaved;
                }
            }
            
            transaction.commit();
        }//if( cleanupStates )

        
      //Make sure there is at least a spectra valid in order to save this state
      std::shared_ptr<const SpecMeas> foreground = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
        
      if( saveState && foreground )
      {
        Wt::log( "info" ) << "Will auto-save state for session-id:" << sessionId() << " at end of session.";

        const int offset = environment().timeZoneOffset();
        WString desc = "End of Session";
        const auto now = chrono::system_clock::now() + chrono::seconds( 60*offset );
        WString name = SpecUtils::to_common_string( chrono::time_point_cast<chrono::microseconds>(now), true ); //"9-Sep-2014 15:02:15"
        
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
                Wt::log("error") << "Could not access current Snapshot State ID while EndOfSession " << e.what();
              }//try / catch
          }
          
        dbstate = sql->session()->add( state );
        
        try
        {
          m_viewer->saveStateToDb( dbstate );
          Wt::log("debug") << "Saved end of session app state to database";
        }catch( std::exception &e )
        {
          Wt::log("error") << "InterSpecApp::prepareForEndOfSession() error: " << e.what();
        }
        transaction.commit();
      }//if( saveState )
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  }catch( std::exception &e )
  {
    Wt::log("error") << "InterSpecApp::prepareForEndOfSession() caught: " << e.what();
  }//try / catch
  
  Wt::log("debug") << "Have prepared for end of session " << sessionId() << ".";
}//void InterSpecApp::prepareForEndOfSession()

#if(  BUILD_AS_WX_WIDGETS_APP )
void InterSpecApp::handleJavaScriptError( const std::string &errorText )
{
  // It doesnt look like we can call wxWidgets here via JS, so we will call into wxWidgets event loop
  doJavaScript( "console.log('Here I am after error');", false );
  
  if( isPrimaryWindowInstance() )
    InterSpecWxUtils::handle_javascript_error( errorText, m_externalToken );
  
  // Default WApplication implementation just logs error, and then calls WApplication:quit()
  WApplication::handleJavaScriptError( errorText );
}//void handleJavaScriptError( const std::string &errorText )
#endif

void InterSpecApp::clearSession()
{
  // Just in case there are any modal dialogs showing - the cover thing can be a little sticky
  doJavaScript( "$('.Wt-dialogcover').hide();" );
  
#if( BUILD_AS_ELECTRON_APP )
  // As a workaround setup a function ElectronUtils::requestNewCleanSession() that
  //  sends websocket message to node land to clear menus, and load a new session
  //  but with argument "restore=no"
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    if( !ElectronUtils::requestNewCleanSession() )
    {
      passMessage( "There was an error requesting a new session - sorry.",
                   WarningWidget::WarningMsgHigh );
    }
  }else
  {
    setupWidgets( false );
  }
#else
  setupWidgets( false );
#endif
}//void clearSession()


void InterSpecApp::miscSignalHandler( const std::string &signal )
{
  if( signal == "clearSession" )
  {
    clearSession();
    return;
  }//if( signal == "clearSession" )
  
  if( SpecUtils::istarts_with( signal, "showRiidAna" ) )
  {
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    if( SpecUtils::icontains(signal, "background") )  //SpecUtils::descriptionText(SpecUtils::SpectrumType::Background)
      type = SpecUtils::SpectrumType::Background;
    else if( SpecUtils::icontains(signal, "secondary") )  //SpecUtils::descriptionText(SpecUtils::SpectrumType::SecondForeground)
      type = SpecUtils::SpectrumType::SecondForeground;
    
    showRiidInstrumentsAna( m_viewer->measurment(type) );
    return;
  }
  
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  if( SpecUtils::istarts_with( signal, "showMap-" ) )
  {
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    if( SpecUtils::icontains(signal, "background") )
      type = SpecUtils::SpectrumType::Background;
    else if( SpecUtils::icontains(signal, "secondary") )
      type = SpecUtils::SpectrumType::SecondForeground;
    
    m_viewer->createMapWindow( type );
    return;
  }//if( SpecUtils::istarts_with( signal, "showMap-" ) )
#endif //#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
 
  if( SpecUtils::istarts_with( signal, "showMsg-" ) )
  {
    string msg = signal.substr(8);
    WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel::WarningMsgInfo;
  
    if( SpecUtils::istarts_with( msg, "info-" ) )
    {
      msg = msg.substr(5);
      level = WarningWidget::WarningMsgLevel::WarningMsgInfo;
    }else if( SpecUtils::istarts_with( msg, "error-" ) )
    {
      msg = signal.substr(6);
      level = WarningWidget::WarningMsgLevel::WarningMsgHigh;
    }
    
    passMessage( msg, level );
    return;
  }//if( SpecUtils::istarts_with( signal, "showMsg-" ) )
  
#if( USE_REMOTE_RID )
  if( signal == "openRemoteRidTool" )
  {
    RemoteRid::handleOpeningRemoteRidTool(m_viewer);
    return;
  }//if( signal == "openRemoteRidTool" )
  
  if( signal == "disableRemoteRid" )
  {
    RemoteRid::disableAutoRemoteRid(m_viewer);
    return;
  }//if( signal == "disableRemoteRid" )
  
  
  if( SpecUtils::istarts_with( signal, "showRemoteRidRefLines-" ) )
  {
    RemoteRid::handleShowingRemoteRidRefLines(m_viewer,signal);
    return;
  }//if( signal == "showRemoteRidRefLines-" )
#endif
  
  if( SpecUtils::istarts_with( signal, "showMultimedia-" ) )
  {
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    if( SpecUtils::icontains(signal, "background") )  //SpecUtils::descriptionText(SpecUtils::SpectrumType::Background)
      type = SpecUtils::SpectrumType::Background;
    else if( SpecUtils::icontains(signal, "secondary") )  //SpecUtils::descriptionText(SpecUtils::SpectrumType::SecondForeground)
      type = SpecUtils::SpectrumType::SecondForeground;
    
    m_viewer->showMultimedia( type );
    
    return;
  }//if( SpecUtils::istarts_with( signal, "showMultimedia" ) )


  // shouldnt ever make it here..
  const string errmsg = "InterSpecApp::miscSignalHandler: unhandled signal '" + signal + "'";
  passMessage( errmsg, WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
  log_developer_error( __func__, errmsg.c_str() );
#endif
  Wt::log("error") << errmsg;
}//void miscSignalHandler( const std::string &signal )


bool InterSpecApp::handleAppUrl( const std::string &url )
{
  try
  {
    m_viewer->handleAppUrl( url );
  }catch( std::exception &e )
  {
    passMessage( "Error handling deep-link: " + string(e.what()), WarningWidget::WarningMsgHigh );
    cerr << "InterSpecApp::handleAppUrl: invalid URL: " << e.what() << endl;
    wApp->log( "error" ) << "InterSpecApp::handleAppUrl: invalid URL: " << e.what();
    
    return false;
  }// try / catch to process th URL
  
  return true;
}//bool handleAppUrl( std::string url )



void InterSpecApp::finalize()
{
  prepareForEndOfSession();
  
  if( m_viewer && m_viewer->m_user )
    Wt::log("info") << "Have finalized session " << sessionId() << " for user "
         << m_viewer->m_user->userName();
}//void InterSpecApp::finalize()


void InterSpecApp::svlog( const WString& message, int priority )
{
  if( m_viewer )
    m_viewer->logMessage( message, priority );
}//void svlog(WString&, int)


void InterSpecApp::svlog( const char *message, int priority )
{
  svlog( Wt::WString::fromUTF8(message), priority );
}

void InterSpecApp::svlog( const std::string& message, int priority )
{
  svlog( Wt::WString::fromUTF8(message), priority );
}


bool InterSpecApp::isMobile() const
{
#if( IOS || ANDROID )
  return true;
#endif

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

bool InterSpecApp::isAndroid() const
{
#if( ANDROID )
  return true;
#elif( BUILD_AS_ELECTRON_APP || IOS || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
  return false;
#endif
  
  const WEnvironment &env = environment();
  const bool isDroid = ( env.agent()==WEnvironment::MobileWebKitAndroid
                      || env.userAgent().find("Android") != std::string::npos
                      || env.userAgent().find("RIM ") != std::string::npos
                      );
  return isDroid;
}


bool InterSpecApp::isPhone() const
{
  const WEnvironment &env = environment();
  
#if( IOS )
  // TODO: we could enable this for all builds, to help with testing; if we do, add this equiv code to isMobile()
  const Http::ParameterMap &parmap = environment().getParameterMap();
  const Http::ParameterMap::const_iterator iter = parmap.find( "isphone" );
  if( (iter != parmap.end()) && iter->second.size() && (iter->second[0] == "1") )
    return true;
#endif
  
  return ( env.userAgent().find("iPhone") != std::string::npos
           || env.userAgent().find("iPod") != std::string::npos
           || (env.userAgent().find("Android") != std::string::npos
               && env.userAgent().find("Mobile") != std::string::npos)
           || env.userAgent().find("RIM ") != std::string::npos);
}//bool InterSpecApp::isPhone()


bool InterSpecApp::isTablet() const
{
  const WEnvironment &env = environment();
  
#if( IOS )
  // TODO: we could enable this for all builds, to help with testing; if we do, add this equiv code to isMobile()
  const Http::ParameterMap &parmap = environment().getParameterMap();
  const Http::ParameterMap::const_iterator iter = parmap.find( "istablet" );
  if( (iter != parmap.end()) && iter->second.size() && (iter->second[0] == "1") )
    return true;
#endif
  
  const string &agent = env.userAgent();
    
  return (agent.find("iPad") != string::npos
            || (agent.find("Android") != string::npos
                && agent.find("Mobile") == string::npos)
            );
}//bool isTablet() const


#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || ANDROID  || BUILD_AS_WX_WIDGETS_APP )
void InterSpecApp::loadSuccesfullCallback()
{
  Wt::log("debug") << "Succesfully loaded session.";

  m_sucessfullyLoadedSignal.reset();
  
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
#if( BUILD_AS_ELECTRON_APP )
    ElectronUtils::notifyNodeJsOfNewSessionLoad();
    ElectronUtils::send_nodejs_message( "SessionFinishedLoading", "" );
#elif( BUILD_AS_OSX_APP )
    macOsUtils::sessionSuccessfullyLoaded();
#elif( ANDROID )
#warning "Need to implement notifying for parent process for Android"
#elif( BUILD_AS_WX_WIDGETS_APP )
    doJavaScript("window.wx.postMessage('SessionFinishedLoading');");
#else
    static_assert( 0, "Something messed up with pre-processor setup" );
#endif
  }//if( InterSpecApp::isPrimaryWindowInstance() )
}//void loadSuccesfullCallback()
#endif
