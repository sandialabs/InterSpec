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

// Get rid of some issues lurking in the boost libraries.
#pragma warning(disable:4244)
#pragma warning(disable:4800)
// Block out some UUID warnings
#pragma warning(disable:4996)

#include <ctime>
#include <tuple>
#include <mutex>
#include <locale>
#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

#include <boost/ref.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WBreak>
#include <Wt/WPoint>
#include <Wt/WServer>
#include <Wt/Dbo/Dbo>
#include <Wt/WSpinBox>
#include <Wt/WTextArea>
#include <Wt/WCheckBox>
#include <Wt/WTemplate>
#include <Wt/WIOService>
#include <Wt/WAnimation>
#include <Wt/WTabWidget>
#include <Wt/WPopupMenu>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WSplitButton>
#include <Wt/WBorderLayout>
#include <Wt/WSelectionBox>
#include <Wt/WCssStyleSheet>
#include <Wt/WSuggestionPopup>
#include <Wt/WContainerWidget>
#include <Wt/WDefaultLoadingIndicator>

#if( USE_CSS_FLEX_LAYOUT )
#include <Wt/WStackedWidget>
#endif

#if( USE_DB_TO_STORE_SPECTRA )
#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>
#include <Wt/Json/Serializer>
#endif

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "InterSpec/MakeDrf.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/FluxTool.h"
#include "InterSpec/PeakEdit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/HelpSystem.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/OneOverR2Calc.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DoseCalcWidget.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/SpecFileSummary.h"
#include "InterSpec/ColorThemeWindow.h"
#include "InterSpec/GammaCountDialog.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/CompactFileManager.h"
#include "InterSpec/UnitsConverterTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/FeatureMarkerWidget.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"
#include "InterSpec/EnergyCalPreserveWindow.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/LicenseAndDisclaimersWindow.h"

#include "InterSpec/D3TimeChart.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

#if( USE_DB_TO_STORE_SPECTRA )
#include "InterSpec/DbStateBrowser.h"
#endif

#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif

#if( IOS )
#include "target/ios/InterSpec/FileHandling.h"
#endif 

#if( INCLUDE_ANALYSIS_TEST_SUITE )
#include "InterSpec/SpectrumViewerTester.h"
#endif

#if( USE_GOOGLE_MAP )
#include "InterSpec/GoogleMap.h"
#endif

#if( USE_SEARCH_MODE_3D_CHART )
#include "InterSpec/SearchMode3DChart.h"
#endif

#if( USE_TERMINAL_WIDGET )
#include "InterSpec/TerminalWidget.h"
#endif

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
#include "InterSpec/SpecFileQueryWidget.h"
#endif

#if( SpecUtils_ENABLE_D3_CHART )
#include "SpecUtils/D3SpectrumExport.h"
#endif

#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#include "js/AppHtmlMenu.js"
#endif

#include "js/InterSpec.js"

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

using namespace Wt;
using namespace std;

std::mutex InterSpec::sm_staticDataDirectoryMutex;
std::string InterSpec::sm_staticDataDirectory = "data";

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
std::mutex InterSpec::sm_writableDataDirectoryMutex;
std::string InterSpec::sm_writableDataDirectory = "";
#endif  //if( not a webapp )



void log_error_message( const std::string &message, const std::string &source, const int priority )
{
  InterSpecApp *app = dynamic_cast<InterSpecApp*>( Wt::WApplication::instance() );
  if( app )
    app->svlog(message,source,priority);
  else
    std::cerr << source << ": " << message << std::endl;
}

namespace
{
  static const char * const PeakInfoTabTitle      = "Peak Manager";
  static const char * const GammaLinesTabTitle    = "Reference Photopeaks";
  static const char * const CalibrationTabTitle   = "Energy Calibration";
  static const char * const NuclideSearchTabTitle = "Nuclide Search";
  static const char * const FileTabTitle          = "Spectrum Files";
  
//#if( !BUILD_FOR_WEB_DEPLOYMENT )
//  const WTabWidget::LoadPolicy TabLoadPolicy = WTabWidget::LazyLoading;
//#else
  const WTabWidget::LoadPolicy TabLoadPolicy = WTabWidget::PreLoading;
//#endif

  void postSvlogHelper( const WString &msg, const int priority )
  {
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    if( app )
      app->svlog( msg, WString(), priority );
  }
  
  //adapted from: http://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
  struct csv_reader: std::ctype<char>
  {
    csv_reader(): std::ctype<char>(get_table()) {}
    static std::ctype_base::mask const* get_table()
    {
      static std::vector<std::ctype_base::mask> rc(table_size, std::ctype_base::mask());
      rc[','] = std::ctype_base::space;
    	rc[' '] = std::ctype_base::space;
      rc['\n'] = std::ctype_base::space;
      return &rc[0];
    }
  };//struct csv_reader
  

  
  //Returns -1 if you shouldnt add the peak to the hint peaks
  int add_hint_peak_pos( const std::shared_ptr<const PeakDef> &peak,
                        const std::deque< std::shared_ptr<const PeakDef> > &existing )
  {
    if( !peak )
      return -1;
    
    if( existing.empty() )
      return 0;
  
    const double sigma_frac = 0.7;
    const double mean = peak->mean();
    const double sigma = peak->sigma();
    
    deque< std::shared_ptr<const PeakDef> >::const_iterator pos;
    pos = lower_bound( existing.begin(), existing.end(),
                       peak, &PeakDef::lessThanByMeanShrdPtr );
    
    bool nearother = false;
    if( pos != existing.begin() )
      nearother |= ((fabs(((*(pos-1))->mean() - mean)/sigma) < sigma_frac));
    if( pos != existing.end() )
      nearother |= ( fabs(((*pos)->mean() - mean)/sigma) < sigma_frac );
    
    return (nearother ? -1 : static_cast<int>(pos - existing.begin()));
  }//int add_hint_peak_pos(...)
  
  
  bool try_update_hint_peak( const std::shared_ptr<const PeakDef> &newpeak,
                             std::shared_ptr<SpecMeas> &meas,
                             const set<int> &samples )
  {
    if( !meas || !newpeak )
      return false;
    
    std::shared_ptr< const SpecMeas::PeakDeque > hintPeaks
                                        = meas->automatedSearchPeaks( samples );
    
    if( !hintPeaks )
      return false;
    
    const int pos = add_hint_peak_pos( newpeak, *hintPeaks );
    if( pos >= 0 )
    {
      std::shared_ptr< SpecMeas::PeakDeque > newHintPeaks
        = std::make_shared<SpecMeas::PeakDeque>( hintPeaks->begin(), hintPeaks->end() );
      newHintPeaks->insert( newHintPeaks->begin() + pos, newpeak );
      meas->setAutomatedSearchPeaks( samples, newHintPeaks );
      return true;
    }//if( pos >= 0 )
    
    return false;
  }//try_update_hint_peak(...)
}//namespace


InterSpec::InterSpec( WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_peakModel( 0 ),
    m_spectrum( 0 ),
    m_timeSeries( 0 ),
    m_detectorToShowMenu( 0 ),
    m_mobileMenuButton(0),
    m_mobileBackButton(0),
    m_mobileForwardButton(0),
    m_notificationDiv(0),
    m_warnings( 0 ),
    m_warningsWindow( 0 ),
    m_fileManager( 0 ),
#if( USE_CSS_FLEX_LAYOUT )
    m_chartResizer( nullptr ),
    m_toolsResizer( nullptr ),
#else
    m_layout( 0 ),
    m_charts( nullptr ),
    m_chartResizer( nullptr ),
    m_toolsLayout( 0 ),
#endif
    m_menuDiv( 0 ),
    m_peakInfoDisplay( 0 ),
    m_peakInfoWindow( 0 ),
    m_peakEditWindow( 0 ),
    m_currentToolsTab( 0 ),
    m_toolsTabs( 0 ),
    m_energyCalTool( 0 ),
    m_energyCalWindow( 0 ),
    m_gammaCountDialog( 0 ),
    m_specFileQueryDialog( 0 ),
    m_shieldingSuggestion( 0 ),
    m_shieldingSourceFit( 0 ),
    m_shieldingSourceFitWindow( 0 ),
    m_materialDB( nullptr ),
    m_nuclideSearchWindow( 0 ),
    m_nuclideSearchContainer(0),
    m_nuclideSearch( 0 ),
    m_fileMenuPopup( 0 ),
    m_toolsMenuPopup( 0 ),
    m_helpMenuPopup( 0 ),
    m_displayOptionsPopupDiv( 0 ),
#if( USE_DB_TO_STORE_SPECTRA )
    m_saveState( 0 ),
    m_saveStateAs( 0 ),
    m_createTag( 0 ),
    m_currentStateID( -1 ),
#endif
    m_rightClickMenu( 0 ),
    m_rightClickEnergy( -DBL_MAX ),
    m_rightClickNuclideSuggestMenu( nullptr ),
    m_rightClickChangeContinuumMenu( nullptr ),
#if( USE_SAVEAS_FROM_MENU )
    m_downloadMenu( 0 ),
  m_downloadMenus{0},
#endif
  m_logYItems{0},
  m_toolTabsVisibleItems{0},
  m_backgroundSubItems{0},
  m_hardBackgroundSub( nullptr ),
  m_verticalLinesItems{0},
  m_horizantalLinesItems{0},
  m_showXAxisSliderItems{ nullptr },
  m_showYAxisScalerItems{ nullptr },
  m_compactXAxisItems{ nullptr },
  m_tabToolsMenuItems{0},
  m_featureMarkersShown{false},
  m_featureMarkers( nullptr ),
  m_featureMarkerMenuItem( nullptr ),
#if( USE_GOOGLE_MAP )
    m_mapMenuItem( 0 ),
#endif
#if( USE_SEARCH_MODE_3D_CHART )
    m_searchMode3DChart( 0 ),
#endif
    m_showRiidResults( nullptr ),
#if( USE_TERMINAL_WIDGET )
    m_terminalMenuItem( 0 ),
    m_terminal( 0 ),
    m_terminalWindow( 0 ),
#endif
    m_clientDeviceType( 0x0 ),
    m_referencePhotopeakLines( 0 ),
    m_referencePhotopeakLinesWindow( 0 ),
    m_licenseWindow( nullptr ),
    m_useInfoWindow( 0 ),
    m_decayInfoWindow( nullptr ),
    m_preserveCalibWindow( 0 ),
    m_renderedWidth( 0 ),
    m_renderedHeight( 0 ),
    m_colorPeaksBasedOnReferenceLines( true ),
    m_findingHintPeaks( false )
{
  //Initialization of the app (this function) takes about 11ms on my 2.6 GHz
  //  Intel Core i7, as of (20150316).
  
  addStyleClass( "InterSpec" );
  
  //Setting setLayoutSizeAware doesnt seem to add any appreciable network
  //  overhead, at least according to the chrome inspectrion panel when running
  //  locally
  setLayoutSizeAware( true );

  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  assert( app );
  //make it so InterSpec::instance() wont return nullptr for calls from within this constructor
  app->m_viewer = this;
  
  //for notification div
  m_notificationDiv = new WContainerWidget();
  m_notificationDiv->setStyleClass("qtipDiv");
  m_notificationDiv->setId("qtip-growl-container");
  
#if( BUILD_AS_ELECTRON_APP && !USE_ELECTRON_NATIVE_MENU )
  m_notificationDiv->addStyleClass( "belowMenu" );
#endif
  
  app->domRoot()->addWidget( m_notificationDiv );
  
  initHotkeySignal();
  
  // Try to grab the username.
  string username = app->getUserNameFromEnvironment();

#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( username == "" )
    username = InterSpecApp::userNameFromOS();
#endif
  
  //If they don't have a username, grab their stored UUID.
  if( username == "" )
  {
    try
    {
      username = wApp->environment().getCookie( "SpectrumViewerUUID" );
    }catch(...)
    {
      stringstream usernamestrm;
      usernamestrm << boost::uuids::random_generator()();
      username = usernamestrm.str();
      wApp->setCookie( "SpectrumViewerUUID", username, 3600*24*365 );
    }//try / catch
  }//if( no username )
  
  detectClientDeviceType();

  
  if( isPhone() )
    username += "_phone";
  else if( isTablet() )
    username += "_tablet";

// Set up the session; open the database.
  m_sql.reset( new DataBaseUtils::DbSession() );
  
  
  // Try to find the user in the database, if not, make a new entry
  {//begin interacting with DB
    DataBaseUtils::DbTransaction transaction( *m_sql );
    m_user = m_sql->session()->find< InterSpecUser >().where( "userName = ?" )
                         .bind( username ).limit(1).resultValue();
    
    if( m_user )
    {
      InterSpecUser::initFromDbValues( m_user, m_sql );
    }else
    {
      InterSpecUser::DeviceType type = InterSpecUser::Desktop;
      if( isPhone() )
        type = InterSpecUser::PhoneDevice;
      else if( isTablet() )
        type = InterSpecUser::TabletDevice;
    
      InterSpecUser *newuser = new InterSpecUser( username, type );
      m_user = m_sql->session()->add( newuser );
    
      InterSpecUser::initFromDefaultValues( m_user, m_sql );
    }//if( m_user ) / else
  
    m_user.modify()->startingNewSession();
    
    transaction.commit();
  }//end interacting with DB
  

  m_peakModel = new PeakModel( this );
  m_spectrum   = new D3SpectrumDisplayDiv();
  m_timeSeries = new D3TimeChart();
  
  m_peakModel->dataChanged().connect( m_spectrum, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->rowsRemoved().connect( m_spectrum, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->rowsInserted().connect( m_spectrum, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->layoutChanged().connect( m_spectrum, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  m_peakModel->modelReset().connect( m_spectrum, &D3SpectrumDisplayDiv::scheduleForegroundPeakRedraw );
  
  
  if( isPhone() )
  {
    //TODO: layoutSizeChanged(...) will trigger the compact axis anyway, but
    //      I should check if doing it here saves a roundtrip
    m_spectrum->setCompactAxis( true );
    m_timeSeries->setCompactAxis( true );
    
    //For iPhoneX we need to:
    //  X - adjust padding for safe area
    //  X - Add padding to mobile menu to cover the notch area.
    //  X - Change position of hamburger menu
    //  - set .Wt-domRoot background color to match the charts
    //  X - Detect orientation changes, and respond appropriately
    //  X - Get the Wt area to fill the entire width (and height)
    //  -Initial AuxWindow width/height correct
    
    LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsDoOrientationChange);
    
    const char *js = INLINE_JAVASCRIPT(
      window.addEventListener("orientationchange", Wt.WT.DoOrientationChange );
      setTimeout( Wt.WT.DoOrientationChange, 0 );
      setTimeout( Wt.WT.DoOrientationChange, 250 );  //JIC - doesnt look necassary
      setTimeout( Wt.WT.DoOrientationChange, 1000 );
    );
    
    doJavaScript( js );
  }//if( isPhone() )
  
  m_spectrum->setPeakModel( m_peakModel );
  
  m_nuclideSearch = new IsotopeSearchByEnergy( this, m_spectrum );
  m_nuclideSearch->setLoadLaterWhenInvisible(true);

  m_warnings = new WarningWidget( this );

  // Set up the energy calibration tool
  m_energyCalTool = new EnergyCalTool( this, m_peakModel );
  m_spectrum->rightMouseDragg().connect( m_energyCalTool, &EnergyCalTool::handleGraphicalRecalRequest );
  displayedSpectrumChanged().connect( m_energyCalTool, &EnergyCalTool::displayedSpecChangedCallback );

  m_fileManager = new SpecMeasManager( this );
  
  m_peakInfoDisplay = new PeakInfoDisplay( this, m_spectrum, m_peakModel );

  initMaterialDbAndSuggestions();
  
#if( BUILD_AS_ELECTRON_APP )
#if( USE_ELECTRON_NATIVE_MENU )
  const bool isAppTitlebar = false;
#else
  const bool isAppTitlebar = InterSpecApp::isPrimaryWindowInstance();
#endif
#else
  const bool isAppTitlebar = false; // !isMobile()
#endif
  
  
  WWidget *menuWidget = NULL;
  if( isMobile() )
  {
    m_mobileMenuButton = new WPushButton( "", wApp->domRoot() );

    m_mobileMenuButton->addStyleClass( "MobileMenuButton btn" );
    //Need to set z-index inline rather than in css so AuxWindows loaded before
    //  the CSS can be forced above the button
    m_mobileMenuButton->setZIndex( 8388635 );
   
    //hamburger
    PopupDivMenu *popup = new PopupDivMenu( m_mobileMenuButton, PopupDivMenu::AppLevelMenu );
    m_mobileMenuButton->removeStyleClass( "Wt-btn" );
    menuWidget = popup;
    
#if( !SpecUtils_ENABLE_D3_CHART )
    m_spectrum->setAvoidMobileMenu( true );
#endif
    
    //Add this transparent overlay when mobile left menu slides in, so that
    // we can capture the click to hide the menu.  If this is not there, the
    // canvas will take over and not propagate the event to close the menu down.
    doJavaScript( "$('body').append('<div class=\"mobilePopupMenuOverlay\" style=\"display: none;\"></div>');" );
    
    m_mobileBackButton = new WContainerWidget( wApp->domRoot() );
    m_mobileBackButton->addStyleClass( "MobilePrevSample btn" );
    m_mobileBackButton->setZIndex( 8388635 );
    m_mobileBackButton->clicked().connect( boost::bind(&InterSpec::handleUserIncrementSampleNum, this, SpecUtils::SpectrumType::Foreground, false) );
    m_mobileBackButton->setHidden(true);
      
    m_mobileForwardButton = new WContainerWidget( wApp->domRoot() );
    m_mobileForwardButton->addStyleClass( "MobileNextSample btn" );
    m_mobileForwardButton->setZIndex( 8388635 );
    m_mobileForwardButton->clicked().connect( boost::bind(&InterSpec::handleUserIncrementSampleNum, this, SpecUtils::SpectrumType::Foreground, true) );
    m_mobileForwardButton->setHidden(true);
  }else  //if( isMobile() )
  {
    m_menuDiv = new WContainerWidget();
    m_menuDiv->addStyleClass( "MenuBar" );
    
    if( !isAppTitlebar )
    {
      menuWidget = m_menuDiv;
      m_menuDiv->addStyleClass( "WebMenuBar" );
      
      WImage *snlLogo = new WImage( "InterSpec_resources/images/SNL_Stacked_Black_Blue_tiny.png", m_menuDiv );
      snlLogo->addStyleClass("SnlWebMenuBarLogo");
    }else
    {
      // If we're here, BUILD_AS_ELECTRON_APP is true, but leaving the compile switch easy to comment out for
      //  development purposes (i.e., there are redundant nested #if statements we can clear up once AppHtmlMenu.js
      //  development is more clear).
      
#if( BUILD_AS_ELECTRON_APP )
      app->useStyleSheet( "InterSpec_resources/AppHtmlMenu.css" );
      m_menuDiv->addStyleClass( "app-titlebar" );
      m_menuDiv->setHeight( 30 );
      
      WContainerWidget *dragRegion = new WContainerWidget( m_menuDiv );
      dragRegion->addStyleClass( "app-titlebar-drag-region" );
      
      //Add InterSpec icon to left side of menubar
      WContainerWidget *appIcon = new WContainerWidget( m_menuDiv );
      appIcon->addStyleClass( "window-appicon" );
      
      menuWidget = new WContainerWidget( m_menuDiv );
      menuWidget->addStyleClass( "AppMenuBtns" );
      menuWidget->setAttributeValue( "role", "menubar" );
      
      WText *menuTitle = new WText( "InterSpec", m_menuDiv );
      menuTitle->setInline( false );
      menuTitle->addStyleClass( "window-title" );
      
      WContainerWidget *titleStretcher = new WContainerWidget( m_menuDiv );
      titleStretcher->addStyleClass( "titlebar-stretcher" );
      
      WContainerWidget *windowControls = new WContainerWidget( m_menuDiv );
      windowControls->addStyleClass( "window-controls-container" );
      
      WContainerWidget *iconDiv = new WContainerWidget( windowControls );
      iconDiv->addStyleClass( "window-icon-bg" );
      WContainerWidget *minimizeIcon = new WContainerWidget( iconDiv );
      minimizeIcon->addStyleClass( "window-icon window-minimize InvertInDark" );
      
      iconDiv = new WContainerWidget( windowControls );
      iconDiv->addStyleClass( "window-icon-bg" );
      WContainerWidget *maximizeIcon = new WContainerWidget( iconDiv );
      maximizeIcon->addStyleClass( "window-icon window-max-restore window-maximize InvertInDark" );
      
      iconDiv = new WContainerWidget( windowControls );
      iconDiv->addStyleClass( "window-icon-bg" );
      WContainerWidget *closeIcon = new WContainerWidget( iconDiv );
      closeIcon->addStyleClass( "window-icon window-close InvertInDark" );
  
      WContainerWidget *resizer = new WContainerWidget( m_menuDiv );
      resizer->addStyleClass( "resizer top" );
      
      resizer = new WContainerWidget( m_menuDiv );
      resizer->addStyleClass( "resizer left" );
      
      LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsSetupAppTitleBar);
      LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsTitleBarChangeMaximized);
      
      doJavaScript( "Wt.WT.SetupAppTitleBar();" );
      
      auto toggleMaximize = [](){
        const char *js =
        "let elem = document.querySelector(\".Wt-domRoot\");"
        "if (!document.fullscreenElement) {"
        "  elem.requestFullscreen().catch(err => {"
        "    console.log( 'Error attempting to enable full-screen mode' );"
        "  });"
        "} else {"
        "  document.exitFullscreen();"
        "}";
        
#if( BUILD_AS_ELECTRON_APP )
        if( InterSpecApp::isPrimaryWindowInstance() )
          ElectronUtils::send_nodejs_message( "ToggleMaximizeWindow", "" );
        else
          wApp->doJavaScript( js );
#else
        wApp->doJavaScript( js );
#endif
      };//doMaximize
      
      maximizeIcon->clicked().connect( std::bind( toggleMaximize ) );
      appIcon->doubleClicked().connect( std::bind( toggleMaximize )  );
      dragRegion->doubleClicked().connect( std::bind( toggleMaximize )  );
      menuTitle->doubleClicked().connect( std::bind( toggleMaximize )  );
      
#if( BUILD_AS_ELECTRON_APP )
      minimizeIcon->clicked().connect( std::bind([](){
        ElectronUtils::send_nodejs_message( "MinimizeWindow", "" );
      }) );
      
      closeIcon->clicked().connect( std::bind([](){
        ElectronUtils::send_nodejs_message( "CloseWindow", "" );
      }) );
#endif //BUILD_AS_ELECTRON_APP
#else //#if( BUILD_AS_ELECTRON_APP - for dev purposes )
      assert( 0 );
#endif //#if( BUILD_AS_ELECTRON_APP - for dev pupropses )
    }//if( !isAppTitlebar ) / else
  }//if( isMobile() ) / else

  addFileMenu( menuWidget, isAppTitlebar );
  addDisplayMenu( menuWidget );
  
  addToolsMenu( menuWidget );
  addAboutMenu( menuWidget );

  
  /* Set the loading indicator so that it's the highest z-index, so always visible  */
  WDefaultLoadingIndicator *indicator = new Wt::WDefaultLoadingIndicator();
  indicator->addStyleClass( "LoadingIndicator" );
  app->setLoadingIndicator( indicator );
   
#if( USE_CSS_FLEX_LAYOUT )
  addStyleClass( "InterSpecFlex" );
  
  if( m_menuDiv )
    addWidget( m_menuDiv );
  
  addWidget( m_spectrum );
  
  m_chartResizer = new WContainerWidget( this );
  m_chartResizer->addStyleClass( "Wt-vsh2" );
  m_chartResizer->setHeight( 5 );
  
  addWidget( m_timeSeries );
  
  m_toolsResizer = new WContainerWidget( this );
  m_toolsResizer->addStyleClass( "Wt-vsh2" );
  m_toolsResizer->setHeight( 5 );
  
  {//begin make tool tabs
    m_nuclideSearchWindow = nullptr;
    m_referencePhotopeakLinesWindow = NULL;
    
    
    m_toolsTabs = new WTabWidget( this );
    
    CompactFileManager *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::LeftToRight );
    m_toolsTabs->addTab( compact, FileTabTitle, TabLoadPolicy );
    
    m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                   boost::placeholders::_1, boost::placeholders::_2 ) );
    
    m_toolsTabs->addTab( m_peakInfoDisplay, PeakInfoTabTitle, TabLoadPolicy );
    
    m_referencePhotopeakLines = new ReferencePhotopeakDisplay( m_spectrum,
                                                              m_materialDB.get(),
                                                              m_shieldingSuggestion,
                                                              this );
    setReferenceLineColors( nullptr );
    
    //PreLoading is necessary on the m_referencePhotopeakLines widget, so that the
    //  "Isotope Search" widget will work properly when a nuclide is clicked
    //  on to display its photpeaks
    //XXX In Wt 3.3.4 at least, the contents of m_referencePhotopeakLines
    //  are not actually loaded to the client until the tab is clicked, and I
    //  cant seem to get this to actually happen.
    //WMenuItem *refPhotoTab =
    m_toolsTabs->addTab( m_referencePhotopeakLines, GammaLinesTabTitle, TabLoadPolicy );
    
    m_toolsTabs->currentChanged().connect( this, &InterSpec::handleToolTabChanged );
    
    m_energyCalTool->setWideLayout();
    m_toolsTabs->addTab( m_energyCalTool, CalibrationTabTitle, TabLoadPolicy );
    
    m_toolsTabs->setHeight( 245 );
    
    assert( !m_nuclideSearchContainer );
    
    m_nuclideSearchContainer = new WContainerWidget();
    WGridLayout *isoSearchLayout = new WGridLayout();
    m_nuclideSearchContainer->setLayout( isoSearchLayout );
    isoSearchLayout->setContentsMargins( 0, 0, 0, 0 );
    isoSearchLayout->addWidget( m_nuclideSearch, 0, 0 );
    m_nuclideSearchContainer->setMargin( 0 );
    m_nuclideSearchContainer->setPadding( 0 );
    isoSearchLayout->setRowStretch( 0, 1 );
    isoSearchLayout->setColumnStretch( 0, 1 );
    
    //WMenuItem *nuclideTab =
    m_toolsTabs->addTab( m_nuclideSearchContainer, NuclideSearchTabTitle, TabLoadPolicy );
    //    const char *tooltip = "Search for nuclides with constraints on energy, "
    //                          "branching ratio, and half life.";
    //    HelpSystem::attachToolTipOn( nuclideTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
    //Make sure the current tab is the peak info display
    m_toolsTabs->setCurrentWidget( m_peakInfoDisplay );
    
    m_toolsTabs->setJavaScriptMember( WT_RESIZE_JS, "function(self,w,h,layout){ console.log( 'wtResize called for tools tab:', w, h, layout ); }" );
    
     // An attempt to call into wtResize from ResizeObserver - not working yet - tool tab contents dont expand...
    const string stackJsRef = m_toolsTabs->contentsStack()->jsRef();
    m_toolsTabs->setJavaScriptMember( "resizeObserver",
      "new ResizeObserver(entries => {"
        "for (let entry of entries) {"
          "if( entry.target && entry.target.wtResize ) {"
            "const w = entry.contentRect.width;"
            "const h = entry.contentRect.height;"
            "console.log( 'Got resize', entry.target.id, 'for {' + w + ',' + h + '}'  );"
            "entry.target.wtResize(entry.target, Math.round(w), Math.round(h), true);"
            "if( (h > 27) && (entry.target.id === '" + m_toolsTabs->id() + "') ){"
              "$('#" + m_toolsTabs->id() + " > .Wt-stack').each( function(i,el){ "
                  "$(el).height( Math.round(h - 27) );"
              "} );"
                                     
            "}"
            "if( (h > 35) && (entry.target.id === '" + m_toolsTabs->id() + "') ){"
              "$('#" + m_toolsTabs->id() + " > .Wt-stack > div').each( function(i,el){ "
                "console.log( 'Setting height to ' + (h - 35) ); "
                "$(el).height( Math.round(h - 35) );"
              "} );"
            "}"
            //"if(" + stackJsRef + " && " + stackJsRef + ".wtResize) {"
            //  + stackJsRef + ".wtResize(" + stackJsRef + ", Math.round(w), Math.round(h-27), true);"
            //  "console.log( 'Will call resize for', " + stackJsRef + " );"
            //"}"
          "}else console.log( 'no wtResize' );"
        "}"
           // "console.log( 'stack=', " + m_toolsTabs->contentsStack()->jsRef() + " );"
           // "console.log( 'stack wtResize=', " + m_toolsTabs->contentsStack()->jsRef() + ".wtResize );"
      "});"
    );
    
    m_toolsTabs->callJavaScriptMember( "resizeObserver.observe", m_toolsTabs->jsRef() );
    //for( int i = 0; i < m_toolsTabs->count(); ++i )
    //  m_toolsTabs->callJavaScriptMember( "resizeObserver.observe", m_toolsTabs->widget(i)->jsRef() );
    
  }//end make tool tabs
  
  // TODO: need to call wtResize of m_toolsTabs so they will get resized correctly
  // TODO: 

  
#else
  

  m_charts = new WContainerWidget();
      
      
  m_charts->addWidget( m_spectrum );
  m_charts->addStyleClass( "charts" );
  m_chartResizer = new WContainerWidget( m_charts );
  m_chartResizer->addStyleClass( "Wt-hrh2" );
  m_chartResizer->setHeight( isMobile() ? 10 : 5 );
  m_charts->addWidget( m_timeSeries );
  
  m_chartResizer->setHidden( m_timeSeries->isHidden() );
  
  LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsInitFlexResizer);
  m_charts->doJavaScript( "Wt.WT.InitFlexResizer('" + m_chartResizer->id() + "','" + m_timeSeries->id() + "');" );
  
  m_layout = new WGridLayout();
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  m_layout->setHorizontalSpacing( 0 );
  m_layout->setVerticalSpacing( 0 );
  
  setLayout( m_layout );

  if( m_menuDiv )
    m_layout->addWidget( m_menuDiv, m_layout->rowCount(), 0 );
  m_layout->addWidget( m_charts, m_layout->rowCount(), 0 );
  m_layout->setRowStretch( m_layout->rowCount() - 1, 1 );
#endif
  
  // No need to updated the default axis titles
  //m_timeSeries->setY1AxisTitle( "Gamma CPS" );
  //m_timeSeries->setY2AxisTitle( "Neutron CPS" );
  //m_timeSeries->setXAxisTitle( "Time of Measurement (seconds)", "Time (seconds)" );
  
  m_timeSeries->chartDragged().connect( this, &InterSpec::timeChartDragged );
  m_timeSeries->chartClicked().connect( this, &InterSpec::timeChartClicked );
  
  m_spectrum->setXAxisTitle( "Energy (keV)" );
  m_spectrum->setYAxisTitle( "Counts" );

  m_spectrum->enableLegend();
  m_spectrum->showHistogramIntegralsInLegend( true );
  m_spectrum->shiftAltKeyDragged().connect( this, &InterSpec::handleShiftAltDrag );

//  m_spectrum->rightClicked().connect( boost::bind( &InterSpec::createPeakEdit, this, boost::placeholders::_1) );
  m_rightClickMenu = new PopupDivMenu( nullptr, PopupDivMenu::TransientMenu );
  m_rightClickMenu->aboutToHide().connect( this, &InterSpec::rightClickMenuClosed );
  
  if( m_rightClickMenu->isMobile() )
  {
    m_rightClickMenu->addPhoneBackItem( nullptr );
  }else
  {
    m_rightClickMenu->setPositionScheme( Wt::Absolute );
    m_rightClickMenu->addStyleClass( " Wt-popupmenu Wt-outset" );
  }
  
  for( RightClickItems i = RightClickItems(0);
       i < kNumRightClickItems; i = RightClickItems(i+1) )
  {
    switch( i )
    {
      case kPeakEdit:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Peak Editor..." );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::peakEditFromRightClick );
        break;
      case kRefitPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Refit Peak" );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::refitPeakFromRightClick );
        break;
      case kRefitROI:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Refit ROI" );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::refitPeakFromRightClick );
        break;
      case kDeletePeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Delete Peak" );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::deletePeakFromRightClick );
        break;
      case kChangeNuclide:
        m_rightClickNuclideSuggestMenu = m_rightClickMenu->addPopupMenuItem( "Change Nuclide" );
        m_rightClickMenutItems[i] = m_rightClickNuclideSuggestMenu->parentItem();
        break;
      case kChangeContinuum:
      {
        m_rightClickChangeContinuumMenu = m_rightClickMenu->addPopupMenuItem( "Change Continuum" );
        m_rightClickMenutItems[i] = m_rightClickChangeContinuumMenu->parentItem();
        
        for( auto type = PeakContinuum::OffsetType(0);
            type <= PeakContinuum::External; type = PeakContinuum::OffsetType(type+1) )
        {
          WMenuItem *item = m_rightClickChangeContinuumMenu->addItem( PeakContinuum::offset_type_label(type) );
          item->triggered().connect( boost::bind( &InterSpec::handleChangeContinuumTypeFromRightClick, this, type ) );
        }//for( loop over PeakContinuum::OffsetTypes )
        break;
      }//case kChangeContinuum:
        
      case kAddPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Add Peak" );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::addPeakFromRightClick );
        break;
      case kShareContinuumWithLeftPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Combine Cont. Left" );
        m_rightClickMenutItems[i]->triggered().connect( boost::bind( &InterSpec::shareContinuumWithNeighboringPeak, this, true ) );
        break;
      case kShareContinuumWithRightPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Combine Cont. Right" );
        m_rightClickMenutItems[i]->triggered().connect( boost::bind( &InterSpec::shareContinuumWithNeighboringPeak, this, false ) );
        break;
        
      case kMakeOwnContinuum:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( "Own Continuum" );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::makePeakFromRightClickHaveOwnContinuum );
        break;
        
      case kNumRightClickItems:
        break;
    }//switch( i )
  }//for( loop over right click menu items )
  
  m_spectrum->rightClicked().connect( boost::bind( &InterSpec::handleRightClick, this,
                                                  boost::placeholders::_1, boost::placeholders::_2,
                                                  boost::placeholders::_3, boost::placeholders::_4 ) );
  m_spectrum->chartClicked().connect( boost::bind( &InterSpec::handleLeftClick, this,
                                                  boost::placeholders::_1, boost::placeholders::_2,
                                                  boost::placeholders::_3, boost::placeholders::_4 ) );
  
//  m_spectrum->controlKeyDragged().connect( boost::bind( &InterSpec::findPeakFromUserRange, this, boost::placeholders::_1, boost::placeholders::_2 ) );
  
  m_spectrum->shiftKeyDragged().connect( boost::bind( &InterSpec::excludePeaksFromRange, this,
                                                     boost::placeholders::_1,
                                                     boost::placeholders::_2 ) );
  m_spectrum->doubleLeftClick().connect( boost::bind( &InterSpec::searchForSinglePeak, this,
                                                     boost::placeholders::_1 ) );

  m_timeSeries->setHidden( true );
  m_chartResizer->setHidden( m_timeSeries->isHidden() );
  
#if( USE_OSX_NATIVE_MENU || USING_ELECTRON_NATIVE_MENU )
  if( InterSpecApp::isPrimaryWindowInstance() )
    m_menuDiv->hide();
#endif
  
  if( !isPhone() )
  {
    setToolTabsVisible( true );
  }else
  {
    //Disallow showing tool tabs when on a phone
    m_toolTabsVisibleItems[1]->setHidden(true);
    
    //Add Ref photopeaks, etc. to tool menu
    addToolsTabToMenuItems();
  }//If( start with tool tabs showing ) / else
 
#if( !ANDROID && !IOS )
  initDragNDrop();
#endif
  
#if( USE_DB_TO_STORE_SPECTRA )
  updateSaveWorkspaceMenu();
#endif
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !IOS && !BUILD_AS_ELECTRON_APP )
  initOsColorThemeChangeDetect();
#endif
  
  applyColorTheme( nullptr );
}//InterSpec( constructor )

InterSpec *InterSpec::instance()
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( !app )
    return nullptr;
  return app->viewer();
}//instance()

void InterSpec::setStaticDataDirectory( const std::string &dir )
{
  std::lock_guard<std::mutex> lock( sm_staticDataDirectoryMutex );
  
  if( !SpecUtils::is_directory(dir) )
    throw runtime_error( "InterSpec::setStaticDataDirectory(): " + dir + " is not a directory." );
  
  sm_staticDataDirectory = dir;

#ifdef _WIN32
  MassAttenuation::set_data_directory( SpecUtils::convert_from_utf8_to_utf16(dir) );
#else
  MassAttenuation::set_data_directory( dir );
#endif
  const string rctn_xml_file = SpecUtils::append_path( dir, "sandia.reactiongamma.xml" );
  if( !SpecUtils::is_file(rctn_xml_file) )
    throw runtime_error( "InterSpec::setStaticDataDirectory(): " + dir + " does not contain a sandia.reactiongamma.xml file." );
  ReactionGammaServer::set_xml_file_location( rctn_xml_file );
  
  const string decay_xml_file = SpecUtils::append_path( dir, "sandia.decay.xml" );
  if( !SpecUtils::is_file(decay_xml_file) )
    throw runtime_error( "InterSpec::setStaticDataDirectory(): " + dir + " does not contain a sandia.decay.xml file." );
  DecayDataBaseServer::setDecayXmlFile( decay_xml_file );
}//void setStaticDataDirectory( const std::string &dir )


std::string InterSpec::staticDataDirectory()
{
  std::lock_guard<std::mutex> lock( sm_staticDataDirectoryMutex );
  return sm_staticDataDirectory;
}

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
void InterSpec::setWritableDataDirectory( const std::string &dir )
{
  std::lock_guard<std::mutex> lock( sm_writableDataDirectoryMutex );
  
  if( !dir.empty() && !SpecUtils::is_directory(dir) )
    throw runtime_error( "InterSpec::setWritableDataDirectory(): " + dir + " is not a directory." );
  
  //Set the serial to module database file, if the user has one in thier
  //  application data directory.
  //Currently this is being done in the target specific code; it should all be
  //  moved here probably.
  //const vector<string> serial_db = SpecUtils::ls_files_in_directory( dir, "serial_to_model.csv" );
  //if( !serial_db.empty() )
  //  SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
  
  sm_writableDataDirectory = dir;
}//setWritableDataDirectory( const std::string &dir )

std::string InterSpec::writableDataDirectory()
{
  std::lock_guard<std::mutex> lock( sm_writableDataDirectoryMutex );
  
  if( sm_writableDataDirectory.empty() )
    throw runtime_error( "writableDataDirectory hasnt been set." );
  
  return sm_writableDataDirectory;
}//string writableDataDirectory()
#endif  //if( not a webapp )



InterSpec::~InterSpec() noexcept(true)
{
  //The deletion of the DOM root node will destroy all the AuxWindows we
  //  have open, but I am manually taking care of them below due to a crash
  //  I have been getting in the WApplication destructor for Wt 3.3.1-rc1

  cerr << "Destructing InterSpec from session '" << (wApp ? wApp->sessionId() : string("")) << "'" << endl;

  if( m_licenseWindow )
  {
    delete m_licenseWindow;
    m_licenseWindow = nullptr;
  }
  
  try
  {
    closeShieldingSourceFitWindow();
  }catch(...)
  {
    cerr << "Caught exception closing shielding source window - shouldnt have happened" << endl;
  }
  
  if( m_peakInfoDisplay )
  {
    if( m_toolsTabs && m_toolsTabs->indexOf(m_peakInfoDisplay)>=0 )
        m_toolsTabs->removeTab( m_peakInfoDisplay );
    if( m_peakInfoWindow )
      m_peakInfoWindow->contents()->removeWidget( m_peakInfoDisplay );
    delete m_peakInfoDisplay;
    m_peakInfoDisplay = nullptr;
  }//if( m_peakInfoDisplay )

  if( m_peakInfoWindow )
  {
    delete m_peakInfoWindow;
    m_peakInfoWindow = nullptr;
  }//if( m_peakInfoWindow )
  
  if( m_energyCalTool )
  {
    if( m_toolsTabs && m_toolsTabs->indexOf(m_energyCalTool)>=0 )
      m_toolsTabs->removeTab( m_energyCalTool );
    
    delete m_energyCalTool;
    m_energyCalTool = nullptr;
  }//if( m_energyCalTool )
  
  if( m_energyCalWindow )
  {
    delete m_energyCalWindow;
    m_energyCalWindow = nullptr;
  }//if( m_energyCalWindow )
    
  if( m_shieldingSourceFitWindow )
  {
    delete m_shieldingSourceFitWindow;
    m_shieldingSourceFitWindow = nullptr;
  }//if( m_shieldingSourceFitWindow )
  
  if( m_nuclideSearch )
  {
    delete m_nuclideSearch;
    m_nuclideSearch = nullptr;
  }//if( m_nuclideSearch )

  if( m_nuclideSearchWindow )
  {
    delete m_nuclideSearchWindow;
    m_nuclideSearchWindow = nullptr;
  }
  
  if( m_referencePhotopeakLines )
  {
    if( m_toolsTabs && m_toolsTabs->indexOf(m_referencePhotopeakLines)>=0 )
      m_toolsTabs->removeTab( m_referencePhotopeakLines );
    
    m_referencePhotopeakLines->clearAllLines();
    delete m_referencePhotopeakLines;
    m_referencePhotopeakLines = nullptr;
  }//if( m_referencePhotopeakLines )
  
  if( m_referencePhotopeakLinesWindow )
  {
    delete m_referencePhotopeakLinesWindow;
    m_referencePhotopeakLinesWindow = nullptr;
  }//if( m_referencePhotopeakLinesWindow )
  
  if( m_warnings )  //WarningWidget isnt necassarily parented, so we do have to manually delete it
  {
    if( m_warningsWindow )
      m_warningsWindow->stretcher()->removeWidget( m_warnings );
    delete m_warnings;
    m_warnings = nullptr;
  }//if( m_warnings )

  if( m_warningsWindow )
  {
    delete m_warningsWindow;
    m_warningsWindow = nullptr;
  }//if( m_warningsWindow )
  
  if( m_peakEditWindow )
  {
    delete m_peakEditWindow;
    m_peakEditWindow = nullptr;
  }//if( m_peakEditWindow )
  
  if( m_mobileMenuButton )
    delete m_mobileMenuButton;

  if (m_mobileBackButton )
    delete m_mobileBackButton;
  
  if (m_mobileForwardButton)
    delete m_mobileForwardButton;
  
  if( m_notificationDiv )
    delete m_notificationDiv;
  
  if( m_fileManager )
    delete m_fileManager;

  if( m_shieldingSuggestion )
    delete m_shieldingSuggestion;
  m_shieldingSuggestion = NULL;

  if( m_menuDiv )
  {
    delete m_menuDiv;
    m_menuDiv = nullptr;
  }//if( m_menuDiv )

  try
  {
    m_user.reset();
  }catch( ... )
  {
    cerr << "Caught unexpected exception doing m_user.reset()" << endl;
  }
  
  try
  {
    m_sql.reset();
  }catch( ... )
  {
    cerr << "Caught unexpected exception doing m_sql.reset()" << endl;
  }
}//InterSpec destructor



#if( SpecUtils_ENABLE_D3_CHART )

string InterSpec::print_d3_reference_gammas() const
{
  string result;
  
  if( m_referencePhotopeakLines )
  {
    result = m_referencePhotopeakLines->jsonReferenceLinesArray();
  }else
  {
    result = "[]";
  }//if( m_referencePhotopeakLines ) / else
  
  return result;
}//string InterSpec::print_d3_reference_gammas() const

//Temporary function (20160224) to aid in development of d3.js spectrum rendering
string InterSpec::print_d3_json() const
{
  std::ostringstream ostr;
  const char *q = "\"";  // for creating valid json format
    
  std::shared_ptr<const SpecUtils::Measurement> foreground = m_spectrum->data();
  std::shared_ptr<const SpecUtils::Measurement> background = m_spectrum->background();
  std::shared_ptr<const SpecUtils::Measurement> secondary  = m_spectrum->secondData();
    
  std::shared_ptr<SpecUtils::Measurement> data = m_spectrum->data();
  std::shared_ptr<SpecUtils::Measurement> back = m_spectrum->background();
  std::shared_ptr<SpecUtils::Measurement> second = m_spectrum->secondData();
    
  typedef deque< PeakModel::PeakShrdPtr > PeakDeque;
  string peakstring;
  
  // Update time
  ostr << "{\n\t" << q << "updateTime" << q << ":" << q << SpecUtils::to_iso_string(boost::posix_time::second_clock::local_time()) << q;
  
  // spectrum values
  ostr << ",\n\t" << q << "spectra" << q << ": [";
  

  if( data )
  {
    string title = Wt::WWebWidget::escapeText(data->title()).toUTF8();
    if( title != data->title() )
      data->set_title( title );  //JIC, proper escaping not implemented in SpecUtils yet.
    
    D3SpectrumExport::D3SpectrumOptions options;
    options.line_color = "black";
    options.display_scale_factor = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::Foreground);
    
    std::shared_ptr<const PeakDeque > peaks = m_dataMeasurement->peaks(m_displayedSamples);
    if( peaks )
    {
      vector<PeakModel::PeakShrdPtr> inpeaks( peaks->begin(), peaks->end() );
      options.peaks_json = PeakDef::peak_json( inpeaks, foreground );
    }
    
    D3SpectrumExport::write_spectrum_data_js( ostr, *data, options, 0, 1 );
  }

  if( back )
  {
    if( data )
      ostr << ",";
    
    string title = Wt::WWebWidget::escapeText(back->title()).toUTF8();
    if( title != back->title() )
      back->set_title( title );  //JIC, proper escaping not implemented in SpecUtils yet.
    
    D3SpectrumExport::D3SpectrumOptions options;
    options.line_color = "steelblue";
    options.display_scale_factor = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::Background);
    
    std::shared_ptr<const PeakDeque > peaks = m_backgroundMeasurement->peaks(m_backgroundSampleNumbers);
    if( peaks )
    {
      vector<PeakModel::PeakShrdPtr> inpeaks( peaks->begin(), peaks->end() );
      options.peaks_json = PeakDef::peak_json( inpeaks, foreground );
    }
    
    D3SpectrumExport::write_spectrum_data_js( ostr, *back, options, 1, -1 );
  }
  
  
  if( second )
  {
    if( data || back )
      ostr << ",";
    
    string title = Wt::WWebWidget::escapeText(second->title()).toUTF8();
    if( title != second->title() )
      second->set_title( title );  //JIC, proper escaping not implemented in SpecUtils yet.
    
    D3SpectrumExport::D3SpectrumOptions options;
    options.line_color = "green";
    options.display_scale_factor = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::SecondForeground);
    
    std::shared_ptr<const PeakDeque > peaks = m_backgroundMeasurement->peaks(m_sectondForgroundSampleNumbers);
    if( peaks )
    {
      vector<PeakModel::PeakShrdPtr> inpeaks( peaks->begin(), peaks->end() );
      options.peaks_json = PeakDef::peak_json( inpeaks, foreground );
    }
    
    D3SpectrumExport::write_spectrum_data_js( ostr, *second, options, 2, 1 );
  }
  
  // end of spectrum json
  ostr << "\n\t]";
  // end of json
  ostr << "\n}\n";
  return ostr.str();
}//std::string InterSpec::print_d3_json() const


D3SpectrumExport::D3SpectrumChartOptions InterSpec::getD3SpectrumOptions() const
{
  double xMin, xMax, yMin, yMax;
  displayedSpectrumRange(xMin, xMax, yMin, yMax);
  
  map<string,string> referc_line_json;
  
  if( m_referencePhotopeakLines )
    referc_line_json = m_referencePhotopeakLines->jsonReferenceLinesMap();
  
  D3SpectrumExport::D3SpectrumChartOptions options(
                                 /* title: */"Interactive Spectrum Development",
                                 /* xAxisTitle: */"Energy", /* yAxisTitle: */"Counts",
                                 /* dataTitle: */(m_spectrum->data() ? m_spectrum->data()->title() :
                                                  m_spectrum->background() ? m_spectrum->background()->title() :
                                                  m_spectrum->secondData() ? m_spectrum->secondData()->title() :
                                                  "Foreground"),
                                 /* useLogYAxis: */m_spectrum->yAxisIsLog(),
                                 /* showVerticalGridLines: */m_spectrum->verticalLinesShowing(),
                                 /* showHorizontalGridLines: */m_spectrum->horizontalLinesShowing(),
                                 /* legendEnabled: */m_spectrum->legendIsEnabled(),
                                 /* compactXAxis: */m_spectrum->isAxisCompacted(),
                                 /* showPeakUserLabels: */m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakUserLabel ),
                                 /* showPeakEnergyLabels: */m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakEnergyLabel ),
                                 /* showPeakNuclideLabels: */m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakNuclideLabel ),
                                 /* showPeakNuclideEnergyLabels: */ m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakNuclideEnergies ),
                                 /* showEscapePeakMarker: */m_featureMarkersShown[static_cast<int>(FeatureMarkerType::EscapePeakMarker)],
                                 /* showComptonPeakMarker: */m_featureMarkersShown[static_cast<int>(FeatureMarkerType::ComptonPeakMarker)],
                                 /* showComptonEdgeMarker: */m_featureMarkersShown[static_cast<int>(FeatureMarkerType::ComptonEdgeMarker)],
                                 /* showSumPeakMarker: */m_featureMarkersShown[static_cast<int>(FeatureMarkerType::SumPeakMarker)],
                                 /* backgroundSubtract: */m_spectrum->backgroundSubtract(),
                                 /* xMin: */xMin, /* xMax: */xMax,
                                 referc_line_json
  );
  

  return options;
}
#endif //#if( SpecUtils_ENABLE_D3_CHART )

/**
 Calls CompactFileManager's refactored and static method handleUserIncrementSampleNum
 */
void InterSpec::handleUserIncrementSampleNum( SpecUtils::SpectrumType type,
                                                      bool increment)
{
    CompactFileManager::handleUserIncrementSampleNum(type, increment, this, m_fileManager->model(), NULL);
}

std::shared_ptr<DataBaseUtils::DbSession> InterSpec::sql()
{
  return m_sql;
};


void InterSpec::layoutSizeChanged( int w, int h )
{
  m_renderedWidth = w;
  m_renderedHeight = h;
  
  const bool comactX = (h <= 420); //Apple iPhone 6+, 6s+, 7+, 8+
  if( (h > 20) && (comactX != m_spectrum->isAxisCompacted()) )
  {
    //Only set axis compact if it isnt already; dont ever set non-compact
    //  ToDo: For case screen is made small, so axis goes compact, then user
    //        makes screen large, currently wont un-compactify axis, even if
    //        user wants this; should evaluate if its worth checking the user
    //        preference about this.
    //const bool makeCompact = InterSpecUser::preferenceValue<bool>( "CompactXAxis", this );
    
    if( comactX && !m_spectrum->isAxisCompacted() )
      m_spectrum->setCompactAxis( comactX );
    
    if( comactX && !m_timeSeries->isAxisCompacted() )
      m_timeSeries->setCompactAxis( comactX );
  }
}//void layoutSizeChanged( int w, int h )


bool InterSpec::isMobile() const
{
  return (m_clientDeviceType & MobileClient);
}

bool InterSpec::isPhone() const
{
  return (m_clientDeviceType & PhoneClient);
}

bool InterSpec::isTablet() const
{
  return (m_clientDeviceType & TabletClient);
}

bool InterSpec::isDesktop() const
{
  return (m_clientDeviceType & DesktopBrowserClient);
}

bool InterSpec::isDedicatedApp() const
{
  return (m_clientDeviceType & DedicatedAppClient);
}

void InterSpec::detectClientDeviceType()
{
  m_clientDeviceType= 0x0;

  InterSpecApp *app= dynamic_cast<InterSpecApp *>( wApp );
  if( !app )
    return;


  const bool phone= app->isPhone();
  const bool tablet= app->isTablet();
  const bool mobile= app->isMobile();

  for( ClientDeviceType type= ClientDeviceType( 0x1 );
       type < NumClientDeviceType; type= ClientDeviceType( type << 1 ) )
  {
    switch( type )
    {
      case DesktopBrowserClient:
        if( !phone && !tablet && !mobile )
          m_clientDeviceType |= type;
        break;

      case PhoneClient:
        if( phone )
          m_clientDeviceType |= type;
        break;

      case TabletClient:
        if( tablet )
          m_clientDeviceType |= type;
        break;

      case MobileClient:
        if( mobile )
          m_clientDeviceType |= type;
        break;

      case HighBandwithClient:
#if( !BUILD_FOR_WEB_DEPLOYMENT )
        m_clientDeviceType|= type;
#else
        if( app->environment().clientAddress().find( "127.0.0" ) != string::npos )
          m_clientDeviceType|= type;
// should probably do some other testing here....
#endif
        //        std::cerr << "env.clientAddress()=" << env.clientAddress() <<
        //        std::endl;
        //        env.clientAddress() //134.252.17.78 (from anothother
        //        computer), 127.0.0.1 from same computer
        break;

      case DedicatedAppClient:
#if( !BUILD_FOR_WEB_DEPLOYMENT )
        m_clientDeviceType|= type;
#endif
        break;

      case NumClientDeviceType:
        break;
    } // switch( type )
  } // for( loop over ClientDeviceType enums )
} // void detectClientDeviceType()


#if( !ANDROID && !IOS )
void InterSpec::initDragNDrop()
{
  LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsFileUploadFcn);
  
  doJavaScript( "$('.Wt-domRoot').data('ForegroundUpUrl','" +
               m_fileManager->foregroundDragNDrop()->url() + "');" );
  
  doJavaScript( "$('.Wt-domRoot').data('BackgroundUpUrl','" +
               m_fileManager->backgroundDragNDrop()->url() + "');" );
  
  doJavaScript( "$('.Wt-domRoot').data('SecondUpUrl','" +
               m_fileManager->secondForegroundDragNDrop()->url() + "');" );
  
  doJavaScript( "Wt.WT.FileUploadFcn();" );
}//void InterSpec::initDragNDrop()
#endif //#if( !ANDROID && !IOS )


void InterSpec::initHotkeySignal()
{
  if( !!m_hotkeySignal )
    return;
  
  //We are specifying for the javascript to not be collected since the response
  //  will change if the tools tabs is shown or not.
  m_hotkeySignal.reset( new JSignal<unsigned int>( this, "hotkey", false ) );
  
  //sender.id was undefined in the following js, so had to work around this a bit
  const char *js = INLINE_JAVASCRIPT(
  function(id,e){
    
    if( !e || !e.key || e.metaKey || e.altKey || e.shiftKey || (typeof e.keyCode === 'undefined') )
      return;
    
    let code = 0;
    if( e.ctrlKey )
    {
      switch( e.key ){
        case '1': case '2': case '3': case '4': case '5': //Shortcuts to switch to the various tabs
        case 'h': // Help dialog
        case 'i': // Info about InterSpec
        case 'k': // Clear showing reference photopeak lines
        case 's': // Store
        case 'l': // Log/Linear
          if( $(".Wt-dialogcover").is(':visible') ) // Dont do shortcut when there is a blocking-dialog showing
            return;
          code = e.key.charCodeAt(0);
          break;
        default:  //Unused - nothing to see here - let the event propagate up
          return;
      }//switch( e.key )
    }else{
      switch( e.key ){
        case "Left":  case "ArrowLeft":  code = 37; break;
        case "Up":    case "ArrowUp":    code = 38; break;
        case "Right": case "ArrowRight": code = 39; break;
        case "Down":  case "ArrowDown":  code = 40; break;
        default:  //Unused - nothing to see here - let the event propagate up
          return;
      }//switch( e.key )
    
      // No menus are active - dont send the signal
      if( $(".MenuLabel.PopupMenuParentButton.active").length === 0 )
        return;
    }//if( e.ctrlKey ) / else
    
    e.preventDefault();
    e.stopPropagation();
    Wt.emit( id, {name:'hotkey'}, code );
  } );
  
  const string jsfcfn = string("function(e){var f=") + js + ";f('" + id() + "',e);}";
  wApp->declareJavaScriptFunction( "appKeyDown", jsfcfn );
  
  const string jsfcn = "document.addEventListener('keydown'," + wApp->javaScriptClass() + ".appKeyDown);";
  doJavaScript( jsfcn );
  
  m_hotkeySignal->connect( boost::bind( &InterSpec::hotKeyPressed, this, boost::placeholders::_1 ) );
}//void initHotkeySignal()


void InterSpec::hotKeyPressed( const unsigned int value )
{
  if( m_toolsTabs )
  {
    string expectedTxt;
    switch( value )
    {
      case '1': expectedTxt = FileTabTitle;          break;
      case '2': expectedTxt = PeakInfoTabTitle;      break;
      case '3': expectedTxt = GammaLinesTabTitle;    break;
      case '4': expectedTxt = CalibrationTabTitle;   break;
      case '5': expectedTxt = NuclideSearchTabTitle; break;
      
      case 'h': case 'H':
        HelpSystem::createHelpWindow( "getting-started" );
        break;
      
      case 'i': case 'I':
        showWelcomeDialog( true );
        break;
      
      case 'k': case 'K':
        // TODO: decide how to handle this if the current reference lines are from the
        //       "Nuclide Search" tool
        if( m_referencePhotopeakLines )
          m_referencePhotopeakLines->clearAllLines();
        break;
      
      case 'l': case 'L':
        setLogY( !m_spectrum->yAxisIsLog() );
        break;
        
      case 's': case 'S':
        stateSave();
        break;
        
      case 37: case 38: case 39: case 40:
        arrowKeyPressed( value );
        break;
    }//switch( value )
  
    if( expectedTxt.empty() )
      return;
  
    for( int i = 0; i < m_toolsTabs->count(); ++i )
    {
      if( m_toolsTabs->tabText(i).toUTF8() == expectedTxt )
      {
        m_toolsTabs->setCurrentIndex( i );
        handleToolTabChanged( i );
        break;
      }
    }//for( int i = 0; i < m_toolsTabs->count(); ++i )
  }else
  {
    switch( value )
    {
      case '1': showCompactFileManagerWindow(); break;
      case '2': showPeakInfoWindow();           break;
      case '3': showGammaLinesWindow();         break;
      case '4': showEnergyCalWindow();          break;
      case '5': showNuclideSearchWindow();      break;
      case 'l':
        setLogY( !m_spectrum->yAxisIsLog() );
        break;
      case 37: case 38: case 39: case 40:
        arrowKeyPressed( value );
        break;
    }//switch( value )
  }//if( tool tabs visible ) / else
}//void hotKeyPressed( const int value )


void InterSpec::arrowKeyPressed( const unsigned int value )
{
#if( USING_ELECTRON_NATIVE_MENU || USE_OSX_NATIVE_MENU || IOS || ANDROID )
  return;
#endif
  
  if( m_mobileMenuButton )
    return;
  
  const vector<PopupDivMenu *> menus{
    m_fileMenuPopup, m_displayOptionsPopupDiv, m_toolsMenuPopup, m_helpMenuPopup
  };
  
  bool foundActive = false;
  size_t activeIndex = 0;
  for( size_t i = 0; i < menus.size(); ++i )
  {
    if( menus[i] && menus[i]->isVisible() && menus[i]->parentButton() )
    {
      foundActive = true;
      activeIndex = i;
      break;
    }
  }
  
  if( !foundActive )
    return;
  
  const bool leftArrow  = (value == 37);
  const bool upArrow    = (value == 38);
  const bool rightArrow = (value == 39);
  const bool downArrow  = (value == 40);
  
  if( leftArrow || rightArrow )
  {
    size_t nextActiveIndex;
    if( leftArrow )
      nextActiveIndex = (activeIndex == 0) ? (menus.size() - 1) : (activeIndex - 1);
    else
      nextActiveIndex = ((activeIndex + 1) % menus.size());
    
    menus[activeIndex]->hide();
    menus[nextActiveIndex]->parentClicked();
  }else if( upArrow || downArrow )
  {
    // TODO: Have the up/down arrow keys select different menu items in the menu
    //   I havent gotten the up/down arrows working to select different menu items - maybe
    //    it needs to be purely JS.
    //    And once we do get the up/down arrows working, then need to listen for 'Enter' and
    //    trigger, again probably all in JS
  }//if( leftArrow || rightArrow ) / else if( upArrow || downArrow )
}//void arrowKeyPressed( const unsigned int value )


void InterSpec::rightClickMenuClosed()
{
  m_rightClickEnergy = -DBL_MAX;
}//void rightClickMenuClosed()


void InterSpec::peakEditFromRightClick()
{
  if( m_rightClickEnergy < -99999.9 )
    return;
  createPeakEdit( m_rightClickEnergy );
}//void peakEditFromRightClick()


std::shared_ptr<const PeakDef> InterSpec::nearestPeak( const double energy ) const
{
  // Why not just call m_peakModel->nearestPeak(energy)
  std::shared_ptr<const PeakDef> nearPeak;
  double minDE = std::numeric_limits<double>::infinity();
  
  const int nrow = m_peakModel->rowCount();
  for( int row = 0; row < nrow; ++row )
  {
    WModelIndex index = m_peakModel->index( row, 0 );
    const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index );
    const double dE = fabs( peak->mean() - energy );
    if( (dE < minDE)
        && ((energy > peak->lowerX()) && (energy < peak->upperX())) )
    {
      minDE = dE;
      nearPeak = peak;
    }//if( dE < minDE )
  }//for( int row = 0; row < nrow; ++row )

  return nearPeak;
}//std::shared_ptr<const PeakDef> nearestPeak() const

void InterSpec::refitPeakFromRightClick()
{
  PeakSearchGuiUtils::refit_peaks_from_right_click( this, m_rightClickEnergy );
}//void refitPeakFromRightClick()


void InterSpec::addPeakFromRightClick()
{
  std::shared_ptr<const SpecUtils::Measurement> dataH = m_spectrum->data();
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak
      || m_rightClickEnergy < peak->lowerX()
      || m_rightClickEnergy > peak->upperX()
      || !m_dataMeasurement
      || !dataH )
  {
    passMessage( "There was no ROI to add peak to",
                 "", WarningWidget::WarningMsgInfo );
    return;
  }//if( !peak )
  
  //get all the peaks that belong to this ROI
  typedef std::shared_ptr<const PeakContinuum> ContPtr;
  typedef map<ContPtr, std::vector<PeakDef > > ContToPeakMap;
  
  
  const auto origContinumm = peak->continuum();
  const std::shared_ptr<const deque<PeakModel::PeakShrdPtr>> allOrigPeaks = m_peakModel->peaks();
  
  ContToPeakMap contToPeaks;
  for( const PeakModel::PeakShrdPtr &thispeak : *allOrigPeaks )
  {
    if( thispeak )
      contToPeaks[thispeak->continuum()].push_back( *thispeak );
  }
  
  std::vector<PeakDef> origRoiPeaks = contToPeaks[peak->continuum()];
  
  //Make sure we dont try to mix data defined and gaussian defined peaks
  for( const PeakDef &p : origRoiPeaks )
  {
    if( !p.gausPeak() )
    {
      passMessage( "Sorry, cant add peaks to a ROI with a data defined peak.",
                  "", WarningWidget::WarningMsgInfo );
      return;
    }
  }//for( const PeakDef &p : origRoiPeaks )
  
  if( origRoiPeaks.empty() )
    throw runtime_error( "Logic error in InterSpec::addPeakFromRightClick()" );
  
  const double x0 = peak->lowerX();
  const double x1 = peak->upperX();
  const size_t lower_channel = dataH->find_gamma_channel( x0 );
  const size_t upper_channel = dataH->find_gamma_channel( x1 );
  const size_t nbin = ((upper_channel > lower_channel) ? (upper_channel-lower_channel) : size_t(0));
  
  double startingChi2, fitChi2;
  
  {//begin codeblock to evaluate startingChi2
    MultiPeakFitChi2Fcn chi2fcn( origRoiPeaks.size(), dataH,
                                peak->continuum()->type(),
                                lower_channel, upper_channel );
    startingChi2 = chi2fcn.evalRelBinRange( 0, chi2fcn.nbin(), origRoiPeaks );
  }//end codeblock to evaluate startingChi2
  
  
  bool inserted = false;
  vector< std::shared_ptr<PeakDef> > answer;
  for( PeakDef p : origRoiPeaks )
  {
    if( fabs(p.mean() - peak->mean()) < 0.01 )
    {
      inserted = true;
      PeakDef newpeak = p;
      newpeak.setMean( m_rightClickEnergy );
      newpeak.setSigma( newpeak.sigma() / sqrt(2.0) );
      newpeak.setAmplitude( 0.25*p.amplitude() );
      p.setAmplitude( 0.75*p.amplitude() );
      p.setSigma( p.sigma() / sqrt(2.0) );
      
      if( newpeak.mean() < p.mean() )
      {
        answer.push_back( std::make_shared<PeakDef>( newpeak ) );
        answer.push_back( std::make_shared<PeakDef>( p ) );
      }else
      {
        answer.push_back( std::shared_ptr<PeakDef>( new PeakDef(p) ) );
        answer.push_back( std::shared_ptr<PeakDef>( new PeakDef(newpeak) ) );
      }//if( newpeak.mean() < p.mean() ) / else
    }else
    {
//      p.setFitFor(PeakDef::Mean, false);
//      p.setFitFor(PeakDef::Sigma, false);
//      p.setFitFor(PeakDef::GaussAmplitude, false);
      answer.push_back( std::make_shared<PeakDef>( p ) );
    }
  }//for( PeakDef p : origRoiPeaks )
  
  if( !inserted )
    throw runtime_error( "Logic error 2 in InterSpec::addPeakFromRightClick()" );

  
  const MultiPeakInitialGuesMethod methods[] = { FromInputPeaks, UniformInitialGuess, FromDataInitialGuess };
  
  for( auto method : methods )
  {
    vector< std::shared_ptr<PeakDef> > orig_answer;
    for( auto p : answer )
      orig_answer.push_back( make_shared<PeakDef>( *p ) );
    
    findPeaksInUserRange( x0, x1, int(answer.size()), method, dataH,
                        m_dataMeasurement->detector(), answer, fitChi2 );
  
    std::vector<PeakDef> newRoiPeaks;
    for( size_t i = 0; i < answer.size(); ++i )
      newRoiPeaks.push_back( *answer[i] );
    
    MultiPeakFitChi2Fcn chi2fcn( newRoiPeaks.size(), dataH,
                                 peak->continuum()->type(),
                                 lower_channel, upper_channel );
    fitChi2 = chi2fcn.evalRelBinRange( 0, chi2fcn.nbin(), newRoiPeaks );
    
    if( fitChi2 < startingChi2 )
    {
      //cout << "Method " << method << " gave fitChi2=" << fitChi2 << ", which is better than initial startingChi2=" << startingChi2 << endl;
      break;
    }
    
    answer = orig_answer;
  }//for( auto method : methods )
  
//  could try to fix all peaks other than the new one, and the nearest one, do the
//  fit, then replace the other peaks to original fitFor state, and refit all of
//  them.
  
  const double dof = (nbin + 3*origRoiPeaks.size() + peak->continuum()->type());
  const double chi2Dof = fitChi2 / dof;
  
  cerr << "m_rightClickEnergy=" << m_rightClickEnergy << endl;
  cerr << "PreChi2=" << startingChi2 << ", PostChi2=" << fitChi2 << endl;
  cerr << "Got chi2Dof=" << chi2Dof << " when fitting for new peak" << endl;
  
  for( size_t i = 0; i < answer.size(); ++i )
  {
    cerr << "Peak " << i << " at " << answer[i]->mean() << " w/ width="
         << answer[i]->sigma() << ", and amp=" << answer[i]->amplitude() << endl;
  }
  
  //Remove old ROI peaks
  if( fitChi2 > startingChi2 )
  {
    passMessage( "Adding a peak did not improve overall ROI fit - not adding one.",
                "", WarningWidget::WarningMsgInfo );
    return;
  }
  
  
  std::map<std::shared_ptr<PeakDef>,PeakModel::PeakShrdPtr> new_to_orig_peaks;
  for( const PeakModel::PeakShrdPtr &thispeak : *allOrigPeaks )
  {
    if( !thispeak || (thispeak->continuum() != origContinumm) )
      continue;
    
    m_peakModel->removePeak( thispeak );
    
    //Associate new peak 
    double dist = 99999.9;
    std::shared_ptr<PeakDef> nearest_new;
    for( std::shared_ptr<PeakDef> newpeak : answer )
    {
      const double d = fabs( newpeak->mean() - thispeak->mean() );
      if( d < dist && !new_to_orig_peaks.count(newpeak) )
      {
        dist = d;
        nearest_new = newpeak;
      }
    }//
    if( nearest_new )
    {
      nearest_new->inheritUserSelectedOptions( *thispeak, true );
      new_to_orig_peaks[nearest_new] = thispeak;
    }
  }//for( const PeakModel::PeakShrdPtr &thispeak : *allOrigPeaks )
  
  
  for( size_t i = 0; i < answer.size(); ++i )
  {
    const bool isNew = !new_to_orig_peaks.count(answer[i]);
    InterSpec::addPeak( *answer[i], isNew );
  }
}//void addPeakFromRightClick()


void InterSpec::makePeakFromRightClickHaveOwnContinuum()
{
  const std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak || !data
      || m_rightClickEnergy < peak->lowerX()
      || m_rightClickEnergy > peak->upperX() )
  {
    passMessage( "There was no ROI to add peak to",
                "", WarningWidget::WarningMsgInfo );
    return;
  }//if( !peak )
  
  std::shared_ptr<const PeakContinuum> oldcont = peak->continuum();
  
  PeakDef newpeak(*peak);
  std::shared_ptr<PeakContinuum> cont
                                = std::make_shared<PeakContinuum>( *oldcont );
  newpeak.setContinuum( cont );
  
  cont->setRange( -1.0, -1.0 );
  size_t lowbin = findROILimit( newpeak, data, false );
  size_t upbin  = findROILimit( newpeak, data, true );
  double minx = data->gamma_channel_lower( lowbin );
  double maxx = data->gamma_channel_upper( upbin );
  cont->setRange( minx, maxx );
  
  m_peakModel->removePeak( peak );
  addPeak( newpeak, true );
  
  refitPeakFromRightClick();
  
  //Now we gotta refit the peaks that have the continuum we took this peak from
  std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks;
  peaks = m_peakModel->peaks();
  if( !peaks )
    return;
  
  //const cast hack!  (I *think* its okay since the peak is refit anyway)
  std::const_pointer_cast<PeakContinuum>(oldcont)->setRange( -1.0, -1.0 );
  
  minx = std::numeric_limits<double>::infinity();
  maxx = -std::numeric_limits<double>::infinity();
  vector<PeakModel::PeakShrdPtr > oldneighbors;
  for( deque< PeakModel::PeakShrdPtr >::const_iterator iter = peaks->begin();
       iter != peaks->end(); ++iter )
  {
    if( (*iter)->continuum() == oldcont )
    {
      oldneighbors.push_back( *iter );
      lowbin = findROILimit( *(*iter), data, false );
      upbin  = findROILimit( *(*iter), data, true );
      const double a = data->gamma_channel_lower( lowbin );
      const double z = data->gamma_channel_upper( upbin );
      minx = std::min( minx, a );
      maxx = std::max( maxx, z );
    }//if( (*iter)->continuum() == oldcont )
  }//for( iterate over peaks )
  
  if( oldneighbors.empty() )
    return;
  
  //const cast hack!  (I *think* its okay since the peak is refit anyway)
  std::const_pointer_cast<PeakContinuum>(oldcont)->setRange( minx, maxx );
  const double oldenergy = m_rightClickEnergy;
  m_rightClickEnergy = oldneighbors[0]->mean();
  refitPeakFromRightClick();
  m_rightClickEnergy = oldenergy;
}//void makePeakFromRightClickHaveOwnContinuum()


void InterSpec::handleChangeContinuumTypeFromRightClick( const int continuum_type )
{
  PeakSearchGuiUtils::change_continuum_type_from_right_click( this, m_rightClickEnergy,
                                                             continuum_type );
}//InterSpec::handleChangeContinuumTypeFromRightClick(...)


void InterSpec::shareContinuumWithNeighboringPeak( const bool shareWithLeft )
{
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak || m_rightClickEnergy < peak->lowerX() || m_rightClickEnergy > peak->upperX() )
  {
    passMessage( "There was no ROI to add peak to", "", WarningWidget::WarningMsgInfo );
    return;
  }//if( !peak )
  
  std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > peaks;
  peaks = m_peakModel->peaks();
  
  if( !peaks )
    return;
  
  std::shared_ptr<const PeakDef> peaktoshare;
  deque<PeakModel::PeakShrdPtr >::const_iterator iter;
  
//  boost::function<bool(const PeakModel::PeakShrdPtr &, const PeakModel::PeakShrdPtr &)> meansort;
//  meansort = boost::bind( &PeakModel::compare, boost::placeholders::_1, boost::placeholders::_2, PeakModel::kMean, Wt::AscendingOrder );
//  iter = lower_bound( peaks->begin(), peaks->end(), peak, meansort );
  iter = std::find( peaks->begin(), peaks->end(), peak );
  
  if( iter==peaks->end() || (*iter)!=peak )
    throw runtime_error( "shareContinuumWithNeighboringPeak: "
                         "error searching for peak I should have found" );
  
  if( shareWithLeft )
  {
    if( iter == peaks->begin() )
      return;
    peaktoshare = *(iter-1);
  }else
  {
    if( (iter+1) == peaks->end() )
      return;
    peaktoshare = *(iter+1);
  }//if( shareWithLeft ) / else
  
  if( peak->continuum() == peaktoshare->continuum() )
    return;
  
  if( !peaktoshare->continuum()->isPolynomial() || !peaktoshare->gausPeak() )
    return;
  
  vector<PeakModel::PeakShrdPtr> leftpeaks, rightpeaks;
  vector<PeakModel::PeakShrdPtr> &pp = (shareWithLeft ? rightpeaks : leftpeaks);
  vector<PeakModel::PeakShrdPtr> &op = (shareWithLeft ? leftpeaks : rightpeaks);
  
  for( auto iter = peaks->begin(); iter != peaks->end(); ++iter )
  {
    if( (*iter)->continuum() == peak->continuum() )
      pp.push_back( *iter );
    else if( (*iter)->continuum() == peaktoshare->continuum() )
      op.push_back( *iter );
  }//for( iterate over peaks )
  
  
  std::shared_ptr<const PeakContinuum> oldcontinuum = peaktoshare->continuum();
  std::shared_ptr<PeakContinuum> continuum = std::make_shared<PeakContinuum>( *oldcontinuum );
  
  const double minx = min( leftpeaks[0]->lowerX(), rightpeaks[0]->lowerX() );
  const double maxx = max( leftpeaks[0]->upperX(), rightpeaks[0]->upperX() );
  continuum->setRange( minx, maxx );
  
  //Make sure the order polynomial is at least the minimum of either side
  if( leftpeaks[0]->continuum()->isPolynomial()
      && continuum->type() < leftpeaks[0]->continuum()->type() )
    continuum->setType( leftpeaks[0]->continuum()->type() );
  
  if( rightpeaks[0]->continuum()->isPolynomial()
     && continuum->type() < rightpeaks[0]->continuum()->type() )
    continuum->setType( rightpeaks[0]->continuum()->type() );
  
  
  for( PeakModel::PeakShrdPtr &p : leftpeaks )
  {
    PeakDef newpeak( *p );
    newpeak.setContinuum( continuum );
    m_peakModel->removePeak( p );
    addPeak( newpeak, false );
  }//for( PeakModel::PeakShrdPtr &p : leftpeaks )
  
  for( PeakModel::PeakShrdPtr &p : rightpeaks )
  {
    PeakDef newpeak( *p );
    newpeak.setContinuum( continuum );
    m_peakModel->removePeak( p );
    addPeak( newpeak, false );
  }//for( PeakModel::PeakShrdPtr &p : leftpeaks )

  
//Instead of the more computationally expensive stuff (lots of signals get
//  emmitted, tables re-rendered, and charts re-drawn - although hopefully
//  mostly lazili) above, you could do the following not very nice const casts
//  std::shared_ptr<PeakContinuum> continuum;
//  continuum = boost::const_pointer_cast<PeakContinuum>( peaktoshare->continuum() );
//  if( shareWithLeft )
//    continuum->setRange( peaktoshare->lowerX(), peak->upperX() );
//  else
//    continuum->setRange( peak->lowerX(), peaktoshare->upperX() );
//  boost::const_pointer_cast<PeakDef>(peak)->setContinuum( continuum );

  refitPeakFromRightClick();
}//void shareContinuumWithNeighboringPeak( const bool shareWithLeft )


void InterSpec::deletePeakFromRightClick()
{
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak )
  {
    passMessage( "There was no peak to delete", "", WarningWidget::WarningMsgInfo );
    return;
  }

  m_peakModel->removePeak( peak );
}//void deletePeakFromRightClick()




void InterSpec::setPeakNuclide( const std::shared_ptr<const PeakDef> peak,
                                     std::string nuclide )
{
  //TODO: should probably add in some error logging or something
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  if( !peak || nuclide.empty() || !db )
  {
    cerr << "InterSpec::setPeakNuclide(): invalid input" << endl;
    return;
  }//if( !peak || nuclide.empty() || !db )
  
  WModelIndex index = m_peakModel->indexOfPeak( peak );
  if( !index.isValid() )
  {
    cerr << "InterSpec::setPeakNuclide(): couldnt find index for peak"
         << endl;
    return;
  }//if( !index.isValid() )
  
  index = m_peakModel->index( index.row(), PeakModel::kIsotope );
  
  m_peakModel->setData( index, boost::any(WString(nuclide)) );
  
  const PeakModel::PeakShrdPtr &p = m_peakModel->peak(index);
  
  //Reactions not implemented yet.
  if( !p || !p->hasSourceGammaAssigned() || p->reaction() )
    return;
  
  //Check if source is same as showing reference line, if so, set peak color to that.
  //ToDo: implement this if the user types into the Peak Manager.
  vector<ReferenceLineInfo> refLines;
  if( m_referencePhotopeakLines )
    refLines = m_referencePhotopeakLines->showingNuclides();
  for( const auto &lines : refLines )
  {
    if( (lines.nuclide || lines.element )
       && lines.nuclide==p->parentNuclide()
       && lines.element==p->xrayElement() )
    {
      index = m_peakModel->index( index.row(), PeakModel::kPeakLineColor );
      m_peakModel->setData( index, boost::any(WString(lines.lineColor.cssText())) );
    }
  }//for( const auto &lines : refLines )
}//void setPeakNuclide(...)


void InterSpec::updateRightClickNuclidesMenu(
                                  const std::shared_ptr<const PeakDef> peak,
                                  std::shared_ptr< std::vector<std::string> > nuclides )
{
  //TODO: should probably add in some error logging or something
  PopupDivMenu *menu = m_rightClickNuclideSuggestMenu;
  
  if( !peak || !nuclides || m_rightClickEnergy < 0.0 ||!menu )
    return;
  
  for( WMenuItem *item : menu->items() )
  {
    if( !item->hasStyleClass("PhoneMenuBack") && !item->hasStyleClass("PhoneMenuClose") )
      delete item;
  }//for( WMenuItem *item : menu->items() )

  
  for( size_t i = 0; i < nuclides->size(); ++i )
  {
    const std::string &nuc = (*nuclides)[i];
    
    if( i==0 && nuc.size()
        && (peak->parentNuclide() || peak->reaction() || peak->xrayElement())  )
    {
      PopupDivMenuItem *item = menu->addMenuItem( nuc, "", true );
      item->setAttributeValue("style", "background: grey; color: white;"
                                       + item->attributeValue("style"));
    }else if( nuc.size() )
    {
      PopupDivMenuItem *item = menu->addMenuItem( nuc, "", true );
      
      if( nuc[0]=='(' && nuc[nuc.size()-1]==')' )
        item->disable();
      else
        item->triggered().connect( boost::bind( &InterSpec::setPeakNuclide, this, peak, nuc ) );
    }else if( i != (nuclides->size()-1) )
    {
      menu->addSeparator();
    }
  }//for( size_t i = 0; i < nuclides->size(); ++i )
  
  //Adjust the menu so its all visible, since it may be larger now, but first
  //  check to make sure the menu is still visible.  We have to do this client
  //  side since the server is unaware whats hidden/visible.
  WMenuItem *parentItem = menu->parentItem();
  WMenu *parentMenu = (parentItem ? parentItem->parentMenu() : (WMenu *)0);
  if( !isMobile() && parentMenu )
  {
    //The positioning margin is copied from WPopupMenu.js and could change in
    //  future versions of Wt (current 3.3.4).  The AdjustTopPos(...) function
    //  is from PopupDiv.cpp, and makes sure Wts positioning is reasonable.
    doJavaScript(  "setTimeout(function(){try{"
                     "var u=" + menu->jsRef() + ", p=" + parentItem->jsRef() + ";"
                     "if(u.style.display==='block'){"
                       "var m=Wt.WT.px(u,'paddingTop')+Wt.WT.px(u,'borderTopWidth');"
                       "Wt.WT.positionAtWidget(u.id,p.id, Wt.WT.Horizontal,-m);"
                       "Wt.WT.AdjustTopPos(u.id);"
                     "}"
                   "}catch(e){console.log('Failed re-positioning menu');}},0);" );
  }//if( !isMobile() && parentMenu )
  
  if( nuclides->empty() )
  {
    PopupDivMenuItem *item = menu->addMenuItem( "No Suggestions", "", true );
    item->setAttributeValue("style", "background: grey; color: white;"
                            + item->attributeValue("style"));
  }
  
  WApplication *app = wApp;
  if( app )
    app->triggerUpdate();
}//updateRightClickNuclidesMenu()


void InterSpec::handleLeftClick( double energy, double counts,
                                      double pageX, double pageY )
{
  // For touch screen non-mobile devices, the right-click menu may be showing
  //  since when it is touch activated (by holding down for >600ms), there is
  //  no mouse-leave signal to cause it to leave.
  //  Note: I dont think m_rightClickMenu is ever null, but we'll check JIC
  if( m_rightClickMenu && !m_rightClickMenu->isHidden() )
    m_rightClickMenu->hide();

  if( (m_toolsTabs && (m_currentToolsTab == m_toolsTabs->indexOf(m_nuclideSearchContainer)))
      || m_nuclideSearchWindow )
  {
    setIsotopeSearchEnergy( energy );
    return;
  }
  
#if( USE_TERMINAL_WIDGET )
  if( m_toolsTabs && m_terminal
      && (m_toolsTabs->currentIndex() == m_toolsTabs->indexOf(m_terminal)) )
  {
    m_terminal->chartClicked(energy,counts,pageX,pageY);
  }
#endif
  
  if( (m_toolsTabs && (m_currentToolsTab == m_toolsTabs->indexOf(m_peakInfoDisplay)))
     || m_peakInfoWindow )
  {
    m_peakInfoDisplay->handleChartLeftClick( energy );
    return;
  }
}//void handleLeftClick(...)




void InterSpec::handleRightClick( double energy, double counts,
                                  double pageX, double pageY )
{
  if( !m_dataMeasurement )
    return;
  
  const std::shared_ptr<const PeakDef> peak = nearestPeak( energy );
  m_rightClickEnergy = energy;
  
  
  typedef std::shared_ptr<const std::deque< PeakModel::PeakShrdPtr > > PeakContainer_t;
  PeakContainer_t peaks = m_peakModel->peaks();
  
  if( !peaks || !peak )
    return;
  
  //We actually need to make a copy of peaks since we will post to outside the
  //  main thread where the deque is expected to remain constant, however, the
  //  PeakModel deque may get changed by like the user deleting a peak.
  peaks = make_shared< std::deque< PeakModel::PeakShrdPtr > >( begin(*peaks), end(*peaks) );
  
  
  //see how many other peaks share ROI
  size_t npeaksInRoi = 0;
  for( const PeakModel::PeakShrdPtr &p : *peaks )
    npeaksInRoi += (p->continuum() == peak->continuum());
  
  std::deque< PeakModel::PeakShrdPtr >::const_iterator iter;
  
  for( RightClickItems i = RightClickItems(0);
      i < kNumRightClickItems; i = RightClickItems(i+1) )
  {
    switch( i )
    {
      case kPeakEdit: case kDeletePeak: case kAddPeak:
//        m_rightClickMenutItems[i]->setHidden( !peak );
      break;
        
      case kChangeContinuum:
      {
        if( !m_rightClickChangeContinuumMenu )
          break;
        
        // Disable current continuum type, enable all others
        const vector<WMenuItem *> items = m_rightClickChangeContinuumMenu->items();
        const char *labelTxt = PeakContinuum::offset_type_label( peak->continuum()->type() );
        for( WMenuItem *item : items )
          item->setDisabled( item->text() == labelTxt );
        
        break;
      }//case kChangeContinuum:
        
      case kChangeNuclide:
      {
        if( m_rightClickNuclideSuggestMenu )
        {
          for( WMenuItem *item : m_rightClickNuclideSuggestMenu->items() )
          {
            if( !item->hasStyleClass("PhoneMenuBack")
                && !item->hasStyleClass("PhoneMenuClose") )
              delete item;
          }//for( WMenuItem *item : m_rightClickNuclideSuggestMenu->items() )
          
          m_rightClickNuclideSuggestMenu->addMenuItem( "...calculating...",
                                                       "", false );
          
          Wt::WServer *server = Wt::WServer::instance();
          assert( server );
          if( server )  //this should always be true
          {
            Wt::WIOService &io = server->ioService();
            std::shared_ptr< vector<string> > candidates
                                       = std::make_shared<vector<string> >();
            boost::function<void(void)> updater = wApp->bind(
                       boost::bind(&InterSpec::updateRightClickNuclidesMenu,
                                   this, peak, candidates ) );
            
            std::shared_ptr<const SpecUtils::Measurement> hist = displayedHistogram( SpecUtils::SpectrumType::Foreground );
            std::shared_ptr<const SpecMeas> meas = measurment( SpecUtils::SpectrumType::Foreground );
            std::shared_ptr<const DetectorPeakResponse> detector;
            if( !!meas )
              detector = meas->detector();
        
            PeakContainer_t hintpeaks;
            if( m_dataMeasurement )
              hintpeaks = m_dataMeasurement->automatedSearchPeaks( m_displayedSamples );
//            if( !hintpeaks || (!!peaks && hintpeaks->size() <= peaks->size()) )
//              hintpeaks = peaks;
            
            //Make a copy of hintpeaks deque since we will be posting to outside
            //  the main thread, where the original deque could actually get
            //  modified by another call, but what we're posting to expected it
            //  to remain constant.
            if( hintpeaks )
              hintpeaks = make_shared< std::deque< PeakModel::PeakShrdPtr > >( begin(*hintpeaks), end(*hintpeaks) );
            
            vector<ReferenceLineInfo> refLines;
            if( m_referencePhotopeakLines )
              refLines = m_referencePhotopeakLines->showingNuclides();
            
            const string session_id = wApp->sessionId();
            
            boost::function<void(void)> worker = [=](){
              IsotopeId::populateCandidateNuclides( hist, peak, hintpeaks,
                      peaks, refLines, detector, session_id, candidates, updater );
            };
            
            io.boost::asio::io_service::post( worker );
          }//if( server )
        }//if( menu )
        else
          cerr << "Serious error getting WMenu" << endl;
        
        break;
      }//case kChangeNuclide:
        
      case kRefitPeak:
        m_rightClickMenutItems[i]->setHidden( !peak->gausPeak() || npeaksInRoi>1 );
      break;
        
      case kRefitROI:
        m_rightClickMenutItems[i]->setHidden( !peak->gausPeak() || npeaksInRoi<2 );
      break;
        
      case kShareContinuumWithLeftPeak:
      {
        iter = std::find( peaks->begin(), peaks->end(), peak );
        if( iter == peaks->begin() )
        {
          m_rightClickMenutItems[i]->setHidden( true );
          break;
        }//if( iter == peaks.begin() )
        
        if( iter==peaks->end() || (*iter)!=peak )
          throw runtime_error( "InterSpec::handleRightClick: "
                               "error searching for peak I should have found 1" );
        
        PeakModel::PeakShrdPtr leftpeak = *(iter - 1);
        if( leftpeak->continuum() == peak->continuum()
            || !leftpeak->continuum()->isPolynomial() )
        {
          m_rightClickMenutItems[i]->setHidden( true );
          break;
        }//if( it alread shares a continuum with the left peak )
        
        const double leftupper = leftpeak->upperX();
        const double rightlower = peak->lowerX();
        const double rightupper = peak->upperX();
        
        //The below 2.0 is arbitrary
        const bool show = ((rightlower-leftupper) < 2.0*(rightupper-rightlower));
        m_rightClickMenutItems[i]->setHidden( !show );
        break;
      }//case kShareContinuumWithLeftPeak:
        
        
      case kShareContinuumWithRightPeak:
      {
        iter = std::find( peaks->begin(), peaks->end(), peak );
        
        if( iter==peaks->end() || (*iter)!=peak )
          throw runtime_error( "InterSpec::handleRightClick: "
                               "error searching for peak I should have found 2" );

        if( (iter+1) == peaks->end() )
        {
          m_rightClickMenutItems[i]->setHidden( true );
          break;
        }//if( iter == peaks.begin() )
        
        
        PeakModel::PeakShrdPtr rightpeak = *(iter + 1);
        if( rightpeak->continuum() == peak->continuum()
            || !rightpeak->continuum()->isPolynomial() )
        {
          m_rightClickMenutItems[i]->setHidden( true );
          break;
        }//if( it alread shares a continuum with the left peak )
        
        const double rightlower = rightpeak->lowerX();
        const double leftlower = peak->lowerX();
        const double leftupper = peak->upperX();
        
        //The below 2.0 is arbitrary
        const bool show = ((rightlower-leftupper) < 2.0*(leftupper-leftlower));
        m_rightClickMenutItems[i]->setHidden( !show );
        break;
      }//case kShareContinuumWithRightPeak:
      
      
      case kMakeOwnContinuum:
      {
        bool shares = false;
        std::shared_ptr<const PeakContinuum> cont = peak->continuum();
        for( PeakModel::PeakShrdPtr p : *peaks )
          shares = (shares || (p!=peak && p->continuum()==cont) );
        m_rightClickMenutItems[i]->setHidden( !shares );
        break;
      }//kMakeOwnContinuum
        
      case kNumRightClickItems:
      break;
    }//switch( i )
  }//for( loop over right click menu items )
  
  if( m_rightClickMenu->isMobile() )
    m_rightClickMenu->showMobile();
  else
    m_rightClickMenu->popup( WPoint(pageX,pageY) );
}//void handleRightClick(...)


void InterSpec::createPeakEdit( double energy )
{
  if( m_peakEditWindow )
  {
    m_peakEditWindow->peakEditor()->changePeak( energy );
  }else
  {
    m_peakEditWindow = new PeakEditWindow( energy, m_peakModel, this );
    m_peakEditWindow->editingDone().connect( this, &InterSpec::deletePeakEdit );
    m_peakEditWindow->finished().connect( this, &InterSpec::deletePeakEdit );
    m_peakEditWindow->resizeToFitOnScreen();
    
    PeakEdit *editor = m_peakEditWindow->peakEditor();
    if( !editor->isEditingValidPeak() )
    {
      delete m_peakEditWindow;
      m_peakEditWindow = 0;
    }//if( !editor->isEditingValidPeak() )
  }//if( m_peakEditWindow ) / else
  

}//void createPeakEdit( Wt::WMouseEvent event )


void InterSpec::deletePeakEdit()
{
  if( m_peakEditWindow )
  {
    delete m_peakEditWindow;
    m_peakEditWindow = NULL;
  }
}//void deletePeakEdit()



void InterSpec::setIsotopeSearchEnergy( double energy )
{
  if( !m_nuclideSearch )
    return;
  
  double sigma = -1.0;

  if( m_toolsTabs )
  {
    if( m_currentToolsTab != m_toolsTabs->indexOf(m_nuclideSearchContainer) )
      return;
  }else if( !m_nuclideSearchWindow )
  {
    return;
  }
  
  //check to see if this is within 3.0 sigma of a peak, and if so, set the
  //  energy to the mean of that peak
  PeakModel::PeakShrdPtr peak = m_peakModel->nearestPeak( energy );
  
  if( !!peak )
  {
    const double width = peak->gausPeak() ? peak->fwhm() : 0.5*peak->roiWidth();
    if( (fabs(peak->mean()-energy) < (3.0/2.634)*width) )
    {
      energy = peak->mean();
      sigma = width;
      
      //For low res spectra make relatively less wide
      const auto spectrum = displayedHistogram(SpecUtils::SpectrumType::Foreground);
      if( m_dataMeasurement && !PeakFitUtils::is_high_res(spectrum) )
        sigma *= 0.35;
    }//if( within 3 sigma of peak )
  }//if( !!peak )
  
  m_nuclideSearch->setNextSearchEnergy( energy, sigma );
}//void setIsotopeSearchEnergy( double energy );



void InterSpec::setFeatureMarkerOption( const FeatureMarkerType option, const bool show )
{
  m_featureMarkersShown[static_cast<int>(option)] = show;
  
  m_spectrum->setFeatureMarkerOption( option, show );
}//setFeatureMarkerOption(...)


bool InterSpec::showingFeatureMarker( const FeatureMarkerType option )
{
  return m_featureMarkersShown[static_cast<int>(option)];
}


void InterSpec::setComptonPeakAngle( const int angle )
{
  m_spectrum->setComptonPeakAngle( angle );
}//void setComptonPeakAngle( const float angle );

void InterSpec::toggleFeatureMarkerWindow()
{
  if( m_featureMarkers )
  {
    deleteFeatureMarkerWindow();
    return;
  }//if( m_featureMarkers )

  m_featureMarkers = new FeatureMarkerWindow( this );
  m_featureMarkers->finished().connect( this, &InterSpec::deleteFeatureMarkerWindow );
  
  m_featureMarkerMenuItem->setText( "Hide Feature Markers" );
}//void toggleFeatureMarkerWindow()


void InterSpec::deleteFeatureMarkerWindow()
{
  if( !m_featureMarkers )
    return;
  
  AuxWindow::deleteAuxWindow( m_featureMarkers );
  m_featureMarkers = nullptr;
  m_featureMarkerMenuItem->setText( "Feature Markers..." );
  
  for( FeatureMarkerType i = FeatureMarkerType(0);
       i < FeatureMarkerType::NumFeatureMarkers;
       i = FeatureMarkerType(static_cast<int>(i)+1) )
  {
    if( m_featureMarkersShown[static_cast<int>(i)] )
      setFeatureMarkerOption( i, false );
  }
}//void deleteFeatureMarkerWindow()


Wt::Signal<SpecUtils::SpectrumType,std::shared_ptr<SpecMeas>, std::set<int>, vector<string> > &
                                      InterSpec::displayedSpectrumChanged()
{
  return m_displayedSpectrumChangedSignal;
}


Wt::Signal<SpecUtils::SpectrumType,double> &InterSpec::spectrumScaleFactorChanged()
{
  return m_spectrumScaleFactorChanged;
}


WModelIndex InterSpec::addPeak( PeakDef peak,
                                    const bool associateShowingNuclideXrayRctn )
{
  if( fabs(peak.mean())<0.1 && fabs(peak.amplitude())<0.1 )
    return WModelIndex();
  
  if( !m_referencePhotopeakLines || !associateShowingNuclideXrayRctn )
    return m_peakModel->addNewPeak( peak );
  
  if( peak.parentNuclide() || peak.xrayElement() || peak.reaction() )
    return m_peakModel->addNewPeak( peak );
  
  const bool showingEscape = showingFeatureMarker(FeatureMarkerType::EscapePeakMarker);
  auto foreground = displayedHistogram(SpecUtils::SpectrumType::Foreground);
  PeakSearchGuiUtils::assign_nuclide_from_reference_lines( peak, m_peakModel,
                         foreground, m_referencePhotopeakLines,
                         m_colorPeaksBasedOnReferenceLines, showingEscape );
  
  WModelIndex newpeakindex = m_peakModel->addNewPeak( peak );
  
  PeakModel::PeakShrdPtr newpeak = m_peakModel->peak(newpeakindex);
  try_update_hint_peak( newpeak, m_dataMeasurement, m_displayedSamples );
  
  return newpeakindex;
}//WModelIndex addPeak( PeakDef peak )



#if( IOS )
void InterSpec::exportSpecFile()
{
  //Need to:
  //  -Build gui to allow selecting foreground/background/secondary as well
  //   selecting which samples/detectors/etc is wanted; this can be made general
  //   for all operating systems or deployments.
  //  -Enable/disable "Export Spectrum" menu item as a spectrum is avaialble or not

  if( !m_dataMeasurement )
    return;
  
  hand_spectrum_to_other_app( m_dataMeasurement );
}//void exportSpecFile()
#endif



#if( USE_DB_TO_STORE_SPECTRA )
void InterSpec::saveStateToDb( Wt::Dbo::ptr<UserState> entry )
{
  if( !entry || entry->user != m_user || !m_user.session() )
    throw runtime_error( "Invalid input to saveStateToDb()" );
  
  if( entry->isWriteProtected() )
    throw runtime_error( "UserState is write protected; "
                         "please clone state to re-save" );
  
//  if( !m_user->preferenceValue<bool>("SaveSpectraToDb") )
//    throw runtime_error( "You must enable saving spectra to database inorder to"
//                         " save the apps state." );

  try
  {
    saveShieldingSourceModelToForegroundSpecMeas();
    
    DataBaseUtils::DbTransaction transaction( *m_sql );
    entry.modify()->serializeTime = WDateTime::currentDateTime();
  
    if( !entry->creationTime.isValid() )
      entry.modify()->creationTime = entry->serializeTime;
    
    const bool forTesting = (entry->stateType == UserState::kForTest);
    
    std::shared_ptr<const SpecMeas> foreground = measurment( SpecUtils::SpectrumType::Foreground );
    std::shared_ptr<const SpecMeas> second = measurment( SpecUtils::SpectrumType::SecondForeground );
    std::shared_ptr<const SpecMeas> background = measurment( SpecUtils::SpectrumType::Background );
    
    Dbo::ptr<UserFileInDb> dbforeground = measurmentFromDb( SpecUtils::SpectrumType::Foreground, true );
    Dbo::ptr<UserFileInDb> dbsecond     = measurmentFromDb( SpecUtils::SpectrumType::SecondForeground, true );
    Dbo::ptr<UserFileInDb> dbbackground = measurmentFromDb( SpecUtils::SpectrumType::Background, true );
  
    //JIC, make sure indices have all been assigned to everything.
    m_sql->session()->flush();
    
    if( (entry->snapshotTagParent || forTesting)
        || (dbforeground && dbforeground->isWriteProtected()) )
    {
      dbforeground = UserFileInDb::makeDeepWriteProtectedCopyInDatabase( dbforeground, *m_sql, true );
      if( dbforeground && !(entry->snapshotTagParent || forTesting) )
        UserFileInDb::removeWriteProtection( dbforeground );
    }
    
    if( (entry->snapshotTagParent || forTesting)
         && dbsecond
         && !dbsecond->isWriteProtected()
         && (foreground != second) )
    {
      dbsecond = UserFileInDb::makeDeepWriteProtectedCopyInDatabase( dbsecond, *m_sql, true );
      if( dbsecond && !(entry->snapshotTagParent || forTesting) )
        UserFileInDb::removeWriteProtection( dbsecond );
    }else if( foreground == second )
      dbsecond = dbforeground;
    
    if( (entry->snapshotTagParent || forTesting)
        && dbbackground
        && !dbbackground->isWriteProtected()
        && (background != foreground)
        && (background != second))
    {
      dbbackground = UserFileInDb::makeDeepWriteProtectedCopyInDatabase( dbbackground, *m_sql, true );
      if( dbbackground && !(entry->snapshotTagParent || forTesting) )
        UserFileInDb::removeWriteProtection( dbbackground );
    }else if( background == foreground )
      dbbackground = dbforeground;
    else if( background == second )
      dbbackground = dbsecond;
    
    if( !dbforeground && foreground )
      throw runtime_error( "Error saving foreground to the database" );
    if( !dbsecond && second )
      throw runtime_error( "Error saving second foreground to the database" );
    if( !dbbackground && background )
      throw runtime_error( "Error saving background to the database" );
    
    entry.modify()->snapshotTagParent = entry->snapshotTagParent;
    entry.modify()->snapshotTags = entry->snapshotTags;
    
    //We need to make sure dbforeground, dbbackground, dbsecond will have been
    //  written to the database, so there id()'s will be not -1.
    m_sql->session()->flush();
    
    entry.modify()->foregroundId = dbforeground.id();
    entry.modify()->backgroundId = dbbackground.id();
    entry.modify()->secondForegroundId = dbsecond.id();
    
    /*
    entry.modify()->shieldSourceModelId = -1;
    if( m_shieldingSourceFit )
    {
      Wt::Dbo::ptr<ShieldingSourceModel> model = m_shieldingSourceFit->modelInDb();
      
      //make sure model has actually been written to the DB, so index wont be 0
      if( model )
        m_sql->session()->flush();
        
      entry.modify()->shieldSourceModelId = model.id();
    }
     */
    
//    entry.modify()->otherSpectraCsvIds;
    {
      string &str = (entry.modify()->foregroundSampleNumsCsvIds);
      str.clear();
      for( const int i : m_displayedSamples )
        str += (str.empty() ? "" : ",") + std::to_string( i );
    }
    
    {
      string &str = (entry.modify()->backgroundSampleNumsCsvIds);
      str.clear();
      for( const int i : m_backgroundSampleNumbers )
        str += (str.empty() ? "" : ",") + std::to_string( i );
    }
    
    {
      string &str = (entry.modify()->secondForegroundSampleNumsCsvIds);
      str.clear();
      for( const int i : m_sectondForgroundSampleNumbers )
        str += (str.empty() ? "" : ",") + std::to_string( i );
    }
    
    {
      const vector<string> dispDet = detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
      const vector<string> det_names = foreground ? foreground->detector_names() : vector<string>();
      const vector<int>  detNums = foreground ? foreground->detector_numbers() : vector<int>();
      
      string &str = (entry.modify()->showingDetectorNumbersCsv);
      str.clear();
      for( const auto &det : dispDet )
      {
        auto pos = std::find( begin(det_names), end(det_names), det );
        const auto index = pos - begin(det_names);
        if( index < static_cast<int>(detNums.size()) )
          str += (str.empty() ? "" : ",") + to_string(detNums[index]);
      }
    }
    
    entry.modify()->energyAxisMinimum = m_spectrum->xAxisMinimum();
    entry.modify()->energyAxisMaximum = m_spectrum->xAxisMaximum();
    entry.modify()->countsAxisMinimum = m_spectrum->yAxisMinimum();
    entry.modify()->countsAxisMaximum = m_spectrum->yAxisMaximum();
    entry.modify()->displayBinFactor = 0;
    
    entry.modify()->shownDisplayFeatures = 0x0;
    if( toolTabsVisible() )
      entry.modify()->shownDisplayFeatures |= UserState::kDockedWindows;
    if( m_spectrum->yAxisIsLog() )
      entry.modify()->shownDisplayFeatures |= UserState::kLogSpectrumCounts;
    
    try
    {
      if( m_user->preferenceValue<bool>( "ShowVerticalGridlines" ) )
        entry.modify()->shownDisplayFeatures |= UserState::kVerticalGridLines;
    }catch(...)
    {
      // We can get here if we loaded from a state that didnt have this preference
    }
    
    try
    {
      if( m_user->preferenceValue<bool>( "ShowHorizontalGridlines" ) )
        entry.modify()->shownDisplayFeatures |= UserState::kHorizontalGridLines;
    }catch(...)
    {
      // We can get here if we loaded from a state that didnt have this preference
    }
    
    if( m_spectrum->legendIsEnabled() )
    {
      cerr << "Legend enabled" << endl;
      entry.modify()->shownDisplayFeatures |= UserState::kSpectrumLegend;
    }else cerr << "Legend NOT enabled" << endl;
    
    if( m_shieldingSourceFit )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingShieldSourceFit;
    
    entry.modify()->backgroundSubMode = UserState::kNoSpectrumSubtract;
    if( m_spectrum->backgroundSubtract() )
      entry.modify()->backgroundSubMode = UserState::kBackgorundSubtract;
    
    entry.modify()->currentTab = UserState::kNoTabs;
    if( m_toolsTabs )
    {
      const WString &txt = m_toolsTabs->tabText( m_toolsTabs->currentIndex() );
      if( txt == PeakInfoTabTitle )
        entry.modify()->currentTab = UserState::kPeakInfo;
      else if( txt == GammaLinesTabTitle )
        entry.modify()->currentTab = UserState::kGammaLines;
      else if( txt == CalibrationTabTitle )
        entry.modify()->currentTab = UserState::kCalibration;
      else if( txt == NuclideSearchTabTitle )
        entry.modify()->currentTab = UserState::kIsotopeSearch;
      else if( txt == FileTabTitle )
        entry.modify()->currentTab = UserState::kFileTab;
    }//if( m_toolsTabs )
    
    entry.modify()->showingMarkers = 0x0;
    
    entry.modify()->disabledNotifications = 0x0;

    for( WarningWidget::WarningMsgLevel level = WarningWidget::WarningMsgLevel(0);
        level <= WarningWidget::WarningMsgHigh;
        level = WarningWidget::WarningMsgLevel(level+1) )
    {
      if( !m_warnings->active( level ) )
        entry.modify()->disabledNotifications |= (0x1<<level);
    }
    
    entry.modify()->showingPeakLabels = 0x0;
    for( SpectrumChart::PeakLabels label = SpectrumChart::PeakLabels(0);
         label < SpectrumChart::kNumPeakLabels;
         label = SpectrumChart::PeakLabels(label+1) )
    {
      if( m_spectrum->showingPeakLabel( label ) )
        entry.modify()->showingPeakLabels |= (0x1<<label);
    }

    entry.modify()->showingWindows = 0x0;
    //m_warningsWindow
    //m_peakInfoWindow
    //m_energyCalWindow
    //m_shieldingSourceFitWindow
    //m_referencePhotopeakLinesWindow
//    m_nuclideSearchWindow
//    enum ShowingWindows
//    {
//      kEnergyCalibration = 0x1, DrfSelectSelect = 0x2 //etc..
//    };
  
    entry.modify()->isotopeSearchEnergiesXml.clear();
    if( m_nuclideSearch )
      m_nuclideSearch->serialize( entry.modify()->isotopeSearchEnergiesXml );
    
    if( !toolTabsVisible() && m_nuclideSearchWindow )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingEnergySearch;
      entry.modify()->currentTab = UserState::kIsotopeSearch;
    }
  
   
    entry.modify()->gammaLinesXml.clear();
    if( m_referencePhotopeakLines )
      m_referencePhotopeakLines->serialize( entry.modify()->gammaLinesXml );
    
  /*
   enum FeatureMarkers
   {
   kEscapePeaks = 0x1, kCompPeak = 0x2, kComptonEdge = 0x4, kSumPeak = 0x8
   };
   */

    Json::Array userOptions;
    const Wt::Dbo::collection< Wt::Dbo::ptr<UserOption> > &prefs
                                                     = m_user->m_dbPreferences;
    
    vector< Dbo::ptr<UserOption> > options;
    std::copy( prefs.begin(), prefs.end(), std::back_inserter(options) );
    
    for( vector< Dbo::ptr<UserOption> >::const_iterator iter = options.begin();
        iter != options.end(); ++iter )
    {
      Dbo::ptr<UserOption> option = *iter;
      Json::Value val( Json::ObjectType );
      Json::Object &obj = val;
      obj["type"]  = int(option->m_type);
      obj["name"]  = WString(option->m_name);
      
      switch( option->m_type )
      {
        case UserOption::String:
          obj["value"] = WString( option->m_value );
        break;
        
        case UserOption::Decimal:
          obj["value"] = std::stod( option->m_value );
          
        break;
        
        case UserOption::Integer:
          obj["value"] = std::stoi( option->m_value );
        break;
        
        case UserOption::Boolean:
          obj["value"] = (option->m_value== "true" || option->m_value=="1");
        break;
      }//switch( m_type )
      
      userOptions.push_back( val );
    }//for( loop over DB entries )
  
    entry.modify()->userOptionsJson = Json::serialize( userOptions, 0 );
    
    //The color theme, if anything other than default is being used.
    //  (I guess we could serialize the default, but whateer for now)
    try
    {
      Wt::Dbo::ptr<ColorThemeInfo> theme;
      const int colorThemIndex = m_user->preferenceValue<int>("ColorThemeIndex", this);
      if( colorThemIndex > 0 )
        theme = m_user->m_colorThemes.find().where( "id = ?" ).bind( colorThemIndex ).resultValue();
      if( theme )
        entry.modify()->colorThemeJson = theme->json_data;
    }catch(...)
    {
      //Probably shouldnt ever happen
      cerr << "saveStateToDb(): Caught exception retrieving color theme from database" << endl;
    }//try / catch
    
    if( forTesting || entry->snapshotTagParent )
      UserState::makeWriteProtected( entry, m_sql->session() );
    
    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "saveStateToDb(...) caught: " << e.what() << endl;
    throw runtime_error( "I'm sorry there was a technical error saving your"
                         " user state to the database." );
    //to do: handle errors
  }//try / catch
}//void saveStateToDb( Wt::Dbo::ptr<UserState> entry )


Dbo::ptr<UserState> InterSpec::serializeStateToDb( const Wt::WString &name,
                                                        const Wt::WString &desc,
                                                        const bool forTesting,
                                                        Dbo::ptr<UserState> parent)
{
  Dbo::ptr<UserState> answer;
  
  UserState *state = new UserState();
  state->user = m_user;
  state->stateType = forTesting ? UserState::kForTest : UserState::kUserSaved;
  state->creationTime = WDateTime::currentDateTime();
  state->name = name;
  state->description = desc;

  
  {//Begin interaction with database
    DataBaseUtils::DbTransaction transaction( *m_sql );
    if( parent )
    {
      //Making sure we get the main parent
      while( parent->snapshotTagParent )
        parent = parent->snapshotTagParent;
      state->snapshotTagParent = parent;
    }//parent
    
    answer = m_user.session()->add( state );
    transaction.commit();
  }//End interaction with database
  
  try
  {
    saveStateToDb( answer );
  }catch( std::exception &e )
  {
    passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
  }
      
  if (parent) {
    WString msg = "Created new tag '";
    msg += name.toUTF8();
    msg += "' under snapshot '";
    msg += parent.get()->name.toUTF8();
    msg += "'.";
    passMessage( msg, "", WarningWidget::WarningMsgInfo );
  }else
  {
    WString msg = "Created new snapshot '";
    msg += name.toUTF8();
    msg += "'";
    passMessage( msg, "", WarningWidget::WarningMsgSave );
  }
    
  return answer;
}//Dbo::ptr<UserState> serializeStateToDb(...)



void InterSpec::loadStateFromDb( Wt::Dbo::ptr<UserState> entry )
{
  //This function takes about 5 ms (to load the previous HPGe state, on new app
  //  construction) on my 2.6 GHz Intel Core i7, as of (20150316).
  if( !entry )
    throw runtime_error( "UserState was invalid" );

  Wt::Dbo::ptr<UserState> parent = entry->snapshotTagParent;
  
  try
  {
    closeShieldingSourceFitWindow();
    
    switch( entry->stateType )
    {
      case UserState::kEndOfSessionTemp:
      case UserState::kUserSaved:
      case UserState::kForTest:
      case UserState::kUndefinedStateType:
      case UserState::kEndOfSessionHEAD:
      break;
    }//switch( entry->stateType )

    DataBaseUtils::DbTransaction transaction( *m_sql );
    
    if( parent )
    {
      while( parent->snapshotTagParent )
        parent = parent->snapshotTagParent;
    }

    
    Dbo::ptr<UserFileInDb> dbforeground, dbsecond, dbbackground;
    if( entry->foregroundId >= 0 )
      dbforeground = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                      .bind( entry->foregroundId );
    if( entry->backgroundId >= 0
        && entry->backgroundId != entry->foregroundId )
      dbbackground = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                      .bind( entry->backgroundId );
    else if( entry->backgroundId == entry->foregroundId )
      dbbackground = dbforeground;
      
    if( entry->secondForegroundId >= 0
        && entry->secondForegroundId != entry->foregroundId
        && entry->secondForegroundId != entry->backgroundId )
      dbsecond = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                  .bind( entry->secondForegroundId );
    else if( entry->secondForegroundId == entry->foregroundId )
      dbsecond = dbforeground;
    else if( entry->secondForegroundId == entry->backgroundId )
      dbsecond = dbbackground;

    if( entry->foregroundId >= 0
        && (!dbforeground || !dbforeground->filedata.size()) )
      throw runtime_error( "Unable to locate foreground in database" );
    if( entry->backgroundId >= 0
        && (!dbbackground || !dbbackground->filedata.size()) )
      throw runtime_error( "Unable to locate background in database" );
    if( entry->secondForegroundId >= 0
        && (!dbsecond || !dbsecond->filedata.size()) )
      throw runtime_error( "Unable to locate second foreground in database" );
    
    //Essentially reset the state of the app
    m_fileManager->removeAllFiles();
    if( m_referencePhotopeakLines )
      m_referencePhotopeakLines->clearAllLines();
    
    deletePeakEdit();
    deleteGammaCountDialog();
    
    if( m_shieldingSourceFitWindow && !m_shieldingSourceFitWindow->isHidden() )
      m_shieldingSourceFitWindow->hide();
    
    bool wasDocked = (entry->shownDisplayFeatures & UserState::kDockedWindows);
    if( toolTabsVisible() != wasDocked )
      setToolTabsVisible( wasDocked );

    //Now start reloading the state
    std::shared_ptr<SpecMeas> foreground, second, background;
    std::shared_ptr<SpecMeas> snapforeground, snapsecond, snapbackground;
    std::shared_ptr<SpectraFileHeader> foregroundheader, backgroundheader,
                                          secondheader;
    
    
    if( parent )
    {
      //If we are loading a snapshot, we actually want to associate the
      //  foreground/second/backgorund entries in the database with the parent
      //  entries, but actually load the snapshots.  This is so when the user
      //  saves changes, the parents saved state recieves the changes, and the
      //  snapshot doesnt get changed.
      //All of this is really ugly and horrible, and should be improved!
      
//      const bool autosave
//        = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", this );
//      if( autosave )
//      {
//      Its actually to late, HEAD will be set to the tag version even if the
//      user immediate selects save as.
//        const char *txt = "Auto saving is currently enabled, meaning the most"
//                          " recent version of this snapshot will get"
//                          " over-written with the current state (plus any"
//                          " changes you make) unless you disable auto saving,"
//                          " or do a &quot;Save As&quot; asap.";
//        passMessage( txt, "", WarningWidget::WarningMsgMedium );
//      }
      
      if( dbforeground && dbforeground->filedata.size() )
        snapforeground = (*dbforeground->filedata.begin())->decodeSpectrum();
    
      if( dbsecond && dbsecond->filedata.size() && dbsecond != dbforeground )
        snapsecond = (*dbsecond->filedata.begin())->decodeSpectrum();
      else if( dbsecond == dbforeground )
        snapsecond = snapforeground;
    
      if( dbbackground && dbbackground->filedata.size()
          && dbbackground != dbforeground
          && dbbackground != dbsecond  )
        snapbackground = (*dbbackground->filedata.begin())->decodeSpectrum();
      else if( dbbackground == dbforeground )
        snapbackground = snapforeground;
      else if( dbbackground == dbsecond )
        snapbackground = snapsecond;
      
      
      if( parent->foregroundId >= 0 )
        dbforeground = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                        .bind( parent->foregroundId );
      if( parent->backgroundId >= 0
          && parent->backgroundId != parent->foregroundId )
      {
        dbbackground = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                        .bind( parent->backgroundId );
      }else if( parent->backgroundId == parent->foregroundId )
      {
        dbbackground = dbforeground;
      }else if( parent->backgroundId < 0 && entry->backgroundId >= 0 )
      {
        parent.modify()->backgroundId = entry->backgroundId;
      }else
      {
        dbbackground.reset();
      }
      
      if( parent->secondForegroundId >= 0
         && parent->secondForegroundId != parent->foregroundId
         && parent->secondForegroundId != parent->backgroundId )
      {
        dbsecond = m_sql->session()->find<UserFileInDb>().where( "id = ?" )
                                    .bind( parent->secondForegroundId );
      }else if( parent->secondForegroundId < 0
                && entry->secondForegroundId >= 0 )
      {
        parent.modify()->secondForegroundId = entry->secondForegroundId;
      }else if( parent->secondForegroundId == parent->foregroundId )
      {
        dbsecond = dbforeground;
      }else if( parent->secondForegroundId == parent->backgroundId )
      {
        dbsecond = dbbackground;
      }else
      {
        dbsecond.reset();
      }
    }//if( parent )
    
    if( dbforeground )
      m_fileManager->setDbEntry( dbforeground, foregroundheader, foreground, 0);
    
    if( !parent || snapbackground )
    {
      if( dbbackground && dbbackground != dbforeground )
        m_fileManager->setDbEntry( dbbackground, backgroundheader, background, 0);
      else if( dbbackground == dbforeground )
        background = foreground;
    }//if( !parent || snapbackground )
    
    if( !parent || snapsecond )
    {
      if( dbsecond && dbsecond != dbforeground && dbsecond != dbbackground )
        m_fileManager->setDbEntry( dbsecond, secondheader, second, 0);
      else if( dbsecond == dbforeground )
        second = foreground;
      else if( dbsecond == dbbackground )
        second = background;
    }//if( !parent || snapsecond )
    
    if( parent )
    {
      if( foregroundheader && snapforeground )
      {
        foregroundheader->setMeasurmentInfo( snapforeground );
        foreground = snapforeground;
      }
      
      if( snapbackground != snapforeground )
      {
        if( backgroundheader && snapbackground )
        {
          backgroundheader->setMeasurmentInfo( snapbackground );
          background = snapbackground;
        }else if( snapbackground && dbbackground )
        {
          dbbackground = UserFileInDb::makeDeepWriteProtectedCopyInDatabase( dbbackground, *m_sql, true );
          if( dbbackground )
            UserFileInDb::removeWriteProtection( dbbackground );
          m_fileManager->setDbEntry( dbbackground, backgroundheader, background, 0);
        }
      }//if( snapbackground != snapforeground )
      
      
      if( (snapsecond != snapforeground) && (snapsecond != snapbackground) )
      {
        if( secondheader && snapsecond )
        {
          secondheader->setMeasurmentInfo( snapsecond );
          second = snapsecond;
        }else if( snapsecond && dbsecond )
        {
          dbsecond = UserFileInDb::makeDeepWriteProtectedCopyInDatabase( dbsecond, *m_sql, true );
          if( dbsecond )
            UserFileInDb::removeWriteProtection( dbsecond );
          m_fileManager->setDbEntry( dbsecond, secondheader, second, 0);
        }
      }
    }//if( parent )
    
    auto csvToInts = []( const string &str ) -> set<int> {
      set<int> samples;
      stringstream strm( str );
      std::copy( istream_iterator<int>( strm ), istream_iterator<int>(),
                std::inserter( samples, end(samples) ) );
      return samples;
    }; //csvToInts lambda
    
    const set<int> foregroundNums = csvToInts( entry->foregroundSampleNumsCsvIds );
    const set<int> secondNums     = csvToInts( entry->secondForegroundSampleNumsCsvIds );
    const set<int> backgroundNums = csvToInts( entry->backgroundSampleNumsCsvIds );
    const set<int> otherSamples   = csvToInts( entry->otherSpectraCsvIds );
    
    
    setSpectrum( foreground, foregroundNums, SpecUtils::SpectrumType::Foreground, 0 );
    if( foreground )
    {
      //If we dont have a foreground, we probably shouldnt be loading the state, but...
      setSpectrum( background, backgroundNums, SpecUtils::SpectrumType::Background, 0 );
      setSpectrum( second, secondNums, SpecUtils::SpectrumType::SecondForeground, 0 );
    }
    
    
    //Load the other spectra the user had opened.  Note that they were not
    //  write protected so they may have been changed or removed
    //...should check to make sure that they are lazily loaded, so as to not
    //  waste to much CPU.
    for( const int id : otherSamples )
    {
      try
      {
        Dbo::ptr<UserFileInDb> dbfile;
        std::shared_ptr<SpecMeas> meas;
        std::shared_ptr<SpectraFileHeader> header;
        dbfile = m_sql->session()->find<UserFileInDb>().where( "id = ?" ).bind( id );
        if( dbfile )
          m_fileManager->setDbEntry( dbfile, header, meas, 0 );
      }catch(...){}
    }//for( const int id : otherSamples )
    
    //Load the ShieldingSourceModel
    Dbo::ptr<ShieldingSourceModel> fitmodel;
    if( entry->shieldSourceModelId >= 0 )
      fitmodel = m_sql->session()->find<ShieldingSourceModel>()
                         .where( "id = ?" ).bind( entry->shieldSourceModelId );
    if( fitmodel )
    {
      cerr << "When loading state from DB, not loading Shielding/Source Fit"
      << " model, due to lame GUI ish" << endl;
//      if( !m_shieldingSourceFit )
//        showShieldingSourceFitWindow();
      if( m_shieldingSourceFit )
        m_shieldingSourceFit->loadModelFromDb( fitmodel );
    }//if( fitmodel )
    
    if( m_shieldingSourceFitWindow )
      m_shieldingSourceFitWindow->hide();
    
    if( m_nuclideSearchWindow )
    {
      m_nuclideSearchContainer->layout()->removeWidget( m_nuclideSearch );
      delete m_nuclideSearchContainer;
      m_nuclideSearchContainer = nullptr;
    }
    
//  std::string showingDetectorNumbersCsv;
//  see:  std::vector<bool> detectorsToDisplay() const;
    
    float xmin = entry->energyAxisMinimum;
    float xmax = entry->energyAxisMaximum;
    std::shared_ptr<const SpecUtils::Measurement> hist = m_spectrum->data();
    
    if( hist && (xmin != xmax) && (xmax > xmin)
        && !IsInf(xmin) && !IsInf(xmax) && !IsNan(xmin) && !IsNan(xmax) )
    {
      const size_t nbin = hist->num_gamma_channels();
      xmin = std::max( xmin, hist->gamma_channel_lower( 0 ) );
      xmax = std::min( xmax, hist->gamma_channel_upper( nbin - 1 ) );
      m_spectrum->setXAxisRange( xmin, xmax );
    }//if( xmin != xmax )
        
    
    if( foreground )
    {
      // If we have a state without a foreground, then for the spectrum chart, changing the legend
      //  enabled status, or setting the horizontal/vertical lines causes the client-side javascript
      //  to crash - I spent some hours trying to figure it out, but no luck.  So for the moment,
      //  we'll just do this work-around.
      //  One hint was for the crashing case D3SpectrumDisplayDiv::defineJavaScript() seems to be
      //  called multiple times, and maybe the id of the parent div changes???
      //  Its really bizarre.
#ifdef _MSC_VER
#pragma message("Not setting show legend or grid lines on charts when loading states that do not have a foreground; should figure out why this causes a client-side error")
#else
#warning "Not setting show legend or grid lines on charts when loading states that do not have a foreground; should figure out why this causes a client-side error"
#endif
      
      const bool vertGridLines = (entry->shownDisplayFeatures & UserState::kVerticalGridLines);
      const bool horizontalGridLines = (entry->shownDisplayFeatures & UserState::kVerticalGridLines);
      m_spectrum->showVerticalLines( vertGridLines );
      m_spectrum->showHorizontalLines( horizontalGridLines );
      m_timeSeries->showVerticalLines( vertGridLines );
      m_timeSeries->showHorizontalLines( horizontalGridLines );
      
      if( (entry->shownDisplayFeatures & UserState::kSpectrumLegend) )
      {
        m_spectrum->enableLegend();
      }else
      {
        m_spectrum->disableLegend();
      }
    }//if( foreground )
    
//  SpectrumSubtractMode backgroundSubMode;
    
    if( wasDocked )
    {
      WString title;
      switch( entry->currentTab )
      {
        case UserState::kPeakInfo:      title = PeakInfoTabTitle;      break;
        case UserState::kGammaLines:    title = GammaLinesTabTitle;    break;
        case UserState::kCalibration:   title = CalibrationTabTitle;   break;
        case UserState::kIsotopeSearch: title = NuclideSearchTabTitle; break;
        case UserState::kFileTab:       title = FileTabTitle;          break;
        case UserState::kNoTabs:                                       break;
      };//switch( entry->currentTab )
      
      for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
      {
        if( m_toolsTabs->tabText( tab ) == title )
        {
          // Note: this causes handleToolTabChanged(...) to be called
          m_toolsTabs->setCurrentIndex( tab );
          m_currentToolsTab = tab;
          break;
        }//if( m_toolsTabs->tabText( tab ) == title )
      }//for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
    }//if( wasDocked )
    
    
    //Should take care of entry->showingWindows here, so we can relie on it below
    if( entry->gammaLinesXml.size() )
    {
      bool toolTabsHidden = false;
      if( !m_referencePhotopeakLines )    //tool tabs must be hidden, so lets show the window
      {
        showGammaLinesWindow();//if entry->showingWindows was implemented
        toolTabsHidden = true;
        assert( m_referencePhotopeakLines );
      }
      
      string data = entry->gammaLinesXml;
      m_referencePhotopeakLines->deSerialize( data );
      
      if( toolTabsHidden )
        closeGammaLinesWindow();
    }//if( entry->gammaLinesXml.size() )
    
    if( m_nuclideSearch && entry->isotopeSearchEnergiesXml.size() )
    {
      string data = entry->isotopeSearchEnergiesXml;
      
      bool display = (entry->currentTab == UserState::kIsotopeSearch);
      
      if( !wasDocked )
      {
        display = (entry->shownDisplayFeatures & UserState::kShowingEnergySearch);
        if( display )
          showNuclideSearchWindow();
      }//if( !wasDocked )
      
      m_nuclideSearch->deSerialize( data, display );
    }//if( m_nuclideSearch && entry->isotopeSearchEnergiesXml.size() )

    for( SpectrumChart::PeakLabels label = SpectrumChart::PeakLabels(0);
        label < SpectrumChart::kNumPeakLabels;
        label = SpectrumChart::PeakLabels(label+1) )
    {
      const bool showing = (entry->showingPeakLabels & (0x1<<label));
      m_spectrum->setShowPeakLabel( label, showing );
    }//for( loop over peak labels )

    
    if( entry->userOptionsJson.size() )
    {
      try
      {
      Json::Value userOptionsVal( Json::ArrayType );
      Json::parse( entry->userOptionsJson, userOptionsVal );
      //cerr << "Parsed User Options, but not doing anything with them" << endl;
        
      const Json::Array &userOptions = userOptionsVal;
      for( size_t i = 0; i < userOptions.size(); ++i )
      {
          const Json::Object &obj = userOptions[i];
          WString name = obj.get("name");
          const int type = obj.get("type");
          const UserOption::DataType datatype = UserOption::DataType( type );
      
          switch( datatype )
          {
            case UserOption::String:
            {
              std::string value = obj.get("value");
              InterSpecUser::pushPreferenceValue( m_user, name.toUTF8(),
                                                 boost::any(value), this, wApp );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Decimal:
            {
              double value = obj.get("value");
              InterSpecUser::pushPreferenceValue( m_user, name.toUTF8(),
                                                 boost::any(value), this, wApp );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Integer:
            {
              int value = obj.get("value");
              InterSpecUser::pushPreferenceValue( m_user, name.toUTF8(),
                                                 boost::any(value), this, wApp );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Boolean:
            {
              bool value = obj.get("value");
              InterSpecUser::pushPreferenceValue( m_user, name.toUTF8(),
                                                  boost::any(value), this, wApp );
            }//case UserOption::Boolean:
          }//switch( datatype )
        }//for( size_t i = 0; i < userOptions.size(); ++i )
      }catch( std::exception &e )
      {
        cerr << "Failed to parse JSON user preference with error: "
             << e.what() << endl;
        //well let this error slide (happens on OSX, Wt::JSON claims boost isnt
        //  new enough)
      }//try / catch
    }//if( entry->userOptionsJson.size() )
    
    if( (entry->shownDisplayFeatures & UserState::kShowingShieldSourceFit) )
      showShieldingSourceFitWindow();
    
    if( entry->colorThemeJson.size() > 32 )  //32 is arbitrary
    {
      try
      {
        auto theme = make_shared<ColorTheme>();
        ColorTheme::fromJson( entry->colorThemeJson, *theme );
        applyColorTheme( theme );
      }catch(...)
      {
        cerr << "loadStateFromDb(): Caught exception to set hte color them" << endl;
      }//try / catch to load ColorTheme.
    }//if( its reasonable we might have color theme information )
    
//  int showingMarkers;        //bitwise or of FeatureMarkers
//  int showingWindows;
    
    m_currentStateID = entry.id();
    if( parent )
      m_currentStateID = parent.id();
    
    transaction.commit();
  }catch( std::exception &e )
  {
    m_currentStateID = -1;
    updateSaveWorkspaceMenu();
    cerr << "\n\n\n\nloadStateFromDb(...) caught: " << e.what() << endl;
    throw runtime_error( e.what() );
  }//try / catch
  
  m_currentStateID = entry.id();
  updateSaveWorkspaceMenu();
}//void loadStateFromDb( Wt::Dbo::ptr<UserState> entry )


//Whenever m_currentStateID changes value, we need to update the menu!
void InterSpec::updateSaveWorkspaceMenu()
{
  m_saveStateAs->setDisabled( !m_dataMeasurement );
  
  //check if Save menu needs to be updated
  if( m_currentStateID >= 0 )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      Dbo::ptr<UserState> state = m_sql->session()
                                       ->find<UserState>( "where id = ?" )
                                       .bind( m_currentStateID );
      transaction.commit();
      
      if( state && (state->stateType==UserState::kEndOfSessionHEAD
                    || state->stateType==UserState::kUserSaved) )
      {
        m_saveState->setDisabled(false);
        m_createTag->setDisabled(false);
      } //UserState exists
      else
      { //no UserState
        m_saveState->setDisabled(true);
        m_createTag->setDisabled(true);
      } //no UserState
      return;
    }catch( std::exception &e )
    {
      cerr << "\nstartStoreStateInDb() caught: " << e.what() << endl;
    }//try / catch
  }

  m_saveState->setDisabled(true);
  m_createTag->setDisabled(true);
}//void updateSaveWorkspaceMenu()


int InterSpec::currentAppStateDbId()
{
  return m_currentStateID;
}//int currentAppStateDbId()


void InterSpec::resetCurrentAppStateDbId()
{
  m_currentStateID = -1;
  updateSaveWorkspaceMenu();
}//void resetCurrentAppStateDbId()
//m_currentStateID

#endif //#if( USE_DB_TO_STORE_SPECTRA )


std::shared_ptr<const ColorTheme> InterSpec::getColorTheme()
{
  if( m_colorTheme )
    return m_colorTheme;
  
  auto theme = make_shared<ColorTheme>();
  
  try
  {
    const int colorThemIndex = m_user->preferenceValue<int>("ColorThemeIndex", this);
    if( colorThemIndex < 0 )
    {
      auto deftheme = ColorTheme::predefinedTheme( ColorTheme::PredefinedColorTheme(-colorThemIndex) );
      if( deftheme )
        theme.reset( deftheme.release() );
    }else
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      auto themeFromDb = m_user->m_colorThemes.find().where( "id = ?" ).bind( colorThemIndex ).resultValue();
      transaction.commit();
   
      theme->setFromDataBase( *themeFromDb );
    }
  }catch(...)
  {
    //We'll get the default color theme here.
  }//try / catch
  
  m_colorTheme = theme;
  
  return theme;
}//Wt::Dbo<ColorTheme> getColorTheme()


void InterSpec::setReferenceLineColors( std::shared_ptr<const ColorTheme> theme )
{
  if( !m_referencePhotopeakLines )
    return;
  
  if( !theme )
    theme = getColorTheme();
  
  m_referencePhotopeakLines->setPeaksGetAssignedRefLineColor( theme->peaksTakeOnReferenceLineColor );
  m_referencePhotopeakLines->setColors( theme->referenceLineColor );
  m_referencePhotopeakLines->setColorsForSpecificSources( theme->referenceLineColorForSources );
}//void setReferenceLineColors( Wt::Dbo<ColorTheme> ptr )


void InterSpec::applyColorTheme( shared_ptr<const ColorTheme> theme )
{
  //20190324: This function, called with a nullptr, takes 30us on my 2018 mbp.
  if( !theme )
    theme = getColorTheme();
  m_colorTheme = theme;

  m_colorPeaksBasedOnReferenceLines = theme->peaksTakeOnReferenceLineColor;
  
  m_spectrum->setForegroundSpectrumColor( theme->foregroundLine );
  m_spectrum->setBackgroundSpectrumColor( theme->backgroundLine );
  m_spectrum->setSecondarySpectrumColor( theme->secondaryLine );
  m_spectrum->setDefaultPeakColor( theme->defaultPeakLine );
  
  m_spectrum->setAxisLineColor( theme->spectrumAxisLines );
  m_spectrum->setChartMarginColor( theme->spectrumChartMargins );
  m_spectrum->setChartBackgroundColor( theme->spectrumChartBackground );
  m_spectrum->setTextColor( theme->spectrumChartText );
  
  m_timeSeries->applyColorTheme( theme );
  
  setReferenceLineColors( theme );
  
  string cssfile;
  if( !theme->nonChartAreaTheme.empty() )
    cssfile = "InterSpec_resources/themes/" + theme->nonChartAreaTheme + "/" + theme->nonChartAreaTheme + ".css";
  
  
  if( !m_currentColorThemeCssFile.empty() && (m_currentColorThemeCssFile != cssfile) )
  {
    wApp->removeStyleSheet( m_currentColorThemeCssFile );
    m_currentColorThemeCssFile.clear();
  }
  
  if( (m_currentColorThemeCssFile != cssfile)
      && !SpecUtils::iequals_ascii(theme->nonChartAreaTheme, "default") )
  {
    if( SpecUtils::is_file(cssfile) )
    {
      wApp->useStyleSheet( cssfile );
      m_currentColorThemeCssFile = cssfile;
    }else
    {
      passMessage( "Could not load the non-peak-area theme '" + theme->nonChartAreaTheme + "'",
                  "", WarningWidget::WarningMsgLow );
    }
  }//if( should load CSS file )
  
  m_colorThemeChanged.emit( theme );
}//void InterSpec::applyColorTheme()


Wt::Signal< std::shared_ptr<const ColorTheme> > &InterSpec::colorThemeChanged()
{
  return m_colorThemeChanged;
}


#if( BUILD_AS_OSX_APP || APPLY_OS_COLOR_THEME_FROM_JS || IOS || BUILD_AS_ELECTRON_APP )
void InterSpec::osThemeChange( std::string name )
{
  cout << "InterSpec::osThemeChange('" << name << "');" << endl;
  
  SpecUtils::to_lower_ascii(name);
  if( name != "light" && name != "dark" && name != "default" && name != "no-preference" && name != "no-support" )
  {
    cerr << "InterSpec::osThemeChange('" << name << "'): unknown OS theme." << endl;
    return;
  }
  
  try
  {
    const int colorThemIndex = m_user->preferenceValue<int>("ColorThemeIndex", this);
    if( colorThemIndex >= 0 )
      return;
    
    const auto oldTheme = ColorTheme::PredefinedColorTheme(-colorThemIndex);
    
    if( oldTheme != ColorTheme::PredefinedColorTheme::DefaultColorTheme )
      return;
    
    std::unique_ptr<ColorTheme> theme;
    
    if( name == "dark" )
    {
      //Check to see if we already have dark applied
      if( m_colorTheme && SpecUtils::icontains( m_colorTheme->theme_name.toUTF8(), "Dark") )
        return;
     
      cerr << "Will set to dark" << endl;
      theme = ColorTheme::predefinedTheme( ColorTheme::PredefinedColorTheme::DarkColorTheme );
    }else if( name == "light" || name == "default" || name == "no-preference" || name == "no-support" )
    {
      //Check to see if we already have light applied
      if( m_colorTheme && SpecUtils::icontains( m_colorTheme->theme_name.toUTF8(), "Default") )
        return;
      
      cerr << "Will set to default" << endl;
      theme = ColorTheme::predefinedTheme( ColorTheme::PredefinedColorTheme::DefaultColorTheme );
    }
    
    if( theme )
      applyColorTheme( make_shared<ColorTheme>(*theme) );
  }catch(...)
  {
    cerr << "InterSpec::osThemeChange() caught exception - not doing anything" << endl;
  }
  
}//void osThemeChange( std::string name )
#endif


void InterSpec::showColorThemeWindow()
{
  //Takes ~5ms (on the server) to create a ColorThemeWindow.
	new ColorThemeWindow(this);
}//void showColorThemeWindow()


AuxWindow *InterSpec::showIEWarningDialog()
{
  try
  {
    wApp->environment().getCookie( "IEWarningDialog" );
    return NULL;
  }catch(...){}


  AuxWindow *dialog = new AuxWindow( "Not Compatible with Internet Explorer", AuxWindowProperties::IsModal );
  dialog->setAttributeValue( "style", "max-width: 80%;" );

  WContainerWidget *contents = dialog->contents();
  contents->setContentAlignment( Wt::AlignCenter );

  WText *instructions = new WText(
          "This webapp has not been tested with Microsoft Internet Explorer, and<br />"
          " as a result some things will either fail to render or may not even work.<br />"
          "<p>Please consider using a reasonably recent version of Firefox, Chrome, or Safari.</p>"
          , XHTMLUnsafeText, contents );
  instructions->setInline( false );

  WPushButton *ok = new WPushButton( "Close", contents );
  ok->clicked().connect( dialog, &AuxWindow::hide );

  WCheckBox *cb = new WCheckBox( "Do not show again", contents );
  cb->setFloatSide( Right );
  cb->checked().connect( boost::bind( &InterSpec::setShowIEWarningDialogCookie, this, false ) );
  cb->unChecked().connect( boost::bind( &InterSpec::setShowIEWarningDialogCookie, this, true ) );

  dialog->centerWindow();
  dialog->show();
  
  dialog->finished().connect( boost::bind( AuxWindow::deleteAuxWindow, dialog ) );
  
  return dialog;
} // AuxWindow *showIEWarningDialog()


void InterSpec::showLicenseAndDisclaimersWindow( const bool is_awk,
                                      std::function<void()> finished_callback )
{
  if( m_licenseWindow )
  {
    m_licenseWindow->show();
    return;
  }
  
  m_licenseWindow = new LicenseAndDisclaimersWindow( is_awk, renderedWidth(), renderedHeight() );
  
  m_licenseWindow->finished().connect( std::bind([this,finished_callback](){
    deleteLicenseAndDisclaimersWindow();
    
    if( finished_callback )
      finished_callback();
  }) );
}//void showLicenseAndDisclaimersWindow()

void InterSpec::startClearSession()
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( !app )
    return;
  
  SimpleDialog *window = new SimpleDialog( "Clear Session?",
                                          "Are you sure you would like to start a clean session?" );
  
  WPushButton *button = window->addButton( "Yes" );
  button->clicked().connect( app, &InterSpecApp::clearSession );
  
  button = window->addButton( "No" );
  button->setFocus();
}//void startClearSession()


void InterSpec::deleteLicenseAndDisclaimersWindow()
{ 
  if( !m_licenseWindow )
    return;
  
  AuxWindow::deleteAuxWindow( m_licenseWindow );
  m_licenseWindow = nullptr;
}//void deleteLicenseAndDisclaimersWindow()


void InterSpec::showWelcomeDialog( bool force )
{
  try
  {
    if( !force && !m_user->preferenceValue<bool>( "ShowSplashScreen" ) )
      return;
  }catch(...)
  {
    //m_user didnt have preference "ShowSplashScreen" for some reason
    if( !force )
      return;
  }
  
  if( m_useInfoWindow )
  {
    m_useInfoWindow->show();
    return;
  }
  
  WServer::instance()->post( wApp->sessionId(), std::bind( [this](){
    /*
     For Android, showing this useInfoWindow at startup causes some exceptions
     in the JavaScript whenever the loading indicator is shown.  I'm pretty
     confused by it.
     I tried setting a new loading indicator in deleteWelcomeCountDialog(), but it didnt work
     
     Havent checked if creating this window via "posting" helps
     */
    if( m_useInfoWindow )
      return;
    
    std::function<void(bool)> dontShowAgainCallback = [this](bool value){
      InterSpecUser::setPreferenceValue<bool>( m_user, "ShowSplashScreen", value, this );
    };
    m_useInfoWindow = new UseInfoWindow( dontShowAgainCallback , this );
  
    m_useInfoWindow->finished().connect( this, &InterSpec::deleteWelcomeCountDialog );
    
    wApp->triggerUpdate();
  } ) );
}//void showWelcomeDialog()


void InterSpec::deleteWelcomeCountDialog()
{
  if( !m_useInfoWindow )
    return;
  
  AuxWindow::deleteAuxWindow( m_useInfoWindow );
  m_useInfoWindow = nullptr;
}//void deleteWelcomeCountDialog()


void InterSpec::deleteEnergyCalPreserveWindow()
{
  if( m_preserveCalibWindow )
  {
    delete m_preserveCalibWindow;
    m_preserveCalibWindow = nullptr;
  }
}//void deleteEnergyCalPreserveWindow()


void InterSpec::setShowIEWarningDialogCookie( bool show )
{
  if( show )
    wApp->setCookie( "IEWarningDialog", "", 0 );
  else
    wApp->setCookie( "IEWarningDialog", "noshow", 3600*24*30 );
} // void InterSpec::setShowIEWarningDialogCookie( bool show )



void InterSpec::showGammaCountDialog()
{
  if( m_gammaCountDialog )
  {
//    m_gammaCountDialog->show();
//    m_gammaCountDialog->expand();
    return;
  }//if( m_gammaCountDialog )

  m_gammaCountDialog = new GammaCountDialog( this );
//  m_gammaCountDialog->rejectWhenEscapePressed();
  m_gammaCountDialog->finished().connect( this, &InterSpec::deleteGammaCountDialog );
//  m_gammaCountDialog->resizeToFitOnScreen();
}//void showGammaCountDialog()


void InterSpec::deleteGammaCountDialog()
{
  if( m_gammaCountDialog )
    delete m_gammaCountDialog;

  m_gammaCountDialog = 0;
}//void deleteGammaCountDialog()




#if( USE_SPECRUM_FILE_QUERY_WIDGET )
void InterSpec::showFileQueryDialog()
{
  if( m_specFileQueryDialog )
    return;
  
  m_specFileQueryDialog = new AuxWindow( "Spectrum File Query Tool", AuxWindowProperties::TabletNotFullScreen );
  //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
  //  size to be the min-size
  m_specFileQueryDialog->setMinimumSize( 640, 480 );
  m_specFileQueryDialog->setResizable( true );
  //m_specFileQueryDialog->disableCollapse();
  
  SpecFileQueryWidget *qw = new SpecFileQueryWidget( this );
  WGridLayout *stretcher = m_specFileQueryDialog->stretcher();
  stretcher->addWidget( qw, 0, 0 );
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  stretcher->setVerticalSpacing( 0 );
  stretcher->setHorizontalSpacing( 0 );
  
  WPushButton *closeButton = m_specFileQueryDialog->addCloseButtonToFooter( "Close", true );
  closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_specFileQueryDialog ) );
  
  m_specFileQueryDialog->finished().connect( this, &InterSpec::deleteFileQueryDialog );
  
  if( m_renderedWidth > 100 && m_renderedHeight > 100 && !isPhone() )
  {
    m_specFileQueryDialog->resize( 0.85*m_renderedWidth, 0.85*m_renderedHeight );
    m_specFileQueryDialog->centerWindow();
  }
  
  AuxWindow::addHelpInFooter( m_specFileQueryDialog->footer(), "spectrum-file-query" );
}//void showFileQueryDialog()


void InterSpec::deleteFileQueryDialog()
{
  if( !m_specFileQueryDialog )
    return;
  delete m_specFileQueryDialog;
  m_specFileQueryDialog = 0;
}//void deleteFileQueryDialog()
#endif




PopupDivMenu * InterSpec::displayOptionsPopupDiv()
{
  return m_displayOptionsPopupDiv;
}//WContainerWidget *displayOptionsPopupDiv()


void InterSpec::logMessage( const Wt::WString& message, const Wt::WString& source, int priority )
{
#if( PERFORM_DEVELOPER_CHECKS )
  static std::mutex s_message_mutex;
  
  {//begin codeblock to logg message
    std::lock_guard<std::mutex> file_gaurd( s_message_mutex );
    ofstream output( "interspec_messages_to_users.txt", ios::out | ios::app );
    const boost::posix_time::ptime now = WDateTime::currentDateTime().toPosixTime();
    output << "Message " << SpecUtils::to_iso_string( now ) << " ";
    if( !source.empty() )
      output << "[" << source.toUTF8() << ", " << priority << "]: ";
    else
      output << "[" << priority << "]: ";
    output << message.toUTF8() << endl << endl;
  }//end codeblock to logg message
#endif
  
  if( wApp )
  {
    m_messageLogged.emit( message, source, priority );
//    wApp->triggerUpdate();
  }else
  {
    const boost::posix_time::ptime now = WDateTime::currentDateTime().toPosixTime();
    cerr << "Message " << SpecUtils::to_iso_string( now ) << " ";
    if( !source.empty() )
      cerr << "[" << source.toUTF8() << ", " << priority << "]: ";
    else
      cerr << "[" << priority << "]: ";
    cerr << message.toUTF8() << endl << endl;
  }
}//void InterSpec::logMessage(...)



Wt::Signal< Wt::WString, Wt::WString, int > &InterSpec::messageLogged()
{
  return m_messageLogged;
}


Wt::WContainerWidget *InterSpec::menuDiv()
{
  return m_menuDiv;
}


void InterSpec::showWarningsWindow()
{
  m_warnings->createContent();
  
  if( !m_warningsWindow )
  {
    m_warningsWindow = new AuxWindow( "Notification/Logs",
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                   | AuxWindowProperties::IsModal
                   | AuxWindowProperties::EnableResize) );
    m_warningsWindow->contents()->setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
    m_warningsWindow->rejectWhenEscapePressed();
    m_warningsWindow->stretcher()->addWidget( m_warnings, 0, 0, 1, 1 );
    m_warningsWindow->stretcher()->setContentsMargins(0,0,0,0);
    //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
    //  size to be the min-size
    m_warningsWindow->setMinimumSize( 640, 480 );
    m_warningsWindow->resizeScaledWindow(0.75, 0.75);
    m_warningsWindow->centerWindow();
    m_warningsWindow->finished().connect( boost::bind( &InterSpec::handleWarningsWindowClose, this, false ) );
        
      
    WPushButton *clearButton = new WPushButton( "Delete Logs", m_warningsWindow->footer() );
    clearButton->clicked().connect( boost::bind( &WarningWidget::clearMessages, m_warnings ) );
    clearButton->addStyleClass( "BinIcon" );
      if (isMobile())
      {
          clearButton->setFloatSide( Right );
      }
      WPushButton *closeButton = m_warningsWindow->addCloseButtonToFooter();
      closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_warningsWindow ) );

  }//if( !m_warningsWindow )
  
  m_warningsWindow->show();
}//void showWarningsWindow()


void InterSpec::handleWarningsWindowClose( bool close )
{
    m_warningsWindow->stretcher()->removeWidget( m_warnings );
    AuxWindow::deleteAuxWindow( m_warningsWindow );
    m_warningsWindow = 0;
}//void handleWarningsWindowClose( bool )


void InterSpec::showPeakInfoWindow()
{
  if( m_toolsTabs )
  {
    cerr << "nterSpec::showPeakInfoWindow()\n\tTemporary hack - we dont want to show the"
         << " peak info window when tool tabs are showing since things get funky"
         << endl;
    m_toolsTabs->setCurrentWidget( m_peakInfoDisplay );
    m_currentToolsTab = m_toolsTabs->currentIndex();
    return;
  }//if( m_toolsTabs )
 
  if( m_toolsTabs && m_peakInfoDisplay )
  {
    m_toolsTabs->removeTab( m_peakInfoDisplay );
    m_toolsTabs->setCurrentIndex( 0 );
    m_currentToolsTab = 0;
  }//if( m_toolsTabs && m_peakInfoDisplay )
  
  if( !m_peakInfoWindow )
  {
    m_peakInfoWindow = new AuxWindow( "Peak Manager" );
    m_peakInfoWindow->rejectWhenEscapePressed();
    WBorderLayout *layout = new WBorderLayout();
    layout->setContentsMargins(0, 0, 15, 0);
    m_peakInfoWindow->contents()->setLayout( layout );
    //m_peakInfoWindow->contents()->setPadding(0);
    //m_peakInfoWindow->contents()->setMargin(0);
    
    layout->addWidget( m_peakInfoDisplay, Wt::WBorderLayout::Center );
    WContainerWidget *buttons = new WContainerWidget();
    layout->addWidget( buttons, Wt::WBorderLayout::South );

    WContainerWidget* footer = m_peakInfoWindow->footer();
    
    WPushButton* closeButton = m_peakInfoWindow->addCloseButtonToFooter("Close",true);
    closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_peakInfoWindow ) );
      
    AuxWindow::addHelpInFooter( m_peakInfoWindow->footer(), "peak-manager" );
    
    WPushButton *b = new WPushButton( CalibrationTabTitle, footer );
    b->clicked().connect( this, &InterSpec::showEnergyCalWindow );
    b->setFloatSide(Wt::Right);
      
    b = new WPushButton( GammaLinesTabTitle, footer );
    b->clicked().connect( this, &InterSpec::showGammaLinesWindow );
    b->setFloatSide(Wt::Right);
      
      
    m_peakInfoWindow->resizeScaledWindow( 0.75, 0.5 );
    m_peakInfoWindow->centerWindow();
    m_peakInfoWindow->setResizable( true );
    m_peakInfoWindow->finished().connect( this, &InterSpec::handlePeakInfoClose );
  }//if( !m_peakInfoWindow )

  m_peakInfoWindow->show();
}//void showPeakInfoWindow()


void InterSpec::handlePeakInfoClose()
{
  if( !m_peakInfoWindow )
    return;

  if( m_toolsTabs )
  {
    m_peakInfoWindow->contents()->removeWidget( m_peakInfoDisplay );
    if( m_toolsTabs->indexOf(m_peakInfoDisplay) < 0 )
    {
      m_toolsTabs->addTab( m_peakInfoDisplay, PeakInfoTabTitle, TabLoadPolicy );
      m_toolsTabs->setCurrentIndex( m_toolsTabs->indexOf(m_peakInfoDisplay) );
    }//if( m_toolsTabs->indexOf(m_peakInfoDisplay) < 0 )
    m_currentToolsTab = m_toolsTabs->currentIndex();
    
    delete m_peakInfoWindow;
    m_peakInfoWindow = NULL;
  }//if( m_toolsTabs )
}//void handlePeakInfoClose()


#if( USE_DB_TO_STORE_SPECTRA )
void InterSpec::finishStoreStateInDb( WLineEdit *nameedit,
                                           WTextArea *descriptionarea,
                                           AuxWindow *window,
                                           const bool forTesting,
                                           Dbo::ptr<UserState> parent )
{
  const WString &name = nameedit->text();
  if( name.empty() )
  {
    passMessage( "You must specify a name", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( name.empty() )
  
  const WString &desc = (!descriptionarea ? "" : descriptionarea->text());
  
  Dbo::ptr<UserState> state = serializeStateToDb( name, desc, forTesting, parent );
  
  if( parent )
    m_currentStateID = parent.id();
  else
    m_currentStateID = state.id();
    
  updateSaveWorkspaceMenu();
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  if( forTesting )
  {
    string filepath = SpecUtils::append_path(TEST_SUITE_BASE_DIR, "analysis_tests");
    
    if( !SpecUtils::is_directory( filepath ) )
      throw runtime_error( "CWD didnt contain a 'analysis_tests' folder as expected" );
    
    const boost::posix_time::ptime localtime = boost::posix_time::second_clock::local_time();
    
    string timestr = SpecUtils::to_iso_string( localtime );
    string::size_type period_pos = timestr.find_last_of( '.' );
    timestr.substr( 0, period_pos );
    
    filepath = SpecUtils::append_path( filepath, (name.toUTF8() + "_" + timestr + ".n42") );
    
#ifdef _WIN32
    const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16(filepath);
    ofstream output( wfilepath.c_str(), ios::binary|ios::out );
#else
    ofstream output( filepath.c_str(), ios::binary|ios::out );
#endif
    if( !output.is_open() )
      throw runtime_error( "Couldnt open test file ouput '"
                           + filepath + "'" );
    
    storeTestStateToN42( output, name.toUTF8(), desc.toUTF8() );
  }//if( forTesting )
#endif
  if (window)
      delete window;
}//void finishStoreStateInDb()


#if( INCLUDE_ANALYSIS_TEST_SUITE )
void InterSpec::storeTestStateToN42( std::ostream &output,
                                          const std::string &name,
                                          const std::string &descr )
{
  try
  {
    if( !m_dataMeasurement )
      throw runtime_error( "No valid foreground file" );
    
    SpecMeas meas;
    meas.uniqueCopyContents( *m_dataMeasurement );
    
    
#if( PERFORM_DEVELOPER_CHECKS )
    try
    {
      SpecMeas::equalEnough( meas, *m_dataMeasurement );
    }catch( std::exception &e )
    {
      log_developer_error( __func__, e.what() );
    }
#endif
    
    const set<int> foregroundsamples = displayedSamples( SpecUtils::SpectrumType::Foreground );
    const set<int> backgroundsamples = displayedSamples( SpecUtils::SpectrumType::Background );
    set<int> newbacksamples;
    
    for( const std::shared_ptr<const SpecUtils::Measurement> &p : meas.measurements() )
    {
      const int samplenum = p->sample_number();
      if( foregroundsamples.count(samplenum) )
        meas.set_source_type( SpecUtils::SourceType::Foreground, p );
      else
        meas.remove_measurement( p, false );
    }//for( const std::shared_ptr<const SpecUtils::Measurement> &p : meas.measurements() )
    
    
    std::shared_ptr< std::deque< std::shared_ptr<const PeakDef> > >
    foregroundpeaks, backgroundpeaks;
    foregroundpeaks = meas.peaks( foregroundsamples );
    
    if( !!m_backgroundMeasurement )
    {
      backgroundpeaks = m_backgroundMeasurement->peaks( backgroundsamples );
      
      std::vector< std::shared_ptr<const SpecUtils::Measurement> > backgrounds = m_backgroundMeasurement->measurements();
      for( const std::shared_ptr<const SpecUtils::Measurement> &p : backgrounds )
      {
        if( !backgroundsamples.count( p->sample_number()) )
          continue;
        
        std::shared_ptr<SpecUtils::Measurement> backmeas = std::make_shared<SpecUtils::Measurement>( *p );
        meas.add_measurement( backmeas, false );
        meas.set_source_type( SpecUtils::SourceType::Background, backmeas );
        newbacksamples.insert( backmeas->sample_number() );
      }//for( const std::shared_ptr<const SpecUtils::Measurement> &p : m_backgroundMeasurement->measurements() )
    }//if( !!m_backgroundMeasurement )
 
    //Now remove all peaks not for the currently displayed samples.
    meas.removeAllPeaks();
    
    //But add back in peaks for just the displayed foreground and background
    if( !!foregroundpeaks && foregroundpeaks->size() )
      meas.setPeaks( *foregroundpeaks, foregroundsamples );
    if( !!backgroundpeaks && backgroundpeaks->size() )
      meas.setPeaks( *backgroundpeaks, newbacksamples );
    
    using namespace ::rapidxml;
    std::shared_ptr< xml_document<char> > n42doc;
    n42doc = meas.create_2012_N42_xml();
    
    if( !n42doc )
      throw runtime_error( "Failed to create 2012 N42 XML document" );
    
    xml_node<char> *RadInstrumentData = n42doc->first_node( "RadInstrumentData", 17 );
    
    if( !RadInstrumentData )
      throw runtime_error( "Logic error trying to retrieve RadInstrumentData node" );
    
    xml_node<char> *InterSpecNode = RadInstrumentData->first_node( "DHS:InterSpec", 13 );
    if( !InterSpecNode )
      throw runtime_error( "Logic error trying to retrieve DHS:InterSpec node" );
    
    //Make a note that this contains information for testing as well
    xml_attribute<char> *attr = n42doc->allocate_attribute( "ForTesting", "true" );
    InterSpecNode->append_attribute( attr );
    
    //Add user options into the XML
    m_user->userOptionsToXml( InterSpecNode, this );
    
    if( newbacksamples.size() )
    {
      vector<int> backsamples( newbacksamples.begin(), newbacksamples.end() );
      stringstream backsamplesstrm;
      for( size_t i = 0; i < backsamples.size(); ++i )
        backsamplesstrm << (i?" ":"") << backsamples[i];
      const char *val = n42doc->allocate_string( backsamplesstrm.str().c_str() );
      xml_node<char> *backgroundsamples_node
                    = n42doc->allocate_node( node_element,
                                      "DisplayedBackgroundSampleNumber", val );
      InterSpecNode->append_node( backgroundsamples_node );
    }//if( newbacksamples.size() )
    
    if( m_shieldingSourceFit
        && m_shieldingSourceFit->userChangedDuringCurrentForeground() )
    {
      m_shieldingSourceFit->serialize( InterSpecNode );
      cerr << "\n\nThe shielding/source fit model WAS changed for current foreground" << endl;
    }else
    {
      cerr << "\n\nThe shielding/source fit model was NOT changed for current foreground" << endl;
    }
    
    const boost::posix_time::ptime localtime = boost::posix_time::second_clock::local_time();
    const string timestr = SpecUtils::to_iso_string( localtime );
    
    const char *val = n42doc->allocate_string( timestr.c_str(), timestr.size()+1 );
    xml_node<char> *node = n42doc->allocate_node( node_element, "TestSaveDateTime", val );
    InterSpecNode->append_node( node );
    
    if( name.size() )
    {
      val = n42doc->allocate_string( name.c_str(), name.size()+1 );
      node = n42doc->allocate_node( node_element, "TestName", val );
      InterSpecNode->append_node( node );
    }//if( name.size() )
    
    if( descr.size() )
    {
      val = n42doc->allocate_string( descr.c_str(), descr.size()+1 );
      node = n42doc->allocate_node( node_element, "TestDescription", val );
      InterSpecNode->append_node( node );
    }//if( name.size() )
    
    string xml_data;
    rapidxml::print( std::back_inserter(xml_data), *n42doc, 0 );
    
    if( !output.write( xml_data.c_str(), xml_data.size() ) )
      throw runtime_error( "" );
  }catch( std::exception &e )
  {
    string msg = "Failed to save test state to N42 file: ";
    msg += e.what();
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void storeTestStateToN42( const std::string &name )



void InterSpec::loadTestStateFromN42( const std::string filename )
{
  using namespace rapidxml;
  try
  {
    std::shared_ptr<SpecMeas> meas = std::make_shared<SpecMeas>();
    
    // TODO: If we call load_from_N42_document(...), the sample numbers get all messed up and we
    //       cant match peaks and displayed samples up right - I'm not sure why this is at the
    //       moment, probably some virtual call issue or something, but for the moment will just
    //       patch through it and re-parse the XML a second time.
    //const bool loaded = meas->load_from_N42_document( doc_node );
    const bool loaded = meas->load_file( filename, SpecUtils::ParserType::N42_2012, ".n42" );
    
    if( !loaded )
      throw runtime_error( "Failed to load SpecMeas from XML document" );
    
    
    std::vector<char> data;
    SpecUtils::load_file_data( filename.c_str(), data );
    
    xml_document<char> doc;
    char * const data_start = &data.front();
    char * const data_end = data_start + data.size();
    doc.parse<rapidxml::parse_trim_whitespace|rapidxml::allow_sloppy_parse>( data_start, data_end );
    xml_node<char> *doc_node = doc.first_node();
    
    const xml_node<char> *RadInstrumentData = doc.first_node( "RadInstrumentData", 17 );
    
    if( !RadInstrumentData )
      throw runtime_error( "Error trying to retrieve RadInstrumentData node" );
    
    const xml_node<char> *InterSpecNode = RadInstrumentData->first_node( "DHS:InterSpec", 13 );
    if( !InterSpecNode )
      throw runtime_error( "Error trying to retrieve DHS:InterSpec node" );

    //meas->decodeSpecMeasStuffFromXml( InterSpecNode );
    
    const xml_attribute<char> *attr = InterSpecNode->first_attribute( "ForTesting" );
    if( !attr || !attr->value()
        || !rapidxml::internal::compare(attr->value(), attr->value_size(), "true", 4, true) )
      throw runtime_error( "Input N42 file doesnt look to be a test state" );
    
    set<int> backgroundsamplenums;
    const xml_node<char> *backsample_node
               = InterSpecNode->first_node( "DisplayedBackgroundSampleNumber" );
    
    if( backsample_node && backsample_node->value_size() )
    {
      std::vector<int> results;
      const bool success = SpecUtils::split_to_ints(
                                      backsample_node->value(),
                                      backsample_node->value_size(), results );
      if( !success )
        throw runtime_error( "Error parsing background sample numbers" );
      backgroundsamplenums.insert( results.begin(), results.end() );
    }//if( backsample_node && backsample_node->value_size() )
    
    
    const xml_node<char> *preferences = InterSpecNode->first_node( "preferences" );
    if( preferences )
      InterSpecUser::restoreUserPrefsFromXml( m_user, preferences, this );
    else
      cerr << "Warning, couldnt find preferences in the XML"
              " when restoring state from XML" << endl;
    
    
//    {
//      std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks;
//      peaks = meas->peaks( meas->displayedSampleNumbers() );
//      if( !!peaks )
//        cerr << "There were NO peaks for displayed sample numbers" << endl;
//      else
//        cerr << "There were " << peaks->size() << " peaks for displayed sample numbers" << endl;
//    
//      cerr << "Displaying sample numbers: ";
//      for( int i : meas->displayedSampleNumbers() )
//        cerr << i << ", ";
//      cerr << endl;
      
//      std::set<std::set<int> > samps = meas->sampleNumsWithPeaks();
//      cerr << "There are samps.size()=" << samps.size() << endl;
//      for( auto a : samps )
//      {
//        cerr << "Set: ";
//        for( auto b : a )
//          cerr << b << " ";
//        cerr << endl;
//      }
//    }
    
    
    
    const xml_node<char> *name = InterSpecNode->first_node( "TestName" );
    const xml_node<char> *descrip = InterSpecNode->first_node( "TestDescription" );
    const xml_node<char> *savedate = InterSpecNode->first_node( "TestSaveDateTime" );
    
    
    std::shared_ptr<SpecMeas> dummy;
    setSpectrum( dummy, {}, SpecUtils::SpectrumType::Background, 0 );
    setSpectrum( dummy, {}, SpecUtils::SpectrumType::SecondForeground, 0 );
    
    string filename = meas->filename();
    if( name && name->value_size() )
      filename = name->value();
    
    std::shared_ptr<SpectraFileHeader> header;
    header = m_fileManager->addFile( filename, meas );
    
    const std::set<int> &dispsamples = meas->displayedSampleNumbers();
    
    setSpectrum( meas, dispsamples, SpecUtils::SpectrumType::Foreground, 0 );
    
    if( backgroundsamplenums.size() )
      setSpectrum( meas, backgroundsamplenums, SpecUtils::SpectrumType::Background, 0 );
    
    const xml_node<char> *sourcefit = InterSpecNode->first_node( "ShieldingSourceFit" );
    if( sourcefit )
    {
      showShieldingSourceFitWindow();
      m_shieldingSourceFit->deSerialize( sourcefit );
      closeShieldingSourceFitWindow();
    }//if( sourcefit )
    
    stringstream msg;
    msg << "Loaded test state";
    if( name && name->value_size() )
      msg << " named '" << name->value() << "'";
    if( savedate && savedate->value_size() )
      msg << " saved " << savedate->value();
    if( descrip && descrip->value_size() )
      msg << " with description " << descrip->value();
    
    passMessage( msg.str(), "", WarningWidget::WarningMsgInfo );
  }catch( std::exception &e )
  {
    string msg = "InterSpec::loadTestStateFromN42: ";
    msg += e.what();
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void InterSpec::loadTestStateFromN42()

namespace
{
  void doTestStateLoad( WSelectionBox *filesbox, AuxWindow *window, InterSpec *viewer )
  {
    const string path_to_tests = SpecUtils::append_path( TEST_SUITE_BASE_DIR, "analysis_tests" );
    const string filename = SpecUtils::append_path( path_to_tests, filesbox->currentText().toUTF8() );
    viewer->loadTestStateFromN42( filename );
    delete window;
  }
}

void InterSpec::startN42TestStates()
{
  const string path_to_tests = SpecUtils::append_path( TEST_SUITE_BASE_DIR, "analysis_tests" );
  const vector<string> files = SpecUtils::recursive_ls( path_to_tests, ".n42" );
  
  if( files.empty() )
  {
    passMessage( "There are no test files to load in the 'analysis_tests' "
                 "directory", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( files.empty() )
  
  AuxWindow *window = new AuxWindow( "Test State N42 Files" );
  window->resizeWindow( 450, 400 );
  
  WGridLayout *layout = window->stretcher();
  WSelectionBox *filesbox = new WSelectionBox();
  layout->addWidget( filesbox, 0, 0, 1, 2 );
  
  //Lets display files by alphabetical name
  vector<string> dispfiles;
  for( const string &p : files )
    dispfiles.push_back( SpecUtils::fs_relative( path_to_tests, p ) );
    
  std::sort( begin(dispfiles), end(dispfiles) );
  
  for( const string &name : dispfiles )
    filesbox->addItem( name );
  
  WPushButton *button = new WPushButton( "Cancel" );
  layout->addWidget( button, 1, 0, AlignCenter );
  button->clicked().connect( boost::bind(&AuxWindow::deleteAuxWindow, window) );
  
  button = new WPushButton( "Load" );
  button->disable();
  button->clicked().connect( boost::bind( &doTestStateLoad, filesbox, window, this ) );
  
  layout->addWidget( button, 1, 1, AlignCenter );
  filesbox->activated().connect( button, &WPushButton::enable );
  
  layout->setRowStretch( 0, 1 );
  
  window->show();
  window->centerWindow();
}//void startN42TestStates()



void InterSpec::startStoreTestStateInDb()
{
  if( m_currentStateID >= 0 )
  {
    AuxWindow *window = new AuxWindow( "Warning",
                                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
                                       | AuxWindowProperties::DisableCollapse) );
    WText *t = new WText( "Overwrite current test state?", window->contents() );
    t->setInline( false );
    
    window->setClosable( false );

    WPushButton *button = new WPushButton( "Overwrite", window->footer() );
    button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    button->clicked().connect( boost::bind( &InterSpec::startStoreStateInDb, this, true, false, true, false ) );
    
    button = new WPushButton( "Create New", window->footer() );
    button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    button->clicked().connect( boost::bind( &InterSpec::startStoreStateInDb, this, true, true, false, false ) );
    
    button = new WPushButton( "Cancel", window->footer() );
    button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    window->centerWindow();
    window->show();
    
    window->rejectWhenEscapePressed();
    window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  }else
  {
    startStoreStateInDb( true, true, true, false );
  }//if( m_currentStateID >= 0 ) / else
}//void startStoreTestStateInDb()
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE )


//Save the snapshot and ALSO the spectra.
void InterSpec::stateSave()
{
  if( m_currentStateID > 0 )
  {
    //TODO: Should check if can save?
    saveShieldingSourceModelToForegroundSpecMeas();
    startStoreStateInDb( false, false, false, false ); //save snapshot
  }else
  {
    //No currentStateID, so just save as new snapshot
    stateSaveAs();
  }//if( m_currentStateID > 0 ) / else
} //void stateSave()


//Copied from startStoreStateInDb()
//This is the new combined Save As method (snapshot/spectra)
void InterSpec::stateSaveAs()
{
  const bool forTesting = false;
    
  AuxWindow *window = new AuxWindow( "Store State As", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::TabletNotFullScreen) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  WGridLayout *layout = window->stretcher();
    
  WLineEdit *edit = new WLineEdit();
  edit->setEmptyText( "(Name to store under)" );
  WText *label = new WText( "Name" );
  layout->addWidget( label, 2, 0 );
  layout->addWidget( edit,  2, 1 );
    
  if( !!m_dataMeasurement && m_dataMeasurement->filename().size() )
      edit->setText( m_dataMeasurement->filename() );
    
  WTextArea *summary = new WTextArea();
  label = new WText( "Desc." );
  summary->setEmptyText( "(Optional description)" );
  layout->addWidget( label, 3, 0 );
  layout->addWidget( summary, 3, 1 );
  
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 3, 1 );
  
  label = new WText( "This will save the current app state to the database - to export<br/> as a separate file, use the &quot;Export File&quot; menu option" );
  label->setInline( false );
  label->setTextAlignment( Wt::AlignCenter );
  layout->addWidget( label, 4, 1 );
  
  WPushButton *closeButton = window->addCloseButtonToFooter("Cancel", false);
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  WPushButton *save = new WPushButton( "Save" , window->footer() );
  save->setIcon( "InterSpec_resources/images/disk2.png" );
  
  save->clicked().connect( boost::bind( &InterSpec::stateSaveAsAction,
                                        this, edit, summary, window, forTesting ) );
    
  window->setMinimumSize(WLength(450), WLength(250));
    
  window->centerWindow();
  window->show();
}//void stateSaveAs()


//New method to save tag for snapshot
void InterSpec::stateSaveTag()
{
  AuxWindow *window = new AuxWindow( "Tag current state",
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal) | AuxWindowProperties::TabletNotFullScreen) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  window->disableCollapse();
  WGridLayout *layout = window->stretcher();
    
  WLineEdit *edit = new WLineEdit();
  if( !isMobile() )
    edit->setFocus(true);
  WText *label = new WText( "Name" );
  layout->addWidget( label, 2, 0 );
  layout->addWidget( edit,  2, 1 );
    
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 2, 1 );
  
  WPushButton *closeButton = window->addCloseButtonToFooter("Cancel", false);
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  WPushButton *save = new WPushButton( "Tag" , window->footer() );
  //save->setIcon( "InterSpec_resources/images/disk2.png" );
    
  save->clicked().connect( boost::bind( &InterSpec::stateSaveTagAction,
                                         this, edit, window ) );
    
  window->setMinimumSize(WLength(450), WLength::Auto);
    
  window->centerWindow();
  window->show();
}//void stateSaveTag()


void InterSpec::stateSaveTagAction( WLineEdit *nameedit, AuxWindow *window )
{
  Dbo::ptr<UserState> state;
    
  if( m_currentStateID >= 0 )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      state = m_sql->session()->find<UserState>( "where id = ?" )
                               .bind( m_currentStateID );
      transaction.commit();
    }catch( std::exception &e )
    {
        cerr << "\nstartStoreStateInDb() caught: " << e.what() << endl;
    }//try / catch
  }//if( m_currentStateID >= 0 )
    
  //Save Snapshot/enviornment
  if( state )
    finishStoreStateInDb( nameedit, 0, 0, false, state );
  else
    passMessage( "You must first save a state before creating a snapshot.",
                 "", WarningWidget::WarningMsgHigh );
      
  //close window
  AuxWindow::deleteAuxWindow( window );
}//stateSaveTagAction()


//Action to Save As Snapshot
void InterSpec::stateSaveAsAction( WLineEdit *nameedit,
                                      WTextArea *descriptionarea,
                                      AuxWindow *window,
                                      const bool forTesting )
{
  //Need to remove the current foregorund/background/second from the file
  //  manager, and then re-add back in the SpecMeas so they will be saved to
  //  new database entries.  Otherwise saving the current state as, will create
  //  a new UserState that points to the current database entries for the
  //  SpecMeas
  for( int i = 0; i < 3; i++ )
  {
    std::shared_ptr<SpecMeas> m = measurment( SpecUtils::SpectrumType(i) );
    if( m )
    {
      m_fileManager->removeSpecMeas( m, false );
      m_fileManager->addFile( m->filename(), m );
      m->setModified();
    }
  }//for( int i = 0; i < 3; i++ )
  
  //Save Snapshot/enviornment
  finishStoreStateInDb( nameedit, descriptionarea, NULL, forTesting, Dbo::ptr<UserState>() );
    
  //close window
  AuxWindow::deleteAuxWindow(window);
}//stateSaveAsAction()



//Note: Used indirectly by new combined snapshot method stateSave()
void InterSpec::startStoreStateInDb( const bool forTesting,
                                          const bool asNewState,
                                          const bool allowOverWrite,
                                          const bool endOfSessionNoDelete )
{
  Dbo::ptr<UserState> state;
  
  if( m_currentStateID >= 0 && !asNewState )
  {
      //If saving a tag, make sure to save to the parent.  We never overwrite TAGS
      try
      {
          DataBaseUtils::DbTransaction transaction( *m_sql );
          Dbo::ptr<UserState> childState  = m_sql->session()->find<UserState>( "where id = ? AND SnapshotTagParent_id >= 0" )
          .bind( m_currentStateID );
          
          if (childState)
          {
              state = m_sql->session()->find<UserState>( "where id = ?" )
              .bind( childState.get()->snapshotTagParent.id() );
          }
          
          transaction.commit();
      }catch( std::exception &e )
      {
          cerr << "\nstartStoreStateInDb() caught: " << e.what() << endl;
      }//try / catch

    if (!state)
    {
        try
        {
            DataBaseUtils::DbTransaction transaction( *m_sql );
            state = m_sql->session()->find<UserState>( "where id = ?" )
                           .bind( m_currentStateID );
            transaction.commit();
        }catch( std::exception &e )
        {
            cerr << "\nstartStoreStateInDb() caught: " << e.what() << endl;
        }//try / catch
    } //!state
  }//if( m_currentStateID >= 0 )
  
  if( state )
  { //Save without dialog
      
    const bool writeProtected = state->isWriteProtected();
    
    if( allowOverWrite && writeProtected )
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      UserState::removeWriteProtection( state, m_sql->session() );
      transaction.commit();
    }//if( allowOverWrite && writeProtected )
    
      if (endOfSessionNoDelete)
      {
          DataBaseUtils::DbTransaction transaction( *m_sql );
          state.modify()->stateType = UserState::kEndOfSessionHEAD;
          transaction.commit();
      }
      
    try
    {
      saveStateToDb( state );
    }catch( std::exception &e )
    {
      passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    }//try / catch
    
    if( writeProtected )
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      UserState::makeWriteProtected( state, m_sql->session() );
      transaction.commit();
    }//if( writeProtected )
    
    
    passMessage( "Saved state '" + state->name + "'", "",
                 WarningWidget::WarningMsgSave );
    return;
  }//if( state )
  
  AuxWindow *window = new AuxWindow( "Create Snapshot",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                      | AuxWindowProperties::TabletNotFullScreen
                                      | AuxWindowProperties::DisableCollapse) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  WGridLayout *layout = window->stretcher();
  WLineEdit *edit = new WLineEdit();
  edit->setEmptyText( "(Name to store under)" );
  WText *label = new WText( "Name" );
  layout->addWidget( label, 0, 0 );
  layout->addWidget( edit,  0, 1 );
  
  if( !!m_dataMeasurement && m_dataMeasurement->filename().size() )
    edit->setText( m_dataMeasurement->filename() );
  
  WTextArea *summary = new WTextArea();
  label = new WText( "Desc." );
  summary->setEmptyText( "(Optional description)" );
  layout->addWidget( label, 1, 0 );
  layout->addWidget( summary, 1, 1 );
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 1, 1 );
  
  
  WPushButton *closeButton = window->addCloseButtonToFooter("Cancel", false);
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );

  WPushButton *save = new WPushButton( "Create" , window->footer() );
  save->setIcon( "InterSpec_resources/images/disk2.png" );
  
  save->clicked().connect( boost::bind( &InterSpec::finishStoreStateInDb,
                                       this, edit, summary, window, forTesting, Wt::Dbo::ptr<UserState>() ) );

  window->setMinimumSize(WLength(450), WLength(250));
  
  window->centerWindow();
  window->show();
}//void InterSpec::startStoreStateInDb()


//Deprecated - no longer used
void InterSpec::browseDatabaseStates( bool allowTests )
{
  new DbStateBrowser( this, allowTests );
}//void InterSpec::browseDatabaseStates()
#endif //#if( USE_DB_TO_STORE_SPECTRA )


#if( USE_DB_TO_STORE_SPECTRA && INCLUDE_ANALYSIS_TEST_SUITE )
void InterSpec::startStateTester()
{
  new SpectrumViewerTesterWindow( this );
}//void startStateTester()
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )




void InterSpec::addFileMenu( WWidget *parent, const bool isAppTitlebar )
{
  if( m_fileMenuPopup )
    return;

  const bool mobile = isMobile();
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addFileMenu(): parent passed in must"
                         " be a PopupDivMenu  or WContainerWidget" );

  string menuname = isAppTitlebar ? "File" : "InterSpec";
#if( !defined(__APPLE__) && USING_ELECTRON_NATIVE_MENU )
  menuname = "File";
#else
  if( mobile )
    menuname = "Spectra";
#endif
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( menuname, menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_fileMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
    
//    button = new WPushButton( "Spectrum", menuDiv );
//    button->addStyleClass( "MenuLabel" );
//    m_fileMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_fileMenuPopup = parentMenu->addPopupMenuItem( menuname );
//    m_fileMenuPopup = parentMenu->addPopupMenuItem( "Spectrum" );
  }//if( menuDiv ) / else

  PopupDivMenuItem *item = (PopupDivMenuItem *)0;
  
#if( defined(__APPLE__) && USING_ELECTRON_NATIVE_MENU )
  PopupDivMenuItem *aboutitem = m_fileMenuPopup->createAboutThisAppItem();
  
  if( aboutitem )
    aboutitem->triggered().connect( boost::bind( &InterSpec::showLicenseAndDisclaimersWindow, this, false, std::function<void()>{} ) );
  
  m_fileMenuPopup->addSeparator();
  m_fileMenuPopup->addRoleMenuItem( PopupDivMenu::MenuRole::Hide );
  m_fileMenuPopup->addRoleMenuItem( PopupDivMenu::MenuRole::HideOthers );
  m_fileMenuPopup->addRoleMenuItem( PopupDivMenu::MenuRole::UnHide );
  m_fileMenuPopup->addRoleMenuItem( PopupDivMenu::MenuRole::Front );
  m_fileMenuPopup->addSeparator();
#elif( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    item = m_fileMenuPopup->addMenuItem( "About InterSpec" );
    item->triggered().connect( boost::bind( &InterSpec::showLicenseAndDisclaimersWindow, this, false, std::function<void()>{} ) );
    m_fileMenuPopup->addSeparator();  //doesnt seem to be showing up for some reason... owe well.
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif
  
  
#if( USE_DB_TO_STORE_SPECTRA )
  //  item = m_fileMenuPopup->addMenuItem( "Test Serialization" );
  //  item->triggered().connect( this, &InterSpec::testLoadSaveState );
  
  //----temporary
  //    PopupDivMenu *backupmenu = m_fileMenuPopup->addPopupMenuItem( "Backup", "InterSpec_resources/images/database_go.png" );
  
  //    m_fileMenuPopup->addSeparator();
  
  if( m_fileManager )
  {
    if( mobile )
    {
      item = m_fileMenuPopup->addMenuItem( "Open File..." );
      item->triggered().connect( boost::bind ( &SpecMeasManager::startQuickUpload, m_fileManager ) );
      
      item = m_fileMenuPopup->addMenuItem( "Loaded Spectra..." );
      item->triggered().connect( this, &InterSpec::showCompactFileManagerWindow );
    }//if( mobile )
    
    
    item = m_fileMenuPopup->addMenuItem( "Manager...", "InterSpec_resources/images/file_manager_small.png" );
    HelpSystem::attachToolTipOn(item, "Manage loaded spectra", showToolTips );
    
    item->triggered().connect( m_fileManager, &SpecMeasManager::startSpectrumManager );
    
    m_fileMenuPopup->addSeparator();
    
    const char *save_txt = "Store";           //"Save"
    const char *save_as_txt = "Store As...";  //"Save As..."
    const char *tag_txt = "Tag...";
    
    // --- new save menu ---
    m_saveState = m_fileMenuPopup->addMenuItem( save_txt );
    m_saveState->triggered().connect( boost::bind( &InterSpec::stateSave, this ) );
    HelpSystem::attachToolTipOn(m_saveState, "Saves the current app state to InterSpecs database", showToolTips );
    
    
    // --- new save as menu ---
    m_saveStateAs = m_fileMenuPopup->addMenuItem( save_as_txt );
    m_saveStateAs->triggered().connect( boost::bind( &InterSpec::stateSaveAs, this ) );
    HelpSystem::attachToolTipOn(m_saveStateAs, "Saves the current app state to a new listing in InterSpecs database", showToolTips );
    m_saveStateAs->setDisabled( true );
    
    // --- new save tag menu ---
    m_createTag = m_fileMenuPopup->addMenuItem( tag_txt , "InterSpec_resources/images/stop_time_small.png");
    m_createTag->triggered().connect( boost::bind(&InterSpec::stateSaveTag, this ));
    
    HelpSystem::attachToolTipOn(m_createTag, "Tags the current Interspec state so you "
                                "can revert to at a later time.  When loading an app state, "
                                "you can choose either the most recent save, or select past tagged versions.", showToolTips );
    
    m_fileMenuPopup->addSeparator();
    
    item = m_fileMenuPopup->addMenuItem( "Previous..." , "InterSpec_resources/images/db_small.png");
    item->triggered().connect( boost::bind( &SpecMeasManager::browseDatabaseSpectrumFiles,
                                            m_fileManager,
                                            SpecUtils::SpectrumType::Foreground ) );
    HelpSystem::attachToolTipOn(item, "Opens previously saved states", showToolTips );
    
    m_fileMenuPopup->addSeparator();
  }//if( m_fileManager )
#endif  //USE_DB_TO_STORE_SPECTRA


  item = m_fileMenuPopup->addMenuItem( "Clear Session..." );
  item->triggered().connect( this, &InterSpec::startClearSession );
  HelpSystem::attachToolTipOn(item, "Starts a clean application state with no spectra loaded", showToolTips );
  m_fileMenuPopup->addSeparator();
  
  
  PopupDivMenu *subPopup = 0;

  subPopup = m_fileMenuPopup->addPopupMenuItem( "Samples" );
  
  const string docroot = wApp->docRoot();
  
  item = subPopup->addMenuItem( "Ba-133 (16k channel)" );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/ba133_source_640s_20100317.n42"),
                                         SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::N42_2006 ) );
  if( mobile )
    item = subPopup->addMenuItem( "Passthrough (16k channel)" );
  else
    item = subPopup->addMenuItem( "Passthrough (16k bin ICD1, 8 det., 133 samples)" );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/passthrough.n42"),
                                         SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::N42_2006 ) );
  
  item = subPopup->addMenuItem( "Background (16k bin N42)" );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/background_20100317.n42"),
                                         SpecUtils::SpectrumType::Background, SpecUtils::ParserType::N42_2006 ) );
  //If its a mobile device, we'll give a few more spectra to play with
  if( mobile )
  {
    item = subPopup->addMenuItem( "Ba-133 (Low Res, No Calib)" );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Ba133LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    
    item = subPopup->addMenuItem( "Co-60 (Low Res, No Calib)" );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Co60LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    
    item = subPopup->addMenuItem( "Cs-137 (Low Res, No Calib)" );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Cs137LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    item = subPopup->addMenuItem( "Th-232 (Low Res, No Calib)" );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Th232LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
  }//if( mobile )
  
  
  if( !mobile )
  {
    string tip = "Creates a dialog to browse for a spectrum file on your file system.";
    if( !mobile )
      tip += " Drag and drop the file directly into the app window as a quicker alternative.";
    
    item = m_fileMenuPopup->addMenuItem( "Open File..." );
    HelpSystem::attachToolTipOn( item, tip, showToolTips );
    item->triggered().connect( boost::bind ( &SpecMeasManager::startQuickUpload, m_fileManager ) );
  }//if( !mobile )
  

#if( USE_SAVEAS_FROM_MENU )
  //only display when not on desktop.
  m_downloadMenu = m_fileMenuPopup->addPopupMenuItem( "Export File" );
  
  for( SpecUtils::SpectrumType i = SpecUtils::SpectrumType(0);
      i <= SpecUtils::SpectrumType::Background;
      i = SpecUtils::SpectrumType(static_cast<int>(i)+1) )
  {
    m_downloadMenus[static_cast<int>(i)] = m_downloadMenu->addPopupMenuItem( descriptionText( i ) );
    
    for( SpecUtils::SaveSpectrumAsType j = SpecUtils::SaveSpectrumAsType(0);
        j < SpecUtils::SaveSpectrumAsType::NumTypes;
        j = SpecUtils::SaveSpectrumAsType(static_cast<int>(j)+1) )
    {
      const string desc = descriptionText( j );
      item = m_downloadMenus[static_cast<int>(i)]->addMenuItem( desc + " File" );
      DownloadCurrentSpectrumResource *resource
      = new DownloadCurrentSpectrumResource( i, j, this, item );
      
#if( USE_OSX_NATIVE_MENU || USING_ELECTRON_NATIVE_MENU )
      //If were using OS X native menus, we obviously cant relly on the
      //  browser responding to a click on an anchor; we also cant click the
      //  link of the PopupDivMenuItem itself in javascript, or we get a
      //  cycle of the server telling the browser to click the link again,
      //  each time it tells it to do it in javascript.  So we will introduce
      //  a false menu item that will actually have the URL of the download
      //  associated with its anchor, and not connect any server side
      //  triggered() events to it, thus breaking the infiinite cycle - blarg.
      PopupDivMenuItem *fakeitem = new PopupDivMenuItem( "", "" );
      m_downloadMenus[static_cast<int>(i)]->WPopupMenu::addItem( fakeitem );
      fakeitem->setLink( WLink( resource ) );
      fakeitem->setLinkTarget( TargetNewWindow );
      
      if( fakeitem->anchor() )
      {
        const string jsclick = "try{document.getElementById('"
        + fakeitem->anchor()->id()
        + "').click();}catch(e){}";
        item->triggered().connect( boost::bind( &WApplication::doJavaScript, wApp, jsclick, true ) );
      }else
      {
        cerr << "Unexpected error accessing the anchor for file downloading" << endl;
      }
#else
      item->setLink( WLink( resource ) );
      item->setLinkTarget( TargetNewWindow );
#endif
      
      const char *tooltip = nullptr;
      switch( j )
      {
        case SpecUtils::SaveSpectrumAsType::Txt:
          tooltip = "A space delimited file will be created with a header of"
          " things like live time, followed by three columns"
          " (channel number, lower channel energy, counts). Each"
          " record in the file will be sequentially recorded with"
          " three line breaks between records.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::Csv:
          tooltip = "A comma delimietd file will be created with a channel"
          " lower energy and a counts column.  All records in the"
          " current spectrum file will be written sequentially,"
          " with an extra line break and &quot;<b>Energy, Data</b>"
          "&quot; header between subsequent records.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::Pcf:
          tooltip = "A binary PCF file will be created which contains all"
          " records of the current file.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::N42_2006:
          tooltip = "A simple spectromiter style 2006 N42 XML file will be"
          " produced which contains all records in the current file.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::N42_2012:
          tooltip = "A 2012 N42 XML document will be produced which contains"
          " all samples of the current spectrum file, as well as"
          " some additional <code>InterSpec</code> information"
          " such as the identified peaks and detector response.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::Chn:
          tooltip = "A binary integer CHN file will be produced containing a"
          " single spectrum matching what is currently shown.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
          tooltip = "A binary floating point SPC file will be produced containing a"
          " single spectrum matching what is currently shown.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
          tooltip = "A binary integer SPC file will be produced containing a"
          " single spectrum matching what is currently shown.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::SpcAscii:
          tooltip = "A ascii based SPC file will be produced containing a"
          " single spectrum matching what is currently shown.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0:
          tooltip = "A binary Exploranium GR130 file will be produced with"
          " each record (spectrum) containing 256 channels.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2:
          tooltip = "A binary Exploranium GR135v2 file (includes neutron info)"
          " will be produced with each record (spectrum) containing"
          " 1024 channels.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::SpeIaea:
          tooltip = "A ASCII based standard file format that will contain a"
          " single spectrum only.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::Cnf:
          tooltip = "A binary Canberra file format that contains only a single spectrum.";
          break;
          
        case SpecUtils::SaveSpectrumAsType::Tka:
          tooltip = "A text based format that provides real time and channel counts only.";
          break;
          
#if( SpecUtils_ENABLE_D3_CHART )
        case SpecUtils::SaveSpectrumAsType::HtmlD3:
          tooltip = "An HTML file using D3.js to generate a spectrum chart"
          " that you can optionally interact with and view offline.";
          break;
#endif //#if( SpecUtils_ENABLE_D3_CHART )
        case SpecUtils::SaveSpectrumAsType::NumTypes:
          break;
      }//switch( j )
      
      assert( tooltip );
      
      const bool showInstantly = true;
      if( tooltip )
        HelpSystem::attachToolTipOn( item, tooltip, showInstantly,
                                    HelpSystem::ToolTipPosition::Right );
    }//for( loop over file types )
    
    m_downloadMenus[static_cast<int>(i)]->disable();
    if( m_downloadMenus[static_cast<int>(i)]->parentItem() )
      m_downloadMenu->setItemHidden( m_downloadMenus[static_cast<int>(i)]->parentItem(), true );
  }//for( loop over spectrums )
  
  m_fileMenuPopup->setItemHidden(m_downloadMenu->parentItem(),true); //set as hidden first, will be visible when spectrum is added
#elif( IOS )
  PopupDivMenuItem *exportFile = m_fileMenuPopup->addMenuItem( "Export Spectrum..." );
  exportFile->triggered().connect( this, &InterSpec::exportSpecFile );
#else
  item = m_fileMenuPopup->addMenuItem( "Export" );
  item->triggered().connect( m_fileManager, &SpecMeasManager::displayQuickSaveAsDialog );
#endif //#if( USE_SAVEAS_FROM_MENU )
  

  
#if( USE_DB_TO_STORE_SPECTRA )
  assert( !m_fileManager || m_saveStateAs );
  assert( !m_fileManager || m_saveState );
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  m_fileMenuPopup->addSeparator();
  PopupDivMenu* testing = m_fileMenuPopup->addPopupMenuItem( "Testing" , "InterSpec_resources/images/testing.png");
  item = testing->addMenuItem( "Store App Test State..." );
  item->triggered().connect( boost::bind( &InterSpec::startStoreTestStateInDb, this ) );
  HelpSystem::attachToolTipOn(item,"Stores app state as part of the automated test suite", showToolTips );
  
  item = testing->addMenuItem( "Restore App Test State..." );
  item->triggered().connect( boost::bind(&InterSpec::browseDatabaseStates, this, true) );
  HelpSystem::attachToolTipOn(item, "Restores InterSpec to a previously stored state.", showToolTips );
  
  item = testing->addMenuItem( "Load N42 Test State..." );
  item->triggered().connect( boost::bind(&InterSpec::startN42TestStates, this) );
  HelpSystem::attachToolTipOn(item, "Loads a InterSpec specific variant of a "
                                    "2012 N42 file that contains Foreground, "
                                    "Background, User Settings, and Shielding/"
                                    "Source model.", showToolTips );
  
  item = testing->addMenuItem( "Show Testing Widget..." );
  item->triggered().connect( boost::bind(&InterSpec::startStateTester, this ) );
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE )


#if(USE_OSX_NATIVE_MENU )
    //Add a separtor before Quit InterSpec
    m_fileMenuPopup->addSeparator();
#endif

#if( BUILD_AS_ELECTRON_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
#if( USING_ELECTRON_NATIVE_MENU )
  m_fileMenuPopup->addSeparator();
  m_fileMenuPopup->addRoleMenuItem( PopupDivMenu::MenuRole::Quit );
#elif( BUILD_AS_ELECTRON_APP )
  m_fileMenuPopup->addSeparator();
  PopupDivMenuItem *exitItem = m_fileMenuPopup->addMenuItem( "Exit" );
  exitItem->triggered().connect( std::bind( []{
    ElectronUtils::send_nodejs_message( "CloseWindow", "" );
  }) );
#endif
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif // BUILD_AS_ELECTRON_APP
} // void InterSpec::addFileMenu( WContainerWidget *menuDiv, bool show_examples )


bool InterSpec::toolTabsVisible() const
{
  return m_toolsTabs;
}

void InterSpec::addToolsTabToMenuItems()
{
  if( !m_toolsMenuPopup )
    return;
  
  if( m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::RefPhotopeaks)] )
    return;
  
  bool showToolTips = false;
  
  try
  {
    showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  }catch(...)
  {
  }
  
  const int index_offest = isMobile() ? 1 : 0;
  
  Wt::WMenuItem *item = nullptr;
  const char *tooltip = nullptr, *icon = nullptr;
  
  icon = "InterSpec_resources/images/reflines.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 0, GammaLinesTabTitle, icon , true );
  item->triggered().connect( boost::bind( &InterSpec::showGammaLinesWindow, this ) );
  tooltip = "Allows user to display x-rays and/or gammas from elements, isotopes, or nuclear reactions."
            " Also provides user with a shortcut to change detector and account for shielding.";
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::RefPhotopeaks)] = item;
  
  icon = "InterSpec_resources/images/peakmanager.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 1, PeakInfoTabTitle, icon, true );
  item->triggered().connect( this, &InterSpec::showPeakInfoWindow );
  tooltip = "Provides shortcuts to search for and identify peaks. Displays parameters of all fit peaks in a sortable table.";
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::PeakManager)] = item;
  
  icon = "InterSpec_resources/images/calibrate.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 2, CalibrationTabTitle, icon, true );
  item->triggered().connect( this, &InterSpec::showEnergyCalWindow );
  tooltip = "Allows user to modify or fit for offset, linear, or quadratic energy calibration terms,"
            " as well as edit non-linear deviation pairs.<br />"
            "Energy calibration can also be done graphically by right-click dragging the spectrum from original to modified energy"
            " (e.g., drag the data peak to the appropriate reference line).";
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::EnergyCal)] = item;
  
  icon = "InterSpec_resources/images/magnifier_black.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 3, NuclideSearchTabTitle, icon, true );
  item->triggered().connect( this, &InterSpec::showNuclideSearchWindow);
  tooltip = "Search for nuclides with constraints on energy, branching ratio, and half life.";
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::NuclideSearch)] = item;

  
  icon = "InterSpec_resources/images/auto_peak_search.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 4, "Auto Peak Search", icon, true );
  item->triggered().connect( boost::bind( &PeakSearchGuiUtils::automated_search_for_peaks, this, true ) );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::AutoPeakSearch)] = item;
  
  
  item = m_toolsMenuPopup->addSeparatorAt( index_offest + 5 );
  
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::Seperator)] = item;
}//void addToolsTabToMenuItems()


void InterSpec::removeToolsTabToMenuItems()
{
  if( !m_toolsMenuPopup )
    return;
  
  for( int i = 0; i < static_cast<int>(ToolTabMenuItems::NumItems); ++i )
  {
    if( !m_tabToolsMenuItems[i] )
      continue;
 
    if( m_tabToolsMenuItems[i]->isSeparator() )
      m_toolsMenuPopup->removeSeperator( m_tabToolsMenuItems[i] );
    
    delete m_tabToolsMenuItems[i];
    m_tabToolsMenuItems[i] = nullptr;
  }
  
}//void removeToolsTabToMenuItems()


void InterSpec::setToolTabsVisible( bool showToolTabs )
{
#if( USE_CSS_FLEX_LAYOUT )
  if( m_toolsTabs->isVisible() == showToolTabs )
    return;
#else
  if( static_cast<bool>(m_toolsTabs) == showToolTabs )
    return;
#endif
  
  m_toolTabsVisibleItems[0]->setHidden( showToolTabs );
  m_toolTabsVisibleItems[1]->setHidden( !showToolTabs );
  
  if( showToolTabs )
    removeToolsTabToMenuItems();
  else
    addToolsTabToMenuItems();

#if( USE_CSS_FLEX_LAYOUT )
  m_toolsTabs->setHidden( !showToolTabs );
  m_toolsResizer->setHidden( !showToolTabs );
#else
  
  const int layoutVertSpacing = isMobile() ? 10 : 5;
  
  if( showToolTabs )
  { //We are showing the tool tabs
    if( m_menuDiv )
      m_layout->removeWidget( m_menuDiv );
  
    m_layout->removeWidget( m_charts );
    m_layout->clear();
    
    //We have to completely replace m_layout or else the vertical spacing
    //  changes wont quite trigger.  Prob a Wt bug, so should investigate again
    //  if we ever upgrade Wt.
    delete m_layout;
    m_layout = new WGridLayout();
    m_layout->setContentsMargins( 0, 0, 0, 0 );
    m_layout->setHorizontalSpacing( 0 );
    setLayout( m_layout );
    
    if( m_peakInfoWindow )
    {
      m_peakInfoWindow->contents()->removeWidget( m_peakInfoDisplay );
      delete m_peakInfoWindow;
      m_peakInfoWindow = NULL;
    }//if( m_peakInfoWindow )
    
    string refNucXmlState;
    if( m_referencePhotopeakLines
        && (!m_referencePhotopeakLines->currentlyShowingNuclide().empty()
             || m_referencePhotopeakLines->persistedNuclides().size()) )
    {
      m_referencePhotopeakLines->serialize( refNucXmlState );
    }
    
    closeGammaLinesWindow();
    
    if( m_energyCalWindow )
    {
      m_energyCalWindow->stretcher()->removeWidget( m_energyCalTool );
      delete m_energyCalWindow;
      m_energyCalWindow = nullptr;
    }//if( m_energyCalWindow )
    
    m_toolsTabs = new WTabWidget();
    //m_toolsTabs->addStyleClass( "ToolsTabs" );
    
    CompactFileManager *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::LeftToRight );
    m_toolsTabs->addTab( compact, FileTabTitle, TabLoadPolicy );
    
    m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                   boost::placeholders::_1, boost::placeholders::_2 ) );
    
    //WMenuItem * peakManTab =
    m_toolsTabs->addTab( m_peakInfoDisplay, PeakInfoTabTitle, TabLoadPolicy );
//    const char *tooltip = "Displays parameters of all identified peaks in a sortable table.";
//    HelpSystem::attachToolTipOn( peakManTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
    if( m_referencePhotopeakLines )
    {
      m_referencePhotopeakLines->clearAllLines();
      delete m_referencePhotopeakLines;
    }
      
    if( m_referencePhotopeakLinesWindow )
      delete m_referencePhotopeakLinesWindow;
    m_referencePhotopeakLinesWindow = NULL;
      
    m_referencePhotopeakLines = new ReferencePhotopeakDisplay( m_spectrum,
                                              m_materialDB.get(),
                                              m_shieldingSuggestion,
                                              this );
    setReferenceLineColors( nullptr );
    
    //PreLoading is necessary on the m_referencePhotopeakLines widget, so that the
    //  "Isotope Search" widget will work properly when a nuclide is clicked
    //  on to display its photpeaks
    //XXX In Wt 3.3.4 at least, the contents of m_referencePhotopeakLines
    //  are not actually loaded to the client until the tab is clicked, and I
    //  cant seem to get this to actually happen.
    //WMenuItem *refPhotoTab =
    m_toolsTabs->addTab( m_referencePhotopeakLines, GammaLinesTabTitle, TabLoadPolicy );
      
//      const char *tooltip = "Allows user to display x-rays and/or gammas from "
//                            "elements, isotopes, or nuclear reactions. Also "
//                            "provides user with a shortcut to change detector "
//                            "and account for shielding.";
//      HelpSystem::attachToolTipOn( refPhotoTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
      
    m_toolsTabs->currentChanged().connect( this, &InterSpec::handleToolTabChanged );
    
    m_energyCalTool->setWideLayout();
    m_toolsTabs->addTab( m_energyCalTool, CalibrationTabTitle, TabLoadPolicy );
    
    m_toolsLayout = new WGridLayout();
    m_toolsLayout->setContentsMargins( 0, 0, 0, 0 );
    m_toolsLayout->setVerticalSpacing( 0 );
    m_toolsLayout->setHorizontalSpacing( 0 );
    m_toolsLayout->addWidget( m_toolsTabs, 0, 0 );
    
    int row = 0;
    if( m_menuDiv )
      m_layout->addWidget( m_menuDiv,  row++, 0 );
    
    m_layout->addWidget( m_charts, row, 0 );
    m_layout->setRowResizable( row, true );
    m_layout->setRowStretch( row, 1 );
    
    m_layout->setVerticalSpacing( layoutVertSpacing );
    if( m_menuDiv && !m_menuDiv->isHidden() )  //get rid of a small amount of space between the menu bar and the chart
      m_charts->setMargin( -layoutVertSpacing, Wt::Top );
    
    //Without using the wrapper below, the tabs widget will change height, even
    //  if it is explicitly set, when different tabs are clicked (unless a
    //  manual resize is performed by the user first)
    //m_layout->addLayout( toolsLayout,   ++row, 0 );
    WContainerWidget *toolsWrapper = new WContainerWidget();
    toolsWrapper->setLayout( m_toolsLayout );
    toolsWrapper->setHeight( 245 );
    m_layout->addWidget( toolsWrapper, ++row, 0 );
    m_layout->setRowStretch( row, 0 );
    
    if( m_nuclideSearchWindow )
    {
      m_nuclideSearch->clearSearchEnergiesOnClient();
      m_nuclideSearchWindow->stretcher()->removeWidget( m_nuclideSearch );
      delete m_nuclideSearchWindow;
      m_nuclideSearchWindow = 0;
    }//if( m_nuclideSearchWindow )
    
    assert( !m_nuclideSearchContainer );
    
    m_nuclideSearchContainer = new WContainerWidget();
    WGridLayout *isoSearchLayout = new WGridLayout();
    m_nuclideSearchContainer->setLayout( isoSearchLayout );
    isoSearchLayout->setContentsMargins( 0, 0, 0, 0 );
    isoSearchLayout->addWidget( m_nuclideSearch, 0, 0 );
    m_nuclideSearchContainer->setMargin( 0 );
    m_nuclideSearchContainer->setPadding( 0 );
    isoSearchLayout->setRowStretch( 0, 1 );
    isoSearchLayout->setColumnStretch( 0, 1 );
    
    //WMenuItem *nuclideTab =
    m_toolsTabs->addTab( m_nuclideSearchContainer, NuclideSearchTabTitle, TabLoadPolicy );
//    const char *tooltip = "Search for nuclides with constraints on energy, "
//                          "branching ratio, and half life.";
//    HelpSystem::attachToolTipOn( nuclideTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
    //nuclideTab->setIcon("InterSpec_resources/images/magnifier.png");
//    if( m_nuclideSearch )
//      m_nuclideSearch->loadSearchEnergiesToClient();
    
    //Make sure the current tab is the peak info display
    m_toolsTabs->setCurrentWidget( m_peakInfoDisplay );
    
    if( m_referencePhotopeakLines && refNucXmlState.size() )
      m_referencePhotopeakLines->deSerialize( refNucXmlState );
  }else
  {
    //We are hiding the tool tabs
    m_layout->removeWidget( m_charts );
    
    if( m_menuDiv )
      m_layout->removeWidget( m_menuDiv );
    m_toolsTabs->removeTab( m_peakInfoDisplay );
    m_toolsTabs->removeTab( m_energyCalTool );
    
    m_nuclideSearch->clearSearchEnergiesOnClient();
    m_nuclideSearchContainer->layout()->removeWidget( m_nuclideSearch );
    m_toolsTabs->removeTab( m_nuclideSearchContainer );
    delete m_nuclideSearchContainer;
    m_nuclideSearchContainer = nullptr;
    
    if( m_referencePhotopeakLines )
      m_referencePhotopeakLines->clearAllLines();
    
    m_toolsTabs = nullptr;
    
    if( !m_referencePhotopeakLinesWindow )
      m_referencePhotopeakLines = nullptr;
    m_toolsLayout = nullptr;
    
    m_layout->clear();
    
    //We have to completely replace m_layout or else the vertical spacing
    //  changes wont quite trigger.  Prob a Wt bug, so should investigate again
    //  if we ever upgrade Wt.
    //delete m_layout;
    m_layout = new WGridLayout();
    m_layout->setContentsMargins( 0, 0, 0, 0 );
    m_layout->setHorizontalSpacing( 0 );
    setLayout( m_layout );
    
    int row = -1;
    if( m_menuDiv )
      m_layout->addWidget( m_menuDiv, ++row, 0 );
    m_layout->addWidget( m_charts, ++row, 0 );
    m_layout->setRowStretch( row, 1 );
  }//if( showToolTabs ) / else
  
  
  // I'm guessing when the charts were temporarily removed from the DOM (or changed or whatever),
  //  the bindings to watch for mousedown and touchstart were removed, so lets re-instate them.
#if( USE_CSS_FLEX_LAYOUT )
#else
  m_charts->doJavaScript( "Wt.WT.InitFlexResizer('" + m_chartResizer->id() + "','" + m_timeSeries->id() + "');" );
#endif
  
  if( m_toolsTabs )
    m_currentToolsTab = m_toolsTabs->currentIndex();
  else
    m_currentToolsTab = -1;
    
    
    if (m_mobileBackButton && m_mobileForwardButton)
    {
        if (!toolTabsVisible() && m_dataMeasurement && !m_dataMeasurement->passthrough()
             && m_dataMeasurement->sample_numbers().size()> 1)
        {
            m_mobileBackButton->setHidden(false);
            m_mobileForwardButton->setHidden(false);
        }
        else
        {
            m_mobileBackButton->setHidden(true);
            m_mobileForwardButton->setHidden(true);
            
        }
    }
  
  //Not sure _why_ this next statement is needed, but it is, or else the
  //  spectrum chart shows up with no data
  m_spectrum->scheduleUpdateForeground();
  m_spectrum->scheduleUpdateBackground();
  m_spectrum->scheduleUpdateSecondData();
  
  m_timeSeries->scheduleRenderAll();
#endif // USE_CSS_FLEX_LAYOUT / else
}//void setToolTabsVisible( bool showToolTabs )


void InterSpec::addDisplayMenu( WWidget *parent )
{
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addDisplayMenu(): parent passed in"
                        " must be a PopupDivMenu  or WContainerWidget" );
 
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( "View", menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_displayOptionsPopupDiv = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_displayOptionsPopupDiv = parentMenu->addPopupMenuItem( "View" );
  }//if( menuDiv ) / else
  
  m_toolTabsVisibleItems[0] = m_displayOptionsPopupDiv->addMenuItem( "Show Tool Tabs" , "InterSpec_resources/images/dock_small.png" );
  m_toolTabsVisibleItems[1] = m_displayOptionsPopupDiv->addMenuItem( "Hide Tool Tabs" , "InterSpec_resources/images/undock_small.png" );
  m_toolTabsVisibleItems[0]->triggered().connect( boost::bind( &InterSpec::setToolTabsVisible, this, true ) );
  m_toolTabsVisibleItems[1]->triggered().connect( boost::bind( &InterSpec::setToolTabsVisible, this, false ) );

  
  if( isPhone() )
  {
    //hide Dock mode on small phone screens
    m_toolTabsVisibleItems[0]->setHidden(true);
  } //isPhone
    
  m_toolTabsVisibleItems[1]->setHidden(true);
  
  m_displayOptionsPopupDiv->addSeparator();

  addDetectorMenu( m_displayOptionsPopupDiv );
  
  m_displayOptionsPopupDiv->addSeparator();
  
  PopupDivMenu *chartmenu = m_displayOptionsPopupDiv->addPopupMenuItem( "Chart Options" , "InterSpec_resources/images/spec_settings_small.png");
  
  bool logypref = true;
  try{ logypref = m_user->preferenceValue<bool>( "LogY" ); }catch(...){}
  
  m_logYItems[0] = chartmenu->addMenuItem( "Log Y Scale" );
  m_logYItems[1] = chartmenu->addMenuItem( "Linear Y Scale" );
  m_logYItems[0]->setHidden( logypref );
  m_logYItems[1]->setHidden( !logypref );
  m_logYItems[0]->triggered().connect( boost::bind( &InterSpec::setLogY, this, true  ) );
  m_logYItems[1]->triggered().connect( boost::bind( &InterSpec::setLogY, this, false ) );
  m_spectrum->setYAxisLog( logypref );
  std::function<void (boost::any)> logy_fcn = [=](boost::any value){
    this->setLogY( boost::any_cast<bool>(value) );
  };
  InterSpecUser::associateFunction( m_user, "LogY", logy_fcn, this );

  
  const bool verticleLines = InterSpecUser::preferenceValue<bool>( "ShowVerticalGridlines", this );
  m_verticalLinesItems[0] = chartmenu->addMenuItem( "Show Vertical Lines" , "InterSpec_resources/images/sc_togglegridvertical.png");
  m_verticalLinesItems[1] = chartmenu->addMenuItem( "Hide Vertical Lines" , "InterSpec_resources/images/sc_togglegridvertical.png");
  m_verticalLinesItems[0]->triggered().connect( boost::bind( &InterSpec::setVerticalLines, this, true ) );
  m_verticalLinesItems[1]->triggered().connect( boost::bind( &InterSpec::setVerticalLines, this, false ) );
  m_verticalLinesItems[0]->setHidden( verticleLines );
  m_verticalLinesItems[1]->setHidden( !verticleLines );
  m_spectrum->showVerticalLines( verticleLines );
  m_timeSeries->showVerticalLines( verticleLines );
  std::function<void (boost::any)> vl_fcn = [=](boost::any value){
    this->setVerticalLines( boost::any_cast<bool>(value) );
  };
  InterSpecUser::associateFunction( m_user, "ShowVerticalGridlines", vl_fcn, this );
  
  const bool horizontalLines = InterSpecUser::preferenceValue<bool>( "ShowHorizontalGridlines", this );
  m_horizantalLinesItems[0] = chartmenu->addMenuItem( "Show Horizontal Lines" , "InterSpec_resources/images/sc_togglegridhorizontal.png");
  m_horizantalLinesItems[1] = chartmenu->addMenuItem( "Hide Horizontal Lines" , "InterSpec_resources/images/sc_togglegridhorizontal.png");
  m_horizantalLinesItems[0]->triggered().connect( boost::bind( &InterSpec::setHorizantalLines, this, true ) );
  m_horizantalLinesItems[1]->triggered().connect( boost::bind( &InterSpec::setHorizantalLines, this, false ) );
  m_horizantalLinesItems[0]->setHidden( horizontalLines );
  m_horizantalLinesItems[1]->setHidden( !horizontalLines );
  m_spectrum->showHorizontalLines( horizontalLines );
  m_timeSeries->showHorizontalLines( horizontalLines );
  std::function<void (boost::any)> hl_fcn = [=](boost::any value){
    this->setHorizantalLines( boost::any_cast<bool>(value) );
  };
  InterSpecUser::associateFunction( m_user, "ShowHorizontalGridlines", hl_fcn, this );
  
  
  if( isPhone() )
  {
    m_compactXAxisItems[0] = m_compactXAxisItems[1] = nullptr;
  }else
  {
    const bool makeCompact = InterSpecUser::preferenceValue<bool>( "CompactXAxis", this );
    m_spectrum->setCompactAxis( makeCompact );
    m_compactXAxisItems[0] = chartmenu->addMenuItem( "Compact X-Axis" , "");
    m_compactXAxisItems[1] = chartmenu->addMenuItem( "Normal X-Axis" , "");
    m_compactXAxisItems[0]->triggered().connect( boost::bind( &InterSpec::setXAxisCompact, this, true ) );
    m_compactXAxisItems[1]->triggered().connect( boost::bind( &InterSpec::setXAxisCompact, this, false ) );
    m_compactXAxisItems[0]->setHidden( makeCompact );
    m_compactXAxisItems[1]->setHidden( !makeCompact );
    std::function<void (boost::any)> cxa_fcn = [=](boost::any value){
      this->setXAxisCompact( boost::any_cast<bool>(value) );
    };
    InterSpecUser::associateFunction( m_user, "CompactXAxis", cxa_fcn, this );
  }
  
  //What we should do here is have a dialog that pops up that lets users  select
  //  colors for foreground, background, secondary, as well as the first number
  //  of reference lines.
  //  Probably also chart background, chart text, chart axis.  Should allow
  //  users to save these as "themes", as well as include a few by default...
  //Probably should also include an option to say whether fitted phtopeaks
  //  should adopt the color of the current refernce photopeak lines
  
  //With the below, it seems the color of the widget never actually gets
  //  changed... also, it wouldnt be supported on native menus.
  //chartmenu->addSeparator();
  //ColorSelect *color = new ColorSelect(ColorSelect::PrefferNative);
  //chartmenu->addWidget( color );
  //color->cssColorChanged().connect(<#const std::string &function#>)
  
  chartmenu->addSeparator();

  PopupDivMenuItem *item = chartmenu->addMenuItem( "Show Spectrum Legend" );
  m_spectrum->legendDisabled().connect( item, &PopupDivMenuItem::show );
  m_spectrum->legendEnabled().connect( item,  &PopupDivMenuItem::hide );
  
  item->triggered().connect( boost::bind( &D3SpectrumDisplayDiv::enableLegend, m_spectrum ) );
  
  item->hide(); //we are already showing the legend
  
  addPeakLabelSubMenu( m_displayOptionsPopupDiv ); //add Peak menu
  
  const bool showSlider = InterSpecUser::preferenceValue<bool>( "ShowXAxisSlider", this );
  m_spectrum->showXAxisSliderChart( showSlider );
  m_showXAxisSliderItems[0] = m_displayOptionsPopupDiv->addMenuItem( "Show Energy Slider" , "");
  m_showXAxisSliderItems[1] = m_displayOptionsPopupDiv->addMenuItem( "Hide Energy Slider" , "");
  m_showXAxisSliderItems[0]->triggered().connect( boost::bind( &InterSpec::setXAxisSlider, this, true ) );
  m_showXAxisSliderItems[1]->triggered().connect( boost::bind( &InterSpec::setXAxisSlider, this, false ) );
  m_showXAxisSliderItems[0]->setHidden( showSlider );
  m_showXAxisSliderItems[1]->setHidden( !showSlider );
  std::function<void (boost::any)> xas_fcn = [=](boost::any value){
    this->setXAxisSlider( boost::any_cast<bool>(value) );
  };
  InterSpecUser::associateFunction( m_user, "ShowXAxisSlider", xas_fcn, this );
  
  
  const bool showScalers = InterSpecUser::preferenceValue<bool>( "ShowYAxisScalers", this );
  m_spectrum->showYAxisScalers( showScalers );
  
  m_showYAxisScalerItems[0] = m_displayOptionsPopupDiv->addMenuItem( "Show Y-Axis Scalers" , "");
  m_showYAxisScalerItems[1] = m_displayOptionsPopupDiv->addMenuItem( "Hide Y-Axis Scalers" , "");
  m_showYAxisScalerItems[0]->triggered().connect( boost::bind( &InterSpec::setShowYAxisScalers, this, true ) );
  m_showYAxisScalerItems[1]->triggered().connect( boost::bind( &InterSpec::setShowYAxisScalers, this, false ) );
  m_showYAxisScalerItems[0]->setHidden( showScalers );
  m_showYAxisScalerItems[1]->setHidden( !showScalers );
  m_showYAxisScalerItems[(showScalers ? 1 : 0)]->disable();
  
  std::function<void (boost::any)> fcnt = [=](boost::any value){
    this->setShowYAxisScalers( boost::any_cast<bool>(value) );
  };
  
  InterSpecUser::associateFunction( m_user, "ShowYAxisScalers", fcnt, this );
  
  m_displayOptionsPopupDiv->addSeparator();
  
  m_backgroundSubItems[0] = m_displayOptionsPopupDiv->addMenuItem( "Background Subtract" );
  m_backgroundSubItems[1] = m_displayOptionsPopupDiv->addMenuItem( "Un-Background Subtract" );
  m_backgroundSubItems[0]->triggered().connect( boost::bind( &InterSpec::setBackgroundSub, this, true ) );
  m_backgroundSubItems[1]->triggered().connect( boost::bind( &InterSpec::setBackgroundSub, this, false ) );
  m_backgroundSubItems[0]->disable();
  m_backgroundSubItems[1]->hide();
  
  
  m_hardBackgroundSub = m_displayOptionsPopupDiv->addMenuItem( "Hard Background Sub..." );
  m_hardBackgroundSub->triggered().connect( this, &InterSpec::startHardBackgroundSub );
  m_hardBackgroundSub->disable();
  
#if( USE_GOOGLE_MAP )
  m_mapMenuItem = m_displayOptionsPopupDiv->addMenuItem( "Map","InterSpec_resources/images/map_small.png" );
  m_mapMenuItem->triggered().connect( boost::bind( &InterSpec::createMapWindow, this, SpecUtils::SpectrumType::Foreground ) );
  m_mapMenuItem->disable();
  HelpSystem::attachToolTipOn( m_mapMenuItem,
                    "Show measurement(s) location on a Google Map. Only enabled"
                    " when GPS coordinates are available.", showToolTips );
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  m_searchMode3DChart = m_displayOptionsPopupDiv->addMenuItem( "Show 3D View","" );
  m_searchMode3DChart->triggered().connect( boost::bind( &InterSpec::create3DSearchModeChart, this ) );
  m_searchMode3DChart->disable();
  HelpSystem::attachToolTipOn( m_searchMode3DChart,
                        "Shows Time vs. Energy vs. Counts view for search mode or RPM data.", showToolTips );
#endif
  
  m_showRiidResults = m_displayOptionsPopupDiv->addMenuItem( "Show RIID Results","" );
  m_showRiidResults->triggered().connect( boost::bind( &InterSpec::showRiidResults, this, SpecUtils::SpectrumType::Foreground ) );
  HelpSystem::attachToolTipOn( m_showRiidResults,
                              "Shows the detectors RIID analysis results included in the spectrum file.", showToolTips );
  m_showRiidResults->disable();
  
  
  {
    m_featureMarkerMenuItem = m_displayOptionsPopupDiv->addMenuItem( "Feature Markers...", "", true );
    HelpSystem::attachToolTipOn( m_featureMarkerMenuItem,
                    "Tool to show single/double escape peaks, Compton peak, Compton Edge, and sum peaks",
                                showToolTips );
    m_featureMarkerMenuItem->triggered().connect( this, &InterSpec::toggleFeatureMarkerWindow );
    
    //If didnt want to use JSlot, could do...
    //    js = can + ".data('compangle',null);";
    //    checkbox->checked().connect( boost::bind( &WApplication::doJavaScript, wApp, js, true ) );
    
    m_displayOptionsPopupDiv->addSeparator();
    auto saveitem = m_displayOptionsPopupDiv->addMenuItem( "Save Spectrum as PNG" );
    saveitem->triggered().connect( boost::bind(&InterSpec::saveChartToImg, this, true, true) );
    saveitem = m_displayOptionsPopupDiv->addMenuItem( "Save Spectrum as SVG" );
    saveitem->triggered().connect( boost::bind(&InterSpec::saveChartToImg, this, true, false) );
  }//if( overlay )
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    m_displayOptionsPopupDiv->addSeparator();
    auto browserItem = m_displayOptionsPopupDiv->addMenuItem( "Use in external browser" );
#if( BUILD_AS_ELECTRON_APP )
    browserItem->triggered().connect( std::bind( [](){
      ElectronUtils::send_nodejs_message("OpenInExternalBrowser", "");
    } ) );
#endif
    
#if( BUILD_AS_OSX_APP )
    // A brief attempt at using javascript to open a browser window failed (probably because I wasnt
    //  doing it right or something), so I just implemented calling back to obj-c; see the
    //  didReceiveScriptMessage method in AppDelegate.mm.
    //  Alternatively we could have added something into macOsUtils.h and then called into there
    //  where we would have the obj-c open a browser window that way, but I wanted to try this
    //  method of communicating between JS and native code (but I shuld check if it introduces any
    //  notable overhead...).
    // A note for the future: should probably have tried to init the javascript:
    //    "try{document.getElementById('" + browserItem->anchor()->id() + "').click();}catch(e){}";
    //    after calling the following c++
    //      browserItem->setLink( WLink("http://localhost:port?restore=no&primary=no") );
    //      browserItem->setLinkTarget( AnchorTarget::TargetNewWindow );
    browserItem->triggered().connect( std::bind([=](){
      doJavaScript( "console.log('Will try to call back to obj-c');"
                    "try{"
                      "window.webkit.messageHandlers.interOp.postMessage({\"action\": \"ExternalInstance\"});"
                     "}catch(error){"
                       "console.warn('Failed to callback to the obj-c: ' + error );"
                     "}" );
    }) );
#endif //BUILD_AS_OSX_APP
  }//if( useNativeMenu )
#endif //BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP
  
  
#if( USING_ELECTRON_NATIVE_MENU )
  m_displayOptionsPopupDiv->addSeparator();
  m_displayOptionsPopupDiv->addRoleMenuItem( PopupDivMenu::MenuRole::ToggleFullscreen );
  m_displayOptionsPopupDiv->addSeparator();
  m_displayOptionsPopupDiv->addRoleMenuItem( PopupDivMenu::MenuRole::ResetZoom );
  m_displayOptionsPopupDiv->addRoleMenuItem( PopupDivMenu::MenuRole::ZoomIn );
  m_displayOptionsPopupDiv->addRoleMenuItem( PopupDivMenu::MenuRole::ZoomOut );
#if( PERFORM_DEVELOPER_CHECKS )
  m_displayOptionsPopupDiv->addSeparator();
  m_displayOptionsPopupDiv->addRoleMenuItem( PopupDivMenu::MenuRole::ToggleDevTools );
#endif
#elif( BUILD_AS_ELECTRON_APP )

  m_displayOptionsPopupDiv->addSeparator();
  PopupDivMenuItem *fullScreenItem = m_displayOptionsPopupDiv->addMenuItem( "Toggle Full Screen" );  //F11
  m_displayOptionsPopupDiv->addSeparator();
  PopupDivMenuItem *resetZoomItem = m_displayOptionsPopupDiv->addMenuItem( "Actual Size" ); //Ctrl+0
  PopupDivMenuItem *zoomInItem = m_displayOptionsPopupDiv->addMenuItem( "Zoom In" ); //Ctrl+Shift+=
  PopupDivMenuItem *zoomOutItem = m_displayOptionsPopupDiv->addMenuItem( "Zoom Out" ); //Ctrl+-

  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsResetPageZoom);
  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsIncreasePageZoom);
  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsDecreasePageZoom);
  
  // Note: the triggered() signal a Wt::Signal, which is C++ only, so we cant just hook it up to
  //       javascript for it to run - we have to make the round-trip JS -> C++ -> JS
  fullScreenItem->triggered().connect( std::bind( []{
    ElectronUtils::send_nodejs_message( "ToggleMaximizeWindow", "" );
  }) );
  
  resetZoomItem->triggered().connect( std::bind( []{
    wApp->doJavaScript( "Wt.WT.ResetPageZoom();" );
  }) );
  
  zoomInItem->triggered().connect( std::bind( []{
    wApp->doJavaScript( "Wt.WT.IncreasePageZoom();" );
  }) );
  
  zoomOutItem->triggered().connect( std::bind( []{
    wApp->doJavaScript( "Wt.WT.DecreasePageZoom();" );
  }) );

#if( BUILD_AS_ELECTRON_APP )
#if( PERFORM_DEVELOPER_CHECKS )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    m_displayOptionsPopupDiv->addSeparator();
    PopupDivMenuItem *devToolItem = m_displayOptionsPopupDiv->addMenuItem( "Toggle Dev Tools" );
    devToolItem->triggered().connect( std::bind( []{
      ElectronUtils::send_nodejs_message( "ToggleDevTools", "" );
    }) );
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif
#endif
  
#endif
}//void addDisplayMenu( menuParentDiv )


void InterSpec::addDetectorMenu( WWidget *menuWidget )
{
  if( m_detectorToShowMenu )
    return;
  
  PopupDivMenu *parentPopup = dynamic_cast<PopupDivMenu *>( menuWidget );

  //TODO: Change "Detectors" to "Detectors Displayed"
  if( parentPopup )
  {
    m_detectorToShowMenu = parentPopup->addPopupMenuItem( "Detectors" /*, "InterSpec_resources/images/detector_small_white.png"*/ );
  }else
  {
    WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>(menuWidget);
    if( !menuDiv )
      throw runtime_error( "addDetectorMenu: serious error" );
    
    WPushButton *button = new WPushButton( "Detectors", menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_detectorToShowMenu = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }//if( parentPopup )
  
  if( m_detectorToShowMenu->parentItem() )
    m_detectorToShowMenu->parentItem()->disable();
//  PopupDivMenuItem *item = m_detectorToShowMenu->addMenuItem( "Energy Calibration" );
//  item->triggered().connect( boost::bind( &WDialog::setHidden, m_energyCalWindow, false, WAnimation() ) );
//  item->triggered().connect( boost::bind( &InterSpec::showEnergyCalWindow, this ) );
}//void addDetectorMenu( WContainerWidget *menuDiv )


void InterSpec::handEnergyCalWindowClose()
{
  if( !m_energyCalWindow || !m_energyCalTool )
    return;
  
  WGridLayout *layout = m_energyCalWindow->stretcher();
  layout->removeWidget( m_energyCalTool );
  
  AuxWindow::deleteAuxWindow( m_energyCalWindow );
  m_energyCalWindow = nullptr;
  
  if( m_toolsTabs )
  {
    m_energyCalTool->setWideLayout();
    if( m_toolsTabs->indexOf(m_energyCalTool) < 0 )
      m_toolsTabs->addTab( m_energyCalTool, CalibrationTabTitle, TabLoadPolicy );
    
    m_currentToolsTab = m_toolsTabs->currentIndex();
  }//if( m_toolsTabs )
}//void handEnergyCalWindowClose()


EnergyCalTool *InterSpec::energyCalTool()
{
  return m_energyCalTool;
}

void InterSpec::showEnergyCalWindow()
{
  if( m_energyCalWindow && !m_toolsTabs )
  {
    m_energyCalWindow->show();
    return;
  }

  const int index = (m_toolsTabs ? m_toolsTabs->indexOf(m_energyCalTool) : -1);
  
  if( index >= 0 )
    m_toolsTabs->removeTab( m_energyCalTool );
  
  if( m_energyCalWindow )
  {
    m_energyCalWindow->stretcher()->removeWidget(m_energyCalTool);
    delete m_energyCalWindow;
  }
    
  m_energyCalWindow = new AuxWindow( "Energy Calibration",
                                WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                    | AuxWindowProperties::TabletNotFullScreen );
  m_energyCalWindow->rejectWhenEscapePressed();
  m_energyCalWindow->stretcher()->addWidget( m_energyCalTool, 0, 0 );
  m_energyCalTool->setTallLayout();
  m_energyCalTool->show();
   
  //m_energyCalWindow->finished().connect(boost::bind( &AuxWindow::deleteAuxWindow, m_energyCalWindow ) );
  m_energyCalWindow->finished().connect( boost::bind( &InterSpec::handEnergyCalWindowClose, this ) );

  if( (m_renderedWidth > 100) && (m_renderedHeight > 100) )
    m_energyCalWindow->setMaximumSize( 0.8*m_renderedWidth, 0.8*m_renderedHeight );
  m_energyCalWindow->setWidth( 380 );
  
  m_energyCalWindow->show();
  m_energyCalWindow->resizeToFitOnScreen();
  m_energyCalWindow->centerWindow();
  
  AuxWindow::addHelpInFooter( m_energyCalWindow->footer(), "energy-calibration" );
  Wt::WPushButton *closeButton = m_energyCalWindow->addCloseButtonToFooter("Close",true);
  closeButton->clicked().connect( boost::bind( &InterSpec::handEnergyCalWindowClose, this ) );
  
  
  if( m_toolsTabs )
    m_currentToolsTab = m_toolsTabs->currentIndex();
}//void showEnergyCalWindow()



void InterSpec::setLogY( bool logy )
{
  InterSpecUser::setPreferenceValue<bool>( m_user, "LogY", logy, this );
  m_logYItems[0]->setHidden( logy );
  m_logYItems[1]->setHidden( !logy );
  m_spectrum->setYAxisLog( logy );
}//void setLogY( bool logy )


void InterSpec::setBackgroundSub( bool subtract )
{
  m_backgroundSubItems[0]->setHidden( subtract );
  m_backgroundSubItems[1]->setHidden( !subtract );
  m_spectrum->setBackgroundSubtract( subtract );
}//void setBackgroundSub( bool sub )


void InterSpec::setVerticalLines( bool show )
{
  InterSpecUser::setPreferenceValue<bool>( m_user, "ShowVerticalGridlines", show, this );
  m_verticalLinesItems[0]->setHidden( show );
  m_verticalLinesItems[1]->setHidden( !show );
  m_spectrum->showVerticalLines( show );
  m_timeSeries->showVerticalLines( show );
}//void setVerticalLines( bool show )


void InterSpec::setHorizantalLines( bool show )
{
  InterSpecUser::setPreferenceValue<bool>( m_user, "ShowHorizontalGridlines", show, this );
  m_horizantalLinesItems[0]->setHidden( show );
  m_horizantalLinesItems[1]->setHidden( !show );
  m_spectrum->showHorizontalLines( show );
  m_timeSeries->showHorizontalLines( show );
}//void setHorizantalLines( bool show )


void InterSpec::startHardBackgroundSub()
{
  const char *msg =
  "<div style=\"text-align: left;\">"
  "<p>The normal background subtraction option only affects display of the data, with operations"
  " like peak-fitting or exporting data are done on the full-statistics original foreground"
  " spectrum.</p>"
  "<p>A &quot;hard background subtraction&quot; creates a modified foreground by doing a bin-by-bin"
  " subtraction of the background from the foreground.</p>"
  "Side effects of doing a hard background subtraction include:"
  "<ul style=\"list-style-type: square; margin-top: 4px;\">"  //list-style-type: none;
    "<li>Variances, i.e. the statistical uncertainty of each channel, will no longer be correct.</li>"
    "<li>Small energy calibration differences between spectra may create artificial features in the data.</li>"
    "<li>If a peak in the foreground overlaps with a peak in the background, the statistical"
         " uncertainty of fit foreground peaks will no longer be correct</li>"
  "</ul>"
  "The primary reasons to choose a hard background subtraction over normal display background subtraction are:"
  "<ul style=\"list-style-type: square; margin-top: 4px;\">"
    "<li>Avoiding artifacts on fit peak continuums.</li>"
    "<li>To simplify background peak subtraction in the <b>Activity/Shielding Fit</b> tool.</li>"
    "<li>You dont care about subtleties this causes (which in practice are minimal).</li>"
  "</ul>"
  "</div>"
  "<br />"
  "<div style=\"text-align: center;\"><b><em>Would you like to proceed?</em></b></div>"
  ;

  
  SimpleDialog *dialog = new SimpleDialog( "Perform Hard Background Subtract?", msg );
  
  // For some reason on Windows Electron version, the dialog does not expand out very wide - so lets
  //  force it
  const int ww = renderedWidth();
  if( ww > 500 )
    dialog->setWidth( std::min(ww/2, 800) );
  
  auto truncate_neg = make_shared<bool>(false);
  auto round_counts = make_shared<bool>(false);
  
  WContainerWidget *optionsDiv = new WContainerWidget( dialog->contents() );
  optionsDiv->setPadding( 40, Wt::Side::Left );
  optionsDiv->setPadding( 20, Wt::Side::Bottom );
  
  WCheckBox *cb = new WCheckBox( "Truncate negative bins at zero", optionsDiv );
  cb->setInline( false );
  cb->checked().connect( std::bind([=](){ *truncate_neg = true; } ) );
  cb->unChecked().connect( std::bind([=](){ *truncate_neg = false; } ) );
  
  cb = new WCheckBox( "Round channel counts to nearest integer", optionsDiv );
  cb->setInline( false );
  cb->checked().connect( std::bind([=](){ *round_counts = true; } ) );
  cb->unChecked().connect( std::bind([=](){ *round_counts = false; } ) );
  
  
  WPushButton *button = dialog->addButton( "Yes" );
  button->setFocus();
  button->clicked().connect( boost::bind( &InterSpec::finishHardBackgroundSub, this, truncate_neg, round_counts ) );
  dialog->addButton( "No" );  //dont need to hook this to anything
}//void startHardBackgroundSub()


void InterSpec::finishHardBackgroundSub( std::shared_ptr<bool> truncate_neg, std::shared_ptr<bool> round_counts )
{
  const auto foreground = m_spectrum->data();
  const auto background = m_spectrum->background();
  const float sf = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::Background);
  
  const bool no_neg = truncate_neg ? *truncate_neg : false;
  const bool do_round = round_counts ? *round_counts : false;
  
  if( !foreground
     || !background
     || !m_dataMeasurement
     || (foreground->num_gamma_channels() < 7)
     || (background->num_gamma_channels() < 7)
     || IsInf(sf) || IsNan(sf) || (sf <= 0.0)
     || !foreground->channel_energies()  //should be covered by num_gamma_channels(), but whatever
     || !background->channel_energies()
     || !background->energy_calibration() //should always be true, but whatever
     || !foreground->energy_calibration()
     || !background->energy_calibration()->valid()
     || !foreground->energy_calibration()->valid()
     )
  {
    passMessage( "Error doing hard background subtraction."
                 " Foreground or background was not available, or background scale factor invalid.",
                 "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  try
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> orig_peaks = m_peakModel->peaks();
    shared_ptr<const vector<float>> fore_counts = foreground->gamma_counts();
    shared_ptr<const vector<float>> back_counts = background->gamma_counts();
    
    // Make sure back_counts has the same energy calibration and fore_counts, so we can subtract
    //  on a bin-by-bin basis
    if( background->energy_calibration() != foreground->energy_calibration()
       && (*background->energy_calibration()) != (*foreground->energy_calibration()) )
    {
      auto new_backchan = make_shared<vector<float>>( fore_counts->size(), 0.0f );
      SpecUtils::rebin_by_lower_edge( *background->channel_energies(), *back_counts,
                                     *foreground->channel_energies(), *new_backchan );
      back_counts = new_backchan;
    }
    
    // Create what will be the background subtracted foreground
    auto back_sub_counts = make_shared<vector<float>>( *fore_counts );
    
    //back_counts and fore_counts should always be the same size, but we'll be a bit verbose anyway
    assert( back_counts->size() == fore_counts->size() );
    const size_t nchann = std::min( back_counts->size(), fore_counts->size() );
    
    // Do the actual background subtraction
    for( size_t i = 0; i < nchann; ++i )
    {
      float &val = (*back_sub_counts)[i];
      val -= sf*(*back_counts)[i];
      
      if( no_neg )
        val = std::max( 0.0f, val );
      
      if( do_round )
        val = std::round( val );
    }//for( size_t i = 0; i < nchann; ++i )
    
    // Create a new Measurement object, based on the old foreground
    auto newspec = make_shared<SpecUtils::Measurement>( *foreground );
    newspec->set_gamma_counts( back_sub_counts, foreground->live_time(), foreground->real_time() );
    vector<string> remarks = foreground->remarks();
    remarks.push_back( "This spectrum has been background subtracted in InterSpec" );
    newspec->set_remarks( remarks );
    newspec->set_sample_number( 1 );
    
    // Create a new spectrum file object, and set new background subtracted Measurement as its only
    //  record
    auto newmeas = make_shared<SpecMeas>();
    
    // Copy just the SpecUtils::SpecFile stuff over to 'newmeas' so we dont copy things like
    //  displayed sample numbers, and uneeded peaks and stuff
    static_cast<SpecUtils::SpecFile &>(*newmeas) = static_cast<SpecUtils::SpecFile &>(*m_dataMeasurement);
    
    newmeas->remove_measurements( newmeas->measurements() );
    
    // Need to make sure UUID will get updated.
    newmeas->set_uuid( "" );
    
    // Update filename
    newmeas->set_filename( "bkgsub_" + newmeas->filename() );
    
    // Actually add the measurement
    newmeas->add_measurement( newspec, true );
    
    // Reset all displayed sample numbers and peaks and stuff
    //newmeas->clearInterSpecDisplayStuff();
    
    // Re-fit peaks and set them
    std::vector<PeakDef> refit_peaks;
    if( orig_peaks && orig_peaks->size() )
    {
      try
      {
        vector<PeakDef> input_peaks;
        for( const auto &i : *orig_peaks )
          input_peaks.push_back( *i );
        
        const double lowE = newspec->gamma_energy_min();
        const double upE = newspec->gamma_energy_max();
      
        refit_peaks = fitPeaksInRange( lowE, upE, 0.0, 0.0, 0.0, input_peaks, newspec, {}, true );
        
        std::deque<std::shared_ptr<const PeakDef> > peakdeque;
        for( const auto &p : refit_peaks )
          peakdeque.push_back( std::make_shared<const PeakDef>(p) );
        
        newmeas->setPeaks( peakdeque, {newspec->sample_number()} );
      }catch( std::exception &e )
      {
        cerr << "Unexpected exception re-fitting peaks doing hard background subtract: "
             << e.what() << endl;
      }//try / catch to fit peaks
    }//if( we need to refit peaks )
    
    // Get rid of the previously displayed background if there was one
    setSpectrum( nullptr, {}, SpecUtils::SpectrumType::Background, 0 );
    
    
    auto header = m_fileManager->addFile( newmeas->filename(), newmeas );
    SpectraFileModel *filemodel = m_fileManager->model();
    auto index = filemodel->index( header );
    m_fileManager->displayFile( index.row(), newmeas,
                                SpecUtils::SpectrumType::Foreground,
                                false, false,
                                SpecMeasManager::VariantChecksToDo::None );
  }catch( std::exception &e )
  {
    passMessage( "There was an error loading the newly created spectrum file, sorry:"
                + string(e.what()),
                "", WarningWidget::WarningMsgHigh );
  }//try / catch
  
}//finishHardBackgroundSub();



void InterSpec::setXAxisSlider( bool show )
{
  InterSpecUser::setPreferenceValue<bool>( m_user, "ShowXAxisSlider", show, this );
  m_showXAxisSliderItems[0]->setHidden( show );
  m_showXAxisSliderItems[1]->setHidden( !show );
  
  if( show )
  {
    //Default to compact x-axis.
    if( m_compactXAxisItems[0] )
      m_compactXAxisItems[0]->setHidden( true );
    if( m_compactXAxisItems[1] )
      m_compactXAxisItems[1]->setHidden( false );
    
     m_spectrum->setCompactAxis( true );
    m_timeSeries->setCompactAxis( true );
  }else
  {
    //Go back to whatever the user wants/selects
    const bool makeCompact = InterSpecUser::preferenceValue<bool>( "CompactXAxis", this );
    
    if( m_compactXAxisItems[0] )
      m_compactXAxisItems[0]->setHidden( makeCompact );
    if( m_compactXAxisItems[1] )
      m_compactXAxisItems[1]->setHidden( !makeCompact );
    
    m_spectrum->setCompactAxis( makeCompact );
    m_timeSeries->setCompactAxis( makeCompact );
  }//show /hide
  
  m_spectrum->showXAxisSliderChart( show );
}//void setXAxisSlider( bool show )


void InterSpec::setXAxisCompact( bool compact )
{
  InterSpecUser::setPreferenceValue<bool>( m_user, "CompactXAxis", compact, this );
  
  if( m_compactXAxisItems[0] )
    m_compactXAxisItems[0]->setHidden( compact );
  if( m_compactXAxisItems[1] )
    m_compactXAxisItems[1]->setHidden( !compact );
  
  m_spectrum->setCompactAxis( compact );
  m_timeSeries->setCompactAxis( compact );
}//void setXAxisCompact( bool compact )


void InterSpec::setShowYAxisScalers( bool show )
{
  const bool hasSecond = (m_secondDataMeasurement || m_backgroundMeasurement);
  
  assert( m_showYAxisScalerItems[0] );
  assert( m_showYAxisScalerItems[1] );
  
  m_showYAxisScalerItems[0]->setHidden( show );
  m_showYAxisScalerItems[0]->setDisabled( !show && !hasSecond );
  
  m_showYAxisScalerItems[1]->setHidden( !show );
  m_showYAxisScalerItems[1]->setDisabled( show && !hasSecond );
  
  m_spectrum->showYAxisScalers( show );
  
  try
  {
    InterSpecUser::setPreferenceValue<bool>( m_user, "ShowYAxisScalers", show, this );
  }catch( std::exception &e )
  {
    cerr << "InterSpec::setShowYAxisScalers: Got exception setting pref: " << e.what() << endl;
  }
}//void setShowYAxisScalers( bool show )



ReferencePhotopeakDisplay *InterSpec::referenceLinesWidget()
{
  return m_referencePhotopeakLines;
}

#if( defined(WIN32) && BUILD_AS_ELECTRON_APP )
  //When users drag files from Outlook on windows into the app
  //  you can call the following functions
void InterSpec::dragEventWithFileContentsStarted()
{
  //Set JS variable to indicate that this isnt a normal file drop on the browser
  doJavaScript( "$('.Wt-domRoot').data('DropFileContents',true);" );
}


void InterSpec::dragEventWithFileContentsFinished()
{ 
  doJavaScript( "$('.Wt-domRoot').data('DropFileContents',false);"
	            "$('#Uploader').remove();"
				"$('.Wt-domRoot').data('IsDragging',false);"); 
}

#endif ///#if( defined(WIN32) && BUILD_AS_ELECTRON_APP )


void InterSpec::addPeakLabelSubMenu( PopupDivMenu *parentWidget )
{
  PopupDivMenu *menu = parentWidget->addPopupMenuItem( "Peak Labels",  "InterSpec_resources/images/tag.svg" );
  
  
//  PopupDivMenuItem *item = menu->addMenuItem( "Show User Labels", "", false );
  //
  WCheckBox *cb = new WCheckBox( "Show User Labels" );
  cb->setChecked(false);
  PopupDivMenuItem *item = menu->addWidget( cb );
  
  cb->checked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakUserLabel, true
  ) );
  cb->unChecked().connect( boost::bind( &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakUserLabel, false
  ) );
  
  cb = new WCheckBox( "Show Peak Energies" );
  cb->setChecked(false);
  item = menu->addWidget( cb );
  
  cb->checked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakEnergyLabel, true
  ) );
  cb->unChecked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakEnergyLabel, false
  ) );
  
  cb = new WCheckBox( "Show Nuclide Names" );
  cb->setChecked(false);
  item = menu->addWidget( cb );
  
  cb->checked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakNuclideLabel, true
  ) );
  cb->unChecked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakNuclideLabel, false
  ) );
  
  cb = new WCheckBox( "Show Nuclide Energies" );
  cb->setChecked(false);
  item = menu->addWidget( cb );
  
  cb->checked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakNuclideEnergies, true
  ) );
  
  cb->unChecked().connect( boost::bind(
          &D3SpectrumDisplayDiv::setShowPeakLabel,
          m_spectrum, SpectrumChart::kShowPeakNuclideEnergies,  false
  ) );
}//void addPeakLabelMenu( Wt::WContainerWidget *menuDiv )



void InterSpec::addAboutMenu( Wt::WWidget *parent )
{
  if( m_helpMenuPopup )
    return;
  
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addAboutMenu(): parent passed in"
                         " must be a PopupDivMenu  or WContainerWidget" );

  if( menuDiv )
  {
    WPushButton *button = new WPushButton( "Help", menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_helpMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_helpMenuPopup = parentMenu->addPopupMenuItem( "Help" );
  }//if( menuDiv ) / else

  PopupDivMenuItem *item = m_helpMenuPopup->addMenuItem( "Welcome..." );
  item->triggered().connect( boost::bind( &InterSpec::showWelcomeDialog, this, true ) );

  m_helpMenuPopup->addSeparator();
  
  item = m_helpMenuPopup->addMenuItem( "Help Contents..." ,  "InterSpec_resources/images/help_minimal.svg");
  
  item->triggered().connect( boost::bind( &HelpSystem::createHelpWindow, string("getting-started") ) );

  Wt::WMenuItem *notifications = m_helpMenuPopup->addMenuItem( "Notification Logs..." , "InterSpec_resources/images/log_file_small.png");
  notifications->triggered().connect( this, &InterSpec::showWarningsWindow );

  m_helpMenuPopup->addSeparator();
  PopupDivMenu *subPopup = m_helpMenuPopup->addPopupMenuItem( "Options", "InterSpec_resources/images/cog_small.png" );
    
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  
  const bool autoStore = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", this );
  WCheckBox *cb = new WCheckBox( " Automatically store session" );
  cb->setChecked( autoStore );
  item = subPopup->addWidget( cb );
  HelpSystem::attachToolTipOn( item, "Automatically stores app state", showToolTips );
  InterSpecUser::associateWidget( m_user, "AutoSaveSpectraToDb", cb, this, false );
  

  if( !isMobile() )
  {
    WCheckBox *checkbox = new WCheckBox( " Show tooltips" );
    item = subPopup->addWidget( checkbox );
    checkbox->setChecked( showToolTips );
    HelpSystem::attachToolTipOn( item,
                                "Show tooltips after hovering over an element for half a second."
                                , true, HelpSystem::ToolTipPosition::Right );
    checkbox->checked().connect( boost::bind( &InterSpec::toggleToolTip, this, true ) );
    checkbox->unChecked().connect( boost::bind( &InterSpec::toggleToolTip, this, false ) );
    InterSpecUser::associateWidget( m_user, "ShowTooltips", checkbox, this, false );
  }//if( !isMobile() )
  
  {//begin add "AskPropagatePeaks" to menu
    const bool doPropogate = InterSpecUser::preferenceValue<bool>( "AskPropagatePeaks", this );
    WCheckBox *checkbox = new WCheckBox( " Ask to Propagate Peaks" );
    checkbox->setChecked( doPropogate );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item,
                                 "When loading spectra from the same detector,"
                                 " ask if should re-fit the same peaks for the"
                                 " new spectrum.  Only applies if new spectrum"
                                 " doesnt have any peaks, but previous"
                                 " foreground did.",
                                 true, HelpSystem::ToolTipPosition::Right );
    checkbox->checked().connect( boost::bind( &InterSpec::toggleToolTip, this, true ) );
    checkbox->unChecked().connect( boost::bind( &InterSpec::toggleToolTip, this, false ) );
    InterSpecUser::associateWidget( m_user, "AskPropagatePeaks", checkbox, this, false );
  }//end add "AskPropagatePeaks" to menu
  
  
  {//begin add "DisplayBecquerel"
    WCheckBox *checkbox = new WCheckBox( " Display in Becquerel" );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, "Display activity in units of becquerel, rather than curie.",
                                 true, HelpSystem::ToolTipPosition::Right );
    InterSpecUser::associateWidget( m_user, "DisplayBecquerel", checkbox, this, false );
  }//end add "DisplayBecquerel"
  
  {//begin add "LoadDefaultDrf"
    WCheckBox *checkbox = new WCheckBox( " Use default DRFs" );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, "When a spectrum is loaded, and you have not previously"
                                " associated a detector response function with the model of"
                                " detector the spectrum is from, whether a default DRF should try"
                                " to be found and loaded automatically.",
                                true, HelpSystem::ToolTipPosition::Right );
    InterSpecUser::associateWidget( m_user, "LoadDefaultDrf", checkbox, this, false );
  }//end add "LoadDefaultDrf"
  
	item = subPopup->addMenuItem("Color Themes...");
	item->triggered().connect(boost::bind(&InterSpec::showColorThemeWindow, this));

#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
  subPopup->addSeparator();
  WCheckBox *promptOnLoad = new WCheckBox( "Prompt to load prev state" );
  item = subPopup->addWidget( promptOnLoad );
  const char *prompttext = "At application start, ask to load previous state.";
  HelpSystem::attachToolTipOn( item, prompttext, showToolTips );
  InterSpecUser::associateWidget( m_user, "PromptStateLoadOnStart", promptOnLoad, this, false );
  
  WCheckBox *doLoad = new WCheckBox( "Load prev state on start" );
  item = subPopup->addWidget( doLoad );
  const char *doloadtext = "At application start, automatically load previous state, if not set to be prompted";
  HelpSystem::attachToolTipOn( item, doloadtext, showToolTips );
  InterSpecUser::associateWidget( m_user, "LoadPrevStateOnStart", doLoad, this, false );
#endif
  
  
#if( BUILD_AS_OSX_APP )
  const bool addAboutInterSpec = !InterSpecApp::isPrimaryWindowInstance();
#else
    const bool addAboutInterSpec = true;
#endif
  
  if( addAboutInterSpec )
  {
    m_helpMenuPopup->addSeparator();
    
    item = m_helpMenuPopup->addMenuItem( "About InterSpec..." );
    item->triggered().connect( boost::bind( &InterSpec::showLicenseAndDisclaimersWindow, this, false, std::function<void()>{} ) );
  }

}//void addAboutMenu( Wt::WContainerWidget *menuDiv )


void InterSpec::toggleToolTip( const bool showToolTips )
{
  //update all existing qtips
  if( showToolTips )
  {
    wApp->doJavaScript( "$('.qtip-rounded.canDisableTt').qtip('option', 'show.event', 'mouseenter focus');" );
  }else
  {
    wApp->doJavaScript( "$('.qtip-rounded.canDisableTt').qtip('option', 'show.event', '');" );
  }
  
}//void toggleToolTip( const bool sticky )



const std::set<int> &InterSpec::displayedSamples( SpecUtils::SpectrumType type ) const
{
  static const std::set<int> empty;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
    {
      if( !m_dataMeasurement )
        return empty;
      return m_displayedSamples;
    }//case SpecUtils::SpectrumType::Foreground:

    case SpecUtils::SpectrumType::SecondForeground:
    {
      if( !m_secondDataMeasurement )
        return empty;
      return m_sectondForgroundSampleNumbers;
    }//case SpecUtils::SpectrumType::SecondForeground:

    case SpecUtils::SpectrumType::Background:
    {
      if( !m_backgroundMeasurement )
        return empty;
      return m_backgroundSampleNumbers;
    }//case SpecUtils::SpectrumType::Background:
  }//switch( type )

  throw runtime_error( "InterSpec::displayedSamples(...) - Serious Badness" );

  return empty;  //keep compiler from complaining
}//const std::set<int> &displayedSamples( SpecUtils::SpectrumType spectrum_type ) const


std::shared_ptr<const SpecMeas> InterSpec::measurment( SpecUtils::SpectrumType type ) const
{
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return m_dataMeasurement;
    case SpecUtils::SpectrumType::SecondForeground:
      return m_secondDataMeasurement;
    case SpecUtils::SpectrumType::Background:
      return m_backgroundMeasurement;
  }//switch( type )

  return std::shared_ptr<const SpecMeas>();
}//measurment(...)


std::shared_ptr<SpecMeas> InterSpec::measurment( SpecUtils::SpectrumType type )
{
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return m_dataMeasurement;
    case SpecUtils::SpectrumType::SecondForeground:
      return m_secondDataMeasurement;
    case SpecUtils::SpectrumType::Background:
      return m_backgroundMeasurement;
  }//switch( type )

  return std::shared_ptr<SpecMeas>();
}//std::shared_ptr<SpecMeas> measurment( SpecUtils::SpectrumType spectrum_type )


#if( USE_DB_TO_STORE_SPECTRA )
Wt::Dbo::ptr<UserFileInDb> InterSpec::measurmentFromDb( SpecUtils::SpectrumType type,
                                                             bool update )
{
  try
  {
    Wt::Dbo::ptr<UserFileInDb> answer;
    std::shared_ptr<SpectraFileHeader> header;
  
    SpectraFileModel *fileModel = m_fileManager->model();
    std::shared_ptr<SpecMeas> meas = measurment( type );
    if( !meas )
      return answer;
  
    WModelIndex index = fileModel->index( meas );
    if( !index.isValid() )
      return answer;
    
    header = fileModel->fileHeader( index.row() );
    if( !header )
      return answer;
  
    answer = header->dbEntry();
    if( answer && !meas->modified() )
      return answer;
  
    bool savePref = false;
    try{ savePref = m_user->preferenceValue<bool>( "AutoSaveSpectraToDb" ); }catch(...){}
    
    if( !savePref )
      return answer;
    
    Dbo::ptr<UserFileInDb> dbback;
    
    if( type == SpecUtils::SpectrumType::Foreground )
    {
      WModelIndex bindex;
      std::shared_ptr<SpectraFileHeader> bheader;
      std::shared_ptr<SpecMeas> background = measurment( SpecUtils::SpectrumType::Background );
      if( background )
         bindex = fileModel->index( background );
      if( bindex.isValid() )
        bheader = fileModel->fileHeader( bindex.row() );
      if( bheader )
      {
        dbback = bheader->dbEntry();
        if( !dbback )
        {
          try
          {
            bheader->saveToDatabase( background );
            dbback = bheader->dbEntry();
          }catch(...){}
        }
      }//if( bheader )
    }//if( type == SpecUtils::SpectrumType::Foreground )
    
    if( update && savePref )
      header->saveToDatabase( meas );
    
    return header->dbEntry();
  }catch( std::exception &e )
  {
    cerr << "\n\nSpectrumViewer::measurmentFromDb(...) caught: " << e.what()
         << endl;
  }//try / catch
  
  return Wt::Dbo::ptr<UserFileInDb>();
}//Wt::Dbo::ptr<UserFileInDb> measurmentFromDb( SpecUtils::SpectrumType type, bool update );
#endif //#if( USE_DB_TO_STORE_SPECTRA )

std::shared_ptr<const SpecUtils::Measurement> InterSpec::displayedHistogram( SpecUtils::SpectrumType spectrum_type ) const
{
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return m_spectrum->data();
    case SpecUtils::SpectrumType::SecondForeground:
      return m_spectrum->secondData();
    case SpecUtils::SpectrumType::Background:
      return m_spectrum->background();
//  m_spectrum->continuum();
  }//switch( spectrum_type )

  throw runtime_error( "InterSpec::displayedHistogram(...): invalid input arg" );

  return std::shared_ptr<const SpecUtils::Measurement>();
}//displayedHistogram(...)


void InterSpec::saveChartToImg( const bool spectrum, const bool asPng )
{
  std::shared_ptr<const SpecMeas> spec = measurment(SpecUtils::SpectrumType::Foreground);
  string filename = (spec ? spec->filename() : string("spectrum"));
  const string ext = SpecUtils::file_extension(filename);
  if( !ext.empty() && (ext.size() <= filename.size()) )
    filename = filename.substr(0,filename.size()-ext.size());
  if( filename.empty() )
    filename = "spectrum";
  if( !spectrum )
    filename += "_timechart";
  std::string timestr = SpecUtils::to_iso_string( boost::posix_time::second_clock::local_time() );
  auto ppos = timestr.find('.');
  if( ppos != string::npos )
    timestr = timestr.substr(0,ppos);
  filename += "_" + timestr + ((!spectrum || asPng) ? ".png" : ".svg");
  
  string illegal_chars = "\\/:?\"<>|";
  SpecUtils::erase_any_character( filename, illegal_chars.c_str() );
  
  if( spectrum )
  {
    m_spectrum->saveChartToImg( filename, asPng );
  }else
  {
    m_timeSeries->saveChartToPng( filename );
  }
}//saveSpectrumToPng()


double InterSpec::displayScaleFactor( SpecUtils::SpectrumType spectrum_type ) const
{
  return m_spectrum->displayScaleFactor( spectrum_type );
}//double displayScaleFactor( SpecUtils::SpectrumType spectrum_type ) const


void InterSpec::setDisplayScaleFactor( const double sf, const SpecUtils::SpectrumType spec_type )
{
  m_spectrum->setDisplayScaleFactor( sf, spec_type );
  m_spectrumScaleFactorChanged.emit( spec_type, sf );
}//void setDisplayScaleFactor( const double sf, SpecUtils::SpectrumType spectrum_type );


float InterSpec::liveTime( SpecUtils::SpectrumType type ) const
{
  if( !measurment(type) )
    return 0.0f;
  
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      return m_spectrum->foregroundLiveTime();
    case SpecUtils::SpectrumType::SecondForeground:
      return m_spectrum->secondForegroundLiveTime();
    case SpecUtils::SpectrumType::Background:
      return m_spectrum->backgroundLiveTime();
  }//switch( type )
  
  return 0.0f;
}//float liveTime( SpecUtils::SpectrumType type ) const


int InterSpec::renderedWidth() const
{
  return m_renderedWidth;
}

int InterSpec::renderedHeight() const
{
  return m_renderedHeight;
}


void InterSpec::createOneOverR2Calculator()
{
//  OneOverR2Calc *calc =
  new OneOverR2Calc();
  
//  if( !toolTabsVisible() )
//  {
//    const int maxHeight = static_cast<int>(0.95*paintedHeight());
//    const int maxWidth = static_cast<int>(0.95*paintedWidth());
//    calc->setMaximumSize( maxWidth, maxHeight );
//    calc->contents()->setOverflow( WContainerWidget::OverflowAuto );
//  }//if( !toolTabsVisible() )
}//void createOneOverR2Calculator()


void InterSpec::createUnitsConverterTool()
{
  new UnitsConverterTool();
}//void createUnitsConverterTool()


void InterSpec::createFluxTool()
{
  new FluxToolWindow( this );
}//void createFluxTool()


void InterSpec::createDecayInfoWindow()
{
  if( !m_decayInfoWindow )
  {
    m_decayInfoWindow = new DecayWindow( this );
    m_decayInfoWindow->finished().connect( boost::bind( &InterSpec::deleteDecayInfoWindow, this ) );
  }
  
  if( m_referencePhotopeakLines )
  {
    const ReferenceLineInfo &nuc
                          = m_referencePhotopeakLines->currentlyShowingNuclide();
    if( nuc.nuclide )
    {
      m_decayInfoWindow->clearAllNuclides();
      
      //\todo We could do a little better and check the Shielding/Source Fit
      //  widget and grab those activities (and ages) if they match
      
      m_decayInfoWindow->addNuclide( nuc.nuclide->atomicNumber,
                         nuc.nuclide->massNumber,
                         nuc.nuclide->isomerNumber,
                         1.0*PhysicalUnits::microCi, true,
                         0.0, 5.0*nuc.age );
    }//if( nuc.nuclide )
  }//if( m_referencePhotopeakLines )
}//void createDecayInfoWindow()


void InterSpec::deleteDecayInfoWindow()
{
  if( m_decayInfoWindow )
    AuxWindow::deleteAuxWindow( m_decayInfoWindow );
  m_decayInfoWindow = nullptr;
}//void deleteDecayInfoWindow()


void InterSpec::createFileParameterWindow()
{
  new SpecFileSummary( this );
}//void createFileParameterWindow()


#if( USE_GOOGLE_MAP )
void InterSpec::displayOnlySamplesWithinView( GoogleMap *map,
                                  const SpecUtils::SpectrumType targetSamples,
                                  const SpecUtils::SpectrumType fromSamples )
{
  float uplat, leftlng, lowerlat, rightlng;
  map->getMapExtent( uplat, leftlng, lowerlat, rightlng );
  
  if( lowerlat > uplat )
    std::swap( lowerlat, uplat );
  if( leftlng > rightlng )
    std::swap( leftlng, rightlng );
  
  std::set<int> sample_numbers;
  std::shared_ptr<SpecMeas> meass = measurment( fromSamples );
  
  if( !meass )
  {
    passMessage( "Could not load spectrum", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !meass )
  
  for( const int sample : meass->sample_numbers() )
  {
    bool samplewithin = false;
    for( const int detnum : meass->detector_numbers() )
    {
      std::shared_ptr<const SpecUtils::Measurement> m = meass->measurement( sample, detnum );
      if( !!m && m->has_gps_info()
         && m->longitude()>=leftlng && m->longitude()<=rightlng
         && m->latitude()>=lowerlat && m->latitude()<=uplat )
      {
        samplewithin = true;
        break;
      }
    }//for( const int detnum : meass->detector_numbers() )
    
    if( samplewithin )
      sample_numbers.insert( sample );
  }//for( int sample : meass->sample_numbers() )
  
  if( sample_numbers.empty() )
  {
    passMessage( "There were no samples in the visible map area.",
                 "", WarningWidget::WarningMsgHigh );
    return;
  }//if( sample_numbers.empty() )
  
  
  if( (fromSamples!=targetSamples) || targetSamples!=SpecUtils::SpectrumType::Foreground )
  {
    setSpectrum( meass, sample_numbers, targetSamples, SetSpectrumOptions::CheckToPreservePreviousEnergyCal );
  }else
  {
    changeDisplayedSampleNums( sample_numbers, targetSamples );
  }
}//displayOnlySamplesWithinView(...)


void InterSpec::createMapWindow( SpecUtils::SpectrumType spectrum_type )
{
  std::shared_ptr<const SpecMeas> meas = measurment( spectrum_type );
  
  if( !meas )
    return;
  
  const set<int> &samples = displayedSamples( spectrum_type );
  
  AuxWindow *window = new AuxWindow( "Map",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::EnableResize)
                                      | AuxWindowProperties::DisableCollapse)
                                    );
  
  int w = 0.66*renderedWidth();
  int h = 0.8*renderedHeight();
  
  const char *label = 0;
  switch( spectrum_type )
  {
    case SpecUtils::SpectrumType::Foreground:       label = "Foreground";        break;
    case SpecUtils::SpectrumType::SecondForeground: label = "Second Foreground"; break;
    case SpecUtils::SpectrumType::Background:       label = "Background";        break;
  }//switch( spectrum_type )
  
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  //window->footer()->setStyleClass( "modal-footer" );
  
  const bool enableLoadingVisible = (meas->sample_numbers().size() > 1);

  GoogleMap *googlemap = new GoogleMap( enableLoadingVisible );
  WGridLayout *layout = window->stretcher();
  googlemap->addMeasurment( meas, label, samples );
  layout->addWidget( googlemap, 0, 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  
  //We need to set the footer height explicitly, or else the window->resize()
  //  messes up.
//  window->footer()->resize( WLength::Auto, WLength(50.0) );
  
  if( enableLoadingVisible )
  {
    WPushButton *button = new WPushButton( "Load Visible Points...", window->footer() );
    WPopupMenu *menu = new WPopupMenu();
    menu->setAutoHide( true );
    button->setMenu( menu );
    WMenuItem *item = menu->addItem( "As Foreground" );
    item->triggered().connect( boost::bind( &InterSpec::displayOnlySamplesWithinView, this, googlemap, SpecUtils::SpectrumType::Foreground, spectrum_type ) );
    item = menu->addItem( "As Background" );
    item->triggered().connect( boost::bind( &InterSpec::displayOnlySamplesWithinView, this, googlemap, SpecUtils::SpectrumType::Background, spectrum_type ) );
    item = menu->addItem( "As Secondary" );
    item->triggered().connect( boost::bind( &InterSpec::displayOnlySamplesWithinView, this, googlemap, SpecUtils::SpectrumType::SecondForeground, spectrum_type ) );
  }//if( meas->measurements().size() > 10 )
  
  
  WPushButton *closeButton = window->addCloseButtonToFooter();
  closeButton->clicked().connect( window, &AuxWindow::hide );
  
  window->resize( WLength(w), WLength(h) );
  window->show();
  window->centerWindow();
  window->rejectWhenEscapePressed();
  
//  window->resizeToFitOnScreen();
}//void createMapWindow()
#endif //#if( USE_GOOGLE_MAP )


#if( USE_SEARCH_MODE_3D_CHART )
void InterSpec::create3DSearchModeChart()
{
  if( !m_dataMeasurement || !m_dataMeasurement->passthrough() )
  {
    passMessage( "The 3D chart is only available for search mode or RPM passthrough data.",
                "", WarningWidget::WarningMsgInfo );
    return;
  }//if( we dont have the proper data to make a 3D chart )
  
  AuxWindow *dialog = new AuxWindow( "3D Data View" );
  dialog->disableCollapse();
  dialog->Wt::WDialog::rejectWhenEscapePressed();
  
  WGridLayout *layout = dialog->stretcher();
  layout->setHorizontalSpacing( 0 );
  layout->setVerticalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  SearchMode3DChart *chart = new SearchMode3DChart( this );
  layout->addWidget( chart, 0, 0 );
  
  dialog->show();
  dialog->setClosable( true );
  dialog->resizeScaledWindow( 0.95, 0.95 );
  dialog->centerWindow();
  dialog->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, dialog ) );
}//void create3DSearchModeChart()
#endif


void InterSpec::showRiidResults( const SpecUtils::SpectrumType type )
{
  showRiidInstrumentsAna( measurment(type) );
}//void showRiidResults( const SpecUtils::SpectrumType type )


#if( USE_TERMINAL_WIDGET )
void InterSpec::createTerminalWidget()
{
  if( m_terminal )
    return;
  
  m_terminal = new TerminalWidget( this );
  m_terminal->focusText();
  
  if( m_toolsTabs )
  {
    WMenuItem *item = m_toolsTabs->addTab( m_terminal, "Terminal" );
    item->setCloseable( true );
    m_toolsTabs->setCurrentWidget( m_terminal );
    const int index = m_toolsTabs->currentIndex();
    m_toolsTabs->setTabToolTip( index, "Numeric, algebraic, and text-based spectrum interaction terminal." );
    m_toolsTabs->tabClosed().connect( this, &InterSpec::handleTerminalWindowClose );
  }else
  {
    m_terminalWindow = new AuxWindow( "Terminal",
                                     (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                      | AuxWindowProperties::EnableResize | AuxWindowProperties::TabletNotFullScreen) );
    
    m_terminalWindow->rejectWhenEscapePressed();
    m_terminalWindow->finished().connect( this, &InterSpec::handleTerminalWindowClose );

    m_terminalWindow->show();
    if( m_renderedWidth > 100 && m_renderedHeight > 100 && !isPhone() )
    {
      m_terminalWindow->resizeWindow( 0.95*m_renderedWidth, 0.25*m_renderedHeight );
      m_terminalWindow->centerWindow();
    }
    
    m_terminalWindow->stretcher()->addWidget( m_terminal, 0, 0 );
  }//if( toolTabsVisible() )
  
  m_terminalMenuItem->disable();
}//void createTerminalWidget()


void InterSpec::handleTerminalWindowClose()
{
  if( !m_terminal )
    return;
 
  m_terminalMenuItem->enable();
  
  if( m_terminalWindow )
  {
    delete m_terminalWindow;
  }else
  {
    delete m_terminal;
    if( m_toolsTabs )
      m_toolsTabs->setCurrentIndex( 2 );
  }
  
  m_terminal = 0;
  m_terminalWindow = 0;
}//void handleTerminalWindowClose()
#endif  //#if( USE_TERMINAL_WIDGET )



void InterSpec::addToolsMenu( Wt::WWidget *parent )
{
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addToolsMenu(): parent passed in"
                        " must be a PopupDivMenu  or WContainerWidget" );
  
  PopupDivMenu *popup = NULL;
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( "Tools", menuDiv );
    button->addStyleClass( "MenuLabel" );
    popup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    popup = parentMenu->addPopupMenuItem( "Tools" );
  }

  m_toolsMenuPopup = popup;
  
  PopupDivMenuItem *item = NULL;

  item = popup->addMenuItem( "Activity/Shielding Fit" );
  HelpSystem::attachToolTipOn( item,"Allows advanced input of shielding material and activity around source isotopes to improve the fit." , showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showShieldingSourceFitWindow, this ) );
  
  item = popup->addMenuItem( "Gamma XS Calc", "" );
  HelpSystem::attachToolTipOn( item,"Allows user to determine the cross section for gammas of arbitrary energy though any material in <code>InterSpec</code>'s library. Efficiency estimates for detection of the gamma rays inside the full energy peak and the fraction of gamma rays that will make it through the material without interacting with it can be provided with the input of additional information.", showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showGammaXsTool, this ) );
    
    
  item = popup->addMenuItem( "Dose Calc", "" );
  HelpSystem::attachToolTipOn( item,
      "Allows you to compute dose, activity, shielding, or distance, given the"
      " other pieces of information.", showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showDoseTool, this ) );
  
//  item = popup->addMenuItem( Wt::WString::fromUTF8("1/r² Calculator") );  // is superscript 2
#if( USE_OSX_NATIVE_MENU  || USING_ELECTRON_NATIVE_MENU )
  item = popup->addMenuItem( Wt::WString::fromUTF8("1/r\x0032 Calculator") );  //works on OS X at least.
#else
  item = popup->addMenuItem( Wt::WString::fromUTF8("1/r<sup>2</sup> Calculator") );
  item->makeTextXHTML();
#endif
  
  HelpSystem::attachToolTipOn( item,"Allows user to use two dose measurements taken at different distances from a source to determine the absolute distance to the source from the nearer measurement.", showToolTips );
  
  item->triggered().connect( this, &InterSpec::createOneOverR2Calculator );

  item = popup->addMenuItem( "Units Converter" );
  HelpSystem::attachToolTipOn( item, "Convert radiation-related units.", showToolTips );
  item->triggered().connect( this, &InterSpec::createUnitsConverterTool );
  
  item = popup->addMenuItem( "Flux Tool" );
  HelpSystem::attachToolTipOn( item,"Converts detected peak counts to gammas emitted by the source.", showToolTips );
  item->triggered().connect( this, &InterSpec::createFluxTool );
  
  item = popup->addMenuItem( "Nuclide Decay Info" );
  HelpSystem::attachToolTipOn( item,"Allows user to obtain advanced information about activities, gamma/alpha/beta production rates, decay chain, and daughter nuclides." , showToolTips );
  item->triggered().connect( this, &InterSpec::createDecayInfoWindow );

  
  item = popup->addMenuItem( "Detector Response Select" );
  HelpSystem::attachToolTipOn( item,"Allows user to change the detector response function.", showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showDrfSelectWindow, this ) );
  
  item = popup->addMenuItem( "Make Detector Response" );
  HelpSystem::attachToolTipOn( item, "Create detector response function from characterization data.", showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showMakeDrfWindow, this ) );

  
  item = popup->addMenuItem( "File Parameters" );
  HelpSystem::attachToolTipOn( item,"Allows user to view/edit the file parameters. If ever the application is unable to render activity calculation, use this tool to provide parameters the original file did not provide; <code>InterSpec</code> needs all parameters for activity calculation.", showToolTips );
  item->triggered().connect( this, &InterSpec::createFileParameterWindow );

  item = popup->addMenuItem( "Energy Range Sum" );
  HelpSystem::attachToolTipOn( item, "Sums the number of gammas in region of interest (ROI). Can also be accessed by left-click dragging over the ROI while holding both the <kbd><b>ALT</b></kbd> and <kbd><b>SHIFT</b></kbd> keys.", showToolTips );
  item->triggered().connect( this, &InterSpec::showGammaCountDialog );
  
#if( USE_SPECRUM_FILE_QUERY_WIDGET )
  
#if( BUILD_AS_OSX_APP )
  const bool addQueryTool = InterSpecApp::isPrimaryWindowInstance();
#else
  const bool addQueryTool = true;
#endif
  
  if( addQueryTool )
  {
    item = popup->addMenuItem( "File Query Tool" );
    HelpSystem::attachToolTipOn( item, "Searches through a directory (recursively) for spectrum files that match specafiable conditions.", showToolTips );
    item->triggered().connect( this, &InterSpec::showFileQueryDialog );
  }//if( addQueryTool )
#endif
  
#if( USE_TERMINAL_WIDGET )
  m_terminalMenuItem = popup->addMenuItem( "Math/Command Terminal" );
  HelpSystem::attachToolTipOn( m_terminalMenuItem, "Creates a terminal that provides numeric and algebraic computations, as well as allowing text based interactions with the spectra.", showToolTips );
  m_terminalMenuItem->triggered().connect( this, &InterSpec::createTerminalWidget );
#endif
}//void InterSpec::addToolsMenu( Wt::WContainerWidget *menuDiv )



void InterSpec::fillMaterialDb( std::shared_ptr<MaterialDB> materialDB,
                                     const std::string sessionid,
                                     boost::function<void(void)> update )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  try
  {
    //materialDB can get destructed if the session ends immediately....
    const string materialfile = SpecUtils::append_path(sm_staticDataDirectory, "MaterialDataBase.txt" );
    materialDB->parseGadrasMaterialFile( materialfile, db, false );
    
    WServer::instance()->post( sessionid, update );
  }catch( std::exception &e )
  {
    WString msg = "Error initializing the material database: " + string(e.what());
    
    WServer::instance()->post( sessionid, boost::bind( &postSvlogHelper, msg, int(WarningWidget::WarningMsgHigh) ) );
    
    return;
  }//try / catch
}//void fillMaterialDb(...)


void InterSpec::pushMaterialSuggestionsToUsers()
{
  if( !m_materialDB || !m_shieldingSuggestion )
    throw runtime_error( "pushMaterialSuggestionsToUsers(): you must"
                        " call initMaterialDbAndSuggestions() first." );
  
  for( const string &name : m_materialDB->names() )
    m_shieldingSuggestion->addSuggestion( name, name );
  
  wApp->triggerUpdate();
}//void pushMaterialSuggestionsToUsers()



void InterSpec::initMaterialDbAndSuggestions()
{
  if( !m_shieldingSuggestion )
  {
    WSuggestionPopup::Options popupOptions;
    popupOptions.highlightBeginTag  = "<b>";          //Open tag to highlight a match in a suggestion.
    popupOptions.highlightEndTag    = "</b>";         //Close tag to highlight a match in a suggestion.
    popupOptions.listSeparator      = ',';            //(char) When editing a list of values, the separator used for different items.
    popupOptions.whitespace         = " \\t()";       //When editing a value, the whitespace characters ignored before the current value.
    popupOptions.wordSeparators     = "-_., ;()";     //To show suggestions based on matches of the edited value with parts of the suggestion.
    popupOptions.appendReplacedText = "";             //When replacing the curr
    m_shieldingSuggestion = new WSuggestionPopup( popupOptions );
    m_shieldingSuggestion->addStyleClass("suggestion");
    m_shieldingSuggestion->setJavaScriptMember("wtNoReparent", "true");
    m_shieldingSuggestion->setFilterLength( 0 );
    m_shieldingSuggestion->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  }//if( !m_shieldingSuggestion )

  if( !m_materialDB )
  {
    m_materialDB = std::make_shared<MaterialDB>();
    
    boost::function<void(void)> success = wApp->bind( boost::bind(&InterSpec::pushMaterialSuggestionsToUsers, this) );
    
    boost::function<void(void)> worker = boost::bind( &fillMaterialDb,
                                    m_materialDB, wApp->sessionId(), success );
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
  }//if( !m_materialDB )
}//void InterSpec::initMaterialDbAndSuggestions()


void InterSpec::showGammaXsTool()
{
  new GammaXsWindow( m_materialDB.get(), m_shieldingSuggestion, this );
} //showGammaXsTool()


void InterSpec::showDoseTool()
{
  new DoseCalcWindow( m_materialDB.get(), m_shieldingSuggestion, this );
}


void InterSpec::showMakeDrfWindow()
{
  MakeDrf::makeDrfWindow( this, m_materialDB.get(), m_shieldingSuggestion );
}//void showDrfSelectWindow()


void InterSpec::showDrfSelectWindow()
{
  std::shared_ptr<DetectorPeakResponse> currentDet;
  if( m_dataMeasurement )
    currentDet = m_dataMeasurement->detector();
  InterSpec *specViewer = this;
  SpectraFileModel *fileModel = m_fileManager->model();

  new DrfSelectWindow( currentDet, specViewer, fileModel );
}//void showDrfSelectWindow()


void InterSpec::showCompactFileManagerWindow()
{
 auto *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::Tabbed );

  m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                 boost::placeholders::_1, boost::placeholders::_2 ) );
  
  AuxWindow *window = new AuxWindow( "Select Opened Spectra to Display", (AuxWindowProperties::TabletNotFullScreen) );
  window->disableCollapse();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  
  WPushButton *closeButton = window->addCloseButtonToFooter();
  closeButton->clicked().connect(window, &AuxWindow::hide);
  
  WGridLayout *layout = window->stretcher();
  layout->addWidget( compact, 0, 0 );
  
  window->show();
//  window->resizeToFitOnScreen();
  window->centerWindow();
}//void showCompactFileManagerWindow()

void InterSpec::closeNuclideSearchWindow()
{
  if( !m_nuclideSearchWindow )
    return;
  
  m_nuclideSearch->clearSearchEnergiesOnClient();
  m_nuclideSearchWindow->stretcher()->removeWidget( m_nuclideSearch );
  
  delete m_nuclideSearchWindow;
  m_nuclideSearchWindow = 0;
  
  if( m_toolsTabs )
  {
    m_nuclideSearchContainer = new WContainerWidget();
    WGridLayout *isotopeSearchGridLayout = new WGridLayout();
    m_nuclideSearchContainer->setLayout( isotopeSearchGridLayout );

    isotopeSearchGridLayout->addWidget( m_nuclideSearch, 0, 0 );
    isotopeSearchGridLayout->setRowStretch( 0, 1 );
    isotopeSearchGridLayout->setColumnStretch( 0, 1 );

    m_toolsTabs->addTab( m_nuclideSearchContainer, NuclideSearchTabTitle, TabLoadPolicy );
    m_currentToolsTab = m_toolsTabs->currentIndex();
  }//if( m_toolsTabs )
}//void closeNuclideSearchWindow()

void InterSpec::showNuclideSearchWindow()
{
  if( m_nuclideSearchWindow )
  {
    m_nuclideSearchWindow->show();
    m_nuclideSearchWindow->resizeToFitOnScreen();
    m_nuclideSearchWindow->centerWindow();
    m_nuclideSearch->loadSearchEnergiesToClient();
    return;
  }//if( m_nuclideSearchWindow )
  
  if( m_toolsTabs && m_nuclideSearchContainer )
  {
    m_nuclideSearchContainer->layout()->removeWidget( m_nuclideSearch );
    m_toolsTabs->removeTab( m_nuclideSearchContainer );
    delete m_nuclideSearchContainer;
    m_nuclideSearchContainer = 0;
  }
  
  
  m_nuclideSearchWindow = new AuxWindow( NuclideSearchTabTitle,
                                        (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                         | AuxWindowProperties::EnableResize) );
  m_nuclideSearchWindow->contents()->setOverflow(Wt::WContainerWidget::OverflowHidden);
  m_nuclideSearchWindow->finished().connect( boost::bind( &InterSpec::closeNuclideSearchWindow, this ) );
  m_nuclideSearchWindow->rejectWhenEscapePressed();
  
  //if( isMobile() )
  //{
  //  m_nuclideSearchWindow->contents()->setPadding( 0 );
  //  m_nuclideSearchWindow->contents()->setMargin( 0 );
  //}//if( isPhone() )
  
  m_nuclideSearchWindow->stretcher()->setContentsMargins( 0, 0, 0, 0 );
  m_nuclideSearchWindow->stretcher()->addWidget( m_nuclideSearch, 0, 0 );
  
  //We need to set the footer height explicitly, or else the window->resize()
  //  messes up.
//  m_nuclideSearchWindow->footer()->resize( WLength::Auto, WLength(50.0) );
  
 
  Wt::WPushButton *closeButton = m_nuclideSearchWindow->addCloseButtonToFooter("Close",true);
  
  closeButton->clicked().connect( boost::bind( &InterSpec::closeNuclideSearchWindow, this ) );
  
  AuxWindow::addHelpInFooter( m_nuclideSearchWindow->footer(), "nuclide-search-dialog" );
  
  m_nuclideSearchWindow->resize( WLength(800,WLength::Pixel), WLength(510,WLength::Pixel));
  m_nuclideSearchWindow->resizeToFitOnScreen();
  m_nuclideSearchWindow->centerWindow();
  m_nuclideSearchWindow->show();
  
  m_nuclideSearch->loadSearchEnergiesToClient(); //clear the isotope search on the canvas
  
  if( m_toolsTabs )
    m_currentToolsTab = m_toolsTabs->currentIndex();
}//void showNuclideSearchWindow()


void InterSpec::showShieldingSourceFitWindow()
{
  if( !m_shieldingSourceFit )
  {
    assert( m_peakInfoDisplay );
    auto widgets = ShieldingSourceDisplay::createWindow( this );
    
    m_shieldingSourceFit = widgets.first;
    m_shieldingSourceFitWindow  = widgets.second;
  }else
  {
    const double windowWidth = 0.95 * renderedWidth();
    const double windowHeight = 0.95 * renderedHeight();
    m_shieldingSourceFitWindow->resizeWindow( windowWidth, windowHeight );
    
    m_shieldingSourceFitWindow->resizeToFitOnScreen();
    m_shieldingSourceFitWindow->show();
    m_shieldingSourceFitWindow->centerWindow();
  }//if( !m_shieldingSourceFit )
}//void showShieldingSourceFitWindow()


void InterSpec::saveShieldingSourceModelToForegroundSpecMeas()
{
  if( !m_shieldingSourceFitWindow || !m_shieldingSourceFit || !m_dataMeasurement )
    return;
  
  string xml_data;
  std::unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
  
  m_shieldingSourceFit->serialize( doc.get() );
  
  //rapidxml::print(std::back_inserter(xml_data), *doc, 0);
  //cout << "\n\nsaveShieldingSourceModelToForegroundSpecMeas Model: " << xml_data << endl << endl;
  
  m_dataMeasurement->setShieldingSourceModel( std::move(doc) );
}//void saveShieldingSourceModelToForegroundSpecMeas()

void InterSpec::closeShieldingSourceFitWindow()
{
  if( !m_shieldingSourceFitWindow || !m_shieldingSourceFit )
    return;
  
#if( USE_DB_TO_STORE_SPECTRA )
  m_shieldingSourceFit->saveModelIfAlreadyInDatabase();
#endif
  
  saveShieldingSourceModelToForegroundSpecMeas();
  
  delete m_shieldingSourceFitWindow;
  m_shieldingSourceFitWindow = nullptr;
  m_shieldingSourceFit = nullptr;
}//void closeShieldingSourceFitWindow()


void InterSpec::showGammaLinesWindow()
{
  if( m_referencePhotopeakLinesWindow )
  {
    m_referencePhotopeakLinesWindow->show();
    return;
  }

  if( m_toolsTabs && m_referencePhotopeakLines )
    m_toolsTabs->removeTab( m_referencePhotopeakLines );

  
  std::string xml_state;
  
  if( m_referencePhotopeakLines )
  {
    if( !m_referencePhotopeakLines->currentlyShowingNuclide().empty()
        || m_referencePhotopeakLines->persistedNuclides().size() )
    m_referencePhotopeakLines->serialize( xml_state );
    
    m_referencePhotopeakLines->clearAllLines();
    delete m_referencePhotopeakLines;
    m_referencePhotopeakLines = NULL;
  }//if( m_referencePhotopeakLines )

  m_referencePhotopeakLinesWindow = new AuxWindow( GammaLinesTabTitle,
                                                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                                   | AuxWindowProperties::EnableResize)
                                                  );
  m_referencePhotopeakLinesWindow->contents()->setOverflow(WContainerWidget::OverflowHidden);
  m_referencePhotopeakLinesWindow->rejectWhenEscapePressed();

  m_referencePhotopeakLines = new ReferencePhotopeakDisplay( m_spectrum,
                                               m_materialDB.get(),
                                               m_shieldingSuggestion,
                                               this );
  setReferenceLineColors( nullptr );
  
  if( xml_state.size() )
    m_referencePhotopeakLines->deSerialize( xml_state );

  Wt::WGridLayout *layout = new Wt::WGridLayout();
  layout->setContentsMargins(5,5,5,5);
  m_referencePhotopeakLinesWindow->contents()->setLayout(layout);
  layout->addWidget( m_referencePhotopeakLines, 0, 0 );

  Wt::WPushButton *closeButton = m_referencePhotopeakLinesWindow->addCloseButtonToFooter("Close",true);
  
  if( isPhone() )
  {
    m_referencePhotopeakLines->displayingNuclide().connect( boost::bind( &WPushButton::setText, closeButton, WString("Show Lines")) );
    m_referencePhotopeakLines->nuclidesCleared().connect( boost::bind( &WPushButton::setText, closeButton, WString("Close")) );
  }
  
  closeButton->clicked().connect( m_referencePhotopeakLinesWindow, &AuxWindow::hide );
  m_referencePhotopeakLinesWindow->finished().connect( boost::bind( &InterSpec::closeGammaLinesWindow, this ) );
  
  AuxWindow::addHelpInFooter( m_referencePhotopeakLinesWindow->footer(),
                              "reference-gamma-lines-dialog" );
  
  double w = renderedWidth();
  if( w < 100.0 )
    w = 800.0;
  w = std::min( w, 800.0 );
  
  m_referencePhotopeakLinesWindow->resize( WLength(w,WLength::Pixel), WLength(310,WLength::Pixel));
  m_referencePhotopeakLinesWindow->resizeToFitOnScreen();
  m_referencePhotopeakLinesWindow->centerWindow();
  m_referencePhotopeakLinesWindow->show();
  
  if( m_toolsTabs )
    m_currentToolsTab = m_toolsTabs->currentIndex();
}//void showGammaLinesWindow()


void InterSpec::closeGammaLinesWindow()
{
  if( !m_referencePhotopeakLinesWindow )
    return;
  
  //When the "back" button is pressed on mobile phones
  if( isPhone() && m_referencePhotopeakLinesWindow->isHidden() )
    return;
  
  if( isPhone() )
  {
    m_referencePhotopeakLinesWindow->hide();
    return;
  }
  
  string xmlstate;
  if( m_referencePhotopeakLines )
  {
    if( !m_referencePhotopeakLines->currentlyShowingNuclide().empty()
         || m_referencePhotopeakLines->persistedNuclides().size() )
      m_referencePhotopeakLines->serialize( xmlstate );
    m_referencePhotopeakLines->clearAllLines();
    if( m_toolsTabs && m_toolsTabs->indexOf( m_referencePhotopeakLines ) >= 0 )
      m_toolsTabs->removeTab( m_referencePhotopeakLines );
    delete m_referencePhotopeakLines;
    m_referencePhotopeakLines = nullptr;
  }//if( m_referencePhotopeakLines )

  delete m_referencePhotopeakLinesWindow;
  m_referencePhotopeakLinesWindow = nullptr;

  if( m_toolsTabs )
  {
    m_referencePhotopeakLines = new ReferencePhotopeakDisplay( m_spectrum,
                                                   m_materialDB.get(),
                                                   m_shieldingSuggestion,
                                                   this );
    setReferenceLineColors( nullptr );
    
    m_toolsTabs->addTab( m_referencePhotopeakLines, GammaLinesTabTitle, TabLoadPolicy );
    
    if( xmlstate.size() )
      m_referencePhotopeakLines->deSerialize( xmlstate );
  }//if( m_toolsTabs )
  
  if( m_toolsTabs )
    m_currentToolsTab = m_toolsTabs->currentIndex();
}//void closeGammaLinesWindow()


void InterSpec::handleToolTabChanged( int tab )
{
  if( !m_toolsTabs )
    return;
  
  const int refTab = m_toolsTabs->indexOf(m_referencePhotopeakLines);
  const int calibtab = m_toolsTabs->indexOf(m_energyCalTool);
  const int searchTab = m_toolsTabs->indexOf(m_nuclideSearchContainer);
  
  if( m_referencePhotopeakLines && (tab == refTab) && !isMobile() )
    m_referencePhotopeakLines->setFocusToIsotopeEdit();
    
  if( m_nuclideSearch && (m_currentToolsTab==searchTab) )
    m_nuclideSearch->clearSearchEnergiesOnClient();
  
  if( m_nuclideSearch && (tab==searchTab) )
    m_nuclideSearch->loadSearchEnergiesToClient();
  
  if( tab == calibtab )
  {
    if( InterSpecUser::preferenceValue<bool>( "ShowTooltips", this ) )
      passMessage( "You can also recalibrate graphically by right-clicking and "
                   "dragging the spectrum to where you want",
                   "", WarningWidget::WarningMsgInfo );
  }//if( tab == calibtab )
  
  m_currentToolsTab = tab;
}//void InterSpec::handleToolTabChanged( int tabSwitchedTo )


SpecMeasManager *InterSpec::fileManager()
{
  return m_fileManager;
}


PeakModel *InterSpec::peakModel()
{
  return m_peakModel;
}


MaterialDB *InterSpec::materialDataBase()
{
  return m_materialDB.get();
}


Wt::WSuggestionPopup *InterSpec::shieldingSuggester()
{
  return m_shieldingSuggestion;
}


Wt::Signal<std::shared_ptr<DetectorPeakResponse> > &InterSpec::detectorChanged()
{
  return m_detectorChanged;
}


Wt::Signal<std::shared_ptr<DetectorPeakResponse> > &InterSpec::detectorModified()
{
  return m_detectorModified;
}


float InterSpec::sample_real_time_increment( const std::shared_ptr<const SpecMeas> &meas,
                                             const int sample,
                                             const std::vector<string> &detector_names )
{
  float realtime = 0.0f;
  
  if( !meas )
    return realtime;
  
  const auto &measurements = meas->sample_measurements( sample );
  
  for( const auto &m : measurements )
  {
    const auto pos = std::find( begin(detector_names), end(detector_names), m->detector_name() );
    if( pos != end(detector_names) )
      realtime = std::max( realtime, m->real_time() );
  }
  return realtime;
  
  
/*
  int nback = 0, nnonback = 0;
  double backtime = 0.0, nonbacktime = 0.0;
  for( const std::shared_ptr<const SpecUtils::Measurement> &m : measurement )
  {
    if( m->source_type() == SpecUtils::SourceType::Background )
    {
      ++nback;
      backtime += m->real_time();
    }else
    {
      ++nnonback;
      nonbacktime += m->real_time();
    }
  }//for( const std::shared_ptr<const SpecUtils::Measurement> &m : measurement )
  
  if( nnonback )
    return nonbacktime/nnonback;
  if( nback )
    return backtime/nback;
  return 0.0;
*/
}//double sample_real_time_increment()


/*
double InterSpec::liveTime( const std::set<int> &samplenums ) const
{
  double time = 0.0;
  
  if( !m_dataMeasurement )
    return 0.0;
  
  const vector<bool> det_use = detectorsToDisplay();
  const set<int> sample_numbers = validForegroundSamples();
  
  const vector<int> &detnums = m_dataMeasurement->detector_numbers();
  const vector<int>::const_iterator detnumbegin = detnums.begin();
  const vector<int>::const_iterator detnumend = detnums.end();
  
  for( int sample : samplenums )
  {
    const vector<std::shared_ptr<const SpecUtils::Measurement>> measurement
    = m_dataMeasurement->sample_measurements( sample );
    
    for( const std::shared_ptr<const SpecUtils::Measurement> &m : measurement )
    {
      const int detn = m->detector_number();
      const size_t detpos = std::find(detnumbegin,detnumend,detn) - detnumbegin;
      
      if( detpos < det_use.size() && det_use[detpos] )
        time += m->live_time();
    }//for( const std::shared_ptr<const SpecUtils::Measurement> &m : measurement )
  }//for( int sample : prev_displayed_samples )

  return time;
}//double liveTime( const std::set<int> &samplenums );
*/



void InterSpec::changeDisplayedSampleNums( const std::set<int> &samples,
                                                const SpecUtils::SpectrumType type )
{
  std::shared_ptr<SpecMeas> meas = measurment( type );
  
  if( !meas )
    return;
  
  std::shared_ptr<const SpecUtils::Measurement> prevhist = displayedHistogram(type);
  
  std::set<int> *sampleset = nullptr;
  
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      sampleset = &m_displayedSamples;
      deletePeakEdit();
    break;
      
    case SpecUtils::SpectrumType::SecondForeground:
      sampleset = &m_sectondForgroundSampleNumbers;
    break;
      
    case SpecUtils::SpectrumType::Background:
      sampleset = &m_backgroundSampleNumbers;
    break;
  }//switch( type )
  
  if( (*sampleset) == samples )
    return;
  
  (*sampleset) = samples;
  
  
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      displayForegroundData( true );
    break;
      
    case SpecUtils::SpectrumType::SecondForeground:
      displaySecondForegroundData();
    break;
      
    case SpecUtils::SpectrumType::Background:
      displayBackgroundData();
    break;
  }//switch( type )
  
  
  //Right now, we will only search for hint peaks for foreground
#if( !ANDROID && !IOS )
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      if( !!m_dataMeasurement
         && !m_dataMeasurement->automatedSearchPeaks(samples) )
        searchForHintPeaks( m_dataMeasurement, samples );
      break;
      
    case SpecUtils::SpectrumType::SecondForeground:
    case SpecUtils::SpectrumType::Background:
      break;
  }//switch( spec_type )
#endif
  
  const auto dets = detectorsToDisplay(type);
  m_displayedSpectrumChangedSignal.emit( type, meas, (*sampleset), dets );
}//void InterSpec::changeDisplayedSampleNums( const std::set<int> &samples )


void InterSpec::timeChartClicked( const int sample_number, Wt::WFlags<Wt::KeyboardModifier> modifiers )
{
  timeChartDragged( sample_number, sample_number, modifiers );
}//void timeChartClicked(...)


void InterSpec::timeChartDragged( const int sample_start_in, const int sample_end_in,
                           Wt::WFlags<Wt::KeyboardModifier> modifiers )
{
  if( !m_dataMeasurement )
    return;
 
  enum class ActionType
  {
    ChangeSamples,
    AddSamples,
    RemoveSamples
  };
  
  ActionType action = ActionType::ChangeSamples;
  if( modifiers.testFlag(KeyboardModifier::ShiftModifier) )
    action = ActionType::AddSamples;
  else if( modifiers.testFlag(KeyboardModifier::ControlModifier) )
    action = ActionType::RemoveSamples;
  
  SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
  if( modifiers.testFlag(KeyboardModifier::AltModifier) )
    type = SpecUtils::SpectrumType::Background;
  else if( modifiers.testFlag(KeyboardModifier::MetaModifier) )
    type = SpecUtils::SpectrumType::SecondForeground;
  
  const int sample_start = std::min( sample_start_in, sample_end_in );
  const int sample_end = std::max( sample_start_in, sample_end_in );
  
  const set<int> &all_samples = m_dataMeasurement->sample_numbers();
  const set<int>::const_iterator start_iter = all_samples.find( sample_start );
  set<int>::const_iterator end_iter = all_samples.find( sample_end );
  
  if( start_iter == end(all_samples) || end_iter == end(all_samples) )
  {
    passMessage( "Received invalid sample number from time chart (" + std::to_string(sample_start)
                 + ", " + std::to_string(sample_end) + ") - this shouldnt have happend"
                 " - not changing sample numbers of spectrum",
                 "", WarningWidget::WarningMsgHigh );
    return;
  }//if( invalid sample numbers )
  
  ++end_iter;
  set<int> interaction_samples;
  for( set<int>::const_iterator iter = start_iter; iter != end_iter; ++iter )
    interaction_samples.insert( *iter );
  
  assert( interaction_samples.size() );
  
  if( measurment(type) != m_dataMeasurement )
  {
    if( action != ActionType::RemoveSamples )
      setSpectrum( m_dataMeasurement, interaction_samples, type, 0 );
    return;
  }//if( the action isnt for the foreground )
  
  std::set<int> dispsamples = displayedSamples(type);
  
  switch( action )
  {
    case ActionType::ChangeSamples:
      dispsamples = interaction_samples;
      changeDisplayedSampleNums( dispsamples, type );
      break;
      
    case ActionType::AddSamples:
      dispsamples.insert( begin(interaction_samples), end(interaction_samples) );
      changeDisplayedSampleNums( dispsamples, type );
      break;
      
    case ActionType::RemoveSamples:
      for( const auto sample : interaction_samples )
        dispsamples.erase(sample);
      
      if( !dispsamples.empty() )
      {
        changeDisplayedSampleNums( dispsamples, type );
      }else
      {
        switch( type )
        {
          case SpecUtils::SpectrumType::Foreground:
            // If user erased all the samples - then lets go back to displaying the default samples.
            //  We dont want to clear the foreground file, because then we'll lose the time chart
            //  and the background/secondary spectra too.
            dispsamples = sampleNumbersForTypeFromForegroundFile(type);
            changeDisplayedSampleNums( dispsamples, type );
            break;
            
          case SpecUtils::SpectrumType::Background:
          case SpecUtils::SpectrumType::SecondForeground:
            setSpectrum( nullptr, std::set<int>{}, type, 0 );
            break;
        }//switch( type )
      }//if( !dispsamples.empty() ) / else
      
      
      break;
  }//switch( action )
}//void timeChartDragged(...)


void InterSpec::findAndSetExcludedSamples( std::set<int> definetly_keep_samples )
{
  m_excludedSamples.clear();

  if( !m_dataMeasurement || !m_dataMeasurement->passthrough() )
    return;

  const set<int> all_samples = m_dataMeasurement->sample_numbers();

  for( const int sample : all_samples )
  {
    vector< std::shared_ptr<const SpecUtils::Measurement> > measurements = m_dataMeasurement->sample_measurements( sample );

    if( definetly_keep_samples.count(sample) > 0 )
      continue;

    if( measurements.empty() )
    {
      m_excludedSamples.insert( sample );
    }else
    {
      const std::shared_ptr<const SpecUtils::Measurement> meas = measurements.front();
      //XXX - Assuming background and calibration statis is the same for all
      //      detectors
      //const bool back = (meas->source_type() == SpecUtils::SourceType::Background);
      const bool calib = (meas->source_type() == SpecUtils::SourceType::Calibration);

      if( /*back ||*/ calib )
        m_excludedSamples.insert( sample );
    }//if( measurements.empty() ) / else
  }//for( const int sample : all_samples )
}//void InterSpec::findAndSetExcludedSamples()


std::set<int> InterSpec::validForegroundSamples() const
{
  set<int> sample_nums;

  if( !m_dataMeasurement )
    return sample_nums;

  sample_nums = m_dataMeasurement->sample_numbers();
  for( const int s : m_excludedSamples )
    sample_nums.erase( s );
  
  // If we have "derived" and non-derived data, then don't let derived data be a valid foreground.
  const bool hasDerivedData = m_dataMeasurement->contains_derived_data();
  const bool hasNonDerivedData = m_dataMeasurement->contains_non_derived_data();
  
  set<int> to_rm;
  for( const int samplenum : sample_nums )
  {
    const auto meass = m_dataMeasurement->sample_measurements(samplenum);
    
    for( const std::shared_ptr<const SpecUtils::Measurement> &m : meass )
    {
      if( hasDerivedData && hasNonDerivedData && m->derived_data_properties() )
        to_rm.insert( samplenum );
      
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
        case SpecUtils::SourceType::Background:
          to_rm.insert( samplenum );
          break;
          
        case SpecUtils::SourceType::Foreground:
        case SpecUtils::SourceType::Unknown:
          break;
      }//switch( m->source_type() )
    }//for( loop over measurements of this sample number )
  }//for( const int s : sample_nums )
  
  for( const int samplenum : to_rm )
    sample_nums.erase( samplenum );

  return sample_nums;
}//std::set<int> validForegroundSamples() const


#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !IOS && !BUILD_AS_ELECTRON_APP )
void InterSpec::initOsColorThemeChangeDetect()
{
  m_osColorThemeChange.reset( new JSignal<std::string>( this, "OsColorThemeChange", true ) );
  m_osColorThemeChange->connect( boost::bind( &InterSpec::osThemeChange, this,
                                             boost::placeholders::_1 ) );
  
  LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsSetupOsColorThemeChangeJs);
  
  doJavaScript( "Wt.WT.SetupOsColorThemeChangeJs('" + id() + "')" );
}//void initOsColorThemeChangeDetect()
#endif


void InterSpec::loadDetectorResponseFunction( std::shared_ptr<SpecMeas> meas,
                                              SpecUtils::DetectorType type,
                                              const std::string serial_number,
                                              const std::string manufacturer,
                                              const std::string model,
                                              const bool tryDefaultDrf,
                                              const std::string sessionId )
{
  if( !meas )
    return;
  
  std::shared_ptr<DetectorPeakResponse> det;

  //First see if the user has opted for a detector for this serial number of
  //  detector model
  det = DrfSelect::getUserPreferredDetector( m_sql, m_user, serial_number, type, model );
  
  if( !det && (type == SpecUtils::DetectorType::Unknown) )
    return;
  
  if( !det && !tryDefaultDrf )
    return;
  
  const bool usingUserDefaultDet = !!det;
  
  if( !det )
  {
    try
    {
      det = DrfSelect::initARelEffDetector( type, manufacturer, model, this );
    }catch( std::exception & )
    {
    }
  }//if( !det )
  
  if( !det )
  {
    try
    {
      det = DrfSelect::initAGadrasDetector( type, this );
    }catch( std::exception & )
    {
    }
  }//if( !det )
  

  if( !det )
    return;
  
  WServer::instance()->post( sessionId, std::bind( [this, meas, det, usingUserDefaultDet](){
    //ToDo: could add button to remove association with DRF in database,
    //      similar to the "Start Fresh Session" button.  Skeleton code to do this
    /*
     std::unique_ptr<Wt::JSignal<> > m_clearDrfAssociation;
     m_clearDrfAssociation.reset( new JSignal<>(this, "removeDrfAssociation", false) );
     m_clearDrfAssociation->connect( this, &InterSpec::... );
     
     WStringStream js;
     js << "File contained RIID analysis results: "
     << riidAnaSummary(meas)
     << "<div onclick="
     "\"Wt.emit('" << wApp->root()->id() << "',{name:'miscSignal'}, 'showRiidAna-" << type << "');"
     //"$('.qtip.jgrowl:visible:last').remove();"
     "try{$(this.parentElement.parentElement).remove();}catch(e){}"
     "return false;\" "
     "class=\"clearsession\">"
     "<span class=\"clearsessiontxt\">Show full RIID results</span></div>";
     
     WarningWidget::displayPopupMessageUnsafe( js.str(), WarningWidget::WarningMsgShowRiid, 20000 );
     */
    
    if( meas != m_dataMeasurement )
    {
      cerr << "Foreground changed by the time DRF was loaded." << endl;
      return;
    }
    
    const bool wasModified = meas->modified();
    const bool wasModifiedSinceDecode = meas->modified_since_decode();
      
    meas->setDetector( det );
      
    if( !wasModified )
      meas->reset_modified();
      
    if( !wasModifiedSinceDecode )
      meas->reset_modified_since_decode();
    
    m_detectorChanged.emit( det );
    
    
    const char *msg = "Using the detector response function you specified to use as default for this detector.";
    if( !usingUserDefaultDet )
      msg = "Have loaded a default detector response function for this detection system.";
    
    passMessage( msg, "", WarningWidget::WarningMsgInfo );
    
    WApplication *app = WApplication::instance();
    if( app )
      app->triggerUpdate();
  }) );
  
}//void InterSpec::loadDetectorResponseFunction( WApplication *app )


void InterSpec::doFinishupSetSpectrumWork( std::shared_ptr<SpecMeas> meas,
                                  vector<boost::function<void(void)> > workers )
{
  if( !meas || workers.empty() )
    return;
  
  bool modified, modifiedSinceDecode;
  
  {//begin codeblock to access meas
    std::lock_guard<std::recursive_mutex> scoped_lock( meas->mutex() );
    modified = meas->modified();
    modifiedSinceDecode = meas->modified_since_decode();
  }//end codeblock to access meas
  
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t i = 0; i < workers.size(); ++i )
      pool.post( workers[i] );
    pool.join();
  }
  
  {//begin codeblock to access meas
    std::lock_guard<std::recursive_mutex> scoped_lock( meas->mutex() );
    if( !modified )
      meas->reset_modified();
    if( !modifiedSinceDecode )
      meas->reset_modified_since_decode();
  }//end codeblock to access meas
}//void InterSpec::doFinishupSetSpectrumWork( boost::function<void(void)> workers )


void InterSpec::setSpectrum( std::shared_ptr<SpecMeas> meas,
                             std::set<int> sample_numbers,
                             const SpecUtils::SpectrumType spec_type,
                             const Wt::WFlags<SetSpectrumOptions> options )
{
  const int spectypeindex = static_cast<int>( spec_type );
  
  vector< boost::function<void(void)> > furtherworkers;
  
  const bool wasModified = (meas ? meas->modified() : false);
  const bool wasModifiedSinceDecode = (meas ? meas->modified_since_decode() : false);
  
  if( m_useInfoWindow && meas )
  {
    // If we are loading a state from the "Welcome To InterSpec" screen, we dont want to delete
    //  m_useInfoWindow because we will still use it, so instead we'll try deleting the window on
    //  the next go around of the event loop.
    auto doDelete = wApp->bind( std::bind([this](){
      WApplication *app = wApp;
      if( !app )
        return;
      deleteWelcomeCountDialog();
      app->triggerUpdate();
    }) );
      
    WServer::instance()->post( wApp->sessionId(), doDelete );
  }//if( meas )
  
  std::shared_ptr<SpecMeas> previous = measurment(spec_type);
  const set<int> prevsamples = displayedSamples(spec_type);
  const vector<string> prevdets = detectorsToDisplay(spec_type);
  const bool sameSpec = (meas==previous);
  std::shared_ptr<const SpecUtils::Measurement> prev_display = m_spectrum->histUsedForXAxis();
  

  if( (spec_type == SpecUtils::SpectrumType::Foreground) && previous && (previous != meas) )
  {
    closeShieldingSourceFitWindow();
    
#if( USE_DB_TO_STORE_SPECTRA )
    //if( m_user->preferenceValue<bool>( "AutoSaveSpectraToDb" ) )
    //{
    //  //We also need to do this in the InterSpec destructor as well.
    //  //   Also maybe change size limitations to only apply to auto saving
    //  if( m_currentStateID >= 0 )
    //  {
    //    //Save to (HEAD) of current state
    //  }else
    //  {
    //    //Create a state
    //    //Handle case where file is to large to be saved
    //  }
    //}
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  }//if( (spec_type == SpecUtils::SpectrumType::Foreground) && !!previous && (previous != meas) )
  
  if( !!meas && isMobile() && !toolTabsVisible()
      /* && options.testFlag(SetSpectrumOptions::CheckToPreservePreviousEnergyCal) */
      && m_referencePhotopeakLines  && (spec_type == SpecUtils::SpectrumType::Foreground) )
  {
    m_referencePhotopeakLines->clearAllLines();
  }
  
  string msg;
  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
#if( USE_DB_TO_STORE_SPECTRA )
      m_currentStateID = -1;
      updateSaveWorkspaceMenu();
#endif
      if( !sameSpec )
        deletePeakEdit();
    break;
    
    case SpecUtils::SpectrumType::SecondForeground:
    case SpecUtils::SpectrumType::Background:
    break;
  }//switch( spec_type )

  if( meas
      /*&& sample_numbers.empty() */
      && !sameSpec
      && (spec_type==meas->displayType())
      && meas->displayedSampleNumbers().size() )
  {
    sample_numbers = meas->displayedSampleNumbers();
  }

#if( USE_SAVEAS_FROM_MENU )
  if( m_downloadMenu && m_downloadMenus[spectypeindex] )
  {
    m_downloadMenus[spectypeindex]->setDisabled( !meas );
    WMenuItem *item = m_downloadMenus[spectypeindex]->parentItem();
    if( item )
      m_downloadMenu->setItemHidden( item, !meas );
      
    bool allhidden=true;
    for( SpecUtils::SpectrumType i = SpecUtils::SpectrumType(0);
         i <= SpecUtils::SpectrumType::Background;
         i = SpecUtils::SpectrumType(static_cast<int>(i)+1) )
    {
      if (!m_downloadMenu->isItemHidden( m_downloadMenus[static_cast<int>(i)]->parentItem()))
      {
        allhidden=false;
        break;
      }//if (!m_downloadMenus[i]->isHidden())
    }//    for( SpecUtils::SpectrumType i = SpecUtils::SpectrumType(0); i<=SpecUtils::SpectrumType::Background; i = SpecUtils::SpectrumType(i+1) )
      
    m_fileMenuPopup->setItemHidden(m_downloadMenu->parentItem(),allhidden);
  }//if( m_downloadMenu && m_downloadMenus[spec_type] )
#endif
  
  
  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      m_detectorChangedConnection.disconnect();
      m_detectorModifiedConnection.disconnect();
      m_displayedSpectrumChanged.disconnect();
      
      if( meas )
      {
        //Lets keep using the same detector if we are loading a new spectrum
        //  with the same number of bins, but doesnt have a detector of its own
        const bool sameNBins = ( m_dataMeasurement && meas
               && (m_dataMeasurement->num_gamma_channels()==meas->num_gamma_channels()) );

        std::shared_ptr<DetectorPeakResponse> old_det;
        if( m_dataMeasurement )
          old_det = m_dataMeasurement->detector();

        if( meas && (!meas->detector() || !meas->detector()->isValid() )
            && old_det && old_det->isValid()
            && sameNBins && meas->num_gamma_channels() )
        {
          meas->setDetector( old_det );
        }

        if( meas->detector() != old_det )
        {
          auto drf = meas->detector();
          m_detectorChanged.emit( meas->detector() );
          
          if( drf )
          {
            auto drfcopy = std::make_shared<DetectorPeakResponse>( *drf );
            boost::function<void(void)> worker = boost::bind( &DrfSelect::updateLastUsedTimeOrAddToDb, drfcopy, m_user.id(), m_sql );
            WServer::instance()->ioService().boost::asio::io_service::post( worker );
          }
        }//if( meas->detector() != old_det )
      
        if( !meas->detector() /* && (detType != SpecUtils::DetectorType::Unknown) */ )
        {
          const bool doLoadDefault = InterSpecUser::preferenceValue<bool>( "LoadDefaultDrf", this );
          
          const string &manufacturer = meas->manufacturer();
          const string &model = meas->instrument_model();
          const string &serial_num = meas->instrument_id();
          SpecUtils::DetectorType type = meas->detector_type();
          if( type == SpecUtils::DetectorType::Unknown )
            type = SpecMeas::guessDetectorTypeFromFileName( meas->filename() );
          
          boost::function<void()> worker
                  = wApp->bind( boost::bind( &InterSpec::loadDetectorResponseFunction,
                          this, meas, type, serial_num, manufacturer, model, doLoadDefault, wApp->sessionId() ) );
          furtherworkers.push_back( worker );
        }//if( we could try to load a detector type )
        
        m_detectorChangedConnection = m_detectorChanged.connect( boost::bind( &SpecMeas::detectorChangedCallback, meas.get(), boost::placeholders::_1 ) );
        m_detectorModifiedConnection = m_detectorModified.connect( boost::bind( &SpecMeas::detectorChangedCallback, meas.get(), boost::placeholders::_1 ) );
        m_displayedSpectrumChanged = m_displayedSpectrumChangedSignal.connect( boost::bind( &SpecMeas::displayedSpectrumChangedCallback, meas.get(), boost::placeholders::_1,
            boost::placeholders::_2, boost::placeholders::_3, boost::placeholders::_4 ) );
      }//if( meas )

      m_dataMeasurement = meas;
      
      findAndSetExcludedSamples( sample_numbers );

      if( !sameSpec && m_shieldingSourceFit )
        m_shieldingSourceFit->newForegroundSet();
    break;

    case SpecUtils::SpectrumType::SecondForeground:
      m_secondDataMeasurement = meas;
      m_sectondForgroundSampleNumbers = sample_numbers;
      if( meas && m_sectondForgroundSampleNumbers.empty() )
        m_sectondForgroundSampleNumbers = meas->sample_numbers();
    break;

    case SpecUtils::SpectrumType::Background:
      m_backgroundMeasurement = meas;
      m_backgroundSampleNumbers = sample_numbers;
      if( meas && m_backgroundSampleNumbers.empty() )
        m_backgroundSampleNumbers = meas->sample_numbers();
    break;
  };//switch( spec_type )

  if( msg.size() )
    passMessage( msg, "", 0 );

  
  //If loading a new foreground that has a different number of channels than
  //  the background/secondary, and differnet number of bins than previous
  //  foreground, get rid of the background/secondary since the user probably
  //  isnt interested in files from a completely different detector anymore.
  if( spec_type == SpecUtils::SpectrumType::Foreground && m_dataMeasurement && !sameSpec )
  {
    //Assume we will use all the detectors (just to determine binning)
    const vector<string> &detectors = m_dataMeasurement->detector_names();
    
    shared_ptr<const vector<float>> binning;
    try
    {
      auto energy_cal = m_dataMeasurement->suggested_sum_energy_calibration( sample_numbers, detectors );
      if( energy_cal )
        binning = energy_cal->channel_energies();
    }catch( std::exception & )
    {
      
    }
    
    shared_ptr<const vector<float>> prev_binning = prev_display ? prev_display->channel_energies() : nullptr;
    
    const bool diff_fore_nchan = ((!prev_binning || !binning) || (prev_binning->size() != binning->size()));
    
    const size_t num_foreground_channels = binning ? binning->size() : 0;
    size_t num_sec_channel = 0, num_back_channel = 0;
    if( m_secondDataMeasurement )
      num_sec_channel = m_secondDataMeasurement->num_gamma_channels();
    if( m_backgroundMeasurement )
      num_back_channel = m_backgroundMeasurement->num_gamma_channels();
    
    if( diff_fore_nchan && num_sec_channel && (num_sec_channel != num_foreground_channels) )
    {
#if( USE_SAVEAS_FROM_MENU )
      m_downloadMenus[static_cast<int>(SpecUtils::SpectrumType::SecondForeground)]->setDisabled( true );
      WMenuItem *item = m_downloadMenus[static_cast<int>(SpecUtils::SpectrumType::SecondForeground)]->parentItem();
      if( item )
        m_downloadMenu->setItemHidden( item, true );
#endif
      
      m_secondDataMeasurement = nullptr;
      m_spectrum->setSecondData( nullptr );
      
      m_displayedSpectrumChangedSignal.emit( SpecUtils::SpectrumType::SecondForeground,
                                             nullptr, {}, {} );
    }//if( num_sec_channel )
    
    if( diff_fore_nchan && num_back_channel && num_foreground_channels && (num_back_channel != num_foreground_channels) )
    {
#if( USE_SAVEAS_FROM_MENU )
      m_downloadMenus[static_cast<int>(SpecUtils::SpectrumType::Background)]->setDisabled( true );
      WMenuItem *item = m_downloadMenus[static_cast<int>(SpecUtils::SpectrumType::Background)]->parentItem();
      if( item )
        m_downloadMenu->setItemHidden( item, true );
#endif
      
      m_backgroundMeasurement = nullptr;
      m_spectrum->setBackground( nullptr );
      m_displayedSpectrumChangedSignal.emit( SpecUtils::SpectrumType::Background, nullptr, {}, {} );
    }//if( nSecondBins )
  }//if( spec_type == SpecUtils::SpectrumType::Foreground )
  
  
  
  //Fall throughs intentional
  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      updateGuiForPrimarySpecChange( sample_numbers );
#if( USE_DB_TO_STORE_SPECTRA )
      m_saveStateAs->setDisabled( !m_dataMeasurement );
#endif
      displayForegroundData( false );
      displayTimeSeriesData();
      
    case SpecUtils::SpectrumType::SecondForeground:
      displaySecondForegroundData();
      
    case SpecUtils::SpectrumType::Background:
      displayBackgroundData();
  };//switch( spec_type )
  
  
  if( !sameSpec )
  {
    const bool enableScaler = (m_secondDataMeasurement || m_backgroundMeasurement);
    const int windex_0 = m_spectrum->yAxisScalersIsVisible() ? 0 : 1;
    const int windex_1 = m_spectrum->yAxisScalersIsVisible() ? 1 : 0;
    m_showYAxisScalerItems[windex_0]->hide();
    m_showYAxisScalerItems[windex_0]->enable();
    m_showYAxisScalerItems[windex_1]->show();
    m_showYAxisScalerItems[windex_1]->setDisabled( !enableScaler );
  }//if( !sameSpec )
  
  //Making fcn call take current data as a argument so that if this way a
  //  recalibration happens (which will change m_spectrum->data()), then the
  //  peak fit routine will get the correct data to use.
  std::function<void(std::shared_ptr<const SpecUtils::Measurement>)> propigate_peaks_fcns;
  
  const bool askToPropigatePeaks
          = InterSpecUser::preferenceValue<bool>( "AskPropagatePeaks", this );
  if( askToPropigatePeaks
     && options.testFlag(InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal)
     && !sameSpec && meas && m_dataMeasurement && previous && m_spectrum->data()
     && spec_type==SpecUtils::SpectrumType::Foreground
     && previous->instrument_id()==meas->instrument_id()
     && previous->num_gamma_channels()==meas->num_gamma_channels() )
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> prevpeak, currpeaks;
    
    //Call const version of peaks so a new deque wont be created if it doesnt exist.
    prevpeak = std::const_pointer_cast<const SpecMeas>(previous)->peaks(prevsamples);
    currpeaks = std::const_pointer_cast<const SpecMeas>(m_dataMeasurement)->peaks(m_displayedSamples);
    
    if( prevpeak && !prevpeak->empty() && (!currpeaks || currpeaks->empty()) )
    {
      std::vector<PeakDef> input_peaks;
      for( const auto &p : *prevpeak )
      {
        if( p )  //Shouldnt be necassary, but JIC
          input_peaks.push_back( *p );
      }
      
      vector<PeakDef> original_peaks;
      
      const std::string sessionid = wApp->sessionId();
      propigate_peaks_fcns = [=]( std::shared_ptr<const SpecUtils::Measurement> data ){
        PeakSearchGuiUtils::fit_template_peaks( this, data, input_peaks, original_peaks,
                       PeakSearchGuiUtils::PeakTemplateFitSrc::PreviousSpectrum,
                       sessionid );
      };
    }//if( prev spec had peaks and new one doesnt )
  }//if( should propogate peaks )
  
  
  
  deleteEnergyCalPreserveWindow();
  
  if( options.testFlag(SetSpectrumOptions::CheckToPreservePreviousEnergyCal)
      && !sameSpec && m_energyCalTool && !!meas && !!m_dataMeasurement )
  {
    switch( spec_type )
    {
      case SpecUtils::SpectrumType::Foreground:
        if( EnergyCalPreserveWindow::candidate(meas,previous) )
          m_preserveCalibWindow = new EnergyCalPreserveWindow( meas, spec_type,
                                         previous, spec_type, m_energyCalTool );
      break;
    
      case SpecUtils::SpectrumType::SecondForeground:
      case SpecUtils::SpectrumType::Background:
        if( EnergyCalPreserveWindow::candidate(meas,m_dataMeasurement) )
          m_preserveCalibWindow = new EnergyCalPreserveWindow( meas, spec_type,
                                                m_dataMeasurement, SpecUtils::SpectrumType::Foreground,
                                                m_energyCalTool );
      break;
    };//switch( spec_type )
  
    if( m_preserveCalibWindow )
    {
      if( propigate_peaks_fcns )
      {
        m_preserveCalibWindow->finished().connect( std::bind( [=](){
          deleteEnergyCalPreserveWindow();
          std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
          WServer::instance()->ioService().boost::asio::io_service::post( std::bind([=](){ propigate_peaks_fcns(data); }) );
        } ) );
        
        propigate_peaks_fcns = nullptr;
      }else
      {
        m_preserveCalibWindow->finished().connect( this, &InterSpec::deleteEnergyCalPreserveWindow );
      }//if( propigate_peaks_fcns ) / else
      
    }
  }//if( !sameSpec && m_energyCalTool && !!meas )
  
  if( propigate_peaks_fcns )
  {
    std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
    WServer::instance()->ioService().boost::asio::io_service::post( std::bind([=](){ propigate_peaks_fcns(data); }) );
    propigate_peaks_fcns = nullptr;
  }
  
  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      if( !!m_dataMeasurement && ((!m_dataMeasurement)!=(!previous)) )
        doJavaScript( "$('.Wt-domRoot').data('HasForeground',1);" );
      else if( !m_dataMeasurement )
        doJavaScript( "$('.Wt-domRoot').data('HasForeground',0);" );
    break;
    
    case SpecUtils::SpectrumType::SecondForeground: case SpecUtils::SpectrumType::Background:
      if( !!meas && !m_dataMeasurement )
      {
        passMessage( "You must load a foreground spectrum before viewing a"
                     " background or second foreground spectrum.",
                    "", WarningWidget::WarningMsgHigh );
      }
    break;
  }//switch( spec_type )
  
  
  // Update the energy calibration tool, as there is new data.
  const auto shownDets = detectorsToDisplay(spec_type);
  m_displayedSpectrumChangedSignal.emit( spec_type, meas, sample_numbers, shownDets );
  
  if( meas )
  {
    if( !wasModified )
      meas->reset_modified();
    if( !wasModifiedSinceDecode )
      meas->reset_modified_since_decode();
  }//if( meas )
  
#if( USE_GOOGLE_MAP )
//  m_secondDataMeasurement && m_dataMeasurement m_backgroundMeasurement
  const bool hasGps = (m_dataMeasurement && m_dataMeasurement->has_gps_info());
  if( m_mapMenuItem )
    m_mapMenuItem->setDisabled( !hasGps );
#endif
  
#if( USE_SEARCH_MODE_3D_CHART )
  const bool isSearchData = (m_dataMeasurement && m_dataMeasurement->passthrough());
  if( m_searchMode3DChart )
    m_searchMode3DChart->setDisabled( !isSearchData );
#endif

  if( m_showRiidResults )
  {
    const bool showRiid = m_dataMeasurement && m_dataMeasurement->detectors_analysis();
    m_showRiidResults->setDisabled( !showRiid );
  }
  
  //Right now, we will only search for hint peaks for foreground
#if( !ANDROID && !IOS )
  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
    {
      if( m_dataMeasurement )
      {
        auto peaks = m_dataMeasurement->automatedSearchPeaks(sample_numbers);
        if( !peaks )
          searchForHintPeaks( m_dataMeasurement, sample_numbers );
      }
      break;
    }
      
    case SpecUtils::SpectrumType::SecondForeground:
    case SpecUtils::SpectrumType::Background:
      break;
  }//switch( spec_type )
#endif
  
  
  //Lets see if there are any parse warnings that we should give to the user.
  if( meas && !sameSpec )
  {
    Wt::WApplication *app = wApp;
    
    auto checkForWarnings = [sample_numbers,app,this,meas](){
      set<string> givenwarnings;
      for( const auto &msg : meas->parse_warnings() )
        givenwarnings.insert( msg );
      
      // @TODO Check if sample_numbers could be empty if the user wants all samples displayd
      for( const auto &m : meas->measurements() )
      {
        if( !m || !sample_numbers.count(m->sample_number()) )
          continue;
        
        for( const auto &msg : m->parse_warnings() )
          givenwarnings.insert( msg );
      }//for( const auto &m : meas->measurements() )
      
      if( !givenwarnings.empty() )
      {
        WApplication::UpdateLock lock(app);
        if( lock )
        {
          for( const auto &msg : givenwarnings )
            passMessage( msg, "", WarningWidget::WarningMsgMedium );
          app->triggerUpdate();
        }//
      }//if( !givenwarnings.empty() )
    };//checkForWarnings lamda
    
    furtherworkers.push_back( checkForWarnings );
  }//if( meas && !sameSpec )
  
  // Check if there are RIID analysis results in the file, and if so let the user know.
  if( !sameSpec && meas && meas->detectors_analysis()
     && options.testFlag(SetSpectrumOptions::CheckForRiidResults) )
  {
    auto ana = meas->detectors_analysis();
    
    // We'll show the analysis results for the foreground, no matter what if its not-empty (e.g., it
    //  might only have algorithm version information).  But for background and second foreground,
    //  we'll only show if there is an identified nuclide - perhaps we shouldnt ever show for the
    //  background/secondary.
    bool worthShowing = true;
    switch( spec_type )
    {
      case SpecUtils::SpectrumType::Foreground:
        worthShowing = !ana->is_empty();
        break;
        
      case SpecUtils::SpectrumType::SecondForeground:
      case SpecUtils::SpectrumType::Background:
        worthShowing = false;
        for( const auto &r : ana->results_ )
          worthShowing = (worthShowing || !r.nuclide_.empty());
        break;
    }//switch( spec_type )
    
    // Also, only show popup if we are only using this file for this display type.
    //  Note though that if someone loads a file as a background, then changes it to a foreground,
    //  they wont get the RIID notification popup
    const int nusedfor = static_cast<int>( meas == m_dataMeasurement )
                         + static_cast<int>( meas == m_backgroundMeasurement )
                         + static_cast<int>( meas == m_secondDataMeasurement );
    
    // Only show notification when we arent already showing the file
    if( worthShowing && (nusedfor == 1) )
    {
      const std::string type = SpecUtils::descriptionText(spec_type);
      WStringStream js;
      js << "File contained RIID analysis results: "
      << riidAnaSummary(meas)
      << "<div onclick="
           "\"Wt.emit('" << wApp->root()->id() << "',{name:'miscSignal'}, 'showRiidAna-" << type << "');"
           //"$('.qtip.jgrowl:visible:last').remove();"
           "try{$(this.parentElement.parentElement).remove();}catch(e){}"
           "return false;\" "
           "class=\"clearsession\">"
         "<span class=\"clearsessiontxt\">Show full RIID results</span></div>";
      
      WarningWidget::displayPopupMessageUnsafe( js.str(), WarningWidget::WarningMsgShowRiid, 20000 );
    }//if( nusedfor == 1 )
  }//if( meas && !sameSpec )
  
  if( meas && furtherworkers.size() )
  {
    boost::function<void(void)> worker = boost::bind(
                                  &InterSpec::doFinishupSetSpectrumWork,
                                  this, meas, furtherworkers );
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
  }//if( meas && furtherworkers.size() )
  
  if( m_mobileBackButton && m_mobileForwardButton )
  {
    if( !toolTabsVisible() && meas && spec_type==SpecUtils::SpectrumType::Foreground
        && !meas->passthrough() && (meas->sample_numbers().size()>1) )
    {
      m_mobileBackButton->setHidden(false);
      m_mobileForwardButton->setHidden(false);
    }else
    {
      m_mobileBackButton->setHidden(true);
      m_mobileForwardButton->setHidden(true);
    }
  }//if( m_mobileBackButton && m_mobileForwardButton )
  
  
  //Check for warnings from parsing the file, and display to user - not
  //  implemented well yet.
  if( m_dataMeasurement && (spec_type == SpecUtils::SpectrumType::Foreground) )
  {
    shared_ptr<const SpecUtils::EnergyCalibration> energy_cal;
    try
    {
      energy_cal = m_dataMeasurement->suggested_sum_energy_calibration( sample_numbers,
                                                             m_dataMeasurement->detector_names() );
    }catch( std::exception & )
    {
    }
    
    if( !energy_cal )
      passMessage( "Warning, no energy calibration for the selected samples",
                   "", WarningWidget::WarningMsgMedium );
  }//if( spec_type == SpecUtils::SpectrumType::Foreground && m_dataMeasurement && !sameSpec )
  

  //Display a notice to the user about how they can select different portions of
  //  passthrough/search-mode data
  /*
  if( spec_type==SpecUtils::SpectrumType::Foreground && !!m_dataMeasurement
      && m_dataMeasurement->passthrough() )
  {
    const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", this );
  
    if( showToolTips )
    {
      const char *tip = "Clicking and dragging on the time-series (bottom)"
      " chart, will change the time range the energy spectrum"
      " is summed over.  You can also shift-click or"
      " shift-drag to add additional time spans."
      " Shift-clicking or shift dragging entirely within a"
      " currently used (highlighted) time span will remove"
      " that portion of the time span. Holding the 'alt' key"
      " while doing the above will perform the same actions,"
      " but for the background if it is the same spectrum"
      " file as the foreground.";
      passMessage( tip, "", WarningWidget::WarningMsgInfo );
    }//if( showToolTips )
  }//if( passthrough foreground )
   */
}//void setSpectrum(...)


void InterSpec::reloadCurrentSpectrum( SpecUtils::SpectrumType spec_type )
{
  std::shared_ptr<SpecMeas> meas;
  std::set<int> sample_numbers;

  switch( spec_type )
  {
    case SpecUtils::SpectrumType::Foreground:
      meas = m_dataMeasurement;
      sample_numbers = m_displayedSamples;
    break;

    case SpecUtils::SpectrumType::SecondForeground:
      meas = m_secondDataMeasurement;
      sample_numbers = m_sectondForgroundSampleNumbers;
    break;

    case SpecUtils::SpectrumType::Background:
      meas = m_backgroundMeasurement;
      sample_numbers = m_backgroundSampleNumbers;
    break;
  }//switch( spec_type )

  setSpectrum( meas, sample_numbers, spec_type, 0 );
}//void reloadCurrentSpectrum( SpecUtils::SpectrumType spec_type )



void InterSpec::finishLoadUserFilesystemOpenedFile(
                                std::shared_ptr<SpecMeas> meas,
                                std::shared_ptr<SpectraFileHeader> header,
                                const SpecUtils::SpectrumType type )
{
  SpectraFileModel *fileModel = m_fileManager->model();
  
  try
  {
    const int row = fileModel->addRow( header );
    m_fileManager->displayFile( row, meas, type, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
  }catch( std::exception & )
  {
    passMessage( "There was an error loading "
                 + (!!meas ? meas->filename() : string("spectrum file")),
                "", WarningWidget::WarningMsgHigh );
  }//try / catch
  
}//finishLoadUserFilesystemOpenedFile(...)


void InterSpec::promptUserHowToOpenFile( std::shared_ptr<SpecMeas> meas,
                                             std::shared_ptr<SpectraFileHeader> header )
{
  const char *msg =
  "This file looks like it's from the same detector as the current foreground."
  "<p>How would you like to open this spectrum file?</p>";
  
  string filename = header->displayName().toUTF8();
  if( filename.size() > 48 )
  {
    SpecUtils::utf8_limit_str_size( filename, 48 );
    filename += "...";
  }
  
  SimpleDialog *dialog = new SimpleDialog( WString::fromUTF8(filename), WString::fromUTF8(msg) );
  WPushButton *button = dialog->addButton( WString::fromUTF8("Foreground") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                          meas, header, SpecUtils::SpectrumType::Foreground ) );
  button->setFocus( true );
  
  button = dialog->addButton( WString::fromUTF8("Background") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                          meas, header, SpecUtils::SpectrumType::Background ) );
  
  button = dialog->addButton( WString::fromUTF8("Secondary") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                         meas, header, SpecUtils::SpectrumType::SecondForeground ) );
}//void promptUserHowToOpenFile(...)



bool InterSpec::userOpenFileFromFilesystem( const std::string path, std::string displayFileName  )
{
  try
  {
    if( displayFileName.empty() )
      displayFileName = SpecUtils::filename(path);
    
    if( !m_fileManager )  //shouldnt ever happen.
      throw runtime_error( "Internal logic error, no valid m_fileManager" );
    
    if( !SpecUtils::is_file(path) )
      throw runtime_error( "Could not access file '" + path + "'" );
  
    auto header = std::make_shared<SpectraFileHeader>( m_user, true, this );
    auto meas = header->setFile( displayFileName, path, SpecUtils::ParserType::Auto );
  
    if( !meas )
      throw runtime_error( "Failed to decode file" );
    
    bool couldBeBackground = true;
    if( !m_dataMeasurement || !meas
        || meas->uuid() == m_dataMeasurement->uuid() )
    {
      couldBeBackground = false;
    }else
    {
      couldBeBackground &= (m_dataMeasurement->instrument_id() == meas->instrument_id());
      couldBeBackground &= (m_dataMeasurement->num_gamma_channels() == meas->num_gamma_channels());
    }//if( !m_dataMeasurement || !meas )

    //Should we check if this meas has the same UUID as the second of background?
    
    if( couldBeBackground )
    {
      promptUserHowToOpenFile( meas, header );
      return true;
    }else
    {
//      return m_fileManager->loadFromFileSystem( path, SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::Auto );
      cout << "Will load file " << path << " requested to be loaded at "
      << WDateTime::currentDateTime().toString(DATE_TIME_FORMAT_STR)
      << endl;
      
      SpectraFileModel *fileModel = m_fileManager->model();
      const int row = fileModel->addRow( header );
      m_fileManager->displayFile( row, meas, SpecUtils::SpectrumType::Foreground, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
      return true;
    }
  }catch( std::exception &e )
  {
    cerr << "Caught exception '" << e.what() << "' when trying to load '"
         << path << "'" << endl;
    SpecMeasManager::displayInvalidFileMsg( displayFileName, e.what() );
    return false;
  }//try / catch
  
  return true;
}//bool userOpenFileFromFilesystem( const std::string filepath )


void InterSpec::handleAppUrl( std::string url )
{
  //Get rid of (optional) leading "interspec://", so URL will look like: 'drf/specify?v=1&n=MyName&"My other Par"'
  const string scheme = "interspec://";
  if( SpecUtils::istarts_with(url, scheme) )
    url = url.substr(scheme.size());
  
  // For the moment, we will require there to be a query string, even if its empty.
  string::size_type q_start_pos = url.find( '?' );
  string::size_type q_end_pos = q_start_pos + 1;
  if( q_start_pos == string::npos )
  {
    // Since '?' isnt a QR alphanumeric code also allow the URL encoded value of it, "%3F"
    q_start_pos = url.find( "%3F" );
    if( q_start_pos == string::npos )
      q_start_pos = url.find( "%3f" );
    
    if( q_start_pos != string::npos )
      q_end_pos = q_start_pos + 3;
  }//
    
  if( q_start_pos == string::npos )
    throw runtime_error( "App URL did not contain the host/path component that specifies the intent"
                         " of the URL (e.g., what tool to use, or what info is contained in the URL)." );
  
  // the q_pos+1 should be safe, even if q_pos is last character in string.
  if( url.find(q_end_pos, 'q') != string::npos )
    throw runtime_error( "App URL contained more than one '?' character, which isnt allowed" );
    
  const string host_path = url.substr( 0, q_start_pos );
  const string::size_type host_end = host_path.find( '/' );
  const string host = (host_end == string::npos) ? host_path : host_path.substr(0,host_end);
  const string path = (host_end == string::npos) ? string("") : host_path.substr(host_end+1);
  const string query_str = url.substr( q_end_pos + 1 );
  
  cout << "host='" << host << "' and path='" << path << "' and query_str='" << query_str << "'" << endl;
  
  if( SpecUtils::iequals_ascii(host,"drf") )
  {
    if( !SpecUtils::iequals_ascii(path,"specify") )
      throw runtime_error( "App 'drf' URL with path '" + path + "' not supported." );
    
    DrfSelect::handle_app_url_drf( query_str );
  }else
  {
    throw runtime_error( "App URL with purpose (host-component) '" + host + "' not supported." );
  }
}//void handleAppUrl( std::string url )



void InterSpec::detectorsToDisplayChanged()
{
  displayBackgroundData();
  displaySecondForegroundData();
  displayForegroundData( true );
  displayTimeSeriesData();
 
  // \TODO: The foreground may be pass-through, but the user could have just de-selected a detector
  //        so that none of the currently selected sample numbers have gamma/neutron data, and in
  //        this case we should fix things up.  Or it could be the case that the remaining detectors
  //        are not time-series (e.g., only have one or two spectra), so we then no longer would
  //        need/want the time chart (and vice versa when adding a detector).
  
  // This function is only called when a checkbox in the "Detectors" sub-menu is changed, so for the
  //  moment, we will only emit that things changed for the foreground.  In the future we should get
  //  rid of this sub-menu and handle things correctly.
  const auto type = SpecUtils::SpectrumType::Foreground;
  const auto meas = measurment(type);
  const auto &samples = displayedSamples(type);
  const auto detectors = detectorsToDisplay(type);
  m_displayedSpectrumChangedSignal.emit(type,meas,samples,detectors);
}//void detectorsToDisplayChanged()


void InterSpec::updateGuiForPrimarySpecChange( std::set<int> display_sample_nums )
{
  m_displayedSamples = display_sample_nums;

  if( m_detectorToShowMenu )
  {
    const vector<WMenuItem *> items = m_detectorToShowMenu->items();
    for( WMenuItem *item : items )
    {
      if( !item->hasStyleClass("PhoneMenuBack") )
      {
//        m_detectorToShowMenu->removeItem( item );  //This seems unecassary, leave commented as test
        delete item;
      }
    }
    
#if( USING_ELECTRON_NATIVE_MENU )
#if( !defined(WIN32) )
    //20190125: hmm, looks like detectors menu is behaving okay - I guess I fixed it somewhere else?
    //#warning "Need to do something to get rid of previous detectors from Electrons menu"
#endif
//  https://github.com/electron/electron/issues/527
    m_detectorToShowMenu->clearElectronMenu();
#endif
  }//if( m_detectorToShowMenu )
  
  m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::Foreground );
  m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::Background );
  m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::SecondForeground );
  
  if( m_displayedSamples.empty() )
    m_displayedSamples = validForegroundSamples();

  const vector<int> det_nums = m_dataMeasurement
                               ? m_dataMeasurement->detector_numbers()
                               : vector<int>();
  const vector<string> &det_names = m_dataMeasurement
                                    ? m_dataMeasurement->detector_names()
                                    : vector<string>();
  const vector<string> &neut_det_names = m_dataMeasurement
                                 ? m_dataMeasurement->neutron_detector_names()
                                 : vector<string>();

  for( size_t index = 0; index < det_nums.size(); ++index )
  {
    const int detnum = det_nums[index];
    
    string name;
    try
    {
      name = det_names.at(index);
    }catch(...)
    {
      cerr << "InterSpec::updateGuiForPrimarySpecChange(...)\n\tSerious logic error = please fix" << endl;
      continue;
    }

    char title[512];
    if( name.size() )
      snprintf( title, sizeof(title), "Det. %i (%s)", detnum, name.c_str() );
    else
      snprintf( title, sizeof(title), "Det. %i", detnum );
    
    vector<string>::const_iterator neut_name_pos;
    neut_name_pos = std::find( neut_det_names.begin(), neut_det_names.end(), name );
    const int isneut = (neut_name_pos != neut_det_names.end());
    //if( isneut )
      //snprintf( title, sizeof(title), "Det. %i (%s)", i, name.c_str() );

    if( m_detectorToShowMenu )
    {
#if( WT_VERSION>=0x3030300 )
#if( USE_OSX_NATIVE_MENU  || USING_ELECTRON_NATIVE_MENU )
      WCheckBox *cb = new WCheckBox( title );
      cb->setChecked( true );
      PopupDivMenuItem *item = m_detectorToShowMenu->addWidget( cb, false );
#else
      PopupDivMenuItem *item = m_detectorToShowMenu->addMenuItem(title,"",false);
      item->setCheckable( true );
      item->setChecked( true );
      Wt::WCheckBox *cb = item->checkBox();
      if( !cb )
        throw runtime_error( "Serious error creating checkbox in menu item" );
#endif
      
#else
      WCheckBox *cb = new WCheckBox( title );
      cb->setChecked(true);
      
      //NOTE: this is necessary to prevent problems in menu state
      //cb->checked().connect( boost::bind(&WCheckBox::setChecked, cb, true) );
      //cb->unChecked().connect( boost::bind(&WCheckBox::setChecked, cb, false) );
      
      PopupDivMenuItem *item = m_detectorToShowMenu->addWidget( cb, false );
#endif
      
      if( isneut )
        item->addStyleClass( "NeutDetCbItem" );
      
      cb->checked().connect( this, &InterSpec::detectorsToDisplayChanged );
      cb->unChecked().connect( this, &InterSpec::detectorsToDisplayChanged );
      
      //item->triggered().connect( boost::bind( &InterSpec::detectorsToDisplayChanged, this, item ) );
      
      if( det_nums.size() == 1 )
        item->disable();
    }//if( m_detectorToShowMenu )
  }//for( gamma detector )
  
  if( m_detectorToShowMenu && m_detectorToShowMenu->parentItem() )
    m_detectorToShowMenu->parentItem()->setDisabled( det_nums.empty() );
}//bool updateGuiForPrimarySpecChange( const std::string &filename )


size_t InterSpec::addHighlightedEnergyRange( const float lowerEnergy,
                                            const float upperEnergy,
                                            const WColor &color )
{
  return m_spectrum->addDecorativeHighlightRegion( lowerEnergy, upperEnergy, color );
}//void setHighlightedEnergyRange( double lowerEnergy, double upperEnergy )


bool InterSpec::removeHighlightedEnergyRange( const size_t regionid )
{
  return m_spectrum->removeDecorativeHighlightRegion( regionid );
}//bool removeHighlightedEnergyRange( const size_t regionid );



void InterSpec::setDisplayedEnergyRange( float lowerEnergy, float upperEnergy )
{
  if( upperEnergy < lowerEnergy )
    std::swap( upperEnergy, lowerEnergy );
  
  m_spectrum->setXAxisRange( lowerEnergy, upperEnergy );
}//void setDisplayedEnergyRange()


bool InterSpec::setYAxisRange( float lower_counts, float upper_counts )
{
  bool success = true;
  if( upper_counts < lower_counts )
    std::swap( lower_counts, upper_counts );
  
  if( (lower_counts <= 1.0E-6) && (upper_counts <= 1.0E-6) )
    return false;
  
  if( (lower_counts < 9.9E-7) && m_spectrum->yAxisIsLog() )
  {
    success = false;
    lower_counts = 1.0E-6;
    if( upper_counts <= lower_counts )
      upper_counts = 10*lower_counts;
  }//if( log y-axis, and specifying approx less than zero counts )
  
  m_spectrum->setYAxisRange( lower_counts, upper_counts );
  
  return success;
}//bool setYAxisRange(...)



void InterSpec::displayedSpectrumRange( double &xmin, double &xmax, double &ymin, double &ymax ) const
{
  m_spectrum->visibleRange( xmin, xmax, ymin, ymax );
}

void InterSpec::handleShiftAltDrag( double lowEnergy, double upperEnergy )
{
  if( !m_gammaCountDialog && upperEnergy<=lowEnergy )
    return;

  showGammaCountDialog();
  m_gammaCountDialog->setEnergyRange( lowEnergy, upperEnergy );
}//void InterSpec::handleShiftAltDrag( double lowEnergy, double upperEnergy )



void InterSpec::searchForSinglePeak( const double x )
{
  if( !m_peakModel )
    throw runtime_error( "InterSpec::searchForSinglePeak(...): "
                        "shoudnt be called if peak model isnt set.");
  
  std::shared_ptr<SpecUtils::Measurement> data = m_spectrum->data();
  
  if( !m_dataMeasurement || !data )
    return;
  
  const double xmin = m_spectrum->xAxisMinimum();
  const double xmax = m_spectrum->xAxisMaximum();
  
  const double specWidthPx = m_spectrum->chartWidthInPixels();
  const double pixPerKeV = (xmax > xmin && xmax > 0.0 && specWidthPx > 10.0) ? std::max(0.001,(specWidthPx/(xmax - xmin))): 0.001;
  
  std::shared_ptr<const DetectorPeakResponse> det = m_dataMeasurement->detector();
  vector< PeakModel::PeakShrdPtr > origPeaks;
  if( !!m_peakModel->peaks() )
  {
    for( const PeakModel::PeakShrdPtr &p : *m_peakModel->peaks() )
    {
      origPeaks.push_back( p );
      
      //Avoid fitting a peak in the same area a data defined peak is.
      if( !p->gausPeak() && (x >= p->lowerX()) && (x <= p->upperX()) )
        return;
    }
  }//if( m_peakModel->peaks() )
  
  pair< PeakShrdVec, PeakShrdVec > foundPeaks;
  foundPeaks = searchForPeakFromUser( x, pixPerKeV, data, origPeaks );
  
  //cerr << "Found " << foundPeaks.first.size() << " peaks to add, and "
  //     << foundPeaks.second.size() << " peaks to remove" << endl;
  
  if( foundPeaks.first.empty()
      || foundPeaks.second.size() >= foundPeaks.first.size() )
  {
    char msg[256];
    snprintf( msg, sizeof(msg), "Couldn't find peak a peak near %.1f keV", x );
    passMessage( msg, "", 0 );
    return;
  }//if( foundPeaks.first.empty() )
  
  
  for( const PeakModel::PeakShrdPtr &p : foundPeaks.second )
    m_peakModel->removePeak( p );

  
  //We want to add all of the previously found peaks back into the model, before
  //  adding the new peak.
  PeakShrdVec peakstoadd( foundPeaks.first.begin(), foundPeaks.first.end() );
  PeakShrdVec existingpeaks( foundPeaks.second.begin(), foundPeaks.second.end() );
  
  //First add all the new peaks that have a nuclide/reaction/xray associated
  //  with them, since we know they are existing peaks
  for( const PeakModel::PeakShrdPtr &p : foundPeaks.first )
  {
    if( p->parentNuclide() || p->reaction() || p->xrayElement() )
    {
      //find nearest previously existing peak, and add new peak, while removing
      //  old one from existingpeaks
      int nearest = -1;
      double smallesdist = DBL_MAX;
      for( size_t i = 0; i < existingpeaks.size(); ++i )
      {
        const PeakModel::PeakShrdPtr &prev = existingpeaks[i];
        if( prev->parentNuclide() != p->parentNuclide() )
          continue;
        if( prev->reaction() != p->reaction() )
          continue;
        if( prev->xrayElement() != p->xrayElement() )
          continue;
        
        const double thisdif = fabs(p->mean() - prev->mean());
        if( thisdif < smallesdist )
        {
          nearest = static_cast<int>( i );
          smallesdist = thisdif;
        }
      }
      
      if( nearest >= 0 )
      {
        addPeak( *p, false );
        existingpeaks.erase( existingpeaks.begin() + nearest );
        peakstoadd.erase( std::find(peakstoadd.begin(), peakstoadd.end(), p) );
      }//if( nearest >= 0 )
    }//if( p->parentNuclide() || p->reaction() || p->xrayElement() )
  }//for( const PeakModel::PeakShrdPtr &p : peakstoadd )
  
  
  //Now go through and add the new versions of the previously existing peaks,
  //  using energy to match the previous to current peak.
  for( const PeakModel::PeakShrdPtr &p : existingpeaks )
  {
    size_t nearest = 0;
    double smallesdist = DBL_MAX;
    for( size_t i = 0; i < peakstoadd.size(); ++i )
    {
      const double thisdif = fabs(p->mean() - peakstoadd[i]->mean());
      if( thisdif < smallesdist )
      {
        nearest = i;
        smallesdist = thisdif;
      }
    }
    
    std::shared_ptr<const PeakDef> peakToAdd = peakstoadd[nearest];
    peakstoadd.erase( peakstoadd.begin() + nearest );
    addPeak( *peakToAdd, false );
  }//for( const PeakModel::PeakShrdPtr &p : existingpeaks )
  
  //Just in case we messed up the associations between the existing peak an
  //  their respective new version, we'll add all peaks that have a
  //  nuclide/reaction/xray associated with them, since they must be previously
  //  existing
  for( const PeakModel::PeakShrdPtr &p : peakstoadd )
    if( p->parentNuclide() || p->reaction() || p->xrayElement() )
      addPeak( *p, false );
  
  //Finally, in principle we will add the new peak here
  for( const PeakModel::PeakShrdPtr &p : peakstoadd )
  {
    if( !p->parentNuclide() && !p->reaction() && !p->xrayElement() )
      addPeak( *p, true );
  }
}//void searchForSinglePeak( const double x )


void InterSpec::automatedPeakSearchStarted()
{
  if( m_peakInfoDisplay )
    m_peakInfoDisplay->enablePeakSearchButton( false );
}//void automatedPeakSearchStarted()


void InterSpec::automatedPeakSearchCompleted()
{
  if( m_peakInfoDisplay )
    m_peakInfoDisplay->enablePeakSearchButton( true );
}//void automatedPeakSearchCompleted()


bool InterSpec::colorPeaksBasedOnReferenceLines() const
{
  return m_colorPeaksBasedOnReferenceLines;
}



void InterSpec::searchForHintPeaks( const std::shared_ptr<SpecMeas> &data,
                                         const std::set<int> &samples )
{
  cerr << "Starting searchForHintPeaks()" << endl;
  
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > origPeaks
                                                        = m_peakModel->peaks();
  if( !!origPeaks )
    origPeaks = std::make_shared<deque<PeakModel::PeakShrdPtr> >( *origPeaks );
  
  std::shared_ptr< vector<std::shared_ptr<const PeakDef> > > searchresults
            = std::make_shared< vector<std::shared_ptr<const PeakDef> > >();
  
  std::weak_ptr<const SpecUtils::Measurement> weakdata = m_spectrum->data();
  std::weak_ptr<SpecMeas> spectrum = data;
  
  boost::function<void(void)> callback = wApp->bind(
                boost::bind(&InterSpec::setHintPeaks,
                this, spectrum, samples, origPeaks, searchresults) );
  
  boost::function<void(void)> worker = boost::bind( &PeakSearchGuiUtils::search_for_peaks_worker,
                                                   weakdata, origPeaks,
                                                   vector<ReferenceLineInfo>(), false,
                                                   searchresults,
                                                   callback, wApp->sessionId(), true );

  if( m_findingHintPeaks )
  {
    m_hintQueue.push_back( worker );
  }else
  {
    Wt::WServer *server = Wt::WServer::instance();
    if( server )  //this should always be true
    {
      m_findingHintPeaks = true;
      server->ioService().boost::asio::io_service::post( worker );
    }//if( server )
  }
}//void searchForHintPeaks(...)


void InterSpec::setHintPeaks( std::weak_ptr<SpecMeas> weak_spectrum,
                  std::set<int> samplenums,
                  std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > existing,
                  std::shared_ptr<std::vector<std::shared_ptr<const PeakDef> > > resultpeaks )
{
  //cerr << "InterSpec::setHintPeaks(...) with "
  //     << (!!resultpeaks ? resultpeaks->size() : size_t(0)) << " peaks." << endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
  if( !wApp )
    log_developer_error( __func__, "setHintPeaks() being called from not within the event loop!" );
#endif

  m_findingHintPeaks = false;
  
  if( m_hintQueue.size() )
  {
    Wt::WServer *server = Wt::WServer::instance();
    if( server )  //this should always be true
    {
      m_findingHintPeaks = true;
      cerr << "InterSpec::setHintPeaks(...): posting queued job" << endl;
      boost::function<void()> worker = m_hintQueue.back();
      m_hintQueue.pop_back();
      server->ioService().boost::asio::io_service::post( worker );
    }//if( server )
  }//if( m_hintQueue.size() )
  
  typedef std::shared_ptr<const PeakDef> PeakPtr;
  typedef deque< PeakPtr > PeakDeque;
  std::shared_ptr<SpecMeas> spectrum = weak_spectrum.lock();
  
  if( !spectrum || !resultpeaks )
  {
    cerr << "InterSpec::setHintPeaks(): invalid SpecMeas" << endl;
    return;
  }//if( !spectrum )
  
//  if( spectrum != m_dataMeasurement && spectrum != m_backgroundMeasurement
//      && spectrum != m_secondDataMeasurement )
//  {
//    cerr << "InterSpec::setHintPeaks(): SpecMeas not current spectrum"
//         << endl;
//    return;
//  }
  
  //we could check to see if the spectrum and sample numbers are still the
  //  current one.  If not the user probably doesnt care about this spectrum,
  //  so why bother storing ther results?
  
  std::shared_ptr< PeakDeque > newpeaks
    = std::make_shared<PeakDeque>( resultpeaks->begin(), resultpeaks->end() );
  
  //See if the user has added any peaks since we did the automated search
  std::shared_ptr<PeakDeque> current_user_peaks = spectrum->peaks( samplenums );
  
  vector< PeakPtr > addedpeaks;
  if( !!current_user_peaks && !existing )
  {
    addedpeaks.insert( addedpeaks.end(), current_user_peaks->begin(), current_user_peaks->end() );
  }else if( !!current_user_peaks && !!existing )
  {
    for( const PeakPtr &p : *current_user_peaks )
      if( std::find(existing->begin(),existing->end(),p) != existing->end() )
        addedpeaks.push_back( p );
  }
  
  if( addedpeaks.size() )
  {
    for( PeakPtr p : addedpeaks )
    {
      const int pos = add_hint_peak_pos( p, *newpeaks );
      if( pos >= 0 )
        newpeaks->insert( newpeaks->begin() + pos, p );
    }
  }//if( addedpeaks.size() )
  
  spectrum->setAutomatedSearchPeaks( samplenums, newpeaks );
//  existing
}//void setHintPeaks(...)




/*
 //Depreciated 20150204 by wcjohns in favor of calling 
 //  InterSpec::findPeakFromControlDrag()
void InterSpec::findPeakFromUserRange( double x0, double x1 )
{
  if( !m_peakModel )
    throw runtime_error( "SpectrumDisplayDiv::findPeakFromUserRange(...): "
                        "shoudnt be called if peak model isnt set.");
  
  std::shared_ptr<SpecUtils::Measurement> data = m_spectrum->data();
  
  if( !data )
  {
    passMessage( "No spectrum data", "", WarningWidget::WarningMsgMedium );
    return;
  }
    
  if( x0 > x1 )
    swap( x0, x1 );
  
  //round to the edges of the bins the drag action was from
  const size_t lowerchannel = data->find_gamma_channel( x0 );
  const size_t upperchannel = data->find_gamma_channel( x1 );
  x0 = data->gamma_channel_lower( lowerchannel );
  x1 = data->gamma_channel_upper( upperchannel );
  
  PeakDef peak( 0.0, 0.0, 0.0 );
  peak.setMean( 0.5*(x0+x1) );
  peak.setSigma( (x1-x0)/8.0 );
  std::shared_ptr<PeakContinuum> continuum = peak.continuum();
  
  continuum->calc_linear_continuum_eqn( data, x0, x1, 1 );
  
  const double data_area = data->gamma_channels_sum( lowerchannel, upperchannel );
  const double continuum_area = peak.offset_integral( x0, x1 );
  peak.setAmplitude( data_area - continuum_area );
  
  std::vector<PeakDef> peakV;
  peakV.push_back( peak );
  
  const double stat_threshold = 0.5;
  const double hypothesis_threshold = 0.5;
  
  //XXX - we promised the user to use a given continuum, however the below only
  //      uses this as a starting point
  peakV = fitPeaksInRange( x0+0.00001, x1-0.00001, 0,
                           stat_threshold, hypothesis_threshold,
                           peakV, data, m_peakModel->peakVec() );
  
  for( size_t i = 0; i < peakV.size(); ++i )
    addPeak( peakV[i], true );
    
  if( peakV.size()==0 )
    passMessage( "No peaks were found. You might try a slighlty changed ROI "
                 "region.", "", WarningWidget::WarningMsgLow );
}//void findPeakFromUserRange( const double x0, const double x1 )
*/

void InterSpec::excludePeaksFromRange( double x0, double x1 )
{
  if( !m_peakModel )
    throw runtime_error( "InterSpec::excludePeaksFromRange(...): "
                        "shoudnt be called if peak model isnt set.");
  if( x0 > x1 )
    swap( x0, x1 );
  
  vector<PeakDef> all_peaks = m_peakModel->peakVec();
//  vector<PeakDef> peaks_in_range = peaksInRange( x0, x1, 4.0, all_peaks );
  vector<PeakDef> peaks_in_range = peaksTouchingRange( x0, x1, all_peaks );
  
  if( peaks_in_range.empty() )
    return;
  
  std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
  if( !data )
    return;
  
  for( const PeakDef &peak : peaks_in_range )
  {
    vector<PeakDef>::iterator iter = std::find( all_peaks.begin(),
                                                all_peaks.end(), peak );
    if( iter != all_peaks.end() )
      all_peaks.erase( iter );
  }//for( peak in range ...)
  
  
  vector<PeakDef> peaks_to_keep;
  double minEffectedPeak = DBL_MAX, maxEffectedPeak = -DBL_MAX;
  
  for( PeakDef peak : peaks_in_range )
  {
    if( peak.mean()>=x0 && peak.mean()<=x1 )
      continue;
    
    std::shared_ptr<PeakContinuum> continuum = peak.continuum();
    
    double lowx = 0.0, upperx = 0.0;
    findROIEnergyLimits( lowx, upperx, peak, data );
    
    minEffectedPeak = std::min( minEffectedPeak, peak.mean() );
    maxEffectedPeak = std::max( maxEffectedPeak, peak.mean() );

    if( x0>=upperx || x1<=lowx )
    {
      peaks_to_keep.push_back( peak );
      continue;
    }
    
    const bool x0InPeak = ((x0>=lowx) && (x0<=upperx));
    const bool x1InPeak = ((x1>=lowx) && (x1<=upperx));
    
    if( !x0InPeak && !x1InPeak )
    {
      peaks_to_keep.push_back( peak );
      continue;
    }else if( x0InPeak && x1InPeak )
    {
      if( x0>peak.mean() )
        upperx = x0;
      else
        lowx = x1;
    }else if( x0InPeak )
    {
      upperx = x0;
    }else //if( x1InPeak )
    {
      lowx = x1;
    }//if / else

    continuum->setRange( lowx, upperx );
    
    peaks_to_keep.push_back( peak );
  }//for( PeakDef peak : peaks_in_range )
  
  
  std::shared_ptr<const DetectorPeakResponse> detector
                                          = measurment(SpecUtils::SpectrumType::Foreground)->detector();
  map<std::shared_ptr<const PeakContinuum>,PeakShrdVec > peaksinroi;
  for( PeakDef peak : peaks_to_keep )
    peaksinroi[peak.continuum()].push_back( std::make_shared<PeakDef>(peak) );
  
  map<std::shared_ptr<const PeakContinuum>, PeakShrdVec >::const_iterator iter;
  for( iter = peaksinroi.begin(); iter != peaksinroi.end(); ++iter )
  {
    //Instead of doing the fitting here, we could just:
    //m_rightClickEnergy = iter->second[0].mean();
    //refitPeakFromRightClick();
    //This would make things more consistent between right clicking to fit, and
    //  erasing part of the ROI (which is what happened if we are here).
    
    PeakShrdVec newpeaks = refitPeaksThatShareROI( data, detector, iter->second, 0.25 );
    if( newpeaks.size() == iter->second.size() )
    {
      for( size_t j = 0; j < newpeaks.size(); ++j )
        all_peaks.push_back( PeakDef(*newpeaks[j]) );
      std::sort( all_peaks.begin(), all_peaks.end(), &PeakDef::lessThanByMean );
    }else
    {
      //if newpeaks.empty(), then the Chi2 of the fit probably did not improve,
      //  so we'll just use the existing peaks.
      //I dont actually know how I feel about this: on one hand I want to
      //  believe in the Chi2, but on the other, we should probably respond to
      //  the user chaning the ROI. refitPeakFromRightClick() will then attempt
      //  a second fitting methos
      if( newpeaks.size() )
        cerr << "Warning: failed to refit peaks for some reason; got "
             << newpeaks.size()  << " with an input of " << iter->second.size()
             << endl;
      for( size_t j = 0; j < iter->second.size(); ++j )
        all_peaks.push_back( PeakDef(*iter->second[j]) );
      
      std::sort( all_peaks.begin(), all_peaks.end(), &PeakDef::lessThanByMean );
      
      /*
      double stat_threshold = 0.5, hypothesis_threshold = 0.5;
      vector<PeakDef> newpeaks;
      newpeaks = fitPeaksInRange( minEffectedPeak,
                                 maxEffectedPeak,
                                 3.0, stat_threshold, hypothesis_threshold,
                                 all_peaks,
                                 m_spectrum->data(),
                                 vector<PeakDef>() );
      if( all_peaks.size() == newpeaks.size() )
        all_peaks = newpeaks;
       */
    }
  }//for( loop over peaksinroi elements )
  
  
  m_peakModel->setPeaks( all_peaks );
}//void excludePeaksFromRange( const double x0, const double x1 )


void InterSpec::guessIsotopesForPeaks( WApplication *app )
{
  InterSpec *viewer = this;
  if( !m_peakModel )
    throw runtime_error( "InterSpec::refitPeakAmplitudes(...): "
                         "shoudnt be called if peak model isnt set.");
  
  std::shared_ptr<const SpecUtils::Measurement> data;
  std::shared_ptr<const DetectorPeakResponse> detector;
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > all_peaks;
  
  {//begin code-block to get input data
    std::unique_ptr<WApplication::UpdateLock> applock;
    if( app )
      applock.reset( new WApplication::UpdateLock(app) );
    
    if( app && m_peakModel->peaks() )
      all_peaks = std::make_shared<deque<PeakModel::PeakShrdPtr> >( *m_peakModel->peaks() );
    else
      all_peaks = m_peakModel->peaks();
    
    if( !all_peaks || all_peaks->empty() )
      return;
    
    data = m_spectrum->data();
  
    if( viewer && viewer->measurment(SpecUtils::SpectrumType::Foreground) )
    {
      detector = viewer->measurment(SpecUtils::SpectrumType::Foreground)->detector();
      
      if( detector )
      {
        DetectorPeakResponse *resp = new DetectorPeakResponse( *detector );
        detector.reset( resp );
      }//if( detector )
    }//if( viewer && viewer->measurment(SpecUtils::SpectrumType::Foreground) )
  
    if( detector && !detector->hasResolutionInfo() )
    {
      try
      {
        std::shared_ptr<DetectorPeakResponse> newdetector = std::make_shared<DetectorPeakResponse>( *detector );
        newdetector->fitResolution( all_peaks, data, DetectorPeakResponse::kGadrasResolutionFcn );
        detector = newdetector;
      }catch( std::exception & )
      {
        detector.reset();
      }
    }
    
    if( !detector || !detector->isValid() )
    {
      DetectorPeakResponse *detPtr = new DetectorPeakResponse();
      detector.reset( detPtr );
    
      string drf_dir;
      if( PeakFitUtils::is_high_res(data) )
        drf_dir = SpecUtils::append_path(sm_staticDataDirectory, "GenericGadrasDetectors/HPGe 40%" );
      else
        drf_dir = SpecUtils::append_path(sm_staticDataDirectory, "GenericGadrasDetectors/NaI 1x1" );
      
      detPtr->fromGadrasDirectory( drf_dir );
    }//if( !detector || !detector->isValid() )
  }//end code-block to get input data
  
  vector<IsotopeId::PeakToNuclideMatch> idd( all_peaks->size() );
  
  size_t peakn = 0, rownum = 0;
  vector<size_t> rownums;
  SpecUtilsAsync::ThreadPool threadpool;
//  vector< boost::function<void()> > workers;
  vector<PeakModel::PeakShrdPtr> inputpeaks;
  vector<PeakDef> modifiedPeaks;
  for( PeakModel::PeakShrdPtr peak : *all_peaks )
  {
    ++rownum;
    inputpeaks.push_back( peak );
    
    if( (peak->parentNuclide()
         && (peak->nuclearTransition()
             || (peak->sourceGammaType()==PeakDef::AnnihilationGamma) ))
        || peak->reaction() || peak->xrayElement() )
    {
      modifiedPeaks.push_back( *peak );
      continue;
    }
    
    threadpool.post( boost::bind( &IsotopeId::suggestNuclides,
                                  boost::ref(idd[peakn]), peak, all_peaks,
                                  data, detector ) );
//    workers.push_back( boost::bind( &suggestNuclides, boost::ref(idd[peakn]),
//                                   peak, all_peaks, data, detector ) ) );
    rownums.push_back( rownum-1 );
    ++peakn;
  }//for( PeakModel::PeakShrdPtr peak : *all_peaks )
  
  threadpool.join();
//  SpecUtils::do_asyncronous_work( workers, false );

    
  for( size_t resultnum = 0; resultnum < rownums.size(); ++resultnum )
  {
    const size_t row = rownums[resultnum];
    const PeakModel::PeakShrdPtr peak = m_peakModel->peakPtr( row ); 
    PeakDef newPeak = *inputpeaks[row];
      
    if( newPeak.parentNuclide() || newPeak.reaction() || newPeak.xrayElement() )
    {
      modifiedPeaks.push_back( newPeak );
      continue;
    }//if( !isotopeData.empty() )
      
    const IsotopeId::PeakToNuclideMatch match = idd[resultnum];
    vector<PeakDef::CandidateNuclide> candidates;
      
    for( size_t j = 0; j < match.nuclideWeightPairs.size(); ++j )
    {
      const IsotopeId::NuclideStatWeightPair &p = match.nuclideWeightPairs[j];
      
      PeakDef::SourceGammaType sourceGammaType = PeakDef::NormalGamma;
      size_t radparticleIndex;
      const SandiaDecay::Transition *transition = NULL;
      const double sigma = newPeak.gausPeak() ? newPeak.sigma() : 0.125*newPeak.roiWidth();
      PeakDef::findNearestPhotopeak( p.nuclide, match.energy,
                                       4.0*sigma, false, transition,
                                       radparticleIndex, sourceGammaType );
      if( j == 0 )
        newPeak.setNuclearTransition( p.nuclide, transition,
                                      int(radparticleIndex), sourceGammaType );
        
      if( transition || (sourceGammaType==PeakDef::AnnihilationGamma) )
      {
        PeakDef::CandidateNuclide candidate;
        candidate.nuclide          = p.nuclide;
        candidate.weight           = p.weight;
        candidate.transition       = transition;
        candidate.sourceGammaType  = sourceGammaType;
        candidate.radparticleIndex = static_cast<int>(radparticleIndex);
        candidates.push_back( candidate );
      }//if( transition )
    }//for( size_t j = 0; j < match.nuclideWeightPairs.size(); ++j )

    newPeak.setCandidateNuclides( candidates );
    modifiedPeaks.push_back( newPeak );
  }//for( size_t i = 0; i < idd.size(); ++i )

  {//begin codeblock to set modified peaks
    std::unique_ptr<WApplication::UpdateLock> applock;
    if( app )
      applock.reset( new WApplication::UpdateLock(app) );
    
    m_peakModel->setPeaks( modifiedPeaks );
    
    if( app )
      app->triggerUpdate();
  }//end codeblock to set modified peaks
}//void guessIsotopesForPeaks()


vector<pair<float,int> > InterSpec::passthroughTimeToSampleNumber() const
{
  if( !m_dataMeasurement )
    return {};
  
  const vector<string> disp_dets = detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  
  /* If we have both "derived" data and non-derived data, we dont want to count the derived data
     spectra as background, because then weird things happen.
   */
  const bool hasDerivedData = m_dataMeasurement->contains_derived_data();
  const bool hasNonDerivedData = m_dataMeasurement->contains_non_derived_data();
  
  // We'll first grab foreground, background, and derived samples separately, then combine them
  //  So background samples will be first, and may be compressed, then foreground, then derived.
  double foretime = 0.0, backtime = 0.0, derivedtime = 0.0;
  vector<pair<float,int> > foreground, background, derived_data;
  
  const set<int> all_sample_nums = m_dataMeasurement->sample_numbers();
  for( const int sample_num : all_sample_nums )
  {
    if( m_excludedSamples.count(sample_num) )
      continue;
    
    const float sampletime = sample_real_time_increment( m_dataMeasurement, sample_num, disp_dets );
    if( sampletime <= 0.0f )
      continue;
    
    bool isFore = false, isback = false, isDerived = false;
    const auto meass = m_dataMeasurement->sample_measurements(sample_num);
    
    for( const std::shared_ptr<const SpecUtils::Measurement> &m : meass )
    {
      if( hasDerivedData && hasNonDerivedData && m->derived_data_properties() )
      {
        isDerived = true;
        continue;
      }
      
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          break;
          
        case SpecUtils::SourceType::Background:
          isback = true;
          break;
          
        case SpecUtils::SourceType::Foreground:
        case SpecUtils::SourceType::Unknown:
          isFore = true;
          break;
      }//switch( m->source_type() )
    }//for( const std::shared_ptr<const SpecUtils::Measurement> &m : meass )
    
    if( isDerived )
    {
      derived_data.push_back( {sampletime, sample_num} );
      derivedtime += sampletime;
    }else if( isFore )
    {
      foreground.push_back( {sampletime, sample_num} );
      foretime += sampletime;
    }else if( isback )
    {
      background.push_back( {sampletime, sample_num} );
      backtime += sampletime;
    }
  }//for( const int sample_num : sample_nums )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
// Note: this check is not always valid.  The code above will put a sample as foreground if
//       and of its measurements are foreground/unknown, but validForegroundSamples() will
//       remove the sample from foreground if any of its measurements are background or cal.
//  const auto prevfore = validForegroundSamples();
//  set<int> newfore;
//  for( auto s : foreground )
//    newfore.insert( s.second );
//  assert( newfore == prevfore );
#endif
  

  double cumulative_time = 0.0;
  vector<pair<float,int> > answer;
  if( background.empty() && foreground.empty() && derived_data.empty() )
    return answer;
  
  answer.reserve( foreground.size() + background.size() + derived_data.size() + 1 );
  
  if( !background.empty() )
  {
    // We want to limit the background samples to take up about 10% of the time chart because we
    //  normally dont care much about the background variation, and a lot of times there can be like
    //  a 5 minute, single spectrum, background, and like 10 seconds foreground, which would make
    //  the foreground not even visible.
    double backscale = 1.0;
    if( !foreground.empty() && (foretime > 0.1) && (backtime > 0.1*foretime) )
      backscale = ( std::ceil(0.1*foretime) ) / backtime;
    
    cumulative_time = -backscale * backtime;
    
    for( const auto &time_sample : background )
    {
      answer.push_back( {static_cast<float>(cumulative_time), time_sample.second} );
      cumulative_time += (backscale * time_sample.first);
    }
  }//if( !background.empty() )
  
  assert( fabs(cumulative_time) < 1.0E-4 );
  
  for( const auto &time_sample : foreground )
  {
    answer.push_back( {static_cast<float>(cumulative_time), time_sample.second} );
    cumulative_time += time_sample.first;
  }
  
  for( const auto &time_sample : derived_data )
  {
    answer.push_back( {cumulative_time, time_sample.second} );
    cumulative_time += time_sample.first;
  }

  // Add in upper edge of last time segment.
  answer.push_back( {cumulative_time, answer.back().second + 1} );
  
  return answer;
}//vector<std::pair<float,int> > passthroughTimeToSampleNumber() const


void InterSpec::displayTimeSeriesData()
{
  if( m_dataMeasurement && m_dataMeasurement->passthrough() )
  {
    if( m_timeSeries->isHidden() )
    {
      m_timeSeries->setHidden( false );
      m_chartResizer->setHidden( m_timeSeries->isHidden() );
    }//if( m_timeSeries->isHidden() )
    
    const vector<string> det_to_use = detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    m_timeSeries->setData( m_dataMeasurement, det_to_use );
    
    const set<int> emptyset;
    const set<int> &fore   = m_displayedSamples;
    const set<int> &back   = (m_backgroundMeasurement == m_dataMeasurement)
                                ? m_backgroundSampleNumbers : emptyset;
    const set<int> &second = (m_secondDataMeasurement == m_dataMeasurement)
                                ? m_sectondForgroundSampleNumbers : emptyset;
  
    m_timeSeries->setHighlightedIntervals( fore, SpecUtils::SpectrumType::Foreground );
    m_timeSeries->setHighlightedIntervals( back, SpecUtils::SpectrumType::Background );
    m_timeSeries->setHighlightedIntervals( second, SpecUtils::SpectrumType::SecondForeground );
  }else
  {
    m_timeSeries->setData( nullptr, {} );
    m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::Foreground );
    m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::Background );
    m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::SecondForeground );
    
    if( !m_timeSeries->isHidden() )
    {
      m_timeSeries->setHidden( true );
      m_chartResizer->setHidden( m_timeSeries->isHidden() );
    }//if( !m_timeSeries->isHidden() )
  }//if( passthrough ) / else
}//void displayTimeSeriesData()


std::set<int> InterSpec::sampleRangeToSet( int start_sample,  int end_sample,
                                std::shared_ptr<const SpecMeas> meas,
                                const std::set<int> &excluded_samples )
{
  std::set<int> display_samples;
  if( meas )
  {
    set<int> sample_numbers = meas->sample_numbers();
    for( int s : excluded_samples )
      sample_numbers.erase( s );

    const int first_sample = (*sample_numbers.begin());
    const int last_sample = (*sample_numbers.rbegin());
  //  const int num_samples = static_cast<int>( meas->sample_numbers().size() );

    if( start_sample < first_sample ) start_sample = first_sample;
    if( end_sample < 0 )              end_sample = last_sample;
    if( end_sample > last_sample )    end_sample = last_sample;

    if( sample_numbers.size() == 1 )
      start_sample = end_sample = first_sample;

    for( int sample : sample_numbers )
      if( sample >= start_sample && sample <= end_sample )
        display_samples.insert( sample );
  }//if( meas )

  return display_samples;
}//sampleRangeToSet



vector<string> InterSpec::detectorsToDisplay( const SpecUtils::SpectrumType type ) const
{
  shared_ptr<const SpecMeas> wanted_meas = measurment(type);
  if( !wanted_meas )
    return vector<string>{};
  
  if( !m_dataMeasurement )
    return wanted_meas->detector_names();
  
  const vector<string> &all_det_names = m_dataMeasurement->detector_names();
  
  if( (type != SpecUtils::SpectrumType::Foreground)
      && (wanted_meas->detector_names() != all_det_names) )
    return wanted_meas->detector_names();
  
  vector<string> detector_names;
  
  if( !m_detectorToShowMenu || all_det_names.empty() )
    return all_det_names;
  
  const vector<WMenuItem *> items = m_detectorToShowMenu->items();
  
  size_t detnum = 0;
  for( size_t i = 0; i < items.size(); ++i )
  {
    if( items[i]->hasStyleClass("PhoneMenuBack") )
      continue;
    
#if( WT_VERSION>=0x3030300 )
    PopupDivMenuItem *item = dynamic_cast<PopupDivMenuItem *>( items[i] );
    WCheckBox *cb = (item ? item->checkBox() : (WCheckBox *)0);
#else
    WCheckBox *cb = 0;
    for( int j = 0; !cb && (j < items[i]->count()); ++j )
      cb = dynamic_cast<WCheckBox *>(items[i]->widget(j));
#endif
    
    if( cb && detnum < all_det_names.size() )
    {
      if( cb->isChecked() )
        detector_names.push_back( all_det_names[detnum] );
    }else
    {
      cerr << "XXX - Failed to get checkbox in menu!" << endl;
    }
    ++detnum;
  }// for( size_t i = 0; i < items.size(); ++i )
  
  if( detnum != all_det_names.size() )
    throw runtime_error( "InterSpec::detectorsToDisplay(): serious logic error" );

  return detector_names;
}//std::vector<string> detectorsToDisplay() const


void InterSpec::refreshDisplayedCharts()
{
  //We want to keep the old display scale factors, but calling displayForegroundData() will
  //  reset them.
  const float backSf = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::Background);
  const float secondSf = m_spectrum->displayScaleFactor(SpecUtils::SpectrumType::SecondForeground);
  
  if( m_dataMeasurement )
    displayForegroundData( true );
  
  if( m_secondDataMeasurement )
  {
    displaySecondForegroundData();
    m_spectrum->setDisplayScaleFactor( secondSf, SpecUtils::SpectrumType::SecondForeground );
  }//if( m_secondDataMeasurement )
  
  if( m_backgroundMeasurement )
  {
    displayBackgroundData();
    m_spectrum->setDisplayScaleFactor( backSf, SpecUtils::SpectrumType::Background );
  }//if( m_backgroundMeasurement )
  
  // Now check if the currently displayed energy range extend past all the ranges of the data.
  //  If so, dont display past where the data is.
  double min_energy = 99999.9, max_energy = -99999.9;

  auto checkEnergyRange = [&min_energy, &max_energy]( std::shared_ptr<const SpecUtils::Measurement> m ){
    auto cal = m ? m->energy_calibration() : nullptr;
    if( !cal || !cal->valid() )
      return;
    min_energy = std::min( static_cast<double>( cal->lower_energy() ), min_energy );
    max_energy = std::max( static_cast<double>( cal->upper_energy() ), max_energy );
  };
  
  checkEnergyRange( m_spectrum->data() );
  checkEnergyRange( m_spectrum->background() );
  checkEnergyRange( m_spectrum->secondData() );
  
  if( (max_energy > min_energy)
     && !IsInf(min_energy) && !IsInf(max_energy)
     && !IsNan(min_energy) && !IsNan(max_energy) )
  {
    const double curr_min = m_spectrum->xAxisMinimum();
    const double curr_max = m_spectrum->xAxisMaximum();
    
    if( (max_energy < curr_max) || (curr_min < min_energy) )
    {
      min_energy = ((curr_min < min_energy) || (curr_min >= max_energy)) ? min_energy : curr_min;
      max_energy = ((max_energy < curr_max) || (curr_max <= min_energy)) ? max_energy : curr_max;
      
      m_spectrum->setXAxisRange( min_energy, max_energy );
    }//if( either min or max range is outside of the data )
  }//if( possible energy range is valid )
}//void refreshDisplayedCharts()


std::set<int> InterSpec::sampleNumbersForTypeFromForegroundFile( const SpecUtils::SpectrumType type ) const
{
  set<int> dispsamples;
  if( !m_dataMeasurement )
    return dispsamples;
  
  for( const auto &m : m_dataMeasurement->measurements() )
  {
    switch( m->source_type() )
    {
      case SpecUtils::SourceType::IntrinsicActivity:
      case SpecUtils::SourceType::Calibration:
        if( type == SpecUtils::SpectrumType::SecondForeground )
          dispsamples.insert( m->sample_number() );
        break;
        
      case SpecUtils::SourceType::Background:
        if( type == SpecUtils::SpectrumType::Background )
          dispsamples.insert( m->sample_number() );
        break;
      case SpecUtils::SourceType::Foreground:
      case SpecUtils::SourceType::Unknown:
        if( type == SpecUtils::SpectrumType::Foreground )
          dispsamples.insert( m->sample_number() );
        break;
    }//switch( m->source_type() )
  }//for( loop over all measurements )
  
  return dispsamples;
}//set<int> sampleNumbersForType( SpectrumType )


void InterSpec::displayForegroundData( const bool current_energy_range )
{
  auto &meas = m_dataMeasurement;
  set<int> &sample_nums = m_displayedSamples;
  const vector<string> detectors = detectorsToDisplay( SpecUtils::SpectrumType::Foreground );
  
  if( meas && !detectors.empty() && sample_nums.empty() )
   {
     sample_nums = validForegroundSamples();
     if( !meas->passthrough() && (sample_nums.size() > 1) )
       sample_nums = { *begin(sample_nums) };
   }
  
  if( !meas || detectors.empty() || sample_nums.empty() )
  {
    m_backgroundSubItems[0]->disable();
    m_backgroundSubItems[0]->show();
    m_backgroundSubItems[1]->hide();
    m_hardBackgroundSub->disable();
    
    m_showYAxisScalerItems[0]->setDisabled( m_showYAxisScalerItems[0]->isVisible() );
    m_showYAxisScalerItems[1]->setDisabled( m_showYAxisScalerItems[1]->isVisible() );
    
    if( m_spectrum->data() )
    {
      m_spectrum->setData( nullptr, false );
      m_peakModel->setPeakFromSpecMeas( nullptr, sample_nums );
    }
    
    return;
  }//if( !meas )


  m_peakModel->setPeakFromSpecMeas( meas, sample_nums );

  const auto energy_cal = meas->suggested_sum_energy_calibration(sample_nums, detectors);
  if( !energy_cal )
  {
    vector<shared_ptr<const SpecUtils::Measurement>> meass;
    bool allNeutron = true, containSpectrum = false;
    for( const auto sample_num : sample_nums )
    {
      for( const auto &m : meas->sample_measurements( sample_num ) )
      {
         if( std::find( begin(detectors), end(detectors), m->detector_name() ) != end(detectors) )
         {
           meass.push_back( m );
           
           const auto cal = m->energy_calibration();
           const bool hasGamma = (m->num_gamma_channels() > 0);
           const bool hasSpectrum = cal && cal->valid() && (cal->num_channels() > 7);
           allNeutron = (allNeutron && !hasGamma);
           containSpectrum = (containSpectrum || hasSpectrum);
         }//if( this Measurement is from a detector we want )
      }//for( loop over Measurements for this sample number )
    }//for( const auto sample_num : sample_nums )
    
    string msg;
    if( meass.empty() )
    {
      if( detectors.size() == meas->detector_names().size() )
      {
        msg = "<p>The spectrum file didn't contain any spectra with the current sample numbers.</p>"
              "<p>Please select different/more sample numbers.</p>";
      }else
      {
        msg = "<p>The spectrum file didn't contain any spectra with the current sample numbers and"
              " detector names.</p>"
              "<p>Please select more detectors or sample numbers.</p>";
      }
    }else if( allNeutron )
    {
      msg = "<p>The current sample numbers and detector names only contain neutron data.</p>"
            "<p>Please select samples that include spectroscopic data.</p>";
    }else if( !containSpectrum )
    {
      msg = "<p>The current sample numbers and detector names do not include any spectroscopic"
            " measurements.</p>"
            "<p>Please select samples that include spectroscopic data.</p>";
    }else
    {
      msg = "<p>I couldn't determine binning to display the spectrum.</p>";
      
      if( meass.size() > 1 )
        msg += "<p>You might try selecting only a single spectra to display in "
        "the <b>File manager</b>, or a sample with gamma counts.</p>";
    }//if(
    
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  }//if( !binning )


  const bool canSub = (m_dataMeasurement && m_backgroundMeasurement);
  const bool isSub = m_spectrum->backgroundSubtract();
  if( m_backgroundSubItems[0]->isEnabled() != canSub )
    m_backgroundSubItems[0]->setDisabled( !canSub );
  if( m_backgroundSubItems[0]->isHidden() != isSub )
  {
    m_backgroundSubItems[0]->setHidden( isSub );
    m_backgroundSubItems[1]->setHidden( !isSub );
  }//if( m_backgroundSubItems[0]->isHidden() != isSub )

  if( m_hardBackgroundSub->isEnabled() != canSub )
    m_hardBackgroundSub->setDisabled( !canSub );
  
  std::shared_ptr<SpecUtils::Measurement> dataH;
  
  if( energy_cal )
    dataH = m_dataMeasurement->sum_measurements( sample_nums, detectors, energy_cal );
  
  if( dataH )
    dataH->set_title( "Foreground" );

  const float lt = dataH ? dataH->live_time() : 0.0f;
  const float rt = dataH ? dataH->real_time() : 0.0f;
  const float sum_neut = dataH ? dataH->neutron_counts_sum() : 0.0f;
  
  m_spectrum->setData( dataH, current_energy_range );
  
  if( !m_timeSeries->isHidden() )
    m_timeSeries->setHighlightedIntervals( sample_nums, SpecUtils::SpectrumType::Foreground );
}//void displayForegroundData()


void InterSpec::displaySecondForegroundData()
{
  const auto &meas = m_secondDataMeasurement;
  std::set<int> &sample_nums = m_sectondForgroundSampleNumbers;
  const auto disp_dets = detectorsToDisplay( SpecUtils::SpectrumType::SecondForeground );
  
  //Note: below will throw exception if 'disp_samples' has any invalid entries
  shared_ptr<const SpecUtils::EnergyCalibration> energy_cal;
  if( meas )
    energy_cal = meas->suggested_sum_energy_calibration( sample_nums, disp_dets );
  
  if( !energy_cal || !m_dataMeasurement )
  {
    //sample_nums.clear();
    if( m_spectrum->secondData() )
      m_spectrum->setSecondData( nullptr );
    
    if( !m_timeSeries->isHidden() )
      m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::SecondForeground );
    
    return;
  }//if( !m_secondDataMeasurement )
  
  if( !meas->num_measurements() )
    throw runtime_error( "Serious logic error in InterSpec::displaySecondForegroundData()" );

  auto histH = meas->sum_measurements( sample_nums, disp_dets, energy_cal );
  if( histH )
    histH->set_title( "Second Foreground" );
    
  m_spectrum->setSecondData( histH );
  
  if( !m_timeSeries->isHidden() )
    m_timeSeries->setHighlightedIntervals( m_sectondForgroundSampleNumbers, SpecUtils::SpectrumType::SecondForeground );
}//void displaySecondForegroundData()


void InterSpec::displayBackgroundData()
{
  const auto &meas = m_backgroundMeasurement;
  set<int> &disp_samples = m_backgroundSampleNumbers;
  
  const vector<string> disp_dets = detectorsToDisplay( SpecUtils::SpectrumType::Background );
  
  //Note: below will throw exception if 'disp_samples' has any invalid entries
  shared_ptr<const SpecUtils::EnergyCalibration> energy_cal;
  if( meas )
    energy_cal = meas->suggested_sum_energy_calibration( disp_samples, disp_dets );
  
  if( !energy_cal || !m_dataMeasurement )
  {
    m_backgroundSubItems[0]->disable();
    m_backgroundSubItems[0]->show();
    m_backgroundSubItems[1]->hide();
    m_hardBackgroundSub->disable();
    //disp_samples.clear();
    if( m_spectrum->background() )
      m_spectrum->setBackground( nullptr );
    
    if( !m_timeSeries->isHidden() )
      m_timeSeries->setHighlightedIntervals( {}, SpecUtils::SpectrumType::Background );
    
    return;
  }//if( !energy_cal || !m_dataMeasurement )
  
  auto backgroundH = meas->sum_measurements( disp_samples, disp_dets, energy_cal );
  if( backgroundH )
    backgroundH->set_title( "Background" );
    
  const float lt = backgroundH ? backgroundH->live_time() : -1.0f;
  const float rt = backgroundH ? backgroundH->real_time() : -1.0f;
  const float neutronCounts = backgroundH ? backgroundH->neutron_counts_sum() : -1.0f;
  m_spectrum->setBackground( backgroundH );
  
  if( !m_timeSeries->isHidden() )
  {
    const auto background = SpecUtils::SpectrumType::Background;
    if( m_backgroundMeasurement != m_dataMeasurement )
      m_timeSeries->setHighlightedIntervals( {}, background );
    else
      m_timeSeries->setHighlightedIntervals( m_backgroundSampleNumbers, background );
  }//if( !m_timeSeries->isHidden() )
  
  const bool canSub = (m_dataMeasurement && m_backgroundMeasurement);
  const bool isSub = m_spectrum->backgroundSubtract();
  if( m_backgroundSubItems[0]->isEnabled() != canSub )
    m_backgroundSubItems[0]->setDisabled( !canSub );
  if( m_backgroundSubItems[0]->isHidden() != isSub )
  {
    m_backgroundSubItems[0]->setHidden( isSub );
    m_backgroundSubItems[1]->setHidden( !isSub );
  }//if( m_backgroundSubItems[0]->isHidden() != isSub )
  
  if( m_hardBackgroundSub->isEnabled() != canSub )
    m_hardBackgroundSub->setDisabled( !canSub );
}//void displayBackgroundData()

