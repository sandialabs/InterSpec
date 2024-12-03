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

#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WPoint>
#include <Wt/WServer>
#include <Wt/Dbo/Dbo>
#include <Wt/WTextArea>
#include <Wt/WCheckBox>
#include <Wt/WIOService>
#include <Wt/WTabWidget>
#include <Wt/WPopupMenu>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WBorderLayout>
#include <Wt/WSelectionBox>
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
#endif

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/MakeDrf.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/FluxTool.h"
#include "InterSpec/PeakEdit.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/DecayWindow.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/OneOverR2Calc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DoseCalcWidget.h"
#include "InterSpec/ExportSpecFile.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "InterSpec/PeakInfoDisplay.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpecFileSummary.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/AddNewPeakDialog.h"
#include "InterSpec/ColorThemeWindow.h"
#include "InterSpec/GammaCountDialog.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/EnterAppUrlWindow.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/MultimediaDisplay.h"
#include "InterSpec/CompactFileManager.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/UnitsConverterTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/FeatureMarkerWidget.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/FileDragUploadResource.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"
#include "InterSpec/EnergyCalPreserveWindow.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/LicenseAndDisclaimersWindow.h"

#include "InterSpec/D3TimeChart.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

#if( IOS )
#include "target/ios/InterSpec/FileHandling.h"
#endif 

#if( INCLUDE_ANALYSIS_TEST_SUITE )
#include "InterSpec/SpectrumViewerTester.h"
#endif

#if( USE_GOOGLE_MAP )
#include "InterSpec/GoogleMap.h"
#elif( USE_LEAFLET_MAP )
#include "InterSpec/LeafletRadMap.h"
#endif

#if( USE_SEARCH_MODE_3D_CHART )
#include "InterSpec/SearchMode3DChart.h"
#endif

#if( USE_TERMINAL_WIDGET )
#include "InterSpec/TerminalWidget.h"
#endif

#if( USE_DETECTION_LIMIT_TOOL )
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/DetectionLimitSimple.h"
#endif

#if( USE_SPECRUM_FILE_QUERY_WIDGET )
#include "InterSpec/SpecFileQueryWidget.h"
#endif

#if( SpecUtils_ENABLE_D3_CHART )
#include "SpecUtils/D3SpectrumExport.h"
#endif

#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#endif

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP)
#include "js/AppHtmlMenu.js"
#endif

#if( BUILD_AS_WX_WIDGETS_APP )
#include "target/wxWidgets/InterSpecWxUtils.h"
#endif 


#if( USE_REMOTE_RID )
#include "InterSpec/RemoteRid.h"
#endif

#if( USE_REL_ACT_TOOL )
#include "InterSpec/RelActAutoGui.h"
#include "InterSpec/RelActManualGui.h"
#endif

#include "js/InterSpec.js"

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

using namespace Wt;
using namespace std;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

namespace
{
  std::mutex ns_staticDataDirectoryMutex;
  bool sm_haveSetStaticDataDirectory = false;
  std::string ns_staticDataDirectory = "data";

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  std::mutex ns_writableDataDirectoryMutex;
  std::string ns_writableDataDirectory = "";
#endif  //if( not a webapp )


  static const string PeakInfoTabTitleKey(       "app-tab-peak-manager" );
  static const string GammaLinesTabTitleKey(     "app-tab-ref-photopeak" );
  static const string CalibrationTabTitleKey(    "app-tab-energy-cal" );
  static const string NuclideSearchTabTitleKey(  "app-tab-nuc-search" );
  static const string FileTabTitleKey(           "app-tab-spec-files" );
#if( USE_TERMINAL_WIDGET )
  static const string TerminalTabTitleKey(       "app-tab-terminal" );
#endif
#if( USE_REL_ACT_TOOL )
  static const string RelActManualTitleKey(      "app-tab-isotopics" );
#endif

//#if( !BUILD_FOR_WEB_DEPLOYMENT )
//  const WTabWidget::LoadPolicy TabLoadPolicy = WTabWidget::LazyLoading;
//#else
  const WTabWidget::LoadPolicy TabLoadPolicy = WTabWidget::PreLoading;
//#endif

  void postSvlogHelper( const WString &msg, const int priority )
  {
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    if( app )
      app->svlog( msg, priority );
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
  
  
  template<class T>
  void del_ptr_set_null( T * &ptr )
  {
    if( ptr )
    {
      delete ptr;
      ptr = nullptr;
    }
  }//void del_ptr_set_null( T * &ptr )
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    /** A simple wrapper for holding the tool-tabs.
     When/if the user resizes the tool-tabs height, this function will call back and let the
     InterSpec class know, so we can restore this same hight when/if the tool tabs are
     toggled.
     */
    class ToolTabContentWrapper : public Wt::WContainerWidget
    {
      int m_height;
      InterSpec *m_interspec;
    
    public:
      ToolTabContentWrapper( InterSpec *interspec, WContainerWidget *parent = nullptr )
      : WContainerWidget( parent ),
      m_height( 0 ),
      m_interspec( interspec )
      {
        assert( m_interspec );
        setLayoutSizeAware( true );
      }
      
      virtual void layoutSizeChanged( int width, int height )
      {
        //cout << "ToolTabContentWrapper: w x h=" << width << " x " << height << endl;
        if( m_height == height )
          return;
        
        m_height = height;
        m_interspec->toolTabContentHeightChanged( height );
      }
    };//class ToolTabContentWrapper
#endif // InterSpec_PHONE_ROTATE_FOR_TABS
}//namespace


InterSpec::InterSpec( WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_user{},
    m_preferences( nullptr ),
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
    m_currentToolsTab( -1 ),
    m_toolsTabs( 0 ),
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    m_toolsTabsContentHeight( 0 ),
#endif
    m_energyCalTool( 0 ),
    m_energyCalWindow( 0 ),
    m_gammaCountDialog( 0 ),
    m_specFileQueryDialog( 0 ),
    m_shieldingSuggestion( 0 ),
    m_shieldingSourceFit( 0 ),
    m_shieldingSourceFitWindow( 0 ),
#if( USE_REL_ACT_TOOL )
    m_relActAutoGui( nullptr ),
    m_relActAutoWindow( nullptr ),
    m_relActAutoMenuItem( nullptr ),
    m_relActManualGui( nullptr ),
    m_relActManualWindow( nullptr ),
    m_relActManualMenuItem( nullptr ),
#endif
    m_materialDB( nullptr ),
    m_nuclideSearchWindow( 0 ),
    m_nuclideSearchContainer(0),
    m_nuclideSearch( 0 ),
    m_fileMenuPopup( 0 ),
    m_editMenuPopup( nullptr ),
    m_toolsMenuPopup( 0 ),
    m_helpMenuPopup( 0 ),
    m_displayOptionsPopupDiv( 0 ),
#if( USE_DB_TO_STORE_SPECTRA )
    m_saveState( 0 ),
    m_saveStateAs( 0 ),
    m_createTag( 0 ),
#endif
    m_languagesSubMenu( nullptr ),
    m_rightClickMenu( 0 ),
    m_rightClickEnergy( -DBL_MAX ),
    m_rightClickNuclideSuggestMenu( nullptr ),
    m_rightClickChangeContinuumMenu( nullptr ),
    m_rightClickChangeSkewMenu( nullptr ),
    m_exportSpecFileMenu{ nullptr },
    m_exportSpecFileWindow{ nullptr },
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
  m_featureMarkersWindow( nullptr ),
  m_featureMarkerMenuItem( nullptr ),
  m_multimedia( nullptr ),
#if( USE_REMOTE_RID )
  m_autoRemoteRidResultDialog( nullptr ),
#endif
  m_gammaXsToolWindow( nullptr ),
  m_doseCalcWindow( nullptr ),
  m_1overR2Calc( nullptr ),
  m_unitsConverter( nullptr ),
  m_fluxTool( nullptr ),
  m_makeDrfTool( nullptr ),
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  m_mapMenuItem( nullptr ),
#if( USE_LEAFLET_MAP )
  m_leafletWarning( nullptr ),
  m_leafletWindow( nullptr ),
#endif
#endif
  m_enterUri( nullptr ),
#if( USE_SEARCH_MODE_3D_CHART )
  m_searchMode3DChart( 0 ),
#endif
  m_showRiidResults( nullptr ),
  m_showMultimedia( nullptr ),
#if( USE_TERMINAL_WIDGET )
  m_terminalMenuItem( 0 ),
  m_terminal( 0 ),
  m_terminalWindow( 0 ),
#endif
#if( USE_REMOTE_RID )
  m_remoteRidMenuItem( nullptr ),
  m_remoteRid( nullptr ),
  m_remoteRidWindow( nullptr ),
#endif
#if( USE_DETECTION_LIMIT_TOOL )
  m_simpleMdaWindow( nullptr ),
  m_detectionLimitWindow( nullptr ),
#endif
  m_clientDeviceType( 0x0 ),
  m_referencePhotopeakLines( 0 ),
  m_referencePhotopeakLinesWindow( 0 ),
  m_helpWindow( nullptr ),
  m_licenseWindow( nullptr ),
  m_useInfoWindow( 0 ),
  m_decayInfoWindow( nullptr ),
  m_addFwhmTool( nullptr ),
  m_preserveCalibWindow( 0 ),
#if( USE_SEARCH_MODE_3D_CHART )
  m_3dViewWindow( nullptr ),
#endif
  m_riidDisplay( nullptr ),
  m_drfSelectWindow( nullptr ),
  m_undo( nullptr ),
  m_renderedWidth( 0 ),
  m_renderedHeight( 0 ),
  m_colorPeaksBasedOnReferenceLines( true ),
  m_findingHintPeaks( false ),
  m_hintQueue{},
  m_infoNotificationsMade{}
{
  //Initialization of the app (this function) takes about 11ms on my 2.6 GHz
  //  Intel Core i7, as of (20150316).
  
  addStyleClass( "InterSpec" );
  
  //Setting setLayoutSizeAware doesnt seem to add any appreciable network
  //  overhead, at least according to the chrome inspection panel when running
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
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
  m_notificationDiv->addStyleClass( "belowMenu" );
#endif
  
  app->domRoot()->addWidget( m_notificationDiv );
  
  app->hotkeySignal().connect( boost::bind( &InterSpec::hotKeyPressed, this, boost::placeholders::_1 ) );
  
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
      // We'll set the user name to be a combination of the initial time, and a hash of the
      //  Wt session ID to make sure its unique-enough, without bringing in anything heavyweight
      auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      username = "user-" + SpecUtils::to_iso_string(now)
                 + "-" + std::to_string( std::hash<std::string>{}(wApp->sessionId()) );
      wApp->setCookie( "SpectrumViewerUUID", username, 3600*24*365 );
    }//try / catch
  }//if( no username )

  
  if( app->isPhone() )
    username += "_phone";
  else if( app->isTablet() )
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
      
    }else
    {
      InterSpecUser::DeviceType type = InterSpecUser::Desktop;
      if( app->isPhone() )
        type = InterSpecUser::PhoneDevice;
      else if( app->isTablet() )
        type = InterSpecUser::TabletDevice;
    
      InterSpecUser *newuser = new InterSpecUser( username, type );
      m_user = m_sql->session()->add( newuser );
    }//if( m_user ) / else
  
    m_preferences = new UserPreferences( this );
    
    m_user.modify()->startingNewSession();
    
    transaction.commit();
  }//end interacting with DB
  
  detectClientDeviceType();

  const string langPref = UserPreferences::preferenceValue<string>("Language", this);
  if( !langPref.empty() )
    wApp->setLocale( WLocale(langPref) );
      
  app->useMessageResourceBundle( "InterSpec" );
    
  // Now that we have m_sql and m_user setup, we can create the undo/redo manager, if we
  //  are using the desktop interface.  We will create this manager before any our widgets
  //
  // We wont enable undo/redo when we are using mobile-menu (i.e., phones, or tablets that dont
  //  have desktop interface enabled).
  std::unique_ptr<UndoRedoManager::BlockUndoRedoInserts> undo_blocker;
  if( (UndoRedoManager::maxUndoRedoSteps() >= 0)
      /* && !app->isPhone()
      && (!app->isTablet() || UserPreferences::preferenceValue<bool>("TabletUseDesktopMenus", this)) */ )
  {
    m_undo = new UndoRedoManager( this );
    undo_blocker = std::unique_ptr<UndoRedoManager::BlockUndoRedoInserts>();
  }//if( desktop interface )
    
  m_peakModel = new PeakModel( this );
  m_spectrum   = new D3SpectrumDisplayDiv();
  m_timeSeries = new D3TimeChart();
  
  if( app->isMobile() )
  {
    //TODO: layoutSizeChanged(...) will trigger the compact axis anyway, but need to check if doing it here saves a roundtrip
    m_spectrum->setCompactAxis( true );
    m_timeSeries->setCompactAxis( true );
    
    LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsDoOrientationChange);
    
    const char *js = INLINE_JAVASCRIPT(
      window.addEventListener("orientationchange", Wt.WT.DoOrientationChange );
      setTimeout( Wt.WT.DoOrientationChange, 0 );
      setTimeout( Wt.WT.DoOrientationChange, 250 );  //JIC - doesnt look necessary
      setTimeout( Wt.WT.DoOrientationChange, 1000 );
    );
    
    doJavaScript( js );
  }//if( isMobile() )
  
  m_spectrum->setPeakModel( m_peakModel );
  m_spectrum->existingRoiEdgeDragUpdate().connect( m_spectrum, &D3SpectrumDisplayDiv::performExistingRoiEdgeDragWork );
  m_spectrum->dragCreateRoiUpdate().connect( m_spectrum, &D3SpectrumDisplayDiv::performDragCreateRoiWork );
  
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
  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
  const bool isAppTitlebar = InterSpecApp::isPrimaryWindowInstance();
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
      
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
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
      
#if( BUILD_AS_WX_WIDGETS_APP )
      if(InterSpecApp::isPrimaryWindowInstance())
      {
        //The Edge WebView2 doesnt respect the "-webkit-app-region: drag" property, so we will
        //  do a workaround to have wxWidgets capture the mouse when titlebar is clicked on
        const char* js = "function(el,evt){"
          // Make sure its left mouse button, and for `m_menuDiv`, make sure element clicked on was not a descendant (like a menu button)
          "if( (evt.button!==0) || (el.id !== evt.target.id)) return;"
          "window.wx.postMessage('MouseDownInTitleBar');"
          "evt.stopPropagation();"
          "evt.preventDefault();"
          "}";

        //titleStretcher doesnt work for this since it is all margin, and no actual width
        dragRegion->mouseWentDown().connect(js); //Doesnt seem to be useful; havent checked out why
        menuTitle->mouseWentDown().connect(js);
        appIcon->mouseWentDown().connect(js);
        m_menuDiv->mouseWentDown().connect(js);

        // Should we also send wxWidgets a mouse up signal, e.g., ` m_menuDiv->mouseWentUp().connect(...)`?
      }//if( primary wxWidgets instance )
#endif
      
      
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP)
      if (InterSpecApp::isPrimaryWindowInstance())
      {
#if( BUILD_AS_ELECTRON_APP )
        auto toggleMaximize = []() {
            ElectronUtils::send_nodejs_message("ToggleMaximizeWindow", "");
        };//doMaximize

        maximizeIcon->clicked().connect(std::bind(toggleMaximize));
        appIcon->doubleClicked().connect(std::bind(toggleMaximize));
        dragRegion->doubleClicked().connect(std::bind(toggleMaximize));
        menuTitle->doubleClicked().connect(std::bind(toggleMaximize));
#elif( BUILD_AS_WX_WIDGETS_APP )
        const string js = "function(){window.wx.postMessage('ToggleMaximizeWindow');}";
        maximizeIcon->clicked().connect(js);
        // None of the below double-click events seem to work
        appIcon->doubleClicked().connect(js);
        dragRegion->doubleClicked().connect(js);
        menuTitle->doubleClicked().connect(js);
        m_menuDiv->doubleClicked().connect(js);
#endif
      }
      else
      {
#endif
        const char* toggleMaxJs =
          "function(){let elem = document.querySelector(\".Wt-domRoot\");"
          "if (!document.fullscreenElement) {"
          "  elem.requestFullscreen().catch(err => {"
          "    console.log( 'Error attempting to enable full-screen mode' );"
          "  });"
          "} else {"
          "  document.exitFullscreen();"
          "}"
          "}";

        maximizeIcon->clicked().connect(toggleMaxJs);
        appIcon->doubleClicked().connect(toggleMaxJs);
        dragRegion->doubleClicked().connect(toggleMaxJs);
        menuTitle->doubleClicked().connect(toggleMaxJs);
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP)
      }
#endif

#if( BUILD_AS_ELECTRON_APP )
      minimizeIcon->clicked().connect( std::bind([](){
        ElectronUtils::send_nodejs_message( "MinimizeWindow", "" );
      }) );
      
      closeIcon->clicked().connect( std::bind([](){
        ElectronUtils::send_nodejs_message( "CloseWindow", "" );
      }) );
#endif //BUILD_AS_ELECTRON_APP

#if( BUILD_AS_WX_WIDGETS_APP )
      minimizeIcon->clicked().connect( "function(){ window.wx.postMessage('MinimizeWindow'); }" );
      closeIcon->clicked().connect( "function(){ window.wx.postMessage('CloseWindow'); }" );
#endif //BUILD_AS_WX_WIDGETS_APP


#else //#if( BUILD_AS_ELECTRON_APP - for dev purposes )
      assert( 0 );
#endif //#if( BUILD_AS_ELECTRON_APP - for dev pupropses )
    }//if( !isAppTitlebar ) / else
  }//if( isMobile() ) / else

  addFileMenu( menuWidget, isAppTitlebar );
  addEditMenu( menuWidget );
  addViewMenu( menuWidget );
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
    m_toolsTabs->addStyleClass( "ToolsTabs" );
    
    CompactFileManager *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::LeftToRight );
    m_toolsTabs->addTab( compact, WString::tr(FileTabTitleKey), TabLoadPolicy );
    
    m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                   boost::placeholders::_1,
                                                   boost::placeholders::_2,
                                                   boost::placeholders::_3 ) );
    
    m_toolsTabs->addTab( m_peakInfoDisplay, WString::tr(PeakInfoTabTitleKey), TabLoadPolicy );
    
    m_referencePhotopeakLines = new ReferencePhotopeakDisplay( m_spectrum,
                                                              m_materialDB.get(),
                                                              m_shieldingSuggestion,
                                                              this );
    setReferenceLineColors( nullptr );
    
    //PreLoading is necessary on the m_referencePhotopeakLines widget, so that the
    //  "Isotope Search" widget will work properly when a nuclide is clicked
    //  on to display its photopeaks
    //XXX In Wt 3.3.4 at least, the contents of m_referencePhotopeakLines
    //  are not actually loaded to the client until the tab is clicked, and I
    //  cant seem to get this to actually happen.
    //WMenuItem *refPhotoTab =
    m_toolsTabs->addTab( m_referencePhotopeakLines, WString::tr(GammaLinesTabTitleKey), TabLoadPolicy );
    
    m_toolsTabs->currentChanged().connect( this, &InterSpec::handleToolTabChanged );
    
    m_energyCalTool->setWideLayout();
    m_toolsTabs->addTab( m_energyCalTool, WString::tr(CalibrationTabTitleKey), TabLoadPolicy );
    
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
    m_toolsTabs->addTab( m_nuclideSearchContainer, WString::tr(NuclideSearchTabTitleKey), TabLoadPolicy );
    //    const char *tooltip = "Search for nuclides with constraints on energy, "
    //                          "branching ratio, and half life.";
    //    HelpSystem::attachToolTipOn( nuclideTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
#if( USE_TERMINAL_WIDGET || USE_REL_ACT_TOOL )
    // Handle when the user closes the tab for the Math/Command terminal and the Manual Relative
    //  Activity tool
    m_toolsTabs->tabClosed().connect( boost::bind( &InterSpec::handleToolTabClosed, this, boost::placeholders::_1 ) );
#endif
    
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
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-peak-editor") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::peakEditFromRightClick );
        break;
      case kRefitPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-refit-peak") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::refitPeakFromRightClick );
        break;
      case kRefitPeakWithDrfFwhm:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-use-drf-fwhm") );
        m_rightClickMenutItems[i]->setToolTip( WString::tr("rclick-mi-tt-use-drf-fwhm") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::refitPeakWithDrfFwhm );
        break;
        
      case kSetMeanToRefPhotopeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-fix-mean") );
        m_rightClickMenutItems[i]->setToolTip( WString::tr("rclick-mi-tt-fix-mean") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::setMeanToRefPhotopeak );
        break;
        
      case kRefitROI:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-refit-roi") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::refitPeakFromRightClick );
        break;
      case kDeletePeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-delete-peak") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::deletePeakFromRightClick );
        break;
      case kChangeNuclide:
        m_rightClickNuclideSuggestMenu = m_rightClickMenu->addPopupMenuItem( WString::tr("rclick-mi-change-nuc") );
        m_rightClickMenutItems[i] = m_rightClickNuclideSuggestMenu->parentItem();
        break;
      case kChangeContinuum:
      {
        m_rightClickChangeContinuumMenu = m_rightClickMenu->addPopupMenuItem( WString::tr("rclick-mi-change-cont") );
        m_rightClickMenutItems[i] = m_rightClickChangeContinuumMenu->parentItem();
        
        for( auto type = PeakContinuum::OffsetType(0);
            type <= PeakContinuum::External; type = PeakContinuum::OffsetType(type+1) )
        {
          WMenuItem *item = m_rightClickChangeContinuumMenu->addItem( WString::tr( PeakContinuum::offset_type_label_tr(type) ) );
          item->triggered().connect( boost::bind( &InterSpec::handleChangeContinuumTypeFromRightClick, this, type ) );
        }//for( loop over PeakContinuum::OffsetTypes )
        break;
      }//case kChangeContinuum:
        
      case kChangeSkew:
      {
        m_rightClickChangeSkewMenu = m_rightClickMenu->addPopupMenuItem( WString::tr("rclick-mi-skew-type") );
        m_rightClickMenutItems[i] = m_rightClickChangeSkewMenu->parentItem();
        for( auto type = PeakDef::SkewType(0);
            type <= PeakDef::SkewType::DoubleSidedCrystalBall; type = PeakDef::SkewType(type+1) )
        {
          WMenuItem *item = m_rightClickChangeSkewMenu->addItem( PeakDef::to_label(type) );
          item->triggered().connect( boost::bind( &InterSpec::handleChangeSkewTypeFromRightClick,
                                                 this, static_cast<int>(type) ) );
        }//for( loop over PeakContinuum::OffsetTypes )
        break;
      }
        
      case kAddPeakToRoi:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-add-peak") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::addPeakFromRightClick );
        break;
      case kShareContinuumWithLeftPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-combine-cont-left") );
        m_rightClickMenutItems[i]->triggered().connect( boost::bind( &InterSpec::shareContinuumWithNeighboringPeak, this, true ) );
        break;
      case kShareContinuumWithRightPeak:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-combine-cont-right") );
        m_rightClickMenutItems[i]->triggered().connect( boost::bind( &InterSpec::shareContinuumWithNeighboringPeak, this, false ) );
        break;
        
      case kMakeOwnContinuum:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-own-cont") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::makePeakFromRightClickHaveOwnContinuum );
        break;
        
#if( USE_DETECTION_LIMIT_TOOL )
      case kFitNewPeakNotInRoi:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-fit-new-peak") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::fitNewPeakNotInRoiFromRightClick );
      break;
      case kAddPeakNotInRoi:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-add-new-peak") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::startAddPeakFromRightClick );
      break;
      case kSearchEnergy:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-mi-search-energy") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::searchOnEnergyFromRightClick );
      break;
      case kSimpleMda:
        m_rightClickMenutItems[i] = m_rightClickMenu->addMenuItem( WString::tr("rclick-simple-mda") );
        m_rightClickMenutItems[i]->triggered().connect( this, &InterSpec::startSimpleMdaFromRightClick );
      break;
#endif
        
      case kNumRightClickItems:
        break;
    }//switch( i )
  }//for( loop over right click menu items )

  // For touch devices, give an obvious way to close the right-click menu since they cant simply
  //  leave with the mouse (although they could just tap somewhere).
  if( app && app->isMobile() && !m_rightClickMenu->isMobile() )
    m_rightClickMenu->addMenuItem( WString::tr("Cancel") );
    
  m_spectrum->rightClicked().connect( boost::bind( &InterSpec::handleRightClick, this,
                                                  boost::placeholders::_1, boost::placeholders::_2,
                                                  boost::placeholders::_3, boost::placeholders::_4 ) );
  m_spectrum->chartClicked().connect( boost::bind( &InterSpec::handleLeftClick, this,
                                                  boost::placeholders::_1, boost::placeholders::_2,
                                                  boost::placeholders::_3, boost::placeholders::_4 ) );
  
  m_spectrum->shiftKeyDragged().connect( boost::bind( &InterSpec::excludePeaksFromRange, this,
                                                     boost::placeholders::_1,
                                                     boost::placeholders::_2 ) );
  m_spectrum->doubleLeftClick().connect( boost::bind( &InterSpec::searchForSinglePeak, this,
                                                     boost::placeholders::_1 ) );
  m_spectrum->xRangeChanged().connect( boost::bind( &InterSpec::handleSpectrumChartXRangeChange, this,
                                                     boost::placeholders::_1,
                                                   boost::placeholders::_2,
                                                   boost::placeholders::_3,
                                                   boost::placeholders::_4,
                                                   boost::placeholders::_5) );
      
  m_spectrum->yAxisScaled().connect(
          boost::bind( &InterSpec::handleDisplayScaleFactorChangeFromSpectrum, this,
                      boost::placeholders::_1, boost::placeholders::_2, boost::placeholders::_3 ) );
    
  m_spectrum->legendDisabled().connect( std::bind([this](){
    if( m_undo && !m_undo->isInUndoOrRedo() )
      m_undo->addUndoRedoStep( [this](){ m_spectrum->enableLegend(); },
                               [this](){ m_spectrum->disableLegend(); },
                              "Hide spectrum legend." );
    } ) );
    
  m_spectrum->legendEnabled().connect( std::bind([this](){
    if( m_undo && !m_undo->isInUndoOrRedo() )
      m_undo->addUndoRedoStep( [this](){ m_spectrum->disableLegend(); },
                               [this](){ m_spectrum->enableLegend(); },
                              "Show spectrum legend." );
    } ) );
    
    
  m_timeSeries->setHidden( true );
  m_chartResizer->setHidden( m_timeSeries->isHidden() );
  
#if( USE_OSX_NATIVE_MENU )
  if( InterSpecApp::isPrimaryWindowInstance() )
    m_menuDiv->hide();
#endif
  
  if( !isPhone() )
  {
    setToolTabsVisible( true );
  }else
  {
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    const WEnvironment &env = wApp->environment();
    const bool isVertical = (env.screenWidth() < env.screenHeight());
    setToolTabsVisible( !isVertical );
#else
    //Disallow showing tool tabs when on a phone
    m_toolTabsVisibleItems[1]->setHidden(true);
    
    //Add Ref photopeaks, etc. to tool menu
    addToolsTabToMenuItems();
#endif //#endif InterSpec_PHONE_ROTATE_FOR_TABS
  }//If( start with tool tabs showing ) / else
 
  initDragNDrop();
  
#if( USE_DB_TO_STORE_SPECTRA )
  updateSaveWorkspaceMenu();
#endif
  
#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !BUILD_AS_ELECTRON_APP )
  initOsColorThemeChangeDetect();
#endif
  
  applyColorTheme( nullptr );
  
#if( APP_MENU_STATELESS_FIX )
  // Make sure the menus get pre-loaded
  PopupDivMenu::pre_render(m_fileMenuPopup);
  PopupDivMenu::pre_render(m_editMenuPopup);
  PopupDivMenu::pre_render(m_toolsMenuPopup);
  PopupDivMenu::pre_render(m_helpMenuPopup);
  PopupDivMenu::pre_render(m_displayOptionsPopupDiv);
#endif
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
  std::lock_guard<std::mutex> lock( ns_staticDataDirectoryMutex );
  
  if( !SpecUtils::is_directory(dir) )
    throw runtime_error( "InterSpec::setStaticDataDirectory(): " + dir + " is not a directory." );
  
  ns_staticDataDirectory = dir;
  sm_haveSetStaticDataDirectory = true;

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
  
  IsotopeId::setDataDirectory( dir );
}//void setStaticDataDirectory( const std::string &dir )


std::string InterSpec::staticDataDirectory()
{
  std::lock_guard<std::mutex> lock( ns_staticDataDirectoryMutex );
  return ns_staticDataDirectory;
}

bool InterSpec::haveSetStaticDataDirectory()
{
  std::lock_guard<std::mutex> lock( ns_staticDataDirectoryMutex );
  return sm_haveSetStaticDataDirectory;
}

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
void InterSpec::setWritableDataDirectory( const std::string &dir )
{
  std::lock_guard<std::mutex> lock( ns_writableDataDirectoryMutex );
  
  if( !dir.empty() && !SpecUtils::is_directory(dir) )
    throw runtime_error( "InterSpec::setWritableDataDirectory(): " + dir + " is not a directory." );
  
  //Set the serial to module database file, if the user has one in thier
  //  application data directory.
  //Currently this is being done in the target specific code; it should all be
  //  moved here probably.
  //const vector<string> serial_db = SpecUtils::ls_files_in_directory( dir, "serial_to_model.csv" );
  //if( !serial_db.empty() )
  //  SerialToDetectorModel::set_detector_model_input_csv( serial_db[0] );
  
  ns_writableDataDirectory = dir;
}//setWritableDataDirectory( const std::string &dir )

std::string InterSpec::writableDataDirectory()
{
  std::lock_guard<std::mutex> lock( ns_writableDataDirectoryMutex );
  
  if( ns_writableDataDirectory.empty() )
    throw runtime_error( "writableDataDirectory hasnt been set." );
  
  return ns_writableDataDirectory;
}//string writableDataDirectory()
#endif  //if( not a webapp )



InterSpec::~InterSpec() noexcept(true)
{
  // The DOM root may not actually being deleted, most likely if the "Clear Session..." option
  //  was invoked.  WPopupWidget descendants, except for AuxWindow, SimpleDialog, and any
  //  WPopupWidget-derived widget explicitly given a parent, will not be deleted, so we'll do
  //  some manual cleanup here (as of 20220917 when AuxWindow and SimpleDialog where explicitly
  //  parented by the current InterSpec instance, we are doing much more cleanup than necessary).

  Wt::log("info") << "Destructing InterSpec from session '" << (wApp ? wApp->sessionId() : string("")) << "'";

  // Get rid of undo/redo, so we dont insert anything into them
  del_ptr_set_null( m_undo );
  del_ptr_set_null( m_licenseWindow );
  
  try
  {
    closeShieldingSourceFit();
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
    
    del_ptr_set_null( m_peakInfoDisplay );
  }//if( m_peakInfoDisplay )

  del_ptr_set_null( m_peakInfoWindow );
  
  if( m_energyCalTool )
  {
    if( m_toolsTabs && m_toolsTabs->indexOf(m_energyCalTool)>=0 )
      m_toolsTabs->removeTab( m_energyCalTool );
    
    del_ptr_set_null( m_energyCalTool );
  }//if( m_energyCalTool )
  
  del_ptr_set_null( m_energyCalWindow );
  del_ptr_set_null( m_shieldingSourceFitWindow );
  del_ptr_set_null( m_nuclideSearch );
  del_ptr_set_null( m_nuclideSearchWindow );
  
  if( m_referencePhotopeakLines )
  {
    if( m_toolsTabs && m_toolsTabs->indexOf(m_referencePhotopeakLines)>=0 )
      m_toolsTabs->removeTab( m_referencePhotopeakLines );
    
    m_referencePhotopeakLines->clearAllLines();
    del_ptr_set_null( m_referencePhotopeakLines );
  }//if( m_referencePhotopeakLines )
  
  del_ptr_set_null( m_referencePhotopeakLinesWindow );
  
  handleWarningsWindowClose();
  del_ptr_set_null( m_warnings ); //WarningWidget isnt necessarily parented, so we do have to manually delete it
  
  deletePeakEdit();
  deleteGammaCountDialog();
  
  // The following may be parented by app->domRoot()
  del_ptr_set_null( m_mobileMenuButton );
  del_ptr_set_null( m_mobileBackButton );
  del_ptr_set_null( m_mobileForwardButton );
  del_ptr_set_null( m_notificationDiv );
  del_ptr_set_null( m_fileManager );
  del_ptr_set_null( m_displayOptionsPopupDiv );
  del_ptr_set_null( m_fileMenuPopup );
  del_ptr_set_null( m_editMenuPopup );
  del_ptr_set_null( m_toolsMenuPopup );
  del_ptr_set_null( m_helpMenuPopup );
  del_ptr_set_null( m_menuDiv );
  
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


D3SpectrumExport::D3SpectrumChartOptions InterSpec::getD3SpectrumOptions() const
{
  double xMin, xMax, yMin, yMax;
  displayedSpectrumRange(xMin, xMax, yMin, yMax);
  
  map<string,string> referc_line_json;
  
  if( m_referencePhotopeakLines )
    referc_line_json = m_referencePhotopeakLines->jsonReferenceLinesMap();

  const string title = "Interactive Spectrum Development";
  const string xAxisTitle = WString::tr("Energy").toUTF8();
  const string yAxisTitle = WString::tr("Counts").toUTF8();
  const string dataTitle = (m_spectrum->data() ? m_spectrum->data()->title() :
                            m_spectrum->background() ? m_spectrum->background()->title() :
                            m_spectrum->secondData() ? m_spectrum->secondData()->title() :
                            WString::tr("Foreground").toUTF8());
  const bool useLogYAxis = m_spectrum->yAxisIsLog();
  const bool showVerticalGridLines = m_spectrum->verticalLinesShowing();
  const bool showHorizontalGridLines = m_spectrum->horizontalLinesShowing();
  const bool legendEnabled = m_spectrum->legendIsEnabled();
  const bool compactXAxis = m_spectrum->isAxisCompacted();
  const bool showPeakUserLabels = m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakUserLabel );
  const bool showPeakEnergyLabels = m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakEnergyLabel );
  const bool showPeakNuclideLabels = m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakNuclideLabel );
  const bool showPeakNuclideEnergyLabels = m_spectrum->showingPeakLabel( SpectrumChart::kShowPeakNuclideEnergies );
  const bool showEscapePeakMarker = (m_featureMarkersWindow && m_featureMarkersShown[static_cast<int>(FeatureMarkerType::EscapePeakMarker)]);
  const bool showComptonPeakMarker = (m_featureMarkersWindow && m_featureMarkersShown[static_cast<int>(FeatureMarkerType::ComptonPeakMarker)]);
  const bool showComptonEdgeMarker = (m_featureMarkersWindow && m_featureMarkersShown[static_cast<int>(FeatureMarkerType::ComptonEdgeMarker)]);
  const bool showSumPeakMarker = (m_featureMarkersWindow && m_featureMarkersShown[static_cast<int>(FeatureMarkerType::SumPeakMarker)]);
  const bool backgroundSubtract = m_spectrum->backgroundSubtract();
  
  
  D3SpectrumExport::D3SpectrumChartOptions options( title, xAxisTitle, yAxisTitle, dataTitle,
      useLogYAxis, showVerticalGridLines, showHorizontalGridLines, legendEnabled, compactXAxis,
      showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels, showPeakNuclideEnergyLabels,
      showEscapePeakMarker, showComptonPeakMarker, showComptonEdgeMarker, showSumPeakMarker,
      backgroundSubtract, xMin, xMax,
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
  
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  if( isPhone() )
  {
    if( (w <= 20) || (h <= 20) )
    {
      w = wApp->environment().screenWidth();
      h = wApp->environment().screenHeight();
    }
    
    if( (w > 20) && (h > 20) )
    {
      const bool isVertical = (h > w);
      
      // If we are changing orientation - close all the open windows
      //  TODO: close are restore all the windows, instead of just closing them
#if( USE_CSS_FLEX_LAYOUT )
      if( m_toolsTabs->isVisible() != isVertical )
#else
      if( static_cast<bool>(m_toolsTabs) != isVertical )
#endif
      {
        //m_energyCalWindow and m_nuclideSearchWindow will get closed by `setToolTabsVisible(...)`
        if( m_gammaCountDialog )
          deleteGammaCountDialog();
        
        if( m_shieldingSourceFitWindow )
          closeShieldingSourceFit();
        
#if( USE_REL_ACT_TOOL )
        if( m_relActAutoWindow )
          m_relActAutoWindow->hide(); //rejects if, calling `handleRelActAutoClose()`
        if( m_relActManualWindow )
          m_relActManualWindow->hide(); //rejects if, calling `handleRelActManualClose()`
#endif
        
        if( m_multimedia )
          programmaticallyCloseMultimediaWindow();
        
#if( USE_REMOTE_RID )
        if( m_autoRemoteRidResultDialog )
          programaticallyCloseAutoRemoteRidResultDialog();
#endif
        if( m_gammaXsToolWindow )
          m_gammaXsToolWindow->hide(); //calls deleteGammaXsTool
        if( m_doseCalcWindow )
          deleteDoseCalcTool();
        if( m_1overR2Calc )
          deleteOneOverR2Calc();
        if( m_unitsConverter )
          deleteUnitsConverterTool();
        if( m_fluxTool )
          deleteFluxTool();
        if( m_makeDrfTool )
          handleCloseMakeDrfWindow( m_makeDrfTool );
#if( USE_LEAFLET_MAP )
        if( m_leafletWindow )
          programmaticallyCloseLeafletMap();
#endif
        if( m_enterUri )
          m_enterUri->accept();
        assert( !m_enterUri );
        if( m_terminalWindow )
          m_terminalWindow->hide();
#if( USE_REMOTE_RID )
        if( m_remoteRidWindow )
          deleteRemoteRidWindow();
#endif
#if( USE_DETECTION_LIMIT_TOOL )
        if( m_simpleMdaWindow )
          handleSimpleMdaWindowClose();
        if( m_detectionLimitWindow )
          handleDetectionLimitWindowClose();
#endif
        if( m_helpWindow )
          closeHelpWindow();
        if( m_licenseWindow )
          deleteLicenseAndDisclaimersWindow();
        if( m_useInfoWindow )
          deleteWelcomeDialog( true );
        if( m_decayInfoWindow )
          deleteDecayInfoWindow();
        if( m_addFwhmTool )
          deleteFwhmFromForegroundWindow();
        if( m_preserveCalibWindow )
          deleteEnergyCalPreserveWindow();
#if( USE_SEARCH_MODE_3D_CHART )
        if( m_3dViewWindow )
          handle3DSearchModeChartClose( m_3dViewWindow );
#endif
        if( m_riidDisplay )
          programmaticallyCloseRiidResults();
        if( m_drfSelectWindow )
          closeDrfSelectWindow();
      }//if( we are changing orientation )
      
      setToolTabsVisible( isVertical );
    }//if( (w > 20) && (h > 20) )
  }//if( isPhone() )
#endif
  
  const bool comactX = (h <= 420); //Apple iPhone 6+, 6s+, 7+, 8+
  if( (h > 20) && (comactX != m_spectrum->isAxisCompacted()) )
  {
    //Only set axis compact if it isnt already; dont ever set non-compact
    //  ToDo: For case screen is made small, so axis goes compact, then user
    //        makes screen large, currently wont un-compactify axis, even if
    //        user wants this; should evaluate if its worth checking the user
    //        preference about this.
    //const bool makeCompact = UserPreferences::preferenceValue<bool>( "CompactXAxis", this );
    
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

bool InterSpec::isAndroid() const
{
  return (m_clientDeviceType & ClientDeviceType::AndroidClient);
}

void InterSpec::detectClientDeviceType()
{
  m_clientDeviceType= 0x0;

  InterSpecApp *app= dynamic_cast<InterSpecApp *>( wApp );
  if( !app )
    return;

  const bool phone  = app->isPhone();
  bool tablet = app->isTablet();
  bool mobile = app->isMobile();
  if( mobile && !phone && tablet )
    tablet = mobile = !UserPreferences::preferenceValue<bool>("TabletUseDesktopMenus", this);

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

      case DedicatedAppClient:
#if( !BUILD_FOR_WEB_DEPLOYMENT )
        m_clientDeviceType |= type;
#endif
        break;

      case AndroidClient:
        if( app->isAndroid() )
          m_clientDeviceType |= type;
        break;
        
      case NumClientDeviceType:
        break;
    }// switch( type )
  }// for( loop over ClientDeviceType enums )
}// void detectClientDeviceType()


void InterSpec::changeLocale( std::string languageCode )
{
  const string prevLangPref = UserPreferences::preferenceValue<string>("Language", this);
  //if( prevLangPref == local )
  //  return;
  
  const set<string> &languages = InterSpecApp::languagesAvailable();
  if( !languages.count(languageCode) && !languageCode.empty() )
  {
    passMessage( WString::tr("err-invalid-language").arg(languageCode),
                WarningWidget::WarningMsgHigh );
    return;
  }//if( !languages.count(languageCode) )
  
  UserPreferences::setPreferenceValue( "Language", languageCode, this );
  
  // Add the style class to the item selected
  if( m_languagesSubMenu )
  {
    for( WMenuItem *item : m_languagesSubMenu->items() )
    {
      item->removeStyleClass( "CurrentLanguage" );
      const string label = item->text().toUTF8();
      if( (label == languageCode)
         || (languageCode.empty() && (label == "Default")) )
      {
        item->addStyleClass( "CurrentLanguage" );
      }
    }//for( WMenuitem *item : m_languagesSubMenu->items() )
  }//if( m_languagesSubMenu )
  
  
  if( languageCode.empty() )
  {
    wApp->setLocale( WLocale() ); // the default locale is chosen.
  }else
  {
    wApp->setLocale( WLocale( languageCode ) );
  }
}//void changeLocale( std::string locale );


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


void InterSpec::hotKeyPressed( const unsigned int value )
{
  if( m_toolsTabs )
  {
    //string expectedTxt;
    switch( value )
    {
      case '1': case '2': case '3': case '4': case '5': case '6': case '7':
      {
        const int tabIndex = value - '1';
        
        if( tabIndex < m_toolsTabs->count() )
        {
          m_toolsTabs->setCurrentIndex( tabIndex );
          handleToolTabChanged( tabIndex );
        }
        
        break;
      }// case '1' through '7'
        
        
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
        
      case 'e': case 'E':
        createExportSpectrumFileDialog();
        break;
        
    // Temporarily add shortcut for showing FAQ window during development
      case 'f': case 'F':
      {
        showWelcomeDialog( true );
        if( m_useInfoWindow )
          m_useInfoWindow->showFaqTab();
        break;
      }
        
        
      case 37: case 38: case 39: case 40:
        arrowKeyPressed( value );
        break;
        
      case 'z':
        if( m_undo )
          m_undo->executeUndo();
        break;
        
      case 'Z':
        if( m_undo )
          m_undo->executeRedo();
        break;
    }//switch( value )
  
    
    /*
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
     */
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
#if( USE_OSX_NATIVE_MENU || IOS || ANDROID )
  return;
#endif
  
  if( m_mobileMenuButton )
    return;
  
  const vector<PopupDivMenu *> menus{
    m_fileMenuPopup, m_editMenuPopup, m_displayOptionsPopupDiv, m_toolsMenuPopup, m_helpMenuPopup
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
    size_t nextActiveIndex = activeIndex;
    do
    {
      if( leftArrow )
        nextActiveIndex = (nextActiveIndex == 0) ? (menus.size() - 1) : (nextActiveIndex - 1);
      else
        nextActiveIndex = ((nextActiveIndex + 1) % menus.size());
    }while( !menus[nextActiveIndex] );
    
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
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  PeakSearchGuiUtils::refit_peaks_from_right_click( this, m_rightClickEnergy );
}//void refitPeakFromRightClick()


void InterSpec::refitPeakWithDrfFwhm()
{
  PeakSearchGuiUtils::refit_peaks_with_drf_fwhm( this, m_rightClickEnergy );
}//void InterSpec::refitPeakWithDrfFwhm()


void InterSpec::setMeanToRefPhotopeak()
{
  PeakSearchGuiUtils::refit_peak_with_photopeak_mean( this, m_rightClickEnergy );
}//void setMeanToRefPhotopeak()


void InterSpec::addPeakFromRightClick()
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  std::shared_ptr<const SpecUtils::Measurement> dataH = m_spectrum->data();
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak
      || m_rightClickEnergy < peak->lowerX()
      || m_rightClickEnergy > peak->upperX()
      || !m_dataMeasurement
      || !dataH )
  {
    passMessage( WString::tr("err-no-roi-to-add-peak"), WarningWidget::WarningMsgInfo );
    return;
  }//if( !peak )
  
  //get all the peaks that belong to this ROI
  typedef std::shared_ptr<const PeakContinuum> ContPtr;
  typedef map<ContPtr, std::vector<PeakDef > > ContToPeakMap;
  
  
  const auto origContinumm = peak->continuum();
  const std::shared_ptr<const deque<PeakModel::PeakShrdPtr>> allOrigPeaksDeque = m_peakModel->peaks();
  
  // We need to make a copy of all the shared pointers because we modify the deque that
  //  allOrigPeaksDeque points at.
  const deque<PeakModel::PeakShrdPtr> allOrigPeaks( begin(*allOrigPeaksDeque), end(*allOrigPeaksDeque) );
  
  ContToPeakMap contToPeaks;
  for( const PeakModel::PeakShrdPtr &thispeak : allOrigPeaks )
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
      passMessage( WString::tr("err-no-add-peak-to-data-roi"), WarningWidget::WarningMsgInfo );
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
    MultiPeakFitChi2Fcn chi2fcn( static_cast<int>(origRoiPeaks.size()),
                                dataH,
                                peak->continuum()->type(),
                                PeakDef::SkewType::NoSkew,
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

  
  const MultiPeakInitialGuessMethod methods[] = { FromInputPeaks, UniformInitialGuess, FromDataInitialGuess };
  
  for( auto method : methods )
  {
    
    vector< std::shared_ptr<PeakDef> > orig_answer;
    try
    {
      for( auto p : answer )
        orig_answer.push_back( make_shared<PeakDef>( *p ) );
      
      findPeaksInUserRange( x0, x1, int(answer.size()), method, dataH,
                           m_dataMeasurement->detector(), answer, fitChi2 );
      
      std::vector<PeakDef> newRoiPeaks;
      for( size_t i = 0; i < answer.size(); ++i )
        newRoiPeaks.push_back( *answer[i] );
      
      MultiPeakFitChi2Fcn chi2fcn( static_cast<int>(newRoiPeaks.size()),
                                  dataH,
                                  peak->continuum()->type(),
                                  PeakDef::SkewType::NoSkew,
                                  lower_channel, upper_channel );
      fitChi2 = chi2fcn.evalRelBinRange( 0, chi2fcn.nbin(), newRoiPeaks );
      
      if( !IsNan(fitChi2) && !IsInf(fitChi2) && (fitChi2 < startingChi2) )
      {
        //cout << "Method " << method << " gave fitChi2=" << fitChi2 << ", which is better than initial startingChi2=" << startingChi2 << endl;
        break;
      }
      
      answer = orig_answer;
    }catch( std::exception & )
    {
      answer = orig_answer;
    }
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
    passMessage( WString::tr("err-add-peak-no-improve"), WarningWidget::WarningMsgInfo );
    return;
  }
  
  
  std::map<std::shared_ptr<PeakDef>,PeakModel::PeakShrdPtr> new_to_orig_peaks;
  for( const PeakModel::PeakShrdPtr &thispeak : allOrigPeaks )
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
  }//for( const PeakModel::PeakShrdPtr &thispeak : allOrigPeaks )
  
  
  for( size_t i = 0; i < answer.size(); ++i )
  {
    const bool isNew = !new_to_orig_peaks.count(answer[i]);
    InterSpec::addPeak( *answer[i], isNew );
  }
}//void addPeakFromRightClick()


void InterSpec::makePeakFromRightClickHaveOwnContinuum()
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  const std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak || !data
      || m_rightClickEnergy < peak->lowerX()
      || m_rightClickEnergy > peak->upperX() )
  {
    passMessage( WString::tr("err-no-roi-to-add-peak"), WarningWidget::WarningMsgInfo );
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
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  PeakSearchGuiUtils::change_continuum_type_from_right_click( this, m_rightClickEnergy,
                                                             continuum_type );
}//InterSpec::handleChangeContinuumTypeFromRightClick(...)


void InterSpec::handleChangeSkewTypeFromRightClick( const int peak_skew_type )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  PeakSearchGuiUtils::change_skew_type_from_right_click( this, m_rightClickEnergy,
                                                        peak_skew_type );
}//void handleChangeSkewTypeFromRightClick( const int peak_continuum_offset_type )


void InterSpec::shareContinuumWithNeighboringPeak( const bool shareWithLeft )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak || m_rightClickEnergy < peak->lowerX() || m_rightClickEnergy > peak->upperX() )
  {
    passMessage( WString::tr("err-no-roi-to-add-peak"), WarningWidget::WarningMsgInfo );
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
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  std::shared_ptr<const PeakDef> peak = nearestPeak( m_rightClickEnergy );
  if( !peak )
  {
    passMessage( WString::tr("err-no-peak-to-del"), WarningWidget::WarningMsgInfo );
    return;
  }

  m_peakModel->removePeak( peak );
}//void deletePeakFromRightClick()




void InterSpec::setPeakNuclide( const std::shared_ptr<const PeakDef> peak,
                                     std::string nuclide )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
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
  if( !p || !p->hasSourceGammaAssigned() )
    return;
  
  //Check if source is same as showing reference line, if so, set peak color to that.
  //ToDo: implement this if the user types into the Peak Manager.
  vector<ReferenceLineInfo> refLines;
  if( m_referencePhotopeakLines )
    refLines = m_referencePhotopeakLines->showingNuclides();
  
  // Try for a simple match - that is if the ref lines were for a nuclide, element, or reaction
  for( const auto &lines : refLines )
  {
    if( (lines.m_nuclide && (lines.m_nuclide==p->parentNuclide()))
       || (lines.m_element && (lines.m_element==p->xrayElement()))
       || (p->reaction() && lines.m_reactions.count(p->reaction()))
      )
    {
      index = m_peakModel->index( index.row(), PeakModel::kPeakLineColor );
      m_peakModel->setData( index, boost::any(WString(lines.m_input.m_color.cssText())) );
      
      return;
    }
  }//for( const auto &lines : refLines )
  
  // Lets try check each of the lines
  for( const auto &lines : refLines )
  {
    if( lines.m_nuclide || lines.m_element || !lines.m_reactions.empty() )
      continue;
    
    for( const ReferenceLineInfo::RefLine &line : lines.m_ref_lines )
    {
      if( (line.m_parent_nuclide && (line.m_parent_nuclide == p->parentNuclide()))
         || (line.m_element && (line.m_element == p->xrayElement()))
         || (line.m_reaction && (line.m_reaction == p->reaction()))
         )
      {
        index = m_peakModel->index( index.row(), PeakModel::kPeakLineColor );
        m_peakModel->setData( index, boost::any(WString(lines.m_input.m_color.cssText())) );
        
        return;
      }
    }//for( const ReferenceLineInfo::RefLine &line : lines.m_ref_lines )
      
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
    PopupDivMenuItem *item = menu->addMenuItem( WString::tr("rclick-mi-no-nuc-suggestions"), "", true );
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
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = m_peakModel->peaks();

#if( !USE_DETECTION_LIMIT_TOOL )
  if( !peaks || !peak )
    return;
#endif
  
  //We actually need to make a copy of peaks since we will post to outside the
  //  main thread where the deque is expected to remain constant, however, the
  //  PeakModel deque may get changed by like the user deleting a peak.
  if( peaks )
    peaks = make_shared< std::deque<shared_ptr<const PeakDef>>>( begin(*peaks), end(*peaks) );
  else
    peaks = make_shared<std::deque<shared_ptr<const PeakDef>>>();
  
  char energy_str[32] = { '\0' };
  snprintf( energy_str, sizeof(energy_str), "%.1f", energy );
  
  //see how many other peaks share ROI
  size_t npeaksInRoi = 0;
  if( peak )
  {
    for( const PeakModel::PeakShrdPtr &p : *peaks )
      npeaksInRoi += (p->continuum() == peak->continuum());
  }//if( peak )
  
  for( RightClickItems i = RightClickItems(0);
      i < kNumRightClickItems; i = RightClickItems(i+1) )
  {
    switch( i )
    {
      case kPeakEdit: 
      case kDeletePeak:
      case kAddPeakToRoi:
        m_rightClickMenutItems[i]->setHidden( !peak );
      break;
        
      case kChangeContinuum:
      {
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peak && m_rightClickChangeContinuumMenu )
        {
          // Disable current continuum type, enable all others
          const vector<WMenuItem *> items = m_rightClickChangeContinuumMenu->items();
          const char *labelTxt = PeakContinuum::offset_type_label_tr( peak->continuum()->type() );
          for( WMenuItem *item : items )
            item->setDisabled( item->text().key() == labelTxt );
        }//if( peak )
        
        break;
      }//case kChangeContinuum:
        
      case kChangeSkew:
      {
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peak && m_rightClickChangeSkewMenu )
        {
          // Disable current skew type, enable all others
          const vector<WMenuItem *> items = m_rightClickChangeSkewMenu->items();
          const char *labelTxt = PeakDef::to_label( peak->skewType() );
          for( WMenuItem *item : items )
            item->setDisabled( item->text() == labelTxt );
        }//if( peak && m_rightClickChangeSkewMenu )
      }//case kChangeSkew:
        
      case kChangeNuclide:
      {
        assert( m_rightClickNuclideSuggestMenu );
        
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peak && m_rightClickNuclideSuggestMenu )
        {
          for( WMenuItem *item : m_rightClickNuclideSuggestMenu->items() )
          {
            if( !item->hasStyleClass("PhoneMenuBack")
                && !item->hasStyleClass("PhoneMenuClose") )
              delete item;
          }//for( WMenuItem *item : m_rightClickNuclideSuggestMenu->items() )
          
          m_rightClickNuclideSuggestMenu->addMenuItem( WString::tr("rclick-mi-calc-status"),
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
        
            shared_ptr<const deque<shared_ptr<const PeakDef>>> hintpeaks;
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
        }//if( peak && m_rightClickNuclideSuggestMenu )
        
        break;
      }//case kChangeNuclide:
        
      case kRefitPeak:
        m_rightClickMenutItems[i]->setHidden( !peak || !peak->gausPeak() || (npeaksInRoi > 1) );
      break;
        
      case kRefitROI:
        m_rightClickMenutItems[i]->setHidden( !peak || !peak->gausPeak() || (npeaksInRoi < 2) );
      break;
        
      case kRefitPeakWithDrfFwhm:
        m_rightClickMenutItems[i]->setHidden( !peak || !peak->gausPeak() );
      break;
        
      case kSetMeanToRefPhotopeak:
      {
        const float energy = peak ? PeakSearchGuiUtils::reference_line_energy_near_peak( this, *peak ) : 0.0;
        const bool hide = (!peak || (energy < 10.0f));
        m_rightClickMenutItems[i]->setHidden( hide );
        if( !hide )
          m_rightClickMenutItems[i]->setText( WString::tr("rclick-mi-fix-energy").arg( energy_str ) );
        
        break;
      }//case kSetMeanToRefPhotopeak:
        
      case kShareContinuumWithLeftPeak:
      {
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peaks && peak )
        {
          auto iter = std::find( peaks->begin(), peaks->end(), peak );
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
        }//if( peak )
        
        break;
      }//case kShareContinuumWithLeftPeak:
        
        
      case kShareContinuumWithRightPeak:
      {
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peaks && peak )
        {
          auto iter = std::find( peaks->begin(), peaks->end(), peak );
          
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
          }//if( it already shares a continuum with the left peak )
          
          const double rightlower = rightpeak->lowerX();
          const double leftlower = peak->lowerX();
          const double leftupper = peak->upperX();
          
          //The below 2.0 is arbitrary
          const bool show = ((rightlower-leftupper) < 2.0*(leftupper-leftlower));
          m_rightClickMenutItems[i]->setHidden( !show );
        }//if( peaks && peak )
        
        break;
      }//case kShareContinuumWithRightPeak:
      
      
      case kMakeOwnContinuum:
      {
        m_rightClickMenutItems[i]->setHidden( !peak );
        
        if( peaks && peak )
        {
          bool shares = false;
          std::shared_ptr<const PeakContinuum> cont = peak->continuum();
          for( PeakModel::PeakShrdPtr p : *peaks )
            shares = (shares || (p!=peak && p->continuum()==cont) );
          m_rightClickMenutItems[i]->setHidden( !shares );
        }//if( peaks && peak )
        
        break;
      }//kMakeOwnContinuum
        
#if( USE_DETECTION_LIMIT_TOOL )
      case kFitNewPeakNotInRoi:
      {
        m_rightClickMenutItems[i]->setHidden( !!peak );
        
        if( !peak )
          m_rightClickMenutItems[i]->setText( WString::tr("rclick-mi-fit-new-peak").arg(energy_str) );
        
        break;
      }//case kFitNewPeakNotInRoi:
        
      case kAddPeakNotInRoi:
      {
        m_rightClickMenutItems[i]->setHidden( !!peak ); //"Add peak..."
        break;
      }//case kAddPeakNotInRoi:
      
      case kSearchEnergy:
      {
        m_rightClickMenutItems[i]->setHidden( !!peak );
        
        if( !peak )
          m_rightClickMenutItems[i]->setText( WString::tr("rclick-mi-search-energy").arg(energy_str) ); //"Search near {1} keV"
        
        break;
      }//case kSearchEnergy:
        
      case kSimpleMda:
      {
        m_rightClickMenutItems[i]->setHidden( !!peak );
        
        if( !peak )
        {
          WString target_txt;
          const tuple<const SandiaDecay::Nuclide *, double, float> near_line
                                = PeakSearchGuiUtils::nuclide_reference_line_near( this, energy );
          const SandiaDecay::Nuclide *ref_nuc = get<0>(near_line);
          const float ref_energy = get<2>(near_line);
            
          if( ref_nuc && (ref_energy > 10.0) )
          {
            char buffer[64] = { '\0' };
            snprintf( buffer, sizeof(buffer), "%s, %.1f keV", ref_nuc->symbol.c_str(), ref_energy );
            target_txt = WString::fromUTF8( buffer );
          }//if( ref_energy > 10.0 )
          
          if( target_txt.empty() )
            target_txt = WString("{1} keV").arg(energy_str);
          
          m_rightClickMenutItems[i]->setText( WString::tr("rclick-simple-mda").arg(target_txt) ); //"Quick MDA {1}..."
        }//if( !peak )
        
        break;
      }//case kSimpleMda:
#endif  //USE_DETECTION_LIMIT_TOOL
        
        
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
  auto create_editor = [this]( double ene ){
    m_peakEditWindow = new PeakEditWindow( ene, m_peakModel, this );
    m_peakEditWindow->editingDone().connect( this, &InterSpec::deletePeakEdit );
    m_peakEditWindow->finished().connect( this, &InterSpec::deletePeakEdit );
    m_peakEditWindow->resizeToFitOnScreen();
    
    PeakEdit *editor = m_peakEditWindow->peakEditor();
    if( !editor->isEditingValidPeak() )
    {
      delete m_peakEditWindow;
      m_peakEditWindow = nullptr;
      return;
    }//if( !editor->isEditingValidPeak() )
  };//
  
  auto change_peak = [this,create_editor]( double ene ){
    if( m_peakEditWindow )
    {
      PeakEdit *editor = m_peakEditWindow->peakEditor();
      editor->changePeak( ene );
    }else
    {
      create_editor( ene );
    }
  };//change_peak
  
  
  if( m_peakEditWindow )
  {
    PeakEdit *editor = m_peakEditWindow->peakEditor();
    const double prev_energy = editor->currentPeakEnergy();
    
    change_peak( energy );
    
    auto undo = [prev_energy,change_peak](){ change_peak( prev_energy ); };
    auto redo = [energy,change_peak](){ change_peak( energy ); };
    if( m_undo )
      m_undo->addUndoRedoStep( undo, redo, "Change peak editor peak." );
  }else
  {
    create_editor( energy );
    
    auto undo = [this](){
      if( m_peakEditWindow )
      {
        delete m_peakEditWindow;
        m_peakEditWindow = nullptr;
      }
    };
    
    auto redo = [energy,change_peak](){ change_peak( energy ); };
    if( m_undo )
      m_undo->addUndoRedoStep( undo, redo, "Open peak editor." );
  }//if( m_peakEditWindow ) / else
}//void createPeakEdit( Wt::WMouseEvent event )


PeakEditWindow *InterSpec::peakEdit()
{
  return m_peakEditWindow;
}


void InterSpec::deletePeakEdit()
{
  if( !m_peakEditWindow )
    return;
  
  PeakEdit *editor = m_peakEditWindow->peakEditor();
  assert( editor );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    const double currentEnergy = editor ? editor->currentPeakEnergy() : 0.0;
    
    auto undo = [this, currentEnergy](){
      m_peakEditWindow = new PeakEditWindow( currentEnergy, m_peakModel, this );
      m_peakEditWindow->editingDone().connect( this, &InterSpec::deletePeakEdit );
      m_peakEditWindow->finished().connect( this, &InterSpec::deletePeakEdit );
      m_peakEditWindow->resizeToFitOnScreen();
    };//auto undo
    
    auto redo = [this](){
      deletePeakEdit();
    };
    
    m_undo->addUndoRedoStep( undo, redo, "Close peak editor." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
  delete m_peakEditWindow;
  m_peakEditWindow = nullptr;
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
    if( (fabs(peak->mean()-energy) < (3.0/2.35482)*width) )
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
  const bool wasShown = m_featureMarkersShown[static_cast<int>(option)];
  
  m_featureMarkersShown[static_cast<int>(option)] = show;
  m_spectrum->setFeatureMarkerOption( option, show );
  
  if( m_undo && m_undo->canAddUndoRedoNow() && (show != wasShown) )
  {
    auto undo = [this,option,show](){
      FeatureMarkerWidget *tool = m_featureMarkersWindow ? m_featureMarkersWindow->tool() : nullptr;
      if( !tool && m_referencePhotopeakLines )
        tool = m_referencePhotopeakLines->featureMarkerTool();
        
      if( tool )
        tool->setFeatureMarkerChecked( option, !show );
      setFeatureMarkerOption( option, !show );
    };
    auto redo = [this,option,show](){
      FeatureMarkerWidget *tool = m_featureMarkersWindow ? m_featureMarkersWindow->tool() : nullptr;
      if( !tool && m_referencePhotopeakLines )
        tool = m_referencePhotopeakLines->featureMarkerTool();
      
      if( tool )
        tool->setFeatureMarkerChecked( option, show );
      setFeatureMarkerOption( option, show );
    };
    
    m_undo->addUndoRedoStep( undo, redo, "Toggle feature marker" );
  }//if( m_undo && (show != wasShown) )
}//setFeatureMarkerOption(...)


bool InterSpec::showingFeatureMarker( const FeatureMarkerType option )
{
  return m_featureMarkersShown[static_cast<int>(option)];
}


void InterSpec::setComptonPeakAngle( const int angle )
{
  const int prev_angle = m_spectrum->comptonPeakAngle();
  m_spectrum->setComptonPeakAngle( angle );
  
  if( m_undo && m_undo->canAddUndoRedoNow() && (prev_angle != angle) )
  {
    auto undo = [this,prev_angle](){
      FeatureMarkerWidget *tool = m_featureMarkersWindow ? m_featureMarkersWindow->tool() : nullptr;
      if( !tool && m_referencePhotopeakLines )
        tool = m_referencePhotopeakLines->featureMarkerTool();
      
      if( tool )
        tool->setDisplayedComptonPeakAngle( prev_angle );
      m_spectrum->setComptonPeakAngle( prev_angle );
    };
    auto redo = [this,angle](){
      FeatureMarkerWidget *tool = m_featureMarkersWindow ? m_featureMarkersWindow->tool() : nullptr;
      if( !tool && m_referencePhotopeakLines )
        tool = m_referencePhotopeakLines->featureMarkerTool();
      
      if( tool )
        tool->setDisplayedComptonPeakAngle( angle );
      m_spectrum->setComptonPeakAngle( angle );
    };
    
    m_undo->addUndoRedoStep( undo, redo, "Change Compton angle" );
  }//if( m_undo && (prev_angle != angle) )
}//void setComptonPeakAngle( const float angle );


void InterSpec::toggleFeatureMarkerWindow()
{
  const bool showing = (m_featureMarkersWindow
                  || (m_referencePhotopeakLines && m_referencePhotopeakLines->featureMarkerTool()));
  
  displayFeatureMarkerWindow( !showing );
}//void toggleFeatureMarkerWindow()


void InterSpec::displayFeatureMarkerWindow( const bool show )
{
  if( m_undo && m_undo->canAddUndoRedoNow() )
    m_undo->addUndoRedoStep( [this,show](){ displayFeatureMarkerWindow(!show); },
                            [this,show](){ displayFeatureMarkerWindow(show); },
                            "Show feature marker window" );
  
  const bool showing = (m_featureMarkersWindow
                  || (m_referencePhotopeakLines && m_referencePhotopeakLines->featureMarkerTool()));
  
  if( showing == show )
  {
    if( m_referencePhotopeakLines && m_referencePhotopeakLines->featureMarkerTool() )
    {
      assert( m_toolsTabs );
      if( m_toolsTabs )
      {
        const int index = m_toolsTabs->indexOf(m_referencePhotopeakLines);
        m_toolsTabs->setCurrentIndex( index );
        m_currentToolsTab = index;
      }
      m_referencePhotopeakLines->emphasizeFeatureMarker();
    }//if( reference photopeak tool is showing feature marker options )
    
    return;
  }//if( showing == show )
  
  if( !show )
  {
    if( m_featureMarkersWindow )
    {
      AuxWindow::deleteAuxWindow( m_featureMarkersWindow );
      m_featureMarkersWindow = nullptr;
    }else
    {
      assert( m_referencePhotopeakLines && m_referencePhotopeakLines->featureMarkerTool() );
      m_referencePhotopeakLines->removeFeatureMarkerTool();
    }
    
    for( FeatureMarkerType i = FeatureMarkerType(0);
         i < FeatureMarkerType::NumFeatureMarkers;
         i = FeatureMarkerType(static_cast<int>(i)+1) )
    {
      if( m_featureMarkersShown[static_cast<int>(i)] )
        m_spectrum->setFeatureMarkerOption( i, false );
    }
    
    m_featureMarkerMenuItem->setText( WString::tr("app-mi-view-feature-markers") );
    
    return;
  }//if( !show )
  
  
  m_featureMarkerMenuItem->setText( WString::tr("app-mi-view-hide-feature-markers") );
  
  // Initing the feature marker tool can add some undo/redo steps, by calling into
  //  `InterSpec::setComptonPeakAngle(...)` or `InterSpec::displayFeatureMarkerWindow(bool)`
  //  - lets block these calls.
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  FeatureMarkerWidget *widget = nullptr;
  if( toolTabsVisible() && m_referencePhotopeakLines )
  {
    widget = m_referencePhotopeakLines->showFeatureMarkerTool();
    assert( m_toolsTabs );
    if( m_toolsTabs )
      m_toolsTabs->setCurrentIndex( m_toolsTabs->indexOf(m_referencePhotopeakLines) );
    m_referencePhotopeakLines->emphasizeFeatureMarker();
  }else
  {
    m_featureMarkersWindow = new FeatureMarkerWindow( this );
    m_featureMarkersWindow->finished().connect( boost::bind(&InterSpec::displayFeatureMarkerWindow, this, false) );
    
    widget = m_featureMarkersWindow->tool();
  }//if( toolTabsVisible() && m_referencePhotopeakLines ) / else
  
  // Restore the widget state to previous opened state.
  for( FeatureMarkerType i = FeatureMarkerType(0);
      i < FeatureMarkerType::NumFeatureMarkers;
      i = FeatureMarkerType(static_cast<int>(i)+1) )
  {
    const bool show = m_featureMarkersShown[static_cast<int>(i)];
    m_spectrum->setFeatureMarkerOption( i, show );
    widget->setFeatureMarkerChecked( i, show );
  }//for( set FeatureMarkers to the last state of the window )
}//void displayFeatureMarkerWindow()



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


#if( USE_DB_TO_STORE_SPECTRA )
void InterSpec::saveStateToDb( Wt::Dbo::ptr<UserState> entry )
{
  if( !entry || (entry->user != m_user) || !m_user.session() )
    throw runtime_error( "Invalid input to saveStateToDb()" );

  try
  {
    saveShieldingSourceModelToForegroundSpecMeas();
#if( USE_REL_ACT_TOOL )
    saveRelActManualStateToForegroundSpecMeas();
    saveRelActAutoStateToForegroundSpecMeas();
#endif
    
    DataBaseUtils::DbTransaction transaction( *m_sql );
    entry.modify()->serializeTime = WDateTime::currentDateTime();
    
    if( !entry->creationTime.isValid() )
      entry.modify()->creationTime = entry->serializeTime;
    
    const bool deepCopy = ((entry->stateType != UserState::kUserSaved) || !!entry->snapshotTagParent);
    
    std::shared_ptr<const SpecMeas> foreground = measurment( SpecUtils::SpectrumType::Foreground );
    std::shared_ptr<const SpecMeas> second = measurment( SpecUtils::SpectrumType::SecondForeground );
    std::shared_ptr<const SpecMeas> background = measurment( SpecUtils::SpectrumType::Background );
    
    Dbo::ptr<UserFileInDb> dbforeground, dbsecond, dbbackground;
    dbforeground = measurementFromDb( SpecUtils::SpectrumType::Foreground, true );
    if( foreground != second )
      dbsecond = measurementFromDb( SpecUtils::SpectrumType::SecondForeground, true );
    if( (background != foreground) && (background != second) )
      dbbackground = measurementFromDb( SpecUtils::SpectrumType::Background, true );
  
    //JIC, make sure indices have all been assigned to everything.
    m_sql->session()->flush();
    
    auto create_copy_or_update_in_db = [this, deepCopy, &entry]( Dbo::ptr<UserFileInDb> &dbfile, shared_ptr<const SpecMeas> &file ){
      if( !dbfile || !file )
        return;
    
      // We dont need to make a copy of the file if it was not part of a save state, but now its
      //  becoming part of one
      if( deepCopy && (dbfile->isPartOfSaveState || dbfile->snapshotParent) )
        dbfile = UserFileInDb::makeDeepCopyOfFileInDatabase( dbfile, *m_sql, true );
      
      dbfile.modify()->isPartOfSaveState = true;
      
      Dbo::ptr<UserFileInDbData> filedata;
      for( auto iter = dbfile->filedata.begin(); iter != dbfile->filedata.end(); ++iter )
      {
        assert( !filedata );
        if( filedata )
          (*iter).remove();
        else
          filedata = *iter;
      }
      
      //assert( filedata );
      if( !filedata )
      {
        // We may get here if we've removed a state (i.e., a `UserState::kEndOfSessionTemp`) from
        //  the database, but the pointer is still in memory.
        UserFileInDbData *newdata = new UserFileInDbData();
        filedata.reset( newdata );
        newdata->fileInfo = dbfile;
        m_sql->session()->add( filedata );
      }
      
      // This next line will throw an exception if the file is too large to store in database
      filedata.modify()->setFileData( file, UserFileInDbData::SerializedFileFormat::k2012N42 );
    };//create_copy_or_update_in_db lambda
    
    create_copy_or_update_in_db( dbforeground, foreground );
    create_copy_or_update_in_db( dbsecond, second );
    create_copy_or_update_in_db( dbbackground, background );

    if( !dbforeground && foreground )
      throw runtime_error( "Error saving foreground to the database" );  //not displayed to user, so not internationalized
    if( !dbsecond && second )
      throw runtime_error( "Error saving second foreground to the database" );
    if( !dbbackground && background )
      throw runtime_error( "Error saving background to the database" );
    
    //We need to make sure dbforeground, dbbackground, dbsecond will have been
    //  written to the database, so there id()'s will be not -1.
    m_sql->session()->flush();
    
    entry.modify()->foregroundId = static_cast<int>( dbforeground.id() );
    entry.modify()->backgroundId = static_cast<int>( dbbackground.id() );
    entry.modify()->secondForegroundId = static_cast<int>( dbsecond.id() );
    
    // If using is clicking save - remove any `kUserStateAutoSavedWork` states, as
    //  they would now be behind the state the user is saving
    if( entry->stateType == UserState::kUserSaved )
    {
      Dbo::collection<Dbo::ptr<UserState>> eosEntries = entry->snapshotTags.find()
        .where( "StateType = ?" )
        .bind(int(UserState::UserStateType::kUserStateAutoSavedWork));
      for( auto iter = eosEntries.begin(); iter != eosEntries.end(); ++iter )
        iter->remove();
    }//if( remove `kUserStateAutoSavedWork` states )
    
    
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
    
    // TODO: should add time chart limits here - if showing
    
    entry.modify()->shownDisplayFeatures = 0x0;
    if( toolTabsVisible() )
      entry.modify()->shownDisplayFeatures |= UserState::kDockedWindows;
    if( m_spectrum->yAxisIsLog() )
      entry.modify()->shownDisplayFeatures |= UserState::kLogSpectrumCounts;
    
    try
    {
      if( UserPreferences::preferenceValue<bool>( "ShowVerticalGridlines", this ) )
        entry.modify()->shownDisplayFeatures |= UserState::kVerticalGridLines;
      
      if( UserPreferences::preferenceValue<bool>( "ShowHorizontalGridlines", this ) )
        entry.modify()->shownDisplayFeatures |= UserState::kHorizontalGridLines;
    }catch(...)
    {
      // We shouldnt get here
      assert( 0 );
    }
    
    if( m_spectrum->legendIsEnabled() )
      entry.modify()->shownDisplayFeatures |= UserState::kSpectrumLegend;
    
    if( m_shieldingSourceFit )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingShieldSourceFit;
    
#if( USE_TERMINAL_WIDGET )
    if( m_terminal )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingTerminalWidget;
#endif
    
#if( USE_REL_ACT_TOOL )
    if( m_relActManualGui )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingRelActManual;
    
    if( m_relActAutoGui )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingRelActAuto;
#endif
    
    if( m_multimedia )
      entry.modify()->shownDisplayFeatures |= UserState::kShowingMultimedia;
    
    if( m_gammaXsToolWindow && m_gammaXsToolWindow->xstool() )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingGammaXsTool;
      entry.modify()->gammaXsToolUri = m_gammaXsToolWindow->xstool()->encodeStateToUrl();
    }
    
    if( m_doseCalcWindow && m_doseCalcWindow->tool() )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingDoseCalcTool;
      entry.modify()->doseCalcToolUri = m_doseCalcWindow->tool()->encodeStateToUrl();
    }
    
    if( m_1overR2Calc )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowing1OverR2Tool;
      entry.modify()->oneOverR2ToolUri = m_1overR2Calc->encodeStateToUrl();
    }
    
    if( m_fluxTool )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingFluxTool;
      entry.modify()->fluxToolUri = m_fluxTool->encodeStateToUrl();
    }//if( m_fluxTool )
    
    if( m_unitsConverter )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingUnitConvertTool;
      entry.modify()->unitsConverterToolUri = m_unitsConverter->encodeStateToUrl();
    }
    
    if( m_decayInfoWindow )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingNucDecayInfo;
      entry.modify()->nucDecayInfoUri = m_decayInfoWindow->encodeStateToUrl();
    }
    
    if( m_gammaCountDialog )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingEnergyRangeSum;
      entry.modify()->energyRangeSumUri = m_gammaCountDialog->encodeStateToUrl();
    }
    
#if( USE_DETECTION_LIMIT_TOOL )
    if( m_detectionLimitWindow )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingDetectionSens;
      entry.modify()->detectionSensitivityToolUri = m_simpleMdaWindow->tool()->encodeStateToUrl();
    }//if( m_detectionLimitWindow )
    
    if( m_simpleMdaWindow )
    {
      entry.modify()->shownDisplayFeatures |= UserState::kShowingSimpleMda;
      entry.modify()->simpleMdaUri = m_simpleMdaWindow->tool()->encodeStateToUrl();
    }//if( m_simpleMdaWindow )
#endif
    
    entry.modify()->backgroundSubMode = UserState::kNoSpectrumSubtract;
    if( m_spectrum->backgroundSubtract() )
      entry.modify()->backgroundSubMode = UserState::kBackgorundSubtract;
    
    entry.modify()->currentTab = UserState::kNoTabs;
    if( m_toolsTabs )
    {
      const WString &txtKey = m_toolsTabs->tabText( m_toolsTabs->currentIndex() ).key();
      if( txtKey == PeakInfoTabTitleKey )
        entry.modify()->currentTab = UserState::kPeakInfo;
      else if( txtKey == GammaLinesTabTitleKey )
        entry.modify()->currentTab = UserState::kGammaLines;
      else if( txtKey == CalibrationTabTitleKey )
        entry.modify()->currentTab = UserState::kCalibration;
      else if( txtKey == NuclideSearchTabTitleKey )
        entry.modify()->currentTab = UserState::kIsotopeSearch;
      else if( txtKey == FileTabTitleKey )
        entry.modify()->currentTab = UserState::kFileTab;
#if( USE_TERMINAL_WIDGET )
      else if( txtKey == TerminalTabTitleKey )
        entry.modify()->currentTab = UserState::kTerminalTab;
#endif
#if( USE_REL_ACT_TOOL )
      else if( txtKey == RelActManualTitleKey )
        entry.modify()->currentTab = UserState::kRelActManualTab;
#endif
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
      const int colorThemIndex = UserPreferences::preferenceValue<int>("ColorThemeIndex", this);
      if( colorThemIndex > 0 )
        theme = m_user->m_colorThemes.find().where( "id = ?" ).bind( colorThemIndex ).resultValue();
      if( theme )
        entry.modify()->colorThemeJson = theme->json_data;
    }catch(...)
    {
      //Probably shouldnt ever happen
      cerr << "saveStateToDb(): Caught exception retrieving color theme from database" << endl;
    }//try / catch
    
    transaction.commit();
  }catch( FileToLargeForDbException &e )
  {
    cerr << "Caught: " << e.what() << endl;
    throw;
  }catch( std::exception &e )
  {
    cerr << "saveStateToDb(...) caught: " << e.what() << endl;
    throw runtime_error( WString::tr("err-save-sate-to-db").toUTF8() );
  }//try / catch
}//void saveStateToDb( Wt::Dbo::ptr<UserState> entry )


void InterSpec::loadStateFromDb( Wt::Dbo::ptr<UserState> entry )
{
  //This function takes about 5 ms (to load the previous HPGe state, on new app
  //  construction) on my 2.6 GHz Intel Core i7, as of (20150316).
  if( !entry )
    throw runtime_error( "UserState was invalid" );

  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  try
  {
    closeShieldingSourceFit();
    
    switch( entry->stateType )
    {
      case UserState::kEndOfSessionTemp:
      case UserState::kUserSaved:
      case UserState::kForTest:
      case UserState::kUndefinedStateType:
      case UserState::kUserStateAutoSavedWork:
      break;
    }//switch( entry->stateType )

    DataBaseUtils::DbTransaction transaction( *m_sql );
    
    Wt::Dbo::ptr<UserState> parent = entry->snapshotTagParent;
    while( parent && parent->snapshotTagParent )
      parent = parent->snapshotTagParent;
    
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
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    bool wasDocked = (entry->shownDisplayFeatures & UserState::kDockedWindows);
    if( isPhone() )
    {
      wasDocked = (m_toolsTabs != nullptr);
    }else if( toolTabsVisible() != wasDocked )
    {
      setToolTabsVisible( wasDocked );
    }
#else
    bool wasDocked = (entry->shownDisplayFeatures & UserState::kDockedWindows);
    if( toolTabsVisible() != wasDocked )
      setToolTabsVisible( wasDocked );
#endif
    
    //Now start reloading the state
    std::shared_ptr<SpecMeas> foreground, second, background;
    std::shared_ptr<SpecMeas> snapforeground, snapsecond, snapbackground;
    std::shared_ptr<SpectraFileHeader> foregroundheader, backgroundheader,
                                          secondheader;
    
    
    if( parent )
    {
      //If we are loading a snapshot, we actually want to associate the
      //  foreground/second/background entries in the database with the parent
      //  entries, but actually load the snapshots.  This is so when the user
      //  saves changes, the parents saved state receives the changes, and the
      //  snapshot doesnt get changed.
      //All of this is really ugly and horrible, and should be improved!
      
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
          dbbackground = UserFileInDb::makeDeepCopyOfFileInDatabase( dbbackground, *m_sql, true );
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
          dbsecond = UserFileInDb::makeDeepCopyOfFileInDatabase( dbsecond, *m_sql, true );
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
    
    setSpectrum( nullptr, {}, SpecUtils::SpectrumType::Background, 0 );
    setSpectrum( nullptr, {}, SpecUtils::SpectrumType::SecondForeground, 0 );
    
    Wt::WFlags<SetSpectrumOptions> options;
#if( USE_REMOTE_RID )
    if( background || second )
      options |= SetSpectrumOptions::SkipExternalRid;
#endif
    
    setSpectrum( foreground, foregroundNums, SpecUtils::SpectrumType::Foreground, options );
    if( foreground )
    {
      //If we dont have a foreground, we probably shouldnt be loading the state, but...
#if( USE_REMOTE_RID )
      if( !second )
        options.clear( SetSpectrumOptions::SkipExternalRid );
#endif
      setSpectrum( background, backgroundNums, SpecUtils::SpectrumType::Background, options );
      
#if( USE_REMOTE_RID )
      options.clear( SetSpectrumOptions::SkipExternalRid );
#endif
      setSpectrum( second, secondNums, SpecUtils::SpectrumType::SecondForeground, options );
    }//if( foreground )
    
    
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
//        shieldingSourceFit();
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
        
    
    if( m_multimedia )
      programmaticallyCloseMultimediaWindow();
    assert( !m_multimedia );
    
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
    
    if( (entry->shownDisplayFeatures & UserState::kShowingShieldSourceFit) )
    {
      shieldingSourceFit();
      
      const double windowWidth = 0.95 * renderedWidth();
      const double windowHeight = 0.95 * renderedHeight();
      m_shieldingSourceFitWindow->resizeWindow( windowWidth, windowHeight );
      
      m_shieldingSourceFitWindow->resizeToFitOnScreen();
      m_shieldingSourceFitWindow->show();
      m_shieldingSourceFitWindow->centerWindow();
    }
    
#if( USE_TERMINAL_WIDGET )
    if( (entry->shownDisplayFeatures & UserState::kShowingTerminalWidget) )
      createTerminalWidget();
#endif
    
#if( USE_REL_ACT_TOOL )
    if( (entry->shownDisplayFeatures & UserState::kShowingRelActManual) )
      createRelActManualWidget();
    
    if( (entry->shownDisplayFeatures & UserState::kShowingRelActAuto) )
      showRelActAutoWindow();
#endif
    
    if( (entry->shownDisplayFeatures & UserState::kShowingMultimedia) )
      showMultimedia( SpecUtils::SpectrumType::Foreground );
    
    
    if( (entry->shownDisplayFeatures & UserState::kShowingGammaXsTool)
       && !entry->gammaXsToolUri.empty() )
    {
      GammaXsWindow *w = showGammaXsTool();
      if( w && w->xstool() )
        w->xstool()->handleAppUrl( entry->gammaXsToolUri );
    }
    
    
    if( (entry->shownDisplayFeatures & UserState::kShowingDoseCalcTool)
       && !entry->doseCalcToolUri.empty() )
    {
      DoseCalcWindow *w = showDoseTool();
      DoseCalcWidget *tool = w ? w->tool() : nullptr;
      if( tool )
      {
        string path, uri = entry->doseCalcToolUri;
        const auto pos = uri.find('?');
        if( pos != string::npos )
        {
          path = uri.substr(0,pos);
          uri = uri.substr(pos+1);
        }
        tool->handleAppUrl( path, uri );
      }//if( tool )
    }
    
    if( (entry->shownDisplayFeatures & UserState::kShowing1OverR2Tool)
       && !entry->oneOverR2ToolUri.empty() )
    {
      OneOverR2Calc *calc = createOneOverR2Calculator();
      if( calc )
        calc->handleAppUrl( entry->oneOverR2ToolUri );
    }
    
    if( (entry->shownDisplayFeatures & UserState::kShowingFluxTool)
       && !entry->fluxToolUri.empty() )
    {
      FluxToolWindow *tool = createFluxTool();
      if( tool )
        tool->handleAppUrl( entry->fluxToolUri );
    }//if( m_fluxTool )
    
    if( (entry->shownDisplayFeatures & UserState::kShowingUnitConvertTool)
       && !entry->unitsConverterToolUri.empty() )
    {
      UnitsConverterTool *converter = createUnitsConverterTool();
      if( converter )
        converter->handleAppUrl( entry->unitsConverterToolUri );
    }
    
    if( (entry->shownDisplayFeatures & UserState::kShowingNucDecayInfo)
       && !entry->nucDecayInfoUri.empty() )
    {
      DecayWindow *decay = createDecayInfoWindow();
      if( decay )
        decay->handleAppUrl( entry->nucDecayInfoUri );
    }
    
    if( (entry->shownDisplayFeatures & UserState::kShowingEnergyRangeSum)
       && !entry->energyRangeSumUri.empty() )
    {
      GammaCountDialog *dialog = showGammaCountDialog();
      if( dialog )
        dialog->handleAppUrl( entry->energyRangeSumUri );
    }
    
    
#if( USE_DETECTION_LIMIT_TOOL )
    if( (entry->shownDisplayFeatures & UserState::kShowingDetectionSens)
       && !entry->detectionSensitivityToolUri.empty() )
    {
      showDetectionLimitTool( entry->detectionSensitivityToolUri );
    }//if( m_detectionLimitWindow )
    
    if( (entry->shownDisplayFeatures & UserState::kShowingSimpleMda)
       && !entry->simpleMdaUri.empty() )
    {
      auto dialog = showSimpleMdaWindow();
      if( dialog )
        dialog->tool()->handleAppUrl( entry->simpleMdaUri );
    }//if( m_simpleMdaWindow )
#endif
    
    
    if( wasDocked )
    {
      WString titleKey;
      switch( entry->currentTab )
      {
        case UserState::kPeakInfo:        titleKey = PeakInfoTabTitleKey;      break;
        case UserState::kGammaLines:      titleKey = GammaLinesTabTitleKey;    break;
        case UserState::kCalibration:     titleKey = CalibrationTabTitleKey;   break;
        case UserState::kIsotopeSearch:   titleKey = NuclideSearchTabTitleKey; break;
        case UserState::kFileTab:         titleKey = FileTabTitleKey;          break;
#if( USE_TERMINAL_WIDGET )
        case UserState::kTerminalTab:     titleKey = TerminalTabTitleKey;      break;
#endif
#if( USE_REL_ACT_TOOL )
        case UserState::kRelActManualTab: titleKey = RelActManualTitleKey;     break;
#endif
        case UserState::kNoTabs:                                               break;
      };//switch( entry->currentTab )
      
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
      if( isPhone() )
      {
        switch( entry->currentTab )
        {
          case UserState::kPeakInfo:
            m_toolsTabs->setCurrentIndex( 1 ); //m_peakInfoDisplay
            break;
          case UserState::kGammaLines:
            m_toolsTabs->setCurrentIndex( 2 ); //m_referencePhotopeakLines
            break;
          case UserState::kCalibration:
            m_toolsTabs->setCurrentIndex( 3 ); //m_energyCalTool
            break;
          case UserState::kIsotopeSearch:  
            m_toolsTabs->setCurrentIndex( 4 ); //m_nuclideSearchContainer
            break;
          case UserState::kFileTab:
            m_toolsTabs->setCurrentIndex( 0 ); //CompactFileManager
            break;
  #if( USE_TERMINAL_WIDGET )
          case UserState::kTerminalTab:
            if( m_terminal )
              m_toolsTabs->setCurrentWidget( m_terminal );
            break;
  #endif
  #if( USE_REL_ACT_TOOL )
          case UserState::kRelActManualTab: 
            if( m_relActManualGui )
              m_toolsTabs->setCurrentWidget( m_relActManualGui );
            break;
  #endif
          case UserState::kNoTabs:  
            // Leave it where it was
            break;
        };//switch( entry->currentTab )
      }else
      {
        for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
        {
          if( m_toolsTabs->tabText(tab).key() == titleKey )
          {
            // Note: this causes handleToolTabChanged(...) to be called
            m_toolsTabs->setCurrentIndex( tab );
            m_currentToolsTab = tab;
            break;
          }//if( m_toolsTabs->tabText( tab ) == title )
        }//for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
      }
#else
      for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
      {
        if( m_toolsTabs->tabText(tab).key() == titleKey )
        {
          // Note: this causes handleToolTabChanged(...) to be called
          m_toolsTabs->setCurrentIndex( tab );
          m_currentToolsTab = tab;
          break;
        }//if( m_toolsTabs->tabText( tab ) == title )
      }//for( int tab = 0; tab < m_toolsTabs->count(); ++tab )
#endif
    }//if( wasDocked )
    
    
    //Should take care of entry->showingWindows here, so we can rely on it below
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

    // We wont change user preferences if this state was an end of session state.
    //  This is mostly for development: when developing we often kill process so
    //  the state doesnt get saved, but then any preferences we've updated will
    //  essentially get lost when auto-loading previous end of session state,
    //  which is a bit annoying
    if( entry->userOptionsJson.size()
       && (entry->stateType != UserState::kEndOfSessionTemp) )
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
          const WString &name_wstr = obj.get("name");
          const string name = name_wstr.toUTF8();
          
          // Users probably wont expect loading a state to mess up preferences for how they
          //  currently like things, so we we'll only set the ones that notably impact presentation
          //  of the data, and are hopefully obvious
          static const std::string prefs_to_load[] = {
            "LogY", "ShowVerticalGridlines", "ShowHorizontalGridlines",
            "ColorThemeIndex", "ShowXAxisSlider", "CompactXAxis",
            "ShowYAxisScalers", "DisplayBecquerel"
          };
          
          if( std::find(begin(prefs_to_load), end(prefs_to_load), name) == end(prefs_to_load) )
            continue;
          
          const int type = obj.get("type");
          const UserOption::DataType datatype = UserOption::DataType( type );
          
          switch( datatype )
          {
            case UserOption::String:
            {
              const std::string value = obj.get("value");
              UserPreferences::setPreferenceValue( name, value, this );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Decimal:
            {
              const double value = obj.get("value");
              UserPreferences::setPreferenceValue( name, value, this );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Integer:
            {
              const int value = obj.get("value");
              UserPreferences::setPreferenceValue( name, value, this );
              break;
            }//case UserOption::Boolean:
              
            case UserOption::Boolean:
            {
              const bool value = obj.get("value");
              UserPreferences::setPreferenceValue( name, value, this );
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
    
    transaction.commit();
    
    const long long int db_index = parent ? parent.id() : entry.id();
    if( m_dataMeasurement && !m_displayedSamples.empty() )
      m_dataMeasurement->setDbStateId( db_index, m_displayedSamples );
    
    updateSaveWorkspaceMenu();
  }catch( std::exception &e )
  {
    if( m_dataMeasurement )
      m_dataMeasurement->setDbStateId( -1, m_displayedSamples );
    
    updateSaveWorkspaceMenu();
    cerr << "\n\n\n\nloadStateFromDb(...) caught: " << e.what() << endl;
    throw runtime_error( e.what() );
  }//try / catch
}//void loadStateFromDb( Wt::Dbo::ptr<UserState> entry )


//Whenever `m_dataMeasurement->dbStateId(db_index, m_displayedSamples)` changes value,
//  we need to update the menu!
void InterSpec::updateSaveWorkspaceMenu()
{
  m_saveStateAs->setDisabled( !m_dataMeasurement );
  
  //check if Save menu needs to be updated
  
  
  const long long int db_index = (m_dataMeasurement && !m_displayedSamples.empty())
                                    ? m_dataMeasurement->dbStateId(m_displayedSamples)
                                    : static_cast<long long int>(-1);
  
  if( db_index >= 0 )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      Dbo::ptr<UserState> state = m_sql->session()
                                       ->find<UserState>( "where id = ?" )
                                       .bind( db_index );
      transaction.commit();
      
      if( state && (state->stateType==UserState::kUserStateAutoSavedWork
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


long long int InterSpec::currentAppStateDbId()
{
  if( !m_dataMeasurement || m_displayedSamples.empty() )
    return -1;
  
  return m_dataMeasurement->dbStateId(m_displayedSamples);
}//int currentAppStateDbId()


void InterSpec::resetCurrentAppStateDbId()
{
  if( m_dataMeasurement )
    m_dataMeasurement->setDbStateId( -1, m_displayedSamples );
  updateSaveWorkspaceMenu();
}//void resetCurrentAppStateDbId()
#endif //#if( USE_DB_TO_STORE_SPECTRA )


std::shared_ptr<const ColorTheme> InterSpec::getColorTheme()
{
  if( m_colorTheme )
    return m_colorTheme;
  
  auto theme = make_shared<ColorTheme>();
  
  try
  {
    const int colorThemIndex = UserPreferences::preferenceValue<int>("ColorThemeIndex", this);
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

  
  m_spectrum->applyColorTheme( theme );
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
    // We'll explicitly check if the CSS file exists on disk... and we'll generous and check both
    //  from the CWD, and the docRoot() directory (e.g., InterSpec_resources), although I guess
    //  we should only do the latter.
    //  Also, we "know" we should have the dark theme, so we wont touch disk for this.
    if( (theme->nonChartAreaTheme == "dark")
        || SpecUtils::is_file( SpecUtils::append_path(wApp->docRoot(), cssfile) )
        || SpecUtils::is_file(cssfile) )
    {
      wApp->useStyleSheet( cssfile );
      m_currentColorThemeCssFile = cssfile;
    }else
    {
      passMessage( "Could not load the non-peak-area theme '" + theme->nonChartAreaTheme + "'",
                  WarningWidget::WarningMsgLow );
    }
  }//if( should load CSS file )
  
  m_colorThemeChanged.emit( theme );
}//void InterSpec::applyColorTheme()


Wt::Signal< std::shared_ptr<const ColorTheme> > &InterSpec::colorThemeChanged()
{
  return m_colorThemeChanged;
}


#if( BUILD_AS_OSX_APP || APPLY_OS_COLOR_THEME_FROM_JS || IOS || BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
void InterSpec::osThemeChange( std::string name )
{
  SpecUtils::to_lower_ascii(name);
  if( name != "light" && name != "dark" && name != "default" && name != "no-preference" && name != "no-support" )
  {
    cerr << "InterSpec::osThemeChange('" << name << "'): unknown OS theme." << endl;
    return;
  }
  
  //TODO: if( name == "no-support" ), then should hide widgets associated with "AutoDarkFromOs"
  
  try
  {
    const bool autoDark = UserPreferences::preferenceValue<bool>( "AutoDarkFromOs", this );
    
    if( autoDark && (name == "dark") )
    {
      //Check to see if we already have dark applied
      if( m_colorTheme && SpecUtils::icontains( m_colorTheme->theme_name.toUTF8(), "Dark") )
        return;
     
      unique_ptr<ColorTheme> theme = ColorTheme::predefinedTheme( ColorTheme::PredefinedColorTheme::DarkColorTheme );
      assert( theme );
      if( theme )
        applyColorTheme( make_shared<ColorTheme>(*theme) );
    }else if( !autoDark || name == "light" || name == "default" || name == "no-preference" || name == "no-support" )
    {
      if( m_colorTheme
         && (m_colorTheme->dbIndex < 0)
         && (ColorTheme::PredefinedColorTheme(-m_colorTheme->dbIndex) == ColorTheme::PredefinedColorTheme::DarkColorTheme) )
      {
        cout << "Will set to default or user specified color theme, from \"Dark\"." << endl;
        
        m_colorTheme = nullptr;
        applyColorTheme( nullptr );
      }else
      {
        cout << "Current theme is not \"Dark\", so not setting color theme." << endl;
      }
    }
  }catch( std::exception &e )
  {
    cerr << "InterSpec::osThemeChange() caught exception - not doing anything:" << e.what() << endl;
  }
}//void osThemeChange( std::string name )
#endif


void InterSpec::showColorThemeWindow()
{
  //Takes ~5ms (on the server) to create a ColorThemeWindow.
	ColorThemeWindow *window = new ColorThemeWindow(this);
  
  new UndoRedoManager::BlockGuiUndoRedo( window ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
}//void showColorThemeWindow()


void InterSpec::showLicenseAndDisclaimersWindow()
{
  if( m_licenseWindow )
  {
    m_licenseWindow->show();
    return;
  }
  
  m_licenseWindow = new LicenseAndDisclaimersWindow( this );
  
  m_licenseWindow->finished().connect( this, &InterSpec::deleteLicenseAndDisclaimersWindow );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ deleteLicenseAndDisclaimersWindow(); };
    auto redo = [this](){ showLicenseAndDisclaimersWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo),
                            "Show disclaimers, credits, and contact window." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void showLicenseAndDisclaimersWindow()


void InterSpec::startClearSession()
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( !app )
    return;
  
  SimpleDialog *window = new SimpleDialog( WString::tr("clear-session"),
                                           WString::tr("clear-session-msg") );
  
  WPushButton *button = window->addButton( WString::tr("Yes") );
  button->clicked().connect( app, &InterSpecApp::clearSession );
  
  button = window->addButton( WString::tr("No") );
  button->setFocus();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto closer = wApp->bind( boost::bind(&WDialog::accept, window) );
    // We wont have redo show the dialog again, because then the closer functionoid wont work, and
    //  it isnt worth making the SimpleDialog a member variable of the InterSpec class.
    m_undo->addUndoRedoStep( closer, nullptr, "Show clear session dialog" );
  }//if( undo )
}//void startClearSession()


UserPreferences *InterSpec::preferences()
{
  return m_preferences;
}


const Wt::Dbo::ptr<InterSpecUser> &InterSpec::user()
{
  return m_user;
}

void InterSpec::reReadUserInfoFromDb()
{
  m_user.reread();
}

void InterSpec::deleteLicenseAndDisclaimersWindow()
{ 
  if( !m_licenseWindow )
    return;
  
  AuxWindow::deleteAuxWindow( m_licenseWindow );
  m_licenseWindow = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ showLicenseAndDisclaimersWindow(); };
    auto redo = [this](){ deleteLicenseAndDisclaimersWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo),
                            "Close disclaimers, credits, and contact window." );
  }//if( add undo step )
}//void deleteLicenseAndDisclaimersWindow()


void InterSpec::showWelcomeDialogWorker( const bool force )
{
  if( m_useInfoWindow )
    return;
  
  std::function<void(bool)> dontShowAgainCallback;
  if( !force )
  {
    dontShowAgainCallback = [this](bool value){
      UserPreferences::setPreferenceValue( "ShowSplashScreen", value, this );
    };
  }//if( !force )
  
  m_useInfoWindow = new UseInfoWindow( dontShowAgainCallback , this );

  m_useInfoWindow->finished().connect( boost::bind( &InterSpec::deleteWelcomeDialog, this, force ) );
  
  if( force && m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ deleteWelcomeDialog(false); };
    auto redo = [this,force](){ showWelcomeDialog(true); };
    m_undo->addUndoRedoStep( undo, redo, "Show welcome dialog." );
  }//if( force && m_undo && m_undo->canAddUndoRedoNow() )
  
  wApp->triggerUpdate();
}//void showWelcomeDialogWorker( bool force );


void InterSpec::showWelcomeDialog( const bool force )
{
  if( m_useInfoWindow )
  {
    m_useInfoWindow->show();
    m_useInfoWindow->resizeToFitOnScreen();
    m_useInfoWindow->centerWindowHeavyHanded();
    return;
  }//if( m_useInfoWindow )
  
  try
  {
    if( !force && !UserPreferences::preferenceValue<bool>("ShowSplashScreen", this) )
      return;
  }catch(...)
  {
    assert(0);
    //m_user didnt have preference "ShowSplashScreen" for some reason
    if( !force )
      return;
  }
  
  if( force )
  {
    showWelcomeDialogWorker( force );
  }else
  {
    /*
     For Android, showing this useInfoWindow at startup causes some exceptions
     in the JavaScript whenever the loading indicator is shown.  I'm pretty
     confused by it.
     I tried setting a new loading indicator in deleteWelcomeDialog(), but it didnt work
     */
    WServer::instance()->post( wApp->sessionId(), std::bind(&InterSpec::showWelcomeDialogWorker, this, force ) );
  }
}//void showWelcomeDialog()


void InterSpec::deleteWelcomeDialog( const bool addUndoRedoStep )
{
  if( !m_useInfoWindow )
    return;
  
  AuxWindow::deleteAuxWindow( m_useInfoWindow );
  m_useInfoWindow = nullptr;
  
  if( addUndoRedoStep && m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ showWelcomeDialog(true); };
    auto redo = [this](){ deleteWelcomeDialog(false); };
    m_undo->addUndoRedoStep( undo, redo, "Close welcome dialog." );
  }//if( force && m_undo && m_undo->canAddUndoRedoNow() )
}//void deleteWelcomeDialog()


void InterSpec::deleteEnergyCalPreserveWindow()
{
  if( m_preserveCalibWindow )
  {
    delete m_preserveCalibWindow;
    m_preserveCalibWindow = nullptr;
  }
}//void deleteEnergyCalPreserveWindow()



ExportSpecFileWindow *InterSpec::createExportSpectrumFileDialog()
{
  assert( !m_exportSpecFileWindow || (m_undo && m_undo->isInUndoOrRedo()) );
  
  if( m_exportSpecFileWindow )
    return m_exportSpecFileWindow;
  
  m_exportSpecFileWindow = new ExportSpecFileWindow( this );
  m_exportSpecFileWindow->finished().connect( this, &InterSpec::handleExportSpectrumFileDialogClose );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){
      if( m_exportSpecFileWindow )
        m_exportSpecFileWindow->accept();
    };
    
    auto redo = [this](){
      createExportSpectrumFileDialog();
    };
    
    m_undo->addUndoRedoStep( undo, redo, "Show spectrum file export tool." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
  return m_exportSpecFileWindow;
}//void createExportSpectrumFileDialog()


void InterSpec::handleExportSpectrumFileDialogClose()
{
  assert( m_exportSpecFileWindow );
  if( !m_exportSpecFileWindow )
    return;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    const string uri = m_exportSpecFileWindow->encodeStateToUrl();
    auto undo = [this,uri](){
      createExportSpectrumFileDialog();
      if( m_exportSpecFileWindow )
        m_exportSpecFileWindow->handleAppUrl( uri );
    };
    
    auto redo = [this](){
      if( m_exportSpecFileWindow )
        m_exportSpecFileWindow->accept();
    };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close spectrum file export tool." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
  m_exportSpecFileWindow = nullptr;
}


GammaCountDialog *InterSpec::showGammaCountDialog()
{
  if( m_gammaCountDialog )
    return m_gammaCountDialog;

  m_gammaCountDialog = new GammaCountDialog( this );
  m_gammaCountDialog->finished().connect( this, &InterSpec::deleteGammaCountDialog );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    m_undo->addUndoRedoStep( [=](){deleteGammaCountDialog();},
                            [=](){showGammaCountDialog();},
                            "Show Energy Range Sum." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
  return m_gammaCountDialog;
}//GammaCountDialog *showGammaCountDialog()


void InterSpec::deleteGammaCountDialog()
{
  if( !m_gammaCountDialog )
    return;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    const string uri = m_gammaCountDialog->encodeStateToUrl();
    auto undo = [this,uri](){
      GammaCountDialog *dialog = showGammaCountDialog();
      if( dialog )
        dialog->handleAppUrl( uri );
    };
    
    auto redo = [this](){ deleteGammaCountDialog(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Energy Range Sum." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
  AuxWindow::deleteAuxWindow( m_gammaCountDialog );
  m_gammaCountDialog = nullptr;
}//void deleteGammaCountDialog()




#if( USE_SPECRUM_FILE_QUERY_WIDGET )
void InterSpec::showFileQueryDialog()
{
  if( m_specFileQueryDialog )
    return;
  
  
  m_specFileQueryDialog = new AuxWindow( WString::tr("window-title-spec-file-query"), 
                                        Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                        | AuxWindowProperties::SetCloseable );
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
  
  WPushButton *closeButton = m_specFileQueryDialog->addCloseButtonToFooter( WString::tr("Close"), true );
  closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_specFileQueryDialog ) );
  
  m_specFileQueryDialog->finished().connect( this, &InterSpec::deleteFileQueryDialog );
  
  if( m_renderedWidth > 100 && m_renderedHeight > 100 && !isPhone() )
  {
    m_specFileQueryDialog->resize( 0.85*m_renderedWidth, 0.85*m_renderedHeight );
    m_specFileQueryDialog->centerWindow();
  }
  
  AuxWindow::addHelpInFooter( m_specFileQueryDialog->footer(), "spectrum-file-query" );
  
  new UndoRedoManager::BlockGuiUndoRedo( m_specFileQueryDialog ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
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


void InterSpec::logMessage( const Wt::WString& message, int priority )
{
#if( PERFORM_DEVELOPER_CHECKS )
  static std::mutex s_message_mutex;
  
  {//begin codeblock to logg message
    std::lock_guard<std::mutex> file_gaurd( s_message_mutex );
    ofstream output( "interspec_messages_to_users.txt", ios::out | ios::app );
    auto now = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    output << "Message " << SpecUtils::to_iso_string( now ) << " ";
    output << "[" << priority << "]: ";
    output << message.toUTF8() << endl << endl;
  }//end codeblock to logg message
#endif
  
  if( wApp )
  {
    m_messageLogged.emit( message, priority );
//    wApp->triggerUpdate();
  }else
  {
    auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
    cerr << "Message " << SpecUtils::to_iso_string( now ) << " ";
    cerr << "[" << priority << "]: ";
    cerr << message.toUTF8() << endl << endl;
  }
}//void InterSpec::logMessage(...)


WarningWidget *InterSpec::warningWidget()
{
  return m_warnings;
}

Wt::Signal< Wt::WString, int > &InterSpec::messageLogged()
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
    m_warningsWindow = new AuxWindow( WString::tr("window-title-notification-log"),
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                   | AuxWindowProperties::DisableCollapse
                   | AuxWindowProperties::EnableResize
                   | AuxWindowProperties::SetCloseable) );
    m_warningsWindow->contents()->setOffsets( WLength(0, WLength::Pixel), Wt::Left | Wt::Top );
    m_warningsWindow->rejectWhenEscapePressed();
    m_warningsWindow->stretcher()->addWidget( m_warnings, 0, 0, 1, 1 );
    m_warningsWindow->stretcher()->setContentsMargins(0,0,0,0);
    //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
    //  size to be the min-size
    m_warningsWindow->setMinimumSize( 640, 480 );
    m_warningsWindow->finished().connect( boost::bind( &InterSpec::handleWarningsWindowClose, this ) );
        
      
    WPushButton *clearButton = new WPushButton( WString::tr("notification-log-clear"),
                                               m_warningsWindow->footer() );
    clearButton->clicked().connect( boost::bind( &WarningWidget::clearMessages, m_warnings ) );
    clearButton->addStyleClass( "BinIcon" );
    if( isMobile() )
      clearButton->setFloatSide( Right );
    WPushButton *closeButton = m_warningsWindow->addCloseButtonToFooter();
    closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_warningsWindow ) );
  }//if( !m_warningsWindow )
  
  m_warningsWindow->show();
  m_warningsWindow->resizeScaledWindow(0.75, 0.75);
  m_warningsWindow->centerWindow();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ handleWarningsWindowClose(); };
    auto redo = [this](){ showWarningsWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show notifications window" );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void showWarningsWindow()


void InterSpec::handleWarningsWindowClose()
{
  if( m_warningsWindow )
  {
    m_warningsWindow->stretcher()->removeWidget( m_warnings );
    AuxWindow::deleteAuxWindow( m_warningsWindow );
    m_warningsWindow = nullptr;
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){ showWarningsWindow(); };
      auto redo = [this](){ handleWarningsWindowClose(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close notifications window" );
    }//if( m_undo && m_undo->canAddUndoRedoNow() )
  }//if( m_warningsWindow )
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
    m_peakInfoWindow = new AuxWindow( WString::tr("window-title-peak-manager"),
                              Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable) );
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
    
    WPushButton* closeButton = m_peakInfoWindow->addCloseButtonToFooter(WString::tr("Close"),true);
    
    closeButton->clicked().connect( boost::bind( &AuxWindow::hide, m_peakInfoWindow ) );
      
    AuxWindow::addHelpInFooter( m_peakInfoWindow->footer(), "peak-manager" );
    
    WPushButton *b = new WPushButton( WString::tr(CalibrationTabTitleKey), footer );
    b->clicked().connect( this, &InterSpec::showEnergyCalWindow );
    b->setFloatSide(Wt::Right);
      
    b = new WPushButton( WString::tr(GammaLinesTabTitleKey), footer );
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
      m_toolsTabs->addTab( m_peakInfoDisplay, WString::tr(PeakInfoTabTitleKey), TabLoadPolicy );
      m_toolsTabs->setCurrentIndex( m_toolsTabs->indexOf(m_peakInfoDisplay) );
    }//if( m_toolsTabs->indexOf(m_peakInfoDisplay) < 0 )
    m_currentToolsTab = m_toolsTabs->currentIndex();
    
    delete m_peakInfoWindow;
    m_peakInfoWindow = NULL;
  }//if( m_toolsTabs )
}//void handlePeakInfoClose()


#if( USE_DB_TO_STORE_SPECTRA )
void InterSpec::finishStoreStateInDb( const WString &name,
                                      const WString &desc,
                                      Dbo::ptr<UserState> parent )
{
  if( name.empty() )
  {
    passMessage( WString::tr("db-error-no-name"), WarningWidget::WarningMsgHigh );
    return;
  }//if( name.empty() )
  
  
  Dbo::ptr<UserState> answer;
  
  UserState *state = new UserState();
  state->user = m_user;
  state->stateType = UserState::kUserSaved;
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
  }catch( FileToLargeForDbException &e )
  {
    passMessage( e.message(), WarningWidget::WarningMsgHigh );
  }catch( std::exception &e )
  {
    passMessage( e.what(), WarningWidget::WarningMsgHigh );
  }
   
#if( PERFORM_DEVELOPER_CHECKS )
  {// Begin check if we are leaving any orphaned spectra that were part of a user state
    DataBaseUtils::DbTransaction transaction( *m_sql );
    
    Dbo::collection< Dbo::ptr<UserFileInDb> > query;
    const char *stateQueryTxt = "A.InterSpecUser_id = ? AND A.StateType = ? AND A.SnapshotTagParent_id IS NULL";
    query = m_sql->session()->query<Dbo::ptr<UserFileInDb>>(
              "SELECT ufd FROM UserFileInDb ufd"
              " LEFT JOIN UserState us ON ufd.id IN (us.ForegroundId, us.BackgroundId, us.SecondForegroundId)"
              " WHERE us.id IS NULL AND ufd.IsPartOfSaveState = 1 AND ufd.InterSpecUser_id = ?" ).bind( m_user.id() );
    
    const size_t num_dangling = query.size();
    if( num_dangling )
    {
      // I think this can happen for end-of-session states - but not really sure
      cerr << "There are " << query.size() << " dangling files" << endl;
      for( auto iter = query.begin(); iter != query.end(); ++iter )
        cerr << (*iter)->filename << endl;
      
      cerr << endl;
      
      //The following query is untested - but I think it would cleanup the database of dangling files.
      //m_sql->session()->execute(
      //  "DELETE FROM UserFileInDb WHERE id IN ("
      //    " SELECT ufd.id FROM UserFileInDb ufd LEFT JOIN UserState us"
      //    " ON ufd.id IN (us.ForegroundId, us.BackgroundId, us.SecondForegroundId)"
      //    " WHERE us.id IS NULL AND ufd.IsPartOfSaveState = 1 AND ufd.InterSpecUser_id = ?"
      //  ")"
      //).bind( m_user.id() );
    }//if( num_dangling )
    transaction.commit();
  }// End check if we are leaving any orphaned spectra that were part of a user state
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  
  if( parent )
  {
    WString msg = WString::tr("db-new-tag");
    passMessage( msg.arg(name).arg(parent.get()->name), WarningWidget::WarningMsgInfo );
  }else
  {
    WString msg = WString::tr("db-new-snapshot");
    passMessage( msg.arg(name), WarningWidget::WarningMsgSave );
  }
  
  const long long int db_index = parent ? parent.id() : answer.id();
  if( m_dataMeasurement )
    m_dataMeasurement->setDbStateId( db_index, m_displayedSamples );
  updateSaveWorkspaceMenu();
}//void finishStoreStateInDb()


#if( INCLUDE_ANALYSIS_TEST_SUITE )
void InterSpec::storeTestStateToN42( const Wt::WString name, const Wt::WString descr  )
{
  try
  {
    if( !m_dataMeasurement )
      throw runtime_error( "No valid foreground file" );
    
    // Open the output stream
    string filepath = SpecUtils::append_path(TEST_SUITE_BASE_DIR, "analysis_tests");
      
    if( !SpecUtils::is_directory( filepath ) )
      throw runtime_error( "CWD didnt contain a 'analysis_tests' folder as expected" );
      
    const int offset = wApp->environment().timeZoneOffset();
    auto localtime = chrono::time_point_cast<chrono::microseconds>(chrono::system_clock::now());
    localtime += std::chrono::seconds(60*offset);
    
    string timestr = SpecUtils::to_iso_string( localtime );
    string::size_type period_pos = timestr.find_last_of( '.' );
    timestr = timestr.substr( 0, period_pos );
    
    string filename = name.toUTF8() + "_" + timestr + ".n42";
    SpecUtils::erase_any_character( filename, "\\/:?\"<>|" );
    
    if( name.empty() )
      throw runtime_error( "Name must NOT be empty." );
    
    filename += "_" + timestr + ".n42";
    
    filepath = SpecUtils::append_path( filepath, filename );
      
  #ifdef _WIN32
    const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16(filepath);
    ofstream output( wfilepath.c_str(), ios::binary|ios::out );
  #else
    ofstream output( filepath.c_str(), ios::binary|ios::out );
  #endif
    if( !output.is_open() )
      throw runtime_error( "Couldnt open test file output '" + filepath + "'" );
    
    
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
    preferences()->userOptionsToXml( InterSpecNode, this );
    
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
    
    
    const char *val = n42doc->allocate_string( timestr.c_str(), timestr.size()+1 );
    xml_node<char> *node = n42doc->allocate_node( node_element, "TestSaveDateTime", val );
    InterSpecNode->append_node( node );
    
    if( !name.empty() )
    {
      const string txt = name.toUTF8();
      val = n42doc->allocate_string( txt.c_str(), txt.size()+1 );
      node = n42doc->allocate_node( node_element, "TestName", val );
      InterSpecNode->append_node( node );
    }//if( name.size() )
    
    if( !descr.empty() )
    {
      const string txt = descr.toUTF8();
      val = n42doc->allocate_string( txt.c_str(), txt.size()+1 );
      node = n42doc->allocate_node( node_element, "TestDescription", val );
      InterSpecNode->append_node( node );
    }//if( name.size() )
    
    string xml_data;
    rapidxml::print( std::back_inserter(xml_data), *n42doc, 0 );
    
    if( !output.write( xml_data.c_str(), xml_data.size() ) )
      throw runtime_error( "writing to file failed." );
  }catch( std::exception &e )
  {
    string msg = "Failed to save test state to N42 file: ";
    msg += e.what();
    passMessage( msg, WarningWidget::WarningMsgHigh );
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
      UserPreferences::restoreUserPrefsFromXml( preferences, this );
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
      shieldingSourceFit();
      m_shieldingSourceFit->deSerialize( sourcefit );
      closeShieldingSourceFit();
    }//if( sourcefit )
    
    stringstream msg;
    msg << "Loaded test state";
    if( name && name->value_size() )
      msg << " named '" << name->value() << "'";
    if( savedate && savedate->value_size() )
      msg << " saved " << savedate->value();
    if( descrip && descrip->value_size() )
      msg << " with description " << descrip->value();
    
    passMessage( msg.str(), WarningWidget::WarningMsgInfo );
  }catch( std::exception &e )
  {
    string msg = "InterSpec::loadTestStateFromN42: ";
    msg += e.what();
    passMessage( msg, WarningWidget::WarningMsgHigh );
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
                 "directory", WarningWidget::WarningMsgHigh );
    return;
  }//if( files.empty() )
  
  AuxWindow *window = new AuxWindow( "Test State N42 Files", 
                                    (WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                      | AuxWindowProperties::DisableCollapse) );
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



void InterSpec::startStoreTestState()
{
  AuxWindow *window = new AuxWindow( "Store app test state to N42",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                     | AuxWindowProperties::TabletNotFullScreen
                                     | AuxWindowProperties::DisableCollapse) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  WGridLayout *layout = window->stretcher();
  WLineEdit *edit = new WLineEdit();
  edit->setEmptyText( WString::tr("db-name-empty-txt") );
  
  edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  edit->setAttributeValue( "autocorrect", "off" );
  edit->setAttributeValue( "spellcheck", "off" );
#endif
  
  WText *label = new WText( WString::tr("Name") );
  layout->addWidget( label, 0, 0 );
  layout->addWidget( edit,  0, 1 );
  
  if( !!m_dataMeasurement && m_dataMeasurement->filename().size() )
    edit->setText( m_dataMeasurement->filename() );
  
  WTextArea *summary = new WTextArea();
  label = new WText( WString::tr("Desc.") );
  summary->setEmptyText( WString::tr("db-dec-empty-txt") );
  layout->addWidget( label, 1, 0 );
  layout->addWidget( summary, 1, 1 );
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 1, 1 );
  
  
  WPushButton *closeButton = window->addCloseButtonToFooter( WString::tr("Cancel"), false);
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  WPushButton *save = new WPushButton( WString::tr("Create"), window->footer() );
  save->setIcon( "InterSpec_resources/images/disk2.png" );
  
  save->clicked().connect( std::bind( [this, edit, summary, window](){
    const WString name = edit ? edit->text() : WString();
    const WString sumtxt = summary ? summary->text() : WString();
    storeTestStateToN42( name, sumtxt );
  } ) );
  
  window->setMinimumSize(WLength(450), WLength(250));
  
  window->centerWindow();
  window->show();
}//void startStoreTestState()
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE )


void InterSpec::stateSave()
{
  Dbo::ptr<UserState> state;
  const long long int current_state_index = currentAppStateDbId();
  if( current_state_index >= 0 )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      state = m_sql->session()->find<UserState>( "where id = ?" )
        .bind( current_state_index );
      
      // find ultimate parent state
      while( state && state->snapshotTagParent )
        state = state->snapshotTagParent;
  
      if( state && (state->stateType != UserState::kUserSaved) )
        state.modify()->stateType = UserState::kUserSaved;
      
      transaction.commit();
    }catch( std::exception &e )
    {
      assert( 0 );
    }//try / catch
  }//if( current_state_index >= 0 )
    
  
  if( state )
  {
    try
    {
      saveStateToDb( state );
    }catch( FileToLargeForDbException &e )
    {
      passMessage( e.message(), WarningWidget::WarningMsgHigh );
    }catch( std::exception &e )
    {
      passMessage( e.what(), WarningWidget::WarningMsgHigh );
    }
  }else //No currentStateID, so just save as new state
  {
    stateSaveAs();
  }
}//void stateSave()


void InterSpec::stateSaveAs()
{
  AuxWindow *window = new AuxWindow( WString::tr("window-title-store-state-as"),
    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
      | AuxWindowProperties::TabletNotFullScreen
      | AuxWindowProperties::DisableCollapse) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  window->disableCollapse(); //Not sure why the property flags dont get this
  WGridLayout *layout = window->stretcher();
    
  WLineEdit *edit = new WLineEdit();
  edit->setEmptyText( WString::tr("db-name-empty-txt") );
  
  edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  edit->setAttributeValue( "autocorrect", "off" );
  edit->setAttributeValue( "spellcheck", "off" );
#endif

  
  WText *label = new WText( WString::tr("Name") );
  layout->addWidget( label, 2, 0 );
  layout->addWidget( edit,  2, 1 );
    
  if( !!m_dataMeasurement && m_dataMeasurement->filename().size() )
      edit->setText( m_dataMeasurement->filename() );
    
  WTextArea *summary = new WTextArea();
  label = new WText( WString::tr("Desc.") );
  summary->setEmptyText( WString::tr("db-dec-empty-txt") );
  layout->addWidget( label, 3, 0 );
  layout->addWidget( summary, 3, 1 );
  
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 3, 1 );
  
  label = new WText( WString::tr("db-dec-export-msg") );
  label->setInline( false );
  label->setTextAlignment( Wt::AlignCenter );
  layout->addWidget( label, 4, 1 );
  
  WPushButton *closeButton = window->addCloseButtonToFooter( WString::tr("Cancel"), false);
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  WPushButton *save = new WPushButton( WString::tr("Save") , window->footer() );
  save->setIcon( "InterSpec_resources/images/disk2.png" );
  
  save->clicked().connect( boost::bind( &InterSpec::stateSaveAsFinish,
                                        this, edit, summary, window ) );
    
  window->setMinimumSize(WLength(450), WLength(250));
    
  window->centerWindow();
  window->show();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto closer = wApp->bind( boost::bind( &AuxWindow::hide, window ) );
    m_undo->addUndoRedoStep( closer, nullptr, "Show save-as dialog." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void stateSaveAs()


//New method to save tag for snapshot
void InterSpec::stateSaveTag()
{
  AuxWindow *window = new AuxWindow( WString::tr("window-title-tag-state"),
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                   | AuxWindowProperties::TabletNotFullScreen) );
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->setClosable( false );
  window->disableCollapse();
  WGridLayout *layout = window->stretcher();
    
  WLineEdit *edit = new WLineEdit();
  
  edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  edit->setAttributeValue( "autocorrect", "off" );
  edit->setAttributeValue( "spellcheck", "off" );
#endif
  
  if( !isMobile() )
    edit->setFocus(true);
  WText *label = new WText( WString::tr("Name") );
  layout->addWidget( label, 2, 0 );
  layout->addWidget( edit,  2, 1 );
    
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 2, 1 );
  
  WPushButton *closeButton = window->addCloseButtonToFooter( WString::tr("Cancel"), false );
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  WPushButton *save = new WPushButton( WString::tr("db-Tag"), window->footer() );
  //save->setIcon( "InterSpec_resources/images/disk2.png" );
    
  save->clicked().connect( boost::bind( &InterSpec::stateSaveTagFinish,
                                         this, edit, window ) );
    
  window->setMinimumSize(WLength(450), WLength::Auto);
    
  window->centerWindow();
  window->show();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto closer = wApp->bind( boost::bind( &AuxWindow::hide, window ) );
    m_undo->addUndoRedoStep( closer, nullptr, "Show save-tag dialog." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void stateSaveTag()


void InterSpec::stateSaveTagFinish( WLineEdit *nameedit, AuxWindow *window )
{
  Dbo::ptr<UserState> state;
    
  const long long int current_state_index = currentAppStateDbId();
  
  if( current_state_index >= 0 )
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_sql );
      state = m_sql->session()->find<UserState>( "where id = ?" )
                               .bind( current_state_index );
      transaction.commit();
    }catch( std::exception &e )
    {
      Wt::log("error") << "startStoreStateInDb() caught: " << e.what();
    }//try / catch
  }//if( current_state_index >= 0 )
    
  //Save Snapshot/environment
  assert( state ); // Should always be able to get the state
  if( state )
  {
    WString name = nameedit ? nameedit->text() : WString();
    finishStoreStateInDb( name, "", state );
  }else
  {
    // Shouldnt ever get here
    passMessage( "You must first save a state before creating a snapshot.",
                WarningWidget::WarningMsgHigh );
  }
      
  //close window
  AuxWindow::deleteAuxWindow( window );
}//stateSaveTagFinish()


void InterSpec::stateSaveAsFinish( WLineEdit *nameedit,
                                      WTextArea *descriptionarea,
                                      AuxWindow *window )
{
  // If the files in the DB are currently part of a saved-state, we need to duplicate
  //  them, otherwise `finishStoreStateInDb` will just assume they are part of the measurement
  //  and overwrite them - so in this case we'll just disconnect
  //  Note that in this case `finishStoreStateInDb(...)` --> `saveStateToDb(...)` will
  //  create a new entry in the database.
  SpectraFileModel *filemodel = m_fileManager->model();
  for( int i = 0; i < 3; i++ )
  {
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
    auto m = measurment( type );
    if( !m )
      continue;
    
    shared_ptr<SpectraFileHeader> hdr = filemodel->fileHeader( m );
    assert( hdr );
    if( !hdr )
      continue;
    
    Wt::Dbo::ptr<UserFileInDb> dbmeas = hdr->dbEntry();
    
    if( dbmeas && (dbmeas->isPartOfSaveState || dbmeas->snapshotParent) )
      hdr->clearDbEntry();
  }//for( int i = 0; i < 3; i++ )
  
  //Save Snapshot
  const WString name = nameedit->text();
  const WString description = descriptionarea ? descriptionarea->text() : WString();
  
  finishStoreStateInDb( name, description, Dbo::ptr<UserState>() );
    
  //close window
  AuxWindow::deleteAuxWindow(window);
}//stateSaveAsFinish()


void InterSpec::saveStateAtForegroundChange( const bool doAsync )
{
  const bool saveState = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", this );
  if( !saveState || !m_dataMeasurement )
  {
    saveShieldingSourceModelToForegroundSpecMeas();
  #if( USE_REL_ACT_TOOL )
    saveRelActManualStateToForegroundSpecMeas();
    saveRelActAutoStateToForegroundSpecMeas();
  #endif
    
    return;
  }//if( !saveState || !m_dataMeasurement )
  
  Wt::log( "info" ) << "Auto-saving state for foreground change for file: '" << m_dataMeasurement->filename() << "'.";

  // TODO: possibly check file size, and if we dont have a hope of saving to DB, i.e.,
  //        if( m_dataMeasurement->memmorysize() > UserFileInDb::sm_maxFileSizeBytes )
  //          return;
  
  const int offset = wApp->environment().timeZoneOffset();
  WString desc = "Working Point";
  const auto now = chrono::system_clock::now() + chrono::seconds( 60*offset );
  WString name = SpecUtils::to_common_string( chrono::time_point_cast<chrono::microseconds>(now), true ); //"9-Sep-2014 15:02:15"
  
  DataBaseUtils::DbTransaction transaction( *m_sql );
  
  try
  {
    const long long int current_state_index = currentAppStateDbId();
    Wt::Dbo::ptr<UserState> parentState;
    if( current_state_index >= 0 )
    {
      parentState = m_sql->session()->find<UserState>( "where id = ?" )
            .bind( current_state_index );
    }
    
    // If we are connected to a state in the database, we'll save a `kUserStateAutoSavedWork` state,
    //  otherwise, we'll save just the files to the database
    if( parentState )
    {
      // Remove previous `kUserStateAutoSavedWork` states for this user-saved state - we will make a new one
      Dbo::collection<Dbo::ptr<UserState>> prevEosStates
        = parentState->snapshotTags.find()
        .where( "StateType = ?" )
        .bind( int(UserState::UserStateType::kUserStateAutoSavedWork) );
      
      for( auto iter = prevEosStates.begin(); iter != prevEosStates.end(); ++iter )
        UserState::removeFromDatabase( *iter, *m_sql );
      
      // We will only save a state if the foreground file has been modified
      if( m_dataMeasurement && m_dataMeasurement->modified() )
      {
        UserState *state = new UserState();
        state->user = m_user;
        state->stateType = UserState::kUserStateAutoSavedWork;
        state->creationTime = WDateTime::currentDateTime();
        state->name = name;
        state->description = desc;
        state->snapshotTagParent = parentState;
        
        Wt::Dbo::ptr<UserState> dbstate = m_sql->session()->add( state );
        
        saveStateToDb( dbstate );
        
        Wt::log("debug") << "saveStateAtForegroundChange(): Saved kUserStateAutoSavedWork"
        " state to database";
      }else
      {
        Wt::log("debug") << "saveStateAtForegroundChange(): not saving kUserStateAutoSavedWork"
        " state to database, since there are no changes since last user save.";
      }//if( we have done work we might want to save )
    }else
    {
      saveShieldingSourceModelToForegroundSpecMeas();
#if( USE_REL_ACT_TOOL )
      saveRelActManualStateToForegroundSpecMeas();
      saveRelActAutoStateToForegroundSpecMeas();
#endif
      
      // Get foreground/background/secondary files, and update them to the database..
      SpectraFileModel *filemodel = m_fileManager->model();
      shared_ptr<SpectraFileHeader> forehdr, backhdr, secohdr;
      
      forehdr = filemodel->fileHeader( m_dataMeasurement );
      if( m_backgroundMeasurement != m_dataMeasurement )
        backhdr = filemodel->fileHeader( m_backgroundMeasurement );
      if( (m_secondDataMeasurement != m_backgroundMeasurement)
         && (m_secondDataMeasurement != m_dataMeasurement) )
        secohdr = filemodel->fileHeader( m_secondDataMeasurement );
      
      auto save_to_db = [doAsync]( std::weak_ptr<SpectraFileHeader> wk_ptr, std::shared_ptr<SpecMeas> meas ){
        auto call_save = [wk_ptr, meas](){
          shared_ptr<SpectraFileHeader> hdr = wk_ptr.lock();
          
          if( !hdr )
          {
            cerr << "SpectraFileHeader no longer in memory - not writing to DB." << endl;
            return;
          }
          
          try
          {
            hdr->saveToDatabase( meas );
          }catch( FileToLargeForDbException &e )
          {
            // Expected problem - dont give an error message
          }catch( std::exception &e )
          {
            Wt::log("error") << "InterSpec::saveStateAtForegroundChange() async error: " << e.what();
          }
        };//call_save lambda
        
        assert( wk_ptr.lock() );
        if( doAsync )
        {
          WServer::instance()->schedule( 25, wApp->sessionId(), call_save, [](){
            cerr << "Failed to save spectrum to DB in worker (session dead?)." << endl;
            assert( 0 );
          } );
        }else
        {
          call_save();
        }
      };//save_to_db lambda
      
      // If we load a spectrum file, and get the dialog that lets us pick back up from a
      //  previous working of the file, and we do, `forehdr` may be nullptr, because
      //  #SpectraFileHeader::setDbEntry has not been called for the original file yet
      if( forehdr )
        save_to_db( forehdr, m_dataMeasurement );
      
      if( backhdr )
        save_to_db( backhdr, m_backgroundMeasurement );
      
      if( secohdr )
        save_to_db( secohdr, m_secondDataMeasurement );
      
      Wt::log("debug") << "saveStateAtForegroundChange(): Instead kUserStateAutoSavedWork state,"
      " saved spectrum files to database since not in a user save-state.";
    }//if( current_state_index >= 0 ) / else
    
    transaction.commit();
  }catch( FileToLargeForDbException &e )
  {
    // Expected problem - dont give an error message - but I dont think we ever end up here
    Wt::log("info") << "InterSpec::saveStateAtForegroundChange() Not saving state: " << e.what();
    transaction.rollback();
  }catch( std::exception &e )
  {
    Wt::log("error") << "InterSpec::saveStateAtForegroundChange() error: " << e.what();
    transaction.rollback();
  }
}//void saveStateAtForegroundChange()


void InterSpec::saveStateForEndOfSession()
{
  const bool saveState = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", this );
  Wt::log( "info" ) << "Will auto-save state for session-id: '" << wApp->sessionId() << "' at end of session.";

  const int offset = wApp->environment().timeZoneOffset();
  WString desc = "End of Session";
  const auto now = chrono::system_clock::now() + chrono::seconds( 60*offset );
  WString name = SpecUtils::to_common_string( chrono::time_point_cast<chrono::microseconds>(now), true ); //"9-Sep-2014 15:02:15"
  
  
  // Remove existing end-of-session state
  DataBaseUtils::DbTransaction transaction( *m_sql );
  try
  {
    Dbo::collection<Dbo::ptr<UserState>> prevEosStates = m_sql->session()->find<UserState>()
      .where( "InterSpecUser_id = ? AND StateType = ?" )
      .bind( m_user.id() )
      .bind( int(UserState::kEndOfSessionTemp) );
    
    //assert( prevEosStates.size() <= 1 );
    for( auto iter = prevEosStates.begin(); iter != prevEosStates.end(); ++iter )
      UserState::removeFromDatabase( *iter, *m_sql );
  }catch( std::exception &e )
  {
    Wt::log("error") << "InterSpec::saveStateForEndOfSession() error removing last state: " << e.what();
    transaction.rollback();
    return;
  }
  
  // I we dont want to save states, or there is no foreground - nothing more to do here
  if( !saveState || !m_dataMeasurement )
  {
    transaction.commit();
    return;
  }
  
  try
  {
    // Check if we are connected with a database state
    const long long int current_state_index = currentAppStateDbId();
    Wt::Dbo::ptr<UserState> parentState;
    if( current_state_index >= 0 )
    {
      parentState = m_sql->session()->find<UserState>( "where id = ?" )
            .bind( current_state_index );
    }
    
    UserState *state = new UserState();
    state->user = m_user;
    state->stateType = UserState::kEndOfSessionTemp;
    state->creationTime = WDateTime::currentDateTime();
    state->name = name;
    state->description = desc;
    state->snapshotTagParent = parentState;
    
    Wt::Dbo::ptr<UserState> dbstate = m_sql->session()->add( state );
    
    saveStateToDb( dbstate );
    
    transaction.commit();
    
    Wt::log("debug") << "Saved end of session app state to database";
  }catch( std::exception &e )
  {
    Wt::log("error") << "InterSpec::saveStateForEndOfSession() error: " << e.what();
    transaction.rollback();
  }
}//void saveStateForEndOfSession()
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
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addFileMenu(): parent passed in must"
                         " be a PopupDivMenu  or WContainerWidget" );

  WString menuname = WString::tr( mobile ? "app-menu-spectra" : (isAppTitlebar ? "app-menu-file" : "interspec") );
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( menuname, menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_fileMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_fileMenuPopup = parentMenu->addPopupMenuItem( menuname );
  }//if( menuDiv ) / else

  PopupDivMenuItem *item = (PopupDivMenuItem *)0;
  
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
    item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-about") );
    item->triggered().connect( this, &InterSpec::showLicenseAndDisclaimersWindow );
    m_fileMenuPopup->addSeparator();  //doesnt seem to be showing up for some reason... owe well.
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif
  
  
#if( USE_DB_TO_STORE_SPECTRA )
  if( m_fileManager )
  {
    if( mobile )
    {
      item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-open") );
#if( ANDROID )
      item->triggered().connect( std::bind([this](){
        doJavaScript( "window.interspecJava.startBrowseToOpenFile();" );
      }) );
#else
      item->triggered().connect( boost::bind ( &SpecMeasManager::startQuickUpload, m_fileManager ) );
#endif
      
      item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-loaded-spec") );
      item->triggered().connect( this, &InterSpec::showCompactFileManagerWindow );
    }//if( mobile )
    
    
    item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-manager"), "InterSpec_resources/images/file_manager_small.png" );
    HelpSystem::attachToolTipOn(item, WString::tr("app-mi-tt-file-manager"), showToolTips );
    
    item->triggered().connect( m_fileManager, &SpecMeasManager::startSpectrumManager );
    
    m_fileMenuPopup->addSeparator();
    
    // --- new save menu ---
    m_saveState = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-store") );
    m_saveState->triggered().connect( boost::bind( &InterSpec::stateSave, this ) );
    HelpSystem::attachToolTipOn(m_saveState, WString::tr("app-mi-tt-file-store"), showToolTips );
    
    
    // --- new save as menu ---
    m_saveStateAs = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-store-as") );
    m_saveStateAs->triggered().connect( boost::bind( &InterSpec::stateSaveAs, this ) );
    HelpSystem::attachToolTipOn(m_saveStateAs, WString::tr("app-mi-tt-file-store-as"), showToolTips );
    m_saveStateAs->setDisabled( true );
    
    // --- new save tag menu ---
    m_createTag = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-tag"), "InterSpec_resources/images/stop_time_small.png");
    m_createTag->triggered().connect( boost::bind(&InterSpec::stateSaveTag, this ));
    
    HelpSystem::attachToolTipOn(m_createTag, WString::tr("app-mi-tt-file-tag"), showToolTips );
    
#if( PERFORM_DEVELOPER_CHECKS || !defined(NDEBUG) )
    // During development the app often gets killed by the IDE, before the state can
    //  be updated - lets add in a way to trigger updating the on-load app state
    item = m_fileMenuPopup->addMenuItem( "Store EoS" );
    HelpSystem::attachToolTipOn(item,
      "For development purposes - Stores current state as your end of session state, "
      "for loading on next launch.", showToolTips );
    item->triggered().connect( this, &InterSpec::saveStateForEndOfSession );
#endif
    
    m_fileMenuPopup->addSeparator();
    
    item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-prev") , "InterSpec_resources/images/db_small.png");
    item->triggered().connect( m_fileManager, &SpecMeasManager::browsePrevSpectraAndStatesDb );
    HelpSystem::attachToolTipOn(item, WString::tr("app-mi-tt-file-prev"), showToolTips );
    
    m_fileMenuPopup->addSeparator();
  }//if( m_fileManager )
#endif  //USE_DB_TO_STORE_SPECTRA


  item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-clear") );
  item->triggered().connect( this, &InterSpec::startClearSession );
  HelpSystem::attachToolTipOn(item, WString::tr("app-mi-tt-file-clear"), showToolTips );
  m_fileMenuPopup->addSeparator();
  
  PopupDivMenu *subPopup = nullptr;
  subPopup = m_fileMenuPopup->addPopupMenuItem( WString::tr("app-mi-samples") );
  
  const string docroot = wApp->docRoot();
  
  item = subPopup->addMenuItem( WString::tr("app-mi-samples-ba133") );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/ba133_source_640s_20100317.n42"),
                                         SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::N42_2006 ) );
  if( mobile )
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-passthrough-mobile") );
  else
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-passthrough") );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/passthrough.n42"),
                                         SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::N42_2006 ) );
  
  item = subPopup->addMenuItem( WString::tr("app-mi-samples-background") );
  item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                         SpecUtils::append_path(docroot, "example_spectra/background_20100317.n42"),
                                         SpecUtils::SpectrumType::Background, SpecUtils::ParserType::N42_2006 ) );
  //If its a mobile device, we'll give a few more spectra to play with
  if( mobile )
  {
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-ba133-lowres") );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Ba133LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-co60-lowres") );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Co60LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-cs137-lowres") );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Cs137LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
    item = subPopup->addMenuItem( WString::tr("app-mi-samples-th232-lowres") );
    item->triggered().connect( boost::bind( &SpecMeasManager::loadFromFileSystem, m_fileManager,
                                           SpecUtils::append_path(docroot, "example_spectra/Th232LowResNoCalib.spe"),
                                           SpecUtils::SpectrumType::Foreground, SpecUtils::ParserType::SpeIaea ) );
  }//if( mobile )
  
  
  if( !mobile )
  { 
    item = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-file-open") );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-file-open"), showToolTips );
    
#if( ANDROID )
    item->triggered().connect( std::bind([this](){
      doJavaScript( "window.interspecJava.startBrowseToOpenFile();" );
    }) );
#else
    item->triggered().connect( boost::bind ( &SpecMeasManager::startQuickUpload, m_fileManager ) );
#endif
  }//if( !mobile )
  
  m_exportSpecFileMenu = m_fileMenuPopup->addMenuItem( WString::tr("app-mi-export-file") );
  m_exportSpecFileMenu->triggered().connect( boost::bind( &InterSpec::createExportSpectrumFileDialog, this ) );
  m_exportSpecFileMenu->disable();
  
#if( USE_DB_TO_STORE_SPECTRA )
  assert( !m_fileManager || m_saveStateAs );
  assert( !m_fileManager || m_saveState );
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  m_fileMenuPopup->addSeparator();
  PopupDivMenu* testing = m_fileMenuPopup->addPopupMenuItem( "Testing" , "InterSpec_resources/images/testing.png");
  item = testing->addMenuItem( "Store App Test State..." );
  item->triggered().connect( boost::bind( &InterSpec::startStoreTestState, this ) );
  HelpSystem::attachToolTipOn(item,"Stores app state as part of the automated test suite", showToolTips );
  
  item = testing->addMenuItem( "Load N42 Test State..." );
  item->triggered().connect( boost::bind(&InterSpec::startN42TestStates, this) );
  HelpSystem::attachToolTipOn(item, "Loads a InterSpec specific variant of a "
                                    "2012 N42 file that contains Foreground, "
                                    "Background, User Settings, and Shielding/"
                                    "Source model.", showToolTips );
  
  item = testing->addMenuItem( "Show Testing Widget..." );
  item->triggered().connect( boost::bind(&InterSpec::startStateTester, this ) );
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE )


#if( USE_OSX_NATIVE_MENU )
    //Add a separtor before Quit InterSpec
    m_fileMenuPopup->addSeparator();
#endif

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
#if( BUILD_AS_ELECTRON_APP )
  m_fileMenuPopup->addSeparator();
  PopupDivMenuItem *exitItem = m_fileMenuPopup->addMenuItem( WString::tr("Exit") );
  exitItem->triggered().connect( std::bind( []{
    ElectronUtils::send_nodejs_message( "CloseWindow", "" );
  }) );
#elif(BUILD_AS_WX_WIDGETS_APP)
    m_fileMenuPopup->addSeparator();
    PopupDivMenuItem* exitItem = m_fileMenuPopup->addMenuItem( WString::tr("Exit") );
    exitItem->triggered().connect(std::bind([] {
      wApp->doJavaScript("window.wx.postMessage('CloseWindow');");
    }));
#endif
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif //  || BUILD_AS_WX_WIDGETS_APP
}//void addFileMenu( WContainerWidget *menuDiv, bool show_examples )


void InterSpec::addEditMenu( Wt::WWidget *parent )
{
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addEditMenu(): parent passed in must"
                         " be a PopupDivMenu or WContainerWidget" );

  const WString menuname = WString::tr( "app-menu-edit" );
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( menuname, menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_editMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_editMenuPopup = parentMenu->addPopupMenuItem( menuname );
  }//if( menuDiv ) / else
  
  int menuindex = -1;
  if( m_undo )
  {
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 0;
#endif
    
    PopupDivMenuItem *undoMenu = m_editMenuPopup->insertMenuItem( menuindex, WString::tr("app-mi-edit-undo"), "", true );
    undoMenu->setDisabled( true );
    undoMenu->triggered().connect( m_undo, &UndoRedoManager::executeUndo );

#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 1;
#endif
    
    PopupDivMenuItem *redoMenu = m_editMenuPopup->insertMenuItem( menuindex, WString::tr("app-mi-edit-redo"), "", true );
    redoMenu->setDisabled( true );
    redoMenu->triggered().connect( m_undo, &UndoRedoManager::executeRedo );
    
    m_undo->undoMenuDisableUpdate().connect( boost::bind(&PopupDivMenuItem::setDisabled, undoMenu, boost::placeholders::_1) );
    m_undo->redoMenuDisableUpdate().connect( boost::bind(&PopupDivMenuItem::setDisabled, redoMenu, boost::placeholders::_1) );
    m_undo->undoMenuToolTipUpdate().connect( boost::bind(&PopupDivMenuItem::setToolTip, undoMenu, boost::placeholders::_1, TextFormat::PlainText) );
    m_undo->redoMenuToolTipUpdate().connect( boost::bind(&PopupDivMenuItem::setToolTip, redoMenu, boost::placeholders::_1, TextFormat::PlainText) );
  }//if( m_undo )
  
#if( !ANDROID )
    // TODO: On Android we dont currently catch the download for these downloads with the content encoded in the URI, so we will just disable these for now there.
  
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 2;
#endif
  m_editMenuPopup->addSeparatorAt( (m_undo ? menuindex : -1) );

#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 3;
#endif
  
  auto saveitem = m_editMenuPopup->insertMenuItem( (m_undo ? 3 : -1), WString::tr("app-mi-edit-save-png"), "", true );
  saveitem->triggered().connect( boost::bind(&InterSpec::saveChartToImg, this, true, true) );
  
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 4;
#endif
  
  saveitem = m_editMenuPopup->insertMenuItem( (m_undo ? 4 : -1), WString::tr("app-mi-edit-save-svg"), "", true );
  saveitem->triggered().connect( boost::bind(&InterSpec::saveChartToImg, this, true, false) );
#endif
  
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 5;
#endif
  m_editMenuPopup->addSeparatorAt( (m_undo ? menuindex : -1) );
  
#if( BUILD_AS_OSX_APP )
  if( InterSpecApp::isPrimaryWindowInstance() )
    menuindex = 6;
#endif
  
  PopupDivMenuItem *uriItem = m_editMenuPopup->insertMenuItem( (m_undo ? menuindex : -1), WString::tr("app-mi-edit-enter-url"), "", true );
  uriItem->triggered().connect( boost::bind(&InterSpec::makeEnterAppUrlWindow, this) );
  
#if( BUILD_AS_OSX_APP )
  // All the macOS native menu stuff (copy/paste/select/etc) will be below our stuff
  if( InterSpecApp::isPrimaryWindowInstance() )
    m_editMenuPopup->addSeparatorAt( (m_undo ? 7 : -1) );
#endif
}//void addEditMenu( Wt::WWidget *menuDiv )


bool InterSpec::toolTabsVisible() const
{
  return m_toolsTabs;
}


#if( InterSpec_PHONE_ROTATE_FOR_TABS )
void InterSpec::toolTabContentHeightChanged( int height )
{
  m_toolsTabsContentHeight = height;
}
#endif


void InterSpec::addToolsTabToMenuItems()
{
  if( !m_toolsMenuPopup )
    return;
  
  if( m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::RefPhotopeaks)] )
    return;
  
  bool showToolTips = false;
  
  try
  {
    showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  }catch(...)
  {
  }
  
  const int index_offest = isMobile() ? 1 : 0;
  
  Wt::WMenuItem *item = nullptr;
  WString tooltip;
  const char *icon = nullptr;
  
  icon = "InterSpec_resources/images/reflines.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 0, WString::tr(GammaLinesTabTitleKey), icon , true );
  item->triggered().connect( boost::bind( &InterSpec::showGammaLinesWindow, this ) );
  
  tooltip = WString::tr("app-tab-tt-ref-photopeak");
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::RefPhotopeaks)] = item;
  
  icon = "InterSpec_resources/images/peakmanager.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 1, WString::tr(PeakInfoTabTitleKey), icon, true );
  item->triggered().connect( this, &InterSpec::showPeakInfoWindow );
  tooltip = WString::tr("app-tab-tt-peak-manager");
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::PeakManager)] = item;
  
  icon = "InterSpec_resources/images/calibrate.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 2, WString::tr(CalibrationTabTitleKey), icon, true );
  item->triggered().connect( this, &InterSpec::showEnergyCalWindow );
  tooltip = WString::tr("app-tab-tt-energy-cal");
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::EnergyCal)] = item;
  
  icon = "InterSpec_resources/images/magnifier_black.svg";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 3, WString::tr(NuclideSearchTabTitleKey), icon, true );
  item->triggered().connect( this, &InterSpec::showNuclideSearchWindow);
  tooltip = WString::tr("app-tab-tt-nuc-search");
  HelpSystem::attachToolTipOn( item, tooltip, showToolTips );
  m_tabToolsMenuItems[static_cast<int>(ToolTabMenuItems::NuclideSearch)] = item;

  icon = "InterSpec_resources/images/auto_peak_search.png";
  item = m_toolsMenuPopup->insertMenuItem( index_offest + 4, WString::tr("app-tab-tt-auto-search"), icon, true );
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
  
  unique_ptr<UndoRedoManager::BlockUndoRedoInserts> undo_sentry;
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this, showToolTabs]{ setToolTabsVisible( !showToolTabs ); };
    auto redo = [this, showToolTabs]{ setToolTabsVisible( showToolTabs ); };
    m_undo->addUndoRedoStep( undo, redo, "Toggle show tool tabs" );
    
    undo_sentry = make_unique<UndoRedoManager::BlockUndoRedoInserts>();
  }//if( m_undo )
  
  // If feature marker window, or column in "Reference Photopeak" tab is open, then we will close
  //  this tool.  We will not open it up at the end of this function, because that seems a little
  //  confusing (and also, it causes a JS exception....)
  const bool showingFeatureMarkers = (m_featureMarkersWindow
                || (m_referencePhotopeakLines && m_referencePhotopeakLines->featureMarkerTool()));
  if( showingFeatureMarkers )
    displayFeatureMarkerWindow( false );
  
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
  
  string refNucXmlState;
  if( m_referencePhotopeakLines
      && ((m_referencePhotopeakLines->currentlyShowingNuclide().m_validity == ReferenceLineInfo::InputValidity::Valid)
           || m_referencePhotopeakLines->persistedNuclides().size()) )
  {
    m_referencePhotopeakLines->serialize( refNucXmlState );
  }
  
  const int layoutVertSpacing = isMobile() ? 10 : 5;

#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    const bool phone = isPhone();
#endif
  
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
    
    closeGammaLinesWindow();
    
    if( m_energyCalWindow )
    {
      m_energyCalWindow->stretcher()->removeWidget( m_energyCalTool );
      delete m_energyCalWindow;
      m_energyCalWindow = nullptr;
    }//if( m_energyCalWindow )
    
    m_toolsTabs = new WTabWidget();
    m_toolsTabs->addStyleClass( "ToolsTabs" );
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    const CompactFileManager::DisplayMode cfmDispMode = phone ? CompactFileManager::Tabbed : CompactFileManager::LeftToRight;
    CompactFileManager *compact = new CompactFileManager( m_fileManager, this, cfmDispMode );
    WString tabTitle = phone ? WString() : WString::tr(FileTabTitleKey);
    WMenuItem *fileTab = m_toolsTabs->addTab( compact, tabTitle, TabLoadPolicy );
    if( phone )
      fileTab->setIcon( "InterSpec_resources/images/spec_files.svg" );
#else
    CompactFileManager *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::LeftToRight );
    m_toolsTabs->addTab( compact, WString::tr(FileTabTitleKey), TabLoadPolicy );
#endif
    
    m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                   boost::placeholders::_1,
                                                   boost::placeholders::_2,
                                                   boost::placeholders::_3 ) );
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    m_peakInfoDisplay->setNarrowPhoneLayout( phone );
    tabTitle = phone ? WString() : WString::tr(PeakInfoTabTitleKey);
    WMenuItem *peakManTab = m_toolsTabs->addTab( m_peakInfoDisplay, tabTitle, TabLoadPolicy );
    if( phone )
      peakManTab->setIcon( "InterSpec_resources/images/peakmanager.svg" ); //
#else
    m_toolsTabs->addTab( m_peakInfoDisplay, WString::tr(PeakInfoTabTitleKey), TabLoadPolicy );
#endif
//    WString tooltip = WString::tr("app-tab-tt-peak-manager");
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
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    m_referencePhotopeakLines->setNarrowPhoneLayout( phone );
    tabTitle = phone ? WString(): WString::tr(GammaLinesTabTitleKey);
    WMenuItem *refPhotoTab = m_toolsTabs->addTab( m_referencePhotopeakLines, tabTitle, TabLoadPolicy );
    if( phone )
      refPhotoTab->setIcon( "InterSpec_resources/images/reflines.svg" );
#else
    m_toolsTabs->addTab( m_referencePhotopeakLines, WString::tr(GammaLinesTabTitleKey), TabLoadPolicy );
#endif
//   WString tooltip = WString::tr("app-tab-tt-ref-photopeak");
//      HelpSystem::attachToolTipOn( refPhotoTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
      
    m_toolsTabs->currentChanged().connect( this, &InterSpec::handleToolTabChanged );
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( phone )
      m_energyCalTool->setTallLayout();
    else
      m_energyCalTool->setWideLayout();
      
    tabTitle = phone ? WString() : WString::tr(CalibrationTabTitleKey);
    WMenuItem *energyCalTab = m_toolsTabs->addTab( m_energyCalTool, tabTitle, TabLoadPolicy );
    if( phone )
      energyCalTab->setIcon( "InterSpec_resources/images/calibrate.svg" );
#else
    m_energyCalTool->setWideLayout();
    m_toolsTabs->addTab( m_energyCalTool, WString::tr(CalibrationTabTitleKey), TabLoadPolicy );
#endif
    
//  WString tooltip = WString::tr("app-tab-tt-energy-cal");
//  HelpSystem::attachToolTipOn( energyCalTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
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
    
#if( OPTIMIZE_D3TimeChart_HIDDEN_LOAD )
    // With out this next line the time chart wont show if we hide tool tabs, and then show
    //  them again.  I think it has something to do with the timechart div getting taken out
    //  of the DOM, and references lost.  Maybe.
    m_timeSeries->refreshJs();
#endif
    
    //Without using the wrapper below, the tabs widget will change height, even
    //  if it is explicitly set, when different tabs are clicked (unless a
    //  manual resize is performed by the user first)
    //m_layout->addLayout( toolsLayout,   ++row, 0 );
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    WContainerWidget *toolsWrapper = new ToolTabContentWrapper( this );
    toolsWrapper->setLayout( m_toolsLayout );
    if( (m_toolsTabsContentHeight <= 100)  //No height set yet (e.g., app is just loading), or user made way to small
       || (m_renderedHeight < 100)         //App is loading
       || ((m_renderedHeight > 0) && (m_toolsTabsContentHeight > (m_renderedHeight/2))) ) //No more than half the screen
    {
      if( phone )
      {
        m_toolsTabsContentHeight = 280;
        if( m_renderedHeight > 100 )
          m_toolsTabsContentHeight = std::max( 245, std::min( m_toolsTabsContentHeight, m_renderedHeight/2 ) );
        else
          m_toolsTabsContentHeight = std::max( 245, std::min( m_toolsTabsContentHeight, wApp->environment().screenHeight()/2 ) );
      }else
      {
        m_toolsTabsContentHeight = 245;
      }
    }
    toolsWrapper->setHeight( m_toolsTabsContentHeight );
#else
    WContainerWidget *toolsWrapper = new WContainerWidget();
    toolsWrapper->setHeight( 245 );
#endif
    
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
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    m_nuclideSearch->setNarrowPhoneLayout( phone );
    tabTitle = phone ? WString() : WString::tr(NuclideSearchTabTitleKey);
    WMenuItem *nuclideTab = m_toolsTabs->addTab( m_nuclideSearchContainer, tabTitle, TabLoadPolicy );
    if( phone )
      nuclideTab->setIcon( "InterSpec_resources/images/magnifier_black.svg" );
#else
    m_toolsTabs->addTab( m_nuclideSearchContainer, WString::tr(NuclideSearchTabTitleKey), TabLoadPolicy );
#endif
//  WString tooltip = WString::tr("app-tab-tt-nuc-search");
//  HelpSystem::attachToolTipOn( nuclideTab, tooltip, showToolTips, HelpSystem::ToolTipPosition::Top );
    
#if( USE_TERMINAL_WIDGET || USE_REL_ACT_TOOL )
    // Handle when the user closes the tab for the Math/Command terminal and the Manual Relative
    //  Activity tool
    m_toolsTabs->tabClosed().connect( boost::bind( &InterSpec::handleToolTabClosed, this, boost::placeholders::_1 ) );
#endif
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( (m_currentToolsTab >= 0) && (m_currentToolsTab < m_toolsTabs->count()) )
    {
      m_toolsTabs->setCurrentIndex( m_currentToolsTab );
    }else
    {
      //These next `setCurrentWidget(...)` lines will cause `handleToolTabChanged(...)` to be called
      if( phone )
        m_toolsTabs->setCurrentWidget( m_referencePhotopeakLines );
      else
        m_toolsTabs->setCurrentWidget( m_peakInfoDisplay );
    }
#else
    //Make sure the current tab is the peak info display
    m_toolsTabs->setCurrentWidget( m_peakInfoDisplay );
#endif
    
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

#if( USE_REL_ACT_TOOL )
    if( m_relActManualGui )
    {
      if( !m_relActManualWindow )
        m_toolsTabs->removeTab( m_energyCalTool );
      handleRelActManualClose();
    }//if( m_relActManualGui )
    
    if( m_relActAutoGui )
      handleRelActAutoClose();
#endif
    
#if( USE_TERMINAL_WIDGET )
    if( m_terminal )
    {
      if( !m_terminalWindow )
        m_toolsTabs->removeTab( m_terminal );
      handleTerminalWindowClose();
    }
#endif
    
    m_nuclideSearch->clearSearchEnergiesOnClient();
    m_nuclideSearchContainer->layout()->removeWidget( m_nuclideSearch );
    m_toolsTabs->removeTab( m_nuclideSearchContainer );
    delete m_nuclideSearchContainer;
    m_nuclideSearchContainer = nullptr;
    
    
    m_toolsTabs = nullptr;
    
    if( m_referencePhotopeakLines )
    {
      m_referencePhotopeakLines->clearAllLines();
      delete m_referencePhotopeakLines;
      m_referencePhotopeakLines = nullptr;
    }
    
    if( m_referencePhotopeakLinesWindow )
      delete m_referencePhotopeakLinesWindow;
    m_referencePhotopeakLinesWindow = nullptr;
    
    if( !refNucXmlState.empty() )
    {
      showGammaLinesWindow();
      if( m_referencePhotopeakLines )
        m_referencePhotopeakLines->deSerialize( refNucXmlState );
      if( phone && m_referencePhotopeakLinesWindow )
        m_referencePhotopeakLinesWindow->hide();
    }//if( !refNucXmlState.empty() )
    
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
  //else
  //  m_currentToolsTab = -1;
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( phone )
    {
      const int w = (m_renderedWidth > 100) ? m_renderedWidth : wApp->environment().screenWidth();
      const int h = (m_renderedWidth > 100) ? m_renderedHeight : wApp->environment().screenHeight();
      if( w < h )
      {
        // On portrait phone, remove the scaler slider, without altering preference
        m_showYAxisScalerItems[0]->show();
        m_showYAxisScalerItems[1]->hide();
        m_spectrum->showYAxisScalers( false );
        
        // TODO: adjust padding in spectrum chart, e.g. in JS:
        //  SpectrumChartD3.padding.labelPad = 5
        //  SpectrumChartD3.padding.left = 5
        //  SpectrumChartD3.padding.right = 10
      }else
      {
        // The phone is landscape, restore showing scaler to user preference
        const bool showScalers = UserPreferences::preferenceValue<bool>( "ShowYAxisScalers", this );
        m_showYAxisScalerItems[0]->setHidden( showScalers );
        m_showYAxisScalerItems[1]->setHidden( !showScalers );
        m_spectrum->showYAxisScalers( showScalers );
      }
    }//if( phone )
#endif
  
  if( m_mobileBackButton && m_mobileForwardButton )
  {
    if( !toolTabsVisible() && m_dataMeasurement && !m_dataMeasurement->passthrough()
        && (m_dataMeasurement->sample_numbers().size() > 1) )
    {
      m_mobileBackButton->setHidden(false);
      m_mobileForwardButton->setHidden(false);
    }else
    {
      m_mobileBackButton->setHidden(true);
      m_mobileForwardButton->setHidden(true);
    }
  }//if( m_mobileBackButton && m_mobileForwardButton )
  
  //Not sure _why_ this next statement is needed, but it is, or else the
  //  spectrum chart shows up with no data
  m_spectrum->scheduleUpdateForeground();
  m_spectrum->scheduleUpdateBackground();
  m_spectrum->scheduleUpdateSecondData();
  
  m_timeSeries->scheduleRenderAll();
#endif // USE_CSS_FLEX_LAYOUT / else
  
  // If we call `displayFeatureMarkerWindow(true);`, we'll get a JS exception - rather than
  //  figure this out, we'll just not re-open it.
  //if( showingFeatureMarkers )
  //  displayFeatureMarkerWindow( true );
}//void setToolTabsVisible( bool showToolTabs )


void InterSpec::addViewMenu( WWidget *parent )
{
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addViewMenu(): parent passed in"
                        " must be a PopupDivMenu  or WContainerWidget" );
 
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( WString::tr("app-menu-view"), menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_displayOptionsPopupDiv = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_displayOptionsPopupDiv = parentMenu->addPopupMenuItem( WString::tr("app-menu-view") );
  }//if( menuDiv ) / else
  
  m_toolTabsVisibleItems[0] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-show-tabs"),
                                                      "InterSpec_resources/images/dock_small.png" );
  m_toolTabsVisibleItems[1] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-hide-tabs"),
                                                    "InterSpec_resources/images/undock_small.png" );
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
  
  PopupDivMenu *chartmenu = m_displayOptionsPopupDiv->addPopupMenuItem( WString::tr("app-mi-view-chart-opts"),
                                              "InterSpec_resources/images/spec_settings_small.png");
  
  bool logypref = true;
  try{ logypref = UserPreferences::preferenceValue<bool>( "LogY", this ); }catch(...){}
  
  m_logYItems[0] = chartmenu->addMenuItem( WString::tr("app-mi-view-logy") );
  m_logYItems[1] = chartmenu->addMenuItem( WString::tr("app-mi-view-liny") );
  m_logYItems[0]->setHidden( logypref );
  m_logYItems[1]->setHidden( !logypref );
  m_logYItems[0]->triggered().connect( boost::bind( &InterSpec::setLogY, this, true  ) );
  m_logYItems[1]->triggered().connect( boost::bind( &InterSpec::setLogY, this, false ) );
  m_spectrum->setYAxisLog( logypref );
  m_spectrum->yAxisLogLinChanged().connect( boost::bind(&InterSpec::setLogY, this, boost::placeholders::_1) );
  m_preferences->addCallbackWhenChanged( "LogY", this, &InterSpec::setLogY );
  
  const bool verticleLines = UserPreferences::preferenceValue<bool>( "ShowVerticalGridlines", this );
  m_verticalLinesItems[0] = chartmenu->addMenuItem( WString::tr("app-mi-view-show-vert"),
                                          "InterSpec_resources/images/sc_togglegridvertical.png");
  m_verticalLinesItems[1] = chartmenu->addMenuItem( WString::tr("app-mi-view-hide-vert"),
                                          "InterSpec_resources/images/sc_togglegridvertical.png");
  m_verticalLinesItems[0]->triggered().connect( boost::bind( &InterSpec::setVerticalLines, this, true ) );
  m_verticalLinesItems[1]->triggered().connect( boost::bind( &InterSpec::setVerticalLines, this, false ) );
  m_verticalLinesItems[0]->setHidden( verticleLines );
  m_verticalLinesItems[1]->setHidden( !verticleLines );
  m_spectrum->showVerticalLines( verticleLines );
  m_timeSeries->showVerticalLines( verticleLines );
  m_preferences->addCallbackWhenChanged( "ShowVerticalGridlines", this, &InterSpec::setVerticalLines );
  
  const bool horizontalLines = UserPreferences::preferenceValue<bool>( "ShowHorizontalGridlines", this );
  m_horizantalLinesItems[0] = chartmenu->addMenuItem( WString::tr("app-mi-view-show-horz"),
                                          "InterSpec_resources/images/sc_togglegridhorizontal.png");
  m_horizantalLinesItems[1] = chartmenu->addMenuItem( WString::tr("app-mi-view-hide-horz"),
                                          "InterSpec_resources/images/sc_togglegridhorizontal.png");
  m_horizantalLinesItems[0]->triggered().connect( boost::bind( &InterSpec::setHorizantalLines, this, true ) );
  m_horizantalLinesItems[1]->triggered().connect( boost::bind( &InterSpec::setHorizantalLines, this, false ) );
  m_horizantalLinesItems[0]->setHidden( horizontalLines );
  m_horizantalLinesItems[1]->setHidden( !horizontalLines );
  m_spectrum->showHorizontalLines( horizontalLines );
  m_timeSeries->showHorizontalLines( horizontalLines );
  m_preferences->addCallbackWhenChanged( "ShowHorizontalGridlines", this, &InterSpec::setHorizantalLines );
  
  
  if( isPhone() )
  {
    m_compactXAxisItems[0] = m_compactXAxisItems[1] = nullptr;
  }else
  {
    const bool makeCompact = UserPreferences::preferenceValue<bool>( "CompactXAxis", this );
    m_spectrum->setCompactAxis( makeCompact );
    m_compactXAxisItems[0] = chartmenu->addMenuItem( WString::tr("app-mi-view-compact-x"), "");
    m_compactXAxisItems[1] = chartmenu->addMenuItem( WString::tr("app-mi-view-normal-x"), "");
    m_compactXAxisItems[0]->triggered().connect( boost::bind( &InterSpec::setXAxisCompact, this, true ) );
    m_compactXAxisItems[1]->triggered().connect( boost::bind( &InterSpec::setXAxisCompact, this, false ) );
    m_compactXAxisItems[0]->setHidden( makeCompact );
    m_compactXAxisItems[1]->setHidden( !makeCompact );
    m_preferences->addCallbackWhenChanged( "CompactXAxis", this, &InterSpec::setXAxisCompact );
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

  PopupDivMenuItem *item = chartmenu->addMenuItem( WString::tr("app-mi-view-show-leg") );
  m_spectrum->legendDisabled().connect( item, &PopupDivMenuItem::show );
  m_spectrum->legendEnabled().connect( item,  &PopupDivMenuItem::hide );
  
  item->triggered().connect( boost::bind( &D3SpectrumDisplayDiv::enableLegend, m_spectrum ) );
  
  item->hide(); //we are already showing the legend
  
  addPeakLabelSubMenu( m_displayOptionsPopupDiv ); //add Peak menu
  
  const bool showSlider = UserPreferences::preferenceValue<bool>( "ShowXAxisSlider", this );
  m_showXAxisSliderItems[0] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-show-slider"), "");
  m_showXAxisSliderItems[1] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-hide-slider"), "");
  m_showXAxisSliderItems[0]->triggered().connect( boost::bind( &InterSpec::setXAxisSlider, this, true, true ) );
  m_showXAxisSliderItems[1]->triggered().connect( boost::bind( &InterSpec::setXAxisSlider, this, false, true ) );
  m_showXAxisSliderItems[0]->setHidden( showSlider );
  m_showXAxisSliderItems[1]->setHidden( !showSlider );
  
  m_spectrum->showXAxisSliderChart( showSlider );
  m_spectrum->xAxisSliderShown().connect( boost::bind(&InterSpec::setXAxisSlider, this, boost::placeholders::_1, true) );
  m_preferences->addCallbackWhenChanged( "ShowXAxisSlider", boost::bind( &InterSpec::setXAxisSlider, this, boost::placeholders::_1, false) );
  
  const bool showScalers = UserPreferences::preferenceValue<bool>( "ShowYAxisScalers", this );
  m_spectrum->showYAxisScalers( showScalers );
  
  m_showYAxisScalerItems[0] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-show-scaler"), "");
  m_showYAxisScalerItems[1] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-hide-scaler"), "");
  m_showYAxisScalerItems[0]->triggered().connect( boost::bind( &InterSpec::setShowYAxisScalers, this, true ) );
  m_showYAxisScalerItems[1]->triggered().connect( boost::bind( &InterSpec::setShowYAxisScalers, this, false ) );
  m_showYAxisScalerItems[0]->setHidden( showScalers );
  m_showYAxisScalerItems[1]->setHidden( !showScalers );
  m_preferences->addCallbackWhenChanged( "ShowYAxisScalers", this, &InterSpec::setShowYAxisScalers );
  
  m_displayOptionsPopupDiv->addSeparator();
  
  m_backgroundSubItems[0] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-back-sub") );
  m_backgroundSubItems[1] = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-back-unsub") );
  m_backgroundSubItems[0]->triggered().connect( boost::bind( &InterSpec::setBackgroundSub, this, true ) );
  m_backgroundSubItems[1]->triggered().connect( boost::bind( &InterSpec::setBackgroundSub, this, false ) );
  m_backgroundSubItems[0]->disable();
  m_backgroundSubItems[1]->hide();
  
  
  m_hardBackgroundSub = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-hard-sub") );
  m_hardBackgroundSub->triggered().connect( this, &InterSpec::startHardBackgroundSub );
  m_hardBackgroundSub->disable();

  
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  m_mapMenuItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-map"),
                                                      "InterSpec_resources/images/map_small.png" );
  m_mapMenuItem->triggered().connect( boost::bind( &InterSpec::createMapWindow, this, 
                                                  SpecUtils::SpectrumType::Foreground, false ) );
  m_mapMenuItem->disable();
  HelpSystem::attachToolTipOn( m_mapMenuItem, WString::tr("app-mi-tt-view-map"), showToolTips );
#endif
  
  /*
   // Example of having a menu-item that triggers opening google maps in an external Window
  auto mapsItem = m_displayOptionsPopupDiv->addMenuItem( "External Google Maps","InterSpec_resources/images/map_small.png" );
  auto dummyItem = m_displayOptionsPopupDiv->addMenuItem( "","" );
  dummyItem->setLink( WLink("#") );
  dummyItem->setLinkTarget( TargetNewWindow );
  
  mapsItem->triggered().connect( std::bind([this,dummyItem](){
    //
    dummyItem->setLink( WLink("https://www.google.com/maps/search/?api=1&query=47.5951518%2C-122.3316393") );
   
    if( dummyItem->anchor() )
    {
      WServer::instance()->schedule( 100, wApp->sessionId(), [dummyItem](){
        const string jsclick = "try{document.getElementById('"+ dummyItem->anchor()->id() + "').click();}catch(e){}";
        wApp->doJavaScript( jsclick );
        wApp->triggerUpdate();
      });
    }else
    {
      cerr << "Unexpected error accessing the anchor for file downloading" << endl;
    }
  }) );
  */
  
  
#if( USE_SEARCH_MODE_3D_CHART )
  m_searchMode3DChart = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-3d"),"" );
  m_searchMode3DChart->triggered().connect( boost::bind( &InterSpec::create3DSearchModeChart, this ) );
  m_searchMode3DChart->disable();
  HelpSystem::attachToolTipOn( m_searchMode3DChart, WString::tr("app-mi-tt-view-3d"), showToolTips );
#endif
  
  m_showRiidResults = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-rid"),"" );
  m_showRiidResults->triggered().connect( boost::bind( &InterSpec::showRiidResults, this, SpecUtils::SpectrumType::Foreground ) );
  HelpSystem::attachToolTipOn( m_showRiidResults, WString::tr("app-mi-tt-view-rid"), showToolTips );
  m_showRiidResults->disable();
  
  m_showMultimedia = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-img"),"" );
  m_showMultimedia->triggered().connect( boost::bind( &InterSpec::showMultimedia, this, SpecUtils::SpectrumType::Foreground ) );
  HelpSystem::attachToolTipOn( m_showMultimedia, WString::tr("app-mi-tt-view-img"), showToolTips );
  m_showMultimedia->disable();
  
  m_featureMarkerMenuItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-feature-markers"), "", true );
  HelpSystem::attachToolTipOn( m_featureMarkerMenuItem, WString::tr("app-mi-tt-view-feature-markers"),
                                showToolTips );
  m_featureMarkerMenuItem->triggered().connect( this, &InterSpec::toggleFeatureMarkerWindow );

  
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
  if (InterSpecApp::isPrimaryWindowInstance())
  {
    m_displayOptionsPopupDiv->addSeparator();

#if( BUILD_AS_ELECTRON_APP )
    auto newWindowItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-new-win") );
    newWindowItem->triggered().connect(std::bind([]() {
      ElectronUtils::send_nodejs_message("NewAppWindow", "");
      }));
#endif

#if( BUILD_AS_WX_WIDGETS_APP )
    auto newWindowItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-new-win") );
    newWindowItem->triggered().connect(std::bind([]() {
      wApp->doJavaScript("window.wx.postMessage('NewAppWindow');");
      }));
#endif

    auto browserItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-webbrowser") );
#if( BUILD_AS_ELECTRON_APP )
    browserItem->triggered().connect(std::bind([]() {
      ElectronUtils::send_nodejs_message("OpenInExternalBrowser", "");
      }));
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
    browserItem->triggered().connect(std::bind([=]() {
      doJavaScript("console.log('Will try to call back to obj-c');"
        "try{"
        "window.webkit.messageHandlers.interOp.postMessage({\"action\": \"ExternalInstance\"});"
        "}catch(error){"
        "console.warn('Failed to callback to the obj-c: ' + error );"
        "}");
      }));
#endif //BUILD_AS_OSX_APP

#if( BUILD_AS_WX_WIDGETS_APP )
    browserItem->triggered().connect(std::bind([=]() {doJavaScript("window.wx.postMessage('OpenInExternalBrowser');"); }));
#endif

  }//if( useNativeMenu )
#endif //BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP
  

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_WX_WIDGETS_APP )
  if (InterSpecApp::isPrimaryWindowInstance())
  {
    m_displayOptionsPopupDiv->addSeparator();
    PopupDivMenuItem *fullScreenItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-full") );  //F11
#if( BUILD_AS_ELECTRON_APP )
    // Note: the triggered() signal a Wt::Signal, which is C++ only, so we cant just hook it up to
    //       javascript for it to run - we have to make the round-trip JS -> C++ -> JS
    fullScreenItem->triggered().connect(std::bind([] {
      ElectronUtils::send_nodejs_message("ToggleMaximizeWindow", "");
    }));
#else
    fullScreenItem->triggered().connect(std::bind([] {
      wApp->doJavaScript("window.wx.postMessage('ToggleMaximizeWindow');");
    }));
#endif
  }//if (InterSpecApp::isPrimaryWindowInstance())

  m_displayOptionsPopupDiv->addSeparator();
  PopupDivMenuItem *resetZoomItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-act-size") ); //Ctrl+0
  PopupDivMenuItem *zoomInItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-zin") ); //Ctrl+Shift+=
  PopupDivMenuItem *zoomOutItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-zout") ); //Ctrl+-

  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsResetPageZoom);
  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsIncreasePageZoom);
  LOAD_JAVASCRIPT(wApp, "js/AppHtmlMenu.js", "AppHtmlMenu", wtjsDecreasePageZoom);
  
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
    PopupDivMenuItem *devToolItem = m_displayOptionsPopupDiv->addMenuItem( WString::tr("app-mi-view-dev-tool") );
    devToolItem->triggered().connect( std::bind( []{
      ElectronUtils::send_nodejs_message( "ToggleDevTools", "" );
    }) );
  }//if( InterSpecApp::isPrimaryWindowInstance() )
#endif
#endif
#endif
}//void addViewMenu( menuParentDiv )


void InterSpec::addDetectorMenu( WWidget *menuWidget )
{
  if( m_detectorToShowMenu )
    return;
  
  PopupDivMenu *parentPopup = dynamic_cast<PopupDivMenu *>( menuWidget );

  //TODO: Change "Detectors" to "Detectors Displayed"
  if( parentPopup )
  {
    m_detectorToShowMenu = parentPopup->addPopupMenuItem( WString::tr("app-mi-view-dets") ); //"InterSpec_resources/images/detector_small_white.png"
  }else
  {
    WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>(menuWidget);
    if( !menuDiv )
      throw runtime_error( "addDetectorMenu: serious error" );
    
    WPushButton *button = new WPushButton( WString::tr("app-mi-view-dets"), menuDiv );
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
      m_toolsTabs->addTab( m_energyCalTool, WString::tr(CalibrationTabTitleKey), TabLoadPolicy );
    
    m_currentToolsTab = m_toolsTabs->currentIndex();
  }//if( m_toolsTabs )
}//void handEnergyCalWindowClose()


EnergyCalTool *InterSpec::energyCalTool()
{
  return m_energyCalTool;
}


UndoRedoManager *InterSpec::undoRedoManager()
{
  return m_undo;
}//UndoRedoManager *undoRedoManager();


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
    
  m_energyCalWindow = new AuxWindow( WString("window-title-energy-cal"),
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
  const bool wasLogY = UserPreferences::preferenceValue<bool>( "LogY", this );
  
  UserPreferences::setPreferenceValue( "LogY", logy, this );
  m_logYItems[0]->setHidden( logy );
  m_logYItems[1]->setHidden( !logy );
  m_spectrum->setYAxisLog( logy );
  
  if( m_undo && (wasLogY != logy) )
  {
    m_undo->addUndoRedoStep( [this,logy](){ setLogY(!logy); },
                            [this,logy](){ setLogY(logy); },
                            "Toggle log-y" );
  }//if( m_undo && (wasLogY != logy) )
}//void setLogY( bool logy )


void InterSpec::setBackgroundSub( bool subtract )
{
  const bool wasBackSub = m_spectrum->backgroundSubtract();
  
  m_backgroundSubItems[0]->setHidden( subtract );
  m_backgroundSubItems[1]->setHidden( !subtract );
  m_spectrum->setBackgroundSubtract( subtract );
  
  if( m_undo && (wasBackSub != subtract) )
  {
    m_undo->addUndoRedoStep( [this,subtract](){ setBackgroundSub(!subtract); },
                            [this,subtract](){ setBackgroundSub(subtract); },
                            "Background subtract" );
  }//if( m_undo && (wasLogY != logy) )
}//void setBackgroundSub( bool sub )


void InterSpec::setVerticalLines( bool show )
{
  const bool wasShow = m_spectrum->verticalLinesShowing();
  
  UserPreferences::setPreferenceValue( "ShowVerticalGridlines", show, this );
  m_verticalLinesItems[0]->setHidden( show );
  m_verticalLinesItems[1]->setHidden( !show );
  m_spectrum->showVerticalLines( show );
  m_timeSeries->showVerticalLines( show );
  
  if( m_undo && (wasShow != show) )
  {
    m_undo->addUndoRedoStep( [this,show](){ setVerticalLines(!show); },
                            [this,show](){ setVerticalLines(show); },
                            "Show vertical lines" );
  }//if( m_undo && (wasLogY != logy) )
}//void setVerticalLines( bool show )


void InterSpec::setHorizantalLines( bool show )
{
  const bool wasShow = m_spectrum->horizontalLinesShowing();
  
  UserPreferences::setPreferenceValue( "ShowHorizontalGridlines", show, this );
  m_horizantalLinesItems[0]->setHidden( show );
  m_horizantalLinesItems[1]->setHidden( !show );
  m_spectrum->showHorizontalLines( show );
  m_timeSeries->showHorizontalLines( show );
  
  if( m_undo && (wasShow != show) )
  {
    m_undo->addUndoRedoStep( [this,show](){ setHorizantalLines(!show); },
                            [this,show](){ setHorizantalLines(show); },
                            "Show horizantal lines" );
  }//if( m_undo && (wasLogY != logy) )
}//void setHorizantalLines( bool show )


void InterSpec::startHardBackgroundSub()
{
  SimpleDialog *dialog = new SimpleDialog( WString::tr("window-title-hard-back-sub"),
                                          WString::tr("window-content-hard-back-sub") );
  
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
  
  WCheckBox *cb = new WCheckBox( WString::tr("window-hard-back-sub-truncate"), optionsDiv );
  cb->setInline( false );
  cb->checked().connect( std::bind([=](){ *truncate_neg = true; } ) );
  cb->unChecked().connect( std::bind([=](){ *truncate_neg = false; } ) );
  
  cb = new WCheckBox( WString::tr("window-hard-back-sub-round"), optionsDiv );
  cb->setInline( false );
  cb->checked().connect( std::bind([=](){ *round_counts = true; } ) );
  cb->unChecked().connect( std::bind([=](){ *round_counts = false; } ) );
  
  
  WPushButton *button = dialog->addButton( WString::tr("Yes") );
  button->setFocus();
  button->clicked().connect( boost::bind( &InterSpec::finishHardBackgroundSub, this, truncate_neg, round_counts ) );
  dialog->addButton( WString::tr("No") );  //dont need to hook this to anything
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
                 WarningWidget::WarningMsgHigh );
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
                WarningWidget::WarningMsgHigh );
  }//try / catch
  
}//finishHardBackgroundSub();



void InterSpec::setXAxisSlider( const bool show, const bool addUndoRedo )
{
  UserPreferences::setPreferenceValue( "ShowXAxisSlider", show, this );
  m_showXAxisSliderItems[0]->setHidden( show );
  m_showXAxisSliderItems[1]->setHidden( !show );
  
  if( show )
  {
    //Default to compact x-axis.
    if( m_compactXAxisItems[0] )
      m_compactXAxisItems[0]->setHidden( true );
    if( m_compactXAxisItems[1] )
      m_compactXAxisItems[1]->setHidden( false );
    
    m_timeSeries->setCompactAxis( true );
  }else
  {
    //Go back to whatever the user wants/selects
    const bool makeCompact = UserPreferences::preferenceValue<bool>( "CompactXAxis", this );
    
    if( m_compactXAxisItems[0] )
      m_compactXAxisItems[0]->setHidden( makeCompact );
    if( m_compactXAxisItems[1] )
      m_compactXAxisItems[1]->setHidden( !makeCompact );
    
    m_spectrum->setCompactAxis( makeCompact ); // This shouldn't be necessary, but JIC
    m_timeSeries->setCompactAxis( makeCompact );
  }//show /hide
  
  m_spectrum->showXAxisSliderChart( show );
  
  
  if( addUndoRedo && m_undo )
  {
    m_undo->addUndoRedoStep( [this,show](){ setXAxisSlider(!show, false); },
                            [this,show](){ setXAxisSlider(show, false); },
                            "Show x-axis slider" );
  }//if( m_undo && (wasLogY != logy) )
}//void setXAxisSlider( bool show )


void InterSpec::setXAxisCompact( bool compact )
{
  const bool wasCompact = m_spectrum->isAxisCompacted();
  
  UserPreferences::setPreferenceValue( "CompactXAxis", compact, this );
  
  if( m_compactXAxisItems[0] )
    m_compactXAxisItems[0]->setHidden( compact );
  if( m_compactXAxisItems[1] )
    m_compactXAxisItems[1]->setHidden( !compact );
  
  m_spectrum->setCompactAxis( compact );
  m_timeSeries->setCompactAxis( compact );
  
  if( m_undo && (wasCompact != compact) )
  {
    m_undo->addUndoRedoStep( [this,compact](){ setXAxisCompact(!compact); },
                            [this,compact](){ setXAxisCompact(compact); },
                            "Set x-axis compact" );
  }//if( m_undo && (wasLogY != logy) )
}//void setXAxisCompact( bool compact )


void InterSpec::setShowYAxisScalers( bool show )
{
  const bool hasSecond = (m_secondDataMeasurement || m_backgroundMeasurement);
  
  assert( m_showYAxisScalerItems[0] );
  assert( m_showYAxisScalerItems[1] );
  
  m_showYAxisScalerItems[0]->setHidden( show );
  m_showYAxisScalerItems[1]->setHidden( !show );
  
  const bool wasShowing = m_spectrum->yAxisScalersIsVisible();
  m_spectrum->showYAxisScalers( show );
  
  try
  {
    UserPreferences::setPreferenceValue( "ShowYAxisScalers", show, this );
  }catch( std::exception &e )
  {
    cerr << "InterSpec::setShowYAxisScalers: Got exception setting pref: " << e.what() << endl;
    assert( 0 );
  }
  
  if( m_undo && (wasShowing != show) )
  {
    m_undo->addUndoRedoStep( [this,show](){ setShowYAxisScalers(!show); },
                            [this,show](){ setShowYAxisScalers(show); },
                            "Show y-axis scalers" );
  }//if( m_undo && (wasLogY != logy) )
}//void setShowYAxisScalers( bool show )



ReferencePhotopeakDisplay *InterSpec::referenceLinesWidget()
{
  return m_referencePhotopeakLines;
}

IsotopeSearchByEnergy *InterSpec::nuclideSearch()
{
  return m_nuclideSearch;
}//IsotopeSearchByEnergy *nuclideSearch();


PeakInfoDisplay *InterSpec::peakInfoDisplay()
{
  return m_peakInfoDisplay;
}//PeakInfoDisplay *peakInfoDisplay();


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
  PopupDivMenu *menu = parentWidget->addPopupMenuItem( WString::tr("app-mi-view-peak-labels"),
                                                        "InterSpec_resources/images/tag.svg" );
  
  // Make a lamda to do work of connecting signals, and undo/redo
  auto setupLabelCbCallbacks = [this]( const SpectrumChart::PeakLabels label, WCheckBox *cb ){
    
    cb->checked().connect( boost::bind( &D3SpectrumDisplayDiv::setShowPeakLabel,
            m_spectrum, label, true
    ) );
    cb->unChecked().connect( boost::bind( &D3SpectrumDisplayDiv::setShowPeakLabel,
            m_spectrum, label, false
    ) );
    
    // FIXME: for some reason undoing things wont change the checkbox state, for first undo, but if you undo/redo/undo it does/  Not totally sure why.
    const auto set_checked = [this,label,cb](){
      cb->setChecked( true );
      // emitting checked() signal, instead of calling `m_spectrum->setShowPeakLabel(label,true);`
      //  to keep `nuc_energy_cb` state consistent.
      cb->checked().emit();
    };
    
    const auto set_unchecked = [this,label,cb](){
      cb->setChecked( false );
      cb->unChecked().emit(); //See comment in `set_checked`
    };
    
    const auto undo_uncheck = [this,set_checked,set_unchecked](){
      if( m_undo )
        m_undo->addUndoRedoStep( set_checked, set_unchecked, "Hide peak label" );
    };
    const auto undo_check = [this,set_checked,set_unchecked](){
      if( m_undo )
        m_undo->addUndoRedoStep( set_unchecked, set_checked, "Show peak label" );
    };
    
    cb->checked().connect( std::bind( undo_check ) );
    cb->unChecked().connect( std::bind( undo_uncheck ) );
  };//auto setupLabelCbCallbacks
  
  WCheckBox *cb = new WCheckBox( WString::tr("app-mi-view-lbl-usr") );
  cb->setChecked(false);
  PopupDivMenuItem *item = menu->addWidget( cb );
  // We will show user labels by default
  item->setChecked( true );
  m_spectrum->setShowPeakLabel( SpectrumChart::PeakLabels::kShowPeakUserLabel, true );
  setupLabelCbCallbacks( SpectrumChart::kShowPeakUserLabel, cb );
  
  cb = new WCheckBox( WString::tr("app-mi-view-lbl-peak-en") );
  cb->setChecked(false);
  item = menu->addWidget( cb );
  setupLabelCbCallbacks( SpectrumChart::kShowPeakEnergyLabel, cb );
  
  cb = new WCheckBox( WString::tr("app-mi-view-lbl-nuc-nm") );
  cb->setChecked(false);
  item = menu->addWidget( cb );
  setupLabelCbCallbacks( SpectrumChart::kShowPeakNuclideLabel, cb );
  
  WCheckBox *nuc_energy_cb = new WCheckBox( WString::tr("app-mi-view-lbl-nuc-en") );
  nuc_energy_cb->setChecked(false);
  item = menu->addWidget( nuc_energy_cb );
  
  nuc_energy_cb->disable();
  cb->checked().connect( nuc_energy_cb, &WCheckBox::enable );
  cb->unChecked().connect( nuc_energy_cb, &WCheckBox::disable );
  cb->unChecked().connect( nuc_energy_cb, &WCheckBox::setUnChecked );
  setupLabelCbCallbacks( SpectrumChart::kShowPeakNuclideEnergies, nuc_energy_cb );
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
    WPushButton *button = new WPushButton( WString::tr("app-menu-help"), menuDiv );
    button->addStyleClass( "MenuLabel" );
    m_helpMenuPopup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    m_helpMenuPopup = parentMenu->addPopupMenuItem( WString::tr("app-menu-help") );
  }//if( menuDiv ) / else
  
  PopupDivMenuItem *item = m_helpMenuPopup->addMenuItem( WString::tr("app-mi-help-welcome") );
  item->triggered().connect( boost::bind( &InterSpec::showWelcomeDialog, this, true ) );

  m_helpMenuPopup->addSeparator();
  
  item = m_helpMenuPopup->addMenuItem( WString::tr("app-mi-help-contents"),
                                      "InterSpec_resources/images/help_minimal.svg");
  
  item->triggered().connect( boost::bind( &HelpSystem::createHelpWindow, string("getting-started") ) );

  Wt::WMenuItem *notifications = m_helpMenuPopup->addMenuItem( WString::tr("app-mi-help-notif-log"),
                                                  "InterSpec_resources/images/log_file_small.png");
  notifications->triggered().connect( this, &InterSpec::showWarningsWindow );

  m_helpMenuPopup->addSeparator();
  PopupDivMenu *subPopup = m_helpMenuPopup->addPopupMenuItem( WString::tr("app-mi-help-opts"),
                                                      "InterSpec_resources/images/cog_small.png" );
    
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  
  // Note: we will associate the WCheckBox's with a option, to set their checked state,
  //       BEFORE adding to the menu - so this way macOS native menus will pickup the
  //       correct values.
  WCheckBox *cb = new WCheckBox( WString::tr("app-mi-help-pref-auto-store") );
  UserPreferences::associateWidget( "AutoSaveSpectraToDb", cb, this );
  item = subPopup->addWidget( cb );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-auto-store"), showToolTips );
  
  
  cb = new WCheckBox( WString::tr("app-mi-help-pref-check-prev") );
  UserPreferences::associateWidget( "CheckForPrevOnSpecLoad", cb, this );
  item = subPopup->addWidget( cb );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-check-prev"), showToolTips );
  
  
  if( !isMobile() )
  {
    WCheckBox *checkbox = new WCheckBox( WString::tr("app-mi-help-pref-show-tt") );
    UserPreferences::associateWidget( "ShowTooltips", checkbox, this );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-show-tt"),
                                true, HelpSystem::ToolTipPosition::Right );
    checkbox->checked().connect( boost::bind( &InterSpec::toggleToolTip, this, true ) );
    checkbox->unChecked().connect( boost::bind( &InterSpec::toggleToolTip, this, false ) );
  }//if( !isMobile() )
  
  {//begin add "AskPropagatePeaks" to menu
    WCheckBox *checkbox = new WCheckBox( WString::tr("app-mi-help-pref-prop-peak") );
    UserPreferences::associateWidget( "AskPropagatePeaks", checkbox, this );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-prop-peak"),
                                 true, HelpSystem::ToolTipPosition::Right );
    checkbox->checked().connect( boost::bind( &InterSpec::toggleToolTip, this, true ) );
    checkbox->unChecked().connect( boost::bind( &InterSpec::toggleToolTip, this, false ) );
  }//end add "AskPropagatePeaks" to menu
  
  
  {//begin add "DisplayBecquerel"
    WCheckBox *checkbox = new WCheckBox( WString::tr("app-mi-help-pref-disp-bq") );
    UserPreferences::associateWidget( "DisplayBecquerel", checkbox, this );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-disp-bq"),
                                 true, HelpSystem::ToolTipPosition::Right );
  }//end add "DisplayBecquerel"
    
  {//begin add "LoadDefaultDrf"
    WCheckBox *checkbox = new WCheckBox( WString::tr("app-mi-help-pref-def-drfs") );
    UserPreferences::associateWidget( "LoadDefaultDrf", checkbox, this );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-def-drfs"),
                                true, HelpSystem::ToolTipPosition::Right );
  }//end add "LoadDefaultDrf"
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  if( app && app->isTablet() )
  {
    WCheckBox *checkbox = new WCheckBox( WString::tr("app-mi-help-pref-desktop") );
    UserPreferences::associateWidget( "TabletUseDesktopMenus", checkbox, this );
    item = subPopup->addWidget( checkbox );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-desktop"),
                                true, HelpSystem::ToolTipPosition::Right );
    
    const WString msg = WString("{1}"
    "<div onclick=\"Wt.emit( $('.specviewer').attr('id'), {name:'miscSignal'}, 'clearSession');"
    "try{$(this.parentElement.parentElement).remove();}catch(e){}return false;\" "
    "class=\"clearsession\"><span class=\"clearsessiontxt\">{2}</span></div>")
      .arg( WString::tr("warn-desktop-disp-switch") )
      .arg( WString::tr("warn-desktop-disp-switch-refresh") );
    
    checkbox->checked().connect( boost::bind( &WarningWidget::addMessageUnsafe,
      m_warnings, msg, WarningWidget::WarningMsgShowOnBoardRiid, 5000 ) );
    checkbox->unChecked().connect( boost::bind( &WarningWidget::addMessageUnsafe,
      m_warnings, msg, WarningWidget::WarningMsgShowOnBoardRiid, 5000 ) );
  }//if( is tablet )
  
  WCheckBox *autoDarkCb = new WCheckBox( WString::tr("app-mi-help-pref-auto-dark") );
  UserPreferences::associateWidget( "AutoDarkFromOs", autoDarkCb, this );
  item = subPopup->addWidget( autoDarkCb );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-auto-dark"),
                              true, HelpSystem::ToolTipPosition::Right );
  
  m_preferences->addCallbackWhenChanged( "AutoDarkFromOs", std::bind([](){
    InterSpec *viewer = InterSpec::instance();
    if( viewer )
      viewer->doJavaScript( "try{ Wt.WT.DetectOsColorThemeJs('" + viewer->id() + "'); }"
                            "catch(e){ console.error('Error with DetectOsColorThemeJs:',e); }" );
  }) );
  
	item = subPopup->addMenuItem( WString::tr("app-mi-help-pref-theme") );
	item->triggered().connect(boost::bind(&InterSpec::showColorThemeWindow, this));

#if( PROMPT_USER_BEFORE_LOADING_PREVIOUS_STATE )
  subPopup->addSeparator();
  WCheckBox *promptOnLoad = new WCheckBox( WString::tr("app-mi-help-pref-prompt-prev") );
  UserPreferences::associateWidget( "PromptStateLoadOnStart", promptOnLoad, this );
  item = subPopup->addWidget( promptOnLoad );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-prompt-prev"), showToolTips );
  
  WCheckBox *doLoad = new WCheckBox( WString::tr("app-mi-help-pref-load-prev") );
  UserPreferences::associateWidget( "LoadPrevStateOnStart", doLoad, this );
  item = subPopup->addWidget( doLoad );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-help-pref-load-prev"), showToolTips );
#endif
  
  
  //Check for multiple languages, and if so, add a menu to let the user select
  const set<string> &languages = InterSpecApp::languagesAvailable();
  if( languages.size() > 1 )
  {
    const string langPref = UserPreferences::preferenceValue<string>("Language", this);
      
    m_languagesSubMenu = m_helpMenuPopup->addPopupMenuItem( WString::tr("app-mi-help-language")  );
    
    item = m_languagesSubMenu->addMenuItem( WString::fromUTF8("Default") );
    if( langPref == "" )
      item->addStyleClass( "CurrentLanguage" );
    item->triggered().connect( boost::bind( &InterSpec::changeLocale, this, "" ) );
    
    
    for( const string &lang : languages )
    {
      string label = lang;
      
      vector<string> parts;
      SpecUtils::split( parts, lang, "-" );
      
      const string &lang_code = !parts.empty() ? parts[0] : lang;
      if( !lang_code.empty() )
      {
        string language;
        if( wApp->messageResourceBundle().resolveKey( "lang-" + lang_code, language ) )
          label = language;
      }//if( !lang_code.empty() )
      
      if( (parts.size() > 1) && !parts[1].empty() )
      {
        string country;
        if( wApp->messageResourceBundle().resolveKey( "country-" + parts[1], country ) )
          label += ", " + country;
      }//if( (parts.size() > 1) && !parts[1].empty() )
      
      item = m_languagesSubMenu->addMenuItem( label );
      if( langPref == lang )
        item->addStyleClass( "CurrentLanguage" );
      item->triggered().connect( boost::bind( &InterSpec::changeLocale, this, lang ) );
    }//for( const string &lang : languages )
  }//if( languages.size() > 1 )
  
#if( BUILD_AS_OSX_APP )
  const bool addAboutInterSpec = !InterSpecApp::isPrimaryWindowInstance();
#else
    const bool addAboutInterSpec = true;
#endif
  
  if( addAboutInterSpec )
  {
    m_helpMenuPopup->addSeparator();
    item = m_helpMenuPopup->addMenuItem( WString::tr("app-mi-help-about") );
    item->triggered().connect( this, &InterSpec::showLicenseAndDisclaimersWindow );
  }

}//void addAboutMenu( Wt::WContainerWidget *menuDiv )


void InterSpec::toggleToolTip( const bool showToolTips )
{
  //update all existing qtips
  if( showToolTips )
  {
    wApp->doJavaScript( "$('.qtip-rounded.canDisableTt').qtip('option', 'show.event', 'mouseenter');" );
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
Wt::Dbo::ptr<UserFileInDb> InterSpec::measurementFromDb( const SpecUtils::SpectrumType type,
                                                          const bool update )
{
  Wt::Dbo::ptr<UserFileInDb> answer;
  
  SpectraFileModel *fileModel = m_fileManager->model();
  std::shared_ptr<SpecMeas> meas = measurment( type );
  if( !meas )
    return answer;
  
  WModelIndex index = fileModel->index( meas );
  if( !index.isValid() )
    return answer;
    
  shared_ptr<SpectraFileHeader> header = fileModel->fileHeader( index.row() );
  if( !header )
    return answer;
  
  answer = header->dbEntry();
  if( answer && !meas->modified() )
    return answer;
    
  if( update )
    header->saveToDatabase( meas );
    
  return header->dbEntry();
}//Wt::Dbo::ptr<UserFileInDb> measurementFromDb( SpecUtils::SpectrumType type, bool update );
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
  
  const int offset = wApp->environment().timeZoneOffset();
  auto localtime = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
  localtime += chrono::seconds(60*offset);
  
  std::string timestr = SpecUtils::to_iso_string( localtime );
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


void InterSpec::setDisplayScaleFactor( const double sf,
                                      const SpecUtils::SpectrumType spec_type,
                                      const bool addUndoRedoStep )
{
  const double prevSF = m_spectrum->displayScaleFactor( spec_type );
  
  m_spectrum->setDisplayScaleFactor( sf, spec_type );
  m_spectrumScaleFactorChanged.emit( spec_type, sf );
  
  if( !addUndoRedoStep || !m_undo || m_undo->isInUndoOrRedo() )
    return;
  
  auto undo = [this, prevSF, spec_type, sf](){
    setDisplayScaleFactor( prevSF, spec_type, false );
    m_spectrum->yAxisScaled().emit( prevSF, sf, spec_type ); // To update the compact file manager
  };
    
  auto redo = [this, sf, spec_type, prevSF](){
    setDisplayScaleFactor( sf, spec_type, false );
    m_spectrum->yAxisScaled().emit( sf, prevSF, spec_type ); // To update the compact file manager
  };
    
  m_undo->addUndoRedoStep( undo, redo, "Change y-axis scale factor" );
}//void setDisplayScaleFactor( const double sf, SpecUtils::SpectrumType spectrum_type );


void InterSpec::handleDisplayScaleFactorChangeFromSpectrum( const double sf, const double prevSF,
                                                           const SpecUtils::SpectrumType spec_type )
{
  // This function is called when the user slides the slider on the spectrum, through the
  //  D3SpectrumDisplayDiv::yAxisScaled() signal.
  
  if( !m_undo || m_undo->isInUndoOrRedo() )
    return;
  
  auto undo = [this, prevSF, spec_type, sf](){
    setDisplayScaleFactor( prevSF, spec_type, false );
    m_spectrum->yAxisScaled().emit( prevSF, sf, spec_type ); // To update the compact file manager
  };
  
  auto redo = [this, sf, spec_type, prevSF](){
    setDisplayScaleFactor( sf, spec_type, false );
    m_spectrum->yAxisScaled().emit( sf, prevSF, spec_type ); // To update the compact file manager
  };
  
  m_undo->addUndoRedoStep( undo, redo, "Change y-axis scale factor" );
}//void handleDisplayScaleFactorChangeFromSpectrum(...)


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


OneOverR2Calc *InterSpec::createOneOverR2Calculator()
{
  if( !m_1overR2Calc )
  {
    m_1overR2Calc = new OneOverR2Calc();
    m_1overR2Calc->finished().connect( boost::bind( &InterSpec::deleteOneOverR2Calc, this ) );
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){ deleteOneOverR2Calc(); };
      auto redo = [this](){ createOneOverR2Calculator(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show 1/r2 calculator" );
    }
  }//if( !m_1overR2Calc )
  
  m_1overR2Calc->show();
  m_1overR2Calc->resizeToFitOnScreen();
  m_1overR2Calc->centerWindowHeavyHanded();
  
  return m_1overR2Calc;
}//void createOneOverR2Calculator()


void InterSpec::deleteOneOverR2Calc()
{
  if( !m_1overR2Calc )
    return;
  
  const bool do_undo = (m_undo && m_undo->canAddUndoRedoNow());
  const string state_uri = do_undo ? m_1overR2Calc->encodeStateToUrl() : string();
  
  AuxWindow::deleteAuxWindow( m_1overR2Calc );
  m_1overR2Calc = nullptr;
  
  if( do_undo )
  {
    auto undo = [this,state_uri](){
      OneOverR2Calc *calc = createOneOverR2Calculator();
      if( calc )
        calc->handleAppUrl( state_uri );
    };
    auto redo = [this](){ deleteOneOverR2Calc(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Hide 1/r2 calculator" );
  }//if( do_undo )
}//void deleteOneOverR2Calc()


UnitsConverterTool *InterSpec::createUnitsConverterTool()
{
  if( !m_unitsConverter )
  {
    m_unitsConverter = new UnitsConverterTool();
    m_unitsConverter->finished().connect( boost::bind( &InterSpec::deleteUnitsConverterTool, this ) );
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){ deleteUnitsConverterTool(); };
      auto redo = [this](){ createUnitsConverterTool(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show units converter" );
    }//if( undo )
  }//if( !m_unitsConverter )
  
  m_unitsConverter->show();
  m_unitsConverter->resizeToFitOnScreen();
  m_unitsConverter->centerWindowHeavyHanded();
  
  return m_unitsConverter;
}//void createUnitsConverterTool()


void InterSpec::deleteUnitsConverterTool()
{
  if( !m_unitsConverter )
    return;
  
  const bool do_undo = (m_undo && m_undo->canAddUndoRedoNow());
  const string state_uri = do_undo ? m_unitsConverter->encodeStateToUrl() : string();
  
  AuxWindow::deleteAuxWindow( m_unitsConverter );
  m_unitsConverter = nullptr;
  
  if( do_undo )
  {
    auto undo = [this,state_uri](){
      UnitsConverterTool *tool = createUnitsConverterTool();
      if( tool )
        tool->handleAppUrl( state_uri );
    };
    auto redo = [this](){ deleteUnitsConverterTool(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Hide units converter" );
  }//if( do_undo )
}//void deleteUnitsConverterTool()


FluxToolWindow *InterSpec::createFluxTool()
{
  if( !m_fluxTool )
  {
    m_fluxTool = new FluxToolWindow( this );
    m_fluxTool->finished().connect( boost::bind( &InterSpec::deleteFluxTool, this ) );
  }
  
  m_fluxTool->show();
  m_fluxTool->resizeToFitOnScreen();
  m_fluxTool->centerWindowHeavyHanded();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ deleteFluxTool(); };
    auto redo = [this](){ createFluxTool(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show flux tool" );
  }//if( undo )
  
  return m_fluxTool;
}//void createFluxTool()


void InterSpec::deleteFluxTool()
{
  if( !m_fluxTool )
    return;
  
  const bool do_undo = (m_undo && m_undo->canAddUndoRedoNow());
  const string state_uri = do_undo ? m_fluxTool->encodeStateToUrl() : string();
  
  AuxWindow::deleteAuxWindow( m_fluxTool );
  m_fluxTool = nullptr;
  
  if( do_undo )
  {
    auto undo = [this,state_uri](){
      FluxToolWindow *flux = createFluxTool();
      if( flux )
        flux->handleAppUrl( state_uri );
    };
    auto redo = [this](){ deleteFluxTool(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Hide flux tool" );
  }//if( do_undo )
}//void deleteFluxTool();


DecayWindow *InterSpec::createDecayInfoWindow()
{
  string initial_uri;
  
  if( m_decayInfoWindow )
  {
    initial_uri = m_decayInfoWindow->encodeStateToUrl();
  }else
  {
    m_decayInfoWindow = new DecayWindow( this );
    m_decayInfoWindow->finished().connect( boost::bind( &InterSpec::deleteDecayInfoWindow, this ) );
  }
  
  if( m_referencePhotopeakLines )
  {
    const ReferenceLineInfo &nuc = m_referencePhotopeakLines->currentlyShowingNuclide();
    if( nuc.m_nuclide )
    {
      m_decayInfoWindow->clearAllNuclides();
      
      // TODO: We could do a little better and check the Shielding/Source Fit widget and grab those activities (and ages) if they match
      
      const bool useBq = UserPreferences::preferenceValue<bool>("DisplayBecquerel", InterSpec::instance());
      const double act = useBq ? PhysicalUnits::MBq : (1.0E-6 * PhysicalUnits::curie);
      const string actStr = useBq ? "1 MBq" : "1 uCi";
      const double hl = nuc.m_nuclide->halfLife;
      double age = hl;
      try
      {
        age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( nuc.m_input.m_age, hl );
      }catch(exception &)
      {}
      
      m_decayInfoWindow->addNuclide( nuc.m_nuclide->atomicNumber,
                         nuc.m_nuclide->massNumber,
                         nuc.m_nuclide->isomerNumber,
                         1.0*PhysicalUnits::microCi, !useBq,
                         0.0, actStr, 5.0*age );
    }//if( nuc.nuclide )
  }//if( m_referencePhotopeakLines )
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this,initial_uri](){
      if( initial_uri.empty() )
      {
        deleteDecayInfoWindow();
      }else
      {
        DecayWindow *dialog = createDecayInfoWindow();
        if( dialog )
          dialog->handleAppUrl( initial_uri );
      }
    };
    auto redo = [this](){ createDecayInfoWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show decay tool" );
  }//if( undo )
  
  return m_decayInfoWindow;
}//void createDecayInfoWindow()


MakeFwhmForDrfWindow *InterSpec::fwhmFromForegroundWindow( const bool use_auto_fit_peaks )
{
  if( m_addFwhmTool )
    return m_addFwhmTool;
  
  m_addFwhmTool = new MakeFwhmForDrfWindow( use_auto_fit_peaks );
  m_addFwhmTool->tool()->updatedDrf().connect( m_addFwhmTool, &AuxWindow::hide );
  m_addFwhmTool->finished().connect( boost::bind( &InterSpec::deleteFwhmFromForegroundWindow, this ) );
  
  if( m_drfSelectWindow )
  {
    m_addFwhmTool->tool()->updatedDrf().connect(
        boost::bind( &DrfSelect::handleFitFwhmFinished,
                    m_drfSelectWindow->widget(), boost::placeholders::_1
    ) );
  }//if( m_drfSelectWindow )
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ deleteFwhmFromForegroundWindow(); };
    auto redo = [this,use_auto_fit_peaks](){ fwhmFromForegroundWindow(use_auto_fit_peaks); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show fit FWHM from foreground tool" );
  }//if( undo )
  
  return m_addFwhmTool;
}//MakeFwhmForDrfWindow *fwhmFromForegroundWindow()


void InterSpec::deleteFwhmFromForegroundWindow()
{
  if( !m_addFwhmTool )
    return;
  
  shared_ptr<MakeFwhmForDrf::ToolState> state = m_addFwhmTool->tool()->currentState();
  
  AuxWindow::deleteAuxWindow( m_addFwhmTool );
  m_addFwhmTool = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this,state](){
      MakeFwhmForDrfWindow *window = fwhmFromForegroundWindow(false);
      if( window )
        window->tool()->setState( state );
      
      shared_ptr<DetectorPeakResponse> drf = m_dataMeasurement ? m_dataMeasurement->detector() : nullptr;
      if( state->m_orig_drf != drf )
      {
        m_detectorChanged.emit(state->m_orig_drf);
      }
    };
    auto redo = [this](){ deleteFwhmFromForegroundWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close fit FWHM from foreground tool" );
  }//if( do_undo )
}//void deleteFwhmFromForegroundWindow()


#if( USE_DETECTION_LIMIT_TOOL )
void InterSpec::showDetectionLimitTool( const std::string &query_str )
{
  if( m_detectionLimitWindow )
    return;
    
  auto tool = createDetectionLimitTool();
  if( tool && !query_str.empty() )
    tool->tool()->handleAppUrl( query_str );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ programmaticallyCloseDetectionLimit(); };
    auto redo = [this,query_str](){ showDetectionLimitTool(query_str); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show detection limit tool" );
  }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
}//void InterSpec::showDetectionLimitTool()


DetectionLimitWindow *InterSpec::createDetectionLimitTool()
{
  if( !m_detectionLimitWindow )
  {
    m_detectionLimitWindow = new DetectionLimitWindow( this, m_materialDB.get(), m_shieldingSuggestion );
    m_detectionLimitWindow->finished().connect( this, &InterSpec::handleDetectionLimitWindowClose );
  }//if( !m_detectionLimitWindow )
  
  return m_detectionLimitWindow;
}//DetectionLimitWindow *createDetectionLimitTool()

     
void InterSpec::handleDetectionLimitWindowClose()
{
  auto *caller = dynamic_cast<DetectionLimitWindow *>( WObject::sender() );
  assert( caller );
  assert( !m_detectionLimitWindow || (caller == m_detectionLimitWindow) );
  
  if( !m_detectionLimitWindow && !caller )
    return;
  
  if( m_detectionLimitWindow && caller && (m_detectionLimitWindow != caller) )
    return;
    
  auto dialog = m_detectionLimitWindow ? m_detectionLimitWindow : caller;
  assert( caller );
  
  if( !m_detectionLimitWindow )
  {
    AuxWindow::deleteAuxWindow( dialog );
    return;
  }
  
  m_detectionLimitWindow = nullptr;
  
  const bool canUndo = (m_undo && m_undo->canAddUndoRedoNow());
  
  string state_uri;
  if( canUndo )
    state_uri = dialog->tool()->encodeStateToUrl();
  
  AuxWindow::deleteAuxWindow( dialog );
  
  if( canUndo )
  {
    // We'll just assume results for the foreground was showing, so we dont have to pass this around
    auto undo = [this,state_uri](){
      DetectionLimitWindow *tool = createDetectionLimitTool();
      assert( tool );
      if( tool && tool->tool() )
        tool->tool()->handleAppUrl(state_uri);
    };
    auto redo = [this](){ programmaticallyCloseSimpleMda(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Detection Limit Tool" );
  }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
}//void handleDetectionLimitWindowClose()


void InterSpec::programmaticallyCloseDetectionLimit()
{
  if( !m_detectionLimitWindow )
    return;
  
  DetectionLimitWindow *dialog = m_detectionLimitWindow;
  m_detectionLimitWindow = nullptr; //Prevents undo/redo
  dialog->hide();
}//void programmaticallyCloseDetectionLimit()


DetectionLimitSimpleWindow *InterSpec::showSimpleMdaWindow()
{
  if( m_simpleMdaWindow )
    return m_simpleMdaWindow;
  
  m_simpleMdaWindow = new DetectionLimitSimpleWindow( m_materialDB.get(), m_shieldingSuggestion, this );
  m_simpleMdaWindow->finished().connect( this, &InterSpec::handleSimpleMdaWindowClose );
  
  return m_simpleMdaWindow;
}//DetectionLimitSimpleWindow *showSimpleMdaWindow()


void InterSpec::programmaticallyCloseSimpleMda()
{
  if( !m_simpleMdaWindow )
    return;
  
  DetectionLimitSimpleWindow *dialog = m_simpleMdaWindow;
  m_simpleMdaWindow = nullptr; //Prevents undo/redo
  dialog->done( WDialog::DialogCode::Accepted );
}//void programmaticallyCloseSimpleMda()


void InterSpec::handleSimpleMdaWindowClose()
{
  auto *caller = dynamic_cast<DetectionLimitSimpleWindow *>( WObject::sender() );
  assert( caller );
  assert( !m_simpleMdaWindow || (caller == m_simpleMdaWindow) );
  
  if( !m_simpleMdaWindow )
    return;
  
  auto dialog = m_simpleMdaWindow;
  m_simpleMdaWindow = nullptr;
  
  const string state_uri = dialog->tool()->encodeStateToUrl();
  
  AuxWindow::deleteAuxWindow( dialog );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    // We'll just assume results for the foreground was showing, so we dont have to pass this around
    auto undo = [this,state_uri](){
      DetectionLimitSimpleWindow *tool = showSimpleMdaWindow();
      assert( tool );
      if( tool && tool->tool() )
        tool->tool()->handleAppUrl(state_uri);
    };
    auto redo = [this](){ programmaticallyCloseSimpleMda(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Simple MDA" );
  }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
}//void handleSimpleMdaWindowClose()


void InterSpec::fitNewPeakNotInRoiFromRightClick()
{
  searchForSinglePeak( m_rightClickEnergy );
}//void fitNewPeakNotInRoiFromRightClick()


void InterSpec::startAddPeakFromRightClick()
{
  // TODO: add AddNewPeakDialog pointer to InterSpec class, like other tools, to fully support undo/redo, and everything.
  const double energy = m_rightClickEnergy;
  
  AddNewPeakDialog *window = new AddNewPeakDialog( energy );
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto redo = [energy](){
      AddNewPeakDialog *window = new AddNewPeakDialog( energy );
      window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    };
    
    auto closer = wApp->bind( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    m_undo->addUndoRedoStep( closer, redo , "Start add peak dialog." );
  }//if( undo )
}//void startAddPeakFromRightClick()


void InterSpec::searchOnEnergyFromRightClick()
{
  if( m_toolsTabs )
  {
    const int prevTab = m_toolsTabs->currentIndex();
    
    assert( m_nuclideSearchContainer );
    m_toolsTabs->setCurrentWidget( m_nuclideSearchContainer );
    
    const int nowTab = m_toolsTabs->currentIndex();
    
    if( prevTab != nowTab )
    {
      if( m_undo && m_undo->canAddUndoRedoNow() )
      {
        auto undo = [this,prevTab](){ m_toolsTabs->setCurrentIndex(prevTab); };
        auto redo = [this](){ m_toolsTabs->setCurrentWidget( m_nuclideSearchContainer ); };
        m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Change to search tab." );
      }//if( m_undo && m_undo->canAddUndoRedoNow() )
    }
  }else
  {
    const bool alreadyShowingWindow = (m_nuclideSearchWindow && m_nuclideSearchWindow->isVisible());
    
    showNuclideSearchWindow();
    
    if( !alreadyShowingWindow && m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){
        if( m_nuclideSearchWindow )
          m_nuclideSearchWindow->hide();
      };
      auto redo = [this](){ showNuclideSearchWindow(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show energy search tool." );
    }
  }//if( m_toolsTabs ) / else
  
  assert( m_nuclideSearch );
  
  vector<IsotopeSearchByEnergy::SearchEnergy *> orig_searchs = m_nuclideSearch->searches();
  for( size_t i = 1; i < orig_searchs.size(); ++i )
    m_nuclideSearch->removeSearchEnergy( orig_searchs[i] );
  
  m_nuclideSearch->setNextSearchEnergy( m_rightClickEnergy ); //Adds its own undo/redo step
}//void searchOnEnergyFromRightClick()


void InterSpec::startSimpleMdaFromRightClick()
{
  string prevState;
  bool wasShowing = false;
  
  if( m_simpleMdaWindow )
  {
    wasShowing = true;
    prevState = m_simpleMdaWindow->tool()->encodeStateToUrl();
  }else
  {
    showSimpleMdaWindow();
  }
  
  assert( m_simpleMdaWindow );
  if( !m_simpleMdaWindow )
    return;
  
  // Check if showing ref photopeak lines for a nuclide
  // Check if near-enough reference photopeak line
  // set nuclide and energy to tool
  if( m_referencePhotopeakLines )
  {
    tuple<const SandiaDecay::Nuclide *, double, float> line
                = PeakSearchGuiUtils::nuclide_reference_line_near( this, m_rightClickEnergy );
    
    const SandiaDecay::Nuclide *nuc = get<0>(line);
    const double age = get<1>(line);
    const float energy = get<2>(line);
    
    if( energy > 10.0f )
      m_simpleMdaWindow->tool()->setNuclide( nuc, age, energy );
    else
      m_simpleMdaWindow->tool()->setNuclide( nullptr, 0.0, m_rightClickEnergy );
  }else
  {
    m_simpleMdaWindow->tool()->setNuclide( nullptr, 0.0, m_rightClickEnergy );
  }//if( m_referencePhotopeakLines )
  
  
  const string currentState = m_simpleMdaWindow->tool()->encodeStateToUrl();
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [=](){
      if( wasShowing )
      {
        if( prevState.empty() )
          programmaticallyCloseSimpleMda();
        showSimpleMdaWindow();
        if( m_simpleMdaWindow && !prevState.empty() )
          m_simpleMdaWindow->tool()->handleAppUrl( prevState );
      }else
      {
        programmaticallyCloseSimpleMda();
      }
    };//undo
    
    auto redo = [=](){
      showSimpleMdaWindow();
      assert( m_simpleMdaWindow );
      if( m_simpleMdaWindow )
        m_simpleMdaWindow->tool()->handleAppUrl( currentState );
    };//redo
    
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Simple MDA Tool." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void startSimpleMdaFromRightClick()
#endif //USE_DETECTION_LIMIT_TOOL


void InterSpec::deleteDecayInfoWindow()
{
  if( m_decayInfoWindow && m_undo && m_undo->canAddUndoRedoNow() )
  {
    const string initial_uri = m_decayInfoWindow->encodeStateToUrl();
    
    auto undo = [this,initial_uri](){
      DecayWindow *dialog = createDecayInfoWindow();
      if( dialog )
        dialog->handleAppUrl( initial_uri );
    };
    auto redo = [this](){ deleteDecayInfoWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close decay tool" );
  }//if( m_decayInfoWindow && m_undo && m_undo->canAddUndoRedoNow() )
  
  
  if( m_decayInfoWindow )
    AuxWindow::deleteAuxWindow( m_decayInfoWindow );
  m_decayInfoWindow = nullptr;
}//void deleteDecayInfoWindow()


void InterSpec::createFileParameterWindow( const SpecUtils::SpectrumType type )
{
  SpecFileSummary *window = new SpecFileSummary( type, this );
  new UndoRedoManager::BlockGuiUndoRedo( window ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
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
    passMessage( "Could not load spectrum", WarningWidget::WarningMsgHigh );
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
                 WarningWidget::WarningMsgHigh );
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

#elif( USE_LEAFLET_MAP )

void InterSpec::createMapWindow( const SpecUtils::SpectrumType spectrum_type, const bool noWarning )
{
  // We will display for all sample numbers, but only include the displayed detectors
  
  const shared_ptr<const SpecMeas> meas = measurment( spectrum_type );
  const vector<string> detectors = detectorsToDisplay( spectrum_type );
  
  if( m_leafletWarning )
    programmaticallyCloseLeafletWarning();
  
  if( m_leafletWindow )
    programmaticallyCloseLeafletMap();
  
  auto onMapCreate = boost::bind( &InterSpec::handleLeafletMapOpen, this, boost::placeholders::_1 );
  
  m_leafletWarning = LeafletRadMap::showForMeasurement( meas, {}, detectors, onMapCreate, noWarning );
    
  if( !m_leafletWarning )
    return;
    
  m_leafletWarning->finished().connect( this, &InterSpec::handleLeafletWarningWindowClose );
    
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    // Note: If the user clicks cancel on the warning screen, hitting undo will then be a blank
    //  step, and it wont bring back up the warning dialog.  We could add another callback to
    //  #LeafletRadMap::showForMeasurement, or have the above case be called with just a nullptr,
    //  but this all is adding more complexity than its worth.
    
    auto undo = [this](){ programmaticallyCloseLeafletWarning(); };
    auto redo = [this, spectrum_type, noWarning](){ createMapWindow( spectrum_type, noWarning ); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show map." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void createMapWindow( SpecUtils::SpectrumType spectrum_type )


void InterSpec::handleLeafletWarningWindowClose()
{
  WObject *signalSender = WObject::sender();
  SimpleDialog *dialogSender = dynamic_cast<SimpleDialog *>( signalSender );
  assert( dialogSender );
  assert( !m_leafletWarning || (dialogSender == m_leafletWarning) );
  
  if( dialogSender == m_leafletWarning )
    m_leafletWarning = nullptr;
}//void handleLeafletWarningWindowClose()


void InterSpec::handleLeafletMapOpen( LeafletRadMapWindow *window )
{
  if( m_leafletWindow )
  {
    AuxWindow::deleteAuxWindow( m_leafletWindow );
    m_leafletWindow = nullptr;
  }
  
  m_leafletWindow = window;
  if( !m_leafletWindow )
    return;
  
  m_leafletWindow->finished().connect( this, &InterSpec::handleLeafletMapClosed );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    using SpecUtils::SpectrumType;
    const shared_ptr<const SpecMeas> meas = m_leafletWindow->map()->measurement();
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    for( auto i : {SpectrumType::Foreground, SpectrumType::Background, SpectrumType::SecondForeground} )
    {
      if( meas == measurment(i) )
      {
        type = i;
        break;
      }
    }
    
    auto undo = [this](){ programmaticallyCloseLeafletMap(); };
    auto redo = [this, type](){ createMapWindow( type, true ); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show map." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void handleLeafletMapOpen()


void InterSpec::programmaticallyCloseLeafletMap()
{
  if( !m_leafletWindow )
    return;
  
  LeafletRadMapWindow *dialog = m_leafletWindow;
  m_leafletWindow = nullptr;
  AuxWindow::deleteAuxWindow( dialog );
}//void programmaticallyCloseLeafletMap();


void InterSpec::programmaticallyCloseLeafletWarning()
{
  if( !m_leafletWarning )
    return;
  
  SimpleDialog *dialog = m_leafletWarning;
  m_leafletWarning = nullptr;
  dialog->done( WDialog::DialogCode::Accepted );
}//void programmaticallyCloseLeafletWarning();


void InterSpec::handleLeafletMapClosed()
{
  WObject *signalSender = WObject::sender();
  LeafletRadMapWindow *mapSender = dynamic_cast<LeafletRadMapWindow *>( signalSender );
  
  assert( mapSender );
  assert( !m_leafletWindow || (m_leafletWindow == mapSender) );
  
  if( m_leafletWindow && (mapSender == m_leafletWindow) && m_undo && m_undo->canAddUndoRedoNow() )
  {
    using SpecUtils::SpectrumType;
    const shared_ptr<const SpecMeas> meas = m_leafletWindow->map()->measurement();
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    for( auto i : {SpectrumType::Foreground, SpectrumType::Background, SpectrumType::SecondForeground} )
    {
      if( meas == measurment(i) )
      {
        type = i;
        break;
      }
    }
    
    auto undo = [this,type](){ createMapWindow(type, true); };
    auto redo = [this](){ programmaticallyCloseLeafletMap(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show map." );
  }//if( we want to add undo/redo step )
  
  m_leafletWindow = nullptr;
}//void handleLeafletMapClosed();

#endif //#if( USE_GOOGLE_MAP ) / #elif( USE_LEAFLET_MAP )


#if( USE_SEARCH_MODE_3D_CHART )
void InterSpec::create3DSearchModeChart()
{
  if( !m_dataMeasurement || !m_dataMeasurement->passthrough() )
  {
    passMessage( "The 3D chart is only available for search mode or RPM passthrough data.",
                WarningWidget::WarningMsgInfo );
    return;
  }//if( we dont have the proper data to make a 3D chart )
  
  if( m_3dViewWindow )
    programmaticallyClose3DSearchModeChart();
  
  m_3dViewWindow = new AuxWindow( WString::tr("window-title-3d"),
                                 (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                  | AuxWindowProperties::EnableResize
                                  | AuxWindowProperties::TabletNotFullScreen) );
  //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
  //  size to be the min-size
  m_3dViewWindow->setMinimumSize( 512, 420 );
  
  m_3dViewWindow->disableCollapse();
  m_3dViewWindow->Wt::WDialog::rejectWhenEscapePressed();
  
  WGridLayout *layout = m_3dViewWindow->stretcher();
  layout->setHorizontalSpacing( 0 );
  layout->setVerticalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  SearchMode3DChart *chart = new SearchMode3DChart( this );
  layout->addWidget( chart, 0, 0 );
  
  Wt::WPushButton *closeButton = m_3dViewWindow->addCloseButtonToFooter( WString::tr("Close"),true);
  closeButton->clicked().connect( m_3dViewWindow, &AuxWindow::hide );
  
  m_3dViewWindow->show();
  m_3dViewWindow->setClosable( true );
  m_3dViewWindow->resizeScaledWindow( 0.7, 0.9 );
  m_3dViewWindow->centerWindow();
  m_3dViewWindow->setResizable( true );
  m_3dViewWindow->finished().connect( boost::bind( &InterSpec::handle3DSearchModeChartClose, this, m_3dViewWindow ) );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ programmaticallyClose3DSearchModeChart(); };
    auto redo = [this](){ create3DSearchModeChart(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show 3D view" );
  }//
}//void create3DSearchModeChart()


void InterSpec::programmaticallyClose3DSearchModeChart()
{
  if( !m_3dViewWindow )
    return;
  
  AuxWindow *dialog = m_3dViewWindow;
  m_3dViewWindow = nullptr;
  AuxWindow::deleteAuxWindow( dialog );
}//void programmaticallyClose3DSearchModeChart()


void InterSpec::handle3DSearchModeChartClose( AuxWindow *window )
{
  assert( window == m_3dViewWindow );
  if( (window == m_3dViewWindow) && window )
  {
    m_3dViewWindow = nullptr;
    AuxWindow::deleteAuxWindow( window );
  }
  
  m_3dViewWindow = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ create3DSearchModeChart(); };
    auto redo = [this](){ programmaticallyClose3DSearchModeChart(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close 3D view" );
  }//
}//void handle3DSearchModeChartClose()
#endif


void InterSpec::showRiidResults( const SpecUtils::SpectrumType type )
{
  if( m_riidDisplay )
    programmaticallyCloseRiidResults();
  
  m_riidDisplay = showRiidInstrumentsAna( measurment(type) );
  m_riidDisplay->finished().connect( this, &InterSpec::handleRiidResultsClose );
  
  if( m_riidDisplay && m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ programmaticallyCloseRiidResults(); };
    auto redo = [this,type](){ showRiidResults(type); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show ID Results" );
  }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
}//void showRiidResults( const SpecUtils::SpectrumType type )


void InterSpec::handleRiidResultsClose()
{
  SimpleDialog *caller = dynamic_cast<SimpleDialog *>( WObject::sender() );
  
  if( !m_riidDisplay || (caller != m_riidDisplay) )
    return;
  
  m_riidDisplay = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    // We'll just assume results for the foreground was showing, so we dont have to pass this around
    auto undo = [this](){ showRiidResults(SpecUtils::SpectrumType::Foreground); };
    auto redo = [this](){ programmaticallyCloseRiidResults(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close ID Results" );
  }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
}//void handleRiidResultsClose()


void InterSpec::programmaticallyCloseRiidResults()
{
  if( !m_riidDisplay )
    return;
  
  SimpleDialog *dialog = m_riidDisplay;
  m_riidDisplay = nullptr;
  dialog->done( WDialog::DialogCode::Accepted );
}//void programmaticallyCloseRiidResults()


void InterSpec::showMultimedia( const SpecUtils::SpectrumType type )
{
  if( m_multimedia )
    programmaticallyCloseMultimediaWindow();
  
  m_multimedia = displayMultimedia( measurment(type) );
  m_multimedia->finished().connect( boost::bind( &InterSpec::handleMultimediaClose, this, m_multimedia ) );
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ programmaticallyCloseMultimediaWindow(); };
    auto redo = [this, type](){ showMultimedia(type); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Multimedia" );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void showMultimedia( const SpecUtils::SpectrumType type )


void InterSpec::handleMultimediaClose( SimpleDialog *dialog )
{
  if( dialog != m_multimedia )
  {
    cerr << "InterSpec::handleMultimediaClose: dialog being closed is not m_multimedia; not doing anything" << endl;
    assert( !m_multimedia );
    return;
  }
  
  m_multimedia = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    // Just assume Foreground so we dont need to copy this info all around
    auto undo = [this](){ showMultimedia(SpecUtils::SpectrumType::Foreground); };
    auto redo = [this](){ programmaticallyCloseMultimediaWindow(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Multimedia" );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
}//void handleMultimediaClose( SimpleDialog *dialog )


void InterSpec::programmaticallyCloseMultimediaWindow()
{
  SimpleDialog *dialog = m_multimedia;
  m_multimedia = nullptr; //Set m_multimedia to nullptr so #handleMultimediaClose wont add undo/redo step
  if( dialog )
    dialog->done( Wt::WDialog::DialogCode::Accepted );
}//void programmaticallyCloseMultimediaWindow()


#if( USE_TERMINAL_WIDGET )
void InterSpec::createTerminalWidget()
{
  if( m_terminal )
    return;
  
  m_terminal = new TerminalWidget( this );
  m_terminal->focusText();
  
  if( m_toolsTabs )
  {
    WMenuItem *item = m_toolsTabs->addTab( m_terminal, WString::tr(TerminalTabTitleKey) );
    item->setCloseable( true );
    m_toolsTabs->setCurrentWidget( m_terminal );
    const int index = m_toolsTabs->currentIndex();
    m_toolsTabs->setTabToolTip( index, WString::tr("app-tab-tt-terminal") );
    
    // Note that the m_toolsTabs->tabClosed() signal has already been hooked up to call
    //  handleToolTabClosed(), which will delete m_terminal when the user closes the tab.
  }else
  {
    m_terminalWindow = new AuxWindow( WString::tr(TerminalTabTitleKey),
                                     (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                      | AuxWindowProperties::EnableResize | AuxWindowProperties::TabletNotFullScreen) );
    
    WPushButton *closeButton = m_terminalWindow->addCloseButtonToFooter();
    closeButton->clicked().connect(m_terminalWindow, &AuxWindow::hide);
    
    AuxWindow::addHelpInFooter( m_terminalWindow->footer(), "terminal-dialog" );
    
    m_terminalWindow->rejectWhenEscapePressed();
    m_terminalWindow->finished().connect( this, &InterSpec::handleTerminalWindowClose );

    m_terminalWindow->show();
    if( m_renderedWidth > 100 && m_renderedHeight > 100 && !isPhone() )
    {
      m_terminalWindow->resizeWindow( 0.95*m_renderedWidth, 0.25*m_renderedHeight );
      m_terminalWindow->centerWindow();
    }
    
    m_terminalWindow->stretcher()->addWidget( m_terminal, 0, 0 );
    
    // We shouldnt block undo/redo any time the terminal window is open
    //new UndoRedoManager::BlockGuiUndoRedo( m_terminalWindow ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
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
    AuxWindow::deleteAuxWindow( m_terminalWindow );
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


#if( USE_REMOTE_RID )
RemoteRid *InterSpec::remoteRid()
{
  return m_remoteRid;
}

void InterSpec::createRemoteRidWindow()
{
  if( m_remoteRid )
    return;
  
  m_remoteRidMenuItem->disable();
  
  auto closeTool = wApp->bind( boost::bind(&InterSpec::deleteRemoteRidWindow, this) );
  auto openTool = wApp->bind( boost::bind(&InterSpec::createRemoteRidWindow, this) );
  
  SimpleDialog *warning = RemoteRid::startRemoteRidDialog( this, [this,closeTool,openTool](AuxWindow *window, RemoteRid *rid ){
    assert( (!window) == (!rid) );
    
    // `window` and `rid` will be nullptr if user canceled the operation
    if( !window || !rid )
    {
      m_remoteRidMenuItem->enable();
      assert( !m_remoteRid && !m_remoteRidWindow );
    }else
    {
      m_remoteRid = rid;
      m_remoteRidWindow = window;
      window->finished().connect( this, &InterSpec::deleteRemoteRidWindow );
      
      if( m_undo && m_undo->canAddUndoRedoNow() )
      {
        auto redo = [this](){
          if( m_remoteRid )
            deleteRemoteRidWindow();
          pair<AuxWindow *, RemoteRid *> res = RemoteRid::createDialog( this );
          if( !res.first )
            return;
          m_remoteRid = res.second;
          m_remoteRidWindow = res.first;
          res.first->finished().connect( this, &InterSpec::deleteRemoteRidWindow );
        };
        
        m_undo->addUndoRedoStep( closeTool, redo, "Show remote RID tool" );
      }//if( undo )
    }
  } );
  
  
  if( warning && m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = wApp->bind( boost::bind(&WDialog::accept, warning) );
    m_undo->addUndoRedoStep( std::move(undo), openTool, "Show remote RID tool" );
  }//if( undo warning )
}//void createRemoteRidWindow()


void InterSpec::deleteRemoteRidWindow()
{
  if( !m_remoteRid )
    return;
  
  m_remoteRidMenuItem->enable();
  
  AuxWindow::deleteAuxWindow( m_remoteRidWindow );
  
  m_remoteRid = nullptr;
  m_remoteRidWindow = nullptr;
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this](){
      if( m_remoteRid )
        deleteRemoteRidWindow();
      pair<AuxWindow *, RemoteRid *> res = RemoteRid::createDialog( this );
      if( !res.first )
        return;
      m_remoteRid = res.second;
      m_remoteRidWindow = res.first;
      res.first->finished().connect( this, &InterSpec::deleteRemoteRidWindow );
    };
    auto redo = wApp->bind( boost::bind(&InterSpec::deleteRemoteRidWindow, this) );
    
    auto openTool = wApp->bind( boost::bind(&InterSpec::createRemoteRidWindow, this) );
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close remote RID tool" );
  }
}//void handleRemoteRidClose()

void InterSpec::setAutoRemoteRidResultDialog( SimpleDialog *dialog )
{
  programaticallyCloseAutoRemoteRidResultDialog();
  m_autoRemoteRidResultDialog = dialog;
}//setAutoRemoteRidResultDialog()


void InterSpec::handleAutoRemoteRidResultDialogClose()
{
  m_autoRemoteRidResultDialog = nullptr;
}//handleAutoRemoteRidResultDialogClose()


void InterSpec::programaticallyCloseAutoRemoteRidResultDialog()
{
  SimpleDialog *dialog = m_autoRemoteRidResultDialog;
  
  //Set m_multimedia to m_autoRemoteRidResultDialog so #handleAutoRemoteRidResultDialogClose
  //  wont add undo/redo step
  m_autoRemoteRidResultDialog = nullptr;
  if( dialog )
    dialog->done( Wt::WDialog::DialogCode::Accepted );
}//programaticallyCloseAutoRemoteRidResultDialog()
#endif  //#if( USE_REMOTE_RID )


#if( USE_REL_ACT_TOOL )
RelActAutoGui *InterSpec::showRelActAutoWindow()
{
  if( !m_relActAutoGui )
  {
    const std::pair<RelActAutoGui *,AuxWindow *> widgets = RelActAutoGui::createWindow( this );
    if( !widgets.first || !widgets.second )
      return m_relActAutoGui;
      
    m_relActAutoGui = widgets.first;
    m_relActAutoWindow  = widgets.second;
    
    m_relActAutoWindow->finished().connect( boost::bind( &InterSpec::handleRelActAutoClose, this ) );
    
    try
    {
      rapidxml::xml_document<char> *relActState = m_dataMeasurement
                                                  ? m_dataMeasurement->relActAutoGuiState()
                                                  : nullptr;
      if( relActState && relActState->first_node() )
      {
        m_relActAutoGui->deSerialize( relActState->first_node() );
        m_relActAutoGui->checkIfInUserConfigOrCreateOne( true );
      }
    }catch( std::exception &e )
    {
      passMessage( "Error setting &quot;Isotopics from nuclides&quot; state to previously used state: "
                  + std::string(e.what()), WarningWidget::WarningMsgHigh );
      
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, ("Error deserializing Rel. Act. GUI state: " + string(e.what())).c_str() );
#endif
      
      assert( 0 );
    }//try / catch
  }else
  {
    const double windowWidth = 0.95 * renderedWidth();
    const double windowHeight = 0.95 * renderedHeight();
    m_relActAutoWindow->resizeWindow( windowWidth, windowHeight );
    
    m_relActAutoWindow->resizeToFitOnScreen();
    m_relActAutoWindow->show();
    m_relActAutoWindow->centerWindow();
  }//if( !m_shieldingSourceFit )
  
  assert( m_relActAutoMenuItem );
  m_relActAutoMenuItem->disable();
  
  
  return m_relActAutoGui;
}//RelActAutoGui *showRelActAutoWindow()


void InterSpec::handleRelActAutoClose()
{
  assert( m_relActAutoMenuItem );
  m_relActAutoMenuItem->enable();
  
  if( m_relActAutoGui )
    saveRelActAutoStateToForegroundSpecMeas();
  
  assert( m_relActAutoWindow );
  if( !m_relActAutoWindow )
    return;
  
  delete m_relActAutoWindow;
  m_relActAutoGui = nullptr;
  m_relActAutoWindow = nullptr;
}//void handleRelActAutoClose()


RelActManualGui *InterSpec::createRelActManualWidget()
{
  assert( m_relActManualMenuItem );
  m_relActManualMenuItem->disable();
  
  if( m_relActManualGui )
    return m_relActManualGui;
 
  const int origTab = m_toolsTabs ? m_toolsTabs->currentIndex() : -1;
  
  m_relActManualGui = new RelActManualGui( this );
  
  if( m_toolsTabs )
  {
    WMenuItem *item = m_toolsTabs->addTab( m_relActManualGui, WString::tr(RelActManualTitleKey) );
    item->setCloseable( true );
    m_toolsTabs->setCurrentWidget( m_relActManualGui );
    const int index = m_toolsTabs->currentIndex();
    m_toolsTabs->setTabToolTip( index, WString::tr("app-tab-tt-rel-eff") );
    
    // Note that the m_toolsTabs->tabClosed() signal has already been hooked up to call
    //  handleToolTabClosed(), which will delete m_relActManualGui when the user closes the tab.
  }else
  {
    m_relActManualWindow = new AuxWindow( WString::tr("window-title-peak-rel-eff"),
                                     (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                      | AuxWindowProperties::EnableResize | AuxWindowProperties::TabletNotFullScreen) );
    
    m_relActManualWindow->rejectWhenEscapePressed();
    m_relActManualWindow->finished().connect( this, &InterSpec::handleRelActManualClose );
    
    m_relActManualWindow->show();
    if( (m_renderedWidth > 100) && (m_renderedHeight > 100) && !isPhone() )
    {
      m_relActManualWindow->resizeWindow( 0.95*m_renderedWidth, 0.80*m_renderedHeight );
      m_relActManualWindow->centerWindow();
    }
    
    m_relActManualWindow->stretcher()->addWidget( m_relActManualGui, 0, 0 );
    
    
    WPushButton *closeButton = m_relActManualWindow->addCloseButtonToFooter();
    closeButton->clicked().connect(m_relActManualWindow, &AuxWindow::hide);
  }//if( toolTabsVisible() )
  
  assert( m_relActManualMenuItem );
  m_relActManualMenuItem->disable();
  
  try
  {
    rapidxml::xml_document<char> *relActSate = m_dataMeasurement
                                                 ? m_dataMeasurement->relActManualGuiState()
                                                 : nullptr;
    if( relActSate && relActSate->first_node() )
      m_relActManualGui->deSerialize( relActSate->first_node() );
  }catch( std::exception &e )
  {
    passMessage( "Error setting &quot;Isotopics from peaks&quot; state to previously used state: "
                 + std::string(e.what()), WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, ("Error deserializing Rel. Act. GUI state: " + string(e.what())).c_str() );
#endif
    
    assert( 0 );
  }//try / catch
  
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [this, origTab](){
      if( !m_relActManualWindow && m_relActManualGui )
        m_toolsTabs->removeTab( m_relActManualGui );
      
      handleRelActManualClose();
      
      if( (origTab >= 0) && m_toolsTabs )
        m_toolsTabs->setCurrentIndex( origTab );
    };
    
    auto redo = [this](){ createRelActManualWidget(); };
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show 'Isotopics from peaks' tool" );
  }//if( m_undo && !m_undo->canAddUndoRedoNow() )
  
  return m_relActManualGui;
}//RelActManualGui *createRelActManualWidget()


void InterSpec::handleRelActManualClose()
{
  assert( m_relActManualGui );
  
  if( m_relActManualGui )
    saveRelActManualStateToForegroundSpecMeas();
  
  if( m_relActManualWindow )
  {
    AuxWindow::deleteAuxWindow( m_relActManualWindow );
  }else
  {
    if( m_relActManualGui )
      delete m_relActManualGui;
    
    if( m_toolsTabs )
      m_toolsTabs->setCurrentIndex( 2 );
  }//if( m_relActManualWindow ) / else
  
  m_relActManualGui = nullptr;
  m_relActManualWindow = nullptr;

  assert( m_relActManualMenuItem );
  m_relActManualMenuItem->enable();
}//void InterSpec::handleRelActManualClose()


void InterSpec::saveRelActManualStateToForegroundSpecMeas()
{
  if( !m_relActManualGui || !m_dataMeasurement )
    return;
  
  string xml_data;
  std::unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
  
  m_relActManualGui->serialize( doc.get() );
  
  m_dataMeasurement->setRelActManualGuiState( std::move(doc) );
}//void saveRelActManualStateToForegroundSpecMeas()


void InterSpec::saveRelActAutoStateToForegroundSpecMeas()
{
  if( !m_relActAutoGui || !m_dataMeasurement )
    return;
  
  string xml_data;
  std::unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
  
  m_relActAutoGui->serialize( doc.get() );
  
  m_dataMeasurement->setRelActAutoGuiState( std::move(doc) );
}//void saveRelActAutoStateToForegroundSpecMeas()

#endif //#if( USE_REL_ACT_TOOL )


#if( USE_TERMINAL_WIDGET || USE_REL_ACT_TOOL )
void InterSpec::handleToolTabClosed( const int tabnum )
{
  assert( m_toolsTabs );
  if( !m_toolsTabs )
    return;
  
  WWidget *w = m_toolsTabs->widget( tabnum );
  
#if( USE_TERMINAL_WIDGET && USE_REL_ACT_TOOL )
  if( w == m_relActManualGui )
  {
    handleRelActManualClose();
  }else if( w == m_terminal )
  {
    handleTerminalWindowClose();
  }else
  {
    assert( 0 );
  }
#elif( USE_TERMINAL_WIDGET )
  handleTerminalWindowClose();
#elif( USE_REL_ACT_TOOL )
  handleRelActManualClose();
  //static_assert( 0, "Need to update handleToolTabClosed logic" );  //20230913 - no updates look to be needed
#endif
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  m_currentToolsTab = m_toolsTabs->currentIndex();
#endif
}//void handleToolTabClosed( const int tabnum )
#endif


void InterSpec::addToolsMenu( Wt::WWidget *parent )
{
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  
  PopupDivMenu *parentMenu = dynamic_cast<PopupDivMenu *>( parent );
  WContainerWidget *menuDiv = dynamic_cast<WContainerWidget *>( parent );
  if( !parentMenu && !menuDiv )
    throw runtime_error( "InterSpec::addToolsMenu(): parent passed in"
                        " must be a PopupDivMenu  or WContainerWidget" );
  
  PopupDivMenu *popup = NULL;
  
  if( menuDiv )
  {
    WPushButton *button = new WPushButton( WString::tr("app-menu-tools"), menuDiv );
    button->addStyleClass( "MenuLabel" );
    popup = new PopupDivMenu( button, PopupDivMenu::AppLevelMenu );
  }else
  {
    popup = parentMenu->addPopupMenuItem( WString::tr("app-menu-tools") );
  }
  
  m_toolsMenuPopup = popup;
  
  PopupDivMenuItem *item = NULL;

  item = popup->addMenuItem( WString::tr("app-mi-tools-act-fit") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-act-fit") , showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::shieldingSourceFit, this ) );
 
#if( USE_REL_ACT_TOOL )
  // The Relative Efficiency tools are not specialized to display on phones yet.
  if( !isPhone() )
  {
    m_relActAutoMenuItem = popup->addMenuItem( WString::tr("app-mi-tools-iso-by-nuc") );
    HelpSystem::attachToolTipOn( m_relActAutoMenuItem, WString::tr("app-mi-tt-tools-iso-by-nuc") , showToolTips );
    m_relActAutoMenuItem->triggered().connect( boost::bind( &InterSpec::showRelActAutoWindow, this ) );
    
    m_relActManualMenuItem = popup->addMenuItem( WString::tr("app-mi-tools-iso-by-peak") );
    HelpSystem::attachToolTipOn( m_relActManualMenuItem, WString::tr("app-mi-tt-tools-iso-by-peak") , showToolTips );
    m_relActManualMenuItem->triggered().connect( boost::bind( &InterSpec::createRelActManualWidget, this ) );
  }//if( !isPhone() )
#endif
  
  item = popup->addMenuItem( WString::tr("app-mi-tools-xs-calc") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-xs-calc"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showGammaXsTool, this ) );
    
  item = popup->addMenuItem( WString::tr("app-mi-tools-dose-calc") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-dose-calc"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showDoseTool, this ) );
  
//  item = popup->addMenuItem( WString::fromUTF8("1/r² Calculator") );  // is superscript 2
#if( USE_OSX_NATIVE_MENU )
  item = popup->addMenuItem( WString::fromUTF8("1/r\x0032 Calculator") );  //works on OS X at least.
#else
  item = popup->addMenuItem( WString::tr("app-mi-tools-1r2") );
  item->makeTextXHTML();
#endif
  
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-1r2"), showToolTips );
  
  item->triggered().connect( boost::bind( &InterSpec::createOneOverR2Calculator, this ) );

  item = popup->addMenuItem( WString::tr("app-mi-tools-unit") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-unit"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::createUnitsConverterTool, this ) );
  
  item = popup->addMenuItem( WString::tr("app-mi-tools-flux") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-flux"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::createFluxTool, this ) );
  
  item = popup->addMenuItem( WString::tr("app-mi-tools-nuc-decay") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-nuc-decay"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::createDecayInfoWindow, this ) );

#if( USE_DETECTION_LIMIT_TOOL )
  if( !isPhone() )
  {
    // TODO: modify formatting of Detection Confidence tool to work on phones
    item = popup->addMenuItem( WString::tr("app-mi-tools-mda") );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-mda"), showToolTips );
    item->triggered().connect( boost::bind( &InterSpec::showDetectionLimitTool, this, string() ) );
  }//if( !isPhone() )
#endif
  
  item = popup->addMenuItem( WString::tr("app-mi-tools-select-drf") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-select-drf"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showDrfSelectWindow, this ) );
  
  if( !isPhone() )
  {
    // The create DRF tool wont layout very good on phones - I dont think this is "phone"
    //  work anyway, so we'll just skip it.
    item = popup->addMenuItem( WString::tr("app-mi-tools-make-drf") );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-make-drf"), showToolTips );
    item->triggered().connect( boost::bind( &InterSpec::showMakeDrfWindow, this ) );
  }
  
  item = popup->addMenuItem( WString::tr("app-mi-tools-file-par") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-file-par"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::createFileParameterWindow, this, SpecUtils::SpectrumType::Foreground ) );

  item = popup->addMenuItem( WString::tr("app-mi-tools-en-sum") );
  HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-en-sum"), showToolTips );
  item->triggered().connect( boost::bind( &InterSpec::showGammaCountDialog, this ) );
  
#if( USE_SPECRUM_FILE_QUERY_WIDGET )
  
#if( BUILD_AS_OSX_APP )
  const bool addQueryTool = InterSpecApp::isPrimaryWindowInstance();
#else
  const bool addQueryTool = true;
#endif
  
  if( addQueryTool )
  {
    item = popup->addMenuItem( WString::tr("app-mi-tools-query") );
    HelpSystem::attachToolTipOn( item, WString::tr("app-mi-tt-tools-query"), showToolTips );
    item->triggered().connect( this, &InterSpec::showFileQueryDialog );
  }//if( addQueryTool )
#endif
  
#if( USE_TERMINAL_WIDGET )
  m_terminalMenuItem = popup->addMenuItem( WString::tr("app-mi-tools-math") );
  HelpSystem::attachToolTipOn( m_terminalMenuItem, WString::tr("app-mi-tt-tools-math"), showToolTips );
  m_terminalMenuItem->triggered().connect( this, &InterSpec::createTerminalWidget );
#endif
  
#if( USE_REMOTE_RID )
  m_remoteRidMenuItem = popup->addMenuItem( WString::tr("app-mi-tools-ext-rid") );
  
  WString extRidTT = WString::tr("app-mi-tt-tools-ext-rid");

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
  extRidTT.arg( WString::tr("app-mi-tt-tools-ext-rid-exe") );
#else
  extRidTT.arg( "" );
#endif
  
  HelpSystem::attachToolTipOn( m_terminalMenuItem, extRidTT, showToolTips );
  m_remoteRidMenuItem->triggered().connect( this, &InterSpec::createRemoteRidWindow );
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
    const string materialfile = SpecUtils::append_path(ns_staticDataDirectory, "MaterialDataBase.txt" );
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
    popupOptions.listSeparator      = '\0';            //(char) When editing a list of values, the separator used for different items.
    popupOptions.whitespace         = "";       //When editing a value, the whitespace characters ignored before the current value.
    //popupOptions.wordSeparators     = "-_., ;()";     //To show suggestions based on matches of the edited value with parts of the suggestion.
    popupOptions.wordStartRegexp = "\\s|^|\\(|\\<";       // Instead of using .wordSeparators, we will use the regex option to start matching at whitespaces, start of line, open-paren, and boundaries of words (probably a bit duplicative).
    popupOptions.appendReplacedText = "";             //
    
    // We may want to parent `m_shieldingSuggestion...
    
    m_shieldingSuggestion = new WSuggestionPopup( popupOptions, this );
    
    m_shieldingSuggestion->addStyleClass("suggestion");
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
    m_shieldingSuggestion->setJavaScriptMember("wtNoReparent", "true");
#endif
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


GammaXsWindow *InterSpec::showGammaXsTool()
{
  if( !m_gammaXsToolWindow )
  {
    m_gammaXsToolWindow = new GammaXsWindow( m_materialDB.get(), m_shieldingSuggestion, this );
    m_gammaXsToolWindow->finished().connect( this, &InterSpec::deleteGammaXsTool );
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){ deleteGammaXsTool(); };
      auto redo = [this](){ showGammaXsTool(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Gamma XS Tool." );
    }//if( m_undo && m_undo->canAddUndoRedoNow() )
  }else
  {
    m_gammaXsToolWindow->show();
    m_gammaXsToolWindow->centerWindowHeavyHanded();
  }
  
  return m_gammaXsToolWindow;
}//showGammaXsTool()


void InterSpec::deleteGammaXsTool()
{
  if( m_gammaXsToolWindow )
  {
    GammaXsGui *tool = m_gammaXsToolWindow->xstool();
    
    if( m_undo && tool )
    {
      const string uri_query = tool->encodeStateToUrl();
      
      auto undo = [this,uri_query](){
        try
        {
          GammaXsWindow *xswindow = showGammaXsTool();
          GammaXsGui *tool = xswindow ? xswindow->xstool() : nullptr;
          if( tool )
            tool->handleAppUrl( uri_query );
        }catch( std::exception &e )
        {
          Wt::log("error") << "Error restoring gamma XS tool state: '"
                           << e.what() << "', from URI='" << uri_query << "'";
          passMessage( "Error restoring Gamma XS Tool State", WarningWidget::WarningMsgInfo );
        }//try /catch
      };//undo
      
      auto redo = [this](){ deleteGammaXsTool(); };
      
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close XS Tool" );
    }//if( m_undo && tool )
    
    AuxWindow::deleteAuxWindow( m_gammaXsToolWindow );
  }//if( m_gammaXsToolWindow )
  
  m_gammaXsToolWindow = nullptr;
}//void deleteGammaXsTool()


DoseCalcWindow *InterSpec::showDoseTool()
{
  if( !m_doseCalcWindow )
  {
    m_doseCalcWindow = new DoseCalcWindow( m_materialDB.get(), m_shieldingSuggestion, this );
    m_doseCalcWindow->finished().connect( this, &InterSpec::deleteDoseCalcTool );
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this](){ deleteDoseCalcTool(); };
      auto redo = [this](){ showDoseTool(); };
      m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Dose Tool." );
    }//if( m_undo && m_undo->canAddUndoRedoNow() )
  }else
  {
    m_doseCalcWindow->show();
    m_doseCalcWindow->centerWindowHeavyHanded();
  }
  
  return m_doseCalcWindow;
}//DoseCalcWindow *showDoseTool()


void InterSpec::deleteDoseCalcTool()
{
  if( m_doseCalcWindow )
  {
    DoseCalcWidget *tool = m_doseCalcWindow->tool();
    string uri = tool ? tool->encodeStateToUrl() : string();
    auto undo = [this,uri](){
      // TODO: the GUI layout can be a bit messed up when undoing the close - not totally sure why.
      DoseCalcWindow *dosewin = showDoseTool();
      DoseCalcWidget *tool = dosewin ? dosewin->tool() : nullptr;
      if( !tool )
        return;
      try
      {
        const size_t pos = uri.find('?');
        if( pos == string::npos )
          throw runtime_error( "Couldnt find path and query components of uri." );
        string path = uri.substr(0, pos);
        string query = uri.substr( pos + 1 );
        tool->handleAppUrl( path, query );
      }catch( std::exception &e )
      {
        passMessage( "Error re-opening Dose Calc tool: " + string(e.what()),
                    WarningWidget::WarningMsgHigh );
      }
      
    };//undo
    auto redo = [this](){ deleteDoseCalcTool(); };
    
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Dose Tool" );
    
    AuxWindow::deleteAuxWindow( m_doseCalcWindow );
  }
  m_doseCalcWindow = nullptr;
}//void deleteDoseCalcTool();


MakeDrfWindow *InterSpec::showMakeDrfWindow()
{
  assert( !m_makeDrfTool );
  if( !m_makeDrfTool )
  {
    // Disable the "Make Detector Response" Tools menu item
    assert( m_toolsMenuPopup );
    for( WMenuItem *item : m_toolsMenuPopup->items() )
    {
      if( item->text().key() == "app-mi-tools-make-drf" )
        item->setDisabled( true );
    }//for( loop over "Tools" menu items )
    
    m_makeDrfTool = new MakeDrfWindow( this, m_materialDB.get(), m_shieldingSuggestion );
    m_makeDrfTool->finished().connect( boost::bind( &InterSpec::handleCloseMakeDrfWindow, this, m_makeDrfTool ) );
    new UndoRedoManager::BlockGuiUndoRedo( m_makeDrfTool ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
  }//if( !m_makeDrfTool )
  
  return m_makeDrfTool;
}//void showDrfSelectWindow()


MakeDrfWindow *InterSpec::makeDrfWindow()
{
  return m_makeDrfTool;
}


void InterSpec::handleCloseMakeDrfWindow( MakeDrfWindow *window )
{
  assert( window == m_makeDrfTool );
  if( !window || (window != m_makeDrfTool) )
    return;
  
  // Enable the "Make Detector Response" Tools menu item
  assert( m_toolsMenuPopup );
  for( WMenuItem *item : m_toolsMenuPopup->items() )
  {
    if( item->text().key() == "app-mi-tools-make-drf" )
      item->setDisabled( false );
  }//for( loop over "Tools" menu items )
  
  m_makeDrfTool = nullptr;
  AuxWindow::deleteAuxWindow( window );
}//void handleCloseMakeDrfWindow( MakeDrfWindow *window )


DrfSelectWindow *InterSpec::showDrfSelectWindow()
{
  auto redu = [this](){
    if( !m_drfSelectWindow )
      m_drfSelectWindow = new DrfSelectWindow( this );
    else
      m_drfSelectWindow->show();
  };
  
  auto undo = [this](){
    if( m_drfSelectWindow )
    {
      DrfSelectWindow *window = m_drfSelectWindow;
      m_drfSelectWindow = nullptr;
      AuxWindow::deleteAuxWindow( window );
    }
  };
  
  redu();
  
  if( m_undo )
    m_undo->addUndoRedoStep( undo, redu, "Show DRF Select Window" );
  
  return m_drfSelectWindow;
}//void showDrfSelectWindow()


void InterSpec::closeDrfSelectWindow()
{
  auto undo = [this](){
    if( !m_drfSelectWindow )
      m_drfSelectWindow = new DrfSelectWindow( this );
    else
      m_drfSelectWindow->show();
  };
  
  auto redo = [this](){
    if( m_drfSelectWindow )
    {
      DrfSelectWindow *window = m_drfSelectWindow;
      m_drfSelectWindow = nullptr;
      AuxWindow::deleteAuxWindow( window );
    }
  };
  
  redo();
  
  if( m_undo )
    m_undo->addUndoRedoStep( undo, redo, "Close DRF Select Window" );
}//void closeDrfSelectWindow()


void InterSpec::showCompactFileManagerWindow()
{
 auto *compact = new CompactFileManager( m_fileManager, this, CompactFileManager::Tabbed );

  m_spectrum->yAxisScaled().connect( boost::bind( &CompactFileManager::handleSpectrumScale, compact,
                                                 boost::placeholders::_1,
                                                 boost::placeholders::_2,
                                                 boost::placeholders::_3 ) );
  
  AuxWindow *window = new AuxWindow( WString::tr("window-title-compact-file"),
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::SetCloseable)
                                     |AuxWindowProperties::TabletNotFullScreen) );
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

    m_toolsTabs->addTab( m_nuclideSearchContainer, WString::tr(NuclideSearchTabTitleKey), TabLoadPolicy );
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
  
  
  m_nuclideSearchWindow = new AuxWindow( WString::tr(NuclideSearchTabTitleKey),
                                        (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                         | AuxWindowProperties::EnableResize
                                         | AuxWindowProperties::SetCloseable) );
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
  
 
  Wt::WPushButton *closeButton = m_nuclideSearchWindow->addCloseButtonToFooter( WString::tr("Close"),true);
  
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


ShieldingSourceDisplay *InterSpec::shieldingSourceFit()
{
  if( m_shieldingSourceFit )
    return m_shieldingSourceFit;
  
  assert( !m_shieldingSourceFitWindow );
  
  assert( m_peakInfoDisplay );
  auto widgets = ShieldingSourceDisplay::createWindow( this );
  
  m_shieldingSourceFit = widgets.first;
  m_shieldingSourceFitWindow = widgets.second;
  
  if( m_undo && !m_undo->isInUndoOrRedo() )
  {
    auto undo = [this](){ closeShieldingSourceFit(); };
    auto redo = [this](){ shieldingSourceFit(); };
    m_undo->addUndoRedoStep( undo, redo, "Open Activity/Shielding Fit tool" );
  }
  
  return m_shieldingSourceFit;
}//ShieldingSourceDisplay *shieldingSourceFit()


void InterSpec::saveShieldingSourceModelToForegroundSpecMeas()
{
  if( !m_shieldingSourceFitWindow || !m_shieldingSourceFit || !m_dataMeasurement )
    return;
  
  string xml_data;
  std::unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
  
  m_shieldingSourceFit->serialize( doc.get() );
  
  m_dataMeasurement->setShieldingSourceModel( std::move(doc) );
}//void saveShieldingSourceModelToForegroundSpecMeas()



void InterSpec::closeShieldingSourceFit()
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
  
  if( m_undo && !m_undo->isInUndoOrRedo() )
  {
    auto undo = [this](){ shieldingSourceFit(); };
    auto redo = [this](){ closeShieldingSourceFit(); };
    m_undo->addUndoRedoStep( undo, redo, "Close Activity/Shielding Fit tool" );
  }
}//void closeShieldingSourceFit()



void InterSpec::showHelpWindow( const std::string &preselect )
{
  if( m_helpWindow )
  {
    m_helpWindow->show();
    m_helpWindow->centerWindow();
    m_helpWindow->setTopic( preselect );
    return;
  }//if( m_helpWindow )
  
  m_helpWindow = new HelpSystem::HelpWindow( preselect );
  m_helpWindow->finished().connect( this, &InterSpec::closeHelpWindow );
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [this](){ closeHelpWindow(); };
    auto redo = [this,preselect](){ showHelpWindow(preselect);};
    
    undoRedo->addUndoRedoStep( undo, redo, "Show help window" );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void InterSpec::showHelpWindow( const std::string &preselect )


void InterSpec::closeHelpWindow()
{
  if( !m_helpWindow )
    return;
  
  const string topic = m_helpWindow->currentTopic();
  
  AuxWindow::deleteAuxWindow( m_helpWindow );
  m_helpWindow = nullptr;
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [this, topic](){ showHelpWindow( topic ); };
    auto redo = [this](){ closeHelpWindow(); };
    
    undoRedo->addUndoRedoStep( undo, redo, "Close help window" );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void closeHelpWindow()


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
    const ReferenceLineInfo &lines = m_referencePhotopeakLines->currentlyShowingNuclide();
    
    if( (lines.m_validity == ReferenceLineInfo::InputValidity::Valid)
        || m_referencePhotopeakLines->persistedNuclides().size() )
    {
      m_referencePhotopeakLines->serialize( xml_state );
    }
    
    m_referencePhotopeakLines->clearAllLines();
    delete m_referencePhotopeakLines;
    m_referencePhotopeakLines = NULL;
  }//if( m_referencePhotopeakLines )

  m_referencePhotopeakLinesWindow = new AuxWindow( WString::tr(GammaLinesTabTitleKey),
                                                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                                                   | AuxWindowProperties::EnableResize
                                                   | AuxWindowProperties::SetCloseable)
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

  Wt::WPushButton *closeButton = m_referencePhotopeakLinesWindow->addCloseButtonToFooter( WString::tr("Close"),true);
  
  if( isPhone() )
  {
    m_referencePhotopeakLines->displayingNuclide().connect( boost::bind( &WPushButton::setText, closeButton, WString::tr("show-lines")) );
    m_referencePhotopeakLines->nuclidesCleared().connect( boost::bind( &WPushButton::setText, closeButton, WString::tr("Close")) );
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
    const ReferenceLineInfo &lines = m_referencePhotopeakLines->currentlyShowingNuclide();
    if( (lines.m_validity == ReferenceLineInfo::InputValidity::Valid)
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
    
    
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
    if( isPhone() && (renderedHeight() > renderedWidth()) )
      m_referencePhotopeakLines->setNarrowPhoneLayout( true );
#endif
    
    m_toolsTabs->addTab( m_referencePhotopeakLines, WString::tr(GammaLinesTabTitleKey), TabLoadPolicy );
    
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
  
  const int prev_tab = m_currentToolsTab;
  
  auto handle_change = [this]( const int current_tab, const bool focus  ){
    const int refTab = m_toolsTabs->indexOf(m_referencePhotopeakLines);
    const int calibtab = m_toolsTabs->indexOf(m_energyCalTool);
    const int searchTab = m_toolsTabs->indexOf(m_nuclideSearchContainer);
    
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    
    if( m_referencePhotopeakLines && focus && (current_tab == refTab) && app && !app->isMobile() )
      m_referencePhotopeakLines->setFocusToIsotopeEdit();
    
    if( m_nuclideSearch && (m_currentToolsTab==searchTab) )
      m_nuclideSearch->clearSearchEnergiesOnClient();
    
    if( m_nuclideSearch && (current_tab==searchTab) )
      m_nuclideSearch->loadSearchEnergiesToClient();
    
    if( focus && (current_tab == calibtab) )
    {
      if( !m_infoNotificationsMade.count("recal-tab")
         && UserPreferences::preferenceValue<bool>( "ShowTooltips", this )
         && !isMobile() )
      {
        m_infoNotificationsMade.insert( "recal-tab" );
        passMessage( WString::tr("info-recal-tab-selected"), WarningWidget::WarningMsgInfo );
      }
    }//if( tab == calibtab )
    
    m_currentToolsTab = current_tab;
  };//handle_change
  
  handle_change( tab, true );
  
  if( !m_undo )
    return;
  
  auto undo = [this, prev_tab, handle_change](){
    if( m_toolsTabs && (prev_tab < m_toolsTabs->count()) )
    {
      //WTabWidget::currentChanged() is emited even when you programmatically change the tab; we
      //  need to block this here, or else the undo/redo history will get all messed up.
      const bool orig_state = m_toolsTabs->currentChanged().isBlocked();
      m_toolsTabs->currentChanged().setBlocked(true);
      m_toolsTabs->setCurrentIndex(prev_tab);
      m_toolsTabs->currentChanged().setBlocked(orig_state);
    }
    handle_change(prev_tab, false);
  };//undo
  
  auto redo = [this, tab, handle_change](){
    if( m_toolsTabs && (tab < m_toolsTabs->count()) )
    {
      // See note above on blocking the WTabWidget::currentChanged() signal
      const bool orig_state = m_toolsTabs->currentChanged().isBlocked();
      m_toolsTabs->currentChanged().setBlocked(true);
      m_toolsTabs->setCurrentIndex(tab);
      m_toolsTabs->currentChanged().setBlocked(orig_state);
    }
    
    handle_change(tab, false);
  };//redo
  
  m_undo->addUndoRedoStep( undo, redo, "Change tool tab." );
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
}//double sample_real_time_increment()


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
                 + ", " + std::to_string(sample_end) + ") - this shouldn't have happened"
                 " - not changing sample numbers of spectrum",
                 WarningWidget::WarningMsgHigh );
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


#if( APPLY_OS_COLOR_THEME_FROM_JS && !BUILD_AS_OSX_APP && !BUILD_AS_ELECTRON_APP )
void InterSpec::initOsColorThemeChangeDetect()
{
  m_osColorThemeChange.reset( new JSignal<std::string>( this, "OsColorThemeChange", false ) );
  m_osColorThemeChange->connect( boost::bind( &InterSpec::osThemeChange, this,
                                             boost::placeholders::_1 ) );
  
  LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsDetectOsColorThemeJs);
  LOAD_JAVASCRIPT(wApp, "js/InterSpec.js", "InterSpec", wtjsSetupOsColorThemeChangeJs);
  
  doJavaScript( "Wt.WT.DetectOsColorThemeJs('" + id() + "')" );
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
      det = DrfSelect::initARelEffDetector( type, manufacturer, model );
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
    //      can be found by searching for "WarningMsgShowOnBoardRiid"
    
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
    
    const char *msg_key = (usingUserDefaultDet ? "info-user-default-drf" : "info-app-default-drf");
    passMessage( WString::tr(msg_key), WarningWidget::WarningMsgInfo );
    
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
  UndoRedoManager::BlockUndoRedoInserts blocker;
  
  const int spectypeindex = static_cast<int>( spec_type );
  
  vector< boost::function<void(void)> > furtherworkers;
  
  const bool wasModified = (meas ? meas->modified() : false);
  const bool wasModifiedSinceDecode = (meas ? meas->modified_since_decode() : false);
  
  if( m_useInfoWindow && meas )
  {
    // If we are loading a state from the "Welcome To InterSpec" screen, we dont want to delete
    //  m_useInfoWindow because we will still use it, so instead we'll try deleting the window on
    //  the next go around of the event loop.
    auto deleter = wApp->bind( boost::bind( &InterSpec::deleteWelcomeDialog, this, false) );
    auto doDelete = [deleter](){ deleter(); wApp->triggerUpdate(); };
    WServer::instance()->post( wApp->sessionId(), doDelete );
  }//if( meas )
  
  std::shared_ptr<SpecMeas> previous = measurment(spec_type);
  const set<int> prevsamples = displayedSamples(spec_type);
  const vector<string> prevdets = detectorsToDisplay(spec_type);
  const bool sameSpecFile = (meas==previous);
  std::shared_ptr<const SpecUtils::Measurement> prev_display = m_spectrum->histUsedForXAxis();
  

  if( (spec_type == SpecUtils::SpectrumType::Foreground) && previous && (previous != meas) )
  {
#if( USE_DB_TO_STORE_SPECTRA )
    saveStateAtForegroundChange( true ); //This function will check "AutoSaveSpectraToDb" pref
#endif //#if( USE_DB_TO_STORE_SPECTRA )
    
    // Close Shielding/Source fit Window
    if( m_shieldingSourceFitWindow )
    {
      delete m_shieldingSourceFitWindow;
      m_shieldingSourceFitWindow = nullptr;
      m_shieldingSourceFit = nullptr;
    }
    
    
    if( m_relActAutoGui )
    {
      // TODO: Should we close this?
    }
    
    if( m_riidDisplay )
      programmaticallyCloseRiidResults();
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
      if( !sameSpecFile )
        deletePeakEdit();
      
      m_exportSpecFileMenu->setDisabled( !meas );
    break;
    
    case SpecUtils::SpectrumType::SecondForeground:
    case SpecUtils::SpectrumType::Background:
    break;
  }//switch( spec_type )

  if( meas
      /*&& sample_numbers.empty() */
      && !sameSpecFile
      && (spec_type==meas->displayType())
      && meas->displayedSampleNumbers().size() )
  {
    sample_numbers = meas->displayedSampleNumbers();
  }
  
  
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
          const bool doLoadDefault = UserPreferences::preferenceValue<bool>( "LoadDefaultDrf", this );
          
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

      if( !sameSpecFile && m_shieldingSourceFit )
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
    passMessage( msg, 0 );

  
  //If loading a new foreground that has a different number of channels than
  //  the background/secondary, and differnet number of bins than previous
  //  foreground, get rid of the background/secondary since the user probably
  //  isnt interested in files from a completely different detector anymore.
  if( spec_type == SpecUtils::SpectrumType::Foreground && m_dataMeasurement && !sameSpecFile )
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
      m_secondDataMeasurement = nullptr;
      m_spectrum->setSecondData( nullptr );
      
      m_displayedSpectrumChangedSignal.emit( SpecUtils::SpectrumType::SecondForeground,
                                             nullptr, {}, {} );
    }//if( num_sec_channel )
    
    if( diff_fore_nchan && num_back_channel && num_foreground_channels && (num_back_channel != num_foreground_channels) )
    {
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
      updateSaveWorkspaceMenu();
#endif
      displayForegroundData( false );
      displayTimeSeriesData();
      
    case SpecUtils::SpectrumType::SecondForeground:
      displaySecondForegroundData();
      
    case SpecUtils::SpectrumType::Background:
      displayBackgroundData();
  };//switch( spec_type )
  
  
  if( !sameSpecFile )
  {
    const bool enableScaler = (m_secondDataMeasurement || m_backgroundMeasurement);
    const int windex_0 = m_spectrum->yAxisScalersIsVisible() ? 0 : 1;
    const int windex_1 = m_spectrum->yAxisScalersIsVisible() ? 1 : 0;
    m_showYAxisScalerItems[windex_0]->hide();
    m_showYAxisScalerItems[windex_1]->show();
  }//if( !sameSpecFile )
  
  //Making fcn call take current data as a argument so that if this way a
  //  recalibration happens (which will change m_spectrum->data()), then the
  //  peak fit routine will get the correct data to use.
  std::function<void(std::shared_ptr<const SpecUtils::Measurement>)> propigate_peaks_fcns;
  
  const bool askToPropigatePeaks
          = UserPreferences::preferenceValue<bool>( "AskPropagatePeaks", this );
  if( askToPropigatePeaks
     && options.testFlag(InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal)
     && !sameSpecFile && meas && m_dataMeasurement && previous && m_spectrum->data()
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
      && !sameSpecFile && m_energyCalTool && !!meas && !!m_dataMeasurement )
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
  }//if( !sameSpecFile && m_energyCalTool && !!meas )
  
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
        passMessage( WString::tr("err-no-foreground"), WarningWidget::WarningMsgHigh );
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
  
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  if( spec_type == SpecUtils::SpectrumType::Foreground )
  {
    const bool hasGps = (meas && meas->has_gps_info());
    if( m_mapMenuItem )
      m_mapMenuItem->setDisabled( !hasGps );
    
    if( hasGps )
    {
      const std::string type = SpecUtils::descriptionText(spec_type);
      
      WString js = WString("{1}"
                 "<div onclick="
                 "\"Wt.emit( $('.specviewer').attr('id'),{name:'miscSignal'}, 'showMap-{2}');"
                 "try{$(this.parentElement.parentElement).remove();}catch(e){}"
                 "return false;\" "
                 "class=\"clearsession\">"
                 "<span class=\"clearsessiontxt\">{3}</span></div>"
                 )
        .arg( WString::tr("info-has-gps-txt") )
        .arg( type )
        .arg( WString::tr("info-has-gps-btn") );
      
      m_warnings->addMessageUnsafe( js, WarningWidget::WarningMsgShowOnBoardRiid, 20000 );
    }//if( hasGps )
  }//if( spec_type == SpecUtils::SpectrumType::Foreground )
#endif //#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  
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
  
  checkEnableViewImageMenuItem();
  
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
  if( meas && !sameSpecFile && !(options & InterSpec::SetSpectrumOptions::SkipParseWarnings) )
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
            passMessage( msg, WarningWidget::WarningMsgMedium );
          app->triggerUpdate();
        }//
      }//if( !givenwarnings.empty() )
    };//checkForWarnings lamda
    
    furtherworkers.push_back( checkForWarnings );
  }//if( meas && !sameSpecFile )
  
  // Check if there are RIID analysis results in the file, and if so let the user know.
  if( !sameSpecFile && meas && meas->detectors_analysis()
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
      
      WString js = WString(
        "{1}: {2}"
        "<div onclick=\"Wt.emit( $('.specviewer').attr('id'),{name:'miscSignal'}, 'showRiidAna-{3}');try{$(this.parentElement.parentElement).remove();}catch(e){} return false;\" class=\"clearsession\">"
          "<span class=\"clearsessiontxt\">{4}</span>"
        "</div>")
      .arg( WString::tr("info-has-rid-txt") )
      .arg( riidAnaSummary(meas) )
      .arg( type )
      .arg( WString::tr("info-has-rid-btn") )
      ;
      
      m_warnings->addMessageUnsafe( js, WarningWidget::WarningMsgShowOnBoardRiid, 20000 );
    }//if( nusedfor == 1 )
  }//if( meas && !sameSpecFile )
  
#if( USE_REMOTE_RID )
  if( meas && !options.testFlag(SetSpectrumOptions::SkipExternalRid) )
  {
    const ExternalRidAuotCallPref pref = RemoteRid::external_rid_call_pref( this );
    
    if( pref != ExternalRidAuotCallPref::DoNotCall )
    {
      Wt::WApplication *app = wApp;
      auto callExternalRid = [app,this,sameSpecFile,pref](){
        WApplication::UpdateLock lock(app);
        if( lock )
        {
          Wt::WFlags<RemoteRid::AnaFileOptions> flags;
          
          // When a search or portal file is loaded as foreground we could be a little smarter
          //  in deciding if we should submit the whole file, or just the displayed sample
          if( sameSpecFile )
            flags |= RemoteRid::AnaFileOptions::OnlyDisplayedSearchSamples;
          
          RemoteRid::startAutomatedOnLoadAnalysis( this, flags );
        }else
        {
          cerr << "Failed to get WApplication::UpdateLock to call external RID." << endl;
        }
      };//callExternalRid lamda
      
      const string appid = wApp->sessionId();
      
      auto postExternalRid = [appid, callExternalRid](){
        WServer::instance()->post( appid, callExternalRid );
      };
      
      furtherworkers.push_back( postExternalRid );
    }//if( user selected to call either external REST API, or external EXE )
  }//if( meas )
#endif //#if( USE_REMOTE_RID )
  
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
      passMessage( WString::tr("warn-no-energy-cal"), WarningWidget::WarningMsgMedium );
  }//if( spec_type == SpecUtils::SpectrumType::Foreground && m_dataMeasurement && !sameSpecFile )
  
  
  if( m_dataMeasurement
     && (spec_type == SpecUtils::SpectrumType::Foreground)
     && !sameSpecFile
     && !m_dataMeasurement->multimedia_data().empty() )
  {
    // Close the multimedia dialog if we open a new spectrum file
    if( m_multimedia )
      programmaticallyCloseMultimediaWindow();
    assert( !m_multimedia );
    
#if( USE_REMOTE_RID )
      if( m_autoRemoteRidResultDialog )
        programaticallyCloseAutoRemoteRidResultDialog();
#endif
    
    bool show_notification = false;
    for( const shared_ptr<const SpecUtils::MultimediaData> &mmd : m_dataMeasurement->multimedia_data() )
    {
      if( mmd && (mmd->data_.size() > 25) )
      {
        //const string mime = Wt::Utils::guessImageMimeTypeData( )
        show_notification = true;
      }
    }//for( const shared_ptr<const SpecUtils::MultimediaData> &mmd : m_dataMeasurement->multimedia_data() )
    
    if( show_notification )
    {
      const bool show = UserPreferences::preferenceValue<bool>( "AutoShowSpecMultimedia", this );
      if( show )
      {
        showMultimedia( SpecUtils::SpectrumType::Foreground );
      }else
      {
        const std::string type = SpecUtils::descriptionText(spec_type);
        WString js = WString(
          "{1}: "
          "<div onclick=\"Wt.emit( $('.specviewer').attr('id'),{name:'miscSignal'}, 'showMultimedia-{2}');try{$(this.parentElement.parentElement).remove();}catch(e){} return false;\" class=\"clearsession\">"
            "<span class=\"clearsessiontxt\">{3}</span>"
          "</div>")
        .arg( WString::tr("info-has-img-txt") )
        .arg( type )
        .arg( WString::tr("info-has-img-btn") )
        ;
        
        m_warnings->addMessageUnsafe( js, WarningWidget::WarningMsgShowOnBoardRiid, 20000 );
      }
    }//if( show_notification )
  }//if( m_dataMeasurement && (spec_type == SpecUtils::SpectrumType::Foreground) )
  

  //Display a notice to the user about how they can select different portions of
  //  passthrough/search-mode data
  /*
  if( spec_type==SpecUtils::SpectrumType::Foreground && !!m_dataMeasurement
      && m_dataMeasurement->passthrough() )
  {
    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", this );
  
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
      passMessage( tip, WarningWidget::WarningMsgInfo );
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
    m_fileManager->displayFile( row, meas, type, true, true, 
            SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
  }catch( std::exception & )
  {
    passMessage( "There was an error loading "
                 + (!!meas ? meas->filename() : string("spectrum file")),
                WarningWidget::WarningMsgHigh );
  }//try / catch
  
}//finishLoadUserFilesystemOpenedFile(...)


void InterSpec::promptUserHowToOpenFile( std::shared_ptr<SpecMeas> meas,
                                             std::shared_ptr<SpectraFileHeader> header )
{
  string filename = header->displayName().toUTF8();
  if( filename.size() > 25 )
  {
    SpecUtils::utf8_limit_str_size( filename, 25 );
    filename += "...";
  }
  
  SimpleDialog *dialog = new SimpleDialog( WString::fromUTF8(filename), WString::tr("prompt-how-open") );
  WPushButton *button = dialog->addButton( WString::tr("Foreground") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                          meas, header, SpecUtils::SpectrumType::Foreground ) );
  button->setFocus( true );
  
  button = dialog->addButton( WString::tr("Background") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                          meas, header, SpecUtils::SpectrumType::Background ) );
  
  button = dialog->addButton( WString::tr("Secondary") );
  button->clicked().connect( boost::bind( &InterSpec::finishLoadUserFilesystemOpenedFile, this,
                                         meas, header, SpecUtils::SpectrumType::SecondForeground ) );
}//void promptUserHowToOpenFile(...)


void InterSpec::userOpenFile( std::shared_ptr<SpecMeas> meas, std::string displayFileName )
{
  assert( meas );
  if( !meas )
    throw runtime_error( "Invalid spectrum file." );
  
  if( !m_fileManager )  //shouldnt ever happen.
    throw runtime_error( "Internal logic error, no valid m_fileManager" );
   
  auto header = std::make_shared<SpectraFileHeader>( true, this );
  header->setFile( displayFileName, meas );
  
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
    return;
  }else
  {
    cout << "Will load file " << displayFileName << " requested to be loaded at "
    << WDateTime::currentDateTime().toString(DATE_TIME_FORMAT_STR)
    << endl;
    
    SpectraFileModel *fileModel = m_fileManager->model();
    const int row = fileModel->addRow( header );
    m_fileManager->displayFile( row, meas, SpecUtils::SpectrumType::Foreground, true, true,
            SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
    return;
  }
}//void InterSpec::userOpenFile( std::shared_ptr<SpecMeas> meas, std::string displayFileName )


bool InterSpec::userOpenFileFromFilesystem( const std::string path, std::string displayFileName  )
{
  try
  {
    if( displayFileName.empty() )
      displayFileName = SpecUtils::filename(path);
    
    assert( m_fileManager );
    if( !m_fileManager )  //shouldnt ever happen.
      throw runtime_error( "Internal logic error, no valid m_fileManager" );
    
    if( !SpecUtils::is_file(path) )
      throw runtime_error( "Could not access file '" + path + "'" );
  
    auto meas = make_shared<SpecMeas>();
    const bool success = meas->load_file( path, SpecUtils::ParserType::Auto, displayFileName );
    
    if( !success )
      throw runtime_error( "Failed to decode file" );
    
    userOpenFile( meas, displayFileName );
    return true;
  }catch( std::exception &e )
  {
    cerr << "Caught exception '" << e.what() << "' when trying to load '"
         << path << "'" << endl;
    SpecMeasManager::displayInvalidFileMsg( displayFileName, e.what() );
    return false;
  }//try / catch
  
  return true;
}//bool userOpenFileFromFilesystem( const std::string filepath )


void InterSpec::handleAppUrl( const std::string &url_encoded_url )
{
  string url = Wt::Utils::urlDecode( url_encoded_url );
    
  // Check if the URL contains "RADDATA://G0/" (dont require start with - incase its a 'emailto:'
  //  URI or something)
  if( SpecUtils::icontains(url, "RADDATA://")
     || SpecUtils::istarts_with(url, "interspec://G0/") )
  {
    m_fileManager->handleSpectrumUrl( std::move(url) );
    return;
  }
  
  string host, path, query_str, frag;
  AppUtils::split_uri( url, host, path, query_str, frag );
  
  if( query_str.empty() && !SpecUtils::iequals_ascii(host, "about") )
    throw runtime_error( "No query string found in URI." );
  
  // I dont think we use the fragment component of URLs anywhere, but maybe we accidentally
  //  included a '#' character somewhere when we shouldnt have.
  if( !frag.empty() )
    query_str += "#" + frag;
  
  //cout << "host='" << host << "' and path='" << path << "' and query_str='" << query_str << "'" << endl;
  
  deleteWelcomeDialog( false );
  deleteEnergyCalPreserveWindow();
  deleteLicenseAndDisclaimersWindow();
  
  if( SpecUtils::iequals_ascii(host,"drf") )
  {
    if( !SpecUtils::iequals_ascii(path,"specify") )
      throw runtime_error( "App 'drf' URL with path '" + path + "' not supported." );
    
    DrfSelect::handle_app_url_drf( query_str );
  }else if( SpecUtils::iequals_ascii(host,"decay") )
  {
    DecayWindow *decay = InterSpec::createDecayInfoWindow();
    if( decay )
      decay->handleAppUrl( path, query_str );
  }else if( SpecUtils::iequals_ascii(host,"specexport") )
  {
    ExportSpecFileWindow *w = createExportSpectrumFileDialog();
    if( w )
      w->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"dose") )
  {
    DoseCalcWindow *dose = showDoseTool();
    if( dose && dose->tool() )
      dose->tool()->handleAppUrl( path, query_str );
  }else if( SpecUtils::iequals_ascii(host,"flux") )
  {
    FluxToolWindow *flux = createFluxTool();
    if( flux )
      flux->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"specsum") )
  {
    GammaCountDialog *gammasum = showGammaCountDialog();
    if( gammasum )
      gammasum->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"gammaxs") )
  {
    GammaXsWindow *xs = showGammaXsTool();
    if( xs && xs->xstool() )
      xs->xstool()->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"1overr2") )
  {
    OneOverR2Calc *calc = createOneOverR2Calculator();
    if( calc )
      calc->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"unit") )
  {
    UnitsConverterTool *converter = createUnitsConverterTool();
    if( converter )
      converter->handleAppUrl( query_str );
  }else if( SpecUtils::iequals_ascii(host,"about") )
  {
    showLicenseAndDisclaimersWindow();
    if( m_licenseWindow )
      m_licenseWindow->handleAppUrlPath( path );
  }else if( SpecUtils::iequals_ascii(host,"welcome") )
  {
    showWelcomeDialog(true);
    if( m_useInfoWindow )
      m_useInfoWindow->handleAppUrl( query_str );
  }
#if( USE_REMOTE_RID )
  else if( SpecUtils::iequals_ascii(host,"remoterid") )
  {
    RemoteRid::handleAppUrl( query_str );
  }
#endif
#if( USE_DETECTION_LIMIT_TOOL )
  else if( SpecUtils::iequals_ascii(host,"detection-limit") )
  {
    showDetectionLimitTool( query_str );
  }else if( SpecUtils::iequals_ascii(host,"simple-mda") )
  {
    string prev_state;
    if( m_simpleMdaWindow )
      prev_state = m_simpleMdaWindow->tool()->encodeStateToUrl();
    
    auto create_window = [this, url](){
      DetectionLimitSimpleWindow *tool = showSimpleMdaWindow();
      if( tool )
        tool->tool()->handleAppUrl( url );
    };//create_window lambda
    
    if( m_undo && m_undo->canAddUndoRedoNow() )
    {
      auto undo = [this,prev_state](){
        if( !prev_state.empty() )
        {
          DetectionLimitSimpleWindow *tool = showSimpleMdaWindow();
          if( tool )
            tool->tool()->handleAppUrl( prev_state );
        }else 
        {
          programmaticallyCloseSimpleMda();
        }
      };//undo lambda
      
      m_undo->addUndoRedoStep( std::move(undo), create_window, "Show Simple MDA URI" );
    }//if( dialog && m_undo && m_undo->canAddUndoRedoNow() )
    
    create_window();
  }
#endif //USE_DETECTION_LIMIT_TOOL
  else
  {
    throw runtime_error( "App URL with purpose (host-component) '" + host + "' not supported." );
  }
}//void handleAppUrl( std::string url )


SimpleDialog *InterSpec::makeEnterAppUrlWindow()
{
  if( m_enterUri )
    return m_enterUri;
  
  m_enterUri = EnterAppUrlWindow::createEntryWindow( this );
  assert( m_enterUri );
  if( m_enterUri )
    m_enterUri->finished().connect( this, &InterSpec::handleAppUrlClosed );
  
  return m_enterUri;
}//SimpleDialog *makeEnterAppUrlWindow()


void InterSpec::handleAppUrlClosed()
{
  m_enterUri = nullptr;
}//void handleAppUrlClosed()


void InterSpec::useMessageResourceBundle( const std::string &name )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
    app->useMessageResourceBundle( name );
}


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
  
  
  if( !m_undo )
    return;
  
  int menu_num_changed = -1;
  bool checked = false;
  WObject *caller = WObject::sender();
  const vector<WMenuItem *> items = m_detectorToShowMenu->items();
  
  for( size_t i = 0; i < items.size(); ++i )
  {
    if( items[i]->hasStyleClass("PhoneMenuBack") )
      continue;
    
    PopupDivMenuItem *item = dynamic_cast<PopupDivMenuItem *>( items[i] );
    WCheckBox *cb = (item ? item->checkBox() : (WCheckBox *)0);
    
    if( cb && ((caller == cb) || (caller == item)) )
    {
      menu_num_changed = static_cast<int>( i );
      assert( cb );
      checked = cb ? cb->isChecked() : false;
      break;
    }
  }//for( size_t i = 0; i < items.size(); ++i )
  
  if( menu_num_changed >= 0 )
  {
    auto undo_redo = [this,menu_num_changed]( const bool set_checked ){
      const vector<WMenuItem *> items = m_detectorToShowMenu->items();
      if( menu_num_changed >= static_cast<int>(items.size()) )
        return;
      
      PopupDivMenuItem *item = dynamic_cast<PopupDivMenuItem *>( items[menu_num_changed] );
      assert( item );
      WCheckBox *cb = (item ? item->checkBox() : (WCheckBox *)0);
      assert( cb );
      if( cb )
      {
        cb->setChecked( set_checked );
        detectorsToDisplayChanged();
      }
    };//auto undo_redo = [...
    
    m_undo->addUndoRedoStep( [undo_redo,checked](){ undo_redo( !checked ); },
                             [undo_redo,checked](){ undo_redo( checked ); },
                            "Toggle display detector" );
  }//if( menu_num_changed >= 0 )
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

    WString title;
    if( name.size() )
      title = WString::tr("app-mi-det-with-name").arg(detnum).arg(name);
    else
      title = WString::tr("app-mi-det-no-name").arg(detnum);
    
    vector<string>::const_iterator neut_name_pos;
    neut_name_pos = std::find( neut_det_names.begin(), neut_det_names.end(), name );
    const int isneut = (neut_name_pos != neut_det_names.end());
    //if( isneut )
      //snprintf( title, sizeof(title), "Det. %i (%s)", i, name.c_str() );

    if( m_detectorToShowMenu )
    {
#if( USE_OSX_NATIVE_MENU )
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
  return m_spectrum->addDecorativeHighlightRegion( lowerEnergy, upperEnergy, color,
                                                  D3SpectrumDisplayDiv::HighlightRegionFill::Full,
                                                  "" );
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


tuple<double,double,Wt::WString> InterSpec::setYAxisRange( double lower_counts, double upper_counts )
{
  return m_spectrum->setYAxisRange( lower_counts, upper_counts );
}//bool setYAxisRange(...)


bool InterSpec::setLogYAxisMin( const double lower_value )
{
  try
  {
    m_spectrum->setLogYAxisMin( lower_value );
  }catch( std::exception & )
  {
    return false;
  }
  return true;
}//setLogYAxisMin( float lower_value )

       
void InterSpec::displayedSpectrumRange( double &xmin, double &xmax, double &ymin, double &ymax ) const
{
  m_spectrum->visibleRange( xmin, xmax, ymin, ymax );
}

void InterSpec::handleShiftAltDrag( double lowEnergy, double upperEnergy )
{
  string prev_uri;
  if( m_gammaCountDialog )
    prev_uri = m_gammaCountDialog->encodeStateToUrl();
  
  {
    //Keep showing dialog and `setEnergyRange` from adding seperate undo states
    UndoRedoManager::BlockUndoRedoInserts blocker;
    
    GammaCountDialog *dialog = showGammaCountDialog();
    if( dialog )
      dialog->setEnergyRange( lowEnergy, upperEnergy );
  }
   
  if( m_undo && m_undo->canAddUndoRedoNow() )
  {
    auto undo = [prev_uri,this](){
      if( prev_uri.empty() )
      {
        deleteGammaCountDialog();
      }else
      {
        GammaCountDialog *dialog = showGammaCountDialog();
        if( dialog )
          dialog->handleAppUrl( prev_uri );
      }
    };
    
    auto redo = [this,lowEnergy,upperEnergy](){
      GammaCountDialog *dialog = showGammaCountDialog();
      if( dialog )
        dialog->setEnergyRange( lowEnergy, upperEnergy );
    };
    
    m_undo->addUndoRedoStep( std::move(undo), std::move(redo), "Set energy range to sum." );
  }//if( m_undo && m_undo->canAddUndoRedoNow() )
  
}//void InterSpec::handleShiftAltDrag( double lowEnergy, double upperEnergy )


void InterSpec::handleSpectrumChartXRangeChange( const double xmin, const double xmax,
                                                 const double oldXmin, const double oldXmax,
                                                 const bool user_interaction )
{
  // Add Undo/Redo step here.
  if( !m_undo || !user_interaction
     || ((fabs(xmin - oldXmin) < 0.5) && (fabs(xmax - oldXmax) < 0.5)) )
    return;
  
  auto undo = [oldXmin, oldXmax, this](){ m_spectrum->setXAxisRange( oldXmin, oldXmax ); };
  auto redo = [xmin, xmax, this](){ m_spectrum->setXAxisRange( xmin, xmax ); };
    
  m_undo->addUndoRedoStep( undo, redo, "Spectrum energy range change." );
}//void handleSpectrumChartXRangeChange(...);


void InterSpec::searchForSinglePeak( const double x )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  if( !m_peakModel )
    throw runtime_error( "InterSpec::searchForSinglePeak(...): "
                        "shoudnt be called if peak model isnt set.");
  
  std::shared_ptr<const SpecUtils::Measurement> data = m_spectrum->data();
  
  if( !m_dataMeasurement || !data )
    return;
  
  const double xmin = m_spectrum->xAxisMinimum();
  const double xmax = m_spectrum->xAxisMaximum();
  
  const double specWidthPx = m_spectrum->chartWidthInPixels();
  const double pixPerKeV = (xmax > xmin && xmax > 0.0 && specWidthPx > 10.0) ? std::max(0.001,(specWidthPx/(xmax - xmin))): 0.001;
  
  std::shared_ptr<const DetectorPeakResponse> det = m_dataMeasurement->detector();
  
  
  PeakSearchGuiUtils::fit_peak_from_double_click( this, x, pixPerKeV, det );
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
  std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > origPeaks
                                                        = m_peakModel->peaks();
  if( !!origPeaks )
    origPeaks = std::make_shared<deque<PeakModel::PeakShrdPtr> >( *origPeaks );
  
  std::shared_ptr< vector<std::shared_ptr<const PeakDef> > > searchresults
            = std::make_shared< vector<std::shared_ptr<const PeakDef> > >();
  
  std::weak_ptr<const SpecUtils::Measurement> weakdata = m_spectrum->data();
  auto drf = data->detector();
  std::weak_ptr<SpecMeas> spectrum = data;
  
  boost::function<void(void)> callback = wApp->bind(
                boost::bind(&InterSpec::setHintPeaks,
                this, spectrum, samples, origPeaks, searchresults) );
  
  boost::function<void(void)> worker = boost::bind( &PeakSearchGuiUtils::search_for_peaks_worker,
                                                   weakdata,
                                                   drf, 
                                                   origPeaks,
                                                   vector<ReferenceLineInfo>(),
                                                   false,
                                                   searchresults,
                                                   callback,
                                                   wApp->sessionId(),
                                                   true );

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


void InterSpec::excludePeaksFromRange( double x0, double x1 )
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
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
        drf_dir = SpecUtils::append_path(ns_staticDataDirectory, "GenericGadrasDetectors/HPGe 40%" );
      else
        drf_dir = SpecUtils::append_path(ns_staticDataDirectory, "GenericGadrasDetectors/NaI 1x1" );
      
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
    
    PopupDivMenuItem *item = dynamic_cast<PopupDivMenuItem *>( items[i] );
    WCheckBox *cb = (item ? item->checkBox() : (WCheckBox *)0);
    
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


void InterSpec::checkEnableViewImageMenuItem()
{
  if( m_showMultimedia )
  {
    const bool showPics = m_dataMeasurement && (m_dataMeasurement->multimedia_data().size() > 0);
    m_showMultimedia->setDisabled( !showPics );
  }
}//void checkEnableViewImageMenuItem()


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
    
    const char *msg = "";
    if( meass.empty() )
    {
      if( detectors.size() == meas->detector_names().size() )
        msg = "err-no-fore-sn";
      else
        msg = "err-no-fore-det-sn";
    }else if( allNeutron )
    {
      msg = "err-no-fore-all-neut";
    }else if( !containSpectrum )
    {
      msg = "err-no-fore-no-spec";
    }else
    {
      msg = (meass.size() > 1) ? "err-no-fore-no-ecal-g1spec" : "err-no-fore-no-ecal-1spec";
    }//if(
    
    passMessage( WString::tr(msg), WarningWidget::WarningMsgHigh );
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
    dataH->set_title( WString::tr("Foreground").toUTF8() );

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
    histH->set_title( WString::tr("second-foreground").toUTF8() );
    
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
    backgroundH->set_title( WString::tr("Background").toUTF8() );
    
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

