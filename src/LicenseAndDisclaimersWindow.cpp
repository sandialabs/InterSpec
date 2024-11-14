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
#include <fstream>

#include <Wt/WText>
#include <Wt/WMenu>
#include <Wt/Utils>
#include <Wt/WBorder>
#include <Wt/WString>
#include <Wt/WTemplate>
#include <Wt/WMenuItem>
#include <Wt/WAnimation>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WCssDecorationStyle>

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
#include <Wt/WServer>
#endif

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/LicenseAndDisclaimersWindow.h"

#ifdef _WIN32
#include "SpecUtils/StringAlgo.h"
#endif

using namespace Wt;
using namespace std;

LicenseAndDisclaimersWindow::LicenseAndDisclaimersWindow( InterSpec *interspec )
: AuxWindow( WString::tr("window-title-license-credit"),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
             | AuxWindowProperties::DisableCollapse
             | AuxWindowProperties::EnableResize
             | AuxWindowProperties::SetCloseable) ),
  m_menu( nullptr )
{
  setClosable( true );
  rejectWhenEscapePressed();
  
  assert( interspec );
  if( interspec )
    interspec->useMessageResourceBundle( "LicenseAndDisclaimersWindow" );
  
  wApp->useStyleSheet( "InterSpec_resources/LicenseAndDisclaimersWindow.css" );
  
  addStyleClass( "LicenseAndDisclaimersWindow" );
  
  const string docroot = wApp->docRoot();
  m_resourceBundle.use( SpecUtils::append_path(docroot,"InterSpec_resources/static_text/copyright_and_about") ,false);
  
  int screen_width = interspec ? interspec->renderedWidth() : 0;
  int screen_height = interspec ? interspec->renderedHeight() : 0;
  if( (screen_width < 100) && interspec && interspec->isMobile() )
  {
    screen_width = wApp->environment().screenWidth();
    screen_height = wApp->environment().screenHeight();
  }
  
  const bool narrow_layout = ((screen_width > 100) 
                              && (screen_width < 500)
                              && (screen_width < screen_height));
  
  double width = 0.5*screen_width;
  double height = 0.8*screen_height;

  if( height < 512.0 )
    height = 1.0*std::min( screen_height, 512 );
  height = std::min( height, 1024.0 );  //1024 not actually tested, could maybye bee 800

  if( !narrow_layout && ((width < 715.0) || (height < 512.0)) )
  {
    setMinimumSize(715,512);
    resize( WLength(50, WLength::FontEm), WLength(80,WLength::Percentage));
  }else
  {
    resize( WLength(width), WLength(height) );
  }
  
  WStackedWidget *stack = new WStackedWidget();
  stack->addStyleClass( "UseInfoStack" );
  stack->setOverflow( WContainerWidget::OverflowAuto );
  
  WAnimation animation(Wt::WAnimation::Fade, Wt::WAnimation::Linear, 200);
  stack->setTransitionAnimation( animation, true );
  
  m_menu = new WMenu( stack, Wt::Vertical );
  if( narrow_layout )
    m_menu->addStyleClass( "VerticalNavMenu HeavyNavMenu HorizontalMenu" );
  else
    m_menu->addStyleClass( "VerticalNavMenu HeavyNavMenu SideMenu" );
  
  //HorizontalMenu
  
  WDialog::contents()->setOverflow( WContainerWidget::OverflowHidden );
  
  //If on phone, need to make text much smaller!
  const bool phone = (interspec && interspec->isPhone());
  //const bool tablet = (app && app->isTablet());
  if( phone )
    WDialog::contents()->addStyleClass( "PhoneCopywriteContent" );
  
  WContainerWidget *topDiv = new WContainerWidget();
  
  WBorder border(WBorder::Solid, WBorder::Explicit, Wt::gray);
  border.setWidth( WBorder::Explicit, WLength(1) );
  
  WGridLayout *layout = stretcher();
  
  if( narrow_layout )
  {
    m_menu->setMargin( 5, Wt::Side::Bottom );
    
    layout->addWidget( topDiv,    0, 0 );
    layout->addWidget( m_menu,    1, 0 );
    layout->addWidget( stack,     2, 0 );
    layout->setRowStretch( 2, 1 );
  }else
  {
    topDiv->decorationStyle().setBorder( border,  Wt::Bottom );
    stack->decorationStyle().setBorder( border,  Wt::Right | Wt::Left );
    m_menu->decorationStyle().setBorder( border,  Wt::Left );
    
    layout->addWidget( topDiv,    0, 0, 1, 2 );
    layout->addWidget( m_menu,    1, 0 );
    layout->addWidget( stack,     1, 1, 1, 1 );
    layout->setRowStretch( 1, 1 );
  }
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  
  
  //Populate topDiv
  //The SNL/NTESS copywrite must appear before any open source software licenses
  string apptitle, copyright;
  m_resourceBundle.resolveKey( "app-title", apptitle );
  m_resourceBundle.resolveKey( "copyright", copyright );
  
  WTemplate *title = new WTemplate( topDiv );
  title->setTemplateText( apptitle );
  title->bindString("build-version", InterSpec_VERSION);
  title->bindString("build-date", std::to_string(AppUtils::compile_date_as_int()) );
  title->bindString("copyright", copyright );
  
  //Add items to the left menu; the contents wont be loaded until shown.
  makeItem( WString::tr("ladw-mi-disclaimer"), "dhs-disclaimer" );
  makeLgplLicenseItem();
  makeItem( WString::tr("ladw-mi-credits"), "credits" );
  makeItem( WString::tr("ladw-mi-contact"), "contact" );
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  makeDataStorageItem();
#endif
  m_menu->select( 0 );
  
  layout->setColumnStretch( 1, 1 );

  //Put in a footer with an Acknowledge that will accept this dialog
  WPushButton *close = addCloseButtonToFooter();
  close->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  
  resizeToFitOnScreen(); //jic
  centerWindow();
}//LicenseAndDisclaimersWindow( InterSpec* viewer ):



LicenseAndDisclaimersWindow::~LicenseAndDisclaimersWindow()
{
}//~LicenseAndDisclaimersWindow()





void LicenseAndDisclaimersWindow::itemCreator( const string &resource, Wt::WContainerWidget *parent, WString title)
{
  std::string resourceContent;
  m_resourceBundle.resolveKey( resource, resourceContent );
  WTemplate *templ = new WTemplate( parent );
  templ->setTemplateText( WString(resourceContent, UTF8), XHTMLUnsafeText );
}//void LicenseAndDisclaimersWindow::itemCreator( const string &resource, Wt::WContainerWidget *parent, WString title)



void LicenseAndDisclaimersWindow::right_select_item( WMenuItem *item )
{
  if( !item )
    return;
  
  m_menu->select( item );
  item->triggered().emit( item ); //doenst look like this is emmitted either
                                  //when body of SideMenuItem is clicked
                                  //stop all players
}//void LicenseAndDisclaimersWindow::select_item(  SideMenuItem *item )


bool LicenseAndDisclaimersWindow::handleAppUrlPath( const std::string &path )
{
  int index = -1;
  if( SpecUtils::istarts_with(path, "discl") )
    index = 0;
  else if( SpecUtils::istarts_with(path, "lic") )
    index = 1;
  else if( SpecUtils::istarts_with(path, "cred") )
    index = 2;
  else if( SpecUtils::istarts_with(path, "cont") )
    index = 3;
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  else if( SpecUtils::istarts_with(path, "data") )
    index = 4;
#endif
  
  if( (index < 0) || (index >= m_menu->count()) )
  {
    passMessage( "Appropriate URI path for '" + Wt::Utils::htmlEncode(path) + "' could not be found.", 3 );
    return false;
  }
  
  right_select_item( m_menu->itemAt(index) );
  return true;
}//bool handleAppUrlPath( const std::string &url )



SideMenuItem *LicenseAndDisclaimersWindow::makeItem( const WString &title, const string &resource)
{
  std::function<void(WContainerWidget *)> f = boost::bind( &LicenseAndDisclaimersWindow::itemCreator,
                                                          this, resource, boost::placeholders::_1,
                                                          title );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( title, w );
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem * LicenseAndDisclaimersWindow::makeItem( const WString &title, const string &resource)


void LicenseAndDisclaimersWindow::lgplLicenseCreator( Wt::WContainerWidget *parent )
{
  //const string approot = wApp->appRoot();  //looks to be blank if CWD
  const string docroot = wApp->docRoot();
  const string license_file = SpecUtils::append_path(docroot, "InterSpec_resources/static_text/lgpl-2.1.txt" );
  
  bool error_reading = false;
  string license_content;
  
  try
  {
#ifdef _WIN32
    const std::wstring wlicense_file = SpecUtils::convert_from_utf8_to_utf16(license_file);
    ifstream stream( wlicense_file.c_str(), ios::in | ios::binary );
#else
    ifstream stream( license_file.c_str(), ios::in | ios::binary );
#endif
    if( !stream )
      throw runtime_error( "Cannot open license file " + license_file );
    
    // Determine file size
    stream.seekg(0, ios::end);
    const size_t file_size = stream.tellg();
    stream.seekg(0);
    
    //The actual file is 26530 bytes - maybe this should just be hardcoded in
    if( file_size > 32768 || file_size < 16384 )
      throw runtime_error( "License file ('" + license_file + "') does not have"
                           " reasonable size contents ("
                           + std::to_string(file_size) + " bytes)." );
    
    license_content.resize( file_size + 1, '\0' );
    if( !stream.read( &(license_content[0]), static_cast<streamsize>(file_size) ) )
      throw runtime_error( "Error reading license file ('" + license_file + "')." );
    license_content[file_size] = '\0';  //should already be zero, but jic
  }catch( std::exception &e )
  {
    //Make sure that even
    error_reading = true;
    license_content = e.what();
    license_content = "<div>" + license_content + "</div>";
    license_content += "<br /><div>InterSpec is licensed under GNU Lesser General Public License, version 2.1 (LGPL v2.1).<br />"
    "A copy of the license should have been included with the software, or it can be obtained by visiting "
    "<a target=\"_blank\" href=\"https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html\">"
    "https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html</a>"
    " or contacting this softwares authors.</div>";
  }//try catch
  
  WText *text = new WText( license_content, XHTMLText, parent );
  
  if( !error_reading )
    text->addStyleClass( "LicenseContent" );
}//void lgplLicenseCreator()


SideMenuItem *LicenseAndDisclaimersWindow::makeLgplLicenseItem()
{
  std::function<void(WContainerWidget *)> f
              = boost::bind( &LicenseAndDisclaimersWindow::lgplLicenseCreator, this,
                            boost::placeholders::_1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( WString::tr("ladw-mi-license"), w );
  
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem *makeItem( const WString &title, const string &resource)

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
void LicenseAndDisclaimersWindow::dataStorageCreator( Wt::WContainerWidget *parent )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  
  InterSpec *viewer = app ? app->viewer() : nullptr;
  
  if( viewer )
  {
    string datadir, userDataDir;
    try
    {
      userDataDir = viewer->writableDataDirectory();
      datadir = Wt::Utils::htmlEncode( userDataDir );
    }catch(...)
    {
      
    }
    
    const string staticDataDir = viewer->staticDataDirectory();
    string staticdir = Wt::Utils::htmlEncode( staticDataDir );
    
    //SpecUtils::is_absolute_path( staticdir )
    try
    {
      SpecUtils::make_canonical_path( datadir );
      SpecUtils::make_canonical_path( staticdir );
    }catch( std::exception &e )
    {
      cerr << "Failed to make_canonical_path for either '" << datadir
           << "' or '" << staticdir << "': " << e.what() << endl;
    }
    
    
#ifdef _WIN32
    SpecUtils::ireplace_all( datadir, "/", "\\" );
    SpecUtils::ireplace_all( staticdir, "/", "\\" );
#endif
    
    
#if( USE_DB_TO_STORE_SPECTRA )
    bool displayStats = false;
    string totalUserTime, totalFilesOpened, totalSessions, firstAccess;
    
    InterSpec *viewer = InterSpec::instance();
    Dbo::ptr<InterSpecUser> user = ((app && viewer) ? viewer->user() : Dbo::ptr<InterSpecUser>());
    if( user )
    {
      try
      {
        totalSessions = std::to_string( user->accessCount() );
        totalFilesOpened = std::to_string( user->numSpectraFilesOpened() );
        chrono::steady_clock::time_point::duration totaltime = user->totalTimeInApp();
        totaltime += app->timeSinceTotalUseTimeUpdated();
        const chrono::seconds numsecs = chrono::duration_cast<chrono::seconds>(totaltime);
        
        totalUserTime = PhysicalUnitsLocalized::printToBestTimeUnits( numsecs.count(), 2, 1.0 );
        
        
        const WDateTime utcStartTime = WDateTime::fromPosixTime( to_ptime(user->firstAccessUTC()) );
        const int offset = app->environment().timeZoneOffset();
        firstAccess = utcStartTime.addSecs(60*offset).toString( "dd-MMM-yyyy" ).toUTF8();
        
        // The alternative of using WLocalDateTime leaves the time in UTC...
        //const WLocalDateTime localStartTime = utcStartTime.toLocalTime();
        //firstAccess = localStartTime.toString("dd-MMM-yyyy").toUTF8();
        
        displayStats = true;
      }catch(...)
      {
        // Dont expect this to ever really happen.
      }//try / catch
    }//if( app && viewer && viewer->m_user )
#endif
    
    auto server = WServer::instance();
    const int httpPort = server ? server->httpPort() : 0;
    
    WText *text = new WText( WString::tr("ladw-data-location-user").arg(datadir), parent );
    text->addStyleClass( "DataLocationSection" );

#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    if( !userDataDir.empty() )
    {
      WPushButton *showBtn = new WPushButton( displOptLower );
#ifdef _WIN32
      const char *txt_key = "ladw-show-data-location-win";
#elif __APPLE__
      const char *txt_key = "ladw-show-data-location-macOS";
#else
      const char *txt_key = "ladw-show-data-location";
#endif
      showBtn->setText( WString::tr(txt_key) );
      showBtn->setStyleClass( "LinkBtn ShowDataLocationBtn" );
      showBtn->clicked().connect( std::bind([userDataDir](){
        AppUtils::showFileInOsFileBrowser(userDataDir);
      }) );
    }
#endif
    
    text = new WText( WString::tr("ladw-data-location-static").arg(staticdir), parent );
    text->addStyleClass( "DataLocationSection" );
    
#if( !ANDROID && !IOS && !BUILD_FOR_WEB_DEPLOYMENT )
    if( !userDataDir.empty() )
    {
      WPushButton *showBtn = new WPushButton( displOptLower );
#ifdef _WIN32
      const char *txt_key = "ladw-show-data-location-win";
#elif __APPLE__
      const char *txt_key = "ladw-show-data-location-macOS";
#else
      const char *txt_key = "ladw-show-data-location";
#endif
      showBtn->setText( WString::tr(txt_key) );
      showBtn->setStyleClass( "LinkBtn ShowDataLocationBtn" );
      showBtn->clicked().connect( std::bind([staticDataDir](){
        AppUtils::showFileInOsFileBrowser(staticDataDir);
      }) );
    }
#endif
    
    text = new WText( WString::tr("ladw-data-location-network").arg(httpPort), parent );
    text->addStyleClass( "DataLocationSection" );
    
#if( USE_DB_TO_STORE_SPECTRA )
    if( displayStats )
    {
      text = new WText( WString::tr("ladw-user-stats")
                          .arg(totalUserTime)
                          .arg(totalFilesOpened)
                          .arg(totalSessions)
                          .arg(firstAccess), parent );
      text->addStyleClass( "DataLocationSection" );
    }//if( displayStats )
#endif //if( USE_DB_TO_STORE_SPECTRA )
    
    WString content("{1}{2}{3}{4}");
    
#if( BUILD_FOR_WEB_DEPLOYMENT )
    // Statement here would depend on web-server policies
    content.arg( "please contact web-administrator." ).arg( "" ).arg( "" );
#else
    
#if( USE_LEAFLET_MAP  || USE_GOOGLE_MAP )
  #if( USE_REMOTE_RID )
    content.arg( WString::tr("ladw-external-remote-rid-and-map") );
  #else
    content.arg( WString::tr("ladw-external-only-map") );
  #endif
#elif( USE_REMOTE_RID )
    content.arg( WString::tr("ladw-external-only-remote-rid") );
#else
    content.arg( WString::tr("ladw-external-none") );
#endif
    
#if( USE_DB_TO_STORE_SPECTRA )
    content.arg( WString::tr("ladw-local-store-info") );
#else
    content.arg( "" );
#endif
    
#if( USE_GOOGLE_MAP )
    content.arg( WString::tr("ladw-google-map") );
#elif( USE_LEAFLET_MAP )
    content.arg( WString::tr("ladw-leaflet-map") );
#else
    content.arg( "" );
#endif
    
#if( USE_REMOTE_RID )
    const string urls = UserPreferences::preferenceValue<string>( "ExternalRidUrl", viewer );
    const string exes = UserPreferences::preferenceValue<string>( "ExternalRidExe", viewer );
    if( urls.empty() && exes.empty() )
      content.arg( WString::tr("ladw-remote-rid-none") );
    
    if( !urls.empty() )
      content.arg( WString::tr("ladw-remote-rid").arg(urls) );
    
    if( !exes.empty() )
      content.arg( WString::tr("ladw-remote-rid").arg(exes) );
#else
    content.arg( "" );
#endif
#endif //#if( BUILD_FOR_WEB_DEPLOYMENT ) / else
    
    text = new WText( content, parent );
    text->addStyleClass( "DataLocationSection" );
  }else
  {
    WText *text = new WText( "Error retrieving directory data", parent );
    text->addStyleClass( "DataLocationSection" );
  }
}//void dataStorageCreator( Wt::WContainerWidget *parent );


SideMenuItem *LicenseAndDisclaimersWindow::makeDataStorageItem()
{
  std::function<void(WContainerWidget *)> f
         = boost::bind( &LicenseAndDisclaimersWindow::dataStorageCreator, this,
                        boost::placeholders::_1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( "Data", w );
  
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem *makeDataStorageItem()
#endif //#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
