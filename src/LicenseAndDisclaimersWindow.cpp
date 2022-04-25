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

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
#include <Wt/WServer>
#endif

#include "SpecUtils/DateTime.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UseInfoWindow.h"
#include "InterSpec/LicenseAndDisclaimersWindow.h"

#ifdef _WIN32
#include "SpecUtils/StringAlgo.h"
#endif

using namespace Wt;
using namespace std;



//The below YEAR MONTH DAY macros are taken from
//http://bytes.com/topic/c/answers/215378-convert-__date__-unsigned-int
//  and I believe to be public domain code
#define YEAR ((((__DATE__ [7] - '0') * 10 + (__DATE__ [8] - '0')) * 10 \
+ (__DATE__ [9] - '0')) * 10 + (__DATE__ [10] - '0'))
#define MONTH (__DATE__ [2] == 'n' && __DATE__ [1] == 'a' ? 0 \
: __DATE__ [2] == 'b' ? 1 \
: __DATE__ [2] == 'r' ? (__DATE__ [0] == 'M' ? 2 : 3) \
: __DATE__ [2] == 'y' ? 4 \
: __DATE__ [2] == 'n' ? 5 \
: __DATE__ [2] == 'l' ? 6 \
: __DATE__ [2] == 'g' ? 7 \
: __DATE__ [2] == 'p' ? 8 \
: __DATE__ [2] == 't' ? 9 \
: __DATE__ [2] == 'v' ? 10 : 11)
#define DAY ((__DATE__ [4] == ' ' ? 0 : __DATE__ [4] - '0') * 10 \
+ (__DATE__ [5] - '0'))

/** \brief Macro to embed compile date (as an decimal int) in ISO format into
 the program
 
 example: Will evaluate to 20120122 for jan 22nd, 2012.
 */
#define COMPILE_DATE_AS_INT (YEAR*10000 + (MONTH+1)*100 + DAY)

//Or could just use: __DATE__


LicenseAndDisclaimersWindow::LicenseAndDisclaimersWindow( const bool is_awk, int screen_width, int screen_height )
: AuxWindow("Disclaimers, Licenses, Credit, and Contact",
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
               | AuxWindowProperties::DisableCollapse | AuxWindowProperties::EnableResize) ),
  m_menu( nullptr )
{
  setClosable( !is_awk );
  if( !is_awk )
    rejectWhenEscapePressed();
  
  const string docroot = wApp->docRoot();
  m_resourceBundle.use( SpecUtils::append_path(docroot,"InterSpec_resources/static_text/copyright_and_about") ,false);
  
  double width = 0.5*screen_width, height = 0.8*screen_height;

  if( height < 512.0 )
    height = 1.0*std::min( screen_height, 512 );
  height = std::min( height, 1024.0 );  //1024 not actually tested, could maybye bee 800

  if( width < 715.0 || height < 512.0 )
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
  m_menu->addStyleClass( "VerticalNavMenu HeavyNavMenu SideMenu" );
  
  WDialog::contents()->setOverflow( WContainerWidget::OverflowHidden );
  
  //If on phone, need to make text much smaller!
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  const bool phone = (app && app->isPhone());
  const bool tablet = (app && app->isTablet());
  if( phone )
    WDialog::contents()->addStyleClass( "PhoneCopywriteContent" );
  
  WContainerWidget *topDiv = new WContainerWidget();
  
  WBorder border(WBorder::Solid, WBorder::Explicit, Wt::gray);
  border.setWidth( WBorder::Explicit, WLength(1) );
  
  topDiv->decorationStyle().setBorder( border,  Wt::Bottom );
  stack->decorationStyle().setBorder( border,  Wt::Right | Wt::Left );
  m_menu->decorationStyle().setBorder( border,  Wt::Left );
  
  WGridLayout *layout = stretcher();
  
  layout->addWidget( topDiv,    0, 0, 1, 2 );
  layout->addWidget( m_menu,    1, 0 );
  layout->addWidget( stack,     1, 1, 1, 1 );
  layout->setRowStretch( 1, 1 );
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
  title->bindString("build-date", std::to_string(COMPILE_DATE_AS_INT) );
  title->bindString("copyright", copyright );
  
  //Add items to the left menu; the contents wont be loaded until shown.
  makeItem( "Disclaimer", "dhs-disclaimer" );
  makeLgplLicenseItem();
  makeItem( "Credits", "credits" );
  makeItem( "Contact", "contact" );
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  makeDataStorageItem();
#endif
  m_menu->select( 0 );
  
  layout->setColumnStretch( 1, 1 );

  //Put in a footer with an Acknowledge that will accept this dialog
  WContainerWidget *bottom = footer();
  
  
  if( phone || (tablet && !is_awk) )
  {
    WPushButton *close = addCloseButtonToFooter();
    close->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  }else
  {
    WPushButton *close = new WPushButton( (is_awk ? "Acknowledge" : "Close"), bottom );
    close->addStyleClass( "CenterBtnInMblAuxWindowHeader" );
    close->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  }
  
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



void LicenseAndDisclaimersWindow::right_select_item(  WMenuItem *item )
{
  m_menu->select( item );
  item->triggered().emit( item ); //doenst look like this is emmitted either
                                  //when body of SideMenuItem is clicked
                                  //stop all players
}//void LicenseAndDisclaimersWindow::select_item(  SideMenuItem *item )


SideMenuItem * LicenseAndDisclaimersWindow::makeItem( const WString &title, const string &resource)
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
    text->setAttributeValue( "style", "display: block; font-family: monospace; white-space: pre; margin: 1em 0; font-size: x-small;" );
}//void lgplLicenseCreator()


SideMenuItem *LicenseAndDisclaimersWindow::makeLgplLicenseItem()
{
  std::function<void(WContainerWidget *)> f
              = boost::bind( &LicenseAndDisclaimersWindow::lgplLicenseCreator, this,
                            boost::placeholders::_1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( "License", w );
  
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem *makeItem( const WString &title, const string &resource)

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
void LicenseAndDisclaimersWindow::dataStorageCreator( Wt::WContainerWidget *parent )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  
  InterSpec *viewer = app ? app->viewer() : nullptr;
  
  string contents;
  if( viewer )
  {
    string datadir;
    try{ datadir = Wt::Utils::htmlEncode( viewer->writableDataDirectory() ); }catch(...){}
    string staticdir = Wt::Utils::htmlEncode( viewer->staticDataDirectory() );
    
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
    Dbo::ptr<InterSpecUser> user = ((app && viewer) ? viewer->m_user : Dbo::ptr<InterSpecUser>());
    if( user )
    {
      try
      {
        totalSessions = std::to_string( user->accessCount() );
        totalFilesOpened = std::to_string( user->numSpectraFilesOpened() );
        boost::posix_time::time_duration totaltime = user->totalTimeInApp();
        // Note that if user has multiple sessions going, this next line wont be exactly correct, but close enough.
        totaltime += app->activeTimeInCurrentSession();
        totalUserTime = PhysicalUnits::printToBestTimeUnits( totaltime.total_seconds() );
        
        const WDateTime utcStartTime = WDateTime::fromPosixTime( user->firstAccessUTC() );
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
    
    const string style = "font-family: monospace;"
                         " background: white;"
                         " color: black;"
                         " white-space: nowrap;"
                         " padding: 4px 5px 4px 10px;"
                         " overflow-x: auto;"
                         " -webkit-user-select: all;"
                         " user-select: all;";
    
    contents =
    "<p>User data is stored in"
    "<div style=\"" + style + "\">" + datadir + "</div>"
    "</p>"
    "<p>The data that comes with InterSpec, such as nuclear decay info,"
    " cross-section, and similar is stored in"
    "<div style=\"" + style + "\">" + staticdir + "</div>"
    "</p>"
    
    "<p>You can also use InterSpec from your browser at"
    " (port number will change when you restart InterSpec):"
    "<div style=\"" + style + "\">http://127.0.0.1:" + std::to_string(httpPort) + "</div>"
    "</p>"
    ;
    
#if( USE_DB_TO_STORE_SPECTRA )
    if( displayStats )
    {
      contents += "<div style=\"margin-top: 10px\"><p>"
      "You have actively used InterSpec for approximately "
      + totalUserTime + ", to open " + totalFilesOpened
      + " files, over " + totalSessions + " sessions since " + firstAccess + "."
      "</p></div>";
    }//if( displayStats )
#endif //if( USE_DB_TO_STORE_SPECTRA )
    
    
#if( BUILD_FOR_WEB_DEPLOYMENT )
    // Statement here would depend on web-server policies
#else
    contents +=
    "<div style=\"margin-top: 10px\"><p>"
    "InterSpec does not send or receive data external to your device"
#if( USE_GOOGLE_MAP )
    ", except when the Google Maps feature is used"
#endif
    "."
    
#if( USE_DB_TO_STORE_SPECTRA )
    "<br />InterSpec does locally"
    " store preferences, spectra you load, saved app states, and use"
    " statistics."
    " This information does not leave your device, and can be deleted by removing the user data"
    " directory shown above."
#endif
    
#if( USE_GOOGLE_MAP )
    "<br />When the Google Maps feature is used, the spectrum file location information is sent to"
    " Google in order to receive maps of the relevant location."
#endif
    "</p></div>";
    
#endif //#if( BUILD_FOR_WEB_DEPLOYMENT ) / else

    
    
    
  }else
  {
    contents = "Error retrieving directory data";
  }
  
  WText *text = new WText( WString::fromUTF8(contents), XHTMLText, parent );
  text->setAttributeValue( "style", "display: block; margin: 1em 0;" );
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
#endif //#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
