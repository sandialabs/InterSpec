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
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WCssDecorationStyle>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/UseInfoWindow.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/LicenseAndDisclaimersWindow.h"

using namespace Wt;
using namespace std;


LicenseAndDisclaimersWindow::LicenseAndDisclaimersWindow( const bool is_awk, int screen_width, int screen_height )
: AuxWindow("Disclaimers, Licenses, Credit, and Contact",
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
               | AuxWindowProperties::DisableCollapse | AuxWindowProperties::EnableResize) ),
  m_menu( nullptr )
{
  setClosable( !is_awk );
  if( !is_awk )
    rejectWhenEscapePressed();
  
  m_resourceBundle.use("InterSpec_resources/static_text/copyright_and_about",false);
  
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
  m_menu->addStyleClass( "SideMenu" );
  
  WDialog::contents()->setOverflow( WContainerWidget::OverflowHidden );
  
  //If on phone, need to make text much smaller!
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  const bool phone = (app && app->isPhone());
  const bool tablet = (app && app->isTablet());
  if( phone )
    WDialog::contents()->addStyleClass( "PhoneCopywriteContent" );
  
  WContainerWidget *topDiv = new WContainerWidget();
  WContainerWidget *bottomDiv = new WContainerWidget();
  
  WBorder border(WBorder::Solid, WBorder::Explicit, Wt::gray);
  border.setWidth( WBorder::Explicit, WLength(1) );
  
  topDiv->decorationStyle().setBorder( border,  Wt::Bottom );
  bottomDiv->decorationStyle().setBorder( border, Wt::Top );
  stack->decorationStyle().setBorder( border,  Wt::Right | Wt::Left );
  m_menu->decorationStyle().setBorder( border,  Wt::Left );
  
  WGridLayout *layout = stretcher();
  
  layout->addWidget( topDiv,    0, 0, 1, 2 );
  layout->addWidget( m_menu,    1, 0 );
  layout->addWidget( stack,     1, 1, 1, 1 );
  layout->addWidget( bottomDiv, 2, 0, 1, 2 );
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
  
  //Populate bottomDiv
  string dhsack;
  m_resourceBundle.resolveKey( "dhs-acknowledgement", dhsack );
  new WText( dhsack, Wt::XHTMLText, bottomDiv );
  bottomDiv->setAttributeValue( "style", "text-align: center; font-size: smaller; padding-top: 0.6em;" );
  
  //Add items to the left menu; the contents wont be loaded until shown.
  makeItem( "Disclaimer", "dhs-disclaimer" );
  makeLgplLicenseItem();
  makeItem( "Credits", "credits" );
  makeItem( "Contact", "contact" );
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__))) )
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
  std::function<void(WContainerWidget *)> f = boost::bind( &LicenseAndDisclaimersWindow::itemCreator, this, resource, _1, title );
  
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
  const string approot = wApp->appRoot();  //looks to be blank if CWD
  const string license_file = UtilityFunctions::append_path(approot, "InterSpec_resources/static_text/lgpl-2.1.txt" );
  
  bool error_reading = false;
  string license_content;
  
  try
  {
#ifdef _WIN32
    const std::wstring wlicense_file = UtilityFunctions::convert_from_utf8_to_utf16(license_file);
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
    text->setAttributeValue( "style", "display: block; font-family: monospace; white-space: pre; margin: 1em 0;" );
}//void lgplLicenseCreator()


SideMenuItem *LicenseAndDisclaimersWindow::makeLgplLicenseItem()
{
  std::function<void(WContainerWidget *)> f = boost::bind( &LicenseAndDisclaimersWindow::lgplLicenseCreator, this, _1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( "License", w );
  
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem *makeItem( const WString &title, const string &resource)

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__))) )
void LicenseAndDisclaimersWindow::dataStorageCreator( Wt::WContainerWidget *parent )
{
  auto app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  
  InterSpec *viewer = app ? app->viewer() : nullptr;
  
  string contents;
  if( viewer )
  {
    string datadir = Wt::Utils::htmlEncode( viewer->writableDataDirectory() );
    string staticdir = Wt::Utils::htmlEncode( viewer->staticDataDirectory() );
    
    //UtilityFunctions::is_absolute_path( staticdir )
    try
    {
      UtilityFunctions::make_canonical_path( datadir );
      UtilityFunctions::make_canonical_path( staticdir );
    }catch( std::exception &e )
    {
      cerr << "Failed to make_canonical_path for either '" << datadir
           << "' or '" << staticdir << "': " << e.what() << endl;
    }
    
    
#ifdef _WIN32
    UtilityFunctions::ireplace_all( datadir, "/", "\\" );
    UtilityFunctions::ireplace_all( staticdir, "/", "\\" );
#endif
    
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
    "</p>";
  }else
  {
    contents = "Error retrieving directory data";
  }
  
  WText *text = new WText( WString::fromUTF8(contents), XHTMLText, parent );
  text->setAttributeValue( "style", "display: block; margin: 1em 0;" );
}//void dataStorageCreator( Wt::WContainerWidget *parent );


SideMenuItem *LicenseAndDisclaimersWindow::makeDataStorageItem()
{
  std::function<void(WContainerWidget *)> f = boost::bind( &LicenseAndDisclaimersWindow::dataStorageCreator, this, _1 );
  
  WWidget *w = deferCreate( f );
  w->addStyleClass( "UseInfoItem" );
  
  SideMenuItem *item = new SideMenuItem( "Data", w );
  
  item->clicked().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  item->mouseWentDown().connect( boost::bind( &LicenseAndDisclaimersWindow::right_select_item, this, item) );
  
  m_menu->addItem( item );
  
  return item;
}//SideMenuItem *makeDataStorageItem()
#endif //#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || (BUILD_AS_LOCAL_SERVER && (defined(WIN32) || defined(__APPLE__))) )
