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

#include <fstream>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WAnchor>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WResource>
#include <Wt/WMenuItem>
#include <Wt/WFileUpload>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WContainerWidget>

#include <Wt/Dbo/Exception>

#include "InterSpec/InterSpec.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/ColorThemeWidget.h"
#include "InterSpec/ColorThemeWindow.h"

using namespace std;
using namespace Wt;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

namespace
{
  //Equiv to SideMenuItem in UseInfoWindow.h, but reimplemented here for potential styling
  class ThemeMenuItem : public Wt::WMenuItem
  {
  public:
    ThemeMenuItem( const Wt::WString &text,
                  std::unique_ptr<ColorTheme> theme,
                  const bool editable )
    : WMenuItem(text, nullptr, WMenuItem::LazyLoading),
      m_editable(editable), m_isEditedSinceApply(false), m_theme( std::move(theme) )
    {
    }
    
    const ColorTheme *theme() { return m_theme.get(); }
    bool editable() const { return m_editable;  }
    
    bool isEditedSinceApply() const { return m_isEditedSinceApply; }
    void setEditedSinceApply( const bool edited ) { m_isEditedSinceApply = edited; }
    
    void setTheme( std::unique_ptr<ColorTheme> theme )
    {
      m_theme = std::move( theme );
    }
  protected:
    
    const bool m_editable;
    bool m_isEditedSinceApply;
    std::unique_ptr<ColorTheme> m_theme;
  };//class ThemeMenuItem
  
  
  class JsonDownloadResource : public Wt::WResource
  {
    ColorThemeWindow *m_display;
    Wt::WApplication *m_app; //it looks like WApplication::instance() will be valid in handleRequest, but JIC
  public:
    JsonDownloadResource( ColorThemeWindow *parent )
    : WResource( parent ), m_display( parent ), m_app( WApplication::instance() )
    {
      assert( m_app );
    }
    
    virtual ~JsonDownloadResource()
    {
      beingDeleted();
    }
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      WApplication::UpdateLock lock( m_app );
      
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in JsonDownloadResource.";
        
        response.out() << "Error grabbing application lock to form JsonDownloadResource resource; please report to InterSpec@sandia.gov.";
        response.setStatus(500);
        assert( 0 );
        
        return;
      }//if( !lock )
      
      if( !m_display )
        return;
      string name = m_display->currentThemeTitle();
      if( name.empty() )
        name = "empty";
      
      //Remove bad filename characters
      const string notallowed = "\\/:?\"<>|*";
      for( auto it = begin(name) ; it < end(name) ; ++it )
      {
        if( notallowed.find(*it) != string::npos )
          *it = ' ';
      }
      
      suggestFileName( name + ".json", WResource::Attachment );
      response.setMimeType( "application/json" );
      response.out() << m_display->currentThemeJson();
    }
  };//class JsonDownloadResource : public Wt::WResource
  
}//namespace



ColorThemeWindow::ColorThemeWindow( InterSpec *interspec )
: AuxWindow( "Color Themes", Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen) | AuxWindowProperties::SetCloseable ),
m_interspec(interspec),
m_menu( nullptr ),
m_edit( nullptr ),
m_removeIcn( nullptr ),
m_addIcn( nullptr ),
m_close( nullptr ),
m_save( nullptr ),
m_apply( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/ColorThemeWindow.css" );
  
  int w = m_interspec->renderedWidth();
  int h = m_interspec->renderedHeight();
  if( w > 0 && h > 0 )
  {
    const int ww = std::min( std::max( w/2, 750 ), (9*w)/10 );
    const int wh = std::min( std::min( std::max( (6*h)/7, 420), (9*h)/10 ), 1024);
    resizeWindow( ww, wh); //Should have no effect on phones
  }else
  {
    resizeScaledWindow(0.9, 0.95);
  }
  
  WContainerWidget *leftMenuDiv = new WContainerWidget();
  leftMenuDiv->addStyleClass( "VerticalColorThemesMenu" );
  
  WContainerWidget *editDiv = new WContainerWidget();
  editDiv->addStyleClass( "ColorThemeEditContainer" );
  
  m_edit = new ColorThemeWidget(editDiv);
  m_edit->edited().connect( boost::bind( &ColorThemeWindow::themEditedCallback, this ) );
  
  m_edit->addStyleClass( "ColorThemeWindowContents" );
  
  WGridLayout *leftMenuDivLayout = new WGridLayout();
  leftMenuDivLayout->setVerticalSpacing( 0 );
  leftMenuDivLayout->setHorizontalSpacing( 0 );
  leftMenuDivLayout->setContentsMargins( 0, 0, 0, 0 );
  
  leftMenuDiv->setLayout( leftMenuDivLayout );
  
  WContainerWidget *menuDiv = new WContainerWidget();
  menuDiv->setOverflow( WContainerWidget::OverflowAuto, Wt::Vertical );
  menuDiv->setOverflow( WContainerWidget::OverflowHidden, Wt::Horizontal );
  
  m_menu = new WMenu(Wt::Vertical,menuDiv);
  leftMenuDivLayout->addWidget( menuDiv, 0, 0 );
  //m_menu->addStyleClass((m_interspec->isMobile() ? "SideMenuPhone" : "SideMenu")); //defined in InterSpec.css
  m_menu->addStyleClass( "SideMenu VerticalNavMenu HeavyNavMenu" );
  
  WContainerWidget *adSubDiv = new WContainerWidget();
  leftMenuDivLayout->addWidget( adSubDiv, 1, 0 );
  adSubDiv->addStyleClass( "AddRemoveThemeButtons" );
  m_removeIcn = new WContainerWidget(adSubDiv); //needed or else button wont show up
  m_removeIcn->setStyleClass( "DeleteColorTheme Wt-icon" );
  m_removeIcn->setInline( true );
  m_removeIcn->clicked().connect( this, &ColorThemeWindow::removeThemeCallback );
  m_removeIcn->setToolTip( "Removes the currently selected theme.", Wt::PlainText );
  
  m_addIcn = new WContainerWidget(adSubDiv); //needed or else button wont show up
  m_addIcn->setStyleClass( "AddColorTheme Wt-icon" );
  m_addIcn->clicked().connect( this, &ColorThemeWindow::cloneThemeCallback );
  m_addIcn->setInline( true );
  m_addIcn->setToolTip( "Creates a new color theme based upon the currently selected theme.", Wt::PlainText );
  
  WPushButton *upload = new WPushButton( "", adSubDiv );
  upload->setIcon( WLink("InterSpec_resources/images/upload_small.svg") );
  upload->addStyleClass( "LinkBtn UploadBtn UploadColorTheme" );
  upload->setToolTip( "Imports a color theme JSON file to use.", Wt::PlainText );
  upload->clicked().connect( this, &ColorThemeWindow::uploadThemeCallback );
  
  JsonDownloadResource *downloadResource = new JsonDownloadResource( this );
  
#if( BUILD_AS_OSX_APP || IOS )
  WAnchor *download = new WAnchor( WLink(downloadResource), adSubDiv );
  download->setTarget( AnchorTarget::TargetNewWindow );
  download->setStyleClass( "LinkBtn DownloadLink DownloadColorTheme" );
  download->setTarget( AnchorTarget::TargetNewWindow );
  download->setText( "&nbsp;" );
#else
  WPushButton *download = new WPushButton( adSubDiv );
  download->setIcon( "InterSpec_resources/images/download_small.svg" );
  download->setLink( WLink(downloadResource) );
  download->setLinkTarget( Wt::TargetNewWindow );
  download->setStyleClass( "LinkBtn DownloadBtn DownloadColorTheme" );
  
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  download->clicked().connect( std::bind([downloadResource](){
    android_download_workaround(downloadResource, "color_theme.xml");
  }) );
#endif //ANDROID

#endif
  
  download->setToolTip( "Exports the currently selected theme to a JSON file.", Wt::PlainText );
  
  leftMenuDivLayout->setRowStretch( 0, 1 );
  
  WGridLayout *layout = stretcher();
  layout->addWidget(leftMenuDiv, 0, 0, AlignLeft);
  layout->addWidget(editDiv, 0, 1 );
  layout->setRowStretch(0, 1);
  layout->setColumnStretch( 1, 1 );
  
  WContainerWidget *foot = footer();
  AuxWindow::addHelpInFooter( foot, "color-theme-dialog" );
  
  const bool phone = isPhone();
  WCheckBox *autoDarkCb = new WCheckBox( "Auto apply \"Dark\" from OS" );
  autoDarkCb->setFloatSide( Wt::Side::Left );
  autoDarkCb->setToolTip( "Apply the \"Dark\" color theme automatically according to"
                          " the operating systems current value, or when it transisitions." );
  
  if( !phone )
    foot->addWidget( autoDarkCb );
  
  InterSpecUser::associateWidget(m_interspec->m_user, "AutoDarkFromOs", autoDarkCb, m_interspec);
  
  m_save = new WPushButton( "Save", foot );
  m_apply = new WPushButton( "Apply", foot );
  m_save->setHiddenKeepsGeometry( true );
  m_apply->setHiddenKeepsGeometry( true );
  
  m_save->clicked().connect( boost::bind( &ColorThemeWindow::saveCallback, this ) );
  m_apply->clicked().connect( boost::bind( &ColorThemeWindow::applyCallback, this ) );
  
  m_close = AuxWindow::addCloseButtonToFooter( "Close", true, foot );
  m_close->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  
  if( phone ) //Keep "Close" the left most item
    foot->addWidget( autoDarkCb );
  
  //To the menu add
  vector<unique_ptr<ColorTheme>> defaultThemes = ColorTheme::predefinedThemes();
  
  //Now grab from the
  vector<unique_ptr<ColorTheme>> dbThemes = userDbThemes();
  
  ThemeMenuItem *currentItem = nullptr;
  const shared_ptr<const ColorTheme> currentColorTheme = m_interspec->getColorTheme();
  //long long colorThemeIndex = -1;
  //if( currentColorTheme )
  //  colorThemeIndex = currentColorTheme->dbIndex;
  //else
  //  colorThemeIndex = InterSpecUser::preferenceValue<int>("ColorThemeIndex", m_interspec);
  
  for( unique_ptr<ColorTheme> &p : defaultThemes )
  {
    const WString name = p->theme_name;
    
    ThemeMenuItem *item = new ThemeMenuItem(name, std::move(p), false );
    m_menu->addItem(item);
    item->clicked().connect( boost::bind(&ColorThemeWindow::selectItem, this, item) );
    
    if( currentColorTheme && (currentColorTheme->theme_name == name) )
      currentItem = item;
  }//for( size_t i = 0; i < defaultThemes.size(); ++i )
  
  
  m_menu->addSeparator();
  
  
  for( unique_ptr<ColorTheme> &p : dbThemes )
  {
    //const auto dbindex = p->dbIndex;
    const WString name = p->theme_name;
    
    ThemeMenuItem *item = new ThemeMenuItem(name, std::move(p), true );
    m_menu->addItem(item);
    item->clicked().connect( boost::bind(&ColorThemeWindow::selectItem, this, item) );
    
    // Note: we cant just check the database for a matching index because app-states save the
    //   color theme JSON to a field in the database, not a link to the database entry.
    //   So we just have to match by name - not perfect, but close enough.
    if( !currentItem && currentColorTheme && (currentColorTheme->theme_name == name) )
      currentItem = item;
  }//for( unique_ptr<ColorTheme> &p : dbThemes )
  
  
  bool currentThemeSelected = true;
  if( !currentItem )
  {
    currentThemeSelected = false;
    currentItem = dynamic_cast<ThemeMenuItem *>( m_menu->itemAt(0) );
  }//if( !currentItem )
  
  if( currentItem )
  {
    m_menu->select( currentItem );
    m_edit->setTheme( currentItem->theme(), currentItem->editable() );
    if( currentThemeSelected )
      m_apply->hide();
  }//if( currentItem )
  
  m_menu->itemSelected().connect(boost::bind(&ColorThemeWindow::themeSelected, this,
                                             boost::placeholders::_1));
  
  finished().connect(boost::bind(&ColorThemeWindow::checkForSavesAndCleanUp, this));
  
  setClosable( true );
  show();
  centerWindow();
}//ColorThemeWindow


ColorThemeWindow::~ColorThemeWindow()
{
  
}//~ColorThemeWindow


void ColorThemeWindow::cloneThemeCallback()
{
  auto item = dynamic_cast<ThemeMenuItem *>(m_menu->currentItem());
  if( !item || !item->theme() ) //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::cloneThemeCallback(): Error getting current theme item" << endl;
    return;
  }
  
  unique_ptr<ColorTheme> newtheme( new ColorTheme( *item->theme() ) );
  newtheme->dbIndex = -1;
  
  newtheme->theme_name += " (cloned)";
  newtheme->creation_time = WDateTime::currentDateTime();
  newtheme->modified_time = newtheme->creation_time;
  auto themename = newtheme->theme_name;
  ThemeMenuItem *newitem = new ThemeMenuItem( themename, std::move(newtheme), true );
  m_menu->addItem( newitem );
  m_menu->select( newitem );
  newitem->clicked().connect( boost::bind(&ColorThemeWindow::selectItem, this, newitem) );
  
  showOrHideApplyButton();
}//void cloneThemeCallback()


void ColorThemeWindow::removeThemeCallback()
{
  auto item = dynamic_cast<ThemeMenuItem *>(m_menu->currentItem());
  if( !item || !item->theme() ) //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::removeThemeCallback(): Error getting current theme item" << endl;
    return;
  }
  
  if( !item->editable() )  //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::removeThemeCallback(): Current theme is not editable" << endl;
    return;
  }
  
  const WString name = item->text();
  const long long dbIndex = item->theme()->dbIndex;
  
  AuxWindow *conf = new AuxWindow( "Delete " + name + "?",
                                  WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse)
                                  | AuxWindowProperties::IsModal | AuxWindowProperties::TabletNotFullScreen );
  
  auto doDelete = [this,dbIndex,item,conf](){
    if( dbIndex >= 0 )
    {
      try
      {
        std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
        
        DataBaseUtils::DbTransaction transaction( *sql );
        Dbo::ptr<ColorThemeInfo> dbentry = m_interspec->m_user->colorThemes().find().where( "id = ?" ).bind( dbIndex );
        dbentry.remove();
        transaction.commit();
      }catch( Wt::Dbo::Exception & )
      {
        m_interspec->logMessage( "Error removing theme from database.", 2 ); //2 = WarningWidget::WarningMsgLevel::WarningMsgMedium
      }
    }//if( in DB )
    
    m_menu->removeItem( item );
    m_menu->select( 0 );
    
    delete item;
    conf->hide();
    
    showOrHideApplyButton();
  };//doDelete
  
  WText *text = new WText( "<span style=\"white-space: nowrap;\">Are you sure you would like to delete <em>" + name + "</em>?</span>" );
  conf->contents()->addWidget( text );
  conf->contents()->setMargin( 15 );

  WContainerWidget *bottom = conf->footer();
  
  WPushButton *cancel = new WPushButton( "No", bottom );
  cancel->clicked().connect( conf, &AuxWindow::hide );
  cancel->setWidth( WLength(47.5,WLength::Percentage) );
  cancel->setFloatSide( Wt::Left );
  
  WPushButton *erase = new WPushButton( "Yes", bottom );
  erase->clicked().connect( std::bind( doDelete ) );
  erase->setWidth( WLength(47.5,WLength::Percentage) );
  erase->setFloatSide( Wt::Right );
  
  conf->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, conf ) );
  
  conf->show();
  conf->rejectWhenEscapePressed();
  conf->resizeToFitOnScreen();
  conf->centerWindow();
  
  showOrHideApplyButton();
}//void removeThemeCallback()


void ColorThemeWindow::uploadThemeCallback()
{
  //ToDo 20181126: This function is currently only minimally implemented, should
  //  deal with errors better, and cleaup the code/logic to be a bit less
  //  spaghetti (e.g., get rid of nested lambdas, etc)

  //1) Make a dialog to allow user to select a file.
  AuxWindow *window = new AuxWindow( "Upload Color Theme JSON",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                     | AuxWindowProperties::TabletNotFullScreen | AuxWindowProperties::DisableCollapse
                                     | AuxWindowProperties::SetCloseable) );
  
  WContainerWidget *contents = window->contents();
  
  WPushButton *close = window->addCloseButtonToFooter( "Cancel" );
  close->clicked().connect( boost::bind( &AuxWindow::hide, window ) );
  
  WFileUpload *upload = new WFileUpload( contents );
  upload->setInline( false );
  
  //2) Upload file and check size.  If to big, reject and display error message.
  upload->uploaded().connect( std::bind( [this,window,upload](){
    const string filename = upload->spoolFileName();
    upload->hide();
    try
    {
      ifstream input( filename.c_str() );
      if( !input )
        throw runtime_error( "Error opening uploaded file." );
      input.unsetf(ios::skipws);
      
      // Determine file size
      input.seekg(0, ios::end);
      const size_t filesize = input.tellg();
      input.seekg(0);
      if( !filesize )
        throw runtime_error( "Uploaded file empty." );
      
      if( filesize > 64*1024 )
        throw runtime_error( "Uploaded file too large; max size 64kb." );
      
      //Read file data
      string data;
      data.resize( filesize );
      if( !input.read( &data.front(), static_cast<streamsize>(filesize) ) )
        throw runtime_error( "Unable to read in uploaded file." );
      
      //3) parse file.  If not valid display error message.
      auto theme = make_shared<ColorTheme>();
      try
      {
        ColorTheme::fromJson(data, *theme);
        
        //4) Ask user if they want to save to database (show disabled version of edit widget)
        new WText( "Would you like to import <em>" + theme->theme_name.toUTF8() + "</em>?", window->contents() );
        WContainerWidget *editParent = new WContainerWidget( window->contents() );
        
        const int w = m_interspec->renderedWidth();
        const int h = m_interspec->renderedHeight();
        
        editParent->setMaximumSize( 0.75*w, 0.75*h );
        ColorThemeWidget *edit = new ColorThemeWidget( editParent );
        edit->setTheme( theme.get(), false );
        WPushButton *yes = window->addCloseButtonToFooter( "Yes" );
        yes->clicked().connect( std::bind( [=](){
          
          //5) Add to database, and set as current theme to display
          std::unique_ptr<ColorTheme> newtheme = this->saveThemeToDb( theme.get() );
          if( newtheme )
          {
            ThemeMenuItem *newitem = new ThemeMenuItem( newtheme->theme_name, std::move( newtheme ), true );
            m_menu->addItem( newitem );
            m_menu->select( newitem );
            newitem->clicked().connect( boost::bind(&ColorThemeWindow::selectItem, this, newitem) );
            m_apply->show();
          }else
          {
             m_interspec->logMessage( "Error saving theme to database - sorry!.", 2 ); //2 = WarningWidget::WarningMsgLevel::WarningMsgInfo
          }
          window->hide();
        }) );
        
        window->resizeToFitOnScreen();
        window->centerWindow();
      }catch( std::exception & )
      {
        throw;
      }catch(...)
      {
        throw runtime_error( "Invalid JSON, or missing required info." );
      }
    }catch( std::exception &e )
    {
      new WText( e.what(), window->contents() );
    }
  }) );
  
  upload->fileTooLarge().connect( std::bind( [window,upload](){
    upload->hide();
    new WText( "Too Large of file.", window->contents() );
  }) );
  upload->changed().connect( upload, &WFileUpload::upload );
  
  
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->show();
  window->resizeToFitOnScreen();
  window->centerWindow();
}//void uploadThemeCallback()


void ColorThemeWindow::selectItem( WMenuItem *item )
{
  m_menu->select( item );
  item->triggered().emit( item ); //doenst look like this is emmitted either
                                  //when body of SideMenuItem is clicked
}//void select_item(  SideMenuItem *item )


void ColorThemeWindow::checkForSavesAndCleanUp()
{
  if( !m_edit->m_currentTheme )
  {
    AuxWindow::deleteAuxWindow( this );
    return;
  }
  
  string content = "<span style=\"white-space: nowrap;\">Save changes to <em>"
                + m_edit->m_currentTheme->theme_name.toUTF8() + "</em>?</span>";
  SimpleDialog *dialog = new SimpleDialog( "Save Changes?", content );
  
  
  WPushButton *discard = dialog->addButton( "Discard" );
  discard->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  WPushButton *save = dialog->addButton( "Save" );
  save->clicked().connect( boost::bind( &ColorThemeWindow::saveAndDelete, this ) );
}//void checkForSavesAndCleanUp()


void ColorThemeWindow::saveAndDelete()
{
  saveThemeToDb( m_edit->m_currentTheme.get() );
  AuxWindow::deleteAuxWindow( this );
}


unique_ptr<ColorTheme> ColorThemeWindow::saveThemeToDb( const ColorTheme *theme )
{
  if( !theme )
    return nullptr;
  
  unique_ptr<ColorTheme> newTheme;
  
  string errormsg;
  Dbo::ptr<InterSpecUser> &user = m_interspec->m_user;
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
  
  if( !theme || !user || !sql )
  {
    m_interspec->logMessage( "Programming Logic Error in ColorThemeWindow::saveCallback(): unable to get database resource", 2 ); //2 = WarningWidget::WarningMsgLevel::WarningMsgMedium
    return newTheme;
  }
  
  if( theme->dbIndex < 0 )
  {
    Wt::Dbo::ptr<ColorThemeInfo> dbtheme( new ColorThemeInfo() );
    
    try
    {
      dbtheme.modify()->user = user;
      dbtheme.modify()->modified_time = dbtheme.modify()->creation_time = WDateTime::currentDateTime();
      
      theme->setToDataBase( *dbtheme.modify() );
      
      DataBaseUtils::DbTransaction transaction( *sql );
      
      dbtheme = sql->session()->add( dbtheme );
      
      transaction.commit();
      
      newTheme.reset( new ColorTheme() );
      newTheme->setFromDataBase( *dbtheme );
      
      newTheme->dbIndex = dbtheme.id();
    } catch( Wt::Dbo::Exception &e )
    {
      errormsg = "Error writing new color theme to database: " + string( e.what() );
      newTheme.reset();
    } catch( std::exception &e )
    {
      errormsg = "Error serializing new color theme: " + string( e.what() );
      newTheme.reset();
    }
  } else
  {
    try
    {
      DataBaseUtils::DbTransaction transaction( *sql );
      
      const Wt::Dbo::collection< Wt::Dbo::ptr<ColorThemeInfo> > &userthemes = user->colorThemes();
      
      Wt::Dbo::ptr<ColorThemeInfo> dbtheme = userthemes.find().where( "id = ?" ).bind( theme->dbIndex );
      
      if( !dbtheme )
        throw std::runtime_error( "Could not find wanted color theme in database." );
      
      theme->setToDataBase( *dbtheme.modify() );
      dbtheme.modify()->modified_time = WDateTime::currentDateTime();
      
      transaction.commit();
      
      newTheme.reset( new ColorTheme() );
      newTheme->setFromDataBase( *dbtheme );
      
      newTheme->dbIndex = dbtheme.id();
    } catch( Wt::Dbo::Exception &e )
    {
      errormsg = "Error writing existing color theme to database: " + string( e.what() );
      newTheme.reset();
    } catch( std::exception &e )
    {
      errormsg = "Error serializing existing color theme: " + string( e.what() );
      newTheme.reset();
    }
  }//if( in not in database ) / else ( already in DB )
  
  if( errormsg.empty() )
    m_interspec->logMessage( theme->theme_name + " color theme saved.", 0 ); //2 = WarningWidget::WarningMsgLevel::WarningMsgInfo
  else
    m_interspec->logMessage( errormsg, 2 ); //2 = WarningWidget::WarningMsgLevel::WarningMsgMedium
  
  return newTheme;
}//void saveThemeToDb( const ColorTheme *theme )


std::string ColorThemeWindow::currentThemeTitle()
{
  const ColorTheme *theme = nullptr;
  if( m_edit->m_currentTheme )
    theme = m_edit->m_currentTheme.get();
  else
    theme = m_edit->m_origTheme.get();
  
  if( !theme )
    return "";
  
  return theme->theme_name.toUTF8();
}//std::string currentThemeTitle()


std::string ColorThemeWindow::currentThemeJson()
{
  const ColorTheme *theme = nullptr;
  if( m_edit->m_currentTheme )
    theme = m_edit->m_currentTheme.get();
  else
    theme = m_edit->m_origTheme.get();
  
  if( !theme )
    return "";
  
  return theme->toJson(*theme);
}//std::string currentThemeJson()


void ColorThemeWindow::saveCallback()
{
  auto item = dynamic_cast<ThemeMenuItem *>( m_menu->currentItem() );
  if( !item ) //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::saveCallback(): Error getting current theme item" << endl;
    return;
  }
  
  if( !item->editable() || !item->theme() )  //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::saveCallback(): Current theme is not editable" << endl;
    return;
  }
  
  if( !m_edit->m_currentTheme )
    return;
  
  auto newTheme = saveThemeToDb( m_edit->m_currentTheme.get() );
  
  if( newTheme )
    item->setTheme( std::move(newTheme) );
  
  m_edit->setTheme( item->theme(), true );
  
  m_save->hide();
}//void saveCallback()


void ColorThemeWindow::applyCallback()
{
  auto item = dynamic_cast<ThemeMenuItem *>(m_menu->currentItem());
  if( !item ) //shouldnt ever happen
  {
    cerr << "ColorThemeWindow::applyCallback(): Error getting current theme item" << endl;
    return;
  }
  
  saveCallback();
  
  const ColorTheme *theme = item->theme();
  assert( theme ); //def shouldnt ever happen
  
  Dbo::ptr<InterSpecUser> &user = m_interspec->m_user;
  InterSpecUser::setPreferenceValue<int>( user, "ColorThemeIndex", static_cast<int>(theme->dbIndex), m_interspec );
  
  m_interspec->applyColorTheme( make_shared<ColorTheme>(*theme) );
  
  item->setEditedSinceApply( false );
  m_apply->hide();
}//void applyCallback()


void ColorThemeWindow::themEditedCallback()
{
  m_save->show();
  m_apply->show();
  
  ThemeMenuItem *themeItem = dynamic_cast<ThemeMenuItem *>(m_menu->currentItem());
  if( themeItem )
    themeItem->setEditedSinceApply( true );  //should always be called
}//void themEditedCallback()


std::vector<std::unique_ptr<ColorTheme>>  ColorThemeWindow::userDbThemes()
{
  std::vector<std::unique_ptr<ColorTheme>> dbthemes;
  
  try
  {
    Dbo::ptr<InterSpecUser> &user = m_interspec->m_user;
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_interspec->sql();
    
    if (!user || !sql)
      throw runtime_error("Invalid database ish");
    
    DataBaseUtils::DbTransaction transaction(*sql);
    user.reread();
    
    const Wt::Dbo::collection< Wt::Dbo::ptr<ColorThemeInfo> > &userthemes = user->colorThemes();
    
    for (auto i = userthemes.begin(); i != userthemes.end(); ++i)
    {
      unique_ptr<ColorTheme> t( new ColorTheme() );
      try
      {
        t->setFromDataBase(*(*i));
        t->dbIndex = (*i).id();
        dbthemes.push_back( std::move(t) );
      }catch( std::exception &e )
      {
        cerr << "Failed to se-serialize theme '" << (*i)->theme_name.toUTF8() << ": " << e.what() << endl;
      }
    }
    
    //ToDo: check for duplicate names, and if so, maybe add the date or something
    
    transaction.commit();
  }catch( Wt::Dbo::Exception &e )
  {
    cerr << "ColorThemeWindow::userDbThemes(): Caught exception searching database for users themes: " << e.what() << " code=" << e.code() << endl;
  }catch( std::exception &e )
  {
    cerr << "ColorThemeWindow::userDbThemes(): Caught error searching database for users themes: " << e.what() << endl;
  }
  
  return dbthemes;
}//userDbThemes()



void ColorThemeWindow::themeSelected( Wt::WMenuItem *item )
{
  auto themeItem = dynamic_cast<ThemeMenuItem *>(item);
  
  if( !themeItem )
  {
    cerr << "ColorThemeWindow::themeSelected(" << item << "): item passed in not a ThemeMenuItem" << endl;
    return;
  }
  
  m_edit->setTheme( themeItem->theme(), themeItem->editable() );
  m_removeIcn->setHidden( !themeItem->editable() );
  m_save->setHidden( true );
  
  showOrHideApplyButton();
}//themeSelected(Wt::WMenuItem *item)


void ColorThemeWindow::showOrHideApplyButton()
{
  Wt::WMenuItem *item = m_menu->currentItem();
  auto themeItem = dynamic_cast<ThemeMenuItem *>(item);
  
  if( !themeItem )
  {
    cerr << "ColorThemeWindow::showOrHideApplyButton(" << item << "): is not a ThemeMenuItem" << endl;
    return;
  }
  
  std::shared_ptr<const ColorTheme> currentTheme = m_interspec->getColorTheme();
  
  bool showApply = false;
  if( !currentTheme || !themeItem->theme() )
    showApply = true;
  else if( currentTheme->dbIndex < 0 || themeItem->theme()->dbIndex < 0 )
    showApply = true;
  else if( currentTheme->dbIndex != themeItem->theme()->dbIndex )
    showApply = true;
  else if( themeItem->editable() && themeItem->isEditedSinceApply() )
    showApply = true;
  
  m_apply->setHidden( !showApply );
}//void showOrHideApplyButton()
