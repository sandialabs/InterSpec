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

#include "InterSpec/DirectorySelector.h"

#include <iostream>
#include <functional>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WServer.h>
#include <Wt/WString.h>
#include <Wt/WLineEdit.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"

#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#endif

#if( BUILD_AS_OSX_APP )
#include "target/osx/macOsUtils.h"
#endif

#if( BUILD_AS_WX_WIDGETS_APP )
namespace
{
  // The wxWidgets directory picker is defined in the wx executable target
  // (where wxApp::ms_appInstance is valid) and registered into LibInterSpec
  // at startup via set_wx_native_directory_picker().  This is a plain
  // function pointer (not std::function) so the value stored here is just an
  // address - no ABI concerns crossing the DLL boundary on Windows /MT
  // builds where the EXE and DLL have separate CRT copies.
  WxNativeDirectoryPickerFn s_wx_native_directory_picker = nullptr;
}//namespace

void set_wx_native_directory_picker( WxNativeDirectoryPickerFn callback )
{
  s_wx_native_directory_picker = callback;
}
#endif // BUILD_AS_WX_WIDGETS_APP


DirectorySelector::DirectorySelector()
  : WContainerWidget(),
    m_currentPath(),
    m_isValid( false ),
    m_label( nullptr ),
    m_pathDisplay( nullptr ),
    m_selectButton( nullptr ),
    m_pathInput( nullptr ),
    m_useNativeDialog( false ),
    m_useEditableField( false )
{
  wApp->useStyleSheet( "InterSpec_resources/DirectorySelector.css" );
  InterSpec::instance()->useMessageResourceBundle( "DirectorySelector" );

  addStyleClass( "DirectorySelector" );

  // Determine if we should use native dialog
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( InterSpecApp::isPrimaryWindowInstance() )
  {
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
    m_useNativeDialog = true;
#endif
  }
#endif

  // The path widget is editable everywhere except macOS-with-native-picker.
  // The macOS native build is sandboxed; only paths returned by NSOpenPanel
  // come with security-scoped access, so a typed/pasted path wouldn't work
  // anyway.  On Linux/Windows (Electron and wxWidgets) both paths work, so
  // expose an editable line edit alongside the picker button.
  m_useEditableField = true;
#if( defined(__APPLE__) )
  if( m_useNativeDialog )
    m_useEditableField = false;
#endif

  setupUI();
}


DirectorySelector::~DirectorySelector()
{
  // Nothing special needed for cleanup
}

void DirectorySelector::setLabelTxt( const Wt::WString &labelTxt )
{
  if( m_label )
    m_label->setText( labelTxt );
}//void setLabelTxt()


std::string DirectorySelector::path() const
{
  return m_currentPath;
}


void DirectorySelector::setPath( const std::string &path, bool validate )
{
  if( m_currentPath == path )
    return;

  m_currentPath = path;

  // Update whichever path widget actually exists.
  if( m_pathDisplay )
    m_pathDisplay->setText( Wt::WString::fromUTF8(path) );
  if( m_pathInput )
    m_pathInput->setText( Wt::WString::fromUTF8(path) );
  
  if( validate )
    validatePath();
  
  // Emit path changed signal
  m_pathChanged.emit( m_currentPath );
}


bool DirectorySelector::isPathValid() const
{
  return m_isValid;
}


void DirectorySelector::setEnabled( bool enabled )
{
  WContainerWidget::setDisabled( !enabled );
  
  if( m_selectButton )
    m_selectButton->setDisabled( !enabled );
    
  if( m_pathInput )
    m_pathInput->setDisabled( !enabled );
}


void DirectorySelector::setupUI()
{
  // Add label
  m_label = addNew<Wt::WLabel>( Wt::WString::tr("ds-path-label") );

  // Path widget: editable line edit, or read-only text on macOS-with-picker.
  if( m_useEditableField )
  {
    m_pathInput = addNew<Wt::WLineEdit>();
    m_pathInput->addStyleClass( "DirectorySelectorPathTxt" );
    m_pathInput->setPlaceholderText( Wt::WString::tr("ds-txt-input-placeholder") );

    m_pathInput->keyWentUp().connect( this, [this](){ handleTextInputChange(); } );
    m_pathInput->enterPressed().connect( this, [this](){ handleTextInputChange(); } );
    m_pathInput->blurred().connect( this, [this](){ handleTextInputChange(); } );
    m_pathInput->changed().connect( this, [this](){ handleTextInputChange(); } );
  }else
  {
    m_pathDisplay = addNew<Wt::WText>( Wt::WString::tr("ds-native-no-path-placeholder") );
    m_pathDisplay->addStyleClass( "DirectorySelectorPathTxt" );
  }

  // Browse button: shown whenever a native picker is available.
  if( m_useNativeDialog )
  {
    m_selectButton = addNew<Wt::WPushButton>( Wt::WString::tr("ds-native-select-path-button") );
    m_selectButton->clicked().connect( this, [this](){ handleNativeDirectorySelection(); } );
  }
}//void setupUI()


void DirectorySelector::validatePath()
{
  const bool wasValid = m_isValid;
  m_isValid = !m_currentPath.empty() && SpecUtils::is_directory( m_currentPath );
  
  // Update validity styling on whichever path widget is actually shown.
  Wt::WWidget * const pathWidget = m_pathInput
    ? static_cast<Wt::WWidget *>(m_pathInput)
    : static_cast<Wt::WWidget *>(m_pathDisplay);
  if( pathWidget )
  {
    if( m_isValid )
    {
      if( pathWidget->hasStyleClass( "DirectorySelectorInvalidPath" ) )
        pathWidget->removeStyleClass( "DirectorySelectorInvalidPath" );
      if( !pathWidget->hasStyleClass( "DirectorySelectorValidPath" ) )
        pathWidget->addStyleClass( "DirectorySelectorValidPath" );
    }else
    {
      if( pathWidget->hasStyleClass( "DirectorySelectorValidPath" ) )
        pathWidget->removeStyleClass( "DirectorySelectorValidPath" );
      if( !m_currentPath.empty() && !pathWidget->hasStyleClass( "DirectorySelectorInvalidPath" ) )
        pathWidget->addStyleClass( "DirectorySelectorInvalidPath" );
    }
  }
  
  // Emit validity changed signal if status changed
  if( wasValid != m_isValid )
    m_pathValidityChanged.emit( m_isValid );
}


void DirectorySelector::handleNativeDirectorySelection()
{
  if( !m_useNativeDialog )
    return;

#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
  auto on_select_callback = [this]( const std::vector<std::string> &paths ){
    assert( paths.empty() || (paths.size() == 1) );
    // If the user cancelled (paths empty), leave the previously selected
    // path unchanged.  setPath() routes the new value to whichever path
    // widget is shown (line edit or read-only text).
    if( !paths.empty() )
      setPath( paths[0], true );
  };

#if( BUILD_AS_ELECTRON_APP )
  auto electron_callback = [on_select_callback]( std::string value ){
    std::vector<std::string> paths;
    if( !value.empty() )
      paths.push_back( value );
    on_select_callback( paths );
  };

  ElectronUtils::browse_for_directory( "Select Directory",
                                       "Select directory.",
                                       electron_callback );
#elif( BUILD_AS_WX_WIDGETS_APP )
  if( !s_wx_native_directory_picker )
  {
    std::cerr << "DirectorySelector: wx native directory picker has not been"
                 " registered - call set_wx_native_directory_picker() at app startup."
              << std::endl;
    return;
  }

  // The wx executable's browse_for_directory invokes its result callback
  // directly on the wx main thread (synchronously, after the dialog closes).
  // We have to bounce that result back into the Wt session ourselves - and
  // crucially, the WServer::post call below has to happen from code compiled
  // into LibInterSpec, because Wt is statically linked separately into the
  // dylib and the executable, so each binary has its own WServer::instance()
  // pointer (the executable's is null).
  Wt::WApplication * const wapp = wApp;
  const std::string session_id = wapp ? wapp->sessionId() : std::string();

  auto wx_result_callback = [on_select_callback, session_id]( const std::vector<std::string> &paths ){
    Wt::WServer * const server = Wt::WServer::instance();
    if( !server || session_id.empty() )
      return;
    server->post( session_id, [on_select_callback, paths](){
      on_select_callback( paths );
      Wt::WApplication * const app = Wt::WApplication::instance();
      if( app )
        app->triggerUpdate();
    } );
  };

  s_wx_native_directory_picker( "Select Directory",
                                "Select directory.",
                                wx_result_callback );
#else
  static_assert( BUILD_AS_OSX_APP, "Need to update preprocessor ifs" );
  const bool canChooseFiles = false;
  const bool canChooseDirectories = true;
  const bool allowsMultipleSelection = false;

  macOsUtils::showFilePicker( "Select Directory",
                              "Select directory.",
                              canChooseFiles,
                              canChooseDirectories,
                              allowsMultipleSelection,
                              on_select_callback );
#endif

#endif // BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP
}


void DirectorySelector::handleTextInputChange()
{
  if( !m_pathInput )
    return;
    
  const std::string newPath = m_pathInput->text().toUTF8();
  
  if( newPath != m_currentPath )
  {
    m_currentPath = newPath;
    validatePath();
    m_pathChanged.emit( m_currentPath );
  }
} 