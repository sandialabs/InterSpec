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
#include <iostream>
#include <filesystem>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WTreeView>
#include <Wt/WLineEdit>
#include <Wt/WMessageBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RefSpectraModel.h"
#include "InterSpec/RefSpectraWidget.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"


using namespace Wt;
using namespace std;

namespace fs = std::filesystem;

namespace 
{
  bool has_spectra_files( const fs::path &path )
  {
    for( const auto &entry : fs::directory_iterator( path ) )
    {
      if( !SpecUtils::likely_not_spec_file( entry.path().string() ) )
        return true;
      if( fs::is_directory( entry.path() ) && has_spectra_files( entry.path() ) )
        return true;
    }
    return false;
  }
}// anonymous namespace


RefSpectraDialog::RefSpectraDialog( const Wt::WString &title )
  : SimpleDialog( title ),
    m_widget( nullptr )
{
  addStyleClass( "RefSpectraDialog" );

  m_widget = new RefSpectraWidget( contents() );
  m_widget->setHeight( 500 );
  
  addButton( WString::tr("Close") );
  
  // Set dialog size based on screen size
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  const bool portrait = ((interspec->renderedWidth() > 100) 
                         && (interspec->renderedHeight() > 100)
                         && (interspec->renderedWidth() < 480));
  if( interspec && (interspec->renderedWidth() > 100) && (interspec->renderedHeight() > 50) )
  {
    double dialogWidth = std::min(800.0, 0.5*interspec->renderedWidth());
    double dialogHeight = std::min(600.0, 0.8*interspec->renderedHeight());
    
    resize( dialogWidth, dialogHeight );
  }
  else
  {
    // Default size for when screen dimensions aren't available
    resize( 640, 480 );
  }
  
  rejectWhenEscapePressed();
}

RefSpectraDialog::~RefSpectraDialog()
{
  // The widget will be automatically deleted by Wt
} 

RefSpectraDialog *RefSpectraDialog::createDialog( const RefSpectraInitialBehaviour initialBehaviour, 
                                                  const SpecUtils::SpectrumType type )
{
  RefSpectraDialog *dialog = new RefSpectraDialog( WString::tr("rs-dialog-title") );
  dialog->m_widget->setLoadSpectrumType( type );

  switch( initialBehaviour )
  {
    case RefSpectraInitialBehaviour::Default:
      break;

    case RefSpectraInitialBehaviour::GuessMatchToForeground:
      dialog->m_widget->selectSimilarToForeground();
      break;
  }//switch( initialBehaviour )

  return dialog;
}

RefSpectraWidget::RefSpectraWidget( Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
  m_treeModel( nullptr ),
  m_treeView( nullptr ),
  m_stack( nullptr ),
  m_spectrum( nullptr ),
  m_dirInfoContainer( nullptr ),
  m_showCurrentForeground( nullptr ),
  m_refBackground( nullptr ),
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  m_addDirButton( nullptr ),
  m_deleteDirButton( nullptr ),
  m_showInExplorerButton( nullptr ),
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
  m_showAsComboBox( nullptr )
{
  setupUI();

  initBaseDirs();
}


RefSpectraWidget::~RefSpectraWidget()
{
}


void RefSpectraWidget::setupUI()
{
  wApp->useStyleSheet( "InterSpec_resources/RefSpectraWidget.css" );

  addStyleClass( "RefSpectraWidget" );
  
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  interspec->useMessageResourceBundle( "RefSpectraWidget" );

  // Update layout so widgets are vertically stacked on portrait phones.
  const bool portrait = ((interspec->renderedWidth() > 100) 
                         && (interspec->renderedHeight() > 100)
                         && (interspec->renderedWidth() < 480));

  WGridLayout *mainLayout = new WGridLayout();
  setLayout( mainLayout );
  
  // Tree view
  m_treeModel = new RefSpectraModel( this );
  m_treeView = new WTreeView();
  m_treeView->setModel( m_treeModel );
  m_treeView->selectionChanged().connect( this, &RefSpectraWidget::handleSelectionChanged );
  
  mainLayout->addWidget( m_treeView, 0, 0 );

  // Stack that will either show intro txt, the spectrum of selected file, or a text box with the readme.txt or readme.xml
  m_stack = new WStackedWidget();
  m_stack->addStyleClass( "RefSpectraStack" );

  const char *directions_key = interspec->isMobile() ? "rs-no-dir-selected-mobile" : "rs-no-dir-selected-desktop";
  WText *noDirSelected = new WText( WString::tr(directions_key) );
  m_stack->addWidget( noDirSelected );

  m_dirInfoContainer = new WContainerWidget();
  m_stack->addWidget( m_dirInfoContainer );

  m_spectrum = new D3SpectrumDisplayDiv();

  // Set up the spectrum display, making it compact, and match current color theme and log y setting
  const bool wasLogY = UserPreferences::preferenceValue<bool>( "LogY", interspec );
  m_spectrum->setCompactAxis( true );
  m_spectrum->applyColorTheme( interspec->getColorTheme() );
  m_spectrum->setXAxisTitle( "" );
  m_spectrum->setYAxisTitle( "", "" );
  m_spectrum->disableLegend();
  m_spectrum->setYAxisLog( wasLogY );
  //m_spectrum->showYAxisScalers( true );
  //m_spectrum->setShowPeakLabel( 0, true ); //SpectrumChart::PeakLabels::kShowPeakUserLabel
  //m_spectrum->setReferncePhotoPeakLines( ReferenceLineInfo() );
  interspec->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, m_spectrum, boost::placeholders::_1 ) );

  m_stack->addWidget( m_spectrum );
  
  if( portrait )
    mainLayout->addWidget( m_stack, 1, 0 );
  else
    mainLayout->addWidget( m_stack, 0, 1 );


  // Options
  WContainerWidget *optionsContainer = new WContainerWidget();
  optionsContainer->addStyleClass( "RefSpectraOptions" );
  if( portrait )
    mainLayout->addWidget( optionsContainer, 2, 0 );
  else
    mainLayout->addWidget( optionsContainer, 1, 0, 1, 2 );
  
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  // Add directory button
  m_addDirButton = new Wt::WContainerWidget( optionsContainer );
  m_addDirButton->clicked().connect( this, &RefSpectraWidget::startAddDirectory );
  m_addDirButton->addStyleClass( "AddDirBtn" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec );
  HelpSystem::attachToolTipOn( m_addDirButton, WString::tr("rs-tt-add-dir-btn"), showToolTips );

  m_deleteDirButton = new Wt::WContainerWidget( optionsContainer );
  m_deleteDirButton->clicked().connect( this, &RefSpectraWidget::removeCurrentlySelectedDirectory );
  m_deleteDirButton->addStyleClass( "DeleteDirBtn" );
  m_deleteDirButton->setHiddenKeepsGeometry( true );
  m_deleteDirButton->setHidden( true );
  HelpSystem::attachToolTipOn( m_deleteDirButton, WString::tr("rs-tt-delete-dir-btn"), showToolTips );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

  WContainerWidget *optionsSpacer = new WContainerWidget( optionsContainer );
  optionsSpacer->addStyleClass( "RefSpectraOptionsSpacer" );
  

  m_showCurrentForeground = new Wt::WCheckBox( WString::tr("rs-show-foreground"), optionsContainer );
  m_showCurrentForeground->addStyleClass( "CbNoLineBreak" );
  m_showCurrentForeground->changed().connect( this, &RefSpectraWidget::updatePreview );
  
  m_refBackground = new Wt::WCheckBox( WString::tr("rs-show-background"), optionsContainer );
  m_refBackground->addStyleClass( "CbNoLineBreak" );
  m_refBackground->changed().connect( this, &RefSpectraWidget::updatePreview );
  m_refBackground->setHiddenKeepsGeometry( true );
  m_refBackground->setHidden( true );
  
  if( portrait )
  {
    mainLayout->setRowStretch( 0, 1 );
    mainLayout->setRowStretch( 1, 1 );
  }else
  {
    mainLayout->setRowStretch( 0, 1 );
  }

  WContainerWidget *showAsContainer = new WContainerWidget();
  showAsContainer->addStyleClass( "ShowAsContainer" );

  if( portrait )
    mainLayout->addWidget( showAsContainer, 3, 0 );
  else 
    mainLayout->addWidget( showAsContainer, 2, 0, 1, 2 );


#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  #ifdef _WIN32
  const char *show_exlorer_key = "rs-show-in-explorer";
  #elif defined( __APPLE__ )
  const char *show_exlorer_key = "rs-show-in-finder";
  #else
  const char *show_exlorer_key = "rs-show-in-linux-file-manager";
  #endif
  m_showInExplorerButton = new Wt::WPushButton( WString::tr(show_exlorer_key), showAsContainer );

  m_showInExplorerButton->addStyleClass( "LinkBtn ShowInExplorerBtn" );
  m_showInExplorerButton->clicked().connect( this, &RefSpectraWidget::showInExplorer );
  m_showInExplorerButton->setHiddenKeepsGeometry( true );
  m_showInExplorerButton->setHidden( true );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

  WContainerWidget *showAsContainerSpacer = new WContainerWidget( showAsContainer );
  showAsContainerSpacer->addStyleClass( "RefSpectraOptionsSpacer" );

  WContainerWidget *showAsGroup = new WContainerWidget(showAsContainerSpacer);
  WLabel *showAsLabel = new WLabel( WString::tr("rs-show-as"), showAsGroup );
  showAsLabel->addStyleClass( "RefSpectraShowAsLabel" );
  m_showAsComboBox = new WComboBox( showAsGroup );
  m_showAsComboBox->addStyleClass( "RefSpectraShowAsComboBox" );
  m_showAsComboBox->addItem( WString::tr("Foreground") );
  m_showAsComboBox->addItem( WString::tr("Background") );
  m_showAsComboBox->addItem( WString::tr("Secondary") );
  m_showAsComboBox->setCurrentIndex( 2 );

  m_treeView->setSelectionMode( Wt::SelectionMode::SingleSelection );
}//void setupUI();

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
void RefSpectraWidget::startAddDirectory()
{
  SimpleDialog *dialog = new SimpleDialog( WString::tr("rs-add-dir-dialog-title") );
  WPushButton *okBtn = dialog->addButton( WString::tr("Save") );
  WPushButton *cancelBtn = dialog->addButton( WString::tr("Cancel") );

  dialog->addStyleClass( "RefSpectraAddDirDialog" );

  WContainerWidget *contents = dialog->contents();
  contents->addStyleClass( "RefSpectraAddDirDialogContents" );

  new WLabel( WString::tr("rs-add-dir-dialog-path"), contents );
  
  
#if( BUILD_AS_ELECTRON_APP )
  WContainerWidget *pathDiv = new WContainerWidget();
  linelayout->addWidget( pathDiv, 0, 1 );

  //ToDo: make a proper path selecting widget for Electron
  WPushButton *pickPath = new WPushButton( "Select Path", contents );
  WText *baseLocation = new WText( "(No Path Selected)", contents );
  baseLocation->addStyleClass( "SpecFileQueryPathTxt" );
  
  
  pickPath->clicked().connect( std::bind( [baseLocation](){
    // We need to use the lifetime management of WObject to safely bind the callback:
    auto bound_callback = boost::bind( &WText::setText, baseLocation, boost::placeholders::_1 );
    
    // Now we need to convert the boost bound thing, to a std::function, so we'll use a lamda as an
    //  intermediary to also capture the bound_function
    auto callback = [bound_callback]( string value ){
      bound_callback( value );
    };
    
    ElectronUtils::browse_for_directory( "Select Directory", "Select base-directory of reference spectra.", callback );
  }) );

  /*
  const string uploadname = dialog->id() + "PathPicker";
  const string uploadhtml = "<input id=\"" + uploadname + "\" type=\"file\" webkitdirectory=\"\" />";
  
  WText *baseLocation = new WText( uploadhtml, XHTMLUnsafeText, contents );
  
  //TODO: put in error handling!
  wApp->doJavaScript( "document.getElementById('" + uploadname + "').onchange = function(event){"
                         "var outputDir = document.getElementById('" + uploadname + "').files[0].path;"
                         "Wt.emit( \"" + id() + "\", { name: 'BaseDirSelected' }, outputDir );"
                     "};"
                     );
  */
#elif( BUILD_AS_OSX_APP )
static_assert( 0, "Need to implement and check" );
  WFileUpload *baseLocation = new WFileUpload( contents );
  baseLocation->changed().connect( std::bind([](){
    
  }) );
  linelayout->addWidget( m_baseLocation, 0, 1 );
#else
  WLineEdit *baseLocation = new WLineEdit( contents );
  auto callback = [baseLocation, okBtn](){
    std::string path = baseLocation->text().toUTF8();
    const bool valid_dir = SpecUtils::is_directory( path );
    okBtn->setEnabled( valid_dir );
  };
  baseLocation->keyPressed().connect( std::bind(callback) );
  baseLocation->setEmptyText( "Path to base directory to search" );

  okBtn->clicked().connect( std::bind([baseLocation, dialog](){
    std::string path = baseLocation->text().toUTF8();
    dialog->accept();
  }) );
#endif
}//


void RefSpectraWidget::addDirectory( const std::string &path )
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;
    
  if( SpecUtils::is_directory( path ) )
  {
    m_treeModel->addBaseDirectory( path );

    std::string normed_path;
    try
    {
      normed_path = SpecUtils::lexically_normalize_path( path );
    }catch( const std::exception &e )
    {
      assert( 0 );
      std::cerr << "RefSpectraWidget::addDirectory() - Invalid directory: '" << path << "'" << std::endl;
      normed_path = path;
    }

    const std::string ref_spectra_dirs = UserPreferences::preferenceValue<std::string>( "RefSpectraDirs", interspec );
    
    std::vector<std::string> base_dirs;
    SpecUtils::split( base_dirs, ref_spectra_dirs, ";" );

    bool is_new = true;
    for( const auto &dir : base_dirs )
    {
      std::string normed_dir;
      try
      {
        normed_dir = SpecUtils::lexically_normalize_path( dir );
      }catch( const std::exception &e )
      {
        std::cerr << "RefSpectraWidget::addDirectory() - Invalid directory: " << dir << std::endl;
        normed_dir = dir;
      }

      if( normed_dir == normed_path )
        is_new = false;
    }
    
    if( is_new )
    {
      base_dirs.push_back( path );
      std::string new_ref_spectra_dirs; 
      for( const auto &dir : base_dirs )
      {
        if( !new_ref_spectra_dirs.empty() )
          new_ref_spectra_dirs += ';';
        new_ref_spectra_dirs += dir;
      }

      try
      {
        UserPreferences::setPreferenceValue( "RefSpectraDirs", new_ref_spectra_dirs, interspec );
      }catch( const std::exception &e )
      {
        std::cerr << "RefSpectraWidget::addDirectory() - Error setting preferences: " << e.what() << std::endl;
        assert( 0 );
      }
    }else
    {
      passMessage( WString::tr("rs-duplicate-dir-path").arg( path ), WarningWidget::WarningMsgMedium );
      std::cerr << "RefSpectraWidget::addDirectory() - Invalid directory: " << path << std::endl;
    }
  }else
  {
    passMessage( WString::tr("rs-invalid-dir-path").arg( path ), WarningWidget::WarningMsgHigh );
    std::cerr << "RefSpectraWidget::addDirectory() - Invalid directory: " << path << std::endl;
  }
}

void RefSpectraWidget::removeCurrentlySelectedDirectory()
{
  const std::set<WModelIndex> selectedIndexes = m_treeView->selectedIndexes();
  assert( selectedIndexes.size() == 1 );
  if( selectedIndexes.empty() )
    return;

  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
    return;

  const std::string path = m_treeModel->getFilePath( index );
  assert( !path.empty() );
  assert( SpecUtils::is_directory( path ) );
  
  if( path.empty() || !SpecUtils::is_directory( path ) )
    return;

  m_treeModel->removeBaseDirectory( path );

  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  const std::string normed_path = ([&path](){
    try
    {
      return SpecUtils::lexically_normalize_path( path );
    }catch( const std::exception &e )
    {
      assert( 0 );
      std::cerr << "RefSpectraWidget::removeCurrentlySelectedDirectory() - Invalid directory: '" << path << "'" << std::endl;
    }
    return path;
  })();

  bool found_path_to_exclude = false;
  std::vector<std::string> new_base_dirs;
  std::vector<std::string> old_base_dirs;
  const std::string ref_spectra_dirs = UserPreferences::preferenceValue<std::string>( "RefSpectraDirs", interspec );
  SpecUtils::split( old_base_dirs, ref_spectra_dirs, ";" ); 
  for( const auto &dir : old_base_dirs )
  {
    std::string normed_dir;
    try
    {
      normed_dir = SpecUtils::lexically_normalize_path( dir );
    }catch( const std::exception &e )
    {
      assert( 0 );
      std::cerr << "RefSpectraWidget::removeCurrentlySelectedDirectory() - Invalid directory: '" << path << "'" << std::endl;
      normed_dir = dir;
    } 

    if( normed_dir != normed_path )
      new_base_dirs.push_back( dir );
    else
      found_path_to_exclude = true;
  }//for( const auto &dir : old_base_dirs )

  if( !found_path_to_exclude )
  {
    passMessage( WString::tr("rs-dir-not-found").arg( path ), WarningWidget::WarningMsgHigh );
    std::cerr << "RefSpectraWidget::removeCurrentlySelectedDirectory() - Directory not found: " << path << std::endl;
    assert( 0 );
    return;
  }//if( !found_path_to_exclude )

  std::string new_ref_spectra_dirs; 
  for( const auto &dir : new_base_dirs )
  {
    if( !new_ref_spectra_dirs.empty() )
      new_ref_spectra_dirs += ';';
    new_ref_spectra_dirs += dir;
  }

  try
  {
    UserPreferences::setPreferenceValue( "RefSpectraDirs", new_ref_spectra_dirs, interspec );
  }catch( const std::exception &e )
  {
    std::cerr << "RefSpectraWidget::removeCurrentlySelectedDirectory() - Error setting preferences: " << e.what() << std::endl;
    assert( 0 );
  }
}//void RefSpectraWidget::removeCurrentlySelectedDirectory()
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT


void RefSpectraWidget::updateTreeView()
{
  m_treeModel->refresh();
}


void RefSpectraWidget::handleSelectionChanged()
{
  const WModelIndexSet selectedIndexes = m_treeView->selectedIndexes();
  assert( selectedIndexes.empty() || (selectedIndexes.size() == 1) );
  
  if( selectedIndexes.empty() )
  {
    m_deleteDirButton->setHidden( true );
    m_showInExplorerButton->setHidden( true );
    m_refBackground->setHidden( true );
    //m_showCurrentForeground->setHidden( true );
    //m_showAsComboBox->setHidden( true );
    m_stack->setCurrentIndex( 0 );
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_showInExplorerButton->hide();
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

    return;
  }//if( selectedIndexes.empty() )


  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
    return;

  const bool isDir = m_treeModel->isDirectory( index );
  const std::string filePath = m_treeModel->getFilePath( index );

  m_dirInfoContainer->clear();

  if( isDir )
  {
    m_treeView->expand( index );
    
    m_deleteDirButton->setHidden( index.parent().parent().isValid() );
    m_showInExplorerButton->setHidden( false );
    m_refBackground->setHidden( true );
    //m_showCurrentForeground->setHidden( true );
    //m_showAsComboBox->setHidden( true );
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_showInExplorerButton->show();
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

    m_stack->setCurrentIndex( 1 );
    
    Wt::WText *title_text = new Wt::WText( m_dirInfoContainer );
    title_text->addStyleClass( "DirInfoTitle" );

    size_t num_spectra_files = 0, num_dirs = 0;
    const std::vector<std::string> files = SpecUtils::ls_files_in_directory( filePath );
    for( const auto &file : files )
    {
      if( SpecUtils::is_directory( file ) )
      {
        num_dirs += 1;
        continue;
      }

      const string fname = SpecUtils::filename( file );
      if( SpecUtils::iequals_ascii( fname, "readme.txt" ) || SpecUtils::iequals_ascii( fname, "readme.xml" ) )
      {
        string file_data;
        if( SpecUtils::file_size( file ) < 10*1024 )
        {
          try
          {
            std::vector<char> data;
            SpecUtils::load_file_data( file.c_str(), data );
            file_data = std::string( data.data(), data.size() );
          }catch( const std::exception &e )
          {
            std::cerr << "RefSpectraWidget::handleSelectionChanged() - Error loading file: " << file << std::endl;
          }
        }//if( SpecUtils::file_size( file ) < 10*1024 )

        Wt::WText *text = new Wt::WText( file_data, m_dirInfoContainer );
        text->addStyleClass( "DirReadmeText" );
      }else
      {
        num_spectra_files += SpecUtils::likely_not_spec_file( file );
      }
    }

    WString title_txt = WString::tr("rs-dir-info")
                          .arg( static_cast<int>( num_spectra_files ) )
                          .arg( static_cast<int>( num_dirs ) );
    title_text->setText( title_txt );

    m_spectrum->setData( nullptr, false );
    m_spectrum->setSecondData( nullptr );
    m_spectrum->setBackground( nullptr );
      
    return;
  }//if( isDir )

  m_refBackground->setHidden( true );
  m_deleteDirButton->setHidden( true );
  m_showInExplorerButton->setHidden( false );
  
  m_showCurrentForeground->setHidden( false );
  m_showAsComboBox->setHidden( false ); //JIC

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_showInExplorerButton->show();
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

  m_stack->setCurrentIndex( 2 );

  updatePreview();
}//void RefSpectraWidget::handleSelectionChanged()


void RefSpectraWidget::updatePreview()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  const WModelIndexSet selectedIndexes = m_treeView->selectedIndexes();
  assert( selectedIndexes.empty() || (selectedIndexes.size() == 1) );
  if( selectedIndexes.empty() )
    return;

  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
    return;

  const bool isDir = m_treeModel->isDirectory( index );
  assert( !isDir );
  if( isDir )
    return;

  const std::string filePath = m_treeModel->getFilePath( index );

  const string parent_path = SpecUtils::parent_path( filePath );
  const string filename = SpecUtils::filename( filePath );
  const string extension = SpecUtils::file_extension( filename );
  const vector<string> background_files = SpecUtils::ls_files_in_directory( parent_path,
    []( const std::string &filename, void *userdata ) -> bool{
      return SpecUtils::starts_with( filename, "background." );
  }, nullptr );

  SpecMeas infile;
  const bool loaded = infile.load_file( filePath, SpecUtils::ParserType::Auto, extension );
  if( !loaded )
  {
    m_stack->setCurrentIndex( 1 );
    WString error_text = WString::tr("rs-error-loading-file").arg( filePath );
    Wt::WText *text = new Wt::WText( error_text, m_dirInfoContainer );
    text->addStyleClass( "ErrorLoadingText" );

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_showInExplorerButton->hide();
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

    return;
  }//if( !loaded )
  
  // Now need to figure out foreground and background spectra.
  shared_ptr<const SpecUtils::Measurement> foreground_meas;
  shared_ptr<const SpecUtils::Measurement> background_meas;

  // TODO: properly check for foreground/background records...
  if( infile.num_measurements() >= 1 )
    foreground_meas = infile.measurement( size_t(0) );
  
  if( !background_meas && !background_files.empty() )
  {
    const string background_file = background_files.front();
    SpecMeas background_infile;
    const bool background_loaded = background_infile.load_file( background_file, SpecUtils::ParserType::Auto, extension );

    // TODO: properly check for background records...
    if( background_loaded )
      background_meas = background_infile.measurement( size_t(0) );
  }//if( !background_meas && !background_files.empty() )
  
  m_refBackground->setHidden( background_files.empty() && !background_meas );
  
  if( !m_refBackground->isChecked() )
    background_meas = nullptr;

  m_spectrum->setData( foreground_meas, false );
  m_spectrum->setBackground( background_meas );
  
  shared_ptr<const SpecUtils::Measurement> disp_foreground;
  if( m_showCurrentForeground->isChecked() )
    disp_foreground = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  m_spectrum->setSecondData( disp_foreground );

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_showInExplorerButton->show();
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT
}//void updatePreview()


void RefSpectraWidget::setLoadSpectrumType( SpecUtils::SpectrumType type )
{
  m_showAsComboBox->setCurrentIndex( static_cast<int>( type ) );
}


void RefSpectraWidget::selectSimilarToForeground()
{
  // TODO: Implement ...
}//void selectSimilarToForeground();


void RefSpectraWidget::initBaseDirs()
{
const std::string static_data_dir = InterSpec::staticDataDirectory();
  const std::string ref_spectra_dir = SpecUtils::append_path( static_data_dir, "reference_spectra" );
  if( SpecUtils::is_directory( ref_spectra_dir ) )
    m_treeModel->addBaseDirectory( ref_spectra_dir );

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  // Check for directories under reference_spectra in the "user_data" directory
  const std::string user_dir = InterSpec::writableDataDirectory();
  const std::string user_ref_spectra_dir = SpecUtils::append_path( user_dir, "reference_spectra" );
  if( SpecUtils::is_directory( user_ref_spectra_dir ) ) 
  {
    bool has_spectra_files_in_base = false;
    vector<string> dirs_to_add;
    for( const auto &entry : fs::directory_iterator( user_ref_spectra_dir ) )
    {
      if( fs::is_directory( entry.path() ) )
      {
        // Only add the directory if it has spectra files, or children directories that have spectra files
        if( has_spectra_files( entry.path() ) )
          dirs_to_add.push_back( entry.path().string() );
      }else{
        if( !SpecUtils::likely_not_spec_file( entry.path().string() ) )
          has_spectra_files_in_base = true;
      }
    }

    // If base dir has spectra files in it, we'll just add base directory.
    // Otherwise, we'll add the sub-directories that do have spectra files in them.
    if( has_spectra_files_in_base )
    {
      m_treeModel->addBaseDirectory( ref_spectra_dir );
    }else
    {
      for( const auto &dir : dirs_to_add )
        m_treeModel->addBaseDirectory( dir );
    }
  }
#endif


#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  
  const std::string ref_spectra_dirs = UserPreferences::preferenceValue<std::string>( "RefSpectraDirs", interspec );
  if( !ref_spectra_dirs.empty() )
  {
    std::vector<std::string> dirs;
    SpecUtils::split( dirs, ref_spectra_dirs, ";" );
    for( const auto &dir : dirs )
    {
      if( SpecUtils::is_directory( dir ) )
        m_treeModel->addBaseDirectory( dir );
      else
        std::cerr << "RefSpectraWidget::initBaseDirs() - Invalid directory: " << dir << std::endl;
    }
  }
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

}//void initBaseDirs();

void RefSpectraWidget::showInExplorer()
{
  // TODO: Implement show in explorer logic
}