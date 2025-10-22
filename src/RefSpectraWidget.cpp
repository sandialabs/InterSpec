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
#include <algorithm>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WServer>
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

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RefSpectraModel.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RefSpectraWidget.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"


#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
#include "InterSpec/DirectorySelector.h"
#endif

#if BUILD_AS_OSX_APP
#include "target/osx/macOsUtils.h"
#endif

using namespace Wt;
using namespace std;

namespace fs = std::filesystem;

namespace 
{
  const SpecUtils::SpectrumType ns_spectrum_types[] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };
  const size_t ns_spectrum_types_size = sizeof( ns_spectrum_types ) / sizeof( ns_spectrum_types[0] );

  bool has_spectra_files( const fs::path &path )
  {
    try
    {
      for( const auto &entry : fs::directory_iterator( path ) )
      {
        if( !SpecUtils::likely_not_spec_file( entry.path().string() ) )
          return true;
        if( fs::is_directory( entry.path() ) && has_spectra_files( entry.path() ) )
          return true;
      }
    }catch( std::exception &e )
    {
      cerr << "has_spectra_files: caught exception iterating directories: " << e.what() << endl;
    }
    
    return false;
  }//bool has_spectra_files( const fs::path &path )
}// anonymous namespace


RefSpectraDialog::RefSpectraDialog( const Wt::WString &title )
  : SimpleDialog( title ),
    m_widget( nullptr ),
    m_loadBtn( nullptr )
{
  addStyleClass( "RefSpectraDialog" );

  WGridLayout *layout = new WGridLayout();
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );

  contents()->setLayout( layout );

  m_widget = new RefSpectraWidget();
  layout->addWidget( m_widget, 0, 0 );

  m_loadBtn = addButton( WString::tr("Load") );
  m_loadBtn->clicked().connect( m_widget, &RefSpectraWidget::loadSelectedSpectrum );
  m_loadBtn->setDisabled( true );
  m_widget->fileSelectionChangedSignal().connect( boost::bind( &RefSpectraDialog::handleSelectionChanged, this, boost::placeholders::_1 ) );

  addButton( WString::tr("Cancel") );
  
  // Set dialog size based on screen size
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  if( interspec && interspec->isPhone() )
  {
    addStyleClass( "RefSpectraDialog-phone" );
#if( IOS )
    addStyleClass( "RefSpectraDialog-iphone" );
#endif
  }

  // We want to target 750px x 500px for normal - portrait phones, we'll take what we can get.
  //  (note: typical, old 7" tablets are at least 600x1024)
  const bool portrait = ((interspec->renderedWidth() > 100) 
                         && (interspec->renderedHeight() > 100)
                         && (interspec->renderedWidth() < 480));
  if( portrait )
  {
    addStyleClass( "RefSpectraDialog-portrait" );
    m_widget->setWidth( 0.95*interspec->renderedWidth() - 30 );
    m_widget->setHeight( 0.95*interspec->renderedHeight() - 90 );
    m_widget->setMinimumSize( 0.95*interspec->renderedWidth() - 30, 0.95*interspec->renderedHeight() - 90 );
    m_widget->setMaximumSize( 0.95*interspec->renderedWidth() - 30, 0.95*interspec->renderedHeight() - 90 );
  }else if( (interspec->renderedWidth() > 100) && (interspec->renderedHeight() > 50) )
  {
    double dialogWidth = std::min(750.0, 0.95*interspec->renderedWidth());
    double dialogHeight = std::min(500.0, 0.95*interspec->renderedHeight());
    
    m_widget->setWidth( dialogWidth - 30 );
    m_widget->setHeight( dialogHeight - 90 );
    m_widget->setMinimumSize( dialogWidth - 30, dialogHeight - 90 );
    m_widget->setMaximumSize( dialogWidth - 30, dialogHeight - 90 );
  }else
  {
    // Default size for when screen dimensions aren't available
    m_widget->setWidth( 750 - 30 );
    m_widget->setHeight( 500 - 90 );
    m_widget->setMinimumSize( 750 - 30, 500 - 90 );
    m_widget->setMaximumSize( 750 - 30, 500 - 90 );
  }
  
  rejectWhenEscapePressed();
}

RefSpectraDialog::~RefSpectraDialog()
{
  // The widget will be automatically deleted by Wt
} 

void RefSpectraDialog::handleSelectionChanged( RefSpectraWidgetSelectionType type )
{
  m_loadBtn->setDisabled( type != RefSpectraWidgetSelectionType::File );
}//void handleSelectionChanged( RefSpectraWidgetSelectionType type );


RefSpectraDialog *RefSpectraDialog::createDialog( const RefSpectraInitialBehaviour initialBehaviour, 
                                                  const SpecUtils::SpectrumType type )
{
  //WString title = WString::tr("rs-dialog-title");
  WString title;
  RefSpectraDialog *dialog = new RefSpectraDialog( title );
  dialog->m_widget->setLoadSpectrumType( type );

  switch( initialBehaviour )
  {
    case RefSpectraInitialBehaviour::Default:
      break;

    case RefSpectraInitialBehaviour::LastUserSelectedSpectra:
      dialog->m_widget->selectLastSelectedPath();
      break;

    //case RefSpectraInitialBehaviour::GuessMatchToForeground:
    //  dialog->m_widget->selectSimilarToForeground();
    //  break;
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
  m_showAsComboBox( nullptr ),
  m_baseDirectories{},
  m_fileSelectionChangedSignal( this ),
  m_lastCollapsedIndex{},
  m_lastExpandedIndex{},
  m_lastExpandTime{}
{
  setupUI();
  initBaseDirs();
}

void RefSpectraWidget::storeLastSelectedPath() noexcept
{
  try
  {
    const std::set<WModelIndex> selectedIndexes = m_treeView->selectedIndexes();
    if( selectedIndexes.size() != 1 )
      return;

    WModelIndex index = *selectedIndexes.begin();
    
    vector<string> selected_paths;
    while( index.isValid() )
    {
      selected_paths.push_back( m_treeModel->getDisplayName( index ) );
      index = m_treeModel->parent( index );
    }
    
    std::reverse( begin(selected_paths), end(selected_paths) );
    string selected_path;
    for( const auto &path : selected_paths )
    {
      if( !selected_path.empty() )
        selected_path += ";";
      selected_path += path;
    }

    InterSpec *interspec = InterSpec::instance();
    if( interspec && !selected_path.empty() )
      UserPreferences::setPreferenceValue( "RefSpectraLastSelection", selected_path, interspec );
  }catch( const std::exception &e )
  {
    std::cerr << "Error setting last selected path: " << e.what() << std::endl;
  }
}//void storeLastSelectedPath()
  
  
void RefSpectraWidget::selectLastSelectedPath()
{
  try
  {
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    if( !interspec )
      return;

    const std::string last_selection = UserPreferences::preferenceValue<std::string>( "RefSpectraLastSelection", interspec );
    if( last_selection.empty() )
      return;

    vector<string> selected_paths;
    SpecUtils::split( selected_paths, last_selection, ";" );
    
    WModelIndex index;
    bool changed_selection = false, current_selection_is_dir = false;
    
    for( const string &display_name : selected_paths )
    {
      index = m_treeModel->indexForDisplayName( display_name, index );
      if( !index.isValid() )
        break;
      m_treeView->expand( index );
      m_treeView->setSelectedIndexes( {index} );
      changed_selection = true;
      current_selection_is_dir = m_treeModel->isDirectory( index );
    }

    if( changed_selection )
    {
      if( current_selection_is_dir )
        m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::Directory );
      else
        m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::File );
    }//if( changed_selection )
  }catch(const std::exception& e)
  {
    std::cerr << "Error selecting last selected path: " << e.what() << std::endl;
  }
}//void selectLastSelectedPath()


RefSpectraWidget::~RefSpectraWidget()
{
  storeLastSelectedPath();
  
#if BUILD_AS_OSX_APP
  // Stop accessing all security scoped resources
  std::vector<std::string> bookmarkKeys = macOsUtils::getSecurityScopedBookmarkKeys("SecBookmark_RefSpectra_");
  for( const auto &key : bookmarkKeys )
  {
    macOsUtils::stopAccessingSecurityScopedBookmark(key);
  }
#endif
}//~RefSpectraWidget()


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
  mainLayout->setContentsMargins( 0, 0, 0, 0 );
  mainLayout->setVerticalSpacing( 0 );
  mainLayout->setHorizontalSpacing( 0 );
  mainLayout->setRowStretch( 0, 1 );
  mainLayout->setColumnStretch( 1, 1 );

  if( portrait )
    addStyleClass( "RefSpectraWidget-portrait" );
  else
    addStyleClass( "RefSpectraWidget-landscape" );

  setLayout( mainLayout );
  
  // Tree view
  m_treeModel = new RefSpectraModel( this );
  m_treeView = new Wt::WTreeView();
  m_treeView->addStyleClass( "RefSpectraTreeView" );
  //m_treeView->setMinimumSize( WLength(250), WLength::Auto );
  m_treeView->setModel( m_treeModel );
  m_treeView->selectionChanged().connect( this, &RefSpectraWidget::handleSelectionChanged );
  m_treeView->collapsed().connect( boost::bind( &RefSpectraWidget::handleCollapsed, this, boost::placeholders::_1 ) );
  m_treeView->expanded().connect( boost::bind( &RefSpectraWidget::handleExpanded, this, boost::placeholders::_1 ) );
  m_treeView->setSelectable( true );
  m_treeView->setSelectionMode( Wt::SelectionMode::SingleSelection );
  m_treeView->setColumnWidth( 0, WLength(250, WLength::Pixel) );
  m_treeView->setWidth( WLength(275, WLength::Pixel) ); 
  m_treeView->setColumnResizeEnabled( false );
  m_treeView->setSortingEnabled( false );

  m_treeModel->layoutAboutToBeChanged().connect( this, &RefSpectraWidget::handleLayoutAboutToBeChanged );
  m_treeModel->layoutChanged().connect( this, &RefSpectraWidget::handleLayoutChanged );

  WContainerWidget *treeViewWrapper = new WContainerWidget();
  treeViewWrapper->addStyleClass( "RefSpectraTreeViewWrapper" );
  WGridLayout *treeViewLayout = new WGridLayout();
  treeViewLayout->setContentsMargins( 0, 0, 0, 0 );
  treeViewLayout->setVerticalSpacing( 0 );
  treeViewLayout->setHorizontalSpacing( 0 );
  treeViewWrapper->setLayout( treeViewLayout );
  treeViewLayout->addWidget( m_treeView, 0, 0 );
  WContainerWidget *treeViewOptions = new WContainerWidget();
  treeViewOptions->addStyleClass( "RefSpectraOptions" );
  treeViewLayout->addWidget( treeViewOptions, 1, 0 );
  treeViewLayout->setRowStretch( 0, 1 );


  mainLayout->addWidget( treeViewWrapper, 0, 0 );

  // Stack that will either show intro txt, the spectrum of selected file, or a text box with the readme.txt or readme.xml
  m_stack = new WStackedWidget();
  m_stack->addStyleClass( "RefSpectraStack" );

  const char *directions_key = interspec->isMobile() ? "rs-no-dir-selected-mobile" : "rs-no-dir-selected-desktop";
  WText *noDirSelected = new WText( WString::tr(directions_key) );
  noDirSelected->addStyleClass( "RefSpectraNoDirSelected" );
  m_stack->addWidget( noDirSelected );

  m_dirInfoContainer = new WContainerWidget();
  m_stack->addWidget( m_dirInfoContainer );

  m_spectrum = new D3SpectrumDisplayDiv();

  // Set up the spectrum display, making it compact, and match current color theme and log y setting
  const bool wasLogY = UserPreferences::preferenceValue<bool>( "LogY", interspec );
  m_spectrum->setCompactAxis( true );
  m_spectrum->setXAxisTitle( "" );
  m_spectrum->setYAxisTitle( "", "" );
  m_spectrum->disableLegend();
  m_spectrum->setYAxisLog( wasLogY );
  //m_spectrum->showYAxisScalers( true );
  //m_spectrum->setShowPeakLabel( 0, true ); //SpectrumChart::PeakLabels::kShowPeakUserLabel
  //m_spectrum->setReferncePhotoPeakLines( ReferenceLineInfo() );

  m_stack->addWidget( m_spectrum );
  
  WContainerWidget *stackWrapper = new WContainerWidget();
  stackWrapper->addStyleClass( "RefSpectraStackWrapper" );
  WGridLayout *stackLayout = new WGridLayout();
  stackLayout->setContentsMargins( 0, 0, 0, 0 );
  stackLayout->setVerticalSpacing( 0 );
  stackLayout->setHorizontalSpacing( 0 );
  stackWrapper->setLayout( stackLayout ); 
  stackLayout->addWidget( m_stack, 0, 0 );

  WContainerWidget *stackOptions = new WContainerWidget();
  stackOptions->addStyleClass( "RefSpectraOptions RefSpecDispOptions" );
  stackLayout->addWidget( stackOptions, 1, 0 );
  stackLayout->setRowStretch( 0, 1 );


  if( portrait )
    mainLayout->addWidget( stackWrapper, 1, 0 );
  else
    mainLayout->addWidget( stackWrapper, 0, 1 );
  
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  // Add directory button
  m_addDirButton = new Wt::WContainerWidget( treeViewOptions );
  m_addDirButton->clicked().connect( this, &RefSpectraWidget::startAddDirectory );
  m_addDirButton->addStyleClass( "AddDirBtn" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec );
  HelpSystem::attachToolTipOn( m_addDirButton, WString::tr("rs-tt-add-dir-btn"), showToolTips );

  m_deleteDirButton = new Wt::WContainerWidget( treeViewOptions );
  m_deleteDirButton->clicked().connect( this, &RefSpectraWidget::removeCurrentlySelectedDirectory );
  m_deleteDirButton->addStyleClass( "DeleteDirBtn" );
  m_deleteDirButton->setHiddenKeepsGeometry( true );
  m_deleteDirButton->setHidden( true );
  HelpSystem::attachToolTipOn( m_deleteDirButton, WString::tr("rs-tt-delete-dir-btn"), showToolTips );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

  m_showCurrentForeground = new Wt::WCheckBox( WString::tr(portrait ? "rs-show-foreground-portrait" : "rs-show-foreground"), stackOptions );
  m_showCurrentForeground->addStyleClass( "CbNoLineBreak" );
  m_showCurrentForeground->changed().connect( this, &RefSpectraWidget::updatePreview );
  

  m_refBackground = new Wt::WCheckBox( WString::tr(portrait ? "rs-show-background-portrait" : "rs-show-background"), stackOptions );
  m_refBackground->addStyleClass( "CbNoLineBreak" );
  m_refBackground->changed().connect( this, &RefSpectraWidget::updatePreview );
  m_refBackground->setHiddenKeepsGeometry( true );
  m_refBackground->setHidden( true );
  
  if( portrait )
  {
    //mainLayout->setRowStretch( 0, 1 );
    //mainLayout->setRowStretch( 1, 1 );
    mainLayout->setRowResizable( 0, true, WLength(50,WLength::Percentage) ); //This is just to make the tree-view and chart-area split the heiht 50/50
  }else
  {
    mainLayout->setRowStretch( 0, 1 );
  }

  WContainerWidget *optionsSpacer = new WContainerWidget( treeViewOptions ); 
  optionsSpacer->addStyleClass( "RefSpectraOptionsSpacer" );

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  #ifdef _WIN32
  const char *show_exlorer_key = "rs-show-in-explorer";
  #elif defined( __APPLE__ )
  const char *show_exlorer_key = "rs-show-in-finder";
  #else
  const char *show_exlorer_key = "rs-show-in-linux-file-manager";
  #endif
  m_showInExplorerButton = new Wt::WPushButton( WString::tr(show_exlorer_key), treeViewOptions );

  m_showInExplorerButton->addStyleClass( "LinkBtn ShowInExplorerBtn" );
  m_showInExplorerButton->clicked().connect( this, &RefSpectraWidget::showInExplorer );
  m_showInExplorerButton->setHiddenKeepsGeometry( true );
  m_showInExplorerButton->setHidden( true );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT


  optionsSpacer = new WContainerWidget( stackOptions ); 
  optionsSpacer->addStyleClass( "RefSpectraOptionsSpacer" );

  WContainerWidget *showAsGroup = new WContainerWidget(stackOptions);
  showAsGroup->addStyleClass( "RefSpectraShowAsGroup" );
  WLabel *showAsLabel = new WLabel( WString::tr("rs-show-as"), showAsGroup );
  showAsLabel->addStyleClass( "RefSpectraShowAsLabel" );
  m_showAsComboBox = new WComboBox( showAsGroup );
  m_showAsComboBox->addStyleClass( "RefSpectraShowAsComboBox" );

  for( const SpecUtils::SpectrumType type : ns_spectrum_types )
  {
    switch( type )
    {
      case SpecUtils::SpectrumType::Foreground:
        m_showAsComboBox->addItem( WString::tr("Foreground") );
        break;
      case SpecUtils::SpectrumType::Background:
        m_showAsComboBox->addItem( WString::tr("Background") );
        break;
      case SpecUtils::SpectrumType::SecondForeground:
        m_showAsComboBox->addItem( WString::tr("Secondary") );
        break;
    }//switch( type )
  }//for( const SpecUtils::SpectrumType type : types )
  
  static_assert( ns_spectrum_types_size == 3, "ns_spectrum_types must have 3 elements" );
  m_showAsComboBox->setCurrentIndex( 2 );
  m_showAsComboBox->setFocus( true );
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
  okBtn->disable();

  DirectorySelector *dirSelector = new DirectorySelector( contents );
  
  dirSelector->pathValidityChanged().connect( 
    boost::bind( &WPushButton::setDisabled, okBtn, 
      boost::bind(std::logical_not<bool>(), boost::placeholders::_1) 
  ) );
  
  okBtn->clicked().connect( std::bind([dirSelector, this](){
    std::string path = dirSelector->path();
    this->addDirectory( path );
  }) );

  dialog->show();
}//void startAddDirectory()


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

#if BUILD_AS_OSX_APP
    // For macOS, create and store a security scoped bookmark
    // Use a deterministic key based on the normalized path (replacing path separators)
    std::string pathKey = normed_path;
    std::replace(pathKey.begin(), pathKey.end(), '/', '_');
    std::replace(pathKey.begin(), pathKey.end(), '\\', '_');
    std::replace(pathKey.begin(), pathKey.end(), ':', '_');
    std::string bookmarkKey = "SecBookmark_RefSpectra_" + pathKey;
    
    // Check if we already have a bookmark for this path
    std::vector<std::string> existingKeys = macOsUtils::getSecurityScopedBookmarkKeys("SecBookmark_RefSpectra_");
    bool bookmarkExists = false;
    for( const auto &key : existingKeys )
    {
      std::string existingPath;
      if( macOsUtils::resolveAndStartAccessingSecurityScopedBookmark(key, existingPath) )
      {
        macOsUtils::stopAccessingSecurityScopedBookmark(key); // Just checking, so stop immediately
        std::string normedExisting;
        try
        {
          normedExisting = SpecUtils::lexically_normalize_path( existingPath );
        }catch( const std::exception &e )
        {
          normedExisting = existingPath;
        }
        
        if( normedExisting == normed_path )
        {
          bookmarkExists = true;
          break;
        }
      }
    }
    
    if( bookmarkExists )
    {
      passMessage( WString::tr("rs-duplicate-dir-path").arg( path ), WarningWidget::WarningMsgMedium );
      return;
    }
    
    // Create the security scoped bookmark
    if( !macOsUtils::createAndStoreSecurityScopedBookmark(normed_path, bookmarkKey) )
    {
      passMessage( WString::tr("rs-error-bookmark-failed").arg( path ), WarningWidget::WarningMsgHigh );
      std::cerr << "RefSpectraWidget::addDirectory() - Failed to create security scoped bookmark for: " << path << std::endl;
      return;
    }
    
    passMessage( WString::tr("rs-dir-added-successfully").arg( path ), WarningWidget::WarningMsgInfo );
    
#else
    // For non-macOS platforms, use the old method
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
      std::cerr << "RefSpectraWidget::addDirectory() - Duplicate directory: " << path << std::endl;
    }
#endif
  }else
  {
    passMessage( WString::tr("rs-invalid-dir-path").arg( path ), WarningWidget::WarningMsgHigh );
    std::cerr << "RefSpectraWidget::addDirectory() - Invalid directory: " << path << std::endl;
  }
}//void RefSpectraWidget::addDirectory( const std::string &path )

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

#if BUILD_AS_OSX_APP
  // For macOS, find and remove the security scoped bookmark
  std::vector<std::string> bookmarkKeys = macOsUtils::getSecurityScopedBookmarkKeys("SecBookmark_RefSpectra_");
  bool foundBookmark = false;
  
  for( const auto &key : bookmarkKeys )
  {
    std::string bookmarkPath;
    if( macOsUtils::resolveAndStartAccessingSecurityScopedBookmark(key, bookmarkPath) )
    {
      macOsUtils::stopAccessingSecurityScopedBookmark(key); // Just checking, so stop immediately
      std::string normedBookmarkPath;
      try
      {
        normedBookmarkPath = SpecUtils::lexically_normalize_path( bookmarkPath );
      }catch( const std::exception &e )
      {
        normedBookmarkPath = bookmarkPath;
      }
      
      if( normedBookmarkPath == normed_path )
      {
        macOsUtils::removeSecurityScopedBookmark(key);
        foundBookmark = true;
        break;
      }
    }
  }
  
  if( !foundBookmark )
  {
    passMessage( WString::tr("rs-dir-not-found").arg( path ), WarningWidget::WarningMsgHigh );
    std::cerr << "RefSpectraWidget::removeCurrentlySelectedDirectory() - Directory bookmark not found: " << path << std::endl;
    return;
  }
  
  passMessage( WString::tr("rs-dir-removed-successfully").arg( path ), WarningWidget::WarningMsgInfo );
  
#else
  // For non-macOS platforms, use the old method
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

    if( (normed_dir != normed_path) || found_path_to_exclude ) //User may have added the same path multiple times, so only remove first instance
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
#endif
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
    m_lastCollapsedIndex = WModelIndex();
    m_lastExpandedIndex = WModelIndex();
    m_lastExpandTime = std::chrono::time_point<std::chrono::steady_clock>{};

#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    m_deleteDirButton->setHidden( true );
    m_showInExplorerButton->setHidden( true );
#endif
    m_refBackground->setHidden( true );
    m_showCurrentForeground->setHidden( true );
    //m_showAsComboBox->setHidden( true );
    m_stack->setCurrentIndex( 0 );

    m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::None );
    return;
  }//if( selectedIndexes.empty() )


  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
  {
    m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::None );
    return;
  }

  const bool isDir = m_treeModel->isDirectory( index );
  const std::string filePath = m_treeModel->getFilePath( index );

  m_dirInfoContainer->clear();

  if( isDir )
  {
    // We'll auto-expand directories when selected, but we'll do it after the current
    //  event loop iteration, incase the user is collapsing the node (in which case we dont want to expand it)
    boost::function<void()> do_expand = wApp->bind( boost::bind(&RefSpectraWidget::tryExpandNode, this, index) );
    const auto doExpand = [do_expand](){ do_expand(); wApp->triggerUpdate(); };
    WServer::instance()->post( wApp->sessionId(), doExpand );
    
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
    int parent_levels = 0;
    WModelIndex parent_index = index;
    while( parent_index.parent().isValid() && (parent_levels < 11) )
    {
      parent_index = parent_index.parent();
      parent_levels += 1;
      assert( parent_levels < 10 );
    }//while( parent_index.parent().isValid() )

    if( parent_index.isValid() )
    {
      const std::string parent_path = m_treeModel->getFilePath( parent_index );
      const std::string parent_path_normed = SpecUtils::lexically_normalize_path( parent_path );
      
      bool is_pref_dir = false;
      const string static_data_dir = SpecUtils::lexically_normalize_path( InterSpec::staticDataDirectory() );
      if( SpecUtils::starts_with( parent_path_normed, static_data_dir.c_str() ) )
        is_pref_dir = true;
  
  #if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
      const std::string user_dir = SpecUtils::lexically_normalize_path( InterSpec::writableDataDirectory() );
      if( SpecUtils::starts_with( parent_path_normed, user_dir.c_str() ) )
        is_pref_dir = true;
  #endif
      m_deleteDirButton->setHidden( is_pref_dir );
    }//if( parent_index.isValid() )

    m_showInExplorerButton->setHidden( false );
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

    m_refBackground->setHidden( true );
    m_showCurrentForeground->setHidden( true );
    //m_showAsComboBox->setHidden( true );

    m_stack->setCurrentIndex( 1 );
    
    Wt::WText *title_text = new Wt::WText( m_dirInfoContainer );
    title_text->addStyleClass( "DirInfoTitle" );

    int num_spectra_files = 0, num_dirs = 0;
    
    try
    {
      // We can get an exception if file path isnt allowed
      for( fs::directory_iterator itr( filePath ); itr != fs::directory_iterator(); ++itr )
      {
        const std::string file = itr->path().string();
        if( fs::is_directory( itr->path() ) )
        {
          num_dirs += 1;
          continue;
        }
        
        const string fname = SpecUtils::filename( file );
        if( SpecUtils::iequals_ascii( fname, "readme.txt" ) || SpecUtils::iequals_ascii( fname, "readme.xml" ) )
        {
          string file_data;
          if( SpecUtils::file_size( file ) < 10240 )
          {
            try
            {
              std::vector<char> data;
              SpecUtils::load_file_data( file.c_str(), data );
              file_data = std::string( data.data(), data.size() ? data.size() - 1 : 0 );
            }catch( const std::exception &e )
            {
              std::cerr << "RefSpectraWidget::handleSelectionChanged() - Error loading file: " << file << std::endl;
            }
          }//if( SpecUtils::file_size( file ) < 10*1024 )
          
          Wt::WText *text = new Wt::WText( file_data, Wt::XHTMLText, m_dirInfoContainer );
          text->addStyleClass( "DirReadmeText" );
        }else
        {
          num_spectra_files += !SpecUtils::likely_not_spec_file( file );
        }
      }
    }catch( std::exception &e )
    {
      cerr << "RefSpectraWidget: caught exception iterating directories: " << e.what() << endl;
    }//try / catch

    WString title_txt;
    if( num_dirs == 0 )
      title_txt = WString::tr("rs-dir-info-only-spectra").arg( num_spectra_files );
    else if( num_spectra_files == 0 )
      title_txt = WString::tr("rs-dir-info-only-subdir").arg( num_dirs );
    else
      title_txt = WString::tr("rs-dir-info-with-subdir-spectra").arg( num_spectra_files ).arg( num_dirs );
    title_text->setText( title_txt );

    m_spectrum->setData( nullptr, false );
    m_spectrum->setSecondData( nullptr );
    m_spectrum->setBackground( nullptr );
      
    m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::Directory );

    return;
  }//if( isDir )

  m_lastCollapsedIndex = WModelIndex();

  m_refBackground->setHidden( true );
#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
  m_deleteDirButton->setHidden( true );
  m_showInExplorerButton->setHidden( false );
#endif
  
  m_showCurrentForeground->setHidden( false );
  m_showAsComboBox->setHidden( false ); //JIC

  m_stack->setCurrentIndex( 2 );

  updatePreview();

  m_fileSelectionChangedSignal.emit( RefSpectraWidgetSelectionType::File );
}//void RefSpectraWidget::handleSelectionChanged()


void RefSpectraWidget::handleCollapsed( const Wt::WModelIndex &index )
{
  // Check if this collapse is almost immediately after an exand of the node, and if so
  //  undo the collapse - this happens when the user uses the "+" icon to expand a node;
  //  not sure if its soemthing we are doing, or a Wt issue, but this works around it.
  if( (index == m_lastExpandedIndex)
     && ((chrono::steady_clock::now() - m_lastExpandTime) < chrono::milliseconds(500)) )
  {
    m_treeView->expand( index );
    m_lastCollapsedIndex = WModelIndex();
    return;
  }//

  m_lastCollapsedIndex = index;

  // The collapse wasnt a near instant collapse after exanding - we can clear ou the following
  m_lastExpandedIndex = WModelIndex();
  m_lastExpandTime = std::chrono::time_point<std::chrono::steady_clock>{};
}//void handleCollapsed( const Wt::WModelIndex &index )


void RefSpectraWidget::handleExpanded( const Wt::WModelIndex &index )
{
  m_lastExpandedIndex = index;
  m_lastExpandTime = std::chrono::steady_clock::now();
}


void RefSpectraWidget::tryExpandNode( const Wt::WModelIndex &index )
{
  if( m_lastCollapsedIndex != index )
  {
    m_treeView->expand( index );
    m_lastCollapsedIndex = WModelIndex();
  }
}


void RefSpectraWidget::handleLayoutAboutToBeChanged()
{
  // Nothing to do here?
}


void RefSpectraWidget::handleLayoutChanged()
{
  const WModelIndexSet selectedIndexes = m_treeView->selectedIndexes();
  if( selectedIndexes.size() != 1 )
    return;

  int parent_levels = 0;
  WModelIndex index = *selectedIndexes.begin();
  std::vector<WModelIndex> nodes_to_expand;

  while( index.isValid() && (parent_levels < 12) )
  {
    nodes_to_expand.push_back( index );
    index = index.parent();
    parent_levels += 1;
    assert( parent_levels < 10 );
  }//while( parent_index.parent().isValid() )

  std::reverse( nodes_to_expand.begin(), nodes_to_expand.end() );
  for( const WModelIndex &node : nodes_to_expand )
    m_treeView->expand( node );

  cout << "handleLayoutChanged() - expanded " << nodes_to_expand.size() << " nodes" << endl;
}


void RefSpectraWidget::updatePreview()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  shared_ptr<const SpecUtils::Measurement> user_foreground = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  m_showCurrentForeground->setHidden( !user_foreground );

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
    []( const std::string &fname, void *userdata ) -> bool{
      return SpecUtils::starts_with( SpecUtils::filename(fname), "background." );
  }, nullptr );

  SpecMeas infile;
  const bool loaded = infile.load_file( filePath, SpecUtils::ParserType::Auto, extension );
  if( !loaded )
  {
    m_stack->setCurrentIndex( 1 );
    WString error_text = WString::tr("rs-error-invalid-file").arg( filename );
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
  {
    foreground_meas = infile.measurements().front();
  }

  if( !background_meas && !background_files.empty() )
  {
    const string background_file = background_files.front();
    const string background_extension = SpecUtils::file_extension( background_file );
    SpecMeas background_infile;
    const bool background_loaded = background_infile.load_file( background_file, SpecUtils::ParserType::Auto, background_extension );

    // TODO: properly check for background records...
    if( background_loaded )
      background_meas = background_infile.measurements().front();
  }//if( !background_meas && !background_files.empty() )
  
  m_refBackground->setHidden( !background_meas );

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
  for( size_t i = 0; i < ns_spectrum_types_size; ++i )
  {
    if( ns_spectrum_types[i] == type )
    {
      m_showAsComboBox->setCurrentIndex( static_cast<int>( i ) );
      return;
    }
  }//for( size_t i = 0; i < ns_spectrum_types_size; ++i )

  assert( 0 ); //shouldnt ever happen
}//void setLoadSpectrumType( SpecUtils::SpectrumType type )


//void RefSpectraWidget::selectSimilarToForeground()
//{
//  // TODO: Implement ...
//}//void selectSimilarToForeground();


void RefSpectraWidget::loadSelectedSpectrum()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  const WModelIndexSet selectedIndexes = m_treeView->selectedIndexes();
  assert( selectedIndexes.empty() || (selectedIndexes.size() == 1) );
  if( selectedIndexes.empty() )
  {
    assert( 0 );
    return;
  }

  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
  {
    assert( 0 );
    return;
  }

  const bool isDir = m_treeModel->isDirectory( index );
  assert( !isDir );
  if( isDir )
  {
    assert( 0 );
    return;
  }

  const std::string filePath = m_treeModel->getFilePath( index );
  const std::string displayName = m_treeModel->getDisplayName( index );
  
  auto meas = make_shared<SpecMeas>();
  const bool success = meas->load_file( filePath, SpecUtils::ParserType::Auto, displayName );
    
  if( !success )
  {
    assert( 0 );
    passMessage( WString::tr("rs-error-loading-file").arg(displayName), WarningWidget::WarningMsgHigh );
    return;
  }


  const int current_index = m_showAsComboBox->currentIndex();
  assert( (current_index >= 0) && (current_index < static_cast<int>( ns_spectrum_types_size )) );
  
  SpecUtils::SpectrumType type = SpecUtils::SpectrumType::SecondForeground;
  if( (current_index >= 0) && (current_index < static_cast<int>( ns_spectrum_types_size )) )
    type = ns_spectrum_types[current_index];
  else
    passMessage( WString::tr("rs-error-invalid-type-selected"), WarningWidget::WarningMsgHigh );
  
 
  shared_ptr<SpectraFileHeader> header = std::make_shared<SpectraFileHeader>( true, interspec );
  header->setFile( displayName, meas );
  SpecMeasManager *fileManager = interspec->fileManager();
  SpectraFileModel *fileModel = fileManager->model();
  const int row = fileModel->addRow( header );

  const bool checkIfPreviouslyOpened = false;
  const bool doPreviousEnergyRangeCheck = false;
  interspec->fileManager()->displayFile( row, meas, type, 
               checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
               SpecMeasManager::VariantChecksToDo::None );
}//void loadSelectedSpectrum()


Wt::Signal<RefSpectraWidgetSelectionType> &RefSpectraWidget::fileSelectionChangedSignal() 
{ 
  return m_fileSelectionChangedSignal; 
}


void RefSpectraWidget::initBaseDirs()
{
  const std::string static_data_dir = InterSpec::staticDataDirectory();
  const std::string ref_spectra_dir = SpecUtils::append_path( static_data_dir, "reference_spectra" );

  const std::string common_nucs_dir = SpecUtils::append_path( ref_spectra_dir, "Common_Field_Nuclides" );

  if( SpecUtils::is_directory( common_nucs_dir ) )
  {
    const WModelIndex index = m_treeModel->addBaseDirectory( common_nucs_dir );
    assert( index.isValid() );
    if( index.isValid() )
    {
      cout << "Expanded common_nucs_dir" << endl;
      m_treeView->expand( index );
    }else
    {
      cout << "Failed to expand base directory: " << common_nucs_dir << endl;
    }
  }//if( SpecUtils::is_directory( common_nucs_dir ) )


#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  // Check for directories under reference_spectra in the "user_data" directory
  const std::string user_dir = InterSpec::writableDataDirectory();
  const std::string user_ref_spectra_dir = SpecUtils::append_path( user_dir, "reference_spectra" );
  if( SpecUtils::is_directory( user_ref_spectra_dir ) ) 
  {
    bool has_spectra_files_in_base = false;
    vector<string> dirs_to_add;
    try
    {
      for( const auto &entry : fs::directory_iterator( user_ref_spectra_dir ) )
      {
        const std::string entry_name = entry.path().filename().string();
        if( entry_name.empty() || entry_name.front() == '.' )
          continue;
        
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
    }catch( std::exception &e )
    {
      cerr << "RefSpectraWidget::initBaseDirs(): caught exception iterating directories: " << e.what() << endl;
    }//try / catc

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
  
#if BUILD_AS_OSX_APP
  // For macOS, resolve security scoped bookmarks
  std::vector<std::string> bookmarkKeys = macOsUtils::getSecurityScopedBookmarkKeys("SecBookmark_RefSpectra_");
  std::vector<std::string> failedBookmarks;
  
  for( const auto &key : bookmarkKeys )
  {
    std::string resolvedPath;
    if( macOsUtils::resolveAndStartAccessingSecurityScopedBookmark(key, resolvedPath) )
    {
      if( SpecUtils::is_directory( resolvedPath ) )
      {
        m_treeModel->addBaseDirectory( resolvedPath );
        // Note: We keep the security scoped resource access active during the app session
      }
      else
      {
        std::cerr << "RefSpectraWidget::initBaseDirs() - Resolved bookmark path is not a directory: " << resolvedPath << std::endl;
        macOsUtils::stopAccessingSecurityScopedBookmark(key);
        failedBookmarks.push_back(key);
      }
    }
    else
    {
      std::cerr << "RefSpectraWidget::initBaseDirs() - Failed to resolve security scoped bookmark: " << key << std::endl;
      failedBookmarks.push_back(key);
    }
  }
  
  // Clean up failed bookmarks
  for( const auto &key : failedBookmarks )
  {
    macOsUtils::removeSecurityScopedBookmark(key);
  }
  
  if( !failedBookmarks.empty() )
  {
    const int numdir = static_cast<int>( failedBookmarks.size() );
    passMessage( WString::tr("rs-some-dirs-removed").arg(numdir), WarningWidget::WarningMsgMedium );
  }
  
#else
  // For non-macOS platforms, use the old method
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
#endif
#endif // !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT

  // Sort it so things are in alphabetical order
  m_treeView->sortByColumn( 0, Wt::SortOrder::AscendingOrder );


  // Expand the "Common Field Nuclides" node
  for( int row = 0; row < m_treeModel->rowCount(); ++row )
  {
    const WModelIndex index = m_treeModel->index( row, 0 );
    std::string name = m_treeModel->getDisplayName( index );
    cout << "name: " << name << endl;
    if( name == "Common Field Nuclides" )
    {
      m_treeView->expand( index );
      m_treeView->setSelectedIndexes( {index} );
      cout << "Expanded and selected Common Field Nuclides" << endl;
      break;
    }
  }//for( int row = 0; row < m_treeModel->rowCount(); ++row )
}//void initBaseDirs();


#if( !IOS && !ANDROID && !BUILD_FOR_WEB_DEPLOYMENT )
void RefSpectraWidget::showInExplorer()
{
  const WModelIndexSet selectedIndexes = m_treeView->selectedIndexes();
  assert( selectedIndexes.empty() || (selectedIndexes.size() == 1) );
  if( selectedIndexes.empty() )
    return;

  const WModelIndex index = *selectedIndexes.begin();
  assert( index.isValid() );
  if( !index.isValid() )
    return;
    
  const std::string filePath = m_treeModel->getFilePath( index );
  AppUtils::showFileInOsFileBrowser(filePath);
}//void showInExplorer()
#endif

/*
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"


void RefSpectraWidget::dev_code()
{
  const string output_path = "/path/to/output/RefSpectra";
  const string base_path = "/path/to/input/peak_area_optimization/peak_fit_accuracy_inject/";
  const string city_time = "Livermore/1800_seconds";

  const vector<pair<string,string>> replacements{
    { "_ShHeavy", "_HeavyShielding" },
    { "_Unsh", "_Unshielded" },
    { "_Sh", "_Shielded" },
  };

  const vector<pair<SpecUtils::DetectorType,string>> det_type_replacements{
    { SpecUtils::DetectorType::DetectiveEx, "Detective-EX" },
    { SpecUtils::DetectorType::DetectiveX, "Detective-X" },
    { SpecUtils::DetectorType::IdentiFinderR500NaI, "IdentiFINDER-R500-NaI" },
    { SpecUtils::DetectorType::RadSeekerLaBr, "Radseeker-LaBr3" },
    { SpecUtils::DetectorType::IdentiFinderNG, "IdentiFINDER-NGH" },
    { SpecUtils::DetectorType::OrtecRadEagleNai, "RadEagle" },
    { SpecUtils::DetectorType::RadiaCode, "Kromek-GR1-CZT" },
    { SpecUtils::DetectorType::Fulcrum40h, "Fulcrum40h" },
    { SpecUtils::DetectorType::Unknown, "CZT_H3D_M400_ORNL_25cm" },
    { SpecUtils::DetectorType::Unknown, "1cm-1cm-1cm" },
    { SpecUtils::DetectorType::Unknown, "LaBr3_1.5x1.5_SNL" },
    { SpecUtils::DetectorType::RadiaCode, "Radiacode-102" },
    { SpecUtils::DetectorType::AvidRsi, "Rack" },
    { SpecUtils::DetectorType::VerifinderNaI, "Verifinder-SN23N" }
  };

  for( const auto &dettype_dirname : det_type_replacements )
  {
    const SpecUtils::DetectorType &det_type = dettype_dirname.first;
    std::string det_type_name = SpecUtils::detectorTypeToString( det_type );
    if( det_type == SpecUtils::DetectorType::Unknown )
      det_type_name = dettype_dirname.second;

    const std::string &dirname = dettype_dirname.second;
    const string in_file_dir = SpecUtils::append_path( SpecUtils::append_path( base_path, dirname ), city_time );
    const string det_out_dir = SpecUtils::append_path( output_path, det_type_name );

    if( !SpecUtils::is_directory(in_file_dir) )
    {
      cerr << "Failed to find directory '" << in_file_dir << "' :(" << endl;
      return;
    }

    if( SpecUtils::create_directory(det_out_dir) == 0 )
    {
      cerr << "Failed to create directory '" << det_out_dir << "' :(" << endl;
      return;
    }


    // loop over all .pcf files in the directory
    bool wrote_back = false;
    vector<string> pcf_files = SpecUtils::ls_files_in_directory(in_file_dir, ".pcf");
    for( const string &filename : pcf_files )
    {
      SpecUtils::SpecFile spec;
      const bool success = spec.load_file( filename, SpecUtils::ParserType::Auto, filename );
      if( !success )
      {
        cerr << "Failed to load file '" << filename << "' :(" << endl;
        return;
      }

      if( spec.num_measurements() < 3 )
      {
        cerr << "File '" << filename << "' has less than 3 measurements :(" << endl;
        return;
      }

      if( det_type == SpecUtils::DetectorType::Unknown )
      {
        spec.set_instrument_model( det_type_name );
      }else
      {
        spec.set_detector_type( det_type );
        spec.set_instrument_model( det_type_name );
      }

      std::shared_ptr<const SpecUtils::Measurement> foreground = spec.measurement( size_t(0) );
      std::shared_ptr<const SpecUtils::Measurement> background = spec.measurement( size_t(1) );
      assert( background && foreground );

      SpecUtils::SpecFile out_foreground = spec, out_background = spec;
      out_foreground.remove_measurements( out_foreground.measurements() );
      out_background.remove_measurements( out_background.measurements() );

      auto fore = std::make_shared<SpecUtils::Measurement>(*foreground);
      auto back = std::make_shared<SpecUtils::Measurement>(*background);

      out_foreground.add_measurement( fore, true );
      out_background.add_measurement( back, true );


      string outname = SpecUtils::filename(filename);
      const string orig_extension = SpecUtils::file_extension(outname);
      outname = outname.substr( 0, outname.size() - orig_extension.size() );
      for( const auto &from_to : replacements )
        SpecUtils::ireplace_all( outname, from_to.first.c_str(), from_to.second.c_str() );

      string title = outname;
      SpecUtils::ireplace_all( title, "_", " " );

      outname += ".txt";
      outname = SpecUtils::append_path( det_out_dir, outname );

      back->set_title( "Background" );
      fore->set_title( title );


      ofstream out_stream( outname.c_str(), ios::out | ios::binary );
      if( !out_stream )
      {
        cerr << "Could not open file " << outname << " for writing" << endl;
        return;
      }

      const bool wrote_fore = out_foreground.write_uri( out_stream, 1, 0 );
      if( !wrote_fore )
      {
        cerr << "Could not write foreground spectrum to file " << outname << endl;
        return;
      }
      cout << "Wrote '" << outname << "'" << endl;

      if( !wrote_back )
      {
        string backname = SpecUtils::append_path( det_out_dir, "background.txt" );
        ofstream back_stream( backname.c_str(), ios::out | ios::binary );
        if( !back_stream )
        {
          cerr << "Could not open file " << backname << " for writing" << endl;
          return;
        }
        const bool wrote_background = out_background.write_uri( back_stream, 1, 0 );
        if( !wrote_background )
        {
          cerr << "Could not write background spectrum to file " << outname << endl;
          return;
        }

        cout << "Wrote '" << backname << "'" << endl;
        wrote_back = true;
      }//if( !wrote_back
    }//for( const string &filename : pcf_files )
  }//for( const auto &dettype_dirname : det_type_replacements )
}//class RefSpectraWidget
*/
