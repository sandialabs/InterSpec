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
#include <Wt/WCheckBox>
#include <Wt/WTreeView>
#include <Wt/WLineEdit>
#include <Wt/WMessageBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/RefSpectraModel.h"
#include "InterSpec/RefSpectraWidget.h"
#include "InterSpec/UserPreferences.h"

using namespace Wt;
using namespace std;

namespace fs = std::filesystem;


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

RefSpectraDialog *RefSpectraDialog::createDialog( SpecUtils::SpectrumType type )
{
  RefSpectraDialog *dialog = new RefSpectraDialog( WString::tr("rs-dialog-title") );
  dialog->m_widget->setLoadSpectrumType( type );
  return dialog;
}

RefSpectraWidget::RefSpectraWidget( Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
  m_treeModel( nullptr ),
  m_treeView( nullptr ),
  m_previewContainer( nullptr ),
  m_overlayCurrentDoc( nullptr ),
  m_overlayOtherDoc( nullptr ),
  m_addDirButton( nullptr ),
  m_deleteDirButton( nullptr ),
  m_loadSpectrumType( SpecUtils::SpectrumType::SecondForeground )
{
  setupUI();
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
  if( interspec )
    interspec->useMessageResourceBundle( "RefSpectraWidget" );

  WGridLayout *mainLayout = new WGridLayout();
  setLayout( mainLayout );
  
  // Tree view
  m_treeModel = new RefSpectraModel( this );
  m_treeView = new WTreeView();
  m_treeView->setModel( m_treeModel );
  m_treeView->clicked().connect( this, &RefSpectraWidget::handleFileSelection );
  
  mainLayout->addWidget( m_treeView, 0, 0 );

  // Options
  WContainerWidget *optionsContainer = new WContainerWidget();
  mainLayout->addWidget( optionsContainer, 1, 0, 1, 2 );

  // Preview area - to be replaced by a stack or something
  m_previewContainer = new WContainerWidget();
  mainLayout->addWidget( m_previewContainer, 0, 1, 1, 1 );

  // Add directory button
  m_addDirButton = new Wt::WContainerWidget( optionsContainer );
  m_addDirButton->clicked().connect( this, &RefSpectraWidget::addDirectory );
  m_addDirButton->addStyleClass( "AddDirBtn" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec );
  HelpSystem::attachToolTipOn( m_addDirButton, WString::tr("rs-tt-add-dir-btn"), showToolTips );

  m_deleteDirButton = new Wt::WContainerWidget( optionsContainer );
  m_deleteDirButton->clicked().connect( this, &RefSpectraWidget::removeDirectory );
  m_deleteDirButton->addStyleClass( "DeleteDirBtn" );
  m_deleteDirButton->setHiddenKeepsGeometry( true );
  m_deleteDirButton->setHidden( true );
  HelpSystem::attachToolTipOn( m_deleteDirButton, WString::tr("rs-tt-delete-dir-btn"), showToolTips );

  m_overlayCurrentDoc = new Wt::WCheckBox( WString::tr("rs-show-foreground"), optionsContainer );
  m_overlayCurrentDoc->addStyleClass( "CbNoLineBreak" );
  m_overlayCurrentDoc->changed().connect( this, &RefSpectraWidget::updatePreview );
  m_overlayCurrentDoc->setFloatSide( Wt::Side::Right );
  
  m_overlayOtherDoc = new Wt::WCheckBox( WString::tr("rs-show-background"), optionsContainer );
  m_overlayOtherDoc->addStyleClass( "CbNoLineBreak" );
  m_overlayOtherDoc->changed().connect( this, &RefSpectraWidget::updatePreview );
  m_overlayOtherDoc->setFloatSide( Wt::Side::Right );
  mainLayout->setRowStretch( 0, 1 );
}

void RefSpectraWidget::addDirectory()
{
  // TODO: Implement directory addition logic
}

void RefSpectraWidget::removeDirectory()
{
  // TODO: Implement directory removal logic
  //       Check if base directory is currently selected, then get directory, and remove from preferences.

  //auto it = std::find( m_baseDirectories.begin(), m_baseDirectories.end(), path );
  //if( it != m_baseDirectories.end() ) {
  //  m_baseDirectories.erase( it );
  //  m_treeModel->removeBaseDirectory( path );
  //}
}

void RefSpectraWidget::updateTreeView()
{
  m_treeModel->refresh();
}

void RefSpectraWidget::handleFileSelection( const Wt::WModelIndex &index )
{
  std::string filePath = m_treeModel->getFilePath( index );
  if( !filePath.empty() ) {
    updatePreview();
  }
}

void RefSpectraWidget::updatePreview()
{
  // TODO: Implement preview update logic based on selected file and overlay options
  m_previewContainer->clear();
  Wt::WText *text = new Wt::WText("Preview will be implemented here");
  m_previewContainer->addWidget( text );
}


void RefSpectraWidget::setLoadSpectrumType( SpecUtils::SpectrumType type )
{
  m_loadSpectrumType = type;
}