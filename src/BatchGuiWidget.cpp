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
#include "InterSpec/BatchGuiWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/FileDragUploadResource.h"


#if( BUILD_AS_ELECTRON_APP )
#include "target/electron/ElectronUtils.h"
#endif

#if( BUILD_AS_OSX_APP )
#include "target/osx/macOsUtils.h"
#endif

using namespace Wt;
using namespace std;

namespace fs = std::filesystem;

namespace 
{
  
}// anonymous namespace


BatchGuiDialog::BatchGuiDialog( FileDragUploadResource *uploadResource, const Wt::WString &title )
  : SimpleDialog( title ),
    m_widget( nullptr ),
    m_processBtn( nullptr )
{
  addStyleClass( "BatchGuiDialog" );

  WGridLayout *layout = new WGridLayout();
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );

  contents()->setLayout( layout );

  m_widget = new BatchGuiWidget( uploadResource );
  layout->addWidget( m_widget, 0, 0 );

  m_processBtn = addButton( WString::tr("Analyze") );
  //m_processBtn->clicked().connect( m_widget, &BatchGuiWidget::doSomething );
  m_processBtn->setDisabled( true );
  m_widget->canDoAnalysis().connect( boost::bind( &WPushButton::setEnabled, m_processBtn, boost::placeholders::_1 ) );

  WPushButton *cancelBtn = addButton( WString::tr("Close") );
  cancelBtn->clicked().connect( m_widget, &BatchGuiWidget::handleClose );
  
  // Set dialog size based on screen size
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  if( interspec && interspec->isPhone() )
  {
    addStyleClass( "BatchGuiDialog-phone" );
#if( IOS )
    addStyleClass( "BatchGuiDialog-iphone" );
#endif
  }

  // We want to target 750px x 500px for normal - portrait phones, we'll take what we can get.
  //  (note: typical, old 7" tablets are at least 600x1024)
  const bool portrait = ((interspec->renderedWidth() > 100) 
                         && (interspec->renderedHeight() > 100)
                         && (interspec->renderedWidth() < 480));
  if( portrait )
  {
    addStyleClass( "BatchGuiDialog-portrait" );
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

BatchGuiDialog::~BatchGuiDialog()
{
  // The widget will be automatically deleted by Wt
} 


BatchGuiDialog *BatchGuiDialog::createDialog( FileDragUploadResource *uploadResource )
{
  //WString title = WString::tr("bgw-dialog-title");
  WString title;
  BatchGuiDialog *dialog = new BatchGuiDialog( uploadResource, title );

  return dialog;
}

BatchGuiWidget::BatchGuiWidget( FileDragUploadResource *uploadResource, Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
    m_uploadResource( uploadResource ),
    m_canDoAnalysis( this ),
    m_inputFiles{}
{
  assert( m_uploadResource );
  
  m_uploadResource->fileDrop().connect( this, &BatchGuiWidget::handleFileDrop );
}//BatchGuiWidget constructor


BatchGuiWidget::~BatchGuiWidget()
{
  for( const tuple<std::string,std::string,bool> &file : m_inputFiles )
  {
    const bool should_delete = std::get<2>(file);
    if( should_delete )
    {
      const string &path_to_file = std::get<1>(file);
      const bool success = SpecUtils::remove_file( path_to_file );
      if( !success )
        cerr << "FileDragUploadResource: Warning, could not delete file '" << path_to_file << "'" << endl;
    }//if( should_delete )
  }//for( const tuple<std::string,std::string,bool> &file : m_inputFiles )
  
  m_inputFiles.clear();
}//~BatchGuiWidget()


void BatchGuiWidget::handleClose()
{
  m_canDoAnalysis.emit( false );
  m_uploadResource->clearSpooledFiles();
}


Wt::Signal<bool> &BatchGuiWidget::canDoAnalysis()
{
  return m_canDoAnalysis;
}


void BatchGuiWidget::handleFileDrop( const std::string &, const std::string & )
{
  addInputFiles( m_uploadResource->takeSpooledFiles() );
}


void BatchGuiWidget::addInputFiles( const std::vector<std::tuple<std::string,std::string,bool>> &files )
{
  for( const tuple<string,string,bool> &file : files )
  {
    m_inputFiles.push_back( file );
  }
  
  updateCanDoAnalysis();
}

void BatchGuiWidget::updateCanDoAnalysis()
{
  m_canDoAnalysis.emit( !m_inputFiles.empty() );
}
