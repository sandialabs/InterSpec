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
#include <chrono>

#include <Wt/WMenu>
#include <Wt/WGroupBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>


#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/BatchGuiWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/BatchGuiAnaWidget.h"
#include "InterSpec/BatchGuiInputFile.h"
#include "InterSpec/DirectorySelector.h"
#include "InterSpec/FileDragUploadResource.h"


using namespace Wt;
using namespace std;


BatchGuiDialog::BatchGuiDialog( FileDragUploadResource *uploadResource, const Wt::WString &title )
: SimpleDialog( title ), m_widget( nullptr ), m_processBtn( nullptr )
{
  addStyleClass( "BatchGuiDialog" );

  WGridLayout *layout = new WGridLayout();
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );

  contents()->setLayout( layout );

  m_widget = new BatchGuiWidget( uploadResource );
  layout->addWidget( m_widget, 0, 0 );

  m_processBtn = new WPushButton( WString::tr( "bgw-analyze-button" ), footer() );
  m_processBtn->setStyleClass( "simple-dialog-btn" );
  m_processBtn->clicked().connect( m_widget, &BatchGuiWidget::performAnalysis );
  m_processBtn->disable();
  m_widget->canDoAnalysis().connect( boost::bind( &WPushButton::setEnabled, m_processBtn, boost::placeholders::_1 ) );

  addButton( WString::tr( "Close" ) );

  // Set dialog size based on screen size
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( !interspec )
    return;

  if( interspec && interspec->isPhone() )
  {
    addStyleClass( "BatchGuiDialog-phone" );
#if ( IOS )
    addStyleClass( "BatchGuiDialog-iphone" );
#endif
  }

  // We want to target 750px x 500px for normal - portrait phones, we'll take what we can get.
  //  (note: typical, old 7" tablets are at least 600x1024)
  const bool portrait = ( ( interspec->renderedWidth() > 100 ) && ( interspec->renderedHeight() > 100 ) &&
                          ( interspec->renderedWidth() < 480 ) );
  if( portrait )
  {
    addStyleClass( "BatchGuiDialog-portrait" );
    m_widget->setWidth( 0.95 * interspec->renderedWidth() - 30 );
    m_widget->setHeight( 0.95 * interspec->renderedHeight() - 90 );
    m_widget->setMinimumSize( 0.95 * interspec->renderedWidth() - 30, 0.95 * interspec->renderedHeight() - 90 );
    m_widget->setMaximumSize( 0.95 * interspec->renderedWidth() - 30, 0.95 * interspec->renderedHeight() - 90 );
  } else if( ( interspec->renderedWidth() > 100 ) && ( interspec->renderedHeight() > 50 ) )
  {
    double dialogWidth = std::min( 750.0, 0.95 * interspec->renderedWidth() );
    double dialogHeight = std::min( 500.0, 0.95 * interspec->renderedHeight() );

    m_widget->setWidth( dialogWidth - 30 );
    m_widget->setHeight( dialogHeight - 90 );
    m_widget->setMinimumSize( dialogWidth - 30, dialogHeight - 90 );
    m_widget->setMaximumSize( dialogWidth - 30, dialogHeight - 90 );
  } else
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
  // WString title = WString::tr("bgw-dialog-title");
  WString title;
  BatchGuiDialog *dialog = new BatchGuiDialog( uploadResource, title );

  return dialog;
}

BatchGuiWidget::BatchGuiWidget( FileDragUploadResource *uploadResource, Wt::WContainerWidget *parent )
: Wt::WContainerWidget( parent ),
  m_uploadResource( uploadResource ),
  m_batch_type_menu( nullptr ),
  m_options_stack( nullptr ),
  m_act_shield_ana_opts( nullptr ),
  m_peak_fit_opts( nullptr ),
  m_input_files_container( nullptr ),
  m_output_dir( nullptr ),
  m_can_do_analysis( false ),
  m_canDoAnalysis( this )
{
  assert( m_uploadResource );
  InterSpec *interspec = InterSpec::instance();
  WApplication *app = WApplication::instance();

  app->useStyleSheet( "InterSpec_resources/BatchGuiWidget.css" );
  app->require( "InterSpec_resources/BatchGuiWidget.js" );
  interspec->useMessageResourceBundle( "BatchGuiWidget" );

  addStyleClass( "BatchGuiWidget" );

  doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', true);" );
  // doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', true);" );

  interspec->saveShieldingSourceModelToForegroundSpecMeas();
#if ( USE_REL_ACT_TOOL )
  interspec->saveRelActManualStateToForegroundSpecMeas();
  interspec->saveRelActAutoStateToForegroundSpecMeas();
#endif
  
  Wt::WGroupBox *options_container = new Wt::WGroupBox( WString::tr( "bgw-type-select-label" ), this );
  options_container->addStyleClass( "TypeSelectContainer" );

  m_options_stack = new Wt::WStackedWidget();
  m_batch_type_menu = new Wt::WMenu( m_options_stack, options_container );
  m_batch_type_menu->addStyleClass( "LightNavMenu VerticalNavMenu" );
  options_container->addWidget( m_options_stack );

  m_act_shield_ana_opts = new BatchGuiActShieldAnaWidget();
  WMenuItem *item = new WMenuItem(
    WString::tr( "bgw-act-shield-ana-opts-label" ), m_act_shield_ana_opts, WMenuItem::LoadPolicy::PreLoading );
  m_act_shield_ana_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
  m_batch_type_menu->addItem( item );

  m_peak_fit_opts = new BatchGuiPeakFitWidget();
  item = new WMenuItem( WString::tr( "bgw-peak-fit-opts-label" ), m_peak_fit_opts, WMenuItem::LoadPolicy::PreLoading );
  m_peak_fit_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
  m_batch_type_menu->addItem( item );

  m_batch_type_menu->select( 0 );

  m_input_files_container = new WGroupBox( WString::tr( "bgw-input-files-label" ), this );
  m_input_files_container->addStyleClass( "InputFilesContainer" );

  m_input_files_container->doJavaScript( "BatchInputDropUploadSetup(" + m_input_files_container->jsRef() +
                                         ", "
                                         " '" +
                                         interspec->fileManager()->batchDragNDrop()->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_input_files_container->id() + "']);" );

  m_output_dir = new DirectorySelector( this );
  m_output_dir->setLabelTxt( WString::tr( "bgw-output-dir-label" ) );
  m_output_dir->pathChanged().connect( this, &BatchGuiWidget::updateCanDoAnalysis );

  m_uploadResource->fileDrop().connect( this, &BatchGuiWidget::handleFileDrop );

  addInputFiles( m_uploadResource->takeSpooledFiles() );
}// BatchGuiWidget constructor

BatchGuiWidget::~BatchGuiWidget()
{
  wApp->doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', null);" );
  // wApp->doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', null);" );
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_input_files_container->id() + "']);" );
}//~BatchGuiWidget()

Wt::Signal<bool> &BatchGuiWidget::canDoAnalysis()
{
  return m_canDoAnalysis;
}

void BatchGuiWidget::handleFileDrop( const std::string &, const std::string & )
{
  addInputFiles( m_uploadResource->takeSpooledFiles() );
}

void BatchGuiWidget::addInputFiles( const std::vector<std::tuple<std::string, std::string, bool>> &files )
{
  for( const tuple<string, string, bool> &file : files )
  {
    const string &display_name = std::get<0>( file );
    const string &path_to_file = std::get<1>( file );
    const bool should_delete = std::get<2>( file );

    BatchGuiInputSpectrumFile *input_file =
      new BatchGuiInputSpectrumFile( display_name, path_to_file, should_delete, m_input_files_container );
    input_file->preview_created_signal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
    input_file->remove_self_request().connect(
      boost::bind( &BatchGuiWidget::handle_remove_input_file, this, boost::placeholders::_1 ) );
  }
}// void BatchGuiWidget::addInputFiles()

void BatchGuiWidget::handle_remove_input_file( BatchGuiInputSpectrumFile *input )
{
  delete input;

  updateCanDoAnalysis();
}// void handle_remove_input_file( BatchGuiInputSpectrumFile *input )

void BatchGuiWidget::updateCanDoAnalysis()
{
  size_t num_input_files = 0;

  for( Wt::WWidget *child : m_input_files_container->children() )
  {
    BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
    num_input_files += ( input_file && input_file->spec_meas() );
  }

  BatchGuiAnaWidget *const batch_ana_widget = dynamic_cast<BatchGuiAnaWidget *>( m_options_stack->currentWidget() );

  const bool can_do_analysis =
    ( num_input_files > 0 ) && m_output_dir->isPathValid() && ( batch_ana_widget && batch_ana_widget->canDoAnalysis() );
  if( can_do_analysis != m_can_do_analysis )
  {
    m_can_do_analysis = can_do_analysis;
    m_canDoAnalysis.emit( can_do_analysis );
  }
}// void updateCanDoAnalysis()

void BatchGuiWidget::performAnalysis()
{
  vector<tuple<string, string, std::shared_ptr<const SpecMeas>>> input_files;
  const string output_dir = m_output_dir->path();

  for( Wt::WWidget *child : m_input_files_container->children() )
  {
    BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
    if( input_file )
    {
      const string display_name = input_file->display_name();
      const string path_to_file = input_file->path_to_file();
      const shared_ptr<const SpecMeas> spec_meas = input_file->spec_meas();

      input_files.push_back( make_tuple( display_name, path_to_file, spec_meas ) );
    }// if( input_file )
  }// for( Wt::WWidget *child : m_input_files_container->children() )

  BatchGuiAnaWidget *batch_ana_widget = dynamic_cast<BatchGuiAnaWidget *>( m_options_stack->currentWidget() );
  assert( batch_ana_widget );
  if( batch_ana_widget )
  {
    batch_ana_widget->performAnalysis( input_files, output_dir );
  } else
  {
    cerr << "BatchGuiWidget::performAnalysis: Warning, no batch analysis widget selected" << endl;
  }
}// void BatchGuiWidget::performAnalysis()
