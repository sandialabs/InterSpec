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
#include <memory>
#include <iostream>
#include <chrono>

#include <Wt/WMenu.h>
#include <Wt/WServer.h>
#include <Wt/WGroupBox.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WStackedWidget.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/Filesystem.h"

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

  WGridLayout *layout = contents()->setLayout( std::make_unique<WGridLayout>() );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );

  m_widget = layout->addWidget( std::make_unique<BatchGuiWidget>( uploadResource ), 0, 0 );

  m_processBtn = footer()->addNew<WPushButton>( WString::tr( "bgw-analyze-button" ) );
  m_processBtn->setStyleClass( "simple-dialog-btn" );
  m_processBtn->clicked().connect( m_widget, &BatchGuiWidget::performAnalysis );
  m_processBtn->disable();
  m_widget->canDoAnalysis().connect( [this]( bool enabled ){ m_processBtn->setEnabled( enabled ); } );

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
    double dialogHeight = std::min( 650.0, 0.95 * interspec->renderedHeight() );

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
  const WString title;
  BatchGuiDialog *dialog = SimpleDialog::make<BatchGuiDialog>( uploadResource, title );

  return dialog;
}

BatchGuiWidget::BatchGuiWidget( FileDragUploadResource *uploadResource )
: Wt::WContainerWidget(),
  m_uploadResource( uploadResource ),
  m_batch_type_menu( nullptr ),
  m_options_stack( nullptr ),
  m_act_shield_ana_opts( nullptr ),
  m_peak_fit_opts( nullptr ),
  m_file_convert_opts( nullptr ),
  m_input_files_container( nullptr ),
  m_output_dir( nullptr ),
  m_input_status_error( nullptr ),
  m_can_do_analysis( false ),
  m_canDoAnalysis()
{
  assert( m_uploadResource );
  InterSpec *interspec = InterSpec::instance();
  WApplication *app = WApplication::instance();

  app->useStyleSheet( "InterSpec_resources/BatchGuiWidget.css" );
  app->require( "InterSpec_resources/BatchGuiWidget.js" );
  interspec->useMessageResourceBundle( "BatchGuiWidget" );

  addStyleClass( "BatchGuiWidget" );

  doJavaScript( "window._IS=window._IS||{};window._IS.BlockFileDrops=true;" );
  // doJavaScript( "window._IS.BatchUploadOnly=true;" );

  interspec->saveShieldingSourceModelToForegroundSpecMeas();
#if ( USE_REL_ACT_TOOL )
  interspec->saveRelActManualStateToForegroundSpecMeas();
  interspec->saveRelActAutoStateToForegroundSpecMeas();
#endif
  
  Wt::WGroupBox *options_container = addNew<Wt::WGroupBox>( WString::tr( "bgw-type-select-label" ) );
  options_container->addStyleClass( "TypeSelectContainer" );

  // Note: In Wt3 menu was added to options_container via constructor, stack added after.
  //  In Wt4 we replicate that order: menu first, then stack.
  {
    auto stack_owner = std::make_unique<Wt::WStackedWidget>();
    m_options_stack = stack_owner.get();
    m_batch_type_menu = options_container->addNew<Wt::WMenu>( m_options_stack );
    m_batch_type_menu->addStyleClass( "LightNavMenu VerticalNavMenu AnaTypeMenu" );
    options_container->addWidget( std::move(stack_owner) );
  }

  {
    auto act_shield_owner = std::make_unique<BatchGuiActShieldAnaWidget>();
    m_act_shield_ana_opts = act_shield_owner.get();
    m_act_shield_ana_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
    m_batch_type_menu->addItem( WString::tr( "bgw-act-shield-ana-opts-label" ),
                                std::move( act_shield_owner ),
                                Wt::ContentLoading::Eager );
  }

  {
    auto peak_fit_owner = std::make_unique<BatchGuiPeakFitWidget>();
    m_peak_fit_opts = peak_fit_owner.get();
    m_peak_fit_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
    m_batch_type_menu->addItem( WString::tr( "bgw-peak-fit-opts-label" ),
                                std::move( peak_fit_owner ),
                                Wt::ContentLoading::Eager );
  }

  {
    auto file_convert_owner = std::make_unique<FileConvertOpts>();
    m_file_convert_opts = file_convert_owner.get();
    m_file_convert_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
    m_batch_type_menu->addItem( WString::tr( "bgw-file-convert-opts-label" ),
                                std::move( file_convert_owner ),
                                Wt::ContentLoading::Eager );
  }
  

  m_batch_type_menu->select( 0 );
  m_batch_type_menu->itemSelected().connect( this, &BatchGuiWidget::updateCanDoAnalysis );

  m_input_files_container = addNew<WGroupBox>( WString::tr( "bgw-input-files-label" ) );
  m_input_files_container->addStyleClass( "InputFilesContainer" );

  m_input_files_container->doJavaScript( "BatchInputDropUploadSetup(" + m_input_files_container->jsRef() +
                                         ", "
                                         " '" +
                                         interspec->fileManager()->batchDragNDrop()->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_input_files_container->id() + "']);" );

  m_output_dir = addNew<DirectorySelector>();
  m_output_dir->setLabelTxt( WString::tr( "bgw-output-dir-label" ) );
  m_output_dir->pathChanged().connect( this, &BatchGuiWidget::updateCanDoAnalysis );

  m_uploadResource->fileDrop().connect( this, &BatchGuiWidget::handleFileDrop );

  m_input_status_error = addNew<WText>();
  m_input_status_error->addStyleClass( "ReasonCantAnalyzeMsg" );
  m_input_status_error->hide();

  
  //const vector<tuple<string,string,bool>> spooled_files = m_uploadResource->takeSpooledFiles();

  // We will load the initial spectrum files, after giving the widget a second to fully load.
  //  I'm not quite sure why, but without doing this, sometimes we can get a JS exception,
  //  maybe because the JS is somehow getting out of order?
  //addInputFiles( spooled_files );

  //std::function<void()> load_files
  //              = wApp->bind( boost::bind( &BatchGuiWidget::addInputFiles, this, spooled_files ) );
  //std::function<void()> worker = [load_files](){
  //  load_files();
  //  wApp->triggerUpdate();
  //};
  
  // Fallback function to clean the files up, incase this session is no longer alive
  //  BUT note that there is a path where if this widget is deleted, before the worker is called,
  //  then the files wont be cleaned up any way.
  //std::function<void()> fall_back = [spooled_files](){
  //  for( const tuple<string, string, bool> &file : spooled_files )
  //  {
  //    const string &path_to_file = std::get<1>( file );
  //    const bool should_delete = std::get<2>( file );
  //    if( should_delete )
  //      SpecUtils::remove_file( path_to_file );
  //  }
  //};//fall_back
  
  //WServer::instance()->schedule( std::chrono::milliseconds(1), wApp->sessionId(), worker, fall_back );

  handleFileDrop( "", "" );

  wApp->triggerUpdate();
}// BatchGuiWidget constructor

BatchGuiWidget::~BatchGuiWidget()
{
  wApp->doJavaScript( "if(window._IS)window._IS.BlockFileDrops=null;" );
  // wApp->doJavaScript( "if(window._IS)window._IS.BatchUploadOnly=null;" );
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_input_files_container->id() + "']);" );
}//~BatchGuiWidget()

Wt::Signal<bool> &BatchGuiWidget::canDoAnalysis()
{
  return m_canDoAnalysis;
}

void BatchGuiWidget::handleFileDrop( const std::string &, const std::string & )
{
  const vector<tuple<string, string, bool>> dropped_files = m_uploadResource->takeSpooledFiles();
  //addInputFiles( dropped_files );
  std::function<void()> worker = [this, dropped_files](){ addInputFiles( dropped_files ); };
  WServer::instance()->schedule( std::chrono::milliseconds(25), wApp->sessionId(), worker );
  wApp->triggerUpdate();
}

void BatchGuiWidget::addInputFiles( const std::vector<std::tuple<std::string, std::string, bool>> &files )
{
  int num_initial_files = m_input_files_container->count();

  for( const tuple<string, string, bool> &file : files )
  {
    const string &display_name = std::get<0>( file );
    const string &path_to_file = std::get<1>( file );
    const bool should_delete = std::get<2>( file );

    const auto show_preview = (num_initial_files < sm_max_spec_file_previews)
                                ? BatchGuiInputSpectrumFile::ShowPreviewOption::Show
                                : BatchGuiInputSpectrumFile::ShowPreviewOption::DontShow;

    BatchGuiInputSpectrumFile *input_file =
      m_input_files_container->addNew<BatchGuiInputSpectrumFile>( display_name, path_to_file, should_delete, show_preview );
    input_file->preview_created_signal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
    input_file->remove_self_request().connect(
      [this]( BatchGuiInputSpectrumFile *f ){ handle_remove_input_file( f ); } );


    num_initial_files += 1;
  }

  wApp->triggerUpdate();
}// void BatchGuiWidget::addInputFiles()

void BatchGuiWidget::handle_remove_input_file( BatchGuiInputSpectrumFile *input )
{
  // removeWidget returns unique_ptr which goes out of scope here, destroying the widget
  m_input_files_container->removeWidget( input );

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

  const pair<bool,WString> ana_status = batch_ana_widget ? batch_ana_widget->canDoAnalysis() : make_pair(false, WString());

  bool can_do_analysis = ana_status.first;
  WString error_msg = ana_status.second;
  if( can_do_analysis && ( num_input_files == 0 ) )
  {
    can_do_analysis = false;
    error_msg = WString::tr("bgw-no-ana-no-input-files");
  }

  if( can_do_analysis && !m_output_dir->isPathValid() )
  {
    can_do_analysis = false;
    error_msg = WString::tr("bgw-no-ana-invalid-output-path");
  }

  m_input_status_error->setHidden( can_do_analysis );
  if( error_msg != m_input_status_error->text() )
    m_input_status_error->setText( error_msg );

  if( (can_do_analysis != m_can_do_analysis) )
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
