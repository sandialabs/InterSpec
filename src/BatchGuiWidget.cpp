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
#include <Wt/WMenu>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WGroupBox>
#include <Wt/WSvgImage>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WFileDropWidget>
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
#include "InterSpec/DirectorySelector.h"
#include "InterSpec/PeakSearchGuiUtils.h"
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
  class InputFileWidget : public Wt::WContainerWidget
  {
    protected:
      const std::string m_filename;
      const std::string m_display_name;
      const bool m_should_cleanup;

      Wt::WContainerWidget *m_preview_container;
      
      bool m_preview_created;
      std::shared_ptr<SpecMeas> m_spec_meas;

      Wt::Signal<bool> m_preview_created_signal;

  public:
    InputFileWidget( const std::string display_name, 
                     const std::string path_to_file, 
                     const bool should_cleanup, 
                     Wt::WContainerWidget *parent )
      : Wt::WContainerWidget( parent ),
        m_filename( path_to_file ),
        m_display_name( display_name ),
        m_should_cleanup( should_cleanup ),
        m_preview_container( nullptr ),
        m_preview_created_signal( this )
    {
      addStyleClass( "InputFileWidget" );

      m_preview_container = new Wt::WContainerWidget( this );
      m_preview_container->addStyleClass( "InputFilePreview" );

      WText *filename_text = new WText( display_name, this );
      filename_text->addStyleClass( "FilenameText" );
      filename_text->setToolTip( "Full path of disk file: '" + path_to_file + "'" );

      WContainerWidget *close_icon = new WContainerWidget( this );
      close_icon->addStyleClass( "closeicon-wtdefault" );
      close_icon->clicked().connect( boost::bind( &InputFileWidget::removeSelf, this ) );

      const string sessionid = wApp->sessionId();
      std::shared_ptr<int> status_ptr = make_shared<int>( 0 );
      std::shared_ptr<SpecMeas> spec_meas = make_shared<SpecMeas>();
      boost::function<void(void)> updateGuiCallback = wApp->bind( boost::bind( &InputFileWidget::set_spectrum, this, spec_meas, status_ptr ) );

      WServer::instance()->ioService().boost::asio::io_service::post( [updateGuiCallback, spec_meas, status_ptr, display_name, path_to_file, sessionid](){

        const bool success = spec_meas->load_file( path_to_file, SpecUtils::ParserType::Auto, display_name );
        if( success )
        {
          *status_ptr = 1;
          spec_meas->set_filename( SpecUtils::filename(display_name) );
        }else
        {
          *status_ptr = 2;
        }

        WServer::instance()->post( sessionid, updateGuiCallback );
      } );
    }//InputFileWidget constructor


    void set_spectrum( std::shared_ptr<SpecMeas> spec_meas, std::shared_ptr<int> status_ptr )
    {
      WText *preview = new WText(  m_preview_container );

      if( !status_ptr || !spec_meas || (spec_meas->num_measurements() == 0) || ((*status_ptr) != 1) )
      {
        preview->setText( WString::tr("bgw-not-spec-preview") );
        preview->setStyleClass( "NotSpectrumFile" );
      }else
      {
        shared_ptr<const deque<shared_ptr<const PeakDef> > > peaks = nullptr;
        vector<shared_ptr<const ReferenceLineInfo>> reflines;

        const bool compact = true;
        const int width_px = 120, height_px = 70;
        const shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();

        //Find most likely foreground.
        shared_ptr<const SpecUtils::Measurement> preview_meas;
        preview_meas = spec_meas->measurement( size_t(0) );

        shared_ptr<WSvgImage> svg = PeakSearchGuiUtils::renderChartToSvg( preview_meas, peaks, reflines,
                                              0.0, 0.0, width_px, height_px, theme, compact );

        if( svg )
        {
          stringstream strm;
          svg->write( strm );
          preview->setText( strm.str() );
        }else
        {
          preview->setText( WString::tr("bgw-preview-error") );
          preview->setStyleClass( "PreviewError" );
        }

        m_spec_meas = spec_meas;
        m_preview_created = true;
      }//if( !status_ptr || !spec_meas || (spec_meas->num_measurements() == 0) || ((*status_ptr) != 1) )

      m_preview_created_signal.emit( !!m_spec_meas );

      wApp->triggerUpdate();
    }//void set_spectrum()
  

    Wt::Signal<bool> &preview_created_signal()
    {
      return m_preview_created_signal;
    }
  
  
    ~InputFileWidget()
    {
      if( m_should_cleanup )
      {
        const bool success = SpecUtils::remove_file( m_filename );
        if( !success )
          cerr << "BatchGuiDialog::InputFileWidget: Warning, could not delete file '" << m_filename << "'" << endl;
      }
    }//~InputFileWidget()


    std::shared_ptr<SpecMeas> spec_meas() const
    {
      return m_spec_meas;
    }

    void removeSelf()
    {
      delete this;
    }//void removeSelf()

    const std::string &display_name() const
    {
      return m_display_name;
    }

    const std::string &path_to_file() const
    {
      return m_filename;
    }
  };//class InputFileWidget
}// anonymous namespace


/** Base class for all batch analysis widgets. */
class BatchGuiAnaWidget : public Wt::WContainerWidget
{
public:
  BatchGuiAnaWidget( Wt::WContainerWidget *parent )
    : Wt::WContainerWidget( parent )
  {
    addStyleClass( "BatchGuiAnaWidget" );
  }

  virtual void performAnalysis( const vector<tuple<string,string,std::shared_ptr<const SpecMeas>>> &input_files, const string &output_dir ) = 0;
};//BatchGuiAnaWidget


class BatchGuiPeakFitWidget : public BatchGuiAnaWidget
{
protected:
  Wt::WContainerWidget *m_peak_container;

public:
  BatchGuiPeakFitWidget( Wt::WContainerWidget *parent = nullptr )
    : BatchGuiAnaWidget( parent )
  {
    addStyleClass( "BatchGuiPeakFitWidget" );

    new WText( "PeakFitAnaOpts", this );


/*
 struct InterSpec_API BatchPeakFitOptions
  {
    bool to_stdout;
    bool refit_energy_cal;
    bool use_exemplar_energy_cal;
    bool write_n42_with_results;
    bool show_nonfit_peaks;
    bool overwrite_output_files;
    bool create_csv_output;
    bool create_json_output;
    std::string background_subtract_file;
    std::set<int> background_subtract_samples;
    bool use_existing_background_peaks;
    bool use_exemplar_energy_cal_for_background;
    double peak_stat_threshold;
    double peak_hypothesis_threshold;
    std::string template_include_dir;
    std::vector<std::string> report_templates;
    std::vector<std::string> summary_report_templates;
  };//struct BatchPeakFitOptions

  void fit_peaks_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                        const ::BatchPeak::BatchPeakFitOptions &options )
or 
use use:
const BatchPeak::BatchPeakFitResult fit_results
                 = fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, nullptr, {}, options );
But duplicate a lot of the logging code.

Or maybe, and probably a better idea, is to over-ride both these functions to take in measurements in memory, as well as thier names; this way we can probably keep things straight in the reporting by using the real names...
*/

  }

  virtual void performAnalysis( const vector<tuple<string,string,std::shared_ptr<const SpecMeas>>> &input_files, const string &output_dir ) override
  {

  }

};//BatchGuiPeakFitWidget

class BatchGuiActShieldAnaWidget : public BatchGuiPeakFitWidget
{
  protected:
    Wt::WContainerWidget *m_act_shield_container;

public:
  BatchGuiActShieldAnaWidget( Wt::WContainerWidget *parent = nullptr )
    : BatchGuiPeakFitWidget( parent ),
    m_act_shield_container( nullptr )
  {
    addStyleClass( "BatchGuiActShieldAnaWidget" );

    m_act_shield_container = new Wt::WContainerWidget( this );
    m_act_shield_container->addStyleClass( "ActShieldOptionsContainer" );

    new WText( "BatchGuiActShieldAnaWidget", m_act_shield_container );

    /*

  struct InterSpec_API BatchActivityFitOptions
    : public BatchPeak::BatchPeakFitOptions
  {
    bool use_bq = false;
    std::shared_ptr<DetectorPeakResponse> drf_override;
    boost::optional<double> distance_override;
    bool hard_background_sub;
  };//struct BatchActivityFitOptions


We need to put in a "display_name" field for the below function call, so reporting will be more readable - or should just pass in a SpecMeas...

BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          const BatchActivityFitOptions &options );
  */

  }

  virtual void performAnalysis( const vector<tuple<string,string,std::shared_ptr<const SpecMeas>>> &input_files, const string &output_dir ) override
  {

  }
};//BatchGuiActShieldAnaWidget





/*
class FileConvertOpts : public Wt::WContainerWidget
{
public:
  FileConvertOpts( Wt::WContainerWidget *parent )
    : Wt::WContainerWidget( parent )
  {
    addStyleClass( "FileConvertOpts" );

    new WText( "FileConvertOpts", this );
  }
};//FileConvertOpts

class RelActAutoOpts : public Wt::WContainerWidget
{
public:
  RelActAutoOpts( Wt::WContainerWidget *parent )
    : Wt::WContainerWidget( parent )
  {
    addStyleClass( "RelActAutoOpts" );

    new WText( "RelActAutoOpts", this );
  }
};//RelActAutoOpts
*/

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

  m_processBtn = addButton( WString::tr("bgw-analyze-button") );
  m_processBtn->clicked().connect( m_widget, &BatchGuiWidget::performAnalysis );
  m_processBtn->disable();
  m_widget->canDoAnalysis().connect( boost::bind( &WPushButton::setEnabled, m_processBtn, boost::placeholders::_1 ) );

  addButton( WString::tr("Close") );
  
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
  wApp->useStyleSheet( "InterSpec_resources/BatchGuiWidget.css" );
  InterSpec::instance()->useMessageResourceBundle( "BatchGuiWidget" );
  
  addStyleClass( "BatchGuiWidget" );
  
  //doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', true);" );
  doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', true);" );

  Wt::WGroupBox *options_container = new WGroupBox( WString::tr("bgw-type-select-label"),  this );
  options_container->addStyleClass( "TypeSelectContainer" );
  
  m_options_stack = new Wt::WStackedWidget();
  m_batch_type_menu = new Wt::WMenu( m_options_stack, options_container );
  m_batch_type_menu->addStyleClass( "LightNavMenu VerticalNavMenu" );
  options_container->addWidget( m_options_stack );
  
  m_act_shield_ana_opts = new BatchGuiActShieldAnaWidget();
  WMenuItem *item = new WMenuItem( WString::tr("bgw-act-shield-ana-opts-label"), m_act_shield_ana_opts, WMenuItem::LoadPolicy::PreLoading );
  m_batch_type_menu->addItem( item );

  m_peak_fit_opts = new BatchGuiPeakFitWidget();
  item = new WMenuItem( WString::tr("bgw-peak-fit-opts-label"), m_peak_fit_opts, WMenuItem::LoadPolicy::PreLoading );
  m_batch_type_menu->addItem( item );

  m_batch_type_menu->select( 0 );

  //WFileDropWidget *file_drop_widget = new WFileDropWidget( this );
  //file_drop_widget->setWidth( WLength( 100, WLength::Percentage ) );
  //file_drop_widget->drop().connect( std::bind([](){
  //  cerr << "Got droped file - not uploaded yet." << endl;
  //}) ); //The signal triggers if one or more files are dropped.
  //Signal<File*>& newUpload() //The signal triggers when the upload of a file is about to begin - e.g., could create widget to show preview - or could do when drop signal is emmittted
  //Signal<File*>& uploaded() //The signal is triggered if any file finished uploading. - could generate file and steal the spooled file
  //Signal<File*>& uploadFailed() //The signal triggers when an upload failed.
  //Signal< File*, ::uint64_t >& tooLarge() //The signal triggers when a file is too large for upload.
  //m_input_files_container = new WGroupBox( WString::tr("bgw-input-files-label"),  file_drop_widget );
  // Should probably put WFileDropWidget inside `m_input_files_container`, and just create another pointer that either points
  
  m_input_files_container = new WGroupBox( WString::tr("bgw-input-files-label"),  this );
  m_input_files_container->addStyleClass( "InputFilesContainer" );

  m_output_dir = new DirectorySelector( this );
  m_output_dir->setLabelTxt( WString::tr("bgw-output-dir-label") );
  m_output_dir->pathChanged().connect( this, &BatchGuiWidget::updateCanDoAnalysis );

  m_uploadResource->fileDrop().connect( this, &BatchGuiWidget::handleFileDrop );

  addInputFiles( m_uploadResource->takeSpooledFiles() );
}//BatchGuiWidget constructor


BatchGuiWidget::~BatchGuiWidget()
{
  //wApp->doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', null);" );
  wApp->doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', null);" );
  
  m_uploadResource->clearSpooledFiles();
}//~BatchGuiWidget()


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
    const string &display_name = std::get<0>(file);
    const string &path_to_file = std::get<1>(file);
    const bool should_delete = std::get<2>(file);

    InputFileWidget *input_file = new InputFileWidget( display_name, path_to_file, should_delete, m_input_files_container );
    input_file->preview_created_signal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
  }
}//void BatchGuiWidget::addInputFiles()

void BatchGuiWidget::updateCanDoAnalysis()
{
  size_t num_input_files = 0;

  for( Wt::WWidget *child : m_input_files_container->children() )
  {
    InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
    num_input_files += (input_file && input_file->spec_meas());
  }

  bool can_do_analysis = (num_input_files > 0) && m_output_dir->isPathValid();
  if( can_do_analysis != m_can_do_analysis )  
  {
    m_can_do_analysis = can_do_analysis;
    m_canDoAnalysis.emit( can_do_analysis );
  }
}//void updateCanDoAnalysis()


void BatchGuiWidget::performAnalysis()
{
  vector<tuple<string,string,std::shared_ptr<const SpecMeas>>> input_files;
  const string output_dir = m_output_dir->path();

  for( Wt::WWidget *child : m_input_files_container->children() )
  {
    InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
    if( input_file )
    {
      const string display_name = input_file->display_name();
      const string path_to_file = input_file->path_to_file();
      const shared_ptr<const SpecMeas> spec_meas = input_file->spec_meas();

      input_files.push_back( make_tuple( display_name, path_to_file, spec_meas ) );
    }//if( input_file )
  }//for( Wt::WWidget *child : m_input_files_container->children() )

  BatchGuiAnaWidget *batch_ana_widget = dynamic_cast<BatchGuiAnaWidget *>( m_options_stack->currentWidget() );
  assert( batch_ana_widget );
  if( batch_ana_widget )
  {
    batch_ana_widget->performAnalysis( input_files, output_dir );
  }else
  {
    cerr << "BatchGuiWidget::performAnalysis: Warning, no batch analysis widget selected" << endl;
  }
}//void BatchGuiWidget::performAnalysis()
