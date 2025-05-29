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

#include <Wt/WText>
#include <Wt/WMenu>
#include <Wt/Utils>
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
#include <Wt/WResource>
#include <Wt/Http/Request>
#include <Wt/Http/Response>

#include <nlohmann/json.hpp>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"


#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/BatchGuiWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DirectorySelector.h"
#include "InterSpec/NativeFloatSpinBox.h"
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

namespace 
{
  // Create a custom resource to serve the HTML content
  class BatchReportResource : public Wt::WResource
  {
    private:
    const std::string m_html_content;
          
    public:
    BatchReportResource( std::string &&html_content, WObject* parent = nullptr )
      : Wt::WResource(parent), m_html_content( std::move(html_content) )
    {
      setDispositionType( DispositionType::Inline );
    }
          
    virtual ~BatchReportResource()
    {
      beingDeleted();
    }
          
    virtual void handleRequest(const Wt::Http::Request& request, Wt::Http::Response& response) override
    {
      response.setMimeType("text/html; charset=utf-8");
      response.addHeader("Cache-Control", "no-cache, no-store, must-revalidate");
      response.addHeader("Pragma", "no-cache");
      response.addHeader("Expires", "0");
            
      response.out() << m_html_content;
    }
  };//class BatchReportResource


class MiscInputFileWidget : public Wt::WContainerWidget
{
protected:
  const std::string m_filename;
  const std::string m_display_name;
  const bool m_should_cleanup;

  Wt::WText *m_txt;

  Wt::Signal<MiscInputFileWidget *> m_remove_self_request_signal;

public:
  MiscInputFileWidget( const std::string display_name,
                  const std::string path_to_file,
                  const bool should_cleanup,
                  Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
  m_filename( path_to_file ),
  m_display_name( display_name ),
  m_should_cleanup( should_cleanup ),
  m_txt( nullptr )
  {
    addStyleClass( "MiscInputFileWidget" );

    m_txt = new WText( display_name, this );
    m_txt->addStyleClass( "FilenameText" );
    m_txt->setToolTip( "Full path of disk file: '" + path_to_file + "'" );

    WContainerWidget *close_icon = new WContainerWidget( this );
    close_icon->addStyleClass( "closeicon-wtdefault" );
    close_icon->clicked().connect( boost::bind( &MiscInputFileWidget::requestRemoveSelf, this ) );
    close_icon->clicked().preventPropagation();
  }//MiscInputFileWidget constructor

  ~MiscInputFileWidget()
  {
    if( m_should_cleanup )
    {
      const bool success = SpecUtils::remove_file( m_filename );
      if( !success )
        cerr << "BatchGuiDialog::MiscInputFileWidget: Warning, could not delete file '" << m_filename << "'" << endl;
    }
  }//~MiscInputFileWidget()


  void requestRemoveSelf()
  {
    m_remove_self_request_signal.emit(this);
  }//void requestRemoveSelf()

  const std::string &display_name() const
  {
    return m_display_name;
  }

  const std::string &path_to_file() const
  {
    return m_filename;
  }

  Wt::Signal<MiscInputFileWidget *> &remove_self_request()
  {
    return m_remove_self_request_signal;
  }
};//class MiscInputFileWidget
}//namespace

/** Represents a spectrum file, with a little thumbnail preview.

 Not in a anonymous namespace because we have to foreward declare it in the header :(
 */
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
  Wt::Signal<InputFileWidget *> m_remove_self_request_signal;

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
    close_icon->clicked().connect( boost::bind( &InputFileWidget::requestRemoveSelf, this ) );
    close_icon->clicked().preventPropagation();

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
    if( !status_ptr || !spec_meas || (spec_meas->num_measurements() == 0) || ((*status_ptr) != 1) )
    {
      WText *preview = new WText( m_preview_container );
      preview->setText( WString::tr("bgw-not-spec-preview") );
      preview->setStyleClass( "NotSpectrumFile" );
      return;
    }else
    {
      shared_ptr<const deque<shared_ptr<const PeakDef> > > peaks = nullptr;
      vector<shared_ptr<const ReferenceLineInfo>> reflines;

      const bool compact = true;
      const int width_px = 120, height_px = 70;
      const shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();

      //Find most likely foreground.  If passthrough, sum everything
      shared_ptr<const SpecUtils::Measurement> preview_meas;

      const std::set<int> &samples = spec_meas->sample_numbers();
      const vector<string> &det_names = spec_meas->detector_names();
      if( spec_meas->passthrough() )
      {
        try
        {
          preview_meas = spec_meas->sum_measurements(samples, det_names, nullptr );
        }catch( std::exception & )
        {
          WText *preview = new WText( m_preview_container );
          preview->setText( WString::tr("bgw-passthrough-sum-error") );
          return;
        }
      }else
      {
        // If not first passthrough, we'll take the first marked Foreground, Unknown,
        //  Background, sample number, in that order
        shared_ptr<const SpecUtils::Measurement> first_back, first_unk, first_fore;

        for( int sample_num : spec_meas->sample_numbers() )
        {
          shared_ptr<const SpecUtils::Measurement> sample_meas;
          try
          {
            if( det_names.size() == 1 )
              sample_meas = spec_meas->measurement( sample_num, det_names.front() );
            else
              sample_meas = spec_meas->sum_measurements({sample_num}, det_names, nullptr );
          }catch( std::exception &e )
          {
          }//try / catch

          if( !sample_meas )
            continue;

          switch( sample_meas->source_type() )
          {
            case SpecUtils::SourceType::IntrinsicActivity:
            case SpecUtils::SourceType::Calibration:
              break;
            case SpecUtils::SourceType::Background:
              if( !first_back )
                first_back = sample_meas;
              break;
            case SpecUtils::SourceType::Foreground:
              if( !first_fore )
                first_fore = sample_meas;
              break;
            case SpecUtils::SourceType::Unknown:
              if( !first_unk )
                first_unk = sample_meas;
              break;
          }//switch( sample_meas->source_type() )

          if( first_fore )
            break;
        }//for( int sample_num : spec_meas->sample_numbers() )

        if( first_fore )
          preview_meas = first_fore;
        else if( first_unk )
          preview_meas = first_unk;
        else
          preview_meas = first_back;
      }//if( spec_meas->passthrough() ) / else

      if( preview_meas )
      {
        D3SpectrumDisplayDiv *spec = new D3SpectrumDisplayDiv( m_preview_container );
        spec->clicked().preventPropagation();
        spec->setThumbnailMode();
        spec->setData( preview_meas, false );
        spec->resize( WLength(100,WLength::Percentage), WLength(100,WLength::Percentage) );

        //We dont currently need to explicitly set the color theme, as all color theme styling for
        //  D3SpectrumDisplayDiv is globablly applied...
        //InterSpec *interspec = InterSpec::instance();
        //spec->applyColorTheme( interspec->getColorTheme() );
        //interspec->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, spec, boost::placeholders::_1 ) );
      }else
      {
        WText *preview = new WText( m_preview_container );
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

  void requestRemoveSelf()
  {
    m_remove_self_request_signal.emit(this);
  }//void requestRemoveSelf()

  const std::string &display_name() const
  {
    return m_display_name;
  }

  const std::string &path_to_file() const
  {
    return m_filename;
  }

  Wt::Signal<InputFileWidget *> &remove_self_request()
  {
    return m_remove_self_request_signal;
  }
};//class InputFileWidget




/** Base class for all batch analysis widgets. */
class BatchGuiAnaWidget : public Wt::WContainerWidget
{
protected:
  Wt::Signal<bool> m_canDoAnalysis;

public:
  BatchGuiAnaWidget( Wt::WContainerWidget *parent )
    : Wt::WContainerWidget( parent )
  {
    addStyleClass( "BatchGuiAnaWidget" );
  }

  virtual void performAnalysis( const vector<tuple<string,string,std::shared_ptr<const SpecMeas>>> &input_files, const string &output_dir ) = 0;

  virtual bool canDoAnalysis() const = 0;


  Wt::Signal<bool> &canDoAnalysisSignal()
  {
    return m_canDoAnalysis;
  }
};//BatchGuiAnaWidget


class BatchGuiPeakFitWidget : public BatchGuiAnaWidget
{
protected:
  Wt::WContainerWidget *m_peak_options_container = nullptr;

  Wt::WContainerWidget *m_exemplar_input = nullptr;
  Wt::WCheckBox *m_use_current_foreground = nullptr;
  Wt::WContainerWidget *m_exemplar_file_drop = nullptr;
  FileDragUploadResource *m_exemplar_file_resource = nullptr;
  std::shared_ptr<SpecMeas> m_uploaded_exemplar;

  Wt::WContainerWidget *m_background_input = nullptr;
  Wt::WCheckBox *m_use_current_background = nullptr;
  Wt::WCheckBox *m_no_background = nullptr;
  Wt::WContainerWidget *m_background_file_drop = nullptr;
  std::shared_ptr<SpecMeas> m_uploaded_background;
  FileDragUploadResource *m_background_file_resource = nullptr;

  // Option widgets
  Wt::WCheckBox *m_fit_all_peaks = nullptr;
  Wt::WCheckBox *m_refit_energy_cal = nullptr;
  Wt::WCheckBox *m_use_exemplar_energy_cal = nullptr;
  Wt::WCheckBox *m_write_n42_with_results = nullptr;
  Wt::WCheckBox *m_show_nonfit_peaks = nullptr;
  Wt::WCheckBox *m_overwrite_output_files = nullptr;
  Wt::WCheckBox *m_create_csv_output = nullptr;
  Wt::WCheckBox *m_create_json_output = nullptr;
  Wt::WCheckBox *m_use_existing_background_peaks = nullptr;
  Wt::WCheckBox *m_use_exemplar_energy_cal_for_background = nullptr;

  // Threshold options with labels
  Wt::WContainerWidget *m_peak_stat_threshold_container = nullptr;
  Wt::WLabel *m_peak_stat_threshold_label = nullptr;
  NativeFloatSpinBox *m_peak_stat_threshold = nullptr;

  Wt::WContainerWidget *m_peak_hypothesis_threshold_container = nullptr;
  Wt::WLabel *m_peak_hypothesis_threshold_label = nullptr;
  NativeFloatSpinBox *m_peak_hypothesis_threshold = nullptr;

  Wt::WContainerWidget *m_reports_container = nullptr;
  Wt::WCheckBox *m_html_report = nullptr;
  Wt::WCheckBox *m_csv_report = nullptr;
  
  // Per-file custom report variables grouped together
  Wt::WCheckBox *m_per_file_custom_report = nullptr;
  Wt::WContainerWidget *m_per_file_custom_report_container = nullptr;
  FileDragUploadResource *m_per_file_custom_report_resource = nullptr;

  // Summary custom report variables grouped together  
  Wt::WCheckBox *m_summary_custom_report = nullptr;
  Wt::WContainerWidget *m_summary_custom_report_container = nullptr;
  FileDragUploadResource *m_summary_custom_report_resource = nullptr;

public:
  BatchGuiPeakFitWidget( Wt::WContainerWidget *parent = nullptr )
    : BatchGuiAnaWidget( parent )
  {
    addStyleClass( "BatchGuiPeakFitWidget" );

    m_exemplar_input = new WGroupBox( WString::tr("bgw-exemplar-grp-title"),  this );
    m_exemplar_input->addStyleClass( "ExemplarToUseOpt" );
    m_use_current_foreground = new WCheckBox( WString::tr("bgw-exemplar-use-current-fore"), m_exemplar_input );
    m_use_current_foreground->addStyleClass( "CbNoLineBreak" );
    m_exemplar_file_drop = new WContainerWidget( m_exemplar_input );
    m_exemplar_file_drop->addStyleClass( "ExemplarFileDrop" );
    const bool have_fore = !!InterSpec::instance()->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    m_use_current_foreground->setChecked( have_fore );
    m_exemplar_file_drop->setHidden( have_fore );
    m_use_current_foreground->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_use_current_foreground->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_exemplar_file_resource = new FileDragUploadResource( m_exemplar_file_drop );
    m_exemplar_file_drop->doJavaScript("BatchInputDropUploadSetup(" + m_exemplar_file_drop->jsRef() + ", "
                                            " '" + m_exemplar_file_resource->url() + "');" );
    doJavaScript("setupOnDragEnterDom(['" + m_exemplar_file_drop->id() + "']);");
    m_exemplar_file_resource->fileDrop().connect( boost::bind( &BatchGuiPeakFitWidget::exemplarUploaded, this,
                                                 boost::placeholders::_1, boost::placeholders::_2) );




    const bool have_back = !!InterSpec::instance()->displayedHistogram(SpecUtils::SpectrumType::Background);
    m_background_input = new WGroupBox( WString::tr("bgw-back-grp-title"),  this );
    m_background_input->addStyleClass( "ExemplarToUseOpt" );

    m_use_current_background = new WCheckBox( WString::tr("bgw-back-use-current"), m_background_input );
    m_use_current_background->addStyleClass( "CbNoLineBreak" );
    m_use_current_background->setChecked( have_back );
    m_use_current_background->setHidden( !have_back );
    m_use_current_background->checked().connect( this, &BatchGuiPeakFitWidget::useCurrentBackgroundChanged );
    m_use_current_background->unChecked().connect( this, &BatchGuiPeakFitWidget::useCurrentBackgroundChanged );

    m_no_background = new WCheckBox( WString::tr("bgw-back-none"), m_background_input );
    m_no_background->addStyleClass( "CbNoLineBreak" );
    m_no_background->setChecked( !have_back );
    m_no_background->checked().connect( this, &BatchGuiPeakFitWidget::useNoBackgroundChanged );
    m_no_background->unChecked().connect( this, &BatchGuiPeakFitWidget::useNoBackgroundChanged );

    m_background_file_drop = new WContainerWidget( m_background_input );
    m_background_file_drop->addStyleClass( "ExemplarFileDrop" );
    m_background_file_drop->setHidden( true );

    m_background_file_resource = new FileDragUploadResource( m_background_file_drop );
    m_background_file_drop->doJavaScript("BatchInputDropUploadSetup(" + m_background_file_drop->jsRef() + ", "
                                         " '" + m_background_file_resource->url() + "');" );
    doJavaScript("setupOnDragEnterDom(['" + m_background_file_drop->id() + "']);");
    m_background_file_resource->fileDrop().connect( boost::bind( &BatchGuiPeakFitWidget::backgroundUploaded, this,
                                                                boost::placeholders::_1, boost::placeholders::_2) );


    m_peak_options_container = new Wt::WContainerWidget( this );
    m_peak_options_container->addStyleClass( "PeakFitOptionsContainer" );


    WContainerWidget *boolOptions = new Wt::WContainerWidget( m_peak_options_container );
    boolOptions->addStyleClass( "PeakFitBoolOptionsContainer" );

    // Create checkbox options
    m_fit_all_peaks = new Wt::WCheckBox( WString::tr("bgw-fit-all-peaks"), boolOptions );
    m_fit_all_peaks->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_fit_all_peaks, WString::tr("bgw-fit-all-peaks-tt"), true );

    m_refit_energy_cal = new Wt::WCheckBox( WString::tr("bgw-refit-energy-cal"), boolOptions );
    m_refit_energy_cal->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_refit_energy_cal, WString::tr("bgw-refit-energy-cal-tt"), true );

    m_use_exemplar_energy_cal = new Wt::WCheckBox( WString::tr("bgw-use-exemplar-energy-cal"), boolOptions );
    m_use_exemplar_energy_cal->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_use_exemplar_energy_cal, WString::tr("bgw-use-exemplar-energy-cal-tt"), true );

    m_write_n42_with_results = new Wt::WCheckBox( WString::tr("bgw-write-n42-with-results"), boolOptions );
    m_write_n42_with_results->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_write_n42_with_results, WString::tr("bgw-write-n42-with-results-tt"), true );
    m_write_n42_with_results->setChecked( true );

    m_show_nonfit_peaks = new Wt::WCheckBox( WString::tr("bgw-show-nonfit-peaks"), boolOptions );
    m_show_nonfit_peaks->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_show_nonfit_peaks, WString::tr("bgw-show-nonfit-peaks-tt"), true );

    m_overwrite_output_files = new Wt::WCheckBox( WString::tr("bgw-overwrite-output-files"), boolOptions );
    m_overwrite_output_files->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_overwrite_output_files, WString::tr("bgw-overwrite-output-files-tt"), true );

    m_create_csv_output = new Wt::WCheckBox( WString::tr("bgw-create-csv-output"), boolOptions );
    m_create_csv_output->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_create_csv_output, WString::tr("bgw-create-csv-output-tt"), true );
    m_create_csv_output->setChecked( true ); // Default to true as per command line
    
    m_create_json_output = new Wt::WCheckBox( WString::tr("bgw-create-json-output"), boolOptions );
    m_create_json_output->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_create_json_output, WString::tr("bgw-create-json-output-tt"), true );

    m_use_existing_background_peaks = new Wt::WCheckBox( WString::tr("bgw-use-existing-background-peaks"), boolOptions );
    m_use_existing_background_peaks->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_use_existing_background_peaks, WString::tr("bgw-use-existing-background-peaks-tt"), true );

    m_use_exemplar_energy_cal_for_background = new Wt::WCheckBox( WString::tr("bgw-use-exemplar-energy-cal-for-background"), boolOptions );
    m_use_exemplar_energy_cal_for_background->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_use_exemplar_energy_cal_for_background, WString::tr("bgw-use-exemplar-energy-cal-for-background-tt"), true );

    // Create threshold options with labels
    WContainerWidget *float_options = new Wt::WContainerWidget( m_peak_options_container );
    float_options->addStyleClass( "PeakFitFloatOptionsContainer" );

    m_peak_stat_threshold_container = new Wt::WContainerWidget( float_options );
    m_peak_stat_threshold_container->addStyleClass( "ThresholdOptionContainer" );
    
    m_peak_stat_threshold_label = new Wt::WLabel( WString::tr("bgw-peak-stat-threshold-label"), m_peak_stat_threshold_container );
    m_peak_stat_threshold_label->setWordWrap( false );
    m_peak_stat_threshold = new NativeFloatSpinBox( m_peak_stat_threshold_container );
    m_peak_stat_threshold_label->setBuddy( m_peak_stat_threshold );
    m_peak_stat_threshold->setValue( 2.0f ); // Default value from command line
    m_peak_stat_threshold->setRange( 0.0f, 10.0f );
    m_peak_stat_threshold->setSpinnerHidden( true );
    m_peak_stat_threshold->setWidth( 40 );
    HelpSystem::attachToolTipOn( m_peak_stat_threshold, WString::tr("bgw-peak-stat-threshold-tt"), true );

    m_peak_hypothesis_threshold_container = new Wt::WContainerWidget( float_options );
    m_peak_hypothesis_threshold_container->addStyleClass( "ThresholdOptionContainer" );
    
    m_peak_hypothesis_threshold_label = new Wt::WLabel( WString::tr("bgw-peak-hypothesis-threshold-label"), m_peak_hypothesis_threshold_container );
    m_peak_hypothesis_threshold_label->setWordWrap( false );
    m_peak_hypothesis_threshold = new NativeFloatSpinBox( m_peak_hypothesis_threshold_container );
    m_peak_hypothesis_threshold_label->setBuddy( m_peak_hypothesis_threshold );
    m_peak_hypothesis_threshold->setValue( 1.0f ); // Default value from command line
    m_peak_hypothesis_threshold->setRange( 0.0f, 10.0f );
    m_peak_hypothesis_threshold->setSpinnerHidden( true );
    m_peak_hypothesis_threshold->setWidth( 40 );
    HelpSystem::attachToolTipOn( m_peak_hypothesis_threshold, WString::tr("bgw-peak-hypothesis-threshold-tt"), true );


    m_reports_container = new WGroupBox( WString::tr("bgw-reports-grp-title"),  this );
    m_reports_container->addStyleClass( "ReportsContainer" );
    m_html_report = new Wt::WCheckBox( WString::tr("bgw-reports-write-html"), m_reports_container );
    m_html_report->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_html_report, WString::tr("bgw-reports-write-html-tooltip"), true );
    m_html_report->setChecked( true );

    m_csv_report = new Wt::WCheckBox( WString::tr("bgw-reports-write-csv"), m_reports_container );
    m_csv_report->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_csv_report, WString::tr("bgw-reports-write-csv-tooltip"), true );
    m_csv_report->setChecked( false );
    m_csv_report->hide(); //This option should not be used for peaks fitting - use `m_create_csv_output`.

    Wt::WContainerWidget *custom_rpt_per_file_opts = new Wt::WContainerWidget( m_reports_container );
    custom_rpt_per_file_opts->addStyleClass( "CustomReportOptions" );
    m_per_file_custom_report = new Wt::WCheckBox( WString::tr("bgw-reports-write-custom-per-file"), custom_rpt_per_file_opts );
    m_per_file_custom_report->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_per_file_custom_report, WString::tr("bgw-reports-write-custom-per-file-tooltip"), true );
    m_per_file_custom_report->setChecked( false );
    m_per_file_custom_report->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_per_file_custom_report->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_per_file_custom_report_container = new Wt::WContainerWidget( custom_rpt_per_file_opts );
    m_per_file_custom_report_container->addStyleClass( "CustomReportUploader EmptyReportUpload" );

    m_per_file_custom_report_resource = new FileDragUploadResource( m_per_file_custom_report_container );
    m_per_file_custom_report_container->doJavaScript("BatchInputDropUploadSetup(" + m_per_file_custom_report_container->jsRef() + ", "
                                            " '" + m_per_file_custom_report_resource->url() + "');" );

    m_per_file_custom_report_container->doJavaScript("setupOnDragEnterDom(['" + m_per_file_custom_report_container->id() + "']);");
    m_per_file_custom_report_resource->fileDrop().connect( boost::bind( &BatchGuiPeakFitWidget::perFileCustomReportUploaded, this,
                                                                boost::placeholders::_1, boost::placeholders::_2) );
    m_per_file_custom_report_container->hide();

    // Group summary custom report elements together
    Wt::WContainerWidget *custom_rpt_summary_opts = new Wt::WContainerWidget( m_reports_container );
    custom_rpt_summary_opts->addStyleClass( "CustomReportOptions" );
    
    m_summary_custom_report = new Wt::WCheckBox( WString::tr("bgw-reports-write-custom-summary"), custom_rpt_summary_opts );
    m_summary_custom_report->addStyleClass( "CbNoLineBreak" );
    HelpSystem::attachToolTipOn( m_summary_custom_report, WString::tr("bgw-reports-write-custom-summary-tooltip"), true );
    m_summary_custom_report->setChecked( false );
    m_summary_custom_report->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_summary_custom_report->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_summary_custom_report_container = new Wt::WContainerWidget( custom_rpt_summary_opts );
    m_summary_custom_report_container->addStyleClass( "CustomReportUploader EmptyReportUpload" );

    m_summary_custom_report_resource = new FileDragUploadResource( m_summary_custom_report_container );
    m_summary_custom_report_container->doJavaScript("BatchInputDropUploadSetup(" + m_summary_custom_report_container->jsRef() + ", "
                                            " '" + m_summary_custom_report_resource->url() + "');" );

    m_summary_custom_report_container->doJavaScript("setupOnDragEnterDom(['" + m_summary_custom_report_container->id() + "']);");
    m_summary_custom_report_resource->fileDrop().connect( boost::bind( &BatchGuiPeakFitWidget::summaryCustomReportUploaded, this,
                                                                boost::placeholders::_1, boost::placeholders::_2) );
    m_summary_custom_report_container->hide();

    // Connect signals to update analysis capability
    m_fit_all_peaks->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_fit_all_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_refit_energy_cal->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_refit_energy_cal->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_use_exemplar_energy_cal->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_use_exemplar_energy_cal->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_write_n42_with_results->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_write_n42_with_results->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_show_nonfit_peaks->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_show_nonfit_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_overwrite_output_files->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_overwrite_output_files->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_create_csv_output->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_create_csv_output->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_create_json_output->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_create_json_output->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_use_existing_background_peaks->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_use_existing_background_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_use_exemplar_energy_cal_for_background->checked().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_use_exemplar_energy_cal_for_background->changed().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_peak_stat_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_peak_stat_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    m_peak_hypothesis_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );
    m_peak_hypothesis_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::peakFitOptionChanged );

    peakFitOptionChanged();
  }//BatchGuiPeakFitWidget constructor


  ~BatchGuiPeakFitWidget()
  {
    wApp->doJavaScript("removeOnDragEnterDom(['" + m_background_file_drop->id() + "']);");
    wApp->doJavaScript("removeOnDragEnterDom(['" + m_exemplar_file_drop->id() + "']);");
    wApp->doJavaScript("removeOnDragEnterDom(['" + m_per_file_custom_report_container->id() + "']);");
    wApp->doJavaScript("removeOnDragEnterDom(['" + m_summary_custom_report_container->id() + "']);");
  }


  void handleFileUpload( WContainerWidget *dropArea, FileDragUploadResource *resource )
  {
    dropArea->clear();

    const vector<tuple<string,string,bool>> spooled = resource->takeSpooledFiles();
    if( spooled.empty() )
      return;

    const tuple<string,string,bool> first_file = spooled.front();

    const string &display_name = std::get<0>(first_file);
    const string &path_to_file = std::get<1>(first_file);
    const bool should_delete = std::get<2>(first_file);

    InputFileWidget *input = new InputFileWidget( display_name, path_to_file, should_delete, dropArea );
    dropArea->removeStyleClass( "EmptyExemplarUpload" );
    input->remove_self_request().connect( boost::bind( &BatchGuiPeakFitWidget::handle_remove_exemplar_upload, this, boost::placeholders::_1 ) );

    // Cleanup all uploads we wont use - dont expect this to actually happen
    for( size_t i = 1; i < spooled.size(); ++i )
    {
      if( get<2>(spooled[i]) )
      {
        const bool success = SpecUtils::remove_file( std::get<1>(spooled[i]) );
        if( !success )
          cerr << "exemplarUploaded: Warning, could not delete file '" << std::get<1>(spooled[i]) << "'" << endl;
      }
    }

    peakFitOptionChanged();
  }//void exemplarUploaded( std::string, std::string )


  void exemplarUploaded( const std::string &, const std::string & )
  {
    handleFileUpload( m_exemplar_file_drop, m_exemplar_file_resource );
  }


  void backgroundUploaded( const std::string &, const std::string & )
  {
    handleFileUpload( m_background_file_drop, m_background_file_resource );
  }


  void perFileCustomReportUploaded( const std::string &, const std::string & )
  {
    const vector<tuple<string,string,bool>> spooled = m_per_file_custom_report_resource->takeSpooledFiles();

    for( const auto &file : spooled )
    {
      const string &display_name = std::get<0>(file);
      const string &path_to_file = std::get<1>(file);
      const bool should_delete = std::get<2>(file);

      MiscInputFileWidget *input = new MiscInputFileWidget( display_name, path_to_file, should_delete, m_per_file_custom_report_container );
      input->remove_self_request().connect( boost::bind( &BatchGuiPeakFitWidget::handle_remove_per_file_custom_report_upload, this, boost::placeholders::_1 ) );
    }//for( const auto &file : spooled )

    peakFitOptionChanged();
  }//void perFileCustomReportUploaded( const std::string &, const std::string & )


  void summaryCustomReportUploaded( const std::string &, const std::string & )
  {
    const vector<tuple<string,string,bool>> spooled = m_summary_custom_report_resource->takeSpooledFiles();

    for( const auto &file : spooled )
    {
      const string &display_name = std::get<0>(file);
      const string &path_to_file = std::get<1>(file);
      const bool should_delete = std::get<2>(file);

      MiscInputFileWidget *input = new MiscInputFileWidget( display_name, path_to_file, should_delete, m_summary_custom_report_container );
      input->remove_self_request().connect( boost::bind( &BatchGuiPeakFitWidget::handle_remove_summary_custom_report_upload, this, boost::placeholders::_1 ) );
    }//for( const auto &file : spooled )

    peakFitOptionChanged();
  }//void summaryCustomReportUploaded( const std::string &, const std::string & )


  void handle_remove_per_file_custom_report_upload( MiscInputFileWidget *input )
  {
    delete input;
    peakFitOptionChanged();
  }

  void handle_remove_summary_custom_report_upload( MiscInputFileWidget *input )
  {
    delete input;
    peakFitOptionChanged();
  }

  void handle_remove_exemplar_upload( InputFileWidget *input )
  {
    delete input;
    peakFitOptionChanged();
  }

  void useCurrentBackgroundChanged()
  {
    if( m_use_current_background->isChecked() )
      m_no_background->setChecked( false );

    peakFitOptionChanged();
  }//void useCurrentBackgroundChanged()


  void useNoBackgroundChanged()
  {
    if( m_no_background->isChecked() )
      m_use_current_background->setChecked( false );

    peakFitOptionChanged();
  }//void useNoBackgroundChanged()


  void peakFitOptionChanged()
  {
    const bool fit_all_peaks = m_fit_all_peaks->isChecked();
    m_refit_energy_cal->setHidden( fit_all_peaks );
    //m_exemplar_input->setHidden( fit_all_peaks );
    //m_use_exemplar_energy_cal->setHidden( fit_all_peaks );
    m_background_input->setHidden( fit_all_peaks );
    m_show_nonfit_peaks->setHidden( fit_all_peaks );
    m_use_existing_background_peaks->setHidden( fit_all_peaks );
    m_use_exemplar_energy_cal_for_background->setHidden( fit_all_peaks );
    m_peak_stat_threshold_container->setHidden( fit_all_peaks );       //Fit thresholds not currently implemented when fitting all peaks
    m_peak_hypothesis_threshold_container->setHidden( fit_all_peaks ); //Fit thresholds not currently implemented when fitting all peaks

    const bool use_current_fore = m_use_current_foreground->isChecked();
    m_exemplar_file_drop->setHidden( use_current_fore );
    if( !use_current_fore )
    {
      if( !m_exemplar_file_drop->isVisible() )
        m_exemplar_file_drop->show();

      const bool has_drop_class = m_exemplar_file_drop->hasStyleClass( "EmptyExemplarUpload" );
      const size_t num_kids = m_exemplar_file_drop->children().size();
      if( (num_kids > 0) && has_drop_class )
        m_exemplar_file_drop->removeStyleClass( "EmptyExemplarUpload" );
      else if( (num_kids == 0) && !has_drop_class )
        m_exemplar_file_drop->addStyleClass( "EmptyExemplarUpload" );
    }else
    {
      if( m_exemplar_file_drop->isVisible() )
        m_exemplar_file_drop->hide();
    }//if( !use_current_fore ) / else


    const bool have_back = !!InterSpec::instance()->displayedHistogram(SpecUtils::SpectrumType::Background);

    const bool use_current_back = m_use_current_background->isChecked();
    const bool use_no_back = m_no_background->isChecked();
    assert( !use_current_back || (use_current_back != use_no_back) ); //shouldnt both be checked at the same time
    auto hide_back_filedrop = [this](){
      const size_t num_kids = m_background_file_drop->children().size();
      if( num_kids )
        m_background_file_drop->clear();
      if( !m_background_file_drop->isHidden() )
        m_background_file_drop->hide();
      const bool has_drop_class = m_background_file_drop->hasStyleClass( "EmptyExemplarUpload" );
      if( !has_drop_class )
        m_background_file_drop->addStyleClass( "EmptyExemplarUpload" );
    };

    m_use_exemplar_energy_cal_for_background->setHidden( use_no_back || fit_all_peaks );

    if( use_no_back )
    {
      assert( !use_current_back );
      if( use_current_back )
        m_use_current_background->setChecked( false );

      hide_back_filedrop();
    }else if( use_current_back )
    {
      assert( !use_no_back );
      hide_back_filedrop();
    }else
    {
      const size_t num_kids = m_background_file_drop->children().size();
      const bool has_drop_class = m_background_file_drop->hasStyleClass( "EmptyExemplarUpload" );

      if( m_background_file_drop->isHidden() )
        m_background_file_drop->show();
      if( num_kids && has_drop_class )
        m_background_file_drop->removeStyleClass( "EmptyExemplarUpload" );
      if( !num_kids && !has_drop_class )
        m_background_file_drop->addStyleClass( "EmptyExemplarUpload" );
    }


    m_per_file_custom_report_container->setHidden( !m_per_file_custom_report->isChecked() );
    m_summary_custom_report_container->setHidden( !m_summary_custom_report->isChecked() );

    const bool have_per_file_custom_report = !m_per_file_custom_report_container->children().empty();
    const bool per_file_has_empty_class = m_per_file_custom_report_container->hasStyleClass( "EmptyReportUpload" );
    if( have_per_file_custom_report && per_file_has_empty_class )
      m_per_file_custom_report_container->removeStyleClass( "EmptyReportUpload" );
    if( !have_per_file_custom_report && !per_file_has_empty_class )
      m_per_file_custom_report_container->addStyleClass( "EmptyReportUpload" );

    const bool have_summary_custom_report = !m_summary_custom_report_container->children().empty();
    const bool summary_has_empty_class = m_summary_custom_report_container->hasStyleClass( "EmptyReportUpload" );
    if( have_summary_custom_report && summary_has_empty_class )
      m_summary_custom_report_container->removeStyleClass( "EmptyReportUpload" );
    if( !have_summary_custom_report && !summary_has_empty_class )
      m_summary_custom_report_container->addStyleClass( "EmptyReportUpload" );

    // Emit signal to parent that analysis capability may have changed
    m_canDoAnalysis.emit( canDoAnalysis() );
  }//void peakFitOptionChanged()


  tuple<shared_ptr<SpecMeas>,string,set<int>> get_exemplar() const
  {
    string filename;
    set<int> samples;
    shared_ptr<SpecMeas> meas;

    if( m_use_current_foreground->isChecked() )
    {
      InterSpec *interspec = InterSpec::instance();

      // We dont need to save the tool states to the SpecMeas here, but we will anyway, just to be safe if anything
      //  changes in the future
      interspec->saveShieldingSourceModelToForegroundSpecMeas();
    #if( USE_REL_ACT_TOOL )
      interspec->saveRelActManualStateToForegroundSpecMeas();
      interspec->saveRelActAutoStateToForegroundSpecMeas();
    #endif

      std::shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground );
      samples = interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );

      if( foreground )
      {
        filename = foreground->filename();

        shared_ptr<SpecMeas> spec_copy = make_shared<SpecMeas>();
        spec_copy->uniqueCopyContents( *foreground );

        meas = spec_copy;
      }//if( foreground )
    }else
    {
      for( Wt::WWidget *child : m_exemplar_file_drop->children() )
      {
        InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
        assert( input_file );
        if( !input_file )
          continue;

        meas = input_file->spec_meas();
        filename = input_file->path_to_file();
        if( meas )
        {
          shared_ptr<SpecMeas> spec_copy = make_shared<SpecMeas>();
          spec_copy->uniqueCopyContents( *meas );
          meas = spec_copy;
        }

        if( meas )
          break;
      }
    }//if( m_use_current_foreground->isChecked() ) / else


    return {meas,filename,samples};
  }//tuple<shared_ptr<const SpecMeas>,string,set<int>> get_exemplar() const


  tuple<shared_ptr<const SpecMeas>,string,set<int>> get_background() const
  {
    string filename;
    set<int> samples;
    shared_ptr<const SpecMeas> meas;

    if( m_no_background->isChecked() )
    {
      // Nothing to do here
    }else if( m_use_current_background->isChecked() )
    {
      InterSpec *interspec = InterSpec::instance();
      meas = interspec->measurment( SpecUtils::SpectrumType::Background );
      samples = interspec->displayedSamples( SpecUtils::SpectrumType::Background );
      if( meas )
        filename = meas->filename();
    }else
    {
      for( Wt::WWidget *child : m_exemplar_file_drop->children() )
      {
        InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
        assert( input_file );
        if( !input_file )
          continue;

        meas = input_file->spec_meas();
        filename = input_file->path_to_file();
        break;
      }
    }//if( m_use_current_foreground->isChecked() ) / else


    return {meas,filename,samples};
  }//tuple<shared_ptr<const SpecMeas>,string,set<int>> get_background() const


  BatchPeak::BatchPeakFitOptions getPeakFitOptions() const
  {
    BatchPeak::BatchPeakFitOptions answer;

    answer.to_stdout = false;
    answer.fit_all_peaks = m_fit_all_peaks->isChecked();
    if( answer.fit_all_peaks )
    {
      answer.refit_energy_cal = false;
      answer.show_nonfit_peaks = false;
    }else
    {
      answer.refit_energy_cal = m_refit_energy_cal->isChecked();
      answer.show_nonfit_peaks = m_show_nonfit_peaks->isChecked();
      answer.peak_stat_threshold = m_peak_stat_threshold->value();
      answer.peak_hypothesis_threshold = m_peak_hypothesis_threshold->value();
    }

    tuple<shared_ptr<const SpecMeas>,string,set<int>> exemplar_info = get_exemplar();

    if( get<0>(exemplar_info) )
      answer.use_exemplar_energy_cal = m_use_exemplar_energy_cal->isChecked();

    answer.write_n42_with_results = m_write_n42_with_results->isChecked();
    answer.overwrite_output_files = m_overwrite_output_files->isChecked();
    answer.create_csv_output = m_create_csv_output->isChecked();
    answer.create_json_output = m_create_json_output->isChecked();

    if( m_no_background->isChecked() )
    {
      answer.use_existing_background_peaks = false;
      answer.use_exemplar_energy_cal_for_background = false;
    }else
    {
      const tuple<shared_ptr<const SpecMeas>,string,set<int>> background_info = get_background();
      const shared_ptr<const SpecMeas> &background_spec = get<0>(background_info);
      if( background_spec )
      {
        shared_ptr<SpecMeas> spec_copy = make_shared<SpecMeas>();
        spec_copy->uniqueCopyContents( *background_spec );
        answer.cached_background_subtract_spec = spec_copy;
      }
      
      answer.background_subtract_file = get<1>(background_info);
      answer.background_subtract_samples = get<2>(background_info);
      answer.use_existing_background_peaks = m_use_existing_background_peaks->isChecked();
      answer.use_exemplar_energy_cal_for_background = m_use_exemplar_energy_cal_for_background->isChecked();
    }

    answer.create_csv_output = m_create_csv_output->isChecked();
    answer.create_json_output = m_create_json_output->isChecked();


    if( m_html_report->isChecked() )
    {
      answer.report_templates.push_back( "html" );
      answer.summary_report_templates.push_back( "html" );
    }
   
    if( m_csv_report->isVisible() && m_csv_report->isChecked() )
    {
      answer.report_templates.push_back( "csv" );
      answer.summary_report_templates.push_back( "csv" );
    }
   
    if( m_per_file_custom_report->isChecked() )
    {
      for( Wt::WWidget *child : m_per_file_custom_report_container->children() )
      {
        MiscInputFileWidget *input_file = dynamic_cast<MiscInputFileWidget *>( child );
        assert( input_file );
        if( !input_file )
          continue;

        answer.report_templates.push_back( input_file->path_to_file() );
      }
    }
   
    if( m_summary_custom_report->isChecked() )
    {
      for( Wt::WWidget *child : m_summary_custom_report_container->children() )
      {
        MiscInputFileWidget *input_file = dynamic_cast<MiscInputFileWidget *>( child );
        assert( input_file );
        if( !input_file )
          continue;

        answer.summary_report_templates.push_back( input_file->path_to_file() );
      }
    }//if( m_summary_custom_report->isChecked() )
   
    //answer.template_include_dir = ...; //std::string - if we allow setting this, then the custom report templates cant be used.
    //answer.output_dir = ...; //std::string

    return answer;
  }//BatchPeak::BatchPeakFitOptions getOptions() const


  virtual void performAnalysis( const vector<tuple<string,string,shared_ptr<const SpecMeas>>> &input_files, const string &output_dir ) override
  {
    BatchPeak::BatchPeakFitOptions options = getPeakFitOptions();
    options.output_dir = output_dir;

    const tuple<shared_ptr<const SpecMeas>,string,set<int>> exemplar_info = get_exemplar();
    const shared_ptr<const SpecMeas> &exemplar = get<0>(exemplar_info);
    const string &exemplar_filename = get<1>(exemplar_info);
    const set<int> &exemplar_samples = get<2>(exemplar_info);

    const tuple<shared_ptr<const SpecMeas>,string,set<int>> background_info = get_background();

    vector<string> file_names;
    vector<shared_ptr<SpecMeas>> spec_files;
    for( size_t i = 0; i < input_files.size(); ++i )
    {
      const string &display_name = get<0>(input_files[i]);
      const string &path_to_file = get<1>(input_files[i]);
      const shared_ptr<const SpecMeas> &spec_meas = get<2>(input_files[i]);
      assert( spec_meas );
      if( !spec_meas )
        continue;

      file_names.push_back( spec_meas->filename() );

      shared_ptr<SpecMeas> spec_copy = make_shared<SpecMeas>(); //we need it non-const (and we dont want to mess with in-memory files, so we get consistent resutls with multiple calls)
      spec_copy->uniqueCopyContents( *spec_meas );
      spec_files.push_back( spec_copy );
    }//for( size_t i = 0; i < input_files.size(); ++i )

    // We need to do the work off the main GUI thread, so we need to post the work to the server.
    const string sessionid = Wt::WApplication::instance()->sessionId();
    auto error_msg = make_shared<string>();
    auto results = make_shared<BatchPeak::BatchPeakFitSummaryResults>();

    SimpleDialog *waiting_dialog = new SimpleDialog( WString::tr("bgw-performing-work-title"), WString::tr("bgw-performing-work-msg") );
    waiting_dialog->addButton( WString::tr("Close") );
    boost::function<void(void)> close_waiting_dialog = wApp->bind( boost::bind( &SimpleDialog::done, waiting_dialog, Wt::WDialog::DialogCode::Accepted ) );


    std::function<void(void)> show_error_dialog = [error_msg, close_waiting_dialog](){
      close_waiting_dialog();
      SimpleDialog *dialog = new SimpleDialog( WString::tr("bgw-error-analysis-title"), WString::tr("bgw-error-analysis-msg").arg(*error_msg) );
      dialog->addStyleClass( "BatchAnalysisErrorDialog" );
      dialog->addButton( WString::tr("Okay") );
      wApp->triggerUpdate();
    };

    std::function<void(void)> update_gui_fcn = [results, options, close_waiting_dialog](){
      close_waiting_dialog();

      SimpleDialog *dialog = new SimpleDialog( WString::tr("bgw-analysis-summary-title") );
      dialog->addStyleClass( "BatchAnalysisResultDialog" );
      
      try
      {
        nlohmann::json summary_json = nlohmann::json::parse(results->summary_json);
        inja::Environment env = BatchInfoLog::get_default_inja_env( options );
        string rpt = BatchInfoLog::render_template("html", env,
                         BatchInfoLog::TemplateRenderType::PeakFitSummary, options, summary_json );
        
        // Create the resource and get its URL
        BatchReportResource *report_resource = new BatchReportResource( std::move(rpt), dialog);
        const string resource_url = report_resource->url();
        
        const string contents = "<iframe "
                     "style=\"width: 100%; height: 100%; border: none;\" "
                     "sandbox=\"allow-scripts\" "
                     "src=\"" + resource_url + "\" "
                     "onerror=\"console.error('Iframe failed to load')\"> "
                     "</iframe>";
        
        dialog->resize( WLength(80.0,WLength::Percentage), WLength(95.0,WLength::Percentage) );
        
        WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
        result_text->addStyleClass( "BatchReportIFrameHolder" );
        result_text->setWidth( WLength(100.0,WLength::Percentage) );
        result_text->setHeight( WLength(100.0,WLength::Percentage) );
      }catch( nlohmann::json::parse_error &e )
      {
        const string contents = "<p>" + WString::tr("bgw-json-parse-error").arg(string(e.what())).toUTF8() + "</p>";
        WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
        result_text->addStyleClass( "BatchReportJsonError" );
      }catch( inja::InjaError &e )
      {
        const string contents = "<p>" + WString::tr("bgw-inja-error").arg(e.message).arg(std::to_string(e.location.line)).arg(std::to_string(e.location.column)).toUTF8() + "</p>";
        
        WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
        result_text->addStyleClass( "BatchReportInjaError" );
      }catch( std::exception &e )
      {
        const string contents = "<p>" + WString::tr("bgw-misc-error").arg(string(e.what())).toUTF8() + "</p>";
        WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
        result_text->addStyleClass( "BatchReportMiscError" );
      }
      
      dialog->addButton( WString::tr("Okay") );
      
      
      if( !results->warnings.empty() )
      {
        SimpleDialog *warnings_dialog = new SimpleDialog( WString::tr("bgw-warning-title") );
        warnings_dialog->addStyleClass( "BatchAnalysisWarningDialog" );
        
        WContainerWidget *contents = warnings_dialog->contents();
        contents->setList( true, false );
        for( const string &warn_msg : results->warnings )
        {
          WContainerWidget *item = new WContainerWidget( contents );
          new WText( warn_msg, item );
        }
        
        warnings_dialog->addButton( WString::tr("Okay") );
      }//if( !results.warnings.empty() )

      wApp->triggerUpdate();
    };//update_gui_fcn lambda
    

    std::function<void(void)> do_work_fcn = [exemplar_filename, exemplar, exemplar_samples, file_names, spec_files, options, results, error_msg, update_gui_fcn, show_error_dialog, sessionid](){
      try
      {
        BatchPeak::fit_peaks_in_files( exemplar_filename, exemplar, exemplar_samples,
                                    file_names, spec_files, options, results.get() );
        
        WServer::instance()->post( sessionid, update_gui_fcn );
      }catch( std::exception &e )
      {
        *error_msg = e.what();
        WServer::instance()->post( sessionid, show_error_dialog );
      }
    };//do_work_fcn lambda

    // Do the peak fitting work in a non-gui thread.
    WServer::instance()->ioService().boost::asio::io_service::post( do_work_fcn );
  }//performAnalysis(...)


  
  virtual bool canDoAnalysis() const override
  {
    bool have_exemplar = false;
    if( m_use_current_foreground->isChecked() )
    {
      have_exemplar = !!(InterSpec::instance()->displayedHistogram(SpecUtils::SpectrumType::Foreground));
    }else
    {
      for( Wt::WWidget *child : m_exemplar_file_drop->children() )
      {
        InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
        have_exemplar = (input_file && input_file->spec_meas());
        if( have_exemplar )
          break;
      }
    }

    if( !have_exemplar )
      return false;

    bool back_status_okay = m_no_background->isChecked();
    if( !back_status_okay && m_use_current_background->isChecked() )
    {
      back_status_okay = !!(InterSpec::instance()->displayedHistogram(SpecUtils::SpectrumType::Background));
    }

    if( !back_status_okay && !m_use_current_background->isChecked() && !m_no_background->isChecked() )
    {
      for( Wt::WWidget *child : m_background_file_drop->children() )
      {
        InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
        have_exemplar = (input_file && input_file->spec_meas());
        if( have_exemplar )
          break;
      }
    }

    if( !back_status_okay )
      return false;

    return true;
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

    m_fit_all_peaks->setChecked( false );
    m_fit_all_peaks->hide();

    m_act_shield_container = new Wt::WContainerWidget();
    m_act_shield_container->addStyleClass( "ActShieldOptionsContainer" );
    insertWidget( 0, m_act_shield_container );



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

  }//void performAnalysis()

  virtual bool canDoAnalysis() const override
  {
    return true;
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

  m_processBtn = new WPushButton( WString::tr("bgw-analyze-button"), footer() );
  m_processBtn->setStyleClass( "simple-dialog-btn" );
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
  InterSpec *interspec = InterSpec::instance();
  WApplication *app = WApplication::instance();

  app->useStyleSheet( "InterSpec_resources/BatchGuiWidget.css" );
  app->require( "InterSpec_resources/BatchGuiWidget.js" );
  interspec->useMessageResourceBundle( "BatchGuiWidget" );

  addStyleClass( "BatchGuiWidget" );
  
  doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', true);" );
  //doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', true);" );

  Wt::WGroupBox *options_container = new Wt::WGroupBox( WString::tr("bgw-type-select-label"),  this );
  options_container->addStyleClass( "TypeSelectContainer" );
  
  m_options_stack = new Wt::WStackedWidget();
  m_batch_type_menu = new Wt::WMenu( m_options_stack, options_container );
  m_batch_type_menu->addStyleClass( "LightNavMenu VerticalNavMenu" );
  options_container->addWidget( m_options_stack );
  
  m_act_shield_ana_opts = new BatchGuiActShieldAnaWidget();
  WMenuItem *item = new WMenuItem( WString::tr("bgw-act-shield-ana-opts-label"), m_act_shield_ana_opts, WMenuItem::LoadPolicy::PreLoading );
  m_act_shield_ana_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
  m_batch_type_menu->addItem( item );

  m_peak_fit_opts = new BatchGuiPeakFitWidget();
  item = new WMenuItem( WString::tr("bgw-peak-fit-opts-label"), m_peak_fit_opts, WMenuItem::LoadPolicy::PreLoading );
  m_peak_fit_opts->canDoAnalysisSignal().connect( this, &BatchGuiWidget::updateCanDoAnalysis );
  m_batch_type_menu->addItem( item );

  m_batch_type_menu->select( 0 );

  m_input_files_container = new WGroupBox( WString::tr("bgw-input-files-label"),  this );
  m_input_files_container->addStyleClass( "InputFilesContainer" );
  
  m_input_files_container->doJavaScript("BatchInputDropUploadSetup(" + m_input_files_container->jsRef() + ", "
                                        " '" + interspec->fileManager()->batchDragNDrop()->url() + "');" );
  doJavaScript("setupOnDragEnterDom(['" + m_input_files_container->id() + "']);");

  m_output_dir = new DirectorySelector( this );
  m_output_dir->setLabelTxt( WString::tr("bgw-output-dir-label") );
  m_output_dir->pathChanged().connect( this, &BatchGuiWidget::updateCanDoAnalysis );

  m_uploadResource->fileDrop().connect( this, &BatchGuiWidget::handleFileDrop );

  addInputFiles( m_uploadResource->takeSpooledFiles() );
}//BatchGuiWidget constructor


BatchGuiWidget::~BatchGuiWidget()
{
  wApp->doJavaScript( "$('.Wt-domRoot').data('BlockFileDrops', null);" );
  //wApp->doJavaScript( "$('.Wt-domRoot').data('BatchUploadOnly', null);" );
  wApp->doJavaScript("removeOnDragEnterDom(['" + m_input_files_container->id() + "']);");
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
    input_file->remove_self_request().connect( boost::bind( &BatchGuiWidget::handle_remove_input_file, this, boost::placeholders::_1 ) );
  }
}//void BatchGuiWidget::addInputFiles()


void BatchGuiWidget::handle_remove_input_file( InputFileWidget *input )
{
  delete input;

  updateCanDoAnalysis();
}//void handle_remove_input_file( InputFileWidget *input )


void BatchGuiWidget::updateCanDoAnalysis()
{
  size_t num_input_files = 0;

  for( Wt::WWidget *child : m_input_files_container->children() )
  {
    InputFileWidget *input_file = dynamic_cast<InputFileWidget *>( child );
    num_input_files += (input_file && input_file->spec_meas());
  }

  BatchGuiAnaWidget *const batch_ana_widget = dynamic_cast<BatchGuiAnaWidget *>( m_options_stack->currentWidget() );

  const bool can_do_analysis = (num_input_files > 0) 
                                && m_output_dir->isPathValid()
                                && (batch_ana_widget && batch_ana_widget->canDoAnalysis());
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
