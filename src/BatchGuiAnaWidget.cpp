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

#include "InterSpec/BatchGuiAnaWidget.h"

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <functional>
#include <cassert>
#include <chrono>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <Wt/WText>
#include <Wt/WMenu>
#include <Wt/Utils>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WGroupBox>
#include <Wt/WResource>
#include <Wt/WSvgImage>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WStackedWidget>
#include <Wt/WFileDropWidget>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>

#include <nlohmann/json.hpp>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/BatchGuiAnaWidget.h"
#include "InterSpec/DirectorySelector.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/FileDragUploadResource.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/BatchGuiInputFile.h"

#if ( USE_REL_ACT_TOOL )
#include "InterSpec/RelActCalc.h"
#endif

#include "rapidxml/rapidxml.hpp"

using namespace Wt;
using namespace std;

namespace
{
  void right_select_item( WMenu *menu, WMenuItem *item )
  {
    menu->select( item );
    item->triggered().emit( item ); //
  }
  
  // This is a duplicate of `ExportSpecFileTool::maxRecordsInCurrentSaveType(...)`.
  uint16_t max_records_in_save_type( const SpecUtils::SaveSpectrumAsType save_format, std::shared_ptr<const SpecMeas> spec )
  {
    switch( save_format )
    {
        // Spectrum file types that can have many spectra in them
      case SpecUtils::SaveSpectrumAsType::Txt:
      case SpecUtils::SaveSpectrumAsType::Csv:
      case SpecUtils::SaveSpectrumAsType::Pcf:
      case SpecUtils::SaveSpectrumAsType::N42_2006:
      case SpecUtils::SaveSpectrumAsType::N42_2012:
      case SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0:
      case SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2:
        return std::numeric_limits<uint16_t>::max();
        
#if( USE_QR_CODES )
        // Spectrum file types that can have two spectra in them (foreground + background)
      case SpecUtils::SaveSpectrumAsType::Uri:
        return (spec && (spec->num_gamma_channels() < 2075)) ? 2 : 1;
#endif
        
#if( SpecUtils_ENABLE_D3_CHART )
      case SpecUtils::SaveSpectrumAsType::HtmlD3:
        return 2;
#endif
        
        // Spectrum file types that can have a single spectra in them
      case SpecUtils::SaveSpectrumAsType::Chn:
      case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
      case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
      case SpecUtils::SaveSpectrumAsType::SpcAscii:
      case SpecUtils::SaveSpectrumAsType::SpeIaea:
      case SpecUtils::SaveSpectrumAsType::Cnf:
      case SpecUtils::SaveSpectrumAsType::Tka:
        return 1;
        
#if( SpecUtils_INJA_TEMPLATES )
      case SpecUtils::SaveSpectrumAsType::Template:
        assert( 0 );
        break;
#endif
        
      case SpecUtils::SaveSpectrumAsType::NumTypes:
        assert( 0 );
        break;
    }//switch( save_format )
    
    assert( 0 );
    
    return 0;
  }//max_records_in_save_type(...)
  
// Create a custom resource to serve the HTML content
class BatchReportResource : public Wt::WResource
{
private:
  const std::string m_html_content;

public:
  BatchReportResource( std::string &&html_content, WObject *parent = nullptr )
  : Wt::WResource( parent ), m_html_content( std::move( html_content ) )
  {
    setDispositionType( DispositionType::Inline );
  }

  virtual ~BatchReportResource()
  {
    beingDeleted();
  }

  virtual void handleRequest( const Wt::Http::Request &request, Wt::Http::Response &response ) override
  {
    response.setMimeType( "text/html; charset=utf-8" );
    response.addHeader( "Cache-Control", "no-cache, no-store, must-revalidate" );
    response.addHeader( "Pragma", "no-cache" );
    response.addHeader( "Expires", "0" );

    response.out() << m_html_content;
  }
};// class BatchReportResource

}// namespace



BatchGuiAnaWidget::BatchGuiAnaWidget( Wt::WContainerWidget *parent )
  : Wt::WContainerWidget( parent ),
  m_canDoAnalysis( this )
{
  addStyleClass( "BatchGuiAnaWidget" );

  WApplication *app = WApplication::instance();
  InterSpec *interspec = InterSpec::instance();

  app->useStyleSheet( "InterSpec_resources/BatchGuiAnaWidget.css" );
  app->require( "InterSpec_resources/BatchGuiWidget.js" ); //Should already be loaded, but JIC
  interspec->useMessageResourceBundle( "BatchGuiAnaWidget" );
}

Wt::Signal<bool,Wt::WString> &BatchGuiAnaWidget::canDoAnalysisSignal()
{
  return m_canDoAnalysis;
}

BatchGuiPeakFitWidget::BatchGuiPeakFitWidget( Wt::WContainerWidget *parent ) : BatchGuiAnaWidget( parent )
{
  addStyleClass( "BatchGuiPeakFitWidget" );

  m_exemplar_input = new WGroupBox( WString::tr( "bgw-exemplar-grp-title" ), this );
  m_exemplar_input->addStyleClass( "ExemplarToUseOpt" );
  m_use_current_foreground = new WCheckBox( WString::tr( "bgw-exemplar-use-current-fore" ), m_exemplar_input );
  m_use_current_foreground->addStyleClass( "CbNoLineBreak" );
  m_exemplar_file_drop = new WContainerWidget( m_exemplar_input );
  m_exemplar_file_drop->addStyleClass( "ExemplarFileDrop" );
  const bool have_fore = !!InterSpec::instance()->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  m_use_current_foreground->setChecked( have_fore );
  m_exemplar_file_drop->setHidden( have_fore );
  m_use_current_foreground->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_use_current_foreground->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_exemplar_file_resource = new FileDragUploadResource( m_exemplar_file_drop );
  m_exemplar_file_drop->doJavaScript( "BatchInputDropUploadSetup(" + m_exemplar_file_drop->jsRef() +
                                      ", "
                                      " '" +
                                      m_exemplar_file_resource->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_exemplar_file_drop->id() + "']);" );
  m_exemplar_file_resource->fileDrop().connect(
    boost::bind( &BatchGuiPeakFitWidget::exemplarUploaded, this, boost::placeholders::_1, boost::placeholders::_2 ) );


  const bool have_back = !!InterSpec::instance()->displayedHistogram( SpecUtils::SpectrumType::Background );
  m_background_input = new WGroupBox( WString::tr( "bgw-back-grp-title" ), this );
  m_background_input->addStyleClass( "ExemplarToUseOpt" );

  m_use_current_background = new WCheckBox( WString::tr( "bgw-back-use-current" ), m_background_input );
  m_use_current_background->addStyleClass( "CbNoLineBreak" );
  m_use_current_background->setChecked( have_back );
  m_use_current_background->setHidden( !have_back );
  m_use_current_background->checked().connect( this, &BatchGuiPeakFitWidget::useCurrentBackgroundChanged );
  m_use_current_background->unChecked().connect( this, &BatchGuiPeakFitWidget::useCurrentBackgroundChanged );

  m_no_background = new WCheckBox( WString::tr( "bgw-back-none" ), m_background_input );
  m_no_background->addStyleClass( "CbNoLineBreak" );
  m_no_background->setChecked( !have_back );
  m_no_background->checked().connect( this, &BatchGuiPeakFitWidget::useNoBackgroundChanged );
  m_no_background->unChecked().connect( this, &BatchGuiPeakFitWidget::useNoBackgroundChanged );

  m_background_file_drop = new WContainerWidget( m_background_input );
  m_background_file_drop->addStyleClass( "ExemplarFileDrop" );
  m_background_file_drop->setHidden( true );

  m_background_file_resource = new FileDragUploadResource( m_background_file_drop );
  m_background_file_drop->doJavaScript( "BatchInputDropUploadSetup(" + m_background_file_drop->jsRef() +
                                        ", "
                                        " '" +
                                        m_background_file_resource->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_background_file_drop->id() + "']);" );
  m_background_file_resource->fileDrop().connect(
    boost::bind( &BatchGuiPeakFitWidget::backgroundUploaded, this, boost::placeholders::_1, boost::placeholders::_2 ) );


  m_peak_options_container = new Wt::WContainerWidget( this );
  m_peak_options_container->addStyleClass( "PeakFitOptionsContainer" );


  WContainerWidget *boolOptions = new Wt::WContainerWidget( m_peak_options_container );
  boolOptions->addStyleClass( "PeakFitBoolOptionsContainer" );

  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );

  // Create checkbox options
  m_fit_all_peaks = new Wt::WCheckBox( WString::tr( "bgw-fit-all-peaks" ), boolOptions );
  m_fit_all_peaks->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_fit_all_peaks, WString::tr( "bgw-fit-all-peaks-tt" ), showToolTips );

  m_refit_energy_cal = new Wt::WCheckBox( WString::tr( "bgw-refit-energy-cal" ), boolOptions );
  m_refit_energy_cal->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_refit_energy_cal, WString::tr( "bgw-refit-energy-cal-tt" ), showToolTips );

  m_use_exemplar_energy_cal = new Wt::WCheckBox( WString::tr( "bgw-use-exemplar-energy-cal" ), boolOptions );
  m_use_exemplar_energy_cal->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn(
    m_use_exemplar_energy_cal, WString::tr( "bgw-use-exemplar-energy-cal-tt" ), showToolTips );

  m_write_n42_with_results = new Wt::WCheckBox( WString::tr( "bgw-write-n42-with-results" ), boolOptions );
  m_write_n42_with_results->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_write_n42_with_results, WString::tr( "bgw-write-n42-with-results-tt" ), showToolTips );
  m_write_n42_with_results->setChecked( true );

  m_show_nonfit_peaks = new Wt::WCheckBox( WString::tr( "bgw-show-nonfit-peaks" ), boolOptions );
  m_show_nonfit_peaks->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_show_nonfit_peaks, WString::tr( "bgw-show-nonfit-peaks-tt" ), showToolTips );

  m_overwrite_output_files = new Wt::WCheckBox( WString::tr( "bgw-overwrite-output-files" ), boolOptions );
  m_overwrite_output_files->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_overwrite_output_files, WString::tr( "bgw-overwrite-output-files-tt" ), showToolTips );

  m_create_json_output = new Wt::WCheckBox( WString::tr( "bgw-create-json-output" ), boolOptions );
  m_create_json_output->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_create_json_output, WString::tr( "bgw-create-json-output-tt" ), showToolTips );

  m_use_existing_background_peaks =
    new Wt::WCheckBox( WString::tr( "bgw-use-existing-background-peaks" ), boolOptions );
  m_use_existing_background_peaks->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn(
    m_use_existing_background_peaks, WString::tr( "bgw-use-existing-background-peaks-tt" ), showToolTips );

  m_use_exemplar_energy_cal_for_background =
    new Wt::WCheckBox( WString::tr( "bgw-use-exemplar-energy-cal-for-background" ), boolOptions );
  m_use_exemplar_energy_cal_for_background->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_use_exemplar_energy_cal_for_background,
                               WString::tr( "bgw-use-exemplar-energy-cal-for-background-tt" ),
                               showToolTips );

  // Create threshold options with labels
  WContainerWidget *float_options = new Wt::WContainerWidget( m_peak_options_container );
  float_options->addStyleClass( "PeakFitFloatOptionsContainer" );

  m_peak_stat_threshold_container = new Wt::WContainerWidget( float_options );
  m_peak_stat_threshold_container->addStyleClass( "ThresholdOptionContainer" );

  m_peak_stat_threshold_label =
    new Wt::WLabel( WString::tr( "bgw-peak-stat-threshold-label" ), m_peak_stat_threshold_container );
  m_peak_stat_threshold_label->setWordWrap( false );
  m_peak_stat_threshold = new NativeFloatSpinBox( m_peak_stat_threshold_container );
  m_peak_stat_threshold_label->setBuddy( m_peak_stat_threshold );
  m_peak_stat_threshold->setValue( 2.0f );// Default value from command line
  m_peak_stat_threshold->setRange( 0.0f, 10.0f );
  m_peak_stat_threshold->setSpinnerHidden( true );
  m_peak_stat_threshold->setWidth( 40 );
  HelpSystem::attachToolTipOn( m_peak_stat_threshold, WString::tr( "bgw-peak-stat-threshold-tt" ), showToolTips );

  m_peak_hypothesis_threshold_container = new Wt::WContainerWidget( float_options );
  m_peak_hypothesis_threshold_container->addStyleClass( "ThresholdOptionContainer" );

  m_peak_hypothesis_threshold_label =
    new Wt::WLabel( WString::tr( "bgw-peak-hypothesis-threshold-label" ), m_peak_hypothesis_threshold_container );
  m_peak_hypothesis_threshold_label->setWordWrap( false );
  m_peak_hypothesis_threshold = new NativeFloatSpinBox( m_peak_hypothesis_threshold_container );
  m_peak_hypothesis_threshold_label->setBuddy( m_peak_hypothesis_threshold );
  m_peak_hypothesis_threshold->setValue( 1.0f );// Default value from command line
  m_peak_hypothesis_threshold->setRange( 0.0f, 10.0f );
  m_peak_hypothesis_threshold->setSpinnerHidden( true );
  m_peak_hypothesis_threshold->setWidth( 40 );
  HelpSystem::attachToolTipOn(
    m_peak_hypothesis_threshold, WString::tr( "bgw-peak-hypothesis-threshold-tt" ), showToolTips );


  m_reports_container = new WGroupBox( WString::tr( "bgw-reports-grp-title" ), this );
  m_reports_container->addStyleClass( "ReportsContainer" );
  m_html_report = new Wt::WCheckBox( WString::tr( "bgw-reports-write-html" ), m_reports_container );
  m_html_report->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_html_report, WString::tr( "bgw-reports-write-html-tooltip" ), showToolTips );
  m_html_report->setChecked( true );

  m_create_csv_output = new Wt::WCheckBox( WString::tr( "bgw-create-csv-output" ), m_reports_container );
  m_create_csv_output->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_create_csv_output, WString::tr( "bgw-create-csv-output-tt" ), showToolTips );
  m_create_csv_output->setChecked( true );// Default to true as per command line


  Wt::WContainerWidget *custom_rpt_per_file_opts = new Wt::WContainerWidget( m_reports_container );
  custom_rpt_per_file_opts->addStyleClass( "CustomReportOptions" );
  m_per_file_custom_report =
    new Wt::WCheckBox( WString::tr( "bgw-reports-write-custom-per-file" ), custom_rpt_per_file_opts );
  m_per_file_custom_report->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn(
    m_per_file_custom_report, WString::tr( "bgw-reports-write-custom-per-file-tooltip" ), showToolTips );
  m_per_file_custom_report->setChecked( false );
  m_per_file_custom_report->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_per_file_custom_report->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_per_file_custom_report_container = new Wt::WContainerWidget( custom_rpt_per_file_opts );
  m_per_file_custom_report_container->addStyleClass( "CustomReportUploader EmptyReportUpload" );

  m_per_file_custom_report_resource = new FileDragUploadResource( m_per_file_custom_report_container );
  m_per_file_custom_report_container->doJavaScript( "BatchInputDropUploadSetup(" +
                                                    m_per_file_custom_report_container->jsRef() +
                                                    ", "
                                                    " '" +
                                                    m_per_file_custom_report_resource->url() + "');" );

  m_per_file_custom_report_container->doJavaScript( "setupOnDragEnterDom(['" +
                                                    m_per_file_custom_report_container->id() + "']);" );
  m_per_file_custom_report_resource->fileDrop().connect( boost::bind(
    &BatchGuiPeakFitWidget::perFileCustomReportUploaded, this, boost::placeholders::_1, boost::placeholders::_2 ) );
  m_per_file_custom_report_container->hide();

  // Group summary custom report elements together
  Wt::WContainerWidget *custom_rpt_summary_opts = new Wt::WContainerWidget( m_reports_container );
  custom_rpt_summary_opts->addStyleClass( "CustomReportOptions" );

  m_summary_custom_report =
    new Wt::WCheckBox( WString::tr( "bgw-reports-write-custom-summary" ), custom_rpt_summary_opts );
  m_summary_custom_report->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn(
    m_summary_custom_report, WString::tr( "bgw-reports-write-custom-summary-tooltip" ), showToolTips );
  m_summary_custom_report->setChecked( false );
  m_summary_custom_report->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_summary_custom_report->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_summary_custom_report_container = new Wt::WContainerWidget( custom_rpt_summary_opts );
  m_summary_custom_report_container->addStyleClass( "CustomReportUploader EmptyReportUpload" );

  m_summary_custom_report_resource = new FileDragUploadResource( m_summary_custom_report_container );
  m_summary_custom_report_container->doJavaScript( "BatchInputDropUploadSetup(" +
                                                   m_summary_custom_report_container->jsRef() +
                                                   ", "
                                                   " '" +
                                                   m_summary_custom_report_resource->url() + "');" );

  m_summary_custom_report_container->doJavaScript( "setupOnDragEnterDom(['" + m_summary_custom_report_container->id() +
                                                   "']);" );
  m_summary_custom_report_resource->fileDrop().connect( boost::bind(
    &BatchGuiPeakFitWidget::summaryCustomReportUploaded, this, boost::placeholders::_1, boost::placeholders::_2 ) );
  m_summary_custom_report_container->hide();

  // Connect signals to update analysis capability
  m_fit_all_peaks->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_fit_all_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_refit_energy_cal->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_refit_energy_cal->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_use_exemplar_energy_cal->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_use_exemplar_energy_cal->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_write_n42_with_results->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_write_n42_with_results->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_show_nonfit_peaks->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_show_nonfit_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_overwrite_output_files->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_overwrite_output_files->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_create_csv_output->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_create_csv_output->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_create_json_output->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_create_json_output->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_use_existing_background_peaks->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_use_existing_background_peaks->unChecked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_use_exemplar_energy_cal_for_background->checked().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_use_exemplar_energy_cal_for_background->changed().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_peak_stat_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_peak_stat_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  m_peak_hypothesis_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::optionsChanged );
  m_peak_hypothesis_threshold->valueChanged().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  optionsChanged();
}// BatchGuiPeakFitWidget constructor

BatchGuiPeakFitWidget::~BatchGuiPeakFitWidget()
{
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_background_file_drop->id() + "']);" );
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_exemplar_file_drop->id() + "']);" );
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_per_file_custom_report_container->id() + "']);" );
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_summary_custom_report_container->id() + "']);" );
}

void BatchGuiPeakFitWidget::handleFileUpload( WContainerWidget *dropArea, FileDragUploadResource *resource )
{
  dropArea->clear();

  const vector<tuple<string, string, bool>> spooled = resource->takeSpooledFiles();
  if( spooled.empty() )
    return;

  const tuple<string, string, bool> first_file = spooled.front();

  const string &display_name = std::get<0>( first_file );
  const string &path_to_file = std::get<1>( first_file );
  const bool should_delete = std::get<2>( first_file );
  const BatchGuiInputSpectrumFile::ShowPreviewOption show_preview = BatchGuiInputSpectrumFile::ShowPreviewOption::Show;

  BatchGuiInputSpectrumFile *input =
    new BatchGuiInputSpectrumFile( display_name, path_to_file, should_delete, show_preview, dropArea );
  dropArea->removeStyleClass( "EmptyExemplarUpload" );
  input->remove_self_request().connect(
    boost::bind( &BatchGuiPeakFitWidget::handle_remove_exemplar_upload, this, boost::placeholders::_1 ) );

  input->preview_created_signal().connect( this, &BatchGuiPeakFitWidget::optionsChanged );

  // Cleanup all uploads we wont use - dont expect this to actually happen
  for( size_t i = 1; i < spooled.size(); ++i )
  {
    if( get<2>( spooled[i] ) )
    {
      const bool success = SpecUtils::remove_file( std::get<1>( spooled[i] ) );
      if( !success )
        cerr << "exemplarUploaded: Warning, could not delete file '" << std::get<1>( spooled[i] ) << "'" << endl;
    }
  }

  optionsChanged();
}// void exemplarUploaded( std::string, std::string )

void BatchGuiPeakFitWidget::exemplarUploaded( const std::string &, const std::string & )
{
  handleFileUpload( m_exemplar_file_drop, m_exemplar_file_resource );
}

void BatchGuiPeakFitWidget::backgroundUploaded( const std::string &, const std::string & )
{
  handleFileUpload( m_background_file_drop, m_background_file_resource );
}

void BatchGuiPeakFitWidget::perFileCustomReportUploaded( const std::string &, const std::string & )
{
  const vector<tuple<string, string, bool>> spooled = m_per_file_custom_report_resource->takeSpooledFiles();

  for( const auto &file : spooled )
  {
    const string &display_name = std::get<0>( file );
    const string &path_to_file = std::get<1>( file );
    const bool should_delete = std::get<2>( file );

    BatchGuiInputFile *input =
      new BatchGuiInputFile( display_name, path_to_file, should_delete, m_per_file_custom_report_container );
    input->remove_self_request().connect( boost::bind(
      &BatchGuiPeakFitWidget::handle_remove_per_file_custom_report_upload, this, boost::placeholders::_1 ) );
  }// for( const auto &file : spooled )

  optionsChanged();
}// void perFileCustomReportUploaded( const std::string &, const std::string & )

void BatchGuiPeakFitWidget::summaryCustomReportUploaded( const std::string &, const std::string & )
{
  const vector<tuple<string, string, bool>> spooled = m_summary_custom_report_resource->takeSpooledFiles();

  for( const auto &file : spooled )
  {
    const string &display_name = std::get<0>( file );
    const string &path_to_file = std::get<1>( file );
    const bool should_delete = std::get<2>( file );

    BatchGuiInputFile *input =
      new BatchGuiInputFile( display_name, path_to_file, should_delete, m_summary_custom_report_container );
    input->remove_self_request().connect( boost::bind(
      &BatchGuiPeakFitWidget::handle_remove_summary_custom_report_upload, this, boost::placeholders::_1 ) );
  }// for( const auto &file : spooled )

  optionsChanged();
}// void summaryCustomReportUploaded( const std::string &, const std::string & )

void BatchGuiPeakFitWidget::handle_remove_per_file_custom_report_upload( BatchGuiInputFile *input )
{
  delete input;
  optionsChanged();
}

void BatchGuiPeakFitWidget::handle_remove_summary_custom_report_upload( BatchGuiInputFile *input )
{
  delete input;
  optionsChanged();
}

void BatchGuiPeakFitWidget::handle_remove_exemplar_upload( BatchGuiInputSpectrumFile *input )
{
  delete input;
  optionsChanged();
}

void BatchGuiPeakFitWidget::useCurrentBackgroundChanged()
{
  if( m_use_current_background->isChecked() )
    m_no_background->setChecked( false );

  optionsChanged();
}// void useCurrentBackgroundChanged()

void BatchGuiPeakFitWidget::useNoBackgroundChanged()
{
  if( m_no_background->isChecked() )
    m_use_current_background->setChecked( false );

  optionsChanged();
}// void useNoBackgroundChanged()

void BatchGuiPeakFitWidget::optionsChanged()
{
  const bool fit_all_peaks = m_fit_all_peaks->isChecked();
  m_refit_energy_cal->setHidden( fit_all_peaks );
  m_exemplar_input->setHidden( fit_all_peaks && !m_use_exemplar_energy_cal->isChecked() );
  // m_use_exemplar_energy_cal->setHidden( fit_all_peaks );
  m_background_input->setHidden( fit_all_peaks );
  m_show_nonfit_peaks->setHidden( fit_all_peaks );
  m_use_existing_background_peaks->setHidden( fit_all_peaks );
  m_use_exemplar_energy_cal_for_background->setHidden( fit_all_peaks );
  m_peak_stat_threshold_container->setHidden( fit_all_peaks );// Fit thresholds not currently implemented when fitting all peaks
  m_peak_hypothesis_threshold_container->setHidden( fit_all_peaks );// Fit thresholds not currently implemented when fitting all peaks

  const bool use_current_fore = m_use_current_foreground->isChecked();
  m_exemplar_file_drop->setHidden( use_current_fore );
  if( !use_current_fore )
  {
    if( !m_exemplar_file_drop->isVisible() )
      m_exemplar_file_drop->show();

    const bool has_drop_class = m_exemplar_file_drop->hasStyleClass( "EmptyExemplarUpload" );
    const size_t num_kids = m_exemplar_file_drop->children().size();
    if( ( num_kids > 0 ) && has_drop_class )
      m_exemplar_file_drop->removeStyleClass( "EmptyExemplarUpload" );
    else if( ( num_kids == 0 ) && !has_drop_class )
      m_exemplar_file_drop->addStyleClass( "EmptyExemplarUpload" );
  } else
  {
    if( m_exemplar_file_drop->isVisible() )
      m_exemplar_file_drop->hide();
  }// if( !use_current_fore ) / else


  const bool have_back = !!InterSpec::instance()->displayedHistogram( SpecUtils::SpectrumType::Background );

  const bool use_current_back = m_use_current_background->isChecked();
  const bool use_no_back = m_no_background->isChecked();
  assert( !use_current_back || ( use_current_back != use_no_back ) );// shouldnt both be checked at the same time
  auto hide_back_filedrop = [this]()
  {
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
  } else if( use_current_back )
  {
    assert( !use_no_back );
    hide_back_filedrop();
  } else
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
  const pair<bool,WString> ana_stat = canDoAnalysis();
  m_canDoAnalysis.emit( ana_stat.first, ana_stat.second );
}// void optionsChanged()

tuple<shared_ptr<SpecMeas>, string, set<int>> BatchGuiPeakFitWidget::get_exemplar() const
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
#if ( USE_REL_ACT_TOOL )
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
    }// if( foreground )
  } else
  {
    for( Wt::WWidget *child : m_exemplar_file_drop->children() )
    {
      BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
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

        //Find the sample with peaks.
        const set<set<int>> samples_with_peaks = meas->sampleNumsWithPeaks();
        if( samples_with_peaks.empty() )
        {
          // No peaks
        }else if( samples_with_peaks.size() == 1 )
        {
          samples = *begin(samples_with_peaks);
        }else
        {
          // Find the foreground with peaks
          for( const set<int> &peak_samples : samples_with_peaks )
          {
            bool is_back = false, is_fore = false;
            for( const int sample : peak_samples )
            {
              vector<shared_ptr<const SpecUtils::Measurement>> meass = meas->sample_measurements(sample);
              for( const shared_ptr<const SpecUtils::Measurement> &m : meass )
              {
                switch( m->source_type() )
                {
                  case SpecUtils::SourceType::IntrinsicActivity:
                  case SpecUtils::SourceType::Calibration:
                    break;

                  case SpecUtils::SourceType::Background:
                    is_back = true;
                    break;

                  case SpecUtils::SourceType::Foreground:
                  case SpecUtils::SourceType::Unknown:
                    is_fore = true;
                    break;
                }//switch( m->source_type() )
              }//for( const shared_ptr<const Measurement> &m : meass )

              if( is_fore && !is_back )
              {
                samples = peak_samples;
                break;
              }
            }//for( const int sample : samples )
          }//for( const set<int> &samples : samples_with_peaks )
        }//for( const set<int> &peak_samples : samples_with_peaks )

        if( samples.empty() )
        {
          //It is ambigous.
          // TODO: Put more effort into trying to divine what set of sample numbers is intended to be the foreground.
          //  Right now we'll just leave `samples` empty, and the user will get an error message that there are now peaks in exemplar...
        }
      }//if( meas )

      if( meas )
        break;
    }//for( Wt::WWidget *child : m_exemplar_file_drop->children() )
  }// if( m_use_current_foreground->isChecked() ) / else


  return { meas, filename, samples };
}// tuple<shared_ptr<const SpecMeas>,string,set<int>> get_exemplar() const

tuple<shared_ptr<const SpecMeas>, string, set<int>> BatchGuiPeakFitWidget::get_background() const
{
  string filename;
  set<int> samples;
  shared_ptr<const SpecMeas> meas;

  if( m_no_background->isChecked() )
  {
    // Nothing to do here
  } else if( m_use_current_background->isChecked() )
  {
    InterSpec *interspec = InterSpec::instance();
    meas = interspec->measurment( SpecUtils::SpectrumType::Background );
    samples = interspec->displayedSamples( SpecUtils::SpectrumType::Background );
    if( meas )
      filename = meas->filename();
  } else
  {
    for( Wt::WWidget *child : m_background_file_drop->children() )
    {
      BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
      assert( input_file );
      if( !input_file )
        continue;

      meas = input_file->spec_meas();
      filename = input_file->path_to_file();
      break;
    }
  }// if( m_use_current_foreground->isChecked() ) / else


  return { meas, filename, samples };
}// tuple<shared_ptr<const SpecMeas>,string,set<int>> get_background() const

BatchPeak::BatchPeakFitOptions BatchGuiPeakFitWidget::getPeakFitOptions() const
{
  BatchPeak::BatchPeakFitOptions answer;

  answer.to_stdout = false;
  answer.fit_all_peaks = m_fit_all_peaks->isChecked();
  if( answer.fit_all_peaks )
  {
    answer.refit_energy_cal = false;
    answer.show_nonfit_peaks = false;
  } else
  {
    answer.refit_energy_cal = m_refit_energy_cal->isChecked();
    answer.show_nonfit_peaks = m_show_nonfit_peaks->isChecked();
    answer.peak_stat_threshold = m_peak_stat_threshold->value();
    answer.peak_hypothesis_threshold = m_peak_hypothesis_threshold->value();
  }

  tuple<shared_ptr<const SpecMeas>, string, set<int>> exemplar_info = get_exemplar();

  if( get<0>( exemplar_info ) )
    answer.use_exemplar_energy_cal = m_use_exemplar_energy_cal->isChecked();

  answer.write_n42_with_results = m_write_n42_with_results->isChecked();
  answer.overwrite_output_files = m_overwrite_output_files->isChecked();
  answer.create_csv_output = m_create_csv_output->isChecked();
  answer.create_json_output = m_create_json_output->isChecked();

  if( m_no_background->isChecked() )
  {
    answer.use_existing_background_peaks = false;
    answer.use_exemplar_energy_cal_for_background = false;
  } else
  {
    const tuple<shared_ptr<const SpecMeas>, string, set<int>> background_info = get_background();
    const shared_ptr<const SpecMeas> &background_spec = get<0>( background_info );
    if( background_spec )
    {
      shared_ptr<SpecMeas> spec_copy = make_shared<SpecMeas>();
      spec_copy->uniqueCopyContents( *background_spec );
      answer.cached_background_subtract_spec = spec_copy;
    }

    answer.background_subtract_file = get<1>( background_info );
    answer.background_subtract_samples = get<2>( background_info );
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

  if( m_create_csv_output->isChecked() )
  {
    // answer.report_templates.push_back( "csv" );
    answer.summary_report_templates.push_back( "csv" );
  }

  if( m_per_file_custom_report->isChecked() )
  {
    for( Wt::WWidget *child : m_per_file_custom_report_container->children() )
    {
      BatchGuiInputFile *input_file = dynamic_cast<BatchGuiInputFile *>( child );
      assert( input_file );
      if( !input_file )
        continue;

      string tmplt_name = input_file->path_to_file();
      if( !input_file->display_name().empty()
         && (SpecUtils::filename(input_file->display_name()) != SpecUtils::filename(input_file->path_to_file())) )
      {
        tmplt_name += BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker;
        tmplt_name += SpecUtils::filename(input_file->display_name());
      }

      answer.report_templates.push_back( tmplt_name );
    }
  }

  if( m_summary_custom_report->isChecked() )
  {
    for( Wt::WWidget *child : m_summary_custom_report_container->children() )
    {
      BatchGuiInputFile *input_file = dynamic_cast<BatchGuiInputFile *>( child );
      assert( input_file );
      if( !input_file )
        continue;


      string tmplt_name = input_file->path_to_file();
      if( !input_file->display_name().empty()
         && (SpecUtils::filename(input_file->display_name()) != SpecUtils::filename(input_file->path_to_file())) )
      {
        tmplt_name += BatchPeak::BatchPeakFitOptions::sm_report_display_name_marker;
        tmplt_name += SpecUtils::filename(input_file->display_name());
      }

      answer.summary_report_templates.push_back( tmplt_name );
    }
  }// if( m_summary_custom_report->isChecked() )

  // answer.template_include_dir = ...; //std::string - if we allow setting this, then the custom report templates
  // cant be used. answer.output_dir = ...; //std::string

  return answer;
}// BatchPeak::BatchPeakFitOptions getOptions() const

void BatchGuiPeakFitWidget::performAnalysis(
  const vector<tuple<string, string, std::shared_ptr<const SpecMeas>>> &input_files,
  const string &output_dir )
{
  BatchPeak::BatchPeakFitOptions options = getPeakFitOptions();
  options.output_dir = output_dir;

  const tuple<shared_ptr<const SpecMeas>, string, set<int>> exemplar_info = get_exemplar();
  const shared_ptr<const SpecMeas> &exemplar =
    get<0>( exemplar_info );// exemplar is nullptr if a peaks CSV file is used.
  const string &exemplar_filename = get<1>( exemplar_info );
  const set<int> &exemplar_samples = get<2>( exemplar_info );

  const tuple<shared_ptr<const SpecMeas>, string, set<int>> background_info = get_background();

  vector<string> file_names;
  vector<shared_ptr<SpecMeas>> spec_files;
  for( size_t i = 0; i < input_files.size(); ++i )
  {
    const string &display_name = get<0>( input_files[i] );
    const string &path_to_file = get<1>( input_files[i] );
    const shared_ptr<const SpecMeas> &spec_meas = get<2>( input_files[i] );
    assert( spec_meas );
    if( !spec_meas )
      continue;

    file_names.push_back( spec_meas->filename() );

    shared_ptr<SpecMeas> spec_copy =
      make_shared<SpecMeas>();// we need it non-const (and we dont want to mess with in-memory files, so we get
                              // consistent resutls with multiple calls)
    spec_copy->uniqueCopyContents( *spec_meas );
    spec_files.push_back( spec_copy );
  }// for( size_t i = 0; i < input_files.size(); ++i )

  // We need to do the work off the main GUI thread, so we need to post the work to the server.
  const string sessionid = Wt::WApplication::instance()->sessionId();
  auto error_msg = make_shared<string>();
  auto results = make_shared<BatchPeak::BatchPeakFitSummary>();

  SimpleDialog *waiting_dialog =
    new SimpleDialog( WString::tr( "bgw-performing-work-title" ), WString::tr( "bgw-performing-work-msg" ) );
  waiting_dialog->addButton( WString::tr( "Close" ) );
  boost::function<void( void )> close_waiting_dialog =
    wApp->bind( boost::bind( &SimpleDialog::done, waiting_dialog, Wt::WDialog::DialogCode::Accepted ) );


  std::function<void( void )> show_error_dialog = [error_msg, close_waiting_dialog]()
  {
    close_waiting_dialog();
    SimpleDialog *dialog = new SimpleDialog( WString::tr( "bgw-error-analysis-title" ),
                                             WString::tr( "bgw-error-analysis-msg" ).arg( *error_msg ) );
    dialog->addStyleClass( "BatchAnalysisErrorDialog" );
    dialog->addButton( WString::tr( "Okay" ) );
    wApp->triggerUpdate();
  };

  std::function<void( void )> update_gui_fcn = [results, options, close_waiting_dialog]()
  {
    close_waiting_dialog();

    SimpleDialog *dialog = new SimpleDialog( WString::tr( "bgw-analysis-summary-title" ) );
    dialog->addStyleClass( "BatchAnalysisResultDialog" );

    try
    {
      nlohmann::json summary_json = nlohmann::json::parse( results->summary_json );
      inja::Environment env = BatchInfoLog::get_default_inja_env( options );
      string rpt = BatchInfoLog::render_template(
        "html", env, BatchInfoLog::TemplateRenderType::PeakFitSummary, options, summary_json );

      // Create the resource and get its URL
      BatchReportResource *report_resource = new BatchReportResource( std::move( rpt ), dialog );
      const string resource_url = report_resource->url();

      const string contents = "<iframe "
                              "style=\"width: 100%; height: 100%; border: none;\" "
                              "sandbox=\"allow-scripts\" "
                              "src=\"" +
                              resource_url +
                              "\" "
                              "onerror=\"console.error('Iframe failed to load')\"> "
                              "</iframe>";

      dialog->resize( WLength( 80.0, WLength::Percentage ), WLength( 95.0, WLength::Percentage ) );

      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportIFrameHolder" );
      result_text->setWidth( WLength( 100.0, WLength::Percentage ) );
      result_text->setHeight( WLength( 100.0, WLength::Percentage ) );
    } catch( nlohmann::json::parse_error &e )
    {
      const string contents = "<p>" + WString::tr( "bgw-json-parse-error" ).arg( string( e.what() ) ).toUTF8() + "</p>";
      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportJsonError" );
    } catch( inja::InjaError &e )
    {
      const string contents = "<p>" +
                              WString::tr( "bgw-inja-error" )
                                .arg( e.message )
                                .arg( std::to_string( e.location.line ) )
                                .arg( std::to_string( e.location.column ) )
                                .toUTF8() +
                              "</p>";

      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportInjaError" );
    } catch( std::exception &e )
    {
      const string contents = "<p>" + WString::tr( "bgw-misc-error" ).arg( string( e.what() ) ).toUTF8() + "</p>";
      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportMiscError" );
    }

    dialog->addButton( WString::tr( "Okay" ) );


    if( !results->warnings.empty() )
    {
      SimpleDialog *warnings_dialog = new SimpleDialog( WString::tr( "bgw-warning-title" ) );
      warnings_dialog->addStyleClass( "BatchAnalysisWarningDialog" );

      InterSpec *interspec = InterSpec::instance();
      const int app_width = interspec->renderedWidth();
      const double warn_width = std::min( 450, (app_width > 100) ? app_width : 450 );
      warnings_dialog->setMinimumSize( WLength(warn_width,WLength::Pixel), WLength::Auto );

      WContainerWidget *contents = warnings_dialog->contents();
      contents->setList( true, false );
      set<string> seen_warnings;

      for( const string &warn_msg : results->warnings )
      {
        // Avoid duplicate warnings
        if( seen_warnings.count(warn_msg) )
          continue;
        seen_warnings.insert( warn_msg );

        WContainerWidget *item = new WContainerWidget( contents );
        new WText( warn_msg, item );
      }

      warnings_dialog->addButton( WString::tr( "Okay" ) );
    }// if( !results.warnings.empty() )

    wApp->triggerUpdate();
  };// update_gui_fcn lambda


  std::function<void( void )> do_work_fcn = [exemplar_filename,
                                             exemplar,
                                             exemplar_samples,
                                             file_names,
                                             spec_files,
                                             options,
                                             results,
                                             error_msg,
                                             update_gui_fcn,
                                             show_error_dialog,
                                             sessionid]()
  {
    try
    {
      BatchPeak::fit_peaks_in_files(
        exemplar_filename, exemplar, exemplar_samples, file_names, spec_files, options, results.get() );

      WServer::instance()->post( sessionid, update_gui_fcn );
    } catch( std::exception &e )
    {
      *error_msg = e.what();
      WServer::instance()->post( sessionid, show_error_dialog );
    }
  };// do_work_fcn lambda

  // Do the peak fitting work in a non-gui thread.
  WServer::instance()->ioService().boost::asio::io_service::post( do_work_fcn );
}// performAnalysis(...)

std::pair<bool,Wt::WString> BatchGuiPeakFitWidget::canDoAnalysis() const
{
  if( !m_fit_all_peaks->isVisible() || !m_fit_all_peaks->isChecked() )
  {
    const tuple<shared_ptr<SpecMeas>, string, set<int>> exemplar_info = get_exemplar();
    const shared_ptr<SpecMeas> &exemplar = get<0>( exemplar_info );
    if( exemplar )
    {
      const set<int> &exemplar_samples = get<2>( exemplar_info );
      std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = exemplar->peaks( exemplar_samples );
      if( !peaks || peaks->empty() )
        return {false, Wt::WString::tr("bgw-no-ana-peak-fit-no-peaks") };
    } else
    {
      // Might be a peaks CSV file
      bool has_peaks_csv = false;
      for( Wt::WWidget *child : m_exemplar_file_drop->children() )
      {
        BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
        assert( input_file );
        if( input_file && input_file->is_peaks_csv() )
        {
          has_peaks_csv = true;
          break;
        }
      }//for( Wt::WWidget *child : m_exemplar_file_drop->children() )

      if( !has_peaks_csv )
        return {false, Wt::WString::tr("bgw-no-ana-peak-fit-need-peaks") };
    }// if( exemplar ) / else
  }// if( !m_fit_all_peaks->isVisible() || !m_fit_all_peaks->isChecked() )

  bool back_status_okay = m_no_background->isChecked();
  if( !back_status_okay && m_use_current_background->isChecked() )
  {
    back_status_okay = !!( InterSpec::instance()->displayedHistogram( SpecUtils::SpectrumType::Background ) );
  }

  if( !back_status_okay && !m_use_current_background->isChecked() && !m_no_background->isChecked() )
  {
    for( Wt::WWidget *child : m_background_file_drop->children() )
    {
      BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
      back_status_okay = ( input_file && input_file->spec_meas() );
      if( back_status_okay )
        break;
    }
  }

  if( !back_status_okay )
    return {false, Wt::WString::tr("bgw-no-ana-peak-fit-bad-background") };

  return {true, WString()};
}//std::pair<bool,Wt::WString> BatchGuiPeakFitWidget::canDoAnalysis() const

BatchGuiActShieldAnaWidget::BatchGuiActShieldAnaWidget( Wt::WContainerWidget *parent )
: BatchGuiPeakFitWidget( parent ),
  m_act_shield_container( nullptr ),
  m_use_bq( nullptr ),
  m_hard_background_sub( nullptr ),
  m_detector_input( nullptr ),
  m_use_detector_override( nullptr ),
  m_detector_file_drop( nullptr ),
  m_detector_file_resource( nullptr ),
  m_uploaded_detector( nullptr ),
  m_override_distance( nullptr ),
  m_distance_input_container( nullptr ),
  m_distance_label( nullptr ),
  m_distance_edit( nullptr ),
  m_csv_report( nullptr )
{
  addStyleClass( "BatchGuiActShieldAnaWidget" );

  m_fit_all_peaks->setChecked( false );
  m_fit_all_peaks->hide();

  m_act_shield_container = new Wt::WContainerWidget();
  m_act_shield_container->addStyleClass( "ActShieldOptionsContainer" );
  const int back_area_index = indexOf( m_background_input );

  insertWidget( back_area_index + 1, m_act_shield_container );


  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  m_csv_report = new Wt::WCheckBox( WString::tr( "bgw-reports-write-csv" ) );
  m_csv_report->addStyleClass( "CbNoLineBreak" );
  HelpSystem::attachToolTipOn( m_csv_report, WString::tr( "bgw-reports-write-csv-tooltip" ), showToolTips );
  m_reports_container->insertWidget( 2, m_csv_report );

  m_use_bq = new WCheckBox( WString::tr( "bgw-use-bq" ) );
  m_peak_options_container->insertWidget( 0, m_use_bq );
  m_use_bq->addStyleClass( "CbNoLineBreak" );
  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  m_use_bq->setChecked( useBq );
  m_use_bq->checked().connect( this, &BatchGuiActShieldAnaWidget::useBqChanged );
  m_use_bq->unChecked().connect( this, &BatchGuiActShieldAnaWidget::useBqChanged );


  m_hard_background_sub = new WCheckBox( WString::tr( "bgw-hard-background-sub" ), m_background_input );
  m_hard_background_sub->addStyleClass( "CbNoLineBreak" );
  m_hard_background_sub->setChecked( false );
  m_hard_background_sub->checked().connect( this, &BatchGuiActShieldAnaWidget::useHardBackgroundSubChanged );
  m_hard_background_sub->unChecked().connect( this, &BatchGuiActShieldAnaWidget::useHardBackgroundSubChanged );

  m_detector_input = new WGroupBox( WString::tr( "bgw-detector-input-label" ), m_act_shield_container );
  m_detector_input->addStyleClass( "DetectorInputContainer" );

  m_use_detector_override = new WCheckBox( WString::tr( "bgw-use-detector-override" ), m_detector_input );
  m_use_detector_override->addStyleClass( "CbNoLineBreak" );
  m_use_detector_override->setChecked( false );
  m_use_detector_override->checked().connect( this, &BatchGuiActShieldAnaWidget::useDetectorOverrideChanged );
  m_use_detector_override->unChecked().connect( this, &BatchGuiActShieldAnaWidget::useDetectorOverrideChanged );

  m_detector_file_drop = new WContainerWidget( m_detector_input );
  m_detector_file_drop->addStyleClass( "DetectorFileDrop EmptyDetectorUpload" );
  m_detector_file_drop->setHidden( true );

  m_detector_file_resource = new FileDragUploadResource( m_detector_file_drop );
  m_detector_file_drop->doJavaScript( "BatchInputDropUploadSetup(" + m_detector_file_drop->jsRef() +
                                      ", "
                                      " '" +
                                      m_detector_file_resource->url() + "');" );
  doJavaScript( "setupOnDragEnterDom(['" + m_detector_file_drop->id() + "']);" );
  m_detector_file_resource->fileDrop().connect( boost::bind(
    &BatchGuiActShieldAnaWidget::detectorUploaded, this, boost::placeholders::_1, boost::placeholders::_2 ) );

  // Distance override
  m_override_distance = new WCheckBox( WString::tr( "bgw-override-distance" ), m_act_shield_container );
  m_override_distance->addStyleClass( "CbNoLineBreak" );
  m_override_distance->setChecked( false );
  m_override_distance->checked().connect( this, &BatchGuiActShieldAnaWidget::overrideDistanceChanged );
  m_override_distance->unChecked().connect( this, &BatchGuiActShieldAnaWidget::overrideDistanceChanged );

  m_distance_input_container = new Wt::WContainerWidget( m_act_shield_container );
  m_distance_input_container->addStyleClass( "DistanceInputContainer" );
  m_distance_input_container->setHidden( true );

  m_distance_label = new Wt::WLabel( WString::tr( "bgw-distance-label" ), m_distance_input_container );
  m_distance_label->setWordWrap( false );
  m_distance_edit = new Wt::WLineEdit( m_distance_input_container );
  m_distance_label->setBuddy( m_distance_edit );
  m_distance_edit->setWidth( 100 );

  WRegExpValidator *distValidator =
    new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, m_distance_edit );
  m_distance_edit->setValidator( distValidator );
  m_distance_edit->changed().connect( this, &BatchGuiActShieldAnaWidget::distanceValueChanged );
  m_distance_edit->enterPressed().connect( this, &BatchGuiActShieldAnaWidget::distanceValueChanged );

  optionsChanged();
}// BatchGuiActShieldAnaWidget()

BatchActivity::BatchActivityFitOptions BatchGuiActShieldAnaWidget::getActivityFitOptions() const
{
  BatchActivity::BatchActivityFitOptions options;

  BatchPeak::BatchPeakFitOptions &peak_options = static_cast<BatchPeak::BatchPeakFitOptions &>( options );
  peak_options = getPeakFitOptions();

  options.use_bq = m_use_bq->isChecked();
  options.hard_background_sub = ( m_hard_background_sub->isVisible() && m_hard_background_sub->isChecked() );

  const shared_ptr<const DetectorPeakResponse> det = detector();
  if( m_use_detector_override->isChecked() && m_uploaded_detector )
    options.drf_override = m_uploaded_detector;

  if( ( !det || !det->isFixedGeometry() ) &&
      ( !m_override_distance->isHidden() && m_override_distance->isChecked() && !m_distance_edit->text().empty() ) )
  {
    try
    {
      std::string distance_str = m_distance_edit->text().toUTF8();

      const double distance = PhysicalUnits::stringToDistance( distance_str );
      if( distance < 0.0 )
        throw runtime_error( "Distance must be positive" );

      options.distance_override = distance;
    } catch( std::exception & )
    {
      assert( 0 );
    }// try / catch to get the distance
  }// if( !m_override_distance->isHidden() && m_override_distance->isChecked() && !m_distance_edit->text().empty() )

  return options;
}// getActivityFitOptions()

void BatchGuiActShieldAnaWidget::performAnalysis(
  const vector<tuple<string, string, std::shared_ptr<const SpecMeas>>> &input_files,
  const string &output_dir )
{
  BatchActivity::BatchActivityFitOptions options = getActivityFitOptions();
  options.output_dir = output_dir;

  const tuple<shared_ptr<SpecMeas>, string, set<int>> exemplar_info = get_exemplar();
  const shared_ptr<SpecMeas> &exemplar = get<0>( exemplar_info );
  const string &exemplar_filename = get<1>( exemplar_info );
  const set<int> &exemplar_samples = get<2>( exemplar_info );

  if( !exemplar )
  {
    SimpleDialog *dialog =
      new SimpleDialog( WString::tr( "bgw-error-analysis-title" ), WString::tr( "bgw-no-exemplar-msg" ) );
    dialog->addStyleClass( "BatchAnalysisErrorDialog" );
    dialog->addButton( WString::tr( "Okay" ) );
    return;
  }

  vector<string> file_names;
  vector<shared_ptr<SpecMeas>> input_files_meas;
  for( size_t i = 0; i < input_files.size(); ++i )
  {
    const string &display_name = get<0>( input_files[i] );
    const string &path_to_file = get<1>( input_files[i] );
    const shared_ptr<const SpecMeas> &spec_meas = get<2>( input_files[i] );
    assert( spec_meas );
    if( !spec_meas )
      continue;

    shared_ptr<SpecMeas> meas_copy = make_shared<SpecMeas>();
    meas_copy->uniqueCopyContents( *spec_meas );
    meas_copy->set_filename( display_name );

    file_names.push_back( display_name );
    input_files_meas.push_back( meas_copy );
  }// for( size_t i = 0; i < input_files.size(); ++i )

  shared_ptr<BatchActivity::BatchActivityFitSummary> summary_results =
    make_shared<BatchActivity::BatchActivityFitSummary>();

  // We need to do the work off the main GUI thread, so we need to post the work to the server.
  const string sessionid = Wt::WApplication::instance()->sessionId();
  auto error_msg = make_shared<string>();

  SimpleDialog *waiting_dialog =
    new SimpleDialog( WString::tr( "bgw-performing-work-title" ), WString::tr( "bgw-performing-work-msg" ) );
  waiting_dialog->addButton( WString::tr( "Close" ) );
  boost::function<void( void )> close_waiting_dialog =
    wApp->bind( boost::bind( &SimpleDialog::done, waiting_dialog, Wt::WDialog::DialogCode::Accepted ) );

  std::function<void( void )> show_error_dialog = [error_msg, close_waiting_dialog]()
  {
    close_waiting_dialog();
    SimpleDialog *dialog = new SimpleDialog( WString::tr( "bgw-error-analysis-title" ),
                                             WString::tr( "bgw-error-analysis-msg" ).arg( *error_msg ) );
    dialog->addStyleClass( "BatchAnalysisErrorDialog" );
    dialog->addButton( WString::tr( "Okay" ) );
    wApp->triggerUpdate();
  };

  std::function<void( void )> update_gui_fcn = [close_waiting_dialog, summary_results, options]()
  {
    close_waiting_dialog();

    SimpleDialog *dialog = new SimpleDialog( WString::tr( "bgw-analysis-summary-title" ) );
    dialog->addStyleClass( "BatchAnalysisResultDialog" );

    try
    {
      nlohmann::json summary_json = nlohmann::json::parse( summary_results->summary_json );
      inja::Environment env = BatchInfoLog::get_default_inja_env( summary_results->options );
      string rpt = BatchInfoLog::render_template(
        "html", env, BatchInfoLog::TemplateRenderType::ActShieldSummary, summary_results->options, summary_json );

      // Create the resource and get its URL
      BatchReportResource *report_resource = new BatchReportResource( std::move( rpt ), dialog );
      const string resource_url = report_resource->url();

      const string contents = "<iframe "
                              "style=\"width: 100%; height: 100%; border: none;\" "
                              "sandbox=\"allow-scripts\" "
                              "src=\"" +
                              resource_url +
                              "\" "
                              "onerror=\"console.error('Iframe failed to load')\"> "
                              "</iframe>";

      dialog->resize( WLength( 80.0, WLength::Percentage ), WLength( 95.0, WLength::Percentage ) );

      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportIFrameHolder" );
      result_text->setWidth( WLength( 100.0, WLength::Percentage ) );
      result_text->setHeight( WLength( 100.0, WLength::Percentage ) );
    } catch( nlohmann::json::parse_error &e )
    {
      const string contents = "<p>" + WString::tr( "bgw-json-parse-error" ).arg( string( e.what() ) ).toUTF8() + "</p>";
      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportJsonError" );
    } catch( inja::InjaError &e )
    {
      const string contents = "<p>" +
                              WString::tr( "bgw-inja-error" )
                                .arg( e.message )
                                .arg( std::to_string( e.location.line ) )
                                .arg( std::to_string( e.location.column ) )
                                .toUTF8() +
                              "</p>";

      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportInjaError" );
    } catch( std::exception &e )
    {
      const string contents = "<p>" + WString::tr( "bgw-misc-error" ).arg( string( e.what() ) ).toUTF8() + "</p>";
      WText *result_text = new WText( contents, Wt::TextFormat::XHTMLUnsafeText, dialog->contents() );
      result_text->addStyleClass( "BatchReportMiscError" );
    }

    dialog->addButton( WString::tr( "Okay" ) );


    if( !summary_results->warnings.empty() )
    {
      SimpleDialog *warnings_dialog = new SimpleDialog( WString::tr( "bgw-warning-title" ) );
      warnings_dialog->addStyleClass( "BatchAnalysisWarningDialog" );

      WContainerWidget *contents = warnings_dialog->contents();
      contents->setList( true, false );
      for( const string &warn_msg : summary_results->warnings )
      {
        WContainerWidget *item = new WContainerWidget( contents );
        new WText( warn_msg, item );
      }

      warnings_dialog->addButton( WString::tr( "Okay" ) );
    }// if( !results.warnings.empty() )

    wApp->triggerUpdate();
  };// update_gui_fcn lambda

  std::function<void( void )> do_work_fcn = [exemplar_filename,
                                             exemplar,
                                             exemplar_samples,
                                             file_names,
                                             input_files_meas,
                                             options,
                                             summary_results,
                                             error_msg,
                                             update_gui_fcn,
                                             show_error_dialog,
                                             sessionid]()
  {
    try
    {
      BatchActivity::fit_activities_in_files(
        exemplar_filename, exemplar, exemplar_samples, file_names, input_files_meas, options, summary_results.get() );
      WServer::instance()->post( sessionid, update_gui_fcn );
    } catch( std::exception &e )
    {
      *error_msg = e.what();
      WServer::instance()->post( sessionid, show_error_dialog );
    }
  };// do_work_fcn lambda

  // Do the activity fitting work in a non-gui thread.
  WServer::instance()->ioService().boost::asio::io_service::post( do_work_fcn );
}// performAnalysis(...)

pair<bool,Wt::WString> BatchGuiActShieldAnaWidget::canDoAnalysis() const
{
  const pair<bool,Wt::WString> peak_fit_status = BatchGuiPeakFitWidget::canDoAnalysis();
  if( !peak_fit_status.first )
    return peak_fit_status;

  const shared_ptr<const DetectorPeakResponse> det = detector();
  if( !det )
    return {false, WString::tr("bgw-no-ana-act-shield-no-drf")};

  // `BatchGuiPeakFitWidget::canDoAnalysis()` should have already checked for the exemplar,
  //   and that it has peaks, but peaks analysis will allow the exemplar to be a peaks CSV file.
  const tuple<shared_ptr<SpecMeas>, string, set<int>> exemplar_info = get_exemplar();
  const shared_ptr<SpecMeas> &exemplar = get<0>( exemplar_info );
  if( !exemplar )
    return {false, WString::tr("bgw-no-ana-act-shield-no-exemplar")};

  const set<int> &exemplar_samples = get<2>( exemplar_info );
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = exemplar->peaks( exemplar_samples );
  if( !peaks || peaks->empty() )
    return {false, WString::tr("bgw-no-ana-act-shield-no-peaks-in-exemplar")};

  const rapidxml::xml_document<char> *shielding_source_model = exemplar->shieldingSourceModel();
  if( !shielding_source_model )
    return {false, WString::tr("bgw-no-ana-act-shield-no-act-fit-model")};

  if( !m_no_background->isChecked() )
  {
    const tuple<shared_ptr<const SpecMeas>, string, set<int>> background_info = get_background();
    const shared_ptr<const SpecMeas> &background = get<0>( background_info );
    if( !background )
      return {false, WString::tr("bgw-no-ana-act-shield-invalid-background")};
  }// if( !m_no_background->isChecked() )

  return {true, WString()};
}//pair<bool,Wt::WString> BatchGuiActShieldAnaWidget::canDoAnalysis() const

void BatchGuiActShieldAnaWidget::optionsChanged()
{
  BatchGuiPeakFitWidget::optionsChanged();

  // Figure out what detector response function we are using, and if a fixed
  //  geometery, hide the distance input.
  const shared_ptr<const DetectorPeakResponse> det = detector();
  if( det && det->isFixedGeometry() )
  {
    m_distance_input_container->setHidden( true );
    m_override_distance->setHidden( true );
  }

  m_hard_background_sub->setHidden( m_no_background->isChecked() );
}// optionsChanged()

void BatchGuiActShieldAnaWidget::useBqChanged()
{
  optionsChanged();
}

void BatchGuiActShieldAnaWidget::useHardBackgroundSubChanged()
{
  optionsChanged();
}

void BatchGuiActShieldAnaWidget::useDetectorOverrideChanged()
{
  const bool use_detector_override = m_use_detector_override->isChecked();
  m_detector_file_drop->setHidden( !use_detector_override );

  if( use_detector_override )
  {
    const bool has_drop_class = m_detector_file_drop->hasStyleClass( "EmptyDetectorUpload" );
    const size_t num_kids = m_detector_file_drop->children().size();
    if( ( num_kids > 0 ) && has_drop_class )
      m_detector_file_drop->removeStyleClass( "EmptyDetectorUpload" );
    else if( ( num_kids == 0 ) && !has_drop_class )
      m_detector_file_drop->addStyleClass( "EmptyDetectorUpload" );
  } else
  {
    // Clear any uploaded files when disabled
    const size_t num_kids = m_detector_file_drop->children().size();
    if( num_kids )
      m_detector_file_drop->clear();
    const bool has_drop_class = m_detector_file_drop->hasStyleClass( "EmptyDetectorUpload" );
    if( !has_drop_class )
      m_detector_file_drop->addStyleClass( "EmptyDetectorUpload" );
  }

  optionsChanged();
}

std::shared_ptr<const DetectorPeakResponse> BatchGuiActShieldAnaWidget::detector() const
{
  if( m_use_detector_override->isChecked() && m_uploaded_detector )
    return m_uploaded_detector;

  if( m_use_current_foreground->isChecked() )
  {
    InterSpec *interspec = InterSpec::instance();
    assert( interspec );
    if( !interspec )
      return nullptr;

    const shared_ptr<SpecMeas> spec_meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
    if( !spec_meas )
      return nullptr;

    return spec_meas->detector();
  }// if( m_use_current_foreground->isChecked() )

  // Try to find a detector in the exemplar file
  for( Wt::WWidget *child : m_exemplar_file_drop->children() )
  {
    BatchGuiInputSpectrumFile *input_file = dynamic_cast<BatchGuiInputSpectrumFile *>( child );
    const shared_ptr<SpecMeas> spec_meas = input_file ? input_file->spec_meas() : nullptr;
    if( spec_meas )
      return spec_meas->detector();
  }

  return nullptr;
}// detector()

void BatchGuiActShieldAnaWidget::detectorUploaded( const std::string &, const std::string & )
{
  const vector<tuple<string, string, bool>> spooled = m_detector_file_resource->takeSpooledFiles();
  if( spooled.empty() )
    return;

  // Clear any existing detector files
  m_detector_file_drop->clear();
  m_uploaded_detector = nullptr;

  for( const auto &file : spooled )
  {
    const string &display_name = std::get<0>( file );
    const string &path_to_file = std::get<1>( file );
    const bool should_delete = std::get<2>( file );

    try
    {
      // Try to load the detector response function
      m_uploaded_detector = BatchActivity::init_drf_from_name( path_to_file, "" );

      if( m_uploaded_detector )
      {
        BatchGuiInputFile *input =
          new BatchGuiInputFile( display_name, path_to_file, should_delete, m_detector_file_drop );
        input->remove_self_request().connect(
          boost::bind( &BatchGuiActShieldAnaWidget::handle_remove_detector_upload, this, boost::placeholders::_1 ) );

        break;// Only use the first valid DRF file
      } else
      {
        // Clean up file if it couldn't be loaded as a DRF
        if( should_delete )
          SpecUtils::remove_file( path_to_file );
      }
    } catch( std::exception &e )
    {
      // Clean up file if it couldn't be loaded as a DRF
      if( should_delete )
        SpecUtils::remove_file( path_to_file );
    }
  }// for( const auto &file : spooled )

  // Clean up any remaining files we didn't use
  for( size_t i = 1; i < spooled.size(); ++i )
  {
    if( get<2>( spooled[i] ) )
    {
      const bool success = SpecUtils::remove_file( std::get<1>( spooled[i] ) );
      if( !success )
        cerr << "detectorUploaded: Warning, could not delete file '" << std::get<1>( spooled[i] ) << "'" << endl;
    }
  }

  const bool has_drop_class = m_detector_file_drop->hasStyleClass( "EmptyDetectorUpload" );
  const size_t num_kids = m_detector_file_drop->children().size();
  if( ( num_kids == 0 ) && !has_drop_class )
    m_detector_file_drop->addStyleClass( "EmptyDetectorUpload" );
  else if( ( num_kids > 0 ) && has_drop_class )
    m_detector_file_drop->removeStyleClass( "EmptyDetectorUpload" );

  optionsChanged();
}

void BatchGuiActShieldAnaWidget::handle_remove_detector_upload( BatchGuiInputFile *input )
{
  delete input;
  m_uploaded_detector.reset();

  const bool has_drop_class = m_detector_file_drop->hasStyleClass( "EmptyDetectorUpload" );
  if( !has_drop_class )
    m_detector_file_drop->addStyleClass( "EmptyDetectorUpload" );

  optionsChanged();
}

void BatchGuiActShieldAnaWidget::overrideDistanceChanged()
{
  m_distance_input_container->setHidden( m_override_distance->isHidden() || !m_override_distance->isChecked() );
  optionsChanged();
}

void BatchGuiActShieldAnaWidget::distanceValueChanged()
{
  string distance_str = m_distance_edit->text().toUTF8();

  // If no units specified, default to cm
  if( distance_str.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
  {
    try
    {
      PhysicalUnits::stringToDistance( distance_str + " cm" );
      distance_str += " cm";
      m_distance_edit->setText( WString::fromUTF8( distance_str ) );
    } catch( const std::exception &e )
    {}
  }// if( distance_str.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )

  try
  {
    const double distance = PhysicalUnits::stringToDistance( distance_str );
    if( distance < 0.0 )
      m_distance_edit->setText( "" );
  } catch( std::exception & )
  {
    m_distance_edit->setText( "" );
  }

  optionsChanged();
}// distanceValueChanged()

BatchGuiActShieldAnaWidget::~BatchGuiActShieldAnaWidget()
{
  wApp->doJavaScript( "removeOnDragEnterDom(['" + m_detector_file_drop->id() + "']);" );
}


FileConvertOpts::FileConvertOpts( Wt::WContainerWidget *parent )
 : BatchGuiAnaWidget( parent ),
  m_format_menu( nullptr ),
  m_overwrite_output( nullptr ),
  m_sum_for_single_output_types( nullptr )
{
  addStyleClass( "FileConvertOpts" );
  
  WContainerWidget *menuHolder = new WContainerWidget( this );
  menuHolder->addStyleClass( "SideMenuHolder" );
  
  m_format_menu = new WMenu( menuHolder );
  m_format_menu->itemSelected().connect( this, &FileConvertOpts::handleFormatChange );
  m_format_menu->addStyleClass( "SideMenu VerticalNavMenu LightNavMenu FileConvertFormatMenu" );
  
  const bool isMobile = false;
  
  Wt::WMessageResourceBundle descrip_bundle;
  if( !isMobile )
  {
    const string docroot = wApp->docRoot();
    const string bundle_file = SpecUtils::append_path(docroot, "InterSpec_resources/app_text/spectrum_file_format_descriptions" );
    descrip_bundle.use(bundle_file,true);
  }//if( !isMobile )
  
  
  auto addFormatItem = [this, &descrip_bundle, isMobile]( const char *label, SpecUtils::SaveSpectrumAsType type ){
    WMenuItem *item = m_format_menu->addItem( label );
    item->clicked().connect( boost::bind(&right_select_item, m_format_menu, item) );
    item->setData( reinterpret_cast<void *>(type) );
    
    if( !isMobile )
    {
      string description;
      if( descrip_bundle.resolveKey(label, description) )
      {
        SpecUtils::trim( description );
        
        description = Wt::Utils::htmlEncode( description, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
        
        WImage *img = new WImage( item );
        img->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
        img->setStyleClass("Wt-icon");
        img->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
        img->setFloatSide( Wt::Side::Right );
        
        HelpSystem::attachToolTipOn( img, description, true,
                                    HelpSystem::ToolTipPosition::Right,
                                    HelpSystem::ToolTipPrefOverride::InstantAlways );
      }//if( we have the description of the file )
    }//if( !isMobile )
  };//addFormatItem lambda
  
  
  addFormatItem( "N42-2012", SpecUtils::SaveSpectrumAsType::N42_2012 );
  addFormatItem( "N42-2006", SpecUtils::SaveSpectrumAsType::N42_2006 );
  addFormatItem( "CHN", SpecUtils::SaveSpectrumAsType::Chn );
  addFormatItem( "IAEA SPE", SpecUtils::SaveSpectrumAsType::SpeIaea );
  addFormatItem( "CSV", SpecUtils::SaveSpectrumAsType::Csv );
  addFormatItem( "TXT", SpecUtils::SaveSpectrumAsType::Txt );
  addFormatItem( "PCF", SpecUtils::SaveSpectrumAsType::Pcf );
  addFormatItem( "CNF", SpecUtils::SaveSpectrumAsType::Cnf );
  addFormatItem( "SPC (int)", SpecUtils::SaveSpectrumAsType::SpcBinaryInt );
  addFormatItem( "SPC (float)", SpecUtils::SaveSpectrumAsType::SpcBinaryFloat );
  addFormatItem( "SPC (ascii)", SpecUtils::SaveSpectrumAsType::SpcAscii );
  addFormatItem( "TKA", SpecUtils::SaveSpectrumAsType::Tka );
  addFormatItem( "GR-130", SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0 );
  addFormatItem( "GR-135", SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2 );
#if( SpecUtils_ENABLE_D3_CHART )
  addFormatItem( "HTML", SpecUtils::SaveSpectrumAsType::HtmlD3 );
#endif
#if( USE_QR_CODES )
  addFormatItem( "URI", SpecUtils::SaveSpectrumAsType::Uri );
#endif
  
  WContainerWidget *opts = new WContainerWidget( this );
  opts->addStyleClass( "FileConvertOptsOptions" );

  WText *txt = new WText( WString::tr("bgw-info-will-convert-file-type"), opts );
  txt->setInline( false );
  txt->addStyleClass( "FileConvertOptsTitle" );


  m_overwrite_output = new WCheckBox( WString::tr("bgw-overwrite-output-cb"), opts );
  m_overwrite_output->addStyleClass( "CbNoLineBreak" );
  m_overwrite_output->checked().connect( this, &FileConvertOpts::optionsChanged );
  
  m_sum_for_single_output_types = new WCheckBox( WString::tr("bgw-sum-multirecord-to-single"), opts );
  m_sum_for_single_output_types->addStyleClass( "CbNoLineBreak" );
  m_sum_for_single_output_types->hide();
  m_sum_for_single_output_types->checked().connect( this, &FileConvertOpts::optionsChanged );

  WText *spacer = new WText( "&nbsp;", opts );
  spacer->addStyleClass( "FileOptSpacer" );

  m_format_menu->select( 0 );
  handleFormatChange();
}//


FileConvertOpts::~FileConvertOpts()
{
  
}

void FileConvertOpts::handleFormatChange()
{
  const SpecUtils::SaveSpectrumAsType save_type = currentSaveType();
  const uint16_t max_records = max_records_in_save_type( save_type, nullptr );

  if( m_sum_for_single_output_types )
    m_sum_for_single_output_types->setHidden( max_records > 3 );
  
  optionsChanged();
}//void FileConvertOpts::handleFormatChange()


SpecUtils::SaveSpectrumAsType FileConvertOpts::currentSaveType() const
{
  const WMenuItem * const currentFormatItem = m_format_menu->currentItem();
  assert( currentFormatItem );
  
  if( currentFormatItem )
  {
    const uint64_t data = reinterpret_cast<uint64_t>( currentFormatItem->data() );
    assert( data < static_cast<int>(SpecUtils::SaveSpectrumAsType::NumTypes) );
    return SpecUtils::SaveSpectrumAsType( data );
  }//if( currentFormatItem )
  
  return SpecUtils::SaveSpectrumAsType::N42_2012;
}//SpecUtils::SaveSpectrumAsType currentSaveType() const;


void FileConvertOpts::performAnalysis( const vector<tuple<string, string, shared_ptr<const SpecMeas>>> &input_files,
                               const string &output_dir )
{
  vector<string> warnings;
  
  const SpecUtils::SaveSpectrumAsType save_type = currentSaveType();
  const bool overwrite = m_overwrite_output->isChecked();
  const bool sum_multi =  m_sum_for_single_output_types->isChecked();
  
  for( size_t index = 0; index < input_files.size(); ++index )
  {
    const string &display_name = get<0>( input_files[index] );
    const string &path_to_file = get<1>( input_files[index] );
    const shared_ptr<const SpecMeas> &spec_meas = get<2>( input_files[index] );
    assert( spec_meas );
    if( !spec_meas )
      continue;
    
    try
    {
      const uint16_t max_records = max_records_in_save_type( save_type, spec_meas );
      
      const string orig_leaf_name = SpecUtils::filename(display_name);
      const string orig_ext = SpecUtils::file_extension(orig_leaf_name);
      const string leaf_name = ((orig_ext.size() > 0) && (orig_ext.size() <= 4) && (orig_ext.size() < orig_leaf_name.size()))
      ? orig_leaf_name.substr( 0, orig_leaf_name.size() - orig_ext.size() )
      : orig_leaf_name;
      const string output_ext = SpecUtils::suggestedNameEnding(save_type);
      const string base_filename = SpecUtils::append_path( output_dir, leaf_name );
      
      const size_t num_records = spec_meas->num_measurements();
      if( num_records > max_records )
      {
        // We will write out every record into a different file.
        if( sum_multi )
        {
          const string outputname = base_filename + "." + output_ext;
          spec_meas->write_to_file( outputname, save_type );
        }else
        {
          vector<shared_ptr<const SpecUtils::Measurement> > meass = spec_meas->measurements();
          for( size_t record_num = 0; record_num < meass.size(); ++record_num )
          {
            const shared_ptr<const SpecUtils::Measurement> &m = meass[record_num];
            
            const string outputname = base_filename + "_" + std::to_string(record_num) + "." + output_ext;
            const bool is_file = SpecUtils::is_file(outputname);
            
            if( SpecUtils::is_file(outputname) && !overwrite )
            {
              warnings.push_back( "Not overwriting existing '" + outputname + "'" );
            }else
            {
              if( is_file && !SpecUtils::remove_file(outputname) )
                warnings.push_back( "Failed to delete existing file '" + outputname + "'" );
              
              spec_meas->write_to_file( outputname, {m->sample_number()}, {m->detector_name()}, save_type );
            }
          }//for( size_t record_num = 0; record_num < meass.size(); ++record_num )
        }//if( sum_multi )
      }else
      {
        const string outputname = base_filename + "." + output_ext;
        const bool is_file = SpecUtils::is_file(outputname);
        
        if( SpecUtils::is_file(outputname) && !overwrite )
        {
          warnings.push_back( "Not overwriting existing '" + outputname + "'" );
        }else
        {
          if( is_file && !SpecUtils::remove_file(outputname) )
            warnings.push_back( "Failed to delete existing file '" + outputname + "'" );
          
          spec_meas->write_to_file( outputname, save_type );
        }
      }//if( num_records > max_records ) / else
      
    }catch( std::exception &e )
    {
      warnings.push_back( "Unexpected error writing file '" + display_name + "': " + string(e.what()) );
    }
  }//for( size_t index = 0; index < input_files.size(); ++index )
  
  
  SimpleDialog *dialog = new SimpleDialog( WString::tr( "bgw-analysis-summary-title" ) );
  dialog->addStyleClass( "BatchAnalysisResultDialog" );
  
  if( warnings.empty() )
  {
    new WText( "File conversion complete.", dialog->contents() );
  }else
  {
    dialog->addStyleClass( "BatchAnalysisWarningDialog" );
    
    WContainerWidget *contents = dialog->contents();
    contents->setList( true, false );
    for( const string &warn_msg : warnings )
    {
      WContainerWidget *item = new WContainerWidget( contents );
      new WText( warn_msg, item );
    }
  }//if( warnings.empty() ) / else
  
  dialog->addButton( WString::tr( "Okay" ) );
}//void FileConvertOpts::performAnalysis(...)

  
pair<bool,Wt::WString> FileConvertOpts::canDoAnalysis() const
{
  return {true, ""};
}//pair<bool,Wt::WString> FileConvertOpts::canDoAnalysis() const


void FileConvertOpts::optionsChanged()
{
  const pair<bool,WString> ana_stat = canDoAnalysis();

  m_canDoAnalysis.emit( ana_stat.first, ana_stat.second );
}//void optionsChanged()
