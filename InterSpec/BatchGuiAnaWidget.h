#ifndef BatchGuiAnaWidget_h
#define BatchGuiAnaWidget_h
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

#include <set>
#include <tuple>
#include <string>
#include <vector>
#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/BatchPeak.h"
#include "InterSpec/BatchActivity.h"

// Forward declarations
class SpecMeas;
class BatchGuiInputFile;
class NativeFloatSpinBox;
class DetectorPeakResponse;
class FileDragUploadResource;
class BatchGuiInputSpectrumFile;

namespace Wt
{
  class WText;
  class WLabel;
  class WCheckBox;
  class WLineEdit;
}

namespace SpecUtils
{
  //enum class SpectrumType : int;
  enum class SaveSpectrumAsType : int;
}


/** Base class for all batch analysis widgets. */
class BatchGuiAnaWidget : public Wt::WContainerWidget
{
protected:
  Wt::Signal<bool> m_canDoAnalysis;

public:
  BatchGuiAnaWidget( Wt::WContainerWidget *parent );

  virtual void performAnalysis( const std::vector<std::tuple<std::string, std::string, std::shared_ptr<const SpecMeas>>> &input_files,
                                const std::string &output_dir ) = 0;

  virtual bool canDoAnalysis() const = 0;

  virtual void optionsChanged() = 0;

  Wt::Signal<bool> &canDoAnalysisSignal();
};

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
  
  // Per-file custom report variables grouped together
  Wt::WCheckBox *m_per_file_custom_report = nullptr;
  Wt::WContainerWidget *m_per_file_custom_report_container = nullptr;
  FileDragUploadResource *m_per_file_custom_report_resource = nullptr;

  // Summary custom report variables grouped together
  Wt::WCheckBox *m_summary_custom_report = nullptr;
  Wt::WContainerWidget *m_summary_custom_report_container = nullptr;
  FileDragUploadResource *m_summary_custom_report_resource = nullptr;

public:
  BatchGuiPeakFitWidget( Wt::WContainerWidget *parent = nullptr );

  ~BatchGuiPeakFitWidget();

  void handleFileUpload( Wt::WContainerWidget *dropArea, FileDragUploadResource *resource );
  void exemplarUploaded( const std::string &, const std::string & );
  void backgroundUploaded( const std::string &, const std::string & );
  void perFileCustomReportUploaded( const std::string &, const std::string & );
  void summaryCustomReportUploaded( const std::string &, const std::string & );

  void handle_remove_per_file_custom_report_upload( BatchGuiInputFile *input );
  void handle_remove_summary_custom_report_upload( BatchGuiInputFile *input );
  void handle_remove_exemplar_upload( BatchGuiInputSpectrumFile *input );

  void useCurrentBackgroundChanged();
  void useNoBackgroundChanged();

  virtual void optionsChanged() override;

  std::tuple<std::shared_ptr<SpecMeas>, std::string, std::set<int>> get_exemplar() const;
  std::tuple<std::shared_ptr<const SpecMeas>, std::string, std::set<int>> get_background() const;

  BatchPeak::BatchPeakFitOptions getPeakFitOptions() const;

  virtual void performAnalysis( const std::vector<std::tuple<std::string, std::string, std::shared_ptr<const SpecMeas>>> &input_files,
                                const std::string &output_dir ) override;

  virtual bool canDoAnalysis() const override;
};

class BatchGuiActShieldAnaWidget : public BatchGuiPeakFitWidget
{
protected:
  Wt::WContainerWidget *m_act_shield_container;

  // Activity/shielding analysis specific options
  Wt::WCheckBox *m_use_bq;
  Wt::WCheckBox *m_hard_background_sub;

  // Detector response function upload
  Wt::WContainerWidget *m_detector_input;
  Wt::WCheckBox *m_use_detector_override;
  Wt::WContainerWidget *m_detector_file_drop;
  FileDragUploadResource *m_detector_file_resource;
  std::shared_ptr<DetectorPeakResponse> m_uploaded_detector;

  // Distance override
  Wt::WCheckBox *m_override_distance;
  Wt::WContainerWidget *m_distance_input_container;
  Wt::WLabel *m_distance_label;
  Wt::WLineEdit *m_distance_edit;

  Wt::WCheckBox *m_csv_report = nullptr;

public:
  BatchGuiActShieldAnaWidget( Wt::WContainerWidget *parent = nullptr );

  BatchActivity::BatchActivityFitOptions getActivityFitOptions() const;

  virtual void performAnalysis( const std::vector<std::tuple<std::string, std::string, std::shared_ptr<const SpecMeas>>> &input_files,
                                const std::string &output_dir ) override;

  virtual bool canDoAnalysis() const override;

  virtual void optionsChanged() override;

  void useBqChanged();
  void useHardBackgroundSubChanged();
  void useDetectorOverrideChanged();

  std::shared_ptr<const DetectorPeakResponse> detector() const;

  void detectorUploaded( const std::string &, const std::string & );
  void handle_remove_detector_upload( BatchGuiInputFile *input );

  void overrideDistanceChanged();
  void distanceValueChanged();

  ~BatchGuiActShieldAnaWidget();
};




class FileConvertOpts : public BatchGuiAnaWidget
{
protected:
  Wt::WMenu *m_format_menu;
  Wt::WCheckBox *m_overwrite_output;
  Wt::WCheckBox *m_sum_for_single_output_types;
  
  void handleFormatChange();
  SpecUtils::SaveSpectrumAsType currentSaveType() const;
  
public:
  FileConvertOpts( Wt::WContainerWidget *parent = nullptr );
  ~FileConvertOpts();
  
  virtual void performAnalysis( const std::vector<std::tuple<std::string, std::string, std::shared_ptr<const SpecMeas>>> &input_files,
                               const std::string &output_dir ) override;
  
  virtual bool canDoAnalysis() const override;
  
  virtual void optionsChanged() override;
};//class FileConvertOpts


/*
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

#endif // BatchGuiAnaWidget_h 
