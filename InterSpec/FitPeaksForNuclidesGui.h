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

#ifndef FitPeaksForNuclidesGui_h
#define FitPeaksForNuclidesGui_h

#include "InterSpec_config.h"

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget>
#include <Wt/WFlags>

#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/SimpleDialog.h"

namespace Wt
{
  class WMenu;
  class WCheckBox;
  class WComboBox;
  class WPushButton;
  class WText;
}//namespace Wt

namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils

class D3SpectrumDisplayDiv;
class NativeFloatSpinBox;
class PeakModel;
class RelEffChart;
class RelActTxtResults;
class ReferenceLineInfo;
class DetectorPeakResponse;

class SimpleDialog;
class FitPeaksAdvancedWidget;


/** Functions in this namespace handle the GUI workflow for fitting peaks
 for nuclides shown in the ReferencePhotopeakDisplay.

 These functions must be called from the Wt application thread.
 */
namespace FitPeaksForNuclidesGui
{
  /** Creates and shows the advanced options dialog for fitting peaks.

   Returns the created dialog so the caller can track it.
   The caller is responsible for managing the dialog lifecycle.

   Must be called from the Wt application thread.

   @return Pointer to the created SimpleDialog, or nullptr if creation failed.
   */
  SimpleDialog *showAdvancedDialog();

  /** Starts the peak fitting workflow for currently displayed nuclides.

   This function performs a 3-stage async operation:
   - Stage A (GUI thread): Peak detection using AnalystChecks
   - Stage B (Background thread): Calls FitPeaksForNuclides::fit_peaks_for_nuclides()
   - Stage C (GUI thread): Results dialog with preview, color assignment, accept/cancel

   @param from_advanced_dialog If true, indicates call is from the advanced dialog
          (currently unused, but reserved for future options).

   Must be called from the Wt application thread.
   */
  void startFitSources( const bool from_advanced_dialog );

}//namespace FitPeaksForNuclidesGui


/** Advanced dialog for fitting peaks for nuclides: options, live preview, status. */
class FitPeaksAdvancedDialog : public SimpleDialog
{
public:
  FitPeaksAdvancedDialog( const Wt::WString &title );
  virtual ~FitPeaksAdvancedDialog();

  FitPeaksAdvancedWidget *widget() { return m_widget; }

protected:

private:
  void onAcceptClicked();
  void onResultUpdated();

  FitPeaksAdvancedWidget *m_widget;
  Wt::WPushButton *m_acceptBtn;
};//class FitPeaksAdvancedDialog


/** Inner content widget for FitPeaksAdvancedDialog: chart, status, warnings, options. */
class FitPeaksAdvancedWidget : public Wt::WContainerWidget
{
public:
  FitPeaksAdvancedWidget( Wt::WContainerWidget *parent = nullptr );
  virtual ~FitPeaksAdvancedWidget();

  /** Called when user clicks Accept; adds peaks to PeakModel and removes old peaks. */
  void acceptResult();

  /** True if a successful result is available and Accept is valid. */
  bool canAccept() const;

  /** Set maximum height of the chart (e.g. 75% of widget width). Call after dialog sets widget size. */
  void setChartMaxHeight( double heightPx );

  /** Emitted when result state changes (calculating finished, result received, error). */
  Wt::Signal<> &resultUpdated();

protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

private:
  enum RenderActions
  {
    UpdateCalculations = 0x01
  };

  void startComputation();
  void updateFromResult( std::shared_ptr<FitPeaksForNuclides::PeakFitResult> result,
                        std::shared_ptr<std::atomic_bool> cancel_flag,
                        size_t calc_number );
  void handleCalcError( std::shared_ptr<std::string> error_msg,
                        std::shared_ptr<std::atomic_bool> cancel_flag );

  void buildOptionsFromConfig();
  void syncConfigFromOptions();
  void scheduleOptionsUpdate();
  void onDontUseRoisChanged();
  void onExistingAsFreeChanged();
  void onDontVaryEnergyCalChanged();
  void onRelEffTypeChanged();
  FitPeaksForNuclides::PeakFitForNuclideConfig currentConfig() const;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> currentOptions() const;
  std::string sourceListTitle() const;

  Wt::WText *m_title;
  Wt::WMenu *m_upper_menu;
  D3SpectrumDisplayDiv *m_chart;
  PeakModel *m_chart_peak_model;
  RelEffChart *m_rel_eff_chart;
  RelActTxtResults *m_txt_results;
  Wt::WText *m_status;
  Wt::WContainerWidget *m_warnings_div;
  Wt::WContainerWidget *m_options_div;

  Wt::WCheckBox *m_opt_dont_use_rois;
  Wt::WCheckBox *m_opt_existing_as_free;
  Wt::WCheckBox *m_opt_dont_vary_energy_cal;
  Wt::WCheckBox *m_opt_dont_refine_energy_cal;
  Wt::WCheckBox *m_opt_fit_bkgnd_peaks;
  Wt::WCheckBox *m_opt_fit_bkgnd_dont_use;
  Wt::WCheckBox *m_opt_use_background;
  NativeFloatSpinBox *m_opt_roi_min_chi2;
  NativeFloatSpinBox *m_opt_roi_min_peak_sig;
  NativeFloatSpinBox *m_opt_obs_initial_sig;
  NativeFloatSpinBox *m_opt_obs_final_sig;
  Wt::WComboBox *m_opt_skew_type;
  Wt::WComboBox *m_opt_fwhm_form;
  Wt::WComboBox *m_opt_rel_eff_type;
  NativeFloatSpinBox *m_opt_rel_eff_order;
  Wt::WContainerWidget *m_rel_eff_order_row;

  Wt::WFlags<RenderActions> m_render_flags;
  bool m_is_calculating;
  size_t m_calc_number;
  std::shared_ptr<std::atomic_bool> m_cancel_calc;
  std::shared_ptr<FitPeaksForNuclides::PeakFitResult> m_current_result;
  Wt::Signal<> m_resultUpdated;

  std::vector<RelActCalcAuto::SrcVariant> m_sources;
  std::vector<RelActCalcAuto::NucInputInfo> m_base_nucs;
  Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> m_base_options;
  std::shared_ptr<SpecUtils::Measurement> m_fg_copy;
  std::shared_ptr<SpecUtils::Measurement> m_bg_copy;
  std::shared_ptr<const DetectorPeakResponse> m_drf;
  bool m_is_hpge;
  std::vector<ReferenceLineInfo> m_ref_lines;
  std::vector<std::shared_ptr<const PeakDef>> m_user_peaks;
  std::shared_ptr<std::vector<std::shared_ptr<const PeakDef>>> m_auto_search_peaks;
};//class FitPeaksAdvancedWidget

#endif //FitPeaksForNuclidesGui_h
