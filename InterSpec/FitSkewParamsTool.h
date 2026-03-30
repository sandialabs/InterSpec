#ifndef FitSkewParamsTool_h
#define FitSkewParamsTool_h
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

#include <atomic>
#include <memory>
#include <vector>
#include <optional>

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/PeakDef.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakFitLM.h"

class PeakDef;
class InterSpec;
class PeakModel;
class NativeFloatSpinBox;
class D3SpectrumDisplayDiv;

namespace Wt
{
  class WText;
  class WComboBox;
  class WCheckBox;
  class WPopupMenu;
  class WPushButton;
}//namespace Wt

namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils


/** A widget for fitting skew parameters from detected peaks across the spectrum.

 Shows a spectrum display with detected peaks, controls for skew type and parameter
 values with per-parameter "fit" checkboxes, and a "Fit" button to run the fit.
 When used inside FitSkewParamsWindow, results can be accepted back to PeakFitDetPrefs.
 */
class FitSkewParamsTool : public Wt::WContainerWidget
{
public:
  FitSkewParamsTool( InterSpec *viewer );
  virtual ~FitSkewParamsTool();

  /** Signal emitted when fit completes or results change; dialog can update Accept button. */
  Wt::Signal<> &resultUpdated();

  /** Whether results are available for acceptance. */
  bool canAccept() const;

  /** Applies the current skew parameter values to the foreground PeakFitDetPrefs,
   and optionally replaces analysis peaks with the fit peaks.
   */
  void acceptResults();

  /** Returns the "update analysis peaks on accept" checkbox, so the dialog can
   place it in the footer if desired.
   */
  Wt::WCheckBox *updatePeaksCb();

protected:
  void initWidgets();

  /** Rebuilds skew parameter rows for the current skew type, with "fit" checkboxes. */
  void updateSkewParamRows();

  /** Runs the fit in a background thread. */
  void doFit();

  /** Called on the GUI thread when the background fit completes. */
  void handleFitResults( const std::shared_ptr<PeakFitLM::FitPeaksResults> &results,
                         const std::shared_ptr<std::atomic_bool> &cancelFlag );

  /** When user edits a skew spin box, apply values to peaks and update the spectrum display. */
  void userEditedSkewValue();

  /** Sets skew type and parameter values on copies of the given peaks and returns them. */
  std::vector<std::shared_ptr<const PeakDef>> applySkewToPeaks(
    const std::vector<std::shared_ptr<const PeakDef>> &peaks ) const;

  /** Handles the user dragging an existing ROI edge on the spectrum chart. */
  void handleRoiDrag( double new_lower, double new_upper, double roi_px,
                      double orig_lower, std::string spec_type, bool is_final );

  /** Handles right-click on the spectrum chart — shows context menu for peak operations. */
  void handleRightClick( double energy, double counts, int pageX, int pageY,
                         std::string ref_line );

  /** Deletes the peak nearest to the given energy from the local model. */
  void deletePeakNearEnergy( double energy );

  /** Changes the continuum type for the peak nearest to the given energy. */
  void changeContinuumTypeNearEnergy( double energy, int continuum_type );

  /** Handles shift+drag on chart: removes peaks whose mean is in [x0, x1]. */
  void handleShiftKeyDrag( double x0, double x1 );


  InterSpec *m_viewer;

  // Spectrum display
  D3SpectrumDisplayDiv *m_chart;
  PeakModel *m_peakModel;
  std::shared_ptr<const SpecUtils::Measurement> m_spectrum;

  // Controls
  Wt::WComboBox *m_skewTypeCombo;
  Wt::WContainerWidget *m_paramsDiv;
  NativeFloatSpinBox *m_lowerSpin[4];
  NativeFloatSpinBox *m_upperSpin[4];
  Wt::WCheckBox *m_fitCb[4];
  Wt::WCheckBox *m_updatePeaksCb;
  Wt::WPushButton *m_fitBtn;
  Wt::WText *m_statusText;

  // Fit state
  std::vector<std::shared_ptr<const PeakDef>> m_detectedPeaks;
  std::vector<std::shared_ptr<const PeakDef>> m_fitPeaks;
  std::optional<PeakFitLM::FitPeaksResults::SkewRelation> m_fitSkewRelation;

  // Right-click context menu
  Wt::WPopupMenu *m_rightClickMenu;
  double m_rightClickEnergy;

  // Background calculation state
  bool m_isCalculating;
  std::shared_ptr<std::atomic_bool> m_cancelCalc;
  Wt::Signal<> m_resultUpdated;
};//class FitSkewParamsTool


/** AuxWindow wrapper around FitSkewParamsTool.

 Lifetime is managed by InterSpec (via showFitSkewParamsWindow / closeFitSkewParamsWindow /
 acceptFitSkewParamsWindow), so undo/redo and spectrum changes can properly manage the window.
 */
class FitSkewParamsWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  virtual ~FitSkewParamsWindow();

  FitSkewParamsTool *tool();

protected:
  // Constructor is protected; use AuxWindow::make<FitSkewParamsWindow>() to create.
  FitSkewParamsWindow( InterSpec *viewer );

private:
  FitSkewParamsTool *m_tool;
  Wt::WPushButton *m_acceptBtn;
};//class FitSkewParamsWindow

#endif //FitSkewParamsTool_h
