#ifndef DetectionLimitDynamic_h
#define DetectionLimitDynamic_h
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
#include <string>
#include <memory>
#include <vector>

#include <Wt/WContainerWidget.h>
#include <Wt/Core/observing_ptr.hpp>

#include "InterSpec/AuxWindow.h"

class PeakModel;
class InterSpec;
class SimpleDialog;
class DetectorDisplay;
class ShieldingSelect;
class NativeFloatSpinBox;
class DetectionLimitDynamic;
class DetectorPeakResponse;
class D3SpectrumDisplayDiv;
class NuclideSourceEnterController;

namespace Wt
{
  class WText;
  class WLabel;
  class WSpinBox;
  class WLineEdit;
  class WComboBox;
  class WTabWidget;
  class WPushButton;
  class WCheckBox;
  class WTableView;
  class WSuggestionPopup;
  class WButtonGroup;
  class WRadioButton;
}//namespace Wt

namespace SandiaDecay
{
  class Nuclide;
}

namespace DetectionLimitCalc
{
  struct DynamicMdaInput;
  struct DynamicMdaResult;
  struct DynamicSearchResult;
}//namespace DetectionLimitCalc


/** AuxWindow wrapping the Dynamic MDA tool.  Constructor is protected; create
 with `AuxWindow::make<DetectionLimitDynamicWindow>(...)`. */
class DetectionLimitDynamicWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  virtual ~DetectionLimitDynamicWindow();

  DetectionLimitDynamic *tool();

protected:
  DetectionLimitDynamicWindow( InterSpec *viewer );

  DetectionLimitDynamic *m_tool;
};//class DetectionLimitDynamicWindow


/** The Dynamic MDA tool widget.  Provides:

  - A "Calculator" tab that solves the inverse problem: given any two of
    {speed, distance-of-closest-approach, source activity}, find the third
    at the user's confidence level and false-positive budget.
  - A "Search" tab that slides a fixed window across a passthrough /
    search-mode spectrum file and flags windows that exceed the per-trial
    decision threshold.
 */
class DetectionLimitDynamic : public Wt::WContainerWidget
{
public:
  DetectionLimitDynamic( InterSpec *viewer );
  virtual ~DetectionLimitDynamic();

  /** Returns the current photopeak energy (keV), or 0 if none selected. */
  float photopeakEnergy() const;
  const SandiaDecay::Nuclide *nuclide() const;

  /** Handles "interspec://dynamic-mda/<path>?<query>" URLs.
   The path indicates the tab ("calculator" or "search"). */
  void handleAppUrl( std::string uri );

  /** Encodes the current tool state as a URL fragment (without the
   "interspec://" prefix or "dynamic-mda" authority). */
  std::string encodeStateToUrl() const;

protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

  void init();

  // ----- input handlers -----
  void handleNuclideChanged();
  void handleGammaChanged();
  void handleConfidenceLevelChanged();
  void handleFpPerHourChanged();
  void handleSolveForChanged();
  void handleSpeedChanged();
  void handleDcaChanged();
  void handleActivityChanged();
  void handleWindowModeChanged();
  void handleWindowLengthChanged();
  void handleStrideModeChanged();
  void handleStrideChanged();
  void handleSnapSamplesChanged();
  void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf );
  void handleSpectrumChanged();
  void handleShieldingChanged();

  // ----- result update + search -----
  void updateCalcResult();
  void runSearch();
  void clearSearchHighlights();
  void showMoreInfoDialog();
  void handleMoreInfoClose( SimpleDialog *dialog );

  /** Returns the confidence level (β probability) selected by the GUI. */
  double currentConfidenceLevel() const;

  /** Pull the median per-sample real-time from the foreground (0 if not a
   passthrough / single-sample file).  Result is cached against the current
   foreground SpecFile pointer; cache is invalidated by handleSpectrumChanged. */
  double currentSamplePeriod() const;

  mutable const void *m_samplePeriodCacheKey = nullptr;  // shared_ptr<...>::get()
  mutable double m_samplePeriodCache = 0.0;

  /** Build a `DynamicMdaInput` from the current widget state.  Throws on
   user-input issues that prevent calculation. */
  DetectionLimitCalc::DynamicMdaInput buildCalcInput() const;

  enum RenderActions
  {
    UpdateLimit = 0x01,
    AddUndoRedoStep = 0x02
  };//enum RenderActions

  Wt::WFlags<RenderActions> m_renderFlags;

  InterSpec *m_viewer;

  // Chart at top: shows the background spectrum with ROI / side-channel
  // highlights, and a peak Gaussian when the MDA calc succeeds.  Cleared
  // of the peak (but still showing data + highlights) on calc failure.
  D3SpectrumDisplayDiv *m_spectrum;
  PeakModel            *m_peakModel;

  Wt::WTabWidget *m_tabs;

  // ----- common (parent-owned) inputs -----
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WLineEdit *m_nuclideAgeEdit;
  NuclideSourceEnterController *m_nucEnterController;
  Wt::WComboBox *m_photoPeakEnergy;
  std::vector<std::pair<double,double>> m_photoPeakEnergiesAndBr;

  Wt::WComboBox *m_confidenceLevel;
  NativeFloatSpinBox *m_fpPerHour;
  DetectorDisplay *m_detectorDisplay;

  NativeFloatSpinBox *m_lowerRoi;
  NativeFloatSpinBox *m_upperRoi;
  Wt::WSpinBox *m_numSideChannel;

  /** Optional shielding directly over the source.  When empty, transmission
   is 1.  The simple `ShieldingSelect( WSuggestionPopup* )` constructor is
   used — no fit-for checkboxes, no trace/self-attenuating sources. */
  ShieldingSelect *m_shieldingSelect;

  // ----- Tab 1: Calculator -----
  enum class SolveForId : int
  {
    Activity = 0,
    Speed = 1,
    Dca = 2
  };
  std::shared_ptr<Wt::WButtonGroup> m_solveForGroup;

  // Unit-aware text inputs.  Parsed via PhysicalUnits helpers and a small
  // local parseSpeed helper.  Each (label, input) pair lives in its own
  // flex sub-row container so the entire row can be hidden when "solve
  // for" points at that quantity.
  Wt::WContainerWidget *m_speedRow;
  Wt::WLabel           *m_speedLabel;
  Wt::WLineEdit        *m_speedInput;        // e.g. "5 m/s", "1 km/h", or bare m/s
  Wt::WContainerWidget *m_dcaRow;
  Wt::WLabel           *m_dcaLabel;
  Wt::WLineEdit        *m_dcaInput;          // e.g. "100 cm", "1 m"
  Wt::WContainerWidget *m_activityRow;
  Wt::WLabel           *m_activityLabel;
  Wt::WLineEdit        *m_activityInput;     // e.g. "1 uCi", "37 kBq"

  enum class WindowModeId : int { Auto = 0, Manual = 1 };
  std::shared_ptr<Wt::WButtonGroup> m_windowModeGroup;
  Wt::WLineEdit *m_windowLengthManual;       // user-entered time (e.g. "1 s", "100 ms")

  enum class StrideModeId : int { Seconds = 0, Samples = 1 };
  std::shared_ptr<Wt::WButtonGroup> m_strideModeGroup;
  Wt::WLineEdit *m_strideSeconds;            // user-entered time (e.g. "1 s")
  Wt::WSpinBox  *m_strideSamples;            // count
  Wt::WText *m_strideDerivedTxt;
  Wt::WCheckBox *m_snapSamplesCb;
  Wt::WText *m_trialsAlphaTxt;

  Wt::WText *m_resultTxt;
  Wt::WText *m_calcErrTxt;
  Wt::WPushButton *m_moreInfoBtn;

  Wt::Core::observing_ptr<SimpleDialog> m_moreInfoWindow;

  // ----- Tab 2: Search -----
  Wt::WText *m_searchInfoTxt;
  Wt::WPushButton *m_runSearchBtn;
  Wt::WTableView *m_hitsView;
  std::shared_ptr<Wt::WStandardItemModel> m_hitsModel;
  Wt::WText *m_searchErrTxt;

  // ----- cached state -----
  std::shared_ptr<DetectionLimitCalc::DynamicMdaInput> m_currentInput;
  std::shared_ptr<DetectionLimitCalc::DynamicMdaResult> m_currentResult;
  std::shared_ptr<DetectionLimitCalc::DynamicSearchResult> m_currentSearch;

  std::string m_stateUri;  // for undo/redo tracking

private:
  void buildCalculatorTab( Wt::WContainerWidget *container );
  void buildSearchTab( Wt::WContainerWidget *container );
  void refreshGammasForCurrentNuclide();
  void updateTrialsDisplay();
};//class DetectionLimitDynamic

#endif //DetectionLimitDynamic_h
