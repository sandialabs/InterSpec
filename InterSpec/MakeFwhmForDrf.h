#ifndef MakeFwhmForDrf_h
#define MakeFwhmForDrf_h
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
#include <string>
#include <vector>
#include <utility>

#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"

class PeakDef;
class AuxWindow;
class InterSpec;
class MakeDrfChart;
class FwhmPeaksModel;
class NativeFloatSpinBox;
class RowStretchTreeView;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WComboBox;
  class WTableView;
  class WPushButton;
}//namespace Wt

/* Adds FWHM functional information to the DRF passed in; resulting
 DRF can be set to the current one for the application
 */
class MakeFwhmForDrf : public Wt::WContainerWidget
{
public:
  /** What the tool does with the DRF and the spectrum when it opens. */
  enum class InitialFit : int
  {
    /** Run the automated peak search, then fit the result. */
    SearchAndFit,

    /** Fit the users existing peaks; the automated search is left to
     #startAutomatedPeakSearch. */
    FitUserPeaks,

    /** Show the FWHM the DRF already has (or "None" when it has none), and fit only while the
     "Fit FWHM" checkbox this mode adds is checked - which it starts out being only when the DRF has
     no FWHM of its own.  For editors of an existing detector, where looking at the FWHM must not
     cost the user the one they already had.  An "Original FWHM" button appears whenever fitting is
     off and the equation has drifted from what the DRF came in with.
     */
    ShowExisting
  };//enum class InitialFit

  /** @param narrow_layout When true the peak table uses abbreviated column headers and turns off
   sorting, so the columns fit where there is little horizontal room (phones, and the
   "Modify Detector Response" dialog, where this tool shares the width with a side menu and the
   equation controls).  Wt only renders the header sort handle for sortable columns, so turning
   sorting off is what reclaims that space.
   */
  MakeFwhmForDrf( const InitialFit initial_fit,
                 InterSpec *viewer,
                 std::shared_ptr<DetectorPeakResponse> drf,
                 const bool narrow_layout );
  
  virtual ~MakeFwhmForDrf() override;

  Wt::Signal<bool> &validationChanged();
  bool isValidFwhm() const;

  /** Emitted whenever a user edit changes this tools state, but only when the tool is embedded in
   an owner that has taken over undo/redo for it (see #setOwnerHandlesUndoRedo).  Connect this, and record
   the step from #currentState / #setState, when this tool is a section of a larger dialog.
   */
  Wt::Signal<> &stateChanged();

  /** Tells this tool that its owner records undo/redo steps for it, so it should stop registering
   its own.  Its own steps resolve their target through `InterSpec::fwhmFromForegroundWindow`, which
   for an embedded tool would pop open the standalone FWHM window instead of reaching this one.
   */
  void setOwnerHandlesUndoRedo( const bool owner_handles );
  
  void setToDrf();
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();

  /** Locks the equation-type selection to `DetectorPeakResponse::kConstantPlusSqrtEnergy`
   (i.e., `FWHM = A0 + A1*sqrt(energy)`) and hides/disables the combo box that would otherwise
   let the user pick a different form.  Intended for callers (e.g. the GENIE CNF export options,
   which can only make use of this one FWHM form) that need this tool restricted to fitting only
   this equation type.
   */
  void restrictToConstantPlusSqrtEnergy();

  /** Starts the (multi-threaded, off-session-thread) automated peak search whose results seed the
   FWHM fit.  Called from the constructor when `auto_fit_peaks` is true; otherwise callers that
   embed this widget can defer the search until the user actually looks at it (see
   `DrfModifyWidget`s FWHM tab), so opening a dialog doesnt cost a full peak search.

   A no-op if a search is already running, or one has already been started.
   */
  void startAutomatedPeakSearch();


public: //Some stuff for undo/redo support
  struct TableRow
  {
    bool m_is_user_peak;
    bool m_use_for_fit;
    std::shared_ptr<const PeakDef> m_peak;
  };//struct TableRow
  
  
  struct ToolState
  {
    int m_fwhm_index;
    int m_sqrt_eqn_index;
    std::vector<MakeFwhmForDrf::TableRow> m_rows;
    std::vector<float> m_parameters, m_uncertainties;
    std::shared_ptr<DetectorPeakResponse> m_orig_drf;
    
    bool operator==( const ToolState &rhs ) const;
  };//struct ToolState

  std::shared_ptr<ToolState> currentState() const;
  void setState( std::shared_ptr<const ToolState> state );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  
  void handleTableDataChange();
  void handleFwhmEqnTypeChange();
  void handleSqrtEqnOrderChange();
  
  void coefficientManuallyChanged( const int coef_num );

  /** The "Fit FWHM" checkbox was toggled. */
  void handleFitCheckChanged();

  /** Puts the equation type, term count and coefficients back to what the DRF came in with. */
  void restoreOriginalFwhm();

  /** Whether the equation should track a fit of the peak table.  Always true where there is no
   "Fit FWHM" checkbox (i.e. everywhere but `InitialFit::ShowExisting`). */
  bool fittingEnabled() const;

  /** Enables the "Fit FWHM" checkbox only when a fit is possible (a real equation form, and a
   foreground to find peaks in), and shows the "Original FWHM" button only when fitting is off and
   the equation differs from the DRFs.  A no-op when neither was created. */
  void updateFitControls();

  void scheduleRefit();
  void doRefitWork();

  /** Pushes #m_parameters / #m_uncertainties into the coefficient edits, the equation text and the
   chart.  Split out of #doRefitWork so the tool can also display coefficients it did not fit -
   e.g. the ones the DRF arrived with.
   */
  void updateEquationDisplay();
  void setEquationToChart();
  
  std::vector<std::shared_ptr<const PeakDef>> get_user_peaks();

  /** Installs the results of #startAutomatedPeakSearch (or of the trivial no-search case) into the
   peak table, and takes the tool out of the "searching" state.

   `error_msg`, when non-empty, is the reason the search failed; the search results will then be
   empty, and the message is shown to the user.  Must be called on the session thread.
   */
  void setPeaksFromAutoSearch( std::vector<std::shared_ptr<const PeakDef>> user_peaks,
                               std::shared_ptr<std::vector<std::shared_ptr<const PeakDef>>> auto_search_peaks,
                               std::string error_msg = std::string() );
  
  void doAddUndoRedoStep();
  void scheduleUndoRedoStep();
  
protected:
  InterSpec *m_interspec;

  /** How the tool started up; `ShowExisting` is what adds the "Fit FWHM" button. */
  const InitialFit m_initial_fit;

  /** See #setOwnerHandlesUndoRedo. */
  bool m_owner_handles_undo_redo;

  bool m_currently_searching;

  /** Set once an automated peak search has been kicked off, so selecting the FWHM tab a second
   time doesnt start another one. */
  bool m_auto_search_started;

  /** Set true by the destructor, so the peak-search completion callback - which is posted from a
   worker thread and runs on the session thread - can tell "this widget is gone" from "this widget
   is alive, but I couldnt find it".  See `WidgetUtils::WidgetHandle`, and the same pattern in
   `SpecFileQueryWidget`.  Only an atomic in a shared control block crosses the thread boundary.
   */
  std::shared_ptr<std::atomic<bool>> m_widget_deleted;

  bool m_refit_scheduled;
  bool m_undo_redo_scheduled;
  
  std::shared_ptr<DetectorPeakResponse> m_orig_drf;
  std::vector<std::shared_ptr<const PeakDef>> m_user_peaks;
  std::vector<std::shared_ptr<const PeakDef>> m_auto_fit_peaks;
  
  MakeDrfChart *m_chart;
  
  Wt::WComboBox *m_fwhmEqnType;
  Wt::WComboBox *m_sqrtEqnOrder;

  /** "Fit FWHM" and "Original FWHM"; only created for `InitialFit::ShowExisting`. */
  Wt::WCheckBox *m_fitCb;
  Wt::WPushButton *m_originalBtn;

  /** What the DRF arrived with, so "Original FWHM" can put it back: the equation-type combo index
   (`kNumResolutionFnctForm` when the DRF had no FWHM), the term-count combo index (-1 when not a
   sqrt-polynomial), and the coefficients. */
  int m_original_form_index;
  int m_original_sqrt_order;
  std::vector<float> m_original_parameters;
  std::vector<NativeFloatSpinBox *> m_parEdits;
  std::vector<float> m_parameters;
  std::vector<float> m_uncertainties;
  
  Wt::WText *m_error;
  Wt::WText *m_equation;
  
  RowStretchTreeView *m_table;
  FwhmPeaksModel *m_model;
  
  Wt::Signal<bool> m_validationChanged;
  Wt::Signal<> m_stateChanged;
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;
  
  std::shared_ptr<const ToolState> m_current_state;
};//class MakeFwhmForDrf


/** Create a AuxWindow with the MakeFwhmForDrf widget as the primary content.
    Returns pointer to the created AuxWindow, but is will already be shown,
    and have the signals to delete it when closed hooked up, so you probably
    wont need the window.
 */
class MakeFwhmForDrfWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  MakeFwhmForDrf *tool();

protected:
  // Constructor is protected; use AuxWindow::make<MakeFwhmForDrfWindow>() to create.
  MakeFwhmForDrfWindow( const bool use_auto_fit_peaks_too );

  MakeFwhmForDrf *m_tool;
};//class MakeFwhmForDrfWindow

#endif //MakeFwhmForDrf_h

