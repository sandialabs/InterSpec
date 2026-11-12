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

#include <map>
#include <atomic>
#include <memory>
#include <vector>
#include <optional>

#include <Wt/WSignal.h>
#include <Wt/WString.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/PeakDef.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakFitLM.h"

class PeakDef;
class InterSpec;
class PeakModel;
class SkewParamsGrid;
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


/** A spectrum, and its peaks, to fit the skew across - e.g., one of the calibration spectra of the
 "Create Detector Response Function" tool; see MakeDrf::skewFitSpectra().
 */
struct FitSkewSpectrumInput
{
  std::shared_ptr<const SpecUtils::Measurement> spectrum;

  /** Plain-text name of the spectrum, for the user to select it by. */
  Wt::WString title;

  /** Background spectra are not used to fit the skew, but their peaks are refit with the new skew on
   accept (if the user asks), so background subtraction uses the same peak shape.
   */
  bool isBackground = false;

  /** If this is the spectrum displayed as the foreground. */
  bool isDisplayedForeground = false;

  /** All the peaks in the spectrum, and the peaks of the ROIs the DRF tool is using. */
  std::vector<std::shared_ptr<const PeakDef>> allPeaks, selectedPeaks;
};//struct FitSkewSpectrumInput


/** A widget for fitting skew parameters from peaks across a spectrum - or, when given the spectra of the
 "Create Detector Response Function" tool, jointly across all of its (non-background) spectra.

 Shows a spectrum display with the peaks, controls for skew type and parameter
 values with per-parameter "fit" checkboxes, and a "Fit" button to run the fit.
 When used inside FitSkewParamsWindow, results can be accepted back to PeakFitDetPrefs.
 */
class FitSkewParamsTool : public Wt::WContainerWidget
{
public:
  /** @param drfSpectra The spectra, and peaks, of the DRF tool; if empty, the foreground, and the peaks
             detected in it, are used.
   */
  FitSkewParamsTool( InterSpec *viewer, std::vector<FitSkewSpectrumInput> drfSpectra = {} );
  virtual ~FitSkewParamsTool();

  /** If fitting across the spectra of the DRF tool. */
  bool isDrfMode() const;

  /** After #acceptResults in DRF mode, with the "refit peaks" checkbox checked: the DRF tool peaks, each
   mapped to its refit with the new skew.  Empty otherwise.
   */
  const std::map<std::shared_ptr<const PeakDef>,std::shared_ptr<const PeakDef>> &drfPeakReplacements() const;

  /** Signal emitted when fit completes or results change; dialog can update Accept button. */
  Wt::Signal<> &resultUpdated();

  /** Whether results are available for acceptance. */
  bool canAccept() const;

  /** Applies the current skew parameter values to the foreground PeakFitDetPrefs, and optionally refits
   the analysis peaks with them - or in DRF mode, the peaks of all the DRF tool spectra (see
   #drfPeakReplacements).
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
  void handleFitResults( const std::shared_ptr<PeakFitLM::MultiSpectrumFitResults> &results,
                         const std::shared_ptr<std::atomic_bool> &cancelFlag );

  /** When user edits a skew spin box, apply values to peaks and update the spectrum display. */
  void userEditedSkewValue();

  /** Sets skew type and parameter values on copies of the given peaks and returns them. */
  std::vector<std::shared_ptr<const PeakDef>> applySkewToPeaks(
    const std::vector<std::shared_ptr<const PeakDef>> &peaks ) const;

  /** Shows the spectrum, and peaks, of `m_spectra[index]`. */
  void showSpectrum( const size_t index );

  /** Sets the displayed peaks from the current spectrum's peaks, with the current skew. */
  void refreshDisplayedPeaks();

  /** In DRF mode, resets each spectrum's peaks to the characterization (selected), or all, peaks of the
   DRF tool.
   */
  void handleOnlyCharPeaksToggled();

  /** Updates the spectrum selector's text, which includes the number of peaks. */
  void updateSpectrumComboText();

  /** Enables the "Fit" button if there is something to fit, and no fit is running. */
  void updateFitButton();

  /** Updates the text describing what will be fit, from which spectra, and what happens on accept. */
  void updateDescription();

  /** Updates everything that depends on the spectra's peaks, after they change. */
  void peaksChanged();

  /** Enables, or disables (while fitting), the choice of peaks and spectra to fit. */
  void setSpectraOptionsEnabled( const bool enabled );

  /** Indices into #m_spectra of the spectra to fit: all of them, or just the displayed one. */
  std::vector<size_t> spectraToFit() const;

  /** The number of peaks in the spectra to fit. */
  size_t numPeaksToFit() const;

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

  /** Replaces the current spectrum's peaks of ROI `continuum` with `roiPeaks`, and updates the display. */
  void replaceRoiPeaks( const std::shared_ptr<const PeakContinuum> &continuum,
                        const std::vector<std::shared_ptr<const PeakDef>> &roiPeaks );


  InterSpec *m_viewer;

  /** A spectrum the skew is fit across. */
  struct SpectrumState
  {
    std::shared_ptr<const SpecUtils::Measurement> spectrum;
    Wt::WString title;

    /** Index into #m_drfInputs this is from; unused if not in DRF mode. */
    size_t drfInputIndex = 0;

    /** The peaks being worked with - edited by the user (ROI drag, continuum change, deletion), and
     replaced by the fit results.  The skew from the GUI is applied to them for display, and fitting.
     */
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };//struct SpectrumState

  /** The spectra from the DRF tool; empty if not in DRF mode. */
  std::vector<FitSkewSpectrumInput> m_drfInputs;

  /** The spectra the skew is fit across; in DRF mode the non-background DRF tool spectra, otherwise
   the foreground (if there is one).
   */
  std::vector<SpectrumState> m_spectra;

  /** Index into #m_spectra of the displayed spectrum. */
  size_t m_currentSpectrum;

  /** The energies the energy-dependent skew values are given at (the lower and upper values in the
   GUI, and PeakFitDetPrefs): the lowest and highest energies of the foreground spectrum.
   */
  double m_anchorLowerEnergy;
  double m_anchorUpperEnergy;

  // Spectrum display
  D3SpectrumDisplayDiv *m_chart;
  std::shared_ptr<PeakModel> m_peakModel;

  // Controls
  Wt::WComboBox *m_spectrumCombo;

  /** In DRF mode: only use the peaks the DRF tool is using (else, all their spectra's peaks). */
  Wt::WCheckBox *m_onlyCharPeaksCb;

  /** In DRF mode, with more than one spectrum: only fit the displayed spectrum. */
  Wt::WCheckBox *m_onlyDisplayedCb;
  Wt::WComboBox *m_skewTypeCombo;
  SkewParamsGrid *m_paramsGrid;
  Wt::WCheckBox *m_updatePeaksCb;
  Wt::WPushButton *m_fitBtn;

  /** The fit status ("Fitting...", the results, errors). */
  Wt::WText *m_statusText;

  /** What will be fit; see #updateDescription. */
  Wt::WText *m_descText;

  // Fit state
  std::optional<PeakFitLM::FitPeaksResults::SkewRelation> m_fitSkewRelation;

  /** Indices into #m_spectra of the spectra the current (or last) fit is of. */
  std::vector<size_t> m_fitSpectrumIndices;

  /** See #drfPeakReplacements. */
  std::map<std::shared_ptr<const PeakDef>,std::shared_ptr<const PeakDef>> m_drfPeakReplacements;

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
  //  See FitSkewParamsTool for `drfSpectra`.
  FitSkewParamsWindow( InterSpec *viewer, std::vector<FitSkewSpectrumInput> drfSpectra = {} );

private:
  FitSkewParamsTool *m_tool;
  Wt::WPushButton *m_acceptBtn;
};//class FitSkewParamsWindow

#endif //FitSkewParamsTool_h
