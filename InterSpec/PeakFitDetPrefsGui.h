#ifndef PeakFitDetPrefsGui_h
#define PeakFitDetPrefsGui_h
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

#include <Wt/WContainerWidget>

class InterSpec;
class NativeFloatSpinBox;
struct PeakFitDetPrefs;

namespace Wt
{
  class WText;
  class WLabel;
  class WComboBox;
  class WCheckBox;
}//namespace Wt


/** A collapsible panel for viewing and editing the PeakFitDetPrefs of the
 current foreground spectrum.

 In compact mode (used in PeakInfoDisplay), the panel has a narrow collapsed
 strip (~25px) that expands to show full controls.  In normal mode (used in
 MakeDrf), the controls are always visible without a collapse mechanism.
 */
class PeakFitDetPrefsGui : public Wt::WContainerWidget
{
public:
  /** Constructor.
   @param viewer The main InterSpec instance.
   @param compactMode If true, show collapsed/expanded strip behavior; if false,
          show full controls directly (for use in MakeDrf).
   @param parent Optional parent widget.
   */
  PeakFitDetPrefsGui( InterSpec *viewer, const bool compactMode,
                       Wt::WContainerWidget *parent = nullptr );

  virtual ~PeakFitDetPrefsGui();

  /** Called externally when the PeakFitDetPrefs on the foreground SpecMeas
   changes (e.g., from auto-loading or DRF change). Updates the GUI to reflect
   the new values.
   */
  void handlePrefsChanged();

  /** Returns the current preferences as set in the GUI.
   Returns nullptr if no foreground is loaded.
   */
  std::shared_ptr<const PeakFitDetPrefs> currentPrefs() const;

protected:
  void init();

  /** Toggles between collapsed and expanded state (compact mode only). */
  void toggleExpanded();

  /** Rebuilds the skew parameter rows for the currently selected skew type. */
  void updateSkewParamRows();

  /** Reads GUI values, builds a PeakFitDetPrefs, and pushes it to the
   foreground SpecMeas. Registers undo/redo step. Called when user changes
   any control.
   */
  void userChangedValue();

  /** Updates all controls from the current foreground SpecMeas prefs.
   Sets m_programmaticUpdate to suppress undo/redo registration.
   */
  void updateFromSpecMeas();

  /** Enables or disables all editing controls. */
  void setControlsEnabled( const bool enabled );

  /** Opens the FitSkewParamsWindow dialog. */
  void showFitSkewDialog();


  InterSpec *m_viewer;
  bool m_compactMode;

  /** Guard to suppress undo/redo when updating controls programmatically. */
  bool m_programmaticUpdate;

  // Compact-mode collapse/expand widgets
  Wt::WContainerWidget *m_collapsedDiv;
  Wt::WContainerWidget *m_expandedDiv;

  // Controls
  Wt::WComboBox *m_detTypeCombo;
  Wt::WComboBox *m_skewTypeCombo;

  /** Container for the dynamically generated skew parameter rows. */
  Wt::WContainerWidget *m_skewParamsDiv;

  /** Lower-energy spin boxes for skew params (index 0..3). May be nullptr if
   the parameter is not used by the current skew type.
   */
  NativeFloatSpinBox *m_lowerSkewSpin[4];

  /** Upper-energy spin boxes for skew params (index 0..3). */
  NativeFloatSpinBox *m_upperSkewSpin[4];

  /** Checkbox for ROI-independent skew. */
  Wt::WCheckBox *m_roiIndepCb;

  /** Link to open the fit skew parameters dialog. */
  Wt::WText *m_fitSkewLink;

  /** Read-only label showing the LoadingSource. */
  Wt::WText *m_sourceLabel;
};//class PeakFitDetPrefsGui

#endif //PeakFitDetPrefsGui_h
