#ifndef DrfModifyWidget_h
#define DrfModifyWidget_h
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

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WFlags.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/MakeMcResponseForDrf.h"

class InterSpec;
class SwitchCheckbox;
class EccUncertOptions;
class NativeFloatSpinBox;
class DetectorPeakResponse;
class DetectorEfficiencyUncert;

namespace Wt
{
  class WText;
  class WTable;
  class WMenu;
  class WMenuItem;
  class WCheckBox;
  class WComboBox;
  class WLineEdit;
  class WTextArea;
  class WPushButton;
  class WStackedWidget;
  class WContainerWidget;
}//namespace Wt

namespace ceelo{ class DetectorResponse; struct GeometryDescriptor; }

/** A single "Modify Detector" editor consolidating the actions that used to
 crowd the Detector Response Select footer: renaming, geometry + Monte-Carlo
 characterization, FWHM fitting, and the efficiency uncertainty.

 Everything is applied to ONE working copy of the DRF when the user accepts,
 so the tabs compose (an MC characterization, an FWHM fit, a rename, and an
 uncertainty edit can all be set in one editing session).  The tool widgets
 are embedded (not their standalone windows), so the Modify dialog owns the
 single Use/Cancel decision.

 A far-field detector is either "Flat Disk" (no geometry-modeled response;
 efficiency is the solid-angle far-field model) or "Geometry Modeled" (a
 Monte-Carlo `ceelo::DetectorResponse` is attached and answers off-axis /
 near-field / uncertainty-aware queries).  The Geom & MC tab carries that
 toggle; fixed-geometry DRFs have no geometry to model and show no such tab.
 The Anchor tab edits the efficiency representation and its uncertainty together, choosing its
 editor by the curve's `EfficiencyFnctForm` (see #AnchorEditor) rather than by whether measured
 points happen to be present - so an ISOCS .ecc, a GADRAS Efficiency.csv and an ANGLE .outx all get
 the same editor, and a fixed-geometry curve is no longer rewritten as far-field on apply.
 */
class DrfModifyWidget : public Wt::WContainerWidget
{
public:
  /** The geometry (for the geometry form and the measured-anchor editor) comes from
   `drf->geometry()`. */
  DrfModifyWidget( InterSpec *viewer,
                   std::shared_ptr<const DetectorPeakResponse> drf );
  virtual ~DrfModifyWidget() override;

  /** Builds the modified DRF from all tabs and emits #updatedDrf. */
  void apply();

  /** The "Use" entry point: applies immediately, unless a confirmation is warranted first -
   detaching a geometry-modeled response when the mode was switched to Flat Disk, or regenerating a
   stale response before use.  See #apply. */
  void requestApply();

  /** Emitted (by #apply) with the modified DRF. */
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();

  /** Emitted with whether the embedded Monte-Carlo tool currently holds a
   generated response.  The window uses it to hold "Use" closed for a DRF that
   arrived with no efficiency of its own (a geometry-only import), where a
   response is the only thing that can give it one.
   */
  Wt::Signal<bool> &mcResponseAvailable();

  /** Whether the DRF being modified needs a Monte-Carlo response before it can
   be used at all, i.e. it came in with no efficiency curve. */
  bool needsMcResponse() const;

  /** The DRF this editor was seeded with, un-modified - what re-creating this dialog takes. */
  std::shared_ptr<const DetectorPeakResponse> originalDrf() const;


  /** Everything the user can change in this dialog, for undo/redo.

   Snapshot/diff at render time, the same shape `RelActAutoGui` and `DetectionLimitSimple` use:
   each edit handler flags `AddUndoRedoStep` and schedules a render; the render compares the new
   snapshot against the previous one and records a step only when they actually differ.
   */
  struct ToolState
  {
    std::string name, description;
    int tabIndex = 0;

    /** Whether the far-field response is Geometry-Modeled (a Monte-Carlo response is/should be
     attached) rather than Flat Disk.  Always false for fixed-geometry DRFs (no Geom & MC tab). */
    bool geometryModeled = false;

    /** The coefficient σ/ρ matrix editor's numeric shadow: the row-major M·M covariance of the
     exp-of-log-power-series coefficients.  Empty for the other two efficiency forms. */
    std::vector<double> coefCovMatrix;

    /** The exp-of-log-power-series coefficient text, one entry per term. */
    std::vector<std::string> coefficients;

    /** The efficiency formula text (`kFunctialEfficienyForm` only). */
    std::string formula;

    /** Energy / efficiency / stat-% / cert-% / source text, one entry per point row.  Cells for
     columns this curve does not show are empty. */
    std::vector<std::array<std::string,5>> anchors;
    std::string anchorRefDistance, anchorDefaultUncert;

    /** The energy-correlation length the correlated column is combined with; `EccUncertOptions`
     `effectiveCorrLength()` semantics, so one field round-trips all three modes.  -1 when the
     active editor has no correlation control. */
    double anchorCorrLength = -1.0;

    /** Selected energy units for the coefficient/formula editors (`PhysicalUnits`). */
    float efficiencyEnergyUnits = 1.0f;

    MakeMcResponseForDrf::State mc;
    std::shared_ptr<const MakeFwhmForDrf::ToolState> fwhm;

    bool operator==( const ToolState &rhs ) const;
    bool operator!=( const ToolState &rhs ) const{ return !((*this) == rhs); }
  };//struct ToolState

  std::shared_ptr<ToolState> currentState() const;

  /** Restores a #currentState snapshot; records no undo/redo step of its own. */
  void setState( const std::shared_ptr<const ToolState> &state );

protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;

  /** Kicks off the FWHM tools automated peak search the first time its tab is opened - it is a
   multi-threaded fit of the whole spectrum, so it is not worth paying for every time this dialog is
   opened, only when the user actually looks at the tab.  (Fitting the FWHM to those peaks is a
   separate, opt-in step - see `MakeFwhmForDrf::InitialFit::ShowExisting`.)  Also records the
   undo/redo step for the section change. */
  void handleTabSelected( Wt::WMenuItem *item );

  enum RenderActions
  {
    AddUndoRedoStep = 0x01,

    /** Rebuild the coefficient σ/ρ covariance table from the numeric shadow; deferred to #render
     so a cell edit does not delete the WLineEdit whose `changed()` is being handled. */
    RebuildCovTable = 0x02
  };//enum RenderActions

  /** Flags this dialogs state as user-edited, so the next render records an undo/redo step. */
  void scheduleUndoRedoStep();

  /** A user edit: records an undo/redo step, marks any generated response stale, and refreshes the
   footer "Generate Response" button.  Every edit handler routes through here. */
  void markEdited();

  /** Re-baselines #m_currentState, and (when flagged) records the step from the old baseline. */
  void doAddUndoRedoStep( const bool add_step );

  /** Flat Disk / Geometry Modeled toggled: greys the Geom & MC tool in Flat Disk, swaps the visible
   Anchor editor, and records the edit. */
  void handleModeToggle();

  /** Whether the far-field response is currently Geometry-Modeled (vs Flat Disk). */
  bool geometryModeled() const;

  /** Which of the Anchor tab's three editors is showing.  Chosen by the efficiency curve's
   `EfficiencyFnctForm`, except that Geometry-Modeled always uses #Points - those rows are the
   Monte-Carlo grounding anchors. */
  enum class AnchorEditor
  {
    /** `kEnergyEfficiencyPairs`: the energy/efficiency/uncertainty table. */
    Points,
    /** `kExpOfLogPowerSeries`: coefficient boxes plus their M·M covariance. */
    Coefficients,
    /** `kFunctialEfficienyForm`: the formula text, a flat default uncertainty, and an optional
     per-energy uncertainty table. */
    Formula
  };//enum class AnchorEditor

  AnchorEditor activeAnchorEditor() const;

  /** Shows exactly one of the three Anchor-tab editors per #activeAnchorEditor, and words the help
   text to match the columns actually rendered. */
  void updateAnchorEditorVisibility();

  /** Appends one point row.  Every column is built; the Efficiency one is hidden by CSS for a
   formula curve (see #updateAnchorEditorVisibility), whose efficiency comes from the formula, so
   its rows describe uncertainty only.  The Source column only exists per #m_anchorHasSourceCol.
   Blank stat cells fall back to the default-uncert on apply. */
  void addAnchorRow( const float energy, const float efficiency,
                     const float fracStatUncert, const float fracCertUncert,
                     const std::string &sourceKey );
  void removeAnchorRow();

  /** Whether every Anchor-tab widget still holds the value it was seeded with.  When it does, the
   apply paths leave the DRF's curve and uncertainty exactly as they were, rather than rebuilding
   an equivalent-but-not-identical one - which matters for a covariance the correlated+diagonal
   model cannot reproduce (a source-blocked `MeasuredDrfPoints` matrix, or one restored from a URL
   with no component split). */
  bool anchorsMatchSeed() const;

  /** Reads the point rows and the correlation control into the uncertainty the Anchor tab means:
   `fromCorrelatedPlusDiagonal` over the rows, or - when there are no rows - a single flat node at
   the default uncert.  Returns nullptr when there is nothing to express. */
  std::shared_ptr<DetectorEfficiencyUncert> buildUncertFromRows();

  /** Rebuilds `working`'s efficiency curve + uncertainty (and, for an absolute reference curve,
   its `MeasuredDrfPoints`) from the point editor.  No-op when that editor was not built, or when
   nothing was edited. */
  void applyAnchorEdits( DetectorPeakResponse &working );

  /** Rebuilds `working`'s exp-of-log-power-series curve and coefficient covariance from the
   coefficient editor. */
  void applyCoefficientEdits( DetectorPeakResponse &working );

  /** Rebuilds `working`'s formula curve and uncertainty from the formula editor. */
  void applyFormulaEdits( DetectorPeakResponse &working );

  /** Seeds the coefficient covariance shadow (#m_coefCovMatrix) from an existing uncertainty's
   `coefficientCovariance()`, then rebuilds the table. */
  void seedCoefCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert );

  /** Re-renders the coefficient covariance table from the numeric shadow: column 0 the static
   `A0..An` labels, diagonal σ (%), editable upper-triangle ρ, disabled lower-triangle mirror. */
  void rebuildCovTable();

  /** Adds/removes an exp-of-log-power-series term, resizing both the coefficient boxes and the
   covariance shadow; the surviving block of the covariance stays bit-exact. */
  void addCoefficient( const float value = 0.0f );
  void removeCoefficient();

  /** Covariance-shadow edit handlers (coefficient indices).  Each mutates the shadow so untouched
   entries stay bit-exact (σ scales its row/col to hold ρ fixed; ρ sets one pair), then rebuilds and
   records the edit. */
  void covSigmaChanged( const std::size_t i, const std::string &text );
  void covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text );

  /** Validates the formula text with a trial `setIntrinsicEfficiencyFormula`, flagging the field
   `Wt-invalid` when it will not parse.  Returns whether it is usable. */
  bool validateFormula();

  /** The `PhysicalUnits` energy unit the equation/formula is written in, per #m_effEnergyUnits.
   Distinct from #m_anchorEnergyUnits - see that member. */
  float equationEnergyUnits() const;

  /** Builds a working DRF from every tab: name/description, the visible Anchor editor, and
   FWHM.  When `includeMcResponse`, attaches the generated (or existing) Monte-Carlo response in
   Geometry-Modeled mode and detaches it in Flat Disk; otherwise (a regeneration seed) always
   detaches, so the manual points/covariance drive grounding. */
  std::shared_ptr<DetectorPeakResponse> buildWorkingDrf( const bool includeMcResponse );

  /** Footer "Generate Response": re-seeds the MC tool with the live edits (detached) and starts a
   generation.  Geometry-Modeled only. */
  void handleGenerateResponse();

  /** Shows/enables the footer generate button per mode, geometry readiness, and pending edits. */
  void updateGenerateButton();

  /** MC tool finished a generation: clears the stale flag (the following #userChanged is not a user
   edit), and, when a regenerate-then-use is pending, applies. */
  void handleResponseGenerated( std::shared_ptr<ceelo::DetectorResponse> response );

  InterSpec *m_interspec;
  std::shared_ptr<const DetectorPeakResponse> m_orig;
  std::shared_ptr<const ceelo::GeometryDescriptor> m_geometry;

  Wt::WMenu *m_tabMenu;
  Wt::WStackedWidget *m_tabStack;

  Wt::WLineEdit *m_name;
  Wt::WTextArea *m_description;

  MakeMcResponseForDrf *m_mcTool;   //null for a fixed-geometry DRF (no Geom & MC tab)
  MakeFwhmForDrf *m_fwhmTool;

  /** The FWHM tabs menu item, so selecting it can kick off the peak search. */
  Wt::WMenuItem *m_fwhmTabItem;

  /** Flat Disk / Geometry Modeled toggle (checked == Geometry Modeled), above the Geom & MC tool;
   null for a fixed-geometry DRF. */
  SwitchCheckbox *m_modeToggle;

  /** Current mode: whether a Monte-Carlo response is (to be) attached.  Always false for a
   fixed-geometry DRF. */
  bool m_geometryModeled;

  /** Whether the seed DRF arrived with efficiency points/pairs. */
  bool m_origHasPoints;

  // --- Anchor tab: the three swappable editors -----------------------------
  /** Help text above the editors; its wording tracks which editor, and which columns, are shown -
      see #updateAnchorEditorVisibility. */
  Wt::WText *m_anchorHelp;
  /** Holds the point-table editor; shown for #AnchorEditor::Points and #AnchorEditor::Formula
   (where its Efficiency column is hidden, the efficiency coming from the formula). */
  Wt::WContainerWidget *m_pointsEditor;
  /** Holds the coefficient boxes and their covariance matrix; shown for
   #AnchorEditor::Coefficients. */
  Wt::WContainerWidget *m_coefEditor;
  /** Holds the efficiency formula text; shown for #AnchorEditor::Formula. */
  Wt::WContainerWidget *m_formulaEditor;

  /** Point editor: one row per point, plus an editable reference distance and a single default
   statistical-uncert %. */
  Wt::WTable *m_anchorTable;
  /** Wraps #m_anchorTable so it can scroll and centre inside the panel. */
  Wt::WContainerWidget *m_anchorTableWrap;
  Wt::WPushButton *m_addAnchor, *m_removeAnchor;
  Wt::WLineEdit *m_anchorRefDistance;
  Wt::WLineEdit *m_anchorDefaultUncert;
  /** `source` is null unless #m_anchorHasSourceCol. */
  struct AnchorRow{ Wt::WLineEdit *energy, *eff, *stat, *cert, *source; };
  std::vector<AnchorRow> m_anchors;

  /** How the correlated column is correlated across energy; null when the points carry a source
   column (see #m_anchorHasSourceCol), whose certificate uncertainty is blocked per source instead. */
  EccUncertOptions *m_uncertOptions;

  /** Whether the point table has a per-source certificate/Source column, i.e. its correlated
   uncertainty is blocked by source key rather than governed by #m_uncertOptions.  Only an absolute
   reference curve has one.  Deliberately distinct from #m_anchorIsAbsolute: conflating the two is
   what hid the correlated column from every fixed-geometry DRF. */
  bool m_anchorHasSourceCol;

  /** The `PhysicalUnits` energy unit the point table's Energy column is in.

   For a pairs curve this is the curve's own unit (a GADRAS CSV may use MeV), so the numbers shown
   are the numbers stored.  For a formula curve the rows are covariance nodes rather than curve
   points, and `DetectorEfficiencyUncert` documents node energies as keV regardless of the curve -
   so it is keV there.  Fixed at construction; deliberately NOT tied to #m_effEnergyUnits, which
   only says what units the equation/formula is written in. */
  float m_anchorEnergyUnits;

  /** Whether the point editor holds ABSOLUTE efficiencies at #m_anchorRefDistance (an ANGLE-style
   reference curve) rather than INTRINSIC ones (a GADRAS Efficiency.csv).  Decides whether the
   reference-distance row and the certificate/source columns are shown, which geometry type
   #applyAnchorEdits writes back, and whether the rows are also recorded as `MeasuredDrfPoints` -
   which are read elsewhere as absolute efficiencies at a distance, so intrinsic points must not
   masquerade as them. */
  bool m_anchorIsAbsolute;

  /** Coefficient editor (`kExpOfLogPowerSeries`): one spin box per term, their M·M σ/ρ covariance
   table, and the energy-units selector the equation is written in. */
  std::vector<NativeFloatSpinBox *> m_coefEdits;
  Wt::WContainerWidget *m_coefParams;
  Wt::WTable *m_covTable;
  Wt::WPushButton *m_addCoef, *m_removeCoef;
  /** Authoritative numeric shadow: the row-major M·M coefficient covariance.  The GUI edits this;
   untouched entries stay bit-exact across a σ or ρ edit. */
  std::vector<double> m_coefCovMatrix;

  /** Formula editor (`kFunctialEfficienyForm`). */
  Wt::WTextArea *m_formulaText;

  /** keV / MeV for the coefficient and formula editors, and the row holding it (hidden for a
   points curve, whose energies are in the column itself). */
  Wt::WComboBox *m_effEnergyUnits;
  Wt::WContainerWidget *m_effUnitsRow;

  /** The Anchor-tab state this dialog opened with, so #anchorsMatchSeed can tell "the user changed
   nothing" from "the user re-typed the same numbers". */
  std::shared_ptr<const ToolState> m_seedState;

  /** Footer "Generate Response" button (Geometry-Modeled only); regenerates the MC response from
   the live edits. */
  Wt::WPushButton *m_generateBtn;

  /** Whether an edit has landed since the last generated response, i.e. the response is stale. */
  bool m_changedSinceGenerate;

  /** Set true just before a regenerate-then-use starts, so #handleResponseGenerated applies once the
   response lands. */
  bool m_applyAfterGenerate;

  /** One-shot: the #userChanged that a generation emits right after #responseGenerated is not a user
   edit, so it must not re-mark the fresh response stale. */
  bool m_suppressNextEditMark;

  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;

  Wt::WFlags<RenderActions> m_renderFlags;

  /** The state the next undo step will restore to; re-baselined on every render. */
  std::shared_ptr<const ToolState> m_currentState;

  /** Set while #setState is restoring, so its widget edits dont look like user edits. */
  bool m_restoringState;
};//class DrfModifyWidget


/** AuxWindow hosting a DrfModifyWidget with Use/Cancel; create via
 AuxWindow::make<DrfModifyWindow>( viewer, drf ).
 */
class DrfModifyWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  DrfModifyWidget *tool();

protected:
  DrfModifyWindow( InterSpec *viewer,
                   std::shared_ptr<const DetectorPeakResponse> drf );

  DrfModifyWidget *m_tool;
};//class DrfModifyWindow

#endif //DrfModifyWidget_h
