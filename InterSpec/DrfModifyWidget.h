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
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WFlags.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/MakeFwhmForDrf.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/MakeMcResponseForDrf.h"

class DrfChart;
class InterSpec;
class SwitchCheckbox;
class EccUncertOptions;
class MeasuredDrfPoints;
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

namespace ceelo{ class DetectorResponse; struct GeometryDescriptor; enum class ResponseProfile : uint8_t; }

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
 editor with `DrfModifyCalc::editorForDrf` - purely from what the DRF carries, never from the
 Flat-Disk / Geometry-Modeled toggle, so flipping that toggle cannot send an apply down a different
 editor than the one the user typed into.

 Two rules keep the emitted DRF self-consistent, and both were learned the hard way:

  - **The rows are seeded from whatever the apply writes back.**  A Create-DRF detector's measured
    points are ABSOLUTE efficiencies at their own source distances while its curve is intrinsic, so
    an apply that took those rows as the curve made the detector ~200x too insensitive.  Now the
    measured-point editors write the points and RE-FIT the equation from them
    (`MakeDrfCalc::refitEfficiencyFromPoints`), which is the only operation that keeps the curve, its
    coefficient covariance, the node covariance and the points all describing one detector.
  - **A response's staleness is derived from content, not tracked with a flag** - see
    #responseStale.  While attached, a `ceelo::DetectorResponse` answers every efficiency and
    covariance query, so an edit it does not reflect is an edit the program ignores.

 `DrfModifyCalc` holds the apply logic itself (it is unit tested; this class is the Wt plumbing).
 */
class DrfModifyWidget : public Wt::WContainerWidget
{
public:
  DrfModifyWidget( InterSpec *viewer,
                   std::shared_ptr<const DetectorPeakResponse> drf );
  virtual ~DrfModifyWidget() override;

  /** Builds the modified DRF from all tabs and emits #updatedDrf.  Assumes #requestApply has already
   validated the edits. */
  void apply();

  /** The "Use" entry point: validates the edits, and applies unless something must be settled first -
   an edit that could not be applied (reported, and nothing is emitted), a geometry-modeled response
   being detached, or a response that no longer reflects the edits.  See #apply.  */
  void requestApply();

  /** Emitted (by #apply) with the modified DRF. */
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();

  /** Emitted with whether the embedded Monte-Carlo tool currently holds a
   generated response.  The window uses it to hold "Use" closed for a DRF that
   arrived with no efficiency of its own (a geometry-only import), where a
   response is the only thing that can give it one.
   */
  Wt::Signal<bool> &mcResponseAvailable();
  
  /** Emitted when a Monte-Carlo run starts or ends, with whether one is now in flight.
   
   For an owner with footer buttons that must not act during a run - "Use" would offer to generate
   a response while the run that would produce it is already going, and then silently do nothing,
   since only one generation may be in flight at a time.
   */
  Wt::Signal<bool> &generatingChanged();
  
  /** Whether a Monte-Carlo run is in flight right now. */
  bool isGenerating() const;

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
    /** One point row: its cells as typed, plus which of the DRF's measured points it describes.

     The seed index has to travel with the state: it is how an edited row keeps its point's
     provenance (peak area, live time, file name, distance uncertainty), and an undo that dropped it
     would silently discard that provenance on the next apply.
     */
    struct RowState
    {
      std::array<std::string,6> cells;
      int seedIndex = -1;

      bool operator==( const RowState &rhs ) const
      {
        return (cells == rhs.cells) && (seedIndex == rhs.seedIndex);
      }
    };//struct RowState

    /** `hashValue()` of the DRF the dialog was opened on.  An undo step outlives the dialog (it is
     re-resolved against whatever tool exists when it runs), so without this a step recorded for one
     detector would be replayed onto a different one - transplanting its cell text and, worse, its
     `RowState::seedIndex` provenance.  #setState ignores a state whose detector is not the one on
     screen. */
    uint64_t drfHash = 0;

    std::string name, description;
    int tabIndex = 0;

    /** Whether the far-field response is Geometry-Modeled (a Monte-Carlo response is/should be
     attached) rather than Flat Disk.  Always false for fixed-geometry DRFs (no Geom & MC tab). */
    bool geometryModeled = false;

    /** The coefficient editor's numeric shadow, kept as 1-sigma per coefficient plus a correlation
     matrix rather than a covariance: typing a 0 sigma then typing it back must not destroy that
     coefficient's correlations, and you cannot recover them from a zeroed row of a covariance. */
    std::vector<double> coefSigmas;
    std::vector<double> coefRho;     //row-major N*N, unit diagonal

    /** The exp-of-log-power-series coefficient text, one entry per term. */
    std::vector<std::string> coefficients;

    /** The efficiency formula text (`kFunctialEfficienyForm` only). */
    std::string formula;

    /** Energy / efficiency / stat-% / cert-% / source / distance-cm text, one entry per point row.
     Cells for columns this curve does not show are empty. */
    std::vector<RowState> anchors;
    std::string anchorRefDistance;

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
  
  /** Selects the "Geom & MC" tab, for an opener that already knows the user came here to
   characterize a detector - e.g. one just imported from a QR code with a shape but no response.
   
   No-op when the tab does not exist (a fixed-geometry DRF has no geometry to model).  Selects by
   item rather than index on purpose: the indices shift with that tab's presence.
   */
  void showGeometryTab();
  
  /** Selects the "Geom & MC" tab, switches to Geometry Modeled, sets the build method, and starts
   a run - so a caller that already knows what the user asked for (e.g. the app-URL import's
   detector-modeling choice) lands them on a characterization already under way, with this tool's
   own progress, ETA and cancel controls, rather than on a form they have to drive themselves.
   
   Returns whether a generation actually started; false when the DRF has no geometry to model, the
   geometry is incomplete, or a run is already in flight.
   */
  bool startMcCharacterization( const MakeMcResponseForDrf::Method method,
                                const ceelo::ResponseProfile profile,
                                const MakeMcResponseForDrf::Precision precision );

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

    /** Rebuild the coefficient σ/ρ table from the numeric shadow; deferred to #render so a cell edit
     does not delete the WLineEdit whose `changed()` is being handled. */
    RebuildCovTable = 0x02,

    /** Re-read the uncertainty this DRF would report, and the response-staleness wording. */
    RefreshSummary = 0x04
  };//enum RenderActions

  /** Flags this dialogs state as user-edited, so the next render records an undo/redo step. */
  void scheduleUndoRedoStep();

  /** A user edit: records an undo/redo step, refreshes the Anchor tab's "what this detector
   reports" line and the General tab's chart and summary, and refreshes the footer "Generate
   Response" button.  Every edit handler routes through here.

   It does NOT set a staleness flag: whether the attached response still describes the edits is
   derived from their content by #responseStale. */
  void markEdited();

  /** A user edit the response does not depend on: the name, the description, or the FWHM (which is
   a separate detector property the Monte Carlo never sees).  Everything #markEdited does except
   marking the response stale - offering to regenerate over one of these is just confusing. */
  void markEditedNoRegen();


  /** Re-baselines #m_currentState, and (when flagged) records the step from the old baseline. */
  void doAddUndoRedoStep( const bool add_step );

  /** Flat Disk / Geometry Modeled toggled: greys the Geom & MC tool in Flat Disk and records the
   edit.  Deliberately does NOT change which Anchor editor is showing - that follows the efficiency
   representation, not how the detector answers off-axis questions. */
  void handleModeToggle();

  /** Whether the far-field response is currently Geometry-Modeled (vs Flat Disk). */
  bool geometryModeled() const;

  /** Whether the point table carries the Source / Distance columns, i.e. its rows are the raw
   measured points rather than the curve's own numbers. */
  bool anchorHasSourceCols() const;

  /** Shows exactly one of the three Anchor-tab editors for #m_editor, and words the help text to
   match. */
  void updateAnchorEditorVisibility();

  /** Appends one point row.  Every column is built; the Efficiency one is hidden by CSS for a
   formula curve (see #updateAnchorEditorVisibility), whose efficiency comes from the formula, so
   its rows describe uncertainty only.  The Source/Distance columns only exist per
   #anchorHasSourceCols.  `seedIndex` is the measured point this row describes, or -1.
   A blank uncertainty cell means zero - there is no default standing in for one. */
  void addAnchorRow( const float energy, const float efficiency,
                     const float fracStatUncert, const float fracCertUncert,
                     const std::string &sourceKey, const float distance = -1.0f,
                     const int seedIndex = -1 );
  void removeAnchorRow();

  /** Per-editor "did the user change anything here" checks, against the state the dialog opened
   with.  One flag for all three editors is what let a mode flip discard the edits of whichever
   editor was not applied; and an untouched editor must leave the DRF bit-identical, which matters
   for a covariance the correlated+diagonal model cannot reproduce (a source-blocked
   `MeasuredDrfPoints` matrix, or one restored from a URL with no component split).  */
  bool pointsEdited() const;
  bool coefficientsEdited() const;
  bool coefCovarianceEdited() const;
  bool formulaEdited() const;

  /** Reads the point table into `DrfModifyCalc::PointRow`s, reporting a malformed cell rather than
   skipping it.  Returns whether every row parsed. */
  bool collectPointRows( std::vector<DrfModifyCalc::PointRow> &rows,
                         std::vector<DrfModifyCalc::Problem> &problems ) const;

  /** The Anchor-tab settings that are not per-row. */
  DrfModifyCalc::AnchorOptions anchorOptions() const;

  /** Applies whichever Anchor editor is showing, when it was edited.  Returns false only when
   something the user typed could not be applied; `problems` always says what. */
  bool applyAnchorTab( DetectorPeakResponse &working,
                       std::vector<DrfModifyCalc::Problem> &problems );

  /** Shows `problems` to the user; blocking ones as errors, the rest as information. */
  void showProblems( const std::vector<DrfModifyCalc::Problem> &problems );

  /** Seeds the coefficient σ/ρ shadow from an existing uncertainty's `coefficientCovariance()`,
   falling back to the legacy per-coefficient uncertainties (flagged as a placeholder - see
   #m_coefCovIsPlaceholder), then rebuilds the table. */
  void seedCoefCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert );

  /** Re-renders the coefficient σ/ρ table from the numeric shadow: column 0 the static
   `A0..An` labels, diagonal σ, editable upper-triangle ρ, disabled lower-triangle mirror.  Also
   updates the "these correlations are impossible" warning. */
  void rebuildCovTable();

  /** Adds/removes an exp-of-log-power-series term, resizing both the coefficient boxes and the σ/ρ
   shadow; the surviving block stays bit-exact. */
  void addCoefficient( const float value = 0.0f );
  void removeCoefficient();

  /** Covariance-shadow edit handlers (coefficient indices).  Each writes one entry of the σ/ρ
   shadow, so every other entry stays bit-exact, then rebuilds and records the edit. */
  void covSigmaChanged( const std::size_t i, const std::string &text );
  void covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text );

  /** Validates the formula text with a trial `setIntrinsicEfficiencyFormula`, flagging the field
   `Wt-invalid` when it will not parse.  Returns whether it is usable. */
  bool validateFormula();

  /** The `PhysicalUnits` energy unit the equation/formula is written in, per #m_effEnergyUnits.
   Distinct from #m_anchorEnergyUnits - see that member. */
  float equationEnergyUnits() const;

  /** Builds a working DRF from every tab: name/description, the visible Anchor editor, the geometry,
   and FWHM.  When `includeMcResponse`, attaches the generated (or existing) Monte-Carlo response in
   Geometry-Modeled mode and detaches it in Flat Disk; otherwise (a regeneration seed) always
   detaches, so the manual points/covariance drive grounding.  In Geometry-Modeled mode with no
   response at all, the geometry typed into the form is recorded on the DRF instead.

   `problems` collects everything that could not be applied; nothing is shown to the user from here
   (#requestApply decides that), so this is safe to call for a seed or a preview. */
  std::shared_ptr<DetectorPeakResponse> buildWorkingDrf( const bool includeMcResponse,
                                      std::vector<DrfModifyCalc::Problem> &problems );

  /** Whether the response that would be attached no longer describes the current edits, derived by
   comparing `DrfModifyCalc::seedFingerprint` against the fingerprint of the seed the attached
   response was built from.

   Content, not a flag: a flag has to be cleared by whoever regenerates, and the paths that forgot to
   (an automatic rebuild from a pre-edit seed; a redo that restored the edits but not the flag) left
   an edited curve behind a response that still answered every query.  It also means the Flat-Disk /
   Geometry-Modeled toggle cannot make a current response look stale - the mode is not content.  */
  bool responseStale();

  /** Rebuilds the General tab's chart and summary table from a quiet preview of the working DRF
   (see #buildWorkingDrf) - what "Use" would produce right now.  Cheap enough for every visit. */
  void refreshGeneralTab();

  /** An edit landed: refresh the General tab now if it is showing, else on its next visit. */
  void markGeneralStale();

  /** The rows of #m_infoTable, in order; #InfoRow::NumInfoRow is the row count. */
  enum InfoRow
  {
    InfoEfficiency, InfoGeometry, InfoSupport, InfoUncert,
    InfoFwhm, InfoTotalEff, InfoRange, InfoDiameter,
    NumInfoRow
  };//enum InfoRow

  /** Builds #m_infoTable's rows once, with their fixed labels.  The values are only ever set with
   `setText` afterwards (see #fillInfoTable) - rebuilding these widgets on each refresh left Wt
   wiring client-side handlers to elements it had already replaced. */
  void buildInfoTable( Wt::WContainerWidget *parent );

  /** Sets each row of #m_infoTable from what `drf` carries: efficiency source, geometry, location
   support, uncertainty, FWHM, total efficiency, energy range, diameter. */
  void fillInfoTable( const std::shared_ptr<const DetectorPeakResponse> &drf );

  /** Footer "Generate Response": re-seeds the MC tool with the live edits (through the seed
   provider) and starts a generation.  Geometry-Modeled only.  Returns whether a generation actually
   started - a caller that wants to apply the result afterwards must not arm itself for a run that
   never began. */
  bool handleGenerateResponse();

  /** Shows the footer generate button whenever the mode is Geometry Modeled, enables it per
   geometry readiness and #responseStale, and puts the reason it is blocked (if any) in
   #m_generateHint. */
  void updateGenerateButton();

  /** MC tool finished a generation: records what it was generated from, and, when a
   regenerate-then-use is pending, applies. */
  void handleResponseGenerated( std::shared_ptr<ceelo::DetectorResponse> response );

  /** Switches to Flat Disk (detaching the response) and applies - the honest alternative to using a
   response that ignores the edits. */
  void detachResponseAndApply();

  /** Shows or hides the note saying the correlation control currently has nothing to act on. */
  void updateCorrelationNote();

  /** Updates the read-only "what this detector reports" line: the efficiency uncertainty the
   analysis would propagate, split into the part this detector's data supports and the part that is
   an ad hoc model envelope.  Also refreshes the response-staleness note. */
  void refreshUncertSummary();

  InterSpec *m_interspec;
  std::shared_ptr<const DetectorPeakResponse> m_orig;

  /** The DRF's raw measured points as the dialog opened, which rows carry provenance from. */
  std::shared_ptr<const MeasuredDrfPoints> m_seedPoints;

  /** Which Anchor editor this DRF gets, from `DrfModifyCalc::editorForDrf`.  Fixed at construction:
   it describes how the DRF represents its efficiency, which this dialog does not change. */
  DrfModifyCalc::AnchorEditor m_editor;

  Wt::WMenu *m_tabMenu;
  Wt::WStackedWidget *m_tabStack;

  Wt::WLineEdit *m_name;
  Wt::WTextArea *m_description;

  /** General tab: a live chart of the efficiency + FWHM, and a summary of what the detector
   carries; both rebuilt from a quiet preview of the working DRF (see #refreshGeneralTab). */
  DrfChart *m_generalChart;
  Wt::WTable *m_infoTable;

  /** The value cell of each #InfoRow; only their text changes after construction. */
  std::array<Wt::WText *,NumInfoRow> m_infoValues;
  Wt::WMenuItem *m_generalTabItem;

  /** An edit landed while another tab was showing; the General tab refreshes on its next visit. */
  bool m_generalStale;

  MakeMcResponseForDrf *m_mcTool;   //null for a fixed-geometry DRF (no Geom & MC tab)
  MakeFwhmForDrf *m_fwhmTool;

  /** The FWHM tabs menu item, so selecting it can kick off the peak search. */
  Wt::WMenuItem *m_fwhmTabItem;

  /** The Geom & MC tabs menu item, so selecting it can time the Monte-Carlo estimate; null for a
   fixed-geometry DRF. */
  Wt::WMenuItem *m_geomTabItem;

  /** Flat Disk / Geometry Modeled toggle (checked == Geometry Modeled), above the Geom & MC tool;
   null for a fixed-geometry DRF. */
  SwitchCheckbox *m_modeToggle;

  /** Current mode: whether a Monte-Carlo response is (to be) attached.  Always false for a
   fixed-geometry DRF. */
  bool m_geometryModeled;

  // --- Anchor tab: the three swappable editors -----------------------------
  /** Help text above the editors; its wording tracks which editor is shown -
      see #updateAnchorEditorVisibility. */
  Wt::WText *m_anchorHelp;

  /** Shown when a geometry-modeled response is attached: that response, not the numbers below,
   answers every efficiency and uncertainty query, and is rebuilt from them on "Use". */
  Wt::WText *m_responseNote;

  /** The uncertainty this detector reports, and how much of it is model envelope - see
   #refreshUncertSummary. */
  Wt::WText *m_uncertSummary;

  /** Holds the point-table editor; shown for every editor except
   #DrfModifyCalc::AnchorEditor::Coefficients (for a formula curve its Efficiency column is hidden,
   the efficiency coming from the formula). */
  Wt::WContainerWidget *m_pointsEditor;
  /** Holds the coefficient boxes and their σ/ρ matrix. */
  Wt::WContainerWidget *m_coefEditor;
  /** Holds the efficiency formula text. */
  Wt::WContainerWidget *m_formulaEditor;

  /** Point editor: one row per point, plus an editable reference distance. */
  Wt::WTable *m_anchorTable;
  /** Wraps #m_anchorTable so it can scroll and centre inside the panel. */
  Wt::WContainerWidget *m_anchorTableWrap;
  Wt::WPushButton *m_addAnchor, *m_removeAnchor;
  Wt::WLineEdit *m_anchorRefDistance;
  /** `source` and `dist` are null unless #anchorHasSourceCols; `seedIndex` is the measured point
   this row describes (-1 for a row the user added). */
  struct AnchorRow
  {
    Wt::WLineEdit *energy, *eff, *stat, *cert, *source, *dist;
    int seedIndex = -1;
  };//struct AnchorRow
  std::vector<AnchorRow> m_anchors;

  /** How the correlated column is correlated across energy; null when the rows are measured points,
   whose certificate uncertainty is blocked per source instead. */
  EccUncertOptions *m_uncertOptions;

  /** Shown while every "Corr. %" cell is blank, i.e. while #m_uncertOptions governs nothing. */
  Wt::WText *m_corrInertNote;

  /** The `PhysicalUnits` energy unit the point table's Energy column is in.

   keV except for a pairs curve, where it is the curve's own unit (a GADRAS CSV may use MeV) so the
   numbers shown are the numbers stored.  Measured points and covariance nodes are keV by contract
   (`MeasuredEffPoint::energy`, `DetectorEfficiencyUncert`) whatever units the curve's equation uses.
   Fixed at construction; deliberately NOT tied to #m_effEnergyUnits, which only says what units the
   equation/formula is written in. */
  float m_anchorEnergyUnits;

  /** Coefficient editor (`kExpOfLogPowerSeries`): one spin box per term, their sigma/rho table, and
   the energy-units selector the equation is written in. */
  std::vector<NativeFloatSpinBox *> m_coefEdits;
  Wt::WContainerWidget *m_coefParams;
  Wt::WTable *m_covTable;
  Wt::WPushButton *m_addCoef, *m_removeCoef;

  /** Authoritative numeric shadow of the coefficient covariance, as 1-sigma per coefficient and a
   row-major N*N correlation matrix (unit diagonal).  The GUI edits these; untouched entries stay
   bit-exact, and a 0 sigma does not erase correlations. */
  std::vector<double> m_coefSigmas;
  std::vector<double> m_coefRho;

  /** Shown when the σ/ρ the user entered cannot describe any real set of errors. */
  Wt::WText *m_covWarning;

  /** Shown while #m_coefCovIsPlaceholder: the matrix on screen is not one this DRF carries. */
  Wt::WText *m_covPlaceholderNote;

  /** Whether the user has edited the σ/ρ table in this session.  Tracked rather than inferred from a
   value diff, because adding or removing a term resizes the shadow and would look like an edit. */
  bool m_coefCovTouched;

  /** True when the σ/ρ shown was manufactured from the legacy per-coefficient uncertainties rather
   than read from a stored covariance.  Such a matrix assumes the coefficients are independent, which
   for a log-power-series fit they emphatically are not, so it is offered as a starting point and
   only written to the DRF if the user actually edits it. */
  bool m_coefCovIsPlaceholder;

  /** Formula editor (`kFunctialEfficienyForm`). */
  Wt::WTextArea *m_formulaText;

  /** keV / MeV for the coefficient and formula editors, and the row holding it (hidden for a
   points curve, whose energies are in the column itself). */
  Wt::WComboBox *m_effEnergyUnits;
  Wt::WContainerWidget *m_effUnitsRow;

  /** The Anchor-tab state this dialog opened with, so the per-editor edited-checks can tell "the
   user changed nothing" from "the user re-typed the same numbers". */
  std::shared_ptr<const ToolState> m_seedState;

  /** Footer "Generate Response" button (Geometry-Modeled only); regenerates the MC response from
   the live edits. */
  Wt::WPushButton *m_generateBtn;

  /** `DrfModifyCalc::seedFingerprint` of the seed the currently held response was built from - the
   DRF as it arrived, or the seed of the last generation.  See #responseStale. */
  std::size_t m_generatedFromFingerprint;

  /** The fingerprint of the seed handed to the generation now running, promoted to
   #m_generatedFromFingerprint when it lands. */
  std::size_t m_pendingSeedFingerprint;

  /** The export tip in the footer, hidden while #m_generateHint has something to say. */
  Wt::WText *m_exportNote;

  /** Why a response cannot be generated right now (an incomplete geometry, or one still guessed
   from the diameter) - shown beside the disabled generate button; hidden when it can. */
  Wt::WText *m_generateHint;

  /** Set true just before a regenerate-then-use starts, so #handleResponseGenerated applies once the
   response lands. */
  /** The MakeMcResponseForDrf generation id that "Generate & use" armed, or -1.  Deliberately an
   id rather than a bool: a generation that fails, is cancelled or is superseded never emits
   `responseGenerated`, so a bool would stay armed and make the *next* generation - including the
   transfer one that rebuilds itself whenever the geometry changes - silently apply the detector
   and close the dialog. */
  int m_applyAfterGenerationId;

  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;

  Wt::Signal<bool> m_generatingChanged;

  /** What #m_generatingChanged last reported, so it only fires on a transition. */
  bool m_wasGenerating;

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
