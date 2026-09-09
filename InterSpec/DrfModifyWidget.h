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
class DetectorPeakResponse;
class DetectorEfficiencyUncert;

namespace Wt
{
  class WText;
  class WTable;
  class WMenu;
  class WMenuItem;
  class WCheckBox;
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
 The Uncertainty tab shows a measured-point editor when the response is
 grounded to measured efficiencies, and a node σ/ρ covariance-matrix editor
 otherwise.
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

    /** The σ/ρ covariance-matrix editor's numeric shadow (display order), used when no measured
     points ground the response.  #covMatrix is the row-major N·N covariance of fractional
     efficiency error at #covEnergies (keV). */
    std::vector<double> covEnergies;
    std::vector<double> covMatrix;

    /** Energy / efficiency / stat-% / cert-% / source text, one entry per measured-point row.
     The last two are empty for an intrinsic (single-uncertainty) curve. */
    std::vector<std::array<std::string,5>> anchors;
    std::string anchorRefDistance, anchorDefaultUncert;

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

    /** Rebuild the σ/ρ covariance table from the numeric shadow; deferred to #render so a cell
     edit does not delete the WLineEdit whose `changed()` is being handled. */
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
   Uncertainty editor, and records the edit. */
  void handleModeToggle();

  /** Whether the far-field response is currently Geometry-Modeled (vs Flat Disk). */
  bool geometryModeled() const;

  /** Whether the Uncertainty tab is showing the measured-points editor (rather than the σ/ρ
   covariance matrix): true when the DRF carries measured points, or the mode is Geometry Modeled. */
  bool pointsEditorVisible() const;

  /** Shows exactly one of the two Uncertainty-tab editors per #pointsEditorVisible. */
  void updateUncertEditorVisibility();

  /** Appends one measured-point row (energy keV / absolute efficiency / statistical-uncert % /
   certificate-uncert % / source-key); the last two are only built for an absolute reference curve.
   Blank stat cells fall back to the default-uncert on apply. */
  void addAnchorRow( const float energy, const float efficiency,
                     const float fracStatUncert, const float fracCertUncert,
                     const std::string &sourceKey );
  void removeAnchorRow();

  /** Rebuilds `working`'s measured points + far-field efficiency + node covariance from the points
   editor (energy/efficiency/stat/cert/source rows, reference distance, and the default uncert %).
   No-op when the editor was not built. */
  void applyAnchorEdits( DetectorPeakResponse &working );

  /** Seeds the σ/ρ covariance shadow (#m_covEnergies / #m_covMatrix) from an existing node
   covariance, then rebuilds the table. */
  void seedCovFromUncert( const std::shared_ptr<const DetectorEfficiencyUncert> &uncert );

  /** Re-renders the covariance table from the numeric shadow: column 0 editable energies, echoed
   column headers, diagonal σ (%), editable upper-triangle ρ, disabled lower-triangle mirror. */
  void rebuildCovTable();

  void addEnergyRow();
  void removeEnergyRow();

  /** Covariance-shadow edit handlers (index into the current display order).  Each mutates the
   shadow so untouched entries stay bit-exact (σ scales its row/col to hold ρ fixed; ρ sets one
   pair), then rebuilds and records the edit. */
  void covSigmaChanged( const std::size_t i, const std::string &text );
  void covRhoChanged( const std::size_t i, const std::size_t j, const std::string &text );
  void covEnergyChanged( const std::size_t i, const std::string &text );

  /** Writes the covariance shadow onto `working` (sorted ascending); clears the uncert when empty.
   No-op when the covariance editor is not the visible one. */
  void applyCovarianceEdits( DetectorPeakResponse &working );

  /** Builds a working DRF from every tab: name/description, the visible Uncertainty editor, and
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

  /** Whether the seed DRF arrived with measured efficiency points/pairs, i.e. the points editor is
   always the right Uncertainty editor for it regardless of mode. */
  bool m_origHasPoints;

  // --- Uncertainty tab: the two swappable editors --------------------------
  /** Help text above the editors; its wording tracks which editor is shown
      (measured points vs. σ/ρ covariance matrix) - see #updateUncertEditorVisibility. */
  Wt::WText *m_uncertHelp;
  /** Holds the measured-points editor; shown when #pointsEditorVisible. */
  Wt::WContainerWidget *m_pointsEditor;
  /** Holds the σ/ρ covariance-matrix editor; shown otherwise. */
  Wt::WContainerWidget *m_covEditor;

  /** Measured-point editor: one row per reference point, plus an editable reference distance and a
   single default statistical-uncert %. */
  Wt::WTable *m_anchorTable;
  Wt::WPushButton *m_addAnchor, *m_removeAnchor;
  Wt::WLineEdit *m_anchorRefDistance;
  Wt::WLineEdit *m_anchorDefaultUncert;
  /** `cert` and `source` are null for an intrinsic (single-uncertainty) curve. */
  struct AnchorRow{ Wt::WLineEdit *energy, *eff, *stat, *cert, *source; };
  std::vector<AnchorRow> m_anchors;

  /** Whether the point editor holds ABSOLUTE efficiencies at #m_anchorRefDistance (an ANGLE-style
   reference curve) rather than INTRINSIC ones (a GADRAS Efficiency.csv).  Decides whether the
   reference-distance row and the certificate/source columns are shown, which geometry type
   #applyAnchorEdits writes back, and whether the rows are also recorded as `MeasuredDrfPoints` -
   which are read elsewhere as absolute efficiencies at a distance, so intrinsic points must not
   masquerade as them. */
  bool m_anchorIsAbsolute;

  /** σ/ρ covariance-matrix editor. */
  Wt::WTable *m_covTable;
  Wt::WPushButton *m_addEnergy, *m_removeEnergy;
  /** Authoritative numeric shadow (display order): node energies (keV) and the row-major N·N
   covariance of fractional efficiency error.  The GUI edits these; untouched entries stay
   bit-exact across a σ or ρ edit. */
  std::vector<double> m_covEnergies;
  std::vector<double> m_covMatrix;

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
