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
class DetectorPeakResponse;

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
}//namespace Wt

namespace ceelo{ struct GeometryDescriptor; }

/** A single "Modify Detector" editor consolidating the actions that used to
 crowd the Detector Response Select footer: renaming, geometry + Monte-Carlo
 characterization, FWHM fitting, and a baseline (energy-range-band) efficiency
 uncertainty.

 Everything is applied to ONE working copy of the DRF when the user accepts,
 so the tabs compose (an MC characterization, an FWHM fit, a rename, and a
 baseline uncertainty can all be set in one editing session).  The tool
 widgets are embedded (not their standalone windows), so the Modify dialog
 owns the single Use/Cancel decision.
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

    /** Lower / upper / fractional-uncert text, one entry per uncertainty-band row. */
    std::vector<std::array<std::string,3>> bands;

    /** Energy / efficiency / uncert text, one entry per measured-anchor row. */
    std::vector<std::array<std::string,3>> anchors;
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
    AddUndoRedoStep = 0x01
  };//enum RenderActions

  /** Flags this dialogs state as user-edited, so the next render records an undo/redo step. */
  void scheduleUndoRedoStep();

  /** Re-baselines #m_currentState, and (when flagged) records the step from the old baseline. */
  void doAddUndoRedoStep( const bool add_step );

  void addBandRow( const float lowerEnergy, const float upperEnergy,
                   const float fracUncert );
  void removeBandRow();

  /** Appends one measured-anchor row (energy keV / absolute efficiency /
   fractional-uncert %); blank cells fall back to the default-uncert on apply. */
  void addAnchorRow( const float energy, const float efficiency,
                     const float fracUncert );
  void removeAnchorRow();

  /** Rebuilds `working`'s measured points + far-field absolute efficiency from
   the anchor editor (energy/efficiency/uncert rows, reference distance, and the
   default uncert %).  No-op when the anchor editor was not built. */
  void applyAnchorEdits( DetectorPeakResponse &working );

  InterSpec *m_interspec;
  std::shared_ptr<const DetectorPeakResponse> m_orig;
  std::shared_ptr<const ceelo::GeometryDescriptor> m_geometry;

  Wt::WMenu *m_tabMenu;
  Wt::WStackedWidget *m_tabStack;

  Wt::WLineEdit *m_name;
  Wt::WTextArea *m_description;

  MakeMcResponseForDrf *m_mcTool;
  MakeFwhmForDrf *m_fwhmTool;

  /** The FWHM tabs menu item, so selecting it can kick off the peak search. */
  Wt::WMenuItem *m_fwhmTabItem;

  /** Baseline uncertainty band editor: one row per energy band. */
  Wt::WTable *m_bandTable;
  Wt::WPushButton *m_addBand, *m_removeBand;
  struct BandRow{ Wt::WLineEdit *lower, *upper, *frac; };
  std::vector<BandRow> m_bands;

  /** Measured-curve anchor editor (only built when seeded with an ANGLE-style
   geometry + measured points): one row per reference point, plus an editable
   reference distance and a single default fractional-uncert %. */
  Wt::WTable *m_anchorTable;
  Wt::WPushButton *m_addAnchor, *m_removeAnchor;
  Wt::WLineEdit *m_anchorRefDistance;
  Wt::WLineEdit *m_anchorDefaultUncert;
  struct AnchorRow{ Wt::WLineEdit *energy, *eff, *uncert; };
  std::vector<AnchorRow> m_anchors;

  /** Whether the anchor editor holds ABSOLUTE efficiencies at #m_anchorRefDistance (an ANGLE-style
   reference curve) rather than INTRINSIC ones (a GADRAS Efficiency.csv).  Decides whether the
   reference-distance row is shown, which geometry type #applyAnchorEdits writes back, and whether
   the rows are also recorded as `MeasuredDrfPoints` - which are read elsewhere as absolute
   efficiencies at a distance, so intrinsic points must not masquerade as them. */
  bool m_anchorIsAbsolute;

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
