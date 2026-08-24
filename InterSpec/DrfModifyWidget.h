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

#include <memory>
#include <vector>

#include <Wt/WContainerWidget.h>

#include "InterSpec/AuxWindow.h"

class InterSpec;
class MakeFwhmForDrf;
class MakeMcResponseForDrf;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WTable;
  class WMenu;
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
  DrfModifyWidget( InterSpec *viewer,
                   std::shared_ptr<const DetectorPeakResponse> drf,
                   std::shared_ptr<const ceelo::GeometryDescriptor> geometry = nullptr );
  virtual ~DrfModifyWidget() override;

  /** Builds the modified DRF from all tabs and emits #updatedDrf. */
  void apply();

  /** Emitted (by #apply) with the modified DRF. */
  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> &updatedDrf();

protected:
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

  Wt::Signal<std::shared_ptr<DetectorPeakResponse>> m_updatedDrf;
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
                   std::shared_ptr<const DetectorPeakResponse> drf,
                   std::shared_ptr<const ceelo::GeometryDescriptor> geometry = nullptr );

  DrfModifyWidget *m_tool;
};//class DrfModifyWindow

#endif //DrfModifyWidget_h
