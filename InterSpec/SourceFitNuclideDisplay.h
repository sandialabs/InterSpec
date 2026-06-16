#ifndef SourceFitNuclideDisplay_h
#define SourceFitNuclideDisplay_h
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
#include <memory>

#include <Wt/WModelIndex.h>
#include <Wt/WContainerWidget.h>

//Forward declarations
class PeakModel;
class InterSpec;
class SourceFitModel;

namespace SandiaDecay
{
  struct Nuclide;
}//namespace SandiaDecay


/** A custom per-nuclide "card" display that replaces the default Wt MVC table view of
 #SourceFitModel in the Activity/Shielding Fit tool.

 #SourceFitModel remains the single source of truth: this widget *reads* through the model's
 data()/flags()/accessors and *writes* through setData() (so undo/redo, cross-isotope age
 propagation, and uncertainty resets all keep working), and it rebuilds/refreshes its cards in
 response to the model's row- and data-changed signals.

 Each card shows the nuclide name + half-life, an editable activity field (with the fit
 uncertainty shown in parentheses when fit), the derived mass and activity-at-T0, an age field
 (only when age is relevant), Fit checkboxes (only when applicable), and a "via progeny" row
 listing daughter nuclides whose emissions are attributed to this parent.
 */
class SourceFitNuclideDisplay : public Wt::WContainerWidget
{
public:
  SourceFitNuclideDisplay( SourceFitModel *model,
                           std::shared_ptr<PeakModel> peakModel,
                           InterSpec *interspec );
  virtual ~SourceFitNuclideDisplay();

  /** Reconciles the displayed cards with the model (adds/removes/reorders cards and refreshes
   each one's contents).  Cheap to call; safe to call any time on the session thread. */
  void reconcileCards();

protected:
  class NuclideCard;  //defined in the .cpp

  /** Connected to the model's dataChanged/rowsInserted/rowsRemoved/layoutChanged signals. */
  void handleModelChanged();

  SourceFitModel *m_model;
  std::shared_ptr<PeakModel> m_peakModel;
  InterSpec *m_interspec;

  /** Scrolling container that holds the per-nuclide cards. */
  Wt::WContainerWidget *m_cardsHolder;

  std::map<const SandiaDecay::Nuclide *, NuclideCard *> m_cards;

  /** Re-entrancy guard: true while applying programmatic widget updates, so the edit handlers
   ignore any spurious change signals and reconcile cannot recurse. */
  bool m_refreshing;

  friend class NuclideCard;
};//class SourceFitNuclideDisplay

#endif //SourceFitNuclideDisplay_h
