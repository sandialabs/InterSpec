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
#include <vector>
#include <string>
#include <functional>

#include <Wt/WModelIndex.h>
#include <Wt/WContainerWidget.h>

//Forward declarations
class PeakModel;
class InterSpec;
class SourceFitModel;

namespace Wt
{
  class WText;
  class WColor;
  class WPushButton;
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
}//namespace SandiaDecay

namespace ShieldingSourceFitCalc
{
  struct BackDecayResults;
}//namespace ShieldingSourceFitCalc


/** A custom per-nuclide "card" display that replaces the default Wt MVC table view of
 #SourceFitModel in the Activity/Shielding Fit tool.

 #SourceFitModel remains the single source of truth: this widget *reads* through the model's
 data()/flags()/accessors and *writes* through setData() (so undo/redo, cross-isotope age
 propagation, and uncertainty resets all keep working), and it rebuilds/refreshes its cards in
 response to the model's row- and data-changed signals.

 Each card shows the nuclide name + half-life, an editable activity field (with the fit
 uncertainty shown in parentheses when fit), the derived mass and - when the "Back decay" option is
 enabled - the activity at the back-decay reference time, an age field (only when age is relevant),
 Fit checkboxes (only when applicable), and a "via progeny" row listing progeny nuclides whose
 emissions are attributed to this parent.  While back-decay is in use, a footer under the cards gives the
 reference time, with a link that opens the decay tool from it.
 */
class SourceFitNuclideDisplay : public Wt::WContainerWidget
{
public:
  /** @param nestPeaks When true (phone "merged" layout), each card becomes an accordion that nests
   the list of *its* peaks (toggle a peak's "use in fit" right inside the card).  When false (the
   default, desktop), the cards behave exactly as before and the peak table stays separate. */
  SourceFitNuclideDisplay( SourceFitModel *model,
                           std::shared_ptr<PeakModel> peakModel,
                           InterSpec *interspec,
                           const bool nestPeaks = false );
  virtual ~SourceFitNuclideDisplay();

  /** Reconciles the displayed cards with the model (adds/removes/reorders cards and refreshes
   each one's contents).  Cheap to call; safe to call any time on the session thread. */
  void reconcileCards();

  /** A compact, fit-ready summary of one source nuclide, used to populate the phone bottom
   status/fit bar (the scrollable activity strip) and the results sheet. */
  struct NucActivitySummary
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;
    std::string symbol;       //!< nuclide id, e.g. "Cs137"
    std::string colorCss;     //!< the nuclide's peak-line color (empty -> use a default)
    std::string activity;     //!< activity formatted in the display units
    std::string activityUncert; //!< uncertainty (empty when not fit / unavailable)
    std::string age;          //!< age (empty when aging is not relevant)
    bool fitActivity = false; //!< whether the activity is being fit for

    /** Activity (± uncertainty) at the back-decay reference time, or a short reason it could not be
     determined; empty when back-decay is not enabled. */
    std::string backDecayActivity;
  };//struct NucActivitySummary

  /** Returns one summary per nuclide currently in the model, in model order. */
  std::vector<NucActivitySummary> activitySummaries() const;

  /** Function giving the activities at the back-decay reference time, for the sources currently in
   the model; it should return nullptr when back-decay is not enabled.  Called once per
   #reconcileCards (each model change), so manually entered activities update without a fit. */
  typedef std::function<std::shared_ptr<const ShieldingSourceFitCalc::BackDecayResults>()> BackDecayProvider;

  /** Sets the back-decay provider, and refreshes the cards (unless clearing it, with nullptr). */
  void setBackDecayProvider( BackDecayProvider provider );

  /** One source as the decay tool should start it: `activity` of the nuclide, with `age`. */
  struct DecayToolSource
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;
    double activity = 0.0;
    double age = 0.0;
  };

  /** The sources at the earlier of the back-decay reference time and the measurement start, as the
   back-decay calculation has them; nuclides whose activity then is not determined, or not
   positive, are left out. */
  static std::vector<DecayToolSource> decayToolStart( const ShieldingSourceFitCalc::BackDecayResults &bd );

protected:
  class NuclideCard;  //defined in the .cpp

  /** Connected to the model's dataChanged/rowsInserted/rowsRemoved/layoutChanged signals. */
  void handleModelChanged();

  /** Connected (only in nestPeaks mode) to the *peak* model's change signals, so toggling a
   peak that does not move a SourceFitModel row still refreshes the nested peak lists. */
  void handlePeakModelChanged();

  /** Looks up the peak-line color for a nuclide (from any of its peaks).  Empty if none. */
  std::string peakColorForNuclide( const SandiaDecay::Nuclide *nuc ) const;

  /** Shows the back-decay footer (reference time and "Show Decay" link) only while back-decay is in
   use, i.e., #m_backDecay is valid; called from #reconcileCards. */
  void updateBackDecayFooter();

  /** Opens the decay tool with #decayToolStart, its chart spanning three times the back-decay time. */
  void showDecayFromRefTime();

  SourceFitModel *m_model;
  std::shared_ptr<PeakModel> m_peakModel;
  InterSpec *m_interspec;

  /** When true, cards are accordions that nest their peak rows (phone layout). */
  bool m_nestPeaks;

  /** See #setBackDecayProvider. */
  BackDecayProvider m_backDecayProvider;

  /** The provider's result, from the most recent #reconcileCards; read by the cards. */
  std::shared_ptr<const ShieldingSourceFitCalc::BackDecayResults> m_backDecay;

  /** Scrolling container that holds the per-nuclide cards. */
  Wt::WContainerWidget *m_cardsHolder;

  /** Footer under the cards, only shown while back-decay is in use; see #updateBackDecayFooter.
   It is where the reference time is shown: the option's own status text only shows errors (except
   on phones, whose options are on a separate tab).  #m_backDecayFootTxt is "T_ref: <date>", and
   #m_backDecayFootSpan, after it on the same line when there is room, how long before or after the
   measurement that is. */
  Wt::WContainerWidget *m_backDecayFoot;
  Wt::WText *m_backDecayFootTxt;
  Wt::WText *m_backDecayFootSpan;
  Wt::WPushButton *m_showDecayBtn;

  std::map<const SandiaDecay::Nuclide *, NuclideCard *> m_cards;

  /** Re-entrancy guard: true while applying programmatic widget updates, so the edit handlers
   ignore any spurious change signals and reconcile cannot recurse. */
  bool m_refreshing;

  friend class NuclideCard;
};//class SourceFitNuclideDisplay

#endif //SourceFitNuclideDisplay_h
