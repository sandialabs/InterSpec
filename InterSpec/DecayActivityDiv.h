#ifndef DecayActivityDiv_h
#define DecayActivityDiv_h
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
#include <string>
#include <vector>

#include <Wt/WContainerWidget.h>
#include <Wt/Core/observing_ptr.hpp>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/DecayActivityCalc.h"

namespace Wt
{
  class WText;
  class WLabel;
  class WLineEdit;
  class WCheckBox;
  class WComboBox;
  class WTabWidget;
  class WPushButton;
  class WMouseEvent;
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
  class NuclideMixture;
}//namespace SandiaDecay

class InterSpec;
struct ColorTheme;
class CsvDownloadGui;
class DecayChainChart;
class DecayCsvResource;
class DecayActivityChart;
class DecayActivityModel;
class DecaySelectNuclide;
struct NuclideSelectedInfo;
class DateLengthCalculator;

// TODO: need to implement undo/redo in DecayActivityDiv as the user changes things.

/** The contents of the "Nuclide Decay Info" tool (see `DecayWindow`): the source nuclides the user
 has entered, and tabs showing their activity over time, their decay chain, and a calculator.

 The decay math and app-url handling live in `DecayActivityCalc`.
 */
class DecayActivityDiv : public Wt::WContainerWidget
{
public:
  DecayActivityDiv( InterSpec *viewer );
  virtual ~DecayActivityDiv();

  /** Adds a source nuclide, resets the time span shown, and redraws.

   If `activityStr` is empty, or does not parse to `activity`, it is formatted from `activity`.
   */
  void addNuclide( const int z, const int a, const int iso,
                   const double activity, const bool useCurie,
                   const double age, const std::string &activityStr );

  void clearAllNuclides();

  /** Sets the time span to show, and redraws. */
  void setDecayChartTimeRange( double dt );

  /** Sets the tool state from an app-url ("interspec://decay/...") - see
   `DecayActivityCalc::fromAppUrl` for the format.  Problems with the URL are shown to the user as
   warnings.

   @param path The path of the URI (e.g., "chart" for "interspec://decay/chart?nuc=...").
   @param query_str The url-decoded query portion of the URI (e.g., "nuc=U238&...").
   */
  void handleAppUrl( std::string path, std::string query_str );

  /** Encodes the current tool state to app-url format; see `DecayActivityCalc::toAppUrl`.  The
   returned string does not include the "interspec://decay/" prefix, and has not been url-encoded.
   */
  std::string encodeStateToUrl();

private:
  friend class CsvDownloadGui;
  friend class DecayCsvResource;
  friend class DateLengthCalculator;

  /** A source nuclide, and the "pill" widget showing it. */
  struct Nuclide : public DecayActivityCalc::SourceNuclide
  {
    Wt::WContainerWidget *display = nullptr;
    Wt::WText *txt = nullptr;

    void updateTxt();
  };//struct Nuclide

  void init();
  void initCharts();
  Wt::WContainerWidget *initDisplayOptionWidgets();

  /** Adds (or, if `replaceIndex` is a valid index, replaces) a source nuclide, without updating
   anything else; callers must follow with #nuclidesChanged (or equivalent).
   */
  void addNuclideNoRefresh( DecayActivityCalc::SourceNuclide src, const int replaceIndex = -1 );
  void clearNuclidesNoRefresh();

  /** Call after the set of source nuclides changes: resets the time span shown, and redraws. */
  void nuclidesChanged();

  /** Adds the nuclide from the nuclide-select dialog, replacing the one being edited, if any. */
  void addTheNuclide( const NuclideSelectedInfo &nuc );

  void removeNuclide( Wt::WContainerWidget *frame );

  /** Asks first when more than one nuclide would be thrown away (nothing here is undoable). */
  void confirmClearAllNuclides();

  /** Shows the placeholder, and enables "remove all", according to whether any nuclide is added. */
  void updateNuclideAreaState();

  /** Shows the nuclide-select dialog to add a new nuclide. */
  void showAddNuclideDialog();

  /** Called when the nuclide-select dialog closes; stops editing any nuclide. */
  void nuclideSelectDialogDone();

  /** Shows the nuclide-select dialog to edit the nuclide the pill `w` displays; when the user
   accepts, #addTheNuclide replaces that nuclide, in place.
   */
  void sourceNuclideDoubleClicked( Wt::WContainerWidget *w );

  /** Stops editing a nuclide (if one is being edited). */
  void clearEditingNuclide();

  /** Half-life of the first source nuclide (which "half-lives" time spans refer to), or -1.0. */
  double firstNuclideHalfLife() const;

  /** Parses a time span, which may be in half-lives of the first nuclide; throws on error. */
  double parseTimeSpan( const std::string &txt ) const;

  /** Whether to show activities of `nuc` in curies: true if the first source nuclide it decays
   from (or is) was entered in curies; the user preference if it is not from any source nuclide.
   */
  bool useCurieFor( const SandiaDecay::Nuclide *nuc ) const;

  /** Returns the time span to show, setting a default into the time span field if it is invalid. */
  double timeToDisplayTill();
  void checkTimeRangeValid();
  void setTimeLimitToDisplay();
  void updateInitialMixture() const;

  /** Applies the log-y check box to the chart axis, without redrawing. */
  void applyYScale();
  void updateYScale();
  void setGridLineStatus();

  /** Redraws the activity chart from #m_nuclides; also updates the calculator results.

   @param syncCalcTimeSpan If true, the calculators time span is set to the charts.
   */
  void refreshDecayDisplay( const bool syncCalcTimeSpan );
  void refreshDecayChart();
  void addDecaySeries();
  void userSetShowSeries( int series, bool show );
  void updateYAxisRange();

  /** Updates the calculator now if its tab is showing, or else when it is next shown. */
  void updateCalculator();
  void chartTabChanged( int index );

  void colorThemeChanged( std::shared_ptr<const ColorTheme> theme );

  void createCsvDownloadGui();
  void deleteCsvDownloadGui();
  void deleteCsvDownloadGuiTriggerUpdate();

  /** Opens the "Batch Decay" tool, seeded with the currently entered nuclides. */
  void openBatchDecayTool();

  void displayMoreInfoPopup( const double time );
  void decayChartClicked( const Wt::WMouseEvent &event );
  void deleteMoreInfoDialog();
  void manageActiveDecayChainNucStyling();

  Wt::WContainerWidget *nuclideInformation( const SandiaDecay::Nuclide *nuclide ) const;
  Wt::WContainerWidget *isotopesSummary( const double time ) const;


  InterSpec *m_viewer;

  std::vector<Nuclide>         m_nuclides;
  Wt::WLineEdit               *m_displayTimeLength;
  Wt::WComboBox               *m_displayActivityUnitsCombo;
  Wt::WLabel                  *m_displayActivityUnitsLabel;
  Wt::WCheckBox               *m_logYScale;
  Wt::WCheckBox               *m_showGridLines;
  Wt::WComboBox               *m_yAxisType;

  //things for the adding of activity
  Wt::WContainerWidget        *m_parentNuclidesDiv;
  Wt::WContainerWidget        *m_nuclidesAddedDiv;

  /** Placeholder shown inside m_nuclidesAddedDiv while no nuclide has been added. */
  Wt::WText                   *m_noNuclidesTxt;
  Wt::WPushButton             *m_createNewNuclideButton;
  Wt::WPushButton             *m_clearNuclidesButton;

  Wt::Core::observing_ptr<AuxWindow> m_nuclideSelectDialog;
  DecaySelectNuclide          *m_nuclideSelect;

  /** The pill of the nuclide being edited in #m_nuclideSelectDialog; null when adding a nuclide. */
  Wt::Core::observing_ptr<Wt::WContainerWidget> m_editingNuclide;

  //Objects to actually display the chart of activities
  Wt::WTabWidget              *m_chartTabWidget;
  DecayActivityChart          *m_decayChart;
  std::shared_ptr<DecayActivityModel> m_decayModel;

  DecayChainChart             *m_decayChainChart;

  Wt::Core::observing_ptr<AuxWindow> m_moreInfoDialog;

  Wt::WContainerWidget        *m_decayLegend;

  DateLengthCalculator        *m_calc;

  /** If the calculator needs updating when it is next shown. */
  bool                         m_calcNeedsUpdate;

  // Wt4: observing_ptr auto-nulls if the dialog is destroyed (e.g. closed by the user), so the raw
  //  pointer can't dangle (matches m_nuclideSelectDialog / m_moreInfoDialog above).
  Wt::Core::observing_ptr<AuxWindow> m_csvDownloadDialog; //null if not displayed

  /** The time unit of the charts x-axis, so clicks can be mapped to a time; negative if the chart
   is empty.
   */
  double m_currentTimeUnits;

  std::unique_ptr<SandiaDecay::NuclideMixture> m_currentMixture;
};//DecayActivityDiv



//DecayActivityDiv_h
#endif
