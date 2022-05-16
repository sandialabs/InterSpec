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

#include <map>
#include <string>
#include <vector>
#include <memory>

#include <Wt/WEvent>
#include <Wt/WConfig.h>
#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"


namespace Wt
{
  class WText;
  class WSlider;
  class WAnchor;
  class WLineEdit;
  class WCheckBox;
  class WComboBox;
  class WTabWidget;
  class WDatePicker;
  class WPushButton;
  class WDoubleSpinBox;
  class WStandardItemModel;

  namespace Chart
  {
    class WCartesianChart;
  }//namespace Chart
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
  class NuclideMixture;
  class SandiaDecayDataBase;
}//namespace SandiaDecay

class InterSpec;
class CsvDownloadGui;
class PeakCsvResource;
class DecayChainChart;
class DecayActivityChart;
class DecayActivityModel;
class DecaySelectNuclide;
class ChartToImageResource;
struct NuclideSelectedInfo;
class DateLengthCalculator;


class DecayActivityDiv : public Wt::WContainerWidget
{
public:
  DecayActivityDiv( InterSpec *viewer, Wt::WContainerWidget *parent = NULL );
  virtual ~DecayActivityDiv();
  
  void addNuclide( const int z, const int a, const int iso,
                  const double activity, const bool useCurrie,
                  const double age );
  
  void clearAllNuclides();
  
  void setDecayChartTimeRange( double dt );
  
  void colorThemeChanged();
  
  /** Handles receiving a "deep-link" url starting with "interspec://decay/...".
   
   Example URIs:
   - "interspec://decay/chain?nuclide=U238&activity=3uCi&initialage=20y&timespan=22y&actunits=ci"
   - "interspec://decay/chart?nuc=Ba133&act=3uCi&nuc=Cs137&act=2uci&actunits=ci&timespan=20y"
   - "interspec://decay/chain?nuc=232-th"
   
   @param path The path specified in the URI; must be one of three values "chart", "chain", or
          "calc".  So for example, in the URI "interspec://decay/chain?nuc...", the path is "chain".
          This value determines which tab of the tool will be shown.
   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://decay/chain?nuclide=U238&...", then this string would be "nuclide=U238&...".
          This string is is in standard URL format of "key1=value1&key2=value2&..." with the
          ordering only mattering if there is more than one nuclide specified.  Possible key values
          are listed below.  Capitalization is not important.
   
   Possible url key-value values:
   - "time", "timespan": The timespan fo display on the chart or calculator.  The value may be
     specified using standard time units, or may be specified in half-lives (of the first nuclide).
     Examples include: "timespan=20y", "timespan=5hl", "timespan=200seconds", "timespan=3half-lives"
   - "actunits": Must have a value of "becquerel", "bq", "curie", or "ci". This is the style of
     units that will be used to display things to the user.
   - "iso", "nuc", "isotope", "nuclide": these are all synonyms; the value must list an isotopes.
      Examples include "nuc=Am241", "iso=Co60m", "iso=Co-60", "iso=60Co"
   - "initialage", "age": The initial age of the nuclide, at the time of starting the age.
   - "act", "activity": this is the activity of the nuclide.  If more than one nuclide is specified,
     then this is the activity for the immediately preceding nuclide (e.g., you specify nuclide,
     then its activity). The value provided must specify magnitude and units of activity.
     Examples include "act=5uci", "activity=10kbq"
   
   TODO: implement display options like Log-Y scale, Y-axis type (Activity, Gamma Rate, etc), and grid-lines.
   TODO: implement showing QR code for the tools current state, and put in the DecayWindow footer. 
   */
  void handleAppUrl( std::string path, std::string query_str );
  
public:
  friend class PeakCsvResource;
  friend class DateLengthCalculator;

  InterSpec *m_viewer;
  
  /** \TODO: make this an actual widget and instead of storing in m_nuclides, and just store only in m_nuclidesAddedDiv */
  struct Nuclide
  {
    int a;
    int z;
    int iso;
    double age;
    double activity;
    bool useCurrie;
    Wt::WContainerWidget *display;
    Wt::WText *txt;
    
    void updateTxt();
  };//struct Nuclide

  std::vector<Nuclide>         m_nuclides;
  Wt::WLineEdit               *m_displayTimeLength;
  Wt::WComboBox               *m_displayActivityUnitsCombo;
  Wt::WLabel                  *m_displayActivityUnitsLabel;
  Wt::WCheckBox               *m_logYScale;
  Wt::WCheckBox               *m_photopeakLogYScale;
  Wt::WCheckBox               *m_showGridLines;

  Wt::WComboBox               *m_yAxisType;
  
  //things for the adding of activity
  Wt::WContainerWidget        *m_parentNuclidesDiv;
  Wt::WContainerWidget        *m_addParentNuclideDiv;
  Wt::WContainerWidget        *m_nuclidesAddedDiv;
  Wt::WPushButton             *m_createNewNuclideButton;

  Wt::WPushButton             *m_clearNuclidesButton;

  AuxWindow                   *m_nuclideSelectDialog;
  DecaySelectNuclide          *m_nuclideSelect;

  //Objects to actually display the chart of activities
  Wt::WTabWidget              *m_chartTabWidget;
  DecayActivityChart          *m_decayChart;
  DecayActivityModel          *m_decayModel;
  
  DecayChainChart           *m_decayChainChart;
    
  AuxWindow                   *m_moreInfoDialog;

  Wt::WContainerWidget        *m_decayLegend;

  DateLengthCalculator *m_calc;
  
  CsvDownloadGui *m_csvDownloadDialog; //null if not displayed
  
  //Below are variable updated by refreshDecayDisplay() that reflect the current
  // stuff being showed by m_decayChart.  They are only used in
  // updateMouseOver(...) and can probably just be gotten rid of
  double           m_currentTimeUnits;
  double           m_currentTimeRange;
  
  int m_currentNumXPoints;
  
  SandiaDecay::NuclideMixture            *m_currentMixture;


  //Functions for dealing with adding new nuclides
  void addSelectedNuclide();
  
  void addTheNuclide( const NuclideSelectedInfo &nuc );
  
  void checkTimeRangeValid();
  
  
  void updateYScale();
  void refreshDecayDisplay( const bool update_calc );
  void addDecaySeries();
  void userSetShowSeries( int series, bool show );
  void updateYAxisRange();
  
  void createCsvDownloadGui();
  void deleteCsvDownloadGui();
  void deleteCsvDownloadGuiTriggerUpdate();
    
  void updateMouseOver( const Wt::WMouseEvent &event );
  void removeNuclide( Wt::WContainerWidget *frame );

  //nuclideSelectDialogDone(): makes sure none of the source Nuclides have the
  //  "EditingNuclide" style class associated with them.
  void nuclideSelectDialogDone();
  
  //sourceNuclideDoubleClicked(...): we are taking a double click by the user
  //  on a source Nuclide to mean they want to edit it.  How we keep track of
  //  this internally is by adding a "EditingNuclide" to the WContainerWidget
  //  rendering the Nuclide.  When the Nuclide div is double clicked, we will
  //  popup the nuclide select window, with the values set to those of the
  //  Nuclide.  Then when the nuclide select widget indicates to add a Nuclide
  //  by calling DecayActivityDiv::addNuclide(...) (through signal/slot), we will
  //  see if any of the source Nuclide divs have the "EditingNuclide" style
  //  class, and if so, place the new Nuclide in its place.
  //  This is a bit hacky, or not as robust as I would like, but its the easiest
  //  relatively okay way to add in this editing capability.
  void sourceNuclideDoubleClicked( Wt::WContainerWidget *w );

  double timeToDisplayTill();

  void displayMoreInfoPopup( const double time );
  void decayChartClicked( const Wt::WMouseEvent& event );
  
  Wt::WContainerWidget *nuclideInformation(
                                    const SandiaDecay::Nuclide *nuclide ) const;
  Wt::WContainerWidget *isotopesSummary( const double time ) const;


  void setTimeLimitToDisplay();
  
  void updateInitialMixture() const;
  
  static double findTimeForActivityFrac( const SandiaDecay::NuclideMixture *mixture,
                                         const double fracT0ActivityWanted,
                                         const double searchStartTime = -1.0,
                                         const double searchEndTime = -1.0 );
  void deleteMoreInfoDialog();
  void manageActiveDecayChainNucStyling();

  void init();
  void initCharts();
  void setGridLineStatus();
  
  void globalKeyPressed( const Wt::WKeyEvent &e );
  
  Wt::WContainerWidget *initDisplayOptionWidgets();

  
};//DecayActivityDiv



//DecayActivityDiv_h
#endif
