#ifndef ReferencePhotopeakDisplay_h
#define ReferencePhotopeakDisplay_h
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

#include <deque>
#include <string>
#include <vector>

#include <Wt/WColor>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/ReactionGamma.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/ReferenceLineInfo.h"

class MaterialDB;
class SpectrumChart;
class ShieldingSelect;
class FeatureMarkerWidget;
class MoreNuclideInfoWindow;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml

namespace Wt
{
  class WText;
  class WCheckBox;
  class WLineEdit;
  class WTreeView;
  class WGridLayout;
  class WPushButton;
  class WDoubleSpinBox;
  class WSuggestionPopup;
}//namespace Wt

class InterSpec;
class ColorSelect;
class DetectorDisplay;
class D3SpectrumDisplayDiv;
class IsotopeNameFilterModel;


class DecayParticleModel : public  Wt::WAbstractItemModel
{
  //Model to display
  // TODO: unify info held by DecayParticleModel::RowData and ReferenceLineInfo, so they can share the same struct, and we can just set the data once

public:
  enum Column
  { kEnergy, kBranchingRatio, kResponsibleNuc,
    kDecayMode, kParticleType, kNumColumn
  };//enum Column

  struct RowData
  {
    static const int XRayDecayMode;
    static const int ReactionToGammaMode;
    static const int NormGammaDecayMode;
    static const int CascadeSumMode;

    float energy, branchRatio;
    int decayMode;
    SandiaDecay::ProductType particle;
    const SandiaDecay::Nuclide *responsibleNuc;
  };//struct RowData

  static bool less_than( const RowData &lhs, const RowData &rhs,
                         const Column r, const Wt::SortOrder order );

  DecayParticleModel( Wt::WObject *parent = 0 );
  virtual ~DecayParticleModel();

  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index,
                                            int role = Wt::DisplayRole ) const;
  boost::any headerData( int section,
                         Wt::Orientation orientation = Wt::Horizontal,
                         int role = Wt::DisplayRole ) const;
  virtual int columnCount( const Wt::WModelIndex &parent
                                                    = Wt::WModelIndex() ) const;
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index) const;
  virtual Wt::WModelIndex index( int row, int column,
                      const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );
  void clear();
  virtual void setRowData( const std::vector<RowData> &newData );

  const std::vector<RowData> &rowData() const;
  
protected:
  Column m_sortColumn;
  Wt::SortOrder m_sortOrder;
  std::vector<RowData> m_data;
};//class DecayParticleModel(...)



class ReferencePhotopeakDisplay : public Wt::WContainerWidget
{
public:
  ReferencePhotopeakDisplay(
                        D3SpectrumDisplayDiv *chart,
                        MaterialDB *materialDB,
                        Wt::WSuggestionPopup *materialSuggest,
                        InterSpec *specViewer,
                        Wt::WContainerWidget *parent = 0 );
  virtual ~ReferencePhotopeakDisplay();

  void setFocusToIsotopeEdit();

  //clearAllLines(): removes lines from both this widget (table and entry form),
  //  as well as on the chart.
  void clearAllLines();

  //fitPeaks(): attempts to fit and add the peaks fir the currently showing
  //  nuclide
  void fitPeaks();
  
  //setIsotope(...): sets the nuclde to specified nuc and age; if age is less
  //  than zero, then current age will not be modified
  void setIsotope( const SandiaDecay::Nuclide *nuc, double age = -1.0 );
  
  //setElement(...): sets to display xrays for the desired element
  void setElement( const SandiaDecay::Element *element );
    
  //setReaction(...): sets to display gammas for the desired reaction
  void setReaction( const ReactionGamma::Reaction *rctn );
  
  //setNuclideAndAge(...): sets the widgets nuclide text to passed in value.
  //  If length of age is non-zero, then will set the age string, otherwise
  //  will defer to auto age setting.
  //  'name' is assumed to be a valid source name, and is not checked.
  void setNuclideAndAge( const std::string &name, const std::string &age );
  
  //setShieldingMaterialAndThickness(...): sets the shielding material as well
  //  as thickness.  The name must be either in the material DB, or a chemical
  //  formula (which if so, will be added to the list of user suggestions for
  //  the current session).  The thickness must be accepted by
  //  PhysicalUnits::stringToDistance(...).  Throws exception if either field
  //  is invalid.
  void setShieldingMaterialAndThickness( const std::string &name,
                                         const std::string &thickness );
  
  /** Returns the current ShieldingSelect.  This is intended as an accessor, not for setting the
   shielding.
   */
  const ShieldingSelect *shieldingSelect();
  
  //currentlyShowingNuclide(): as its name implies; returns NULL if not showing
  //  a nuclide.
  const ReferenceLineInfo &currentlyShowingNuclide() const;

  //persistedNuclides(): nuclides the user as asked to be persisted, and are
  //  showing.  Does not include currentlyShowingNuclide(), however the user
  //  may have previously persisted the same nuclide as
  //  currentlyShowingNuclide() so a nuclide could be in both places.
  const std::vector<ReferenceLineInfo> &persistedNuclides() const;
  
  //Returns all nuclides currently showing
  std::vector<ReferenceLineInfo> showingNuclides() const;

  /** Simple accessor for the detectorDisplay widget; used for CSV download of displayed info.  */
  const DetectorDisplay *detectorDisplay() const;
  
  /** Accessor to the shielding select; used for CSV download of displayed info. */
  const ShieldingSelect *shieldingSelect() const;
  
  const DecayParticleModel *particleModel() const;
  
  //serialize(...): Serializes the current state of this widget to the XML
  //  document or string.  Currently, all persisted lines will be assumed to
  //  have the same shielding as the currently saved lines.
  void serialize( std::string &xml_data  ) const;
  void serialize( rapidxml::xml_node<char> *parent_node ) const;

  //deSerialize(...): de-serializes the widget from the XML string created by
  //  serialize(...).
  void deSerialize( std::string &xml_data  );
  
  /** returns a JSON array corresponding to currently showing reference lines.
   Will be in the format: "[{},{},...]", where even if there are no showing
  lines you will get back "[]"
   */
  //std::string jsonReferenceLinesArray();
  
  /** Provides a map from the reference line label (e.x., "Ba133") to the
   actual lines JSON (ex. {...}).
   "*/
  std::map<std::string,std::string> jsonReferenceLinesMap();
  
  //persistCurentLines(): persists current lines
  void persistCurentLines();
  
  /** Set #m_lineColors.  This will be used to find the first available color
  that should be used for the current ref lines.  The behaviour is influenced by
   #m_peaksGetAssignedRefLineColor, with #m_specificSourcelineColors and
   #m_previouslyPickedSourceColors being consulted to determine color to use
   for ref lines.
   */
  void setColors( const std::vector<Wt::WColor> &referenceLineColor );
  
  /** Sets the ref line colors to use for specific sources.  Will be overridden
   by #m_previouslyPickedSourceColors, if it has an entry for the specific
   sources.
   */
  void setColorsForSpecificSources( const std::map<std::string,Wt::WColor> &referenceLineColorForSources );
  
  /** Setting this determines if we will always use the first non-currently
   displayed reference line color in #m_lineColors (set by #setColors), or if
   we will look at the colors of all of the current peaks, as well as currently
   displayed ref lines to search through #m_lineColors for the first non-used
   color.
   */
  void setPeaksGetAssignedRefLineColor( const bool theydo );
  
  /** Sets the nuclide ID from an external service (a seperate exectuable or web-service); e.g.
   if user has it setup to automatically get nuclide ID on file load.
   
   @param algo_name The name to display as the header under "Suggestions" for these results.
          Should be short (less than ~12 characters); if empty will use "External RID".
   @param iso_descrips The pairs of nuclide names, and thier descriptions.  For each pair,
          if the name is a valid nuclide (or rather reference line that can be displayed), then
          the first element should be nuclide name, and second element its description
          (e.g., "Industrial", "SNM", etc).  If not a nuclide, then the second element should
          be empty, which will make it so user cant click on this result to show ref. lines.
          (I know, not a great system, but it will probably need re-vamped anyway once
          we get a little more use-cases and experience).
   
   \sa RemoteRid::startAutomatedOnLoadAnalysis
   */
  void setExternalRidResults( const std::string &algo_name,
                             const std::vector<std::pair<std::string,std::string>> &iso_descrips );
  
  /** Signal emmitted whenever the user selects a new nuclide to be shown. */
  Wt::Signal<> &displayingNuclide();
  
  /** Signal emmitted whenever user clicks "Clear All" button, or goes from a
   valid nuclide, to no cuclide (and there are none persisted).
   */
  Wt::Signal<> &nuclidesCleared();
  
  /** Return the material database this widget uses */
  const MaterialDB *materialDB() const;
  
  /** Returns a shared pointer to `m_undo_redo_sentry` - as long as there are any shared ptrs to
   the sentry alive, an undo/redo step wont be inserted.
   */
  std::shared_ptr<void> getDisableUndoRedoSentry();
  
  /** Will return the showing FeatureMarkerWidget, or null if not currently showing. */
  FeatureMarkerWidget *featureMarkerTool();
  
  /** Returns valid FeatureMarkerWidget ptr. Creates a FeatureMarkerWidget, if not currently showing. */
  FeatureMarkerWidget *showFeatureMarkerTool();
  
  /** Removes FeatureMarkerWidget, if currently showing. */
  void removeFeatureMarkerTool();
  
  /** Called when user checks or unchecks to show feature markers. */
  void featureMarkerCbToggled();
  
  /** Blinks the feature marker widget a bit. */
  void emphasizeFeatureMarker();
  
  void programmaticallyCloseMoreInfoWindow();
  
  MoreNuclideInfoWindow *moreInfoWindow();
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  /** Currently just hides or shows the table with gamma info. */
  void setNarrowPhoneLayout( const bool narrow );
#endif
  
  /** Pass in a source string (U235, Fe(n,n), Pb, etc), and get back the color that has been previously used
   by the user.
   Checks `currentlyUsedPeakColors()`, then  `m_previouslyPickedSourceColors`, then
   `m_specificSourcelineColors`.  If no history of it can be found, returns default/empty color.
   
   Note that `src` passed in must exactly match nuclide/reaction/element.
   
   This function is intended to be called by other widgets.
   */
  Wt::WColor suggestColorForSource( const std::string &src ) const;
  
  /** Updates `m_previouslyPickedSourceColors` for the specified source.
   
   If source string is empty, no action is taken.
   If color is default, then the source in the cache is removed.
   */
  void updateColorCacheForSource( const std::string &source, const Wt::WColor &color );
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  
  void addUndoRedoPoint( const ReferenceLineInfo &starting_showing,
                        const std::vector<ReferenceLineInfo> &starting_persisted,
                        const std::deque<RefLineInput> &starting_prev_nucs,
                        const bool starting_user_color,
                        const RefLineInput &current_user_input );
  
  void updateDisplayChange();
  void updateDisplayFromInput( RefLineInput user_input );
 
  RefLineInput userInput() const;

  static std::vector<DecayParticleModel::RowData> createTableRows( const ReferenceLineInfo &refLine );
  
  Wt::WColor colorForNewSource( const std::string &src );

  void handleIsotopeChange( const bool useCurrentAge );
  
  void handleDrfChange( std::shared_ptr<DetectorPeakResponse> det );

  //refreshLinesDisplayedToGui(): makes setting and re-sends to client the lines
  //  that should be displayed, based on m_currentlyShowingNuclide and
  //  m_persisted objects
  void refreshLinesDisplayedToGui();
  
  void userColorSelectCallback( const Wt::WColor &color );
  
  /** Returns a map from source name (ex "U235", "Be(n,n')", etc) to colors used
   for peaks with that as a source, as well as for persisted reference lines.
   */
  std::map<std::string,std::vector<Wt::WColor>> currentlyUsedPeakColors() const;
  
  void toggleShowOptions();

  void updateOtherNucsDisplay();
  void updateAssociatedNuclides();
  void showMoreInfoWindow();
  void handleMoreInfoWindowClose( MoreNuclideInfoWindow *window );

  /** Function that will be called whenever any displayed spectrum
  gets changed (different file, or sample numbers).

  If foreground spectrum _file_ changes, then the RIID analysis results
  of the siggested nuclides may need updating.  However, for simplicity, 
   we'll update suggested nuclides whenever the foreground gets updated
  */
  void handleSpectrumChange(SpecUtils::SpectrumType type);

  D3SpectrumDisplayDiv *m_chart;

  InterSpec *m_spectrumViewer;
  
  /** If we are setting the GUI state using a RefLineInput object, or serializing from XML,
   we dont want to react to signals from shielding widget, or DRF, or whatever that they
   have changed, otherwise we can get into an infinite recursive state.
   */
  bool m_currently_updating;
  
  /** if `m_undo_redo_sentry.lock()` yields a valid pointer, than an undo/redo step wont be inserted.
   \sa getDisableUndoRedoSentry();
   */
  std::weak_ptr<void> m_undo_redo_sentry;
  
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WSuggestionPopup *m_nuclideSuggest;

  Wt::WLineEdit *m_ageEdit;
  Wt::WDoubleSpinBox *m_lowerBrCuttoff;

  Wt::WCheckBox *m_promptLinesOnly;

  Wt::WText *m_halflife;
  Wt::WPushButton *m_moreInfoBtn;
    
  Wt::WPushButton *m_persistLines;
  Wt::WPushButton *m_clearLines;
  //Wt::WPushButton *m_fitPeaks;

  Wt::WPushButton *m_options_icon;
  Wt::WContainerWidget *m_options;
  Wt::WContainerWidget *m_optionsContent;

  Wt::WCheckBox *m_showGammas;
  Wt::WCheckBox *m_showXrays;
  Wt::WCheckBox *m_showAlphas;
  Wt::WCheckBox *m_showBetas;
  Wt::WCheckBox *m_showCascadeSums;
  Wt::WText* m_cascadeWarn;
  Wt::WCheckBox *m_showEscapes;
  Wt::WCheckBox *m_showRiidNucs;
  Wt::WCheckBox *m_showPrevNucs;
  Wt::WCheckBox *m_showAssocNucs;
  Wt::WCheckBox *m_showFeatureMarkers;

  Wt::WContainerWidget *m_otherNucsColumn;
  Wt::WContainerWidget *m_otherNucs;

  const size_t m_max_prev_nucs = 8; //arbitrary
  
  std::deque<RefLineInput> m_prevNucs;
  
  /** Name of "external" RID algorithm used, as set by #setExternalRidResults. */
  std::string m_external_algo_name;
  /** "External" RID results, as set by #setExternalRidResults. */
  std::vector<std::pair<std::string,std::string>> m_external_ids;
  
  Wt::WContainerWidget *m_featureMarkerColumn;
  
  DetectorDisplay *m_detectorDisplay;
  MaterialDB *m_materialDB;                 //not owned by this object
  Wt::WSuggestionPopup *m_materialSuggest;  //not owned by this object
  ShieldingSelect *m_shieldingSelect;

  Wt::WTreeView *m_particleView;
  DecayParticleModel *m_particleModel;

  std::vector<ReferenceLineInfo> m_persisted;
  ReferenceLineInfo m_currentlyShowingNuclide;
  
  ColorSelect *m_colorSelect;
  Wt::WInteractWidget *m_csvDownload;
  
  /** Indicates if the user has picked a color for the current source.
   */
  bool m_userHasPickedColor;
  
  /** Only effects how colors are chosen for reference line, does not influence
      coloring of peaks.
   */
  bool m_peaksGetAssignedRefLineColor;
  
  /** A cache of the most recent color the user has chosen for a particular
      source.  This is consulted before #m_specificSourcelineColors.
   */
  std::map<std::string,Wt::WColor> m_previouslyPickedSourceColors;
  
  /** List of colors that will be used for reference lines.  */
  std::vector<Wt::WColor> m_lineColors;
  
  /** Colors of referfence lines to use for specific sources.
   Will be consulted after #m_previouslyPickedSourceColors.
   */
  std::map<std::string,Wt::WColor> m_specificSourcelineColors;
  
  Wt::Signal<> m_displayingNuclide;
  
  Wt::Signal<> m_nuclidesCleared;
  
  MoreNuclideInfoWindow *m_nucInfoWindow;
  
  FeatureMarkerWidget *m_featureMarkers;
  
  static const int sm_xmlSerializationVersion;
};//class ReferencePhotopeakDisplay


#endif  //ifndef ReferencePhotopeakDisplay_h
