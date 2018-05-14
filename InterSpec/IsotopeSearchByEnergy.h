#ifndef IsotopeSearchByEnergy_h
#define IsotopeSearchByEnergy_h
/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <vector>
#include <memory>

#include <boost/any.hpp>

#include <Wt/WString>
#include <Wt/WModelIndex>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>
#include <Wt/WAbstractItemDelegate>

#include "InterSpec/AuxWindow.h"
#include "sandia_decay/SandiaDecay.h"
#include "InterSpec/ReactionGamma.h"

#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif

namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
  struct Transition;
  struct RadParticle;
}//namespace SandiaDecay

namespace Wt
{
  class WCheckBox;
  class WLineEdit;
  class WPushButton;
  class WDoubleSpinBox;
  class WTreeView;
}//namespace Wt

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


class InterSpec;

class CustomAbstractItemDelegate : public Wt::WAbstractItemDelegate
{
  
public:
  CustomAbstractItemDelegate(WObject *parent = 0)
  : WAbstractItemDelegate(parent)
  { }
protected:
  Wt::WWidget* update	(	Wt::WWidget * 	widget,
                                               const Wt::WModelIndex & 	index,
                       Wt::WFlags< Wt::ViewItemRenderFlag > 	flags
                                               );
};

class IsotopeSearchByEnergyModel : public Wt::WAbstractItemModel
{
public:
  enum Column
  {
    ParentIsotope, Energy, Distance, BranchRatio,
    SpecificIsotope, ParentHalfLife, AssumedAge, NumColumns
  };//enum Column
  
  enum RadSource
  {
    kGamma    = 0x1,
    kXRay     = 0x2,
    kReaction = 0x4
  };//enum RadSource
  
  struct IsotopeMatch
  {
    IsotopeMatch();
    
    double m_distance;    //sum of distance over all energies searched
    
    double m_age;         //age assumed for listing things
    double m_branchRatio; //branching ratio
    
    //Only one of the following will be valid: m_nuclide, m_element, m_reaction
    const SandiaDecay::Nuclide     *m_nuclide;
    const SandiaDecay::Transition  *m_transition;
    const SandiaDecay::RadParticle *m_particle;
    PeakDef::SourceGammaType m_sourceGammaType;
    
    const SandiaDecay::Element *m_element;
    const SandiaDecay::EnergyIntensityPair *m_xray;
    
    const ReactionGamma::Reaction *m_reaction;
    ReactionGamma::EnergyAbundance m_reactionEnergy;
    
    Wt::WString m_displayData[NumColumns];
  };//struct IsotopeMatch
  
public:
  IsotopeSearchByEnergyModel( Wt::WObject *parent = 0 );
  virtual ~IsotopeSearchByEnergyModel();
  
  void clearResults();
  
  
  //SearchWorkingSpace: this is a place to store search results, and what led to
  //  them, before changing the data of the IsotopeSearchByEnergyModel.  This
  //  struct is necasary since the search results are computed in a background
  //  thread, and then the results are posted to the WApplications event loop,
  //  of which the call to (including alll input parameters) must be bound
  //  before posting the search to the background thread (for safety against
  //  calling updateSearchResults() after *this has been delted).
  struct SearchWorkingSpace
  {
    std::vector<double> energies;
    std::vector<double> windows;
    std::vector< std::vector<IsotopeMatch> > matches;
    Column sortColumn;
    Wt::SortOrder sortOrder;
    boost::function< void(void) > searchdoneCallback;
  };//struct SearchWorkingSpace

  void updateSearchResults( std::shared_ptr<SearchWorkingSpace> workingspace );
  static void setSearchEnergies( std::shared_ptr<SearchWorkingSpace> workingspace,
                          const double minbr,
                          const double minHalfLife,
                          Wt::WFlags<RadSource> radiation,
                          const std::string appid,
                          boost::function< void(void) > updatefcn );
  
  
  virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;
  virtual boost::any data( const Wt::WModelIndex &index,
                          int role = Wt::DisplayRole) const;
  
  //nuclide(...): returns the nuclide IF it is a nuclide defined row, otherwise
  //  NULL
  const SandiaDecay::Nuclide *nuclide( const Wt::WModelIndex &index ) const;

  //nuclide(...): returns the element IF it is a xray defined row, otherwise
  //  NULL
  const SandiaDecay::Element *xrayElement( const Wt::WModelIndex &index ) const;
  
  //nuclide(...): returns the reaction IF it is a reaction defined row,
  //  otherwise NULL
  const ReactionGamma::Reaction *reaction( const Wt::WModelIndex &index ) const;
  
  double assumedAge( const Wt::WModelIndex &index ) const;
  
  virtual boost::any headerData( int section,
                                Wt::Orientation orientation = Wt::Horizontal,
                                int role = Wt::DisplayRole ) const;
  
  virtual Wt::WModelIndex index( int row, int column,
                                const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;
  
  virtual void sort( int column, Wt::SortOrder order = Wt::AscendingOrder );

  virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const;

  

  static void sortData( std::vector< std::vector<IsotopeMatch> > &input,
                       const std::vector<double> &energies,
                       int column, Wt::SortOrder order );
  
  typedef std::map<const SandiaDecay::Nuclide *, std::set<double> > NucToEnergiesMap;
  typedef std::vector< std::vector<IsotopeSearchByEnergyModel::IsotopeMatch> > SearchResults;
  
  //nuclidesWithAllEnergies(...): allows x-rays to be considered gamma rays
  static void nuclidesWithAllEnergies( const NucToEnergiesMap &candidates,
                                       const std::vector<double> &energies,
                                       const std::vector<double> &windows,
                                       const double minbr,
                                       SearchResults &answer );
  
  //xraysWithAllEnergies(...): fairly inefficient
  static void xraysWithAllEnergies( const std::vector<double> &energies,
                                    const std::vector<double> &windows,
                                    SearchResults &answer );
  
  //reactionsWithAllEnergies(...): not very well yet.  Does not take into
  //  account the minimum branching ratio desired
  static void reactionsWithAllEnergies( const std::vector<double> &energies,
                                        const std::vector<double> &windows,
                                        SearchResults &answer );
  
  Column sortColumn() const;
  Wt::SortOrder sortOrder() const;
  
protected:
  Column m_sortColumn;
  Wt::SortOrder m_sortOrder;
  std::vector<double> m_windows;
  std::vector<double> m_energies;
  std::vector< std::vector<IsotopeMatch> > m_matches;
};//IsotopeSearchByEnergyModel


class IsotopeSearchByEnergy : public Wt::WContainerWidget
{
protected:
  class SearchEnergy : public Wt::WContainerWidget
  {
    //A class that displays a user input search energy, and window, as well as
    //  notifies IsotopeSearchByEnergy of any changes.
  public:
    SearchEnergy( Wt::WContainerWidget *parent = 0 );
    virtual ~SearchEnergy(){}
  
    double energy() const;
    double window() const;
    void setEnergy( double energy );
    void setWindow( double window );
    
    void enableRemove();
    void disableRemove();
    
    void enableAddAnother();
    void disableAddAnother();
    
    Wt::Signal<> &enter();
    Wt::Signal<> &remove();
    Wt::Signal<> &changed();
    Wt::Signal<> &gotFocus();
    Wt::Signal<> &addAnother();
    
  protected:
    void emitEnter();
    void emitRemove();
    void emitChanged();
    void emitGotFocus();
    void emitAddAnother();
    
    Wt::WText *m_removeIcn;
    Wt::WText *m_addAnotherIcn;
    Wt::WDoubleSpinBox *m_energy;
    Wt::WDoubleSpinBox *m_window;
    Wt::Signal<> m_enter, m_changed, m_remove, m_focus, m_addAnother;
  };//class SearchEnergy
  
  
public:
  IsotopeSearchByEnergy(  InterSpec *viewer,
#if ( USE_SPECTRUM_CHART_D3 )
                          D3SpectrumDisplayDiv *chart,
#else
                          SpectrumDisplayDiv *chart,
#endif
                          Wt::WContainerWidget *parent = 0 );
  
  virtual ~IsotopeSearchByEnergy();
  
  //startSearch():  starts a job (in a thread outside of the event loop) that
  //  updates the IsotopeSearchByEnergyModel results with the displayed
  //  SearchEnergy widgets, min BR, and min HL.
  void startSearch( const bool refreshBr );
  
  //addNewSearchEnergy(): adds a SearchEnergy widget; iniatial energy 0.0, and
  //  returns created widget
  SearchEnergy *addNewSearchEnergy();
  
  //addSearchEnergy(): a convience function for calling addNewSearchEnergy()
  //  inside signals/slots
  void addSearchEnergy();
  
  //setNextSearchEnergy(...): intended to be called when the user clicks on the
  //  chart; set the value of the SearchEnergy widget pointed to by
  //  m_nextSearchEnergy to be the value passed in, then increments
  //  m_nextSearchEnergy
  void setNextSearchEnergy( double energy, double window = -1.0 );
  
  //searchEnergyRecievedFocus(...): called when the user clicks on a search
  //  energy widget, so that if the user then clicks on the spectrum to set
  //  the energy, the SearchEnergy passed into this function will be set.
  void searchEnergyRecievedFocus( SearchEnergy *enrgy );
  
  //removeSearchEnergy(...): delets the passed in SearchEnergy widget, and
  //  updates m_nextSearchEnergy if necessary
  void removeSearchEnergy( SearchEnergy *energy );
  
  //loadSearchEnergiesToClient(): sets the search energies to the client canvas
  //  for dragging (e.g. the spectrum) will show the bands representing the
  //  search range
  void loadSearchEnergiesToClient();
  
  //clearSearchEnergiesOnClient(): clears the search energies displayed on the
  //  spectrum
  void clearSearchEnergiesOnClient();
  
  //minBrOrHlChanged(): intiialtes a job in the thread pool to update the
  //  EnergyToNuclideServer to new min BR and HL
  void minBrOrHlChanged();

  //resultSelectionChanged(): called when user selects a row of the results
  void resultSelectionChanged();
  
  //serialize(): serializes current search energies to xml
  void serialize( std::string &xml_data ) const;
  
  //deSerialize(): de-serializes search eneries from the xml passed in
  //  If renderOnChart is specified, the search energies/ranges will be
  //  loaded to the chart, and search performed.
  void deSerialize( std::string &xml_data, const bool renderOnChart );
  
  //searches(): returns the current search widgets
  std::vector<SearchEnergy *> searches();
  
protected:
  //hideSearchingTxt(): hides the searching text only if
  //  searchNum==m_currentSearch.  This is incase multiple searches are posted
  //  to the thread pool, only the most recent one will hide the "Searching"
  //  text.
  void hideSearchingTxt( const int searchNum );
  
  InterSpec *m_viewer;
#if ( USE_SPECTRUM_CHART_D3 )
  D3SpectrumDisplayDiv *m_chart;
#else
  SpectrumDisplayDiv *m_chart;
#endif
  Wt::WContainerWidget *m_searchEnergies;
  int m_currentSearch;
  Wt::WText *m_searching;
  Wt::WTreeView *m_results;
  Wt::WDoubleSpinBox *m_minBranchRatio;
  Wt::WLineEdit *m_minHalfLife;
  IsotopeSearchByEnergyModel *m_model;
  
  Wt::WCheckBox *m_gammas;
  Wt::WCheckBox *m_xrays;
  Wt::WCheckBox *m_reactions;
  
  size_t m_nextSearchEnergy;
  double m_minBr, m_minHl;
  
  static const int sm_xmlSerializationVersion;
};//class IsotopeSearchByEnergy


#endif  //ifndef IsotopeSearchByEnergy_h
