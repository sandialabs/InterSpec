#ifndef IsotopeSearchByEnergy_h
#define IsotopeSearchByEnergy_h
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

#include <vector>
#include <memory>

#include <boost/any.hpp>

#include <Wt/WString>
#include <Wt/WContainerWidget>

#if ( USE_SPECTRUM_CHART_D3 )
class D3SpectrumDisplayDiv;
#else
class SpectrumDisplayDiv;
#endif


namespace Wt
{
  class WText;
  class WCheckBox;
  class WLineEdit;
  class WTreeView;
  class WPushButton;
  class WDoubleSpinBox;
}//namespace Wt

class InterSpec;
class IsotopeSearchByEnergyModel;


/** A class that displays a user input search energy, and window, as well as
  notifies IsotopeSearchByEnergy of any changes.
 
  \sa IsotopeSearchByEnergyModel
*/
class IsotopeSearchByEnergy : public Wt::WContainerWidget
{
protected:
  /** Class that represents an energy (and its +-range) that a user wants
      to search on.  The IsotopeSearchByEnergy widget will always have at least
      one of these, but the user may add more, to require search results to
      match additional energies.
   */
  class SearchEnergy : public Wt::WContainerWidget
  {
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
    
    Wt::WContainerWidget *m_removeIcn;
    Wt::WContainerWidget *m_addAnotherIcn;
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
