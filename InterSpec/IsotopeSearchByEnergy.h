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

#include "InterSpec/ReactionGamma.h" //Because we cant forward declare ReactionGamma::Reaction

class D3SpectrumDisplayDiv;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
  class WSplitButton;
}//namespace Wt

namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
}//namespace SandiaDecay


namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml

class InterSpec;
class RowStretchTreeView;
class NativeFloatSpinBox;
class IsotopeSearchByEnergyModel;


/** A class that displays a user input search energy, and window, as well as
  notifies IsotopeSearchByEnergy of any changes.
 
  \sa IsotopeSearchByEnergyModel
*/
class IsotopeSearchByEnergy : public Wt::WContainerWidget
{
public:
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
    NativeFloatSpinBox *m_energy;
    NativeFloatSpinBox *m_window;
    Wt::Signal<> m_enter, m_changed, m_remove, m_focus, m_addAnother;
  };//class SearchEnergy
  
  
public:
  IsotopeSearchByEnergy(  InterSpec *viewer,
                          D3SpectrumDisplayDiv *chart,
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

  /** `m_search_catagory_select` is changed. */
  void catagoryChanged( const bool update_results );
  
  //resultSelectionChanged(): called when user selects a row of the results
  void resultSelectionChanged();
  
  /** Returns how many search energies are reasonably intended to be on peaks. */
  int numSearchEnergiesOnPeaks();
  
  /** Returns how many peaks plausibly line up with a current reference line energy.
   
   @param require_peaks_with_no_id If true, only return count of peaks with no current nuclide ID
   
   Currently doesnt do anything fancy (like Rel Eff analysis, etc) to determine this, but just niavely
   does stupid matching.
   */
  int numCurrentNuclideLinesOnPeaks( const bool require_peaks_with_no_id );
  
  /** Assigns peaks that are reasonably near search energies, to the selected rows nuclide. */
  void assignSearchedOnPeaksToSelectedNuclide();
  
  /** Assigns peaks near currently selected reference lines, to the current nuclide.
   
   @param require_no_peak_id If true, only consider peaks that dont currently have a nuclide/x-ray/reaction assigned.
   */
  void assignPeaksNearReferenceLinesToSelectedNuclide( const bool require_no_peak_id );
  
  /** Clears current reference lines, and selection in table. */
  void clearSelectionAndRefLines();
  
  /** Shows/hides the clear Ref Lines/selection button. */
  void updateClearSelectionButton();
  
  /** Currently just calls updateClearSelectionButton(); the below describes some commented out cuntionality that can likely be removed.
   
   Makes sure the highlighted row in the results matches the displayed reference lines,
   and also calls #updateClearSelectionButton.  If the selected row didnt previously match
   the ref. lines, then the first matching result will be selected.
   
   This function is primarily to facilitate the undo/redo functionality (not perfectly, but close-enough
   to allow us to rely on the ref line widget to handle changing selection rows).
   */
  void handleRefLinesUpdated();
  
  //serialize(): serializes current search energies to xml
  void serialize( std::string &xml_data ) const;
  
  //deSerialize(): de-serializes search eneries from the xml passed in
  //  If renderOnChart is specified, the search energies/ranges will be
  //  loaded to the chart, and search performed.
  void deSerialize( std::string &xml_data, const bool renderOnChart );
  
  //searches(): returns the current search widgets
  std::vector<SearchEnergy *> searches();
  
  /** Returns a shared pointer to `m_undo_redo_sentry` - as long as there are any shared ptrs to
   the sentry alive, an undo/redo step wont be inserted.
   */
  std::shared_ptr<void> getDisableUndoRedoSentry();
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  void setResultTableColumnWidths( const bool narrow );
  
  /** Currently just repositions the results table, and adjusts its row sizes and visible columns. */
  void setNarrowPhoneLayout( const bool narrow );
#endif
protected:
  /** Populates `m_filters` and  `m_typeFilter`.  */
  void initFilterTypes();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

  
  //hideSearchingTxt(): hides the searching text only if
  //  searchNum==m_currentSearch.  This is incase multiple searches are posted
  //  to the thread pool, only the most recent one will hide the "Searching"
  //  text.
  void hideSearchingTxt( const int searchNum );
  
  InterSpec *m_viewer;
  D3SpectrumDisplayDiv *m_chart;
  Wt::WContainerWidget *m_searchConditionsColumn;
  Wt::WContainerWidget *m_searchEnergies;
  Wt::WContainerWidget *m_assignBtnRow;
  Wt::WPushButton *m_clearRefLines;
  Wt::Signals::connection m_refLineUpdateConnection;
  Wt::Signals::connection m_refLineClearConnection;
  Wt::WSplitButton *m_assignPeakToSelected;
  int m_currentSearch;
  Wt::WText *m_searching;
  RowStretchTreeView *m_results;
  
  Wt::WContainerWidget *m_minBranchRatioDiv;
  NativeFloatSpinBox *m_minBranchRatio;
  Wt::WContainerWidget *m_minHalfLiveDiv;
  Wt::WLineEdit *m_minHalfLife;
  IsotopeSearchByEnergyModel *m_model;
  
  struct NucSearchCatagory
  {
    /** Name to display to user in the WComboBox */
    Wt::WString m_name;
    /** A description for the filter type - currently set as tool-top on the WComboBox. */
    Wt::WString m_description;
    
    /** Values of min BR and min HL.  Users can still change these afterwards. */
    double m_minBr, m_minHl;
    
    /** If should search on nuclide gammas and decay x-rays. */
    bool m_nuclides;
    /** If should search on element flourescense x-rays. */
    bool m_fluorescence_xrays;
    /** If should include reactions.  */
    bool m_reactions;
    /** If should search for alpha energy of nuclides. If true, then betas, x-rays, reactions, and nulcides must be false.  */
    bool m_alphas;
    /** If should search for beta end-point energy of nuclides. If true, then alphas, x-rays, reactions, and nulcides must be false. */
    bool m_beta_endpoint;
    
    /** Specfic elements to search for fluorescence x-rays.  If empty, and fluorescence x-rays are allowed, will use all elements. */
    std::vector<const SandiaDecay::Element *> m_specific_elements;
    /** Specific nuclides to search.  If empty and decay x-rays, nuc gammas, alphas, or betas are allowed, will use all nuclides. */
    std::vector<const SandiaDecay::Nuclide *> m_specific_nuclides;
    /** Specific reactions to search.  If empty and reactions are allowed, then all reacitions will be used. */
    std::vector<const ReactionGamma::Reaction *> m_specific_reactions;
    
    /** Sets the information from an XML element.
     
     The passed-in xml node must have name "NucSearchCatagory".
     Throws exception on error.
     
     <NucSearchCatagory>
       <Name>My Fav Nucs</Name>
       <Description>A really good description</Description>
       <MinBr>0.0</MinBr>                    <!-- Defaults to zero if not specified -->
       <MinHl>60 s</MinHl>                   <!-- Defaults to zero if not specified -->
       <AllowNucGammas>1</AllowNucGammas>    <!-- If should search on nuclide gammas and decay x-rays. -->
       <AllowFluorXrays>1</AllowFluorXrays>  <!-- If should search on element flourescense x-rays. -->
       <AllowReactions>1</AllowReactions>    <!-- If should include reactions. -->
       <AllowAlphas>0</AllowAlphas>          <!-- If should search for alpha energy of nuclides. If true, then betas, x-rays, reactions, and nulcides must be false. -->
       <AllowBetas>0</AllowBetas>            <!-- If should search for beta end-point energy of nuclides. If true, then alphas, x-rays, reactions, and nulcides must be false. -->
     
       <Elements>                            <!-- If empty or not specified, will allow all elements. -->
         <Element>U</Element>                 <!-- Use same values you would type into Ref. Photopeak -->
         <Element>Pu</Element>                <!-- Invalid values will be silently discarded -->
       </Elements>
     
       <Nuclides>                            <!-- If empty or not specified, will allow all nuclides. -->
         <Nuclide>Co60</Nuclide>             <!-- Use same values you would type into Ref. Photopeak -->
         <Nuclide>Am241</Nuclide>            <!-- Invalid values will be silently discarded -->
       </Nuclides>
    
       <Reactions>                           <!-- If empty or not specified, will allow all nuclides. -->
         <Reaction>Fe(n,g)</Reaction>        <!-- Use same values you would type into Ref. Photopeak -->
         <Reaction>U(n,n)</Reaction>         <!-- Invalid values will be silently discarded -->
       </Reactions>
     </NucSearchCatagory>
     */
    void deSerialize( const rapidxml::xml_node<char> * const xml_data );
  };//struct NucSearchCatagory
  
  std::vector<NucSearchCatagory> m_search_catagories;
  
  Wt::WComboBox *m_search_catagory_select;
  
  size_t m_nextSearchEnergy;
  double m_minBr, m_minHl;
  
  
  /** A struct that represents the GUI state, using basic types.
   
   We will use this struct as an intermediary to/from XML, as well as for undo/redo state.
   */
  struct WidgetState
  {
    Wt::WString MinBranchRatio;
    Wt::WString MinHalfLife;
    
    size_t NextSearchEnergy;
    size_t CatagoryIndex;
    std::vector<std::pair<double,double>> SearchEnergies;
    
    void serialize( std::string &xml_data ) const;
    
    /** Throws exception on error. */
    void deSerialize( std::string &xml_data, const std::vector<NucSearchCatagory> &catagories );
    
    bool operator==(const WidgetState &rhs) const;
  };//struct WidgetState

  WidgetState guiState() const;
  void setGuiState( const WidgetState &state, const bool renderOnChart );
  
  /** Updates #m_state (based on actual GUI), and adds a undo/redo point, if necassary and appropriate. */
  void addUndoRedoPoint();
  
  
  /** if `m_undo_redo_sentry.lock()` yeilds a valid pointer, than an undo/redo step wont be inserted.
   \sa getDisableUndoRedoSentry();
   */
  std::weak_ptr<void> m_undo_redo_sentry;
  
  WidgetState m_state;
  int m_selected_row;
  
  static const int sm_xmlSerializationVersion;
};//class IsotopeSearchByEnergy


#endif  //ifndef IsotopeSearchByEnergy_h
