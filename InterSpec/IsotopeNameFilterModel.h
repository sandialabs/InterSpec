#ifndef IsotopeNameFilterModel_h
#define IsotopeNameFilterModel_h
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

#include <Wt/WSignal>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/ReactionGamma.h"

namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
  class SandiaDecayDataBase;
}//namespace SandiaDecay

namespace Wt
{
  class WLineEdit;
  class WSuggestionPopup;
}


class IsotopeNameFilterModel : public  Wt::WAbstractItemModel
{
  /*
    IsotopeNameFilterModel is intended to filter user input to guess
    which isotope they want to select.  Does not suggest stable isotopes, or 
    isotopes which do not have any transitions.
    */
public:
  IsotopeNameFilterModel( Wt::WObject *parent );
  virtual ~IsotopeNameFilterModel();

  virtual Wt::WModelIndex index( int row, int column,
                                 const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;

  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;

  virtual int rowCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;

  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const;

  void addCustomSuggestPossibility( const std::string &str );
  
  //determineAndRemoveIsoLevel(...) looks for patterns such as 'Co60m', 'Co60meta',
  //  'Co 60 m', etc. to determine if the user is inputting a metts stable state.
  //  The function returns the iso level (right no just 0, 1, or 2), and
  //  removes the portion of the text indicating the meta level, from the input.
  //--taken from IsotopeNameFilterModel 20121014
  static int determineAndRemoveIsoLevel( std::string &label );

  //All returned alpha strings will be lowercase
  //--taken from IsotopeNameFilterModel 20121014, but a slightly different implementation now
  static void getAlphaAndNumericSubStrs( std::string label,
                                         std::vector<std::string> &alphastrs,
                                         std::vector<std::string> &numericstrs );

  //XXX - filter(...) may be fairly inefficient, and should be cleaned up
  virtual void filter( const Wt::WString &text );
//  Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex & index ) const;
  
  //originalText: what user has typed in
  //
  static std::vector<const ReactionGamma::Reaction *> suggestReactions(
                                            const Wt::WString &originalText,
            const std::vector<const SandiaDecay::Nuclide *> &suggestion_nucs,
            const std::set<const SandiaDecay::Element *> &candidate_elements );
  
  
  static void suggestNuclides( const std::vector<std::string> &alphastrs,
                              const std::vector<std::string> &numericstrs,
                              const int metalevel,
                    const std::set<const SandiaDecay::Element *> &candidate_elements,
                    std::vector<const SandiaDecay::Nuclide *> &suggestions,
                std::vector< const SandiaDecay::Element * > &suggest_elements );
  
  static std::set<const SandiaDecay::Element *> possibleElements( const std::vector<std::string> &alphastrs );
  
  static void replacerJs( std::string &js );
  static void nuclideNameMatcherJs( std::string &js );

  //WSuggestionPopup refilter popup can fail when the user types in realy fast,
  //  at least for Wt 3.3.4 (I think this wa fixed in later versions) so here
  //  we will do a hack to catch this.
  static void setQuickTypeFixHackjs( Wt::WSuggestionPopup *popup );
  
  //Fixes multiple issues with WSuggestionPopup:
  //  1. Race condition: when user types quickly and hits Enter, the suggestion popup
  //     may not have filtered yet, causing wrong suggestion to be selected. This
  //     intercepts Enter key and checks if typed text exactly matches any suggestion,
  //     using that match instead of the highlighted one.
  //  2. Arrow key navigation: ensures arrow key selection works correctly when hitting Enter.
  //  3. Mouse click handling: prevents premature change events when clicking suggestions,
  //     ensuring only the final suggestion value triggers a change event.
  //  4. Popup hiding: automatically hides the suggestion popup when values are entered.
  //  5. Special character input: allows special characters like "(" to be entered without
  //     being blocked by the suggestion popup.
  static void setEnterKeyMatchFixJs( Wt::WSuggestionPopup *popup, Wt::WLineEdit *edit );
  
  /** By default nuclides are included in suggestions, this function allows
   * disabling/enabling this.
   */
  void excludeNuclides( const bool exclude = true );
  
  /** By default xrays (just the element name) are included in suggestions, this
   * function allows disabling/enabling this. 
   */
  void excludeXrays( const bool exclude = true );
  
  /** By default escapes (S.E., D.E.( are included in suggestions, this function
   * allows disabling/enabling that functionality.
   */
  void excludeEscapes( const bool exclude = true );
  
  /** By default nuclear reactions are included in suggestions, this allows 
   * disabling. 
   */
  void excludeReactions( const bool exclude = true );
  
protected:
  const double m_minHalfLife; //seconds
  
  bool m_includeXray;
  bool m_includeEscape;
  bool m_includeNuclides;
  bool m_includeReactions;
  
  Wt::WString m_filter;
  std::vector<std::string> m_customPotentials;
  
  std::vector< const SandiaDecay::Element * > m_candidatesElements;
  std::vector< const SandiaDecay::Nuclide * > m_candidatesNuclides;  //what about isomeric state
  std::vector< const ReactionGamma::Reaction * > m_candidatesReactions;
  std::vector<Wt::WString> m_customSuggests;
  
  //m_typePrefix: Is empty for most user inputs, but if the user is inputting
  //  single or double escape peaks, then may be "S.E. " or "D.E. ", and is
  //  only returned as part of data(...) for nuclides and reactions.
  std::string m_typePrefix;
};//class IsotopeNameFilterModel


#endif
