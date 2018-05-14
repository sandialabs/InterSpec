#ifndef DecaySelectNuclideDiv_h
#define DecaySelectNuclideDiv_h
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

#include <map>
#include <string>
#include <vector>

#include <Wt/WSignal>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

#include "InterSpec/AuxWindow.h"

namespace Wt
{
  class WText;
  class WLineEdit;
  class WComboBox;
  class WPushButton;
  class WSelectionBox;
  class WSuggestionPopup;
}//namespace Wt


namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
  class SandiaDecayDataBase;
}//namespace SandiaDecay

struct NuclideSelectedInfo
{
  int z, a, metasable;
  double activity, initialAge;  //in units of PhysicalUnits
  bool useCurrie;  //used purely for display purposes
};//struct NuclideSelectedInfo

//IsotopeNameFilterModel

class SimpleIsotopeNameFilterModel;

class SimpleIsotopeNameFilterModel : public  Wt::WAbstractItemModel
{
  /*
    SimpleIsotopeNameFilterModel is intended to filter user input to guess
    which isotope they want to select.
    */
public:
  SimpleIsotopeNameFilterModel( DecaySelectNuclide *parent );
  virtual ~SimpleIsotopeNameFilterModel();

  virtual Wt::WModelIndex index( int row, int column,
                                 const Wt::WModelIndex &parent = Wt::WModelIndex() ) const;

  virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const;

  virtual int rowCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;
  virtual int columnCount( const Wt::WModelIndex & parent = Wt::WModelIndex() ) const;

  virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const;

  //determineAndRemoveIsoLevel(...) looks for patterns such as 'Co60m', 'Co60meta',
  //  'Co 60 m', etc. to determine if the user is inputting a metts stable state.
  //  The function returns the iso level (right no just 0, 1, or 2), and
  //  removes the portion of the text indicating the meta level, from the input.
  //--taken from IsotopeNameFilterModel 20121014
  static int determineAndRemoveIsoLevel( std::string &label );

  //All returned alpha strings will be lowercase
  //--taken from IsotopeNameFilterModel 20121014
  static void getAlphaAndNumericSubStrs( std::string label,
                                         std::vector<std::string> &alphastrs,
                                         std::vector<std::string> &numericstrs );

  void filter( const Wt::WString &text );
//  Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex & index ) const;


  static void replacerJs( std::string &js );
  static void nuclideNameMatcherJs( std::string &js );

protected:
  DecaySelectNuclide *m_parent;
  const double m_minHalfLife; //seconds
  Wt::WString m_filter;
  std::vector< const SandiaDecay::Element * > m_candidatesElements;
  std::vector< const SandiaDecay::Nuclide * > m_candidatesNuclides;  //what about isomeric state
};//class SimpleIsotopeNameFilterModel



class DecaySelectNuclide : public Wt::WContainerWidget
{
protected:
  const SandiaDecay::SandiaDecayDataBase *m_nuclideDB;
  Wt::WContainerWidget         *m_footer;
  AuxWindow                    *m_auxWindow;
  Wt::WSelectionBox            *m_elementSelection;
  Wt::WSelectionBox            *m_massSelection;
  Wt::WPushButton              *m_acceptButton;
  Wt::WLineEdit                *m_nuclideActivityEdit;
  Wt::WLineEdit                *m_nuclideAgeEdit;
  Wt::WText                    *m_selectedIsotopeHalfLife;
  Wt::WLineEdit                *m_isotopeSearch;
  Wt::WSuggestionPopup         *m_isotopeSuggestions;
  SimpleIsotopeNameFilterModel *m_isoSearchFilterModel;

  Wt::Signal<void> m_doneSignal;
  Wt::Signal< NuclideSelectedInfo > m_selectedSignal;

  void init();
  void initActivityAgeSelects();
  void makeMassList();
  void currentlySelectedIsotope( int &a, int &z, int &meta );
  void updateSelectedHalfLife();
  void updateNuclideSuggestBox();
  void emitDone();
  void emitAccepted();
  void enableAcceptButton();
  int getZ( const std::string &symbol ) const;

public:
  DecaySelectNuclide( const SandiaDecay::SandiaDecayDataBase *decayDatabase,
                    Wt::WContainerWidget *parent = NULL , AuxWindow* window = NULL);
  virtual ~DecaySelectNuclide();

  Wt::Signal<void> &done();
//  Wt::Signal<int,int,int,double,std::string,double> &selected();
  Wt::Signal<NuclideSelectedInfo> &selected();

  void setNuclideSearchToFocus();
  
  void setCurrentInfo( int a, int z, int iso, double age, double activity, bool useCurrie );
  
  //The "Add" button can have either the text "Add" or "Accept" depending
  //  if the intention is to add, or to edit a displayed nuclide.  The following
  //  are convience functions to do this
  void setAddButtonToAdd();
  void setAddButtonToAccept();
  
  friend class SimpleIsotopeNameFilterModel;
};//class DecaySelectNuclide


#endif //DecaySelectNuclideDiv_h
