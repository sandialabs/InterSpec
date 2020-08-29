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
#include <set>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WSignal>
#include <Wt/WSpinBox>
#include <Wt/WComboBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSelectionBox>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>
#include <Wt/WDoubleValidator>

#include "SpecUtils/StringAlgo.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayActivityDiv.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DecaySelectNuclideDiv.h"
#include "InterSpec/IsotopeNameFilterModel.h"

using namespace Wt;
using namespace std;

// See also: http://www.webtoolkit.eu/wt/blog/2010/03/02/javascript_that_is_c__
#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
  template<class T> struct index_compare_assend
  {
    index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] < arr[b];
    }
    const T arr;
  };//struct index_compare
}



DecaySelectNuclide::DecaySelectNuclide( const bool phone, Wt::WContainerWidget *parent, AuxWindow *auxWindow )
  : WContainerWidget( parent ),
    m_phone( phone ),
    m_footer(auxWindow->footer()),
    m_auxWindow(auxWindow),
    m_elementSelection( NULL ),
    m_massSelection( NULL ),
    m_acceptButton( NULL ),
    m_nuclideActivityEdit( NULL ),
    m_nuclideAgeEdit( NULL ),
    m_selectedIsotopeHalfLife( NULL ),
    m_isotopeSearch( NULL ),
    m_isotopeSuggestions( NULL ),
    m_isoSearchFilterModel( NULL )
{
  addStyleClass( "DecaySelectNuclide" );
  init();
}

DecaySelectNuclide::~DecaySelectNuclide()
{
  if( m_isotopeSuggestions )
    delete m_isotopeSuggestions;
}

//Wt::Signal<int,int,int,double,std::string,double> &DecaySelectNuclide::selected()
Wt::Signal<NuclideSelectedInfo> &DecaySelectNuclide::selected()
{
  return m_selectedSignal;
}

Wt::Signal<void> &DecaySelectNuclide::done()
{
  return m_doneSignal;
}



void DecaySelectNuclide::setNuclideSearchToFocus()
{
  m_isotopeSearch->setFocus();
}


void DecaySelectNuclide::setAddButtonToAdd()
{
  m_acceptButton->setText( "Add" );
  m_acceptButton->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
}

void DecaySelectNuclide::setAddButtonToAccept()
{
  m_acceptButton->setText( "Accept" );
  m_acceptButton->setIcon( "InterSpec_resources/images/accept.png" );
}


void DecaySelectNuclide::setCurrentInfo( int a, int z, int iso,
                                       double age, double activity, bool useCurrie )
{
  using PhysicalUnits::bestActivityUnit;
  using PhysicalUnits::UnitNameValuePair;
  using PhysicalUnits::bestTimeUnitLongHtml;
  using PhysicalUnits::sm_timeUnitNameValues;
  using PhysicalUnits::sm_activityUnitNameValues;
  using PhysicalUnits::sm_timeUnitHtmlNameValues;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  const SandiaDecay::Element *el = db->element( z );
  if( !el )
    return;

  const string agestr = PhysicalUnits::printToBestTimeUnits( age );
  m_nuclideAgeEdit->setText( agestr );

  const string actstr = PhysicalUnits::printToBestActivityUnits( activity, 2, useCurrie );
  m_nuclideActivityEdit->setText( actstr );
  
  
  bool foundElement = false;
  const int nelements = m_elementSelection->count();
  for( int i = 0; i < nelements; ++i )
  {
    const string itemtxt = m_elementSelection->itemText( i ).toUTF8();
    if( itemtxt == el->name )
    {
      foundElement = true;
      m_elementSelection->setCurrentIndex( i );
      break;
    }
  }//for( int i = 0; i < nelements; ++i )
  
  if( !foundElement )  //shouldnt ever happen
  {
    cerr << "DecaySelectNuclide::setCurrentInfo(...): serious logic error, "
            "couldnt find element for z=" << z << endl;
    return;
  }//if( !foundElement )

  makeMassList();
  
  const int nmasses = m_massSelection->count();
  for( int i = 0; i < nmasses; ++i )
  {
    const string masstxt = m_massSelection->itemText(i).toUTF8();
    const string nucsymbol = el->symbol + masstxt;
    const SandiaDecay::Nuclide *nuc = db->nuclide( nucsymbol );
    
    if( !nuc ) //shouldnt ever happen
    {
      cerr << "DecaySelectNuclide::setCurrentInfo(...): couldnt find nuclide for "
           << nucsymbol << endl;
      continue;
    }//if( !nuc )
    
    if( (nuc->massNumber == a) && (nuc->isomerNumber == iso) )
    {
      m_massSelection->setCurrentIndex( i );
      m_isotopeSearch->setText( nuc->symbol );
      break;
    }//if( (nuc->massNumber == a) && (nuc->isomerNumber == iso) )
  }//for( int i = 0; i < nmasses; ++i )
  
  updateSelectedHalfLife();
  enableAcceptButton();
}//void setCurrentInfo(...)


void DecaySelectNuclide::init()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  WContainerWidget::clear();
  
  WLabel *label = 0;
  m_elementSelection         = new WSelectionBox();
  m_massSelection            = new WSelectionBox();
  m_nuclideActivityEdit      = new WLineEdit();//WSpinBox() try WDoubleSpinBox?
  m_nuclideAgeEdit           = new WLineEdit();//WSpinBox(); try WDoubleSpinBox?
  m_selectedIsotopeHalfLife  = new WText( "&lambda;<sub>&frac12;</sub>=",
                                          Wt::XHTMLUnsafeText );
  m_isotopeSearch            = new WLineEdit();
  m_isoSearchFilterModel     = new SimpleIsotopeNameFilterModel( this );
  string matcherJS, replaceJS;
  SimpleIsotopeNameFilterModel::nuclideNameMatcherJs( matcherJS );
  SimpleIsotopeNameFilterModel::replacerJs( replaceJS );
  m_isotopeSuggestions = new WSuggestionPopup( matcherJS, replaceJS );
//  m_isotopeSuggestions = new WSuggestionPopup( m_isotopeSearch );
  wApp->domRoot()->addWidget( m_isotopeSuggestions );
  m_isotopeSuggestions->setMaximumSize( WLength::Auto,
                                        WLength(15, WLength::FontEm) );

  m_isoSearchFilterModel->filter( "" );
  m_isotopeSuggestions->setFilterLength( -1 );
  m_isotopeSuggestions->setModel( m_isoSearchFilterModel );
  m_isotopeSuggestions->filterModel().connect( m_isoSearchFilterModel,
                                      &SimpleIsotopeNameFilterModel::filter );

  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_isotopeSuggestions );
  
  
  m_isotopeSuggestions->forEdit( m_isotopeSearch,
                   WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );
  m_isotopeSearch->enterPressed().connect( this,
                                           &DecaySelectNuclide::emitAccepted );
  
  initActivityAgeSelects();

  m_elementSelection
      ->activated().connect( this, &DecaySelectNuclide::enableAcceptButton );
  m_massSelection
      ->activated().connect( this, &DecaySelectNuclide::enableAcceptButton );
  m_elementSelection
      ->activated().connect( this, &DecaySelectNuclide::updateSelectedHalfLife );
  m_elementSelection
      ->activated().connect( this, &DecaySelectNuclide::makeMassList );
  m_elementSelection
      ->activated().connect( this, &DecaySelectNuclide::updateNuclideSuggestBox );
  m_massSelection
      ->activated().connect( this, &DecaySelectNuclide::updateSelectedHalfLife );
  m_massSelection
      ->activated().connect( this, &DecaySelectNuclide::updateNuclideSuggestBox );
  m_nuclideActivityEdit
      ->changed().connect( this, &DecaySelectNuclide::enableAcceptButton );
  m_nuclideAgeEdit
      ->changed().connect( this, &DecaySelectNuclide::enableAcceptButton );

  
  if( m_phone )
  {
    m_massSelection->setVerticalSize( 10 );
    m_elementSelection->setVerticalSize( 10 );
  }else
  {
    m_massSelection->setVerticalSize( 20 );
    m_elementSelection->setVerticalSize( 20 );
  }//if( !m_phone )
  
  m_massSelection->setSelectionMode( SingleSelection );
  m_elementSelection->setSelectionMode( SingleSelection );
  
  
  
  m_acceptButton = new WPushButton( "Add" ,m_footer);
  m_acceptButton->setFloatSide(Wt::Right);
  m_acceptButton->clicked().connect( this, &DecaySelectNuclide::emitAccepted );
  WPushButton *cancelButton = m_auxWindow->addCloseButtonToFooter("Close");
  cancelButton->clicked().connect( this, &DecaySelectNuclide::emitDone );

  m_elementSelection->addStyleClass( "m_elementSelection" );
  m_acceptButton->addStyleClass("AddIcon");
    

  m_nuclideActivityEdit->addStyleClass( "m_nuclideActivityEdit" );

  m_nuclideAgeEdit->addStyleClass( "m_nuclideAgeEdit" );

  m_selectedIsotopeHalfLife->addStyleClass( "m_selectedIsotopeHalfLife" );

  for( int z = 1; z < 119; ++z )
  {
    const SandiaDecay::Element *el = db->element( z );
    if( !el )
      continue;
    if( db->nuclides( el ).empty() )
      continue;
    m_elementSelection->addItem( el->name );
  }//for each( const NameToZMap::value_type &name_z, m_nameToZMap )

  m_elementSelection->setCurrentIndex( 0 );
  makeMassList();

  WGridLayout *layout = new WGridLayout();
  layout->addWidget( m_elementSelection,        0, 0, 1, 1 );
  layout->addWidget( m_massSelection,           0, 1, 1, 1 );
  
  label = new WLabel( "Nuclide: " );
  layout->addWidget( label, 1, 0, 1, 1 );
  layout->addWidget( m_isotopeSearch, 1, 1, 1, 1 );
  layout->addWidget( m_selectedIsotopeHalfLife, 2, 0, 1, 2, Wt::AlignCenter );
  
  label = new WLabel( "Activity: " );
  layout->addWidget( label,                      3, 0, 1, 1 );
  layout->addWidget( m_nuclideActivityEdit,      3, 1, 1, 1 );
  
  label = new WLabel( "Initial Age: " );
  layout->addWidget( label,                   4, 0, 1, 1 );
  layout->addWidget( m_nuclideAgeEdit,        4, 1, 1, 1 );
  
//  WContainerWidget *buttonDiv = new WContainerWidget();
//  buttonDiv->addWidget( m_acceptButton );
//  buttonDiv->addWidget( cancelButton );
//  layout->addWidget( buttonDiv,                 5, 0, 1, 3, Wt::AlignLeft );
  
  layout->setRowStretch( 0, 1 );
  WContainerWidget::setLayout( layout );

  enableAcceptButton();
}//void initNuclideMenu()



void DecaySelectNuclide::initActivityAgeSelects()
{
  WRegExpValidator *actvalidator = new WRegExpValidator( PhysicalUnits::sm_activityRegex, m_nuclideActivityEdit );
  actvalidator->setFlags(Wt::MatchCaseInsensitive);
  m_nuclideActivityEdit->setValidator( actvalidator );
  m_nuclideActivityEdit->setTextSize( 10 );
  m_nuclideActivityEdit->setText( "1 uCi" );

  m_nuclideAgeEdit->setText( "0.0 us" );
  m_nuclideAgeEdit->setTextSize( 10 );

  WRegExpValidator *agevalidator = new WRegExpValidator( PhysicalUnits::sm_timeDurationHalfLiveOptionalRegex, m_nuclideAgeEdit );
  agevalidator->setFlags(Wt::MatchCaseInsensitive);
  m_nuclideAgeEdit->setValidator(agevalidator);
}//void initActivityAgeSelects()




void DecaySelectNuclide::emitDone()
{
  m_doneSignal.emit();
}


void DecaySelectNuclide::currentlySelectedIsotope( int &a, int &z, int &meta )
{
  try
  {
    const int massIndex    = m_massSelection->currentIndex();
    const int nuclideIndex = m_elementSelection->currentIndex();

    if( (massIndex<0) || (nuclideIndex<0) )
      throw std::runtime_error( "Nuclide or Mass not selected" );

    const WString massWStr = m_massSelection->itemText( massIndex );
    const WString nucWStr  = m_elementSelection->itemText( nuclideIndex );
    const string nucStr    = nucWStr.toUTF8();
    string massStr         = massWStr.toUTF8();

    if( SpecUtils::iends_with( massStr, "m" ) )
       meta = 1;
    else if( SpecUtils::iends_with( massStr, "m2" ) )
       meta = 2;
    else if( SpecUtils::iends_with( massStr, "m3" ) )
      meta = 3;
    else
      meta = 0;

    if( meta==1 )
      massStr = massStr.substr( 0, massStr.length() - 1 );
    else if( meta )
      massStr = massStr.substr( 0, massStr.length() - 2 );

    a = std::stoi( massStr );
    z = getZ( nucStr );
  }catch( exception & /*e*/ )
  {
    //I think all boost/geant4/Wt exceptions inherit from this std::exception
//    cerr << "DecaySelectNuclide::currentlySelectedIsotope(int&, int&): "
//         << "exception caught!" << endl << e.what() << endl;
    meta = 0;
    a = z = -1;
  }//try / catch
}//void currentlySelectedIsotope( int &a, int &z )



void DecaySelectNuclide::emitAccepted()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  try
  {
    int a, z, meta;
    currentlySelectedIsotope( a, z, meta );

    bool validAct  = false, validAge = false;
    double activity = 0.0, age = 0.0;
    const string activityTxt = m_nuclideActivityEdit->text().narrow();
    const string ageText = m_nuclideAgeEdit->text().narrow();
    
    try
    {
      
      activity = PhysicalUnits::stringToActivity( activityTxt );
      validAct = true;
    }catch(...){ }
    
    try
    {
      double hl = 0.0;
      const SandiaDecay::Nuclide *nuc = db->nuclide( z, a, meta );
      if( nuc )
        hl = nuc->halfLife;
      
      age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( ageText, hl );
      validAge = true;
    }catch(...){ }
    
    
    if( (a <= 0) || (z <= 0) || !validAct || !validAge )
    {
      cerr << "DecaySelectNuclide::emitAccepted: shouldnt be here!" << endl;
      return;
    }//if( not correctly selected )

    enableAcceptButton();

    NuclideSelectedInfo selected;
    selected.z = z;
    selected.a = a;
    selected.metasable = meta;
    selected.activity = activity;
    selected.initialAge = age;

    //units activity regex can accept: (bq|becquerel|ci|cu|curie|c)
    const string::size_type unitpos = activityTxt.find_first_of( "CcBb" );
    if( unitpos == string::npos || activityTxt[unitpos]=='C' || activityTxt[unitpos]=='c' )
      selected.useCurrie = true;
    else
      selected.useCurrie = false;
    
    m_selectedSignal.emit( selected );
    m_doneSignal.emit();
  }catch(...)
  {
    cerr << "DecaySelectNuclide::emitAccepted(): exception caught!" << endl;
  }//try / catch
}//void DecaySelectNuclide::emitAccepted()



void DecaySelectNuclide::enableAcceptButton()
{
  try
  {
    const int massIndex    = m_massSelection->currentIndex();
    const int nuclideIndex = m_elementSelection->currentIndex();

    const bool activity  = m_nuclideActivityEdit->validate()==WValidator::Valid;
    const bool age       = m_nuclideAgeEdit->validate()==WValidator::Valid;
    
    if( (massIndex<0) || (nuclideIndex<0) || !activity || !age )
    {
       m_acceptButton->disable();
//      m_acceptButton->setStyleClass( "m_acceptButton Wt-btn with-label Wt-disabled" );
    }else
    {
      m_acceptButton->enable();
//      m_acceptButton->setStyleClass( "m_acceptButton Wt-btn with-label" );
    }
  }catch(...)
  {
    m_acceptButton->disable();
    cerr << "DecaySelectNuclide::enableAcceptButton(): exception caught!" << endl;
  }//try / catch
}//void DecaySelectNuclide::enableAcceptButton()


int DecaySelectNuclide::getZ( const std::string &symbol ) const
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element *el = db->element( symbol );
  
  if( !el )
  {
    stringstream msg;
    msg << "DecaySelectNuclide::getZ() Error: Could not find Z cooresponding to "
        << "symbol '" << symbol << "'";
    throw std::runtime_error( msg.str() );
  }//if( z_name_iter == m_nameToZMap.end() )

  return el->atomicNumber;
}//int getZ( const std::string &symbol )


void DecaySelectNuclide::updateNuclideSuggestBox()
{
  const int massIndex    = m_massSelection->currentIndex();
  const int nuclideIndex = m_elementSelection->currentIndex();

  if( nuclideIndex<0 )
  {
    m_isotopeSearch->setText( "" );
    m_isoSearchFilterModel->filter( "" );
    return;
  }//if( nuclideIndex<0 )

  WString nucWStr = m_elementSelection->itemText( nuclideIndex );
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Element *el = db->element( nucWStr.toUTF8() );

  WString searchstr;
  if( el )
    searchstr = el->symbol;

  if( massIndex >= 0 )
  {
    const WString massWStr = m_massSelection->itemText( massIndex );
    searchstr += massWStr;
  }

  m_isotopeSearch->setText( searchstr );
}//void DecaySelectNuclide::updateNuclideSuggestBox()


void DecaySelectNuclide::updateSelectedHalfLife()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  m_selectedIsotopeHalfLife->setText( "&lambda;<sub>&frac12;</sub>=" );

  int a, z, meta;
  currentlySelectedIsotope( a, z, meta );

  if( (a<=0) || (z<=0) )
    return;

  const SandiaDecay::Nuclide *nuclide = db->nuclide( z, a, meta );

  if( !nuclide )
  {
    cerr << "DecaySelectNuclide::updateSelectedHalfLife() couldnt find z="
         << z << " a=" << a << " iso=" << meta << endl;
    return;
  }//if( !nuclide )

  const string text = "&lambda;<sub>&frac12;</sub>="
                       + PhysicalUnits::printToBestTimeUnits(nuclide->halfLife,2);
  m_selectedIsotopeHalfLife->setText( text );
}//void updateSelectedHalfLife()



void DecaySelectNuclide::makeMassList()
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  const int index = m_elementSelection->currentIndex();

  m_massSelection->clear();

  if( index < 0 )
    return;

  const WString selectedWStr = m_elementSelection->itemText( index );
  const string selectedStr = selectedWStr.toUTF8();
  const SandiaDecay::Element *el = db->element( selectedStr );
  
  if( !el )
  {
    cerr << "DecaySelectNuclide::makeMassList(): invalid element" << endl;
    return;
  }
  
  const vector<const SandiaDecay::Nuclide *> nucs = db->nuclides( el );  
  for( const SandiaDecay::Nuclide *nuc : nucs )
  {
    if( (nuc->halfLife<=0.0) || IsInf(nuc->halfLife) )
      continue;
    
    const int mass = nuc->massNumber;
    const int meta = nuc->isomerNumber;

    stringstream subname;
    subname << mass;

    if( meta )
      subname << "m";
    if( meta > 1 )
      subname << meta;

    m_massSelection->addItem( subname.str() );
  }//for( each atomic mass )

  m_massSelection->setCurrentIndex( -1 );
}//void DecaySelectNuclide::makeMassList()




SimpleIsotopeNameFilterModel::SimpleIsotopeNameFilterModel( DecaySelectNuclide *parent )
  : WAbstractItemModel( parent ),
    m_parent( parent ),
    m_minHalfLife( 0.0 )
{
}

SimpleIsotopeNameFilterModel::~SimpleIsotopeNameFilterModel()
{
}

Wt::WModelIndex SimpleIsotopeNameFilterModel::index( int row, int column,
                                                     const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() || (column!=0) || (row<0) )
    return WModelIndex();

  const int nrows = static_cast<int>( m_candidatesNuclides.size() + m_candidatesElements.size() );

  if( row >= nrows )
    return WModelIndex();

  return createIndex( row, column, (void *)NULL );
}//index(...)


Wt::WModelIndex SimpleIsotopeNameFilterModel::parent( const Wt::WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent(...)


int SimpleIsotopeNameFilterModel::rowCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_candidatesNuclides.size() + m_candidatesElements.size() );
}

int SimpleIsotopeNameFilterModel::columnCount( const Wt::WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return 1;
}

boost::any SimpleIsotopeNameFilterModel::data( const Wt::WModelIndex &index, int role ) const
{
  const int row = index.row();
  const int column = index.column();

  const int nnuc = static_cast<int>( m_candidatesNuclides.size() );
  const int nel = static_cast<int>( m_candidatesElements.size() );
  const int nrows = nnuc + nel;

  if( row<0 || row>=nrows || column!=0 )
    return boost::any();

//  if( role == Wt::ToolTipRole )
//    return WString( "the_tool_tip" );

  if( row < nnuc )
  {
    const SandiaDecay::Nuclide *nuc = m_candidatesNuclides[row];
    //we could customize tool tip or display roles here to give more information...
    return boost::any( WString( nuc->symbol ) );
  }else
  {
    const int elnum = row - nnuc;
    const SandiaDecay::Element *el = m_candidatesElements[elnum];
    //we could customize tool tip or display roles here to give more information...
    return boost::any( WString( el->symbol ) );
  }
}//boost::any SimpleIsotopeNameFilterModel::data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const

int SimpleIsotopeNameFilterModel::determineAndRemoveIsoLevel( std::string &label )
{
  // taken from IsotopeNameFilterModel 20121014
  //determineAndRemoveIsoLevel(...) looks for patterns such as 'Co60m', 'Co60meta',
  //  'Co 60 m', etc. to determine if the user is inputting a metts stable state.
  //  The function returns the iso level (right no just 0, 1, or 2), and
  //  removes the portion of the text indicating the meta level, from the input.

  int user_input_is_meta = 0;
  size_t len = label.size();
  size_t mpos = label.find("meta");

  if( mpos != string::npos )
  {
    if( label.find("meta2") == string::npos )
    {
      label.erase( mpos, 4 );
      user_input_is_meta = 1;
    }else
    {
      label.erase( mpos, 5 );
      user_input_is_meta = 2;
    }
  }else if( label.find('m') != string::npos )
  {
    mpos = 0;
    do
    {
      if( len==0 )
        break;

      mpos = label.find('m', mpos);

      if( mpos == string::npos )
        break;

      if( mpos == (len-1) )
      {
        size_t prev = ((mpos > 0) ? mpos-1 : mpos);
        while( prev > 0 && !isalpha(label[prev]) && !isdigit(label[prev]) )
          --prev;

        if( prev > 0 )
        {
          user_input_is_meta = (isdigit( label[prev] ) ? 1 : 0);
          if( user_input_is_meta )
          {
            label.erase( mpos, 1 );
//            len = label.length();
            break;
          }//if( user_input_is_meta )
        }//if( prev > 0 )
      }else if( mpos > 0 )
      {
        if( (!isdigit(label[mpos+1]) && !isalpha( label[mpos+1]))
            || (label[mpos+1]=='2' && ( (mpos+1)==(len-1) || (!isdigit(label[mpos+2]) && !isalpha( label[mpos+2])) ) ) )
        {
          const bool is2m = (label[mpos+1]=='2'
                             && ( (mpos+1)==(len-1)
                                  || (!isdigit(label[mpos+2]) && !isalpha( label[mpos+2])) ));
          size_t prev = mpos - 1;
          while( prev > 0 && !isalpha(label[prev]) && !isdigit(label[prev]) )
            --prev;
          if( prev > 0 )
          {
            user_input_is_meta = (isdigit( label[prev] ) ? 1 : 0);
            if( user_input_is_meta )
            {
              if( !is2m )
                label.erase( mpos, 1 );
              else
              {
                user_input_is_meta = 2;
                label.erase( mpos, 2 );
              }
//              len = label.length();
              break;
            }//if( user_input_is_meta )
          }//if( prev > 0 )
        }//if( !isdigit(label[mpos+1]) && !isalpha( label[mpos+1]) )
      }//if( mpos == (len-1) ) / else

      mpos += 1;
    }while( mpos != string::npos && mpos < len );
  }//if( has "meta" ) / else has 'm'

  return user_input_is_meta;
}//int determineAndRemoveIsoLevel( std::string &label )


void SimpleIsotopeNameFilterModel::getAlphaAndNumericSubStrs( std::string label,
                                       std::vector<std::string> &alphastrs,
                                       std::vector<std::string> &numericstrs )
{
  SpecUtils::erase_any_character( label, " -_,\t<>/?[]{}\\|!@#$%^&*();:\"'~`+=" );
  SpecUtils::to_lower_ascii( label );

  const string::size_type len = label.length();

  for( size_t i = 0; i < len; )
  {
    string strStr, numStr;
    while( (i<len) && isalpha( label[i] ) )
      strStr += label[i++];
    while( (i<len) && isdigit( label[i] ) )
      numStr += label[i++];

    if( strStr.size() )
      alphastrs.push_back( strStr );
    if( numStr.size() )
      numericstrs.push_back( numStr );
    if( numStr.empty() && strStr.empty() )
      ++i;
  }//for( size_t i = 0; i < len; )

}//void IsotopeNameFilterModel::getAlphaAndNumericSubStrs(...)



void SimpleIsotopeNameFilterModel::filter( const Wt::WString &text )
{
  beginRemoveRows( WModelIndex(), 0, rowCount() );
  m_candidatesElements.clear();
  m_candidatesNuclides.clear();
  endRemoveRows();

  string testTxt = text.toUTF8();
  const int metalevel = determineAndRemoveIsoLevel( testTxt );

  vector<string> alphastrs, numericstrs;
  getAlphaAndNumericSubStrs( testTxt, alphastrs, numericstrs );

  if( text.empty() )
    return;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  vector<const SandiaDecay::Nuclide *> nuclide_suggestions;

  const SandiaDecay::Nuclide *nuc = db->nuclide( text.toUTF8() );


  if( nuc )
  { //First lets see if there is an exact match for a nuclide, if so select it
    m_candidatesNuclides.push_back( nuc );

    const SandiaDecay::Element *element = db->element( nuc->atomicNumber );
    int index = -1;
    if( element )
      index = m_parent->m_elementSelection->findText( element->name );
    if( index >= 0 )
    {
      m_parent->m_elementSelection->setCurrentIndex( index );
      m_parent->makeMassList();
      stringstream massnumberstrm;
      massnumberstrm << nuc->massNumber;
      if( nuc->isomerNumber )
        massnumberstrm << "m";
      if( nuc->isomerNumber >= 2 )
        massnumberstrm << nuc->isomerNumber;
      index = m_parent->m_massSelection->findText( massnumberstrm.str() );
      if( index >= 0 )
      {
        m_parent->m_massSelection->setCurrentIndex( index );
        m_parent->enableAcceptButton();
      }
    }//if( index >= 0 )
  }else
  { //see if user has typed in a legitamate element, if so select that one
    const SandiaDecay::Element *element = db->element( text.toUTF8() );
    int index = -1;
    if( element )
      index = m_parent->m_elementSelection->findText( element->name );
    if( index >= 0 )
    {
      m_parent->m_elementSelection->setCurrentIndex( index );
      m_parent->makeMassList();
    }
  }//if( nuc ) / else

  m_parent->updateSelectedHalfLife();
  
  //suggest based of of alphastrs

  const std::vector<const SandiaDecay::Element *> &elements = db->elements();

  set<const SandiaDecay::Element *> candidate_elements;

  for( const SandiaDecay::Element *el : elements )
  {
    const string name   = SpecUtils::to_lower_ascii_copy( el->name );
    const string symbol = SpecUtils::to_lower_ascii_copy( el->symbol );

    for( string str : alphastrs )
    {
      if( SpecUtils::starts_with( symbol, str.c_str() )
          || SpecUtils::starts_with( name, str.c_str() ) )
        candidate_elements.insert( el );
    }//for( const string &str : alphastrs )
  }//for( const SandiaDecay::Element *el : elements )

  
  
  //Go through and find the candidate Nuclides, but make sure they are
  //  compatible with numericstrs
  if( !numericstrs.empty() )
  {
    for( const SandiaDecay::Element *el : candidate_elements )
    {
      const vector<const SandiaDecay::Nuclide *> nuclides = db->nuclides( el );
      for( const SandiaDecay::Nuclide *nuc : nuclides )
      {
        if( (nuc->halfLife<=0.0) || IsInf(nuc->halfLife) )
          continue;
        
        bool numeric_compat = false;
        for( const string &str : numericstrs )
          numeric_compat |= SpecUtils::contains( std::to_string(nuc->massNumber), str.c_str() );
        
        if( metalevel > 0 )
          numeric_compat = (numeric_compat && (nuc->isomerNumber==metalevel));
        
        if( numeric_compat )
          nuclide_suggestions.push_back( nuc );
      }//for( const SandiaDecay::Nuclide *nuc : nuclides )
    }//for( const SandiaDecay::Element *el : candidate_elements )
  }//if( !numericstrs.empty() )
  
  
  if( numericstrs.empty() )
  {
    //If the user hasnt typed numbers, only suggest elements
    nuclide_suggestions.clear();
  }else
  {
    //If the user has typed numbers, only suggest nuclides
    candidate_elements.clear();
  }
  
  
  //Sort suggestions by Levenshtein distance, so that way the word closest to
  //  what the user has typed in, will be at the top of the list.  If we dont
  //  do this, then if the user types "Ra226" and hits enter, then "Rn226"
  //  will actually be selected.
  const std::string usertxt = text.toUTF8();
    
  //Sort suggested nuclides
  vector<unsigned int> distance( nuclide_suggestions.size() );
  vector<size_t> sort_indices( nuclide_suggestions.size() );
  for( size_t i = 0; i < nuclide_suggestions.size(); ++i )
  {
    sort_indices[i] = i;
    const string &symbol = nuclide_suggestions[i]->symbol;
    distance[i] = SpecUtils::levenshtein_distance( usertxt, symbol );
  }//for( size_t i = 0; i < nuclide_suggestions.size(); ++i )
  std::stable_sort( sort_indices.begin(), sort_indices.end(),
                   index_compare_assend<vector<unsigned int>&>(distance) );
  
  //Sort suggested elements
  vector<const SandiaDecay::Element *> suggest_elements( begin(candidate_elements), end(candidate_elements) );
  
  distance.resize( suggest_elements.size() );
  vector<size_t> element_sort_indices( suggest_elements.size() );
  for( size_t i = 0; i < suggest_elements.size(); ++i )
  {
    element_sort_indices[i] = i;
    const string &symbol = suggest_elements[i]->symbol;
    const string &name = suggest_elements[i]->name;
    distance[i] = std::min( SpecUtils::levenshtein_distance( usertxt, symbol ),
                           SpecUtils::levenshtein_distance( usertxt, name ) );
  }//for( size_t i = 0; i < suggest_elements.size(); ++i )
  std::stable_sort( element_sort_indices.begin(), element_sort_indices.end(),
                   index_compare_assend<vector<unsigned int>&>(distance) );
  
  
  const int nrow = static_cast<int>( candidate_elements.size() + nuclide_suggestions.size() );
  
  beginInsertRows( WModelIndex(), 0, nrow );
  
  m_candidatesElements.resize( suggest_elements.size() );
  for( size_t i = 0; i < suggest_elements.size(); ++i )
    m_candidatesElements[i] = suggest_elements[ element_sort_indices[i] ];
  
  m_candidatesNuclides.resize( nuclide_suggestions.size() );
  for( size_t i = 0; i < nuclide_suggestions.size(); ++i )
    m_candidatesNuclides[i] = nuclide_suggestions[ sort_indices[i] ];
  
  endInsertRows();
}//void SimpleIsotopeNameFilterModel::filter( const Wt::WString &text )


//  Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex & index ) const;


void SimpleIsotopeNameFilterModel::nuclideNameMatcherJs( std::string &js )
{
  /*
   * INLINE_JAVASCRIPT is a macro which allows entry of JavaScript
   * directly in a C++ file - awesome!
   */
  js = INLINE_JAVASCRIPT
  (
    function( edit )
    {
      var value = edit.value;

      return function( suggestion )
      {
        if( !suggestion )
          return value;

        //XXX - TODO - right now the below ignores the case of meta-stable elements
        //             and could possibly cause results like <b><b>M</b>g13<b>m</b></b>
        //             ex. Co60m, U235m, etc
        //XXX - TODO - right now the gamma lines in parenthesis, are also
        //             matched in the regex, this should be avoided

       //replace the special characters the user may have typed in
       value = value.replace(/[.*+?|()\\[\\]{}\\\\$^]/g, "\\$&");

        //Highlight the matching numeric parts of the answer
        var numericregex = new RegExp('(\\d+)');
        var numericmatch = numericregex.exec(value);
        if( numericmatch && numericmatch.length > 1 )
          for( var i = 1; i < numericmatch.length; ++i )
          {
            var regex = new RegExp( '(' + numericmatch[i] + ')', 'gi' );
            suggestion = suggestion.replace( regex, "<b>$1</b>" );
          }

        //Highlight the matching alpha parts of the answer
        var alpharegex = new RegExp('([a-zA-Z]+)');
        var alphamatch = alpharegex.exec(value);
        if( alphamatch && alphamatch.length > 1 )
          for( var i = 1; i < alphamatch.length; ++i )
            if( (alphamatch[i] != 'm' && alphamatch[i] != 'M')
                || !numericmatch || numericmatch.length<2 )
            {
              var regex = new RegExp( '(' + alphamatch[i] + ')', 'gi' );
              suggestion = suggestion.replace( regex, "<b>$1</b>" );
            }

         //For right now we will filter everything server side, so all
         //  suggestions will be a match.
        return { match : true, suggestion : suggestion };
      }
    }
  );
}//void nuclideNameMatcherJs( std::string &js )


void SimpleIsotopeNameFilterModel::replacerJs( std::string &js )
{
  js = INLINE_JAVASCRIPT
  (
      function (edit, suggestionText, suggestionValue)
      {
        edit.value = suggestionValue;

        if (edit.selectionStart)
          edit.selectionStart = edit.selectionEnd = suggestionValue.length;
      }//function (edit, suggestionText, suggestionValue)
   );//js = INLINE_JAVASCRIPT( ... )
}//void replacerJs( std::string &js )
