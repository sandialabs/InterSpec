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

#include <set>
#include <vector>
#include <sstream>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WIOService>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSplitButton>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>


#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PopupDiv.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"

using namespace Wt;
using namespace std;


const int IsotopeSearchByEnergy::sm_xmlSerializationVersion = 0;

const std::string IsotopeSearchByEnergy::sm_medical_category_key    = "isbe-category-medical";
const std::string IsotopeSearchByEnergy::sm_industrial_category_key = "isbe-category-industrial";
const std::string IsotopeSearchByEnergy::sm_snm_category_key        = "isbe-category-snm";
const std::string IsotopeSearchByEnergy::sm_norm_category_key       = "isbe-category-norm";
const std::string IsotopeSearchByEnergy::sm_common_category_key     = "isbe-category-common";
const std::string IsotopeSearchByEnergy::sm_fission_category_key    = "isbe-category-fission";

namespace
{
  const WString ActiveSearchEnergyClass = "ActiveSearchEnergy";
  
  /** Returns peaks that were likely used to enter search energies from. */
  vector<PeakModel::PeakShrdPtr> peaks_searched( const PeakModel * const pmodel,
                                                 const vector<IsotopeSearchByEnergy::SearchEnergy *> &searches )
  {
    if( !pmodel )
      return {};
    
    set<PeakModel::PeakShrdPtr> peaks;
    for( const IsotopeSearchByEnergy::SearchEnergy *search : searches )
    {
      const double energy = search->energy();
      const double window = search->window();
      
      PeakModel::PeakShrdPtr peak = pmodel->nearestPeak( energy );
      if( !peak )
        continue;
      
      const double mean = peak->mean();
      const double fwhm = peak->fwhm();
      const double diff = fabs( mean - energy );
      
      if( (diff < 0.5*fwhm) && (diff < window) )
        peaks.insert( peak );
    }//for( SearchEnergy *search : searches() )
    
    vector<PeakModel::PeakShrdPtr> answer( begin(peaks), end(peaks) );
    std::sort( begin(answer), end(answer),  &PeakDef::lessThanByMeanShrdPtr );
    
    return answer;
  }//peaks_searched()
}//namespace


void IsotopeSearchByEnergy::SearchEnergy::emitRemove()
{
  m_remove.emit();
}


void IsotopeSearchByEnergy::SearchEnergy::emitChanged()
{
  m_changed.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitEnter()
{
  m_enter.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitGotFocus()
{
  m_focus.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitAddAnother()
{
  m_addAnother.emit();
}

IsotopeSearchByEnergy::SearchEnergy::SearchEnergy( Wt::WContainerWidget *p )
: WContainerWidget( p ),
  m_removeIcn( 0 ),
  m_addAnotherIcn( 0 ),
  m_energy( 0 ),
  m_window( 0 )
{
  addStyleClass( "SearchEnergy" );
  
  WLabel *label = new WLabel( "Energy", this );
  label->addStyleClass( "SearchEnergyLabel" );
  
  m_energy = new NativeFloatSpinBox( this );
  label->setBuddy( m_energy );
  m_energy->setMinimum( 0.0f );
  m_energy->setMaximum( 1000000.0f );
  m_energy->enterPressed().connect( this, &SearchEnergy::emitEnter );
  
  label = new WLabel( "+/-", this );
  label->addStyleClass( "SearchEnergyWindowLabel" );
  
  m_window = new NativeFloatSpinBox( this );
  label->setBuddy( m_window );
  
  m_window->setMinimum( 0.0f );
  m_window->setMaximum( 1000000.0f );
  m_window->setValue( 10.0f );
  m_window->enterPressed().connect( this, &SearchEnergy::emitEnter );
  
  label = new WLabel( "keV", this );
  label->addStyleClass( "KeVLabel" );
  
  m_energy->valueChanged().connect( this, &SearchEnergy::emitChanged );
  m_window->valueChanged().connect( this, &SearchEnergy::emitChanged );
  //m_energy->keyPressed().connect( this, &SearchEnergy::emitChanged );
  
  m_energy->focussed().connect( this, &SearchEnergy::emitGotFocus );
  m_window->focussed().connect( this, &SearchEnergy::emitGotFocus );
  clicked().connect( this, &SearchEnergy::emitGotFocus );
    
  // Add a spacer incase we get wider than we could reasonable want to grow the text inputs, so
  //  there will be a space between the window input and the add/remove buttons.
  WContainerWidget *spacer = new WContainerWidget( this );
  spacer->addStyleClass( "SearchEnergySpacer" );
  
  m_removeIcn = new WContainerWidget( this ); //needed or else button wont show up
  m_removeIcn->setStyleClass( "DeleteSearchEnergy Wt-icon" );
  m_removeIcn->clicked().connect( this, &SearchEnergy::emitRemove );
  m_removeIcn->clicked().preventPropagation();  //make it so we wont emit gotFocus()
  
  m_addAnotherIcn = new WContainerWidget( this ); //needed or else button wont show up
  m_addAnotherIcn->setStyleClass( "AddSearchEnergy Wt-icon" );
  m_addAnotherIcn->clicked().connect( this, &SearchEnergy::emitAddAnother );
  m_addAnotherIcn->clicked().preventPropagation(); //make it so we wont emit gotFocus(), which would
                                                   // keep new search energy from getting focus due
                                                   // order of handling signal callbacks
}//SearchEnergy constructor


void IsotopeSearchByEnergy::SearchEnergy::enableAddAnother()
{
  m_addAnotherIcn->enable();
}//void enableAddAnother()


void IsotopeSearchByEnergy::SearchEnergy::disableAddAnother()
{
  m_addAnotherIcn->disable();
}//void disableAddAnother()


void IsotopeSearchByEnergy::SearchEnergy::enableRemove()
{
  m_removeIcn->enable();
}


void IsotopeSearchByEnergy::SearchEnergy::disableRemove()
{
  m_removeIcn->disable();
}


double IsotopeSearchByEnergy::SearchEnergy::energy() const
{
  if( m_energy->validate() == WValidator::Valid )
    return m_energy->value();
  return 0.0;
}

void IsotopeSearchByEnergy::SearchEnergy::setEnergy( double energy )
{
  m_energy->setValue( energy );
}

void IsotopeSearchByEnergy::SearchEnergy::setWindow( double window )
{
  window = floor( 100.0*window + 0.5 ) / 100.0;
  m_window->setValue( window );
}

double IsotopeSearchByEnergy::SearchEnergy::window() const
{
  if( m_window->validate() == WValidator::Valid )
    return m_window->value();
  return 0.0;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::enter()
{
  return m_enter;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::changed()
{
  return m_changed;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::gotFocus()
{
  return m_focus;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::remove()
{
  return m_remove;
}


Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::addAnother()
{
  return m_addAnother;
}


IsotopeSearchByEnergy::IsotopeSearchByEnergy( InterSpec *viewer,
                                              D3SpectrumDisplayDiv *chart,
                                              Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_viewer( viewer ),
  m_chart( chart ),
  m_searchConditionsColumn( nullptr ),
  m_searchEnergies( NULL ),
  m_assignBtnRow( nullptr ),
  m_clearRefLines( nullptr ),
  m_assignPeakToSelected( nullptr ),
  m_currentSearch( 0 ),
  m_searching( NULL ),
  m_results( NULL ),
  m_minBranchRatioDiv( nullptr ),
  m_minBranchRatio( NULL ),
  m_minHalfLiveDiv( nullptr ),
  m_minHalfLife( NULL ),
  m_model( NULL ),
  m_search_categories{},
  m_search_category_select( nullptr ),
  m_nextSearchEnergy( 0 ),
  m_minBr( 0.0 ), m_minHl( 6000.0 * PhysicalUnits::second ),
  m_undo_redo_sentry{},
  m_state{},
  m_selected_row( -1 )
{
  wApp->useStyleSheet( "InterSpec_resources/IsotopeSearchByEnergy.css" );
  
  addStyleClass( "IsotopeSearchByEnergy" );
    
  viewer->useMessageResourceBundle( "IsotopeSearchByEnergy" );
  
  shared_ptr<void> undo_sentry = getDisableUndoRedoSentry();
  
  m_searchConditionsColumn = new WContainerWidget( this );
  m_searchConditionsColumn->setStyleClass( "IsotopeSearchConditions" );
  
  m_searchEnergies = new WContainerWidget( m_searchConditionsColumn );
  m_searchEnergies->setStyleClass( "IsotopeSearchEnergies" );
  
  m_assignBtnRow = new WContainerWidget( m_searchConditionsColumn );
  m_assignBtnRow->addStyleClass( "AssignToSelectedRow" );
  
  m_clearRefLines = new WPushButton( WString::tr("isbe-clear-ref"), m_assignBtnRow );
  m_clearRefLines->addStyleClass( "LightButton" );
  m_clearRefLines->clicked().connect( this, &IsotopeSearchByEnergy::clearSelectionAndRefLines );
  m_clearRefLines->hide();
  
  
  m_assignPeakToSelected = new WSplitButton( "&nbsp;", m_assignBtnRow ); //Space is needed so Wt will add the ".with-label" style class
  m_assignPeakToSelected->addStyleClass( "LightButton" );
  m_assignPeakToSelected->actionButton()->addStyleClass( "LightButton" );
  m_assignPeakToSelected->dropDownButton()->addStyleClass( "LightButton" );
  m_assignPeakToSelected->actionButton()->clicked().connect( this, &IsotopeSearchByEnergy::assignSearchedOnPeaksToSelectedNuclide );
  m_assignPeakToSelected->hide();
  
  WPopupMenu *assignPeakMenu = nullptr;
  if( m_viewer->isMobile() )
    assignPeakMenu = new WPopupMenu();
  else
    assignPeakMenu = new PopupDivMenu( nullptr, PopupDivMenu::MenuType::TransientMenu);
  m_assignPeakToSelected->setMenu( assignPeakMenu );
  // We will add relevant menu items to this button when the time comes
  
  m_search_category_select = new WComboBox( m_searchConditionsColumn );
  m_search_category_select->setStyleClass( "IsotopeSourceTypes" );
  m_search_category_select->activated().connect( boost::bind( &IsotopeSearchByEnergy::categoryChanged, this, true ) );
  
  WContainerWidget *searchOptions = new WContainerWidget( m_searchConditionsColumn );
  searchOptions->setStyleClass( "IsotopeSearchMinimums" );

  
  auto helpBtn = new WContainerWidget( searchOptions );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "nuclide-search-dialog" ) );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  m_minBranchRatioDiv = new WContainerWidget( searchOptions );
  m_minBranchRatioDiv->setHiddenKeepsGeometry( true );
  WLabel *label = new WLabel( WString::tr("isbe-min-br"), m_minBranchRatioDiv );
//  HelpSystem::attachToolTipOn( label,"Toggle or type minimum branching ratio.", showToolTips , HelpSystem::ToolTipPosition::Top);

  m_minBranchRatio = new NativeFloatSpinBox( m_minBranchRatioDiv );
  HelpSystem::attachToolTipOn( m_minBranchRatio, WString::tr("isbe-tt-min-br"),
                              showToolTips , HelpSystem::ToolTipPosition::Top);
  
  m_minBranchRatio->setValue( m_minBr );
  m_minBranchRatio->setRange( 0.0f, 1.0f );
  m_minBranchRatio->setSingleStep( 0.1f );
  label->setBuddy( m_minBranchRatio );
  
  m_minHalfLiveDiv = new WContainerWidget( searchOptions );
  m_minHalfLiveDiv->setHiddenKeepsGeometry( true );
  label = new WLabel( WString::tr("isbe-min-hl"), m_minHalfLiveDiv ); //"Min. T\xc2\xbd"
 
  m_minHalfLife = new WLineEdit( "6000 s", m_minHalfLiveDiv );
  m_minHl = 6000.0 * PhysicalUnits::second;
  
  m_minHalfLife->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_minHalfLife->setAttributeValue( "autocorrect", "off" );
  m_minHalfLife->setAttributeValue( "spellcheck", "off" );
#endif
    
  HelpSystem::attachToolTipOn( {label,m_minHalfLife}, WString::tr("isbe-tt-min-hl"),
                              showToolTips , HelpSystem::ToolTipPosition::Top );

  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationRegex(), this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_minHalfLife->setValidator(validator);
  
  label->setBuddy( m_minHalfLife );

  m_minBranchRatio->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->valueChanged().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->enterPressed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->blurred().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );

  m_minHalfLife->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minHalfLife->enterPressed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minHalfLife->blurred().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  
  m_model = new IsotopeSearchByEnergyModel( this );
  
  // Even though we are using a CSS flex layout to control the size of the m_results table, the JS
  //  wtResize function will be called on this RowStretchTreeView.
  m_results = new RowStretchTreeView( this );
  m_results->setRootIsDecorated(false); //makes the tree look like a table! :)
  
  m_results->setModel( m_model );
  m_results->addStyleClass( "IsotopeSearchResultTable ToolTabSection" );
  m_results->setAlternatingRowColors( true );
  m_results->sortByColumn( IsotopeSearchByEnergyModel::Column::ProfileDistance, Wt::DescendingOrder );
  
#if( InterSpec_PHONE_ROTATE_FOR_TABS )
  setResultTableColumnWidths( false );
#else
    for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
        col < IsotopeSearchByEnergyModel::NumColumns;
        col = IsotopeSearchByEnergyModel::Column(col+1) )
    {
      m_results->setColumnHidden( col, false );
      m_results->setSortingEnabled( col, true );

      switch( col )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::Energy:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::Distance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::BranchRatio:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::NumColumns:
        break;
      }//switch( col )
    }//for( loop over peak columns )
#endif //InterSpec_PHONE_ROTATE_FOR_TABS / else
  
  m_results->setSelectionMode( Wt::SingleSelection );
  m_results->setSelectionBehavior( Wt::SelectRows );
  m_results->selectionChanged().connect( this, &IsotopeSearchByEnergy::resultSelectionChanged );
  
  m_searching = new WText( WString::tr("isbe-searching"), this );
  m_searching->addStyleClass( "IsotopeSearchInProgress" );
  m_searching->setInline( false );
  m_searching->hide();
  
  
  // Add in one non-removable search energy
  SearchEnergy *enrgy = new SearchEnergy( m_searchEnergies );
  enrgy->enter().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->addAnother().connect( this, &IsotopeSearchByEnergy::addSearchEnergy );
  enrgy->gotFocus().connect( boost::bind( &IsotopeSearchByEnergy::searchEnergyRecievedFocus,
                                        this, enrgy ) );
  enrgy->remove().connect( boost::bind( &IsotopeSearchByEnergy::removeSearchEnergy, this, enrgy) );
  enrgy->addStyleClass( ActiveSearchEnergyClass );
  enrgy->disableRemove();
  enrgy->changed().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );

  // We could initialize filter types here, but instead we'll wait until first render, to slightly
  //  help out initial app render
  //initFilterTypes();
    
  // Previous to 20230928 we called `minBrOrHlChanged()` here, but this had the side-effect of
  //  also clearing the ReferencePhotopeak nuclide, which if we are restoring a state, we
  //  dont want to do.  So instead we will pre-initialize the mapping from energy to
  //  nuclide explicitly, in a separate thread, once the `render()` function is called
  //  the first time (which is immediately after this since we preload the tool-tabs, but whatever).
  //minBrOrHlChanged();
    
  // During normal desktop construction of this widget, the ReferencePhotopeakDisplay is still
  //  nullptr at this point, but we'll check here, JIC
  updateClearSelectionButton();
  
  ReferencePhotopeakDisplay *display = m_viewer ? m_viewer->referenceLinesWidget() : nullptr;
  if( display )
  {
    m_refLineUpdateConnection
      = display->displayingNuclide().connect( this, &IsotopeSearchByEnergy::handleRefLinesUpdated );
    m_refLineClearConnection
      = display->nuclidesCleared().connect( this, &IsotopeSearchByEnergy::handleRefLinesUpdated );
  }//if( display )
}//IsotopeSearchByEnergy constructor


void IsotopeSearchByEnergy::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WContainerWidget::render( flags );
  
  if( flags.testFlag(Wt::RenderFlag::RenderFull) )
  {
    initFilterTypes();
    
    // Initialize the mapping from energies to nuclides when we render this widget the first
    //  time, so this way it will be ever so slightly quicker when the user does the first search.
    const double minHl = m_minHl, minBr = m_minBr;
    WServer::instance()->ioService().boost::asio::io_service::post( [minHl,minBr](){
      EnergyToNuclideServer::setLowerLimits( minHl, minBr );
      EnergyToNuclideServer::energyToNuclide();
    } );
  }
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void IsotopeSearchByEnergy::searchEnergyRecievedFocus( SearchEnergy *enrgy )
{
  std::vector<SearchEnergy *> searchv = searches();
  for( size_t index = 0; index < searchv.size(); ++index )
  {
    if( searchv[index] == enrgy )
    {
      m_nextSearchEnergy = index;
      if( !enrgy->hasStyleClass( ActiveSearchEnergyClass ) )
        enrgy->addStyleClass( ActiveSearchEnergyClass );
    }else
    {
      searchv[index]->removeStyleClass( ActiveSearchEnergyClass );
    }
  }//for( WWidget *kid : kids )
}//void searchEnergyRecievedFocus( SearchEnergy *enrgy )


IsotopeSearchByEnergy::SearchEnergy *IsotopeSearchByEnergy::addNewSearchEnergy()
{
  std::vector<SearchEnergy *> searchv = searches();
  for( SearchEnergy *enrgy : searchv )
    enrgy->disableAddAnother();
  
  SearchEnergy *enrgy = new SearchEnergy( m_searchEnergies );
  enrgy->changed().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->enter().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->remove().connect(
                      boost::bind( &IsotopeSearchByEnergy::removeSearchEnergy,
                                   this, enrgy) );
  enrgy->addAnother().connect( this, &IsotopeSearchByEnergy::addSearchEnergy );
  enrgy->gotFocus().connect(
                boost::bind( &IsotopeSearchByEnergy::searchEnergyRecievedFocus,
                             this, enrgy ) );
//  m_search->changed().connect( this,
//                          &IsotopeSearchByEnergy::loadSearchEnergiesToClient );
//  m_search->enable();
  
  searchv = searches();
  if( searchv.size() > 1 )
    searchv[0]->enableRemove();
  
  for( size_t i = 0; i < searchv.size(); ++i )
  {
    if( searchv[i]->energy() < 0.1 )
    {
      searchEnergyRecievedFocus( searchv[i] );
      // Note: if there are more callbacks queued, like gotFocus(), then they may happen after here, so we
      //    need to either be careful of the ordering of callbacks, or use something like a WTimer here.
      //WTimer::singleShot( 1, boost::bind( &IsotopeSearchByEnergy::searchEnergyRecievedFocus, this, searchv[i]) );
      break;
    }//if( searchv[i]->energy() < 0.01 )
  }//for( size_t i = 0; i < searchv.size(); ++i )
  
  return enrgy;
}//SearchEnergy *addNewSearchEnergy()


void IsotopeSearchByEnergy::addSearchEnergy()
{
  addNewSearchEnergy();
  
  //For sake of undo/redo
  addUndoRedoPoint();
}//void addSearchEnergy()


vector<IsotopeSearchByEnergy::SearchEnergy *> IsotopeSearchByEnergy::searches()
{
  vector<SearchEnergy *> searchEnergies;
  
  const vector<WWidget *> &children = m_searchEnergies->children();
  for( size_t index = 0; index < children.size(); ++index )
  {
    SearchEnergy *ww = dynamic_cast<SearchEnergy *>( children[index] );
    if( ww )
      searchEnergies.push_back( ww );
  }
  
  return searchEnergies;
}//vector<SearchEnergy *> searches()


#if( InterSpec_PHONE_ROTATE_FOR_TABS )
void IsotopeSearchByEnergy::setResultTableColumnWidths( const bool narrow )
{
  if( narrow )
  {
    m_results->setLineHeight( 13 );
    m_results->setHeaderHeight( 15 );
    
    for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
        col < IsotopeSearchByEnergyModel::NumColumns;
        col = IsotopeSearchByEnergyModel::Column(col+1) )
    {
      m_results->setSortingEnabled( col, true );
      
      switch( col )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::Energy:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::Distance:
          m_results->setColumnWidth( col, WLength(35,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::BranchRatio:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
          break;
          
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::NumColumns:
        break;
      }//switch( col )
    }//for( loop over peak columns )
  }else
  {
    m_results->setLineHeight( 20 );
    m_results->setHeaderHeight( 20 );
    
    for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
        col < IsotopeSearchByEnergyModel::NumColumns;
        col = IsotopeSearchByEnergyModel::Column(col+1) )
    {
      m_results->setColumnHidden( col, false );
      m_results->setSortingEnabled( col, true );
      
      switch( col )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::Energy:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::Distance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::BranchRatio:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::NumColumns:
          break;
      }//switch( col )
    }//for( loop over peak columns )
  }//if( narrow ) / else
}//void setResultTableColumnWidths( const bool narrow )


void IsotopeSearchByEnergy::setNarrowPhoneLayout( const bool narrow )
{
  //WWidget *p = m_results->parent();
  //if( p )
  //  p->removeChild( m_results );
  
  if( narrow )
  {
    addStyleClass( "NarrowNuclideSearch" );
    m_searchConditionsColumn->insertBefore( m_results, m_assignBtnRow );
    m_results->setLineHeight( 13 );
    m_results->setHeaderHeight( 15 );
    
    for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
        col < IsotopeSearchByEnergyModel::NumColumns;
        col = IsotopeSearchByEnergyModel::Column(col+1) )
    {
      switch( col )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::Energy:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::Distance:
          m_results->setColumnWidth( col, WLength(35,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::BranchRatio:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
        break;
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          m_results->setColumnWidth( col, WLength(50,WLength::Pixel) );
          m_results->setColumnHidden( col, false );
          break;
          
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
          m_results->setColumnHidden( col, true );
        break;
          
        case IsotopeSearchByEnergyModel::NumColumns:
        break;
      }//switch( col )
    }//for( loop over peak columns )
  }else
  {
    removeStyleClass( "NarrowNuclideSearch" );
    addWidget( m_results );
    m_results->setLineHeight( 20 );
    m_results->setHeaderHeight( 20 );
    
    for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
        col < IsotopeSearchByEnergyModel::NumColumns;
        col = IsotopeSearchByEnergyModel::Column(col+1) )
    {
      m_results->setColumnHidden( col, false );
      m_results->setSortingEnabled( col, true );

      switch( col )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::Energy:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::Distance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::BranchRatio:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
          break;
          
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
        break;
          
        case IsotopeSearchByEnergyModel::NumColumns:
        break;
      }//switch( col )
    }//for( loop over peak columns )
  }
}//void setNarrowPhoneLayout( const bool narrow )
#endif //InterSpec_PHONE_ROTATE_FOR_TABS


const vector<IsotopeSearchByEnergy::NucSearchCategory> &IsotopeSearchByEnergy::search_categories()
{
  if( m_search_categories.empty() )
    initFilterTypes();
  
  return m_search_categories;
}//const vector<NucSearchCategory> search_categories() const


const IsotopeSearchByEnergy::NucSearchCategory &IsotopeSearchByEnergy::get_category_info( const string &cat_key, 
                                                      const vector<NucSearchCategory> &categories )
{
  assert( !categories.empty() );
  
  if( categories.empty() )
    throw runtime_error( "Categories not inited." );
  
  const NucSearchCategory *cat = nullptr;
  for( const NucSearchCategory &nsc : categories )
  {
    string key = nsc.m_name.key();
    if( key.empty() )
      key = nsc.m_name.toUTF8();
    
    if( key == cat_key )
    {
      cat = &nsc;
      break;
    }
  }//for( loop over categories )
  
  if( !cat )
  {
    vector<string> fields;
    SpecUtils::split( fields, cat_key, "-" ); //cat_key is like "isbe-category-medical"
    throw runtime_error( "Couldnt find " + (fields.size()>2 ? fields[2] : cat_key) + " category." );
  }//if( !cat )
  
  return *cat;
}//const NucSearchCategory &get_category_info( const string &cat_key )


bool IsotopeSearchByEnergy::is_in_category( const SandiaDecay::Nuclide *nuc,
                           const std::string &cat_name,
                           const std::vector<NucSearchCategory> &category )
{
  const NucSearchCategory &cat = get_category_info( cat_name, category );
  const auto pos = std::find( begin(cat.m_specific_nuclides), end(cat.m_specific_nuclides), nuc );
  return ((cat.m_nuclides && cat.m_specific_nuclides.empty())
          || (pos != end(cat.m_specific_nuclides)));
}


bool IsotopeSearchByEnergy::is_in_category( const SandiaDecay::Element *el,
                           const std::string &cat_name,
                           const std::vector<NucSearchCategory> &category )
{
  const NucSearchCategory &cat = get_category_info( cat_name, category );
  const auto pos = std::find( begin(cat.m_specific_elements), end(cat.m_specific_elements), el );
  return ((cat.m_fluorescence_xrays && cat.m_specific_elements.empty())
          || (pos != end(cat.m_specific_elements)));
}


bool IsotopeSearchByEnergy::is_in_category( const ReactionGamma::Reaction *rctn,
                           const std::string &cat_name,
                           const std::vector<NucSearchCategory> &category )
{
  const NucSearchCategory &cat = get_category_info( cat_name, category );
  const auto pos = std::find( begin(cat.m_specific_reactions), end(cat.m_specific_reactions), rctn );
  return ((cat.m_reactions && cat.m_specific_reactions.empty())
          || (pos != end(cat.m_specific_reactions)));
}


std::shared_ptr<void> IsotopeSearchByEnergy::getDisableUndoRedoSentry()
{
  shared_ptr<void> answer = m_undo_redo_sentry.lock();
  if( answer )
    return answer;
  
  int *dummy = new int(0);
  auto deleter = []( void *obj ){
    int *sentry = (int *)obj;
    if( sentry )
      delete sentry;
  };
  
  answer = shared_ptr<void>( dummy, deleter );
  m_undo_redo_sentry = answer;
  return answer;
}//std::shared_ptr<void> getDisableUndoRedoSentry()


void IsotopeSearchByEnergy::loadSearchEnergiesToClient()
{
  // We get here when user changes the current tools tab to this tool
  
  vector<pair<double,double>> searchRegions;
  
  for( auto sw : searches() )
  {
    if( sw->energy() > 0.1 )
      searchRegions.push_back( make_pair(sw->energy(), sw->window()) );
  }
  
  m_chart->setSearchEnergies( searchRegions );
  
  // Update the "Clear Ref. Lines" button
  updateClearSelectionButton();
  
  // Listen for changes to display reference lines, and update "Clear Ref. Lines" button for those
  ReferencePhotopeakDisplay *display = m_viewer ? m_viewer->referenceLinesWidget() : nullptr;
  if( display )
  {
    if( !m_refLineUpdateConnection.connected() )
      m_refLineUpdateConnection
        = display->displayingNuclide().connect( this, &IsotopeSearchByEnergy::handleRefLinesUpdated );
    if( !m_refLineClearConnection.connected() )
      m_refLineClearConnection
        = display->nuclidesCleared().connect( this, &IsotopeSearchByEnergy::handleRefLinesUpdated );
  }//if( display )
}//void loadSearchEnergiesToClient()


IsotopeSearchByEnergyModel::Column IsotopeSearchByEnergyModel::sortColumn() const
{
  return m_sortColumn;
}

Wt::SortOrder IsotopeSearchByEnergyModel::sortOrder() const
{
  return m_sortOrder;
}


void IsotopeSearchByEnergy::clearSearchEnergiesOnClient()
{
  // Called when the user clicks off this tab.
  m_chart->setSearchEnergies( vector<pair<double,double>>() );
  
  if( m_refLineUpdateConnection.connected() )
    m_refLineUpdateConnection.disconnect();
  
  if( m_refLineClearConnection.connected() )
    m_refLineClearConnection.disconnect();
}//void clearSearchEnergiesOnClient()


void IsotopeSearchByEnergy::setNextSearchEnergy( double energy, double sigma )
{
  const vector<SearchEnergy *> searchW = searches();
  if( searchW.empty() )  //shouldnt ever happen, bu JIC
    return;
  
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  
  
  for( size_t index = 0; index < searchW.size(); ++index )
  {
    searchW[index]->removeStyleClass( ActiveSearchEnergyClass );
    if( fabs(searchW[index]->energy()-energy) < 0.01*PhysicalUnits::keV )
      m_nextSearchEnergy = index;
  }
  
  searchW[m_nextSearchEnergy]->setEnergy( energy );
  if( sigma > 0.0 )
    searchW[m_nextSearchEnergy]->setWindow( sigma );
  
  ++m_nextSearchEnergy;
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  searchW[m_nextSearchEnergy]->addStyleClass( ActiveSearchEnergyClass );
  
  startSearch( false );
}//void setNextSearchEnergy( double energy )


void IsotopeSearchByEnergy::removeSearchEnergy( IsotopeSearchByEnergy::SearchEnergy *energy )
{
  vector<SearchEnergy *> searchW = searches();
  
  if( searchW.size() < 2 )
    return;
  
  vector<SearchEnergy *>::iterator pos
                     = std::find( searchW.begin(), searchW.end(), energy );
  
  if( pos == searchW.end() )  //shouldnt ever happen, but JIC
    return;
  
  const size_t index = pos - searchW.begin();
  if( index < m_nextSearchEnergy )
    --m_nextSearchEnergy;
  
  searchW.erase( pos );
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  
  for( SearchEnergy *s : searchW )
    s->removeStyleClass( ActiveSearchEnergyClass );
  
  searchW[m_nextSearchEnergy]->addStyleClass( ActiveSearchEnergyClass );
  searchW.back()->enableAddAnother();
  if( searchW.size() == 1 )
    searchW[0]->disableRemove();
  else if( !searchW.empty() )
    searchW[0]->enableRemove();
  
  delete energy;
  startSearch( false );
}//void removeSearchEnergy( SearchEnergy *energy )


void IsotopeSearchByEnergy::init_category_info( std::vector<NucSearchCategory> &results )
{
  if( !results.empty() )
  {
    cout << "IsotopeSearchByEnergy: filters already initiated." << endl;
    return;
  }
  
  auto append_categories = [&results]( const std::string &xml_file_name ) -> int {
#ifdef _WIN32
    const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(xml_file_name);
    std::ifstream input( wfilename.c_str() );
#else
    std::ifstream input( xml_file_name.c_str() );
#endif
      
    if( !input.is_open() )
      return -1;
    
    int num_loaded = 0;
    
    try
    {
      rapidxml::file<char> input_file( input );
      
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_trim_whitespace>( input_file.data() );
      
      const auto categories_node = doc.first_node( "NucSearchCategories" );
      if( !categories_node )
        throw runtime_error( "No NucSearchCategories found." );
      
      XML_FOREACH_CHILD( category_node, categories_node, "NucSearchCategory" )
      {
        try
        {
          NucSearchCategory cat;
          cat.deSerialize( category_node );
          results.push_back( std::move(cat) );
          
          num_loaded += 1;
        }catch( std::exception &e )
        {
          cerr << "Failed to deserialize NucSearchCategory node: " << e.what() << endl;
        }
      }//for( loop over NucSearchCategory nodes )
    }catch( rapidxml::parse_error &e )
    {
      string msg = "Failed parsing '" + xml_file_name + "': " + string(e.what());
      const char * const position = e.where<char>();
      if( position && *position )
      {
        const char *end_pos = position;
        for( size_t i = 0; (*end_pos) && (i < 80); ++i )
          end_pos += 1;
        msg += "\n" + std::string(position, end_pos);
      }//if( position )
      
      cerr << msg << endl;
    }catch( std::exception &e )
    {
      cerr << "Failed parsing '" << xml_file_name << "': " << e.what() << endl;
    }//try / catch to parse XML and load them
    
    return num_loaded;
  };//append_categories(...)
  
  {// Begin add a default search category, so even if we dont find default XML, things wont be completely borked
    NucSearchCategory def_category;
    def_category.m_name = WString::tr("isbe-category-nuc-xray");
    def_category.m_description = WString::tr("isbe-category-nuc-xray-desc");
    def_category.m_minBr = 0.0;
    def_category.m_minHl = 6000*PhysicalUnits::second;
    def_category.m_nuclides = true;
    def_category.m_fluorescence_xrays = true;
    def_category.m_reactions = false;
    def_category.m_alphas = false;
    def_category.m_beta_endpoint = false;
    def_category.m_no_progeny = false;
    
    results.push_back( std::move(def_category) );
  }// End add a default search category
  
  
  // Load the default category info
  const string static_data_dir = InterSpec::staticDataDirectory();
  const string def_xml = SpecUtils::append_path(static_data_dir, "NuclideSearchCatagories.xml");
  
  if( append_categories(def_xml) <= 0 )
    throw runtime_error( "Failed to load default NuclideSearchCatagories.xml" );

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP || BUILD_AS_UNIT_TEST_SUITE )
  // Load user specific category info
  const string user_data_dir = InterSpec::writableDataDirectory();
  const std::string user_xml = SpecUtils::append_path(def_xml, "NuclideSearchCatagories.xml");
  if( append_categories(user_xml) == 0 )
  {
    if( SpecUtils::is_file(user_xml) )
      throw runtime_error( "Error loading user NuclideSearchCatagories.xml" );
  }
#endif
}//static void init_category_info( std::vector<NucSearchCategory> &results );


void IsotopeSearchByEnergy::initFilterTypes()
{
  if( !m_search_categories.empty() )
  {
    cout << "IsotopeSearchByEnergy: filters already initiated." << endl;
    return;
  }
  
  try
  {
    init_category_info( m_search_categories );
  }catch( std::exception &e )
  {
    cerr << "IsotopeSearchByEnergy::initFilterTypes(): " << e.what() << endl;
  }
  
  for( const NucSearchCategory &cat : m_search_categories )
    m_search_category_select->addItem( cat.m_name );
    
  m_search_category_select->setCurrentIndex( 0 );
  
  categoryChanged(false);
}//void initFilterTypes()


void IsotopeSearchByEnergy::minBrOrHlChanged()
{
  if( m_minBranchRatio->validate() == WValidator::Valid )
    m_minBr = m_minBranchRatio->value();
  else
    m_minBranchRatio->setValue( m_minBr );
  
  try
  {
    const string hltxt = m_minHalfLife->valueText().toUTF8();
    m_minHl = PhysicalUnitsLocalized::stringToTimeDuration( hltxt );
  }catch(...)
  {
    m_minHalfLife->setText( PhysicalUnitsLocalized::printToBestTimeUnits(m_minHl,2) );
  }//try / catch

  startSearch( true );
}//void IsotopeSearchByEnergy::minBrOrHlChanged()


void IsotopeSearchByEnergy::categoryChanged( const bool update_results )
{
  assert( !m_search_categories.empty() );
  if( m_search_categories.empty() )
    return;
  
  int current_index = m_search_category_select->currentIndex();
  assert( (current_index >= 0) && (current_index < static_cast<int>(m_search_categories.size())) );
  
  if( (current_index < 0) || (current_index > static_cast<int>(m_search_categories.size())) )
  {
    current_index = 0;
    m_search_category_select->setCurrentIndex( current_index );
  }
  
  addUndoRedoPoint();
  
  const NucSearchCategory &cat = m_search_categories[current_index];
  m_search_category_select->setToolTip( cat.m_description );
  m_minHl = cat.m_minHl;
  m_minHalfLife->setValueText( PhysicalUnitsLocalized::printToBestTimeUnits(m_minHl, 1) );
  
  const bool refreshBr = (m_minBr != cat.m_minBr);
  m_minBr = cat.m_minBr;
  m_minBranchRatio->setValue( m_minBr );
  
  m_minBranchRatioDiv->setHidden( !cat.m_nuclides );
  m_minHalfLiveDiv->setHidden( !cat.m_nuclides );
  
  if( update_results )
  {
    // User changed things here, so if applicable, get rid of showing alpha/beta
    //  ref. lines.  Its possible the user explicitly
    ReferencePhotopeakDisplay *display = m_viewer->referenceLinesWidget();
    if( display )
    {
      UndoRedoManager::BlockUndoRedoInserts undo_blocker;
      if( !display->showingGammaLines() )
        display->setShowGammaLines( cat.m_nuclides || cat.m_reactions );
      if( !display->showingXrayLines() )
        display->setShowXrayLines( cat.m_nuclides || cat.m_fluorescence_xrays );
      display->setShowAlphaLines( cat.m_alphas );
      display->setShowBetaLines( cat.m_beta_endpoint );
    }//if( display )
    
    startSearch( refreshBr );
  }
}//void categoryChanged()


void IsotopeSearchByEnergy::resultSelectionChanged()
{
  if( !m_viewer )
    return;
  
  // This is a bit of a hack, but if the tool tabs are not visible, then the
  //  ReferencePhotopeakDisplay pointer will be nullptr, so we will create a reference gamma lines
  //  window - havent tested this on phones, as of 20201030
  ReferencePhotopeakDisplay *display = m_viewer->referenceLinesWidget();
  if( !display && !m_viewer->toolTabsVisible() )
  {
    m_viewer->showGammaLinesWindow();
    if( m_viewer->isPhone() )
      m_viewer->closeGammaLinesWindow();
    display = m_viewer->referenceLinesWidget();
  }
  
  shared_ptr<void> ref_line_undo_sentry = display ? display->getDisableUndoRedoSentry() : nullptr;
  
  string orig_state_xml;
  if( display )
    display->serialize( orig_state_xml );
  
  if( display )
  {
    int current_index = m_search_category_select->currentIndex();
    assert( (current_index >= 0) && (current_index < static_cast<int>(m_search_categories.size())) );
    
    if( (current_index < 0) || (current_index > static_cast<int>(m_search_categories.size())) )
      current_index = 0;
    assert( current_index < static_cast<int>(m_search_categories.size()) );
    const NucSearchCategory &cat = m_search_categories[current_index];
    
    UndoRedoManager::BlockUndoRedoInserts undo_blocker;
    if( !display->showingGammaLines() )
      display->setShowGammaLines( cat.m_nuclides || cat.m_reactions );
    if( !display->showingXrayLines() )
      display->setShowXrayLines( cat.m_nuclides || cat.m_fluorescence_xrays );
    display->setShowAlphaLines( cat.m_alphas );
    display->setShowBetaLines( cat.m_beta_endpoint );
  }//if( display )
  
  const int orig_selected_row = m_selected_row;
  
  const int nPeaksSearched = numSearchEnergiesOnPeaks();
 
  UndoRedoManager *undoManager = m_viewer->undoRedoManager();
  auto undo = [orig_state_xml, orig_selected_row](){
    InterSpec *viewer = InterSpec::instance();
    ReferencePhotopeakDisplay *display = viewer ? viewer->referenceLinesWidget() : nullptr;
    IsotopeSearchByEnergy *search = viewer ? viewer->nuclideSearch() : nullptr;
    assert( display && search );
    if( !display || !search )
      return;
    
    shared_ptr<void> ref_line_undo_sentry = display->getDisableUndoRedoSentry();
    shared_ptr<void> search_undo_sentry = search->getDisableUndoRedoSentry();
    
    search->m_selected_row = orig_selected_row;
    if( orig_selected_row >= 0 )
      search->m_results->setSelectedIndexes( {search->m_model->index(orig_selected_row, 0)} );
    else
      search->m_results->setSelectedIndexes( {} );
    search->resultSelectionChanged();
    
    try
    {
      string orig_state_xml_copy = orig_state_xml;
      display->deSerialize( orig_state_xml_copy );
    }catch( std::exception & )
    {
      assert( 0 );
    }
  };//undo
  
  
  WModelIndexSet selected = m_results->selectedIndexes();
  if( selected.empty() )
  {
    m_selected_row = -1;
    
    updateClearSelectionButton();
    
    m_assignPeakToSelected->hide();
    
    if( display )
      display->setIsotope( nullptr );
    
    auto redo = [](){
      InterSpec *viewer = InterSpec::instance();
      ReferencePhotopeakDisplay *display = viewer ? viewer->referenceLinesWidget() : nullptr;
      IsotopeSearchByEnergy *search = viewer ? viewer->nuclideSearch() : nullptr;
      assert( display && search );
      if( !display || !search )
        return;
      
      shared_ptr<void> ref_line_undo_sentry = display->getDisableUndoRedoSentry();
      shared_ptr<void> search_undo_sentry = search->getDisableUndoRedoSentry();
      display->setIsotope( nullptr );
      search->m_results->setSelectedIndexes( {} );
      search->m_selected_row = -1;
      search->resultSelectionChanged();
    };
    
    if( undoManager && (orig_selected_row >= 0) && !m_undo_redo_sentry.lock() )
      undoManager->addUndoRedoStep( undo, redo, "Clear nuclide search selection." );
    
    return;
  }//if( selected.empty() )
  
  const WModelIndex index = *selected.begin();
  m_selected_row = index.row();
  
  const SandiaDecay::Nuclide *nuc = m_model->nuclide( index );
  const SandiaDecay::Element *el = m_model->xrayElement( index );
  const ReactionGamma::Reaction *rctn = m_model->reaction( index );
  
  if( display )
  {
    if( nuc )
      display->setIsotope( nuc, m_model->assumedAge(index) );
    else if( el )
      display->setElement( el );
    else if( rctn )
      display->setReaction( rctn );
    
    string final_state_xml;
    display->serialize( final_state_xml );
    
    auto redo = [final_state_xml, index](){
      InterSpec *viewer = InterSpec::instance();
      ReferencePhotopeakDisplay *display = viewer ? viewer->referenceLinesWidget() : nullptr;
      IsotopeSearchByEnergy *search = viewer ? viewer->nuclideSearch() : nullptr;
      assert( display && search );
      if( !display || !search )
        return;
      
      shared_ptr<void> ref_line_undo_sentry = display->getDisableUndoRedoSentry();
      shared_ptr<void> search_undo_sentry = search->getDisableUndoRedoSentry();
      
      search->m_results->setSelectedIndexes( {index} );
      search->resultSelectionChanged();
      
      //try
      //{
      //  string final_state_xml_copy = final_state_xml;
      //  display->deSerialize( final_state_xml_copy );
      //}catch( std::exception & )
      //{
      //  assert( 0 );
      //}
    };//redo
    
    
    if( undoManager && !m_undo_redo_sentry.lock() )
      undoManager->addUndoRedoStep( undo, redo, "Change search row" );
  }//if( display )
  
  
  const bool showBtn = ((nuc || el || rctn) && (nPeaksSearched > 0));
  m_assignPeakToSelected->setHidden( !showBtn );
  if( showBtn )
  {
    const string symbol = (nuc ? nuc->symbol : string())
                          + (el ? (el->symbol + " x-ray") : string())
                          + (rctn ? rctn->name() : string());
    WString btntxt = WString::trn( "isbe-assign-to-btn", nPeaksSearched ).arg( symbol );
    
    m_assignPeakToSelected->actionButton()->setText( btntxt );
    
    const int nPeaksOnNuc = numCurrentNuclideLinesOnPeaks(false);
    const int nPeaksNoIdOnNuc = numCurrentNuclideLinesOnPeaks(true);
    
    const bool hide_menu_btn = (nPeaksOnNuc <= nPeaksSearched);
    m_assignPeakToSelected->dropDownButton()->setHidden( hide_menu_btn );
    
    WPopupMenu *menu = m_assignPeakToSelected->menu();
    assert( menu );
    if( menu )
    {
      // Clear the menu
      for( const auto item : menu->items() )
        menu->removeItem( item );
      
      if( !hide_menu_btn )
      {
        WString txt = WString::trn( "isbe-assign-searched-energies", nPeaksSearched ).arg( symbol );
        WMenuItem *item = menu->addItem( txt );
        item->triggered().connect( this, &IsotopeSearchByEnergy::assignSearchedOnPeaksToSelectedNuclide );
        
        if( nPeaksNoIdOnNuc > nPeaksSearched )
        {
          txt = WString::tr("isbe-assign-matching-peaks-no-id").arg( symbol );
          item = menu->addItem( txt );
          item->triggered().connect( boost::bind( &IsotopeSearchByEnergy::assignPeaksNearReferenceLinesToSelectedNuclide, this, true) );
        }//if( nPeaksNoIdOnNuc > nPeaksSearched )
        
        if( nPeaksOnNuc > nPeaksSearched )
        {
          txt = WString::tr("isbe-assign-all-matching-peaks").arg( symbol );
          item = menu->addItem( txt );
          item->triggered().connect( boost::bind( &IsotopeSearchByEnergy::assignPeaksNearReferenceLinesToSelectedNuclide, this, false) );
        }//if( nPeaksOnNuc > nPeaksSearched )
      }//if( !hide_menu_btn )
    }//if( menu )
    
    m_assignPeakToSelected->dropDownButton()->setHidden( hide_menu_btn );
  }//if( showBtn )
  
  updateClearSelectionButton();
}//void resultSelectionChanged()


int IsotopeSearchByEnergy::numSearchEnergiesOnPeaks()
{
  PeakModel *pmodel = m_viewer ? m_viewer->peakModel() : nullptr;
  assert( pmodel );
  
  const vector<PeakModel::PeakShrdPtr> peaks = peaks_searched( pmodel, searches() );
  
  return static_cast<int>( peaks.size() );
}//int numSearchEnergiesOnPeaks()


int IsotopeSearchByEnergy::numCurrentNuclideLinesOnPeaks( const bool require_peaks_with_no_id )
{
  PeakModel *pmodel = m_viewer ? m_viewer->peakModel() : nullptr;
  if( !pmodel )
    return 0;
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks = pmodel->peaks();
  if( !all_peaks )
    return 0;
  
  vector<shared_ptr<const PeakDef>> candidate_peaks;
  for( const shared_ptr<const PeakDef> &p : *all_peaks )
  {
    if( !require_peaks_with_no_id || !p->hasSourceGammaAssigned() )
      candidate_peaks.push_back( p );
  }//for( const shared_ptr<const PeakDef> &p : all_peaks )
  
  
  const WModelIndexSet selected = m_results->selectedIndexes();
  if( selected.empty() )
    return 0;
  
  const WModelIndex row_index = *selected.begin();
  
  const SandiaDecay::Nuclide *nuc = m_model->nuclide( row_index );
  const SandiaDecay::Element *el = m_model->xrayElement( row_index );
  const ReactionGamma::Reaction *rctn = m_model->reaction( row_index );
  
  if( !nuc && !el && !rctn )
    return 0;
  
  ReferencePhotopeakDisplay *refline_widget = m_viewer->referenceLinesWidget();
  if( !refline_widget )
    return 0;
  
  const ReferenceLineInfo &reflines = refline_widget->currentlyShowingNuclide();
  
  if( reflines.m_validity != ReferenceLineInfo::InputValidity::Valid )
    return 0;
  
  assert( (nuc && (reflines.m_nuclide == nuc))
         || (el && (reflines.m_element == el))
         || (rctn && reflines.m_reactions.count(rctn)) );
  
  if( !(nuc && (reflines.m_nuclide == nuc))
     && !(el && (reflines.m_element == el))
     && !(rctn && reflines.m_reactions.count(rctn)) )
  {
    return 0;
  }
  
  // Get the reference line energy and normalized intensities; we will sort them
  //  to speed (worst case) matching to peaks
  vector<pair<double,double>> ref_lines_energy_br;
  ref_lines_energy_br.reserve( reflines.m_ref_lines.size() );
  
  for( const ReferenceLineInfo::RefLine &line : reflines.m_ref_lines )
  {
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
      case ReferenceLineInfo::RefLine::Particle::Beta:
        continue;
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
        break;
    }//switch( line.m_particle_type )
    
    // I guess we *could* be interested in Rel. Eff. line of sum or escape peaks,
    //  but for the moment, we'll just use normal gammas and x-rays
    switch( line.m_source_type )
    {
      case ReferenceLineInfo::RefLine::RefGammaType::Normal:
      case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
        break;
        
      case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
      case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
      case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
      case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
        continue;
        break;
    }//switch( line.m_source_type )
    
    if( line.m_normalized_intensity > 0.0 ) // TODO: have some reasonable threshold, or logic to actually consider this gamma
      ref_lines_energy_br.emplace_back( line.m_energy, line.m_normalized_intensity );
  }//for( const ReferenceLineInfo::RefLine &line : refline.m_ref_lines )
  
  std::sort( begin(ref_lines_energy_br), end(ref_lines_energy_br) );
  
  
  set<PeakModel::PeakShrdPtr> peaks;
  for( const shared_ptr<const PeakDef> &peak : candidate_peaks )
  {
    const double mean = peak->mean();
    const double width = peak->gausPeak() ? peak->fwhm() : 0.5*peak->roiWidth();
    //What is used in `InterSpec::setIsotopeSearchEnergy(energy)` to match search energy to peak
    const double match_tol = (3.0/2.35482)*width;
    
    const auto lb = std::lower_bound(begin(ref_lines_energy_br), end(ref_lines_energy_br), make_pair(mean - match_tol, 0.0) );
    const auto ub = std::upper_bound(lb, end(ref_lines_energy_br), make_pair(mean + match_tol, 0.0) );
    
    for( auto iter = lb; iter != ub; ++iter )
    {
      if( fabs(iter->first - mean) <= match_tol )
      {
        peaks.insert( peak );
        break;
      }
    }//for( auto iter = lb; iter != ub; ++iter )
  }//for( const shared_ptr<const PeakDef> &peak : candidate_peaks )
  
  return static_cast<int>( peaks.size() );
}//int numCurrentNuclideLinesOnPeaks( const bool require_peaks_with_no_id )


void IsotopeSearchByEnergy::assignSearchedOnPeaksToSelectedNuclide()
{
  UndoRedoManager::PeakModelChange peak_undo_creator;
  
  PeakModel *pmodel = m_viewer ? m_viewer->peakModel() : nullptr;
  assert( pmodel );
  if( !pmodel )
    return;
  
  const vector<PeakModel::PeakShrdPtr> peaks = peaks_searched( pmodel, searches() );
  if( peaks.empty() )
    return;
  
  const WModelIndexSet selected = m_results->selectedIndexes();
  if( selected.empty() )
    return;
  
  const WModelIndex row_index = *selected.begin();
  
  const SandiaDecay::Nuclide *nuc = m_model->nuclide( row_index );
  const SandiaDecay::Element *el = m_model->xrayElement( row_index );
  const ReactionGamma::Reaction *rctn = m_model->reaction( row_index );
   
  assert( nuc || el || rctn );
  
  WString nucstr;
  if( nuc )
    nucstr = WString::fromUTF8( nuc->symbol );
  else if( el )
    nucstr = WString( "{1} {2}" ).arg(el->symbol).arg( WString::tr("x-ray") );
  else if( rctn )
    nucstr = WString::fromUTF8( rctn->name() );
  else
    return;
  
  for( const auto &peak : peaks )
  {
    // Dont set the nuclide, if it is already assigned to this nuc/el/rct
    //  (the user may have set a specific gamma energy, and we dont want
    //  change this by guessing which one it should be).
    if( (nuc && (nuc == peak->parentNuclide()))
       || (el && (el == peak->xrayElement()))
       || (rctn && (rctn == peak->reaction())) )
    {
      continue;
    }
    
    const WModelIndex peak_index = pmodel->indexOfPeak(peak);
    assert( peak_index.isValid() );
    const int peak_row = peak_index.row();
    
    const WModelIndex iso_index = pmodel->index( peak_row, PeakModel::Columns::kIsotope );
    pmodel->setData( iso_index, boost::any(nucstr) );
    
    // PeakModel doesnt know about reference line color, so lets get the color of the reference
    //  lines (which, if previous peaks have been assigned current nuc/el/rctn, then ref lines
    //  color should match them), and change the peaks color to that.
    ReferencePhotopeakDisplay *refLines = m_viewer->referenceLinesWidget();
    if( refLines )
    {
      // We'll double-check that the current reference lines are the selected nuc/el/rctn,
      //  but this should always be the case
      const ReferenceLineInfo &current = refLines->currentlyShowingNuclide();
      if( (nuc && (current.m_nuclide == nuc))
         || (el && (current.m_element == el))
         || (rctn && current.m_reactions.count(rctn)) )
      {
        const WColor &color = current.m_input.m_color;
        if( !color.isDefault() )
        {
          const WModelIndex color_index = pmodel->index( peak_row, PeakModel::Columns::kPeakLineColor );
          pmodel->setData( color_index, boost::any( WString(color.cssText()) ) );
        }
      }else
      {
        assert( 0 );
      }
    }
  }//for( const auto &peak : peaks )
  
  // Hide the button, so the user knows something has happened.
  if( m_assignPeakToSelected )
    m_assignPeakToSelected->hide();
}//void assignSearchedOnPeaksToSelectedNuclide()


void IsotopeSearchByEnergy::assignPeaksNearReferenceLinesToSelectedNuclide( const bool require_no_peak_id )
{
  PeakSearchGuiUtils::assign_peak_nuclides_from_reference_lines( m_viewer, require_no_peak_id, true );
}//void assignPeaksNearReferenceLinesToSelectedNuclide( const bool require_no_peak_id )


void IsotopeSearchByEnergy::clearSelectionAndRefLines()
{
  m_results->setSelectedIndexes( {} );
  ReferencePhotopeakDisplay *refLines = m_viewer->referenceLinesWidget();
  if( refLines )
    refLines->clearAllLines();
  
  m_clearRefLines->hide();
}//void clearSelectionAndRefLines()


void IsotopeSearchByEnergy::updateClearSelectionButton()
{
  if( !m_clearRefLines )
    return;
  
  WModelIndexSet selected = m_results->selectedIndexes();
  if( !selected.empty() )
  {
    m_clearRefLines->show();
    return;
  }
  
  ReferencePhotopeakDisplay *refLines = m_viewer->referenceLinesWidget();
  if( !refLines )
  {
    m_clearRefLines->hide();
    return;
  }
  
  m_clearRefLines->setHidden( refLines->showingNuclides().empty() );
}//void updateClearSelectionButton();



void IsotopeSearchByEnergy::handleRefLinesUpdated()
{
  updateClearSelectionButton();
  
  /*
  ReferencePhotopeakDisplay *refLines = m_viewer->referenceLinesWidget();
  if( !refLines )
  {
    updateClearSelectionButton();
    return;
  }
  
  const ReferenceLineInfo &current = refLines->currentlyShowingNuclide();
  
  if( current.m_validity != ReferenceLineInfo::InputValidity::Valid )
  {
    m_results->setSelectedIndexes( {} );
    updateClearSelectionButton();
    return;
  }
  
  const WModelIndexSet selected = m_results->selectedIndexes();
  
  if( !selected.empty() )
  {
    const WModelIndex row_index = *selected.begin();
    
    const SandiaDecay::Nuclide *nuc = m_model->nuclide( row_index );
    const SandiaDecay::Element *el = m_model->xrayElement( row_index );
    const ReactionGamma::Reaction *rctn = m_model->reaction( row_index );
    
    assert( nuc || el || rctn );
    
    if( (nuc && (current.m_nuclide == nuc))
         || (el && (current.m_element == el))
         || (rctn && current.m_reactions.count(rctn)) )
    {
      updateClearSelectionButton();
      return;
    }
  }//if( !selected.empty() )
  
  const int nrows = m_model->rowCount();
  for( int row = 0; row < nrows; ++row )
  {
    const auto row_index = m_model->index(row, 0);
    const SandiaDecay::Nuclide *nuc = m_model->nuclide( row_index );
    const SandiaDecay::Element *el = m_model->xrayElement( row_index );
    const ReactionGamma::Reaction *rctn = m_model->reaction( row_index );
    
    if( (nuc && (current.m_nuclide == nuc))
         || (el && (current.m_element == el))
         || (rctn && current.m_reactions.count(rctn)) )
    {
      m_results->setSelectedIndexes( {row_index} );
      updateClearSelectionButton();
      return;
    }
  }//for( int row = 0; row < nrows; ++row )
  
  m_results->setSelectedIndexes( {} );
  updateClearSelectionButton();
   */
}//void IsotopeSearchByEnergy::handleRefLinesUpdated()



IsotopeSearchByEnergy::WidgetState IsotopeSearchByEnergy::guiState() const
{
  WidgetState state;
  
  state.MinBranchRatio = m_minBranchRatio->valueText();
  state.MinHalfLife = m_minHalfLife->valueText();
  state.NextSearchEnergy = m_nextSearchEnergy;
  const int cat_index = m_search_category_select->currentIndex();
  state.CategoryIndex = static_cast<int>( std::max( 0, cat_index ) );
  
  const vector<WWidget *> children = m_searchEnergies->children();
  
  for( const WWidget *w : children )
  {
    const SearchEnergy *ww = dynamic_cast<const SearchEnergy *>( w );
    //if( ww && (ww->energy() > 0.000001) )
    state.SearchEnergies.push_back( {ww->energy(), ww->window()} );
  }//for( const WWebWidget *w : children )
  
  return state;
}//guiState()


void IsotopeSearchByEnergy::setGuiState( const WidgetState &state, const bool renderOnChart )
{
  vector<SearchEnergy *> origSearches;
  for( WWidget *w : m_searchEnergies->children() )
  {
    SearchEnergy *ww = dynamic_cast<SearchEnergy *>( w );
    if( ww )
      origSearches.push_back( ww );
  }//for( const WWebWidget *w : children )
  
  for( size_t i = 1; i < origSearches.size(); ++i )
    removeSearchEnergy( origSearches[i] );
  
  m_minBranchRatio->setValueText( state.MinBranchRatio );
  m_minHalfLife->setValueText( state.MinHalfLife );
  m_nextSearchEnergy = state.NextSearchEnergy;
  
  size_t cat_index = state.CategoryIndex;
  const size_t num_cats = m_search_categories.size();
  if( cat_index >= num_cats )
    cat_index = 0;
  m_search_category_select->setCurrentIndex( static_cast<int>(cat_index) );
  if( cat_index < num_cats )
    m_search_category_select->setToolTip( m_search_categories[cat_index].m_description );
  
  int nnode = 0;
  for( const pair<double,double> &ene : state.SearchEnergies )
  {
    SearchEnergy *searcher = (SearchEnergy *)0;
    if( nnode++ )
    {
      searcher = addNewSearchEnergy();
    }else
    {
      for( WWidget *w : m_searchEnergies->children() )
      {
        searcher = dynamic_cast<SearchEnergy *>( w );
        if( searcher )
          break;
      }//for( WWidget *w : m_searchEnergies->children() )
      
      assert( searcher );
      if( !searcher )
        searcher = addNewSearchEnergy();
    }//if( nnode++ ) / else
    
    searcher->setEnergy( ene.first );
    searcher->setWindow( ene.second );
  }//for( const pair<double,double> &ene : state.SearchEnergies )
  
  
  if( renderOnChart )
  {
    const vector<SearchEnergy *> vals = searches();
    if( vals.size() && (fabs(vals[0]->energy())>0.01 || vals.size()>2) )
      startSearch( true );
    loadSearchEnergiesToClient();
  }else
  {
    //This next call appears to be necassary or else search energies will show
    //  but I'm not certain why...
    clearSearchEnergiesOnClient();
  }//if( renderOnChart )
}//void setGuiState( const WidgetState &state );


void IsotopeSearchByEnergy::WidgetState::serialize( std::string &xml_data ) const
{
  rapidxml::xml_document<char> doc;
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *searchEnergiesNode, *searchnode;
  rapidxml::xml_attribute<> *attr;
  
  name = "IsotopeSearchByEnergy";
  base_node = doc.allocate_node( rapidxml::node_element, name );
  doc.append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  value = doc.allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc.allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  name = "MinBranchRatio";
  value = doc.allocate_string( MinBranchRatio.toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "MinHalfLife";
  value = doc.allocate_string( MinHalfLife.toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "NextSearchEnergy";
  value = doc.allocate_string( std::to_string(NextSearchEnergy).c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "CategoryIndex";
  value = doc.allocate_string( std::to_string(CategoryIndex).c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "SearchEnergies";
  searchEnergiesNode = doc.allocate_node( rapidxml::node_element, name );
  base_node->append_node( searchEnergiesNode );
  
  
  for( const std::pair<double,double> &w : SearchEnergies )
  {
    const double &energy = w.first;
    const double &window = w.second;
    
    if( energy > 0.000001 )
    {
      searchnode = doc.allocate_node( rapidxml::node_element, "SearchEnergy" );
      searchEnergiesNode->append_node( searchnode );
      
      value = doc.allocate_string( std::to_string(energy).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Energy", value );
      searchnode->append_node( node );
      
      value = doc.allocate_string( std::to_string(window).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Window", value );
      searchnode->append_node( node );
    }//if( ww && (ww->energy() > 0.000001) )
  }//for( const WWebWidget *w : children )
  
  xml_data.clear();
  rapidxml::print(std::back_inserter(xml_data), doc, 0);
}//void IsotopeSearchByEnergy::WidgetState::serialize( std::string &xml_data ) const


void IsotopeSearchByEnergy::WidgetState::deSerialize( std::string &xml_data,
                                                  const std::vector<NucSearchCategory> &categories )
{
  rapidxml::xml_document<char> doc;
  const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
  if( xml_data.size() )
    doc.parse<flags>( &(xml_data[0]) );
    
  rapidxml::xml_attribute<char> *attr;
  rapidxml::xml_node<char> *base_node, *node, *search_nodes, *value_node;
    
  base_node = doc.first_node( "IsotopeSearchByEnergy", 21 );
  if( !base_node )
    throw runtime_error( "Couldnt get base node, IsotopeSearchByEnergy" );
    
  int version = 0;
  attr = base_node->first_attribute( "version", 7 );
  if( !attr || !attr->value()
      || !(stringstream(attr->value()) >> version)
      || (version != sm_xmlSerializationVersion) )
    throw runtime_error( "Mising or invalid NuclideSearchByEnergy version" );

  node = base_node->first_node( "MinBranchRatio", 14 );
  if( !node || !node->value() )
    throw runtime_error( "Missing MinBranchRatio node" );
  MinBranchRatio = node->value();
    
  node = base_node->first_node( "MinHalfLife", 11 );
  if( !node || !node->value() )
    throw runtime_error( "Missing MinHalfLife node" );
  MinHalfLife = node->value();
    
  node = base_node->first_node( "NextSearchEnergy", 16 );
  if( !node || !node->value()
      || !(stringstream(node->value()) >> NextSearchEnergy) )
    throw runtime_error( "Missing/invalid NextSearchEnergy node" );

  node = base_node->first_node( "CategoryIndex", 13 );
  if( node )
  {
    CategoryIndex = 0;
    int cat_index =  0;
    const bool success = SpecUtils::parse_int(node->value(), node->value_size(), cat_index);
    assert( success );
    if( !success || (cat_index < 0) || (cat_index >= categories.size()) )
      cat_index = 0;
    CategoryIndex = static_cast<size_t>( cat_index );
  }else
  {
    // 20241114 (but not merged into main until after v1.0.13), we changed from some checkboxes
    //  for gamma/xray/reaction options, to a selection option to select things.
    //  For could just default to use the zeroth category.
    // However, we can basically remain backward compatible, if we do whats commented out here,
    //  and basically find the option in the drop-down that matches the options selected.
    
    CategoryIndex = 0;
    
     
    bool IncludeGammas = true, IncludeXRays = true, IncludeReactions = false;
    node = base_node->first_node( "IncludeGammas", 13 );
    if( node && node->value() && strlen(node->value()))
      IncludeGammas = (node->value()[0] != '0');

    node = base_node->first_node( "IncludeXRays", 12 );
    if( node && node->value() && strlen(node->value()))
      IncludeXRays = (node->value()[0] != '0');
      
    node = base_node->first_node( "IncludeReactions", 16 );
    if( node && node->value() && strlen(node->value()))
      IncludeReactions = (node->value()[0] != '0');
    
    for( size_t i = 0; i < categories.size(); ++i )
    {
      const NucSearchCategory &cat = categories[i];
      if( (cat.m_nuclides == IncludeGammas)
         && (cat.m_fluorescence_xrays == IncludeXRays)
         && (cat.m_reactions == IncludeReactions)
         && !cat.m_alphas
         && !cat.m_beta_endpoint
         && cat.m_specific_elements.empty()
         && cat.m_specific_nuclides.empty()
         && cat.m_specific_reactions.empty() )
      {
        CategoryIndex = i;
        break;
      }
    }//for( loop over m_search_categories )
  }//if( CategoryIndex node ) / else
  
    
  search_nodes = base_node->first_node( "SearchEnergies", 14 );
  if( !search_nodes )
    throw runtime_error( "Missing SearchEnergies node" );

  SearchEnergies.clear();
  for( node = search_nodes->first_node( "SearchEnergy", 12 );
        node; node = node->next_sibling( "SearchEnergy", 12 ) )
  {
    double energy, window;

    value_node = node->first_node( "Energy", 6 );
    if( !value_node || !value_node->value()
        || !(stringstream(value_node->value()) >> energy) )
      throw runtime_error( "Missing/invalid SearchEnergy energy entry" );
      
    value_node = node->first_node( "Window", 6 );
    if( !value_node || !value_node->value()
        || !(stringstream(value_node->value()) >> window) )
      throw runtime_error( "Missing/invalid SearchEnergy window entry" );
      
    SearchEnergies.push_back( {energy, window} );
  }//for( loop over nodes )
}//void IsotopeSearchByEnergy::WidgetState::deSerialize( std::string &xml_data )


bool IsotopeSearchByEnergy::WidgetState::operator==(const WidgetState &rhs) const
{
  return ((MinBranchRatio == rhs.MinBranchRatio)
          && (MinHalfLife == rhs.MinHalfLife)
          && (NextSearchEnergy == rhs.NextSearchEnergy)
          && (CategoryIndex == rhs.CategoryIndex)
          && (SearchEnergies == rhs.SearchEnergies)
  );
}//bool WidgetState::operator==(const WidgetState &rhs) const;


void IsotopeSearchByEnergy::NucSearchCategory::deSerialize( const rapidxml::xml_node<char> * const xml_data )
{
  if( !xml_data || !xml_data->name_size() )
    throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory::deSerialize: invalid input." );
    
  if( !XML_NAME_ICOMPARE(xml_data, "NucSearchCategory") )
    throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory:"
                        " input must be a NucSearchCategory elelement." );
  
  auto get_bool_val = [xml_data]( const rapidxml::xml_base<char> *el ) -> bool {
    return (el && el->value_size() && (el->value()[0] == '1'));
  };
  
  const rapidxml::xml_node<char> *name_el = XML_FIRST_NODE(xml_data, "Name");
  const string name = SpecUtils::xml_value_str( name_el );
  if( name.empty() )
    throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory: No valid name category." );
  
  const rapidxml::xml_attribute<char> *trans_ettrib = XML_FIRST_ATTRIB(name_el, "internationalize");
  bool translate = get_bool_val( trans_ettrib );
  
  m_name = translate ? Wt::WString::tr(name) : Wt::WString::fromUTF8(name);
  
  const rapidxml::xml_node<char> *desc_el = XML_FIRST_NODE(xml_data, "Description");
  const string desc = SpecUtils::xml_value_str( desc_el );
  const rapidxml::xml_node<char> *trans_el = XML_FIRST_NODE(desc_el, "internationalize");
  translate = get_bool_val( trans_el );
  m_description = translate ? Wt::WString::tr(desc) : Wt::WString::fromUTF8(desc);
  
  m_minBr = 0.0;
  const rapidxml::xml_node<char> *min_br_el = XML_FIRST_NODE(xml_data, "MinBr");
  if( min_br_el )
  {
    if( !SpecUtils::parse_double(min_br_el->value(), min_br_el->value_size(), m_minBr) )
      throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory: min BR invalid" );
    if( (m_minBr < 0.0) || (m_minBr > 1.0) )
      throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory: min BR must be between 0 and 1" );
  }//if( min_br_el )
  
  m_minHl = 0.0;
  const rapidxml::xml_node<char> *min_hl_el = XML_FIRST_NODE(xml_data, "MinHl");
  if( min_hl_el )
  {
    const string val = SpecUtils::xml_value_str(min_hl_el);
    m_minHl = PhysicalUnits::stringToTimeDuration(val);
    if( m_minHl < 0.0 )
      m_minHl = 0.0;
  }//if( min_hl_el )
  
  m_nuclides = get_bool_val( XML_FIRST_NODE(xml_data, "AllowNucGammas") );
  m_fluorescence_xrays = get_bool_val( XML_FIRST_NODE(xml_data, "AllowFluorXrays") );
  m_reactions = get_bool_val( XML_FIRST_NODE(xml_data, "AllowReactions") );
  m_alphas = get_bool_val( XML_FIRST_NODE(xml_data, "AllowAlphas") );
  m_beta_endpoint = get_bool_val( XML_FIRST_NODE(xml_data, "AllowBetas") );
  m_no_progeny = get_bool_val( XML_FIRST_NODE(xml_data, "NoProgeny") );
  
  if( (m_alphas || m_beta_endpoint)
      && ((m_alphas == m_beta_endpoint) || m_reactions || m_fluorescence_xrays || m_nuclides) )
  {
    throw runtime_error( "IsotopeSearchByEnergy::NucSearchCategory: If alpha or betas are"
                        " specified, then no other source may be specified." );
  }
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Couldnt get DecayDataBaseServer." );
  
  m_specific_elements.clear();
  m_specific_nuclides.clear();
  m_specific_reactions.clear();
  const rapidxml::xml_node<char> *elements_el = XML_FIRST_NODE(xml_data, "Elements");
  const rapidxml::xml_node<char> *nuclides_el = XML_FIRST_NODE(xml_data, "Nuclides");
  const rapidxml::xml_node<char> *reactions_el = XML_FIRST_NODE(xml_data, "Reactions");
  if( elements_el )
  {
    XML_FOREACH_CHILD( element, elements_el, "Element" )
    {
      const string el_str = SpecUtils::xml_value_str(element);
      const SandiaDecay::Element * const el = db->element(el_str);
      if( el )
        m_specific_elements.push_back( el );
      else
        cerr << "Couldnt find element '" << el_str << "'." << endl;
    }//for( loop over <Elements> )
  }//if( elements_el )
  
  if( nuclides_el )
  {
    XML_FOREACH_CHILD( nuclide, nuclides_el, "Nuclide" )
    {
      const string nuc_str = SpecUtils::xml_value_str(nuclide);
      const SandiaDecay::Nuclide * const nuc = db->nuclide(nuc_str);
      if( nuc )
        m_specific_nuclides.push_back( nuc );
      else
        cerr << "Couldnt find nuclide '" << nuc_str << "'." << endl;
    }//for( loop over <Nuclide> )
  }//if( nuclides_el )
  
  if( reactions_el )
  {
    const ReactionGamma *reactionDb = nullptr;
   
    set<const ReactionGamma::Reaction *> reactions_seen;
      
    XML_FOREACH_CHILD( reaction, reactions_el, "Reaction" )
    {
      if( !reactionDb )
      {
        try
        {
          reactionDb = ReactionGammaServer::database();
        }catch(...)
        {
          cerr << "NucSearchCategory::deSerialize(): Failed to open gamma reactions XML file" << endl;
        }

        if( !reactionDb )
          break;
      }//if( !reactionDb )

      const string rectn_str = SpecUtils::xml_value_str(reaction);
        
      try
      {
        vector<ReactionGamma::ReactionPhotopeak> photopeaks;
        reactionDb->gammas( rectn_str, photopeaks);
          
        for( const ReactionGamma::ReactionPhotopeak &gamma : photopeaks )
        {
          if( gamma.reaction && !reactions_seen.count(gamma.reaction) )
          {
            reactions_seen.insert( gamma.reaction );
            m_specific_reactions.push_back( gamma.reaction );
              // We're not exiting the loop here, because some of the element reactions
              // may be composed of specific nuclide reactions
          }//if( gamma.reaction && !reactions_seen.count(gamma.reaction) )
        }//for( const ReactionGamma::ReactionPhotopeak &gamma : photopeaks )
      }catch( std::exception &e )
      {
        cerr << "NucSearchCategory::deSerialize(): Failed to get reaction '" << rectn_str << "'" << endl;
      }
    }//for( loop over <Reactions> )
  }//if( reactions_el )
  
  if( !m_specific_nuclides.empty() && !m_nuclides && !m_alphas && !m_beta_endpoint )
    throw runtime_error( "Specific nuclides specified, but AllowNucGammas, AllowAlphas, and AllowBetas is false" );
  
  if( !m_specific_elements.empty() && !m_fluorescence_xrays )
    throw runtime_error( "Specific elements specified, but AllowFluorXrays is false" );
  
  if( !m_specific_reactions.empty() && !m_reactions )
    throw runtime_error( "Specific reactions specified, but AllowReactions is false" );
}//void IsotopeSearchByEnergy::NucSearchCategory::deSerialize(...)


void IsotopeSearchByEnergy::serialize( std::string &xml_data ) const
{
  const WidgetState state = guiState();
  state.serialize( xml_data );
}//void serialize( std::string &xmlOutput ) const


void IsotopeSearchByEnergy::deSerialize( std::string &xml_data,
                                         const bool renderOnChart )
{
  std::shared_ptr<void> search_undo_sentry = getDisableUndoRedoSentry();
  
  try
  {
    WidgetState state;
    state.deSerialize( xml_data, m_search_categories );
    setGuiState( state, renderOnChart );
  }catch( std::exception &e )
  {
    cerr << "IsotopeSearchByEnergy::deSerialize(...) caught: " << e.what() << endl;
    stringstream msg;
    msg << "Error opening displayed photopeaks from database for search: " << e.what();
    passMessage( msg.str(), WarningWidget::WarningMsgHigh );
  }// try / catch
}//void IsotopeSearchByEnergy::deSerialize(...)


void IsotopeSearchByEnergy::addUndoRedoPoint()
{
  const WidgetState prev_state = m_state;
  m_state = guiState();
 
  // Dont update state if we are in an undo/redo
  UndoRedoManager *undoManager = m_viewer->undoRedoManager();
  if( undoManager && !undoManager->canAddUndoRedoNow() )
    return;
  
  if( !m_undo_redo_sentry.lock() )
  {
    if( undoManager && !(m_state == prev_state) )
    {
      auto undo = [prev_state](){
        InterSpec *viewer = InterSpec::instance();
        ReferencePhotopeakDisplay *display = viewer ? viewer->referenceLinesWidget() : nullptr;
        IsotopeSearchByEnergy *search = viewer ? viewer->nuclideSearch() : nullptr;
        assert( display && search );
        if( !display || !search )
          return;
        
        shared_ptr<void> ref_line_undo_sentry = display->getDisableUndoRedoSentry();
        shared_ptr<void> search_undo_sentry = search->getDisableUndoRedoSentry();
        
        try
        {
          search->setGuiState( prev_state, true );
        }catch(std::exception &)
        {
          assert(0);
        }
      };//undo
      
      const WidgetState current_state = m_state;
      auto redo = [current_state](){
        InterSpec *viewer = InterSpec::instance();
        ReferencePhotopeakDisplay *display = viewer ? viewer->referenceLinesWidget() : nullptr;
        IsotopeSearchByEnergy *search = viewer ? viewer->nuclideSearch() : nullptr;
        assert( display && search );
        if( !display || !search )
          return;
        
        shared_ptr<void> ref_line_undo_sentry = display->getDisableUndoRedoSentry();
        shared_ptr<void> search_undo_sentry = search->getDisableUndoRedoSentry();
        
        try
        {
          search->setGuiState( current_state, true );
        }catch( std::exception &)
        {
          assert(0);
        }
      };//redo
      
      undoManager->addUndoRedoStep( undo, redo, "Update nuclide search." );
    }//if( undoManager )
  }//if( !m_undo_redo_sentry.lock() )
}//void addUndoRedoPoint()


void IsotopeSearchByEnergy::startSearch( const bool refreshBr )
{
  addUndoRedoPoint();
  
  if( refreshBr )
    EnergyToNuclideServer::setLowerLimits( m_minHl, m_minBr );
  
  loadSearchEnergiesToClient();
  
  vector<double> energies, windows;
  const vector<WWidget *> children = m_searchEnergies->children();
  
  for( const WWidget *w : children )
  {
    const SearchEnergy *ww = dynamic_cast<const SearchEnergy *>( w );
    if( ww && (ww->energy() > 0.000001) )
    {
      energies.push_back( ww->energy() );
      windows.push_back( ww->window() );
    }
  }//for( const WWebWidget *w : children )
  
  WFlags<IsotopeSearchByEnergyModel::RadSource> srcs;
  
  std::vector<const SandiaDecay::Element *> elements;
  std::vector<const SandiaDecay::Nuclide *> nuclides;
  std::vector<const ReactionGamma::Reaction *> reactions;
  
  const int cat_index = std::max(0, m_search_category_select->currentIndex() );
  if( cat_index <= static_cast<int>(m_search_categories.size()) )
  {
    const NucSearchCategory &cat = m_search_categories[cat_index];
    
    if( cat.m_nuclides )
      srcs |= IsotopeSearchByEnergyModel::RadSource::NuclideGammaOrXray;
    if( cat.m_fluorescence_xrays )
      srcs |= IsotopeSearchByEnergyModel::RadSource::kFluorescenceXRay;
    if( cat.m_reactions )
      srcs |= IsotopeSearchByEnergyModel::RadSource::kReaction;
    if( cat.m_alphas )
      srcs |= IsotopeSearchByEnergyModel::RadSource::kAlpha;
    if( cat.m_beta_endpoint )
      srcs |= IsotopeSearchByEnergyModel::RadSource::kBetaEndpoint;
    if( cat.m_no_progeny )
      srcs |= IsotopeSearchByEnergyModel::RadSource::kNoNuclideProgeny;
    
    elements = cat.m_specific_elements;
    nuclides = cat.m_specific_nuclides;
    reactions = cat.m_specific_reactions;
  }//if( cat_index <= static_cast<int>(m_search_categories.size()) )
  
  WApplication *app = wApp;
  auto workingspace = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>();
  workingspace->energies = energies;
  workingspace->windows = windows;
  workingspace->sortColumn = m_model->sortColumn();
  workingspace->sortOrder = m_model->sortOrder();
  workingspace->isHPGe = PeakFitUtils::is_likely_high_res( m_viewer );
  workingspace->undoSentry = getDisableUndoRedoSentry(); //m_undo_redo_sentry.lock();
  
  std::shared_ptr<SpecMeas> foreground = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( foreground )
  {
    const set<int> &samplenums = m_viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
    auto userpeaks = foreground->peaks( samplenums );
    auto autopeaks = foreground->automatedSearchPeaks( samplenums );
    
    workingspace->foreground = foreground;
    workingspace->foreground_samplenums = samplenums;
    
    if( userpeaks )
      workingspace->user_peaks.insert( end(workingspace->user_peaks),
                                       begin(*userpeaks), end(*userpeaks) );
    if( autopeaks )
      workingspace->automated_search_peaks.insert(
                        end(workingspace->automated_search_peaks),
                        begin(*autopeaks), end(*autopeaks) );
    workingspace->detector_response_function = foreground->detector();
    workingspace->displayed_measurement = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  }//if( foreground )
  
  ++m_currentSearch;
  workingspace->searchdoneCallback = app->bind(
                          boost::bind( &IsotopeSearchByEnergy::hideSearchingTxt,
                                       this, m_currentSearch ) );
  
  //Verified below is safe if the WApplication instance is terminated before
  //  search results are completed, as well as if the WApplication isnt
  //  terminated but m_model is deleted.
  boost::function< void(void) > updatefcnt = app->bind( boost::bind(
                            &IsotopeSearchByEnergyModel::updateSearchResults,
                            m_model, workingspace ) );
  boost::function< void(void) > worker = boost::bind(
                                &IsotopeSearchByEnergyModel::setSearchEnergies,
                                workingspace, m_minBr, m_minHl, srcs,
                                std::move(elements), std::move(nuclides), std::move(reactions),
                                app->sessionId(), updatefcnt );
  WServer::instance()->ioService().boost::asio::io_service::post( worker );
  
  m_searching->show();
}//void startSearch()


void IsotopeSearchByEnergy::hideSearchingTxt( const int searchNum )
{
  if( searchNum == m_currentSearch )
    m_searching->hide();
}

IsotopeSearchByEnergy::~IsotopeSearchByEnergy()
{
  
}//IsotopeSearchByEnergy destructor



