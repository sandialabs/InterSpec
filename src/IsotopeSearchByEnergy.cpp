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
#include <Wt/WIOService>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>


#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"

#if( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#else
#include "InterSpec/SpectrumDisplayDiv.h"
#endif

using namespace Wt;
using namespace std;


const int IsotopeSearchByEnergy::sm_xmlSerializationVersion = 0;

namespace
{
  const WString ActiveSearchEnergyClass = "ActiveSearchEnergy";
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
#if ( USE_SPECTRUM_CHART_D3 )
                                              D3SpectrumDisplayDiv *chart,
#else
                                              SpectrumDisplayDiv *chart,
#endif
                                              Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_viewer( viewer ),
  m_chart( chart ),
  m_searchEnergies( NULL ),
  m_currentSearch( 0 ),
  m_searching( NULL ),
  m_results( NULL ),
  m_minBranchRatio( NULL ),
  m_minHalfLife( NULL ),
  m_model( NULL ),
  m_gammas( NULL ),
  m_xrays( NULL ),
  m_reactions( NULL ),
  m_nextSearchEnergy( 0 ),
  m_minBr( 0.0 ), m_minHl( 6000.0 * PhysicalUnits::second )
{
  wApp->useStyleSheet( "InterSpec_resources/IsotopeSearchByEnergy.css" );
  
  addStyleClass( "IsotopeSearchByEnergy" );
  
  WContainerWidget *searchConditions = new WContainerWidget( this );
  searchConditions->setStyleClass( "IsotopeSearchConditions" );
  
  m_searchEnergies = new WContainerWidget( searchConditions );
  m_searchEnergies->setStyleClass( "IsotopeSearchEnergies" );
  
  
  WContainerWidget *sourceTypes = new WContainerWidget( searchConditions );
  sourceTypes->setStyleClass( "IsotopeSourceTypes" );
  m_gammas = new WCheckBox( "Gammas", sourceTypes );
  m_gammas->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_gammas->setChecked();

  m_xrays = new WCheckBox( "X-rays", sourceTypes );
  m_xrays->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_xrays->setChecked();
  
  m_reactions = new WCheckBox( "Reactions", sourceTypes );
  m_reactions->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
//  m_reactions->setChecked();

  
  WContainerWidget *searchOptions = new WContainerWidget( searchConditions );
  searchOptions->setStyleClass( "IsotopeSearchMinimums" );

  
  auto helpBtn = new WContainerWidget( searchOptions );
  helpBtn->addStyleClass( "Wt-icon ContentHelpBtn" );
  helpBtn->clicked().connect( boost::bind( &HelpSystem::createHelpWindow, "nuclide-search-dialog" ) );
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  WContainerWidget *optionDiv = new WContainerWidget( searchOptions );
  WLabel *label = new WLabel( "Min. BR", optionDiv );
//  HelpSystem::attachToolTipOn( label,"Toggle or type minimum branching ratio.", showToolTips , HelpSystem::ToolTipPosition::Top);

  m_minBranchRatio = new NativeFloatSpinBox( optionDiv );
  string tip = "Minimum branching ratio.";
  HelpSystem::attachToolTipOn( m_minBranchRatio, tip, showToolTips , HelpSystem::ToolTipPosition::Top);
  
  m_minBranchRatio->setValue( m_minBr );
  m_minBranchRatio->setRange( 0.0f, 1.0f );
  m_minBranchRatio->setSingleStep( 0.1f );
  label->setBuddy( m_minBranchRatio );
  
  optionDiv = new WContainerWidget( searchOptions );
  label = new WLabel( "Min. HL", optionDiv );
 
  tip = "Minimum half life of nuclides to be searched.<br />"
    "<div>Age can be specified using a combination of time units, "
    "similar to '<b>5.3y 8d 22m</b>'.</div>"
    "<div>"
    "Acceptible time units: <b>year</b>, <b>yr</b>, <b>y</b>, <b>day</b>, <b>d</b>, <b>hrs</b>, <b>hour</b>, <b>h</b>, <b>minute</b>, "
    "<b>min</b>, <b>m</b>, <b>second</b>, <b>s</b>, <b>ms</b>, <b>microseconds</b>, <b>us</b>, <b>nanoseconds</b>, <b>ns</b>, or "
    "you can specify time period by <b>hh:mm:ss</b>. "
    "</div>"
    "<div>"
    "When multiple time periods are "
    "specified, they are summed, e.x. '1y6months 3m' is interpreted as "
    "18 months and 3 minutes"
    "</div>";
//    HelpSystem::attachToolTipOn( label, tip, showToolTips , HelpSystem::ToolTipPosition::Top);

  m_minHalfLife = new WLineEdit( "6000 s", optionDiv );
  HelpSystem::attachToolTipOn( m_minHalfLife, tip, showToolTips , HelpSystem::ToolTipPosition::Top );

  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationRegex, this );
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
  m_results->addStyleClass( "IsotopeSearchResultTable" );
  m_results->setAlternatingRowColors( true );
  //m_results->sortByColumn( IsotopeSearchByEnergyModel::Distance, Wt::AscendingOrder );
  m_results->sortByColumn( IsotopeSearchByEnergyModel::Column::ProfileDistance, Wt::DescendingOrder );
  
  
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
  
  m_results->setSelectionMode( Wt::SingleSelection );
  m_results->setSelectionBehavior( Wt::SelectRows );
  m_results->selectionChanged().connect( this, &IsotopeSearchByEnergy::resultSelectionChanged );
  
  
  m_searching = new WText( "Searching", this );
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


  minBrOrHlChanged();
}//IsotopeSearchByEnergy constuctor


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


void IsotopeSearchByEnergy::loadSearchEnergiesToClient()
{
#if( USE_SPECTRUM_CHART_D3 )
  
  vector<pair<double,double>> searchRegions;
  
  for( auto sw : searches() )
  {
    if( sw->energy() > 0.1 )
      searchRegions.push_back( make_pair(sw->energy(), sw->window()) );
  }
  
  m_chart->setSearchEnergies( searchRegions );
#else
  CanvasForDragging *can = m_chart->overlayCanvas();
  if( !can )
    return;
  
  string js;
  const vector<SearchEnergy *> searchW = searches();
  
  if( !searchW.empty() )
  {
    js += "[";
    char buffer[120];
    for( size_t index = 0; index < searchW.size(); ++index )
    {
      if( searchW[index]->energy() > 0.1 )
      {
        if( js.size() > 2 )
          js += ",";
        snprintf( buffer, sizeof(buffer), "[%.2f,%.2f]",
                 searchW[index]->energy(), searchW[index]->window() );
        js += buffer;
      }//if( searchW[index]->energy() > 0.1 )
    }//for( size_t index = 0; index < searchW.size(); ++index )
    js += "]";
  }//if( searchW.empty() ) / else
  
  if( js.size() < 2 )
    js = "null";
  doJavaScript( "$('#c"+can->id()+"').data('SearchEnergies'," + js + ");"
                "Wt.WT.DrawGammaLines('c" + can->id() + "',true);" );
#endif
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
#if( USE_SPECTRUM_CHART_D3 )
  m_chart->setSearchEnergies( vector<pair<double,double>>() );
#else
  //Calling wApp->doJavaScript(...) rather than this->doJavaScript(...) to
  //  ensure the command is done, even if this widget is deleted
  CanvasForDragging *can = m_chart->overlayCanvas();
  if( can )
    wApp->doJavaScript( "$('#c"+can->id()+"').data('SearchEnergies',null);"
                       "Wt.WT.DrawGammaLines('c" + can->id() + "',true);" );
#endif
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



void IsotopeSearchByEnergy::minBrOrHlChanged()
{
  if( m_minBranchRatio->validate() == WValidator::Valid )
    m_minBr = m_minBranchRatio->value();
  else
    m_minBranchRatio->setValue( m_minBr );
  
  try
  {
    const string hltxt = m_minHalfLife->valueText().narrow();
    m_minHl = PhysicalUnits::stringToTimeDuration( hltxt );
  }catch(...)
  {
    m_minHalfLife->setText( PhysicalUnits::printToBestTimeUnits(m_minHl,2) );
  }//try / catch

  startSearch( true );
}//void IsotopeSearchByEnergy::minBrOrHlChanged()


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
    display = m_viewer->referenceLinesWidget();
  }
  
  if( !display )
    return;
  
  
  WModelIndexSet selected = m_results->selectedIndexes();
  if( selected.empty() )
  {
    display->setIsotope( NULL );
    return;
  }//if( selected.empty() )
  

  if( m_model->nuclide( *selected.begin() ) )
  {
    const SandiaDecay::Nuclide *nuc = m_model->nuclide( *selected.begin() );
    const double age = m_model->assumedAge( *selected.begin() );
    
    display->setIsotope( nuc, age );
  }else if( m_model->xrayElement( *selected.begin() ) )
  {
    const SandiaDecay::Element *el = m_model->xrayElement( *selected.begin() );
    
    display->setElement( el );
  }else if( m_model->reaction( *selected.begin() ) )
  {
    const ReactionGamma::Reaction *rctn = m_model->reaction( *selected.begin() );
    
    display->setReaction( rctn );
  }//if( m_model->nuclide( *selected.begin() ) ) / else ...
}//void resultSelectionChanged()


void IsotopeSearchByEnergy::deSerialize( std::string &xml_data,
                                         const bool renderOnChart )
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
  
  try
  {
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace
    | rapidxml::parse_trim_whitespace;
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
    m_minBranchRatio->setValueText( node->value() );
    
    node = base_node->first_node( "MinHalfLife", 11 );
    if( !node || !node->value() )
      throw runtime_error( "Missing MinHalfLife node" );
    m_minHalfLife->setValueText( node->value() );
    
    node = base_node->first_node( "NextSearchEnergy", 16 );
    if( !node || !node->value()
        || !(stringstream(node->value()) >> m_nextSearchEnergy) )
      throw runtime_error( "Missing/invalid NextSearchEnergy node" );

    node = base_node->first_node( "IncludeGammas", 13 );
    if( node && node->value() && strlen(node->value()))
      m_gammas->setChecked( (node->value()[0] != '0') );
    else
      throw runtime_error( "Missing/invalid IncludeGammas node" );

    node = base_node->first_node( "IncludeXRays", 12 );
    if( node && node->value() && strlen(node->value()))
      m_xrays->setChecked( (node->value()[0] != '0') );
    else
      throw runtime_error( "Missing/invalid IncludeXRays node" );
    
    node = base_node->first_node( "IncludeReactions", 16 );
    if( node && node->value() && strlen(node->value()))
      m_reactions->setChecked( (node->value()[0] != '0') );
    else
      throw runtime_error( "Missing/invalid IncludeXRays node" );
    
    search_nodes = base_node->first_node( "SearchEnergies", 14 );
    if( !search_nodes )
      throw runtime_error( "Missing SearchEnergies node" );

    int nnode = 0;
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
        
        if( !searcher )
          throw runtime_error( "Serious logic error in "
                               "IsotopeSearchByEnergy::deSerialize(...)" );
      }//if( nnode++ ) / else
      
      searcher->setEnergy( energy );
      searcher->setWindow( window );
    }//for( loop over nodes )
    
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
  }catch( std::exception &e )
  {
    cerr << "IsotopeSearchByEnergy::deSerialize(...) caught: " << e.what() << endl;
    stringstream msg;
    msg << "Error opening displayed photopeaks from database for search: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void IsotopeSearchByEnergy::deSerialize( std::string &xml_data )


void IsotopeSearchByEnergy::serialize( std::string &xml_data  ) const
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
  value = doc.allocate_string( m_minBranchRatio->valueText().toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "MinHalfLife";
  value = doc.allocate_string( m_minHalfLife->valueText().toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "NextSearchEnergy";
  value = doc.allocate_string( std::to_string(m_nextSearchEnergy).c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeGammas";
  value = (m_gammas->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeXRays";
  value = (m_xrays->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeReactions";
  value = (m_reactions->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "SearchEnergies";
  searchEnergiesNode = doc.allocate_node( rapidxml::node_element, name );
  base_node->append_node( searchEnergiesNode );
  
  const vector<WWidget *> children = m_searchEnergies->children();
  
  for( const WWidget *w : children )
  {
    const SearchEnergy *ww = dynamic_cast<const SearchEnergy *>( w );
    if( ww && (ww->energy() > 0.000001) )
    {
      searchnode = doc.allocate_node( rapidxml::node_element, "SearchEnergy" );
      searchEnergiesNode->append_node( searchnode );
      
      value = doc.allocate_string( std::to_string(ww->energy()).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Energy", value );
      searchnode->append_node( node );
      
      value = doc.allocate_string( std::to_string(ww->window()).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Window", value );
      searchnode->append_node( node );
    }//if( ww && (ww->energy() > 0.000001) )
  }//for( const WWebWidget *w : children )
  
  xml_data.clear();
  rapidxml::print(std::back_inserter(xml_data), doc, 0);
}//void serialize( std::string &xmlOutput ) const


void IsotopeSearchByEnergy::startSearch( const bool refreshBr )
{
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
  
  if( m_gammas->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kGamma;
  if( m_xrays->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kXRay;
  if( m_reactions->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kReaction;
  
  
  WApplication *app = wApp;
  auto workingspace = make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>();
  workingspace->energies = energies;
  workingspace->windows = windows;
  workingspace->sortColumn = m_model->sortColumn();
  workingspace->sortOrder = m_model->sortOrder();
  
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
                                app->sessionId(), updatefcnt );
  WServer::instance()->ioService().post( worker );
  
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



