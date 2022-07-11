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

#include <tuple>
#include <limits>
#include <memory>
#include <vector>
#include <utility>

#include <boost/optional.hpp>

#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WPushButton>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/InterSpec.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/D3TimeChart.h"

using namespace Wt;
using namespace std;

using SpecUtils::SpecFile;
using SpecUtils::Measurement;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

namespace
{
  const std::string &jsbool( bool val )
  {
    static const std::string t = "true";
    static const std::string f = "false";
    
    return val ? t : f;
  };
}//namespace


#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WMenuItem>
#include <Wt/WTabWidget>
#include <Wt/WSelectionBox>
#include <Wt/WStackedWidget>
#include <Wt/WStackedWidget>

#include "InterSpec/NativeFloatSpinBox.h"

class D3TimeChartFilters : public WContainerWidget
{
  D3TimeChart *m_parentChart;
  
  enum InteracModeIndex
  {
    IM_Normal, IM_Zoom, IM_Pan, IM_Select, IM_Add, IM_Remove, IM_NumMode
  };
  
  WMenu *m_interactModeMenu;
  
  WContainerWidget *m_specTypeDiv;
  WMenu *m_specTypeSelect;
  
  NativeFloatSpinBox *m_lowerEnergy;
  NativeFloatSpinBox *m_upperEnergy;
  WPushButton *m_clearEnergyFilterBtn;
  WCheckBox *m_dontRebin;
  WCheckBox *m_hideNeutrons;
  NativeFloatSpinBox *m_gammaNeutRelEmphasis;
  
public:
  D3TimeChartFilters( D3TimeChart *parent )
    : WContainerWidget( parent ),
      m_parentChart( parent ),
      m_interactModeMenu( nullptr ),
      m_specTypeDiv( nullptr ),
      m_specTypeSelect( nullptr ),
      m_lowerEnergy( nullptr ),
      m_upperEnergy( nullptr ),
      m_clearEnergyFilterBtn( nullptr ),
      m_dontRebin( nullptr ),
      m_hideNeutrons( nullptr ),
      m_gammaNeutRelEmphasis( nullptr )
  {
    assert( parent );
    
    addStyleClass( "D3TimeChartFilters" );
    
    
    WContainerWidget *closeIcon = new WContainerWidget( this );
    closeIcon->addStyleClass( "closeicon-wtdefault" );
    closeIcon->clicked().connect( boost::bind( &D3TimeChart::showFilters, parent, false ) );
    
    
    WTabWidget *tabs = new WTabWidget( this );
    tabs->addStyleClass( "D3TimeFiltersTab" );
    
    WContainerWidget *interact = new WContainerWidget();
    interact->addStyleClass( "D3TimeInteractTab" );
    tabs->addTab(interact, "Interact"); //Returns a
    
    WStackedWidget *instructionsStack = new WStackedWidget();
    instructionsStack->addStyleClass( "D3TimeInteractInst" );
    
    WText *title = new WText( "Interaction mode:", interact );
    title->addStyleClass( "D3TimeInteractModeText" );
    
    m_interactModeMenu = new WMenu( instructionsStack, interact );
    m_interactModeMenu->addStyleClass( "D3TimeInteractMenu D3TimeInteractModeMenu LightNavMenu" );
    
    WMenuItem *item = nullptr;
    WText *instructions = nullptr;

    for( InteracModeIndex index = InteracModeIndex(0); index < IM_NumMode;
        index = InteracModeIndex(index + 1) )
    {
      item = nullptr;
      
      switch( index )
      {
        case InteracModeIndex::IM_Normal:
          instructions = new WText( "Left mouse: select<br />"
                                    "&nbsp;&nbsp;&nbsp;&nbsp;no modifiers: foreground<br />"
                                    "&nbsp;&nbsp;&nbsp;&nbsp;Option-key: background<br />"
                                    "&nbsp;&nbsp;&nbsp;&nbsp;shift-key: add times<br />"
                                    "&nbsp;&nbsp;&nbsp;&nbsp;Ctrl-key: remove times<br />"
                                    "Right mouse: zoom in/out<br />"
                                    "Wheel: zoom and left/right<br />" );
          item = m_interactModeMenu->addItem( "Normal", instructions );
          break;
          
        case InteracModeIndex::IM_Zoom:
          instructions = new WText( "Click drag right zoom in<br />"
                                    "Click drag left zoom out<br />"
                                    "Drag x-axis to pan<br />" );
          item = m_interactModeMenu->addItem( "Zoom", instructions );
          break;
          
        case InteracModeIndex::IM_Pan:
          instructions = new WText( "Click drag chart left/right<br />");
          item = m_interactModeMenu->addItem( "Pan", instructions );
          break;
          
        case InteracModeIndex::IM_Select:
          instructions = new WText( "Click drag to select<br />"
                                    "to sum for spectrum.");
          item = m_interactModeMenu->addItem( "Select", instructions );
          break;
          
        case InteracModeIndex::IM_Add:
          instructions = new WText( "Click drag to add additional<br />"
                                    "time periods to spectrum sum.");
          item = m_interactModeMenu->addItem( "Add", instructions );
          break;
          
        case InteracModeIndex::IM_Remove:
          instructions = new WText( "Click drag to remove time<br />"
                                    "periods from spectrum sum." );
          item = m_interactModeMenu->addItem( "Remove", instructions );
          break;
          
        case InteracModeIndex::IM_NumMode:
          assert(0);
          break;
      }//switch( index )
      
      assert( item );
      item->clicked().connect( boost::bind(&WMenuItem::select, item) );
      item->triggered().connect( this, &D3TimeChartFilters::handleInteractionModeChange );
    }//for( loop over InteracModeIndex index )
    
    
    
    // Add menu to allow choose Foreground/Background/Secondary for the "Select", "Add", and "Remove" options
    m_specTypeDiv = new WContainerWidget( interact );
    m_specTypeDiv->setHidden( true );
    
    title = new WText( "Spectrum Type:", m_specTypeDiv );
    title->addStyleClass( "D3TimeInteractSpecTypeText" );
    
    m_specTypeSelect = new WMenu( m_specTypeDiv );
    m_specTypeSelect->addStyleClass( "D3TimeInteractMenu D3TimeInteractSpecTypeMenu LightNavMenu" );
    item = m_specTypeSelect->addItem( "Fore." );
    item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    item = m_specTypeSelect->addItem( "Back." );
    item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    item = m_specTypeSelect->addItem( "Sec." );
    item->clicked().connect( boost::bind(&WMenuItem::select, item) );
    m_specTypeSelect->select( m_specTypeSelect->itemAt(0) );
    m_specTypeSelect->itemSelected().connect( this, &D3TimeChartFilters::handleInteractionModeChange );
    
    
    interact->addWidget( instructionsStack );
    
    
    
    WContainerWidget *filterContents = new WContainerWidget();
    filterContents->addStyleClass( "D3TimeFilterTab" );
    WMenuItem *filterTabItem = tabs->addTab(filterContents, "Filter");
    
    title = new WText( "Lower energy limit:", filterContents );
    title->addStyleClass( "D3TimeInteractModeText" );
    m_lowerEnergy = new NativeFloatSpinBox( filterContents );
    
    title = new WText( "Upper energy limit:", filterContents );
    title->addStyleClass( "D3TimeInteractModeText" );
    m_upperEnergy = new NativeFloatSpinBox( filterContents );
    
    m_lowerEnergy->setSpinnerHidden( true );
    m_upperEnergy->setSpinnerHidden( true );
    
    m_lowerEnergy->setPlaceholderText( "Lower keV" );
    m_upperEnergy->setPlaceholderText( "Upper keV" );
    
    m_lowerEnergy->setText( "" );
    m_upperEnergy->setText( "" );
    
    //m_lowerEnergy->valueChanged().connect( this, &D3TimeChartFilters::energyFilterChanged );
    //m_upperEnergy->valueChanged().connect( this, &D3TimeChartFilters::energyFilterChanged );
    
    m_lowerEnergy->valueChanged().connect( m_parentChart, &D3TimeChart::userChangedEnergyRangeFilterCallback );
    m_lowerEnergy->valueChanged().connect( this, &D3TimeChartFilters::energyFilterChangedCallback );
    m_upperEnergy->valueChanged().connect( m_parentChart, &D3TimeChart::userChangedEnergyRangeFilterCallback );
    m_upperEnergy->valueChanged().connect( this, &D3TimeChartFilters::energyFilterChangedCallback );
    
    
    m_clearEnergyFilterBtn = new WPushButton( "Clear Energies", filterContents );
    m_clearEnergyFilterBtn->addStyleClass( "D3TimeFilterClear" );
    m_clearEnergyFilterBtn->clicked().connect( this, &D3TimeChartFilters::clearFilterEnergiesCallback );
    m_clearEnergyFilterBtn->hide();
    m_clearEnergyFilterBtn->setHiddenKeepsGeometry( true );
    
    
    m_dontRebin = new WCheckBox( "Don't rebin", filterContents );
    m_dontRebin->addStyleClass( "D3TimeDontRebin" );
    m_dontRebin->setToolTip( "Disable combining time samples on the chart when there are more time samples than pixels." );
    m_dontRebin->checked().connect( this, &D3TimeChartFilters::dontRebinChanged );
    m_dontRebin->unChecked().connect( this, &D3TimeChartFilters::dontRebinChanged );
    
    
    m_hideNeutrons = new WCheckBox( "Hide neutrons", filterContents );
    m_hideNeutrons->addStyleClass( "D3TimeHideNeutrons" );
    m_hideNeutrons->setToolTip( "Do not show neutrons on chart" );
    m_hideNeutrons->checked().connect( this, &D3TimeChartFilters::hideNeutronsChanged );
    m_hideNeutrons->unChecked().connect( this, &D3TimeChartFilters::hideNeutronsChanged );
    
    
    WContainerWidget *sfDiv = new WContainerWidget( filterContents );
    sfDiv->addStyleClass( "D3TimeYAxisRelScale" );
    
    sfDiv->setToolTip( "This value allows emphasizing either the gamma or neutron data by scaling"
                       " one of the y-axis's range to be a multiple of the maximum data displayed.\n"
                       " A value greater than 1.0 will scale the neutron axis by the entered"
                       " value, while not affecting the gamma y-range.  A value less than one"
                       " will scale gamma axis by one over the entered value.\n "
                       "i.e., to emphasize gamma chart enter value >1, for neutron <1.\n"
                       "Must be between 0.04 and 25." );
    WLabel *label = new WLabel( "Rel. y-max:" , sfDiv );
    
    m_gammaNeutRelEmphasis = new NativeFloatSpinBox( sfDiv );
    label->setBuddy( m_gammaNeutRelEmphasis );
    m_gammaNeutRelEmphasis->setText( "1.0" );
    m_gammaNeutRelEmphasis->setSingleStep( 0.2f );
    m_gammaNeutRelEmphasis->valueChanged().connect( this, &D3TimeChartFilters::handleGammaNeutRelEmphasisChanged );
  }//D3TimeChartFilters
  
  
  void handleInteractionModeChange()
  {
    int index = m_interactModeMenu->currentIndex();
    
    if( index < 0 ) //probably shouldnt happen, but JIC
    {
      WMenuItem *item = m_interactModeMenu->itemAt(0);
      if( !item )
        return;
    
      index = 0;
      item->select();
    }//if( index < 0 )
    
    assert( index >= 0 );
    assert( index < InteracModeIndex::IM_NumMode );
    
    bool hideSpecTypeMenu = true;
    
    switch( InteracModeIndex(index) )
    {
      case InteracModeIndex::IM_Normal:
      case InteracModeIndex::IM_Zoom:
      case InteracModeIndex::IM_Pan:
        hideSpecTypeMenu = true;
        break;
        
      case InteracModeIndex::IM_Select:
      case InteracModeIndex::IM_Add:
      case InteracModeIndex::IM_Remove:
      case InteracModeIndex::IM_NumMode:
        hideSpecTypeMenu = false;
        break;
    }//switch( index )
    
    if( hideSpecTypeMenu != m_specTypeDiv->isHidden() )
      m_specTypeDiv->setHidden( hideSpecTypeMenu );
    
    if( m_parentChart )
      m_parentChart->setUserInteractionMode( interactionMode() );
  }//void handleInteractionModeChange()
  
  
  D3TimeChart::UserInteractionMode interactionMode()
  {
    const int index = m_interactModeMenu->currentIndex();
    if( index < 0 )
    {
      auto item = m_interactModeMenu->itemAt(0);
      if( item )
        item->select();
      return D3TimeChart::UserInteractionMode::Default;
    }
    
    assert( index >= 0 );
    assert( index < InteracModeIndex::IM_NumMode );
    
    switch( index )
    {
      case InteracModeIndex::IM_Normal:
        return D3TimeChart::UserInteractionMode::Default;
        
      case InteracModeIndex::IM_Zoom:
        return D3TimeChart::UserInteractionMode::Zoom;
        
      case InteracModeIndex::IM_Pan:
        return D3TimeChart::UserInteractionMode::Pan;
      
      case InteracModeIndex::IM_Select:
        switch( m_specTypeSelect->currentIndex() )
        {
          case 0: return D3TimeChart::UserInteractionMode::SelectForeground;
          case 1: return D3TimeChart::UserInteractionMode::SelectBackground;
          case 2: return D3TimeChart::UserInteractionMode::SelectSecondary;
        }//switch( m_specTypeSelect->currentIndex() )
        break;
        
      case InteracModeIndex::IM_Add:
        switch( m_specTypeSelect->currentIndex() )
        {
          case 0: return D3TimeChart::UserInteractionMode::AddForeground;
          case 1: return D3TimeChart::UserInteractionMode::AddBackground;
          case 2: return D3TimeChart::UserInteractionMode::AddSecondary;
        }//switch( m_specTypeSelect->currentIndex() )
        break;
        
      case InteracModeIndex::IM_Remove:
        switch( m_specTypeSelect->currentIndex() )
        {
          case 0: return D3TimeChart::UserInteractionMode::RemoveForeground;
          case 1: return D3TimeChart::UserInteractionMode::RemoveBackground;
          case 2: return D3TimeChart::UserInteractionMode::RemoveSecondary;
        }//switch( m_specTypeSelect->currentIndex() )
        break;
        
      case InteracModeIndex::IM_NumMode:
        assert( 0 );
        break;
    }//switch( index )
    
    assert( 0 );
    
    return D3TimeChart::UserInteractionMode::Default;
  }//D3TimeChart::UserInteractionMode interactionMode()
  
  void clearFilterEnergiesCallback()
  {
    m_clearEnergyFilterBtn->hide();
    
    if( m_lowerEnergy->text().empty() && m_upperEnergy->text().empty() )
      return;
    
    m_lowerEnergy->setText( "" );
    m_upperEnergy->setText( "" );
    
    m_parentChart->userChangedEnergyRangeFilterCallback();
  }
  
  void resetFilterEnergies()
  {
    m_lowerEnergy->setText( "" );
    m_upperEnergy->setText( "" );
    m_clearEnergyFilterBtn->hide();
  }
  
  
  void energyFilterChangedCallback()
  {
    if( !m_lowerEnergy->text().empty() && !m_upperEnergy->text().empty() )
    {
      if( fabs(m_lowerEnergy->value() - m_upperEnergy->value()) < 0.1 )
      {
        m_lowerEnergy->setText( "" );
        m_upperEnergy->setText( "" );
        
        auto interspec = InterSpec::instance();
        const char *msg = "Energy sum range must be larger than 0.1 keV.";
        if( interspec )
          interspec->logMessage( msg, "", 3 );
      }//
    }//if( both lower and upper energies are specified )
    
    const bool emptyInput = (m_lowerEnergy->text().empty() && m_upperEnergy->text().empty());
    //m_clearEnergyFilterBtn->setEnabled( notEmpty );
    m_clearEnergyFilterBtn->setHidden( emptyInput );
    
    //const auto ranges = energyRangeFilters();
    //m_parentChart->userChangedEnergyRangeFilterCallback( ranges.first, ranges.second );
  }//void energyFilterChanged()
  
  
  pair<boost::optional<float>,boost::optional<float>> energyRangeFilters() const
  {
    pair<boost::optional<float>,boost::optional<float>> answer;
    if( !m_lowerEnergy->text().empty() )
      answer.first = m_lowerEnergy->value();
    
    if( !m_upperEnergy->text().empty() )
      answer.second = m_upperEnergy->value();
    
    if( answer.first && answer.second && ((*answer.first) > (*answer.second)) )
      std::swap( answer.first, answer.second );
    
    return answer;
  }//energyRangeFilters()
  
  
  void dontRebinChanged()
  {
    const bool dontRebin = m_dontRebin->isChecked();
    m_parentChart->setDontRebin( dontRebin );
  }
  
  void hideNeutronsChanged()
  {
    m_parentChart->setNeutronsHidden( m_hideNeutrons->isChecked() );
  }
  
  
  void handleGammaNeutRelEmphasisChanged()
  {
    float value = 1.0f;
    switch( m_gammaNeutRelEmphasis->validate() )
    {
      case Wt::WValidator::Invalid:
      case Wt::WValidator::InvalidEmpty:
        m_gammaNeutRelEmphasis->setText( "1.0" );
        break;
        
      case Wt::WValidator::Valid:
        value = m_gammaNeutRelEmphasis->value();
        if( value < 0.04f )
        {
          value = 0.04f;
          m_gammaNeutRelEmphasis->setText( "0.04" );
        }
        
        if( value > 25.0f )
        {
          value = 25.0f;
          m_gammaNeutRelEmphasis->setText( "25" );
        }
        break;
    }//switch( m_gammaNeutRelEmphasis->validate() )
    
    const string js = m_parentChart->m_jsgraph
                      +  ".setYAxisGammaNeutronRelMaxSf(" + std::to_string(value) + ");";
    m_parentChart->doJavaScript( js );
  }//handleGammaNeutRelEmphasisChanged()
  
  
  void setDontRebin( const bool dontRebin )
  {
    m_dontRebin->setChecked( dontRebin );
  }
  
  void setNeutronsHidden( const bool hide )
  {
    m_hideNeutrons->setChecked( hide );
  }
  
  void setHideNeutronsOptionVisible( const bool visible )
  {
    m_hideNeutrons->setHidden( !visible );
    m_gammaNeutRelEmphasis->setHidden( !visible );
    // Should we also reset m_gammaNeutRelEmphasis to 1.0 if we are hiding the option?
  }
};//class D3TimeChartOptions


D3TimeChart::D3TimeChart( Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_renderFlags( 0 ),
  m_layoutWidth( 0 ),
  m_layoutHeight( 0 ),
  m_chartWidthPx( 0.0 ),
  m_compactXAxis( false ),
  m_showVerticalLines( false ),
  m_showHorizontalLines( false ),
  m_dontRebin( false ),
  m_hideNeutrons( false ),
  m_spec( nullptr ),
  m_detectors_to_display(),
  m_highlights(),
  m_xAxisTitle( "Time of Measurement (seconds)"),
  m_compactXAxisTitle( "Time (seconds)" ),
  m_y1AxisTitle( "Gamma CPS" ),
  m_y2AxisTitle( "Neutron CPS"),
  m_chartClicked( this ),
  m_chartDragged( this ),
  m_displayedXRangeChange( this ),
  //m_energyRangeFilterChanged( this ),
  m_chartClickedJS( nullptr ),
  m_chartDraggedJS( nullptr ),
  m_displayedXRangeChangeJS( nullptr ),
  m_jsgraph( jsRef() + ".chart" ),
  m_chart( nullptr ),
  m_options( nullptr ),
  m_showOptionsIcon( nullptr ),
  m_gammaLineColor( 0x00, 0x00, 0x00 ),
  m_neutronLineColor( 0x00, 0x00, 0x00 ),
  m_foregroundHighlightColor( 0x00, 0x00, 0x00 ),
  m_backgroundHighlightColor( 0x00, 0x00, 0x00 ),
  m_secondaryHighlightColor( 0x00, 0x00, 0x00 ),
  m_occLineColor(),
  m_textColor( 0x00, 0x00, 0x00 ),
  m_axisColor( 0x00, 0x00, 0x00 ),
  m_chartMarginColor( 0x00, 0x00, 0x00 ),
  m_chartBackgroundColor( 0x00, 0x00, 0x00 )
{
  addStyleClass( "D3TimeChartParent" );
    
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/D3TimeChart.js" );
  wApp->useStyleSheet( "InterSpec_resources/D3TimeChart.css" );
  initChangeableCssRules();
  
  m_chart = new WContainerWidget( this );
  m_chart->addStyleClass( "D3TimeChart" );
  
  // Cancel right-click events for the div, we handle it all in JS
  m_chart->setAttributeValue( "oncontextmenu",
                    "event.cancelBubble = true; event.returnValue = false; return false;" );
  
  m_options = new D3TimeChartFilters( this );
  m_options->hide();
  
  m_showOptionsIcon = new WContainerWidget( this );
  m_showOptionsIcon->setStyleClass( "ShowD3TimeChartFilters Wt-icon" );
  m_showOptionsIcon->clicked().connect( boost::bind( &D3TimeChart::showFilters, this, true) );
}//D3TimeChart(...)


D3TimeChart::~D3TimeChart()
{
  // Do we need to cleanup any JS objects here?
}


void D3TimeChart::defineJavaScript()
{
  //WWebWidget::doJavaScript( "$(" + m_chart->jsRef() + ").append('<svg width=\"400\" height=\"75\" style=\"box-sizing: border-box; \"><rect width=\"300\" height=\"75\" style=\"fill:rgb(0,0,255);stroke-width:3;stroke:rgb(0,0,0)\" /></svg>' ); console.log( 'Added rect' );" );
  
  //WWebWidget::doJavaScript( "$(" + m_chart->jsRef() + ").append('<div style=\"box-sizing: border-box; position: absolute; border: 1px solid blue; width: 100px; height: 75px;\" />' );" );
  //cerr << "\n\nFor debugging - skipping D3TimeChart::defineJavaScript(...) - you need to fix this\n\n";
  //return;
  
  string options = "{";
  options += "xtitle: '" + (m_compactXAxis ? m_compactXAxisTitle : m_xAxisTitle) + "'";
  options += ", y1title: '" + m_y1AxisTitle + "'";
  options += ", y2title: '" + m_y2AxisTitle + "'";
  options += ", compactXAxis: " + jsbool(m_compactXAxis);
  options += ", gridx: " + jsbool(m_showVerticalLines);
  options += ", gridy: " + jsbool(m_showHorizontalLines);
  options += ", chartLineWidth: 1.0";  //ToDo: Let this be specified in C++
  options += ", dontRebin: false"; 
  options += "}";
  
  setJavaScriptMember( "chart", "new D3TimeChart(" + m_chart->jsRef() + "," + options + ");");
  
  //setJavaScriptMember( WT_RESIZE_JS,
  //                     "function(self, w, h, layout){"
  //                     " setTimeout( function(){" + m_jsgraph + ".handleResize();},0); "
  //                     "}" );
  
  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + m_chart->id() + "') )"
          + m_jsgraph + ".handleResize();"
      "}"
    "});"
  );
  
  callJavaScriptMember( "resizeObserver.observe", m_chart->jsRef() );
  
  if( !m_chartClickedJS )
  {
    m_chartClickedJS.reset( new Wt::JSignal<int,int>(m_chart, "timeclicked", false) );
    m_chartDraggedJS.reset( new Wt::JSignal<int,int,int>(m_chart, "timedragged", false) );
    m_displayedXRangeChangeJS.reset( new Wt::JSignal<int,int,int>(m_chart,"timerangechange",false) );
    
    m_chartClickedJS->connect( this, &D3TimeChart::chartClickedCallback );
    m_chartDraggedJS->connect( this, &D3TimeChart::chartDraggedCallback );
    m_displayedXRangeChangeJS->connect( this, &D3TimeChart::displayedXRangeChangeCallback );
  }//if( !m_xRangeChangedJS )
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void D3TimeChart::doJavaScript( const std::string& js )
{
  //cerr << "\n\nFor debugging - skipping D3TimeChart::doJavaScript(...) - you need to fix this\n\n";
  //return;
  
  if( isRendered() )
    WContainerWidget::doJavaScript( js );
  else
    m_pendingJs.push_back( std::move(js) );
}//doJavaScript(...)


void D3TimeChart::setData( std::shared_ptr<const SpecUtils::SpecFile> data,
                           vector<string> det_to_display )
{
  if( !data && !det_to_display.empty() )
    throw runtime_error( "D3TimeChart::setData: detectors specified with null SpecFile." );
  
  
  if( data )
  {
    if( det_to_display.empty() )
    {
      // Display all detectors if none were specified
      det_to_display = m_spec->detector_names();
    }else
    {
      // Check that all specified detectors are valid names
      const vector<string> &valid_names = data->detector_names();
      for( const string &n : det_to_display )
      {
        const auto pos = std::find( begin(valid_names), end(valid_names), n );
        if( pos == end(valid_names) )
          throw runtime_error( "D3TimeChart::setData: invalid detector name ('"
                               + n + "') specified." );
      }
    }//if( no dets specified ) / else
    
    std::sort( begin(det_to_display), end(det_to_display) );
    if( det_to_display != m_detectors_to_display )
      scheduleRenderAll();
  }//if( data )
  
  
  if( !m_highlights.empty() )
  {
    m_highlights.clear();
    scheduleHighlightRegionRender();
  }//if( !m_highlights.empty() )
  
  if( m_spec != data )
  {
    // Reset energy range summed whenever we load a new file.
    if( m_options )
      m_options->resetFilterEnergies();
    
    // Dont hide neutrons until
    if( m_hideNeutrons )
      setNeutronsHidden( false );
    
    // Schedule updating data to client
    scheduleRenderAll();
  }//if( this is a different spectrum file )
  
  m_spec = data;
  m_detectors_to_display = det_to_display;
}//void setData(...)
  


void D3TimeChart::setDataToClient()
{
  if( !m_spec )
  {
    doJavaScript( m_jsgraph +  ".setData( null );" );
    return;
  }//if( !m_spec )
  
  /** Description of JSON format sent to client JS charting
   {
     // All arrays of numbers (realTimes, sampleNumbers, and various counts) will be the same length
     // Time in seconds of each time interval being plotted.  These numbers specify the length of
     //  each time interval on the x-axis.
     realTimes: [302.4,0.1,0.1,0.9,1,...],
   
     
     // The start times of each sample is represented as the number of milliseconds since the
     //  UNIX epoch, like in JavaScript, however, to save space we will define a 'startTimeOffset'
     //  that will need to be added to each of the actual 'startTimes' to get the final time.  Also,
     //  note that JS Date assumes UTC time, so make sure when you display the time, its in UTC so
     //  no timezone offsets get added to it
     //  'startTimeOffset' and 'startTimes' may not be present if start times are not known.  Also,
     //  entries in 'startTimes' may be null if start time wasnt available for that sample.
     startTimeOffset: 1264291704078,
     startTimes: [0,604800,100,100,900, ...],
   
     // The sample number for each time interval.  These usually be monotonically increasing, but
     //  this isnt guaranteed, and the starting value isnt guaranteed either.  These values link-up
     //  to the SampleNumber of the SpecFile the data is being loaded from.
     sampleNumbers: [-1,2,3,4,5,...],
   
     // Array, that if present, will be same length as reaTimes and sampleNumbers, and gives the
     //  source-type of each sample.  The values in this array correspond to the numerical values
     //  of the SpecUtils::SourceType enum, specifically:
     //    IntrinsicActivity : 0, Calibration : 1, Background : 2, Foreground : 3, Unknown : 4.
     sourceTypes: [2,3,3,3,2,...],
   
     // The GPS coordinates of the timeslice.
     //  This array will be missing if no GPS coordinates are available, and if there isnt
     //  coordinates for a given sample, null will be provided.
     //  In the future may add timestamp as third element of array.
     gpsCoordinates: [[37.675821,-121.709863],[37.675821,-121.709863],null,[37.675821,-121.709863],
                      [37.675821,-121.709863],...],
   
     // The gamma counts to plot.  Usually we will only plot one line, but if it makes sense to
     //  plan ahead now, being able to plot multiple gamma lines would be useful (either stacked, or
     //  all lines just drawn on top of each other); if there are multiple lines than the mouse going
     //  over a line could maybe pop-up the counts value for that time interval, and detector name
     //  (unless only a single detector, then name isnt necessary).
     //  The live-time array may not be present; if it isnt, then you can assume live-times are
     //  equal to real times, and use the realTimes array.  If live-time array is present, it will
     //  be the same length as the counts array; use the live-time to divide the counts by to get
     //  counts per second for each time interval (live-time is typically a little less than
     //  real-time, but may be up to like ~90% less in extreme situations).
     //  Entries in 'counts' and/or 'liveTimes' arrays may be null if that information is not
     //  available for the corresponding sample number.
     gammaCounts: [{detName: 'Aa1', color: 'rgb(13,55,19)', counts: [1022,974,333,827,921,...],
                    liveTimes: [300,0.1,0.1,0.99,0.9,...]},
                   {detName: 'Ba1', color: 'green', counts: [55,60,18,99,1000,...],
                    liveTimes: [299,0.092,0.1,1,1,...]}
                  ],
   
     // The neutron counts to plot, similar to gamma, but the field may not exist, or be null, in
     //  which case there will be no neutron counts to plot, and the Y2-axis should disapear.
     //  Notes under gammaCounts above for the live-times array also apply here
     neutronCounts: [{detName: 'Aa1n', color: '#145366', counts: [1,0,4,10,12,...],
                     liveTimes: [300,0.1,0.1,0.1,0.1,...]},
                     {detName: 'VD2n', color: '#ff0000', counts: [8,12,13,19,0,...]}],
   
     // The sample number ranges for "occupied" and non-occupied times of the data.
     //  This field may not be present, may be null, may be empty array, or have entries as
     //  described.  This field is typically used for describing radiation portal monitor data, and
     //  will be absent for portable detection systems.
     occupancies: [{
                    // Says whether this sample range corresponds to being "occupied" (corresponding
                    //  to (SpecUtils::OccupancyStatus::Occupied) or if false, "not occupied"
                    //  corresponding to SpecUtils::OccupancyStatus::NotOccupied).  If occupancy
                    //  status for a sample is unknown or undefined, then that sample will not be
                    //  included in any of these ranges.
                    //  Typically non-occupied intervals will also be marked as Background in the
                    //  'sourceTypes' array, and occupied intervals will be marked as Foreground.
                    status: true,
   
                    // The sample number corresponding to the first occupied sample
                    // In the previous implementation a vertical line was drawn, on left side of
                    // sample with the text "occ. start"
                    startSample: 4,
        
                    // The last occupied sample number.  Was previously drawn as a vertical line
                    //  with text "occ. end"
                    endSample: 10,
                   
                    // Color to draw vertical line and text
                    color: 'grey'
                   },
                   {...}
                  ]
   }
   */
  
  //Variable to control if we will only plot a single gamma and neutron line, or if we will plot
  //  each detector separately.  This may become a user option at some point, or the idea of more
  //  then one line for gamma/neutron may get scrapped if it is to confusing or unhelpful.
  const bool plotDetectorsSeparate = false;
  static_assert( !plotDetectorsSeparate, "Plotting detectors separate JS functionality has not been tested." );
  
  vector<double> realTimes;
  vector<int> sampleNumbers;
  vector<SpecUtils::SourceType> sourceTypes;
  vector<std::tuple<double,double,boost::posix_time::ptime>> gpsCoordinates;
  int64_t startTimesOffset = 0;
  vector<int64_t> startTimes;
  bool allUnknownSourceType = true, anyStartTimeKnown = false;
  bool haveAnyGps = false, plottingNeutrons = false;
  
  map<string,bool> hasGamma, hasNuetron;
  map<string,vector<double>> gammaCounts, neutronCounts, liveTimes;  //maps from detector name to counts
  
  const set<int> &sample_numbers = m_spec->sample_numbers();
  //const vector<string> &detNames = m_spec->detector_names();
  const vector<string> &detNames = m_detectors_to_display;
  
#if( PERFORM_DEVELOPER_CHECKS )
  {// begin check that detector names to plot are valid
    const vector<string> &all_det_names = m_spec->detector_names();
    for( const string &name : detNames )
    {
      const auto pos = std::find( begin(all_det_names), end(all_det_names), name );
      if( pos == end(all_det_names) )
      {
        log_developer_error( __func__, ("Found invalid detector name, '" + name + "'").c_str() );
      }
      
      assert( pos != end(all_det_names) );
    }//for( const string &name : detNames )
  }// end check that detector names to plot are valid
#endif
  

  boost::optional<float> lowerEnergy, upperEnergy;
  if( m_options )
  {
    const auto energyRange = m_options->energyRangeFilters();
    lowerEnergy = energyRange.first;
    upperEnergy = energyRange.second;
  }//if( m_options )
  
#define Q_DBL_NaN std::numeric_limits<double>::quiet_NaN()
  
  for( const int sample_num : sample_numbers )
  {
    /// \TODO: check that all Measurements have same real times; right now we'll just use the max
    ///        value for each sample
    float realTime = 0.0f;
    boost::posix_time::ptime startTime;
    SpecUtils::SourceType sourcetype = SpecUtils::SourceType::Unknown;
    tuple<double,double,boost::posix_time::ptime> coords{ Q_DBL_NaN, Q_DBL_NaN, {} };
    
    // Lets make sure this sample number has data from any of the currently selected data, and if
    //  not, we wont include this sample.
    bool haveAnyDataThisSample = false;
    for( const string &detName : detNames )
    {
      if( m_spec->measurement( sample_num, detName ) )
      {
        haveAnyDataThisSample = true;
        break;
      }
    }//for( const string &detName : detNames )
    
    if( !haveAnyDataThisSample )
      continue;
    
    sampleNumbers.push_back( sample_num );
    
    for( const string &detName : detNames )
    {
      const auto m = m_spec->measurement( sample_num, detName );
      
      if( !m )
      {
        liveTimes[detName].push_back( Q_DBL_NaN );
        gammaCounts[detName].push_back( Q_DBL_NaN );
        neutronCounts[detName].push_back( Q_DBL_NaN );
        
        if( !hasGamma.count(detName) )
          hasGamma[detName] = false;
        
        if( !hasNuetron.count(detName) )
          hasNuetron[detName] = false;
        
        continue;
      }//if( !m )
      
      realTime = std::max( realTime, m->real_time() );
      if( startTime.is_special() )
        startTime = m->start_time();
      
      auto gammas = m->gamma_counts();
      if( gammas && gammas->size() )
        hasGamma[detName] = true;
      else if( !hasGamma.count(detName) )
        hasGamma[detName] = false;
               
      if( m->contained_neutron() )
        hasNuetron[detName] = true;
      else if( !hasNuetron.count(detName) )
        hasNuetron[detName] = false;
               
      if( m->source_type() != SpecUtils::SourceType::Unknown )
        sourcetype = m->source_type();
      
      if( lowerEnergy || upperEnergy )
      {
        const float specMinEnergy = m->gamma_energy_min();
        const float specMaxEnergy = m->gamma_energy_max();
        
        double gamma_sum = m->gamma_count_sum();
        
        if( (!lowerEnergy || (lowerEnergy.get() < specMinEnergy) )
           && (!upperEnergy || (upperEnergy.get() > specMaxEnergy)) )
        {
          // gamma_sum = m->gamma_count_sum();
        }else if( lowerEnergy && upperEnergy )
        {
          gamma_sum = m->gamma_integral(*lowerEnergy, *upperEnergy);
        }else if( lowerEnergy )
        {
          gamma_sum = m->gamma_integral(*lowerEnergy, specMaxEnergy + 1000);
        }else if( upperEnergy )
        {
          gamma_sum = m->gamma_integral( specMinEnergy - 1000, *upperEnergy);
        }
        
        gammaCounts[detName].push_back( gamma_sum );
      }else
      {
        gammaCounts[detName].push_back( m->gamma_count_sum() );
      }
      neutronCounts[detName].push_back( m->neutron_counts_sum() );
      liveTimes[detName].push_back( m->live_time() );
      
      if( m->has_gps_info() )
      {
        haveAnyGps = true;
        std::get<0>(coords) = m->latitude();
        std::get<1>(coords) = m->longitude();
        std::get<2>(coords) = m->position_time();
      }
    }//for( const string &name : m_spec->detector_names() )
    
    realTimes.push_back( realTime );
    sourceTypes.push_back( sourcetype );
    allUnknownSourceType = (allUnknownSourceType && (sourcetype == SpecUtils::SourceType::Unknown));
    
    if( startTime.is_special() )
    {
      startTimes.push_back( std::numeric_limits<int64_t>::min() );
    }else
    {
      anyStartTimeKnown = true;
       
      static const boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
      
      // std::time_t (what boost uses for seconds type) may be 32-bit, so convert to int64_t to
      //  avoid overflow.
      int64_t nmilli = (startTime - epoch).total_seconds();
      nmilli *= 1000;
      
      const double frac = startTime.time_of_day().fractional_seconds()
                        / static_cast<double>(boost::posix_time::time_duration::ticks_per_second());
      nmilli += static_cast<int64_t>( std::round(1000*frac) );
      
      if( startTimesOffset <= 0 && nmilli > 0 )
        startTimesOffset = nmilli;
      
      startTimes.push_back( nmilli - startTimesOffset );
    }
    
    gpsCoordinates.push_back( coords );
  }//for( loop over sample numbers )
  
  const size_t numSamples = sampleNumbers.size();
  //A quick sanity check all the arrays will be the same length.
  assert( numSamples == realTimes.size() );
  assert( numSamples == sourceTypes.size() );
  assert( numSamples == startTimes.size() );
  assert( numSamples == gpsCoordinates.size() );
  
  for( const auto &p : gammaCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  for( const auto &p : neutronCounts )
  {
    assert( p.second.size() == numSamples );
  }
  
  for( const auto &p : liveTimes )
  {
    assert( p.second.size() == numSamples );
  }
  
  if( !numSamples )
  {
    doJavaScript( m_jsgraph +  ".setData( null );" );
    return;
  }//if( !numSamples )
  
  const char jssep[] = { '[', ',' };
  
  auto printNumberArray = [&jssep]( WStringStream &js, const vector<double> &arr ){
    
    if( arr.empty() )
    {
      js << "[]";
      return;
    }//if( !arr.empty() )
    
    for( size_t i = 0; i < arr.size(); ++i )
    {
      if( IsNan(arr[i]) )
        js << jssep[i ? 1 : 0] << "null";
      else
        js << jssep[i ? 1 : 0] << arr[i];
    }//for( loop over realTimes )
    
    js << "]";
  };//printNumberArray
  
  
  WStringStream js;
  js << m_jsgraph <<  ".setData( {\n";
  
  js << "\t\"realTimes\": ";
  printNumberArray( js, realTimes );
  
  if( anyStartTimeKnown )
  {
    //doubles have 53 bits of integer precision - should be fine
    // long long int should be at least 64 bits, and otherwise int64_t may not be supported
    // by WStringStream.
    assert( static_cast<long long>(startTimesOffset) == (startTimesOffset) );
    js << ",\n\t\"startTimeOffset\": " << static_cast<long long>(startTimesOffset);
    js << ",\n\t\"startTimes\": ";
    for( size_t i = 0; i < startTimes.size(); ++i )
    {
      js << jssep[i ? 1 : 0];
      if( startTimes[i] == std::numeric_limits<int64_t>::min() )
        js << "null";
      else
        js << static_cast<long long>(startTimes[i]);
    }
    js << "]";
  }//if( anyStartTimeKnown )
  
  js << ",\n\t\"sampleNumbers\": ";
  for( size_t i = 0; i < sampleNumbers.size(); ++i )
    js << jssep[i ? 1 : 0] << sampleNumbers[i];
  js << "]";
  
  if( !allUnknownSourceType )
  {
    js << ",\n\t\"sourceTypes\": ";
    for( size_t i = 0; i < sourceTypes.size(); ++i )
      js << jssep[i ? 1 : 0] << static_cast<int>(sourceTypes[i]);
    js << "]";
  }//if( !allUnknownSourceType )
  
  if( lowerEnergy )
    js << ",\n\t\"filterLowerEnergy\": " << static_cast<double>(*lowerEnergy);
  
  if( upperEnergy )
    js << ",\n\t\"filterUpperEnergy\": " << static_cast<double>(*upperEnergy);
    
  if( haveAnyGps )
  {
    js << ",\n\t\"gpsCoordinates\": ";
    for( size_t i = 0; i < gpsCoordinates.size(); ++i )
    {
      js << jssep[i ? 1 : 0];
      
      const auto &coords = gpsCoordinates[i];
      if( IsNan(std::get<0>(coords)) )
        js << "null";
      else
        js << "[" << std::get<0>(coords) << "," << std::get<1>(coords) << "]";
    }//for( size_t i = 0; i < gpsCoordinates.size(); ++i )
    
    js << "]";
  }//if( haveAnyGps )
  
  
  if( plotDetectorsSeparate )
  {
    int nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasGamma[detName] )
        continue;
    
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"gammaCounts\": [" )
         << "\n\t\t{\"detName\": \"" << detName << "\", \"color\": \""
         << (m_gammaLineColor.isDefault() ? string("#cfced2") :  m_gammaLineColor.cssText())
         << "\", \"counts\": ";
      
      const auto &counts = gammaCounts[detName];
      printNumberArray( js, counts );
      
      
      const auto &lts = liveTimes[detName];
      js << ",\n\t\"liveTimes\": ";
      printNumberArray( js, lts );
      
      js << "}";
    }//for( const auto &detName : detNames )
    
    if( nwrote )
      js << "]";
    
    nwrote = 0;
    for( const auto &detName : detNames )
    {
      if( !hasNuetron[detName] )
        continue;
    
      plottingNeutrons = true;
      
      //For now we'll just make all detector lines the same color (if we are even doing multiple lines)
      js << string(nwrote++ ? "," : ",\n\t\"neutronCounts\": [" ) << "\n\t\t{\"detName\": \""
         << detName << "\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("rgb(0,128,0)") :  m_neutronLineColor.cssText())
         << "\", \"counts\": ";

      const auto &counts = neutronCounts[detName];
      printNumberArray( js, counts );
      
      js << "}";
    }//for( const auto &detName : detNames )
    
    if( nwrote )
      js << "]";
  }else
  {
    bool haveAnyGamma = false, haveAnyNeutron = false;
    for( const auto &p : hasGamma )
      haveAnyGamma |= p.second;
    for( const auto &p : hasNuetron )
      haveAnyNeutron |= p.second;
      
    if( haveAnyGamma )
    {
      js << ",\n\t\"gammaCounts\": [{\"detName\": \"\", \"color\": \""
         << (m_gammaLineColor.isDefault() ? string("#cfced2") :  m_gammaLineColor.cssText())
         << "\",\n\t\t\"counts\": [";
      
      for( size_t i = 0; i < numSamples; ++i )
      {
        double sum = 0.0;
        bool anyNonNan = false;
        for( const auto &p : gammaCounts )
        {
          assert( p.second.size() == numSamples );
          
          if( !IsNan(p.second[i]) )
          {
            anyNonNan = true;
            sum += p.second[i];
          }
        }//for( const auto &p : gammaCounts )
        
        double liveTime = 0.0;
        bool haveLiveTime = false;
        for( const auto &p : liveTimes )
        {
          if( !IsNan(p.second[i]) )
          {
            haveLiveTime = true;
            liveTime += p.second[i];
          }
        }
        
        
        js << string(i ? "," : "");
        if( anyNonNan )
          js << sum;
        else
          js << "null";
      }//for( size_t i = 0; i < numSamples; ++i )
      js << "]";
      
      js << ",\n\t\t\"liveTimes\": [";
      for( size_t i = 0; i < numSamples; ++i )
      {
        double liveTime = 0.0;
        bool haveLiveTime = false;
        for( const auto &p : liveTimes )
        {
          if( !IsNan(p.second[i]) )
          {
            haveLiveTime = true;
            liveTime += p.second[i];
          }
        }
        
        js << string(i ? "," : "");
        if( haveLiveTime )
          js << liveTime;
        else
          js << "null";
      }//for( size_t i = 0; i < numSamples; ++i )
      js << "]";
      
      js << "\n\t}]";
    }//if( haveAnyGamma )
    
    if( haveAnyNeutron )
    {
      plottingNeutrons = true;
      
      js << ",\n\t\"neutronCounts\": [{\"detName\": \"\", \"color\": \""
         << (m_neutronLineColor.isDefault() ? string("#cfced2") :  m_neutronLineColor.cssText())
         << "\", \"counts\": [";
      
      vector<double> neutronLiveTimes( numSamples, std::numeric_limits<double>::quiet_NaN() );
      for( size_t i = 0; i < numSamples; ++i )
      {
        bool haveNonNan = false;
        double sum = 0.0;
        for( const auto &p : neutronCounts )
        {
          if( !IsNan(p.second[i]) )
          {
            haveNonNan = true;
            sum += p.second[i];
            
            if( IsNan(neutronLiveTimes[i]) )
              neutronLiveTimes[i] = 0;
            if( !IsNan(realTimes[i]) )
              neutronLiveTimes[i] += realTimes[i];
          }//if( we have neutron counts )
        }//for( loop over neutron detectors to sum their counts )
        js << string(i ? "," : "");
        if( haveNonNan )
          js << sum;
        else
          js << "null";
      }//for( loop over samples )
      
      js << "],\n\t\t\"liveTimes\": ";
      printNumberArray( js, neutronLiveTimes );
      
      js << "\n\t}]";
    }//if( haveAnyNeutron )
  }//if( plotDetectorsSeparate ) / else
  
  auto occRanges = sampleNumberRangesWithOccupancyStatus( SpecUtils::OccupancyStatus::Occupied, m_spec );
  if( occRanges.size() )
  {
    js << ",\n\t\"occupancies\": [";
    for( size_t i = 0; i < occRanges.size(); ++i )
    {
      js << string(i ? ",\n\t\t" : "\n\t\t")
         << "{ startSample: " << occRanges[i].first
         << ", endSample: " << occRanges[i].second
         << ", color: \""
         << (m_occLineColor.isDefault() ? string("rgb(128,128,128)") :  m_occLineColor.cssText())
         << "\" }";
    }//for( size_t i = 0; i < occRanges.size(); ++i )
    js << "\n\t]";
  }//if( occRanges.size() )
  
  js << "\n\t} );";
  
  // cout << "\n\nWill set time data with JSON=" + js.str() + "\n\n" << endl;

  doJavaScript( js.str() );
  
  
  
  // If neutrons will be displayed, we need to update the location of the filter icon so it
  //  doesnt overlap with the neutron y-axis (the y-axis on right of chart) too bad.
  const WString neutOpSc( "HasNeutrons" );
  const bool showNeutron = plottingNeutrons && !m_hideNeutrons;
  if( m_showOptionsIcon && (showNeutron != m_showOptionsIcon->hasStyleClass(neutOpSc)) )
  {
    if( showNeutron )
      m_showOptionsIcon->addStyleClass(neutOpSc);
    else
      m_showOptionsIcon->removeStyleClass(neutOpSc);
  }
  
  if( m_options )
    m_options->setHideNeutronsOptionVisible( plottingNeutrons );
}//void setDataToClient()


void D3TimeChart::setHighlightedIntervals( const std::set<int> &sample_numbers,
                                           const SpecUtils::SpectrumType type )
{
  const size_t npre = m_highlights.size();
  m_highlights.erase( std::remove_if( begin(m_highlights), end(m_highlights),
    [type](const D3TimeChart::HighlightRegion &region) -> bool {
      return (region.type == type);
  }), end(m_highlights) );
 
  if( sample_numbers.empty() )
  {
    if( npre != m_highlights.size() )
      scheduleHighlightRegionRender();
    return;
  }
  
  auto interspec = InterSpec::instance();
  if( !m_spec )
  {
    const char *msg = "Time chart passed sample numbers when no data is set - should not happen";
    if( interspec )
      interspec->logMessage( msg, "", 3 );
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg );
#endif
    return;
  }//if( !m_spec )

#if( PERFORM_DEVELOPER_CHECKS )
  for( const int sample : sample_numbers )
  {
    if( !m_spec->sample_numbers().count(sample) )
    {
      log_developer_error( __func__, ("Was passed an invalid samplenumber: "
                                       + std::to_string(sample)).c_str() );
    }//if( spec file doesnt have the input sample number )
  }//for( loop over all incomming sample numbers )
#endif
  
  const set<int> allSamplesThisType = interspec->sampleNumbersForTypeFromForegroundFile(type);
  if( sample_numbers == allSamplesThisType )
  {
    scheduleHighlightRegionRender();
    return;
  }//
  
  shared_ptr<const ColorTheme> theme = (interspec ? interspec->getColorTheme() : nullptr);
  
  D3TimeChart::HighlightRegion region;
  region.type = type;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      region.color = theme ? theme->timeHistoryForegroundHighlight : Wt::WColor(255,255,0,155);
      break;
    case SpecUtils::SpectrumType::SecondForeground:
      region.color = theme ? theme->timeHistorySecondaryHighlight : Wt::WColor(0,255,255,75);
      break;
    case SpecUtils::SpectrumType::Background:
      region.color = theme ? theme->timeHistoryBackgroundHighlight : Wt::WColor(0,128,0,75);
      break;
  }//switch( type )
  
  int firstInRange = *(sample_numbers.begin());
  int previous = firstInRange;
  
  for( auto iter = begin(sample_numbers); iter != end(sample_numbers); ++iter )
  {
    const int thisval = *iter;
      
    if( (thisval > (previous+1)) )
    {
      region.start_sample_number = firstInRange;
      region.end_sample_number = previous;
      m_highlights.push_back( region );
      
      firstInRange = thisval;
    }//if( thisval > (previous+1) )
      
    previous = thisval;
  }//for( loop over smaple_numbers )
    
  region.start_sample_number = firstInRange;
  region.end_sample_number = previous;
  m_highlights.push_back( region );
  
  scheduleHighlightRegionRender();
}//setHighlightedIntervals(...)


Wt::Signal<int,Wt::WFlags<Wt::KeyboardModifier>> &D3TimeChart::chartClicked()
{
  return m_chartClicked;
}


Wt::Signal<int,int,Wt::WFlags<Wt::KeyboardModifier>> &D3TimeChart::chartDragged()
{
  return m_chartDragged;
}


Wt::Signal<int,int,int> &D3TimeChart::displayedXRangeChange()
{
  return m_displayedXRangeChange;
}


//Wt::Signal<boost::optional<float>,boost::optional<float>> &D3TimeChart::energyRangeFilterChanged()
//{
//  return m_energyRangeFilterChanged;
//}


vector<pair<int,int>> D3TimeChart::sampleNumberRangesWithOccupancyStatus(
                                                const SpecUtils::OccupancyStatus wanted_occ_status,
                                                std::shared_ptr<const SpecUtils::SpecFile> spec )
{
  vector< pair<int,int> > occ_ranges;
  
  if( !spec )
    return occ_ranges;
  
  std::set<int> wantedSamples;
  const set<int> &allsamples = spec->sample_numbers();
  const vector<string> &detnames = spec->detector_names();
  
  bool inWantedRegion = false;
  int startSample = numeric_limits<int>::min(), prevSample = numeric_limits<int>::min();
    
  for( const int sample : allsamples )
  {
    bool wantSample = false;
    for( size_t i = 0; !wantSample && i < detnames.size(); ++i )
    {
      auto m = spec->measurement( sample, detnames[i] );
      wantSample = (m && (m->occupied() == wanted_occ_status));
    }
    
    if( wantSample )
      wantedSamples.insert( sample );
    
    if( wantSample && !inWantedRegion )
      startSample = sample;
    else if( inWantedRegion && !wantSample )
      occ_ranges.push_back( {startSample,prevSample} );
    
    prevSample = sample;
    inWantedRegion = wantSample;
  }//for( const int sample : allsamples )
  
  if( inWantedRegion )
    occ_ranges.push_back( {startSample,prevSample} );

  if( wantedSamples == spec->sample_numbers() )
    return {};
  
  return occ_ranges;
}//sampleNumberRangesWithOccupancyStatus(...)


void D3TimeChart::setHighlightRegionsToClient()
{
  if( !m_spec || m_highlights.empty() )
  {
    // Dont set the highlight regions if no data.  We probably shouldnt ever get here anyway.
    doJavaScript( m_jsgraph + ".setHighlightRegions(null);" );
    return;
  }//if( !m_spec )
  
  auto type_to_str = []( const SpecUtils::SpectrumType type ) -> std::string {
    switch (type )
    {
      case SpecUtils::SpectrumType::Foreground:       return "FOREGROUND";
      case SpecUtils::SpectrumType::SecondForeground: return "SECONDARY";
      case SpecUtils::SpectrumType::Background:       return "BACKGROUND";
    }
    return "";
  };//type_to_str
  
  
  WStringStream js;
  js << m_jsgraph << ".setHighlightRegions( [";
  
  for( size_t i = 0; i < m_highlights.size(); ++i )
  {
    const auto &region = m_highlights[i];
    js << string(i ? "," : "")
       << "{startSample: " << region.start_sample_number
       << ", endSample: " << region.end_sample_number
       << ", fillColor: \"" << region.color.cssText() << "\""
       << ", type: \"" << type_to_str(region.type) << "\"}";
  }//for( loop over highlight regions )
  js << "] );";
  
  cout << "\n\nWill set highlight regions with JSON=" + js.str() + "\n\n" << endl;
  doJavaScript( js.str() );
}//setHighlightRegionsToClient()


void D3TimeChart::saveChartToPng( const std::string &filename )
{
  // \TODO: implement - see D3SpectrumDisplayDiv for a starting point, although it isnt finished
}//void saveChartToPng( const std::string &filename )


void D3TimeChart::scheduleRenderAll()
{
  m_renderFlags |= TimeRenderActions::UpdateData;
  m_renderFlags |= TimeRenderActions::UpdateHighlightRegions;
  
  scheduleRender();
}//scheduleRenderAll()


void D3TimeChart::scheduleHighlightRegionRender()
{
  m_renderFlags |= TimeRenderActions::UpdateHighlightRegions;
  
  scheduleRender();
}//void scheduleHighlightRegionRender();


void D3TimeChart::applyColorTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return;
  
  setGammaLineColor( theme->timeChartGammaLine );
  setNeutronLineColor( theme->timeChartNeutronLine );
  setAxisLineColor( theme->timeAxisLines );
  setChartBackgroundColor( theme->timeChartBackground );
  setChartMarginColor( theme->timeChartMargins );
  setTextColor( theme->timeChartText );
  m_occLineColor = theme->occupancyIndicatorLines;
  m_foregroundHighlightColor = theme->timeHistoryForegroundHighlight;
  m_backgroundHighlightColor = theme->timeHistoryBackgroundHighlight;
  m_secondaryHighlightColor = theme->timeHistorySecondaryHighlight;
  
  scheduleRenderAll();
}//void applyColorTheme( std::shared_ptr<const ColorTheme> theme );


void D3TimeChart::setGammaLineColor( const Wt::WColor &color )
{
  m_gammaLineColor = color.isDefault() ? WColor( 0x00, 0x00, 0x00 ) : color;
  scheduleRenderAll();
}

void D3TimeChart::setNeutronLineColor( const Wt::WColor &color )
{
  m_neutronLineColor = color.isDefault() ? WColor(0x00,0xff,0xff) : color;
  scheduleRenderAll();
}

void D3TimeChart::setTextColor( const Wt::WColor &color )
{
  m_textColor = color.isDefault() ? WColor(0,0,0) : color;
  const string c = m_textColor.cssText();
  
  const string rulename = "TimeTextColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ".D3TimeChart .xaxistitle,"
                                       " .D3TimeChart .yaxistitle,"
                                       " .D3TimeChart .yaxis,"
                                       " .D3TimeChart .yaxislabel,"
                                       " .D3TimeChart .xaxis",
                                       "fill: " + c );
}//setTextColor(...)


void D3TimeChart::setAxisLineColor( const Wt::WColor &color )
{
  m_axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  string rulename = "TimeAxisColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ".D3TimeChart .xaxis > .domain,"
                                        " .D3TimeChart .yaxis > .domain,"
                                        " .D3TimeChart .xaxis > .tick > line,"
                                        " .D3TimeChart .yaxis > .tick > line,"
                                        " .D3TimeChart .yaxistick",
                                       "stroke: " + m_axisColor.cssText() );
}//setAxisLineColor(...)


void D3TimeChart::setChartMarginColor( const Wt::WColor &color )
{
  m_chartMarginColor = color;
  
  const string rulename = "TimeMarginColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  if( color.isDefault() )
  {
    if( m_cssRules.count(rulename) )
    {
      style.removeRule( m_cssRules[rulename] );
      m_cssRules.erase( rulename );
    }
    
    return;
  }//if( color.isDefault() )
  
  //Actualkly this will set the background for the entire chart...
  m_cssRules[rulename] = style.addRule( "#" + id() + " > svg", "background: " + color.cssText() );
}//setChartMarginColor(...)


void D3TimeChart::setChartBackgroundColor( const Wt::WColor &color )
{
  m_chartBackgroundColor = color;
  const string c = color.isDefault() ? "rgba(0,0,0,0)" : color.cssText();
  
  const string rulename = "TimeBackgroundColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  if( color.isDefault() )
  {
    if( m_cssRules.count(rulename) )
      style.removeRule( m_cssRules[rulename] );
  }//if( color.isDefault() )
  
  m_cssRules[rulename] = style.addRule( "#chartarea" + id(), "fill: " + c );
}//setChartBackgroundColor(...)


void D3TimeChart::setXAxisTitle( const std::string &title, const std::string &compactTitle )
{
  m_xAxisTitle = title;
  SpecUtils::ireplace_all( m_xAxisTitle, "'", "&#39;" );
  
  m_compactXAxisTitle = compactTitle;
  SpecUtils::ireplace_all( m_compactXAxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle('"
                  + (m_compactXAxis ? m_compactXAxisTitle : m_xAxisTitle) + "');" );
}//void setXAxisTitle( const std::string &title )


void D3TimeChart::setY1AxisTitle( const std::string &title )
{
  m_y1AxisTitle = title;
  SpecUtils::ireplace_all( m_y1AxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setY1AxisTitle('" + m_y1AxisTitle + "');" );
}//void setY1AxisTitle( const std::string &title )


void D3TimeChart::setY2AxisTitle( const std::string &title )
{
  m_y2AxisTitle = title;
  SpecUtils::ireplace_all( m_y2AxisTitle, "'", "&#39;" );
  
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setY2AxisTitle('" + m_y2AxisTitle + "');" );
}//void setY2AxisTitle( const std::string &title )


void D3TimeChart::showFilters( const bool show )
{
  if( !m_options || !m_showOptionsIcon || (m_showOptionsIcon->isHidden() == show) )
    return;
  
  m_showOptionsIcon->setHidden( show );
  m_options->setHidden( !show );
  
  if( !show )
  {
    //I'm not the fence if we should remove the filters when hiding the filter widget...
    //const auto filterRange = m_options->energyRangeFilters();
    //if( filterRange.first || filterRange.second )
    //{
    //  m_options->resetFilterEnergies();
    //  scheduleRenderAll();
    //}
  }//if( !show )
  
  //doJavaScript( "setTimeout( function(){" + m_jsgraph + ".handleResize();},0); " );
  
  const auto mode = show ? m_options->interactionMode() : UserInteractionMode::Default;
  setUserInteractionMode( mode );
}//void showFilters( const bool show )


void D3TimeChart::setUserInteractionMode( const UserInteractionMode mode )
{
  string modestr = "";
  switch( mode )
  {
    case UserInteractionMode::Default:          modestr = "Default";          break;
    case UserInteractionMode::Zoom:             modestr = "Zoom";             break;
    case UserInteractionMode::Pan:              modestr = "Pan";              break;
    case UserInteractionMode::SelectForeground: modestr = "SelectForeground"; break;
    case UserInteractionMode::SelectBackground: modestr = "SelectBackground"; break;
    case UserInteractionMode::SelectSecondary:  modestr = "SelectSecondary";  break;
    case UserInteractionMode::AddForeground:    modestr = "AddForeground";    break;
    case UserInteractionMode::AddBackground:    modestr = "AddBackground";    break;
    case UserInteractionMode::AddSecondary:     modestr = "AddSecondary";     break;
    case UserInteractionMode::RemoveForeground: modestr = "RemoveForeground"; break;
    case UserInteractionMode::RemoveBackground: modestr = "RemoveBackground"; break;
    case UserInteractionMode::RemoveSecondary:  modestr = "RemoveSecondary";  break;
  }//switch( mode )
  
  doJavaScript( m_jsgraph + ".setUserInteractionMode('" + modestr + "');"  );
}//void setUserInteractionMode( const UserInteractionMode mode )


//void D3TimeChart::userChangedEnergyRangeFilterCallback( const boost::optional<float> lowerEnergy,
//                                                const boost::optional<float> upperEnergy )
void D3TimeChart::userChangedEnergyRangeFilterCallback()
{
  //m_energyRangeFilterChanged.emit( lowerEnergy, upperEnergy );
  
  scheduleRenderAll();
}//void userChangedEnergyRangeFilterCallback( const float lowerEnergy, const float upperEnergy )


//pair<boost::optional<float>,boost::optional<float>> D3TimeChart::energyRangeFilters()
//{
//  if( m_options )
//    return m_options->energyRangeFilters();
//  return pair<boost::optional<float>,boost::optional<float>>();
//}


void D3TimeChart::setCompactAxis( const bool compact )
{
  if( isRendered() )
  {
    if( m_compactXAxis != compact )
      doJavaScript( m_jsgraph + ".setXAxisTitle('"
                   + (compact ? m_compactXAxisTitle : m_xAxisTitle) + "');" );
    
    doJavaScript( m_jsgraph + ".setCompactXAxis(" + jsbool(compact) + ");" );
  }//if( isRendered() )
  
  m_compactXAxis = compact;
}//void setCompactAxis( compact )


bool D3TimeChart::isAxisCompacted() const
{
  return m_compactXAxis;
}

void D3TimeChart::showGridLines( bool show )
{
  m_showVerticalLines = show;
  m_showHorizontalLines = show;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(show) + ");"
                  + m_jsgraph + ".setGridY(" + jsbool(show) + ");" );
}

void D3TimeChart::showVerticalLines( const bool draw )
{
  m_showVerticalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridX(" + jsbool(draw) + ");" );
}

void D3TimeChart::showHorizontalLines( const bool draw )
{
  m_showHorizontalLines = draw;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setGridY(" + jsbool(draw) + ");" );
}

bool D3TimeChart::verticalLinesShowing() const
{
  return m_showVerticalLines;
}

bool D3TimeChart::horizontalLinesShowing() const
{
  return m_showHorizontalLines;
}


void D3TimeChart::setDontRebin( const bool dontRebin )
{
  m_dontRebin = dontRebin;
  if( m_options )
    m_options->setDontRebin( dontRebin );
  doJavaScript( m_jsgraph + ".setDontRebin(" + jsbool(dontRebin) +  ");" );
}


bool D3TimeChart::dontRebin() const
{
  return m_dontRebin;
}


void D3TimeChart::setNeutronsHidden( const bool hide )
{
  m_hideNeutrons = hide;
  if( m_options )
    m_options->setNeutronsHidden( hide );
  doJavaScript( m_jsgraph + ".setNeutronsHidden(" + jsbool(hide) +  ");" );
}


bool D3TimeChart::neutronsHidden() const
{
  return m_hideNeutrons;
}

void D3TimeChart::setXAxisRangeSamples( const int min_sample_num, const int max_sample_num )
{
  doJavaScript( m_jsgraph + ".setXAxisZoomSamples("
                + std::to_string(min_sample_num) + "," + std::to_string(max_sample_num) + ");" );
}


void D3TimeChart::initChangeableCssRules()
{
  WCssStyleSheet &style = wApp->styleSheet();

  m_cssRules["TimeGridColor"] = style.addRule( " .D3TimeChart .xgrid > .tick,  .D3TimeChart .ygrid > .tick", "stroke: #b3b3b3" );
  m_cssRules["TimeMinorGridColor"] = style.addRule( " .D3TimeChart .minorgrid", "stroke: #e6e6e6" );
}//void initChangeableCssRules()


void D3TimeChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  //const bool renderUpdate = (flags & Wt::RenderFlag::RenderUpdate);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
  
  if( m_renderFlags.testFlag(TimeRenderActions::UpdateData) )
    setDataToClient();
  
  if( m_renderFlags.testFlag(TimeRenderActions::UpdateHighlightRegions) )
    setHighlightRegionsToClient();
  
  m_renderFlags = 0;
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void D3TimeChart::chartClickedCallback( int sample_number, int modifier_keys )
{
  m_chartClicked.emit(sample_number, Wt::WFlags<Wt::KeyboardModifier>(Wt::KeyboardModifier(modifier_keys)) );
}//chartClickedCallback(...)


void D3TimeChart::chartDraggedCallback( int first_sample_number, int last_sample_number, int modifier_keys )
{
  m_chartDragged.emit( first_sample_number, last_sample_number, Wt::WFlags<Wt::KeyboardModifier>(Wt::KeyboardModifier(modifier_keys)) );
}//chartDraggedCallback(...)


void D3TimeChart::displayedXRangeChangeCallback( int first_sample_number, int last_sample_number, int samples_per_channel )
{
  m_displayedXRangeChange.emit( first_sample_number, last_sample_number, samples_per_channel );
}//displayedXRangeChangeCallback(...)
