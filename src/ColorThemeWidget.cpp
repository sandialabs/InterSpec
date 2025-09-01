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

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WComboBox>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WSuggestionPopup>
#include <Wt/WContainerWidget>

#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/ColorThemeWidget.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

ColorThemeWidget::ColorThemeWidget(WContainerWidget *parent)
	: WContainerWidget(parent),
	m_themeTitle( nullptr ),
	m_themeDescription( nullptr ),
  m_nonChartAreaCssTheme( nullptr ),
	m_colorSelects{ nullptr },
	m_specMarginSameAsBackground( nullptr ),
	m_timeMarginSameAsBackground( nullptr ),
  m_noSpectrumBackground( nullptr ),
  m_noSpectrumMargin( nullptr ),
  m_noTimeBackground( nullptr ),
  m_noTimeMargin( nullptr ),
	m_peaksTakeRefLineColor( nullptr ),
  m_peakLabelFontSize( nullptr ),
  m_peakLabelAngle( nullptr ),
  m_logYAxisMin( nullptr ),
  m_referenceLineColor{ nullptr },
  m_specificRefLineName{ nullptr },
  m_specificRefLineColor{ nullptr }
{
  wApp->useStyleSheet( "InterSpec_resources/ColorThemeWidget.css" );
  
  addStyleClass( "ColorThemeWidget" );
  
	WTable *table = new WTable( this );
  table->addStyleClass( "ColorThemeSelectTable" );

  const bool nativeColorSelect = ColorSelect::willUseNativeColorPicker();
  
	int row = 0;
	WTableCell *cell = table->elementAt(row, 0);
	WLabel *label = new WLabel("Theme Name", cell);
	cell = table->elementAt(row, 1);
	cell->setColumnSpan(2);
	m_themeTitle = new WLineEdit("Title not assigned", cell);
  m_themeTitle->setAttributeValue( "ondragstart", "return false" );
	label->setBuddy(m_themeTitle);
	m_themeTitle->setWidth(WLength(95, WLength::Percentage));

	++row;
	cell = table->elementAt(row, 0);
	label = new WLabel("Description",cell);
	cell = table->elementAt(row, 1);
	cell->setColumnSpan(2);
  m_themeDescription = new WLineEdit("Desc not assigned", cell);
  m_themeDescription->setAttributeValue( "ondragstart", "return false" );
	label->setBuddy(m_themeDescription);
	m_themeDescription->setWidth(WLength(95, WLength::Percentage));
	
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->setColumnSpan(3);
  cell->setContentAlignment(Wt::AlignmentFlag::AlignCenter);
  cell->addStyleClass( "ThemeColorsTitle" );
  new WText("Non-Chart Area", cell);
  
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Backdrop Theme", cell);
  cell = table->elementAt(row, 1);
  
  //ToDo: right m_nonChartAreaCssTheme only lists names given by ColorTheme::predefinedThemeName(),
  //      but what we should *really* do is list options in InterSpec_resources/themes
  m_nonChartAreaCssTheme = new WComboBox( cell );
  m_nonChartAreaCssTheme->setNoSelectionEnabled( true );
  
  for( ColorTheme::PredefinedColorTheme t = ColorTheme::DefaultColorTheme;
      t < ColorTheme::NumPredefinedColorTheme;
      t = ColorTheme::PredefinedColorTheme(t+1) )
  {
    m_nonChartAreaCssTheme->addItem( ColorTheme::predefinedThemeName(t) );
  }
  
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("Theme for all non-chart area", cell);
  

	++row;
	cell = table->elementAt(row, 0);
	cell->setColumnSpan(3);
	cell->setContentAlignment(Wt::AlignmentFlag::AlignCenter);
  cell->addStyleClass( "ThemeColorsTitle" );
	new WText("Chart Colors", cell);
  

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spec. Foreground", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[ForegroundLine] = new ColorSelect(0,cell); 
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Line color for the foreground spectrum", cell);
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spec. Background", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[BackgroundLine] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Line color for the background spectrum", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spec. Secondary", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SecondaryLine] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Line color for the secondary spectrum", cell);
	


	++row;
	cell = table->elementAt(row, 0);
	cell->setColumnSpan(2);
	m_peaksTakeRefLineColor = new WCheckBox("Peaks Use Reference Line Color", cell);
	cell = table->elementAt(row, 2);
	new WText("Determines if newly fitted peaks take on color of reference line for that peak.", cell);


	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Default Peak Color", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DefaultPeakLine] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color used for peaks that dont have an assigned color", cell);
	

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line Medical", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineMedical] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from medical sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line Industrial", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineIndustrial] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from industrial sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line NORM", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineNorm] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from NORM sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line SNM", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineSnm] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from SNM sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line Common", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineCommon] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from common sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Dynamic Ref Line Other", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DynamicRefLineOther] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for dynamic reference lines from other/uncategorized sources", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spectrum Chart Axis", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumAxis] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Line color for spectrum X and Y axis", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spectrum Chart Text", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumChartText] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Spectrum chart text (y and x axis labels) color", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spectrum Area", cell);
	cell = table->elementAt(row, 1);
  
  if( nativeColorSelect )
  {
    m_noSpectrumBackground = new WCheckBox( "None", cell );
    m_noSpectrumBackground->setInline( false );
  }
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumChartBackground] = new ColorSelect(ColorSelect::AllowNoColor,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Background color for spectrum chart", cell);
    
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Spectrum Chart Bckgrnd", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noSpectrumMargin = new WCheckBox( "None", cell );
    m_noSpectrumMargin->setInline( false );
  }
	m_colorSelects[SpectrumChartMargins] = new ColorSelect(ColorSelect::AllowNoColor,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for spectrum chart margins", cell);
	m_specMarginSameAsBackground = new WCheckBox("Same as spectrum area", cell);
	m_specMarginSameAsBackground->setInline(false);


  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Time Chart Gamma Counts", cell);
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[TimeChartGammas] = new ColorSelect(0,cell);
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("Line color for the gamma counts in time chart", cell);
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Time Chart Neutron Counts", cell);
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[TimeChartNeutrons] = new ColorSelect(0,cell);
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("Line color for the neutron counts in time chart", cell);
  
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time Chart Axis", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeAxisLines] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Line color for time chart X and Y axis", cell);

  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time Chart Background", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noTimeBackground = new WCheckBox( "None", cell );
    m_noTimeBackground->setInline( false );
  }
  m_colorSelects[TimeChartBackground] = new ColorSelect(ColorSelect::AllowNoColor,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Background color for time chart", cell);
  

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time Chart Margins", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noTimeMargin = new WCheckBox( "None", cell );
    m_noTimeMargin->setInline( false );
  }
  m_colorSelects[TimeChartMargins] = new ColorSelect(ColorSelect::AllowNoColor,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color for spectrum chart margins", cell);
	m_timeMarginSameAsBackground = new WCheckBox("Same as background", cell);
	m_timeMarginSameAsBackground->setInline(false);
  
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time Chart Text", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeChartText] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Text color for time chart (x and y-axis labels)", cell);

  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Occ. Indicator Lines", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[OccupancyIndicatorLines] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Lines that indicate the beginning and end of an occupancy (portal data) or time period of interest (search mode data)", cell);

	
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time foreground highlight", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistoryForegroundHighlight] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color used to indicate on the time chart, which segments of time are being summed for the foreground spectrum.", cell);

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time background highlight", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistoryBackgroundHighlight] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color used to indicate on the time chart, which segments of time are being summed for the background spectrum.", cell);
		
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	new WLabel("Time secondary highlight", cell);
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistorySecondaryHighlight] = new ColorSelect(0,cell);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	new WText("Color used to indicate on the time chart, which segments of time are being summed for the secondary spectrum.", cell);
	
  
  WContainerWidget *refLineContainer = new WContainerWidget( this );
  refLineContainer->addStyleClass( "RefLinesArea" );
  WGridLayout *refLayout = new WGridLayout();
  refLayout->setContentsMargins( 0, 0, 0, 0 );
  refLayout->setVerticalSpacing( 0 );
  refLayout->setHorizontalSpacing( 0 );
  refLineContainer->setLayout(refLayout);
  
  WContainerWidget *genericRef = new WContainerWidget();
  genericRef->addStyleClass( "GenericRefLineColors" );
  refLayout->addWidget( genericRef, 0, 0 );
  
  WText *genericTitle = new WText( "Reference Line Colors", genericRef );
  genericTitle->addStyleClass( "ColorRefLineTitle" );
  genericTitle->setInline( false );
  
  WText *genericDesc = new WText( "The ordered colors to use for reference lines", genericRef );
  genericDesc->addStyleClass( "ColorRefLineDesc" );
  genericDesc->setInline( false );
  
  table = new WTable( genericRef );
  table->addStyleClass( "ReferenceLinColorTable" );
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    cell = table->elementAt(i,0);
    new WText( std::to_string(i+1), cell );
    cell->addStyleClass( "GenericRefLineNumber" );
    cell = table->elementAt(i,1);
    m_referenceLineColor[i] = new ColorSelect( 0, cell );
    m_referenceLineColor[i]->cssColorChanged().connect( boost::bind( &ColorThemeWidget::genericRefLineColorChangedCallback, this, i ) );
  }
  
  WContainerWidget *specificRef = new WContainerWidget();
  specificRef->addStyleClass( "SpecificRefLineColors" );
  refLayout->addWidget( specificRef, 0, 1 );
  
  WText *specificTitle = new WText( "Specific Source Reference Line Colors", specificRef );
  specificTitle->addStyleClass( "ColorRefLineTitle" );
  specificTitle->setInline( false );
  
  WText *specificDesc = new WText( "The colors to use for reference lines for given sources", specificRef );
  specificDesc->addStyleClass( "ColorRefLineDesc" );
  specificDesc->setInline( false );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->addCustomSuggestPossibility( "background" );
  WSuggestionPopup *nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  nuclideSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );

  //See ReferencePhotopeakDisplay.cpp for note on this next hack section of code
  string js = INLINE_JAVASCRIPT(
    var addTryCatch = function( elid ){
      var dofix = function(elid){
        var el = Wt.WT.getElement(elid);
        var self = el ? jQuery.data(el, 'obj') : null;
        if( !self ) self = el ? el.wtObj : null;; //Wt 3.7.1
        if( !self ){ setTimeout( function(){dofix(elid);}, 100 ); return; }
        var oldfcn = self.refilter;
        self.refilter = function(value){ try{ oldfcn(value); }catch(e){ console.log('My refilter caught: ' + e ); } };
      };
      dofix(elid);
    };
  );
  
  nuclideSuggest->doJavaScript( js + " addTryCatch('" + nuclideSuggest->id() + "');" );
  
  isoSuggestModel->filter( "" );
  nuclideSuggest->setFilterLength( -1 );
  nuclideSuggest->setModel( isoSuggestModel );
  nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  
  
  table = new WTable( specificRef );
  table->addStyleClass( "SpecificRefLineColorTable" );
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    cell = table->elementAt(i,0);
    m_specificRefLineName[i] = new WLineEdit( cell );
    m_specificRefLineName[i]->setAttributeValue( "ondragstart", "return false" );
    nuclideSuggest->forEdit( m_specificRefLineName[i], WSuggestionPopup::Editing );
    m_specificRefLineName[i]->changed().connect( boost::bind( &ColorThemeWidget::specificRefLineSourceChangedCallback, this, i ) );
    m_specificRefLineName[i]->enterPressed().connect( boost::bind( &ColorThemeWidget::specificRefLineSourceChangedCallback, this, i ) );
    
    m_specificRefLineName[i]->addStyleClass( "SpecificRefLineInput" );
    cell = table->elementAt(i,1);
    cell->addStyleClass( "SpecificRefLineColorCell" );
    m_specificRefLineColor[i] = new ColorSelect( ColorSelect::AllowNoColor, cell );
    m_specificRefLineColor[i]->cssColorChanged().connect( boost::bind( &ColorThemeWidget::specificRefLineColorChangedCallback, this, i ) );
    m_specificRefLineColor[i]->setColor( WColor() );
  }
  
  
  
  
  row = 0;
  table = new WTable( this );
  table->addStyleClass( "ColorThemeOtherOptionsTable" );
  
  cell = table->elementAt(row, 0);
  cell->setColumnSpan(3);
  cell->setContentAlignment(Wt::AlignmentFlag::AlignCenter);
  cell->addStyleClass( "ThemePeakLabelTitle" );
  new WText("Other Options", cell);
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Peak label size", cell);
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_peakLabelFontSize = new WComboBox( cell);
  m_peakLabelFontSize->addItem( "Default" );
  m_peakLabelFontSize->addItem( "smaller" );
  m_peakLabelFontSize->addItem( "larger" );
  m_peakLabelFontSize->addItem( "xx-small" );
  m_peakLabelFontSize->addItem( "x-small" );
  m_peakLabelFontSize->addItem( "small" );
  m_peakLabelFontSize->addItem( "medium" );
  m_peakLabelFontSize->addItem( "large" );
  m_peakLabelFontSize->addItem( "x-large" );
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("Specifies font-size of peak labels", cell);
  
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Peak label angle", cell);
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_peakLabelAngle = new NativeFloatSpinBox( cell);
  m_peakLabelAngle->setSpinnerHidden( true );
  m_peakLabelAngle->setRange( -180, 270 );
  m_peakLabelAngle->setWidth( 50 );
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("Rotation, in degrees, of peak labels", cell);
  
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  new WLabel("Log-Y minimum", cell);
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_logYAxisMin = new NativeFloatSpinBox( cell);
  m_logYAxisMin->setSpinnerHidden( true );
  m_logYAxisMin->setRange( 1.0E-8, 1000 );
  m_logYAxisMin->setWidth( 50 );
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  new WText("The minimum value displayed on the logarithmic y-axis, when there is a channel"
            " with zero or negative counts in the displayed x-range.", cell);
  
  
  m_themeTitle->changed().connect( boost::bind( &ColorThemeWidget::titleChangedCallback, this ) );
  m_themeDescription->changed().connect( boost::bind( &ColorThemeWidget::descriptionChangedCallback, this ) );
  m_peaksTakeRefLineColor->checked().connect( boost::bind( &ColorThemeWidget::peaksTakeRefLineColorChangedCallback, this ) );
  m_peaksTakeRefLineColor->unChecked().connect( boost::bind( &ColorThemeWidget::peaksTakeRefLineColorChangedCallback, this ) );

  m_peakLabelFontSize->activated().connect( boost::bind( &ColorThemeWidget::peakLabelFontSizeChanged, this ) );
  m_peakLabelAngle->valueChanged().connect( boost::bind( &ColorThemeWidget::peakLabelRotationChanged, this ) );
  m_logYAxisMin->valueChanged().connect( boost::bind( &ColorThemeWidget::logYAxisMinValueChanged, this ) );
  
  m_nonChartAreaCssTheme->changed().connect( boost::bind( &ColorThemeWidget::nonChartAreaThemeChanged, this ) );
  
  m_specMarginSameAsBackground->checked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, SpectrumChartMargins ) );
  m_specMarginSameAsBackground->unChecked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, SpectrumChartMargins ) );
  
  m_timeMarginSameAsBackground->checked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartMargins ) );
  m_timeMarginSameAsBackground->unChecked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartMargins ) );

  if( m_noSpectrumBackground )
  {
    m_noSpectrumBackground->checked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, SpectrumChartBackground ) );
    m_noSpectrumBackground->unChecked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, SpectrumChartBackground ) );
  }
  
  if( m_noSpectrumMargin )
    m_noSpectrumMargin->changed().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartMargins ) );
  
  if( m_noTimeBackground )
  {
    m_noTimeBackground->checked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartBackground ) );
    m_noTimeBackground->unChecked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartBackground ) );
  }
  
  if( m_noTimeMargin )
  {
    m_noTimeMargin->checked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartMargins ) );
    m_noTimeMargin->unChecked().connect( boost::bind( &ColorThemeWidget::newColorSelectedCallback, this, TimeChartMargins ) );
  }


	for( SelectableColor color = SelectableColor(0); color < NumSelectableColors; color = SelectableColor(color + 1) )
	{
		m_colorSelects[color]->cssColorChanged().connect(boost::bind(&ColorThemeWidget::newColorSelectedCallback, this, color));
	}
}//ColorThemeWidget( constructor )


ColorThemeWidget::~ColorThemeWidget()
{

}//ColorThemeWidget( destructor )


void ColorThemeWidget::setTheme(const ColorTheme *theme, const bool modifieable)
{
	m_currentTheme.reset();

	if( !theme )
	{
		m_themeTitle->setText("(empty)");
		m_themeDescription->setText("");
		m_origTheme.reset();
    setDisabled( true );
		return;
	}
  
  const bool wasEnabled = isEnabled();
  if( wasEnabled != modifieable )
  {
    setDisabled( !modifieable );
    
    //We have to explicitly enable the spectrum.js based color pickers if they
    //  were created with the parent disabled (due to using a default color
    //  theme)
    for( auto i : m_colorSelects )
      i->setDisabled( !modifieable );
    for( auto i : m_referenceLineColor )
      i->setDisabled( !modifieable );
    for( auto i : m_specificRefLineColor )
      i->setDisabled( !modifieable );
  }//if( wasEnabled != modifieable )
  
  m_themeModifiable = modifieable;
  
	m_origTheme.reset( new ColorTheme(*theme) );
	m_themeTitle->setText(theme->theme_name);
	m_themeDescription->setText(theme->theme_description);

  
  m_peakLabelFontSize->setCurrentIndex(0);
  for( int i = 1; i < m_peakLabelFontSize->count(); ++i )
  {
    if( m_peakLabelFontSize->itemText(i).toUTF8() == theme->spectrumPeakLabelSize )
    {
      m_peakLabelFontSize->setCurrentIndex(i);
      break;
    }
  }//
  
  m_peakLabelAngle->setValue( -static_cast<float>(theme->spectrumPeakLabelRotation) );
  if( theme->spectrumLogYAxisMin > 0.0 )
    m_logYAxisMin->setValue( static_cast<float>(theme->spectrumLogYAxisMin) );
  else
    m_logYAxisMin->setValue( 0.1f );
  
  if( theme->nonChartAreaTheme.empty() || theme->nonChartAreaTheme=="default" )
  {
    m_nonChartAreaCssTheme->setCurrentIndex( 0 );
  }else
  {
    m_nonChartAreaCssTheme->setCurrentIndex( -1 );
    for( ColorTheme::PredefinedColorTheme t = ColorTheme::DefaultColorTheme;
        t < ColorTheme::NumPredefinedColorTheme;
        t = ColorTheme::PredefinedColorTheme(t+1) )
    {
      //ToDo: right m_nonChartAreaCssTheme only lists names given by ColorTheme::predefinedThemeName(),
      //      but what we should *really* do is list options in InterSpec_resources/themes
      if( theme->nonChartAreaTheme == ColorTheme::predefinedThemeName(t) )
      {
        m_nonChartAreaCssTheme->setCurrentIndex( t - ColorTheme::DefaultColorTheme );
        break;
      }
    }
  }//if( default ) / else
  
	m_peaksTakeRefLineColor->setChecked(theme->peaksTakeOnReferenceLineColor);

	m_colorSelects[ForegroundLine]->setColor(theme->foregroundLine);
	m_colorSelects[BackgroundLine]->setColor(theme->backgroundLine);
	m_colorSelects[SecondaryLine]->setColor(theme->secondaryLine);
	m_colorSelects[DefaultPeakLine]->setColor(theme->defaultPeakLine);
	m_colorSelects[SpectrumAxis]->setColor(theme->spectrumAxisLines);
	m_colorSelects[SpectrumChartText]->setColor(theme->spectrumChartText);
	m_colorSelects[SpectrumChartBackground]->setColor(theme->spectrumChartBackground);
  const bool noSpecBack = theme->spectrumChartBackground.isDefault();
  if( m_noSpectrumBackground )
  {
    m_noSpectrumBackground->setChecked( noSpecBack );
    m_colorSelects[SpectrumChartBackground]->setHidden( noSpecBack );
  }

  //No way this logic to control checking/unchecking, and showing.hiding the checkboxes is correct in the first go.
  const bool specMarginSame = theme->spectrumChartMargins.isDefault();
  m_colorSelects[SpectrumChartMargins]->setColor( theme->spectrumChartMargins );
  m_specMarginSameAsBackground->setChecked( specMarginSame );
  if( m_noSpectrumMargin )
  {
    m_colorSelects[SpectrumChartMargins]->setHidden( (specMarginSame && !noSpecBack) );
    m_noSpectrumMargin->setChecked( specMarginSame );
    m_noSpectrumMargin->setHidden( specMarginSame );
  }

  
  m_colorSelects[TimeChartGammas]->setColor(theme->timeChartGammaLine);
  m_colorSelects[TimeChartNeutrons]->setColor(theme->timeChartNeutronLine);
  
  
	m_colorSelects[TimeAxisLines]->setColor(theme->timeAxisLines);
	m_colorSelects[TimeChartText]->setColor(theme->timeChartText);


  const bool noTimeBack = theme->timeChartBackground.isDefault();
  if( m_noTimeBackground )
  {
    m_noTimeBackground->setChecked( noTimeBack );
    m_colorSelects[TimeChartBackground]->setHidden( noTimeBack );
  }
  m_colorSelects[TimeChartBackground]->setColor( theme->timeChartBackground );
  

  //No way this logic to control checking/unchecking, and showing.hiding the checkboxes is correct in the first go.
  const bool timeMarginSame = theme->timeChartMargins.isDefault();
  m_colorSelects[TimeChartMargins]->setColor( theme->timeChartMargins );
  m_timeMarginSameAsBackground->setChecked( timeMarginSame );
  m_timeMarginSameAsBackground->setHidden( !timeMarginSame );
  if( m_noTimeMargin )
  {
    m_colorSelects[TimeChartMargins]->setHidden( (timeMarginSame && !noTimeBack) );
    m_noTimeMargin->setChecked( timeMarginSame );
    m_noTimeMargin->setHidden( timeMarginSame );
  }


		
	m_colorSelects[OccupancyIndicatorLines]->setColor(theme->occupancyIndicatorLines);
	m_colorSelects[TimeHistoryForegroundHighlight]->setColor(theme->timeHistoryForegroundHighlight);
	m_colorSelects[TimeHistoryBackgroundHighlight]->setColor(theme->timeHistoryBackgroundHighlight);
	m_colorSelects[TimeHistorySecondaryHighlight]->setColor(theme->timeHistorySecondaryHighlight);
	m_colorSelects[DynamicRefLineMedical]->setColor( theme->dynamicRefLineMedicalColor.isDefault() ? 
                                                    Wt::WColor( ColorTheme::sm_dynamic_ref_line_medical_color ) : 
                                                    theme->dynamicRefLineMedicalColor );
	m_colorSelects[DynamicRefLineIndustrial]->setColor( theme->dynamicRefLineIndustrialColor.isDefault() ? 
                                                       Wt::WColor( ColorTheme::sm_dynamic_ref_line_industrial_color ) : 
                                                       theme->dynamicRefLineIndustrialColor );
	m_colorSelects[DynamicRefLineNorm]->setColor( theme->dynamicRefLineNormColor.isDefault() ? 
                                                 Wt::WColor( ColorTheme::sm_dynamic_ref_line_norm_color ) : 
                                                 theme->dynamicRefLineNormColor );
	m_colorSelects[DynamicRefLineSnm]->setColor( theme->dynamicRefLineSnmColor.isDefault() ? 
                                                Wt::WColor( ColorTheme::sm_dynamic_ref_line_snm_color ) : 
                                                theme->dynamicRefLineSnmColor );
	m_colorSelects[DynamicRefLineCommon]->setColor( theme->dynamicRefLineCommonColor.isDefault() ? 
                                                   Wt::WColor( ColorTheme::sm_dynamic_ref_line_common_color ) : 
                                                   theme->dynamicRefLineCommonColor );
	m_colorSelects[DynamicRefLineOther]->setColor( theme->dynamicRefLineOtherColor.isDefault() ? 
                                                  Wt::WColor( ColorTheme::sm_dynamic_ref_line_other_color ) : 
                                                  theme->dynamicRefLineOtherColor );
  
  
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    m_referenceLineColor[i]->setColor( WColor() );
    m_specificRefLineName[i]->setText( "" );
    m_specificRefLineColor[i]->setColor( WColor() );
  }
  
  for( size_t i = 0; i < theme->referenceLineColor.size(); ++i )
  {
    if( static_cast<int>(i) < sm_numRefLineColors )
      m_referenceLineColor[i]->setColor( theme->referenceLineColor[i] );
  }
  

  int position = 0;
  for( const auto &c : theme->referenceLineColorForSources )
  {
    if( !c.first.empty() && !c.second.isDefault() && position < sm_numRefLineColors )
    {
      m_specificRefLineName[position]->setText( c.first );
      m_specificRefLineColor[position]->show();
      m_specificRefLineColor[position]->setColor( c.second );
      ++position;
    }else
    {
      m_specificRefLineColor[position]->hide();
    }
  }
  
}//setTheme(...)


Wt::Signal<> &ColorThemeWidget::edited()
{
  return m_edited;
}


void ColorThemeWidget::newColorSelectedCallback(const ColorThemeWidget::SelectableColor color)
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
	switch( color )
	{
	  case ForegroundLine:
      m_currentTheme->foregroundLine = m_colorSelects[color]->color();
      break;
      
	  case BackgroundLine:
      m_currentTheme->backgroundLine = m_colorSelects[color]->color();
      break;
      
    case SecondaryLine:
      m_currentTheme->secondaryLine = m_colorSelects[color]->color();
      break;
      
	  case DefaultPeakLine:
      m_currentTheme->defaultPeakLine = m_colorSelects[color]->color();
      break;
      
	  case SpectrumAxis:
      m_currentTheme->spectrumAxisLines = m_colorSelects[color]->color();
      break;
      
	  case SpectrumChartText:
      m_currentTheme->spectrumChartText = m_colorSelects[color]->color();
      break;
      
	  case SpectrumChartBackground:
      if( m_noSpectrumBackground )
      {
        m_colorSelects[color]->setHidden( m_noSpectrumBackground->isChecked() );
        if( m_noSpectrumBackground->isChecked() )
          m_currentTheme->spectrumChartBackground = WColor();
        else
          m_currentTheme->spectrumChartBackground = m_colorSelects[color]->color();
      }else
      {
        m_currentTheme->spectrumChartBackground = m_colorSelects[color]->color();
      }
      break;
      
	  case SpectrumChartMargins:
      if( m_noSpectrumMargin )
      {
        m_colorSelects[color]->setHidden( m_specMarginSameAsBackground->isChecked() || m_noSpectrumMargin->isChecked() );
        m_noSpectrumMargin->setHidden( m_specMarginSameAsBackground->isChecked() );
        if( m_specMarginSameAsBackground->isChecked() )
          m_noSpectrumMargin->setChecked( false );
      }
      
      if( m_noSpectrumMargin )
      {
        if( m_specMarginSameAsBackground->isChecked() || m_noSpectrumMargin->isChecked() )
          m_currentTheme->spectrumChartMargins = WColor();
        else
          m_currentTheme->spectrumChartMargins = m_colorSelects[color]->color();
      }else
      {
        if( m_specMarginSameAsBackground->isChecked() )
          m_currentTheme->spectrumChartMargins = WColor();
        else
          m_currentTheme->spectrumChartMargins = m_colorSelects[color]->color();
      }
      break;
      
    case TimeChartGammas:
      m_currentTheme->timeChartGammaLine = m_colorSelects[color]->color();
      break;
      
    case TimeChartNeutrons:
      m_currentTheme->timeChartNeutronLine = m_colorSelects[color]->color();
      break;
      
	  case TimeAxisLines:
      m_currentTheme->timeAxisLines = m_colorSelects[color]->color();
      break;
      
	  case TimeChartText:
      m_currentTheme->timeChartText = m_colorSelects[color]->color();
      break;
      
	  case TimeChartBackground:
      if( m_noTimeBackground )
      {
        m_colorSelects[color]->setHidden( m_noTimeBackground->isChecked() );
        if( m_noTimeBackground->isChecked() )
          m_currentTheme->timeChartBackground = WColor();
        else
          m_currentTheme->timeChartBackground = m_colorSelects[color]->color();
      }else
      {
        m_currentTheme->timeChartBackground = m_colorSelects[color]->color();
      }
      
      break;
      
	  case TimeChartMargins:
      if( m_noTimeMargin )
      {
        m_colorSelects[color]->setHidden( m_noTimeMargin->isChecked() || m_timeMarginSameAsBackground->isChecked() );
      
        m_noTimeMargin->setHidden( m_timeMarginSameAsBackground->isChecked() );
        if( m_timeMarginSameAsBackground->isChecked() )
          m_noTimeMargin->setChecked( false );
      
        if( m_timeMarginSameAsBackground->isChecked() || m_noTimeMargin->isChecked() )
          m_currentTheme->timeChartBackground = WColor();
        else
          m_currentTheme->timeChartBackground = m_colorSelects[color]->color();
      }else
      {
        if( m_timeMarginSameAsBackground->isChecked() )
          m_currentTheme->timeChartBackground = WColor();
        else
          m_currentTheme->timeChartBackground = m_colorSelects[color]->color();
      }
      break;
      
	  case OccupancyIndicatorLines:
      m_currentTheme->occupancyIndicatorLines = m_colorSelects[color]->color();
      break;
      
	  case TimeHistoryForegroundHighlight:
      m_currentTheme->timeHistoryForegroundHighlight = m_colorSelects[color]->color();
      break;
      
	  case TimeHistoryBackgroundHighlight:
      m_currentTheme->timeHistoryBackgroundHighlight = m_colorSelects[color]->color();
      break;
      
	  case TimeHistorySecondaryHighlight:
      m_currentTheme->timeHistorySecondaryHighlight = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineMedical:
      m_currentTheme->dynamicRefLineMedicalColor = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineIndustrial:
      m_currentTheme->dynamicRefLineIndustrialColor = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineNorm:
      m_currentTheme->dynamicRefLineNormColor = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineSnm:
      m_currentTheme->dynamicRefLineSnmColor = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineCommon:
      m_currentTheme->dynamicRefLineCommonColor = m_colorSelects[color]->color();
      break;
      
	  case DynamicRefLineOther:
      m_currentTheme->dynamicRefLineOtherColor = m_colorSelects[color]->color();
      break;
      
	  case NumSelectableColors:
		  break;
	}//switch (color)

  m_edited.emit();
}//void newColorSelectedCallback(color)


void ColorThemeWidget::titleChangedCallback()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    if( m_origTheme->theme_name == m_themeTitle->text() )
      return;
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  m_currentTheme->theme_name = m_themeTitle->text();
  m_edited.emit();
}//void titleChangedCallback()


void ColorThemeWidget::descriptionChangedCallback()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    if( m_origTheme->theme_description == m_themeDescription->text() )
      return;
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  m_currentTheme->theme_description = m_themeDescription->text();
  m_edited.emit();
}//void descriptionChangedCallback()


void ColorThemeWidget::nonChartAreaThemeChanged()
{
  const ColorTheme::PredefinedColorTheme theme = ColorTheme::PredefinedColorTheme(ColorTheme::DefaultColorTheme + m_nonChartAreaCssTheme->currentIndex());
  string themeName = ColorTheme::predefinedThemeName(theme);
  if( theme == ColorTheme::DefaultColorTheme )
    themeName = "";
  
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    if( m_origTheme->nonChartAreaTheme == themeName )
      return;
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  m_currentTheme->nonChartAreaTheme = themeName;
  
  m_edited.emit();
}//void nonChartAreaThemeChanged()


void ColorThemeWidget::peaksTakeRefLineColorChangedCallback()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    if( m_origTheme->peaksTakeOnReferenceLineColor == m_peaksTakeRefLineColor->isChecked() )
      return;
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  m_currentTheme->peaksTakeOnReferenceLineColor = m_peaksTakeRefLineColor->isChecked();
  m_edited.emit();
}//void peaksTakeRefLineColorChangedCallback()


void ColorThemeWidget::peakLabelFontSizeChanged()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  if( m_peakLabelFontSize->currentIndex() == 0 )
    m_currentTheme->spectrumPeakLabelSize = "";
  else
    m_currentTheme->spectrumPeakLabelSize = m_peakLabelFontSize->currentText();
  
  m_edited.emit();
}//void peakLabelFontSizeChanged()


void ColorThemeWidget::peakLabelRotationChanged()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  m_currentTheme->spectrumPeakLabelRotation = -m_peakLabelAngle->value();
  m_edited.emit();
}//void peakLabelRotationChanged()


void ColorThemeWidget::logYAxisMinValueChanged()
{
  if( !m_currentTheme )
  {
    assert( m_origTheme );
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  }
  
  float val = m_logYAxisMin->value();
  if( val <= 0.0f )
  {
    val = 0.1f;
    m_logYAxisMin->setValue( val );
  }
  
  m_currentTheme->spectrumLogYAxisMin = val;
  m_edited.emit();
}//void logYAxisMinValueChanged()


void ColorThemeWidget::genericRefLineColorChangedCallback( const int num )
{
  if( num < 0 || num >= sm_numRefLineColors ) //shouldnt ever happen
    return;
 
  //ToDo: This funtion is only minimially hacked - make better
  if( !m_currentTheme )
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  
  m_currentTheme->referenceLineColor.clear();
  for( int i = 0; i < sm_numRefLineColors; ++i )
    m_currentTheme->referenceLineColor.push_back( m_referenceLineColor[i]->color() );
  
  m_edited.emit();
}//void genericRefLineColorChangedCallback( const int num );


void ColorThemeWidget::specificRefLineColorChangedCallback( const int num )
{
  if( num < 0 || num >= sm_numRefLineColors ) //shouldnt ever happen
    return;
  
  //ToDo: This funtion is only minimially hacked - make better
  if( !m_currentTheme )
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  
  m_currentTheme->referenceLineColorForSources.clear();
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    const string symbol = m_specificRefLineName[i]->text().toUTF8();
    if( !symbol.empty() )
    {
      m_currentTheme->referenceLineColorForSources[symbol] = m_specificRefLineColor[i]->color();
    }else
    {
      m_specificRefLineColor[i]->setColor( WColor() );
    }
  }
  
  m_edited.emit();
}//void specificRefLineColorChangedCallback( const int num );


void ColorThemeWidget::specificRefLineSourceChangedCallback( const int num )
{
  if( num < 0 || num >= sm_numRefLineColors ) //shouldnt ever happen
    return;
  
  //ToDo: This funtion is only minimially hacked - make better
  if( !m_currentTheme )
    m_currentTheme.reset( new ColorTheme(*m_origTheme) );
  
  m_currentTheme->referenceLineColorForSources.clear();
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    const string symbol = m_specificRefLineName[i]->text().toUTF8();
    if( !symbol.empty() )
    {
      m_currentTheme->referenceLineColorForSources[symbol] = m_specificRefLineColor[i]->color();
      if( m_currentTheme->referenceLineColorForSources[symbol].isDefault() )
      {
        m_specificRefLineColor[i]->setColor( WColor("#000000") );
        m_currentTheme->referenceLineColorForSources[symbol] = WColor("#000000");
      }
    }else
    {
      m_specificRefLineColor[i]->setColor( WColor() );
    }
  }
  
  m_edited.emit();
}//void specificRefLineSourceChangedCallback( const int num );


