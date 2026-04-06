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

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WComboBox.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WSuggestionPopup.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/ColorThemeWidget.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/IsotopeNameFilterModel.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__

ColorThemeWidget::ColorThemeWidget()
	: WContainerWidget(),
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
  InterSpec::instance()->useMessageResourceBundle( "ColorThemeWidget" );
  wApp->useStyleSheet( "InterSpec_resources/ColorThemeWidget.css" );
  
  addStyleClass( "ColorThemeWidget" );
  
	WTable *table = addNew<WTable>();
  table->addStyleClass( "ColorThemeSelectTable" );

  const bool nativeColorSelect = ColorSelect::willUseNativeColorPicker();
  
	int row = 0;
	WTableCell *cell = table->elementAt(row, 0);
	WLabel *label = cell->addNew<WLabel>(WString::tr("ctwidget-theme-name"));
	cell = table->elementAt(row, 1);
	cell->setColumnSpan(2);
	m_themeTitle = cell->addNew<WLineEdit>(WString::tr("ctwidget-title-placeholder"));
  m_themeTitle->setAttributeValue( "ondragstart", "return false" );
	label->setBuddy(m_themeTitle);
	m_themeTitle->setWidth(WLength(95, WLength::Unit::Percentage));

	++row;
	cell = table->elementAt(row, 0);
	label = cell->addNew<WLabel>(WString::tr("ctwidget-description"));
	cell = table->elementAt(row, 1);
	cell->setColumnSpan(2);
  m_themeDescription = cell->addNew<WLineEdit>(WString::tr("ctwidget-desc-placeholder"));
  m_themeDescription->setAttributeValue( "ondragstart", "return false" );
	label->setBuddy(m_themeDescription);
	m_themeDescription->setWidth(WLength(95, WLength::Unit::Percentage));
	
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->setColumnSpan(3);
  cell->setContentAlignment(Wt::AlignmentFlag::Center);
  cell->addStyleClass( "ThemeColorsTitle" );
  cell->addNew<WText>(WString::tr("ctwidget-non-chart-area"));
  
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-backdrop-theme"));
  cell = table->elementAt(row, 1);

  //ToDo: right m_nonChartAreaCssTheme only lists names given by ColorTheme::predefinedThemeName(),
  //      but what we should *really* do is list options in InterSpec_resources/themes
  m_nonChartAreaCssTheme = cell->addNew<WComboBox>();
  m_nonChartAreaCssTheme->setNoSelectionEnabled( true );
  
  for( ColorTheme::PredefinedColorTheme t = ColorTheme::DefaultColorTheme;
      t < ColorTheme::NumPredefinedColorTheme;
      t = ColorTheme::PredefinedColorTheme(t+1) )
  {
    m_nonChartAreaCssTheme->addItem( ColorTheme::predefinedThemeName(t) );
  }
  
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-backdrop-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-background"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppBackground] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-background-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-text"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppText] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-text-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-border"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppBorder] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-border-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-link"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppLink] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-link-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-label"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppLabel] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-label-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-input-bg"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppInputBackground] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-input-bg-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-button-bg"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppButtonBackground] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-button-bg-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-button-border"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppButtonBorder] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-button-border-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-button-text"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppButtonText] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-button-text-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-menubar-bg"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppMenuBarBackground] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-menubar-bg-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-menubar-active"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppMenuBarActiveColor] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-menubar-active-desc"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-app-menubar-hover"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[AppMenuBarHoverColor] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-app-menubar-hover-desc"));

	++row;
	cell = table->elementAt(row, 0);
	cell->setColumnSpan(3);
	cell->setContentAlignment(Wt::AlignmentFlag::Center);
  cell->addStyleClass( "ThemeColorsTitle" );
	cell->addNew<WText>(WString::tr("ctwidget-chart-colors"));
  

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spec-foreground"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[ForegroundLine] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spec-foreground-desc"));
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spec-background"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[BackgroundLine] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spec-background-desc"));

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spec-secondary"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SecondaryLine] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spec-secondary-desc"));
	


	++row;
	cell = table->elementAt(row, 0);
	cell->setColumnSpan(2);
	m_peaksTakeRefLineColor = cell->addNew<WCheckBox>(WString::tr("ctwidget-peaks-use-ref-color"));
	cell = table->elementAt(row, 2);
	cell->addNew<WText>(WString::tr("ctwidget-peaks-use-ref-desc"));


	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-default-peak-color"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[DefaultPeakLine] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-default-peak-desc"));

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spectrum-axis"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumAxis] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spectrum-axis-desc"));

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spectrum-text"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumChartText] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spectrum-text-desc"));

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spectrum-area"));
	cell = table->elementAt(row, 1);

  if( nativeColorSelect )
  {
    m_noSpectrumBackground = cell->addNew<WCheckBox>( WString::tr("None") );
    m_noSpectrumBackground->setInline( false );
  }
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[SpectrumChartBackground] = cell->addNew<ColorSelect>(ColorSelect::AllowNoColor);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spectrum-area-desc"));
    
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-spectrum-margin"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noSpectrumMargin = cell->addNew<WCheckBox>( WString::tr("None") );
    m_noSpectrumMargin->setInline( false );
  }
	m_colorSelects[SpectrumChartMargins] = cell->addNew<ColorSelect>(ColorSelect::AllowNoColor);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-spectrum-margin-desc"));
	m_specMarginSameAsBackground = cell->addNew<WCheckBox>(WString::tr("ctwidget-same-as-spectrum"));
	m_specMarginSameAsBackground->setInline(false);


  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-time-gamma"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[TimeChartGammas] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-time-gamma-desc"));
  
  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-time-neutron"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_colorSelects[TimeChartNeutrons] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-time-neutron-desc"));
  
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-axis"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeAxisLines] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-axis-desc"));

  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-background"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noTimeBackground = cell->addNew<WCheckBox>( WString::tr("None") );
    m_noTimeBackground->setInline( false );
  }
  m_colorSelects[TimeChartBackground] = cell->addNew<ColorSelect>(ColorSelect::AllowNoColor);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-background-desc"));
  

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-margins"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  if( nativeColorSelect )
  {
    m_noTimeMargin = cell->addNew<WCheckBox>( WString::tr("None") );
    m_noTimeMargin->setInline( false );
  }
  m_colorSelects[TimeChartMargins] = cell->addNew<ColorSelect>(ColorSelect::AllowNoColor);
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-margins-desc"));
	m_timeMarginSameAsBackground = cell->addNew<WCheckBox>(WString::tr("ctwidget-same-as-background"));
	m_timeMarginSameAsBackground->setInline(false);
  
  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-text"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeChartText] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-text-desc"));

  
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-occ-lines"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[OccupancyIndicatorLines] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-occ-lines-desc"));

	
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-fore-highlight"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistoryForegroundHighlight] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-fore-highlight-desc"));

	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-back-highlight"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistoryBackgroundHighlight] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-back-highlight-desc"));
		
	++row;
	cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
	cell->addNew<WLabel>(WString::tr("ctwidget-time-sec-highlight"));
	cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
	m_colorSelects[TimeHistorySecondaryHighlight] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
	cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
	cell->addNew<WText>(WString::tr("ctwidget-time-sec-highlight-desc"));
	
  
  WContainerWidget *refLineContainer = addNew<WContainerWidget>();
  refLineContainer->addStyleClass( "RefLinesArea" );
  auto refLayoutOwned = std::make_unique<WGridLayout>();
  WGridLayout *refLayout = refLayoutOwned.get();
  refLayout->setContentsMargins( 0, 0, 0, 0 );
  refLayout->setVerticalSpacing( 0 );
  refLayout->setHorizontalSpacing( 0 );
  refLineContainer->setLayout( std::move(refLayoutOwned) );

  auto genericRefOwned = std::make_unique<WContainerWidget>();
  WContainerWidget *genericRef = genericRefOwned.get();
  genericRef->addStyleClass( "GenericRefLineColors" );
  refLayout->addWidget( std::move(genericRefOwned), 0, 0 );

  WText *genericTitle = genericRef->addNew<WText>( WString::tr("ctwidget-ref-line-colors") );
  genericTitle->addStyleClass( "ColorRefLineTitle" );
  genericTitle->setInline( false );

  WText *genericDesc = genericRef->addNew<WText>( WString::tr("ctwidget-ref-line-colors-desc") );
  genericDesc->addStyleClass( "ColorRefLineDesc" );
  genericDesc->setInline( false );

  table = genericRef->addNew<WTable>();
  table->addStyleClass( "ReferenceLinColorTable" );
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    cell = table->elementAt(i,0);
    cell->addNew<WText>( std::to_string(i+1) );
    cell->addStyleClass( "GenericRefLineNumber" );
    cell = table->elementAt(i,1);
    m_referenceLineColor[i] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
    m_referenceLineColor[i]->cssColorChanged().connect( this, [this, i](){ genericRefLineColorChangedCallback( i ); } );
  }

  auto specificRefOwned = std::make_unique<WContainerWidget>();
  WContainerWidget *specificRef = specificRefOwned.get();
  specificRef->addStyleClass( "SpecificRefLineColors" );
  refLayout->addWidget( std::move(specificRefOwned), 0, 1 );

  WText *specificTitle = specificRef->addNew<WText>( WString::tr("ctwidget-specific-ref-colors") );
  specificTitle->addStyleClass( "ColorRefLineTitle" );
  specificTitle->setInline( false );

  WText *specificDesc = specificRef->addNew<WText>( WString::tr("ctwidget-specific-ref-desc") );
  specificDesc->addStyleClass( "ColorRefLineDesc" );
  specificDesc->setInline( false );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  auto isoSuggestModel = std::make_shared<IsotopeNameFilterModel>();
  isoSuggestModel->addCustomSuggestPossibility( "background" );
  WSuggestionPopup *nuclideSuggest = addNew<WSuggestionPopup>( matcherJs, replacerJs );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  nuclideSuggest->addStyleClass( "nuclide-suggest" );

  IsotopeNameFilterModel::setQuickTypeFixHackjs( nuclideSuggest );

  isoSuggestModel->filter( "" );
  nuclideSuggest->setFilterLength( -1 );
  nuclideSuggest->setModel( isoSuggestModel );
  nuclideSuggest->filterModel().connect( isoSuggestModel.get(), &IsotopeNameFilterModel::filter );


  table = specificRef->addNew<WTable>();
  table->addStyleClass( "SpecificRefLineColorTable" );
  for( int i = 0; i < sm_numRefLineColors; ++i )
  {
    cell = table->elementAt(i,0);
    m_specificRefLineName[i] = cell->addNew<WLineEdit>();
    m_specificRefLineName[i]->setAttributeValue( "ondragstart", "return false" );
    nuclideSuggest->forEdit( m_specificRefLineName[i], PopupTrigger::Editing );
    IsotopeNameFilterModel::setEnterKeyMatchFixJs( nuclideSuggest, m_specificRefLineName[i] );
    m_specificRefLineName[i]->changed().connect( this, [this, i](){ specificRefLineSourceChangedCallback( i ); } );
    m_specificRefLineName[i]->enterPressed().connect( this, [this, i](){ specificRefLineSourceChangedCallback( i ); } );

    m_specificRefLineName[i]->addStyleClass( "SpecificRefLineInput" );
    cell = table->elementAt(i,1);
    cell->addStyleClass( "SpecificRefLineColorCell" );
    m_specificRefLineColor[i] = cell->addNew<ColorSelect>( ColorSelect::AllowNoColor );
    m_specificRefLineColor[i]->cssColorChanged().connect( this, [this, i](){ specificRefLineColorChangedCallback( i ); } );
    m_specificRefLineColor[i]->setColor( WColor() );
  }
  
  
  
  {
    auto dynamicRefCatsOwned = std::make_unique<WContainerWidget>();
    WContainerWidget *dynamicRefCats = dynamicRefCatsOwned.get();
    dynamicRefCats->addStyleClass( "DynamicRefLineColors" );
    refLayout->addWidget( std::move(dynamicRefCatsOwned), 1, 0 );

    WText *dynamicTitle = dynamicRefCats->addNew<WText>( WString::tr("ctwidget-dynamic-ref-lines") );
    dynamicTitle->addStyleClass( "ColorRefLineTitle" );
    dynamicTitle->setInline( false );

    WText *genericDesc = dynamicRefCats->addNew<WText>( WString::tr("ctwidget-dynamic-ref-desc") );
    genericDesc->addStyleClass( "ColorRefLineDesc" );
    genericDesc->setInline( false );

    table = dynamicRefCats->addNew<WTable>();
    table->addStyleClass( "ReferenceLinColorTable" );

    int dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-medical"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineMedical] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});


    dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-industrial"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineIndustrial] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});


    dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-norm"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineNorm] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});

    dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-snm"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineSnm] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});

    dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-common"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineCommon] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});

    dynamic_ref_row = table->rowCount();
    cell = table->elementAt(dynamic_ref_row, 0);
    cell->addStyleClass( "CTRowLabel" );
    cell->addNew<WLabel>(WString::tr("ctwidget-other"));
    cell = table->elementAt(dynamic_ref_row, 1);
    cell->addStyleClass( "CTSelect" );
    m_colorSelects[DynamicRefLineOther] = cell->addNew<ColorSelect>(Wt::WFlags<ColorSelect::ColorSelectOptions>{});
  }
  
  
  
  row = 0;
  table = addNew<WTable>();
  table->addStyleClass( "ColorThemeOtherOptionsTable" );

  cell = table->elementAt(row, 0);
  cell->setColumnSpan(3);
  cell->setContentAlignment(Wt::AlignmentFlag::Center);
  cell->addStyleClass( "ThemePeakLabelTitle" );
  cell->addNew<WText>(WString::tr("ctwidget-other-options"));

  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-peak-label-size"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_peakLabelFontSize = cell->addNew<WComboBox>();
  m_peakLabelFontSize->addItem( WString::tr("Default") );
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
  cell->addNew<WText>(WString::tr("ctwidget-peak-label-size-desc"));


  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-peak-label-angle"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_peakLabelAngle = cell->addNew<NativeFloatSpinBox>();
  m_peakLabelAngle->setSpinnerHidden( true );
  m_peakLabelAngle->setRange( -180, 270 );
  m_peakLabelAngle->setWidth( 50 );
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-peak-label-angle-desc"));


  ++row;
  cell = table->elementAt(row, 0);
  cell->addStyleClass( "CTRowLabel" );
  cell->addNew<WLabel>(WString::tr("ctwidget-log-y-min"));
  cell = table->elementAt(row, 1);
  cell->addStyleClass( "CTSelect" );
  m_logYAxisMin = cell->addNew<NativeFloatSpinBox>();
  m_logYAxisMin->setSpinnerHidden( true );
  m_logYAxisMin->setRange( 1.0E-8, 1000 );
  m_logYAxisMin->setWidth( 50 );
  cell = table->elementAt(row, 2);
  cell->addStyleClass( "CTRowDesc" );
  cell->addNew<WText>(WString::tr("ctwidget-log-y-min-desc"));
  
  
  m_themeTitle->changed().connect( this, [this](){ titleChangedCallback(); } );
  m_themeDescription->changed().connect( this, [this](){ descriptionChangedCallback(); } );
  m_peaksTakeRefLineColor->checked().connect( this, [this](){ peaksTakeRefLineColorChangedCallback(); } );
  m_peaksTakeRefLineColor->unChecked().connect( this, [this](){ peaksTakeRefLineColorChangedCallback(); } );

  m_peakLabelFontSize->activated().connect( this, [this]( int ){ peakLabelFontSizeChanged(); } );
  m_peakLabelAngle->valueChanged().connect( this, [this]( double ){ peakLabelRotationChanged(); } );
  m_logYAxisMin->valueChanged().connect( this, [this]( double ){ logYAxisMinValueChanged(); } );

  m_nonChartAreaCssTheme->changed().connect( this, [this](){ nonChartAreaThemeChanged(); } );

  m_specMarginSameAsBackground->checked().connect( this, [this](){ newColorSelectedCallback( SpectrumChartMargins ); } );
  m_specMarginSameAsBackground->unChecked().connect( this, [this](){ newColorSelectedCallback( SpectrumChartMargins ); } );

  m_timeMarginSameAsBackground->checked().connect( this, [this](){ newColorSelectedCallback( TimeChartMargins ); } );
  m_timeMarginSameAsBackground->unChecked().connect( this, [this](){ newColorSelectedCallback( TimeChartMargins ); } );

  if( m_noSpectrumBackground )
  {
    m_noSpectrumBackground->checked().connect( this, [this](){ newColorSelectedCallback( SpectrumChartBackground ); } );
    m_noSpectrumBackground->unChecked().connect( this, [this](){ newColorSelectedCallback( SpectrumChartBackground ); } );
  }

  if( m_noSpectrumMargin )
    m_noSpectrumMargin->changed().connect( this, [this](){ newColorSelectedCallback( TimeChartMargins ); } );

  if( m_noTimeBackground )
  {
    m_noTimeBackground->checked().connect( this, [this](){ newColorSelectedCallback( TimeChartBackground ); } );
    m_noTimeBackground->unChecked().connect( this, [this](){ newColorSelectedCallback( TimeChartBackground ); } );
  }

  if( m_noTimeMargin )
  {
    m_noTimeMargin->checked().connect( this, [this](){ newColorSelectedCallback( TimeChartMargins ); } );
    m_noTimeMargin->unChecked().connect( this, [this](){ newColorSelectedCallback( TimeChartMargins ); } );
  }


	for( SelectableColor color = SelectableColor(0); color < NumSelectableColors; color = SelectableColor(color + 1) )
	{
		m_colorSelects[color]->cssColorChanged().connect( this, [this, color](){ newColorSelectedCallback( color ); } );
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
		m_themeTitle->setText(WString::tr("ctwidget-empty"));
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

  m_colorSelects[AppBackground]->setColor( theme->appBackgroundColor );
  m_colorSelects[AppText]->setColor( theme->appTextColor );
  m_colorSelects[AppBorder]->setColor( theme->appBorderColor );
  m_colorSelects[AppLink]->setColor( theme->appLinkColor );
  m_colorSelects[AppLabel]->setColor( theme->appLabelColor );
  m_colorSelects[AppInputBackground]->setColor( theme->appInputBackground );
  m_colorSelects[AppButtonBackground]->setColor( theme->appButtonBackground );
  m_colorSelects[AppButtonBorder]->setColor( theme->appButtonBorderColor );
  m_colorSelects[AppButtonText]->setColor( theme->appButtonTextColor );
  m_colorSelects[AppMenuBarBackground]->setColor( theme->appMenuBarBackground );
  m_colorSelects[AppMenuBarActiveColor]->setColor( theme->appMenuBarActiveColor );
  m_colorSelects[AppMenuBarHoverColor]->setColor( theme->appMenuBarHoverColor );

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

    case AppBackground:
      m_currentTheme->appBackgroundColor = m_colorSelects[color]->color();
      break;

    case AppText:
      m_currentTheme->appTextColor = m_colorSelects[color]->color();
      break;

    case AppBorder:
      m_currentTheme->appBorderColor = m_colorSelects[color]->color();
      break;

    case AppLink:
      m_currentTheme->appLinkColor = m_colorSelects[color]->color();
      break;

    case AppLabel:
      m_currentTheme->appLabelColor = m_colorSelects[color]->color();
      break;

    case AppInputBackground:
      m_currentTheme->appInputBackground = m_colorSelects[color]->color();
      break;

    case AppButtonBackground:
      m_currentTheme->appButtonBackground = m_colorSelects[color]->color();
      break;

    case AppButtonBorder:
      m_currentTheme->appButtonBorderColor = m_colorSelects[color]->color();
      break;

    case AppButtonText:
      m_currentTheme->appButtonTextColor = m_colorSelects[color]->color();
      break;

    case AppMenuBarBackground:
      m_currentTheme->appMenuBarBackground = m_colorSelects[color]->color();
      break;

    case AppMenuBarActiveColor:
      m_currentTheme->appMenuBarActiveColor = m_colorSelects[color]->color();
      break;

    case AppMenuBarHoverColor:
      m_currentTheme->appMenuBarHoverColor = m_colorSelects[color]->color();
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


