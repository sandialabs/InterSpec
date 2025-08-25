#ifndef ColorThemeWidget_h
#define ColorThemeWidget_h
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

#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

//Forward declarations
namespace Wt
{
  class WCheckBox;
  class WLineEdit;
  class WComboBox;
}

struct ColorTheme;
class ColorSelect;
class NativeFloatSpinBox;

/** Widget to display, and allow editing a ColorTheme.
    Does not handle saving or loading ColorThemes to/from database; see
    #ColorThemeWindow
 */
class ColorThemeWidget : public Wt::WContainerWidget
{
public:
	ColorThemeWidget( Wt::WContainerWidget *parent = nullptr );
	virtual ~ColorThemeWidget();

	/** A copy of the theme passed in will be made, and modified. */
  void setTheme( const ColorTheme *theme, const bool modifieable );

  Wt::Signal<> &edited();
  
	enum SelectableColor 
	{  
		ForegroundLine,
    BackgroundLine,
		SecondaryLine,
		DefaultPeakLine,
		SpectrumAxis,
		SpectrumChartText,
		SpectrumChartBackground,
		SpectrumChartMargins,
    TimeChartGammas,
    TimeChartNeutrons,
		TimeAxisLines,
		TimeChartText,
		TimeChartBackground,
		TimeChartMargins,
		OccupancyIndicatorLines,
		TimeHistoryForegroundHighlight,
		TimeHistoryBackgroundHighlight,
		TimeHistorySecondaryHighlight,
		KineticRefLineDefault,
		NumSelectableColors
	};//enum class EditableColor

	/** */
	void newColorSelectedCallback(const SelectableColor color);
  void titleChangedCallback();
  void descriptionChangedCallback();
  void nonChartAreaThemeChanged();
  void peaksTakeRefLineColorChangedCallback();
  void peakLabelFontSizeChanged();
  void peakLabelRotationChanged();
  void logYAxisMinValueChanged();

  void genericRefLineColorChangedCallback( const int num );
  void specificRefLineSourceChangedCallback( const int num );
  void specificRefLineColorChangedCallback( const int num );
  
  std::unique_ptr<const ColorTheme> m_origTheme;
  bool m_themeModifiable;
  std::unique_ptr<ColorTheme> m_currentTheme; //only created when there is a modification

	Wt::WLineEdit       *m_themeTitle;
	Wt::WLineEdit       *m_themeDescription;

  //ToDo: right this combo box only lists names given by ColorTheme::predefinedThemeName(),
  //      but what we should *really* do is list options in InterSpec_resources/themes
  Wt::WComboBox *m_nonChartAreaCssTheme;
  
	ColorSelect *m_colorSelects[NumSelectableColors];

	Wt::WCheckBox *m_specMarginSameAsBackground;
	Wt::WCheckBox *m_timeMarginSameAsBackground;

  Wt::WCheckBox *m_noSpectrumBackground;
  Wt::WCheckBox *m_noSpectrumMargin;

  Wt::WCheckBox *m_noTimeBackground;
  Wt::WCheckBox *m_noTimeMargin;

	Wt::WCheckBox *m_peaksTakeRefLineColor;

  Wt::WComboBox *m_peakLabelFontSize;
  NativeFloatSpinBox *m_peakLabelAngle;
  NativeFloatSpinBox *m_logYAxisMin;
  
  static const int sm_numRefLineColors = 12;
  ColorSelect *m_referenceLineColor[sm_numRefLineColors];
  
  Wt::WLineEdit *m_specificRefLineName[sm_numRefLineColors];
  ColorSelect *m_specificRefLineColor[sm_numRefLineColors];
  
  Wt::Signal<> m_edited;
};//class ColorThemeWidget


#endif //ColorThemeWidget_h
