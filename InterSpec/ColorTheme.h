#ifndef ColorTheme_h
#define ColorTheme_h
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
#include <memory>
#include <vector>
#include <string>

#include <Wt/WColor>
#include <Wt/WString>
#include <Wt/WDateTime>


class ColorThemeInfo;

/** Color theme of the app.
  Currently only the coloring of the charts is supported.
 
 ToDo:
 -Allow specifying peak fill and line colors seperately.
 -Allow user to pick alpha (now that spectrum.js can be used) for chart things.
 
 Note: Alpha channel is not supported currently, due to issue with Wt<4.0.4
 https://redmine.webtoolkit.eu/issues/6702
 (could implement a workaround for now...)
 */
struct ColorTheme
{
  /** The version of JSON format currently written by #toJson. Increment this
   value when introducing breaking (ie, non backwards compatible) changes to the
   json encoding of #ColorTheme.
   */
  static const char * const sm_color_theme_json_version;  //Currently "1"
  
  /** A listing of the predifeined color themes.
      Note that a userPreference for "ColorThemeIndex" that has a negative value
      indicates to use one of these predifined values, where the absolute value
      tells which one.
   */
  enum PredefinedColorTheme
  {
    DefaultColorTheme = 1,
    DarkColorTheme = 2,
    NumPredefinedColorTheme = 3
  };//enum PredefinedColorTheme
  
  /** */
  static std::string predefinedThemeName( const PredefinedColorTheme theme );
  
  /** Returns the predifined color theme.
   Returns nullptr if asked for theme is invalid.
   */
  static std::unique_ptr<ColorTheme> predefinedTheme( const PredefinedColorTheme theme );
  
  //ToDo: include mechanism to use JSON files to add to this.
  static std::vector<std::unique_ptr<ColorTheme>> predefinedThemes();

  /** Initialized to InterSpec defaults. */
  ColorTheme();
  
  /** Sets the information from an instance of a ColorThemeInfo.
      Throws exception if error.
   */
  void setFromDataBase( const ColorThemeInfo &db );
  
  /** Serializes this class to a ColorThemeInfo. */
  void setToDataBase( ColorThemeInfo &db ) const;
  
  /** Write the contents of a #ColorTheme object to JSON.
   A different #ColorTheme object can be intiated from the JSON using #fromJson
   to get an identical #ColorTheme instance.
   */
  static std::string toJson( const ColorTheme &info );
  
  /** Sets the passed in #ColorTheme object to values specified in the JSON
     passed in.  If a value is not specified in the JSON string, then the
     default value for that member variable will be used.
   
     Throws exception if error.
   */
  static void fromJson( const std::string &json, ColorTheme &info );
  
  /** Not serialized to/from JSON, but used for saving to database.  Will be -1
      if not saved to DB.
   */
  long long dbIndex;
  
  Wt::WString theme_name;
  Wt::WString theme_description;
  
  Wt::WDateTime creation_time;
  Wt::WDateTime modified_time;
  
  /** CSS theme to styla the rest of the app besides the charts.
      To create one of these themes you need to at least creat a CSS file at:
        InterSpec_resources/themes/<theme name>/<theme name>.css
   */
  std::string nonChartAreaTheme;
  
  /** Line color of foreground spectrum. */
  Wt::WColor foregroundLine;
  
  /** Line color of background spectrum. */
  Wt::WColor backgroundLine;
  
  /** Line color of secondary spectrum. */
  Wt::WColor secondaryLine;
  
  /** The axis lines, and ticks for the spectrum chart. */
  Wt::WColor spectrumAxisLines;
  
  /** The background color for the spectrum chart.  If this is specified, but
   #spectrumChartMargins is not, then the entire chart area will be set to
   this color.  However, if #spectrumChartMargins is specified, then,
   #spectrumChartBackground will be used to color the square area falling
   inside the chart axis.
   If #spectrumChartBackground is not given, then the background will be white.
   */
  Wt::WColor spectrumChartBackground;
  
  /** If specified, this color will be used for the area of the chart outside
   the axis.  If this quantitiy is not specified, but #spectrumChartBackground
   is, then this area will be colored according to #spectrumChartBackground.
   */
  Wt::WColor spectrumChartMargins;
  
  /** The color for all text on the chart, besides peak labels.  So axis tick
   values, axis titles, and text that may appear in the chart.
   */
  Wt::WColor spectrumChartText;
  
  /** The peak label size.
   Example values: "8px", "smaller", "12", "10px", "x-small", etc.
   Empty string gives default size
   */
  Wt::WString spectrumPeakLabelSize;
  
  /** Peak label rotation angle, in degrees.
   
   A negative value rotates it the direction you probably want.
   A value of 0 is horizontal, a value of -90 is vertical (i.e. up-and-down).  Only tested [0,-90]
   */
  double spectrumPeakLabelRotation;
  
  /** The minimum value the y-axis will go down to, if there there is a channel count
   with zero or less counts.
   Must be greater than zero.
   Untested with anything besides powers of 10, e.g., 0.1, 0.001, etc
   */
  double spectrumLogYAxisMin;
  
  /** The color for the line that indicates the gamma counts on the time chart
   */
  Wt::WColor timeChartGammaLine;
  
  /** The color for the line that indicates the neutron counts on the time chart
   */
  Wt::WColor timeChartNeutronLine;
  
  /** Equivalent of #spectrumAxisLines, but for the time chart (passthrough and
   search-mode data).
  */
  Wt::WColor timeAxisLines;
  
  /** Equivalent of #spectrumChartBackground, but for the time chart
      (passthrough and search-mode data).
   */
  Wt::WColor timeChartBackground;
  
  /** Equivalent of #spectrumChartMargins, but for the time chart (passthrough
   and search-mode data).
   */
  Wt::WColor timeChartMargins;
  
  /** Equivalent of #spectrumChartText, but for the time chart (passthrough and
   search-mode data).
   */
  Wt::WColor timeChartText;
  
  /** Vertical lines (and text) that indicate the begining and end of occupied
   regions of the time history chart; used for some RPM detection systems.
   */
  Wt::WColor occupancyIndicatorLines;
  
  //ToDo: range count sum highlight area, and erase peaks highlight area, etc.
  
  
  Wt::WColor timeHistoryForegroundHighlight;
  Wt::WColor timeHistoryBackgroundHighlight;
  Wt::WColor timeHistorySecondaryHighlight;
  
  /** The color for peaks not associated with any reference lines, or if
   #peaksTakeOnReferenceLineColor is false, all peaks.
   */
  Wt::WColor defaultPeakLine;
  
  /** If false, then peaks will always be set to #defaultPeakLine.
      If true, then peaks associated with a reference line will be set to that
      reference lines color, and if not associated with a refline, then will
      be the #defaultPeakLine color.
   */
  bool peaksTakeOnReferenceLineColor;
  
  /** For reference lines, whose source is not in #referenceLineColorForSources,
   and whos source is not already represented by a peak in the spectrum, will be
   selected from the first availble #referenceLineColor entry not already used.
   */
  std::vector<Wt::WColor> referenceLineColor;
  
  /** Allows mapping of a specific nuclide/x-ray/reaction to a color always.
      The key will be something like 'U235', 'Co65', 'Fe', etc.,
   */
  std::map<std::string,Wt::WColor> referenceLineColorForSources;
};//struct ColorTheme


#endif //ColorTheme_h
