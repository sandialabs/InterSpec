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
#include <regex>
#include <chrono>
#include <algorithm>
#include <vector>
#include <string>
#include <typeinfo>

#include <Wt/WColor.h>
#include <Wt/WString.h>
#include <Wt/WDateTime.h>
#include <Wt/Json/Value.h>
#include <Wt/Json/Array.h>
#include <Wt/Json/Parser.h>
#include <Wt/Json/Object.h>
#include <Wt/Json/Serializer.h>

#include "InterSpec/AppUtils.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecUser.h"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

using namespace std;
using namespace Wt;

const char * const ColorTheme::sm_color_theme_json_version = "1";
const char * const ColorTheme::sm_dynamic_ref_line_medical_color = "#FF0000";
const char * const ColorTheme::sm_dynamic_ref_line_industrial_color = "#7e7e80";
const char * const ColorTheme::sm_dynamic_ref_line_norm_color = "#39a004";
const char * const ColorTheme::sm_dynamic_ref_line_snm_color = "#ff7600";
const char * const ColorTheme::sm_dynamic_ref_line_common_color = "#0033ff";
const char * const ColorTheme::sm_dynamic_ref_line_other_color = "#8055ef";

namespace
{
  /*
  void minimal_test()
  {
    ColorTheme theme;
    theme.backgroundLine = Wt::WColor(10,11,12,13);
    theme.foregroundLine = Wt::WColor(99,98,97,96);
    theme.referenceLineColorForSources["Th999"] = Wt::WColor(21,22,23,24);
    const std::string themestr = ColorTheme::toJson(theme);
    ColorTheme duplicate;
    ColorTheme::fromJson(themestr, duplicate);
    std::cout << "Original:\n" << themestr << "\n" << std::endl;
    std::cout << "Duplicate:\n" << ColorTheme::toJson(duplicate) << std::endl;
    assert( themestr == ColorTheme::toJson(duplicate) );
  }//void minimal_test()
   */

  /** Adds the `--interspec-<name>: <value>;` declarations found in the CSS file at `filePath`
   to `tokens` (later files override earlier ones).  Does nothing if the file does not exist.
   */
  void add_css_token_declarations( const string &filePath, map<string,string> &tokens )
  {
    if( !SpecUtils::is_file(filePath) )
      return;

    string css;
    try
    {
      css = AppUtils::file_contents( filePath );
    }catch( const std::exception &e )
    {
      // Only cosmetic (the editor shows the theme's default beside each colour), so an unreadable
      //  file must not take the dialog down with it.
      cerr << "Failed to read theme CSS '" << filePath << "': " << e.what() << endl;
      return;
    }

    // Strip comments, so a commented-out declaration is not picked up
    css = std::regex_replace( css, std::regex("/\\*[\\s\\S]*?\\*/"), "" );

    const std::regex declaration( "--interspec-([a-z0-9-]+)\\s*:\\s*([^;{}]+);" );
    const std::sregex_iterator end_it;
    for( std::sregex_iterator it( begin(css), end(css), declaration ); it != end_it; ++it )
    {
      string value = (*it)[2].str();
      SpecUtils::trim( value );
      tokens[(*it)[1].str()] = value;
    }
  }//void add_css_token_declarations(...)
}//namespace


const std::vector<ColorTheme::AppColorToken> &ColorTheme::appColorTokens()
{
  // Same order and grouping as InterSpec_resources/themes/default/default.css.
  static const vector<AppColorToken> tokens{
    // token name,                  legacy JSON key,       i18n stem,            editor group
    { "background-color",          "backgroundColor",     "background",         "surfaces" },
    { "input-background",          "inputBackground",     "input-bg",           "surfaces" },
    { "panel-background",          nullptr,               "panel-bg",           "surfaces" },
    { "header-background",         nullptr,               "header-bg",          "surfaces" },
    { "alt-row-background",        nullptr,               "alt-row-bg",         "surfaces" },
    { "popup-background",          nullptr,               "popup-bg",           "surfaces" },
    { "tooltip-background",        nullptr,               "tooltip-bg",         "surfaces" },
    { "tooltip-text-color",        nullptr,               "tooltip-text",       "surfaces" },
    { "tooltip-border-color",      nullptr,               "tooltip-border",     "surfaces" },
    { "dialog-cover-color",        nullptr,               "dialog-cover",       "surfaces" },
    { "titlebar-background",       nullptr,               "titlebar-bg",        "surfaces" },
    { "titlebar-hover-background", nullptr,               "titlebar-hover-bg",  "surfaces" },
    { "text-color",                "textColor",           "text",               "text" },
    { "label-color",               "labelColor",          "label",              "text" },
    { "secondary-text-color",      nullptr,               "secondary-text",     "text" },
    { "fainter-text-color",        nullptr,               "fainter-text",       "text" },
    { "disabled-text-color",       nullptr,               "disabled-text",      "text" },
    { "disabled-background",       nullptr,               "disabled-bg",        "text" },
    { "border-color",              "borderColor",         "border",             "lines" },
    { "outline-color",             nullptr,               "outline",            "lines" },
    { "link-color",                "linkColor",           "link",               "interactive" },
    { "accent-color",              nullptr,               "accent",             "interactive" },
    { "highlight-background",      nullptr,               "highlight-bg",       "interactive" },
    { "selected-background",       nullptr,               "selected-bg",        "interactive" },
    { "selected-text-color",       nullptr,               "selected-text",      "interactive" },
    { "hover-background",          nullptr,               "hover-bg",           "interactive" },
    { "button-background",         "buttonBackground",    "button-bg",          "interactive" },
    { "button-hover-background",   nullptr,               "button-hover-bg",    "interactive" },
    { "button-border-color",       "buttonBorderColor",   "button-border",      "interactive" },
    { "button-text-color",         "buttonTextColor",     "button-text",        "interactive" },
    { "menubar-background",        "menuBarBackground",   "menubar-bg",         "interactive" },
    { "menubar-active-color",      "menuBarActiveColor",  "menubar-active",     "interactive" },
    { "menubar-hover-color",       "menuBarHoverColor",   "menubar-hover",      "interactive" },
    { "warning-color",             nullptr,               "warning",            "semantic" },
    { "warning-background",        nullptr,               "warning-bg",         "semantic" },
    { "warning-border-color",      nullptr,               "warning-border",     "semantic" },
    { "error-color",               nullptr,               "error",              "semantic" },
    { "error-background",          nullptr,               "error-bg",           "semantic" },
    { "error-border-color",        nullptr,               "error-border",       "semantic" },
    { "success-color",             nullptr,               "success",            "semantic" },
    { "success-background",        nullptr,               "success-bg",         "semantic" },
    { "success-border-color",      nullptr,               "success-border",     "semantic" },
  };

  return tokens;
}//appColorTokens()


const ColorTheme::AppColorToken *ColorTheme::appColorToken( const std::string &name )
{
  for( const AppColorToken &token : appColorTokens() )
  {
    if( name == token.name )
      return &token;
  }

  return nullptr;
}//appColorToken( name )


std::map<std::string,std::string> ColorTheme::cssThemeTokenDefaults( const std::string &resourcesDir,
                                                                    const std::string &cssTheme )
{
  map<string,string> tokens;
  const string themesDir = SpecUtils::append_path( resourcesDir, "themes" );

  const string baseDir = SpecUtils::append_path( themesDir, "default" );
  add_css_token_declarations( SpecUtils::append_path( baseDir, "default.css" ), tokens );

  if( !cssTheme.empty() && !SpecUtils::iequals_ascii(cssTheme, "default") )
  {
    const string themeDir = SpecUtils::append_path( themesDir, cssTheme );
    add_css_token_declarations( SpecUtils::append_path( themeDir, cssTheme + ".css" ), tokens );
  }

  return tokens;
}//cssThemeTokenDefaults(...)


std::vector<std::string> ColorTheme::availableCssThemes( const std::string &resourcesDir )
{
  vector<string> themes{ "default" };

  const string themesDir = SpecUtils::append_path( resourcesDir, "themes" );
  vector<string> dirs = SpecUtils::ls_directories_in_directory( themesDir );
  std::sort( begin(dirs), end(dirs) );

  for( const string &dir : dirs )
  {
    const string name = SpecUtils::filename( dir );
    const string cssFile = SpecUtils::append_path( SpecUtils::append_path(themesDir, name), name + ".css" );
    if( (name != "default") && SpecUtils::is_file(cssFile) )
      themes.push_back( name );
  }

  return themes;
}//availableCssThemes(...)


std::string ColorTheme::predefinedThemeName( const PredefinedColorTheme theme )
{
  switch( theme )
  {
    case DefaultColorTheme: return "default";
    case DarkColorTheme: return "dark";
    case NumPredefinedColorTheme: return "NumPredefinedColorTheme";
  }
  return "InvalidPredefinedColorTheme";
}//std::string ColorTheme::predefinedThemeName( const PredefinedColorTheme theme )


std::unique_ptr<ColorTheme> ColorTheme::predefinedTheme( const PredefinedColorTheme theme )
{
  unique_ptr<ColorTheme> themePtr;
  
  switch( theme )
  {
    case DefaultColorTheme:
      themePtr.reset( new ColorTheme() );
      break;
      
    case DarkColorTheme:
    {
      const string darkJson = R"delim(
      {
        "JsonVersion" : "1",
        "created" : "2018-12-24T07:02:39+0000",
        "defaultPeakLineColor" : "#cecfd0",
        "description" : "Dark InterSpec color scheme.",
        "modified" : "2018-12-26T02:58:03+0000",
        "name" : "Dark",
        "nonChartArea": { "cssTheme" : "dark" },
        "peaksTakeOnReferenceLineColor" : true,
        "referenceLines" : {
          "lineColors" : ["#c0c0c0", "#ffff99", "#b8d9f0", "#9933FF", "#FF66FF", "#CC3333", "#FF6633", "#FFFF99", "#CCFFCC", "#0000CC", "#666666", "#003333"],
          "specificSources" : {
            "Ba131" : "rgb(192,192,192)",
            "background" : "#967f55"
          },
          "dynamicRefLineMedicalColor" : "#FF0000",
          "dynamicRefLineIndustrialColor" : "#D7D7DA", 
          "dynamicRefLineNormColor" : "#39a004",
          "dynamicRefLineSnmColor" : "#ff7600",
          "dynamicRefLineCommonColor" : "#9da2f2",
          "dynamicRefLineOtherColor" : "#8055ef"
        },
        "spectrum" : {
          "axisLines" : "#cfced2",
          "backgroundColor" : "#a9ebdd",
          "chartText" : "#cfced2",
          "foregroundColor" : "#cfced2",
          "secondaryColor" : "#fd8273"
        },
        "timeChart" : {
          "axisLines" : "#cfced2",
          "backgroundHighlightColor" : "rgb(0,255,255)",
          "chartText" : "#cfced2",
          "foregroundHighlightColor" : "#f6f6ad",
          "gammaLine" : "#cfced2",
          "neutronLine" : "rgb(0,128,0)",
          "occIndicatorLines" : "rgb(128,128,128)",
          "secondaryHighlightColor" : "rgb(0,128,0)"
        }
      }
      )delim";
      
      themePtr.reset( new ColorTheme() );
      fromJson( darkJson, *themePtr );
      
      break;
    }//case DarkColorTheme:
      
    case NumPredefinedColorTheme:
      break;
  }
  
  if( themePtr )
    themePtr->dbIndex = -static_cast<long long>( theme );
  
  return themePtr;
}//predefinedTheme(...)


vector<unique_ptr<ColorTheme>> ColorTheme::predefinedThemes()
{
	vector<unique_ptr<ColorTheme>> answer;
  
  for( PredefinedColorTheme t = DefaultColorTheme;
      t < NumPredefinedColorTheme;
      t = PredefinedColorTheme(t+1) )
	answer.push_back( predefinedTheme(t) );
  
	return answer;
}

ColorTheme::ColorTheme()
{
  dbIndex = -1;

  theme_name = "Default";
  theme_description = "Default InterSpec color scheme.";
  creation_time = WDateTime( WDate(2018,11,1), WTime(22,56) );  //when I originally created this class
  modified_time = WDateTime( WDate(2022,8,3), WTime(14,00) );   //when I slightly modified it

  // The base (light) CSS theme is always loaded, and `appColors` stays empty: the app colours
  //  come from InterSpec_resources/themes/default/default.css unless a theme overrides them.
  nonChartAreaTheme = "";

  foregroundLine = Wt::WColor(0x00, 0x00, 0x00); //Wt::Wt::StandardColor::Black
  backgroundLine = Wt::WColor(0x00, 0xff, 0xff); //Wt::GlobalColor::cyan
  secondaryLine = Wt::WColor(0x00, 0x80, 0x00);  //Wt::Wt::StandardColor::DarkGreen
  
  timeHistoryForegroundHighlight = Wt::WColor( 255, 255, 0, 155 );
  timeHistoryBackgroundHighlight = Wt::WColor( 0, 255, 255, 75 );
  timeHistorySecondaryHighlight = Wt::WColor( 0, 128, 0, 75 );
  
  defaultPeakLine = WColor( 0, 51, 255 );
  
  spectrumAxisLines = WColor( Wt::StandardColor::Black );
  spectrumChartBackground = WColor();
  spectrumChartMargins = WColor();
  spectrumChartText = WColor( Wt::StandardColor::Black );
  
  spectrumPeakLabelSize = WString();
  spectrumPeakLabelRotation = 0.0;
  spectrumLogYAxisMin = 0.1;
  
  timeChartGammaLine = WColor( Wt::StandardColor::Black );
  timeChartNeutronLine = Wt::WColor( Wt::StandardColor::DarkGreen );
  
  timeAxisLines = WColor( Wt::StandardColor::Black );
  timeChartBackground = WColor();
  timeChartMargins = WColor();
  timeChartText = WColor( Wt::StandardColor::Black );
  
  occupancyIndicatorLines = WColor( 128, 128, 128 ); //alpha channel of 75
  
  peaksTakeOnReferenceLineColor = true;
  referenceLineColor = std::vector<Wt::WColor>{
    {"#00b9ff"},
    {"#006600"},
    {"#cc3333"},
    {"#9933FF"},
    {"#FF66FF"},
    {"#830808"},
    {"#FF6633"},
    {"#F1C232"},
    {"#CCFFCC"},
    {"#0000CC"},
    {"#666666"},
    {"#003333"}
  };
  
  referenceLineColorForSources["Ba131"] = WColor( "#800080" ); //Just an example
  referenceLineColorForSources["background"] = WColor( "#967f55" );  //brownish
  
  dynamicRefLineMedicalColor = WColor( sm_dynamic_ref_line_medical_color );
  dynamicRefLineIndustrialColor = WColor( sm_dynamic_ref_line_industrial_color );
  dynamicRefLineNormColor = WColor( sm_dynamic_ref_line_norm_color );
  dynamicRefLineSnmColor = WColor( sm_dynamic_ref_line_snm_color );
  dynamicRefLineCommonColor = WColor( sm_dynamic_ref_line_common_color );
  dynamicRefLineOtherColor = WColor( sm_dynamic_ref_line_other_color );
}//ColorTheme() constructor


std::string ColorTheme::toJson( const ColorTheme &info )
{
  Json::Object base;
  base["JsonVersion"] = WString(sm_color_theme_json_version);
  base["name"] = info.theme_name;
  base["description"] = info.theme_description;
  base["created"] = info.creation_time.toString( "yyyy-MM-ddThh:mm:ssZ" ); //"2018-11-01T22:56:00+0000"
  base["modified"] = info.modified_time.toString( "yyyy-MM-ddThh:mm:ssZ" );
  
  
  Json::Object &nonChartArea = base["nonChartArea"] = Json::Value(Json::Type::Object);
  if( info.nonChartAreaTheme.empty() )
    nonChartArea["cssTheme"] = WString("default");
  else
    nonChartArea["cssTheme"] = WString( info.nonChartAreaTheme );

  // Only explicit overrides are written, keyed by CSS token name.  Alpha is kept, since the
  //  hover/highlight colours are translucent.
  for( const auto &nameAndColor : info.appColors )
  {
    if( !nameAndColor.second.isDefault() )
      nonChartArea[nameAndColor.first] = WString( nameAndColor.second.cssText(true) );
  }

  Json::Object &spectrum = base["spectrum"] = Json::Value(Json::Type::Object);
  if( !info.foregroundLine.isDefault() )
    spectrum["foregroundColor"] = WString( info.foregroundLine.cssText(false) );
  if( !info.backgroundLine.isDefault() )
    spectrum["backgroundColor"] = WString( info.backgroundLine.cssText(false) );
  if( !info.secondaryLine.isDefault() )
    spectrum["secondaryColor"] = WString( info.secondaryLine.cssText(false) );
  if( !info.spectrumAxisLines.isDefault() )
   spectrum["axisLines"] = WString( info.spectrumAxisLines.cssText(false) );
  if( !info.spectrumChartBackground.isDefault() )
    spectrum["chartBackground"] = WString( info.spectrumChartBackground.cssText(false) );
  if( !info.spectrumChartMargins.isDefault() )
    spectrum["chartMargins"] = WString( info.spectrumChartMargins.cssText(false) );
  if( !info.spectrumChartText.isDefault() )
    spectrum["chartText"] = WString( info.spectrumChartText.cssText(false) );
  if( !info.spectrumPeakLabelSize.empty() )
    spectrum["peakLabelSize"] = info.spectrumPeakLabelSize;
  if( info.spectrumPeakLabelRotation != 0.0 )
    spectrum["peakLabelRotation"] = info.spectrumPeakLabelRotation;
  if( info.spectrumLogYAxisMin > 0.0 )
    spectrum["logYAxisMin"] = info.spectrumLogYAxisMin;
  
  Json::Object &timechart = base["timeChart"] = Json::Value(Json::Type::Object);
  if( !info.timeHistoryForegroundHighlight.isDefault() )
    timechart["foregroundHighlightColor"] = WString( info.timeHistoryForegroundHighlight.cssText(false) );
  if( !info.timeHistoryBackgroundHighlight.isDefault() )
    timechart["backgroundHighlightColor"] = WString( info.timeHistoryBackgroundHighlight.cssText(false) );
  if( !info.timeHistorySecondaryHighlight.isDefault() )
    timechart["secondaryHighlightColor"] = WString( info.timeHistorySecondaryHighlight.cssText(false) );
  if( !info.timeChartGammaLine.isDefault() )
    timechart["gammaLine"] = WString( info.timeChartGammaLine.cssText(false) );
  if( !info.timeChartNeutronLine.isDefault() )
    timechart["neutronLine"] = WString( info.timeChartNeutronLine.cssText(false) );
  if( !info.timeAxisLines.isDefault() )
    timechart["axisLines"] = WString( info.timeAxisLines.cssText(false) );
  if( !info.timeChartBackground.isDefault() )
    timechart["chartBackground"] = WString( info.timeChartBackground.cssText(false) );
  if( !info.timeChartMargins.isDefault() )
    timechart["chartMargins"] = WString( info.timeChartMargins.cssText(false) );
  if( !info.timeChartText.isDefault() )
    timechart["chartText"] = WString( info.timeChartText.cssText(false) );
  if( !info.occupancyIndicatorLines.isDefault() )
    timechart["occIndicatorLines"] = WString( info.occupancyIndicatorLines.cssText(false) );
  
  
  base["peaksTakeOnReferenceLineColor"] = info.peaksTakeOnReferenceLineColor;
  if( !info.defaultPeakLine.isDefault() )
    base["defaultPeakLineColor"] = WString( info.defaultPeakLine.cssText(false) );
  
  // Always create the referenceLines object since we may need it for dynamic ref line colors
  Json::Object &refLines = base["referenceLines"] = Json::Value(Json::Type::Object);
  
  // Add dynamic reference line colors to the referenceLines section
  if( !info.dynamicRefLineMedicalColor.isDefault() )
    refLines["dynamicRefLineMedicalColor"] = WString( info.dynamicRefLineMedicalColor.cssText(false) );
  if( !info.dynamicRefLineIndustrialColor.isDefault() )
    refLines["dynamicRefLineIndustrialColor"] = WString( info.dynamicRefLineIndustrialColor.cssText(false) );
  if( !info.dynamicRefLineNormColor.isDefault() )
    refLines["dynamicRefLineNormColor"] = WString( info.dynamicRefLineNormColor.cssText(false) );
  if( !info.dynamicRefLineSnmColor.isDefault() )
    refLines["dynamicRefLineSnmColor"] = WString( info.dynamicRefLineSnmColor.cssText(false) );
  if( !info.dynamicRefLineCommonColor.isDefault() )
    refLines["dynamicRefLineCommonColor"] = WString( info.dynamicRefLineCommonColor.cssText(false) );
  if( !info.dynamicRefLineOtherColor.isDefault() )
    refLines["dynamicRefLineOtherColor"] = WString( info.dynamicRefLineOtherColor.cssText(false) );
  
  if( info.referenceLineColor.size() )
  {
    Json::Array &lineColors = refLines["lineColors"] = Json::Value(Json::Type::Array);
    for( const auto &c : info.referenceLineColor )
      lineColors.push_back( WString(c.cssText(false)) );
  }
  
  if( info.referenceLineColorForSources.size() )
  {
    Json::Object &srcs = refLines["specificSources"] = Json::Value(Json::Type::Object);
    for( const auto &p : info.referenceLineColorForSources )
      srcs[p.first] = WString( p.second.cssText(true) );
  }
  
  return Json::serialize( base );
}//string toJson( const ColorThemeInfo &info )


void ColorTheme::fromJson( const std::string &json, ColorTheme &info )
{
  info = ColorTheme();
  
  Json::Value baseValue;
  try
  {
    Json::parse( json, baseValue );  //Throws ParseError on failure
  }catch( Wt::Json::ParseError &e )
  {
    throw runtime_error( "Invalid JSON: " + string(e.what()) );
  }
  
  Json::Object &base = baseValue;
  
  auto &nameval = base.get("name");
  if( nameval.type() != Json::Type::String )//required!
    throw runtime_error( "JSON didnt contain string required property 'name'" );
  
  info.theme_name = nameval;
  if( base.contains("description") )
    info.theme_description = base.get("description");
  
  if( base.contains("created") )
  {
    const WString &createdstr = static_cast<const WString &>( base["created"] );
    //For some reason
    //info.creation_time = WDateTime::fromString( createdstr, "yyyy-MM-ddThh:mm:ssZ" );
    const SpecUtils::time_point_t created_time = SpecUtils::time_from_string( createdstr.toUTF8() );
    info.creation_time = WDateTime::fromTime_t( chrono::system_clock::to_time_t( created_time ) );
  }
  if( base.contains("modified") )
  {
    const WString &modifiedstr = static_cast<const WString &>( base["modified"] );
    //info.modified_time = WDateTime::fromString( modifiedstr, "yyyy-MM-ddThh:mm:ssZ" );
    const SpecUtils::time_point_t mod_time = SpecUtils::time_from_string( modifiedstr.toUTF8() );
    info.modified_time = WDateTime::fromTime_t( chrono::system_clock::to_time_t( mod_time ) );
  }
  
  
  info.nonChartAreaTheme = "";
  info.appColors.clear();
  const Json::Value &nonChartArea = base.get("nonChartArea");
  if( nonChartArea.type()==Json::Type::Object
      && static_cast<const Json::Object &>(nonChartArea).contains("cssTheme") )
  {
    const Json::Object &nonChartAreaObj = nonChartArea;
    const Json::Value &cssThemeVal = nonChartAreaObj.get("cssTheme");

    string val;
    if( cssThemeVal.type() == Wt::Json::Type::String )
    {
      if( cssThemeVal.hasType( typeid(Wt::WString) ) )
       val = static_cast<const WString &>(cssThemeVal).toUTF8();
      else if( cssThemeVal.hasType(typeid(std::string)) )
        val = static_cast<const std::string &>(cssThemeVal);
      else
        cout << "Json type is not string" << endl;
    }else
      cout << "CssThemeVal is type " << static_cast<int>(cssThemeVal.type()) << endl;

    // "default" is the always-loaded base theme, which is represented by an empty name
    if( SpecUtils::iequals_ascii(val, "default") )
      val = "";

    info.nonChartAreaTheme = val;

    // Only the colours the JSON lists become overrides; the rest keep the CSS theme's values.
    //  Older versions of InterSpec wrote camelCase keys, so those are accepted too.
    for( const AppColorToken &token : appColorTokens() )
    {
      const char *key = nullptr;
      if( nonChartAreaObj.contains(token.name) )
        key = token.name;
      else if( token.legacyJsonKey && nonChartAreaObj.contains(token.legacyJsonKey) )
        key = token.legacyJsonKey;

      if( !key )
        continue;

      const string colorStr = static_cast<const WString &>(nonChartAreaObj.get(key)).toUTF8();
      if( !colorStr.empty() )
        info.appColors[token.name] = WColor( colorStr );
    }//for( const AppColorToken &token : appColorTokens() )
  }//if( nonChartArea object with a cssTheme )
  
  if( base.contains("spectrum") /*&& base["spectrum"].type()==Json::Type::Object*/ )
  {
    Json::Object &spectrum = base["spectrum"];
    if( spectrum.contains("foregroundColor") )
      info.foregroundLine = WColor( static_cast<const WString &>(spectrum["foregroundColor"]) );
    if( spectrum.contains("backgroundColor") )
      info.backgroundLine = WColor( static_cast<const WString &>(spectrum["backgroundColor"]) );
    if( spectrum.contains("secondaryColor") )
      info.secondaryLine = WColor( static_cast<const WString &>(spectrum["secondaryColor"]) );
    if( spectrum.contains("axisLines") )
      info.spectrumAxisLines = WColor( static_cast<const WString &>(spectrum["axisLines"]) );
    if( spectrum.contains("chartBackground") )
      info.spectrumChartBackground = WColor( static_cast<const WString &>(spectrum["chartBackground"]) );
    if( spectrum.contains("chartMargins") )
      info.spectrumChartMargins = WColor( static_cast<const WString &>(spectrum["chartMargins"]) );
    if( spectrum.contains("chartText") )
      info.spectrumChartText = WColor( static_cast<const WString &>(spectrum["chartText"]) );
    if( spectrum.contains("peakLabelSize") )
      info.spectrumPeakLabelSize = spectrum["peakLabelSize"];
    if( spectrum.contains("peakLabelRotation") )
      info.spectrumPeakLabelRotation = spectrum["peakLabelRotation"];
    if( spectrum.contains("logYAxisMin") )
    {
      info.spectrumLogYAxisMin = spectrum["logYAxisMin"];
      if( info.spectrumLogYAxisMin <= 0.0 || IsNan(info.spectrumLogYAxisMin) || IsInf(info.spectrumLogYAxisMin) )
        info.spectrumLogYAxisMin = 0.1;
    }//
  }//if( spectrum node )
  
  if( base.contains("timeChart") )
  {
    Json::Object &timechart = base["timeChart"];
    if( timechart.contains("foregroundHighlightColor") )
      info.timeHistoryForegroundHighlight = WColor( static_cast<const WString &>(timechart["foregroundHighlightColor"]) );
    if( timechart.contains("backgroundHighlightColor") )
      info.timeHistoryBackgroundHighlight = WColor( static_cast<const WString &>(timechart["backgroundHighlightColor"]) );
    if( timechart.contains("secondaryHighlightColor") )
      info.timeHistorySecondaryHighlight = WColor( static_cast<const WString &>(timechart["secondaryHighlightColor"]) );
    if( timechart.contains("gammaLine") )
      info.timeChartGammaLine = WColor( static_cast<const WString &>(timechart["gammaLine"]) );
    if( timechart.contains("neutronLine") )
      info.timeChartNeutronLine = WColor( static_cast<const WString &>(timechart["neutronLine"]) );
    if( timechart.contains("axisLines") )
      info.timeAxisLines = WColor( static_cast<const WString &>(timechart["axisLines"]) );
    if( timechart.contains("chartBackground") )
      info.timeChartBackground = WColor( static_cast<const WString &>(timechart["chartBackground"]) );
    if( timechart.contains("chartMargins") )
      info.timeChartMargins = WColor( static_cast<const WString &>(timechart["chartMargins"]) );
    if( timechart.contains("chartText") )
      info.timeChartText = WColor( static_cast<const WString &>(timechart["chartText"]) );
    if( timechart.contains("occIndicatorLines") )
      info.occupancyIndicatorLines = WColor( static_cast<const WString &>(timechart["occIndicatorLines"]) );
  }//if( timeChart node )
  
  
  if( base.contains("peaksTakeOnReferenceLineColor") )
    info.peaksTakeOnReferenceLineColor = base["peaksTakeOnReferenceLineColor"];
  
  if( base.contains("defaultPeakLineColor") )
    info.defaultPeakLine = WColor( static_cast<const WString &>(base["defaultPeakLineColor"]) );
  
  // Set default values for dynamic reference line colors
  info.dynamicRefLineMedicalColor = WColor( sm_dynamic_ref_line_medical_color );
  info.dynamicRefLineIndustrialColor = WColor( sm_dynamic_ref_line_industrial_color );
  info.dynamicRefLineNormColor = WColor( sm_dynamic_ref_line_norm_color );
  info.dynamicRefLineSnmColor = WColor( sm_dynamic_ref_line_snm_color );
  info.dynamicRefLineCommonColor = WColor( sm_dynamic_ref_line_common_color );
  info.dynamicRefLineOtherColor = WColor( sm_dynamic_ref_line_other_color );
  
  if( base.contains("referenceLines") )
  {
    Json::Object &refLines = base["referenceLines"];
    
    // Check for dynamic reference line colors
    if( refLines.contains("dynamicRefLineMedicalColor") )
      info.dynamicRefLineMedicalColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineMedicalColor"]) );
    if( refLines.contains("dynamicRefLineIndustrialColor") )
      info.dynamicRefLineIndustrialColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineIndustrialColor"]) );
    if( refLines.contains("dynamicRefLineNormColor") )
      info.dynamicRefLineNormColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineNormColor"]) );
    if( refLines.contains("dynamicRefLineSnmColor") )
      info.dynamicRefLineSnmColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineSnmColor"]) );
    if( refLines.contains("dynamicRefLineCommonColor") )
      info.dynamicRefLineCommonColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineCommonColor"]) );
    if( refLines.contains("dynamicRefLineOtherColor") )
      info.dynamicRefLineOtherColor = WColor( static_cast<const WString &>(refLines["dynamicRefLineOtherColor"]) );
    
    if( refLines.contains("lineColors") )
    {
      info.referenceLineColor.clear();
      const Json::Array &lineColors = refLines["lineColors"];
      
      for( const auto &c : lineColors )
        info.referenceLineColor.push_back( WColor(static_cast<const WString &>(c)) );
    }
    
    if( refLines.contains("specificSources") )
    {
      info.referenceLineColorForSources.clear();
      const Json::Object &srcs = refLines["specificSources"];
      for( const auto &p : srcs )
        info.referenceLineColorForSources[p.first] = WColor(static_cast<const WString &>(p.second));
    }
  }//if( referenceLines node )
}//void fromJson( const std::string &json, ColorThemeInfo &info )


void ColorTheme::setFromDataBase( const ColorThemeInfo &db )
{
  ColorTheme::fromJson( db.json_data, *this );
  
  //this->dbIndex = [this isnt correct] db.user.id();
  
  //Should check if values from JSON match from ColorThemeInfo, but whatever for now.
  this->theme_name = db.theme_name;
  this->theme_description = db.theme_description;
  this->creation_time = db.creation_time;
  this->modified_time = db.modified_time;
}//void setFromDataBase( const ColorThemeInfo &db )
  

void ColorTheme::setToDataBase( ColorThemeInfo &db ) const
{
  db.json_data = ColorTheme::toJson( *this );
  db.theme_name = this->theme_name;
  db.theme_description = this->theme_description;
  db.creation_time = this->creation_time;
  db.modified_time = this->modified_time;
}//void setToDataBase( ColorThemeInfo &db ) const
