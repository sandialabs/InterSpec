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
#include <chrono>
#include <vector>
#include <string>
#include <typeinfo>

#include <Wt/WColor>
#include <Wt/WString>
#include <Wt/WDateTime>
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/Json/Serializer>

#include "InterSpec/ColorTheme.h"
#include "InterSpec/InterSpecUser.h"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"

using namespace std;
using namespace Wt;

const char * const ColorTheme::sm_color_theme_json_version = "1";
const char * const ColorTheme::sm_kinetic_ref_line_default_color = "#FFA500"; //Orange

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
}//namespace


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
        "kineticRefLineDefaultColor" : "#ff8c00",
        "modified" : "2018-12-26T02:58:03+0000",
        "name" : "Dark",
        "nonChartArea": { "cssTheme" : "dark" },
        "peaksTakeOnReferenceLineColor" : true,
        "referenceLines" : {
          "lineColors" : ["#c0c0c0", "#ffff99", "#b8d9f0", "#9933FF", "#FF66FF", "#CC3333", "#FF6633", "#FFFF99", "#CCFFCC", "#0000CC", "#666666", "#003333"],
          "specificSources" : {
            "Ba133" : "rgb(192,192,192)",
            "Th232" : "#4abfb9",
            "U235" : "rgb(128,0,128)",
            "background" : "#967f55"
          }
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
  
  foregroundLine = Wt::WColor(0x00, 0x00, 0x00); //Wt::GlobalColor::black
  backgroundLine = Wt::WColor(0x00, 0xff, 0xff); //Wt::GlobalColor::cyan
  secondaryLine = Wt::WColor(0x00, 0x80, 0x00);  //Wt::GlobalColor::darkGreen
  
  timeHistoryForegroundHighlight = Wt::WColor( 255, 255, 0, 155 );
  timeHistoryBackgroundHighlight = Wt::WColor( 0, 255, 255, 75 );
  timeHistorySecondaryHighlight = Wt::WColor( 0, 128, 0, 75 );
  
  defaultPeakLine = WColor( 0, 51, 255, 155 );
  
  spectrumAxisLines = WColor( GlobalColor::black );
  spectrumChartBackground = WColor();
  spectrumChartMargins = WColor();
  spectrumChartText = WColor( GlobalColor::black );
  
  spectrumPeakLabelSize = WString();
  spectrumPeakLabelRotation = 0.0;
  spectrumLogYAxisMin = 0.1;
  
  timeChartGammaLine = WColor( GlobalColor::black );
  timeChartNeutronLine = Wt::WColor( GlobalColor::darkGreen );
  
  timeAxisLines = WColor( GlobalColor::black );
  timeChartBackground = WColor();
  timeChartMargins = WColor();
  timeChartText = WColor( GlobalColor::black );
  
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
  
  referenceLineColorForSources["U235"] = WColor( "#800080" );
  referenceLineColorForSources["background"] = WColor( "#967f55" );  //brownish
  
  kineticRefLineDefaultColor = WColor( sm_kinetic_ref_line_default_color );
}//ColorTheme() constructor


std::string ColorTheme::toJson( const ColorTheme &info )
{
  Json::Object base;
  base["JsonVersion"] = WString(sm_color_theme_json_version);
  base["name"] = info.theme_name;
  base["description"] = info.theme_description;
  base["created"] = info.creation_time.toString( "yyyy-MM-ddThh:mm:ssZ" ); //"2018-11-01T22:56:00+0000"
  base["modified"] = info.modified_time.toString( "yyyy-MM-ddThh:mm:ssZ" );
  
  
  Json::Object &nonChartArea = base["nonChartArea"] = Json::Value(Json::ObjectType);
  if( info.nonChartAreaTheme.empty() )
    nonChartArea["cssTheme"] = WString("default");
  else
    nonChartArea["cssTheme"] = WString( info.nonChartAreaTheme );
  
  Json::Object &spectrum = base["spectrum"] = Json::Value(Json::ObjectType);
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
  
  Json::Object &timechart = base["timeChart"] = Json::Value(Json::ObjectType);
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
  
  if( !info.kineticRefLineDefaultColor.isDefault() )
    base["kineticRefLineDefaultColor"] = WString( info.kineticRefLineDefaultColor.cssText(false) );
  
  if( info.referenceLineColor.empty() && info.referenceLineColorForSources.empty() )
    return Json::serialize( base );
  
  Json::Object &refLines = base["referenceLines"] = Json::Value(Json::ObjectType);
  
  if( info.referenceLineColor.size() )
  {
    Json::Array &lineColors = refLines["lineColors"] = Json::Value(Json::ArrayType);
    for( const auto &c : info.referenceLineColor )
      lineColors.push_back( WString(c.cssText(false)) );
  }
  
  if( info.referenceLineColorForSources.size() )
  {
    Json::Object &srcs = refLines["specificSources"] = Json::Value(Json::ObjectType);
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
  if( nameval.type() != Json::Type::StringType )//required!
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
  const Json::Value &nonChartArea = base.get("nonChartArea");
  if( nonChartArea.type()==Json::ObjectType
      && static_cast<const Json::Object &>(nonChartArea).contains("cssTheme") )
  {
    const Json::Object &nonChartAreaObj = nonChartArea;
    const Json::Value &cssThemeVal = nonChartAreaObj.get("cssTheme");
    
    string val;
    if( cssThemeVal.type() == Wt::Json::StringType )
    {
      if( cssThemeVal.hasType( typeid(Wt::WString) ) )
       val = static_cast<const WString &>(cssThemeVal).toUTF8();
      else if( cssThemeVal.hasType(typeid(std::string)) )
        val = static_cast<const std::string &>(cssThemeVal);
      else
        cout << "Json type is not string" << endl;
    }else
      cout << "CssThemeVal is type " << cssThemeVal.type() << endl;
    
    if( SpecUtils::iequals_ascii(info.nonChartAreaTheme, "default") )
      val = "";
    
    info.nonChartAreaTheme = val;
  }
  
  if( base.contains("spectrum") /*&& base["spectrum"].type()==Json::ObjectType*/ )
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
  
  if( base.contains("kineticRefLineDefaultColor") )
    info.kineticRefLineDefaultColor = WColor( static_cast<const WString &>(base["kineticRefLineDefaultColor"]) );
  else
    info.kineticRefLineDefaultColor = WColor( sm_kinetic_ref_line_default_color );  // Default color if not in JSON
  
  if( base.contains("referenceLines") )
  {
    Json::Object &refLines = base["referenceLines"];
    
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
