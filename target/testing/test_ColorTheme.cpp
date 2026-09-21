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

// Must be defined before Windows.h (or any header that includes it) is included; see CLAUDE.md.
#ifdef _WIN32
  #define WIN32_LEAN_AND_MEAN
  #include <winsock2.h>
  #include <windows.h>
#endif

#include "InterSpec_config.h"

#include <map>
#include <set>
#include <regex>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#define BOOST_TEST_MODULE TestColorTheme
#include <boost/test/included/unit_test.hpp>

#include <Wt/WColor.h>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/ColorTheme.h"

using namespace std;
using namespace boost::unit_test;


/** Keeps the `--interspec-*` colour-token system consistent.

 The tokens are declared in CSS (InterSpec_resources/themes/default/default.css is the light
 palette, themes/<name>/<name>.css re-declares them for other themes), consumed by every other
 stylesheet with no fallback, and published as user overrides from the C++ `ColorTheme` table.
 Nothing at build time otherwise ties those three together: a token used in a stylesheet but
 declared in neither theme silently falls back to whatever the browser inherits, and a token missing
 from one theme silently keeps the other theme's value.  These checks catch that drift.
 */

namespace
{
  string g_source_dir;

  /** Stylesheets that are also used outside the app (embedded in standalone HTML reports, or
   shared with SpecUtils), where every `var()` keeps a fallback normalized to the light value.
   */
  const set<string> sm_keep_fallback_files{
    "SpectrumChartD3.css", "ShieldingSourceFitPlot.css", "RelEffPlot.css", "DrfChart.css",
    "DecayChainChart.css"
  };

  /** Tokens declared in CSS that are theme "mode" rather than colours; they are deliberately not
   user-overridable, so they are not in the ColorTheme table.
   */
  const set<string> sm_css_only_tokens{ "icon-filter", "color-scheme" };


  /** Derives the source tree from `--datadir` (the `data` directory sits at its root). */
  void set_source_dir()
  {
    if( !g_source_dir.empty() )
      return;

    const int argc = framework::master_test_suite().argc;
    const char * const * const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }
    SpecUtils::ireplace_all( datadir, "%20", " " );

    if( !datadir.empty() )
      g_source_dir = SpecUtils::parent_path( datadir );

    if( g_source_dir.empty() )
    {
      for( const char * const d : { ".", "..", "../..", "../../..", "../../../.." } )
      {
        if( SpecUtils::is_directory( SpecUtils::append_path( d, "InterSpec_resources" ) ) )
        {
          g_source_dir = d;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( !g_source_dir.empty(), "Could not find the source tree; pass --datadir=" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( SpecUtils::append_path( g_source_dir, "InterSpec_resources" ) ),
                          "'" << g_source_dir << "' has no InterSpec_resources directory" );
  }//void set_source_dir()


  string resources_dir()
  {
    set_source_dir();
    return SpecUtils::append_path( g_source_dir, "InterSpec_resources" );
  }


  string file_contents( const string &path )
  {
    ifstream input( path.c_str(), ios::in | ios::binary );
    BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open '" << path << "'" );
    stringstream strm;
    strm << input.rdbuf();
    return strm.str();
  }


  string strip_comments( string css )
  {
    return regex_replace( css, regex("/\\*[\\s\\S]*?\\*/"), "" );
  }


  /** The `--interspec-<name>: <value>;` declarations of one stylesheet. */
  map<string,string> declared_tokens( const string &path )
  {
    map<string,string> tokens;
    const string css = strip_comments( file_contents(path) );
    const regex declaration( "--interspec-([a-z0-9-]+)\\s*:\\s*([^;{}]+);" );
    const sregex_iterator end_it;
    for( sregex_iterator it( begin(css), end(css), declaration ); it != end_it; ++it )
    {
      string value = (*it)[2].str();
      SpecUtils::trim( value );
      tokens[(*it)[1].str()] = value;
    }
    return tokens;
  }


  string normalize_value( string value )
  {
    SpecUtils::to_lower_ascii( value );
    SpecUtils::erase_any_character( value, " \t\r\n" );
    return value;
  }


  /** Every `var(--interspec-NAME[, FALLBACK])` use in `text`, as (name, fallback) pairs; fallback
   is empty when there is none.
   */
  vector<pair<string,string>> token_uses( const string &text )
  {
    vector<pair<string,string>> uses;
    const string opener = "var(--interspec-";
    size_t pos = text.find( opener );
    while( pos != string::npos )
    {
      const size_t paren = pos + 3;  // the "(" of "var("
      int depth = 0;
      size_t close = string::npos;
      for( size_t i = paren; i < text.size(); ++i )
      {
        if( text[i] == '(' )
          depth += 1;
        else if( text[i] == ')' )
        {
          depth -= 1;
          if( depth == 0 )
          {
            close = i;
            break;
          }
        }
      }
      if( close == string::npos )
        break;

      const string inner = text.substr( paren + 1, close - paren - 1 );  // "--interspec-x, fallback"
      const size_t comma = inner.find( ',' );
      string name = inner.substr( string("--interspec-").size(),
                                  (comma == string::npos) ? string::npos : comma - string("--interspec-").size() );
      string fallback = (comma == string::npos) ? string() : inner.substr( comma + 1 );
      SpecUtils::trim( name );
      SpecUtils::trim( fallback );
      uses.emplace_back( name, fallback );

      pos = text.find( opener, pos + opener.size() );
    }
    return uses;
  }


  /** All first-party stylesheets (third-party ones under assets/ excluded). */
  vector<string> first_party_css()
  {
    vector<string> files;
    for( const string &f : SpecUtils::recursive_ls( resources_dir(), ".css" ) )
    {
      if( f.find( "/assets/" ) == string::npos && f.find( "\\assets\\" ) == string::npos )
        files.push_back( f );
    }
    BOOST_REQUIRE_MESSAGE( files.size() > 50, "Found only " << files.size() << " stylesheets" );
    return files;
  }
}//namespace


BOOST_AUTO_TEST_CASE( TokenBlocksAgree )
{
  const string themes = SpecUtils::append_path( resources_dir(), "themes" );
  const map<string,string> light = declared_tokens( SpecUtils::append_path( themes, "default/default.css" ) );
  const map<string,string> dark = declared_tokens( SpecUtils::append_path( themes, "dark/dark.css" ) );

  BOOST_REQUIRE( light.size() > 30 );

  for( const auto &nameAndValue : light )
  {
    BOOST_CHECK_MESSAGE( !nameAndValue.second.empty(), "--interspec-" << nameAndValue.first << " is empty in default.css" );
    BOOST_CHECK_MESSAGE( dark.count(nameAndValue.first), "--interspec-" << nameAndValue.first << " is declared in default.css but not dark.css" );
  }
  for( const auto &nameAndValue : dark )
  {
    BOOST_CHECK_MESSAGE( !nameAndValue.second.empty(), "--interspec-" << nameAndValue.first << " is empty in dark.css" );
    BOOST_CHECK_MESSAGE( light.count(nameAndValue.first), "--interspec-" << nameAndValue.first << " is declared in dark.css but not default.css" );
  }

  // Every other CSS theme on disk must carry the full list too
  for( const string &theme : ColorTheme::availableCssThemes( resources_dir() ) )
  {
    if( theme == "default" )
      continue;
    const map<string,string> other = declared_tokens( SpecUtils::append_path( SpecUtils::append_path(themes, theme), theme + ".css" ) );
    for( const auto &nameAndValue : light )
      BOOST_CHECK_MESSAGE( other.count(nameAndValue.first), "--interspec-" << nameAndValue.first << " missing from theme '" << theme << "'" );
  }
}//BOOST_AUTO_TEST_CASE( TokenBlocksAgree )


BOOST_AUTO_TEST_CASE( TableMatchesCss )
{
  const map<string,string> light = declared_tokens( SpecUtils::append_path( resources_dir(), "themes/default/default.css" ) );

  set<string> table_names;
  for( const ColorTheme::AppColorToken &token : ColorTheme::appColorTokens() )
  {
    BOOST_CHECK_MESSAGE( table_names.insert( token.name ).second, "Duplicate table entry " << token.name );
    BOOST_CHECK_MESSAGE( light.count( token.name ), "Table token '" << token.name << "' is not declared in default.css" );
    BOOST_CHECK( ColorTheme::appColorToken( token.name ) == &token );
  }

  for( const auto &nameAndValue : light )
  {
    BOOST_CHECK_MESSAGE( table_names.count(nameAndValue.first) || sm_css_only_tokens.count(nameAndValue.first),
                        "--interspec-" << nameAndValue.first << " is declared in CSS but is neither in the"
                        " ColorTheme table nor in the CSS-only list" );
  }

  BOOST_CHECK( !ColorTheme::appColorToken( "no-such-token" ) );
}//BOOST_AUTO_TEST_CASE( TableMatchesCss )


BOOST_AUTO_TEST_CASE( EveryUsedTokenIsDeclared )
{
  const map<string,string> light = declared_tokens( SpecUtils::append_path( resources_dir(), "themes/default/default.css" ) );

  vector<string> files = first_party_css();
  for( const string &f : SpecUtils::ls_files_in_directory( resources_dir(), ".js" ) )
    files.push_back( f );
  for( const string &f : SpecUtils::ls_files_in_directory( SpecUtils::append_path( g_source_dir, "src" ), ".cpp" ) )
    files.push_back( f );

  size_t total_uses = 0;
  for( const string &path : files )
  {
    const string text = file_contents( path );
    for( const pair<string,string> &use : token_uses( text ) )
    {
      total_uses += 1;
      BOOST_CHECK_MESSAGE( light.count( use.first ),
                          SpecUtils::filename(path) << " uses undeclared token --interspec-" << use.first );
    }
    BOOST_CHECK_MESSAGE( text.find("--interspec-warning-text-color") == string::npos,
                        SpecUtils::filename(path) << " still uses the retired --interspec-warning-text-color" );
  }

  BOOST_CHECK_MESSAGE( total_uses > 300, "Only found " << total_uses << " token uses; the scan is probably broken" );
}//BOOST_AUTO_TEST_CASE( EveryUsedTokenIsDeclared )


BOOST_AUTO_TEST_CASE( FallbackPolicy )
{
  const map<string,string> light = declared_tokens( SpecUtils::append_path( resources_dir(), "themes/default/default.css" ) );

  for( const string &path : first_party_css() )
  {
    const string name = SpecUtils::filename( path );
    const bool keep = sm_keep_fallback_files.count( name );
    const string css = strip_comments( file_contents( path ) );

    for( const pair<string,string> &use : token_uses( css ) )
    {
      if( !keep )
      {
        BOOST_CHECK_MESSAGE( use.second.empty(), name << ": in-app stylesheets carry no fallback, but found"
                            " var(--interspec-" << use.first << ", " << use.second << ")" );
        continue;
      }

      const map<string,string>::const_iterator pos = light.find( use.first );
      if( pos == end(light) )
        continue;  // reported by EveryUsedTokenIsDeclared
      BOOST_CHECK_MESSAGE( normalize_value(use.second) == normalize_value(pos->second),
                          name << ": fallback for --interspec-" << use.first << " is '" << use.second
                          << "' but the light value is '" << pos->second << "'" );
    }
  }
}//BOOST_AUTO_TEST_CASE( FallbackPolicy )


BOOST_AUTO_TEST_CASE( JsonRoundTrip )
{
  ColorTheme theme;
  theme.theme_name = "Round trip";
  theme.nonChartAreaTheme = "dark";
  theme.appColors["text-color"] = Wt::WColor( "white" );
  theme.appColors["background-color"] = Wt::WColor( 44, 45, 48 );
  theme.appColors["menubar-hover-color"] = Wt::WColor( 18, 101, 200, 64 );  // alpha must survive

  const string json = ColorTheme::toJson( theme );
  ColorTheme parsed;
  ColorTheme::fromJson( json, parsed );

  BOOST_CHECK_EQUAL( parsed.nonChartAreaTheme, "dark" );
  BOOST_REQUIRE_EQUAL( parsed.appColors.size(), theme.appColors.size() );
  for( const auto &nameAndColor : theme.appColors )
  {
    BOOST_REQUIRE_MESSAGE( parsed.appColors.count(nameAndColor.first), "Lost " << nameAndColor.first );
    BOOST_CHECK_EQUAL( parsed.appColors[nameAndColor.first].cssText(true), nameAndColor.second.cssText(true) );
  }
  BOOST_CHECK_EQUAL( parsed.appColors["menubar-hover-color"].alpha(), 64 );

  // A theme with no overrides writes only the CSS theme name, and reads back empty
  ColorTheme bare;
  bare.theme_name = "Bare";
  ColorTheme parsed_bare;
  ColorTheme::fromJson( ColorTheme::toJson( bare ), parsed_bare );
  BOOST_CHECK( parsed_bare.appColors.empty() );
  BOOST_CHECK( parsed_bare.nonChartAreaTheme.empty() );
}//BOOST_AUTO_TEST_CASE( JsonRoundTrip )


BOOST_AUTO_TEST_CASE( LegacyJsonLoads )
{
  // The `nonChartArea` block older versions wrote (camelCase keys), from the old Dark literal
  const string legacy = R"json({
    "JsonVersion" : "1",
    "name" : "Legacy dark",
    "nonChartArea": {
      "cssTheme" : "dark",
      "backgroundColor": "rgb(44,45,48)",
      "textColor": "white",
      "borderColor": "rgb(136,136,136)",
      "linkColor": "rgb(18,101,200)",
      "labelColor": "rgb(197,198,201)",
      "inputBackground": "rgb(59,60,63)",
      "buttonBackground": "rgb(86,88,90)",
      "buttonBorderColor": "rgb(106,107,109)",
      "buttonTextColor": "rgb(230,230,230)"
    }
  })json";

  ColorTheme theme;
  ColorTheme::fromJson( legacy, theme );
  BOOST_CHECK_EQUAL( theme.nonChartAreaTheme, "dark" );
  BOOST_CHECK_EQUAL( theme.appColors.size(), 9 );
  BOOST_CHECK_EQUAL( theme.appColors["background-color"].cssText(), "rgb(44,45,48)" );
  BOOST_CHECK_EQUAL( theme.appColors["text-color"].cssText(), "white" );
  BOOST_CHECK_EQUAL( theme.appColors["button-text-color"].cssText(), "rgb(230,230,230)" );
  BOOST_CHECK( !theme.appColors.count("menubar-background") );  // absent stays absent (no fallback ladder)

  // "default" is the base theme, represented by an empty name
  const string light = R"json({ "name" : "Light", "nonChartArea": { "cssTheme" : "default", "text-color": "#123456" } })json";
  ColorTheme lightTheme;
  ColorTheme::fromJson( light, lightTheme );
  BOOST_CHECK( lightTheme.nonChartAreaTheme.empty() );
  BOOST_CHECK_EQUAL( lightTheme.appColors.size(), 1 );
  BOOST_CHECK_EQUAL( lightTheme.appColors["text-color"].cssText(), "#123456" );

  // The predefined dark theme now leaves the palette to dark.css
  const unique_ptr<ColorTheme> dark = ColorTheme::predefinedTheme( ColorTheme::DarkColorTheme );
  BOOST_REQUIRE( dark );
  BOOST_CHECK_EQUAL( dark->nonChartAreaTheme, "dark" );
  BOOST_CHECK( dark->appColors.empty() );
}//BOOST_AUTO_TEST_CASE( LegacyJsonLoads )


BOOST_AUTO_TEST_CASE( FuzzCorpusParses )
{
  // The corpus files are fuzz seeds, not all of them valid themes (minimal.json has no "name"), so
  //  the requirement is only that parsing either succeeds or reports failure by exception.
  set_source_dir();
  const string corpus = SpecUtils::append_path( g_source_dir, "target/fuzzing/corpus/color_theme" );
  const vector<string> files = SpecUtils::ls_files_in_directory( corpus, ".json" );
  BOOST_REQUIRE_MESSAGE( !files.empty(), "No corpus files in " << corpus );
  for( const string &path : files )
  {
    ColorTheme theme;
    try
    {
      ColorTheme::fromJson( file_contents(path), theme );
      BOOST_TEST_MESSAGE( SpecUtils::filename(path) << " parsed" );
    }catch( const std::exception &e )
    {
      BOOST_TEST_MESSAGE( SpecUtils::filename(path) << " rejected: " << e.what() );
    }
  }

  ColorTheme theme;
  BOOST_CHECK_NO_THROW( ColorTheme::fromJson( file_contents( SpecUtils::append_path(corpus, "default_theme.json") ), theme ) );
}//BOOST_AUTO_TEST_CASE( FuzzCorpusParses )


BOOST_AUTO_TEST_CASE( CssThemeHelpers )
{
  const vector<string> themes = ColorTheme::availableCssThemes( resources_dir() );
  BOOST_REQUIRE( !themes.empty() );
  BOOST_CHECK_EQUAL( themes.front(), "default" );
  BOOST_CHECK( std::find( begin(themes), end(themes), "dark" ) != end(themes) );

  const map<string,string> light = ColorTheme::cssThemeTokenDefaults( resources_dir(), "" );
  const map<string,string> dark = ColorTheme::cssThemeTokenDefaults( resources_dir(), "dark" );
  BOOST_CHECK_EQUAL( light.at("text-color"), "black" );
  BOOST_CHECK_EQUAL( dark.at("text-color"), "white" );
  BOOST_CHECK_EQUAL( dark.at("color-scheme"), "dark" );
  BOOST_CHECK_EQUAL( light.size(), dark.size() );

  // An unknown theme name just yields the base values
  const map<string,string> unknown = ColorTheme::cssThemeTokenDefaults( resources_dir(), "no-such-theme" );
  BOOST_CHECK( unknown == light );
}//BOOST_AUTO_TEST_CASE( CssThemeHelpers )
