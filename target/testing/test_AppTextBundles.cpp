/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
 Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain
 rights in this software.
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

#include <set>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#define BOOST_TEST_MODULE TestAppTextBundles
#include <boost/test/included/unit_test.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

using namespace std;
using namespace boost::unit_test;


/** Structural checks over every localization bundle in `InterSpec_resources/app_text`.

 These files are data, so nothing else in the build looks at them: CMake only globs them for IDE
 visibility, and a broken one fails at *runtime* by being silently ignored - with no message at all
 on Windows (see `app_text/README.md`).  That combination is why 18 dead keys, three malformed
 `<message id="dls-">` stubs, and a drifting set of files survived unnoticed.

 The checks are deliberately structural rather than linguistic - nothing here can tell whether a
 translation is any good.  What it can tell is that a file will load, that an id a translator added
 still corresponds to something the C++ asks for, and that a translated string did not lose a
 `{1}` and so silently drop a number out of a sentence.
 */

namespace
{
  /** Message ids and the argument placeholders each one uses. */
  typedef map<string,set<string>> BundleContents;

  string g_app_text_dir;


  /** Locates `InterSpec_resources/app_text`.

   Accepts `--apptextdir=`, and otherwise derives it from `--datadir` (the two are siblings in the
   source tree) before falling back to a search upward from the working directory - so the test
   runs both from ctest and by hand, with the same arguments the other tests take.
   */
  void set_app_text_dir()
  {
    if( !g_app_text_dir.empty() )
      return;

    const int argc = framework::master_test_suite().argc;
    const char * const * const argv = framework::master_test_suite().argv;

    string datadir;
    for( int i = 1; i < argc; ++i )
    {
      const string arg = argv[i];
      if( SpecUtils::istarts_with( arg, "--apptextdir=" ) )
        g_app_text_dir = arg.substr( 13 );
      if( SpecUtils::istarts_with( arg, "--datadir=" ) )
        datadir = arg.substr( 10 );
    }

    SpecUtils::ireplace_all( g_app_text_dir, "%20", " " );
    SpecUtils::ireplace_all( datadir, "%20", " " );

    // `data` and `InterSpec_resources` sit beside each other, so --datadir already says where the
    //  source tree is; this keeps the test runnable with exactly the arguments the others take.
    if( g_app_text_dir.empty() && !datadir.empty() )
    {
      const string parent = SpecUtils::parent_path( datadir );
      const string candidate = SpecUtils::append_path(
                                 SpecUtils::append_path( parent, "InterSpec_resources" ),
                                 "app_text" );
      if( SpecUtils::is_directory( candidate ) )
        g_app_text_dir = candidate;
    }

    if( g_app_text_dir.empty() )
    {
      for( const char * const d : { ".", "..", "../..", "../../..", "../../../.." } )
      {
        const string candidate = SpecUtils::append_path(
                                   SpecUtils::append_path( d, "InterSpec_resources" ), "app_text" );
        if( SpecUtils::is_file( SpecUtils::append_path( candidate, "README.md" ) ) )
        {
          g_app_text_dir = candidate;
          break;
        }
      }
    }

    BOOST_REQUIRE_MESSAGE( !g_app_text_dir.empty(),
                          "Could not find InterSpec_resources/app_text; pass --apptextdir=" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_directory( g_app_text_dir ),
                          "'" << g_app_text_dir << "' is not a directory" );
  }//void set_app_text_dir()


  /** The named character entities Wt will accept.

   Wt parses these files with its own rapidxml variant, whose `translate_xhtml_entity` table is the
   whole vocabulary available (Wt 3.7.1, `src/3rdparty/rapidxml/rapidxml_xhtml.hpp`; the version is
   pinned, see CLAUDE.md).  Anything outside it is a parse error that discards the entire file at
   startup, which is exactly the failure that produces no message on Windows - so it is worth
   catching here rather than in a bug report.  Numeric references (`&#160;` / `&#xA0;`) are always
   fine and are checked separately.
   */
  const set<string> &allowed_entities()
  {
    static const set<string> s_allowed = {
      "AElig", "Aacute", "Acirc", "Agrave", "Alpha", "Aring", "Atilde", "Auml", "Beta", "Ccedil",
      "Chi", "Dagger", "Delta", "Dstrok", "ETH", "Eacute", "Ecirc", "Egrave", "Epsilon", "Eta",
      "Euml", "Gamma", "Iacute", "Icirc", "Igrave", "Iota", "Iuml", "Kappa", "Lambda", "Mu",
      "Ntilde", "Nu", "OElig", "Oacute", "Ocirc", "Ograve", "Omega", "Omicron", "Oslash", "Otilde",
      "Ouml", "Phi", "Pi", "Prime", "Psi", "Rho", "Scaron", "Sigma", "THORN", "Tau", "Theta",
      "Uacute", "Ucirc", "Ugrave", "Upsilon", "Uuml", "Xi", "Yacute", "Yuml", "Zeta", "aacute",
      "acirc", "acute", "aelig", "agrave", "alefsym", "alpha", "amp", "and", "ang", "apos",
      "aring", "asymp", "atilde", "auml", "bdquo", "beta", "brkbar", "brvbar", "bull", "cap",
      "ccedil", "cedil", "cent", "chi", "circ", "clubs", "cong", "copy", "crarr", "cup", "curren",
      "dArr", "dagger", "darr", "deg", "delta", "diams", "die", "divide", "eacute", "ecirc",
      "egrave", "empty", "emsp", "ensp", "epsilon", "equiv", "eta", "eth", "euml", "euro", "exist",
      "fnof", "forall", "frac12", "frac14", "frac34", "frasl", "gamma", "ge", "gt", "hArr", "harr",
      "hearts", "hellip", "hibar", "iacute", "icirc", "iexcl", "igrave", "image", "infin", "int",
      "iota", "iquest", "isin", "iuml", "kappa", "lArr", "lambda", "lang", "laquo", "larr",
      "lceil", "ldquo", "le", "lfloor", "lowast", "loz", "lrm", "lsaquo", "lsquo", "lt", "macr",
      "mdash", "micro", "middot", "minus", "mu", "nabla", "nbsp", "ndash", "ne", "ni", "not",
      "notin", "nsub", "ntilde", "nu", "oacute", "ocirc", "oelig", "ograve", "oline", "omega",
      "omicron", "oplus", "or", "ordf", "ordm", "oslash", "otilde", "otimes", "ouml", "para",
      "part", "permil", "perp", "phi", "pi", "piv", "plusmn", "pound", "prime", "prod", "prop",
      "psi", "quot", "rArr", "radic", "rang", "raquo", "rarr", "rceil", "rdquo", "real", "reg",
      "rfloor", "rho", "rlm", "rsaquo", "rsquo", "sbquo", "scaron", "sdot", "sect", "shy", "sigma",
      "sigmaf", "sim", "spades", "sub", "sube", "sum", "sup", "sup1", "sup2", "sup3", "supe",
      "szlig", "tau", "there4", "theta", "thetasym", "thinsp", "thorn", "tilde", "times", "trade",
      "uArr", "uacute", "uarr", "ucirc", "ugrave", "uml", "upsih", "upsilon", "uuml", "weierp",
      "xi", "yacute", "yen", "yuml", "zeta", "zwj", "zwnj"
    };

    return s_allowed;
  }//const set<string> &allowed_entities()


  /** True if \p text is well-formed UTF-8.  A bundle saved in Latin-1 loads as mojibake rather
   than failing, so nothing downstream would report it.
   */
  bool is_valid_utf8( const string &text )
  {
    size_t i = 0;
    while( i < text.size() )
    {
      const unsigned char lead = static_cast<unsigned char>( text[i] );
      size_t extra = 0;

      if( lead < 0x80 )
        extra = 0;
      else if( (lead & 0xE0) == 0xC0 )
        extra = 1;
      else if( (lead & 0xF0) == 0xE0 )
        extra = 2;
      else if( (lead & 0xF8) == 0xF0 )
        extra = 3;
      else
        return false;

      if( (i + extra) >= text.size() )
        return false;

      for( size_t j = 1; j <= extra; ++j )
      {
        if( (static_cast<unsigned char>( text[i + j] ) & 0xC0) != 0x80 )
          return false;
      }

      i += 1 + extra;
    }//while( i < text.size() )

    return true;
  }//bool is_valid_utf8( const string &text )


  /** The `{1}`, `{2}`, ... argument slots a message uses, as a set.

   A set rather than a count: `WString::arg()` fills slots by number, so a translation that keeps
   three slots but renumbers one still drops a value.  Order does not matter - languages reorder
   the sentence, which is the point of numbered slots.
   */
  set<string> placeholders_in( const string &text )
  {
    set<string> answer;
    for( size_t open = text.find( '{' ); open != string::npos; open = text.find( '{', open + 1 ) )
    {
      const size_t close = text.find( '}', open + 1 );
      if( close == string::npos )
        break;

      const string inner = text.substr( open + 1, close - open - 1 );
      if( !inner.empty()
         && (inner.find_first_not_of( "0123456789" ) == string::npos) )
      {
        answer.insert( inner );
      }
    }

    return answer;
  }//set<string> placeholders_in( const string &text )


  /** Named entity references used in \p text that Wt would reject. */
  set<string> disallowed_entities_in( const string &text )
  {
    set<string> answer;
    for( size_t amp = text.find( '&' ); amp != string::npos; amp = text.find( '&', amp + 1 ) )
    {
      const size_t semi = text.find( ';', amp + 1 );
      if( (semi == string::npos) || ((semi - amp) > 12) )
        continue;

      const string name = text.substr( amp + 1, semi - amp - 1 );
      if( name.empty() )
        continue;

      if( name[0] == '#' )   // numeric reference; always valid
        continue;

      if( name.find_first_not_of( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789" )
         != string::npos )
      {
        continue;   // not an entity reference at all, just a stray '&'
      }

      if( !allowed_entities().count( name ) )
        answer.insert( name );
    }

    return answer;
  }//set<string> disallowed_entities_in( const string &text )


  /** Concatenated text of a node and everything under it, entities left as written. */
  string node_text( const rapidxml::xml_node<char> * const node )
  {
    string answer;
    if( !node )
      return answer;

    if( node->type() == rapidxml::node_data || node->type() == rapidxml::node_cdata )
      answer.append( node->value(), node->value_size() );

    for( const rapidxml::xml_node<char> *child = node->first_node();
        child; child = child->next_sibling() )
    {
      answer += node_text( child );
    }

    return answer;
  }//string node_text( const rapidxml::xml_node<char> *node )


  /** Parses one bundle, checking everything that can be judged from the file alone.

   Returns the ids it defines and the placeholders each uses; `BOOST_ERROR`s and returns empty if
   the file could not be parsed at all.
   */
  BundleContents check_one_bundle( const string &path )
  {
    BundleContents answer;

    vector<char> data;
    {
      ifstream input( path.c_str(), ios::in | ios::binary );
      BOOST_REQUIRE_MESSAGE( input.is_open(), "Failed to open '" << path << "'" );
      data.assign( istreambuf_iterator<char>( input ), istreambuf_iterator<char>() );
    }

    BOOST_REQUIRE_MESSAGE( !data.empty(), "'" << path << "' is empty" );

    const string as_text( data.begin(), data.end() );
    BOOST_CHECK_MESSAGE( is_valid_utf8( as_text ), "'" << path << "' is not valid UTF-8" );

    data.push_back( '\0' );

    rapidxml::xml_document<char> document;
    try
    {
      // Entity translation is off on purpose: the entity vocabulary is Wt's, not rapidxml's, so
      //  the set is checked explicitly below rather than by whichever parser happens to run here.
      document.parse<rapidxml::parse_no_entity_translation>( &data.front() );
    }catch( std::exception &e )
    {
      BOOST_ERROR( "'" << path << "' is not well-formed XML: " << e.what() );
      return answer;
    }

    const rapidxml::xml_node<char> * const root = document.first_node( "messages", 8 );
    if( !root )
    {
      BOOST_ERROR( "'" << path << "' has no root <messages> element" );
      return answer;
    }

    for( const rapidxml::xml_node<char> *message = root->first_node( "message", 7 );
        message; message = message->next_sibling( "message", 7 ) )
    {
      const rapidxml::xml_attribute<char> * const id_attrib
                                                 = message->first_attribute( "id", 2 );
      if( !id_attrib || (id_attrib->value_size() == 0) )
      {
        BOOST_ERROR( "'" << path << "' has a <message> with no id" );
        continue;
      }

      const string id( id_attrib->value(), id_attrib->value_size() );

      // An id that is only the file's prefix is the shape the malformed stubs took; they are
      //  invisible at runtime (nothing ever looks them up) and so survive indefinitely.
      if( (id.size() < 2) || (id.back() == '-') )
        BOOST_ERROR( "'" << path << "' has the placeholder-looking id '" << id << "'" );

      if( answer.count( id ) )
        BOOST_ERROR( "'" << path << "' defines id '" << id << "' more than once" );

      const string text = node_text( message );

      const set<string> bad = disallowed_entities_in( text );
      for( const string &entity : bad )
      {
        BOOST_ERROR( "'" << path << "' message '" << id << "' uses entity '&" << entity
                    << ";', which Wt cannot parse - use the UTF-8 character instead" );
      }

      answer[id] = placeholders_in( text );
    }//for( loop over <message> elements )

    BOOST_CHECK_MESSAGE( !answer.empty(), "'" << path << "' defines no messages" );

    return answer;
  }//BundleContents check_one_bundle( const string &path )


  /** Splits "DetectionLimitTool_es.xml" into ("DetectionLimitTool", "es").

   Locale suffixes here are always a short alphabetic tag; base names legitimately contain
   underscores, so the split is on the *last* one, and only when what follows looks like a locale.
   */
  bool split_bundle_name( const string &filename, string &base, string &locale )
  {
    string stem = filename;
    if( !SpecUtils::iends_with( stem, ".xml" ) )
      return false;
    stem = stem.substr( 0, stem.size() - 4 );

    const size_t underscore = stem.rfind( '_' );
    if( underscore == string::npos )
    {
      base = stem;
      locale.clear();
      return true;
    }

    const string suffix = stem.substr( underscore + 1 );
    const bool looks_like_locale
        = !suffix.empty() && (suffix.size() <= 5)
          && (suffix.find_first_not_of( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-" )
              == string::npos);

    if( !looks_like_locale )
    {
      base = stem;
      locale.clear();
      return true;
    }

    base = stem.substr( 0, underscore );
    locale = suffix;
    return true;
  }//bool split_bundle_name(...)
}//namespace


BOOST_AUTO_TEST_CASE( EveryBundleIsWellFormed )
{
  set_app_text_dir();

  const vector<string> files
      = SpecUtils::recursive_ls( g_app_text_dir, ".xml" );

  BOOST_REQUIRE_MESSAGE( files.size() > 50,
                        "Only found " << files.size() << " bundles in '" << g_app_text_dir
                        << "' - the directory is probably wrong" );

  size_t num_messages = 0;
  for( const string &path : files )
    num_messages += check_one_bundle( path ).size();

  BOOST_TEST_MESSAGE( "Checked " << files.size() << " bundles, " << num_messages << " messages" );
}//BOOST_AUTO_TEST_CASE( EveryBundleIsWellFormed )


BOOST_AUTO_TEST_CASE( LocaleBundlesTrackTheirEnglishBase )
{
  set_app_text_dir();

  const vector<string> files = SpecUtils::recursive_ls( g_app_text_dir, ".xml" );
  BOOST_REQUIRE( !files.empty() );

  // base name -> locale ("" for the English file) -> contents
  map<string,map<string,BundleContents>> bundles;
  for( const string &path : files )
  {
    string base, locale;
    if( !split_bundle_name( SpecUtils::filename( path ), base, locale ) )
      continue;

    bundles[base][locale] = check_one_bundle( path );
  }

  size_t num_locale_files = 0;
  for( const auto &base_entry : bundles )
  {
    const string &base = base_entry.first;
    const map<string,BundleContents> &by_locale = base_entry.second;

    const map<string,BundleContents>::const_iterator english = by_locale.find( string() );
    if( english == by_locale.end() )
    {
      // A locale file with no English base means either the base was deleted or the split above
      //  mistook part of the name for a locale; both are worth saying out loud.
      BOOST_ERROR( "Bundle '" << base << "' has locale files but no English '" << base
                  << ".xml'" );
      continue;
    }

    for( const auto &locale_entry : by_locale )
    {
      if( locale_entry.first.empty() )
        continue;

      ++num_locale_files;

      const string what = base + "_" + locale_entry.first + ".xml";

      for( const auto &message : locale_entry.second )
      {
        const BundleContents::const_iterator base_message
                                             = english->second.find( message.first );

        // Subset, not equality.  Partial translation is the documented and intended steady state
        //  (`app_text/README.md`): a missing key falls back to English, and many bundles ship that
        //  way on purpose.  An *extra* key is the actual defect - it is either a typo, so the
        //  translation never appears, or a leftover of a key the C++ stopped asking for.
        if( base_message == english->second.end() )
        {
          BOOST_ERROR( "'" << what << "' defines '" << message.first
                      << "', which no longer exists in '" << base
                      << ".xml' - retired key, or a misspelled id whose translation is dead" );
          continue;
        }

        // An *extra* placeholder is always broken: `WString::arg()` is called as many times as the
        //  English text has slots, so a slot beyond that is never filled and the user is shown a
        //  literal "{4}".
        vector<string> extra_args;
        set_difference( message.second.begin(), message.second.end(),
                        base_message->second.begin(), base_message->second.end(),
                        back_inserter( extra_args ) );

        BOOST_CHECK_MESSAGE( extra_args.empty(),
                            "'" << what << "' message '" << message.first << "' uses placeholder {"
                            << (extra_args.empty() ? string() : extra_args.front())
                            << "}, which the English text does not - it will render literally" );

        // Using *fewer* is only usually a defect.  Some slots carry a plural marker that a
        //  language does not need (French writes "gamma(s)" rather than "gamma{1}"), so this warns
        //  rather than fails - but it does catch a translation that quietly dropped a value.
        vector<string> unused_args;
        set_difference( base_message->second.begin(), base_message->second.end(),
                        message.second.begin(), message.second.end(),
                        back_inserter( unused_args ) );

        BOOST_WARN_MESSAGE( unused_args.empty(),
                           "'" << what << "' message '" << message.first << "' never uses {"
                           << (unused_args.empty() ? string() : unused_args.front())
                           << "}, which the English text does - check a value was not dropped" );
      }//for( loop over the locale's messages )
    }//for( loop over locales of this bundle )
  }//for( loop over bundle base names )

  BOOST_TEST_MESSAGE( "Compared " << num_locale_files << " locale bundles against "
                     << bundles.size() << " English bases" );
}//BOOST_AUTO_TEST_CASE( LocaleBundlesTrackTheirEnglishBase )


BOOST_AUTO_TEST_CASE( DetectionLimitBundlesAreFullyTranslated )
{
  set_app_text_dir();

  // The detection-limit tools are the reason this test exists, and unlike the rest of `app_text`
  //  they were swept to completeness deliberately.  Pinning that here means the next person to add
  //  a `dl-` key finds out immediately, rather than shipping it in English to every locale.
  const vector<string> bases = { "DetectionLimit", "DetectionLimitSimple", "DetectionLimitTool" };
  const vector<string> locales = { "ar", "bn", "es", "fr", "hi", "hu", "id", "ja", "pt", "ru",
                                   "zh" };

  for( const string &base : bases )
  {
    const string english_path = SpecUtils::append_path( g_app_text_dir, base + ".xml" );
    BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( english_path ),
                          "Missing '" << english_path << "'" );

    const BundleContents english = check_one_bundle( english_path );
    BOOST_REQUIRE( !english.empty() );

    for( const string &locale : locales )
    {
      const string path = SpecUtils::append_path( g_app_text_dir,
                                                  base + "_" + locale + ".xml" );
      if( !SpecUtils::is_file( path ) )
      {
        BOOST_ERROR( "Missing translation '" << base << "_" << locale << ".xml'" );
        continue;
      }

      const BundleContents translated = check_one_bundle( path );

      vector<string> missing;
      for( const auto &message : english )
      {
        if( !translated.count( message.first ) )
          missing.push_back( message.first );
      }

      BOOST_CHECK_MESSAGE( missing.empty(),
                          "'" << base << "_" << locale << ".xml' is missing "
                          << missing.size() << " of " << english.size() << " messages, starting "
                          << "with '" << (missing.empty() ? string() : missing.front()) << "'" );
    }//for( loop over locales )
  }//for( loop over bundle base names )
}//BOOST_AUTO_TEST_CASE( DetectionLimitBundlesAreFullyTranslated )
