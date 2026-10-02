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

#include <cctype>
#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>

#include <Wt/WLength.h>
#include <Wt/WString.h>
#include <Wt/WSuggestionPopup.h>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/MaterialDB.h"
#include "InterSpec/ShieldMaterialSuggestion.h"


namespace
{
  /** The popup matching options used historically by InterSpec for shielding
   material name suggestions; see the original `m_shieldingSuggestion` setup
   in `InterSpec::initMaterialDbAndSuggestions()`.
   */
  Wt::WSuggestionPopup::Options buildMaterialPopupOptions()
  {
    Wt::WSuggestionPopup::Options options;
    options.highlightBeginTag  = "<b>";
    options.highlightEndTag    = "</b>";
    options.listSeparator      = '\0';
    options.whitespace         = "";
    options.wordStartRegexp    = "\\s|^|\\(|\\<";
    options.appendReplacedText = "";
    return options;
  }


  /** How well `text` matches the typed `filter` (both lower-case), lower being better, or -1 if
   `text` does not contain it: 0 the whole text, 1 a whole word of it, 2 the start of it, 3 the start
   of a word, 4 anywhere.  Words are separated by whitespace, brackets, commas and the like, but not
   hyphens: "al" is a whole word of "al (aluminum)" - no material is named just "Al", yet typing it
   should rank aluminum first - while "u" is not a whole word of the alloy "u-al".
   */
  int match_rank( const std::string &text, const std::string &filter )
  {
    if( filter.empty() || (text == filter) )
      return 0;

    const auto is_separator = []( const char c ) -> bool {
      return std::isspace( static_cast<unsigned char>(c) )
             || ((c != '\0') && (std::strchr( "()[]{},;/&", c ) != nullptr));
    };

    int best = -1;
    for( size_t pos = text.find( filter ); pos != std::string::npos; pos = text.find( filter, pos + 1 ) )
    {
      const size_t end = pos + filter.size();
      const bool word_start = (pos == 0) || is_separator( text[pos - 1] );
      const bool word_end = (end == text.size()) || is_separator( text[end] );

      const int rank = (word_start && word_end) ? 1 : ((pos == 0) ? 2 : (word_start ? 3 : 4));
      if( (best < 0) || (rank < best) )
        best = rank;
    }//for( each occurrence of filter in text )

    return best;
  }//match_rank(...)
}//namespace


ShieldMaterialSuggestion::ShieldMaterialSuggestion()
  : Wt::WSuggestionPopup( buildMaterialPopupOptions() )
{
  addStyleClass( "suggestion" );

  // -1 means: re-emit `filterModel()` on every keystroke (so server-side
  //   results stay authoritative).  Same convention used by
  //   `IsotopeNameFilterModel` in `DecaySelectNuclideDiv.cpp`.
  setFilterLength( -1 );

  setMaximumSize( Wt::WLength::Auto,
                  Wt::WLength( 15, Wt::LengthUnit::FontEm ) );

  setDropDownIconUnfiltered( true );

  filterModel().connect( this, &ShieldMaterialSuggestion::handleFilter );

  // Don't seed eagerly: with `setFilterLength(-1)`, Wt's `modelRowsInserted`
  //   bails out unless an active filter is running, so calling `addSuggestion`
  //   here would set model data on rows that have no backing DOM widget and
  //   would crash inside `modelDataChanged`.  The first user-driven filter
  //   event (typing or dropdown-icon click) will populate the popup.
}//ShieldMaterialSuggestion::ShieldMaterialSuggestion()


void ShieldMaterialSuggestion::addFormulaMaterial( const std::string &name )
{
  if( name.empty() )
    return;

  for( const std::string &existing : m_formulaMaterials )
  {
    if( SpecUtils::iequals_ascii( existing, name ) )
      return;
  }

  if( MaterialDB::initialized() )
  {
    const std::shared_ptr<const MaterialDB> matDb = MaterialDB::instance();
    for( const std::string &existing : matDb->names() )
    {
      if( SpecUtils::iequals_ascii( existing, name ) )
        return;
    }
  }//if( MaterialDB::initialized() )

  m_formulaMaterials.push_back( name );
}//void addFormulaMaterial( const std::string &name )


void ShieldMaterialSuggestion::handleFilter( const Wt::WString &filter )
{
  clearSuggestions();

  const std::string filterStr = SpecUtils::to_lower_ascii_copy( filter.toUTF8() );

  // The popup lists rows in the order added and pre-selects the first one shown (what Enter/Tab
  //  takes), so the closest matches go first: typing "Al" should give "Al (aluminum)", not whichever
  //  "Al..." material the database happens to list first.
  struct Candidate
  {
    int rank;
    std::string text, lower;
  };
  std::vector<Candidate> candidates;

  const auto consider = [&candidates, &filterStr]( const std::string &text ){
    // Never add a blank row (the built-in "void" material has no description): WSuggestionPopup's
    //  matcher treats it as matching any input, shows it as "undefined", and picking it wipes the edit.
    if( SpecUtils::trim_copy( text ).empty() )
      return;

    std::string lower = SpecUtils::to_lower_ascii_copy( text );
    const int rank = match_rank( lower, filterStr );
    if( rank >= 0 )
      candidates.push_back( Candidate{ rank, text, std::move(lower) } );
  };

  if( MaterialDB::initialized() )
  {
    const std::shared_ptr<const MaterialDB> matDb = MaterialDB::instance();
    const std::vector<std::shared_ptr<const Material>> &mats = matDb->materials();
    for( const std::shared_ptr<const Material> &mat : mats )
    {
      if( !mat )
        continue;

      const std::string &name = mat->name;
      const std::string &desc = mat->description;

      // Filter out Pu-enrichment variants (e.g. "PuO2 - X.X% Pu240 Plutonium
      //  dioxide"); InterSpec does not use these and they only clutter the
      //  suggestion list.  Same filter the old code applied.
      if( name.find( "% Pu" ) != std::string::npos )
        continue;
      if( desc.find( "% Pu" ) != std::string::npos )
        continue;

      consider( name );

      if( SpecUtils::iequals_ascii( name, desc ) )
        continue;

      // Only add the description as a separate suggestion if it isn't
      //  already represented by the name string.
      if( SpecUtils::ifind_substr_ascii( name, desc.c_str() ) == std::string::npos )
        consider( desc );
    }//for( const std::shared_ptr<const Material> &mat : mats )
  }//if( MaterialDB::initialized() )

  for( const std::string &name : m_formulaMaterials )
    consider( name );

  // Best rank first; within a rank the shortest (fewest extra characters, so the closest
  //  completion), then alphabetical.  With nothing typed (the dropdown icon) keep the database order.
  if( !filterStr.empty() )
  {
    std::stable_sort( begin(candidates), end(candidates),
      []( const Candidate &lhs, const Candidate &rhs ) -> bool {
        if( lhs.rank != rhs.rank )
          return lhs.rank < rhs.rank;
        if( lhs.lower.size() != rhs.lower.size() )
          return lhs.lower.size() < rhs.lower.size();
        return lhs.lower < rhs.lower;
    } );
  }//if( !filterStr.empty() )

  for( const Candidate &candidate : candidates )
    addSuggestion( candidate.text, candidate.text );
}//void handleFilter( const Wt::WString &filter )
