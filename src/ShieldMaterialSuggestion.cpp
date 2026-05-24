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
#include <string>
#include <vector>

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


  /** True if `filter` is empty, OR `name` contains `filter` (case-insensitive). */
  bool name_matches( const std::string &name, const std::string &filter )
  {
    if( filter.empty() )
      return true;
    return SpecUtils::ifind_substr_ascii( name, filter.c_str() ) != std::string::npos;
  }
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

  const std::string filterStr = filter.toUTF8();

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

      const bool nameHit = name_matches( name, filterStr );
      const bool descHit = name_matches( desc, filterStr );
      if( !nameHit && !descHit )
        continue;

      addSuggestion( name, name );

      if( SpecUtils::iequals_ascii( name, desc ) )
        continue;

      // Only add the description as a separate suggestion if it isn't
      //  already represented by the name string.
      const size_t sub_pos = SpecUtils::ifind_substr_ascii( name, desc.c_str() );
      if( sub_pos == std::string::npos )
        addSuggestion( desc, desc );
    }//for( const std::shared_ptr<const Material> &mat : mats )
  }//if( MaterialDB::initialized() )

  for( const std::string &name : m_formulaMaterials )
  {
    if( name_matches( name, filterStr ) )
      addSuggestion( name, name );
  }
}//void handleFilter( const Wt::WString &filter )
