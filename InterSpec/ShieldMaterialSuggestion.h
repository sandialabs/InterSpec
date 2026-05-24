#ifndef ShieldMaterialSuggestion_h
#define ShieldMaterialSuggestion_h
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

#include <string>
#include <vector>

#include <Wt/WSuggestionPopup.h>

namespace Wt
{
  class WString;
}


/** A self-contained `Wt::WSuggestionPopup` that dynamically populates its
 suggestions from the `MaterialDB` singleton in response to user typing.

 Each instance owns its own state (no shared global popup), so each widget
 with a material-name edit field can simply create one and call
 `forEdit( m_materialEdit, ... )`.  This avoids pre-loading the full material
 list into the DOM at page load and removes the need to plumb a shared
 `WSuggestionPopup *` through every widget constructor.

 Material names entered as chemical formulas (which are not part of the
 `MaterialDB` singleton) can be remembered for the lifetime of this popup
 via `addFormulaMaterial(...)`; that cache is per-instance only.
 */
class ShieldMaterialSuggestion : public Wt::WSuggestionPopup
{
public:
  ShieldMaterialSuggestion();

  ~ShieldMaterialSuggestion() override = default;

  /** Adds a user-entered material name (typically from
   `MaterialDB::materialFromChemicalFormula(...)`) to this popup's local
   cache, so it will appear as a suggestion on subsequent filter events.

   The cache is per-instance and is not shared across other
   `ShieldMaterialSuggestion` objects.  Duplicate entries (matched
   case-insensitively against the cache and against `MaterialDB` names)
   are silently ignored.
   */
  void addFormulaMaterial( const std::string &name );

private:
  /** Slot connected to `filterModel()`; rebuilds the suggestion list from
   `MaterialDB::instance()` plus the local formula cache, filtered against
   the user's current input.
   */
  void handleFilter( const Wt::WString &filter );

  std::vector<std::string> m_formulaMaterials;
};//class ShieldMaterialSuggestion

#endif //ShieldMaterialSuggestion_h
