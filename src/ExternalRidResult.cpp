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

#include "InterSpec/ExternalRidResult.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/IsotopeSearchByEnergy.h"

const SandiaDecay::Nuclide *ExternalRidIsotope::nuclide() const
{
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(source) )
    return std::get<const SandiaDecay::Nuclide *>(source);
  return nullptr;
}

const SandiaDecay::Element *ExternalRidIsotope::element() const
{
  if( std::holds_alternative<const SandiaDecay::Element *>(source) )
    return std::get<const SandiaDecay::Element *>(source);
  return nullptr;
}

const ReactionGamma::Reaction *ExternalRidIsotope::reaction() const
{
  if( std::holds_alternative<const ReactionGamma::Reaction *>(source) )
    return std::get<const ReactionGamma::Reaction *>(source);
  return nullptr;
}

bool ExternalRidIsotope::is_null() const
{
  return (!nuclide() && !element() && !reaction());
}

std::string ExternalRidIsotope::source_name() const
{
  if( const SandiaDecay::Nuclide * const nuc = nuclide() )
    return nuc->symbol;
  else if( const SandiaDecay::Element * const el = element() )
    return el->symbol;
  else if( const ReactionGamma::Reaction * const rctn = reaction() )
    return rctn->name();
  return name;
}


void ExternalRidIsotope::init()
{
  // If source is already set, we might still need to fill in type
  if( std::holds_alternative<std::monostate>(source) && !name.empty() )
  {
    // Try to resolve name to a Nuclide first
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    if( !db )
      throw std::runtime_error( "Failed to get SandiaDecayDataBase!?!" );
    
    const SandiaDecay::Nuclide * const nuc = db->nuclide(name);
    if( nuc )
      source = nuc;
    
    if( !nuc )
    {
      if( const SandiaDecay::Element * const el = db->element(name) )
      {
        source = el;
      }else if( const ReactionGamma * const rctn_db = ReactionGammaServer::database() )
      {
        std::vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
        rctn_db->gammas( name, possible_rctns );
        
        if( !possible_rctns.empty() && possible_rctns[0].reaction ) // Take the first reaction found
          source = possible_rctns[0].reaction;
      }
    }//if( !nuc )
  }//if( std::holds_alternative<std::monostate>(source) && !name.empty() )
  
  
  // If type is empty and we have a source, try to categorize it
  if( type.empty() && !is_null() )
  {
    InterSpec * const interspec = InterSpec::instance();
    IsotopeSearchByEnergy * const search = interspec ? interspec->nuclideSearch() : nullptr;
    if( search )
    {
      const std::vector<IsotopeSearchByEnergy::NucSearchCategory> &categories = search->search_categories();
      
      auto get_type = [&]( auto src ) -> std::string {
        if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_medical_category_key, categories) )
          return "Medical";
        else if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_industrial_category_key, categories) )
          return "Industrial";
        else if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_norm_category_key, categories) )
          return "NORM";
        else if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_snm_category_key, categories) )
          return "SNM";
        return "";
      };//get_type
      
      if( const SandiaDecay::Nuclide * const nuc = nuclide() )
        type = get_type( nuc );
      else if( const SandiaDecay::Element * const el = element() )
        type = get_type( el );
      else if( const ReactionGamma::Reaction * const rctn = reaction() )
        type = get_type( rctn );
    }//if( search )
  }//if( type.empty() && !is_null() )
}//void init()


bool ExternalRidIsotope::operator==(const ExternalRidIsotope& other) const
{
  return (name == other.name && 
          type == other.type &&
          confidenceStr == other.confidenceStr &&
          countRate == other.countRate &&
          confidence == other.confidence &&
          source == other.source);
}
