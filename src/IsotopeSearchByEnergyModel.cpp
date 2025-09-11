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

#include <set>
#include <map>
#include <deque>
#include <vector>
#include <sstream>

#include <Wt/WServer>
#include <Wt/WApplication>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"


using namespace Wt;
using namespace std;

namespace
{
  struct Sorter
  {
    const vector<double> &m_energies;
    IsotopeSearchByEnergyModel::Column m_column;
    Wt::SortOrder m_order;
    
    Sorter( const vector<double> &energies, int column, Wt::SortOrder order )
    : m_energies( energies ),
    m_column( IsotopeSearchByEnergyModel::Column(column) ),
    m_order( order )
    {}
    
    bool operator()( const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &lhsin,
                    const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &rhsin )
    {
      const bool asscending = (m_order==AscendingOrder);
      
      assert( lhsin.size()==rhsin.size() && lhsin.size()==m_energies.size() );
      
      if( m_energies.empty() )
        return asscending;
      
      const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &lhs = (asscending ? lhsin : rhsin);
      const vector<IsotopeSearchByEnergyModel::IsotopeMatch> &rhs = (asscending ? rhsin : lhsin);
      
      switch( m_column )
      {
        case IsotopeSearchByEnergyModel::ParentIsotope:
        case IsotopeSearchByEnergyModel::SpecificIsotope:
          return (lhs[0].m_displayData[m_column] < rhs[0].m_displayData[m_column] );
          
        case IsotopeSearchByEnergyModel::Distance:
          return (lhs[0].m_distance < rhs[0].m_distance);
          
        case IsotopeSearchByEnergyModel::Energy:
        {
          return (std::stod(lhs[0].m_displayData[m_column].narrow())
                  < std::stod(rhs[0].m_displayData[m_column].narrow()));
        }
          
        case IsotopeSearchByEnergyModel::ProfileDistance:
          return (lhs[0].m_profileDistance < rhs[0].m_profileDistance);
          
        case IsotopeSearchByEnergyModel::BranchRatio:
        {
          double lhsval = 0.0, rhsval = 0.0;
          for( size_t i = 0; i < lhs.size(); ++i )
            lhsval += std::stod(lhs[i].m_displayData[m_column].narrow());
          for( size_t i = 0; i < rhs.size(); ++i )
            rhsval += std::stod(rhs[i].m_displayData[m_column].narrow());
          return lhsval < rhsval;
        }
          
        case IsotopeSearchByEnergyModel::ParentHalfLife:
          if( lhs[0].m_nuclide && rhs[0].m_nuclide )
            return (lhs[0].m_nuclide->halfLife < rhs[0].m_nuclide->halfLife);
          return (lhs[0].m_nuclide != NULL);
          
        case IsotopeSearchByEnergyModel::AssumedAge:
          return (lhs[0].m_age < rhs[0].m_age);
          
        case IsotopeSearchByEnergyModel::NumColumns:
          return false;
      }//switch( col )
      
      return false;
    }//bool operator()
  };//struct Sorter
  
  
  /** Conviencince typedef from (pointer to nuclide) to the set of energies of
     that nuclide
   */
  typedef map<const SandiaDecay::Nuclide *, set<double> > NuclideMatches;
  
  /** Returns nuclides filtered such that they have gammas/x-rays in the
      specified energy ranges, with at least the minimum branching ratio and
      half-lives specified.
   */
  NuclideMatches filter_nuclides( const double min_relative_br,
                                  const double minHalfLife,
                                  const bool no_progeny,
                                  const vector<double> &energies,
                                  const vector<double> &windows )
  {
    using namespace SandiaDecay;
    
    assert( energies.size() == windows.size() );
    
    NuclideMatches filteredNuclides;
    
    //Get isotopes with gammas in all ranges...
    EnergyToNuclideServer::setLowerLimits( minHalfLife, min_relative_br );

    auto nucnuc = EnergyToNuclideServer::energyToNuclide();
    if( !nucnuc )
      throw runtime_error( "Couldnt get EnergyToNuclideServer" );

    for( size_t i = 0; i < energies.size(); ++i )
    {
      const float minenergy = static_cast<float>(energies[i] - windows[i]);
      const float maxenergy = static_cast<float>(energies[i] + windows[i]);
      
      const bool canBeAnnih = ((510.99891f >= minenergy) && (510.99891f <= maxenergy));
      

      //nucnuc should be sorted by energy (EnergyNuclidePair::operator<)
      const auto begin = lower_bound( nucnuc->begin(), nucnuc->end(),
                                   EnergyToNuclideServer::EnergyNuclidePair{minenergy, nullptr} );
      const auto end = upper_bound( begin, nucnuc->end(),
                                   EnergyToNuclideServer::EnergyNuclidePair{maxenergy, nullptr} );
      
      for( auto pos = begin; pos != end; ++pos )
      {
        //nucnuc actually only contains the nuclides that they themselves give
        //  off the requested energies, meaning we have to go through and inspect
        //  all the nuclides that could decay to the nuclide in nucnuc to see
        //  if they are compatible with the search criteria, if we are considering progeny.
        int stop = 0;
        double transbr = 1.0;
        for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
        {
          for( const SandiaDecay::RadParticle &r : t->products )
          {
            if( ((r.type==SandiaDecay::GammaParticle) || (r.type==SandiaDecay::XrayParticle))
               && fabs(r.energy - pos->energy) < 0.001 )
            {
              transbr = t->branchRatio * r.intensity;
              stop = 1;
            }else if( canBeAnnih && r.type==SandiaDecay::PositronParticle )
            {
              //XXX - there is a slight issue with nuclides that have positrons and gammas near 511 keV, but whatever
              if( !stop )
                transbr = 0.0;
              transbr += 2.0 * t->branchRatio * r.intensity;
              stop = 2;
            }
          }//for( const SandiaDecay::RadParticle &r : t->products )
          
          if( stop == 1 )
            break;
        }//for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
        
        if( no_progeny )
        {
          filteredNuclides[pos->nuclide].insert( energies[i] );
        }else
        {
          const vector<const Nuclide *> forebearers = pos->nuclide->forebearers();
          for( const Nuclide *nuc : forebearers )
          {
            if( nuc->halfLife < minHalfLife )
              continue;
            
            filteredNuclides[nuc].insert( energies[i] );
          }
        }//if( !no_progeny )
      }//for( pos = begin; pos != end; ++pos )
    }//for( size_t i = 0; i < energies.size(); ++i )

    if( min_relative_br > 0.0 )
    {
      //Age each nuclide to default nuclide age; find its max rel intensity line, then divide each gamma by that
      //  and re-due the above.  This pathway is a lot slower since we are doing the full decay calculation, so we'll
      //  try to use multiple threads.
      SpecUtilsAsync::ThreadPool pool;
      for( NuclideMatches::value_type &nuc_matches : filteredNuclides )
      {
        const SandiaDecay::Nuclide * const nuc = nuc_matches.first;
        set<double> * const matched_energies = &(nuc_matches.second);

        pool.post( [min_relative_br, nuc, no_progeny, matched_energies, &energies, &windows](){
          matched_energies->clear();
          const double age = no_progeny ? 0.0 : PeakDef::defaultDecayTime( nuc );
          SandiaDecay::NuclideMixture mixture;
          mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuc,1.0E4) );
          const vector<SandiaDecay::EnergyRatePair> photons
                                        = mixture.photons(age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy);
          double max_intensity = 0.0;
          for( const SandiaDecay::EnergyRatePair &energy_rate : photons )
            max_intensity = std::max( max_intensity, energy_rate.numPerSecond );

          if( (max_intensity <= 0.0) || IsNan(max_intensity) || IsInf(max_intensity) )
            max_intensity = 1.0; //JIC

          for( size_t window_index = 0; window_index < energies.size(); ++window_index )
          {
            const double min_energy = energies[window_index] - fabs(windows[window_index]);
            const double max_energy = energies[window_index] + fabs(windows[window_index]);

            const auto photon_begin = lower_bound( begin(photons), end(photons), min_energy,
              []( const SandiaDecay::EnergyRatePair &el, const double value ){
                return el.energy < value;
            } );

            const auto photon_end = upper_bound( begin(photons), end(photons), max_energy,
              []( const double value, const SandiaDecay::EnergyRatePair &el ){
                return value < el.energy;
            } );

            auto max_in_range_iter = photon_end;
            for( auto iter = photon_begin; iter != photon_end; ++iter )
            {
              if( (max_in_range_iter == photon_end) || (iter->numPerSecond > max_in_range_iter->numPerSecond) )
                max_in_range_iter = iter;
            }
            if( (max_in_range_iter != photon_end) && ( (max_in_range_iter->numPerSecond/max_intensity) >= min_relative_br) )
              matched_energies->insert( max_in_range_iter->energy );
          }//for( size_t window_index = 0; window_index < energies.size(); ++window_index )
        } );
      }//for( const auto &nuc_matches : filteredNuclides )
      pool.join();
    }//if( we are filtering on BRs )

    return filteredNuclides;
  }//filter_nuclides( )


  void alphaOrBetasWithAllEnergies( const bool isAlpha,
                                 const vector<double> &energies,
                                 const vector<double> &windows,
                                 const shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                 const vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                 const vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                 const bool no_progeny,
                                 const vector<const SandiaDecay::Nuclide *> &nuclides,
                                   IsotopeSearchByEnergyModel::SearchResults &answer )
  {
    if( energies.empty() )
      return;
    
    assert( energies.size() == windows.size() );
    
    char buffer[32] = { '\0' };
    
    const SandiaDecay::ProductType wantedType = isAlpha ? SandiaDecay::ProductType::AlphaParticle
    : SandiaDecay::ProductType::BetaParticle;
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    
    const vector<const SandiaDecay::Nuclide *> &all_nuclides = db->nuclides();
    for( const SandiaDecay::Nuclide * nuc : all_nuclides )
    {
      if( !nuclides.empty() )
      {
        const auto pos = std::find( begin(nuclides), end(nuclides), nuc );
        if( pos == end(nuclides) )
          continue;
      }
      
      // We'll require the parent nuclide to match at least one energy Window, so we dont have
      //  to decay everything - I'm not really sure if this is valid, or what, so it might get
      //  deleted later.
      bool parent_matches_one = false;
      for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &part : trans->products )
        {
          if( part.type != wantedType )
            continue;
          
          for( size_t windex = 0; !parent_matches_one && (windex < energies.size()); windex += 1 )
          {
            const double min_energy = energies[windex] - windows[windex];
            const double max_energy = energies[windex] + windows[windex];
            parent_matches_one = ((part.energy >= min_energy) && (part.energy <= max_energy));
          }
          
          if( parent_matches_one )
            break;
        }//for( const RadParticle &part : trans->products )
        
        if( parent_matches_one )
          break;
      }//for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
      
      if( !parent_matches_one )
        continue;
      
      bool matches_all_windows = true;
      
      //We will track the Transition, and particle index that best each of the windows;
      //  upon sucess there will be the same number of elements in this vector as search
      //  windows, with them ordered the same.
      vector<tuple<const SandiaDecay::Transition *,size_t,double>> matching_trans_part;
      
      // If we dont want progeny, we will use an age of 0.0, which should make activity of all
      //  progeny 0, which we will then skip over.
      const double age = no_progeny ? 0.0 : PeakDef::defaultDecayTime( nuc );
      
      SandiaDecay::NuclideMixture mixture;
      mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuc,1.0) );
      const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity(age);
      
      for( size_t windex = 0; matches_all_windows && (windex < energies.size()); windex += 1 )
      {
        const double min_energy = energies[windex] - windows[windex];
        const double max_energy = energies[windex] + windows[windex];
        
        bool match_this_window = false;
        for( const SandiaDecay::NuclideActivityPair &nap : activities )
        {
          if( !nap.nuclide || (nap.activity <= 0.0) )
            continue;
          
          for( const SandiaDecay::Transition *trans : nap.nuclide->decaysToChildren )
          {
            for( size_t trans_index = 0; trans_index < trans->products.size(); ++trans_index )
            {
              const SandiaDecay::RadParticle &part = trans->products[trans_index];
              
              const bool matches = ((part.type == wantedType)
                                    && (part.energy >= min_energy)
                                    && (part.energy <= max_energy)
                                    && (part.intensity > 0.0) );
              if( matches )
              {
                match_this_window = true;
                
                if( matching_trans_part.size() <= windex )
                {
                  matching_trans_part.emplace_back(trans, trans_index, nap.activity);
                }else
                {
                  assert( (windex + 1) == matching_trans_part.size() );
                  
                  // We already found another match - let see if this is a better one
                  // TODO: we are keeping the one closest - but we should consider amplitude
                  const SandiaDecay::Transition *prev_trans = std::get<0>(matching_trans_part[windex]);
                  const size_t prev_part_index = std::get<1>(matching_trans_part[windex]);
                  assert( prev_trans );
                  assert( prev_part_index < prev_trans->products.size() );
                  const double prev_de = fabs( energies[windex] - prev_trans->products[prev_part_index].energy );
                  const double this_de = fabs( energies[windex] - part.energy );
                  if( this_de < prev_de )
                  {
                    std::get<0>(matching_trans_part[windex]) = trans;
                    std::get<1>(matching_trans_part[windex]) = trans_index;
                    std::get<2>(matching_trans_part[windex]) = nap.activity;
                  }
                }//if( this is first particle matching this search energy ) / else (see if this is better match)
              }//if( matches )
            }//for( const RadParticle &part : trans->products )
          }//for( const SandiaDecay::Transition *trans : nap.nuclide->decaysToChildren )
        }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
        
        if( !match_this_window )
          matches_all_windows = false;
      }//for( size_t windex = 0; !matches_all_windows && (windex < energies.size()); windex += 1 )
      
      if( !matches_all_windows )
        continue;
      
      const vector<SandiaDecay::EnergyRatePair> particle_rates
      = isAlpha ? mixture.alphas(age) : mixture.betas(age);
      
      double max_rate = 0;
      for( const SandiaDecay::EnergyRatePair &erp : particle_rates )
        max_rate = std::max( max_rate, erp.numPerSecond );
      
      assert( matching_trans_part.size() == energies.size() );
      
      double total_dist = 0.0;
      for( size_t i = 0; i < matching_trans_part.size(); ++i )
      {
        const double energy = energies[i];
        
        const SandiaDecay::Transition *const trans = std::get<0>(matching_trans_part[i]);
        const size_t trans_index = std::get<1>(matching_trans_part[i]);
        assert( trans );
        assert( trans_index < trans->products.size() );
        const SandiaDecay::RadParticle &part = trans->products[trans_index];
        total_dist += fabs( energy - part.energy );
      }
      
      // Use `matching_trans_part` to generate results
      vector<IsotopeSearchByEnergyModel::IsotopeMatch> nucmatches;
      for( size_t i = 0; i < matching_trans_part.size(); ++i )
      {
        const double energy = energies[i];
        
        const SandiaDecay::Transition *const trans = std::get<0>(matching_trans_part[i]);
        const size_t trans_index = std::get<1>(matching_trans_part[i]);
        const double activity = std::get<2>(matching_trans_part[i]);
        assert( trans );
        assert( trans_index < trans->products.size() );
        
        const SandiaDecay::RadParticle &part = trans->products[trans_index];
        
        
        IsotopeSearchByEnergyModel::IsotopeMatch match;
        match.m_distance = total_dist;
        match.m_age = age;
        match.m_branchRatio = activity * trans->branchRatio * part.intensity / max_rate;
        match.m_profileDistance = 0.0; //We'll calc this a little later, if i==0
        match.m_nuclide = nuc;
        match.m_transition = trans;
        match.m_particle = &part;
        match.m_element = nullptr;
        match.m_xray = nullptr;
        match.m_reaction = nullptr;
        
        match.m_displayData[IsotopeSearchByEnergyModel::ParentIsotope] = nuc->symbol;
        snprintf( buffer, sizeof(buffer), "%.2f", part.energy );
        match.m_displayData[IsotopeSearchByEnergyModel::Energy] = buffer;
        
        snprintf( buffer, sizeof(buffer), "%.2f", total_dist );
        match.m_displayData[IsotopeSearchByEnergyModel::Distance] = buffer;
        
        snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
        match.m_displayData[IsotopeSearchByEnergyModel::BranchRatio] = buffer;
        
        if( trans->parent && trans->child )
        {
          match.m_displayData[IsotopeSearchByEnergyModel::SpecificIsotope]
                    = trans->parent->symbol + "&rarr;" + trans->child->symbol;
        }else if( trans->parent )
        {
          match.m_displayData[IsotopeSearchByEnergyModel::SpecificIsotope]
                    = SandiaDecay::to_str(trans->mode) + string(" of ") + trans->parent->symbol;
        }
        
        if( !i )
        {
          match.m_displayData[IsotopeSearchByEnergyModel::ParentHalfLife]
                        = PhysicalUnitsLocalized::printToBestTimeUnits(match.m_nuclide->halfLife);
          match.m_displayData[IsotopeSearchByEnergyModel::AssumedAge]
                                      = PhysicalUnitsLocalized::printToBestTimeUnits(match.m_age);
          
          if( isAlpha )
          {
            const double weight = IsotopeId::profile_weight( nullptr, displayed_measurement,
                                                 user_peaks, automated_search_peaks,
                                                 particle_rates, energies, windows,
                                                 26, 0.0,
                                                 match.m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope] );
            
            match.m_profileDistance = weight;
            snprintf( buffer, sizeof(buffer), "%.2f", weight );
            match.m_displayData[IsotopeSearchByEnergyModel::ProfileDistance] = buffer;
          }else
          {
            // Since betas are endpoint energies, it makes no sense to use profile_weight
            match.m_profileDistance = -total_dist;
            snprintf( buffer, sizeof(buffer), "%.2f", -total_dist );
            match.m_displayData[IsotopeSearchByEnergyModel::ProfileDistance] = buffer;
          }//if( isAlpha ) / else
        }//if( !i )
        
        
        nucmatches.push_back( match );
      }//for( size_t i = 0; i < matching_trans_part.size(); ++i )
      
      assert( nucmatches.size() == energies.size() );
      answer.push_back( nucmatches );
    }//for( const SandiaDecay::Nuclide * nuc : all_nuclides )
  }//void alphaOrBetasWithAllEnergies(...)

}//namespace


IsotopeSearchByEnergyModel::IsotopeMatch::IsotopeMatch()
: m_distance(0.0), m_age(0.0), m_branchRatio(0.0), m_profileDistance(-1.0),
m_nuclide(0), m_transition(0), m_particle(0),
m_sourceGammaType(PeakDef::NormalGamma),
m_element(0), m_xray(0), m_reaction(0)
{
}


IsotopeSearchByEnergyModel::IsotopeSearchByEnergyModel( Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_sortColumn( IsotopeSearchByEnergyModel::Column::ProfileDistance ),
  m_sortOrder( Wt::AscendingOrder )
{
  
}//IsotopeSearchByEnergyModel( constuctor )


IsotopeSearchByEnergyModel::~IsotopeSearchByEnergyModel()
{
  
}//IsotopeSearchByEnergyModel destructor


void IsotopeSearchByEnergyModel::nuclidesWithAllEnergies(
            const IsotopeSearchByEnergyModel::NucToEnergiesMap &filteredNuclides,
            const vector<double> &energies,
            const vector<double> &windows,
            const double min_rel_br,
            const bool no_progeny,
            const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
            const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
            const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
            const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
            std::vector< vector<IsotopeSearchByEnergyModel::IsotopeMatch> > &answer )
{
  // Note: most parts of this function are wrapped in try/catch blocks, but exceptions are really
  //       not expected.
  
  assert( energies.size() == windows.size() );
  
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  // TODO: currently only taking 'most likely' combination of matching of
  //       energies - should implement getting all permutations!
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  for( const NucToEnergiesMap::value_type &nm : filteredNuclides )
  {
    //check to see if this nuclide has gammas for each energy, not strictly
    //  necessary, but probably computationally faster (do we care about this
    //  here though)
    bool hasAll = true;
    vector<size_t> energy_xray_is_for;
    vector<const SandiaDecay::EnergyIntensityPair *> matching_xrays;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      try
      {
        const double energy = energies[i];
        const double de = fabs( windows[i] );
        set<double>::const_iterator lb, up, iter;
        lb = nm.second.lower_bound( energy - de );
        up = nm.second.upper_bound( energy + de );
        bool found = false;
        for( iter = lb; iter != up; ++iter )
          found |= (fabs((*iter)-energy) <= de);

        if( !found && (energy < 125.0*PhysicalUnits::keV) )
        {//lets look through the fluorescent x-rays for this element
          double minxraydist = 999.9;
          SandiaDecay::EnergyIntensityPair closest_xray(0.0,0.0);
          const SandiaDecay::Element *el = db->element( nm.first->atomicNumber );
          const vector<SandiaDecay::EnergyIntensityPair> &xrays = el->xrays;
          for( size_t j = 0; j < xrays.size(); ++j )
          {
            const double dist = fabs(xrays[j].energy-energy);
            if( dist <= minxraydist )
            {
              minxraydist = dist;
              closest_xray = xrays[j];
            }//if( dist <= minxraydist )
            
            if( minxraydist <= de )
            {
              found = true;
              energy_xray_is_for.push_back( i );
              matching_xrays.push_back( &(xrays[j]) );
            }//if( minxraydist <= de )
          }//for( loop iver x-rays )
        }//if( !found && (energy < 115.0*PhysicalUnits::keV) )
        
        if( !found )
        {
          hasAll = false;
          break;
        }//if( !found )
      }catch( std::exception &e )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        stringstream msg;
        msg << "Unexpected exception (1): " << e.what();
        log_developer_error( __func__, msg.str().c_str() );
#endif
      }//try / catch
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( !hasAll )
      continue;

    double dist = 0.0;
    vector<IsotopeMatch> nucmatches;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      try
      {
        vector<size_t>::const_iterator xraypos = find( energy_xray_is_for.begin(),
                                                      energy_xray_is_for.end(), i );
        
        if( xraypos != energy_xray_is_for.end() )
        {
          const size_t pos = xraypos - energy_xray_is_for.begin();
          const double energy = energies[i];
          const SandiaDecay::EnergyIntensityPair *xray = matching_xrays[pos];
          
          IsotopeMatch match;
          match.m_distance = fabs(energy - xray->energy);    //sum of distance over all energies searched
          dist += match.m_distance;
          
          match.m_age = 0.0;         //age assumed for listing things
          match.m_branchRatio = xray->intensity;
          
          //Only one of the following will be valid: m_nuclide, m_element, m_reaction
          match.m_nuclide = NULL;
          match.m_transition = NULL;
          match.m_particle = NULL;
          match.m_element = db->element( nm.first->atomicNumber );
          match.m_xray = xray;
          match.m_reaction = NULL;
          //        match.m_reactionEnergy;
          
          //        match.m_displayData[ParentHalfLife] = "";
          match.m_displayData[AssumedAge] = PhysicalUnitsLocalized::printToBestTimeUnits( 0.0 );
          //        match.m_displayData[SpecificIsotope] = "";
          if( match.m_element )
            match.m_displayData[ParentIsotope] = match.m_element->symbol;
          
          if( xray )
          {
            snprintf( buffer, sizeof(buffer), "%.2f", xray->energy );
            match.m_displayData[Energy] = buffer;
          }//if( xray )
          
          snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
          match.m_displayData[Distance] = buffer;
          
          snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
          match.m_displayData[BranchRatio] = buffer;
          
          nucmatches.push_back( match );
        }else//if( xraypos != energy_xray_is_for.end() )
        {
          size_t trans_index = 0;
          const SandiaDecay::Transition *transition = NULL;
          IsotopeMatch match;
          
          PeakDef::SourceGammaType sourceGammaType = PeakDef::NormalGamma;
          
          const double nucAge = no_progeny ? 0.0 : -1.0;
          PeakDef::findNearestPhotopeak( nm.first, energies[i], windows[i],
                                          false, nucAge, transition, trans_index, sourceGammaType );
          
          if( !transition && (sourceGammaType != PeakDef::AnnihilationGamma) )
            continue;
          
          match.m_nuclide = nm.first;
          match.m_transition = transition;
          if( transition )
            match.m_particle = &(transition->products[trans_index]);
          match.m_sourceGammaType = sourceGammaType;
          
          match.m_distance = 99999999.9;
          switch( sourceGammaType )
          {
            case PeakDef::NormalGamma:
            case PeakDef::XrayGamma:
              if( match.m_particle )
                match.m_distance = fabs(energies[i] - match.m_particle->energy);
              break;
              
            case PeakDef::AnnihilationGamma:
              match.m_distance = fabs(energies[i] - 510.99891*SandiaDecay::keV );
              break;
              
            case PeakDef::SingleEscapeGamma:
              if( match.m_particle )
                match.m_distance = fabs(energies[i] - (match.m_particle->energy - 510.99891));
              break;
              
            case PeakDef::DoubleEscapeGamma:
              if( match.m_particle )
                match.m_distance = fabs(energies[i] - (match.m_particle->energy - 2.0*510.99891) );
              break;
          }//switch( sourceGammaType )
          
          match.m_age = no_progeny ? 0.0 : PeakDef::defaultDecayTime( nm.first );
          SandiaDecay::NuclideMixture mixture;
          mixture.addNuclide( SandiaDecay::NuclideActivityPair(nm.first,1.0) );
          const vector<SandiaDecay::EnergyRatePair> photons
                  = mixture.photons( match.m_age, SandiaDecay::NuclideMixture::OrderByAbundance );
          
          double nearestEnergy = 999999.9, nearestAbun = 0.0, maxAbund = -999.9;
          for( const SandiaDecay::EnergyRatePair &aep : photons )
          {
            double d = 999999.9;
            
            if( match.m_sourceGammaType == PeakDef::AnnihilationGamma )
              d = fabs( aep.energy - 510.99891*SandiaDecay::keV );
            else if( match.m_particle )
              d = fabs( aep.energy - match.m_particle->energy );
            
            if( d < nearestEnergy )
            {
              nearestEnergy = d;
              nearestAbun = aep.numPerSecond;
            }//if( d < nearestEnergy )
            
            maxAbund = std::max( maxAbund, aep.numPerSecond );
          }//for( const SandiaDecay::AbundanceEnergyPair &aep : photons )
          
          match.m_branchRatio = nearestAbun / maxAbund;
          
          if( (match.m_branchRatio < min_rel_br) || (match.m_branchRatio <= 0.0) )
            continue;

          match.m_displayData[ParentIsotope] = match.m_nuclide->symbol;
          
          if( match.m_sourceGammaType == PeakDef::AnnihilationGamma )
            snprintf( buffer, sizeof(buffer), "510.99" );
          else if( match.m_particle )
            snprintf( buffer, sizeof(buffer), "%.2f", match.m_particle->energy );
          
          match.m_displayData[Energy] = buffer;
          
          snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
          match.m_displayData[Distance] = buffer;
          
          snprintf( buffer, sizeof(buffer), "%.2f", match.m_branchRatio );
          match.m_displayData[BranchRatio] = buffer;
          
          string xstn_str;
          if( match.m_transition && match.m_transition->parent && match.m_transition->child )
          {
            xstn_str = match.m_transition->parent->symbol + "&rarr;"
                        + match.m_transition->child->symbol;
          }else if( match.m_transition && match.m_transition->parent )
          {
            xstn_str = string(SandiaDecay::to_str(match.m_transition->mode)) + " of "
                        + match.m_transition->parent->symbol;
          }else if( !match.m_transition )
          {
            xstn_str = "Annih. Gamma";
          }//if( match.m_transition->parent... ) / else
          
          if( match.m_transition && (match.m_particle->type == SandiaDecay::XrayParticle) )
            xstn_str = " xray";
          
          match.m_displayData[SpecificIsotope] = std::move(xstn_str);
          
          if( !i )
          {
            match.m_displayData[ParentHalfLife]
            = PhysicalUnitsLocalized::printToBestTimeUnits(match.m_nuclide->halfLife);
            match.m_displayData[AssumedAge]
            = PhysicalUnitsLocalized::printToBestTimeUnits(match.m_age);
          }//if( !i )
          
          dist += match.m_distance;
          
          nucmatches.push_back( match );
        }//if( xraypos != energy_xray_is_for.end() ) / else
      }catch( std::exception &e )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        stringstream msg;
        msg << "Unexpected exception (2): " << e.what();
        log_developer_error( __func__, msg.str().c_str() );
#endif
      }//try/catch
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( nucmatches.size() != energies.size() )
      continue;
    
    try
    {
      if( !nucmatches[0].m_nuclide )
      {
        for( IsotopeMatch &match : nucmatches )
        {
          if( match.m_nuclide )
          {
            nucmatches[0].m_nuclide = match.m_nuclide;
            nucmatches[0].m_displayData[ParentIsotope] = match.m_nuclide->symbol;
            nucmatches[0].m_displayData[Energy]
                            = nucmatches[0].m_displayData[Energy].toUTF8() + " (xray)";
            break;
          }//if( match.m_nuclide )
        }//for( IsotopeMatch &match : nmagicucmatches )
      }//if( !nucmatches[0].m_nuclide )
      
      nucmatches[0].m_distance = dist;
      
      const double gcm2 = PhysicalUnits::g / PhysicalUnits::cm2;
      const double atomic_nums[]   = { 1.0, 26.0, 74.0 };
      const double areal_density[] = { 0.0*gcm2, 10.0*gcm2, 25.0*gcm2 };
      static_assert( sizeof(atomic_nums) == 3*sizeof(atomic_nums[0]), "" );
      static_assert( sizeof(areal_density) == 3*sizeof(areal_density[0]), "" );
      
      double mw = -999.9;
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity( nucmatches[0].m_nuclide, 0.001*SandiaDecay::curie );
      vector<SandiaDecay::EnergyRatePair> srcgammas = mix.photons( nucmatches[0].m_age );
      
      for( size_t i = 0; i < 3; ++i )
      {
        const double weight = IsotopeId::profile_weight( detector_response_function,
                                        displayed_measurement,
                                        user_peaks,
                                        automated_search_peaks, srcgammas,
                                        energies, windows,
                                        atomic_nums[i], areal_density[i],
                                        nucmatches[0].m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope] );
        mw = std::max(mw,weight);
      }
      nucmatches[0].m_profileDistance = mw;
      snprintf( buffer, sizeof(buffer), "%.2f", mw );
      nucmatches[0].m_displayData[ProfileDistance] = buffer;
      
      
      snprintf( buffer, sizeof(buffer), "%.2f", dist );
      nucmatches[0].m_displayData[Distance] = buffer;
      answer.push_back( nucmatches );
    }catch( std::exception &e )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      stringstream msg;
      msg << "Unexpected exception (3): " << e.what();
      log_developer_error( __func__, msg.str().c_str() );
#endif
    }//try / catch
  }//for( const NuclideMatches::value_type &nm : filteredNuclides )
}//void nuclidesWithAllEnergies


void IsotopeSearchByEnergyModel::xraysWithAllEnergies(
                                                      const std::vector<double> &energies,
                                                      const std::vector<double> &windows,
                                                      const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                                      const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                                      const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                                      const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                                      const std::vector<const SandiaDecay::Element *> &limited_elements,
                                                      SearchResults &answer )
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  for( size_t i = 0; i < energies.size(); ++i )
    if( (energies[i]-windows[i]) > 120*PhysicalUnits::keV )
      return;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const vector<const SandiaDecay::Element *> &elements = db->elements();
  
  for( const SandiaDecay::Element *el : elements )
  {
    if( !limited_elements.empty() )
    {
      auto pos = std::find( begin(limited_elements), end(limited_elements), el );
      if( pos == end(limited_elements) )
        continue;
    }//if( !limited_elements.empty() )
    
    vector<IsotopeMatch> nucmatches;
    vector<const SandiaDecay::EnergyIntensityPair *> xray_matches;
    const vector<SandiaDecay::EnergyIntensityPair> &xrays = el->xrays;
    if( xrays.size() < energies.size() )
      continue;
    
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      const double de = fabs( windows[i] );
      double minEnergy = 1000.0 * de;
      const SandiaDecay::EnergyIntensityPair *xray = NULL;
      
      for( const SandiaDecay::EnergyIntensityPair &e : xrays )
      {
        const double delta = fabs((e.energy)-energy);
        if( delta < minEnergy )
        {
          minEnergy = delta;
          xray = &e;
        }//if( delta < minEnergy )
      }//for( const SandiaDecay::EnergyIntensityPair &e : xrays )
      
      if( minEnergy < de )
        xray_matches.push_back( xray );
      else
        break;
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( xray_matches.size() != energies.size() )
      continue;
    
    double dist = 0.0;
    for( size_t i = 0; i < xray_matches.size(); ++i )
    {
      const double energy = energies[i];
      const SandiaDecay::EnergyIntensityPair *xray = xray_matches[i];
      
      IsotopeMatch match;
      match.m_distance = fabs(energy - xray->energy);    //sum of distance over all energies searched
      dist += match.m_distance;
      
      match.m_age = 0.0;         //age assumed for listing things
      match.m_branchRatio = xray->intensity;
      
      //Only one of the following will be valid: m_nuclide, m_element, m_reaction
      match.m_nuclide = NULL;
      match.m_transition = NULL;
      match.m_particle = NULL;
      match.m_element = el;
      match.m_xray = xray;
      match.m_reaction = NULL;
      match.m_displayData[AssumedAge] = PhysicalUnitsLocalized::printToBestTimeUnits( 0.0 );
      match.m_displayData[ParentIsotope] = match.m_element->symbol;
      
      snprintf( buffer, sizeof(buffer), "%.2f", xray->energy );
      match.m_displayData[Energy] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
      match.m_displayData[Distance] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
      match.m_displayData[BranchRatio] = buffer;
      
      nucmatches.push_back( match );
    }//for( const SandiaDecay::EnergyIntensityPair *xray : xray_matches )
    
    nucmatches[0].m_distance = dist;
    
    {//begin code to get profile distance
      vector<SandiaDecay::EnergyRatePair> srcxrays;
      for( const auto &x : xrays )
        srcxrays.emplace_back( x.intensity, x.energy );
      
      nucmatches[0].m_profileDistance
         = IsotopeId::profile_weight( detector_response_function, displayed_measurement,
                          user_peaks, automated_search_peaks, srcxrays,
                          energies, windows, 1.0, 0.0,
                          nucmatches[0].m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope] );
      
      snprintf( buffer, sizeof(buffer), "%.2f", nucmatches[0].m_profileDistance );
      nucmatches[0].m_displayData[ProfileDistance] = buffer;
    }//end code to get profile distance
    
    snprintf( buffer, sizeof(buffer), "%.2f", dist );
    nucmatches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( nucmatches );
  }//for( const SandiaDecay::Element *el : elements )
}//void xraysWithAllEnergies(...)



void IsotopeSearchByEnergyModel::reactionsWithAllEnergies(
                                                          const std::vector<double> &energies,
                                                          const std::vector<double> &windows,
                                                          const std::shared_ptr<const DetectorPeakResponse> detector_response_function,
                                                          const std::shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                                                          const std::vector<std::shared_ptr<const PeakDef>> &user_peaks,
                                                          const std::vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                                                          const std::vector<const ReactionGamma::Reaction *> &limited_reactions,
                                                          SearchResults &answer )
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  const ReactionGamma *db = ReactionGammaServer::database();
  set<const ReactionGamma::Reaction *> reactions;
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    const float energy = static_cast<float>(energies[i]);
    const float de = static_cast<float>( fabs( windows[i] ) );
    vector<const ReactionGamma::Reaction *> thesereaction;
    db->reactions( energy-de, energy+de, thesereaction );
    
    if( i == 0 )
    {
      reactions.insert( thesereaction.begin(), thesereaction.end() );
    }else
    {
      set<const ReactionGamma::Reaction *> surviving_reactions;
      
      for( const ReactionGamma::Reaction *rctn : thesereaction )
      {
        if( reactions.count( rctn ) )
          surviving_reactions.insert( rctn );
      }//for( const ReactionGamma::Reaction *rctn : thesereaction )
      
      reactions.swap( surviving_reactions );
      
      if( reactions.empty() )
        return;
    }//if( i == 0 ) / else
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  
  if( !limited_reactions.empty() )
  {
    set<const ReactionGamma::Reaction *> surviving_reactions;
    for( const ReactionGamma::Reaction *rctn : reactions )
    {
      auto pos = std::find(begin(limited_reactions), end(limited_reactions), rctn);
      if( pos != end(limited_reactions) )
        surviving_reactions.insert(rctn);
    }
    reactions.swap( surviving_reactions );
  }//if( !limited_reactions.empty() )
  
  
  for( const ReactionGamma::Reaction *rctn : reactions )
  {
    double dist = 0.0;
    vector<IsotopeMatch> matches;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      double smallestDelta = 999999999.9;
      ReactionGamma::Reaction::EnergyYield nearesteA;
      
      for( const ReactionGamma::Reaction::EnergyYield &ea : rctn->gammas )
      {
        const double delta = fabs( ea.energy - energy );
        if( delta < smallestDelta  )
        {
          nearesteA = ea;
          smallestDelta = delta;
        }//if( delta < nearestEnergy  )
      }//for( const ReactionGamma::EnergyAbundance &ea : rctn->gammas )
      
      IsotopeMatch match;
      match.m_distance = energy - nearesteA.energy;
      dist += fabs(match.m_distance);
      match.m_age = 0.0;
      match.m_reaction = rctn;
      match.m_branchRatio = nearesteA.abundance;
      match.m_reactionEnergy = nearesteA;
      
      match.m_displayData[AssumedAge] = PhysicalUnitsLocalized::printToBestTimeUnits( 0.0 );
      match.m_displayData[ParentIsotope] = rctn->name();
      
      snprintf( buffer, sizeof(buffer), "%.2f", nearesteA.energy );
      match.m_displayData[Energy] = buffer;
      
      //      snprintf( buffer, sizeof(buffer), "%.2f", match.m_distance );
      //      match.m_displayData[Distance] = buffer;
      
      snprintf( buffer, sizeof(buffer), "%.2g", match.m_branchRatio );
      match.m_displayData[BranchRatio] = buffer;
      
      //      if( rctn->targetNuclide )
      //        match.m_displayData[SpecificIsotope] = rctn->targetNuclide->symbol;
      
      matches.push_back( match );
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    matches[0].m_distance = dist;
    
    {//begin code to get profile distance
      vector<SandiaDecay::EnergyRatePair> srcgammas;
      for( const ReactionGamma::Reaction::EnergyYield &ea : rctn->gammas )
        srcgammas.emplace_back( ea.abundance, ea.energy );
      
      matches[0].m_profileDistance
        = IsotopeId::profile_weight( detector_response_function, displayed_measurement,
                       user_peaks, automated_search_peaks, srcgammas,
                       energies, windows, 1.0, 0.0,
                       matches[0].m_displayData[IsotopeSearchByEnergyModel::Column::ParentIsotope] );
      
      snprintf( buffer, sizeof(buffer), "%.2f", matches[0].m_profileDistance );
      matches[0].m_displayData[ProfileDistance] = buffer;
    }//end code to get profile distance
    
    snprintf( buffer, sizeof(buffer), "%.2g", matches[0].m_distance );
    matches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( matches );
  }//for( const ReactionGamma::Reaction *rctn : reactions )
  
}//void reactionsWithAllEnergies(...)


void IsotopeSearchByEnergyModel::alphasWithAllEnergies( const vector<double> &energies,
                              const vector<double> &windows,
                              const shared_ptr<const SpecUtils::Measurement> displayed_measurement,
                              const vector<std::shared_ptr<const PeakDef>> &user_peaks,
                              const vector<std::shared_ptr<const PeakDef>> &automated_search_peaks,
                              const bool no_progeny,
                              const vector<const SandiaDecay::Nuclide *> &nuclides,
                              SearchResults &answer )
{
  alphaOrBetasWithAllEnergies( true, energies, windows, displayed_measurement,
                              user_peaks, automated_search_peaks, no_progeny, nuclides, answer );
}


void IsotopeSearchByEnergyModel::betaEndpointWithAllEnergies( const std::vector<double> &energies,
                                  const std::vector<double> &windows,
                                  const bool no_progeny,
                                  const std::vector<const SandiaDecay::Nuclide *> &nuclides,
                                  SearchResults &answer )
{
  alphaOrBetasWithAllEnergies( false, energies, windows, nullptr, {}, {}, no_progeny, nuclides, answer );
}//void betaEndpointWithAllEnergies(...)


void IsotopeSearchByEnergyModel::clearResults()
{
  beginRemoveRows( WModelIndex(), 0, rowCount()-1 );
  m_matches.clear();
  endRemoveRows();
}//void clearResults();


void IsotopeSearchByEnergyModel::updateSearchResults(
                                                     std::shared_ptr<IsotopeSearchByEnergyModel::SearchWorkingSpace> workingspace )
{
  const vector<double> &energies = workingspace->energies;
  const vector<double> &windows = workingspace->windows;
  
  vector< vector<IsotopeMatch> > &matches = workingspace->matches;
  
  m_sortColumn = workingspace->sortColumn;
  m_sortOrder = workingspace->sortOrder;
  
  if( m_matches.size() )
  {
    beginRemoveRows( WModelIndex(), 0, rowCount()-1 );
    m_matches.clear();
    endRemoveRows();
  }//if( m_matches.size() )
  
  if( matches.size() )
  {
    const int ninsert = static_cast<int>( matches.size() * energies.size() );
    beginInsertRows( WModelIndex(), 0, ninsert - 1 );
    m_windows = windows;
    m_energies = energies;
    m_matches.swap( matches );
    endInsertRows();
  }//if( matches.size() )
  
  // There is an apparent bug in Wt 3.7.1 (at least), that if we search on an energy, then change
  //  the window or something, then scroll down in the results table, we will then just get
  //  the loading indicator, and the rows will never render.
  //  I couldnt easily figure out the underlying cause of this, but calling `reset()` or
  //  layoutChanged() in this next line keeps this from happening; we are re-rendering the whole
  //  table anyway, so I guess this isnt that heavy-handed of a work-around...
  //  \sa PeakModel::setPeakFromSpecMeas, PeakModel::setPeaks, ...
  //reset();
  layoutAboutToBeChanged().emit();
  layoutChanged().emit();
  
  workingspace->searchdoneCallback();
  
  
  if( !workingspace->error_msg.empty() )
  {
    // Probably wont ever get here - but JIC.
    passMessage( workingspace->error_msg, 3 );
  }
  
  
  
  wApp->triggerUpdate();
}//void updateSearchResults()



void IsotopeSearchByEnergyModel::setSearchEnergies(
                                                   std::shared_ptr<SearchWorkingSpace> workingspace,
                                                   const double min_rel_br,
                                                   const double minHalfLife,
                                                   Wt::WFlags<IsotopeSearchByEnergyModel::RadSource> srcs,
                                                   const std::vector<const SandiaDecay::Element *> &elements,
                                                   const std::vector<const SandiaDecay::Nuclide *> &nuclides,
                                                   const std::vector<const ReactionGamma::Reaction *> &specific_reactions,
                                                   const std::string appid,
                                                   boost::function< void(void) > updatefcn )
{
  try
  {
    if( !workingspace )
      throw runtime_error( "setSearchEnergies(...): invalid workingspace" );
    
    const vector<double> &energies = workingspace->energies;
    const vector<double> &windows = workingspace->windows;
    
    vector< vector<IsotopeMatch> > &matches = workingspace->matches;
    matches.clear();
    
    if( energies.size() != windows.size() )
      throw runtime_error( "setSearchEnergies(...): input error" );
    
    if( energies.empty() )
    {
      if( updatefcn )
        WServer::instance()->post( appid, updatefcn );
      return;
    }//if( energies.empty() )
    
    
    using SandiaDecay::Element;
    using SandiaDecay::Nuclide;
    using SandiaDecay::EnergyIntensityPair;
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const vector<const SandiaDecay::Element *> &elements = db->elements();
    
    //Get x-rays with at least one of the energies
    map<const SandiaDecay::Element *, vector<EnergyIntensityPair> > filteredXrays;
    for( const Element *el : elements )
      for( const EnergyIntensityPair &xray : el->xrays )
        filteredXrays[el].push_back( xray );
    
    //Get reactions with at least one of the energies
    const ReactionGamma *rctnDb = ReactionGammaServer::database();
    vector<ReactionGamma::ReactionPhotopeak> reactions;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const float minenergy = static_cast<float>(energies[i] - windows[i]);
      const float maxenergy = static_cast<float>(energies[i] + windows[i]);
      rctnDb->reactions( minenergy, maxenergy, reactions );
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    
    //Time to make all the pairings
    auto &user_peaks = workingspace->user_peaks;
    std::vector<std::shared_ptr<const PeakDef>> auto_peaks = workingspace->automated_search_peaks;
    
    if( auto_peaks.empty() && workingspace->foreground && workingspace->displayed_measurement )
    {
      //iOS/Android may not have auto-search peaks yet.  Also recently loaded
      //  spectra and it looks like sometimes when previous states were loaded
      const auto data = workingspace->displayed_measurement;
      auto userpeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( begin(user_peaks), end(user_peaks) );
      const bool singleThreaded = false;
      const bool isHPGe = workingspace->isHPGe;
      auto_peaks = ExperimentalAutomatedPeakSearch::search_for_peaks( data, workingspace->detector_response_function, userpeaksdeque, singleThreaded, isHPGe );
      
      std::shared_ptr<SpecMeas> foreground = workingspace->foreground;
      const set<int> samplenums = workingspace->foreground_samplenums;
      auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( begin(auto_peaks), end(auto_peaks) );
      
      WServer::instance()->post( appid, std::bind( [foreground,samplenums,autopeaksdeque](){
        foreground->setAutomatedSearchPeaks( samplenums, autopeaksdeque );
      }) );
    }//
    
    const auto &meas = workingspace->displayed_measurement;
    
    auto &drf = workingspace->detector_response_function;
    
    // TODO: if searching for more than one source type (e.g., gamma, xray, reaction) use separate
    //       threads and then combine results
    
    const bool no_progeny = (srcs & kNoNuclideProgeny);
    
    //Nuclides that match all energies
    if( srcs & RadSource::NuclideGammaOrXray )
    {
      // nuclidesWithAllEnergies probably wont throw - it will discard any sub-results that cause
      //  unexpected (and there really are none expected) exceptions.
      NuclideMatches filteredNuclides = filter_nuclides( min_rel_br, minHalfLife, no_progeny, energies, windows );

      if( !nuclides.empty() )
      {
        NuclideMatches valid_matches;
        for( const auto &match : filteredNuclides )
        {
          const SandiaDecay::Nuclide * const nuc = match.first;
          const set<double> &energies = match.second;
          const auto wanted_pos = std::find( begin(nuclides), end(nuclides), nuc );
          if( wanted_pos != end(nuclides) )
            valid_matches[nuc] = energies;
        }
        valid_matches.swap( filteredNuclides );
      }//if( !nuclides.empty() )
      
      nuclidesWithAllEnergies( filteredNuclides, energies, windows, min_rel_br, no_progeny,
                              drf, meas, user_peaks, auto_peaks, matches );
    }//if( srcs & RadSource::NuclideGammaOrXray )
    
    //Get elements with x-rays which match all energies. Note: not affected by `no_progeny`
    if( srcs & RadSource::kFluorescenceXRay )
      xraysWithAllEnergies( energies, windows, drf, meas, user_peaks, auto_peaks, elements, matches );
    
    //Get elements with reactions which match all energies. Note: not affected by `no_progeny`
    if( srcs & RadSource::kReaction )
      reactionsWithAllEnergies( energies, windows, drf, meas, user_peaks, auto_peaks, specific_reactions, matches );
    
    if( srcs & RadSource::kAlpha )
      alphasWithAllEnergies( energies, windows, meas, user_peaks, auto_peaks, no_progeny, nuclides, matches );
    
    if( srcs & RadSource::kBetaEndpoint )
      betaEndpointWithAllEnergies( energies, windows, no_progeny, nuclides, matches );
    
    //Get elements with gamma+xrays which match all energies
    
    //Get elements with xrays+reactions which match all energies
    
    //Get elements with gamma+reactions which match all energies
    
    //sort the data
    sortData( matches, energies, workingspace->sortColumn, workingspace->sortOrder );
  }catch( std::exception &e )
  {
    workingspace->error_msg = "The was an unexpected exception, and search results may not be complete.";
    
    stringstream msg;
    msg << "IsotopeSearchByEnergyModel::setSearchEnergies: caught exception: " << e.what();
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.str().c_str() );
#endif
    
    cerr << msg.str() << endl;
  }//try / catch

  if( updatefcn )
    WServer::instance()->post(  appid, updatefcn );
}//void setSearchEnergies( const vector<double> &energies, const double window )


int IsotopeSearchByEnergyModel::columnCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return NumColumns;
}//int columnCount( const Wt::WModelIndex &parent ) const


int IsotopeSearchByEnergyModel::rowCount( const WModelIndex &parent ) const
{
  if( parent.isValid() )
    return 0;
  return static_cast<int>( m_matches.size() * m_energies.size() );
}//int rowCount( const WModelIndex &parent ) const


WModelIndex IsotopeSearchByEnergyModel::parent( const WModelIndex & ) const
{
  return WModelIndex();
}//WModelIndex parent( const Wt::WModelIndex &index ) const;


const SandiaDecay::Nuclide *IsotopeSearchByEnergyModel::nuclide( const WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_nuclide;
}//const SandiaDecay::Nuclide *nuclide( const Wt::WModelIndex &index ) const


const SandiaDecay::Element *IsotopeSearchByEnergyModel::xrayElement( const WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_element;
}//xrayElement( const Wt::WModelIndex &index ) const


const ReactionGamma::Reaction *IsotopeSearchByEnergyModel::reaction( const WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_reaction;
}//reaction( const Wt::WModelIndex &index ) const


double IsotopeSearchByEnergyModel::assumedAge( const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return -1.0;
  
  return m_matches[matchNum].at(0).m_age;
}//double assumedAge( const Wt::WModelIndex &index ) const;


boost::any IsotopeSearchByEnergyModel::data( const WModelIndex &index, int role ) const
{
  const int row = index.row();
  const Column col = Column( index.column() );
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  const size_t energyNum = static_cast<size_t>( row % m_energies.size() );
  
  
  if( (role!=Wt::DisplayRole) || (matchNum>m_matches.size())
     || (col>=NumColumns) )
    return boost::any();
  
  const vector<IsotopeMatch> &match = m_matches[matchNum];
  const IsotopeMatch &iso = match[energyNum];
  
  
  switch( col )
  {
    case ParentIsotope:
    case Distance:
      if( energyNum )
        return boost::any();
      //fallthrough intentional (my first intentional use in like 5 years)
    case Energy: case BranchRatio: case ProfileDistance:
    case SpecificIsotope: case ParentHalfLife: case AssumedAge:
      return iso.m_displayData[col];
      break;
      
    case NumColumns:
      break;
  }//switch( col )
  
  return boost::any();
}//boost::any data(...)


boost::any IsotopeSearchByEnergyModel::headerData( int section,
                                                  Wt::Orientation orientation,
                                                  int role ) const
{
  if( (orientation == Wt::Orientation::Horizontal) && (role == Wt::ItemDataRole::LevelRole) )
    return 0;
  
  if( role == Wt::ItemDataRole::DisplayRole )
  {
    switch( section )
    {
      case IsotopeSearchByEnergyModel::Column::ParentIsotope:
        return WString::tr("isbem-parent-nuc");
      case IsotopeSearchByEnergyModel::Column::Distance:
        return WString::tr("isbem-diff");
      case IsotopeSearchByEnergyModel::Column::Energy:
      {
        if( InterSpec::instance()->isPhone() )
          return WString::tr("Energy");
        return WString::tr("Energy (keV)");
      }
      case IsotopeSearchByEnergyModel::Column::BranchRatio:
        return WString::tr("isbem-rel-br");
      case IsotopeSearchByEnergyModel::Column::ProfileDistance:
        return WString::tr("isbem-profile");
      case IsotopeSearchByEnergyModel::Column::SpecificIsotope:
        return WString::tr("isbem-decay");
      case IsotopeSearchByEnergyModel::Column::ParentHalfLife:
        return WString::tr("isbem-part-hl");
      case IsotopeSearchByEnergyModel::Column::AssumedAge:
        return WString::tr("isbem-assumed-age");
      case IsotopeSearchByEnergyModel::Column::NumColumns:
        break;
    }//switch( col )
  }else if( role == Wt::ItemDataRole::ToolTipRole )
  {
    switch( section )
    {
      case IsotopeSearchByEnergyModel::Column::ParentIsotope:
        return WString::tr("isbem-tt-parent-nuc");
      case IsotopeSearchByEnergyModel::Column::Distance:
        return WString::tr("isbem-tt-diff");
      case IsotopeSearchByEnergyModel::Column::Energy:
        return WString::tr("isbem-tt-energy");
      case IsotopeSearchByEnergyModel::Column::BranchRatio:
        return WString::tr("isbem-tt-rel-bre");
      case IsotopeSearchByEnergyModel::Column::ProfileDistance:
        return WString::tr("isbem-tt-profile");
      case IsotopeSearchByEnergyModel::Column::SpecificIsotope:
        return boost::any();
      case IsotopeSearchByEnergyModel::Column::ParentHalfLife:
        return WString::tr("isbem-tt-parent-hl");
      case IsotopeSearchByEnergyModel::Column::AssumedAge:
        return WString::tr("isbem-tt-assumed-age");
      case IsotopeSearchByEnergyModel::Column::NumColumns:
        break;
    }//switch( col )
  }//ToolTipRole
  
  return boost::any();
}//boost::any headerData( int section, Orientation orientation, int role ) const


WModelIndex IsotopeSearchByEnergyModel::index( int row, int column,
                                              const WModelIndex &parent ) const
{
  if( (column>=NumColumns) || (row>=rowCount()) )
    return WModelIndex();
  
  //const size_t num_energies = m_energies.size();
  //const size_t match_num = row / num_energies;
  //const size_t sub_val = row % num_energies;
  //assert( match_num < m_matches.size() );
  //assert( sub_val < m_matches[match_num].size() );
  //void *ptr = (void *)&(m_matches[match_num][sub_val]);
  
  void *ptr = nullptr;  // I dont think we ever use this
  return createIndex( row, column, ptr );
}//WModelIndex index( int row, int column, const WModelIndex &parent ) const



void IsotopeSearchByEnergyModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();
  m_sortOrder = order;
  m_sortColumn = IsotopeSearchByEnergyModel::Column(column);
  sortData( m_matches, m_energies, column, order );
  layoutChanged().emit();
}//void sort(...)


WFlags<ItemFlag> IsotopeSearchByEnergyModel::flags( const WModelIndex &index ) const
{
  const Column col = Column( index.column() );
  const int level = (index.row() % m_energies.size());
  
  if( (col==ParentIsotope) && (level==0) )
    return WFlags<ItemFlag>( ItemIsSelectable | ItemIsXHTMLText );
  
  return WFlags<ItemFlag>( ItemIsXHTMLText );
}//WFlags<ItemFlag> flags( const Wt::WModelIndex &index ) const


void IsotopeSearchByEnergyModel::sortData( vector< vector<IsotopeMatch> > &data,
                                          const vector<double> &energies,
                                          int column, Wt::SortOrder order )
{
  std::stable_sort( data.begin(), data.end(), Sorter(energies, column, order) );
}//sortData(...)

