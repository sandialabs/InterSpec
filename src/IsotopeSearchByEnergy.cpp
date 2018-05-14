/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#include <vector>
#include <sstream>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WServer>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WTreeView>
#include <Wt/WIOService>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WHBoxLayout>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "sandia_decay/SandiaDecay.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/PhotopeakLineDisplay.h"
#include "InterSpec/IsotopeSearchByEnergy.h"

#if ( USE_SPECTRUM_CHART_D3 )
#include "InterSpec/D3SpectrumDisplayDiv.h"
#endif

using namespace Wt;
using namespace std;


const int IsotopeSearchByEnergy::sm_xmlSerializationVersion = 0;

namespace
{
  const WString ActiveSearchEnergyClass = "ActiveSearchEnergy";
}


IsotopeSearchByEnergyModel::IsotopeMatch::IsotopeMatch()
: m_distance(0.0), m_age(0.0), m_branchRatio(0.0), m_nuclide(0),
  m_transition(0), m_particle(0), m_sourceGammaType(PeakDef::NormalGamma),
  m_element(0), m_xray(0), m_reaction(0)
{
}


IsotopeSearchByEnergyModel::IsotopeSearchByEnergyModel( Wt::WObject *parent )
: WAbstractItemModel( parent ),
  m_sortColumn( Distance ),
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
                                              const double minBR,
       std::vector< vector<IsotopeSearchByEnergyModel::IsotopeMatch> >  &answer)
{
  if( energies.empty() )
    return;
  
  char buffer[32];
  
  //XXX - currently only taking 'most likely' combination of matching of
  //      energies - should implement getting all permutations!
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  for( const NucToEnergiesMap::value_type &nm : filteredNuclides )
  {
    //chack to see if this nuclide has gammas for each energy, not strictly
    //  necassarry, but probably computationally faster (do we care about this
    //  here though)
    bool hasAll = true;
    vector<size_t> energy_xray_is_for;
    vector<const SandiaDecay::EnergyIntensityPair *> matching_xrays;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      const double de = fabs( windows[i] );
      set<double>::const_iterator lb, up, iter;
      lb = nm.second.lower_bound( energy - de );
      up = nm.second.upper_bound( energy + de );
      bool found = false;
      for( iter = lb; iter != up; ++iter )
        found |= (fabs((*iter)-energy) <= de);
      

      if( !found && (energy < 115.0*PhysicalUnits::keV) )
      {//lets look through the x-rays for this
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
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( !hasAll )
      continue;
    
    double dist = 0.0;
    vector<IsotopeMatch> nucmatches;
    for( size_t i = 0; i < energies.size(); ++i )
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
        match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
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
        
        PeakDef::findNearestPhotopeak( nm.first, energies[i],
                         windows[i], transition, trans_index, sourceGammaType );
        if( !transition && (sourceGammaType!=PeakDef::AnnihilationGamma) )
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
        
        match.m_age = PeakDef::defaultDecayTime( nm.first );
      
        SandiaDecay::NuclideMixture mixture;
        mixture.addNuclide( SandiaDecay::NuclideActivityPair(nm.first,1.0) );
        const vector<SandiaDecay::EnergyRatePair> gammas
               = mixture.gammas( match.m_age,
                                 SandiaDecay::NuclideMixture::OrderByAbundance, true );
      
        double nearestEnergy = 999999.9, nearestAbun = 0.0, maxAbund = -999.9;
        for( const SandiaDecay::EnergyRatePair &aep : gammas )
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
        }//for( const SandiaDecay::AbundanceEnergyPair &aep : gammas )
      
        match.m_branchRatio = nearestAbun / maxAbund;
      
        if( match.m_branchRatio < minBR )
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
        
        stringstream trnsitionstrm;
        if( match.m_transition && match.m_transition->parent && match.m_transition->child )
        {
          trnsitionstrm << match.m_transition->parent->symbol << "&rarr;"
                        << match.m_transition->child->symbol;
        }else if( match.m_transition && match.m_transition->parent )
        {
          using namespace SandiaDecay;
          trnsitionstrm << match.m_transition->mode
                        << " of " << match.m_transition->parent->symbol;
        }else if( !match.m_transition )
        {
          trnsitionstrm << "Annih. Gamma";
        }//if( match.m_transition->parent... ) / else
        
        if( match.m_transition && match.m_particle->type == SandiaDecay::XrayParticle )
          trnsitionstrm << " xray";
        
        match.m_displayData[SpecificIsotope] = trnsitionstrm.str();
      
        if( !i )
        {
          match.m_displayData[ParentHalfLife]
                 = PhysicalUnits::printToBestTimeUnits(match.m_nuclide->halfLife);
          match.m_displayData[AssumedAge]
                 = PhysicalUnits::printToBestTimeUnits(match.m_age);
        }//if( !i )
      
        dist += match.m_distance;
        nucmatches.push_back( match );
      }//if( xraypos != energy_xray_is_for.end() ) / else
    }//for( size_t i = 0; i < energies.size(); ++i )
    
    if( nucmatches.size() != energies.size() )
      continue;
    
    if( !nucmatches[0].m_nuclide )
    {
      for( IsotopeMatch &match : nucmatches )
      {
        if( match.m_nuclide )
        {
          nucmatches[0].m_nuclide = match.m_nuclide;
          nucmatches[0].m_displayData[ParentIsotope] = match.m_nuclide->symbol;
          nucmatches[0].m_displayData[Energy]
                    = nucmatches[0].m_displayData[Energy].narrow() + " (xray)";
          break;
        }//if( match.m_nuclide )
      }//for( IsotopeMatch &match : nucmatches )
    }//if( !nucmatches[0].m_nuclide )
    
    nucmatches[0].m_distance = dist;
    
    snprintf( buffer, sizeof(buffer), "%.2f", dist );
    nucmatches[0].m_displayData[Distance] = buffer;
    answer.push_back( nucmatches );
  }//for( const NuclideMatches::value_type &nm : filteredNuclides )
}//void nuclidesWithAllEnergies


void IsotopeSearchByEnergyModel::xraysWithAllEnergies(
                                 const std::vector<double> &energies,
                                 const std::vector<double> &windows,
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
      match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
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
    
    snprintf( buffer, sizeof(buffer), "%.2f", dist );
    nucmatches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( nucmatches );
  }//for( const SandiaDecay::Element *el : elements )
}//void xraysWithAllEnergies(...)



void IsotopeSearchByEnergyModel::reactionsWithAllEnergies(
                                            const std::vector<double> &energies,
                                            const std::vector<double> &windows,
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
  

  for( const ReactionGamma::Reaction *rctn : reactions )
  {
    double dist = 0.0;
    vector<IsotopeMatch> matches;
    for( size_t i = 0; i < energies.size(); ++i )
    {
      const double energy = energies[i];
      double smallestDelta = 999999999.9;
      ReactionGamma::EnergyAbundance nearesteA;
      
      for( const ReactionGamma::EnergyAbundance &ea : rctn->gammas )
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
      
      match.m_displayData[AssumedAge] = PhysicalUnits::printToBestTimeUnits( 0.0 );
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
    
    snprintf( buffer, sizeof(buffer), "%.2g", matches[0].m_distance );
    matches[0].m_displayData[Distance] = buffer;
    
    answer.push_back( matches );
  }//for( const ReactionGamma::Reaction *rctn : reactions )
  
}//void reactionsWithAllEnergies(...)



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
  
  workingspace->searchdoneCallback();
  
  wApp->triggerUpdate();
}//void updateSearchResults()


void IsotopeSearchByEnergyModel::setSearchEnergies(
                            std::shared_ptr<SearchWorkingSpace> workingspace,
                                                   const double minbr,
                                                   const double minHalfLife,
                  Wt::WFlags<IsotopeSearchByEnergyModel::RadSource> srcs,
                                              const std::string appid,
                                              boost::function< void(void) > updatefcn )
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
    WServer::instance()->post(  appid, updatefcn );
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
  
  //Get isotopes with gammas in all ranges...
  EnergyToNuclideServer::setLowerLimits( minHalfLife, minbr );
  std::shared_ptr< const EnergyToNuclideServer::EnergyNuclidePairVec > nucnuc
                                     = EnergyToNuclideServer::energyToNuclide();
  if( !nucnuc )
    throw runtime_error( "Couldnt get EnergyToNuclideServer" );
  
  typedef map<const Nuclide *, set<double> > NuclideMatches;
  NuclideMatches filteredNuclides;
  for( size_t i = 0; i < energies.size(); ++i )  
  {
    const float minenergy = static_cast<float>(energies[i] - windows[i]);
    const float maxenergy = static_cast<float>(energies[i] + windows[i]);
    
    const bool canBeAnnih = (510.99891f>=minenergy && 510.99891f<=maxenergy);
    
    EnergyToNuclideServer::EnergyNuclidePair enPair( minenergy, NULL );
    vector<EnergyToNuclideServer::EnergyNuclidePair>::const_iterator begin, end, pos;
  
    //nucnuc should be sorted by energy (EnergyNuclidePair::operator<)
    begin = lower_bound( nucnuc->begin(), nucnuc->end(), enPair );
    enPair.energy = maxenergy;
    end = upper_bound( begin, nucnuc->end(), enPair );
    for( pos = begin; pos != end; ++pos )
    {
      //nucnuc actually only contians the nuclides that they themselves give
      //  off the requested energies, meaning we have to go through and inspect
      //  all the nuclides that could decay to the nuclide in nucnuc to see
      //  if they are compatible with the search criteria.  For this we
      //  need to
      int stop = 0;
      double transbr = 1.0;
      for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &r : t->products )
        {
          if( r.type==SandiaDecay::GammaParticle && fabs(r.energy - pos->energy) < 0.001 )
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
        }
        
        if( stop==1 )
          break;
      }//for( const SandiaDecay::Transition *t : pos->nuclide->decaysToChildren )
      
      
      const vector<const Nuclide *> forebearers = pos->nuclide->forebearers();
      for( const Nuclide *nuc : forebearers )
      {
        if( nuc->halfLife < minHalfLife )
          continue;
        
        if( minbr<=0.0 || nuc->branchRatioToDecendant(pos->nuclide)*transbr>minbr )
          filteredNuclides[nuc].insert( energies[i] );
      }
    }
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  //Time to make all the pairings
  
  //Nuclides that match all energies
  if( srcs & kGamma )
    nuclidesWithAllEnergies( filteredNuclides, energies, windows, minbr, matches );

  //Get elements with x-rays which match all energies
  if( srcs & kXRay )
    xraysWithAllEnergies( energies, windows, matches );
  
  //Get elements with reactions which match all energies
  if( srcs & kReaction )
    reactionsWithAllEnergies( energies, windows, matches );
  
  //Get elements with gamma+xrays which match all energies
  
  //Get elements with xrays+reactions which match all energies
  
  //Get elements with gamma+reactions which match all energies
  
  //sort the data
  sortData( matches, energies, workingspace->sortColumn, workingspace->sortOrder );
  
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


const SandiaDecay::Nuclide *IsotopeSearchByEnergyModel::nuclide(
                                           const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;

  return m_matches[matchNum].at(0).m_nuclide;
}//const SandiaDecay::Nuclide *nuclide( const Wt::WModelIndex &index ) const


const SandiaDecay::Element *IsotopeSearchByEnergyModel::xrayElement(
                                            const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_element;
}//xrayElement( const Wt::WModelIndex &index ) const


const ReactionGamma::Reaction *IsotopeSearchByEnergyModel::reaction(
                                            const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return NULL;
  
  return m_matches[matchNum].at(0).m_reaction;
}//reaction( const Wt::WModelIndex &index ) const


double IsotopeSearchByEnergyModel::assumedAge(
                                            const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const size_t matchNum = static_cast<size_t>( row / m_energies.size() );
  
  if( matchNum>m_matches.size() )
    return -1.0;

  return m_matches[matchNum].at(0).m_age;
}//double assumedAge( const Wt::WModelIndex &index ) const;


boost::any IsotopeSearchByEnergyModel::data( const WModelIndex &index,
                                             int role ) const
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
    case Energy: case BranchRatio:
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
  if( orientation==Wt::Horizontal && role==Wt::LevelRole )
    return 0;
  
  if (role==Wt::DisplayRole)
  {
    switch( section )
    {
      case ParentIsotope:
        return WString("Parent");
      case Distance:
        return WString("Difference");
      case Energy:
        return WString("Energy (keV)");
      case BranchRatio:
        return WString("Branch Ratio");
      case SpecificIsotope:
        return WString("Decay");
      case ParentHalfLife:
        return WString("Parent H.L.");
      case AssumedAge:
        return WString("Assumed Age");
      case NumColumns:
        break;
    }//switch( col )
  }//DisplayRole
  else if (role==Wt::ToolTipRole)
  {
    switch( section )
    {
      case ParentIsotope:
        return WString("Parent nuclide");
      case Distance:
        return WString("Difference between selected nuclide's energy level and searched energy level");
      case Energy:
        return WString("Nuclide energy");
      case BranchRatio:
        return WString("Branching ratio of selected nuclide");
      case SpecificIsotope:
        return boost::any();
      case ParentHalfLife:
        return WString("Parent half life");
      case AssumedAge:
        return WString("Assumed age of nuclide");
      case NumColumns:
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
  
  void *ptr = (void *)&(m_matches[row]);  //ah, whatever
  return createIndex( row, column, ptr );
}//WModelIndex index( int row, int column, const WModelIndex &parent ) const



void IsotopeSearchByEnergyModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();
  m_sortOrder = order;
  m_sortColumn = Column(column);
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
}//namespace

void IsotopeSearchByEnergyModel::sortData( vector< vector<IsotopeMatch> > &data,
                                           const vector<double> &energies,
                                           int column, Wt::SortOrder order )
{
  std::stable_sort( data.begin(), data.end(), Sorter(energies, column, order) );
}//sortData(...)


void IsotopeSearchByEnergy::SearchEnergy::emitRemove()
{
  m_remove.emit();
}


void IsotopeSearchByEnergy::SearchEnergy::emitChanged()
{
  m_changed.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitEnter()
{
  m_enter.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitGotFocus()
{
  m_focus.emit();
}

void IsotopeSearchByEnergy::SearchEnergy::emitAddAnother()
{
  m_addAnother.emit();
}

IsotopeSearchByEnergy::SearchEnergy::SearchEnergy( Wt::WContainerWidget *p )
: WContainerWidget( p ),
  m_removeIcn( 0 ),
  m_addAnotherIcn( 0 ),
  m_energy( 0 ),
  m_window( 0 )
{
  addStyleClass( "SearchEnergy" );
  WHBoxLayout* layout = new WHBoxLayout();

  layout->setContentsMargins(0, 0, 0, 0);
  setLayout(layout);
  WLabel *label = new WLabel( "Energy " );
  layout->addWidget(label);
  label->addStyleClass( "SearchEnergyLabel" );
  m_energy = new WDoubleSpinBox(  );
  layout->addWidget(m_energy);
  m_energy->setMinimum( 0.0 );
  m_energy->setMaximum( 1000000.0 );
  m_energy->setTextSize( 5 );
  m_energy->enterPressed().connect( this, &SearchEnergy::emitEnter );
  
  label = new WLabel( "+/-" );
  layout->addWidget(label);
  label->addStyleClass( "SearchEnergyWindowLabel" );
  m_window = new WDoubleSpinBox(  );
  layout->addWidget(m_window);
  m_window->setMinimum( 0.0 );
  m_window->setMaximum( 1000000.0 );
  m_window->setValue( 10.0 );
  m_window->setTextSize( 7 );
  m_window->enterPressed().connect( this, &SearchEnergy::emitEnter );

//  InterSpecApp *app= dynamic_cast<InterSpecApp *>( wApp );
//  if( !app && app->isMobile())
//  {
//20150123: on android at least, calling setNativeControl() causes
//  a javascript exception (having to do with the validate)
//    m_energy->setNativeControl(true); //mobile should not show spinner
//    m_window->setNativeControl(true); //mobile should not show spinner
//  } //(isMobile())
  
  label = new WLabel( "keV" );
  layout->addWidget(label);
  
  m_energy->valueChanged().connect( this, &SearchEnergy::emitChanged );
  m_window->valueChanged().connect( this, &SearchEnergy::emitChanged );
  m_energy->keyPressed().connect( this, &SearchEnergy::emitChanged );
  
  m_energy->focussed().connect( this, &SearchEnergy::emitGotFocus );
  m_window->focussed().connect( this, &SearchEnergy::emitGotFocus );
  clicked().connect( this, &SearchEnergy::emitGotFocus );
    
  WContainerWidget *adsubdiv = new WContainerWidget(   );

  m_removeIcn = new WText( adsubdiv );
    
  m_removeIcn->setStyleClass( "DeleteIcon" );
  m_removeIcn->addStyleClass( "energyAddDelete" );
  m_removeIcn->setHiddenKeepsGeometry(true);
  m_removeIcn->clicked().connect( this, &SearchEnergy::emitRemove );
  m_removeIcn->clicked().preventPropagation();
  m_addAnotherIcn = new WText( adsubdiv );

  m_addAnotherIcn->setStyleClass( "AddIcon" );
  m_addAnotherIcn->setHiddenKeepsGeometry(true);
  m_addAnotherIcn->addStyleClass( "energyAddDelete" );
  m_addAnotherIcn->clicked().connect( this, &SearchEnergy::emitAddAnother );
  m_addAnotherIcn->clicked().preventPropagation();
  layout->addWidget(adsubdiv);
}//SearchEnergy constructor


void IsotopeSearchByEnergy::SearchEnergy::enableAddAnother()
{
  m_addAnotherIcn->show();
}//void enableAddAnother()


void IsotopeSearchByEnergy::SearchEnergy::disableAddAnother()
{
  m_addAnotherIcn->hide();
}//void disableAddAnother()


void IsotopeSearchByEnergy::SearchEnergy::enableRemove()
{
  m_removeIcn->show();
}


void IsotopeSearchByEnergy::SearchEnergy::disableRemove()
{
  m_removeIcn->hide();
}


double IsotopeSearchByEnergy::SearchEnergy::energy() const
{
  if( m_energy->validate() == WValidator::Valid )
    return m_energy->value();
  return 0.0;
}

void IsotopeSearchByEnergy::SearchEnergy::setEnergy( double energy )
{
  m_energy->setValue( energy );
}

void IsotopeSearchByEnergy::SearchEnergy::setWindow( double window )
{
  window = floor( 100.0*window + 0.5 ) / 100.0;
  m_window->setValue( window );
}

double IsotopeSearchByEnergy::SearchEnergy::window() const
{
  if( m_window->validate() == WValidator::Valid )
    return m_window->value();
  return 0.0;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::enter()
{
  return m_enter;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::changed()
{
  return m_changed;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::gotFocus()
{
  return m_focus;
}

Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::remove()
{
  return m_remove;
}


Wt::Signal<> &IsotopeSearchByEnergy::SearchEnergy::addAnother()
{
  return m_addAnother;
}


IsotopeSearchByEnergy::IsotopeSearchByEnergy( InterSpec *viewer,
#if ( USE_SPECTRUM_CHART_D3 )
                                              D3SpectrumDisplayDiv *chart,
#else
                                              SpectrumDisplayDiv *chart,
#endif
                                              Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_viewer( viewer ),
  m_chart( chart ),
  m_searchEnergies( NULL ),
  m_currentSearch( 0 ),
  m_searching( NULL ),
  m_results( NULL ),
  m_minBranchRatio( NULL ),
  m_minHalfLife( NULL ),
  m_model( NULL ),
  m_gammas( NULL ),
  m_xrays( NULL ),
  m_reactions( NULL ),
  m_nextSearchEnergy( 0 ),
  m_minBr( 0.0 ), m_minHl( 6000.0 * PhysicalUnits::second )
{
  addStyleClass( "IsotopeSearchByEnergy" );
  WContainerWidget *searchEnergiesDiv = new WContainerWidget();
    
  //This is needed to prevent it from collapsing!
  searchEnergiesDiv->setMinimumSize(370, Wt::WLength::Auto);
  
  WGridLayout* searchEnergiesDivLayout = new WGridLayout();
  searchEnergiesDiv->setLayout(searchEnergiesDivLayout);
 
  searchEnergiesDivLayout->setContentsMargins( 0,0,0,0 );
//  searchEnergiesDiv->setInline( true );
//  searchEnergiesDiv->addStyleClass( "SearchEnergiesLeftDiv" );
  
  m_searchEnergies = new WContainerWidget(  );
  m_searchEnergies->setOverflow( WContainerWidget::OverflowAuto );
  searchEnergiesDivLayout->addWidget(m_searchEnergies,0,0);
  searchEnergiesDivLayout->setRowStretch(0, 1);
  searchEnergiesDivLayout->setColumnStretch(0, 1);
//  m_searchEnergies->addStyleClass( "SearchEnergies" );

  WContainerWidget *buttonDiv = new WContainerWidget(  );

//  WCssDecorationStyle boxxyDeco;
//  boxxyDeco.setBackgroundColor( WColor( color, color, color ) );
//  
//  WString msgstr = prefix + body + /*timestamp +*/ suffix;
//  
////  WText *boxxy = new WText( msgstr, Wt::XHTMLUnsafeText );
//  buttonDiv->setDecorationStyle( boxxyDeco );
  
  searchEnergiesDivLayout->addWidget(buttonDiv,1,0);
  WGridLayout *buttonDivLayout = new WGridLayout();
  buttonDiv->setLayout(buttonDivLayout);
  buttonDivLayout->setContentsMargins(3, 3, 3, 3);
//  buttonDiv->addStyleClass( "SearchEnergiesButtonDiv" );
  
//  WContainerWidget *srcDiv = new WContainerWidget( );

//  srcDiv->addStyleClass( "IsoSearchSrcDiv" );
  
  //XXX - the srcDiv and the searchEnergiesDiv can overlap in the case the
  //      window/tool-bar height isnt very much.  Since this will only happen
  //      when a bunch of energies are searched for, I wont worry about it
  m_gammas = new WCheckBox( "Gammas" );
  buttonDivLayout->addWidget(m_gammas,0,0,1,2, AlignMiddle);
//  m_gammas->addStyleClass( "IsoSearchSrcCB" );
  //  m_gammas->setFloatSide(Wt::Right);
  m_gammas->setChecked();
  m_gammas->checked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_gammas->unChecked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );

  m_xrays = new WCheckBox( "X-rays" );
  buttonDivLayout->addWidget(m_xrays,0,2,1,2, AlignMiddle);

//  m_xrays->addStyleClass( "IsoSearchSrcCB" );
  //  m_xrays->setFloatSide(Wt::Right);
  m_xrays->setChecked();
  m_xrays->checked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_xrays->unChecked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  
  m_reactions = new WCheckBox( "Reactions" );
    buttonDivLayout->addWidget(m_reactions,0,4, AlignMiddle);
//  m_reactions->addStyleClass( "IsoSearchSrcCB" );
//  m_reactions->setFloatSide(Wt::Right);
  //  m_reactions->setChecked();
  m_reactions->checked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_reactions->unChecked().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );

  SearchEnergy *enrgy = new SearchEnergy( m_searchEnergies );
  enrgy->enter().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->gotFocus().connect(
                boost::bind( &IsotopeSearchByEnergy::searchEnergyRecievedFocus,
                                        this, enrgy ) );
  enrgy->remove().connect(
                        boost::bind( &IsotopeSearchByEnergy::removeSearchEnergy,
                                      this, enrgy) );
  enrgy->addStyleClass( ActiveSearchEnergyClass );
  enrgy->addAnother().connect( this, &IsotopeSearchByEnergy::addSearchEnergy );
  enrgy->disableRemove();
  
  
  const bool showToolTipInstantly
  = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  WLabel *label = new WLabel( "Min. BR" );
//  HelpSystem::attachToolTipOn( label,"Toggle or type minimum branching ratio.", showToolTipInstantly , HelpSystem::Top);
   buttonDivLayout->addWidget(label,1,0);
//  label->addStyleClass( "MinBrLabel" );
  m_minBranchRatio = new WDoubleSpinBox();
  string tip = "Minimum branching ratio.";
  HelpSystem::attachToolTipOn( m_minBranchRatio, tip, showToolTipInstantly , HelpSystem::Top);
  buttonDivLayout->addWidget(m_minBranchRatio,1,1);

  m_minBranchRatio->setWidth( 35 );
  m_minBranchRatio->setValue( m_minBr );
  m_minBranchRatio->setRange( 0.0, 1.0 );
  m_minBranchRatio->setSingleStep( 0.1 );
  label = new WLabel( "Min. HL" );
 
  tip = "Minimum half life of nuclides to be searched.<br>"
    "<div>Age can be specified using a combination of time units, "
    "similar to '<b>5.3y 8d 22m</b>'.</div>"
    "<div>"
    "Acceptible time units: <b>year</b>, <b>yr</b>, <b>y</b>, <b>day</b>, <b>d</b>, <b>hrs</b>, <b>hour</b>, <b>h</b>, <b>minute</b>, "
    "<b>min</b>, <b>m</b>, <b>second</b>, <b>s</b>, <b>ms</b>, <b>microseconds</b>, <b>us</b>, <b>nanoseconds</b>, <b>ns</b>, or "
    "you can specify time period by <b>hh:mm:ss</b>. "
    "</div>"
    "<div>"
    "When multiple time periods are "
    "specified, they are summed, e.x. '1y6months 3m' is interpreted as "
    "18 months and 3 minutes"
    "</div>";
//    HelpSystem::attachToolTipOn( label, tip, showToolTipInstantly , HelpSystem::Top);

  buttonDivLayout->addWidget(label,1,2);
  m_minHalfLife = new WLineEdit( "6000 s" );
  HelpSystem::attachToolTipOn( m_minHalfLife, tip, showToolTipInstantly , HelpSystem::Top );

  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_timeDurationRegex, this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_minHalfLife->setValidator(validator);
    
    
    buttonDivLayout->addWidget(m_minHalfLife,1,3);
  m_minHalfLife->setWidth( 55 );
  //XXX should set WRegExpValidator for m_minHalfLife here.

  m_minBranchRatio->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->valueChanged().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->enterPressed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minBranchRatio->blurred().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );

  m_minHalfLife->changed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minHalfLife->enterPressed().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );
  m_minHalfLife->blurred().connect( this, &IsotopeSearchByEnergy::minBrOrHlChanged );

  
  m_searching = new WText( "Searching" );
  buttonDivLayout->addWidget( m_searching, 1, 4, AlignCenter | AlignMiddle );
  m_searching->addStyleClass( "SearchingEnergiesTxt" );
  enrgy->changed().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  
  m_model = new IsotopeSearchByEnergyModel( this );
  
  m_results = new RowStretchTreeView();
  m_results->setRootIsDecorated(false); //makes the tree look like a table! :)
  
  m_results->setModel( m_model );
  m_results->addStyleClass( "IsotopeSearchResultTable" );
  m_results->setAlternatingRowColors( true );
  m_results->sortByColumn( IsotopeSearchByEnergyModel::Distance, Wt::AscendingOrder );
  
  
  for( IsotopeSearchByEnergyModel::Column col = IsotopeSearchByEnergyModel::Column(0);
      col < IsotopeSearchByEnergyModel::NumColumns;
      col = IsotopeSearchByEnergyModel::Column(col+1) )
  {
    m_results->setColumnHidden( col, false );
    m_results->setSortingEnabled( col, true );

    switch( col )
    {
      case IsotopeSearchByEnergyModel::ParentIsotope:
        m_results->setColumnWidth( col, WLength(5,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::Energy:
        m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::Distance:
        m_results->setColumnWidth( col, WLength(6,WLength::FontEm) );
//        m_results->setItemDelegateForColumn(col, new CustomAbstractItemDelegate());
      break;
        
      case IsotopeSearchByEnergyModel::BranchRatio:
        m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::SpecificIsotope:
        m_results->setColumnWidth( col, WLength(8,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::ParentHalfLife:
        m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::AssumedAge:
        m_results->setColumnWidth( col, WLength(7,WLength::FontEm) );
      break;
        
      case IsotopeSearchByEnergyModel::NumColumns:
      break;
    }//switch( col )
  }//for( loop over peak columns )
  
  m_results->setSelectionMode( Wt::SingleSelection );
  m_results->setSelectionBehavior( Wt::SelectRows );
  m_results->selectionChanged().connect( this, &IsotopeSearchByEnergy::resultSelectionChanged );
  
  m_searching->hide();
    
  //XXX - I would like to just add searchEnergiesDiv and tableHolder to *this,
  //  however I cant quite get the CSS right for both Webkit and FireFox, so
  //  I am temporarily using a WGridLayout.
//  addWidget( searchEnergiesDiv );
//  addWidget( tableHolder );
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  layout->addWidget( searchEnergiesDiv, 0, 0 );
  layout->addWidget( m_results, 0, 1 );
  layout->setColumnStretch( 1, 1 );
  layout->setHorizontalSpacing( 0 );
  layout->setVerticalSpacing( 0 );
  
  minBrOrHlChanged();
}//IsotopeSearchByEnergy constuctor

WWidget* CustomAbstractItemDelegate::update	(	Wt::WWidget * 	widget,
                                                                         const Wt::WModelIndex & 	index,
                                                 Wt::WFlags< Wt::ViewItemRenderFlag > 	flags
                                                                         )
{
  WText *item;
  
  if (widget) {
    item = dynamic_cast<WText *>(widget);
  } else {
    item = new WText();
    widget = item;
  }
  
  string text = asString(index.data(DisplayRole)).toUTF8();
  if (text.length()!=0)
  {
    const double val = atof( text.c_str() );
    if (val<0.4)
      item->setAttributeValue("style", "background-color:#C4DF9B"); //green
    else  if (val<1.0)
      item->setAttributeValue("style", "background-color:#FDC68A"); //orange
    else
      item->setAttributeValue("style", "background-color:#F49AC2"); //red
  }
  item->setText(text);
  
  return widget;
} //WWidget* Wt::CustomAbstractItemDelegate::update	(	WWidget * 	widget,const WModelIndex & 	index,WFlags< ViewItemRenderFlag > 	flags)

void IsotopeSearchByEnergy::searchEnergyRecievedFocus( SearchEnergy *enrgy )
{
  std::vector<SearchEnergy *> searchv = searches();
  for( size_t index = 0; index < searchv.size(); ++index )
  {
    if( searchv[index] == enrgy )
    {
      m_nextSearchEnergy = index;
      if( !enrgy->hasStyleClass( ActiveSearchEnergyClass ) )
        enrgy->addStyleClass( ActiveSearchEnergyClass );
    }else
    {
      searchv[index]->removeStyleClass( ActiveSearchEnergyClass );
    }
  }//for( WWidget *kid : kids )
}//void searchEnergyRecievedFocus( SearchEnergy *enrgy )


IsotopeSearchByEnergy::SearchEnergy *IsotopeSearchByEnergy::addNewSearchEnergy()
{
  std::vector<SearchEnergy *> searchv = searches();
  for( SearchEnergy *enrgy : searchv )
    enrgy->disableAddAnother();
  
  SearchEnergy *enrgy = new SearchEnergy( m_searchEnergies );
  enrgy->changed().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->enter().connect( boost::bind( &IsotopeSearchByEnergy::startSearch, this, false ) );
  enrgy->remove().connect(
                      boost::bind( &IsotopeSearchByEnergy::removeSearchEnergy,
                                   this, enrgy) );
  enrgy->addAnother().connect( this, &IsotopeSearchByEnergy::addSearchEnergy );
  enrgy->gotFocus().connect(
                boost::bind( &IsotopeSearchByEnergy::searchEnergyRecievedFocus,
                             this, enrgy ) );
//  m_search->changed().connect( this,
//                          &IsotopeSearchByEnergy::loadSearchEnergiesToClient );
//  m_search->enable();
  
  searchv = searches();
  if( searchv.size() > 1 )
    searchv[0]->enableRemove();
  
  for( size_t i = 0; i < searchv.size(); ++i )
  {
    if( searchv[i]->energy() < 0.1 )
    {
      searchEnergyRecievedFocus( searchv[i] );
      break;
    }//if( searchv[i]->energy() < 0.01 )
  }//for( size_t i = 0; i < searchv.size(); ++i )
  
  return enrgy;
}//SearchEnergy *addNewSearchEnergy()


void IsotopeSearchByEnergy::addSearchEnergy()
{
  addNewSearchEnergy();
}//void addSearchEnergy()


vector<IsotopeSearchByEnergy::SearchEnergy *> IsotopeSearchByEnergy::searches()
{
  vector<SearchEnergy *> searchEnergies;
  
  const vector<WWidget *> &children = m_searchEnergies->children();
  for( size_t index = 0; index < children.size(); ++index )
  {
    SearchEnergy *ww = dynamic_cast<SearchEnergy *>( children[index] );
    if( ww )
      searchEnergies.push_back( ww );
  }
  
  return searchEnergies;
}//vector<SearchEnergy *> searches()


void IsotopeSearchByEnergy::loadSearchEnergiesToClient()
{
  CanvasForDragging *can = m_chart->overlayCanvas();
  if( !can )
    return;
  
  string js;
  const vector<SearchEnergy *> searchW = searches();
  
  if( !searchW.empty() )
  {
    js += "[";
    char buffer[120];
    for( size_t index = 0; index < searchW.size(); ++index )
    {
      if( searchW[index]->energy() > 0.1 )
      {
        if( js.size() > 2 )
          js += ",";
        snprintf( buffer, sizeof(buffer), "[%.2f,%.2f]",
                 searchW[index]->energy(), searchW[index]->window() );
        js += buffer;
      }//if( searchW[index]->energy() > 0.1 )
    }//for( size_t index = 0; index < searchW.size(); ++index )
    js += "]";
  }//if( searchW.empty() ) / else
  
  if( js.size() < 2 )
    js = "null";
  doJavaScript( "$('#c"+can->id()+"').data('SearchEnergies'," + js + ");"
                "Wt.WT.DrawGammaLines('c" + can->id() + "',true);" );
}//void loadSearchEnergiesToClient()


IsotopeSearchByEnergyModel::Column IsotopeSearchByEnergyModel::sortColumn() const
{
  return m_sortColumn;
}

Wt::SortOrder IsotopeSearchByEnergyModel::sortOrder() const
{
  return m_sortOrder;
}


void IsotopeSearchByEnergy::clearSearchEnergiesOnClient()
{
  //Calling wApp->doJavaScript(...) rather than this->doJavaScript(...) to
  //  ensure the command is done, even if this widget is deleted
  CanvasForDragging *can = m_chart->overlayCanvas();
  if( can )
    wApp->doJavaScript( "$('#c"+can->id()+"').data('SearchEnergies',null);"
                       "Wt.WT.DrawGammaLines('c" + can->id() + "',true);" );
}//void clearSearchEnergiesOnClient()


void IsotopeSearchByEnergy::setNextSearchEnergy( double energy, double sigma )
{
  const vector<SearchEnergy *> searchW = searches();
  if( searchW.empty() )  //shouldnt ever happen, bu JIC
    return;
  
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  
  
  for( size_t index = 0; index < searchW.size(); ++index )
  {
    searchW[index]->removeStyleClass( ActiveSearchEnergyClass );
    if( fabs(searchW[index]->energy()-energy) < 0.01*PhysicalUnits::keV )
      m_nextSearchEnergy = index;
  }
  
  searchW[m_nextSearchEnergy]->setEnergy( energy );
  if( sigma > 0.0 )
    searchW[m_nextSearchEnergy]->setWindow( sigma );
  
  ++m_nextSearchEnergy;
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  searchW[m_nextSearchEnergy]->addStyleClass( ActiveSearchEnergyClass );
  
  startSearch( false );
}//void setNextSearchEnergy( double energy )


void IsotopeSearchByEnergy::removeSearchEnergy( IsotopeSearchByEnergy::SearchEnergy *energy )
{
  vector<SearchEnergy *> searchW = searches();
  
  if( searchW.size() < 2 )
    return;
  
  vector<SearchEnergy *>::iterator pos
                     = std::find( searchW.begin(), searchW.end(), energy );
  
  if( pos == searchW.end() )  //shouldnt ever happen, but JIC
    return;
  
  const size_t index = pos - searchW.begin();
  if( index < m_nextSearchEnergy )
    --m_nextSearchEnergy;
  
  searchW.erase( pos );
  if( m_nextSearchEnergy >= searchW.size() )
    m_nextSearchEnergy = 0;
  
  for( SearchEnergy *s : searchW )
    s->removeStyleClass( ActiveSearchEnergyClass );
  
  searchW[m_nextSearchEnergy]->addStyleClass( ActiveSearchEnergyClass );
  searchW.back()->enableAddAnother();
  if( searchW.size() == 1 )
    searchW[0]->disableRemove();
  else if( !searchW.empty() )
    searchW[0]->enableRemove();
  
  delete energy;
  startSearch( false );
}//void removeSearchEnergy( SearchEnergy *energy )



void IsotopeSearchByEnergy::minBrOrHlChanged()
{
  if( m_minBranchRatio->validate() == WValidator::Valid )
    m_minBr = m_minBranchRatio->value();
  else
    m_minBranchRatio->setValue( m_minBr );
  
  try
  {
    const string hltxt = m_minHalfLife->valueText().narrow();
    m_minHl = PhysicalUnits::stringToTimeDuration( hltxt );
  }catch(...)
  {
    m_minHalfLife->setText( PhysicalUnits::printToBestTimeUnits(m_minHl,2) );
  }//try / catch

  startSearch( true );
}//void IsotopeSearchByEnergy::minBrOrHlChanged()


void IsotopeSearchByEnergy::resultSelectionChanged()
{
  if( !m_viewer || !m_viewer->isotopeLinesWidget() )
    return;
 
  PhotopeakLineDisplay *display = m_viewer->isotopeLinesWidget();
  
  WModelIndexSet selected = m_results->selectedIndexes();
  if( selected.empty() )
  {
    display->setIsotope( NULL );
    return;
  }//if( selected.empty() )
  

  if( m_model->nuclide( *selected.begin() ) )
  {
    const SandiaDecay::Nuclide *nuc = m_model->nuclide( *selected.begin() );
    const double age = m_model->assumedAge( *selected.begin() );
    
    display->setIsotope( nuc, age );
  }else if( m_model->xrayElement( *selected.begin() ) )
  {
    const SandiaDecay::Element *el = m_model->xrayElement( *selected.begin() );
    
    display->setElement( el );
  }else if( m_model->reaction( *selected.begin() ) )
  {
    const ReactionGamma::Reaction *rctn = m_model->reaction( *selected.begin() );
    
    display->setReaction( rctn );
  }//if( m_model->nuclide( *selected.begin() ) ) / else ...
}//void resultSelectionChanged()


void IsotopeSearchByEnergy::deSerialize( std::string &xml_data,
                                         const bool renderOnChart )
{
  vector<SearchEnergy *> origSearches;
  for( WWidget *w : m_searchEnergies->children() )
  {
    SearchEnergy *ww = dynamic_cast<SearchEnergy *>( w );
    if( ww )
      origSearches.push_back( ww );
  }//for( const WWebWidget *w : children )
  
  for( size_t i = 1; i < origSearches.size(); ++i )
    removeSearchEnergy( origSearches[i] );
  
  try
  {
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace
    | rapidxml::parse_trim_whitespace;
    if( xml_data.size() )
      doc.parse<flags>( &(xml_data[0]) );
    
    rapidxml::xml_attribute<char> *attr;
    rapidxml::xml_node<char> *base_node, *node, *search_nodes, *value_node;
    
    base_node = doc.first_node( "IsotopeSearchByEnergy", 21 );
    if( !base_node )
      throw runtime_error( "Couldnt get base node, IsotopeSearchByEnergy" );
    
    int version = 0;
    attr = base_node->first_attribute( "version", 7 );
    if( !attr || !attr->value()
       || !(stringstream(attr->value()) >> version)
       || (version != sm_xmlSerializationVersion) )
      throw runtime_error( "Mising or invalid NuclideSearchByEnergy version" );

    node = base_node->first_node( "MinBranchRatio", 14 );
    if( !node || !node->value() )
      throw runtime_error( "Missing MinBranchRatio node" );
    m_minBranchRatio->setValueText( node->value() );
    
    node = base_node->first_node( "MinHalfLife", 11 );
    if( !node || !node->value() )
      throw runtime_error( "Missing MinHalfLife node" );
    m_minHalfLife->setValueText( node->value() );
    
    node = base_node->first_node( "NextSearchEnergy", 16 );
    if( !node || !node->value()
        || !(stringstream(node->value()) >> m_nextSearchEnergy) )
      throw runtime_error( "Missing/invalid NextSearchEnergy node" );

    node = base_node->first_node( "IncludeGammas", 13 );
    if( node && node->value() && strlen(node->value()))
      m_gammas->setChecked( (node->value()[0] == '1') );
    else
      throw runtime_error( "Missing/invalid IncludeGammas node" );

    node = base_node->first_node( "IncludeXRays", 12 );
    if( node && node->value() && strlen(node->value()))
      m_xrays->setChecked( (node->value()[0] == '1') );
    else
      throw runtime_error( "Missing/invalid IncludeXRays node" );
    
    node = base_node->first_node( "IncludeReactions", 16 );
    if( node && node->value() && strlen(node->value()))
      m_reactions->setChecked( (node->value()[0] == '1') );
    else
      throw runtime_error( "Missing/invalid IncludeXRays node" );
    
    search_nodes = base_node->first_node( "SearchEnergies", 14 );
    if( !search_nodes )
      throw runtime_error( "Missing SearchEnergies node" );

    int nnode = 0;
    for( node = search_nodes->first_node( "SearchEnergy", 12 );
         node; node = node->next_sibling( "SearchEnergy", 12 ) )
    {
      double energy, window;

      value_node = node->first_node( "Energy", 6 );
      if( !value_node || !value_node->value()
          || !(stringstream(value_node->value()) >> energy) )
        throw runtime_error( "Missing/invalid SearchEnergy energy entry" );
      
      value_node = node->first_node( "Window", 6 );
      if( !value_node || !value_node->value()
         || !(stringstream(value_node->value()) >> window) )
        throw runtime_error( "Missing/invalid SearchEnergy window entry" );
      
      SearchEnergy *searcher = (SearchEnergy *)0;
      if( nnode++ )
      {
        searcher = addNewSearchEnergy();
      }else
      {
        for( WWidget *w : m_searchEnergies->children() )
        {
          searcher = dynamic_cast<SearchEnergy *>( w );
          if( searcher )
            break;
        }//for( WWidget *w : m_searchEnergies->children() )
        
        if( !searcher )
          throw runtime_error( "Serious logic error in "
                               "IsotopeSearchByEnergy::deSerialize(...)" );
      }//if( nnode++ ) / else
      
      searcher->setEnergy( energy );
      searcher->setWindow( window );
    }//for( loop over nodes )
    
    if( renderOnChart )
    {
      const vector<SearchEnergy *> vals = searches();
      if( vals.size() && (fabs(vals[0]->energy())>0.01 || vals.size()>2) )
        startSearch( true );
      loadSearchEnergiesToClient();
    }else
    {
      //This next call appears to be necassary or else search energies will show
      //  but I'm not certain why...
      clearSearchEnergiesOnClient();
    }//if( renderOnChart )
  }catch( std::exception &e )
  {
    cerr << "IsotopeSearchByEnergy::deSerialize(...) caught: " << e.what() << endl;
    stringstream msg;
    msg << "Error opening displayed photopeaks from database for search: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
  }//try / catch
}//void IsotopeSearchByEnergy::deSerialize( std::string &xml_data )


void IsotopeSearchByEnergy::serialize( std::string &xml_data  ) const
{
  rapidxml::xml_document<char> doc;
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *searchEnergiesNode, *searchnode;
  rapidxml::xml_attribute<> *attr;
  
  name = "IsotopeSearchByEnergy";
  base_node = doc.allocate_node( rapidxml::node_element, name );
  doc.append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  value = doc.allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc.allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  name = "MinBranchRatio";
  value = doc.allocate_string( m_minBranchRatio->valueText().toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "MinHalfLife";
  value = doc.allocate_string( m_minHalfLife->valueText().toUTF8().c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "NextSearchEnergy";
  value = doc.allocate_string( std::to_string(m_nextSearchEnergy).c_str() );
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeGammas";
  value = (m_gammas->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeXRays";
  value = (m_xrays->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "IncludeReactions";
  value = (m_reactions->isChecked() ? "1": "0");
  node = doc.allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  name = "SearchEnergies";
  searchEnergiesNode = doc.allocate_node( rapidxml::node_element, name );
  base_node->append_node( searchEnergiesNode );
  
  const vector<WWidget *> children = m_searchEnergies->children();
  
  for( const WWidget *w : children )
  {
    const SearchEnergy *ww = dynamic_cast<const SearchEnergy *>( w );
    if( ww && (ww->energy() > 0.000001) )
    {
      searchnode = doc.allocate_node( rapidxml::node_element, "SearchEnergy" );
      searchEnergiesNode->append_node( searchnode );
      
      value = doc.allocate_string( std::to_string(ww->energy()).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Energy", value );
      searchnode->append_node( node );
      
      value = doc.allocate_string( std::to_string(ww->window()).c_str() );
      node = doc.allocate_node( rapidxml::node_element, "Window", value );
      searchnode->append_node( node );
    }//if( ww && (ww->energy() > 0.000001) )
  }//for( const WWebWidget *w : children )
  
  xml_data.clear();
  rapidxml::print(std::back_inserter(xml_data), doc, 0);
}//void serialize( std::string &xmlOutput ) const


void IsotopeSearchByEnergy::startSearch( const bool refreshBr )
{
  if( refreshBr )
    EnergyToNuclideServer::setLowerLimits( m_minHl, m_minBr );
  
  loadSearchEnergiesToClient();
  
  vector<double> energies, windows;
  const vector<WWidget *> children = m_searchEnergies->children();
  
  for( const WWidget *w : children )
  {
    const SearchEnergy *ww = dynamic_cast<const SearchEnergy *>( w );
    if( ww && (ww->energy() > 0.000001) )
    {
      energies.push_back( ww->energy() );
      windows.push_back( ww->window() );
    }
  }//for( const WWebWidget *w : children )
  
  
  WFlags<IsotopeSearchByEnergyModel::RadSource> srcs;
  
  if( m_gammas->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kGamma;
  if( m_xrays->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kXRay;
  if( m_reactions->isChecked() )
    srcs |= IsotopeSearchByEnergyModel::kReaction;
  
  
  WApplication *app = wApp;
  std::shared_ptr<IsotopeSearchByEnergyModel::SearchWorkingSpace>
          workingspace( new IsotopeSearchByEnergyModel::SearchWorkingSpace() );
  workingspace->energies = energies;
  workingspace->windows = windows;
  workingspace->sortColumn = m_model->sortColumn();
  workingspace->sortOrder = m_model->sortOrder();
  
  
  ++m_currentSearch;
  workingspace->searchdoneCallback = app->bind(
                          boost::bind( &IsotopeSearchByEnergy::hideSearchingTxt,
                                       this, m_currentSearch ) );
  
  //Verified bellow is safe if the WApplication instance is terminated before
  //  search results are completed, as well as if the WApplication isnt
  //  terminated but m_model is deleted.
  boost::function< void(void) > updatefcnt = app->bind( boost::bind(
                            &IsotopeSearchByEnergyModel::updateSearchResults,
                            m_model, workingspace ) );
  boost::function< void(void) > worker = boost::bind(
                                &IsotopeSearchByEnergyModel::setSearchEnergies,
                                workingspace, m_minBr, m_minHl, srcs,
                                app->sessionId(), updatefcnt );
  WServer::instance()->ioService().post( worker );
  
  m_searching->show();
}//void startSearch()


void IsotopeSearchByEnergy::hideSearchingTxt( const int searchNum )
{
  if( searchNum == m_currentSearch )
    m_searching->hide();
}

IsotopeSearchByEnergy::~IsotopeSearchByEnergy()
{
  
}//IsotopeSearchByEnergy destuctor



