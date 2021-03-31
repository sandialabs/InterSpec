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
#include <iostream>
#include <boost/filesystem.hpp>

#include "SpecUtils/Filesystem.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

string DecayDataBaseServer::sm_decayXrayXmlLocation = "data/sandia.decay.xml";

std::mutex DecayDataBaseServer::sm_dataBaseMutex;
SandiaDecay::SandiaDecayDataBase DecayDataBaseServer::sm_dataBase;

std::mutex EnergyToNuclideServer::sm_mutex;
std::shared_ptr<const EnergyToNuclideServer::EnergyNuclidePairVec> EnergyToNuclideServer::sm_energyToNuclide;



const SandiaDecay::SandiaDecayDataBase *DecayDataBaseServer::database()
{
  std::lock_guard<std::mutex> lock( sm_dataBaseMutex );

  if( !sm_dataBase.initialized() )
    sm_dataBase.initialize( sm_decayXrayXmlLocation );
  
  return &sm_dataBase;
}//const SandiaDecayDataBase *database()


void DecayDataBaseServer::initialize()
{
  std::lock_guard<std::mutex> lock( sm_dataBaseMutex );

  if( sm_dataBase.initialized() )
    return;

  try
  {
    sm_dataBase.initialize( sm_decayXrayXmlLocation );
    
    if( !sm_dataBase.xmlContainedDecayXRayInfo() )
      throw std::runtime_error( "InterSpec requires nuclear decay XML file to contain decay x-ray info" );
    
    if( !sm_dataBase.xmlContainedElementalXRayInfo() )
      throw std::runtime_error( "InterSpec requires nuclear decay XML file to contain flouresnce x-ray info" );
  }catch(...)
  {
    sm_dataBase.reset();
    
    cerr << "DecayDataBaseServer::initialize()\n\tError initializing SandiaDecayDataBase!"
         << endl << endl;
  }//try/caatch
}//void initialize()


bool DecayDataBaseServer::initialized()
{
  std::lock_guard<std::mutex> lock( sm_dataBaseMutex );
  return sm_dataBase.initialized();
}//bool initialized()

//The below functions must be called before first call to database(), and
//  if called afterwards with a different location, will throw an exception
void DecayDataBaseServer::setDecayXmlFile( const std::string &path_and_file )
{
  std::lock_guard<std::mutex> lock( sm_dataBaseMutex );
  
  if( sm_dataBase.initialized() && path_and_file!=sm_decayXrayXmlLocation )
    throw runtime_error( "You can not call DecayDataBaseServer::setDecayXmlFile"
                        "(...) after already initializing the nuclide "
                        "database" );
  sm_decayXrayXmlLocation = path_and_file;
}//void setDecayXmlFile( const std::string &path_and_file )


void DecayDataBaseServer::setXmlFileDirectory( const std::string &dir )  //assumes file named sandia.decay.xml
{
  std::lock_guard<std::mutex> lock( sm_dataBaseMutex );

  const string file = SpecUtils::append_path( dir, "sandia.decay.xml" );

  if( sm_dataBase.initialized() )
  {
    try{
      if( boost::filesystem::equivalent( file, sm_decayXrayXmlLocation ) )
        return;
    }catch(...){}

    throw runtime_error( "You can not call "
                         "DecayDataBaseServer::setXmlFileDirectory(...) "
                         "after already initializing the nuclide database" );
  }//if( sm_dataBase.initialized() )

  sm_decayXrayXmlLocation = boost::filesystem::path(file).make_preferred().string();
}//void setXmlFileDirectory( const std::string &dir )



double EnergyToNuclideServer::sm_halfLife = 6000.0*SandiaDecay::second;
double EnergyToNuclideServer::sm_branchRatio = 0.0;

const std::shared_ptr< const EnergyToNuclideServer::EnergyNuclidePairVec > &
                                        EnergyToNuclideServer::energyToNuclide()
{
  std::lock_guard<std::mutex> lock( sm_mutex );

  if( !sm_energyToNuclide )
  {
    auto result = make_shared<EnergyNuclidePairVec>();
    result->reserve( 79264 );
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    EnergyToNuclideServer::initGammaToNuclideMatches( db, *result, sm_halfLife, sm_branchRatio );
    sm_energyToNuclide = result;
  }//if( sm_energyToNuclide.empty() )

  return sm_energyToNuclide;
}//const EnergyNuclidePairVec &energyToNuclide()

bool EnergyToNuclideServer::initialized()
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  return static_cast<bool>( sm_energyToNuclide );
}//bool EnergyToNuclideServer::initialized()

void EnergyToNuclideServer::unintialize()
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  sm_energyToNuclide.reset();
}//void unintialize()

double EnergyToNuclideServer::minHalfLife()
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  return sm_halfLife;
}

double EnergyToNuclideServer::minBranchingRatio()
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  return sm_branchRatio;
}


void EnergyToNuclideServer::setLowerLimits( const double halfLife,
                                            const double branchRatio )
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  if( (sm_branchRatio==branchRatio) && (sm_halfLife==halfLife) )
    return;
  sm_energyToNuclide.reset();
  sm_branchRatio = branchRatio;
  sm_halfLife    = halfLife;
}//setLowerLimits(...)


EnergyToNuclideServer::EnergyNuclidePair::EnergyNuclidePair( float e, const SandiaDecay::Nuclide *n )
: energy( e ), nuclide( n )
{
}

bool EnergyToNuclideServer::EnergyNuclidePair::operator<( const EnergyToNuclideServer::EnergyNuclidePair &rhs ) const
{
  return energy < rhs.energy;
}


std::vector<const SandiaDecay::Nuclide *> EnergyToNuclideServer::nuclidesWithGammaInRange( const float lowE,
                                                      const float highE,
                                                      const std::vector<const SandiaDecay::Nuclide *> &candidates,
                                                      const bool allowaging
                                                      )
{
  using namespace SandiaDecay;
  
  std::vector<const Nuclide *> answer;
  
  if( allowaging )
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      NuclideMixture mix;
      mix.addNuclideByActivity( nuc, 1.0E6 * becquerel );
      const vector<EnergyRatePair> gammas = mix.gammas( nuc->halfLife, NuclideMixture::OrderByAbundance, true );
      
      for( size_t i = 0; i < gammas.size(); ++i )
      {
        const EnergyRatePair &a = gammas[i];
        if( a.energy >= lowE && a.energy <= highE )
          if( find(answer.begin(), answer.end(), nuc) == answer.end() )
            answer.push_back( nuc );
      }//for( size_t i = 0; i < gammas.size(); ++i )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }else
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      const std::vector<const Transition *> &trans = nuc->decaysToChildren;
      for( size_t j = 0; j < trans.size(); ++j )
      {
        const Transition * const tran = trans[j];
        const std::vector<RadParticle> &particles = tran->products;
        for( size_t k = 0; k < particles.size(); ++k )
        {
          const RadParticle &particle = particles[k];
          if( particle.type != GammaParticle )
            continue;
          
          if( particle.energy >= lowE && particle.energy <= highE )
            if( find(answer.begin(),answer.end(),nuc) == answer.end() )
              answer.push_back( nuc );
        }//for( size_t k = 0; k < particles.size(); ++k )
      }//for( size_t j = 0; j < trans.size(); ++j )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }//if( allowaging ) / else
  
  return answer;
}//nuclidesWithGammaInRange(...)



vector<const SandiaDecay::Nuclide *> EnergyToNuclideServer::nuclidesWithGammaInRange( float lowE,
                                                 float highE,
                                                 const EnergyNuclidePairVec &gammaToNuc )
{
  if( highE < lowE )
    swap( highE, lowE );
  const EnergyNuclidePair lowerE( lowE, NULL ), upperE( highE, NULL );
  
  const EnergyNuclidePairVec::const_iterator begin = gammaToNuc.begin();
  const EnergyNuclidePairVec::const_iterator end   = gammaToNuc.end();
  
  if( highE < lowE )
    swap( highE, lowE );
  
  EnergyNuclidePairVec::const_iterator lower, upper;
  lower = lower_bound( begin, end, lowerE );
  upper = upper_bound( begin, end, upperE );
  
  vector<const SandiaDecay::Nuclide *> answer;
  
  for( ; lower != upper; ++lower )
    answer.push_back( lower->nuclide );
  
  return answer;
}//nuclidesWithGammaInRange



void EnergyToNuclideServer::initGammaToNuclideMatches( const SandiaDecay::SandiaDecayDataBase *database,
                               EnergyNuclidePairVec &results,
                               const double min_halflife,
                               const double min_gamma_intensity )
{
  if( !results.empty() )
    results.clear();
  
  vector<const SandiaDecay::Nuclide *>::const_iterator nuclideIter = database->nuclides().begin();
  const vector<const SandiaDecay::Nuclide *>::const_iterator end = database->nuclides().end();
  
  
  for( ; nuclideIter != end; ++nuclideIter )
  {
    const SandiaDecay::Nuclide * const nuclide = (*nuclideIter);
    const vector<const SandiaDecay::Transition *> &transitions = nuclide->decaysToChildren;
    
    if( nuclide->halfLife < min_halflife )
    {
      const vector<const SandiaDecay::Nuclide *> forbears = nuclide->forebearers();
      if( forbears.size() < 2 )
        continue;
      bool keep = false;
      for( size_t i = 0; i < forbears.size(); ++i )
      {
        if( forbears[i]->halfLife > min_halflife )
          keep = true;
      }
      
      if( !keep )
        continue;
    }//if( nuclide->halfLife < min_halflife )
    
    float gamma_br_sum = 0.0;
    if( min_gamma_intensity > 0.0 )
    {
      for( size_t trans = 0; trans < transitions.size(); ++trans )
      {
        const SandiaDecay::Transition *transition = transitions[trans];
        const vector<SandiaDecay::RadParticle> &products = transition->products;
        for( size_t part = 0; part < products.size(); ++part )
        {
          if( products[part].type == SandiaDecay::GammaParticle )
            gamma_br_sum += products[part].intensity * transition->branchRatio;
          else if( products[part].type == SandiaDecay::PositronParticle )
            gamma_br_sum += 2.0f * products[part].intensity * transition->branchRatio;
        }//for( size_t part = 0; part < products.size(); ++part )
      }//for( size_t trans = 0; trans < transitions.size(); ++trans )
    }//if( min_gamma_intensity > 0.0 )
    
    
    for( size_t trans = 0; trans < transitions.size(); ++trans )
    {
      const SandiaDecay::Transition *transition = transitions[trans];
      const vector<SandiaDecay::RadParticle> &products = transition->products;
      for( size_t part = 0; part < products.size(); ++part )
      {
        if( products[part].type == SandiaDecay::GammaParticle || products[part].type == SandiaDecay::PositronParticle )
        {
          if( min_gamma_intensity > 0.0 )
          {
            const float mult = 1.0f + static_cast<float>(products[part].type==SandiaDecay::PositronParticle);
            const float intensity = mult * products[part].intensity * transition->branchRatio / gamma_br_sum;
            if( intensity < min_gamma_intensity )
              continue;
          }//if( min_gamma_intensity > 0.0 )
          
          const double energy = (products[part].type==SandiaDecay::GammaParticle ? products[part].energy : static_cast<float>(510.99891*SandiaDecay::keV) );
          const EnergyNuclidePair enPair( energy, nuclide );
          std::vector<EnergyNuclidePair>::iterator begin, end, pos;
          
          begin = results.begin();
          end = results.end();
          pos = upper_bound(  begin, end, enPair );
          
          results.insert( pos, enPair );
        }//if( products[part].type == GammaParticle )
      }//for( loop over RadParticles, part )
    }//for( loop over transitions, trans )
  }//for( loop over nuclides in database )
  
}//void NuclidePeakMatcher() constructor
