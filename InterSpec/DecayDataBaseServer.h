#ifndef DecayDataBaseServer_h
#define DecayDataBaseServer_h
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

#include <mutex>
#include <memory>
#include <string>
#include <vector>

#include "SandiaDecay/SandiaDecay.h"

//Since all operations on the database are const, with the exception of
//  initialization, there is no reason to not share a single copy between all
//  sessions/threads/users, so we'll do this through the following class.
class InterSpec_API DecayDataBaseServer
{
public:
  //database() will throw if it is unable to initialse the database
  //  - otherwise will always return a valid pointer
  static const SandiaDecay::SandiaDecayDataBase *database();

  static void initialize();
  static bool initialized();

  //The below functions must be called before first call to database(), and
  //  if called afterwards with a different location, will throw an exception
  static void setDecayXmlFile( const std::string &path_and_file );
  
  static void setXmlFileDirectory( const std::string &dir );  //assumes files named sandia.decay.xml

private:
  
  static std::string sm_decayXrayXmlLocation; //defaults to ./data/sandia.decay.xml
  
  static std::mutex sm_dataBaseMutex;
  static SandiaDecay::SandiaDecayDataBase sm_dataBase;
};//class DecayDataBaseServer


class EnergyToNuclideServer
{
  public:
  /** A convience data struct to allow associating nuclides with specific
   energies gammas for sorting and matching.
   */
  struct EnergyNuclidePair
  {
    //A simple struct to allow sorting in the m_gammaEnergyToNuclideMap object
    float energy;
    const SandiaDecay::Nuclide *nuclide;
    EnergyNuclidePair( float _energy, const SandiaDecay::Nuclide *_nuclide );
    bool operator<( const EnergyNuclidePair &rhs ) const;
  };//struct NuclideEnergyPair
  
  // - TODO - should consider switching from a vector to a deque
  typedef std::vector<EnergyNuclidePair> EnergyNuclidePairVec;
  
public:
  //matcher() will throw if it is unable to initialse the database
  //  so you may wish to call DecayDataBaseServer::setDecayXmlFile(...) or
  //  DecayDataBaseServer::setDecayXmlFile(...) before calling energyToNuclide()
  //  - otherwise will always return a valid pointer
  static const std::shared_ptr< const EnergyNuclidePairVec > energyToNuclide();

  static bool initialized();

  static void unintialize();
  
  static void setLowerLimits( const double halfLife, const double minRelativeBranchRatio );
  static double minHalfLife();
  static double minBranchingRatio();
  
private:
  static std::mutex sm_mutex;
  static double sm_halfLife;
  static double sm_minRelativeBranchRatio;
  static std::shared_ptr< const EnergyNuclidePairVec > sm_energyToNuclide;


public:
  //initGammaToNuclideMatches() times measured 20120514 on my mid 2011 MacBook Pro:
  //  No cuts:
  //          takes 1.6s (cpu and wall), and creates 79259 matches
  //  Regecting isotopes with have half lives less than 600s:
  //          takes 0.47s  (cpu and wall), and creates 44429 matches
  //  Regecting isotopes with no parents that have half lives less than 600s:
  //          takes 1.3s (cpu and wall), and creates 70451 matches
  //  Regecting isotopes who them and all parents have halfLife less than 600s
  //          takes 0.53s (cpu and wall), and creates 47197 matches
  //  Regecting isotopes who them and all parents have halfLife less than 6000s
  //          takes 0.21s (cpu and wall), and creates 29783 matches
  //  Regecting isotopes who them and all parents have halfLife less than 6000s,
  //          and photopeaks which have less than 0.1% of gammas for that decay:
  //          takes 0.05s (cpu and wall), and creates 13744 matches
  //
  //min_halflife specifies not just the minimum halflife of the nuclide the
  //  transition is from, but of all of its forbeares (eg if a grandparent
  //  nuclide has halfLife > min_halflife, then transition is kept)
  //min_gamma_rel_br ranges from 0 to 1.0 and specifies the minimum relative
  //  branching ratio (e.g., the gammas br, divided by max br of that nuclide)
  //  that should be included.  Note that this applies to an individual decay
  //  of a nuclide, not relative to a parent nuclide - so let more things through
  //  than you will want for a long decay chain, so think of it as a pre-filter.
  static void initGammaToNuclideMatches( const SandiaDecay::SandiaDecayDataBase *database,
                                 EnergyNuclidePairVec &results,
                                 const double min_halflife,
                                 const double min_gamma_rel_br );

  //nuclidesWithGammaInRange(...) returns nuclides with gammas in the specified
  //  range.  Note that the range are inclusive, and it is assumed gammaToNuc
  //  is sorted.
  static std::vector<const SandiaDecay::Nuclide *> nuclidesWithGammaInRange( float lowE,
                                                        float highE,
                                                        const EnergyNuclidePairVec &gammaToNuc );
  
  //nuclidesWithGammaInRange(...) returns nuclides with gammas in the specified
  //  range by looping over all input candidates.
  //  Can be computationally very slow if aging is allowed.
  static std::vector<const SandiaDecay::Nuclide *> nuclidesWithGammaInRange( const float lowE,
                                                        const float highE,
                                                        const std::vector<const SandiaDecay::Nuclide *> &candidates,
                                                        const bool allowaging = false  );
};//class NuclidePeakMatcherServer

#endif //DecayDataBaseServer_h
