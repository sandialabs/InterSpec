#ifndef ReactionGamma_h
#define ReactionGamma_h
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
#include <string>
#include <vector>
#include <memory>

enum ReactionType
{
  AlphaNeutron,    //gammas produced by alpha,n sources; X(a,n).  nucleus + alpha -> (nucleus+2p+1n) + n + gammas
  NeutronAlpha,    //gammas produced by n,alpha sources; X(n,a).  ex. B10 + neutron -> Li7 + alpha + gammas
  AlphaProton,     //gammas produced by alpha,p sources; X(a,p).  ex. N14 + alpha -> O17 + proton + gammas
  NeutronCapture,  //neutron capture reactions; X(a,n).           nucleus + neutron -> (nucleus+1n) + gammas
  NeutronInelasticScatter, //inelastic scatters; X(a,n).          nucleus + neutron -> nucleus + neutron + gammas
  AlphaInelasticScatter,   //inelastic scatters; X(a,a).          nucleus + alpha -> nucleus + alpha + gammas
  AnnihilationReaction,
  NumReactionType
};//enum ReactionType


class ReactionGamma;
namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  class SandiaDecayDataBase;
}//namespace SandiaDecay

namespace rapidxml
{
  template<class Ch> class xml_node;
}//namespace rapidxml


//Since all operations on the database are const, with the exception of
//  initialization, there is no reason to not share a single copy between all
//  sessions/threads/users, so we'll do this through the following class.
//Note
class InterSpec_API ReactionGammaServer
{
public:
  static const ReactionGamma *database();

  /** Sets the XML file name to use; intended to be called near start of program
      execution.  Must be called before ReactionGammaServer::database() is ever
      called, or else an exception will be thrown.
      Default is "data/sandia.reactiongamma.xml"
   */
  static void set_xml_file_location( const std::string &filename );
  
private:
  ReactionGammaServer();

  static std::string sm_xmlFileLocation; //defaults to data/sandia.reactiongamma.xml
  static std::mutex sm_dataBaseMutex;
  static std::unique_ptr<ReactionGamma> sm_dataBase;
};


class ReactionGamma
{
public:

  enum ReactionParticle
  {
    NeutronParticle,
    ProtonParticle,
    GammaParticle,
    AlphaParticle
  };//enum ReactionParticle

  //to_particle(...): accepts p, n, a, g, proton, neutron, alpha, gamma, n'
  //  throws on invalid input; case dependant.
  static ReactionParticle to_particle( const std::string &label );
  static ReactionParticle ingoing_particle( ReactionType type );
  static ReactionParticle outgoing_particle( ReactionType type );

  struct EnergyAbundance
  {
    EnergyAbundance();

    float energy;
    float abundance;

    //less_than(...): sort by energy
    static bool less_than( const EnergyAbundance &lhs,
                           const EnergyAbundance &rhs );
  };//struct EnergyAbundance


  struct Reaction
  {
    // TODO: make constructor private so only ReactionGamma can create it
    Reaction();

    // TODO: pre-populate the name when reacting the Reaction object.
    std::string name() const;

    ReactionType type;
    
    std::string remark;
    
    struct EnergyYield
    {
      float energy;
      float abundance;
      std::string remark; //will be empty if not given
    };//struct EnergyYield
    
    std::vector<EnergyYield> gammas;  //sorted by energy
    const SandiaDecay::Nuclide *targetNuclide;  //may be null if rates given for natural abundance
    const SandiaDecay::Element *targetElement;  //may be null if rate is given for individual nuclide
    const SandiaDecay::Nuclide *productNuclide; //may be null if targetNuclide is null

    //less_than(...): sorts by atomic number, then by mass number (if available)
    //  or else
    static bool less_than( const Reaction &lhs, const Reaction &rhs );
  };//struct Reaction


  struct ReactionPhotopeak
  {
    float energy;    //in units of PhysicalUnits (e.g. keV)
    float abundance; //fractional abundance (e.g. from 0.0 to 1.0)
    const Reaction *reaction;  //only valid as long as the ReactionGamma object
                               //  used to create this is; not owned by this
                               //  ReactionPhotopeak
    std::string remark;
  };//struct ReactionPhotopeak

  //ReactionGamma(...): will throw exception on mis-formed XML or invalid
  //  xml file name or database pointer
  ReactionGamma( const std::string &sandia_reaction_xml,
                 const SandiaDecay::SandiaDecayDataBase *database );
  ~ReactionGamma();

  //gammas(...): fills 'answer' out w/ reactions described by names similar to
  //  Fe(n,g), Al(a,p), Al(n,g), Ge(n,n'), O(a,n), Annihilation (also Ann, Ann.)
  //Returns a comma separated list of the reaction names.
  //Note that if you specify only the element, and not specific isotope, then
  //  you probably would want intensities normalized according to natural
  //  material fractions - however the code does not explicitly due this and
  //  instead relies on the fact data has at most one isotope per element, so
  //  things _currently_ work out okay.
  //An exception will be thrown if the input is ill-formed (eg invalid element
  //  or isotope, non enclosing parenthesis, or invalid reaction types).
  //(a,n) is alphaneutron   : Be9, C13, O17, O18, F19,
  //(n,a) is neutronalpha   : B10
  //(a,p) is alphaproton    : N14, Al27
  //(n,g) is neutroncapture : Al27, Be9, Bi, Ca, C, Cl, Cr, Cu, Fe, Ge, H, Mg, N, Ni, O, Pb, Si, Ta,
  //(n,n') is neutroninelasticscatter : Al27, Fe
  //(a,a') is alphainelasticscatter
  std::string gammas( const std::string &name,
                      std::vector<ReactionPhotopeak> &answer ) const;

  // TODO: need to implement retrieving a Reaction by name
  
  //reactions(...): gives all reactions which initial element is
  //  specified by first argument; note reaction may be defined for specific
  //  isotopes as well.
  //Note that 'answer' is not cleared, but appended to.
  void reactions( const SandiaDecay::Element *el,
                  std::vector<const Reaction *> &answer ) const;

  //reactions(...): gives all reactions which initial nuclide is
  //  specified by first argument
  //Note that 'answer' is not cleared, but appended to.
  void reactions( const SandiaDecay::Nuclide *nuc,
                  std::vector<const Reaction *> &answer ) const;

  //Access reactions of a given type
  const std::vector<const Reaction *> &reactions( ReactionType type ) const;

  //reactions(...): looks through all reactions and returns results in 'answer'
  //  for every reaction with a gamma in the energy range specified.
  //  Note that 'answer' is appended to (not replaced).
  //  There is no garuntee of sort order of 'answer'.
  //  Also note that energy range is inclusive.
  void reactions( float lower_energy,
                  float upper_energy,
                  std::vector<ReactionPhotopeak> &answer ) const;

  //reactions(...): looks through all reactions and returns results in 'answer'
  //  for every reaction with a gamma in the energy range specified.
  //  Note that 'answer' is appended to (not replaced).
  //  There is no garuntee of sort order of 'answer'.
  //  Also note that energy range is inclusive.
  void reactions( float lower_energy,
                  float upper_energy,
                  std::vector<const Reaction *> &answer ) const;

  //reactions(...): filters the 'input' and places results into 'answer'.
  //  Note that 'answer' is appended to (not replaced).
  //  There is no garuntee of sort order of 'answer'.
  //  Also note that energy range is inclusive.
  static void reactions( float lower_energy,
                         float upper_energy,
                         const std::vector<const Reaction *> &input,
                         std::vector<ReactionPhotopeak> &answer );





protected:
  void init( const std::string &sandia_reaction_xml );
  void populate_reaction( const rapidxml::xml_node<char> *node,
                          ReactionType type,
                          std::vector<const Reaction *>  &results );

protected:
  const SandiaDecay::SandiaDecayDataBase *m_decayDatabase;
  std::vector<const Reaction *> m_reactions[NumReactionType];
};//class ReactionGamma

#endif   //ReactionGamma
