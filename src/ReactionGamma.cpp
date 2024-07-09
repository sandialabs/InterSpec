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
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

string ReactionGammaServer::sm_xmlFileLocation = "data/sandia.reactiongamma.xml";
std::mutex ReactionGammaServer::sm_dataBaseMutex;
std::unique_ptr<ReactionGamma> ReactionGammaServer::sm_dataBase;

/*
void print_reaction( std::string name )
{
  try
  {
    const ReactionGamma *db = ReactionGammaServer::database();
    std::vector<ReactionGamma::EnergyAbundance> g = db->gammas( name );
    cerr << "Reaction " << name << " has " << g.size() << " gammas:" << endl;
    for( size_t i = 0; i < g.size(); ++i )
      cerr << "\t" << g[i].energy << "keV, I=" << g[i].abundance << endl;
  }catch( std::exception &e )
  {
    cerr << e.what() << endl;
  }
}//void print_reaction( std::string name )
*/

void ReactionGammaServer::set_xml_file_location( const std::string &filename )
{
  std::unique_lock<std::mutex> lock( sm_dataBaseMutex );
  if( sm_dataBase )
    throw std::runtime_error( "ReactionGammaServer::set_xml_file_location(...)"
                              " can not be called after ReactionGammaServer::database()"
                              " has been called" );
  sm_xmlFileLocation = filename;
}//void set_xml_file_location()

const ReactionGamma *ReactionGammaServer::database()
{
  //is it thread safe to test if a scoped_ptr is valid?
  //  -if so we could avoid getting the thread lock here, and only get it if
  //   we have to create the ReactionGamma object
  std::unique_lock<std::mutex> lock( sm_dataBaseMutex );

  if( sm_dataBase )
    return sm_dataBase.get();

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  sm_dataBase.reset( new ReactionGamma( sm_xmlFileLocation, db ) );

  return sm_dataBase.get();
}//database()


ReactionGammaServer::ReactionGammaServer()
{
  throw runtime_error( "A ReactionGammaServer object "
                            "shouldnt ever be created" );
}


ReactionGamma::EnergyAbundance::EnergyAbundance()
  : energy( 0.0 ),
    abundance( 0.0 )
{
}

bool ReactionGamma::EnergyAbundance::less_than( const EnergyAbundance &lhs,
                                                const EnergyAbundance &rhs )
{
  return lhs.energy < rhs.energy;
}

ReactionGamma::Reaction::Reaction()
  : type( NumReactionType ),
    targetNuclide( NULL ),
    targetElement( NULL ),
    productNuclide( NULL )
{
}


string ReactionGamma::Reaction::name() const
{
  string answer;

  if( targetElement )
    answer = targetElement->symbol;
  else if( targetNuclide )
    answer = targetNuclide->symbol;

  switch( type )
  {
    case AlphaNeutron:            answer += "(a,n)";        break;
    case NeutronAlpha:            answer += "(n,a)";        break;
    case AlphaProton:             answer += "(a,p)";        break;
    case NeutronCapture:          answer += "(n,g)";        break;
    case NeutronInelasticScatter: answer += "(n,n)";        break; // TODO: Add a prime to the second n, the nuclide suggestion code will then need to be updated - the character to add is probably: "\xE2\x80\xB2", "&prime;", "&#8242;" 
    case AlphaInelasticScatter:   answer += "(a,a)";        break; // TODO: Add a prime to the second a...
    case AnnihilationReaction:    answer =  "Annihilation"; break;
    case NumReactionType:         break;
  }//switch( type )

  return answer;
}//string ReactionGamma::Reaction::name() const

//less_than(...): sorts by atomic number, then by mass number (if available)
//  or else
bool ReactionGamma::Reaction::less_than( const Reaction &lhs,
                                         const Reaction &rhs )
{
  int rhs_an = 0, lhs_an = 0, lhs_mn = 0, rhs_mn = 0;
  if( lhs.targetNuclide )
  {
    lhs_an = lhs.targetNuclide->atomicNumber;
    lhs_mn = lhs.targetNuclide->massNumber;
  }else if( lhs.targetElement )
    lhs_an = lhs.targetElement->atomicNumber;

  if( rhs.targetNuclide )
  {
    rhs_an = rhs.targetNuclide->atomicNumber;
    rhs_mn = rhs.targetNuclide->massNumber;
  }else if( rhs.targetElement )
    rhs_an = rhs.targetElement->atomicNumber;

  if( lhs_an == rhs_an )
    return lhs_mn < rhs_mn;
  return lhs_an < rhs_an;
}//less_than(...)



ReactionGamma::ReactionGamma( const string &sandia_reaction_xml,
                              const SandiaDecay::SandiaDecayDataBase *database )
  : m_decayDatabase( database )
{
  init( sandia_reaction_xml );
}//ReactionGamma constructor


ReactionGamma::~ReactionGamma()
{
  for( ReactionType rt = ReactionType(0);
       rt < NumReactionType;
       rt = ReactionType(rt+1) )
  {
    vector<const Reaction *> &rctns = m_reactions[rt];
    for( size_t i = 0; i < rctns.size(); ++i )
      delete rctns[i];
    rctns.clear();
  }//for( loop over reactions )
}//~ReactionGamma()


ReactionGamma::ReactionParticle ReactionGamma::to_particle( const std::string &label )
{
  if( label=="n" || label=="neutron" || label=="n'" )
    return NeutronParticle;
  else if( label=="p" || label=="proton" )
    return ProtonParticle;
  else if( label=="g" || label=="gamma" )
    return GammaParticle;
  else if( label=="a" || label=="alpha" )
    return AlphaParticle;
  throw runtime_error( "Particle, " + label + ", is an invalid particle type"
                       " (valid types: a,g,p,n,alpha,gamma,proton,neutron)" );
}//ReactionParticle to_particle( const std::string &label )


ReactionGamma::ReactionParticle ReactionGamma::ingoing_particle( ReactionType type )
{
  switch( type )
  {
    case AlphaNeutron:            return AlphaParticle;
    case NeutronAlpha:            return NeutronParticle;
    case AlphaProton:             return AlphaParticle;
    case NeutronCapture:          return NeutronParticle;
    case NeutronInelasticScatter: return NeutronParticle;
    case AlphaInelasticScatter:   return AlphaParticle;
    case AnnihilationReaction:    return GammaParticle;
    case NumReactionType: break;
  }//switch( type )

  throw runtime_error( "ReactionGamma::ingoing_particle(...): invalid input" );
}//ReactionParticle ingoing_particle( ReactionType type )


ReactionGamma::ReactionParticle ReactionGamma::outgoing_particle( ReactionType type )
{
  switch( type )
  {
    case AlphaNeutron:            return NeutronParticle;
    case NeutronAlpha:            return AlphaParticle;
    case AlphaProton:             return ProtonParticle;
    case NeutronCapture:          return GammaParticle;
    case NeutronInelasticScatter: return NeutronParticle;
    case AlphaInelasticScatter:   return AlphaParticle;
    case AnnihilationReaction:    return GammaParticle;
    case NumReactionType: break;
  }//switch( type )

  throw runtime_error( "ReactionGamma::outgoing_particle(...): invalid input" );
}//ReactionParticle outgoing_particle( ReactionType type )

//gammas(...): returns the reactions described by names similar to
//  Fe(n,g), Al(a,p), Al(n,g), Ge(n,n'), O(a,n)
//Note that if you specify only the element, and not specific isotope, then
//  resulting abundances will be normalized according to natural material
//  fractions.
//An exception will be thrown if the input is mal-formed (eg invalid element
//  or isotope, non enclosing parenthesis, or invalid reaction types).
//(a,n) is alphaneutron   : Be9, C13, O17, O18, F19,
//(n,a) is neutronalpha   : B10
//(a,p) is alphaproton    : N14, Al27
//(n,g) is neutroncapture : Al27, Be9, Bi, Ca, C, Cl, Cr, Cu, Fe, Ge, H, Mg, N, Ni, O, Pb, Si, Ta,
//(n,n') is neutroninelasticscatter : Al27, Fe
//(a,a') is alphainelasticscatter : Al(27(a,a'g),
std::string ReactionGamma::gammas( const string &name,
                         vector<ReactionGamma::ReactionPhotopeak> &answer ) const
{
  string reactionNames;
  answer.clear();

  if( SpecUtils::icontains(name,"ann") )
  {
    assert( m_reactions[AnnihilationReaction].size() > 0 );
    const Reaction *rcnt = m_reactions[AnnihilationReaction][0];
    assert( rcnt->gammas.size() > 0 );
    
    ReactionPhotopeak rcnanswer;
    rcnanswer.energy = rcnt->gammas[0].energy;
    rcnanswer.abundance = rcnt->gammas[0].abundance;
    rcnanswer.reaction = rcnt;
    answer.push_back( rcnanswer );
    
    return "Annihilation";
  }//if( SpecUtils::icontains(name,"ann") )
  
  const string::size_type open_paren = name.find( '(' );
  const string::size_type close_paren = name.find( ')', open_paren );
  const string::size_type comma_pos = name.find( ',', open_paren );

  if( open_paren==string::npos || close_paren==string::npos )
    throw runtime_error( "ReactionGamma::gammas(...): must have a opening "
                         "and closing parenthisis" );

  if( comma_pos==string::npos || comma_pos>close_paren )
    throw runtime_error( "ReactionGamma::gammas(...): expected a comma in the"
                         " paranthesis" );

  string label = name.substr( 0, open_paren );
  string first_rctn = name.substr( open_paren+1, comma_pos-open_paren-1 );
  string second_rctn = name.substr( comma_pos+1, close_paren-comma_pos-1 );

  SpecUtils::trim( label );
  SpecUtils::trim( first_rctn );
  SpecUtils::trim( second_rctn );
  SpecUtils::to_lower_ascii( first_rctn );
  SpecUtils::to_lower_ascii( second_rctn );

  const bool is_isotope = (label.find_first_of("0123456789")!=string::npos);

  const SandiaDecay::Element *el = NULL;
  const SandiaDecay::Nuclide *nuc = NULL;

  if( !is_isotope )
  {
    el = m_decayDatabase->element( label );
    if( !el )
      throw runtime_error( "Invalid element: " + label );
  }else
  {
    nuc = m_decayDatabase->nuclide( label );
    if( !nuc )
      throw runtime_error( "Invalid nuclide: " + label );
  }//if( is_isotope )


  const ReactionParticle first_prtcl = to_particle( first_rctn );
  const ReactionParticle second_prtcl = to_particle( second_rctn );

  ReactionType reaction_type;
  if( first_prtcl==AlphaParticle && second_prtcl==NeutronParticle )
    reaction_type = AlphaNeutron;
  else if( first_prtcl==NeutronParticle && second_prtcl==AlphaParticle )
    reaction_type = NeutronAlpha;
  else if( first_prtcl==AlphaParticle && second_prtcl==ProtonParticle )
    reaction_type = AlphaProton;
  else if( first_prtcl==NeutronParticle && second_prtcl==GammaParticle )
    reaction_type = NeutronCapture;
  else if( first_prtcl==NeutronParticle && second_prtcl==NeutronParticle )
    reaction_type = NeutronInelasticScatter;
  else if( first_prtcl==AlphaParticle && second_prtcl==AlphaParticle )
    reaction_type = AlphaInelasticScatter;
  else
    throw runtime_error( "Invalid reaction: " + first_rctn + "->"+second_rctn );

  const SandiaDecay::SandiaDecayDataBase *db = m_decayDatabase;

  const vector<const Reaction *> &reactions = m_reactions[reaction_type];
  for( const Reaction *rctn : reactions )
  {
    const SandiaDecay::Element *rctn_el = rctn->targetElement;
    const SandiaDecay::Nuclide *rctn_nuc = rctn->targetNuclide;

    //user asked for reactions cooresponding to this element or nuclide
    //else if user asked for this element, this reaction is for one its isotope
    //XXX - right now we're reying on the reaction data having at most one
    //      isotope per atom, if this changes
    if( (rctn_el && (el==rctn_el)) || (rctn_nuc && nuc && (nuc==rctn_nuc)) )
    {
      for( const Reaction::EnergyYield &ea : rctn->gammas )
      {
        ReactionPhotopeak rp;
        rp.energy = ea.energy;
        rp.abundance = ea.abundance;
        rp.remark = ea.remark;
        rp.reaction = rctn;
        answer.push_back( rp );
      }//for( const EnergyAbundance &ea : rctn->gammas )

//      answer.insert( answer.end(), rctn->gammas.begin(), rctn->gammas.end() );
      if( !reactionNames.empty() )
        reactionNames += ",";
      reactionNames += rctn->name();
    }else if( rctn_nuc && el && (el==db->element(rctn_nuc->atomicNumber)) )
    {
//      answer.insert( answer.end(), rctn->gammas.begin(), rctn->gammas.end() );
      for( const Reaction::EnergyYield &ea : rctn->gammas )
      {
        ReactionPhotopeak rp;
        rp.energy = ea.energy;
        rp.abundance = ea.abundance;
        rp.remark = ea.remark;
        rp.reaction = rctn;
        answer.push_back( rp );
      }//for( const EnergyAbundance &ea : rctn->gammas )

      if( !reactionNames.empty() )
        reactionNames += ",";
      reactionNames += rctn->name();
    }

    //THere is some unecessary ;ooping here - in principle
    // 'if( answer.empty() ) break;' would be fine _with_current_data_
  }//for( const Reaction *rctn : reactions )

  if( answer.empty() )
    throw runtime_error( "Reaction data doesnt contain: " + name );

  return reactionNames;
}//std::vector<ReactionPhotopeak> gammas( const std::string &name ) const


void ReactionGamma::reactions( const SandiaDecay::Element *el,
                std::vector<const ReactionGamma::Reaction *> &answer ) const
{
  if( !el )
    return;

  for( ReactionType type = ReactionType(0);
       type < NumReactionType;
       type = ReactionType(type+1) )
  {
    for( const Reaction *rctn : m_reactions[type] )
    {
      const SandiaDecay::Nuclide *nuc = rctn->targetNuclide;
      if( (rctn->targetElement==el)
          || (nuc && (nuc->atomicNumber==el->atomicNumber)) )
        answer.push_back( rctn );
    }//for( const Reaction *rctn : rctns )
  }//for( loop over ReactionType )
}//reactions(...)


void ReactionGamma::reactions( const SandiaDecay::Nuclide *nuc,
                std::vector<const ReactionGamma::Reaction *> &answer ) const
{
  if( !nuc )
    return;

  for( ReactionType type = ReactionType(0);
       type < NumReactionType;
       type = ReactionType(type+1) )
  {
    for( const Reaction *rctn : m_reactions[type] )
    {
      if( rctn->targetNuclide == nuc )
        answer.push_back( rctn );
    }//for( const Reaction *rctn : rctns )
  }//for( loop over ReactionType )
}//reactions(...)


const vector<const ReactionGamma::Reaction *> &ReactionGamma::reactions(
                                                       ReactionType type ) const
{
  if( type >= NumReactionType )
    throw runtime_error( "ReactionGamma::reactions(): invalid reaction type" );

  return m_reactions[type];
}//reactions( ReactionType type )


//reactions(...): looks through all reactions and returns results in 'answer'
//  for every reaction with a gamma in the energy range specified.
//  Note that 'answer' is appended to (not replaced).
//  There is no garuntee of sort order of 'answer'.
void ReactionGamma::reactions( float lower_energy,
                               float upper_energy,
                       vector<ReactionGamma::ReactionPhotopeak> &answer ) const
{
  for( ReactionType type = ReactionType(0);
       type < NumReactionType;
       type = ReactionType(type+1) )
    reactions( lower_energy, upper_energy, reactions(type), answer );
}//reactions(...)


void ReactionGamma::reactions( float lower_energy,
                               float upper_energy,
                        vector<const ReactionGamma::Reaction *> &answer ) const
{
  for( ReactionType type = ReactionType(0);
       type < NumReactionType;
       type = ReactionType(type+1) )
  {
    const vector<const Reaction *> &rctns = m_reactions[type];

    for( size_t i = 0; i < rctns.size(); ++i )
    {
      const Reaction *reaction = rctns[i];
      const vector<Reaction::EnergyYield> &gammas = reaction->gammas;
      for( size_t j = 0; j < gammas.size(); ++j )
      {
        const Reaction::EnergyYield &ea = gammas[j];
        if( ea.energy >= lower_energy && ea.energy <= upper_energy )
        {
          answer.push_back( reaction );
          break;
        }//if( ea.energy >= lower_energy && ea.energy <= upper_energy )
      }//for( size_t j = 0; j < gammas.size(); ++j )
    }//for( size_t i = 0; i < answer.size(); ++i )
  }//for( loop over ReactionType )
}//reactions(...)


//reactions(...): filters the 'input' and places results into 'answer'.
//  Note that 'answer' is appended to (not replaced).
//  There is no garuntee of sort order of 'answer'.
void ReactionGamma::reactions( float lower_energy,
                       float upper_energy,
                       const vector<const Reaction *> &input,
                       vector<ReactionGamma::ReactionPhotopeak> &answer )
{
  for( size_t i = 0; i < input.size(); ++i )
  {
    const Reaction *reaction = input[i];

    const vector<Reaction::EnergyYield> &gammas = reaction->gammas;
    for( size_t j = 0; j < gammas.size(); ++j )
    {
      const Reaction::EnergyYield &ea = gammas[j];
      if( ea.energy >= lower_energy && ea.energy <= upper_energy )
      {
        ReactionPhotopeak photopeak;
        photopeak.abundance = ea.abundance;
        photopeak.energy = ea.energy;
        photopeak.remark = ea.remark;
        photopeak.reaction = reaction;
        answer.push_back( photopeak );
      }//if( ea.energy >= lower_energy && ea.energy <= upper_energy )
    }//for( size_t j = 0; j < gammas.size(); ++j )
  }//for( size_t i = 0; i < answer.size(); ++i )
}//void ReactionGamma::reactions(...)


void ReactionGamma::populate_reaction( const rapidxml::xml_node<char> *parent,
                        ReactionType type,
                        vector<const ReactionGamma::Reaction *>  &results )
{
  if( !results.empty() )
    throw runtime_error( "ReactionGamma::populate_reaction(...) expects empty "
                         "results passed in" );

  if( !m_decayDatabase )
    throw runtime_error( "ReactionGamma::populate_reaction(...) invalid nuclide"
                         " database ptr" );

  if( !parent )
  {
    cerr << "populate_reaction(...): Warning, null node passed in for reaction="
         << type << endl;
    return;
  }//if( !parent )

  try
  {
    for( const rapidxml::xml_node<char> *node = parent->first_node( "target", 6 );
        node;
        node = node->next_sibling( "target", 6 ) )
    {
      const rapidxml::xml_node<char> *nuclide_node  = node->first_node( "nuclide", 7 );
      const rapidxml::xml_node<char> *product_node  = node->first_node( "product", 7 );
      const rapidxml::xml_node<char> *energies_node = node->first_node( "energies", 8 );
      const rapidxml::xml_node<char> *yeilds_node   = node->first_node( "yields", 6 );
      const rapidxml::xml_node<char> *remark_node   = node->first_node( "remark", 6 );
      
      Reaction rctn;
      rctn.type = type;

      if( !nuclide_node )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Could not "
                             "find nuclide node for a target" );

      if( !nuclide_node->value() )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Empty "
                             "nuclide node for a target" );

      //not checking product_node, since this is optional
      if( !energies_node )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Could not "
                             "find energies node for a target" );
      if( !yeilds_node )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Could not "
                             "find yields node for a target" );

      if( remark_node && remark_node->value_size() )
        rctn.remark.append( remark_node->value(), remark_node->value_size() );
      
      const SandiaDecay::SandiaDecayDataBase *db = m_decayDatabase;

      const string nuclide_name = nuclide_node->value();
      const bool nucNameHasNumbers = (nuclide_name.find_first_of("0123456789") != string::npos);
      
      rctn.targetNuclide = nucNameHasNumbers ? db->nuclide(nuclide_name) : nullptr;
      if( rctn.targetNuclide )
      {
        rctn.targetElement = db->element( rctn.targetNuclide->atomicNumber );
      }else
      {
        rctn.targetElement = db->element(nuclide_name);
//        if( type == NeutronCapture )
//        check if only one natural nuclide and if so, fill that out
      }//if( the nuclide is listed ) / else

      if( !rctn.targetNuclide && !rctn.targetElement )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Couldn't "
                             "find a isotope or element by the name '" + nuclide_name
                              + "', for a target." );

      if( product_node && product_node->value() && product_node->value_size() )
      {
        const string name = product_node->value();
        if( name.find_first_of("0123456789") != string::npos )
          rctn.productNuclide = db->nuclide( name );
      }

      for( const rapidxml::xml_node<char> *enode = energies_node->first_node( "energy", 6 );
           enode;
           enode = enode->next_sibling( "energy", 6 ) )
      {
        if( enode->value() )
        {
          const string str = enode->value();
          Reaction::EnergyYield val;
          if( !(stringstream(str) >> val.energy) )
            throw runtime_error( "ReactionGamma::populate_reaction(...) "
                                 "Couldnt convert " + str + " to a float" );
          val.energy *= static_cast<float>(PhysicalUnits::keV);
          
          const rapidxml::xml_attribute<char> *remark = enode->first_attribute( "remark", 6 );
          if( remark && remark->value_size() )
            val.remark.append( remark->value(), remark->value_size() );
          
          rctn.gammas.push_back( val );
        }//if( enode->value() )
      }//for( loop over energy nodes )

      size_t yeild_num = 0;
      for( const rapidxml::xml_node<char> *ynode = yeilds_node->first_node( "yield", 5 );
           ynode;
           ynode = ynode->next_sibling( "yield", 5 ) )
      {
        if( ynode->value() )
        {
          if( yeild_num >= rctn.gammas.size() )
            throw runtime_error( "ReactionGamma::populate_reaction(...) Less "
                                 "energies than yields" );

          const string str = ynode->value();
          if( !(stringstream(str) >> rctn.gammas[yeild_num].abundance) )
            throw runtime_error( "ReactionGamma::populate_reaction(...) Couldnt"
                                 " convert '" + str + "' to a float" );
          rctn.gammas[yeild_num].abundance /= 100.0;
          
          const rapidxml::xml_attribute<char> *remark = ynode->first_attribute( "remark", 6 );
          if( remark && remark->value_size() )
          {
            if( !rctn.gammas[yeild_num].remark.empty() )
              rctn.gammas[yeild_num].remark += ". ";
            rctn.gammas[yeild_num].remark.append( remark->value(), remark->value_size() );
          }//if( remark && remark->value_size() )
          
          ++yeild_num;
        }//if( ynode->value() )
      }//for( loop over yeild nodes )

      if( yeild_num != rctn.gammas.size() )
        throw runtime_error( "ReactionGamma::populate_reaction(...) Less yields"
                             " than energies" );

      Reaction *ea = new Reaction();
      (*ea) = rctn;
      results.push_back( ea );
    }//for( loop over targets )
  }catch( exception &e )
  {
    for( size_t i = 0; i < results.size(); ++i )
      delete results[i];
    results.clear();
    throw;
  }//try / catch
}//void populate_reaction(...)




void ReactionGamma::init( const string &input )
{
  using namespace rapidxml;

  vector<char> inputdata;
  SpecUtils::load_file_data( input.c_str(), inputdata );

  xml_document<char> document;
  document.parse<parse_normalize_whitespace | parse_trim_whitespace>( &inputdata.front() );
  const rapidxml::xml_node<char> *doc_node = document.first_node();

  if( !doc_node )
    throw runtime_error( "No document node in " + input );

  const rapidxml::xml_node<char> *alphaneutron_node = doc_node->first_node( "alphaneutron", 12 );
  const rapidxml::xml_node<char> *neutronalpha_node = doc_node->first_node( "neutronalpha", 12 );
  const rapidxml::xml_node<char> *alphaproton_node = doc_node->first_node( "alphaproton", 11 );
  const rapidxml::xml_node<char> *capture_node = doc_node->first_node( "neutroncapture", 14 );
  const rapidxml::xml_node<char> *scatter_node = doc_node->first_node( "neutroninelasticscatter", 23 );
  const rapidxml::xml_node<char> *alpha_scatter_node = doc_node->first_node( "alphainelasticscatter", 21 );

  populate_reaction( alphaneutron_node, AlphaNeutron, m_reactions[AlphaNeutron] );
  populate_reaction( neutronalpha_node, NeutronAlpha, m_reactions[NeutronAlpha] );
  populate_reaction( alphaproton_node, AlphaProton, m_reactions[AlphaProton] );
  populate_reaction( capture_node, NeutronCapture, m_reactions[NeutronCapture] );
  populate_reaction( scatter_node, NeutronInelasticScatter, m_reactions[NeutronInelasticScatter] );
  populate_reaction( alpha_scatter_node, AlphaInelasticScatter, m_reactions[AlphaInelasticScatter] );
  
  Reaction *annrctn = new Reaction();
  annrctn->type = AnnihilationReaction;
  annrctn->gammas.resize( 1 );
  annrctn->gammas.back().energy = 510.99891f;
  annrctn->gammas.back().abundance = 1.0f;
  m_reactions[AnnihilationReaction].push_back( annrctn );
}//void ReactionGamma::init(...)
