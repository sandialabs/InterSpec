/* SandiaDecay: a library that provides nuclear decay info and calculations.
 
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

#include <set>
#include <map>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <locale>
#include <limits>
#include <cstring>
#include <sstream>
#include <fstream>
#include <ctype.h>
#include <float.h>
#include <iostream>
#include <assert.h>
#include <stdexcept>
#include <algorithm>


#include "rapidxml/rapidxml.hpp"

#include "SandiaDecay.h"

// Ignore this warning. It shouldn't be important, and seems to be Windows-only.
#pragma warning(disable:4244) // warning C4244: 'initializing' : conversion from 'std::streamoff' to 'size_t', possible loss of data

using namespace std;

// Workaround for Windows portability. Theoretically stable.
#if ( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) )
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define fpu_error(x) (isinf(x) || isnan(x))

#define snprintf c99_snprintf

namespace
{
//from http://stackoverflow.com/questions/2915672/snprintf-and-visual-studio-2010
inline int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
  int count = -1;
  if( size != 0 )
    count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
  if( count == -1 )
    count = _vscprintf(format, ap);
  return count;
}

inline int c99_snprintf(char* str, size_t size, const char* format, ...)
{
  int count;
  va_list ap;
  va_start( ap, format );
  count = c99_vsnprintf( str, size, format, ap );
  va_end( ap );
  return count;
}
}//namespace
#endif

#if __cplusplus > 199711L
#define IsInf(x) std::isinf(x)
#define IsNan(x) std::isnan(x)
#else
#define IsInf(x) isinf(x)
#define IsNan(x) isnan(x)
#endif



namespace
{
//a struct so that we can sort an array using an extanal array of indexes,
//  with out modifying the data array.
template<class T> struct index_compare_descend
{
  index_compare_descend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] > arr[b];
  }
  const T arr;
};//struct index_compare


//a struct so that we can sort by index
template<class T> struct index_compare_assend
{
  index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }
  const T arr;
};//struct index_compare

bool convert( const std::string& arg, SandiaDecay::ProductType &s )
{
  if( arg == "beta" )                 s = SandiaDecay::BetaParticle;
  else if( arg == "gamma" )           s = SandiaDecay::GammaParticle;
  else if( arg == "alpha" )           s = SandiaDecay::AlphaParticle;
  else if( arg == "positron" )        s = SandiaDecay::PositronParticle;
  else if( arg == "electronCapture" ) s = SandiaDecay::CaptureElectronParticle;
  else if( arg == "xray" )            s = SandiaDecay::XrayParticle;
  else return false;
  return true;
}//bool convert( const std::string& arg, ProductType &s )

SandiaDecay::DecayMode decaymode_fromstr( const std::string &decayMode )
{
  if( decayMode == "b-" )        return SandiaDecay::BetaDecay;
  else if( decayMode == "it" )   return SandiaDecay::IsometricTransitionDecay;
  else if( decayMode == "b-n" )  return SandiaDecay::BetaAndNeutronDecay;
  else if( decayMode == "b+" )   return SandiaDecay::BetaPlusDecay;
  else if( decayMode == "ec" )   return SandiaDecay::ElectronCaptureDecay;
  else if( decayMode == "b-a" )  return SandiaDecay::BetaAndAlphaDecay;
  else if( decayMode == "b-2n" ) return SandiaDecay::BetaAndTwoNeutronDecay;
  else if( decayMode == "ecp" )  return SandiaDecay::ElectronCaptureAndProtonDecay;
  else if( decayMode == "eca" )  return SandiaDecay::ElectronCaptureAndAlphaDecay;
  else if( decayMode == "b+p" )  return SandiaDecay::BetaPlusAndProtonDecay;
  else if( decayMode == "ec2p" ) return SandiaDecay::ElectronCaptureAndTwoProtonDecay;
  else if( decayMode == "b+2p" ) return SandiaDecay::BetaPlusAndTwoProtonDecay;
  else if( decayMode == "b+3p" ) return SandiaDecay::BetaPlusAndThreeProtonDecay;
  else if( decayMode == "b+a" )  return SandiaDecay::BetaPlusAndAlphaDecay;
  else if( decayMode == "2b-" )  return SandiaDecay::DoubleBetaDecay;
  else if( decayMode == "2ec" )  return SandiaDecay::DoubleElectronCaptureDecay;
  else if( decayMode == "a" )    return SandiaDecay::AlphaDecay;
  else if( decayMode == "p" )    return SandiaDecay::ProtonDecay;
  else if( decayMode == "14c" )  return SandiaDecay::Carbon14Decay;
  else if( decayMode == "sf" )   return SandiaDecay::SpontaneousFissionDecay;
  else if( decayMode == "Undefined" ) return SandiaDecay::UndefinedDecay;
  return SandiaDecay::UndefinedDecay;
}//DecayMode decaymode_fromstr( const char *str )
  
SandiaDecay::ForbiddennessType forbiddenness_fromstr( const char *str )
{
  const size_t len = strlen(str);
  if( len == 1 )
  {
    if( str[0]=='1' )
      return SandiaDecay::FirstForbidden;
    if( str[0]=='2' )
      return SandiaDecay::SecondForbidden;
    if( str[0]=='3' )
      return SandiaDecay::ThirdForbidden;
    if( str[0]=='4' )
      return SandiaDecay::FourthForbidden;
  }else if( len==2 && str[1]=='u' )
  {
    if( str[0]=='1' )
      return SandiaDecay::FirstUniqueForbidden;
    if( str[0]=='2' )
      return SandiaDecay::SecondUniqueForbidden;
    if( str[0]=='3' )
      return SandiaDecay::ThirdUniqueForbidden;
  }
  
  return SandiaDecay::NoForbiddenness;
}//ForbiddennessType forbiddenness_fromstr( const char *str )
  
  
SandiaDecay::RadParticle::RevisionType revision_fromstr( const char *str )
{
  const size_t len = strlen(str);
  if( !len )
    return SandiaDecay::RadParticle::NoRevision;
  if( !strcmp(str,"edited") )
    return SandiaDecay::RadParticle::EditedRevision;
  if( !strcmp(str,"inserted") )
    return SandiaDecay::RadParticle::InsertedRevision;
  if( !strcmp(str,"deleted") )
    return SandiaDecay::RadParticle::DeletedRevision;
  return SandiaDecay::RadParticle::NoRevision;
}//RevisionType revision_fromstr( const char *str )

  
bool lessThanBySymbol( const SandiaDecay::Nuclide *lhs, const std::string &rhs )
{
  if( !lhs )
    return false;
  return (lhs->symbol < rhs);
}//bool lessThanBySymbol( const Nuclide *lhs, const std::string &rhs )

  
bool pointerLessThanBySymbol( const SandiaDecay::Nuclide *lhs, const SandiaDecay::Nuclide *rhs )
{
  if( !lhs )
    return false;
  if( !rhs )
    return true;
  return (lhs->symbol < rhs->symbol);
}//bool pointerLessThanBySymbol( const Nuclide *lhs, const Nuclide *lhs )
  
  
vector<const SandiaDecay::Nuclide *>::const_iterator find(
                            const vector<const SandiaDecay::Nuclide *> &data,
                            const std::string &symbol )
{
  vector<const SandiaDecay::Nuclide *>::const_iterator pos
          = lower_bound( data.begin(), data.end(), symbol, &lessThanBySymbol );
    
  if( (pos == data.end()) || ((*pos)->symbol != symbol) )
    return data.end();
    
  return pos;
}//find(...)

  
typedef double DecayCoefType;
struct BatemanWorkingSpace
{
  std::vector< std::vector<const SandiaDecay::Nuclide *> > nuclide_path;
  std::vector< std::vector< std::vector<DecayCoefType> > > decay_coeffs;
};//struct BatemanWorkingSpace
  
//The following nuclides are known to give nan or inf coeficients in
//  BatemanWorkingSpace::decay_coeffs :  Sm136, Pt183, Ir172, Hg179, Eu136
//  (this issue was found 20121011, but not investigated furhter - I think
//  its been fixed (20180402), but I didnt check that)
void calc_bateman_coef( std::vector< vector<DecayCoefType> > &coefs,
                        std::vector<const SandiaDecay::Transition *> decay_path,
                        BatemanWorkingSpace &ws )
{
  const SandiaDecay::Transition *trans = decay_path.back();
  
  const SandiaDecay::Nuclide *child = trans->child;
  assert( child );
  
  const size_t row = coefs.size();
  coefs.push_back( vector<DecayCoefType>(row+1, 0.0) );  //13.5% of time
  
  for( size_t col = 0; col < row; ++col )
  {
    assert( decay_path.at(row-1)->child );
    const double lambda_iminus1 = decay_path.at(row-1)->parent->decayConstant();
    const double lambda_i = decay_path.at(row-1)->child->decayConstant();
    const double lambda_j = decay_path.at(col)->parent->decayConstant();
    const double br = decay_path[row-1]->branchRatio;
    
    coefs[row][col] = br * (lambda_iminus1/(lambda_i - lambda_j)) * coefs[row-1][col];
    coefs[row][row] -= coefs[row][col];
  }//for( loop over 'col' of matrix A )
  
  
  //We have to get rid of any transition without a daughter for efficiency sake
  vector<const SandiaDecay::Transition *> decays;
  decays.reserve( child->decaysToChildren.size() );
  
  for( size_t t = 0; t < child->decaysToChildren.size(); ++t )
    if( child->decaysToChildren[t]->child )
      decays.push_back( child->decaysToChildren[t] );
  
  
  if( decays.empty() )
  {
    const size_t ndecaypath = decay_path.size();
    ws.nuclide_path.push_back( vector<const SandiaDecay::Nuclide *>() );
    vector<const SandiaDecay::Nuclide *> &decay_nuclides = ws.nuclide_path.back();
    decay_nuclides.resize( ndecaypath + 1 );
    for( size_t i = 0; i < ndecaypath; ++i )
      decay_nuclides[i] = decay_path[i]->parent;
    
    decay_nuclides[ndecaypath] = child;
    ws.decay_coeffs.push_back( coefs );  //17.2% of time
    
    coefs.resize( row );
    
    return;
  }//if( decays.empty() )
  
  const size_t ndecays = decays.size();
  for( size_t t = 0; t < ndecays; ++t )
  {
    const size_t ndecaypath = decay_path.size();
    vector<const SandiaDecay::Transition *> new_path( ndecaypath + 1 );
    for( size_t i = 0; i < ndecaypath; ++i )
      new_path[i] = decay_path[i];
    new_path[ndecaypath] = decays[t];
    
    calc_bateman_coef( coefs, new_path, ws );
  }//for( loop over decaysToChildren )
  
  coefs.resize( row );
}//void calc_bateman_coef(...)
}//namespace


namespace SandiaDecay
{

const char *to_str( SandiaDecay::ProductType s )
{
  switch( s )
  {
    case SandiaDecay::BetaParticle:
      return "beta";
    case SandiaDecay::GammaParticle:
      return "gamma";
    case SandiaDecay::AlphaParticle:
      return "alpha";
    case SandiaDecay::PositronParticle:
      return "positron";
    case SandiaDecay::CaptureElectronParticle:
      return "electronCapture";
    case SandiaDecay::XrayParticle:
      return "xray";
  }//switch( s )

  return "InvalidProductType";
}//const char *to_str( ProductType s )

const char *to_str( const DecayMode &s )
{
  switch( s )
  {
    case AlphaDecay:
      return "Alpha";
    case BetaDecay:
      return "Beta";
    case BetaPlusDecay:
      return "Beta+";
    case ProtonDecay:
      return "Proton";
    case IsometricTransitionDecay:
      return "Isometric";
    case BetaAndNeutronDecay:
      return "Beta and Neutron";
    case BetaAndTwoNeutronDecay:
      return "Beta and 2 Neutron";
    case ElectronCaptureDecay:
      return "Electron Capture";
    case ElectronCaptureAndProtonDecay:
      return "Electron Capture and Proton";
    case ElectronCaptureAndAlphaDecay:
      return "Electron Capture and Alpha";
    case ElectronCaptureAndTwoProtonDecay:
      return "Electron Capture and 2 Proton";
    case BetaAndAlphaDecay:
      return "Beta and Alpha";
    case BetaPlusAndProtonDecay:
      return "Beta Plus and Proton";
    case BetaPlusAndTwoProtonDecay:
      return "Beta Plus and 2 Proton";
    case BetaPlusAndThreeProtonDecay:
      return "Beta Plus and 3 Proton";
    case BetaPlusAndAlphaDecay:
      return "Beta Plus and Alpha";
    case DoubleBetaDecay:
      return "2 Beta+";
    case DoubleElectronCaptureDecay:
      return "2 Electron Capture";
    case Carbon14Decay: //Carbon-14 ?
      return "Carbon-14";
    case SpontaneousFissionDecay:   //I think this is spontaneous fission
      return "Spontaneous Fission";
    case ClusterDecay:
      return "Cluster Emission";
    case UndefinedDecay:
      return "Undefined";
    break;
  };//

  return "";
}//const char *to_str( DecayMode s )
  
std::istream& operator>>( std::istream& input, ProductType &s )
{
  std::string arg;
  input >> arg;

  if( !convert( arg, s ) )
    input.setstate( ios::failbit );

  return input;
}//istream & operator>>( istream &, ProductType & )


std::string human_str_summary( const Nuclide &obj )
{
  stringstream ostr;
    
  ostr << obj.symbol << " Atomic Number " << obj.atomicNumber <<", Atomic Mass "
    << obj.massNumber << ", Isomer Number " << obj.isomerNumber << " "
    << obj.atomicMass << " AMU, HalfLife=" << obj.halfLife << " seconds";
  const size_t nParents = obj.decaysFromParents.size();
  if( nParents )
    ostr << "\n  Parent";
  if( nParents > 1 )
    ostr << "s";
  if( nParents )
    ostr << ": ";
    
  for( size_t i = 0; i < nParents; ++i )
  {
    if( i )
      ostr << ", ";
    if( obj.decaysFromParents[i] && obj.decaysFromParents[i]->parent )
      ostr << obj.decaysFromParents[i]->parent->symbol;
  }//for( loop over parents )
    
  const size_t nChilds = obj.decaysToChildren.size();
    
  for( size_t i = 0; i < nChilds; ++i )
    ostr << "\n    " << human_str_summary(*obj.decaysToChildren[i]);
    
  return ostr.str();
}//std::ostream& operator<<( std::ostream& ostr, const Nuclide &obj )


std::string human_str_summary( const Transition &obj )
{
  stringstream ostr;
    
  ostr << "Transition ";
  if( obj.parent )
    ostr << obj.parent->symbol;
  ostr << "&rarr;";
  if( obj.child )
    ostr << obj.child->symbol;
  else
    ostr << "various";
  
  ostr << ": mode=" << obj.mode << " branchRatio=" << obj.branchRatio;
  if( obj.products.size() )
    ostr << "; Products:";
  for( size_t i = 0; i < obj.products.size(); ++i )
    ostr << "\n    " << SandiaDecay::human_str_summary( obj.products[i] );
    
  return ostr.str();
}//ostream& operator<<( ostream& ostr, const Transition &obj )

  
std::string human_str_summary( const RadParticle &obj )
{
  stringstream ostr;
  
  switch( obj.type )
  {
    case BetaParticle:   ostr << "beta:"; break;
    case GammaParticle:  ostr << "gamma:"; break;
    case AlphaParticle:           ostr << "alpha:"; break;
    case PositronParticle:        ostr << "positron:"; break;
    case CaptureElectronParticle: ostr << "electronCapture:"; break;
    case XrayParticle:           ostr << "xray:"; break;
  };//switch( obj.type )
    
  ostr << " energy=" << obj.energy << "keV intensity=" << obj.intensity;
  switch( obj.type )
  {
    case GammaParticle:
      break;
      
    case AlphaParticle:
      ostr << " hinderence=" << obj.hindrance;
      break;
      
    case BetaParticle: case PositronParticle: case CaptureElectronParticle:
      ostr << " forbiddenness=" << obj.forbiddenness << " logFT=" << obj.logFT;
      break;

    case XrayParticle:
      break;
  };//switch( obj.type )
    
  return ostr.str();
}//std::string human_str_summary( const RadParticle &obj )

  
  


RadParticle::RadParticle( const ::rapidxml::xml_node<char> *node )
{
  set( node );
}

RadParticle::RadParticle(ProductType p_type, float p_energy, float p_intensity)
  : type( p_type ),
    energy( p_energy ),
    intensity( p_intensity ),
    hindrance( 0.0f ),
    logFT( 0.0f ),
    forbiddenness( NoForbiddenness ),
    revision( RadParticle::NoRevision )
{
}


void RadParticle::set( const ::rapidxml::xml_node<char> *node )
{
  //Note: assumes rapidxml::parse_no_string_terminators has NOT been used
  using ::rapidxml::internal::compare;
  typedef ::rapidxml::xml_attribute<char> Attribute;

  const std::string name( node->name(), node->name()+node->name_size() );
  if( !convert( name, type ) )
  {
    const string msg = "Failed converting '" + name + "' to ProductType";
    throw runtime_error( msg );
  }//if( !convert( name, type ) )

  const Attribute *energy_att        = node->first_attribute( "energy", 6 );
  const Attribute *intensity_att     = node->first_attribute( "intensity", 9 );
  const Attribute *hindrance_att     = node->first_attribute( "hindrance", 9 );
  const Attribute *forbiddenness_att = node->first_attribute( "forbiddenness", 13 );
  const Attribute *logFT_att         = node->first_attribute( "logFT", 5 );
  const Attribute *revision_att      = node->first_attribute( "revision", 8 );

  if( energy_att ) energy = (float)atof( energy_att->value() );
  else             energy = 0.0;

  if( intensity_att ) intensity = (float)atof( intensity_att->value() );
  else                intensity = 0.0f;

  if( hindrance_att ) hindrance = (float)atof( hindrance_att->value() );
  else                hindrance = 0.0f;

  if( forbiddenness_att )
    forbiddenness = forbiddenness_fromstr( forbiddenness_att->value() );
  else forbiddenness = NoForbiddenness;

  if( logFT_att ) logFT = static_cast<float>( atof( logFT_att->value() ) );
  else logFT = 0.0f;

  if( revision_att )
    revision = revision_fromstr( revision_att->value() );
  else revision = RadParticle::NoRevision;

  // get any coincident gammas
  // this causes memory corruption somehow?
  for( const ::rapidxml::xml_node<char> *coincidence = node->first_node( "coincidentgamma", 15, true );
      coincidence;
      coincidence = coincidence->next_sibling("coincidentgamma", 15, true) )
  {
    const Attribute *cenergy_att    = coincidence->first_attribute( "energy", 6 );
    const Attribute *cintensity_att = coincidence->first_attribute( "intensity", 9 );
    if( cenergy_att && cenergy_att->value()
        && cintensity_att && cintensity_att->value() )
    {
      const float cenergy = (float)atof( cenergy_att->value() );
      const float cintensity = (float)atof( cintensity_att->value() );
      coincidentRadiation.push_back( RadParticle(GammaParticle, cenergy, cintensity) );
//      std::cout << "Coincident e=" << cenergy << " i=" << cintensity << endl;
    }else
    {
      std::cerr << "RadParticle::set(...): Warning - invalid coincident entry, "
                   "please check the XML file for errors" << endl;
    }
  }//for( loop over coincident gammas )
}//void set( ::rapidxml::xml_node<char> *node )






//Transition::Transition()
//{
//  parent = child = NULL;
//  mode = UndefinedDecay;
//  branchRatio = 0.0;
//}//default Transition constructor


Transition::Transition( const ::rapidxml::xml_node<char> *node,
                        const vector<const Nuclide *> &nuclides )
{
  set( node, nuclides );
}


void Transition::set( const ::rapidxml::xml_node<char> *node,
                      const vector<const Nuclide *> &nuclides )
{
  products.clear();

  //Note: assumes ::rapidxml::parse_no_string_terminators has NOT been used
  using ::rapidxml::internal::compare;
  typedef ::rapidxml::xml_attribute<char> Attribute;

  if( !compare( node->name(), node->name_size(), "transition", 10, true ) )
    throw runtime_error( "Transition::set(xml_node *node): node must be \"transition\"" );

  const Attribute *parent_att      = node->first_attribute( "parent", 6 );
  const Attribute *child_att       = node->first_attribute( "child", 5 );
  const Attribute *mode_att        = node->first_attribute( "mode", 4 );
  const Attribute *branchRatio_att = node->first_attribute( "branchRatio", 11 );

  if( parent_att )
  {
    const std::string parentstr = parent_att->value();
    vector<const Nuclide *>::const_iterator iter = find( nuclides, parentstr );
    if( (parentstr != "") && (iter != nuclides.end()) )
      parent = (*iter);
    else parent = NULL;
  }else parent = NULL;

  if( child_att )
  {
    const std::string childstr = child_att->value();
    vector<const Nuclide *>::const_iterator iter = find( nuclides, childstr );
    if( (childstr != "") && (iter != nuclides.end()) )
      child = (*iter);
    else child = NULL;
  }else child = NULL;


  if( mode_att )
    mode = decaymode_fromstr( mode_att->value() );
  else mode = UndefinedDecay;


  if( branchRatio_att ) branchRatio = (float)atof( branchRatio_att->value() );
  else                  branchRatio = (float)0.0;


  for( const ::rapidxml::xml_node<char> *particle = node->first_node();
      particle;
      particle = particle->next_sibling() )
  {
    if( compare( particle->name(), particle->name_size(), "edit", 4, true )
        || compare( particle->name(), particle->name_size(), "insertion", 9, true ) )
      continue;

    try
    {
      RadParticle product( particle );

      if( product.revision != RadParticle::DeletedRevision )
        products.push_back( product );
    }catch( std::exception &e )
    {
      cerr << "Caught: " << e.what()
//           << " on transition from "
//           << (parent ? parent->symbol : string("null"))
//           << " to "
//           << (child ? child->symbol : string("null"))
           << endl;
    }
  }//for( loop over decaysToChildren )

}//void set( ::rapidxml::xml_node<char> *node )






// XXX only kinda tested
bool Nuclide::operator<( const Nuclide &rhs ) const
{
  //is 'this' a daughter of 'rhs' ? If not, is 'this' lighter than 'rhs
  if( massNumber != rhs.massNumber )
    return (massNumber < rhs.massNumber);
  if( atomicNumber == rhs.atomicNumber )
    return (isomerNumber < rhs.isomerNumber);

// XXX
//The below commented out code makes this function not act correctly,
// (eg allow duplicate nuclides in a map<Nuclide,double> object) but I am unsure why...
//
//  const int an_diff = atomicNumber - rhs.atomicNumber;
//  if( abs( an_diff ) > 1 )
//    return (atomicNumber<rhs.atomicNumber);

  const float lhsToRhsBr = this->branchRatioToDecendant( &rhs );
  const float rhsToLhsBr = rhs.branchRatioToDecendant( this );

  if( lhsToRhsBr==0.0 && rhsToLhsBr==0.0 )
    return (atomicMass<rhs.atomicMass);

  return (lhsToRhsBr<rhsToLhsBr);
}//operator<


bool Nuclide::lessThan( const Nuclide *lhs, const Nuclide *rhs )
{
  return ((*lhs) < (*rhs));
}

bool Nuclide::greaterThan( const Nuclide *lhs, const Nuclide *rhs )
{
  return lessThan( rhs, lhs );
}

bool  Nuclide::operator==( const Nuclide &rhs ) const
{
  return ((atomicNumber==rhs.atomicNumber) && (massNumber==rhs.massNumber) && (isomerNumber==rhs.isomerNumber));
}

bool  Nuclide::operator!=( const Nuclide &rhs ) const
{
  return ((atomicNumber!=rhs.atomicNumber) || (massNumber!=rhs.massNumber) || (isomerNumber!=rhs.isomerNumber));
}

size_t Nuclide::memsize() const
{
  size_t size = sizeof(Nuclide);
  size += symbol.capacity()*sizeof(symbol[0]);
  size += decaysToChildren.capacity()*sizeof(decaysToChildren[0]);
  size += decaysFromParents.capacity()*sizeof(decaysFromParents[0]);

  return size;
}//size_t memsize() const


size_t RadParticle::memsize() const
{
  size_t size = sizeof(RadParticle) + coincidentRadiation.capacity()*sizeof(RadParticle);
  
  for( size_t i = 0; i < coincidentRadiation.size(); ++i )
    size += coincidentRadiation[i].memsize();
  
  return size;
}//size_t memsize() const
  
  
size_t Transition::memsize() const
{
  size_t size = sizeof(Transition) + products.capacity()*sizeof(RadParticle);
  
  for( size_t i = 0; i < products.size(); ++i )
    size += (products[i].memsize()-sizeof(products[i]));
  
  return size;
}
  
size_t Element::memsize() const
{
  size_t size = sizeof(Element);
  size += symbol.capacity()*sizeof(symbol[0]);
  size += name.capacity()*sizeof(name[0]);
  size += isotopes.capacity()*sizeof(isotopes[0]);
  size += xrays.capacity()*sizeof(xrays[0]);
  
  return size;
}

size_t SandiaDecayDataBase::memsize() const
{
  size_t size = sizeof(SandiaDecayDataBase);
  
  size += m_nuclides.capacity()*sizeof(m_nuclides[0]);
  for( size_t i = 0; i < m_nuclideStore.size(); ++i )
    size += (m_nuclideStore[i].memsize()-sizeof(m_nuclideStore[i]));
  size += m_nuclideStore.capacity()*sizeof(m_nuclideStore[0]);
  
  size += m_elements.capacity()*sizeof(m_elements[0]);
  for( size_t i = 0; i < m_elements.size(); ++i )
    size += (m_elements[i]->memsize()-sizeof(m_elements[i]));
  size += m_elementStore.capacity()*sizeof(m_elementStore[0]);
  
  for( size_t i = 0; i < m_transitionStore.size(); ++i )
    size += (m_transitionStore[i].memsize()-sizeof(m_transitionStore[i]));
  size += m_transitionStore.capacity()*sizeof(m_transitionStore[0]);
  
  return size;
}//size_t memsize() const
  
  
size_t Nuclide::numDecays() const
{
  return decaysToChildren.size();
}

double Nuclide::decayConstant() const
{
  //0.5 = exp( -decay_const*halfLife )
  //IsInf() dosnt seem to be available in Microsoft VC++
  if( !(halfLife <= DBL_MAX && halfLife >= -DBL_MAX) || (halfLife<=0.0) ) //is this check really neaded?
    return 0.0;

  return 0.693147180559945286 / halfLife;
}//double decayConstant()

double Nuclide::partialDecayConstant( const size_t transition_num ) const
{
  //0.5 = exp( -(decay_const_1 + decay_const_2 + ...) * halfLife )
  const double total_decay_constant = decayConstant();
  return total_decay_constant * decaysToChildren.at(transition_num)->branchRatio;
}//double partialDecayConstant( const size_t transition_num )

float Nuclide::branchingRatio( const size_t transition_num ) const
{
  return decaysToChildren.at( transition_num )->branchRatio;
}//double branchingRatio( const size_t transition_num ) const


const Nuclide *Nuclide::child( const size_t transition_num ) const
{
  return decaysToChildren.at( transition_num )->child;
}



double Nuclide::numAtomsToActivity( const double num_atoms ) const
{
  return num_atoms * decayConstant();
}//double numAtomsToActivity( const double num_atoms ) const

double Nuclide::activityToNumAtoms( const double activity ) const
{
  return activity / decayConstant();
}//double activityToNumAtoms( const double activity ) const


// XXX - Barely tested
// XXX - If there is a cycle: Nuclide_1 --> Nuclide_2 --> Nuclide_1, then
//       this will cause an infinite recurtion
float Nuclide::branchRatioToDecendant( const Nuclide *descendant ) const
{
  if( this == descendant )
    return 1.0f;

  if( massNumber < descendant->massNumber )
    return 0.0f;

  float br_sum = 0.0;
  for( size_t i = 0; i < decaysToChildren.size(); ++i )
  {
    const Nuclide *child = decaysToChildren[i]->child;
    if( child )
    {
      const float br = decaysToChildren[i]->branchRatio;
      br_sum += br * child->branchRatioToDecendant( descendant );
    }//if( child )
  }//for( loop over decay paths )

  return br_sum;
}//double branchRatioToDecendant( const Nuclide *descendant ) const


// XXX - barels tested
// XXX - If there is a cylcle: Nuclide_1 --> Nuclide_2 --> Nuclide_1, then
//       this will cause an infinite recurtion
float Nuclide::branchRatioFromForebear( const Nuclide *ancestor ) const
{
  //I think the below would work, but
  // return ancestor->branchRatioToDecendant( this );

  if( this == ancestor )
    return 1.0f;

  if( massNumber > ancestor->massNumber )
    return 0.0f;

  float br_sum = 0.0;
  for( size_t i = 0; i < decaysFromParents.size(); ++i )
  {
    const Nuclide *parent = decaysFromParents[i]->parent;
    const float br = decaysFromParents[i]->branchRatio;
    br_sum += br * parent->branchRatioFromForebear( ancestor );
  }//for( loop over decay paths )

  return br_sum;
}//double branchRatioFromForebear( const Nuclide *ancestor ) const



vector<const Nuclide *> Nuclide::descendants() const
{
  vector<const Nuclide *> children( 1, this );

  for( size_t decN = 0; decN < decaysToChildren.size(); ++decN )
  {
    const Transition * const transition = decaysToChildren[decN];
    const Nuclide * const child = transition->child;

    if( !child )
      continue;

    const vector<const Nuclide *> results = child->descendants();

    for( size_t resN = 0; resN < results.size(); ++resN )
    {
      const Nuclide * const result = results[resN];
      const vector<const Nuclide *>::iterator begin = children.begin();
      const vector<const Nuclide *>::iterator end = children.end();
      vector<const Nuclide *>::iterator pos;
      pos = lower_bound( begin, end, result, &Nuclide::greaterThan );

      if( (pos==end) || ((*pos)!=result) )
        children.insert( pos, result );
    }//for( size_t resN = 0; resN < results.size(); ++resN )
  }//for( size_t decN = 0; decN < decaysToChildren.size(); ++decN )


  return children;
}//std::vector<const Nuclide *> descendants() const


std::vector<const Nuclide *> Nuclide::forebearers() const
{
  std::vector<const Nuclide *> parents( 1, this );

  //
  for( size_t decN = 0; decN < decaysFromParents.size(); ++decN )
  {
    const Transition * const transition = decaysFromParents[decN];
    const Nuclide * const parent = transition->parent;
    const vector<const Nuclide *> results = parent->forebearers();

    for( size_t resN = 0; resN < results.size(); ++resN )
    {
      const Nuclide * const result = results[resN];
      const vector<const Nuclide *>::iterator begin = parents.begin();
      const vector<const Nuclide *>::iterator end = parents.end();

      vector<const Nuclide *>::iterator pos;
      pos = lower_bound( begin, end, result, &Nuclide::lessThan );

      if( (pos==end) || ((*pos)!=result) )
        parents.insert( pos, result );
    }//for( size_t resN = 0; resN < results.size(); ++resN )
  }//for( size_t decN = 0; decN < decaysToChildren.size(); ++decN )

  return parents;
}//std::vector<const Nuclide *> forebearers() const



double Nuclide::secularEquilibriumHalfLife() const
{
  //20111012 Verified gives correct correct for U235, U238, Pu240, Pu242, Th232
  if( IsInf(halfLife) )
    return numeric_limits<double>::infinity();

  double secHalflife = 0.0;
  for( size_t i = 0; i < decaysToChildren.size(); ++i )
  {
    const Nuclide *child = decaysToChildren[i]->child;
    if( child && (child != this) )
    {
      if( !IsInf(child->halfLife) )
        secHalflife = max( secHalflife,
                         max( child->halfLife,
                              child->secularEquilibriumHalfLife() ) );
    }
  }//for( size_t i = 0; i < decaysToChildren.size(); ++i )

  return secHalflife;
}//double Nuclide::secularEquilibriumHalfLife() const


bool Nuclide::canObtainSecularEquilibrium() const
{
  return (secularEquilibriumHalfLife() < halfLife);
}//bool canObtainSecularEquilibrium() const


double Nuclide::promptEquilibriumHalfLife() const
{
  //20111012 Verified gives correct correct for U235, U238, Pu240, Pu242, Th232

  //XXX
  //  shouldnt this function just return the max 1/2 life of the daughter elements?
  //  Right now I'll leave below implementation, since this implinetation roughly
  //  follows some notes I found from by a fellow sandian, but I dont necassarily
  //  trust.
  if( IsInf(halfLife) )
    return 0.0;

  double promtHalflife = 0.0;

  for( size_t i = 0; i < decaysToChildren.size(); ++i )
  {
    const Nuclide *child = decaysToChildren[i]->child;
    if( child && (child != this) )
    {
      if( !IsInf(child->halfLife)
          && child->halfLife < this->halfLife )
        promtHalflife = max( promtHalflife,
                         max( child->halfLife,
                              child->promptEquilibriumHalfLife() ) );
    }
  }//for( size_t i = 0; i < decaysToChildren.size(); ++i )

  return promtHalflife;
}//double Nuclide::promptEquilibriumHalfLife() const


bool Nuclide::canObtainPromptEquilibrium() const
{
  return (promptEquilibriumHalfLife() > DBL_MIN);
}//bool Nuclide::canObtainPromptEquilibrium() const


double Nuclide::atomsPerGram() const
{
  //1 mole of C12 weighs 12 grams
  const double mole = 6.02214179E+23;
  return mole/atomicMass;
}//double Nuclide::atomsPerGram() const


double Nuclide::activityPerGram() const
{
  const double num_atoms_in_gram = atomsPerGram();
  return numAtomsToActivity( num_atoms_in_gram );
}//double Nuclide::activityPerGram() const


bool Nuclide::decaysToStableChildren() const
{
  for( size_t i = 0; i < decaysToChildren.size(); ++i )
  {
    const Transition *trans = decaysToChildren[i];
    if( trans->child && !IsInf(trans->child->halfLife) )
      return false;
  }//for( size_t i = 0; i < decaysToChildren.size(); ++i )

  return true;
}//bool decaysToStableChildren() const


Nuclide::Nuclide( const ::rapidxml::xml_node<char> *node )
{
  set( node );
}

//Nuclide::Nuclide()
//{
//};

Nuclide::~Nuclide()
{
}//~Nuclide()



void Nuclide::set( const ::rapidxml::xml_node<char> *node )
{
  //Note: assumes ::rapidxml::parse_no_string_terminators has NOT been used
  using ::rapidxml::internal::compare;
  typedef ::rapidxml::xml_attribute<char> Attribute;

  if( !compare( node->name(), node->name_size(), "nuclide", 7, true ) )
    throw runtime_error( "Nuclide::set(xml_node *node): node must be \"nuclide\"" );

  decaysToChildren.clear();

  const Attribute *symbol_att, *atomicNumber_att, *massNumber_att;
  const Attribute *isomerNumber_att, *atomicMass_att, *halfLife_att;

  symbol_att       = node->first_attribute( "symbol", 6 );
  atomicNumber_att = node->first_attribute( "atomicNumber", 12 );
  massNumber_att   = node->first_attribute( "massNumber", 10 );
  isomerNumber_att = node->first_attribute( "isomerNumber", 12 );
  atomicMass_att   = node->first_attribute( "atomicMass", 10 );
  halfLife_att     = node->first_attribute( "halfLife", 8 );

  if( symbol_att ) symbol = symbol_att->value();
  else             symbol = "";

  if( atomicNumber_att ) atomicNumber = (short int)atoi( atomicNumber_att->value() );
  else                   atomicNumber = 0;

  if( massNumber_att )   massNumber = (short int)atoi( massNumber_att->value() );
  else                   massNumber = 0;

  if( isomerNumber_att ) isomerNumber = (short int)atoi( isomerNumber_att->value() );
  else                   isomerNumber = 0;

  if( atomicMass_att ) atomicMass = (float)atof( atomicMass_att->value() );
  else                 atomicMass = 0.0;

  if( halfLife_att )
  {
    if( sscanf( halfLife_att->value(), "%lf", &halfLife ) != 1 )
    {
      if( !compare( halfLife_att->value(), halfLife_att->value_size(), "INF", 3, true )  )
        cerr << "Non-INF, non-numeric halfLife: " << halfLife_att->value() << endl;
      halfLife = numeric_limits<double>::infinity();
    }//if( ss.fail() )
  }else
  {
    halfLife = numeric_limits<double>::infinity();
  }
}//void set( ::rapidxml::xml_node<char> *node )




Element::Element( ::rapidxml::xml_node<char> *node,
                  const vector<const Nuclide *> &nuclides )
{
  set( node, nuclides );
}


void Element::set( ::rapidxml::xml_node<char> *node,
                   const vector<const Nuclide *> &nuclides )
{
  //Note: assumes ::rapidxml::parse_no_string_terminators has NOT been used
  using ::rapidxml::internal::compare;
  typedef ::rapidxml::xml_attribute<char> Attribute;

  if( !compare( node->name(), node->name_size(), "element", 7, true ) )
    throw runtime_error( "Element::set(xml_node *node): node must be \"element\"" );

  isotopes.clear();

  Attribute *symbol_att       = node->first_attribute( "symbol", 6 );
  Attribute *name_att         = node->first_attribute( "name", 4 );
  Attribute *atomicNumber_att = node->first_attribute( "atomicNumber", 12 );

  if( symbol_att ) symbol = symbol_att->value();
  else             symbol = "";

  if( name_att ) name = name_att->value();
  else           name = "";

  if( atomicNumber_att ) atomicNumber = (short int)atoi( atomicNumber_att->value() );
  else                   atomicNumber = 0;


  for( ::rapidxml::xml_node<char> *isotope = node->first_node( "isotope", 7 );
      isotope;
      isotope = isotope->next_sibling( "isotope", 7 ) )
  {
    Attribute *symbol_att = isotope->first_attribute( "symbol", 6 );
    Attribute *abundance_att = isotope->first_attribute( "abundance", 9 );

    if( symbol_att && abundance_att )
    {
      NuclideAbundancePair iso;
      const string sym = symbol_att->value();
      vector<const Nuclide *>::const_iterator pos = find( nuclides, sym );
      if( pos == nuclides.end() )
      {
        const string msg = "Could not find the symbol '"
                            + sym + "' in <element>-><isotope> node";
        throw runtime_error( msg );
      }//if( pos == nuclides.end() )

      iso.nuclide = (*pos);
      iso.abundance = atof( abundance_att->value() );

      isotopes.push_back( iso );
    }else throw runtime_error( "<element>\n        <isotope> node ill formed" );
  }//for( loop over isotopes )

  for( ::rapidxml::xml_node<char> *xray = node->first_node( "xray", 4 );
       xray;
       xray = xray->next_sibling( "xray", 4 ) )
  {
    Attribute *energy_att = xray->first_attribute( "energy", 6 );
    Attribute *intensity_att = xray->first_attribute( "relintensity", 12 );
    
    if( energy_att && intensity_att )
    {
      xrays.push_back( EnergyIntensityPair( atof(energy_att->value()), atof(intensity_att->value()) ) );
    }
  }
}//void Element::set(...)


double Element::atomicMass() const
{
  double atomicMass = 0.0;
  for( size_t i = 0; i < isotopes.size(); ++i )
  {
    const NuclideAbundancePair &p = isotopes[i];
    atomicMass += p.abundance * p.nuclide->atomicMass;
  }//for( loop over isotopes )


  //If the isotope is not natrally occuring, then atomicMass==0.0, so lets 
  //  supliment this.  Note that I got these atomic masses from CERNs ROOT
  //  using the code:
  // TGeoElementTable *t = new TGeoElementTable();
  // t->BuildDefaultElements()
  // for( int i = 1; i < 112; ++i){ cout << t->GetElement(i)->A() << ", ";}
  
  static const double atomic_masses[] = 
  { 
    0.0,  1.00794, 4.0026, 6.941, 9.01218, 10.811, 12.0107, 14.0067, 15.9994, 
    18.9984, 20.1797, 22.9898, 24.305, 26.9815, 28.0855, 30.9738, 32.066, 
    35.4527, 39.948, 39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938, 
    55.845, 58.9332, 58.6934, 63.546, 65.39, 69.723, 72.61, 74.9216, 78.96, 
    79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224, 92.9064, 95.94, 98, 
    101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76, 
    127.6, 126.904, 131.29, 132.905, 137.327, 138.905, 140.116, 140.908, 144.24, 
    145, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.26, 168.934, 
    173.04, 174.967, 178.49, 180.948, 183.84, 186.207, 190.23, 192.217, 
    195.078, 196.967, 200.59, 204.383, 207.2, 208.98, 209, 210, 222, 223, 226, 
    227, 232.038, 231.036, 238.029, 237, 244, 243, 247, 247, 251, 252, 257, 
    258, 259, 262, 261, 262, 263, 262, 265, 266, 269, 272 
  };//static const float atomic_masses[] = {...}

  if( atomicMass<1.0 )
    atomicMass = static_cast<float>( atomic_masses[atomicNumber] );

/*
  if( atomicMass<1.0 )
  {
  //XXX - isotopes is actually empty
    double maxHalfLife = 0.0;
    for( size_t i = 0; i < isotopes.size(); ++i )
    {
      const NuclideAbundancePair &p = isotopes[i];
      if( p.nuclide->halfLife > maxHalfLife )
      {
        maxHalfLife = p.nuclide->halfLife;
        atomicMass = p.nuclide->atomicMass;
      }
    }//for( loop over isotopes )
  }//if( atomicMass<1.0 )
*/
  
  return atomicMass;
}//double atomicMass() const


bool Element::lessThan( const Element *lhs, const Element *rhs )
{
  if( !lhs || !rhs )
    return false;

  return (lhs->atomicNumber < rhs->atomicNumber);
}//bool lessThan( const Element *lhs, const Element *rhs )


bool Element::lessThanForAtomicNum( const Element *lhs, const int atomicNum )
{
  if( !lhs )
    return false;

  return (lhs->atomicNumber < atomicNum);
}//bool lessThanForAtomicNum( const Element *lhs, const int atomicNum )


NuclideMixture::NuclideMixture()
{
}

NuclideMixture::~NuclideMixture()
{
}


void NuclideMixture::clear()
{
  m_originalNuclides.clear();
  m_decayedToNuclides.clear();
}//void NuclideMixture::clear()

void NuclideMixture::addNuclideByAbundance( const Nuclide *nuclide, double num_initial_atom )
{
  addNuclide( NuclideNumAtomsPair( nuclide, num_initial_atom )  );
}


void NuclideMixture::addNuclideByActivity( const Nuclide *nuclide, double activity )
{
  const double num_initial_atom = nuclide->activityToNumAtoms( activity );
  addNuclide( NuclideNumAtomsPair( nuclide, num_initial_atom ) );
}


void NuclideMixture::addAgedNuclideByActivity( const Nuclide *nuclide,
                                               double activity, double age )
{
  NuclideMixture ager;
  ager.addNuclideByActivity( nuclide, activity );
  const std::vector<NuclideActivityPair> aged_activities = ager.activity( age );

  double age_sf = 1.0;
  for( size_t i = 0; i < aged_activities.size(); ++i )
  {
    const NuclideActivityPair &red = aged_activities[i];
    if( red.nuclide == nuclide )
    {
      age_sf = activity / red.activity;
      i = aged_activities.size();
    }
  }//for( size_t i = 0; i < aged_activities.size(); ++i )

  for( size_t i = 0; i < aged_activities.size(); ++i )
  {
    const NuclideActivityPair &red = aged_activities[i];
    addNuclideByActivity( red.nuclide, age_sf * red.activity );
  }
}//void addAgedNuclideByActivity(...)


void NuclideMixture::addNuclide( const NuclideNumAtomsPair &nuclide )
{
  m_decayedToNuclides.clear();
  m_originalNuclides.push_back( nuclide );
}


void NuclideMixture::addNuclide( const NuclideActivityPair &nucidePair )
{
  const Nuclide *nuclide = nucidePair.nuclide;
  const double num_atoms = nuclide->activityToNumAtoms( nucidePair.activity );
  addNuclide( NuclideNumAtomsPair( nuclide, num_atoms ) );
}



bool NuclideMixture::addNuclideInSecularEquilibrium( const Nuclide *parent,
                                                     double parent_activity )
{
  const double secEqHL = parent->secularEquilibriumHalfLife();
  if( (secEqHL >= parent->halfLife) || (secEqHL<=0.0)  )
  {
    cerr << "NuclideMixture::addNuclideInSecularEquilibrium( " << parent->symbol
         << ", " << parent_activity/becquerel << "bq ): Nuclide "
         << parent->symbol << " can not reach secular equilibrium! I am "
         << "cowardly refusing to add it to the NuclideMixture." << endl;
    return false;
  }//if( (secEqHL >= parent->halfLife) || (secEqHL<=0.0)  )

  const vector<const Transition *> &transitions = parent->decaysToChildren;

  addNuclideByActivity( parent, parent_activity );

  for( size_t trans = 0; trans < transitions.size(); ++trans )
  {
    const Nuclide *child = transitions[trans]->child;
    const float branchRatio = transitions[trans]->branchRatio;

    if( child && (child!=parent) && !IsInf(child->halfLife))
      addNuclideInSecularEquilibrium( child, branchRatio*parent_activity );
  }//for( size_t trans = 0; trans < transitions.size(); ++trans )
  
  return true;
}//void NuclideMixture::addNuclideInSecularEquilibrium(...)


void NuclideMixture::addNuclideInPromptEquilibrium( const Nuclide *parent,
                                                   double activity )
{
  if( IsInf(parent->halfLife) )
    return;

  addNuclideByActivity( parent, activity );

  for( size_t i = 0; i < parent->decaysToChildren.size(); ++i )
  {
    const Transition *trans = parent->decaysToChildren[i];
    const Nuclide *child = trans->child;
    if( child && (child->halfLife < parent->halfLife) )
      addNuclideInPromptEquilibrium( child, trans->branchRatio*activity );
  }//for( size_t i = 0; i < decaysToChildren.size(); ++i )
}//void addNuclideInPromptEquilibrium(...)



int NuclideMixture::numInitialNuclides() const
{
  return static_cast<int>( m_originalNuclides.size() );
}

int NuclideMixture::numSolutionNuclides() const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();
  return static_cast<int>( m_decayedToNuclides.size() );
}


double NuclideMixture::activity( double time_in_seconds, const Nuclide *nuclide )    const
{
  return activity( time_in_seconds, internalIndexNumber( nuclide ) );
}

double NuclideMixture::activity( double time_in_seconds, const std::string &symbol ) const
{
  return activity( time_in_seconds, internalIndexNumber( symbol ) );
}

double NuclideMixture::activity( double time_in_seconds, int z, int atomic_mass, int iso )    const
{
  return activity( time_in_seconds, internalIndexNumber( z, atomic_mass, iso ) );
}

double NuclideMixture::activity( double time_in_seconds, int index ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  const int n_nuclides = static_cast<int>( m_decayedToNuclides.size() );
  if( (index < 0) || (index >= n_nuclides) )
    throw std::runtime_error( "NuclideMixture::activity(...): hecka ish" );

  return m_decayedToNuclides[index].activity( time_in_seconds );
}//double activity( double time_in_seconds, int internal_index_number ) const


double NuclideMixture::totalActivity( double time ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  double activity = 0.0;
  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
    activity += m_decayedToNuclides[i].activity( time );

  return activity;
}//double NuclideMixture::totalActivity( double time ) const


double NuclideMixture::totalMassInGrams( double time ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  double mass_amu = 0.0;
  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
  {
    const NuclideTimeEvolution &evo = m_decayedToNuclides[i];
    const Nuclide *nuc = evo.nuclide;
    mass_amu += evo.numAtoms( time ) * nuc->atomicMass;
  }//for( loop over m_decayedToNuclides )

  const double mole = 6.02214179E+23;  //AMU/gram
  return mass_amu/mole;
}//double totalMassInGrams( double time )


double NuclideMixture::numAtoms( double time_in_seconds, const Nuclide *nuclide )    const
{
  return numAtoms( time_in_seconds, internalIndexNumber( nuclide ) );
}


double NuclideMixture::numAtoms( double time_in_seconds, const std::string &symbol ) const
{
  return numAtoms( time_in_seconds, internalIndexNumber( symbol ) );
}


double NuclideMixture::numAtoms( double time_in_seconds, int z, int atomic_mass, int iso )    const
{
  return numAtoms( time_in_seconds, internalIndexNumber( z, atomic_mass, iso ) );
}//double numAtoms( double time_in_seconds, int z, int atomic_mass )  const;


double NuclideMixture::numAtoms( double time_in_seconds, int index ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  const int n_nuclides = static_cast<int>( m_decayedToNuclides.size() );
  if( (index < 0) || (index >= n_nuclides) )
    throw std::runtime_error( "NuclideMixture::numAtoms(...): Invalid Index" );

  return m_decayedToNuclides[index].numAtoms( time_in_seconds );
}//double numAtoms( double time_in_seconds, int internal_index_number ) const;


vector<NuclideActivityPair> NuclideMixture::activity( double time_in_seconds )  const
{
  vector<NuclideActivityPair> answer;

  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
  {
    const Nuclide *nuclide = m_decayedToNuclides[i].nuclide;
    const double activity = m_decayedToNuclides[i].activity( time_in_seconds );
    answer.push_back( NuclideActivityPair( nuclide, activity )  );
  }//for( loop over atoms decayed too, i )

  return answer;
}//vector<NuclideActivityPair> activity( double time_in_seconds )  const


vector<NuclideNumAtomsPair> NuclideMixture::numAtoms( double time_in_seconds )  const
{
  vector<NuclideNumAtomsPair> answer;
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
  {
    const Nuclide *nuclide = m_decayedToNuclides[i].nuclide;
    const double numAtoms = m_decayedToNuclides[i].numAtoms( time_in_seconds );
    answer.push_back( NuclideNumAtomsPair( nuclide, numAtoms )  );
  }//for( loop over atoms decayed too, i )

  return answer;
}//vector<NuclideNumAtomsPair> numAtoms( double time_in_seconds )  const


bool includesNuclide( const std::vector<NuclideNumAtomsPair> &lhs, const Nuclide *nuc )
{
  for( size_t i = 0; i < lhs.size(); ++i )
  {
    if( lhs[i].nuclide == nuc )
      return true;
  }
  return false;
}


vector<EnergyRatePair> NuclideMixture::decayParticle( double time,
                                                           ProductType type,
                                                           HowToOrder sortType ) const
{
  //TODO: from profiling it turns out this function is incredibly slow!
  //      It appears most of this innefieceincy is associated with using
  //      the map<double,double>

  //Should make a vector of energies first, sort them, then a second vector
  //  of abundances made to coordespond to the vector of energies

  //(note that all times given for debug build, release build is ~10x faster)
  //For original version of code, on my MacBookPro, for the Ba133 spectrum to associate IDs
  //  For 5880 gamma calls:  5.221995s wall, 5.260000s user + 0.000000s system = 5.260000s CPU (100.7%)
  //
  //For version of code written 20120928 (eg one used when USE_NEW_VERSION is true)
  //  for 5880 gamma calls:  4.237621s wall, 4.210000s user + 0.000000s system = 4.210000s CPU (99.3%)

  //Pre 20121010 I had thought the call to activity was harmless (1.8% of time),
  //  but using "Instruments" in OSX, instead of "Activity Monitor" profiling
  //  shows it taking up close to 100% of the time of this function, at least
  //  for the case fitting shielding and source activities in the spectrum
  //  viewer - the inneficient function looks to be SandiaDecayDataBase::getTimeEvolution(...)
  const vector<NuclideActivityPair> activities = activity( time ); //takes up 1.8% of the time for this function call

#define USE_NEW_VERSION 1


#if( USE_NEW_VERSION )
  vector<double> energies, branchratios;
//  energies.reserve( 2048 );
//  branchratios.reserve( 2048 );

//If we only use RadParticles which have
//  activity * particle.intensity>0.000001*total_activity
//  we bring time down from 4.21 seconds to 3.65 second
//  double total_activity = 0.0;
//  for( size_t i = 0; i < activities.size(); ++i )
//    total_activity += activities[i].activity;


  for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
  {
    const Nuclide *nuclide = activities[nucIndex].nuclide;
    const double activity = activities[nucIndex].activity;

    const size_t n_decaysToChildren = nuclide->decaysToChildren.size();

    for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    {
      const Transition *transition = nuclide->decaysToChildren[decayIndex];
      const size_t n_products = transition->products.size();

      for( size_t productNum = 0; productNum < n_products; ++productNum )
      {
        const RadParticle &particle = transition->products[productNum];
        if( particle.type == type )
        {

  //Comment from Greg Thoreson in his code to read tori.db3:
  // The ToRI x-ray database will return x-rays for GS->isomeric state transitions,
  // even if the x-rays originate from the isomeric de-excitation
  // so, we don't want to repeat different isomeric states of the same daughter
  //
  // Gregs protection against this issue didnt appear to be consistent with the
  //  comment, so it is _very_ likely wcjohns is mis-interpreting things, and
  //  the protection he made bellow is incorrect!
  //
  //20160114: wcjohns implemented not adding in contribution of of xrays from
  //  isomeric decays to a lower iso state, if the iso state was not the
  //  originating nuclide
  //
  //Note: wcjohns does not understand atoms with A>1.  For instance, he doesnt
  //  get why isomeric transforms would produce x-rays, as the change of
  //  potential on the electrons caused by the reoganization of the nucleus
  //  seems like it should at best result in very low energy x-rays, but this
  //  is not the case, as shown by:
  //    http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=850504
          if( type == XrayParticle )
          {
            //XXX: Note bellow neglects situations where you add a nuclide to
            //  the NuclideMixture that is a iso>0 state, and its a decendant of
            //  another nuclide that you added
            if( transition->parent->isomerNumber
//               && !transition->child->isomerNumber
               && transition->parent && transition->child
               && transition->parent->atomicNumber == transition->child->atomicNumber
               && transition->parent->massNumber == transition->child->massNumber
               && !includesNuclide( m_originalNuclides, transition->parent )
               )
            {
              continue;
            }
          }//if( type == XrayParticle )
          
          energies.push_back( particle.energy );
          branchratios.push_back( activity * particle.intensity * transition->branchRatio );
        }//if( particle.type == type )
      }//for( size_t productNum = 0; productNum < n_products; ++productNum )
    }//for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
  }//for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )

  const size_t nenergies = energies.size();
  vector<size_t> sort_indices( nenergies );
  for( size_t i = 0; i < nenergies; ++i )
    sort_indices[i] = i;

  std::sort( sort_indices.begin(), sort_indices.end(),
             index_compare_assend<vector<double>&>(energies) );

  double current_energy = -1.0;
  vector<EnergyRatePair> answer;
  answer.reserve( energies.size() );

  for( size_t i = 0; i < nenergies; ++i )
  {
    const double energy = energies[sort_indices[i]];
    const double branchratio = branchratios[sort_indices[i]];

    if( energy != current_energy )
    {
      current_energy = energy;
      answer.push_back( EnergyRatePair(branchratio,energy) );
    }else
    {
      answer.back().numPerSecond += branchratio;
    }
  }//for( size_t i = 0; i < energies.size(); ++i )

#else

  map<double,double> energyToRateMap;
  typedef map<double,double>::value_type MapDD_t;

  for( size_t nucIndex = 0; nucIndex < activities.size(); ++nucIndex )
  {
    const Nuclide *nuclide = activities[nucIndex].nuclide;
    const double activity = activities[nucIndex].activity;

    const size_t n_decaysToChildren = nuclide->decaysToChildren.size();

    for( size_t decayIndex = 0; decayIndex < n_decaysToChildren; ++decayIndex )
    {
      const Transition *transition = nuclide->decaysToChildren[decayIndex];
      const size_t n_products = transition->products.size();

      for( size_t productNum = 0; productNum < n_products; ++productNum )
      {
        const RadParticle &particle = transition->products[productNum];
        if( particle.type == type )
        {
          map<double,double>::iterator pos = energyToRateMap.find( particle.energy ); //takes up about 11% of the time for this funciton call

          if( pos == energyToRateMap.end() )
            pos = energyToRateMap.insert( MapDD_t( particle.energy, 0.0 ) ).first;  //takes up about 35% of the time for this function call

          pos->second += activity * particle.intensity * transition->branchRatio;
        }//if( this is a gamma )
      }//for( loop over decay products, productNum )
    }//for( loop over types of decays, decayIndex )
  }//for( loop over nuclides decayed to, nucIndex )


  vector<AbundanceEnergyPair> answer;
  answer.reserve( energyToRateMap.size() );

  map<double,double>::iterator iter;
  for( iter = energyToRateMap.begin(); iter != energyToRateMap.end(); ++iter )
    answer.push_back( AbundanceEnergyPair( iter->second, iter->first ) );     //~1% of the time for this function call

/*
  cerr << "answer.size()=" << answer.size() << ", answer2.size()=" << answer2.size() << endl;
  for( size_t i = 0; i < min(answer.size(),answer2.size()); ++i )
    cerr << "\t" << i << ": (" << answer[i].energy << ", "  << answer[i].numPerSecond << ") - (" << answer2[i].energy << ", "  << answer2[i].numPerSecond << ")" << endl;
  assert( answer.size() == answer2.size() );
*/
#endif //#if( USE_NEW_VERSION )

  switch( sortType )
  {
    case OrderByAbundance:
      sort( answer.begin(), answer.end(),
           &EnergyRatePair::moreThanByNumPerSecond );
      break;

    case OrderByEnergy:
      //already sorted by energy
      break;
  }//switch( sortType )

  return answer;
}//vector<AbundanceEnergyPair> decayParticle( double time_in_seconds, ProductType type, SortType sortType = kByAbundance ) const



vector<EnergyRatePair> NuclideMixture::gammas( const double time_in_seconds,
                                                    const NuclideMixture::HowToOrder sortType,
                                                    const bool includeAnnihilation ) const
{
  if( !includeAnnihilation )
    return decayParticle( time_in_seconds, GammaParticle, sortType );

  vector<EnergyRatePair> gams = decayParticle( time_in_seconds, GammaParticle, sortType );
  const vector<EnergyRatePair> poss = decayParticle( time_in_seconds, PositronParticle, OrderByEnergy );

  EnergyRatePair possbr( 0.0, 510.998910 );
  
  for( size_t i = 0; i < poss.size(); ++i )
    possbr.numPerSecond += 2.0*poss[i].numPerSecond; //an annihaltion event creates two photons
  
  vector<EnergyRatePair>::iterator pos;
  switch( sortType )
  {
    case OrderByAbundance:
      pos = lower_bound( gams.begin(), gams.end(), possbr, &EnergyRatePair::moreThanByNumPerSecond );
    break;
        
    case OrderByEnergy:
      pos = lower_bound( gams.begin(), gams.end(), possbr, &EnergyRatePair::lessThanByEnergy );
    break;
  }//switch( sortType )
    
  gams.insert( pos, possbr );
  
  return gams;
}//vector<AbundanceEnergyPair> gammas( double time_in_seconds, SortType sortType ) const


std::vector<EnergyRatePair> NuclideMixture::alphas( double time_in_seconds, NuclideMixture::HowToOrder sortType ) const
{
  return decayParticle( time_in_seconds, AlphaParticle, sortType );
}//vector<AbundanceEnergyPair> gammas( double time_in_seconds, SortType sortType ) const


std::vector<EnergyRatePair> NuclideMixture::betas( double time_in_seconds,  NuclideMixture::HowToOrder sortType ) const
{
  return decayParticle( time_in_seconds, BetaParticle, sortType );
}


std::vector<EnergyRatePair> NuclideMixture::betaPlusses( double time_in_seconds, NuclideMixture::HowToOrder sortType ) const
{
  return decayParticle( time_in_seconds, PositronParticle, sortType );
}


std::vector<EnergyRatePair> NuclideMixture::xrays( double time_in_seconds, NuclideMixture::HowToOrder sortType ) const
{
  return decayParticle( time_in_seconds, XrayParticle, sortType );
}

std::vector<EnergyRatePair> NuclideMixture::photons( const double time_in_seconds, const NuclideMixture::HowToOrder sortType ) const
{
  vector<EnergyRatePair> gam = gammas( time_in_seconds, sortType, true );
  const vector<EnergyRatePair> x = xrays( time_in_seconds, sortType );
  
  gam.insert( gam.begin(), x.begin(), x.end() );
  
  switch( sortType )
  {
    case OrderByAbundance:
      sort( gam.begin(), gam.end(), &EnergyRatePair::moreThanByNumPerSecond );
      break;
      
    case OrderByEnergy:
      sort( gam.begin(), gam.end(), &EnergyRatePair::lessThanByEnergy );
      break;
  }//switch( sortType )
  
  return gam;
}


//internal_index_number runs from 0 to numSolutionNuclides()
int NuclideMixture::internalIndexNumber( const Nuclide *nuclide ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
    if( m_decayedToNuclides[i].nuclide->symbol == nuclide->symbol ) //we could just compare pointer values...
      return static_cast<int>( i );

  throw std::runtime_error( nuclide->symbol + " is not in this NuclideMixture" );
  return -1;
}//int internalIndexNumber( const Nuclide *nuclide )


int NuclideMixture::internalIndexNumber( const std::string &symbol ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
    if( m_decayedToNuclides[i].nuclide->symbol == symbol )
      return static_cast<int>( i );

  const string msg = symbol + " is not in this NuclideMixture.  Did you perhaps"
                              " not put symbol into its normailzed form (ex. "
                              "Co60, U238m, etc)?";
  throw std::runtime_error( msg );
  return -1;
}//int internalIndexNumber( const std::string &symbol ) const


int NuclideMixture::internalIndexNumber(int z, int atomic_mass, int iso ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  for( size_t i = 0; i < m_decayedToNuclides.size(); ++i )
  {
    const Nuclide *nuclide = m_decayedToNuclides[i].nuclide;
    if( nuclide->atomicMass == atomic_mass
        && nuclide->atomicNumber == z
        && nuclide->isomerNumber == iso )
      return static_cast<int>( i );
  }

  char buffer[256];
  snprintf( buffer, sizeof(buffer),
            "Z=%i AN=%i, ISO=%i is not in this NuclideMixture.",
            z, atomic_mass, iso );
  throw std::runtime_error( buffer );
  return -1;
}//int internalIndexNumber(int z, int atomic_mass ) const


const Nuclide *NuclideMixture::initialNuclide( int index ) const
{
  return m_originalNuclides.at( index ).nuclide;
}//const Nuclide *initialNuclide( int index ) const


const Nuclide *NuclideMixture::solutionNuclide( int internal_index_number ) const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  return m_decayedToNuclides.at( internal_index_number ).nuclide;
}//const Nuclide *solutionNuclide( int internal_index_number ) const



std::string NuclideMixture::info( const double time ) const
{
  if( m_originalNuclides.empty() )
    return "";

  const vector<NuclideActivityPair> activities = activity( time );
  const vector<EnergyRatePair> gammasAbundances = gammas( time, NuclideMixture::OrderByAbundance, true );
  const vector<EnergyRatePair> alphasAbundances = alphas( time );
  const vector<EnergyRatePair> betasAbundances = betas( time );
  const vector<EnergyRatePair> betaPlussesAbundances = betaPlusses( time );

  stringstream infostrm;

  infostrm << "Starting from:";

  {//begin codeblock to put starting nuclides information into the stream
    stringstream tempstrm;
    tempstrm << "\n    ";

    for( size_t i = 0; i < m_originalNuclides.size(); ++i )
    {
      const NuclideNumAtomsPair &pair = m_originalNuclides[i];
      if( i ) tempstrm << ",";
      if( tempstrm.str().length() > 78 )
      {
        infostrm << tempstrm.str();
        tempstrm.clear();
        tempstrm << "\n    ";
      }//if( we've reached our length limit )

      tempstrm << pair.nuclide->numAtomsToActivity( pair.numAtoms ) / becquerel
               << " Bq " << pair.nuclide->symbol;
    }//for( loop over orignal nuclides )

    infostrm << tempstrm.str();
  }//end codeblock to put starting nuclides information into the stream


  infostrm << "\nThe following nuclides are present at " << time/second
           << " seconds:";


  {//begin codeblock to put decayed to nuclides information into the stream
    stringstream tempstrm;
    tempstrm << "\n    ";

    for( size_t i = 0; i < activities.size(); ++i )
    {
      const NuclideActivityPair &pair = activities[i];
      if( i ) tempstrm << ",";
      if( tempstrm.str().length() > 78 )
      {
        infostrm << tempstrm.str();
        tempstrm.clear();
        tempstrm << "\n    ";
      }//if( we've reached our length limit )

      tempstrm << " " << pair.activity / becquerel << " Bq " << pair.nuclide->symbol;
    }//for( loop over activities )

    infostrm << tempstrm.str();
  }//end codeblock to put decayed to nuclides information into the stream

  if( gammasAbundances.size() )
  {
    infostrm << "\nGammas Present:";
    for( size_t i = 0; i < gammasAbundances.size(); ++i )
      infostrm << "\n    " << gammasAbundances[i].energy/keV << " keV\t-\t"
               << gammasAbundances[i].numPerSecond << " per second";
  }//if( gammasAbundances.size() )

  if( alphasAbundances.size() )
  {
    infostrm << "\nAlphas Present:";
    for( size_t i = 0; i < alphasAbundances.size(); ++i )
      infostrm << "\n    " << alphasAbundances[i].energy/keV << " keV\t-\t"
              << alphasAbundances[i].numPerSecond << " per second";
  }//if( alphasAbundances.size() )

  if( betasAbundances.size() )
  {
    infostrm << "\nBetas Present:";
    for( size_t i = 0; i < betasAbundances.size(); ++i )
      infostrm << "\n    " << betasAbundances[i].energy/keV << " keV\t-\t"
               << betasAbundances[i].numPerSecond << " per second";
  }//if( alphasAbundances.size() )


  if( betaPlussesAbundances.size() )
  {
    infostrm << "\nBeta+'s Present:";
    for( size_t i = 0; i < betaPlussesAbundances.size(); ++i )
      infostrm << "\n    " << betaPlussesAbundances[i].energy/keV << " keV\t-\t"
               << betaPlussesAbundances[i].numPerSecond << " per second";
  }//if( betaPlussesAbundances.size() )


  return infostrm.str();
}//string info() const

const std::vector<NuclideTimeEvolution> &NuclideMixture::decayedToNuclidesEvolutions() const
{
  if( m_decayedToNuclides.empty() )
    performTimeEvolution();

  return m_decayedToNuclides;
}//const std::vector<NuclideTimeEvolution> &NuclideMixture::decayedToNuclidesEvolutions() const



void NuclideMixture::performTimeEvolution() const
{
  m_decayedToNuclides = SandiaDecayDataBase::getTimeEvolution( m_originalNuclides );
}//void performTimeEvolution()


SandiaDecayDataBase::SandiaDecayDataBase( const std::string &filename )
  : m_xmlFileContainedDecayXrays( false ),
    m_xmlFileContainedElementalXrays( false )
{
  parseSandiaDecayXml( filename, m_nuclides, m_nuclideStore,
                        m_elements, m_elementStore, m_transitionStore );
  
  checkIfContainedXrays();
}

  
SandiaDecayDataBase::SandiaDecayDataBase()
  : m_xmlFileContainedDecayXrays( false ),
    m_xmlFileContainedElementalXrays( false )
{
}

  
void SandiaDecayDataBase::checkIfContainedXrays()
{
  m_xmlFileContainedDecayXrays = false;
  for( size_t i = 0; !m_xmlFileContainedDecayXrays && (i < m_transitionStore.size()); ++i )
    for( size_t j = 0; !m_xmlFileContainedDecayXrays && (j < m_transitionStore[i].products.size()); ++j )
      m_xmlFileContainedDecayXrays = (m_transitionStore[i].products[j].type == XrayParticle);
  
  m_xmlFileContainedElementalXrays = ((m_elements.size()>=93) && !(m_elements[92]->xrays.empty()));
}//void checkIfContainedXrays()
  
  
bool SandiaDecayDataBase::initialized() const
{
  return !m_nuclides.empty();
}//bool initialized() const


void SandiaDecayDataBase::initialize( const std::string &filename )
{
  if( !m_nuclides.empty() )
    throw std::runtime_error( "SandiaDecayDataBase::initialize(...): Database "
                              "already initialezed form XML file" );
  parseSandiaDecayXml( filename, m_nuclides, m_nuclideStore,
                       m_elements, m_elementStore, m_transitionStore );
  
  checkIfContainedXrays();
}//void initialize( const std::string &filename )

  
void SandiaDecayDataBase::reset()
{
  m_nuclides.clear();
  m_nuclideStore.clear();
  
  m_elements.clear();
  m_elementStore.clear();

  m_transitionStore.clear();
  
  //Reclaim memory
  {
    vector<Element> tempElementStore;
    tempElementStore.swap( m_elementStore );
  }
  
  {
    vector<Nuclide> tempNuclideStore;
    tempNuclideStore.swap( m_nuclideStore );
  }
  
  {
    vector<Transition> tempTransitionStore;
    tempTransitionStore.swap( m_transitionStore );
  }
  
//#if __cplusplus <= 19971  1L
//#else
//  m_elementStore.shrink_to_fit();
//  m_nuclideStore.shrink_to_fit();
//  m_transitionStore.shrink_to_fit();
//#endif

  m_xmlFileContainedDecayXrays = false;
  m_xmlFileContainedElementalXrays = false;
}
  
  
SandiaDecayDataBase::~SandiaDecayDataBase()
{
  //unInitialize();
}//~SandiaDecayDataBase()



const vector<const Nuclide *> &SandiaDecayDataBase::nuclides() const
{
  return m_nuclides;
}

const std::vector<const Element *> &SandiaDecayDataBase::elements() const
{
  return m_elements;
}

const std::vector<Transition> &SandiaDecayDataBase::transitions() const
{
  return m_transitionStore;
}
  
const Element *SandiaDecayDataBase::element( const int atomicNumber ) const
{
  typedef vector<const Element *>::const_iterator ElementIter;

  const ElementIter begin = m_elements.begin();
  const ElementIter end = m_elements.end();

#ifdef _WIN32
  ElementIter pos = begin;
  while( pos != end && Element::lessThanForAtomicNum( *pos, atomicNumber ) )
    ++pos;
#else
  const ElementIter pos = lower_bound( begin, end, atomicNumber,
                                       &Element::lessThanForAtomicNum );
#endif

  if( pos == end )
    return NULL;

  return *pos;
}//const Element *SandiaDecayDataBase::element( const int massNumber ) const


const Element *SandiaDecayDataBase::element( const std::string &input ) const
{

  string el;
  for( size_t i = 0; i < input.size(); ++i )
    if( isalpha( input[i] ) )
      el += tolower( input[i] );

  if( el.empty() )
    return NULL;

  const bool isSymbol = ( (el.length() <= 2)
                        || (el.length()==3 && el.substr(0,2)=="uu" ) );

  if( isSymbol )
   el[0] = toupper( el[0] );

  for( size_t i = 0; i < m_elements.size(); ++i )
  {
    const Element *testEl = m_elements[i];
    if( testEl->symbol == el || testEl->name == el )
      return testEl;
  }//for( loop over m_elements to match the current string to, i )

  return NULL;
}//const Element *element( std::string label ) const



void SandiaDecayDataBase::parseSandiaDecayXml( const std::string &filename,
                                               vector<const Nuclide *> &nuclides,
                                               std::vector<Nuclide> &stored_nuclides,
                                               std::vector<const Element *> &elements,
                                               std::vector<Element> &stored_elements,
                                               std::vector<Transition> &transitionStore )
{
  vector<char> data;

  {//Begin codeblock for reading the data from disk
    //this method for reading the file was taken from rapidxml_utils.hpp
    //  and is much faster than other methods tried
    basic_ifstream<char> stream(filename.c_str(), ios::binary);
    if (!stream)
      throw runtime_error(string("Cannot open file ") + filename);
    stream.unsetf(ios::skipws);

    // Determine stream size
    stream.seekg(0, ios::end);
    size_t size = stream.tellg();
    stream.seekg(0);

    // Load data and add terminating 0
    data.resize(size + 1);
    stream.read(&data.front(), static_cast<streamsize>(size));
    data[size] = 0;
  }//end codeblock for reading the data from disk

  elements.clear();
  stored_elements.clear();

  try
  {
    using namespace ::rapidxml;

    xml_document<char> doc;
    doc.parse<parse_normalize_whitespace | parse_trim_whitespace>( &data.front() );

    const xml_node<char> *document_node = doc.first_node();

    nuclides.clear();
    stored_nuclides.clear();
    
    stored_nuclides.reserve( 3720 );
    
    for( const xml_node<char> *node = document_node->first_node( "nuclide", 7 );
        node;
        node = node->next_sibling( "nuclide", 7 ) )
    {
#if __cplusplus > 199711L
      stored_nuclides.emplace_back( Nuclide(node) );
#else
      stored_nuclides.push_back( Nuclide(node) );
#endif
    }//for( loop over "nuclide" nodes )

    const size_t nnuc = stored_nuclides.size();
    nuclides.resize( nnuc );
    for( size_t i = 0; i < nnuc; ++i )
      nuclides[i] = &stored_nuclides[i];
    std::sort( nuclides.begin(), nuclides.end(), &pointerLessThanBySymbol );
    
    transitionStore.clear();
    transitionStore.reserve(3400); //actually 3376
    
    for( const xml_node<char> *node = document_node->first_node( "transition", 10 );
        node;
        node = node->next_sibling( "transition", 10 ) )
    {
      transitionStore.push_back( Transition(node, nuclides) );
      Transition *trans = &transitionStore.back();

      //XXX 20111012: the sandia.decay.xml file list transitions for isotopes with
      //    halfLife==inf, which can screw up some of the decay logic, since if
      //    the half lilfe is inf, there should be no decays, right?
      const bool use = ( (trans->parent != NULL)
                         && (trans->branchRatio > 0.0)
                         && !IsInf(trans->parent->halfLife) );

      //Perhaps the 'Transition' should insert itself into it's parent 'Nuclide' ?
      if( use )
      {
        Nuclide *child = const_cast<Nuclide *>( trans->child );
        Nuclide *parent = const_cast<Nuclide *>( trans->parent );
        parent->decaysToChildren.push_back( trans );
        if( child )
          child->decaysFromParents.push_back( trans );
      }else
        transitionStore.pop_back();
    }//for( loop over "transition" nodes )

    
    vector<Element *> els; //we need non-const pointers to the Elements to add
                           //  xrays to them

    stored_elements.reserve( 119 );
    for( xml_node<char> *node = document_node->first_node( "element", 7 );
       node;
       node = node->next_sibling( "element", 7 ) )
    {
      stored_elements.push_back( Element( node, nuclides ) );
    }
    
    const size_t nelements = stored_elements.size();
    els.resize( nelements );
    elements.reserve( nelements );
    
    for( size_t i = 0; i < nelements; ++i )
      els[i] = &stored_elements[i];
    std::sort( els.begin(), els.end(), &Element::lessThan );
    elements.insert( elements.end(), els.begin(), els.end() );
    
    //Check to make sure the xray data was actually in the file.
    //  This could probably be skipped to save a little time.
    //if( nelements<93 || elements[92]->xrays.empty() )
    //  throw runtime_error( filename + " did not contain elemental xray data: (" + elements[92]->symbol + ")" );
  }catch (::rapidxml::parse_error &e)
  {
    long line = static_cast<long>( std::count(&data.front(), e.where<char>(), '\n') + 1 );
    cerr << "Error parsing file on line " << line << endl;
    stored_nuclides.clear();
    nuclides.clear();
  }catch( std::exception &e )
  {
    std::cerr << "Exception " << e.what()
              << " caught when reading in SandiaDecay XML" << std::endl;
    stored_nuclides.clear();
    nuclides.clear();
  }
}//void SandiaDecayDataBase::parseSandiaDecayXml(...)


bool SandiaDecayDataBase::xmlContainedDecayXRayInfo() const
{
  return m_xmlFileContainedDecayXrays;
}
  
bool SandiaDecayDataBase::xmlContainedElementalXRayInfo() const
{
  return m_xmlFileContainedElementalXrays;
}
  
const Nuclide *SandiaDecayDataBase::nuclide( int z, int massNum, int iso ) const
{
  //This function can be optimized in the future!

  vector<const Nuclide *>::const_iterator pos;
  for( pos = m_nuclides.begin(); pos != m_nuclides.end(); ++pos )
  {
    const Nuclide *nuc = (*pos);
    if( nuc->atomicNumber==z && nuc->massNumber==massNum && nuc->isomerNumber==iso )
      return nuc;
  }

  return NULL;
}//const Nuclide *nuclide( int z, int atomicNumber ) const


const Nuclide *SandiaDecayDataBase::nuclide( const std::string &label ) const
{
  try
  {
    const string symbol = toNormalSymbolForm( label );
    // This will do a binary search to locate the nuclide with the symbol, if
    // such a nuclide exists.
    vector<const Nuclide *>::const_iterator pos = find( m_nuclides, symbol );

    if( pos == m_nuclides.end() )
      return NULL;
    return (*pos);
  }catch(...)
  {
  }

  return NULL;
}//const Nuclide *nuclide( const std::string &symbol ) const


const Nuclide *SandiaDecayDataBase::nuclideUnchecked(
                                                    const string &symbol ) const
{
  vector<const Nuclide *>::const_iterator pos = find( m_nuclides, symbol );

  if( pos == m_nuclides.end() )
    return NULL;

  return (*pos);
}//nuclideUnchecked(..)



//Get the nuclides corresponding to a given element
std::vector<const Nuclide *> SandiaDecayDataBase::nuclides(
                                                const Element *element ) const
{
  //  XXX - untested!

  if( !element )
    return std::vector<const Nuclide *>();

  vector<const Nuclide *>::const_iterator start = lower_bound( m_nuclides.begin(),
                                                  m_nuclides.end(),
                                                  element->symbol,
                                                  &lessThanBySymbol );
  vector<const Nuclide *>::const_iterator end = start;
  while( end != m_nuclides.end() )
  {
    string nucSym = (*end)->symbol;
    string::size_type pos = nucSym.find_first_of( "0123456789" );
    nucSym = nucSym.substr( 0, pos );

    if( nucSym != element->symbol )
      break;
    ++end;
  }//while( end != m_nuclides.end() )

  return std::vector<const Nuclide *>( start, end );
}//nuclides( const Element * )


std::vector<const Nuclide *> SandiaDecayDataBase::nuclides(
                                             const std::string &input ) const
{
  //  XXX - untested!

  string el;
  for( size_t i = 0; i < input.size(); ++i )
    if( isalpha( input[i] ) )
      el += tolower( input[i] );

  if( el.empty() )
    throw runtime_error( "SandiaDecayDataBase::nuclides( const string &):"
                         " String passed in contained no alphabetical chars" );

  const bool isSymbol = ( (el.length() <= 2)
                        || (el.length()==3 && el.substr(0,2)=="uu" ) );

  if( isSymbol )
   el[0] = toupper( el[0] );

  const Element *element = NULL;
  for( size_t i = 0; !element && (i < m_elements.size()); ++i )
  {
    const Element *testEl = m_elements[i];
    if( testEl->symbol == el || testEl->name == el )
      element = testEl;
  }//for( loop over m_elements to match the current string to, i )

  if( !element )
  {
    stringstream msg;
    msg << "SandiaDecayDataBase::nuclides( \"" << input << "\" ):"
        << " A element with name or symbol '" << el << "' couldnt be found";
    throw runtime_error( msg.str() );
  }//if( !foundSymbol )

  return nuclides( element );
}//nuclides( const std::string &element ) const




std::string SandiaDecayDataBase::toNormalSymbolForm( const std::string &symbol ) const
{
  vector<int> numbers;
  vector<string> elementSymbols;

  int meta = 0;

  const size_t len = symbol.length();

  for( size_t i = 0; i < len; )
  {
    string strStr, numStr;
    while( (i<len) && isalpha( symbol[i] ) )
      strStr += tolower( symbol[i++] );

    if( strStr.length() )
    {
      //Lets check if this text could be a meta stable state indicator
      if( (strStr.length()==1 && strStr[0]=='m')
         || (strStr.length()>=4 && strStr.substr(0,4)=="meta") )
      {
        meta = 1;
        if( i < len )
        {
          const int metaLevel = atoi( symbol.c_str()+i );
          if( metaLevel>0 && metaLevel<5 )
          {
            ++i;
            meta = metaLevel;
          }//if( it was a valid metalevel following the 'M' or 'META'
        }//if( i < len )
      }else elementSymbols.push_back( strStr );
    }//if( we found a alphabetical substring )

    while( (i<len) && isdigit( symbol[i] ) )
      numStr += symbol[i++];

    if( numStr.length() )
      numbers.push_back( atoi( numStr.c_str() ) );

    if( numStr.empty() && strStr.empty() )
      ++i;
  }//for(...)


  int atomicMass = -1;

  //The user may have specified both the atomicNumber and atomicMass, so if
  //  there was more than one number specified (excluding the isomeric level,
  //  which should not have made it into 'numbers'), then assume the largest
  //  number is the atomicMass
  if( !numbers.empty() )
    atomicMass = *max_element( numbers.begin(), numbers.end() );


  if( elementSymbols.size()==1 )
  {
    string &el = elementSymbols[0];
    const bool isSymbol = ( (el.length() <= 2)
                          || (el.length()==3 && el.substr(0,2)=="uu" ) );

    if( isSymbol )
     el[0] = toupper( el[0] );

    bool foundSymbol = false;
    for( size_t i = 0; i < m_elements.size(); ++i )
    {
      const Element *element = m_elements[i];
      if( element->symbol == el || element->name == el )
      {
        foundSymbol = true;
        el = element->symbol;
        break;
      }//
    }//for( loop over m_elements to match the current string to, i )


    if( !foundSymbol )
    {
      stringstream msg;
      msg << "SandiaDecayDataBase::toNormalSymbolForm( \"" << symbol << "\" ):"
          << " A element with name or symbol '" << el << "' couldnt be found";
      throw runtime_error( msg.str() );
    }//if( !foundSymbol )

    char buffer[32];
    if( meta<=0 )
      snprintf( buffer, sizeof(buffer), "%s%i", el.c_str(), atomicMass );
    else if( meta==1 )
      snprintf( buffer, sizeof(buffer), "%s%im", el.c_str(), atomicMass );
    else
      snprintf( buffer, sizeof(buffer), "%s%im%i", el.c_str(), atomicMass, meta );
    
    const std::string nuffstr = buffer;
    vector<const Nuclide *>::const_iterator pos = find( m_nuclides, nuffstr );

    if( pos == m_nuclides.end() )
    {
      stringstream msg;
      msg << "SandiaDecayDataBase::toNormalSymbolForm( \"" << symbol << "\" ): "
          << el << " doesnt have an isotope with atmic mass " << atomicMass
          << " isomerNumber=" << meta;
      throw runtime_error( msg.str() );
    }//if( it was an invalid isotope mass )

    return nuffstr;
  }//if( elementSymbols.size()==1 && numbers.size()==1 )

  vector<string> candiateSymbls;

  for( size_t i = 0; i < m_elements.size(); ++i )
  {
    const Element *element = m_elements[i];
    const string &name = element->name;
    string symbol = element->symbol;
    symbol[0] = tolower( symbol[0] );

    for( size_t j = 0; j < elementSymbols.size(); ++j )
      if( (elementSymbols[j] == symbol)
          || (elementSymbols[j] == name) )
        candiateSymbls.push_back( element->symbol );
  }//for( loop over m_elements to match the current string to, i )

  if( candiateSymbls.size()!=1 )
  {
    stringstream msg;
    msg << "SandiaDecayDataBase::toNormalSymbolForm( \"" << symbol << "\" ): "
        << "Found " << candiateSymbls.size() << " candidate elements";
    throw runtime_error( msg.str() );
  }//if( candiateSymbls.size()!=1 )

  char buffer[32];
  if( meta<=0 )
    snprintf( buffer, sizeof(buffer), "%s%i", candiateSymbls[0].c_str(), atomicMass );
  else if( meta==1 )
    snprintf( buffer, sizeof(buffer), "%s%im", candiateSymbls[0].c_str(), atomicMass );
  else
    snprintf( buffer, sizeof(buffer), "%s%im%i", candiateSymbls[0].c_str(), atomicMass, meta );
  
  const std::string nuffstr = buffer;

  vector<const Nuclide *>::const_iterator pos = find( m_nuclides, nuffstr );

  if( (pos == m_nuclides.end()) )
  {
    stringstream msg;
    msg << "SandiaDecayDataBase::toNormalSymbolForm( \"" << symbol << "\" ): "
        << candiateSymbls[0] << " doesnt have an isotope with atomic mass " << atomicMass
        << " isomerNumber=" << meta;
    throw runtime_error( msg.str() );
  }//if( it was an invalid isotope mass )

  return nuffstr;
}//std::string SandiaDecayDataBase::toNormalSymbolForm( const std::string &symbol ) const



std::vector<NuclideTimeEvolution> SandiaDecayDataBase::getTimeEvolution( const Nuclide *parent,
                                                 const double orignal_activity  )
{
  const double n_atoms = parent->activityToNumAtoms( orignal_activity );
  vector<NuclideNumAtomsPair> parentV( 1, NuclideNumAtomsPair(parent, n_atoms) );
  return getTimeEvolution( parentV );
}//NuclideMixture getTimeEvolution(...)



std::vector<NuclideActivityPair> SandiaDecayDataBase::decay( const Nuclide *parent,
                                                             const double orignal_activity,
                                                             const double time_in_seconds )
{
  NuclideMixture mixture;
  mixture.addNuclideByActivity( parent, orignal_activity );
  return mixture.activity( time_in_seconds );
}//decay(...)


vector<NuclideActivityPair> SandiaDecayDataBase::decay( const vector<NuclideNumAtomsPair> &parents, const double time_in_seconds )
{
  NuclideMixture mixture;
  for( size_t i = 0; i < parents.size(); ++i )
    mixture.addNuclide( parents[i] );
  return mixture.activity( time_in_seconds );
}//vector<NuclideActivityPair> decay( const vector<NuclideNumAtomsPair> &parents, const double time_in_seconds )


vector<NuclideActivityPair> SandiaDecayDataBase::decay( const std::vector<NuclideActivityPair> &parents, const double time_in_seconds )
{
  NuclideMixture mixture;
  for( size_t i = 0; i < parents.size(); ++i )
    mixture.addNuclide( parents[i] );
  return mixture.activity( time_in_seconds );
}//vector<NuclideActivityPair> decay( const vector<NuclideActivityPair> &parents, const double time_in_seconds )





std::vector<NuclideTimeEvolution> SandiaDecayDataBase::getTimeEvolution( const std::vector<NuclideActivityPair> &parents )
{
  std::vector<NuclideNumAtomsPair> parents_by_num_atoms;
  for( size_t i =0; i < parents.size(); ++i )
  {
    const Nuclide *nuclide = parents[i].nuclide;
    const double n_atoms = nuclide->activityToNumAtoms( parents[i].activity );
    parents_by_num_atoms.push_back( NuclideNumAtomsPair(nuclide, n_atoms) );
  }

  return getTimeEvolution( parents_by_num_atoms );
}//vector<NuclideTimeEvolution> getTimeEvolution( const vector<NuclideActivityPair> & )


std::vector<NuclideTimeEvolution> SandiaDecayDataBase::getTimeEvolution( const std::vector<NuclideNumAtomsPair> &input )
{
  //This function needs serious cleanup!
  //Also, profiling done 20121010 shows that this function is what should really
  //  be optimized, which makes sense given its implementation

  BatemanWorkingSpace ws;
  ws.decay_coeffs.reserve( 65 );  //64 is largest in database
  ws.nuclide_path.reserve( 65 );

  //Lets make sure if the user inputed a vector with multiple entries for
  //  the same nuclide, we combine these entries
  vector<NuclideNumAtomsPair> parents;
  for( size_t i = 0; i < input.size(); ++i )
  {
    const Nuclide *inputNuclide = input[i].nuclide;

    bool added = false;
    for( size_t p = 0; !added && (p < parents.size()); ++p )
    {
      if( (*(parents[p].nuclide)) == (*inputNuclide) ) //note: could just use pointer comparison...
      {
        added = true;
        parents[p].numAtoms += input[i].numAtoms;
      }//if( parents has this element )
    }//for( loop over current parents, p )

    if( !added )
      parents.push_back( input[i] );
  }//for( loop over input, i )


  for( size_t parentNum = 0; parentNum < parents.size(); ++parentNum )
  {
    const Nuclide *parent = parents[parentNum].nuclide;
    const double n_original = parents[parentNum].numAtoms;

    for( size_t decayNum = 0; decayNum < parent->decaysToChildren.size(); ++decayNum )
    {
      const Transition *trans = parent->decaysToChildren[decayNum];
      if( trans->child )
      {
        const vector<const Transition *> decay_path( 1, trans );
        vector< vector<DecayCoefType> > coefs(1, vector<DecayCoefType>(1, n_original) );
        coefs.reserve( 128 );
        calc_bateman_coef( coefs, decay_path, ws );
      }//if( trans->child )
    }//for( loop over transition, decayNum )
  }//for( loop over parental nucildes, parentNum )


  //Now we need to stitch the paths back together, and add up contributions
  //  to each nuclide from each decay path, but we have to be careful we
  //  dont double count contributions from identical sub paths.
  typedef vector<const Nuclide *> DecayPath;
  typedef DecayPath::const_iterator DecayPathIter;
  typedef vector< DecayPath > DecayPathVec;
  typedef DecayPathVec::const_iterator DecayPathVecIter;
  typedef map<const Nuclide *,double> NuclideToMagnitudeMap;
  typedef map<const Nuclide *, NuclideToMagnitudeMap, bool(*)(const Nuclide *,const Nuclide *) > NuclToCoefMapMap;

  DecayPathVec sub_paths_used;
  NuclToCoefMapMap nuc_to_coef_map( &Nuclide::lessThan );

  const DecayPathVec  &nuclide_paths = ws.nuclide_path;
  const vector< vector< vector<DecayCoefType> > > &decay_coeffs = ws.decay_coeffs;

//  cerr << "There were " << decay_coeffs.size() << " paths found" << endl;

  for( size_t solution = 0; solution < decay_coeffs.size(); ++solution )
  {
    for( size_t i = 0; i < decay_coeffs[solution].size(); ++i )
    {
      DecayPathIter path_begin = nuclide_paths[solution].begin();
      DecayPathIter path_end   = path_begin + 1 + i;
      const DecayPath sub_path( path_begin, path_end );

      DecayPathVecIter sub_paths_begin = sub_paths_used.begin();
      DecayPathVecIter sub_paths_end   = sub_paths_used.end();

      if( find( sub_paths_begin, sub_paths_end, sub_path ) != sub_paths_end )  //20121010: This find takes up about 60% of the time of this function
        continue;

      sub_paths_used.push_back( sub_path );

      const Nuclide *nuclide = sub_path[i];

      assert( sub_path.size() > i );

      NuclideToMagnitudeMap &coeff_mapp = nuc_to_coef_map[nuclide];

      for( size_t j = 0; j <= i; ++j )
      {
        const float coeff = decay_coeffs[solution][i][j];
        const Nuclide *orig_nuc = sub_path[j];

        if( !coeff_mapp.count( orig_nuc ) )
          coeff_mapp[orig_nuc] = 0.0;

        coeff_mapp[orig_nuc] += coeff;
      }
    }//for( loop over nuclides, i, that we are forming solutions for )
  }//for( loop over solution )


  std::vector<NuclideTimeEvolution> answer;

  for( NuclToCoefMapMap::const_iterator coef_map_iter = nuc_to_coef_map.begin();
      coef_map_iter != nuc_to_coef_map.end();
      ++coef_map_iter )
  {
    const Nuclide *nuclide = coef_map_iter->first;
    answer.push_back( NuclideTimeEvolution( nuclide ) );
    NuclideTimeEvolution &evolution = answer.back();

    for( NuclideToMagnitudeMap::const_iterator coefs_iter = coef_map_iter->second.begin();
        coefs_iter != coef_map_iter->second.end();
        ++coefs_iter )
    {
      evolution.addEvolutionTerm( coefs_iter->second,  coefs_iter->first );
    }
  }//for( each nuc_to_coef_map entry )



  return answer;
}//std::vector<NuclideTimeEvolution> getTimeEvolution( std::vector<NuclideNumAtomsPair> & ) const







NuclideNumAtomsPair::NuclideNumAtomsPair( const Nuclide *_nuclide, double _numAtoms )
  :   nuclide( _nuclide ), numAtoms( _numAtoms )
{}



  NuclideActivityPair::NuclideActivityPair( const Nuclide *_nuclide, double _activity )
    : nuclide( _nuclide ), activity( _activity )
  {}


  TimeEvolutionTerm::TimeEvolutionTerm( double mag,  const Nuclide *nuc )
    : termCoeff( mag ), exponentialCoeff( nuc->decayConstant() )
  {}

  double TimeEvolutionTerm::eval( double time_in_seconds ) const
  {
    return termCoeff * exp( -time_in_seconds*exponentialCoeff );
  }


 EnergyRatePair::EnergyRatePair( double _numPerSecond, double _energy )
    : energy( _energy ), numPerSecond( _numPerSecond )
  {}


  bool EnergyRatePair::moreThanByNumPerSecond( const EnergyRatePair &lhs,
                                       const EnergyRatePair &rhs )
  {
      return (lhs.numPerSecond > rhs.numPerSecond);
  }

  bool EnergyRatePair::lessThanByEnergy( const EnergyRatePair &lhs,
                                              const EnergyRatePair &rhs )
  {
    return lhs.energy < rhs.energy;
  }


  NuclideTimeEvolution::NuclideTimeEvolution( const Nuclide *_nuclide )
     : nuclide(_nuclide){}

  double NuclideTimeEvolution::numAtoms( double time_in_seconds ) const
  {
    double n_atoms = 0.0;
    for( size_t i = 0; i < evolutionTerms.size(); ++i )
      n_atoms += evolutionTerms[i].eval( time_in_seconds );
    return max( 0.0, n_atoms );
  }

  double NuclideTimeEvolution::activity( double time_in_seconds ) const
  {
    const double n_atoms = numAtoms( time_in_seconds );
    return nuclide->numAtomsToActivity( n_atoms );
  }

  void NuclideTimeEvolution::addEvolutionTerm( double mag,  const Nuclide *nuc )
  {
    evolutionTerms.push_back( TimeEvolutionTerm( mag, nuc ) );
  }

}//namespace SandiaDecay
