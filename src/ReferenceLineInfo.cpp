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

#include <vector>
#include <string>
#include <stdexcept>

#include <Wt/WWebWidget>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"


#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "sandia_decay/SandiaDecay.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;

const int ReferenceLineInfo::sm_xmlSerializationVersion = 1;


namespace
{
  const size_t ns_numColors = 13;
  const char *ns_colors[ns_numColors]
  = { "#0000FF",
    "#006600", "#006666", "#0099FF",
    "#9933FF", "#FF66FF", "#CC3333", "#FF6633",
    "#FFFF99", "#CCFFCC", "#000CC", "#666666", "#003333"
  };
  
  template<class T> struct index_compare_assend
  {
    index_compare_assend(const T arr) : arr(arr) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      return arr[a] < arr[b];
    }
    const T arr;
  };//struct index_compare
  
  std::string jsQuote( const std::string &str )
  {
    return Wt::WWebWidget::jsStringLiteral(str,'\'');
  }
  
}//namespace


//I need practice/help memorizing background lines, so I typed the following in
//  from "Practical gamma-ray spectrometry" pg 361
const BackgroundLine BackgroundReactionLines[28] =
{
  BackgroundLine( 53.44f,   10.34f, "Ge(n,g)", BackgroundReaction, "Ge72(n,g), Ge74(n,2n)" ),
  BackgroundLine( 68.75f,   0.1f,   "Ge(n,n)", BackgroundReaction, "Ge73(n,n) broad antisymetric peak" ),
  BackgroundLine( 139.68f,  39.0f,  "Ge75m",   BackgroundReaction, "Ge74(n,g), Ge76(n,2n)" ),
  BackgroundLine( 159.7f,   10.33f, "Ge(n,g)", BackgroundReaction, "Ge76(n,g)" ),
  BackgroundLine( 174.95f,  0.0f,   "Ge(n,g)", BackgroundReaction, "Ge70(n,g) activation" ),
  BackgroundLine( 198.39f,  0.0f,   "Ge71m",   BackgroundReaction, "Sum peak Ge70(n,g)" ),
  BackgroundLine( 278.26f,  0.0f,   "Cu64",    BackgroundReaction, "Cu63(n,g), Cu65(n,2n) prompt gamma" ),
  BackgroundLine( 336.24f,  45.9f,  "Cd115m,In115m", BackgroundReaction, "Activation of Cd (daughter of Cd115)" ),
  BackgroundLine( 416.86f,  27.7f,  "In116m",  BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  BackgroundLine( 527.90f,  27.5f,  "Cd115",   BackgroundReaction, "Cd114(n,g) activation" ),
  BackgroundLine( 558.46f,  0.0f,   "Cd114",   BackgroundReaction, "Cd113(n,g) prompt gamma" ),
  BackgroundLine( 569.7f,   97.89f, "Pb207m",  BackgroundReaction, "Pb207(n,n)" ),
  BackgroundLine( 579.2f,   0.0f,   "Pb207",   BackgroundReaction, "Pb207(n,n) prompt gamma" ),
  BackgroundLine( 595.85f,  0.0f,   "Ge74",    BackgroundReaction, "Ge74(n,n) broad asymmetric peak" ),
  BackgroundLine( 669.62f,  0.0f,   "Cu63",    BackgroundReaction, "Cu63(n,n) prompt gamma" ),
  BackgroundLine( 689.6f,   0.0f,   "Ge72",    BackgroundReaction, "Ge72(n,n) broad asymetric peak" ),
  BackgroundLine( 803.06f,  0.0f,   "Pb206",   BackgroundReaction, "Pb206(n,n) prompt gamma" ),
  BackgroundLine( 843.76f,  71.8f,  "Mg27",    BackgroundReaction, "Mg26(n,g) or Al27(n,p) of encapsilation" ),
  BackgroundLine( 846.77f,  0.0f,   "Fe56",    BackgroundReaction, "Fe56(n,n)" ),
  BackgroundLine( 962.06f,  0.0f,   "Cu63",    BackgroundReaction, "Cu63(n,n) prompt gamma" ),
  BackgroundLine( 1014.44f, 28.0f,  "Mg27",    BackgroundReaction, "Mg26(n,g) or Al27(n,p) of encapsilation" ),
  BackgroundLine( 1063.66f, 88.5f,  "Pb207m",  BackgroundReaction, "Pb207(n,n)" ),
  BackgroundLine( 1097.3f,  56.2f,  "In116",   BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  BackgroundLine( 1115.56f, 0.0f,   "Cu65",    BackgroundReaction, "Cu65(n,n)" ),
  BackgroundLine( 1173.23f, 99.85f, "Co60",    BackgroundReaction, "Activation" ),
  BackgroundLine( 1293.54f, 84.4f,  "In116",   BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  BackgroundLine( 1332.49f, 99.98f, "Co60",    BackgroundReaction, "Activation" ),
  BackgroundLine( 2224.57f, 0.0f,   "H2",      BackgroundReaction, "H(n,g)" )
};//BackgroundReactionLines[]


//BackgroundLines looks to be taking up ~14 kb of executable size on Win7
const BackgroundLine BackgroundLines[89] =
{
  /*xrays bellow 46 keV not inserted*/
  BackgroundLine( 46.54f,   4.25f,  "Pb210",   Ra226Series,        "" ),
  BackgroundLine( 53.23f,   1.060f, "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 63.28f,   4.8f,   "Th234",   U238Series,         "" ),
  BackgroundLine( 72.81f,   27.7f,  "Pb xray", BackgroundXRay,     "Flourescence and Tl208 decay" ),
  BackgroundLine( 74.82f,   27.7f,  "Bi xray", BackgroundXRay,     "Pb212, Pb214 decay" ),
  BackgroundLine( 74.97f,   46.2f,  "Pb xray", BackgroundXRay,     "Flourescence and Tl208 decay" ),
  BackgroundLine( 77.11f,   46.2f,  "Bi xray", BackgroundXRay,     "Pb212, Pb214 decay" ),
  BackgroundLine( 79.29f,   46.1f,  "Po xray", BackgroundXRay,     "Fluorescence and Bi212, Bi214 decay" ),
  BackgroundLine( 81.23f,   0.9f,   "Th231",   U235Series,         "" ),
  BackgroundLine( 84.94f,   10.7f,  "Pb xray", BackgroundXRay,     "Flourescence and Tl208 decay" ),
  BackgroundLine( 87.3f,    3.91f,  "Pb xray", BackgroundXRay,     "Flourescence and Tl208 decay" ),
  BackgroundLine( 87.35f,   10.7f,  "Bi xray", BackgroundXRay,     "Pb212, Pb214 decay" ),
  BackgroundLine( 89.78f,   3.93f,  "Bi xray", BackgroundXRay,     "Pb212, Pb214 decay" ),
  BackgroundLine( 89.96f,   28.1f,  "Th xray", BackgroundXRay,     "U235 and Ac228 decay" ),
  BackgroundLine( 92.58f,   5.58f,  "Th234",   U238Series,         " - doublet" ),
  BackgroundLine( 93.35f,   45.4f,  "Th xray", BackgroundXRay,     "U235 and Ac228 decay" ),
  BackgroundLine( 105.6f,   10.7f,  "Th xray", BackgroundXRay,     "U235 and Ac228 decay" ),
  BackgroundLine( 109.16f,  1.54f,  "U235",    U235Series,         "" ),
  BackgroundLine( 112.81f,  0.28f,  "Th234",   U238Series,         "" ),
  BackgroundLine( 122.32f,  1.192f, "Ra223",   U235Series,         "" ),
  BackgroundLine( 129.06f,  2.42f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 143.76f,  10.96f, "U235",    U235Series,         "" ),
  BackgroundLine( 163.33f,  5.08f,  "U235",    U235Series,         "" ),
  BackgroundLine( 185.72f,  57.2f,  "U235",    U235Series,         "" ),
  BackgroundLine( 186.21f,  3.555f, "Ra226",   U238Series,         "" ),
  BackgroundLine( 205.31f,  5.01f,  "U235",    U235Series,         "" ),
  BackgroundLine( 209.26f,  3.89f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 238.63f,  43.6f,  "Pb212",   Th232Series,        "" ),
  BackgroundLine( 240.89f,  4.12f,  "Ra224",   Th232Series,        "" ),
  BackgroundLine( 242.0f,   7.268f, "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 269.49f,  13.7f,  "Ra223",   U235Series,         "" ),
  BackgroundLine( 270.24f,  3.46f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 277.37f,  2.37f,  "Tl208",   Th232Series,        "" ),
  BackgroundLine( 295.22f,  18.5f,  "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 299.98f,  2.16f,  "Th227",   U235Series,         "" ),
  BackgroundLine( 300.07f,  2.47f,  "Pa231",   U235Series,         "" ),
  BackgroundLine( 300.09f,  3.18f,  "Pb212",   Th232Series,        "" ),
  BackgroundLine( 328.0f,   2.95f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 338.28f,  2.79f,  "Ra223",   U235Series,         "" ),
  BackgroundLine( 338.32f,  11.27f, "Ac228",   Th232Series,        "" ),
  BackgroundLine( 351.06f,  12.91f, "Bi211",   U235Series,         "" ),
  BackgroundLine( 351.93f,  35.60f, "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 409.46f,  1.92f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 447.6f,   10.44f, "Be7",     OtherBackground,    "Cosmic" ),
  BackgroundLine( 462.0f,   0.213f, "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 463.0f,   4.40f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 510.7f,   6.29f,  "Tl208",   Th232Series,        "" ),
  BackgroundLine( 511.0f,   1.0f,   "",        OtherBackground,    "Annihilation radiation (beta+)" ),
  BackgroundLine( 570.82f,  0.182f, "Ac228",   Th232Series,        "" ),
  BackgroundLine( 583.19f,  30.6f,  "Tl208",   Th232Series,        "" ),
  BackgroundLine( 609.31f,  45.49f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 661.66f,  40.0f, "Cs137",   OtherBackground,    "Fission" ),
  BackgroundLine( 726.86f,  0.62f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 727.33f,  6.74f,  "Bi212",   Th232Series,        "" ),
  BackgroundLine( 755.31f,  1.0f,   "Ac228",   Th232Series,        "" ),
  BackgroundLine( 768.36f,  4.891f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 794.95f,  4.25f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 806.17f,  1.262f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 832.01f,  3.52f,  "Pb211",   U235Series,         "" ),
  BackgroundLine( 835.71f,  1.61f,  "AC228",   Th232Series,        "" ),
  BackgroundLine( 839.04f,  0.587f, "Pb214",   Ra226Series,        "" ),
  BackgroundLine( 860.56f,  4.48f,  "Tl208",   Th232Series,        "" ),
  BackgroundLine( 911.20f,  25.8f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 934.06f,  3.096f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 964.77f,  4.99f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 968.97f,  15.8f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 1001.03f, 1.021f, "Pa234m",  U238Series,         "..." ),
  BackgroundLine( 1120.29f, 14.91f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1155.19f, 1.635f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1238.11f, 5.827f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1377.67f, 3.967f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1407.98f, 2.389f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1459.14f, 0.83f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 1460.82f, 10.66f, "K40",     K40Background,      "" ),
  BackgroundLine( 1588.2f,  3.22f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 1592.51f, 1.0f,   "",        OtherBackground,    "Double Escape, 2614 keV" ),
  BackgroundLine( 1620.74f, 1.51f,  "Bi212",   Th232Series,        "" ),
  BackgroundLine( 1630.63f, 1.51f,  "Ac228",   Th232Series,        "" ),
  BackgroundLine( 1686.09f, 0.095f, "Ac228",   Th232Series,        "" ),
  BackgroundLine( 1661.28f,  1.12f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1729.6f,  2.843f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1764.49f, 15.28f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1838.36f, 0.346f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 1847.42f, 2.023f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 2103.51f, 2.0f,   "",        OtherBackground,    "Single Escape, 2614 keV" ),
  BackgroundLine( 2118.55f, 0.454f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 2204.21f, 4.913f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 2447.86f, 1.518f, "Bi214",   Ra226Series,        "" ),
  BackgroundLine( 2614.51f, 35.85f, "Tl208",   Th232Series,        "Th232 series; Pb208(n,p)" ),
};//BackgroundLines


ReferenceLineInfo::ReferenceLineInfo()
{
  reset();
}

void ReferenceLineInfo::reset()
{
  nuclide = NULL;
  element = NULL;
  energies.clear();
  intensities.clear();
  particlestrs.clear();
  decaystrs.clear();
  elementstrs.clear();
  particle_sf.clear();
  reactionGammas.clear();
  backgroundLines.clear();
  labelTxt.clear();
  reactionsTxt.clear();
  showGammas = showXrays = showAlphas = showBetas = false;
  promptLinesOnly = showLines = isBackground = displayLines = false;
  lowerBrCuttoff = age = 0.0;
}//void ReferenceLineInfo::reset()


bool ReferenceLineInfo::empty() const
{
  return (!nuclide && !element
          && reactionGammas.empty() && backgroundLines.empty());
}//bool empty() const


bool ReferenceLineInfo::operator==( const ReferenceLineInfo &rhs ) const
{
  return (nuclide==rhs.nuclide
          && element==rhs.element
          && labelTxt==rhs.labelTxt );
}//ReferenceLineInfo::operator==


void ReferenceLineInfo::toJson( string &json, size_t num ) const
{
  const size_t color_index = (num % ns_numColors);
  
  std::stringstream jsonstrm;
  jsonstrm << "{color:'" << ns_colors[color_index] << "',"
  << "parent:'" << parentLabel() << "',";
  
  if( promptLinesOnly )
    jsonstrm << "prompt:true,";
  if( !backgroundLines.empty() )
    jsonstrm << "age:'Primordial',";
  else
    jsonstrm << "age:" << jsQuote( PhysicalUnits::printToBestTimeUnits(age) ) << ",";
  
  if( !detectorName.empty() )
    jsonstrm << "detector:" << jsQuote( detectorName ) << ",";
  
  if( !shieldingName.empty() )
  {
    const string name = jsQuote( shieldingName );
    string thickness = PhysicalUnits::printToBestLengthUnits( shieldingThickness );
    thickness = jsQuote( thickness );
    
    jsonstrm << "shielding:" << name << ",";
    jsonstrm << "shieldingThickness:" << thickness << ",";
  }//if( !shieldingName.empty() )
  
  jsonstrm << "particleSf:{";
  
  for( map<string,double>::const_iterator iter = particle_sf.begin();
       iter != particle_sf.end(); ++iter )
  {
    jsonstrm << ((iter==particle_sf.begin()) ? "" : ",") << jsQuote(iter->first)
             << ":" << iter->second;
  }
  
  jsonstrm << "},lines:[";
  
  bool printed = false;
  for( size_t i = 0; i < energies.size(); ++i )
  {
    if( intensities[i] <= 0.0 )
      continue;
    
    if( printed )
      jsonstrm << ",";
    
    printed = true;
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.3g", intensities[i] );
    
    jsonstrm << "{e:" << energies[i]
    << ",h:" << buffer;
    if( !particlestrs[i].empty() )
      jsonstrm << ",particle:'" << particlestrs[i] << "'";
    if( !decaystrs[i].empty() )
      jsonstrm << ",decay:'" << decaystrs[i] << "'";
    if( !elementstrs[i].empty() )
      jsonstrm << ",el:'" << elementstrs[i] << "'";
    jsonstrm << "}";
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  jsonstrm <<"]}";
  
  json += jsonstrm.str();
}//std::string toJson( const ReferenceLineInfo &displnuc )


void ReferenceLineInfo::sortByEnergy()
{
  const size_t len = energies.size();
  
  if( len != energies.size() || len != intensities.size()
     || len != particlestrs.size() || len != decaystrs.size()
     || len != elementstrs.size() )
    throw runtime_error( "ReferenceLineInfo::sortByEnergy(): inconsistent input size" );
  
  vector<size_t> sort_indices( energies.size() );
  for( size_t i = 0; i < len; ++i )
    sort_indices[i] = i;
  
  std::sort( sort_indices.begin(), sort_indices.end(),
            index_compare_assend<vector<double>&>(energies) );
  
  ReferenceLineInfo tmp = *this;
  for( size_t i = 0; i < len; ++i )
  {
    const size_t index = sort_indices[i];
    energies[i]     = tmp.energies[index];
    intensities[i]  = tmp.intensities[index];
    particlestrs[i] = tmp.particlestrs[index];
    decaystrs[i]    = tmp.decaystrs[index];
    elementstrs[i]  = tmp.elementstrs[index];
  }//for( size_t i = 0; i < len; ++i )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < len; ++i )
  {
    if( energies[i-1] > energies[i] )
    {
      log_developer_error( BOOST_CURRENT_FUNCTION, "Error sorting reference photopeaks" );
      break;
    }
  }//for( size_t i = 0; i < len; ++i )
#endif
}//void PhotopeakLineDisplay::sortByEnergy()


std::string ReferenceLineInfo::parentLabel() const
{
  if( nuclide )
    return nuclide->symbol;
  else if( element )
    return element->symbol;
  else if( !reactionGammas.empty() )
    return reactionsTxt;
  else if( !backgroundLines.empty() )
    return "Background";
  return "";
}//std::string parentLabel() const


void ReferenceLineInfo::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  reset();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  rapidxml::xml_node<char> *node, *nuclide_node, *element_node, *lines_node,
  *energy_node, *br_node, *str_node, *sfnode;
  rapidxml::xml_attribute<char> *attr, *typeattr, *valattr;
  
  if( !base_node || (base_node->name() != string("DisplayedSource")) )
    throw runtime_error( "Invalid base node for DisplayedSource" );
  
  assert( sm_xmlSerializationVersion == 1 );
  
  int version;
  attr = base_node->first_attribute( "version", 7 );
  if( !attr || !attr->value()
     || !(stringstream(attr->value()) >> version)
     || (version != 0 && version != 1) )
    throw runtime_error( "Mising or invalid DisplayedSource version" );
  
  nuclide_node = base_node->first_node( "Nuclide", 7 );
  element_node = base_node->first_node( "Element", 7 );
  
  if( nuclide_node && nuclide_node->value() )
    nuclide = db->nuclide( nuclide_node->value() );
  if( element_node && element_node->value() )
    element = db->element( element_node->value() );
  
  node = base_node->first_node( "DisplayLines", 12 );
  if( node && node->value() && strlen(node->value()))
    displayLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowGammas", 10 );
  if( node && node->value() && strlen(node->value()))
    showGammas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowXrays", 9 );
  if( node && node->value() && strlen(node->value()))
    showXrays = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowAlphas", 10 );
  if( node && node->value() && strlen(node->value()))
    showAlphas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowBetas", 9 );
  if( node && node->value() && strlen(node->value()))
    showBetas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowLines", 9 );
  if( node && node->value() && strlen(node->value()))
    showLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "PromptLinesOnly", 15 );
  if( node && node->value() && strlen(node->value()))
    promptLinesOnly = (node->value()[0] == '1');
  
  node = base_node->first_node( "IsBackground", 12 );
  if( node && node->value() && strlen(node->value()))
    isBackground = (node->value()[0] == '1');
  
  node = base_node->first_node( "IsReaction", 10 );
  if( node && node->value() && strlen(node->value()))
    isReaction = (node->value()[0] == '1');
  
  node = base_node->first_node( "Age", 3 );
  if( node && node->value() && !(stringstream(node->value()) >> age) )
    throw runtime_error( "Ivalid age parameter" );
  
  node = base_node->first_node( "LowestBranchRatio", 17 );
  if( node && node->value() && !(stringstream(node->value())>> lowerBrCuttoff) )
    throw runtime_error( "Ivalid LowerBrCuttoff parameter" );
  
  node = base_node->first_node( "ShieldingName", 13 );
  if( node && node->value() && node->value_size() )
    shieldingName = node->value();

  node = base_node->first_node( "ShieldingThickness", 18 );
  if( node && node->value() && !(stringstream(node->value())>> shieldingThickness) )
    throw runtime_error( "Ivalid ShieldingThickness parameter" );

  node = base_node->first_node( "DetectorName", 12 );
  if( node && node->value() && node->value_size() )
    detectorName = node->value();
  
  node = base_node->first_node( "Nuclide", 7 );
  if( node && node->value() )
    labelTxt = node->value();
  
  if( isBackground )
  {
    for( const BackgroundLine &bl : BackgroundLines )
    {
      if( std::get<3>(bl)!=BackgroundXRay || showXrays )
        backgroundLines.push_back( &bl );
    }//for( const BackgroundLine &bl : BackgroundLines )
  }//if( isBackground )
  
  if( isReaction )
  {
    try
    {
      const ReactionGamma *rctnDb = ReactionGammaServer::database();
      if( rctnDb )
        reactionsTxt = rctnDb->gammas( labelTxt, reactionGammas );
    }catch( std::exception &e )
    {
      cerr << "ReferenceLineInfo::deSerialize() caught: " << e.what() << endl;
    }//try / catch
  }//if( isReaction )
  
  lines_node = base_node->first_node( "Lines", 5 );
  if( lines_node )
  {
    for( node = lines_node->first_node( "Line", 4 );
        node; node = node->next_sibling( "Line", 4 ) )
    {
      double energy, intensity;
      energy_node = node->first_node( "Energy", 6 );
      br_node = node->first_node( "Intensity", 9 );
      if( !energy_node || !br_node
         || !energy_node->value() || !br_node->value()
         || !(stringstream(energy_node->value()) >> energy)
         || !(stringstream(br_node->value()) >> intensity) )
        throw runtime_error( "Invalid line element" );
      
      string partstr, decaystr, elstr;
      str_node = node->first_node( "Particle", 8 );
      if( str_node && str_node->value() )
        partstr = str_node->value();
      str_node = node->first_node( "Decay", 5 );
      if( str_node && str_node->value() )
        decaystr = str_node->value();
      str_node = node->first_node( "Element", 7 );
      if( str_node && str_node->value() )
        elstr = str_node->value();
      
      energies.push_back(     energy );
      intensities.push_back(  intensity );
      particlestrs.push_back( partstr );
      decaystrs.push_back(    decaystr );
      elementstrs.push_back(  elstr );
    }//for( loop over line elements )
  }//if( lines_node )
  
  
  if( version >= 1 && particlestrs.size() )
  {
    sfnode = base_node->first_node( "ParticleScaleFactors", 20 );
    if( !sfnode )
      throw runtime_error( "Missing ParticleScaleFactors node" );
    
    for( node = sfnode->first_node( "ParticleSF", 10 );
        node; node = node->next_sibling( "ParticleSF", 10 ) )
    {
      double val;
      typeattr = node->first_attribute("type",4);
      valattr = node->first_attribute("value",5);
      if( !typeattr || !typeattr->value() || !valattr || !valattr->value() || !(stringstream(valattr->value()) >> val) )
        throw runtime_error( "Invalid ParticleSF" );
      
      particle_sf[typeattr->value()] = val;
    }
    
    for( const string &str : particlestrs )
    {
      if( !particle_sf.count(str) )
        throw runtime_error( "Particle type does not have a SF: " + str );
    }
  }else
  {
    for( size_t i = 0; i < particlestrs.size(); ++i )
      particle_sf[particlestrs[i]] = 1.0;
  }//if( version >= 1 ) / else
}//void deSerialize( const rapidxml::xml_node<char> *node );


void ReferenceLineInfo::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  char buffer[64];
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *lines_node, *sf_node,
  *part_sf_node, *energy_node, *br_node, *str_node;
  rapidxml::xml_attribute<char> *attr;
  
  name = "DisplayedSource";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  assert( sm_xmlSerializationVersion == 1 );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  value = doc->allocate_string( buffer );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  /*
   if( nuclide )
   {
   name = "Nuclide";
   //    value = nuclide->symbol.c_str(); never read
   node = doc->allocate_node( rapidxml::node_element, name );
   base_node->append_node( node );
   }//if( nuclide )
   
   if( element )
   {
   name = "Element";
   //    value = element->symbol.c_str(); //never read
   node = doc->allocate_node( rapidxml::node_element, name );
   base_node->append_node( node );
   }//if( element )
   
   if( reactionGammas.size() && reactionGammas[0].reaction )
   {
   name = "Reaction";
   //    value = reactionGammas[0].reaction->name().c_str(); //never read
   node = doc->allocate_node( rapidxml::node_element, name );
   base_node->append_node( node );
   }//if( reactionGammas.size() && reactionGammas[0].reaction )
   */
  
  //  if( backgroundLines.size() )
  //  {
  //    name = "BackgroundLines";
  //    switch( backgroundLines[0]->get<3>() )
  //    {
  //      case U238Series:  value = "U238 series"; break;
  //      case U235Series:  value = "U235 series"; break;
  //      case Th232Series: value = "Th232 series"; break;
  //      case Ra226Series: value = "U238 (Ra226) series"; break;
  //      case K40Background: value = "Primordial"; break;
  //      case OtherBackground: case BackgroundXRay: case BackgroundReaction:
  //        value = backgroundLines[0]->get<2>().c_str();
  //      break;
  //    }//switch( backgroundLines[i]->get<3>() )
  //    node = doc->allocate_node( rapidxml::node_element, name );
  //    base_node->append_node( node );
  //  }//if( backgroundLines.size() )
  
  name = "Nuclide";
  value = labelTxt.c_str();
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  if( age >= 0.0 )
  {
    name = "Age";
    const string agestrval = PhysicalUnits::printToBestTimeUnits(age);
    value = doc->allocate_string( agestrval.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( age >= 0.0 )
  
  name = "LowestBranchRatio";
  snprintf( buffer, sizeof(buffer), "%g", lowerBrCuttoff );
  value = doc->allocate_string( buffer );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "DisplayLines";
  value = (displayLines ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowGammas";
  value = (showGammas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowXrays";
  value = (showXrays ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowAlphas";
  value = (showAlphas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowBetas";
  value = (showBetas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "PromptLinesOnly";
  value = (promptLinesOnly ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowLines";
  value = (showLines ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "IsBackground";
  value = (isBackground ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "IsReaction";
  value = (isReaction ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  if( shieldingName.size() )
  {
    name = "ShieldingName";
    value = shieldingName.c_str();
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
    
    if( shieldingThickness > 0.0 )
    {
      name = "ShieldingThickness";
      snprintf( buffer, sizeof(buffer), "%g", shieldingThickness );
      value = doc->allocate_string( buffer );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
    }
  }//if( shieldingName.size() )
  
  if( detectorName.size() )
  {
    name = "DetectorName";
    value = detectorName.c_str();
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( detectorName.size() )
  
  
  if( energies.size() )
  {
    if( (energies.size() != intensities.size())
       || (energies.size() != particlestrs.size())
       || (energies.size() != decaystrs.size())
       || (energies.size() != elementstrs.size()) )
      throw runtime_error( "energies.size() != intensities.size() or related" );
    
    name = "Lines";
    lines_node = doc->allocate_node( rapidxml::node_element, name );
    base_node->append_node( lines_node );
    
    for( size_t i = 0; i < energies.size(); ++i )
    {
      using rapidxml::node_element;
      
      node = doc->allocate_node( node_element, "Line" );
      lines_node->append_node( node );
      
      snprintf( buffer, sizeof(buffer), "%g", energies[i] );
      value = doc->allocate_string( buffer );
      energy_node = doc->allocate_node( node_element, "Energy", value );
      node->append_node( energy_node );
     
       snprintf( buffer, sizeof(buffer), "%g", intensities[i] );
      value = doc->allocate_string( buffer );
      br_node = doc->allocate_node( node_element, "Intensity", value );
      node->append_node( br_node );
      
      if( particlestrs[i].size() )
      {
        value = doc->allocate_string( particlestrs[i].c_str() );
        str_node = doc->allocate_node( node_element, "Particle", value );
        node->append_node( str_node );
      }//if( particlestrs[i].size() )
      
      if( decaystrs[i].size() )
      {
        value = doc->allocate_string( decaystrs[i].c_str() );
        str_node = doc->allocate_node( node_element, "Decay", value );
        node->append_node( str_node );
      }//if( decaystrs[i].size() )
      
      if( elementstrs[i].size() )
      {
        value = doc->allocate_string( elementstrs[i].c_str() );
        str_node = doc->allocate_node( node_element, "Element", value );
        node->append_node( str_node );
      }//if( elementstrs[i].size() )
    }//for( size_t i = 0; i < energies.size(); ++i )
  }//if( energies.size() )
  
  if( particlestrs.size() )
  {
    assert( sm_xmlSerializationVersion >= 1 );
    sf_node = doc->allocate_node( rapidxml::node_element, "ParticleScaleFactors" );
    base_node->append_node( sf_node );
    
    for( std::map<std::string,double>::const_iterator iter = particle_sf.begin();
        iter != particle_sf.end(); ++iter )
    {
      snprintf( buffer, sizeof(buffer), "%g", iter->second );
      part_sf_node = doc->allocate_node( rapidxml::node_element, "ParticleSF" );
      sf_node->append_node( part_sf_node );
      value = doc->allocate_string( iter->first.c_str(), iter->first.size() );
      attr = doc->allocate_attribute( "type", value, 4, iter->first.size() );
      part_sf_node->append_attribute( attr );
      
      value = doc->allocate_string( buffer );
      attr = doc->allocate_attribute( "value", value );
      part_sf_node->append_attribute( attr );
    }
  }//if( particlestrs.size() )
}//void ReferenceLineInfo::serialize( rapidxml::xml_node<char> *parent_node )



